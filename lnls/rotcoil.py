"""library to read rotation coil data and create excitation data file."""

import os as _os
import numpy as _np

from siriuspy import envars as _envars
from siriuspy import util as _util
from siriuspy.magnet import util as _mutil
from siriuspy.ramp import util as _rutil
from matplotlib import gridspec as _gridspec


class RotCoilData:
    """Rotating coild data."""

    _del = ('(C)', '(rps)', '(rps^2)', '(A)', '(V)', '(ohm)', '(m)', '(um)')
    _params = (
        'file',
        'date',
        'hour',
        'operator',
        'software_version',
        'bench',
        'temperature(C)',
        'integrator_gain',
        'n_integration_points',
        'velocity(rps)',
        'acceleration(rps^2)',
        'n_collections',
        'n_turns',
        'analysis_interval',
        'rotation',
        'main_coil_current_avg(A)',
        'main_coil_current_std(A)',
        'main_coil_volt_avg(V)',
        'main_coil_volt_std(V)',
        'magnet_resistance_avg(ohm)',
        'magnet_resistance_std(ohm)',
        'ch_coil_current_avg(A)',
        'ch_coil_current_std(A)',
        'cv_coil_current_avg(A)',
        'cv_coil_current_std(A)',
        'qs_coil_current_avg(A)',
        'qs_coil_current_std(A)',
        'trim_coil_current_avg(A)',
        'trim_coil_current_std(A)',
        'rotating_coil_name',
        'rotating_coil_type',
        'measurement_type',
        'pulse_start_collect',
        'n_turns_main_coil',
        'main_coil_internal_radius(m)',
        'main_coil_external_radius(m)',
        'n_turns_bucked_coil',
        'bucked_coil_internal_radius(m)',
        'bucked_coil_external_radius(m)',
        'magnetic_center_x(um)',
        'magnetic_center_y(um)', )

    def __init__(self, path, conv_mpoles_sign):
        """Init."""
        self.path = path
        self._read_data(conv_mpoles_sign)

    def _read_data(self, conv_mpoles_sign):
        # read all text
        with open(self.path, 'r', encoding='latin1') as f:
            text = f.read()
        lines = text.splitlines()

        # process header file
        self._process_header(lines)

        # process data
        self._process_data(lines, conv_mpoles_sign)

    def _process_header(self, lines):
        # process header file
        for line in lines:
            param, *data = line.replace('\t', ' ').strip().split(' ')
            for p in RotCoilData._params:
                if param == p:

                    # delete unwanted substrings
                    param = self._del_unwanted(param)

                    # try to convert data to int or float
                    data = ' '.join(data).strip()
                    try:
                        data = int(data)
                    except ValueError:
                        try:
                            data = float(data)
                        except ValueError:
                            pass

                    # finally create attribute
                    setattr(self, param, data)

    def _process_data(self, lines, conv_mpoles_sign):
        self.harmonics = list()
        self.intmpole_normal_avg = list()
        self.intmpole_skew_avg = list()
        n1 = lines.index('##### Reading Data #####') + 3
        try:
            n2 = lines.index('##### Raw Data Stored(V.s) [1e-12] #####') - 5
        except ValueError:
            n2 = lines.index('##### Raw Data Stored(V.s) #####') - 5
        for i in range(n1, n2):
            words = lines[i].replace('\t', ' ').strip().split()
            n = int(words[0])
            multipoles = [conv_mpoles_sign * float(w) for w in words[1:]]
            self.harmonics.append(n)
            self.intmpole_normal_avg.append(multipoles[0])
            self.intmpole_skew_avg.append(multipoles[2])
            # print(n, multipoles)

    @staticmethod
    def _del_unwanted(param):
        for r in RotCoilData._del:
            param = param.replace(r, '')
        return param


class RotCoilMeas_Cor:
    """Rotation coil measurement of corrector magnets."""

    main_harmonic = 1  # 1: dipole, 2:quadrupole, etc...
    pwrsupply_polarity = 'bipolar'


class RotCoilMeas_HCor(RotCoilMeas_Cor):
    """Rotation coil measurement of horizontal corrector magnets."""

    main_harmonic_type = 'normal'


class RotCoilMeas_VCor(RotCoilMeas_Cor):
    """Rotation coil measurement of vertical corrector magnets."""

    main_harmonic_type = 'skew'


class RotCoilMeas_Quad:
    """Rotation coil measurement of quadrupole magnets."""

    main_harmonic = 2  # 1: dipole, 2:quadrupole, etc...
    main_harmonic_type = 'normal'
    pwrsupply_polarity = 'monopolar'


class RotCoilMeas_Sext:
    """Rotation coil measurement of sextupole magnets."""

    main_harmonic = 3  # 1: dipole, 2:quadrupole, etc...
    main_harmonic_type = 'normal'
    pwrsupply_polarity = 'monopolar'


class RotCoilMeas:
    """Rotation coil measurement of SI magnets."""

    excitation_type = 'main'
    family_folder = ''
    lnls_ima_path = _envars.folder_lnls_ima

    _excdata_obs = (
        '# POLARITY TABLE',
        '# ==============',
        '#',
        ('# Magnet function         | IntStrength(1) | IntField(2) |'
         ' ConvSign(3) | Current(4)'),
        ('# ------------------------|----------------|-------------|'
         '-------------|-----------'),
        ('# dipole                  | Angle > 0      | BYL  < 0    |'
         ' -1.0        | I > 0'),
        ('# corrector-horizontal    | HKick > 0      | BYL  > 0    |'
         ' +1.0        | I > 0'),
        ('# corrector-vertical      | VKick > 0      | BXL  < 0    |'
         ' -1.0        | I > 0'),
        ('# quadrupole (focusing)   | KL    > 0      | D1NL < 0    |'
         ' -1.0        | I > 0'),
        ('# quadrupole (defocusing) | KL    < 0      | D1NL > 0    |'
         ' -1.0        | I > 0'),
        ('# quadrupole (skew)       | KL    > 0      | D1SL > 0    |'
         ' +1.0        | I > 0'),
        ('# sextupole  (focusing)   | SL    > 0      | D2NL < 0    |'
         ' -1.0        | I > 0'),
        ('# sextupole  (defocusing) | SL    < 0      | D2NL > 0    |'
         ' -1.0        | I > 0'),
        '#',
        '# Defs:',
        '# ----',
        '# BYL   := \\int{dz By|_{x=y=0}}.',
        '# BXL   := \\int{dz Bx|_{x=y=0}}.',
        '# D1NL  := \\int{dz \\frac{dBy}{dx}_{x=y=0}}',
        '# D2NL  := (1/2!) \\int{dz \\frac{d^2By}{dx^2}_{x=y=0}}',
        '# D1SL  := \\int{dz \\frac{dBx}{dx}_{x=y=0}}',
        '# Brho  := magnetic rigidity.',
        '# Angle := ConvSign * BYL / abs(Brho)',
        '# HKick := ConvSign * BYL / abs(Brho)',
        '# VKick := ConvSign * BXL / abs(Brho)',
        '# KL    := ConvSign * D1NL / abs(Brho)',
        '# SL    := ConvSign * D2NL / abs(Brho)',
        '#',
        '# Obs:',
        '# ---',
        '# (1) Parameter definition.',
        ('#     IntStrength values correspond to integrated PolynomA and '
         'PolynomB parameters'),
        ('#     of usual beam tracking codes, with the exception that VKick '
         'has its sign'),
        '#     reversed with respecto to its corresponding value in PolynomA.',
        '# (2) Sirius coordinate system and Lorentz force.',
        '# (3) Conversion sign for IntField <-> IntStrength',
        ('# (4) Convention of magnet excitation polarity, so that when'
         ' I > 0 the strength'),
        '#     of the magnet has the expected conventional sign.',
        '',
        '# STATIC DATA FILE FORMAT',
        '# =======================',
        '#',
        ('# These static data files should comply with the following '
         'formatting rules:'),
        ('# 1. If the first alphanumeric character of the line is not the'
         ' pound sign'),
        '#    then the lines is a comment.',
        '# 2. If the first alphanumeric character is "#" then if',
        ('#    a) it is followed by "[<parameter>] <value>" a parameter '
         'names <parameter>'),
        ('#       is define with value <value>. if the string <value> has '
         'spaces in it'),
        '#       it is split as a list of strings.',
        '#    b) otherwise the line is ignored as a comment line.',)

    def __init__(self, serial_number):
        """Init."""
        self.serial_number = serial_number
        self._read_rotcoil_data()
        self._calc_magnetic_center()

    @property
    def data_sets(self):
        """Return list of data set."""
        return self._get_data_sets()

    @property
    def size(self):
        """Return number of current measurements."""
        data = self._rotcoildata[self.data_sets[0]]
        return len(data)

    def get_nominal_main_intmpole_values(self, energy):
        """Nominal integrated main multipole."""
        brho, *_ = _util.beam_rigidity(energy)
        intmpole = dict()
        for fam, strength in self.nominal_KL_values.items():
            intmpole[fam] = - strength * brho
        return intmpole

    def get_data_set_measurements(self, data_set):
        """."""
        return self._rotcoildata[data_set]

    def get_max_current_index(self):
        """Return max current index."""
        max_c_i = []
        for data_set in self.data_sets:
            # data = self._rotcoildata[data_set]
            currents = self.get_currents(data_set)
            max_c = max(currents)
            # print(self.serial_number, data_set, max_c)
            max_c_i.append(currents.index(max_c))
        umaxci = _np.unique(max_c_i)
        # print(umaxci)
        if len(umaxci) > 1:
            raise ValueError('Inconsistent current values in data sets')
        return umaxci[0]

    def get_rampup(self, data_set):
        """Rampup data."""
        c = self.get_currents(data_set)
        if self.main_harmonic_type == 'normal':
            gl = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            gl = self.get_intmpole_skew_avg(data_set, self.main_harmonic)
        ind = self.get_rampup_indices()
        return [c[i] for i in ind], [gl[i] for i in ind]

    def get_rampdown_hysteresis(self, data_set):
        """Rampdown hysteresis."""
        c = self.get_currents(data_set)
        if self.main_harmonic_type == 'normal':
            gl = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            gl = self.get_intmpole_skew_avg(data_set, self.main_harmonic)

        ind = self.get_rampup_indices()
        c_lin, gl_lin = zip(*[(c[i], gl[i]) for i in ind])
        # i_max = self.get_max_current_index()
        # c_lin = c[0:i_max+1]
        # gl_lin = gl[0:i_max+1]
        gl_int = _np.interp(c, c_lin, gl_lin)
        gl_dif = [gl_int[i] - gl[i] for i in range(len(c))]
        # gl, c, h = gl[i_max:], c[i_max:], gl_dif[i_max:]
        ind = self.get_rampdown_indices()
        gl, c, h = zip(*[(gl[i], c[i], gl_dif[i]) for i in ind])

        area = -_np.trapz(h, c)
        return gl, c, h, area

    def get_currents(self, data_set):
        """Return currents of a data set."""
        return self.get_currents_avg(data_set)

    def get_currents_avg(self, data_set):
        """Return currents of a data set."""
        data = self._rotcoildata[data_set]
        return [d.main_coil_current_avg for d in data]

    def get_currents_std(self, data_set):
        """Return currents of a data set."""
        data = self._rotcoildata[data_set]
        return [d.main_coil_current_std for d in data]

    def get_intmpole_normal_avg(self, data_set, n):
        """Return average integrated normal multipole."""
        i = self.harmonics.index(n)
        data = self._rotcoildata[data_set]
        p = []
        for datum in data:
            p.append(datum.intmpole_normal_avg[i])
        return p

    def get_intmpole_normal_avg_current(self, data_set, idx):
        """."""
        data = self._rotcoildata[data_set][idx]
        p = []
        for i in range(len(data.harmonics)):
            p.append(data.intmpole_normal_avg[i])
        return p

    def get_intmpole_skew_avg(self, data_set, n):
        """Return average integrated skew multipole."""
        i = self.harmonics.index(n)
        data = self._rotcoildata[data_set]
        p = []
        for datum in data:
            p.append(datum.intmpole_skew_avg[i])
        return p

    def get_intmpole_skew_avg_current(self, data_set, idx):
        """."""
        data = self._rotcoildata[data_set][idx]
        p = []
        for i in range(len(data.harmonics)):
            p.append(data.intmpole_skew_avg[i])
        return p

    def get_magnetic_center_x(self, data_set):
        """List with horizontal position of magnetic center."""
        data = self._rotcoildata[data_set]
        x = [d.magnetic_center_x for d in data]
        return x

    def get_magnetic_center_y(self, data_set):
        """List with vertical position of magnetic center."""
        data = self._rotcoildata[data_set]
        y = [d.magnetic_center_y for d in data]
        return y

    def rampup_curr_2_main_mpole(self, data_set, current):
        """Interpolate."""
        c, gl = self.get_rampup(data_set)
        gl_interp = _np.interp(current, c, gl)
        return gl_interp

    def rampup_main_mpole_2_curr(self, data_set, main_mpole):
        """Interpolate."""
        c, gl = self.get_rampup(data_set)
        c, gl = _np.array(c), _np.array(gl)
        arr1inds = gl.argsort()
        sorted_c = c[arr1inds]
        sorted_gl = gl[arr1inds]
        current_interp = _np.interp(main_mpole, sorted_gl, sorted_c)
        return current_interp

    def save_excdata(self, data_set, harmonics=None):
        """Save data."""
        lines = self._excitation_text(data_set, harmonics)
        filename = self.magnet_type_name + '-' + self.serial_number
        # save data to file
        with open(filename + '.txt', 'w') as fp:
            for line in lines:
                fp.write(line + '\n')

    def multipoles_kicks_spec_sys(self, data_set, current_index,
                                  energy, nrpoints=301):
        """."""
        brho, *_ = _util.beam_rigidity(energy)
        r0 = self.spec_r0 / 1000.0
        main_harm = self.main_harmonic - 1
        if self.main_harmonic_type == 'normal':
            mpole = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            mpole = self.get_intmpole_skew_avg(data_set, self.main_harmonic)
        main_mpole = mpole[current_index]
        normal_harms = self.spec_normal_sys_harms - 1
        normal_mpoles = self.spec_normal_sys_mpoles
        skew_harms = self.spec_skew_sys_harms - 1
        skew_mpoles = self.spec_skew_sys_mpoles
        return RotCoilMeas._get_kick(brho, r0,
                                     main_harm,
                                     main_mpole,
                                     normal_harms,
                                     normal_mpoles,
                                     skew_harms,
                                     skew_mpoles,
                                     nrpoints,
                                     True)

    def multipoles_kicks_spec_rms(self, data_set, current_index,
                                  energy, nrpoints=101, nrmpoles=1000):
        """."""
        brho, *_ = _util.beam_rigidity(energy)
        r0 = self.spec_r0 / 1000.0
        main_harm = self.main_harmonic - 1
        if self.main_harmonic_type == 'normal':
            mpole = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            mpole = self.get_intmpole_skew_avg(data_set, self.main_harmonic)
        main_mpole = mpole[current_index]

        kickx_min, kickx_max = None, None
        kicky_min, kicky_max = None, None
        for j in range(500):

            # --- normal multipoles
            nharms = set()
            nharms.update(self.spec_normal_sys_harms)
            nharms.update(self.spec_normal_rms_harms)
            nharms = _np.array(sorted(nharms))
            normal_mpoles = _np.zeros(nharms.shape)
            # add sys multipoles
            for i in range(len(nharms)):
                if nharms[i] in self.spec_normal_sys_harms:
                    idx = _np.argwhere(self.spec_normal_sys_harms == nharms[i])
                    idx = idx[0][0]
                    normal_mpoles[i] += self.spec_normal_sys_mpoles[idx]
            # add rms multipoles
            for i in range(len(nharms)):
                if nharms[i] in self.spec_normal_rms_harms:
                    idx = _np.argwhere(self.spec_normal_rms_harms == nharms[i])
                    idx = idx[0][0]
                    normal_mpoles[i] += self.spec_normal_rms_mpoles[idx] * \
                        2.0 * (_np.random.random() - 0.5)

            # --- skew multipoles
            sharms = set()
            sharms.update(self.spec_skew_sys_harms)
            sharms.update(self.spec_skew_rms_harms)
            sharms = _np.array(sorted(nharms))
            skew_mpoles = _np.zeros(sharms.shape)
            # add sys multipoles
            for i in range(len(sharms)):
                if sharms[i] in self.spec_skew_sys_harms:
                    idx = _np.argwhere(self.spec_skew_sys_harms == sharms[i])
                    idx = idx[0][0]
                    skew_mpoles[i] += self.spec_skew_sys_mpoles[idx]
            # add rms multipoles
            for i in range(len(sharms)):
                if sharms[i] in self.spec_skew_rms_harms:
                    idx = _np.argwhere(self.spec_skew_rms_harms == sharms[i])
                    idx = idx[0][0]
                    skew_mpoles[i] += self.spec_skew_rms_mpoles[idx] * \
                        2.0 * (_np.random.random() - 0.5)

            x, y, kickx, kicky = RotCoilMeas._get_kick(brho, r0,
                                                       main_harm,
                                                       main_mpole,
                                                       nharms - 1,
                                                       normal_mpoles,
                                                       sharms - 1,
                                                       skew_mpoles,
                                                       nrpoints,
                                                       True)
            if kickx_min is None:
                kickx_min, kickx_max = kickx, kickx
                kicky_min, kicky_max = kicky, kicky
            else:
                kickx_min = _np.min([kickx_min, kickx], axis=0)
                kickx_max = _np.max([kickx_max, kickx], axis=0)
                kicky_min = _np.min([kicky_min, kicky], axis=0)
                kicky_max = _np.max([kicky_max, kicky], axis=0)

        return x, y, [kickx_min, kickx_max], [kicky_min, kicky_max]

    def multipoles_kicks_residual_old(self, data_set, current_index,
                                      energy,
                                      include_dipole=False,
                                      include_quadrupole=True,
                                      nrpoints=301):
        """."""
        brho, *_ = _util.beam_rigidity(energy)
        r0 = self.spec_r0 / 1000.0
        main_harm = self.main_harmonic - 1
        if self.main_harmonic_type == 'normal':
            mpole = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            mpole = self.get_intmpole_skew_avg(data_set, self.main_harmonic)
        main_mpole = mpole[current_index]

        normal_harms = _np.array(self.harmonics) - 1
        skew_harms = _np.array(self.harmonics) - 1
        normal_mpoles, skew_mpoles = [], []
        for i in range(len(self.harmonics)):
            h = self.harmonics[i]
            if h != self.main_harmonic:
                nmpole = self.get_intmpole_normal_avg(data_set, h)
                smpole = self.get_intmpole_skew_avg(data_set, h)
                if (not include_dipole and h == 1) or \
                   (not include_quadrupole and h == 2):  # NOTE: TEST !!!!
                    # does not include dipolar error
                    normal_mpoles.append(0.0)
                    skew_mpoles.append(0.0)
                else:
                    normal_mpoles.append(nmpole[current_index])
                    skew_mpoles.append(smpole[current_index])
            else:
                if self.main_harmonic_type == 'normal':
                    smpole = self.get_intmpole_skew_avg(data_set, h)
                    normal_mpoles.append(0.0)
                    skew_mpoles.append(smpole[current_index])
                else:
                    nmpole = self.get_intmpole_normal_avg(data_set, h)
                    normal_mpoles.append(nmpole[current_index])
                    skew_mpoles.append(0.0)
        return RotCoilMeas._get_kick(brho, r0,
                                     main_harm,
                                     main_mpole,
                                     normal_harms,
                                     normal_mpoles,
                                     skew_harms,
                                     skew_mpoles,
                                     nrpoints,
                                     False)

    def multipoles_kicks_residual(self, data_set, current_index,
                                  energy,
                                  excluded_monomials_norm=None,
                                  excluded_monomials_skew=None,
                                  nrpoints=301):
        """."""
        if excluded_monomials_norm is None:
            excluded_monomials_norm = []
        if excluded_monomials_skew is None:
            excluded_monomials_skew = []
        brho, *_ = _util.beam_rigidity(energy)
        r0 = self.spec_r0 / 1000.0
        main_harm = self.main_harmonic - 1
        if self.main_harmonic_type == 'normal':
            mpole = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            mpole = self.get_intmpole_skew_avg(data_set, self.main_harmonic)
        main_mpole = mpole[current_index]

        normal_harms = _np.array(self.harmonics) - 1
        skew_harms = _np.array(self.harmonics) - 1
        normal_mpoles, skew_mpoles = [], []
        for i in range(len(self.harmonics)):
            h = self.harmonics[i]
            nmpole = self.get_intmpole_normal_avg(data_set, h)
            smpole = self.get_intmpole_skew_avg(data_set, h)
            if h in excluded_monomials_norm:
                normal_mpoles.append(0.0)
            else:
                normal_mpoles.append(nmpole[current_index])
            if h in excluded_monomials_skew:
                skew_mpoles.append(0.0)
            else:
                skew_mpoles.append(smpole[current_index])
        return RotCoilMeas._get_kick(brho, r0,
                                     main_harm,
                                     main_mpole,
                                     normal_harms,
                                     normal_mpoles,
                                     skew_harms,
                                     skew_mpoles,
                                     nrpoints,
                                     False)

    @staticmethod
    def _get_kick(brho, r0,
                  main_harm,
                  main_mpole,
                  normal_harms,
                  normal_mpoles,
                  skew_harms,
                  skew_mpoles,
                  nrpoints,
                  denormalize=True):
        """."""
        # print('brho: {}'.format(brho))
        # print('main_harm: {}'.format(brho))

        x = _np.linspace(-r0, r0, nrpoints)
        y = 0 * x
        z = x + 1j * y
        b = 0 * z

        h = normal_harms
        if denormalize:
            m = (main_mpole * r0**main_harm / r0**h) * normal_mpoles
        else:
            m = _np.array(normal_mpoles)
        for i in range(len(h)):
            b += m[i] * (z**h[i])

        h = skew_harms
        if denormalize:
            m = (main_mpole * r0**main_harm / r0**h) * skew_mpoles * 1j
        else:
            m = _np.array(skew_mpoles) * 1j
        for i in range(len(h)):
            b += m[i] * (z**h[i])

        kickx = - (_np.real(b) / brho)
        kicky = + (_np.imag(b) / brho)

        return x, y, kickx, kicky

    @staticmethod
    def get_excdata_text(
            pwrsupply_polarity,
            magnet_type_label,
            magnet_serial_number,
            data_set,
            main_harmonic,
            main_harmonic_type,
            harmonics,
            currents,
            mpoles_n,
            mpoles_s,
            filename):
        """Return excitation data text."""
        # harmonics = ' '.join([str(h) for h in sorted(harmonics)])
        # main_harmonic = (main_harmonic-1, main_harmonic_type)
        units = ''
        for h in harmonics:
            unit = _mutil.get_multipole_si_units(h-1)
            units += unit + ' ' + unit + '  '
        units = units.strip()
        harms = ' '.join((str(h-1) for h in harmonics))

        lines = list()
        a = lines.append
        # HEADER
        a('# HEADER')
        a('# ======')
        a('# label           {}'.format(filename))
        a('# harmonics       {}'.format(harms))
        a('# main_harmonic   {} {}'.format(main_harmonic-1,
                                           main_harmonic_type))
        a('# units           Ampere  {}'.format(units))
        a('')
        # EXCITATION DATA
        a('# EXCITATION DATA')
        a('# ===============')
        # excdata for bipolar pwrsupplies replicate data with
        # negative currents and multipoles.
        if pwrsupply_polarity == 'bipolar':
            for i in reversed(range(len(currents))):
                v = '{:+010.4f}  '.format(-currents[i])
                for j in range(len(harmonics)):
                    v += '{:+11.4e} {:+11.4e}  '.format(-mpoles_n[i, j],
                                                        -mpoles_s[i, j])
                a(v.strip())
        for i in range(len(currents)):
            v = '{:+010.4f}  '.format(currents[i])
            for j in range(len(harmonics)):
                v += '{:+11.4e} {:+11.4e}  '.format(mpoles_n[i, j],
                                                    mpoles_s[i, j])
            a(v.strip())
        a('')
        # COMMENTS
        a('# COMMENTS')
        a('# ========')
        a('# 1. file generated automatically from rotating coil '
          'measurement data')
        a('# 2. timestamp: {}'.format(_util.get_timestamp()))
        a('# 3. magnet_type_label: {}'.format(magnet_type_label))
        if magnet_serial_number is not None:
            a('# 4. magnet_serial_number: {}'.format(magnet_serial_number))
        else:
            a('# 4. magnet_serial_number: {}'.format('AVERAGE OF MAGNETS'))
        if data_set is not None:
            a('# 5. data_set: {}'.format(data_set))
        else:
            a('# 5. data_set: {}'.format('UNDEFINED'))
        a('')
        # OBS
        for line in RotCoilMeas._excdata_obs:
            a(line)
        return lines

    def get_files(self, data_set):
        """Return list of data files in a data set."""
        return self._get_files(data_set)

    def get_rampup_indices(self):
        """."""
        if self.magnet_type_label == 'Q30' and self.serial_number == '011':
            return list(self._specialized_rampupind_Q30_011())
        elif self.magnet_type_label == 'BC':
            return list(self._specialized_rampupind_BC())
        elif self.magnet_type_label == 'TBC':
            return list(self._specialized_rampupind_TBC())
        elif self.magnet_type_label == 'TBQ':
            return list(self._specialized_rampupind_TBQ())
        else:
            idx = self.get_max_current_index()
            return list(range(idx+1))

    def get_rampdown_indices(self):
        """."""
        if self.magnet_type_label == 'Q30' and self.serial_number == '011':
            return tuple(range(37, 49+1))
        elif self.magnet_type_label == 'BC':
            return self._specialized_rampdownind_BC()
        elif self.magnet_type_label == 'TBC':
            return self._specialized_rampdownind_TBC()
        elif self.magnet_type_label == 'TBQ':
            return self._specialized_rampdownind_TBQ()
        else:
            idx = self.get_max_current_index()
            return tuple(range(idx, self.size))

    def _calc_magnetic_center(self):
        # B = D + Q*z + S*z**2
        #
        # B = By + Bx * 1j
        # z = x + y * 1j
        #
        # Dipolar root for quadrupoles:
        # B(z0) = 0 => z0 = -D/Q
        #
        # Quadrupolar root for sextupoles:
        # B = (D - S*z0**2) + S*(z - z0)**2
        # z0 = -Q/(2S)
        for data_set in self._rotcoildata:
            for d in self._rotcoildata[data_set]:
                idx_dip = d.harmonics.index(1)
                idx_quad = d.harmonics.index(2)
                a0 = d.intmpole_skew_avg[idx_dip]
                b0 = d.intmpole_normal_avg[idx_dip]
                a1 = d.intmpole_skew_avg[idx_quad]
                b1 = d.intmpole_normal_avg[idx_quad]
                D = b0 + a0 * 1j
                Q = b1 + a1 * 1j
                if isinstance(self, RotCoilMeas_Quad):
                    z0 = -D/Q
                elif isinstance(self, RotCoilMeas_Sext):
                    idx_sext = d.harmonics.index(3)
                    a2 = d.intmpole_skew_avg[idx_sext]
                    b2 = d.intmpole_normal_avg[idx_sext]
                    S = b2 + a2 * 1j
                    z0 = -Q/S/2.0
                    # B = D - S*z0**2
                    # d.magnetic_center_intbx = B.imag
                    # d.magnetic_center_intby = B.real
                elif isinstance(self, RotCoilMeas_Cor):
                    z0 = 0 + 0 * 1j
                else:
                    raise NotImplementedError()
                d.magnetic_center_x = 1e6 * z0.real
                d.magnetic_center_y = 1e6 * z0.imag

    def _excitation_text(self, data_set, harmonics):

        pwrsupply_polarity = self.pwrsupply_polarity
        main_harmonic = int(self.main_harmonic)
        main_harmonic_type = self.main_harmonic_type
        if harmonics is None:
            harmonics = sorted([int(h) for h in self.harmonics])
        magnet_type_label = self.magnet_type_label
        magnet_serial_number = self.serial_number
        filename = self.magnet_type_name + '-' + magnet_serial_number
        units = ''
        for h in harmonics:
            unit = _mutil.get_multipole_si_units(h-1)
            units += unit + ' ' + unit + '  '
        units = units.strip()

        currents, _ = self.get_rampup(data_set)
        shape = (len(currents), len(harmonics))
        mpoles_n = _np.zeros(shape)
        mpoles_s = _np.zeros(shape)
        idx = self.get_rampup_indices()
        for j in range(len(harmonics)):
            h = harmonics[j]
            n = self.get_intmpole_normal_avg(data_set, h)
            n = [n[i] for i in idx]
            s = self.get_intmpole_skew_avg(data_set, h)
            s = [s[i] for i in idx]
            mpoles_n[:, j] = n
            mpoles_s[:, j] = s

        lines = RotCoilMeas.get_excdata_text(
            pwrsupply_polarity,
            magnet_type_label,
            magnet_serial_number,
            data_set,
            main_harmonic,
            main_harmonic_type,
            harmonics,
            currents,
            mpoles_n,
            mpoles_s,
            filename)
        return lines

    def _get_data_path(self):
        mag_type_name = self.magnet_type_name
        mag_type_name = mag_type_name.replace('quadrupole', 'quadrupoles')
        mag_type_name = mag_type_name.replace('sextupole-sf', 'sextupole')
        mag_type_name = mag_type_name.replace('sextupole', 'sextupoles')
        mag_type_name = mag_type_name.replace('corrector-ch', 'correctors')
        mag_type_name = mag_type_name.replace('corrector-cv', 'correctors')
        data_path = \
            self.lnls_ima_path + '/' + mag_type_name + '/' + \
            self.model_version + '/measurement/magnetic/rotcoil/' + \
            self.family_folder + \
            self.magnet_type_label + '-' + \
            self.serial_number + '/' + self.excitation_type
        return data_path

    def _get_data_sets(self):
        data_path = self._get_data_path()
        if self.magnet_type_label == 'BQF' and self.serial_number == '053':
            fs = self._specialized_data_sets_BQF_053()
        elif self.magnet_type_label == 'Q20' and \
                self.serial_number in ('006', '093', '008', '004'):
            fs = self._specialized_data_sets_Q20()
        elif self.magnet_type_label == 'BC' and \
                self.serial_number in ('059', '060', '062'):
            fs = self._specialized_data_sets_BC()
        elif self.magnet_type_label == 'BS' and \
                self.serial_number in ('007', ):
            fs = self._specialized_data_sets_BS_007()
        elif self.magnet_type_label == 'S15' and \
                self.excitation_type == 'main_qs':
            fs = self._specialized_data_sets_S15_SKEW()
        else:
            fs = _os.listdir(data_path)
        files = []
        for f in fs:
            if _os.path.isdir(data_path + '/' + f):
                files.append(f)
        return files

    def _get_files(self, data_set):
        data_path = self._get_data_path()
        data_path += '/' + data_set
        files = _os.listdir(data_path)
        return files

    def _read_rotcoil_data(self):
        # read rotcoildata
        self._rotcoildata = dict()
        for data_set in self.data_sets:
            tstamps, mdata = [], []
            files = self.get_files(data_set)
            for file in files:
                path = self._get_data_path()
                path += '/' + data_set + '/' + file
                try:
                    meas = RotCoilData(path, self.conv_mpoles_sign)
                except Exception:
                    print('Error while trying to read {}'.format(path))
                    raise
                tstamps.append(meas.hour)
                mdata.append(meas)
            if self.magnet_type_label == 'Q14' and self.serial_number == '060':
                dataset_datum = self._specialized_sort_Q14_060(mdata)
            else:
                # sort by timestamp
                dataset_datum = [d for _, d in sorted(zip(tstamps, mdata))]
            self._rotcoildata[data_set] = dataset_datum
        # check consistency of meas data
        self._check_measdata()

    def _check_measdata(self):
        # harmonics
        self.harmonics = None
        for data_set, datum in self._rotcoildata.items():
            for d in datum:
                if self.harmonics is None:
                    self.harmonics = list(d.harmonics)
                else:
                    if d.harmonics != self.harmonics:
                        raise ValueError('Inconsistent parameter harmonics')

    def _specialized_data_sets_BQF_053(self):
        # M4 and M5 are incomplete
        return ['M1', 'M2', 'M3']

    def _specialized_data_sets_Q20(self):
        # M1_ferromag is incomplete
        return ['M1']

    def _specialized_data_sets_S15_SKEW(self):
        return ['M1']

    def _specialized_data_sets_BC(self):
        return ['M1']

    def _specialized_data_sets_BS_007(self):
        return ['M2', 'M3']

    def _specialized_sort_Q14_060(self, mdata):
        files = (
            'Q14-060_Q_BOA_000.0A_180407_095148.dat',
            'Q14-060_Q_BOA_002.0A_180407_095212.dat',
            'Q14-060_Q_BOA_004.0A_180407_095236.dat',
            'Q14-060_Q_BOA_006.0A_180407_095259.dat',
            'Q14-060_Q_BOA_008.0A_180407_095323.dat',
            'Q14-060_Q_BOA_010.0A_180407_095347.dat',
            'Q14-060_Q_BOA_030.0A_180407_095412.dat',
            'Q14-060_Q_BOA_050.0A_180407_095437.dat',
            'Q14-060_Q_BOA_070.0A_180407_095502.dat',
            'Q14-060_Q_BOA_090.0A_180407_095528.dat',
            'Q14-060_Q_BOA_110.0A_180407_095553.dat',
            'Q14-060_Q_BOA_130.0A_180407_095618.dat',
            'Q14-060_Q_BOA_148.0A_180407_100443.dat',
            'Q14-060_Q_BOA_130.0A_180407_095708.dat',
            'Q14-060_Q_BOA_110.0A_180407_095734.dat',
            'Q14-060_Q_BOA_090.0A_180407_095759.dat',
            'Q14-060_Q_BOA_070.0A_180407_095824.dat',
            'Q14-060_Q_BOA_050.0A_180407_095849.dat',
            'Q14-060_Q_BOA_030.0A_180407_095914.dat',
            'Q14-060_Q_BOA_010.0A_180407_095939.dat',
            'Q14-060_Q_BOA_008.0A_180407_100003.dat',
            'Q14-060_Q_BOA_006.0A_180407_100027.dat',
            'Q14-060_Q_BOA_004.0A_180407_100050.dat',
            'Q14-060_Q_BOA_002.0A_180407_100114.dat',
            'Q14-060_Q_BOA_000.0A_180407_100137.dat',
        )
        return self._sort(mdata, files)

    def _specialized_rampupind_Q30_011(self):
        return tuple(range(25, 37+1))

    def _specialized_rampupind_BC(self):
        return [18, 19, 20, 21, 22, 23, 24, 1, 2, 3, 4, 5, 6]

    def _specialized_rampdownind_BC(self):
        return [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]

    def _specialized_rampupind_TBQ(self):
        return [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]

    def _specialized_rampdownind_TBQ(self):
        return [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42]

    def _specialized_rampupind_TBC(self):
        return [15, 16, 17, 18, 19, 20, 1, 2, 3, 4, 5]

    def _specialized_rampdownind_TBC(self):
        return [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

    def _sort(self, mdata, files):
        dataset_datum = []
        dfiles = [d.file for d in mdata]
        for file in files:
            idx = dfiles.index(file)
            data = mdata[idx]
            dataset_datum.append(data)
        return dataset_datum


class RotCoilMeas_SI(RotCoilMeas):
    """Rotation coil measurement of SI magnets."""

    # used in case meas was taken with opposite current polarity
    conv_mpoles_sign = +1.0


class RotCoilMeas_BO(RotCoilMeas):
    """Rotation coil measurement of BO magnets."""

    # used in case meas was taken with opposite current polarity
    conv_mpoles_sign = +1.0


class RotCoilMeas_TB(RotCoilMeas):
    """Rotation coil measurement of TB magnets."""

    # used in case meas was taken with opposite current polarity
    conv_mpoles_sign = +1.0


class RotCoilMeas_SIQuadQ14(RotCoilMeas_SI, RotCoilMeas_Quad):
    """Rotation coil measurement of SI quadrupole magnets Q14."""

    conv_mpoles_sign = -1.0  # meas with opposite current polarity!
    magnet_type_label = 'Q14'
    magnet_type_name = 'si-quadrupole-q14'
    model_version = 'model-04'
    magnet_hardedge_length = 0.14  # [m]
    nominal_KL_values = {
        'SI-Fam:MA-QDA': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDA'],
        'SI-Fam:MA-QDB1': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDB1'],
        'SI-Fam:MA-QDB2': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDB2'],
        'SI-Fam:MA-QDP1': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDP1'],
        'SI-Fam:MA-QDP2': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDP2'],
    }
    spec_main_intmpole_rms_error = 0.05  # [%]
    spec_main_intmpole_max_value = 5.2116053477732  # [T] (spec in wiki-sirius)
    spec_magnetic_center_x = 40.0  # [um]
    spec_magnetic_center_y = 40.0  # [um]
    spec_roll = 0.3  # [mrad]

    spec_r0 = 12.0  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13, 17]) + 1
    spec_normal_sys_mpoles = _np.array([-3.9e-4, +1.7e-3, -8.0e-4, +8.5e-5])
    spec_normal_rms_harms = _np.array([2, 3, 4, 5]) + 1
    spec_normal_rms_mpoles = _np.array([1.5, 1.5, 1.5, 1.5])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5]) + 1
    spec_skew_rms_mpoles = _np.array([0.5, 0.5, 0.5, 0.5])*1e-4


class RotCoilMeas_SIQuadQ30(RotCoilMeas_SI, RotCoilMeas_Quad):
    """Rotation coil measurement of SI quadrupole magnets Q30."""

    conv_mpoles_sign = +1.0  # meas with default current polarity!
    magnet_type_label = 'Q30'
    magnet_type_name = 'si-quadrupole-q30'
    model_version = 'model-06'
    magnet_hardedge_length = 0.30  # [m]
    nominal_KL_values = {
        'SI-Fam:MA-QFB': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QFB'],
        'SI-Fam:MA-QFP': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QFP'],
    }
    spec_main_intmpole_rms_error = 0.05  # [%]
    spec_main_intmpole_max_value = 13.62942873208  # [T] (spec in wiki-sirius)
    spec_magnetic_center_x = 40.0  # [um]
    spec_magnetic_center_y = 40.0  # [um]
    spec_roll = 0.3  # [mrad]

    spec_r0 = 12.0  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13, 17]) + 1
    spec_normal_sys_mpoles = _np.array([-4.3e-4, +1.8e-3, -8.1e-4, +7.2e-5])
    spec_normal_rms_harms = _np.array([2, 3, 4, 5]) + 1
    spec_normal_rms_mpoles = _np.array([1.5, 1.5, 1.5, 1.5])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5]) + 1
    spec_skew_rms_mpoles = _np.array([0.5, 0.5, 0.5, 0.5])*1e-4


class RotCoilMeas_SIQuadQ20(RotCoilMeas_SI, RotCoilMeas_Quad):
    """Rotation coil measurement of SI quadrupole magnets Q20."""

    family_folder = 'family_1/'

    conv_mpoles_sign = +1.0  # meas with default current polarity!
    magnet_type_label = 'Q20'
    magnet_type_name = 'si-quadrupole-q20'
    model_version = 'model-05'
    magnet_hardedge_length = 0.20  # [m]
    nominal_KL_values = {
        'SI-Fam:MA-QFA': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QFA'],
        'SI-Fam:MA-Q1': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-Q1'],
        'SI-Fam:MA-Q2': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-Q2'],
        'SI-Fam:MA-Q3': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-Q3'],
        'SI-Fam:MA-Q4': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-Q4'],
    }
    spec_main_intmpole_rms_error = 0.05  # [%]
    spec_main_intmpole_max_value = 9.0862858213864  # [T] (spec in wiki-sirius)
    spec_magnetic_center_x = 40.0  # [um]
    spec_magnetic_center_y = 40.0  # [um]
    spec_roll = 0.3  # [mrad]

    spec_r0 = 12.0  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13, 17]) + 1
    spec_normal_sys_mpoles = _np.array([-4.1e-4, +1.7e-3, -7.7e-4, +5.9e-5])
    spec_normal_rms_harms = _np.array([2, 3, 4, 5]) + 1
    spec_normal_rms_mpoles = _np.array([1.5, 1.5, 1.5, 1.5])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5]) + 1
    spec_skew_rms_mpoles = _np.array([0.5, 0.5, 0.5, 0.5])*1e-4


class RotCoilMeas_SISextS15(RotCoilMeas_SI, RotCoilMeas_Sext):
    """Rotation coil measurement of SI sextupole S15."""

    family_folder = 'family_1/'

    conv_mpoles_sign = +1.0  # meas with default current polarity!
    magnet_type_label = 'S15'
    magnet_type_name = 'si-sextupole-s15'
    model_version = 'model-07'
    magnet_hardedge_length = 0.15  # [m]
    nominal_KL_values = {
        'SI-Fam:MA-SFA0': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFA0'],
        'SI-Fam:MA-SFB0': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFB0'],
        'SI-Fam:MA-SFP0': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFP0'],
        'SI-Fam:MA-SFA1': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFA1'],
        'SI-Fam:MA-SFB1': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFB1'],
        'SI-Fam:MA-SFP1': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFP1'],
        'SI-Fam:MA-SFA2': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFA2'],
        'SI-Fam:MA-SFB2': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFB2'],
        'SI-Fam:MA-SFP2': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SFP2'],
        # using nominal integrated strengths
        'SI-Fam:MA-SDA0': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDA0'],
        'SI-Fam:MA-SDB0': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDB0'],
        'SI-Fam:MA-SDP0': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDP0'],
        'SI-Fam:MA-SDA1': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDA1'],
        'SI-Fam:MA-SDB1': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDB1'],
        'SI-Fam:MA-SDP1': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDP1'],
        'SI-Fam:MA-SDA2': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDA2'],
        'SI-Fam:MA-SDB2': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDB2'],
        'SI-Fam:MA-SDP2': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDP2'],
        'SI-Fam:MA-SDA3': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDA3'],
        'SI-Fam:MA-SDB3': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDB3'],
        'SI-Fam:MA-SDP3': - _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-SDP3'],
    }
    spec_main_intmpole_rms_error = 0.05  # [%]
    spec_main_intmpole_max_value = 360.24921758801  # [T/m] (spec in wiki)
    spec_magnetic_center_x = 40.0  # [um]
    spec_magnetic_center_y = 40.0  # [um]
    spec_roll = 0.3  # [mrad]

    spec_r0 = 12.0  # [mm]
    spec_normal_sys_harms = _np.array([4, 6, 8, 14]) + 1
    spec_normal_sys_mpoles = _np.array([-7.0e-5, -1.4e-4, -2.4e-3, +1.4e-3])
    spec_normal_rms_harms = _np.array([3, 4, 5, 6]) + 1
    spec_normal_rms_mpoles = _np.array([7.0, 5.0, 4.0, 2.0])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([3, 4, 5, 6]) + 1
    spec_skew_rms_mpoles = _np.array([5.0, 5.0, 0.5, 0.5])*1e-4


class RotCoilMeas_BOSext(RotCoilMeas_BO, RotCoilMeas_Sext):
    """Rotation coil measurement of BO sextupoles."""

    conv_mpoles_sign = +1.0  # meas with default current polarity!
    magnet_type_label = 'BS'
    magnet_type_name = 'bo-sextupole-sf'
    model_version = 'model-03'
    magnet_hardedge_length = 0.105  # [m]
    nominal_KL_values = {
        'BO-Fam:MA-SF': _rutil.NOMINAL_STRENGTHS['BO-Fam:MA-SF'],
        'BO-Fam:MA-SD': _rutil.NOMINAL_STRENGTHS['BO-Fam:MA-SD'],
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = 21.014537692633  # [T/m] (spec wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    pwrsupply_polarity = 'bipolar'
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    spec_normal_sys_harms = _np.array([8, 14]) + 1
    spec_normal_sys_mpoles = _np.array([-2.7e-2, -1.4e-2])
    spec_normal_rms_harms = _np.array([3, 4, 5, 6, 7, 8, 9]) + 1
    spec_normal_rms_mpoles = _np.array([4, 4, 4, 4, 4, 4, 4])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([3, 4, 5, 6, 7, 8, 9]) + 1
    spec_skew_rms_mpoles = _np.array([1, 1, 1, 1, 1, 1, 1])*1e-4


class RotCoilMeas_BOQuadQD(RotCoilMeas_BO, RotCoilMeas_Quad):
    """Rotation coil measurement of BO quadrupole magnets QD."""

    conv_mpoles_sign = +1.0  # meas with default current polarity!
    magnet_type_label = 'BQD'
    magnet_type_name = 'bo-quadrupole-qd'
    model_version = 'model-02'
    magnet_hardedge_length = 0.10  # [m]
    nominal_KL_values = {
        'BO-Fam:MA-QD': _rutil.NOMINAL_STRENGTHS['BO-Fam:MA-QD'],
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = 0.52536344231582  # [T] (spec wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    pwrsupply_polarity = 'bipolar'
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13]) + 1
    spec_normal_sys_mpoles = _np.array([-4.7e-3, +1.2e-3, +5.4e-7])
    spec_normal_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_normal_rms_mpoles = _np.array([7, 4, 4, 4, 4, 4, 4])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_skew_rms_mpoles = _np.array([10, 5, 1, 1, 1, 1, 1])*1e-4


class RotCoilMeas_BOQuadQF(RotCoilMeas_BO, RotCoilMeas_Quad):
    """Rotation coil measurement of BO quadrupole magnets QF."""

    conv_mpoles_sign = +1.0  # meas with opposite current polarity!
    magnet_type_label = 'BQF'
    magnet_type_name = 'bo-quadrupole-qf'
    model_version = 'model-06'
    magnet_hardedge_length = 0.228  # [m]
    nominal_KL_values = {
        'BO-Fam:MA-QF': _rutil.NOMINAL_STRENGTHS['BO-Fam:MA-QF'],
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = 4.2554438827581  # [T] (spec wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13]) + 1
    spec_normal_sys_mpoles = _np.array([-1.0e-3, +1.1e-3, +8.0e-5])
    spec_normal_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_normal_rms_mpoles = _np.array([7, 4, 4, 4, 4, 4, 4])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_skew_rms_mpoles = _np.array([10, 5, 1, 1, 1, 1, 1])*1e-4


class RotCoilMeas_BOCorH(RotCoilMeas_BO, RotCoilMeas_HCor):
    """Rotation coil measurement of BO horizontal correctors."""

    conv_mpoles_sign = -1.0  # meas with opposite current polarity!
    magnet_type_label = 'BC'
    magnet_type_name = 'bo-corrector-ch'
    model_version = 'model-02'
    magnet_hardedge_length = 0.15018  # [m]
    nominal_KL_values = {
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = 0.003102146040341  # [T.m] (wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    # there is not multipole spec. using one based on prototype meas and
    # dynapt testes.
    spec_normal_sys_harms = _np.array([1, 2, 3, 4, 5, 6]) + 1
    spec_normal_sys_mpoles = _np.array(
        [-3.0e-4, +3.0e-3, +1.3e-4, -3.3e-3, +6.2e-4, -3.2e-3])
    spec_normal_rms_harms = _np.array([])
    spec_normal_rms_mpoles = _np.array([])
    spec_skew_sys_harms = _np.array([0, ]) + 1
    spec_skew_sys_mpoles = _np.array([0, ])
    spec_skew_rms_harms = _np.array([])
    spec_skew_rms_mpoles = _np.array([])


class RotCoilMeas_TBCorH(RotCoilMeas_TB, RotCoilMeas_HCor):
    """Rotation coil measurement of TB horizontal correctors."""

    excitation_type = 'CH'

    # conv_mpoles_sign = -1.0  # meas with opposite current polarity!
    magnet_type_label = 'TBC'
    magnet_type_name = 'tb-corrector-ch'
    model_version = 'model-03'
    magnet_hardedge_length = 0.081  # [m]
    nominal_KL_values = {
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = -0.00125  # [T.m] (wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    # there is not multipole spec. using one based on prototype meas and
    # dynapt testes.
    spec_normal_sys_harms = _np.array([0, 1, 2, 4, 6, 8]) + 1
    spec_normal_sys_mpoles = _np.array(
        [+1.0000e+00, -2.0131e-06, -3.8712e-01, -2.0729e-02,
         -1.1201e-02, -9.0639e-03])  # from model-03
    spec_normal_rms_harms = _np.array([])
    spec_normal_rms_mpoles = _np.array([])
    spec_skew_sys_harms = _np.array([0, ]) + 1
    spec_skew_sys_mpoles = _np.array([0, ])
    spec_skew_rms_harms = _np.array([])
    spec_skew_rms_mpoles = _np.array([])


class RotCoilMeas_TBCorV(RotCoilMeas_TB, RotCoilMeas_VCor):
    """Rotation coil measurement of TB vertical correctors."""

    excitation_type = 'CV'

    # conv_mpoles_sign = -1.0  # meas with opposite current polarity!
    magnet_type_label = 'TBC'
    magnet_type_name = 'tb-corrector-cv'
    model_version = 'model-03'
    magnet_hardedge_length = 0.081  # [m]
    nominal_KL_values = {
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = -0.00125  # [T.m] (wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    # there is not multipole spec. using one based on prototype meas and
    # dynapt testes.
    spec_normal_sys_harms = _np.array([0, ]) + 1
    spec_normal_sys_mpoles = _np.array([0, ])
    spec_normal_rms_harms = _np.array([])
    spec_normal_rms_mpoles = _np.array([])
    spec_skew_sys_harms = _np.array([0, 1, 2, 4, 6, 8]) + 1
    spec_skew_sys_mpoles = _np.array(
        [+1.0000e+00, -1.0154e-06, +3.8711e-01, -1.1342e-02,
         +6.6516e-03, +1.2167e-02])  # from model-03
    spec_skew_rms_harms = _np.array([])
    spec_skew_rms_mpoles = _np.array([])


class RotCoilMeas_TBQuad(RotCoilMeas_TB, RotCoilMeas_Quad):
    """Rotation coil measurement of TB quadrupole magnets."""

    conv_mpoles_sign = +1.0  # meas with default current polarity!
    magnet_type_label = 'TBQ'
    magnet_type_name = 'tb-quadrupole'
    model_version = 'model-01'
    magnet_hardedge_length = 0.10  # [m]
    nominal_KL_values = {
        'TB-01:MA-QD1': _rutil.NOMINAL_STRENGTHS['TB-01:MA-QD1'],
        'TB-01:MA-QF1': _rutil.NOMINAL_STRENGTHS['TB-01:MA-QF1'],
        'TB-02:MA-QD2A': _rutil.NOMINAL_STRENGTHS['TB-02:MA-QD2A'],
        'TB-02:MA-QF2A': _rutil.NOMINAL_STRENGTHS['TB-02:MA-QF2A'],
        'TB-02:MA-QD2B': _rutil.NOMINAL_STRENGTHS['TB-02:MA-QD2B'],
        'TB-02:MA-QF2B': _rutil.NOMINAL_STRENGTHS['TB-02:MA-QF2B'],
        'TB-03:MA-QD3': _rutil.NOMINAL_STRENGTHS['TB-03:MA-QD3'],
        'TB-03:MA-QF3': _rutil.NOMINAL_STRENGTHS['TB-03:MA-QF3'],
        'TB-04:MA-QD4': _rutil.NOMINAL_STRENGTHS['TB-04:MA-QD4'],
        'TB-04:MA-QF4': _rutil.NOMINAL_STRENGTHS['TB-04:MA-QF4'],
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = 0.8  # [T] (spec wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    pwrsupply_polarity = 'bipolar'
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13]) + 1  # from model-01
    spec_normal_sys_mpoles = _np.array([-4.7e-3, +1.2e-3, -4.0e-6])  # model-01
    # booster quads
    spec_normal_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_normal_rms_mpoles = _np.array([7, 4, 4, 4, 4, 4, 4])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_skew_rms_mpoles = _np.array([10, 5, 1, 1, 1, 1, 1])*1e-4


class MagnetsAnalysis:
    """Measurements of a magnet type magnets."""

    def __init__(self, rotcoilmeas_cls, serial_numbers):
        """Init."""
        if isinstance(serial_numbers, dict):
            # serial_numbers is a dict (for various magnet families)
            self.serials = []
            for fam, sn in serial_numbers.items():
                rotcoilmeas_cls.family_folder = fam
                for s in sn:
                    self._magnetsdata[s] = rotcoilmeas_cls(s)
                self.serials.extend(sn)
        else:
            # serial_numbers is a list (for single magnet family)
            self.serials = serial_numbers
            self._magnetsdata = dict()
            for s in serial_numbers:
                self._magnetsdata[s] = rotcoilmeas_cls(s)
        self._average = dict()

    def init(self):
        """Init."""
        # Load all data
        self.tmpl = self._magnetsdata[self.serials[0]]
        self.max_i = self.tmpl.get_max_current_index()
        self.spec_max = - self.tmpl.conv_mpoles_sign * \
            self.tmpl.spec_main_intmpole_max_value

    def print_info(self):
        """Print info."""
        if self.tmpl.conv_mpoles_sign != 1.0:
            print(('WARNING: rotating coil measurements were taken with '
                   'opposite polarity.'))
            print('')
            print(('positive currents of monopolar power supply used '
                   'generated field with opposite sign.'))
            print(('signs of all multipole values will be therefore '
                   'inverted so as to generate default'))
            print(('excitation data tables: positive currents correspond '
                   'to nominal focusing or defocusing field'))
            print('properties.')
            print('')
        fmtstr = 'index: {:02d}, serial_number: {}, data sets: {}'
        for i in range(len(self.serials)):
            sn = self.serials[i]
            print(fmtstr.format(i, sn, self._magnetsdata[sn].data_sets))

    def main_intmpole_at_max_current(self, data_set):
        """."""
        fmtstr = ('index:{:02d}, serial:{}, idx:{:02d}, max_current: '
                  '{:+10.4f} [A], diff_spec: {:+.2f} [%]')
        self.max_mpole = []
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            c = d.get_currents(data_set)
            idx = d.get_max_current_index()
            if self.tmpl.main_harmonic_type == 'normal':
                mpoles = d.get_intmpole_normal_avg(data_set, d.main_harmonic)
            else:
                mpoles = d.get_intmpole_skew_avg(data_set, d.main_harmonic)
            self.max_mpole.append(mpoles[idx])
            diff_spec = 100*(mpoles[idx] - self.spec_max)/self.spec_max
            print(fmtstr.format(i, d.serial_number, idx, c[idx], diff_spec))

    def main_intmpole_at_max_current_plot(self, plt):
        """."""
        y = (self.spec_max, ) * 2
        plt.plot([0, len(self.max_mpole)-1], y, '--k')
        plt.plot(self.max_mpole, 'og')
        plt.grid(True)
        plt.legend(('Spec', 'Data'))
        plt.xlabel('Serial Number Index')
        if isinstance(self.tmpl, RotCoilMeas_Quad):
            plt.ylabel('Integrated Quadrupole [T]')
            plt.title(('Comparison of Integrated Quadrupole at Maximum '
                       'Current x Specification'))
        elif isinstance(self.tmpl, RotCoilMeas_Sext):
            plt.ylabel('Integrated Sextupole [T/m]')
            plt.title(('Comparison of Integrated Sextupole at Maximum '
                       'Current x Specification'))
        elif isinstance(self.tmpl, RotCoilMeas_Cor):
            plt.ylabel('Integrated Dipole [T.m]')
            plt.title(('Comparison of Integrated Dipole at Maximum '
                       'Current x Specification'))
        else:
            raise NotImplementedError

    def magnetic_center_direction_plot(self, data_set, direction, plt):
        """."""
        v = []
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            c = d.get_currents(data_set)
            if direction in ('x', 'X', 'h', 'H'):
                u = d.get_magnetic_center_x(data_set)
                plt.plot(c, u, 'b')
            elif direction in ('y', 'Y', 'v', 'V'):
                u = d.get_magnetic_center_y(data_set)
                plt.plot(c, u, 'r')
            else:
                raise NotImplementedError()
            idx = d.get_max_current_index()
            v.append(u[idx])

        v = _np.array(v)
        if direction in ('x', 'X', 'h', 'H'):
            dstr = 'Horizontal '
            specp = (+self.tmpl.spec_magnetic_center_x, ) * 2
            specn = (-self.tmpl.spec_magnetic_center_x, ) * 2
        elif direction in ('y', 'Y', 'v', 'V'):
            dstr = 'Vertical '
            specp = (+self.tmpl.spec_magnetic_center_y, ) * 2
            specn = (-self.tmpl.spec_magnetic_center_y, ) * 2

        fmtstr = dstr + 'center at maximum current [um]: {:+.2f} ± {:.2f}'
        print(fmtstr.format(_np.mean(v), _np.std(v)))

        plt.plot([min(c), max(c)], specp, '--k')
        plt.plot([min(c), max(c)], specn, '--k')
        plt.xlabel('Current [A]')
        plt.ylabel(dstr + 'position [um]')
        plt.title(dstr + 'center of magnets fields x current')
        plt.grid(True)

    def magnetic_center_plot(self, data_set, plt):
        """."""
        xv, yv = [], []
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            idx = d.get_max_current_index()
            x = d.get_magnetic_center_x(data_set)
            y = d.get_magnetic_center_y(data_set)
            xv.append(x[idx])
            yv.append(y[idx])
        specp = (+self.tmpl.spec_magnetic_center_x, ) * 2
        specn = (-self.tmpl.spec_magnetic_center_x, ) * 2
        plt.plot([0, len(self)-1], specp, '--k')
        plt.plot([0, len(self)-1], specn, '--k')
        plt.plot(xv, 'ob')
        plt.plot(yv, 'or')
        plt.xlabel('Serial Number Index')
        plt.ylabel('Position [um]')
        plt.legend(('Spec', 'Spec', 'X', 'Y'))
        plt.title('Magnetic Centers at Maximum Current')

    def magnetic_center_transverse_plot(self, data_set, plt):
        """."""
        xv, yv = [], []
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            idx = d.get_max_current_index()
            x = d.get_magnetic_center_x(data_set)
            y = d.get_magnetic_center_y(data_set)
            xv.append(x[idx])
            yv.append(y[idx])

        # plot
        sx = self.tmpl.spec_magnetic_center_x
        sy = self.tmpl.spec_magnetic_center_y
        plt.plot([-sx, -sx], [-sy, sy], '--k')
        plt.plot([-sx, sx], [sy, sy], '--k')
        plt.plot([sx, sx], [sy, -sy], '--k')
        plt.plot([sx, -sx], [-sy, -sy], '--k')
        for x, y in zip(xv, yv):
            plt.plot([x], [y], 'o', color=[1, 0, 1])
        plt.xlabel('Horizontal Position [um]')
        plt.ylabel('Vertical Position [um]')
        plt.title('Magnetic Centers at Maximum Current')

    def rotation_error_plot(self, data_set, plt, curr_idx):
        """."""
        rot_error = []
        ind = self.tmpl.get_rampup_indices()
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            spec_roll = d.spec_roll
            n = _np.array(d.get_intmpole_normal_avg(data_set, d.main_harmonic))
            s = _np.array(d.get_intmpole_skew_avg(data_set, d.main_harmonic))
            n, s = n[ind], s[ind]
            if self.tmpl.main_harmonic_type == 'normal':
                r = s[curr_idx]/n[curr_idx]
            else:
                r = n[curr_idx]/s[curr_idx]
            theta = _np.arctan(r)/d.main_harmonic
            rot_error.append(theta)
        rot_error = 1000*_np.array(rot_error)
        avg = _np.mean(rot_error)
        std = _np.std(rot_error)
        if plt is not None:
            plt.plot([0, len(rot_error)-1], [+spec_roll, ]*2, 'k--')
            plt.plot([0, len(rot_error)-1], [-spec_roll, ]*2, 'k--')
            plt.plot([0, len(rot_error)-1], [avg, avg], 'r-')
            plt.plot([0, len(rot_error)-1], [avg-std, avg-std], 'r--')
            plt.plot([0, len(rot_error)-1], [avg+std, avg+std], 'r--')
            plt.plot(rot_error, 'ro')
            plt.xlabel('Magnet Index')
            plt.ylabel('Rotation Error [mrad]')
            plt.title(('Main Multipole Rotation Error '
                       '(Current Index {})').format(curr_idx))
        return spec_roll, avg, std

    def rotation_error_vs_current_plot(self, data_set, energy, plt):
        """."""
        vec_curr, *_ = self.tmpl.get_rampup(data_set)
        vec_curr = _np.array(vec_curr)
        vec_avg = _np.zeros(len(vec_curr))
        vec_std = _np.zeros(len(vec_curr))
        for i in range(len(vec_curr)):
            spec, vec_avg[i], vec_std[i] = \
                self.rotation_error_plot(data_set, None, i)
        d = self.tmpl.get_nominal_main_intmpole_values(energy)
        nom_curr, nom_fams = [], []
        for fam, gl in d.items():
            c = self.tmpl.rampup_main_mpole_2_curr(data_set, gl)
            nom_curr.append(c)
            nom_fams.append(fam)
        nom_curr = _np.array(nom_curr)
        nom_fams = _np.array(nom_fams)
        arr1inds = nom_curr.argsort()
        nom_curr = nom_curr[arr1inds]
        nom_fams = nom_fams[arr1inds]
        plt.plot([min(vec_curr), max(vec_curr)], [+spec, ]*2, 'k--')
        plt.plot([min(vec_curr), max(vec_curr)], [-spec, ]*2, 'k--')
        plt.plot(vec_curr, vec_avg, '-r')
        plt.plot(vec_curr, vec_avg - vec_std, '--r')
        plt.plot(vec_curr, vec_avg + vec_std, '--r')
        y = [min(-spec, min(vec_avg-vec_std)),
             max(+spec, max(vec_avg+vec_std))]
        print('Currents for nominal strengths:')
        for i in range(len(nom_curr)):
            plt.plot([nom_curr[i], ]*2, y, '-.g')
            print('{:<10s}: {:.1f} A'.format(nom_fams[i], nom_curr[i]))
        plt.xlabel('RampUp Current [A]')
        plt.ylabel('Rotation Error [mrad]')
        plt.title(('Rotation Error as Defined by Skew and Normal '
                   'of Main Multipole'))
        plt.grid(True)

    def rampup_excitation_curve_plot(self, data_set, energy, plt):
        """."""
        c_min, c_max = 1, 1
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            c, gl = d.get_rampup(data_set)
            c_min, c_max = min(c_min, min(c)), max(c_max, max(c))
            plt.plot(c, gl, 'og')

        y = (self.spec_max, ) * 2
        plt.plot([c_min, c_max], y, '--k')

        if isinstance(self.tmpl, RotCoilMeas_Quad):
            sstr = 'Integrated Quadrupole [T]'
        elif isinstance(self.tmpl, RotCoilMeas_Sext):
            sstr = 'Integrated Sextupole [T/m]'
        elif isinstance(self.tmpl, RotCoilMeas_Cor):
            sstr = 'Integrated Dipole [T.m]'
        else:
            raise NotImplementedError()

        print('Nominal ' + sstr + ':')
        nom = self.tmpl.get_nominal_main_intmpole_values(energy)
        for fam, v in nom.items():
            print('{:<16s}: {:+.6f}'.format(fam, v))
            plt.plot([c_min, c_max], [v, v], '--', color=[0, 0.5, 0])

        plt.xlabel('Current [A]')
        plt.ylabel(sstr)
        plt.title('Ramp Up of All Magnets')
        plt.grid(True)

    def rampup_excitation_curve_dispersion_plot(self, data_set, plt):
        """."""
        ind = self.tmpl.get_rampup_indices()
        shape = (len(self.serials), len(ind))
        # shape = (len(self.serials), 1+self.tmpl.get_max_current_index())
        print(shape)
        c, g = _np.zeros(shape), _np.zeros(shape)
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            ct, gt = d.get_rampup(data_set)
            c[i, :] = ct
            g[i, :] = gt

        c_avg = _np.mean(c, axis=0)
        g_avg = _np.mean(g, axis=0)
        g_std = _np.std(g, axis=0)

        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            gl_interp = d.rampup_curr_2_main_mpole(data_set, c_avg)
            g_dif = gl_interp - g_avg
            # print(d.serial_number, max(abs(g_dif)))
            plt.plot(c_avg, g_dif)
        plt.plot(c_avg, +g_std, '--k', linewidth=4)
        plt.plot(c_avg, -g_std, '--k', linewidth=4)

        if isinstance(self.tmpl, RotCoilMeas_Quad):
            sstr = 'Integrated Quadrupole'
        elif isinstance(self.tmpl, RotCoilMeas_Sext):
            sstr = 'Integrated Sextupole'
        elif isinstance(self.tmpl, RotCoilMeas_Cor):
            sstr = 'Integrated Dipole'
        else:
            raise NotImplementedError()

        plt.xlabel('Current [A]')
        plt.ylabel(sstr + ' [T]')
        plt.title('Difference of Magnets ' + sstr + ' from Average')

    def rampup_excitation_curve_rms_error_print(self, data_set):
        """."""
        def get_gl_set(current_index):
            c, g = [], []
            for d in self._magnetsdata.values():
                ct, gt = d.get_rampup(data_set)
                c.append(ct[current_index])
                g.append(gt[current_index])
            g_avg = _np.mean(g)
            g_std = _np.std(g)
            return g_avg, g_std, c, g

        currents, _ = self.tmpl.get_rampup(data_set)
        errors, cs, gs = [], [], []
        for i in range(len(currents)):
            g_avg, g_std, c, g = get_gl_set(i)
            error = [100*(gv - g_avg)/g_avg for gv in g]
            fmtstr = ('current {:02d}: {:+8.3f} [A], rms_error: {:7.4f} [%], '
                      'max_error: {:7.4f} [%]')
            print(fmtstr.format(i, _np.mean(c),
                                abs(100*g_std/g_avg), max(_np.abs(error))))
            errors.append(error)
            cs.append(c)
            gs.append(g)
        self.errors = errors

    def rampup_excitation_curve_rms_error_plot(self, plt):
        """."""
        dat = self.errors[-1]
        spec_rms = self.tmpl.spec_main_intmpole_rms_error
        # avg, std = _np.mean(dat), _np.std(dat)
        plt.plot(dat, 'og')
        plt.plot((1, len(dat)), (spec_rms, spec_rms), '--k')
        plt.plot((1, len(dat)), (-spec_rms, -spec_rms), '--k')
        plt.title('Magnets Integrated Main Multipole at Maximum Current')
        plt.xlabel('Serial Number Index')
        plt.ylabel('Difference from average [%]')

    def hysteresis_absolute_plot(self, data_set, plt):
        """."""
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            gl, c, h, area = d.get_rampdown_hysteresis(data_set)
            plt.plot(c, h)
        plt.title('Absolute Rampdown Hysteresis of All Magnets')
        plt.xlabel('Current [A]')
        if isinstance(self.tmpl, RotCoilMeas_Quad):
            plt.ylabel('Quadrupole Hysteresis [T]')
        elif isinstance(self.tmpl, RotCoilMeas_Sext):
            plt.ylabel('Sextupole Hysteresis [T/m]')
        elif isinstance(self.tmpl, RotCoilMeas_Cor):
            plt.ylabel('Sextupole Hysteresis [T.m]')
        else:
            raise NotImplementedError()
        plt.grid(True)

    def hysteresis_relative_plot(self, data_set, plt):
        """."""
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            gl, c, h, area = d.get_rampdown_hysteresis(data_set)
            r = [100*h[i]/gl[i] for i in range(len(h))]
            plt.plot(c, r)
            plt.title('Relative Rampdown Hysteresis of All Magnets')
            plt.xlabel('Current [A]')
            if isinstance(self.tmpl, RotCoilMeas_Quad):
                plt.ylabel('Quadrupole Hysteresis [%]')
            elif isinstance(self.tmpl, RotCoilMeas_Sext):
                plt.ylabel('Sextupole Hysteresis [%]')
            elif isinstance(self.tmpl, RotCoilMeas_Cor):
                plt.ylabel('Dipole Hysteresis [%]')
            else:
                raise NotImplementedError()
            plt.grid(True)

    # def save_excdata_average(self, data_set, harmonics=None):
    #     """Save excitation data."""
    #     pwrsupply_polarity = self.tmpl.pwrsupply_polarity
    #     main_harmonic = int(self.tmpl.main_harmonic)
    #     main_harmonic_type = self.tmpl.main_harmonic_type
    #     if harmonics is None:
    #         harmonics = sorted([int(h) for h in self.tmpl.harmonics])
    #     magnet_type_label = self.tmpl.magnet_type_label
    #     magnet_serial_number = None
    #     filename = self.tmpl.magnet_type_name + '-fam'
    #     units = ''
    #     for h in harmonics:
    #         unit = _mutil.get_multipole_si_units(h-1)
    #         units += unit + ' ' + unit + '  '
    #     units = units.strip()
    #
    #     if data_set not in self._average:
    #         self._create_average_mpoles(data_set)
    #
    #     # build text lines
    #     lines = RotCoilMeas.get_excdata_text(
    #         pwrsupply_polarity,
    #         magnet_type_label,
    #         magnet_serial_number,
    #         data_set,
    #         main_harmonic,
    #         main_harmonic_type,
    #         harmonics,
    #         self._average[data_set]['currents'],
    #         self._average[data_set]['mpoles_n'],
    #         self._average[data_set]['mpoles_s'],
    #         filename)
    #
    #     # save data to file
    #     with open(filename + '.txt', 'w') as fp:
    #         for line in lines:
    #             fp.write(line + '\n')

    def conv_current_2_mpoles(self, data_set, current):
        """."""
        if data_set not in self._average:
            self._create_average_mpoles(data_set)
        currents = self._average[data_set]['currents']
        mpoles_n = self._average[data_set]['mpoles_n']
        mpoles_s = self._average[data_set]['mpoles_s']
        interp_mpoles_n = _np.zeros(mpoles_n.shape)
        interp_mpoles_s = _np.zeros(mpoles_s.shape)
        return None

    # def save_excdata_individuals(self, data_set, harmonics=None):
    #     """Save excdata of all individual magnets."""
    #     for data in self._magnetsdata.values():
    #         data.save_excdata(data_set, harmonics)

    def multipole_errors_kickx_plot(self, data_set, plt, curr_idx=None,
                                    excluded_monomials_plot1=None,
                                    excluded_monomials_plot2=None,
                                    energy=3.0):
        """."""
        excmon1, excmon2 = None, None
        nr_plots = 0
        if excluded_monomials_plot1 is not None:
            nr_plots += 1
            excmon1 = excluded_monomials_plot1
        if excluded_monomials_plot2 is not None:
            nr_plots += 1
            if excmon1 is None:
                excmon1 = excluded_monomials_plot2
            else:
                excmon2 = excluded_monomials_plot2

        if nr_plots == 0:
            print('Invalid excluded monomials argument!')
            return

        # convert rampup index to general index
        ind = self.tmpl.get_rampup_indices()
        idx = ind[curr_idx]

        gs = _gridspec.GridSpec(1, nr_plots)
        gs.update(left=0.10, right=0.98, hspace=0, wspace=0.30, top=0.98)
        ax1 = plt.subplot(gs[0, 0])
        if excmon2 is not None:
            ax2 = plt.subplot(gs[0, 1], sharex=ax1, sharey=ax1)
        for s in self._magnetsdata:
            mdata = self._magnetsdata[s]
            if excmon1 is not None:
                x, y, kx, ky = mdata.multipoles_kicks_residual(
                    data_set, idx,
                    excluded_monomials_norm=excmon1,
                    excluded_monomials_skew=excmon1,
                    energy=energy)
                ax1.plot(1e3*x, 1e6*kx, color=[0, 0, 0.5])
            if excmon2 is not None:
                x, y, kx, ky = mdata.multipoles_kicks_residual(
                    data_set, idx,
                    excluded_monomials_norm=excmon2,
                    excluded_monomials_skew=excmon2,
                    energy=energy)
                ax2.plot(1e3*x, 1e6*kx, color=[0, 0, 0.5])
        x, y, kx, ky = self.tmpl.multipoles_kicks_spec_sys(
            data_set, idx, energy)
        if excmon1 is not None:
            ax1.plot(1e3*x, 1e6*kx, 'r')
        if excmon2 is not None:
            ax2.plot(1e3*x, 1e6*kx, 'r')
        x, y, kx, ky = self.tmpl.multipoles_kicks_spec_rms(
            data_set, idx, energy)
        if excmon1 is not None:
            ax1.plot(1e3*x, 1e6*kx[0], '--r')
            ax1.plot(1e3*x, 1e6*kx[1], '--r')
            ax1.set_xlabel('X [mm]')
            ax1.set_ylabel('Residual kick [urad]')
            ax1.grid()
        if excmon2 is not None:
            ax2.plot(1e3*x, 1e6*kx[0], '--r')
            ax2.plot(1e3*x, 1e6*kx[1], '--r')
            ax2.set_xlabel('X [mm]')
            ax2.set_ylabel('Residual kick [urad]')
            ax2.grid()

    def multipole_errors_kicky_plot(self, data_set, plt, curr_idx=None,
                                    excluded_monomials_plot1=None,
                                    excluded_monomials_plot2=None,
                                    energy=3.0):
        """."""
        excmon1, excmon2 = None, None
        nr_plots = 0
        if excluded_monomials_plot1 is not None:
            nr_plots += 1
            excmon1 = excluded_monomials_plot1
        if excluded_monomials_plot2 is not None:
            nr_plots += 1
            if excmon1 is None:
                excmon1 = excluded_monomials_plot2
            else:
                excmon2 = excluded_monomials_plot2

        if nr_plots == 0:
            print('Invalid excluded monomials argument!')
            return

        # convert rampup index to general index
        ind = self.tmpl.get_rampup_indices()
        idx = ind[curr_idx]

        gs = _gridspec.GridSpec(1, nr_plots)
        gs.update(left=0.10, right=0.98, hspace=0, wspace=0.30, top=0.98)
        ax1 = plt.subplot(gs[0, 0])
        if excmon2 is not None:
            ax2 = plt.subplot(gs[0, 1], sharex=ax1, sharey=ax1)
        for s in self._magnetsdata:
            mdata = self._magnetsdata[s]
            if excmon1 is not None:
                x, y, kx, ky = mdata.multipoles_kicks_residual(
                    data_set, idx,
                    excluded_monomials_norm=excmon1,
                    excluded_monomials_skew=excmon1,
                    energy=energy)
                ax1.plot(1e3*x, 1e6*ky, color=[0, 0, 0.5])
            if excmon2 is not None:
                x, y, kx, ky = mdata.multipoles_kicks_residual(
                    data_set, idx,
                    excluded_monomials_norm=excmon2,
                    excluded_monomials_skew=excmon2,
                    energy=energy)
                ax2.plot(1e3*x, 1e6*ky, color=[0, 0, 0.5])
        x, y, kx, ky = self.tmpl.multipoles_kicks_spec_sys(
            data_set, idx, energy)
        if excmon1 is not None:
            ax1.plot(1e3*x, 1e6*ky, 'r')
        if excmon2 is not None:
            ax2.plot(1e3*x, 1e6*ky, 'r')
        x, y, kx, ky = self.tmpl.multipoles_kicks_spec_rms(
            data_set, idx, energy)
        if excmon1 is not None:
            ax1.plot(1e3*x, 1e6*ky[0], '--r')
            ax1.plot(1e3*x, 1e6*ky[1], '--r')
            ax1.set_xlabel('X [mm]')
            ax1.set_ylabel('Residual kick [urad]')
            ax1.grid()
        if excmon2 is not None:
            ax2.plot(1e3*x, 1e6*ky[0], '--r')
            ax2.plot(1e3*x, 1e6*ky[1], '--r')
            ax2.set_xlabel('X [mm]')
            ax2.set_ylabel('Residual kick [urad]')
            ax2.grid()

    def readme_print(self, data_set, curr_idx):
        """."""
        c, gl = self.tmpl.get_rampup(data_set)
        cidx = curr_idx
        line1 = self.tmpl.magnet_type_label + \
            ' Magnetic Center and Integrated Main Multipole'
        print(line1)
        print(len(line1)*'='+'\n')
        print('As measured with rotcoil for I = {0:3.0f}A\n'.format(c[cidx]))
        print('{0:7s} |{1:^29s} |'.format('Magnet', 'M1'))
        print('{0:7s} |'.format(''), end='')
        if isinstance(self.tmpl, RotCoilMeas_Quad):
            st = '{:>8s} {:>8s} {:>10s} |'.format(
                'x0 [mm]', 'y0 [mm]', 'GL/I [T/mA]')
        elif isinstance(self.tmpl, RotCoilMeas_Sext):
            st = '{:>8s} {:>8s} {:>10s} |'.format(
                'x0 [mm]', 'y0 [mm]', 'SL/I [T/m/mA]')
        elif isinstance(self.tmpl, RotCoilMeas_Cor):
            st = '{:>8s} {:>8s} {:>10s} |'.format(
                'x0 [mm]', 'y0 [mm]', 'BL/I [T.m/mA]')
        print(1*st)
        for s in self._magnetsdata:
            print('{}-{} |'.format(self.tmpl.magnet_type_label, s), end='')
            for med in self.tmpl.data_sets:
                x = self._magnetsdata[s].get_magnetic_center_x(med)[cidx]
                y = self._magnetsdata[s].get_magnetic_center_y(med)[cidx]
                c, gl = self._magnetsdata[s].get_rampup(med)
                print('{:+8.1f} {:+8.1f} {:+10.4f}  |'.format(
                    x, y, 1000*gl[cidx]/c[cidx]), end='')
            print()

    def readme_multipoles_print(self, data_set, curr_idx):
        """."""
        # print comment lines
        print(('# multipoles are divided by excitation current and units are'
               ' [T], [m] and [A]'))
        print('# harmonics (dipole n=0): ', end='')
        for h in self.tmpl.harmonics:
            print('{:02d} '.format(h-1), end='')
        print('')
        print(('# MAG_LABEL   CURRENT[A]   '
               '<NORMAL_MULTIPOLES/CURRENT>[T/m^(n-1)/A]   '
               '<SKEW_MULTIPOLES/CURRENT>[T/m^(n-1)/A]'))

        # print multipoles
        for mag, mdata in self._magnetsdata.items():
            print('{:9<s}   '.format(self.tmpl.magnet_type_label + '-' + mag),
                  end='')
            c, _ = mdata.get_rampup(data_set)
            current = c[curr_idx]
            if current == 0:
                current = 1e-4  # eventual to avoid division by zero
            print('{:+.6e}   '.format(current), end='')
            ind = mdata.get_rampup_indices()
            idx = ind[curr_idx]
            nmpoles = mdata.get_intmpole_normal_avg_current(data_set, idx)
            nmpoles = _np.array(nmpoles) / current
            for i in range(len(nmpoles)):
                print('{:+.6e} '.format(nmpoles[i]), end='')
            print('  ', end='')
            smpoles = mdata.get_intmpole_skew_avg_current(data_set, idx)
            smpoles = _np.array(smpoles) / current
            for i in range(len(smpoles)):
                print('{:+.6e} '.format(smpoles[i]), end='')
            print('  ', end='')
            print('')

    def _create_average_mpoles(self, data_set):
        # average current
        currents = list()
        for data in self._magnetsdata.values():
            c, _ = data.get_rampup(data_set)
            currents.append(c)
        currents = _np.mean(_np.array(currents), axis=0)

        # calc average integrated multipoles
        shape = (len(currents), len(self.tmpl.harmonics))
        mpoles_n = _np.zeros(shape)
        mpoles_s = _np.zeros(shape)
        idx = self.tmpl.get_rampup_indices()
        for j in range(len(self.tmpl.harmonics)):
            h = self.tmpl.harmonics[j]
            for data in self._magnetsdata.values():
                n = data.get_intmpole_normal_avg(data_set, h)
                n = [n[i] for i in idx]
                s = data.get_intmpole_skew_avg(data_set, h)
                s = [s[i] for i in idx]
                mpoles_n[:, j] += n
                mpoles_s[:, j] += s
            mpoles_n[:, j] /= len(self._magnetsdata)
            mpoles_s[:, j] /= len(self._magnetsdata)
            self._average[data_set] = {
                'currents': currents,
                'mpoles_n': mpoles_n,
                'mpoles_s': mpoles_s}

    def __getitem__(self, key):
        """Return magnet data."""
        return self._magnetsdata[key]

    def __len__(self):
        """Return number of magnets."""
        return len(self._magnetsdata)

    def __iter__(self):
        """Iter."""
        return iter(self._magnetsdata)
