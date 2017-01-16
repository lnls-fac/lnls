#!/usr/bin/env python3

#import lnls

#import shutil
import os as _os
import numpy as _numpy
import math as _math
import mathphys as _mp
import matplotlib.pyplot as _plt
import matplotlib.gridspec as _gridspec
import matplotlib.ticker as _ticker
import os as _os
import lnls.utils as _utils
#from . import utils as _utils

_colors = ['blue', 'red', 'green', 'orange', 'black']
_labels = {1:'Integrated dipole [T.m]',
          2:'Integrated quadrupole [T]',
          3:'Integrated sextupole [T/m]'}

class AnalysisParameters():
    def __init(self):
        pass
    def __str__(self):
        r = ''
        r += '{0:<30s} {1:s}'.format('label', self.label)
        r += '\n{0:<30s} {1:s}'.format('main_harmonic', '{0:d} ({1:s})'.format(self.main_multipole_harmonic, get_harmonic_label(self.main_multipole_harmonic)))
        r += '\n{0:<30s} {1:s}'.format('main_harmonic_is_skew', str(self.main_multipole_is_skew))
        r += '\n{0:<30s} {1:f}'.format('reference_radius[mm]', 1000*self.ref_radius)
        r += '\n{0:<30s} {1:s}'.format('harmonics', str(self.harmonics))
        return r

class BSAnalysisParameters(AnalysisParameters):

    def __init__(self):
        self.label = 'Parameters for booster sextupoles'
        self.main_multipole_harmonic = 3  # [1: dipole, 2:quadrupole, ...]
        self.main_multipole_is_skew  = False
        self.harmonics = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] # [1: dipole, 2:quadrupole, ...]
        self.ref_radius = 0.0175 # [m]
        self.x_misalignment_spec = 160     #[um]
        self.y_misalignment_spec = 160     #[um]
        self.roll_rotation_spec  = 0.8     #[mrad]
        self.max_integ_mult_spec = -21.05  #[T/m]
        self.excitation_rms_spec = 0.3     #[%]
        self.multipoles_spec = {
            # (normal_sys, skew_sys) (normal_std, skew_std)
            4: ((0,0),       (+4.0e-4,1e-4)),
            5: ((0,0),       (+4.0e-4,1e-4)),
            6: ((0,0),       (+4.0e-4,1e-4)),
            7: ((0,0),       (+4.0e-4,1e-4)),
            8: ((0,0),       (+4.0e-4,1e-4)),
            9: ((-2.5e-2,0), (+4.0e-4,1e-4)),
            10:((0,0),       (+4.0e-4,1e-4)),
            15:((-1.5e-2,0), (0,0)),
        }

class BQAnalysisParameters(AnalysisParameters):

    def __init__(self):
        self.label = 'Parameters for booster quadrupoles'
        self.main_multipole_harmonic = 2  # [1: dipole, 2:quadrupole, ...]
        self.main_multipole_is_skew  = False
        self.harmonics = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] # [1: dipole, 2:quadrupole, ...]
        self.ref_radius = 0.0175 # [m]
        self.x_misalignment_spec = 160     #[um]
        self.y_misalignment_spec = 160     #[um]
        self.roll_rotation_spec  = 0.8     #[mrad]
        self.max_integ_mult_spec = -4.255  #[T]
        self.excitation_rms_spec = 0.3     #[%]
        self.multipoles_spec = {
            # (normal_sys, skew_sys) (normal_std, skew_std)
            3: ((0,0),       (+7.0e-4,1e-3)),
            4: ((0,0),       (+4.0e-4,5e-4)),
            5: ((0,0),       (+4.0e-4,1e-4)),
            6: ((-1.0e-3,0), (+4.0e-4,1e-4)),
            7: ((0,0),       (+4.0e-4,1e-4)),
            8: ((0,0),       (+4.0e-4,1e-4)),
            9: ((0,0),       (+4.0e-4,1e-4)),
            10:((+1.1e-3,0), (0,0)),
            14:((+8.0e-5,0), (0,0)),
        }

class SI_Q14_AnalysisParameters(AnalysisParameters):

    def __init__(self):
        self.label = 'Parameters for SI Q14 quadrupoles'
        self.main_multipole_harmonic = 2  # [1: dipole, 2:quadrupole, ...]
        self.main_multipole_is_skew  = False
        self.harmonics = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18] # [1: dipole, 2:quadrupole, ...]
        self.ref_radius = 0.012 # [m]
        self.x_misalignment_spec = 40      #[um]
        self.y_misalignment_spec = 40      #[um]
        self.roll_rotation_spec  = 0.3     #[mrad]
        self.max_integ_mult_spec = 5.2116  #[T]
        self.excitation_rms_spec = 0.05    #[%]
        self.multipoles_spec = {
            # (normal_sys, skew_sys) (normal_std, skew_std)
            3: ((0,0),               (1.5e-4,0.5e-4)),
            4: ((0,0),               (1.5e-4,0.5e-4)),
            5: ((0,0),               (1.5e-4,0.5e-4)),
            6: ((-3.9e-4,0),         (+1.5e-4,0.5e-4)),
            10:((+1.7e-3,0),         (0,0)),
            14:((-8.0e-4,0),         (0,0)),
            18:((+8.5e-5,0),         (0,0)),
        }

class SI_Q20_AnalysisParameters(AnalysisParameters):

    def __init__(self):
        self.label = 'Parameters for SI Q20 quadrupoles'
        self.main_multipole_harmonic = 2  # [1: dipole, 2:quadrupole, ...]
        self.main_multipole_is_skew  = False
        self.harmonics = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18] # [1: dipole, 2:quadrupole, ...]
        self.ref_radius = 0.012 # [m]
        self.x_misalignment_spec = 40      #[um]
        self.y_misalignment_spec = 40      #[um]
        self.roll_rotation_spec  = 0.3     #[mrad]
        self.max_integ_mult_spec = 9.08629 #[T]
        self.excitation_rms_spec = 0.05    #[%]
        self.multipoles_spec = {
            # (normal_sys, skew_sys) (normal_std, skew_std)
            3: ((0,0),               (1.5e-4,0.5e-4)),
            4: ((0,0),               (1.5e-4,0.5e-4)),
            5: ((0,0),               (1.5e-4,0.5e-4)),
            6: ((-4.1e-4,0),         (1.5e-4,0.5e-4)),
            10:((+1.7e-3,0),         (0,0)),
            14:((-7.7e-4,0),         (0,0)),
            18:((+5.9e-5,0),         (0,0)),
        }

class RadialRotatingCoil:
    def __init__(self, lines):
        for line in lines:
            words = line.strip().split()
            if not words: continue
            first = words[0]
            if first[0] == '#':
                if len(words)>3 and words[3] == 'Armazenados(V.s)':
                    self.data_conversion_factor = float(words[4].replace('[','').replace(']',''))
                continue
            elif first == 'nome_bobina_girante':
                self.label = ' '.join(words[1:])
            elif first == 'tipo_bobina_girante':
                self.type = ' '.join(words[1:])
            elif first == 'velocidade(rps)':
                self.rotation_velocity = float(words[1])
            elif first == 'aceleracao(rps^2)':
                self.rotation_acceleration = float(words[1])
            elif first == 'sentido_de_rotacao':
                self.rotation_wise = ' '.join(words[1:])
            elif first == 'ganho_integrador':
                self.integrator_gain = float(words[1])
            elif first == 'n_espiras_bobina_principal':
                self.nr_coil_turns = int(words[1])
            elif first == 'raio_interno_bobina_princip(m)':
                self.inner_radius = float(words[1])
            elif first == 'raio_externo_bobina_princip(m)':
                self.outer_radius = float(words[1])
            elif first == 'tipo_medicao':
                self.measurement_type = ' '.join(words[1:])
            elif first == 'pulso_start_coleta':
                self.init_tick = int(words[1])
            elif first == 'nr_pontos_integracao':
                self.nr_points = int(words[1])
            elif first == 'n_espiras_bobina_bucked':
                self.nr_bucked_coil_turns_bucked = int(words[1])
            elif first == 'raio_interno_bobina_bucked(m)':
                self.bucked_inner_radius = float(words[1])
            elif first == 'raio_externo_bobina_bucked(m)':
                self.bucked_outer_radius = float(words[1])
            else:
                pass

    def __str__(self):
        r = ''
        r +=   '{0:<30s} {1:s}'.format('label', self.label)
        r += '\n{0:<30s} {1:s}'.format('type', self.type)
        r += '\n{0:<30s} {1:s}'.format('rotation_direction', self.rotation_wise)
        r += '\n{0:<30s} {1:f}'.format('rotation_velocity[rps]', self.rotation_velocity)
        r += '\n{0:<30s} {1:f}'.format('rotation_acceleration[rps^2]', self.rotation_acceleration)
        r += '\n{0:<30s} {1:f}'.format('inner_radius[m]', self.inner_radius)
        r += '\n{0:<30s} {1:f}'.format('outer_radius[m]', self.outer_radius)
        r += '\n{0:<30s} {1:f}'.format('integrator_gain', self.integrator_gain)
        r += '\n{0:<30s} {1:d}'.format('nr_coil_turns', self.nr_coil_turns)
        r += '\n{0:<30s} {1:s}'.format('measurement_type', self.measurement_type)
        r += '\n{0:<30s} {1:d}'.format('init_tick', self.init_tick)
        r += '\n{0:<30s} {1:d}'.format('nr_points_integrator', self.nr_points)
        return r

class Measurement:

    def __init__(self, filename):
        lines = [line.strip() for line in open(filename, encoding='latin-1')]
        self.rotating_coil = RadialRotatingCoil(lines)
        self.filename = filename
        self.timestamp = ''
        for i in range(len(lines)):
            line = lines[i]
            if not line: continue
            if 'Volta_1' in line:
                raw_data_txt = lines[i+1:]
                break
            if line[0] == '#': continue
            words = line.split()
            first = words[0]
            if first == 'nr_voltas':
                self.nr_turns = int(words[1])
            elif first == 'nr_pontos_integracao':
                self.nr_points = int(words[1])
            elif first == 'corrente_alim_principal_avg(A)':
                self.current1_avg = float(words[1])
            elif first == 'corrente_alim_principal_std(A)':
                self.current1_std = float(words[1])
            elif first == 'corrente_alim_secundaria_avg(A)':
                self.current2_avg = float(words[1])
            elif first == 'corrente_alim_secundaria_std(A)':
                self.current2_std = float(words[1])
            elif first == 'data':
                self.timestamp += words[1] + ' '
            elif first == 'hora':
                self.timestamp += words[1] + ' '
            elif first == 'temperatura_ima(C)':
                self.temperature = words[1]

        self.timestamp = self.timestamp.strip()

        # reads raw data section of file and converts data into numpy matrix
        raw_data = _numpy.zeros((self.nr_turns, self.nr_points))

        idx = 0
        for line in raw_data_txt:
            raw_data[:,idx] = [self.rotating_coil.data_conversion_factor*float(word) for word in line.split(' ')]
            idx += 1

        # subtracts linear drift
        linear_data = _numpy.reshape(raw_data,(1,-1))
        raw_data = _numpy.cumsum(_numpy.roll(linear_data, 0))
        avg_raw_data = _numpy.mean(linear_data)
        linear_data_corrected = linear_data - avg_raw_data
        raw_data_corrected = _numpy.cumsum(_numpy.roll(linear_data_corrected, 0))
        raw_data_corrected = _numpy.reshape(raw_data_corrected, (self.nr_turns,-1))

        self.flux = raw_data_corrected

        self.magnet_label = _os.path.basename(self.filename).split('_')[0]
        #self.magnet_label = _os.path.basename(self.filename)[:5]

    def __str__(self):
        r = ''
        r += '\n\n--- measurement ---\n\n'
        r += '{0:<30s} {1:s}'.format('data_filename', _os.path.basename(self.filename))
        r += '\n{0:<30s} {1:s}'.format('time_stamp', self.timestamp)
        r += '\n{0:<30s} {1:+.3f}'.format('main_current_avg[A]', self.current1_avg)
        r += '\n{0:<30s} {1:.3f}'.format('main_current_std[A]', self.current1_std)
        r += '\n{0:<30s} {1:s}'.format('temperature', self.temperature)
        r += '\n{0:<30s} {1:d}'.format('nr_turns', self.nr_turns)
        r += '\n\n--- rotating coil ---\n\n'
        r += self.rotating_coil.__str__()
        return r

class Analysis:

    def __init__(self, measurement, parameters):
        self.measurement = measurement
        self.parameters = parameters
        self._run_analysis()

    def _run_analysis(self):
        if self.measurement.rotating_coil.type == 'Bobina Radial':
            self._calc_multipoles_for_radial_rotcoil()
            self._calc_relative_multipoles()
            self._calc_statistics()
        else:
            Exception('rotating coil type not defined!')

    def _calc_multipoles_for_radial_rotcoil(self):

        flux_data = self.measurement.flux
        harmonics = self.parameters.harmonics
        rotcoil = self.measurement.rotating_coil
        nr_coil_turns = rotcoil.nr_coil_turns
        r1, r2 = rotcoil.inner_radius, rotcoil.outer_radius
        nr_turns = self.measurement.nr_turns

        self.polynom_a = _numpy.zeros((nr_turns,len(harmonics)))
        self.polynom_b = _numpy.zeros((nr_turns,len(harmonics)))
        for i in range(nr_turns):
            flux = flux_data[i,:]
            angle = _numpy.linspace(2*_math.pi/len(flux),2*_math.pi,len(flux))
            polynom_a, polynom_b, polynom = [],[],[]
            for harmonic in harmonics:
                #n = harmonic + 1
                n = harmonic
                integ_cos = _numpy.trapz(flux * _numpy.cos(n*angle), angle)
                integ_sin = _numpy.trapz(flux * _numpy.sin(n*angle), angle)
                #if (n == 3):
                #    print(integ_cos, integ_sin)
                b_n = +n * integ_cos / _math.pi / nr_coil_turns / (r2**n - r1**n)
                a_n = -n * integ_sin / _math.pi / nr_coil_turns / (r2**n - r1**n)
                polynom_b.append(b_n)
                polynom_a.append(a_n)
            self.polynom_a[i,:] = polynom_a
            self.polynom_b[i,:] = polynom_b

    def _calc_relative_multipoles(self):

        ref_radius = self.parameters.ref_radius
        main_harmonic = self.parameters.main_multipole_harmonic
        main_is_skew = self.parameters.main_multipole_is_skew
        harmonics = self.parameters.harmonics

        idx = harmonics.index(main_harmonic)
        if main_is_skew:
            main_multipole_at_ref_radius = self.polynom_a[:,idx] * (ref_radius ** (main_harmonic-1))
        else:
            main_multipole_at_ref_radius = self.polynom_b[:,idx] * (ref_radius ** (main_harmonic-1))
        self.polynom_a_relative = 0 * self.polynom_a
        self.polynom_b_relative = 0 * self.polynom_b
        for i in range(len(harmonics)):
            polynom_a_at_rotcoil_radius = self.polynom_a[:,i] * (ref_radius ** (harmonics[i]-1))
            polynom_b_at_rotcoil_radius = self.polynom_b[:,i] * (ref_radius ** (harmonics[i]-1))
            self.polynom_a_relative[:,i] = polynom_a_at_rotcoil_radius / main_multipole_at_ref_radius
            self.polynom_b_relative[:,i] = polynom_b_at_rotcoil_radius / main_multipole_at_ref_radius

    def _calc_statistics(self):
        self.polynom_a_relative_avg = _numpy.mean(self.polynom_a_relative, axis=0)
        self.polynom_a_relative_std = _numpy.std(self.polynom_a_relative, axis=0)
        self.polynom_b_relative_avg = _numpy.mean(self.polynom_b_relative, axis=0)
        self.polynom_b_relative_std = _numpy.std(self.polynom_b_relative, axis=0)
        self.polynom_a_avg = _numpy.mean(self.polynom_a, axis=0)
        self.polynom_a_std = _numpy.std(self.polynom_a, axis=0)
        self.polynom_b_avg = _numpy.mean(self.polynom_b, axis=0)
        self.polynom_b_std = _numpy.std(self.polynom_b, axis=0)

    def __str__(self):
        r = ''
        r += self.measurement.__str__()
        r += self.parameters.__str__()
        r += '\n\n--- analysis ---\n\n'
        r += 'integrated multipoles relative to main multipole at r0 = {0:.1f} mm\n'.format(1000*self.parameters.ref_radius)
        r += '{0:2s}  {1:<11s}  {2:<11s}  {3:<11s}  {4:<11s}    {5:<10s}\n'.format('n', 'avg (Nn/Mp)','std (Nn/Mp)','avg (Sn/Mp)','std (Sn/Mp)', 'label')
        for i in range(len(self.parameters.harmonics)):
            h = self.parameters.harmonics[i]
            r += '{0:02d}  {1:+.4e}  {2:.4e}   {3:+.4e}  {4:.4e}    {5:<10s}\n'.format(h, self.polynom_b_relative_avg[i], self.polynom_b_relative_std[i], self.polynom_a_relative_avg[i], self.polynom_a_relative_std[i], get_harmonic_label(h))
        r += '\n'
        r += 'integrated multipoles\n'
        r += '{0:2s}  {1:<11s}  {2:<11s}  {3:<11s}  {4:<11s}    {5:<10s}    units\n'.format('n', '  avg (Nn)',' std (Nn)','  avg (Sn)',' std (Sn)', 'label')
        for i in range(len(self.parameters.harmonics)):
            h = self.parameters.harmonics[i]
            r += '{0:02d}  {1:+.4e}  {2:.4e}   {3:+.4e}  {4:.4e}    {5:<10s}    {6:s}\n'.format(h, self.polynom_b_avg[i], self.polynom_b_std[i], self.polynom_a_avg[i], self.polynom_a_std[i], get_harmonic_label(h), get_units(h))
        r += '\nlegend:'
        r += '\nNn   : integrated normal 2n-polar field'
        r += '\nSn   : integrated skew 2n-polar field'
        r += '\nNn/Mp: integrated normal 2n-polar field relative ar r0 by integral of main multipole at r0'
        r += '\nSn/Mp: integrated skew 2n-polar field relative at r0 by integral of main multipole at r0'
        r += '\navg(): measurement average'
        r += '\navg(): measurement stddev'
        return r

class AnalysisFromSet:
    pass

def get_harmonic_label(n):
    dic = {1:'dipole',2:'quadrupole',3:'sextupole',4:'octupole',5:'decapole',6:'dudecapole'}
    if n > max(dic.keys()):
        return str(n)
    else:
        return dic[n]

def get_units(harmonic):
    try:
        u = {1:'T.m',2:'T',3:'T/m'}[harmonic]
    except:
        u = 'T/m^{0:d}'.format(harmonic-2)
    return u

def run_analysis(parms, fnames, print_flag=True):
    '''reads data and does multipolar analysis, store it in "analysis"'''
    analysis = []
    for i in range(len(fnames)):
        m = Measurement(fnames[i])
        a = Analysis(m,parms)
        *first, fname = _os.path.split(fnames[i])
        if print_flag: print('{3:02d} - {0}, current: {1:+9.4f} +/- {2:.4f}'.format(fname, a.measurement.current1_avg, a.measurement.current1_std, i))
        analysis.append(a)
    return analysis

def select_ramp_up(current_data_set, current_threshold=0.0):
    '''returns only ramp up points from a data set'''
    c = -float('inf')
    new_data = []
    for data in current_data_set:
        current = data.measurement.current1_avg
        if abs(current) > current_threshold:
            if current < c: break
            c = current
            new_data.append(data)
    return new_data

def get_remanent_field(magnet_data_set, current_threshold=0.0):
    '''does a remanent field calculation from a data set'''
    currents, polya, polyb = [], [], []
    for d1 in magnet_data_set:
        for d2 in d1:
            for d3 in d2:
                current = d3.measurement.current1_avg
                if abs(current) < current_threshold:
                    currents.append(current)
                    polya.append(d3.polynom_a_avg)
                    polyb.append(d3.polynom_b_avg)
    polya_avg = _numpy.mean(_numpy.array(polya), axis=0)
    polya_std = _numpy.std (_numpy.array(polya), axis=0)
    polyb_avg = _numpy.mean(_numpy.array(polyb), axis=0)
    polyb_std = _numpy.std (_numpy.array(polyb), axis=0)
    return ((polya_avg,polya_std),(polyb_avg,polyb_std))

def get_maximum_main_multipole(current_data_set, parms, current_threshold=0.0):
    '''returns maximum value of main multipole'''
    selection = select_ramp_up(current_data_set, current_threshold=current_threshold)
    idx = parms.harmonics.index(parms.main_multipole_harmonic)
    max_current = selection[-1].measurement.current1_avg
    if parms.main_multipole_is_skew:
        return selection[-1].polynom_a_avg[idx], max_current
    else:
        return selection[-1].polynom_b_avg[idx], max_current

def get_multipole_from_data_set(current_data_set, parms, h, mtype='normal', relative=False):
    '''exctracts from data set current and multipole'''
    idx = parms.harmonics.index(h)
    c, multipole = [],[]
    for d in current_data_set:
        c.append(d.measurement.current1_avg)
        if mtype == 'normal':
            if relative:
                multipole.append(d.polynom_b_relative_avg[idx])
            else:
                multipole.append(d.polynom_b_avg[idx])
        else:
            if relative:
                multipole.append(d.polynom_a_relative_avg[idx])
            else:
                multipole.append(d.polynom_a_avg[idx])
    return c, multipole

def find_current(meas_data_set, parms, multi_norm, energy, mtype='normal', current_threshold=0.0):
    brho, *_ = _mp.beam_optics.beam_rigidity(energy = energy)
    multi = - multi_norm * brho
    current = []
    for d1 in meas_data_set:
        d2 = select_ramp_up(d1, current_threshold)
        c,s = get_multipole_from_data_set(d2, parms, parms.main_multipole_harmonic, mtype='normal', relative=False)
        tc = _numpy.interp(multi, s, c, left=float('nan'), right=float('nan'))
        current.append(tc)
    return _numpy.mean(current), _numpy.std(current)

def plot_relative_multipoles(magnet_data_set, parms, h, mtype='normal', current_threshold=0.0, currents = None, label = None, show=True, save=False, xlim=None, ylim=None, ax=None):

    if currents is None: currents = []
    idx = parms.harmonics.index(h)

    try:
        r = parms.multipoles_spec[h]
    except:
        r = ((0,0),(0,0))
    if mtype == 'normal':
        sys, rms = r[0][0], r[1][0]
    else:
        sys, rms = r[0][1], r[1][1]

    if not ax:
        f = _plt.figure()
        ax = _plt.axes()
    f = ax.get_figure()

    #_plt.clf()
    c = -1
    minc, maxc = float('inf'), -float('inf')
    minm, maxm = float('inf'), -float('inf')
    for d1 in magnet_data_set:
        c = (c+1) % len(_colors)
        first = True
        for d2 in d1:
            d3 = select_ramp_up(d2, current_threshold)
            if not label:
                label = d3[0].measurement.magnet_label
            cu, mu = get_multipole_from_data_set(d3, parms, h=h, mtype=mtype, relative=True)
            if first:
                ax.plot(cu, mu, _colors[c], label = label)
                first = False
            else:
                ax.plot(cu, mu, _colors[c])
            minc, maxc = min([minc,min(cu)]), max([maxc,max(cu)])
            minm, maxm = min([minm,min(mu)]), max([maxm,max(mu)])
    ax.plot([minc,maxc],[sys-rms,sys-rms], 'k--')
    ax.plot([minc,maxc],[sys,sys], 'k')
    ax.plot([minc,maxc],[sys+rms,sys+rms], 'k--')
    for c in currents:
        ax.plot([c,c],[minm,maxm], 'k--')
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xlabel('Current [A]')
    ax.set_ylabel('Relative multipole @ r = ' + str(parms.ref_radius * 1000) + ' mm')
    ax.grid('on')
    ax.legend(loc='best')
    if mtype == 'normal':
        ax.set_title('Relative normal multipole of order h = ' + str(h))
        if save: f.savefig('relative_normal_multipole_{0:02d}.png'.format(h))
    else:
        ax.set_title('Relative skew multipole of order h = ' + str(h))
        if save: f.savefig('relative_skew_multipole_{0:02d}.png'.format(h))
    if show: f.show()

def plot_magnetic_center(magnet_data_set, parms, mtype='normal', currents=None, show=True, save=False, xlim=None, ylim=None, ax=None):

    if not ax:
        f = _plt.figure()
        ax = _plt.axes()
    f = ax.get_figure()

    if currents is None: currents = []
    sys = 0.0
    if mtype == 'normal':
        rms = parms.x_misalignment_spec
    else:
        rms = parms.y_misalignment_spec

    c = -1
    minc, maxc = float('inf'), -float('inf')
    minm, maxm = float('inf'), -float('inf')
    for d1 in magnet_data_set:
        c = (c+1) % len(_colors)
        first = True
        for d2 in d1:
            d3 = select_ramp_up(d2, 0.5)
            cu,mu = [],[]
            for d4 in d3:
                cu.append(d4.measurement.current1_avg)
                if parms.main_multipole_harmonic == 3:
                    D = d4.polynom_b_avg[0] + d4.polynom_a_avg[0] * 1j
                    Q = d4.polynom_b_avg[1] + d4.polynom_a_avg[1] * 1j
                    S = d4.polynom_b_avg[2] + d4.polynom_a_avg[2] * 1j
                    z = _numpy.roots([2*S,Q])
                elif parms.main_multipole_harmonic == 2:
                    D = d4.polynom_b_avg[0] + d4.polynom_a_avg[0] * 1j
                    Q = d4.polynom_b_avg[1] + d4.polynom_a_avg[1] * 1j
                    z = _numpy.roots([Q,D])
                if mtype == 'normal':
                    mu.append(1e6*z.real)
                else:
                    mu.append(1e6*z.imag)
            if first:
                ax.plot(cu, mu, _colors[c], label=d4.measurement.magnet_label)
                first = False
            else:
                ax.plot(cu, mu, _colors[c])
            minc, maxc = min([minc,min(cu)]), max([maxc,max(cu)])
            minm, maxm = min([minm,min(mu)]), max([maxm,max(mu)])
    ax.plot([minc,maxc],[sys-rms,sys-rms], 'k--')
    ax.plot([minc,maxc],[sys,sys], 'k')
    ax.plot([minc,maxc],[sys+rms,sys+rms], 'k--')
    ax.legend(loc='best')
    for c in currents:
        if ylim:
            ax.plot([c,c],ylim,'k--')
        else:
            ax.plot([c,c],[mimm,maxm],'k--')
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xlabel('Current [A]')
    ax.set_ylabel('Position [um]')
    ax.grid('on')
    if mtype == 'normal':
        if parms.main_multipole_harmonic == 3:
            ax.set_title('Horizontal position where quadrupole vanishes')
        else:
            ax.set_title('Horizontal position where field vanishes')
        if save: f.savefig('magnetic_center_horizontal_pos.png')
    else:
        if parms.main_multipole_harmonic == 3:
            ax.set_title('Vertical position where field vanishes')
        else:
            ax.set_title('Vertical position where field vanishes')
        if save: f.savefig('magnetic_center_vertical_pos.png')
    if show: f.show()

def plot_rotation_angle(magnet_data_set, parms, currents = None, show=True, save=False, ax=None):

    if not ax:
        f = _plt.figure()
        ax = _plt.axes()
    f = ax.get_figure()

    sys = 0
    rms = parms.roll_rotation_spec  #[mrad]
    idx = parms.harmonics.index(parms.main_multipole_harmonic)

    if currents is None: currents = []

    c = -1
    minc, maxc = float('inf'), -float('inf')
    minm, maxm = float('inf'), -float('inf')
    for d1 in magnet_data_set:
        c = (c+1) % len(_colors)
        first = True
        for d2 in d1:
            d3 = select_ramp_up(d2, 0.5)
            cu,mu = [],[]
            for d4 in d3:
                cu.append(d4.measurement.current1_avg)
                angle = _math.atan(d4.polynom_a_avg[idx]/d4.polynom_b_avg[idx])/(parms.main_multipole_harmonic-1)
                mu.append(1e3*angle)
            if first:
                ax.plot(cu, mu, _colors[c], label = d4.measurement.magnet_label)
                first = False
            else:
                ax.plot(cu, mu, _colors[c])
            minc, maxc = min([minc,min(cu)]), max([maxc,max(cu)])
            minm, maxm = min([minm,min(mu)]), max([maxm,max(mu)])
    ax.plot([minc,maxc],[sys-rms,sys-rms], 'k--')
    ax.plot([minc,maxc],[sys,sys], 'k')
    ax.plot([minc,maxc],[sys+rms,sys+rms], 'k--')
    for c in currents:
        ax.plot([c,c],[minm,maxm], 'k--')
    ax.legend(loc='best')
    for c in currents:
        ax.plot([c,c],[minm,maxm],'k--')
    ax.set_xlabel('Current [A]')
    ax.set_ylabel('Rotation angle [mrad]')
    ax.grid('on')
    ax.set_title('Rotation angle from main multipole component')
    if save: f.savefig('rotation_angle.png')
    if show: f.show()

def print_multipoles_single_magnet(meas_data_set, parms, current, current_threshold=0.0):

    pa_avg, pa_std = [], []
    pb_avg, pb_std = [], []
    strheader = '{0:<3s}: {1:^24s}  |  {2:^24s}'.format('h','(B_n/B_2)@r_0', '(A_n/B_2)@r_0')
    print(strheader)
    print('-'*len(strheader))
    for i in range(len(parms.harmonics)):
        polya, polyb = [], []
        for d2 in meas_data_set: # over different measurements for one magnet
            d = select_ramp_up(d2, current_threshold)
            curr, tpolya, tpolyb = [], [], []
            for d3 in d: # over different excitation currents
                curr.append(d3.measurement.current1_avg)
                tpolya.append(d3.polynom_a_relative_avg[i])
                tpolyb.append(d3.polynom_b_relative_avg[i])
            polya.append(_numpy.interp(current, curr, tpolya))
            polyb.append(_numpy.interp(current, curr, tpolyb))
        polya_avg, polya_std = _numpy.mean(polya), _numpy.std(polya)
        polyb_avg, polyb_std = _numpy.mean(polyb), _numpy.std(polyb)
        pa_avg.append(polya_avg), pa_std.append(polya_std)
        pb_avg.append(polyb_avg), pb_std.append(polyb_std)
        print('{0:02d} : {1:+.3e} +/- {2:.3e}  |  {3:+.3e} +/- {4:.3e}'.format(parms.harmonics[i], polyb_avg, polyb_std, polya_avg, polya_std))

    print()
    print('harmonics: ', _numpy.array(parms.harmonics)-1)
    print('polyb_avg: [', end='')
    for d in pb_avg: print('{0:+.1e},'.format(d), end='')
    print(']')
    print('polyb_std: [', end='')
    for d in pb_std: print('{0:+.1e},'.format(d), end='')
    print(']')
    print('polya_avg: [', end='')
    for d in pa_avg: print('{0:+.1e},'.format(d), end='')
    print(']')
    print('polya_std: [', end='')
    for d in pa_std: print('{0:+.1e},'.format(d), end='')
    print(']')

def get_excitation_curve(current_data_set, parms):
    harmonics = parms.harmonics
    main_harmonic_idx = harmonics.index(parms.main_multipole_harmonic)
    data = select_ramp_up(current_data_set, current_threshold=0.0)
    currents, exctable = [], _numpy.zeros((len(data),len(harmonics)))
    for i in range(len(data)):
        currents.append(data[i].measurement.current1_avg)
        exctable[i,:] = data[i].polynom_b_avg
    exccurve = exctable[:,main_harmonic_idx]
    return exctable, exccurve, currents, harmonics, parms.main_multipole_harmonic

def get_average_excitation_curve(meas_data_set, parms, currents):
    multipoles_sum1, multipoles_sum2 = 0,0
    for current_data_set in meas_data_set:
        exctable, exccurve, tcurrents, harmonics, main_multipole_harmonic = get_excitation_curve(current_data_set, parms)
        multipoles = _numpy.zeros((len(currents),len(harmonics)))
        for i in range(len(harmonics)):
            multipoles[:,i] = _numpy.interp(currents, tcurrents, exctable[:,i], left=float('nan'), right=float('nan'))
        multipoles_sum1 += multipoles
        multipoles_sum2 += multipoles**2
    idx = harmonics.index(main_multipole_harmonic)
    exctable_avg = multipoles_sum1/len(meas_data_set)
    exctable_std = _numpy.sqrt(multipoles_sum2/len(meas_data_set) - exctable_avg**2)
    exccurve_avg = exctable_avg[:,idx]
    exccurve_std = exctable_std[:,idx]
    return (exctable_avg, exccurve_avg, currents, harmonics, parms.main_multipole_harmonic), (exctable_std, exccurve_std)
    #return currents, multipoles_avg, multipoles_std, harmonics, main_multipole_harmonic

def plot_excitation_curve(meas_data_set, parms, currents, show=True, save=False, ax=None):

    excitation_curve, std = get_average_excitation_curve(meas_data_set, parms, currents)
    exctable, exccurve, currents, harmonics, main_multipole_harmonic = excitation_curve
    exctable_std, exccurve_std = std
    #currents, multipoles_avg, multipoles_std, harmonics, main_multipole_harmonic = get_average_excitation_curve(meas_data_set, parms, currents)

    if not ax:
        f = _plt.figure()
        ax = _plt.axes()
    f = ax.get_figure()

    idx = harmonics.index(main_multipole_harmonic)
    ax.plot(currents, exctable[:,idx] - exctable_std[:,idx], '--', color='blue')
    ax.plot(currents, exctable[:,idx], color = 'blue')
    ax.plot(currents, exctable[:,idx] + exctable_std[:,idx], '--', color = 'blue')

    ax.set_xlabel('Current [A]')
    ax.set_ylabel(_labels[main_multipole_harmonic])
    ax.grid('on')
    ax.set_title('Excitation Curve')
    if save: f.savefig('excitation_curve.png')
    if show: f.show()

def calc_excitation_curve_nonlinearity(meas_data_set, parms, currents, show=True, save=False, ax=None):

    if not ax:
        f = _plt.figure()
        ax = _plt.axes()
    f = ax.get_figure()
    currents = _numpy.array(currents)

    excitation_curve, std = get_average_excitation_curve(meas_data_set, parms, currents)
    exctable, exccurve, currents, harmonics, main_multipole_harmonic = excitation_curve
    exctable_std, exccurve_std = std

    idx = harmonics.index(main_multipole_harmonic)
    pfit = _numpy.poly1d(_numpy.polyfit(currents, exctable[:,idx], idx));
    fit_error = _numpy.array([(pfit(currents[i]) - exctable[i,idx])/exctable[i,idx] for i in range(len(currents))])

    c_neg, f_neg = currents[fit_error < 0], fit_error[fit_error < 0]
    c_pos, f_pos = currents[fit_error > 0], fit_error[fit_error > 0]

    ax.set_yscale('log')
    ax.plot(c_pos, 100*f_pos, 'bo');
    ax.plot(c_neg, 100*abs(f_neg), 'ro');

    ax.set_ylim([1e-2,100])
    ax.set_xlabel('Current [A]')
    ax.set_ylabel('Nonlinearity [%]')
    ax.grid('on')
    ax.set_title('Nonlinearity of Excitation Curve')
    if save: f.savefig('excitation_curve_non_linearity.png')
    if show: f.show()

    return fit_error, currents, exctable

def plot_hysteresis(meas_data_set, parms, excitation_curve, legends=None, show=True, save=False, ax=None):

    if not ax:
        f = _plt.figure()
        ax = _plt.axes()
    f = ax.get_figure()

    exctable, exccurve, currents, harmonics, main_multipole_harmonic = excitation_curve

    currents = _numpy.array(currents)
    idx = harmonics.index(main_multipole_harmonic)
    pfit = _numpy.poly1d(_numpy.polyfit(currents, exctable[:,idx], idx));
    for current_data_set in meas_data_set:
        c,m = [],[]
        for data in current_data_set:
            main_multipole = data.polynom_b_avg[idx]
            current = data.measurement.current1_avg
            multipole = main_multipole - pfit(current)
            c.append(current)
            m.append(multipole)
        ax.plot(c,m)
    ax.set_xlabel('Current [A]')
    ax.set_ylabel(_labels[main_multipole_harmonic])
    ax.grid('on')
    ax.set_title('Histeresis of Main Multipole')
    if legends: ax.legend(legends)
    if save: f.savefig('histeresis.png')
    if show: f.show()

def plot_trim_coil_excitation_curves(meas_data_set, parms, excitation_curve, legends=None, show=True, save=False, ax=None):

    if not ax:
        f = _plt.figure()
        ax = _plt.axes()
    f = ax.get_figure()

    harmonics = parms.harmonics
    idx = harmonics.index(parms.main_multipole_harmonic)


    exctable, exccurve, currents, harmonics, main_multipole_harmonic = excitation_curve

    pfit = _numpy.poly1d(_numpy.polyfit(currents, exctable[:,idx], idx));

    for data in meas_data_set:
        x,y = [],[]
        main_current = data[0].measurement.current1_avg
        for i in range(len(data)):
            x.append(data[i].measurement.current2_avg)
            y.append(data[i].polynom_b_avg[idx] - pfit(main_current))
            #y.append(100*(data[i].polynom_b_avg[idx] - pfit(main_current))/pfit(main_current))
        ax.plot(x,y)
    ax.set_xlabel('Trim Current [A]')
    ax.set_ylabel('Variation of ' + _labels[main_multipole_harmonic])
    if legends: ax.legend(legends, loc='best')
    if save: f.savefig('trim_coil_excitation.png')
    if show: f.show()
