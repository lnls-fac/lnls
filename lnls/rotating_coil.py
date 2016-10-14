#!/usr/bin/env python3

#import lnls

#import shutil
import os as _os
import numpy as _numpy
import math as _math
import mathphys as _mp
import matplotlib.pyplot as _plt
import matplotlib.ticker as _ticker
import os as _os
import lnls.utils as _utils
#from . import utils as _utils

_colors = ['blue', 'red', 'green', 'orange', 'black']

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
        self.max_integ_mult_spec = 5.2116 #[T]
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
            r += '{0:02d}  {1:+.4e}  {2:.4e}   {3:+.4e}  {4:.4e}    {5:<10s}\n'.format(h, self.polynom_b_relative_avg[i], self.polynom_b_relative_std[i], self.polynom_a_relative_avg[i], self.polynom_a_relative_std[i], get_label(h))
        r += '\n'
        r += 'integrated multipoles\n'
        r += '{0:2s}  {1:<11s}  {2:<11s}  {3:<11s}  {4:<11s}    {5:<10s}    units\n'.format('n', '  avg (Nn)',' std (Nn)','  avg (Sn)',' std (Sn)', 'label')
        for i in range(len(self.parameters.harmonics)):
            h = self.parameters.harmonics[i]
            r += '{0:02d}  {1:+.4e}  {2:.4e}   {3:+.4e}  {4:.4e}    {5:<10s}    {6:s}\n'.format(h, self.polynom_b_avg[i], self.polynom_b_std[i], self.polynom_a_avg[i], self.polynom_a_std[i], get_label(h), get_units(h))
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

def run_analysis(parms, fnames):
    '''reads data and does multipolar analysis, store it in "analysis"'''
    analysis = []
    for i in range(len(fnames)):
        m = Measurement(fnames[i])
        a = Analysis(m,parms)
        *first, fname = _os.path.split(fnames[i])
        print('{3:02d} - {0}, current: {1:+9.4f} +/- {2:.4f}'.format(fname, a.measurement.current1_avg, a.measurement.current1_std, i))
        analysis.append(a)
    return analysis

def select_ramp_up(data, current_threshold=0.0):
    '''returns only ramp up points from a data set'''
    c = -float('inf')
    new_data = []
    for d in data:
        current = d.measurement.current1_avg
        if abs(current) > current_threshold:
            if current < c: break
            c = current
            new_data.append(d)
    return new_data

def get_remanent_field(data, current_threshold=0.0):
    '''does a remanent field calculation from a data set'''
    currents, polya, polyb = [], [], []
    for d1 in data:
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

def get_maximum_main_multipole(data, parms, current_threshold=0.0):
    '''returns maximum value of main multipole'''
    selection = select_ramp_up(data, current_threshold=current_threshold)
    idx = parms.harmonics.index(parms.main_multipole_harmonic)
    if parms.main_multipole_is_skew:
        return selection[-1].polynom_a_avg[idx]
    else:
        return selection[-1].polynom_b_avg[idx]

def get_multipole_from_data_set(data, parms, h, mtype='normal', relative=False):
    '''exctracts from data set current and multipole'''
    idx = parms.harmonics.index(h)
    c, multipole = [],[]
    for d in data:
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

def plot_excitation_curve(data, parms, current_grid, current_threshold=0.0, show=False):
    '''plots excitation curves'''

    labels = {1:'Integrated dipole [T.m]',
              2:'Integrated quadrupole [T]',
              3:'Integrated sextupole [T/m]'}

    c = current_grid
    I = _numpy.linspace(c[0], c[1], c[2])
    _plt.clf()
    c = -1
    multi_set = []
    for meas_set in data:
        c = c+1 % len(_colors)
        for d2 in meas_set:
            d = select_ramp_up(d2, current_threshold)
            currents, multipoles = get_multipole_from_data_set(d, parms, h=parms.main_multipole_harmonic, mtype='normal', relative=False)
            interpMulti = _numpy.interp(I, currents, multipoles)
            multi_set.append(interpMulti)
            _plt.plot(currents, multipoles, _colors[c])

    _plt.xlabel('Current [A]')
    _plt.ylabel(labels[parms.main_multipole_harmonic])
    _plt.grid('on')
    _plt.title('Excitation Curves')
    _plt.savefig('excitation_curves.png')
    if show: _plt.show()

    _plt.clf()
    multi_avg = _numpy.mean(_numpy.array(multi_set), axis=0)
    multi_std = _numpy.std(_numpy.array(multi_set), axis=0)
    for d in multi_set:
        difrel = 100*(_numpy.array(d) - multi_avg)/multi_avg
        _plt.plot(I, difrel, 'b')
    mult_rel = 100*multi_std/abs(multi_avg)
    _plt.plot(I,mult_rel, 'r')
    rms_spec = parms.excitation_rms_spec
    _plt.plot([min(I), max(I)], [+rms_spec, +rms_spec], 'k--')
    _plt.plot([min(I), max(I)], [-rms_spec, -rms_spec], 'k--')
    _plt.grid('on')
    _plt.xlabel('Current [A]')
    _plt.ylabel('excitation curve rms [%]')
    _plt.title('Excitation Curve R.M.S')
    _plt.savefig('excitation_curve_rms.png')
    if show: _plt.show()

    return I, mult_rel

def plot_relative_multipoles(data, parms, h, mtype='normal', current_threshold=0.0, currents = None, show=False, xlim=None, ylim=None):

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

    _plt.clf()
    c = -1
    minc, maxc = float('inf'), -float('inf')
    minm, maxm = float('inf'), -float('inf')
    for d1 in data:
        c = (c+1) % len(_colors)
        first = True
        for d2 in d1:
            d3 = select_ramp_up(d2, current_threshold)
            cu, mu = get_multipole_from_data_set(d3, parms, h=h, mtype=mtype, relative=True)
            if first:
                _plt.plot(cu, mu, _colors[c], label = d3[0].measurement.magnet_label)
                first = False
            else:
                _plt.plot(cu, mu, _colors[c])
            minc, maxc = min([minc,min(cu)]), max([maxc,max(cu)])
            minm, maxm = min([minm,min(mu)]), max([maxm,max(mu)])
    _plt.plot([minc,maxc],[sys-rms,sys-rms], 'k--')
    _plt.plot([minc,maxc],[sys,sys], 'k')
    _plt.plot([minc,maxc],[sys+rms,sys+rms], 'k--')
    for c in currents:
        _plt.plot([c,c],[minm,maxm], 'k--')
    if xlim: _plt.xlim(xlim)
    if ylim: _plt.ylim(ylim)
    _plt.xlabel('Current [A]')
    _plt.ylabel('Relative multipole @ r = ' + str(parms.ref_radius * 1000) + ' mm')
    _plt.grid('on')
    _plt.legend(loc='best')
    if mtype == 'normal':
        _plt.title('Relative normal multipole of order h = ' + str(h))
        _plt.savefig('relative_normal_multipole_{0:02d}.png'.format(h))
    else:
        _plt.title('Relative skew multipole of order h = ' + str(h))
        _plt.savefig('relative_skew_multipole_{0:02d}.png'.format(h))
    if show: _plt.show()

def find_current(data, parms, multi_norm, energy, mtype='normal', current_threshold=0.0):
    brho, *_ = _mp.beam_optics.beam_rigidity(energy = energy)
    multi = - multi_norm * brho
    current = []
    for d1 in data:
        d2 = select_ramp_up(d1, current_threshold)
        c,s = get_multipole_from_data_set(d2, parms, parms.main_multipole_harmonic, mtype='normal', relative=False)
        tc = _numpy.interp(multi, s[::-1], c[::-1], left=float('nan'), right=float('nan'))
        current.append(tc)
    return _numpy.mean(current), _numpy.std(current)

def plot_magnetic_center(data, parms, mtype='normal', currents=None, show=False):

    if currents is None: currents = []
    sys = 0.0
    if mtype == 'normal':
        rms = parms.x_misalignment_spec
    else:
        rms = parms.y_misalignment_spec

    _plt.clf()
    c = -1
    minc, maxc = float('inf'), -float('inf')
    minm, maxm = float('inf'), -float('inf')
    for d1 in data:
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
                _plt.plot(cu, mu, _colors[c], label=d4.measurement.magnet_label)
                first = False
            else:
                _plt.plot(cu, mu, _colors[c])
            minc, maxc = min([minc,min(cu)]), max([maxc,max(cu)])
            minm, maxm = min([minm,min(mu)]), max([maxm,max(mu)])
    _plt.plot([minc,maxc],[sys-rms,sys-rms], 'k--')
    _plt.plot([minc,maxc],[sys,sys], 'k')
    _plt.plot([minc,maxc],[sys+rms,sys+rms], 'k--')
    _plt.legend(loc='best')
    for c in currents:
        _plt.plot([c,c],[minm,maxm],'k--')
    _plt.xlabel('Current [A]')
    _plt.ylabel('Position [um]')
    _plt.grid('on')
    if mtype == 'normal':
        _plt.title('Horizontal position where quadrupole vanishes')
        _plt.savefig('magnetic_center_horizontal_pos.png')
    else:
        _plt.title('Vertical position where quadrupole vanishes')
        _plt.savefig('magnetic_center_vertical_pos.png')
    if show: _plt.show()

def plot_rotation_angle(data, parms, currents = None, show=False):

    sys = 0
    rms = parms.roll_rotation_spec  #[mrad]
    idx = parms.harmonics.index(parms.main_multipole_harmonic)

    if currents is None: currents = []

    _plt.clf()
    c = -1
    minc, maxc = float('inf'), -float('inf')
    minm, maxm = float('inf'), -float('inf')
    for d1 in data:
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
                _plt.plot(cu, mu, _colors[c], label = d4.measurement.magnet_label)
                first = False
            else:
                _plt.plot(cu, mu, _colors[c])
            minc, maxc = min([minc,min(cu)]), max([maxc,max(cu)])
            minm, maxm = min([minm,min(mu)]), max([maxm,max(mu)])
    _plt.plot([minc,maxc],[sys-rms,sys-rms], 'k--')
    _plt.plot([minc,maxc],[sys,sys], 'k')
    _plt.plot([minc,maxc],[sys+rms,sys+rms], 'k--')
    for c in currents:
        _plt.plot([c,c],[minm,maxm], 'k--')
    _plt.legend(loc='best')
    for c in currents:
        _plt.plot([c,c],[minm,maxm],'k--')
    _plt.xlabel('Current [A]')
    _plt.ylabel('Rotation angle [mrad]')
    _plt.grid('on')
    _plt.title('Rotation angle from sextupole components')
    _plt.savefig('rotation_angle.png')
    if show: _plt.show()

def print_multipoles(data, parms, current, current_threshold=0.0):

    pa_avg, pa_std = [], []
    pb_avg, pb_std = [], []
    strheader = '{0:<3s}: {1:^24s}  |  {1:^24s}'.format('h','(B_n/B_2)@r_0', '(A_n/B_2)@r_0')
    print(strheader)
    print('-'*len(strheader))
    for i in range(len(parms.harmonics)):
        polya, polyb = [], []
        for d1 in data:
            t1polya, t1polyb = [], []
            for d2 in d1:
                d = select_ramp_up(d2, current_threshold)
                curr, tpolya, tpolyb = [], [], []
                for d3 in d:
                    curr.append(d3.measurement.current1_avg)
                    tpolya.append(d3.polynom_a_relative_avg[i])
                    tpolyb.append(d3.polynom_b_relative_avg[i])
                t1polya.append(_numpy.interp(current, curr, tpolya))
                t1polyb.append(_numpy.interp(current, curr, tpolyb))
            polya.append(_numpy.mean(t1polya))
            polyb.append(_numpy.mean(t1polyb))

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

def plot_magnetic_center2(data, parms, current, mtype='normal', current_threshold=0.0, show=False):

    if currents is None: currents = []
    sys = 0.0
    if mtype == 'normal':
        rms = parms.x_misalignment_spec
    else:
        rms = parms.y_misalignment_spec

    magnet,center = [],[]
    for d1 in data:
        z = 0 + 0*1j
        for d2 in d1:
            d3 = select_ramp_up(d2, 0.5)
            if parms.main_multipole_harmonic == 2:
                c,d = get_multipole_from_data_set(d3, parms, h=1, mtype='normal', relative=False)
                DB = _numpy.interp(current, c, d, left=float('nan'), right=float('nan'))
                c,d = get_multipole_from_data_set(d3, parms, h=1, mtype='skew', relative=False)
                DA = _numpy.interp(current, c, d, left=float('nan'), right=float('nan'))
                c,q = get_multipole_from_data_set(d3, parms, h=2, mtype='normal', relative=False)
                QB = _numpy.interp(current, c, q, left=float('nan'), right=float('nan'))
                c,q = get_multipole_from_data_set(d3, parms, h=2, mtype='skew', relative=False)
                QA = _numpy.interp(current, c, q, left=float('nan'), right=float('nan'))
                D = DB + 1j * DA
                Q = QB + 1j * QA
                z += _numpy.roots([Q,D])
        z /= len(d1)
        magnet.append(d3[0].measurement.magnet_name)
        if mtype == 'normal':
            center.append(1e6*z.real)
        else:
            center.append(1e6*z.imag)


# def get_units(harmonic):
#     try:
#         u = {1:'T.m',2:'T',3:'T/m'}[harmonic]
#     except:
#         u = 'T/m^{0:d}'.format(harmonic-2)
#     return u
#
# def get_label(harmonic):
#     try:
#         m = {1:'dipole',2:'quadrupole',3:'sextupole',4:'octupole',5:'decapole',6:'dudecapole'}[harmonic]
#     except:
#         m = '{0:d}-polar'.format(2*harmonic)
#     return m
#
# def get_multipole_name(harmonic,mult_type=None,ref_radius=None):
#     m = get_label(harmonic)
#     if ref_radius:
#         r = 'Normalized ' + mult_type.lower() + ' ' + m + ' (at r = {0:.1f} mm)'.format(1000*ref_radius)
#     else:
#         u = get_units(harmonic)
#         r = mult_type.lower().capitalize() + ' ' + m + ' / ({0:s})'.format(u)
#     return r
#
# def plot_harmonic_analysis(analysis, ylim=None, show_fig=False, save_fig=True):
#
#     if ylim is None: ylim=(1e-6,1e-1)
#     filename = _os.path.splitext(_os.path.basename(analysis.measurement.filename))[0]
#     harmonics = analysis.parameters.harmonics
#     ref_radius = analysis.parameters.ref_radius
#     # Normal
#     _plt.clf()
#     multipoles = analysis.polynom_b_relative_avg
#     multipoles_err = analysis.polynom_b_relative_std
#     colors = ['red' if value < 0 else 'blue' for value in multipoles]
#     _plt.bar(harmonics, abs(multipoles), yerr=multipoles_err, log=True, color=colors, error_kw=dict(elinewidth=2,ecolor='black'), align='center')
#     _plt.ylim(ylim)
#     _plt.title(filename)
#     _plt.xlabel('harmonics (dipole:n=1)')
#     _plt.ylabel('relative normal multipole (at r = {0:.1f} mm)'.format(1000*ref_radius))
#     _plt.grid('on')
#     if show_fig: _plt.show()
#     if save_fig: _plt.savefig('normal_' + filename + '.png')
#     # Skew
#     _plt.clf()
#     multipoles = analysis.polynom_a_relative_avg
#     multipoles_err = analysis.polynom_a_relative_std
#     colors = ['red' if value < 0 else 'blue' for value in multipoles]
#     _plt.bar(harmonics, abs(multipoles), yerr=multipoles_err, log=True, color=colors, error_kw=dict(elinewidth=2,ecolor='black'), align='center')
#     _plt.ylim(ylim)
#     _plt.title(filename)
#     _plt.xlabel('harmonics (dipole:n=1)')
#     _plt.ylabel('relative skew multipole (at r = {0:.1f} mm)'.format(1000*ref_radius))
#     _plt.grid('on')
#     if show_fig: _plt.show()
#     if save_fig: _plt.savefig('skew_' + filename + '.png')
#
# def plot_skew_rotation_vs_current(analysis_set, show_fig=False, save_fig=True):
#     current, angle, angle_err = [], [], []
#     for analysis in analysis_set:
#         parameters = analysis.parameters
#         data = analysis.measurement if hasattr(analysis,'measurement') else analysis
#         magnet_name = data.magnet_label
#         main_harmonic = analysis.parameters.main_multipole_harmonic
#         idx = analysis.parameters.harmonics.index(main_harmonic)
#         current.append(data.current1_avg)
#         a = 1000 * _math.atan(analysis.polynom_a_relative_avg[idx])/main_harmonic
#         angle.append(a)
#         #angle_err.append(1000 * )
#
#     limcurr = [min(current),max(current)]
#     maxy = [1.1*max(1e3*parameters.roll_rotation_spec, max(angle)), 1.1*min(-1e3*parameters.roll_rotation_spec, min(angle))]
#
#     _plt.clf()
#     _plt.plot(current, angle, '-o', color='darkgreen')
#     _plt.plot(limcurr,[1e3*parameters.roll_rotation_spec]*2, '--', color='black', label='Spec')
#     _plt.plot(limcurr,[-1e3*parameters.roll_rotation_spec]*2, '--', color='black')
#     _plt.ylim(maxy)
#     _plt.xlabel('current / A')
#     _plt.ylabel('Angle / mrad')
#     _plt.grid('on')
#     _plt.title('Residual skew of main multipolar field for ' + magnet_name)
#     _plt.legend(loc=0)
#     if show_fig: _plt.show()
#     if save_fig: _plt.savefig('skew_rotation_vs_current.png')
#
# def plot_dislocation_vs_current(analysis_set, show_fig=False, save_fig=True):
#     current, posx, posy = [], [], []
#     posx_err, posy_err = [], []
#     for analysis in analysis_set:
#         parameters = analysis.parameters
#         data = analysis.measurement if hasattr(analysis,'measurement') else analysis
#         magnet_name = data.magnet_label
#         main_harmonic = analysis.parameters.main_multipole_harmonic
#         idx = analysis.parameters.harmonics.index(main_harmonic)
#         current.append(data.current1_avg)
#         y = - 1e6*analysis.polynom_a_avg[0] / analysis.polynom_b_avg[1]
#         x = - 1e6*analysis.polynom_b_avg[0] / analysis.polynom_b_avg[1]
#         a,b   = 1e6*analysis.polynom_a_avg[0], 1e6*analysis.polynom_b_avg[1]
#         da,db = 1e6*analysis.polynom_a_std[0], 1e6*analysis.polynom_b_std[1]
#         y_err = ((a/b)**2*((da/a)**2+(db/b)**2))**0.5
#         a,b   = 1e6*analysis.polynom_b_avg[0], 1e6*analysis.polynom_b_avg[1]
#         da,db = 1e6*analysis.polynom_b_std[0], 1e6*analysis.polynom_b_std[1]
#         x_err = ((a/b)**2*((da/a)**2+(db/b)**2))**0.5
#         posx.append(x); posy.append(y)
#         posx_err.append(x_err); posy_err.append(y_err);
#
#     limcurr = [min(current),max(current)]
#     _plt.clf()
#     _plt.plot(current, posx, '-o', color='blue', label='Horizontal')
#     _plt.errorbar(x=current, y=posx, ecolor='black', fmt='none', yerr=posx_err)
#     _plt.plot(current, posy, '-s', color='red', label='Vertical')
#     _plt.errorbar(x=current, y=posy, ecolor='black', fmt='none', yerr=posy_err)
#     _plt.plot(limcurr,[1e6*parameters.x_misalignment_spec]*2, '--', color='blue', label='Horizontal spec')
#     _plt.plot(limcurr,[-1e6*parameters.x_misalignment_spec]*2, '--', color='blue')
#     _plt.plot(limcurr,[1e6*parameters.y_misalignment_spec]*2, '--', color='red', label='Vertical spec')
#     _plt.plot(limcurr,[-1e6*parameters.y_misalignment_spec]*2, '--', color='red')
#     _plt.legend(loc=0)
#     _plt.xlabel('current / A')
#     _plt.ylabel('Pos / um')
#     _plt.grid('on')
#     _plt.title('Residual dipolar fields for ' + magnet_name)
#     if show_fig: _plt.show()
#     if save_fig: _plt.savefig('dislocation_vs_current.png')
#
# def plot_multipole_vs_current(analysis_set, title=None, harmonic=None, relative=False, show_fig=False, save_fig=True):
#     polynom_a_avg, polynom_b_avg = [], []
#     polynom_a_std, polynom_b_std = [], []
#     current = []
#     for analysis in analysis_set:
#         parameters = analysis.parameters
#         if not harmonic:  harmonic=parameters.main_multipole_harmonic
#         data = analysis.measurement if hasattr(analysis, 'measurement') else analysis
#         current.append(data.current1_avg)
#         idx = analysis.parameters.harmonics.index(harmonic)
#         if relative:
#             polynom_a_avg.append(analysis.polynom_a_relative_avg[idx])
#             polynom_a_std.append(analysis.polynom_a_relative_std[idx])
#             polynom_b_avg.append(analysis.polynom_b_relative_avg[idx])
#             polynom_b_std.append(analysis.polynom_b_relative_std[idx])
#         else:
#             polynom_a_avg.append(analysis.polynom_a_avg[idx])
#             polynom_a_std.append(analysis.polynom_a_std[idx])
#             polynom_b_avg.append(analysis.polynom_b_avg[idx])
#             polynom_b_std.append(analysis.polynom_b_std[idx])
#     if relative:
#         if not title:
#             normal_title = data.magnet_label + '_relative_normal-{0:02d}-pole-field_vs_current'.format(2*harmonic)
#             skew_title = data.magnet_label + '_relative_skew-{0:02d}-pole-field_vs_current'.format(2*harmonic)
#         else:
#             normal_title = title+'_normal'; skew_title = title+'_skew'
#         normal_ylabel = get_multipole_name(harmonic,mult_type='normal',ref_radius=analysis.parameters.ref_radius)
#         skew_ylabel = get_multipole_name(harmonic,mult_type='skew',ref_radius=analysis.parameters.ref_radius)
#     else:
#         if not title:
#             normal_title = data.magnet_label + '_normal-{0:02d}-pole-field_vs_current'.format(2*harmonic)
#             skew_title = data.magnet_label + '_skew-{0:02d}-pole-field_vs_current'.format(2*harmonic)
#         else:
#             normal_title = title+'_normal'; skew_title = title+'_skew'
#         normal_ylabel = get_multipole_name(harmonic,mult_type='normal',ref_radius=None)
#         skew_ylabel = get_multipole_name(harmonic,mult_type='skew',ref_radius=None)
#
#     _plt.clf()
#     _plt.errorbar(current, polynom_b_avg, yerr=polynom_b_std)
#     _plt.xlabel('current / A')
#     _plt.ylabel(normal_ylabel)
#     _plt.grid('on')
#     _plt.title(normal_title)
#     try:
#         spec, _ = parameters.multipoles_spec[harmonic]
#         if spec[0]:
#             spec_sys = spec[0]
#             _plt.plot([min(current),max(current)], [spec_sys]*2, '-', color='black', label='spec_sys')
#             if spec[1]:
#                 spec_rms = [spec_sys-spec[1],spec_sys+spec[1]]
#                 _plt.plot([min(current),max(current)], [spec_rms[0]]*2, '--', color='black', label='spec_sys +/- spec_rms')
#                 _plt.plot([min(current),max(current)], [spec_rms[1]]*2, '--', color='black')
#             _plt.legend()
#         else:
#             if spec[1]:
#                 spec_rms = [-spec_rms,+spec_rms]
#                 _plt.plot([min(current),max(current)], [spec_rms[0]]*2, '--', color='black', label='+/- spec_rms')
#                 _plt.plot([min(current),max(current)], [spec_rms[1]]*2, '--', color='black')
#                 _plt.legend()
#     except KeyError:
#         pass
#
#     #_plt.gca().get_yaxis().set_major_formatter(_ticker.FormatStrFormatter('%+0.1E'))
#     if show_fig: _plt.show()
#     if save_fig: _plt.savefig(normal_title+'.png')
#
#     _plt.clf()
#     _plt.errorbar(current, polynom_a_avg, yerr=polynom_a_std)
#     _plt.xlabel('current / A')
#     _plt.ylabel(skew_ylabel)
#     _plt.grid('on')
#     _plt.title(skew_title)
#     try:
#         _, spec = parameters.multipoles_spec[harmonic]
#         if spec[0]:
#             spec_sys = spec[0]
#             _plt.plot([min(current),max(current)], [spec_sys]*2, '-', color='black', label='spec_sys')
#             if spec[1]:
#                 spec_rms = [spec_sys-spec[1],spec_sys+spec[1]]
#                 _plt.plot([min(current),max(current)], [spec_rms[0]]*2, '--', color='black', label='spec_sys +/- spec_rms')
#                 _plt.plot([min(current),max(current)], [spec_rms[1]]*2, '--', color='black')
#             _plt.legend()
#         else:
#             if spec[1]:
#                 spec_rms = [-spec_rms,+spec_rms]
#                 _plt.plot([min(current),max(current)], [spec_rms[0]]*2, '--', color='black', label='+/- spec_rms')
#                 _plt.plot([min(current),max(current)], [spec_rms[1]]*2, '--', color='black')
#                 _plt.legend()
#     except KeyError:
#         pass
#     if show_fig: _plt.show()
#     if save_fig: _plt.savefig(skew_title+'.png')
#
# def merge_measurements_with_same_current(analysis_set, magnet_label=None):
#
#     harmonics = analysis_set[0].parameters.harmonics
#
#     # partition data into a set of measurements with distinct currents
#     dset = []
#     for a in analysis_set:
#         t_avg, t_std = a.measurement.current1_avg, a.measurement.current1_std
#         matching_element = None
#         for element in dset:
#             avg, std = _numpy.mean(element[0]), max(element[1])
#             if abs(t_avg - avg) < 2*max(t_std,std):
#                 matching_element = element; break
#         if matching_element:
#             matching_element[0].append(t_avg)
#             matching_element[1].append(t_std)
#             matching_element[2].append(a)
#         else:
#             dset.append([[t_avg],[t_std],[a]])
#
#     new_analysis_set = []
#     for current_avg_set, current_std_set, analysis_set in dset:
#         analysis = AnalysisFromSet()
#         if magnet_label: analysis.magnet_label = magnet_label
#         analysis.current1_avg = _numpy.mean(current_avg_set)
#         analysis.current1_std = max(_numpy.std(current_avg_set), max(current_std_set))
#         polynom_a = _numpy.zeros((len(analysis_set), len(harmonics)))
#         polynom_b = _numpy.zeros((len(analysis_set), len(harmonics)))
#         polynom_a_relative = _numpy.zeros((len(analysis_set), len(harmonics)))
#         polynom_b_relative = _numpy.zeros((len(analysis_set), len(harmonics)))
#         for i in range(len(analysis_set)):
#             analysis.parameters = analysis_set[i].parameters
#             polynom_a[i,:] = analysis_set[i].polynom_a_avg
#             polynom_b[i,:] = analysis_set[i].polynom_b_avg
#             polynom_a_relative[i,:] = analysis_set[i].polynom_a_relative_avg
#             polynom_b_relative[i,:] = analysis_set[i].polynom_b_relative_avg
#         analysis.polynom_a_avg = _numpy.mean(polynom_a, axis=0)
#         analysis.polynom_a_std = _numpy.std(polynom_a, axis=0)
#         analysis.polynom_b_avg = _numpy.mean(polynom_b, axis=0)
#         analysis.polynom_b_std = _numpy.std(polynom_b, axis=0)
#         #sel = (max(analysis_set[i].polynom_a_std, axis=0) > analysis.polynom_a_std)
#         #analysis.polynom_a_std(sel) =
#         analysis.polynom_a_relative_avg = _numpy.mean(polynom_a_relative, axis=0)
#         analysis.polynom_a_relative_std = _numpy.std(polynom_a_relative, axis=0)
#         analysis.polynom_b_relative_avg = _numpy.mean(polynom_b_relative, axis=0)
#         analysis.polynom_b_relative_std = _numpy.std(polynom_b_relative, axis=0)
#         for i in range(len(analysis_set)):
#             # checks if intrinsic errors are larger than dispersion amongst data
#             sel = analysis_set[i].polynom_a_std > analysis.polynom_a_std
#             analysis.polynom_a_std[sel] = analysis_set[i].polynom_a_std[sel]
#             sel = analysis_set[i].polynom_a_relative_std > analysis.polynom_a_relative_std
#             analysis.polynom_a_relative_std[sel] = analysis_set[i].polynom_a_relative_std[sel]
#             sel = analysis_set[i].polynom_b_std > analysis.polynom_b_std
#             analysis.polynom_b_std[sel] = analysis_set[i].polynom_b_std[sel]
#             sel = analysis_set[i].polynom_b_relative_std > analysis.polynom_b_relative_std
#             analysis.polynom_b_relative_std[sel] = analysis_set[i].polynom_b_relative_std[sel]
#
#         new_analysis_set.append(analysis)
#     return new_analysis_set
#
# def run_analysis(local_parameters, analysis_type='excitation'):
#
#     parameters = local_parameters.parameters
#
#     # loads default parameters for analysis
#     print('Analysis parameters: {0:s}'.format(parameters.label))
#
#     # gets list of files (excluding I=0A data)
#     top_folder = _os.path.join(_os.getcwd(), local_parameters.data_folder)
#     all_files = _utils.files_get_matches(folder=top_folder, recursive=True, strs_in=['_BOB_','.dat'], strs_out='+0000A')
#     print('Analysis of files in {1:s}...'.format(len(all_files), top_folder))
#
#     # analyzes all data
#     analysis_set = []
#     for fname in all_files:
#         print('  ' + _os.path.basename(fname) + '...')
#         m = Measurement(fname)
#         if local_parameters.only_ramp_up and analysis_set and m.current1_avg < analysis_set[-1].measurement.current1_avg: break
#         a = Analysis(m, parameters)
#         txt_fname = _os.path.basename(fname).replace('.dat','_analysis.txt')
#         with open(txt_fname,'w') as tf:
#             tf.write(a.__str__())
#         analysis_set.append(a)
#         if local_parameters.plot_each_data_file: plot_harmonic_analysis(a,
#             ylim=local_parameters.mult_ylim,
#             show_fig=local_parameters.show_fig,
#             save_fig=local_parameters.save_fig)
#
#     magnet_label = analysis_set[0].measurement.magnet_label
#     harmonics = analysis_set[0].parameters.harmonics
#
#     if analysis_type == 'excitation':
#         # merge different measurements at the same current
#         print('Merge measurement data with the same current')
#         analysis_set = merge_measurements_with_same_current(analysis_set, magnet_label=magnet_label)
#
#     # plots excitation curve
#     print('Plot of excitation curve')
#     plot_multipole_vs_current(analysis_set,
#         title=magnet_label + '_excitation_curve',
#         harmonic=None,
#         relative=False,
#         show_fig=local_parameters.show_fig,
#         save_fig=local_parameters.save_fig)
#
#     print('Plots magnet off-center dislocation')
#     plot_dislocation_vs_current(analysis_set,
#                                 show_fig=local_parameters.show_fig,
#                                 save_fig=local_parameters.save_fig)
#
#     print('Plots magnet equivalent roll rotation')
#     plot_skew_rotation_vs_current(analysis_set,
#                                   show_fig=local_parameters.show_fig,
#                                   save_fig=local_parameters.save_fig)
#
#     # plots all relative multipoles as function of current
#     print('Plot of current-dependent relative multipoles')
#     for h in harmonics:
#         plot_multipole_vs_current(analysis_set,
#             title=None,
#             harmonic=h,
#             relative=True,
#             show_fig=local_parameters.show_fig,
#             save_fig=local_parameters.save_fig)
#
#     if local_parameters.move_files:
#         print('Moving analysis file to analysis folder')
#         # moves file to analysis folder
#         plot_files = _utils.files_get_matches(folder=_os.getcwd(), recursive=False, strs_in='.png', strs_out=None)
#         txt_files = _utils.files_get_matches(folder=_os.getcwd(), recursive=False, strs_in='analysis.txt', strs_out=None)
#         analysis_folder = _os.path.join(top_folder, 'analysis')
#         if not _os.path.exists(analysis_folder): _os.makedirs(analysis_folder)
#         for plot_file in plot_files + txt_files:
#             new_path = _os.path.join(top_folder, 'analysis',_os.path.basename(plot_file))
#             _os.rename(plot_file, new_path)
#
# def run_excitation_curve_analysis(local_parameters):
#
#     local_parameters.data_folder = local_parameters.excitation_data_folder
#     run_analysis(local_parameters, analysis_type='excitation')
#
# def run_histeresis_curve_analysis(local_parameters):
#
#     local_parameters.data_folder = local_parameters.histeresis_data_folder
#     run_analysis(local_parameters, analysis_type='histeresis')
