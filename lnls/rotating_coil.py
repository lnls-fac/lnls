#!/usr/bin/env python3

#import lnls
#import os
#import shutil
import numpy as _numpy
import math as _math
#import matplotlib.pyplot as plt
#import re


class BQAnalysisParameters:
    def __init__(self):
        self.main_multipole_harmonic = 2  # [1: dipole, 2:quadrupole, ...]
        self.main_multipole_is_skew  = False
        self.harmonics = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] # [1: dipole, 2:quadrupole, ...]
        self.multipoles_ylim = [1e-6,1e1]
        self.ref_radius = 0.0175 # [m]
        self.multipoles_spec = {
            # (normal_sys, normal_std) (skew_sys, skew_std)
            1: ((None,None),    (None,None)),
            2: ((None,None),    (None,None)),
            3: ((None,None),    (+7.0e-4,1e-3)),
            4: ((None,None),    (+4.0e-4,5e-4)),
            5: ((None,None),    (+4.0e-4,1e-4)),
            6: ((-1.0e-3,None), (+4.0e-4,1e-4)),
            7: ((None,None),    (+4.0e-4,1e-4)),
            8: ((None,None),    (+4.0e-4,1e-4)),
            9: ((None,None),    (+4.0e-4,1e-4)),
            10:((+1.1e-3,None), (None,None)),
            14:((+8.0e-5,None), (None,None)),
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


    def __str__(self):
        r = ''
        r +=   '{0:<30s} {1:s}'.format('data filename', self.filename)
        r += '\n{0:<30s} {1:s}'.format('timestamp', self.timestamp)
        r += '\n{0:<30s} {1:+.3f}'.format('main_current_avg[A]', self.current1_avg)
        r += '\n{0:<30s} {1:.3f}'.format('main_current_std[A]', self.current1_std)
        r += '\n{0:<30s} {1:s}'.format('temperature', self.temperature)
        r += '\n{0:<30s} {1:d}'.format('nr_turns', self.nr_turns)


        r += '\n\n--- rotating coil ---\n\n'
        r += self.rotating_coil.__str__()
        return r


class Analysis:

    def __init__(self, measurement, parameters, run_analysis=True):
        self.measurement = measurement
        self.parameters = parameters
        if run_analysis: self.run_analysis()

    def run_analysis(self):
        if self.measurement.rotating_coil.type == 'Bobina Radial':
            self._calc_multipoles_for_radial_rotcoil()
            self._calc_normalized_multipoles()
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


    def _calc_normalized_multipoles(self):

        ref_radius = self.parameters.ref_radius
        main_harmonic = self.parameters.main_multipole_harmonic
        main_is_skew = self.parameters.main_multipole_is_skew
        harmonics = self.parameters.harmonics

        idx = harmonics.index(main_harmonic)
        if main_is_skew:
            main_multipole_at_ref_radius = self.polynom_a[:,idx] * (ref_radius ** (main_harmonic-1))
        else:
            main_multipole_at_ref_radius = self.polynom_b[:,idx] * (ref_radius ** (main_harmonic-1))
        self.polynom_a_normalized = 0 * self.polynom_a
        self.polynom_b_normalized = 0 * self.polynom_b
        for i in range(len(harmonics)):
            polynom_a_at_rotcoil_radius = self.polynom_a[:,i] * (ref_radius ** (harmonics[i]-1))
            polynom_b_at_rotcoil_radius = self.polynom_b[:,i] * (ref_radius ** (harmonics[i]-1))
            self.polynom_a_normalized[:,i] = polynom_a_at_rotcoil_radius / main_multipole_at_ref_radius
            self.polynom_b_normalized[:,i] = polynom_b_at_rotcoil_radius / main_multipole_at_ref_radius




#
# def load_data_from_file(file_name, rotcoil_factor, plot=True):
#
#     _parameters = {
#         # rotating coil parameters
#         'n_espiras_bobina_principal':'nr_coil_turns',
#         'raio_interno_bobina_princip(m)':'r1',
#         'raio_externo_bobina_princip(m)':'r2',
#         # measurement parameters
#         'corrente_alim_principal_avg(A)':'current1_avg',
#         'corrente_alim_principal_std(A)':'current1_std',
#         'corrente_alim_secundaria_avg(A)':'current2_avg',
#         'corrente_alim_secundaria_std(A)':'current2_std',
#         'nr_pontos_integracao':'nr_pts',
#         'nr_voltas':'nr_turns',
#         }
#
#     # reads raw data from file
#     lines = [line.strip() for line in open(file_name, encoding='latin-1')]
#
#     parameters = dict()
#     for i in range(len(lines)):
#         line = lines[i]
#         if 'Volta_1' in line:
#             raw_data_txt = lines[i+1:]
#             break
#         if not line: continue
#         if not line or line[0] == '#': continue
#         for parameter in _parameters.keys():
#             if parameter in line:
#                 data = str.replace(line, parameter, '').strip()
#                 parameters[_parameters[parameter]] = eval(data) #float(data)
#                 #print(_parameters[parameter] + ': ' + data)
#
#     raw_data = numpy.zeros((parameters['nr_turns'],parameters['nr_pts']))
#     idx = 0
#     for line in raw_data_txt:
#         raw_data[:,idx] = [rotcoil_factor*float(word) for word in line.split(' ')]
#         idx += 1
#
#     # subtracts linear drift
#     linear_data = numpy.reshape(raw_data,(1,-1))
#     raw_data = numpy.cumsum(numpy.roll(linear_data, 0))
#     avg_raw_data = numpy.mean(linear_data)
#     linear_data_corrected = linear_data - avg_raw_data
#     raw_data_corrected = numpy.cumsum(numpy.roll(linear_data_corrected, 0))
#
#     if plot:
#         plt.plot(raw_data)
#         plt.plot(raw_data_corrected)
#         plt.show()
#     raw_data_corrected = numpy.reshape(raw_data_corrected, (parameters['nr_turns'],-1))
#     parameters['flux'] = raw_data_corrected
#     return parameters
#
# def get_all_data_files(folder=None, recursive=True, sub_strs=None):
#     if folder is None:
#         folder = os.getcwd()
#     if sub_strs is None:
#         sub_strs = ('.dat','BOB_')
#     elif isinstance(tokens, str):
#         tokens = (tokens,)
#     files = []
#     local_files = os.listdir(folder)
#     for fname in local_files:
#         path = os.path.join(folder, fname)
#         if os.path.isdir(path):
#             files.extend(get_all_data_files(path, recursive, sub_strs))
#         else:
#             is_substr = [1 if token in path else 0 for token in sub_strs]
#             if sum(is_substr): files.append(path)
#     return files
#
# def select_data_files(all_files, magnet_str, current_str, non_permitted_str='NON_PERMITTED'):
#     files = []
#     for filename in all_files:
#         if magnet_str in filename and current_str in filename and non_permitted_str not in filename:
#             files.append(filename)
#     return files
#
# def calc_multipoles_from_radial_coil_flux(flux, harmonics, r1, r2, nr_coil_turns):
#     polynom_a, polynom_b = [],[]
#     angle = numpy.linspace(2*math.pi/len(flux),2*math.pi,len(flux))
#     for harmonic in harmonics:
#         #n = harmonic + 1
#         n = harmonic
#         integ_cos = numpy.trapz(flux * numpy.cos(n*angle), angle)
#         integ_sin = numpy.trapz(flux * numpy.sin(n*angle), angle)
#         #if (n == 3):
#         #    print(integ_cos, integ_sin)
#         b_n = +n * integ_cos / math.pi / nr_coil_turns / (r2**n - r1**n)
#         a_n = -n * integ_sin / math.pi / nr_coil_turns / (r2**n - r1**n)
#         polynom_b.append(b_n)
#         polynom_a.append(a_n)
#     return polynom_a, polynom_b
#
# def calc_multipoles_from_radial_coil_data(parameters, plot=True):
#     flux = parameters['flux']
#     r1, r2 = parameters['r1'], parameters['r2']
#     nr_coil_turns = int(parameters['nr_coil_turns'])
#     nr_turns = int(parameters['nr_turns'])
#     harmonics = parameters['harmonics']
#     polynom_a, polynom_b = numpy.zeros((nr_turns, len(harmonics))), numpy.zeros((nr_turns, len(harmonics)))
#     for i in range(nr_turns):
#          polynom_a[i,:], polynom_b[i,:] = calc_multipoles_from_radial_coil_flux(flux[i,:], harmonics, r1, r2, nr_coil_turns)
#     parameters['polynom_a'] = polynom_a
#     parameters['polynom_b'] = polynom_b
#     parameters['polynom'] = (polynom_a**2 + polynom_b**2)**0.5
#     if plot:
#         for f in flux:
#             plt.plot(f)
#         plt.show()
#
# def run_radial_coil_analysis(file_name, harmonics, main_harmonic, rotcoil_radius, rotcoil_factor, plot_flag=True, print_flag=True):
#
#     # reads data from file
#     parameters  = load_data_from_file(file_name, rotcoil_factor, plot=plot_flag)
#
#     parameters['file_name'] = file_name
#
#     # calcs multipoles
#     parameters['harmonics'] = harmonics
#     calc_multipoles_from_radial_coil_data(parameters, plot=plot_flag)
#
#     # calcs normalized multipoles
#     calc_normalized_multipoles(main_harmonic=main_harmonic,
#                                parameters=parameters,
#                                rotcoil_radius=rotcoil_radius)
#     # calcs stats
#     calc_stats_unique_excitation(parameters)
#
#     # prints data
#     if print_flag:
#         print_analysis(parameters)
#
#     return parameters
#
# def calc_stats_unique_excitation(parameters):
#     try:
#         polynom_a = parameters['polynom_a']
#         polynom_b = parameters['polynom_b']
#         polynom   = parameters['polynom']
#         parameters['avg_polynom_b'] = numpy.mean(polynom_b,axis=0)
#         parameters['avg_polynom_a'] = numpy.mean(polynom_a,axis=0)
#         parameters['std_polynom_b'] = numpy.std(polynom_b,axis=0)
#         parameters['std_polynom_a'] = numpy.std(polynom_a,axis=0)
#         parameters['avg_polynom']   = numpy.mean(polynom, axis=0)
#         parameters['std_polynom']   = numpy.std(polynom, axis=0)
#     except KeyError:
#         pass
#     try:
#         norm_polynom_a = parameters['norm_polynom_a']
#         norm_polynom_b = parameters['norm_polynom_b']
#         norm_polynom   = parameters['norm_polynom']
#         parameters['norm_avg_polynom_b'] = numpy.mean(norm_polynom_b,axis=0)
#         parameters['norm_avg_polynom_a'] = numpy.mean(norm_polynom_a,axis=0)
#         parameters['norm_std_polynom_b'] = numpy.std(norm_polynom_b,axis=0)
#         parameters['norm_std_polynom_a'] = numpy.std(norm_polynom_a,axis=0)
#         parameters['norm_avg_polynom']   = numpy.mean(norm_polynom, axis=0)
#         parameters['norm_std_polynom']   = numpy.std(norm_polynom, axis=0)
#     except KeyError:
#         pass
#
# def calc_stats_variable_excitation(analysis_list, label='sample_average'):
#
#     norm_polynom_a, norm_polynom_b, norm_polynom = [], [], []
#     n = len(analysis_list)
#     for i in range(n):
#         analysis = analysis_list[i]
#         norm_polynom_a.extend(list(analysis['norm_polynom_a']))
#         norm_polynom_b.extend(list(analysis['norm_polynom_b']))
#         norm_polynom.extend(list(analysis['norm_polynom']))
#     parameters = {}
#     parameters['harmonics'] = analysis['harmonics']
#     parameters['current1_avg'] = float('nan')
#     parameters['file_name'] = 'N.D.'
#     parameters['rotcoil_radius'] = analysis['rotcoil_radius']
#     parameters['norm_polynom_a'] = numpy.array(norm_polynom_a)
#     parameters['norm_polynom_b'] = numpy.array(norm_polynom_b)
#     parameters['norm_polynom']   = numpy.array(norm_polynom)
#
#     #d = [norm_polynom_b[0] for norm_polynom_b in parameters['norm_polynom_b']]
#     #print(d)
#     calc_stats_unique_excitation(parameters)
#     return parameters
#
# def merge_unique_excitation_data(analysis_list):
#     polynom_a, polynom_b = [], []
#     n = len(analysis_list)
#     for i in range(n):
#         analysis = analysis_list[i]
#         polynom_a.extend(list(analysis['polynom_a']))
#         polynom_b.extend(list(analysis['polynom_b']))
#     parameters = {}
#     parameters.update(analysis_list[0])
#     parameters['file_name'] = 'MERGED_DATA'
#     parameters['polynom_a'] = numpy.array(polynom_a)
#     parameters['polynom_b'] = numpy.array(polynom_b)
#     calc_normalized_multipoles(main_harmonic=parameters['main_harmonic'],
#                                parameters=parameters,
#                                rotcoil_radius=parameters['rotcoil_radius'])
#     calc_stats_unique_excitation(parameters)
#     return parameters
#
# def calc_normalized_multipoles(main_harmonic,parameters,rotcoil_radius):
#
#     parameters['rotcoil_radius'] = rotcoil_radius
#     parameters['main_harmonic'] = main_harmonic
#     harmonics = parameters['harmonics']
#     polynom_a = parameters['polynom_a']
#     polynom_b = parameters['polynom_b']
#     idx = harmonics.index(main_harmonic[0])
#     if main_harmonic[1]:
#         main_multipole_at_rotcoil_radius = polynom_a[:,idx] * (rotcoil_radius ** (main_harmonic[0]-1))
#     else:
#         main_multipole_at_rotcoil_radius = polynom_b[:,idx] * (rotcoil_radius ** (main_harmonic[0]-1))
#     norm_polynom_a = 0 * polynom_a
#     norm_polynom_b = 0 * polynom_b
#     for i in range(len(harmonics)):
#         # polynom_a_at_rotcoil_radius = polynom_a[:,i] * (rotcoil_radius ** harmonics[i])
#         # polynom_b_at_rotcoil_radius = polynom_b[:,i] * (rotcoil_radius ** harmonics[i])
#         polynom_a_at_rotcoil_radius = polynom_a[:,i] * (rotcoil_radius ** (harmonics[i]-1))
#         polynom_b_at_rotcoil_radius = polynom_b[:,i] * (rotcoil_radius ** (harmonics[i]-1))
#         norm_polynom_a[:,i] = polynom_a_at_rotcoil_radius / main_multipole_at_rotcoil_radius
#         norm_polynom_b[:,i] = polynom_b_at_rotcoil_radius / main_multipole_at_rotcoil_radius
#     parameters['norm_polynom_a'] = norm_polynom_a
#     parameters['norm_polynom_b'] = norm_polynom_b
#     parameters['norm_polynom'] = (norm_polynom_a**2 + norm_polynom_b**2)**0.5
#
# def print_analysis(parameters, normalized=False):
#     harmonics = parameters['harmonics']
#     #polynom_a = parameters['polynom_a']
#     #polynom_b = parameters['polynom_b']
#     if normalized:
#         try:
#             avg_polynom_a = parameters['norm_avg_polynom_a']
#             avg_polynom_b = parameters['norm_avg_polynom_b']
#             std_polynom_a = parameters['norm_std_polynom_a']
#             std_polynom_b = parameters['norm_std_polynom_b']
#             avg_polynom   = parameters['norm_avg_polynom']
#             std_polynom   = parameters['norm_std_polynom']
#         except KeyError:
#             print('normalized data missing!')
#             return
#     else:
#         try:
#             avg_polynom_a = parameters['avg_polynom_a']
#             avg_polynom_b = parameters['avg_polynom_b']
#             std_polynom_a = parameters['std_polynom_a']
#             std_polynom_b = parameters['std_polynom_b']
#             avg_polynom   = parameters['avg_polynom']
#             std_polynom   = parameters['std_polynom']
#         except KeyError:
#             print('non-normalized data missing!')
#             return
#     header = '{0:2s} | {1:11s} {2:10s} | {3:11s} {4:10s} | {5:10s} {6:10s}'.format('n','AvgPolyB','StdPolyB','AvgPolyA','StdPolyA', 'AvgAbsPoly','StdAbsPoly')
#     print(header)
#     print('-'*len(header))
#     for i in range(len(harmonics)):
#         print('{0:02d} | {1:+.4e} {2:.4e} | {3:+.4e} {4:.4e} | {5:.4e} {6:.4e}'.format(1+harmonics[i], avg_polynom_b[i], std_polynom_b[i], avg_polynom_a[i], std_polynom_a[i], avg_polynom[i], std_polynom[i]))
#     print('-'*len(header))
#
# def plot_excitation_curve(analysis_list,
#                           harmonic,
#                           plot_title = '',
#                           poly_type = 'polynom_b',
#                           normalized=False,
#                           current_flag=True,
#                           show_plot=False,
#                           save_figs=False,
#                           save_fname=None,
#                           spec=(None,None)
#                           ):
#     if not analysis_list: return
#     plt.clf()
#     if normalized:
#         rotcoil_radius = analysis_list[0]['rotcoil_radius']
#         multipole_labels = [
#             'normalized dipolar field at $r_0$ = {0:.1f} mm'.format(1000*rotcoil_radius),
#             'normalized quadrupolar field at $r_0$ = {0:.1f} mm'.format(1000*rotcoil_radius),
#             'normalized sextupolar field at $r_0$ = {0:.1f} mm'.format(1000*rotcoil_radius),
#             'normalized octupolar field at $r_0$ = {0:.1f} mm'.format(1000*rotcoil_radius),
#             'normalized decapolar field at $r_0$ = {0:.1f} mm'.format(1000*rotcoil_radius),
#             'normalized duodecapolar field at $r_0$ = {0:.1f} mm'.format(1000*rotcoil_radius),
#         ]
#         default_label = 'normalized multipolar (n={0:d}) field at $r_0$ = {1:.1f} mm'.format(harmonic, 1000*rotcoil_radius)
#     else:
#         multipole_labels = [
#             'integrated dipolar field component [T.m]',
#             'integrated quadrupolar field component [T]',
#             'integrated sextupolar field component [T/m]',
#             'integrated octupolar field component [T/m^2]',
#             'integrated decapolar field component [T/m^3]',
#             'integrated duodecapolar field component [T/m^4]',
#         ]
#         default_label = 'integrated (n='+str(harmonic)+') multipolar component [T/m^{0}]'.format(harmonic-2)
#     current, multipole_avg, multipole_std = [], [], []
#     for i in range(len(analysis_list)):
#         analysis = analysis_list[i]
#         this_current = analysis['current1_avg'] if current_flag else i
#         current.append(this_current)
#         idx = analysis['harmonics'].index(harmonic)
#         if normalized:
#             avg_polynom = analysis['norm_avg_' + poly_type]
#             std_polynom = analysis['norm_std_' + poly_type]
#         else:
#             avg_polynom = analysis['avg_' + poly_type]
#             std_polynom = analysis['std_' + poly_type]
#             spec = (None,None)
#         multipole_avg.append(avg_polynom[idx])
#         multipole_std.append(std_polynom[idx])
#     try:
#         ylabel = multipole_labels[harmonic-1]
#     except:
#         ylabel = default_label
#
#     #colors = ['red' if value < 0 else 'blue' for value in multipole_avg]
#     #plt.bar(current, numpy.array(multipole_avg), yerr=multipole_std, log=False, color=colors, error_kw=dict(elinewidth=2,ecolor='black'), align='center')
#     plt.errorbar(current, multipole_avg, yerr=multipole_std)
#     if not spec[0]:
#         if spec[1]:
#             plt.plot(current, -spec[1]*len(current), color='black')
#             plt.plot(current, +spec[1]*len(current), color='black')
#     else:
#         if not spec[1]:
#             plt.plot(current, [spec[0]]*len(current), color='black')
#         else:
#             plt.plot(current, [spec[0]-spec[1]]*len(current), color='black')
#             plt.plot(current, [spec[0]+spec[1]]*len(current), color='black')
#     plt.ylabel(ylabel)
#     plt.grid('on')
#     plt.title(plot_title)
#     if current_flag:
#         plt.xlabel('current [A]')
#     else:
#         plt.xlabel('index')
#     if show_plot:plt.show()
#     if save_figs: plt.savefig(save_fname)
#
# def plot_multipoles(parameters,
#                     poly_type='polynom_b',
#                     title='',
#                     ylim=None,
#                     normalized=False,
#                     show_figs=True,
#                     save_fname=None):
#     plt.clf()
#     harmonics  = parameters['harmonics']
#     if normalized:
#         avg_polynom = parameters['norm_avg_'+poly_type]
#         std_polynom = parameters['norm_std_'+poly_type]
#         rotcoil_radius = parameters['rotcoil_radius']
#         ylabel = 'normalized multipoles @ $r_0$ = {0:.1f} mm'.format(1000*rotcoil_radius)
#     else:
#         avg_polynom = parameters['avg_'+poly_type]
#         std_polynom = parameters['std_'+poly_type]
#         ylabel = 'absolute multipoles'
#     # if poly_type == 'polynom_b':
#     #     title = 'Normal Multipoles'
#     # else:
#     #     title = 'Skew Multipoles'
#
#     colors = ['red' if value < 0 else 'blue' for value in avg_polynom]
#     plt.bar(harmonics, abs(avg_polynom), yerr=std_polynom, log=True, color=colors, error_kw=dict(elinewidth=2,ecolor='black'), align='center')
#     plt.xlabel('harmonics')
#     plt.ylabel(ylabel)
#     plt.grid('on')
#     plt.title(title)
#     if ylim is not None: plt.ylim(ylim)
#     if show_figs: plt.show()
#     if save_fname: plt.savefig(save_fname)
#
# def run_analysis(files, analysis_parms, variable_excitation=True):
#     if not files: return None
#     alist = []
#     for file in files:
#         if analysis_parms['rotcoil_type'] == 'radial':
#             analysis = run_radial_coil_analysis(file,
#                                                 harmonics=analysis_parms['harmonics'],
#                                                 main_harmonic=analysis_parms['main_harmonic'],
#                                                 rotcoil_radius=analysis_parms['rotcoil_radius'],
#                                                 rotcoil_factor=analysis_parms['rotcoil_factor'],
#                                                 plot_flag=False,
#                                                 print_flag=False)
#         else:
#             Exception('analysis of this kind of rotating coil not implemented!')
#
#         alist.append(analysis)
#     if variable_excitation:
#         avg = [calc_stats_variable_excitation(alist)]
#     else:
#         parameters = merge_unique_excitation_data(alist)
#         calc_stats_unique_excitation(parameters)
#         avg = [parameters]
#     return avg
#
# def save_multipoles_fig(analysis_list, label1='sample_average', label2='', ylim=None):
#     if ylim is None:
#         ylim=[1e-5,1e2]
#     n = len(analysis_list)
#     for i in range(n):
#         analysis = analysis_list[i]
#         try:
#             current = round(analysis['current1_avg'])
#         except:
#             label2 = ''
#         if label2 == '':
#             plot_title_a = 'Skew multipoles for '+label1
#             plot_title_b = 'Normal multipoles for '+label1
#             fname_a = label1+'_harmonic_skew_'+label2+'_{0:02d}.png'.format(i+1) if n > 1 else label1+'_harmonics_skew.png'
#             fname_b = label1+'_harmonics_normal_'+label2+'_{0:02d}.png'.format(i+1) if n > 1 else label1+'_harmonics_normal.png'
#         else:
#             plot_title_a = 'Skew multipoles for '+label1+', I={0}A'.format(current)
#             plot_title_b = 'Normal multipoles for '+label1+', I={0}A'.format(current)
#             fname_a = label1+'_harmonics_skew_'+label2+'_{0:02d}.png'.format(i+1) if n > 1 else label1+'_harmonics_skew_'+label2+'.png'
#             fname_b = label1+'_harmonics_normal_'+label2+'_{0:02d}.png'.format(i+1) if n > 1 else label1+'_harmonics_normal_'+label2+'.png'
#         plot_multipoles(analysis,
#                         title=plot_title_a,
#                         ylim=ylim,
#                         normalized=True,
#                         show_figs=False,
#                         poly_type='polynom_a',
#                         save_fname=fname_a)
#         plot_multipoles(analysis,
#                         title=plot_title_b,
#                         ylim=ylim,
#                         normalized=True,
#                         show_figs=False,
#                         poly_type='polynom_b',
#                         save_fname=fname_b)
#
# def run_analysis_separate_excitations(parameters, magnet_name, all_files):
#     print('<<< analysis separate excitations for '+magnet_name+' >>>')
#     current_values = parameters['current_values']
#     label1 = magnet_name
#     data = {}
#     for i in range(len(current_values)):
#         key = '{0:s}'.format(current_values[i])
#         data[key] = select_data_files(all_files, magnet_name, current_values[i])
#     for current_label in data.keys():
#         if not data[current_label]: continue
#         string = 'current value '+current_label
#         print(string); #print('='*len(string))
#         avg = run_analysis(data[current_label], parameters, variable_excitation=False)
#         save_multipoles_fig(avg, label1=label1, label2=current_label, ylim=parameters['multipoles_ylim'])
#         #print('absolute multipoles:'); print_analysis(avg[0], normalized=False)
#         #print('normalized multipoles:'); print_analysis(avg[0], normalized=True)
#         #print('')
#
# def run_analysis_all_excitations(parameters, magnet_name, all_files):
#     print('<<< analysis all excitations for '+magnet_name+' >>>')
#     files = select_data_files(all_files, magnet_name, '', '+0000A')
#     avg = run_analysis(files, parameters, variable_excitation=True)
#     save_multipoles_fig(avg, label1=magnet_name, label2='', ylim=parameters['multipoles_ylim'])
#     print('absolute multipoles:'); print_analysis(avg[0], normalized=False)
#     print('normalized multipoles:'); print_analysis(avg[0], normalized=True)
#
# def run_analysis_excitation_curve(analysis_parms, all_files):
#     for magnet in analysis_parms['magnets']:
#         print('<<< analysis separate excitations curve for '+magnet+' >>>')
#         excitation = []
#         for current_name in analysis_parms['current_values']:
#             files = select_data_files(all_files, magnet, current_name)
#             if not files: continue
#             avg = run_analysis(files, analysis_parms, variable_excitation=False)
#             excitation.append(avg[0])
#         for harmonic in analysis_parms['harmonics']:
#             print('magnet {0:s}, harmonic {1:02d}'.format(magnet, harmonic))
#             try:
#                 spec = analysis_parms['multipoles_spec'][harmonic]
#             except:
#                 spec = ((None,None),)*2
#             plot_excitation_curve(excitation,
#                                   harmonic=harmonic,
#                                   plot_title='Normal multipoles for ' + magnet,
#                                   normalized=False,
#                                   show_plot=False,
#                                   save_figs=True,
#                                   poly_type='polynom_b',
#                                   spec=spec[0],
#                                   save_fname=magnet+'_excitation_normal_absolute_h={0:02d}'.format(harmonic))
#             plot_excitation_curve(excitation,
#                                  harmonic=harmonic,
#                                  plot_title='Normal multipoles for ' + magnet,
#                                  normalized=True,
#                                  show_plot=False,
#                                  save_figs=True,
#                                  poly_type='polynom_b',
#                                  spec=spec[0],
#                                  save_fname=magnet+'_excitation_normal_relative_h={0:02d}'.format(harmonic))
#             plot_excitation_curve(excitation,
#                                 harmonic=harmonic,
#                                 plot_title='Skew multipoles for ' + magnet,
#                                 normalized=False,
#                                 show_plot=False,
#                                 save_figs=True,
#                                 poly_type='polynom_a',
#                                 spec=spec[1],
#                                 save_fname=magnet+'_excitation_skew_absolute_h={0:02d}'.format(harmonic))
#             plot_excitation_curve(excitation,
#                                harmonic=harmonic,
#                                plot_title='Skew multipoles for ' + magnet,
#                                normalized=True,
#                                show_plot=False,
#                                save_figs=True,
#                                poly_type='polynom_a',
#                                spec=spec[1],
#                                save_fname=magnet+'_excitation_skew_relative_h={0:02d}'.format(harmonic))
#
# def store_analysis_files(analysis_parms):
#     files = os.listdir(analysis_parms['top_folder'])
#     for file in files:
#         if 'png' in file:
#             old_path = os.path.join(analysis_parms['top_folder'],file)
#             new_dir  = os.path.join(analysis_parms['top_folder'],analysis_parms['analysis_subfolder'])
#             if not os.path.isdir(new_dir):
#                 os.makedirs(new_dir)
#             new_path = os.path.join(new_dir,file)
#             try:
#                 os.remove(new_path)
#             except:
#                 pass
#             shutil.move(old_path, new_path)
#
# def run_magnet_analysis(magnet_object, magnet_name):
#
#     # --- Loads default global analysis parameters
#     parameters = dict(magnet_object.__dict__)
#
#     # --- Sets specifics analysis parameters
#     parameters['magnet']             = magnet_name
#     parameters['magnets']            = [parameters['magnet']]
#     parameters['top_folder']         = os.getcwd()
#     parameters['analysis_subfolder'] = 'analysis'
#
#     # --- Gets list of all data file within top_folder and subfolders
#     all_files = get_all_data_files(parameters['top_folder'])
#
#     # --- Analyzes multipoles for each separate excitation current
#     run_analysis_separate_excitations(parameters, parameters['magnet'], all_files)
#
#     # --- Analyzes multipole as a function of excitation curve for each magnet
#     run_analysis_excitation_curve(parameters, all_files)
#
#     # --- Analyzes multipoles for all currents (excluding I=0A)
#     run_analysis_all_excitations(parameters, parameters['magnet'], all_files)
#
#     # --- Stores analysis file into appropriate directories
#     store_analysis_files(parameters)
#
# def process_analysis_cmd(argv, magnet_object, magnet_name):
#
#     def run(magnet_object, magnet_name):
#         run_magnet_analysis(magnet_object, magnet_name)
#
#     def clean():
#         all_files = get_all_data_files(folder=os.getcwd(), recursive=True, tokens='.png')
#         for path in all_files:
#             print('removing ' + path + '...')
#             os.remove(path)
#
#     # processes input arguments
#     if len(argv) == 2:
#         cmd = argv[1].lower()
#         if cmd == 'clean':
#             clean()
#         elif cmd == 'run':
#             run(magnet_object, magnet_name)
#         else:
#             print(argv[0] + ': invalid arguments')
#     else:
#         run(magnet_object, magnet_name)
