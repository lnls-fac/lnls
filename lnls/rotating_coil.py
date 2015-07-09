#!/usr/bin/env python3

import lnls
import os
import numpy
import math
import matplotlib.pyplot as plt

def load_data_from_file(file_name, conversion, plot=True):

    _parameters = {
        'corrente_alim_principal_avg(A)':'current1_avg',
        'corrente_alim_principal_std(A)':'current1_std',
        'corrente_alim_secundaria_avg(A)':'current2_avg',
        'corrente_alim_secundaria_std(A)':'current2_std',
        'n_espiras_bobina_principal':'nr_coil_turns',
        'raio_interno_bobina_princip(m)':'r1',
        'raio_externo_bobina_princip(m)':'r2',
        'nr_pontos_integracao':'nr_pts',
        'nr_voltas':'nr_turns',
        }

    # reads raw data from file
    lines = [line.strip() for line in open(file_name, encoding='latin-1')]
    parameters = dict()
    for i in range(len(lines)):
        line = lines[i]
        if 'Volta_1' in line:
            raw_data_txt = lines[i+1:]
            break
        if not line: continue
        if not line or line[0] == '#': continue
        for parameter in _parameters.keys():
            if parameter in line:
                data = str.replace(line, parameter, '').strip()
                parameters[_parameters[parameter]] = eval(data) #float(data)
                #print(_parameters[parameter] + ': ' + data)

    raw_data = numpy.zeros((parameters['nr_turns'],parameters['nr_pts']))
    idx = 0
    for line in raw_data_txt:
        raw_data[:,idx] = [conversion*float(word) for word in line.split(' ')]
        idx += 1

    # subtracts linear drift
    linear_data = numpy.reshape(raw_data,(1,-1))
    raw_data = numpy.cumsum(numpy.roll(linear_data, 0))
    avg_raw_data = numpy.mean(linear_data)
    linear_data_corrected = linear_data - avg_raw_data
    raw_data_corrected = numpy.cumsum(numpy.roll(linear_data_corrected, 0))

    if plot:
        plt.plot(raw_data)
        plt.plot(raw_data_corrected)
        plt.show()
    raw_data_corrected = numpy.reshape(raw_data_corrected, (parameters['nr_turns'],-1))
    parameters['flux'] = raw_data_corrected
    return parameters

def calc_multipoles_from_radial_coil_flux(flux, monomials, r1, r2, nr_coil_turns):
    polynom_a, polynom_b = [],[]
    angle = numpy.linspace(2*math.pi/len(flux),2*math.pi,len(flux))
    for monomial in monomials:
        n = monomial + 1
        integ_cos = numpy.trapz(flux * numpy.cos(n*angle), angle)
        integ_sin = numpy.trapz(flux * numpy.sin(n*angle), angle)
        #if (n == 3):
        #    print(integ_cos, integ_sin)
        b_n = +n * integ_cos / math.pi / nr_coil_turns / (r2**n - r1**n)
        a_n = -n * integ_sin / math.pi / nr_coil_turns / (r2**n - r1**n)
        polynom_b.append(b_n)
        polynom_a.append(a_n)
    return polynom_a, polynom_b

def calc_multipoles_from_radial_coil_data(parameters, plot=True):
    flux = parameters['flux']
    r1, r2 = parameters['r1'], parameters['r2']
    nr_coil_turns = int(parameters['nr_coil_turns'])
    nr_turns = int(parameters['nr_turns'])
    monomials = parameters['monomials']
    polynom_a, polynom_b = numpy.zeros((nr_turns, len(monomials))), numpy.zeros((nr_turns, len(monomials)))
    for i in range(nr_turns):
         polynom_a[i,:], polynom_b[i,:] = calc_multipoles_from_radial_coil_flux(flux[i,:], monomials, r1, r2, nr_coil_turns)
    parameters['polynom_a'] = polynom_a
    parameters['polynom_b'] = polynom_b
    parameters['polynom'] = (polynom_a**2 + polynom_b**2)**0.5
    if plot:
        for f in flux:
            plt.plot(f)
        plt.show()

def calc_stats_unique_excitation(parameters):
    try:
        polynom_a = parameters['polynom_a']
        polynom_b = parameters['polynom_b']
        polynom   = parameters['polynom']
        parameters['avg_polynom_b'] = numpy.mean(polynom_b,axis=0)
        parameters['avg_polynom_a'] = numpy.mean(polynom_a,axis=0)
        parameters['std_polynom_b'] = numpy.std(polynom_b,axis=0)
        parameters['std_polynom_a'] = numpy.std(polynom_a,axis=0)
        parameters['avg_polynom']   = numpy.mean(polynom, axis=0)
        parameters['std_polynom']   = numpy.std(polynom, axis=0)
    except KeyError:
        pass
    try:
        norm_polynom_a = parameters['norm_polynom_a']
        norm_polynom_b = parameters['norm_polynom_b']
        norm_polynom   = parameters['norm_polynom']
        parameters['norm_avg_polynom_b'] = numpy.mean(norm_polynom_b,axis=0)
        parameters['norm_avg_polynom_a'] = numpy.mean(norm_polynom_a,axis=0)
        parameters['norm_std_polynom_b'] = numpy.std(norm_polynom_b,axis=0)
        parameters['norm_std_polynom_a'] = numpy.std(norm_polynom_a,axis=0)
        parameters['norm_avg_polynom']   = numpy.mean(norm_polynom, axis=0)
        parameters['norm_std_polynom']   = numpy.std(norm_polynom, axis=0)
    except KeyError:
        pass

def calc_stats_variable_excitation(analysis_list, label='sample_average'):

    norm_polynom_a, norm_polynom_b, norm_polynom = [], [], []
    n = len(analysis_list)
    for i in range(n):
        analysis = analysis_list[i]
        norm_polynom_a.extend(list(analysis['norm_polynom_a']))
        norm_polynom_b.extend(list(analysis['norm_polynom_b']))
        norm_polynom.extend(list(analysis['norm_polynom']))
    parameters = {}
    parameters['monomials'] = analysis['monomials']
    parameters['current1_avg'] = float('nan')
    parameters['file_name'] = 'N.D.'
    parameters['r0'] = analysis['r0']
    parameters['norm_polynom_a'] = numpy.array(norm_polynom_a)
    parameters['norm_polynom_b'] = numpy.array(norm_polynom_b)
    parameters['norm_polynom']   = numpy.array(norm_polynom)
    calc_stats_unique_excitation(parameters)
    return parameters

def merge_unique_excitation_data(analysis_list):
    polynom_a, polynom_b = [], []
    n = len(analysis_list)
    for i in range(n):
        analysis = analysis_list[i]
        polynom_a.extend(list(analysis['polynom_a']))
        polynom_b.extend(list(analysis['polynom_b']))
    parameters = {}
    parameters.update(analysis_list[0])
    parameters['file_name'] = 'MERGED_DATA'
    parameters['polynom_a'] = numpy.array(polynom_a)
    parameters['polynom_b'] = numpy.array(polynom_b)
    calc_normalized_multipoles(main_monomial=parameters['main_monomial'],
                               parameters=parameters,
                               r0=parameters['r0'],
                               is_skew=parameters['is_skew'])
    calc_stats_unique_excitation(parameters)
    return parameters

def calc_normalized_multipoles(main_monomial,parameters,r0, is_skew=False):

    parameters['r0'] = r0
    parameters['main_monomial'] = main_monomial
    parameters['is_skew'] = is_skew
    monomials = parameters['monomials']
    polynom_a = parameters['polynom_a']
    polynom_b = parameters['polynom_b']
    idx = monomials.index(main_monomial)
    if is_skew:
        main_multipole_at_r0 = polynom_a[:,idx] * (r0 ** main_monomial)
    else:
        main_multipole_at_r0 = polynom_b[:,idx] * (r0 ** main_monomial)
    norm_polynom_a = 0 * polynom_a
    norm_polynom_b = 0 * polynom_b
    for i in range(len(monomials)):
        polynom_a_at_r0 = polynom_a[:,i] * (r0 ** monomials[i])
        polynom_b_at_r0 = polynom_b[:,i] * (r0 ** monomials[i])
        norm_polynom_a[:,i] = polynom_a_at_r0 / main_multipole_at_r0
        norm_polynom_b[:,i] = polynom_b_at_r0 / main_multipole_at_r0
    parameters['norm_polynom_a'] = norm_polynom_a
    parameters['norm_polynom_b'] = norm_polynom_b
    parameters['norm_polynom'] = (norm_polynom_a**2 + norm_polynom_b**2)**0.5

def run_radial_coil_analysis(file_name, monomials, main_monomial, r0, conversion, is_skew=False, plot_flag=True, print_flag=True):

    # reads data from file
    parameters  = load_data_from_file(file_name, conversion, plot=plot_flag)

    parameters['file_name'] = file_name

    # calcs multipoles
    parameters['monomials'] = monomials
    calc_multipoles_from_radial_coil_data(parameters, plot=plot_flag)

    # calcs normalized multipoles
    calc_normalized_multipoles(main_monomial=main_monomial,
                               parameters=parameters,
                               r0=r0,
                               is_skew=is_skew)
    # calcs stats
    calc_stats_unique_excitation(parameters)

    # prints data
    if print_flag:
        print_analysis(parameters)



    return parameters

def print_analysis(parameters, normalized=False):
    monomials = parameters['monomials']
    #polynom_a = parameters['polynom_a']
    #polynom_b = parameters['polynom_b']
    if normalized:
        try:
            avg_polynom_a = parameters['norm_avg_polynom_a']
            avg_polynom_b = parameters['norm_avg_polynom_b']
            std_polynom_a = parameters['norm_std_polynom_a']
            std_polynom_b = parameters['norm_std_polynom_b']
            avg_polynom   = parameters['norm_avg_polynom']
            std_polynom   = parameters['norm_std_polynom']
        except KeyError:
            print('normalized data missing!')
            return
    else:
        try:
            avg_polynom_a = parameters['avg_polynom_a']
            avg_polynom_b = parameters['avg_polynom_b']
            std_polynom_a = parameters['std_polynom_a']
            std_polynom_b = parameters['std_polynom_b']
            avg_polynom   = parameters['avg_polynom']
            std_polynom   = parameters['std_polynom']
        except KeyError:
            print('non-normalized data missing!')
            return
    header = '{0:2s} | {1:11s} {2:10s} | {3:11s} {4:10s} | {5:10s} {6:10s}'.format('n','AvgPolyB','StdPolyB','AvgPolyA','StdPolyA', 'AvgAbsPoly','StdAbsPoly')
    print(header)
    print('-'*len(header))
    for i in range(len(monomials)):
        print('{0:02d} | {1:+.4e} {2:.4e} | {3:+.4e} {4:.4e} | {5:.4e} {6:.4e}'.format(1+monomials[i], avg_polynom_b[i], std_polynom_b[i], avg_polynom_a[i], std_polynom_a[i], avg_polynom[i], std_polynom[i]))
    print('-'*len(header))

def plot_excitation_curve(analysis_list, monomial, plot_title = '', poly_type = 'polynom_b', current_flag=True):

    current, multipole = [], []
    for i in range(len(analysis_list)):
        analysis = analysis_list[i]
        this_current = analysis['current1_avg'] if current_flag else i
        current.append(this_current)
        idx = analysis['monomials'].index(monomial)
        avg_polynom = analysis['avg_' + poly_type]
        multipole.append(avg_polynom[idx])

    plt.plot(current, multipole)
    plt.ylabel('integrated multipole (n='+str(monomial)+') [T.m/m^{0}]'.format(monomial))
    plt.grid('on')
    plt.title(plot_title)
    if current_flag:
        plt.xlabel('current [A]')
    else:
        plt.xlabel('index')
    plt.show()

def plot_multipoles(parameters, poly_type='polynom_b', title='', ylim=None, normalized=False, plot_flag=True, save_fname=None):
    plt.clf()
    monomials  = parameters['monomials']
    if normalized:
        avg_polynom = parameters['norm_avg_'+poly_type]
        std_polynom = parameters['norm_std_'+poly_type]
        r0 = parameters['r0']
    else:
        avg_polynom = parameters['avg_'+poly_type]
        std_polynom = parameters['std_'+poly_type]

    colors = ['red' if value < 0 else 'blue' for value in avg_polynom]
    plt.bar(monomials, abs(avg_polynom), yerr=std_polynom, log=True, color=colors, error_kw=dict(elinewidth=2,ecolor='black'), align='center')
    plt.xlabel('harmonics')
    plt.ylabel('multipoles')
    plt.grid('on')
    if ylim is not None:
        plt.ylim(ylim)
    if normalized:
        plt.title('multipoles of ' + poly_type + ' (normalized at r0 = {0} mm)'.format(1000*r0))
    plt.title(title)
    if plot_flag:
        plt.show()
    if save_fname:
        plt.savefig(save_fname)

def run_analysis(files, analysis_parms, variable_excitation=True):
    alist = []
    for file in files:
        analysis = run_radial_coil_analysis(file,
                                            monomials=analysis_parms['monomials'],
                                            main_monomial=analysis_parms['main_monomial'],
                                            is_skew=analysis_parms['is_skew'],
                                            r0=analysis_parms['r0'],
                                            conversion=analysis_parms['conversion'],
                                            plot_flag=False,
                                            print_flag=False)
        alist.append(analysis)
    if variable_excitation:
        avg = [calc_stats_variable_excitation(alist)]
    else:
        parameters = merge_unique_excitation_data(alist)
        calc_stats_unique_excitation(parameters)
        avg = [parameters]
    return avg

def save_multipoles_fig(analysis_list, label1='sample_average', label2='', ylim=None):
    if ylim is None:
        ylim=[1e-5,1e2]
    n = len(analysis_list)
    for i in range(n):
        analysis = analysis_list[i]
        try:
            current = round(analysis['current1_avg'])
        except:
            label2 = ''
        if label2 == '':
            plot_title_a = 'normalized polynom_a for '+label1
            plot_title_b = 'normalized polynom_b for '+label1
            fname_a = label1+'_polynom_a_'+label2+'_{0:02d}.png'.format(i+1) if n > 1 else label1+'_polynom_a.png'
            fname_b = label1+'_polynom_b_'+label2+'_{0:02d}.png'.format(i+1) if n > 1 else label1+'_polynom_b.png'
        else:
            plot_title_a = 'normalized polynom_a for '+label1+', I={0}A'.format(current),
            plot_title_b = 'normalized polynom_b for '+label1+', I={0}A'.format(current),
            fname_a = label1+'_polynom_a_'+label2+'_{0:02d}.png'.format(i+1) if n > 1 else label1+'_polynom_a_'+label2+'.png'
            fname_b = label1+'_polynom_b_'+label2+'_{0:02d}.png'.format(i+1) if n > 1 else label1+'_polynom_b_'+label2+'.png'
        plot_multipoles(analysis,
                        title=plot_title_a,
                        ylim=ylim,
                        normalized=True,
                        plot_flag=False,
                        poly_type='polynom_a',
                        save_fname=fname_a)
        plot_multipoles(analysis,
                        title=plot_title_b,
                        ylim=ylim,
                        normalized=True,
                        plot_flag=False,
                        poly_type='polynom_b',
                        save_fname=fname_b)
