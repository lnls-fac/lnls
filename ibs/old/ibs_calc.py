#!/usr/bin/env python3

import numpy as _np

#import other modules for calculation of IBS, lifetime and Landau Cavity tracking and Input File reading
#from Read_input import *
import Read_input as _read_input
import CIMP_3HC as _cimp_3hc
import Lifetime as _lifetime
import sirius_phase1_parameters as sirius_phase1_parameters
import os as _os

def calc_ibs(param, twiss, I0, print_flag=False):

    #Definition of mA
    #mA=1e-3/param['Ccoul']*param['C']/param['cluz']

    #Define new vectors
    Np=_np.zeros(1)
    exi=_np.zeros(1)
    eyi=_np.zeros(1)
    spi=_np.zeros(1)
    ssi=_np.zeros(1)
    LFtous=_np.zeros(1)
    LFine=_np.zeros(1)
    LFelas=_np.zeros(1)

    if print_flag:
        print ('-----------------------------------------------')
        print ('Calculates IBS effects and Lifetime results')
    #Calculates IBS effects on the emittances
    (exi[0],eyi[0],spi[0],ssi[0])=_cimp_3hc.Iterate_emittances(twiss,param)
    #Uses the IBS results to calculates lifetimes
    (LFtous[0],LFine[0],LFelas[0])=_lifetime.Calc_Lifetime(param,I0,twiss,exi[0],eyi[0],spi[0],ssi[0])
    if print_flag:
        print ('Bunch number = {0}'.format(0+1))
        print ('Ib = {0:4.1f} mA' .format(I0[0]))
        print ('ex_fim = {0:0.3f} nm rad' .format(exi[0]*1e9))
        print ('ey_fim = {0:0.3f} pm rad' .format(eyi[0]*1e12))
        print ('sp_fim = {0:0.3f} %' .format(spi[0]*100))
        print ('ss_fim = {0:0.3f} mm or {1:0.3f} ps' .format(ssi[0]*1e3,ssi[0]/param['cluz']*1e12))
        print ('Lifetime results [h]: Touschek = {0:0.2f}, Elastic = {1:0.2f} and Inelastic  = {2:0.2f}'.format(LFtous[0],LFelas[0],LFine[0]))
    if print_flag:
        print ('-----------------------------------------------')
        print ('\n')

    param['Ib'] = I0
    param['ex_fim'] = exi
    param['ey_fim'] = eyi
    param['sp_fim'] = spi
    param['ss_fim'] = ssi
    param['lft_tous'] = LFtous
    param['lft_elas'] = LFelas
    param['lft_ine'] = LFine

    return param

def saves_results(param):

    #Saves final data
    f=open('Lifetime.txt','w')
    f.write('Initial parameters: \n')
    f.write('ex[m rad] \t ey[m rad] \t sp \t\t ss[m]\n')
    f.write(str('{0:0.5e}'.format(param['ex0'])) + '\t')
    f.write(str('{0:0.5e}'.format(param['ey0'])) + '\t')
    f.write(str('{0:0.5e}'.format(param['sp0'])) + '\t')
    f.write(str('{0:0.5e}'.format(param['ss0'])) + '\t')
    f.write('\n\n')
    f.write('Ib[mA] \t LFtous[h] \t LFelas[h] \t LFine [h] \t ex[nm rad] \t ey[pm rad] \t sp[%] \t ss[mm]\n')
    for j in range(Npoints):
        f.write(str('{0:0.4f}'.format(I0[j])) + '\t')
        f.write(str('{0:0.2f}'.format(LFtous[j])) + '\t')
        f.write(str('{0:0.2f}'.format(LFelas[j])) + '\t')
        f.write(str('{0:0.2f}'.format(LFine[j])) + '\t')
        f.write(str('{0:0.5e}'.format(exi[j]*1e9)) + '\t')
        f.write(str('{0:0.5e}'.format(eyi[j]*1e12)) + '\t')
        f.write(str('{0:0.5e}'.format(spi[j]*100)) + '\t')
        f.write(str('{0:0.5e}'.format(ssi[j]*1e3)) + '\n')
    f.close()
    print('Data saved in file: Lifetime.txt')

def calc_ibs_from_rms_folders(folder, sirius_phase_name):

    # selection of input parameters based on sirius phase
    if sirius_phase_name == 'phase1':
        print("using parameters defined in 'sirius_phase1_parameters.py'")
        input_parameters = sirius_phase1_parameters.parameters
        input_parameters['gamma'] = (input_parameters['En']/1e9)/input_parameters['me']
    else:
        raise Exception('sirius phase not defined')

    print('running ibs calculation on random machines in folder ' + folder)
    files = _os.listdir(folder)
    tau_inel, tau_elas = [], []
    tau_tota, tau_tous = [], []
    emitt_x, emitt_y = [], []
    emitt_p, emitt_s = [], []
    for file in files:
        if 'rms' in file:
            print(file+': ', end='')
            parameters  =  input_parameters.copy()
            twiss_fname = _os.path.join(folder, file,'twiss.txt')
            eaccp_fname = _os.path.join(folder, file,'dynap_ma_out.txt')
            parameters, twiss  = _read_input.Read_data(twiss_fname, eaccp_fname, parameters, 0, print_flag=False)
            I0=parameters['I0']
            param = calc_ibs(parameters, twiss, I0, print_flag=False)
            tau_tous.append(param['lft_tous'][0])
            tau_elas.append(param['lft_elas'][0])
            tau_inel.append(param['lft_ine'][0])
            tau_tota.append(1.0/(1.0/param['lft_tous'][0] + 1.0/param['lft_elas'][0] + 1.0/param['lft_ine'][0]))
            emitt_x.append(param['ex_fim']), emitt_y.append(param['ey_fim'])
            emitt_p.append(param['sp_fim']), emitt_s.append(param['ss_fim'])
            print('tau_tous: {0:04.1f}, '.format(tau_tous[-1]), end='')
            print('tau_inel: {0:04.1f}, '.format(tau_inel[-1]), end='')
            print('tau_elas: {0:04.1f}, '.format(tau_elas[-1]), end='')
            print('tau_tota: {0:04.1f}, '.format(tau_tota[-1]), end='')
            print('')
    print('emittance_x  : {0:0.3f} +/- {1:.3f} (from {2:0.3f}) nm.rad'.format(1e9*_np.mean(emitt_x), 1e9*_np.std(emitt_x), 1e9*parameters['ex0']))
    print('emittance_y  : {0:0.3f} +/- {1:.3f} (from {2:0.3f}) pm.rad'.format(1e12*_np.mean(emitt_y), 1e12*_np.std(emitt_y), 1e12*parameters['ey0']))
    print('sigma_e      : {0:0.3f} +/- {1:.3f} (from {2:0.3f}) %'.format(1e2*_np.mean(emitt_p), 1e2*_np.std(emitt_p), 1e2*parameters['sp0']))
    print('sigma_s      : {0:0.3f} +/- {1:.3f} (from {2:0.3f}) mm'.format(1e3*_np.mean(emitt_p), 1e3*_np.std(emitt_p), 1e3*parameters['sp0']))
    print('tau_inelastic: {0:04.1f} +/- {1:3.1f}'.format(_np.mean(tau_inel), _np.std(tau_inel)))
    print('tau_elastic  : {0:04.1f} +/- {1:3.1f}'.format(_np.mean(tau_elas), _np.std(tau_elas)))
    print('tau_touscheck: {0:04.1f} +/- {1:3.1f}'.format(_np.mean(tau_tous), _np.std(tau_tous)))
    print('tau_total    : {0:04.1f} +/- {1:3.1f}'.format(_np.mean(tau_tota), _np.std(tau_tota)))


if __name__ == "__main__":
    #calc_ibs_from_rms_folders('/home/fac_files/data/sirius/si/beam_dynamics/oficial/v07/c05/multi.cod.tune.coup/trackcpp/', sirius_phase_name='phase1')
    calc_ibs_from_rms_folders('/home/fac_files/data/sirius/si/beam_dynamics/oficial/v07/c05/multi.cod.tune.coup/trackcpp/', sirius_phase_name='phase1')
