#!/usr/bin/env python-sirius

import numpy as np
import matplotlib.pyplot as plt

from lnls.rotcoil import RotCoilMeas_SIFCH

RotCoilMeas_SIFCH.lnls_ima_path = '/home/ximenes/repos-dev/'
# RotCoilMeas_SIFCH.excitation_type = 'FC1'
# RotCoilMeas_SIFCH.excitation_type = 'FC2'
RotCoilMeas_SIFCH.excitation_type = 'NOT_DEFINED'

# NOTE: these are fast correctors with slow QS coils installed previously
installation_fc1_previous = {
    '030': '01C2', '066': '02C2',
    '065': '03C2', '078': '04C2',
    '079': '05C2', '070': '06C2',
    '063': '07C2', '010': '08C2',
    '074': '09C2', '068': '10C2',
    '032': '11C2', '052': '12C2',
    '061': '13C2', '064': '14C2',
    '038': '15C2', '021': '16C2',
    '048': '17C2', '054': '18C2',
    '071': '19C2', '009': '20C2',
}

# NOTE: these dicts have information on serial numbers of all 
# recently installed fast correctors
installation_fc1 = {
    '081': '01M1', '040': '11M1',
    '075': '01M2', '051': '11M2',
    '057': '02M1', '029': '12M1',
    '083': '02M2', '004': '12M2',
    '056': '03M1', '015': '13M1',
    '073': '03M2', '027': '13M2',
    '060': '04M1', '023': '14M1',
    '072': '04M2', '008': '14M2',
    '047': '05M1', '005': '15M1',
    '053': '05M2', '058': '15M2',
    '044': '06M1', '059': '16M1',
    '011': '06M2', '076': '16M2',
    '084': '07M1', '086': '17M1',
    '069': '07M2', '024': '17M2',
    '085': '08M1', '046': '18M1',
    '039': '08M2', '020': '18M2',
    '030': '09M1', '012': '19M1',
    '077': '09M2', '006': '19M2',
    '067': '10M1', '049': '20M1',
    '062': '10M2', '014': '20M2',
}

installation_fc2 = {
    '026': '01C3', '052': '11C3',
    '036': '02C3', '055': '12C3',
    '019': '03C3', '084': '13C3',
    '045': '04C3', '043': '14C3',
    '033': '05C3', '050': '15C3',
    '018': '06C3', '037': '16C3',
    '034': '07C3', '028': '17C3',
    '022': '08C3', '042': '18C3',
    '031': '09C3', '035': '19C3',
    '016': '10C3', '013': '20C3',
}


def get_rotcoil_multipoles(type, serial):
    rc = RotCoilMeas_SIFCH(serial)
    data = rc.get_data_set_measurements(type)
    ids = np.zeros(len(data))
    currs = np.zeros(len(data))
    harmonics = data[0].harmonics
    norms = np.zeros((len(data), len(harmonics)))
    skews = np.zeros((len(data), len(harmonics)))
    for i, rcd in enumerate(data):
        ids[i] = rcd.id
        currs[i] = rcd.main_coil_current_avg
        norms[i, :] = rcd.intmpole_normal_avg
        skews[i, :] = rcd.intmpole_skew_avg
    # sort
    inds = np.argsort(ids)
    ids = ids[inds]
    currs = currs[inds]
    norms = norms[inds, :]
    skews = skews[inds, :]
    return currs, harmonics, norms, skews, rc


def calc_all_data(type, serials):
    all_currs = []
    all_norms = []
    all_skews = []
    for i, serial in enumerate(serials):
        currs, harmonics, norms, skews, _ = get_rotcoil_multipoles(type, serial)
        all_currs += list(currs)
        for multp in norms:
            all_norms.append(multp)
        for multp in skews:
            all_skews.append(multp)
    all_norms = np.array(all_norms)
    all_skews = np.array(all_skews)
    return all_currs, all_norms, all_skews, harmonics


def calc_fitting(all_currs, all_norms, all_skews, currs_fit=None, order=1):
    if currs_fit is None:
        currs_fit = np.linspace(0,1,10)
    currs_fit_tile = np.tile(currs_fit, (15,1))
    p, *_ = np.polyfit(all_currs, all_norms, 4, full=True)
    norms_fit = np.polyval(p, currs_fit_tile.T)
    p, *_ = np.polyfit(all_currs, all_skews, 4, full=True)
    skews_fit = np.polyval(p, currs_fit_tile.T)
    return currs_fit, norms_fit, skews_fit


def calc_avg_multipoles(type, serials, currs_fit=None, order=1):
    all_currs, all_norms, all_skews, harmonics = calc_all_data(type, serials)
    calc_fitting(all_currs, all_norms, all_skews, currs_fit, order=order)
    currs_fit, norms_fit, skews_fit = calc_fitting(all_currs, all_norms, all_skews, currs_fit, order=order)
    return currs_fit, norms_fit, skews_fit, all_currs, all_norms, all_skews, harmonics


def plot_excitation_curves(fc_dict, exc_type, plot_type):

    if fc_dict == installation_fc1:
        mag_type = 'FC1'
    else:
        mag_type = 'FC2'
    RotCoilMeas_SIFCH.excitation_type = mag_type
    serials = list(fc_dict.keys())
    for serial in serials:
        currs, harmonics, norms, skews, _ = get_rotcoil_multipoles(exc_type, serial)
        if plot_type == 'ch':
            plt.plot(currs, norms[:, 0]*1e5, label=serial)
        else:
            plt.plot(currs, skews[:, 0]*1e5, label=serial)
        # currs_fit = np.linspace(0,1.5,10)
        # currs_fit, norms_fit, skews_fit, all_currs, all_norms, all_skews = calc_avg_multipoles('ch', [serial, ], currs_fit)
        # plt.plot(all_currs, all_norms[:, 0], 'o')
        # plt.plot(currs_fit, norms_fit[:, 0])
    plt.xlabel('Current [A]')
    plt.ylabel('Integrated ' + plot_type + ' dipolar field / Brho@3GeV [urad]')
    serials = [el[1:] for el in serials]
    serials = ' '.join(sorted(serials))
    plt.title(exc_type.upper() + ' excitation curves for ' + mag_type + ' magnets\n' + 'serials: ' + serials)
    plt.grid()
    plt.tight_layout()
    plt.show()


def plot_average_excitation_curve(fc_dict, exc_type, plot_type):
    serials = list(fc_dict.keys())
    if fc_dict == installation_fc1:
        mag_type = 'FC1'
    else:
        mag_type = 'FC2'
    RotCoilMeas_SIFCH.excitation_type = mag_type
    currs_fit = np.linspace(0,1.5,10)
    currs_fit, norms_fit, skews_fit, all_currs, all_norms, all_skews, harmonics = calc_avg_multipoles(exc_type, serials, currs_fit)
    conv = 1e6/10  # T.m -> urad # 3 GeV
    if plot_type == 'ch':
        ylabel = 'Integrated By field / Brho@3GeV [urad]'
        plt.plot(all_currs, conv * all_norms[:, 0], 'o')
        plt.plot(currs_fit, conv * norms_fit[:, 0])
    else:
        ylabel = 'Integrated Bx field / Brho@3GeV [urad]'
        plt.plot(all_currs, conv * all_skews[:, 0], 'o')
        plt.plot(currs_fit, conv * skews_fit[:, 0])
    
    plt.xlabel('Current [A]')
    plt.ylabel(ylabel)
    serials = [el[1:] for el in serials]
    serials = ' '.join(sorted(serials))
    plt.title(exc_type.upper() + ' excitation curves for ' + mag_type + ' magnets\n' + 'serials: ' + serials)
    plt.grid()
    plt.tight_layout()
    plt.show()


def create_excdata_file(fc_dict, type, sign):
    serials = list(fc_dict.keys())
    if fc_dict == installation_fc1:
        mag_type = 'FC1'
    else:
        mag_type = 'FC2'
    RotCoilMeas_SIFCH.excitation_type = mag_type
    currs_fit = np.linspace(0, 1.0, 5)
    currs_fit, norms_fit, skews_fit, all_currs, all_norms, all_skews, harmonics = calc_avg_multipoles(type, serials, currs_fit, order=2)
    norms_fit *= sign
    skews_fit *= sign
    norms_fit -= norms_fit[0, :]
    skews_fit -= skews_fit[0, :]
    harms = [1, 3, 9]

    for i, curr in enumerate(currs_fit):
        stg = f'{curr:+08.4f}  '
        for idx, j in enumerate(harmonics):
            if j in harms:
                stg += f'{norms_fit[i, idx]:+.4e} '
                stg += f'{skews_fit[i, idx]:+.4e} '
        print(stg)

    currs_fit = -1 * np.flip(currs_fit)
    norms_fit = -1 * np.flip(norms_fit, 0)
    skews_fit = -1 * np.flip(skews_fit, 0)

    for i, curr in enumerate(currs_fit):
        stg = f'{curr:+08.4f}  '
        for idx, j in enumerate(harmonics):
            if j in harms:
                stg += f'{norms_fit[i, idx]:+.4e} '
                stg += f'{skews_fit[i, idx]:+.4e} '
        print(stg)
    

def run():

    installation_fc1.update(installation_fc1_previous)
    # plot_excitation_curves(fc_dict=installation_fc1, exc_type='ch', plot_type='ch')
    # plot_excitation_curves(fc_dict=installation_fc1, exc_type='ch', plot_type='cv')
    # plot_excitation_curves(fc_dict=installation_fc1, exc_type='cv', plot_type='cv')
    # plot_excitation_curves(fc_dict=installation_fc1, exc_type='cv', plot_type='ch')
    # plot_excitation_curves(fc_dict=installation_fc2, exc_type='ch', plot_type='ch')
    # plot_excitation_curves(fc_dict=installation_fc2, exc_type='ch', plot_type='cv')
    # plot_excitation_curves(fc_dict=installation_fc2, exc_type='cv', plot_type='cv')
    # plot_excitation_curves(fc_dict=installation_fc2, exc_type='cv', plot_type='ch')
    # plot_average_excitation_curve(fc_dict=installation_fc1, exc_type='ch', plot_type='ch')
    # plot_average_excitation_curve(fc_dict=installation_fc1, exc_type='ch', plot_type='cv')
    # plot_average_excitation_curve(fc_dict=installation_fc1, type='cv')
    # plot_average_excitation_curve(fc_dict=installation_fc2, type='ch')
    # plot_average_excitation_curve(fc_dict=installation_fc2, type='cv')

    # create_excdata_file(installation_fc1, 'ch', sign=+1)
    # create_excdata_file(installation_fc1, 'cv', sign=-1)
    # create_excdata_file(installation_fc2, 'ch', sign=+1)
    create_excdata_file(installation_fc2, 'cv', sign=-1)


if __name__ == '__main__':
    run()
