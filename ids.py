#!/usr/bin/env python3

from tkinter import Tk as _Tk
from tkinter.filedialog import askopenfilename as _askopenfilename
import lnls as _app
import os as _os
import numpy as _np
import matplotlib.pyplot as _plt
import mathphys as _mathphys

def read_kicktable(fname):

    # reads raw data from file
    lines = [line.strip() for line in open(fname)]

    # reads header
    data, idx = [], []
    for i in range(len(lines)):
        if lines[i] == 'START':
            idx.append(i+1)
            continue
        try:
            data.append(float(lines[i]))
        except ValueError:
            pass
    id_length, id_nrpts_x, id_nrpts_y = data[0], int(data[1]), int(data[2])

    # reads kickx
    id_posx, id_posy = [float(word) for word in lines[idx[0]].split()], []
    id_kickx = _np.zeros((id_nrpts_y,id_nrpts_x))
    data = lines[idx[0]+1:idx[1]-1]
    idx_y = 0
    for line in data:
        try:
            datum = [float(word) for word in line.split()]
            id_posy.append(datum[0])
            id_kickx[idx_y,:] = datum[1:]
            idx_y += 1
        except ValueError:
            pass

    # reads kicky
    id_kicky = _np.zeros((id_nrpts_y,id_nrpts_x))
    data = lines[idx[1]+1:]
    idx_y = 0
    for line in data:
        try:
            datum = [float(word) for word in line.split()]
            id_kicky[idx_y,:] = datum[1:]
            idx_y += 1
        except ValueError:
            pass

    return (id_length, id_posx, id_posy, id_kickx, id_kicky)

def select_kicktable_file():

    default_folder = _os.path.join(_app.folder_mml_sirius_ids,'id_modelling')
    opt = {'initialdir':default_folder, 'title':'select kicktable file', 'defaultextension':'.txt', 'filetypes':[('text files', '*.txt')]}
    _Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    fname = _askopenfilename(**opt) # show an "Open" dialog box and return the path to the selected file
    return fname

def plot_kicktable(fname=None, energy = 3.0, print_flag=True, savefigs_flag=True, display_flag=True):
    """Plot ID kicktable stored in a file

    Accepts name of file witk ID kicktable. If not provided a dialogbox will
    show for users input.

    Keyword arguments:
    fname          -- Name of the file with ID kicktable
    energy         -- Beam energy [GeV]
    print_flag     -- True/False. If true (default), prints kicktable info in stdout.
    savefigs_flag  -- True/False.
    display_flag   -- True/False.

    Returns:
    id_length --
    id_posx   --
    id_posy   --
    id_kickx  --
    id_kicky  --

    """
    #fname = _os.path.join(_app.folder_mml_sirius_ids,'id_modelling','U25','U25_kicktable_4meters.txt')
    # if fname is missing opens dialogbox for users selection of filename
    if fname is None:
        fname = select_kicktable_file()
        if not fname:
            return None

    # reads kicktable from file
    id_length, id_posx, id_posy, id_kickx, id_kicky = read_kicktable(fname)

    if print_flag:
        print('filename : {0}'.format(fname))
        print('length[m]: {0}'.format(id_length))
        print('nrpts_x  : {0}'.format(len(id_posx)))
        print('nrpts_y  : {0}'.format(len(id_posy)))
        print('posx[mm] : {0} ... {1}'.format(1000*id_posx[0], 1000*id_posx[-1]))
        print('posy[mm] : {0} ... {1}'.format(1000*id_posy[0], 1000*id_posy[-1]))

    brho,_,_,_,_ = _mathphys.beam_optics.beam_rigidity(energy = energy * 1e9)
    _os.path.basename(fname)

    # kickx
    plot_idx = [int(len(id_posy)/2), int(len(id_posy)/4), 0]
    leg = []
    print(plot_idx)
    for i in plot_idx:
        _plt.plot(1000*_np.array(id_posx), (1e6/brho**2)*id_kickx[i,:])
        leg.append('{0:+.2f} mm'.format(1000*id_posy[i]))
    _plt.xlabel('posx [mm]'), _plt.ylabel('kickx [um]')
    _plt.grid(), _plt.suptitle('Insertion Device Horizontal Kick')
    _plt.legend(leg)
    if savefigs_flag:
        _plt.savefig('kickx.svg')
    if display_flag:
        _plt.show()
    _plt.clf()

    # kicky
    plot_idx = [int(len(id_posx)/2), int(len(id_posx)/4), 0]
    leg = []
    for i in plot_idx:
        _plt.plot(1000*_np.array(id_posy), (1e6/brho**2)*id_kicky[:,i])
        leg.append('{0:+.2f} mm'.format(1000*id_posx[i]))
    _plt.xlabel('posy [mm]'), _plt.ylabel('kicky [um]')
    _plt.grid(), _plt.suptitle('Insertion Device Vertical Kick')
    _plt.legend(leg)
    if savefigs_flag:
        _plt.savefig('kicky.svg')
    if display_flag:
        _plt.show()
    _plt.clf()

    return (id_length, id_posx, id_posy, id_kickx, id_kicky)
