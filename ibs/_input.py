import numpy as _np
import copy as _copy

def read_energy_acceptance_file(fname, eRF):

    # reads raw data from file
    lines = [line.strip() for line in open(fname)]

    # processes raw data
    accp, accn = [], []
    for line in lines:
        if not line or line[0] == '#':
            continue
        values = [float(word) for word in line.split()]
        pos, e_ac = values[4], values[7]
        if e_ac > 0.0:
            accp.append([pos,min(abs(e_ac),eRF)])
        else:
            accn.append([pos,min(abs(e_ac),eRF)])

    accp = _np.array(accp)
    accn = _np.array(accn)
    return (accp,accn)


def read_twiss_file(fname, orig_parameters):

    # reads raw data from file
    lines = [line.strip() for line in open(fname)]

    parameters = _copy.deepcopy(orig_parameters)

    # processes raw data into twiss and element structures
    twiss, elements = [], []
    for line in lines:
        words = line.split()
        if not words or words[0][0] == '*':
            continue
        if words[0][0] == '#':
            if words[0] == '#MCF':
                parameters.mcf = float(words[1])
            elif words[0] == '#I1':
                parameters.latt_i1 = float(words[1])
            elif words[0] == '#I2':
                parameters.latt_i2 = float(words[1])
            elif words[0] == '#I3':
                parameters.latt_i3 = float(words[1])
            elif words[0] == '#I4':
                parameters.latt_i4 = float(words[1])
            elif words[0] == '#I5':
                parameters.latt_i5 = float(words[1])
            elif words[0] == '#I6':
                parameters.latt_i6 = float(words[1])
            else:
                pass
            continue
        else:
            if float(words[3]) > 0:
                values = [float(word) for word in words[2:]]
                values = values + [0,0] # for acceptances insertion latter on
                #print(values)
                twiss.append(values)
                elements.append(words[0])

    twiss = _np.array(twiss)
    elements = _np.array(elements)

    return (elements, twiss, parameters)
