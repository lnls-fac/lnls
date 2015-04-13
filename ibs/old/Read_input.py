import numpy as _np
import math as _math

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

def read_twiss_file(fname, param):

    # reads raw data from file
    lines = [line.strip() for line in open(fname)]

    # processes raw data into twiss and element structures
    twiss, elements = [], []
    for line in lines:
        words = line.split()
        if not words or words[0][0] == '*':
            continue
        if words[0][0] == '#':
            if words[0] == '#I1':
                param['I1'] += float(words[1])
            elif words[0] == '#I2':
                param['I2'] += float(words[1])
            elif words[0] == '#I3':
                param['I3'] += float(words[1])
            elif words[0] == '#I4':
                param['I4'] += float(words[1])
            elif words[0] == '#I5':
                param['I5'] += float(words[1])
            elif words[0] == '#I6':
                param['I6'] += float(words[1])
            else:
                pass
            continue
        if words[0][0] == '@':
            if words[1] == 'K_beta':
                param['k_beta'] = float(words[3])
            elif words[1] == 'K_dw':
                param['k_dw'] = float(words[3])-param['k_beta']
            elif words[1] == 'EX':
                param['ex0'] = float(words[3])
        else:
            if float(words[3]) > 0:
                values = [float(word) for word in words[2:]]
                values = values + [0,0] # for acceptances insertion latter on
                #print(values)
                twiss.append(values)
                elements.append(words[0])

    twiss = _np.array(twiss)
    elements = _np.array(elements)

    #print('Is: ', param['I1'], param['I2'], param['I3'], param['I4'])

    return (elements, twiss)

def Read_data(ftwiss,facc,param,HC,print_flag=False):


    # read twiss file
    elements, twiss = read_twiss_file(ftwiss, param)

    #Calcultaes emittances
    param = calc_emitt(param,HC,print_flag)

    #Calculates RF acceptance
    U0  = param['Cgamma']/(2*_math.pi)*(param['En']/1e+9)**4*param['I2']*1e+9
    q   = param['Vrf']/U0
    eRF = _math.sqrt(2*U0/(_math.pi*param['ap']*param['hh']*param['En'])*(_math.sqrt(q**2-1)-_math.acos(1/q)))

    # reads file with energy acceptance
    taccp, taccn = read_energy_acceptance_file(facc, eRF)


    # propagates acceptance data to all machine periods
    accp, accn = _np.array(taccp), _np.array(taccn)
    for i in range(10-1):
        taccp2 = _np.array(taccp)
        taccp2[:,0] += (i+1) * param['C']/10.0
        accp = _np.vstack((accp, taccp2))
        taccn2 = _np.array(taccn)
        taccn2[:,0] += (i+1) * param['C']/10.0
        accn = _np.vstack((accn, taccn2))

    #interpolate acceptance in along the twiss positions
    accptwiss=_np.zeros(len(twiss))
    accptwiss=_np.interp(twiss[:,0],accp[:,0],accp[:,1])
    accntwiss=_np.zeros(len(twiss))
    accntwiss=_np.interp(twiss[:,0],accn[:,0],accn[:,1])

    #concatenate the acceptance columns in the twiss matrix
    twiss[:,12]=accptwiss
    twiss[:,13]=accntwiss

    return param, twiss

def MGrowthRate(fname):

    #Read file with microwave gorth rates and currents
    tfile = open(fname,'rU')
    sigS=[]
    Curr=[]
    GT=[]
    for aline in tfile:
        values=aline.split()
        sigS.append(float(values[0])*1e-3)
        Curr.append(float(values[1]))
        GT.append(float(values[2]))

    tfile.close()
    return (sigS,Curr,GT)


def calc_emitt(param,HC,print_flag=False):


    mA=1e-3/param['Ccoul']*param['C']/param['cluz']
    param['Np'] = mA*param['I0']

    #Calculate energy spread
    sigp=_math.sqrt(param['Cq']*param['gamma']**2*(param['I3']/(2*param['I2']+param['I4'])))
    param['sp0']=sigp

    #Calculate bunch length
    U0=param['Cgamma']/(2*_math.pi)*(param['En']/1e+9)**4*param['I2']*1e+9
    if (HC==0):
        param['ss0']=param['sp0']*param['C']*_math.sqrt(param['ap']*param['En']/(2*_math.pi*param['hh']*(param['Vrf']**2-U0**2)**0.5))
    elif (HC==1):
        Qs0=_math.sqrt(param['ap']*param['hh']*_math.sqrt(param['Vrf']**2-U0**2)/(2*_math.pi*param['En']))
        param['ss0']=0.432343807*(param['ap']*param['hh']*_math.pi*sigp/Qs0)**(0.5)*param['C']/(2*pi*param['hh'])


    #Calculate ex
    Jx=1-param['I4']/param['I2']
    ex=param['Cq']*param['gamma']**2*param['I5']/(Jx*param['I2'])
    param['ex0']=ex

    #Calculate ey
    ey=(param['k_beta']+param['k_dw'])*ex
    param['ey0']=ey

    mA=1e-3/param['Ccoul']*param['C']/param['cluz']

    if print_flag:
        for j in range(len(I0)):
            print ('-----------------------------------------------')
            print ('Initial beam parameters')
            print ('Itot   = {0:4.1f} mA' .format(param['Np'][j]*param['hh']/mA))
            print ('ex_ini = {0:0.3f} nm rad' .format(param['ex0']*1e9))
            print ('ey_ini = {0:0.3f} pm rad' .format(param['ey0']*1e12))
            print ('sp_ini = {0:0.3f} %' .format(param['sp0']*100))
            print ('ss_ini = {0:0.3f} mm' .format(param['ss0']*1e3))
            print ('-----------------------------------------------')
            print ('\n')

    return param
