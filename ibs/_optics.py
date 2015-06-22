import mathphys as _mp
import math as _math
import mathphys.beam_optics as _bo

_twiss_dict = {'s':0,'len':1,'mux':2,'betax':3,'alphax':4,'etax':5,'etapx':6,'muy':7,'betay':8,'alphay':9,'etay':10,'etapy':11,'eaccp_pos':12,'eaccp_neg':13}

def get_twiss(twiss, label):
    return twiss[:,_twiss_dict[label]]

def calc_rf_acceptance(parameters):

    parameters.U0 = _bo.calc_U0(parameters.beam_energy, parameters.latt_i2 + 1*parameters.ids_i2)
    parameters.q  = _bo.calc_overvoltage(parameters.Vrf, parameters.U0)
    parameters.rf_acceptance = _bo.calc_rf_acceptance(parameters.U0, parameters.mcf, parameters.harmonic_number, parameters.q, parameters.beam_energy)

def calc_emittances(parameters):

    p = parameters
    p.nr_electrons_per_bunch = _bo.calc_number_of_electrons(energy=p.beam_energy, circumference=p.circumference, current=p.beam_current_per_bunch)
    p.natural_sigmae = _bo.calc_natural_energy_spread(p.beam_energy, p.latt_i2, p.latt_i3, p.latt_i4)
    #p.natural_sigmal = _bo.calc_natural_bunch_length(p.beam_energy, p.circumference, p.natural_sigmae, p.U0, p.mcf, p.harmonic_number, p.Vrf, p.hcavities)

    print(p)

# #Calculate bunch length
# U0=param['Cgamma']/(2*_math.pi)*(param['En']/1e+9)**4*param['I2']*1e+9
# if (HC==0):
#     param['ss0']=param['sp0']*param['C']*_math.sqrt(param['ap']*param['En']/(2*_math.pi*param['hh']*(param['Vrf']**2-U0**2)**0.5))
# elif (HC==1):
#     Qs0=_math.sqrt(param['ap']*param['hh']*_math.sqrt(param['Vrf']**2-U0**2)/(2*_math.pi*param['En']))
#     param['ss0']=0.432343807*(param['ap']*param['hh']*_math.pi*sigp/Qs0)**(0.5)*param['C']/(2*pi*param['hh'])


    return


    parameters.Np = parameters.beam_current_per_bunch * mA

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
