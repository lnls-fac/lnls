import numpy as _np

#Calculates the final bunch legth given by the HHC (no IBS included)

    f0=param['cluz']/param['C']    #Revolution frequency
    m=param['mharm']             #harmonic number
    Vrf=param['Vrf']            #main RF voltage [V]
    #Energy loss per turn [eV]
    U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+9

    #Frequency in rad/s:


    # print important values
    #print "Ideal case:"
    sigCen=_np.zeros(nbunches)
    Ratio=_np.zeros(nbunches)
    posz=_np.zeros(5000)
    pot=_np.zeros(5000)

    #Define longitudinal scale array
    posz=_np.arange(0,5000.)/10-250. # in milimiters

    sigma_t=sigma_mm/1000*1/param['cluz']
    print ('Natural bunch length= {0:0.2f} ps or {1:0.2f} mm' .format(sigma_t*1e12,sigma_mm))

    #Calculate bunch length according to the theoretical RF potential (Main RF + 3HC)
    pot=1/(param['En']*param['C'])*param['cluz']/w_rf*(Vrf*(cos(Phi_sync)-_np.cos(posz/1000*w_rf/param['cluz']+Phi_sync))+Vhideal/m*(cos(m*pi-Phi_synch)-_np.cos(m*posz/1000*w_rf/param['cluz']+m*pi-Phi_synch)))-1/(param['En']*param['C'])*U0*posz/1000
    perfil=_np.exp(-pot/(param['ap']*sigma_e**2))
    sigma_t=sigma_mm/1000*1/param['cluz']
    print ('Ideal bunch legnth with 3HC = {0:0.2f} ps or {1:0.2f} mm' .format(sigma_t*1e12, sigma_mm))
        pot=1/(param['En']*param['C'])*param['cluz']/w_rf*(Vmain[i]*1e3*(cos(Phi_sync_nat-phis[i])-_np.cos(posz/1000*w_rf/param['cluz']+Phi_sync_nat-phis[i]))+Vh[i]*1e3/m*(cos(m*pi-phih[i])-_np.cos(m*posz/1000*w_rf/param['cluz']+m*pi-phih[i])))-1/(param['En']*param['C'])*U0*posz/1000
        perfil=_np.exp(-pot/(param['ap']*sigma_e**2))
        (pos0,sigma_mm)=calc_sigma(posz,perfil)




    f=open('BunchLength.txt','w')
    f.write('bunch number \t sigmaz [ps] \t centroid [ps]\n')
    for j in range(nbunches):
        f.write(str('{0:4d}'.format(j)) + '\t')
        f.write(str('{0:0.5e}'.format(sigSh[j])) + '\t')
        f.write(str('{0:0.5e}'.format(sigCen[j])) + '\t')
        f.write('\n')
    f.close()
    print('Data saved in file: BunchLength.txt')
    print("-----------------------------------------------")
    print('\n')
    A, mu, sigma = p
    return A*_np.exp(-(x-mu)**2/(2.*sigma**2))

def calc_sigma(pos,profile):

    max_val=max(profile)
    yvar=sqrt(_np.sum(aux)/_np.sum(profile))

    return (ycm,yvar)