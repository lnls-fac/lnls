import numpy as _np

#Calculates the final bunch legth given by the HHC (no IBS included)def Calc_Sigma(Vh,phih,Vmain,phis,I0,nbunches,param):    print ("-----------------------------------------------")    print ('Calculation of Bunch length across the bunch train')

    f0=param['cluz']/param['C']    #Revolution frequency    f_rf=param['hh']*f0            #Ressoant frequency of the accelerating cavity
    m=param['mharm']             #harmonic number    h=param['hh']                #harmoniy number of the main RF
    Vrf=param['Vrf']            #main RF voltage [V]    sigma_e=param['sp0']        #Natural energy spread
    #Energy loss per turn [eV]
    U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+9

    #Frequency in rad/s:    w_rf=2*pi*f_rf
    #Necessary quantities    #1)Active System (acc. cavities)    #Modified Synchronous phase    Phi_sync=asin(m**2/(m**2-1)*U0/Vrf)    #Natural Synchronous phase    Phi_sync_nat=asin(U0/Vrf)    #2) Passive System (3rd harmonic cavities)    #Harmonic Phase    Phi_synch=-atan2((m*U0/Vrf),sqrt((m**2-1)**2-(m**2*U0/Vrf)**2))    #Harmonic cavity gap voltage    Vhideal=Vrf*sqrt(1/m**2-(U0/Vrf)**2/(m**2-1))

    # print important values
    #print "Ideal case:"    #print "Vg1=",Vrf/1e6, "MV, Vg2=", Vhideal/1e3," kV"    #print "Phi_sync_nat =", Phi_sync_nat    #print "Phi_sync =", Phi_sync    #print "Phi_synch =", Phi_synch    #print "U0 =", U0/1e3,' keV'    sigSh=_np.zeros(nbunches)
    sigCen=_np.zeros(nbunches)
    Ratio=_np.zeros(nbunches)    #Creates arrays for 3HC calculation
    posz=_np.zeros(5000)    perfil=_np.zeros(5000)
    pot=_np.zeros(5000)

    #Define longitudinal scale array
    posz=_np.arange(0,5000.)/10-250. # in milimiters
    #Calculation of the natural bunch length    pot=param['ap']/(param['En']*param['C'])*param['cluz']*Vrf/w_rf*(cos(Phi_sync_nat)-_np.cos(posz/1000*w_rf/param['cluz']+Phi_sync_nat))-param['ap']/(param['En']*param['C'])*U0*posz/1000    perfil=_np.exp(-pot/(param['ap']**2*sigma_e**2))    (pos0,sigma_mm)=calc_sigma(posz,perfil)
    sigma_t=sigma_mm/1000*1/param['cluz']    perfil=perfil**2    R0=_np.trapz(perfil,posz)
    print ('Natural bunch length= {0:0.2f} ps or {1:0.2f} mm' .format(sigma_t*1e12,sigma_mm))    sigS0=sigma_t

    #Calculate bunch length according to the theoretical RF potential (Main RF + 3HC)
    pot=1/(param['En']*param['C'])*param['cluz']/w_rf*(Vrf*(cos(Phi_sync)-_np.cos(posz/1000*w_rf/param['cluz']+Phi_sync))+Vhideal/m*(cos(m*pi-Phi_synch)-_np.cos(m*posz/1000*w_rf/param['cluz']+m*pi-Phi_synch)))-1/(param['En']*param['C'])*U0*posz/1000
    perfil=_np.exp(-pot/(param['ap']*sigma_e**2))    (pos0,sigma_mm)=calc_sigma(posz,perfil)
    sigma_t=sigma_mm/1000*1/param['cluz']    perfil=perfil**2    Rh=_np.trapz(perfil,posz)
    print ('Ideal bunch legnth with 3HC = {0:0.2f} ps or {1:0.2f} mm' .format(sigma_t*1e12, sigma_mm))    print ('Ratio of Touscheck lifetime increase = {0:0.2f}' .format(Rh/R0))    #Calculation of bunch with 3rd harmonic cavity results from tracking    for i in range(nbunches):
        pot=1/(param['En']*param['C'])*param['cluz']/w_rf*(Vmain[i]*1e3*(cos(Phi_sync_nat-phis[i])-_np.cos(posz/1000*w_rf/param['cluz']+Phi_sync_nat-phis[i]))+Vh[i]*1e3/m*(cos(m*pi-phih[i])-_np.cos(m*posz/1000*w_rf/param['cluz']+m*pi-phih[i])))-1/(param['En']*param['C'])*U0*posz/1000
        perfil=_np.exp(-pot/(param['ap']*sigma_e**2))
        (pos0,sigma_mm)=calc_sigma(posz,perfil)        sigSh[i]=sigma_mm/1000*1/param['cluz']        sigCen[i]=pos0/param['cluz']*1/1000

    print ('Total phase shift= {0:0.3f} ps' .format(fabs(phis[0]-phis[nbunches-1])/w_rf*1e12))    print ('Average bunch length = {0:0.2f} +/- {1:0.2f} ps' .format(_np.mean(sigSh)*1e12,_np.std(sigSh)*1e12))    print ('Average bunch length = {0:0.2f} +/- {1:0.2f} mm' .format(_np.mean(sigSh)*1e3*param['cluz'], _np.std(sigSh)*1e3*param['cluz']))    print ('Average bunch length increase = {0:0.2f}' .format(_np.mean(sigSh)/sigS0))

    #Saves calculated bunch lengths and centroid
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
    print('\n')def gauss(x, *p):
    A, mu, sigma = p
    return A*_np.exp(-(x-mu)**2/(2.*sigma**2))

def calc_sigma(pos,profile):

    max_val=max(profile)    min_val=min(profile)    profile=(profile-min_val)/max_val    aux=profile*pos    ycm=_np.sum(aux)/_np.sum(profile)    aux=profile*(pos-ycm)**2
    yvar=sqrt(_np.sum(aux)/_np.sum(profile))

    return (ycm,yvar)