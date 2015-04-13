import numpy as _np
import mathphys as _mp

#Simulations of the phase deviation caused by the transient beam loading  in a 3rd harmonic cavity system#The 3rd harmonic system is passive and so the beam is the responsible for the build up of voltagedef Phase_Shift_Sirius(n_voltas,I_tot,gap,hres,param,n_voltas_max):# n_voltas - number of turns,I_tot - total current, hres (floats) - detune# gap - matrix with logitudinal filling pattern
# param - list of ring parameters    C=param['C']            # Ring Circunference (m)    f0=param['cluz']/C        # Revolution frequency    alfa=param['ap']        # momentum compaction factor    f_rf=param['hh']*f0        # Ressoant frequency of the accelerating cavity
    f_res=f_rf-param['Detune0'] #Generator Frequency    f_resn=hres*f0            # Ressoant frequency of the 3rd harmonic cavity    m=param['mharm']         # harmonic number    Vrf=param['Vrf']*1e-3    # RF voltage in keV    E=param['En']*1e-3        # Energy in keV
    h=int(param['hh'])        # harmonic number of the main RF    # Beam properties
    U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+06
    Js=2+param['I4']/param['I2']
    taus=(2*E*param['C'])/(Js*U0*param['cluz'])
    lambdas = 2/taus


    # Theoretical 3rd harmonic cavity voltage
    Vgh=Vrf*sqrt(1/m**2-(U0/Vrf)**2/(m**2-1)) #(kV)
    #Frequencies in rad/s:    w_rf=2*pi*f_rf
    w_res=2*pi*f_res    w_resn=2*pi*f_resn

    #Vectors definition
    phi=_np.zeros((n_voltas_max+1,h))    delta=_np.zeros((n_voltas_max+1,h))    I0=_np.zeros(h)
    Vtot=_np.zeros(h)
    Vharm=_np.zeros(h)
    Vmain=_np.zeros(h)
    phimain=_np.zeros(h)
    phih=_np.zeros(h)
    phi_final=_np.zeros(h)
    delta_final=_np.zeros(h)    Vb=_np.zeros((h),dtype=complex)
    Vbn=_np.zeros((h),dtype=complex)
    VmainC=_np.zeros((h),dtype=complex)

    #Calculate the current distribution and the energy and phase of the bunches    ini=0
    tot=0
    for i in range(len(gap)):
        I0[ini:(gap[i][0])]=gap[i][1]
        ini=gap[i][0]

    n_bunches=int(_np.sum(I0))    I_tot=I_tot/1000.0 #Ampere - total current    I_b=I_tot/n_bunches #Ampere - current per bunch    I0=I0*I_b    print ("Total stored current : {0:0.3f} A" .format(I_tot))
    print ("Total number of filled bunches : {0:3d}" .format(n_bunches))

    #Necessary quantities    #1)Active System (acc. cavities)    #Loaded Quality Factor    QL=param['Q0']/(1+param['betac'])    #Synchronous phase    Phi_sync=asin(U0/Vrf)    #Detune    Psi=atan2(2*QL*(w_res-w_rf),w_res)
    #Loss Current    I_0=(Vrf*1000)/(param['Rs0']/(1+param['betac'])) #Ampere
    #Loaded Angle    #PhiL=atan2((tan(Psi)-2*I_tot/I_0*cos(pi-Phi_sync)),(1+2*I_tot/I_0*sin(pi-Phi_sync)))
    PhiL=atan2((tan(Psi)-2*I_tot/I_0*cos(Phi_sync)),(1+2*I_tot/I_0*sin(Phi_sync)))    #Gap voltage given by the generator    Vg0 = param['Rs0']/(1+param['betac'])*cos(Psi)/cos(PhiL)*(I_0+2*I_tot*sin(pi-Phi_sync))/1000
    #LossFactor    k0=w_res*param['Rs0']/(2*param['Q0'])    #2) Passive System (3rd harmonic cavities)    #Harmonic Phase    Phi_syncn=-atan2((m*U0/Vrf),sqrt((m**2-1)**2-(m**2*U0/Vrf)**2))/m
    #Detune    Psi_n=atan2(2*param['Q0n']*(w_resn-m*w_rf),w_resn)    #LossFactor    kn=w_resn*param['Rsn']/(2*param['Q0n'])    #Loss Current    I_0n=(Vgh*1000)/param['Rsn'] #Ampere    #Loaded Angle    PhiL_n=atan2((tan(Psi_n)-2*I_tot/I_0n*cos(pi/2-m*Phi_syncn)),(1+2*I_tot/I_0n*sin(pi/2-m*Phi_syncn)))    #Harmonic cavity gap voltage    Vg0n=param['Rsn']*2*I_tot*sin(pi/2-m*Phi_syncn)*cos(Psi_n)/1000
    #Print important values
    print ('-----------------------------------------------')
    print ('Main parameters for calculation')    print ('Ideal Vharm = {0:0.3f} kV' .format(Vrf*sqrt(1/m**2-1/(m**2-1)*(U0/Vrf)**2)))    print ('Induced beam current in main RF = {0:0.3f} kV' .format(2*I_tot*param['Rs0']/(1+param['betac'])*cos(Psi)/1000))
    print ('Optimum Detune = {0:0.3f} kHz' .format(2*I_tot/(Vrf*1000)*param['Rs0']/(1+param['betac'])*cos(pi-Phi_sync)*f_rf/(2*QL)*1e-3))    print ('Actual Detune = {0:0.3f} kHz' .format(tan(Psi)*f_rf/(2*QL)*1e-3))    print ('RF phase shift = {0:0.3f}' .format(-w_res/(2*QL*f0)*tan(Psi)))    print ('Ideal HC shunt impedance = {0:0.3f} MOhms' .format(Vrf**2/(2*I_tot*U0)*(m**2-1)/m**2*1/1000))    print ('Ideal HC phase = {0:0.3f}' .format(Psi_n-pi/2))    print ('Acc. cavitiy phase: Phi_sync = {0:0.3f}' .format(Phi_sync))    print ('HHC cavitiy phase: Phi_syncn = {0:0.3f}' .format(m*Phi_syncn))
    print ('Generator voltages (1) Acc. and (2) HHC:')    print ('        Vg1 = {0:0.3f} kV and Vg2 = {1:0.3f} kV' .format(Vg0,Vg0n))    print ('Cavities loaded angles (1) Acc. and (2) HHC:')    print ('        Psi1 = {0:0.3f} and PhiL1 = {1:0.3f}' .format(Psi,PhiL))    print ('        Psi2 = {0:0.3f} and PhiL2 = {1:0.3f}' .format(Psi_n,PhiL_n))    print ('Relative number of turns: turn/n_turn_max = {0:2.0f}' .format(floor(n_voltas/n_voltas_max)))
    print ('Total number of filled bunches = {0:3.0f}' .format(n_bunches))    print ('-----------------------------------------------\n')

    #Tracking    for l in range((n_voltas/n_voltas_max)):        print ('turn # {0:2d} - Dt_max = {1:0.3e}, V3HC = {2:0.2f} kV, Vmain = {3:0.2f} kV'.format(l+1,(phi[n_voltas_max][0]-phi[n_voltas_max][h-1])/w_rf,_np.average(_np.absolute(Vbn)),_np.average(_np.absolute(Vb))))        if(l>0):
            for p in range(h):
                phi[0][p]=phi[n_voltas_max][p]                delta[0][p]=delta[n_voltas_max][p]
        print ('        phi(first bucket) = {1:0.3e} and phi(last bucket) = {1:0.3e}'.format(phi[0][1],phi[0][n_bunches-1]))
        for n in range(1,n_voltas_max+1): #loop in turns
            for p in range(h):# loop in bunches
                phi[n][p]=phi[n-1][p]+2*pi* alfa*h*delta[n-1][p]
                if (p==0 and n==1 and l==0): # first turn of the first particle
                    Vb[p]=-(k0*I0[p]/f0)/1000                    Vbn[p]=-(kn*I0[p]/f0)/1000
                elif (p==0):# other turns for first particles
                    Dt=(phi[n][0]-phi[n-1][h-1])/w_rf+2*pi/w_rf
                    Vb[0]=(Vb[h-1]-(k0*I0[h-1]/f0)/1000)*complex(cos(w_res*Dt),sin(w_res*Dt))*exp(-w_res*Dt/(2*QL))-(k0*I0[0]/f0)/1000
                    Vbn[0]=(Vbn[h-1]-(kn*I0[h-1]/f0)/1000)*complex(cos(w_resn*Dt),sin(w_resn*Dt))*exp(-w_resn*Dt/(2*param['Q0n']))-(kn*I0[0]/f0)/1000
                else:# other turns for all other particles
                    Dt=(phi[n][p]-phi[n][p-1])/w_rf+2*pi/w_rf
                    Vb[p]=(Vb[p-1]-(k0*I0[p-1]/f0)/1000)*complex(cos(w_res*Dt),sin(w_res*Dt))*exp(-w_res*Dt/(2*QL))-(k0*I0[p]/f0)/1000
                    Vbn[p]=(Vbn[p-1]-(kn*I0[p-1]/f0)/1000)*complex(cos(w_resn*Dt),sin(w_resn*Dt))*exp(-w_resn*Dt/(2*param['Q0n']))-(kn*I0[p]/f0)/1000

                Vtot[p]=Vg0*_np.sin(pi-Phi_sync+phi[n][p]+Psi-PhiL)+_np.real(Vb[p]+Vbn[p])                delta[n][p]=(1-2*lambdas/f0)*delta[n-1][p]+(Vtot[p]-U0)/E
            if (n==n_voltas_max):
                phi_final=phi[n]
                delta_final=delta[n]
    #Harmonic voltage and phase
    Vharm=_np.absolute(Vbn)
    phih=pi/2+_np.angle(Vbn)
    #Main voltage and final synchronous phase
    VmainC=Vg0*complex(sin(pi-Phi_sync+phi_final[p]+Psi-PhiL),cos(pi-Phi_sync+phi_final[p]+Psi-PhiL))+Vb    Vmain=_np.absolute(VmainC)
    phimain=pi/2-_np.angle(VmainC)
    #Saves raw data
    f=open('Landau.txt','w')
    f.write('I0[mA] \t Vharm[kV] \t phi_harm [rad] \t Vmain[kV] \t phi_main [rad] \t phi_final [rad] \t delta_final\n')
    for j in range(h):
        f.write(str('{0:0.5e}'.format(I0[j])) + '\t')
        f.write(str('{0:0.5e}'.format(Vharm[j])) + '\t')
        f.write(str('{0:0.5e}'.format(phih[j])) + '\t')
        f.write(str('{0:0.5e}'.format(Vmain[j])) + '\t')
        f.write(str('{0:0.5e}'.format(phimain[j])) + '\t')
        f.write(str('{0:0.5e}'.format(phi_final[j])) + '\t')
        f.write(str('{0:0.5e}'.format(delta_final[j])) + '\t')
        f.write(str('{0:0.5e}'.format(abs(Vb[j]))) + '\t')
        f.write(str('{0:0.5e}'.format(phase(Vb[j]))) + '\t')
        f.write('\n')
    f.close()
    print('Data saved in file: Landau.txt \n')    return (Vharm,Vmain,phih,phi_final,I0)