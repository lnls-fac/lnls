#!/usr/bin/python

import sys
import time
import threading
import numpy
import string
import copy
from scipy.optimize import curve_fit
from math import sqrt,exp,log,pi,acos,atan,cos,asin


def g_CIMP(x):

	x=x[0,:]
	g=numpy.zeros(len(x))
	#g=2.691*(1-0.2288964/x)*1/((1+0.16*x)*(1+1.35*numpy.exp(-x/0.2)))
	
	if (max(x)<1):
		g=-18.743261164767357+101.6507221241339*x**(0.33333333)-104.59646433814892*numpy.sqrt(x)+33.73393945878933*x-10.325598001906716*x**(1.5)
	elif (min(x)>1):
		g=1.1976693536243692*(1-0.2660904859953754/x)*1/((1+0.04920104690300144*x)*(1-0.5874697493344921*numpy.exp(-x*0.09913039025775051)))
	else:
		g=2.691*(1-0.2288964/x)*1/((1+0.16*x)*(1+1.35*numpy.exp(-x/0.2)))
			
	return g
	

def gauss_function(x, a, x0, sigma):
    return a*numpy.exp(-(x-x0)**2/(2*sigma**2))

def calc_sigma(pos,profile):

	max_val=max(profile)
	min_val=min(profile)
	profile=(profile-min_val)/max_val
	aux=profile*pos
	ycm=numpy.sum(aux)/numpy.sum(profile)
	aux=profile*(pos-ycm)**2
	yvar=sqrt(numpy.sum(aux)/numpy.sum(profile))	

	return (ycm,yvar)

def Calc_Growth(twiss,param):

	#Define parameters
	brel=sqrt(1-1/param['gamma']**2)	
	ex=param['exi']
	ey=param['eyi']
	ss=param['ssi']
	sp=param['spi']
	
	#Define twiss arrays
	
	s=numpy.zeros(len(twiss))
	betax=numpy.zeros(len(twiss))
	alphax=numpy.zeros(len(twiss))
	betay=numpy.zeros(len(twiss))
	alphay=numpy.zeros(len(twiss))
	Dx=numpy.zeros(len(twiss))
	Dpx=numpy.zeros(len(twiss))
	Dy=numpy.zeros(len(twiss))
	Dpy=numpy.zeros(len(twiss))

	s=twiss[:,0]
	betax=twiss[:,2]
	alphax=twiss[:,3]
	betay=twiss[:,6]
	alphay=twiss[:,7]
	Dx=twiss[:,4]
	Dpx=twiss[:,5]
	Dy=twiss[:,8]
	Dpy=twiss[:,9]

	#Calculate the parameters
	A=param['cluz']*param['Np']*param['r0']**2/(64*numpy.pi**2*brel**3*param['gamma']**4*ex*ey*ss*sp)
	logCIMP=numpy.log(param['gamma']**2*ex*numpy.sqrt(betay*ey)/(param['r0']*betax))
	Hx=1/betax*[Dx**2+(betax*Dpx+alphax*Dx)**2]
	Hy=1/betay*[Dy**2+(betay*Dpy+alphay*Dy)**2]
	SigH=numpy.sqrt(1/sp**2+Hx/ex+Hy/ey)**(-1)
	aCIMP=SigH/param['gamma']*numpy.sqrt(betax/ex)
	bCIMP=SigH/param['gamma']*numpy.sqrt(betay/ey)	
	
	#Calculate Function g
	g_ab=g_CIMP(aCIMP/bCIMP)
	g_ba=g_CIMP(bCIMP/aCIMP)
		
	#Saves values for the ration a/b and b/a
	#f=open('RatioAB.txt','w')
	#for j in range(len(aCIMP[0,:])):
	#	f.write(str(aCIMP[0,j]/bCIMP[0,j])+'\t\t'+str(bCIMP[0,j]/aCIMP[0,j])+'\n')
	#f.close()	
	
	#Calculate Growth Rates
	fp=A*logCIMP*(SigH**2/sp**2)*(g_ba/aCIMP+g_ab/bCIMP)
	fx=A*logCIMP*(-aCIMP*g_ba+Hx*SigH**2/ex*(g_ba/aCIMP+g_ab/bCIMP))
	fy=A*logCIMP*(-bCIMP*g_ab+Hy*SigH**2/ey*(g_ba/aCIMP+g_ab/bCIMP))
	
	#Integrate along the s coordinate
	invTp=2*pi**(3.0/2.0)*numpy.trapz(fp,s)
	invTx=2*pi**(3.0/2.0)*numpy.trapz(fx,s)
	invTy=2*pi**(3.0/2.0)*numpy.trapz(fy,s)
	
	#Calculate growth
	Tp=invTp
	Tx=invTx
	Ty=invTy
	
	return (Tx,Ty,Tp)

# Fucntion that iterates emittances for the case with no harmonic system (simple calculation of the bunch length)	
def Iterate_emittances(twiss,param):

	#Define differences
	i=1
	time=0
	diff1=1
	diff2=1
	diff3=1
	diff4=1
	difftot=diff1+diff2+diff3+diff4
	
	#Calculate U0
	U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+9
	#print U0
	
	#Calculate damping partition numbers
	Jx=1-param['I4']/param['I2']
	Jy=1
	Jp=2+param['I4']/param['I2']
	#print Jx,Jy,Jp
	
	# Caluclate damping times
	taux=(2*param['En']*param['C'])/(Jx*U0*param['cluz'])
	tauy=(2*param['En']*param['C'])/(Jy*U0*param['cluz'])
	taup=(2*param['En']*param['C'])/(Jp*U0*param['cluz'])
	#print taux,tauy,taup
	
	#Define step for iteration
	tt=taux/5
	
	# Synchrotron tune
	Qs0=sqrt(param['ap']*param['hh']*sqrt(param['Vrf']**2-U0**2)/(2*pi*param['En']))

	
	#Cretaes an array that's a subgroup of param
	inter={}
	inter['exi']=param['ex0']
	inter['eyi']=(param['k_dw']+param['k_beta'])*param['ex0']
	inter['ssi']=param['ss0']
	inter['spi']=param['sp0']
	inter['gamma']=param['gamma']
	inter['r0']=param['r0']
	inter['Np']=param['Np']
	inter['cluz']=param['cluz']
		
	while (difftot>10**(-7)):
		(Tx,Ty,Tp)=Calc_Growth(twiss,inter)
		Tx=float(Tx)/param['C']
		Ty=float(Ty)/param['C']
		Tp=float(Tp)/param['C']
		#print Tx,Ty,Tp

		exx=(-param['ex0']+exp(2*tt*(Tx-1/taux))*(param['ex0']+inter['exi']*(-1+Tx*taux)))/(-1+Tx*taux)

		eyy=(-(param['k_dw']*param['ex0']+param['k_beta']*exx*(1-tauy/Ty))+exp(2*tt*(Ty-1/tauy))*((param['k_dw']*param['ex0']+param['k_beta']*exx*(1-tauy/Ty))+inter['eyi']*(-1+Ty*tauy)))/(-1+Ty*tauy)

		spp=(-param['sp0']+exp(tt*(Tp-1/taup))*(param['sp0']+inter['spi']*(-1+Tp*taup)))/(-1+Tp*taup)
		# Accelerating cavity system only

		sss=inter['spi']*param['C']*sqrt(param['ap']*param['En']/(2*pi*param['hh']*(param['Vrf']**2-U0**2)**0.5));
		#print exx,eyy,spp,sss

		diff1=abs(exx-inter['exi'])/inter['exi']
		diff2=abs(eyy-inter['eyi'])/inter['eyi']
		diff3=abs(spp-inter['spi'])/inter['spi']
		diff4=abs(sss-inter['ssi'])/inter['ssi']
		difftot=diff1+diff2+diff3+diff4
		#print difftot
				
		inter['exi']=exx;

		inter['eyi']=eyy;

		inter['spi']=spp;

		inter['ssi']=sss;

		time=i*tt;
		i=i+1

	return (exx,eyy,spp,sss)

# Function that iterates emittances using the results from tracking to calculate bunch length
def Iterate_emittances3HC(twiss,param,phimain,Vmain,phiharm,Vharm):

	#Define differences
	i=1
	time=0
	diff1=1
	diff2=1
	diff3=1
	diff4=1
	difftot=diff1+diff2+diff3+diff4
	
	#Calculate U0
	U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+9
	
	#Calculate synchronous phase
	Phi_sync_nat=asin(U0/param['Vrf'])

	
	#Calculate damping partition numbers
	Jx=1-param['I4']/param['I2']
	Jy=1 
	Jp=2+param['I4']/param['I2']
	#print Jx,Jy,Jp
	
	# Caluclate damping times
	taux=(2*param['En']*param['C'])/(Jx*U0*param['cluz'])
	tauy=(2*param['En']*param['C'])/(Jy*U0*param['cluz'])
	taup=(2*param['En']*param['C'])/(Jp*U0*param['cluz'])
	#print taux,tauy,taup
	
	#Define step for iteration
	tt=taux/5
	
	#RF frequency
	w_rf =2*pi*(param['hh']*param['cluz']/param['C']-param['Detune0'])	#Generator Frequency

	
	#Creates arrays for 3HC calculation
	posz=numpy.zeros(5000)
	perfil=numpy.zeros(5000)
	pot=numpy.zeros(5000)
	
	#Define longitudinal scale array
	posz=numpy.arange(0,5000.)/10-250 # in milimiters
	
	#Cretaes an array that's a subgroup of param
	inter={}
	inter['exi']=param['ex0']
	inter['eyi']=(param['k_dw']+param['k_beta'])*param['ex0']
	inter['spi']=param['sp0']
	inter['gamma']=param['gamma']
	inter['r0']=param['r0']
	inter['Np']=param['Np']
	inter['cluz']=param['cluz']
	
	pot=1/(param['En']*param['C'])*param['cluz']/w_rf*(Vmain*1e3*(cos(Phi_sync_nat-phimain)-numpy.cos(posz/1000*w_rf/param['cluz']+Phi_sync_nat-phimain))+Vharm*1e3/param['mharm']*(cos(param['mharm']*pi-phiharm)-numpy.cos(param['mharm']*posz/1000*w_rf/param['cluz']+param['mharm']*pi-phiharm)))-1/(param['En']*param['C'])*U0*posz/1000
	perfil=numpy.exp(-pot/(param['ap']*param['sp0']**2))
	(pos0,sigma_mm)=calc_sigma(posz,perfil)
	inter['ssi']=sigma_mm/1000			
				
	while (difftot>10**(-7)):
		(Tx,Ty,Tp)=Calc_Growth(twiss,inter)
		Tx=float(Tx)/param['C']
		Ty=float(Ty)/param['C']
		Tp=float(Tp)/param['C']
		#print Tx,Ty,Tp
		
		exx=(-param['ex0']+exp(2*tt*(Tx-1/taux))*(param['ex0']+inter['exi']*(-1+Tx*taux)))/(-1+Tx*taux)

		eyy=(-(param['k_dw']*param['ex0']+param['k_beta']*exx*(1-tauy/Ty))+exp(2*tt*(Ty-1/tauy))*((param['k_dw']*param['ex0']+param['k_beta']*exx*(1-tauy/Ty))+inter['eyi']*(-1+Ty*tauy)))/(-1+Ty*tauy)

		spp=(-param['sp0']+exp(tt*(Tp-1/taup))*(param['sp0']+inter['spi']*(-1+Tp*taup)))/(-1+Tp*taup)


		
		#Calculate bunch length according to the RF potential (Main RF + 3HC)
		pot=1/(param['En']*param['C'])*param['cluz']/w_rf*(Vmain*1e3*(cos(Phi_sync_nat-phimain)-numpy.cos(posz/1000*w_rf/param['cluz']+Phi_sync_nat-phimain))+Vharm*1e3/param['mharm']*(cos(param['mharm']*pi-phiharm)-numpy.cos(param['mharm']*posz/1000*w_rf/param['cluz']+param['mharm']*pi-phiharm)))-1/(param['En']*param['C'])*U0*posz/1000
		perfil=numpy.exp(-pot/(param['ap']*spp**2))
		(pos0,sigma_mm)=calc_sigma(posz,perfil)
		sss=sigma_mm/1000			
		
		#print exx,eyy,spp,sss

		diff1=abs(exx-inter['exi'])/inter['exi']
		diff2=abs(eyy-inter['eyi'])/inter['eyi']
		diff3=abs(spp-inter['spi'])/inter['spi']
		diff4=abs(sss-inter['ssi'])/inter['ssi']
		difftot=diff1+diff2+diff3+diff4
		#print difftot
				
		inter['exi']=exx;

		inter['eyi']=eyy;

		inter['spi']=spp;

		inter['ssi']=sss;

		time=i*tt;
		i=i+1

	return (exx,eyy,spp,sss)

# Function that iterates emittances for the case with no harmonic system (simple calculation of the bunch length) but
# takes into account the longitudinal growth rate due to microwave instability	
def Iterate_emittancesMW(twiss,param,sigS,Curr,GT):

	#Define differences
	i=1
	time=0
	diff1=1
	diff2=1
	diff3=1
	diff4=1
	difftot=diff1+diff2+diff3+diff4
	
	#Calculate U0
	U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+9
	#print U0
	
	#Calculate damping partition numbers
	Jx=1-param['I4']/param['I2']
	Jy=1
	Jp=2+param['I4']/param['I2']
	#print Jx,Jy,Jp
	
	# Caluclate damping times
	taux=(2*param['En']*param['C'])/(Jx*U0*param['cluz'])
	tauy=(2*param['En']*param['C'])/(Jy*U0*param['cluz'])
	taup=(2*param['En']*param['C'])/(Jp*U0*param['cluz'])
	#print taux,tauy,taup
	
	#Define step for iteration
	tt=taux/5
	
	# Synchrotron tune
	Qs0=sqrt(param['ap']*param['hh']*sqrt(param['Vrf']**2-U0**2)/(2*pi*param['En']))

	# Define the interpolation function for Microwave Instability
	microwave=interp2d(sigS,Curr,GT,kind='linear')
	
	#Cretaes an array that's a subgroup of param
	inter={}
	inter['exi']=param['ex0']
	inter['eyi']=(param['k_dw']+param['k_beta'])*param['ex0']
	inter['ssi']=param['ss0']
	inter['spi']=param['sp0']
	inter['gamma']=param['gamma']
	inter['r0']=param['r0']
	inter['Np']=param['Np']
	inter['cluz']=param['cluz']
	sss=param['ss0']
		
	while (difftot>10**(-7)):
		#Add the Microwave growth rate to the longitudinal plane
		DTp=microwave(sss,param['Np'])
		#print DTp

		(Tx,Ty,Tp)=Calc_Growth(twiss,inter)
		Tx=float(Tx)/param['C']
		Ty=float(Ty)/param['C']
		Tp=float(Tp)/param['C']+DTp
		
		exx=(-param['ex0']+exp(2*tt*(Tx-1/taux))*(param['ex0']+inter['exi']*(-1+Tx*taux)))/(-1+Tx*taux)

		#eyy=(-param['ey0']+exp(2*tt*(Ty-1/tauy))*(param['ey0']+inter['eyi']*(-1+Ty*tauy)))/(-1+Ty*tauy)
		eyy=(-(param['k_dw']*param['ex0']+param['k_beta']*exx*(1-tauy/Ty))+exp(2*tt*(Ty-1/tauy))*((param['k_dw']*param['ex0']+param['k_beta']*exx*(1-tauy/Ty))+inter['eyi']*(-1+Ty*tauy)))/(-1+Ty*tauy)

		spp=(-param['sp0']+exp(tt*(Tp-1/taup))*(param['sp0']+inter['spi']*(-1+Tp*taup)))/(-1+Tp*taup)
		# Accelerating cavity system only

		sss=inter['spi']*param['C']*sqrt(param['ap']*param['En']/(2*pi*param['hh']*(param['Vrf']**2-U0**2)**0.5));

		diff1=abs(exx-inter['exi'])/inter['exi']
		diff2=abs(eyy-inter['eyi'])/inter['eyi']
		diff3=abs(spp-inter['spi'])/inter['spi']
		diff4=abs(sss-inter['ssi'])/inter['ssi']
		difftot=diff1+diff2+diff3+diff4
		#print difftot
				
		inter['exi']=exx;

		inter['eyi']=eyy;

		inter['spi']=spp;

		inter['ssi']=sss;

		time=i*tt;
		i=i+1

	return (exx,eyy,spp,sss)


# Function that iterates emittances using the results from tracking to calculate bunch length
# and also takes into account the longitudinal growth rate due to microwave instability	
def Iterate_emittances3HC_MW(twiss,param,phimain,Vmain,phiharm,Vharm,sigS,Curr,GT):

	#Definde differences
	i=1
	time=0
	diff1=1
	diff2=1
	diff3=1
	diff4=1
	difftot=diff1+diff2+diff3+diff4
	
	#Calculate U0
	U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+9
	
	#Calculate synchronous phase
	Phi_sync_nat=asin(U0/param['Vrf'])
	
	#Calculate damping partition numbers
	Jx=1-param['I4']/param['I2']
	Jy=1
	Jp=2+param['I4']/param['I2']
	#print Jx,Jy,Jp
	
	# Caluclate damping times
	taux=(2*param['En']*param['C'])/(Jx*U0*param['cluz'])
	tauy=(2*param['En']*param['C'])/(Jy*U0*param['cluz'])
	taup=(2*param['En']*param['C'])/(Jp*U0*param['cluz'])
	#print taux,tauy,taup
	
	#Define step for iteration
	tt=taux/5
	
	# Synchrotron tune
	Qs0=sqrt(param['ap']*param['hh']*sqrt(param['Vrf']**2-U0**2)/(2*pi*param['En']))

	# Define the interpolation function for Microwave Instability
	microwave=interp2d(sigS,Curr,GT,kind='linear')	

	#RF frequency
	w_rf =2*pi*(param['hh']*param['cluz']/param['C']-param['Detune0'])	#Generator Frequency

	
	#Creates arrays for 3HC calculation
	posz=numpy.zeros(5000)
	perfil=numpy.zeros(5000)
	pot=numpy.zeros(5000)
	
	#Define longitudinal scale array
	posz=numpy.arange(0,5000.)/10-250 # in milimiters
	
	#Cretaes an array that's a subgroup of param
	inter={}
	inter['exi']=param['ex0']
	inter['eyi']=(param['k_dw']+param['k_beta'])*param['ex0']
	inter['spi']=param['sp0']
	inter['gamma']=param['gamma']
	inter['r0']=param['r0']
	inter['Np']=param['Np']
	inter['cluz']=param['cluz']
	
	pot=1/(param['En']*param['C'])*param['cluz']/w_rf*(Vmain*1e3*(cos(Phi_sync_nat-phimain)-numpy.cos(posz/1000*w_rf/param['cluz']+Phi_sync_nat-phimain))+Vharm*1e3/param['mharm']*(cos(param['mharm']*pi-phiharm)-numpy.cos(param['mharm']*posz/1000*w_rf/param['cluz']+param['mharm']*pi-phiharm)))-1/(param['En']*param['C'])*U0*posz/1000
	perfil=numpy.exp(-pot/(param['ap']*param['sp0']**2))
	(pos0,sigma_mm)=calc_sigma(posz,perfil)
	inter['ssi']=sigma_mm/1000			
		
	while (difftot>10**(-7)):
		#Add the Microwave growth rate to the longitudinal plane
		DTp=microwave(sss,param['Np'])
		#print DTp

		(Tx,Ty,Tp)=Calc_Growth(twiss,inter)
		Tx=float(Tx)/param['C']
		Ty=float(Ty)/param['C']
		Tp=float(Tp)/param['C']+DTp
		
		exx=(-param['ex0']+exp(2*tt*(Tx-1/taux))*(param['ex0']+inter['exi']*(-1+Tx*taux)))/(-1+Tx*taux)

		eyy=(-(param['k_dw']*param['ex0']+param['k_beta']*exx*(1-tauy/Ty))+exp(2*tt*(Ty-1/tauy))*((param['k_dw']*param['ex0']+param['k_beta']*exx*(1-tauy/Ty))+inter['eyi']*(-1+Ty*tauy)))/(-1+Ty*tauy)
		spp=(-param['sp0']+exp(tt*(Tp-1/taup))*(param['sp0']+inter['spi']*(-1+Tp*taup)))/(-1+Tp*taup)
		
		#Calculate bunch length according to the RF potential (Main RF + 3HC)
		pot=1/(param['En']*param['C'])*param['cluz']/w_rf*(Vmain*1e3*(cos(Phi_sync_nat-phimain)-numpy.cos(posz/1000*w_rf/param['cluz']+Phi_sync_nat-phimain))+Vharm*1e3/param['mharm']*(cos(param['mharm']*pi-phiharm)-numpy.cos(param['mharm']*posz/1000*w_rf/param['cluz']+param['mharm']*pi-phiharm)))-1/(param['En']*param['C'])*U0*posz/1000
		perfil=numpy.exp(-pot/(param['ap']*spp**2))
		(pos0,sigma_mm)=calc_sigma(posz,perfil)
		sss=sigma_mm/1000			

		
		#print exx,eyy,spp,sss

		diff1=abs(exx-inter['exi'])/inter['exi']
		diff2=abs(eyy-inter['eyi'])/inter['eyi']
		diff3=abs(spp-inter['spi'])/inter['spi']
		diff4=abs(sss-inter['ssi'])/inter['ssi']
		difftot=diff1+diff2+diff3+diff4
		#print difftot
				
		inter['exi']=exx;

		inter['eyi']=eyy;

		inter['spi']=spp;

		inter['ssi']=sss;

		time=i*tt;
		i=i+1

	return (exx,eyy,spp,sss)


