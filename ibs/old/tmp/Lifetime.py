#!/usr/bin/python

import sys
import time
import threading
import numpy
import string
from math import sqrt,exp,log,pi,atan,sin
from scipy.integrate import quad


def Calc_Lifetime(Pmed,param,I0,twiss,ex,ey,sigP,sigS):
# float Pmed = residual gas pressure (nTorr)
# Dicionary param = basic machine parameters
# array I0,twiss,acc = current distribution(A), twiss parameters and momentum acceptance along the ring (%)
# array ex,ey,sigP,sigS = Calculated values for emittances, energy spread and bunch length for each bunch

	#Twiss parameters
	s=numpy.zeros(len(twiss))
	betax=numpy.zeros(len(twiss))
	alphax=numpy.zeros(len(twiss))
	betay=numpy.zeros(len(twiss))
	alphay=numpy.zeros(len(twiss))
	Dx=numpy.zeros(len(twiss))
	Dpx=numpy.zeros(len(twiss))
	Dy=numpy.zeros(len(twiss))
	Dpy=numpy.zeros(len(twiss))
	accp=numpy.zeros(len(twiss))
	accn=numpy.zeros(len(twiss))

	s=twiss[:,0]
	betax=twiss[:,2]
	alphax=twiss[:,3]
	betay=twiss[:,6]
	alphay=twiss[:,7]
	Dx=twiss[:,4]
	Dpx=twiss[:,5]
	Dy=twiss[:,8]
	Dpy=twiss[:,9]
	accp=twiss[:,10]
	accn=twiss[:,11]
	Ds=numpy.zeros(len(twiss))
	acc=numpy.zeros(len(twiss))
	
	Ds=s-numpy.roll(s,1)
	Ds[0]=0
	acc=numpy.minimum(accp,accn)

	#Calculate average beta functions
	betax_avg=numpy.average(betax,weights=Ds)
	betay_avg=numpy.average(betay,weights=Ds)
	#print "<betax> = ",betax_avg
	#print "<betay> = ",betay_avg
	#print "<acc> = ",numpy.average(acc,weights=Ds)#, "<accp> = ",numpy.average(accp,weights=Ds), "<accn> = ",numpy.average(accn,weights=Ds) 

	#Machine parameters
	C = param['C'] #Circunference (m)
	frev=param['C']/param['cluz'] #Rev. freq (Hz)		
	theta_x = sqrt(param['Ax']/betax_avg)
	theta_y = sqrt(param['Ay']/betay_avg)
	R=theta_y/theta_x
	FR=pi+(R**2+1.0)*sin(2*atan(R))+2.0*(R**2.0-1.0)*atan(R)
	
	#Elastic Scattering Lifetime
	Telas=10.25*2.0*pi/FR*(param['En']/1.0e+09)**2*param['Ay']/(betay_avg*Pmed)
	
	#Inelastic Scattering Lifetime
	Tine=1/(0.0065*Pmed*log(1./(numpy.average(acc,weights=Ds))))
	
	#Touschek Lifetime Calculation
	sigx=numpy.zeros(len(twiss))
	sigy=numpy.zeros(len(twiss))
	epsilon=numpy.zeros(len(twiss))
	Tv=numpy.zeros(len(twiss))
	dsdT=numpy.zeros(len(twiss)-1)
		
	sigx=numpy.sqrt(betax*ex+(Dx*sigP)**2)
	sigy=numpy.sqrt(betay*ey+(Dy*sigP)**2)
	epsilonp=(accp)**2*(betax/ex)/(1957.0*param['En']/1.0e+09)**2
	epsilonn=(accn)**2*(betax/ex)/(1957.0*param['En']/1.0e+09)**2
	Tvp=(5.39e17*(param['En']/1.0e+09)**2*(accp)**3*sigx*sigy*sigS/(De(epsilonp)*C))**(-1.0)
	Tvn=(5.39e17*(param['En']/1.0e+09)**2*(accn)**3*sigx*sigy*sigS/(De(epsilonn)*C))**(-1.0)
 	aux=0.5*(numpy.average(Tvp,weights=Ds)+numpy.average(Tvn,weights=Ds))
 	Ttous=1.0/(aux*I0)
	 
	return (Ttous,Tine,Telas)
	
def Calc_Lifetime_Matlab(Pmed,param,twiss,ex,ey,sigP,sigS):
# float Pmed = residual gas pressure (nTorr)
# Dicionary param = basic machine parameters
# array I0,twiss,acc = current distribution(A), twiss parameters and momentum acceptance along the ring (%)
# array ex,ey,sigP,sigS = Calculated values for emittances, energy spread and bunch length for each bunch

	#Twiss parameters
	s=numpy.zeros(len(twiss))
	betax=numpy.zeros(len(twiss))
	alphax=numpy.zeros(len(twiss))
	betay=numpy.zeros(len(twiss))
	alphay=numpy.zeros(len(twiss))
	Dx=numpy.zeros(len(twiss))
	Dpx=numpy.zeros(len(twiss))
	Dy=numpy.zeros(len(twiss))
	Dpy=numpy.zeros(len(twiss))
	accp=numpy.zeros(len(twiss))
	accn=numpy.zeros(len(twiss))

	s=twiss[:,0]
	betax=twiss[:,2]
	alphax=twiss[:,3]
	betay=twiss[:,6]
	alphay=twiss[:,7]
	Dx=twiss[:,4]
	Dpx=twiss[:,5]
	Dy=twiss[:,8]
	Dpy=twiss[:,9]
	accp=twiss[:,10]
	accn=twiss[:,11]
	Ds=numpy.zeros(len(twiss))
	acc=numpy.zeros(len(twiss))
	
	Ds=s-numpy.roll(s,1)
	Ds[0]=0
	acc=numpy.minimum(accp,accn)

	#Calculate average beta functions
	betax_avg=numpy.average(betax,weights=Ds)
	betay_avg=numpy.average(betay,weights=Ds)
	#print "<betax> = ",betax_avg
	#print "<betay> = ",betay_avg
	#print "<acc> = ",numpy.average(acc,weights=Ds)#, "<accp> = ",numpy.average(accp,weights=Ds), "<accn> = ",numpy.average(accn,weights=Ds) 

	#Machine parameters
	C = param['C'] #Circunference (m)
	frev=param['C']/param['cluz'] #Rev. freq (Hz)	
	gamma=param['gamma']
	Np=param['Np']	
	theta_x = sqrt(param['Ax']/betax_avg)
	theta_y = sqrt(param['Ay']/betay_avg)
	R=theta_y/theta_x
	FR=pi+(R**2+1.0)*sin(2*atan(R))+2.0*(R**2.0-1.0)*atan(R)
	
	#Elastic Scattering Lifetime
	Telas=10.25*2.0*pi/FR*(param['En']/1.0e+09)**2*param['Ay']/(betay_avg*Pmed)
	
	#Inelastic Scattering Lifetime
	Tine=1/(0.0065*Pmed*log(1./(numpy.average(acc,weights=Ds))))
	
	#Touschek Lifetime Calculation
	sigx=numpy.zeros(len(twiss))
	sigy=numpy.zeros(len(twiss))
	epsilon=numpy.zeros(len(twiss))
	Tv=numpy.zeros(len(twiss))
	dsdT=numpy.zeros(len(twiss)-1)
	
	#bunch size and volume
	sigx=numpy.sqrt(betax*ex+(Dx*sigP)**2)
	sigy=numpy.sqrt(betay*ey+(Dy*sigP)**2)
	V=sigS*sigx*sigy
	
	# parameters
	Sx2=ex*betax
	factor=betax*Dpx+alphax*Dx
	A1=1.0/(4.0*sigP**2)+(Dx**2+factor**2)/(4.0*Sx2)
	B1=betax*factor/(2.0*Sx2)
	C1=betax**2/(4.0*Sx2)-B1**2/(4.0*A1)
	
	#Epsilon factors	
	#epsilonp=(2.0*numpy.sqrt(C1)/gamma*accp)**2
	#epsilonn=(2.0*numpy.sqrt(C1)/gamma*accn)**2
	epsilonp=(2.0*numpy.sqrt(C1)*accp)**2/(1957.0*param['En']/1.0e+09)**2
	epsilonn=(2.0*numpy.sqrt(C1)*accn)**2/(1957.0*param['En']/1.0e+09)**2
	
	
	#Lifetime
	Tvp=9.4718e-23*Np/(gamma**2)*1/(accp**3)*De(epsilonp)/V
	Tvn=9.4718e-23*Np/(gamma**2)*1/(accn**3)*De(epsilonn)/V
 	aux=0.5*(numpy.average(Tvp,weights=Ds)+numpy.average(Tvn,weights=Ds))*3600
 	Ttous=1.0/aux
	 
	return (Ttous,Tine,Telas)	
	
def De(e):
	
	out=numpy.zeros(len(e))
	for j in range(len(e)):
		(int1,err1)=quad(integrand1,e[j],numpy.inf)
		(int2,err2)=quad(integrand2,e[j],numpy.inf)
		out[j]=0.5*sqrt(e[j])*(-3.0*exp(-e[j])+e[j]*int1+int2*(3.0*e[j]-e[j]*log(e[j])+2.0))							
	return out
	
def integrand1(x):		
	return log(x)/x*exp(-x)

def integrand2(x):
	return exp(-x)/x	