#!/usr/bin/python

import sys
import time
import threading
import numpy
import string
import copy
from scipy.optimize import curve_fit
from math import sqrt,exp,log,pi,acos

#import other modules for calculation of IBS, lifetime and Landau Cavity tracking and Input File reading
from Read_input import *
from CavLandau import *
from BunchLength import *
from CIMP_3HC import *
from Lifetime import *

#Load and define parameters 
#reads the parameter file and creates a dictionary with the parameters
param=Parameters_setup('Parameters_Sirius_Phase2.txt') 
#reads the twiss and the momentum acceptance file and interpolate the momentum aperture (in case the number of points don;t match the twiss file)
#using the parameter list calculates the initial conditions of the beam (intial bunch curent and 6D emittance)
# The order of the parameters in the twiss files is: name,KEYWORD,s,L,betx,alfx,dx,dpx,bety,alfy,dy,dpy
# The order in the acceptance file is: s, positive acceptance values, begative acceptance values
twiss=Read_data('twiss_Sirius.txt','acc_Sirius_V08.txt',param,1) 

#Definition of mA
mA=1e-3/param['Ccoul']*param['C']/param['cluz']

# Total stored current
Itot=param['Np']*param['hh']/mA

#this quatity is a multidimensional array to create the filling patter for the machine
#for example an evenlly filled machine has gap=[[number of buckets,1]], meaning all the buncketshave the same intensity
#a machine swith a gap of 50 buckets and a single bunch in the middle is gap=[[buckets-50,1],[25,0],[1,5],[24,0]], meaning we have (buckets-50) bunches
#with even intensity, 24 empty buckets, a single bunch with a intensity 5 times higher than the other buckets and## finally 24 empty buckets.
gap=[[864,1]]
#Harmonic cavity resont harmonic 
frf=param['cluz']/param['C']*param['hh']
hres=param['mharm']*param['hh']+param['DetuneHC']*param['hh']/frf

#Traking calculation, gives back 5 vectors (or arrays) of numbers with
#	I0 - current per bunch (A)
# 	Vh - harminc voltage (kV)
# 	ph - harminc phase (rad)
# 	Vm - accelerating voltage (kV)
# 	pm - accelerating pahse phase (rad), difference to the theoretical sincrotron phase
#You can continue to use these vectors or load everything form the saved file: Landau.txt
#(Vh,Vm,ph,pm,I0)=Phase_Shift_Sirius(500000,Itot,gap,hres,param,10000)

Npoints=1
I0=numpy.zeros(Npoints)
Vh=numpy.zeros(Npoints)
ph=numpy.zeros(Npoints)
Vm=numpy.zeros(Npoints)
pm=numpy.zeros(Npoints)

#Open file with bunch's properties (results from tracking)
#	I0 - current per bunch (A)
# 	Vh - harminc voltage (kV)
# 	ph - harmonic phase (rad)
# 	Vm - accelerating voltage (kV)
# 	pm - accelerating pahse phase (rad), difference to the theoretical sincrotron phase
bfile=open('Landau_nogap.txt','r')
values=bfile.readline().split('\t')#read title's line
for i in range(Npoints):
	values=bfile.readline().split('\t')
#	I0[i]=float(values[0])*1e3 #convert to mA
#	Vh[i]=float(values[1])
#	ph[i]=float(values[2])
#	Vm[i]=float(values[3])
#	pm[i]=float(values[5])
	I0[i]=Itot/864#float(values[2])*1e3 #convert to mA
	Vh[i]=float(values[1])
	ph[i]=float(values[5])
	Vm[i]=float(values[0])
	pm[i]=float(values[4])
bfile.close()


#Define new vectors
Np=numpy.zeros(Npoints)
exi=numpy.zeros(Npoints)
eyi=numpy.zeros(Npoints)
spi=numpy.zeros(Npoints)
ssi=numpy.zeros(Npoints)
LFtous=numpy.zeros(Npoints)
LFine=numpy.zeros(Npoints)
LFelas=numpy.zeros(Npoints)

#Define number of particles/bunch (or filling pattern array)
Np=I0*mA

#Calculates the bunch length and saves data in the file BunchLength.txt
#Also prints some statistics at the end
Calc_Sigma(Vh,ph,Vm,pm,I0,Npoints,param)

#Average pressure in the machine
Pmed=1 #[nTorr]

#Check places where I0 differnte then and get the indexes
not_zero=numpy.flatnonzero(I0)

print '-----------------------------------------------'
print 'Calculates IBS effects and Lifetime results'
for j in range(len(not_zero)):
	#Redefines parameter Np
	param['Np']=I0[not_zero[j]]*mA
	#Calculates IBS effects on the emittances
	(exi[j],eyi[j],spi[j],ssi[j])=Iterate_emittances3HC(twiss,param,pm[not_zero[j]],Vm[not_zero[j]],ph[not_zero[j]],Vh[not_zero[j]])	
	#Uses the IBS results to calculates lifetimes
	(LFtous[j],LFine[j],LFelas[j])=Calc_Lifetime(Pmed,param,I0[not_zero[j]],twiss,exi[j],eyi[j],spi[j],ssi[j])
	print 'Bunch number = ', not_zero[j]+1
	print 'Ib = {0:4.1f} mA' .format(I0[not_zero[j]])
	print 'ex_fim = {0:0.3f} nm rad' .format(exi[j]*1e9)
	print 'ey_fim = {0:0.3f} pm rad' .format(eyi[j]*1e12)
	print 'sp_fim = {0:0.3f} %' .format(spi[j]*100)
	print 'ss_fim = {0:0.3f} mm or {1:0.3f} ps' .format(ssi[j]*1e3,ssi[j]/param['cluz']*1e12)
	print 'Lifetime results [h]: Touschek = {0:0.2f}, Elastic = {1:0.2f} and Inelastic  = {2:0.2f}'.format(LFtous[j],LFelas[j],LFine[j])

print '-----------------------------------------------'
print '\n'


#Saves final data
f=open('Lifetime.txt','w')
f.write('Initial parameters: \n')
f.write('ex[m rad] \t ey[m rad] \t sp \t\t ss[m]\n')
f.write(str('{0:0.5e}'.format(param['ex0'])) + '\t')
f.write(str('{0:0.5e}'.format(param['ey0'])) + '\t')
f.write(str('{0:0.5e}'.format(param['sp0'])) + '\t')
f.write(str('{0:0.5e}'.format(param['ss0'])) + '\t')
f.write('\n\n')
f.write('Ib[mA] \t LFtous[h] \t LFelas[h] \t LFine [h] \t ex[nm rad] \t ey[pm rad] \t sp[%] \t ss[mm]\n')
for j in range(Npoints):
	f.write(str('{0:0.4f}'.format(I0[j])) + '\t')
	f.write(str('{0:0.2f}'.format(LFtous[j])) + '\t')
	f.write(str('{0:0.2f}'.format(LFelas[j])) + '\t')
	f.write(str('{0:0.2f}'.format(LFine[j])) + '\t')
	f.write(str('{0:0.5e}'.format(exi[j]*1e9)) + '\t')
	f.write(str('{0:0.5e}'.format(eyi[j]*1e12)) + '\t')
	f.write(str('{0:0.5e}'.format(spi[j]*100)) + '\t')
	f.write(str('{0:0.5e}'.format(ssi[j]*1e3)) + '\n')
f.close()
print('Data saved in file: Lifetime.txt')		
		



