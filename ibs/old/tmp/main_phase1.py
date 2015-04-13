#!/usr/bin/env python

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
param=Parameters_setup('Parameters_Sirius_Phase1.txt') 
#reads the twiss and the momentum acceptance file and interpolate the momentum aperture (in case the number of points don;t match the twiss file)
#using the parameter list calculates the initial conditions of the beam (intial bunch curent and 6D emittance)
# The order of the parameters in the twiss files is: name,KEYWORD,s,L,betx,alfx,dx,dpx,bety,alfy,dy,dpy
# The order in the acceptance file is: s, positive acceptance values, begative acceptance values
twiss=Read_data('twiss_Sirius.txt','acc_Sirius_V08.txt',param,0) 

#Definition of mA
mA=1e-3/param['Ccoul']*param['C']/param['cluz']


I0=[0.116]#numpy.arange(0.1,4.0,0.05)
Npoints=len(I0)

#Define new vectors
Np=numpy.zeros(Npoints)
exi=numpy.zeros(Npoints)
eyi=numpy.zeros(Npoints)
spi=numpy.zeros(Npoints)
ssi=numpy.zeros(Npoints)
LFtous=numpy.zeros(Npoints)
LFine=numpy.zeros(Npoints)
LFelas=numpy.zeros(Npoints)

#Average pressure in the machine
Pmed=1 #[nTorr]

print '-----------------------------------------------'
print 'Calculates IBS effects and Lifetime results'
for j in range(len(I0)):
	#Redefines parameter Np
	param['Np']=I0[j]*mA
	#Calculates IBS effects on the emittances
	(exi[j],eyi[j],spi[j],ssi[j])=Iterate_emittances(twiss,param)	
	#Uses the IBS results to calculates lifetimes
	(LFtous[j],LFine[j],LFelas[j])=Calc_Lifetime(Pmed,param,I0[j],twiss,exi[j],eyi[j],spi[j],ssi[j])
	print 'Bunch number = ', j+1
	print 'Ib = {0:4.1f} mA' .format(I0[j])
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
		



