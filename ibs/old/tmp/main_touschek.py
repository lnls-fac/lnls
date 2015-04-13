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


exi=numpy.arange(0.05,2.0,0.05)
exi=exi*1e-9
Npoints=len(exi)

#Define new vectors
eyi=param['ey0']
spi=param['sp0']
ssi=param['ss0']
I0=param['Np']/mA
LFtous=numpy.zeros(Npoints)
LFine=numpy.zeros(Npoints)
LFelas=numpy.zeros(Npoints)

#Average pressure in the machine
Pmed=1 #[nTorr]

print '-----------------------------------------------'
print 'Calculates IBS effects and Lifetime results'
for j in range(Npoints):
	#Uses the IBS results to calculates lifetimes
	(LFtous[j],LFine[j],LFelas[j])=Calc_Lifetime(Pmed,param,I0,twiss,exi[j],eyi,spi,ssi)
	print 'Bunch number = ', j+1
	print 'Ib = {0:4.1f} mA' .format(I0)
	print 'ex_fim = {0:0.3f} nm rad' .format(exi[j]*1e9)
	print 'ey_fim = {0:0.3f} pm rad' .format(eyi*1e12)
	print 'sp_fim = {0:0.3f} %' .format(spi*100)
	print 'ss_fim = {0:0.3f} mm or {1:0.3f} ps' .format(ssi*1e3,ssi/param['cluz']*1e12)
	print 'Lifetime results [h]: Touschek = {0:0.2f}, Elastic = {1:0.2f} and Inelastic  = {2:0.2f}'.format(LFtous[j],LFelas[j],LFine[j])

print '-----------------------------------------------'
print '\n'


#Saves final data
f=open('Touschek_scan.txt','w')
f.write('Initial parameters: \n')
f.write('I0[mA] \t ey[m rad] \t sp[%] \t ss[m]\n')
f.write(str('{0:0.5e}'.format(I0)) + '\t')
f.write(str('{0:0.5e}'.format(eyi)) + '\t')
f.write(str('{0:0.5e}'.format(spi)) + '\t')
f.write(str('{0:0.5e}'.format(ssi)) + '\t')
f.write('\n\n')
f.write('ex[mA] \t LFtous[h] \t LFelas[h] \t LFine [h]\n')
for j in range(Npoints):
	f.write(str('{0:0.4f}'.format(exi[j]*1e9)) + '\t')
	f.write(str('{0:0.2f}'.format(LFtous[j])) + '\t')
	f.write(str('{0:0.2f}'.format(LFelas[j])) + '\t')
	f.write(str('{0:0.2f}'.format(LFine[j])) + '\n')
f.close()
print('Data saved in file: Lifetime.txt')		
		



