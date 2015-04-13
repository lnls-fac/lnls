#!/usr/bin/python

import sys
import time
import threading
import numpy
import string
import copy
from scipy.interpolate import interp2d
from scipy.optimize import curve_fit
from math import sqrt,exp,log,pi,acos,atan,cos

def Read_data(ftwiss,facc,param,HC):

	#read twiss file
	tfile = open(ftwiss,'rU')
	twiss=numpy.zeros(12)
	newrow=numpy.zeros(12)
	element=numpy.zeros(1)
	for aline in tfile:
		values=aline.split()
		if (values[0]=='@'):
			if (values[1]=='K_beta'):
				param['k_beta']= float(values[3])
			elif (values[1]=='K_dw'):
				param['k_dw']=float(values[3])-param['k_beta']
			elif (values[1]=='EX'):
				param['ex0']=float(values[3])
		elif (values[0]<>'@' and values[0]<>'*'):
			if float(values[3])>0:
				for i in range(10):
					newrow[i]=float(values[i+2])
				twiss=numpy.vstack((twiss,newrow))
				element=numpy.vstack((element,values[0]))
	tfile.close()
	
	twiss=numpy.delete(twiss,0,axis=0)
	element=numpy.delete(element,0,axis=0)
	
	calc_emitt(param,HC)
	
	
	#Calculates RF acceptance
	U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+9
	q=param['Vrf']/U0
	eRF=sqrt(2*U0/(pi*param['ap']*param['hh']*param['En'])*(sqrt(q**2-1)-acos(1/q)))
		
	#read Acceptance file
	tfile = open(facc,'rU')
	accp=numpy.zeros(2)
	accn=numpy.zeros(2)
	newrow=numpy.zeros(2)
	for aline in tfile:
		values=aline.split()
		newrow[0]=float(values[0])
		newrow[1]=min(abs(float(values[1])),eRF)
		accp=numpy.vstack((accp,newrow))
		newrow[1]=min(abs(float(values[2])),eRF)
		accn=numpy.vstack((accn,newrow))
	tfile.close()
	
	accp=numpy.delete(accp,0,axis=0)
	accn=numpy.delete(accn,0,axis=0)
	
	#increase acc for the ring length
	newrow2=numpy.zeros((len(accp),2))
	newrow2=copy.copy(accp)
	last=accp[-1,0]
	rowmax=len(accp)

	for i in range(1,10):
		newrow2[:,0]=i*last+copy.copy(accp[0:rowmax,0])
		accp=numpy.vstack((accp,newrow2))
		
	newrow2=numpy.zeros((len(accn),2))
	newrow2=copy.copy(accn)
	last=accn[-1,0]
	rowmax=len(accn)

	for i in range(1,10):
		newrow2[:,0]=i*last+copy.copy(accn[0:rowmax,0])
		accn=numpy.vstack((accn,newrow2))
	
	#interpolate acceptance in along the twiss positions
	accptwiss=numpy.zeros(len(twiss))
	accptwiss=numpy.interp(twiss[:,0],accp[:,0],accp[:,1])
	accntwiss=numpy.zeros(len(twiss))
	accntwiss=numpy.interp(twiss[:,0],accn[:,0],accn[:,1])
	
	#concatenate the acceptance columns in the twiss matrix
	twiss[:,10]=accptwiss
	twiss[:,11]=accntwiss
		
	return twiss

def Parameters_setup(fname):

	#Read file with machine parameters
	tfile = open(fname,'rU')
	param={}
	for aline in tfile:
		if aline[0]<>'#':
			values=aline.split()
			param[values[0]]=float(values[1])
		
	tfile.close()
	return param
	
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


def calc_emitt(param,HC):
	
	
	#Calculate energy spread
	sigp=sqrt(param['Cq']*param['gamma']**2*(param['I3']/(2*param['I2']+param['I4'])))
	param['sp0']=sigp
	
	#Calculate bunch length
	U0=param['Cgamma']/(2*pi)*(param['En']/1e+9)**4*param['I2']*1e+9
	if (HC==0):
		param['ss0']=param['sp0']*param['C']*sqrt(param['ap']*param['En']/(2*pi*param['hh']*(param['Vrf']**2-U0**2)**0.5))
	elif (HC==1):
		Qs0=sqrt(param['ap']*param['hh']*sqrt(param['Vrf']**2-U0**2)/(2*pi*param['En']))
		param['ss0']=0.432343807*(param['ap']*param['hh']*pi*sigp/Qs0)**(0.5)*param['C']/(2*pi*param['hh'])

			
	#Calculate ex
	Jx=1-param['I4']/param['I2']
	ex=param['Cq']*param['gamma']**2*param['I5']/(Jx*param['I2'])
	param['ex0']=ex		
	
	#Calculate ey
	ey=(param['k_beta']+param['k_dw'])*ex	
	param['ey0']=ey	
	
	mA=1e-3/param['Ccoul']*param['C']/param['cluz']
	
	print '-----------------------------------------------'

	print 'Initial beam parameters'
	print 'Itot = {0:4.1f} mA' .format(param['Np']*param['hh']/mA)
	print 'ex_ini = {0:0.3f} nm rad' .format(param['ex0']*1e9)
	print 'ey_ini = {0:0.3f} pm rad' .format(param['ey0']*1e12)
	print 'sp_ini = {0:0.3f} %' .format(param['sp0']*100)
	print 'ss_ini = {0:0.3f} mm' .format(param['ss0']*1e3)
	print '-----------------------------------------------'
	print '\n'


