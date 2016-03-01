#!/usr/bin/python2

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm, sigmaclip
import sys
import pyaneti as pti

#----- HERE THE MAIN PROGRAM STARTS -----#

fname = "corot_transit.tbl"
 
skrow = 0
#Read the data fid
dummyd,dummyf,flag = np.loadtxt(fname,usecols=(2,9,10), comments='\\',unpack=True,skiprows=skrow)

#Let us take the good data with the flag
wflux = []
hdate = []
for i in range(0,len(flag)):
	if ( flag[i] == 0):
		wflux.append(dummyf[i])
		hdate.append(dummyd[i])

errs = np.sqrt(wflux)

def find_transits(x,y):

	newy = sigmaclip(y,low=3.0,high=3.0)

	transits, miny, maxy = y - newy	

#----------------------------------------------

ntr = 2

tls = [None]*ntr 

tls[0] = [3217,3218]
tls[1] = [3232.1,3233.1] 

def parabola(x,a,b,c):
	y = a + b*x +c*x*x 
	return y 

def normalize(x,y,err,limits):
	dummyx  = []
	dummyy  = []
	dummyerr= []
	newx = []
	newy = []
	newerr = []
	for i in range(0,len(x)):
		if ( x[i] > limits[0] and x[i] < limits[1] ):
			dummyy.append(y[i])
			dummyx.append(x[i])
			dummyerr.append(err[i])

	newy = sigmaclip(dummyy,low=3.0,high=3.0)
	newy = sigmaclip(newy[0],low=2.0,high=2.0)
	for i in range(0,len(dummyx)):
		if ( dummyy[i] > newy[1] and dummyy[i] < newy[2] ):
			newx.append(dummyx[i])
			newerr.append(dummyerr[i])

	popt, cov = curve_fit(parabola,newx,newy[0],sigma=newerr)
	dummyx = np.array(dummyx)
	par = parabola(dummyx,popt[0],popt[1],popt[2])
	dummyy = dummyy / par
	dummyerr = dummyerr / par
	return dummyx, dummyy, dummyerr
	
	#Now we have  the values around the transit
	
xt= [None]*ntr	
yt= [None]*ntr	
et= [None]*ntr	
#
xt[0],yt[0], et[0] = normalize(hdate,wflux,errs,tls[0])
xt[1],yt[1], et[1] = normalize(hdate,wflux,errs,tls[1])

if ( ntr < 2):
	print "you do not have enought transit data!"
	sys.exit("I could crash!")

#Let us start to make a good priors calculation
#T0 = [None]*ntr
#for i in range(0,ntr):
#	T0[i] = min(xt[i]) + 0.5*(max(xt[i])-min(xt[i]))
T0 = min(xt[0]) + 0.5*(max(xt[0])-min(xt[0]))

#
nx = [None]*ntr
ny = [None]*ntr
ne = [None]*ntr
dx = []
dy = []
de = []
for i in range(0,ntr):
	for j in range(0,len(xt[i])):
		if ( j % 1 == 0):
			dx.append(xt[i][j])
			dy.append(yt[i][j])
			de.append(et[i][j])
	nx[i] = dx
	ny[i] = dy
	ne[i] = de
	dx = []
	dy = []
	de = []

xt = []
yt = []
et = []

xt = nx
yt = ny
et = ne

#

tlimits = [None]*(ntr+1)

tlimits[0] = int(0)
for i in range(1,ntr+1):
	tlimits[i] = tlimits[i-1] + len(xt[i-1])

megax = np.concatenate(xt)
megay = np.concatenate(yt)
megae = np.concatenate(et)

P = max(megax) - min(megax)
P = 15.
e = 0.3
w = np.pi / 2.0
i = np.pi/2
a = 16.0
u1 = 0.42
u2 = 0.25
pz = 0.2
prec = 2.5e-3
maxi = int(1e8)
chi2_toler = 0.3
thin_factor = int(2e3)
ics = False
nconv = 100

fit_e = False
fit_w = True
fit_i = True
fit_a = True
fit_u1 = False
fit_u2 = False
fit_pz = True
fit_t0 = True
fit_P = True

what_fit = [int(fit_e),int(fit_w),int(fit_i),int(fit_a), \
					  int(fit_u1),int(fit_u2),int(fit_pz),int(fit_t0), \
						int(fit_P)]

#pti.metropolis_hastings_tr(megax, megay,megae, T0,P,e, w, i, a, u1, u2, pz,tlimits,prec, maxi, thin_factor, chi2_toler, ics, what_fit, nconv)

vari,chi2,chi2red,eo, wo, io, ao, u1o, u2o, pzo,t0o,Po = np.loadtxt('mh_trfit.dat', comments='#',unpack=True)

#t0o = [None]*ntr

#for j in range(0,ntr):
#	n = [10+j]
#	t0o[j] = np.loadtxt('mh_trfit.dat', comments='#',unpack=True, usecols=(n))



def find_errb(x,nconv):
	#Let us take only the converging part
	iout = len(x) - nconv
	xnew = x[iout:]
	mu,std = norm.fit(xnew)
	return mu, std

#t0_val = [None]*ntr
#t0_err = [None]*ntr
#for i in range(0,ntr):
#	t0_val[i],t0_err[i] = find_errb(t0o[i],nconv)



e_val,e_err = find_errb(eo,nconv)
w_val,w_err = find_errb(wo,nconv)
i_val,i_err = find_errb(io,nconv)
a_val,a_err = find_errb(ao,nconv)
u1_val,u1_err = find_errb(u1o,nconv)
u2_val,u2_err = find_errb(u2o,nconv)
pz_val,pz_err = find_errb(pzo,nconv)
t0_val,t0_err = find_errb(t0o,nconv)
P_val, P_err  = find_errb(Po,nconv)

print ('The planet parameters are:')
print ('T0= %4.4e +/- %4.4e' %(t0_val,t0_err))
print ('e = %4.4e +/- %4.4e' %(e_val,e_err))
print ('w = %4.4e +/- %4.4e' %(w_val,w_err))
print ('i = %4.4e +/- %4.4e' %(i_val,i_err))
print ('a = %4.4e +/- %4.4e' %(a_val,a_err))
print ('u1= %4.4e +/- %4.4e' %(u1_val,u1_err))
print ('u2= %4.4e +/- %4.4e' %(u2_val,u2_err))
print ('pz= %4.4e +/- %4.4e' %(pz_val,pz_err))
print ('P = %4.4e +/- %4.4e' %(P_val,P_err))

#x_dum = [None]*ntr
#for i in range(0,ntr):
#	x_dum[i] = xt[i] - t0_val[i]

#New megax
#megax = np.concatenate(x_dum)

#Plot the initial conditions
zvec = pti.find_z(megax,t0_val,P_val,e_val,w_val,i_val,a_val, tlimits)

mud, mu0 = pti.occultquad(zvec,u1_val,u2_val,pz_val)

res = megay - mud

plt.figure(2,figsize=(10,10))
plt.subplot(211)
plt.xlim(min(xt[0]),max(xt[0]))
plt.plot(linewidth=2.0)
plt.plot(megax,megay,'bo',alpha=0.3)
plt.plot(megax,mud,'k')
plt.subplot(212)
plt.xlim(min(xt[0]),max(xt[0]))
plt.plot(linewidth=2.0)
plt.plot(megax,res,'bo',alpha=0.3)
plt.plot(megax,np.zeros(len(megax)),'k--')
plt.show()

