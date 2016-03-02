#!/usr/bin/python2
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm, sigmaclip
import sys
import pyaneti as pti

#Read the file with all the python functions
execfile('todo-py.py')


#PREPATARION RV DATA

ylab = 'RV (km/s)'
units_ms = False
ktom = 1.0
Porb = 1.0
T0   = 0.0
mstar= 1.0
imcmc = 5000000
is_circular = False
chi_toler = 0.099
prec = 1e-5
tfactor = 10000
e0 = 0.3
w0 = np.pi

fname_rv = "corot1.dat"

#Read the input_rv.py to know the 
execfile('input_rv.py')

#Read the data file
time,rv,err,tspe = np.loadtxt(fname_rv,usecols=(0,1,2,3), \
  dtype={'names': ('time', 'rv', 'err','telescope'), \
	'formats': ('float', 'float', 'float', 'S1')}, \
	comments='#',unpack=True)

#Transform rv from km/s to m/s
if(units_ms):
	ktom = 1000
	rv=rv*ktom
	err=err*ktom
	ylab = 'RV (m/s)'


#These lists have lists with data for the different telescopes
time_all=[]
rv_all=[]
errs_all=[]

#Number of telescopes
nt = len(telescopes)

#tspe is the data with the telescopes read from the file
for i in range(0,nt):
	time_dum =[]
	rv_dum =[]
	errs_dum =[]
	if (len(telescopes[i]) == 0):
		print ("There is no data for %s"% telescopes[i])
	else:
		for j in range(0,len(tspe)):
			if (tspe[j] == telescopes[i]):
				time_dum.append(time[j])
				rv_dum.append(rv[j])
				errs_dum.append(err[j])
	#The *all variables are lists of lists, each list constains
	# a list with the data of each telescope
		time_all.append(time_dum)
		rv_all.append(rv_dum)
		errs_all.append(errs_dum)


#The mega* variables contains all telescope data
mega_rv = []
mega_time = []
mega_err  = []
tlab = []
for i in range(0,nt): 
	for j in range(0,len(rv_all[i])):
		tlab.append(i)
		mega_rv.append(rv_all[i][j])
		mega_time.append(time_all[i][j])
		mega_err.append(errs_all[i][j])


#PREPARATION TRANSIT DATA

fname_tr = "corot_transit.tbl"

skrow = 0
#Read the data file
dummyd,dummyf,flag = np.loadtxt(fname_tr,usecols=(2,9,10), comments='\\',unpack=True,skiprows=skrow)

#Let us take the good data with the flag
nobin_wflux = []
nobin_hdate = []
for i in range(0,len(flag)):
	if ( flag[i] == 0):
		nobin_wflux.append(dummyf[i])
		nobin_hdate.append(dummyd[i])

nbin = 16

#bin the data to do fastest test
hdate, err_hdate = bin_data(nobin_hdate,nbin)
wflux, errs = bin_data(nobin_wflux,nbin)

print errs[0], np.sqrt(wflux[0])

#THIS HAS TO BE DONE AUTOMATICALLY

plt.xlim(3217,3218)
plt.errorbar(hdate,wflux,errs)
plt.show()

ntr = 2

tls = [None]*ntr 

tls[0] = [3217.,3218.]
tls[1] = [3232.1,3233.1] 

#

#crash if you do not have more than one transit
if ( ntr < 2):
	print "you do not have enought transit data!"
	sys.exit("I crashed because I want more data!")
	
#Each element of these vectors will have the information
#of a given transit
xt= [None]*ntr	
yt= [None]*ntr	
et= [None]*ntr	

#Normalize all the transit independently
for i in range(0,ntr):
	xt[i],yt[i],et[i] = normalize_transit(hdate,wflux,errs,tls[i])

#Let us put together the information of all the arrays
megax = np.concatenate(xt)
megay = np.concatenate(yt)
megae = np.concatenate(et)


#PRIORS SECTION

#Let us initialize the MCMC Metropolis-Hasting algorithm
#Let us try to do a guess for the init values
k0vecmax = [None]*nt
k0vecmin = [None]*nt
k0 = 0.0
for i in range(0,nt):
	k0vecmin[i] = min(rv_all[i])
	k0vecmax[i] = max(rv_all[i])
	k0 = k0 + ( k0vecmax[i] - k0vecmin[i] ) / 2

k0 = k0 / nt
v0 = np.zeros(nt)

for i in range(0,nt):
	v0[i] = ( k0vecmin[i] + k0vecmax[i] ) / 2.0


T0 = min(xt[0]) + 0.5*(max(xt[0])-min(xt[0]))
P = 15.
e = 0.3
w = np.pi / 2.0
i = np.pi/2
a = 13.0
u1 = 0.42
u2 = 0.25
pz = 0.5
prec = 5.e-5
maxi = int(1e8)
chi2_toler = 0.3
thin_factor = int(2e3)
ics = False
nconv = 100

fit_t0 	= True
fit_P 	= True
fit_e 	= True
fit_w 	= True
fit_i 	= True
fit_a 	= True
fit_u1 	= False
fit_u2 	= False
fit_pz 	= True
fit_k		= True
fit_v0	= True


#Transform the logical values to integers to select what
# to fit
what_fit = [int(fit_t0), int(fit_P), int(fit_e), int(fit_w), \
            int(fit_i),int(fit_a), int(fit_u1),int(fit_u2),  \
            int(fit_pz), int(fit_k), int(fit_v0)]

dummy = [T0,P,e,w,i,a,u1,u2,pz,k0]
params = np.concatenate((dummy,v0))

#Start the powerful mcmc routine to do the fit
#pti.metropolis_hastings_tr(megax, megay, megae, T0, P, e, \
#	 w, i, a, u1, u2, pz, prec, maxi, thin_factor, ics,     \
#	 what_fit, nconv)

pti.metropolis_hastings_tr_improved(megax, megay, megae,  \
	params, prec, maxi, thin_factor, ics, what_fit, nconv)

vari, chi2,chi2red,t0o,Po,eo,wo,io,ao,u1o,u2o,pzo = \
			 np.loadtxt('mh_trfit.dat', comments='#',unpack=True)

#Read the output values from the file mh_trfit.dat
#vari,chi2,chi2red,eo, wo, io, ao, u1o, u2o, pzo, t0o, Po = \
#			 np.loadtxt('mh_trfit.dat', comments='#',unpack=True)

#Find the values and errors
t0_val,t0_err = find_vals_gauss(t0o,nconv)
P_val, P_err  = find_vals_gauss(Po,nconv)
e_val,e_err 	= find_vals_gauss(eo,nconv)
w_val,w_err 	= find_vals_gauss(wo,nconv)
i_val,i_err 	= find_vals_gauss(io,nconv)
a_val,a_err 	= find_vals_gauss(ao,nconv)
u1_val,u1_err = find_vals_gauss(u1o,nconv)
u2_val,u2_err = find_vals_gauss(u2o,nconv)
pz_val,pz_err = find_vals_gauss(pzo,nconv)

w_deg = w_val * 180. / np.pi
w_deg_err = w_err * 180. / np.pi
i_deg = i_val * 180. / np.pi
i_deg_err = i_err * 180. / np.pi

#Print the best fit values values
print ('The best fit planet parameters are:')
print ('T0    = %4.4e +/- %4.4e days'%(t0_val,t0_err))
print ('e     = %4.4e +/- %4.4e'			%(e_val,e_err))
print ('w     = %4.4e +/- %4.4e deg'	%(w_deg,w_deg_err))
print ('i     = %4.4e +/- %4.4e deg' %(i_val,i_err))
print ('a/r*  = %4.4e +/- %4.4e' 		%(a_val,a_err))
print ('u1	  = %4.4e +/- %4.4e' 		%(u1_val,u1_err))
print ('u2    = %4.4e +/- %4.4e' 		%(u2_val,u2_err))
print ('rp/r* = %4.4e +/- %4.4e' 		%(pz_val,pz_err))
print ('P     = %4.4e +/- %4.4e' 		%(P_val,P_err))


#Let us obtain the residuals to do a nice plot

#Move all the points to T0
for i in range(0,ntr):
	xt[i] = xt[i] - P_val * i

#Redefine megax with the new xt valies
megax = np.concatenate(xt)
z_val = pti.find_z(megax,t0_val,P_val,e_val,w_val,i_val,a_val)
mud_val, mu0_val = pti.occultquad(z_val,u1_val,u2_val,pz_val)
#Residuals
res = megay - mud_val

#Get the model data to do the plot
nvec = int(1e5)
dx = ( max(megax) - min(megax) ) / nvec
xvec = np.zeros(nvec)
xvec[0] = min(megax) 
for i in range(1,nvec):
	xvec[i] = xvec[i-1] + dx
zvec = pti.find_z(xvec,t0_val,P_val,e_val,w_val,i_val,a_val)
mud, mu0 = pti.occultquad(zvec,u1_val,u2_val,pz_val)
#Now we have data to plot a nice model

#Do the plot
plt.figure(2,figsize=(10,10))
#Plot the transit light curve
plt.subplot(211)
plt.xlim(min(xt[0]),max(xt[0]))
plt.errorbar(megax,megay,megae,fmt='o',alpha=0.3)
plt.plot(xvec,mud,'k',linewidth=2.0)
#Plot the residuals
plt.subplot(212)
plt.xlim(min(xt[0]),max(xt[0]))
plt.errorbar(megax,res,megae,fmt='o',alpha=0.3)
plt.plot(megax,np.zeros(len(megax)),'k--',linewidth=2.0)
plt.show()

#The and of transit_fmcmc.py
