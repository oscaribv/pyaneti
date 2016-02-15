
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
import sys
import fmcmc

#-------------------------------------
	#This functions gives the rv functions
	#it assumes a circular orbit
def rv_circular(t,rv0,k):
	rv = rv0 - k * np.sin( 2.0 * np.pi * ( t - T0 ) / Porb )
	return rv
#------------------------------------
def rv_curve_rv0(t,rv0,k0,ecc,omega):
	#man is the mean anomally
	man = 2 * np.pi * (t - T0) / Porb
	#call to find_anomaly function to find the eccentric anomaly
	anomaly = fmcmc.find_anomaly(man,ecc,1e-4,1e6)
	#The general equation for rv motion is
	rv = rv0 + k0 * (np.cos( anomaly + omega ) + ecc*np.cos( omega ) ) / np.sqrt(1.0 - ecc*ecc)
	return rv
#Linearized model from:
#http://adsabs.harvard.edu/abs/2009ApJS..182..205W
def rv_curve_wh(t,h,c,ecc,v0,tp):
  man = 2 * np.pi * ( t - tp ) / Porb
  anomaly = fmcmc.find_anomaly(man,ecc,1e-4,1e6)
  ut = h * np.cos(anomaly) + c * np.sin(anomaly) + v0
  return ut
#------------------------------------
def rv_curve(t,k0,ecc,omega):
	#man is the mean anomally
	man = 2 * np.pi * ( t - T0 ) / Porb
	#call to find_anomaly function to find the eccentric anomaly
	anomaly = fmcmc.find_anomaly(man,ecc,1e-4,1e6)
	#The general equation for rv motion is
	rv = + k0 * (np.cos( anomaly + omega ) + ecc*np.cos( omega ) )  / np.sqrt(1.0 - ecc*ecc)
	return rv
#------------------------------------
#Extract systemic velocities is highly important
#so let us take care with this
def find_rv0(time,fase,err,tpe):
	#Let us first fit assuming a circular orbit 
	#(to first estimate v0 and k0)
	popt,pcov = curve_fit(rv_circular,time,fase,sigma=err)
	#These results will be our first guesses for v0 and k0
	rvc = popt[0]
	sdrvc= np.sqrt(pcov[0][0])
	kc  = popt[1]
	#Now let us fit again now with the guesses as input
	popt,pcov = curve_fit(rv_curve_rv0,time,fase,sigma=err,p0=[rvc,kc,0,0])
	#Let us check if the function is "ellipsable"
	#If we have two points, try to fit an elliptical orbit is not correct
	#If this is the case, the errors wll go to infinite, let us check this
	if ( np.isinf(pcov[0][0]) == True ): #Then the function is not ellipsable
		rv0 = rvc
		sdrv0= sdrvc
	else: #The function is ellipsable, let us keep playing with the fit!
		rv0 = popt[0]
		sdrv0 = pcov[0][0]
		k0 = popt[1]
		w = popt[2]
		#Now we have guesses for rv0, k0 and w, let us re do the fit, now to improve e
		popt,pcov = curve_fit(rv_curve_rv0,time,fase,sigma=err,p0=[rv0,k0,0,w])
		rv0= popt[0]
		sdrv0= np.sqrt(pcov[0][0])
	print ('for %1s -> rv0 = %5.5f +/- %5.5f m/s'%(tpe,rv0,sdrv0))
	return rv0
#-----------------------------
def scale_period(jd,T0,P):
	x = [None]*len(jd)
	for i in range(len(jd)):
		x[i] = ( ( jd[i] - T0 ) % P ) /  P
	return x
#-----------------------------
def planet_mass(mstar,k,P,ecc):
	#Gravitational costant
	Gc = 6.67408e-11 #m^3 / (kgs^2)
	Gc = Gc * 1.989e30 # m^3 / (Msun s^2)
	#period in seconds
	P = P * 24. * 3600
	mpsin = k * ( 2. * np.pi * Gc / P)**(-1./3.)  * \
	mstar**(2./3.) * np.sqrt( 1.0 - ecc*ecc )
	return mpsin
#-----------------------------

#----- HERE THE MAIN PROGRAM STARTS -----#

ylab = 'RV (km/s)'
ktom = 1.0
Porb = 1.0
T0   = 0.0
mstar= 1.0
extract_rv = False
units_ms = False
imcmc = 5000000

#Read the input_rv.py to know the 
execfile('input_rv.py')

#Read the data file
time,fase,err,tspe = np.loadtxt(fname,usecols=(0,1,2,3), \
  dtype={'names': ('time', 'fase', 'err','telescope'), \
	'formats': ('float', 'float', 'float', 'S1')}, \
	comments='#',unpack=True)

#Transform fase from km/s to m/s
if(units_ms):
	ktom = 1000
	fase=fase*ktom
	err=err*ktom
	ylab = 'RV (m/s)'
	

#These lists have lists with data for the different telescopes
time_all=[]
fase_all=[]
errs_all=[]

#Number of telescopes
nt = len(telescopes)

#tspe is the data with the telescopes read from the file
for i in range(0,nt):
	time_dum =[]
	fase_dum =[]
	errs_dum =[]
	if (len(telescopes[i]) == 0):
		print ("There is no data for %s"% telescopes[i])
	else:
		for j in range(0,len(tspe)):
			if (tspe[j] == telescopes[i]):
				time_dum.append(time[j])
				fase_dum.append(fase[j])
				errs_dum.append(err[j])
	#The *all variables are lists of lists, each list constains
	# a list with the data of each telescope
		time_all.append(time_dum)
		fase_all.append(fase_dum)
		errs_all.append(errs_dum)

if (extract_rv):
	print ("Extracting systemic velocities for the %i telescopes"%nt)
	#Find all the systematic velocities and put them in rv0_all
	rv0_all=[None]*nt
	for i in range(0,nt):
		rv0_all[i] = find_rv0(time_all[i],fase_all[i],errs_all[i],telescopes[i])

	#Let us take off the offset for each telescope
	for i in range(0,nt): #each i is a different telescope
		for j in range(0,len(fase_all[i])):
			fase_all[i][j] = fase_all[i][j] - rv0_all[i]


#The mega* variables contains all telescope data
mega_fase = []
mega_time = []
mega_err  = []
tlab = []
for i in range(0,nt): 
	for j in range(0,len(fase_all[i])):
		tlab.append(i)
		mega_fase.append(fase_all[i][j])
		mega_time.append(time_all[i][j])
		mega_err.append(errs_all[i][j])

#----------------------------------------------
#
# MCMC calculation starts here
#
#----------------------------------------------

def find_errb(x):
	iout = len(x)/5 * 4
	#Let us take only the converging part
	xnew = x[iout:]
	mu,std = norm.fit(xnew)
	return mu, std 


#Let us initialize the MCMC Metropolis-Hasting algorithm
#Let us try to do a guess for the init values
v0 = np.zeros(nt)
kmin = min(mega_fase)
kmax = max(mega_fase)
k0 = (kmax - kmin) /  2.0
e0 = 0.1
w0 = np.pi
prec = 1e-3
chi_limit = 1.0

for i in range(0,nt):
	v0[i] = k0 + kmin

#Start the FORTRAN calling

ndata = 1000

#vmc = np.empty(ndata,nt)
kmc = np.empty(ndata)
emc = np.empty(ndata)
wmc = np.empty(ndata)
t0mc= np.empty(ndata)
pmc = np.empty(ndata)

print tlab

#sys.exit()

vmc, kmc, emc, wmc, t0mc, pmc = fmcmc.mcmc_rv(mega_time,mega_fase,mega_err,tlab,v0,k0,e0,w0,T0,Porb,prec,imcmc,ndata,chi_limit)
#fmcmc.mcmc_rv(mega_time,mega_fase,mega_err,v0,k0,e0,w0,T0,Porb,prec,imcmc,ndata)

#end fortran calling	

#plt.figure(1,figsize=(8,4))
#plt.xlabel("Iteration")
#plt.ylabel("Value")
#plt.plot(vmc,label='v0')
#plt.plot(kmc,label='k')
#plt.plot(emc,label='e')
#plt.plot(wmc,label='w')
#plt.plot(t0mc,label='w')
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
#           ncol=4, mode="expand", borderaxespad=0.)
#plt.show()

kval,sigk  = find_errb(kmc) 
ecval,sige = find_errb(emc)
wval,sigw  = find_errb(wmc)
tval,sigt  = find_errb(t0mc)
pval,sigp  = find_errb(pmc)

vval = [None]*nt
sigv = [None]*nt

for i in range(0,nt):
	vval[i], sigv[i]  = find_errb(vmc[i])

#if eccentricity is negative, there is a shift in omega
if (ecval < 0.0 ):
	ecval = - ecval
	wval  = wval + np.pi

#----------------------------------------------
#
# MCMC calculation finishes here
#
#----------------------------------------------

#Now let us calculate the plantary mass
#e = 0.0
k = kval

#ecval = 0.0

mpsin = planet_mass(mstar,k*ktom,Porb,ecval)
mjup = 0.0009543
mearth = 0.000003003 

#Print the results
print (' k = %4.4e +/- %4.4e m/s' %(k,sigk))
print (' w = %4.4e +/- %4.4e rad' %(wval,sigw))
print (' e = %4.4e +/- %4.4e    ' %(ecval,sige))
print (' t0= %4.4e +/- %4.4e    ' %(tval,sigt))
print (' P = %4.4e +/- %4.4e    ' %(pval,sigp))
for i in range(0,nt):
	print ('For %s, v = %4.4e +/- %4.4e m/s' %(tspe[i],vval[i],sigv[i]))

print ('Planet mass of %1.4e M_j (for a %2.2f solar mass star)' % (mpsin/mjup, mstar))

#Let us do a nice plot

T0 = tval

#Create the RV fitted curve
n = 500
xmin = T0
xmax = T0 + Porb
dn = (xmax - xmin) /  n  
rvx = np.empty([n])
rvx[0] = xmin 
for i in range(1,n):
	rvx[i] = rvx[i-1] + dn
#rvy = rv_curve(rvx,k,e,w)
#rvy = rv_circular(rvx,vval,kval)
#rvy = rv_curve_rv0(rvx,vval,kval,ecval,wval)
rvy = fmcmc.rv_curve(rvx,0.0,T0,kval,Porb,ecval,wval)
#rvy = rv_curve_wh(rvx,hval,cval,ecval,vval,tval)

res = [None]*nt
for i in range(0,nt):
	res[i] = fmcmc.rv_curve(time_all[i],0.0,T0,kval,Porb,ecval,wval)
	fase_all[i] = fase_all[i] - vval[i]
	res[i] = fase_all[i] - res[i]

#Residuals

p_rv = scale_period(rvx,T0,Porb)
p_all = [None]*nt
for i in range(0,nt):
	p_all[i] = scale_period(time_all[i],T0,Porb)

#error bars -> http://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html
plt.figure(2,figsize=(8,8))
plt.subplot(211)
plt.xlabel("Phase")
plt.ylabel(ylab)
#plt.ylim(1.5*min(rvy),max(rvy)*1.5)
plt.plot(p_rv,rvy,'k',label=('k=%2.2f m/s'%k ))
mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'v']
for i in range(0,nt):
	plt.errorbar(p_all[i],fase_all[i],errs_all[i],label=telescopes[i],fmt=mark[i],alpha=0.7)
	#plt.errorbar(time_all[i],fase_all[i],errs_all[i],label=telescopes[i],fmt=mark[i])
#plt.legend()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
           ncol=4, mode="expand", borderaxespad=0.)
plt.subplot(212)
plt.xlabel("Phase")
plt.ylabel(ylab)
plt.plot([0.,1.],[0.,0.],'k--')
mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'v']
for i in range(0,nt):
	plt.errorbar(p_all[i],res[i],errs_all[i],label=telescopes[i],fmt=mark[i],alpha=0.7)
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
#           ncol=4, mode="expand", borderaxespad=0.)
plt.savefig('rv_fit.png')
plt.show()
