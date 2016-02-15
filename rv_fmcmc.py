
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
import sys
import fmcmc

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

ndata = 10000

#vmc = np.empty(ndata,nt)
kmc = np.empty(ndata)
emc = np.empty(ndata)
wmc = np.empty(ndata)
t0mc= np.empty(ndata)
pmc = np.empty(ndata)

vmc, kmc, emc, wmc, t0mc, pmc = fmcmc.mcmc_rv(mega_time,mega_fase,mega_err,tlab,v0,k0,e0,w0,T0,Porb,prec,imcmc,ndata,chi_limit)

#end fortran calling	


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
k = kval

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
	print ('For %s, v = %4.4e +/- %4.4e m/s' %(telescopes[i],vval[i],sigv[i]))

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
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
           ncol=4, mode="expand", borderaxespad=0.)
plt.subplot(212)
plt.xlabel("Phase")
plt.ylabel(ylab)
plt.plot([0.,1.],[0.,0.],'k--')
mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'v']
for i in range(0,nt):
	plt.errorbar(p_all[i],res[i],errs_all[i],label=telescopes[i],fmt=mark[i],alpha=0.7)
plt.savefig('rv_fit.png')
plt.show()
