#More information about curve_fit in the following link
#http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.curve_fit.html

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

#-----------------------------------
#This function find the eccentry anomaly by using the
#Newton-Rapshon method
def find_anomaly(man,ecc,delta=1.e-5):
	anomaly = [0.0]*len(man)
	f  = anomaly - ecc * np.sin(anomaly) - man
	df =	   1.0 - ecc * np.cos(anomaly)
	#print f,df
	#sys.exit()
	for i in range(0,len(man)):
		while ( np.absolute(f[i]) >= delta):
			dum = anomaly[i] - f[i] / df[i]
			anomaly[i] = dum
			f[i] = anomaly[i] - ecc * np.sin(anomaly[i]) - man[i]
			df[i]=	1 - ecc * np.cos(anomaly[i])
	return anomaly	
#-------------------------------------
	#This functions assumes that the orbit is circular
def rv_circular(t,rv0,k):
	rv = rv0 + k * np.sin( 2.0 * np.pi * ( t - T0 ) / Porb )
	#rv = rv0 + k * np.sin( 2.0 * np.pi * ( t ) / Porb )
	return rv
#------------------------------------
def rv_curve_rv0(t,rv0,k0,ecc,omega):
	#The general equation for rv motion is
	man = 2 * np.pi * (t - T0) / Porb
	#man = 2 * np.pi * (t) / Porb
	anomaly = find_anomaly(man,ecc)
	rv = rv0 + k0 * (np.cos( anomaly + omega ) + ecc*np.cos( omega ) ) / np.sqrt(1.0 - ecc*ecc)
	return rv
#------------------------------------
def rv_curve(t,k0,ecc,omega):
	#The general equation for rv motion is
	man = 2 * np.pi * ( t - T0 ) / Porb
	#man = 2 * np.pi * ( t ) / Porb
	anomaly = find_anomaly(man,ecc)
	rv = + k0 * (np.cos( anomaly + omega ) + ecc*np.cos( omega ) )  / np.sqrt(1.0 - ecc*ecc)
	return rv
#------------------------------------
#Extract systemic velocities is highly important, so let us take care with this
def find_rv0(time,fase,err,tpe):
	#Let us first fit assuming a circular orbit (to first estimate v0 and k0)
	popt,pcov = curve_fit(rv_circular,time,fase)
	#This results will be our first guesses for v0 and k0
	rv0 = popt[0]
	k0  = popt[1]
	#Now let us fit again now with your guesses as input
	#popt,pcov = curve_fit(rv_circular,time,fase,sigma=err,p0=[rv0,k0])
	#rv0 = popt[0]
	#k0  = popt[1]
	#Now rv0 and kv0 are really good guesses to try to fit an eccentric orbit
	popt,pcov = curve_fit(rv_curve_rv0,time,fase,sigma=err,p0=[rv0,k0,0,0])
	rv0 = popt[0]
	k0  = popt[1]
	w   = popt[3]
	#Now we have guesses for rv0, k0 and w, let us re do the fit, now to improve e
	popt,pcov = curve_fit(rv_curve_rv0,time,fase,sigma=err,p0=[rv0,k0,0,w])
	rvrv0= popt[0]
	sdrv0= np.sqrt(pcov[0][0])
	print ('for %1s -> rv0 = %5.5f +/- %5.5f m/s'%(tpe,rvrv0,sdrv0))
	return rvrv0
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

#Read the input parameters from input_rv.txt
idata = np.loadtxt('input_rv.txt',comments='#',unpack=True,dtype='S')

#Let us use all the read data
Porb = np.float(idata[0])
T0	 = np.float(idata[1])
tpes = idata[2]
fname= idata[3]
telescopes = [None]*len(idata[2])
for i in range(0,len(idata[2])):
	telescopes[i] = idata[2][i]
#

#Read the data file
time,fase,err,tspe = np.loadtxt(fname,usecols=(0,1,2,3), \
  dtype={'names': ('time', 'fase', 'err','telescope'), \
	'formats': ('float', 'float', 'float', 'S1')}, \
	comments='#',unpack=True)

#Transform fase from km/s to m/s
fase=fase*1000.
err=err*1000.

#These lists have lists with data for the different telescopes
time_all=[]
fase_all=[]
errs_all=[]

#tspe is the data with the telescopes read from the file
for i in range(0,len(telescopes)):
	time_dum =[]
	fase_dum =[]
	errs_dum =[]
	for j in range(0,len(tspe)):
		if (tspe[j] == telescopes[i]):
			time_dum.append(time[j])
			fase_dum.append(fase[j])
			errs_dum.append(err[j])
	time_all.append(time_dum)
	fase_all.append(fase_dum)
	errs_all.append(errs_dum)

rv0_all=[None]*len(telescopes)

print ("Extracting systemic velocities for the %i telescopes"%len(telescopes))
#Find all the systematic velocities and put them in rv0_all
for i in range(0,len(telescopes)):
	rv0_all[i] = find_rv0(time_all[i],fase_all[i],errs_all[i],telescopes[i])

#Let us take off the offset for each telescope
for i in range(0,len(telescopes)):
	for j in range(0,len(fase_all[i])):
		fase_all[i][j] = fase_all[i][j] - rv0_all[i]

mega_fase = []
mega_time = []
mega_err  = []
#Let us create a mega array with all the data
for i in range(0,len(telescopes)):
	for j in range(0,len(fase_all[i])):
		mega_fase.append(fase_all[i][j])
		mega_time.append(time_all[i][j])
		mega_err.append(errs_all[i][j])

#It is time to fit the curve for k for all the data
popt,pcov = curve_fit(rv_curve,mega_time,mega_fase,sigma=mega_err,p0=[1,0,0])
kg = popt[0]
wg = popt[2]
popt,pcov = curve_fit(rv_curve,mega_time,mega_fase,sigma=mega_err,p0=[kg,0,wg])

#Let us store the values of k and its corresponding sigma
k = popt[0]
sigk = np.sqrt(pcov[0][0])
e = popt[1]
sige = np.sqrt(pcov[1][1])
w = popt[2]
sigw = np.sqrt(pcov[2][2])


#Now let us calculate the plantary mass

mstar = 1.0
mpsin = planet_mass(mstar,k,Porb,e)
mjup = 0.0009543
mearth = 0.000003003 

#Print the results
print ('semi-amplitude k= %8f +- %8f m/s' %(k,sigk))
print ('semi-amplitude e= %8f +- %8f m/s' %(e,sige))
print ('semi-amplitude w= %8f +- %8f m/s' %(w,sigw))

print ('This corresponds to a planet mass of %1.4e M_j for a %2.2f solar mass star' % (mpsin/mjup, mstar))

#Let us do a nice plot

#Create the RV fitted curve
n = 500
xmin = T0
xmax = T0 + Porb
dn = (xmax - xmin) /  n  
rvx = np.empty([n])
rvx[0] = xmin 
for i in range(1,n):
	rvx[i] = rvx[i-1] + dn
rvy = rv_curve(rvx,k,e,w)

p_rv = scale_period(rvx,T0,Porb)
p_all = [None]*len(telescopes)
for i in range(0,len(telescopes)):
	p_all[i] = scale_period(time_all[i],T0,Porb)

#error bars -> http://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html
plt.xlabel("Phase")
plt.ylabel("k (m/s)")
plt.ylim(1.5*min(rvy),max(rvy)*3)
plt.plot(p_rv,rvy,'k',label=('RV fit with k=%2.2f m/s'%k ))
mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'v']
for i in range(0,len(telescopes)):
	plt.errorbar(p_all[i],fase_all[i],errs_all[i],label=telescopes[i],fmt=mark[i])
plt.legend()
plt.savefig('rv_fit.png')
plt.show()

