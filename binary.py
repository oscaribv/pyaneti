
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import lmfit as lf
import sys

#-----------------------------------
def find_anomaly(man,ecc,delta=1.e-4,imax=5000):
	anomaly = [0.0]*len(man)
	anomaly = np.array(anomaly)
	f  = anomaly - ecc * np.sin(anomaly) - man
	df =	   1.0 - ecc * np.cos(anomaly)
	counter = 0
	for i in range(0,len(man)):
		counter = 0
		while ( np.absolute(f[i]) >= delta):
			dum = anomaly[i] - f[i] / df[i]
			anomaly[i] = dum
			f[i] = anomaly[i] - ecc * np.sin(anomaly[i]) - man[i]
			df[i]=	1 - ecc * np.cos(anomaly[i])
			counter = counter + 1
			if (counter > imax):
				sys.exit("I am tired!")
	anomaly = np.sqrt(1+ecc/(1-ecc)) * np.tan(anomaly/2.0)
	anomaly = 2. * np.arctan(anomaly)
	return anomaly	
#-------------------------------------
#Linearized model from:
#http://adsabs.harvard.edu/abs/2009ApJS..182..205W
def rv_curve(t,h,c,ecc,v0,tp):
	man = 2 * np.pi * ( t - tp ) / Porb
	anomaly = find_anomaly(man,ecc)
	ut = h * np.cos(anomaly) + c * np.sin(anomaly) + v0
	return ut
#------------------------------------
def scale_period(jd,T0,P):
	x = [None]*len(jd)
	for i in range(len(jd)):
		x[i] = ( ( jd[i] - T0 ) % P ) /  P
	return x
#-----------------------------
def planet_mass(mstar,k,P,ecc):
	#Gravitational costant
	k = k*1e3
	Gc = 6.67408e-11 #m^3 / (kgs^2)
	Gc = Gc * 1.989e30 # m^3 / (Msun s^2)
	#period in seconds
	P = P * 24. * 3600
	mpsin = k * ( 2. * np.pi * Gc / P)**(-1./3.)  * \
	mstar**(2./3.) * np.sqrt( 1.0 - ecc*ecc )
	return mpsin
#-----------------------------

#----- HERE THE MAIN PROGRAM STARTS -----#

#Read the input parameters from input_rv.txt
idata = np.loadtxt('input_bin.txt',comments='#',unpack=True,dtype='S')

#Let us use all the read data
Porb = np.float(idata[0])
T0	 = np.float(idata[1])
tpes = idata[2]
fname= idata[3]
mstar= np.float(idata[4])
telescopes = [None]*len(idata[2])
for i in range(0,len(idata[2])):
	telescopes[i] = idata[2][i]
#

#Read the data file
time,fase,err,tspe = np.loadtxt(fname,usecols=(0,1,2,3), \
  dtype={'names': ('time', 'fase', 'err','telescope'), \
	'formats': ('float', 'float', 'float', 'S1')}, \
	comments='#',unpack=True)


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
		time_all.append(time_dum)
		fase_all.append(fase_dum)
		errs_all.append(errs_dum)

#The mega* variables contains all telescope data
mega_fase = []
mega_time = []
mega_err  = []
for i in range(0,nt):
	for j in range(0,len(fase_all[i])):
		mega_fase.append(fase_all[i][j])
		mega_time.append(time_all[i][j])
		mega_err.append(errs_all[i][j])


rvm = lf.Model(rv_curve)
rvm.set_param_hint('ecc',value=0.0)
rvm.set_param_hint('h',  value=1.0)
rvm.set_param_hint('c',  value=1.0)
rvm.set_param_hint('v0', value=1.0)
rvm.set_param_hint('tp', value=0.0)
result = rvm.fit(mega_fase,t=mega_time,weights=mega_err,verbose=False)

h  = result.values['h']
c  = result.values['c']
e  = result.values['ecc']
v0 = result.values['v0']

rvm.set_param_hint('ecc',value=e)
rvm.set_param_hint('h',value=h)
rvm.set_param_hint('c',value=c)
rvm.set_param_hint('v0',value=v0)
rvm.set_param_hint('tp',value=0)
result = rvm.fit(mega_fase,t=mega_time,weights=mega_err)

h  = result.values['h']
c  = result.values['c']
e  = result.values['ecc']
v0 = result.values['v0']
T0 = result.values['tp']

k = np.sqrt(h*h + c*c)
w = np.arctan(-c/h)
vp = v0 - k * e * np.cos(w)

#if ( e < 0.0 ):
#	e = - e
#	w = w + np.pi/2.0

print(result.fit_report())

#Now let us calculate the plantary mass

mpsin = planet_mass(mstar,k,Porb,e)
mjup = 0.0009543
mearth = 0.000003003 

#Print the results
print (' k = %4.4e' %(k))
print (' w = %4.4e' %(w))
print (' e = %4.4e' %(e))
print (' vp= %4.4e' %(vp))
print (' tp= %4.4e' %(T0))

print ('Planet mass of %1.4e M_j (for a %2.2f solar mass star)' % (mpsin/mjup, mstar))

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
rvy = rv_curve(rvx,h,c,e,v0,T0)

p_rv = scale_period(rvx,T0,Porb)
p_all = [None]*nt
for i in range(0,nt):
	p_all[i] = scale_period(time_all[i],T0,Porb)

#error bars -> http://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html
plt.figure(1,figsize=(10,5))
plt.xlabel("Phase")
plt.ylabel("k (m/s)")
#plt.ylim(1.5*min(rvy),max(rvy)*1.5)
plt.plot(p_rv,rvy,'k',label=('k=%2.2f m/s'%k ))
mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'v']
for i in range(0,nt):
	plt.errorbar(p_all[i],fase_all[i],errs_all[i],label=telescopes[i],fmt=mark[i])
#plt.legend()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
           ncol=4, mode="expand", borderaxespad=0.)
plt.savefig('rv_fit.png')
plt.show()

