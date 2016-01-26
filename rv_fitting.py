import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

#-----------------------------------

def rv_curve(t,rv0,k):
	rv = rv0 - k * np.sin( 2.0 * np.pi * ( t - T0 ) / Porb )
	return rv

#------------------------------------

def rv_curve_k(t,k):
	rv = - k * np.sin( 2.0 * np.pi * ( t - T0 ) / Porb )
	return rv

#------------------------------------
			
def find_k_rv0(time,fase,err,tpe):
	#How is the fit made? 
	#http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.curve_fit.html
	popt,pcov = curve_fit(rv_curve,time,fase,sigma=err)

	rvrv0= popt[0]
	rvk  = popt[1]
	sdrv0= np.sqrt(pcov[0][0])
	sdk  = np.sqrt(pcov[1][1])

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
	rv0_all[i] = find_k_rv0(time_all[i],fase_all[i],errs_all[i],telescopes[i])

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
popt,pcov = curve_fit(rv_curve_k,mega_time,mega_fase,sigma=mega_err)

#Let us store the values of k and its corresponding sigma
k = popt[0]
sigk = np.sqrt(pcov[0][0])

#Now let us calculate the plantary mass
ecc = 0.0

mstar = 1.0
mpsin = planet_mass(mstar,k,Porb,ecc)
mjup = 0.0009543
mearth = 0.000003003 

#Print the results
print ('semi-amplitude k= %8f +- %8f m/s' %(k,sigk))

print ('This corresponds to a planet mass of %1.4e M_j for a %2.2f solar mass star' % (mpsin/mjup, mstar))

#Let us do a nice plot

#Create the RV fitted curve
n = 100
rvx = [None]*n
rvy = [None]*n
xmin = T0
xmax = T0 + Porb
dn = (xmax - xmin) /  n  
rvx[0] = xmin 
rvy[0] = rv_curve_k(rvx[0],k)
for i in range(1,n):
	newx   = rvx[i-1] + dn
	rvx[i] = newx 
	rvy[i] = rv_curve_k(rvx[i],k)

p_rv = scale_period(rvx   ,T0,Porb)
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

