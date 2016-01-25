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
			
def find_k_rv0(time,fase,err):
	#How is the fit made? 
	#http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.curve_fit.html
	popt,pcov = curve_fit(rv_curve,time,fase,sigma=err)

	rvrv0= popt[0]
	rvk  = popt[1]
	sdrv0= np.sqrt(pcov[0][0])
	sdk  = np.sqrt(pcov[1][1])

	print ('k   = %5.5f +- %5.5f'%(rvk  , sdk  ))
	print ('rv0 = %5.5f +- %5.5f'%(rvrv0, sdrv0))

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

#INPUT PARAMETERS
Porb = 3.258907
T0 = 7064.43314 
telescopes = ['F','H','M']
fname='prueba.dat'

#Read the data file
time,fase,err,tspe = np.loadtxt(fname,usecols=(0,1,2,3), \
  dtype={'names': ('time', 'fase', 'err','telescope'), \
	'formats': ('float', 'float', 'float', 'S1')}, \
	comments='#',unpack=True)

#Transform fase from km/s to m/s
fase=fase*1000.
err=err*1000.


#Let us separate the data for the different telescopes
time_F=[]
fase_F=[]
err_F =[]
time_H=[]
fase_H=[]
err_H =[]
time_M=[]
fase_M=[]
err_M =[]

#Let us fill this data with the corresponding telescope
for i in range(0,len(tspe)):
	if ( tspe[i] == 'F' ):
		time_F.append(time[i])
		fase_F.append(fase[i])
		err_F.append(err[i])
	if ( tspe[i] == 'H' ):
		time_H.append(time[i])
		fase_H.append(fase[i])
		err_H.append(err[i])
	if ( tspe[i] == 'M' ):
		time_M.append(time[i])
		fase_M.append(fase[i])
		err_M.append(err[i])

#Let us find the fit for each telescope
#In this way can take off the offset of each one
rv0_F = find_k_rv0(time_F,fase_F,err_F)
rv0_H = find_k_rv0(time_H,fase_H,err_H)
rv0_M = find_k_rv0(time_M,fase_M,err_M)
#Let us take off the offset for each telescope
for i in range(0,len(fase_F)):
	fase_F[i] = fase_F[i] - rv0_F

for i in range(0,len(fase_H)):
	fase_H[i] = fase_H[i] - rv0_H

for i in range(0,len(fase_M)):
	fase_M[i] = fase_M[i] - rv0_M


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

#Find all the systematic velocities and put them in rv0_all
for i in range(0,len(telescopes)):
	rv0_all[i] = find_k_rv0(time_all[i],fase_all[i],errs_all[i])

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

#Let us create a giant array with the data without offset
#mega_fase = np.concatenate((fase_F,fase_M,fase_H))
#mega_time = np.concatenate((time_F,time_M,time_H))
#mega_err  = np.concatenate((err_F,err_M,err_H))

#It is time to fit the curve for k for all the data
popt,pcov = curve_fit(rv_curve_k,mega_time,mega_fase,sigma=mega_err)

#Let us store the values of k and its corresponding sigma
k = popt[0]
sigk = np.sqrt(pcov[0][0])

#Print the result
print ('k= %8f +- %8f' %(k,sigk))

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

#Let us scale the time with the Period
p_F  = scale_period(time_F,T0,Porb)
p_H  = scale_period(time_H,T0,Porb)
p_M  = scale_period(time_M,T0,Porb)
p_rv = scale_period(rvx   ,T0,Porb)

p_all = [None]*len(telescopes)
for i in range(0,len(telescopes)):
	p_all[i] = scale_period(time_all[i],T0,Porb)

#error bars -> http://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html
plt.xlabel("Phase")
plt.ylabel("k (m/s)")
plt.ylim(1.5*min(rvy),max(rvy)*3)
plt.plot(p_rv,rvy,'k',label=('RV fit with k=%2.2f m/s'%k ))

for i in range(0,len(telescopes)):
	plt.errorbar(p_all[i],fase_all[i],errs_all[i],label=telescopes[i])
#plt.errorbar(p_F, fase_F, err_F, fmt='d', color='b', label='FIES')
#plt.errorbar(p_H, fase_H, err_H, fmt='o', color='r', label='HARPS')
#plt.errorbar(p_M, fase_M, err_M, fmt='s', color='g', label='McDonald')
plt.legend()
plt.savefig('rv_fit.png')
plt.show()


#Now let us calculate the plantary mass
ecc = 0.0

mpsin = planet_mass(1.0,k,Porb,ecc)
mjupiter = 0.0009543
mearth = 0.000003003 

print 'planet mass = '
print ('%3.5e solar masses' %mpsin)
print ('%3.5e jupiter masses' %(mpsin/mjupiter))
print ('%3.5e earth masses' %(mpsin/mearth))
