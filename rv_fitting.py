import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def rv_curve(t,rv0,k):
	rv = rv0 - k * np.sin( 2.0 * np.pi * ( t - T0 ) / Porb )
	return rv

#Planet parameters
Porb = 3.258907
T0 = 7064.43314 

fname='prueba.dat'

#Read the data file
data = np.loadtxt(fname,usecols=(0,1,2),comments='#')
#data is a complicated python object, let us split it
time = []
fase = []
err  = []
for i in range(0,len(data)-1):
	a,b,c=data[i][:]
	time.append(a)
	fase.append(b)
	err.append(c)

#How does the fit is made? 
#http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.curve_fit.html
popt,pcov = curve_fit(rv_curve,time,fase,sigma=err)

rvrv0= popt[0]
rvk  = popt[1]
sdrv0= np.sqrt(pcov[0][0])
sdk  = np.sqrt(pcov[1][1])

n = 1000
xmin = min(time) - Porb/2.
xmax = max(time) + Porb/2.

rvx = []
rvy = []
dn = ( xmax - xmin ) / n
rvx.append(xmin)
rvy.append(rv_curve(rvx[0],rvrv0,rvk))
for	i in range(1,n-1):
	xnew = rvx[i-1] + dn
	rvx.append( xnew )
	rvy.append( rv_curve(xnew,rvrv0,rvk) )

print ('k   = %5.5f +- %5.5f'%(rvk  , sdk  ))
print ('rv0 = %5.5f +- %5.5f'%(rvrv0, sdrv0))


#error bars -> http://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html
plt.errorbar(time, fase, yerr=err, fmt='o', color='b')
plt.plot(rvx,rvy,'r-')
plt.show()
