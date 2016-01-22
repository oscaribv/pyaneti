import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

def rv_curve(t,rv0,k):
	rv = rv0 - k * np.sin( 2.0 * np.pi * ( t - T0 ) / Porb )
	return rv

def rv_curve_k(t,k):
	rv = - k * np.sin( 2.0 * np.pi * ( t - T0 ) / Porb )
	return rv
#------------------------------------
			
def find_k_rv0(time,fase,err):
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
	for	i in range(1,n):
		xnew = rvx[i-1] + dn
		rvx.append( xnew )
		rvy.append( rv_curve(xnew,rvrv0,rvk) )

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

#Planet parameters
Porb = 3.258907
T0 = 7064.43314 

fname='prueba.dat'

#Read the data file
time,fase,err,tspe = np.loadtxt(fname,usecols=(0,1,2,3), dtype={'names': ('time', 'fase', 'err','telescope'),'formats': ('float', 'float', 'float', 'S1')}, comments='#',unpack=True)

telescopes = ['F','H','M']

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

	
rv0_F = find_k_rv0(time_F,fase_F,err_F)
rv0_H = find_k_rv0(time_H,fase_H,err_H)
rv0_M = find_k_rv0(time_M,fase_M,err_M)

th_F=[None]*len(fase_F)
for i in range(0,len(fase_F)):
	fase_F[i] = fase_F[i] - rv0_F

th_H=[None]*len(fase_H)
for i in range(0,len(fase_H)):
	fase_H[i] = fase_H[i] - rv0_H

th_M=[None]*len(fase_M)
for i in range(0,len(fase_M)):
	fase_M[i] = fase_M[i] - rv0_M

mega_fase = np.concatenate((fase_F,fase_M,fase_H))
mega_time = np.concatenate((time_F,time_M,time_H))
mega_err  = np.concatenate((err_F,err_M,err_H))

popt,pcov = curve_fit(rv_curve_k,mega_time,mega_fase,sigma=mega_err)

k = popt[0]
sigk = np.sqrt(pcov[0][0])

print ('k= %8f +- %8f' %(k,sigk))

#Let us do a nice plot

#Create the RV fitted curve
n = 100
rvx = [None]*n
rvy = [None]*n
xmin = 0
xmax = Porb
dn = (xmax - xmin) /  n 
rvx[0] = xmin 
rvy[0] = rv_curve_k(rvx[0],k)
for i in range(1,n):
	newx = rvx[i-1] + dn
	rvx[i] = newx 
	rvy[i] = rv_curve_k(rvx[i],k)

p_F = scale_period(time_F,T0,Porb)
p_H = scale_period(time_H,T0,Porb)
p_M = scale_period(time_M,T0,Porb)
p_rv= scale_period(rvx,T0,Porb)


#error bars -> http://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html
plt.plot(p_rv,rvy,'k')
plt.errorbar(p_F, fase_F, err_F, fmt='o', color='g')
plt.errorbar(p_H, fase_H, err_H, fmt='o', color='r')
plt.errorbar(p_M, fase_M, err_M, fmt='o', color='y')
plt.show()
