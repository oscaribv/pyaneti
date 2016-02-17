import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm, sigmaclip
import sys
import frv

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

ntr = 2

tls = [None]*ntr 

tls[0] = [3217,3218]
#tls[1] = 

plt.xlim(tls[0])
plt.plot(hdate,wflux)
plt.show()
#Let us try to fit a parabola to the curve

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

	#print len(x),len(newx), len(newy[0])	
	popt, cov = curve_fit(parabola,newx,newy[0],sigma=newerr)
	dummyx = np.array(dummyx)
	par = parabola(dummyx,popt[0],popt[1],popt[2])
	#plt.plot(newx,newy[0],'b',dummyx,par,'k')
	#plt.show()
	#print popt

	dummyy = dummyy / par
	dummyerr = dummyerr / par
	#plt.plot(dummyx,dummyy,'k')
	#plt.show()

	return dummyx, dummyy, dummyerr
	
	#Now we have  the values around the transit
	
xt= [None]*ntr	
yt= [None]*ntr	
et= [None]*ntr	
xt[0],yt[0], et[0] = normalize(hdate,wflux,errs,tls[0])

plt.plot(xt[0],yt[0],'k')
plt.show()

def mandel_agol():
	#Write here the function
	#This function has to be written in fortran	

sys.exit('I am here!')

plt.figure(2,figsize=(8,8))
plt.subplot(211)
plt.xlabel("Phase")
plt.ylabel(ylab)
plt.plot(p_rv,rvy,'k',label=('k=%2.2f m/s'%k ))
mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'v']
for i in range(0,nt):
	plt.errorbar(p_all[i],rv_all[i],errs_all[i],label=telescopes[i],fmt=mark[i],alpha=0.6)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
           ncol=4, mode="expand", borderaxespad=0.)
plt.subplot(212)
plt.xlabel("Phase")
plt.ylabel(ylab)
plt.plot([0.,1.],[0.,0.],'k--')
mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'v']
for i in range(0,nt):
	plt.errorbar(p_all[i],res[i],errs_all[i],label=telescopes[i],fmt=mark[i],alpha=0.6)
plt.savefig('rv_fit.png')
plt.show()
