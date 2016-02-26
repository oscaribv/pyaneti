import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm, sigmaclip
import sys
import pyaneti as pti

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
tls[1] = [3232.1,3233.1] 

plt.xlim(tls[0])
#plt.plot(hdate,wflux)
#plt.show()
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
#
#
#TAKE CARE WITH THE ARRAY SIZES
##
xt[0],yt[0], et[0] = normalize(hdate,wflux,errs,tls[0])
xt[1],yt[1], et[1] = normalize(hdate,wflux,errs,tls[1])

#plt.plot(xt[0],yt[0])
#plt.show()
#plt.plot(xt[1],yt[1])
#plt.show()

if ( ntr < 2):
	print "you do not have enought transit data!"
	sys.exit("I could crash!")

#Let us start to make a good priors calculation
T0 = [None]*ntr
for i in range(0,ntr):
	T0[i] = min(xt[i]) + 0.5*(max(xt[i])-min(xt[i]))

#
nx = [None]*ntr
ny = [None]*ntr
ne = [None]*ntr
dx = []
dy = []
de = []
for i in range(0,ntr):
	for j in range(0,len(xt[i])):
		if ( j % 1 == 0):
			dx.append(xt[i][j])
			dy.append(yt[i][j])
			de.append(et[i][j])
	nx[i] = dx
	ny[i] = dy
	ne[i] = de
	dx = []
	dy = []
	de = []

xt = []
yt = []
et = []

xt = nx
yt = ny
et = ne

#


tlimits = [None]*(ntr+1)

tlimits[0] = int(0)
for i in range(1,ntr+1):
	tlimits[i] = tlimits[i-1] + len(xt[i-1])

print tlimits

megax = np.concatenate(xt)
megay = np.concatenate(yt)
megae = np.concatenate(et)


e = 0.3
w = np.pi / 2.0
i = np.pi/2
a = 30.0
u1 = 0.42
u2 = 0.25
pz = 0.1
prec = 5e-5
maxi = int(1e8)
chi2_toler = 0.89
thin_factor = int(1e3)
ics = False

fit_e = False
fit_w = True
fit_i = True
fit_a = True
fit_u1 = False
fit_u2 = False
fit_pz = True
fit_t0 = True

what_fit = [int(fit_e),int(fit_w),int(fit_i),int(fit_a), \
					  int(fit_u1),int(fit_u2),int(fit_pz),int(fit_t0)]

print what_fit

P = T0[1] - T0[0]

print megax

#Plot the initial conditions
zvec = pti.find_z(megax,T0,e,w,P, i, a, tlimits)

print len(zvec)
print zvec

mud, mu0 = pti.occultquad(zvec,u1,u2,pz)

#plt.plot(megax,mud,'bo',megax,mu0,'rs')
#plt.show()


#sys.exit('test here')

t0o = [None]*ntr

t0o, eo, wo, io, ao, u1o, u2o, pzo = pti.metropolis_hastings_tr(megax, megay, megae, T0, e, w, i, a, u1, u2, pz,tlimits,prec, maxi, thin_factor, chi2_toler, ics, what_fit)


Po = pti.find_porb(t0o)


print t0o, Po, eo, wo, io, ao, u1o, u2o, pzo

#Plot the initial conditions
zvec = pti.find_z(megax,t0o,eo,wo,Po, io, ao, tlimits)

mud, mu0 = pti.occultquad(zvec,u1,u2,pzo)


plt.figure(1,figsize=(8,4))
plt.xlim(min(xt[0]),max(xt[0]),'r')
plt.plot(megax,megay,'bo',alpha=0.5)
plt.plot(megax,mud,'k')
plt.plot(linewidth=2.0)
plt.show()



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
