#-----------------------------------------------------------
#											todo-py.py
# This file contains a lot of useful of python functions.
#					Oscar Barragan, March 1st, 2016   
#-----------------------------------------------------------


# Useful libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
import sys

#-----------------------------------------------------------
#scale_period -> To calculate the periodogram of a RV curve
#input: jd -> time vector in julian date to be escaled (days)
#				T0 -> time zero of the transit (days)
#				P  -> planet orbital period (days)
#output: x -> vector with the scaled values (days)
#-----------------------------------------------------------
def scale_period(jd,Tp,P):
  x = [None]*len(jd)
  for i in range(len(jd)):
    x[i] = ( ( jd[i] - Tp ) % P ) /  P
  return x

#-----------------------------------------------------------
#planet_mass -> gives the planet mass from RV parameters
#input: mstar -> mass of the orbited star (solar masses)
#				k			-> semi-amplitude of the RV (m/s)
#				P  		-> planet orbital period (days)
#				ecc   -> orbit eccentricity (no unit)
#output:mpsin -> planet mass (stellar masses) times sin(i)
#-----------------------------------------------------------
def planet_mass(mstar,k,P,ecc):
  #Gravitational costant
  Gc = 6.67408e-11 #m^3 / (kgs^2)
  Gc = Gc * 1.989e30 # m^3 / (Msun s^2)
  P = P * 24. * 3600 # s
  mpsin = k * ( 2. * np.pi * Gc / P)**(-1./3.)  * \
  mstar**(2./3.)
  return mpsin

#-----------------------------------------------------------
#bin_data - bin a vector x each nbin points
#input: x	 		-> vector of a given len (n) 
#				nbin  -> number of points to bin
#output: nx		-> binned vector with len int(n/nbin)	 
#WARNING!
#if the len of your vector is not a multiple of nbin, you 
#will lost the last residual points
#-----------------------------------------------------------
def bin_data(x,nbin):
	nx = []
	sd = []
	dx = [None]*nbin
	for i in range(0,len(x)):
		dx[int(i%nbin)] = x[i]
		if ( (i+1)%nbin == 0 ):
			nx.append(np.mean(dx))
			sd.append(np.std(dx,ddof=1))
			dx = [None]*nbin
	
	return nx,sd
#----------------------------------------------

def find_transits(x,y):

  newy = sigmaclip(y,low=3.0,high=3.0)

  transits, miny, maxy = y - newy

#-----------------------------------------------------------
#parabola -> a parabola
#input: x	 -> vector of x values 
#				a  -> zero order term 
#				b  -> first order term 
#				c  -> second order term 
#output: y -> your parabola vector
#-----------------------------------------------------------
def parabola(x,a,b,c):
  y = a + b*x +c*x*x
  return y

#-----------------------------------------------------------
# normalize_transit -> normalizes a planet transit curve
#input: x	     -> vector with the temporal values
#				y  		 -> vector with the flux values
#				err		 -> vector with the flux error values
#				limits -> 2d vector with the start and ending date
#									of a planetary transit
#output: dummyx, dummyy, dummyerr 
#				 the normalized x, y and err vectors, respectively
#WARNING!
# This routine assumes that your transit depth is 
#	above 3 sigma of the data				 
#-----------------------------------------------------------
def normalize_transit(x,y,err,limits):
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

  popt, cov = curve_fit(parabola,newx,newy[0],sigma=newerr)
  dummyx = np.array(dummyx)
  par = parabola(dummyx,popt[0],popt[1],popt[2])
  dummyy = dummyy / par
  dummyerr = dummyerr / par

  return dummyx, dummyy, dummyerr

#-----------------------------------------------------------
# find_vals_gauss -> find the mean and standard deviation
#							  of the last nconv points of a given array
#input: x 		-> vector with a minimum size nconv
#				nconv -> the last nconv points to be taken account
#								 in the gaussian fit
#output: mu -> peak position of the gaussian fit
#				 std-> standard deviation of the gaussian fit
#-----------------------------------------------------------
def find_vals_gauss(x,nconv):
  #let us take only the converging part
  iout = len(x) - nconv
  xnew = x[iout:]
  mu,std = norm.fit(xnew)
  return mu, std

#-----------------------------------------------------------
# find_vals_perc -> find the mean and standard deviation
#							  of the last nconv points of a given array
#input: x 		-> vector with a minimum size nconv
#				nconv -> the last nconv points to be taken account
#								 in the gaussian fit
#output: mu -> peak position of the gaussian fit
#				 std-> standard deviation of the gaussian fit
#-----------------------------------------------------------
def find_vals_perc(x,nconv):
  #let us take only the converging part
  iout = len(x) - nconv
  xnew = x[iout:]
	#With a 68% confidence interval
  med, mine,maxe = np.percentile(xnew,[16,50,84])
  return med, mine, maxe


#Let us do the plots here

def plot_transit():
	  #Move all the points to T0
  for i in range(0,ntr):
    xt[i] = xt[i] - P_val * i

  #Redefine megax with the new xt values
  megax = np.concatenate(xt)
  z_val = pti.find_z(megax,t0_val,P_val,e_val,w_val\
          ,i_val,a_val)
  mud_val, mu0_val = pti.occultquad(z_val,u1_val,u2_val\
          ,pz_val)
  #Residuals
  res = megay - mud_val

  #Get the model data to do the plot
  nvec = int(1e5)
  dx = ( max(megax) - min(megax) ) / nvec
  xvec = np.zeros(nvec)
  xvec[0] = min(megax)
  for i in range(1,nvec):
    xvec[i] = xvec[i-1] + dx
  zvec = pti.find_z(xvec,t0_val,P_val,e_val,w_val,i_val,a_val)
  mud, mu0 = pti.occultquad(zvec,u1_val,u2_val,pz_val)
  #Now we have data to plot a nice model

  #Do the plot
  plt.figure(2,figsize=(10,10))
  #Plot the transit light curve
  plt.subplot(211)
  plt.xlim(min(xt[0]),max(xt[0]))
  plt.errorbar(megax,megay,megae,fmt='o',alpha=0.3)
  plt.plot(xvec,mud,'k',linewidth=2.0)
  #Plot the residuals
  plt.subplot(212)
  plt.xlim(min(xt[0]),max(xt[0]))
  plt.errorbar(megax,res,megae,fmt='o',alpha=0.3)
  plt.plot(megax,np.zeros(len(megax)),'k--',linewidth=2.0)
  plt.show()


#Plot RV
def plot_rv():
  #Create the RV fitted curve
  n = 5000
  xmin = t0_val
  xmax = t0_val + P_val
  dn = (xmax - xmin) /  n
  rvx = np.empty([n])
  rvx[0] = xmin
  for i in range(1,n):
    rvx[i] = rvx[i-1] + dn
  if ( is_circular ):
    rvy = pti.rv_circular(rvx,0.0,t0_val,k_val,P_val)
  else:
    rvy = pti.rv_curve(rvx,0.0,t0_val,k_val,P_val,e_val,w_val)

  res = [None]*nt
  for i in range(0,nt):
   if (is_circular):
      res[i] = pti.rv_circular(time_all[i],0.0,t0_val,\
               k_val,P_val)
   else:
      res[i] = pti.rv_curve(time_all[i],0.0,t0_val,k_val,\
               P_val,e_val,w_val)
   rv_all[i] = rv_all[i] - v_val[i]
   res[i] = rv_all[i] - res[i]

  p_rv = scale_period(rvx,t0_val,P_val)
  p_all = [None]*nt
  #tp_val = pti.find_tp(t0_val,e_val,w_val,P_val)
  for i in range(0,nt):
    p_all[i] = scale_period(time_all[i],t0_val,P_val)

  plt.figure(3,figsize=(10,20))
  plt.subplot(311)
  plt.xlabel("Phase")
  plt.ylabel(ylab)
  plt.ylim(-1.4*k_val,1.4*k_val)
  plt.plot(p_rv,rvy,'k',label=('k=%2.2f m/s'%k_val ))
  mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*']
  for i in range(0,nt):
    plt.errorbar(p_all[i],rv_all[i],errs_all[i],\
    label=telescopes[i],fmt=mark[i],alpha=0.6)
  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
             ncol=4, mode="expand", borderaxespad=0.)
  plt.subplot(312)
  plt.xlabel("Phase")
  plt.ylabel(ylab)
  plt.plot([0.,1.],[0.,0.],'k--')
  for i in range(0,nt):
    plt.errorbar(p_all[i],res[i],errs_all[i],\
    label=telescopes[i],fmt=mark[i],alpha=0.6)
  plt.show()


#PRINT INITIAL CONFIGURATION
def print_init():
	print ''
	print '------------------------------'
	print "INITIAL CONFIGURATION"
	print '------------------------------'
	print 'is_circular    = ', is_circular
	print 'iter max       = ', maxi
	print 'step precision = ', prec
	print 'thin factor    =', thin_factor
	print 'nconv          =',nconv
	print 'fit RV         =', fit_rv
	print 'fit Transit    =', fit_tr
	print '------------------------------'
	print 'Priors'
	print '------------------------------'
	print 'T_0         = ', T0
	print 'Period      = ', P
	print 'eccentriciy =', e
	print 'periastron  =', w
	if (fit_tr):
	  print 'sin(i)    = ', np.sin(ii)
	  print 'a/r*      = ', a
	  print 'u1        = ', u1
	  print 'u2        = ', u2
	  print 'rp/r*     = ', pz
	print '------------------------------'
	print 'What am I fitting?'
	print '------------------------------'
	print 'fit T0= ', fit_t0
	print 'fit P = ', fit_P
	print 'fit e = ', fit_e
	print 'fit w = ', fit_w
	if (fit_tr):
	  print 'fit i = ', fit_i
	  print 'fit a = ', fit_a
	  print 'fit u1= ', fit_u1
	  print 'fit u2= ', fit_u2
	  print 'fit pz= ', fit_pz
	if (fit_rv):
	  print 'fit k = ', fit_k
	  print 'fit v0= ', fit_v0
	print '------------------------------'
	print ''
	
