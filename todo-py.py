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
#       i     -> orbit inclination (radians)
#output:mp    -> planet mass (solar masses)
# WARNING: the defalut value for i is pi/2 (90 degrees),
# if you do not know the orbit inclination, this function
# computes the planet mass times sin(i)
#-----------------------------------------------------------
def planet_mass(mstar,k,P,ecc,i=np.pi/2.):

  #Gravitational costant
	#Gc = 6.67408e-11 #m^3 / (kgs^2)
	Gc = 6.67259e-11
	#Gc = Gc * 1.989e30 # m^3 / (Msun s^2)
	Gc = Gc * 1.98855e30 # m^3 / (Msun s^2)
	P = P * 24. * 3600 # s

	unoe	 = np.sqrt(1.-ecc*ecc) 

	mpsin = k * ( 2. * np.pi * Gc / P)**(-1./3.)  * \
	mstar**(2./3.) * unoe

	mp = mpsin / np.sin(i)

	return mp

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
			sd.append(np.std(dx,ddof=1)/np.sqrt(nbin-1))
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
# find_vals_perc -> find the median and the errors within
#									 a 68% confidence interval
#input: x 		-> vector with a minimum size nconv
#				nconv -> the last nconv points to be taken account
#								 in the gaussian fit
#output: med 	-> median value
#				 mine	-> left error (50% - 16%)
#				 maxe	-> right error (84% - 50%)
#-----------------------------------------------------------
def find_vals_perc(x,nconv,sf=1.0):
	iout = len(x) - nconv
	xnew = x[iout:]
	#With a 68% confidence interval
	mine, med, maxe = np.percentile(xnew,[16,50,84])
	maxe = ( maxe - med ) / sf
	mine = ( med - mine ) / sf
	return med, mine, maxe


#-----------------------------------------------------------
# PRINT INITIAL CONFIGURATION
#-----------------------------------------------------------
def print_init():
	print ''
	print '=============================='
	print '------------------------------'
	print "    INITIAL CONFIGURATION     "
	print '------------------------------'
	print 'is_circular    = ', is_circular
	print 'iter max       = ', maxi
	print 'step precision = ', prec
	print 'thin factor    =', thin_factor
	print 'nconv          =',nconv
	print 'fit RV         =', fit_rv
	print 'fit Transit    =', fit_tr
	print '------------------------------'
	print '          Priors              '
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
	print '     What am I fitting?       '
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
	print '=============================='
	print ''


