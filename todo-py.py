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
def scale_period(jd,T0,P):
  x = [None]*len(jd)
  for i in range(len(jd)):
    x[i] = ( ( jd[i] - T0 ) % P ) /  P
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
			sd.append(np.std(dx))
			dx = [None]*nbin
	
	sd = sd / np.sqrt(nbin)

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
# find_vals_gauss -> Find the mean and standard deviation
#							  of the last nconv points of a given array
#input: x 		-> vector with a minimum size nconv
#				nconv -> the last nconv points to be taken account
#								 in the gaussian fit
#output: mu -> peak position of the gaussian fit
#				 std-> standard deviation of the gaussian fit
#-----------------------------------------------------------
def find_vals_gauss(x,nconv):
  #Let us take only the converging part
  iout = len(x) - nconv
  xnew = x[iout:]
  mu,std = norm.fit(xnew)
  return mu, std

