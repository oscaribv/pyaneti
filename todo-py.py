import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
import sys
import pyaneti as pti

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
  mstar**(2./3.)
  return mpsin
#-----------------------------
def bin_data(x,y,e,nbin):
  nx = []
  ny = []
  ne = []
  dx = 0.0
  dy = 0.0
  de = 0.0
  for i in range(0,len(x)):
    dx = dx + x[i]
    dy = dy + y[i]
    de = de + e[i]
    if ( (i+1)%nbin == 0 ):
      nx.append(dx/nbin)
      ny.append(dy/nbin)
      ne.append(de/nbin)
      dx = 0.0
      dy = 0.0
      de = 0.0

  return nx, ny, ne
#----------------------------------------------

def find_nsits(x,y):

  newy = sigmaclip(y,low=3.0,high=3.0)

  transits, miny, maxy = y - newy

#----------------------------------------------

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

  popt, cov = curve_fit(parabola,newx,newy[0],sigma=newerr)
  dummyx = np.array(dummyx)
  par = parabola(dummyx,popt[0],popt[1],popt[2])
  dummyy = dummyy / par
  dummyerr = dummyerr / par
  return dummyx, dummyy, dummyerr


def find_errb(x,nconv):
  #Let us take only the converging part
  iout = len(x) - nconv
  xnew = x[iout:]
  mu,std = norm.fit(xnew)
  return mu, std



