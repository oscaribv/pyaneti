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

