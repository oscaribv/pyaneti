#-----------------------------------------------------------
#                    todo-py.py
# This file contains a lot of useful of python functions.
#	    Oscar Barragan, March, 2016   
#-----------------------------------------------------------


# Useful libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
import sys

#Calculated as of http://arxiv.org/abs/1501.05415
def get_BIC(chi2tot_val):

  #Get the number of data and parameters
  if (fit_rv and fit_tr ):
    ndata = len(megax) + len(mega_rv)
    npars = sum(what_fit) + nt - 1
  elif(fit_rv and not fit_tr):
    ndata = len(mega_rv)
    npars = sum(what_fit) + nt - nplanets
  elif(not fit_rv and fit_tr):
    ndata = len(megax)
    npars = sum(what_fit)

#  if ( fit_rv and not fit_tr ):
#    rv_dum = []
#    for j in range(0,nt):
#      rv_dum.append(rv_all[j])
#    res = [None]*nt
#    for j in range(0,nt):
#	  #This is the model of the actual planet
#	    res[j] = pti.rv_curve_mp(time_all[j],0.0,t0_val,k_val,\
#	    P_val,e_val,w_val)
#	    #the actual value, minus the systemic velocity
#	    rv_dum[j] = rv_dum[j] - v_val[j] 
#	    res[j] = rv_dum[j] - res[j]

  #variance = np.concatenate(res)
  #variance = sum(abs(variance*variance)) / ndata

  #BIC = ndata * np.log(variance) + npars * np.log(ndata)

  BIC = chi2tot_val + npars * np.log(ndata)  

  return BIC


#-----------------------------------------------------------
#scale_period -> To calculate the periodogram of a RV curve
#input: jd -> time vector in julian date to be escaled (days)
# T0 -> time zero of the transit (days)
# P  -> planet orbital period (days)
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
# k    -> semi-amplitude of the RV (m/s)
# P    -> planet orbital period (days)
# ecc  -> orbit eccentricity (no unit)
# i    -> orbit inclination (radians)
#output: mp -> planet mass (solar masses)
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
  P = P * 24. * 3600. # s

  unoe = np.sqrt(1.-ecc*ecc) 

  mpsin = k * ( 2. * np.pi * Gc / P)**(-1./3.)  * \
  mstar**(2./3.) * unoe

  mp = mpsin / np.sin(i)

  #find the mass by solving the mass function, this is useful for 
  #stars orbited by other stars

  f = [1.0]*len(P)
  cte = - unoe**3 * P * k**3 / 2. / np.pi / Gc 
  sini = np.sin(i)
  flag = True
  while ( flag ):
    f = cte + (mp * sini )**3 / ( mstar + mp )**2
    #df= 3. * mp**2 * sini**3 / ( mstar + mp )**2 - ( 2. * mp**3 * sini**3 / ( mstar + mp )**3 )
    df= mp**2 * sini**3 / ( mstar + mp )**2 * ( 3. - 2. * mp / (mstar + mp ) )
    mp = mp - f/df
    for j in range(0,len(P)):
      if ( f[j] > 1.e-8 ):
        flag = True
        break
      else:
        flag = False

  return mp

#-----------------------------------------------------------

def get_rhostar(P,a):

  P = P * 24. * 3600. # s
  rho = 3. * np.pi * a**3 / ( G_cgs * P * P)

  return rho

#-----------------------------------------------------------

#Get ranges, assumes that the consecutive 
#data is almost equally spaced
def get_transit_ranges(times):

  xlimits = []
  xlimits.append(times[0])
  dx = times[3] - times[0]
  for i in range(1,len(times)-1):
    if ( times[i+1] - times[i] > dx ): # We are in other transit
      xlimits.append(times[i])
      xlimits.append(times[i+1])
  xlimits.append(times[len(times)-1]) 

  print xlimits
  
  lims = [None]*(len(xlimits)/2)
  for i in range(0,len(xlimits),2):
    lims[i/2] = [xlimits[i],xlimits[i+1]]

  ntr = len(lims)

  return lims, ntr

#-----------------------------------------------------------
#bin_data - bin a vector x each nbin points
# input: x -> vector of a given len (n) 
# nbin  -> number of points to bin
#output: nx -> binned vector with len int(n/nbin)	 
#WARNING!
#if the len of your vector is not a multiple of nbin, you 
#will lost the last residual points
#-----------------------------------------------------------
def bin_data_mean(x,nbin):
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
#-----------------------------------------------------------
def bin_data_median(x,nbin):
  nx = []
  sd = []
  dx = [None]*nbin
  for i in range(0,len(x)):
    dx[int(i%nbin)] = x[i]
    if ( (i+1)%nbin == 0 ):
      nx.append(np.median(dx))
      sd.append(np.std(dx,ddof=1)/np.sqrt(nbin-1))
      dx = [None]*nbin
  
  return nx, sd
#----------------------------------------------

def find_transits(x,y):

  newy = sigmaclip(y,low=3.0,high=3.0)

  transits, miny, maxy = y - newy

#-----------------------------------------------------------
#parabola -> a parabola
#input: x -> vector of x values 
#       a  -> zero order term 
#       b  -> first order term 
#       c  -> second order term 
#output: y -> your parabola vector
#-----------------------------------------------------------
def parabola(x,a,b,c):
  y = a + b*x +c*x*x
  y =  y - y + 1.0
  return y

#-----------------------------------------------------------
# normalize_transit -> normalizes a planet transit curve
#input: x -> vector with the temporal values
#       y -> vector with the flux values
#     err -> vector with the flux error values
#  limits -> 2d vector with the start and ending date
#	     of a planetary transit
#output: dummyx, dummyy, dummyerr 
#     the normalized x, y and err vectors, respectively
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

def separate_transits(x,y,err,limits):
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

  return dummyx, dummyy, dummyerr


#-----------------------------------------------------------
# find_vals_gauss -> find the mean and standard deviation
#               of the last nconv points of a given array
#input: x -> vector with a minimum size nconv
#   nconv -> the last nconv points to be taken account
#	     in the gaussian fit
#output: mu -> peak position of the gaussian fit
#        std-> standard deviation of the gaussian fit
#-----------------------------------------------------------
def find_vals_gauss(x,nconv):
  #let us take only the converging part
  iout = len(x) - nconv
  xnew = x[iout:]
  mu,std = norm.fit(xnew)
  return mu, std

#-----------------------------------------------------------
# find_vals_perc -> find the median and the errors within
#  a 68% confidence interval
#input: x -> vector with a minimum size nconv
#   nconv -> the last nconv points to be taken account
#            in the gaussian fit
#output: med -> median value
#	mine -> left error (50% - 16%)
#	maxe -> right error (84% - 50%)
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
  print 'iter max       = ', maxi
  print 'thin factor    = ', thin_factor
  print 'nconv          = ', nconv
  print 'nwalkers       = ', nwalkers
  print '------------------------------'
  print 'fit RV         = ', fit_rv
  print 'fit Transit    = ', fit_tr
  print '------------------------------'
  if (fit_tr):
    print 'LC data        = ', lc_data
    print ('cadence time   =  %2.3f min'%(t_cad*60.*24))
    print 'n rebinning    = ', n_cad
  print '------------------------------'
  print 'fitting T0     = ', fit_t0
  print 'fitting P      = ', fit_P
  print 'fitting e      = ', fit_e
  print 'fitting w      = ', fit_w
  if (fit_tr):
    print 'fitting i      = ', fit_i
    print 'fitting a      = ', fit_a
    print 'fitting q1     = ', fit_q1
    print 'fitting q2     = ', fit_q2
    print 'fitting pz     = ', fit_pz
  if (fit_rv):
    print 'fitting k      = ', fit_k
    print 'fitting v0     = ', fit_v0
  print '------------------------------'
  print '        PRIOR RANGES          '
  print '------------------------------'
  if ( min_t0.__class__ == float ):
    print ('T0 = [ %4.4f , %4.4f ]' %(min_t0,max_t0))
    print ('P  = [ %4.4f , %4.4f ]' %(min_P,max_P))
    print ('e  = [ %4.4f , %4.4f ]' %(min_e,max_e))
    print ('w  = [ %4.4f , %4.4f ]' %(min_w,max_w))
    if (fit_tr):
      print ('i  = [ %4.4f , %4.4f ]' %(min_i,max_i))
      print ('a  = [ %4.4f , %4.4f ]' %(min_a,max_a))
      print ('pz = [ %4.4f , %4.4f ]' %(min_pz,max_pz))
      print ('q1 = [ %4.4f , %4.4f ]' %(min_q1,max_q1))
      print ('q2 = [ %4.4f , %4.4f ]' %(min_q2,max_q2))
    if (fit_rv):
      print ('K  = [ %4.4f , %4.4f ]' %(min_k,max_k))
      print ('rv0= [ %4.4f , %4.4f ]' %(min_rv0,max_rv0))
  else:
    for j in range(0,nplanets):
      print 'Planet ', j + 1
      print ('T0 = [ %4.4f , %4.4f ]' %(min_t0[j],max_t0[j]))
      print ('P  = [ %4.4f , %4.4f ]' %(min_P[j],max_P[j]))
      print ('e  = [ %4.4f , %4.4f ]' %(min_e[j],max_e[j]))
      print ('w  = [ %4.4f , %4.4f ]' %(min_w[j],max_w[j]))
      if (fit_tr):
        print ('i  = [ %4.4f , %4.4f ]' %(min_i[j],max_i[j]))
        print ('a  = [ %4.4f , %4.4f ]' %(min_a[j],max_a[j]))
        print ('pz = [ %4.4f , %4.4f ]' %(min_pz[j],max_pz[j]))
        print ('q1 = [ %4.4f , %4.4f ]' %(min_q1[j],max_q1[j]))
        print ('q2 = [ %4.4f , %4.4f ]' %(min_q2[j],max_q2[j]))
      if (fit_rv):
        print ('K  = [ %4.4f , %4.4f ]' %(min_k[j],max_k[j]))
        print ('rv0= [ %4.4f , %4.4f ]' %(min_rv0,max_rv0))
  print '------------------------------'
  print '     PHYSICAL LIMITS          '
  print '------------------------------'
  if ( min_t0.__class__ == float ):
    print ('T0 = [ %4.4f , %4.4f ]' %(min_phys_t0,max_phys_t0))
    print ('P  = [ %4.4f , %4.4f ]' %(min_phys_P,max_phys_P))
    print ('e  = [ %4.4f , %4.4f ]' %(min_phys_e,max_phys_e))
    print ('w  = [ %4.4f , %4.4f ]' %(min_phys_w,max_phys_w))
    if (fit_tr):
      print ('i  = [ %4.4f , %4.4f ]' %(min_phys_i,max_phys_i))
      print ('a  = [ %4.4f , %4.4f ]' %(min_phys_a,max_phys_a))
      print ('pz = [ %4.4f , %4.4f ]' %(min_phys_pz,max_phys_pz))
      print ('q1 = [ %4.4f , %4.4f ]' %(min_phys_q1,max_phys_q1))
      print ('q2 = [ %4.4f , %4.4f ]' %(min_phys_q2,max_phys_q2))
    if (fit_rv):
      print ('K  = [ %4.4f , %4.4f ]' %(min_phys_k,max_phys_k))
      print ('rv0= [ %4.4f , %4.4f ]' %(min_rv0,max_rv0))
  print '------------------------------'
  print '=============================='
