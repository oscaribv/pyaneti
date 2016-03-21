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
def planet_mass(mstar,sdm,k,sdk,P,sdP,ecc,sdecc):
  #Gravitational costant
	Gc = 6.67408e-11 #m^3 / (kgs^2)
	Gc = Gc * 1.989e30 # m^3 / (Msun s^2)
	P = P * 24. * 3600 # s

	mstar3 = mstar**(1./3.)
	const  = (2*np.pi*Gc)**( - 1./3.)
	unoe	 = np.sqrt(1.-ecc*ecc) 
	P3		 = P**(1./3.)

	mpsin = k * ( 2. * np.pi * Gc / P)**(-1./3.)  * \
	mstar**(2./3.) * unoe

	#dmdP = k/3.0 * P**(-2./3.) * const * mstar3**2 * unoe

	#dmdk = mpsin / k

	#dmdm = k * P3 * const * 2./3. * mstar(-1./3.) * unoe

	#dmde = - mpsin / unoe**2 * ecc

	return mpsin#, sdmpsin

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
  mine, med, maxe = np.percentile(xnew,[16,50,84])
	#maxe = maxe - med
	#mine = med - mine
  return med, mine, maxe



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

#This function works for more than one planet
def print_errors_planets():

	for j in range(0,nt):
		if ( is_log_rv0 ):
			vo[j] = np.power(10.,vo[j])

	for l in range(0,nplanets):

		if (is_log_P):
			Po[l] = np.power(10.,Po[l])

		if (is_ew):
			dummy_e = eo[l]
			eo[l] = eo[l] * eo[l] + wo[l] * wo[l]
			wo[l] = np.arctan2(dummy_e,wo[l])
	
		if ( fit_rv ):
			if ( is_log_k ):
				ko[l] = np.power(10.,ko[l])
	
		if ( errores == 'gauss' ):
	
			chi2_val[l], chi2_errs[l] = find_vals_gauss(chi2red[l],nconv)
			t0_val[l],t0_err[l] = find_vals_gauss(t0o[l],nconv)
			P_val[l], P_err[l]  = find_vals_gauss(Po[l],nconv)
			e_val[l],e_err[l] 	= find_vals_gauss(eo[l],nconv)
			w_val[l],w_err[l] 	= find_vals_gauss(wo[l],nconv)
			if (w_val[l] < 0.0 ):
				w_val[l] = w_val[l] + 2 * np.pi	
			w_deg[l] 		= w_val[l] * 180. / np.pi
			w_deg_err[l] = w_err[l] * 180. / np.pi
			if ( fit_rv ):
				k_val[l], k_err[l]  = find_vals_gauss(ko[l],nconv)
				v_val[l] = [None]*nt
				v_err[l] = [None]*nt
				for j in range(0,nt):
					v_val[j], v_err[j] = find_vals_gauss(vo[j],nconv)
		
			#Print the best fit values values
			print ('chi2_red = %1.4f +/- %1.4f' %(chi2_val[l],chi2_errs[l]))
			print ('The best fit planet parameters are:')
			print ('T0    = %4.4f +/- %4.4f days'%(t0_val[l],t0_err[l]))
			print ('P     = %4.4f +/- %4.4f days' 		%(P_val[l],P_err[l]))
			print ('e     = %4.4f +/- %4.4f     '			%(e_val[l],e_err[l]))
			print ('w     = %4.4f +/- %4.4f deg '	%(w_deg[l],w_deg_err[l]))
			if (fit_rv):
				print ('K    = %4.4f +/- %4.4f m/s' 		%(k_val[l]/1e-3,k_err[l]/1e-3))
				for i in range(0,nt):
					print ('%s v0 = %4.4f +/- %4.4f km/s' 	%(telescopes[i], \
					v_val[i],v_err[i]))
			
		#Percentile errors
		
		if ( errores == 'perc' ):	
	
			chi2_val[l], chi2_errl[l], chi2_errr[l] = find_vals_perc(chi2red[l],nconv)
		
			t0_val[l], t0_errl[l], t0_errr[l] = find_vals_perc(t0o[l],nconv)
			P_val[l], P_errl[l], P_errr[l]  = find_vals_perc(Po[l],nconv)
			e_val[l],e_errl[l], e_errr[l] 	= find_vals_perc(eo[l],nconv)
			w_val[l],w_errl[l], w_errr[l] 	= find_vals_perc(wo[l],nconv)
			if (w_val[l] < 0.0 ):
				w_val[l] = w_val[l] + 2 * np.pi	
			w_deg[l] 		= w_val[l] * 180. / np.pi
			w_deg_errl[l] = w_errl[l] * 180. / np.pi
			w_deg_errr[l] = w_errr[l] * 180. / np.pi
			if ( fit_rv ):
				k_val[l], k_errl[l], k_errr[l]  = find_vals_perc(ko[l],nconv)
				for j in range(0,nt):
					v_val[j], v_errl[j], v_errr[j] = find_vals_perc(vo[j],nconv)
			
		
			#Print the best fit values values
			print ('chi2_red = %1.4f + %1.4f - %1.4f' %(chi2_val[l],chi2_errr[l]-chi2_val[l], chi2_val[l] - chi2_errl[l]))
			print ('The best fit planet parameters are:')
			print ('T0    = %4.4f + %4.4f - %4.4f days'%(t0_val[l],t0_errr[l]-t0_val[l], t0_val[l]-t0_errl[l]))
			print ('P     = %4.4f + %4.4f - %4.4f days' 		%(P_val[l], P_errr[l] - P_val[l], P_val[l] - P_errl[l]))
			print ('e     = %4.4f + %4.4f - %4.4f     '			%(e_val[l], e_errr[l] - e_val[l], e_val[l] - e_errl[l]))
			print ('w     = %4.4f + %4.4f - %4.4f deg '	%(w_deg[l],w_deg_errr[l] - w_deg[l], w_deg[l] - w_deg_errl[l]))
			if (fit_rv):
				print ('K     = %4.4f + %4.4f - %4.4f m/s' 		%(k_val[l]/1.e-3,(k_errr[l] - k_val[l])/1.e-3, (k_val[l] - k_errl[l])/1e-3))
				for i in range(0,nt):
					print ('%s v0  = %4.4f + %4.4f - %4.4f km/s' 	%(telescopes[i], \
						v_val[i],v_errr[i] - v_val[i], v_val[i] - v_errl[i]))
		
