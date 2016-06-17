#!/usr/bin/python2.7

#-----------------------------------------------------------
#                         pyaneti.py
#                        DESCRIPTION
#                   Barragan O, March 2016
#-----------------------------------------------------------

#Load libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm, sigmaclip
import sys
import os
import pyaneti as pti #FORTRAN module

#Load the input file
#You have to run the program as ./pyaneti star_name
star = str(sys.argv[1])

inf_name = 'inpy/'+star+'/input_fit.py'

if ( not os.path.isfile( inf_name ) ):
  print 'You have not created', inf_name
  sys.exit()

#Read the file with all the python functions
execfile('src/todo-py.py')

#Read the file with the default values
execfile('src/default.py')

#Read input file
execfile(inf_name)

#Prepare data
execfile('src/prepare_data.py')

#Create ouput directory
outdir = 'outpy/' + star + '_out'
if not os.path.exists(outdir):
  os.makedirs(outdir)

#PRIORS SECTION

#Let us try to do a guess for the init values
if (fit_rv):
  #Estimate systemic velocity priors and limits from data
  min_rv0 = min(mega_rv)
  max_rv0 = max(mega_rv)
  v0 = [min_rv0]*nt
  #Estimate k priors and limits from data
  if ( P.__class__ == float  ):
    max_phys_k = (max_rv0 - min_rv0) / 2.0 
    max_k = min( [max_k,max_phys_k] )
  else:
    max_phys_k = [(max_rv0 - min_rv0) / 2.0]*nplanets
    max_k = [(max_rv0 - min_rv0) / 2.0]*nplanets
  #P
  max_phys_P = max(mega_time) - min(mega_time)
  #T0
  min_phys_t0 = min(mega_time)
  max_phys_t0 = max(mega_time)
if (fit_tr):
  #Let us estimate limits for planet size from data 
  min_flux = min(megay)  
  max_flux = max(megay)  
  max_phys_pz = max_flux - min_flux 
  max_phys_pz = np.sqrt(max_phys_pz)
  min_phys_pz = max_phys_pz*0.5
  max_pz = min([max_pz,max_phys_pz])
  min_pz = max([min_pz,min_phys_pz])
  #P
  #tls is the list with the limits of the transits
  max_phys_P = tls[1][1] - tls[0][0]
  min_phys_P = tls[1][0] - tls[0][1]
  if ( min_P < min_phys_P ):
    min_P = min_phys_P
  if ( max_P > max_phys_P ):
    max_P = max_phys_P
  #t0
  min_phys_t0 = tls[0][0]
  max_phys_t0 = tls[0][1]
  #i
  if ( fit_e == False ):
    min_phys_i = ( 1. + max_phys_pz ) / min_phys_a
    min_phys_i = np.arccos(min_phys_i)
    if ( min_i < min_phys_i ):
      min_i = min_phys_i

if ( is_circular ):
  fit_e = False
  fit_w = False
  is_ew = False
  e = 1e-5
  w = np.pi / 2.0

#Print intial configuration
print_init()

#-------------------------------------------------------------
#                   FITTING ROUTINES
#-------------------------------------------------------------

#FIT TRANSIT AND RV CURVES
if (fit_rv and fit_tr ):

  if ( a_from_kepler ):
    is_log_a = False
    fit_a = False
  pstar = [mstar_mean,rstar_mean]
  lpstar = [mstar_sigma,rstar_sigma]

  flag = [is_log_P,is_ew,is_sini,is_log_a,is_log_k,is_log_rv0]

  what_fit = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w), \
              int(fit_i),int(fit_a),int(fit_q1),int(fit_q1),\
              int(fit_pz), int(fit_k), int(fit_v0)]

  dummy = [T0,P,e,w,ii,a,q1,q2,pz,k0]
  params = np.concatenate((dummy,v0))

  #Call the fit routine

  if ( method == 'sm' ):

    min_phys_t0 = min_t0
    max_phys_t0 = max_t0

    vec_rv0_limits = []
    for m in range(0,nt):
      vec_rv0_limits.append(min_rv0) 
      vec_rv0_limits.append(max_rv0) 

    dummy_lims = \
    [ min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w \
    , min_i, max_i, min_a, max_a, min_q1, max_q1, min_q1, \
      max_q1, min_pz, max_pz, min_k, max_k]

    dummy_lims_physical = \
    [min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w \
    , min_phys_i, max_phys_i, min_phys_a, max_phys_a, min_phys_q1, max_phys_q1, min_phys_q1, \
    max_phys_q1, min_phys_pz, max_phys_pz,min_phys_k,max_phys_k]

    limits = np.concatenate((dummy_lims,vec_rv0_limits)) 
    limits_p = np.concatenate((dummy_lims_physical,vec_rv0_limits)) 


    pti.stretch_move(mega_time,mega_rv,mega_err,tlab \
    ,megax, megay, megae, params,pstar,lpstar,limits, limits_p , nwalkers,a_factor, maxi, thin_factor, \
    n_cad,t_cad,what_fit, flag,a_from_kepler, nconv)

  elif ( method == 'plot' ):
    print 'I will only print the values and generate the plot'

  else:
    print 'You did not choose a method!'
    print 'method = sm   -> Stretch move'
    print 'method = plot -> Plot of a previous run'
    sys.exit('choose your favorite.')

  print 'Reading the data file, wait a bit!'

  newfile = outdir+'/'+star+'_rv-tr.dat'
  if ( os.path.isfile('mh_fit.dat') ):
    os.rename('mh_fit.dat',newfile)
        
  #Read the data
  vari,chi2,chi2red,t0o,Po,eo,wo,io,ao,q1o,q2o,pzo,ko =  \
  np.loadtxt(newfile, comments='#',unpack=True, \
  usecols=range(0,13))
  vo = [None]*nt
  for j in range(0,nt):
    n = [13+j]
    a = np.loadtxt(newfile, comments='#', \
    unpack=True, usecols=(n))
    vo[j] = a

#-----------------------------------------------------------
#                   FIT TRANSIT CURVE ONLY
#-----------------------------------------------------------

elif ( not fit_rv and fit_tr ):

  min_phys_t0 = min_t0
  max_phys_t0 = max_t0

  flag = [is_log_P, is_ew, is_sini, is_log_a]

  what_fit = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w),  \
                int(fit_i),int(fit_a), int(fit_q1),int(fit_q1),\
                int(fit_pz)]

  params = [T0,P,e,w,ii,a,q1,q2,pz]

  if ( method == 'sm' ):
    #The transit time should be in the first window
    limits = \
    [ min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w \
    , min_i, max_i, min_a, max_a, min_q1, max_q1, \
    min_q1, max_q1, min_pz, max_pz]
               
    limits_physical = \
    [min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w \
    , min_phys_i, max_phys_i, min_phys_a, max_phys_a, min_phys_q1, max_phys_q1, min_phys_q1, \
    max_phys_q1, min_phys_pz, max_phys_pz]

    pti.stretch_move_tr(megax, megay, megae,  \
    params,limits, limits_physical, nwalkers,a_factor,maxi, thin_factor,n_cad,t_cad, what_fit,flag,nconv)

  elif ( method == 'plot' ):
    print 'I will only print the values and generate the plot'

  else:
    print 'You did not choose a method!'
    print 'method = sm   -> Stretch move'
    print 'method = plot -> Plot of a previous run'
    sys.exit('choose your favorite.')

  print 'Reading the data file, wait a bit!'

  newfile = outdir+'/'+star+'_tr.dat'
  if ( os.path.isfile('mh_trfit.dat') ):
    os.rename('mh_trfit.dat',newfile)
  #Read the data
  vari, chi2,chi2red,t0o,Po,eo,wo,io,ao,q1o,q2o,pzo = \
  np.loadtxt(newfile, comments='#',unpack=True)

#-------------------------------------------------------------
#                   FIT RV CURVE ONLY
#-------------------------------------------------------------

elif ( fit_rv and not fit_tr ):

  flag = [is_log_P,is_ew,is_log_k,is_log_rv0]

  if ( P.__class__ == float ):
    what_fit = [fit_t0, fit_P, fit_e, fit_w, fit_k, fit_v0 ]
    dparams = [T0, P, e, w, k0]
    params = np.concatenate((dparams,v0))
	
    vec_rv0_limits = []
    for m in range(0,nt):
      vec_rv0_limits.append(min_rv0) 
      vec_rv0_limits.append(max_rv0) 
	
    dummy_lims = \
    [ min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w, \
    min_k, max_k]

    dummy_lims_physical = \
    [ min_t0, max_t0, min_P, 1000, 1e-10, 0.99, 0.0, 2*np.pi, \
    1e-3,1e4]

    limits = np.concatenate((dummy_lims,vec_rv0_limits)) 
    limits_p = np.concatenate((dummy_lims_physical,vec_rv0_limits)) 
		
  else:
    what_fit = [None]*6*nplanets
    params   = [None]*(5+nt)*nplanets	
    limits   = [None]*(5+nt)*2*nplanets
    #Let us fill the input variables for 
    #all the number of planets
    for m in range(0,nplanets):
      #What to fit from the input lists	
      what_fit[0+6*m] = int(fit_t0[m]) 
      what_fit[1+6*m] = int(fit_P[m]) 
      what_fit[2+6*m] = int(fit_e[m]) 
      what_fit[3+6*m] = int(fit_w[m]) 
      what_fit[4+6*m] = int(fit_k[m]) 
      what_fit[5+6*m] = int(fit_v0[m]) 
      #fill the parameters vector
      params[0+(5+nt)*m] = T0[m]
      params[1+(5+nt)*m] = P[m]
      params[2+(5+nt)*m] = e[m]
      params[3+(5+nt)*m] = w[m]
      params[4+(5+nt)*m] = k[m]
      #fill the systemic velocities
      for j in range(0,nt):
        params[(5+j)+(5+nt)*m] = v0[j]
	#fill the limits
      limits[0+(5+nt)*2*m] = min_t0[m]
      limits[1+(5+nt)*2*m] = max_t0[m]
      limits[2+(5+nt)*2*m] = min_P[m]
      limits[3+(5+nt)*2*m] = max_P[m]
      limits[4+(5+nt)*2*m] = min_e[m]
      limits[5+(5+nt)*2*m] = max_e[m]
      limits[6+(5+nt)*2*m] = min_w[m]
      limits[7+(5+nt)*2*m] = max_w[m]
      limits[8+(5+nt)*2*m] = min_k[m]
      limits[9+(5+nt)*2*m] = max_k[m]
      for j in range(0,nt):
        limits[(10+j*2)+(5+nt)*2*m] = min_rv0
        limits[(11+j*2)+(5+nt)*2*m] = max_rv0

  if ( method == 'sm' ):

    pti.stretch_move_rv(mega_time,mega_rv,mega_err,tlab,\
    params, limits, nwalkers, a_factor, maxi, thin_factor, \
    what_fit,flag,nconv,datas=len(mega_time), \
    nt=nt,npl=nplanets)

  elif ( method == 'plot' ):
    print 'I will only print the values and generate the plot'

  else:
    print 'You did not choose a method!'
    print 'method = sm   -> Stretch move'
    print 'method = plot -> Plot of a previous run'
    sys.exit('choose your favorite.')

  print 'Reading the data file, wait a bit!'

  if ( nplanets == 1):
    out_file = 'planet1.dat'
    newfile = outdir+'/'+star+'_rv.dat'
    if ( os.path.isfile(out_file) ):
      os.rename(out_file,newfile)
  elif ( nplanets > 1):
    out_file = [None]*nplanets
    newfile = [None]*nplanets
    for m in range(0,nplanets):
      out_file[m] = 'planet' + str(m+1) + '.dat'
      newfile[m] = outdir+'/'+star+'_rv'+str(m+1)+'.dat'
      if ( os.path.isfile(out_file[m]) ):
        os.rename(out_file[m],newfile[m])

  if ( nplanets == 1 ):
    vari,chi2,chi2red,t0o,Po,eo,wo,ko = \
    np.loadtxt(newfile, comments='#', unpack=True,\
    usecols=range(0,8))
    #Read the systemic velocities
    vo = [None]*nt
    for j in range(0,nt):
      n = [8+j]
      a = np.loadtxt(newfile, comments='#',unpack=True,\
      usecols=(n))
      vo[j] = a
  else:
    #Create all the variables, list of lists
    vari = [[]]*nplanets
    chi2 = [[]]*nplanets
    chi2red = [[]]*nplanets
    t0o = [[]]*nplanets
    Po = [[]]*nplanets
    eo = [[]]*nplanets
    wo = [[]]*nplanets
    ko = [[]]*nplanets
    #each l index is for a different planet
    for l in range(0,nplanets):
      vari[l],chi2[l],chi2red[l],t0o[l],Po[l],eo[l], \
      wo[l],ko[l] = np.loadtxt(newfile[l], comments='#', \
      unpack=True, usecols=range(0,8))
    #The  systemic velocities are the same for all the planets
    vo = [None]*nt
    for j in range(0,nt):
      n = [8+j]
      a = np.loadtxt(newfile[0], comments='#', \
      unpack=True, usecols=(n))
      vo[j] = a
		
#Nothing to fit!
else:
  sys.exit("Nothing to fit!")

#-------------------------------------------------------------
#             	END FITTING ROUTINES
#-------------------------------------------------------------

#Print the values
execfile('src/print_values.py')

#Create plots
execfile('src/plot_data.py')

if ( nplanets == 1):

  plot_histogram()

  plot_correlations()

  #PLOT TRANSIT
  if ( fit_tr ):
    plot_transit()

  #PLOT RV CURVE
  if ( fit_rv ):
    plot_rv_one()
	
elif (nplanets > 1):

  hist_mp_rv()

  #plot_correlations()

  #PLOT THE RV curves
  plot_rv_mp()

