#-----------------------------------------------------------
#                    todo-py.py
# This file contains a lot of useful of python functions.
#	    Oscar Barragan, March, 2016   
#-----------------------------------------------------------

# Useful libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

#-----------------------------------------------------------
#This suborutine calculates the Bayesian Information Criteria
#Calculated as of http://arxiv.org/abs/1501.05415
# Input  -> chi2 (not reduced)
# output -> BIC
#-----------------------------------------------------------
def get_BIC(chi2tot_val):

  npars = sum(wtf_all) + sum(wtf_ldc) + sum(wtf_rvs)
  ndata = len(megax) + len(mega_rv)

  BIC = chi2tot_val + npars * np.log(ndata)

  return BIC

#-----------------------------------------------------------
#This function returns the phase of a temporal array given
#the period
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

  P = P * 24. * 3600. # s

  #Let us make a guess by assuming mstar >> mp
  unoe = np.sqrt(1.-ecc*ecc) 
  mpsin = k * ( 2. * np.pi * S_GM_SI / P)**(-1./3.)  * \
  mstar**(2./3.) * unoe
  mp = mpsin / np.sin(i)

  #find the mass by solving the mass function, 
  #this is useful for stars orbited by other stars

  f = [1.0]*len(P)
  cte = - unoe**3 * P * k**3 / 2. / np.pi / S_GM_SI 
  sini = np.sin(i)
  flag = True
  #Start Newton-Raphson algorithm
  while ( flag ):
    f = cte + (mp * sini )**3 / ( mstar + mp )**2
    df= mp**2 * sini**3 / ( mstar + mp )**2 * ( 3. - 2. * mp / (mstar + mp ) )
    mp = mp - f/df
    for j in range(0,len(P)):
      #check that all the array elemets have converged
      if ( f[j] > 1.e-8 ):
        flag = True
        break
      else:
        flag = False

  return mp

#-----------------------------------------------------------
# This routine calculates the stellar density
# Based on eq. 30 from Winn., 2014
# Assuming the companion is too small
# Input:
# P -> Period
# a -> semi-major axis
# Output:
# rho -> stellar density
#-----------------------------------------------------------
def get_rhostar(P,a):
  P = P * 24. * 3600. # s
  rho = 3. * np.pi * a**3 / ( G_cgs * P * P)
  return rho


#-----------------------------------------------------------

def get_teq(Tstar,albedo,rstar,a):
  Tp = Tstar*( 1.0 - albedo)**(0.25)
  Tp = (rstar/2.0/a)**(0.5) * Tp
  return Tp

#-----------------------------------------------------------
# check if we force a circular orbit
#-----------------------------------------------------------
def check_circular():
 
  if (nplanets == 1 ):
    if ( is_circular ):
      fit_e = False
      fit_w = False
      is_ew = False
      e = 1e-5
      w = np.pi / 2.0
  else:
      fit_e = [False]*nplanets
      fit_w = [False]*nplanets
      is_ew = False
      e = [1e-5]*nplanets
      w = [np.pi / 2.0]*nplanets

#-----------------------------------------------------------
#  Smart priors, get the best values of the physical and 
#  priors limits 
#-----------------------------------------------------------
def smart_priors():
  #We are using global variables
  global fit_tr, fit_rv
  global tota_rv_fit, total_tr_fit
  global min_rv0, max_rv0, v0, min_k, max_k, min_phys_k, max_phys_k
  global min_P, max_P, min_phys_P, max_phys_P, min_t0, max_t0, \
         min_phys_t0, max_phys_t0, min_pz, max_pz, min_phys_pz, \
         max_phys_pz, min_i, max_i, min_phys_i, max_phys_i, min_phys_rv0, max_phys_rv0

  #Let us try to do a guess for the init values
  if ( total_rv_fit ):

    #Estimate systemic velocity priors and limits from data
    #The systemic velocity value of all the telescope should
    #be between the smallest and larger RV datapoint
    min_rv0 = [None]*nt
    max_rv0 = [None]*nt
    min_phys_rv0 = [None]*nt
    max_phys_rv0 = [None]*nt
    for o in range(0,nt):
        min_rv0[o] = min(rv_all[o]) - 0.1
        max_rv0[o] = max(rv_all[o]) + 0.1
        min_phys_rv0[o] = min(rv_all[o]) - 1.
        max_phys_rv0[o] = max(rv_all[o]) + 1.

    v0 = list(min_rv0)

  if ( total_tr_fit ):

    min_flux = min(megay)
    max_flux = max(megay)
    min_phys_pz =[None]*len(xt)
    max_phys_pz =[None]*len(xt)
    for o in range(0,nplanets):
      min_phys_pz[o] = 0.0
      max_phys_pz[o] = max_flux - min_flux
      max_phys_pz[o] = 2.*np.sqrt(max_phys_pz[o])
      #Let us assume that the smallest planet in the data
      #is around 10% of the maximum depth
      #If we gave a worst prior for planet size, take a better
      max_pz[o] = min([max_pz[o],max_phys_pz[o]])

    min_phys_P =[None]*len(xt)
    max_phys_P =[None]*len(xt)
    min_phys_t0=[None]*len(xt)
    max_phys_t0=[None]*len(xt)
    for o in range(0,nplanets):
      #Max and minimum Period for each planet from the data
      #----          ----|min_phys_P|----         ------
      #    -       -                     -       -
      #     -------                       -------
      #|                  max_phys_P                    |
      #|min_phys_t0
      #     max_phys_t0  |
      max_phys_P[o] = xt[o][1][len(xt[o][1])-1] - xt[o][0][0]
      min_phys_P[o] = xt[o][1][0] - xt[o][0][len(xt[o][0])-1]
      min_P[o] = max(min_P[o],min_phys_P[o])
      max_P[o] = min(max_P[o],max_phys_P[o])
      min_phys_t0[o] = xt[o][0][0]
      max_phys_t0[o] = xt[o][0][len(xt[o][0])-1]

    #If we are ussing spectroscopic priors, let us estimante the limits on a
    ms_min = mstar_mean - 5.*mstar_sigma
    ms_max = mstar_mean + 5.*mstar_sigma
    rs_min = rstar_mean - 5.*rstar_sigma
    rs_max = rstar_mean + 5.*rstar_sigma
    for o in range(0,nplanets):
      if ( a_from_kepler[o] ):
        p_min = min_P[o]
        p_max = max_P[o]
        min_phys_a[o] = pti.get_a_scaled(ms_min,rs_max,p_min)
        max_phys_a[o] = pti.get_a_scaled(ms_max,rs_min,p_max)
        print min_phys_a[o], max_phys_a[o]

  #sys.exit()


#-----------------------------------------------------------
#Get transit ranges, assumes that the consecutive
#data is almost equally spaced
# it assumes the times vector has a gap between transit sets
#   ----          ---       (GAP)      ----          ----
#       -        -                         -        -
#        --------                           --------
#Input:
# times -> vector with the time stamps
#Output:
# lims  -> transit limits list size (ntr), each list has
#          two values, the first and last transit datapoint
# ntr   -> number of transits
#-----------------------------------------------------------
def get_transit_ranges(times,gap):

  xlimits = []
  xlimits.append(times[0])
  dx = times[gap] - times[0] #the between transit gap
  for i in range(1,len(times)-1):
    if ( times[i+1] - times[i] > dx ):
      # We are in other transit
      xlimits.append(times[i])
      xlimits.append(times[i+1])
  xlimits.append(times[len(times)-1])
  #xlims has all the data points when the transit data sets
  #start and end

  #Let us put them in lims
  lims = [None]*(len(xlimits)/2)
  for i in range(0,len(xlimits),2):
    lims[i/2] = [xlimits[i],xlimits[i+1]]

  ntr = len(lims)
  print 'Number of transits = ', ntr
  print 'Transit limits ='
  print lims

  return lims, ntr

#-----------------------------------------------------------
# This subroutine separate the transit data in different
# data sets, each data set has a transit signal
# Input:
#   x      -> time array
#   y      -> flux array
#   err    -> error array
#   limits -> time limits for a given transit
# Output:
#   dummyx   ->  time vector with the time data between
#                limits[0] and limits[1]
#   dummyy   ->  flux vector with the time data between 
#                limits[0] and limits[1]
#   dummyerr ->  errors vector with the time data between
#                limits[0] and limits[1]
#-----------------------------------------------------------
def separate_transits(x,y,err,limits):
  dummyx  = []
  dummyy  = []
  dummyerr= []
  for i in range(0,len(x)):
    if ( x[i] > limits[0] and x[i] < limits[1] ):
      dummyx.append(x[i])
      dummyy.append(y[i])
      dummyerr.append(err[i])

  return dummyx, dummyy, dummyerr

#-----------------------------------------------------------
# find_vals_perc -> find the median and the errors within
#  a 68% credible interval
#input: 
#       x -> vector with a minimum size nconv
#   nconv -> the last nconv points to be taken account
#            in the gaussian fit
#output:
#        med -> median value
#	mine -> left error (50% - 16%)
#	maxe -> right error (84% - 50%)
#-----------------------------------------------------------
def find_vals_perc(x,sf=1.0):
  #With a 68% confidence interval
  mine, med, maxe = np.percentile(x,[16.0,50.0,84.0])
  maxe = ( maxe - med ) / sf
  mine = ( med - mine ) / sf
  
  return med, mine, maxe


#-----------------------------------------------------------

def good_clustering(chi2,chain_lab,nconv,nwalkers):
  #Let us find the good indixes for the cluster
  #We have n walkers

  print 'STARTING CHAIN CLUSTERING'
  print 'Initial number of chains:', nwalkers

  #Extract all the chains
  chi2_walkers = [None]*nwalkers
  chi2_mean = [None]*nwalkers
  walk_dummy = []
  for i in range(0,nwalkers):
    for j in range (0,len(chain_lab)):
      if (chain_lab[j] == i ):
        walk_dummy.append(chi2[j])
    chi2_walkers[i] = walk_dummy
    walk_dummy = []


  #The mean of each walker
  for i in range (0,nwalkers):
    chi2_mean[i] = np.mean(chi2_walkers[i])

  #get the minimum chi2
  total_min = min(chi2_mean)

  good_chain = []
  #Let us kill all the walkers 5 times the minimum
  for i in range(0,nwalkers):
    if ( chi2_mean[i]/total_min < 1.0 + clustering_delta ):
      #We are saving the good chain labels
      good_chain.append(i)

  #Now we know how many good chains we have
  new_nwalkers = len(good_chain)

  print 'Final number of chains:', new_nwalkers

  #Let us save the good index
  good_index = []
  for i in range(0, len(chain_lab)):
      for j in range(0,len(good_chain)):
          if ( chain_lab[i] == good_chain[j] ):
              good_index.append(i)

  return good_index, new_nwalkers

#-----------------------------------------------------------

def clustering(par,good_index):

  cluster_par = np.zeros(len(good_index))
  for i in range(0,len(good_index)):
    n = good_index[i]
    cluster_par[i] = par[n]

  return cluster_par

#-----------------------------------------------------------
#         FIT JOINT RV-TRANSIT DATA
#-----------------------------------------------------------
def joint_fit():
  global wtf_all, wtf_ldc, wtf_rvs, nt
  global a_from_kepler, mstar_mean, rstar_mean, mstar_sigma_rstar_sigma
  global is_log_P, is_ew, is_b_factor, is_log_k, is_log_rv0
  global fit_t0, fit_P, fit_e, fit_w, fit_i, fit_a,fit_q1, fit_q2, fit_pz, fit_k,fit_v0
  global T0,P,e,w,ii,a,q1,q2,pz,k0,alpha,beta, v0
  global min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w, min_i, max_i, min_a,\
         max_a, min_q1, max_q1, min_q1, max_q1, min_pz, max_pz, min_k, max_k, min_alpha, max_alpha, \
         min_beta, max_beta, min_rv0, max_rv0
  global min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w, \
         min_phys_i, max_phys_i, min_phys_a, max_phys_a, min_phys_q1, max_phys_q1, min_phys_q1, \
         max_phys_q1, min_phys_pz, max_phys_pz,min_phys_k,max_phys_k, min_phys_alpha, max_phys_alpha, \
         min_phys_beta, max_phys_beta, min_phys_rv0, max_phys_rv0
  global vari,chi2,chi2red,t0o,Po,eo,wo,io,ao,q1o,q2o,pzo,ko,alphao,betao,vo, what_fit
  global new_nwalkers, good_index
  global jrvo, jtro

  wtf_all = [None]*8*nplanets
  for o in range(0,nplanets):
    wtf_all[o*8:(o+1)*8] = [int(fit_t0[o]),int(fit_P[o]),int(fit_e[o]),int(fit_w[o]), \
                            int(fit_i[o]),int(fit_a[o]), int(fit_pz[o]), int(fit_k[o]) ]

  wtf_rvs = []
  for o in range(0,nt):
    wtf_rvs.append(int(fit_v0))

  wtf_ldc = [int(fit_q1), int(fit_q2)]

  #Let us check what do we want to fit
  total_fit_flag = [ total_rv_fit, total_tr_fit ]
  pars = [None]*8*nplanets
  for i in range(0,nplanets):
    pars[i*8:(i+1)*8]= [T0[i],P[i],e[i],w[i],ii[i],a[i],pz[i],k0[i]]

  rvs = v0
  ldc   = [q1,q2]
  flags = [is_log_P,is_ew,is_b_factor,is_log_a,is_log_k,is_log_rv0]

  if ( method == 'mcmc' ):

    vec_rv0_limits = []
    vec_rv0_phys_limits = []
    for m in range(0,nt):
      vec_rv0_limits.append(min_rv0[m])
      vec_rv0_limits.append(max_rv0[m])
      vec_rv0_phys_limits.append(min_phys_rv0[m])
      vec_rv0_phys_limits.append(max_phys_rv0[m])

    dummy_lims = [None]*8*2*nplanets
    dummy_lims_physical = [None]*8*2*nplanets
    for o in range(0,nplanets):
      dummy_lims[o*8*2:(o+1)*8*2 ] = \
      [ min_t0[o], max_t0[o], min_P[o], max_P[o], min_e[o], max_e[o], min_w[o], max_w[o] \
      , min_i[o], max_i[o], min_a[o], max_a[o], min_pz[o], max_pz[o], min_k[o], max_k[o] ]
      dummy_lims_physical[o*8*2:(o+1)*8*2] = \
      [min_phys_t0[o], max_phys_t0[o], min_phys_P[o], max_phys_P[o], min_phys_e[o], max_phys_e[o], min_phys_w[o], max_phys_w[o] \
      , min_phys_i[o], max_phys_i[o], min_phys_a[o], max_phys_a[o], min_phys_pz[o], max_phys_pz[o],min_phys_k[o],max_phys_k[o] ]

    limits = dummy_lims
    limits_p = dummy_lims_physical

    limits_rvs = vec_rv0_limits
    limits_p_rvs = vec_rv0_phys_limits

    limits_ldc = [ min_q1, max_q1, min_q2, max_q2]
    limits_p_ldc = [ min_phys_q1, max_phys_q1, min_phys_q2, max_phys_q2]

    stellar_pars = [mstar_mean,mstar_sigma,rstar_mean,rstar_sigma]
    is_jitter = [is_jitter_rv, is_jitter_tr]

    pti.multi_all_stretch_move(\
    mega_time,mega_rv,megax,megay,mega_err,megae, \
    tlab,megap,pars,rvs,ldc,stellar_pars,a_from_kepler,\
    flags,total_fit_flag,is_jitter,wtf_all,wtf_rvs,wtf_ldc, \
    nwalkers,maxi,thin_factor,nconv, limits, limits_rvs, \
    limits_ldc,limits_p, limits_p_rvs, limits_p_ldc, \
    n_cad, t_cad, nplanets, nt)


  elif ( method == 'plot' ):
    print 'I will only print the values and generate the plot'

  else:
    print 'You did not choose a method!'
    print 'method = mcmc   -> Run the MCMC code'
    print 'method = plot   -> Plot of a previous run'
    sys.exit('choose your favorite.')

  newfile = outdir+'/'+star+'_all_data.dat'
  if ( os.path.isfile('all_data.dat') ):
    os.rename('all_data.dat',newfile)

  newfile_jitter = outdir+'/'+star+'_jitter_data.dat'
  if ( os.path.isfile('jitter_data.dat') ):
    os.rename('jitter_data.dat',newfile_jitter)

#-----------------------------------------------------------
#          PRINT INITIAL CONFIGURATION
#-----------------------------------------------------------
def print_init():
  out_init_file = outdir+'/'+star+'_init.dat'
  oif = open(out_init_file,'w')
  oif.write ('\n')
  oif.write ('==============================\n')
  oif.write ('------------------------------\n')
  oif.write ("    INITIAL CONFIGURATION     \n")
  oif.write ('------------------------------\n')
  oif.write ('Star           = %s\n'%star)
  oif.write ('No. planets    = %d\n'%nplanets)
  oif.write ('------------------------------\n')
  oif.write ('iter max       = %d\n' %maxi)
  oif.write ('thin factor    = %d\n' %thin_factor)
  oif.write ('nconv          = %d\n' %nconv)
  oif.write ('nwalkers       = %d\n' %nwalkers)
  oif.write ('------------------------------\n')
  oif.write ('fit RV         = %s\n' %fit_rv)
  oif.write ('fit Transit    = %s\n' %fit_tr)
  oif.write ('------------------------------\n')
  if ( total_tr_fit ):
    oif.write ('LC data        = %s\n' %lc_data)
    oif.write ('cadence time   =  %2.3f min\n'%(t_cad*60.*24))
    oif.write ('n rebinning    = %d\n' %n_cad)
    oif.write ('Stellar priors = %s\n' %a_from_kepler)
  oif.write ('------------------------------\n')
  oif.write ('fitting T0     = %s\n'% fit_t0)
  oif.write ('fitting P      = %s\n'% fit_P)
  oif.write ('fitting e      = %s\n'% fit_e)
  oif.write ('fitting w      = %s\n'% fit_w)
  if (fit_tr):
    oif.write ('fitting i      = %s\n'% fit_i)
    oif.write ('fitting a      = %s\n'% fit_a)
    oif.write ('fitting q1     = %s\n'% fit_q1)
    oif.write ('fitting q2     = %s\n'% fit_q2)
    oif.write ('fitting pz     = %s\n'% fit_pz)
  if (fit_rv):
    oif.write ('fitting k      = %s\n'% fit_k)
    oif.write ('fitting v0     = %s\n'% fit_v0)
  for j in range(0,nplanets):
    oif.write ('------------------------------\n')
    oif.write ('  PLANET %s \n' %(star + plabels[j]))
    oif.write ('------------------------------\n')
    oif.write ('        PRIOR RANGES          \n')
    oif.write ('------------------------------\n')
    oif.write ('T0 = [ %4.4f , %4.4f ]\n' %(min_t0[j],max_t0[j]))
    oif.write ('P  = [ %4.4f , %4.4f ]\n' %(min_P[j],max_P[j]))
    oif.write ('e  = [ %4.4f , %4.4f ]\n' %(min_e[j],max_e[j]))
    oif.write ('w  = [ %4.4f , %4.4f ]\n' %(min_w[j],max_w[j]))
    oif.write ('i  = [ %4.4f , %4.4f ]\n' %(min_i[j],max_i[j]))
    oif.write ('a  = [ %4.4f , %4.4f ]\n' %(min_a[j],max_a[j]))
    oif.write ('pz = [ %4.4f , %4.4f ]\n' %(min_pz[j],max_pz[j]))
    oif.write ('K  = [ %4.4f , %4.4f ]\n' %(min_k[j],max_k[j]))
    oif.write ('------------------------------\n')
    oif.write ('     PHYSICAL LIMITS          \n')
    oif.write ('------------------------------\n')
    oif.write ('T0 = [ %4.4f , %4.4f ]\n' %(min_phys_t0[j],max_phys_t0[j]))
    oif.write ('P  = [ %4.4f , %4.4f ]\n' %(min_phys_P[j],max_phys_P[j]))
    oif.write ('e  = [ %4.4f , %4.4f ]\n' %(min_phys_e[j],max_phys_e[j]))
    oif.write ('w  = [ %4.4f , %4.4f ]\n' %(min_phys_w[j],max_phys_w[j]))
    oif.write ('i  = [ %4.4f , %4.4f ]\n' %(min_phys_i[j],max_phys_i[j]))
    oif.write ('a  = [ %4.4f , %4.4f ]\n' %(min_phys_a[j],max_phys_a[j]))
    oif.write ('pz = [ %4.4f , %4.4f ]\n' %(min_phys_pz[j],max_phys_pz[j]))
    oif.write ('K  = [ %4.4f , %4.4f ]\n' %(min_phys_k[j],max_phys_k[j]))
  oif.write ('------------------------------\n')
  oif.write ('   Other parameters limits \n')
  oif.write ('------------------------------\n')
  oif.write ('q1 = [ %4.4f , %4.4f ]\n' %(min_q1,max_q1))
  oif.write ('q2 = [ %4.4f , %4.4f ]\n' %(min_q2,max_q2))
  for m in range(0,nt):
    oif.write ('%s = [ %4.4f , %4.4f ]\n' %(telescopes_labels[m],min_rv0[m],max_rv0[m]))
  oif.write ('------------------------------\n')
  oif.write ('Other parameters physical limits\n ')
  oif.write ('------------------------------\n')
  oif.write ('q1 = [ %4.4f , %4.4f ]\n' %(min_phys_q1,max_phys_q1))
  oif.write ('q2 = [ %4.4f , %4.4f ]\n' %(min_phys_q2,max_phys_q2))
  for m in range(0,nt):
    oif.write ('%s = [ %4.4f , %4.4f ]\n' %(telescopes_labels[m],min_phys_rv0[m],max_phys_rv0[m]))
  oif.write ('==============================\n')

  oif.close()
  dummy_file = open(out_init_file)
  for line in dummy_file:
    print line,
  dummy_file.close()
