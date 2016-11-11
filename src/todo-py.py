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

  #Get the number of data and parameters
#  if (fit_rv and fit_tr ):
#    ndata = len(megax) + len(mega_rv)
###    npars = sum(what_fit) + nt - 1
#  elif(fit_rv and not fit_tr):
#    ndata = len(mega_rv)
#    npars = sum(what_fit) + nt - nplanets
#  elif(not fit_rv and fit_tr):
#    ndata = len(megax)
#    npars = sum(what_fit)
#

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
  global min_rv0, max_rv0, v0, min_k, max_k, min_phys_k, max_phys_k
  global min_P, max_P, min_phys_P, max_phys_P, min_t0, max_t0, \
         min_phys_t0, max_phys_t0, min_pz, max_pz, min_phys_pz, \
         max_phys_z, min_i, max_i, min_phys_i, max_phys_i, min_phys_rv0, max_phys_rv0

  #Let us try to do a guess for the init values
  if (fit_rv):

    #Estimate systemic velocity priors and limits from data
    #The systemic velocity value of all the telescope should
    #be between the smallest and larger RV datapoint
    #min_rv0 = min(mega_rv)
    #max_rv0 = max(mega_rv)
    #min_phys_rv0 = min(mega_rv)
    #max_phys_rv0 = max(mega_rv)

    min_rv0 = [None]*nt
    max_rv0 = [None]*nt
    min_phys_rv0 = [None]*nt
    max_phys_rv0 = [None]*nt
    for o in range(0,nt):
        min_rv0[o] = min(rv_all[o])
        max_rv0[o] = max(rv_all[o])
        min_phys_rv0[o] = min(rv_all[o])
        max_phys_rv0[o] = max(rv_all[o])


    if ( fit_v0 ):
      v0 = list(min_rv0)
    #Estimate k priors and limits from data
    if ( P.__class__ == float  ):
      if ( fit_alpha or fit_beta ):
        min_rv0 = min(min_rv0/2.,2.*min_rv0)
        max_rv0 = max(max_rv0/2.,2.*max_rv0)
        min_phys_rv0 = min(min_rv0/2.,2.*min_rv0)
        max_phys_rv0 = max(max_rv0/2.,2.*max_rv0)
      max_phys_k = 0.0
      for i in range(0,nt):
        # The maximum value ok K should be 
        max_phys_k = max_phys_k + max(rv_all[i]) - min(rv_all[i]) 
      max_phys_k = max_phys_k / nt
     
      # if we gave a worst prior for k, let us take a better
      max_k = min( [max_k,max_phys_k] )
    else:
      # The maximum value ok K should be 
    #  max_phys_k = (max_rv0 - min_rv0) 
      # if we gave a worst prior for k, let us take a better
    #  max_k = [(max_rv0 - min_rv0) ]*nplanets
    #P
      ppp = 0
    #The period should not be larger than our obsrvational run
    #max_phys_P = max(mega_time) - min(mega_time)
    #T0
    #The T0 should not be larger than our obsrvational run
    #min_phys_t0 = min(mega_time)
    #max_phys_t0 = max(mega_time)

  if (fit_tr):

    #Let us estimate limits for planet size from data 
    #The maximum depth of the data is max_flux - min_flux
    min_flux = min(megay)  
    max_flux = max(megay)  
    max_phys_pz = max_flux - min_flux   
    #Then, our maximum planet size should be (eq. 23 Winn) 
    max_phys_pz = np.sqrt(max_phys_pz)
    #Let us assume that the smallest planet in the data
    #is around 10% of the maximum depth
    min_phys_pz = 0.1 * max_phys_pz
    #If we gave a worst prior for planet size, take a better
    max_pz = min([max_pz,max_phys_pz])
    min_pz = max([min_pz,min_phys_pz])
    #P
    #tls is the list with the limits of the transits
    #The period cannot be larger (smaller)  than the outer
    # (inner) boudarie of two consecutive transit sets
    max_phys_P = tls[1][1] - tls[0][0]
    min_phys_P = tls[1][0] - tls[0][1]
    min_P = max(min_P,min_phys_P)
    max_P = min(max_P,max_phys_P)
    #t0
    #T0 should be in the first transit data set
    min_phys_t0 = tls[0][0]
    max_phys_t0 = tls[0][1]

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
  mine, med, maxe = np.percentile(x,[16,50,84])
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
def fit_new_method():
  global wtf_all, wtf_ldc, wtf_rvs
  global a_from_kepler, mstar_mean, rstar_mean, mstar_sigma_rstar_sigma
  global is_log_P, is_ew, is_b_factor, is_log_k, is_log_rv0
  global fit_t0, fit_P, fit_e, fit_w, fit_i, fit_a,fit_q1, fit_q2, fit_pz, fit_k,fit_v0
  global T0,P,e,w,ii,a,q1,q2,pz,k0,alpha,beta, v0
  global min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w, min_i, max_i, min_a,\
         max_a, min_q1, max_q1, min_q1, max_q1, min_pz, max_pz, min_k, max_k, min_alpha, max_alpha, \
         min_beta, max_beta
  global min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w, \
         min_phys_i, max_phys_i, min_phys_a, max_phys_a, min_phys_q1, max_phys_q1, min_phys_q1, \
         max_phys_q1, min_phys_pz, max_phys_pz,min_phys_k,max_phys_k, min_phys_alpha, max_phys_alpha, \
         min_phys_beta, max_phys_beta, min_phys_rv0, max_phys_rv0
  global vari,chi2,chi2red,t0o,Po,eo,wo,io,ao,q1o,q2o,pzo,ko,alphao,betao,vo, what_fit
  global new_nwalkers, good_index
  global jrvo, jtro


  if ( a_from_kepler ):
    k_log_a = False
    fit_a = False
  else:
    k_log_a = is_log_a
    fit_a = fit_a


  wtf_all = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w), \
            int(fit_i),int(fit_a), int(fit_pz), int(fit_k) ]
  wtf_rvs = [int(fit_v0)]*nt
  wtf_ldc = [fit_q1, fit_q2]

  if ( method == 'new' ):


    pars  = [T0,P,e,w,ii,a,pz,k0]
    rvs   = v0
    ldc   = [q1,q2]
    flags = [is_log_P,is_ew,is_b_factor,k_log_a,is_log_k,is_log_rv0]


    vec_rv0_limits = []
    vec_rv0_phys_limits = []
    for m in range(0,nt):
      vec_rv0_limits.append(min_rv0[m])
      vec_rv0_limits.append(max_rv0[m])
      vec_rv0_phys_limits.append(min_phys_rv0[m])
      vec_rv0_phys_limits.append(max_phys_rv0[m])

    dummy_lims = \
    [ min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w \
    , min_i, max_i, min_a, max_a, min_pz, max_pz, min_k, max_k ]

    dummy_lims_physical = \
    [min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w \
    , min_phys_i, max_phys_i, min_phys_a, max_phys_a, min_phys_pz, max_phys_pz,min_phys_k,max_phys_k ]

    limits = dummy_lims
    limits_p = dummy_lims_physical

    limits_rvs = vec_rv0_limits
    limits_p_rvs = vec_rv0_phys_limits

    limits_ldc = [ min_q1, max_q1, min_q2, max_q2]
    limits_p_ldc = [ min_phys_q1, max_phys_q1, min_phys_q2, max_phys_q2]

    pti.multi_all_stretch_move(\
                      mega_time,mega_rv,megax,megay,mega_err,megae, \
                      tlab,pars,rvs,ldc,flags,wtf_all,wtf_rvs,wtf_ldc, \
                      nwalkers,maxi,thin_factor,nconv, limits, limits_rvs, \
                      limits_ldc,limits_p, limits_p_rvs, limits_p_ldc, \
                      n_cad, t_cad, nplanets, nt)


  elif ( method == 'plot' ):
    print 'I will only print the values and generate the plot'

  else:
    print 'You did not choose a method!'
    print 'method = sm   -> Stretch move'
    print 'method = plot -> Plot of a previous run'
    sys.exit('choose your favorite.')

  print 'Reading the data file, wait a bit!'

  newfile = outdir+'/'+star+'_all_data.dat'
  if ( os.path.isfile('all_data.dat') ):
    os.rename('all_data.dat',newfile)



#-----------------------------------------------------------
#         FIT JOINT RV-TRANSIT DATA
#-----------------------------------------------------------
def fit_joint():

  global a_from_kepler, mstar_mean, rstar_mean, mstar_sigma_rstar_sigma
  global is_log_P, is_ew, is_b_factor, is_log_k, is_log_rv0
  global fit_t0, fit_P, fit_e, fit_w, fit_i, fit_a,fit_q1, fit_q2, fit_pz, fit_k,fit_v0
  global T0,P,e,w,ii,a,q1,q2,pz,k0,alpha,beta, v0
  global min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w, min_i, max_i, min_a,\
         max_a, min_q1, max_q1, min_q1, max_q1, min_pz, max_pz, min_k, max_k, min_alpha, max_alpha, \
         min_beta, max_beta
  global min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w, \
         min_phys_i, max_phys_i, min_phys_a, max_phys_a, min_phys_q1, max_phys_q1, min_phys_q1, \
         max_phys_q1, min_phys_pz, max_phys_pz,min_phys_k,max_phys_k, min_phys_alpha, max_phys_alpha, \
         min_phys_beta, max_phys_beta, min_phys_rv0, max_phys_rv0
  global vari,chi2,chi2red,t0o,Po,eo,wo,io,ao,q1o,q2o,pzo,ko,alphao,betao,vo, what_fit
  global new_nwalkers, good_index
  global jrvo, jtro


  if ( a_from_kepler ):
    k_log_a = False
    fit_a = False
  else:
    k_log_a = is_log_a
    fit_a = fit_a

  pstar = [mstar_mean,rstar_mean]
  lpstar = [mstar_sigma,rstar_sigma]

  flag = [is_log_P,is_ew,is_b_factor,k_log_a,is_log_k,is_log_rv0]

  what_fit = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w), \
              int(fit_i),int(fit_a),int(fit_q1),int(fit_q2),\
              int(fit_pz), int(fit_k),int(fit_alpha,), int(fit_beta), int(fit_v0)]

  dummy = [T0,P,e,w,ii,a,q1,q2,pz,k0,alpha,beta]
  params = np.concatenate((dummy,v0))

  #Call the fit routine

  if ( method == 'sm' ):

    vec_rv0_limits = []
    vec_rv0_phys_limits = []
    for m in range(0,nt):
      vec_rv0_limits.append(min_rv0) 
      vec_rv0_limits.append(max_rv0)
      vec_rv0_phys_limits.append(min_phys_rv0)
      vec_rv0_phys_limits.append(max_phys_rv0) 

    dummy_lims = \
    [ min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w \
    , min_i, max_i, min_a, max_a, min_q1, max_q1, min_q1, \
      max_q1, min_pz, max_pz, min_k, max_k,min_alpha, max_alpha, \
         min_beta, max_beta ]

    dummy_lims_physical = \
    [min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w \
    , min_phys_i, max_phys_i, min_phys_a, max_phys_a, min_phys_q1, max_phys_q1, min_phys_q1, \
    max_phys_q1, min_phys_pz, max_phys_pz,min_phys_k,max_phys_k, min_phys_alpha, max_phys_alpha, \
         min_phys_beta, max_phys_beta ]

    limits = np.concatenate((dummy_lims,vec_rv0_limits)) 
    limits_p = np.concatenate((dummy_lims_physical,vec_rv0_phys_limits)) 

    is_jitter = [ is_jitter_rv, is_jitter_tr ]

    pti.stretch_move(mega_time,mega_rv,mega_err,tlab \
    ,megax, megay, megae, params,pstar,lpstar,limits, limits_p , nwalkers,a_factor, maxi, thin_factor, \
    n_cad,t_cad,what_fit, flag,is_jitter,a_from_kepler, nconv,nt=nt,npars=12)

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
  dvari,dchain_lab,dchi2,djrvo,djtro,dt0o,dPo,deo,dwo,dio,dao,dq1o,dq2o,dpzo,dko, dalphao, dbetao =  \
  np.loadtxt(newfile, comments='#',unpack=True, \
  usecols=range(0,17))
  dvo = [None]*nt
  for j in range(0,nt):
    n = [17+j]
    a = np.loadtxt(newfile, comments='#', \
    unpack=True, usecols=(n))
    dvo[j] = a

  #Starting clustering
  good_index, new_nwalkers = good_clustering(dchi2,dchain_lab,nconv,nwalkers)
  vari = clustering(dvari,good_index)
  chi2 = clustering(dchi2,good_index)
  jrvo = clustering(djrvo,good_index)
  jtro = clustering(djtro,good_index)
  t0o = clustering(dt0o,good_index)
  Po = clustering(dPo,good_index)
  eo = clustering(deo,good_index)
  wo = clustering(dwo,good_index)
  io = clustering(dio,good_index)
  ao = clustering(dao,good_index)
  q1o = clustering(dq1o,good_index)
  q2o = clustering(dq2o,good_index)
  pzo = clustering(dpzo,good_index)
  ko = clustering(dko,good_index)
  alphao = clustering(dalphao,good_index)
  betao = clustering(dbetao,good_index)
  vo = [None]*nt
  for j in range(0,nt):
    vo[j] = clustering(dvo[j],good_index)


#-----------------------------------------------------------
#         FIT TRANSIT DATA
#-----------------------------------------------------------
def fit_transit():

  global a_from_kepler, mstar_mean, rstar_mean, mstar_sigma_rstar_sigma
  global is_log_P, is_ew, is_b_factor, _is_log_a
  global fit_t0, fit_P, fit_e, fit_w, fit_i, fit_a,fit_q1, fit_q2, fit_pz
  global T0,P,e,w,ii,a,q1,q2,pz
  global min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w, min_i, max_i, min_a,\
         max_a, min_q1, max_q1, min_q1, max_q1, min_pz, max_pz
  global min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w, \
         min_phys_i, max_phys_i, min_phys_a, max_phys_a, min_phys_q1, max_phys_q1, min_phys_q1, \
         max_phys_q1, min_phys_pz, max_phys_pz
  global vari,chi2,jtro,t0o,Po,eo,wo,io,ao,q1o,q2o,pzo, what_fit
  global new_nwalkers, good_index


  if ( a_from_kepler ):
    k_log_a = False
    fit_a = False
  else:
    k_log_a = is_log_a
    fit_a = fit_a
 
  pstar = [mstar_mean,rstar_mean]
  lpstar = [mstar_sigma,rstar_sigma]

  flag = [is_log_P, is_ew, is_b_factor, is_log_a]

  what_fit = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w),  \
                int(fit_i),int(fit_a), int(fit_q1),int(fit_q2),\
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
    params,pstar,lpstar,limits, limits_physical, nwalkers,a_factor,maxi, thin_factor,n_cad,t_cad, what_fit \
    ,flag,a_from_kepler,nconv)

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
  dvari,dchain_lab,dchi2,djtro,dt0o,dPo,deo,dwo,dio,dao,dq1o,dq2o,dpzo = \
  np.loadtxt(newfile, comments='#',unpack=True)

  #Starting clustering
  good_index, new_nwalkers = good_clustering(dchi2,dchain_lab,nconv,nwalkers)
  jtro = clustering(djtro,good_index)
  chi2 = clustering(dchi2,good_index)
  vari = clustering(dvari,good_index)
  t0o = clustering(dt0o,good_index)
  Po = clustering(dPo,good_index)
  eo = clustering(deo,good_index)
  wo = clustering(dwo,good_index)
  io = clustering(dio,good_index)
  ao = clustering(dao,good_index)
  q1o = clustering(dq1o,good_index)
  q2o = clustering(dq2o,good_index)
  pzo = clustering(dpzo,good_index)


#-----------------------------------------------------------
#                 FIT RV DATA
#-----------------------------------------------------------
def fit_radial_velocity():

  global a_from_kepler, mstar_mean, rstar_mean, mstar_sigma_rstar_sigma
  global is_log_P, is_ew, is_log_k, is_log_rv0
  global fit_t0, fit_P, fit_e, fit_w,fit_k,fit_v0
  global T0,P,e,w,k0, v0
  global min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w,\
         min_k, max_k, min_alpha, max_alpha, min_beta, max_beta
  global min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w, \
         min_phys_k,max_phys_k, min_phys_alpha, max_phys_alpha, min_phys_beta, max_phys_beta, min_phys_rv0, max_phys_rv0
  global vari,chi2,chi2red,t0o,Po,eo,wo,ko,alphao,betao,vo, what_fit
  global new_nwalkers, good_index
  global jrvo

  flag = [is_log_P,is_ew,is_log_k,is_log_rv0]

  if ( P.__class__ == float ):
    what_fit = [fit_t0, fit_P, fit_e, fit_w, fit_k, fit_alpha, fit_beta, fit_v0 ]
    dparams = [T0, P, e, w, k0, alpha, beta]
    params = np.concatenate((dparams,v0))
	
    vec_rv0_limits = []
    vec_rv0_phys_limits = []
    for m in range(0,nt):
      vec_rv0_limits.append(min_rv0) 
      vec_rv0_limits.append(max_rv0) 
      vec_rv0_phys_limits.append(min_phys_rv0) 
      vec_rv0_phys_limits.append(max_phys_rv0) 
	
    dummy_lims = \
    [ min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w, \
    min_k, max_k, min_alpha, max_alpha, min_beta, max_beta]

    dummy_lims_physical = \
    [ min_t0, max_t0, min_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w, \
    min_phys_k,max_phys_k,  min_alpha, max_alpha, min_beta, max_beta]

    limits = np.concatenate((dummy_lims,vec_rv0_limits)) 
    limits_p = np.concatenate((dummy_lims_physical,vec_rv0_phys_limits)) 
		
  else:
    what_fit = [None]*(8*nplanets)
    params   = [None]*((7+nt)*nplanets)	
    limits   = [None]*((7+nt)*2*nplanets)
    #Let us fill the input variables for 
    #all the number of planets
    for m in range(0,nplanets):
      #What to fit from the input lists	
      what_fit[0+8*m] = int(fit_t0[m]) 
      what_fit[1+8*m] = int(fit_P[m]) 
      what_fit[2+8*m] = int(fit_e[m]) 
      what_fit[3+8*m] = int(fit_w[m]) 
      what_fit[4+8*m] = int(fit_k[m]) 
      what_fit[5+8*m] = int(fit_alpha[m]) 
      what_fit[6+8*m] = int(fit_beta[m]) 
      what_fit[7+8*m] = int(fit_v0[m]) 
      #fill the parameters vector
      params[0+(7+nt)*m] = T0[m]
      params[1+(7+nt)*m] = P[m]
      params[2+(7+nt)*m] = e[m]
      params[3+(7+nt)*m] = w[m]
      params[4+(7+nt)*m] = k[m]
      params[5+(7+nt)*m] = alpha[m]
      params[6+(7+nt)*m] = beta[m]
      #fill the systemic velocities
      for j in range(0,nt):
        params[(7+j)+(7+nt)*m] = v0[j]
	#fill the limits
      limits[0+(7+nt)*2*m] = min_t0[m]
      limits[1+(7+nt)*2*m] = max_t0[m]
      limits[2+(7+nt)*2*m] = min_P[m]
      limits[3+(7+nt)*2*m] = max_P[m]
      limits[4+(7+nt)*2*m] = min_e[m]
      limits[5+(7+nt)*2*m] = max_e[m]
      limits[6+(7+nt)*2*m] = min_w[m]
      limits[7+(7+nt)*2*m] = max_w[m]
      limits[8+(7+nt)*2*m] = min_k[m]
      limits[9+(7+nt)*2*m] = max_k[m]
      limits[10+(7+nt)*2*m] = min_alpha[m]
      limits[11+(7+nt)*2*m] = max_alpha[m]
      limits[12+(7+nt)*2*m] = min_beta[m]
      limits[13+(7+nt)*2*m] = max_beta[m]
      for j in range(0,nt):
        limits[(14+j*2)+(7+nt)*2*m] = min_rv0
        limits[(15+j*2)+(7+nt)*2*m] = max_rv0

    vec_rv0_limits = []
    vec_rv0_phys_limits = []
    for m in range(0,nt):
      vec_rv0_limits.append(min_rv0) 
      vec_rv0_limits.append(max_rv0) 
      vec_rv0_phys_limits.append(min_phys_rv0) 
      vec_rv0_phys_limits.append(max_phys_rv0) 
	
    dummy_lims_physical = \
    [ min_phys_t0, max_phys_t0, min_phys_P, max_phys_P, min_phys_e, max_phys_e, min_phys_w, max_phys_w, \
    min_phys_k,max_phys_k,  min_phys_alpha, max_phys_alpha, min_phys_beta, max_phys_beta]

    limits_p = np.concatenate((dummy_lims_physical,vec_rv0_phys_limits)) 
	
  if ( method == 'sm' ):

    pti.stretch_move_rv(mega_time,mega_rv,mega_err,tlab,\
    params, limits,limits_p, nwalkers, a_factor, maxi, thin_factor, \
     what_fit,flag, nconv, datas=len(mega_time), nt=nt, \
    npl=nplanets, npars=7)

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
    dvari,dchain_lab,dchi2,djrvo,dt0o,dPo,deo,dwo,dko,dalphao, dbetao = \
    np.loadtxt(newfile, comments='#', unpack=True,\
    usecols=range(0,11))
    #Read the systemic velocities
    dvo = [None]*nt
    for j in range(0,nt):
      n = [11+j]
      a = np.loadtxt(newfile, comments='#',unpack=True,\
      usecols=(n))
      dvo[j] = a

    #Cluster variables
    good_index, new_nwalkers = good_clustering(dchi2,dchain_lab,nconv,nwalkers)
    vari = clustering(dvari,good_index)
    chi2 = clustering(dchi2,good_index)
    jrvo = clustering(djrvo,good_index)
    t0o = clustering(dt0o,good_index)
    Po = clustering(dPo,good_index)
    eo = clustering(deo,good_index)
    wo = clustering(dwo,good_index)
    ko = clustering(dko,good_index)
    alphao = clustering(dalphao,good_index)
    betao = clustering(dbetao,good_index)
    vo = [None]*nt
    for j in range(0,nt):
      vo[j] = clustering(dvo[j],good_index)

   
  else:
    new_nwalkers = nwalkers
    #Create all the variables, list of lists
    vari = [[]]*nplanets
    chi2 = [[]]*nplanets
    dvari = [[]]*nplanets
    t0o = [[]]*nplanets
    Po = [[]]*nplanets
    eo = [[]]*nplanets
    wo = [[]]*nplanets
    ko = [[]]*nplanets
    alphao = [[]]*nplanets
    betao = [[]]*nplanets
    #each l index is for a different planet


    for l in range(0,nplanets):
      dvari,dchain_lab,dchi2,dt0o,dPo,deo, \
      dwo,dko, dalphao, dbetao = np.loadtxt(newfile[l], comments='#', \
      unpack=True, usecols=range(0,10))
    #Cluster variables
      good_index, new_nwalkers = good_clustering(dchi2,dchain_lab,nconv,nwalkers)
      vari[l] = clustering(dvari,good_index)
      chi2[l] = clustering(dchi2,good_index)
      t0o[l] = clustering(dt0o,good_index)
      Po[l] = clustering(dPo,good_index)
      eo[l] = clustering(deo,good_index)
      wo[l] = clustering(dwo,good_index)
      ko[l] = clustering(dko,good_index)
      alphao[l] = clustering(dalphao,good_index)
      betao[l] = clustering(dbetao,good_index)

    #The  systemic velocities are the same for all the planets
    dvo = [None]*nt
    vo = [None]*nt
    for j in range(0,nt):
      n = [10+j]
      a = np.loadtxt(newfile[0], comments='#', \
      unpack=True, usecols=(n))
      dvo[j] = a

    for j in range(0,nt):
      vo[j] = clustering(dvo[j],good_index)

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
      for m in range(0,nt):
        print ('rv0= [ %4.4f , %4.4f ]' %(min_rv0[m],max_rv0[m]))
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
        for m in range(0,nt):
          print ('rv0= [ %4.4f , %4.4f ]' %(min_rv0[m],max_rv0[m]))
  print '------------------------------'
  print '     PHYSICAL LIMITS          '
  print '------------------------------'
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
    for m in range(0,nt):
      print ('rv0= [ %4.4f , %4.4f ]' %(min_rv0[m],max_rv0[m]))
  print '------------------------------'
  print '=============================='
