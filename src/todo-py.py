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
#Return the equilibrium temeprature given the stellar temperature
#albedo, stellar radius and distance to the star
def get_teq(Tstar,albedo,rstar,a):
  Tp = Tstar*( 1.0 - albedo)**(0.25)
  Tp = (rstar/2.0/a)**(0.5) * Tp
  return Tp

#Sigma clipping functions copied from exotrending
#x and y are the original arrays, z is the vector with the residuals
def sigma_clip(x,y,z,limit_sigma=5,is_plot=False):
  control = True
  new_y = list(y)
  new_x = list(x)
  new_z = list(z)
  dummy_x = []
  dummy_y = []
  dummy_z = []
  n = 1
  while ( control ):
    sigma = np.std(new_z)
    for i in range(0,len(new_z)):
      if ( np.abs(new_z[i]) < limit_sigma*sigma ):
        dummy_x.append(new_x[i])
        dummy_y.append(new_y[i])
        dummy_z.append(new_z[i])
    if ( len(dummy_x) == len(new_x) ): #We did not cut, so the sigma clipping is done
      control = False
    new_y = list(dummy_y)
    new_x = list(dummy_x)
    new_z = list(dummy_z)
    dummy_x = []
    dummy_y = []
    dummy_z = []
    n = n + 1

  if ( is_plot ):
    plt.plot(x,y,'or',new_x,new_y,'ob')
    plt.show()

  return new_x, new_y

#-----------------------------------------------------------
#  Smart priors, get the best values of the physical and
#  priors limits
#-----------------------------------------------------------
def smart_priors():
  #We are using global variables
  global fit_tr, fit_rv
  global tota_rv_fit, total_tr_fit
  global min_rv0, max_rv0, v0, min_k, max_k
  global min_P, max_P, min_t0, max_t0, \
         min_rp, max_rp, \
         min_i, max_i

  #Let us try to do a guess for the init values
  if ( total_rv_fit ):

    #Estimate systemic velocity priors and limits from data
    #The systemic velocity value of all the telescope should
    #be between the smallest and larger RV datapoint
    min_rv0 = [None]*nt
    max_rv0 = [None]*nt
    for o in range(0,nt):
        if (fit_v0 == 'u'):
          min_rv0[o] = min(rv_all[o]) - 1.0e-1
          max_rv0[o] = max(rv_all[o]) + 1.0e-1
          if is_linear_trend == 'u':
            min_rv0[o] = min(rv_all[o]) - 1.0
            max_rv0[o] = max(rv_all[o]) + 1.0
        else:
          min_rv0[o] = 0.0
          max_rv0[o] = 0.0

#Create function to create xt vectors
def create_transit_data(time,flux,errs,planet=0,span=0.0):
  global trt_vec, P_vec, T0_vec #vector with the transit durations calculated in print_values.py

  P  = best_value(P_vec[planet],maxloglike,get_value)
  T0 = best_value(T0_vec[planet],maxloglike,get_value)
  tt = best_value(trt_vec[planet],maxloglike,get_value)
  tt = tt/24.0

  if ( span < 1e-5 ):
    #We have to calculate things
    span = 3*tt
#  else:
    #We have to use the span given by the user

  #Let us fold the data to the period
  t_inicial = T0 - span/2

  folded_t = list(time)
  for o in range(0,len(time)):
    folded_t[o] = int( ( time[o] - t_inicial ) / P )
    folded_t[o] = time[o] - folded_t[o] * P
    folded_t[o] = folded_t[o] - T0

  lt = []
  xt = []
  yt = []
  et = []

  for o in range(0,len(time)):
      if ( folded_t[o] > - span/2 and folded_t[o] < span/2 ):
          lt.append(time[o])
          xt.append(folded_t[o])
          yt.append(flux[o])
          et.append(errs[o])

  #let us separate the transits
  lt_out = []
  xt_out = []
  yt_out = []
  et_out = []
  lt_d = []
  xt_d = []
  yt_d = []
  et_d = []
  lt_d.append(lt[0])
  xt_d.append(xt[0])
  yt_d.append(yt[0])
  et_d.append(et[0])
  for o in range(1,len(xt)):
      if ( xt[o] > xt[o-1] ):
        lt_d.append(lt[o])
        xt_d.append(xt[o])
        yt_d.append(yt[o])
        et_d.append(et[o])
      else:
        lt_out.append(lt_d)
        xt_out.append(xt_d)
        yt_out.append(yt_d)
        et_out.append(et_d)
        lt_d = []
        xt_d = []
        yt_d = []
        et_d = []
        lt_d.append(lt[o])
        xt_d.append(xt[o])
        yt_d.append(yt[o])
        et_d.append(et[o])


  lt_out.append(lt_d)
  xt_out.append(xt_d)
  yt_out.append(yt_d)
  et_out.append(et_d)

  return lt_out, xt_out, yt_out, et_out

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
def find_vals_perc(x,sf=1.0,prob=68.3):
  #With a 68% confidence interval
  mnval = 50.0 - prob/2.0
  mxval = 50.0 + prob/2.0
  mine, med, maxe = np.percentile(x,[mnval,50.0,mxval])
  maxe = ( maxe - med ) / sf
  mine = ( med - mine ) / sf

  return med, mine, maxe


#-----------------------------------------------------------
def best_value(vector,loglike,cual):
    if ( cual == 'median'):
        result = np.median(vector)
    elif( cual == 'mode' ):
        result = my_mode(vector)
    elif( cual == 'maxlike' ):
        maxindex = np.argmax(loglike)
        result = vector[maxindex]

    return result


#-----------------------------------------------------------
#This routine calculates the mode of a vector
#The vector in divided in bins and count the maximum value
def my_mode(vector,bins=50):
  dx = np.max(vector) - np.min(vector)
  dx = dx / bins
  b = np.sort(vector)
  i = 0
  j = 0
  o = 0
  limite = np.min(vector) + dx
  if ( dx > 1e-10 ): #the parameter is fixed
    while(o < len(b)):
        if ( b[o] < limite ):
            i = i + 1
            if ( i > j ):
                j = i
                maximo = limite - dx/2.0
            o = o + 1
        else:
            i = 0
            limite = limite + dx
  else:
      maximo = np.median(vector)

  return maximo

#-----------------------------------------------------------
def mode_and_99(vector):
    a = my_mode(vector)
    d, b, c = find_vals_perc(vector, sf=1.0, prob=99)

    return a,b,c

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
def good_clustering_fast(chi2,nconv,nwalkers):
  #Let us find the good indixes for the cluster
  #We have n walkers

  print 'STARTING CHAIN CLUSTERING'
  print 'Initial number of chains:', nwalkers

  #This variable will have each walker information
  chi2_walkers = [None]*nwalkers
  chi2_mean = [None]*nwalkers
  for i in range (0,nwalkers):
   chi2_walkers[i] = chi2[i::nconv]

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

  new_nwalkers = len(good_chain)

  print 'Final number of chains:', new_nwalkers
  return good_chain, new_nwalkers


def good_clustering_likelihood(like,nconv,nwalkers):
  #Let us find the good indixes for the cluster
  #We have n walkers

  print 'STARTING CHAIN CLUSTERING'
  print 'Initial number of chains:', nwalkers

  #This variable will have each walker information
  like_walkers = [None]*nwalkers
  like_mean = [None]*nwalkers
  for i in range (0,nwalkers):
   like_walkers[i] = like[i::nconv]

 #The mean of each walker
  for i in range (0,nwalkers):
    like_mean[i] = np.mean(like_walkers[i])

  #get the minimum chi2
  total_max = max(like_mean)

  good_chain = []
  #Let us kill all the walkers 5 times the minimum
  for i in range(0,nwalkers):
    if ( like_mean[i] > total_max * 0.9 ):
      #We are saving the good chain labels
      good_chain.append(i)

  new_nwalkers = len(good_chain)

  print 'Final number of chains:', new_nwalkers

  return good_chain, new_nwalkers

#-----------------------------------------------------------

def clustering(par,good_index):

  cluster_par = np.zeros(len(good_index))
  for i in range(0,len(good_index)):
    n = good_index[i]
    cluster_par[i] = par[n]

  return cluster_par

#-----------------------------------------------------------

def clustering_fast(par,good_index,nconv):

  dummy_par = []
  for i in good_index:
    dummy_par.append(par[i::nwalkers])

  cluster_par = np.ndarray(len(good_index)*nconv)

  n = 0
  for i in range(0,len(dummy_par)):
    for j in range(0,len(dummy_par[i])):
      cluster_par[n] = dummy_par[i][j]
      n = n + 1

  return cluster_par

#-----------------------------------------------------------
def print_values(vector,var,vartex,unit,unittex):
#fname is the variable where we are writting the numbers
#vector is the posterior vectors
#var is the label of the variable to save
#unit is the label of the variable to save
#vals do we want to print median or mode
  medv, minv, maxv = find_vals_perc(vector,s_factor)
  nd = 1
  if ( abs(minv) > 1e-20 and abs(maxv) > 1e-20 ):
    nd = int(np.log10(max(1./minv,1./maxv))) + 2
  opars.write('%10s = %4.7f - %4.7f + %4.7f %8s \n'%(var,medv,minv,maxv,unit))
  otex.write('\\newcommand{\\'+vartex+'}[1]['+unittex+'] \
  {'+str(round(medv,nd))+' _{ - '+str(round(minv,nd))+' } ^ { + '+str(round(maxv,nd))+' }$~#1} \n')
  if ( is_print_mode ):
    medv, minv, maxv = mode_and_99(vector)
    opars.write('%10s  %4.7f , %4.7f , %4.7f %8s \n'%('',medv,minv,maxv,unit))

#-----------------------------------------------------------
#         FIT JOINT RV-TRANSIT DATA
#-----------------------------------------------------------
def joint_fit():
  global fit_all, fit_ldc, fit_rvs, nt
  global a_from_kepler, mstar_mean, rstar_mean, mstar_sigma_rstar_sigma
  global is_log_P, is_ew, is_b_factor, is_log_k, is_log_rv0
  global fit_t0, fit_P, fit_e, fit_w, fit_i, fit_a,fit_q1, fit_q2, fit_rp, fit_k,fit_v0
  global T0,P,e,w,ii,a,q1,q2,rp,k0,alpha,beta, v0
  global min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w, min_i, max_i, min_a,\
         max_a, min_q1, max_q1, min_q1, max_q1, min_rp, max_rp, min_k, max_k, min_alpha, max_alpha, \
         min_beta, max_beta, min_rv0, max_rv0
  global vari,chi2,chi2red,t0o,Po,eo,wo,io,ao,q1o,q2o,rpo,ko,alphao,betao,vo, what_fit
  global new_nwalkers, good_index, nwalkers
  global jrvo, jtro, total_fit_flag, flags
  global limits, priorf, priorl, limits_ldc, limits_rvs


  if ( is_ew ):
    min_e = min_ew1
    max_e = max_ew1
    min_w = min_ew2
    max_w = max_ew2
    fit_e = fit_ew1
    fit_w = fit_ew2

  if ( is_b_factor ):
    min_i = min_b
    max_i = max_b
    fit_i = fit_b

  fit_all = [None]*8*nplanets
  for o in range(0,nplanets):
    fit_all[o*8:(o+1)*8] = [fit_t0[o],fit_P[o],fit_e[o],fit_w[o], \
                            fit_i[o],fit_a[o], fit_rp[o], fit_k[o] ]

  fit_rvs = []
  for o in range(0,nt):
    fit_rvs.append(fit_v0)

  fit_ldc = [fit_q1, fit_q2]

  fit_trends = [is_linear_trend,is_quadratic_trend]

  #Let us check what do we want to fit
  total_fit_flag = [ total_rv_fit, total_tr_fit ]

  flags = [is_log_P,is_ew,is_b_factor,is_den_a,is_log_k,is_log_rv0]


  vec_rv0_limits = []
  for m in range(0,nt):
    vec_rv0_limits.append(min_rv0[m])
    vec_rv0_limits.append(max_rv0[m])

  dummy_lims = [None]*8*2*nplanets

  for o in range(0,nplanets):

    dummy_lims[o*8*2:(o+1)*8*2 ] = \
    [ min_t0[o], max_t0[o], min_P[o], max_P[o], min_e[o], max_e[o], min_w[o], max_w[o] \
    , min_i[o], max_i[o], min_a[o], max_a[o], min_rp[o], max_rp[o], min_k[o], max_k[o] ]

  limits = dummy_lims

  limits_rvs = vec_rv0_limits

  limits_ldc = [ min_q1, max_q1, min_q2, max_q2]

  stellar_pars = [mstar_mean,mstar_sigma,rstar_mean,rstar_sigma]
  is_jitter = [is_jitter_rv, is_jitter_tr]

  if ( method == 'mcmc' ):

    #Ensure nwalkers is divisible by 2
    if ( nwalkers%2 != 0):
         nwalkers = nwalkers + 1

    pti.mcmc_stretch_move(\
    mega_time,mega_rv,megax,megay,mega_err,megae, \
    tlab,jrvlab,\
#    tlab,jrvlab,stellar_pars,a_from_kepler,\
    flags,total_fit_flag,is_jitter,fit_all,fit_rvs,fit_ldc,fit_trends, \
    nwalkers,maxi,thin_factor,nconv, limits, limits_rvs, \
    limits_ldc,n_cad, t_cad, npl=nplanets,n_tel=nt,n_jrv=n_jrv)

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

  newfile_trends = outdir+'/'+star+'_trends_data.dat'
  if ( os.path.isfile('trends_data.dat') ):
    os.rename('trends_data.dat',newfile_trends)

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
  for j in range(0,nplanets):
    oif.write ('------------------------------\n')
    oif.write ('  PLANET %s \n' %(star + plabels[j]))
    oif.write ('------------------------------\n')
    oif.write ('        PRIOR RANGES          \n')
    oif.write ('------------------------------\n')
    oif.write ('T0 = %s[ %4.4f , %4.4f ]\n' %(fit_t0[j],min_t0[j],max_t0[j]))
    oif.write ('P  = %s[ %4.4f , %4.4f ]\n' %(fit_P[j],min_P[j],max_P[j]))
    if ( is_ew ):
      oif.write ('ew1= %s[ %4.4f , %4.4f ]\n' %(fit_ew1[j],min_ew1[j],max_ew1[j]))
      oif.write ('ew2= %s[ %4.4f , %4.4f ]\n' %(fit_ew2[j],min_ew2[j],max_ew2[j]))
    else:
      oif.write ('e  = %s[ %4.4f , %4.4f ]\n' %(fit_e[j],min_e[j],max_e[j]))
      oif.write ('w  = %s[ %4.4f , %4.4f ]\n' %(fit_w[j],min_w[j],max_w[j]))
    if ( is_b_factor ):
      oif.write ('b  = %s[ %4.4f , %4.4f ]\n' %(fit_b[j],min_b[j],max_b[j]))
    else:
      oif.write ('i  = %s[ %4.4f , %4.4f ]\n' %(fit_i[j],min_i[j],max_i[j]))
    oif.write ('a  = %s[ %4.4f , %4.4f ]\n' %(fit_a[j],min_a[j],max_a[j]))
    oif.write ('rp = %s[ %4.4f , %4.4f ]\n' %(fit_rp[j],min_rp[j],max_rp[j]))
    oif.write ('K  = %s[ %4.4f , %4.4f ]\n' %(fit_k[j],min_k[j],max_k[j]))
    oif.write ('------------------------------\n')
  oif.write (' Other parameter priors \n')
  oif.write ('------------------------------\n')
  oif.write ('q1 = %s[ %4.4f , %4.4f ]\n' %(fit_q1,min_q1,max_q1))
  oif.write ('q2 = %s[ %4.4f , %4.4f ]\n' %(fit_q2,min_q2,max_q2))
  for m in range(0,nt):
    oif.write ('%s = %s[ %4.4f , %4.4f ]\n' %(telescopes_labels[m],fit_v0,min_rv0[m],max_rv0[m]))
  oif.write ('==============================\n')

  oif.close()
  dummy_file = open(out_init_file)
  for line in dummy_file:
    print line,
  dummy_file.close()

#------------------------------------------------------------------------#
#            Automatic creation of input for TANGO                       #
#------------------------------------------------------------------------#

def tango_params(param,vec,parhs=True):
    vlen = len(vec)
    letra = param + ' = '
    if ( parhs ):
     letra = letra + ' [ '
    for o in range(0,vlen):
        letra = letra + str(np.median(vec[o]))
        if (o < vlen - 1):
          letra = letra + ','

    if ( parhs ):
      letra = letra + ' ]'
    letra = letra + '\n'

    return letra

#This routine create an input file to create animations using tango
#github.com/oscaribv/tango
def create_tango_input():
  tangof = outdir+'/'+star+'_tango_input.py'
  tf = open(tangof,'w')

  tf.write('#Input file for tango\n')
  tf.write('#system:'+ star+'\n' )
  tf.write('#Created automatically with pyaneti\n')

  tf.write('\n')

  tf.write('#Data file with the flattened light curve\n')
  tf.write('lcname = \''+star+'_new_lc.dat\'\n')

  tf.write('#--------------------------------------------------------------------\n')
  tf.write('#                 Planet and orbit parameters\n')
  tf.write('# Each parameter is a list in which each element\n')
  if ( nplanets == 1 ):
    tf.write('# correspond to a planet. For this case, there is '+str(nplanets)+' planet\n')
  else:
    tf.write('# correspond to a planet. For this case, there are '+str(nplanets)+' planets\n')
  tf.write('#--------------------------------------------------------------------\n')

  tf.write('\n')

  #Orbital period (days)
  tf.write(tango_params('P',P_vec))
  tf.write(tango_params('T0',T0_vec))
  tf.write(tango_params('e',e_vec))
  tf.write(tango_params('w',w_vec))
  tf.write(tango_params('a',ar_vec))
  tf.write(tango_params('inclination',i_vec))
  tf.write(tango_params('rp',rr_vec))
  tf.write(tango_params('u1',[u1_vec],False))
  tf.write(tango_params('u2',[u2_vec],False))

##Integration time of the data
  tf.write('t_cad = ' + str(t_cad) +' \n')
  tf.write('n_cad = ' + str(n_cad) +' \n')

  tf.write('\n')

  #Calculate the transit duration for planet b to create the time ranges
  tfull = np.median(trt_vec[0])
  tfull = tfull/24.
  mit0  = np.median(T0_vec[0])

  tf.write('#--------------------------------------------------------------------\n')
  tf.write('#              Animation controls \n')
  tf.write('#--------------------------------------------------------------------\n')
  tf.write('#Window size to show the data (days)\n')
  tf.write('size_time = 0.5\n')
  tf.write('#1./(photograms per day) in this case the code will create a photogram each 7.2 min\n')
  tf.write('vel_time  = 1./200.\n')
  tf.write('#Animation minimum time (Be sure that you are using the same units as in your data file)\n')
  tf.write('tmin = '+str(mit0 - 2*tfull)+'\n')
  tf.write('#Animation maximum time (Be sure that you are using the same units as in your data file)\n')
  tf.write('tmax = '+str(mit0 + 2*tfull)+'\n')

  tf.write('\n')

  tf.write('#--------------------------------------------------------------------\n')
  tf.write('#                     Plot controls\n')
  tf.write('#--------------------------------------------------------------------\n')

  tf.write('\n')

  tf.write('#Control if we overplot the light curve model\n')
  tf.write('#You need to have installed pyaneti in your computer to use it\n')
  tf.write('is_plot_model = False\n')

  tf.write('\n')

  tf.write('#-----------------------------------------------------------------\n')
  tf.write('#                         END\n')
  tf.write('#-----------------------------------------------------------------\n')

