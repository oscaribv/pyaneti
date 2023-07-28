##################################################################
#                         K2_233_final
#                    Oscar Barragan, Jul 2023
# This input file reproduces the final pyaneti joint RV and transit
# model of K2-233 as in Barragan et al., 2023, MNRAS, 522, 3458
# If you use this file to create the setup of your analysis
# please cite the paper.
##################################################################

#Indicate the RV and transit files
fname_rv = ['spec_time_series.dat']
fname_tr = ['lc.dat']

#Specify the number of planets, in our case 3
nplanets = 3


#MCMC controls
#the thin factor for the chains
thin_factor = 10
#The number of iterations to be taken into account
#The TOTAL number of iterations for the burn-in phase is thin_factor*niter
niter       = 500
#Number of independent Markov chains for the ensemble sampler
nchains     = 250

#Choose the method that we want to use
# mcmc -> runs the mcmc fit program
# plot -> this option create the plots only if a previus run was done
method = 'mcmc'
#method = 'plot'

#We are working with Kepler long-cadence data, integrated over 30 min
#So we need to integrate the transit model as well
#We integrate over 30 steps. See Kipping 2013 for more details.
n_cad = [30]             #30 steps
t_cad = [30./60./24.]    #30 min
#Typically TESS, CHEOPS or ground-base data are short-cadence so they do not need
#the model to be resampled, but always check!

#Give the stellar parameters are given by Lillo-Box et al., (2020)
mstar_mean  = 0.79
mstar_sigma = 0.01
rstar_mean  = 0.71
rstar_sigma = 0.01
tstar_mean  = 5033
tstar_sigma = 50.

#Give the J mag to compute the TSM (Kempton 2018)
mag_j = 8.968


#Print the derived parameters in Earth units
unit_mass = 'earth'

#We have three planets, and we want to fit RVs and transit for all of them
#Let's tell pyaneti that we want to fill rv and tr for all
fit_rv = [True,True,True]
fit_tr = [True,True,True]

##################################################################
#Prior section
##################################################################

#The priors are set as in Table 3, Barragán et al., 2023
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
# b -> Beta distribution priors
fit_t0 = ['u','u','u']   #We fit for t0 with uniform priors
fit_P  = ['u','u','u']   #We fit for P with uniform priors
fit_e  = ['f','f','b']   #We fix e, it works only if is_ew = False
fit_w  = ['f','f','u']   #We fix w, it works only if is_ew = False
fit_b  = ['u','u','u']   #We fix the impact factor
fit_a  = ['u','u','u']   #We fit a with gaussian priors (given by the stellar parameters)
fit_rp = ['u','u','u']   #We fit rp with uniform priors
fit_k  = ['u','u','u']   #We fit k with uniform priors
fit_v0 = 'u'             #We fit systemic velicities with uniform priors
fit_q1 = ['u']     #We fit q1 with uniform priors
fit_q2 = ['u']     #We fit q2 with uniform priors

##################################################################
#Prior values
##################################################################

#q1 and q2 have uniform priors between 0 and 1
min_q1 = [0]
max_q1 = [1]
min_q2 = [0]
max_q2 = [1]

#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
#if 'b' is selected for the parameter A, sets a beta distribution with shape parameters Beta(min_A,max_A)

#Set the priors for all the parameters for the 3 planets

#time of mid-transit, days
min_t0  = [7991.6745,7586.8425,8005.5640]
max_t0  = [7991.7077,7586.9105,8005.6008]
#orbital period, days
min_P   = [2.4667,7.0547,24.3509]
max_P   = [2.4683,7.0655,24.3781]
#impact parameter
min_b   = [0.0]*nplanets
max_b   = [1.0]*nplanets
#Dopler semi-amplitude, km/s
min_k   = [0.0]*nplanets
max_k   = [0.01]*nplanets
#Scaled planetary radii
min_rp  = [0.0]*nplanets
max_rp  = [0.1]*nplanets
#Eccentricity, See Barragán et al., 2023 to understand these priors
min_e   = [0,0,1.52]
max_e   = [1,1,29.0]
#Angle of periastron
min_w = [0.0]*nplanets
max_w = [2*np.pi]*nplanets

#We are sampling for the stellar density instead of the
#scaled semi-major axes of the three planets
sample_stellar_density = True

if sample_stellar_density:
    #If we want to set a Gaussian prior based on the stellar parameters
    #we could so something like this
    #fit_a = ['g']*nplanets
    #min_a   = [3.11]*nplanets
    #max_a   = [0.14]*nplanets
    #But in our case we are giving an uniform prior on the stellar density
    #So we choose
    fit_a = ['u']*nplanets
    min_a   = [1.]*nplanets  #g/cm^3
    max_a   = [5.]*nplanets  #g/cm^3

#We sample for a jitter term for RV and transit data
is_jitter_rv = True
is_jitter_tr = True

##################################################################
#Multi-GP part
##################################################################

#This vector has to be filled with the label that we use for each telescope in the RV data file
telescopes_labels = ['RV','FWHM','BIS']
#This vector has to be filled with the name of each telescope telescopes[i]
#In the file with spectroscopic data 0 is RV, 1 is FWHM and 2 is BIS
telescopes = ['0','1','2']

#3-dimensional GPs with quasi-periodic kernel
kernel_rv = 'MQ3'
#We need 9 parameters, 6 amplitudes + 3 QP hyper-parameters
#See https://github.com/oscaribv/pyaneti/blob/master/inpy/example_timeseries_k2100/input_fit.py
#For more details on how to turn on a multi-GP in pyaneti
fit_krv = ['f']*9
fit_krv[0] = 'u'
fit_krv[1] = 'u'
fit_krv[2] = 'u'
fit_krv[3] = 'f'
fit_krv[4] = 'u'
fit_krv[5] = 'u'
krv_priors = [-0.1,0.1,-0.1,0.1,0.0,0.5,0.0,0.5,-0.5,0.5,-0.5,0.5]

fit_krv[6] = 'u'
fit_krv[7] = 'u'
fit_krv[8] = 'u'
QP_priors = [5.,50.0,0.1,1.,8,11]

krv_priors = np.concatenate([krv_priors,QP_priors])


##################################################################
# END
##################################################################
