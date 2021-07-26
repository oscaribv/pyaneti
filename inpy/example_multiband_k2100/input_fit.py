#Input file for multi-band fitting of K2-100b
#As it appears in Barragan et al., 2019 (https://academic.oup.com/mnras/article/490/1/698/5569669) 
#Created by Oscar Barragan, March 2021.

#The multi-band aproach is similar to the multi-instrument RV fit, where we have to include the band (instrument)
#label in the fourth column in the input file

#File name with the light curve data
#Note that there are four columns ordered as:
#time, normalised flux, normalised flux error, instrument/band label
fname_tr = ['multi_band_K2-100.dat']

#Define a vector with the bands 
#In this case
#K2 - K2 C5 long cadence data
#G  - ARCTIC data
#SC - K2 C18 short cadence data
#band1 - MUSCAT i data
#band2 - MUSCAT r data
#band3 - MUSCAT z data
bands = ['K2','G','SC','band1','band2','band3']

#Pyaneti is able to deal with multi band with a different cadence for each band
#In this example K2 C5 data was observed in long cadence (29.425 min) we need to resample the model
#All other bands do not need resampling
#We need to define t_cad, a vector indicating the integration time (in units of days), 
#each element has to correspond to the bands defined in bands variable
#Note that when no resample of the model is needed, the value that we put in t_cad is not important, but we still need to give a value
t_cad = [29.425 / 60. / 24.0,1.5 / 60. / 24.0,1.5 / 60. / 24.0,1.5 / 60. / 24.0,1.5 / 60. / 24.0,1.5 / 60. / 24.0]
#Vector indicating the number of steps to integrate the model, each element has to correspond to the bands defined in bands variable
#We only integrate the model for the K2 C5 long cadence with 10 steps
n_cad = [10,1,1,1,1,1]

#This variable controls if we want to fit a single radius for all bands or a radius fit for each band
#If False (default) pyaneti will fit a single planet radius for all bands
#If True  pyaneti will fit a planet radius for EACH band, this might be useful to check if the planet radius is consistent within all bands
is_multi_radius = True

#MCMC controls
#the thin factor for the chains
thin_factor = 10
#The number of iterations to be taken into account
#The TOTAL number of iterations for the burn-in phase is thin_factor*niter
niter       = 500
#Number of independent Markov chains for the ensemble sampler
nchains     = 100

#Choose the method that we want to use
# mcmc -> runs the mcmc fit program
# plot -> this option create the plots only if a previus run was done
method = 'mcmc'
#method = 'plot'

#If you want a plot with the seaborn library, is_seaborn_plot has to be True
is_seaborn_plot = True

#K2-100 parameters as in Barragan et al., 2019
mstar_mean  = 1.15
mstar_sigma = 0.05
rstar_mean  = 1.24
rstar_sigma = 0.05
tstar_mean  = 5945.
tstar_sigma = 110.
mag_j =  9.463

#What units do you prefer for your planet parameters?
# earth, jupiter or solar
unit_mass = 'earth'

fit_rv = [False]
fit_tr = [True]

n_columns_posterior = 5

#Prior section
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
fit_t0 = ['u']   #We fit for t0 with uniform priors
fit_P  = ['u']   #We fit for P with uniform priors
fit_e  = ['f']   #We fix e, it works only if is_ew = False
fit_w  = ['f']   #We fix w, it works only if is_ew = False
fit_ew1= ['f']   #We fix e, it works only if is_ew = False
fit_ew2= ['f']   #We fix w, it works only if is_ew = False
fit_b  = ['u']   #We fix the impact factor
fit_a  = ['u']   #We fit a with gaussian priors (given by the stellar parameters)
fit_rp = ['u']   #We fit rp with uniform priors
fit_k  = ['u']   #We fit k with uniform priors
fit_v0 = 'u'     #We fit systemc velicities with uniform priors

#Now we have to fit LDC for each band, in this case we will set uniform priors for all bands
#This makes the trick to create a uniform prior for all 6 bands
fit_q1 = ['u']*6     #We fit q1 with gaussian priors
fit_q2 = ['u']*6     #We fit q2 with gaussian priors

#Set the prior limits for all bands
min_q1 = [0]*6
max_q1 = [1]*6
min_q2 = [0]*6
max_q2 = [1]*6

#We fit for a jitter term for each band
is_jitter_tr = True

#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
min_t0  = [2307.70]
max_t0  = [2307.75]
min_P   = [1.6737]
max_P   = [1.6740]

min_a   = [1.1]
max_a   = [15.]
min_b   = [0.0]
max_b   = [1.0]
min_k   = [0.0]
max_k   = [0.05]
min_rp  = [0.0]
max_rp  = [0.05]
min_w   = [np.pi/2]
max_w   = [np.pi/2]
