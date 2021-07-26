#Input file for multi-band modelling with pyaneti
#Created by Oscar Barragan, July 2021

#The multi-band aproach is similar to the multi-instrument RV fit, where we have to include the band (instrument)
#label in the fourth column in the input file

#File name with the light curve data
#Note that there are four columns ordered as:
#time, normalised flux, normalised flux error, instrument/band label
fname_tr = ['multiband.dat']

#Define a vector with the bands 
#In this case
bands = ['b1','b2']

#Pyaneti is able to deal with multi band with a different cadence for each band
#In this example both  bands are short candece, and we do not need to resample the model
#We need to define t_cad, a vector indicating the integration time (in units of days), 
#Note that when no resample of the model is needed, the value that we put in t_cad is not important, but we still need to give a value
t_cad = [1.5 / 60. / 24.0,1.5 / 60. / 24.0]
#Vector indicating the number of steps to integrate the model, each element has to correspond to the bands defined in bands variable
n_cad = [1,1]

#Label for the time-series plot
tr_xlabel = "Time (days)"

#This variable controls if we want to fit a single radius for all bands or a radius fit for each band
#If False (default) pyaneti will fit a single planet radius for all bands
#If True  pyaneti will fit a planet radius for EACH band, this might be useful to check if the planet radius is consistent within all bands
is_multi_radius = False

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

plot_binned_data = True


#What units do you prefer for your planet parameters?
# earth, jupiter or solar
unit_mass = 'jupiter'

#Let us indicate to pyaneti to fit model two planets
nplanets = 2
#this is only a transit fit, so we turn off the RV modelling
fit_rv = [False,False]
#and we turn on the transit modelling
fit_tr = [True,True]


#Prior section
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
fit_t0 = ['u','u']   #We fit for t0 with uniform priors
fit_P  = ['u','u']   #We fit for P with uniform priors
fit_e  = ['f','f']   #We fix e, it works only if is_ew = False
fit_w  = ['f','f']   #We fix w, it works only if is_ew = False
fit_ew1= ['f','f']   #We fix e, it works only if is_ew = False
fit_ew2= ['f','f']   #We fix w, it works only if is_ew = False
fit_b  = ['u','u']   #We fix the impact factor
fit_a  = ['u','u']   #We fit a with gaussian priors (given by the stellar parameters)
fit_rp = ['u','u']   #We fit rp with uniform priors
fit_k  = ['f','f']   #We fit k with uniform priors
fit_v0 = 'f'     #We fit systemc velicities with uniform priors

#We will sample for the stellar density instead of sampling for the scaled semi-major axis
#We recover the semi-major axis for each planet using the approximation in Winn 2010
sample_stellar_density = True

#Now we have to fit LDC for each band, 
#in this case we will set uniform priors for all bands
#First element corresponds to band 1, second to  band 2, and so on
fit_q1 = ['u','u']    #We fit q1 with gaussian priors
fit_q2 = ['u','u']    #We fit q2 with gaussian priors

#Set the prior limits for all bands
#First element corresponds to band 1, second to  band 2, and so on
min_q1 = [0,0]
max_q1 = [1,1]
min_q2 = [0,0]
max_q2 = [1,1]

#We do not fit for a jitter term in this example
is_jitter_tr = False

#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
min_t0  = [3.95,2.95]
max_t0  = [4.05,3.05]
min_P   = [2.95,9.95]
max_P   = [3.05,10.05]
min_b   = [0.0]*nplanets
max_b   = [1.0]*nplanets
min_k   = [0.0]*nplanets
max_k   = [0.05]*nplanets
min_rp  = [0.0]*nplanets
max_rp  = [0.15]*nplanets
min_e   = [0]*nplanets
max_e   = [1]*nplanets
min_w   = [np.pi/2]*nplanets
max_w   = [np.pi/2]*nplanets

#We are sampling for the stellar density, so let us indicate the priors
#between 0.1 and 5 g/cm^3
min_a = [0.01]*nplanets
max_a = [5.]*nplanets

