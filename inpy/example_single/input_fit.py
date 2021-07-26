#Input file for single transit fit
#Created by Barragan O., March 2021

#Input files
#fname_tr contains the transit data
fname_tr = ['lc_single.dat']

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


#number of planets to model, in tis  case it is only one
nplanets = 1

#We do not sample for a jitter term in this case
is_jitter_tr = False

#What units do you prefer for your planet parameters?
# earth, jupiter or solar
unit_mass = 'earth'

#If we want posterior, correlation and/or chain plots these options have to be set True
is_plot_posterior    = True
is_plot_correlations = True
is_plot_chains       = True

#We want to fit transit and RV 
#For a pure RV fit, fit_tr has to be False
#For a pure TR fit, fit_rv has to be False
#For multi-planet fits fit_rv and fit_tr have the form [True,True,False,...]
#one element for each planet.
#In this case we are only interested in modelling transit
fit_rv = [False]
fit_tr = [True]


#THIS IS THE MOST IMPORTANT LINE TO ADD IN A SINGLE TRANSIT MODEL
#We need to tell pyaneti that we will be fitting only one transit
#This will tell to the code how to deal with the period and the dummy semi-amplitude
is_single_transit = True

#Prior section
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
fit_t0 = ['u']   #We fit for t0 with uniform priors
fit_P  = ['u']   #We fit for P with uniform priors
fit_e  = ['f']   #We fix e, it works only if is_ew = False
fit_w  = ['f']   #We fix w, it works only if is_ew = False
fit_b  = ['u']   #We fix the impact factor
fit_a  = ['u']   #We fit a with gaussian priors (given by the stellar parameters)
fit_rp = ['u']   #We fit rp with uniform priors
fit_k  = ['f']   #We fit k with uniform priors
fit_v0 = 'f'     #We fit systemc velicities with uniform priors
fit_q1 = 'u'     #We fit q1 with gaussian priors
fit_q2 = 'u'     #We fit q2 with gaussian priors

#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
min_t0  = [9.9] 
max_t0  = [10.1]  
min_P   = [18.]
max_P   = [35.]
min_a   = [1.1]
max_a   = [50.]
min_b   = [0.0]
max_b   = [1.1]
min_k   = [0.0]
max_k   = [0.05]
min_rp  = [0.0]
max_rp  = [0.03]
min_e = [0.]
min_w = [np.pi/2.]
max_e = [1.0]
max_w = [1.0]
