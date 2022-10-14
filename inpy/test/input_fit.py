#Input file for test problem for pyaneti
#Created by Barragan O.
#This file contains the basic parameters to change for simple pyaneti run
#There are more flags which can control the input/output of the code
#They can be found inside src/default.py, all the flags inside this file
#can be changed inside the input file.

#Telescope labels
#This vector has to be filled with the label that we use for each telescope in the RV data file
telescopes = ['S']
#This vector has to be filled with the name of each telescope telescopes[i]
telescopes_labels = ['Super_telescope']

#Input files
#fname_rv contains the RV data
fname_rv = ['earth_rv.dat']
#fname_tr contains the transit data
fname_tr = ['earth_lc.dat']

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

#Define the star parameters to calculate the planet parameters
mstar_mean  = 1.0
mstar_sigma = 0.1
rstar_mean  = 1.0
rstar_sigma = 0.1
tstar_mean  = 5772.0
tstar_sigma = 50.

#What units do you prefer for your planet parameters?
# earth, jupiter or solar
unit_mass = 'earth'

#If we want posterior, correlation and/or chain plots these options have to be set True
is_plot_posterior    = True
is_plot_correlations = True
is_plot_chains       = False

#Control the labels of the time-series plots
rv_xlabel = "BJD - 2440000 (days)"
tr_xlabel = "BJD - 2440000 (days)"

#We want to fit transit and RV 
#For a pure RV fit, fit_tr has to be False
#For a pure TR fit, fit_rv has to be False
#For multi-planet fits fit_rv and fit_tr have the form [True,True,False,...]
#one element for each planet.
fit_rv = [True]
fit_tr = [True]

#is_ew controls the parametrization sqrt(e)sin w and sqrt(e) cos w
#if True we fit for the parametrization parameters
#if False we fit for e and w
#Default is True
is_ew = True

#Prior section
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
fit_t0 = ['u']   #We fit for t0 with uniform priors
fit_P  = ['u']   #We fit for P with uniform priors
fit_ew1= ['f']   #We fix sqrt(e) sin w, it works only if is_ew = True
fit_ew2= ['f']   #We fix sqrt(e) cos w, it works only if is_ew = True
fit_e  = ['f']   #We fix e, it works only if is_ew = False
fit_w  = ['f']   #We fix w, it works only if is_ew = False
fit_b  = ['f']   #We fix the impact factor
fit_a  = ['u']   #We fit a with gaussian priors (given by the stellar parameters)
fit_rp = ['u']   #We fit rp with uniform priors
fit_k  = ['u']   #We fit k with uniform priors
fit_v0 = 'u'     #We fit systemc velicities with uniform priors
fit_q1 = 'g'     #We fit q1 with gaussian priors
fit_q2 = 'g'     #We fit q2 with gaussian priors

#if you are in my_test uncomment the next line
#fit_b  = ['u']   #We fit the impact factor

#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
min_t0  = [8285.05] 
max_t0  = [8285.15]  
min_P   = [365.206]
max_P   = [365.306]
min_ew1 = [0.0]
min_ew2 = [0.0]
max_ew1 = [1.0]
max_ew2 = [1.0]
min_a   = [200.]
max_a   = [250.]
min_b   = [0.0]
max_b   = [1.0]
min_k   = [0.0]
max_k   = [0.001]
min_rp  = [0.0]
max_rp  = [0.1]
min_q1  = 0.3464
max_q1  = 0.05
min_q2  = 0.2839
max_q2  = 0.05
