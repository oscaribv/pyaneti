#----------------------------------------------------------
#Example to run a multi-dimensional GP analysis with pyaneti
#Barragan O., 2021
#This analysis reproduces the full analysis of K2-100 presented in
#http://dx.doi.org/10.1093/mnras/stz2569
#----------------------------------------------------------

#----------------------------------------------------------
#This part is exactly as any other pyaneti run
#----------------------------------------------------------

nplanets = 1
#Filename including all the RV and activity indicators time-series
fname_rv = ['timeseries.dat']
#File name with the light curve data
#Note that there are four columns ordered as:
#time, normalised flux, normalised flux error, instrument/band label
fname_tr = ['multi_band_K2-100.dat']

#MCMC controls
thin_factor = 10
niter = 500
nchains = 100

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
is_multi_radius = False

method = 'mcmc'
#method = 'plot'

#K2-100 parameters as in Barragan et al., 2019
mstar_mean  = 1.15
mstar_sigma = 0.05
rstar_mean  = 1.24
rstar_sigma = 0.05
tstar_mean  = 5945.
tstar_sigma = 110.

unit_mass = 'earth'

fit_rv = [True]
fit_tr = [True]

is_jitter_tr = True
is_jitter_rv = True

#Prior section
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
fit_t0 = ['u']   
fit_P  = ['u']   
fit_e  = ['f']   
fit_w  = ['f']   
fit_ew1= ['f']   
fit_ew2= ['f']   
fit_b  = ['u']   
fit_a  = ['u']   
fit_rp = ['u']   
fit_k  = ['u']   
fit_v0 = 'u'    

#Now we have to fit LDC for each band, in this case we will set uniform priors for all bands
#This makes the trick to create a uniform prior for all 6 bands
fit_q1 = ['u']*6     #We fit q1 with gaussian priors
fit_q2 = ['u']*6     #We fit q2 with gaussian priors

#Set the prior limits for all bands
min_q1 = [0]*6
max_q1 = [1]*6
min_q2 = [0]*6
max_q2 = [1]*6

#There is a difference in the photometric and spectroscopic time-series of 4833. days
#We can solve this by addig 4833. to the light curve as
textra = 4833.

#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
min_t0  = [7141.70]
max_t0  = [7142.75]
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

#----------------------------------------------------------
# Here ends the part that is as any other pyaneti run
#----------------------------------------------------------

#----------------------------------------------------------
#In this part is where we start to add our GPs magic
#----------------------------------------------------------

#we are going to reproduce the result in Barragan et al., 2019, where the multi-gp approach as the form
# RV    = A_0 G + A_1 dG
# R_hk  = A_2 G + A_3 dG
# BIS   = A_4 G + A_5 dG
# with A_3 = 0 to recover original approach in Rajpaul et al. (2015)

#The Quasi-periodic kernel implemented in pyaneti is

#G(ti,tj) = A * exp[ sin(pi (ti - tj)/P_GP)**2 / (2*lambda_p**2) - (ti-tj)**2/(2*lamda_e**2) ]

#where A, lambda_e, lambda_p, and P_GP are the hyperparameters

#The first thing to note is how pyaneti deal with the RVs and the activity indicators
#The file timeseries.dat contains all data to be used in this example
#The data has to be provided as if each activity indicator were an extra instrument, i.e., 
#The activity indicator level has to be given in the fourth column of the data file, in this case timeseries.dat
#We are going to use RVs, log_RHK and Bisector span, they are labelled as INST, RHK, and BIS, respectively
#We indicate to pyaneti to read the data 
telescopes = ['INST','RHK','BIS']
#This vector has to be filled with the name of each telescope telescopes[i]
telescopes_labels = ['HARPS RV','HARPS log $R_{HK}$','HARPS BIS SPAN']


#We have to tell to pyaneti what correlation matrix we want for our RV analysis, this is done by setting
kernel_rv = 'MQ3'
# where 'MQ3' is the keyword for the quasiperiodic kernel in the multi-dimensional GP framework, with three time-series. 
#Pyaneti can deal with a different number of time-series by changing to MQX, where X is the number of time-series.

#now we have to treat the hyper parameters as "normal" pyaneti parameters
#i.e., we need to say what kind of prior we will use and the ranges

#First we need to tell to pyaneti what kind of priors we want, equivalent to the fit_* parameters previously used in this file
#The variable that pass this to pyaneti is called fit_krv. In this case we are going to deal with 9 parameters, the six A_i amplitudes
#plus the lambda_e, lambda_p and P_GP from the QP Kernel.
#The first elements of the fit_krv vector are always the A_i's and the last are the kernel hyperparameters.
#The fit_krv vector would be then
#kernel_rv = 'None'
fit_krv = [None]*9 #Define the list with 9 elements
fit_krv[0] = 'u' #A_0
fit_krv[1] = 'u' #A_1
fit_krv[2] = 'u' #A_2
fit_krv[3] = 'f' #A_3, we fix it to zero to recover Rajpaul et al., 2015 approach
fit_krv[4] = 'u' #A_4
fit_krv[5] = 'u' #A_5
fit_krv[6] = 'u' #lambda_e
fit_krv[7] = 'u' #lambda_p
fit_krv[8] = 'u' #P_GP

#We have already indicated the kind of prior that we want, the next step is to set the prior ranges
#The prior ranges are stored in a list called krv_priors, this list lenght is two times the lenght of fit_krv
#It follows the same logic as fit_krv, the first elements correspond to the ranges of the A_i's variables and then the kernel variables
#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
krv_priors = [     
        0.0,0.5,   #ranges for A_0
        0.0,0.5,   #ranges for A_1
        0.0,0.5,   #ranges for A_2
        0.0,0.0,   #ranges for A_3, note that this one is fixed to zero
        -1.5,1.5,  #ranges for A_4
        -1.5,1.5,  #ranges for A_5
        1,80,      #ranges for lambda_e
        0.1,2.0,   #ranges for lambda_p
        4.,5.1     #ranges for P_GP
        ]

#----------------------------------------------------------
#END
#----------------------------------------------------------
