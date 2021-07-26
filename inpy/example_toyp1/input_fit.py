#----------------------------------------------------------
#Example to run a multi-dimensional GP analysis with pyaneti
#Barragan O., 2021
#This file reproduces the analysis presented in Sect. 5.1 of the
#pyaneti II paper
#----------------------------------------------------------


nplanets = 1
fname_rv = ['data_3mdgp.dat']

thin_factor = 10
niter = 500
nchains = 100

method = 'mcmc'

is_plot_correlations = True

rv_xlabel = "Time (days)"

#To run a multi-dimensional GP approach in the spectroscopic time-series
#we need to turn on the RV fit
fit_rv = [True]
fit_tr = [False]

#Let us fit for a jitter term for the RV data
is_jitter_rv = True

#Prior section for the planet parameters
#In this example we are not modelling a planetary signal
#therefore, we fix all the parameter values
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
fit_t0 = ['f']   
fit_P  = ['f']   
fit_e  = ['f']   
fit_w  = ['f']   
fit_ew1= ['f']   
fit_ew2= ['f']   
fit_b  = ['f']   
fit_a  = ['f']   
fit_rp = ['f']   
#We fix the amplitude of the planetary signal to zero
#in this way pyaneti will not add a planetary signal in the first time-series
fit_k  = ['f']   

#This has to be activated to fit for the offsets
fit_v0 = 'u'    

#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
min_t0  = [7140.71934] 
max_t0  = [0.00027]  
min_P   = [1.6739038]
max_P   = [0.0000004]
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

#we are going to reproduce the result in Barragan et al., 2021, where the multi-gp approach as the form
# RV    = A_0 G + B_0 dG
# R_hk  = A_1 G + B_1 dG
# BIS   = A_2 G + B_2 dG
# with A_3 = 0 to recover original approach in Rajpaul et al. (2015)

#The Quasi-periodic kernel implemented in pyaneti is

#G(ti,tj) = A * exp[ sin(pi (ti - tj)/P_GP)**2 / (2*lambda_p**2) - (ti-tj)**2/(2*lamda_e**2) ]

#where A, lambda_e, lambda_p, and P_GP are the hyperparameters

#The first thing to note is how pyaneti deal with the RVs and the activity indicators
#The file data_3mdgp.dat contains all data to be used in this example
#The data has to be provided as if each activity indicator were an extra instrument, i.e., 
#The activity indicator level has to be given in the fourth column of the data file
#We are going to use RVs, s2 and s3 as given in the data_3mdgp.dat file
#We indicate to pyaneti to read the data, every time-series has to be passed as a telescope label 
telescopes = ['rvs','s2','s3']
#This vector has to be filled with the label for each instrument, in this case S1, S2 and S3
telescopes_labels = ['$S_1$','$S_2$','$S_3$']


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
fit_krv[1] = 'u' #B_0
fit_krv[2] = 'u' #A_1
fit_krv[3] = 'u' #B_1
fit_krv[4] = 'u' #A_2
fit_krv[5] = 'u' #B_2
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
        0.0,0.5,   #ranges for B_0
        -0.0,0.5,  #ranges for A_1
        -0.5,0.5,  #ranges for B_1
        -0.5,0.5,  #ranges for A_2
        -0.5,0.5,  #ranges for B_2
        1,80,      #ranges for lambda_e
        0.01,2.0,  #ranges for lambda_p
        4.,6.      #ranges for P_GP
        ]

#----------------------------------------------------------
#END
#----------------------------------------------------------
