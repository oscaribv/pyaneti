#----------------------------------------------------------
#Example to run a multi-dimensional GP analysis with pyaneti
#Barragan O., 2021
#This file reproduces the analysis presented in Sect. 5.2 of the
#pyaneti II paper
#----------------------------------------------------------

fname_rv = ['2inst_data.dat']

thin_factor = 10
niter = 500
nchains = 100

method = 'mcmc'
#method = 'plot'

plot_rv_std = True

rv_xlabel = "Time (days)"

unit_mass = 'earth'

#For this example we are modelling two planetary signals
nplanets = 2
#To run a multi-dimensional GP approach in the spectroscopic time-series
#we need to turn on the RV fit
fit_rv = [True,True]
fit_tr = [False,False]

#We do not fit for a jitter term for this example
is_jitter_rv = False

#WE will sample the eccentricity using the 
#sqrt(e) cos w and sqrt(e) sin w parametrisations
#This is turned on with is_ew = True
is_ew = True

#Prior section
#Set the priors type for both Keplerian signals
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
fit_t0 = ['g','g']   
fit_P  = ['g','g']   
fit_ew1= ['u','u']   
fit_ew2= ['u','u']   
fit_b  = ['f','f']   
fit_a  = ['f','f']   
fit_rp = ['f','f']   
fit_k  = ['u','u']   

#This has to be activated to fit for the offsets for all the time-series
#Pyaneti selects the prior ranges automatically based on the data
fit_v0 = 'u'    

#Set the priors values for both Keplerian signals
#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
min_t0  = [1,2] 
max_t0  = [1e-3,1e-3]  
min_P   = [3,10]
max_P   = [1e-3,1e-3]
min_a   = [1.1]*nplanets
max_a   = [15.]*nplanets
min_b   = [0.0]*nplanets
max_b   = [1.0]*nplanets
min_k   = [0.0]*nplanets
max_k   = [0.05]*nplanets
min_rp  = [0.0]*nplanets
max_rp  = [0.05]*nplanets
min_ew1 = [-1,-1]*nplanets
max_ew1 = [1,1]*nplanets
min_ew2 = [-1,-1]*nplanets
max_ew2 = [1,1]*nplanets
min_e = [0,0]*nplanets
max_e = [1,1]*nplanets
min_w = [0,0]*nplanets
max_w = [2*np.pi,2*np.pi]*nplanets

#----------------------------------------------------------
# Here ends the part that is as any other pyaneti run
#----------------------------------------------------------

#----------------------------------------------------------
#In this part is where we start to add our GPs magic
#----------------------------------------------------------

#we are going to reproduce the result in Barragan et al., 2021, where the multi-gp approach as the form
# S1    = A_1 G + B_1 dG
# S2    = A_2 G + B_2 dG

#The Quasi-periodic kernel implemented in pyaneti is

#G(ti,tj) = A * exp[ sin(pi (ti - tj)/P_GP)**2 / (2*lambda_p**2) - (ti-tj)**2/(2*lamda_e**2) ]

#where A, lambda_e, lambda_p, and P_GP are the hyperparameters

#The first thing to note is how pyaneti deal with the RVs and the activity indicators
#The file 2inst_data.dat contains all data to be used in this example
#The data has to be provided as if each activity indicator were an extra instrument, i.e., 
#The activity indicator level has to be given in the fourth column of the data file
#We are going to use 'rv_i1','rv_i2','s2_i1','s2_i2'
#We indicate to pyaneti to read the data 
telescopes = ['rv_i1','rv_i2','s2_i1','s2_i2']
#This vector has to be filled with the name of each telescope telescopes[i]
telescopes_labels = ['RV I1','RV I2','S2 I1', 'S2 I2']


#We have to tell to pyaneti what correlation matrix we want for our RV analysis, this is done by setting
kernel_rv = 'MQ2'
# where 'MQ2' is the keyword for the quasiperiodic kernel in the multi-dimensional GP framework, with two time-series. 
#Pyaneti can deal with a different number of time-series by changing to MQX, where X is the number of time-series.

#WE NOTE THAT WE HAVE 4 DIFFERENT LABELS IN TELESCOPES, BUT PYANETI KNOWS HOW MANY TIME-SERIES TO USE IN THE
#MULTI-GP REGRESSION VIA THE kernel_rv VARIABLE.

#now we have to treat the hyper parameters as "normal" pyaneti parameters
#i.e., we need to say what kind of prior we will use and the ranges

#First we need to tell to pyaneti what kind of priors we want, equivalent to the fit_* parameters previously used in this file
#The variable that pass this to pyaneti is called fit_krv. In this case we are going to deal with 7 parameters, the four A_i amplitudes
#plus the lambda_e, lambda_p and P_GP from the QP Kernel.
#The first elements of the fit_krv vector are always the amplitudes and the last are the kernel hyperparameters.
#The fit_krv vector would be then
fit_krv = [None]*7 #Define the list with 9 elements
fit_krv[0] = 'u' #A_1
fit_krv[1] = 'u' #B_1
fit_krv[2] = 'u' #A_2
fit_krv[3] = 'u' #B_2
fit_krv[4] = 'u' #lambda_e
fit_krv[5] = 'u' #lambda_p
fit_krv[6] = 'u' #P_GP

#We have already indicated the kind of prior that we want, the next step is to set the prior ranges
#The prior ranges are stored in a list called krv_priors, this list lenght is two times the lenght of fit_krv
#It follows the same logic as fit_krv, the first elements correspond to the ranges of the amplitude variables and then the kernel variables
#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A
krv_priors = [     
        0.0,0.5,    #ranges for A_1
        0.0,0.5,    #ranges for B_1
        -0.0,0.5,   #ranges for A_2
        -0.5,0.5,   #ranges for B_2
        1,200,      #ranges for lambda_e
        0.01,2.0,   #ranges for lambda_p
        4.8,5.2     #ranges for P_GP
        ]

#----------------------------------------------------------
#END
#----------------------------------------------------------
