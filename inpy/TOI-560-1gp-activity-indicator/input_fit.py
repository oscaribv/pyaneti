#----------------------------------------------------------------------------------------
#Example on how to fit 1D GP models with a constant mean function
#This means that the code will not fit for planetary signals in the first dimension
#The trick is to fit for a planet with Doppler semi-amplitude of zero
#This can be useful to fit 1D GP models of activity indicators
#This example shows how to use it for the FWHM but it can be easily changed to any time series
#Oscar Barragan, March, 2025
#----------------------------------------------------------------------------------------

#Filename that includes all the time series with the corresponding label
fname_rv = ['timeseries.dat']

#nplanets = 1, if = 0 the code does not work
nplanets = 1

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

#We still need to specify that we are fitting "one planet" 
fit_rv = [True]

#We fix all the planetary parameters
#The values of the orbital parameters are not relevant if 
#The Doppler semi-amplitude = 0
fit_t0 = ['f']   #We fit for t0 with uniform priors
fit_P  = ['f']   #We fit for P with uniform priors
fit_e  = ['f']   #We fix e, it works only if is_ew = False
fit_w  = ['f']   #We fix w, it works only if is_ew = False
fit_a  = ['f']   #We fit a with gaussian priors (given by the stellar parameters)
fit_k  = ['f']   #We fit k with uniform priors

#WE NEED TO SAMPLE FOR OFFSETS
fit_v0 = 'u'    

#We are also sampling for a jitter term
is_jitter_rv = True

#Fixing the Doppler semi-amplitude to zero
min_k   = [0.0]
max_k   = [0.05]


#In our input file, the FWHM time series has the label 2
telescopes = ['2']
#We re-name it for a clearer view in the plot
telescopes_labels = ['FWHM']

#In this context where we are working with the multi-GP
#it is better to work with the multi-GP keywords, so even for a 1D GP
#The label is MQ1, in this way pyaneti knows that there is a dG term in
#the modelling, even if fixed to zero
kernel_rv = 'MQ1'
fit_krv = ['f']*5
fit_krv[0] = 'u' #We allow to sample for the first amplitude of the GP
#There is no point on sampling for a dG term in 1D GP regresions
fit_krv[1] = 'f' #We fix the amplitude of the derivative to zero
#We set the range for the first amplitude and fix the second one to zero
krv_priors = [0,0.1,0.0,0.1]

#We now set the priors for the Kernel hyperparameters, in this case the QP kernel
fit_krv[2] = 'u'
fit_krv[3] = 'u'
fit_krv[4] = 'u'
QP_priors = [1.,500.0,0.1,2.0,11.5,13]

krv_priors = np.concatenate([krv_priors,QP_priors])

#The setup is ready to run
