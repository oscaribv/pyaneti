#Input file for test_2p
#30, Jun, 2016
#Barragan, O.

#telescopes labels
telescopes = ['S']
telescopes_labels = ['Telescope']

#Input files
fname_rv = 'two_planets.dat'

#chain parameters
thin_factor = 10
nconv = 500
nwalkers = 250

method = 'sm'
#method = 'plot'

#Output controls
unit_mass = 'jupiter'

#Star parameters
mstar_mean = 1.0
rstar_mean = 1.0
inclination = np.pi/2

#We want to fit transit and RV 
fit_rv = True

nplanets = 2

#What do we want to fit?
fit_t0 = [True,True]
fit_P  = [True,True]
fit_e  = [True,True]
fit_w  = [True,True]
fit_k  = [True,True]
fit_v0 = [True,True]
fit_alpha = [False,False]
fit_beta = [False,False]

#Priors
T0 = [7303.023611,7567.400694]
P  = [176.251,264.377]
e  = [1e-3,1e-3]
w  = [np.pi/2,np.pi/2]
k  = [0.1,0.1]
alpha = [0.0,0.0]
beta = [0.,0.]

#Prior ranges
min_t0  = [T0[0] - 0.5, T0[1] - 0.5]
max_t0  = [T0[0] + 0.5, T0[1] + 0.5 ]
min_P   = [P[0] - 0.5,P[1] - 0.5]
max_P   = [P[0] + 0.5,P[1] + 0.5]
min_e   = [0.0,0.0]
max_e   = [1.0,1.0]
min_w   = [0.0,0.0]
max_w   = [2*np.pi,2*np.pi]
min_k   = [1e-5,1e-5]
max_k   = [1,1]
min_alpha   = [-1,-1]
max_alpha   = [1,1]
min_beta   = [-1,-1]
max_beta   = [1,1]
