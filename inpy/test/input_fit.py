#Input file for test
#Mon 27, Jun, 2016
#Barragan, O.

#telescopes labels
telescopes = ['S']
telescopes_labels = ['Super_telescope']

#Input files
fname_rv = 'earth_rv.dat'
fname_tr = 'earth_lc.dat'

#chain parameters
thin_factor = 1
nconv = 500
nwalkers = 100

#photometric data
lc_data = 'free'
n_cad = 1
t_cad = 30. / 60. / 24.

method = 'sm'
method = 'plot'

mstar_mean = 1.0
mstar_sigma = 0.1
rstar_mean = 1.0
rstar_sigma = 0.1

#Output controls
unit_mass = 'earth'

is_ew = False
#Star parameters
a_from_kepler = True

#We want to fit transit and RV 
#fit_rv = True
fit_tr = True

#What do we want to fit?
fit_t0 = True
fit_P  = True
fit_e  = False
fit_w  = False
fit_i  = True
fit_a  = True
fit_q1 = True
fit_q2 = True
fit_pz = True
fit_k  = True
fit_v0 = True

#Priors
T0 = 2448285.2569
P  = 365.256
e  = 1e-5
w  = np.pi/2
k  = 0.001
ii = 0.999999*np.pi/2
q1 = 1e-5
q2 = 1e-5

#Prior ranges
min_t0  = 2448285.05  
max_t0  = 2448285.15  
min_P   = P - 0.05
max_P   = P + 0.05
min_e   = 0.0
max_e   = 1.0
min_w   = 0.0
max_w   = 2*np.pi
min_i   = 85*np.pi/180
max_i   = np.pi/2
min_a   = 200
min_phys_a = 200
max_a   = 250
max_phys_a = 250
min_k   = 1e-5
min_pz  = 0.0085
max_pz  = 0.0095
