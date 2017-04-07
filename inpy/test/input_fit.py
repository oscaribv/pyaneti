#Input file for test
#Barragan, O.

#telescopes labels
telescopes = ['S']
telescopes_labels = ['Super_telescope']

#Input files
fname_rv = ['earth_rv.dat']
fname_tr = ['earth_lc.dat']

#chain parameters
thin_factor = 1
nconv = 500
nwalkers = 100

#photometric data
lc_data = 'free'
n_cad = 1
t_cad = 30. / 60. / 24.

method = 'mcmc'
#method = 'plot'
is_seaborn_plot = True

mstar_mean = 1.0
mstar_sigma = 0.1
rstar_mean = 1.0
rstar_sigma = 0.1

#Output controls
unit_mass = 'earth'

is_plot_histogram = False
is_plot_correlations = False
is_plot_chains = False

#Star parameters
a_from_kepler = [True]

#We want to fit transit and RV 
fit_rv = [True]
fit_tr = [True]

#What do we want to fit?
fit_t0 = ['u']
fit_P  = ['u']
fit_e  = ['f']
fit_ew1= ['f']
fit_ew2= ['f']
fit_w  = ['f']
fit_b  = ['f']
fit_a  = ['u']
fit_q1 = 'g'
fit_q2 = 'g'
fit_rp = ['u']
fit_k  = ['u']
fit_v0 = 'u'

#Priors
T0 = [2448285.2569]
P  = [365.256]
e  = [1e-5]
w  = [np.pi/2]
k  = [0.001]
ii = [0.999999*np.pi/2]
q1 = [1e-5]
q2 = [1e-5]

#is_ew = False

#Prior ranges
min_t0  = [2448285.05] 
max_t0  = [2448285.15]  
min_P   = [P[0] - 0.05]
max_P   = [P[0] + 0.05]
min_e   = [0.0]
min_w   = [0.0]
min_ew1 = [0.0]
min_ew2 = [0.0]
min_i   = [0.0]
max_i   = [1.0]
min_a   = [200]
max_a   = [250]
min_k   = [0.0]
max_k   = [0.001]
min_rp  = [0.0]
max_rp  = [0.1]
min_q1  = 0.3464
max_q1  = 0.05
min_q2  = 0.2839
max_q2  = 0.05
