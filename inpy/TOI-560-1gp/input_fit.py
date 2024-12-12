#----------------------------------------------------------
# 1D GP analysis with pyaneti
# For a multi-planet multi-instrument case
#----------------------------------------------------------

# Input RV data file
fname_rv = ['toi-560-rvs.dat']

# Number of planets in the system
nplanets = 2

#----------------------------------------------------------
# MCMC Configuration
#----------------------------------------------------------

# Thin factor for the Markov chains
thin_factor = 1

# Number of iterations per chain (effective iterations = thin_factor * niter)
niter = 500

# Number of independent Markov chains for the ensemble sampler
nchains = 100

# Select the method for analysis
# 'mcmc' -> Perform the MCMC fit
# 'plot' -> Generate plots using previously obtained results
method = 'mcmc'
#method = 'plot'

# Unit selection for planetary parameters
# Options: 'earth', 'jupiter', or 'solar'
unit_mass = 'earth'

# Flags to indicate which datasets are fitted
# True: Fit RV data; False: Do not fit
fit_rv = [True, True]
fit_tr = [False, False]

# Stellar parameters (Gaussian priors)
mstar_mean = 0.75  # Stellar mass mean in solar masses
mstar_sigma = 0.09  # Stellar mass uncertainty
rstar_mean = 0.63  # Stellar radius mean in solar radii
rstar_sigma = 0.05  # Stellar radius uncertainty
tstar_mean = 4695.  # Stellar temperature mean in K
tstar_sigma = 138.  # Stellar temperature uncertainty

# Flag for eccentricity and argument of periastron fitting (default: False)
is_ew = False

#----------------------------------------------------------
# Prior Section
#----------------------------------------------------------

# Prior types for each parameter:
# 'f' -> Fixed value
# 'u' -> Uniform prior
# 'g' -> Gaussian prior
fit_t0 = ['g', 'g']   # Mid-transit times (Gaussian priors)
fit_P = ['g', 'g']    # Orbital periods (Gaussian priors)
fit_e = ['f', 'f']    # Eccentricity (fixed)
fit_w = ['f', 'f']    # Argument of periastron (fixed)
fit_ew1 = ['u', 'u']  # Sine of eccentricity and argument of periastron (uniform)
fit_ew2 = ['u', 'u']  # Cosine of eccentricity and argument of periastron (uniform)
fit_b = ['u', 'u']    # Impact parameter (uniform)
fit_a = ['u', 'u']    # Semi-major axis (Gaussian priors, derived from stellar parameters)
fit_rp = ['u', 'u']   # Planetary radii (uniform)
fit_k = ['u', 'u']    # RV semi-amplitudes (uniform)
fit_v0 = 'u'          # Systemic velocities (uniform)
fit_q1 = ['u']        # Quadratic limb darkening parameter q1 (uniform)
fit_q2 = ['u']        # Quadratic limb darkening parameter q2 (uniform)

# Plotting flags
is_plot_correlations = False  # Do not generate correlation plots
is_plot_posterior = True      # Generate posterior distribution plots
is_plot_chains = True         # Generate chain plots

# Instrument jitter settings
is_jitter_rv = True  # Include jitter in RV fits
is_jitter_tr = True  # Include jitter in transit fits

#----------------------------------------------------------
# Prior Ranges for Parameters
#----------------------------------------------------------

# Parameter A:
# If 'f' -> Fixed to min_A
# If 'u' -> Uniform prior between min_A and max_A
# If 'g' -> Gaussian prior with mean=min_A, sigma=max_A
min_t0 = [8517.6901364, 9232.1708224]  # Mid-transit times (Gaussian prior means)
max_t0 = [0.0008259, 0.0016128]        # Mid-transit times (Gaussian prior sigmas)
min_P = [6.3980436, 18.8795529]        # Orbital periods (Gaussian prior means)
max_P = [0.0000132, 0.0008066]         # Orbital periods (Gaussian prior sigmas)
min_a = [1.1] * nplanets               # Semi-major axes (uniform prior lower bounds)
max_a = [40.] * nplanets               # Semi-major axes (uniform prior upper bounds)
min_b = [0.0] * nplanets               # Impact parameters (uniform lower bounds)
max_b = [1.0] * nplanets               # Impact parameters (uniform upper bounds)
min_e = [0.0] * nplanets               # Eccentricity (uniform lower bounds)
max_e = [1.0] * nplanets               # Eccentricity (uniform upper bounds)
min_w = [0.0] * nplanets               # Argument of periastron (uniform lower bounds)
max_w = [1.0] * nplanets               # Argument of periastron (uniform upper bounds)
min_k = [0.0] * nplanets               # RV semi-amplitudes (uniform lower bounds)
max_k = [0.05] * nplanets              # RV semi-amplitudes (uniform upper bounds)
min_rp = [0.0] * nplanets              # Planetary radii (uniform lower bounds)
max_rp = [0.05] * nplanets             # Planetary radii (uniform upper bounds)
min_ew1 = [-1.0] * nplanets            # Sine of eccentricity and argument of periastron (lower bounds)
max_ew1 = [1.0] * nplanets             # Sine of eccentricity and argument of periastron (upper bounds)
min_ew2 = [-1.0] * nplanets            # Cosine of eccentricity and argument of periastron (lower bounds)
max_ew2 = [1.0] * nplanets             # Cosine of eccentricity and argument of periastron (upper bounds)

#----------------------------------------------------------
# Telescope and Instrument Information
#----------------------------------------------------------

# Telescope identifiers used in the RV data
telescopes = ['HARPS', 'HIRES', 'PFS']

# Labels corresponding to each telescope
telescopes_labels = ['HARPS RV', 'HIRES RV', 'PFS RV']

#----------------------------------------------------------
# Gaussian Process (GP) Kernel Setup
#----------------------------------------------------------

# Kernel type for RV analysis
kernel_rv = 'QPK'  # Quasi-periodic kernel for RV analysis

# Priors for GP parameters
fit_krv = ['f'] * 4  # Initialize priors for GP parameters
fit_krv[0] = 'u'     # Uniform prior for GP parameter A_0
krv_priors = [0, 0.1]  # Prior range for A_0

fit_krv[1] = 'u'  # Uniform prior for GP parameter A_1
fit_krv[2] = 'u'  # Uniform prior for GP parameter A_2
fit_krv[3] = 'u'  # Uniform prior for GP parameter A_3

# Prior ranges for quasi-periodic kernel hyperparameters
QP_priors = [1., 100.0, 0.1, 2.0, 10, 14]  # Hyperparameters: Amplitude, lambda_e, lambda_p, etc.

# Combine GP parameter priors
krv_priors = np.concatenate([krv_priors, QP_priors])

#----------------------------------------------------------
# END
#----------------------------------------------------------

