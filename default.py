#-----------------------------------------------------------
#                       default.py
#  This file contains the defalut values to run pyaneti.
#  You can change all the values here. You can also control
#	 the values in input_file.py
#                 O. Barragan, March 2016
#-----------------------------------------------------------

#Plot parameters
units_ms = False
ylab = 'RV (km/s)'
ktom = 1.0
mstar= 1.0

#No telescopes
telescopes = []

#MCMC controls
method = 'sm'
is_circular = False
ics = False
prec = 1.e-4
maxi = int(1e8)
thin_factor = int(5e3)
nconv = 100
errores='perc'
nwalkers = 100
unit_mass = 'solar'
scale_error_bars = True
s_factor = 1.0

#Default priors
P = 15.
e = 0.1
w = 3.1415926 / 2.0
ii = 3.1415026 / 2.0
a = 13.0
u1 = 0.42
u2 = 0.25
pz = 0.5


#Transit fit
nbin = 16

#The defaul number of planets is 1
nplanets = 1

#Fit nothing
fit_rv = False
fit_tr = False

#All the parameters are fitted
#by default
fit_t0  = True
fit_P   = True
fit_e   = True
fit_w   = True
#transit fit
fit_i   = True
fit_a   = True
fit_u1  = True
fit_u2  = True
fit_pz  = True
#rv fit
fit_k   = True
fit_v0  = True


#flags
is_log_P     = True
is_ew        = True
is_sini      = True
is_log_a     = True
is_log_k     = True
is_log_rv0   = True

#Limits
min_t0  = 0.0          #days
max_t0  = 1e6          #days
min_P   = 0.1 	       #days
max_P   = 1e4 	       #days
min_e   = 1.e-10       #zero
max_e   = 0.999	       #one
min_w   = 0.0	       #rad
max_w   = 2*np.pi      #rad
#transit fit
min_i   = 1.22173      # 70 degrees
max_i   = np.pi / 2.0  # 90 degrees
min_a   = 1.5	       # The planet is outside the star
max_a   = 1.e8	       # The planet is really far
min_u1  = 0.0	       #
max_u1  = 0.5          #
min_u2  = 0.0	       #
max_u2  = 0.5	       #
min_pz  = 1.e-2	       # Earth size planet / sun
max_pz  = 0.99	       # a really big planet
#rv fit
min_k   = 1.e-3	       # m/s amplitudes
max_k   = 30.	       # a really big planet
min_rv0 = 1.	       #Systemic velocities
max_rv0 = 100.	       #systemic velocities

