#Plot parameters
ylab = 'RV (km/s)'
units_ms = False
ktom = 1.0
mstar= 1.0

#MCMC controls
is_circular = False
tfactor = 10000
ics = False
prec = 1.e-4
maxi = int(1e8)
thin_factor = int(5e3)
nconv = 100

#Default priors
P = 15.
e = 0.1
w = 3.1415926 / 2.0
ii = 3.1415026 / 2
a = 13.0
u1 = 0.42
u2 = 0.25
pz = 0.5

#All the parameters are fitted
#by default
fit_t0  = True
fit_P   = True
fit_e   = True
fit_w   = True
fit_i   = True
fit_a   = True
fit_u1  = True
fit_u2  = True
fit_pz  = True
fit_k   = True
fit_v0  = True

#Fit nothing 
fit_rv = False
fit_tr = False

