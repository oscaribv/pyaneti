#-----------------------------------------------------------
#                       default.py
#  This file contains the defalut values to run pyaneti.
#  You can change all the values here. You can also control
#	 the values inside the input_file.py
#                 O. Barragan, March 2016
#-----------------------------------------------------------

#Constants
#Solar and planetary constants according to IAU 2015
# http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1605.09788
#Sun
S_radius_SI    = 6.957e8       	  #m
S_radius_cgs   = S_radius_SI*1e2  #cm
S_GM_SI        = 1.3271244e20     #m^3 s^{-1} 
S_GM_cgs       = S_GM_SI*1e6      #cm^3 s^{-1} 
G_SI           = 6.67408e-11      #m^3 kg^{-1} s^{-2}
G_cgs          = G_SI*1e3         #cm^3 g^{-1} s^{-2}
S_vol_SI       = 4./3.*np.pi*S_radius_SI**3  #m^3 
S_vol_cgs      = 4./3.*np.pi*S_radius_cgs**3 #cm^3
S_den_SI       = S_GM_SI/G_SI/S_vol_SI    #kg/m^3
S_den_cgs      = S_GM_cgs/G_cgs/S_vol_cgs #g/cm^3
S_Teff         = 5772.0 # K
#Earth
E_GM_SI        = 3.986004e14      #m^3 s^{-2}
E_radius_e_SI  = 6.3781e6         # ecuatorial radius [m]
E_radius_p_SI  = 6.3568e6         # polar radius [m]
#Jupiter
J_GM_SI        = 1.2668653e17     #m^3 s^{-1} 
J_radius_e_SI  = 7.1492e7         # ecuatorial radius [m]
J_radius_p_SI  = 6.6854e7         # polar radius [m]
#Other constants
AU_SI = 1.4960e11 # m

#Default stellar parameters
mstar_mean = 1.0
rstar_mean = 1.0
tstar_mean = S_Teff
mstar_sigma = 1.e-1
rstar_sigma = 1.e-1
tstar_sigma = 100.

#Default transit data cadence
n_cad = 1
t_cad = 30. / 60. / 24.
#Default input columns for transit fit
columns_tr = [0,1,2]

#No telescopes
telescopes = []
telescopes_labels = ['']

#rv jitter controls
is_special_jitter = False
jrvvec = []

#MCMC controls
method = 'mcmc'
is_circular = False
maxi = int(1e8)
thin_factor = 10
niter   = 100
nchains = 100
n_transits = 1
nt = 1
errores='perc'
unit_mass = 'solar'
scale_error_bars = False
s_factor = 1.0
lc_data = 'free'
a_from_kepler = [False]
a_factor = 2.0
is_plot_posterior = True
is_plot_correlations = True
is_plot_chains = True
is_corner_plot = False
is_smart_priors = True
is_print_mode = False
is_clustering = True
latex_values = True
is_plot_all_tr = [False]
is_special_plot_tr = False #Plot percentiles in the phase pholded transit plots
is_special_plot_rv = False #Plot percentiles in the phase pholded RV plots
is_fit_tr_errorbars = True
plot_unbinned_model = False #If True, overplots an unbinned model in the light curve
is_jitter_rv = False
is_jitter_tr = False
plot_parameters = []
is_seaborn_plot = False
seaborn_palette = 'deep'
plot_tr_errorbars = False
is_plot_std_tr = False
rv_markersize  = 8
rv_fillstyle = 'full'
outdir = 'outpy/'
resize_rv = False
resize_tr = True
span_tr  = [0.0]
is_TTV = True
textra = 0.0
my_tr_err = 0
select_y_tr = False
y_lim_max = 1.05
y_lim_min = 0.95
select_y_rv = False
rv_lim_max = -100 #m/s
rv_lim_min = 100  #m/s
#The defaul number of planets is 1
nplanets = 1
#The default planet labels
plabels = ['b','c','d','e','f','g','h','i']
#Figure size default
figure_size_x = 23./2.56 #25cm
figure_size_y = 23./2.56/1.618
font_size_label = 18
clustering_delta = 0.5
get_value = 'median'
sigma_clean = 10
#mark for plots
rv_xlabel = "BJD - 2450000 (days)"
tr_xlabel = "BJD - 2450000 (days)"
is_rv_legend = True
mark = ['o', 'D', 's', 'p', 'h', '8', '^', '<', '*', 'v','>','.', 'H', 'd','+']
rv_colors = ['b','r','g','y','m','c']
tr_colors = 'r'

#Fit nothing
fit_rv = False
fit_tr = False

#All the parameters are fitted
#by default
fit_t0  = ['f']
fit_P   = ['f']
fit_e   = ['f']
fit_w   = ['f']
fit_ew1 = ['f']
fit_ew2 = ['f']
#transit fit
fit_i   = ['f']
fit_b   = ['f']
fit_a   = ['f']
fit_q1  = 'f'
fit_q2  = 'f'
fit_rp  = ['f']
#rv fit
fit_k   = ['f']
fit_v0  = ['f']

is_linear_trend = 'f'
is_quadratic_trend = 'f'

#flags
is_log_P     = False #Allows sampling for the log10 of the period
is_ew        = True  #allows sampling for sqrt(e)sinw and sqrt(e)cosw
is_b_factor  = True  #allows sampling for b instead of i
is_den_a     = False #allows sampling for stellar density instead of a
is_log_k     = False #allows sampling for log of K
is_log_rv0   = False

#Default priors rages (wide)
min_t0  = [0.0]          #days
max_t0  = [1e6]         #days
min_P   = [0.1]       #days
max_P   = [1e4]      #days
min_e   = [1.e-10]      #zero
max_e   = [0.999]       #one
min_w   = [0.0]      #rad
max_w   = [2.0*np.pi]      #rad
min_ew1  = [-1.0]
max_ew1  = [1.0]
min_ew2  = [-1.0]
max_ew2  = [1.0]
#transit fit
min_i   = [0.0 ]    # 70 degrees
max_i   = [np.pi/2.0]  # 90 degrees
min_b   = [0.0]
max_b   = [1.0]
min_a   = [1.5]       # The planet is outside the star
max_a   = [1.e3]       # The planet is really far
min_q1  = 0.0   #
max_q1  = 1.0     #
min_q2  = 0.0     #
max_q2  = 1.0      #
min_rp  = [0.0]       # Earth size planet / sun
max_rp  = [0.99]      # a really big planet
#rv fit
min_k   = [1.e-6]      # m/s amplitudes
max_k   = [30.]     # a really big planet
