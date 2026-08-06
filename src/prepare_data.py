# -----------------------------------------------------------
#                       prepare_data.py
#  This file contains all the variable initializations,
#  both for RV and Transit fittings.
#                   O. Barragan, March 2016
# -----------------------------------------------------------
import numpy as np
import pandas as pd
import pickle

nconv = niter
nwalkers = nchains

# -----------------------------------------------------------
#                         RV DATA
# -----------------------------------------------------------

# Let us check the kind of variable
nplanets_rv = 0
for o in range(0, len(fit_rv)):
    nplanets_rv += int(fit_rv[o])

#Pyaneti will do FCO method automatically if is_run_fco = True
if nplanets_rv > 0:
    if is_run_fco:
        file_path = f'inpy/{star}/{fname_rv[0]}'
        out_file_path = f'inpy/{star}/fco-{fname_rv[0]}'
        create_fco_file(file_path,out_file_path,fco_window)

if nplanets_rv > 0:
    # Read the data file and determine the number of columns
    file_path = f'inpy/{star}/{fname_rv[0]}'
    if is_run_fco: file_path = out_file_path
    with open(file_path, 'r') as file:
    # Loop to skip comment lines that start with #
        while True:
            first_line = file.readline().strip()
            if not first_line:  # End of file with no non-comment lines
                raise ValueError("The file contains only comments or is empty.")
            if not first_line.startswith('#'):
                break

    # Split the first non-comment line and count the number of columns
    n_cols = len(first_line.split())

    # Define the column names and usecols based on the number of columns
    col_names = ['time', 'rv', 'err', 'telescope', 'sjitter', 'dim_lab']
    usecols = list(range(min(n_cols, len(col_names))))

    # Load the data with the appropriate columns
    data = pd.read_csv(file_path,sep='\s+', comment='#', header=None,
                       names=col_names[:n_cols], usecols=usecols, engine='python')


    # Fill missing columns with default values
    if 'telescope' not in data.columns:
        data['telescope'] = '0'
    if 'sjitter' not in data.columns:
        data['sjitter'] = None
    if 'dim_lab' not in data.columns:
        data['dim_lab'] = None

    # Initialize lists for storing data for each telescope/offset
    time_all = []
    rv_all = []
    errs_all = []
    sjitter_all = []
    dim_lab_all = []

    if is_run_fco:
        jrvvec = list(set(data.sjitter))
        telescopes = list(set(data.telescope))
        telescopes = list(map(str,telescopes))
        is_special_jitter = True
        telescopes_labels_jitter = list(telescopes_labels)
        telescopes_labels = ["offset_" + x for x in telescopes]
        #Pyaneti has a limit of 10 default colours and symbols
        #Since now we are using more than 10 offsets (instruments), then pyaneti will not be able to deal with the colours
        #and symbols automatically, but we can tell to pyaneti the symbols and colours that we want
        #In this example, all our symbols will be 'o' and the colours will be generated from random hexagesimal values
        #The following line tells to pyaneti to create a mark for each offset as 'o', by assingning a list to the variable mark
        mark = ['o']*len(telescopes)
        #The next lines generate random hexagesimal values that are used as colour codes for matplotlib
        #the vector that passes the colours to pyaneti is rv_colors
        from random import randint, seed
        rv_colors = [] #Initiante rv_colors as an empty list
        seed(1)        #Give a seed to the random number generator
        #Add a different colour, element by element to rv_colors
        for i in range(len(telescopes)):
            rv_colors.append('#%06X' % randint(0, 0xFFFFFF))
        #Let us remove the legend that indicates each offset from the RV plots
        is_rv_legend = False
        #We search for only one planet in a FCO run
        nplanets = 1


    # Number of telescopes/offsets
    nt = len(telescopes)
    #Let's be sure that the telescope data are strings, they can be given as numbers sometimes
    data.telescope = data.telescope.astype(str)


    # Check if telescope labels are provided
    if not telescopes:
        print('Please, indicate the telescope labels!')
        sys.exit('')


    # Separate the data for each telescope
    for tel in telescopes:
        mask = (data['telescope'] == tel)
        time_all.append(data['time'][mask].tolist())
        rv_all.append(data['rv'][mask].tolist())
        errs_all.append(data['err'][mask].tolist())
        sjitter_all.append(data['sjitter'][mask].tolist())
        dim_lab_all.append(data['dim_lab'][mask].tolist())

    # Create mega variables containing all the data
    rv_vals = []
    rv_time = []
    rv_errs = []
    sjitter_vals = []
    dim_lab_vals = []
    tlab = []

    # Aggregate data from all telescopes
    for i, telescope in enumerate(telescopes):
        rv_vals.extend(rv_all[i])
        rv_time.extend(time_all[i])
        rv_errs.extend(errs_all[i])
        sjitter_vals.extend(sjitter_all[i])
        dim_lab_vals.extend(dim_lab_all[i])
        tlab.extend([i] * len(rv_all[i]))

    tlab = np.array(tlab)


    # If special jitter handling is required
    if is_special_jitter:
        rv_jitter_labels = list(set(sjitter_vals))
        n_jrv = len(jrvvec)
        jrvlab = np.zeros(shape=len(tlab))
        seen = set()
        unique_elements = []
        for element in sjitter_vals:
            if element not in seen:
                seen.add(element)
                unique_elements.append(element)

        for i in range(len(jrvlab)):
            for j in range(len(unique_elements)):
                if (sjitter_vals[i] == unique_elements[j]):
                    jrvlab[i] = j
                    break
    else:
        rv_jitter_labels = telescopes_labels
        jrvlab = tlab
        n_jrv = len(set(tlab))


    # Handle multidimensional Gaussian Processes dimension labels if not all data are simultaneous
    if kernel_rv[0] == 'M':
        if data['dim_lab'][0] != None: #if we gave the dimension labels in the input data file
            dim_lab = np.zeros(shape=len(tlab))
            ndim = int(kernel_rv[2])
            seen = set()
            unique_elements = []
            for element in dim_lab_vals:
                if element not in seen:
                    seen.add(element)
                    unique_elements.append(element)
            for i in range(len(dim_lab)):
                for j in range(ndim):
                    if (dim_lab_vals[i] == unique_elements[j]):
                        dim_lab[i] = j
                        break
        else: #In this case we assume all data are simultaneous
            ndim = int(kernel_rv[2])
            nobs_per_d = len(rv_time) // ndim
            dim_lab = np.concatenate([[o] * nobs_per_d for o in range(ndim)])
    else:
         dim_lab = np.zeros(shape=len(tlab))

    dim_lab = np.array(dim_lab)
    total_rv_fit = True

else:
    # If no planets to fit, initialize variables with default values
    n_jrv = 1
    nt = 1
    tlab = [0]
    jrvlab = tlab
    dim_lab = tlab
    rv_jitter_labels = tlab
    rv_vals = [1]
    rv_time = [1]
    rv_errs = [1]
    total_rv_fit = False
    is_jitter_rv = False



# -----------------------------------------------------------
# RV DATA READY
# -----------------------------------------------------------

# -----------------------------------------------------------
#                     TRANSIT DATA
# -----------------------------------------------------------

# Let us check the kind of variable
nplanets_tr = 0
for o in range(len(fit_tr)):
    nplanets_tr += int(fit_tr[o])

if (nplanets_tr > 0):

    file_path = f'inpy/{star}/{fname_tr[0]}'
    if (my_tr_err == 0):  # The error bars come from the input file
        lc_time, lc_flux, lc_errs = np.loadtxt(file_path, usecols=[0,1,2],
                                               comments='#', unpack=True)
    else:  # The error bars are given in the input
        lc_time, lc_flux = np.loadtxt(file_path, usecols=[0, 1],
                                       comments='#', unpack=True)
        lc_errs = np.array([my_tr_err]*len(lc_time))

    #Add time in order to match with the spectroscopic time-series
    lc_time += textra

    # Create the label vectors for each instrument and jitter
    if (bands[0] == None):
        bands[0] = 'lc'
        nbands = 1
        nradius = 1
        trlab = [0]*len(lc_time)
        jtrlab = [0]*len(lc_time)
        try:
            n_cad = n_cad[:1]
            t_cad = t_cad[:1]
        except:
            n_cad = [n_cad]
            t_cad = [t_cad]
    else:
        trlab = []
        jtrlab = []
        nbands = len(bands)
        nradius = 1
        if is_multi_radius:
            nradius = nbands
        instrument = np.loadtxt(file_path, usecols=[3], dtype=str, unpack=True)
        good_index = []
        for o in range(len(lc_time)):
            for m in range(nbands):
                if (instrument[o] == bands[m]):
                    trlab.append(m)
                    jtrlab.append(m)
                    good_index.append(o)

        #Now let us store the data with the bands that we indicated in the input file
        lc_time = lc_time[good_index]
        lc_flux = lc_flux[good_index]
        lc_errs = lc_errs[good_index]
        n_cad = n_cad[:nbands]
        t_cad = t_cad[:nbands]

    n_jtr = nbands

    total_tr_fit = True

    if bin_lc > 0:
        lc_time, lc_flux, lc_errs = bin_data(lc_time, lc_flux, lc_errs, bin_lc)
        nbands = 1
        nradius = 1
        trlab = [0]*len(lc_time)
        jtrlab = [0]*len(lc_time)

else:
    lc_time = [1.]
    lc_flux = [1.]
    lc_errs = [1.]
    total_tr_fit = False
    is_jitter_tr = False
    fit_q1 = 'f'
    fit_q2 = 'f'
    bands[0] = ' '

# TRANSIT DATA READY

# CHECK WHAT WE HAVE TO FIT
# If we are not going to fit RV or TR data, let us turn off the variables
# for the given case
for o in range(0, nplanets):
    if (fit_tr[o] == False):
        fit_rp[o] = 'f'
        min_rp[o] = 0.0
        if (not is_b_factor):
            fit_i[o] = 'f'
            min_i[o] = np.pi/2.0
        else:
            fit_b[o] = 'f'
            min_b[o] = 0.0
        fit_a[o] = 'f'
    if (fit_rv[o] == False):
        fit_k[o] = 'f'

if is_single_transit:
    for o in range(0, nplanets):
        fit_P[o] = 'f'
        min_P[o] = 10*(max(lc_time)-min(lc_time))
        fit_a[o] = 'u'
        min_a[o] = 1.1
        max_a[o] = 1e3
        plot_binned_data = False

# Let us turn off velocity offset for a pure TR fit
if (not total_rv_fit):
    fit_v0 = 'f'
    nt = 1
    min_rv0 = [0.0]
    max_rv0 = [1.]
    telescopes = ['telescope']
    telescopes_labels = ['telescope']

# This ensures that previous 1-band pyaneti input files work
if (min_q1.__class__ != list):
    min_q1 = [min_q1]
if (min_q2.__class__ != list):
    min_q2 = [min_q2]
if (max_q1.__class__ != list):
    max_q1 = [max_q1]
if (max_q2.__class__ != list):
    max_q2 = [max_q2]
if (fit_q1.__class__ == str):
    fit_q1 = [fit_q1]
if (fit_q2.__class__ == str):
    fit_q2 = [fit_q2]

#Let us make the code complatible with the old pyaneti files that include is_den_a
if "is_den_a" in locals():
    sample_stellar_density = is_den_a

if sample_stellar_density:  # For a multiplanet system the density has to be the same
        for o in range(1, nplanets):
            fit_a[o] = 'f'

#Check if the variable unit_mass exists to make it compatible with older versions of pyaneti
#unit_mass was replaced by planetary_units
if "unit_mass" in locals(): planetary_units = unit_mass

# -----------------------------------------------------------
#  Smart priors, get the best values of the physical and
#  priors limits
# -----------------------------------------------------------

# Let us try to do a guess for the init values
if (total_rv_fit):

    # Estimate systemic velocity priors and limits from data
    # The systemic velocity value of all the telescope should
    # be between the smallest and larger RV datapoint
    if fit_v0.__class__ == str:
        min_rv0 = [None]*nt
        max_rv0 = [None]*nt
        fit_v0_val = fit_v0
        fit_v0 = [None]*nt
        for o in range(nt):
            if (fit_v0_val == 'u'):
                fit_v0[o] = 'u'
                min_rv0[o] = min(rv_all[o]) - 0.1
                max_rv0[o] = max(rv_all[o]) + 0.1
                if is_linear_trend == 'u':
                    min_rv0[o] = min(rv_all[o]) - 1.0
                    max_rv0[o] = max(rv_all[o]) + 1.0
            else:
                fit_v0[o] = 'f'
                min_rv0[o] = 0.0
                max_rv0[o] = 0.0
    else:
        if (len(fit_v0) != nt):
            print('fit_v0 in your input file does not match the number of telescopes')
            sys.exit()

if (is_ew):
        min_e = min_ew1
        max_e = max_ew1
        min_w = min_ew2
        max_w = max_ew2
        fit_e = fit_ew1
        fit_w = fit_ew2

if (is_b_factor):
        min_i = min_b
        max_i = max_b
        fit_i = fit_b

# Let us check what do we want to fit
total_fit_flag = [total_rv_fit, total_tr_fit]

flags = [is_log_P, is_ew, is_b_factor, sample_stellar_density, is_log_k,is_log_rv0]

# Parameter priors
pars_prior_flag = [None]*7*nplanets
for o in range(0, nplanets):
        pars_prior_flag[o*7:(o+1)*7] = [fit_t0[o], fit_P[o], fit_e[o], fit_w[o],
                                        fit_i[o], fit_a[o], fit_k[o]]

pars_prior_vals = [None]*7*2*nplanets
for o in range(0, nplanets):
        pars_prior_vals[o*7*2:(o+1)*7*2] = \
            [min_t0[o], max_t0[o], min_P[o], max_P[o], min_e[o], max_e[o], min_w[o],
                max_w[o], min_i[o], max_i[o], min_a[o], max_a[o], min_k[o], max_k[o]]

# Planet radius priors
prs_prior_flag = [None]*nplanets
prs_prior_values = [None]*nplanets
for o in range(0, nplanets):
        prs_prior_flag[o] = [fit_rp[o]]*nradius
        prs_prior_values[o] = [min_rp[o], max_rp[o]]*nradius

prs_prior_flag = np.concatenate(prs_prior_flag)
prs_prior_values = np.concatenate(prs_prior_values)

# LDC priors
ldc_prior_flag = []
ldc_prior_values = []
for o in range(nbands):
        ldc_prior_flag.append(fit_q1[o])
        ldc_prior_flag.append(fit_q2[o])
        ldc_prior_values.append(min_q1[o])
        ldc_prior_values.append(max_q1[o])
        ldc_prior_values.append(min_q2[o])
        ldc_prior_values.append(max_q2[o])

# Offsets priors
RVS_prior_flag = []
RVS_prior_vals = []
for o in range(0, nt):
        RVS_prior_flag.append(fit_v0[o])
        RVS_prior_vals.append(min_rv0[o])
        RVS_prior_vals.append(max_rv0[o])

# RV jitter priors
if is_jitter_rv:
        if jrv_prior_flag == None: #We set authomatic priors
            jrv_prior_flag = ['u']*n_jrv
            jrv_prior_vals = [0,1e-1]*n_jrv
        else: # We gave the priors in the input file
            if (len(jrv_prior_flag) != n_jrv):
                print('jrv_prior_flag in your input file does not match the number of jitter terms to sample!')
                sys.exit()
            if (len(jrv_prior_vals) != 2*n_jrv):
                print('jrv_prior_vals in your input file does not contain the expected number of prior values!')
                sys.exit()
else:
        jrv_prior_flag = ['f']*n_jrv
        jrv_prior_vals = [0., 0.5]*n_jrv

    # Transit jitter priors
if is_jitter_tr:
        if jtr_prior_flag == None: #We set authomatic priors
            jtr_prior_flag = ['u']*n_jtr
            jtr_prior_vals = [0,1e-2]*n_jtr
        else: # We gave the priors in the input file
            if (len(jtr_prior_flag) != n_jtr):
                print('jtr_prior_flag in your input file does not match the number of jitter terms to sample!')
                sys.exit()
            if (len(jtr_prior_vals) != 2*n_jtr):
                print('jtr_prior_vals in your input file does not contain the expected number of prior values!')
                sys.exit()
else:
        n_jtr = 0
        trlab = [0]*len(lc_time)
        jtrlab = [0]*len(lc_time)
        jtr_prior_flag = ['f']*n_jtr
        jtr_prior_vals = [0., 1.e-3]*n_jtr

   # Trends priors
trends_prior_flag = [is_linear_trend, is_quadratic_trend]
val_linear = 0.1
val_quadratic = 0.1
if (is_linear_trend == 'f'):
        val_linear = 0.0
if (is_quadratic_trend == 'f'):
        val_quadratic = 0.0
trends_prior_vals = [-val_linear, val_linear, -val_quadratic, val_quadratic]

# RV Kernels
if kernel_rv == 'None':
    krv_prior_flag = []
    krv_prior_vals = []
    krv_labels = []
else:
    krv_prior_flag = fit_krv
    krv_prior_vals = krv_priors

    if kernel_rv == 'QP2':
        krv_labels = ['A', 'Gamma', 'metric', 'P_GP']
    elif kernel_rv in ['M32', 'M52']:
        krv_labels = ['A', 'metric']
    elif kernel_rv == 'QPK':
        krv_labels = ['A', 'le', 'lp', 'P_GP']
    elif kernel_rv[0:2] in ['MQ', 'SQ', 'ME', 'MM']:
        krv_labels = [f'A{m}' for m in range(len(fit_krv) - (1 if kernel_rv[0:2] in ['ME', 'MM'] else 3))]

        if kernel_rv[0:2] == 'MQ':
            krv_labels += ['lambdae', 'lambdap', 'PGP']
        elif kernel_rv[0:2] == 'SQ':
            krv_labels += ['lambdae', 'lambdap1', 'PGP1', 'lambdap2', 'PGP2']
        elif kernel_rv[0:2] in ['ME', 'MM']:
            krv_labels += ['lambdae']

#If we give a list with labels to be used, we can use it instead of the authomatically created label
if len(gprv_labels) < 2 or len(gprv_labels) != len(krv_labels):
    gprv_labels = krv_labels

np_rv = len(krv_prior_flag)

# TR Kernels
if kernel_tr == 'None':
    ktr_prior_flag = []
    ktr_prior_vals = []
    ktr_labels = []
else:
    ktr_prior_flag = fit_ktr
    ktr_prior_vals = ktr_priors

    if kernel_tr == 'QPK':
        ktr_labels = ['A', '$\lambda_e$', '$\lambda_p$', 'P']
    elif kernel_tr == 'QP2':
        ktr_labels = ['A', 'metric', '$\Gamma$', 'P']
    elif kernel_tr == 'SEK':
        ktr_labels = ['A', '$\lambda$']


np_tr = len(ktr_prior_flag)

kernels = kernel_rv[0:3]+kernel_tr[0:3]

prior_vals = np.concatenate([pars_prior_vals, prs_prior_values, ldc_prior_values,
                                 RVS_prior_vals, jrv_prior_vals, jtr_prior_vals, trends_prior_vals,
                                 krv_prior_vals, ktr_prior_vals])
prior_flags = np.concatenate([pars_prior_flag, prs_prior_flag, ldc_prior_flag,
                                  RVS_prior_flag, jrv_prior_flag, jtr_prior_flag, trends_prior_flag,
                                  krv_prior_flag, ktr_prior_flag])

model_int = [nplanets, nt, nbands, nradius,
                 nldc, n_jrv, n_jtr, np_rv, np_tr]
model_int = np.concatenate([model_int, n_cad])
model_double = t_cad

# Define the labels of the sampled parameters
labs = []
elab = 'e'
wlab = 'omega'
ilab = 'i'
klab = 'k'
alab = 'arstar'
if (is_ew):
    elab = 'esinomega'
    wlab = 'ecosomega'
if (is_b_factor):
    ilab = 'b'
if (is_log_k):
    klab = 'log10k'
if (sample_stellar_density):
    alab = 'rhostar'
# planet parameter labels
for o in range(nplanets):
    etiquetas = ['t0'+plabels[o], 'p'+plabels[o], elab+plabels[o],
                    wlab+plabels[o], ilab+plabels[o], alab+plabels[o],
                    klab+plabels[o]]
    labs.append(etiquetas)
# planet radius labels
for o in range(nplanets):
    for m in range(nradius):
        if nradius == 1:
            labs.append(['rprstar'+plabels[o]])
        else:
            labs.append(['rprstar'+bands[m]+plabels[o]])
# LDC labels
if nbands == 1:
    labs.append(['q1', 'q2'])
else:
    for m in range(nbands):
        labs.append(['q1'+bands[m], 'q2'+bands[m]])
    # RV instrument labels
labs.append(telescopes_labels)
# jitter labels
rv_jitter_names = []
for o in range(n_jrv):
    labs.append(['rv_jitter'+str(rv_jitter_labels[o])])
    rv_jitter_names.append(['rv_jitter'+str(rv_jitter_labels[o])])
rv_jitter_names = list(np.concatenate(rv_jitter_names))
for o in range(n_jtr):
    labs.append(['tr_jitter'+str(bands[o])])
# trends labels
labs.append(['linear_trend'])
labs.append(['quadratic_trend'])
labs.append(krv_labels)
labs.append(ktr_labels)
labs = np.concatenate(labs)
param_names = list(labs)
# Total labels vector)

# Initialize the unit labels for parameters
units = []

elab_unit = ''  # Eccentricity has no unit
wlab_unit = 'radians'  # Argument of periastron in radians
ilab_unit = 'radians'  # Inclination in radians
klab_unit = 'km/s'  # RV semi-amplitude in km/s
alab_unit = ''  # Scaled semi-major axis, no unit

if is_ew:
    elab_unit = ''  # Eccentricity components have no unit
    wlab_unit = ''  # No unit for eccentricity components
if is_b_factor:
    ilab_unit = ''  # Impact parameter, no unit
if is_log_k:
    klab_unit = ''  # Logarithmic RV semi-amplitude, no unit
if sample_stellar_density:
    alab_unit = 'g/cm^3'  # Stellar density in g/cm^3

# Planet parameter units
for o in range(nplanets):
    etiquetas_units = ['days', 'days', elab_unit, wlab_unit, ilab_unit, alab_unit, klab_unit]
    units.append(etiquetas_units)

# Planet radius units
for o in range(nplanets):
    for m in range(nradius):
        if nradius == 1:
            units.append([''])  # Planet-to-star radius ratio, no unit
        else:
            units.append([''])  # Planet-to-star radius ratio for different bands, no unit

# LDC units
for m in range(nbands):
    if nradius == 1:
        units.append(['', ''])  # Limb darkening coefficients, no units
    else:
        units.append(['', ''])  # Limb darkening coefficients for different bands, no units

# RV instrument units
units.append(['km/s'] * len(telescopes_labels))  # Radial velocity measurements in km/s

# Jitter units
for o in range(n_jrv):
    units.append(['km/s'])  # RV jitter in km/s
for o in range(n_jtr):
    units.append([''])  # Transit jitter, no units

# Trends units
units.append(['km/s/day'])  # Linear trend in RV measurements
units.append(['km/s/day^2'])  # Quadratic trend in RV measurements
#Kernel units
units.append(krv_labels)
units.append(ktr_labels)
#
units = np.concatenate(units)
sampled_param_units = list(units)
# param_units now contains all the units corresponding to the parameters



#Save the pickle file that needs to be used in likelihood calls
# Collect all variables into a dictionary
all_data_dict = {
    'rv_time': rv_time,
    'rv_vals': rv_vals,
    'rv_errs': rv_errs,
    'lc_time': lc_time,
    'lc_flux': lc_flux,
    'lc_errs': lc_errs,
    'prior_flags': prior_flags,
    'prior_vals': prior_vals,
    'tlab': tlab,
    'jrvlab': jrvlab,
    'dim_lab': dim_lab,
    'trlab': trlab,
    'jtrlab': jtrlab,
    'total_fit_flag': total_fit_flag,
    'flags': flags,
    'kernels': kernels,
    'model_int': model_int,
    'model_double': model_double,
    'param_names': param_names
}

# Save dictionary to a pickle file
fpickle = f'{outdir}/posterior_parameters.pkl'
with open(fpickle, 'wb') as pickle_file:
    pickle.dump(all_data_dict, pickle_file)


# -----------------------------------------------------------
#          PRINT INITIAL CONFIGURATION
# -----------------------------------------------------------



if True:
    out_init_file = outdir+'/'+star+'_init.dat'
    oif = open(out_init_file, 'w')
    oif.write('\n')
    oif.write('==============================\n')
    oif.write('------------------------------\n')
    oif.write("    INITIAL CONFIGURATION     \n")
    oif.write('------------------------------\n')
    oif.write('Star           = %s\n' % star)
    oif.write('No. planets    = %d\n' % nplanets)
    oif.write('------------------------------\n')
    oif.write('iter max       = %d\n' % maxi)
    oif.write('thin factor    = %d\n' % thin_factor)
    oif.write('nconv          = %d\n' % nconv)
    oif.write('nwalkers       = %d\n' % nwalkers)
    oif.write('------------------------------\n')
    oif.write('fit RV         = %s\n' % fit_rv)
    oif.write('fit Transit    = %s\n' % fit_tr)
    oif.write('------------------------------\n')
    oif.write('Parameters and priors   \n')
    oif.write('------------------------------\n')
    max_length = max(len(name) for name in param_names)
    for i in range(len(prior_flags)):
        oif.write(f"Parameter {param_names[i]:<{max_length}}: {prior_flags[i]} [{prior_vals[2*i]:.4f},{prior_vals[2*i+1]:.4f}]\n")
    oif.write('------------------------------\n')
    oif.close()

    dummy_file = open(out_init_file)
    for line in dummy_file:
        print(line,end='')
    dummy_file.close()



