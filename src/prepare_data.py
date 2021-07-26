# -----------------------------------------------------------
#                       prepare_data.py
#  This file contains all the variable initializations,
#  both for RV and Transit fittings.
#                   O. Barragan, March 2016
# -----------------------------------------------------------

nconv = niter
nwalkers = nchains

# -----------------------------------------------------------
#                         RV DATA
# -----------------------------------------------------------

# Let us check the kind of variable
nplanets_rv = 0
for o in range(0, len(fit_rv)):
    nplanets_rv = nplanets_rv + int(fit_rv[o])

# Let us ensure that we want to fit rv data
if (nplanets_rv > 0):

    # Read the data file
    # time, RV, errors, and Telescope label
    try:
        time, rv, err, tspe = np.loadtxt('inpy/'+star+'/'+fname_rv[0], usecols=(0, 1, 2, 3),
                                     dtype={'names': ('time', 'rv', 'err', 'telescope'),
                                            'formats': ('float', 'float', 'float', 'S10')},
                                     comments='#', unpack=True)
    except:
        time, rv, err = np.loadtxt('inpy/'+star+'/'+fname_rv[0], usecols=(0, 1, 2),
                                     dtype={'names': ('time', 'rv', 'err'),
                                            'formats': ('float', 'float', 'float')},
                                     comments='#', unpack=True)
        tspe = ['0']*len(time)
        telescopes = ['0']
        telescopes_labels = ['0']
    if (is_special_jitter):
        sjitter = np.loadtxt(
            'inpy/'+star+'/'+fname_rv[0], usecols=(4), dtype=str, unpack=True)

    # These lists have lists with data for
    # the different telescopes
    time_all = []
    rv_all = []
    errs_all = []

    # Number of telescopes
    nt = len(telescopes)

    if (nt < 1):
        print('Please, indicate the telescope labels!')
        sys.exit('')

    # Separate the data for each telescope and create labels
    for i in range(0, nt):
        time_dum = []
        rv_dum = []
        errs_dum = []
        for j in range(0, len(tspe)):
            if (b(tspe[j]) == telescopes[i]):
                time_dum.append(time[j])
                rv_dum.append(rv[j])
                errs_dum.append(err[j])
        # The *all variables are lists of lists, each list constains
        # a list with the data of each telescope
        time_all.append(time_dum)
        rv_all.append(rv_dum)
        errs_all.append(errs_dum)

    # The mega* variables contains all the data
    # All this is neccesary because you do not have
    # the same number of data for each telescope
    rv_vals = []
    rv_time = []
    rv_errs = []
    tlab = []
    jrvlab = []
    # create mega with data of all telescopes
    for i in range(0, nt):
        # fill the mega variable with all the data of the
        # telescope i
        for j in range(0, len(rv_all[i])):
            # tlab has the label of the telescope (an integer)
            # this is useful because matches with the index of
            # the mega variables
            tlab.append(i)
            rv_vals.append(rv_all[i][j])
            rv_time.append(time_all[i][j])
            rv_errs.append(errs_all[i][j])

    if (is_special_jitter):
        n_jrv = len(jrvvec)
        for o in range(0, len(tlab)):
            for ndum in range(0, len(jrvvec)):
                if ((sjitter[o]) == jrvvec[ndum]):
                    jrvlab.append(ndum)
                    break
    else:
        jrvlab = list(tlab)
        n_jrv = nt

    total_rv_fit = True

    if bin_rv > 0:
        rv_time, rv_vals, rv_errs = bin_data(rv_time, rv_vals, rv_errs, bin_rv)
        nbands = 1
        nradius = 1
        tlab = [0]*len(rv_time)
        jrvlab = [0]*len(rv_time)

else:
    n_jrv = 1
    nt = 1
    tlab = [0]
    jrvlab = list(tlab)
    rv_vals = [None]
    rv_time = [None]
    rv_errs = [None]
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
for o in range(0, len(fit_tr)):
    nplanets_tr = nplanets_tr + int(fit_tr[o])

if (nplanets_tr > 0):

    # Each transit planet hasa different file
    myn = len(fname_tr)
    xt = [None]*myn
    yt = [None]*myn
    et = [None]*myn

    for o in range(0, myn):

        filename = 'inpy/'+star+'/'+fname_tr[o]
        if (my_tr_err == 0):  # The error bars come from the input file
            dummyd, dummyf, dummye = np.loadtxt(filename, usecols=columns_tr,
                                                comments='#', unpack=True)
        else:  # The error bars are given in the input
            dummyd, dummyf = np.loadtxt(filename, usecols=[0, 1],
                                        comments='#', unpack=True)
            dummye = [my_tr_err]*len(dummyd)

        dummyd = dummyd + textra

        hdate = dummyd
        wflux = dummyf
        errs = dummye

        # Each element of these lists will have the information
        # of a given transit
        xt[o] = hdate
        yt[o] = wflux
        et[o] = errs

    # Let us put together the information of all the arrays
    # the mega* lists have the data of all the transits
    # in 1D array
    lc_time = np.concatenate(xt)
    lc_flux = np.concatenate(yt)
    lc_errs = np.concatenate(et)
    megap = [0]*len(lc_time)

    # Create the label vectors for each instrument and jitter
    if (len(bands) == 1):
        nbands = 1
        nradius = 1
        trlab = [0]*len(lc_time)
        jtrlab = [0]*len(lc_time)
    else:
        trlab = []
        jtrlab = []
        nbands = len(bands)
        nradius = 1
        if is_multi_radius:
            nradius = nbands
        instrument = np.loadtxt(filename, usecols=[3], dtype=str, unpack=True)
        for o in range(0, len(lc_time)):
            for m in range(0, nbands):
                if (instrument[o] == bands[m]):
                    trlab.append(m)
                    jtrlab.append(m)

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
    megap = [0]
    total_tr_fit = False
    is_jitter_tr = False
    fit_q1 = 'f'
    fit_q2 = 'f'

# TRANSIT DATA READY
# Take care with span_tr
if (len(span_tr) == 1 and nplanets > 1):  # The user did not change this option in input_file.py
    span_tr = [0.0]*nplanets

if (len(a_from_kepler) == 1 and nplanets > 1):
    a_from_kepler = [False]*nplanets

for o in range(0, nplanets):
    if a_from_kepler[o]:
        fit_a[o] = 'g'
        min_a[o], max_a[o] = pti.get_a_err(
            mstar_mean, mstar_sigma, rstar_mean, mstar_sigma, (max_P[o]+min_P[o])/2.)

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


# What transit data are we fitting
if n_cad.__class__ == int or t_cad.__class__ == float:
    if (lc_data == 'kepler_lc'):
        n_cad = [10]*nbands
        t_cad = [29.425 / 60. / 24.0]*nbands  # days
    elif (lc_data == 'kepler_sc'):
        n_cad = [1]*nbands
        t_cad = [1.5 / 60. / 24.0]*nbands  # days
    elif (lc_data == 'tess_sc'):
        n_cad = [10]*nbands
        t_cad = [2.0 / 60. / 24.0]*nbands  # days
    elif (lc_data == 'free'):
        # values given by the user
        n_cad = [n_cad]*nbands
        t_cad = [t_cad]*nbands
    else:  # The n_cad and t_cad vectors are defined in the input file
        n_cad = n_cad
        t_cad = t_cad

#Let us make the code complatible with the old pyaneti files that include is_den_a
if is_den_a:
    sample_stellar_density = is_den_a

if sample_stellar_density:  # For a multiplanet system the density has to be the same
        for o in range(1, nplanets):
            fit_a[o] = 'f'
