# -----------------------------------------------------------
#                    todo-py.py
# This file contains a lot of useful of python functions.
#	    Oscar Barragan, March, 2016
# -----------------------------------------------------------

# Useful libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Solve compatibility issues between bytes and str for python > 3
if sys.version_info < (3,):
    def b(x):
        return x
else:
    def b(x):
        return x.decode()

# -----------------------------------------------------------
# This function returns the phase of a temporal array given
# the period
# input: jd -> time vector in julian date to be escaled (days)
# T0 -> time zero of the transit (days)
# P  -> planet orbital period (days)
# output: x -> vector with the scaled values (days)
# -----------------------------------------------------------


def scale_period(jd, Tp, P):
    x = [None]*len(jd)
    for i in range(len(jd)):
        x[i] = ((jd[i] - Tp) % P) / P
    return x

# -----------------------------------------------------------
# planet_mass -> gives the planet mass from RV parameters
# input: mstar -> mass of the orbited star (solar masses)
# k    -> semi-amplitude of the RV (m/s)
# P    -> planet orbital period (days)
# ecc  -> orbit eccentricity (no unit)
# i    -> orbit inclination (radians)
# output: mp -> planet mass (solar masses)
# WARNING: the defalut value for i is pi/2 (90 degrees),
# if you do not know the orbit inclination, this function
# computes the planet mass times sin(i)
# -----------------------------------------------------------


def planet_mass(mstar, k, P, ecc, i=np.pi/2.):

    P = P * 24. * 3600.  # s

    # Let us make a guess by assuming mstar >> mp
    unoe = np.sqrt(1.-ecc*ecc)
    mpsin = k * (2. * np.pi * S_GM_SI / P)**(-1./3.) * \
        mstar**(2./3.) * unoe
    mp = mpsin / np.sin(i)

    # find the mass by solving the mass function,
    # this is useful for stars orbited by other stars

    f = [1.0]*len(P)
    cte = - unoe**3 * P * k**3 / 2. / np.pi / S_GM_SI
    sini = np.sin(i)
    flag = True
    # Start Newton-Raphson algorithm
    while (flag):
        f = cte + (mp * sini)**3 / (mstar + mp)**2
        df = mp**2 * sini**3 / (mstar + mp)**2 * (3. - 2. * mp / (mstar + mp))
        mp = mp - f/df
        for j in range(0, len(P)):
            # check that all the array elemets have converged
            if (f[j] > 1.e-8):
                flag = True
                break
            else:
                flag = False

    return mp

# -----------------------------------------------------------
# This routine calculates the stellar density
# Based on eq. 30 from Winn., 2014
# Assuming the companion is too small
# Input:
# P -> Period
# a -> semi-major axis
# Output:
# rho -> stellar density
# -----------------------------------------------------------


def get_rhostar(P, a):
    P = P * 24. * 3600.  # s
    rho = 3. * np.pi * a**3 / (G_cgs * P * P)
    return rho


# -----------------------------------------------------------
# Return the equilibrium temeprature given the stellar temperature
# albedo, stellar radius and distance to the star
def get_teq(Tstar, albedo, rstar, a):
    Tp = Tstar*(1.0 - albedo)**(0.25)
    Tp = (rstar/2.0/a)**(0.5) * Tp
    return Tp

# -----------------------------------------------------------
#  Smart priors, get the best values of the physical and
#  priors limits
# -----------------------------------------------------------


def smart_priors():
    # We are using global variables
    global fit_tr, fit_rv, fit_v0
    global tota_rv_fit, total_tr_fit
    global min_rv0, max_rv0, v0, min_k, max_k
    global min_P, max_P, min_t0, max_t0, \
        min_rp, max_rp, \
        min_i, max_i

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
            for o in range(0, nt):
                if (fit_v0_val == 'u'):
                    fit_v0[o] = 'u'
                    min_rv0[o] = min(rv_all[o]) - 0.5
                    max_rv0[o] = max(rv_all[o]) + 0.5
                    if is_linear_trend == 'u':
                        min_rv0[o] = min(rv_all[o]) - 1.0
                        max_rv0[o] = max(rv_all[o]) + 1.0
                else:
                    fit_v0[o] = 'f'
                    min_rv0[o] = 0.0
                    max_rv0[o] = 0.0
        else:
            if (len(fit_v0) != nt):
                print('fit_v0 does not match the number of telescope labels')
                sys.exit()

# -----------------------------------------------------------
# find_vals_perc -> find the median and the errors within
#  a 68% credible interval
# input:
#       x -> vector with a minimum size nconv
#   nconv -> the last nconv points to be taken account
#            in the gaussian fit
# output:
#        med -> median value
#	mine -> left error (50% - 16%)
#	maxe -> right error (84% - 50%)
# -----------------------------------------------------------


def find_vals_perc(x, sf=1.0, prob=68.3):
    # With a 68% confidence interval
    mnval = 50.0 - prob/2.0
    mxval = 50.0 + prob/2.0
    mine, med, maxe = np.percentile(x, [mnval, 50.0, mxval])
    maxe = (maxe - med) / sf
    mine = (med - mine) / sf

    return med, mine, maxe


# -----------------------------------------------------------
def best_value(vector, loglike, cual):
    if (cual == 'median'):
        result = np.median(vector)
    elif(cual == 'mode'):
        result = my_mode(vector)
    elif(cual == 'map'):
        maxindex = np.argmax(loglike)
        result = vector[maxindex]

    return result


# -----------------------------------------------------------
# This routine calculates the mode of a vector
# The vector in divided in bins and count the maximum value
def my_mode(vector, bins=50):
    dx = np.max(vector) - np.min(vector)
    dx = dx / bins
    b = np.sort(vector)
    i = 0
    j = 0
    o = 0
    limite = np.min(vector) + dx
    if (dx > 1e-10):  # the parameter is fixed
        while(o < len(b)):
            if (b[o] < limite):
                i = i + 1
                if (i > j):
                    j = i
                    maximo = limite - dx/2.0
                o = o + 1
            else:
                i = 0
                limite = limite + dx
    else:
        maximo = np.median(vector)

    return maximo

# -----------------------------------------------------------


def mode_and_99(vector):
    a = my_mode(vector)
    d, b, c = find_vals_perc(vector, sf=1.0, prob=98)
    b = d - b
    c = c + d

    return a, b, c

# -----------------------------------------------------------


def good_clustering(chi2, chain_lab, nconv, nwalkers):
    # Let us find the good indixes for the cluster
    # We have n walkers

    print('STARTING CHAIN CLUSTERING')
    print('Initial number of chains:', nwalkers)

    # Extract all the chains
    chi2_walkers = [None]*nwalkers
    chi2_mean = [None]*nwalkers
    walk_dummy = []
    for i in range(0, nwalkers):
        for j in range(0, len(chain_lab)):
            if (chain_lab[j] == i):
                walk_dummy.append(chi2[j])
        chi2_walkers[i] = walk_dummy
        walk_dummy = []

    # The mean of each walker
    for i in range(0, nwalkers):
        chi2_mean[i] = np.mean(chi2_walkers[i])

    # get the minimum chi2
    total_min = min(chi2_mean)

    good_chain = []
    # Let us kill all the walkers 5 times the minimum
    for i in range(0, nwalkers):
        if (chi2_mean[i]/total_min < 1.0 + clustering_delta):
            # We are saving the good chain labels
            good_chain.append(i)

    # Now we know how many good chains we have
    new_nwalkers = len(good_chain)

    print('Final number of chains:', new_nwalkers)

    # Let us save the good index
    good_index = []
    for i in range(0, len(chain_lab)):
        for j in range(0, len(good_chain)):
            if (chain_lab[i] == good_chain[j]):
                good_index.append(i)

    return good_index, new_nwalkers


# -----------------------------------------------------------
def good_clustering_fast(chi2, nconv, nwalkers):
    # Let us find the good indixes for the cluster
    # We have n walkers

    print('STARTING CHAIN CLUSTERING')
    print('Initial number of chains:', nwalkers)

    # This variable will have each walker information
    chi2_walkers = [None]*nwalkers
    chi2_mean = [None]*nwalkers
    for i in range(0, nwalkers):
        chi2_walkers[i] = chi2[i::nconv]

      # The mean of each walker
    for i in range(0, nwalkers):
        chi2_mean[i] = np.mean(chi2_walkers[i])

    # get the minimum chi2
    total_min = min(chi2_mean)

    good_chain = []
    # Let us kill all the walkers 5 times the minimum
    for i in range(0, nwalkers):
        if (chi2_mean[i]/total_min < 1.0 + clustering_delta):
            # We are saving the good chain labels
            good_chain.append(i)

    new_nwalkers = len(good_chain)

    print('Final number of chains:', new_nwalkers)
    return good_chain, new_nwalkers


def good_clustering_likelihood(pos, nconv, nwalkers):
    # Let us find the good indixes for the cluster
    # We have n walkers

    print('STARTING CHAIN CLUSTERING')
    print('Initial number of chains:', nwalkers)

    # This variable will have each walker information
    pos_walkers = [None]*nwalkers
    pos_mean = np.zeros(nwalkers)
    for i in range(0, nwalkers):
        pos_walkers[i] = pos[i::nconv]

    # The mean of each walker
    for i in range(0, nwalkers):
        pos_mean[i] = np.mean(pos_walkers[i])


    sorted_indices = np.argsort(pos_mean)

    good_chain = []
    # Let us kill all the walkers 5 times the minimum
    for m,i in enumerate(sorted_indices):
        otros = i != sorted_indices
        otros[0:m] = False
        otros_walkers = np.mean(pos_mean[otros]) - 3*np.std(pos_mean[otros])
        if (pos_mean[i] > otros_walkers):
            # We are saving the good chain labels
            good_chain.append(i)
        else:
            print("removing chain {}, iteration {}".format(i,m))

    new_nwalkers = len(good_chain)

    print('Final number of chains:', new_nwalkers)

    return good_chain, new_nwalkers

# -----------------------------------------------------------


def clustering(par, good_index):

    cluster_par = np.zeros(len(good_index))
    for i in range(0, len(good_index)):
        n = good_index[i]
        cluster_par[i] = par[n]

    return cluster_par

# -----------------------------------------------------------


def clustering_fast(par, good_index, nconv):

    dummy_par = []
    for i in good_index:
        dummy_par.append(par[i::nwalkers])

    cluster_par = np.ndarray(len(good_index)*nconv)

    n = 0
    for i in range(0, len(dummy_par)):
        for j in range(0, len(dummy_par[i])):
            cluster_par[n] = dummy_par[i][j]
            n = n + 1

    return cluster_par

# -----------------------------------------------------------


def bin_data(tvector, fvector, rvector, tbin):
    leftt = min(tvector)
    right = leftt + tbin
    xbined = []
    fbined = []
    rbined = []
    while (leftt < max(tvector)):
        fdummy = []
        rdummy = []
        for i in range(0, len(tvector)):
            if (tvector[i] > leftt and tvector[i] < right):
                fdummy.append(fvector[i])
                rdummy.append(rvector[i])
        fbined.append(np.mean(fdummy))
        rbined.append(np.mean(rdummy))
        xbined.append((leftt + tbin/2.))
        leftt = leftt + tbin
        right = right + tbin
    xbined = np.asarray(xbined)
    fbined = np.asarray(fbined)
    rbined = np.asarray(rbined)
    xbined = xbined[~np.isnan(rbined)]
    fbined = fbined[~np.isnan(rbined)]
    rbined = rbined[~np.isnan(rbined)]
    return xbined, fbined, rbined



# -----------------------------------------------------------


def print_values(vector, var, vartex, unit, unittex):
    # fname is the variable where we are writting the numbers
    # vector is the posterior vectors
    # var is the label of the variable to save
    # unit is the label of the variable to save
    # vals do we want to print median or mode
    medv, minv, maxv = find_vals_perc(vector, s_factor)
    nd = 1
    if (abs(minv) > 1e-20 and abs(maxv) > 1e-20):
        try:
            nd = int(np.log10(max(1./minv, 1./maxv))) + 2
        except:
            pass
    opars.write('{:10s} = {:4.7f} - {:4.7f} + {:4.7f} {:8s} \n'.format(var, medv, minv, maxv, unit))
    otex.write('\\newcommand{\\'+vartex+'}[1]['+unittex+'] \
  {$'+str(round(medv, nd))+' _{ - '+str(round(minv, nd))+' } ^ { + '+str(round(maxv, nd))+' }$~#1} \n')
    if (is_print_mode):
        medv, minv, maxv = mode_and_99(vector)
        opars.write('%10s  %4.7f , %4.7f , %4.7f %8s (mode, 1 percent, 99 percent) \n' % (
            '', medv, minv, maxv, unit))


def joint_fit():
    global fit_all, fit_ldc, fit_rvs, nt
    global mstar_mean, rstar_mean, mstar_sigma_rstar_sigma
    global is_log_P, is_ew, is_b_factor, is_log_k, is_log_rv0
    global fit_t0, fit_P, fit_e, fit_w, fit_i, fit_a, fit_q1, fit_q2, fit_rp, fit_k, fit_v0
    global T0, P, e, w, ii, a, q1, q2, rp, k0, alpha, beta, v0
    global min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w, min_i, max_i, min_a,\
        max_a, min_q1, max_q1, min_q1, max_q1, min_rp, max_rp, min_k, max_k, min_alpha, max_alpha, \
        min_beta, max_beta, min_rv0, max_rv0
    global vari, chi2, chi2red, t0o, Po, eo, wo, io, ao, q1o, q2o, rpo, ko, alphao, betao, vo, what_fit
    global new_nwalkers, good_index, nwalkers
    global jrvo, jtro, total_fit_flag, flags
    global limits, priorf, priorl, limits_ldc, limits_rvs
    global prior_flags, prior_vals, model_int, model_double
    global n_jtr, trlab, jtrlab, jtr_prior_flag, jtr_prior_vals, jrv_prior_flag, jrv_prior_vals
    global kernels, krv_labels, krv_prior_flag, krv_prior_vals
    global ktr_labels, ktr_prior_flag, ktr_prior_vals, np_rv, np_tr

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

    flags = [is_log_P, is_ew, is_b_factor, sample_stellar_density, is_log_k, is_log_rv0]

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
    for o in range(0, nbands):
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
        jrv_prior_flag = ['m']*n_jrv
        jrv_prior_vals = [1e-3,1]*n_jrv
        #jrv_prior_vals = [0.0, 0.05]*n_jrv
    else:
        jrv_prior_flag = ['f']*n_jrv
        jrv_prior_vals = [0., 0.5]*n_jrv

    # Transit jitter priors
    if is_jitter_tr:
        if len(jtr_prior_flag) < 1:
            jtr_prior_flag = ['m']*n_jtr
            jtr_prior_vals = [1e-3, 1.e-2]*n_jtr
    else:
        n_jtr = 0
        trlab = [0]*len(lc_time)
        jtrlab = [0]*len(lc_time)
        jtr_prior_flag = ['u']*n_jtr
        jtr_prior_vals = [0., 1.e-3]*n_jtr

    # Trends priors
    trends_prior_flag = [is_linear_trend, is_quadratic_trend]
    val_linear = 0.1
    val_quadratic = 0.1
    if (is_linear_trend == 'f'):
        val_linear = 0.0
    if (is_quadratic_trend == 'f'):
        val_quadratic = 0.0
    trends_prior_vals = [-val_linear,
                         val_linear, -val_quadratic, val_quadratic]

    # RV Kernels
    if kernel_rv == 'None':
        krv_prior_flag = []
        krv_prior_vals = []
        krv_labels = []
    elif kernel_rv == 'QP2':
        # fit_krv has to be a four-dimensional vector (A,Gamma1,Gamma2,P)
        krv_prior_flag = fit_krv
        # This has to be a 4-dimensional vector with the prior limits
        krv_prior_vals = krv_priors
        krv_labels = ['A', 'Gamma', 'metric', 'P_GP']
    elif kernel_rv == 'M32' or kernel_rv == 'M52':
        krv_prior_flag = fit_krv
        krv_prior_vals = krv_priors
        krv_labels = ['A','metric']
    elif kernel_rv == 'QPK':
        # fit_krv has to be a four-dimensional vector (A,Gamma1,Gamma2,P)
        krv_prior_flag = fit_krv
        # This has to be a 4-dimensional vector with the prior limits
        krv_prior_vals = krv_priors
        krv_labels = ['A', 'le', 'lp', 'P_GP']
    # Add remaining kernels
    elif kernel_rv[0:2] == 'MQ':
        krv_prior_flag = fit_krv
        krv_prior_vals = krv_priors
        krv_labels = [None]*len(fit_krv)
        for m in range(len(fit_krv) - 3):
            krv_labels[m] = 'A'+str(m)
        krv_labels[len(fit_krv)-3] = 'lambdae'
        krv_labels[len(fit_krv)-2] = 'lambdap'
        krv_labels[len(fit_krv)-1] = 'PGP'
    elif kernel_rv[0:2] == 'ME' or kernel_rv[0:2] == 'MM':
        krv_prior_flag = fit_krv
        krv_prior_vals = krv_priors
        krv_labels = [None]*len(fit_krv)
        for m in range(len(fit_krv)- 1):
            krv_labels[m] = 'A'+str(m)
        krv_labels[len(fit_krv)-1] = 'lambdae'

    np_rv = len(krv_prior_flag)

    # TR Kernels
    if kernel_tr == 'None':
        ktr_prior_flag = []
        ktr_prior_vals = []
        ktr_labels = []
    elif kernel_tr == 'QPK':
        # fit_krv has to be a four-dimensional vector (A,Gamma1,Gamma2,P)
        ktr_prior_flag = fit_ktr
        # This has to be a 4-dimensional vector with the prior limits
        ktr_prior_vals = ktr_priors
        ktr_labels = ['A', '$\lambda_e$', '$\lambda_p$', 'P']
    elif kernel_tr == 'QP2':
        # fit_krv has to be a four-dimensional vector (A,Gamma1,Gamma2,P)
        ktr_prior_flag = fit_ktr
        # This has to be a 4-dimensional vector with the prior limits
        ktr_prior_vals = ktr_priors
        ktr_labels = ['A','metric', '$\Gamma$', 'P']
    # Add remaining kernels

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

    if (method == 'mcmc'):

        # Ensure nwalkers is divisible by 2
        if (nwalkers % 2 != 0):
            nwalkers = nwalkers + 1

        pti.mcmc_stretch_move(
            rv_time, rv_vals, lc_time, lc_flux, rv_errs, lc_errs,
            tlab, jrvlab, trlab, jtrlab,
            flags, total_fit_flag, prior_flags, prior_vals,
            kernels, model_int,
            model_double,
            nwalkers, maxi, thin_factor, nconv)

        #Change fortran file to a lighter csv file
        data = np.loadtxt('all_data.dat',unpack=True)
        np.savetxt('all_data.dat',data.T,delimiter=',',fmt='%1.12e')

    elif (method == 'plot'):
        print('I will only print the values and generate the plot')

    else:
        print('You did not choose a method!')
        print('method = mcmc   -> Run the MCMC code')
        print('method = plot   -> Plot of a previous run')
        sys.exit('choose your favorite.')



#
    newfile = outdir+'/'+star+'_all_data.dat'
    if (os.path.isfile('all_data.dat')):
        os.rename('all_data.dat', newfile)

        # Define the labels to be used in the plots
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
        for o in range(0, nplanets):
            etiquetas = ['t0'+plabels[o], 'p'+plabels[o], elab+plabels[o],
                         wlab+plabels[o], ilab+plabels[o], alab+plabels[o],
                         klab+plabels[o]]
            labs.append(etiquetas)
        # planet radius labels
        for o in range(0, nplanets):
            for m in range(0, nradius):
                labs.append(['rprstar'+plabels[o]+bands[m]])
        # LDC labels
        for m in range(0, nbands):
            labs.append(['q1'+bands[m], 'q2'+bands[m]])
        # RV instrument labels
        labs.append(telescopes_labels)
        # jitter labels
        for o in range(0, n_jrv):
            labs.append(['rv_jitter'+str(telescopes_labels[o])])
        for o in range(0, n_jtr):
            labs.append(['vv_jitter'+str(bands[o])])
        # trends labels
        labs.append(['linear_trend'])
        labs.append(['quadratic_trend'])
        labs.append(krv_labels)
        labs.append(ktr_labels)
        # Total labels vector

        def line_prepender(filename, text):
            with open(filename, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                for o in range(len(text)):
                    f.write(text[o]+'   ')
                f.write('\n' + content)

        labels = np.concatenate([['#i'], ['log_posterior'], ['chi2_rv'], [
                                'chi2_tr'], np.concatenate(labs)])
        line_prepender(newfile, labels)
# -----------------------------------------------------------
#          PRINT INITIAL CONFIGURATION
# -----------------------------------------------------------


def print_init():
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
#  if ( total_tr_fit ):
#    oif.write ('LC data        = %s\n' %lc_data)
#    oif.write ('cadence time   =  %2.3f min\n'%(t_cad*60.*24))
#    oif.write ('n rebinning    = %d\n' %n_cad)
    for j in range(0, nplanets):
        oif.write('------------------------------\n')
        oif.write('  PLANET %s \n' % (star + plabels[j]))
        oif.write('------------------------------\n')
        oif.write('        PRIOR RANGES          \n')
        oif.write('------------------------------\n')
        oif.write('T0 = %s[ %4.4f , %4.4f ]\n' %
                  (fit_t0[j], min_t0[j], max_t0[j]))
        oif.write('P  = %s[ %4.4f , %4.4f ]\n' %
                  (fit_P[j], min_P[j], max_P[j]))
        if (is_ew):
            oif.write('ew1= %s[ %4.4f , %4.4f ]\n' %
                      (fit_ew1[j], min_ew1[j], max_ew1[j]))
            oif.write('ew2= %s[ %4.4f , %4.4f ]\n' %
                      (fit_ew2[j], min_ew2[j], max_ew2[j]))
        else:
            oif.write('e  = %s[ %4.4f , %4.4f ]\n' %
                      (fit_e[j], min_e[j], max_e[j]))
            oif.write('w  = %s[ %4.4f , %4.4f ]\n' %
                      (fit_w[j], min_w[j], max_w[j]))
        if (is_b_factor):
            oif.write('b  = %s[ %4.4f , %4.4f ]\n' %
                      (fit_b[j], min_b[j], max_b[j]))
        else:
            oif.write('i  = %s[ %4.4f , %4.4f ]\n' %
                      (fit_i[j], min_i[j], max_i[j]))
        oif.write('a  = %s[ %4.4f , %4.4f ]\n' %
                  (fit_a[j], min_a[j], max_a[j]))
        oif.write('rp = %s[ %4.4f , %4.4f ]\n' %
                  (fit_rp[j], min_rp[j], max_rp[j]))
        oif.write('K  = %s[ %4.4f , %4.4f ]\n' %
                  (fit_k[j], min_k[j], max_k[j]))
        oif.write('------------------------------\n')
    oif.write(' Other parameter priors \n')
    oif.write('------------------------------\n')
    for m in range(0, nbands):
        oif.write('q1 %s = %s[ %4.4f , %4.4f ]\n' %
                  (bands[m], fit_q1[m], min_q1[m], max_q1[m]))
        oif.write('q2 %s = %s[ %4.4f , %4.4f ]\n' %
                  (bands[m], fit_q2[m], min_q2[m], max_q2[m]))
    for m in range(0, nt):
        oif.write('%s = %s[ %4.4f , %4.4f ]\n' % (
            telescopes_labels[m], fit_v0[m], min_rv0[m], max_rv0[m]))
    oif.write('==============================\n')

    oif.close()
    dummy_file = open(out_init_file)
    for line in dummy_file:
        print(line)
    dummy_file.close()

#------------------------------------------------------------------------#
#            Automatic creation of input for TANGO                       #
#------------------------------------------------------------------------#


def tango_params(param, vec, parhs=True):
    vlen = len(vec)
    letra = param + ' = '
    if (parhs):
        letra = letra + ' [ '
    for o in range(0, vlen):
        letra = letra + str(np.median(vec[o]))
        if (o < vlen - 1):
            letra = letra + ','

    if (parhs):
        letra = letra + ' ]'
    letra = letra + '\n'

    return letra

# This routine create an input file to create animations using tango
# github.com/oscaribv/tango


def create_tango_input():
    tangof = outdir+'/'+star+'_tango_input.py'
    tf = open(tangof, 'w')

    tf.write('#Input file for tango\n')
    tf.write('#system:' + star+'\n')
    tf.write('#Created automatically with pyaneti\n')

    tf.write('\n')

    tf.write('#Data file with the flattened light curve\n')
    tf.write('lcname = \''+star+'_new_lc.dat\'\n')

    tf.write('#--------------------------------------------------------------------\n')
    tf.write('#                 Planet and orbit parameters\n')
    tf.write('# Each parameter is a list in which each element\n')
    if (nplanets == 1):
        tf.write('# correspond to a planet. For this case, there is ' +
                 str(nplanets)+' planet\n')
    else:
        tf.write('# correspond to a planet. For this case, there are ' +
                 str(nplanets)+' planets\n')
    tf.write('#--------------------------------------------------------------------\n')

    tf.write('\n')

    # Orbital period (days)
    tf.write(tango_params('P', P_vec))
    tf.write(tango_params('T0', T0_vec))
    tf.write(tango_params('e', e_vec))
    tf.write(tango_params('w', w_vec))
    tf.write(tango_params('a', ar_vec))
    tf.write(tango_params('inclination', i_vec))
    tf.write(tango_params('rp', rr_vec))
    tf.write(tango_params('u1', [u1_vec], False))
    tf.write(tango_params('u2', [u2_vec], False))

# Integration time of the data
    tf.write('t_cad = ' + str(t_cad) + ' \n')
    tf.write('n_cad = ' + str(n_cad) + ' \n')

    tf.write('\n')

    # Calculate the transit duration for planet b to create the time ranges
    tfull = np.median(trt_vec[0])
    tfull = tfull/24.
    mit0 = np.median(T0_vec[0])

    tf.write('#--------------------------------------------------------------------\n')
    tf.write('#              Animation controls \n')
    tf.write('#--------------------------------------------------------------------\n')
    tf.write('#Window size to show the data (days)\n')
    tf.write('size_time = 0.5\n')
    tf.write(
        '#1./(photograms per day) in this case the code will create a photogram each 7.2 min\n')
    tf.write('vel_time  = 1./200.\n')
    tf.write('#Animation minimum time (Be sure that you are using the same units as in your data file)\n')
    tf.write('tmin = '+str(mit0 - 2*tfull)+'\n')
    tf.write('#Animation maximum time (Be sure that you are using the same units as in your data file)\n')
    tf.write('tmax = '+str(mit0 + 2*tfull)+'\n')

    tf.write('\n')

    tf.write('#--------------------------------------------------------------------\n')
    tf.write('#                     Plot controls\n')
    tf.write('#--------------------------------------------------------------------\n')

    tf.write('\n')

    tf.write('#Control if we overplot the light curve model\n')
    tf.write('#You need to have installed pyaneti in your computer to use it\n')
    tf.write('is_plot_model = False\n')

    tf.write('\n')

    tf.write('#-----------------------------------------------------------------\n')
    tf.write('#                         END\n')
    tf.write('#-----------------------------------------------------------------\n')
