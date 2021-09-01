# Read the data
# Dummy params vector contains
# [0] -> i
# [1] -> chain label
# [2] -> chi2_rv
# [3] -> chi2_tr
# [4 - 7*nplanets] -> Planet parameters(T0,P,e,w,b,a/R*,K)
# [4+7*nplanets - nplanets*nbands ] -> rps_ij (i to planet, j to band)
# [] LDC bands
# [] Offsets
# [] jitter
# [] trends
srp = 7*nplanets  # start of planet radius
sldc = srp + nplanets*nradius  # start of LDC
srv = sldc + nldc*nbands  # start of RV offsets
sjitrv = srv + nt  # start of RV jitter
sjittr = sjitrv + n_jrv  # start of TR jitter
strends = sjittr + n_jtr  # start of RV trends
skrv = strends + 2  # start of RV kernel parameters
sktr = skrv + np_rv  # start of TR kernel parameters

newfile = outdir+'/'+star+'_all_data.dat'

#If we are using the new pyaneti, then the file is in csv format
try:
    dparams = np.loadtxt(newfile, comments='#', unpack=True,delimiter=',')
except: #use the old pyaneti files
    dparams = np.loadtxt(newfile, comments='#', unpack=True)

# Let us do the clustering
params = list(dparams)
#par_likelihood = list(ldparams[2])
new_nwalkers = nwalkers
# The maximum log(likelihood)
dmaxloglike = dparams[1]
maxloglike = list(dmaxloglike)
if (is_clustering):
    good_index, new_nwalkers = good_clustering_likelihood(
        dparams[1], nconv, nwalkers)
    for o in range(0, len(dparams)):
        params[o] = clustering_fast(dparams[o], good_index, nconv)
    maxloglike = params[1]


# Create the stellar data
mstar = np.random.normal(
    loc=mstar_mean, scale=mstar_sigma, size=new_nwalkers*nconv)
rstar = np.random.normal(
    loc=rstar_mean, scale=rstar_sigma, size=new_nwalkers*nconv)
tstar = np.random.normal(
    loc=tstar_mean, scale=tstar_sigma, size=new_nwalkers*nconv)
vsini = np.random.normal(
    loc=vsini_mean, scale=vsini_sigma, size=new_nwalkers*nconv)

# Calculate the BIC
ndata = 0
if (total_rv_fit):
    ndata = ndata + len(rv_vals)
if (total_tr_fit):
    ndata = ndata + len(lc_time)

npars = 0

if (len(plot_parameters) < 2):

    # plot_parameters stores the flag of each fitted parameter
    plot_parameters = []

    for o in range(0, len(prior_flags)):
        if (prior_flags[o] != 'f'):
            npars = npars + 1
            plot_parameters.append(o)

dummy_pars = [0.0]*len(params)
for o in range(0, len(params)):
    dummy_pars[o] = best_value(params[o], maxloglike, get_value)

log_like_total, chi2tot_val_rv, chi2tot_val_tr = \
    pti.get_loglike(rv_time, rv_vals, lc_time, lc_flux, rv_errs, lc_errs,
                    tlab, jrvlab, trlab, jtrlab, total_fit_flag, flags, kernels,
                    dummy_pars[4:], model_int, model_double)

bic_from_loglikelihood = np.log(ndata)*npars - 2.0*log_like_total
aic_from_loglikelihood = 2.0*npars - 2.0*log_like_total

chi2tot_val = chi2tot_val_rv + chi2tot_val_tr
chi2_val = chi2tot_val / (ndata - npars)

if (scale_error_bars):
    s_factor = np.sqrt(chi2_val)
    if (chi2_val > 1.0):
        s_factor = 1.0 / s_factor
else:
    s_factor = 1.0


if (method == 'mcmc' or method == 'plot'):

    base = 4  # Where do the parameters start?
# Fitted parameters
    T0_vec = [None]*nplanets
    P_vec = [None]*nplanets
    e_vec = [None]*nplanets
    w_vec = [None]*nplanets
    b_vec = [None]*nplanets
    ar_vec = [None]*nplanets
    k_vec = [None]*nplanets
    rr_vec = [None]*nplanets*nradius
# Derived parameters
    Teq_vec = [None]*nplanets  # Planet temperature
    r_vec = [None]*nplanets*nradius
    a_vec = [None]*nplanets  # semi-major axis
    m_vec = [None]*nplanets  # planet mass
    i_vec = [None]*nplanets  # orbit inclination
    Tpe_vec = [None]*nplanets  # Periastron passage time
    ds_vec = [None]*nplanets  # stellar density
    dp_vec = [None]*nplanets  # planet density
    gp_vec = [None]*nplanets  # planet surface gravity
    trt_vec = [None]*nplanets  # Total transit duration
    tri_vec = [None]*nplanets  # Ingress/egress duration


# Print the summary
    out_params_file = outdir+'/'+star+'_params.dat'
    out_tex_file = outdir+'/'+star+'_params.tex'
    opars = open(out_params_file, 'w')
    otex = open(out_tex_file, 'w')
    opars.write('\n')
    opars.write(
        '--------------------------------------------------------------\n')
    opars.write('Summary:\n')
    opars.write('N_chains         = %8i \n' % nwalkers)
    opars.write('N_iter           = %8i \n' % nconv)
    opars.write('thin_factor      = %8i \n' % thin_factor)
    opars.write('N_rv_data        = %8i \n' % len(rv_time))
    opars.write('N_tr_data        = %8i \n' % len(lc_time))
    opars.write('N_data           = %8i \n' % ndata)
    opars.write('N_pars           = %8i \n' % npars)
    opars.write('chi2_rv          = %4.4f\n' % (chi2tot_val_rv))
    opars.write('chi2_tr          = %4.4f\n' % (chi2tot_val_tr))
    opars.write('chi2             = %4.4f\n' % (chi2tot_val))
    opars.write('dof              = %8i \n' % (ndata - npars))
    opars.write('chi2/dof         = %4.4f \n' % chi2_val)
    opars.write('ln likelihood    = %4.4f\n' % (log_like_total))
    opars.write('BIC              = %4.4f\n' % (bic_from_loglikelihood))
    opars.write('AIC              = %4.4f\n' % (aic_from_loglikelihood))
    opars.write(
        '--------------------------------------------------------------\n')
    opars.write('             INPUT STELLAR PARAMETERS\n')
    opars.write(
        '--------------------------------------------------------------\n')
    opars.write('M_*     = %4.4f - %4.4f + %4.4f solar masses\n' %
                (mstar_mean, mstar_sigma, mstar_sigma))
    opars.write('R_*     = %4.4f - %4.4f + %4.4f solar radii\n' %
                (rstar_mean, rstar_sigma, rstar_sigma))
    opars.write('T_*     = %4.4f - %4.4f + %4.4f K\n' %
                (tstar_mean, tstar_sigma, tstar_sigma))
    opars.write('vsini     = %4.4f - %4.4f + %4.4f km/s (for RM effect estimation)\n' %
                (vsini_mean, vsini_sigma, vsini_sigma))
    opars.write('J mag     = %4.4f  (for TSM estimation)\n' %
                (mag_j))
    # tex
    #otex.write ('--------------------------------------------------------------\n')
    #otex.write ('             INPUT STELLAR PARAMETERS\n')
    #otex.write ('--------------------------------------------------------------\n')
    otex.write('\\newcommand{\smass}[1][$M_{\odot}$]{ $ %4.7f _{- %4.7f}^{ + %4.7f} $ #1} \n' % (
        mstar_mean, mstar_sigma, mstar_sigma))
    otex.write('\\newcommand{\sradius}[1][$R_{\odot}$]{ $%4.7f _{ - %4.7f}^{ + %4.7f} $ #1}\n' % (
        rstar_mean, rstar_sigma, rstar_sigma))
    otex.write('\\newcommand{\stemp}[1][$\mathrm{K}$]{ $ %4.7f _{- %4.7f}^{ + %4.7f} $ #1 }\n' % (
        tstar_mean, tstar_sigma, tstar_sigma))

    # Print the data for all the planets
    for o in range(0, nplanets):
        T0_vec[o] = np.asarray(list(params[base + 0]))
        P_vec[o] = np.asarray(list(params[base + 1]))
        e_vec[o] = np.asarray(list(params[base + 2]))
        w_vec[o] = np.asarray(list(params[base + 3]))
        b_vec[o] = np.asarray(list(params[base + 4]))
        ar_vec[o] = np.asarray(list(params[base + 5]))
        k_vec[o] = np.asarray(list(params[base + 6]))
        for m in range(0, nradius):
            rr_vec[o*nradius +
                   m] = list(np.asarray(list(params[4+srp+nradius*o+m])))
            rr_vec[o*nradius+m] = np.asarray(rr_vec[o*nradius+m])

# STARTING CALCULATIONS

        if (is_log_P):
            P_vec[o] = 10.0**(P_vec[o])
        if (sample_stellar_density):
            if (o == 0):
                miden = list(params[base+5])
            for m in range(0, len(miden)):
                ar_vec[o][m] = pti.rhotoa(miden[m], P_vec[o][m], 1)
        if (is_log_k):
            k_vec[o] = 10.0**(k_vec[o])

        if (is_ew):
            e_dum = list(e_vec[o])
            w_dum = list(w_vec[o])
            e_vec[o], w_vec[o] = pti.ewto(e_dum, w_dum)
            w_vec[o] = w_vec[o]  # % (2*np.pi)

    # Change between b and i
        if (is_b_factor):
            i_vec[o] = pti.btoi(b_vec[o], ar_vec[o], e_vec[o], w_vec[o])
        else:
            # calculate the impact parameter (eq. 7 Winn 2014)
            i_vec[o] = list(b_vec[o])
            b_vec[o] = ar_vec[o] * np.cos(b_vec[o]) * ((1. - e_vec[o]**2)
                                                       / (1.0 + e_vec[o]*np.sin(w_vec[o])))
            i_vec[o] = np.array(i_vec[o])

        # Calculate equilibrium temperature
        # assuming albedo=0
        Teq_vec[o] = get_teq(tstar, 0.0, 1.0, ar_vec[o])

        # Get the star periastron pasage
        w_s_deg, w_s_deg_l, w_s_deg_r = find_vals_perc(
            w_vec[o]*180./np.pi, s_factor)
        # planet periastron passage
        w_p_deg = (w_s_deg + 180.) % 360

    # Transit durations aproximations (eq. 14, 15, 16 from Winn 2014)
        ec_factor = np.sqrt(
            (1. - e_vec[o]*e_vec[o])) / (1.0 + e_vec[o]*np.sin(w_vec[o]))
        #total transit duration (from contact 1 to 4)
        trt_vec[o] = np.sqrt((1. + rr_vec[o])**2 - b_vec[o]
                             ** 2) / (ar_vec[o] * np.sin(i_vec[o]))
        trt_vec[o] = P_vec[o] / np.pi * \
            np.arcsin(trt_vec[o]) * ec_factor * 24.0
        #full transit duration (from contract 2 to 3)
        tri_vec[o] = np.sqrt((1. - rr_vec[o])**2 - b_vec[o]
                             ** 2) / (ar_vec[o] * np.sin(i_vec[o]))
        tri_vec[o] = P_vec[o] / np.pi * \
            np.arcsin(tri_vec[o]) * ec_factor * 24.0
        # tri_vec[o] = ( trt_vec[o] - tri_vec[o] ) / 2.0 #ingress egress time

        # Calculate the star density from transit data
        # Eq. (30) Winn 2014
        ds_vec[o] = get_rhostar(P_vec[o], ar_vec[o])  # cgs

        # Time of periastron passage
        Tpe_vec[o] = list(T0_vec[o])
        for m in range(0, len(Tpe_vec[o])):
            Tpe_vec[o][m] = pti.find_tp(
                T0_vec[o][m], e_vec[o][m], w_vec[o][m], P_vec[o][m])

        # Density from the input stellar parameters
        irho_vec = mstar/rstar**3 * 1.411

        # Get planet mass, radius and orbit semi-major axis in real units
        for m in range(0, nradius):
            r_vec[o*nradius+m] = rr_vec[o*nradius+m] * rstar

        a_vec[o] = ar_vec[o] * rstar * S_radius_SI / AU_SI
        m_vec[o] = planet_mass(mstar, k_vec[o]*1.e3,
                               P_vec[o], e_vec[o], i_vec[o])

        if is_single_transit:
#            # Let us make the magic if we are dealing with single transits
#            # we have the scaled radius rr_vec and the impact parameter b_vec, we need the transit durations
            # Let us estimate the transit duration from all the points where the planet is crossing the star
            z_vec = [None]*len(Tpe_vec[o])
            tr_dur_vec = []
            for m in range(len(Tpe_vec[o])):
                #pars = [np.mean(Tpe_vec[o]),np.mean(P_vec[o]),np.mean(e_vec[o]),np.mean(w_vec[o]),np.mean(i_vec[o]),np.mean(ar_vec[o])]
                pars = [Tpe_vec[o][m], P_vec[o][m], e_vec[o][m], w_vec[o][m], i_vec[o][m], ar_vec[o][m]]
                # This works for a single planet and for a single band
                rplanet = np.mean(rr_vec[o][m])
                tiempo = np.arange(min(lc_time), max( lc_time), 1./60./24.)  # minute precision
                z_vec[m] = pti.find_z(tiempo, pars)
                # Z contains is the vector containing the distance between the star and planet centres in stellar radius units
                tr_index = z_vec[m] < 1 + rplanet
                tr_time = tiempo[tr_index]
                # Now we have the times where the transit happens
                tr_dur_vec.append(max(tr_time) - min(tr_time))
            tr_dur_vec = np.array(tr_dur_vec)
            # Now we have a vector with the estimated transit duration, time to estimate planet velocity
            # Using eq. (1) of Osborn et al., 2016, MNRAS, 457, 2273, https://academic.oup.com/mnras/article/457/3/2273/2588921
            vpl_vec = np.array(2.*np.sqrt((1.+np.array(rr_vec))**2-np.array(b_vec)**2)/(tr_dur_vec*24.*60.*60.))
            vpl_vec = vpl_vec[0]
            # Now we can estimate planetary orbit period assuming the orbit is circular
            # Using eq. (2) of Osborn et al., 2016, MNRAS, 457, 2273
            P_circ_vec = 8. * np.pi**2 * G_cgs * irho_vec / 3. / vpl_vec**3
            # Convert P_circ_vec to days
            P_circ_vec = P_circ_vec/24./3600.

        # Jeans scape parameter (Fossati et al., 2017, A&A, 598, A90)
        jsp = k_btz * np.asarray(Teq_vec[o]) * \
            S_radius_cgs * np.asarray(r_vec[o*nradius])
        jsp = S_GM_cgs * np.asarray(m_vec[o]) * m_hyd / jsp

        # Planet-star distance at the time of eclipse
        true_anomaly_vec = [None]*len(T0_vec[o])
        for l in range(0, len(true_anomaly_vec)):
            true_anomaly_vec[l] = pti.find_anomaly(
                T0_vec[o][l], T0_vec[o][l], e_vec[o][l], w_vec[o][l], P_vec[o][l])
        dummy_anomaly = np.concatenate(true_anomaly_vec)
        psd_vec = ar_vec[o] * (1. - e_vec[o]**2) / \
            (1. + e_vec[o]*np.cos(dummy_anomaly))
        psd_vec_units = psd_vec * rstar * S_radius_SI / AU_SI

        # Kepler cociente
        pa_vec = (P_vec[o]*3600.*24.0)**2 * S_GM_SI * (mstar + m_vec[o])
        pa_vec = pa_vec / (4.*np.pi**2 * (a_vec[o]*AU_SI)**3)

        # stimate planet gravity and density
        pden_vec = m_vec[o] / r_vec[o*nradius]**3  # solar units
        pden_vec = pden_vec * S_den_cgs  # g/cm^3

        # We can stimate planet surface gravity (eq. (31) Winn)
        pgra_vec = (P_vec[o]*24.*3600.) * \
            (rr_vec[o*nradius]/ar_vec[o])**2 * np.sin(i_vec[o])
        pgra_vec = 2. * np.pi * \
            np.sqrt(1. - e_vec[o]**2) * (k_vec[o]*1.e5) / pgra_vec  # cm/s^2

        # Estimate surface gravity from the derived parameters
        pgra_vec2 = m_vec[o] / r_vec[o*nradius]**2  # in solar units
        pgra_vec2 = pgra_vec2 * 28.02 * 981.  # cm/s^2

        # Stellar luminosity in solar units
        Ls = (rstar)**2*(tstar/S_Teff)**4
        # planet insolation in Flux received at Earth
        Fp = Ls/a_vec[o]**2

        #Estimate the maximum magnitude difference of a star that can generate a transit
        #Based on eq. 4 of Vanderburg et al., (2019), ApJL, 881, L19
        #First lest estimate t12 and t_13
        t12_vec = (trt_vec[o] - tri_vec[o])/2. #assuming ingress and egress are the same
        t13_vec = tri_vec[o] + t12_vec
        delta_mag_vec = t12_vec**2/(t13_vec**2*rr_vec[o*nradius]**2)
        delta_mag_vec = 2.5*np.log10(delta_mag_vec)

        # Convert units
        usymbol = '{\odot}'
        if (unit_mass == 'earth'):
            usymbol = '{\oplus}'
            if (fit_rv[o]):
                m_vec[o] = m_vec[o] * S_GM_SI / E_GM_SI
            if (fit_tr[o]):
                for m in range(0, nradius):
                    r_vec[o*nradius+m] = r_vec[o*nradius+m] * \
                        S_radius_SI / E_radius_e_SI
        elif (unit_mass == 'jupiter'):
            usymbol = '\mathrm{J}'
            if (fit_rv[o]):
                m_vec[o] = m_vec[o] * S_GM_SI / J_GM_SI
            if (fit_tr[o]):
                for m in range(0, nradius):
                    r_vec[o*nradius+m] = r_vec[o*nradius+m] * \
                        S_radius_SI / J_radius_e_SI

        # Estimate Transmission Spectroscopy Metric (TSM)
        # https://iopscience.iop.org/article/10.1088/1538-3873/aadf6f
        # Check the scaled factor
        if (np.mean(r_vec[o*nradius]) < 1.5):
            TSM = 0.190
        elif (np.mean(r_vec[o*nradius]) > 1.5 and np.mean(r_vec[o*nradius]) < 2.75):
            TSM = 1.26
        elif (np.mean(r_vec[o*nradius]) > 2.75 and np.mean(r_vec[o*nradius]) < 4.0):
            TSM = 1.28
        else:
            TSM = 1.15
        #
        TSM = TSM * r_vec[o*nradius]**3*Teq_vec[o] / \
            m_vec[o] / rstar**2 * 10.**(-mag_j/5.)

        # Print the parameters
        # Fitted parameters
        pl = plabels[o]
        opars.write(
            '--------------------------------------------------------------\n')
        opars.write('                   Parameters %s\n' %
                    (star + ' ' + plabels[o]))
        opars.write(
            '-------------------------Fitted-------------------------------\n')
        print_values(T0_vec[o], 'T0', 'Tzero'+pl, 'days', 'days')
        print_values(P_vec[o], 'P', 'P'+pl, 'days', 'days')
        if (is_ew):
            print_values(e_dum, 'ew 1', 'esin'+pl, ' ', ' ')
            print_values(w_dum, 'ew 2', 'ecos'+pl, ' ', ' ')
        else:
            print_values(e_vec[o], 'e', 'e'+pl, ' ', ' ')
            print_values(w_vec[o]*180./np.pi, 'w', 'w'+pl, 'deg', 'deg')
        if (fit_tr[o]):
            print_values(b_vec[o], 'b', 'b'+pl, ' ', ' ')
            if (sample_stellar_density):
                if (o == 0):
                    print_values(params[base+5], 'rho*', 'dentrhee'+pl,
                                 'g^{1/3}/cm', '${\\rm g^{1/3}\,cm^{-1}}$')
            else:
                print_values(ar_vec[o], 'a/R*', 'ar'+pl, ' ', ' ')
            for m in range(0, nradius):
                print_values(rr_vec[o*nradius+m], 'rp/R*' +
                             bands[m], 'rr'+pl+bands[m], ' ', ' ')
        if (fit_rv[o]):
            print_values(k_vec[o]*1e3, 'K', 'k'+pl,
                         'm/s', '${\\rm m\,s^{-1}}$')
        opars.write(
            '-------------------------Derived------------------------------\n')
        if (fit_rv[o]):
            print_values(m_vec[o], 'Mp', 'mp'+pl, 'M_' +
                         unit_mass, '$M_'+usymbol+'$')
        if (fit_tr[o]):
            for m in range(0, nradius):
                print_values(r_vec[o*nradius+m], 'Rp'+bands[m],
                             'rp'+pl+bands[m], 'R_'+unit_mass, '$R_'+usymbol+'$')
        print_values(Tpe_vec[o], 'Tperi', 'Tperi'+pl, 'days', 'days')
        if (is_ew):
            print_values(e_vec[o], 'e', 'e'+pl, ' ', ' ')
            print_values(w_vec[o]*180./np.pi, 'w', 'w'+pl, 'deg', 'deg')
        if (fit_tr[o]):
            print_values(i_vec[o]*180./np.pi, 'i', 'i'+pl, 'deg', 'deg')
            if (sample_stellar_density):
                print_values(ar_vec[o], 'a/R*', 'ar'+pl, ' ', ' ')
            print_values(a_vec[o], 'a', 'a'+pl, 'AU', 'AU')
            for m in range(0, nradius):
                print_values(
                    1.e6*rr_vec[o*nradius+m]**2, 'depth'+bands[m], 'depth'+pl+bands[m], 'ppm', 'ppm')
                print_values((np.sqrt(1.-b_vec[o]**2)*rr_vec[o*nradius+m]**2*vsini)
                             * 1.e3, 'RM'+bands[m], 'RM'+pl+bands[m], 'm/s', '${\\rm m\,s^{-1}}$')
            print_values(Fp, 'Insolation', 'insolation'+pl,
                         'F_Earth', '${\\rm F_{\\oplus}}$')
            print_values(TSM, 'TSM', 'tsm'+pl, 'TSM (Kempton et al 2018)', ' ')
            print_values(ds_vec[o], 'rho*', 'denstr'+pl,
                         'g/cm^3 (transit)', '${\\rm g\,cm^{-3}}$')
            print_values(irho_vec, 'rho*', 'denssp'+pl,
                         'g/cm^3 (stellar paramters)', '${\\rm g\,cm^{-3}}$')
            print_values(Teq_vec[o], 'Teq', 'Teq'+pl, 'K (albedo=0)', 'K')
            print_values(trt_vec[o], 'T_tot', 'ttot'+pl, 'hours', 'hours')
            print_values(tri_vec[o], 'T_full', 'tful'+pl, 'hours', 'hours')
            print_values((trt_vec[o] - tri_vec[o])/2.0,
                         'T_in/eg', 'teg'+pl, 'hours', 'hours')
            print_values(delta_mag_vec,
                         'Delta mag', 'deltamag'+pl, 'eq. (4) Vanderburg et al., 2019', '')
            if is_single_transit:
                print_values(tr_dur_vec*24., 'T_duration',
                             't_dur'+pl, 'hours', 'hours')
                print_values(
                    P_circ_vec, 'Possible period (if curcular orbit)', 'p_cir'+pl, 'days', 'days')
        if fit_rv[o] and fit_tr[o]:
            print_values(pden_vec, 'rho_p', 'denp'+pl,
                         'g/cm^3', '${\\rm g\,cm^{-3}}$')
            print_values(pgra_vec, 'g_p', 'grap'+pl,
                         'cm/s^2 (K and Rp/R*)', '${\\rm cm\,s^{-2}}$')
            print_values(pgra_vec2, 'g_p', 'grappars'+pl,
                         'cm/s^2 (planet parameters)', '${\\rm cm\,s^{-2}}$')
            print_values(jsp, 'Jeans', 'jsp'+pl,
                         'Jeans Scape Parameter (Fossati et al 2017)', ' ')
        opars.write(
            '--------------------------------------------------------------\n')

        # Let us change to the next planet
        base = base + 7


# The other parameters
q1_vec = [None]*nbands
q2_vec = [None]*nbands
u1_vec = [None]*nbands
u2_vec = [None]*nbands
for o in range(0, nbands):
    q1_vec[o] = params[4+sldc+2*o]
    q2_vec[o] = params[4+sldc+1+2*o]

    u1_vec[o], u2_vec[o] = pti.get_us(q1_vec[o], q2_vec[o])

rv_vec = [None]*nt
for o in range(0, nt):
    rv_vec[o] = params[4+srv+o]

opars.write('--------------------  Other parameters -----------------------\n')
if (total_tr_fit):
    for o in range(0, nbands):
        print_values(q1_vec[o], 'q1'+bands[o], 'qone'+bands[o], '', '')
        print_values(q2_vec[o], 'q2'+bands[o], 'qtwo'+bands[o], '', '')
        print_values(u1_vec[o], 'u1'+bands[o], 'uone'+bands[o], '', '')
        print_values(u2_vec[o], 'u2'+bands[o], 'utwo'+bands[o], '', '')
if (total_rv_fit):
    for o in range(0, nt):
        print_values(rv_vec[o], 'Sys. vel. '+telescopes_labels[o],
                     telescopes_labels[o], 'km/s', '${\\rm km\,s^{-1}}$')
opars.write('--------------------------------------------------------------\n')


if (is_jitter_rv or is_jitter_tr):
    if (total_rv_fit):
        for o in range(0, n_jrv):
            print_values(params[4+sjitrv+o]*1.e3, telescopes_labels[o] +
                         ' jitter', 'j'+telescopes_labels[o], 'm/s', '${\\rm m\,s^{-1}}$')
    if (total_tr_fit):
        for o in range(0, n_jtr):
            print_values(params[4+sjittr+o], 'tr jitter' +
                         bands[o], 'jtr'+bands[o], '', '')
    opars.write(
        '--------------------------------------------------------------\n')

if (is_linear_trend != 'f' or is_quadratic_trend != 'f'):
    if (total_rv_fit):
        print_values(params[4+strends]*1.e3, 'linear trend',
                     'ltrend', 'm/s/days', '${\\rm m\,s^{-1}\,d^{-1}}$')
        print_values(params[4+strends+1]*1.e3, 'quadratic trend',
                     'qtrend', 'm/s/days^2', '${\\rm m\,s^{-1}\,d^{-2}}$')
    opars.write(
        '--------------------------------------------------------------\n')
opars.write('\n')


# Hyper-parameters
# RV
for o in range(0, np_rv):
    print_values(params[4+skrv+o], krv_labels[o], 'j'+krv_labels[o], '', '')

# TR
for o in range(0, np_tr):
    print_values(params[4+sktr+o], ktr_labels[o], 'j'+ktr_labels[o], '', '')


# RESIZE TRANSIT ERROR BARS
# if ( is_jitter_tr and resize_tr ):
#  jit_tr = best_value(params_jitter[n_jrv],maxloglike,get_value)
#  for o in range(0,len(et)):
#      for m in range(0,len(et[o])):
#              et[o][m] = np.sqrt(et[o][m]**2 + jit_tr**2)
#  for o in range(0,len(lc_errs)):
#    lc_errs[o] = np.sqrt( lc_errs[o]**2 + jit_tr**2)

if (total_rv_fit):
    new_errs_all = [None]*len(errs_all)
    for o in range(0, len(errs_all)):
        new_errs_all[o] = list(errs_all[o])

if (is_jitter_rv and resize_rv):
    for j in range(0, n_jrv):
        jit_rv = best_value(params[4+sjitrv+j], maxloglike, get_value)
        for o in range(0, len(errs_all[j])):
            new_errs_all[j][o] = 1.e3*np.sqrt(errs_all[j][o]**2 + jit_rv**2)

opars.close()
otex.close()

# Print the output in the screen
with open(out_params_file) as dummy_file:
    for l in dummy_file:
        print(l)
