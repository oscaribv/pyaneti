# -----------------------------------------------------------------
#       TRANSIT PARAMERS TO BE USED TO GENERATE PLOTS
# -----------------------------------------------------------------

# Create parameters vector
pars_tr = np.zeros(shape=(nplanets, 6))
for m in range(0, nplanets):
    pars_tr[m][0] = tp_val[m]
    pars_tr[m][1] = P_val[m]
    pars_tr[m][2] = e_val[m]
    pars_tr[m][3] = w_val[m]
    pars_tr[m][4] = i_val[m]
    pars_tr[m][5] = a_val[m]


tbin = 1.

# ===========================================================
# ===========================================================

def create_folded_tr_plots():

    for o in range(nplanets):

        tr_vector = [None]*nbands

        if (fit_tr[o]):

            for m in range(nbands):

                #compute the models only for the current label
                indices = m == np.asarray(trlab)
                localx = lc_time[indices]
                localy = lc_flux[indices]
                locale = lc_errs[indices]
                localt = np.zeros(len(localx),dtype=int)

                try:
                    tr_vector[m] = create_tr_vector(
                        localx, localy, locale, localt, pars_tr, rp_val, o, m)
                except:
                    tr_vector[m] = np.array([0.]),np.array([1.]),np.array([1e-6]), \
                                    np.array([0.]),np.array([1.]),np.array([1.]), \
                                    np.array([0.]),np.array([1.])

            #Transposing the list to feed fancy_tr_plot as requiered
            #note that trasposing as numpy arrays was creating some issues
            #solved using https://stackoverflow.com/questions/6473679/transpose-list-of-lists
            transpose_tr = list(map(list, zip(*tr_vector)))

            fancy_tr_plot(transpose_tr, o)

# these lists contains all the information for a given label
def fancy_tr_plot(tr_vector, pnumber):

    fname = outdir+'/'+star+plabels[pnumber]+'_tr.pdf'
    print('Creating ', fname)

    # Do the plot
    tfc = 24.  # time factor conversion to hours
    local_T0 = 0.

    # Extract the vectors to be plotted from tr_vector
    xtime = tr_vector[0]
    yflux = tr_vector[1]
    eflux = tr_vector[2]
    xmodel = tr_vector[3]
    xmodel_res = tr_vector[4]
    fd_ub = tr_vector[5]
    res_res = tr_vector[6]
    fd_ub_unbinned = tr_vector[7]

    # Start the plot
    plt.figure(1, figsize=(fsx, fsy+(nbands-1)*0.75*fsy))
    # Plot the transit light curve
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[
                           3.0, 1./(1+(nbands-1)*0.75)])
    gs.update(hspace=0.00)
    ax1 = plt.subplot(gs[0],rasterized=is_rasterized)
    plt.tick_params(labelsize=fos, direction='in')
    x_lim = (min(np.concatenate(xtime))-local_T0)*tfc
    plt.xlim(x_lim, -x_lim)
    if (select_y_tr):
        plt.ylim(y_lim_min, y_lim_max)
    min_val_model = max(np.concatenate(fd_ub)) - min(np.concatenate(fd_ub))
    deltay = 0.
    dy = 0.
    for m in range(0, nbands):
        if (plot_tr_errorbars):
            plt.errorbar((xtime[m]-local_T0)*tfc, yflux, errors,
                         color=tr_colors[m], marker='.', alpha=1.0)
        else:
            local_color = tr_colors[m]
            if plot_binned_data:
                local_color = '#C0C0C0'
            plt.plot((xtime[m]-local_T0)*tfc, yflux[m]-deltay, color=local_color,
                     ms=7, marker=mark_tr[m], alpha=0.75, linewidth=0.)
            plt.plot((xmodel[m]-local_T0)*tfc, fd_ub[m] -
                     deltay, 'k', linewidth=2.0, alpha=1.0)
            plt.errorbar(-x_lim*(0.95), min(yflux[m]-deltay), eflux[m]
                         [0], color=local_color, ms=7, marker=mark_tr[m], alpha=1.0)
            plt.annotate(
                'Error bar '+bands[m], xy=(-x_lim*(0.70), min(yflux[m]-deltay)), fontsize=fos*0.7)
            if plot_binned_data:
                xbined, fbined, rbined = bin_data(xtime[m]*tfc,yflux[m],res_res[m]*1e6,tbin)
                plt.plot(xbined, fbined-deltay, 'o',color=tr_colors[m])
            if (m < nbands - 1):
                dy = max(yflux[m+1])-min(yflux[m+1])
            deltay = deltay + dy

        # save the data
        model_matrxi = np.asarray([(xmodel[m]-local_T0)*tfc, fd_ub[m]])
        np.savetxt(
            outdir+'/'+star+plabels[pnumber]+'-trmodel'+bands[m]+'.txt', model_matrxi.T)
        data_matrxi = np.asarray(
            [(xtime[m]-local_T0)*tfc, yflux[m], eflux[m], res_res[m]])
        np.savetxt(
            outdir+'/'+star+plabels[pnumber]+'-trdata'+bands[m]+'.txt', data_matrxi.T)

    if (nbands == 1):
        plt.ylabel('Flux', fontsize=fos)
    if (nbands > 1):
        plt.ylabel('Flux + offset', fontsize=fos)
    # Calculate the optimal step for the plot
    step_plot = int(abs(x_lim))  # the value of the x_axis
    # now we ensure the result is par
    step_plot = step_plot + int(step_plot % 2)
    step_plot = int(step_plot / 8.) + 1  # The size of the jump depends
    # let us get the new limit
    nuevo = np.arange(0, int(abs(x_lim)) + step_plot, step_plot)
    mxv = np.max(nuevo)
#  plt.xticks( np.arange(int(x_lim),int(-x_lim)+1,step_plot))
    plt.xticks(np.arange(-mxv, mxv+step_plot, step_plot))
    plt.minorticks_on()
    plt.ticklabel_format(useOffset=False, axis='y')
    plt.xlim(x_lim, -x_lim)
    plt.tick_params(axis='x', which='both', direction='in', labelbottom=False)
    plt.tick_params(axis='y', which='both', direction='in')
    # ------------------------------------------------------------
    # Plot the residuals
    # ------------------------------------------------------------
    ax0 = plt.subplot(gs[1],rasterized=is_rasterized)
    plt.tick_params(labelsize=fos, direction='in')
    deltay = 0.
    dy = 0.
    for m in range(nbands):
        if (plot_tr_errorbars):
            plt.errorbar((xmodel_res[m]-local_T0)*tfc, res_res[m]
                         * 1e6, eflux[m]*1e6, fmt='.', alpha=0.5)
        else:
            local_color = tr_colors[m]
            if plot_binned_data:
                local_color = '#C0C0C0'
            plt.plot((xmodel_res[m]-local_T0)*tfc, res_res[m]*1e6-deltay,mark_tr[m],
                     color=local_color, ms=7, alpha=0.5)
            if plot_binned_data:
                xbined, fbined, rbined = bin_data(xtime[m]*tfc,yflux[m],res_res[m]*1e6,tbin)
                plt.plot(xbined,rbined-deltay,'o',color=tr_colors[m])
            if (m < nbands - 1):
                dy = max(res_res[m])-min(res_res[m])
                dy = dy*1e6
            deltay = deltay + dy
    #plt.plot([x_lim, -x_lim], [0.0, 0.0], 'k--', linewidth=1.0, alpha=1.0)
    plt.xticks(np.arange(-mxv, mxv+step_plot, step_plot))
    plt.xlim(x_lim, -x_lim)
    yylims = ax0.get_ylim()
    miy = (max(abs(yylims[0]), abs(yylims[1])))
    plt.yticks(range(-int(miy), int(miy), int(miy/2.)))
    #plt.ylim(-miy, miy)
    # Calcualte the rms
    if (is_plot_std_tr):
        trsigma = np.std(res_res*1e6, ddof=1)
        trsstr = ('%4.0f ppm' % (trsigma))
        y0, yyyy = ax0.get_ylim()
        plt.annotate('$\sigma = $'+trsstr, xy=(x_lim *
                                               (0.80), y0 + 1.8*miy), fontsize=fos*0.7)
#  if ( select_y_tr ):
#    plt.ylim( - ( y_lim_max - 1.0),y_lim_max - 1.0 )
    # Plot the residuals
    plt.minorticks_on()
    plt.tick_params(axis='x', which='both', direction='in')
    plt.tick_params(axis='y', which='both', direction='in')
    plt.ylabel('Residuals (ppm)', fontsize=fos*0.75)
    plt.xlabel("T - T0 (hours)", fontsize=fos)
    plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=300)
    plt.savefig(fname[:-3]+'png', format='png', bbox_inches='tight', dpi=300)
    plt.close()


def create_tr_vector(time, flujo, eflujo, trlab, pars_tr, rp, plabel, bandlab):


    time = np.asarray(time)
    flujo = np.asarray(flujo)
    eflujo = np.asarray(eflujo)

    P = best_value(P_vec[plabel], maxloglike, get_value)
    T0 = best_value(T0_vec[plabel], maxloglike, get_value)
    tt = best_value(trt_vec[plabel], maxloglike, get_value)
    tt = tt/24.

    #span has to be given in units of days
    span = 2*tt
    if ( span_tr[0] > 0.):
        span = span_tr[plabel]

    #indices = []
    phase = abs(((time-T0)%P)/P)
    phase[phase>0.5] -= 1
    indices = abs(phase) <= span/P

    #Let us use only data that appears close to the transits
    #Local time is useful to compute the transits of other planets in the light curve
    local_time = time[indices]
    #xtime, yflux, and e flux contain the phase, flux and e flux of the data
    xtime = phase[indices]*P
    yflux = flujo[indices]
    eflux = eflujo[indices]

    xmodel_res = np.asarray(xtime)
    mimax = span
    xmodel = np.arange(-mimax, mimax, 1.0/60./24.)
    newtrlab = [0]*len(xmodel)

    # Let us create the model
    control = 1
    if nradius == 1:
        control = 0
    # The model has T0 = 0
    dumtp = pti.find_tp(0.0, e_val[plabel], w_val[plabel], P_val[plabel])
    tr_p = np.concatenate([[dumtp], pars_tr[plabel][1:]])
    flux_model = pti.flux_tr(xmodel, newtrlab, tr_p, rp_val[plabel*nradius+bandlab*control],
                        my_ldc[bandlab*2:bandlab*2+2], n_cad[bandlab], t_cad[bandlab], nradius=1)
    # Let us create an unbinned model plot
    flux_model_unbinned = pti.flux_tr(
        xmodel, newtrlab, tr_p, rp_val[plabel*nradius+bandlab*control], my_ldc[bandlab*2:bandlab*2+2], [1], t_cad[bandlab], nradius=1)
    # Calculate the flux to copute the residuals
    newtrlab = [0]*len(xmodel_res)
    flux_model_res = pti.flux_tr(xmodel_res, newtrlab, tr_p, rp_val[plabel*nradius+bandlab*control],
                            my_ldc[bandlab*2:bandlab*2+2], n_cad[bandlab], t_cad[bandlab], nradius=1)


   #############################################################################
   # Let us calculate the flux caused by the other planets
    # Define a vector which will contain the data of other planers for multi fits
    flux_other_planets = np.ones(shape=len(flux_model_res))
    for p in range(nplanets):
        if (p != plabel):
            # flux_other_planets stores the flux of a star for each independent
            newtrlab = [0]*len(local_time)
            flux_other_planets = flux_other_planets * pti.flux_tr(local_time, newtrlab, pars_tr[p], rp_val[p*nradius+bandlab*control],
                                                    my_ldc[bandlab*2:bandlab*2+2], n_cad[bandlab], t_cad[bandlab], nradius=1)

    # Remove extra planets from the data
    yflux_local = yflux / flux_other_planets
    # The flux has been corrected for the other planets

    # Get the residuals
    res_res = yflux_local - flux_model_res
    # are we plotting a GP together with the RV curve
    if kernel_tr[0:2] != 'No':
        xvec = local_time
        yvec = res_res
        evec = eflux
        m, C = pti.pred_gp(kernel_tr, pk_tr, xvec, yvec, evec,
                           xvec, jtr[bandlab], [0]*len(xvec))
        yflux_local = yflux_local - m
        res_res = res_res - m

    return xtime, yflux_local, eflux, xmodel, xmodel_res, flux_model, res_res, flux_model_unbinned

# ===========================================================
#              plot all transits
# ===========================================================

# Now this functions works only with one band


def plot_all_transits():
    global plot_tr_errorbars

    # Create the plot of the whole light
    model_flux = pti.flux_tr(
        lc_time, [0]*len(lc_time), pars_tr, rp_val, my_ldc, n_cad, t_cad)
    res_flux = lc_flux - model_flux

    for i in range(0, nplanets):
        if (fit_tr[i]):
            if (m < nbands - 1):
                dy = max(yflux[m+1])-min(yflux[m+1])
            if (m < nbands - 1):
                dy = max(yflux[m+1])-min(yflux[m+1])

            xt, dt, yt, et = create_transit_data(
                lc_time, lc_flux, lc_errs, i, span_tr[i])
            xt2, dt2, yt2, et2 = create_transit_data(
                lc_time, res_flux, lc_errs, i, span_tr[i])

            if (is_plot_all_tr[i]):
                for j in range(0, len(xt)):
                    xtm = np.arange(min(xt[j]), max(xt[j]), 1./20./24.)
                    ytm = pti.flux_tr(xtm, trlab, pars_tr,
                                      rp_val, my_ldc, n_cad, t_cad)

                    fname = outdir+'/'+star+plabels[i]+'_transit'+str(j)+'.pdf'
                    n = xt[j][len(xt[j])-1] - xt[0][0]
                    n = int(n/P_val[i])
                    #is_err = plot_tr_errorbars
                    #plot_tr_errorbars = True
                    fancy_tr_plot(t0_val[i]+P_val[i]*n, xt[j], yt[j], et[j],
                                  xtm, xt2[j], ytm, np.array(yt2[j]), ytm, fname)
                    #plot_tr_errorbars = is_err


def plot_lightcurve_timeseries():
    ''' This function creates a light curve time-series plot '''

    # Here I need to create a special trlab in order to separate the different colors
    # Now let us imagine that it works with 1-band
    xmodel = np.arange(min(lc_time), max(lc_time), 5./60./24.)
    my_trlab = [0]*len(xmodel)
    ymodel = pti.flux_tr(xmodel, my_trlab, pars_tr.transpose(),
                         rp_val, my_ldc[0:2], n_cad[0], t_cad[0],nradius)
    #ymodel = pti.flux_tr(xmodel,my_trlab,pars_tr,rp_val,my_ldc,n_cad,t_cad)

    # Calcualte the residuals
    trres = pti.flux_tr(lc_time, trlab, pars_tr.transpose(),
                        rp_val, my_ldc, n_cad, t_cad, nradius)
    trres = lc_flux - trres

    # are we plotting a GP together with the RV curve
    if kernel_tr[0:2] != 'No':
        xvec = lc_time
        yvec = lc_flux
        evec = lc_errs
        m, C = pti.pred_gp(kernel_tr, pk_tr, xvec, trres,
                           evec, xmodel, jtr, jtrlab)
        tr_mvec = [xmodel, ymodel, 1.+m, (ymodel+m)]
        model_labels = ['Planetary signal', 'GP', 'P+GP']
        tr_mvec = np.array(tr_mvec)
        mcolors = ['r', 'b', 'k']
        malpha = [0.7, 0.7, 0.9]
    else:
        tr_mvec = np.zeros(shape=(2,len(xmodel)))
        tr_mvec[0] = xmodel
        tr_mvec[1] = ymodel
        model_labels = ['Planetary signal']
        mcolors = ['k']
        malpha = [1.]


    indices = abs(trres) < ( np.mean(trres) + 3.5*np.std(trres) )
    x_c = lc_time[indices]
    y_c = lc_flux[indices]
    e_c = lc_errs[indices]
    r_c = trres[indices]
    local_trlab = np.asarray(trlab)
    t_c = local_trlab[indices]

    tr_dvec_file = np.asarray([x_c, y_c, e_c, r_c, t_c],dtype=object)

    # save the data
    header_1 = 'Time  Flux_model'
    np.savetxt(
        outdir+'/'+star+'-trmodel_lightcurve.txt',tr_mvec.T,header=header_1,fmt='%1.7e')
    header_2 = 'This light curve has been cleaned with a 3.5-sigma clipling \n Time  Flux eFlux residuals telescope_label'
    np.savetxt(
         outdir+'/'+star+'-trdata_lightcurve.txt',tr_dvec_file.T,header=header_2,fmt='%1.7e %1.7e %1.7e %1.7e %i')


    # Name of plot file
    fname = outdir+'/'+star+'_lightcurve.pdf'
    tr_dvec = np.asarray([lc_time, lc_flux, lc_errs, lc_errs, trres, trlab],dtype=object)
    plot_labels_tr = [rv_xlabel, 'Flux', 'Residuals']
    # Create the RV timeseries plot
    create_nice_plot(tr_mvec, tr_dvec, plot_labels_tr, model_labels, bands, fname,
                     plot_residuals=False, fsx=2*fsx, model_colors=mcolors, model_alpha=malpha,colors=tr_colors)
