# ===========================================================
#                   RV PLOTS
# ===========================================================
# ===========================================================
#        Plot RV timeseries
# ===========================================================


def plot_rv_timeseries():
    # -------------------------------------------------------------------------------
    # Convert factor from km to m
    cfactor = np.float(1.e3)
    data_rv_inst = [None]*nt
    data_erv_inst = [None]*nt
    for i in range(0, nt):
        data_rv_inst[i] = list(rv_all[i])
        data_erv_inst[i] = list(errs_all[i])
# -------------------------------------------------------------------------------
    # CREATE MODEL TO PLOT
    # MODEL LENGHT IS THE WHOLE RANGE +/- 10%
    xmin = min(np.concatenate(time_all))
    xmax = max(np.concatenate(time_all))
    total_tt = int(xmax - xmin)
    add = total_tt*0.1
    xmax = xmax + add
    xmin = xmin - add
    n = total_tt*100
    if (n > 5000):
        n = 5000
    # CREATE TIME VECTOR
    rvx = np.arange(xmin, xmax, (xmax-xmin)/n)
    # COMPUTE MODEL WITH ALL THE PLANETS
    rvy = pti.rv_curve_mp(rvx, 0.0, t0_val, k_val, P_val,
                          e_val, w_val, alpha_val, beta_val)
# -------------------------------------------------------------------------------
    # CORRECT DATA FOR EACH SPECTROGRAPH OFFSET
    rv_no_offset = [None]*nt
    rv_residuals = [None]*nt
    for j in range(0, nt):
        # Remove the compute offset for instrument j
        rv_no_offset[j] = np.asarray(data_rv_inst[j] - v_val[j])
        # Store the model with the planets to compute the residuals for instrument j
        rv_residuals[j] = pti.rv_curve_mp(time_all[j], 0.0, t0_val,
                                          k_val, P_val, e_val, w_val, alpha_val, beta_val)
        # Compute the residuals for the instrument j
        rv_residuals[j] = np.asarray(rv_no_offset[j] - rv_residuals[j])
# -------------------------------------------------------------------------------
    if kernel_rv[0:2] == 'MQ' or kernel_rv[0:2] == 'ME' or kernel_rv[0:2] == 'MM' or  kernel_rv[0:2] == 'SQ':
        # How many timeseries do we have?
        #ns = int((len(fit_krv) - 3)/2)
        ns = int(kernel_rv[2])
        # This vector contains all the timeseries (a TxN vector)
        xvec = rv_time
        # This vector contains the predicted timeseries for N cases
        rvx_tot = np.concatenate([rvx]*ns)
        # Let us create our vector with the residuals for the N timeseries
        yvec = [None]*len(rv_vals)
        # Our first chunk of timeseries is always the RV, in this case, we have to
        # Create the residuals once removed the planet signals
        yvec_noplanet = np.concatenate(rv_no_offset)
        yvec_planet = np.concatenate(rv_residuals)
        #The first elements of the array correspond to the GP, so we want the vector with no planets and no offset
        for i in range(int(len(rv_vals)/ns)):
            yvec[i] = yvec_planet[i]
        # Now, let us store the ancilliary data vectors, so we want to remove only the offsets
        for i in range(int(len(rv_vals)/ns), int(len(rv_vals))):
            yvec[i] = yvec_noplanet[i]
        evec = rv_errs
        # Now, predict the GP for all the timeseries
        m, C = pti.pred_gp(kernel_rv, pk_rv, xvec, yvec,
                           evec, rvx_tot, jrv, jrvlab)
        # Let us remove the GP model from the data
        m_gp, C_gp = pti.pred_gp(
            kernel_rv, pk_rv, xvec, yvec, evec, xvec, jrv, jrvlab)
        #Now we can remove the GP model to the data
        yvec = yvec - m_gp
        # -----------------------------------------------------------------------------------
        # Let us create random samples of data
        if False:
            nsamples = 1000
            for j in range(0, nsamples):
                # Create the Gaussian samples
                y_dummy = np.random.normal(yvec, evec, len(yvec))
                out_f = outdir+'/'+star+'_rv_'+str(j)+'.dat'
                opars = open(out_f, 'w')
                for k in range(0, len(yvec)):
                    opars.write('%4.7f %4.7f %4.7f %s \n' % (
                        xvec[k], y_dummy[k], evec[k], telescopes[tlab[k]]))
                opars.close()
        # -----------------------------------------------------------------------------------
        #
        # Let us create the vectors that we will use for the plots
        plot_vector = [None]*ns
        nts = len(rvx)
        # This corresponds to the RV timeseries
        plot_vector[0] = [rvx, (rvy)*cfactor, m[0:nts]
                          * cfactor, (rvy+m[0:nts])*cfactor,np.sqrt(np.matrix.diagonal(C[0:nts,0:nts]))*cfactor]
        for i in range(1, ns):
            plot_vector[i] = [rvx, m[i*nts:(i+1)*nts]*cfactor,np.sqrt(np.matrix.diagonal(C[i*nts:(i+1)*nts,i*nts:(i+1)*nts]))*cfactor]
    # are we plotting a GP together with the RV curve
    elif kernel_rv[0:2] != 'No':
        xvec = rv_time
        yvec = np.concatenate(rv_residuals)
        evec = rv_errs
        m, C = pti.pred_gp(kernel_rv, pk_rv, xvec,
                           yvec, evec, rvx, jrv, jrvlab)
        rv_mvec = [rvx, (rvy)*cfactor, m*cfactor, (rvy+m)*cfactor]
        model_labels = ['Planetary signal', 'GP', 'P+GP']
        mcolors = ['r', 'b', 'k']
        malpha = [0.7, 0.7, 0.9]
    else:
        rv_mvec = [rvx, (rvy)*cfactor]
        model_labels = ['Full model']
        mcolors = ['k']
        malpha = [1.]

    # Name of plot file
    fname = outdir+'/'+star+'_rv_timeseries.pdf'
    # Name of residuals file
    out_f = outdir+'/'+star+'_rv_residuals.dat'

    # Create the vectors to be used in create_nice_plot()
    vec_x = np.concatenate(time_all)
    #yvec_noplanet = np.concatenate(rv_residuals)
    vec_y = np.concatenate(rv_residuals)
    if kernel_rv[0:2] == 'MQ' or kernel_rv[0:2] == 'ME' or kernel_rv[0:2] == 'MM' or kernel_rv[0:2] == 'SQ':
        vec_y = np.array(yvec)
    vec_z = np.concatenate(new_errs_all)/cfactor
    #
    xdata = vec_x
    ydata = np.asarray(np.concatenate(rv_no_offset))*cfactor
    edata = np.asarray(np.concatenate(errs_all))*cfactor
    ejdata = np.asarray(np.concatenate(new_errs_all))
    res = np.asarray(vec_y)*cfactor

    if kernel_rv[0:2] == 'MQ' or kernel_rv[0:2] == 'ME' or kernel_rv[0:2] == 'MM' or kernel_rv[0:2] == 'SQ' :
        ns = int(kernel_rv[2])
        ts_len = int(len(xdata)/ns)
        # Create the vector with the data needed in the create_nice_plot function
        for o in range(0, ns):
            rv_dvec = [xdata[o*ts_len:(o+1)*ts_len], ydata[o*ts_len:(o+1)*ts_len], edata[o*ts_len:(o+1)*ts_len],
                       ejdata[o*ts_len:(o+1)*ts_len], res[o*ts_len:(o+1)*ts_len], tlab[o*ts_len:(o+1)*ts_len]]
            rv_dvecnp = np.asarray(rv_dvec)
            mvec = plot_vector[o]
            mvecnp = np.asarray(mvec)
            np.savetxt(outdir+'/timeseries_model'+str(o) +
                       '.dat', mvecnp.T, fmt='%8.8f')
            np.savetxt(outdir+'/timeseries_data'+str(o)+'.dat', rv_dvecnp.T, fmt='%8.8f %8.8f %8.8f %8.8f %8.8f %i',
                       header='time rv erv jitter rvnoplanet tlab', comments='#')
            if o == 0:
                model_labels = ['Planetary signal', 'GP', 'P+GP']
                mcolors = ['r', 'b', 'k']
                malpha = [0.7, 0.7, 0.9]
                fname = outdir+'/'+star+'_rv_timeseries.pdf'
            else:
                model_labels = ['GP timeseries'+str(o+1)]
                mcolors = ['k']
                malpha = [0.9]
                fname = outdir+'/'+star+'_timeseries'+str(o+1)+'.pdf'
            if (len(rv_labels) == 1):
                plot_labels_rv = [rv_xlabel, 'RV (m/s)', 'Residuals (m/s)']
            else:
                plot_labels_rv = [
                    rv_xlabel[o*ns:(o+1)*ns], rv_labels[o], rv_res[o]]
            #The mvec[:-1] is to ignore the extra dimension added to create the variance of the P
            create_nice_plot(mvec[:-1], rv_dvec, plot_labels_rv, model_labels, telescopes_labels, fname, std_model=mvec[-1],
                             plot_residuals=False, fsx=2*fsx, model_colors=mcolors, model_alpha=malpha,colors=rv_colors)
    else:
        rv_dvec = [xdata, ydata, edata, ejdata, res, tlab]
        rv_dvecnp = np.asarray(rv_dvec)
        mvecnp = np.asarray(rv_mvec)
        np.savetxt(outdir+'/timeseries_model_rv.dat', mvecnp.T, fmt='%8.8f')
        np.savetxt(outdir+'/timeseries_data_rv.dat', rv_dvecnp.T, fmt='%8.8f %8.8f %8.8f %8.8f %8.8f %i',
                       header='time rv erv jitter rvnoplanet tlab', comments='#')
        plot_labels_rv = [rv_xlabel, 'RV (m/s)', 'Residuals (m/s)']
        # Create the RV timeseries plot
        create_nice_plot(rv_mvec, rv_dvec, plot_labels_rv, model_labels, telescopes_labels, fname,
                         plot_residuals=False, fsx=2*fsx, model_colors=mcolors, model_alpha=malpha,colors=rv_colors)

    # Create residuals file
    of = open(out_f, 'w')
    for i in range(0, len(vec_x)):
        of.write(' %8.8f   %8.8f  %8.8f  %s \n' % (vec_x[i], res[i]*1e-3, vec_z[i], telescopes_labels[tlab[i]]))

    of.close()

    return rv_residuals

# ===========================================================
#                RV multi-planet fit
# ===========================================================


def plot_rv_phasefolded():

    cfactor = np.float(1.e3)

    for i in range(0, nplanets):

        # Create the RV fitted model for the planet i
        rvx = np.arange(t0_val[i], t0_val[i]+P_val[i]*0.999, P_val[i]/4999.)
        rvy = pti.rv_curve_mp(
            rvx, 0.0, t0_val[i], k_val[i], P_val[i], e_val[i], w_val[i], 0.0, 0.0)
        # rvx and rvy are the model timeseries for planet i

        #Let us compute the shadow region for the RV plots
        rv_std = []
        if plot_rv_std:
            #Compute 1000 random models from the samples
            rvy_vec = [None]*1000
            for j in range(1000):
                my_j = np.random.randint(len(T0_vec[i]))
                rvy_vec[j] = pti.rv_curve_mp(
                    rvx, 0.0, T0_vec[i][j], k_vec[i][j], P_vec[i][j], e_vec[i][j], w_vec[i][j], 0.0, 0.0)
            #Compute the standard deviation of the models sample
            rv_std = np.std(rvy_vec,axis=0)
            rv_std *= cfactor



        # Now it is time to remove the planets j != i from the data
        rv_pi = pti.rv_curve_mp(
            rv_time, 0.0, t0_val[i], k_val[i], P_val[i], e_val[i], w_val[i], 0., 0.)
        # This variable has all the planets
        rv_pall = pti.rv_curve_mp(
            rv_time, 0.0, t0_val, k_val, P_val, e_val, w_val, 0.0, 0.0)

        # Let us remove all the signals from the data
        res = np.zeros(shape=len(rv_vals))
        for m in range(0, len(rv_vals)):
            res[m] = rv_vals[m] - v_val[tlab[m]] - rv_pall[m]  \
                - alpha_val*(rv_time[m] - t0_val[0]) - \
                beta_val*(rv_time[m] - t0_val[0])**2
        # Let us add the  signal of the planet i to the data
        rv_planet_i = res + rv_pi

        # Did we fit for a GP?
        evec = np.asarray(rv_errs)
        if kernel_rv[0:2] != 'No':
            xvec = rv_time
            yvec = res
            kernel_val, C = pti.pred_gp(
                kernel_rv, pk_rv, xvec, yvec, evec, xvec, jrv, jrvlab)
            res = res - kernel_val
            rv_planet_i = rv_planet_i - kernel_val

        rvy = np.asarray(rvy)*cfactor
        res = np.asarray(res)*cfactor
        rv_planet_i = np.asarray(rv_planet_i)*cfactor
        evec = evec*cfactor
        ejvec = np.concatenate(new_errs_all)

        p_rv = scale_period(rvx, t0_val[i], P_val[i])
        p_all = scale_period(rv_time, t0_val[i], P_val[i])

        fname = outdir+'/'+star+plabels[i]+'_rv.pdf'
#      plot_rv_fancy([p_rv,rvy,p_all,rv_planet_i,evec,ejvec,res,tlab],fname)
        rv_dvec = np.array([p_all, rv_planet_i, evec, ejvec, res, tlab])
        tellabs = telescopes_labels
        if kernel_rv[0:2] == 'MQ' or kernel_rv[0:2] == 'ME' or kernel_rv[0:2] == 'MM' or kernel_rv[0:2] == 'SQ':
            ns = int(kernel_rv[2])
            nd = int(len(p_all)/ns)
            rv_dvec = [p_all[0:nd], rv_planet_i[0:nd],
                       evec[0:nd], ejvec[0:nd], res[0:nd], tlab[0:nd]]
            tellabs = [0]
            tellabs = telescopes_labels[0:int((len(telescopes_labels))/ns)]
        rv_dvec = np.array(rv_dvec)
        rv_mvec = np.array([p_rv, rvy])
        plot_labels_rv = ['Orbital phase', 'RV (m/s)', 'Residuals (m/s)']
        model_labels = ['']
        np.savetxt(fname[:-4]+'-data.dat',rv_dvec.T,header='phase rv_planet'+plabels[i]+'(m/s) eRV(m/s)  eRV_with_jitter(m/s)  residuals(m/s)   instrument',
        fmt='%4.8f  %4.8f  %4.8f  %4.8f  %4.8f  %i')
        np.savetxt(fname[:-4]+'-model.dat',rv_mvec.T,header='phase rv_planet'+plabels[i]+'(m/s)',fmt='%4.8f  %8.8f')
        create_nice_plot(rv_mvec, rv_dvec, plot_labels_rv,
                         model_labels, tellabs, fname,colors=rv_colors,std_model=rv_std)
