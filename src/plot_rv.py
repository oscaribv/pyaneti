import scipy.stats as stats

#--------------------------------------------------------------------------------
#                  Radial velocity timeseries
#--------------------------------------------------------------------------------

def compute_full_rv_model(x,df_vec,planet_labels=plabels[:nplanets]):

    #Get the parameters
    f = np.zeros(shape=(len(x)))
    for p in planet_labels:
        f += pti.rv_curve_mp(x,0.0, df_vec[f't0{p}'],df_vec[f'p{p}'],df_vec[f'e{p}'],df_vec[f'omega{p}'],df_vec[f'k{p}'],0,0)
    f += (x-x[0])*df_vec['linear_trend'] + (x-x[0])**2*df_vec['quadratic_trend']
    return f

def compute_rv_model_notrends(x,df_vec,planet_labels=plabels[:nplanets]):

    #Get the parameters
    f = np.zeros(shape=(len(x)))
    for p in planet_labels:
        f += pti.rv_curve_mp(x,0.0, df_vec[f't0{p}'],df_vec[f'p{p}'],df_vec[f'e{p}'],df_vec[f'omega{p}'],df_vec[f'k{p}'],0,0)
    return f

def compute_gp_model(x,df_vec,pred_dl):

    jitter_values = np.array(df_vec[rv_jitter_names])

    if kernel_rv[0:2] in ['MQ', 'SQ', 'ME', 'MM']: #multi-GP run
        m, C = pti.pred_gp_multigp(kernel_rv,df_vec[krv_labels],rv_time,rv_noplanets,
                       rv_errs,x,jitter_values,jrvlab,dim_lab,pred_dl)
    else: #Normal GP run
        m, C = pti.pred_gp(kernel_rv,df_vec[krv_labels],rv_time,rv_noplanets,
                       rv_errs,x,jitter_values,jrvlab)

    return m, C


def standardize_data(data, errors):
    """
    Standardize the data by subtracting the mean and dividing by the error bars.
    """
    standardized_data = data #- np.mean(data) / errors
    return standardized_data

def normality_test(data, errors,ftw):
    """
    Perform tests on the standardized data to check for consistency with a standard normal distribution.
    """
    standardized_data = standardize_data(data, errors)

    #ks_statistic, p_value = stats.kstest(data, 'norm')
    #print(f'p_value ks_test {p_value}')
    ks_statistic, p_value = stats.shapiro(data)
    #print(f'p_value sh_test {p_value}')

    #print(f"KS Statistic: {ks_statistic}")
    #print(f"P-value: {p_value}")

    if p_value > 0.05:
        ftw.write(f"#Residuals are consistent with normal with p_value Shapiro test of {p_value:2.6e} \n")
    else:
        ftw.write(f"#Residuals are NOT consistent with normal with p_value Shapiro test of {p_value:2.6e} \n")

    return ks_statistic, p_value

#Let us compute the error bars with jitter term
rv_errs_jit = np.array(rv_errs)
for i,t in enumerate(rv_jitter_labels):
    ix = i == np.array(jrvlab)
    rv_errs_jit[ix] = np.sqrt(df_medians[f'rv_jitter{rv_jitter_labels[i]}']**2 + rv_errs_jit[ix]**2)

# check the rv_colors vector is long enough, lengthen it if not
rv_labels_length = len(telescopes_labels)
if len(rv_colors)<rv_labels_length:
    print('Insufficient RV colours, consider adding more in src/default.py; colours will cycle')
    for i in range(rv_labels_length - len(rv_colors)):
        rv_colors.append(rv_colors[i])

#Let's do plots for the nomical case with no GPs
if kernel_rv == 'None':

    plt.figure(1, figsize=(fsx, fsy/2))


    #Remove offsets
    rv_no_offset = np.array(rv_vals)
    #Let us remove the offsets
    for i,t in enumerate(telescopes_labels):
        ix = i == np.array(tlab)
        rv_no_offset[ix] -= df_medians[f'{t}']


    local_rv = np.array(rv_no_offset)
    local_time = np.array(rv_time)
    local_errs = np.array(rv_errs)

    #Compute planet models
    npoints = int(local_time.max() - local_time.min())*24
    if npoints > 5000: npoints = 5000
    xmodel = np.linspace(min(rv_time)-3, max(rv_time)+3,npoints)
    rvmodel = compute_full_rv_model(xmodel,df_medians)
    rvsamples = []
    for i in range(n_samples):
        rvsamples.append(compute_full_rv_model(xmodel,df_samples.iloc[i]))

    #Remove planets
    rv_noplanets = local_rv - compute_rv_model_notrends(local_time,df_medians)
    #residuals to be used in the phase folded plots, remove any remaining trend
    rv_presiduals = rv_noplanets - (local_time-local_time[0])*df_medians['linear_trend'] - (local_time-local_time[0])**2*df_medians['quadratic_trend']

    for i,t in enumerate(telescopes_labels):
        ix = i == np.array(tlab)
        plt.errorbar(local_time[ix],rv_no_offset[ix]*1e3,rv_errs_jit[ix]*1e3,fmt='.',color=rv_colors[i],alpha=0.5,zorder=5)
        plt.errorbar(local_time[ix],rv_no_offset[ix]*1e3,local_errs[ix]*1e3,fmt=mark[i],mfc='w',label=t,color=rv_colors[i],zorder=5)
    plt.plot(xmodel,rvmodel*1e3,'k-')
    for i in range(n_samples):
        plt.plot(xmodel,rvsamples[i]*1e3,'k',alpha=0.05,lw=0.5,zorder=4)
    if is_rv_legend: plt.legend()
    plt.xlim(xmodel.min(),xmodel.max())
    plt.xlabel(rv_xlabel)
    plt.ylabel(r'RV ($\mathrm{m\,s^{-1}}$)')

    #Save the file with the models
    dname = f'{outdir}/timeseries_rv_model.dat'
    with open(dname,'w') as f:
        f.write(f'#Time rv_model   \n')
        for k in range(len(xmodel)):
            f.write(f'{xmodel[k]}  {rvmodel[k]}  \n')

    #Save the file with the data
    dname = f'{outdir}/timeseries_rv_data.dat'
    with open(dname,'w') as f:
        f.write(f'#time rv_no_offset rv_err rv_jit_err rv_res label \n')
        #Perform KS-test to check if the residuals are consistent with a Gaussian distribution
        #ks_statistic, p_value = normality_test(rv_presiduals,rv_errs_jit,f)
        for k in range(len(local_time)):
            f.write(f'{local_time[k]} {rv_no_offset[k]} {local_errs[k]}  {rv_errs_jit[k]} {rv_presiduals[k]} {telescopes_labels[tlab[k]]}  \n')


    fname = outdir+'/'+star+'_rv_timeseries.pdf'
    print(f'Creating  {fname}')
    plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=300)
    plt.savefig(fname[:-3]+'png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

else:

    #Let us find how many dimensions are we dealing with
    ndim = len(set(dim_lab))

    #Remove offsets for all the timeseries
    rv_no_offset = np.array(rv_vals)
    #Let us remove the offsets
    for i,t in enumerate(telescopes_labels):
        ix = i == np.array(tlab)
        rv_no_offset[ix] -= df_medians[f'{t}']

    local_time = np.array(rv_time)
    local_errs = np.array(rv_errs)

    #Remove planets for dimension 0
    id0 = dim_lab == 0
    rv_noplanets = np.array(rv_no_offset)
    rv_noplanets[id0] -= compute_full_rv_model(local_time[id0],df_medians)
    #Now we have a vector of residials with no planets and no offsets for activity indicators

    #Let us compute the GP model at for the observations
    m_data, C_data = compute_gp_model(local_time,df_medians,dim_lab)
    #Now lete us compute the predictive continuous model
    #Compute planet models
    npoints = int(local_time.max() - local_time.min())*24
    if npoints > 5000: npoints = 5000
    xmodel = np.linspace(min(rv_time)-3, max(rv_time)+3,npoints)
    rvmodel = compute_full_rv_model(xmodel,df_medians)
    #Now let us create the time vector that we will pass to the function to compute the multi-GP
    xgp = []
    lgp = []
    for i in range(ndim):
        xgp.extend(xmodel)
        lgp.extend([i]*len(xmodel))
    #time containing timeseries stamps for each dimension xgp = xmodel*ndim
    xgp = np.array(xgp)
    #dimension labels that correspond to the xgp stamps
    lgp = np.array(lgp)
    #Compute the full model
    m_model, C_model = compute_gp_model(xgp,df_medians,lgp)

    #plotting time
    plt.figure(1, figsize=(fsx, ndim*fsy/2))
    gs = gridspec.GridSpec(nrows=ndim, ncols=1)
    gs.update(hspace=0.00)
    ax = [None]*ndim
    for i in range(ndim):
        ax[i] = plt.subplot(gs[i],rasterized=is_rasterized)
        #dimension indices for the data
        idx = dim_lab == i
        #dimension indices for the model
        idm = lgp == i
        #Let us do a plot for each time-series
        #Plot the data, for each telescope
        for j,t in enumerate(telescopes_labels):
            ix = j == np.array(tlab) #Let us extract the indices for the current telescope
            ix = np.logical_and(ix,idx) #let us extract the index of the telescope and the valid dimension index
            if all(ix == False): continue #if there is no telescope in currect dimension we avoid the plot
            ax[i].errorbar(local_time[ix],rv_no_offset[ix]*1e3,rv_errs_jit[ix]*1e3,fmt='.',color=rv_colors[j],alpha=0.5,zorder=5)
            ax[i].errorbar(local_time[ix],rv_no_offset[ix]*1e3,local_errs[ix]*1e3,fmt=mark[j],mfc='w',label=t,color=rv_colors[j],zorder=5)
        ymodel = m_model[idm]
        if i == 0: #We add the planetary model
            ymodel +=  rvmodel
            ax[i].plot(xmodel,rvmodel*1e3,'-',color='r',alpha=1,label='Keplerians')
            rv_presiduals = rv_noplanets[idx] - m_data[idx]
        #Plot the model
        ax[i].plot(xgp[idm],ymodel*1e3,'-',color='k',label='GP model')
        std_model = np.sqrt(np.diagonal(C_model)[idm])
        ax[i].fill_between(xgp[idm],(ymodel-1*std_model)*1e3,(ymodel+1*std_model)*1e3,color='k',alpha=0.1,lw=0,zorder=1)
        ax[i].fill_between(xgp[idm],(ymodel-2*std_model)*1e3,(ymodel+2*std_model)*1e3,color='k',alpha=0.1,lw=0,zorder=1)
        #
        ax[i].set_xlim(xmodel.min(),xmodel.max())
        ax[i].set_ylabel(rv_ylabel[i])
        if is_rv_legend: ax[i].legend()
        if i == ndim - 1:
            ax[i].set_xlabel(rv_xlabel)
        else:
            ax[i].set_xticks([])
        #
        #Save the file with the models
        dname = f'{outdir}/timeseries_dim{i}_model.dat'
        with open(dname,'w') as f:
            if i == 0:
                f.write(f'#Time Full_model std_model GP_model rv_model   \n')
                for k in range(len(xgp[idm])):
                    f.write(f'{xgp[idm][k]} {ymodel[k]} {std_model[k]} {m_model[idm][k]} {rvmodel[k]}  \n')
            else:
                f.write(f'#Time GP_model std_model residuals  \n')
                for k in range(len(xgp[idm])):
                    f.write(f'{xgp[idm][k]} {ymodel[k]} {std_model[k]} \n')
        #Save the file with the data
        dname = f'{outdir}/timeseries_dim{i}_data.dat'
        with open(dname,'w') as f:
            if i == 0:
                f.write(f'#time rv_no_offset rv_err  rv_jit_err rv_activity rv_res label \n')
                #Perform KS-test to check if the residuals are consistent with a Gaussian distribution
                #ks_statistic, p_value = normality_test(rv_presiduals,rv_errs_jit[idx],f)
                for k in range(len(local_time[idx])):
                    f.write(f'{local_time[idx][k]} {rv_no_offset[idx][k]} {local_errs[idx][k]}  {rv_errs_jit[idx][k]} {rv_noplanets[idx][k]} {rv_presiduals[k]} {telescopes_labels[tlab[idx][k]]}  \n')
            else:
                f.write(f'#Time no_offset err jit_err res label \n')
                #Perform KS-test to check if the residuals are consistent with a Gaussian distribution
                local_res = rv_no_offset[idx]-m_data[idx]
                #ks_statistic, p_value = normality_test(local_res,rv_errs_jit[idx],f)
                for k in range(len(local_time[idx])):
                    f.write(f'{local_time[idx][k]} {rv_no_offset[idx][k]} {local_errs[idx][k]}  {rv_errs_jit[idx][k]}  {local_res[k]} {telescopes_labels[tlab[idx][k]]}  \n')

    #
    fname = f'{outdir}/{star}_timeseries_gp.pdf'
    print(f'Creating  {fname}')
    plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=300)
    plt.savefig(fname[:-3]+'png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

    exec(open('src/mgic_rv.py').read())

#--------------------------------------------------------------------------------
#                   Phase folded RV plots
#--------------------------------------------------------------------------------

#Create phase-folded plots for all planets
for k,p in enumerate(plabels[:nplanets]):

    if not fit_rv[k]: continue

    plt.figure(1, figsize=(fsx, fsy))

    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3,1])
    gs.update(hspace=0.00)
    ax0 = plt.subplot(gs[0],rasterized=is_rasterized) #models and data plots
    ax1 = plt.subplot(gs[1],rasterized=is_rasterized) #residuals

    id0 = dim_lab == 0
    #We just computed the residuals with no planets or other signals, let us use to create the plots
    local_rv = np.array(rv_presiduals)
    local_time = np.array(rv_time)[id0]
    local_errs = np.array(rv_errs)[id0]
    local_errs_jit = np.array(rv_errs_jit)[id0]
    local_tel_labs = np.array(telescopes_labels)[tlab[id0]]


    #Compute the model for the current planet
    xmodel = np.linspace(df_medians[f't0{p}'],df_medians[f't0{p}']+df_medians[f'p{p}'],100)
    fmodel = compute_rv_model_notrends(xmodel,df_medians,p)
    #Let us add the signal of the planet we want to plot
    local_rv += compute_rv_model_notrends(local_time,df_medians,p)
    #other samples
    #compute random samples
    fsamples = []
    for j in range(n_samples):
        fsamples.append(compute_rv_model_notrends(xmodel,df_samples.iloc[j],p))
    xmodel -= df_medians[f't0{p}']
    pmodel = xmodel/df_medians[f'p{p}']

    #Let us phase fold the data
    phase = abs(((local_time-df_medians[f't0{p}'])%df_medians[f'p{p}'])/df_medians[f'p{p}'])

    #Store the RV files to recreate the plots
    frvm = f'{outdir}/{star}{p}_rv_model.dat'
    with open(frvm,'w') as f:
        f.write('#Time  RV \n')
        for j in range(len(pmodel)):
            f.write(f'{pmodel[j]} {fmodel[j]*1e3} \n')

    frvd = f'{outdir}/{star}{p}_rv_data.dat'
    with open(frvd,'w') as f:
        f.write('#Time  RV err jit_err residuals  telescope_label\n')
        for j in range(len(phase)):
            f.write(f'{phase[j]} {local_rv[j]*1e3} {local_errs[j]*1e3}  {local_errs_jit[j]*1e3} {rv_presiduals[j]*1e3}  {local_tel_labs[j]}\n')

    for i,t in enumerate(telescopes_labels):
        ix = i == np.array(tlab)[id0]
        if all(ix == False): continue
        #Models
        ax0.errorbar(phase[ix],local_rv[ix]*1e3,local_errs_jit[ix]*1e3,fmt='.',color=rv_colors[i],alpha=0.5)
        ax0.errorbar(phase[ix],local_rv[ix]*1e3,local_errs[ix]*1e3,mfc='w',fmt=mark[i],color=rv_colors[i],label=t,zorder=5)
        #Residuals
        ax1.errorbar(phase[ix],rv_presiduals[ix]*1e3,local_errs_jit[ix]*1e3,fmt='.',color=rv_colors[i],alpha=0.5)
        ax1.errorbar(phase[ix],rv_presiduals[ix]*1e3,local_errs[ix]*1e3,mfc='w',fmt=mark[i],color=rv_colors[i],zorder=5)

    ax0.plot(pmodel,fmodel*1e3,'k',zorder=5)
    for i in range(n_samples):
        ax0.plot(pmodel,fsamples[i]*1e3,'k',alpha=0.05,lw=0.5,zorder=4)

    ax0.set_xlim(0,1)
    if is_rv_legend: ax0.legend()
    ax0.set_xticks([])
    ax0.set_ylabel(r'RV ($\mathrm{m\,s^{-1}}$)')

    ax1.set_xlim(0,1)
    ax1.set_xlabel('Phase')
    ax1.set_ylabel(r'Residuals ($\mathrm{m\,s^{-1}}$)')

    fname = f'{outdir}/{star}{p}_rv.pdf'
    print(f'Creating  {fname}')
    plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=300)
    plt.savefig(fname[:-3]+'png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

