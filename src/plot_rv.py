#===========================================================
#              plot rv fancy function
#===========================================================
def plot_rv_fancy(rv_vector,fname):

  print 'Creating ', fname
  tmodel  = rv_vector[0] #time vector for the model to plot
  rvmodel = rv_vector[1] #rv vector for the model to plot
  tdata   = rv_vector[2]
  rvdata  = rv_vector[3]
  edata   = rv_vector[4]
  ejdata  = rv_vector[5]
  res     = rv_vector[6]
  tel_lab = rv_vector[7]
  #
  plt.figure(3,figsize=(fsx,fsy))
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3., 1.])
  gs.update(hspace=0.00)
  #
  #RV curve plot
  #
  ax0 = plt.subplot(gs[0])
  plt.tick_params(labelsize=fos,direction='in')
  plt.minorticks_on()
  plt.xlabel("")
  plt.ylabel("RV (m/s)",fontsize=fos)
  #PLOT DATA
  for j in range(0,len(tdata)):
    #
    plt.errorbar(tdata[j],rvdata[j],ejdata[j],fmt=mark[tel_lab[j]],\
    alpha=1.0 ,color='#C0C0C0',markersize=rv_markersize,fillstyle=rv_fillstyle)
    #
    plt.errorbar(tdata[j],rvdata[j],edata[j],\
    fmt=mark[tel_lab[j]],alpha=1.0,color=rv_colors[tel_lab[j]], \
    markersize=rv_markersize,fillstyle=rv_fillstyle)
  for j in range(0,len(telescopes_labels)):
    plt.errorbar(-1,0,1e-5,fmt=mark[j],label=telescopes_labels[j],\
    alpha=1.0,color=rv_colors[j], \
    markersize=rv_markersize)
  #PLOT MODEL
  plt.plot(tmodel,rvmodel,'k',linewidth=2.0,alpha=0.9,zorder=3)
  #
  if ( is_rv_legend ): plt.legend(loc=2, ncol=1,scatterpoints=1,numpoints=1,frameon=True,fontsize=fos*0.7)
  #
  plt.xticks(np.arange(0.,1.01,0.1))
  plt.tick_params( axis='x',which='both',direction='in',labelbottom=False)
  plt.tick_params( axis='y',which='both',direction='in')
  yylims = ax0.get_ylim()
  miy = int(max(abs(yylims[0]),abs(yylims[1])))
  if ( miy == 0): #To avoid errors for really small RV variations
    miy = max(abs(yylims[0]),abs(yylims[1]))
  plt.ylim(-miy,miy)
  if ( select_y_rv ):
    plt.ylim(rv_lim_min,rv_lim_max)
  plt.xlim(0.,1.)
  #NEW SUBPLOT: RESIDUALS
  ax1 = plt.subplot(gs[1])
  plt.tick_params(labelsize=fos,direction='in')
  plt.xlabel("Orbital phase",fontsize=fos)
  plt.tick_params( axis='x',which='minor',direction='in',bottom=True,left=True,right=True,top=True)
  plt.tick_params( axis='y',which='both',direction='in')
  plt.xticks(np.arange(0.,1.01,0.1))
  plt.ylabel('Residuals (m/s)',fontsize=fos*0.75)
  #PLOT DATA
  for j in range(0,len(tdata)):
    #
    plt.errorbar(tdata[j],res[j],edata[j],fmt=mark[tel_lab[j]],\
    alpha=1.0 ,color='#C0C0C0',markersize=rv_markersize,fillstyle=rv_fillstyle)
    #
    plt.errorbar(tdata[j],res[j],edata[j],label=telescopes_labels[tel_lab[j]],\
    fmt=mark[tel_lab[j]],alpha=1.0,color=rv_colors[tel_lab[j]], \
    markersize=rv_markersize,fillstyle=rv_fillstyle)
  plt.plot([0.,1.],[0.,0.],'k--',linewidth=1.0)
  #
  yylims = ax1.get_ylim()
  miy = int(max(abs(yylims[0]),abs(yylims[1])))
  if ( miy == 0): #To avoid errors for really small RV variations
    miy = max(abs(yylims[0]),abs(yylims[1]))
  plt.yticks(np.arange(-miy,miy,2.*miy/4.))
  plt.ylim(-miy,miy)
  plt.minorticks_on()
  plt.xlim(0.,1.)
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight',dpi=300)
  plt.close()


#===========================================================
#                   RV PLOTS
#===========================================================


def pars_rv_chain(params,nchain):

  #v_vec_val = [None]*nt
  #v_val = [None]*nt
  ##3 + npars + ldc
  #v_vec_val[:] = params[4+8*nplanets+2:4+8*nplanets+2+nt]
  #for o in range(0,nt):
    #v_val[o] = best_value(v_vec_val[o],maxloglike,get_value)
    #if ( is_log_rv0 ):
      #v_val[o] = 10.0**(v_val[o])

  alpha_val = 0.0
  beta_val = 0.0
  if ( is_linear_trend != 'f' or is_quadratic_trend != 'f' ):
    alpha_val = params_trends[0][nchain]
    beta_val  = params_trends[1][nchain]

  base = 4
  t0_val = np.ndarray(nplanets)
  P_val  = np.ndarray(nplanets)
  e_val  = np.ndarray(nplanets)
  w_val  = np.ndarray(nplanets)
  k_val  = np.ndarray(nplanets)

  for o in range(0,nplanets):
    t0_val[o] = params[base + 0][nchain]
    P_val[o]  = params[base + 1][nchain]
    e_val[o]  = params[base + 2][nchain]
    w_val[o]  = params[base + 3][nchain]
    k_val[o]  = params[base + 7][nchain]
    if ( is_log_P ):
      P_val[o] = 10.0**(P_val[o])
    if ( is_log_k ):
      k_val[o] = 10.0**(k_val[o])
    if ( is_ew ):
      edum_val = e_val[o]
      e_val[o] = e_val[o]**2 + w_val[o]**2
      w_val[o] = np.arctan2(edum_val,w_val[o])

    base = base + 8

  return t0_val, P_val, e_val, w_val, k_val

if ( nplanets > 0 ):
  #Plot without fold the data
  #===========================================================
  #      Plot the light curve with all the parameters
  #===========================================================
  def plot_rv_all_data():

    cfactor = np.float(1.e3)
    rv_datas = [None]*nt
    errs_datas = [None]*nt
    for i in range(0,nt):
      rv_datas[i] = list(rv_all[i])
      errs_datas[i] = list(errs_all[i])

    #Let us save all the RV data in rv_dum
    n = 5000
    xmin = min(np.concatenate(time_all))
    xmax = max(np.concatenate(time_all))
    total_tt = xmax - xmin
    agregar = total_tt*0.1
    xmax = xmax + agregar
    xmin = xmin - agregar
    rvx = np.arange(xmin,xmax,(xmax-xmin)/n)

    #Model curve
    rvy = pti.rv_curve_mp(rvx,0.0,t0_val,k_val,P_val,e_val,w_val,alpha_val,beta_val)

    rv_dum = [None]*nt
    res_dum_all = [None]*nt
    for j in range(0,nt):
      res_dum_all[j] = pti.rv_curve_mp(time_all[j],0.0,t0_val,\
      k_val,P_val,e_val,w_val,alpha_val,beta_val)
      #This is the model of the actual planet
      #the actual value, minus the systemic velocity
      rv_dum[j] = np.asarray(rv_datas[j] - v_val[j])
      res_dum_all[j] = np.asarray(rv_dum[j] - res_dum_all[j])

    #start the plot
    plt.figure(1,figsize=(2*fsx,fsy))

    #are we plotting a GP together with the RV curve
    if kernel_rv[0:2] != 'No':
        xvec = mega_time
        yvec = np.concatenate(res_dum_all)
        evec = mega_err
        m, C =pti.pred_gp(kernel_rv,pk_rv,xvec,yvec,evec,rvx,jrv,jrvlab)
        plt.plot(rvx,rvy*cfactor,'r',alpha=0.7,label='Planetary signal')
        plt.plot(rvx,m*cfactor,'b',alpha=0.7,label='Gaussian Process')
        plt.plot(rvx,(rvy+m)*cfactor,'k',label='Planetary signal + GP')
        sig = np.sqrt(np.diag(C))
        plt.fill_between(rvx,cfactor*(m+sig),cfactor*(m-sig),color='b',alpha=0.1)
    else:
        plt.plot(rvx,(rvy)*cfactor,'k',label='Planetary signal')

    plt.minorticks_on()
    plt.xlabel(rv_xlabel,fontsize=fos)
    plt.ylabel('RV (m/s)',fontsize=fos)
    plt.xlim(xmin,xmax)
    plt.tick_params(labelsize=fos,direction='in')
    for j in range(0,nt):
      plt.errorbar(time_all[j],rv_dum[j]*cfactor,np.asarray(new_errs_all[j]),color='#C0C0C0',\
      fmt=mark[j],alpha=1.0,markersize=rv_markersize)
      #
      plt.errorbar(time_all[j],rv_dum[j]*cfactor,np.asarray(errs_datas[j])*cfactor,color=rv_colors[j],\
      label=telescopes_labels[j],fmt=mark[j],alpha=1.0,markersize=rv_markersize)
    if ( is_rv_legend ): plt.legend(ncol=1,scatterpoints=1,numpoints=1,frameon=True,fontsize=fos*0.8)
    fname = outdir+'/'+star+'_rv_all.pdf'
    print 'Creating ', fname
    plt.savefig(fname,format='pdf',bbox_inches='tight')
    plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight',dpi=300)
    plt.close()

    #Let us create or detrended file
    out_f = outdir+'/'+star+'_rv_residuals.dat'
    vec_x = np.concatenate(time_all)
    vec_y = np.concatenate(res_dum_all)
    vec_z = np.concatenate(errs_all)
    of = open(out_f,'w')
    #of.write('#This detrended light curve was created with pyaneti/lunas\n')
    for i in range(0,len(vec_x)):
      of.write(' %8.8f   %8.8f  %8.8f \n'%(vec_x[i],vec_y[i],vec_z[i]*cfactor))

    of.close()

    return res_dum_all

#===========================================================
#                RV multi-planet fit
#===========================================================
  def plot_rv_mp():

    cfactor = np.float(1.e3)

    for i in range(0,nplanets):

      #Create the RV fitted model for the planet i
      rvx = np.arange(t0_val[i],t0_val[i]+P_val[i]*0.999,P_val[i]/4999.)
      rvy = pti.rv_curve_mp(rvx,0.0,t0_val[i],k_val[i],P_val[i],e_val[i],w_val[i],0.0,0.0)
      #rvx and rvy are the model timeseries for planet i

      #Now it is time to remove the planets j != i from the data
      rv_pi  = pti.rv_curve_mp(mega_time,0.0,t0_val[i],k_val[i],P_val[i],e_val[i],w_val[i],0.,0.)
      #This variable has all the planets
      rv_pall = pti.rv_curve_mp(mega_time,0.0,t0_val,k_val,P_val,e_val,w_val,0.0,0.0)


      #Let us remove all the signals from the data
      res = np.zeros(shape=len(mega_rv))
      for m in range(0,len(mega_rv)):
        res[m] = mega_rv[m] - v_val[tlab[m]] - rv_pall[m] - alpha_val - beta_val
      #Let us add the  signal of the planet i to the data
      rv_planet_i = res + rv_pi

      #Did we fit for a GP?
      evec = np.asarray(mega_err)
      if kernel_rv[0:2] != 'No':
         xvec = mega_time
         yvec = res
         kernel_val, C = pti.pred_gp(kernel_rv,pk_rv,xvec,yvec,evec,xvec,jrv,jrvlab)
         res = res - kernel_val
         rv_planet_i = rv_planet_i - kernel_val

      rvy = np.asarray(rvy)*cfactor
      res = np.asarray(res)*cfactor
      rv_planet_i = np.asarray(rv_planet_i)*cfactor
      evec = evec*cfactor
      ejvec = np.concatenate(new_errs_all)

      p_rv  = scale_period(rvx,t0_val[i],P_val[i])
      p_all = scale_period(mega_time,t0_val[i],P_val[i])

      fname = outdir+'/'+star+plabels[i]+'_rv.pdf'
      plot_rv_fancy([p_rv,rvy,p_all,rv_planet_i,evec,ejvec,res,tlab],fname)
