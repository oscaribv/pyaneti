#===========================================================
#              plot rv fancy function
#===========================================================

def plot_rv_fancy(p_rv,rvy,p_all,rv_dum,errs_all,res,telescopes_labels,fname,is_special=False):
  print 'Creating ', fname
  if ( not is_special ):
      rv_model = rvy
  else:
      rv_model = rvy[0]
  #
  plt.figure(3,figsize=(fsx,fsy))
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3., 1.])
  gs.update(hspace=0.00)
  ax0 = plt.subplot(gs[0])
  plt.tick_params(labelsize=fos,direction='in')
  plt.minorticks_on()
  #plt.subplot(311)
  plt.xlabel("")
  plt.ylabel("RV (m/s)",fontsize=fos)
  #plt.plot([0.,1.],[0.,0.],'k--')
#  if ( is_special ):
#    [plt.fill_between(p_rv,*rvy[i:i+2,:],alpha=0.3,facecolor='k') for i in range(1,6,2)]
#    plt.fill_between(p_rv,*rvy[5:7,:],alpha=0.5,facecolor='k')
  for j in range(0,nt):
    #
    plt.errorbar(p_all[j],rv_dum[j],new_errs_all[j],\
    fmt=mark[j],\
    alpha=1.0 ,color='#C0C0C0',\
    markersize=rv_markersize,fillstyle=rv_fillstyle)
    #
    plt.errorbar(p_all[j],rv_dum[j],errs_all[j],\
    label=telescopes_labels[j],\
    fmt=mark[j],\
    alpha=1.0 ,color=rv_colors[j],\
    markersize=rv_markersize,fillstyle=rv_fillstyle)
  #
  plt.plot(p_rv,rv_model,'k',linewidth=2.0,alpha=0.9,zorder=3)
  #
  if ( is_rv_legend ): plt.legend(loc=2, ncol=1,scatterpoints=1,numpoints=1,frameon=True,fontsize=fos*0.7)
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
  #plt.yticks(np.arange(-miy,miy,2.*miy/8.))
#  plt.yticks(np.arange(yylims[0],yylims[1],(yylims[1]-yylims[0])/7.))
  #NEW SUBPLOT
  ax1 = plt.subplot(gs[1])
  plt.tick_params(labelsize=fos,direction='in')
  plt.xlabel("Orbital phase",fontsize=fos)
  plt.tick_params( axis='x',which='minor',direction='in',bottom=True,left=True,right=True,top=True)
  plt.tick_params( axis='y',which='both',direction='in')
  plt.xticks(np.arange(0.,1.01,0.1))
  plt.ylabel('Residuals (m/s)',fontsize=fos*0.75)
  plt.plot([0.,1.],[0.,0.],'k--',linewidth=1.0)
  for j in range(0,nt):
    #
    plt.errorbar(p_all[j],res[j],new_errs_all[j],\
    fmt=mark[j],\
    alpha=1.0 ,color='#C0C0C0',\
    markersize=rv_markersize,fillstyle=rv_fillstyle)
    #
    plt.errorbar(p_all[j],res[j],errs_all[j],\
    label=telescopes_labels[j],fmt=mark[j],color=rv_colors[j], \
    alpha=1.0,markersize=rv_markersize,fillstyle=rv_fillstyle)
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
    for i in range(0,nt):
      for j in range(0,len(rv_all[i])):
        rv_datas[i][j] = cfactor*rv_datas[i][j]
        errs_datas[i][j] = cfactor*errs_datas[i][j]

    #Let us save all the RV data in rv_dum
    n = 5000
    xmin = min(np.concatenate(time_all))
    xmax = max(np.concatenate(time_all))
    total_tt = xmax - xmin
    agregar = total_tt*0.1
    xmax = xmax + agregar
    xmin = xmin - agregar
    dn = (xmax - xmin) /  n
    rvx = np.empty([n])
    rvx[0] = xmin
    for j in range(1,n):
      rvx[j] = rvx[j-1] + dn

    #Model curve
    rvy = pti.rv_curve_mp(rvx,0.0,t0_val,\
    k_val*cfactor,P_val,e_val,w_val,alpha_val*cfactor,beta_val*cfactor)
    rv_dum = [None]*nt
    res_dum_all = [None]*nt
    for j in range(0,nt):
      res_dum_all[j] = pti.rv_curve_mp(time_all[j],0.0,t0_val,\
      k_val*cfactor,P_val,e_val,w_val,alpha_val*cfactor,beta_val*cfactor)
      #This is the model of the actual planet
      #the actual value, minus the systemic velocity
      rv_dum[j] = rv_datas[j] - v_val[j]*cfactor
      res_dum_all[j] = rv_dum[j] - res_dum_all[j]

    plt.figure(1,figsize=(2*fsx,fsy))
    plt.plot(rvx,rvy,'k')

    plt.minorticks_on()
    plt.xlabel(rv_xlabel,fontsize=fos)
    plt.ylabel('RV (m/s)',fontsize=fos)
    plt.xlim(xmin,xmax)
    plt.tick_params(labelsize=fos,direction='in')
    for j in range(0,nt):
      plt.errorbar(time_all[j],rv_dum[j],new_errs_all[j],color='#C0C0C0',\
      fmt=mark[j],alpha=1.0,markersize=rv_markersize)
      #
      plt.errorbar(time_all[j],rv_dum[j],errs_datas[j],color=rv_colors[j],\
      label=telescopes_labels[j],fmt=mark[j],alpha=1.0,markersize=rv_markersize)
    if ( is_rv_legend ): plt.legend(loc=2, ncol=1,scatterpoints=1,numpoints=1,frameon=True,fontsize=fos*0.8)
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
    k_dum = [None]*nplanets
    for i in range(0,nplanets):
      k_dum[i] = cfactor*k_val[i]
    for i in range(0,nt):
      v_val[i] = cfactor*v_val[i]
      for j in range(0,len(rv_all[i])):
        rv_all[i][j] = cfactor*rv_all[i][j]
        errs_all[i][j] = cfactor*errs_all[i][j]

    for i in range(0,nplanets):
      rv_dum = []
      for j in range(0,nt):
        rv_dum.append(list(rv_all[j]))
      #Create the RV fitted model for the planet i
      rvx = np.arange(t0_val[i],t0_val[i]+P_val[i]*0.999,P_val[i]/4999.)
      rvy = pti.rv_curve_mp(rvx,0.0,t0_val[i],\
      k_dum[i],P_val[i],e_val[i],w_val[i],0.0 ,0.0)

      #If we want to plot the percentiles
      if ( is_special_plot_rv ):

        #len of the chain vector
        len_chain = len(params[0])

        nc = len_chain/50

        rv_vector = [None]*nc
        rv_vector_res = [None]*nc
        for l in range(0,nc):
          #generate a random chain number
          mi_chain = int(np.random.uniform(0,len_chain))
          #Call the parameters
          lt0, lp, le, lw, lk = pars_rv_chain(params,mi_chain)
          rv_vector[l] = pti.rv_curve_mp(rvx,0.0,lt0[i],\
          lk[i]*cfactor,lp[i],le[i],lw[i],0.0 ,0.0)
          #rv_vector_res[l] = pti.rv_curve_mp(time_all,0.0,lt0,\
          #lk*cfactor,lp,le,lw,0.0 ,0.0)

        rv_vector = np.array(rv_vector)
        #rv_vector_res = np.array(rv_vector_res)
        rvy = np.percentile(rv_vector, [50.,0.15,99.85, 2.5,97.5,16.,84.], 0)
        #rvy_res = np.percentile(res_vector_res, [50,0.15,99.85, 2.5,97.5, 16,84], 0)
        #fd_ub_res = res_pc_res[0]

      dt0_val = []
      dk_dum = []
      dP_val = []
      de_val = []
      dw_val = []

      j = 0
      while ( j < nplanets ):
        if ( j != i and nplanets>1 ):
          dt0_val.append(t0_val[j])
          dk_dum.append(k_dum[j])
          dP_val.append(P_val[j])
          de_val.append(e_val[j])
          dw_val.append(w_val[j])
        j = j + 1

      res = [None]*nt
      drvy = [None]*nt
      for j in range(0,nt): #Remove the signal for each telescope
        #This is the model of the actual planet
        res[j] = pti.rv_curve_mp(time_all[j],0.0,t0_val[i],k_dum[i],\
        P_val[i],e_val[i],w_val[i], 0.0, 0.0)

        #This variable has all the others planets
        if ( nplanets > 1 ):
          drvy[j] = pti.rv_curve_mp(time_all[j],0.0,dt0_val,dk_dum \
          ,dP_val,de_val,dw_val, 0.0 \
          ,0.0)
        else:
          drvy[j] = 0.0

        alpha_time = [0.0]*len(time_all[j])
        beta_time = [0.0]*len(time_all[j])
        if ( is_linear_trend or is_quadratic_trend ):
          for m in range(0,len(time_all[j])):
              alpha_time[m] = (time_all[j][m]-t0_val[0])**1 * alpha_val * cfactor
              beta_time[m]  = (time_all[j][m]-t0_val[0])**2 * beta_val  * cfactor

        #the actual value, minus the systemic velocity, minus the other planets
        for o in range(len(time_all[j])):
          if ( nplanets > 1 ):
            rv_dum[j][o] = rv_dum[j][o] - v_val[j] - drvy[j][o] - alpha_time[o] - beta_time[o]
          else:
            rv_dum[j][o] = rv_dum[j][o] - v_val[j] - drvy[j] - alpha_time[o] - beta_time[o]
          res[j][o] = rv_dum[j][o] - res[j][o]

      p_rv = scale_period(rvx,t0_val[i],P_val[i])
      p_all = [None]*nt
      for j in range(0,nt):
        p_all[j] = scale_period(time_all[j],t0_val[i],P_val[i])

      fname = outdir+'/'+star+plabels[i]+'_rv.pdf'
      plot_rv_fancy(p_rv,rvy,p_all,rv_dum,errs_all,res,telescopes_labels,fname,is_special_plot_rv)

