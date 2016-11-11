#Let us do the plots here

from matplotlib import gridspec
from matplotlib.colors import LogNorm

#The size of our plots follow an aurea rectangle 
fsx = figure_size_x
fsy = figure_size_y
fos = font_size_label

vari = params[0]
chi2 = params[2]

def plot_chains():
  plt.xlabel('iteration')
  plt.ylabel('$\chi^2$')
  if (nplanets == 1):
   #plt.plot(vari,chi2,'b.')
   plt.hist2d(vari,chi2,bins=100,norm=LogNorm())
  else:
   #plt.plot(vari[0],chi2[0],'b.')
   plt.hist2d(vari[0],chi2[0],bins=100,norm=LogNorm())
  fname = outdir+'/'+star+plabels[0]+'_chains.pdf'
  print 'Creating ', fname
  plt.savefig(fname,bbox_inches='tight')
  plt.close()

  #=========================#
  #       Transit plot      #
  #=========================#

#Ntransit is the number of the transit that we want to plot
def fancy_tr_plot(xtime,yflux,errors,pars,ldc,flag,fname):

  t0_val = pars[0]
  P_val  = pars[1]
  e_val = pars[2]
  w_val = pars[3]
  i_val = pars[4]
  a_val = pars[5]
  pz_val = pars[6]

  q1_val = ldc[0]
  q2_val = ldc[1]

  print 'Creating ', fname
  #Let us create the model
  xmodel_res = list(xtime)
  xmodel = np.arange(min(xtime), max(xtime), (max(xtime)-min(xtime))/1.e3 )

  xd_ub = np.ndarray(shape=(n_cad,len(xmodel)))
  xd_ub_res = np.ndarray(shape=(n_cad,len(xmodel)))
  zd_ub = [None]*n_cad
  zd_ub_res = [None]*n_cad
  fd_ub = [None]*n_cad
  fd_ub_res = [None]*n_cad
  #Use the long cadence data
  for m in range(0,n_cad):
    #This vector has the model fit
    for n in range(0,len(xmodel)):
      xd_ub[m][n] = xmodel[n] + t_cad * ( (m+1) - 0.5 * (n_cad + 1 ))/n_cad
    #This vector has the residuals
    for n in range(0,len(xmodel_res)):
      xd_ub_res[m][n] = xmodel_res[n] + t_cad * ( (m+1) - 0.5 * (n_cad + 1))/n_cad

  #Calculate the transit curve for all the data
  for m in range(0,n_cad):
    zd_ub[m] = pti.find_z(xd_ub[m][:],[t0_val,P_val,e_val,w_val,i_val,a_val],flag)
    fd_ub[m], dummm = pti.occultquad(zd_ub[m],q1_val,q2_val,pz_val)
    zd_ub_res[m] = pti.find_z(xd_ub_res[m][:],[t0_val,P_val,e_val,w_val,i_val,a_val],flag)
    fd_ub_res[m], dummm = pti.occultquad(zd_ub_res[m],q1_val,q2_val,pz_val)

  #Bin the data
  fd_reb = [0.0]*len(xmodel)
  fd_reb_res = [0.0]*len(xmodel_res)
  for m in range(0,len(xmodel)):
    for n in range(0,n_cad):
      fd_reb[m] = fd_reb[m] + fd_ub[n][m]/n_cad
  for m in range(0,len(xmodel_res)):
    for n in range(0,n_cad):
      fd_reb_res[m] = fd_reb_res[m] + fd_ub_res[n][m]/n_cad

  #Residuals
  res_res = yflux - fd_reb_res

  #Do the plot
  tfc = 24. # time factor conversion to hours
  local_T0 = t0_val
  plt.figure(1,figsize=(fsx,fsy))
  #Plot the transit light curve
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3.0, 1.])
  gs.update(hspace=0.00)
  ax1 = plt.subplot(gs[0])
  x_lim = (min(xtime)-local_T0)*tfc
  plt.xlim(x_lim,-x_lim)
  min_val_model = max(fd_reb) -  min(fd_reb)
  plt.errorbar((xtime-local_T0)*tfc,yflux,errors,fmt='r.',alpha=1.0)
  plt.plot((xmodel-local_T0)*tfc,fd_reb,'k',linewidth=1.0)
  plt.ylabel('Relative flux',fontsize=fos)
  plt.xticks( np.arange(int(x_lim),int(-x_lim)+1,1))
  plt.minorticks_on()
  plt.ticklabel_format(useOffset=False, axis='y')
  plt.tick_params( axis='x',which='both',labelbottom='off')
  #Plot the residuals
  dplot = plt.subplot(gs[1])
  plt.plot([x_lim,-x_lim],[0.0,0.0],'k--',linewidth=1.0)
  plt.errorbar((xmodel_res-local_T0)*tfc,res_res,errors,fmt='r.',alpha=1.0)
  yylims = dplot.get_ylim()
  plt.yticks(np.arange(yylims[0],yylims[1],(yylims[1]-yylims[0])/4.))
  plt.xticks( np.arange(int(x_lim),int(-x_lim)+1,1))
  plt.xlim(x_lim,-x_lim)
  #Plot the residuals
  plt.ylabel('Residuals',fontsize=fos*0.75)
  plt.xlabel("T - T0 (hours)",fontsize=fos)
  plt.minorticks_on()
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight')
  plt.close()


def plot_transit_nice():
#Move all the points to T0
  base = 3
  for o in range(0,nplanets):

    pars_vec = params[base:base+7]
    pars = [None]*7
    for m in range(0,7):
      pars[m] = np.median(pars_vec[m])

    ldc = [ np.mean(params[3+8*nplanets]), np.median(params[3+9*nplanets]) ]

    xt_dummy = list(xt)
    for i in range(0,ntr):
      P_val = pars[1]
      n = xt_dummy[i][len(xt_dummy[i])-1] - xt_dummy[0][0]
      n = int(n/P_val)
      xt_dummy[i] = xt_dummy[i] - P_val * n

    #Redefine megax with the new xt values
    xtime = np.concatenate(xt_dummy)
    flag = [False, False, is_b_factor, False]
    fname = outdir+'/'+star+plabels[o]+'_tr.pdf'
    fancy_tr_plot(xtime, megay, megae,pars,ldc,flag, fname)
    base = base + 8


def plot_all_transits():
  flag = [False, False, False, False]
  xt_dummy = list(xt)
  for i in range(0,ntr):
    xt_dummy[i] = xt_dummy[i] - P_val * i
  for i in range(0,ntr):
    fname = outdir+'/'+star+plabels[0]+'_transit'+str(i)+'.pdf'
    fancy_tr_plot(np.array(xt_dummy[i]),np.array(yt[i]),np.array(et[i]),flag,fname)

def plot_rv_fancy(p_rv,rvy,p_all,rv_dum,errs_all,res,telescopes_labels,fname):
  print 'Creating ', fname
  plt.figure(3,figsize=(fsx,fsy))
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3., 1.])
  gs.update(hspace=0.00)
  ax0 = plt.subplot(gs[0])
  plt.minorticks_on()
  #plt.subplot(311)
  ax0 = plt.xlabel("")
  ax0 = plt.ylabel("RV (m/s)",fontsize=fos)
  ax0 = plt.plot([0.,1.],[0.,0.],'k--')
  ax0 = plt.plot(p_rv,rvy,'k',linewidth=1.0)
  for j in range(0,nt):
    ax0 = plt.errorbar(p_all[j],rv_dum[j],errs_all[j],\
    label=telescopes_labels[j],\
    fmt=mark[j],\
    alpha=1.0 ,\
    markersize=4)
  plt.legend(loc=0, ncol=1,scatterpoints=1,numpoints=1,frameon=False,fontsize='small')
  plt.xticks(np.arange(0.,1.01,0.1))
  plt.tick_params( axis='x',which='both',labelbottom='off')
  #plt.subplot(312)
  ax1 = plt.subplot(gs[1])
  plt.xlabel("Orbital phase",fontsize=fos)
  plt.tick_params( axis='x',which='minor',bottom='on',left='on',right='on',top='on')
  plt.xticks(np.arange(0.,1.01,0.1))
  plt.ylabel('Residuals (m/s)',fontsize=fos*0.75)
  plt.plot([0.,1.],[0.,0.],'k--',linewidth=1.0)
  for j in range(0,nt):
    plt.errorbar(p_all[j],res[j],errs_all[j],\
    label=telescopes_labels[j],fmt=mark[j],alpha=1.0,markersize=4)
  yylims = ax1.get_ylim()
  plt.yticks(np.arange(yylims[0],yylims[1],(yylims[1]-yylims[0])/4.))
  plt.minorticks_on()
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight')
  plt.close()

  #=========================#
  #        RV plot          #
  #=========================#

#===========================================================
#                   One planet plots
#===========================================================

v_vec_val = [None]*nt
v_val = [None]*nt
#3 + npars + ldc
v_vec_val[:] = params[3+10*nplanets:3+10*nplanets+nt]
for o in range(0,nt):
  v_val[o] = np.median(v_vec_val[o])

base = 3
t0_val = np.ndarray(nplanets)
P_val  = np.ndarray(nplanets)
e_val  = np.ndarray(nplanets)
w_val  = np.ndarray(nplanets)
k_val  = np.ndarray(nplanets)
alpha_val = 0.0
beta_val = 0.0

for o in range(0,nplanets):
  t0_val[o] = np.median(params[base + 0])
  P_val[o]  = np.median(params[base + 1])
  e_val[o]  = np.median(params[base + 2])
  w_val[o]  = np.median(params[base + 3])
  k_val[o]  = np.median(params[base + 7])
  base = base + 8

if ( nplanets > 0 ):

  #Plot without fold the data
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
    for j in range(0,nt):
      #This is the model of the actual planet
      #the actual value, minus the systemic velocity
      rv_dum[j] = rv_datas[j] - v_val[j]*cfactor

    plt.figure(1,figsize=(fsx,fsy))
    plt.plot(rvx,rvy,'k')
    plt.minorticks_on()
    plt.xlabel("JD (days)",fontsize=fos)
    plt.ylabel('RV (m/s)',fontsize=fos)
    plt.xlim(xmin,xmax)
    for j in range(0,nt):
      plt.errorbar(time_all[j],rv_dum[j],errs_datas[j],\
      label=telescopes_labels[j],fmt=mark[j],alpha=1.0)
    plt.legend(loc=0, ncol=1,scatterpoints=1,numpoints=1,frameon=False,fontsize='small')
    fname = outdir+'/'+star+plabels[0]+'_rv_all.pdf'
    print 'Creating ', fname
    plt.savefig(fname,format='pdf',bbox_inches='tight')
    plt.close()

#===========================================================

#Plot RV for one planet
  def plot_rv_one():
    cfactor = np.float(1.e3)
    rv_dat2 = [None]*len(rv_all)
    errs_dat2 = [None]*len(errs_all)
    for i in range(0,nt):
      rv_dat2[i] = list(rv_all[i])
      errs_dat2[i] = list(errs_all[i])
    for i in range(0,nt):
      for j in range(0,len(rv_dat2[i])):
        rv_dat2[i][j] = cfactor*rv_dat2[i][j]
        errs_dat2[i][j] = cfactor*errs_dat2[i][j]

    #Let us save all the RV data in rv_dum
    n = 5000
    xmin = t0_val
    xmax = t0_val + P_val[0]
    dn = (xmax - xmin) /  n
    rvx = np.empty([n])
    rvx[0] = xmin
    for j in range(1,n):
      rvx[j] = rvx[j-1] + dn

    #Model curve
    rvy = pti.rv_curve_mp(rvx,0.0,t0_val[0],\
                          k_val[0]*cfactor,P_val[0],e_val[0],w_val[0],0.0,0.0)
    res = [None]*nt
    rv_dum = [None]*nt
    for j in range(0,nt):
        #This is the model of the actual planet
        res[j] = pti.rv_curve_mp(time_all[j],0.0,t0_val[0],k_val[0]*cfactor,\
                                 P_val[0],e_val[0],w_val[0],0.0,0.0)
        alpha_time = [None]*len(time_all[j])
        beta_time = [None]*len(time_all[j])
        for m in range(0,len(time_all[j])):
            alpha_time[m] = (time_all[j][m]-t0_val)**1 * alpha_val[0] * cfactor
            beta_time[m]  = (time_all[j][m]-t0_val)**2 * beta_val[0]  * cfactor

        rv_dum[j] = rv_dat2[j] - v_val[j]*cfactor - alpha_time - beta_time
        res[j] = rv_dum[j] - res[j]

        p_rv = scale_period(rvx,t0_val,P_val)
        p_all = [None]*nt
        for j in range(0,nt):
            p_all[j] = scale_period(time_all[j],t0_val,P_val)

    fname = outdir+'/'+star+plabels[0]+'_rv.pdf'
    plot_rv_fancy(p_rv,rvy,p_all,rv_dum,errs_dat2,res,telescopes_labels,fname)


#===========================================================
#                   Multi-planet plots
#===========================================================

#else:
  #Plot RV for multiplanet
  def plot_rv_mp():
    cfactor = np.float(1.e3)
    rvy = [None]*nplanets
    p_rv = [None]*nplanets
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
        rv_dum.append(rv_all[j])
      #Create the RV fitted model for the planet i
      n = 5000
      xmin = t0_val[i]
      xmax = t0_val[i] + P_val[i]
      dn = (xmax - xmin) /  n
      rvx = np.empty([n])
      rvx[0] = xmin
      for j in range(1,n):
        rvx[j] = rvx[j-1] + dn
      rvy[i] = pti.rv_curve_mp(rvx,0.0,t0_val[i],\
      k_dum[i],P_val[i],e_val[i],w_val[i],alpha_val*cfactor \
      , beta_val*cfactor)

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
      for j in range(0,nt):
        #This is the model of the actual planet
        res[j] = pti.rv_curve_mp(time_all[j],0.0,t0_val[i],k_dum[i],\
        P_val[i],e_val[i],w_val[i], cfactor*alpha_val, cfactor*beta_val)

        #This variable has all the others planets
        if ( nplanets > 1 ):
          drvy[j] = pti.rv_curve_mp(time_all[j],0.0,dt0_val,dk_dum \
          ,dP_val,de_val,dw_val, alpha_val*cfactor \
          ,beta_val*cfactor)
        else:
          drvy[j] = 0.0

        #the actual value, minus the systemic velocity, minus the other planets
        rv_dum[j] = rv_dum[j] - v_val[j] - drvy[j]
        res[j] = rv_dum[j] - res[j]

      p_rv[i] = scale_period(rvx,t0_val[i],P_val[i])
      p_all = [None]*nt
      for j in range(0,nt):
        p_all[j] = scale_period(time_all[j],t0_val[i],P_val[i])

      fname = outdir+'/'+star+plabels[i]+'_rv.pdf'
      plot_rv_fancy(p_rv[i],rvy[i],p_all,rv_dum,errs_all,res,telescopes_labels,fname)

#===========================================================
#                   Histogram plots
#===========================================================

def create_plot_histogram(params,plabs,cbars='red',nb=50):
  n = len(params)
  plt.figure(1,figsize=(15,4*(n)/2))
  gs = gridspec.GridSpec(nrows=(n+1)/2,ncols=2)
  for i in range(0,n):
    plt.subplot(gs[i])
    vpar, lpar, rpar = find_vals_perc(params[i],1.0)
    #best_val = params[i][minchi2_index]
    #plt.axvline(x=best_val,c='yellow')
    plt.axvline(x=vpar,c=cbars)
    plt.axvline(x=vpar-lpar,c=cbars,ls='--')
    plt.axvline(x=vpar+rpar,c=cbars,ls='--')
    plt.xlabel(plabs[i])
    plt.hist(params[i],normed=True,bins=nb)

  fname = outdir+'/'+star+plabels[0]+'_histogram.pdf'
  print 'Creating ', fname
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.close()


def plot_histogram_2():
    plabs = ['dummy']*len(params[3:])
    create_plot_histogram(params[3:], plabs, cbars='red', nb=50)

def plot_histogram(rf=1):

  if ( fit_tr and fit_rv ):
    dparams = [t0o[0::rf],Po[0::rf],eo[0::rf],wo[0::rf],bo[0::rf],ao[0::rf],q1o[0::rf],q2o[0::rf],pzo[0::rf],ko[0::rf]]
    dplabs = ['$T0$','$P$','$e$','$\omega$','$b$','$a/R_*$','$q_1$','$q_2$','$R_p/R_*$','$k$']

    vlabs = [None]*nt
    dumvo = [None]*nt
    for i in range(0,nt):
      vlabs[i] = 'rv0 ' + telescopes_labels[i]
      dumvo[i] = vo[i][0::rf]

    params = np.concatenate([dparams,dumvo])
    labs = np.concatenate([dplabs,vlabs])

    create_plot_histogram(params,labs)


  if ( fit_tr and not fit_rv ):
    params = [t0o[0::rf],Po[0::rf],eo[0::rf],wo[0::rf],bo[0::rf],ao[0::rf],q1o[0::rf],q2o[0::rf],pzo[0::rf]]
    labs = ['$T0$','$P$','$e$','$\omega$','$b$','$a/R_*$','$q_1$','$q_2$','$R_p/R_*$']

    create_plot_histogram(params,labs)

  if ( not fit_tr and fit_rv ):
    if (nplanets == 1 ):
      dparams = [t0o[0::rf],Po[0::rf],eo[0::rf],wo[0::rf],ko[0::rf],alphao[0::rf],betao[0::rf]]
      dplabs = ['$T0$','$P$','$e$','$\omega$','$k$','alpha','beta']
    else:
      dparams = [None]*(7*nplanets)
      dplabs = [None]*(7*nplanets)
      for i in range(0,nplanets):
        dparams[0+7*i] = t0o[i][0::rf]
        dparams[1+7*i] = Po[i][0::rf]
        dparams[2+7*i] = eo[i][0::rf]
        dparams[3+7*i] = wo[i][0::rf]
        dparams[4+7*i] = ko[i][0::rf]
        dparams[5+7*i] = alphao[i][0::rf]
        dparams[6+7*i] = betao[i][0::rf]
        dplabs[0+7*i] = 'T0'+plabels[i]
        dplabs[1+7*i] = 'P'+plabels[i]
        dplabs[2+7*i] = 'e'+plabels[i]
        dplabs[3+7*i] = '$\omega$'+plabels[i]
        dplabs[4+7*i] = 'k'+plabels[i]
        dplabs[5+7*i] = '$\\alpha$'+plabels[i]
        dplabs[6+7*i] = '$\\beta$'+plabels[i]

    vlabs = [None]*nt
    dvo = [None]*nt
    for i in range(0,nt):
      vlabs[i] = 'rv0 ' + telescopes_labels[i]
      dvo[i] = vo[i][0::rf]

    params = np.concatenate([dparams,dvo])
    labs = np.concatenate([dplabs,vlabs])

    create_plot_histogram(params,labs)

#===========================================================
#                   Correlation plots
#===========================================================

def create_plot_correlation(params,plabs,col='red',mark='.'):
  n = len(params)
  plt.figure(1,figsize=(2*n,2*n))
  gs = gridspec.GridSpec(nrows=n,ncols=n)
  for i in range(0,n):
    for j in range(0,i):
      plt.subplot(gs[i*n+j])
      if ( j == 0 ):
        plt.ylabel(plabs[i])
      elif ( j == i - 1 ):
        #plt.colorbar()
        plt.tick_params( axis='y',which='both',labelleft='off')
      else:
        plt.tick_params( axis='y',which='both',labelleft='off')
      if ( i == n-1):
        plt.xlabel(plabs[j])
      else:
        plt.tick_params( axis='x',which='both',labelbottom='off')

      plt.hist2d(params[j],params[i],bins=100,norm=LogNorm())

  fname = outdir+'/'+star+plabels[0]+'_correlations.pdf'
  print 'Creating ', fname
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.close()

def plot_correlations_2():
  labs = ['dummy']*len(params[3:])
  create_plot_correlation(params[3:],labs,col='blue')

def plot_correlations(rf=1):

  if ( fit_tr and fit_rv ):
    dparams = [t0o[0::rf],Po[0::rf],eo[0::rf],wo[0::rf],bo[0::rf],ao[0::rf],q1o[0::rf],q2o[0::rf],pzo[0::rf],ko[0::rf]]
    dplabs = ['$T0$','$P$','$e$','$\omega$','$b$','$a/R_*$','$q_1$','$q_2$','$R_p/R_*$','$k$']

    vlabs = [None]*nt
    dvo = [None]*nt
    for i in range(0,nt):
      vlabs[i] = 'rv0 ' + telescopes_labels[i]
      dvo[i] = vo[i][0::rf]

    params = np.concatenate([dparams,dvo])
    labs = np.concatenate([dplabs,vlabs])

    create_plot_correlation(params,labs,col='blue')

  if ( fit_tr and not fit_rv ):

    params = [t0o[1::rf],Po[1::rf],eo[1::rf],wo[1::rf],bo[1::rf],ao[1::rf],q1o[1::rf],q2o[1::rf],pzo[1::rf]]
    labs = ['$T0$','$P$','$e$','$\omega$','$b$','$a/R_*$','$q_1$','$q_2$','$R_p/R_*$']

    create_plot_correlation(params,labs,col='blue')

  #Now it works only for RV fit
  if ( not fit_tr and fit_rv ):
    if ( nplanets == 1 ):
      dparams = [t0o[1::rf],Po[1::rf],eo[1::rf],wo[1::rf],ko[1::rf],alphao[1::rf],betao[1::rf]]
      dplabs = ['$T0$','$P$','$e$','$\omega$','$k$','alpha','beta']
    else:
      dparams = [None]*(7*nplanets)
      dplabs = [None]*(7*nplanets)
      for i in range(0,nplanets):
        dparams[0+7*i] = t0o[i][1::rf]
        dparams[1+7*i] = Po[i][1::rf]
        dparams[2+7*i] = eo[i][1::rf]
        dparams[3+7*i] = wo[i][1::rf]
        dparams[4+7*i] = ko[i][1::rf]
        dparams[5+7*i] = alphao[i][1::rf]
        dparams[6+7*i] = betao[i][1::rf]
        dplabs[0+7*i] = 'T0'+ plabels[i]
        dplabs[1+7*i] = 'P'+ plabels[i]
        dplabs[2+7*i] = 'e'+ plabels[i]
        dplabs[3+7*i] = '$\omega$'+ plabels[i]
        dplabs[4+7*i] = '$k$'+ plabels[i]
        dplabs[5+7*i] = '$\\alpha$'+ plabels[i]
        dplabs[6+7*i] = '$\\beta$'+ plabels[i]

    vlabs = [None]*nt
    dvo = [None]*nt
    for i in range(0,nt):
      vlabs[i] = 'rv0 ' + telescopes_labels[i]
      dvo[i] = vo[i][1::rf]

    params = np.concatenate([dparams,dvo])
    labs = np.concatenate([dplabs,vlabs])

    create_plot_correlation(params,labs,col='blue')

