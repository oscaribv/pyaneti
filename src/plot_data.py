#Let us do the plots here

from matplotlib import gridspec
from matplotlib.colors import LogNorm

if ( is_seaborn_plot ):
  import seaborn as sns
  sns.set(style='ticks')
  sns.set_color_codes()

fsx = figure_size_x
fsy = figure_size_y
fos = font_size_label

vari = params[0]
chi2 = params[2] + params[3]


#===========================================================
#              plot chains
#===========================================================

def plot_chains():
  plt.xlabel('iteration')
  plt.ylabel('Reduced $\chi^2$')
  plt.hist2d(vari,chi2/(ndata-npars),bins=50)
  fname = outdir+'/'+star+'_chains.pdf'
  print 'Creating ', fname
  plt.savefig(fname,bbox_inches='tight')
  plt.close()

def plot_likelihood():
  plt.xlabel('iteration')
  plt.ylabel('$\ln \mathcal{L}$')
  plt.hist2d(vari,par_likelihood,bins=50)
  fname = outdir+'/'+star+'_likelihood.pdf'
  print 'Creating ', fname
  plt.savefig(fname,bbox_inches='tight')
  plt.close()
##===========================================================
#              plot tr fancy function
#===========================================================

#Ntransit is the number of the transit that we want to plot
def fancy_tr_plot(t0_val,xtime,yflux,errors,xmodel,xmodel_res,fd_reb,res_res,fname):

  print 'Creating ', fname
  #Do the plot
  tfc = 24. # time factor conversion to hours
  local_T0 = t0_val
  plt.figure(1,figsize=(fsx,fsy))
  #Plot the transit light curve
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3.0, 1.])
  gs.update(hspace=0.00)
  ax1 = plt.subplot(gs[0])
  plt.tick_params(labelsize=fos)
  x_lim = (min(xtime)-local_T0)*tfc
  plt.xlim(x_lim,-x_lim)
  if ( select_y_tr ):
    plt.ylim(y_lim_min,y_lim_max)
  min_val_model = max(fd_reb) -  min(fd_reb)
  if ( plot_tr_errobars ):
    plt.errorbar((xtime-local_T0)*tfc,yflux,errors,fmt='ro',alpha=0.8)
  else:
    plt.plot((xtime-local_T0)*tfc,yflux,'ro',ms=8,alpha=1.0)
  plt.plot((xmodel-local_T0)*tfc,fd_reb,'k',linewidth=2.0,alpha=1.0)
  plt.ylabel('Relative flux',fontsize=fos)
  plt.xticks( np.arange(int(x_lim),int(-x_lim)+1,1))
  plt.minorticks_on()
  plt.ticklabel_format(useOffset=False, axis='y')
  plt.tick_params( axis='x',which='both',labelbottom='off')
  #Plot the residuals
  dplot = plt.subplot(gs[1])
  plt.tick_params(labelsize=fos)
  if ( plot_tr_errobars ):
    plt.errorbar((xmodel_res-local_T0)*tfc,res_res,errors,fmt='ro',alpha=0.8)
  else:
    plt.plot((xmodel_res-local_T0)*tfc,res_res,'ro',ms=8,alpha=1.0)
  plt.plot([x_lim,-x_lim],[0.0,0.0],'k--',linewidth=1.0,alpha=1.0)
  yylims = dplot.get_ylim()
  plt.yticks(np.arange(yylims[0],yylims[1],(yylims[1]-yylims[0])/4.))
  plt.xticks( np.arange(int(x_lim),int(-x_lim)+1,1))
  plt.xlim(x_lim,-x_lim)
  if ( select_y_tr ):
    plt.ylim( - ( y_lim_max - 1.0),y_lim_max - 1.0 )
  #Plot the residuals
  plt.ylabel('Residuals',fontsize=fos*0.75)
  plt.xlabel("T - T0 (hours)",fontsize=fos)
  plt.minorticks_on()
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight')
  plt.close()


#===========================================================
#              plot the folded data for each transit
#===========================================================
def plot_transit_nice():
#Move all the points to T0
  ldc = [ best_value(params[4+8*nplanets],get_value), best_value(params[5+8*nplanets],get_value) ]
  q1_val = ldc[0]
  q2_val = ldc[1]
  u1_val = np.sqrt(q1_val)
  u2_val = u1_val * (1.0 -2.0*q2_val)
  u1_val = 2.0*u1_val*q2_val
  flag = [is_log_P, is_ew, is_b_factor, is_log_a]

  t0_val = [None]*nplanets
  P_val  = [None]*nplanets
  e_val  = [None]*nplanets
  w_val  = [None]*nplanets
  i_val  = [None]*nplanets
  a_val  = [None]*nplanets
  rp_val = [None]*nplanets

  base = 4
  for m in range(0,nplanets):
    t0_val[m] = best_value(params[base + 0],get_value)
    P_val[m]  = best_value(params[base + 1],get_value)
    e_val[m]  = best_value(params[base + 2],get_value)
    w_val[m]  = best_value(params[base + 3],get_value)
    i_val[m]  = best_value(params[base + 4],get_value)
    a_val[m]  = best_value(params[base + 5],get_value)
    rp_val[m] = best_value(params[base + 6],get_value)
    base = base + 8

  for o in range(0,nplanets):

    if ( fit_tr[o] ):


      if ( len(span_tr) < 1 ):
          local_span = 0.0
      else:
          local_span = span_tr[o]

      local_time_d, xtime_d, yflux_d, eflux_d = create_transit_data(megax,megay,megae,o,local_span)

      local_time = np.concatenate(local_time_d)
      xtime = np.concatenate(xtime_d)
      yflux = np.concatenate(yflux_d)
      eflux = np.concatenate(eflux_d)

      xmodel_res = xtime
      xmodel = np.arange(min(xtime), max(xtime),1.0/20./24.)
      #Let us create the model

      xd_ub = np.ndarray(shape=(len(xmodel),n_cad))
      xd_ub_res = np.ndarray(shape=(len(xmodel_res),n_cad))
      zd_ub = [None]*len(xmodel)
      zd_ub_res = [None]*len(xmodel_res)
      fd_ub = [None]*len(xmodel)
      fd_ub_res = [None]*len(xmodel_res)
      #Use the long cadence data
      for m in range(0,len(xmodel)):
        #This vector has the model fit
        for n in range(0,n_cad):
          xd_ub[m][n] = xmodel[m] + t_cad * ( (n+1) - 0.5 * (n_cad + 1 ))/n_cad
        #This vector has the residuals
      for m in range(0,len(xmodel_res)):
        for n in range(0,n_cad):
          xd_ub_res[m][n] = xmodel_res[m] + t_cad * ( (n+1) - 0.5 * (n_cad + 1))/n_cad

    #Calculate the transit curve for all the data
      for m in range(0,len(xmodel)):
        zd_ub[m] = pti.find_z(xd_ub[m][:],[0.0,P_val[o],e_val[o],w_val[o],i_val[o],a_val[o]],flag)
        fd_ub[m], dummm = pti.occultquad(zd_ub[m],u1_val,u2_val,rp_val[o])
      for m in range(0,len(xmodel_res)):
        zd_ub_res[m] = pti.find_z(xd_ub_res[m][:],[0.0,P_val[o],e_val[o],w_val[o],i_val[o],a_val[o]],flag)
        fd_ub_res[m], dummm = pti.occultquad(zd_ub_res[m],u1_val,u2_val,rp_val[o])

      #Bin the data
      fd_reb = [0.0]*len(xmodel)
      fd_reb_res = [0.0]*len(xmodel_res)
      for m in range(0,len(xmodel)):
        for n in range(0,n_cad):
          fd_reb[m] = fd_reb[m] + fd_ub[m][n]/n_cad
      for m in range(0,len(xmodel_res)):
        for n in range(0,n_cad):
          #This is the flux caused by the time stams which come from the data
          fd_reb_res[m] = fd_reb_res[m] + fd_ub_res[m][n]/n_cad

     #############################################################################
     # Let us calculate the flux caused by the other planets
#      xmodel_vec = np.concatenate(list(local_time))
      xmodel_vec = local_time
      xd_ub_vec = np.ndarray(shape=(len(xmodel_vec),n_cad))
      for m in range(0,len(xmodel_vec)):
        for n in range(0,n_cad):
          xd_ub_vec[m][n] = xmodel_vec[m] + t_cad * ( (n+1) - 0.5 * (n_cad + 1))/n_cad

      fd_ub_dum = [None]*len(xmodel_vec)
      fd_ub_total = [0.0]*len(xmodel_vec)
      zd_ub_dum = [None]*len(xmodel_vec)
      for p in range(0,nplanets):
        fd_reb_dum  = [0.0]*len(xmodel_vec)
        if ( p != o ):
          for m in range(0,len(xmodel_vec)):
            zd_ub_dum[m] = pti.find_z(xd_ub_vec[m][:],[t0_val[p],P_val[p],e_val[p],w_val[p],i_val[p],a_val[p]],flag)
            fd_ub_dum[m], dummm = pti.occultquad(zd_ub_dum[m],u1_val,u2_val,rp_val[p])
            for n in range(0,n_cad):
              #This is the flux caused by the time stams which come from the data
              fd_reb_dum[m] = fd_reb_dum[m] + fd_ub_dum[m][n]/n_cad
            fd_ub_total[m] = fd_ub_total[m] + fd_reb_dum[m]

      yflux_local = yflux - fd_ub_total
      yflux_local = yflux_local - 1 + nplanets
      #The flux has been corrected for the other planets

      res_res = yflux_local - fd_reb_res

     #############################################################################

      fname = outdir+'/'+star+plabels[o]+'_tr.pdf'
      #xtime is the folded time
      #yflux is the data flux
      #eflux is the error related to yflux
      #xmodel is the light curve model timestamps
      #xmodel_res is the residuals time_stamps
      #fd_reb is the modeled light cuve
      #res_res are the residuals
      fancy_tr_plot(0.0,xtime,yflux_local,eflux,xmodel,xmodel_res,fd_reb,res_res,fname)

#===========================================================
#              plot all transits
#===========================================================

def plot_all_transits():

  ldc = [ best_value(params[4+8*nplanets],get_value), best_value(params[5+8*nplanets],get_value) ]
  q1_val = ldc[0]
  q2_val = ldc[1]
  u1_val = np.sqrt(q1_val)
  u2_val = u1_val * (1.0 -2.0*q2_val)
  u1_val = 2.0*u1_val*q2_val
  flag = [is_log_P, is_ew, is_b_factor, is_log_a]

  t0_val = [None]*nplanets
  P_val  = [None]*nplanets
  e_val  = [None]*nplanets
  w_val  = [None]*nplanets
  i_val  = [None]*nplanets
  a_val  = [None]*nplanets
  rp_val = [None]*nplanets

  base = 4
  for m in range(0,nplanets):
    t0_val[m] = best_value(params[base + 0],get_value)
    P_val[m]  = best_value(params[base + 1],get_value)
    e_val[m]  = best_value(params[base + 2],get_value)
    w_val[m]  = best_value(params[base + 3],get_value)
    i_val[m]  = best_value(params[base + 4],get_value)
    a_val[m]  = best_value(params[base + 5],get_value)
    rp_val[m] = best_value(params[base + 6],get_value)
    base = base + 8

  xt = [None]*nplanets
  dt = [None]*nplanets
  yt = [None]*nplanets
  et = [None]*nplanets
  for i in range(0,nplanets):
    if ( fit_tr[i] ):

      xt[i], dt[i], yt[i], et[i] = create_transit_data(megax,megay,megae,i)
      if ( is_plot_all_tr[i] ):
        for j in range(0,len(xt[i])):

          xvec = np.array(xt[i][j])
          xvec_model = np.arange(min(xvec),max(xvec),1./20./24.)

          xd_ub = np.ndarray(shape=(len(xvec),n_cad))
          xd_ub_res = np.ndarray(shape=(len(xvec_model),n_cad))
          zd_ub = [None]*len(xvec)
          zd_ub_res = [None]*len(xvec_model)
          fd_ub = [None]*len(xvec)
          fd_ub_res = [None]*len(xvec_model)
          fd_ub_total = [0.0]*len(xvec)
          fd_ub_res_total = [0.0]*len(xvec_model)
          fd_reb = [0.0]*len(xvec)
          fd_reb_res = [0.0]*len(xvec_model)

          for m in range(0,len(xvec)):
          #This vector has the model fit
            for n in range(0,n_cad):
              xd_ub[m][n] = xvec[m] + t_cad * ( (n+1) - 0.5 * (n_cad + 1 ))/n_cad
              #This vector has the residuals
          for m in range(0,len(xvec_model)):
            for n in range(0,n_cad):
              xd_ub_res[m][n] = xvec_model[m] + t_cad * ( (n+1) - 0.5 * (n_cad + 1))/n_cad

          #Calculate the transit curve for all the planets
          for o in range(0,nplanets):
            for m in range(0,len(xvec)):
              zd_ub[m] = pti.find_z(xd_ub[m][:],[t0_val[o],P_val[o],e_val[o],w_val[o],i_val[o],a_val[o]],flag)
              fd_ub[m], dummm = pti.occultquad(zd_ub[m],u1_val,u2_val,rp_val[o])
              for n in range(0,n_cad):
                fd_reb[m] = fd_reb[m] + fd_ub[m][n]/n_cad
              fd_ub_total[m] = fd_ub_total[m] + fd_reb[m]
              fd_reb[m] = 0.0
            for m in range(0,len(xvec_model)):
              zd_ub_res[m] = pti.find_z(xd_ub_res[m][:],[t0_val[o],P_val[o],e_val[o],w_val[o],i_val[o],a_val[o]],flag)
              fd_ub_res[m], dummm = pti.occultquad(zd_ub_res[m],u1_val,u2_val,rp_val[o])
              for n in range(0,n_cad):
                fd_reb_res[m] = fd_reb_res[m] + fd_ub_res[m][n]/n_cad
              fd_ub_res_total[m] = fd_ub_res_total[m] + fd_reb_res[m]
              fd_reb_res[m] = 0.0

          yflux_local  = list(yt[i][j])
          fd_ub_total = np.array(fd_ub_total) - nplanets + 1
          fd_ub_res_total = np.array(fd_ub_res_total) - nplanets + 1
          res_res = np.array(yflux_local) - np.array(fd_ub_total)

          fname = outdir+'/'+star+plabels[i]+'_transit'+str(j)+'.pdf'
          #xtime is the folded time
          #yflux is the data flux
          #eflux is the error related to yflux
          #xmodel is the light curve model timestamps
          #xmodel_res is the residuals time_stamps
          #fd_reb is the modeled light cuve
          #res_res are the residuals
          n = xvec[len(xvec)-1] - xt[i][0][0]
          n = int(n/P_val[i])
          fancy_tr_plot(t0_val[i]+P_val[i]*n,xvec,yflux_local,et[i][j],xvec_model,xvec,fd_ub_res_total,res_res,fname)

#===========================================================
#              plot rv fancy function
#===========================================================

def plot_rv_fancy(p_rv,rvy,p_all,rv_dum,errs_all,res,telescopes_labels,fname):
  print 'Creating ', fname
  plt.figure(3,figsize=(fsx,fsy))
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3., 1.])
  gs.update(hspace=0.00)
  ax0 = plt.subplot(gs[0])
  plt.tick_params(labelsize=fos)
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
    markersize=rv_markersize,fillstyle=rv_fillstyle)
  plt.legend(loc=0, ncol=1,scatterpoints=1,numpoints=1,frameon=False,fontsize='small')
  plt.xticks(np.arange(0.,1.01,0.1))
  plt.tick_params( axis='x',which='both',labelbottom='off')
  #plt.subplot(312)
  ax1 = plt.subplot(gs[1])
  plt.tick_params(labelsize=fos)
  plt.xlabel("Orbital phase",fontsize=fos)
  plt.tick_params( axis='x',which='minor',bottom='on',left='on',right='on',top='on')
  plt.xticks(np.arange(0.,1.01,0.1))
  plt.ylabel('Residuals (m/s)',fontsize=fos*0.75)
  plt.plot([0.,1.],[0.,0.],'k--',linewidth=1.0)
  for j in range(0,nt):
    plt.errorbar(p_all[j],res[j],errs_all[j],\
    label=telescopes_labels[j],fmt=mark[j],alpha=1.0,markersize=rv_markersize,fillstyle=rv_fillstyle)
  yylims = ax1.get_ylim()
  plt.yticks(np.arange(yylims[0],yylims[1],(yylims[1]-yylims[0])/4.))
  plt.minorticks_on()
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight')
  plt.close()


#===========================================================
#                   RV PLOTS
#===========================================================

v_vec_val = [None]*nt
v_val = [None]*nt
#3 + npars + ldc
v_vec_val[:] = params[4+8*nplanets+2:4+8*nplanets+2+nt]
for o in range(0,nt):
  v_val[o] = best_value(v_vec_val[o],get_value)
  if ( is_log_rv0 ):
    v_val[o] = 10.0**(v_val[o])


alpha_val = 0.0
beta_val = 0.0
if ( is_linear_trend != 'f' or is_quadratic_trend != 'f' ):
  alpha_val = best_value(params_trends[0],get_value)
  beta_val  = best_value(params_trends[1],get_value)


base = 4
t0_val = np.ndarray(nplanets)
P_val  = np.ndarray(nplanets)
e_val  = np.ndarray(nplanets)
w_val  = np.ndarray(nplanets)
k_val  = np.ndarray(nplanets)

for o in range(0,nplanets):
  t0_val[o] = best_value(params[base + 0],get_value)
  P_val[o]  = best_value(params[base + 1],get_value)
  e_val[o]  = best_value(params[base + 2],get_value)
  w_val[o]  = best_value(params[base + 3],get_value)
  k_val[o]  = best_value(params[base + 7],get_value)
  if ( is_log_P ):
    P_val[o] = 10.0**(P_val)
  if ( is_log_k ):
    k_val[o] = 10.0**(k_val)
  if ( is_ew ):
    edum_val = e_val[o]
    e_val[o] = e_val[o]**2 + w_val[o]**2
    w_val[o] = np.arctan2(edum_val,w_val[o])

  base = base + 8

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
    plt.xlabel("BJD - 2450000 (days)",fontsize=fos)
    plt.ylabel('RV (m/s)',fontsize=fos)
    plt.xlim(xmin,xmax)
    for j in range(0,nt):
      plt.errorbar(time_all[j],rv_dum[j],errs_datas[j],\
      label=telescopes_labels[j],fmt=mark[j],alpha=1.0,markersize=4)
    plt.legend(loc=0, ncol=1,scatterpoints=1,numpoints=1,frameon=False,fontsize='small')
    fname = outdir+'/'+star+'_rv_all.pdf'
    print 'Creating ', fname
    plt.savefig(fname,format='pdf',bbox_inches='tight')
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
        rv_dum.append(list(rv_all[j]))
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
      k_dum[i],P_val[i],e_val[i],w_val[i],0.0 \
      ,0.0)

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

      p_rv[i] = scale_period(rvx,t0_val[i],P_val[i])
      p_all = [None]*nt
      for j in range(0,nt):
        p_all[j] = scale_period(time_all[j],t0_val[i],P_val[i])

      fname = outdir+'/'+star+plabels[i]+'_rv.pdf'
      plot_rv_fancy(p_rv[i],rvy[i],p_all,rv_dum,errs_all,res,telescopes_labels,fname)

#===========================================================
#                   Histogram plots
#===========================================================

def create_plot_histogram(params,plabs,cbars='red',nb=50,num=[]):
  if ( len(num) < 2 ):
    n = range(0,len(params))
  else:
    n = num
  plt.figure(1,figsize=(12,4*(len(n))/2))
  gs = gridspec.GridSpec(nrows=(len(n)+1)/2,ncols=2)
  j = 0
  for i in n:
    plt.subplot(gs[j])
    vpar, lpar, rpar = find_vals_perc(params[i],1.0)
    moda = my_mode(params[i])
    #best_val = params[i][minchi2_index]
    #plt.axvline(x=best_val,c='yellow')
    plt.axvline(x=vpar,c=cbars)
    plt.axvline(x=moda,c='y',ls='-.')
    plt.axvline(x=vpar-lpar,c=cbars,ls='--')
    plt.axvline(x=vpar+rpar,c=cbars,ls='--')
    plt.xlabel(plabs[i])
    plt.hist(params[i],normed=True,bins=nb)
    j = j + 1

  fname = outdir+'/'+star+'_histogram.pdf'
  print 'Creating ', fname
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.close()

def plot_histogram():
    labs = []
    for o in range(0,nplanets):
      etiquetas = ['$T0$'+plabels[o]+' (days)','$P$'+plabels[o]+' (days)','$e$'+plabels[o], \
                 '$\omega$'+plabels[o],'$b$'+plabels[o],'$a/R_\star$'+plabels[o], \
                 '$R_{\mathrm{p}}/R_\star$'+plabels[o],'$k$'+plabels[o]+' (kms$^{-1}$)']
      labs.append(etiquetas)
    labs.append(['$q_1$','$q_2$'])
    labs.append(telescopes_labels)
    labels = np.concatenate(labs)
    create_plot_histogram(params[4:],labels, cbars='red', nb=50, num=plot_parameters)

#===========================================================
#                   Correlation plots
#===========================================================

def create_plot_correlation(params,plabs,col='red',mark='.',num=[]):
  if ( len(num) < 1 ):
    n = range(0,len(params))
  else:
    n = num
  plt.figure(1,figsize=(4*len(n),4*len(n)))
  gs = gridspec.GridSpec(nrows=len(n),ncols=len(n))
  o = 0
  for i in n:
    p = 0
    for j in n:
      if ( j < i ):
        plt.subplot(gs[o*len(n)+p])
        plt.tick_params( axis='y',which='both',labelleft='off')
        plt.tick_params( axis='x',which='both',labelbottom='off')
        plt.ticklabel_format(useOffset=False, axis='both')
        if ( j == n[0] ):
           plt.ylabel(plabs[i],fontsize=25)
        elif ( j == i - 1 ):
          plt.tick_params( axis='y',which='both',labelleft='off')
          plt.tick_params( axis='x',which='both',labelbottom='off')
        else:
          plt.tick_params( axis='y',which='both',labelleft='off')
          plt.tick_params( axis='x',which='both',labelbottom='off')
        if ( i == n[len(n)-1]):
          plt.xlabel(plabs[j],fontsize=25)
        else:
          plt.tick_params( axis='y',which='both',labelleft='off')
          plt.tick_params( axis='x',which='both',labelbottom='off')
        plt.hist2d(params[j],params[i],bins=100,norm=LogNorm())
        p = p + 1
    o = o + 1

  fname = outdir+'/'+star+'_correlations.pdf'
  print 'Creating ', fname
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.close()

def plot_correlations():
  labs = []
  for o in range(0,nplanets):
    etiquetas = ['$T0$'+plabels[o],'$P$'+plabels[o],'$e$'+plabels[o], \
                 '$\omega$'+plabels[o],'$b$'+plabels[o],'$a/R_\star$'+plabels[o], \
                 '$R_{\mathrm{p}}/R_\star$'+plabels[o],'$k$'+plabels[o]]
    labs.append(etiquetas)
  labs.append(['$q_1$','$q_2$'])
  labs.append(telescopes_labels)
  labels = np.concatenate(labs)
  create_plot_correlation(params[4:],labels,col='blue',num=plot_parameters)

