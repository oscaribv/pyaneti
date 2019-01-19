#-----------------------------------------------------------------
#       TRANSIT PARAMERS TO BE USED TO GENERATE PLOTS
#-----------------------------------------------------------------

#Create parameters vector
pars_tr = np.zeros(shape=(nplanets,6))
for m in range(0,nplanets):
    pars_tr[m][0] = tp_val[m]
    pars_tr[m][1] = P_val[m]
    pars_tr[m][2] = e_val[m]
    pars_tr[m][3] = w_val[m]
    pars_tr[m][4] = i_val[m]
    pars_tr[m][5] = a_val[m]

#Returns the pars_tr array for a given chain number
def pars_tr_chain(params,nchain):
  ldc = [ params[4+8*nplanets][nchain], params[5+8*nplanets][nchain]]
  q1_val = ldc[0]
  q2_val = ldc[1]
  u1_val = np.sqrt(q1_val)
  u2_val = u1_val * (1.0 -2.0*q2_val)
  u1_val = 2.0*u1_val*q2_val
  flag = [False]*4
  my_ldc = [u1_val,u2_val]

  t0_val = np.ndarray(nplanets)
  tp_val = np.ndarray(nplanets)
  P_val  = np.ndarray(nplanets)
  e_val  = np.ndarray(nplanets)
  w_val  = np.ndarray(nplanets)
  i_val  = np.ndarray(nplanets)
  a_val  = np.ndarray(nplanets)
  rp_val = np.ndarray(nplanets)
  tp_val = np.ndarray(nplanets)


  base = 4
  for m in range(0,nplanets):
    t0_val[m] = params[base + 0][nchain]
    P_val[m]  = params[base + 1][nchain]
    e_val[m]  = params[base + 2][nchain]
    w_val[m]  = params[base + 3][nchain]
    i_val[m]  = params[base + 4][nchain]
    a_val[m]  = params[base + 5][nchain]
    rp_val[m] = params[base + 6][nchain]
    base = base + 8


    if ( is_den_a ):
      a_val[m] = a_val[0]*(P_val[m]*P_val[m]*7464960000.*G_cgs/3.0/np.pi)**(1./3.)
      if ( m > 0):
        a_val[m] = (P_val[m]/P_val[0])**(2./3.)*a_val[0]

    #Check flags
    #Change between b and i
    if ( is_b_factor ):
      i_val[m] = np.arccos( i_val[m] / a_val[m] * \
              ( 1.0 + e_val[m] * np.sin(w_val[m] ) / ( 1.0 - e_val[m]**2 ) ) )

    if ( is_ew ):
      e_dum = e_val[m]
      w_dum = w_val[m]
      e_val[m] = e_val[m]**2 + w_val[m]**2
      w_val[m] = np.arctan2(e_dum,w_val[m])
      w_val[m] = w_val[m] % (2*np.pi)

    tp_val[m] = pti.find_tp(t0_val[m],e_val[m],w_val[m],P_val[m])

    #Create parameters vector
    pars_tr = np.zeros(shape=(7,nplanets))
    for m in range(0,nplanets):
      pars_tr[0][m] = tp_val[m]
      pars_tr[1][m] = P_val[m]
      pars_tr[2][m] = e_val[m]
      pars_tr[3][m] = w_val[m]
      pars_tr[4][m] = i_val[m]
      pars_tr[5][m] = a_val[m]
      pars_tr[6][m] = rp_val[m]

  return pars_tr


#===========================================================
#===========================================================

def create_folded_tr_plots():

  for o in range(0,nplanets):

    if (fit_tr[o]):

      tr_vector = [None]*nbands

      for m in range(0,nbands):
        localx = []
        localy = []
        locale = []
        localt = []
        for n in range(0,len(megax)):
          if ( trlab[n] == m ):
              localx.append(megax[n])
              localy.append(megay[n])
              locale.append(megae[n])
              localt.append(0)
        tr_vector[m] = plot_parameters_tr(localx,localy,locale,localt,pars_tr,rp_val,o,m)

      transpose_tr = np.asarray(tr_vector)
      transpose_tr = transpose_tr.transpose()
      fancy_tr_plot(transpose_tr, o)


#Ntransit is the number of the transit that we want to plot
#tr_vector contains:
#xtime,yflux,eflux,xmodel,xmodel_res,fd_ub,res_res,fd_ub_unbinned
#these lists contains all the information for a given label
def fancy_tr_plot(tr_vector,pnumber):

  fname = outdir+'/'+star+plabels[pnumber]+'_tr.pdf'
  print 'Creating ', fname
  #Do the plot
  tfc = 24. # time factor conversion to hours
  local_T0 = 0.

  #Extract the vectors to be plotted from tr_vector
  xtime = tr_vector[0]
  yflux = tr_vector[1]
  eflux = tr_vector[2]
  xmodel= tr_vector[3]
  xmodel_res = tr_vector[4]
  fd_ub = tr_vector[5]
  res_res = tr_vector[6]
  fd_ub_unbinned = tr_vector[7]

  #Start the plot
  plt.figure(1,figsize=(fsx,fsy))
  #Plot the transit light curve
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3.0, 1.])
  gs.update(hspace=0.00)
  ax1 = plt.subplot(gs[0])
  plt.tick_params(labelsize=fos,direction='in')
  x_lim = (min(np.concatenate(xtime))-local_T0)*tfc
  plt.xlim(x_lim,-x_lim)
  if ( select_y_tr ):
    plt.ylim(y_lim_min,y_lim_max)
  min_val_model = max(np.concatenate(fd_ub)) -  min(np.concatenate(fd_ub))
  deltay = 0.
  dy = max(rp_val[pnumber*nbands:(pnumber+1)*nbands])**2
  for m in range(0,nbands):
    if ( plot_tr_errorbars  ):
      plt.errorbar((xtime-local_T0)*tfc,yflux,errors,color=tr_colors,fmt='.',alpha=1.0)
    else:
      plt.plot((xtime[m]-local_T0)*tfc,yflux[m]-deltay,'o',ms=7,alpha=0.8)
      ###plt.errorbar(-x_lim*(0.95),y0 +eflux[m][0]*1.5-deltay,eflux[m][0],color=tr_colors,ms=7,fmt='o',alpha=1.0)
      ###plt.annotate('Error bar',xy=(-x_lim*(0.70),y0 +eflux[m][0]*1.65-deltay),fontsize=fos*0.7)
      #if ( is_special ):
      #  [plt.fill_between((xmodel-local_T0)*tfc,*flux_model[i:i+2,:],alpha=0.3,facecolor='b') for i in range(1,6,2)]
      #  plt.fill_between((xmodel-local_T0)*tfc,*flux_model[5:7,:],alpha=0.5,facecolor='k')
      #if (plot_unbinned_model):
      #  plt.plot((xmodel-local_T0)*tfc,fd_ub_unbinned,'b',linewidth=2.0,alpha=1.0)
      plt.plot((xmodel[m]-local_T0)*tfc,fd_ub[m]-deltay,'k',linewidth=2.0,alpha=1.0)
      deltay = deltay + dy
  y0,yyyy = ax1.get_ylim()
  ##plot binned data
  #plt.plot(np.asarray(xbined)*tfc,fbined,'ko')
  ##
  if (nbands == 1 ): plt.ylabel('Flux',fontsize=fos)
  if (nbands > 1 ): plt.ylabel('Flux + offset',fontsize=fos)
  #Calculate the optimal step for the plot
  step_plot = int(abs(x_lim)) #the value of the x_axis
  step_plot = step_plot + int( step_plot % 2 ) # now we ensure the result is par
  step_plot = int ( step_plot / 8. ) + 1 #The size of the jump depends
  #let us get the new limit
  nuevo = np.arange(0,int(abs(x_lim)) + step_plot ,step_plot)
  mxv = np.max(nuevo)
#  plt.xticks( np.arange(int(x_lim),int(-x_lim)+1,step_plot))
  plt.xticks( np.arange(-mxv,mxv+step_plot,step_plot))
  plt.minorticks_on()
  plt.ticklabel_format(useOffset=False, axis='y')
  plt.xlim(x_lim,-x_lim)
  plt.tick_params( axis='x',which='both',direction='in',labelbottom=False)
  plt.tick_params( axis='y',which='both',direction='in')
  #------------------------------------------------------------
  #Plot the residuals
  #------------------------------------------------------------
  ax0 = plt.subplot(gs[1])
  plt.tick_params(labelsize=fos,direction='in')
  for m in range(nbands):
    if ( plot_tr_errorbars  ):
      plt.errorbar((xmodel_res[m]-local_T0)*tfc,res_res[m]*1e6,eflux[m]*1e6,fmt='.',alpha=0.5)
    else:
      plt.plot((xmodel_res[m]-local_T0)*tfc,res_res[m]*1e6,'o',ms=7,alpha=0.5)
  plt.plot([x_lim,-x_lim],[0.0,0.0],'k--',linewidth=1.0,alpha=1.0)
  plt.xticks( np.arange(-mxv,mxv+step_plot,step_plot))
  plt.xlim(x_lim,-x_lim)
  yylims = ax0.get_ylim()
  miy = (max(abs(yylims[0]),abs(yylims[1])))
  plt.yticks(np.arange(-miy,miy,miy/2.))
  plt.ylim(-miy,miy*1.25)
  #Calcualte the rms
  if ( is_plot_std_tr ):
    trsigma = np.std(res_res*1e6,ddof=1)
    trsstr = ('%4.0f ppm'%(trsigma))
    y0,yyyy = ax0.get_ylim()
    plt.annotate('$\sigma = $'+trsstr,xy=(x_lim*(0.80),y0 + 1.8*miy),fontsize=fos*0.7)
#  if ( select_y_tr ):
#    plt.ylim( - ( y_lim_max - 1.0),y_lim_max - 1.0 )
  #Plot the residuals
  plt.minorticks_on()
  plt.tick_params( axis='x',which='both',direction='in')
  plt.tick_params( axis='y',which='both',direction='in')
  plt.ylabel('Residuals (ppm)',fontsize=fos*0.75)
  plt.xlabel("T - T0 (hours)",fontsize=fos)
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight',dpi=300)
  plt.close()

def plot_parameters_tr(time,flujo,eflujo,trlab,pars_tr,rp,plabel,bandlab):

    if ( len(span_tr) < 1 ):
        local_span = 0.0
    else:
        local_span = span_tr[plabel]

    local_time_d, xtime_d, yflux_d, eflux_d = create_transit_data(time,flujo,eflujo,plabel,local_span)

    local_time = np.concatenate(local_time_d)
    xtime = np.concatenate(xtime_d)
    yflux = np.concatenate(yflux_d)
    eflux = np.concatenate(eflux_d)

    xmodel_res = xtime
    mimax = abs(max(abs(min(xtime)),abs(max(xtime))))
    xmodel = np.arange(-mimax, mimax,1.0/20./24.)
    newtrlab = [0]*len(xmodel)
    #Let us create the model

    #The model has T0 = 0
    dumtp = pti.find_tp(0.0,e_val[plabel],w_val[plabel],P_val[plabel])
    dparstr = np.concatenate([[dumtp],pars_tr[plabel][1:]])
    #fd_ub = pti.flux_tr(xmodel,dparstr,my_ldc,n_cad,t_cad)
    fd_ub = pti.flux_tr(xmodel,newtrlab,dparstr,rp[plabel*nbands+bandlab],my_ldc[bandlab*2:bandlab*2+2],n_cad,t_cad)
    #Let us create an unbinned model plot
    #fd_ub_unbinned = pti.flux_tr(xmodel,dparstr,my_ldc,1,t_cad)
    fd_ub_unbinned = pti.flux_tr(xmodel,newtrlab,dparstr,rp[plabel*nbands+bandlab],my_ldc[bandlab*2:bandlab*2+2],1,t_cad)
    #Calculate the flux to copute the residuals
    #fd_ub_res = pti.flux_tr(xmodel_res,dparstr,my_ldc,n_cad,t_cad)
    newtrlab=[0]*len(xmodel_res)
    fd_ub_res = pti.flux_tr(xmodel_res,newtrlab,dparstr,rp[plabel*nbands+bandlab],my_ldc[bandlab*2:bandlab*2+2],n_cad,t_cad)

    #Define a vector which will contain the data of other planers for multi fits
    fd_ub_total = list(fd_ub_res)
    fd_ub_total = np.zeros(shape=len(fd_ub_res))

   #############################################################################
   # Let us calculate the flux caused by the other planets
    for p in range(0,nplanets):
      if ( p != plabel ):
        #fd_ub_total stores the flux of a star for each independent
        #fd_ub_total = fd_ub_total + pti.flux_tr(local_time,pars_tr[:,p],my_ldc,n_cad,t_cad)
        newtrlab=[0]*len(local_time)
        fd_ub_total = fd_ub_total + pti.flux_tr(local_time,newtrlab,pars_tr[p],rp[p*nbands+bandlab],my_ldc[bandlab*2:bandlab*2+2],n_cad,t_cad)

    #Remove extra planets from the data
    yflux_local = yflux - fd_ub_total
    yflux_local = yflux_local - 1 + nplanets
    #The flux has been corrected for the other planets

    #Get the residuals
    res_res = yflux_local - fd_ub_res

   #############################################################################


#    #Let us make a binning
#    nbins  = 20
#    fbined = [0.]*nbins
#    xbined = [0.]*nbins
#    deltax = (max(xtime) - min(xtime))/nbins
#    mim = min(xtime)
#    mix = mim + deltax
#    mindex = []
#    for l in range(0,20):
#        for k in range(0,len(xtime)):
#            if (xtime[k] > mim and xtime[k] < mix):
#                mindex.append(k)
#        xbined[l] = np.mean(xtime[mindex])
#        fbined[l] = np.mean(yflux_local[mindex])
#        mindex = []
#        mim = mix
#        mix = mix + deltax


#    fname = outdir+'/'+star+plabels[o]+'_tr.pdf'
    #xtime is the folded time
    #yflux_local is the data flux
    #eflux is the error related to yflux
    #xmodel is the light curve model timestamps
    #xmodel_res is the residuals time_stamps
    #fd_reb is the modeled light cuve
    #res_res are the residuals
    #fd_ub_unbinned is the unbinned model
    return xtime,yflux_local,eflux,xmodel,xmodel_res,fd_ub,res_res,fd_ub_unbinned

#===========================================================
#              plot all transits
#===========================================================

#Now this functions works only with one band
def plot_all_transits():
  global plot_tr_errorbars

  #Create the plot of the whole light
  model_flux = pti.flux_tr(megax,pars_tr,my_ldc,n_cad,t_cad)
  res_flux = megay - model_flux

  for i in range(0,nplanets):
    if ( fit_tr[i] ):

      xt, dt, yt, et = create_transit_data(megax,megay,megae,i,span_tr[i])
      xt2, dt2, yt2, et2 = create_transit_data(megax,res_flux,megae,i,span_tr[i])

      if ( is_plot_all_tr[i] ):
        for j in range(0,len(xt)):
          xtm = np.arange(min(xt[j]),max(xt[j]),1./20./24.)
          ytm = pti.flux_tr(xtm,trlab,pars_tr,rps,my_ldc,n_cad,t_cad)

          fname = outdir+'/'+star+plabels[i]+'_transit'+str(j)+'.pdf'
          n = xt[j][len(xt[j])-1] - xt[0][0]
          n = int(n/P_val[i])
          #is_err = plot_tr_errorbars
          #plot_tr_errorbars = True
          fancy_tr_plot(t0_val[i]+P_val[i]*n,xt[j],yt[j],et[j],xtm,xt2[j],ytm,np.array(yt2[j]),ytm,fname)
          #plot_tr_errorbars = is_err


#===========================================================
#              clean transits
#  This function cleans the light curve with a N-sigma aogorithm
#===========================================================

def clean_transits(sigma=10):

    #Now we are ready to call the function in fortran
    #All the data is in megax, megay and megae
    model_flux = pti.flux_tr(megax,trlab,pars_tr,rps,my_ldc,n_cad,t_cad)
    xvec_model = np.arange(min(megax),max(megax),1./20./24.)
    solution_flux = pti.flux_tr(xvec_model,trlab,pars_tr,rps,my_ldc,n_cad,t_cad)

    #Calcualte the residuals
    res_flux = megay - model_flux

    #Call the sigma clipping functions
    new_t, new_f = sigma_clip(megax,megay,res_flux,limit_sigma=sigma)

    #Recalculate the error bars
    new_model_flux = pti.flux_tr(new_t,trlab,pars_tr,rps,my_ldc,n_cad,t_cad)
    #New residuals
    new_res_flux = new_f - new_model_flux
    #Recompute the error bars from the std of the residuals
    new_err = np.std(new_res_flux,ddof=1)


    plt.figure(3,figsize=(5*fsx,fsy))
    fname = outdir+'/'+star+'_lightcurve.pdf'
    print 'Creating ', fname
    plt.plot(megax,megay,'ro')
    plt.plot(xvec_model,solution_flux,'k-')
    plt.plot(new_t,new_f,'bo')
    plt.ylabel('Flux',fontsize=fos*0.75)
    plt.xlabel(tr_xlabel,fontsize=fos)
    plt.xlim(min(megax),max(megax))
    plt.minorticks_on()
    plt.savefig(fname,format='pdf',bbox_inches='tight')
    plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight',dpi=300)
    plt.close()

    #Write the cleaned light curve into a file
    #Let us create or detrended file
    out_f = outdir+'/'+star+'_new_lc.dat'
    of = open(out_f,'w')
    for i in range(0,len(new_t)):
      of.write(' %8.8f   %8.8f  %8.8f \n'%(new_t[i],new_f[i],new_err))

    of.close()

