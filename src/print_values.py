from matplotlib import gridspec

#Read the data
#Dummy params vector contains
#[0] -> i
#[1] -> chain label
#[2] -> chi2_rv
#[3] -> chi2_tr
#[4-8*nplanets] -> parameters
#[4+8*nplanets-+2] -> ldc
#[4+8*nplanets+2-+ntelescopes] -> rvs
newfile = outdir+'/'+star+'_all_data.dat'
dparams = np.loadtxt(newfile, comments='#',unpack=True)


#Let us do the clustering
params = list(dparams)
#par_likelihood = list(ldparams[2])
params_jitter = [0.0]*2
new_nwalkers = nwalkers
#The maximum log(likelihood)
dmaxloglike = dparams[1]
maxloglike = dmaxloglike
if ( is_clustering ):
  #Starting clustering
  good_index, new_nwalkers = good_clustering_fast(dparams[2]+dparams[3],nconv,nwalkers)
  maxloglike = clustering_fast(dmaxloglike,good_index,nconv)
  for o in range(0,len(dparams)):
    params[o] = clustering_fast(dparams[o],good_index,nconv)

if ( is_jitter_rv or is_jitter_tr ):
  newfile_jitter = outdir+'/'+star+'_jitter_data.dat'
  dparams_jitter = np.loadtxt(newfile_jitter, comments='#',unpack=True)
  params_jitter = list(dparams_jitter)
  if ( is_clustering ):
    for o in range(0,n_jrv+1):
      params_jitter[o] = clustering_fast(dparams_jitter[o],good_index,nconv)

if ( is_linear_trend != 'f' or is_quadratic_trend != 'f' ):
  newfile_trends = outdir+'/'+star+'_trends_data.dat'
  dparams_trends = np.loadtxt(newfile_trends, comments='#',unpack=True)
  params_trends = list(dparams_trends)
  if ( is_clustering ):
    for o in range(0,2):
      params_trends[o] = clustering_fast(dparams_trends[o],good_index,nconv)

#Create the stellar data
mstar = np.random.normal(loc=mstar_mean,scale=mstar_sigma,size=new_nwalkers*nconv)
rstar = np.random.normal(loc=rstar_mean,scale=rstar_sigma,size=new_nwalkers*nconv)
tstar = np.random.normal(loc=tstar_mean,scale=tstar_sigma,size=new_nwalkers*nconv)

#Calculate the BIC
ndata = 0
if ( total_rv_fit ):
 ndata = ndata + len(mega_rv)
if ( total_tr_fit ):
 ndata = ndata + len(megax)

npars = 0

if ( len(plot_parameters) < 2):

  #plot_parameters stores the flag of each fitted parameter
  plot_parameters = []

  for o in range(0,len(fit_all)):
      if ( fit_all[o] != 'f' ):
          npars = npars + 1
          plot_parameters.append(o)

  for o in range(0,len(fit_ldc)):
      if ( fit_ldc[o] != 'f' ):
          npars = npars + 1
          plot_parameters.append(o+8*nplanets)

  for o in range(0,len(telescopes)):
      if ( fit_v0 != 'f' ):
          npars = npars + 1
          plot_parameters.append(o+2+8*nplanets)

if ( is_linear_trend != 'f' ):
  npars = npars + 1
if ( is_quadratic_trend != 'f' ):
  npars = npars + 1
if ( is_jitter_rv ):
  npars = npars + 1*n_jrv
if ( is_jitter_tr ):
  npars = npars + 1

dummy_pars = [0.0]*len(params)
for o in range(0,len(params)):
    dummy_pars[o] = best_value(params[o],maxloglike,get_value)

fit_pars = dummy_pars[4:4+8*nplanets]
ldc_pars = dummy_pars[4+8*nplanets:6+8*nplanets]
rvs_pars = dummy_pars[6+8*nplanets:6+8*nplanets+nt]

fit_jrv = [0.0]*nt
fit_jtr = 0.0

if ( is_jitter_rv or is_jitter_tr ) :
  for o in range(0,n_jrv):
    fit_jrv[o] = best_value(params_jitter[o],maxloglike,get_value)
  fit_jtr = best_value(params_jitter[n_jrv],maxloglike,get_value)

fit_trends = [0.0]*2
if ( is_linear_trend != 'f' or is_quadratic_trend != 'f' ):
 fit_trends = [best_value(params_trends[0],maxloglike,get_value),best_value(params_trends[1],maxloglike,get_value)]

#Calculate the final chi2 for each case
log_like_rv, chi2tot_val_rv, dummy = \
  pti.get_loglike(mega_time,mega_rv,megax,megay,mega_err,megae,\
                     tlab,jrvlab,[True,False],flags,t_cad, n_cad, \
                     fit_pars,rvs_pars,ldc_pars,fit_trends,fit_jrv,fit_jtr \
                     )

log_like_tr, dummy, chi2tot_val_tr = \
  pti.get_loglike(mega_time,mega_rv,megax,megay,mega_err,megae,\
                     tlab,jrvlab,[False,True],flags,t_cad, n_cad, \
                     fit_pars,rvs_pars,ldc_pars,fit_trends,fit_jrv,fit_jtr \
                     )

log_like_total = log_like_rv + log_like_tr

bic_from_loglikelihood = np.log(ndata)*npars - 2.0*log_like_total
aic_from_loglikelihood = 2.0*npars - 2.0*log_like_total

chi2tot_val  = chi2tot_val_rv + chi2tot_val_tr
chi2_val = chi2tot_val / ( ndata - npars )

if ( scale_error_bars ):
  s_factor = np.sqrt( chi2_val )
  if ( chi2_val > 1.0 ):
    s_factor = 1.0 / s_factor
else:
  s_factor = 1.0


if ( method == 'mcmc' or method == 'plot' ):

  base = 4 #Where do the parameters start?
#Fitted parameters
  T0_vec = [None]*nplanets
  P_vec  = [None]*nplanets
  e_vec  = [None]*nplanets
  w_vec  = [None]*nplanets
  b_vec  = [None]*nplanets
  ar_vec = [None]*nplanets
  rr_vec = [None]*nplanets
  k_vec  = [None]*nplanets
#Derived parameters
  Teq_vec= [None]*nplanets #Planet temperature
  r_vec  = [None]*nplanets #planet radius
  a_vec  = [None]*nplanets #semi-major axis
  m_vec  = [None]*nplanets #planet mass
  i_vec  = [None]*nplanets #orbit inclination
  Tpe_vec= [None]*nplanets #Periastron passage time
  ds_vec = [None]*nplanets #stellar density
  dp_vec = [None]*nplanets #planet density
  gp_vec = [None]*nplanets #planet surface gravity
  trt_vec= [None]*nplanets #Total transit duration
  tri_vec= [None]*nplanets #Ingress/egress duration


#Print the summary
  out_params_file = outdir+'/'+star+'_params.dat'
  out_tex_file = outdir+'/'+star+'_params.tex'
  opars = open(out_params_file,'w')
  otex  = open(out_tex_file,'w')
  opars.write('\n')
  opars.write ('--------------------------------------------------------------\n')
  opars.write('Summary:\n')
  opars.write('N_chains         = %8i \n'%nwalkers)
  opars.write('N_iter           = %8i \n'%nconv)
  opars.write('thin_factor      = %8i \n'%thin_factor)
  opars.write('N_rv_data        = %8i \n'%len(mega_time))
  opars.write('N_tr_data        = %8i \n'%len(megax))
  opars.write('N_data           = %8i \n'%ndata)
  opars.write('N_pars           = %8i \n'%npars)
  opars.write('chi2_rv          = %4.4f\n' %(chi2tot_val_rv))
  opars.write('chi2_tr          = %4.4f\n' %(chi2tot_val_tr))
  opars.write('chi2             = %4.4f\n' %(chi2tot_val))
  opars.write('dof              = %8i \n' %(ndata - npars))
  opars.write('chi2/dof         = %4.4f \n' %chi2_val)
  opars.write('ln likelihood_rv = %4.4f\n' %(log_like_rv))
  opars.write('ln likelihood_tr = %4.4f\n' %(log_like_tr))
  opars.write('ln likelihood    = %4.4f\n' %(log_like_total))
  opars.write('BIC              = %4.4f\n' %(bic_from_loglikelihood))
  opars.write('AIC              = %4.4f\n' %(aic_from_loglikelihood))
  opars.write ('--------------------------------------------------------------\n')
  opars.write ('             INPUT STELLAR PARAMETERS\n')
  opars.write ('--------------------------------------------------------------\n')
  opars.write ('M_*     = %4.7f - %4.7f + %4.7f solar masses\n'%(mstar_mean,mstar_sigma,mstar_sigma))
  opars.write ('R_*     = %4.7f - %4.7f + %4.7f solar radii\n'%(rstar_mean,rstar_sigma,rstar_sigma))
  opars.write ('T_*     = %4.7f - %4.7f + %4.7f K\n'%(tstar_mean,tstar_sigma,tstar_sigma))
  #tex
  #otex.write ('--------------------------------------------------------------\n')
  #otex.write ('             INPUT STELLAR PARAMETERS\n')
  #otex.write ('--------------------------------------------------------------\n')
  otex.write ('\\newcommand{\smass}[1][$M_{\odot}$]{ $ %4.7f _{- %4.7f}^{ + %4.7f} $ #1} \n'%(mstar_mean,mstar_sigma,mstar_sigma))
  otex.write ('\\newcommand{\sradius}[1][$R_{\odot}$]{ $%4.7f _{ - %4.7f}^{ + %4.7f} $ #1}\n'%(rstar_mean,rstar_sigma,rstar_sigma))
  otex.write ('\\newcommand{\stemp}[1][$\mathrm{K}$]{ $ %4.7f _{- %4.7f}^{ + %4.7f} $ #1 }\n'%(tstar_mean,tstar_sigma,tstar_sigma))

  #Print the data for all the planets
  for o in range(0,nplanets):
    T0_vec[o] = params[base + 0]
    P_vec[o]  = params[base + 1]
    e_vec[o]  = params[base + 2]
    w_vec[o]  = params[base + 3]
    b_vec[o]  = params[base + 4]
    ar_vec[o] = params[base + 5]
    rr_vec[o] = params[base + 6]
    k_vec[o]  = params[base + 7]

#STARTING CALCULATIONS

    if ( is_log_P ):
      P_vec[o] = 10.0**(P_vec[o])
    if ( is_den_a ):
      ar_vec[o] = ar_vec[o]*(P_vec[o]*P_vec[o]*7464960000.*G_cgs/3.0/np.pi)**(1./3.)
      if ( o > 0 ):
        ar_vec[o] = (P_vec[o]/P_vec[0])**(2./3.)*ar_vec[0]
    if ( is_log_k ):
      k_vec[o] = 10.0**(k_vec[o])

    if ( is_ew ):
      e_dum = list(e_vec[o])
      w_dum = list(w_vec[o])
      e_vec[o] = e_vec[o]**2 + w_vec[o]**2
      w_vec[o] = np.arctan2(e_dum,w_vec[o])
      w_vec[o] = w_vec[o] % (2*np.pi)

  #Change between b and i
    if ( is_b_factor ):
      i_vec[o] = list(b_vec[o])
      i_vec[o] = np.arccos( b_vec[o] / ar_vec[o] * \
              ( 1.0 + e_vec[o] * np.sin(w_vec[o] ) / ( 1.0 - e_vec[o]**2 ) ) )
    else:
      #calculate the impact parameter (eq. 7 Winn 2014)
      #wo is the star periastron, add pi to have the planet one
      i_vec[o] = list(b_vec[o])
      b_vec[o] =  ar_vec[o] * np.cos(b_vec[o]) * ( ( 1. - e_vec[o]**2 ) \
               / ( 1.0 + e_vec[o]*np.sin(w_vec[o] )))
      i_vec[o] = np.array(i_vec[o])

    #Calculate equilibrium temperature
    #assuming albedo=0
    Teq_vec[o] = get_teq(tstar,0.0,1.0,ar_vec[o])

    #Get the star periastron pasage
    w_s_deg, w_s_deg_l, w_s_deg_r = find_vals_perc(w_vec[o]*180./np.pi,s_factor)
    #planet periastron passage
    w_p_deg = (w_s_deg + 180.) % 360

  #Transit durations aproximations (eq. 14, 15, 16 from Winn 2014)
    ec_factor = np.sqrt(( 1. - e_vec[o]*e_vec[o] )) / ( 1.0 + e_vec[o]*np.sin(w_vec[o] ))
    trt_vec[o] = np.sqrt( (1. + rr_vec[o])**2 - b_vec[o]**2 ) / ( ar_vec[o] * np.sin(i_vec[o]))
    trt_vec[o] = P_vec[o] / np.pi * np.arcsin(trt_vec[o]) * ec_factor * 24.0
    tri_vec[o] = np.sqrt( (1. - rr_vec[o])**2 - b_vec[o]**2 ) / ( ar_vec[o] * np.sin(i_vec[o]))
    tri_vec[o] = P_vec[o] / np.pi * np.arcsin(tri_vec[o]) * ec_factor * 24.0
    #tri_vec[o] = ( trt_vec[o] - tri_vec[o] ) / 2.0 #ingress egress time
    #Calculate the star density from transit data
    #Eq. (30) Winn 2014
    ds_vec[o] = get_rhostar(P_vec[o],ar_vec[o]) #cgs

    #Time of periastron passage
    Tpe_vec[o] = list(T0_vec[o])
    for m in range(0,len(Tpe_vec[o])):
      Tpe_vec[o][m] = pti.find_tp(T0_vec[o][m],e_vec[o][m],w_vec[o][m],P_vec[o][m])

    #Density from the input stellar parameters
    irho_vec = mstar/rstar**3 * 1.411

    #Get planet mass, radius and orbit semi-major axis in real units
    r_vec[o] = rr_vec[o] * rstar
    a_vec[o] = ar_vec[o] * rstar * S_radius_SI / AU_SI
    m_vec[o] = planet_mass(mstar,k_vec[o]*1.e3,P_vec[o],e_vec[o],i_vec[o])

    #Planet-star distance at the time of eclipse
    true_anomaly_vec = [None]*len(T0_vec[o])
    for l in range(0,len(true_anomaly_vec)):
      true_anomaly_vec[l] = pti.find_anomaly(T0_vec[o][l],T0_vec[o][l],e_vec[o][l],w_vec[o][l],P_vec[o][l])
    dummy_anomaly = np.concatenate(true_anomaly_vec)
    psd_vec = ar_vec[o] * ( 1. - e_vec[o]**2) / ( 1. + e_vec[o]*np.cos(dummy_anomaly))
    psd_vec_units = psd_vec * rstar * S_radius_SI / AU_SI

    #Kepler cociente
    pa_vec = (P_vec[o]*3600.*24.0)**2 * S_GM_SI * (mstar + m_vec[o])
    pa_vec = pa_vec / ( 4.*np.pi**2 * (a_vec[o]*AU_SI)**3 )

    #stimate planet gravity and density
    pden_vec = m_vec[o] / r_vec[o]**3 #solar units
    pden_vec = pden_vec * S_den_cgs   #g/cm^3
    #We can stimate planet surface gravity (eq. (31) Winn)
    pgra_vec = (P_vec[o]*24.*3600.) * (rr_vec[o]/ar_vec[o])**2 * np.sin(i_vec[o])
    pgra_vec = 2. * np.pi * np.sqrt(1. - e_vec[o]**2) * (k_vec[o]*1.e5) / pgra_vec #cm/s^2

    #Estimate surface gravity from the derived parameters
    pgra_vec2 = m_vec[o] / r_vec[o]**2   #in solar units
    pgra_vec2 = pgra_vec2 * 28.02 * 981. #cm/s^2

    #Stellar luminosity in solar units
    Ls = (rstar)**2*(tstar/S_Teff)**4
    #planet insolation in Flux received at Earth
    Fp = Ls/a_vec[o]**2

    #Estimate the stellar mass assuming the surface gravity of the planet is true
    #Planet mass from surface gravity of the planet
    mpgra = pgra_vec*(r_vec[o]*S_radius_cgs)**2/G_cgs #g
    #Stellar mass from surface gravity of the planet
    msgra = (k_vec[o]*1.e5)*np.sqrt(1. - e_vec[o]**2)/np.sin(i_vec[o])
    msgra = msgra * ( (P_vec[o]*24.*3600.) / 2. / np.pi / G_cgs )**(1./3.)
    msgra = (mpgra / msgra)**(3./2.)
    msgra = msgra - mpgra
    msgra = msgra / S_GM_cgs * G_cgs

    #Convert units
    usymbol = '{\odot}'
    if ( unit_mass == 'earth'):
      usymbol = '{\oplus}'
      if ( fit_rv ):
        m_vec[o] = m_vec[o] * S_GM_SI / E_GM_SI
      if ( fit_tr ):
        r_vec[o] = r_vec[o] * S_radius_SI / E_radius_e_SI
    elif ( unit_mass == 'jupiter'):
      usymbol = '\mathrm{J}'
      if ( fit_rv ):
        m_vec[o] = m_vec[o] * S_GM_SI / J_GM_SI
      if ( fit_tr ):
        r_vec[o] = r_vec[o] * S_radius_SI / J_radius_e_SI

    #Print the parameters
    #Fitted parameters
    opars.write ('--------------------------------------------------------------\n')
    opars.write ('                   Parameters %s\n' %( star +' '+ plabels[o]))
    opars.write ('-------------------------Fitted-------------------------------\n')
    pl = plabels[o]
    print_values(T0_vec[o],'T0','Tzero'+pl,'days','days')
    print_values(P_vec[o],'P','P'+pl,'days','days')
    if ( is_ew ):
      print_values(e_dum,'ew 1','esin'+pl,' ',' ')
      print_values(w_dum,'ew 2','ecos'+pl,' ',' ')
    else:
      print_values(e_vec[o],'e','e'+pl,' ',' ')
      print_values(w_vec[o]*180./np.pi,'w','w'+pl,'deg','deg')
    if ( fit_tr[o] ):
      print_values(b_vec[o],'b','b'+pl,' ',' ')
      if ( is_den_a ):
        print_values(params[4+5],'rho*^1/3','dentrhee'+pl,'g^{1/3}/cm','${\\rm g^{1/3}\,cm^{-1}}$')
      else:
        print_values(ar_vec[o],'a/R*','ar'+pl,' ',' ')
      print_values(rr_vec[o],'rp/R*','rr'+pl,' ',' ')
    if ( fit_rv[o] ):
      print_values(k_vec[o]*1e3,'K','k'+pl,'m/s','${\\rm m\,s^{-1}}$')
    opars.write ('-------------------------Derived------------------------------\n')
    if (fit_rv[o]): print_values(m_vec[o],'Mp','mp'+pl,'M_'+unit_mass,'$M_'+usymbol+'$')
    if (fit_tr[o]): print_values(r_vec[o],'Rp','rp'+pl,'R_'+unit_mass,'$R_'+usymbol+'$')
    print_values(Tpe_vec[o],'Tperi','Tperi'+pl,'days','days')
    if ( is_ew ):
      print_values(e_vec[o],'e','e'+pl,' ',' ')
      print_values(w_vec[o]*180./np.pi,'w','w'+pl,'deg','deg')
    if ( fit_tr[o]):
      print_values(i_vec[o]*180./np.pi,'i','i'+pl,'deg','deg')
      if ( is_den_a ):
        print_values(ar_vec[o],'a/R*','ar'+pl,' ',' ')
      print_values(a_vec[o],'a','a'+pl,'AU','AU')
      print_values(Fp,'Insolation','insolation'+pl,'F_Earth','${\\rm F_{\\oplus}}$')
      print_values(ds_vec[o],'rho*','denstr'+pl,'g/cm^3 (transit)','${\\rm g\,cm^{-3}}$')
      print_values(irho_vec,'rho*','denssp'+pl,'g/cm^3 (stellar paramters)','${\\rm g\,cm^{-3}}$')
      print_values(pden_vec,'rho_p','denp'+pl,'g/cm^3','${\\rm g\,cm^{-3}}$')
      print_values(pgra_vec,'g_p','grap'+pl,'cm/s^2 (K and Rp/R*)','${\\rm cm\,s^{-2}}$')
      print_values(pgra_vec2,'g_p','grappars'+pl,'cm/s^2 (planet parameters)','${\\rm cm\,s^{-2}}$')
      print_values(msgra,'M_*','mspars'+pl,' solar masses (scaled parameters)','$M_\odot$')
      print_values(Teq_vec[o],'Teq','Teq'+pl,'K (albedo=0)','K')
      print_values(trt_vec[o],'T_tot','ttot'+pl,'hours','hours')
      print_values(tri_vec[o],'T_full','tful'+pl,'hours','hours')
    opars.write ('--------------------------------------------------------------\n')

    #Let us change to the next planet
    base = base + 8

#The other parameters
q1_vec = params[base]
q2_vec = params[base+1]

u1_vec = np.sqrt(q1_vec)
u2_vec = u1_vec * (1. - 2.*q2_vec)
u1_vec = 2.*u1_vec*q2_vec

rv_vec = [None]*nt
for o in range(0,nt):
  rv_vec[o] = params[base+2+o]

opars.write ('--------------------  Other parameters -----------------------\n')
if ( total_tr_fit ):
  print_values(q1_vec,'q1','qone','','')
  print_values(q2_vec,'q2','qtwo','','')
  print_values(u1_vec,'u1','uone','','')
  print_values(u2_vec,'u2','utwo','','')
if ( total_rv_fit ):
  for o in range(0,nt):
    print_values(rv_vec[o],'Sys. vel. '+telescopes_labels[o],telescopes_labels[o],'m/s','${\\rm m\,s^{-1}}$')
opars.write ('--------------------------------------------------------------\n')

if ( is_jitter_rv or is_jitter_tr ):
  if ( total_rv_fit ):
    for o in range(0,n_jrv):
      print_values(params_jitter[o]*1.e3,telescopes_labels[o]+' jitter','j'+telescopes_labels[o],'m/s','${\\rm m\,s^{-1}}$')
  if ( total_tr_fit ):
    print_values(params_jitter[n_jrv],'tr jitter','jtr','','')
  opars.write ('--------------------------------------------------------------\n')
if ( is_linear_trend != 'f' or is_quadratic_trend != 'f' ):
  if ( total_rv_fit ):
    print_values(params_trends[0]*1.e3,'linear trend','ltrend','m/s/days','${\\rm m\,s^{-1}\,d^{-1}}$')
    print_values(params_trends[1]*1.e3,'quadratic trend','qtrend','m/s/days^2','${\\rm m\,s^{-1}\,d^{-2}}$')
  opars.write ('--------------------------------------------------------------\n')
opars.write('\n')

#RESIZE TRANSIT ERROR BARS
if ( is_jitter_tr and resize_tr ):
  jit_tr = best_value(params_jitter[n_jrv],maxloglike,get_value)
  for o in range(0,len(et)):
      for m in range(0,len(et[o])):
              et[o][m] = np.sqrt(et[o][m]**2 + jit_tr**2)
  for o in range(0,len(megae)):
    megae[o] = np.sqrt( megae[o]**2 + jit_tr**2)

if ( total_rv_fit ):
  new_errs_all = [None]*len(errs_all)
  for o in range(0,len(errs_all)):
    new_errs_all[o] = list(errs_all[o])

if ( is_jitter_rv and resize_rv ):
    for j in range(0,n_jrv):
      jit_rv = best_value(params_jitter[j],maxloglike,get_value)
      for o in range(0,len(errs_all[j])):
          new_errs_all[j][o] = 1.e3*np.sqrt(errs_all[j][o]**2 + jit_rv**2)

opars.close()
otex.close()

#Print the output in the screen
dummy_file = open(out_params_file)
for line in dummy_file:
  print line,
dummy_file.close()
