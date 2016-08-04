from matplotlib import gridspec

#what is the minimum chi2
minchi2_index = np.argmin(chi2)

#Create the stellar data
mstar = np.random.normal(loc=mstar_mean,scale=mstar_sigma,size=new_nwalkers*nconv)
rstar = np.random.normal(loc=rstar_mean,scale=rstar_sigma,size=new_nwalkers*nconv)
tstar = np.random.normal(loc=tstar_mean,scale=tstar_sigma,size=new_nwalkers*nconv)
inclination = np.random.normal(loc=inclination_mean,scale=inclination_sigma,size=new_nwalkers*nconv)

if ( nplanets == 1 ):

  #If some parameters are transformed let us go back

  if (is_log_P):
    Po = np.power(10.,Po)

  if (is_ew):
    dummy_e = eo
    eo = eo * eo + wo * wo
    wo = np.arctan2(dummy_e,wo)

  if ( fit_tr ):

    if (is_log_a):
      ao = np.power(10.0,ao)
    if (is_sini):
      io = np.arccos(io/ao)
    inclination = io


    #Calculate equilibrium temperature
    #assuming albedo=0
    Teqo = get_teq(tstar,0.0,1.0,ao)

    #Take back the u1 and u2 values, Kipping 2013
    u1o = 2*np.sqrt(q1o)*q2o
    u2o = np.sqrt(q1o)*(1.-2.*q2o) 
 
    #Calculate the planet size in solar radii
    #rstar must be given in solar radius           
    rpo = pzo * rstar
    
    #calculate the impact parameter (eq. 7 Winn 2014)
    bo =  ao * np.cos(io) * ( ( 1. - eo*eo ) / ( 1.0 + eo*np.sin(wo)))
    #Transit durations aproximations (eq. 14, 15, 16 from Winn 2014)
    ec_factor = np.sqrt(( 1. - eo*eo )) / ( 1.0 + eo*np.sin(wo))
    tto = np.sqrt( (1. + pzo)**2 - bo**2 ) / ( ao * np.sin(io))
    tto = Po / np.pi * np.arcsin(tto) * ec_factor * 24.0
    tfo = np.sqrt( (1. - pzo)**2 - bo**2 ) / ( ao * np.sin(io))
    tfo = Po / np.pi * np.arcsin(tfo) * ec_factor * 24.0
    tfo = ( tto - tfo ) / 2.0 #ingress egress time
    #Calculate the star density from transit data
    #Eq. (30) Winn 2014
    rhoo = get_rhostar(Po,ao) #cgs
    #physical semi-mayor axis
    aphyo = ao * rstar * S_radius_SI / AU_SI

  if ( fit_rv ):

    if ( is_log_k ):
      ko = np.power(10.,ko)

    for j in range(0,nt):
      if ( is_log_rv0 ):
        vo[j] = np.power(10.,vo[j])

    #Get tp and mass from the the values
    tpo = [None]*len(t0o) 
    masso = [None]*len(t0o) 
    for m in range(0,len(t0o)):
      tpo[m] = pti.find_tp(t0o[m],eo[m],wo[m],Po[m])

    masso = planet_mass(mstar,ko*1.e3,Po,eo,inclination)

    if ( fit_tr ):               
      #We can estimate the planet density
      rho_p = masso / ( rpo**3 ) #all is given in solar units
      rho_p = rho_p * S_den_cgs  #sun density [g/cm^3]
      #We can stimate planet surface gravity (eq. (31) Winn)
      gpo = (Po*24.*3600.) * (pzo/ao)**2 * np.sin(io)
      gpo = 2.*np.pi * np.sqrt(1. - eo*eo) * (ko*1.e5) / gpo # cm/s^2

  #Change to earth or jupiter parameters
  #based on IAU resolution http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1605.09788

  usymbol = '{\odot}'
  if ( unit_mass == 'earth'):
    usymbol = '{\oplus}'
    if ( fit_rv ):
      masso = masso * S_GM_SI / E_GM_SI
    if ( fit_tr ):
      rpo = rpo * S_radius_SI / E_radius_e_SI
  elif ( unit_mass == 'jupiter'):
    usymbol = '\mathrm{J}'
    if ( fit_rv ):
      masso = masso * S_GM_SI / J_GM_SI
    if ( fit_tr ):
      rpo = rpo * S_radius_SI / J_radius_e_SI
	
	
  #Calculate the BIC
  if (fit_rv and fit_tr ):
    ndata = len(megax) + len(mega_rv)
    npars = sum(what_fit) + nt - 1
  elif(fit_rv and not fit_tr):
    ndata = len(mega_rv)
    npars = sum(what_fit) + nt - 1
  elif(not fit_rv and fit_tr):
    ndata = len(megax)
    npars = sum(what_fit)
			

  #Let us estimate the errors and print them
  if ( errores == 'perc' ):

    #Obtain the minimin chi2 and calculate the reduced chi2
    chi2tot_val  = np.amin(chi2)
    chi2_val = chi2tot_val / ( ndata - npars )

    if ( scale_error_bars ):
      s_factor = np.sqrt( chi2_val )
      if ( chi2_val > 1.0 ):
        s_factor = 1.0 / s_factor
    else:
      s_factor = 1.0

    #Input stellar parameters
    if (mstar.__class__ == float ):
      if ( fit_rv ):
        ms_val = mstar
        ms_errr = 0.0
        ms_errl = 0.0
      if ( fit_tr ):
        rs_val = rstar
        rs_errr = 0.0
        rs_errl = 0.0
    else:
      if ( fit_rv ):
        ms_val, ms_errl, ms_errr = find_vals_perc(mstar,1.0)
      if ( fit_tr ):
        rs_val, rs_errl, rs_errr = find_vals_perc(rstar,1.0)
        ts_val, ts_errl, ts_errr = find_vals_perc(tstar,1.0)


    t0_val, t0_errl, t0_errr = find_vals_perc(t0o,s_factor)
    T0 = t0_val
    P_val, P_errl, P_errr = find_vals_perc(Po,s_factor)
    e_val,e_errl, e_errr = find_vals_perc(eo,s_factor)
    w_val,w_errl, w_errr = find_vals_perc(wo,s_factor)
    if (w_val < 0.0 ):
      w_val = w_val + 2.e0 * np.pi	
    #Transform periastron from rad to degrees
    w_deg = w_val * 180. / np.pi
    w_deg_errl = w_errl * 180. / np.pi
    w_deg_errr = w_errr * 180. / np.pi

    if ( fit_tr ):
      i_val,i_errl, i_errr = find_vals_perc(io,s_factor)
      a_val,a_errl, a_errr = find_vals_perc(ao,s_factor)
      aphy_val,aphy_errl, aphy_errr = find_vals_perc(aphyo,s_factor)
      q1_val,q1_errl, q1_errr = find_vals_perc(q1o,s_factor)
      q2_val,q2_errl, q2_errr = find_vals_perc(q2o,s_factor)
      u1_val,u1_errl, u1_errr = find_vals_perc(u1o,s_factor)
      u2_val,u2_errl, u2_errr = find_vals_perc(u2o,s_factor)
      pz_val,pz_errl, pz_errr = find_vals_perc(pzo,s_factor)
      rp_val,rp_errl, rp_errr = find_vals_perc(rpo,s_factor)
      b_val , b_errl, b_errr  = find_vals_perc(bo,s_factor)
      tt_val , tt_errl, tt_errr  = find_vals_perc(tto,s_factor)
      tf_val , tf_errl, tf_errr  = find_vals_perc(tfo,s_factor)
      rho_val , rho_errl, rho_errr  = find_vals_perc(rhoo,s_factor)
      Teq_val,Teq_errl, Teq_errr = find_vals_perc(Teqo,s_factor)
      if ( fit_rv ):
        rhop_val,rhop_errl, rhop_errr = find_vals_perc(rho_p,s_factor)
        gp_val,gp_errl, gp_errr = find_vals_perc(gpo,s_factor)
      i_deg = i_val * 180. / np.pi
      i_deg_errl = i_errl * 180. / np.pi
      i_deg_errr = i_errr * 180. / np.pi
    if ( fit_rv ):
      k_val, k_errl, k_errr = find_vals_perc(ko,s_factor)
      alpha_val, alpha_errl, alpha_errr = find_vals_perc(alphao,s_factor)
      beta_val, beta_errl, beta_errr = find_vals_perc(betao,s_factor)
      m_val, m_errl, m_errr  	= find_vals_perc(masso,s_factor)
      tp_val, tp_errl, tp_errr= find_vals_perc(tpo,s_factor)
      #Systemic velocities are stored in v_val
      v_val = [None]*nt
      v_errl = [None]*nt
      v_errr = [None]*nt
      for j in range(0,nt):
        v_val[j], v_errl[j], v_errr[j] = find_vals_perc(vo[j],s_factor)
      if ( not fit_tr ):
        inc_deg = inclination * 180 / np.pi
        if ( inclination.__class__ == float ):
          iinp_val = inc_deg
          iinp_errl = 0.0
          iinp_errr = 0.0 
        else:
          iinp_val,iinp_errl, iinp_errr = find_vals_perc(inc_deg,s_factor)
	

    bic2 = get_BIC(chi2tot_val)

    #Print the best fit values values
    print ('')
    print 'Summary:'
    print 'N_chains    = ', nwalkers
    print 'N_conv      = ', nconv
    print 'thin_factor = ', thin_factor
    print 'N_data      = ', ndata
    print 'N_pars      = ', npars
    print ('chi2       = %1.4f' %(chi2tot_val))
    print 'DOF         = ', ndata - npars
    print ('chi2_red   = %1.4f' %(chi2_val))
    print 'scale factor= ', s_factor
    print ('BIC        = %1.4f' %(bic2))
    print 'Input parameters'
    if ( fit_rv ):
      print ('M_*     = %4.7f + %4.7f - %4.7f solar masses'%(ms_val,ms_errr,ms_errl))
    if ( fit_tr ):
      print ('R_*     = %4.7f + %4.7f - %4.7f solar radii'%(rs_val, rs_errr , rs_errl))
      print ('T_*     = %4.7f + %4.7f - %4.7f solar radii'%(ts_val, ts_errr , ts_errl))
    if ( fit_rv and not fit_tr ):
      print ('i       = %4.7f + %4.7f - %4.7f deg'%(iinp_val, iinp_errr , iinp_errl))
    print ('')
    print ('The best fit planet parameters are:')
    print ('T0    = %4.7f + %4.7f - %4.7f days'%(t0_val,t0_errr,t0_errl))
    print ('P     = %4.7f + %4.7f - %4.7f days'%(P_val, P_errr , P_errl))
    print ('e     = %4.4f + %4.4f - %4.4f     '%(e_val, e_errr , e_errl))
    print ('w     = %4.4f + %4.4f - %4.4f deg '%(w_deg,w_deg_errr, w_deg_errl))
    if (fit_tr):
      print ('Transit fit parameters:')
      print ('i     = %4.4f + %4.4f - %4.4f deg' %(i_deg,i_deg_errr, i_deg_errl))
      print ('a/r*  = %4.4f + %4.4f - %4.4f    '%(a_val, a_errr , a_errl))
      print ('rp/r* = %4.4f + %4.4f - %4.4f    '%(pz_val,pz_errr, pz_errl))
      print ('q_1    = %4.4f + %4.4f - %4.4f    '%(q1_val,q1_errr, q1_errl))
      print ('q_2    = %4.4f + %4.4f - %4.4f    '%(q2_val,q2_errr, q2_errl))
    if (fit_rv):
      print ('RV fit parameters:')
      print ('alpha = %4.4e + %4.4e - %4.4e     '%(alpha_val, alpha_errr , alpha_errl))
      print ('beta  = %4.4e + %4.4e - %4.4e     '%(beta_val, beta_errr , beta_errl))
      print ('K     = %4.4f + %4.4f - %4.4f m/s'%(k_val/1.e-3,(k_errr)/1.e-3, (k_errl)/1.e-3))
      for i in range(0,nt):
        print ('%s v0  = %4.4f + %4.4f - %4.4f km/s'%(telescopes[i], \
              v_val[i],v_errr[i],v_errl[i]))
      print ('')
    
    print ('Derived parameters:')
    if (fit_tr):
      print ('r_p   = %4.4f + %4.4f - %4.4f R_%s' 	        %(rp_val,rp_errr, rp_errl,unit_mass))
      print ('a   = %4.4f + %4.4f - %4.4f  AU ' 		%(aphy_val, aphy_errr , aphy_errl))  
      print ('b r*  = %4.4f + %4.4f - %4.4f' 	         	%(b_val,b_errr, b_errl))
      print ('t_total = %4.4f + %4.4f - %4.4f hours' 		%(tt_val,tt_errr, tt_errl))
      print ('t_in/eg = %4.4f + %4.4f - %4.4f hours' 		%(tf_val,tf_errr, tf_errl))
      print ('rho_* = %4.4f + %4.4f - %4.4f g/cm^3' 		%(rho_val,rho_errr, rho_errl))
      print ('u_1    = %4.4f + %4.4f - %4.4f    ' 		%(u1_val,u1_errr, u1_errl))
      print ('u_2    = %4.4f + %4.4f - %4.4f    ' 		%(u2_val,u2_errr, u2_errl))
      print ('T_eq   = %4.4f + %4.4f - %4.4fK   '               %(Teq_val,Teq_errr, Teq_errl))
    if (fit_rv):
      print ('Tp    = %4.4f + %4.4f - %4.4f days' 		%(tp_val,tp_errr, tp_errl))
      print ('mp    = %4.4f + %4.4f - %4.4f %s masses' 		%(m_val,m_errr, m_errl,unit_mass))
      if ( fit_tr ):
        print ('rho_p = %4.4f + %4.4f - %4.4f g/cm^3' 		%(rhop_val,rhop_errr, rhop_errl))
        print ('g_p = %4.4f + %4.4f - %4.4f cm/s^2' 		%(gp_val,gp_errr, gp_errl))
    
    if ( latex_values ):
      print '%LaTeX commands of the parameters'
      if ( fit_rv ):
        print ('\\newcommand{\smass}[1][$ %s $]{$ %4.7f ^{+ %4.7f}_{ - %4.7f}$ #1}'%('M_{\odot}',ms_val,ms_errr,ms_errl))
      if ( fit_tr ):
        print ('\\newcommand{\sradius}[1][$ %s $]{$ %4.7f ^{+ %4.7f} _{ - %4.7f} $ #1}'%('R_{\odot}',rs_val, rs_errr , rs_errl))
        print ('\\newcommand{\stemp}[1][K]{$ %4.2f ^{+ %4.2f} _{ - %4.2f} $ #1}'%(ts_val, ts_errr , ts_errl))
      if ( fit_rv and not fit_tr ):
        print ('\\newcommand{\inclination}[1][deg]{$%4.7f^{+%4.7f }_{- %4.7f} $ #1}'%(iinp_val, iinp_errr , iinp_errl))
      print ('%')
      print ('%The best fit planet parameters are:')
      print ('\\newcommand{\\tzero}[1][days]{ $%4.7f^{ + %4.7f}_{ - %4.7f}$ #1}'%(t0_val,t0_errr,t0_errl))
      print ('\\newcommand{\porb}[1][days]{$%4.7f^{ + %4.7f}_{ - %4.7f}$ #1}'%(P_val, P_errr , P_errl))
      print ('\\newcommand{\ec}[1][]{ $%4.4f^{ + %4.4f}_{ - %4.4f}$ #1 }'%(e_val, e_errr , e_errl))
      print ('\\newcommand{\periastron}[1][deg]{$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1} '%(w_deg,w_deg_errr, w_deg_errl))
      if (fit_tr):
        print ('%Transit fit parameters:')
        print ('\\newcommand{\inclination}[1][deg]{$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1}' %(i_deg,i_deg_errr, i_deg_errl))
        print ('\\newcommand{\saxis}[1][]{$ %4.4f^{+ %4.4f}_{ - %4.4f}$ #1} '%(a_val, a_errr , a_errl))
        print ('\\newcommand{\spradius}[1][]{$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1} '%(pz_val,pz_errr, pz_errl))
        print ('\\newcommand{\qone}[1][]{$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1}    '%(q1_val,q1_errr, q1_errl))
        print ('\\newcommand{\qtwo}[1][]{$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1}    '%(q2_val,q2_errr, q2_errl))
      if (fit_rv):
        print ('%RV fit parameters:')
        print ('\\newcommand{\krv}[1][m\, s$^{-2}$]{ $ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1}'%(k_val/1.e-3,(k_errr)/1.e-3, (k_errl)/1.e-3))
        print ('\\newcommand{\\alpha}[1][]{%4.4e^{ + %4.4e }_{- %4.4e} #1}    '%(alpha_val, alpha_errr , alpha_errl))
        print ('\\newcommand{\\beta}[1][]{%4.4e^{ + %4.4e }_{- %4.4e} #1}     '%(beta_val, beta_errr , beta_errl))
        for i in range(0,nt):
          print ('\\newcommand{\\vel%s}[1][km\, s$^{-2}$]{$ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1}'%(telescopes[i], \
                v_val[i],v_errr[i],v_errl[i]))
      print ('')
    
      print ('%Derived parameters:')
      if (fit_tr):
        if (unit_mass == 'earth'):
          print ('\\newcommand{\pradius}[1][$R_{%s}$]{$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1}' 	        %(usymbol,rp_val,rp_errr, rp_errl))
        print ('\\newcommand{\\axis}[1][AU]{$ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1} ' 		%(aphy_val, aphy_errr , aphy_errl))  
        print ('\\newcommand{\impactp}[1][]{$ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1} ' 	         	%(b_val,b_errr, b_errl))
        print ('\\newcommand{\\ttotal}[1][hours] {$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1 }' 		%(tt_val,tt_errr, tt_errl))
        print ('\\newcommand{\\tineg }[1][hours] {$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1 }' 		%(tf_val,tf_errr, tf_errl))
        print ('\\newcommand{\sden }[1][g\, cm$^{-3}$]{$ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1}' 		%(rho_val,rho_errr, rho_errl))
        print ('\\newcommand{\uone}[1][]{ $ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1 } ' 		%(u1_val,u1_errr, u1_errl))
        print ('\\newcommand{\utwo}[1][]{ $ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1 }  ' 		%(u2_val,u2_errr, u2_errl))
        print ('\\newcommand{\\tequi}[1][K]{ $ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1}  ' 		%(Teq_val,Teq_errr, Teq_errl))
      if (fit_rv):
        print ('\\newcommand{\\tp}[1][days]{$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1 }' 		%(tp_val,tp_errr, tp_errl))
        print ('\\newcommand{\pmass}[1][$M_{%s}$]{ $ %4.4f^{ + %4.4f}_{ - %4.4f} $ #1}' 		%(usymbol,m_val,m_errr, m_errl))
        if ( fit_tr ):
          print ('\\newcommand{\pden}[1][g\, cm$^{-3}$]{$ %4.4f^{ + %4.4f}_{ - %4.4f}$ #1}' 		%(rhop_val,rhop_errr, rhop_errl))
          print ('\\newcommand{\pg}[1][cm\, s$^-2$]{$  %4.4f^{ + %4.4f}_{ - %4.4f} $ #1}' 		%(gp_val,gp_errr, gp_errl))




#Multiplanet fit
else:

  #Define global variables
  mass_val = [None]*nplanets
  mass_errr = [None]*nplanets
  mass_errl = [None]*nplanets
  chi2_val = [None]*nplanets
  chi2_errr = [None]*nplanets
  chi2_errl = [None]*nplanets
  chi2tot_val = [None]*nplanets
  chi2tot_errr = [None]*nplanets
  chi2tot_errl = [None]*nplanets
  t0_val = [None]*nplanets
  t0_errr = [None]*nplanets
  t0_errl = [None]*nplanets
  tp_val = [None]*nplanets
  tp_errr = [None]*nplanets
  tp_errl = [None]*nplanets
  P_val = [None]*nplanets
  P_errr = [None]*nplanets
  P_errl = [None]*nplanets
  e_val = [None]*nplanets
  e_errr = [None]*nplanets
  e_errl = [None]*nplanets
  w_val = [None]*nplanets
  w_errr = [None]*nplanets
  w_errl = [None]*nplanets
  w_deg = [None]*nplanets
  w_deg_errr = [None]*nplanets
  w_deg_errl = [None]*nplanets
  k_val = [None]*nplanets
  k_errr = [None]*nplanets
  k_errl = [None]*nplanets
  alpha_val = [None]*nplanets
  alpha_errr = [None]*nplanets
  alpha_errl = [None]*nplanets
  beta_val = [None]*nplanets
  beta_errr = [None]*nplanets
  beta_errl = [None]*nplanets
  v_val = [None]*nt
  v_errr = [None]*nt
  v_errl = [None]*nt
  #end global variables

  #This function works for more than one planet
  def print_errors_planets():

    ndata = len(mega_rv)
    npars = sum(what_fit) + nt - nplanets
    npln = npars * np.log(ndata)	
	
    for j in range(0,nt):
      if ( is_log_rv0 ):
        vo[j] = np.power(10.,vo[j])

    for l in range(0,nplanets):

      if (is_log_P):
        Po[l] = np.power(10.,Po[l])

      if (is_ew):
        dummy_e = eo[l]
        eo[l] = eo[l] * eo[l] + wo[l] * wo[l]
        wo[l] = np.arctan2(dummy_e,wo[l])
	
      if ( fit_rv ):
        if ( is_log_k ):
          ko[l] = np.power(10.,ko[l])

      #Get tp and mass from the the values
      tpo = [None]*len(t0o[l]) 
      masso = [None]*len(t0o[l]) 
      for m in range(0,len(t0o[l])):
        tpo[m] = pti.find_tp(t0o[l][m],eo[l][m],wo[l][m],Po[l][m])

      #Change to earth or jupiter parameters
      #based on IAU resolution http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1605.09788

      masso = planet_mass(mstar,ko[l]*1.e3,Po[l],eo[l],inclination)

      if ( unit_mass == 'earth'):
        if ( fit_rv ):
          masso = masso * S_GM_SI / E_GM_SI
        if ( fit_tr ):
          rpo = rpo * S_radius_SI / E_radius_e_SI
      elif ( unit_mass == 'jupiter'):
        if ( fit_rv ):
          masso = masso * S_GM_SI / J_GM_SI
        if ( fit_tr ):
          rpo = rpo * S_radius_SI / J_radius_e_SI
	
      #Percentile errors
		
      if ( errores == 'perc' ):	

        s_factor = 1.0
	
        chi2tot_val[l] = np.amin(chi2[l])
        chi2_val = chi2tot_val[l] / ( ndata - npars )

        if ( scale_error_bars ):
          s_factor = np.sqrt( chi2_val )
          if ( chi2_val > 1.0 ):
            s_factor = 1.0 / s_factor
        else:
          s_factor = 1.0

        #Input stellar parameters
        if (mstar.__class__ == float ):
          if ( fit_rv ):
            ms_val = mstar
            ms_errr = 0.0
            ms_errl = 0.0
          if ( fit_tr ):
            rs_val = rstar
            rs_errr = 0.0
            rs_errl = 0.0
        else:
          if ( fit_rv ):
            ms_val, ms_errl, ms_errr = find_vals_perc(mstar,1.0)
          if ( fit_tr ):
            rs_val, rs_errl, rs_errr = find_vals_perc(rstar,1.0)

	
        t0_val[l], t0_errl[l], t0_errr[l] = find_vals_perc(t0o[l],s_factor)
        mass_val[l], mass_errl[l], mass_errr[l] = find_vals_perc(masso,s_factor)
        tp_val[l], tp_errl[l], tp_errr[l] = find_vals_perc(tpo,s_factor)
        P_val[l], P_errl[l], P_errr[l]  = find_vals_perc(Po[l],s_factor)
        e_val[l],e_errl[l], e_errr[l] 	= find_vals_perc(eo[l],s_factor)
        w_val[l],w_errl[l], w_errr[l] 	= find_vals_perc(wo[l],s_factor)
        if (w_val[l] < 0.0 ):
          w_val[l] = w_val[l] + 2 * np.pi	
        w_deg[l] = w_val[l] * 180. / np.pi
        w_deg_errl[l] = w_errl[l] * 180. / np.pi
        w_deg_errr[l] = w_errr[l] * 180. / np.pi
        if ( fit_rv ):
          k_val[l], k_errl[l], k_errr[l]  = find_vals_perc(ko[l],s_factor)
          alpha_val[l], alpha_errl[l], alpha_errr[l]  = find_vals_perc(alphao[l],s_factor)
          beta_val[l], beta_val[l], beta_val[l]  = find_vals_perc(betao[l],s_factor)
          for j in range(0,nt):
            v_val[j], v_errl[j], v_errr[j] = find_vals_perc(vo[j],s_factor)
            inc_deg = inclination * 180 / np.pi
          if ( inclination.__class__ == float ):
            iinp_val = inc_deg
            iinp_errl = 0.0
            iinp_errr = 0.0 
          else:
            iinp_val,iinp_errl,iinp_errr = find_vals_perc(inc_deg,1.0)
	
        bic2 = get_BIC(chi2tot_val[l])

        #Print the best fit values values
        if (l == 0):
          print 'N_chains    = ', nwalkers
          print 'N_conv      = ', nconv
          print 'thin_factor = ', thin_factor
          print 'N_data      = ', ndata
          print 'N_pars      = ', npars
          print ('chi2       = %1.4f' %(chi2tot_val[l]))
          print 'DOF         = ', ndata - npars
          print ('chi2_red   = %1.4f' %(chi2_val))
          print 'scale factor= ', s_factor
          print ('BIC        = %1.4f' %(bic2))
          print 'Input parameters'
          if ( fit_rv ):
            print ('M_*     = %4.7f + %4.7f - %4.7f solar masses'%(ms_val,ms_errr,ms_errl))
          if ( fit_tr ):
            print ('R_*     = %4.7f + %4.7f - %4.7f solar radii'%(rs_val, rs_errr , rs_errl))
          if ( fit_rv and not fit_tr ):
            print ('i       = %4.7f + %4.7f - %4.7f deg'%(iinp_val, iinp_errr , iinp_errl))
        print ('')
        print ('The best fit planet parameters are:')
        print ('T0    = %4.4f + %4.4f - %4.4f days'%(t0_val[l],t0_errr[l],t0_errl[l]))
        print ('Tp    = %4.4f + %4.4f - %4.4f days'%(tp_val[l],tp_errr[l],tp_errl[l]))
        print ('P     = %4.4f + %4.4f - %4.4f days'%(P_val[l], P_errr[l], P_errl[l]))
        print ('e     = %4.4f + %4.4f - %4.4f     '%(e_val[l], e_errr[l],e_errl[l]))
        print ('w     = %4.4f + %4.4f - %4.4f deg ' %(w_deg[l],w_deg_errr[l], w_deg_errl[l]))
        if (fit_rv):
          print ('K     = %4.4f + %4.4f - %4.4f m/s'%(k_val[l]/1.e-3,(k_errr[l])/1.e-3, (k_errl[l])/1e-3))
          print ('mpsin = %4.4f + %4.4f - %4.4f (%s masses) '%(mass_val[l],mass_errr[l], mass_errl[l], unit_mass))
          for i in range(0,nt):
            print ('%s v0  = %4.4f + %4.4f - %4.4f km/s'%(telescopes[i], \
          v_val[i],v_errr[i],v_errl[i]))


  #Run the previous function
  print_errors_planets()	

