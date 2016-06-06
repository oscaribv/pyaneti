from matplotlib import gridspec


#what is the minimum chi2
minchi2_index = np.argmin(chi2)

if ( nplanets == 1 ):

  #If some parameters are transformed let us go back

  if (is_log_P):
    Po = np.power(10.,Po)

  if (is_ew):
    dummy_e = eo
    eo = eo * eo + wo * wo
    wo = np.arctan2(dummy_e,wo)

  if ( fit_tr ):

    if (is_sini):
      io = np.arcsin(io)
    inclination = io
    if (is_log_a):
      ao = np.power(10.,ao)

    #Take back the u1 and u2 values, Kipping 2013
    q1o = u1o
    q2o = u2o
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

    t0_val, t0_errl, t0_errr = find_vals_perc(t0o,nconv,s_factor)
    T0 = t0_val
    P_val, P_errl, P_errr = find_vals_perc(Po,nconv,s_factor)
    e_val,e_errl, e_errr = find_vals_perc(eo,nconv,s_factor)
    w_val,w_errl, w_errr = find_vals_perc(wo,nconv,s_factor)
    if (w_val < 0.0 ):
      w_val = w_val + 2 * np.pi	
    #Transform periastron from rad to degrees
    w_deg = w_val * 180. / np.pi
    w_deg_errl = w_errl * 180. / np.pi
    w_deg_errr = w_errr * 180. / np.pi

    if ( fit_tr ):
      i_val,i_errl, i_errr = find_vals_perc(io,nconv,s_factor)
      a_val,a_errl, a_errr = find_vals_perc(ao,nconv,s_factor)
      aphy_val,aphy_errl, aphy_errr = find_vals_perc(aphyo,nconv,s_factor)
      q1_val,q1_errl, q1_errr = find_vals_perc(q1o,nconv,s_factor)
      q2_val,q2_errl, q2_errr = find_vals_perc(q2o,nconv,s_factor)
      u1_val,u1_errl, u1_errr = find_vals_perc(u1o,nconv,s_factor)
      u2_val,u2_errl, u2_errr = find_vals_perc(u2o,nconv,s_factor)
      pz_val,pz_errl, pz_errr = find_vals_perc(pzo,nconv,s_factor)
      rp_val,rp_errl, rp_errr = find_vals_perc(rpo,nconv,s_factor)
      b_val , b_errl, b_errr  = find_vals_perc(bo,nconv,s_factor)
      tt_val , tt_errl, tt_errr  = find_vals_perc(tto,nconv,s_factor)
      tf_val , tf_errl, tf_errr  = find_vals_perc(tfo,nconv,s_factor)
      rho_val , rho_errl, rho_errr  = find_vals_perc(rhoo,nconv,s_factor)
      rhop_val,rhop_errl, rhop_errr = find_vals_perc(rho_p,nconv,s_factor)
      gp_val,gp_errl, gp_errr = find_vals_perc(gpo,nconv,s_factor)
      i_deg = i_val * 180. / np.pi
      i_deg_errl = i_errl * 180. / np.pi
      i_deg_errr = i_errr * 180. / np.pi
    if ( fit_rv ):
      k_val, k_errl, k_errr = find_vals_perc(ko,nconv,s_factor)
      m_val, m_errl, m_errr  	= find_vals_perc(masso,nconv,s_factor)
      tp_val, tp_errl, tp_errr= find_vals_perc(tpo,nconv,s_factor)
      #Systemic velocities are stored in v_val
      v_val = [None]*nt
      v_errl = [None]*nt
      v_errr = [None]*nt
      for j in range(0,nt):
        v_val[j], v_errl[j], v_errr[j] = find_vals_perc(vo[j],nconv,s_factor)
	

    npln = npars * np.log(ndata)	

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
    print ('BIC        = %1.4f' %(chi2tot_val + npln))
    print ('')
    print ('The best fit planet parameters are:')
    print ('T0    = %4.7f + %4.7f - %4.7f days'%(t0_val,t0_errr,t0_errl))
    print ('P     = %4.7f + %4.7f - %4.7f days'%(P_val, P_errr , P_errl))
    print ('e     = %4.4f + %4.4f - %4.4f     '%(e_val, e_errr , e_errl))
    print ('w     = %4.4f + %4.4f - %4.4f deg '%(w_deg,w_deg_errr, w_deg_errl))
    if (fit_tr):
      print ('Transit fit parameters:')
      print ('i     = %4.4f + %4.4f - %4.4f deg' %(i_deg,i_deg_errr, i_deg_errl))
      print ('a/r*  = %4.4f + %4.4f - %4.4f    ' 		%(a_val, a_errr , a_errl))
      print ('rp/r* = %4.4f + %4.4f - %4.4f    ' 		%(pz_val,pz_errr, pz_errl))
      print ('q_1    = %4.4f + %4.4f - %4.4f    ' 		%(q1_val,q1_errr, q1_errl))
      print ('q_2    = %4.4f + %4.4f - %4.4f    ' 		%(q2_val,q2_errr, q2_errl))
    if (fit_rv):
      print ('RV fit parameters:')
      print ('K     = %4.4f + %4.4f - %4.4f m/s' 		%(k_val/1.e-3,(k_errr)/1.e-3, (k_errl)/1.e-3))
      for i in range(0,nt):
        print ('%s v0  = %4.4f + %4.4f - %4.4f km/s' 	%(telescopes[i], \
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
    if (fit_rv):
      print ('Tp    = %4.4f + %4.4f - %4.4f days' 		%(tp_val,tp_errr, tp_errl))
      print ('mp    = %4.4f + %4.4f - %4.4f %s masses' 		%(m_val,m_errr, m_errl,unit_mass))
      if ( fit_tr ):
        print ('rho_p = %4.4f + %4.4f - %4.4f g/cm^3' 		%(rhop_val,rhop_errr, rhop_errl))
        print ('g_p = %4.4f + %4.4f - %4.4f cm/s^2' 		%(gp_val,gp_errr, gp_errl))


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

			masso = planet_mass(mstar,ko[l]*1.e3,Po[l],eo[l],inclination)
			if ( unit_mass == 'earth'):
				#masso[m] = 332946 * masso[m]
				masso = 332967.750577677 * masso
			elif ( unit_mass == 'jupiter'):
				#masso[m] = 1047.56 * masso[m]
				masso = 1047.353299069 * mass
	
			if ( errores == 'gauss' ):
	
				chi2_val[l], chi2_errs[l] = find_vals_gauss(chi2red[l],nconv)
				t0_val[l],t0_err[l] = find_vals_gauss(t0o[l],nconv)
				tp_val[l],tp_err[l] = find_vals_gauss(tpo[l],nconv)
				mass_val[l],mass_err[l] = find_vals_gauss(masso[l],nconv)
				P_val[l], P_err[l]  = find_vals_gauss(Po[l],nconv)
				e_val[l],e_err[l] 	= find_vals_gauss(eo[l],nconv)
				w_val[l],w_err[l] 	= find_vals_gauss(wo[l],nconv)
				if (w_val[l] < 0.0 ):
					w_val[l] = w_val[l] + 2 * np.pi	
				w_deg[l] 		= w_val[l] * 180. / np.pi
				w_deg_err[l] = w_err[l] * 180. / np.pi
				if ( fit_rv ):
					k_val[l], k_err[l]  = find_vals_gauss(ko[l],nconv)
					v_val[l] = [None]*nt
					v_err[l] = [None]*nt
					for j in range(0,nt):
						v_val[j], v_err[j] = find_vals_gauss(vo[j],nconv)
		
				#Print the best fit values values
				print ('chi2_red = %1.4f +/- %1.4f' %(chi2_val[l],chi2_errs[l]))
				print ('The best fit planet parameters are:')
				print ('T0    = %4.4f +/- %4.4f days'%(t0_val[l],t0_err[l]))
				print ('Tp    = %4.4f +/- %4.4f days'%(tp_val[l],tp_err[l]))
				print ('P     = %4.4f +/- %4.4f days' 		%(P_val[l],P_err[l]))
				print ('e     = %4.4f +/- %4.4f     '			%(e_val[l],e_err[l]))
				print ('w     = %4.4f +/- %4.4f deg '	%(w_deg[l],w_deg_err[l]))
				if (fit_rv):
					print ('K    = %4.4f +/- %4.4f m/s' 		%(k_val[l]/1e-3,k_err[l]/1e-3))
					print ('mpsin = %4.4f +/- %4.4f m/s' 		%(mass_val[l],mass_err[l]))
					for i in range(0,nt):
						print ('%s v0 = %4.4f +/- %4.4f km/s' 	%(telescopes[i], \
						v_val[i],v_err[i]))
			
			#Percentile errors
		
			if ( errores == 'perc' ):	

				s_factor = 1.0
	
				#chi2tot_val[l], chi2tot_errl[l], chi2tot_errr[l] = find_vals_perc(chi2[l],nconv,s_factor)
				#chi2_val[l], chi2_errl[l], chi2_errr[l] = find_vals_perc(chi2red[l],nconv,s_factor)
				chi2tot_val[l] = np.amin(chi2[l])
                                chi2_val = chi2tot_vali[l] / ( ndata - npars )
				#chi2_val[l]    = np.amin(chi2red[l])

				if ( scale_error_bars ):
					s_factor = np.sqrt( chi2_val[l] )
					if ( chi2_val > 1.0 ):
						s_factor = 1.0 / s_factor
				else:
					s_factor = 1.0

	
				t0_val[l], t0_errl[l], t0_errr[l] = find_vals_perc(t0o[l],nconv,s_factor)
				mass_val[l], mass_errl[l], mass_errr[l] = find_vals_perc(masso,nconv,s_factor)
				tp_val[l], tp_errl[l], tp_errr[l] = find_vals_perc(tpo,nconv,s_factor)
				P_val[l], P_errl[l], P_errr[l]  = find_vals_perc(Po[l],nconv,s_factor)
				e_val[l],e_errl[l], e_errr[l] 	= find_vals_perc(eo[l],nconv,s_factor)
				w_val[l],w_errl[l], w_errr[l] 	= find_vals_perc(wo[l],nconv,s_factor)
				if (w_val[l] < 0.0 ):
					w_val[l] = w_val[l] + 2 * np.pi	
				w_deg[l] 		= w_val[l] * 180. / np.pi
				w_deg_errl[l] = w_errl[l] * 180. / np.pi
				w_deg_errr[l] = w_errr[l] * 180. / np.pi
				if ( fit_rv ):
					k_val[l], k_errl[l], k_errr[l]  = find_vals_perc(ko[l],nconv,s_factor)
					for j in range(0,nt):
						v_val[j], v_errl[j], v_errr[j] = find_vals_perc(vo[j],nconv,s_factor)
			
		
				#Print the best fit values values
				print 'N_data      = ', ndata
				print 'N_pars      = ', npars
				print 'DOF         = ', ndata - npars
				print 'scale factor= ', s_factor
				print ('chi2 = %1.4f' %(chi2tot_val[l]))
				print ('chi2_red = %1.4f' %(chi2_val[l]))
				print ('BIC        = %1.4f ' %(chi2tot_val[l] + npln))
				print ('The best fit planet parameters are:')
				print ('T0    = %4.4f + %4.4f - %4.4f days'%(t0_val[l],t0_errr[l],t0_errl[l]))
				print ('Tp    = %4.4f + %4.4f - %4.4f days'%(tp_val[l],tp_errr[l],tp_errl[l]))
				print ('P     = %4.4f + %4.4f - %4.4f days' 		%(P_val[l], P_errr[l], P_errl[l]))
				print ('e     = %4.4f + %4.4f - %4.4f     '			%(e_val[l], e_errr[l],e_errl[l]))
				print ('w     = %4.4f + %4.4f - %4.4f deg '	%(w_deg[l],w_deg_errr[l], w_deg_errl[l]))
				if (fit_rv):
					print ('K     = %4.4f + %4.4f - %4.4f m/s' 		%(k_val[l]/1.e-3,(k_errr[l])/1.e-3, (k_errl[l])/1e-3))
					print ('mpsin = %4.4f + %4.4f - %4.4f (%s masses) '	%(mass_val[l],mass_errr[l], mass_errl[l], unit_mass))
					for i in range(0,nt):
						print ('%s v0  = %4.4f + %4.4f - %4.4f km/s' 	%(telescopes[i], \
							v_val[i],v_errr[i],v_errl[i]))

	#Run the previous function
	print_errors_planets()	

