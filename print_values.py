from matplotlib import gridspec


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
	  if (is_log_a):
	    ao = np.power(10.,ao)

	if ( fit_rv ):
	  if ( is_log_k ):
	    ko = np.power(10.,ko)
	  for j in range(0,nt):
	    if ( is_log_rv0 ):
	      vo[j] = np.power(10.,vo[j])


	if ( errores == 'gauss' ):

		chi2_val, chi2_errs = find_vals_gauss(chi2red,nconv)
		t0_val,t0_err = find_vals_gauss(t0o,nconv)
		P_val, P_err  = find_vals_gauss(Po,nconv)
		e_val,e_err 	= find_vals_gauss(eo,nconv)
		w_val,w_err 	= find_vals_gauss(wo,nconv)
		if (w_val < 0.0 ):
			w_val = w_val + 2 * np.pi	
		w_deg 		= w_val * 180. / np.pi
		w_deg_err = w_err * 180. / np.pi
		if ( fit_tr ):
			i_val,i_err 	= find_vals_gauss(io,nconv)
			a_val,a_err 	= find_vals_gauss(ao,nconv)
			u1_val,u1_err = find_vals_gauss(u1o,nconv)
			u2_val,u2_err = find_vals_gauss(u2o,nconv)
			pz_val,pz_err = find_vals_gauss(pzo,nconv)
			i_deg 		= i_val * 180. / np.pi
			i_deg_err = i_err * 180. / np.pi
		if ( fit_rv ):
			k_val, k_err  = find_vals_gauss(ko,nconv)
			v_val = [None]*nt
			v_err = [None]*nt
			for j in range(0,nt):
				v_val[j], v_err[j] = find_vals_gauss(vo[j],nconv)

		t0_errl = t0_err
		t0_errr = t0_err
		P_errl = P_err
		P_errr = P_err
		e_errl = e_err
		e_errr = e_err 
		w_errl = w_err
		w_errr = w_err
		k_errl = k_err
		k_errr = k_err
		v_errl = [None]*nt
		v_errr = [None]*nt
		for m in range(0,nt):
			v_errl[m] = v_err[m]
			v_errr[m] = v_err[m]


		#Print the best fit values values
		print ('chi2_red = %1.4f +/- %1.4f' %(chi2_val,chi2_errs))
		print ('The best fit planet parameters are:')
		print ('T0    = %4.4f +/- %4.4f days'%(t0_val,t0_err))
		print ('P     = %4.4f +/- %4.4f days' 		%(P_val,P_err))
		print ('e     = %4.4f +/- %4.4f     '			%(e_val,e_err))
		print ('w     = %4.4f +/- %4.4f deg '	%(w_deg,w_deg_err))
		if (fit_tr):
			print ('i     = %4.4f +/- %4.4f deg' %(i_deg,i_deg_err))
			print ('a/r*  = %4.4f +/- %4.4f    ' 		%(a_val,a_err))
			print ('u1    = %4.4f +/- %4.4f    ' 		%(u1_val,u1_err))
			print ('u2    = %4.4f +/- %4.4f    ' 		%(u2_val,u2_err))
			print ('rp/r* = %4.4f +/- %4.4f    ' 		%(pz_val,pz_err))
		if (fit_rv):
			print ('K    = %4.4f +/- %4.4f m/s ' 		%(k_val/1.e-3,k_err/1e-3))
			for i in range(0,nt):
				print ('%s v0 = %4.4f +/- %4.4f km/s' 	%(telescopes[i], \
				v_val[i],v_err[i]))
	
	#Percentile errors

	if ( errores == 'perc' ):

		chi2_val, chi2_errl, chi2_errr = find_vals_perc(chi2red,nconv)

		t0_val, t0_errl, t0_errr = find_vals_perc(t0o,nconv)
		P_val, P_errl, P_errr  = find_vals_perc(Po,nconv)
		e_val,e_errl, e_errr 	= find_vals_perc(eo,nconv)
		w_val,w_errl, w_errr 	= find_vals_perc(wo,nconv)
		if (w_val < 0.0 ):
			w_val = w_val + 2 * np.pi	
		w_deg 		= w_val * 180. / np.pi
		w_deg_errl = w_errl * 180. / np.pi
		w_deg_errr = w_errr * 180. / np.pi
		if ( fit_tr ):
			i_val,i_errl, i_errr 	= find_vals_perc(io,nconv)
			a_val,a_errl, a_errr 	= find_vals_perc(ao,nconv)
			u1_val,u1_errl, u1_errr = find_vals_perc(u1o,nconv)
			u2_val,u2_errl, u2_errr = find_vals_perc(u2o,nconv)
			pz_val,pz_errl, pz_errr = find_vals_perc(pzo,nconv)
			i_deg 		= i_val * 180. / np.pi
			i_deg_errl = i_errl * 180. / np.pi
			i_deg_errr = i_errr * 180. / np.pi
		if ( fit_rv ):
			k_val, k_errl, k_errr  = find_vals_perc(ko,nconv)
			v_val = [None]*nt
			v_errl = [None]*nt
			v_errr = [None]*nt
			for j in range(0,nt):
				v_val[j], v_errl[j], v_errr[j] = find_vals_perc(vo[j],nconv)
	


		#Print the best fit values values
		print ('chi2_red = %1.4f + %1.4f - %1.4f' %(chi2_val,chi2_errr,chi2_errl))
		print ('The best fit planet parameters are:')
		print ('T0    = %4.4f + %4.4f - %4.4f days'%(t0_val,t0_errr,t0_errl))
		print ('P     = %4.4f + %4.4f - %4.4f days'%(P_val, P_errr , P_errl))
		print ('e     = %4.4f + %4.4f - %4.4f     '%(e_val, e_errr , e_errl))
		print ('w     = %4.4f + %4.4f - %4.4f deg '%(w_deg,w_deg_errr, w_deg_errl))
		if (fit_tr):
			print ('i     = %4.4f + %4.4f - %4.4f deg' %(i_deg,i_deg_errr, i_deg_errl))
			print ('a/r*  = %4.4f + %4.4f - %4.4f    ' 		%(a_val, a_errr , a_errl))
			print ('u1    = %4.4f + %4.4f - %4.4f    ' 		%(u1_val,u1_errr, u1_errl))
			print ('u2    = %4.4f + %4.4f - %4.4f    ' 		%(u2_val,u2_errr, u2_errl))
			print ('rp/r* = %4.4f + %4.4f - %4.4f    ' 		%(pz_val,pz_errr, pz_errl))
	if (fit_rv):
			print ('K     = %4.4f + %4.4f - %4.4f m/s' 		%(k_val/1.e-3,(k_errr)/1.e-3, (k_errl)/1.e-3))
			for i in range(0,nt):
				print ('%s v0  = %4.4f + %4.4f - %4.4f km/s' 	%(telescopes[i], \
				v_val[i],v_errr[i],v_errl[i]))

#Multiplanet fit
else:


	#Define global variables
	mass_val = [None]*nplanets
	mass_errr = [None]*nplanets
	mass_errl = [None]*nplanets
	chi2_val = [None]*nplanets
	chi2_errr = [None]*nplanets
	chi2_errl = [None]*nplanets
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
				masso[m] = planet_mass(mstar,ko[l][m]*1.e3,Po[l][m],eo[l][m])
				masso[m] = 332946 * masso[m]
				
	
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
	
				chi2_val[l], chi2_errl[l], chi2_errr[l] = find_vals_perc(chi2red[l],nconv)
		
				t0_val[l], t0_errl[l], t0_errr[l] = find_vals_perc(t0o[l],nconv)
				mass_val[l], mass_errl[l], mass_errr[l] = find_vals_perc(masso,nconv)
				tp_val[l], tp_errl[l], tp_errr[l] = find_vals_perc(tpo,nconv)
				P_val[l], P_errl[l], P_errr[l]  = find_vals_perc(Po[l],nconv)
				e_val[l],e_errl[l], e_errr[l] 	= find_vals_perc(eo[l],nconv)
				w_val[l],w_errl[l], w_errr[l] 	= find_vals_perc(wo[l],nconv)
				if (w_val[l] < 0.0 ):
					w_val[l] = w_val[l] + 2 * np.pi	
				w_deg[l] 		= w_val[l] * 180. / np.pi
				w_deg_errl[l] = w_errl[l] * 180. / np.pi
				w_deg_errr[l] = w_errr[l] * 180. / np.pi
				if ( fit_rv ):
					k_val[l], k_errl[l], k_errr[l]  = find_vals_perc(ko[l],nconv)
					for j in range(0,nt):
						v_val[j], v_errl[j], v_errr[j] = find_vals_perc(vo[j],nconv)
			
		
				#Print the best fit values values
				print ('chi2_red = %1.4f + %1.4f - %1.4f' %(chi2_val[l],chi2_errr[l],chi2_errl[l]))
				print ('The best fit planet parameters are:')
				print ('T0    = %4.4f + %4.4f - %4.4f days'%(t0_val[l],t0_errr[l],t0_errl[l]))
				print ('Tp    = %4.4f + %4.4f - %4.4f days'%(tp_val[l],tp_errr[l],tp_errl[l]))
				print ('P     = %4.4f + %4.4f - %4.4f days' 		%(P_val[l], P_errr[l], P_errl[l]))
				print ('e     = %4.4f + %4.4f - %4.4f     '			%(e_val[l], e_errr[l],e_errl[l]))
				print ('w     = %4.4f + %4.4f - %4.4f deg '	%(w_deg[l],w_deg_errr[l], w_deg_errl[l]))
				if (fit_rv):
					print ('K     = %4.4f + %4.4f - %4.4f m/s' 		%(k_val[l]/1.e-3,(k_errr[l])/1.e-3, (k_errl[l])/1e-3))
					print ('mpsin = %4.4f + %4.4f - %4.4f (Earth masses) '	%(mass_val[l],mass_errr[l], mass_errl[l]))
					for i in range(0,nt):
						print ('%s v0  = %4.4f + %4.4f - %4.4f km/s' 	%(telescopes[i], \
							v_val[i],v_errr[i],v_errl[i]))

	#Run the previous function
	print_errors_planets()	

