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

	#Print the best fit values values
	print ('chi2_red = %1.4e +/- %1.4e' %(chi2_val,chi2_errs))
	print ('The best fit planet parameters are:')
	print ('T0    = %4.4e +/- %4.4e days'%(t0_val,t0_err))
	print ('P     = %4.4e +/- %4.4e' 		%(P_val,P_err))
	print ('e     = %4.4e +/- %4.4e'			%(e_val,e_err))
	print ('w     = %4.4e +/- %4.4e deg'	%(w_deg,w_deg_err))
	if (fit_tr):
		print ('i     = %4.4e +/- %4.4e deg' %(i_deg,i_deg_err))
		print ('a/r*  = %4.4e +/- %4.4e' 		%(a_val,a_err))
		print ('u1    = %4.4e +/- %4.4e' 		%(u1_val,u1_err))
		print ('u2    = %4.4e +/- %4.4e' 		%(u2_val,u2_err))
		print ('rp/r* = %4.4e +/- %4.4e' 		%(pz_val,pz_err))
	if (fit_rv):
		print ('K    = %4.4e +/- %4.4e' 		%(k_val,k_err))
		for i in range(0,nt):
			print ('%s v0 = %4.4e +/- %4.4e' 	%(telescopes[i], \
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
	print ('chi2_red = %1.4e + %1.4e - %1.4e' %(chi2_val,chi2_errr-chi2_val, chi2_val - chi2_errl))
	print ('The best fit planet parameters are:')
	print ('T0    = %4.4e + %4.4e - %4.4e days'%(t0_val,t0_errr-t0_val, t0_val-t0_errl))
	print ('P     = %4.4e + %4.4e - %4.4e' 		%(P_val, P_errr - P_val, P_val - P_errl))
	print ('e     = %4.4e + %4.4e - %4.4e'			%(e_val, e_errr - e_val, e_val - e_errl))
	print ('w     = %4.4e + %4.4e - %4.4e deg'	%(w_deg,w_deg_errr - w_deg, w_deg - w_deg_errl))
	if (fit_tr):
		print ('i     = %4.4e + %4.4e - %4.4e deg' %(i_deg,i_deg_errr-i_deg, i_deg - i_deg_errl))
		print ('a/r*  = %4.4e + %4.4e - %4.4e' 		%(a_val, a_errr - a_val ,a_val - a_errl))
		print ('u1    = %4.4e + %4.4e - %4.4e' 		%(u1_val,u1_errr - u1_val, u1_val - u1_errl))
		print ('u2    = %4.4e + %4.4e - %4.4e' 		%(u2_val,u2_errr - u2_val, u2_val - u2_errl))
		print ('rp/r* = %4.4e + %4.4e - %4.4e' 		%(pz_val,pz_errr - pz_val, pz_val - pz_errl))
if (fit_rv):
		print ('K     = %4.4e + %4.4e - %4.4e' 		%(k_val,k_errr - k_val, k_val - k_errl))
		for i in range(0,nt):
			print ('%s v0  = %4.4e + %4.4e - %4.4e' 	%(telescopes[i], \
			v_val[i],v_errr[i] - v_val[i], v_val[i] - v_errl[i]))

