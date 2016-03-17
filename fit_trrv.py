#!/usr/bin/python2

#Load libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm, sigmaclip
import sys
import pyaneti as pti

#Read the file with all the python functions
execfile('todo-py.py')

#Read the file with the default values
execfile('default.py')

#Read input file
execfile('input_fit.py')

#PREPATARION RV DATA
if (fit_rv):

	#Read the data file
	time,rv,err,tspe = np.loadtxt(fname_rv,usecols=(0,1,2,3), \
  	dtype={'names': ('time', 'rv', 'err','telescope'), \
		'formats': ('float', 'float', 'float', 'S1')}, \
		comments='#',unpack=True)

	#Transform rv from km/s to m/s
	if(units_ms):
		ktom = 1000
		rv=rv*ktom
		err=err*ktom
		ylab = 'RV (m/s)'
	
	#These lists have lists with data for the different telescopes
	time_all=[]
	rv_all=[]
	errs_all=[]

	#Number of telescopes
	nt = len(telescopes)

	#tspe is the data with the telescopes read from the file
	for i in range(0,nt):
		time_dum =[]
		rv_dum =[]
		errs_dum =[]
		if (len(telescopes[i]) == 0):
			sys.exit("There is no data for %s"% telescopes[i])
		else:
			for j in range(0,len(tspe)):
				if (tspe[j] == telescopes[i]):
					time_dum.append(time[j])
					rv_dum.append(rv[j])
					errs_dum.append(err[j])
		#The *all variables are lists of lists, each list constains
		# a list with the data of each telescope
			time_all.append(time_dum)
			rv_all.append(rv_dum)
			errs_all.append(errs_dum)

	#The mega* variables contains all telescope data
	mega_rv = []
	mega_time = []
	mega_err  = []
	tlab = []
	for i in range(0,nt): 
		for j in range(0,len(rv_all[i])):
			tlab.append(i)
			mega_rv.append(rv_all[i][j])
			mega_time.append(time_all[i][j])
			mega_err.append(errs_all[i][j])

#The RV data is ready

#PREPARATION TRANSIT DATA

if (fit_tr):

	#Read the data file
	dummyd,dummyf,flag = np.loadtxt(fname_tr,usecols=(2,9,10), \
	comments='\\',unpack=True)

	#Let us take the good data with the flag
	nobin_wflux = []
	nobin_hdate = []
	for i in range(0,len(flag)):
		if ( flag[i] == 0):
			nobin_wflux.append(dummyf[i])
			nobin_hdate.append(dummyd[i])

	nbin = 16
	#nbin = 8

	#bin the data to do fastest test
	hdate, err_hdate = bin_data(nobin_hdate,nbin)
	wflux, errs = bin_data(nobin_wflux,nbin)

	#THIS HAS TO BE DONE AUTOMATICALLY!	

	ntr = 2

	tls = [None]*ntr

	tls[0] = [3217.,3218.]
	tls[1] = [3232.1,3233.1]

	#crash if you do not have more than one transit
	if ( ntr < 2):
		print "you do not have enought transit data!"
		sys.exit("I crashed because I want more data!")
	
	#Each element of these vectors will have the information
	#of a given transit
	xt= [None]*ntr	
	yt= [None]*ntr	
	et= [None]*ntr	

	#Normalize all the transit independently
	for i in range(0,ntr):
		xt[i],yt[i],et[i] = normalize_transit(hdate,wflux,errs,tls[i])

	#Let us put together the information of all the arrays
	megax = np.concatenate(xt)
	megay = np.concatenate(yt)
	megae = np.concatenate(et)

#TRANSIT DATA READY

#PRIORS SECTION

#Let us try to do a guess for the init values
if (fit_rv):
	k0vecmax = [None]*nt
	k0vecmin = [None]*nt
	k0 = 0.0
	for i in range(0,nt):
		k0vecmin[i] = min(rv_all[i])
		k0vecmax[i] = max(rv_all[i])
		k0 = k0 + ( k0vecmax[i] - k0vecmin[i] ) / 2

	k0 = k0 / nt
	v0 = np.zeros(nt)

	for i in range(0,nt):
		v0[i] = ( k0vecmin[i] + k0vecmax[i] ) / 2.0

if (fit_tr):
	T0 = min(xt[0]) + 0.5*(max(xt[0])-min(xt[0]))

if ( is_circular ):
	fit_e = False
	fit_w = False
	e = 0.0
	w = np.pi / 2.0

#Print intial configuration
print_init()

#-------------------------------------------------------------
#			FITTING ROUTINES
#-------------------------------------------------------------

#FIT TRANSIT AND RV CURVES
if (fit_rv and fit_tr ):

	flag = [is_log_P,is_ew,is_sini,is_log_a,is_log_k,is_log_rv0]

	what_fit = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w), \
              int(fit_i),int(fit_a),int(fit_u1),int(fit_u2),\
              int(fit_pz), int(fit_k), int(fit_v0)]
	dummy = [T0,P,e,w,ii,a,u1,u2,pz,k0]
	params = np.concatenate((dummy,v0))

	#Call the fit routine
	#pti.metropolis_hastings(mega_time,mega_rv,mega_err,tlab \
	#,megax, megay, megae, params, prec, maxi, thin_factor, \
	#is_circular, what_fit, flag, nconv)

	min_t0	= min(xt[0])
	max_t0 	= max(xt[0])
	min_P	 	= 13.
	max_P	 	= 16.
	min_e		= 1.e-8		
	max_e		= 0.5
	min_w		= 0.0
	max_w		= 2*np.pi
	min_i		= 0.
	max_i		= np.pi / 2.0
	min_a		= 5.0
	max_a		= 20.0
	min_u1	= 0.0
	max_u1	= 0.5
	min_u2	= 0.0
	max_u2	= 0.5
	min_pz	= 1e-3
	max_pz	= 0.5
	min_k		= 5
	max_k		= 20
	min_rv0	= 30
	max_rv0 = 70

	vec_rv0_limits = []
	for m in range(0,nt):
		vec_rv0_limits.append(min_rv0) 
		vec_rv0_limits.append(max_rv0) 

	dummy_lims = \
	[	min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w \
		, min_i, max_i, min_a, max_a, min_u1, max_u1, min_u2, \
		max_u2, min_pz, max_pz, min_k, max_k]

	limits = np.concatenate((dummy_lims,vec_rv0_limits)) 
	
	nwalks = 20*len(params)

	pti.stretch_move(mega_time,mega_rv,mega_err,tlab \
	,megax, megay, megae, params,limits, nwalks, prec, maxi, thin_factor, \
	is_circular, what_fit, flag, nconv)

	#Read the data
	vari,chi2,chi2red,t0o,Po,eo,wo,io,ao,u1o,u2o,pzo,ko =  \
	np.loadtxt('mh_fit.dat', comments='#',unpack=True, \
	usecols=range(0,13))
	vo = [None]*nt
	for j in range(0,nt):
		n = [13+j]
		a = np.loadtxt('mh_fit.dat', comments='#', \
		unpack=True, usecols=(n))
		vo[j] = a

#FIT TRANSIT CURVE ONLY
elif ( not fit_rv and fit_tr ):

	flag = [is_log_P, is_ew, is_sini, is_log_a]

	what_fit = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w),  \
              int(fit_i),int(fit_a), int(fit_u1),int(fit_u2),\
            int(fit_pz)]

	params = [T0,P,e,w,ii,a,u1,u2,pz]

	#Call fit routine
	#pti.metropolis_hastings_tr(megax, megay, megae,  \
	#params, prec, maxi, thin_factor, is_circular, what_fit,flag,nconv)

	min_t0	= min(xt[0])
	max_t0 	= max(xt[0])
	min_P	 	= 13.
	max_P	 	= 16.
	min_e		= 1.e-8		
	max_e		= 0.5
	min_w		= 0.0
	max_w		= 2*np.pi
	min_i		= 0.
	max_i		= 1*np.pi
	min_a		= 5.0
	max_a		= 20.0
	min_u1	= 0.0
	max_u1	= 0.5
	min_u2	= 0.0
	max_u2	= 0.5
	min_pz	= 1e-3
	max_pz	= 0.5

	limits = \
	[	min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w \
		, min_i, max_i, min_a, max_a, min_u1, max_u1, \
		min_u2, max_u2, min_pz, max_pz]

	nwalks = 20*len(params)
	nwalks = 50
	
	pti.stretch_move_tr(megax, megay, megae,  \
	params,limits, nwalks, prec, maxi, thin_factor, is_circular, what_fit,flag,nconv)

	#Read the data
	vari, chi2,chi2red,t0o,Po,eo,wo,io,ao,u1o,u2o,pzo = \
       np.loadtxt('mh_trfit.dat', comments='#',unpack=True)

#FIT RV CURVE ONLY
elif ( fit_rv and not fit_tr ):

	flag = [is_log_P,is_ew,is_log_k,is_log_rv0]

	if ( is_circular ):
		fit_e = False
		fit_w = False

	what_fit_p1 = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w), \
					    int(fit_k), int(fit_v0)]
	what_fit_p2 = [ 0, 0,int(fit_e),int(fit_w), \
					    int(fit_k), int(fit_v0)]
	dummy_p1 = [6813.38345,7.919454,e,w,k0]
	dummy_p2 = [6817.2759,11.90701,e,w,k0]
	params_p1 = np.concatenate((dummy_p1,v0))
	params_p2 = np.concatenate((dummy_p2,v0))

	min_t0	= 6812.
	max_t0 	= 6820.
	min_P	 	= 7.5
	max_P	 	= 8.5
	min_e		= 1.e-10		
	max_e		= 0.999
	min_w		= 0.0
	max_w		= 2*np.pi
	min_k		= 0.00001
	max_k		= 2.0
	min_rv0	= 1.
	max_rv0 = 10.

	vec_rv0_limits = []
	for m in range(0,nt):
		vec_rv0_limits.append(min_rv0) 
		vec_rv0_limits.append(max_rv0) 

	dummy_lims_p2 = \
	[	6800., 6830., min_P, 3*max_P, min_e, max_e, min_w, max_w \
		, min_k, max_k]

	dummy_lims_p1 = \
	[	min_t0, max_t0, min_P, max_P, min_e, max_e, min_w, max_w \
		, min_k, max_k]

	limits_p1 = np.concatenate((dummy_lims_p1,vec_rv0_limits)) 
	limits_p2 = np.concatenate((dummy_lims_p2,vec_rv0_limits)) 

	nplanets = 1
	
	#Call fit routine
	#pti.metropolis_hastings_rv(mega_time,mega_rv,mega_err,tlab,\
  #params, prec, maxi, thin_factor, is_circular, what_fit,flag,nconv)

	#list_pars = np.concatenate((params_p1,params_p2))
	#list_lims = np.concatenate((limits_p1,limits_p2))
	#list_wtfs = np.concatenate((what_fit_p1,what_fit_p2))

	list_pars = (params_p1)
	list_lims = (limits_p1)
	list_wtfs = (what_fit_p1)

	#print list_pars

	#print list_lims
	
	#print list_wtfs	

	#sys.exit('O')

	nwalkers = 20 * len(params_p1)
	nwalkers = 20

	out_file = 'planet1.dat'
	#out_file = ['planet1.dat','planet2.dat']

	pti.stretch_move_rv(mega_time,mega_rv,mega_err,tlab,\
  list_pars, list_lims, nwalkers, prec, maxi, thin_factor, \
	is_circular, list_wtfs,flag,nconv,datas=len(mega_time), \
	nt=nt,npl=nplanets)
	#Read the data
	nconv = nconv * (nwalkers-1)

	if ( nplanets == 1 ):
		vari,chi2,chi2red,t0o,Po,eo,wo,ko = \
		np.loadtxt(out_file, comments='#', unpack=True,\
		usecols=range(0,8))
		vo = [None]*nt
		for j in range(0,nt):
			n = [8+j]
			a = np.loadtxt(out_file, comments='#',unpack=True, usecols=(n))
			vo[j] = a
	else:
		#Create all the variables, list of lists
		vari = [[]]*nplanets
		chi2 = [[]]*nplanets
		chi2red = [[]]*nplanets
		t0o = [[]]*nplanets
		Po = [[]]*nplanets
		eo = [[]]*nplanets
		wo = [[]]*nplanets
		ko = [[]]*nplanets
		#Define global variables
		chi2_val = [None]*nplanets
		chi2_errr = [None]*nplanets
		chi2_errl = [None]*nplanets
		t0_val = [None]*nplanets
		t0_errr = [None]*nplanets
		t0_errl = [None]*nplanets
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
		for l in range(0,nplanets):
			vari[l],chi2[l],chi2red[l],t0o[l],Po[l],eo[l],wo[l],ko[l] = \
			np.loadtxt(out_file[l], comments='#', unpack=True,\
			usecols=range(0,8))
		#The  systemic velocities are the same for all the planets
		vo = [None]*nt
		for j in range(0,nt):
			n = [8+j]
			a = np.loadtxt(out_file[0], comments='#',unpack=True, usecols=(n))
			vo[j] = a
		

#Nothing to fit!
else:
	sys.exit("Nothing to fit!")

#-------------------------------------------------------------
#			END FITTING ROUTINES
#-------------------------------------------------------------
	
if ( nplanets == 1 ):

	errores = 'perc'
	execfile('print_errors.py')

	#PLOT TRANSIT
	if ( fit_tr ):
		plot_transit()

	#PLOT RV CURVE
	if ( fit_rv ):
		plot_rv()

else:
	errores = 'perc'
	print_errors_planets()	


