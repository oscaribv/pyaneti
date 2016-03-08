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
	pti.metropolis_hastings(mega_time,mega_rv,mega_err,tlab \
	,megax, megay, megae, params, prec, maxi, thin_factor, \
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
	pti.metropolis_hastings_tr(megax, megay, megae,  \
	params, prec, maxi, thin_factor, is_circular, what_fit,flag,nconv)

	#Read the data
	vari, chi2,chi2red,t0o,Po,eo,wo,io,ao,u1o,u2o,pzo = \
       np.loadtxt('mh_trfit.dat', comments='#',unpack=True)

#FIT RV CURVE ONLY
elif ( fit_rv and not fit_tr ):

	flag = [is_log_P,is_ew,is_log_k,is_log_rv0]

	if ( is_circular ):
		fit_e = False
		fit_w = False

	what_fit = [int(fit_t0),int(fit_P),int(fit_e),int(fit_w), \
					    int(fit_k), int(fit_v0)]
	dummy = [T0,P,e,w,k0]
	params = np.concatenate((dummy,v0))
	
	#Call fit routine
	pti.metropolis_hastings_rv(mega_time,mega_rv,mega_err,tlab,\
  params, prec, maxi, thin_factor, is_circular, what_fit,flag,nconv)

	#Read the data
	vari,chi2,chi2red,t0o,Po,eo,wo,ko = \
	np.loadtxt('mh_rvfit.dat', comments='#', unpack=True,\
	usecols=range(0,8))
	vo = [None]*nt
	for j in range(0,nt):
		n = [8+j]
		a = np.loadtxt('mh_rvfit.dat', comments='#',unpack=True, usecols=(n))
		vo[j] = a

#Nothing to fit!
else:
	sys.exit("Nothing to fit!")

#-------------------------------------------------------------
#			END FITTING ROUTINES
#-------------------------------------------------------------

#Find the values and errors
t0_val,t0_err = find_vals_gauss(t0o,nconv)
if (is_log_P):
	Po = np.power(10.,Po)
P_val, P_err  = find_vals_gauss(Po,nconv)
if (is_ew):
	dummy_e = eo
	eo = eo * eo + wo * wo
	wo = np.arctan2(dummy_e,wo)
e_val,e_err 	= find_vals_gauss(eo,nconv)
w_val,w_err 	= find_vals_gauss(wo,nconv)
if (w_val < 0.0 ):
	w_val = w_val + 2 * np.pi	
w_deg 		= w_val * 180. / np.pi
w_deg_err = w_err * 180. / np.pi
if ( fit_tr ):
	if (is_sini):
		io = np.arcsin(io)
	i_val,i_err 	= find_vals_gauss(io,nconv)
	if (is_log_a):
		ao = np.power(10.,ao)
	a_val,a_err 	= find_vals_gauss(ao,nconv)
	u1_val,u1_err = find_vals_gauss(u1o,nconv)
	u2_val,u2_err = find_vals_gauss(u2o,nconv)
	pz_val,pz_err = find_vals_gauss(pzo,nconv)
	i_deg 		= i_val * 180. / np.pi
	i_deg_err = i_err * 180. / np.pi
if ( fit_rv ):
	if ( is_log_k ):
		ko = np.power(10.,ko)
	k_val, k_err  = find_vals_gauss(ko,nconv)
	v_val = [None]*nt
	v_err = [None]*nt
	for j in range(0,nt):
		if ( is_log_rv0 ):
			vo[j] = np.power(10.,vo[j])
		v_val[j], v_err[j] = find_vals_gauss(vo[j],nconv)


#Print the best fit values values
print ('The best fit planet parameters are:')
print ('T0    = %4.4e +/- %4.4e days'%(t0_val,t0_err))
print ('P     = %4.4e +/- %4.4e' 		%(P_val,P_err))
print ('e     = %4.4e +/- %4.4e'			%(e_val,e_err))
print ('w     = %4.4e +/- %4.4e deg'	%(w_deg,w_deg_err))
if (fit_tr):
	print ('i     = %4.4e +/- %4.4e deg' %(i_def,i_deg_err))
	print ('a/r*  = %4.4e +/- %4.4e' 		%(a_val,a_err))
	print ('u1    = %4.4e +/- %4.4e' 		%(u1_val,u1_err))
	print ('u2    = %4.4e +/- %4.4e' 		%(u2_val,u2_err))
	print ('rp/r* = %4.4e +/- %4.4e' 		%(pz_val,pz_err))
if (fit_rv):
	print ('K    = %4.4e +/- %4.4e' 		%(k_val,k_err))
	for i in range(0,nt):
		print ('%s v0 = %4.4e +/- %4.4e' 	%(telescopes[i], \
		v_val[i],v_err[i]))


#PLOT TRANSIT
if ( fit_tr ):
	plot_transit()

#PLOT RV CURVE
if ( fit_rv ):
	plot_rv()
