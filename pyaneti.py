#!/usr/bin/python2.7

#-----------------------------------------------------------
#                         pyaneti.py
#                        DESCRIPTION
#                   Barragan O, March 2016
#-----------------------------------------------------------

#Load libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm, sigmaclip
import sys
import pyaneti as pti #FORTRAN module

#Read the file with all the python functions
execfile('todo-py.py')

#Read the file with the default values
execfile('default.py')

#Read input file
execfile('input_fit.py')

#Prepare data
execfile('prepare_data.py')

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
#			              FITTING ROUTINES
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
	if ( method == 'mh' ):
		pti.metropolis_hastings(mega_time,mega_rv,mega_err,tlab \
		,megax, megay, megae, params, prec, maxi, thin_factor, \
		is_circular, what_fit, flag, nconv)

	elif ( method == 'sm' ):

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

	elif ( method == 'plot' ):
		print 'I will only print the values and generate the plot'

	else:
		print 'You did not choose a method!'
		print 'method = mh   -> Metropolis-Hasting'
		print 'method = sm   -> Stretch move'
		print 'method = plot -> Plot of a previous run'
		sys.exit('choose your favorite.')

	print 'Reading the data file, wait a bit!'

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
	if ( method == 'mh' ):
		pti.metropolis_hastings_tr(megax, megay, megae,  \
		params, prec, maxi, thin_factor, is_circular, what_fit,flag,nconv)

	elif ( method == 'sm' ):
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

	elif ( method == 'plot' ):
		print 'I will only print the values and generate the plot'

	else:
		print 'You did not choose a method!'
		print 'method = mh   -> Metropolis-Hasting'
		print 'method = sm   -> Stretch move'
		print 'method = plot -> Plot of a previous run'
		sys.exit('choose your favorite.')

	print 'Reading the data file, wait a bit!'

		#Read the data
	vari, chi2,chi2red,t0o,Po,eo,wo,io,ao,u1o,u2o,pzo = \
       np.loadtxt('mh_trfit.dat', comments='#',unpack=True)

#FIT RV CURVE ONLY
elif ( fit_rv and not fit_tr ):

	flag = [is_log_P,is_ew,is_log_k,is_log_rv0]

	if ( is_circular ):
		fit_e = False
		fit_w = False

	#All this should be in the input file
	fit_t0 = [fit_t0,0]
	fit_P  = [fit_P, 0]
	fit_e  = [fit_e,fit_e]
	fit_w  = [fit_w,fit_w]
	fit_k  = [fit_k,fit_k]  
	fit_v0 = [fit_v0,fit_v0]

	T0 = [6813.38345,6817.2759]
	P  = [7.919454,11.90701]
	e  = [e,e]
	w  = [w,w]
	k  = [k0,k0]

	min_t0	= [6812.,6800.]
	max_t0 	= [6820.,6830.]
	min_P	 	= [7.5,7.5]
	max_P	 	= [8.5,15.0]
	min_e		= [1.e-10,1e-10]
	max_e		= [0.8,0.8]
	min_w		= [0.0,0.0]
	max_w		= [2*np.pi,2*np.pi]
	min_k		= [0.00001,0.00001]
	max_k		= [2.0,2.0]
	min_v0	= 1.
	max_v0  = 10.
	#All this should be in the input file


	what_fit = [None]*6*nplanets
	params   = [None]*(5+nt)*nplanets	
	limits   = [None]*(5+nt)*2*nplanets
	#Let us fill the input variables for 
	#all the number of planets
	for m in range(0,nplanets):
		#What to fit from the input lists	
		what_fit[0+6*m] = int(fit_t0[m]) 
		what_fit[1+6*m] = int(fit_P[m]) 
		what_fit[2+6*m] = int(fit_e[m]) 
		what_fit[3+6*m] = int(fit_w[m]) 
		what_fit[4+6*m] = int(fit_k[m]) 
		what_fit[5+6*m] = int(fit_v0[m]) 
		#fill the parameters vector
		params[0+(5+nt)*m] = T0[m]
		params[1+(5+nt)*m] = P[m]
		params[2+(5+nt)*m] = e[m]
		params[3+(5+nt)*m] = w[m]
		params[4+(5+nt)*m] = k[m]
		#fill the systemic velocities
		for j in range(0,nt):
			params[(5+j)+(5+nt)*m] = v0[j]
		#fill the limits
		limits[0+(5+nt)*2*m] = min_t0[m]
		limits[1+(5+nt)*2*m] = max_t0[m]
		limits[2+(5+nt)*2*m] = min_P[m]
		limits[3+(5+nt)*2*m] = max_P[m]
		limits[4+(5+nt)*2*m] = min_e[m]
		limits[5+(5+nt)*2*m] = max_e[m]
		limits[6+(5+nt)*2*m] = min_w[m]
		limits[7+(5+nt)*2*m] = max_w[m]
		limits[8+(5+nt)*2*m] = min_k[m]
		limits[9+(5+nt)*2*m] = max_k[m]
		for j in range(0,nt):
			limits[(10+j*2)+(5+nt)*2*m] = min_v0
			limits[(11+j*2)+(5+nt)*2*m] = max_v0


#	print limits
#	sys.exit('')

	nwalkers = 20 * len(params)
#	nwalkers = 1000
#	nwalkers = 500

	#out_file = 'planet1.dat'
	out_file = ['planet1.dat','planet2.dat']


	if ( method == 'mh' ):
		pti.metropolis_hastings_rv(mega_time,mega_rv,mega_err,tlab,\
  	params, prec, maxi, thin_factor, is_circular, what_fit,flag,nconv)

	elif ( method == 'sm' ):

		pti.stretch_move_rv(mega_time,mega_rv,mega_err,tlab,\
  	params, limits, nwalkers, prec, maxi, thin_factor, \
		is_circular, what_fit,flag,nconv,datas=len(mega_time), \
		nt=nt,npl=nplanets)

	elif ( method == 'plot' ):
		print 'I will only print the values and generate the plot'

	else:
		print 'You did not choose a method!'
		print 'method = mh   -> Metropolis-Hasting'
		print 'method = sm   -> Stretch move'
		print 'method = plot -> Plot of a previous run'
		sys.exit('choose your favorite.')

	print 'Reading the data file, wait a bit!'

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

#Print the values
execfile('print_values.py')

#Create plots
execfile('plot_data.py')

