


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

