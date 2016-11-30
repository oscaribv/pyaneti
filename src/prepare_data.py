#-----------------------------------------------------------
#                       prepare_data.py
#  This file contains all the variable initializations,
#  both for RV and Transit fittings.
#                   O. Barragan, March 2016
#-----------------------------------------------------------

#Create the stellar data
#mstar = np.random.normal(loc=mstar_mean,scale=mstar_sigma,size=nwalkers*nconv)
#rstar = np.random.normal(loc=rstar_mean,scale=rstar_sigma,size=nwalkers*nconv)

#What transit data are we fitting
if ( lc_data == 'kepler_lc' ):
  n_cad = 10
  t_cad = 29.425 / 60. / 24.0 #days
elif ( lc_data == 'kepler_sc' ):
  n_cad = 1
  t_cad = 1.5 / 60. / 24.0 #days
elif ( lc_data == 'free' ):
  #values given by the user
  n_cad = n_cad
  t_cad = t_cad

#-----------------------------------------------------------
#                         RV DATA
#-----------------------------------------------------------

#Let us check the kind of variable
nplanets_rv = 0
for o in range(0,len(fit_rv)):
  nplanets_rv = nplanets_rv + int(fit_rv[o])

#Let us ensure that we want to fit rv data
if ( nplanets_rv > 0 ):

	#Read the data file
	#time, RV, errors, and Telescope label
	time,rv,err,tspe = np.loadtxt('inpy/'+star+'/'+fname_rv[0],usecols=(0,1,2,3), \
  	dtype={'names': ('time', 'rv', 'err','telescope'), \
		'formats': ('float', 'float', 'float', 'S10')}, \
		comments='#',unpack=True)

	#These lists have lists with data for
	#the different telescopes
	time_all = []
	rv_all   = []
	errs_all = []

	#Number of telescopes
	nt = len(telescopes)

	if ( nt < 1 ):
		print 'Please, indicate the telescope labels!'
		sys.exit('')

	#Separate the data for each telescope and create labels
	for i in range(0,nt):
		time_dum = []
		rv_dum   = []
		errs_dum = []
		for j in range(0,len(tspe)):
			if (tspe[j][0] == telescopes[i][0]):
				time_dum.append(time[j])
				rv_dum.append(rv[j])
				errs_dum.append(err[j])
	#The *all variables are lists of lists, each list constains
	# a list with the data of each telescope
		time_all.append(time_dum)
		rv_all.append(rv_dum)
		errs_all.append(errs_dum)

	#The mega* variables contains all the data 
	#All this is neccesary because you do not have
	#the same number of data for each telescope
	mega_rv   = []
	mega_time = []
	mega_err  = []
	tlab      = []
	#create mega with data of all telescopes
	for i in range(0,nt): 
		#fill the mega variable with all the data of the
		#telescope i
		for j in range(0,len(rv_all[i])):
			#tlab has the label of the telescope (an integer)
			#this is useful because matches with the index of 
			#the mega variables
			tlab.append(i)
			mega_rv.append(rv_all[i][j])
			mega_time.append(time_all[i][j])
			mega_err.append(errs_all[i][j])

	total_rv_fit = True

else:
  tlab = [0]
  mega_rv = [None]
  mega_time = [None]
  mega_err  = [None]
  total_rv_fit = False

#RV DATA READY

#-----------------------------------------------------------
#                     TRANSIT DATA
#-----------------------------------------------------------

#Let us check the kind of variable
nplanets_tr = 0
for o in range(0,len(fit_tr)):
  nplanets_tr = nplanets_tr + int(fit_tr[o])

if ( nplanets_tr > 0 ):

  #Each transit planet hasa different file
  xt= [None]*nplanets
  yt= [None]*nplanets
  et= [None]*nplanets
  for o in range(0,nplanets):

    filename = 'inpy/'+star+'/'+fname_tr[o]
    dummyd,dummyf,dummye = np.loadtxt(filename,usecols=columns_tr, \
    comments='#',unpack=True)

    dummyd = dummyd + textra

    hdate = dummyd
    wflux = dummyf
    errs  = dummye

    if ( my_tr_ranges == True ):
      if (ntr < 2 or len(tls) < 2 ):
        print 'You selected my_tr_ranges = True\n'
        print 'Please, define ntr and ranges'
          #Get the transit ranges
    else:      #This assumes that the input file has the different transits separated
      tls, ntr = get_transit_ranges(hdate,gap_between_transits[o])
    #print tls
    #sys.exit()


    #crash if you do not have more than one transit
    if ( ntr < 2):
      print "you do not have enought transit data!"
      sys.exit("I crashed because I want more data!")

    #Each element of these lists will have the information
    #of a given transit
    xt[o]= [None]*ntr
    yt[o]= [None]*ntr
    et[o]= [None]*ntr

    #Normalize all the transit independently
    #the transit data is inside the limits tls
    for i in range(0,ntr):
      xt[o][i],yt[o][i],et[o][i] = separate_transits(hdate,wflux,errs,tls[i])


  #Let us put together the information of all the arrays
  #the mega* lists have the data of all the transits
  #in 1D array
  megax = []
  megay = []
  megae = []
  megap = []
  for i in range(0,nplanets):
      for j in range(0,len(xt[i])):
        for k in range(0,len(xt[i][j])):
            megap.append(i)
            megax.append(xt[i][j][k])
            megay.append(yt[i][j][k])
            megae.append(et[i][j][k])

  total_tr_fit = True

else:
  megax = [1.]
  megay = [1.]
  megae = [1.]
  megap = [0]
  total_tr_fit = False

#TRANSIT DATA READY

#CHECK WHAT WE HAVE TO FIT
#If we are not going to fit RV or TR data, let us turn off the variables
#for the given case
for o in range(0,nplanets):
    if (fit_tr[o] == False ):
      fit_pz[o] = False
      pz[o] = 1e-10
      fit_i[o] = False
      ii[o] = 1e-10
      fit_a[o] = False
      a[o] = 1e-10
    if (fit_rv[o] == False ):
      fit_k[o] = False
      k0[o] = 1e-10

#Let us turn off velocity offset for a pure TR fit
if ( not total_rv_fit ):
  fit_v0 = False
  nt = 1
  min_rv0 = [-1.]
  max_rv0 = [1.]
  min_phys_rv0 = [-1.]
  max_phys_rv0 = [1.]
  rvs = [0.0]
  telescopes = ['O']
  telescopes_labels = ['']

