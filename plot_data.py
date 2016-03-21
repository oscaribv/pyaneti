	#Let us do the plots here

#one planet case
if ( nplanets == 1 ):

	def plot_transit():
		  #Move all the points to T0
	  for i in range(0,ntr):
	    xt[i] = xt[i] - P_val * i

	  #Redefine megax with the new xt values
	  megax = np.concatenate(xt)
	  z_val = pti.find_z(megax,t0_val,P_val,e_val,w_val\
		  ,i_val,a_val)
	  mud_val, mu0_val = pti.occultquad(z_val,u1_val,u2_val\
		  ,pz_val)
	  #Residuals
	  res = megay - mud_val

	  #Get the model data to do the plot
	  nvec = int(1e5)
	  dx = ( max(megax) - min(megax) ) / nvec
	  xvec = np.zeros(nvec)
	  xvec[0] = min(megax)
	  for i in range(1,nvec):
	    xvec[i] = xvec[i-1] + dx
	  zvec = pti.find_z(xvec,t0_val,P_val,e_val,w_val,i_val,a_val)
	  mud, mu0 = pti.occultquad(zvec,u1_val,u2_val,pz_val)
	  #Now we have data to plot a nice model

	  #Do the plot
	  plt.figure(2,figsize=(10,10))
	  #Plot the transit light curve
	  plt.subplot(211)
	  plt.xlim(min(xt[0]),max(xt[0]))
	  plt.errorbar(megax,megay,megae,fmt='o',alpha=0.3)
	  plt.plot(xvec,mud,'k',linewidth=2.0)
	  #Plot the residuals
	  plt.subplot(212)
	  plt.xlim(min(xt[0]),max(xt[0]))
	  plt.errorbar(megax,res,megae,fmt='o',alpha=0.3)
	  plt.plot(megax,np.zeros(len(megax)),'k--',linewidth=2.0)
	  plt.show()


	#Plot RV
	def plot_rv():
	  #Create the RV fitted curve
	  n = 5000
	  xmin = t0_val
	  xmax = t0_val + P_val
	  dn = (xmax - xmin) /  n
	  rvx = np.empty([n])
	  rvx[0] = xmin
	  for i in range(1,n):
	    rvx[i] = rvx[i-1] + dn
	  if ( is_circular ):
	    rvy = pti.rv_circular(rvx,0.0,t0_val,k_val,P_val)
	  else:
	    rvy = pti.rv_curve(rvx,0.0,t0_val,k_val,P_val,e_val,w_val)

	  res = [None]*nt
	  for i in range(0,nt):
	   if (is_circular):
	      res[i] = pti.rv_circular(time_all[i],0.0,t0_val,\
		       k_val,P_val)
	   else:
	      res[i] = pti.rv_curve(time_all[i],0.0,t0_val,k_val,\
		       P_val,e_val,w_val)
	   rv_all[i] = rv_all[i] - v_val[i]
	   res[i] = rv_all[i] - res[i]

	  p_rv = scale_period(rvx,t0_val,P_val)
	  p_all = [None]*nt
	  #tp_val = pti.find_tp(t0_val,e_val,w_val,P_val)
	  for i in range(0,nt):
	    p_all[i] = scale_period(time_all[i],t0_val,P_val)

	  plt.figure(3,figsize=(10,10))
	  plt.subplot(311)
	  plt.xlabel("Phase")
	  plt.ylabel(ylab)
	  plt.ylim(-1.4*k_val,1.4*k_val)
	  plt.plot(p_rv,rvy,'k',label=('k=%2.2f m/s'%k_val ))
	  mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*']
	  for i in range(0,nt):
	    plt.errorbar(p_all[i],rv_all[i],errs_all[i],\
	    label=telescopes[i],fmt=mark[i],alpha=0.6)
	  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
		     ncol=4, mode="expand", borderaxespad=0.)
	  plt.subplot(312)
	  plt.xlabel("Phase")
	  plt.ylabel(ylab)
	  plt.plot([0.,1.],[0.,0.],'k--')
	  for i in range(0,nt):
	    plt.errorbar(p_all[i],res[i],errs_all[i],\
	    label=telescopes[i],fmt=mark[i],alpha=0.6)
	  plt.show()

	#PLOT TRANSIT
	  if ( fit_tr ):
	  	plot_transit()

  #PLOT RV CURVE
	  if ( fit_rv ):
	  	plot_rv()

else:
	#Plot RV for multiplanet
	def plot_rv_mp():
		rvy = [None]*nplanets
		p_rv = [None]*nplanets
		k_dum = [None]*nplanets
		for i in range(0,nplanets):
			k_dum[i] = 1e3*k_val[i]
		for i in range(0,nt):
			v_val[i] = 1e3*v_val[i]
			for j in range(0,len(rv_all[i])):
				rv_all[i][j] = 1e3*rv_all[i][j]
				errs_all[i][j] = 1e3*errs_all[i][j]


		for i in range(0,nplanets):
			rv_dum = []
			for j in range(0,nt):
				rv_dum.append(rv_all[j])
			#Create the RV fitted model for the planet i
			n = 5000
			xmin = t0_val[i]
			xmax = t0_val[i] + P_val[i]
			dn = (xmax - xmin) /  n
			rvx = np.empty([n])
			rvx[0] = xmin
			for j in range(1,n):
				rvx[j] = rvx[j-1] + dn
			rvy[i] = pti.rv_curve_mp(rvx,0.0,t0_val[i],\
						k_dum[i],P_val[i],e_val[i],w_val[i])

			dt0_val = []
			dk_dum = []
			dP_val = []
			de_val = []
			dw_val = []

			j = 0
			while ( j < nplanets ):
				if ( j != i ):
					dt0_val.append(t0_val[j])
					dk_dum.append(k_dum[j])
					dP_val.append(P_val[j])
					de_val.append(e_val[j])
					dw_val.append(w_val[j])
				j = j + 1


			res = [None]*nt
			drvy = [None]*nt
			for j in range(0,nt):
				#This is the model of the actual planet
				res[j] = pti.rv_curve_mp(time_all[j],0.0,t0_val[i],k_dum[i],\
	      	         P_val[i],e_val[i],w_val[i])
				#This variable has all the others planets 
				drvy[j] = pti.rv_curve_mp(time_all[j],0.0,dt0_val,dk_dum \
									,dP_val,de_val,dw_val)
				#the actual value, minus the systemic velocity, minus the other planets
				rv_dum[j] = rv_dum[j] - v_val[j] - drvy[j]
				res[j] = rv_dum[j] - res[j]

			p_rv[i] = scale_period(rvx,t0_val[i],P_val[i])
			p_all = [None]*nt
			for j in range(0,nt):
				p_all[j] = scale_period(time_all[j],t0_val[i],P_val[i])

			plt.figure(3,figsize=(10,10))
			plt.subplot(311)
			plt.xlabel("")
			plt.ylabel(ylab)
			#plt.ylim(-1.4*k_dum[i],1.4*k_dum[i])
			plt.plot(p_rv[i],rvy[i],'k',label=('k=%2.2f m/s'%k_dum[i] ))
			mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*']
			for j in range(0,nt):
				plt.errorbar(p_all[j],rv_dum[j],errs_all[j],\
				label=telescopes[j],fmt=mark[j],alpha=0.6)
				plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
		     ncol=4, mode="expand", borderaxespad=0.)
			plt.subplot(312)
			plt.xlabel("Phase")
			plt.ylabel('Residuals (m/s)')
			plt.plot([0.,1.],[0.,0.],'k--')
			for j in range(0,nt):
				plt.errorbar(p_all[j],res[j],errs_all[j],\
				label=telescopes[j],fmt=mark[j],alpha=0.6)
			fname = 'planet' + str(i+1) + '.pdf'
			plt.savefig(fname)
			plt.show()

	#call the function
	plot_rv_mp()
