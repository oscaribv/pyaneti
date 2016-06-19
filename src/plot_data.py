#Let us do the plots here

from matplotlib import gridspec
from matplotlib.colors import LogNorm


def plot_rv_fancy(p_rv,rvy,p_all,rv_dum,errs_all,res,telescopes_labels,fname):
  plt.figure(3,figsize=(7,6))
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3., 1.])
  gs.update(hspace=0.00) 
  ax0 = plt.subplot(gs[0])
  plt.minorticks_on()
  #plt.subplot(311)
  ax0 = plt.xlabel("")
  ax0 = plt.ylabel("RV (m/s)")
  ax0 = plt.plot([0.,1.],[0.,0.],'k--')
  ax0 = plt.plot(p_rv,rvy,'k',linewidth=1.0)
  mark = ['o', 'd', '^', '<', '>', '8', 's', 'p', '*']
  for j in range(0,nt):
    ax0 = plt.errorbar(p_all[j],rv_dum[j],errs_all[j],\
    label=telescopes_labels[j],fmt=mark[j],alpha=0.8)

  #ax0 = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4, ncol=nt, mode="expand", borderaxespad=0.)
  plt.legend(loc=0, ncol=1,scatterpoints=1,numpoints=1,frameon=False,fontsize='small')

  plt.xticks(np.arange(0.,1.01,0.1)) 
  plt.tick_params( axis='x',which='both',labelbottom='off') 
  #plt.subplot(312)
  ax1 = plt.subplot(gs[1])
  plt.xlabel("Orbital phase")
  plt.tick_params( axis='x',which='minor',bottom='on',left='on',right='on',top='on') 
  plt.xticks(np.arange(0.,1.01,0.1)) 
  plt.ylabel('Residuals (m/s)')
  plt.plot([0.,1.],[0.,0.],'k--',linewidth=1.0)
  for j in range(0,nt):
    plt.errorbar(p_all[j],res[j],errs_all[j],\
    label=telescopes_labels[j],fmt=mark[j],alpha=0.8)
  yylims = ax1.get_ylim()
  plt.yticks(np.arange(yylims[0],yylims[1],(yylims[1]-yylims[0])/4.))
  plt.minorticks_on()
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.show()


#===========================================================
#                   One planet plots
#===========================================================

if ( nplanets == 1 ):

	#=========================#
	#       Transit plot      #
	#=========================#
	def plot_transit():
		  #Move all the points to T0
	  for i in range(0,ntr):
	    xt[i] = xt[i] - P_val * i 

	  #Redefine megax with the new xt values
	  megax = np.concatenate(xt)
          flag = [False, False, False, False]
	  z_val = pti.find_z(megax,[t0_val,P_val,e_val,w_val
		  ,i_val,a_val],flag)
	  mud_val, mu0_val = pti.occultquad(z_val,q1_val,q2_val\
		  ,pz_val)
	  #Residuals
	  res = megay - mud_val

          #Let us plot the binned model 

	  #Redefine megax with the new xt values
	  megax = np.concatenate(xt)
	  smegax = sorted(megax)
          flag = [False, False, False, False]

          xd_ub = np.ndarray(shape=(n_cad,len(smegax)))
          xd_ub_res = np.ndarray(shape=(n_cad,len(smegax)))
          zd_ub = [None]*n_cad
          zd_ub_res = [None]*n_cad
          fd_ub = [None]*n_cad
          fd_ub_res = [None]*n_cad
          for m in range(0,n_cad):
            for n in range(0,len(smegax)):
                xd_ub[m][n] = smegax[n] + t_cad * ( (m+1) - 0.5 * (n_cad + 1 )  ) / n_cad
                xd_ub_res[m][n] = megax[n] + t_cad * ( (m+1) - 0.5 * (n_cad + 1 )  ) / n_cad

          for m in range(0,n_cad):
             zd_ub[m] = pti.find_z(xd_ub[m][:],[t0_val,P_val,e_val,w_val,i_val,a_val],flag)
	     fd_ub[m], dummm = pti.occultquad(zd_ub[m],q1_val,q2_val,pz_val)
             zd_ub_res[m] = pti.find_z(xd_ub_res[m][:],[t0_val,P_val,e_val,w_val,i_val,a_val],flag)
	     fd_ub_res[m], dummm = pti.occultquad(zd_ub_res[m],q1_val,q2_val,pz_val)

	  fd_reb = [0.0]*len(megax)
	  fd_reb_res = [0.0]*len(megax)
          for m in range(0,len(megax)):
	    for n in range(0,n_cad):
             fd_reb[m] = fd_reb[m] + fd_ub[n][m]/n_cad
             fd_reb_res[m] = fd_reb_res[m] + fd_ub_res[n][m]/n_cad

          res_res = megay - fd_reb_res

	  #Get the model data to do the plot
	  nvec = int(1e5)
	  dx = ( max(megax) - min(megax) ) / nvec
	  xvec = np.zeros(nvec)
	  xvec[0] = min(megax)
	  for i in range(1,nvec):
	    xvec[i] = xvec[i-1] + dx
	  zvec = pti.find_z(xvec,[t0_val,P_val,e_val,w_val,i_val,a_val],flag)
	  mud, mu0 = pti.occultquad(zvec,q1_val,q2_val,pz_val)
	  #Now we have data to plot a nice model

	  #Do the plot
          tfc = 24. # time factor conversion to hours
	  plt.figure(1,figsize=(7,6))
	  #Plot the transit light curve
	  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3.0, 1.])
          gs.update(hspace=0.00) 
	  plt.subplot(gs[0])
	  x_lim = (min(xt[0])-T0)*tfc
	  plt.xlim(x_lim,-x_lim)
          min_val_model = max(fd_reb) -  min(fd_reb)
	  plt.errorbar((megax-T0)*tfc,megay,megae,fmt='r.',alpha=0.8)
	  plt.plot((smegax-T0)*tfc,fd_reb,'k',linewidth=1.0)
          plt.ylabel('Relative flux')
          plt.xticks( np.arange(int(x_lim),int(-x_lim)+1,1))
          plt.minorticks_on()
          plt.tick_params( axis='x',which='both',labelbottom='off') 
	  #Plot the residuals
	  dplot = plt.subplot(gs[1])
	  plt.plot((smegax-T0)*tfc,np.zeros(len(smegax)),'k--',linewidth=1.0)
	  plt.errorbar((megax-T0)*tfc,res_res,megae,fmt='r.',alpha=0.8)
          yylims = dplot.get_ylim()
          plt.yticks(np.arange(yylims[0],yylims[1],(yylims[1]-yylims[0])/5.))
          plt.xticks( np.arange(int(x_lim),int(-x_lim)+1,1))
	  plt.xlim(x_lim,-x_lim)
	  #plt.errorbar(megax,res,megae,fmt='o',alpha=0.8)
	  #Plot the residuals
          plt.ylabel('Residuals')
          plt.xlabel("T - T0 (hours)")
          plt.minorticks_on()
	  plt.savefig(outdir+'/'+star+plabels[0]+'_tr.pdf',format='pdf',bbox_inches='tight')
	  plt.show()


	#=========================#
	#        RV plot          #
	#=========================#

  #Plot RV for one planet
	def plot_rv_one():
		k_dum = 1e3*k_val
		for i in range(0,nt):
			v_val[i] = 1e3*v_val[i]
			for j in range(0,len(rv_all[i])):
				rv_all[i][j] = 1e3*rv_all[i][j]
				errs_all[i][j] = 1e3*errs_all[i][j]

		rv_dum = []
		for j in range(0,nt):
			rv_dum.append(rv_all[j])
		n = 5000
		xmin = t0_val
		xmax = t0_val + P_val
		dn = (xmax - xmin) /  n
		rvx = np.empty([n])
		rvx[0] = xmin
		for j in range(1,n):
			rvx[j] = rvx[j-1] + dn

		rvy = pti.rv_curve_mp(rvx,0.0,t0_val,\
						         k_dum,P_val,e_val,w_val)

		res = [None]*nt
		for j in range(0,nt):
			#This is the model of the actual planet
			res[j] = pti.rv_curve_mp(time_all[j],0.0,t0_val,k_dum,\
	      	         P_val,e_val,w_val)
			#the actual value, minus the systemic velocity
			rv_dum[j] = rv_dum[j] - v_val[j] 
			res[j] = rv_dum[j] - res[j]

			p_rv = scale_period(rvx,t0_val,P_val)
			p_all = [None]*nt
			for j in range(0,nt):
				p_all[j] = scale_period(time_all[j],t0_val,P_val)

		fname = outdir+'/'+star+plabels[0]+'_rv.pdf'
                plot_rv_fancy(p_rv,rvy,p_all,rv_dum,errs_all,res,telescopes_labels,fname)


#===========================================================
#                   Multi-planet plots
#===========================================================

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

                        fname = outdir+'/'+star+plabels[i]+'_rv.pdf'
                        plot_rv_fancy(p_rv[i],rvy[i],p_all,rv_dum,errs_all,res,telescopes_labels,fname)

#===========================================================
#                   Histogram plots
#===========================================================

def create_plot_histogram(params,plabs,cbars='red',nb=50):
	n = len(params)
	plt.figure(1,figsize=(15,4*(n)/2))
	gs = gridspec.GridSpec(nrows=(n+1)/2,ncols=2)
	for i in range(0,n):
		plt.subplot(gs[i])
		vpar, lpar, rpar = find_vals_perc(params[i],nconv,1.0)
                minchi2_val = params[i][minchi2_index]
		plt.axvline(x=vpar,c=cbars)
		plt.axvline(x=minchi2_val,c='yellow')
		plt.axvline(x=vpar-lpar,c=cbars,ls='--')
		plt.axvline(x=vpar+rpar,c=cbars,ls='--')
		plt.xlabel(plabs[i])
		plt.hist(params[i],normed=True,bins=nb)

	plt.savefig(outdir+'/'+star+plabels[0]+'_histogram.pdf',format='pdf',bbox_inches='tight')
	plt.show()

def plot_histogram(rf=1):

	if ( fit_tr and fit_rv ):
		dparams = [ chi2red[0::rf],t0o[0::rf],Po[0::rf],eo[0::rf],wo[0::rf],io[0::rf],ao[0::rf],q1o[0::rf],q2o[0::rf],pzo[0::rf],ko[0::rf]]
		dplabs = ['$\chi^2_{red}$','$T0$','$P$','$e$','$\omega$','$i$','$a/R_*$','$q_1$','$q_2$','$R_p/R_*$','$k$']

		vlabs = [None]*nt
		dvo = [None]*nt
		for i in range(0,nt):
			vlabs[i] = 'rv0 ' + telescopes_labels[i]
			dvo[i] = vo[i][0::rf]


		params = np.concatenate([dparams,dvo])
		labs = np.concatenate([dplabs,vlabs])
	
		create_plot_histogram(params,labs)


	if ( fit_tr and not fit_rv ):
		params = [t0o[0::rf],Po[0::rf],eo[0::rf],wo[0::rf],io[0::rf],ao[0::rf],q1o[0::rf],q2o[0::rf],pzo[0::rf]]
		labs = ['$T0$','$P$','$e$','$\omega$','$i$','$a/R_*$','$q_1$','$q_2$','$R_p/R_*$']
	
		create_plot_histogram(params,labs)

	if ( not fit_tr and fit_rv ):
		if (nplanets == 1 ):
			dparams = [t0o[0::rf],Po[0::rf],eo[0::rf],wo[0::rf],ko[0::rf]]
			dplabs = ['$T0$','$P$','$e$','$\omega$','$k$']
		else:
			dparams = [None]*5*nplanets
			dplabs = [None]*5*nplanets
			for i in range(0,nplanets):
				dparams[0+5*nplanets] = t0o[i][0::rf]
				dparams[1+5*nplanets] = Po[i][0::rf]
				dparams[2+5*nplanets] = eo[i][0::rf]
				dparams[3+5*nplanets] = wo[i][0::rf]
				dparams[4+5*nplanets] = ko[i][0::rf]
				dplabs[0+5*nplanets] = 'T0'+str(i)
				dplabs[1+5*nplanets] = 'P'+str(i)
				dplabs[2+5*nplanets] = 'e'+str(i)
				dplabs[3+5*nplanets] = '$\omega$'+str(i)
				dplabs[4+5*nplanets] = 'k'+str(i)

		vlabs = [None]*nt
		dvo = [None]*nt
		for i in range(0,nt):
			vlabs[i] = 'rv0 ' + telescopes_labels[i]
			dvo[i] = vo[i][0::rf]


		params = np.concatenate([dparams,dvo])
		labs = np.concatenate([dplabs,vlabs])
	
		create_plot_histogram(params,labs)


def hist_mp_rv(cbars='red',nb=50):

	for l in range(0,nplanets):
		plt.figure(1,figsize=(10,3*(5+nt)/2))
		gs = gridspec.GridSpec(nrows=int((5+nt)/2.+0.5),ncols=2)
		ax0 = plt.subplot(gs[0])
		ax0 = plt.axvline(x=t0_val[l],c=cbars)
		ax0 = plt.axvline(x=t0_val[l]-t0_errl[l],c=cbars,ls='--')
		ax0 = plt.axvline(x=t0_val[l]+t0_errr[l],c=cbars,ls='--')
		ax0 = plt.xlabel('T0')
		ax0 = plt.axvline(x=t0_val[l])
		ax0 = plt.hist(t0o[l],normed=True,bins=nb)
		ax1 = plt.subplot(gs[1])
		ax1 = plt.axvline(x=P_val[l],c=cbars)
		ax1 = plt.axvline(x=P_val[l]-P_errl[l],c=cbars,ls='--')
		ax1 = plt.axvline(x=P_val[l]+P_errr[l],c=cbars,ls='--')
		ax1 = plt.xlabel('P')
		ax1 = plt.hist(Po[l],normed=True,bins=nb)
		ax2 = plt.subplot(gs[2])
		ax2 = plt.axvline(x=e_val[l],c=cbars)
		ax2 = plt.axvline(x=e_val[l]-e_errl[l],c=cbars,ls='--')
		ax2 = plt.axvline(x=e_val[l]+e_errr[l],c=cbars,ls='--')
		ax2 = plt.xlabel('$e$')
		ax2 = plt.hist(eo[l],normed=True,bins=nb)
		ax3 = plt.subplot(gs[3])
		ax3 = plt.axvline(x=w_val[l],c=cbars)
		ax3 = plt.axvline(x=w_val[l]-w_errl[l],c=cbars,ls='--')
		ax3 = plt.axvline(x=w_val[l]+w_errr[l],c=cbars,ls='--')
		ax3 = plt.xlabel('$\omega$')
		ax3 = plt.hist(wo[l],normed=True,bins=nb)
		ax4 = plt.subplot(gs[4])
		ax4 = plt.axvline(x=k_val[l],c=cbars)
		ax4 = plt.axvline(x=k_val[l]-k_errl[l],c=cbars,ls='--')
		ax4 = plt.axvline(x=k_val[l]+k_errr[l],c=cbars,ls='--')
		ax4 = plt.xlabel('k')
		ax4 = plt.hist(ko[l],normed=True,bins=nb)
		for m in range(0,nt):
			plt.subplot(gs[5+m])
			plt.axvline(x=v_val[m],c=cbars)
			plt.axvline(x=v_val[m]-v_errl[m],c=cbars,ls='--')
			plt.axvline(x=v_val[m]+v_errr[m],c=cbars,ls='--')
			plt.xlabel('%s rv0'%(telescopes_labels[m]))
			plt.hist(vo[m],normed=True,bins=nb)
		plt.savefig(outdir+'/'+star+'_hist_params'+str(l)+'.pdf',format='pdf',bbox_inches='tight')
		plt.show()

#===========================================================
#                   Correlation plots
#===========================================================

def create_plot_correlation(params,plabs,col='red',mark='.'):
	n = len(params)
	plt.figure(1,figsize=(2*n,2*n))
	gs = gridspec.GridSpec(nrows=n,ncols=n)
	for i in range(0,n):
		for j in range(0,i):
			plt.subplot(gs[i*n+j])
			if ( j == 0 ):
				plt.ylabel(plabs[i])
			elif ( j == i - 1 ):
				#plt.colorbar() 
				plt.tick_params( axis='y',which='both',labelleft='off') 
			else:
				plt.tick_params( axis='y',which='both',labelleft='off') 
			if ( i == n-1):
				plt.xlabel(plabs[j])
			else:
				plt.tick_params( axis='x',which='both',labelbottom='off') 
			#plt.plot(params[j],params[i],c=col,marker=mark,ls='',alpha=0.5)

                        plt.hist2d(params[j],params[i],bins=100,norm=LogNorm())
	plt.savefig(outdir+'/'+star+plabels[0]+'_correlations.pdf',format='pdf',bbox_inches='tight')
	plt.show()

def plot_correlations(rf=1):

	if ( fit_tr and fit_rv ):
		dparams = [chi2red[0::rf],t0o[0::rf],Po[0::rf],eo[0::rf],wo[0::rf],io[0::rf],ao[0::rf],q1o[0::rf],q2o[0::rf],pzo[0::rf],ko[0::rf]]
		dplabs = ['$\chi^2_{red}$','$T0$','$P$','$e$','$\omega$','$i$','$a/R_*$','$q_1$','$q_2$','$R_p/R_*$','$k$']

		vlabs = [None]*nt
		dvo = [None]*nt
		for i in range(0,nt):
			vlabs[i] = 'rv0 ' + telescopes_labels[i]
			dvo[i] = vo[i][0::rf]

		params = np.concatenate([dparams,dvo])
		labs = np.concatenate([dplabs,vlabs])
	
		create_plot_correlation(params,labs,col='blue')


	if ( fit_tr and not fit_rv ):

		params = [t0o[1::rf],Po[1::rf],eo[1::rf],wo[1::rf],io[1::rf],ao[1::rf],q1o[1::rf],q2o[1::rf],pzo[1::rf]]
		labs = ['$T0$','$P$','$e$','$\omega$','$i$','$a/R_*$','$q_1$','$q_2$','$R_p/R_*$']
	
		create_plot_correlation(params,labs,col='blue')

	#Now it works only for RV fit	
	if ( not fit_tr and fit_rv ):
		if ( nplanets == 1 ):
			dparams = [t0o[1::rf],Po[1::rf],eo[1::rf],wo[1::rf],ko[1::rf]]
			dplabs = ['$T0$','$P$','$e$','$\omega$','$k$']
		else:
			dparams = [None]*5*nplanets
			dplabs = [None]*5*nplanets
			for i in range(0,nplanets):
				dparams[0+5*nplanets] = t0o[i][1::rf]
				dparams[1+5*nplanets] = Po[i][1::rf]
				dparams[2+5*nplanets] = eo[i][1::rf]
				dparams[3+5*nplanets] = wo[i][1::rf]
				dparams[4+5*nplanets] = ko[i][1::rf]
				dplabs[0+5*nplanets] = 'T0'
				dplabs[1+5*nplanets] = 'P'
				dplabs[2+5*nplanets] = 'e'
				dplabs[3+5*nplanets] = '$\omega$'
				dplabs[4+5*nplanets] = 'k'

		vlabs = [None]*nt
		dvo = [None]*nt
		for i in range(0,nt):
			vlabs[i] = 'rv0 ' + telescopes_labels[i]
			dvo[i] = vo[i][1::rf]


		params = np.concatenate([dparams,dvo])
		labs = np.concatenate([dplabs,vlabs])
	
		create_plot_correlation(params,labs,col='blue')

