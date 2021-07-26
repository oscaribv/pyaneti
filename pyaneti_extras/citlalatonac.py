################################################################################
#                              citlalatonac.py                                 #
#   Citlalatonac is the name of the Aztec god who created the stars            #
#                           Oscar Barragan, March 2021                         #
################################################################################

import numpy as np
import sys
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
import pyaneti as pti
from numpy.random import multivariate_normal
import seaborn as sns
from scipy.interpolate import interp1d
from matplotlib import gridspec
sns.set(style='ticks')
sns.set_color_codes('deep')

#Details on how astroplan works here -> https://astroplan.readthedocs.io/en/latest/tutorials/summer_triangle.html
from astroplan import Observer
from astropy.coordinates import SkyCoord
from astroplan import FixedTarget
from astropy.time import Time

#Brute force function to create times in which the target is observed from a given observatory
#This function create a times vector of lenght ndata in which the star is observed at the observatory
#in the interval between tmin and tmax
def create_real_times(tmin,tmax,ndata=50,air_mass_limit=1.5,tformat='mjd',star='K2-100',observatory='lapalma'):
    times = []
    #Create a observatory instance with astroplan
    observatory = Observer.at_site(observatory)
    #Create star object from astroplan, the name has to be a valid simbad name for the target
    star = FixedTarget.from_name(star)
    while len(times) < ndata:
        #Draw random time between tmin and tmax
        drt = np.random.uniform(tmin,tmax)
        #Get the time object from astropy.time
        time = Time(drt,format=tformat)
        #Compute the air_mass that the target has at time t
        air_mass = observatory.altaz(time,star).secz
        #If the target is observable at time drt and the air_mass < air_mass_limit then
        #We accept the dummy random time
        if observatory.target_is_up(time,star) and air_mass < air_mass_limit:
            times.append(drt)

    #Sort thet times
    times = sorted(times)

    #Return a numpy array with all the times
    return np.array(times)


class citlali():
    """
    Class to create simulated RVs and acitivity/simmetry indicators assuming they all can be described
    by the same Gaussian Process (following Rajpaul et al., 2015, MNRAS, 442, 2269).
    """


    def __init__(self,time=[],tmin=0,tmax=60,nseries=3,\
                 amplitudes=[0.0058,0.0421,0.024,0.0,0.02,-0.86],
                 kernel_parameters=[31.2,0.55,4.315], kernel='QPK',
                 time_series=[],
                 points_per_day=10,seed=123):
        """
        Input parameters:
        tmin        -> minimum temporal value
        tmax        -> maximum temporal value
        amplitudes  -> Amplitudes of the mutli-GP approach, each time-series has to amplitudes following Rajpaul et al., 2015.
        kernel_parameter -> Hyper-parameters for the given kernel
        kernel      -> Chose the kernel,Quasi-Periodic Kernel 'QPK', Matern 5/2 'M52', Exponential 'EXP'
        time_series -> time_series for the symmetry/activity indicators, the RV time-series is labelled by default as 'rvs'
        seed        -> Set a seed for random numbers

        The instance has the following attributes:
        amplitudes  -> Amplitudes that relate the time-series (see Rajpaul et al., 2015)
        time        -> time-stamps between tmin and tmax
        seed        -> Seed used to create the time-series
        gp_rvs      -> GP induced RV time-series
        time_series -> Name of the time-series
        rvs         -> RV time-series
        """

        #Attributes to be stored in the class
        self.kernel = kernel
        self.kernel_parameters = kernel_parameters
        #Multi-GP amplitudes as attribute
        self.amplitudes = amplitudes[:]

        if len(time) > 0:
            self.time = time
        else:
            #Create vector time
            self.time = np.linspace(tmin,tmax,int(points_per_day*(tmax-tmin)))

        #Initiate random seed
        self.seed = seed
        np.random.seed(self.seed)

        #Create timeseries to be used in case they were not given as input
        if len(time_series) == 0: time_series = [ 'a'+str(i) for i in range(1,nseries)]

        #Generate time_series attribute with the name of all the time-series contained in the instance
        self.time_series = np.concatenate([['rvs'],time_series])

        #Create vector with nseries repetions for self.time (required by pyaneti)
        self.bigtime = np.array(np.concatenate([list(self.time)]*nseries))

        #Create vector with hyper parameters as required in covfunc input
        gp_parameters = np.concatenate([amplitudes,kernel_parameters])

        #Get the kernel label that we need to compute the correlation matrix with pyaneti
        mi_kernel = ''
        if self.kernel == 'QPK':
            mi_kernel = 'MQ'+str(nseries)
        elif self.kernel == 'M52':
            mi_kernel = 'MM'+str(nseries)
        elif self.kernel == 'EXP':
            mi_kernel = 'ME'+str(nseries)
        else:
            print("The input kernel is not supported. Supported kernels are: \n")
            print("QPK -> Quasi-Periodic Kernel\n")
            print("M52 -> Matern 5/2 Kernel\n")
            print("EXP -> Exponential Kernel\n")
            sys.exit()

        #Compute the covariance matrix
        cov = pti.covfunc(mi_kernel,gp_parameters,self.bigtime,self.bigtime)

        #Draw a sample
        ceros = np.zeros(len(self.bigtime))
        samples = multivariate_normal(ceros,cov,size=1)
        #samples contains a big vector with the RVs and all the activity indicators

        #Lenght of the time vector
        lt = len(self.time)

        #Separate the big sample in multi-dimensional samples

        #Extract the RVs and put them in the gp_rvs atribute
        self.gp_rvs = np.array(samples[0][0:lt])

        #Generate all the attributes with the time-series
        #each time-series is labeled
        for i,label in enumerate(self.time_series):
            setattr(self, label, np.array(samples[0][i*lt:(i+1)*lt]))


    def add_planet(self,planet_params=[0,0.011,1.67,0,np.pi/2],planet_name='planet_b'):
        """
        This method adds a RV planet signal following a Keplerian orbit where
        T0 = time of minimum conjunction
        K  = Semi-amplitude of the planet induced RV signal
        P  = planetary orbital period
        e  = orbit eccentricity
        w  = angle of periastron
        Input parameters:
        planet_params -> has to be given as T0, K, P, e, w
        planet_name   -> Name of the planet, default is 'planet_b'
        """

        #Create planet_name atrribute with the exoplanet parameters (planet_params)
        setattr(self,planet_name,planet_params)

        #Compute the Keplerian curve with the input planet parameters using pyaneti
        prv = pti.rv_curve_mp(self.time,0.0,*planet_params,0.0,0.0)

        #Set the rv_planet_name attribute with the induced RV signal for the actual planet
        setattr(self,'rv_'+planet_name,prv)

        #Create the planet_names attribute (an empty list) if has not beed previously defined
        if not hasattr(self,'planet_names'):
            self.planet_names = []

        #The following lines will run only if the planet has not beed added previously
        #This avoids to add the same planet more than once, to add more than one planet
        #each planet has to have a different name
        if planet_name not in self.planet_names:
            self.planet_names.append(planet_name)
            self.rvs = self.rvs + getattr(self,'rv_'+planet_name)
        #the planet signal has been added

    def remove_planet(self,planet_name='planet_b'):
        """
        Remove the planet planet_name from the RV signal (self.rvs)
        """

        if planet_name  in self.planet_names:
            #Remove the planet signal
            self.rvs = self.rvs - getattr(self,'rv_'+planet_name)
            #Remove the attributes corresponding to the planet
            delattr(self,'rv_'+planet_name)
            delattr(self,planet_name)
        else:
            #There is no planet to remove
            print('There is no {} to remove'.format(planet_name))


    def create_data(self,t=[],ndata=0):
        """
        This method creates the *_data attributes
        *_data atributes mimic the observations that we want to analyse with pyaneti
        #these attributes can be modified with red and white noise
        There are three ways in which data are created
        if the t vector is provided, then the data points are created at each element of t
        if ndata = X where X is larger than zero, then the code will create X random observations
        if none is specified, then the code assumes that all of the model points are observations
        """

        #First check if we are given an input vector of times to create the data stamps
        if len(t) > 0:
            #If yes, create the time-stamps doing interpolation
            self.ndata = len(t)
            self.time_data = t
            for label in self.time_series:
                f = interp1d(self.time,getattr(self,label), kind='cubic')
                setattr(self,label+'_data',f(self.time_data))
        #We can also say how many ramdom points we want to extract from the sample
        elif (ndata > 0):
            if ndata > len(self.time):
                self.ndata = len(self.time)
            else:
                self.ndata = ndata
            #extract the indices
            #in this option is selected randomly from an uniform distribution
            indices = np.random.random_integers(0,len(self.time),self.ndata)
            #The *_data attributes are filled using indices
            self.time_data = self.time[indices]
            for label in self.time_series:
                setattr(self,label+'_data',getattr(self,label)[indices])
        #If user does not specify, the code assumes that all the sample will be used as data
        else:
            self.ndata = len(self.time)
            self.time_data = self.time[:]
            for label in self.time_series:
                setattr(self,label+'_data',getattr(self,label))


    def add_white_noise(self,err=[]):
        """
        Add white noise to the simulated data
        You need to run get_data_from_sample before
        #err can be a list in which element is rather,
        a float indicating the error in each time_series or
        a list with individual errors for each time_series
        """

        #First check if the data attributes have been created, if not, create them
        if not hasattr(self,'time_data'):
            print("You have not created the *_data atributes yet!")
            return

        #Avoids the code to crash in case the input is not a list
        if err.__class__ == float:
            err = [err]

        #If the list is not as requested, the values are set randomly
        if len(err) == 0 or len(err) < len(self.time_series):
            #If user does not specify the error, the error is computed randomly
            err = np.random.uniform(1e-5,1e-3,len(self.time_series))

        #Add white noise to each time-series
        for i, label in enumerate(self.time_series):
            #Create the label+'_err'
            setattr(self,label+'_err',err[i])
            #Modify the label+'_data' by adding the white noise
            setattr(self,label+'_data',np.random.normal(getattr(self,label+'_data'),getattr(self,label+'_err')))


    def add_red_noise(self,se_parameters=[1e-4,1]):
        """
        Add red noise to each time-series using a square-exponential kernel
        It allows for a different red noise set of hyperparameters per each time-series
        """

        #First check if the data attributes have been created, if not, create them
        if not hasattr(self,'time_data'):
            print("You have not created the *_data atributes yet!")
            return

        kernel_parameters = []
        #Check if the user is providing enough number of square exponential parameters per each time-series
        if len(se_parameters) < len(self.time_series):
            #We add the same hyperparameters for all the time-series
            #this is not recommended!
            for i in self.time_series:
                kernel_parameters.append(se_parameters)
        else:
            kernel_parameters = se_parameters[:]
        #Now kernel_parameters is a list where each element is a two-element list with hyperparameters of a square-exponential kernel

        #Add red noise to each time-series
        for i,label in enumerate(self.time_series):
            #Let's start by adding the same level of noise to all time-series
            #Add red noise by adding a sample of a quasi periodic kernel
            cov = pti.covfunc('SEK',kernel_parameters[i],self.time_data,self.time_data)
            #Draw a sample
            ceros = np.zeros(len(self.time_data))
            sample = multivariate_normal(ceros,cov)
            #Add the sample to the corresponding time-series
            setattr(self, label+'_data', getattr(self,label+'_data') + sample)


    def periodogram_rvs(self,freqs=np.linspace(0.01,10,100)):
        from scipy.signal import lombscargle
        self.pgram = lombscargle(self.time_data,self.rvs_data,freqs,normalize=True)
        plt.ylabel('Power')
        plt.xlabel('Period')
        plt.plot(1/freqs,self.pgram)
        plt.show()

    def plot(self,fsx=10,fsy=10./4.,save=False,fname='timeseries.pdf',show=True):

        nseries = len(self.time_series)

        plt.figure(figsize=(fsx,nseries*fsy))
        gs = gridspec.GridSpec(nrows=nseries,ncols=1)
        gs.update(hspace=0.002)
        plt.subplot(gs[0])
        plt.plot(self.time,self.rvs,'-',label='RV')
        if hasattr(self,'planet_names'):
            for planet in self.planet_names:
                plt.plot(self.time,getattr(self,'rv_'+planet),'-',label=planet)
        if hasattr(self,'time_data'):
                plt.plot(self.time_data,self.rvs_data,'o',label='simulated data')
        plt.ylabel('RV [km/s]')
        plt.legend(ncol=1,scatterpoints=1,numpoints=1,frameon=True)

        try:
            for i,label in enumerate(self.time_series[1:]):
                plt.subplot(gs[i+1])
                plt.ylabel(label)
                plt.plot(self.time,getattr(self,label),'-',label=label)
                try:
                    plt.plot(self.time_data,getattr(self,label+'_data'),'o',label='simulated data')
                except:
                    pass
                plt.legend(ncol=1,scatterpoints=1,numpoints=1,frameon=True)
        except:
            pass

        plt.xlabel('Time [days]')

        if save: plt.savefig(fname,format='pdf',bbox_inches='tight')

        if show: plt.show()


    def save_data(self,fname='multigp_data.dat'):
        """
        Save the data in the format that pyaneti needs to run
        Time    Value    Value_error   label
        """

        #First check if the data attributes have been created, if not, create them
        if not hasattr(self,'time_data'):
            print("You have not created the *_data atributes yet!")
            return

        with open(fname,'w') as f:
            f.write('# {} \n'.format(fname))
            if hasattr(self,'planet_names'): f.write('# Number of planets = {} \n'.format(len(self.planet_names)))
            f.write('# Number of simulated observations = {} \n'.format(self.ndata))
            f.write('# Kernel used : {}, with parameters: {} \n'.format(self.kernel,self.kernel_parameters))
            f.write('# Random number seed = {} \n'.format(self.seed))
            f.write('# Time   Value   Value_error  time_series_label \n')
            for label in self.time_series:
                for i in range(len(self.time_data)):
                    f.write("{:1.7e}  {:1.7e}  {:1.7e}  {}\n".format(self.time_data[i],getattr(self,label+'_data')[i], \
                                                    getattr(self,label+'_err'),label))
