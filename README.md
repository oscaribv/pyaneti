# __pyaneti__
#### Written by Barragán O. 2016
##### email: oscaribv@gmail.com
##### Updated Jun 21, 2016

### __Introduction__

* _Pianeti_ is the Italian word for planets.
* It is a _python/fortran_ software suite which finds the best fitting solution  by $\chi^2$ minimization using Marcov Chain Monte Carlo (MCMC) methods.
* Ensemble sampler with affine invariance algorithm
for a major coverage of the parameter space
([Godman & Weare, 2010](http://msp.org/camcos/2010/5-1/p04.xhtml)).
* _Python_ does the nice things: Plots, call to functions, prints, input files.
* _Fortran_ does the hard work: MCMC evolution, $\chi^2$ calculation, ensemble sampler evolution.
* Open-source code (GPL v 3.0). No need to pay!
* ** Free and fast code with the robustness of _Fortran_ and the versatility of _Python_ **.

## __Power of pyaneti__

* Multiple independent Marcov Chains to sample the space parameters.
* Easy-to-use: it runs by providing only one input_fit.py file.
* Automatic creation of histograms, correlations, and ready-to-publish plots.
* Circular and elliptical orbits.
* RV multi-planet fitting.
* Systemic velocities for multiple instruments in RV fittings.
* Stellar limb darkening [(Mandel & Agol, 2002)](http://iopscience.iop.org/article/10.1086/345520/meta#artAbst).
* Correct treatment of short and long cadence data ([Kipping, 2010](http://mnras.oxfordjournals.org/content/408/3/1758)).
* Single and joint RV-transit fits.


## ** Start to use it now! **

Ensure you have installed the gfortran compiler and the python libraries numpy, scipy and matplotlib.

You do not need install anything, just clone or download pyaneti, type make and it is done!

```
make
```

if you have all the packages, you are ready to run pyaneti.

```
./pyaneti.py test
```

or

```
python pyaneti.py test
```

The program will start. Wait a couple of minutes and then you should see something like:

```
Summary:
N_chains    =  100
N_conv      =  500
thin_factor =  1
N_data      =  174
N_pars      =  7
chi2       = 190.2928
DOF         =  167
chi2_red   = 1.1395
scale factor=  0.936800476854
BIC        = 226.4061
Input parameters
M_*     = 1.0000000 + 0.0000099 - 0.0000101 solar masses
R_*     = 1.0000000 + 0.0000100 - 0.0000099 solar radii

The best fit planet parameters are:
T0    = 2448285.2568223 + 0.0001310 - 0.0001186 days
P     = 365.2558537 + 0.0017473 - 0.0019224 days
e     = 0.0000 + 0.0000 - 0.0000     
w     = 90.0000 + 0.0000 - 0.0000 deg
Transit fit parameters:
i     = 89.9646 + 0.0090 - 0.0100 deg
a/r*  = 213.2180 + 0.8937 - 1.1412    
rp/r* = 0.0091 + 0.0000 - 0.0000    
q_1    = 0.0000 + 0.0000 - 0.0000    
q_2    = 0.0000 + 0.0000 - 0.0000    
RV fit parameters:
K     = 0.0847 + 0.0013 - 0.0014 m/s
S v0  = 0.0000 + 0.0000 - 0.0000 km/s

Derived parameters:
r_p   = 0.9960 + 0.0021 - 0.0020 R_earth
a   = 0.9915 + 0.0042 - 0.0053  AU
b r*  = 0.1318 + 0.0364 - 0.0332
t_total = 13.0915 + 0.0069 - 0.0075 hours
t_in/eg = 0.1206 + 0.0014 - 0.0010 hours
rho_* = 1.3745 + 0.0174 - 0.0220 g/cm^3
u_1    = 0.0000 + 0.0000 - 0.0000    
u_2    = 0.0032 + 0.0000 - 0.0000    
Tp    = 2448285.2568 + 0.0001 - 0.0001 days
mp    = 0.9469 + 0.0143 - 0.0154 earth masses
rho_p = 5.2639 + 0.0944 - 0.0928 g/cm^3
g_p = 918.1049 + 17.9905 - 17.4897 cm/s^2
```
Once see this, you will see some plots similar to



## Documentation

_I am working on it!_

## Science  with pyaneti
