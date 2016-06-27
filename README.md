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

## ** Dependencies **

You need to have in your computer:
* gfortran
* numpy
* matplotlib

## ** Start to use it now! **

You do not need install anything, just clone or download pyaneti.

```
git clone https://github.com/oscaribv/pyaneti
```

The advantage about cloning the repository is the possibility to follow the changes to this package easily with git pull (learn more about git
in [https://git-scm.com/](https://git-scm.com/)).
Or if you want

```
wget https://github.com/oscaribv/pyaneti/archive/master.zip
unzip master.zip
mv pyaneti_master pyaneti
```

if you choose this option, you should repeat it every time the code is updated.

The next step is to get inside the pyaneti folder and see what we can find inside it

```
cd pyaneti
ls
  LICENSE   outpy  src inpy  makefile  pyaneti.py  README.md
```

 type ``make``, you would see a lot of text apearing in the terminal, after it finishes, you are done!

```
make
```

If you have all the dependencies installed, the make procces ended without any error.
Now you are ready to run the code for the first time.

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
N_chains    =  50
N_conv      =  500
thin_factor =  1
N_data      =  174
N_pars      =  7
chi2       = 168.4746
DOF         =  167
chi2_red   = 1.0088
scale factor=  0.995614042083
BIC        = 204.5880
Input parameters
M_*     = 1.0005287 + 0.0498459 - 0.0501476 solar masses
R_*     = 0.9996468 + 0.0502219 - 0.0496047 solar radii

The best fit planet parameters are:
T0    = 2448285.0956823 + 0.0017752 - 0.0019077 days
P     = 365.2533004 + 0.0014328 - 0.0014945 days
e     = 0.0000 + 0.0000 - 0.0000     
w     = 90.0000 + 0.0000 - 0.0000 deg
Transit fit parameters:
i     = 89.9999 + 0.0000 - 0.0000 deg
a/r*  = 213.9276 + 11.9416 - 9.1588    
rp/r* = 0.0090 + 0.0000 - 0.0000    
q_1    = 0.5510 + 0.0821 - 0.0759    
q_2    = 0.3672 + 0.0792 - 0.0945    
RV fit parameters:
K     = 0.0896 + 0.0007 - 0.0007 m/s
S v0  = 22.0720 + 0.0000 - 0.0000 km/s

Derived parameters:
r_p   = 0.9821 + 0.0498 - 0.0491 R_earth
a   = 0.9976 + 0.0732 - 0.0676  AU
b r*  = 0.0003 + 0.0000 - 0.0000
t_total = 13.1609 + 0.5886 - 0.6962 hours
t_in/eg = 0.1175 + 0.0053 - 0.0062 hours
rho_* = 1.3882 + 0.2457 - 0.1708 g/cm^3
u_1    = 0.5471 + 0.0727 - 0.1145    
u_2    = 0.1969 + 0.1628 - 0.1210    
Tp    = 2448285.0957 + 0.0018 - 0.0019 days
mp    = 1.0017 + 0.0347 - 0.0352 earth masses
rho_p = 5.8022 + 0.9945 - 0.8179 g/cm^3
g_p = 1006.6511 + 115.7284 - 85.6492 cm/s^2
```
Once see this, you will see some plots similar to

![Transit first fit](/outpy/test_out/transit_test.png)
![Radial-Velocity first fit](/outpy/test_out/rv_test.png)

Let me explain you briefly what this is:  
> If you were an advanced alien civilization with really high technology, and "lucky" enough to see the Earth crossing in front of the Sun, **this is how the Earth would look like to you**.

Look at those well known parameters:    
* 1 Earth Mass
* 1 Earth radii
* Period of 365 days
* 1 AU semi-major axis
* Density of 5.5-6.0 g/cm^2,
* Gravity of 10 m/s^2.  

Of course you would need a spectograph with a precision of a few cm/s and also a very nice photometer.

> If you are at this point, you learned two things. First, with good data you can obtain really nice planet parameters and second, you learned how to run pyaneti.  


## ** Create your own setup **


## Documentation

_I am working on it!_

## Science  with pyaneti
