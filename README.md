

<p align="center">
  <img width = "500" src="./src/images/logo_pyaneti.png"/>
</p>

# __pyaneti__
#### Written by Barragán O. et al. 2017
##### email: oscaribv@gmail.com
##### Updated April 19, 2017

### __Introduction__

* _Pianeti_ is the Italian word for planets.
* It is a _python/fortran_ software suite which finds the best fitting solution using Marcov Chain Monte Carlo (MCMC) methods.
* Ensemble sampler with affine invariance algorithm
for a major coverage of the parameter space
([Godman & Weare, 2010](http://msp.org/camcos/2010/5-1/p04.xhtml)).
* _Python_ does the nice things: Plots, call to functions, prints, input files.
* _Fortran_ does the hard work: MCMC evolution, Likelihood calculation, ensemble sampler evolution.
* Open-source code (GPL v 3.0). No need to pay!
* **Free and fast code with the robustness of _Fortran_ and the versatility of _Python_**.

## __Power of pyaneti__

* Multiple independent Marcov Chains to sample the space parameters.
* Easy-to-use: it runs by providing only one input_fit.py file.
* Parallel computing with OpenMP.
* Automatic creation of histograms, correlations, and ready-to-publish plots.
* Circular and elliptical orbits.
* Multi-planet fitting.
* Systemic velocities for multiple instruments.
* Stellar limb darkening [(Mandel & Agol, 2002)](http://iopscience.iop.org/article/10.1086/345520/meta#artAbst).
* Correct treatment of short and long cadence data ([Kipping, 2010](http://mnras.oxfordjournals.org/content/408/3/1758)).
* Single and joint RV and transit fits.

## Dependencies

You need to install in your computer:
* gfortran
* OpenMP
* numpy
* matplotlib
* seaborn (optional)

## Start to use it now!

Just clone or download pyaneti.

```
git clone https://github.com/oscaribv/pyaneti
```

The advantage about cloning the repository is the possibility to follow the changes to this package easily with git pull (learn more about git
in [https://git-scm.com/](https://git-scm.com/)).

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

The program will start. You will see a lot of things appearing on your screen, ignore them now. Wait a couple of minutes and then you should see something like:

```
--------------------------------------------------------------
Summary:
N_chains         =      100
N_conv           =      500
thin_factor      =        1
N_rv_data        =       36
N_tr_data        =      108
N_data           =      144
N_pars           =        8
chi2_rv          = 7.2743
chi2_tr          = 121.4930
chi2             = 128.7674
dof              =      136
chi2/dof         = 0.9468
ln likelihood_rv = 377.7464
ln likelihood_tr = 1090.7782
ln likelihood    = 1468.5246
BIC              = -2897.2906
--------------------------------------------------------------
             INPUT STELLAR PARAMETERS
--------------------------------------------------------------
M_*     = 1.0000000 - 0.1000000 + 0.1000000 solar masses
R_*     = 1.0000000 - 0.1000000 + 0.1000000 solar radii
T_*     = 5600.0000000 - 100.0000000 + 100.0000000 K
--------------------------------------------------------------
                   Parameters test b
-------------------------Fitted-------------------------------
T0   = 2448285.0952767 - 0.0030782 + 0.0032988  days
P    = 365.2509089 - 0.0044430 + 0.0046877  days
ew 1 = 0.0000000 - 0.0000000 + 0.0000000
ew 2 = 0.0000000 - 0.0000000 + 0.0000000
b    = 0.0000000 - 0.0000000 + 0.0000000
a/R* = 215.5628668 - 2.1954734 + 2.2406256
Rp/R*= 0.0093297 - 0.0000787 + 0.0000751
K    = 0.0880342 - 0.0023864 + 0.0022963  m/s
-------------------------Derived------------------------------
e    = 0.0000000 - 0.0000000 + 0.0000000
w*   = 0.0000000 - 0.0000000 + 0.0000000  deg
i    = 90.0000000 - 0.0000000 + 0.0000000  deg
a    = 1.0026256 - 0.1019944 + 0.1003731  AU
rho* = 1.4203407 - 0.0429408 + 0.0447797  g/cm^3 (transit light curve)
rho* = 1.4057819 - 0.3710566 + 0.5574300  g/cm^3 (input stellar parameters)
Mp   = 0.9825596 - 0.0710232 + 0.0711111 M_earth
Rp   = 1.0175344 - 0.1033253 + 0.1011605 R_earth
rho_p= 5.1140865 - 1.3040095 + 1.9831153  g/cm^3
g_p  = 935.2212496 - 31.5878768 + 33.8424277  cm/s^2
Tperi= 2448193.7825041 - 0.0038257 + 0.0041448  days
Teq  = 269.7042101 - 4.9940052 + 5.0331832  K (albedo=0)
T_tot= 13.0651372 - 0.1342224 + 0.1340438  hours
T_i/e= 0.1207591 - 0.0013183 + 0.0012596  hours
--------------------------------------------------------------
--------------------  Other parameters -----------------------
q1    = 0.3727471 - 0.0410572 + 0.0434470
q2    = 0.3208062 - 0.0434040 + 0.0488237
u1    = 0.3918669 - 0.0568551 + 0.0597566
u2    = 0.2177398 - 0.0608878 + 0.0574922
Super_telescope = 22.0719863 - 0.0000017 + 0.0000018 km/s
--------------------------------------------------------------
```
If you see this output it means that pyaneti ended succesfully!

Now let us check the plots.

```
evince testb_rv.pdf testb_tr.pdf
```

You will see something like this


<img src="./src/images/testb_tr.png" style="width: 250px;"/>
<img src="./src/images/testb_rv.png" style="width: 250px;"/>

Let me explain you briefly what this is:
> If you were an advanced alien civilization with really high technology, and "lucky" enough to see the Earth crossing in front of the Sun, **this is how the Earth would look like to you**.

Look at those well known parameters:
* 1 Earth Mass
* 1 Earth radii
* Period of 365 days
* 1 AU semi-major axis
* Density of 5.5 g/cm^2,
* Gravity of 10 m/s^2.

Of course you would need a spectograph with a precision of a few cm/s and also a very nice photometer.

> If you are at this point, you learned two things. First, with good data you can obtain really nice planet parameters and second, you learned how to run pyaneti.


## Documentation

There is no a manual by now. But there are some examples on how to run different
systems. The cases are:
* test - Joint radial velocity and light curve fit.
* kepler-10 - two planets radial velocity fit.

## Science  with pyaneti
* Barragán et al., 2017, _EPIC 218916923 b: a low-mass warm Jupiter on a 29-day orbit transiting an active K0V star_,
MNRAS, submitted (https://arxiv.org/abs/1702.00691).
* Nespral et al., 2017, _Mass determination of K2-19b and K2-19c from radial velocities and transit timing variations_,
A&A, in press (http://arxiv.org/abs/1604.01265).
* Nowak et al., 2017, _EPIC 219388192 b - an inhabitant of the brown dwarf desert in the Ruprecht 147 open cluster_,
 AJ, 153, 131 (https://arxiv.org/abs/1610.08571).
* Eigmüller et al., 2017, _K2-60b and K2-107b. A sub-jovian and a jovian planet from the k2 mission_,
AJ, 153, 130 (https://arxiv.org/abs/1611.03704).
* Barragán et al, 2016, _K2-98b: A 32-M⊕ Neptune-sized planet in a 10-day orbit transiting an F8 star_,
 AJ, 152, 6 (http://arxiv.org/abs/1608.01165).

## Acknowledgements
* to Davide Gandolfi, for his support and motivation to create this code.
* to Hannu Parviainen, to help me to interpret the first result of the PDF of the MCMC chains. I learned a lot!
* to Mabel Valerdi, to help me with as the first pyaneti user, to detect bugs and errors in this manual. As well for her idea for pyaneti's logo.
* to Salvador Curiel, for his suggestions to parallelize the code.

**THANKS!**
