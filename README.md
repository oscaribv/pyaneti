

<p align="center">
  <img width = "500" src="./src/images/logo_pyaneti.png"/>
</p>

# __pyaneti__

<a href="https://academic.oup.com/mnras/advance-article/doi/10.1093/mnras/sty2472/5094600"><img src="https://img.shields.io/badge/MNRAS-2019,482,1017-purple.svg" alt="MNRAS" /></a>
<a href="https://arxiv.org/abs/1809.04609"><img src="https://img.shields.io/badge/arXiv-1809.04609-green.svg" alt="arXiv:1809.04609" /></a>
<a href="http://ascl.net/1707.003"><img src="https://img.shields.io/badge/ascl-1707.003-green.svg" alt="ascl:1707.003" /></a>
<a href="https://github.com/oscaribv/pyaneti/wiki"><img src="https://img.shields.io/badge/wiki-building-yellow.svg" alt="pyaneti wiki" /></a>

### Written by Barragán O., Gandolfi D. & Antoniciello G.
#### email: oscaribv@gmail.com
#### Updated December 09, 2018


# The code

* _Pianeti_ is the Italian word for planets.
* Multi-planet fitting of radial velocity and transit data.
* It uses Markov chain Monte Carlo (MCMC) methods with a Bayesian approach.
* Ensemble sampler with affine invariance algorithm
([Godman & Weare, 2010](http://msp.org/camcos/2010/5-1/p04.xhtml)).
* _Python_ does the nice things: plots, call functions, printing, in/output files.
* _Fortran_ does the hard work: MCMC evolution, likelihood calculation, ensemble sampler evolution.
* Open-source code (GPL v 3.0).

**Free and fast code with the robustness of _Fortran_ and the versatility of _Python_**.

## Power of pyaneti

* Multiple independent Markov chains to sample the parameter space.
* Easy-to-use: it runs by providing only one input_fit.py file.
* Parallel computing with OpenMP.
* Automatic creation of posteriors, correlations, and ready-to-publish plots.
* Circular and eccentric orbits.
* Multi-planet fitting.
* Inclusion of RV and photometry jitter.
* Systemic velocities for multiple instruments.
* Stellar limb darkening [(Mandel & Agol, 2002)](http://iopscience.iop.org/article/10.1086/345520/meta#artAbst).
* Correct treatment of short and long cadence data ([Kipping, 2010](http://mnras.oxfordjournals.org/content/408/3/1758)).
* Single joint RV + transit fitting.


## Learn how to install and use pyaneti [here](https://github.com/oscaribv/pyaneti/wiki/Start-to-use-pyaneti-now!)


### Check pyaneti wiki to learn how to use it:

* How to run pyaneti in parallel [here](https://github.com/oscaribv/pyaneti/wiki/Parallel-run)
* How to fit a single RV signal with 51 Peg b [here](https://github.com/oscaribv/pyaneti/wiki/RV-fit-for-a-single-planet)

## Citing

If you use pyaneti in your research, please cite it as

```
Barragán, O., Gandolfi, D., & Antoniciello, G., 2019, MNRAS, 482, 1017
```

you can use the bibTeX entry

```
@ARTICLE{pyaneti,
       author = {Barrag\'an, O. and Gandolfi, D. and Antoniciello, G.},
        title = "{PYANETI: a fast and powerful software suite for multiplanet radial
        velocity and transit fitting}",
      journal = {\mnras},
     keywords = {methods: numerical, techniques: photometric, techniques: spectroscopic,
        planets and satellites: general, Astrophysics - Earth and
        Planetary Astrophysics, Astrophysics - Instrumentation and
        Methods for Astrophysics, Physics - Data Analysis, Statistics
        and Probability},
         year = 2019,
        month = Jan,
       volume = {482},
        pages = {1017-1030},
          doi = {10.1093/mnras/sty2472},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/#abs/2019MNRAS.482.1017B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## What will come next?

* Gaussian process.
* TTV.
* Multiband transit photometry fitting.
* Graphical User Interface.


**If you have any comments, requests, suggestions or just need any help, please don't think twice, just contact us!**

#### Warning: This code is under developement and it may contain bugs. If you find something please contact us at oscaribv@gmail.com

## Acknowledgements
* Hannu Parviainen, thank you for helping us to interpret the first result of the PDF of the MCMC chains. We learned a lot from you!
* Salvador Curiel, thank you for  suggestions to parallelize the code.
* Mabel Valerdi, thank you for being the first _pyaneti_ user, for spotting typos and errors in this document. And thank you much for the awesome idea for pyaneti's logo.
* Lauren Flor, thank you for testing the code before release.
* Jorge Prieto-Arranz, thank you for all the suggestions which have helped to improve the code.

**THANKS A LOT!**

