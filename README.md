<p align="center">
  <img width = "500" src="./src/images/logo_pyaneti.png"/>
</p>

# __pyaneti__*
*From the Italian word _pianeti_, which means planets

#### email:  oscar.barragan_at_physics.ox.ac.uk
#### Updated November 2021

## Paper I
### Written by Barragán O., Gandolfi D. & Antoniciello G.

<a href="https://academic.oup.com/mnras/advance-article/doi/10.1093/mnras/sty2472/5094600"><img src="https://img.shields.io/badge/MNRAS-2019,482,1017-blueviolet.svg" alt="MNRAS" /></a>
<a href="https://arxiv.org/abs/1809.04609"><img src="https://img.shields.io/badge/arXiv-1809.04609-green.svg" alt="arXiv:1809.04609" /></a>
<a href="http://ascl.net/1707.003"><img src="https://img.shields.io/badge/ascl-1707.003-green.svg" alt="ascl:1707.003" /></a>
<a href="https://github.com/oscaribv/pyaneti/wiki"><img src="https://img.shields.io/badge/wiki-building-yellow.svg" alt="pyaneti wiki" /></a>

## Paper II
### Written by Barragán O., Aigrain S., Rajpaul V. M., & Zicher N.

<a href="https://academic.oup.com/mnras/advance-article-abstract/doi/10.1093/mnras/stab2889/6383008"><img src="https://img.shields.io/badge/MNRAS-2022, 509, 866-blueviolet.svg" alt="MNRAS" /></a>
<a href="https://arxiv.org/abs/2109.14086"><img src="https://img.shields.io/badge/arXiv-2109.140860-green.svg" alt="arXiv:2109.140860" /></a>


#### Brief description on pyaneti:

* The code runs in _python 3_.
* Transit fits for single transits.
* Multi-band fits.
* Gaussian Process (GP) and multidimensional GP regressions.
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
* Joint RV + transit fitting.

If you want to see the cool stuff that pyaneti can do check
[these papers](https://ui.adsabs.harvard.edu/search/fq=database%3A%20astronomy&fq=%7B!type%3Daqp%20v%3D%24fq_database%7D&fq_database=database%3A%20astronomy&p_=0&q=%20full%3A%22pyaneti%22&sort=date%20desc%2C%20bibcode%20desc).

#### Check pyaneti [wiki](https://github.com/oscaribv/pyaneti/wiki) to learn how to use it

### Citing

If you use pyaneti in your research, please cite it as

```
Barragán, O., Gandolfi, D., & Antoniciello, G., 2019, MNRAS, 482, 1017
```

you can use this bibTeX entry

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

If you also use the new routines of pyaneti (multi-band modelling, single transit modelling, or Gaussian Process regression), please cite also this paper

```
Barragán,  O.,  Aigrain,  S.,  Rajpaul,  V.  M.,  &  Zicher,  N.,  2022, MNRAS. 509, 866
```

you can use this bibTeX entry

```
@ARTICLE{pyaneti2,
       author = {{Barrag{\'a}n}, Oscar and {Aigrain}, Suzanne and {Rajpaul}, Vinesh M. and {Zicher}, Norbert},
        title = "{PYANETI - II. A multidimensional Gaussian process approach to analysing spectroscopic time-series}",
      journal = {\mnras},
     keywords = {methods: numerical, techniques: photometry, techniques: spectroscopy, planets and satellites: general, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2022,
        month = jan,
       volume = {509},
       number = {1},
        pages = {866-883},
          doi = {10.1093/mnras/stab2889},
archivePrefix = {arXiv},
       eprint = {2109.14086},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..866B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```


**If you have any comments, requests, suggestions or just need any help, please don't think twice, just contact us!**

#### Warning: This code is under developement and it may contain bugs. If you find something please contact us at oscar.barragan_at_physics.ox.ac.uk

## Acknowledgements
* Hannu Parviainen, thank you for helping us to interpret the first result of the PDF of the MCMC chains. We learned a lot from you!
* Salvador Curiel, thank you for  suggestions to parallelize the code.
* Mabel Valerdi, thank you for being the first _pyaneti_ user, for spotting typos and errors in this document. And thank you much for the awesome idea for pyaneti's logo.
* Lauren Flor, thank you for testing the code before release.
* Jorge Prieto-Arranz, thank you for all the suggestions which have helped to improve the code.

**THANKS A LOT!**

### For a beta version of pyaneti check [here](https://github.com/oscaribv/pyaneti-dev)
