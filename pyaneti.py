#!/usr/bin/python2.7

#-----------------------------------------------------------
#                         pyaneti.py
#                        DESCRIPTION
#                   Barragan O, March 2016
#-----------------------------------------------------------

#Load libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pyaneti as pti #FORTRAN module

#-------------------------------------------------------------
#                   INITIALIZATION
#-------------------------------------------------------------

#Load the input file
#You have to run the program as ./pyaneti star_name
star = str(sys.argv[1])

#Create path to the input_fit.py file
inf_name = 'inpy/'+star+'/input_fit.py'

#Did you create an input_fit.py file?
if ( not os.path.isfile( inf_name ) ):
  print 'You have not created', inf_name
  sys.exit()

#Read the file with all the python functions
execfile('src/todo-py.py')

#Read the file with the default values
execfile('src/default.py')

#Read input file
execfile(inf_name)

#Prepare data
execfile('src/prepare_data.py')

#Create ouput directory
outdir = 'outpy/' + star + '_out'
if not os.path.exists(outdir):
  os.makedirs(outdir)

#Obtain smart priors based on iput data
smart_priors()

#if is_circular is turned on, we ensure a circular orbit
check_circular()

#Print intial configuration
print_init()

#-------------------------------------------------------------
#                   FITTING ROUTINES
#-------------------------------------------------------------

#Joint fit
if (fit_rv and fit_tr ):
  fit_joint()

#Transit fit
elif ( not fit_rv and fit_tr ):
  fit_transit()

#Radial velocity fit
elif ( fit_rv and not fit_tr ):
  fit_radial_velocity()

#Nothing to fit!
else:
  sys.exit("Nothing to fit!")

#-------------------------------------------------------------
#             	PRINT AND PLOT ROUTINES
#-------------------------------------------------------------

#Print the values
execfile('src/print_values.py')

#Create plots
execfile('src/plot_data.py')

if ( nplanets == 1):

  if (plot_histogram):
    plot_histogram()

  if (plot_correlations):
    plot_correlations()

  #PLOT TRANSIT
  if ( fit_tr ):
    plot_transit()

  #PLOT RV CURVE
  if ( fit_rv ):
    plot_rv_one()
	
elif (nplanets > 1):

  if (plot_histogram):
    hist_mp_rv()

  #plot_correlations()

  #PLOT THE RV curves
  plot_rv_mp()

