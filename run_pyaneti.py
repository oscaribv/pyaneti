# -----------------------------------------------------------
#                         pyaneti.py
#                     Main pyaneti file
#                   Barragan O, March 2016
# -----------------------------------------------------------

# Load libraries
from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import pyaneti as pti  # FORTRAN module

# -------------------------------------------------------------
#                   INITIALIZATION
# -------------------------------------------------------------

# Load the input file
# You have to run the program as ./pyaneti star_name
star = str(sys.argv[1])

# Create path to the input_fit.py file
inf_name = 'inpy/'+star+'/input_fit.py'

# Did you create an input_fit.py file?
if (not os.path.isfile(inf_name)):
    print('You have not created', inf_name)
    sys.exit()

# Read the file with all the python functions
exec(open('src/todo-py.py').read())

# Read the file with the default values
exec(open('src/default.py').read())

# Read input file
exec(open(inf_name).read())

# Prepare data
exec(open('src/prepare_data.py').read())

# Create ouput directory
outdir = outdir + star + '_out'
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Obtain smart priors based on iput data
if is_smart_priors:
    smart_priors()

print_init()

# -------------------------------------------------------------
#                   FITTING ROUTINES
# -------------------------------------------------------------

joint_fit()

# -------------------------------------------------------------
#             	PRINT AND PLOT ROUTINES
# -------------------------------------------------------------

exec(open('src/output.py').read())

# -------------------------------------------------------------
#             	 END pyaneti.py FILE
# -------------------------------------------------------------
