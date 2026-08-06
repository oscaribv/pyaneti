# -----------------------------------------------------------
#                    todo-py.py
# This file contains a lot of useful of python functions.
#	    Oscar Barragan, March, 2016
# -----------------------------------------------------------

# Useful libraries
import numpy as np
import matplotlib.pyplot as plt
import sys


# -----------------------------------------------------------
# planet_mass -> gives the planet mass from RV parameters
# input: mstar -> mass of the orbited star (solar masses)
# k    -> semi-amplitude of the RV (m/s)
# P    -> planet orbital period (days)
# ecc  -> orbit eccentricity (no unit)
# i    -> orbit inclination (radians)
# output: mp -> planet mass (solar masses)
# WARNING: the defalut value for i is pi/2 (90 degrees),
# if you do not know the orbit inclination, this function
# computes the planet mass times sin(i)
# -----------------------------------------------------------

def planet_mass(mstar, k, P, ecc, i):
    # Broadcast/scalar handling: work with 1d arrays internally
    mstar_a = np.atleast_1d(mstar).astype(float)
    k_a     = np.atleast_1d(k).astype(float)
    P_a     = np.atleast_1d(P).astype(float)
    ecc_a   = np.atleast_1d(ecc).astype(float)
    inc_a   = np.atleast_1d(i).astype(float)   # avoid shadowing name i

    # Convert P to seconds (work on a copy)
    P_s = P_a * 24.0 * 3600.0  # s

    # small-planet analytic approximation used as initial guess
    unoe = np.sqrt(1.0 - ecc_a * ecc_a)
    # note: S_GM_SI ~ G * M_sun (SI); mstar is in solar masses so keep consistency:
    # mstar**(2/3.) here expects solar-mass units scaled with S_GM_SI usage below,
    # we will solve in *solar mass units* consistently to preserve original spirit:
    # compute m_p*sin i (in solar mass units) approx:
    # convert factor: (2*pi*S_GM_SI / P)^{1/3} is SI; to keep units consistent we will use SI -> solar mass later
    # To keep minimal changes, follow original algebra but be careful in Newton below.
    # Compute an initial m_p*sin i in SI-kg units:
    # mpsin_si = k * (2*pi*S_GM_SI/P_s)**(-1/3) * (mstar_a * M_SUN)**(2/3) * unoe
    mpsin_si = k_a * (2.0 * np.pi * S_GM_SI / P_s)**(-1.0/3.0) * (mstar_a**(2.0/3.0)) * unoe

    # convert mpsin_si (which currently uses S_GM_SI with mstar in solar units) into solar-mass units:
    # since S_GM_SI = G * M_sun, the algebra above returns mpsin in units of solar masses already
    # (this matches the original code expectation). So we can treat mpsin as solar-mass units:
    mpsin = mpsin_si

    # initial mp guess (solar mass units) : divide by sin(i) unless sini==0
    sini = np.sin(inc_a)
    with np.errstate(divide='ignore', invalid='ignore'):
        mp = mpsin / sini

    # prepare constants for exact mass-function solve:
    # mass function RHS (solar-mass units): mf = (P * K^3 * (1-e^2)^{3/2}) / (2*pi*S_GM_SI)
    # Note: with S_GM_SI defined as G*M_sun and mstar in solar units the units are consistent with mp in solar units.
    mf = (P_s * k_a**3 * unoe**3) / (2.0 * np.pi * S_GM_SI)

    # elementwise Newton-Raphson to solve:
    #   F(m) = (m^3 * sin^3 i) / (mstar + m)^2 - mf = 0
    # derivative:
    #   dF/dm = sin^3 i * m^2 * (3*mstar + m) / (mstar + m)^3
    tol = 1e-8
    max_iter = 200

    # ensure mp array has same shape and is writable
    mp = np.array(mp, dtype=float, copy=True)
    mstar_a = np.array(mstar_a, dtype=float, copy=False)

    # loop over elements (small arrays expected in typical usage)
    for idx in range(mp.size):
        # local values (solar-mass units)
        mstar_i = mstar_a.flat[idx]
        mp_i = mp.flat[idx]
        mf_i = mf.flat[idx]
        s_i = sini.flat[idx]
        s3 = s_i**3

        # invalid mf or zero sine => cannot solve, leave mp as NaN
        if not np.isfinite(mf_i) or mf_i <= 0 or (s_i == 0.0):
            mp.flat[idx] = np.nan
            continue

        # Newton iterations
        converged = False
        for it in range(max_iter):
            denom = (mstar_i + mp_i)
            if denom <= 0.0:
                mp_i = np.nan
                break

            F = (s3 * mp_i**3) / (denom**2) - mf_i
            # stable derivative
            dF = s3 * mp_i**2 * (3.0 * mstar_i + mp_i) / (denom**3)

            # safeguard small/zero derivative
            if not np.isfinite(dF) or dF == 0.0:
                break

            delta = F / dF
            mp_new = mp_i - delta

            # damp if step leads to negative mass or NaN
            if not np.isfinite(mp_new) or mp_new <= 0.0:
                mp_new = mp_i - 0.5 * delta

            # check convergence (relative)
            if np.abs(mp_new - mp_i) < tol * max(1.0, mp_i):
                mp_i = mp_new
                converged = True
                break

            mp_i = mp_new

        mp.flat[idx] = mp_i

        # if not converged and mp_i is nan or not finite, keep initial approx (mpsin/sin i)
        if (not converged) or (not np.isfinite(mp_i)):
            # fallback: use small-planet approx (mpsin/sin i)
            # ensure we don't divide by zero (s_i already checked)
            mp.flat[idx] = (mpsin.flat[idx] / s_i) if s_i != 0.0 else np.nan

    # Return shape handling: if original inputs were scalars, return scalar
    if np.ndim(mstar) == 0 and np.ndim(k) == 0 and np.ndim(P) == 0 and np.ndim(ecc) == 0 and np.ndim(i) == 0:
        return float(mp[0])
    return mp

def planet_mass_old(mstar, k, P, ecc,i):


    mstar = np.array(mstar)
    k = np.array(k)
    P = np.array(P)
    ecc = np.array(ecc)
    i = np.array(i)

    # Convert the orbital period P to seconds
    P = P * 24. * 3600.  # s

    # Initial estimate of mass (assuming mstar >> mp)
    unoe = np.sqrt(1. - ecc * ecc)
    mpsin = k * (2. * np.pi * S_GM_SI / P)**(-1./3.) * mstar**(2./3.) * unoe
    mp = mpsin / np.sin(i)

    # Newton-Raphson algorithm to solve the mass function
    cte = - unoe**3 * P * k**3 / (2. * np.pi * S_GM_SI)
    sini = np.sin(i)

    # Initialize flag and convergence condition for each element in the array
    tolerance = 1e-8     # Tolerance for convergence
    # Start Newton-Raphson iteration
    for i in range(len(P)):
        flag = True
        while flag:
            # Compute f and df for Newton-Raphson step
            f = cte[i] + (mp[i] * sini[i])**3 / (mstar[i] + mp[i])**2
            df = mp[i]**2 * sini[i]**3 / (mstar[i] + mp[i])**2 * (3. - 2. * mp[i] / (mstar[i] + mp[i]))

            # Update the mass using the Newton-Raphson step
            mp[i] = mp[i] - f / df

            # Check for convergence
            flag = np.abs(f) > tolerance

    #Remove all the nans
    mp = mp[~np.isnan(mp)]

    return mp


# -----------------------------------------------------------
# This routine calculates the stellar density
# Based on eq. 30 from Winn., 2014
# Assuming the companion is too small
# Input:
# P -> Period
# a -> semi-major axis
# Output:
# rho -> stellar density
# -----------------------------------------------------------


def get_rhostar(P, a):
    P = P * 24. * 3600.  # s
    rho = 3. * np.pi * a**3 / (G_cgs * P * P)
    return rho


# -----------------------------------------------------------
# Return the equilibrium temeprature given the stellar temperature
# albedo, stellar radius and distance to the star
def get_teq(Tstar, albedo, rstar, a):
    Tp = Tstar*(1.0 - albedo)**(0.25)
    Tp = (rstar/2.0/a)**(0.5) * Tp
    return Tp

# -----------------------------------------------------------
# find_vals_perc -> find the median and the errors within
#  a 68% credible interval
# input:
#       x -> vector with a minimum size nconv
#   nconv -> the last nconv points to be taken account
#            in the gaussian fit
# output:
#        med -> median value
#	mine -> left error (50% - 16%)
#	maxe -> right error (84% - 50%)
# -----------------------------------------------------------


def find_vals_perc(x, sf=1.0, prob=68.3):
    # With a 68% confidence interval
    mnval = 50.0 - prob/2.0
    mxval = 50.0 + prob/2.0
    mine, med, maxe = np.percentile(x, [mnval, 50.0, mxval])
    maxe = (maxe - med) / sf
    mine = (med - mine) / sf

    return med, mine, maxe


# -----------------------------------------------------------
def best_value(vector, loglike, cual):
    if (cual == 'median'):
        result = np.median(vector)
    elif(cual == 'mode'):
        result = my_mode(vector)
    elif(cual == 'map'):
        maxindex = np.argmax(loglike)
        result = vector[maxindex]

    return result


# -----------------------------------------------------------
# This routine calculates the mode of a vector
# The vector in divided in bins and count the maximum value
def my_mode(vector, bins=50):
    dx = np.max(vector) - np.min(vector)
    dx = dx / bins
    b = np.sort(vector)
    i = 0
    j = 0
    o = 0
    limite = np.min(vector) + dx
    if (dx > 1e-10):  # the parameter is fixed
        while(o < len(b)):
            if (b[o] < limite):
                i = i + 1
                if (i > j):
                    j = i
                    maximo = limite - dx/2.0
                o = o + 1
            else:
                i = 0
                limite = limite + dx
    else:
        maximo = np.median(vector)

    return maximo

# -----------------------------------------------------------


def mode_and_99(vector):
    a = my_mode(vector)
    d, b, c = find_vals_perc(vector, sf=1.0, prob=98)
    b = d - b
    c = c + d

    return a, b, c

#Function that returns the bad chains
def categorise_chains(pos, nchains):
    # Let us find the good indixes for the cluster
    # We have n walkers

    nwalkers = max(nchains) + 1

    print('STARTING CHAIN CLUSTERING')
    print('Initial number of chains:', nwalkers)

    # This variable will have each walker information
    pos_walkers = [None]*nwalkers
    pos_mean = np.zeros(nwalkers)
    for i in range(nwalkers):
        ichain = nchains == i
        pos_walkers[i] = pos[ichain]

    # The mean of each walker
    for i in range(0, nwalkers):
        #pos_mean[i] = np.mean(pos_walkers[i])
        #Let's see if instaed of working with the mean of the chains, we check the value at the beggining
        #pos_mean[i] = pos_walkers[i][0]
        pos_mean[i] = np.mean(pos_walkers[i][:10])


    sorted_indices = np.argsort(pos_mean)

    good_chain = []
    bad_chain = []
    otros = np.ones(shape=(len(sorted_indices)))
    # Let us kill all the walkers 5 times the minimum
    for i,chain in enumerate(sorted_indices[:-1]):
        otros = np.logical_and(otros,(chain != sorted_indices))
        q1 = np.percentile(pos_mean[otros], 25)
        q3 = np.percentile(pos_mean[otros], 75)
        iqr = q3 - q1
        otros_walkers = np.median(pos_mean[otros]) - clustering_sigma*iqr
        if (pos_mean[chain] > otros_walkers):
            # We are saving the good chain labels
            good_chain.append(chain)
        else:
            print(f"removing chain {chain}, iteration {i}")
            # We are saving the good chain labels
            bad_chain.append(chain)
    #add last index
    good_chain.append(sorted_indices[-1])

    new_nwalkers = len(good_chain)

    print('Final number of chains:', new_nwalkers)

    return good_chain, bad_chain

#
#This function calcualtes the time of periastron given
#Time of minimum conjunction, eccentricit, angle of periastron, period.
def find_tp(T0,e,w,P):
  ereal = e + np.cos(np.pi/2. - w)
  eimag = np.sqrt(1. - e*e) * np.sin( np.pi / 2. - w )
  theta_p = np.arctan2(eimag,ereal)
  theta_p = theta_p - e * np.sin(theta_p)
  Tp = T0 - theta_p * P / 2. / np.pi
  return Tp

# -----------------------------------------------------------


def bin_data_old(tvector, fvector, rvector, tbin):
    tvector = np.asarray(tvector)
    fvector = np.asarray(fvector)
    rvector = np.asarray(rvector)

    # Define bin edges and centers
    leftt = min(tvector)
    right = max(tvector)

    # Generate bins
    bin_edges = np.arange(leftt, right + tbin, tbin)
    bin_centers = bin_edges[:-1] + tbin / 2.  # Compute bin centers

    # Initialize arrays for binned data
    fbined = []
    rbined = []

    # Loop through bins and calculate mean within each bin
    for i in range(len(bin_edges) - 1):
        # Use boolean indexing to find data within the current bin
        bin_mask = (tvector >= bin_edges[i]) & (tvector < bin_edges[i+1])

        # Safeguard for empty bins
        if np.any(bin_mask):
            fbined.append(np.mean(fvector[bin_mask]))
            rbined.append(np.mean(rvector[bin_mask]))
        else:
            fbined.append(np.nan)
            rbined.append(np.nan)

    # Convert to arrays
    fbined = np.asarray(fbined)
    rbined = np.asarray(rbined)
    bin_centers = np.asarray(bin_centers)

    # Remove NaN bins
    valid_mask = ~np.isnan(rbined)
    bin_centers = bin_centers[valid_mask]
    fbined = fbined[valid_mask]
    rbined = rbined[valid_mask]

    return bin_centers, fbined, rbined

def bin_data(tvector,fvector,rvector,tbin=0.1):
    leftt = min(tvector)
    right = leftt + tbin
    xbined = []
    fbined = []
    rbined = []
    while ( leftt < max(tvector) - tbin/2 ):
        fdummy = []
        rdummy = []
        for i in range(0,len(tvector)):
            if ( tvector[i] > leftt and tvector[i] < right ):
                fdummy.append(fvector[i])
                rdummy.append(rvector[i])
        fbined.append(np.mean(fdummy))
        rbined.append(np.mean(rdummy)/np.sqrt(len(rdummy)))
        #fbined.append(np.average(fdummy,weights=rdummy))
        xbined.append((leftt + tbin/2.))
        leftt = leftt + tbin
        right = right + tbin
    fbined = np.asarray(fbined)
    rbined = np.asarray(rbined)
    return xbined, fbined, rbined

def bin_data_tr(tvector,fvector,rvector,tbin=10./60/24.):
    leftt = min(tvector)
    right = leftt + tbin
    xbined = []
    fbined = []
    rbined = []
    while ( leftt < max(tvector) - tbin/2 ):
        fdummy = []
        rdummy = []
        for i in range(0,len(tvector)):
            if ( tvector[i] > leftt and tvector[i] < right ):
                fdummy.append(fvector[i])
                rdummy.append(rvector[i])
        fbined.append(np.mean(fdummy))
        rbined.append(np.mean(rdummy))
        #fbined.append(np.average(fdummy,weights=rdummy))
        xbined.append((leftt + tbin/2.)*24.)
        leftt = leftt + tbin
        right = right + tbin
    fbined = np.asarray(fbined)
    rbined = np.asarray(rbined)
    return xbined, fbined, rbined

# -----------------------------------------------------------
# Define a mapping from common units to LaTeX-friendly units
unit_to_latex = {
    'd': r'$\mathrm{days}$',
    'days': r'$\mathrm{days}$',
    'radians': r'$\mathrm{radians}$',
    'deg': r'$\mathrm{deg}$',
    'K': r'$\mathrm{K}$',
    'g/cm^3': r'$\mathrm{g\,cm^{-3}}$',
    'km/s': r'$\mathrm{km\,s^{-1}}$',
    'km/s/day': r'$\mathrm{km\,s^{-1}\,day^{-1}}$',
    'km/s/day^2': r'$\mathrm{km\,s^{-1}\,day^{-2}}$',
    'm/s': r'$\mathrm{m\,s^{-1}}$',
    'm/s/day': r'$\mathrm{m\,s^{-1}\,day^{-1}}$',
    'm/s/day^2': r'$\mathrm{m\,s^{-1}\,day^{-2}}$',
    'h': r'$\mathrm{hours}$',
    'cm/s^2': r'$\mathrm{cm\,s^{-2}}$',
    'AU': r'$\mathrm{AU}$',
    'R_solar': r'$\mathrm{R_{\odot}}$',
    'R_Earth': r'$\mathrm{R_{\oplus}}$',
    'R_Jupiter': r'$\mathrm{R_{J}}$',
    'M_solar': r'$\mathrm{M_{\odot}}$',
    'M_Earth': r'$\mathrm{M_{\oplus}}$',
    'M_Jupiter': r'$\mathrm{M_{J}}$',
    'F_Earth': r'$\mathrm{F_{\oplus}}$',
    '': ''  # For parameters that have no units
}

# Dictionary to convert numbers to their string equivalents
number_to_word = {
    '0': 'zero',
    '1': 'one',
    '2': 'two',
    '3': 'three',
    '4': 'four',
    '5': 'five',
    '6': 'six',
    '7': 'seven',
    '8': 'eight',
    '9': 'nine'
}

# Function to convert param names to LaTeX commands
def to_latex_command(param_name):
    # Remove underscores and replace numbers with their string equivalents
    param_name_command = param_name.replace('_', '').replace(' ', '')

    # Replace numbers with their word equivalents
    for num, word in number_to_word.items():
        param_name_command = param_name_command.replace(num, word)

    # Prepend the LaTeX command symbol
    return f'{param_name_command}'

def print_values(vector, var, unit):
    # vector is the posterior vectors
    # var is the label of the variable to save
    # unit is the label of the variable to save
    medv, minv, maxv = find_vals_perc(vector)
    nd = 1
    if (abs(minv) > 1e-20 and abs(maxv) > 1e-20):
        try:
            nd = int(np.log10(max(1./minv, 1./maxv))) + 2
        except:
            pass
    # Calculate dynamic widths based on the size of the variables
    opars.write('{:10s} = {:4.7f} - {:4.7f} + {:4.7f} {:8s} \n'.format(var, medv, minv, maxv, unit))
    nd=max(0,nd) #Avoid negative numbers
    rmedv = (float(round(medv,nd)))
    rminv = (float(round(minv,nd)))
    rmaxv = (float(round(maxv,nd)))
    vartex = to_latex_command(var)
    try:
        unittex = unit_to_latex[unit]
    except:
        unittex = ''
    if (rminv == rmaxv ):
        otex.write('\\newcommand{{\\{:s}}}[1][{:s}]{{ $ {:.{width}f} \pm {:.{width}f} $~#1 }} \n'.format(vartex,unittex,rmedv,rmaxv,width=int(nd)))
    else:
        otex.write('\\newcommand{{\\{:s}}}[1][{:s}]{{ $ {:.{width}f}_{{-{:.{width}f}}}^{{+{:.{width}f}}} $~#1 }} \n'.format(vartex,unittex,rmedv,rminv,rmaxv,width=int(nd)))
    if (is_print_mode):
        medv, minv, maxv = mode_and_99(vector)
        opars.write('%10s  %4.7f , %4.7f , %4.7f %8s (mode, 1 percent, 99 percent) \n' % (
            '', medv, minv, maxv, unit))


def line_prepender(filename, text):
            with open(filename, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                for o in range(len(text)):
                    f.write(text[o])
                f.write('\n' + content)






def create_fco_file(fname,outfname,chunk_size=0.5):
    """
    This function creates a RV file ready to be run using Float Chunk Offset method with pyaneti
    The input file fname has to contain four columns
    time RV RV_error instrument_label
    chunk_size delimits the lenght of the chunk where the code will search for at least two points
    """
    #Read the file
    print(f"Transforming {fname} file to be used in the FCO method")
    t, rv, e = np.loadtxt(fname,usecols=(0,1,2),unpack=True)
    tel = np.loadtxt(fname,usecols=(3),dtype=str)

    #Detect the telescope labels
    tl = []
    for l in tel:
        if l not in tl:
            tl.append(l)
    print("Detected {} telescope labels:".format(len(tl)))
    print(tl)

    good_index = []
    label = []
    n = 0
    t_f  = []
    rv_f = []
    e_f  = []
    c_f  = []
    l_f  = []

    #Check the chunks per each instrument
    for j in range(len(tl)):
        #print("Checking for {}".format(tl[j]))
        o = 0
        #Let us find all the times for a given instrument
        ins = tel == tl[j]
        tinst = t[ins]
        #print("{} observations for the current telescope".format(len(tinst)))
        lmin = min(tinst) - 1e-2
        lmax = min(tinst) + chunk_size
        i = np.arange(len(tinst))
        while o < len(tinst) :
            index =  (tinst > lmin) & (tinst < lmax)
            index = i[index]
            if len(index) > 1:
                t_f.append(tinst[index])
                rv_f.append(rv[ins][index])
                e_f.append(e[ins][index])
                c_f.append([str(n)]*len(index))
                l_f.append([tl[j]]*len(index))
                n += 1
                o += len(index)
            else:
                o += 1
            if o >= len(tinst):
                break
            lmin = tinst[o] - 1e-2
            lmax = tinst[o] + chunk_size


    #now we have the chunks for every instrument, let us put everytthing in a big array
    t_f = np.concatenate(t_f)
    rv_f = np.concatenate(rv_f)
    e_f = np.concatenate(e_f)
    c_f = np.concatenate(c_f)
    l_f = np.concatenate(l_f)

    #Write the output file
    with open(outfname,'w') as f:
        f.write('#time RV eRV chunk_label instrument_label \n')
        for i in range(len(t_f)):
                f.write('{} {} {} {} {} \n'.format(t_f[i],rv_f[i],e_f[i],c_f[i],l_f[i]))

#------------------------------------------------------------------------#
#            Automatic creation of input for TANGO                       #
#------------------------------------------------------------------------#


def tango_params(param, vec, parhs=True):
    vlen = len(vec)
    letra = param + ' = '
    if (parhs):
        letra = letra + ' [ '
    for o in range(0, vlen):
        letra = letra + str(np.median(vec[o]))
        if (o < vlen - 1):
            letra = letra + ','

    if (parhs):
        letra = letra + ' ]'
    letra = letra + '\n'

    return letra

# This routine create an input file to create animations using tango
# github.com/oscaribv/tango


def create_tango_input():
    tangof = outdir+'/'+star+'_tango_input.py'
    tf = open(tangof, 'w')

    tf.write('#Input file for tango\n')
    tf.write('#system:' + star+'\n')
    tf.write('#Created automatically with pyaneti\n')

    tf.write('\n')

    tf.write('#Data file with the flattened light curve\n')
    tf.write('lcname = \''+star+'_new_lc.dat\'\n')

    tf.write('#--------------------------------------------------------------------\n')
    tf.write('#                 Planet and orbit parameters\n')
    tf.write('# Each parameter is a list in which each element\n')
    if (nplanets == 1):
        tf.write('# correspond to a planet. For this case, there is ' +
                 str(nplanets)+' planet\n')
    else:
        tf.write('# correspond to a planet. For this case, there are ' +
                 str(nplanets)+' planets\n')
    tf.write('#--------------------------------------------------------------------\n')

    tf.write('\n')

    # Orbital period (days)
    tf.write(tango_params('P', P_vec))
    tf.write(tango_params('T0', T0_vec))
    tf.write(tango_params('e', e_vec))
    tf.write(tango_params('w', w_vec))
    tf.write(tango_params('a', ar_vec))
    tf.write(tango_params('inclination', i_vec))
    tf.write(tango_params('rp', rr_vec))
    tf.write(tango_params('u1', [u1_vec], False))
    tf.write(tango_params('u2', [u2_vec], False))

# Integration time of the data
    tf.write('t_cad = ' + str(t_cad) + ' \n')
    tf.write('n_cad = ' + str(n_cad) + ' \n')

    tf.write('\n')

    # Calculate the transit duration for planet b to create the time ranges
    tfull = np.median(trt_vec[0])
    tfull = tfull/24.
    mit0 = np.median(T0_vec[0])

    tf.write('#--------------------------------------------------------------------\n')
    tf.write('#              Animation controls \n')
    tf.write('#--------------------------------------------------------------------\n')
    tf.write('#Window size to show the data (days)\n')
    tf.write('size_time = 0.5\n')
    tf.write(
        '#1./(photograms per day) in this case the code will create a photogram each 7.2 min\n')
    tf.write('vel_time  = 1./200.\n')
    tf.write('#Animation minimum time (Be sure that you are using the same units as in your data file)\n')
    tf.write('tmin = '+str(mit0 - 2*tfull)+'\n')
    tf.write('#Animation maximum time (Be sure that you are using the same units as in your data file)\n')
    tf.write('tmax = '+str(mit0 + 2*tfull)+'\n')

    tf.write('\n')

    tf.write('#--------------------------------------------------------------------\n')
    tf.write('#                     Plot controls\n')
    tf.write('#--------------------------------------------------------------------\n')

    tf.write('\n')

    tf.write('#Control if we overplot the light curve model\n')
    tf.write('#You need to have installed pyaneti in your computer to use it\n')
    tf.write('is_plot_model = False\n')

    tf.write('\n')

    tf.write('#-----------------------------------------------------------------\n')
    tf.write('#                         END\n')
    tf.write('#-----------------------------------------------------------------\n')
