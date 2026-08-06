#Read pyaneti output file with all posteriors
newfile = f'{outdir}/{star}_all_data.dat'
df = pd.read_csv(newfile)
df = df.iloc[:, :-1] #Remove the last column that is empty

if is_clustering:
    #We have the posteriors, let us perform clustering to remove the bad chains
    good_chains, bad_chains = categorise_chains(df.log_likelihood,df.chain_number)
    #Let us remove the bad chains from the data frame
    df = df[~df.chain_number.isin(bad_chains)]
    df = df.copy()
    #let us save the file
    newfile = f'{outdir}/{star}_all_data_clustered.dat'
    df.to_csv(newfile)

#Extract the names of the sampled parameters
sampled_parameters = list(df.iloc[:,3:])


#Create list that will have the units of each column in the data frame
units = []
#Let us create columns with the stellar parameters to the computations
mstar = np.random.normal(loc=mstar_mean, scale=mstar_sigma,size=len(df.log_likelihood))
tstar = np.random.normal(loc=tstar_mean, scale=tstar_sigma,size=len(df.log_likelihood))
rstar = np.random.normal(loc=rstar_mean, scale=rstar_sigma,size=len(df.log_likelihood))
vsini = np.random.normal(loc=vsini_mean, scale=vsini_sigma,size=len(df.log_likelihood))
# Density from the input stellar parameters
irho = mstar / rstar**3 * 1.411  # Stellar density in cgs units


#Let us compute all the derived parameters
for i,planet in enumerate(plabels[:nplanets]):
# STARTING CALCULATIONS
    # Logarithmic transformation of the period
    if is_log_P:
        df[f'p{planet}'] = 10.0 ** df[f'p{planet}']
        units.append(['days'])

    # Stellar density and semi-major axis calculation
    if sample_stellar_density:
        # Using the stellar density of the first planet as a reference
        if planet == plabels[0]:
            miden = df[f'rhostar{planet}'].copy()  # Copy the stellar density for safety
        # Calculate the semi-major axis and store it in the corresponding column
        for index, row in df.iterrows():
            df.at[index, f'arstar{planet}'] = pti.rhotoa(miden[index], row[f'p{planet}'])[0]
        units.append([''])

    # Logarithmic transformation of the semi-amplitude K
    if is_log_k:
        df[f'k{planet}'] = 10.0 ** df[f'k{planet}']
        units.append(['km/s'])

    #Let us transform the doppler semi-amplitude to meters per second
    df[f'k{planet}_ms'] = df[f'k{planet}']*1e3
    units.append(['m/s'])


    # Conversion of eccentricity and argument of periastron if combined
    if is_ew:
        e_dum = df[f'esinomega{planet}'].copy()
        omega_dum = df[f'ecosomega{planet}'].copy()
        df[f'e{planet}'], df[f'omega{planet}'] = pti.ewto(e_dum, omega_dum)
        units.append(['','radians'])

    # Impact parameter and inclination calculation
    if is_b_factor:
        df[f'i{planet}'] = pti.btoi(df[f'b{planet}'], df[f'arstar{planet}'], df[f'e{planet}'], df[f'omega{planet}'])
        units.append(['radians'])
    else:
        # Calculate the impact parameter according to Winn 2014 equation (7)
        df[f'i{planet}'] = np.arccos(df[f'b{planet}'])  # Initial guess for inclination
        df[f'b{planet}'] = df[f'arstar{planet}'] * np.cos(df[f'i{planet}']) * \
                          ((1. - df[f'e{planet}']**2) / (1.0 + df[f'e{planet}'] * np.sin(df[f'omega{planet}'])))
        df[f'i{planet}'] = np.array(df[f'i{planet}'])  # Ensuring that the inclination is an array
        units.append([''])

    # Convert inclination from radians to degrees and store in i_deg column
    df[f'i_deg{planet}'] = np.degrees(df[f'i{planet}'])  # Converts radians to degrees
    units.append(['deg'])

    # Time of periastron passage
    df[f'Tpe{planet}'] = np.array(df[f't0{planet}'])
    for index, row in df.iterrows():
        df.at[index, f'Tpe{planet}'] = pti.find_tp(row[f't0{planet}'], row[f'e{planet}'], row[f'omega{planet}'], row[f'p{planet}'])
    units.append(['d'])

    #Let us just compute the values for the planet if we sampled transit parameters for this planet
    if fit_tr[i]:

        # Calculate equilibrium temperature assuming albedo = 0
        df[f'Teq{planet}'] = get_teq(tstar, 0.0, 1.0, df[f'arstar{planet}'])
        units.append(['K'])

        # Calculate the stellar density from transit data (Eq. 30, Winn 2014)
        df[f'ds{planet}'] = get_rhostar(df[f'p{planet}'], df[f'arstar{planet}'])
        units.append(['g/cm^3'])

        # Compute the semi-major axis in real units (SI units, in AU)
        df[f'a{planet}'] = df[f'arstar{planet}'] * rstar * S_radius_SI / AU_SI  # Result in AU
        units.append(['AU'])

        if nradius == 1:
            df[f'rp{planet}'] = df[f'rprstar{planet}'] * rstar
            if planetary_units == 'solar':
                units.append(['R_solar'])
            elif planetary_units == 'earth':
                df[f'rp{planet}'] *=  S_radius_SI / E_radius_e_SI  # Convert radius to Earth units
                units.append(['R_Earth'])
            elif planetary_units == 'jupiter':
                df[f'rp{planet}'] *=  S_radius_SI / J_radius_e_SI   # Convert radius to Jupiter units
                units.append(['R_Jupiter'])
        else:
            for band in bands:
                df[f'rp{band}{planet}'] =  df[f'rprstar{band}{planet}'] * rstar
                if planetary_units == 'solar':
                    units.append(['R_solar'])
                elif planetary_units == 'earth':
                    df[f'rp{band}{planet}'] *=  S_radius_SI / E_radius_e_SI  # Convert radius to Earth units
                    units.append(['R_Earth'])
                elif planetary_units == 'jupiter':
                    df[f'rp{band}{planet}'] *=  S_radius_SI / J_radius_e_SI   # Convert radius to Jupiter units
                    units.append(['R_Jupiter'])

        #Let's find the planet radius that we will use for the next calculations
        #We assume that small differences in radii for different bands are not important and we use the value of the first band
        radius_pattern = f"rp{planet}"
        # Extract the first element that contains the string pattern
        radius_name = next((x for x in list(df) if radius_pattern in x), None)
        radius_pattern = f"rprstar{planet}"
        # Extract the first element that contains the string pattern
        radius_rstar_name = next((x for x in list(df) if radius_pattern in x), None)
        #radius_name contains the df name of the first element of radius rp for {planet}

        # Total transit duration (from contact 1 to 4), equation from Winn 2014
        # Calculate the eccentricity factor for transit durations
        ec_factor = np.sqrt(1.0 - df[f'e{planet}']**2) / (1.0 + df[f'e{planet}'] * np.sin(df[f'omega{planet}']))
        df[f'trt{planet}'] = np.sqrt((1.0 + df[radius_rstar_name])**2 - df[f'b{planet}']**2) / (df[f'arstar{planet}'] * np.sin(df[f'i{planet}']))
        df[f'trt{planet}'] = df[f'p{planet}'] / np.pi * np.arcsin(np.clip(df[f'trt{planet}'], -1, 1)) * ec_factor * 24.0  # Convert to hours
        units.append(['h'])

        # Full transit duration (from contact 2 to 3)
        df[f'tri{planet}'] = np.sqrt((1.0 - df[radius_rstar_name])**2 - df[f'b{planet}']**2) / (df[f'arstar{planet}'] * np.sin(df[f'i{planet}']))
        df[f'tri{planet}'] = df[f'p{planet}'] / np.pi * np.arcsin(np.clip(df[f'tri{planet}'], -1, 1)) * ec_factor * 24.0  # Convert to hours
        units.append(['h'])

        #Estimate planet insolation
        # Stellar luminosity in solar units
        Ls = (rstar) ** 2 * (tstar / S_Teff) ** 4
        # Planet insolation flux received at Earth
        df[f'Fp{planet}'] = Ls / df[f'a{planet}'] ** 2
        units.append(['F_Earth'])

        # Estimate the maximum magnitude difference of a star that can generate a transit
        # Based on eq. 4 of Vanderburg et al., 2019, ApJL, 881, L19
        t12 = (df[f'trt{planet}'] - df[f'tri{planet}']) / 2.  # assuming ingress and egress are the same
        t13 = df[f'tri{planet}'] + t12
        # Calculate the maximum magnitude difference
        df[f'delta_mag{planet}'] = t12 ** 2 / (t13 ** 2 * df[radius_rstar_name] ** 2)
        df[f'delta_mag{planet}'] = 2.5 * np.log10(df[f'delta_mag{planet}'])
        units.append([''])


    if fit_rv[i]:

        # Compute the planet mass in real units (solar masses)
        df[f'mp{planet}'] = planet_mass(
            mstar,                          # Stellar mass in solar masses
            df[f'k{planet}']*1e3,              # Semi-amplitude K converted to m/s
            df[f'p{planet}'],                     # Orbital period in days
            df[f'e{planet}'],                     # Eccentricity
            df[f'i{planet}']                      # Inclination in radians
        )
        if planetary_units == 'solar':
            units.append(['M_solar'])
        elif planetary_units == 'earth':
            df[f'mp{planet}'] *=  S_GM_SI / E_GM_SI  # Convert mass to Earth units
            units.append(['M_Earth'])
        elif planetary_units == 'jupiter':
            df[f'mp{planet}'] *=  S_GM_SI / J_GM_SI   # Convert mass to Jupiter units
            units.append(['M_Jupiter'])

    #Values that can only be computed for transiting and RV planets
    if fit_tr[i] and fit_rv[i]:

        # Estimate planet density
        df[f'pden{planet}'] = df[f'mp{planet}'] / df[radius_name] ** 3  # in solar units
        if planetary_units == 'solar': df[f'pden{planet}'] = df[f'pden{planet}'] * S_den_cgs  # Convert to g/cm^3
        if planetary_units == 'earth': df[f'pden{planet}'] = df[f'pden{planet}'] * E_den_cgs  # Convert to g/cm^3
        if planetary_units == 'jupiter': df[f'pden{planet}'] = df[f'pden{planet}'] * J_den_cgs  # Convert to g/cm^3
        units.append(['g/cm^3'])

        # Estimate surface gravity from period and semi-major axis (eq. (31) Winn)
        df[f'pgra{planet}'] = (
            (df[f'p{planet}'] * 24. * 3600.) * (df[radius_rstar_name] / df[f'arstar{planet}']) ** 2 * np.sin(df[f'i{planet}'])
        )
        df[f'pgra{planet}'] = (
            2. * np.pi * np.sqrt(1. - df[f'e{planet}'] ** 2) * (df[f'k{planet}'] * 1.e5) / df[f'pgra{planet}']  # in cm/s^2
        )
        units.append(['cm/s^2'])

        # Estimate surface gravity from derived parameters (in solar units)
        df[f'pgra2{planet}'] = df[f'mp{planet}'] / df[radius_name] ** 2
        if planetary_units == 'solar': df[f'pgra2{planet}'] *=   28.02 * 981.  # in cm/s^2
        if planetary_units == 'earth': df[f'pgra2{planet}'] *=  1 * 981.  # in cm/s^2
        if planetary_units == 'jupiter': df[f'pgra2{planet}'] *= 2.528 * 981.  # in cm/s^2
        units.append(['cm/s^2'])

        # Assign TSM scaling factor based on the mean radius
        mean_radius =  df[radius_name].mean()
        conditions = [
        mean_radius < 1.5,
        (mean_radius > 1.5) & (mean_radius < 2.75),
        (mean_radius > 2.75) & (mean_radius < 4.0),
        mean_radius >= 4.0
        ]
        scaling_factors = [0.190, 1.26, 1.28, 1.15]
        TSM_factor = np.select(conditions, scaling_factors)

        # Perform the final TSM calculation
        df[f'TSM{planet}'] = (
        TSM_factor * df[radius_name]**3 * df[f'Teq{planet}'] /
        df[f'mp{planet}'] / rstar**2 * 10.**(-mag_j/5.)
        )
        units.append([''])

    if is_single_transit:
        # Initialize a column in the DataFrame for transit durations
        df[f'tr_dur{planet}'] = np.nan  # Set default values to NaN

        # Initialize arrays to store z_vec
        z_vec = [None] * len(df)  # Initialize a list for z vectors

        for idx, row in df.iterrows():
            # Create parameter list for find_z function
            pars = [
                row[f'Tpe{planet}'],  # Time of periastron passage
                row[f'p{planet}'],    # Orbital period
                row[f'e{planet}'],    # Eccentricity
                row[f'omega{planet}'], # Argument of periastron
                row[f'i{planet}'],    # Inclination
                row[f'arstar{planet}']  # Semi-major axis
            ]

            # Use mean scaled radius for this planet
            rplanet = np.mean(row[f'rprstar{planet}'])

            # Define the time array with minute precision
            tiempo = np.arange(min(lc_time), max(lc_time), 1./60./24.)  # Minute precision

            # Find the z vector (distance between planet and star center in stellar radius units)
            z_vec[idx] = pti.find_z(tiempo, pars)

            # Identify the points where the planet is transiting (z < 1 + rplanet)
            tr_index = z_vec[idx] < (1 + rplanet)
            tr_time = tiempo[tr_index]

            # Calculate the transit duration (time between the first and last contact)
            if len(tr_time) > 0:  # Ensure that we have transit times
                df.at[idx, f'tr_dur{planet}'] = max(tr_time) - min(tr_time)
            else:
                df.at[idx, f'tr_dur{planet}'] = np.nan  # Handle cases where there is no valid transit


        # Estimate planet velocity at transit time using Eq. 1 from Osborn et al., 2016
        vpl_vec = np.array(2. * np.sqrt((1. + df[f'rprstar{planet}'])**2 - df[f'b{planet}']**2) / (df[f'tr_dur{planet}'] * 24. * 60. * 60.))

        # Store the velocity in the DataFrame
        df[f'vpl{planet}'] = vpl_vec

        # Now estimate planetary orbital period assuming the orbit is circular using Eq. 2 from Osborn et al., 2016
        P_circ_vec = 8. * np.pi**2 * G_cgs * irho / 3. / vpl_vec**3

        # Convert P_circ_vec to days and store in DataFrame
        df[f'P_circ{planet}'] = P_circ_vec / 24. / 3600.

if any(fit_tr) == True:
    #Parameters independent of the number of planets
    if nbands == 1:
        df[f'u1'], df[f'u2'] = pti.get_us(df[f'q1'],df[f'q2'])
        units.append(['']*2)
    else:
        for band in bands:
            df[f'u1{band}'], df[f'u2{band}'] = pti.get_us(df[f'q1{band}'],df[f'q2{band}'])
            units.append(['']*2)

if is_linear_trend == 'u':
    df['linear_trend_ms'] = df['linear_trend'] * 1e3
    units.append(['m/s/day'])

if is_quadratic_trend == 'u':
    df['quadratic_trend_ms'] = df['quadratic_trend'] * 1e3
    units.append(['m/s/day^2'])

derived_param_units = np.concatenate(units)
derived_parameters = [col for col in df.columns[3:] if col not in sampled_parameters]

#Create the output files
out_params_file = outdir+'/'+star+'_params.dat'
out_tex_file = outdir+'/'+star+'_params.tex'
opars = open(out_params_file, 'w')
otex = open(out_tex_file, 'w')

#Let us extract the indices of the parameters that we are sampl
sampled_indices = prior_flags != 'f'
#Dimension of the sampled space
npars = len(prior_flags[sampled_indices])

ndata = len(rv_time) + len(lc_time)

#Compute the median of all the parameters to be used to compute the plots
df_medians = df[sampled_parameters].median()

log_like_total = pyaneti_likelihood_parallel(df_medians[sampled_indices])

#log_like_total = max(df.log_likelihood)
bic_from_loglikelihood = np.log(ndata)*npars - 2.0*log_like_total
aic_from_loglikelihood = 2.0*npars - 2.0*log_like_total

separator = '-' * 62

# Start summary
opars.write('\n')
opars.write(f'{separator}\n')
opars.write('Summary:\n')
opars.write(f'N_chains         = {nwalkers:>8d} \n')
opars.write(f'N_iter           = {nconv:>8d} \n')
opars.write(f'thin_factor      = {thin_factor:>8d} \n')
opars.write(f'N_rv_data        = {len(rv_time):>8d} \n')
opars.write(f'N_tr_data        = {len(lc_time):>8d} \n')
opars.write(f'N_data           = {ndata:>8d} \n')
opars.write(f'N_pars           = {npars:>8d} \n')
opars.write(f'dof              = {ndata - npars:>8d} \n')
opars.write(f'ln likelihood    = {log_like_total:>8.4f}\n')
opars.write(f'BIC              = {bic_from_loglikelihood:>8.4f}\n')
opars.write(f'AIC              = {aic_from_loglikelihood:>8.4f}\n')
opars.write(f'{separator}\n')

# Stellar parameters
opars.write('             INPUT STELLAR PARAMETERS\n')
opars.write(f'{separator}\n')
opars.write(f'M_*     = {mstar_mean:>8.4f} - {mstar_sigma:<8.4f} + {mstar_sigma:<8.4f} solar masses\n')
opars.write(f'R_*     = {rstar_mean:>8.4f} - {rstar_sigma:<8.4f} + {rstar_sigma:<8.4f} solar radii\n')
opars.write(f'T_*     = {tstar_mean:>8.4f} - {tstar_sigma:<8.4f} + {tstar_sigma:<8.4f} K\n')
opars.write(f'vsini   = {vsini_mean:>8.4f} - {vsini_sigma:<8.4f} + {vsini_sigma:<8.4f} km/s (for RM effect estimation)\n')
opars.write(f'J mag   = {mag_j:>8.4f}  (for TSM estimation)\n')
opars.write(f'{separator}\n')

# Sampled parameters
opars.write('-------------       Sampled parameters   ---------------------\n')
opars.write(f'{separator}\n')


# Print the formatted rows
for i, p in enumerate(sampled_parameters):
    print_values(df[p],p,sampled_param_units[i])

opars.write(f'{separator}\n')
opars.write('---------------     Derived parameters   ---------------------\n')
opars.write(f'{separator}\n')

# Print the formatted rows
for i, p in enumerate(derived_parameters):
    print_values(df[p],p,derived_param_units[i])

opars.write(f'{separator}\n')
opars.write(f'{separator}\n')

opars.close()
otex.close()

# Print the output in the screen
with open(out_params_file) as dummy_file:
    for l in dummy_file:
        print(l,end='')


