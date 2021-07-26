# Print the values
exec(open('src/print_values.py').read())

# Create plots
exec(open('src/plot_data.py').read())
exec(open('src/plot_tr.py').read())
exec(open('src/plot_rv.py').read())

if (is_plot_chains):
    plot_chains()

if (is_corner_plot):
    create_corner_plot()
else:
    if (is_plot_posterior):
        plot_posterior()

if (is_plot_correlations):
    plot_correlations()

 # PLOT TRANSIT
if (total_tr_fit):
    plot_lightcurve_timeseries()
    create_folded_tr_plots()
#  if (nbands == 1):
#    plot_all_transits()
    # clean_transits(sigma_clean)
    # create_tango_input()

# PLOT RV CURVE
if (total_rv_fit):
    plot_rv_timeseries()
    plot_rv_phasefolded()
