
#Print the values
execfile('src/print_values.py')

#Create plots
execfile('src/plot_data.py')

if ( is_plot_chains ):
  plot_chains()
#  plot_postiter()

if ( is_corner_plot ):
  create_corner_plot()
else:
  if ( is_plot_posterior ):
    plot_posterior()

  if ( is_plot_correlations ):
    plot_correlations()

 #PLOT TRANSIT
if ( total_tr_fit ):
  plot_transit_nice()
  plot_all_transits()
  clean_transits(sigma_clean)
  create_tango_input()

#PLOT RV CURVE
if ( total_rv_fit ):
  plot_rv_all_data()
  plot_rv_mp()
