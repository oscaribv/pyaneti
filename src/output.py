
#Print the values
execfile('src/print_values.py')

#Create plots
execfile('src/plot_data.py')

plot_chains()

if ( is_plot_histogram ):
  #plot_histogram()
  plot_histogram_2()

if ( is_plot_correlations ):
  #plot_correlations()
  plot_correlations_2()

 #PLOT TRANSIT
if ( fit_tr ):
  plot_transit_nice()
  if ( plot_all_tr ):
    plot_all_transits()

#PLOT RV CURVE
if ( fit_rv ):
  if ( nplanets == 1 ):
    plot_rv_all_data()
    #plot_rv_one()
    plot_rv_mp()
  else:
  #PLOT THE RV curves
    plot_rv_mp()
