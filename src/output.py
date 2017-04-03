
#Print the values
execfile('src/print_values.py')

#Create plots
execfile('src/plot_data.py')

plot_chains()

#plot_likelihood()

if ( is_plot_histogram ):
  plot_histogram()

if ( is_plot_correlations ):
  plot_correlations()

 #PLOT TRANSIT
if ( total_tr_fit ):
  plot_transit_nice()
  plot_all_transits()

#PLOT RV CURVE
if ( total_rv_fit ):
  plot_rv_all_data()
  plot_rv_mp()
