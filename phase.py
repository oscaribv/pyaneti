import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import jdcal as jd
import sys

#
def find_night(days,fase,sn,en):
	night=[]
	fnight=[]
	#This assumes that the night ends the next day that starts
	sndf = 24 % sn #star night day fraction
	n1 = sndf / 24.0
	n2 = en / 24.0
	#In JD format midnight is at 0.5
	#Then the night will start at 0.5-n1 and end at 0.5+n2
	start = 0.5 - n1
	end   = 0.5 + n2
	#Then all the rages between these values will be night
	for i in range(0,len(days)-1):
		if ( days[i] > int(days[i])+start and days[i] < int(days[i])+end ):
			#If our range falls in the night, let us store it in the
			#night vector (we save date and phase of that date)
			night.append(days[i])
			fnight.append(fase[i])
	return night, fnight

def find_eclipse(T0,P,D,jds,jdf):
	eclipses = []
	#How many days since T0?
	hmd = jds - T0
	if ( hmd < 0 ):
		print "You have not observed this planet yet!"
		sys.exit("Error message")
	#How much time do we need for the next eclipse?
	#The residual will give us the elapsed period
	elapsed = np.float64(hmd) % np.float64(P)
	remaining =  P - elapsed
	#The following eclipse will be on
	nextec = jds + remaining
	#Now we have the midpoint of the next eclise
	#Let us find the starting and ending points (by using D)
	while ( nextec < jdf ):
		eclipses.append(nextec - 0.5*D)
		eclipses.append(nextec + 0.5*D)
		nextec = nextec + P
	return eclipses

#Here the main program starts

#Read the input parameters from input_phase.txt
idata = np.loadtxt('input_phase.txt',comments='#',unpack=True,dtype='S')

#Put the data in the correct variables
T0 = np.float(idata[0])
P  = np.float(idata[1])
D  = np.float(idata[2])
objeto = idata[3]
sd = [np.int(idata[4][0:4]),np.int(idata[4][5:7]),np.int(idata[4][8:10])]
fd = [np.int(idata[5][0:4]),np.int(idata[5][5:7]),np.int(idata[5][8:10])]
sn = np.float(idata[6][0:2])
en = np.float(idata[6][3:5])


#THE MAGIC BEGINS HERE

#Let us obtain starting date (sd) and final date (fd) in JD
jdsd = jd.gcal2jd(sd[0],sd[1],sd[2])
jdfd = jd.gcal2jd(fd[0],fd[1],fd[2])
sd_float = jdsd[0] + jdsd[1] + 0.5 
fd_float = jdfd[0] + jdfd[1] + 1.5 
#Number of points to plot (48 per day)
m = 48 * int(fd_float - sd_float)
#Let us create a vector with the days and fase
days = [None]*m
fase = [None]*m
rv   = [None]*m
#jump size
dm = (fd_float - sd_float) / m
#The calculation start at the midnight of the first night
days[0] = sd_float
rv[0]   = - np.sin(2.0*np.pi *(days[0] - T0) / P )
fase[0] = (days[0] - T0) / P 
fase[0] = (days[0] - T0) / P - int(fase[0])
#Now let us do it for all the days (by doing dm jumps)
for i in range(1,m):
	days[i] = days[i-1] + dm 
	rv[i]   = - np.sin(2.0*np.pi *(days[i] - T0) / P )
	fase[i] = (days[i] - T0) / P
	fase[i] = (days[i] - T0) / P - int(fase[i])

#When will it be night?
nightrv,fnightrv = find_night(days,rv,sn,en)
nightfa,fnightfa = find_night(days,fase,sn,en)

#Let us find the eclipses
eclipses = find_eclipse(T0,P,D,sd_float,fd_float)
#A dummy vector to have a y value for the plot
dummy = [0]*len(eclipses)

#Let us do a nice plot 

#Let us change the zero to the starting date
for i in range(0,m):
	days[i] = days[i] - sd_float

for i in range(0,len(nightrv)):
	nightrv[i] = nightrv[i] - sd_float
	nightfa[i] = nightfa[i] - sd_float

for i in range(0,len(eclipses)):
	eclipses[i] = eclipses[i] - sd_float

fname = objeto + '_' + idata[4] + '_to_' + idata[5] 

#Let us create some labels
plt.figure(1,figsize=(8,6))
plt.subplot(211)
plt.xlabel( ' T (JD - %7f)'%sd_float )
plt.ylabel( 'RV' )
plt.ylim(-1.2,1.2)
plt.plot(days,rv,'r-', label="RV/Phase")
plt.plot(nightrv,fnightrv,'bo', label="night")
plt.plot(eclipses,dummy,'gs', label="transit")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)
plt.subplot(212)
plt.xlabel( ' T (JD - %7f)'%sd_float )
plt.ylabel( 'phase' )
plt.ylim(0,1.1)
plt.plot(days,fase,'r-', label="Theoretical RV")
plt.plot(nightfa,fnightfa,'bo', label="night")
plt.savefig(fname + '.pdf')
plt.show()

#Let us make a nice print in the screen with hours and date

hdate = [None]*len(nightfa)

for i in range(0,len(hdate)):
	hdate[i] = jd.jd2gcal(sd_float,np.array(nightfa[i]))

dataf = open(fname+'.dat','w')	

dataf.write('#Object: %s\n'%objeto)
dataf.write('#Phases between %s and %s\n'%(idata[4],idata[5]))
dataf.write('# hour	phase\n')
dataf.write( '#------------------------#\n')
dataf.write('# %4i %2i %2i\n'%(hdate[0][0],hdate[0][1],hdate[0][2]))
dataf.write('#------------------------#\n')
for i in range(0,len(nightfa)):
	dataf.write('%2.2f	%1.4f\n'%(hdate[i][3]*24.0,fnightfa[i]))
	if (hdate[i][2] != hdate[i-1][2] and i > 1):
		dataf.write('#------------------------#\n')
		dataf.write(('# %4i %2i %2i\n'%(hdate[i][0],hdate[i][1],hdate[i][2])))
		dataf.write('#------------------------#\n')

dataf.close()

