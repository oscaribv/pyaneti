import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import jdcal as jd
import sys

#
def find_night(days,fase):
	night=[]
	fnight=[]
	sn = 20 #night start in hours
	en = 7 	#night end in hours
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

#Here starts the main program

#USER INPUT VALUES 
#T0 of the eclise (JD)
T0 = 2307.72237 + 2454833.0
#Period (days)
P  = 1.673774
#Duration (hours)
D = 2
#Starting date (year,month,day)
sd = [2016,1,29]
#Ending date (year,month,day)
fd = [2016,2,3]
#Number of points to plot
m = 1000

#THE MAGIC BEGINS HERE

#durations in day fraction
D = D / 24.0

#Let us obtain starting date (sd) and final date (fd) in JD
jdsd = jd.gcal2jd(sd[0],sd[1],sd[2])
jdfd = jd.gcal2jd(fd[0],fd[1],fd[2])
sd_float = jdsd[0] + jdsd[1]
fd_float = jdfd[0] + jdfd[1]
#Let us create a vector with the days and fase
days = [None]*m
fase = [None]*m
#jump size
dm = (jdfd[1] - jdsd[1]) / m
#The calculation start at the midnight of the first night
days[0] = sd_float
fase[0] = - np.sin(2.0*np.pi *(days[0] - T0) / P )
#Now let us do it for all the days (by doing dm jumps)
for i in range(1,m):
	days[i] = days[i-1] + dm 
	fase[i] = - np.sin(2.0*np.pi *(days[i] - T0) / P )

#When will it be night?
night,fnight = find_night(days,fase)

#Let us find the eclipses
eclipses = find_eclipse(T0,P,D,sd_float,fd_float)
#A dummy vector to have a y value for the plot
dummy = [0]*len(eclipses)

#Let us do a nice plot 

#Let us change the zero to the starting date
for i in range(0,m):
	days[i] = days[i] - sd_float

for i in range(0,len(night)):
	night[i] = night[i] - sd_float

for i in range(0,len(eclipses)):
	eclipses[i] = eclipses[i] - sd_float

#Let us create some labels
plt.xlabel( ' T (JD - %7f)'%sd_float )
plt.ylabel( 'Phase' )
plt.ylim(-1.2,2.2)
plt.plot(days,fase,'r--', label="Theoretical RV")
plt.plot(night,fnight,'bo', label="night")
plt.plot(eclipses,dummy,'gs', label="eclipse")
plt.legend()
plt.savefig('test.png')
plt.show()
