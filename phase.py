import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import jdcal as jd
import sys

#Ask for the values
#print "Please enter the T0 and Period of the planet"
#T0 = input("T0 (JD):")
#P  = input("P (days):")

def find_night(days,fase):
	#Define your night (start and end)
	night=[]
	fnight=[]
	sn = 20 #hours
	en = 6 #hours
	#In JD the midnight is at 0.5
	#This assumes that the night ends the next day that starts
	sndf = 24 % sn #star night day fraction
	endf = 24 % ( 24 - en ) #end night day fraction
	n1 = sndf / 24.0
	n2 = endf / 24.0
	#Then the night will be at
	start = 0.5 - n1
	end   = 0.5 + n2
	#Then all the day rages between these values will be night
	for i in range(0,len(days)-1):
		if ( days[i] > int(days[i])+start and days[i] < int(days[i])+end ):
			night.append(days[i])
			fnight.append(fase[i])
	return night, fnight

#T0, period, duration, startday, endday
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
	#eclipses.append(nextec - 0.5*D)
	#eclipses.append(nextec + 0.5*D)
	while ( nextec < jdf ):
		eclipses.append(nextec - 0.5*D)
		eclipses.append(nextec + 0.5*D)
		nextec = nextec + P
	return eclipses

T0 = 2307.72237 + 2454833.0
P  = 1.673774
#duration
D = 2. / 24.

sd = [2016,1,29]
fd = [2016,2,2]

jdsd = jd.gcal2jd(sd[0],sd[1],sd[2])
jdfd = jd.gcal2jd(fd[0],fd[1],fd[2])

sd_float = jdsd[0] + jdsd[1]
fd_float = jdfd[0] + jdfd[1]

#
#Let us add the timezone
tz = 0
sd_float = sd_float + tz / 24.0
fd_float = fd_float + tz / 24.0
#
m = 1000
#Let us create a vector with the days and fase
days = [None]*m
fase = [None]*m
dm = (jdfd[1] - jdsd[1]) / m
days[0] = sd_float
fase[0] = - np.sin(2.0*np.pi *(days[0] - T0) / P )
for i in range(1,m):
	days[i] = days[i-1] + dm 
	fase[i] = - np.sin(2.0*np.pi *(days[i] - T0) / P )

#When will it be night?
night,fnight = find_night(days,fase)

#Let us find the eclipses
eclipses = find_eclipse(T0,P,D,sd_float,fd_float)
dummy = [0]*len(eclipses)

#Let us do a nice plot 
#for i in range(0,m):
#	days[i] = days[i] - (jdsd[1] + jdsd[0])

#for i in range(0,len(night)):
#	night[i] = night[i] - (jdsd[1] + jdsd[0])

plt.plot(days,fase,'r--',night,fnight,'bo',eclipses,dummy,'gs')
