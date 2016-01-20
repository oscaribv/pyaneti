import numpy as np
import datetime as dt
import jdcal as jd
import sys

#Ask for the values
print "Please enter the T0 and Period of the planet"
T0 = input("T0 (JD):")
P  = input("P (days):")

#What day is today?
now = dt.datetime.now()
year=now.year
month=now.month
day=now.day
print "Today is (y, m, d)", year, month, day

#Let us transform today to Julian Date
nowjd = jd.gcal2jd(year,month,day)
todayjd= np.float64(nowjd[0]) + np.float64(nowjd[1])
print "today JD is (days)", todayjd

#How many days since T0?
hmd = np.float64(todayjd) - np.float64(T0)

#If dum is negative, then you have not oberved the planet!
if ( hmd < 0 ): 
	print "You have not observed this planet yet!"
	sys.exit("Error message")

#How much time do we need for the next eclipse?
#The residual will give us the elapsed period
nextec = np.float64(hmd) % np.float64(P)


print "Today phase is: ", nextec * 360 / P


#The next eclipse will be in P - nextec days
#Let us add this to todayjd
eclipsejd = np.float64(todayjd) + np.float64(P) - np.float64(nextec)
#Now we know when the next eclipse will happen (in JD)

#We want to know when the next n eclispes will happen
n=input("Number of following eclipses you want to know?") 

for i in range(0,n):
	#The first next eclipse will be on eclipsejd
	#Then we just need to add the period i times
	nexteclipse = np.float64(eclipsejd) + i*np.float64(P)
	#Put the eclipse date in a human friendly way
	eclipse=jd.jd2gcal(2400000.5,nexteclipse-2400000.5) 
	hour = int(eclipse[3] * 24.0)
	minu = (eclipse[3] * 24.0 - hour) * 60
	print ("eclipse %3d --> %2d %2d %4d - %2d : %2d "  \
	% (i+1,eclipse[2], eclipse[1], eclipse[0],hour,int(minu)))
