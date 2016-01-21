import numpy as np
import matplotlib.pyplot as plt

fname='prueba.dat'

#Read the data file
data = np.loadtxt(fname,usecols=(0,1,2),comments='#')
#data is a complicated python object, let us split it
time = []
fase = []
err  = []
for i in range(0,len(data)-1):
	a,b,c=data[i][:]
	time.append(a)
	fase.append(b)
	err.append(c)

plt.plot(time,fase,'ro')
