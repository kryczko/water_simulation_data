#!/usr/bin/env python

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata

nooa = 64

f = open('hbonds_histogram.dat', 'r')

xcoords = []
ycoords = []
zcoords = []
hbonds = []

for line in f.readlines():
    x = float(line.split()[0])
    y = float(line.split()[1])
    z = float(line.split()[2])
    hb = float(line.split()[3])

    xcoords.append(x)
    ycoords.append(y)
    zcoords.append(z)
    hbonds.append(hb)	

plt.figure()

xlist = xcoords
ylist = ycoords               #np.linspace(-2., 1., 100)
zlist = zcoords               #np.linspace(-1., 1., 100)
hblist = hbonds

yi = np.linspace(0., 12.42, 100)
zi = np.linspace(0., 12.42, 100)

hbi = griddata (ylist, zlist, hblist, yi, zi)

CS = plt.contour(yi,zi,hbi,15,linewidths=0.5,colors='k')
CS = plt.contourf(yi,zi,hbi,15,cmap=plt.cm.rainbow)
plt.colorbar()


plt.title('Number of hydrogen bonds per molecule wrt y and z axis')
plt.xlabel('y-axis [Angstroms]')
plt.ylabel('z-axis [Angstroms]')

#plt.clabel(CS, inline=1, fontsize=10)
plt.show()

