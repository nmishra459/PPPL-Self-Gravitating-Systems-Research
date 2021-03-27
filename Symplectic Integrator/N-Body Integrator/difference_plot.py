#!/usr/bin/python
import numpy as np
from matplotlib import pylab
from pylab import *
from mpl_toolkits import mplot3d
import subprocess

num_bodies = 3

ax = plt.plot()

##ax.set_xlabel('x (meters)')
##ax.set_ylabel('y (meters)')

#ax.grid() 

rcParams["axes.grid"] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 18

with open("time.txt") as f:
    time = [double(d) for d in f]

distance_difference = [None] * len(time)


figure()
for n in range(1,num_bodies+1):

    x_string = "x" + str(n) + ".txt"
    y_string = "y" + str(n) + ".txt"
    z_string = "z" + str(n) + ".txt" 

    with open(x_string) as f:
        x_array = [double(a) for a in f]
    
    with open(y_string) as f:
        y_array = [double(a) for a in f]

    with open(z_string) as f:
        z_array = [double(a) for a in f]
    
    ghost_x_string = "gx" + str(n) + ".txt"
    ghost_y_string = "gy" + str(n) + ".txt"
    ghost_z_string = "gz" + str(n) + ".txt" 

    with open(ghost_x_string) as f:
        ghost_x_array = [double(a) for a in f]
    
    with open(ghost_y_string) as f:
        ghost_y_array = [double(a) for a in f]

    with open(ghost_z_string) as f:
        ghost_z_array = [double(a) for a in f]


for b in range(0, len(time)):
    distance_difference[b] = sqrt((ghost_x_array[b]-x_array[b])**2 + (ghost_y_array[b]-y_array[b])**2 + (ghost_z_array[b]-z_array[b])**2) 

plot (time, distance_difference)

tight_layout()




show()

print("Graphing Code (Python) Termination: Complete\n")

