#!/usr/bin/python
import numpy as np
from matplotlib import pylab
from pylab import *
from mpl_toolkits import mplot3d
import subprocess

num_bodies = 30

subprocess.call(["gcc", "N-Body Computation.c"]) ##compiles C code (integration)

print("Integration Code (C) Compilation: Success\n");

tmp=subprocess.call("./a.exe") ##initates code execution

print("Graphing Code (Python) Initiation: Success\n")

ax = plt.axes(projection='3d')

ax.set_xlabel('x (meters)')
ax.set_ylabel('y (meters)')
ax.set_zlabel('z (meters)')

ax.grid() 

rcParams["axes.grid"] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 18

##figure()

##with open("ET.txt") as f:b
##    ET = [double(d) for d in f]
with open("time.txt") as f:
    time = [double(d) for d in f]

for n in range(1,num_bodies+1):

    x_string = "x" + str(n) + ".txt"
    y_string = "y" + str(n) + ".txt"
    z_string = "z" + str(n) + ".txt" 
    E_string = "E" + str(n) + ".txt" 

    with open(x_string) as f:
        x_array = [double(a) for a in f]
    
    with open(y_string) as f:
        y_array = [double(a) for a in f]

    with open(z_string) as f:
        z_array = [double(a) for a in f]
    
    with open(E_string) as f:
        E_array = [double(a) for a in f]
    if (n%1 == 0):
        ax.plot3D(x_array, y_array, z_array)
        plot(time, E_array)


##plot(time, ET)

ylabel("Energy (J)")
xlabel("Time (Sec)")

tight_layout()

print("Graphing Code (Python) Termination: Complete\n")

show()
