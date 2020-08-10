#!/usr/bin/python
import numpy as np
from matplotlib import pylab
from pylab import *
from mpl_toolkits import mplot3d
import subprocess

num_bodies = 3

subprocess.call(["gcc", "-lm", "Ghost-Body-Energy.c"]) ##compiles C code (integration)

print("Integration Code (C) Compilation: Success\n");

tmp=subprocess.call("./c.exe") ##initates code execution

print("Graphing Code (Python) Initiation: Success\n")

ax = plt.axes(projection='3d')

ax.set_xlabel('x (meters)')
ax.set_ylabel('y (meters)')
ax.set_zlabel('z (meters)')

ax.grid() 

rcParams["axes.grid"] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 18


with open("time.txt") as f:
    time = [double(d) for d in f]

figure()
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

tight_layout()
ylabel("Energy (J)")
xlabel("Time (Sec)")

# figure()
# for n in range(1,num_bodies+1):
#     KE_string = "KE" + str(n) + ".txt" 
#     with open(KE_string) as f:
#         KE_array = [double(a) for a in f]
#         plot(time, KE_array)



# tight_layout()
# ylabel("Kinetic Energy (J)")
# xlabel("Time (Sec)")



# figure()
# for n in range(1,num_bodies+1):
#     PE_string = "PE" + str(n) + ".txt" 
#     with open(PE_string) as f:
#         PE_array = [double(a) for a in f]
#         plot(time, PE_array)

# tight_layout()
# ylabel("Potential Energy (J)")
# xlabel("Time (Sec)")



figure()
with open("ET.txt") as f:
    ET = [double(d) for d in f]
# with open("KET.txt") as f:
#     KET = [double(d) for d in f]
# with open("PET.txt") as f:
#     PET = [double(d) for d in f]

plot(time, ET)
##plot(time, KET)
##plot(time, PET)

tight_layout()
ylabel("Total Energy (J)")
xlabel("Time (Sec)")

show()
show()
show()
show()


print("Graphing Code (Python) Termination: Complete\n")

