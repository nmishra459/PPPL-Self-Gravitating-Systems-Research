#!/usr/bin/python
import numpy as np
from matplotlib import pylab
from pylab import *
from mpl_toolkits import mplot3d
import subprocess

subprocess.call(["gcc", "integrate.c"]) ##compiles C code (integration)

print("Integration Code (C) Compilation: Success\n");

tmp=subprocess.call("./b.exe") ##initates code execution

print("Graphing Code (Python) Initiation: Success")

x1 = []
y1 = []
z1 = []

x2 = []
y2 = []
z2 = []

x3 = []
y3 = []
z3 = []

E1 = []
E2 = []
E3 = []
ET = []

time = []

with open("x1.txt") as f:
    x1 = [ double(i) for i in f]
with open("y1.txt") as f:
    y1 = [ double(i) for i in f]
with open("z1.txt") as f:
    z1 = [ double(i) for i in f]

with open("x2.txt") as f:
    x2 = [ double(i) for i in f]
with open("y2.txt") as f:
    y2 = [ double(i) for i in f]
with open("z2.txt") as f:
    z2 = [ double(i) for i in f]

with open("x3.txt") as f:
    x3 = [ double(i) for i in f]
with open("y3.txt") as f:
    y3 = [ double(i) for i in f]
with open("z3.txt") as f:
    z3 = [ double(i) for i in f]

with open("E1.txt") as f:
    E1 = [ double(i) for i in f]
with open("E2.txt") as f:
    E2 = [ double(i) for i in f]
with open("E3.txt") as f:
    E3 = [ double(i) for i in f]
with open("ET.txt") as f:
    ET = [ double(i) for i in f]

with open("time.txt") as f:
    time = [ double(i) for i in f]


ax = plt.axes(projection='3d')

ax.set_xlabel('x (meters)')
ax.set_ylabel('y (meters)')
ax.set_zlabel('z (meters)')

# Data for three-dimensional scattered points
ax.plot3D(x1, y1, z1, 'red')
ax.plot3D(x2, y2, z2, 'green')
ax.plot3D(x3, y3, z3, 'blue')

ax.grid()

rcParams["axes.grid"] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 18

figure()
plot(time, E1)
plot(time, E2)
plot(time, E3)
plot(time, ET)
ylabel("Energy (J)")
xlabel("Time (Sec)")

tight_layout()

print("Graphing Code (Python) Termination: Complete\n")

show()


