#!/usr/bin/python
import numpy as np
from matplotlib import pylab
from pylab import *

ax.set_xlabel('Vertices (n)')
ax.set_ylabel('Avg. Dist.')

ax.grid() 

rcParams["axes.grid"] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 18


with open("time.txt") as f:
    time = [double(d) for d in f]

figure()

n = [128,256,512,1024,2048,4096,8192,16384,32768,65536,131072]
dist = [7.581495, 10.678034, 14.967085, 21.056606, 29.551414, 41.910149,58.824887,83.197649,117.469283,166.030158,234.722987]

tight_layout()
ylabel("Energy (J)")
xlabel("Time (Sec)")

plot(n, dist)
##plot(time, KET)
##plot(time, PET)

show()


print("Graphing Code (Python) Termination: Complete\n")

