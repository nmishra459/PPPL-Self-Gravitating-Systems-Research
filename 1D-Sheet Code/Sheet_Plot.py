#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import subprocess

num_bodies = 10000

subprocess.call(["gcc", "-lm", "Sheet_N_Body.c"]) ##compiles C code (integration)
print("Integration Code (C) Compilation Completed\n")
tmp=subprocess.call("./a.exe") ##initates code execution

x = []
figure()
for n in range(1,num_bodies+1):

    x_string = "x" + str(n) + ".txt"

    v_string = "v" + str(n) + ".txt" 

    with open(x_string) as f:
        x_array = [double(a) for a in f]
        x.append(x_array)
    with open(v_string) as f:
        v_array = [double(a) for a in f]
    
    if (n%1 == 0):
        plt.plot(x_array, v_array, 'o', color = 'black')

# merged_list = []

# for l in x:
#     merged_list += l

# hist, bins = np.histogram(merged_list, bins=100, normed=True)
# bin_centers = (bins[1:]+bins[:-1])*0.5
# plt.plot(bin_centers, hist)



# #plt.hist(merged_list, bins=100)

tight_layout()
xlabel('x (meters)')
ylabel('v (meters/seconds)')
rcParams["axes.grid"] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 18

show()

print("Graphing Code (Python) Completed\n")

