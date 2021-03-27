from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import numpy as np
from matplotlib import pylab
from pylab import *
from mpl_toolkits import mplot3d
import subprocess


num_bodies = 8
start_value = 1

np.set_printoptions(threshold=np.inf)
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')

rcParams["axes.grid"] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 18


with open("time.txt") as f:
    time = [double(d) for d in f]

my_line = ax1.plot([],[],[])
xdata , ydata, zdata = [], [], []

def init():
    global x_array,y_array,z_array,E_array
    global xdata, ydata, zdata
    #plots[num_bodies]

    for n in range(start_value,num_bodies+1):
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

    my_line = ax1.plot([],[],[])
    xdata , ydata, zdata = [], [], []



def animate(i):  
    for n in range(start_value,num_bodies+1):
        if (n%1 == 0):
            #xdata.append(x_array[:i])
            #ydata.append(y_array[:i])
            #zdata.append(z_array[:i])
 
            #my_line.set_data(xdata, ydata, zdata)

            ax1.plot3D(x_array[:i], y_array[:i], z_array[:i], color = 'green')
            tight_layout()


anim = animation.FuncAnimation(fig, animate, init_func=init, interval=.001, save_count=10000000)

ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_ylabel('z')

anim.save('trasic_animation8.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

#plt.show()