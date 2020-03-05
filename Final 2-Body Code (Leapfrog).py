#!/usr/bin/python

from matplotlib import pylab
from pylab import *

###############################################################################
N = 10000
G = 1
t = linspace(0,10,N)
dt = t[1] - t[0]

###############################################################################
# functions

def Force(mass1, mass2, radius, coord1, coord2):
    return (G*mass1*mass2)/(radius**2)*(coord2-coord1)/(radius)

def Energy(mass1, mass2, xvelocity, yvelocity, radius):
    return (1/2)*mass1*(xvelocity**2+yvelocity**2) + (-G*mass1*mass2)/radius

def radius(coord1x, coord1y, coord2x, coord2y):
    return ((coord2y - coord1y)**2 + (coord2x-coord1x)**2)**(1/2)

def integrate(F1,F2,m1,m2,x10,y10,vx10,vy10,x20,y20,vx20,vy20):
    ###########################################################################
    # arrays are allocated and filled with zeros
    r12 = zeros(N)

    x1 = zeros(N)
    y1 = zeros(N)

    x2 = zeros(N)
    y2 = zeros(N)
 
    vx1 = zeros(N)
    vx2 = zeros(N)

    vy1 = zeros(N)
    vy2 = zeros(N)
    
    E1 = zeros(N)    
    E2 = zeros(N)

    ###########################################################################    
    # initial conditions          
    x1[0] = x10                   
    y1[0] = y10                   

    x2[0] = x20
    y2[0] = y20

    r12[0] = radius(x1[0], y1[0], x2[0], y2[0])
    #r12[-1] = radius(x1[-1], y1[-1], x2[-1], y2[-1])
        
    vx1[0] = vx10
    vx2[0] = vx20
    
    vy1[0] = vy10
    vy2[0] = vy20

    E1[0] = Energy(m1, m2, vx1[0], vy1[0], r12[0])
    E2[0] = Energy(m2, m1, vx1[0], vy1[0], r12[0])


    ###########################################################################
    # integration

    for i in range(N-1):

         r12[i] = radius(x1[i], y1[i], x2[i], y2[i])
        

         #HALF TIME STEP UPDATE
         x1half = vx1[i] * (dt/2) + x1[i]
         y1half = vy1[i] * (dt/2) + y1[i]
         
         x2half = vx2[i] * (dt/2) + x2[i]
         y2half = vy2[i] * (dt/2) + y2[i]

         rhalf = radius(x1half, y1half, x2half, y2half)
            #print (Force(m1, m2, rhalf, x1half, x2half)/m1)
         vx1[i+1] = vx1[i] + dt*Force(m1, m2, rhalf, x1half, x2half)/m1
         vx2[i+1] = vx2[i] + dt*Force(m1, m2, rhalf, x2half, x1half)/m2

         #print((G*m1*m2/(rhalf)**2)*(x2half-x1half)/(rhalf))

         vy1[i+1] = vy1[i] + dt*Force(m1, m2, rhalf, y2half, y1half)/m1
         vy2[i+1] = vy2[i] + dt*Force(m1, m2, rhalf, y1half, y2half)/m2  

         x1[i+1] = x1half + dt/2 * vx1[i+1]
         y1[i+1] = y1half + dt/2 * vy1[i+1]
    
         x2[i+1] = x2half + dt/2 * vx2[i+1]
         y2[i+1] = y2half + dt/2 * vy2[i+1]

         r12[i+1] = radius(x1[i+1], y1[i+1], x2[i+1], y2[i+1])
                                                                                      
         E1[i+1] += Energy(m1, m2, vx1[i+1], vy1[i+1], r12[i+1])
         E2[i+1] += Energy(m2, m1, vx1[i+1], vy1[i+1], r12[i+1])


    E1[-1] = Energy(m1, m2, vx1[-1], vy1[-1], r12[-1])
    E2[-1] = Energy(m2, m1, vx1[-1], vy1[-1], r12[-1])

    # return solution
    return m1,m2,x1,y1,vx1,vy1,x2,y2,vx2,vy2,E1,E2

###############################################################################
# numerical integration
F1 = zeros(N)
F2 = zeros(N)
m1,m2,x1,y1,vx1,vy1,x2,y2,vx2,vy2,E1,E2 = integrate(  F1,   F2,   1,   1,  -.1,  0,    -.1,     .1,    .1,   0,   .1,    -.1) 
                                                  #   F1,   F2,   m1,  m2, x10,  y10,  vx10,  vy10,  x20,  y20, vx20,  vy20
###############################################################################
rcParams["axes.grid"] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 18

figure()
plot(x1,y1)
plot(x2,y2)
ylabel("y (m)")  
xlabel("x (m)")


figure()
plot(t,E1)
plot(t,E2)
ylabel("Energy (J)")
xlabel("Time (Sec)")

tight_layout()


###############################################################################

show()

###############################################################################




