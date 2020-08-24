# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:48:46 2020

@author: Felipe Bravo
"""
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#UNIDADES
cm = 0.01 #metro
inch = 2.54*cm 

#variables FISICAS
rho = 1.225
Cd = 0.47
D = 8.5*inch
r = D/2
A=sp.pi*r**2
CD = 0.5*rho*Cd*A
g=9.81 #m/s^2
m=15 #kg

#Viento en contra
V1 = 0. #m/s
V2 = 10.
V3 = 20.

#Funcion a Integrar
# Variable estado Z = [x, y , xp, yp]

def bala1(z,t):
    zp=sp.zeros(4)
    zp[0] = z[2]
    zp[1] = z[3]
    v = z[2:4]
    v[0] = v[0] - V1 
    v2 = sp.dot(v,v)
    vnorm = sp.sqrt(sp.dot(v,v))
    FD = -CD*v2*(v/vnorm)
    zp[2] = FD[0]/m
    zp[3] = FD[1]/m -g
    
    return zp

def bala2(z,t):
    zp=sp.zeros(4)
    zp[0] = z[2]
    zp[1] = z[3]
    v = z[2:4]
    v[0] = v[0] - V2 
    v2 = sp.dot(v,v)
    vnorm = sp.sqrt(sp.dot(v,v))
    FD = -CD*v2*(v/vnorm)
    zp[2] = FD[0]/m
    zp[3] = FD[1]/m -g
    
    return zp

def bala3(z,t):
    zp=sp.zeros(4)
    zp[0] = z[2]
    zp[1] = z[3]
    v = z[2:4]
    v[0] = v[0] - V3 
    v2 = sp.dot(v,v)
    vnorm = sp.sqrt(sp.dot(v,v))
    FD = -CD*v2*(v/vnorm)
    zp[2] = FD[0]/m
    zp[3] = FD[1]/m -g
    
    return zp

t = sp.linspace(0,5.5,1001)
#condiciones iniciales 
vi = 100*1000./3600.
z0 = sp.array([0,0,vi, vi])
sol1 = odeint(bala1, z0, t)
sol2 = odeint(bala2, z0, t)
sol3 = odeint(bala3, z0, t)

x1= sol1[:,0]
y1 = sol1[:,1]

x2= sol2[:,0]
y2 = sol2[:,1]

x3= sol3[:,0]
y3 = sol3[:,1]

fig = plt.figure(1)
plt.plot(x1,y1, label='V=0 m/s')
plt.plot(x2,y2, label='V=10.0 m/s')
plt.plot(x3,y3, label='V=20.0 m/s')

plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.title('Trayectoria para distintos vientos')
plt.axis([0, 150, 0, 50])
plt.yticks([0, 10, 20, 30, 40, 50])
plt.legend()
plt.grid()
plt.show()
plt.tight_layout()
fig.savefig('balistica.png')