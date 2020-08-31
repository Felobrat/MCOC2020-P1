# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:41:12 2020

@author: Felipe Bravo
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
# PARAMETROS

pi = np.pi
m = 1       #[kg]
f = 1       #[Hz]
chi = 0.2
omega = 2*pi*f
k = m*omega**2 
c = 2*chi*omega*m


# z = [x xp]
# zp =[xp, xpp]



def zpunto(z,t):
    zp = np.zeros(2)
    zp[0] = z[1]
    zp[1] = -(c*z[1]+k*z[0])/m
    return zp

def eulerint(zp,z0,t,Nsub=1):
    Nt=len(t)
    Ndim=len(z0)
    z = np.zeros((Nt, Ndim))
    z[0,:] = z0
    
    for i in range(1,Nt):
        t_ant = t[i-1]
        dt=(t[i]-t[i-1])/Nsub
        z_temp = z[i-1,:]*1.0
        for k in range(Nsub):
            z_temp+=dt*zp(z_temp,t_ant+k*dt)
        z[i,:]=z_temp
    return z
    

z0 = np.array([0,1])
t = np.linspace(0,4.,1000)
ureal = 1/omega*np.exp(-c/(2*m)*t)*np.sin(omega*t)
sol=odeint(zpunto,z0,t)
soleu0 = eulerint(zpunto,z0,t,1)
soleu10 = eulerint(zpunto,z0,t,10)
soleu100 = eulerint(zpunto,z0,t,100)
x = sol[:,0]
x0= soleu0[:,0]
x1=soleu10[:,0]
x2=soleu100[:,0]
fig = plt.figure(figsize=(7,7)) 
ax = fig.add_subplot(1,1,1)
ax.plot(t,ureal, color="k", linewidth=2, label="u real")
ax.plot(t,x, color="b", linewidth=2, label="u odeint")
ax.plot(t,x0, color="green", linestyle="--", linewidth=2, label="u Euler N=1")
ax.plot(t,x1, color="red", linestyle="--", linewidth=2, label="u Euler N=10")
ax.plot(t,x2, color="orange", linestyle="--", linewidth=2, label="u Euler N=100")
ax.set_xlabel("tiempo [s]")
ax.set_ylabel("u [m]")

plt.legend()
plt.show()

fig.savefig('EDM_DUMPU_homogenea.png')