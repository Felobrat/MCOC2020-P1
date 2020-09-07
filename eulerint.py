# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 21:55:50 2020

@author: Felipe Bravo
"""
import numpy as np

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