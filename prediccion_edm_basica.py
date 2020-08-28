# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 18:41:43 2020

@author: Felipe Bravo
"""

import scipy as sp
import numpy as np
from scipy.integrate import odeint
from numpy import linalg as LA
from datetime import datetime

#UNIDADES
km = 1000 #metros
cm = 0.01 #metro
inch = 2.54*cm 
h = 3600. #s
pi = np.pi

#variables FISICAS
Mt = 5.972*10**24        # Masa de la tierra [kg]
Rt = 6371.*km            # Radio de la tierra [m]
omegat=-2*pi/86400.        # Rapidez angular de la tierra [rad/s]
G = 6.673*10**(-11)      # Constate Gravitacional [N*m**2/kg**2]
ms = 2300.               # masa satelite [kg]
HG = 700*km              # Altura de orbita del satelite

FGmax = -G*Mt/Rt**2     #Fuerza gravitacional maxima en la superficie de la tierra

zp=np.zeros(6)

def zpunto(z,t):
    c = np.cos(omegat*t)
    s = np.sin(omegat*t)
    
    R   = np.array([[ c, s,0],[-s, c,0],[0,0,1]])
    Rp  = omegat*np.array([[-s, c,0],[-c,-s,0],[0,0,0]])
    Rpp = omegat**2*np.array([[-c,-s,0],[ s,-c,0],[0,0,0]])
    
    z1 = z[0:3]
    z2 = z[3:6]
    
    r2 = np.dot(z1,z1)
    r  = np.sqrt(r2)
    Fg = (-G*Mt/Rt**2)*(R@(z1/r))
    
    zp[0:3] = z2
    zp[3:6] = R.T@(Fg - (2*(Rp@z2) + (Rpp@z1)))
    return zp


# TIEMPOS INICIALES

      #<TAI>TAI=2020-08-05T23:00:19.000000</TAI>
      #<UTC>UTC=2020-08-05T22:59:42.000000</UTC>
      #<UT1>UT1=2020-08-05T22:59:41.796108</UT1>
      #<Absolute_Orbit>+33778</Absolute_Orbit>

            
ti="2020-08-05T22:59:42.000000"
ti=ti.split("T")
ti="{} {}".format(ti[0],ti[1])
ti=datetime.strptime(ti,"%Y-%m-%d %H:%M:%S.%f")

# TIEMPOS FINALES

      #<TAI>TAI=2020-08-07T01:00:19.000000</TAI>
      #<UTC>UTC=2020-08-07T00:59:42.000000</UTC>
      #<UT1>UT1=2020-08-07T00:59:41.796594</UT1>
      #<Absolute_Orbit>+33793</Absolute_Orbit>

tf="2020-08-07T00:59:42.000000"
tf=tf.split("T")
tf="{} {}".format(tf[0],tf[1])
tf=datetime.strptime(tf,"%Y-%m-%d %H:%M:%S.%f")

deltaT=(tf-ti).seconds

# DATOS INICIALES

      #<X unit="m">792018.848343</X>
      #<Y unit="m">-5352421.287231</Y>
      #<Z unit="m">4551102.773278</Z>
      #<VX unit="m/s">-2467.932652</VX>
      #<VY unit="m/s">4443.919325</VY>
      #<VZ unit="m/s">5640.241199</VZ>

# POSICION 
x_i=792018.848343        #[m]
y_i=-5352421.287231      #[m]
z_i=4551102.773278       #[m]

# VELOCIDAD
vx_i=-2467.932652        #[m/s]
vy_i=4443.919325         #[m/s]
vz_i=5640.241199          #[m/s]

# DATOS FINALES

      #<X unit="m">-1000730.204814</X>
      #<Y unit="m">-5946686.700579</Y>
      #<Z unit="m">-3710650.426014</Z>
      #<VX unit="m/s">-2311.357909</VX>
      #<VY unit="m/s">-3536.726521</VY>
      #<VZ unit="m/s">6300.485353</VZ>
      
# POSICION  
x_f=-1000730.204814      #[m]
y_f=-5946686.700579      #[m]
z_f=-3710650.426014      #[m]

# VELOCIDAD
vx_f=-2311.357909        #[m/s]
vy_f=-3536.726521        #[m/s]
vz_f=6300.485353         #[m/s]


t=np.linspace(0,deltaT,1000000)
z0=np.array([x_i,y_i,z_i,vx_i,vy_i,vz_i])
sol=odeint(zpunto,z0,t)
x=sol[:,:]

pos_final=np.array([x_f,y_f,z_f,vx_f,vy_f,vz_f])-sol[-1]

dif_dist = np.array([pos_final[0],pos_final[1],pos_final[2]])

dist = LA.norm(dif_dist)
print(f"{dist} [m]")

