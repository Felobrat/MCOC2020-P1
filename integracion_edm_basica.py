# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 08:09:41 2020

@author: Felipe Bravo
"""

import scipy as sp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#UNIDADES
km = 1000 #metros
cm = 0.01 #metro
inch = 2.54*cm 
pi = np.pi

#variables FISICAS
Mt = 5.972*10**24        # Masa de la tierra [kg]
Rt = 6371.*km            # Radio de la tierra [m]
omegat=2*pi/86400.        # Rapidez angular de la tierra [rad/s]
G = 6.673*10**(-11)      # Constate Gravitacional [N*m**2/kg**2]
ms = 2300.               # masa satelite [kg]
print(f"Masa de la tierra = {Mt} [kg]")
print(f"Radio de la tierra = {Rt} [m]")
print(f"Vel. Ang. de la tierra = {omegat} [rad/s] ")
print(f"Masa del Satelite = {ms} [kg]")
print(f"Constante Gravitacional = {G} [Nm^2 / kg^2]")
# por simplicidad del problema, se usan coordenadas cilindricas con vectores unitarios (er,etheta,ez)
# y coordenadas cilindricas con coordenadas  (r, theta , z)
# con derivadas temporales de primer orden (rt, thetat, zt) 
# con derivadas temporales de segundo orden (rtt, thetatt, ztt) 

# el cambio de coordenadas cilindricas a cartesianas es de la forma: 
    # [  cos(theta) sen(theta) 0] [ex]   [er]
    # [- sen(theta) cos(theta) 0] [ey] = [etheta]
    # [      0          0      1] [ez]   [ez]
    
# la ecuacion queda de forma 
# ms*a(er,etheta,ez)=FG(-er,0,0)

# ms*[(rtt-r*thetat**2)er+(2*rt*thetat+r*thetatt)*etheta+ztt*ez]=-G*mt*ms/(r**2)(-er)
# tenemos 3 ecuaciones:
    # (1) er:        ms*((rtt-r*thetat**2)) = G*Mt*ms/(r**2)
    # (2) etheta:       ms*(2*rt*thetat+r*thetatt) = 0
    # (3) ez:                        ms*ztt = 0

# z = [r, theta , rt, thetat]
# zt = [rt , thetat, rtt, thetatt]
def sentinel1A(z,t):
    zt=np.zeros(4)
    zt[0] = z[2]
    zt[1] = z[3]
    zt[2] = -G*Mt/((z[0])*(z[0]))+(z[0])*(z[3])*(z[3])
    zt[3] = -2*z[2]*(z[3])/(z[0])
    
    return zt

t = np.linspace(0,86400*2,1001)
#condiciones iniciales

#Dado que se lanza desde 80km sobre la superficie de la tierra,
# Usando la ecuacion asociada a la direccion radial en la superficie de la tierra
#con una aceleracion radial nula, se tiene
# r*thetat**2 = G*Mt*ms/(r**2)
# despejando thetat
# thetat = sqrt(G*Mt*ms/(r**3))
# r0 = 6371 [km]
# Mt = 5.972e+24 [kg]
# se obtiene que thetat0 = 1.241e-3 [rad/s]
z0 = np.array([Rt+80*km,0,0,1.241*10**(-3)])
sol = odeint(sentinel1A, z0, t)
x1= sol[:,0]
y1 = sol[:,1]+omegat*t
y2 = sol[:,1]
print(f"H. max. del satelite lanzado con vel ang de la tierra = {max(sol[:,0])-Rt} [m]")

# se observa que a esa altura, se lleva a 6953531.685610919
# si se quisiera que esté a 700km sobre la superficie, se puede dar velocidad adicional para que llegue a 700 km.
# haciendo el calculo, la velocodad maxima para que se mantenga entre 80km y 700km sobre es 1.241e-3 [rad/s]
z1 = np.array([Rt+80*km,0,0,1.246*10**(-3)])
sol2 = odeint(sentinel1A, z1, t)
x3 = sol2[:,0]
y3 = sol2[:,1]+omegat*t
print(f"H. max. del satelite con vel ang mayor que de la tierra = {max(sol2[:,0])-Rt}")
theta = np.linspace(0,2*pi,1001)
theta2 = np.linspace(0,60*pi,1001)

radiotierra = np.linspace(Rt,Rt,1001)
atmosfera = np.linspace(Rt+1000*km, Rt+1000*km, 1001)
lininfsat = np.linspace(Rt+80*km, Rt+80*km, 1001)
linsupsat = np.linspace(Rt+700*km, Rt+700*km, 1001)


fig = plt.figure(figsize=(15,15))  
plt.subplot
ax = fig.add_subplot(2,3,2, projection="polar")
ax.plot(y1,x1,color="red",linewidth=1, label='SATELITE')
ax.plot(theta,radiotierra,color="green",linewidth=1, label='Superficie de la Tierra')
ax.plot(theta,atmosfera,color="blue",linewidth=1, label='Límite de la atmosfera')
ax.plot(theta,lininfsat,color="purple",linewidth=1, label='Limite inferior satelite')
ax.plot(theta,linsupsat,color="purple",linewidth=1, label='Limite superior satelite')
ax.legend()

ax.set_title('Movimiento Relativo a la tierra')
graf = fig.add_subplot(2,3,5)
graf.set_yscale("log")
graf.plot(y1,x1,color="red",linewidth=1, label='SATELITE')
graf.plot(theta2,radiotierra,color="green",linewidth=1, label='Superficie de la Tierra')
graf.plot(theta2,atmosfera,color="blue",linewidth=1, label='Límite de la atmosfera')
graf.plot(theta2,lininfsat,color="purple",linewidth=1, label='Limite inferior satelite')
graf.plot(theta2,linsupsat,color="purple",linewidth=1, label='Limite superior satelite')
graf.set_xlim([0, 100])
graf.set_xlabel("ANGULO THETHA [rad]")
graf.set_ylabel("ALTURA [m]")

# ejes del satelite
ax2 = fig.add_subplot(2,3,1, projection="polar")
ax2.plot(y2,x1,color="red",linewidth=1, label='SATELITE')
ax2.plot(theta,radiotierra,color="green",linewidth=1, label='Superficie de la Tierra')
ax2.plot(theta,atmosfera,color="blue",linewidth=1, label='Límite de la atmosfera')
ax2.plot(theta,lininfsat,color="purple",linewidth=1, label='Limite inferior satelite')
ax2.plot(theta,linsupsat,color="purple",linewidth=1, label='Limite superior satelite')
ax2.legend()
ax2.set_title('Movimiento EJES DEL SATELITE')


graf2 = fig.add_subplot(2,3,4)
graf2.set_yscale("log")
graf2.plot(y3,x1,color="red",linewidth=1, label='SATELITE')
graf2.plot(theta2,radiotierra,color="green",linewidth=1, label='Superficie de la Tierra')
graf2.plot(theta2,atmosfera,color="blue",linewidth=1, label='Límite de la atmosfera')
graf2.plot(theta2,lininfsat,color="purple",linewidth=1, label='Limite inferior satelite')
graf2.plot(theta2,linsupsat,color="purple",linewidth=1, label='Limite superior satelite')
graf2.set_xlim([0, 100])
graf2.set_xlabel("ANGULO THETHA [rad]")
graf2.set_ylabel("ALTURA [m]")

# hasta 700 km
ax3 = fig.add_subplot(2,3,3, projection="polar")
ax3.plot(y3,x3,color="red",linewidth=1, label='SATELITE')
ax3.plot(theta,radiotierra,color="green",linewidth=1, label='Superficie de la Tierra')
ax3.plot(theta,atmosfera,color="blue",linewidth=1, label='Límite de la atmosfera')
ax3.plot(theta,lininfsat,color="purple",linewidth=1, label='Limite inferior satelite')
ax3.plot(theta,linsupsat,color="purple",linewidth=1, label='Limite superior satelite')
ax3.legend()
ax3.set_title('Movimiento relativo a la tierra hasta 700 km')


graf3 = fig.add_subplot(2,3,6)
graf3.set_yscale("log")
graf3.plot(y3,x3,color="red",linewidth=1, label='SATELITE')
graf3.plot(theta2,radiotierra,color="green",linewidth=1, label='Superficie de la Tierra')
graf3.plot(theta2,atmosfera,color="blue",linewidth=1, label='Límite de la atmosfera')
graf3.plot(theta2,lininfsat,color="purple",linewidth=1, label='Limite inferior satelite')
graf3.plot(theta2,linsupsat,color="purple",linewidth=1, label='Limite superior satelite')
graf3.set_xlim([0, 100])
graf3.set_xlabel("ANGULO THETHA [rad]")
graf3.set_ylabel("ALTURA [m]")

plt.tight_layout()
fig.savefig('SATELITE.png')