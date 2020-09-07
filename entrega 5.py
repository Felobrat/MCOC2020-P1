# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 01:04:52 2020

@author: Felipe Bravo
"""

import datetime as dt
from time import perf_counter
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from orbita_predicha import satelite
from leer_eof import leer_eof 
from eulerint import eulerint
from tiempo import utc2time
from orbita_predichaJ import sateliteJ
from numpy import linalg as LA

km = 1000
Rt = 6378.1363*km 
km = 1000
h = 3600
t ,x ,y ,z ,vx ,vy ,vz = leer_eof("S1A_OPER_AUX_POEORB_OPOD_20200826T121220_V20200805T225942_20200807T005942.EOF")

pos = np.array([x,y,z])

z_i= np.array([x[0], y[0], z[0], vx[0], vy[0], vz[0]]) #vector de condicion inicial
z_f=np.array([x[-1], y[-1], z[-1],vx[-1], vy[-1], vz[-1] ]) #vector final

t1 = perf_counter()
sol = odeint(satelite, z_i, t)
t2 = perf_counter()
dt = t2 - t1 
print (f"Para odeint se demoró: {dt}")

t3 = perf_counter()
soleu = eulerint(satelite, z_i, t, 1)
t4 = perf_counter()
dte = t4 - t3
print (f"Para eulerint se demoró: {dte}")
#soleu1 = eulerint(satelite, z_i, t, 10)
#soleu2= eulerint(satelite, z_i, t, 20)
#soleu3 = eulerint(satelite, z_i, t, 50)
#soleu4 = eulerint(satelite, z_i, t, 100)
#t5 = perf_counter()
#soleu5 = eulerint(satelite, z_i, t, 500)
#t6 = perf_counter()
#dte2 = t6 - t5

dif = sol - soleu
#dif1 = sol - soleu1
#dif2 = sol - soleu2
#dif3 = sol - soleu3
#dif4 = sol - soleu4
#dif5 = sol - soleu5
posder = np.array(dif)
#posder1 = np.array(dif1)
#posder2 = np.array(dif2)
#posder3 = np.array(dif3)
#posder4 = np.array(dif4)
#posder5 = np.array(dif5)
#net = np.sqrt(np.sum(sol[:,:3]**2,axis=1))
deriva = np.sqrt(np.sum(posder[:,:3]**2,axis=1))
#deriva1 = np.sqrt(np.sum(posder1[:,:3]**2,axis=1))
#deriva2 = np.sqrt(np.sum(posder2[:,:3]**2,axis=1))
#deriva3 = np.sqrt(np.sum(posder3[:,:3]**2,axis=1))
#deriva4 = np.sqrt(np.sum(posder4[:,:3]**2,axis=1))
#deriva5 = np.sqrt(np.sum(posder5[:,:3]**2,axis=1))
DELT = max(deriva)/km
#DELT1 = max(deriva1)/km
#DELT2 = max(deriva2)/km
#DELT3 = max(deriva3/net)
#DELT4 = max(deriva4/net)
#DELT5 = max(deriva5/net*100)
#print (f'deriva maxima {DELT}')




fig = plt.figure(figsize=(7,7)) 

ax = fig.add_subplot(3,1,1)
ax.plot(t/h,x/km, linestyle= '--',label='REAL')
ax.plot(t/h,sol[:,0]/km ,label='PREDICHO')
ax.set_ylabel("x [km]")
plt.title('Posición')

ay = fig.add_subplot(3,1,2)
ay.plot(t/h,y/km, linestyle= '--', label='REAL')
ay.plot(t/h,sol[:,1]/km, label='PREDICHO')
ay.set_ylabel("y [km]")

az = fig.add_subplot(3,1,3)
az.plot(t/h,z/km, linestyle= '--', label='REAL')
az.plot(t/h,sol[:,2]/km, label='PREDICHO')
az.set_xlabel("Tiempo, t [horas]")
az.set_ylabel("z [km]")
plt.legend()
plt.show()
fig.savefig('SATELITE1.png')


fig1 = plt.figure(figsize=(7,7)) 
graf = fig1.add_subplot(1,1,1)
graf.plot(t/h,deriva/km, linestyle= '-',label='REAL')
graf.set_xlabel("Tiempo, t [horas]")
graf.set_ylabel("Deriva [km]")
plt.title('Deriva entre posiciones predichas por ODEINT y EULERINT, $\delta_{max}$' f'= {round(DELT)} [km]')
fig1.savefig('SATELITE2.png')
#fig2 = plt.figure(figsize=(7,7)) 
#graf = fig2.add_subplot(1,1,1)
#graf.plot(t/h,deriva/net*100, linestyle= '-',label='1')
#graf.plot(t/h,deriva1/net*100, linestyle= '-',label='10')
##graf.plot(t/h,deriva2/net*100, linestyle= '-',label='20')
#graf.plot(t/h,deriva3/net*100, linestyle= '-',label='50')
#graf.plot(t/h,deriva4/net*100, linestyle= '-',label='100')
#graf.plot(t/h,deriva5/net*100, linestyle= '-',label='500')
#graf.set_xlabel("Tiempo, t [horas]")
#graf.set_ylabel("%")
#plt.title('Deriva entre posiciones predichas por ODEINT y EULERINT, $\delta_{max}$' f'={round(DELT5)} [%]')
#plt.legend()

sol = odeint(satelite, z_i, t)

pos_final=np.array([x[-1], y[-1], z[-1],vx[-1], vy[-1], vz[-1] ])-sol[-1]


dif_dist = np.array([pos_final[0],pos_final[1],pos_final[2]])
dist = LA.norm(dif_dist)
print(f"Diferencia con Modelo {round(dist/km)} [m]")

t7=perf_counter()
solJ = odeint(sateliteJ, z_i, t)
t8=perf_counter()
tj = t8-t7

print (f'Para odeint con J2 y J3 se demoro: {tj}')
pos_finalJ=np.array([x[-1], y[-1], z[-1],vx[-1], vy[-1], vz[-1] ])-solJ[-1]


dif_dist = np.array([pos_finalJ[0],pos_finalJ[1],pos_finalJ[2]])
dist = LA.norm(dif_dist)
print(f"Diferencia con modelo Gepotencial = {round(dist/km)} [m]")


radiotierra = np.linspace(Rt,Rt,1001)
figu = plt.figure(figsize=(15,5)) 
n = figu.add_subplot(1,3,1, projection="3d", adjustable='box')
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x2 = Rt*np.cos(u)*np.sin(v)
y2 = Rt*np.sin(u)*np.sin(v)
z2 = Rt*np.cos(v)
#n.plot_wireframe(x2, y2, z2, color="b")
n.plot(x/km,y/km,z/km, color="blue",linewidth=1, label='REAL')
n.plot(sol[:,0]/km,sol[:,1]/km,sol[:,2]/km, color="orange",linewidth=1, label='PREDICCION')
n.plot(solJ[:,0]/km,solJ[:,1]/km,solJ[:,2]/km, color="r",linewidth=1, label='PREDICCION J')
plt.title('Posición 3D')
n.set_xlabel("x [km]")
n.set_ylabel("y [km]")
n.set_zlabel("z [km]")
plt.legend()

a = figu.add_subplot(1,3,2)
a.plot(x/km,y/km, color="b",linewidth=1, label='REAL')
a.plot(sol[:,0]/km, sol[:,1]/km, color="orange",linewidth=1, label='PREDICCION')
a.plot(solJ[:,0]/km, solJ[:,1]/km, color="r",linewidth=1, label='PREDICCION J')
plt.title('Posición X-Y')
a.set_xlabel("x [km]")
a.set_ylabel("y [km]")
plt.legend()

b = figu.add_subplot(1,3,3)
b.plot(x/km,z/km, color="b",linewidth=1, label='REAL')
b.plot(sol[:,0]/km, sol[:,2]/km, color="orange",linewidth=1, label='PREDICCION')
b.plot(solJ[:,0]/km, solJ[:,2]/km, color="r",linewidth=1, label='PREDICCION J')
plt.title('Posición X-Z')
b.set_xlabel("x [km]")
b.set_ylabel("z [km]")
plt.legend()
figu.savefig('SATELITE4.png')


fig3 = plt.figure(figsize=(7,7)) 

ax3 = fig3.add_subplot(3,1,1)
ax3.plot(t/h,x/km, linestyle= '--',label='REAL')
ax3.plot(t/h,solJ[:,0]/km ,label='PREDICCION')
ax3.set_ylabel("x [km]")
plt.title('Posición')

ay3 = fig3.add_subplot(3,1,2)
ay3.plot(t/h,y/km, linestyle= '--', label='REAL')
ay3.plot(t/h,solJ[:,1]/km, label='PREDICCION')
ay3.set_ylabel("y [km]")

az3 = fig3.add_subplot(3,1,3)
az3.plot(t/h,z/km, linestyle= '--', label='REAL')
az3.plot(t/h,solJ[:,2]/km, label='PREDICCION')
az3.set_xlabel("Tiempo, t [horas]")
az3.set_ylabel("z [km]")

plt.legend()
fig3.savefig('SATELITE5.png')
t6= perf_counter()
tf=t6-t1
print (f'Tiempo total = {tf}')
plt.tight_layout()