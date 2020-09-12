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
import xml.etree.ElementTree as ET
from numpy import zeros
from scipy import special

pi=np.pi
km = 1000
Rt = 6378.1363*km 
km = 1000
h = 3600
ncoef = 8
G = 6.673*10**(-11)      # Constate Gravitacional [m**3/kg s**2]
Mt = 5.972*10**24        # Masa de la tierra [kg]    

""""FUNCIONES
"""

def utc2time(utc, ut1, EOF_datetime_format = "%Y-%m-%dT%H:%M:%S.%f"):
	t1 = dt.datetime.strptime(ut1,EOF_datetime_format)
	t2 = dt.datetime.strptime(utc,EOF_datetime_format)
	return (t2 - t1).total_seconds()


def leer_eof(fname):
	tree = ET.parse(fname)
	root = tree.getroot()

	Data_Block = root.find("Data_Block")		
	List_of_OSVs = Data_Block.find("List_of_OSVs")

	count = int(List_of_OSVs.attrib["count"])

	t = zeros(count)
	x = zeros(count)
	y = zeros(count)
	z = zeros(count)
	vx = zeros(count)
	vy = zeros(count)
	vz = zeros(count)

	set_ut1 = False
	for i, osv in enumerate(List_of_OSVs):
		UTC = osv.find("UTC").text[4:]
		
		x[i] = osv.find("X").text   #conversion de string a double es implicita
		y[i] = osv.find("Y").text
		z[i] = osv.find("Z").text
		vx[i] = osv.find("VX").text
		vy[i] = osv.find("VY").text
		vz[i] = osv.find("VZ").text

		if not set_ut1:
			ut1 = UTC
			set_ut1 = True

		t[i] = utc2time(UTC, ut1)

	return t, x, y, z, vx, vy, vz


def potencial(x, y, z, ncoef=ncoef):
    r = np.sqrt(x**2 + y**2 + z**2)
    θ = np.arctan2(z, np.sqrt(x**2 + y**2))
    φ = np.arctan2(y, x)
    
    # calcula el potencial
    U, dU = special.lpmn(ncoef+1,ncoef+1, np.sin(θ))
    
    u = 1 / r
    for n in range(0,ncoef):
        for m in range(n):
            u += U[m,n] * np.sin(m*φ) / r**(n+1)
            u += U[m,n] * np.cos(m*φ) / r**(n+1)
    
    return G*Mt * u
    

def fuerza(Z, t, ncoef=ncoef, omega= 2*pi/86200., dx=1):
    """Calcula la fuerza que siente un cuerpo
    
    coef: coeficientes de Legendre
    """
    x, y, z = Z[0], Z[1], Z[2]
    vx, vy, vz = Z[3], Z[4], Z[5]
    
    c=np.cos(omega*t)
    s=np.sin(omega*t)
    
    R=np.array([[c,-s,0],
                [s,c,0],
                [0,0,1]])
    Rp=omega*np.array([[-s,-c,0],
                       [c,-s,0],
                       [0,0,0]])
    Rpp=(omega**2)*np.array([[-c,s,0],
                             [-s,-c,0],
                             [0,0,0]])
    
    # TODO: calcula la fuerza a partir del potencial.
    fx = 0.5 * (potencial(x+dx,y,z,ncoef) - potencial(x-dx,y,z,ncoef) ) / dx
    fy = 0.5 * (potencial(x,y+dx,z,ncoef) - potencial(x,y-dx,z,ncoef) ) / dx
    fz = 0.5 * (potencial(x,y,z+dx,ncoef) - potencial(x,y,z-dx,ncoef) ) / dx
    pos = np.array([x,y,z])
    vel = np.array([vx,vy,vz])
    
    f = np.array([fx, fy, fz]) - R.T@ (2*Rp @ vel + Rpp @ pos)

    return np.concatenate([[vx, vy, vz], f])



#------------------------------------------------------------------------------------------------------
#                                   CODIGO
#----------------------------------------------------------------------------------------------------
t ,x ,y ,z ,vx ,vy ,vz = leer_eof("1.EOF")
pos = np.array([x,y,z])
print (t[-1])
z_i= np.array([x[0], y[0], z[0], vx[0], vy[0], vz[0]]) #vector de condicion inicial
z_f=np.array([x[-1], y[-1], z[-1],vx[-1], vy[-1], vz[-1] ]) #vector final

t1 = perf_counter()
sol = odeint(fuerza, z_i, t)
t2 = perf_counter()
dts = t2 - t1 
print (f"Para odeint se demoró: {dts} [s]")
px, py, pz, pvx, pvy, pvz = sol[-1]

DELT = np.sqrt((px-z_f[0])**2+(py-z_f[1])**2+(pz-z_f[2])**2)/km
print (f'deriva maxima {DELT} [km]')
print (f'{px-z_f[0]}, {py-z_f[1]}, {pz-z_f[2]}')
print (f'{pvx-z_f[3]}, {pvy-z_f[4]}, {pvz-z_f[5]}')

#--------------------------------------------------------------------------------------------------------
#                              GRAFICOS
#--------------------------------------------------------------------------------------------------------
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
#n.plot(solJ[:,0]/km,solJ[:,1]/km,solJ[:,2]/km, color="r",linewidth=1, label='PREDICCION J')
plt.title('Posición 3D')
n.set_xlabel("x [km]")
n.set_ylabel("y [km]")
n.set_zlabel("z [km]")
plt.legend()

a = figu.add_subplot(1,3,2)
a.plot(x/km,y/km, color="b",linewidth=1, label='REAL')
a.plot(sol[:,0]/km, sol[:,1]/km, color="orange",linewidth=1, label='PREDICCION')
#a.plot(solJ[:,0]/km, solJ[:,1]/km, color="r",linewidth=1, label='PREDICCION J')
plt.title('Posición X-Y')
a.set_xlabel("x [km]")
a.set_ylabel("y [km]")
plt.legend()

b = figu.add_subplot(1,3,3)
b.plot(x/km,z/km, color="b",linewidth=1, label='REAL')
b.plot(sol[:,0]/km, sol[:,2]/km, color="orange",linewidth=1, label='PREDICCION')
#b.plot(solJ[:,0]/km, solJ[:,2]/km, color="r",linewidth=1, label='PREDICCION J')
plt.title('Posición X-Z')
b.set_xlabel("x [km]")
b.set_ylabel("z [km]")
plt.legend()
plt.tight_layout()
figu.savefig('SATELITE.png')
