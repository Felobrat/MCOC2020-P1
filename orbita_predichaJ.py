import numpy as np

#Sistema metrico
km = 1000 #metros
cm = 0.01 #metro
inch = 2.54*cm 
h = 3600. #s
pi = np.pi

#variables FISICAS

Rt = 6371.*km            # Radio de la tierra [m]
G = 6.673*10**(-11)      # Constate Gravitacional [m**3/kg s**2]
Mt = 5.972*10**24        # Masa de la tierra [kg]
omega= 2*pi/86200. 
ms = 2300.               # masa satelite [kg]

J2 = 1.75553*10**10*km**5   #[m**5*s**-2]
J3 = -2.6193*10**11*km**6  #[m**6*s**-2]

def sateliteJ(z,t):
    
    c=np.cos(omega*t)
    s=np.sin(omega*t)
    
    R=np.array([[c,-s,0],[s,c,0],[0,0,1]])
    Rp=omega*np.array([[-s,-c,0],[c,-s,0],[0,0,0]])
    Rpp=(omega**2)*np.array([[-c,s,0],[-s,-c,0],[0,0,0]])
    zp=np.zeros(6)
    
    zp[0] = z[3]
    zp[1] = z[4]
    zp[2] = z[5]
    
    r_3 = (np.sqrt(z[0]**2 + z[1]**2 + z[2]**2))**3
    r_7 = (np.sqrt(z[0]**2 + z[1]**2 + z[2]**2))**7
    r_9 = (np.sqrt(z[0]**2 + z[1]**2 + z[2]**2))**9
    z1 = z[0:3]
    z2 = z[3:6]
    F2 = J2*np.array([z[0]/r_9*(6*z[2]**2-3/2*(z[0]**2+z[1]**2)),z[1]/r_9*(6*z[2]**2-3/2*(z[0]**2+z[1]**2)),z[2]/r_9*(3*z[2]**2-9/2*(z[0]**2+z[1]**2))])
    F3 = J3*np.array([z[0]*z[2]/r_9*(10*z[2]**2-15/2*(z[0]**2+z[1]**2)),z[1]*z[2]/r_9*(10*z[2]**2-15/2*(z[0]**2+z[1]**2)),1/r_9*(4*z[2]**2*(z[2]**2-3*(z[0]**2)+z[1]**2)+3/2*(z[1]**2+z[1]**2)**2)])
    
    zp[3:6] =-G*Mt*z1/(r_3)+ F2 + F3 - R.T@ (2*Rp @ z2 + Rpp @ z1)
    
    return zp

