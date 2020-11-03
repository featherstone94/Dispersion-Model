# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:11:13 2020

@author: tjfea
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

H=15   #Effective height of stack in m
Q=0.01    #contaminant emission rate in Kg/s
u=5   #wind velocity in m/s



posx=200
posy=340
#define turbulance perameters 
Ry=[0.443,0.324,0.216,0.141,0.105,0.071]
ry=[0.894,0.894,0.894,0.894,0.894,0.894]
Iy=[-1.0104,-1.634,-2.054,-2.555,-2.754,-3.243]
Jy=[0.9878,1.350,1.02331,1.0423,1.0106,1.0148]
Ky=[-0.0076,-0.0096,-0.0076,-0.0087,-0.0064,-0.0070]
Iz=[4.678,-1.999,-2.341,-3.186,-3.783,-4.490]
Jz=[-1.7172,0.8752,0.9477,1.1737,1.2010,1.4024]
Kz=[0.2770,0.0136,-0.002,-0.0316,-0.0450,-0.0540]

Pasquill_classification = 1 #Plasquill classification: A=0, B=1, C=2... ect

area=500
resolution=1

x=np.arange(0,area,resolution)
y=np.arange(-area/2,area/2,resolution)
z=np.arange(0,100)

X,Y=np.meshgrid(x,y)


def particulate_concentration(x,y,z,Q,u,H,sigma_y,sigma_z):
    f=np.exp(-y**2/(2*sigma_y**2))
    g1=np.exp(-(z-H)**2/(2*sigma_z**2))
    g2=np.exp(-(z+H)**2/(2*sigma_z**2))
    C=Q/(2*np.pi*u*sigma_y*sigma_z)*f*(g1+g2)
    return C

def standard_dev(x,I,J,K):
    sigma=np.exp(I+J*np.log(x)+K*(np.log(x))**2)
    return sigma

vertical_sd=np.zeros(len(x))
horizontal_sd=np.zeros(len(x))

for n in range(0,len(x)):
    vertical_sd[n]=standard_dev(x[n],Iz[Pasquill_classification],Jz[Pasquill_classification],Kz[Pasquill_classification])
    horizontal_sd[n]=standard_dev(x[n],Iy[Pasquill_classification],Jy[Pasquill_classification],Ky[Pasquill_classification])

results=particulate_concentration(X,Y,0,Q,u,H,horizontal_sd,vertical_sd)
plt.figure(0)
plt.pcolormesh(x,y,results,cmap="YlOrBr")
plt.colorbar(label="ground level contamination (kg/$m^3$)")
plt.scatter(0,0,c="r")
plt.xlabel("down wind position (m)")
plt.ylabel("cross wind position (m)")

plt.show()
