# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:11:13 2020

@author: tjfea
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

H=[15,15,14,15]   #Effective height of stacks in m
Q=[0.04,0.04,0.03,0.03]    #contaminant emission rates in Kg/s
u=5   #wind velocity in m/s (assumed constant for all stacks)

posx=[100,150,250,250]
posy=[30,-100,0,200]


#define turbulance perameters from Pasquill's atmospheric stability class. 
#Classifications A through F from indecies 0 through 5
Ry=[0.443,0.324,0.216,0.141,0.105,0.071]
ry=[0.894,0.894,0.894,0.894,0.894,0.894]
Iy=[-1.0104,-1.634,-2.054,-2.555,-2.754,-3.243]
Jy=[0.9878,1.350,1.02331,1.0423,1.0106,1.0148]
Ky=[-0.0076,-0.0096,-0.0076,-0.0087,-0.0064,-0.0070]
Iz=[4.678,-1.999,-2.341,-3.186,-3.783,-4.490]
Jz=[-1.7172,0.8752,0.9477,1.1737,1.2010,1.4024]
Kz=[0.2770,0.0136,-0.002,-0.0316,-0.0450,-0.0540]

Pasquill_classification = 2 #Plasquill classification: A=0, B=1, C=2... ect

area=2000    #area over which to calculate in m^2 
resolution=1    #resolution of calculations in m

#create arrays for coordinates for calculations
x=np.arange(0,area,resolution)
y=np.arange(-area/2,area/2,resolution)
z=np.arange(0,100)

#deine gaussian plume equation
def particulate_concentration(x,y,z,Q,u,H,sigma_y,sigma_z):
    f=np.exp(-y**2/(2*sigma_y**2))
    g1=np.exp(-(z-H)**2/(2*sigma_z**2))
    g2=np.exp(-(z+H)**2/(2*sigma_z**2))
    C=Q/(2*np.pi*u*sigma_y*sigma_z)*f*(g1+g2)
    return C

#define standard dev for 
def standard_dev(x,I,J,K):
    sigma=np.exp(I+J*np.log(x)+K*(np.log(x))**2)
    return sigma

#create meshgrid
X,Y=np.meshgrid(x,y)

#prime arrays
concentration=[]

#loop through the different sources
for i in range(0,len(H)):

    #prime sd arrays
    vertical_sd=np.zeros(len(x))
    horizontal_sd=np.zeros(len(x))
    
    #calculate vertical and horizontal standard deviations
    for n in range(0,len(x)):
        vertical_sd[n]=standard_dev(x[n],Iz[Pasquill_classification],Jz[Pasquill_classification],Kz[Pasquill_classification])
        horizontal_sd[n]=standard_dev(x[n],Iy[Pasquill_classification],Jy[Pasquill_classification],Ky[Pasquill_classification])
    
    #calculated concentration at ground level 
    results=particulate_concentration(X,Y,0,Q[i],u,H[i],horizontal_sd,vertical_sd)
    #convert to pandas dataframe for shift function
    results2=pd.DataFrame(results)
    #shift x and y position to correct for source position
    results3=results2.shift(periods=posy[i],fill_value=0)
    results4=results3.shift(periods=posx[i],axis="columns",fill_value=0)
    concentration.append(results4)


final=pd.DataFrame(np.array([np.zeros(area)]*area).T)
for j in range (0,len(H)):
    final=final+concentration[j]

plt.figure(0)
#plot hat map
plt.pcolormesh(x,y,final,cmap="YlOrBr")
#create color bar key
plt.colorbar(label="ground level contcentration (kg/$m^3$)")
#plot positions of sources as scatter
plt.scatter(posx,posy,c="Black",marker="x")
plt.ylim(np.amax(posy)-(area/2),np.amin(posy)+(area/2))
plt.xlabel("down wind position (m)")
plt.ylabel("cross wind position (m)")

plt.show()
