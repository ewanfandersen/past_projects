# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 12:18:29 2022

@author: Ewan Andersen
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from os.path import isfile, join

directory = 'C:/ewana/Gaia/'

data = np.loadtxt(directory + 'ASTR301_NGC_2169_cleaned.txt', skiprows=2)

iso_dir = directory + 'Isochrones/'
iso = np.sort([iso_dir+f for f in os.listdir(iso_dir) if isfile(join(iso_dir,f))])

RA = data[:,0] #right ascension
DE = data[:,1] #declination
plx = data[:,2] #parallax
e_plx = data[:,3] #error in parallax
pmRA = data[:,4] #proper motion right ascension
e_pmRA = data[:,5] #error in proper motion right ascension
pmDE = data[:,6] #proper motion declination
e_pmDE = data[:,7] #error in proper motion declination
Gmag = data[:,8] #G magnitude
e_Gmag = data[:,9] #error in G magnitude
BPmag = data[:,10] #BP magnitude
e_BPmag = data[:,11] #error in BP magnitude
RPmag = data[:,12] #RP magnitude
e_RPmag = data[:,13] #error in RP magnitude

plt.figure()
plt.hist(plx, bins=200, range=[0.5,6])#histogram of the raw data from VizieR
plt.show()

plt.figure()
plt.plot(pmRA, pmDE, '.', color='green', markersize=0.5)#raw data plot of proper motion right ascension and declination

filter1 = np.where((plx>1.0)&(plx<1.3)&(e_plx<0.2))#isolates the parallaxes of the star cluster that we are interested in, or any row with a parallax between 1.0 and 1.3
filter2 = np.where((pmRA>1.9)&(pmRA<2.7)&(pmDE>-3.3)&(pmDE<-2.37))#isolates the proper motion RA and Dec of the cluster
filter3 = np.where((plx>1.0)&(plx<1.3)&(pmRA>1.9)&(pmRA<2.7)&(pmDE>-3.25)&(pmDE<-2.38))#combination of the above two filters
filter4 = np.where((plx>1.0)&(plx<1.3)&(pmRA>1.9)&(pmRA<2.7)&(pmDE>-3.25)&(pmDE<-2.38)&(e_plx<0.2))#combination of all filters below with an extra filter for excluding the error in parallax greater than 20%.

plt.figure()
plt.hist(plx, bins=200, range=[0.5,6],zorder=1, label='all sources')#base histogram of raw parallaxes from VizieR
plt.hist(plx[filter1], bins=200, range=[0.5,6], zorder=2, label='filtered parallax')#parallax filtered between plx 1 and 1.3
plt.hist(plx[filter4], bins=200, range=[0.5,6], zorder=3, label='filtered parallax and proper motion')#parallax filtered by cluster parallax as well as proper motion of RA and Dec
plt.xlabel('parallax (mas)')
plt.ylabel('number of sources')
plt.title('Parallax counts')
plt.legend()
plt.show()     


plt.figure()
plt.plot(pmRA, pmDE, '.', color='green', markersize=0.5, zorder=1,label='all sources') #plots the proper motion of the RA and Dec
plt.plot(pmRA[filter1], pmDE[filter1], '.', color='orange', markersize=0.5, zorder=2,label='filtered parallax') #plots the proper motion of the RA and Dec with the parallax filter
plt.plot(pmRA[filter2], pmDE[filter2], '.', color='red', markersize=1, zorder=3, label='filtered parallax and proper motion of RA and Dec')#plots the same as above but with filter 3
plt.xlabel('proper motion in right ascension (mas/yr)')
plt.ylabel('proper motion in declination (mas/yr)')
plt.title('Area of Cluster')
plt.legend(markerscale=10)
plt.show()

color=BPmag-RPmag #color is the difference between the BP and the RP magnitudes
e_color=np.sqrt((e_BPmag)**2+(e_RPmag)**2) #adding the errors of these two magnitudes in quadriture to find the error in color

plt.figure()
plt.plot(color, Gmag, '.',  color='blue', markersize=0.5, alpha=0.25, zorder=1, label='all sources')
plt.errorbar(color[filter3], Gmag[filter3], xerr=e_color[filter3],fmt='.', alpha=0.5, color='purple', markersize=1, zorder=2, label='NGC 2168')
plt.errorbar(color[filter4], Gmag[filter4], xerr=e_color[filter4],fmt='.', color='red', markersize=1, zorder=3, label='NGC 2168 < 20%  err in plx')
plt.xlabel('BP-RP')
plt.ylabel('G magnitude')
plt.gca().invert_yaxis()
plt.title('Color-Magnitude Diagram w/ filters')
plt.legend(markerscale=10)
plt.show()

plt.figure()
plt.errorbar(color[filter4], Gmag[filter4], xerr=e_color[filter4],fmt='.', color='purple', markersize=1, zorder=3, label='NGC 2168')
plt.xlabel('BP-RP')
plt.ylabel('G magnitude')
plt.gca().invert_yaxis()
('Color-Magnitude Diagram w/ isolated cluster')
plt.legend(markerscale=10)
plt.show()

RA_range = np.max(RA[filter3]) - np.min(RA[filter3])
Dec_range = np.max(DE[filter3]) - np.min(DE[filter3])
plx_range = (1/(np.min(plx[filter3])*10**-3))-(1/(np.max(plx[filter3])*10**-3))

Dec_range_arcsecs = Dec_range * 3600 * (816/206265)
RA_range_arcsecs = RA_range * 3600 * (816/206265) * np.cos(np.deg2rad(Dec_range_arcsecs))

isochrones = []
for ii in range(len(iso)):
    data = np.loadtxt(iso[ii])
    isochrones.append(data)

def app(num):
   return num + (5*np.log10(816/10))#converts magnitudes from absolute to relative.

mags = [] #g, gb, gr at 100, 150, 175, 200, 250 Myrs. example: [[g100,gb100,gr100], [g150,gb150,gr150]]
for jj in range(len(isochrones)):
    filter_mags = [] 
    for kk in range(28,31):
        filter_col = isochrones[jj][:, kk]
        filtered_values = app(filter_col)
        filter_mags.append(filtered_values)
    mags.append(filter_mags)

colors = []
for ii in range(len(mags)):
    color_ii = mags[ii][1] - mags[ii][2]
    colors.append(color_ii)

#BP-RP magnitudes to get the x-axis portion of the CMD (the 'G-magnitude' parts cancel leaving only BP and RP)
plt.figure()
plt.plot(color[filter4], Gmag[filter4], '.', color='blue', markersize=2, zorder=1, label='NGC 2168')
for i in range(len(colors)):
    plt.plot(colors[i] + 0.4, mags[i][0] + 0.8, color=['pink', 'green', 'red', 'orange', 'purple'][i], markersize=2, zorder=i+2, label=f'{[100, 150, 175, 200, 250][i]} Myr isochrone')
plt.xlabel('BP-RP')
plt.ylabel('G magnitude')
plt.gca().invert_yaxis()                                                                    
plt.legend(markerscale=10)
plt.show()











