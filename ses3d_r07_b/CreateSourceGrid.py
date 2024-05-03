#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

## Read in boundaries of computational domain
fid = open('INPUT/setup').readlines()
for line in fid:
    if line.find('theta_min') >= 0: colat_min = float(line.split()[0]); lat_max = 90 - colat_min
    if line.find('theta_max') >= 0: colat_max = float(line.split()[0]); lat_min = 90 - colat_max
    if line.find('phi_min') >= 0: lon_min = float(line.split()[0])
    if line.find('phi_max') >= 0: lon_max = float(line.split()[0])

## Read source location from file
fid = open('INPUT/source.txt').readlines()
for line in fid:
    if line.find('latitude') >= 0: src_lat = float(line.split()[0])
    if line.find('longitude') >= 0: src_lon = float(line.split()[0])
    if line.find('depth') >= 0: src_dep = float(line.split()[0])
    
## Read receliver locations from file
rec_name = np.loadtxt('INPUT/stations.txt', usecols=[0], dtype=str)
rec_lat = np.loadtxt('INPUT/stations.txt', usecols=[1])
rec_lon = np.loadtxt('INPUT/stations.txt', usecols=[2])

## Create mesh of source
src_glat = src_lat + np.arange(-.2, .21, .1)
src_glon = src_lon + np.arange(-.2, .21, .1)
src_gdep = src_dep + np.arange(0, 15001, 1000)
src_mlat, src_mlon, src_mdep = np.meshgrid(src_glat, src_glon, src_gdep)
count = 0
fid = open('INPUT/source_grid.txt', 'w')
for (lat, lon, dep) in zip(src_mlat.flatten(), src_mlon.flatten(), src_mdep.flatten()):
    fid.write('GR%05d %10.5f %10.5f %8.1f\n' % (count, lat, lon, dep))
    count += 1

## Plot map and station
m = Basemap(projection='merc',llcrnrlat=lat_min,urcrnrlat=lat_max,\
            llcrnrlon=lon_min,urcrnrlon=lon_max,lat_ts=.5*(lat_max+lat_min),resolution='f')
m.etopo(scale=1)
m.drawcoastlines()
m.drawparallels(np.arange(-90.,91.,2.), labels=(1, 1, 0, 0))
m.drawmeridians(np.arange(-180.,181.,2.), labels=(0, 0, 1, 1))

xx, yy = m(src_mlon[:, :, 0], src_mlat[:, :, 0])
#plt.plot(xx, yy, 'or', ms=1)
plt.plot(xx, yy, 'vb', ms=1)
xx, yy = m(rec_lon, rec_lat)
#plt.plot(xx, yy, 'vb', ms=5)
plt.plot(xx, yy, 'or', ms=5)
for x, y, n in zip(xx, yy, rec_name):
    plt.text(x, y +10e3, n)

plt.tight_layout()
plt.savefig('EventMap1.pdf')
plt.show()
