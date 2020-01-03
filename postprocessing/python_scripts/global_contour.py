#!/usr/bin/env python
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

varname = 'cloud_frac'

months = ["0000-04", "0000-05", "0000-06", "0000-07", "0000-08", "0000-09", "0000-10", "0000-11", "0000-12", "0001-01", "0001-02"]

# Information
nc = netCDF4.Dataset('/glade/scratch/raut/tuning_iter1/run/tuning_iter1.cam.h0.0000-01.nc')
nlat = len(nc.dimensions['lat'])
nlon = len(nc.dimensions['lon'])
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
nc.close()

term = np.empty( (nlat,nlon) )

for month in months:
    print month
    nc = netCDF4.Dataset('/glade/scratch/raut/tuning_iter1/run/tuning_iter1.cam.h0.'+month+'.nc')
    var = nc.variables[varname]
    term = np.maximum(term, np.amax(var[0,:,:,:], axis=0))
    nc.close()

# Averaging
#term /= len(months)

# Plot
m = Basemap(llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution='c')
lonmesh, latmesh = np.meshgrid(lon,lat)
x, y = m(lonmesh,latmesh)
m.drawcoastlines()
m.drawmapboundary()
m.contourf(x,y,term)
plt.colorbar()
plt.title(varname)
plt.show()
