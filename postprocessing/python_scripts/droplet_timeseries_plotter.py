#!/usr/bin/env python

import netCDF4
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Pick a column
lat = 63
lon = 92

months = ["0000-04", "0000-05", "0000-06", "0000-07", "0000-08", "0000-09", "0000-10", "0000-11", "0000-12", "0001-01", "0001-02"]

terms = np.empty((len(months),2))

varname = "AREL"
units = None

def getvar(nc):
    var = nc.variables[varname]
    if units == None:
        units = var.units
    if varname == "CDNUMC":
        return var[0,lat,lon]
    else:
        return np.sum(var[0,:,lat,lon])

for imon in range(len(months)):
    month = months[imon]
    nc = netCDF4.Dataset('/glade/scratch/raut/repository_r77834/run/repository_r77834.cam.h0.'+month+'.nc')
    terms[imon,0] = getvar(nc)
    nc.close()

    nc = netCDF4.Dataset('/glade/scratch/raut/cam5_base/run/cam5_base.cam.h0.'+month+'.nc')
    terms[imon,1] = getvar(nc)
    nc.close()

plt.plot(range(len(months)), terms[:,0], label="CAM-CLUBB-SILHS")
plt.plot(range(len(months)), terms[:,1], label="CAM-BASE")

plt.legend()
plt.xlabel("Months since 0000-04")
plt.ylabel("Value ["+units+"]")
plt.title(varname)
plt.show()
