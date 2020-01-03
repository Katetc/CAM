#!/usr/bin/env python

import netCDF4
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Pick a column (there's only one)
lat = 0
lon = 0

nc = netCDF4.Dataset('/home/raut/work/cam_ticket_58/subcol_SILHS_UWM/run_plots/runs/camclubb707_L30_T1200/arm97/run/camclubb707_L30_T1200.cam.h0.0097-06-18-84585.nc')
nts = len(nc.dimensions['time'])

terms = np.empty((nts,2))

varname = "AREL"
units = None

def getvar(nc):
    var = nc.variables[varname]
    global units
    if units == None:
        units = var.units
    if varname == "CDNUMC":
        return var[:,lat,lon]
    else:
        return np.sum(var[:,:,lat,lon], axis=1)

terms[:,0] = getvar(nc)
nc.close()

nc = netCDF4.Dataset('/home/raut/work/cam_ticket_58/subcol_SILHS_UWM/run_plots/runs/base707_L30_T1200/arm97/run/base707_L30_T1200.cam.h0.0097-06-18-84585.nc')
terms[:,1] = getvar(nc)
nc.close()

print("Averages:")
print("CAM-CLUBB-SILHS: " + str(np.average(terms[:,0])))
print("CAM-BASE: " + str(np.average(terms[:,1])))

plt.plot(range(nts), terms[:,0], label="CAM-CLUBB-SILHS")
plt.plot(range(nts), terms[:,1], label="CAM-BASE")

plt.legend()
plt.xlabel("Timestep")
plt.ylabel("Value ["+units+"]")
plt.title(varname)
plt.show()
