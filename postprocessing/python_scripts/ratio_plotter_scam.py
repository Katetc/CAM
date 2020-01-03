#!/usr/bin/env python

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

ncfalse = netCDF4.Dataset('/home/raut/work/cam_ticket_58/subcol_SILHS_UWM/l_const_false/camclubb707_L30_T1200.cam.h0.0097-06-18-84585.nc')
nctrue = netCDF4.Dataset('/home/raut/work/cam_ticket_58/subcol_SILHS_UWM/l_const_true/camclubb707_L30_T1200.cam.h0.0097-06-18-84585.nc')

ncs = [ncfalse,nctrue]

avgs = [np.mean(nc.variables['QCRAT'], axis=0)[:,0,0] for nc in ncs]

heights = ncfalse.variables['lev'][:]

labels = ['l_const_false', 'l_const_true']

for i in range(len(avgs)):
    plt.plot(avgs[i], heights, label=labels[i])

plt.gca().invert_yaxis()
plt.legend()
plt.xlabel("Ratio")
plt.ylabel("Altitude [hPa]")
plt.show()
