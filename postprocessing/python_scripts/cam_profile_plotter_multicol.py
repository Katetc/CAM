#!/usr/bin/env python

import netCDF4
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os

case = sys.argv[1]

prefix='/glade/scratch/'+os.environ['USER']+'/'+case+'/run'

# An array of dictionaries of plot options.
plots = []

# Add variables
varnames = ['cloud_frac', 'SILHS_CLUBB_PRECIP_FRAC', 'SILHS_CLUBB_ICE_SS_FRAC', 'thlm', 'rtm', 'rcm', 'wp2', 'wp3', 'wpthlp', 'wprtp', 'thlp2', 'rtp2', 'rtpthlp', 'OMEGA', 'AQRAIN', 'ANRAIN', 'AWNC', 'CLDICE', 'CLDLIQ', 'CLOUD', 'AWNI', 'AQSNOW', 'ANSNOW', 'wpthvp', 'T', 'C7_Skw_fnc', 'um', 'vm', 'upwp', 'vpwp', 'C6rt_Skw_fnc', 'C6thl_Skw_fnc', 'Lscale', 'C11_Skw_fnc', 'up2', 'vp2', 'Kh_zm', 'Q', 'SL', 'RHW', 'QRS', 'QRL', 'HR', 'FDL']

upper_altitude_hpa = 150.0

for varname in varnames:
    plots.append({'vname': varname})

### Generate plots ###
months = ["0000-04", "0000-05", "0000-06", "0000-07", "0000-08", "0000-09", "0000-10", "0000-11", "0000-12", "0001-01", "0001-02"]

# Pick columns: lat-lon pairs
columns = [(20,190), (27,240), (-20,275), (-20,285), (-5,355), (-1,259), (-60,340), (60,180), (2,140), (9,229), (56,311), (76,320), (45,180), (-0,295)]

# Index pairs
pairs = []
nc = netCDF4.Dataset(prefix+'/'+case+".cam.h0."+months[0]+".nc")
nclat = nc.variables['lat'][:]
nclon = nc.variables['lon'][:]
for column in columns:
    latidx = np.abs(nclat-column[0]).argmin()
    lonidx = np.abs(nclon-column[1]).argmin()
    pairs.append( (latidx,lonidx) )
nc.close()

os.mkdir(case+"_profiles")
for ipair in range(len(pairs)):
    print("Column "+str(ipair+1))
    pair = pairs[ipair]
    lat = pair[0]
    lon = pair[1]
    path = case+'_profiles/col_'+str(ipair+1)
    os.mkdir(path)

    # Another list of dictionaries
    terms = []

    firstMonth = True
    for month in months:
        print(month)
        nc = netCDF4.Dataset(prefix+'/'+case+'.cam.h0.'+month+'.nc')
        for iplot in range(len(plots)):
            plot = plots[iplot]
            ncvar = nc.variables[plot['vname']]
            var = ncvar[0,:,lat,lon]
            if firstMonth == True:
                levvar = nc.variables[ncvar.dimensions[1]]
                terms.append({'array': var, 'units': ncvar.units, 'long_name': ncvar.long_name, 'levels': levvar[:], 'levelunits': levvar.units})
            else:
                terms[iplot]['array'] += var
        firstMonth = False

    # Averaging
    for term in terms:
        term['array'] /= len(months)

    # Plotting
    for iplot in range(len(plots)):
        plot = plots[iplot]
        term = terms[iplot]
        plt.figure()
        upper_index = np.argmax(term['levels'] > upper_altitude_hpa)
        plt.plot(term['array'][upper_index:], term['levels'][upper_index:])
        plt.xlabel(plot['vname'] + ' [' + term['units'] + ']')
        plt.ylabel('Altitude [' + term['levelunits'] + ']')
        plt.title('col_'+str(ipair+1)+': '+plot['vname'] + ': ' + term['long_name'])
        plt.gca().invert_yaxis()
        # Use scientific notation for numbers outside (10^-3, 10+3) (in absolute value??)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(-3,3))
        plt.savefig(path+'/'+plot['vname']+'.svg')
        plt.close()
