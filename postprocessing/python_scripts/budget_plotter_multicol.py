#!/usr/bin/env python

import netCDF4
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import errno

case = sys.argv[1]

prefix = "/glade/scratch/"+os.environ['USER']+"/"+case+"/run"

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

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

for ipair in range(len(pairs)):
    print("Column "+str(ipair+1))
    path=case+"_budgets/col_"+str(ipair+1)
    mkdir_p(path)
    pair = pairs[ipair]
    lat = pair[0]
    lon = pair[1]


    def make_plot(budget_terms, filename):
        for month in months:
            nc = netCDF4.Dataset(prefix+"/"+case+".cam.h0."+month+".nc")
            # First month only:
            if month == months[0]:
                ilev = nc.variables['ilev'][:]
                lev  = nc.variables['lev'] [:]
                nzm = len(nc.dimensions['ilev'])
                nzt = len(nc.dimensions['lev'])

                if nc.variables[budget_terms[0]].shape[1] == nzm:
                    levvar = ilev
                    nz = nzm
                else:
                    levvar = lev
                    nz = nzt
                term_sums = np.zeros((nz,len(budget_terms)))
            for i in range(len(budget_terms)):
                var = nc.variables[budget_terms[i]]
                term_sums[:,i] += var[0,:,lat,lon]
            nc.close()

        terms = term_sums[:,:] / len(months)

        residual = np.zeros(nz)
        residual += terms[:,0]
        for i in range(1,len(budget_terms)):
            residual -= terms[:,i]

        # Plotting
        plt.figure(figsize=(16,16))
        # matplotlib cycles through 7 colors, so when it cycles we want to change
        # the line style.
        linestyles = ['-', '--', '-.']
        linestyle = linestyles.pop(0)
        plot_num = 0;
        for i in range(len(budget_terms)):
            plt.plot(terms[:,i], levvar, linestyle, label=budget_terms[i])
            plot_num += 1
            if plot_num % 7 == 0:
                linestyle = linestyles.pop(0)

        plt.plot(residual, levvar, linestyle, label='residual')

        plt.gca().invert_yaxis()

        plt.legend()
        plt.xlabel('Budget')
        plt.ylabel('Altitude [hPa]')
        plt.title('col_'+str(ipair+1))
        plt.savefig(path+"/"+filename)
        plt.close()

    print("rtp2")
    term = 'rtp2'
    budget_ends = ['_bt', '_ma', '_ta', '_tp', '_dp1', '_dp2', '_cl', '_pd', '_sf', '_forcing']
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "rtp2.svg")

    print("thlp2")
    term = 'thlp2'
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "thlp2.svg")

    print("wprtp")
    term = 'wprtp'
    budget_ends = ['_bt', '_ma', '_ta', '_tp', '_ac', '_bp', '_pr1', '_pr2', '_pr3', '_dp1', '_mfl', '_cl', '_sicl', '_pd', '_forcing']
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "wprtp.svg")

    print("wpthlp")
    term = 'wpthlp'
    budget_ends = ['_bt', '_ma', '_ta', '_tp', '_ac', '_bp', '_pr1', '_pr2', '_pr3', '_dp1', '_mfl', '_cl', '_sicl', '_forcing']
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "wpthlp.svg")

    print("rtpthlp")
    term = 'rtpthlp'
    budget_ends = ['_bt', '_ma', '_ta', '_tp1', '_tp2', '_dp1', '_dp2', '_cl', '_sf', '_forcing']
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "rtpthlp.svg")

    print("wp2")
    term = 'wp2'
    budget_ends = ['_bt', '_ma', '_ta', '_ac', '_bp', '_pr1', '_pr2', '_pr3', '_dp1', '_dp2', '_cl', '_pd', '_sf']
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "wp2.svg")

    print("wp3")
    term = 'wp3'
    budget_ends = ['_bt', '_ma', '_ta', '_tp', '_ac', '_bp1', '_bp2', '_pr1', '_pr2', '_dp1', '_cl']
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "wp3.svg")

    print("up2")
    term = 'up2'
    budget_ends = ['_bt', '_ma', '_ta', '_tp', '_dp1', '_dp2', '_pr1', '_pr2', '_cl', '_pd', '_sf']
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "up2.svg")

    print("vp2")
    term = 'vp2'
    budget_ends = ['_bt', '_ma', '_ta', '_tp', '_dp1', '_dp2', '_pr1', '_pr2', '_cl', '_pd', '_sf']
    budget_terms = [term + budget_ends_i for budget_ends_i in budget_ends]
    make_plot(budget_terms, "vp2.svg")
