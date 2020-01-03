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

## Constants
Lv = 2.501e6
Cp = 1.00464e3

def T_in_K2thl( T, exner, rc ):
    #--------
    thl = ( T - Lv/Cp * rc ) / exner
    #--------
    return thl

secs_per_month = 86400.0 * 30.0

ncs = [netCDF4.Dataset(prefix+'/'+case+'.cam.h0.'+month+'.nc') for month in months]
ncprev = netCDF4.Dataset(prefix+'/'+case+'.cam.h0.0000-03.nc')
lev = ncs[0].variables['lev'][:]
nz = len(ncs[0].dimensions['lev'])

for ipair in range(len(pairs)):
    print("Column "+str(ipair+1))
    path=case+"_budgets/col_"+str(ipair+1)
    mkdir_p(path)
    pair = pairs[ipair]
    lat = pair[0]
    lon = pair[1]

    def read_all(vars, var_vals):
        print('Reading vars ' + str(vars))
        for var in vars:
            var_vals[var] = np.zeros(nz)

        for i in range(len(months)):
            nc = ncs[i]
            for var in vars:
                if var == 'TFIX':
                    var_vals[var] += nc.variables[var][0,lat,lon]
                else:
                    var_vals[var] += nc.variables[var][0,:nz,lat,lon]
        for var in vars:
            var_vals[var] /= len(months)

    def read_all_tend(vars, var_vals):
        print('Reading vars ' + str(vars))
        for var in vars:
            var_vals[var] = np.zeros(nz)

        for i in range(len(months)):
            nc = ncs[i]
            if i == 0:
                var_vals_prev = { }
                for var in vars:
                    var_vals_prev[var] = ncprev.variables[var][0,:,lat,lon]
            for var in vars:
                var_val_curr = nc.variables[var][0,:,lat,lon]
                var_vals[var] += (var_val_curr - var_vals_prev[var]) / secs_per_month
                var_vals_prev[var] = var_val_curr

        for var in vars:
            var_vals[var] /= len(months)

    def plot_wrapper(x, y, label):
        global style,cnt
        plt.plot(x, y, linestyle=style, label=label)
        cnt = cnt + 1
        if cnt%7 == 0:
            style = styles.pop(0)

    # Temperature budget
    var_vals_T = {}
    vars_T = ['STEND_CLUBB', 'MPDT', 'QRS', 'QRL', 'TTGWORO', 'TFIX', 'PTTEND']
    read_all(vars_T, var_vals_T)
    vars_T = ['TAP']
    read_all_tend(vars_T, var_vals_T)
    var_vals_T['TTEND_FV'] = var_vals_T.pop('TAP')

    ### Compute vars
    print('Generating T budget')
    T_terms_minus = ['STEND_CLUBB', 'MPDT', 'QRS', 'QRL', 'TTGWORO', 'TFIX']
    T_terms_minus_fact = []
    for term in T_terms_minus:
        if term in ('STEND_CLUBB', 'MPDT'):
            T_terms_minus_fact.append(1.0/Cp)
        else:
            T_terms_minus_fact.append(1.0)
    ADIAB_FV_calc1 = np.array(var_vals_T['TTEND_FV'])
    for iterm in range(len(T_terms_minus)):
        term = T_terms_minus[iterm]
        ADIAB_FV_calc1 -= (var_vals_T[term] * T_terms_minus_fact[iterm])
    ADIAB_FV_calc2 = (var_vals_T['TTEND_FV'] - (var_vals_T['PTTEND'] + var_vals_T['TFIX']))

    resid_T = (ADIAB_FV_calc2 - ADIAB_FV_calc1)

    ### Create plot
    plt.figure(figsize=(16,16))
    styles = ['-', '--', '-.']; cnt = 0; style = styles.pop(0)

    plot_wrapper(var_vals_T['TTEND_FV'], lev, label='TTEND_FV')
    plot_wrapper(ADIAB_FV_calc2, lev, label='ADIAB_FV')
    for iterm in range(len(T_terms_minus)):
        term = T_terms_minus[iterm]
        plot_wrapper(var_vals_T[term] * T_terms_minus_fact[iterm], lev, label=term)
    plot_wrapper(-resid_T, lev, label='resid')

    plt.gca().invert_yaxis()
    plt.legend()
    plt.savefig(path+'/T.svg')
    plt.close()

    # Cloud water budget
    var_vals_qc = {}
    vars_qc = ['TACLDLIQ', 'MPDLIQ', 'RCMTEND_CLUBB', 'DMECLDLIQ', 'CLDLIQTGW', 'PTECLDLIQ']
    read_all(vars_qc, var_vals_qc)

    ### Compute vars
    print('Generating CLDLIQ budget')
    terms = list(vars_qc)
    terms.pop(terms.index('PTECLDLIQ'))
    TECLDLIQ_FV_calc1 = np.zeros(nz)
    for term in terms:
        TECLDLIQ_FV_calc1 += var_vals_qc[term]
    TECLDLIQ_FV_calc2 = (var_vals_qc['TACLDLIQ'] + var_vals_qc['PTECLDLIQ'])
    resid_CLDLIQ = (TECLDLIQ_FV_calc2 - TECLDLIQ_FV_calc1)

    ### Create plot
    plt.figure(figsize=(16,16))
    styles = ['-', '--', '-.']; cnt = 0; style = styles.pop(0)
    plot_wrapper(TECLDLIQ_FV_calc2, lev, label='TECLDLIQ')
    for term in terms:
        plot_wrapper(var_vals_qc[term], lev, label=term)
    plot_wrapper(resid_CLDLIQ, lev, label='resid')

    plt.gca().invert_yaxis()
    plt.legend()
    plt.savefig(path+'/CLDLIQ.svg')
    plt.close()

    # Vapor budget
    var_vals_q = {}
    vars_q = ['TAQ', 'MPDQ', 'RVMTEND_CLUBB' ,'DMEQ', 'QTGW', 'CT_H2O', 'PTEQ']
    read_all(vars_q, var_vals_q)

    ### Compute vars
    print('Generating Q budget')
    terms = list(vars_q)
    terms.pop(terms.index('PTEQ'))
    TEQ_FV_calc1 = np.zeros(nz)
    for term in terms:
        TEQ_FV_calc1 += var_vals_q[term]
    TEQ_FV_calc2 = (var_vals_q['TAQ'] + var_vals_q['PTEQ'])
    resid_Q = (TEQ_FV_calc2 - TEQ_FV_calc1)

    ### Create plot
    plt.figure(figsize=(16,16))
    styles = ['-', '--', '-.']; cnt = 0; style = styles.pop(0)
    plot_wrapper(TEQ_FV_calc2, lev, label='TEQ')
    for term in terms:
        plot_wrapper(var_vals_q[term], lev, label=term)
    plot_wrapper(resid_Q, lev, label='resid')

    plt.gca().invert_yaxis()
    plt.legend()
    plt.savefig(path+'/Q.svg')
    plt.close()

    # Total water budget
    print('Generating QT budget')
    terms_qt = ['TAQT', 'MPDQT', 'RTMTEND_CLUBB', 'DMEQT', 'QTTGW', 'CT_H2O']
    var_vals_qt = { }
    TEQT_FV = TECLDLIQ_FV_calc2 + TEQ_FV_calc2
    var_vals_qt['TAQT'] = var_vals_qc['TACLDLIQ'] + var_vals_q['TAQ']
    var_vals_qt['MPDQT'] = var_vals_qc['MPDLIQ'] + var_vals_q['MPDQ']
    var_vals_qt['RTMTEND_CLUBB'] = var_vals_qc['RCMTEND_CLUBB'] + var_vals_q['RVMTEND_CLUBB']
    var_vals_qt['DMEQT'] = var_vals_qc['DMECLDLIQ'] + var_vals_q['DMEQ']
    var_vals_qt['QTTGW'] = var_vals_qc['CLDLIQTGW'] + var_vals_q['QTGW']
    var_vals_qt['CT_H2O'] = var_vals_q['CT_H2O']

    resid_QT = (TECLDLIQ_FV_calc2 - TECLDLIQ_FV_calc1) + (TEQ_FV_calc2 - TEQ_FV_calc1)

    ### Create plot
    plt.figure(figsize=(16,16))
    styles = ['-', '--', '-.']; cnt = 0; style = styles.pop(0)
    plot_wrapper(TEQT_FV, lev, label='TEQT')
    for term in terms_qt:
        plot_wrapper(var_vals_qt[term], lev, label=term)
    plot_wrapper(resid_QT, lev, label='resid')

    plt.gca().invert_yaxis()
    plt.legend()
    plt.savefig(path+'/QT.svg')
    plt.close()

    # Theta-l budget
    exner_dict = { }
    read_all( ['exner'], exner_dict )
    exner = exner_dict['exner']

    ### Compute vars
    print('Generating THL budget')
    terms_thl = ['THLTEND_FV','ADIAB4_FV','THLTEND_CLUBB','MPDTHL','QRS','QRL','THLTGW','TFIX','TATHL','DMETHL','resid']
    var_vals_thl = { }
    var_vals_thl['THLTEND_FV'] = T_in_K2thl( var_vals_T['TTEND_FV'], exner, TECLDLIQ_FV_calc2 )
    var_vals_thl['ADIAB4_FV'] = T_in_K2thl( ADIAB_FV_calc2, exner, 0.0 )
    var_vals_thl['THLTEND_CLUBB'] = T_in_K2thl( var_vals_T['STEND_CLUBB']/Cp, exner, var_vals_qc['RCMTEND_CLUBB'] )
    var_vals_thl['MPDTHL'] = T_in_K2thl( var_vals_T['MPDT']/Cp, exner, var_vals_qc['MPDLIQ'] )
    var_vals_thl['QRS']    = T_in_K2thl( var_vals_T['QRS'], exner, 0.0 )
    var_vals_thl['QRL']    = T_in_K2thl( var_vals_T['QRL'], exner, 0.0 )
    var_vals_thl['THLTGW'] = T_in_K2thl( var_vals_T['TTGWORO'], exner, var_vals_qc['CLDLIQTGW'] )
    var_vals_thl['TFIX']   = T_in_K2thl( var_vals_T['TFIX'], exner, 0.0 )
    var_vals_thl['TATHL']  = T_in_K2thl( 0.0, exner, var_vals_qc['TACLDLIQ'] )
    var_vals_thl['DMETHL'] = T_in_K2thl( 0.0, exner, var_vals_qc['DMECLDLIQ'] )
    var_vals_thl['resid']  = T_in_K2thl( -resid_T, exner, resid_CLDLIQ )

    ### Create plot
    plt.figure(figsize=(16,16))
    upper_altitude_hpa = 200.0
    upper_index = np.argmax(lev > upper_altitude_hpa)
    styles = ['-', '--', '-.']; cnt = 0; style = styles.pop(0)
    for var in terms_thl:
        plot_wrapper(var_vals_thl[var][upper_index:], lev[upper_index:], label=var)

    plt.gca().invert_yaxis()
    plt.legend()
    plt.savefig(path+'/THL.svg')
    plt.close()
