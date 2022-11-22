#!/usr/bin/python3
 
import subprocess
import argparse
import netCDF4
import sys
import numpy

parser = argparse.ArgumentParser(description='Run a test')

parser.add_argument("files", nargs=2,
                    help="need two files to diff")
                    
args = parser.parse_args()


# Create diff command, using ncdiff and -O to overwrite file if it exists
ncdiffCommand = "ncdiff -O " + args.files[0] + " " + args.files[1] + " diff.nc"

# Run diff command, printing output (there should be no output unless there's an error)
process = subprocess.Popen(ncdiffCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

# Create dataset from nc file
dset = netCDF4.Dataset("diff.nc")

all_zero = True
for var in dset.variables:
    
    # Let's only check 3D variables, there's some 2D ones that might be worth
    # checking, but I've never seen those differ when the 3D ones dont, and 
    # there's some ones that look like they're 1D but are for some reason
    # 2D (date_written, mdimnames, time_bnds, time_written, zlong_bnds).
    # We could also change this to instead specifically ignore only those,
    # but that feels less robust and checking only 3D variables seems like
    # a fine check.
    if dset[var].ndim == 3:
        if not numpy.all(dset[var][:] == 0):
            print(var + " is NON-ZERO")
            all_zero = False
            
if all_zero:
    print("PASS: All 3D output vars are zero.")
    sys.exit(0)
else:
    print("FAIL: There were some non-zero 3D fields.")
    sys.exit(1)
