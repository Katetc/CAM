#!/usr/bin/python3

import subprocess
import argparse
import sys
import numpy
import os
import shutil
from pathlib import Path
import fileinput
import re

parser = argparse.ArgumentParser(description='Run a test')

parser.add_argument("repoPath", nargs=1,
                    help="cam git repo path")
                    
parser.add_argument("-m", "--mach", nargs=1,
                    help="number of subcolumns")
                    
parser.add_argument("-c", "--caseroot", nargs=1,
                    help="number of subcolumns")
                    
parser.add_argument("-s", "--silhs", nargs=1,
                    help="number of subcolumns")
                    

parser.add_argument("--reuse_build", action='store_true',
                    help="number of subcolumns")
                    
args = parser.parse_args()


repo = Path(args.repoPath[0])

if not repo.exists():
    print("Path does not exist. Cloning cam into: "+args.repoPath)
    #clone repo
    
os.chdir(repo)

print("Checking out externals...")
checkoutExternalsCmd = "python3 manage_externals/checkout_externals -e Externals.cfg"
process = subprocess.Popen(checkoutExternalsCmd.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


print("Copying cime config files")
machineOptions = "cime_config_cesm_machines_files/config_machines.xml"
compilerOptions = "cime_config_cesm_machines_files/config_compilers.xml"
cimeConfigLocation = Path("cime/config/cesm/machines")
shutil.copy(machineOptions, cimeConfigLocation)
shutil.copy(compilerOptions, cimeConfigLocation)


for line in fileinput.input("run_scripts/run_cesm_uwm_coarse_res.sh", inplace=True):
    
    if args.caseroot is not None:
        line = re.sub('CASEROOT=".*"','CASEROOT="'+args.caseroot[0]+'"',line)
        
    if args.mach is not None:
        line = re.sub('MACH=".*"','MACH="'+args.mach[0]+'"',line)
        
    if args.silhs is not None:
        line = re.sub('NUMSC=.*','NUMSC='+args.silhs[0],line)
        
    if args.reuse_build:
        line = re.sub('create_newcase','create_newcase --handle-preexisting-dirs u',line)
        
        
    print(line.rstrip("\r\n"))




    
    




