#!/bin/bash

##########################################################################
# This is a convenience script that runs the AMWG diagnostics, budget
# plotter (in python_scripts), and profile plotter (also in
# python_scripts), and stores all the output in a ZIP archive.
#
# Usage: ./all_diagnostics.bash CASE
#
# where CASE is the name of a directory containing a CAM run, stored
# at /glade/scratch/$USER/CASE.
##########################################################################

scratch=/glade/scratch/$USER

( echo "Running AMWG Diagnostics"
./amwg_cmd.sh $1 &&
./CAM_amwgplots2pdf.bash -f $scratch/amwg/$1/*.tar &&
mv $scratch/amwg/$1/*/UWM_diag/*.pdf . ) &

( echo "Running budget plotter"
python_scripts/budget_plotter_multicol.py $1 ) &

( echo "Running CAM profile plotter"
python_scripts/cam_profile_plotter_multicol.py $1 ) &

( echo "Running CAM budget plotter"
python_scripts/cam_budget_plotter.py $1 ) &

wait

zip -r $1.zip "$1"*
