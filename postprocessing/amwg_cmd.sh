#!/bin/bash

# This script simply sets test_casename in diag140804.csh equal
# to the command line argument, then runs that script.
#
# Usage: ./amwg_cmd.sh CASE
#
# where CASE is the name of a directory containing a CAM run, stored
# at /glade/scratch/$USER/CASE.

copyfile="`mktemp`"
cp --preserve diag140804.csh "$copyfile"
sed "s/^.*set.*test_casename.*=.*\$/set test_casename  = $1/" -i "$copyfile"
"$copyfile"
rm "$copyfile"

