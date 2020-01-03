#!/bin/bash

# Variables
CASE="trunk_tag_cam5_4_175"
CASEROOT="/glade/scratch/$USER/$CASE"
MACH="cheyenne"
COMPSET="F2000_DEV"
#RES="f09_f09"
RES="f19_f19"
QUEUE="regular"
#WALL_TIME="08:00:00" # H:MM:SS or HH:MM:SS
WALL_TIME="03:00:00" # H:MM:SS or HH:MM:SS

# Configuration parameters
NUMSC=4
MGVER=2 # Currently "1" and "2" are allowed

# Obtain the CAM source directory. This solution is from Stack Overflow.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
CAM_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

# Set up case
echo "----- Case Setup -----"
cd "$CAM_DIR/cime/scripts"
./create_newcase --case "$CASEROOT" --mach "$MACH" --compset "$COMPSET" --res "$RES" --run-unsupported || {
    echo "Error creating new case" > /dev/stderr
    exit 1
}
cd "$CASEROOT"

#-----------------------------------------
# Change the number of tasks used
#-----------------------------------------
for component in ATM LND ICE OCN CPL GLC ROF WAV
do
  #./xmlchange NTASKS_$component=864 || exit 1
  ./xmlchange NTASKS_$component=720 || exit 1
done

./xmlchange JOB_QUEUE="$QUEUE"
./xmlchange JOB_WALLCLOCK_TIME="$WALL_TIME"

./case.setup || exit 1

# Compile configuration
echo "----- Compile Configuration -----"

# Compile
echo "----- Compile -----"
./case.build || { echo "Error building CAM" > /dev/stderr; exit 1; }

# Run configuration
echo "----- Run configuration -----"
./xmlchange RUN_STARTDATE="2000-01-01"
#./xmlchange STOP_OPTION=nyears
#./xmlchange STOP_N=4
./xmlchange STOP_OPTION=nmonths
./xmlchange STOP_N=14
./xmlchange DOUT_S="FALSE" # Turn on/off short-term archiving

# Run submission
echo "----- Run Submission -----"
./case.submit || { echo "Error submitting run" > /dev/stderr ; exit 1; }

# Success?!
echo "Success"'!'
exit 0
