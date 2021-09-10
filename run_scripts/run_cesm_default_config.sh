#!/bin/bash

# Variables
CASE="test_CLUBB_merge_2"
CASEROOT="/glade/scratch/$USER/$CASE"
MACH="cheyenne"
COMPSET="F2000climo"
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
./create_newcase --case "$CASEROOT" --mach "$MACH" --compset "$COMPSET" --res "$RES" --run-unsupported 
if [ $? -gt 0 ]
then 
    echo "Error creating new case" > /dev/stderr
    exit 1
fi

cd "$CASEROOT"

#-----------------------------------------
# Change the number of tasks used
#-----------------------------------------
for component in ATM LND ICE OCN CPL GLC ROF WAV
do
  #./xmlchange NTASKS_$component=864 || exit 1
  ./xmlchange NTASKS_$component=720 
if [ $? -gt 0 ]
then
    echo "Failed at the xmlchange step" > /dev/stderr
    exit 1
fi
done

./xmlchange JOB_QUEUE="$QUEUE"
./xmlchange JOB_WALLCLOCK_TIME="$WALL_TIME"

./case.setup
if [ $? -gt 0 ]
then 
   echo "Failed at the case.setup step" > /dev/stderr
   exit 1
fi

# Compile configuration
echo "----- Compile Configuration -----"

# Compile
echo "----- Compile -----"
./case.build 
if [ $? -gt 0 ]
then
   echo "Error building CAM" > /dev/stderr
   exit 1
fi

# Run configuration
echo "----- Run configuration -----"
./xmlchange RUN_STARTDATE="2000-01-01"
#./xmlchange STOP_OPTION=nyears
#./xmlchange STOP_N=4
./xmlchange STOP_OPTION=nmonths
./xmlchange STOP_N=14
./xmlchange DOUT_S="FALSE" # Turn on/off short-term archiving

cat >> user_nl_cam << EOF
nhtfrq = 0,-24,-6,-3
mfilt = 1,5000,5000
history_amwg = .true.
fincl1 = 'U:A','PS:A','T:A','V:A','OMEGA:A','Z3:A','PRECT:A',
'CLDLIQ:A', 'CLDICE:A', 'LWCF:A', 'SWCF:A', 'FLUT:A',
'TMQ:A', 'PRECC:A', 'PRECL:A', 'CME:A', 'PRODPREC:A',
'EVAPPREC:A','EVAPSNOW:A','ICWMRST:A','ICIMRST:A','PRAO:A',
'PRCO:A','QCSEVAP:A','QISEVAP:A','QVRES:A','CMEIOUT:A','VTRMI:A',
'VTRMC:A','QCSEDTEN:A','QISEDTEN:A','MNUCCCO:A','MNUCCTO:A'
fincl2 = 'CLDTOT', 'CLDST','CDNUMC','CLDLIQ','CLDICE','FLUT',
'LWCF','SWCF','PRECT'
EOF

# Run submission
echo "----- Run Submission -----"
./case.submit 
if [ $? -gt 0 ]
then
    echo "Error submitting run" > /dev/stderr
    exit 1
fi

# Success?!
echo "Success"'!'
exit 0
