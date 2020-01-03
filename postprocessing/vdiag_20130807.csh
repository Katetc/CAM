#!/bin/tcsh -f

# to use this package, set the following namelist variables as
# fincl2 = 'PRECT','U200','V200'
# fincl3 = 'OMEGA500','PRECT','U200','U850','FLUT'
# fincl4 = 'PRECT'
#
# avgflag_pertape(2) = 'A'
# avgflag_pertape(3) = 'A'
# avgflag_pertape(4) = 'A'
#
# nhtfrq(2) = -24  #h1 daily
# nhtfrq(3) = -6   #h2 6-hourly
# nhtfrq(4) = -3   #h3 3-hourly
#
# mfilt(2) = 30    # about 1 file per month for each type
# mfilt(3) = 120
# mfilt(4) = 240
# 
# these will produce
# monthly  avg h0
# daily    avg h1  'PRECT','U200','V200'
# 6-hourly avg h2  'OMEGA500','PRECT','U200','U850','FLUT'
# 3-hourly avg h3  'PRECT'


set echo verbose

# ======================================================================
# set up
# ======================================================================
#
# which packages to run
#
set do_monthly = 0    # 0=off, 1=on
set do_daily   = 0    # 0=off, 1=on
set do_hourly  = 1    # 0=off, 1=on    # Wavenumber/Power-spectra Wheeler/Kilidas, frequency of 1,3,6,24 hrs possible
set do_diurnal = 1    # 0=off, 1=on    # diurnal cycle of precipitation
set make_variab_tar   = 1    # 1 = create tar file of results, send to webdir if set


setenv CLEAN 1

# -------------------------
# Case info
#
setenv CASENAME r87242

# Declared the simulation to be non-AMIP
# because this needs to be known by 
# the diurnal cycle plotting script (precip_dc.csh)
# V. Larson 20170211
setenv amip 0
#
# Year stamp of data,  re-written below in 4-digit yr
#
# if the time range is different for different frequencies of model output
# set here
#
set firstyr_monthly = 0
set lastyr_monthly  = 1

set firstyr_daily = 0
set lastyr_daily  = 0

set firstyr_hourly = $firstyr_daily
set lastyr_hourly  = $lastyr_daily

setenv firstyr_diurnal $firstyr_monthly
setenv lastyr_diurnal  $lastyr_monthly

# 
# wkdir (ncl will generate .ps files here)
#
setenv WKDIR /glade/scratch/vlarson/vdiag/$CASENAME
echo "using WKDIR $WKDIR"


# 
# Path to input data. The different types of files can be
#
# EITHER
# 1) in separate directories DIRPATH_HOURLY, DIRPATH_DAILY, DIRPATH_MONTHLY 
#               e.g. setenv DIRPATH_HOURLY  $DATADIR/h2
#                    setenv DIRPATH_DAILY   $DATADIR/h1
#                    setenv DIRPATH_MONTHLY $DATADIR/h0
# OR
# 2) provide the FILESTR_X tags to search for
#               e.g. setenv FILESTR_HOURLY  h2;  setenv DIRPATH_HOURLY $DATADIR/atm/hist
#
# Okay to do both (1) and (2)           
#
# Used like this (hourly.csh): ls $DIRPATH_HOURLY/*${filestr}*.nc > $LIST_HOURLY_IN
#

setenv DATADIR  /glade/scratch/vlarson/$CASENAME

setenv DIRPATH_HOURLY  $DATADIR/run_h2
setenv DIRPATH_DAILY   $DATADIR/run
setenv DIRPATH_MONTHLY $DATADIR/run
setenv DIRPATH_DIURNAL $DATADIR/run   # 3-hourly output

setenv FILESTR_HOURLY   h2
setenv FILESTR_DAILY    h1
setenv FILESTR_MONTHLY  h0
setenv FILESTR_DIURNAL  h3

#
# LANDFRAC is required for daily.csh (prect_hist.ncl )
#    Some resolutions have the field provided in var_data/landmasks
#    Others need to provide the path to a file with the LANDFRAC field *in the right resolution*
#
#    To make your own, use CAM monthly output or initial conditions file (input.nc)
#       ncks -v LANDFRAC input.nc output.nc                                                                                                               
#       setenv LANDFRAC_FILE output.nc 
#       setenv LANDFRAC_NAME LANDFRAC      # name of the field in the nc file 
#
setenv LANDFRAC_FILE /glade/p/work/cchen/VDIAG_20130807/var_data/landmasks/FV_LANDFRAC_1.9x2.5.nc                                                      
setenv LANDFRAC_NAME LANDFRAC             
#setenv LANDFRAC_FILE /glade/u/home/bundy/diag/variab_amwg_yaga_mine/var_data/landmasks/USGS-gtopo30_0.9x1.25_remap_c051027_landmask.nc
#setenv LANDFRAC_NAME LANDFRAC


# 
# If WEBDIR is set, this script will copy and untar the file there. 
# Otherwise, leave unset
#
#unsetenv WEBDIR
#
setenv WEBDIR  /glade/scratch/vlarson/vdiag/web


#----------
# settings for *HOURLY* analysis
#----------

# hourly rate of time sampling
setenv RATE 6             # time-sample frequency [hours]  (1,2,3,4,6,8,12,24)

# usually leave the following alone
@ spd = ( 24 / $RATE )     # calculate samples per day
setenv SPD $spd            
setenv TSTRIDE 1           # TSTRIDE > 1 skips data samples in data set (e.g. = 2 uses every other sample)
setenv LATBND 15           # define tropics (-S,N)


#----------
# Software 
#----------
#
# Diagnostic package location and settings
#
# The environment variable DIAG_HOME must be set to run this script
#    It indicates where the variability package source code lives and should
#    contain the directories var_code and var_data although these can be 
#    located elsewhere by specifying below.

setenv DIAG_HOME /glade/p/work/cchen/VDIAG
setenv VARCODE $DIAG_HOME/var_code
setenv VARDATA $DIAG_HOME/var_data

setenv QUIT_ON_ERROR 1   # 0-try to continue, 1-quit
setenv OVERWRITE_EXISTING_FILES n  # y/n

#-----------------------------------------------------------------------------------------------------
#---------------------------- end of user modification -----------------------------------------------
#-----------------------------------------------------------------------------------------------------


# Check directories
if ( ! -d $DIAG_HOME ) then   # check that it exists
    echo "DIAG_HOME directory does not exist $DIAG_HOME"
    exit -1
else 
    echo DIAG_HOME is set to $DIAG_HOME
endif

if ( ! $?VARCODE ) setenv VARCODE $DIAG_HOME/var_code
if ( ! -d $VARCODE ) then
    echo "VARCODE directory $VARCODE does not exist "
    echo "The default setting is \$DIAG_HOME/var_code but user settings can override this"
    exit -1
else 
    echo VARCODE is set to $VARCODE
endif

if ( ! $?VARDATA ) setenv VARDATA $DIAG_HOME/../var_data
if ( ! -d $VARDATA ) then
    echo "VARDATA directory $VARDATA does not exist "
    echo "The default setting is \$DIAG_HOME/../var_data but user settings can override this"
    exit -1
else 
    echo VARDATA is set to $VARDATA
endif

if ( ! $?RGB ) setenv RGB $DIAG_HOME/rgb
if ( ! -d $RGB ) then
    echo "RGB directory $RGB does not exist "
    echo "The default setting is \$DIAG_HOME/rgb but user settings can override this"
    exit -1
else 
    echo RGB is set to $RGB
endif

#
# NCL / NCARG
#
if (! $?NCARG_ROOT) then
  echo ERROR: You do not have the environment
  echo variable NCARG_ROOT defined, which is used by NCL
    exit -1
endif

if ( ! -d $NCARG_ROOT ) then
    echo "NCARG_ROOT directory does not exist $NCARG_ROOT"
    exit -1
else
    echo NCARG_ROOT is set to $NCARG_ROOT
endif

set whichncl = `which ncl`
if ( $#whichncl > 0 ) then
    setenv NCL $whichncl
    echo "using ncl $NCL"
else
    echo "ERROR: ncl not found"
    exit -1
endif


# This is done for a reason: make sure all years are 4-digit 
#set evans_year = 

setenv FIRSTYR_MONTHLY $firstyr_monthly
setenv LASTYR_MONTHLY  $lastyr_monthly
setenv FIRSTYR_DAILY   $firstyr_daily
setenv LASTYR_DAILY    $lastyr_daily
setenv FIRSTYR_HOURLY  $firstyr_hourly
setenv LASTYR_HOURLY   $lastyr_hourly
setenv FIRSTYR_DIURNAL $firstyr_diurnal
setenv LASTYR_DIURNAL  $lastyr_diurnal
#setenv FIRSTYR_MONTHLY `printf "%04d" {$firstyr_monthly}`
#setenv LASTYR_MONTHLY  `printf "%04d" {$lastyr_monthly}`
#setenv FIRSTYR_DAILY   `printf "%04d" {$firstyr_daily}`
#setenv LASTYR_DAILY    `printf "%04d" {$lastyr_daily}`
#setenv FIRSTYR_HOURLY  `printf "%04d" {$firstyr_hourly}`
#setenv LASTYR_HOURLY   `printf "%04d" {$lastyr_hourly}`
#setenv FIRSTYR_DIURNAL `printf "%04d" {$firstyr_diurnal}`
#setenv LASTYR_DIURNAL  `printf "%04d" {$lastyr_diurnal}`

if ( $do_monthly ) echo "monthly: using date ranges $FIRSTYR_MONTHLY $LASTYR_MONTHLY"
if ( $do_daily   ) echo "daily  : using date ranges $FIRSTYR_DAILY   $LASTYR_DAILY"
if ( $do_hourly  ) echo "hourly : using date ranges $FIRSTYR_HOURLY  $LASTYR_HOURLY"
if ( $do_diurnal ) echo "diurnal: using date ranges $FIRSTYR_DIURNAL $LASTYR_DIURNAL"


#
# Directory to store figures, structured for web page
# DON'T CHANGE! This name is (annoyingly) hardcoded in some of the underlying scripts.
#
setenv variab_dir  $WKDIR/VDIAG_$CASENAME      


# make directories that don't already exist
if ( ! -d $WKDIR ) mkdir $WKDIR
if ( ! -d $variab_dir ) mkdir $variab_dir
if ( $?WEBDIR ) if ( ! -d $WEBDIR  ) mkdir -p $WEBDIR 


setenv OS `uname`

#  env vars still used in ncl code 
 setenv DIRPATH $WKDIR:h

# --------------------------------------------------------------------------- 
# Time to do the work
#

#cat README
cd $WKDIR

if ( $do_diurnal ) then
   $VARCODE/precip_dc.csh
endif

if ( $do_hourly  )  then
    $VARCODE/hourly.csh
endif	


if ( $do_monthly ) then

    if (! -d $WKDIR/VDIAG_{$CASENAME}/monthly/$CASENAME)  mkdir -p $WKDIR/VDIAG_{$CASENAME}/monthly/$CASENAME
    cd $WKDIR/VDIAG_{$CASENAME}/monthly/$CASENAME
   $VARCODE/monthly.csh 
endif
    

if ( $do_daily ) then
   $VARCODE/daily.csh 
endif
    
if ( $make_variab_tar ) then 

    if ( -e $variab_dir.tar ) mv -f $variab_dir.tar $variab_dir.tar_old
    cd $WKDIR
    tar cf VDIAG_{$CASENAME}.tar VDIAG_{$CASENAME}
    if ( $status != 0 ) then
	echo "ERROR $0"
	echo "trying to do:     tar cf $variab_dir.tar $variab_dir"
        exit -1
    endif

    if ( $?WEBDIR ) then
	mv $variab_dir.tar $WEBDIR
        if ( $status == 0 ) then 
	    cd $WEBDIR
	    tar xf VDIAG_{$CASENAME}.tar
	    if ( $status == 0 ) then
		ls $variab_dir.tar
    	    endif # tar xf status == 0 
        else
            ls -l $cwd/$variab_dir.tar
    	endif     # mv tar files status == 0 
    else # $?WEBDIR
        echo "no WEBDIR set so leaving tar file here"
        ls -l $cwd/$variab_dir.tar
    endif # $?WEBDIR

endif # make_variab_tar


exit 0
#
#
# ======================================================================
# end of succesful script
# ======================================================================



