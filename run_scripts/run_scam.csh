#!/bin/csh -f
#From John Truesdale. 
#use modified version for scam-clubb. Pete Bogenschutz

setenv VERSION 20091023

set silhs_enabled = "true" # Change this to "true" to use SILHS (only valid when running camclubb)

set MGVER = 2
set NUMSC = 4

alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'  # do not modify this

#CASES
# 0 = arm97
# 3 = ATEX_48hr
# 4 = BOMEX_5day
# 5 = DYCOMSrf01_4day
# 6 = DYCOMSrf02_06hr
# 7 = RICO_3day
# 8 = ARM_CC
# 12 = IOPCASE
# 13 = TWP_ICE

setenv USER_ROOT `pwd`  #change this, should be directory where
                        # this script is

#foreach num ('1012') #change this, this is case identifier in case you run a case more than once, 

foreach num ('707') # DO NOT COMMIT ANY CHANGE TO THIS NUMBER TO THE REPOSITORY! IT WILL BREAK THE NIGHTLY TESTS!

                  #so it won't overwrite the previous cases.  For example, if num = '100' and 		  
		  #simtype is 'camclubb' with 30 levels and dt = 1800 s then the case name will
		  #be called "camclubb100_L30_T1800".
		  
# series 1010 is for marine Sc
# series 1020 is for trade cu
# series 1030 is for storm tracks  
# series 1040 is for deep convection

#NOTE, below you need to set some important stuff if running IOPCASE (option 12), 
# if you are not running case 12, then just ignore below.
set dolat = -5 #change this (only if case 12), latitude to simulate (from -90 to 90)   ! 15, 180 good trade cu
set dolon = 355 #change this (only if case 12), longitude to simulation (from 0 to 360)  ! -15, 275 Sc
setenv domonth 'jan' #change this (only if running case 12), month to simulation (either 'jan' or 'july')

# Bora... lat -2 lon 118
set dodeep = 'no' # yes = want ZM deep convection, no = do not want ZM deep convection
set advmoms = 'no' # yes = advect CLUBB's moments, no = do not advect moments

MATH num2 = 1 + $num

#NOTE: Do you want to run CAM-BASE or CAM-CLUBB (or both?)
foreach simtype ('camclubb') # change this, setting to 'base' will run standard
                                    # CAM5, setting 'camclubb' will run CAM-CLUBB.  If 
				    # both are set then both will be run back-to-back, 
				    # ie. "('base' 'camclubb')

if ($simtype != 'camclubb' && $silhs_enabled == "true") then
  echo "Cannot run SILHS without SCAM-CLUBB enabled."
  exit 1
endif

# NOTE: Below, set caseid, levarr, tarray.  can be more than one values, if so will loop 
foreach caseid (0 3 4 5 6 7 8 12 13) # change this, see above (depending on case you want to simulate)
#foreach caseid (12) # IOPCASE only

foreach levarr(30)  # change this, number of levels to run
                    # options are 30, 60, 60, 90, 120, 150, 180, 210, 240
  foreach tarray(1200) # change this, the host model timestep 
                       # options are 60, 300, 600, 900, 1200, 1500, 1800, 2100, 2400 s 

setenv SCAM_TEST {$simtype}{$num}_L{$levarr}_T{$tarray}
setenv CAM_TAG cam5_2_06    # change this, this is the CAM tag you are using
#setenv CAM_TAG cam5_2_09-TRUNK
#setenv CAM_TAG clubb65_cam5_1_38
#setenv CAM_TAG clubb52_radiation

if ($CAM_TAG == clubb51_cam5_1_26_adv) set advmoms = 'yes'

set TEST_MDIR = $SCAM_TEST 
setenv SCAM_TLEVS $levarr       

######## IOP
if ($caseid == 0) then
setenv IOP_NAME  'arm97' ; setenv SCAM_IOP 'ARM 1997 (18 Jun 1997 - 18 Jul 1997)'
endif

if ($caseid == 1) then
setenv IOP_NAME  'gate'  ; setenv SCAM_IOP 'GATE IOP III (20 Aug 1974 - 18 Sep 1974)'
endif

if ($caseid == 2) then
setenv IOP_NAME  'toga' ; setenv SCAM_IOP  'TOGA-COARE IOP II (12 Dec 1992 - 1 Jan 1993)'
endif

if ($caseid == 3) then
setenv IOP_NAME  'ATEX_48hr' ; setenv SCAM_IOP 'ATEX IOP February 1969'
endif

if ($caseid == 4) then
setenv IOP_NAME  'BOMEX_5day' ; setenv SCAM_IOP 'BOMEX IOP June 1969'
endif

if ($caseid == 5) then
setenv IOP_NAME  'DYCOMSrf01_4day' ; setenv SCAM_IOP 'DYCOMS2-RF01 July 1999'
endif

if ($caseid == 6) then
setenv IOP_NAME  'DYCOMSrf02_06hr' ; setenv SCAM_IOP 'DYCOMS2-RF02 July 1999'
endif

if ($caseid == 7) then
setenv IOP_NAME  'RICO_3day' ; setenv SCAM_IOP 'RICO July 1995'
endif

if ($caseid == 8) then
setenv IOP_NAME  'ARM_CC' ; setenv SCAM_IOP 'ARM June 21 1997'
endif

if ($caseid == 9) then
setenv IOP_NAME  'ASTEX'  ; setenv SCAM_IOP 'ASTEX June 13 1992'
endif

if ($caseid == 10) then
setenv IOP_NAME  'mpace' ; setenv SCAM_IOP 'MPACE October 9 2004' 
endif

if ($caseid == 11) then
setenv IOP_NAME  'VOCALS_exper' ; setenv SCAM_IOP 'VOCALS January 1900' 
endif

if ($caseid == 12) then
setenv IOP_NAME  'IOPCASE' ; setenv SCAM_IOP 'Taken IOP from global simulation January 1900'
endif

if ($caseid == 13) then
setenv IOP_NAME  'TWP06'   ; setenv SCAM_IOP 'TWP-ICE case'
endif

#set input_name = $IOP_NAME
## Generic users' directories

cd ..

setenv SCAM_ROOT   $USER_ROOT/../run_plots                #change this, this is where output will be placed
setenv IOP_ROOT /home/pub/cam/cam_inputdata

#NOTE, Below is where code is stored
setenv CAM_DIR   $USER_ROOT/..    #change this, this is where CAM code is

######  Copy over the nessecary mods to the right directories

### testing only
#setenv MODDIR /project/amp/bogensch/scam5/mods_v2/$SCAM_TEST && if (! -e $MODDIR) mkdir $MODDIR
#rm $MODDIR/*.F90

#cp /project/amp/bogensch/scam5/mods_v2/consdrop/*.F90 $MODDIR

#cp /project/amp/bogensch/scam5/mods_38/pbuf_times/*.F90 $MODDIR

#cp /project/amp/bogensch/scam5/mods_38/energy_balance/*.F90 $MODDIR
#cp /project/amp/bogensch/scam5/mods_38/pbl_mods/*.F90 $MODDIR
#cp /project/amp/bogensch/scam5/mods_38/sfcflx_mods2/*.F90 $MODDIR
#cp /project/amp/bogensch/scam5/mods_38/ndrop_mods/*.F90 $MODDIR
#cp /project/amp/bogensch/scam5/mods_38/energy_balance2/*.F90 $MODDIR

#cp /project/amp/bogensch/scam5/mods_38/sfcflx_mods3/*.F90 $MODDIR
#cp /project/amp/bogensch/scam5/mods_38/test_b4b/*.F90 $MODDIR

#########################################################################
#### Set some case specific parameters here
#### Here the boundary layer cases use prescribed aerosols while the deep convection
#### and mixed phase cases use prognostic aerosols.  This is because the boundary layer
#### cases are so short that the aerosols do not have time to spin up.
if ($IOP_NAME  == 'arm97' || $IOP_NAME == 'gate' || $IOP_NAME == 'toga' || $IOP_NAME == 'mpace' || $IOP_NAME == 'IOPCASE' || $IOP_NAME  == 'TWP06') then
  set aero_mode = 'modal'
else
  set aero_mode = 'bulk'
endif

if ($advmoms == 'yes') then
  if ($aero_mode == 'modal') set    numadv = 34
  if ($aero_mode == 'bulk') set     numadv = 14
endif

if ($advmoms == 'no') then
  if ($aero_mode == 'modal') set    numadv = 25
  if ($aero_mode == 'bulk') set     numadv = 5
endif

#########################################################################
## Computer Stuff

# Select compiler+libraries depending on machine.

if (! $?OS_NAME)   set OS_NAME   = `uname -s`
if (! $?PROCESSOR) set PROCESSOR = `uname -m`
setenv CSMDATA $IOP_ROOT
if ($OS_NAME == "Linux") then # Linux
    if ($PROCESSOR == "i686") then # 32-bit machines (NCAR-CGD: bangkok,calgary)
	set FORTC = /usr/local/lf9562/bin/lf95
#	set FORTC = /usr/local/lf95
	setenv USER_FC lf95
	setenv NCHOME /usr/local/netcdf3-gcc-lf95
    else if ($PROCESSOR == "x86_64") then # 64-bit machines (NCAR: mineral,hurricane,(salina?))

#	set FORTC = /usr/local/lf95
#	setenv USER_FC lf95
#	setenv NCHOME /usr/local/netcdf3-gcc-lf95
 
#	set DBUG = "-debug"
        set DBUG = " "
#	set FORTC = /usr/local/pgi/linux86-64
#        setenv FORTC /usr/local/pgi/linux86-64
#	setenv USER_FC pgf90
#	setenv NCHOME /usr/local/netcdf-pgi

setenv FORTC /opt/intel/composer_xe_2013.5.192/bin/intel64/
setenv USER_FC ifort
setenv NCHOME /usr/local/netcdf-intel64/

    endif
endif    

setenv PATH $FORTC/bin:$PATH
setenv LD_LIBRARY_PATH $FORTC/lib:$LD_LIBRARY_PATH   
setenv INC_NETCDF ${NCHOME}/include
setenv LIB_NETCDF ${NCHOME}/lib
setenv NCARG_ROOT /contrib/ncarg
setenv PATH $NCARG_ROOT/bin:$PATH
setenv LD_LIBRARY_PATH ${NCHOME}/lib:$LD_LIBRARY_PATH 
set    NCL = /usr/local/bin/ncl
#set    LDFLAGS = "-llapack -lblas"

echo "***** scam_iop *****"
echo "***** Version - $VERSION *****"
echo "***** $SCAM_IOP *****"

### Run scam
# Directory structures
    
    setenv WKDIR     $SCAM_ROOT                 && if (! -e  $WKDIR) mkdir $WKDIR
    setenv WKDIR     $SCAM_ROOT/runs            && if (! -e  $WKDIR) mkdir $WKDIR
    setenv SCAM_SCR  $CAM_DIR/components/cam/bld

# CONTROL #

    setenv WKDIR     $SCAM_ROOT/runs/$SCAM_TEST && if (! -e  $WKDIR) mkdir $WKDIR
    setenv WKDIR     $WKDIR/$IOP_NAME           && if (! -e  $WKDIR) mkdir $WKDIR
    setenv SCAM_BLD  $WKDIR/bld                 && if (! -e  $SCAM_BLD) mkdir $SCAM_BLD
    setenv SCAM_EXE  $WKDIR/run                 && if (! -e  $SCAM_EXE) mkdir $SCAM_EXE
    setenv SCAM_MODS $WKDIR/mods                && if (! -e  $SCAM_MODS) mkdir $SCAM_MODS  
    setenv SCAM_LEVS $SCAM_TLEVS 
    
#    rm $WKDIR/run/*.nc

##--------------------------
## Configure
##--------------------------

    set silhs_config_option = ""
    if ($silhs_enabled == "true") then
      set silhs_config_option = "-psubcols $NUMSC"
    endif

    cd $SCAM_BLD

    # Set options
if ($aero_mode == 'bulk') then
    set chem_opt = "-chem none"
else
    set chem_opt = ""
endif
if ($simtype == 'camclubb') then
    set clubb_opt="-clubb_sgs"
else
    set clubb_opt=""
endif
if ($silhs_enabled == "true") then
    set cppdefs = "-DDISABLE_TIMERS -DUWM_MISC -DSILHS"
else
    set cppdefs = "-DDISABLE_TIMERS -DUWM_MISC"
endif
    $SCAM_SCR/configure -phys cam5 $chem_opt -nlev $SCAM_TLEVS -dyn eul -res 64x128 $silhs_config_option $clubb_opt -nospmd -nosmp -cppdefs "$cppdefs" -scam -usr_src $SCAM_MODS -fc $USER_FC -cc gcc $DBUG -ldflags "-mkl -Mnobounds" -microphys mg$MGVER

##--------------------------
## compile
##--------------------------

   echo " -- Compile"
   gmake -j6 >&! MAKE.out || echo 'ERROR: Compile failed for' $SCAM_TEST '- exiting scam2web'&& exit 1

##--------------------------
## Build the namelist
##---------------------------

# Input something for L26/27/30/31 runs
if ($SCAM_LEVS == 26) set icfile = $CSMDATA/atm/cam2/inic/gaus/cami_0000-09-01_64x128_L26_c030918.nc
#if ($SCAM_LEVS == 27) set icfile = /project/amp/rneale/data/inic_4cam/cami_0000-09-01_64x128_L27_c030918.nc

if ($SCAM_LEVS == 30) set icfile = $CSMDATA/csm/cami_0000-01-01_64x128_L30_c090102.nc
#if ($SCAM_LEVS == 30) set icfile = /project/convection/sungsup/scam/data_I3D/profile/cami_0000-01-01_64x128_L30_c090102.nc
if ($SCAM_LEVS == 60) set icfile = $CSMDATA/csm/cami_0000-01-01_64x128_L60_c090102.nc
if ($SCAM_LEVS == 90) set icfile = $CSMDATA/csm/cami_0000-01-01_64x128_L90_c090102.nc
if ($SCAM_LEVS == 120) set icfile = $CSMDATA/csm/cami_0000-01-01_64x128_L120_c090102.nc
if ($SCAM_LEVS == 150) set icfile = $CSMDATA/csm/cami_0000-01-01_64x128_L150_c090102.nc
if ($SCAM_LEVS == 180) set icfile = $CSMDATA/csm/cami_0000-01-01_64x128_L180_c090102.nc
if ($SCAM_LEVS == 210) set icfile = $CSMDATA/csm/cami_0000-01-01_64x128_L210_c090102.nc
if ($SCAM_LEVS == 240) set icfile = $CSMDATA/csm/cami_0000-01-01_64x128_L240_c090102.nc

#if ($IOP_NAME == 'STORM') set icfile = /project/amp/bogensch/scam4alpha/iop/cam5_1_13_hourly.cam.i.0001-01-01-00000.nc

if ($IOP_NAME == 'IOPCASE') then
  if ($domonth == 'july') set icfile = $IOP_ROOT/iop/camclubb_iop2.cam.i.0001-07-01-00000.nc
  if ($domonth == 'july' && $SCAM_LEVS == 60) set icfile = $IOP_ROOT/iop/july_60.nc
  if ($domonth == 'jan') set icfile = $IOP_ROOT/iop/camclubb_iop.cam.i.0001-01-01-00000.nc
endif

echo icfile $icfile

if ($IOP_NAME == 'IOPCASE') then
  if ($domonth == 'july') set iopfile = $IOP_ROOT/iop/camclubb_iop.cam.h1.0001-07-01-00000.nc
  if ($domonth == 'jan') set iopfile = $IOP_ROOT/iop/camclubb_iop.cam.h1.0001-01-01-00000.nc
endif

#  Data extracted from the IOP file.

if ($IOP_NAME != 'IOPCASE') then

# Set this to the directory containing the IOP files you want to use.
set iopfile = $IOP_ROOT/iop/${IOP_NAME}_4scam.nc

set lat       = `ncks -H -s "%5.2f\n" -v lat    ${iopfile}`
set lon       = `ncks -H -s "%5.2f\n" -v lon    ${iopfile}`

set lon_360 = `echo $lon | cut -c1`                 # Switcheroo if lon is -ve.
echo lon_360 $lon_360

#if ($lon_360 == -) set lon = 360. + lon
if ($lon_360 == -) then
  MATH lon = 360 + $lon
endif
echo aerosols $aero_mode
echo lon $lon

else

set lat = $dolat   # STORM = -60
set lon = $dolon  # STORM = 250

endif

#set timestep, default is 1200
set dtime = $tarray

if ($IOP_NAME != 'IOPCASE') then

set basedate  = `ncks -H -s "%i\n"    -v bdate  ${iopfile}`
set basesecs  = `ncks -H -s "%i\n"    -v tsec   ${iopfile} | head -1`
set steplen   = `ncks -H -s "%i\n"    -v tsec   ${iopfile} | sed -n '2p'`
@ steplen = $steplen - $basesecs
set endstep   = `ncks -H -s "%i\n"    -v tsec   ${iopfile} | sed -n '$='`

#set model timesteps based on size of IOP times and stride (step length) in iop file
echo timestep $dtime
echo basedate $basedate
echo basesecs $basesecs
echo steplen $steplen
echo endstep $endstep

@ stop_n  = (($endstep - 3) * $steplen) / $dtime

# basedate (bdate in nc iop file) may or may not include the century.
# A bit of logic to find out what it is.

#if ($basedate < 10000000) then # Does not contain century.
#    if ($basedate >= 500000) set start_ymd = 19{$basedate}       # (19)50-(19)99.
#    if ($basedate <  500000) set start_ymd = 20`printf "%06d" {$basedate}`  # 20{00)-(20)50.  
#    if ($basedate <  500000) set start_ymd = 190{$basedate}   #hack for MPACE (assume 1904? not 2004)
#else
   set start_ymd = {$basedate} # Contains century info - just copy.
#endif

else

if ($domonth == 'jan') set start_ymd = '10101'
if ($domonth == 'july') set start_ymd = '10701'
set basesecs = 0
set stop_n=2000

endif

set clmin = 'scam_01.clm2.r.0001-07-03-80400.nc'

if ($dtime >= 900) set nhtfrqval = 1
if ($dtime == 600) set nhtfrqval = 2
if ($dtime == 300) set nhtfrqval = 4
if ($dtime == 150) set nhtfrqval = 8
if ($dtime == 60) set nhtfrqval = 20

#echo basedate $basedate
echo scmlat $lat
echo scmlon $lon
echo startymd $start_ymd
echo starttod $basesecs
echo stopn $stop_n
echo nhtfrq $nhtfrqval

#set start_ymd = '19010701'
#set stopn = 1000
#endif

if ($simtype == 'camclubb') then
  set shallow_in = 'CLUBB_SGS'
  set macrop_in = 'CLUBB_SGS'
  set eddy_scheme_in = 'CLUBB_SGS'
else
  set shallow_in = 'UW'
  set macrop_in = 'park'
  set eddy_scheme_in = 'diag_TKE'
endif

if ($dodeep == 'yes') then 
  set deep_in = 'ZM'
else
  set deep_in = 'off'
endif

if ($IOP_NAME == 'IOPCASE') then
  set mfilt_in=2100
else
  set mfilt_in = 5760
endif 

#set stop_n=3

if ($simtype == "camclubb") then

  ## A list of CLUBB variables

  set clubb_vars_zt_list = "'thlm', 'thvm', 'rtm', 'rcm', 'rvm', 'um', 'vm', 'um_ref','vm_ref','ug', 'vg', 'cloud_frac', 'cloud_cover', 'rcm_in_layer', 'rcm_in_cloud', 'p_in_Pa', 'exner', 'rho_ds_zt', 'thv_ds_zt', 'Lscale', 'Lscale_pert_1', 'Lscale_pert_2', 'T_in_K', 'rel_humidity', 'wp3', 'wpthlp2', 'wp2thlp', 'wprtp2', 'wp2rtp', 'Lscale_up', 'Lscale_down', 'tau_zt', 'Kh_zt', 'wp2thvp', 'wp2rcp', 'wprtpthlp', 'sigma_sqd_w_zt', 'rho', 'radht', 'radht_LW', 'radht_SW', 'Ncm', 'Nc_in_cloud', 'Nc_activated', 'snowslope', 'sed_rcm', 'rsat', 'rsati', 'diam', 'mass_ice_cryst', 'rcm_icedfs', 'u_T_cm', 'rtm_bt', 'rtm_ma', 'rtm_ta', 'rtm_mfl', 'rtm_tacl', 'rtm_cl', 'rtm_forcing', 'rtm_sdmp','rtm_mc', 'rtm_pd', 'rvm_mc', 'rcm_mc', 'rcm_sd_mg_morr', 'thlm_bt', 'thlm_ma', 'thlm_ta', 'thlm_mfl', 'thlm_tacl', 'thlm_cl', 'thlm_forcing', 'thlm_sdmp','thlm_mc', 'thlm_old', 'thlm_without_ta', 'thlm_mfl_min', 'thlm_mfl_max', 'thlm_enter_mfl', 'thlm_exit_mfl', 'rtm_old', 'rtm_without_ta', 'rtm_mfl_min', 'rtm_mfl_max', 'rtm_enter_mfl', 'rtm_exit_mfl', 'um_bt', 'um_ma', 'um_gf', 'um_cf', 'um_ta', 'um_f', 'um_sdmp', 'um_ndg', 'vm_bt', 'vm_ma', 'vm_gf', 'vm_cf', 'vm_ta', 'vm_f', 'vm_sdmp', 'vm_ndg', 'wp3_bt', 'wp3_ma', 'wp3_ta', 'wp3_tp', 'wp3_ac', 'wp3_bp1', 'wp3_pr_turb', 'wp3_pr_dfsn', 'wp3_pr1', 'wp3_pr2', 'wp3_dp1', 'wp3_cl', 'mixt_frac', 'w_1', 'w_2', 'varnce_w_1', 'varnce_w_2', 'thl_1', 'thl_2', 'varnce_thl_1', 'varnce_thl_2', 'rt_1', 'rt_2', 'varnce_rt_1', 'varnce_rt_2', 'rc_1', 'rc_2', 'rsatl_1', 'rsatl_2', 'cloud_frac_1', 'cloud_frac_2', 'a3_coef_zt', 'wp3_on_wp2_zt', 'chi_1', 'chi_2', 'stdev_chi_1', 'stdev_chi_2', 'stdev_eta_1', 'stdev_eta_2', 'covar_chi_eta_1', 'covar_chi_eta_2', 'corr_chi_eta_1', 'corr_chi_eta_2', 'corr_rt_thl_1', 'crt_1', 'crt_2', 'cthl_1', 'cthl_2', 'precip_frac', 'precip_frac_1', 'precip_frac_2', 'Ncnm', 'wp2_zt', 'thlp2_zt', 'wpthlp_zt', 'wprtp_zt', 'rtp2_zt', 'rtpthlp_zt', 'up2_zt', 'vp2_zt', 'upwp_zt', 'vpwp_zt'"

  set clubb_vars_zm_list = "'wp2', 'rtp2', 'thlp2', 'rtpthlp', 'wprtp', 'wpthlp', 'wp4', 'up2', 'vp2', 'wpthvp', 'rtpthvp', 'thlpthvp', 'tau_zm', 'Kh_zm', 'wprcp', 'rc_coef', 'wm_zm', 'thlprcp', 'rtprcp', 'rcp2', 'upwp', 'vpwp', 'rho_zm', 'sigma_sqd_w', 'Skw_velocity', 'gamma_Skw_fnc', 'C6rt_Skw_fnc', 'C6thl_Skw_fnc', 'C7_Skw_fnc', 'C1_Skw_fnc', 'a3_coef', 'wp3_on_wp2', 'rcm_zm', 'rtm_zm', 'thlm_zm', 'cloud_frac_zm', 'rho_ds_zm', 'thv_ds_zm', 'em', 'mean_w_up', 'mean_w_down', 'shear', 'wp3_zm', 'Frad', 'Frad_LW', 'Frad_SW', 'Frad_LW_up', 'Frad_SW_up', 'Frad_LW_down', 'Frad_SW_down', 'Fprec', 'Fcsed', 'wp2_bt', 'wp2_ma', 'wp2_ta', 'wp2_ac', 'wp2_bp', 'wp2_pr1', 'wp2_pr2', 'wp2_pr3', 'wp2_dp1', 'wp2_dp2', 'wp2_cl', 'wp2_pd', 'wp2_sf', 'vp2_bt', 'vp2_ma', 'vp2_ta', 'vp2_tp', 'vp2_dp1', 'vp2_dp2', 'vp2_pr1', 'vp2_pr2', 'vp2_cl', 'vp2_pd', 'vp2_sf', 'up2_bt', 'up2_ma', 'up2_ta', 'up2_tp', 'up2_dp1', 'up2_dp2', 'up2_pr1', 'up2_pr2', 'up2_cl', 'up2_pd', 'up2_sf', 'wprtp_bt', 'wprtp_ma', 'wprtp_ta', 'wprtp_tp', 'wprtp_ac', 'wprtp_bp', 'wprtp_pr1', 'wprtp_pr2', 'wprtp_pr3', 'wprtp_dp1', 'wprtp_mfl', 'wprtp_cl', 'wprtp_sicl', 'wprtp_pd', 'wprtp_forcing', 'wprtp_mc', 'wpthlp_bt', 'wpthlp_ma', 'wpthlp_ta', 'wpthlp_tp', 'wpthlp_ac', 'wpthlp_bp', 'wpthlp_pr1', 'wpthlp_pr2', 'wpthlp_pr3', 'wpthlp_dp1', 'wpthlp_mfl', 'wpthlp_cl', 'wpthlp_sicl', 'wpthlp_forcing', 'wpthlp_mc', 'rtp2_bt', 'rtp2_ma', 'rtp2_ta', 'rtp2_tp', 'rtp2_dp1', 'rtp2_dp2', 'rtp2_cl', 'rtp2_pd', 'rtp2_sf', 'rtp2_forcing', 'rtp2_mc', 'thlp2_bt', 'thlp2_ma', 'thlp2_ta', 'thlp2_tp', 'thlp2_dp1', 'thlp2_dp2', 'thlp2_cl', 'thlp2_pd', 'thlp2_sf', 'thlp2_forcing', 'thlp2_mc', 'rtpthlp_bt', 'rtpthlp_ma', 'rtpthlp_ta', 'rtpthlp_tp1', 'rtpthlp_tp2', 'rtpthlp_dp1', 'rtpthlp_dp2', 'rtpthlp_cl', 'rtpthlp_sf', 'rtpthlp_forcing', 'rtpthlp_mc', 'wpthlp_entermfl', 'wpthlp_exit_mfl', 'wprtp_enter_mfl', 'wprtp_exit_mfl', 'wpthlp_mfl_min', 'wpthlp_mfl_max', 'wprtp_mfl_min', 'wprtp_mfl_max'"

  if ($silhs_enabled == "true") then
    set subcol_scheme_nl_entry = "subcol_scheme = 'SILHS',"
    set use_subcol_microp_nl = ".true."
  else
    set subcol_scheme_nl_entry = ""
    set use_subcol_microp_nl = ".false."
  endif
else
  set clubb_vars_zt_list = "'RHO'"
  set clubb_vars_zm_list = ""
  set subcol_scheme_nl_entry = ""
  set use_subcol_microp_nl = ".false."
endif

## temporary namelist
cd $SCAM_BLD                      || echo "cd $blddir failed" && exit 1
cat <<EOF >! tmp_namelistfile
&camexp 
    case_name            = '$SCAM_TEST', 
    stop_option          = 'nsteps',
    stop_n               = $stop_n, 
    scmlat               = $lat,
    scmlon               = $lon, 
    single_column        = .true., 
    history_budget       = .false.,
    history_aerosol      = .true.,
    ndens                = 1,1,1,1,1,1,
    scm_iop_srf_prop     = .true.,
    print_energy_errors  = .true., 
    ncdata               = '$icfile',
    iopfile              = '$iopfile',
    start_ymd            = $start_ymd, 
    start_tod            = $basesecs,
    dtime                = $dtime,
    prescribed_ozone_file = '$IOP_ROOT/atm/cam2/ozone/ozone_1.9x2.5_L26_2000clim_c090803.nc' ,
    prescribed_ozone_cycle_yr = 2000 ,
    prescribed_ozone_type = 'CYCLICAL',
    prescribed_ozone_name = 'O3' ,
    prescribed_aero_model = '$aero_mode',
    nhtfrq               = $nhtfrqval,$nhtfrqval 
    mfilt                = $mfilt_in,$mfilt_in
    deep_scheme          = '$deep_in',
    use_subcol_microp    = $use_subcol_microp_nl,
    $subcol_scheme_nl_entry
    microp_scheme        = 'MG',
    micro_mg_version     = $MGVER,
    micro_mg_sub_version = 0,
    ice_supersat         = .true.,
    macrop_scheme        = '$macrop_in',
    shallow_scheme       = '$shallow_in',
    eddy_scheme          = '$eddy_scheme_in',
    fincl1               = 'CLDST',
       'ICWMR','ICIMR','FREQL','FREQI','LANDFRAC','CDNUMC','FICE','WSUB','CCN3','ICLDIWP',
       'CDNUMC', 'AQSNOW',  'WSUB', 'CCN3', 'FREQI', 'FREQL', 'FREQR', 'FREQS', 'CLDLIQ', 'CLDICE',
       'FSDS', 'FLDS','AREL','AREI','NSNOW','QSNOW','DSNOW',
       'FLNT','FLNTC','FSNT','FSNTC','FSNS','FSNSC','FLNT','FLNTC','QRS','QRSC','QRL','QRLC',
       'LWCF','SWCF', 'NCAI', 'NCAL', 'NIHF','NIDEP','NIIMM','NIMEY','ICLDIWP','ICLDTWP', 'CONCLD',
       'QCSEVAP', 'QISEVAP', 'QVRES', 'CMELIQ', 'CMEIOUT', 'EVAPPREC', 'EVAPSNOW', 'TAQ',
       'ICLMRCU', 'ICIMRCU' ,'ICWMRSH' ,'ICWMRDP', 'ICLMRTOT' , 'ICIMRTOT' , 'SH_CLD' ,  'DP_CLD',
       'ICWMRST', 'ICIMRST', 'ADRAIN','ADSNOW','WSUBI','DIVQ3D','DIVT3D','PRODPREC','PRECT',
       'INEGCLPTEND', 'LNEGCLPTEND', 'VNEGCLPTEND', $clubb_vars_zt_list, $clubb_vars_zm_list
    scm_clubb_iop_name   =  '$IOP_NAME'
/
EOF

if ($simtype == "camclubb") then
  cat <<EOF >>tmp_namelistfile
&clubb_his_nl
    clubb_history        = .true.,
    clubb_rad_history    = .false.,
    clubb_vars_zt        = $clubb_vars_zt_list
    clubb_vars_zm        = $clubb_vars_zm_list 
/
EOF
  if ($silhs_enabled == "true") then
    cat <<EOF >>tmp_namelistfile
&subcol_SILHS_nl
  subcol_SILHS_weight         = .true.,
  subcol_SILHS_numsubcol      = $NUMSC,
  subcol_SILHS_corr_file_path = ' $CAM_DIR/components/cam/src/physics/silhs/',
  subcol_SILHS_corr_file_name = 'arm_97',
  subcol_silhs_use_clear_col = .false.,
  subcol_silhs_q_to_micro = .true.,
  subcol_silhs_n_to_micro = .true.,
  subcol_silhs_constrainmn = .false., 
  subcol_silhs_ncnp2_on_ncnm2 = 0.55,
  hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip = 1.0,
  hmp2_ip_on_hmm2_ip_ratios%nrp2_ip_on_nrm2_ip = 1.0,
  hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip = 1.0,
  hmp2_ip_on_hmm2_ip_ratios%nip2_ip_on_nim2_ip = 1.0,
  hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip = 1.0,
  hmp2_ip_on_hmm2_ip_ratios%nsp2_ip_on_nsm2_ip = 1.0,
/
EOF
  endif
endif

#  subcol_silhs_c7 = 0.5,
#  subcol_silhs_c8 = 3.0,
#  subcol_silhs_c11 = 0.65,
#  subcol_silhs_c11b = 0.53,
#  subcol_silhs_gamma_coef = 0.32,
#  subcol_silhs_mult_coef = 1.5,
#  subcol_silhs_mu = 0.001

$SCAM_SCR/build-namelist -test -v -runtype startup -infile tmp_namelistfile \
    || echo "build-namelist failed" && exit 1

#### RUN

mv -f $SCAM_BLD/docn.stream.txt $SCAM_EXE

/bin/mv *_in $SCAM_EXE
cd $SCAM_EXE
echo " -- Run SCAM in $SCAM_EXE"
$SCAM_BLD/cam >&! scam_output.txt
echo "*** Finished SCAM for - $SCAM_TEST"
echo ""

echo "** FINISHED **"

  end
  end

  end
  
#  if ($simtype == 'camclubb') then
#    MATH num = 1 + $num
#  endif
  
end

end 

exit 0

