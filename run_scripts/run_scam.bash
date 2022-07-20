#!/bin/bash
#
# How to use:
#  Pass the case name you wish to run as the last argument (this is = $1)
#  e.g. "./run_scam.bash atex"
#
#  To run with SILHS sampling enabled, add -s or --silhs before the case name
#  e.g "./run_scam.bash -s atex"
# 
##########################################################
# Variables
CASE="scam_${@: -1}"
CASEROOT="/home/$USER/projects/scratch/$CASE"
MACH="nelson"
COMPSET="FSCAM"
RES="T42_T42"
QUEUE="regular"  
WALL_TIME="00:15:00" # H:MM:SS

# Configuration parameters
NUMSC=4
MGVER=2 # Currently "1" and "2" are allowed
SILHS_ENABLED=false
subcol_scheme_namelist_value="off"
use_subcol_microp_namelist_value=".false."
microp_uniform_namelist_value=".false."

# Note that we use `"$@"' to let each command-line parameter expand to a
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
# This also expects gnu-getopt as opposed to BSD getopt.
# Make sure you have gnu-getopt installed and it is before BSD getopt in your PATH.
TEMP=`getopt -o s --long silhs -n 'compile.bash' -- "$@"`

if [ $? != 0 ] ; then echo "Run with -h for help." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
	case "$1" in
		-s|--silhs) # Indicates we want to run with SILHS enabled

      # Set silhs flag to true
      SILHS_ENABLED=true

      # Ensure use_subcol_microp and microp_uniform true in namelist, required for silhs
      use_subcol_microp_namelist_value=".true."
      microp_uniform_namelist_value=".true."
      
      # Set subcol_scheme in the namelist to SILHS
      subcol_scheme_namelist_value="SILHS"

			shift ;;
		--) shift ; break ;;
		*) echo "Something bad happened!" ; exit 1 ;;
	esac
done

# Obtain the CAM source directory. This solution is from Stack Overflow.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
CAM_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

# Set up case
echo "----- Case Setup -----"
cd "$CAM_DIR/cime/scripts"
./create_newcase --case "test_$CASE" --mach "$MACH" --compset "$COMPSET" --res "$RES" --user-mods-dir "$CASE" --run-unsupported || {
    echo "Error creating new case" >> /dev/stderr
    exit 1
}
#cd "$CASEROOT"
cd "test_$CASE"

#./xmlchange JOB_QUEUE="$QUEUE"
#./xmlchange JOB_WALLCLOCK_TIME="$WALL_TIME"

./case.setup || exit 1

# Compile configuration 
echo "----- Compile Configuration -----"

if [ "$SILHS_ENABLED" = true ] ; then
  # Setup includes -silhs, -psubcols $NUMSC, and -DSILHS to enable SILHS
  ./xmlchange CAM_CONFIG_OPTS="-clubb_sgs -rad rrtmg -silhs -silent -microphys mg$MGVER -psubcols $NUMSC -cppdefs '-DUWM_MISC -DSILHS'"
else
  ./xmlchange CAM_CONFIG_OPTS="-clubb_sgs -rad rrtmg -silent -microphys mg$MGVER -cppdefs '-DUWM_MISC'"
fi

#./xmlchange DEBUG="FALSE" # Set to TRUE for run-time debugging

# Run configuration
echo "----- Run configuration -----"
#./xmlchange RUN_STARTDATE="0000-01-01"
#./xmlchange STOP_OPTION=nmonths
#./xmlchange STOP_N=600
./xmlchange DOUT_S="FALSE" # Turn on/off short-term archiving

# Compile
echo "----- Compile -----"
./case.build || { echo "Error building CAM" >> /dev/stderr; exit 1; }

# A list of CLUBB variables
clubb_vars_zt_list="'thlm', 'thvm', 'rtm', 'rcm', 'rvm', 'um', 'vm', 'um_ref','vm_ref','ug', 'vg', 'cloud_frac', 'cloud_cover', 'rcm_in_layer', 'rcm_in_cloud', 'p_in_Pa', 'exner', 'rho_ds_zt', 'thv_ds_zt', 'Lscale', 'Lscale_pert_1', 'Lscale_pert_2', 'T_in_K', 'rel_humidity', 'wp3', 'wpthlp2', 'wp2thlp', 'wprtp2', 'wp2rtp', 'Lscale_up', 'Lscale_down', 'tau_zt', 'Kh_zt', 'wp2thvp', 'wp2rcp', 'wprtpthlp', 'sigma_sqd_w_zt', 'rho', 'radht', 'radht_LW', 'radht_SW', 'Ncm', 'Nc_in_cloud', 'Nc_activated', 'snowslope', 'sed_rcm', 'rsat', 'rsati', 'diam', 'mass_ice_cryst', 'rcm_icedfs', 'u_T_cm', 'rtm_bt', 'rtm_ma', 'rtm_ta', 'rtm_mfl', 'rtm_tacl', 'rtm_cl', 'rtm_forcing', 'rtm_sdmp','rtm_mc', 'rtm_pd', 'rvm_mc', 'rcm_mc', 'rcm_sd_mg_morr', 'thlm_bt', 'thlm_ma', 'thlm_ta', 'thlm_mfl', 'thlm_tacl', 'thlm_cl', 'thlm_forcing', 'thlm_sdmp','thlm_mc', 'thlm_old', 'thlm_without_ta', 'thlm_mfl_min', 'thlm_mfl_max', 'thlm_enter_mfl', 'thlm_exit_mfl', 'rtm_old', 'rtm_without_ta', 'rtm_mfl_min', 'rtm_mfl_max', 'rtm_enter_mfl', 'rtm_exit_mfl', 'um_bt', 'um_ma', 'um_gf', 'um_cf', 'um_ta', 'um_f', 'um_sdmp', 'um_ndg', 'vm_bt', 'vm_ma', 'vm_gf', 'vm_cf', 'vm_ta', 'vm_f', 'vm_sdmp', 'vm_ndg', 'wp3_bt', 'wp3_ma', 'wp3_ta', 'wp3_tp', 'wp3_ac', 'wp3_bp1', 'wp3_pr_turb', 'wp3_pr_dfsn', 'wp3_pr1', 'wp3_pr2', 'wp3_dp1', 'wp3_cl', 'mixt_frac', 'w_1', 'w_2', 'varnce_w_1', 'varnce_w_2', 'thl_1', 'thl_2', 'varnce_thl_1', 'varnce_thl_2', 'rt_1', 'rt_2', 'varnce_rt_1', 'varnce_rt_2', 'rc_1', 'rc_2', 'rsatl_1', 'rsatl_2', 'cloud_frac_1', 'cloud_frac_2', 'a3_coef_zt', 'wp3_on_wp2_zt', 'chi_1', 'chi_2', 'stdev_chi_1', 'stdev_chi_2', 'stdev_eta_1', 'stdev_eta_2', 'covar_chi_eta_1', 'covar_chi_eta_2', 'corr_chi_eta_1', 'corr_chi_eta_2', 'corr_rt_thl_1', 'crt_1', 'crt_2', 'cthl_1', 'cthl_2', 'precip_frac', 'precip_frac_1', 'precip_frac_2', 'Ncnm', 'wp2_zt', 'thlp2_zt', 'wpthlp_zt', 'wprtp_zt', 'rtp2_zt', 'rtpthlp_zt', 'up2_zt', 'vp2_zt', 'upwp_zt', 'vpwp_zt', 'C11_Skw_fnc'"
clubb_vars_zm_list="'wp2', 'rtp2', 'thlp2', 'rtpthlp', 'wprtp', 'wpthlp', 'wp4', 'up2', 'vp2', 'wpthvp', 'rtpthvp', 'thlpthvp', 'tau_zm', 'Kh_zm', 'wprcp', 'wm_zm', 'thlprcp', 'rtprcp', 'rcp2', 'upwp', 'vpwp', 'rho_zm', 'sigma_sqd_w', 'Skw_velocity', 'gamma_Skw_fnc', 'C6rt_Skw_fnc', 'C6thl_Skw_fnc', 'C7_Skw_fnc', 'C1_Skw_fnc', 'a3_coef', 'wp3_on_wp2', 'rcm_zm', 'rtm_zm', 'thlm_zm', 'cloud_frac_zm', 'rho_ds_zm', 'thv_ds_zm', 'em', 'mean_w_up', 'mean_w_down', 'shear', 'wp3_zm', 'Frad', 'Frad_LW', 'Frad_SW', 'Frad_LW_up', 'Frad_SW_up', 'Frad_LW_down', 'Frad_SW_down', 'Fprec', 'Fcsed', 'wp2_bt', 'wp2_ma', 'wp2_ta', 'wp2_ac', 'wp2_bp', 'wp2_pr1', 'wp2_pr2', 'wp2_pr3', 'wp2_dp1', 'wp2_dp2', 'wp2_cl', 'wp2_pd', 'wp2_sf', 'vp2_bt', 'vp2_ma', 'vp2_ta', 'vp2_tp', 'vp2_dp1', 'vp2_dp2', 'vp2_pr1', 'vp2_pr2', 'vp2_cl', 'vp2_pd', 'vp2_sf', 'up2_bt', 'up2_ma', 'up2_ta', 'up2_tp', 'up2_dp1', 'up2_dp2', 'up2_pr1', 'up2_pr2', 'up2_cl', 'up2_pd', 'up2_sf', 'wprtp_bt', 'wprtp_ma', 'wprtp_ta', 'wprtp_tp', 'wprtp_ac', 'wprtp_bp', 'wprtp_pr1', 'wprtp_pr2', 'wprtp_pr3', 'wprtp_dp1', 'wprtp_mfl', 'wprtp_cl', 'wprtp_sicl', 'wprtp_pd', 'wprtp_forcing', 'wprtp_mc', 'wpthlp_bt', 'wpthlp_ma', 'wpthlp_ta', 'wpthlp_tp', 'wpthlp_ac', 'wpthlp_bp', 'wpthlp_pr1', 'wpthlp_pr2', 'wpthlp_pr3', 'wpthlp_dp1', 'wpthlp_mfl', 'wpthlp_cl', 'wpthlp_sicl', 'wpthlp_forcing', 'wpthlp_mc', 'rtp2_bt', 'rtp2_ma', 'rtp2_ta', 'rtp2_tp', 'rtp2_dp1', 'rtp2_dp2', 'rtp2_cl', 'rtp2_pd', 'rtp2_sf', 'rtp2_forcing', 'rtp2_mc', 'thlp2_bt', 'thlp2_ma', 'thlp2_ta', 'thlp2_tp', 'thlp2_dp1', 'thlp2_dp2', 'thlp2_cl', 'thlp2_pd', 'thlp2_sf', 'thlp2_forcing', 'thlp2_mc', 'rtpthlp_bt', 'rtpthlp_ma', 'rtpthlp_ta', 'rtpthlp_tp1', 'rtpthlp_tp2', 'rtpthlp_dp1', 'rtpthlp_dp2', 'rtpthlp_cl', 'rtpthlp_sf', 'rtpthlp_forcing', 'rtpthlp_mc', 'wpthlp_enter_mfl', 'wpthlp_exit_mfl', 'wprtp_enter_mfl', 'wprtp_exit_mfl', 'wpthlp_mfl_min', 'wpthlp_mfl_max', 'wprtp_mfl_min', 'wprtp_mfl_max', 'Richardson_num', 'shear_sqd'"


# To shut off SILHS, set use_subcol_microp = .false. and microp_uniform = .false..
# To turn on the Zhang-McFarland scheme, set deep_scheme = 'ZM'.
# After shutting off SILHS, to retune the high clouds, adjust 
# micro_mg_dcs, micro_mg_berg_eff_factor, cldfrc2m_rhmini, cldfrc2m_rhmaxi.
cat >> user_nl_cam << EOF
dtime = 300
nhtfrq = 1,1
mfilt = 5000,5000,5000
ndens = 2,2,2,2,2,2
history_budget = .true.
microp_scheme = 'MG'
micro_mg_version = $MGVER
micro_mg_sub_version = 0
!micro_mg_num_steps = 1
!micro_mg_dcs = 390e-6
macrop_scheme = 'CLUBB_SGS'
eddy_scheme = 'CLUBB_SGS'
shallow_scheme = 'CLUBB_SGS'
deep_scheme = 'ZM'
subcol_scheme = '$subcol_scheme_namelist_value'
use_subcol_microp = $use_subcol_microp_namelist_value
microp_uniform = $microp_uniform_namelist_value
history_amwg = .true.
history_vdiag = .false.
clubb_do_adv = .false.
clubb_expldiff = .false.
clubb_rainevap_turb = .false.
clubb_cloudtop_cooling = .false.
clubb_l_use_C7_Richardson = .true.
clubb_l_use_C11_Richardson = .true.
clubb_l_brunt_vaisala_freq_moist = .true.
clubb_l_use_thvm_in_bv_freq = .true.
clubb_l_rcm_supersat_adj = .false.
clubb_l_stability_correct_tau_zm = .false.
clubb_C14 = 1.0
fincl1 = 'U:A','PS:A','T:A','V:A','OMEGA:A','Z3:A','PRECT:A',
'CLDLIQ:A', 'CLDICE:A', 'LWCF:A', 'SWCF:A', 'FLUT:A',
'TMQ:A', 'PRECC:A', 'PRECL:A', 'CME:A', 'PRODPREC:A',
'EVAPPREC:A','EVAPSNOW:A','ICWMRST:A','ICIMRST:A','PRAO:A',
'PRCO:A','QCSEVAP:A','QISEVAP:A','QVRES:A','CMEIOUT:A','VTRMI:A',
'VTRMC:A','QCSEDTEN:A','QISEDTEN:A','MNUCCCO:A','MNUCCTO:A',
'MNUCCDO:A','MNUCCDOhet:A','MSACWIO:A','PSACWSO:A','BERGSO:A',
'BERGO:A','MELTO:A','HOMOO:A','QCRESO:A','PRCIO:A','PRAIO:A',
'MELTSDT:A','FRZRDT:A','ADRAIN:A','ADSNOW:A','FREQR:A','FREQS:A',
'PE:A','APRL:A','PEFRAC:A','VPRCO:A','VPRAO:A','RACAU:A',
'QIRESO:A','QCRESO:A','PRACSO:A','MPDT:A','MPDQ:A','MPDLIQ:A',
'MPDICE:A','INEGCLPTEND', 'LNEGCLPTEND', 'VNEGCLPTEND',
'QCRAT:A', $clubb_vars_zt_list,$clubb_vars_zm_list,
'SL', 'Q', 'RHW', 'QRS', 'QRL', 'HR', 'FDL'
!, 'SILHS_CLUBB_PRECIP_FRAC','SILHS_CLUBB_ICE_SS_FRAC'
ncdata='/home/pub/cam_inputdata/atm/cam/inic/gaus/cami_0000-01-01_64x128_L32_c170510.nc' 
fincl2 = 'CLDTOT', 'CLDST','CDNUMC','CLDLIQ','CLDICE','FLUT',
'LWCF','SWCF','PRECT'
subcol_SILHS_weight = .true.
subcol_SILHS_numsubcol = $NUMSC
subcol_SILHS_corr_file_name = 'arm_97'
subcol_silhs_q_to_micro = .true. ! if .false. gridbox means are used instead of sample points
subcol_silhs_n_to_micro = .true. ! if .false. gridbox means are used instead of sample points
subcol_silhs_use_clear_col = .false.
subcol_SILHS_constrainmn = .false.
subcol_silhs_ncnp2_on_ncnm2 = 0.05,
subcol_silhs_hmp2_ip_on_hmm2_ip_intrcpt%rr = 1.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_intrcpt%Nr = 1.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_intrcpt%rs = 1.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_intrcpt%Ns = 1.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_intrcpt%ri = 1.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_intrcpt%Ni = 1.0
subcol_silhs_hmp2_ip_on_hmm2_ip_slope%rr = 0.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_slope%Nr = 0.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_slope%rs = 0.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_slope%Ns = 0.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_slope%ri = 0.0,
subcol_silhs_hmp2_ip_on_hmm2_ip_slope%Ni = 0.0
l_lh_importance_sampling = .true.
l_Lscale_vert_avg = .false.
l_lh_straight_mc = .false.
l_lh_clustered_sampling = .true.
l_rcm_in_cloud_k_lh_start = .true.
l_random_k_lh_start = .false.
l_max_overlap_in_cloud = .true.
l_lh_instant_var_covar_src = .true.
l_lh_limit_weights = .true.
l_lh_var_frac = .false.
l_lh_normalize_weights = .true.
sol_facti_cloud_borne = 1.0D0
dust_emis_fact = 0.3D0
nucleate_ice_subgrid = 1.0
seasalt_emis_scale = 0.6

! CLUBB history!!!
clubb_history = .true.
clubb_rad_history = .false.
clubb_vars_zt = $clubb_vars_zt_list
clubb_vars_zm = $clubb_vars_zm_list
EOF

# Run submission
echo "----- Run Submission -----"
./case.submit --no-batch 2>&1 | tee -a log

if [[ $(grep -c "ERROR" log) -eq 0 ]]; then
  rm log
  echo "Success!"
  exit 0
else
  rm log
  echo "Error submitting run" >> /dev/stderr
  exit 1
fi
