module micro_mg_utils

!--------------------------------------------------------------------------
!
! This module contains process rates and utility functions used by the MG
! microphysics.
!
! Original MG authors: Andrew Gettelman, Hugh Morrison
! Contributions from: Peter Caldwell, Xiaohong Liu and Steve Ghan
!
! Separated from MG 1.5 by B. Eaton.
! Separated module switched to MG 2.0 and further changes by S. Santos.
!
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!
!--------------------------------------------------------------------------
!
! List of required external functions that must be supplied:
!   gamma --> standard mathematical gamma function (if gamma is an
!       intrinsic, define HAVE_GAMMA_INTRINSICS)
!
!--------------------------------------------------------------------------
!
! Constants that must be specified in the "init" method (module variables):
!
! kind            kind of reals (to verify correct linkage only) -
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                   J kg-1 K-1
! rh2o            gas constant for water vapor                   J kg-1 K-1
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! tmelt           temperature of melting point for water         K
! latvap          latent heat of vaporization                    J kg-1
! latice          latent heat of fusion                          J kg-1
!
!--------------------------------------------------------------------------

! shr_spfn_gamma is not currently compatible with OpenACC. So compilation with 
! OpenACC requires intrinsic gamma.
#ifndef _OPENACC
#ifndef HAVE_GAMMA_INTRINSICS
use shr_spfn_mod, only: gamma => shr_spfn_gamma
#endif
#endif

implicit none
private
save

public :: &
     micro_mg_utils_init, &
     size_dist_param_liq, &
     size_dist_param_liq_2D, &
     size_dist_param_basic, &
     size_dist_param_basic_2D, &
     avg_diameter, &
     avg_diameter_2D, &
     rising_factorial, &
     ice_deposition_sublimation, &
     ice_deposition_sublimation_2D, &
     sb2001v2_liq_autoconversion,&
     sb2001v2_accre_cld_water_rain,&       
     sb2001v2_accre_cld_water_rain_2D, &
     kk2000_liq_autoconversion, &
     ice_autoconversion, &
     ice_autoconversion_2D, &
     immersion_freezing, &
     immersion_freezing_2D, &
     contact_freezing, &
     contact_freezing_2D, &
     snow_self_aggregation, &
     snow_self_aggregation_2D, &
     accrete_cloud_water_snow, &
     accrete_cloud_water_snow_2D, &
     secondary_ice_production, &
     secondary_ice_production_2D, &
     accrete_rain_snow, &
     accrete_rain_snow_2D, &
     heterogeneous_rain_freezing, &
     heterogeneous_rain_freezing_2D, &
     accrete_cloud_water_rain, &
     accrete_cloud_water_rain_2D, &
     self_collection_rain, &
     self_collection_rain_2D, &
     accrete_cloud_ice_snow, &
     accrete_cloud_ice_snow_2D, &
     evaporate_sublimate_precip, &
     evaporate_sublimate_precip_2D, &
     bergeron_process_snow, &
     bergeron_process_snow_2D, &
     graupel_collecting_snow, &
     graupel_collecting_snow_2D, &
     graupel_collecting_rain, &
     graupel_collecting_rain_2D, &
     graupel_collecting_cld_water, &
     graupel_collecting_cld_water_2D, &
     graupel_riming_liquid_snow, &
     graupel_riming_liquid_snow_2D, &
     graupel_rain_riming_snow, &
     graupel_rain_riming_snow_2D, &
     graupel_rime_splintering, &
     graupel_rime_splintering_2D, &
     evaporate_sublimate_precip_graupel, &
     evaporate_sublimate_precip_graupel_2D

! 8 byte real and integer
integer, parameter, public :: r8 = selected_real_kind(12)
integer, parameter, public :: i8 = selected_int_kind(18)

public :: MGHydrometeorProps

type :: MGHydrometeorProps
   ! Density (kg/m^3)
   real(r8) :: rho
   ! Information for size calculations.
   ! Basic calculation of mean size is:
   !     lambda = (shape_coef*nic/qic)^(1/eff_dim)
   ! Then lambda is constrained by bounds.
   real(r8) :: eff_dim
   real(r8) :: shape_coef
   real(r8) :: lambda_bounds(2)
   ! Minimum average particle mass (kg).
   ! Limit is applied at the beginning of the size distribution calculations.
   real(r8) :: min_mean_mass
end type MGHydrometeorProps

interface MGHydrometeorProps
   module procedure NewMGHydrometeorProps
end interface

type(MGHydrometeorProps), public :: mg_liq_props
type(MGHydrometeorProps), public :: mg_ice_props
type(MGHydrometeorProps), public :: mg_rain_props
type(MGHydrometeorProps), public :: mg_snow_props
type(MGHydrometeorProps), public :: mg_graupel_props
type(MGHydrometeorProps), public :: mg_hail_props

interface size_dist_param_liq
  module procedure size_dist_param_liq_vect
  module procedure size_dist_param_liq_line
end interface
interface size_dist_param_basic
  module procedure size_dist_param_basic_vect
  module procedure size_dist_param_basic_line
end interface

!=================================================
! Public module parameters (mostly for MG itself)
!=================================================

! Pi to 20 digits; more than enough to reach the limit of double precision.
real(r8), parameter, public :: pi = 3.14159265358979323846_r8

! "One minus small number": number near unity for round-off issues.
real(r8), parameter, public :: omsm   = 1._r8 - 1.e-5_r8

! Smallest mixing ratio considered in microphysics.
real(r8), parameter, public :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter, public :: mincld = 0.0001_r8

real(r8), parameter, public :: rhosn = 250._r8  ! bulk density snow
real(r8), parameter, public :: rhoi = 500._r8   ! bulk density ice
real(r8), parameter, public :: rhow = 1000._r8  ! bulk density liquid
real(r8), parameter, public :: rhows = 917._r8  ! bulk density water solid

!Hail and Graupel (set in MG3)
real(r8), parameter, public :: rhog = 500._r8 
real(r8), parameter, public :: rhoh = 900._r8 

! fall speed parameters, V = aD^b (V is in m/s)
! droplets
real(r8), parameter, public :: ac = 3.e7_r8
real(r8), parameter, public :: bc = 2._r8
! snow
real(r8), parameter, public :: as = 11.72_r8
real(r8), parameter, public :: bs = 0.41_r8
! cloud ice
real(r8), parameter, public :: ai = 700._r8
real(r8), parameter, public :: bi = 1._r8
! small cloud ice (r< 10 um) - sphere, bulk density
real(r8), parameter, public :: aj = ac*((rhoi/rhows)**(bc/3._r8))*rhows/rhow
real(r8), parameter, public :: bj = bc
! rain
real(r8), parameter, public :: ar = 841.99667_r8
real(r8), parameter, public :: br = 0.8_r8
! graupel
real(r8), parameter, public :: ag = 19.3_r8
real(r8), parameter, public :: bg = 0.37_r8
! hail
real(r8), parameter, public :: ah = 114.5_r8 
real(r8), parameter, public :: bh = 0.5_r8

! mass of new crystal due to aerosol freezing and growth (kg)
! Make this consistent with the lower bound, to support UTLS and
! stratospheric ice, and the smaller ice size limit.
real(r8), parameter, public :: mi0 = 4._r8/3._r8*pi*rhoi*(1.e-6_r8)**3

! mass of new graupel particle  (assume same as mi0 for now, may want to make bigger?)
!real(r8), parameter, public :: mg0 = 4._r8/3._r8*pi*rhoi*(1.e-6_r8)**3
!or set based on M2005:
real(r8), parameter, public :: mg0 = 1.6e-10_r8
! radius of contact nuclei
real(r8), parameter, public :: mmult = 4._r8/3._r8*pi*rhoi*(5.e-6_r8)**3

!=================================================
! Private module parameters
!=================================================

! Signaling NaN bit pattern that represents a limiter that's turned off.
integer(i8), parameter :: limiter_off = int(Z'7FF1111111111111', i8)

! alternate threshold used for some in-cloud mmr
real(r8), parameter :: icsmall = 1.e-8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d
! exponent
real(r8), parameter :: dsph = 3._r8

! Bounds for mean diameter for different constituents.
real(r8), parameter :: lam_bnd_rain(2) = 1._r8/[500.e-6_r8, 20.e-6_r8]
real(r8), parameter :: lam_bnd_snow(2) = 1._r8/[2000.e-6_r8, 10.e-6_r8]

! Minimum average mass of particles.
real(r8), parameter :: min_mean_mass_liq = 1.e-20_r8
real(r8), parameter :: min_mean_mass_ice = 1.e-20_r8

! ventilation parameters
! for snow
real(r8), parameter :: f1s = 0.86_r8
real(r8), parameter :: f2s = 0.28_r8
! for rain
real(r8), parameter :: f1r = 0.78_r8
real(r8), parameter :: f2r = 0.308_r8

! collection efficiencies
! aggregation of cloud ice and snow
real(r8), parameter :: eii = 0.5_r8
! collection efficiency, ice-droplet collisions
real(r8), parameter, public :: ecid = 0.7_r8
! collection efficiency between droplets/rain and snow/rain
real(r8), parameter, public :: ecr = 1.0_r8

! immersion freezing parameters, bigg 1953
real(r8), parameter :: bimm = 100._r8
real(r8), parameter :: aimm = 0.66_r8

! Mass of each raindrop created from autoconversion.
real(r8), parameter :: droplet_mass_25um = 4._r8/3._r8*pi*rhow*(25.e-6_r8)**3
real(r8), parameter :: droplet_mass_40um = 4._r8/3._r8*pi*rhow*(40.e-6_r8)**3

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)

real(r8) :: ra        ! dry air gas constant

! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

! additional constants to help speed up code
real(r8) :: gamma_bs_plus3
real(r8) :: gamma_half_br_plus5
real(r8) :: gamma_half_bs_plus5
real(r8) :: gamma_2bs_plus2

!=========================================================
! Utilities that are cheaper if the compiler knows that
! some argument is an integer.
!=========================================================

interface rising_factorial
   module procedure rising_factorial_r8
   module procedure rising_factorial_2D_r8
   module procedure rising_factorial_integer
   module procedure rising_factorial_2D_integer
end interface rising_factorial

interface var_coef
   module procedure var_coef_r8
   module procedure var_coef_integer
end interface var_coef

!==========================================================================
contains
!==========================================================================

! Initialize module variables.
!
! "kind" serves no purpose here except to check for unlikely linking
! issues; always pass in the kind for a double precision real.
!
! "errstring" is the only output; it is blank if there is no error, or set
! to a message if there is an error.
!
! Check the list at the top of this module for descriptions of all other
! arguments.
subroutine micro_mg_utils_init( kind, rair, rh2o, cpair, tmelt_in, latvap, &
     latice, dcs, errstring)

  integer,  intent(in)  :: kind
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice
  real(r8), intent(in)  :: dcs

  character(128), intent(out) :: errstring

  ! Name this array to workaround an XLF bug (otherwise could just use the
  ! expression that sets it).
  real(r8) :: ice_lambda_bounds(2)

  !-----------------------------------------------------------------------

  errstring = ' '

  if( kind .ne. r8 ) then
     errstring = 'micro_mg_init: KIND of reals does not match'
     return
  endif

  ! declarations for MG code (transforms variable names)

  rv= rh2o                  ! water vapor gas constant
  cpp = cpair               ! specific heat of dry air
  tmelt = tmelt_in

  ra = rair     !dry air gas constant

  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! Define constants to help speed up code (this limits calls to gamma function)
  gamma_bs_plus3=gamma(3._r8+bs)
  gamma_half_br_plus5=gamma(5._r8/2._r8+br/2._r8)
  gamma_half_bs_plus5=gamma(5._r8/2._r8+bs/2._r8)
  gamma_2bs_plus2=gamma(2._r8*bs+2._r8)  

  ! Don't specify lambda bounds for cloud liquid, as they are determined by
  ! pgam dynamically.
  mg_liq_props = MGHydrometeorProps(rhow, dsph, &
       min_mean_mass=min_mean_mass_liq)

  ! Mean ice diameter can not grow bigger than twice the autoconversion
  ! threshold for snow.
  ice_lambda_bounds = 1._r8/[2._r8*dcs, 1.e-6_r8]

  mg_ice_props = MGHydrometeorProps(rhoi, dsph, &
       ice_lambda_bounds, min_mean_mass_ice)

  mg_rain_props = MGHydrometeorProps(rhow, dsph, lam_bnd_rain)
  mg_snow_props = MGHydrometeorProps(rhosn, dsph, lam_bnd_snow)
  mg_graupel_props = MGHydrometeorProps(rhog, dsph, lam_bnd_snow)
  mg_hail_props = MGHydrometeorProps(rhoh, dsph, lam_bnd_snow)


end subroutine micro_mg_utils_init

! Constructor for a constituent property object.
function NewMGHydrometeorProps(rho, eff_dim, lambda_bounds, min_mean_mass) &
     result(res)
  real(r8), intent(in) :: rho, eff_dim
  real(r8), intent(in), optional :: lambda_bounds(2), min_mean_mass
  type(MGHydrometeorProps) :: res

  res%rho = rho
  res%eff_dim = eff_dim
  if (present(lambda_bounds)) then
     res%lambda_bounds = lambda_bounds
  else
     res%lambda_bounds = no_limiter()
  end if
  if (present(min_mean_mass)) then
     res%min_mean_mass = min_mean_mass
  else
     res%min_mean_mass = no_limiter()
  end if

  res%shape_coef = rho*pi*gamma(eff_dim+1._r8)/6._r8

end function NewMGHydrometeorProps

!========================================================================
!FORMULAS
!========================================================================

! Use gamma function to implement rising factorial extended to the reals.
pure function rising_factorial_r8(x, n) result(res)
  real(r8), intent(in) :: x, n
  real(r8) :: res

  res = gamma(x+n)/gamma(x)

end function rising_factorial_r8

! General rising factorial for 2 dimensions 
pure function rising_factorial_2D_r8(mgncol, nlev, x, n) result(res)
  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: x
  real(r8), intent(in) :: n
  real(r8), dimension(mgncol,nlev) :: res
  
  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      res(i,k) = gamma(x(i,k)+n)/gamma(x(i,k))
    end do
  end do

end function rising_factorial_2D_r8

! Rising factorial can be performed much cheaper if n is a small integer.
pure function rising_factorial_integer(x, n) result(res)
  real(r8), intent(in) :: x
  integer, intent(in) :: n
  real(r8) :: res

  integer :: i
  real(r8) :: factor

  res = 1._r8
  factor = x

  do i = 1, n
     res = res * factor
     factor = factor + 1._r8
  end do

end function rising_factorial_integer

! Rising factorial can be performed much cheaper if n is a small integer.
pure function rising_factorial_2D_integer(mgncol, nlev, x, n) result(res)
  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: x
  integer, intent(in) :: n
  real(r8), dimension(mgncol,nlev) :: res

  integer :: i, k, m
  real(r8) :: factor

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      res(i,k) = 1._r8
      factor = x(i,k)

      do m = 1, n
         res(i,k) = res(i,k) * factor
         factor = factor + 1._r8
      end do
    end do
  end do

end function rising_factorial_2D_integer

! Calculate correction due to latent heat for evaporation/sublimation
elemental function calc_ab(t, qv, xxl) result(ab)
  
  real(r8), intent(in) :: t     ! Temperature
  real(r8), intent(in) :: qv    ! Saturation vapor pressure
  real(r8), intent(in) :: xxl   ! Latent heat

  real(r8) :: ab

  real(r8) :: dqsdt

  dqsdt = xxl*qv / (rv * t**2)
  ab = 1._r8 + dqsdt*xxl/cpp

end function calc_ab

! Calculate correction due to latent heat for evaporation/sublimation
function calc_ab_2D(mgncol, nlev, t, qv, xxl) result(ab)
  
  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: t     ! Temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: qv    ! Saturation vapor pressure
  real(r8), intent(in) :: xxl   ! Latent heat

  real(r8), dimension(mgncol,nlev) :: ab

  real(r8) :: dqsdt
  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
      dqsdt = xxl*qv(i,k) / (rv * t(i,k)**2)
      ab(i,k) = 1._r8 + dqsdt*xxl/cpp
    end do
  end do

end function calc_ab_2D

! get cloud droplet size distribution parameters
elemental subroutine size_dist_param_liq_line(props, qcic, ncic, rho, pgam, lamc)
  type(MGHydrometeorProps), intent(in) :: props
  real(r8), intent(in) :: qcic
  real(r8), intent(inout) :: ncic
  real(r8), intent(in) :: rho

  real(r8), intent(out) :: pgam
  real(r8), intent(out) :: lamc

  type(MGHydrometeorProps) :: props_loc

  if (qcic > qsmall) then

     ! Local copy of properties that can be modified.
     ! (Elemental routines that operate on arrays can't modify scalar
     ! arguments.)
     props_loc = props

     ! Get pgam from fit Rotstayn and Liu 2003 (changed from Martin 1994 for CAM6)
     pgam = 1.0_r8 - 0.7_r8 * exp(-0.008_r8*1.e-6_r8*ncic*rho)
     pgam = 1._r8/(pgam**2) - 1._r8
     pgam = max(pgam, 2._r8)

     ! Set coefficient for use in size_dist_param_basic.
     ! The 3D case is so common and optimizable that we specialize it:
     if (props_loc%eff_dim == 3._r8) then
        props_loc%shape_coef = pi / 6._r8 * props_loc%rho * &
             rising_factorial(pgam+1._r8, 3)
     else
        props_loc%shape_coef = pi / 6._r8 * props_loc%rho * &
             rising_factorial(pgam+1._r8, props_loc%eff_dim)
     end if

     ! Limit to between 2 and 50 microns mean size.
     props_loc%lambda_bounds = (pgam+1._r8)*1._r8/[50.e-6_r8, 2.e-6_r8]

     call size_dist_param_basic(props_loc, qcic, ncic, lamc)

  else
     ! pgam not calculated in this case, so set it to a value likely to
     ! cause an error if it is accidentally used
     ! (gamma function undefined for negative integers)
     pgam = -100._r8
     lamc = 0._r8
  end if

end subroutine size_dist_param_liq_line

! get cloud droplet size distribution parameters

subroutine size_dist_param_liq_vect(props, qcic, ncic, rho, pgam, lamc, mgncol)

  type(mghydrometeorprops), intent(in) :: props
  integer,                          intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: qcic
  real(r8), dimension(mgncol), intent(inout) :: ncic
  real(r8), dimension(mgncol), intent(in) :: rho
  real(r8), dimension(mgncol), intent(out) :: pgam
  real(r8), dimension(mgncol), intent(out) :: lamc
  type(mghydrometeorprops) :: props_loc
  integer :: i

  do i=1,mgncol
     if (qcic(i) > qsmall) then
        ! Local copy of properties that can be modified.
        ! (Elemental routines that operate on arrays can't modify scalar
        ! arguments.)
        props_loc = props
        ! Get pgam from fit by Rotstayn and Liu 2003 (changed from Martin 1994 for CAM6)
        pgam(i) = 1.0_r8 - 0.7_r8 * exp(-0.008_r8*1.e-6_r8*ncic(i)*rho(i))
        pgam(i) = 1._r8/(pgam(i)**2) - 1._r8
        pgam(i) = max(pgam(i), 2._r8)
     endif
  enddo
  do i=1,mgncol
     if (qcic(i) > qsmall) then
        ! Set coefficient for use in size_dist_param_basic.
        ! The 3D case is so common and optimizable that we specialize
        ! it:
        if (props_loc%eff_dim == 3._r8) then
           props_loc%shape_coef = pi / 6._r8 * props_loc%rho * &
                rising_factorial(pgam(i)+1._r8, 3)
        else
           props_loc%shape_coef = pi / 6._r8 * props_loc%rho * &
                rising_factorial(pgam(i)+1._r8, props_loc%eff_dim)
        end if
        ! Limit to between 2 and 50 microns mean size.
        props_loc%lambda_bounds(1) = (pgam(i)+1._r8)*1._r8/50.e-6_r8
        props_loc%lambda_bounds(2) = (pgam(i)+1._r8)*1._r8/2.e-6_r8
        call size_dist_param_basic(props_loc, qcic(i), ncic(i), lamc(i))
     endif
  enddo
  do i=1,mgncol
     if (qcic(i) <= qsmall) then
        ! pgam not calculated in this case, so set it to a value likely to
        ! cause an error if it is accidentally used
        ! (gamma function undefined for negative integers)
        pgam(i) = -100._r8
        lamc(i) = 0._r8
     end if
  enddo

end subroutine size_dist_param_liq_vect

subroutine size_dist_param_liq_2D(mgncol, nlev, props, qcic, ncic, rho, &
                                  pgam, lamc)

  integer, intent(in) :: mgncol
  integer, intent(in) :: nlev
  
  type(mghydrometeorprops), intent(in) :: props
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic
  real(r8), dimension(mgncol,nlev), intent(inout) :: ncic
  real(r8), dimension(mgncol,nlev), intent(in) :: rho
  
  real(r8), dimension(mgncol,nlev), intent(out) :: pgam
  real(r8), dimension(mgncol,nlev), intent(out) :: lamc
  
  real(r8) :: eff_dim
              
  real(r8), dimension(mgncol,nlev) :: lower_bound, &
                                      upper_bound, &
                                      shape_coef, &
                                      pgam_tmp, &
                                      rising_factorial_term
  
  integer :: i, k
  
  !$acc data create(lower_bound, upper_bound, shape_coef, pgam_tmp, rising_factorial_term) &
  !$acc&     copyin(props) async(1)
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qcic(i,k) > qsmall) then
        ! Get pgam from fit by Rotstayn and Liu 2003 (changed from Martin 1994 for CAM6)
        pgam(i,k) = 1.0_r8 - 0.7_r8 * exp(-0.008_r8*1.e-6_r8*ncic(i,k)*rho(i,k))
        pgam(i,k) = 1._r8/(pgam(i,k)**2) - 1._r8
        pgam(i,k) = max(pgam(i,k), 2._r8)
        pgam_tmp(i,k) = pgam(i,k) + 1.0_r8
      else
        ! pgam not calculated in this case, so set it to a value likely to
        ! cause an error if it is accidentally used
        ! (gamma function undefined for negative integers)
        pgam(i,k) = -100._r8
        pgam_tmp(i,k) = 0.0_r8
      endif
    end do
  end do
  
  ! Set coefficient for use in size_dist_param_basic.
  ! The 3D case is so common and optimizable that we specialize it
  if (props%eff_dim == 3._r8) then
    rising_factorial_term(:,:) = rising_factorial(mgncol, nlev, pgam_tmp, 3)
  else
    rising_factorial_term(:,:) = rising_factorial(mgncol, nlev, pgam_tmp, props%eff_dim)
  end if
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qcic(i,k) > qsmall) then
        shape_coef(i,k) = pi / 6._r8 * props%rho * rising_factorial_term(i,k)
      end if
    end do
  end do
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      ! Limit to between 2 and 50 microns mean size.
      lower_bound(i,k) = pgam_tmp(i,k) * 1._r8 / 50.e-6_r8
      upper_bound(i,k) = pgam_tmp(i,k) * 1._r8 / 2.e-6_r8
    end do
  end do
        
  ! add upper limit to in-cloud number concentration to prevent
  ! numerical error
  if (limiter_is_on(props%min_mean_mass)) then
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k = 1, nlev
      do i = 1, mgncol
       ncic(i,k) = min(ncic(i,k), qcic(i,k) / props%min_mean_mass)
     end do
    end do
  end if

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qcic(i,k) > qsmall) then
        ! lambda = (c n/q)^(1/d)
        lamc(i,k) = (shape_coef(i,k) * ncic(i,k)/qcic(i,k))**(1._r8/props%eff_dim)
      else
        lamc(i,k) = 0._r8
      endif
    end do
  end do

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qcic(i,k) > qsmall) then
        ! check for slope
        ! adjust vars
        if ( lamc(i,k) < lower_bound(i,k) ) then
           lamc(i,k) = lower_bound(i,k)
           ncic(i,k) = lower_bound(i,k)**props%eff_dim * qcic(i,k)/shape_coef(i,k)
        else if ( lamc(i,k) > upper_bound(i,k) ) then
           lamc(i,k) = upper_bound(i,k)
           ncic(i,k) = upper_bound(i,k)**props%eff_dim * qcic(i,k)/shape_coef(i,k)
        end if
      end if
    end do
  end do
  
  !$acc end data

end subroutine size_dist_param_liq_2D

! Basic routine for getting size distribution parameters.
elemental subroutine size_dist_param_basic_line(props, qic, nic, lam, n0)
  type(MGHydrometeorProps), intent(in) :: props
  real(r8), intent(in) :: qic
  real(r8), intent(inout) :: nic

  real(r8), intent(out) :: lam
  real(r8), intent(out), optional :: n0

  if (qic > qsmall) then

     ! add upper limit to in-cloud number concentration to prevent
     ! numerical error
     if (limiter_is_on(props%min_mean_mass)) then
        nic = min(nic, qic / props%min_mean_mass)
     end if

     ! lambda = (c n/q)^(1/d)
     lam = (props%shape_coef * nic/qic)**(1._r8/props%eff_dim)

     ! check for slope
     ! adjust vars
     if (lam < props%lambda_bounds(1)) then
        lam = props%lambda_bounds(1)
        nic = lam**(props%eff_dim) * qic/props%shape_coef
     else if (lam > props%lambda_bounds(2)) then
        lam = props%lambda_bounds(2)
        nic = lam**(props%eff_dim) * qic/props%shape_coef
     end if

  else
     lam = 0._r8
  end if

  if (present(n0)) n0 = nic * lam

end subroutine size_dist_param_basic_line

subroutine size_dist_param_basic_vect(props, qic, nic, lam, mgncol, n0)

  type (mghydrometeorprops), intent(in) :: props
  integer,                          intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: qic
  real(r8), dimension(mgncol), intent(inout) :: nic
  real(r8), dimension(mgncol), intent(out) :: lam
  real(r8), dimension(mgncol), intent(out), optional :: n0
  integer :: i
  do i=1,mgncol

     if (qic(i) > qsmall) then

        ! add upper limit to in-cloud number concentration to prevent
        ! numerical error
        if (limiter_is_on(props%min_mean_mass)) then
           nic(i) = min(nic(i), qic(i) / props%min_mean_mass)
        end if

        ! lambda = (c n/q)^(1/d)
        lam(i) = (props%shape_coef * nic(i)/qic(i))**(1._r8/props%eff_dim)

        ! check for slope
        ! adjust vars
        if (lam(i) < props%lambda_bounds(1)) then
           lam(i) = props%lambda_bounds(1)
           nic(i) = lam(i)**(props%eff_dim) * qic(i)/props%shape_coef
        else if (lam(i) > props%lambda_bounds(2)) then
           lam(i) = props%lambda_bounds(2)
           nic(i) = lam(i)**(props%eff_dim) * qic(i)/props%shape_coef
        end if

     else
        lam(i) = 0._r8
     end if

  enddo

  if (present(n0)) n0 = nic * lam

end subroutine size_dist_param_basic_vect

subroutine size_dist_param_basic_2D(mgncol, nlev, props, qic, nic, lam, n0)

  integer,                   intent(in) :: mgncol, nlev
  type (mghydrometeorprops), intent(in) :: props
  
  real(r8), dimension(mgncol,nlev), intent(in) :: qic
  
  real(r8), dimension(mgncol,nlev), intent(inout) :: nic
  
  real(r8), dimension(mgncol,nlev), intent(out) :: lam
  real(r8), dimension(mgncol,nlev), intent(out), optional :: n0
  
  integer :: i, k
  
  !$acc data copyin(props) async(1)
  
  ! add upper limit to in-cloud number concentration to prevent
  ! numerical error
  if (limiter_is_on(props%min_mean_mass)) then
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k = 1, nlev
      do i = 1, mgncol
        if (qic(i,k) > qsmall) then
          nic(i,k) = min(nic(i,k), qic(i,k) / props%min_mean_mass)
        end if
      end do
    end do
  end if
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol

      if (qic(i,k) > qsmall) then

        ! lambda = (c n/q)^(1/d)
        lam(i,k) = (props%shape_coef * nic(i,k)/qic(i,k))**(1._r8/props%eff_dim)

        ! check for slope
        ! adjust vars
        if (lam(i,k) < props%lambda_bounds(1)) then
           lam(i,k) = props%lambda_bounds(1)
           nic(i,k) = lam(i,k)**(props%eff_dim) * qic(i,k)/props%shape_coef
        else if (lam(i,k) > props%lambda_bounds(2)) then
           lam(i,k) = props%lambda_bounds(2)
           nic(i,k) = lam(i,k)**(props%eff_dim) * qic(i,k)/props%shape_coef
        end if

      else
        lam(i,k) = 0._r8
      end if

    end do
  end do
  
  !$acc end data
  
  if (present(n0)) then
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k = 1, nlev
      do i = 1, mgncol
        n0(i,k) = nic(i,k) * lam(i,k)
      end do
    end do
  end if

end subroutine size_dist_param_basic_2D


real(r8) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  real(r8), intent(in) :: q         ! mass mixing ratio
  real(r8), intent(in) :: n         ! number concentration (per volume)
  real(r8), intent(in) :: rho_air   ! local density of the air
  real(r8), intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (pi * rho_sub * n/(q*rho_air))**(-1._r8/3._r8)

end function avg_diameter

function avg_diameter_2D(mgncol, nlev, q, n, rho_air, rho_sub) result(res)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  integer, intent(in) :: mgncol, nlev
  real(r8), intent(in), dimension(mgncol,nlev) :: q         ! mass mixing ratio
  real(r8), intent(in), dimension(mgncol,nlev) :: n         ! number concentration (per volume)
  real(r8), intent(in), dimension(mgncol,nlev) :: rho_air   ! local density of the air
  real(r8), intent(in) :: rho_sub   ! density of the particle substance
  real(r8), dimension(mgncol,nlev) :: res   ! local density of the air
  

  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
      if ( (q(i,k) * rho_air(i,k)) > 0._r8 ) then
        res(i,k) = (pi * rho_sub * n(i,k)/(q(i,k)*rho_air(i,k)))**(-1._r8/3._r8)
      else
        res(i,k) = 0._r8
      end if
    end do
  end do

end function avg_diameter_2D

elemental function var_coef_r8(relvar, a) result(res)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(r8), intent(in) :: relvar
  real(r8), intent(in) :: a
  real(r8) :: res

  res = rising_factorial(relvar, a) / relvar**a

end function var_coef_r8

elemental function var_coef_integer(relvar, a) result(res)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(r8), intent(in) :: relvar
  integer, intent(in) :: a
  real(r8) :: res

  res = rising_factorial(relvar, a) / relvar**a

end function var_coef_integer

!========================================================================
!MICROPHYSICAL PROCESS CALCULATIONS
!========================================================================
!========================================================================
! Initial ice deposition and sublimation loop.
! Run before the main loop
! This subroutine written by Peter Caldwell

subroutine ice_deposition_sublimation(t, qv, qi, ni, &
                                      icldm, rho, dv,qvl, qvi, &
                                      berg, vap_dep, ice_sublim, mgncol)

  !INPUT VARS:
  !===============================================
  integer,  intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: t
  real(r8), dimension(mgncol), intent(in) :: qv
  real(r8), dimension(mgncol), intent(in) :: qi
  real(r8), dimension(mgncol), intent(in) :: ni
  real(r8), dimension(mgncol), intent(in) :: icldm
  real(r8), dimension(mgncol), intent(in) :: rho
  real(r8), dimension(mgncol), intent(in) :: dv
  real(r8), dimension(mgncol), intent(in) :: qvl
  real(r8), dimension(mgncol), intent(in) :: qvi

  !OUTPUT VARS:
  !===============================================
  real(r8), dimension(mgncol), intent(out) :: vap_dep !ice deposition (cell-ave value)
  real(r8), dimension(mgncol), intent(out) :: ice_sublim !ice sublimation (cell-ave value)
  real(r8), dimension(mgncol), intent(out) :: berg !bergeron enhancement (cell-ave value)

  !INTERNAL VARS:
  !===============================================
  real(r8) :: ab
  real(r8) :: epsi
  real(r8) :: qiic
  real(r8) :: niic
  real(r8) :: lami
  real(r8) :: n0i
  integer :: i

  do i=1,mgncol
     if (qi(i)>=qsmall) then

        !GET IN-CLOUD qi, ni
        !===============================================
        qiic = qi(i)/icldm(i)
        niic = ni(i)/icldm(i)

        !Compute linearized condensational heating correction
        ab=calc_ab(t(i), qvi(i), xxls)
        !Get slope and intercept of gamma distn for ice.
        call size_dist_param_basic(mg_ice_props, qiic, niic, lami, n0i)
        !Get depletion timescale=1/eps
        epsi = 2._r8*pi*n0i*rho(i)*Dv(i)/(lami*lami)

        !Compute deposition/sublimation
        vap_dep(i) = epsi/ab*(qv(i) - qvi(i))

        !Make this a grid-averaged quantity
        vap_dep(i)=vap_dep(i)*icldm(i)

        !Split into deposition or sublimation.
        if (t(i) < tmelt .and. vap_dep(i)>0._r8) then
           ice_sublim(i)=0._r8
        else
        ! make ice_sublim negative for consistency with other evap/sub processes
           ice_sublim(i)=min(vap_dep(i),0._r8)
           vap_dep(i)=0._r8
        end if

        !sublimation occurs @ any T. Not so for berg.
        if (t(i) < tmelt) then

           !Compute bergeron rate assuming cloud for whole step.
           berg(i) = max(epsi/ab*(qvl(i) - qvi(i)), 0._r8)
        else !T>frz
           berg(i)=0._r8
        end if !T<frz

     else !where qi<qsmall
        berg(i)=0._r8
        vap_dep(i)=0._r8
        ice_sublim(i)=0._r8
     end if !qi>qsmall
  enddo
end subroutine ice_deposition_sublimation

subroutine ice_deposition_sublimation_2D(mgncol, nlev, t, qv, qi, ni, &
                                         icldm, rho, dv,qvl, qvi, &
                                         berg, vap_dep, ice_sublim )

  !INPUT VARS:
  !===============================================
  integer,  intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: t
  real(r8), dimension(mgncol,nlev), intent(in) :: qv
  real(r8), dimension(mgncol,nlev), intent(in) :: qi
  real(r8), dimension(mgncol,nlev), intent(in) :: ni
  real(r8), dimension(mgncol,nlev), intent(in) :: icldm
  real(r8), dimension(mgncol,nlev), intent(in) :: rho
  real(r8), dimension(mgncol,nlev), intent(in) :: dv
  real(r8), dimension(mgncol,nlev), intent(in) :: qvl
  real(r8), dimension(mgncol,nlev), intent(in) :: qvi

  !OUTPUT VARS:
  !===============================================
  real(r8), dimension(mgncol,nlev), intent(out) :: vap_dep !ice deposition (cell-ave value)
  real(r8), dimension(mgncol,nlev), intent(out) :: ice_sublim !ice sublimation (cell-ave value)
  real(r8), dimension(mgncol,nlev), intent(out) :: berg !bergeron enhancement (cell-ave value)

  !INTERNAL VARS:
  !===============================================
  real(r8), dimension(mgncol,nlev) :: ab
  real(r8) :: epsi
  real(r8), dimension(mgncol,nlev) :: qiic
  real(r8), dimension(mgncol,nlev) :: niic
  real(r8), dimension(mgncol,nlev) :: lami
  real(r8), dimension(mgncol,nlev) :: n0i
  integer :: i, k
  
  !$acc data create(ab, qiic, niic, lami, n0i) async(1)
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      !GET IN-CLOUD qi, ni
      !===============================================
      qiic(i,k) = qi(i,k) / icldm(i,k)
      niic(i,k) = ni(i,k) / icldm(i,k)
    end do
  end do
  
  !Compute linearized condensational heating correction
  ab(:,:) = calc_ab_2D(mgncol, nlev, t(:,:), qvi(:,:), xxls)
        
  !Get slope and intercept of gamma distn for ice.
  call size_dist_param_basic_2D(mgncol, nlev, mg_ice_props, qiic(:,:), niic(:,:), lami(:,:), n0i(:,:))

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qi(i,k)>=qsmall) then

        !Get depletion timescale=1/eps
        epsi = 2._r8*pi*n0i(i,k)*rho(i,k)*Dv(i,k)/(lami(i,k)*lami(i,k))

        !Compute deposition/sublimation
        vap_dep(i,k) = epsi/ab(i,k)*(qv(i,k) - qvi(i,k))

        !Make this a grid-averaged quantity
        vap_dep(i,k)=vap_dep(i,k)*icldm(i,k)

        !Split into deposition or sublimation.
        if (t(i,k) < tmelt .and. vap_dep(i,k)>0._r8) then
           ice_sublim(i,k)=0._r8
        else
        ! make ice_sublim negative for consistency with other evap/sub processes
           ice_sublim(i,k)=min(vap_dep(i,k),0._r8)
           vap_dep(i,k)=0._r8
        end if

        !sublimation occurs @ any T. Not so for berg.
        if (t(i,k) < tmelt) then

           !Compute bergeron rate assuming cloud for whole step.
           berg(i,k) = max(epsi/ab(i,k)*(qvl(i,k) - qvi(i,k)), 0._r8)
        else !T>frz
           berg(i,k)=0._r8
        end if !T<frz

      else !where qi<qsmall
        berg(i,k)=0._r8
        vap_dep(i,k)=0._r8
        ice_sublim(i,k)=0._r8
      end if !qi>qsmall
    end do
  end do
  
  !$acc end data
  
end subroutine ice_deposition_sublimation_2D

!========================================================================
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
! minimum qc of 1 x 10^-8 prevents floating point error

subroutine kk2000_liq_autoconversion(mgncol, nlev, microp_uniform, qcic, &
                                     ncic, rho, relvar, &
                                     prc, nprc, nprc1 )

  integer, intent(in) :: mgncol, nlev
  logical, intent(in) :: microp_uniform

  real(r8), dimension(mgncol,nlev), intent(in) :: qcic
  real(r8), dimension(mgncol,nlev), intent(in) :: ncic
  real(r8), dimension(mgncol,nlev), intent(in) :: rho

  real(r8), dimension(mgncol,nlev), intent(in) :: relvar

  real(r8), dimension(mgncol,nlev), intent(out) :: prc
  real(r8), dimension(mgncol,nlev), intent(out) :: nprc
  real(r8), dimension(mgncol,nlev), intent(out) :: nprc1

  real(r8), dimension(mgncol,nlev) :: prc_coef
  
  integer :: i, k
  
  !$acc data create(prc_coef) async(1)
  
  ! Take variance into account, or use uniform value.
  if (.not. microp_uniform) then
    !$acc update host(relvar)
    do k = 1, nlev
      do i = 1, mgncol
        prc_coef(i,k) = var_coef(relvar(i,k), 2.47_r8)
      end do
    end do
    !$acc update device(prc_coef)
  else
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k = 1, nlev
      do i = 1, mgncol
        prc_coef(i,k) = 1._r8
      end do
    end do
  end if

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1,mgncol
      if (qcic(i,k) >= icsmall) then

        ! nprc is increase in rain number conc due to autoconversion
        ! nprc1 is decrease in cloud droplet conc due to autoconversion

        ! assume exponential sub-grid distribution of qc, resulting in additional
        ! factor related to qcvar below
        ! switch for sub-columns, don't include sub-grid qc

        prc(i,k) = prc_coef(i,k) * &
             0.01_r8 * 1350._r8 * qcic(i,k)**2.47_r8 * (ncic(i,k)*1.e-6_r8*rho(i,k))**(-1.1_r8)
        nprc(i,k) = prc(i,k) * (1._r8/droplet_mass_25um)
        nprc1(i,k) = prc(i,k)*ncic(i,k)/qcic(i,k)

      else
        prc(i,k)=0._r8
        nprc(i,k)=0._r8
        nprc1(i,k)=0._r8
      end if
    end do
  end do
  
  !$acc end data
  
end subroutine kk2000_liq_autoconversion
  
  !========================================================================
subroutine sb2001v2_liq_autoconversion(mgncol,nlev,pgam,qc,nc,qr,rho,relvar,au,nprc,nprc1)
  !
  ! ---------------------------------------------------------------------
  ! AUTO_SB:  calculates the evolution of mass- and number mxg-ratio for
  ! drizzle drops due to autoconversion. The autoconversion rate assumes
  ! f(x)=A*x**(nu_c)*exp(-Bx) in drop MASS x. 

  ! Code from Hugh Morrison, Sept 2014

  ! autoconversion
  ! use simple lookup table of dnu values to get mass spectral shape parameter
  ! equivalent to the size spectral shape parameter pgam
    
  integer, intent(in) :: mgncol, nlev
  
  real(r8), dimension(mgncol,nlev), intent (in)    :: pgam
  real(r8), dimension(mgncol,nlev), intent (in)    :: qc  ! = qc (cld water mixing ratio)
  real(r8), dimension(mgncol,nlev), intent (in)    :: nc  ! = nc (cld water number conc /kg)    
  real(r8), dimension(mgncol,nlev), intent (in)    :: qr  ! = qr (rain water mixing ratio)
  real(r8), dimension(mgncol,nlev), intent (in)    :: rho ! = rho : density profile
  real(r8), dimension(mgncol,nlev), intent (in)    :: relvar 
  
  real(r8), dimension(mgncol,nlev), intent (out)   :: au ! = prc autoconversion rate
  real(r8), dimension(mgncol,nlev), intent (out)   :: nprc1 ! = number tendency
  real(r8), dimension(mgncol,nlev), intent (out)   :: nprc ! = number tendency fixed size for rain
 
  ! parameters for droplet mass spectral shape, 
  !used by Seifert and Beheng (2001)                             
  ! warm rain scheme only (iparam = 1)                                                                        
  real(r8), parameter :: dnu(16) = [0._r8,-0.557_r8,-0.430_r8,-0.307_r8, & 
     -0.186_r8,-0.067_r8,0.050_r8,0.167_r8,0.282_r8,0.397_r8,0.512_r8, &
     0.626_r8,0.739_r8,0.853_r8,0.966_r8,0.966_r8]

  ! parameters for Seifert and Beheng (2001) autoconversion/accretion                                         
  real(r8), parameter :: kc = 9.44e9_r8
  real(r8), parameter :: kr = 5.78e3_r8
  real(r8) :: dum, dum1, nu
  integer :: dumi, i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol

      if (qc(i,k) > qsmall) then
        dumi=int(pgam(i,k))
        nu=dnu(dumi)+(dnu(dumi+1)-dnu(dumi))* &
               (pgam(i,k)-dumi)

        dum = 1._r8-qc(i,k)/(qc(i,k)+qr(i,k))
        dum1 = 600._r8*dum**0.68_r8*(1._r8-dum**0.68_r8)**3

        au(i,k) = kc/(20._r8*2.6e-7_r8)* &
          (nu+2._r8)*(nu+4._r8)/(nu+1._r8)**2._r8* &
          (rho(i,k)*qc(i,k)/1000._r8)**4._r8/(rho(i,k)*nc(i,k)/1.e6_r8)**2._r8* &
          (1._r8+dum1/(1._r8-dum)**2)*1000._r8 / rho(i,k)

        nprc1(i,k) = au(i,k)*2._r8/2.6e-7_r8*1000._r8
        nprc(i,k) = au(i,k)/droplet_mass_40um
      else
       au(i,k) = 0._r8
       nprc1(i,k) = 0._r8
       nprc(i,k)=0._r8
      end if
    end do
  end do

  end subroutine sb2001v2_liq_autoconversion 
  
!========================================================================
!SB2001 Accretion V2

subroutine sb2001v2_accre_cld_water_rain(qc,nc,qr,rho,relvar,pra,npra,mgncol)
  !
  ! ---------------------------------------------------------------------
  ! ACCR_SB calculates the evolution of mass mxng-ratio due to accretion
  ! and self collection following Seifert & Beheng (2001).  
  !
  
  integer, intent(in) :: mgncol
  
  real(r8), dimension(mgncol), intent (in)    :: qc  ! = qc (cld water mixing ratio)
  real(r8), dimension(mgncol), intent (in)    :: nc  ! = nc (cld water number conc /kg)    
  real(r8), dimension(mgncol), intent (in)    :: qr  ! = qr (rain water mixing ratio)
  real(r8), dimension(mgncol), intent (in)    :: rho ! = rho : density profile
  real(r8), dimension(mgncol), intent (in)    :: relvar

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: pra  ! MMR
  real(r8), dimension(mgncol), intent(out) :: npra ! Number

  ! parameters for Seifert and Beheng (2001) autoconversion/accretion                                         
  real(r8), parameter :: kc = 9.44e9_r8
  real(r8), parameter :: kr = 5.78e3_r8

  real(r8) :: dum, dum1
  integer :: i

  ! accretion

  do i =1,mgncol

    if (qc(i) > qsmall) then
      dum = 1._r8-qc(i)/(qc(i)+qr(i))
      dum1 = (dum/(dum+5.e-4_r8))**4._r8
      pra(i) = kr*rho(i)*0.001_r8*qc(i)*qr(i)*dum1
      npra(i) = pra(i)*rho(i)*0.001_r8*(nc(i)*rho(i)*1.e-6_r8)/ &
           (qc(i)*rho(i)*0.001_r8)*1.e6_r8 / rho(i)
    else
      pra(i) = 0._r8
      npra(i) = 0._r8
    end if 
  
  enddo
 
end subroutine sb2001v2_accre_cld_water_rain   
  
subroutine sb2001v2_accre_cld_water_rain_2D(mgncol,nlev,qc,nc,qr,rho,relvar,pra,npra)
  !
  ! ---------------------------------------------------------------------
  ! ACCR_SB calculates the evolution of mass mxng-ratio due to accretion
  ! and self collection following Seifert & Beheng (2001).  
  !
  
  integer, intent(in) :: mgncol, nlev
  
  real(r8), dimension(mgncol,nlev), intent (in)    :: qc  ! = qc (cld water mixing ratio)
  real(r8), dimension(mgncol,nlev), intent (in)    :: nc  ! = nc (cld water number conc /kg)    
  real(r8), dimension(mgncol,nlev), intent (in)    :: qr  ! = qr (rain water mixing ratio)
  real(r8), dimension(mgncol,nlev), intent (in)    :: rho ! = rho : density profile
  real(r8), dimension(mgncol,nlev), intent (in)    :: relvar

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: pra  ! MMR
  real(r8), dimension(mgncol,nlev), intent(out) :: npra ! Number

  ! parameters for Seifert and Beheng (2001) autoconversion/accretion                                         
  real(r8), parameter :: kc = 9.44e9_r8
  real(r8), parameter :: kr = 5.78e3_r8

  real(r8) :: dum, dum1
  integer :: i, k

  ! accretion
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol

      if (qc(i,k) > qsmall) then
        dum = 1._r8-qc(i,k)/(qc(i,k)+qr(i,k))
        dum1 = (dum/(dum+5.e-4_r8))**4._r8
        pra(i,k) = kr*rho(i,k)*0.001_r8*qc(i,k)*qr(i,k)*dum1
        npra(i,k) = pra(i,k)*rho(i,k)*0.001_r8*(nc(i,k)*rho(i,k)*1.e-6_r8)/ &
             (qc(i,k)*rho(i,k)*0.001_r8)*1.e6_r8 / rho(i,k)
      else
        pra(i,k) = 0._r8
        npra(i,k) = 0._r8
      end if 
  
    end do
  end do
 
end subroutine sb2001v2_accre_cld_water_rain_2D

!========================================================================
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)

subroutine ice_autoconversion(t, qiic, lami, n0i, dcs, prci, nprci, mgncol)

  integer, intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: t
  real(r8), dimension(mgncol), intent(in) :: qiic
  real(r8), dimension(mgncol), intent(in) :: lami
  real(r8), dimension(mgncol), intent(in) :: n0i
  real(r8),                    intent(in) :: dcs

  real(r8), dimension(mgncol), intent(out) :: prci
  real(r8), dimension(mgncol), intent(out) :: nprci

  ! Assume autoconversion timescale of 180 seconds.
  real(r8), parameter :: ac_time = 180._r8

  ! Average mass of an ice particle.
  real(r8) :: m_ip
  ! Ratio of autoconversion diameter to average diameter.
  real(r8) :: d_rat
  integer :: i

  do i=1,mgncol
     if (t(i) <= tmelt .and. qiic(i) >= qsmall) then

        d_rat = lami(i)*dcs

        ! Rate of ice particle conversion (number).
        nprci(i) = n0i(i)/(lami(i)*ac_time)*exp(-d_rat)

        m_ip = (rhoi*pi/6._r8) / lami(i)**3

        ! Rate of mass conversion.
        ! Note that this is:
        ! m n (d^3 + 3 d^2 + 6 d + 6)
        prci(i) = m_ip * nprci(i) * &
             (((d_rat + 3._r8)*d_rat + 6._r8)*d_rat + 6._r8)

     else
        prci(i) = 0._r8
        nprci(i) = 0._r8
     end if
  enddo
end subroutine ice_autoconversion

subroutine ice_autoconversion_2D(mgncol, nlev, t, qiic, lami, n0i, dcs, prci, nprci)

  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: t
  real(r8), dimension(mgncol,nlev), intent(in) :: qiic
  real(r8), dimension(mgncol,nlev), intent(in) :: lami
  real(r8), dimension(mgncol,nlev), intent(in) :: n0i
  real(r8),                         intent(in) :: dcs

  real(r8), dimension(mgncol,nlev), intent(out) :: prci
  real(r8), dimension(mgncol,nlev), intent(out) :: nprci

  ! Assume autoconversion timescale of 180 seconds.
  real(r8), parameter :: ac_time = 180._r8

  ! Average mass of an ice particle.
  real(r8) :: m_ip
  ! Ratio of autoconversion diameter to average diameter.
  real(r8) :: d_rat
  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
    
      if (t(i,k) <= tmelt .and. qiic(i,k) >= qsmall) then

        d_rat = lami(i,k)*dcs

        ! Rate of ice particle conversion (number).
        nprci(i,k) = n0i(i,k)/(lami(i,k)*ac_time)*exp(-d_rat)

        m_ip = (rhoi*pi/6._r8) / lami(i,k)**3

        ! Rate of mass conversion.
        ! Note that this is:
        ! m n (d^3 + 3 d^2 + 6 d + 6)
        prci(i,k) = m_ip * nprci(i,k) * &
             (((d_rat + 3._r8)*d_rat + 6._r8)*d_rat + 6._r8)

      else
        prci(i,k) = 0._r8
        nprci(i,k) = 0._r8
      end if
    end do
  end do
  
end subroutine ice_autoconversion_2D

! immersion freezing (Bigg, 1953)
!===================================

subroutine immersion_freezing(microp_uniform, t, pgam, lamc, &
     qcic, ncic, relvar, mnuccc, nnuccc, mgncol)

  integer, intent(in) :: mgncol
  logical, intent(in) :: microp_uniform

  ! Temperature
  real(r8), dimension(mgncol), intent(in) :: t

  ! Cloud droplet size distribution parameters
  real(r8), dimension(mgncol), intent(in) :: pgam
  real(r8), dimension(mgncol), intent(in) :: lamc

  ! MMR and number concentration of in-cloud liquid water
  real(r8), dimension(mgncol), intent(in) :: qcic
  real(r8), dimension(mgncol), intent(in) :: ncic

  ! Relative variance of cloud water
  real(r8), dimension(mgncol), intent(in) :: relvar

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: mnuccc ! MMR
  real(r8), dimension(mgncol), intent(out) :: nnuccc ! Number

  ! Coefficients that will be omitted for sub-columns
  real(r8), dimension(mgncol) :: dum
  integer :: i

  if (.not. microp_uniform) then
     dum(:) = var_coef(relvar, 2)
  else
     dum(:) = 1._r8
  end if
  do i=1,mgncol

     if (qcic(i) >= qsmall .and. t(i) < 269.15_r8) then

        nnuccc(i) = &
             pi/6._r8*ncic(i)*rising_factorial(pgam(i)+1._r8, 3)* &
             bimm*(exp(aimm*(tmelt - t(i)))-1._r8)/lamc(i)**3

        mnuccc(i) = dum(i) * nnuccc(i) * &
             pi/6._r8*rhow* &
             rising_factorial(pgam(i)+4._r8, 3)/lamc(i)**3

     else
        mnuccc(i) = 0._r8
        nnuccc(i) = 0._r8
     end if ! qcic > qsmall and t < 4 deg C
  enddo

end subroutine immersion_freezing

subroutine immersion_freezing_2D(mgncol, nlev, microp_uniform, t, pgam, lamc, &
                                 qcic, ncic, relvar, mnuccc, nnuccc )

  integer, intent(in) :: mgncol, nlev
  logical, intent(in) :: microp_uniform

  ! Temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: t

  ! Cloud droplet size distribution parameters
  real(r8), dimension(mgncol,nlev), intent(in) :: pgam
  real(r8), dimension(mgncol,nlev), intent(in) :: lamc

  ! MMR and number concentration of in-cloud liquid water
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic
  real(r8), dimension(mgncol,nlev), intent(in) :: ncic

  ! Relative variance of cloud water
  real(r8), dimension(mgncol,nlev), intent(in) :: relvar

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: mnuccc ! MMR
  real(r8), dimension(mgncol,nlev), intent(out) :: nnuccc ! Number

  ! Coefficients that will be omitted for sub-columns
  real(r8), dimension(mgncol,nlev) :: dum, rising_factorial_term, pgam_tmp
  
  integer :: i, k

  !$acc data create(dum, rising_factorial_term, pgam_tmp) async(1)

  if (.not. microp_uniform) then
    !$acc update host(relvar)
    do k = 1, nlev
      do i = 1, mgncol
        dum(i,k) = var_coef(relvar(i,k), 2)
      end do
    end do
    !$acc update device(dum)
  else
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k = 1, nlev
      do i = 1, mgncol
        dum(i,k) = 1._r8
      end do
    end do
  end if
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      pgam_tmp(i,k) = pgam(i,k) + 1.0_r8
    end do
  end do
  
  rising_factorial_term(:,:) = rising_factorial(mgncol, nlev, pgam_tmp, 3)
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol

      if (qcic(i,k) >= qsmall .and. t(i,k) < 269.15_r8) then

        nnuccc(i,k) = &
             pi/6._r8*ncic(i,k)* rising_factorial_term(i,k) * &
             bimm*(exp(aimm*(tmelt - t(i,k)))-1._r8)/lamc(i,k)**3

      else
        nnuccc(i,k) = 0._r8
      end if ! qcic > qsmall and t < 4 deg C
    end do
  end do
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      pgam_tmp(i,k) = pgam(i,k) + 4.0_r8
    end do
  end do
  
  rising_factorial_term(:,:) = rising_factorial(mgncol, nlev, pgam_tmp, 3)
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol

      if (qcic(i,k) >= qsmall .and. t(i,k) < 269.15_r8) then
        mnuccc(i,k) = dum(i,k) * nnuccc(i,k) * &
             pi/6._r8*rhow* &
             rising_factorial_term(i,k)/lamc(i,k)**3

      else
        mnuccc(i,k) = 0._r8
      end if ! qcic > qsmall and t < 4 deg C
    end do
  end do
  
  !$acc end data

end subroutine immersion_freezing_2D

! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
!===================================================================
! dust size and number in multiple bins are read in from companion routine

subroutine contact_freezing (microp_uniform, t, p, rndst, nacon, &
     pgam, lamc, qcic, ncic, relvar, mnucct, nnucct, mgncol, mdust)

  logical, intent(in) :: microp_uniform

  integer, intent(in) :: mgncol
  integer, intent(in) :: mdust

  real(r8), dimension(mgncol), intent(in) :: t            ! Temperature
  real(r8), dimension(mgncol), intent(in) :: p            ! Pressure
  real(r8), dimension(mgncol, mdust), intent(in) :: rndst ! Radius (for multiple dust bins)
  real(r8), dimension(mgncol, mdust), intent(in) :: nacon ! Number (for multiple dust bins)

  ! Size distribution parameters for cloud droplets
  real(r8), dimension(mgncol), intent(in) :: pgam
  real(r8), dimension(mgncol), intent(in) :: lamc

  ! MMR and number concentration of in-cloud liquid water
  real(r8), dimension(mgncol), intent(in) :: qcic
  real(r8), dimension(mgncol), intent(in) :: ncic

  ! Relative cloud water variance
  real(r8), dimension(mgncol), intent(in) :: relvar

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: mnucct ! MMR
  real(r8), dimension(mgncol), intent(out) :: nnucct ! Number

  real(r8) :: tcnt                  ! scaled relative temperature
  real(r8) :: viscosity             ! temperature-specific viscosity (kg/m/s)
  real(r8) :: mfp                   ! temperature-specific mean free path (m)

  ! Dimension these according to number of dust bins, inferred from rndst size
  real(r8) :: nslip(size(rndst,2))  ! slip correction factors
  real(r8) :: ndfaer(size(rndst,2)) ! aerosol diffusivities (m^2/sec)

  ! Coefficients not used for subcolumns
  real(r8) :: dum, dum1

  ! Common factor between mass and number.
  real(r8) :: contact_factor

  integer  :: i

  do i = 1,mgncol

     if (qcic(i) >= qsmall .and. t(i) < 269.15_r8) then

        if (.not. microp_uniform) then
           dum = var_coef(relvar(i), 4._r8/3._r8)
           dum1 = var_coef(relvar(i), 1._r8/3._r8)
        else
           dum = 1._r8
           dum1 = 1._r8
        endif

        tcnt=(270.16_r8-t(i))**1.3_r8
        viscosity = 1.8e-5_r8*(t(i)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
        mfp = 2.0_r8*viscosity/ &                         ! Mean free path (m)
                     (p(i)*sqrt( 8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i)) ))

        ! Note that these two are vectors.
        nslip = 1.0_r8+(mfp/rndst(i,:))*(1.257_r8+(0.4_r8*exp(-(1.1_r8*rndst(i,:)/mfp))))! Slip correction factor

        ndfaer = 1.381e-23_r8*t(i)*nslip/(6._r8*pi*viscosity*rndst(i,:))  ! aerosol diffusivity (m2/s)

        contact_factor = dot_product(ndfaer,nacon(i,:)*tcnt) * pi * &
             ncic(i) * (pgam(i) + 1._r8) / lamc(i)

        mnucct(i) = dum * contact_factor * &
             pi/3._r8*rhow*rising_factorial(pgam(i)+2._r8, 3)/lamc(i)**3

        nnucct(i) =  dum1 * 2._r8 * contact_factor

     else

        mnucct(i)=0._r8
        nnucct(i)=0._r8

     end if ! qcic > qsmall and t < 4 deg C
  end do

end subroutine contact_freezing

subroutine contact_freezing_2D(mgncol, nlev, mdust, microp_uniform, t, p, rndst, nacon, &
                               pgam, lamc, qcic, ncic, relvar, mnucct, nnucct)

  logical, intent(in) :: microp_uniform

  integer, intent(in) :: mgncol, nlev, mdust

  real(r8), dimension(mgncol,nlev), intent(in) :: t            ! Temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: p            ! Pressure
  real(r8), dimension(mgncol,nlev, mdust), intent(in) :: rndst ! Radius (for multiple dust bins)
  real(r8), dimension(mgncol,nlev, mdust), intent(in) :: nacon ! Number (for multiple dust bins)

  ! Size distribution parameters for cloud droplets
  real(r8), dimension(mgncol,nlev), intent(in) :: pgam
  real(r8), dimension(mgncol,nlev), intent(in) :: lamc

  ! MMR and number concentration of in-cloud liquid water
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic
  real(r8), dimension(mgncol,nlev), intent(in) :: ncic

  ! Relative cloud water variance
  real(r8), dimension(mgncol,nlev), intent(in) :: relvar

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: mnucct ! MMR
  real(r8), dimension(mgncol,nlev), intent(out) :: nnucct ! Number

  real(r8), dimension(mgncol,nlev) :: tcnt        ! scaled relative temperature
  real(r8), dimension(mgncol,nlev) :: viscosity   ! temperature-specific viscosity (kg/m/s)
  real(r8), dimension(mgncol,nlev) :: mfp         ! temperature-specific mean free path (m)

  ! Dimension these according to number of dust bins, inferred from rndst size
  real(r8) :: nslip  ! slip correction factors
  real(r8), dimension(mgncol,nlev, mdust) :: ndfaer ! aerosol diffusivities (m^2/sec)

  ! Coefficients not used for subcolumns
  real(r8), dimension(mgncol,nlev) :: dum, dum1, dot_product_factor, &
                                      rising_factorial_term, pgam_tmp

  ! Common factor between mass and number.
  real(r8) :: contact_factor

  integer  :: i, k, m
  
  !$acc data create(dum, dum1, dot_product_factor, tcnt, viscosity, mfp, ndfaer, &
  !$acc&            rising_factorial_term, pgam_tmp) async(1)
  
  if (.not. microp_uniform) then
    !$acc update host(relvar)
    do k = 1, nlev
      do i = 1, mgncol
        dum(i,k) = var_coef(relvar(i,k), 4._r8/3._r8)
        dum1(i,k) = var_coef(relvar(i,k), 1._r8/3._r8)
      end do
    end do
    !$acc update device(dum,dum1)
  else
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k = 1, nlev
      do i = 1, mgncol
        dum(i,k) = 1._r8
        dum1(i,k) = 1._r8
      end do
    end do
  endif

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol

      tcnt(i,k) = (270.16_r8-t(i,k))**1.3_r8
      
      viscosity(i,k) = 1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
      
      mfp(i,k) = 2.0_r8*viscosity(i,k)/ &                         ! Mean free path (m)
                   (p(i,k)*sqrt( 8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k)) ))
    end do
  end do

  !$acc parallel loop collapse(2) default(present) async(1) 
  do m = 1, mdust
    do k = 1, nlev
      do i = 1, mgncol
        
        ! Slip correction factor
        nslip = 1.0_r8+(mfp(i,k)/rndst(i,k,m)) &
                        * (1.257_r8+(0.4_r8*exp(-(1.1_r8*rndst(i,k,m)/mfp(i,k)))))

        ! aerosol diffusivity (m2/s)
        ndfaer(i,k,m) = 1.381e-23_r8*t(i,k)*nslip &
                          /(6._r8*pi*viscosity(i,k)*rndst(i,k,m))  
      end do
    end do
  end do

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      dot_product_factor(i,k) = 0.0_r8
      do m = 1, mdust
        dot_product_factor(i,k) = dot_product_factor(i,k) &
                                  + ndfaer(i,k,m) * nacon(i,k,m) * tcnt(i,k)
      end do
    end do
  end do
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      pgam_tmp(i,k) = pgam(i,k) + 2.0_r8
    end do
  end do
  
  rising_factorial_term(:,:) = rising_factorial(mgncol, nlev, pgam_tmp, 3)

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      contact_factor = dot_product_factor(i,k) * pi * ncic(i,k) * (pgam(i,k) + 1._r8) / lamc(i,k)

      mnucct(i,k) = dum(i,k) * contact_factor * pi/3._r8 * rhow &
                    * rising_factorial_term(i,k) / lamc(i,k)**3

      nnucct(i,k) =  dum1(i,k) * 2._r8 * contact_factor
    end do
  end do
  
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qcic(i,k) < qsmall .or. t(i,k) >= 269.15_r8) then
        mnucct(i,k)=0._r8
        nnucct(i,k)=0._r8
      end if
    end do
  end do
  
  !$acc end data

end subroutine contact_freezing_2D

! snow self-aggregation from passarelli, 1978, used by reisner, 1998
!===================================================================
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

subroutine snow_self_aggregation(t, rho, asn, rhosn, qsic, nsic, nsagg, mgncol)

  integer,                          intent(in) :: mgncol

  real(r8), dimension(mgncol), intent(in) :: t     ! Temperature
  real(r8), dimension(mgncol), intent(in) :: rho   ! Density
  real(r8), dimension(mgncol), intent(in) :: asn   ! fall speed parameter for snow
  real(r8),                    intent(in) :: rhosn ! density of snow

  ! In-cloud snow
  real(r8), dimension(mgncol), intent(in) :: qsic ! MMR
  real(r8), dimension(mgncol), intent(in) :: nsic ! Number

  ! Output number tendency
  real(r8), dimension(mgncol), intent(out) :: nsagg

  integer :: i

  do i=1,mgncol
     if (qsic(i) >= qsmall .and. t(i) <= tmelt) then
        nsagg(i) = -1108._r8*eii/(4._r8*720._r8*rhosn)*asn(i)*qsic(i)*nsic(i)*rho(i)*&
             ((qsic(i)/nsic(i))*(1._r8/(rhosn*pi)))**((bs-1._r8)/3._r8)
     else
        nsagg(i)=0._r8
     end if
  enddo
end subroutine snow_self_aggregation

subroutine snow_self_aggregation_2D(mgncol, nlev, t, rho, asn, rhosn, qsic, nsic, nsagg)

  integer, intent(in) :: mgncol, nlev

  real(r8), dimension(mgncol,nlev), intent(in) :: t     ! Temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: rho   ! Density
  real(r8), dimension(mgncol,nlev), intent(in) :: asn   ! fall speed parameter for snow
  real(r8),                    intent(in) :: rhosn ! density of snow

  ! In-cloud snow
  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: nsic ! Number

  ! Output number tendency
  real(r8), dimension(mgncol,nlev), intent(out) :: nsagg

  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qsic(i,k) >= qsmall .and. t(i,k) <= tmelt) then
        nsagg(i,k) = -1108._r8*eii/(4._r8*720._r8*rhosn)*asn(i,k)*qsic(i,k)*nsic(i,k)*rho(i,k)*&
             ((qsic(i,k)/nsic(i,k))*(1._r8/(rhosn*pi)))**((bs-1._r8)/3._r8)
      else
        nsagg(i,k)=0._r8
      end if
    end do
  end do
  
end subroutine snow_self_aggregation_2D

! accretion of cloud droplets onto snow/graupel
!===================================================================
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice
! since minimum size ice particle for accretion is 50 - 150 micron

subroutine accrete_cloud_water_snow(t, rho, asn, uns, mu, qcic, ncic, qsic, &
     pgam, lamc, lams, n0s, psacws, npsacws, mgncol)

  integer, intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: t   ! Temperature
  real(r8), dimension(mgncol), intent(in) :: rho ! Density
  real(r8), dimension(mgncol), intent(in) :: asn ! Fallspeed parameter (snow)
  real(r8), dimension(mgncol), intent(in) :: uns ! Current fallspeed   (snow)
  real(r8), dimension(mgncol), intent(in) :: mu  ! Viscosity

  ! In-cloud liquid water
  real(r8), dimension(mgncol), intent(in) :: qcic ! MMR
  real(r8), dimension(mgncol), intent(in) :: ncic ! Number

  ! In-cloud snow
  real(r8), dimension(mgncol), intent(in) :: qsic ! MMR

  ! Cloud droplet size parameters
  real(r8), dimension(mgncol), intent(in) :: pgam
  real(r8), dimension(mgncol), intent(in) :: lamc

  ! Snow size parameters
  real(r8), dimension(mgncol), intent(in) :: lams
  real(r8), dimension(mgncol), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: psacws  ! Mass mixing ratio
  real(r8), dimension(mgncol), intent(out) :: npsacws ! Number concentration

  real(r8) :: dc0 ! Provisional mean droplet size
  real(r8) :: dum
  real(r8) :: eci ! collection efficiency for riming of snow by droplets

  ! Fraction of cloud droplets accreted per second
  real(r8) :: accrete_rate
  integer :: i

  ! ignore collision of snow with droplets above freezing

  do i=1,mgncol
     if (qsic(i) >= qsmall .and. t(i) <= tmelt .and. qcic(i) >= qsmall) then

        ! put in size dependent collection efficiency
        ! mean diameter of snow is area-weighted, since
        ! accretion is function of crystal geometric area
        ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

        dc0 = (pgam(i)+1._r8)/lamc(i)
        dum = dc0*dc0*uns(i)*rhow*lams(i)/(9._r8*mu(i))
        eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

        eci = max(eci,0._r8)
        eci = min(eci,1._r8)

        ! no impact of sub-grid distribution of qc since psacws
        ! is linear in qc
        accrete_rate = pi/4._r8*asn(i)*rho(i)*n0s(i)*eci*gamma_bs_plus3 / lams(i)**(bs+3._r8)
        psacws(i) = accrete_rate*qcic(i)
        npsacws(i) = accrete_rate*ncic(i)
     else
        psacws(i) = 0._r8
        npsacws(i) = 0._r8
     end if
  enddo
end subroutine accrete_cloud_water_snow

subroutine accrete_cloud_water_snow_2D(mgncol, nlev, t, rho, asn, uns, mu, qcic, ncic, qsic, &
                                       pgam, lamc, lams, n0s, psacws, npsacws )

  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: t   ! Temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: rho ! Density
  real(r8), dimension(mgncol,nlev), intent(in) :: asn ! Fallspeed parameter (snow)
  real(r8), dimension(mgncol,nlev), intent(in) :: uns ! Current fallspeed   (snow)
  real(r8), dimension(mgncol,nlev), intent(in) :: mu  ! Viscosity

  ! In-cloud liquid water
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic ! MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: ncic ! Number

  ! In-cloud snow
  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! MMR

  ! Cloud droplet size parameters
  real(r8), dimension(mgncol,nlev), intent(in) :: pgam
  real(r8), dimension(mgncol,nlev), intent(in) :: lamc

  ! Snow size parameters
  real(r8), dimension(mgncol,nlev), intent(in) :: lams
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: psacws  ! Mass mixing ratio
  real(r8), dimension(mgncol,nlev), intent(out) :: npsacws ! Number concentration

  real(r8) :: dc0 ! Provisional mean droplet size
  real(r8) :: dum
  real(r8) :: eci ! collection efficiency for riming of snow by droplets

  ! Fraction of cloud droplets accreted per second
  real(r8) :: accrete_rate
  integer :: i, k

  ! ignore collision of snow with droplets above freezing

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qsic(i,k) >= qsmall .and. t(i,k) <= tmelt .and. qcic(i,k) >= qsmall) then

        ! put in size dependent collection efficiency
        ! mean diameter of snow is area-weighted, since
        ! accretion is function of crystal geometric area
        ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

        dc0 = (pgam(i,k)+1._r8)/lamc(i,k)
        dum = dc0*dc0*uns(i,k)*rhow*lams(i,k)/(9._r8*mu(i,k))
        eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

        eci = max(eci,0._r8)
        eci = min(eci,1._r8)

        ! no impact of sub-grid distribution of qc since psacws
        ! is linear in qc
        accrete_rate = pi/4._r8*asn(i,k)*rho(i,k)*n0s(i,k)*eci*gamma_bs_plus3 / lams(i,k)**(bs+3._r8)
        psacws(i,k) = accrete_rate*qcic(i,k)
        npsacws(i,k) = accrete_rate*ncic(i,k)
      else
        psacws(i,k) = 0._r8
        npsacws(i,k) = 0._r8
      end if
    end do
  end do
  
end subroutine accrete_cloud_water_snow_2D

! add secondary ice production due to accretion of droplets by snow
!===================================================================
! (Hallet-Mossop process) (from Cotton et al., 1986)

subroutine secondary_ice_production(t, psacws, msacwi, nsacwi, mgncol)

  integer, intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: t ! Temperature

  ! Accretion of cloud water to snow tendencies
  real(r8), dimension(mgncol), intent(inout) :: psacws ! MMR

  ! Output (ice) tendencies
  real(r8), dimension(mgncol), intent(out) :: msacwi ! MMR
  real(r8), dimension(mgncol), intent(out) :: nsacwi ! Number
  integer :: i

  do i=1,mgncol
     if((t(i) < 270.16_r8) .and. (t(i) >= 268.16_r8)) then
        nsacwi(i) = 3.5e8_r8*(270.16_r8-t(i))/2.0_r8*psacws(i)
     else if((t(i) < 268.16_r8) .and. (t(i) >= 265.16_r8)) then
        nsacwi(i) = 3.5e8_r8*(t(i)-265.16_r8)/3.0_r8*psacws(i)
     else
        nsacwi(i) = 0.0_r8
     endif
  enddo

  do i=1,mgncol
     msacwi(i) = min(nsacwi(i)*mi0, psacws(i))
     psacws(i) = psacws(i) - msacwi(i)
  enddo
end subroutine secondary_ice_production

subroutine secondary_ice_production_2D(mgncol, nlev, t, psacws, msacwi, nsacwi)

  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: t ! Temperature

  ! Accretion of cloud water to snow tendencies
  real(r8), dimension(mgncol,nlev), intent(inout) :: psacws ! MMR

  ! Output (ice) tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: msacwi ! MMR
  real(r8), dimension(mgncol,nlev), intent(out) :: nsacwi ! Number
  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if((t(i,k) < 270.16_r8) .and. (t(i,k) >= 268.16_r8)) then
        nsacwi(i,k) = 3.5e8_r8*(270.16_r8-t(i,k))/2.0_r8*psacws(i,k)
      else if((t(i,k) < 268.16_r8) .and. (t(i,k) >= 265.16_r8)) then
        nsacwi(i,k) = 3.5e8_r8*(t(i,k)-265.16_r8)/3.0_r8*psacws(i,k)
      else
        nsacwi(i,k) = 0.0_r8
      endif
    end do
  end do

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      msacwi(i,k) = min(nsacwi(i,k)*mi0, psacws(i,k))
      psacws(i,k) = psacws(i,k) - msacwi(i,k)
    end do
  end do
end subroutine secondary_ice_production_2D

! accretion of rain water by snow
!===================================================================
! formula from ikawa and saito, 1991, used by reisner et al., 1998

subroutine accrete_rain_snow(t, rho, umr, ums, unr, uns, qric, qsic, &
     lamr, n0r, lams, n0s, pracs, npracs, mgncol)

  integer,                          intent(in) :: mgncol

  real(r8), dimension(mgncol), intent(in) :: t   ! Temperature
  real(r8), dimension(mgncol), intent(in) :: rho ! Density

  ! Fallspeeds
  ! mass-weighted
  real(r8), dimension(mgncol), intent(in) :: umr ! rain
  real(r8), dimension(mgncol), intent(in) :: ums ! snow
  ! number-weighted
  real(r8), dimension(mgncol), intent(in) :: unr ! rain
  real(r8), dimension(mgncol), intent(in) :: uns ! snow

  ! In cloud MMRs
  real(r8), dimension(mgncol), intent(in) :: qric ! rain
  real(r8), dimension(mgncol), intent(in) :: qsic ! snow

  ! Size distribution parameters
  ! rain
  real(r8), dimension(mgncol), intent(in) :: lamr
  real(r8), dimension(mgncol), intent(in) :: n0r
  ! snow
  real(r8), dimension(mgncol), intent(in) :: lams
  real(r8), dimension(mgncol), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: pracs  ! MMR
  real(r8), dimension(mgncol), intent(out) :: npracs ! Number

  ! Collection efficiency for accretion of rain by snow
  real(r8), parameter :: ecr = 1.0_r8

  ! Ratio of average snow diameter to average rain diameter.
  real(r8) :: d_rat
  ! Common factor between mass and number expressions
  real(r8) :: common_factor
  integer :: i

  do i=1,mgncol
     if (qric(i) >= icsmall .and. qsic(i) >= icsmall .and. t(i) <= tmelt) then

        common_factor = pi*ecr*rho(i)*n0r(i)*n0s(i)/(lamr(i)**3 * lams(i))

        d_rat = lamr(i)/lams(i)

        pracs(i) = common_factor*pi*rhow* &
             sqrt((1.2_r8*umr(i)-0.95_r8*ums(i))**2 + 0.08_r8*ums(i)*umr(i)) / lamr(i)**3 * &
             ((0.5_r8*d_rat + 2._r8)*d_rat + 5._r8)

        npracs(i) = common_factor*0.5_r8* &
             sqrt(1.7_r8*(unr(i)-uns(i))**2 + 0.3_r8*unr(i)*uns(i)) * &
             ((d_rat + 1._r8)*d_rat + 1._r8)

     else
        pracs(i) = 0._r8
        npracs(i) = 0._r8
     end if
  enddo
end subroutine accrete_rain_snow

subroutine accrete_rain_snow_2D(mgncol, nlev, t, rho, umr, ums, unr, uns, qric, qsic, &
                                lamr, n0r, lams, n0s, pracs, npracs)

  integer, intent(in) :: mgncol, nlev

  real(r8), dimension(mgncol,nlev), intent(in) :: t   ! Temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: rho ! Density

  ! Fallspeeds
  ! mass-weighted
  real(r8), dimension(mgncol,nlev), intent(in) :: umr ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: ums ! snow
  ! number-weighted
  real(r8), dimension(mgncol,nlev), intent(in) :: unr ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: uns ! snow

  ! In cloud MMRs
  real(r8), dimension(mgncol,nlev), intent(in) :: qric ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! snow

  ! Size distribution parameters
  ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: lamr
  real(r8), dimension(mgncol,nlev), intent(in) :: n0r
  ! snow
  real(r8), dimension(mgncol,nlev), intent(in) :: lams
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: pracs  ! MMR
  real(r8), dimension(mgncol,nlev), intent(out) :: npracs ! Number

  ! Collection efficiency for accretion of rain by snow
  real(r8), parameter :: ecr = 1.0_r8

  ! Ratio of average snow diameter to average rain diameter.
  real(r8) :: d_rat
  ! Common factor between mass and number expressions
  real(r8) :: common_factor
  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qric(i,k) >= icsmall .and. qsic(i,k) >= icsmall .and. t(i,k) <= tmelt) then

        common_factor = pi*ecr*rho(i,k)*n0r(i,k)*n0s(i,k)/(lamr(i,k)**3 * lams(i,k))

        d_rat = lamr(i,k)/lams(i,k)

        pracs(i,k) = common_factor*pi*rhow* &
             sqrt((1.2_r8*umr(i,k)-0.95_r8*ums(i,k))**2 + 0.08_r8*ums(i,k)*umr(i,k)) / lamr(i,k)**3 * &
             ((0.5_r8*d_rat + 2._r8)*d_rat + 5._r8)

        npracs(i,k) = common_factor*0.5_r8* &
             sqrt(1.7_r8*(unr(i,k)-uns(i,k))**2 + 0.3_r8*unr(i,k)*uns(i,k)) * &
             ((d_rat + 1._r8)*d_rat + 1._r8)

      else
        pracs(i,k) = 0._r8
        npracs(i,k) = 0._r8
      end if
    end do
  end do
  
end subroutine accrete_rain_snow_2D

! heterogeneous freezing of rain drops
!===================================================================
! follows from Bigg (1953)

subroutine heterogeneous_rain_freezing(t, qric, nric, lamr, mnuccr, nnuccr, mgncol)

  integer,                          intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: t    ! Temperature

  ! In-cloud rain
  real(r8), dimension(mgncol), intent(in) :: qric ! MMR
  real(r8), dimension(mgncol), intent(in) :: nric ! Number
  real(r8), dimension(mgncol), intent(in) :: lamr ! size parameter

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: mnuccr ! MMR
  real(r8), dimension(mgncol), intent(out) :: nnuccr ! Number
  integer :: i

  do i=1,mgncol

     if (t(i) < 269.15_r8 .and. qric(i) >= qsmall) then
        nnuccr(i) = pi*nric(i)*bimm* &
             (exp(aimm*(tmelt - t(i)))-1._r8)/lamr(i)**3

        mnuccr(i) = nnuccr(i) * 20._r8*pi*rhow/lamr(i)**3

     else
        mnuccr(i) = 0._r8
        nnuccr(i) = 0._r8
     end if
  enddo
end subroutine heterogeneous_rain_freezing

subroutine heterogeneous_rain_freezing_2D(mgncol, nlev, t, qric, nric, lamr, mnuccr, nnuccr )

  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: t    ! Temperature

  ! In-cloud rain
  real(r8), dimension(mgncol,nlev), intent(in) :: qric ! MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: nric ! Number
  real(r8), dimension(mgncol,nlev), intent(in) :: lamr ! size parameter

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: mnuccr ! MMR
  real(r8), dimension(mgncol,nlev), intent(out) :: nnuccr ! Number
  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol

      if (t(i,k) < 269.15_r8 .and. qric(i,k) >= qsmall) then
        nnuccr(i,k) = pi*nric(i,k)*bimm* &
             (exp(aimm*(tmelt - t(i,k)))-1._r8)/lamr(i,k)**3

        mnuccr(i,k) = nnuccr(i,k) * 20._r8*pi*rhow/lamr(i,k)**3

      else
        mnuccr(i,k) = 0._r8
        nnuccr(i,k) = 0._r8
      end if
    end do
  end do
  
end subroutine heterogeneous_rain_freezing_2D

! accretion of cloud liquid water by rain
!===================================================================
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected

subroutine accrete_cloud_water_rain(microp_uniform, qric, qcic, &
     ncic, relvar, accre_enhan, pra, npra, mgncol)

  logical, intent(in) :: microp_uniform
  integer, intent(in) :: mgncol
  ! In-cloud rain
  real(r8), dimension(mgncol), intent(in) :: qric ! MMR

  ! Cloud droplets
  real(r8), dimension(mgncol), intent(in) :: qcic ! MMR
  real(r8), dimension(mgncol), intent(in) :: ncic ! Number

  ! SGS variability
  real(r8), dimension(mgncol), intent(in) :: relvar
  real(r8), dimension(mgncol), intent(in) :: accre_enhan

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: pra  ! MMR
  real(r8), dimension(mgncol), intent(out) :: npra ! Number

  ! Coefficient that varies for subcolumns
  real(r8), dimension(mgncol) :: pra_coef

  integer :: i

  if (.not. microp_uniform) then
    pra_coef(:) = accre_enhan * var_coef(relvar(:), 1.15_r8)
  else
    pra_coef(:) = 1._r8
  end if

  do i=1,mgncol

    if (qric(i) >= qsmall .and. qcic(i) >= qsmall) then

      ! include sub-grid distribution of cloud water
      pra(i) = pra_coef(i) * 67._r8*(qcic(i)*qric(i))**1.15_r8

      npra(i) = pra(i)*ncic(i)/qcic(i)

    else
      pra(i) = 0._r8
      npra(i) = 0._r8
    end if
  end do
end subroutine accrete_cloud_water_rain

subroutine accrete_cloud_water_rain_2D(mgncol, nlev, microp_uniform, qric, qcic, &
                                       ncic, relvar, accre_enhan, pra, npra )

  logical, intent(in) :: microp_uniform
  integer, intent(in) :: mgncol, nlev
  ! In-cloud rain
  real(r8), dimension(mgncol,nlev), intent(in) :: qric ! MMR

  ! Cloud droplets
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic ! MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: ncic ! Number

  ! SGS variability
  real(r8), dimension(mgncol,nlev), intent(in) :: relvar
  real(r8), dimension(mgncol,nlev), intent(in) :: accre_enhan

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: pra  ! MMR
  real(r8), dimension(mgncol,nlev), intent(out) :: npra ! Number

  ! Coefficient that varies for subcolumns
  real(r8), dimension(mgncol,nlev) :: pra_coef

  integer :: i, k
  
  !$acc data create(pra_coef) async(1)

  if (.not. microp_uniform) then
    !$acc update host(relvar, accre_enhan)
    pra_coef(:,:) = accre_enhan(:,:) * var_coef(relvar(:,:), 1.15_r8)
    !$acc update device(pra_coef)
  else
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k = 1, nlev
      do i = 1, mgncol
        pra_coef(i,k) = 1._r8
      end do
    end do
  end if

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol

      if (qric(i,k) >= qsmall .and. qcic(i,k) >= qsmall) then

        ! include sub-grid distribution of cloud water
        pra(i,k) = pra_coef(i,k) * 67._r8*(qcic(i,k)*qric(i,k))**1.15_r8

        npra(i,k) = pra(i,k)*ncic(i,k)/qcic(i,k)

      else
        pra(i,k) = 0._r8
        npra(i,k) = 0._r8
      end if
    end do
  end do
  
  !$acc end data
  
end subroutine accrete_cloud_water_rain_2D

! Self-collection of rain drops
!===================================================================
! from Beheng(1994)

subroutine self_collection_rain(rho, qric, nric, nragg, mgncol)

  integer,                          intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: rho  ! Air density

  ! Rain
  real(r8), dimension(mgncol), intent(in) :: qric ! MMR
  real(r8), dimension(mgncol), intent(in) :: nric ! Number

  ! Output number tendency
  real(r8), dimension(mgncol), intent(out) :: nragg

  integer :: i

  do i=1,mgncol
     if (qric(i) >= qsmall) then
        nragg(i) = -8._r8*nric(i)*qric(i)*rho(i)
     else
        nragg(i) = 0._r8
     end if
  enddo
end subroutine self_collection_rain

subroutine self_collection_rain_2D(mgncol, nlev, rho, qric, nric, nragg)

  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: rho  ! Air density

  ! Rain
  real(r8), dimension(mgncol,nlev), intent(in) :: qric ! MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: nric ! Number

  ! Output number tendency
  real(r8), dimension(mgncol,nlev), intent(out) :: nragg

  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qric(i,k) >= qsmall) then
        nragg(i,k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
      else
        nragg(i,k) = 0._r8
      end if
    end do
  end do
  
end subroutine self_collection_rain_2D


! Accretion of cloud ice by snow
!===================================================================
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection

subroutine accrete_cloud_ice_snow(t, rho, asn, qiic, niic, qsic, &
     lams, n0s, prai, nprai, mgncol)

  integer,                          intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: t    ! Temperature
  real(r8), dimension(mgncol), intent(in) :: rho   ! Density

  real(r8), dimension(mgncol), intent(in) :: asn  ! Snow fallspeed parameter

  ! Cloud ice
  real(r8), dimension(mgncol), intent(in) :: qiic ! MMR
  real(r8), dimension(mgncol), intent(in) :: niic ! Number

  real(r8), dimension(mgncol), intent(in) :: qsic ! Snow MMR

  ! Snow size parameters
  real(r8), dimension(mgncol), intent(in) :: lams
  real(r8), dimension(mgncol), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: prai ! MMR
  real(r8), dimension(mgncol), intent(out) :: nprai ! Number

  ! Fraction of cloud ice particles accreted per second
  real(r8) :: accrete_rate

  integer :: i

  do i=1,mgncol
     if (qsic(i) >= qsmall .and. qiic(i) >= qsmall .and. t(i) <= tmelt) then

        accrete_rate = pi/4._r8 * eii * asn(i) * rho(i) * n0s(i) * gamma_bs_plus3/ &
             lams(i)**(bs+3._r8)

        prai(i) = accrete_rate * qiic(i)
        nprai(i) = accrete_rate * niic(i)

     else
        prai(i) = 0._r8
        nprai(i) = 0._r8
     end if
  enddo
end subroutine accrete_cloud_ice_snow

subroutine accrete_cloud_ice_snow_2D(mgncol, nlev, t, rho, asn, qiic, niic, qsic, &
     lams, n0s, prai, nprai)

  integer, intent(in) :: mgncol, nlev
  real(r8), dimension(mgncol,nlev), intent(in) :: t    ! Temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: rho   ! Density

  real(r8), dimension(mgncol,nlev), intent(in) :: asn  ! Snow fallspeed parameter

  ! Cloud ice
  real(r8), dimension(mgncol,nlev), intent(in) :: qiic ! MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: niic ! Number

  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! Snow MMR

  ! Snow size parameters
  real(r8), dimension(mgncol,nlev), intent(in) :: lams
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: prai ! MMR
  real(r8), dimension(mgncol,nlev), intent(out) :: nprai ! Number

  ! Fraction of cloud ice particles accreted per second
  real(r8) :: accrete_rate

  integer :: i, k

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qsic(i,k) >= qsmall .and. qiic(i,k) >= qsmall .and. t(i,k) <= tmelt) then

        accrete_rate = pi/4._r8 * eii * asn(i,k) * rho(i,k) * n0s(i,k) * gamma_bs_plus3/ &
             lams(i,k)**(bs+3._r8)

        prai(i,k) = accrete_rate * qiic(i,k)
        nprai(i,k) = accrete_rate * niic(i,k)

      else
        prai(i,k) = 0._r8
        nprai(i,k) = 0._r8
      end if
    end do
  end do
  
end subroutine accrete_cloud_ice_snow_2D

! calculate evaporation/sublimation of rain and snow
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

subroutine evaporate_sublimate_precip(t, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, precip_frac, arn, asn, qcic, qiic, qric, qsic, lamr, n0r, lams, n0s, &
     pre, prds, am_evp_st, mgncol)

  integer,  intent(in) :: mgncol

  real(r8), dimension(mgncol), intent(in) :: t    ! temperature
  real(r8), dimension(mgncol), intent(in) :: rho  ! air density
  real(r8), dimension(mgncol), intent(in) :: dv   ! water vapor diffusivity
  real(r8), dimension(mgncol), intent(in) :: mu   ! viscosity
  real(r8), dimension(mgncol), intent(in) :: sc   ! schmidt number
  real(r8), dimension(mgncol), intent(in) :: q    ! humidity
  real(r8), dimension(mgncol), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), dimension(mgncol), intent(in) :: qvi  ! saturation humidity (ice)
  real(r8), dimension(mgncol), intent(in) :: lcldm  ! liquid cloud fraction
  real(r8), dimension(mgncol), intent(in) :: precip_frac ! precipitation fraction (maximum overlap)

  ! fallspeed parameters
  real(r8), dimension(mgncol), intent(in) :: arn  ! rain
  real(r8), dimension(mgncol), intent(in) :: asn  ! snow

  ! In-cloud MMRs
  real(r8), dimension(mgncol), intent(in) :: qcic ! cloud liquid
  real(r8), dimension(mgncol), intent(in) :: qiic ! cloud ice
  real(r8), dimension(mgncol), intent(in) :: qric ! rain
  real(r8), dimension(mgncol), intent(in) :: qsic ! snow

  ! Size parameters
  ! rain
  real(r8), dimension(mgncol), intent(in) :: lamr
  real(r8), dimension(mgncol), intent(in) :: n0r
  ! snow
  real(r8), dimension(mgncol), intent(in) :: lams
  real(r8), dimension(mgncol), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: pre
  real(r8), dimension(mgncol), intent(out) :: prds
  real(r8), dimension(mgncol), intent(out) :: am_evp_st ! Fractional area where rain evaporates.

  real(r8) :: qclr   ! water vapor mixing ratio in clear air
  real(r8) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  real(r8), dimension(mgncol) :: dum

  integer :: i

  am_evp_st = 0._r8
  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise
  do i=1,mgncol
     if (qcic(i)+qiic(i) < 1.e-6_r8) then
        dum(i) = 0._r8
     else
        dum(i) = lcldm(i)
     end if
  enddo
  do i=1,mgncol
  ! only calculate if there is some precip fraction > cloud fraction

     if (precip_frac(i) > dum(i)) then

        if (qric(i) >= qsmall .or. qsic(i) >= qsmall) then
           am_evp_st(i) = precip_frac(i) - dum(i)

           ! calculate q for out-of-cloud region
           qclr=(q(i)-dum(i)*qvl(i))/(1._r8-dum(i))
        end if

        ! evaporation of rain
        if (qric(i) >= qsmall) then

           ab = calc_ab(t(i), qvl(i), xxlv)
           eps = 2._r8*pi*n0r(i)*rho(i)*Dv(i)* &
                (f1r/(lamr(i)*lamr(i))+ &
                f2r*(arn(i)*rho(i)/mu(i))**0.5_r8* &
                sc(i)**(1._r8/3._r8)*gamma_half_br_plus5/ &
                (lamr(i)**(5._r8/2._r8+br/2._r8)))

           pre(i) = eps*(qclr-qvl(i))/ab

           ! only evaporate in out-of-cloud region
           ! and distribute across precip_frac
           pre(i)=min(pre(i)*am_evp_st(i),0._r8)
           pre(i)=pre(i)/precip_frac(i)
        else
           pre(i) = 0._r8
        end if

        ! sublimation of snow
        if (qsic(i) >= qsmall) then
           ab = calc_ab(t(i), qvi(i), xxls)
           eps = 2._r8*pi*n0s(i)*rho(i)*Dv(i)* &
                (f1s/(lams(i)*lams(i))+ &
                f2s*(asn(i)*rho(i)/mu(i))**0.5_r8* &
                sc(i)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
                (lams(i)**(5._r8/2._r8+bs/2._r8)))
           prds(i) = eps*(qclr-qvi(i))/ab

           ! only sublimate in out-of-cloud region and distribute over precip_frac
           prds(i)=min(prds(i)*am_evp_st(i),0._r8)
           prds(i)=prds(i)/precip_frac(i)
        else
           prds(i) = 0._r8
        end if

     else
        prds(i) = 0._r8
        pre(i) = 0._r8
     end if
  enddo

end subroutine evaporate_sublimate_precip

subroutine evaporate_sublimate_precip_2D(mgncol, nlev, t, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, precip_frac, arn, asn, qcic, qiic, qric, qsic, lamr, n0r, lams, n0s, &
     pre, prds, am_evp_st)

  integer,  intent(in) :: mgncol, nlev

  real(r8), dimension(mgncol,nlev), intent(in) :: t    ! temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: rho  ! air density
  real(r8), dimension(mgncol,nlev), intent(in) :: dv   ! water vapor diffusivity
  real(r8), dimension(mgncol,nlev), intent(in) :: mu   ! viscosity
  real(r8), dimension(mgncol,nlev), intent(in) :: sc   ! schmidt number
  real(r8), dimension(mgncol,nlev), intent(in) :: q    ! humidity
  real(r8), dimension(mgncol,nlev), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), dimension(mgncol,nlev), intent(in) :: qvi  ! saturation humidity (ice)
  real(r8), dimension(mgncol,nlev), intent(in) :: lcldm  ! liquid cloud fraction
  real(r8), dimension(mgncol,nlev), intent(in) :: precip_frac ! precipitation fraction (maximum overlap)

  ! fallspeed parameters
  real(r8), dimension(mgncol,nlev), intent(in) :: arn  ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: asn  ! snow

  ! In-cloud MMRs
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic ! cloud liquid
  real(r8), dimension(mgncol,nlev), intent(in) :: qiic ! cloud ice
  real(r8), dimension(mgncol,nlev), intent(in) :: qric ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! snow

  ! Size parameters
  ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: lamr
  real(r8), dimension(mgncol,nlev), intent(in) :: n0r
  ! snow
  real(r8), dimension(mgncol,nlev), intent(in) :: lams
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: pre
  real(r8), dimension(mgncol,nlev), intent(out) :: prds
  real(r8), dimension(mgncol,nlev), intent(out) :: am_evp_st ! Fractional area where rain evaporates.

  real(r8), dimension(mgncol,nlev) :: qclr   ! water vapor mixing ratio in clear air
  real(r8), dimension(mgncol,nlev) :: ab     ! correction to account for latent heat
  real(r8), dimension(mgncol,nlev) :: eps    ! 1/ sat relaxation timescale

  real(r8), dimension(mgncol,nlev) :: dum

  integer :: i, k
  
  !$acc data create(qclr, ab, eps, dum) async(1)

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
      if (qcic(i,k)+qiic(i,k) < 1.e-6_r8) then
        dum(i,k) = 0._r8
      else
        dum(i,k) = lcldm(i,k)
      end if
    end do
  end do
  
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
     ! calculate q for out-of-cloud region
     if (precip_frac(i,k) > dum(i,k) .and. &
         (qric(i,k) >= qsmall .or. qsic(i,k) >= qsmall)) then
         
         am_evp_st(i,k) = precip_frac(i,k) - dum(i,k)

         ! calculate q for out-of-cloud region
         qclr(i,k) = (q(i,k)-dum(i,k)*qvl(i,k))/(1._r8-dum(i,k))
      else
        am_evp_st(i,k) = 0._r8
      end if
    end do
  end do
  
  ab(:,:) = calc_ab_2D(mgncol, nlev, t(:,:), qvl(:,:), xxlv)
  
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
      if (qric(i,k) >= qsmall) then
        eps(i,k) = 2._r8*pi*n0r(i,k)*rho(i,k)*Dv(i,k)* &
             (f1r/(lamr(i,k)*lamr(i,k))+ &
             f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
             sc(i,k)**(1._r8/3._r8)*gamma_half_br_plus5/ &
             (lamr(i,k)**(5._r8/2._r8+br/2._r8)))
      end if
    end do
  end do
  
  ! evaporation of rain
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
      ! only calculate if there is some precip fraction > cloud fraction
      if (precip_frac(i,k) > dum(i,k) .and. qric(i,k) >= qsmall) then
         pre(i,k) = eps(i,k)*(qclr(i,k)-qvl(i,k))/ab(i,k)

         ! only evaporate in out-of-cloud region
         ! and distribute across precip_frac
         pre(i,k)=min(pre(i,k)*am_evp_st(i,k),0._r8)
         pre(i,k)=pre(i,k)/precip_frac(i,k)
      else
        pre(i,k) = 0._r8
      end if
    end do
  end do
  
  ab(:,:) = calc_ab_2D(mgncol, nlev, t(:,:), qvi(:,:), xxls)
  
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
      if (qsic(i,k) >= qsmall) then
        eps(i,k) = 2._r8*pi*n0s(i,k)*rho(i,k)*Dv(i,k)* &
             (f1s/(lams(i,k)*lams(i,k))+ &
             f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
             sc(i,k)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
             (lams(i,k)**(5._r8/2._r8+bs/2._r8)))
      end if
    end do
  end do

  ! sublimation of snow
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
      if (precip_frac(i,k) > dum(i,k) .and. qsic(i,k) >= qsmall) then
         prds(i,k) = eps(i,k)*(qclr(i,k)-qvi(i,k))/ab(i,k)

         ! only sublimate in out-of-cloud region and distribute over precip_frac
         prds(i,k)=min(prds(i,k)*am_evp_st(i,k),0._r8)
         prds(i,k)=prds(i,k)/precip_frac(i,k)
      else
        prds(i,k) = 0._r8
      end if
    end do
  end do
  
  !$acc end data

end subroutine evaporate_sublimate_precip_2D

! evaporation/sublimation of rain, snow and graupel
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

subroutine evaporate_sublimate_precip_graupel(t, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, precip_frac, arn, asn, agn, bg, qcic, qiic, qric, qsic, qgic, lamr, n0r, lams, n0s, lamg, n0g, &
     pre, prds, prdg, am_evp_st, mgncol)

  integer,  intent(in) :: mgncol

  real(r8), dimension(mgncol), intent(in) :: t    ! temperature
  real(r8), dimension(mgncol), intent(in) :: rho  ! air density
  real(r8), dimension(mgncol), intent(in) :: dv   ! water vapor diffusivity
  real(r8), dimension(mgncol), intent(in) :: mu   ! viscosity
  real(r8), dimension(mgncol), intent(in) :: sc   ! schmidt number
  real(r8), dimension(mgncol), intent(in) :: q    ! humidity
  real(r8), dimension(mgncol), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), dimension(mgncol), intent(in) :: qvi  ! saturation humidity (ice)
  real(r8), dimension(mgncol), intent(in) :: lcldm  ! liquid cloud fraction
  real(r8), dimension(mgncol), intent(in) :: precip_frac ! precipitation fraction (maximum overlap)

  ! fallspeed parameters
  real(r8), dimension(mgncol), intent(in) :: arn  ! rain
  real(r8), dimension(mgncol), intent(in) :: asn  ! snow
  real(r8), dimension(mgncol), intent(in) :: agn  ! graupel
  real(r8),                    intent(in) :: bg 

  ! In-cloud MMRs
  real(r8), dimension(mgncol), intent(in) :: qcic ! cloud liquid
  real(r8), dimension(mgncol), intent(in) :: qiic ! cloud ice
  real(r8), dimension(mgncol), intent(in) :: qric ! rain
  real(r8), dimension(mgncol), intent(in) :: qsic ! snow
  real(r8), dimension(mgncol), intent(in) :: qgic ! graupel

  ! Size parameters
  ! rain
  real(r8), dimension(mgncol), intent(in) :: lamr
  real(r8), dimension(mgncol), intent(in) :: n0r
  ! snow
  real(r8), dimension(mgncol), intent(in) :: lams
  real(r8), dimension(mgncol), intent(in) :: n0s
  ! graupel
  real(r8), dimension(mgncol), intent(in) :: lamg
  real(r8), dimension(mgncol), intent(in) :: n0g

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: pre
  real(r8), dimension(mgncol), intent(out) :: prds
  real(r8), dimension(mgncol), intent(out) :: prdg
  real(r8), dimension(mgncol), intent(out) :: am_evp_st ! Fractional area where rain evaporates.

  real(r8) :: qclr   ! water vapor mixing ratio in clear air
  real(r8) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  real(r8), dimension(mgncol) :: dum

  integer :: i

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise
  am_evp_st = 0._r8
  do i=1,mgncol
     if (qcic(i)+qiic(i) < 1.e-6_r8) then
        dum(i) = 0._r8
     else
        dum(i) = lcldm(i)
     end if
  enddo
  do i=1,mgncol
  ! only calculate if there is some precip fraction > cloud fraction

     if (precip_frac(i) > dum(i)) then

        if (qric(i) >= qsmall .or. qsic(i) >= qsmall .or. qgic(i) >= qsmall) then
           am_evp_st(i) = precip_frac(i) - dum(i)

           ! calculate q for out-of-cloud region
           qclr=(q(i)-dum(i)*qvl(i))/(1._r8-dum(i))
        end if

        ! evaporation of rain
        if (qric(i) >= qsmall) then

           ab = calc_ab(t(i), qvl(i), xxlv)
           eps = 2._r8*pi*n0r(i)*rho(i)*Dv(i)* &
                (f1r/(lamr(i)*lamr(i))+ &
                f2r*(arn(i)*rho(i)/mu(i))**0.5_r8* &
                sc(i)**(1._r8/3._r8)*gamma_half_br_plus5/ &
                (lamr(i)**(5._r8/2._r8+br/2._r8)))

           pre(i) = eps*(qclr-qvl(i))/ab

           ! only evaporate in out-of-cloud region
           ! and distribute across precip_frac
           pre(i)=min(pre(i)*am_evp_st(i),0._r8)
           pre(i)=pre(i)/precip_frac(i)
        else
           pre(i) = 0._r8
        end if

        ! sublimation of snow
        if (qsic(i) >= qsmall) then
           ab = calc_ab(t(i), qvi(i), xxls)
           eps = 2._r8*pi*n0s(i)*rho(i)*Dv(i)* &
                (f1s/(lams(i)*lams(i))+ &
                f2s*(asn(i)*rho(i)/mu(i))**0.5_r8* &
                sc(i)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
                (lams(i)**(5._r8/2._r8+bs/2._r8)))
           prds(i) = eps*(qclr-qvi(i))/ab

          ! only sublimate in out-of-cloud region and distribute over precip_frac
           prds(i)=min(prds(i)*am_evp_st(i),0._r8)
           prds(i)=prds(i)/precip_frac(i)
        else
           prds(i) = 0._r8
        end if

        ! add graupel, do the Same with prdg.
        if (qgic(i).ge.qsmall) then
           ab = calc_ab(t(i), qvi(i), xxls)
           
           eps = 2._r8*pi*n0g(i)*rho(i)*Dv(i)*                    &
                (f1s/(lamg(i)*lamg(i))+                           &
                f2s*(agn(i)*rho(i)/mu(i))**0.5_r8*                &
                sc(i)**(1._r8/3._r8)*gamma(5._r8/2._r8+bg/2._r8)/ &
                (lamg(i)**(5._r8/2._r8+bs/2._r8)))
           prdg(i) = eps*(qclr-qvi(i))/ab
           
           ! only sublimate in out-of-cloud region and distribute over precip_frac
           prdg(i)=min(prdg(i)*am_evp_st(i),0._r8)
           prdg(i)=prdg(i)/precip_frac(i)
        else
           prdg(i) = 0._r8
        end if

     else
        prds(i) = 0._r8
        pre(i) = 0._r8
        prdg(i) = 0._r8
     end if
  enddo

end subroutine evaporate_sublimate_precip_graupel

subroutine evaporate_sublimate_precip_graupel_2D(mgncol, nlev, t, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, precip_frac, arn, asn, agn, bg, qcic, qiic, qric, qsic, qgic, lamr, n0r, lams, n0s, lamg, n0g, &
     pre, prds, prdg, am_evp_st)

  integer,  intent(in) :: mgncol, nlev

  real(r8), dimension(mgncol,nlev), intent(in) :: t    ! temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: rho  ! air density
  real(r8), dimension(mgncol,nlev), intent(in) :: dv   ! water vapor diffusivity
  real(r8), dimension(mgncol,nlev), intent(in) :: mu   ! viscosity
  real(r8), dimension(mgncol,nlev), intent(in) :: sc   ! schmidt number
  real(r8), dimension(mgncol,nlev), intent(in) :: q    ! humidity
  real(r8), dimension(mgncol,nlev), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), dimension(mgncol,nlev), intent(in) :: qvi  ! saturation humidity (ice)
  real(r8), dimension(mgncol,nlev), intent(in) :: lcldm  ! liquid cloud fraction
  real(r8), dimension(mgncol,nlev), intent(in) :: precip_frac ! precipitation fraction (maximum overlap)

  ! fallspeed parameters
  real(r8), dimension(mgncol,nlev), intent(in) :: arn  ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: asn  ! snow
  real(r8), dimension(mgncol,nlev), intent(in) :: agn  ! graupel
  real(r8),                    intent(in) :: bg 

  ! In-cloud MMRs
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic ! cloud liquid
  real(r8), dimension(mgncol,nlev), intent(in) :: qiic ! cloud ice
  real(r8), dimension(mgncol,nlev), intent(in) :: qric ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! snow
  real(r8), dimension(mgncol,nlev), intent(in) :: qgic ! graupel

  ! Size parameters
  ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: lamr
  real(r8), dimension(mgncol,nlev), intent(in) :: n0r
  ! snow
  real(r8), dimension(mgncol,nlev), intent(in) :: lams
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s
  ! graupel
  real(r8), dimension(mgncol,nlev), intent(in) :: lamg
  real(r8), dimension(mgncol,nlev), intent(in) :: n0g

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: pre
  real(r8), dimension(mgncol,nlev), intent(out) :: prds
  real(r8), dimension(mgncol,nlev), intent(out) :: prdg
  real(r8), dimension(mgncol,nlev), intent(out) :: am_evp_st ! Fractional area where rain evaporates.

  real(r8), dimension(mgncol,nlev) :: qclr   ! water vapor mixing ratio in clear air
  real(r8), dimension(mgncol,nlev) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  real(r8), dimension(mgncol,nlev) :: dum

  integer :: i, k
  
  !$acc data create(ab, qclr,dum)

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i=1,mgncol
      am_evp_st(i,k) = 0._r8
    end do
  end do
  
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i=1,mgncol
      if (qcic(i,k)+qiic(i,k) < 1.e-6_r8) then
        dum(i,k) = 0._r8
      else
        dum(i,k) = lcldm(i,k)
      end if
    end do
  end do
  
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i=1,mgncol
      ! only calculate if there is some precip fraction > cloud fraction
      if (precip_frac(i,k) > dum(i,k)) then

        if (qric(i,k) >= qsmall .or. qsic(i,k) >= qsmall .or. qgic(i,k) >= qsmall) then
           am_evp_st(i,k) = precip_frac(i,k) - dum(i,k)

           ! calculate q for out-of-cloud region
           qclr(i,k)=(q(i,k)-dum(i,k)*qvl(i,k))/(1._r8-dum(i,k))
        end if
      end if
    end do
  end do

  ab(:,:) = calc_ab_2D(mgncol, nlev, t(:,:), qvl(:,:), xxlv)

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i=1,mgncol
      ! evaporation of rain
      if (precip_frac(i,k) > dum(i,k) .and. qric(i,k) >= qsmall) then
           eps = 2._r8*pi*n0r(i,k)*rho(i,k)*Dv(i,k)* &
                (f1r/(lamr(i,k)*lamr(i,k))+ &
                f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                sc(i,k)**(1._r8/3._r8)*gamma_half_br_plus5/ &
                (lamr(i,k)**(5._r8/2._r8+br/2._r8)))

           pre(i,k) = eps*(qclr(i,k)-qvl(i,k))/ab(i,k)

           ! only evaporate in out-of-cloud region
           ! and distribute across precip_frac
           pre(i,k)=min(pre(i,k)*am_evp_st(i,k),0._r8)
           pre(i,k)=pre(i,k)/precip_frac(i,k)
      else
         pre(i,k) = 0._r8
      end if
    end do
  end do
  
  ab(:,:) = calc_ab_2D(mgncol, nlev, t(:,:), qvi(:,:), xxls)

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i=1,mgncol
      ! sublimation of snow
      if (precip_frac(i,k) > dum(i,k) .and. qsic(i,k) >= qsmall) then
           eps = 2._r8*pi*n0s(i,k)*rho(i,k)*Dv(i,k)* &
                (f1s/(lams(i,k)*lams(i,k))+ &
                f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                sc(i,k)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
                (lams(i,k)**(5._r8/2._r8+bs/2._r8)))
           prds(i,k) = eps*(qclr(i,k)-qvi(i,k))/ab(i,k)

          ! only sublimate in out-of-cloud region and distribute over precip_frac
           prds(i,k)=min(prds(i,k)*am_evp_st(i,k),0._r8)
           prds(i,k)=prds(i,k)/precip_frac(i,k)
      else
         prds(i,k) = 0._r8
      end if
    end do
  end do

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i=1,mgncol
      ! add graupel, do the Same with prdg.
      if (precip_frac(i,k) > dum(i,k) .and. qgic(i,k).ge.qsmall) then
           eps = 2._r8*pi*n0g(i,k)*rho(i,k)*Dv(i,k)*                    &
                (f1s/(lamg(i,k)*lamg(i,k))+                           &
                f2s*(agn(i,k)*rho(i,k)/mu(i,k))**0.5_r8*                &
                sc(i,k)**(1._r8/3._r8)*gamma(5._r8/2._r8+bg/2._r8)/ &
                (lamg(i,k)**(5._r8/2._r8+bs/2._r8)))
           prdg(i,k) = eps*(qclr(i,k)-qvi(i,k))/ab(i,k)
           
           ! only sublimate in out-of-cloud region and distribute over precip_frac
           prdg(i,k)=min(prdg(i,k)*am_evp_st(i,k),0._r8)
           prdg(i,k)=prdg(i,k)/precip_frac(i,k)
      else
        prdg(i,k) = 0._r8
      end if
    end do
  end do
  
  !$acc end data

end subroutine evaporate_sublimate_precip_graupel_2D


! bergeron process - evaporation of droplets and deposition onto snow
!===================================================================

subroutine bergeron_process_snow(t, rho, dv, mu, sc, qvl, qvi, asn, &
     qcic, qsic, lams, n0s, bergs, mgncol)

  integer, intent(in) :: mgncol

  real(r8), dimension(mgncol), intent(in) :: t    ! temperature
  real(r8), dimension(mgncol), intent(in) :: rho  ! air density
  real(r8), dimension(mgncol), intent(in) :: dv   ! water vapor diffusivity
  real(r8), dimension(mgncol), intent(in) :: mu   ! viscosity
  real(r8), dimension(mgncol), intent(in) :: sc   ! schmidt number
  real(r8), dimension(mgncol), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), dimension(mgncol), intent(in) :: qvi  ! saturation humidity (ice)

  ! fallspeed parameter for snow
  real(r8), dimension(mgncol), intent(in) :: asn

  ! In-cloud MMRs
  real(r8), dimension(mgncol), intent(in) :: qcic ! cloud liquid mixing ratio
  real(r8), dimension(mgncol), intent(in) :: qsic ! snow mixing ratio

  ! Size parameters for snow
  real(r8), dimension(mgncol), intent(in) :: lams
  real(r8), dimension(mgncol), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: bergs

  real(r8) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  integer :: i

  do i=1,mgncol
     if (qsic(i) >= qsmall.and. qcic(i) >= qsmall .and. t(i) < tmelt) then
        ab = calc_ab(t(i), qvi(i), xxls)
        eps = 2._r8*pi*n0s(i)*rho(i)*Dv(i)* &
             (f1s/(lams(i)*lams(i))+ &
             f2s*(asn(i)*rho(i)/mu(i))**0.5_r8* &
             sc(i)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
             (lams(i)**(5._r8/2._r8+bs/2._r8)))
        bergs(i) = eps*(qvl(i)-qvi(i))/ab
     else
        bergs(i) = 0._r8
     end if
  enddo
end subroutine bergeron_process_snow

subroutine bergeron_process_snow_2D(mgncol, nlev, t, rho, dv, mu, sc, qvl, qvi, asn, &
     qcic, qsic, lams, n0s, bergs)

  integer, intent(in) :: mgncol, nlev

  real(r8), dimension(mgncol,nlev), intent(in) :: t    ! temperature
  real(r8), dimension(mgncol,nlev), intent(in) :: rho  ! air density
  real(r8), dimension(mgncol,nlev), intent(in) :: dv   ! water vapor diffusivity
  real(r8), dimension(mgncol,nlev), intent(in) :: mu   ! viscosity
  real(r8), dimension(mgncol,nlev), intent(in) :: sc   ! schmidt number
  real(r8), dimension(mgncol,nlev), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), dimension(mgncol,nlev), intent(in) :: qvi  ! saturation humidity (ice)

  ! fallspeed parameter for snow
  real(r8), dimension(mgncol,nlev), intent(in) :: asn

  ! In-cloud MMRs
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic ! cloud liquid mixing ratio
  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! snow mixing ratio

  ! Size parameters for snow
  real(r8), dimension(mgncol,nlev), intent(in) :: lams
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: bergs

  real(r8), dimension(mgncol,nlev) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  integer :: i, k

  !$acc data create(ab) async(1)
  
  ab(:,:) = calc_ab_2D(mgncol, nlev, t(:,:), qvi(:,:), xxls)

  !$acc parallel loop collapse(2) default(present) async(1) 
  do k = 1, nlev
    do i = 1, mgncol
      if (qsic(i,k) >= qsmall.and. qcic(i,k) >= qsmall .and. t(i,k) < tmelt) then
        
        eps = 2._r8*pi*n0s(i,k)*rho(i,k)*Dv(i,k)* &
             (f1s/(lams(i,k)*lams(i,k))+ &
             f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
             sc(i,k)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
             (lams(i,k)**(5._r8/2._r8+bs/2._r8)))
        bergs(i,k) = eps*(qvl(i,k)-qvi(i,k))/ab(i,k)
      else
        bergs(i,k) = 0._r8
      end if
    end do
  end do
  
  !$acc end data
  
end subroutine bergeron_process_snow_2D

!========================================================================
! Collection of snow by rain to form graupel
!========================================================================

subroutine graupel_collecting_snow(qsic,qric,umr,ums,rho,lamr,n0r,lams,n0s, &
     psacr, mgncol)
  
  integer, intent(in) :: mgncol

  ! In-cloud MMRs
  real(r8), dimension(mgncol), intent(in) :: qsic ! snow
  real(r8), dimension(mgncol), intent(in) :: qric ! rain

  ! mass-weighted fall speeds
  real(r8), dimension(mgncol), intent(in) :: umr ! rain
  real(r8), dimension(mgncol), intent(in) :: ums ! snow
 
  real(r8), dimension(mgncol), intent(in) :: rho  ! air density


  ! Size parameters for rain
  real(r8), dimension(mgncol), intent(in) :: lamr
  real(r8), dimension(mgncol), intent(in) :: n0r

  ! Size parameters for snow
  real(r8), dimension(mgncol), intent(in) :: lams
  real(r8), dimension(mgncol), intent(in) :: n0s

  real(r8), dimension(mgncol), intent(out) :: psacr ! conversion due to coll of snow by rain

  real(r8) :: cons31
  integer :: i

  cons31=pi*pi*ecr*rhosn
	
  do i=1,mgncol

     if (qsic(i).ge.0.1e-3_r8 .and. qric(i).ge.0.1e-3_r8) then
        psacr(i) = cons31*(((1.2_r8*umr(i)-0.95_r8*ums(i))**2+              &
             0.08_r8*ums(i)*umr(i))**0.5_r8*rho(i)*                     &
             n0r(i)*n0s(i)/lams(i)**3*                               &
             (5._r8/(lams(i)**3*lamr(i))+                    &
             2._r8/(lams(i)**2*lamr(i)**2)+                  &
             0.5_r8/(lams(i)*lamr(i)**3)))            
     else
        psacr(i) = 0._r8
     end if

  end do

end subroutine graupel_collecting_snow

subroutine graupel_collecting_snow_2D(mgncol, nlev, qsic,qric,umr,ums, &
                                      rho,lamr,n0r,lams,n0s, psacr)
  
  integer, intent(in) :: mgncol, nlev

  ! In-cloud MMRs
  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! snow
  real(r8), dimension(mgncol,nlev), intent(in) :: qric ! rain

  ! mass-weighted fall speeds
  real(r8), dimension(mgncol,nlev), intent(in) :: umr ! rain
  real(r8), dimension(mgncol,nlev), intent(in) :: ums ! snow
 
  real(r8), dimension(mgncol,nlev), intent(in) :: rho  ! air density


  ! Size parameters for rain
  real(r8), dimension(mgncol,nlev), intent(in) :: lamr
  real(r8), dimension(mgncol,nlev), intent(in) :: n0r

  ! Size parameters for snow
  real(r8), dimension(mgncol,nlev), intent(in) :: lams
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s

  real(r8), dimension(mgncol,nlev), intent(out) :: psacr ! conversion due to coll of snow by rain

  real(r8) :: cons31
  integer :: i, k

  cons31=pi*pi*ecr*rhosn
	
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol

      if (qsic(i,k).ge.0.1e-3_r8 .and. qric(i,k).ge.0.1e-3_r8) then
        psacr(i,k) = cons31*(((1.2_r8*umr(i,k)-0.95_r8*ums(i,k))**2+              &
             0.08_r8*ums(i,k)*umr(i,k))**0.5_r8*rho(i,k)*                     &
             n0r(i,k)*n0s(i,k)/lams(i,k)**3*                               &
             (5._r8/(lams(i,k)**3*lamr(i,k))+                    &
             2._r8/(lams(i,k)**2*lamr(i,k)**2)+                  &
             0.5_r8/(lams(i,k)*lamr(i,k)**3)))            
      else
        psacr(i,k) = 0._r8
      end if

    end do
  end do

end subroutine graupel_collecting_snow_2D

!========================================================================
! Collection of cloud water by graupel
!========================================================================

subroutine graupel_collecting_cld_water(qgic,qcic,ncic,rho,n0g,lamg,bg,agn, &
     psacwg, npsacwg, mgncol)

  integer, intent(in) :: mgncol

  ! In-cloud MMRs
  real(r8), dimension(mgncol), intent(in) :: qgic ! graupel
  real(r8), dimension(mgncol), intent(in) :: qcic ! cloud water

  real(r8), dimension(mgncol), intent(in) :: ncic ! cloud water number conc

  real(r8), dimension(mgncol), intent(in) :: rho  ! air density

  ! Size parameters for graupel
  real(r8), dimension(mgncol), intent(in) :: lamg
  real(r8), dimension(mgncol), intent(in) :: n0g

  ! fallspeed parameters for graupel
  ! Set AGN and BG  as input (in micro_mg3_0.F90)  (select hail or graupel as appropriate)
  real(r8),                    intent(in) :: bg
  real(r8), dimension(mgncol), intent(in) :: agn

  ! Output tendencies
  real(r8), dimension(mgncol), intent(out) :: psacwg
  real(r8), dimension(mgncol), intent(out) :: npsacwg

  real(r8) :: cons
  integer :: i 

  cons = gamma(bg + 3._r8)*pi/4._r8 * ecid

  do i=1,mgncol

        if (qgic(i).ge.1.e-8_r8 .and. qcic(i).ge.qsmall) then

           psacwg(i) = cons*agn(i)*qcic(i)*rho(i)*               &
                  n0g(i)/                        &
                  lamg(i)**(bg+3._r8)
           npsacwg(i) = cons*agn(i)*ncic(i)*rho(i)*              &
                  n0g(i)/                        &
                  lamg(i)**(bg+3._r8)
        else
           psacwg(i)=0._r8
           npsacwg(i)=0._r8
        end if
     enddo
end subroutine graupel_collecting_cld_water

subroutine graupel_collecting_cld_water_2D(mgncol, nlev, qgic, qcic, ncic, rho, &
                                           n0g, lamg, bg, agn, psacwg, npsacwg)

  integer, intent(in) :: mgncol, nlev

  ! In-cloud MMRs
  real(r8), dimension(mgncol,nlev), intent(in) :: qgic ! graupel
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic ! cloud water

  real(r8), dimension(mgncol,nlev), intent(in) :: ncic ! cloud water number conc

  real(r8), dimension(mgncol,nlev), intent(in) :: rho  ! air density

  ! Size parameters for graupel
  real(r8), dimension(mgncol,nlev), intent(in) :: lamg
  real(r8), dimension(mgncol,nlev), intent(in) :: n0g

  ! fallspeed parameters for graupel
  ! Set AGN and BG  as input (in micro_mg3_0.F90)  (select hail or graupel as appropriate)
  real(r8),                    intent(in) :: bg
  real(r8), dimension(mgncol,nlev), intent(in) :: agn

  ! Output tendencies
  real(r8), dimension(mgncol,nlev), intent(out) :: psacwg
  real(r8), dimension(mgncol,nlev), intent(out) :: npsacwg

  real(r8) :: cons
  integer :: i , k

  cons = gamma(bg + 3._r8)*pi/4._r8 * ecid

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol

      if (qgic(i,k).ge.1.e-8_r8 .and. qcic(i,k).ge.qsmall) then

         psacwg(i,k) = cons*agn(i,k)*qcic(i,k)*rho(i,k)*               &
                n0g(i,k)/                        &
                lamg(i,k)**(bg+3._r8)
         npsacwg(i,k) = cons*agn(i,k)*ncic(i,k)*rho(i,k)*              &
                n0g(i,k)/                        &
                lamg(i,k)**(bg+3._r8)
      else
         psacwg(i,k)=0._r8
         npsacwg(i,k)=0._r8
      end if
    end do
  end do
 
end subroutine graupel_collecting_cld_water_2D

!========================================================================
! Conversion of rimed cloud water onto snow to graupel/hail
!========================================================================

subroutine graupel_riming_liquid_snow(psacws,qsic,qcic,nsic,rho,rhosn,rhog,asn,lams,n0s,dtime, &
     pgsacw,nscng,mgncol)

  integer, intent(in) :: mgncol

  ! Accretion of cloud water to snow tendency (modified)
  real(r8), dimension(mgncol), intent(inout) :: psacws

  real(r8), dimension(mgncol), intent(in) :: qsic ! snow mixing ratio
  real(r8), dimension(mgncol), intent(in) :: qcic ! cloud liquid mixing ratio
  real(r8), dimension(mgncol), intent(in) :: nsic ! snow number concentration

  real(r8), dimension(mgncol), intent(in) :: rho   ! air density
  real(r8),                    intent(in) :: rhosn ! snow density
  real(r8),                    intent(in) :: rhog ! graupel density

  real(r8), dimension(mgncol), intent(in) :: asn   ! fall speed parameter for snow

  ! Size parameters for snow
  real(r8), dimension(mgncol), intent(in) :: lams
  real(r8), dimension(mgncol), intent(in) :: n0s

  real(r8),                    intent(in) :: dtime

  !Output process rates
  real(r8), dimension(mgncol), intent(out) :: pgsacw  ! dQ graupel due to collection droplets by snow
  real(r8), dimension(mgncol), intent(out) :: nscng   ! dN graupel due to collection droplets by snow

  real(r8) :: cons
  real(r8) :: rhosu
  real(r8) :: dum
  integer :: i 

!........................................................................
!Input: PSACWS,qs,qc,n0s,rho,lams,rhos,rhog
!Output:PSACWS,PGSACW,NSCNG

  rhosu = 85000._r8/(ra * tmelt)    ! typical air density at 850 mb

  do i=1,mgncol

      cons=4._r8 *2._r8 *3._r8 *rhosu*pi*ecid*ecid*gamma_2bs_plus2/(8._r8*(rhog-rhosn))
      
! Only allow conversion if qni > 0.1 and qc > 0.5 g/kg following Rutledge and Hobbs (1984)
     if (psacws(i).gt.0._r8 .and. qsic(i).GE.0.1e-3_r8 .AND. qcic(i).GE.0.5E-3_r8) then

! portion of riming converted to graupel (Reisner et al. 1998, originally IS1991)
! dtime here is correct.
        pgsacw(i) = min(psacws(i),cons*dtime*n0s(i)*qcic(i)*qcic(i)* &
             asn(i)*asn(i)/ &
             (rho(i)*lams(i)**(2._r8*bs+2._r8))) 

! Mix rat converted into graupel as embryo (Reisner et al. 1998, orig M1990)
        dum= max(rhosn/(rhog-rhosn)*pgsacw(i),0._r8) 

! Number concentraiton of embryo graupel from riming of snow
        nscng(i) = dum/mg0*rho(i)
! Limit max number converted to snow number  (dtime here correct? We think yes)
        nscng(i) = min(nscng(i),nsic(i)/dtime)

! Portion of riming left for snow
        psacws(i) = psacws(i) - pgsacw(i)
     else
        pgsacw(i) = 0._r8
        nscng(i) = 0._r8
     end if
     
  enddo

end subroutine graupel_riming_liquid_snow

subroutine graupel_riming_liquid_snow_2D(mgncol,nlev,psacws,qsic,qcic,nsic,rho, &
                                         rhosn,rhog,asn,lams,n0s,dtime,pgsacw,nscng)

  integer, intent(in) :: mgncol, nlev

  ! Accretion of cloud water to snow tendency (modified)
  real(r8), dimension(mgncol,nlev), intent(inout) :: psacws

  real(r8), dimension(mgncol,nlev), intent(in) :: qsic ! snow mixing ratio
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic ! cloud liquid mixing ratio
  real(r8), dimension(mgncol,nlev), intent(in) :: nsic ! snow number concentration

  real(r8), dimension(mgncol,nlev), intent(in) :: rho   ! air density
  real(r8),                    intent(in) :: rhosn ! snow density
  real(r8),                    intent(in) :: rhog ! graupel density

  real(r8), dimension(mgncol,nlev), intent(in) :: asn   ! fall speed parameter for snow

  ! Size parameters for snow
  real(r8), dimension(mgncol,nlev), intent(in) :: lams
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s

  real(r8),                    intent(in) :: dtime

  !Output process rates
  real(r8), dimension(mgncol,nlev), intent(out) :: pgsacw  ! dQ graupel due to collection droplets by snow
  real(r8), dimension(mgncol,nlev), intent(out) :: nscng   ! dN graupel due to collection droplets by snow

  real(r8) :: cons
  real(r8) :: rhosu
  real(r8) :: dum
  integer :: i, k

!........................................................................
!Input: PSACWS,qs,qc,n0s,rho,lams,rhos,rhog
!Output:PSACWS,PGSACW,NSCNG

  rhosu = 85000._r8/(ra * tmelt)    ! typical air density at 850 mb

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i=1,mgncol

      cons=4._r8 *2._r8 *3._r8 *rhosu*pi*ecid*ecid*gamma_2bs_plus2/(8._r8*(rhog-rhosn))
      
! Only allow conversion if qni > 0.1 and qc > 0.5 g/kg following Rutledge and Hobbs (1984)
      if (psacws(i,k).gt.0._r8 .and. qsic(i,k).GE.0.1e-3_r8 .AND. qcic(i,k).GE.0.5E-3_r8) then

! portion of riming converted to graupel (Reisner et al. 1998, originally IS1991)
! dtime here is correct.
        pgsacw(i,k) = min(psacws(i,k),cons*dtime*n0s(i,k)*qcic(i,k)*qcic(i,k)* &
             asn(i,k)*asn(i,k)/ &
             (rho(i,k)*lams(i,k)**(2._r8*bs+2._r8))) 

! Mix rat converted into graupel as embryo (Reisner et al. 1998, orig M1990)
        dum= max(rhosn/(rhog-rhosn)*pgsacw(i,k),0._r8) 

! Number concentraiton of embryo graupel from riming of snow
        nscng(i,k) = dum/mg0*rho(i,k)
! Limit max number converted to snow number  (dtime here correct? We think yes)
        nscng(i,k) = min(nscng(i,k),nsic(i,k)/dtime)

! Portion of riming left for snow
        psacws(i,k) = psacws(i,k) - pgsacw(i,k)
      else
        pgsacw(i,k) = 0._r8
        nscng(i,k) = 0._r8
      end if
      
    end do
  end do

end subroutine graupel_riming_liquid_snow_2D

!========================================================================
!CHANGE IN Q,N COLLECTION RAIN BY GRAUPEL
!========================================================================

subroutine graupel_collecting_rain(qric,qgic,umg,umr,ung,unr,rho,n0r,lamr,n0g,lamg,&
     pracg,npracg,mgncol)

  integer, intent(in) :: mgncol

  !MMR
  real(r8), dimension(mgncol), intent(in) :: qric  !rain mixing ratio
  real(r8), dimension(mgncol), intent(in) :: qgic  !graupel mixing ratio

  !Mass weighted Fall speeds
  real(r8), dimension(mgncol), intent(in) :: umg ! graupel fall speed
  real(r8), dimension(mgncol), intent(in) :: umr ! rain fall speed

  !Number weighted fall speeds
  real(r8), dimension(mgncol), intent(in) :: ung ! graupel fall speed
  real(r8), dimension(mgncol), intent(in) :: unr ! rain fall speed

  real(r8), dimension(mgncol), intent(in) :: rho   ! air density

 ! Size parameters for rain
  real(r8), dimension(mgncol), intent(in) :: n0r
  real(r8), dimension(mgncol), intent(in) :: lamr

 ! Size parameters for graupel
  real(r8), dimension(mgncol), intent(in) :: n0g
  real(r8), dimension(mgncol), intent(in) :: lamg


  !Output process rates
  real(r8), dimension(mgncol), intent(out) :: pracg   ! Q collection rain by graupel
  real(r8), dimension(mgncol), intent(out) :: npracg  ! N collection rain by graupel

! Add collection of graupel by rain above freezing
! assume all rain collection by graupel above freezing is shed
! assume shed drops are 1 mm in size

  integer :: i 
  real(r8) :: cons41
  real(r8) :: cons32
  real(r8) :: dum

  cons41=pi*pi*ecr*rhow
  cons32=pi/2._r8*ecr

  do i=1,mgncol

     if (qric(i).ge.1.e-8_r8.and.qgic(i).ge.1.e-8_r8) then

! pracg is mixing ratio of rain per sec collected by graupel/hail
        pracg(i) = cons41*(((1.2_r8*umr(i)-0.95_r8*umg(i))**2._r8+                   &
             0.08_r8*umg(i)*umr(i))**0.5_r8*rho(i)*                      &
             n0r(i)*n0g(i)/lamr(i)**3._r8*                              &
             (5._r8/(lamr(i)**3._r8*lamg(i))+                    &
             2._r8/(lamr(i)**2._r8*lamg(i)**2._r8)+				   &
             0.5_r8/(lamr(i)*lamg(i)**3._r8)))

! assume 1 mm drops are shed, get number shed per sec

        dum = pracg(i)/5.2e-7_r8
        
        npracg(i) = cons32*rho(i)*(1.7_r8*(unr(i)-ung(i))**2._r8+            &
             0.3_r8*unr(i)*ung(i))**0.5_r8*n0r(i)*n0g(i)*              &
             (1._r8/(lamr(i)**3._r8*lamg(i))+                      &
             1._r8/(lamr(i)**2._r8*lamg(i)**2._r8)+                   &
             1._r8/(lamr(i)*lamg(i)**3._r8))
        
! hm 7/15/13, remove limit so that the number of collected drops can smaller than 
! number of shed drops
!            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
        npracg(i)=npracg(i)- dum
     else
        npracg(i)=0._r8
        pracg(i)=0._r8
     end if
     
  enddo

end subroutine graupel_collecting_rain

subroutine graupel_collecting_rain_2D(mgncol,nlev,qric,qgic,umg,umr,ung,unr,rho, &
                                      n0r,lamr,n0g,lamg,pracg,npracg)

  integer, intent(in) :: mgncol, nlev

  !MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: qric  !rain mixing ratio
  real(r8), dimension(mgncol,nlev), intent(in) :: qgic  !graupel mixing ratio

  !Mass weighted Fall speeds
  real(r8), dimension(mgncol,nlev), intent(in) :: umg ! graupel fall speed
  real(r8), dimension(mgncol,nlev), intent(in) :: umr ! rain fall speed

  !Number weighted fall speeds
  real(r8), dimension(mgncol,nlev), intent(in) :: ung ! graupel fall speed
  real(r8), dimension(mgncol,nlev), intent(in) :: unr ! rain fall speed

  real(r8), dimension(mgncol,nlev), intent(in) :: rho   ! air density

 ! Size parameters for rain
  real(r8), dimension(mgncol,nlev), intent(in) :: n0r
  real(r8), dimension(mgncol,nlev), intent(in) :: lamr

 ! Size parameters for graupel
  real(r8), dimension(mgncol,nlev), intent(in) :: n0g
  real(r8), dimension(mgncol,nlev), intent(in) :: lamg


  !Output process rates
  real(r8), dimension(mgncol,nlev), intent(out) :: pracg   ! Q collection rain by graupel
  real(r8), dimension(mgncol,nlev), intent(out) :: npracg  ! N collection rain by graupel

! Add collection of graupel by rain above freezing
! assume all rain collection by graupel above freezing is shed
! assume shed drops are 1 mm in size

  integer :: i, k
  real(r8) :: cons41
  real(r8) :: cons32
  real(r8) :: dum

  cons41=pi*pi*ecr*rhow
  cons32=pi/2._r8*ecr

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i=1,mgncol

      if (qric(i,k).ge.1.e-8_r8.and.qgic(i,k).ge.1.e-8_r8) then

! pracg is mixing ratio of rain per sec collected by graupel/hail
        pracg(i,k) = cons41*(((1.2_r8*umr(i,k)-0.95_r8*umg(i,k))**2._r8+                   &
             0.08_r8*umg(i,k)*umr(i,k))**0.5_r8*rho(i,k)*                      &
             n0r(i,k)*n0g(i,k)/lamr(i,k)**3._r8*                              &
             (5._r8/(lamr(i,k)**3._r8*lamg(i,k))+                    &
             2._r8/(lamr(i,k)**2._r8*lamg(i,k)**2._r8)+				   &
             0.5_r8/(lamr(i,k)*lamg(i,k)**3._r8)))

! assume 1 mm drops are shed, get number shed per sec

        dum = pracg(i,k)/5.2e-7_r8
        
        npracg(i,k) = cons32*rho(i,k)*(1.7_r8*(unr(i,k)-ung(i,k))**2._r8+            &
             0.3_r8*unr(i,k)*ung(i,k))**0.5_r8*n0r(i,k)*n0g(i,k)*              &
             (1._r8/(lamr(i,k)**3._r8*lamg(i,k))+                      &
             1._r8/(lamr(i,k)**2._r8*lamg(i,k)**2._r8)+                   &
             1._r8/(lamr(i,k)*lamg(i,k)**3._r8))
        
! hm 7/15/13, remove limit so that the number of collected drops can smaller than 
! number of shed drops
!            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
        npracg(i,k)=npracg(i,k)- dum
      else
        npracg(i,k)=0._r8
        pracg(i,k)=0._r8
      end if
      
    end do
  end do

end subroutine graupel_collecting_rain_2D

!========================================================================
! Rain riming snow to graupel
!========================================================================
! Conversion of rimed rainwater onto snow converted to graupel

subroutine graupel_rain_riming_snow(pracs,npracs,psacr,qsic,qric,nric,nsic,n0s,lams,n0r,lamr,dtime,&
     pgracs,ngracs,mgncol)

  integer, intent(in) :: mgncol
  
  ! Accretion of rain by snow
  real(r8), dimension(mgncol), intent(inout) :: pracs
  real(r8), dimension(mgncol), intent(inout) :: npracs
  real(r8), dimension(mgncol), intent(inout) :: psacr  ! conversion due to coll of snow by rain

  !MMR
  real(r8), dimension(mgncol), intent(in) :: qsic  !snow mixing ratio
  real(r8), dimension(mgncol), intent(in) :: qric  !rain mixing ratio

  real(r8), dimension(mgncol), intent(in) :: nric ! rain number concentration
  real(r8), dimension(mgncol), intent(in) :: nsic ! snow number concentration

 ! Size parameters for snow
  real(r8), dimension(mgncol), intent(in) :: n0s
  real(r8), dimension(mgncol), intent(in) :: lams

 ! Size parameters for rain
  real(r8), dimension(mgncol), intent(in) :: n0r
  real(r8), dimension(mgncol), intent(in) :: lamr

  real(r8),                    intent(in) :: dtime

  !Output process rates
  real(r8), dimension(mgncol), intent(out) :: pgracs  ! Q graupel due to collection rain by snow
  real(r8), dimension(mgncol), intent(out) :: ngracs  ! N graupel due to collection rain by snow

!Input: PRACS,NPRACS,PSACR,qni,qr,lams,lamr,nr,ns  Note: No PSACR in MG2
!Output:PGRACS,NGRACS,PRACS,PSACR

  integer :: i 
  real(r8) :: cons18
  real(r8) :: cons19
  real(r8) :: dum,fmult

  cons18=rhosn*rhosn
  cons19=rhow*rhow

  do i=1,mgncol
     
     fmult=0._r8

! only allow conversion if qs > 0.1 and qr > 0.1 g/kg following rutledge and hobbs (1984)
     if (pracs(i).gt.0._r8.and.qsic(i).ge.0.1e-3_r8.and.qric(i).ge.0.1e-3_r8) then

          ! portion of collected rainwater converted to graupel (reisner et al. 1998)
          dum = cons18*(4._r8/lams(i))**3*(4._r8/lams(i))**3 &    
               /(cons18*(4._r8/lams(i))**3*(4._r8/lams(i))**3+ &  
               cons19*(4._r8/lamr(i))**3*(4._r8/lamr(i))**3)
          dum=min(dum,1._r8)
          dum=max(dum,0._r8)
          pgracs(i) = (1._r8-dum)*pracs(i)
          ngracs(i) = (1._r8-dum)*npracs(i)
          ! limit max number converted to min of either rain or snow number concentration
          ngracs(i) = min(ngracs(i),nric(i)/dtime)
          ngracs(i) = min(ngracs(i),nsic(i)/dtime)
          
          ! amount left for snow production
          pracs(i) = pracs(i) - pgracs(i)
          npracs(i) = npracs(i) - ngracs(i)
          
          ! conversion to graupel due to collection of snow by rain
          psacr(i)=psacr(i)*(1._r8-dum)
       else
          pgracs(i) = 0._r8
          ngracs(i) = 0._r8
        end if
     enddo 

end subroutine graupel_rain_riming_snow

subroutine graupel_rain_riming_snow_2D(mgncol,nlev,pracs,npracs,psacr,qsic,qric,nric,nsic,n0s, &
                                       lams,n0r,lamr,dtime, pgracs,ngracs)

  integer, intent(in) :: mgncol, nlev
  
  ! Accretion of rain by snow
  real(r8), dimension(mgncol,nlev), intent(inout) :: pracs
  real(r8), dimension(mgncol,nlev), intent(inout) :: npracs
  real(r8), dimension(mgncol,nlev), intent(inout) :: psacr  ! conversion due to coll of snow by rain

  !MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: qsic  !snow mixing ratio
  real(r8), dimension(mgncol,nlev), intent(in) :: qric  !rain mixing ratio

  real(r8), dimension(mgncol,nlev), intent(in) :: nric ! rain number concentration
  real(r8), dimension(mgncol,nlev), intent(in) :: nsic ! snow number concentration

 ! Size parameters for snow
  real(r8), dimension(mgncol,nlev), intent(in) :: n0s
  real(r8), dimension(mgncol,nlev), intent(in) :: lams

 ! Size parameters for rain
  real(r8), dimension(mgncol,nlev), intent(in) :: n0r
  real(r8), dimension(mgncol,nlev), intent(in) :: lamr

  real(r8),                         intent(in) :: dtime

  !Output process rates
  real(r8), dimension(mgncol,nlev), intent(out) :: pgracs  ! Q graupel due to collection rain by snow
  real(r8), dimension(mgncol,nlev), intent(out) :: ngracs  ! N graupel due to collection rain by snow

!Input: PRACS,NPRACS,PSACR,qni,qr,lams,lamr,nr,ns  Note: No PSACR in MG2
!Output:PGRACS,NGRACS,PRACS,PSACR

  integer :: i, k
  real(r8) :: cons18
  real(r8) :: cons19
  real(r8) :: dum,fmult

  cons18=rhosn*rhosn
  cons19=rhow*rhow

  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol
     
      fmult=0._r8

! only allow conversion if qs > 0.1 and qr > 0.1 g/kg following rutledge and hobbs (1984)
      if (pracs(i,k).gt.0._r8.and.qsic(i,k).ge.0.1e-3_r8.and.qric(i,k).ge.0.1e-3_r8) then

          ! portion of collected rainwater converted to graupel (reisner et al. 1998)
          dum = cons18*(4._r8/lams(i,k))**3*(4._r8/lams(i,k))**3 &    
               /(cons18*(4._r8/lams(i,k))**3*(4._r8/lams(i,k))**3+ &  
               cons19*(4._r8/lamr(i,k))**3*(4._r8/lamr(i,k))**3)
          dum=min(dum,1._r8)
          dum=max(dum,0._r8)
          pgracs(i,k) = (1._r8-dum)*pracs(i,k)
          ngracs(i,k) = (1._r8-dum)*npracs(i,k)
          ! limit max number converted to min of either rain or snow number concentration
          ngracs(i,k) = min(ngracs(i,k),nric(i,k)/dtime)
          ngracs(i,k) = min(ngracs(i,k),nsic(i,k)/dtime)
          
          ! amount left for snow production
          pracs(i,k) = pracs(i,k) - pgracs(i,k)
          npracs(i,k) = npracs(i,k) - ngracs(i,k)
          
          ! conversion to graupel due to collection of snow by rain
          psacr(i,k)=psacr(i,k)*(1._r8-dum)
      else
          pgracs(i,k) = 0._r8
          ngracs(i,k) = 0._r8
      end if
      
    end do
  end do 

end subroutine graupel_rain_riming_snow_2D

!========================================================================
! Rime Splintering
!========================================================================
subroutine graupel_rime_splintering(t,qcic,qric,qgic,psacwg,pracg,&
     qmultg,nmultg,qmultrg,nmultrg,mgncol)

  integer, intent(in) :: mgncol
  
  real(r8), dimension(mgncol), intent(in) :: t  !temperature

  !MMR
  real(r8), dimension(mgncol), intent(in) :: qcic  !liquid mixing ratio
  real(r8), dimension(mgncol), intent(in) :: qric  !rain mixing ratio
  real(r8), dimension(mgncol), intent(in) :: qgic  !graupel mixing ratio

  ! Already calculated terms for collection 
  real(r8), dimension(mgncol), intent(inout) :: psacwg ! collection droplets by graupel
  real(r8), dimension(mgncol), intent(inout) :: pracg  ! collection rain by graupel

  !Output process rates for splintering
  real(r8), dimension(mgncol), intent(out) :: qmultg  ! Q ice mult of droplets/graupel
  real(r8), dimension(mgncol), intent(out) :: nmultg  ! N ice mult of droplets/graupel
  real(r8), dimension(mgncol), intent(out) :: qmultrg  ! Q ice mult of rain/graupel
  real(r8), dimension(mgncol), intent(out) :: nmultrg  ! N ice mult of rain/graupel


!Input: qg,qc,qr, PSACWG,PRACG,T
!Output: NMULTG,QMULTG,NMULTRG,QMULTRG,PSACWG,PRACG

  integer :: i 
  real(r8) :: fmult
  real(r8) :: tm_3,tm_5,tm_8

  tm_3 = tmelt - 3._r8
  tm_5 = tmelt - 5._r8
  tm_8 = tmelt - 8._r8


!nmultg,qmultg                                                                             .
!========================================================================
  do i=1,mgncol

     nmultrg(i)=0._r8
     qmultrg(i)=0._r8
     nmultg(i)=0._r8
     qmultg(i)=0._r8

     if (qgic(i).ge.0.1e-3_r8) then
        if (qcic(i).ge.0.5e-3_r8.or.qric(i).ge.0.1e-3_r8) then
           if (psacwg(i).gt.0._r8.or.pracg(i).gt.0._r8) then
              if (t(i).lt.tm_3 .and. t(i).gt.tm_8) then
                 if (t(i).gt.tm_3) then
                    fmult = 0._r8
                 else if (t(i).le.tm_3.and.t(i).gt.tm_5)  then
                    fmult = (tm_3-t(i))/2._r8
                 else if (t(i).ge.tm_8.and.t(i).le.tm_5)   then
                    fmult = (t(i)-tm_8)/3._r8
                 else if (t(i).lt.tm_8) then
                    fmult = 0._r8
                 end if

! 1000 is to convert from kg to g  
! splintering from droplets accreted onto graupel

                 if (psacwg(i).gt.0._r8) then
                    nmultg(i) = 35.e4_r8*psacwg(i)*fmult*1000._r8
                    qmultg(i) = nmultg(i)*mmult
                    
! constrain so that transfer of mass from graupel to ice cannot be more mass
! than was rimed onto graupel
                    
                    qmultg(i) = min(qmultg(i),psacwg(i))
                    psacwg(i) = psacwg(i)-qmultg(i)
                    
                 end if


!nmultrg,qmultrg                                                                             .
!========================================================================

! riming and splintering from accreted raindrops

                 if (pracg(i).gt.0._r8) then
                    nmultrg(i) = 35.e4_r8*pracg(i)*fmult*1000._r8
                    qmultrg(i) = nmultrg(i)*mmult
                    
! constrain so that transfer of mass from graupel to ice cannot be more mass
! than was rimed onto graupel

                    qmultrg(i) = min(qmultrg(i),pracg(i))
                    pracg(i) = pracg(i)-qmultrg(i)

                 end if

              end if
           end if
        end if
     end if
  enddo

end subroutine graupel_rime_splintering

subroutine graupel_rime_splintering_2D(mgncol,nlev,t,qcic,qric,qgic,psacwg,pracg,&
     qmultg,nmultg,qmultrg,nmultrg)

  integer, intent(in) :: mgncol, nlev
  
  real(r8), dimension(mgncol,nlev), intent(in) :: t  !temperature

  !MMR
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic  !liquid mixing ratio
  real(r8), dimension(mgncol,nlev), intent(in) :: qric  !rain mixing ratio
  real(r8), dimension(mgncol,nlev), intent(in) :: qgic  !graupel mixing ratio

  ! Already calculated terms for collection 
  real(r8), dimension(mgncol,nlev), intent(inout) :: psacwg ! collection droplets by graupel
  real(r8), dimension(mgncol,nlev), intent(inout) :: pracg  ! collection rain by graupel

  !Output process rates for splintering
  real(r8), dimension(mgncol,nlev), intent(out) :: qmultg  ! Q ice mult of droplets/graupel
  real(r8), dimension(mgncol,nlev), intent(out) :: nmultg  ! N ice mult of droplets/graupel
  real(r8), dimension(mgncol,nlev), intent(out) :: qmultrg  ! Q ice mult of rain/graupel
  real(r8), dimension(mgncol,nlev), intent(out) :: nmultrg  ! N ice mult of rain/graupel


!Input: qg,qc,qr, PSACWG,PRACG,T
!Output: NMULTG,QMULTG,NMULTRG,QMULTRG,PSACWG,PRACG

  integer :: i, k
  real(r8) :: fmult
  real(r8) :: tm_3,tm_5,tm_8

  tm_3 = tmelt - 3._r8
  tm_5 = tmelt - 5._r8
  tm_8 = tmelt - 8._r8


!nmultg,qmultg                                                                             .
!========================================================================
  !$acc parallel loop collapse(2) default(present) async(1)
  do k = 1, nlev
    do i = 1, mgncol

      nmultrg(i,k)=0._r8
      qmultrg(i,k)=0._r8
      nmultg(i,k)=0._r8
      qmultg(i,k)=0._r8

      if (qgic(i,k).ge.0.1e-3_r8) then
        if (qcic(i,k).ge.0.5e-3_r8.or.qric(i,k).ge.0.1e-3_r8) then
           if (psacwg(i,k).gt.0._r8.or.pracg(i,k).gt.0._r8) then
              if (t(i,k).lt.tm_3 .and. t(i,k).gt.tm_8) then
                 if (t(i,k).gt.tm_3) then
                    fmult = 0._r8
                 else if (t(i,k).le.tm_3.and.t(i,k).gt.tm_5)  then
                    fmult = (tm_3-t(i,k))/2._r8
                 else if (t(i,k).ge.tm_8.and.t(i,k).le.tm_5)   then
                    fmult = (t(i,k)-tm_8)/3._r8
                 else if (t(i,k).lt.tm_8) then
                    fmult = 0._r8
                 end if

! 1000 is to convert from kg to g  
! splintering from droplets accreted onto graupel

                 if (psacwg(i,k).gt.0._r8) then
                    nmultg(i,k) = 35.e4_r8*psacwg(i,k)*fmult*1000._r8
                    qmultg(i,k) = nmultg(i,k)*mmult
                    
! constrain so that transfer of mass from graupel to ice cannot be more mass
! than was rimed onto graupel
                    
                    qmultg(i,k) = min(qmultg(i,k),psacwg(i,k))
                    psacwg(i,k) = psacwg(i,k)-qmultg(i,k)
                    
                 end if


!nmultrg,qmultrg                                                                             .
!========================================================================

! riming and splintering from accreted raindrops

                 if (pracg(i,k).gt.0._r8) then
                    nmultrg(i,k) = 35.e4_r8*pracg(i,k)*fmult*1000._r8
                    qmultrg(i,k) = nmultrg(i,k)*mmult
                    
! constrain so that transfer of mass from graupel to ice cannot be more mass
! than was rimed onto graupel

                    qmultrg(i,k) = min(qmultrg(i,k),pracg(i,k))
                    pracg(i,k) = pracg(i,k)-qmultrg(i,k)

                 end if

              end if
           end if
        end if
      end if
    end do
  end do

end subroutine graupel_rime_splintering_2D

!========================================================================
!UTILITIES
!========================================================================

pure function no_limiter()
  real(r8) :: no_limiter

  no_limiter = transfer(limiter_off, no_limiter)

end function no_limiter

pure function limiter_is_on(lim)
  real(r8), intent(in) :: lim
  logical :: limiter_is_on

  limiter_is_on = transfer(lim, limiter_off) /= limiter_off

end function limiter_is_on

end module micro_mg_utils
