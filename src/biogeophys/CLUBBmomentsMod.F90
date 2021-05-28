module CLUBBmomentsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module computes the surface higher order moments that the CLUBB
  ! scheme in CAM needs. These are computed accounting for land surface 
  ! heterogeneity, so all moments are first computed at the patch level 
  ! before scaling up to gridcell levels, according to Machulskaya and
  ! Mironov (2018), which is based on the original Andre et al. 
  ! 1978 equations. 
  ! 
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use clm_varctl             , only : iulog
  use clm_varctl             , only : compute_CLUBB_HMG, compute_CLUBB_HTG 
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varcon             , only : spval
  use decompMod              , only : bounds_type
  !
  use FrictionvelocityMod    , only : frictionvel_type
  use EnergyFluxType         , only : energyflux_type
  use atm2lndType            , only : atm2lnd_type
  use lnd2atmType            , only : lnd2atm_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use TemperatureType        , only : temperature_type
  use clm_varcon             , only : cpair,grav,hvap,vkc
  !
  use subgridAveMod          , only : p2g
  use GridcellType           , only : grc
  use LandunitType           , only : lun
  use ColumnType             , only : col
  use PatchType              , only : patch
  ! 
  ! use histFileMod            , only : hist_addfld1d
  
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CLUBBmoments            ! Calculate moments 
  !
  type, public :: clubbmoments_type

     ! Moments
     real(r8), pointer :: wp2_grid   (:)   ! Vertical velocity varaince (m**2/s**2)
     real(r8), pointer :: thlp2_grid (:)   ! Temperature variance (K**2)
     real(r8), pointer :: qp2_grid   (:)   ! Specific humidity variance (kg**2/kg**2)
     real(r8), pointer :: wpthlp_grid(:)   ! Covariance of vertical velocity and temperature (m/s * K)
     real(r8), pointer :: wpqp_grid(:)     ! Covariance of vertical velocity and humidity (m/s * kg/kg)
     real(r8), pointer :: thlpqp_grid(:)   ! Covariance of temperature and humidity (K * kg/kg)
     real(r8), pointer :: wp3_grid(:)      ! Third order moment of vertical velocity (m**3/s**3)
     real(r8), pointer :: up2_grid(:)      ! Horizontal velocity variance (m**2/s**2)
     real(r8), pointer :: wp4_grid(:)      ! Fourth order moment of vertical velocity (m**4/s**4)
     real(r8), pointer :: wp2thlp_grid(:)  ! Covariance of vertical velocity variance and temperature (m**2/s**2 * K)
     real(r8), pointer :: wp2qp_grid(:)    ! Covariance of vertical velocity variance and humidity (m**2/s**2 * kg/kg)
     real(r8), pointer :: wpqp2_grid(:)    ! Covariance of vertical velocity and humidity variance (m/s * kg**2/kg**2)
     real(r8), pointer :: wpthlp2_grid(:)  ! Covariance of vertical velocity and temperature variance (m/s * K**2)
     real(r8), pointer :: wpthlpqp_grid(:) ! Product of surface vertical velocity, humidity, and temperature (m/s kg/kg K)
     real(r8), pointer :: upwp_grid(:)     ! Zonal momentum flux 
     real(r8), pointer :: vpwp_grid(:)     ! Meridional momentum flux 

   contains
     procedure, public  :: Init            ! Public initialization method
     procedure, private :: InitAllocate    ! initialize/allocate
     procedure, private :: InitHistory     ! setup history fields
     procedure, private :: InitCold        ! initialize for cold start
     !procedure, public  :: Restart         ! setup restart fields

  end type clubbmoments_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    !    Allocate and initialize the data type and setup history, and initialize
    !    for cold-start.
    ! !USES:
    implicit none
    ! !ARGUMENTS:
    class(clubbmoments_type)       :: this
    type(bounds_type) , intent(in) :: bounds
    !real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )
    !logical           , intent(in) :: is_simple_buildtemp        ! If using
    !simple building temp method
    !logical           , intent(in) :: is_prog_buildtemp          ! If using
    !prognostic building temp method

    !SHR_ASSERT_ALL_FL((ubound(t_grnd_col) == (/bounds%endc/)),
    !sourcefile,__LINE__)

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds )

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak
    implicit none
    !
    ! !ARGUMENTS:
    class(clubbmoments_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !integer :: begp, endp
    !integer :: begc, endc
    !integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    !begp = bounds%begp; endp= bounds%endp
    !begc = bounds%begc; endc= bounds%endc
    !begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate( this%wp2_grid      (begg:endg))          ; this%wp2_grid(:)     = nan
    allocate( this%thlp2_grid    (begg:endg))          ; this%thlp2_grid(:)   = nan
    allocate( this%qp2_grid      (begg:endg))          ; this%qp2_grid(:)     = nan
    allocate( this%wpthlp_grid   (begg:endg))          ; this%wpthlp_grid(:)  = nan
    allocate( this%wpqp_grid     (begg:endg))          ; this%wpqp_grid(:)    = nan
    allocate( this%thlpqp_grid   (begg:endg))          ; this%thlpqp_grid(:)  = nan
    allocate( this%wp3_grid      (begg:endg))          ; this%wp3_grid(:)     = nan
    allocate( this%up2_grid      (begg:endg))          ; this%up2_grid(:)     = nan
    allocate( this%wp4_grid      (begg:endg))          ; this%wp4_grid(:)     = nan
    allocate( this%wp2thlp_grid  (begg:endg))          ; this%wp2thlp_grid(:) = nan
    allocate( this%wp2qp_grid    (begg:endg))          ; this%wp2qp_grid(:)   = nan
    allocate( this%wpqp2_grid    (begg:endg))          ; this%wpqp2_grid(:)   = nan
    allocate( this%wpthlp2_grid  (begg:endg))          ; this%wpthlp2_grid(:) = nan
    allocate( this%wpthlpqp_grid (begg:endg))          ; this%wpthlpqp_grid(:) = nan
    allocate( this%upwp_grid     (begg:endg))          ; this%upwp_grid(:)    = nan
    allocate( this%vpwp_grid     (begg:endg))          ; this%vpwp_grid(:)    = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Setup fields that can be output to history files
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !use clm_varpar     , only : nlevsno, nlevgrnd
    !use clm_varctl     , only : use_cn, use_hydrstress
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    use ncdio_pio      , only : ncd_inqvdlen
    implicit none
    !
    ! !ARGUMENTS:
    class(clubbmoments_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !integer           :: begp, endp
    !integer           :: begc, endc
    !integer           :: begl, endl
    integer           :: begg, endg
    integer           :: dimlen
    integer           :: err_code
    logical           :: do_io
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    !begp = bounds%begp; endp= bounds%endp
    !begc = bounds%begc; endc= bounds%endc
    !begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg


    this%wp2_grid(begg:endg) = spval
    call hist_addfld1d (fname='WP2_CLUBB',  units='m^2/s^2',  &
         avgflag='A', long_name='Surface vertical velocity variance for CLUBB',&
         ptr_gcell=this%wp2_grid, default='inactive')

    this%thlp2_grid(begg:endg) = spval
    call hist_addfld1d (fname='THLP2_CLUBB',  units='K^2',  &
         avgflag='A', long_name='Surface temperature variance for CLUBB',&
         ptr_gcell=this%thlp2_grid, default='inactive')

    this%qp2_grid(begg:endg) = spval
    call hist_addfld1d (fname='QP2_CLUBB',  units='kg^2/kg^2',  &
         avgflag='A', long_name='Surface specific humidity variance for CLUBB',&
         ptr_gcell=this%qp2_grid, default='inactive')

    this%wpthlp_grid(begg:endg) = spval
    call hist_addfld1d (fname='WPTHLP_CLUBB',  units='m/s K',  &
         avgflag='A', long_name='Surface covariance of vertical velocity and temperature for CLUBB',&
         ptr_gcell=this%wpthlp_grid, default='inactive')

    this%wpqp_grid(begg:endg) = spval
    call hist_addfld1d (fname='WPQP_CLUBB',  units='m/s kg/kg',  &
         avgflag='A', long_name='Surface covariance of vertical velocity and humidity for CLUBB',&
         ptr_gcell=this%wpqp_grid, default='inactive')

    this%thlpqp_grid(begg:endg) = spval
    call hist_addfld1d (fname='THLPQP_CLUBB',  units='kg/kg K',  &
         avgflag='A', long_name='Surface covariance of temperature and humidity for CLUBB',&
         ptr_gcell=this%thlpqp_grid, default='inactive')

    this%wp3_grid(begg:endg) = spval
    call hist_addfld1d (fname='WP3_CLUBB',  units='m^3/s^3',  &
         avgflag='A', long_name='Surface third order moment of vertical velocity for CLUBB',&
         ptr_gcell=this%wp3_grid, default='inactive')

    this%up2_grid(begg:endg) = spval
    call hist_addfld1d (fname='UP2_CLUBB',  units='m^2/s^2',  &
         avgflag='A', long_name='Surface horizontal velocity variance for CLUBB',&
         ptr_gcell=this%up2_grid, default='inactive')

    this%wp4_grid(begg:endg) = spval
    call hist_addfld1d (fname='WP4_CLUBB',  units='m^4/s^4',  &
         avgflag='A', long_name='Surface fourth order moment of vertical velocity for CLUBB',&
         ptr_gcell=this%wp4_grid, default='inactive')

    this%wp2thlp_grid(begg:endg) = spval
    call hist_addfld1d (fname='WP2THLP_CLUBB',  units='m^2/s^2 K',  &
         avgflag='A', long_name='Surface covariance of vertical velocity variance and temperature for CLUBB',&
         ptr_gcell=this%wp2thlp_grid, default='inactive')

   this%wp2qp_grid(begg:endg) = spval
    call hist_addfld1d (fname='WP2QP_CLUBB',  units='m^2/s^2 kg/kg',  &
         avgflag='A', long_name='Surface covariance of vertical velocity variance and humidity for CLUBB',&
         ptr_gcell=this%wp2qp_grid, default='inactive')

   this%wpqp2_grid(begg:endg) = spval
    call hist_addfld1d (fname='WPQP2_CLUBB',  units='m/s kg^2/kg^2',  &
         avgflag='A', long_name='Surface covariance of vertical velocity and humidity variance for CLUBB',&
         ptr_gcell=this%wpqp2_grid, default='inactive')

   this%wpthlp2_grid(begg:endg) = spval
    call hist_addfld1d (fname='WPTHLP2_CLUBB',  units='m/s K^2',  &
         avgflag='A', long_name='Surface covariance of vertical velocity and temperature variance for CLUBB',&
         ptr_gcell=this%wpthlp2_grid, default='inactive')

   this%wpthlpqp_grid(begg:endg) = spval
    call hist_addfld1d (fname='WPTHLPQP_CLUBB',  units='m/s K kg/kg',  &
         avgflag='A', long_name='Surface product of vertical velocity, humidity, and temperature for CLUBB',&
         ptr_gcell=this%wpthlpqp_grid, default='inactive')

    this%upwp_grid(begg:endg) = spval
    call hist_addfld1d (fname='UPWP_CLUBB',  units='m^2/s^2',  &
         avgflag='A', long_name='Surface zonal momentum flux for CLUBB',&
         ptr_gcell=this%upwp_grid, default='inactive')

    this%vpwp_grid(begg:endg) = spval
    call hist_addfld1d (fname='VPWP_CLUBB',  units='m^2/s^2',  &
         avgflag='A', long_name='Surface meridional momentum flux for CLUBB',&
         ptr_gcell=this%vpwp_grid, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    !use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(clubbmoments_type), intent(in) :: this
    type(bounds_type)       , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !integer :: c,l
    !-----------------------------------------------------------------------

    this%wp2_grid(bounds%begg:bounds%endg)       = 0.0_r8
    this%thlp2_grid(bounds%begg:bounds%endg)     = 0.0_r8
    this%qp2_grid(bounds%begg:bounds%endg)       = 0.0_r8
    this%wpthlp_grid(bounds%begg:bounds%endg)    = 0.0_r8
    this%wpqp_grid(bounds%begg:bounds%endg)      = 0.0_r8
    this%thlpqp_grid(bounds%begg:bounds%endg)    = 0.0_r8
    this%wp3_grid(bounds%begg:bounds%endg)       = 0.0_r8
    this%up2_grid(bounds%begg:bounds%endg)       = 0.0_r8
    this%wp4_grid(bounds%begg:bounds%endg)       = 0.0_r8
    this%wp2thlp_grid(bounds%begg:bounds%endg)   = 0.0_r8
    this%wp2qp_grid(bounds%begg:bounds%endg)     = 0.0_r8
    this%wpqp2_grid(bounds%begg:bounds%endg)     = 0.0_r8
    this%wpthlp2_grid(bounds%begg:bounds%endg)   = 0.0_r8
    this%wpthlpqp_grid(bounds%begg:bounds%endg)  = 0.0_r8
    this%upwp_grid(bounds%begg:bounds%endg)      = 0.0_r8
    this%vpwp_grid(bounds%begg:bounds%endg)      = 0.0_r8

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine CLUBBmoments(bounds,clubbmoments_inst,frictionvel_inst,&
       energyflux_inst,atm2lnd_inst,lnd2atm_inst,&
       waterdiagnosticbulk_inst,temperature_inst)
    !
    ! !DESCRIPTION:
    !
    ! Compute moments that CLUBB uses at the surface:
    !    Start with wp2 (vertical variance) 
    !    Next up is thlp2 (temeprature variance).
    !
    ! ARGUMENTS:
    !class(clubbmoments_type)               , intent(inout)         :: this
    type(bounds_type)                      , intent(in)            :: bounds
    type(clubbmoments_type)                , intent(inout)         :: clubbmoments_inst
    type(frictionvel_type)                 , intent(in)            :: frictionvel_inst
    type(energyflux_type)                  , intent(in)            :: energyflux_inst
    type(atm2lnd_type)                     , intent(in)            :: atm2lnd_inst
    type(lnd2atm_type)                     , intent(in)            :: lnd2atm_inst
    type(waterdiagnosticbulk_type)         , intent(in)            :: waterdiagnosticbulk_inst
    type(temperature_type)                 , intent(in)            :: temperature_inst
    !  !LOCAL VARIABLES:
    integer  :: p                                        ! Patch indices
    integer  :: c,g,l                                    ! Column/gridcell/landunit indices
    real(r8) :: rhoair                                   ! Air density
    real(r8) :: sfcP                                     ! Surface atm pressure (Pa) 
    real(r8) :: theta                                    ! Potential temperature (K) 
    real(r8) :: eflx_sh_pr_conversionPatch               ! Patch level eflx_sh correction term 
    real(r8) :: eflx_sh_il_conversionPatch               ! Patch level eflx_sh correction term
    real(r8) :: eflx_dynbal_grcPatch                     ! Patch level eflx_sh correction term
    real(r8) :: eflx_sh_patchCorrected                   ! Total patch eflx_sh, including correction terms
    real(r8) :: Sk                                       ! Simplified budget ignoring <w'w'T'> (accoring to Nate Chaney's script)
    real(r8) :: c8                                       ! Closure constant from Cheng et al. (2005) 
    real(r8) :: p1                                       ! Selected to match limit in free convection (according to Nate Chaney's script)
    real(r8) :: phiuu,phivv,phiww                        ! Stability corrections for velocity variances 
    real(r8) :: phieps                                   ! Stability correction for TKE dissipation rate
    real(r8) :: phitau                                   ! Stability correction for relaxation time scale
    real(r8) :: wbar,swt,sw1,sw2,Skw,a,w1,w2             ! Variables needed to compute wp4

    real(r8) :: KH(bounds%begp:bounds%endp)              ! Kinematic heat flux 
    real(r8) :: KL(bounds%begp:bounds%endp)              ! Humidity flux
    real(r8) :: wstar(bounds%begp:bounds%endp)           ! Convective velocity scale [m/s]
    real(r8) :: wp2(bounds%begp:bounds%endp)             ! Vertical velocityvariance (m2/s2)
    real(r8) :: thlp2(bounds%begp:bounds%endp)           ! Temperature variance (K2)
    real(r8) :: upwp(bounds%begp:bounds%endp)            ! Zonal momentum flux (m2/s2)
    real(r8) :: vpwp(bounds%begp:bounds%endp)            ! Meridional momentum flux (m2/s2)
    real(r8) :: qp2(bounds%begp:bounds%endp)             ! Specific humidity variance (kg2/kg2)
    real(r8) :: thlpqp(bounds%begp:bounds%endp)          ! Covariance of temperature and humidity (K * kg/kg)
    real(r8) :: Tv(bounds%begp:bounds%endp)              ! Virtual temperature at the surface (column=patch level)
    real(r8) :: thlp2_grid(bounds%begg:bounds%endg)      ! Gridcell weighted mean of THLP2
    real(r8) :: theta_grid(bounds%begg:bounds%endg)      ! Gridcell weighted mean of THETA (Tv)
    real(r8) :: theta_sqr(bounds%begp:bounds%endp)       ! Patch level difference in THETA from grid mean 
    real(r8) :: theta_sqr_grid(bounds%begg:bounds%endg)  ! Gridcell weighted mean of patch THETA differences 
    real(r8) :: qp2_grid(bounds%begg:bounds%endg)        ! Gridcell weighted mean of QP2
    real(r8) :: q_grid(bounds%begg:bounds%endg)          ! Gridcell weighted mean of Q
    real(r8) :: q_sqr(bounds%begp:bounds%endp)           ! Patch level difference in Q from grid mean 
    real(r8) :: q_sqr_grid(bounds%begg:bounds%endg)      ! Gridcell weighted mean of patch Q differences 
    real(r8) :: thlpqp_grid(bounds%begg:bounds%endg)     ! Gridcell weighted mean of thlpqp
    real(r8) :: thlpqp_sqr(bounds%begp:bounds%endp)      ! Patch level difference in theta and Q from grid mean
    real(r8) :: thlpqp_sqr_grid(bounds%begg:bounds%endg) ! Gridcell weighted mean of patch theta and Q differences
    real(r8) :: wp3(bounds%begp:bounds%endp)             ! Vertical velocity skew (m3/s3)
    real(r8) :: up2(bounds%begp:bounds%endp)             ! Horizontal velocityvariance (m2/s2)
    real(r8) :: wp4(bounds%begp:bounds%endp)             ! Fourth order moment of vertical velocity (m4/s4) 
    real(r8) :: wp2thlp(bounds%begp:bounds%endp)         ! Covariance of vertical velocity variance and temperature (m2/s2 * K) 
    real(r8) :: wp2qp(bounds%begp:bounds%endp)           ! Covariance of vertical velocity variance and humidity (m2/s2 * kg/kg)
    real(r8) :: wpqp2(bounds%begp:bounds%endp)           ! Covariance of vertical velocity and humidity variance (m/s * kg2/kg2)
    real(r8) :: wpqp2_grid(bounds%begg:bounds%endg)      ! Gridcell weighted mean of wpqp2 
    real(r8) :: wpthlp2(bounds%begp:bounds%endp)         ! Covariance of vertical velocity and temperature variance (m/s * kg2/kg2)
    real(r8) :: wpthlp2_grid(bounds%begg:bounds%endg)    ! Gridcell weighted mean of wpthlp2
    real(r8) :: wpthlpqp_term1(bounds%begp:bounds%endp)         ! wpthlpqp
    real(r8) :: wpthlpqp_term1_grid(bounds%begg:bounds%endg)    ! Gridcell weighted mean of wpthlpqp_term1
    real(r8) :: wpthlpqp_term2(bounds%begp:bounds%endp)         ! wpthlpqp
    real(r8) :: wpthlpqp_term2_grid(bounds%begg:bounds%endg)    ! Gridcell weighted mean of wpthlpqp_term2

    !------------------------------------------------------------------------

    associate(             &
         ustar                 => frictionvel_inst%ustar_patch                 , &  ! Input: [real(r8) (:)   ]  friction velocity [m/s]  
         !wstar                 => frictionvel_inst%wstar_patch                 , &  ! Input: [real(r8) (:)   ]  convective velocity scale [m/s]
         zeta                  => frictionvel_inst%zeta_patch                  , &  ! Input: [real(r8) (:)   ]  dimensionless stability parameter 
         forc_rho              => atm2lnd_inst%forc_rho_downscaled_col         , &  ! Input: [real(r8) (:)   ]  density (kg/m**3)   
         forc_pbot             => atm2lnd_inst%forc_pbot_downscaled_col        , &  ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         q_ref2m               => waterdiagnosticbulk_inst%q_ref2m_patch       , &  ! Input: [real(r8) (:)   ]  2 m height surface specific humidity (kg/g)
         t_ref2m               => temperature_inst%t_ref2m_patch               , &  ! Input: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)
         eflx_sh_tot           => energyflux_inst%eflx_sh_tot_patch            , &  ! Input: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_pr_conversion => energyflux_inst%eflx_sh_precip_conversion_col, &  ! Input: [real(r8) (:)   ]  correction to eflx_sh
         eflx_sh_il_conversion => lnd2atm_inst%eflx_sh_ice_to_liq_col          , &  ! Input: [real(r8) (:)   ]  correction to eflx_sh  
         eflx_dynbal_grc       => energyflux_inst%eflx_dynbal_grc              , &  ! Input: [real(r8) (:)   ]  correction eflx_sh 
         eflx_lh_tot           => energyflux_inst%eflx_lh_tot_patch            , &  ! Input: [real(r8) (:)   ]  total latent heat flux (W/m**2)  [+ to atm]
         taux                  => energyflux_inst%taux_patch                   , &  ! Input: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)                                  
         tauy                  => energyflux_inst%tauy_patch                     &  ! Input: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2) 
         ) 


       do p = bounds%begp,bounds%endp
          if (patch%active(p)) then

               ! write(iulog,*)'MDF: patch type... ',patch%itype(p)
               ! write(iulog,*)'MDF: patch weight: ',patch%wtgcell(p)

               ! Compute the variance of vertical velocity for each patch 
               ! --------------------------------------------------------
               if (zeta(p) >= 0._r8) then  
                  wp2(p)=1.75_r8*(ustar(p)**2.0_r8)
               else
                  wp2(p)=(ustar(p)**2.0_r8)*(1.75_r8+(2.0_r8*(-zeta(p))**(2.0_r8/3.0_r8)))
               end if  

               ! Compute the variance of potential temperature for each patch
               ! ------------------------------------------------------------ 
               c = patch%column(p) 
               g = patch%gridcell(p)
               l = patch%landunit(p)    ! Only defining l for writing out log file -- DELETE later
   
               ! write(iulog,*)'MDF: col type...   ',col%itype(c)
               ! write(iulog,*)'MDF: landunit type... ',lun%itype(l)
 
               ! Use column or gridcell values for each patch 
               rhoair = forc_rho(c)    
               sfcP   = forc_pbot(c)

               ! Compute zonal and meridional momentum fluxes 
               upwp(p) = taux(p)/rhoair 
               vpwp(p) = tauy(p)/rhoair

               ! Correct total SHFLX as in lnd2atmMod
               eflx_sh_pr_conversionPatch = eflx_sh_pr_conversion(c) 
               eflx_sh_il_conversionPatch = eflx_sh_il_conversion(c)
               eflx_dynbal_grcPatch       = eflx_dynbal_grc(g) 
               eflx_sh_patchCorrected = eflx_sh_tot(p) + eflx_sh_pr_conversionPatch + eflx_sh_il_conversionPatch - eflx_dynbal_grcPatch
              
               ! Compute kinematic heat flux based on eflx_sh 
               KH(p) = eflx_sh_patchCorrected / rhoair / cpair 

               if (zeta(p)>=0.0_r8) then 
                  thlp2(p) = (KH(p)**2.0_r8/ustar(p)**2.0_r8)*4.0_r8
               else 
                  thlp2(p) = (KH(p)**2.0_r8/ustar(p)**2.0_r8)*(4.0_r8*((1.0_r8 - 8.3_r8*zeta(p))**(-2.0_r8/3.0_r8)) )
               end if 

               ! Also compute virtual temperature, which we'll use as 'theta' 
               theta       = t_ref2m(p)*(100000.0_r8/sfcP)**0.286_r8
               !Tv(p)       = (1.0_r8 + 0.61_r8*q_ref2m(p))*t_ref2m(p)
               Tv(p)       = (1.0_r8 + 0.61_r8*q_ref2m(p))*theta

               ! Compute the variance of total water specific humidity
               ! -----------------------------------------------------

               ! Compute humidity flux based on eflx_lh_tot
               KL(p) = eflx_lh_tot(p) / rhoair / hvap
               
               ! Compute qp2 for patch
               if (zeta(p)>=0.0_r8) then
                  qp2(p) = (KL(p)**2.0_r8/ustar(p)**2.0_r8)*4.0_r8
               else
                  qp2(p) = (KL(p)**2.0_r8/ustar(p)**2.0_r8)*(4.0_r8*((1.0_r8 -8.3_r8*zeta(p))**(-2.0_r8/3.0_r8)) )
               end if

               ! Compute the covariance of temperature and humidity 
               ! -------------------------------------------------
               thlpqp(p) = (thlp2(p)**0.5_r8) * (qp2(p)**0.5_r8)

               ! Compute the third order moment of vertical velocity 
               !   See Cheng et al. (2005) and Hsieh and Katul (1997) 
               !   for more information on this formaultion. 
               ! ---------------------------------------------------
               c8 = 5.0_r8 
               p1 = 4.0_r8 * 5.0_r8

               if (zeta(p)<0.0_r8) then    !Only valid for convective ABL (zeta<0)
                  phiuu = (2.7_r8  * (1.0_r8-3.0_r8*zeta(p))**(1.0_r8/3.0_r8) )**2.0_r8
                  phivv = (2.1_r8  * (1.0_r8-3.0_r8*zeta(p))**(1.0_r8/3.0_r8) )**2.0_r8
                  phiww = (1.25_r8 * (1.0_r8-3.0_r8*zeta(p))**(1.0_r8/3.0_r8) )**2.0_r8
                  phieps = vkc * (10.0_r8 + 7.5_r8*(-zeta(p)) + 6.25_r8*zeta(p)**2.0_r8 )/(4.0_r8 + 2.5_r8*(-zeta(p)))
                  phitau = (phiuu+phivv+phiww)/phieps
   
                  Sk = ( (3.0_r8/(c8+p1))*1.25_r8*vkc*zeta(p) / ((1.0_r8-3.0_r8*zeta(p))**(2.0_r8/3.0_r8)) )*phitau
               else
                  Sk = 0.0_r8
               end if

               wp3(p) = Sk * (wp2(p)**(3.0_r8/2.0_r8))

               ! Compute the variance of horizontal velocity 
               ! -------------------------------------------
                if (KH(p)>=0.0_r8) then
                  wstar(p) = ((grav/Tv(p))*KH(p))**(1.0_r8/3.0_r8) 
                  up2(p)   = (4.0_r8 * ustar(p)**2.0_r8) + (0.3_r8 * wstar(p)**2.0_r8)
               else 
                  up2(p) = 4.0_r8 * ustar(p)**2.0_r8 
               end if 

               ! Compute the fourth order moment of vertical velocity
               ! ---------------------------------------------------
               wbar = 0.0_r8       ! MOST constraint 
               swt  = 0.4 

               sw1 = swt*wp2(p)**0.5_r8
               sw2 = swt*wp2(p)**0.5_r8
               Skw = wp3(p)/(wp2(p)**(3.0_r8/2.0_r8))
               ! QC from Nate Chaney:
               if (Skw<-0.8_r8) then 
                  Skw=-0.8
               end if  

               a = 0.5_r8*(1.0_r8 - Skw*(1.0_r8/(4.0_r8*(1-swt**2.0_r8)**3.0_r8) + Skw**2.0_r8)**0.5_r8)            
               w1 = wp2(p)**0.5_r8*((1.0_r8-a)/a)**0.5_r8*(1.0_r8 - swt**2.0_r8)**0.5_r8 + wbar
               w2 = -wp2(p)**0.5_r8*((1.0_r8-a)/a)**0.5_r8*(1.0_r8 - swt**2.0_r8)**0.5_r8 + wbar

               wp4(p) = a*((w1 - wbar)**4.0_r8 + 6.0_r8*(w1 - wbar)**2.0_r8*sw1**2.0_r8 + 3.0_r8*sw1**4.0_r8) + (1.0_r8-a)*((w2 - wbar)**4.0_r8 + 6.0_r8*(w2 - wbar)**2.0_r8*sw2**2.0_r8 + 3.0_r8*sw2**4.0_r8)

         end if 
       end do      

    ! Get weighted means over the gridcell 
    ! Note: Some of these moments don't require any special treatment 
    !       in that the grid mean is computed stragith from p2g. 
    !       But others represent SGS heterogeneity in way that needs 
    !       this grid mean as well as patch deviations from it.
    !       See Machulskaya and Mironov (2018) for details of approach. 
    ! -----------------
    call p2g(bounds, &
         wp2(bounds%begp:bounds%endp), & 
         clubbmoments_inst%wp2_grid (bounds%begg:bounds%endg), & 
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')    

    call p2g(bounds, &
         wp3(bounds%begp:bounds%endp), &
         clubbmoments_inst%wp3_grid (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')

    call p2g(bounds, &
         up2(bounds%begp:bounds%endp), &
         clubbmoments_inst%up2_grid (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')

    call p2g(bounds, &
         upwp(bounds%begp:bounds%endp), &
         clubbmoments_inst%upwp_grid (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')

    call p2g(bounds, &
         vpwp(bounds%begp:bounds%endp), &
         clubbmoments_inst%vpwp_grid (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')

    call p2g(bounds, &
         wp4(bounds%begp:bounds%endp), &
         clubbmoments_inst%wp4_grid (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')

    call p2g(bounds, &
         KH(bounds%begp:bounds%endp), &
         clubbmoments_inst%wpthlp_grid (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')

    call p2g(bounds, &
         KL(bounds%begp:bounds%endp), &
         clubbmoments_inst%wpqp_grid (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')

    if (compute_CLUBB_HMG) then 
 
       call  p2g(bounds, &
             thlp2(bounds%begp:bounds%endp), &
             clubbmoments_inst%thlp2_grid (bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             qp2(bounds%begp:bounds%endp), &
             clubbmoments_inst%qp2_grid (bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             thlpqp(bounds%begp:bounds%endp), &
             clubbmoments_inst%thlpqp_grid (bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')
        
       do g = bounds%begg,bounds%endg
          clubbmoments_inst%wp2thlp_grid(g)  = 0.0_r8
          clubbmoments_inst%wp2qp_grid(g)    = 0.0_r8
          clubbmoments_inst%wpqp2_grid(g)    = 0.0_r8
          clubbmoments_inst%wpthlp2_grid(g)  = 0.0_r8 
          clubbmoments_inst%wpthlpqp_grid(g) = 0.0_r8
       end do 

    else if (compute_CLUBB_HTG) then

       call  p2g(bounds, & 
             thlp2(bounds%begp:bounds%endp), &
             thlp2_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             Tv(bounds%begp:bounds%endp), &
             theta_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             qp2(bounds%begp:bounds%endp), &
             qp2_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             q_ref2m(bounds%begp:bounds%endp), &
             q_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             thlpqp(bounds%begp:bounds%endp), &
             thlpqp_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       ! Some moments need a second term to represent SGS heterogeneity
       !   (The difference in patch values from gridcell mean) 
       do p = bounds%begp,bounds%endp
             if (patch%active(p)) then
                g = patch%gridcell(p)
             
                ! Compute squared difference between patch and grid mean temperature 
                theta_sqr(p) = (Tv(p) - theta_grid(g))**2.0_r8

                ! Compute squared diff between patch and grid mean humidity
                q_sqr(p) = (q_ref2m(p) - q_grid(g))**2.0_r8

                ! Compute diff between patch and grid mean humidity and temperature 
                thlpqp_sqr(p) = (Tv(p) - theta_grid(g)) * (q_ref2m(p) - q_grid(g))
    
                ! Compute patch level covariances between variances and other variables
                wp2thlp(p) = (Tv(p) - theta_grid(g))*clubbmoments_inst%wp2_grid(g)
                wp2qp(p)   = (q_ref2m(p) - q_grid(g))*clubbmoments_inst%qp2_grid(g)
                wpqp2(p)   = (q_ref2m(p) - q_grid(g))*KL(p)                
                wpthlp2(p) = (Tv(p) - theta_grid(g))*KH(p)
           
                ! Compute terms that go into means of wpthlpqp
                wpthlpqp_term1(p) = (Tv(p) - theta_grid(g))*KL(p)
                wpthlpqp_term2(p) = (q_ref2m(p) - q_grid(g))*KH(p)
 
             end if
       end do 

       call  p2g(bounds, &
             wp2thlp(bounds%begp:bounds%endp), &
             clubbmoments_inst%wp2thlp_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             wp2qp(bounds%begp:bounds%endp), &
             clubbmoments_inst%wp2qp_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             theta_sqr(bounds%begp:bounds%endp), &
             theta_sqr_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             q_sqr(bounds%begp:bounds%endp), &
             q_sqr_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             thlpqp_sqr(bounds%begp:bounds%endp), &
             thlpqp_sqr_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             wpqp2(bounds%begp:bounds%endp), &
             wpqp2_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             wpthlp2(bounds%begp:bounds%endp), &
             wpthlp2_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             wpthlpqp_term1(bounds%begp:bounds%endp), &
             wpthlpqp_term1_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       call  p2g(bounds, &
             wpthlpqp_term2(bounds%begp:bounds%endp), &
             wpthlpqp_term2_grid(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

       do g = bounds%begg,bounds%endg
          clubbmoments_inst%thlp2_grid(g)  = thlp2_grid(g)  + theta_sqr_grid(g)
          clubbmoments_inst%qp2_grid(g)    = qp2_grid(g)    + q_sqr_grid(g)
          clubbmoments_inst%thlpqp_grid(g) = thlpqp_grid(g) + thlpqp_sqr_grid(g)
          clubbmoments_inst%wpqp2_grid(g)  = 2.0_r8*wpqp2_grid(g) 
          clubbmoments_inst%wpthlp2_grid(g)  = 2.0_r8*wpthlp2_grid(g)
          clubbmoments_inst%wpthlpqp_grid(g) = wpthlpqp_term1_grid(g) + wpthlpqp_term2_grid(g)
       end do
    end if 
    end associate

  end subroutine CLUBBmoments

end module CLUBBmomentsMod
