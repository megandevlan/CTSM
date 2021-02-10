module CLUBBmomentsType


#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! CLUBB moments data structure
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use PatchType      , only : patch
  !use AnnualFluxDribbler, only : annual_flux_dribbler_type,annual_flux_dribbler_gridcell
  !
  implicit none
  save
  private
  !
  type, public :: clubbmoments_type

     ! Moments
     real(r8), pointer :: wp2_grid   (:)   ! Vertical velocity varaince (m**2/s**2)
     real(r8), pointer :: thlp2_grid (:)   ! Temperature variance (K**2)

   contains
     procedure, public  :: Init            ! Public initialization method
     procedure, private :: InitAllocate    ! initialize/allocate
     procedure, private :: InitHistory     ! setup history fields
     procedure, private :: InitCold        ! initialize for cold start
     !procedure, public  :: Restart         ! setup restart fields

  end type clubbmoments_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

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
    !logical           , intent(in) :: is_simple_buildtemp        ! If using simple building temp method
    !logical           , intent(in) :: is_prog_buildtemp          ! If using prognostic building temp method

    !SHR_ASSERT_ALL_FL((ubound(t_grnd_col) == (/bounds%endc/)), sourcefile,__LINE__)

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

    allocate( this%wp2_grid      (begg:endg))          ; this%wp2_grid        (:)   = nan
    allocate( this%thlp2_grid    (begg:endg))          ; this%thlp2_grid      (:)   = nan

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
         avgflag='A', long_name='Surface vertical velocity variance for CLUBB', &
         ptr_gcell=this%wp2_grid, default='inactive')

    this%thlp2_grid(begg:endg) = spval
    call hist_addfld1d (fname='THLP2_CLUBB',  units='K^2',  &
         avgflag='A', long_name='Surface temperature variance for CLUBB',&
         ptr_gcell=this%thlp2_grid, default='inactive')

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
    this%thlp2_grid(bounds%begg:bounds%endg)       = 0.0_r8

  end subroutine InitCold

  !------------------------------------------------------------------------


end module CLUBBmomentsType
