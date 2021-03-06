module test_mkncdio
! Module for testing mkncdio

  use mkncdio
  use test_mod

  implicit none
  private

  public :: test_get_dim_lengths
  public :: test_get_nonexisting_var

  character(len=*), parameter :: modname = 'test_mkncdio'
  
contains
  
!------------------------------------------------------------------------------
  subroutine test_get_dim_lengths

    implicit none
    
    character(len=128) :: testname
    integer            :: ncid
    character(len=128) :: varname
    integer            :: ndims, ndims_t
    integer            :: dim_lengths(nf_max_var_dims), dim_lengths_t(nf_max_var_dims)

    character(len=*), parameter :: filename = 'unit_testers/inputs/test_lookup_2d_netcdf.nc'

    character(len=*), parameter :: subname = 'test_get_dim_lengths'

    ! Open netcdf file that will be used for most tests
    call check_ret(nf_open(filename, 0, ncid), subname)

    testname = '3d variable'
    varname = 'lookup3d'
    ndims_t = 3
    dim_lengths_t = 0
    dim_lengths_t(1) = 2
    dim_lengths_t(2) = 3
    dim_lengths_t(3) = 4
    call get_dim_lengths(ncid, varname, ndims, dim_lengths)
    call check_results

    call check_ret(nf_close(ncid), subname)

  contains
    subroutine check_results
      call test_is(ndims, ndims_t, modname//' -- '//subname//' -- '//trim(testname)//' -- ndims')
      call test_is(dim_lengths(1:ndims), dim_lengths_t(1:ndims_t), &
                   modname//' -- '//subname//' -- '//trim(testname)//' -- dim_lengths')
    end subroutine check_results

  end subroutine test_get_dim_lengths

!------------------------------------------------------------------------------
  subroutine test_get_nonexisting_var

    implicit none
    
    character(len=128) :: testname
    integer            :: ncid
    character(len=128) :: varname
    integer            :: varid
    logical            :: varexists

    character(len=*), parameter :: filename = 'unit_testers/inputs/test_lookup_2d_netcdf.nc'

    character(len=*), parameter :: subname = 'test_get_nonexiting_var'

    testname = 'check if variables exist'
    varname = 'lookup3d'
    ! Open netcdf file that will be used for most tests
    call check_ret(nf_open(filename, 0, ncid), subname)
    call check_ret(nf_inq_varid(ncid, "zztop", varid), subname, varexists=varexists)
    call test_is(.not.varexists, modname//' -- '//subname//' -- '//trim(testname)//' -- non existing var')
    call check_ret(nf_inq_varid(ncid, varname, varid), subname, varexists=varexists)
    call test_is(varexists, modname//' -- '//subname//' -- '//trim(testname)//' -- existing var')

  end subroutine test_get_nonexisting_var

end module test_mkncdio
