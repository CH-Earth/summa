!-----------------------------------------------------------------------
! Compatibility shim for mizuRoute.
! Provides the NetCDF data-type identifiers required by popMetadat
! without introducing a dependency on the PIO library.
!-----------------------------------------------------------------------
module pio_utils

  use nrtype, only: i4b
  use netcdf, only: nf90_float, nf90_double, nf90_int

  implicit none
  private

  public :: ncd_float
  public :: ncd_double
  public :: ncd_int

  integer(i4b), parameter :: ncd_float  = nf90_float
  integer(i4b), parameter :: ncd_double = nf90_double
  integer(i4b), parameter :: ncd_int    = nf90_int

end module pio_utils
