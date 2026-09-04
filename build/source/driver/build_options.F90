module build_options

  implicit none
  private

#ifdef NGEN_ACTIVE
  logical, parameter, public :: ngen_active = .true.
#else
  logical, parameter, public :: ngen_active = .false.
#endif

end module build_options
