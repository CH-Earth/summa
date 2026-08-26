module build_options

  implicit none
  private

#ifdef NGEN_ACTIVE
  logical, parameter, public :: ngen_active = .true.
#else
  logical, parameter, public :: ngen_active = .false.
#endif

#ifdef NGEN_FORCING_ACTIVE
  logical, parameter, public :: ngen_forcing_active = .true.
#else
  logical, parameter, public :: ngen_forcing_active = .false.
#endif

#ifdef MIZUROUTE_ACTIVE
  logical, parameter, public :: mizuroute_active = .true.
#else
  logical, parameter, public :: mizuroute_active = .false.
#endif

end module build_options
