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

#ifdef NGEN_OUTPUT_ACTIVE
  logical, parameter, public :: ngen_output_active = .true.
#else
  logical, parameter, public :: ngen_output_active = .false.
#endif

#ifdef SUNDIALS_ACTIVE
  logical, parameter, public :: sundials_active = .true.
#else
  logical, parameter, public :: sundials_active = .false.
#endif

#ifdef OPENWQ_ACTIVE
  logical, parameter, public :: openwq_active = .true.
#else
  logical, parameter, public :: openwq_active = .false.
#endif

#ifdef MIZUROUTE_ACTIVE
  logical, parameter, public :: mizuroute_active = .true.
#else
  logical, parameter, public :: mizuroute_active = .false.
#endif

#ifdef ACTORS_ACTIVE
  logical, parameter, public :: actors_active = .true.
#else
  logical, parameter, public :: actors_active = .false.
#endif

end module build_options
