module parameter_search

  USE nr_type, only: i4b, rkind, lgt

  implicit none
  private


  ! **************************************************************************************************
  ! Information describing one model parameter.
  !
  ! The model-specific layer populates this structure. Parameters with sampled=.true. are included
  ! in the search vector. Parameters with sampled=.false. may still participate in parameter
  ! dependencies using their scalar trial value.
  !
  ! The current implementation assumes that all parameters represented by this structure are
  ! spatially uniform. Spatially varying parameter fields are not currently supported.
  ! **************************************************************************************************

  type, public :: parameter_info

    character(len=64)  :: name
    character(len=128) :: long_name = ''
    character(len=32)  :: units     = '-'

    real(rkind)        :: trial_value
    real(rkind)        :: lower
    real(rkind)        :: upper

    logical(lgt)       :: sampled = .false.

    character(len=16)  :: transformation = 'none'

  end type parameter_info


  ! **************************************************************************************************
  ! Ordered parameter dependency.
  !
  ! param_index contains indices into parameter_spec%params. Adjacent parameters must be separated
  ! by at least gap_fraction times the complete range from the lower bound of the first parameter
  ! to the upper bound of the last parameter.
  ! **************************************************************************************************

  type, public :: ordered_constraint

    integer(i4b), allocatable :: param_index(:)

    real(rkind)               :: gap_fraction = 0._rkind

  end type ordered_constraint


  ! **************************************************************************************************
  ! Model-provided specification of the parameter search problem.
  !
  ! Any model can use the parameter-search routines by populating this structure with scalar trial
  ! values, bounds, sampled flags, transformations, and optional ordered dependencies.
  ! **************************************************************************************************

  type, public :: parameter_spec

    type(parameter_info),     allocatable :: params(:)
    type(ordered_constraint), allocatable :: ordered(:)

  end type parameter_spec


  ! **************************************************************************************************
  ! Resolved information for one ordered constraint.
  ! **************************************************************************************************

  type :: ordered_search_info

    integer(i4b), allocatable :: param_index(:)

    real(rkind)               :: gap

  end type ordered_search_info


  ! **************************************************************************************************
  ! Resolved information used during parameter sampling.
  ! **************************************************************************************************

  type, public :: parameter_search_info

    type(parameter_info), allocatable :: params(:)

    ! Master parameter index -> search-vector index; zero for non-sampled parameters
    integer(i4b), allocatable :: search_index(:)

    character(len=64), allocatable :: param_names(:)
    character(len=16), allocatable :: transformation(:)

    ! Physical/model-space bounds
    real(rkind), allocatable :: lower(:)
    real(rkind), allocatable :: upper(:)

    ! Bounds in transformed search coordinates
    real(rkind), allocatable :: search_lower(:)
    real(rkind), allocatable :: search_upper(:)

    type(ordered_search_info), allocatable :: ordered(:)

  end type parameter_search_info


  public :: initialize_parameter_search
  public :: sample_parameters
  public :: perturb_parameters
  public :: perturb_parameters_dds

contains


  ! **************************************************************************************************
  ! Initialize the parameter search.
  !
  ! Validates the model-provided parameter specification, constructs the sampled parameter vector,
  ! resolves parameter transformations and ordered dependencies, and verifies that the constrained
  ! parameter region is feasible.
  ! **************************************************************************************************

  subroutine initialize_parameter_search(spec, search, err, message)

    implicit none

    type(parameter_spec),        intent(in)  :: spec
    type(parameter_search_info), intent(out) :: search
    integer(i4b),                intent(out) :: err
    character(*),                intent(out) :: message

    character(len=1024) :: cmessage

    err = 0
    message = 'initialize_parameter_search/'

    call validate_parameter_spec(spec,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    call build_search_info(spec,search,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    call initialize_constraints(spec,search,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

  end subroutine initialize_parameter_search


  ! **************************************************************************************************
  ! Validate the model-provided parameter specification.
  ! **************************************************************************************************

  subroutine validate_parameter_spec(spec, err, message)

    implicit none

    type(parameter_spec), intent(in)  :: spec
    integer(i4b),         intent(out) :: err
    character(*),         intent(out) :: message

    integer(i4b) :: i
    integer(i4b) :: j
    integer(i4b) :: nSampled

    err = 0
    message = 'validate_parameter_spec/'

    ! parameter specification is required
    if(.not.allocated(spec%params))then
      message=trim(message)//'parameter specification is not defined'
      err=20; return
    endif

    if(size(spec%params) == 0)then
      message=trim(message)//'parameter specification is empty'
      err=20; return
    endif

    nSampled = 0

    do i=1,size(spec%params)

      ! parameter name is required
      if(len_trim(spec%params(i)%name) == 0)then
        message=trim(message)//'parameter name is empty'
        err=20; return
      endif

      ! duplicate parameter names are not permitted
      do j=1,i-1

        if(trim(spec%params(i)%name) == trim(spec%params(j)%name))then
          message=trim(message)//'duplicate parameter: '//trim(spec%params(i)%name)
          err=20; return
        endif

      enddo

      ! valid physical parameter bounds are required
      if(spec%params(i)%lower > spec%params(i)%upper)then
        message=trim(message)//'invalid bounds for parameter: '//trim(spec%params(i)%name)
        err=20; return
      endif

      ! trial values of non-sampled parameters must lie within their physical bounds
      if(.not.spec%params(i)%sampled)then

        if(spec%params(i)%trial_value < spec%params(i)%lower .or. &
           spec%params(i)%trial_value > spec%params(i)%upper)then

          message=trim(message)//'trial value outside bounds for parameter: '// &
                  trim(spec%params(i)%name)
          err=20; return

        endif

      endif

      ! validate parameter transformation
      select case(trim(spec%params(i)%transformation))

        case ('none')

          ! no additional requirements

        case ('log')

          if(spec%params(i)%sampled)then

            if(spec%params(i)%lower <= 0._rkind .or. &
               spec%params(i)%upper <= 0._rkind)then

              message=trim(message)//'log transformation requires positive bounds for parameter: '// &
                      trim(spec%params(i)%name)
              err=20; return

            endif

          endif

        case default

          message=trim(message)//'unsupported parameter transformation: '// &
                  trim(spec%params(i)%transformation)
          err=20; return

      end select

      if(spec%params(i)%sampled) nSampled = nSampled + 1

    enddo

    if(nSampled == 0)then
      message=trim(message)//'no parameters are included in the parameter search'
      err=20; return
    endif

  end subroutine validate_parameter_spec


  ! **************************************************************************************************
  ! Construct the sampled parameter vector and master-to-search index mapping.
  !
  ! Stores both physical parameter bounds and transformed search-space bounds. Sampling routines
  ! operate in transformed coordinates and convert proposed values back to physical model space.
  ! **************************************************************************************************

  subroutine build_search_info(spec, search, err, message)

    implicit none

    type(parameter_spec),        intent(in)    :: spec
    type(parameter_search_info), intent(inout) :: search
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message

    integer(i4b) :: i
    integer(i4b) :: iSearch
    integer(i4b) :: nSampled

    character(len=256) :: cmessage

    err = 0
    message = 'build_search_info/'

    nSampled = 0

    do i=1,size(spec%params)
      if(spec%params(i)%sampled) nSampled = nSampled + 1
    enddo

    search%params = spec%params

    allocate(search%search_index(size(spec%params)), &
             search%param_names(nSampled),           &
             search%transformation(nSampled),        &
             search%lower(nSampled),                 &
             search%upper(nSampled),                 &
             search%search_lower(nSampled),          &
             search%search_upper(nSampled),          &
             stat=err)

    if(err/=0)then
      message=trim(message)//'unable to allocate parameter search information'
      return
    endif

    search%search_index = 0

    iSearch = 0

    do i=1,size(spec%params)

      if(.not.spec%params(i)%sampled) cycle

      iSearch = iSearch + 1

      search%search_index(i) = iSearch

      search%param_names(iSearch)     = spec%params(i)%name
      search%transformation(iSearch) = spec%params(i)%transformation

      ! physical/model-space bounds
      search%lower(iSearch) = spec%params(i)%lower
      search%upper(iSearch) = spec%params(i)%upper

      ! transformed search-space bounds
      call transform_parameter(spec%params(i)%lower,         &
                               spec%params(i)%transformation, &
                               search%search_lower(iSearch),  &
                               err,cmessage)

      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif

      call transform_parameter(spec%params(i)%upper,         &
                               spec%params(i)%transformation, &
                               search%search_upper(iSearch),  &
                               err,cmessage)

      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif

    enddo

  end subroutine build_search_info


  ! **************************************************************************************************
  ! Resolve and validate ordered parameter constraints.
  !
  ! Converts the model-provided gap fractions to absolute physical-space gaps and verifies that each
  ! ordered list contains valid parameter indices and admits at least one feasible parameter
  ! combination.
  ! **************************************************************************************************

  subroutine initialize_constraints(spec, search, err, message)

    implicit none

    type(parameter_spec),        intent(in)    :: spec
    type(parameter_search_info), intent(inout) :: search
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message

    integer(i4b) :: i
    integer(i4b) :: j
    integer(i4b) :: iConstraint
    integer(i4b) :: ixParam
    integer(i4b) :: nOrdered

    integer(i4b), allocatable :: constraint_owner(:)

    real(rkind) :: lower_feasible
    real(rkind) :: previous_value

    err = 0
    message = 'initialize_constraints/'

    ! no ordered constraints
    if(.not.allocated(spec%ordered))then
      allocate(search%ordered(0))
      return
    endif

    allocate(search%ordered(size(spec%ordered)),  &
             constraint_owner(size(spec%params)), &
             stat=err)

    if(err/=0)then
      message=trim(message)//'unable to allocate ordered parameter constraints'
      return
    endif

    constraint_owner = 0

    do iConstraint=1,size(spec%ordered)

      if(.not.allocated(spec%ordered(iConstraint)%param_index))then
        write(message,'(A,I0)') trim(message)// &
          'parameter indices not defined for ordered constraint = ',iConstraint
        err=20; return
      endif

      nOrdered = size(spec%ordered(iConstraint)%param_index)

      if(nOrdered < 2)then
        write(message,'(A,I0)') trim(message)// &
          'ordered constraint must contain at least two parameters, constraint = ',iConstraint
        err=20; return
      endif

      if(spec%ordered(iConstraint)%gap_fraction < 0._rkind)then
        write(message,'(A,I0)') trim(message)// &
          'gap_fraction must be non-negative, constraint = ',iConstraint
        err=20; return
      endif

      allocate(search%ordered(iConstraint)%param_index(nOrdered),stat=err)

      if(err/=0)then
        message=trim(message)//'unable to allocate ordered constraint information'
        return
      endif

      search%ordered(iConstraint)%param_index = &
        spec%ordered(iConstraint)%param_index

      ! ---------------------------------------------------------------------------------------------
      ! Validate parameter indices and constraint membership.
      ! ---------------------------------------------------------------------------------------------

      do i=1,nOrdered

        ixParam = search%ordered(iConstraint)%param_index(i)

        if(ixParam < 1 .or. ixParam > size(search%params))then
          write(message,'(A,I0,A,I0)') trim(message)// &
            'invalid parameter index in ordered constraint = ',iConstraint, &
            ', index = ',ixParam
          err=20; return
        endif

        ! duplicate parameters within the same ordered list are not permitted
        do j=1,i-1

          if(ixParam == search%ordered(iConstraint)%param_index(j))then
            message=trim(message)//'duplicate parameter in ordered constraint: '// &
                    trim(search%params(ixParam)%name)
            err=20; return
          endif

        enddo

        ! for now, parameters cannot occur in multiple ordered constraints
        if(constraint_owner(ixParam) > 0)then
          message=trim(message)//'parameter occurs in more than one ordered constraint: '// &
                  trim(search%params(ixParam)%name)
          err=20; return
        endif

        constraint_owner(ixParam) = iConstraint

      enddo

      ! ---------------------------------------------------------------------------------------------
      ! Compute absolute minimum gap in physical parameter space.
      ! ---------------------------------------------------------------------------------------------

      search%ordered(iConstraint)%gap = &
        spec%ordered(iConstraint)%gap_fraction * &
        ( search%params(search%ordered(iConstraint)%param_index(nOrdered))%upper - &
          search%params(search%ordered(iConstraint)%param_index(1))%lower )

      if(search%ordered(iConstraint)%gap < 0._rkind)then
        write(message,'(A,I0)') trim(message)// &
          'invalid total parameter range in ordered constraint = ',iConstraint
        err=20; return
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Check that at least one feasible parameter combination exists.
      !
      ! Sampled parameters use their smallest feasible physical value. Non-sampled parameters use
      ! their scalar trial value.
      ! ---------------------------------------------------------------------------------------------

      do i=1,nOrdered

        ixParam = search%ordered(iConstraint)%param_index(i)

        if(i == 1)then

          lower_feasible = search%params(ixParam)%lower

        else

          lower_feasible = max(search%params(ixParam)%lower, &
                               previous_value + search%ordered(iConstraint)%gap)

        endif

        if(search%params(ixParam)%sampled)then

          if(lower_feasible > search%params(ixParam)%upper)then
            write(message,'(A,A,A,ES12.5,A,ES12.5,A,A,A,ES12.5)') trim(message)// &
              'ordered constraint is infeasible near parameter: ', &
              trim(search%params(ixParam)%name), &
              ', minimum feasible value = ',lower_feasible, &
              ', upper bound = ',search%params(ixParam)%upper, &
              '; must exceed ', &
              trim(search%params(search%ordered(iConstraint)%param_index(i-1))%name), &
              ' by at least ',search%ordered(iConstraint)%gap
            err=20
            return
          endif

          previous_value=lower_feasible

        else

         if(search%params(ixParam)%trial_value < lower_feasible)then
           write(message,'(A,A,A,ES12.5,A,A,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') trim(message)// &
             'fixed parameter violates ordered constraint: ', &
             trim(search%params(ixParam)%name), &
             ', trial value = ',search%params(ixParam)%trial_value, &
             '; previous parameter = ', &
             trim(search%params(search%ordered(iConstraint)%param_index(i-1))%name), &
             ', previous value = ',previous_value, &
             ', gap = ',search%ordered(iConstraint)%gap, &
             ', parameter lower bound = ',search%params(ixParam)%lower, &
             ', minimum feasible value = ',lower_feasible
           err=20
           return
         endif

         previous_value=search%params(ixParam)%trial_value

        endif

      enddo

    enddo

  end subroutine initialize_constraints


  ! **************************************************************************************************
  ! Sample parameters uniformly over the feasible search space.
  !
  ! Each sampled parameter is drawn uniformly within its transformed search bounds, converted back
  ! to physical model space, and combined with the other parameters. Complete parameter vectors that
  ! violate an ordered physical-space constraint are rejected and resampled.
  ! **************************************************************************************************

  subroutine sample_parameters(search, param_values, err, message)

    implicit none

    type(parameter_search_info), intent(in)  :: search
    real(rkind),                 intent(out) :: param_values(:)
    integer(i4b),                intent(out) :: err
    character(*),                intent(out) :: message

    integer(i4b)            :: i
    integer(i4b)            :: ntry
    integer(i4b), parameter :: maxtry = 100000

    real(rkind) :: u
    real(rkind) :: search_value

    character(len=256) :: cmessage

    err = 0
    message = 'sample_parameters/'

    if(size(param_values) /= size(search%param_names))then
      message=trim(message)//'incorrect parameter vector size'
      err=20; return
    endif

    do ntry=1,maxtry

      do i=1,size(search%param_names)

        call random_number(u)

        search_value = search%search_lower(i) + &
                       u*(search%search_upper(i)-search%search_lower(i))

        call inverse_transform_parameter(search_value,            &
                                         search%transformation(i), &
                                         param_values(i),          &
                                         err,cmessage)

        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

      enddo

      ! ordered constraints are always evaluated in physical model space
      if(check_ordered_constraints(search,param_values)) return

    enddo

    message=trim(message)// &
            'unable to generate a parameter vector satisfying ordered constraints'
    err=20

  end subroutine sample_parameters


  ! **************************************************************************************************
  ! Perturb parameters around a current best parameter vector.
  !
  ! Best parameter values are transformed to search coordinates, perturbed using a truncated normal
  ! distribution, and converted back to physical model space. Complete parameter vectors violating
  ! an ordered physical-space constraint are rejected and resampled.
  ! **************************************************************************************************

  subroutine perturb_parameters(search, best_values, step_fraction, param_values, err, message)

    implicit none

    type(parameter_search_info), intent(in)  :: search
    real(rkind),                 intent(in)  :: best_values(:)
    real(rkind),                 intent(in)  :: step_fraction
    real(rkind),                 intent(out) :: param_values(:)
    integer(i4b),                intent(out) :: err
    character(*),                intent(out) :: message

    integer(i4b)            :: i
    integer(i4b)            :: ntry
    integer(i4b), parameter :: maxtry = 10000

    real(rkind) :: sigma
    real(rkind) :: best_search_value
    real(rkind) :: search_value

    character(len=256) :: cmessage

    err = 0
    message = 'perturb_parameters/'

    if(size(best_values) /= size(search%param_names) .or. &
       size(param_values) /= size(search%param_names))then
      message=trim(message)//'incorrect parameter vector size'
      err=20; return
    endif

    if(step_fraction <= 0._rkind)then
      message=trim(message)//'step_fraction must be greater than zero'
      err=20; return
    endif

    do ntry=1,maxtry

      do i=1,size(search%param_names)

        ! transform current best parameter to search coordinates
        call transform_parameter(best_values(i),            &
                                 search%transformation(i),   &
                                 best_search_value,          &
                                 err,cmessage)

        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

        ! perturbation scale is a fraction of transformed parameter range
        sigma = step_fraction * &
                (search%search_upper(i)-search%search_lower(i))

        call sample_truncated_normal(best_search_value,      &
                                     sigma,                  &
                                     search%search_lower(i), &
                                     search%search_upper(i), &
                                     search_value,           &
                                     err,cmessage)

        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

        ! return proposed value to physical model space
        call inverse_transform_parameter(search_value,            &
                                         search%transformation(i), &
                                         param_values(i),          &
                                         err,cmessage)

        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

      enddo

      if(check_ordered_constraints(search,param_values)) return

    enddo

    message=trim(message)// &
            'unable to generate a parameter vector satisfying ordered constraints'
    err=20

  end subroutine perturb_parameters

  ! **************************************************************************************************
  ! Generate a Dynamically Dimensioned Search (DDS) candidate solution.
  !
  ! Randomly selects a subset of the nD decision variables for inclusion in the DDS neighborhood {N}.
  ! The probability P(i) of selecting each decision variable decreases logarithmically with function
  ! evaluation i. Selected variables are perturbed around the current best solution x_best using a
  ! normal random variable with standard deviation r times the decision-variable range.
  !
  ! This implementation differs slightly from the original DDS algorithm of Tolson and Shoemaker
  ! (2007): perturbations that fall outside the parameter bounds are rejected and the full candidate
  ! is regenerated rather than reflected back into the feasible range. Ordered constraints are
  ! checked after transforming the complete candidate solution back to physical model space.
  ! **************************************************************************************************
  
  subroutine perturb_parameters_dds(search,x_best,i,m,r,x_new,err,message)
  
    implicit none
  
    type(parameter_search_info), intent(in)  :: search        ! parameter-search information
    real(rkind),                 intent(in)  :: x_best(:)     ! current best decision-variable vector
    integer(i4b),                intent(in)  :: i             ! current function-evaluation number
    integer(i4b),                intent(in)  :: m             ! maximum number of function evaluations
    real(rkind),                 intent(in)  :: r             ! DDS neighborhood perturbation size
    real(rkind),                 intent(out) :: x_new(:)      ! new candidate decision-variable vector
    integer(i4b),                intent(out) :: err           ! error code
    character(*),                intent(out) :: message       ! error message
  
    integer(i4b)            :: d                              ! decision-variable index
    integer(i4b)            :: nD                             ! total number of decision variables
    integer(i4b)            :: J                              ! number of dimensions in neighborhood {N}
    integer(i4b)            :: d_random                       ! random dimension selected if {N} is empty
    integer(i4b)            :: ntry                           ! attempt number for generating valid candidate
    integer(i4b), parameter :: maxtry=10000                   ! maximum attempts to generate valid candidate
  
    real(rkind) :: P_i                                       ! probability dimension d is included in {N}
    real(rkind) :: u                                         ! uniform random number
    real(rkind) :: z                                         ! standard normal random variable N(0,1)
    real(rkind) :: sigma_d                                   ! perturbation standard deviation for dimension d
    real(rkind) :: x_best_d                                  ! x_best(d) in transformed search space
    real(rkind) :: x_new_d                                   ! x_new(d) in transformed search space
    real(rkind) :: x_min_d                                   ! lower bound for dimension d in search space
    real(rkind) :: x_max_d                                   ! upper bound for dimension d in search space
  
    logical(lgt), dimension(size(x_best)) :: N               ! DDS neighborhood {N}; true if d is perturbed
  
    character(len=256) :: cmessage                           ! message returned by called routines
  
  
    err=0
    message='perturb_parameters_dds/'
  
    nD=size(x_best)                                          ! determine number of decision variables
  
    ! -----------------------------------------------------------------------------------------------
    ! Check input arguments.
    ! -----------------------------------------------------------------------------------------------
  
    if(size(x_new) /= nD .or. size(search%param_names) /= nD)then
      message=trim(message)//'incorrect decision-variable vector size'
      err=20
      return
    endif
  
    if(m <= 1)then
      message=trim(message)//'maximum number of function evaluations must be greater than one'
      err=20
      return
    endif
  
    if(i < 1 .or. i > m)then
      message=trim(message)//'function-evaluation index is outside valid DDS range'
      err=20
      return
    endif
  
    if(r <= 0._rkind)then
      message=trim(message)//'DDS neighborhood perturbation size must be greater than zero'
      err=20
      return
    endif
  
    ! -----------------------------------------------------------------------------------------------
    ! Compute probability that each decision variable is included in the DDS neighborhood.
    ! -----------------------------------------------------------------------------------------------
  
    P_i=1._rkind-log(real(i,rkind))/log(real(m,rkind))       ! P(i) = 1 - ln(i)/ln(m)
  
    ! -----------------------------------------------------------------------------------------------
    ! Generate a candidate solution satisfying all parameter constraints.
    ! -----------------------------------------------------------------------------------------------
  
    candidate_loop: do ntry=1,maxtry                         ! repeatedly generate candidates until one is feasible
  
      x_new=x_best                                           ! initialize candidate at current best solution
      N=.false.                                              ! initialize DDS neighborhood {N}
      J=0                                                    ! initialize number of perturbed dimensions
  
      ! ---------------------------------------------------------------------------------------------
      ! Randomly select dimensions for inclusion in neighborhood {N}.
      ! ---------------------------------------------------------------------------------------------
  
      do d=1,nD                                              ! loop through all decision variables
  
        call random_number(u)                                ! sample dimension-selection random number
  
        if(u < P_i)then
          N(d)=.true.                                        ! include dimension d in DDS neighborhood
          J=J+1                                              ! increment neighborhood dimension
        endif
  
      enddo
  
      if(J == 0)then                                         ! ensure at least one dimension is perturbed
        call random_number(u)                                ! sample one dimension uniformly
        d_random=1+int(u*real(nD,rkind))                     ! convert random number to dimension index
        d_random=min(d_random,nD)                            ! guard against upper-end roundoff
        N(d_random)=.true.                                   ! include selected dimension in neighborhood
        J=1                                                  ! neighborhood now contains one dimension
      endif
  
      ! ---------------------------------------------------------------------------------------------
      ! Perturb the selected dimensions.
      ! ---------------------------------------------------------------------------------------------
  
      do d=1,nD                                              ! loop through all decision variables
  
        if(.not.N(d)) cycle                                  ! skip decision variables outside DDS neighborhood {N}
  
        call transform_parameter(x_best(d),                  & ! transform x_best(d) from physical to search space
                                 search%transformation(d),   & ! transformation for decision variable d
                                 x_best_d,                   & ! transformed value of x_best(d)
                                 err,cmessage)                 ! error code and message
  
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif
  
        x_min_d=search%search_lower(d)                        ! lower decision-variable bound in search space
        x_max_d=search%search_upper(d)                        ! upper decision-variable bound in search space
        sigma_d=r*(x_max_d-x_min_d)                          ! sigma_d = r [x_max(d) - x_min(d)]
  
        call random_normal(z)                                ! generate standard normal random variable N(0,1)
  
        x_new_d=x_best_d+sigma_d*z                           ! perturb x_best(d) to generate x_new(d)
  
        if(x_new_d < x_min_d .or. x_new_d > x_max_d) &       ! reject candidate if perturbation violates bounds
          cycle candidate_loop
  
        call inverse_transform_parameter(x_new_d,            & ! transform x_new(d) from search to physical space
                                         search%transformation(d), & ! transformation for decision variable d
                                         x_new(d),           & ! proposed physical-space value of dimension d
                                         err,cmessage)         ! error code and message
  
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif
  
      enddo
  
      if(check_ordered_constraints(search,x_new)) return      ! accept candidate if ordered constraints are satisfied
  
    enddo candidate_loop
  
    message=trim(message)//'unable to generate DDS candidate satisfying parameter constraints'
    err=20
  
  end subroutine perturb_parameters_dds
  
  ! **************************************************************************************************
  ! Transform a parameter from physical model space to parameter-search space.
  ! **************************************************************************************************

  subroutine transform_parameter(param_value, transformation, search_value, err, message)

    implicit none

    real(rkind),  intent(in)  :: param_value
    character(*), intent(in)  :: transformation
    real(rkind),  intent(out) :: search_value
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    err = 0
    message = 'transform_parameter/'

    select case(trim(transformation))

      case ('none')

        search_value = param_value

      case ('log')

        if(param_value <= 0._rkind)then
          message=trim(message)//'log transformation requires a positive parameter value'
          err=20; return
        endif

        search_value = log(param_value)

      case default

        message=trim(message)//'unsupported parameter transformation: '// &
                trim(transformation)
        err=20; return

    end select

  end subroutine transform_parameter


  ! **************************************************************************************************
  ! Transform a parameter from parameter-search space back to physical model space.
  ! **************************************************************************************************

  subroutine inverse_transform_parameter(search_value, transformation, param_value, err, message)

    implicit none

    real(rkind),  intent(in)  :: search_value
    character(*), intent(in)  :: transformation
    real(rkind),  intent(out) :: param_value
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    err = 0
    message = 'inverse_transform_parameter/'

    select case(trim(transformation))

      case ('none')

        param_value = search_value

      case ('log')

        param_value = exp(search_value)

      case default

        message=trim(message)//'unsupported parameter transformation: '// &
                trim(transformation)
        err=20; return

    end select

  end subroutine inverse_transform_parameter


  ! **************************************************************************************************
  ! Sample from a truncated normal distribution.
  !
  ! Generates a normally distributed random sample with the specified mean and standard deviation.
  ! Samples outside the specified lower and upper bounds are rejected until a valid sample is found.
  ! **************************************************************************************************

  subroutine sample_truncated_normal(mean, sigma, lower, upper, sample, err, message)

    implicit none

    real(rkind),  intent(in)  :: mean
    real(rkind),  intent(in)  :: sigma
    real(rkind),  intent(in)  :: lower
    real(rkind),  intent(in)  :: upper
    real(rkind),  intent(out) :: sample
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    integer(i4b)            :: ntry
    integer(i4b), parameter :: maxtry = 10000

    real(rkind) :: z

    err = 0
    message = 'sample_truncated_normal/'

    if(lower > upper)then
      message=trim(message)//'lower bound exceeds upper bound'
      err=20; return
    endif

    if(mean < lower .or. mean > upper)then
      message=trim(message)//'mean is outside truncated interval'
      err=20; return
    endif

    if(sigma < 0._rkind)then
      message=trim(message)//'standard deviation must be non-negative'
      err=20; return
    endif

    if(sigma == 0._rkind .or. lower == upper)then
      sample = mean
      return
    endif

    do ntry=1,maxtry

      call random_normal(z)

      sample = mean + sigma*z

      if(sample >= lower .and. sample <= upper) return

    enddo

    message=trim(message)//'unable to generate truncated normal sample'
    err=20

  end subroutine sample_truncated_normal


  ! **************************************************************************************************
  ! Generate a standard normal random variate.
  !
  ! Uses the Box-Muller transform to generate a normally distributed random sample with zero mean
  ! and unit variance from two independent uniform random numbers.
  ! **************************************************************************************************

  subroutine random_normal(sample)

    implicit none

    real(rkind), intent(out) :: sample

    real(rkind) :: u1
    real(rkind) :: u2

    real(rkind), parameter :: pi = acos(-1._rkind)

    call random_number(u1)
    call random_number(u2)

    u1 = max(u1,tiny(1._rkind))

    sample = sqrt(-2._rkind*log(u1))*cos(2._rkind*pi*u2)

  end subroutine random_normal


  ! **************************************************************************************************
  ! Check ordered parameter constraints.
  !
  ! Uses sampled physical values for parameters included in the search vector and scalar trial
  ! values for non-sampled parameters. Returns .true. only when every ordered constraint and minimum
  ! gap is satisfied.
  ! **************************************************************************************************

  logical(lgt) function check_ordered_constraints(search,param_values)

    implicit none

    type(parameter_search_info), intent(in) :: search
    real(rkind),                 intent(in) :: param_values(:)

    integer(i4b) :: i
    integer(i4b) :: iConstraint
    integer(i4b) :: ixParam
    integer(i4b) :: ixSearch

    real(rkind) :: param_value
    real(rkind) :: previous_value

    check_ordered_constraints = .true.

    do iConstraint=1,size(search%ordered)

      do i=1,size(search%ordered(iConstraint)%param_index)

        ixParam  = search%ordered(iConstraint)%param_index(i)
        ixSearch = search%search_index(ixParam)

        if(ixSearch > 0)then
          param_value = param_values(ixSearch)
        else
          param_value = search%params(ixParam)%trial_value
        endif

        if(i > 1)then

          if(param_value < previous_value + search%ordered(iConstraint)%gap)then
            check_ordered_constraints = .false.
            return
          endif

        endif

        previous_value = param_value

      enddo

    enddo

  end function check_ordered_constraints


end module parameter_search
