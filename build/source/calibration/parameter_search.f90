module parameter_search

  USE nr_type, only: i4b, rkind, lgt
  USE, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  ! **************************************************************************************************
  ! Parameterization types used to construct physical parameter fields from calibration variables.
  !
  ! The calibration variable x can represent:
  !   1. A physical parameter value
  !   2. A spatial multiplier applied to an initialized parameter field
  !   3. A parameter ratio applied to another parameter's effective physical values
  !
  ! The effective physical parameter values are:
  !
  !   Physical parameter:
  !     P(i) = x
  !
  !   Spatial multiplier:
  !     P(i) = x * P_init(i)
  !
  !   Parameter ratio:
  !     P(i) = x * P_ref(i)
  !
  ! Here, i indexes the spatial elements of the parameter field.
  ! Physical parameters are assigned a uniform value across all elements.
  ! Spatial multipliers preserve the variability of the initialized field,
  ! while parameter ratios use the effective values of a reference parameter.
  ! **************************************************************************************************
  integer(i4b), parameter, public :: param_physical = 1
  integer(i4b), parameter, public :: param_spatial_multiplier = 2
  integer(i4b), parameter, public :: param_ratio = 3

  ! **************************************************************************************************
  ! State information used to generate parameter samples during a search.
  !
  ! The elite individual contains the parameter vector with the best fitness identified during the
  ! search. The population contains a collection of parameter vectors and their associated fitness
  ! values. Different search methods use different components of this structure. For example, DDS
  ! generates new samples by perturbing the elite individual, whereas population-based methods such
  ! as genetic algorithms and differential evolution generate new samples from the population.
  ! **************************************************************************************************
  type, public :: search_state_type
    real(rkind), allocatable :: elite_individual(:)
    real(rkind)              :: elite_fitness
    real(rkind), allocatable :: population(:,:)
    real(rkind), allocatable :: fitness(:)
  end type search_state_type

  ! **************************************************************************************************
  ! Information describing one model parameter.
  !
  ! The model-specific layer populates this structure. Parameters with sampled=.true. are included
  ! in the search vector, with one scalar calibration variable per parameter. The param_type
  ! determines how this variable is used to construct the effective physical parameter field:
  !
  !   - Physical parameter: assign a spatially uniform value
  !   - Spatial multiplier: scale the initialized parameter field
  !   - Parameter ratio: scale the effective field of a reference parameter
  !
  ! Parameters with sampled=.false. retain their initialized values but may participate in
  ! ordered constraints or serve as reference parameters. The reference_name identifies the
  ! parameter used in a ratio, without assuming how its physical values are stored.
  !
  ! The lower and upper bounds apply to the scalar calibration variable, not necessarily
  ! to the resulting physical parameter values.
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
    integer(i4b)       :: param_type = param_physical   ! Physical value, spatial multiplier, or parameter ratio
    character(len=64)  :: reference_name = ''           ! Name of the reference parameter for a parameter ratio
    character(len=16)  :: reference_reduction = 'none'  ! none: elementwise; min: local reference minimum
  end type parameter_info

  ! **************************************************************************************************
  ! Ordered parameter dependency.
  !
  ! param_index contains indices into parameter_spec%params. Adjacent parameters must be separated
  ! by at least gap_fraction times the complete range from the lower bound of the first parameter
  ! to the upper bound of the last parameter.
  !
  ! Constraint evaluation is owned by parameter_constraints, not by the optimizers.
  ! **************************************************************************************************
  type, public :: ordered_constraint
    integer(i4b), allocatable :: param_index(:)
    real(rkind)               :: gap_fraction = 0._rkind
  end type ordered_constraint

  ! **************************************************************************************************
  ! Multiplicative parameter dependency.
  !
  ! The dependent parameter is sampled as a multiplier of the reference parameter:
  !
  !   physical(dependent) = physical(reference) * sampled(dependent)
  !
  ! Indices refer to parameter_spec%params.
  ! **************************************************************************************************
  type, public :: multiplier_constraint
    integer(i4b) :: reference_index = 0
    integer(i4b) :: dependent_index = 0
  end type multiplier_constraint

  ! **************************************************************************************************
  ! Model-provided specification of the parameter search problem.
  !
  ! Any model can use the parameter-search routines by populating this structure with scalar trial
  ! values, bounds, sampled flags, transformations, and optional ordered dependencies.
  ! **************************************************************************************************
  type, public :: parameter_spec
    type(parameter_info),        allocatable :: params(:)
    type(ordered_constraint),    allocatable :: ordered(:)
    type(multiplier_constraint), allocatable :: multiplier(:)
  end type parameter_spec

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
  end type parameter_search_info

  public :: initialize_parameter_search
  public :: generate_search_sample
  public :: sample_parameters
  public :: perturb_parameters
  public :: perturb_parameters_dds
  public :: transform_parameter
  public :: inverse_transform_parameter
  ! The optimizer knows only the candidate and the yes/no/error result.
  abstract interface
    subroutine feasibility_callback(context,x,feasible,err,message)
      import :: rkind,i4b,lgt
      class(*), intent(inout) :: context
      real(rkind), intent(in) :: x(:)
      logical(lgt), intent(out) :: feasible
      integer(i4b), intent(out) :: err
      character(*), intent(out) :: message
    end subroutine feasibility_callback
  end interface
  public :: feasibility_callback

contains

  ! **************************************************************************************************
  ! Initialize the parameter search.
  !
  ! Validates the model-provided parameter specification, constructs the sampled parameter vector,
  ! resolves parameter transformations, and retains the configured search bounds.
  ! Constraint initialization and feasibility checking are separate operations.
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
      if(.not.ieee_is_finite(spec%params(i)%lower) .or. &
         .not.ieee_is_finite(spec%params(i)%upper))then
        message=trim(message)//'non-finite bounds for parameter: '//trim(spec%params(i)%name)
        err=20; return
      endif
      if(spec%params(i)%lower > spec%params(i)%upper)then
        message=trim(message)//'invalid bounds for parameter: '//trim(spec%params(i)%name)
        err=20; return
      endif

      ! Non-sampled parameters retain their initialized model values.
      ! Warn if these values fall outside the nominal parameter bounds.
      if(.not.spec%params(i)%sampled)then
        if(spec%params(i)%trial_value < spec%params(i)%lower .or. &
           spec%params(i)%trial_value > spec%params(i)%upper)then
          write(*,'(A,A,A,ES12.5,A,ES12.5,A,ES12.5)') &
            'WARNING: initialized value outside nominal bounds for parameter ', &
            trim(spec%params(i)%name), &
            ', value = ',spec%params(i)%trial_value, &
            ', lower = ',spec%params(i)%lower, &
            ', upper = ',spec%params(i)%upper
        endif
      endif

      ! validate parameter transformation
      select case(trim(spec%params(i)%transformation))

        case ('none')

          ! no additional requirements
        case ('log')
          if(spec%params(i)%sampled)then
            if(spec%params(i)%lower <= 0._rkind .or. spec%params(i)%upper <= 0._rkind)then
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
    character(len=1024) :: cmessage

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
      call transform_parameter(spec%params(i)%lower,spec%params(i)%transformation, &
                               search%search_lower(iSearch),err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif
      call transform_parameter(spec%params(i)%upper,spec%params(i)%transformation, &
                               search%search_upper(iSearch),err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif
    enddo

  end subroutine build_search_info

  ! **************************************************************************************************
  ! Generate the next parameter sample for the selected search method.
  !
  ! The search method determines how the next individual is generated from the current search state.
  ! Random sampling draws independently from the feasible parameter space. DDS perturbs the current
  ! elite individual. If no elite individual has yet been established, DDS generates an independent
  ! random sample to initialize the search.
  ! **************************************************************************************************
  subroutine generate_search_sample(method,search,state,i,m,param_values,is_feasible,context,err,message)
    implicit none
    procedure(feasibility_callback) :: is_feasible
    class(*), intent(inout) :: context

    character(*),                intent(in)    :: method
    type(parameter_search_info), intent(in)    :: search
    type(search_state_type),     intent(in)    :: state
    integer(i4b),                intent(in)    :: i
    integer(i4b),                intent(in)    :: m
    real(rkind),                 intent(out)   :: param_values(:)
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message

    real(rkind), parameter :: dds_r=0.2_rkind
    character(len=1024) :: cmessage

    err=0
    message='generate_search_sample/'

    select case(trim(method))

      case ('random')
        call sample_parameters(search,param_values,is_feasible,context,err,cmessage)
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

      case ('dds')
        ! Initialize DDS with independent random samples until an elite individual exists.
        if(.not.allocated(state%elite_individual))then
          call sample_parameters(search,param_values,is_feasible,context,err,cmessage)
          if(err/=0)then
            message=trim(message)//trim(cmessage)
            return
          endif

        else
          call perturb_parameters_dds(search,state%elite_individual,i,m, &
                                      dds_r,param_values,is_feasible,context,err,cmessage)
          if(err/=0)then
            message=trim(message)//trim(cmessage)
            return
          endif

        endif

      case default
        message=trim(message)//'unknown parameter sampling method: '//trim(method)
        err=20
        return

    end select

  end subroutine generate_search_sample


  ! **************************************************************************************************
  ! Sample parameters uniformly over the feasible search space.
  !
  ! Each sampled parameter is drawn uniformly within its transformed search bounds, converted back
  ! to physical model space, and combined with the other parameters. The supplied feasibility predicate is
  ! the acceptance decision; an infeasible complete vector is regenerated.
  ! **************************************************************************************************
  subroutine sample_parameters(search, param_values,is_feasible,context,err,message)
    implicit none
    procedure(feasibility_callback) :: is_feasible
    class(*), intent(inout) :: context
    type(parameter_search_info), intent(in)  :: search
    real(rkind),                 intent(out) :: param_values(:)
    integer(i4b),                intent(out) :: err
    character(*),                intent(out) :: message
    integer(i4b)            :: i
    integer(i4b)            :: ntry
    integer(i4b), parameter :: maxtry = 100000
    real(rkind) :: u
    real(rkind) :: search_value
    character(len=1024) :: cmessage

    logical(lgt) :: feasible

    err = 0
    message = 'sample_parameters/'
    if(size(param_values) /= size(search%param_names))then
      message=trim(message)//'incorrect parameter vector size'
      err=20; return
    endif
    do ntry=1,maxtry
      do i=1,size(search%param_names)
        call random_number(u)
        search_value = search%search_lower(i) + u*(search%search_upper(i)-search%search_lower(i))
        call inverse_transform_parameter(search_value,search%transformation(i),param_values(i), &
                                         err,cmessage)
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif
      enddo

      ! The supplied checker owns all model feasibility rules.
      call is_feasible(context,param_values,feasible,err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif
      if(feasible) return
    enddo
    message=trim(message)// 'candidate retry limit reached; last rejection: '//trim(cmessage)
    err=20

  end subroutine sample_parameters

  ! **************************************************************************************************
  ! Perturb parameters around a current best parameter vector.
  !
  ! Best parameter values are transformed to search coordinates, perturbed using a truncated normal
  ! distribution, and converted back to physical model space. Complete parameter vectors violating
  ! an ordered physical-space constraint are rejected and resampled.
  ! **************************************************************************************************
  subroutine perturb_parameters(search, best_values, step_fraction, param_values,is_feasible,context,err,message)
    implicit none
    procedure(feasibility_callback) :: is_feasible
    class(*), intent(inout) :: context
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
    character(len=1024) :: cmessage

    logical(lgt) :: feasible

    err = 0
    message = 'perturb_parameters/'
    if(size(best_values) /= size(search%param_names) .or. size(param_values) /= size(search%param_names))then
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
        call transform_parameter(best_values(i), search%transformation(i), best_search_value, err,cmessage)
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

        ! perturbation scale is a fraction of transformed parameter range
        sigma = step_fraction * (search%search_upper(i)-search%search_lower(i))
        call sample_truncated_normal(best_search_value,sigma,search%search_lower(i), &
                                     search%search_upper(i),search_value,err,cmessage)
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

        ! return proposed value to physical model space
        call inverse_transform_parameter(search_value,search%transformation(i), &
                                         param_values(i),err,cmessage)
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif
      enddo
      call is_feasible(context,param_values,feasible,err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif
      if(feasible) return
    enddo
    message=trim(message)// 'candidate retry limit reached; last rejection: '//trim(cmessage)
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
  subroutine perturb_parameters_dds(search,x_best,i,m,r,x_new,is_feasible,context,err,message)
    implicit none
    procedure(feasibility_callback) :: is_feasible
    class(*), intent(inout) :: context
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
    character(len=1024) :: cmessage                           ! message returned by called routines

    logical(lgt) :: feasible
    character(len=1024) :: last_rejection

    err=0
    message='perturb_parameters_dds/'
    last_rejection='proposal outside configured bounds'
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
        if(x_new_d < x_min_d .or. x_new_d > x_max_d)then
          last_rejection='proposal outside configured bounds'
          cycle candidate_loop
        endif
        call inverse_transform_parameter(x_new_d,            & ! transform x_new(d) from search to physical space
                                         search%transformation(d), & ! transformation for decision variable d
                                         x_new(d),           & ! proposed physical-space value of dimension d
                                         err,cmessage)         ! error code and message
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif
      enddo
      call is_feasible(context,x_new,feasible,err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif
      if(feasible) return
      last_rejection=cmessage
    enddo candidate_loop
    message=trim(message)//'DDS candidate retry limit reached; last rejection: '//trim(last_rejection)
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
        message=trim(message)//'unsupported parameter transformation: '// trim(transformation)
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
        message=trim(message)//'unsupported parameter transformation: '// trim(transformation)
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



end module parameter_search
