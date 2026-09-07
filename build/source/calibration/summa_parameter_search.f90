module summa_parameter_search

  USE nr_type,    only: i4b, rkind, lgt
  USE summa_type, only: config_info

  implicit none
  private

  ! Resolved information for one ordered constraint
  type, public :: ordered_search_info

    character(len=64), allocatable :: param_names(:)     ! Parameter names in increasing order
    integer(i4b),      allocatable :: calib_index(:)     ! Index in calibration vector; 0 = fixed
    logical(lgt),      allocatable :: sampled(:)         ! .true. if parameter is calibrated

    real(rkind),       allocatable :: lower(:)           ! Native lower bounds
    real(rkind),       allocatable :: upper(:)           ! Native upper bounds
    real(rkind),       allocatable :: fixed_value(:)     ! Default value for fixed parameters
    real(rkind),       allocatable :: upper_feasible(:)  ! Maximum feasible value after constraints

    real(rkind)                    :: gap_fraction        ! Minimum gap as fraction of total range
    real(rkind)                    :: gap                 ! Absolute minimum gap

  end type ordered_search_info


  ! Resolved information for the complete parameter search
  type, public :: parameter_search_info

    character(len=64), allocatable        :: param_names(:)   ! Calibration parameter names
    real(rkind),       allocatable        :: lower(:)         ! Calibration lower bounds
    real(rkind),       allocatable        :: upper(:)         ! Calibration upper bounds
    logical(lgt),      allocatable        :: constrained(:)   ! Parameter occurs in ordered constraint

    type(ordered_search_info), allocatable :: ordered(:)      ! Resolved ordered constraints

  end type parameter_search_info


  public :: get_parameter_bounds
  public :: initialize_parameter_search
  public :: sample_parameters
  public :: perturb_parameters

contains


  ! **************************************************************************************************
  ! Initialize the calibration parameter search.
  !
  ! Resolves parameter names against SUMMA parameter metadata, obtains parameter bounds and default
  ! values, identifies parameters participating in ordered constraints, computes minimum gaps, and
  ! verifies that all ordered constraints have a feasible parameter space.
  !
  ! This routine is intended to be called once before repeated parameter sampling.
  ! **************************************************************************************************

  subroutine initialize_parameter_search(config, search, err, message)

    implicit none

    type(config_info),          intent(in)  :: config
    type(parameter_search_info),intent(out) :: search
    integer(i4b),               intent(out) :: err
    character(*),               intent(out) :: message

    integer(i4b) :: i
    integer(i4b) :: j
    integer(i4b) :: k
    integer(i4b) :: iConstraint
    integer(i4b) :: nConstraint
    integer(i4b) :: nOrdered
    integer(i4b) :: ixCalib

    real(rkind) :: lower_feasible
    real(rkind) :: previous_value
    real(rkind) :: upper_candidate

    character(len=256) :: cmessage

    err = 0
    message = 'initialize_parameter_search/'

    ! calibration parameter list is required
    if(.not.allocated(config%calib%param_list))then
      message=trim(message)//'calibration parameter list is not defined'
      err=20; return
    endif

    ! check that there are actually parameters to calibrate
    if(size(config%calib%param_list) == 0)then
      message=trim(message)//'calibration parameter list is empty'
      err=20; return
    endif

    ! save calibration parameter names
    search%param_names = config%calib%param_list

    ! check for duplicate calibration parameters
    do i=1,size(search%param_names)-1
      do j=i+1,size(search%param_names)
    
        if(trim(search%param_names(i)) == trim(search%param_names(j)))then
          message=trim(message)//'duplicate calibration parameter: '// &
                  trim(search%param_names(i))
          err=20; return
        endif
    
      enddo
    enddo

    ! obtain native bounds for calibration parameters
    call get_parameter_bounds(search%param_names, search%lower, search%upper, &
                              err, cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    ! identify parameters participating in ordered constraints
    allocate(search%constrained(size(search%param_names)), stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate constrained parameter mask'
      return
    endif

    search%constrained = .false.

    ! no ordered constraints
    if(.not.allocated(config%calib%ordered))then
      allocate(search%ordered(0))
      return
    endif

    nConstraint = size(config%calib%ordered)

    allocate(search%ordered(nConstraint), stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate ordered parameter constraints'
      return
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Resolve each ordered constraint.
    ! -----------------------------------------------------------------------------------------------

    do iConstraint=1,nConstraint

      if(.not.allocated(config%calib%ordered(iConstraint)%parameters))then
        write(message,'(A,I0)') trim(message)// &
          'parameter list not defined for ordered constraint = ',iConstraint
        err=20; return
      endif

      nOrdered = size(config%calib%ordered(iConstraint)%parameters)

      if(nOrdered < 2)then
        write(message,'(A,I0)') trim(message)// &
          'ordered constraint must contain at least two parameters, constraint = ',iConstraint
        err=20; return
      endif

      if(config%calib%ordered(iConstraint)%gap_fraction < 0._rkind)then
        write(message,'(A,I0)') trim(message)// &
          'gap_fraction must be non-negative, constraint = ',iConstraint
        err=20; return
      endif

      allocate(search%ordered(iConstraint)%param_names(nOrdered),     &
               search%ordered(iConstraint)%calib_index(nOrdered),     &
               search%ordered(iConstraint)%sampled(nOrdered),         &
               search%ordered(iConstraint)%lower(nOrdered),           &
               search%ordered(iConstraint)%upper(nOrdered),           &
               search%ordered(iConstraint)%fixed_value(nOrdered),     &
               search%ordered(iConstraint)%upper_feasible(nOrdered),  &
               stat=err)

      if(err/=0)then
        message=trim(message)//'unable to allocate ordered constraint information'
        return
      endif

      search%ordered(iConstraint)%param_names = &
        config%calib%ordered(iConstraint)%parameters

      search%ordered(iConstraint)%gap_fraction = &
        config%calib%ordered(iConstraint)%gap_fraction

      ! ---------------------------------------------------------------------------------------------
      ! Resolve each parameter name.
      ! ---------------------------------------------------------------------------------------------

      do i=1,nOrdered

        ! reject duplicate parameters within the same constraint
        do j=1,i-1
          if(trim(search%ordered(iConstraint)%param_names(i)) == &
             trim(search%ordered(iConstraint)%param_names(j)))then

            message=trim(message)//'duplicate parameter in ordered constraint: '// &
                    trim(search%ordered(iConstraint)%param_names(i))
            err=20; return
          endif
        enddo

        ! for the first implementation, do not allow overlapping ordered constraints
        do j=1,iConstraint-1
          do k=1,size(search%ordered(j)%param_names)

            if(trim(search%ordered(iConstraint)%param_names(i)) == &
               trim(search%ordered(j)%param_names(k)))then

              message=trim(message)//'parameter occurs in more than one ordered constraint: '// &
                      trim(search%ordered(iConstraint)%param_names(i))
              err=20; return
            endif

          enddo
        enddo

        ! obtain native bounds and default parameter value
        call get_parameter_info(                                      &
               trim(search%ordered(iConstraint)%param_names(i)),      &
               search%ordered(iConstraint)%fixed_value(i),            &
               search%ordered(iConstraint)%lower(i),                  &
               search%ordered(iConstraint)%upper(i),                  &
               err,cmessage)

        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

        ! determine whether this parameter is included in the calibration vector
        ixCalib = find_parameter(                                    &
                    trim(search%ordered(iConstraint)%param_names(i)), &
                    search%param_names)

        search%ordered(iConstraint)%calib_index(i) = ixCalib
        search%ordered(iConstraint)%sampled(i)     = (ixCalib > 0)

        ! mark this calibration parameter as constrained
        if(ixCalib > 0) search%constrained(ixCalib) = .true.

      enddo

      ! ---------------------------------------------------------------------------------------------
      ! Compute the absolute minimum gap.
      !
      ! gap_fraction is relative to the complete admissible range from the lower bound of the first
      ! member to the upper bound of the last member.
      ! ---------------------------------------------------------------------------------------------

      search%ordered(iConstraint)%gap = &
        search%ordered(iConstraint)%gap_fraction * &
        (search%ordered(iConstraint)%upper(nOrdered) - &
         search%ordered(iConstraint)%lower(1))

      if(search%ordered(iConstraint)%gap < 0._rkind)then
        write(message,'(A,I0)') trim(message)// &
          'invalid parameter range in ordered constraint = ',iConstraint
        err=20; return
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Propagate feasible upper bounds backwards.
      !
      ! Fixed parameters are represented by their default value. This prevents an earlier sampled
      ! parameter from being chosen so large that a later fixed or sampled parameter cannot satisfy
      ! the required minimum gap.
      ! ---------------------------------------------------------------------------------------------

      i = nOrdered

      if(search%ordered(iConstraint)%sampled(i))then
        search%ordered(iConstraint)%upper_feasible(i) = &
          search%ordered(iConstraint)%upper(i)
      else
        search%ordered(iConstraint)%upper_feasible(i) = &
          search%ordered(iConstraint)%fixed_value(i)
      endif

      do i=nOrdered-1,1,-1

        if(search%ordered(iConstraint)%sampled(i))then

          upper_candidate = search%ordered(iConstraint)%upper(i)

          search%ordered(iConstraint)%upper_feasible(i) = &
            min(upper_candidate, &
                search%ordered(iConstraint)%upper_feasible(i+1) - &
                search%ordered(iConstraint)%gap)

        else

          ! fixed parameter must leave enough room for the following parameter
          if(search%ordered(iConstraint)%fixed_value(i) > &
             search%ordered(iConstraint)%upper_feasible(i+1) - &
             search%ordered(iConstraint)%gap)then

            message=trim(message)//'fixed parameter violates ordered constraint: '// &
                    trim(search%ordered(iConstraint)%param_names(i))
            err=20; return
          endif

          search%ordered(iConstraint)%upper_feasible(i) = &
            search%ordered(iConstraint)%fixed_value(i)

        endif

      enddo

      ! ---------------------------------------------------------------------------------------------
      ! Forward feasibility check using the smallest feasible values.
      ! ---------------------------------------------------------------------------------------------

      do i=1,nOrdered

        if(i == 1)then
          lower_feasible = search%ordered(iConstraint)%lower(i)
        else
          lower_feasible = max(search%ordered(iConstraint)%lower(i), &
                               previous_value + search%ordered(iConstraint)%gap)
        endif

        if(search%ordered(iConstraint)%sampled(i))then

          if(lower_feasible > search%ordered(iConstraint)%upper_feasible(i))then
            message=trim(message)//'ordered constraint is infeasible near parameter: '// &
                    trim(search%ordered(iConstraint)%param_names(i))
            err=20; return
          endif

          ! use the lowest possible value when checking feasibility downstream
          previous_value = lower_feasible

        else

          ! fixed parameter must satisfy both the lower constraint and propagated upper constraint
          if(search%ordered(iConstraint)%fixed_value(i) < lower_feasible .or. &
             search%ordered(iConstraint)%fixed_value(i) > &
             search%ordered(iConstraint)%upper_feasible(i))then

            message=trim(message)//'fixed parameter violates ordered constraint: '// &
                    trim(search%ordered(iConstraint)%param_names(i))
            err=20; return
          endif

          previous_value = search%ordered(iConstraint)%fixed_value(i)

        endif

      enddo

    enddo

  end subroutine initialize_parameter_search

  ! **************************************************************************************************
  ! Sample calibration parameters uniformly over the feasible parameter space.
  !
  ! Samples each calibration parameter independently from its native lower and upper bounds.
  ! Parameter vectors that violate an ordered constraint are rejected and resampled. The resulting
  ! accepted samples are uniformly distributed over the feasible region defined by the parameter
  ! bounds and ordered constraints.
  ! **************************************************************************************************
  
  subroutine sample_parameters(search, param_value, err, message)
  
    implicit none
  
    type(parameter_search_info), intent(in)  :: search
    real(rkind),                 intent(out) :: param_value(:)
    integer(i4b),                intent(out) :: err
    character(*),                intent(out) :: message
  
    integer(i4b)            :: i
    integer(i4b)            :: ntry
    integer(i4b), parameter :: maxtry = 100000
  
    real(rkind) :: u
  
    err = 0
    message = 'sample_parameters/'
  
    ! check dimensions
    if(size(param_value) /= size(search%param_names))then
      message=trim(message)//'incorrect parameter vector size'
      err=20; return
    endif
  
    ! generate parameter vectors until all constraints are satisfied
    do ntry=1,maxtry
  
      ! sample uniformly from native parameter bounds
      do i=1,size(search%param_names)
  
        call random_number(u)
  
        param_value(i) = search%lower(i) + &
                         u*(search%upper(i)-search%lower(i))
  
      enddo
  
      ! accept the complete parameter vector if all constraints are satisfied
      if(check_ordered_constraints(search,param_value)) return
  
    enddo
  
    message=trim(message)// &
            'unable to generate a parameter vector satisfying ordered constraints'
    err=20
  
  end subroutine sample_parameters


  ! **************************************************************************************************
  ! Perturb calibration parameters around a current best parameter vector.
  !
  ! Each calibration parameter is perturbed independently using a normal distribution centered
  ! on the current best value. The standard deviation is defined as a fraction of the native
  ! parameter range. Proposed parameter vectors that violate an ordered parameter constraint
  ! are rejected and resampled.
  ! **************************************************************************************************
  
  subroutine perturb_parameters(search, best_value, step_fraction, param_value, err, message)
  
    implicit none
  
    type(parameter_search_info), intent(in)  :: search
    real(rkind),                 intent(in)  :: best_value(:)
    real(rkind),                 intent(in)  :: step_fraction
    real(rkind),                 intent(out) :: param_value(:)
    integer(i4b),                intent(out) :: err
    character(*),                intent(out) :: message
  
    integer(i4b)            :: i
    integer(i4b)            :: ntry
    integer(i4b), parameter :: maxtry = 10000
  
    real(rkind) :: sigma
  
    character(len=256) :: cmessage
  
    err = 0
    message = 'perturb_parameters/'
  
    ! check dimensions
    if(size(best_value) /= size(search%param_names) .or. &
       size(param_value) /= size(search%param_names))then
      message=trim(message)//'incorrect parameter vector size'
      err=20; return
    endif
  
    if(step_fraction <= 0._rkind)then
      message=trim(message)//'step_fraction must be greater than zero'
      err=20; return
    endif
  
    ! -----------------------------------------------------------------------------------------------
    ! Generate parameter vectors until all constraints are satisfied.
    ! -----------------------------------------------------------------------------------------------
  
    do ntry=1,maxtry
  
      ! perturb every calibration parameter independently
      do i=1,size(search%param_names)
  
        sigma = step_fraction*(search%upper(i)-search%lower(i))
  
        call sample_truncated_normal(best_value(i),   &
                                     sigma,           &
                                     search%lower(i), &
                                     search%upper(i), &
                                     param_value(i),  &
                                     err,cmessage)
  
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif
  
      enddo
  
      ! check ordered parameter constraints
      if(check_ordered_constraints(search,param_value)) return
  
    enddo
  
    message=trim(message)//'unable to generate a parameter vector satisfying ordered constraints'
    err=20
  
  end subroutine perturb_parameters

  ! **************************************************************************************************
  ! Sample from a truncated normal distribution.
  !
  ! Generates a normally distributed random value with the specified mean and standard deviation.
  ! Values outside the specified lower and upper bounds are rejected and resampled until a value
  ! within the interval is obtained.
  ! **************************************************************************************************

  subroutine sample_truncated_normal(mean, sigma, lower, upper, value, err, message)

    implicit none
  
    real(rkind), intent(in)  :: mean
    real(rkind), intent(in)  :: sigma
    real(rkind), intent(in)  :: lower
    real(rkind), intent(in)  :: upper
    real(rkind), intent(out) :: value
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
      value = mean
      return
    endif
  
    do ntry=1,maxtry
  
      call random_normal(z)
  
      value = mean + sigma*z
  
      if(value >= lower .and. value <= upper) return
  
    enddo
  
    message=trim(message)//'unable to generate truncated normal sample'
    err=20
  
  end subroutine sample_truncated_normal

  ! **************************************************************************************************
  ! Generate a standard normal random variate.
  !
  ! Uses the Box-Muller transform to generate a normally distributed random value with zero mean
  ! and unit variance from two independent uniform random numbers.
  ! **************************************************************************************************

  subroutine random_normal(z)
  
    implicit none
  
    real(rkind), intent(out) :: z
  
    real(rkind) :: u1
    real(rkind) :: u2
    real(rkind), parameter :: pi = acos(-1._rkind)
  
    call random_number(u1)
    call random_number(u2)
  
    u1 = max(u1,tiny(1._rkind))
  
    z = sqrt(-2._rkind*log(u1))*cos(2._rkind*pi*u2)
  
  end subroutine random_normal


  ! **************************************************************************************************
  ! Check ordered parameter constraints.
  !
  ! Evaluates each ordered parameter list using sampled values for calibration parameters and
  ! default values for parameters that are not included in the calibration vector. Returns
  ! .true. only when all ordered constraints and minimum gaps are satisfied.
  ! **************************************************************************************************
  
  logical(lgt) function check_ordered_constraints(search, param_value)
  
    implicit none
  
    type(parameter_search_info), intent(in) :: search
    real(rkind),                 intent(in) :: param_value(:)
  
    integer(i4b) :: i
    integer(i4b) :: iConstraint
    integer(i4b) :: ixCalib
  
    real(rkind) :: value
    real(rkind) :: previous_value
  
    check_ordered_constraints = .true.
  
    do iConstraint=1,size(search%ordered)
  
      do i=1,size(search%ordered(iConstraint)%param_names)
  
        if(search%ordered(iConstraint)%sampled(i))then
          ixCalib = search%ordered(iConstraint)%calib_index(i)
          value   = param_value(ixCalib)
        else
          value = search%ordered(iConstraint)%fixed_value(i)
        endif
  
        if(i > 1)then
          if(value < previous_value + search%ordered(iConstraint)%gap)then
            check_ordered_constraints = .false.
            return
          endif
        endif
  
        previous_value = value
  
      enddo
  
    enddo
  
  end function check_ordered_constraints


  ! **************************************************************************************************
  ! Get parameter bounds for calibration.
  !
  ! Maps each requested parameter name to the corresponding SUMMA local or basin parameter and
  ! extracts the lower and upper limits defined in the parameter metadata. Returns an error if a
  ! requested parameter cannot be found.
  ! **************************************************************************************************

  subroutine get_parameter_bounds(param_names, lower, upper, err, message)

    implicit none

    character(len=*), intent(in)          :: param_names(:)
    real(rkind), allocatable, intent(out) :: lower(:)
    real(rkind), allocatable, intent(out) :: upper(:)

    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    integer(i4b) :: i
    real(rkind)  :: value

    character(len=256) :: cmessage

    err = 0
    message = 'get_parameter_bounds/'

    allocate(lower(size(param_names)), upper(size(param_names)), stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate parameter bounds'
      return
    endif

    do i=1,size(param_names)

      call get_parameter_info(trim(param_names(i)), value, lower(i), upper(i), &
                              err,cmessage)

      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif

      ! check parameter bounds
      if(lower(i) > upper(i))then
        message=trim(message)//'invalid bounds for parameter: '//trim(param_names(i))
        err=20; return
      endif

    enddo

  end subroutine get_parameter_bounds


  ! **************************************************************************************************
  ! Get SUMMA parameter metadata.
  !
  ! Returns the default value and native lower and upper limits for a local or basin parameter.
  ! **************************************************************************************************

  subroutine get_parameter_info(param_name, value, lower, upper, err, message)

    USE get_ixname_module, only: get_ixParam
    USE get_ixname_module, only: get_ixBpar

    USE globalData, only: localParFallback
    USE globalData, only: basinParFallback

    implicit none

    character(*), intent(in)  :: param_name
    real(rkind),  intent(out) :: value
    real(rkind),  intent(out) :: lower
    real(rkind),  intent(out) :: upper
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    integer(i4b) :: ixParam
    integer(i4b) :: ixBasin

    err = 0
    message = 'get_parameter_info/'

    ! local parameter
    ixParam = get_ixParam(trim(param_name))

    if(ixParam > 0)then

      value = localParFallback(ixParam)%default_val
      lower = localParFallback(ixParam)%lower_limit
      upper = localParFallback(ixParam)%upper_limit
      
      if(lower > upper)then
        message=trim(message)//'invalid bounds for parameter: '//trim(param_name)
        err=20; return
      endif
      
      return

    endif

    ! basin parameter
    ixBasin = get_ixBpar(trim(param_name))

    if(ixBasin > 0)then

      value = basinParFallback(ixBasin)%default_val
      lower = basinParFallback(ixBasin)%lower_limit
      upper = basinParFallback(ixBasin)%upper_limit
      
      if(lower > upper)then
        message=trim(message)//'invalid bounds for parameter: '//trim(param_name)
        err=20; return
      endif
      
      return

    endif

    message=trim(message)//'parameter not found: '//trim(param_name)
    err=20

  end subroutine get_parameter_info


  ! **************************************************************************************************
  ! Find a parameter in a parameter-name vector.
  !
  ! Returns the one-based parameter index, or zero if the parameter is not present.
  ! **************************************************************************************************

  integer(i4b) function find_parameter(param_name, param_names)

    implicit none

    character(*), intent(in) :: param_name
    character(*), intent(in) :: param_names(:)

    integer(i4b) :: i

    find_parameter = 0

    do i=1,size(param_names)

      if(trim(param_names(i)) == trim(param_name))then
        find_parameter = i
        return
      endif

    enddo

  end function find_parameter


end module summa_parameter_search
