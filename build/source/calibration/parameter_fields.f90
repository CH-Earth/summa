
module parameter_fields

  USE nr_type, only: i4b,rkind,lgt
  USE parameter_search, only: parameter_spec,parameter_search_info
  USE parameter_search, only: param_physical,param_spatial_multiplier,param_ratio
  USE, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  integer(i4b), parameter, public :: err_ordered_violation=252000

  ! **************************************************************************************************
  ! Physical values of one model parameter.
  !
  ! initial contains the parameter values obtained from model initialization.
  ! effective contains the values used for the current calibration trial.
  !
  ! The two arrays represent the same spatial elements in the same order.
  ! **************************************************************************************************
  type, public :: parameter_field
    real(rkind), allocatable :: initial(:)
    real(rkind), allocatable :: effective(:)
  end type parameter_field

  public :: reset_parameter_fields
  public :: apply_parameter_sample
  public :: compute_feasible_interval
  public :: validate_effective_ordered_constraints

contains

  ! **************************************************************************************************
  ! Reset the effective parameter fields to their initialized values.
  !
  ! Called before applying a new calibration trial. This prevents spatial
  ! multipliers and parameter ratios from accumulating changes across trials.
  !
  ! All initial fields must have been allocated by the calling model.
  ! **************************************************************************************************
  subroutine reset_parameter_fields(field)
    implicit none

    type(parameter_field), intent(inout) :: field(:)
    integer(i4b) :: i

    do i=1,size(field)
      field(i)%effective=field(i)%initial
    enddo

  end subroutine reset_parameter_fields

  ! **************************************************************************************************
  ! Construct effective physical parameter fields for one calibration trial.
  !
  ! Each parameter has one scalar calibration variable x. The parameterization
  ! determines how x is applied to the physical parameter field:
  !
  !   Physical parameter:
  !     P(i) = x
  !
  !   Spatial multiplier:
  !     P(i) = x * P_init(i)
  !
  !   Parameter ratio:
  !     P(i) = x * P_reference(i)
  !
  ! For a parameter ratio, P_reference contains the EFFECTIVE physical values
  ! of the reference parameter, including any changes made during calibration.
  !
  ! The sampled_values array contains the untransformed scalar calibration
  ! variables, ordered according to search%param_names. The conversion from
  ! transformed search coordinates has already been performed by the search.
  !
  ! Non-sampled parameters retain their initialized physical fields. They may
  ! participate in ordered constraints or serve as fixed reference parameters.
  !
  ! Parameter ratios are evaluated after their references, so references may
  ! themselves be calibrated. Missing references, incompatible field lengths,
  ! and cyclic dependencies are rejected.
  !
  ! Ordered constraints are handled separately. The model-specific adapter
  ! must also verify that related vectors represent corresponding spatial
  ! elements, since equal vector lengths alone do not establish alignment.
  ! **************************************************************************************************
  subroutine apply_parameter_sample(spec,search,sampled_values,field,err,message)
    implicit none

    type(parameter_spec),        intent(in)    :: spec
    type(parameter_search_info), intent(in)    :: search
    real(rkind),                 intent(in)    :: sampled_values(:)
    type(parameter_field),       intent(inout) :: field(:)
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message

    integer(i4b) :: i,j,k
    integer(i4b) :: nParam
    integer(i4b) :: ixSearch
    integer(i4b) :: iPass
    integer(i4b), allocatable :: reference_index(:)
    logical(lgt), allocatable :: complete(:)
    logical(lgt) :: changed
    real(rkind) :: x

    err=0
    message='apply_parameter_sample/'

    ! -----------------------------------------------------------------------------------------------
    ! Validate the parameter specification and search vector.
    ! -----------------------------------------------------------------------------------------------
    if(.not.allocated(spec%params))then
      message=trim(message)//'parameter specification is not allocated'
      err=20; return
    endif

    if(.not.allocated(search%search_index) .or. .not.allocated(search%param_names))then
      message=trim(message)//'parameter search information is not initialized'
      err=20; return
    endif

    nParam=size(spec%params)

    if(size(field)/=nParam .or. size(search%search_index)/=nParam .or. &
       size(sampled_values)/=size(search%param_names))then
      message=trim(message)//'inconsistent parameter registry or search-vector sizes'
      err=20; return
    endif

    allocate(reference_index(nParam),complete(nParam),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate dependency workspace'
      return
    endif

    reference_index=0

    ! -----------------------------------------------------------------------------------------------
    ! Verify that all initialized parameter fields are available.
    ! -----------------------------------------------------------------------------------------------
    do i=1,nParam

      if(.not.allocated(field(i)%initial))then
        message=trim(message)//'initial field not allocated for '//trim(spec%params(i)%name)
        err=20; return
      endif

      if(size(field(i)%initial)==0)then
        message=trim(message)//'initial field is empty for '//trim(spec%params(i)%name)
        err=20; return
      endif

    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Validate the sampled variables and resolve parameter-ratio references.
    ! -----------------------------------------------------------------------------------------------
    do i=1,nParam

      if(.not.spec%params(i)%sampled) cycle

      ixSearch=search%search_index(i)

      if(ixSearch<1 .or. ixSearch>size(sampled_values))then
        message=trim(message)//'invalid search index for '//trim(spec%params(i)%name)
        err=20; return
      endif

      if(trim(search%param_names(ixSearch))/=trim(spec%params(i)%name))then
        message=trim(message)//'search-vector mismatch for '//trim(spec%params(i)%name)
        err=20; return
      endif

      if(.not.ieee_is_finite(sampled_values(ixSearch)))then
        message=trim(message)//'non-finite calibration variable for '//trim(spec%params(i)%name)
        err=20; return
      endif

      select case(spec%params(i)%param_type)

      case(param_physical,param_spatial_multiplier)

        ! These parameterizations do not require reference parameters.

      case(param_ratio)

        if(len_trim(spec%params(i)%reference_name)==0)then
          message=trim(message)//'missing ratio reference for '//trim(spec%params(i)%name)
          err=20; return
        endif

        ! Resolve the reference parameter by name.
        do j=1,nParam
          if(trim(spec%params(j)%name)==trim(spec%params(i)%reference_name))then
            reference_index(i)=j
            exit
          endif
        enddo

        j=reference_index(i)

        if(j==0)then
          message=trim(message)//'unknown ratio reference for '//trim(spec%params(i)%name)//': '// &
                  trim(spec%params(i)%reference_name)
          err=20; return
        endif

        if(j==i)then
          message=trim(message)//'parameter cannot reference itself: '//trim(spec%params(i)%name)
          err=20; return
        endif

        ! Related physical fields must have compatible dimensions.
        if(size(field(i)%initial)/=size(field(j)%initial))then
          message=trim(message)//'incompatible field lengths for ratio '// &
                  trim(spec%params(i)%name)//' / '//trim(spec%params(j)%name)
          err=20; return
        endif

      case default

        message=trim(message)//'unsupported parameterization for '//trim(spec%params(i)%name)
        err=20; return

      end select

    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Restore initialized parameter fields before constructing this trial.
    ! -----------------------------------------------------------------------------------------------
    call reset_parameter_fields(field)

    ! Non-sampled parameters are already complete and can serve as references.
    do i=1,nParam
      complete(i)=.not.spec%params(i)%sampled
    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Construct effective parameter fields.
    !
    ! A parameter ratio is deferred until its reference parameter has been
    ! evaluated. Repeated passes support dependencies in any registry order.
    ! -----------------------------------------------------------------------------------------------
    do iPass=1,nParam

      changed=.false.

      do i=1,nParam

        if(complete(i)) cycle

        j=reference_index(i)

        ! Wait until a sampled reference parameter has been evaluated.
        if(j>0)then
          if(.not.complete(j)) cycle
        endif

        k=search%search_index(i)
        x=sampled_values(k)

        select case(spec%params(i)%param_type)

        case(param_physical)

          field(i)%effective=x

        case(param_spatial_multiplier)

          field(i)%effective=x*field(i)%initial

        case(param_ratio)

          field(i)%effective=x*field(j)%effective

        end select

        complete(i)=.true.
        changed=.true.

      enddo

      if(all(complete)) return
      if(.not.changed) exit

    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Unresolved dependencies indicate a cycle among parameter ratios.
    ! -----------------------------------------------------------------------------------------------
    do i=1,nParam

      if(.not.complete(i))then
        message=trim(message)//'cyclic parameter-ratio dependency involving '// &
                trim(spec%params(i)%name)
        err=20; return
      endif

    enddo

  end subroutine apply_parameter_sample



  ! **************************************************************************************************
  ! Determine the allowable scalar calibration interval for one parameter.
  !
  ! An ordered chain requires P_b(i) >= P_a(i) + gap at every spatial element i.
  ! The physical fields are represented as P_j(i) = x_j*S_j(i), where S_j is
  !  1                  for a physical parameter,
  !  P_j,init(i)        for a spatial multiplier, or
  !  P_ref,effective(i) for a parameter ratio.
  !
  ! Previously assigned parameters and unsampled parameters are held fixed.
  ! Other unassigned sampled parameters retain their configured scalar bounds.
  ! Bounds are propagated through all ordered chains connected to ixParam, so
  ! unassigned neighbours reserve room for subsequent samples.
  !
  ! is_assigned(j) is true only after the effective field of sampled parameter j
  ! has been constructed. Ratio references must be assigned (or unsampled) before
  ! this routine is called. Call again after each newly assigned parameter.
  !
  ! Returned bounds are necessary bounds. Because several unassigned variables
  ! may be coupled, a nonempty interval does not guarantee every choice can be
  ! completed; validate all physical fields when the complete sample is formed.
  ! **************************************************************************************************
  subroutine compute_feasible_interval(spec,field,is_assigned,ixParam,x_lower,x_upper,err,message)
    implicit none

    type(parameter_spec), intent(in)  :: spec            ! Parameter definitions and ordered constraints
    type(parameter_field),intent(in)  :: field(:)        ! Initialized and already-assigned physical fields
    logical(lgt),         intent(in)  :: is_assigned(:)  ! Whether each sampled parameter has been assigned
    integer(i4b),         intent(in)  :: ixParam         ! Master-registry index of the requested parameter
    real(rkind),          intent(out) :: x_lower         ! Feasible lower bound of scalar calibration variable
    real(rkind),          intent(out) :: x_upper         ! Feasible upper bound of scalar calibration variable
    integer(i4b),         intent(out) :: err             ! Error code
    character(*),         intent(out) :: message         ! Error message

    type :: scale_vector
      real(rkind), allocatable :: values(:)
    end type scale_vector

    type(scale_vector), allocatable :: scale(:)          ! Scales for current intervals (assigned = fixed)
    type(scale_vector), allocatable :: base(:)           ! Scales for computing physical chain gaps
    real(rkind), allocatable        :: lower(:),upper(:) ! Scalar bounds for constraint propagation
    real(rkind), allocatable        :: gap(:)            ! Physical minimum gap per ordered chain
    logical(lgt), allocatable       :: active(:)         ! Parameters connected to the requested one

    integer(i4b) :: nParam,nConstraint
    integer(i4b) :: i,j,k,a,b,ref,ic,ipass,iloc,nloc,maxpass
    real(rkind)  :: lo0,hi0,range_low,range_high
    real(rkind)  :: sa,sb,min_a,max_b,candidate,tol,tolx
    logical(lgt) :: changed,has_constraint

    err=0
    message='compute_feasible_interval/'
    x_lower=0._rkind
    x_upper=0._rkind

    if(.not.allocated(spec%params))then
      message=trim(message)//'parameter specification is not allocated'
      err=20; return
    endif

    nParam=size(spec%params)
    if(size(field)/=nParam .or. size(is_assigned)/=nParam)then
      message=trim(message)//'inconsistent field and assignment-array sizes'
      err=20; return
    endif

    if(ixParam<1 .or. ixParam>nParam)then
      message=trim(message)//'invalid parameter index'
      err=20; return
    endif

    if(.not.spec%params(ixParam)%sampled .or. is_assigned(ixParam))then
      message=trim(message)//'requested parameter must be sampled and not yet assigned'
      err=20; return
    endif

    x_lower=spec%params(ixParam)%lower
    x_upper=spec%params(ixParam)%upper
    if(.not.ieee_is_finite(x_lower) .or. .not.ieee_is_finite(x_upper))then
      message=trim(message)//'non-finite scalar bounds'
      err=20; return
    endif
    if(x_lower>x_upper)then
      message=trim(message)//'inverted scalar bounds'
      err=20; return
    endif
    if(.not.allocated(spec%ordered)) return

    nConstraint=size(spec%ordered)
    if(nConstraint==0) return

    ! Validate the ordered chains before indexing the parameter registry.
    has_constraint=.false.
    do ic=1,nConstraint
      if(.not.allocated(spec%ordered(ic)%param_index))then
        message=trim(message)//'unallocated ordered-constraint indices'
        err=20; return
      endif
      if(size(spec%ordered(ic)%param_index)<2)then
        message=trim(message)//'ordered constraint must contain at least two parameters'
        err=20; return
      endif
      if(.not.ieee_is_finite(spec%ordered(ic)%gap_fraction))then
        message=trim(message)//'non-finite gap_fraction'
        err=20; return
      endif
      if(spec%ordered(ic)%gap_fraction<0._rkind)then
        message=trim(message)//'negative gap_fraction'
        err=20; return
      endif
      do k=1,size(spec%ordered(ic)%param_index)
        j=spec%ordered(ic)%param_index(k)
        if(j<1 .or. j>nParam)then
          message=trim(message)//'invalid parameter index in ordered constraint'
          err=20; return
        endif
        if(j==ixParam) has_constraint=.true.
      enddo
    enddo
    if(.not.has_constraint) return

    allocate(active(nParam),lower(nParam),upper(nParam),scale(nParam),base(nParam),gap(nConstraint))
    active=.false.
    active(ixParam)=.true.

    ! Include all ordered chains transitively connected through shared parameters.
    do
      changed=.false.
      do ic=1,nConstraint
        if(.not.any(active(spec%ordered(ic)%param_index))) cycle
        do k=1,size(spec%ordered(ic)%param_index)
          j=spec%ordered(ic)%param_index(k)
          if(active(j)) cycle
          active(j)=.true.
          changed=.true.
        enddo
      enddo
      if(.not.changed) exit
    enddo

    ! Build the physical scaling field for every parameter in the connected chains.
    do j=1,nParam
      if(.not.active(j)) cycle
      if(.not.allocated(field(j)%initial))then
        message=trim(message)//'missing initialized field for '//trim(spec%params(j)%name)
        err=20; return
      endif
      if(size(field(j)%initial)==0)then
        message=trim(message)//'empty initialized field for '//trim(spec%params(j)%name)
        err=20; return
      endif
      if(any(.not.ieee_is_finite(field(j)%initial)))then
        message=trim(message)//'non-finite initialized field for '//trim(spec%params(j)%name)
        err=20; return
      endif

      ! The baseline scale is independent of whether this parameter is assigned.
      ! This keeps the gap based on configured physical bounds, not sampled values.
      if(.not.spec%params(j)%sampled)then
        base(j)%values=field(j)%initial
        lo0=1._rkind
        hi0=1._rkind
      else
        lo0=spec%params(j)%lower
        hi0=spec%params(j)%upper
        if(.not.ieee_is_finite(lo0) .or. .not.ieee_is_finite(hi0))then
          message=trim(message)//'non-finite calibration bounds for '//trim(spec%params(j)%name)
          err=20; return
        endif
        if(lo0>hi0)then
          message=trim(message)//'inverted calibration bounds for '//trim(spec%params(j)%name)
          err=20; return
        endif
        select case(spec%params(j)%param_type)
        case (param_physical)
          allocate(base(j)%values(size(field(j)%initial)))
          base(j)%values=1._rkind
        case (param_spatial_multiplier)
          base(j)%values=field(j)%initial
        case (param_ratio)
          ref=0
          do i=1,nParam
            if(trim(spec%params(i)%name)==trim(spec%params(j)%reference_name))then
              ref=i
              exit
            endif
          enddo
          if(ref==0 .or. ref==j)then
            message=trim(message)//'invalid ratio reference for '//trim(spec%params(j)%name)
            err=20; return
          endif
          if(spec%params(ref)%sampled)then
            if(.not.is_assigned(ref))then
              message=trim(message)//'assign ratio reference first: '//trim(spec%params(ref)%name)
              err=20; return
            endif
            if(.not.allocated(field(ref)%effective))then
              message=trim(message)//'missing effective reference field for '//trim(spec%params(j)%name)
              err=20; return
            endif
            base(j)%values=field(ref)%effective
          else
            if(.not.allocated(field(ref)%initial))then
              message=trim(message)//'missing initial reference field for '//trim(spec%params(j)%name)
              err=20; return
            endif
            base(j)%values=field(ref)%initial
          endif
        case default
          message=trim(message)//'unsupported parameterization for '//trim(spec%params(j)%name)
          err=20; return
        end select
      endif
      if(size(base(j)%values)/=size(field(j)%initial))then
        message=trim(message)//'misaligned parameter scale for '//trim(spec%params(j)%name)
        err=20; return
      endif
      if(any(.not.ieee_is_finite(base(j)%values)))then
        message=trim(message)//'non-finite parameter scale for '//trim(spec%params(j)%name)
        err=20; return
      endif

      if(.not.spec%params(j)%sampled)then
        scale(j)%values=field(j)%initial
        lower(j)=1._rkind
        upper(j)=1._rkind
      else if(is_assigned(j))then
        if(.not.allocated(field(j)%effective))then
          message=trim(message)//'missing assigned effective field for '//trim(spec%params(j)%name)
          err=20; return
        endif
        if(size(field(j)%effective)/=size(field(j)%initial))then
          message=trim(message)//'misaligned effective field for '//trim(spec%params(j)%name)
          err=20; return
        endif
        scale(j)%values=field(j)%effective
        lower(j)=1._rkind
        upper(j)=1._rkind
      else
        scale(j)%values=base(j)%values
        lower(j)=lo0
        upper(j)=hi0
      endif
      if(any(.not.ieee_is_finite(scale(j)%values)))then
        message=trim(message)//'non-finite effective field for '//trim(spec%params(j)%name)
        err=20; return
      endif
    enddo

    ! A chain's gap is computed from the configured PHYSICAL endpoint ranges.
    do ic=1,nConstraint
      if(.not.any(active(spec%ordered(ic)%param_index))) cycle
      k=size(spec%ordered(ic)%param_index)
      a=spec%ordered(ic)%param_index(1)
      b=spec%ordered(ic)%param_index(k)
      lo0=1._rkind
      hi0=1._rkind
      if(spec%params(a)%sampled)then
        lo0=spec%params(a)%lower
        hi0=spec%params(a)%upper
      endif
      range_low=minval(min(lo0*base(a)%values,hi0*base(a)%values))
      if(spec%params(b)%sampled)then
        lo0=spec%params(b)%lower
        hi0=spec%params(b)%upper
      else
        lo0=1._rkind
        hi0=1._rkind
      endif
      range_high=maxval(max(lo0*base(b)%values,hi0*base(b)%values))
      if(.not.ieee_is_finite(range_low) .or. .not.ieee_is_finite(range_high))then
        message=trim(message)//'non-finite physical bounds in ordered constraint'
        err=20; return
      endif
      if(range_high<range_low)then
        message=trim(message)//'inverted physical range in ordered constraint'
        err=20; return
      endif
      gap(ic)=spec%ordered(ic)%gap_fraction*(range_high-range_low)
      if(.not.ieee_is_finite(gap(ic)))then
        message=trim(message)//'non-finite physical gap'
        err=20; return
      endif
      do k=2,size(spec%ordered(ic)%param_index)
        a=spec%ordered(ic)%param_index(k-1)
        b=spec%ordered(ic)%param_index(k)
        if(size(scale(a)%values)/=size(scale(b)%values))then
          message=trim(message)//'different field lengths in ordered constraint'
          err=20; return
        endif
      enddo
    enddo

    ! Propagate bounds on the scalar calibration variables through every chain.
    maxpass=max(30,10*nParam)
    do ipass=1,maxpass
      changed=.false.
      do ic=1,nConstraint
        if(.not.any(active(spec%ordered(ic)%param_index))) cycle
        do k=2,size(spec%ordered(ic)%param_index)
          a=spec%ordered(ic)%param_index(k-1)
          b=spec%ordered(ic)%param_index(k)
          nloc=size(scale(a)%values)
          do iloc=1,nloc
            sa=scale(a)%values(iloc)
            sb=scale(b)%values(iloc)
            min_a=min(lower(a)*sa,upper(a)*sa)
            max_b=max(lower(b)*sb,upper(b)*sb)
            tol=128._rkind*epsilon(1._rkind)*max(tiny(1._rkind),abs(min_a),abs(max_b),abs(gap(ic)))

            ! B - A >= gap: bound B using the smallest attainable A.
            if(sb>0._rkind)then
              candidate=(min_a+gap(ic))/sb
              if(.not.ieee_is_finite(candidate))then
                message=trim(message)//'non-finite propagated lower bound'
                err=20; return
              endif
              if(candidate>lower(b))then
                tolx=128._rkind*epsilon(1._rkind)*max(tiny(1._rkind),abs(candidate),abs(lower(b)))
                if(candidate-lower(b)>tolx) changed=.true.
                lower(b)=candidate
              endif
            else if(sb<0._rkind)then
              candidate=(min_a+gap(ic))/sb
              if(.not.ieee_is_finite(candidate))then
                message=trim(message)//'non-finite propagated upper bound'
                err=20; return
              endif
              if(candidate<upper(b))then
                tolx=128._rkind*epsilon(1._rkind)*max(tiny(1._rkind),abs(candidate),abs(upper(b)))
                if(upper(b)-candidate>tolx) changed=.true.
                upper(b)=candidate
              endif
            else if(min_a+gap(ic)>tol)then
              message=trim(message)//'zero scale makes ordered constraint infeasible'
              err=20; return
            endif

            ! B - A >= gap: bound A using the largest attainable B.
            if(sa>0._rkind)then
              candidate=(max_b-gap(ic))/sa
              if(.not.ieee_is_finite(candidate))then
                message=trim(message)//'non-finite propagated upper bound'
                err=20; return
              endif
              if(candidate<upper(a))then
                tolx=128._rkind*epsilon(1._rkind)*max(tiny(1._rkind),abs(candidate),abs(upper(a)))
                if(upper(a)-candidate>tolx) changed=.true.
                upper(a)=candidate
              endif
            else if(sa<0._rkind)then
              candidate=(max_b-gap(ic))/sa
              if(.not.ieee_is_finite(candidate))then
                message=trim(message)//'non-finite propagated lower bound'
                err=20; return
              endif
              if(candidate>lower(a))then
                tolx=128._rkind*epsilon(1._rkind)*max(tiny(1._rkind),abs(candidate),abs(lower(a)))
                if(candidate-lower(a)>tolx) changed=.true.
                lower(a)=candidate
              endif
            else if(max_b-gap(ic)<-tol)then
              message=trim(message)//'zero scale makes ordered constraint infeasible'
              err=20; return
            endif

          enddo
        enddo
      enddo

      do j=1,nParam
        if(.not.active(j)) cycle
        tol=128._rkind*epsilon(1._rkind)*max(tiny(1._rkind),abs(lower(j)),abs(upper(j)))
        if(lower(j)>upper(j)+tol)then
          message=trim(message)//'infeasible ordered constraints for '//trim(spec%params(j)%name)
          err=20; return
        endif
        if(lower(j)>upper(j))then
          lower(j)=0.5_rkind*(lower(j)+upper(j))
          upper(j)=lower(j)
        endif
      enddo

      if(.not.changed) exit
    enddo

    if(changed)then
      message=trim(message)//'constraint bound propagation did not converge'
      err=20; return
    endif

    x_lower=lower(ixParam)
    x_upper=upper(ixParam)

  end subroutine compute_feasible_interval

  ! **************************************************************************************************
  ! Check ordered constraints requiring physical fields rather than scalar values.
  !
  ! Call after apply_parameter_sample and before writing fields into the model.
  ! Scalar-only ordered constraints remain the responsibility of parameter_search.
  ! **************************************************************************************************
  subroutine validate_effective_ordered_constraints(spec,field,only_constraint,err,message)
    implicit none

    type(parameter_spec),   intent(in)  :: spec
    type(parameter_field),  intent(in)  :: field(:)
    integer(i4b), optional, intent(in)  :: only_constraint
    integer(i4b),           intent(out) :: err
    character(*),           intent(out) :: message

    integer(i4b) :: ic,j,k,a,b,nOrdered
    real(rkind)  :: first_lower,last_upper,unused,gap,tol
    logical(lgt) :: field_based
    character(len=1024) :: cmessage

    err=0
    message='validate_effective_ordered_constraints/'

    if(.not.allocated(spec%ordered)) return

    if(present(only_constraint))then
      if(only_constraint<1 .or. only_constraint>size(spec%ordered))then
        message=trim(message)//'invalid ordered constraint index'
        err=20; return
      endif
    endif

    if(.not.allocated(spec%params))then
      message=trim(message)//'parameter specification is not allocated'
      err=20; return
    endif

    if(size(field)/=size(spec%params))then
      message=trim(message)//'parameter field size does not match parameter registry'
      err=20; return
    endif

    do ic=1,size(spec%ordered)

      if(present(only_constraint))then
        if(ic/=only_constraint) cycle
      endif

      if(.not.allocated(spec%ordered(ic)%param_index))then
        message=trim(message)//'ordered constraint has no parameter indices'
        err=20; return
      endif

      nOrdered=size(spec%ordered(ic)%param_index)

      if(nOrdered<2)then
        message=trim(message)//'ordered constraint must contain two or more parameters'
        err=20; return
      endif

      ! Check whether this constraint needs physical-field validation.
      field_based=.false.

      do j=1,nOrdered
        k=spec%ordered(ic)%param_index(j)

        if(k<1 .or. k>size(field))then
          message=trim(message)//'invalid ordered-constraint parameter index'
          err=20; return
        endif

        if(.not.spec%params(k)%sampled .or. &
          spec%params(k)%param_type/=param_physical) field_based=.true.
      enddo

      if(.not.field_based) cycle

      ! ---------------------------------------------------------------------------------------------
      ! Compute a consistent minimum gap in physical units.
      ! ---------------------------------------------------------------------------------------------
      a=spec%ordered(ic)%param_index(1)
      b=spec%ordered(ic)%param_index(nOrdered)

      call physical_parameter_limits(spec,field,a,first_lower,unused,err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif

      call physical_parameter_limits(spec,field,b,unused,last_upper,err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif

      if(.not.ieee_is_finite(spec%ordered(ic)%gap_fraction) .or. &
         spec%ordered(ic)%gap_fraction<0._rkind .or. last_upper<first_lower)then
        message=trim(message)//'invalid physical range or gap_fraction in ordered constraint'
        err=20; return
      endif

      gap=spec%ordered(ic)%gap_fraction*(last_upper-first_lower)

      if(.not.ieee_is_finite(gap))then
        message=trim(message)//'non-finite physical gap in ordered constraint'
        err=20; return
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Enforce each adjacent inequality at every spatial element.
      ! ---------------------------------------------------------------------------------------------
      do j=2,nOrdered

        a=spec%ordered(ic)%param_index(j-1)
        b=spec%ordered(ic)%param_index(j)

        if(.not.allocated(field(a)%effective) .or. .not.allocated(field(b)%effective))then
          message=trim(message)//'missing effective field in ordered constraint'
          err=20; return
        endif

        if(size(field(a)%effective)/=size(field(b)%effective))then
          message=trim(message)//'different field lengths in ordered constraint'
          err=20; return
        endif

        if(any(.not.ieee_is_finite(field(a)%effective)) .or. &
           any(.not.ieee_is_finite(field(b)%effective)))then
          message=trim(message)//'non-finite field values in ordered constraint'
          err=20; return
        endif

        do k=1,size(field(a)%effective)

          ! Scale numerical tolerance to the physical parameter magnitudes.
          tol=128._rkind*epsilon(1._rkind)* &
              max(tiny(1._rkind), &
                  abs(field(a)%effective(k)), &
                  abs(field(b)%effective(k)), &
                  abs(gap))

          if(field(b)%effective(k)-field(a)%effective(k)<gap-tol)then

            write(cmessage,'(A,I0,3(A,ES12.5))') &
                  ' at spatial element ',k, &
                  ', lower=',field(a)%effective(k), &
                  ', upper=',field(b)%effective(k), &
                  ', gap=',gap

            message=trim(message)//'violated ordered constraint: '// &
                    trim(spec%params(a)%name)//' < '//trim(spec%params(b)%name)//trim(cmessage)

            err=err_ordered_violation; return
          endif

        enddo

      enddo

    enddo

  end subroutine validate_effective_ordered_constraints
  
  ! **************************************************************************************************
  ! Return the attainable physical range of a parameter under its scalar bounds.
  ! Used to convert gap_fraction into a gap with physical parameter units.
  ! **************************************************************************************************
  subroutine physical_parameter_limits(spec,field,ix,physical_lower,physical_upper,err,message)
    implicit none

    type(parameter_spec),  intent(in)  :: spec
    type(parameter_field), intent(in)  :: field(:)
    integer(i4b),          intent(in)  :: ix
    real(rkind),           intent(out) :: physical_lower,physical_upper
    integer(i4b),          intent(out) :: err
    character(*),          intent(out) :: message

    real(rkind), allocatable :: scale(:)
    integer(i4b) :: j,ixReference
    real(rkind)  :: lower,upper

    err=0
    message='physical_parameter_limits/'
    physical_lower=0._rkind
    physical_upper=0._rkind

    if(.not.allocated(field(ix)%initial))then
      message=trim(message)//'missing initialized field for '//trim(spec%params(ix)%name)
      err=20; return
    endif

    if(size(field(ix)%initial)==0)then
      message=trim(message)//'empty initialized field for '//trim(spec%params(ix)%name)
      err=20; return
    endif

    ! A non-sampled parameter retains its initialized field.
    if(.not.spec%params(ix)%sampled)then

      if(any(.not.ieee_is_finite(field(ix)%initial)))then
        message=trim(message)//'non-finite initialized field'
        err=20; return
      endif

      physical_lower=minval(field(ix)%initial)
      physical_upper=maxval(field(ix)%initial)
      return
    endif

    lower=spec%params(ix)%lower
    upper=spec%params(ix)%upper

    if(.not.ieee_is_finite(lower) .or. .not.ieee_is_finite(upper) .or. lower>upper)then
      message=trim(message)//'invalid scalar bounds for '//trim(spec%params(ix)%name)
      err=20; return
    endif

    select case(spec%params(ix)%param_type)

    case(param_physical)
      physical_lower=lower
      physical_upper=upper
      return

    case(param_spatial_multiplier)
      scale=field(ix)%initial

    case(param_ratio)

      ixReference=0
      do j=1,size(spec%params)
        if(trim(spec%params(j)%name)==trim(spec%params(ix)%reference_name))then
          ixReference=j
          exit
        endif
      enddo

      if(ixReference==0 .or. ixReference==ix)then
        message=trim(message)//'invalid ratio reference for '//trim(spec%params(ix)%name)
        err=20; return
      endif

      if(.not.allocated(field(ixReference)%effective))then
        message=trim(message)//'reference field is not evaluated for '//trim(spec%params(ix)%name)
        err=20; return
      endif

      if(size(field(ixReference)%effective)/=size(field(ix)%initial))then
        message=trim(message)//'misaligned ratio reference for '//trim(spec%params(ix)%name)
        err=20; return
      endif

      scale=field(ixReference)%effective

    case default
      message=trim(message)//'unsupported parameterization for '//trim(spec%params(ix)%name)
      err=20; return

    end select

    if(any(.not.ieee_is_finite(scale)))then
      message=trim(message)//'non-finite parameter scale for '//trim(spec%params(ix)%name)
      err=20; return
    endif

    ! Allow positive, zero, or negative spatial scaling factors.
    physical_lower=minval(min(lower*scale,upper*scale))
    physical_upper=maxval(max(lower*scale,upper*scale))

    if(.not.ieee_is_finite(physical_lower) .or. .not.ieee_is_finite(physical_upper))then
      message=trim(message)//'non-finite physical range for '//trim(spec%params(ix)%name)
      err=20; return
    endif

  end subroutine physical_parameter_limits

end module parameter_fields
