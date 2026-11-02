
module parameter_fields

  USE nr_type, only: i4b,rkind,lgt
  USE parameter_search, only: parameter_spec,parameter_search_info
  USE parameter_search, only: param_physical,param_spatial_multiplier,param_ratio
  USE, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

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

end module parameter_fields
