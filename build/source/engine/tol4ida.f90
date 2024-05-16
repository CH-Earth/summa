module tol4ida_module

!======= Inclusions ===========
use, intrinsic :: iso_c_binding
use nrtype
use type4ida

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas        ! named variables for canopy air space
USE globalData,only:iname_veg        ! named variables for vegetation canopy
USE globalData,only:iname_snow       ! named variables for snow
USE globalData,only:iname_soil       ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair  ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy  ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy  ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy  ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer   ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer   ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer   ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer   ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:iname_watAquifer ! named variable defining the water storage in the aquifer

! metadata for information in the data structures
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (rkind)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookDIAG             ! named variables for structure elements
USE var_lookup,only:iLookPROG             ! named variables for structure elements
USE var_lookup,only:iLookDERIV            ! named variables for structure elements
USE var_lookup,only:iLookPARAM            ! named variables for structure elements
USE var_lookup,only:iLookINDEX            ! named variables for structure elements


! privacy
implicit none
private
public::computWeight4ida
public::popTol4ida


contains

! **********************************************************************************************************
! public function computWeight4ida: compute w_i = 1 / ( rtol_i * y_i + atol_i )
! **********************************************************************************************************
! Return values:
!    0 = success,
!   -1 = non-recoverable error, NaN or negative values
! ----------------------------------------------------------------
integer(c_int) function computWeight4ida(sunvec_y, sunvec_ewt, user_data) &
  result(ierr) bind(C,name='computWeight4ida')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_serial_mod
  use nrtype
  use type4ida

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)          :: sunvec_y  ! solution N_Vector    y
  type(N_Vector)          :: sunvec_ewt ! derivative N_Vector W
  type(c_ptr), value      :: user_data ! user-defined data

  ! pointers to data in SUNDIALS vectors
  type(data4ida), pointer :: tol_data ! equations data
  real(rkind), pointer    :: stateVec(:)
  real(rkind), pointer    :: weightVec(:)
  integer(c_int)          :: iState

  !======= Internals ============

  ! get equations data from user-defined data
  call c_f_pointer(user_data, tol_data)

  ! get data arrays from SUNDIALS vectors
  stateVec(1:tol_data%nState)  => FN_VGetArrayPointer(sunvec_y)
  weightVec(1:tol_data%nState)  => FN_VGetArrayPointer(sunvec_ewt)

  do iState = 1,tol_data%nState
    weightVec(iState) = tol_data%rtol(iState) * abs( stateVec(iState) ) + tol_data%atol(iState)
    weightVec(iState) = 1._rkind / weightVec(iState)
  end do

  ierr = 0
  return

end function computWeight4ida


! **********************************************************************************************************
! public subroutine popTol4ida: populate tolerances for state vectors
! **********************************************************************************************************
subroutine popTol4ida(&
                      ! input: data structures
                      nState,                  & ! intent(in):    number of desired state variables
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                      indx_data,               & ! intent(in):    indices defining model states and layers
                      mpar_data,               & ! intent(in)
                      ! output
                      absTol,                  & ! intent(out):   model state vector
                      relTol,                  &
                      err,message)               ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: data structures
  integer(i4b),intent(in)         :: nState                 ! number of desired state variables
  type(var_dlength),intent(in)    :: prog_data              ! prognostic variables for a local HRU
  type(var_dlength),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
  type(var_ilength),intent(in)    :: indx_data              ! indices defining model states and layers
  type(var_dlength),intent(in)    :: mpar_data              ! model parameters
  ! output
  real(rkind),intent(out)         :: absTol(:)            ! model state vector (mixed units)
  real(rkind),intent(out)         :: relTol(:)            ! model state vector (mixed units)
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! state subsets
  integer(i4b)                       :: iState                    ! index of state within the snow+soil domain
  integer(i4b)                       :: iLayer                    ! index of layer within the snow+soil domain
  integer(i4b)                       :: ixStateSubset             ! index within the state subset
  logical(lgt),dimension(nState)     :: tolFlag                   ! flag to denote that the state is populated
  real(rkind)                        :: absTolTempCas
  real(rkind)                        :: relTolTempCas
  real(rkind)                        :: absTolTempVeg
  real(rkind)                        :: relTolTempVeg
  real(rkind)                        :: absTolWatVeg
  real(rkind)                        :: relTolWatVeg
  real(rkind)                        :: absTolTempSoilSnow
  real(rkind)                        :: relTolTempSoilSnow
  real(rkind)                        :: absTolWatSnow
  real(rkind)                        :: relTolWatSnow
  real(rkind)                        :: absTolMatric
  real(rkind)                        :: relTolMatric
  real(rkind)                        :: absTolAquifr
  real(rkind)                        :: relTolAquifr

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  fixedLength: associate(&
    scalarCanairTemp    => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in) : [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp    => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in) : [dp]     temperature of the vegetation canopy (K)
    scalarCanopyWat     => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in) : [dp]     mass of total water on the vegetation canopy (kg m-2)
    scalarCanopyLiq     => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in) : [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    ! model state variable vectors for the snow-soil layers
    mLayerTemp          => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in) : [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracWat    => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in) : [dp(:)]  volumetric fraction of total water (-)
    mLayerVolFracLiq    => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in) : [dp(:)]  volumetric fraction of liquid water (-)
    mLayerMatricHead    => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in) : [dp(:)]  matric head (m)
    mLayerMatricHeadLiq => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(in) : [dp(:)]  matric potential of liquid water (m)
    ! model state variables for the aquifer
    scalarAquiferStorage=> prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)   ,& ! intent(in) : [dp]     storage of water in the aquifer (m)
    ! indices defining specific model states
    ixCasNrg            => indx_data%var(iLookINDEX%ixCasNrg)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy air space energy state variable
    ixVegNrg            => indx_data%var(iLookINDEX%ixVegNrg)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy energy state variable
    ixVegHyd            => indx_data%var(iLookINDEX%ixVegHyd)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy hydrology state variable (mass)
    ixAqWat             => indx_data%var(iLookINDEX%ixAqWat)%dat                  ,& ! intent(in) : [i4b(:)] [length=1] index of aquifer storage state variable
    ! vector of energy and hydrology indices for the snow and soil domains
    ixSnowSoilNrg       => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowSoilHyd       => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    nSnowSoilNrg        => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd        => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
    ! type of model state variabless
    ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
    ixHydType           => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in) : [i4b(:)] index of the type of hydrology states in snow+soil domain
    ! number of layers
    nSnow               => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in) : [i4b]    number of snow layers
    nSoil               => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in) : [i4b]    number of soil layers
    nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in) : [i4b]    total number of layers
    )  ! end association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='popTol4ida/'

    absTolTempCas      = mpar_data%var(iLookPARAM%absTolTempCas)%dat(1)
    relTolTempCas      = mpar_data%var(iLookPARAM%relTolTempCas)%dat(1)
    absTolTempVeg      = mpar_data%var(iLookPARAM%absTolTempVeg)%dat(1)
    relTolTempVeg      = mpar_data%var(iLookPARAM%relTolTempVeg)%dat(1)
    absTolWatVeg       = mpar_data%var(iLookPARAM%absTolWatVeg)%dat(1)
    relTolWatVeg       = mpar_data%var(iLookPARAM%relTolWatVeg)%dat(1)
    absTolTempSoilSnow = mpar_data%var(iLookPARAM%absTolTempSoilSnow)%dat(1)
    relTolTempSoilSnow = mpar_data%var(iLookPARAM%relTolTempSoilSnow)%dat(1)
    absTolWatSnow      = mpar_data%var(iLookPARAM%absTolWatSnow)%dat(1)
    relTolWatSnow      = mpar_data%var(iLookPARAM%relTolWatSnow)%dat(1)
    absTolMatric       = mpar_data%var(iLookPARAM%absTolMatric)%dat(1)
    relTolMatric       = mpar_data%var(iLookPARAM%relTolMatric)%dat(1)
    absTolAquifr       = mpar_data%var(iLookPARAM%absTolAquifr)%dat(1)
    relTolAquifr       = mpar_data%var(iLookPARAM%relTolAquifr)%dat(1)
 
    ! -----
    ! * initialize state vectors...
    ! -----------------------------

    ! initialize flags
    tolFlag(:) = .false.

    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
    do concurrent (iState=1:size(ixCasNrg),ixCasNrg(iState)/=integerMissing)
      absTol( ixCasNrg(iState) )  = absTolTempCas            ! transfer canopy air temperature to the state vector
      relTol( ixCasNrg(iState) )  = relTolTempCas
      tolFlag( ixCasNrg(iState) ) = .true.                 ! flag to denote that tolerances are populated
    end do

    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
    do concurrent (iState=1:size(ixVegNrg),ixVegNrg(iState)/=integerMissing)
      absTol( ixVegNrg(iState) )  = absTolTempVeg      ! transfer vegetation temperature to the state vector
      relTol( ixVegNrg(iState) )  = relTolTempVeg     ! transfer vegetation temperature to the state vector
      tolFlag( ixVegNrg(iState) ) = .true.                 ! flag to denote that tolerances are populated
    end do

    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
    do concurrent (iState=1:size(ixVegHyd),ixVegHyd(iState)/=integerMissing)
      tolFlag( ixVegHyd(iState) ) = .true.                 ! flag to denote that tolerances are populated
      select case(ixStateType_subset( ixVegHyd(iState) ))
        case(iname_watCanopy); absTol( ixVegHyd(iState) )  = absTolWatVeg ; relTol( ixVegHyd(iState) )  = relTolWatVeg
        case(iname_liqCanopy); absTol( ixVegHyd(iState) )  = absTolWatVeg ; relTol( ixVegHyd(iState) )  = relTolWatVeg       ! transfer liquid canopy water to the state vector
        case default; tolFlag( ixVegHyd(iState) ) = .false. ! flag to denote that tolerances are populated
      end select
    end do

    ! tolerance for tempreture of the snow and soil domain
    if(nSnowSoilNrg>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
        ixStateSubset            = ixSnowSoilNrg(iLayer)  ! index within the state vector
        absTol(ixStateSubset)  = absTolTempSoilSnow    ! transfer temperature from a layer to the state vector
        relTol(ixStateSubset)  = relTolTempSoilSnow
        tolFlag(ixStateSubset) = .true.                 ! flag to denote that tolerances are populated
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
    if(nSnowSoilHyd>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
        ixStateSubset            = ixSnowSoilHyd(iLayer)   ! index within the state vector
        tolFlag(ixStateSubset) = .true.                  ! flag to denote that tolerances are populated
        select case( ixHydType(iLayer) )
          case(iname_watLayer); absTol(ixStateSubset) = absTolWatSnow ;  relTol(ixStateSubset) = relTolWatSnow
          case(iname_liqLayer); absTol(ixStateSubset) = absTolWatSnow ;  relTol(ixStateSubset) = relTolWatSnow
          case(iname_matLayer); absTol(ixStateSubset) = absTolMatric ;  relTol(ixStateSubset) = relTolMatric
          case(iname_lmpLayer); absTol(ixStateSubset) = absTolMatric ;  relTol(ixStateSubset) = relTolMatric
          case default; tolFlag(ixStateSubset) = .false.  ! flag to denote that tolerances are populated
        end select
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! build the state vector for the aquifer storage
    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer aquifer
    do concurrent (iState=1:size(ixAqWat),ixAqWat(iState)/=integerMissing)
      absTol( ixAqWat(iState) )  = absTolAquifr    ! transfer aquifer storage to the state vector
      relTol( ixAqWat(iState) )  = relTolAquifr
      tolFlag( ixAqWat(iState) ) = .true.                  ! flag to denote that tolerances are populated
    end do

    ! check that we specified tolerances for all state variables
    if(count(tolFlag)/=nState)then
      print*, 'tolFlag = ', tolFlag
      message=trim(message)//'tolerances not specified for some state variables'
      err=20; return
    endif

  end associate fixedLength      ! end association to variables in the data structure where vector length does not change
end subroutine popTol4ida


end module tol4ida_module
