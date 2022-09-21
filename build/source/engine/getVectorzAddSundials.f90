

module getVectorzAddSundials_module

! data types
USE nrtype

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

! constants
USE multiconst,only:&
                    gravity,      & ! acceleration of gravity              (m s-2)
                    Tfreeze,      & ! temperature at freezing              (K)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

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

! provide access to routines to update states
USE updatState_module,only:updateSnow     ! update snow states
USE updatState_module,only:updateSoil     ! update soil states

! provide access to functions for the constitutive functions and derivatives
USE snow_utils_module,only:fracliquid     ! compute the fraction of liquid water (snow)
USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
USE soil_utils_module,only:dTheta_dTk     ! differentiate the freezing curve w.r.t. temperature (soil)
USE soil_utils_module,only:dTheta_dPsi    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:dPsi_dTheta    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric water content
USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
USE soil_utils_module,only:liquidHead     ! compute the liquid water matric potential

implicit none
private
public::residDiscontinuity
public::countDiscontinuity
contains



 ! **********************************************************************************************************
 ! public subroutine residDiscontinuity:
 ! **********************************************************************************************************
 subroutine residDiscontinuity(&
                       ! input
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       diag_data,                                 & ! intent(in):    model diagnostic variables for a local HRU
                       prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       ! output
                       resid,                                     & ! intent(out)
                       err,message)                                 ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input
 real(rkind),intent(in)             :: stateVec(:)                     ! model state vector (mixed units)
 type(var_dlength),intent(in)    :: diag_data                       ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)    :: prog_data                       ! prognostic variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data                       ! indices defining model states and layers
 real(qp),intent(out)            :: resid(:)
 ! output: error control
 integer(i4b),intent(out)        :: err                             ! error code
 character(*),intent(out)        :: message                         ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iLayer                          ! index of layer within the snow+soil domain
 integer(i4b)                    :: iCount
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of model layers, and layer type
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
 ! indices defining model states and layers
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
 ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in):  [i4b]    index of the squifer storage state variable
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for energy states in the snow+soil subdomain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for hydrology states in the snow+soil subdomain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):  [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):  [i4b]    number of hydrology variables in the snow+soil domain
 ! indices defining type of model state variables
 ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in):  [i4b(:)] [state subset] type of desired model state variables
 ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                & ! intent(in):  [i4b(:)] index of the type of hydrology states in snow+soil domain
)! association with variables in the data structures

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='residDiscontinuity/'

 iCount = 1


 ! check if computing the vegetation flux
 if(ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing)then

  ! temperature of the canopy air space
  if(ixCasNrg/=integerMissing)then
    resid(iCount) = stateVec(ixCasNrg) - Tfreeze    !  scalarCanairTempTrial - Tfreeze
    iCount = iCount + 1
  endif

  ! canopy temperature
  if(ixVegNrg/=integerMissing)then
    resid(iCount) = stateVec(ixVegNrg) - Tfreeze    !  scalarCanopyTempTrial - Tfreeze
    iCount = iCount + 1
  endif


 endif  ! not computing the vegetation flux

 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   resid( iCount ) = stateVec( ixSnowSoilNrg(iLayer) ) - Tfreeze          ! mLayerTempTrial(iLayer) - Tfreeze
   iCount = iCount + 1
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 end associate


 end subroutine residDiscontinuity

  ! **********************************************************************************************************
 ! public subroutine countDiscontinuity:
 ! **********************************************************************************************************
 subroutine countDiscontinuity(&
                       ! input
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       diag_data,                                 & ! intent(in):    model diagnostic variables for a local HRU
                       prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       ! output
                       countD,                                     & ! intent(out)
                       err,message)                                 ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input
 real(rkind),intent(in)             :: stateVec(:)                     ! model state vector (mixed units)
 type(var_dlength),intent(in)    :: diag_data                       ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)    :: prog_data                       ! prognostic variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data                       ! indices defining model states and layers
 integer(i4b),intent(out)        :: countD
 ! output: error control
 integer(i4b),intent(out)        :: err                             ! error code
 character(*),intent(out)        :: message                         ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iLayer                          ! index of layer within the snow+soil domain
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of model layers, and layer type
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
 ! indices defining model states and layers
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
 ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in):  [i4b]    index of the squifer storage state variable
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for energy states in the snow+soil subdomain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for hydrology states in the snow+soil subdomain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):  [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):  [i4b]    number of hydrology variables in the snow+soil domain
 ! indices defining type of model state variables
 ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in):  [i4b(:)] [state subset] type of desired model state variables
 ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                & ! intent(in):  [i4b(:)] index of the type of hydrology states in snow+soil domain
)! association with variables in the data structures

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='countDiscontinuity/'

 ! *** extract state variables for the vegetation canopy

 countD = 0
 ! check if computing the vegetation flux
 if(ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing)then

  ! temperature of the canopy air space
  if(ixCasNrg/=integerMissing)  countD = countD + 1

  ! canopy temperature
  if(ixVegNrg/=integerMissing)  countD = countD + 1


 endif  ! not computing the vegetation flux

 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   countD = countD + 1
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 end associate

 end subroutine countDiscontinuity



end module getVectorzAddSundials_module
