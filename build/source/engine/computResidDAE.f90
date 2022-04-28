

module computResidDAE_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (rkind)

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! define access to state variables to print
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
! privacy
implicit none
private::printResidDAE
public::computResidDAE
contains

 ! **********************************************************************************************************
 ! public subroutine computResidDAE: compute the residual vector
 ! **********************************************************************************************************
 subroutine computResidDAE(&
                        ! input: model control
                        nSnow,                     & ! intent(in):    number of snow layers
                        nSoil,                     & ! intent(in):    number of soil layers
                        nLayers,                   & ! intent(in):    total number of layers
                        ! input: flux vectors
                        sMul,                      & ! intent(in):    state vector multiplier (used in the residual calculations)
                        fVec,                      & ! intent(in):    flux vector
                        ! input: state variables (already disaggregated into scalars and vectors)
                        scalarCanopyTempTrial,     & ! intent(in):
                        mLayerTempTrial,           & ! intent(in)
                        scalarCanairTempPrime,     & ! intent(in):    trial value for the temperature of the canopy air space (K)
                        scalarCanopyTempPrime,     & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                        scalarCanopyWatPrime,      & !
                        mLayerTempPrime,           & ! intent(in):    trial value for the temperature of each snow and soil layer (K) water content
                        scalarAquiferStoragePrime, & ! intent(in):    trial value of storage of water in the aquifer (m)
                        ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                        scalarCanopyIcePrime,      & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                        scalarCanopyLiqPrime,      & ! intent(in):
                        mLayerVolFracIcePrime,     & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                        mLayerVolFracWatPrime,     &
                        mLayerVolFracLiqPrime,     &
                        scalarCanopyCmTrial,       &
                        mLayerCmTrial,             & ! intent(in)
                        ! input: data structures
                        prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,                 & ! intent(in):    model diagnostic variables for a local HRU
                        flux_data,                 & ! intent(in):    model fluxes for a local HRU
                        indx_data,                 & ! intent(in):    index data
                        ! output
                        rAdd,                      & ! intent(out):   additional (sink) terms on the RHS of the state equation
                        rVec,                      & ! intent(out):   residual vector
                        err,message)                 ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 integer(i4b),intent(in)         :: nSnow                     ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                   ! total number of layers in the snow+soil domain
 ! input: flux vectors
 real(qp),intent(in)             :: sMul(:)   ! NOTE: qp      ! state vector multiplier (used in the residual calculations)
 real(rkind),intent(in)             :: fVec(:)                   ! flux vector
 ! input: state variables (already disaggregated into scalars and vectors)
  real(rkind),intent(in)            :: scalarCanopyTempTrial
  real(rkind),intent(in)            :: mLayerTempTrial(:)
 real(rkind),intent(in)             :: scalarCanairTempPrime     ! trial value for temperature of the canopy air space (K)
 real(rkind),intent(in)             :: scalarCanopyTempPrime     ! trial value for temperature of the vegetation canopy (K)
 real(rkind),intent(in)             :: scalarCanopyWatPrime      ! derivative value for liquid water storage in the canopy (kg m-2)
 real(rkind),intent(in)             :: mLayerTempPrime(:)        ! trial value for temperature of each snow/soil layer (K) content
 real(rkind),intent(in)             :: scalarAquiferStoragePrime ! trial value of aquifer storage (m)
 ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
 real(rkind),intent(in)             :: scalarCanopyIcePrime      ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(rkind),intent(in)             :: scalarCanopyLiqPrime
 real(rkind),intent(in)             :: mLayerVolFracIcePrime(:)  ! trial value for volumetric fraction of ice (-)
 real(rkind),intent(in)             :: mLayerVolFracLiqPrime(:)
 real(rkind),intent(in)             ::  mLayerVolFracWatPrime(:)
 real(qp),intent(in)             ::  scalarCanopyCmTrial
 real(qp),intent(in)             ::  mLayerCmTrial(:)
 ! input: data structures
 type(var_dlength),intent(in)    :: prog_data                 ! prognostic variables for a local HRU
 type(var_dlength),intent(in)    :: diag_data                 ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)    :: flux_data                 ! model fluxes for a local HRU
 type(var_ilength),intent(in)    :: indx_data                 ! indices defining model states and layers
 ! output
 real(rkind),intent(out)            :: rAdd(:)                   ! additional (sink) terms on the RHS of the state equation
 real(qp),intent(out)            :: rVec(:)   ! NOTE: qp      ! residual vector
 integer(i4b),intent(out)        :: err                       ! error code
 character(*),intent(out)        :: message                   ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 integer(i4b)                    :: iLayer                    ! index of layer within the snow+soil domain
 integer(i4b),parameter          :: ixVegVolume=1             ! index of the desired vegetation control volumne (currently only one veg layer)
 real(rkind),dimension(nLayers)     :: mLayerVolFracHydPrime          ! vector of volumetric water content (-), either liquid water content or total water content
  real(rkind)                       :: scalarCanopyHydPrime      ! trial value for canopy water (kg m-2), either liquid water content or total water content
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! link to the necessary variables for the residual computations
 associate(&
  ! canopy and layer depth
  canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in): [dp]      canopy depth (m)
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
  ! model fluxes (sink terms in the soil domain)
  mLayerTranspire         => flux_data%var(iLookFLUX%mLayerTranspire)%dat           ,& ! intent(in): [dp]     transpiration loss from each soil layer (m s-1)
  mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,& ! intent(in): [dp(:)]  baseflow from each soil layer (m s-1)
  mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,& ! intent(in): [dp(:)]  change in storage associated with compression of the soil matrix (-)
  ! number of state variables of a specific type
  nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
  nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
  nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
  ! model indices
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy energy state variable
  ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
  ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in): [i4b]    index of water storage in the aquifer
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
  ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
  ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
  ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
  ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in): [i4b(:)] named variables defining the type of hydrology states in snow+soil domain
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                 & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
 ) ! association to necessary variables for the residual computations
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="computResidDAE/"

 ! ---
 ! * compute sink terms...
 ! -----------------------

 ! intialize additional terms on the RHS as zero
 rAdd(:) = 0._rkind

 ! compute energy associated with melt freeze for the vegetation canopy
 if(ixVegNrg/=integerMissing) rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*scalarCanopyIcePrime/canopyDepth   ! energy associated with melt/freeze (J m-3)
 ! compute energy associated with melt/freeze for snow
 ! NOTE: allow expansion of ice during melt-freeze for snow; deny expansion of ice during melt-freeze for soil
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   select case( layerType(iLayer) )
    case(iname_snow); rAdd( ixSnowSoilNrg(iLayer) ) = rAdd( ixSnowSoilNrg(iLayer) ) + LH_fus*iden_ice * mLayerVolFracIcePrime(iLayer)
    case(iname_soil); rAdd( ixSnowSoilNrg(iLayer) ) = rAdd( ixSnowSoilNrg(iLayer) ) + LH_fus*iden_water * mLayerVolFracIcePrime(iLayer)
   end select
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! sink terms soil hydrology (-)
 ! NOTE 1: state variable is volumetric water content, so melt-freeze is not included
 ! NOTE 2: ground evaporation was already included in the flux at the upper boundary
 ! NOTE 3: rAdd(ixSnowOnlyWat)=0, and is defined in the initialization above
 ! NOTE 4: same sink terms for matric head and liquid matric potential
 if(nSoilOnlyHyd>0)then
  do concurrent (iLayer=1:nSoil,ixSoilOnlyHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
   rAdd( ixSoilOnlyHyd(iLayer) ) = rAdd( ixSoilOnlyHyd(iLayer) ) + (mLayerTranspire(iLayer) - mLayerBaseflow(iLayer) )/mLayerDepth(iLayer+nSnow) - mLayerCompress(iLayer)
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! ---
 ! * compute the residual vector...
 ! --------------------------------

 ! compute the residual vector for the vegetation canopy
 ! NOTE: sMul(ixVegHyd) = 1, but include as it converts all variables to quadruple precision
 ! --> energy balance
 if(ixCasNrg/=integerMissing) rVec(ixCasNrg) = sMul(ixCasNrg)*scalarCanairTempPrime - ( fVec(ixCasNrg) + rAdd(ixCasNrg) )
 if(ixVegNrg/=integerMissing) rVec(ixVegNrg) = sMul(ixVegNrg) * scalarCanopyTempPrime + scalarCanopyCmTrial * scalarCanopyWatPrime/canopyDepth  - ( fVec(ixVegNrg) + rAdd(ixVegNrg) )
 ! --> mass balance
 if(ixVegHyd/=integerMissing)then
  scalarCanopyHydPrime = merge(scalarCanopyWatPrime, scalarCanopyLiqPrime, (ixStateType( ixHydCanopy(ixVegVolume) )==iname_watCanopy) )
   rVec(ixVegHyd)  = sMul(ixVegHyd)*scalarCanopyHydPrime  - ( fVec(ixVegHyd) + rAdd(ixVegHyd) )
 endif

 ! compute the residual vector for the snow and soil sub-domains for energy
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   rVec( ixSnowSoilNrg(iLayer) ) = sMul( ixSnowSoilNrg(iLayer) ) * mLayerTempPrime(iLayer) + mLayerCmTrial(iLayer) * mLayerVolFracWatPrime(iLayer) - ( fVec( ixSnowSoilNrg(iLayer) ) + rAdd( ixSnowSoilNrg(iLayer) ) )
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! compute the residual vector for the snow and soil sub-domains for hydrology
 ! NOTE: residual depends on choice of state variable
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
   ! (get the correct state variable)
   mLayerVolFracHydPrime(iLayer)      = merge(mLayerVolFracWatPrime(iLayer), mLayerVolFracLiqPrime(iLayer), (ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_matLayer) )
   ! (compute the residual)
 rVec( ixSnowSoilHyd(iLayer) ) = mLayerVolFracHydPrime(iLayer) - ( fVec( ixSnowSoilHyd(iLayer) ) + rAdd( ixSnowSoilHyd(iLayer) ) )
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! compute the residual vector for the aquifer
  if(ixAqWat/=integerMissing)  rVec(ixAqWat) = sMul(ixAqWat)*scalarAquiferStoragePrime - ( fVec(ixAqWat) + rAdd(ixAqWat) )

 ! print result
 if(globalPrintFlag)then
  write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
  write(*,'(a,1x,100(e12.5,1x))') 'fVec = ', fVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
  !print*, 'PAUSE:'; read(*,*)
 endif

 !call printResidDAE(nSnow,nSoil,nLayers,indx_data,rAdd,rVec)

 ! check
 if(any(isNan(rVec)))then
  call printResidDAE(nSnow,nSoil,nLayers,indx_data,rAdd,rVec)
  message=trim(message)//'we found NaN'
  err=20; return
 endif

 ! end association with the necessary variabiles for the residual calculations
 end associate

 end subroutine computResidDAE

 ! **********************************************************************************************************
 ! private subroutine printResidDAE: print the residual vector mainly for debugging
 ! **********************************************************************************************************
 subroutine printResidDAE( &
                          ! input: model control
                          nSnow,                     & ! intent(in):    number of snow layers
                          nSoil,                     & ! intent(in):    number of soil layers
                          nLayers,                   & ! intent(in):    total number of layers
                          ! input: data structures
                          indx_data,                 & ! intent(in):    index data
                          ! output
                          rAdd,                      & ! intent(out):   additional (sink) terms on the RHS of the state equation
                          rVec)                        ! intent(out):   residual vector

! --------------------------------------------------------------------------------------------------------------------------------
implicit none
! input: model control
integer(i4b),intent(in)         :: nSnow                     ! number of snow layers
integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
integer(i4b),intent(in)         :: nLayers                   ! total number of layers in the snow+soil domain
type(var_ilength),intent(in)    :: indx_data                 ! indices defining model states and layers
! output
real(rkind),intent(in)          :: rAdd(:)                   ! additional (sink) terms on the RHS of the state equation
real(qp),intent(in)             :: rVec(:)   ! NOTE: qp      ! residual vector
! --------------------------------------------------------------------------------------------------------------------------------
! local variables
! --------------------------------------------------------------------------------------------------------------------------------
integer(i4b)                    :: iLayer                    ! index of layer within the snow+soil domain
! --------------------------------------------------------------------------------------------------------------------------------
! link to the necessary variables for the residual computations
associate(&
! number of state variables of a specific type
nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
! model indices
ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy energy state variable
ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in): [i4b]    index of water storage in the aquifer
ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in): [i4b(:)] named variables defining the type of hydrology states in snow+soil domain
layerType               => indx_data%var(iLookINDEX%layerType)%dat                 & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
) ! association to necessary variables for the residual computations
! --------------------------------------------------------------------------------------------------------------------------------

if(ixVegNrg/=integerMissing) print *, 'rAdd(ixVegNrg) = ', rAdd(ixVegNrg)

if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)
    select case( layerType(iLayer) )
      case(iname_snow)
        print *, 'rAdd( ixSnowSoilNrg(iLayer) ) = ', rAdd( ixSnowSoilNrg(iLayer) )
      case(iname_soil); print *, 'rAdd( ixSnowSoilNrg(iLayer) ) = ', rAdd( ixSnowSoilNrg(iLayer) )
    end select
  end do
endif

if(nSoilOnlyHyd>0)then
  do concurrent (iLayer=1:nSoil,ixSoilOnlyHyd(iLayer)/=integerMissing)
    print *, 'rAdd( ixSoilOnlyHyd(iLayer) ) = ', rAdd( ixSoilOnlyHyd(iLayer) )
  end do
endif

if(ixCasNrg/=integerMissing) print *, 'rVec(ixCasNrg) = ', rVec(ixCasNrg)
if(ixVegNrg/=integerMissing) print *, 'rVec(ixVegNrg) = ', rVec(ixVegNrg)
if(ixVegHyd/=integerMissing)then
  print *, 'rVec(ixVegHyd) = ', rVec(ixVegHyd)
endif

if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)
    print *, 'rVec( ixSnowSoilNrg(iLayer) ) = ', rVec( ixSnowSoilNrg(iLayer) )
  end do
endif

if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)
    print *, 'rVec( ixSnowSoilHyd(iLayer) ) = ', rVec( ixSnowSoilHyd(iLayer) )
  end do
endif

if(ixAqWat/=integerMissing)  print *, ' rVec(ixAqWat) = ', rVec(ixAqWat)

end associate

end subroutine printResidDAE

end module computResidDAE_module
