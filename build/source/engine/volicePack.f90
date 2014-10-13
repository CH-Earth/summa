module volicePack_module
! numerical recipes data types
USE nrtype
! named variables for snow and soil
USE data_struc,only:ix_soil,ix_snow
! physical constants
USE multiconst,only:&
                    Tfreeze,  & ! freezing point              (K)
                    LH_fus,   & ! latent heat of fusion       (J kg-1)
                    LH_vap,   & ! latent heat of vaporization (J kg-1)
                    LH_sub,   & ! latent heat of sublimation  (J kg-1)
                    iden_air, & ! intrinsic density of air    (kg m-3)
                    iden_ice, & ! intrinsic density of ice    (kg m-3)
                    iden_water  ! intrinsic density of water  (kg m-3)
! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,   & ! number of snow layers  
                    nSoil,   & ! number of soil layers  
                    nLayers    ! total number of layers
implicit none
private
public::volicePack

contains

 ! ************************************************************************************************
 ! new subroutine: account for changes in volumetric ice content in the snowpack
 ! (add snowfall; combine and sub-divide layers if necessary)
 ! ************************************************************************************************
 subroutine volicePack(&
                       ! input: model control
                       dt,                          & ! intent(in): time step (seconds)
                       ! input/output: model data structures
                       mpar_data,                   & ! intent(in):    model parameters
                       indx_data,                   & ! intent(inout): type of each layer
                       mvar_data,                   & ! intent(inout): model variables for a local HRU
                       model_decisions,             & ! intent(in):    model decisions 
                       ! output: error control
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_d,            & ! data vector (dp)
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure
 ! external subroutine
 USE var_derive_module,only:calcHeight   ! module to calculate height at layer interfaces and layer mid-point
 USE layerMerge_module,only:layerMerge   ! merge snow layers if they are too thin
 USE layerDivide_module,only:layerDivide ! sub-divide layers if they are too thick
 implicit none
 ! ------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                  ! time step (seconds)
 ! input/output: model data structures
 type(var_d),intent(in)          :: mpar_data           ! model parameters
 type(var_ilength),intent(inout) :: indx_data           ! type of each layer
 type(var_dlength),intent(inout) :: mvar_data           ! model variables for a local HRU
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 ! output: error control
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! ------------------------------------------------------------------------------------------------
 ! variables in the data structures
 ! input: model parameters
 real(dp)                        :: snowfrz_scale              ! freeezing curve parameter for snow (K-1)
 ! input: diagnostic scalar variables
 real(dp)                        :: scalarSnowfallTemp         ! computed temperature of fresh snow (K)
 real(dp)                        :: scalarNewSnowDensity       ! computed density of new snow (kg m-3)
 real(dp)                        :: scalarThroughfallSnow      ! throughfall of snow through the canopy (kg m-2 s-1)
 real(dp)                        :: scalarCanopySnowUnloading  ! unloading of snow from the canopy (kg m-2 s-1)
 ! input/output: scalar state variables
 real(dp)                        :: scalarSWE                  ! SWE (kg m-2)
 real(dp)                        :: scalarSnowDepth            ! total snow depth (m)
 ! input/output: state variables
 real(dp),pointer                :: mLayerTemp(:)              ! temperature of each snow layer (K)
 real(dp),pointer                :: mLayerDepth(:)             ! depth of each snow layer (m) 
 real(dp),pointer                :: mLayerVolFracIce(:)        ! volumetric fraction of ice in each snow layer (-)
 real(dp),pointer                :: mLayerVolFracLiq(:)        ! volumetric fraction of liquid water in each snow layer (-)
 ! ------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 real(dp)                        :: volSub              ! volumetric sublimation (kg m-3)
 real(dp),parameter              :: verySmall=tiny(1._dp) ! a very small number
 ! initialize error control
 err=0; message='volicePack/'

 ! associate variables in the data structures
 associate(&

 ! input: model parameters
 snowfrz_scale             => mpar_data%var(iLookPARAM%snowfrz_scale),                   & ! freeezing curve parameter for snow (K-1)

 ! input: diagnostic scalar variables
 scalarSnowfallTemp        => mvar_data%var(iLookMVAR%scalarSnowfallTemp)%dat(1),        & ! computed temperature of fresh snow (K) 
 scalarNewSnowDensity      => mvar_data%var(iLookMVAR%scalarNewSnowDensity)%dat(1),      & ! computed density of new snow (kg m-3)
 scalarThroughfallSnow     => mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1),     & ! throughfall of snow through the canopy (kg m-2 s-1)
 scalarCanopySnowUnloading => mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1), & ! unloading of snow from the canopy (kg m-2 s-1)

 ! input/output: state variables
 scalarSWE                 => mvar_data%var(iLookMVAR%scalarSWE)%dat(1),                 & ! SWE (kg m-2)
 scalarSnowDepth           => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)            & ! total snow depth (m)
 
 ) ! end association to variables in the data structures

 ! point to variables in the data structures
 ! NOTE: use pointers because the dimension length changes
 mLayerTemp                => mvar_data%var(iLookMVAR%mLayerTemp)%dat                      ! temperature of each snow layer (K)
 mLayerDepth               => mvar_data%var(iLookMVAR%mLayerDepth)%dat                     ! depth of each snow layer (m)
 mLayerVolFracIce          => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat                ! volumetric fraction of ice in each snow layer (-)
 mLayerVolFracLiq          => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat                ! volumetric fraction of liquid water in each snow layer (-)

 ! *****
 ! * add new snowfall...
 ! *********************
 call newsnwfall(&
                 ! input: model control
                 dt,                        & ! time step (seconds)
                 snowfrz_scale,             & ! freeezing curve parameter for snow (K-1)
                 ! input: diagnostic scalar variables
                 scalarSnowfallTemp,        & ! computed temperature of fresh snow (K) 
                 scalarNewSnowDensity,      & ! computed density of new snow (kg m-3)
                 scalarThroughfallSnow,     & ! throughfall of snow through the canopy (kg m-2 s-1)
                 scalarCanopySnowUnloading, & ! unloading of snow from the canopy (kg m-2 s-1)
                 ! input/output: state variables
                 scalarSWE,                 & ! SWE (kg m-2)
                 scalarSnowDepth,           & ! total snow depth (m)
                 mLayerTemp(1),             & ! temperature of the top snow or soil layer (K)
                 mLayerDepth(1),            & ! depth of the top snow or soil layer (m)
                 mLayerVolFracIce(1),       & ! volumetric fraction of ice in the top snow or soil layer (-)
                 mLayerVolFracLiq(1),       & ! volumetric fraction of liquid water in the top snow or soil layer (-)
                 ! output: error control
                 err,cmessage)                ! error control
 if(err/=0)then; err=30; message=trim(message)//trim(cmessage); return; endif

 ! re-compute snow depth and SWE
 if(nSnow > 0)then
  scalarSnowDepth = sum(mLayerDepth(1:nSnow))
  scalarSWE       = sum( (mLayerVolFracLiq(1:nSnow)*iden_water + mLayerVolFracIce(1:nSnow)*iden_ice)*mLayerDepth(1:nSnow) )
 endif

 ! update coordinate variables
 call calcHeight(&
                 ! input/output: data structures
                 indx_data,   & ! intent(in): layer type
                 mvar_data,   & ! intent(inout): model variables for a local HRU
                 ! output: error control
                 err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! * combine and sub-divide layers, if necessary...
 ! ************************************************

 ! divide snow layers if too thick
 call layerDivide(&
                  ! input/output: model data structures
                  model_decisions,             & ! intent(in):    model decisions
                  mpar_data,                   & ! intent(in):    model parameters
                  indx_data,                   & ! intent(inout): type of each layer
                  mvar_data,                   & ! intent(inout): model variables for a local HRU
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=65; message=trim(message)//trim(cmessage); return; endif

 ! merge snow layers if they are too thin
 call layerMerge(&
                 ! input/output: model data structures
                 model_decisions,             & ! intent(in):    model decisions
                 mpar_data,                   & ! intent(in):    model parameters
                 indx_data,                   & ! intent(inout): type of each layer
                 mvar_data,                   & ! intent(inout): model variables for a local HRU
                 ! output: error control
                 err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=65; message=trim(message)//trim(cmessage); return; endif

 ! recompute the number of snow and soil layers
 nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
 nSoil = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)

 ! re-compute snow depth and SWE
 if(nSnow > 0)then
  mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
  mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                          mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                          * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
 endif

 !write(*,'(a,1x,i4,f20.10)') ' after combine; mvar_data%var(iLookMVAR%scalarSWE)%dat(1) = ', nSnow, mvar_data%var(iLookMVAR%scalarSWE)%dat(1)

 ! end association to variables in the data structures
 end associate


 end subroutine volicePack



 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! ***** PRIVATE SUBROUTINES ***************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************

 ! ************************************************************************************************
 ! new subroutine: add new snowfall to the system
 ! ************************************************************************************************
 subroutine newsnwfall(&
                       ! input: model control
                       dt,                        & ! time step (seconds)
                       fc_param,                  & ! freeezing curve parameter for snow (K-1)
                       ! input: diagnostic scalar variables
                       scalarSnowfallTemp,        & ! computed temperature of fresh snow (K) 
                       scalarNewSnowDensity,      & ! computed density of new snow (kg m-3)
                       scalarThroughfallSnow,     & ! throughfall of snow through the canopy (kg m-2 s-1)
                       scalarCanopySnowUnloading, & ! unloading of snow from the canopy (kg m-2 s-1)
                       ! input/output: state variables
                       scalarSWE,                 & ! SWE (kg m-2)
                       scalarSnowDepth,           & ! total snow depth (m)
                       surfaceLayerTemp,          & ! temperature of each layer (K)
                       surfaceLayerDepth,         & ! depth of each layer (m)
                       surfaceLayerVolFracIce,    & ! volumetric fraction of ice in each layer (-)
                       surfaceLayerVolFracLiq,    & ! volumetric fraction of liquid water in each layer (-)
                       ! output: error control
                       err,message                ) ! error control
 ! computational modules
 USE snow_utils_module,only:fracliquid,templiquid                  ! functions to compute temperature/liquid water
 ! add new snowfall to the system
 implicit none
 ! input: model control
 real(dp),intent(in)                 :: dt                         ! time step (seconds)
 real(dp),intent(in)                 :: fc_param                   ! freeezing curve parameter for snow (K-1)
 ! input: diagnostic scalar variables
 real(dp),intent(in)                 :: scalarSnowfallTemp         ! computed temperature of fresh snow (K) 
 real(dp),intent(in)                 :: scalarNewSnowDensity       ! computed density of new snow (kg m-3)
 real(dp),intent(in)                 :: scalarThroughfallSnow      ! throughfall of snow through the canopy (kg m-2 s-1)
 real(dp),intent(in)                 :: scalarCanopySnowUnloading  ! unloading of snow from the canopy (kg m-2 s-1)
 ! input/output: state variables
 real(dp),intent(inout)              :: scalarSWE                  ! SWE (kg m-2)
 real(dp),intent(inout)              :: scalarSnowDepth            ! total snow depth (m)
 real(dp),intent(inout)              :: surfaceLayerTemp           ! temperature of each layer (K)
 real(dp),intent(inout)              :: surfaceLayerDepth          ! depth of each layer (m)
 real(dp),intent(inout)              :: surfaceLayerVolFracIce     ! volumetric fraction of ice in each layer (-)
 real(dp),intent(inout)              :: surfaceLayerVolFracLiq     ! volumetric fraction of liquid water in each layer (-)
 ! output: error control
 integer(i4b),intent(out)            :: err                        ! error code
 character(*),intent(out)            :: message                    ! error message
 ! define local variables
 real(dp)                            :: newSnowfall                ! new snowfall -- throughfall and unloading (kg m-2 s-1)
 real(dp)                            :: newSnowDepth               ! new snow depth (m)
 real(dp),parameter                  :: densityCanopySnow=200._dp  ! density of snow on the vegetation canopy (kg m-3)
 real(dp)                            :: totalMassIceSurfLayer      ! total mass of ice in the surface layer (kg m-2)
 real(dp)                            :: totalDepthSurfLayer        ! total depth of the surface layer (m)
 real(dp)                            :: volFracWater               ! volumetric fraction of total water, liquid and ice (-)
 real(dp)                            :: fracLiq                    ! fraction of liquid water (-)
 real(dp)                            :: SWE                        ! snow water equivalent after snowfall (kg m-2)
 real(dp)                            :: tempSWE0                   ! temporary SWE before snowfall, used to check mass balance (kg m-2)
 real(dp)                            :: tempSWE1                   ! temporary SWE after snowfall, used to check mass balance (kg m-2)
 real(dp)                            :: xMassBalance               ! mass balance check (kg m-2)
 real(dp),parameter                  :: verySmall=1.e-8_dp         ! a very small number -- used to check mass balance
 ! initialize error control
 err=0; message="newsnwfall/"

 ! compute the new snowfall
 newSnowfall = scalarThroughfallSnow + scalarCanopySnowUnloading

 ! early return if there is no snowfall
 if(newSnowfall < tiny(dt)) return

 ! compute depth of new snow
 newSnowDepth     = dt*(scalarThroughfallSnow/scalarNewSnowDensity + scalarCanopySnowUnloading/densityCanopySnow)  ! new snow depth (m)

 ! process special case of "snow without a layer"
 if(nSnow==0)then
  ! increment depth and water equivalent
  scalarSnowDepth = scalarSnowDepth + newSnowDepth
  scalarSWE       = scalarSWE + dt*newSnowfall

 ! add snow to the top layer (more typical case where snow layers already exist)
 else

  ! get SWE in the upper layer (used to check mass balance)
  tempSWE0 = (surfaceLayerVolFracIce*iden_ice + surfaceLayerVolFracLiq*iden_water)*surfaceLayerDepth

  ! get the total mass of liquid water and ice (kg m-2)
  totalMassIceSurfLayer  = iden_ice*surfaceLayerVolFracIce*surfaceLayerDepth + newSnowfall*dt
  ! get the total snow depth
  totalDepthSurfLayer    = surfaceLayerDepth + newSnowDepth
  ! compute the new temperature
  surfaceLayerTemp       = (surfaceLayerTemp*surfaceLayerDepth + scalarSnowfallTemp*newSnowDepth) / totalDepthSurfLayer
  ! compute new SWE for the upper layer (kg m-2)
  SWE = totalMassIceSurfLayer + iden_water*surfaceLayerVolFracLiq*surfaceLayerDepth
  ! compute new volumetric fraction of liquid water and ice (-)
  volFracWater = (SWE/totalDepthSurfLayer)/iden_water
  fracLiq      = fracliquid(surfaceLayerTemp,fc_param)                           ! fraction of liquid water
  surfaceLayerVolFracIce = (1._dp - fracLiq)*volFracWater*(iden_water/iden_ice)  ! volumetric fraction of ice (-)
  surfaceLayerVolFracLiq =          fracLiq *volFracWater                        ! volumetric fraction of liquid water (-)
  ! update new layer depth (m)
  surfaceLayerDepth      = totalDepthSurfLayer

  ! get SWE in the upper layer (used to check mass balance)
  tempSWE1 = (surfaceLayerVolFracIce*iden_ice + surfaceLayerVolFracLiq*iden_water)*surfaceLayerDepth

  ! check SWE
  xMassBalance = tempSWE1 - (tempSWE0 + newSnowfall*dt)
  if (abs(xMassBalance) > verySmall)then
   write(*,'(a,1x,f20.10)') 'SWE mass balance = ', xMassBalance
   message=trim(message)//'mass balance problem'
   err=20; return
  endif

 endif  ! if snow layers already exist

 end subroutine newsnwfall


end module volicePack_module
