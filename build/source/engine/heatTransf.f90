module heatTransf_module
USE nrtype
! physical constants
USE multiconst,only:&
                    sb,          & ! Stefan Boltzman constant      (W m-2 K-4)
                    Em_Sno,      & ! emissivity of snow            (-)
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    LH_fus,      & ! latent heat of fusion         (J kg-1)
                    LH_vap,      & ! latent heat of vaporization   (J kg-1)
                    LH_sub,      & ! latent heat of sublimation    (J kg-1)
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water     ! intrinsic density of water    (kg m-3)
! named variables for snow and soil
USE data_struc,only:ix_soil,ix_snow                        ! names variables for snow and soil
! provide access to look-up values for model decisions
USE mDecisions_module,only:      &
 ! look-up values for the numerical method
 iterative,                      & ! iterative
 nonIterative,                   & ! non-iterative
 iterSurfEnergyBal,              & ! iterate only on the surface energy balance
 ! look-up values for method used to compute derivative
 numerical,                      & ! numerical solution
 analytical,                     & ! analytical solution
 ! look-up values for choice of boundary conditions for thermodynamics
 prescribedTemp,                 & ! prescribed temperature
 energyFlux,                     & ! energy flux
 zeroFlux,                       & ! zero flux
 ! look-up values for choice of boundary conditions for soil hydrology
 prescribedHead,                 & ! prescribed head
 ! look-up values for the choice of groundwater parameterization
 equilWaterTable,                & ! equilibrium water table
 pseudoWaterTable,               & ! pseudo water table
 bigBucket,                      & ! a big bucket (lumped aquifer model)
 noExplicit                        ! no explicit groundwater parameterization
! -------------------------------------------------------------------------------------------------
implicit none
private
public::heatTransf
! local parameters
real(dp),parameter            :: RHsurf=1._dp             ! relative humidity of the surface (-)
! number of soil and snow layers
integer(i4b)                  :: nSoil                    ! number of soil layers
integer(i4b)                  :: nSnow                    ! number of snow layers
integer(i4b)                  :: nLayers                  ! total number of layers
contains


 ! ************************************************************************************************
 ! new subroutine: compute change in temperature over the time step
 ! ************************************************************************************************
 subroutine heatTransf(&
                       ! input
                       dt,                       & ! intent(in): time step (seconds)
                       iter,                     & ! intent(in): current iteration count
                       firstSubstep,             & ! intent(in): flag to indicate if we are processing the first sub-step
                       scalarCanopyTempIter,     & ! intent(in): trial temperature of the vegetation canopy at the current iteration (K)
                       scalarCanopyIceIter,      & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                       scalarCanopyLiqIter,      & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       mLayerTempIter,           & ! intent(in): trial temperature of each model layer at the current iteration (K)
                       mLayerVolFracIceIter,     & ! intent(in): trial volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqIter,     & ! intent(in): trial volumetric fraction of liquid water at the current iteration (-)
                       mLayerMatricHeadIter,     & ! intent(in): trial matric head at the current iteration (m)
                       ! output
                       mLayerTempDiff,           & ! intent(out): iteration increment for temperature (K)
                       scalarCanopyTempNew,      & ! intent(out): new temperature of the vegetation canopy (K)
                       mLayerTempNew,            & ! intent(out): new temperature of each model layer (K)
                       mLayerVolFracIceNew,      & ! intent(out): new volumetric fraction of ice (-)
                       mLayerVolFracLiqNew,      & ! intent(out): new volumetric fraction of liquid water (-)
                       mLayerMatricHeadNew,      & ! intent(out): new matric head (m)
                       err,message)                ! intent(out): error control
 ! model decisions
 USE data_struc,only:model_decisions                            ! model decision structure
 USE var_lookup,only:iLookDECISIONS                             ! named variables for elements of the decision structure
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 ! input
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 integer(i4b),intent(in)       :: iter                     ! iteration count
 logical(i4b),intent(in)       :: firstSubStep             ! flag to indicate if we are processing the first sub-step
 real(dp),intent(in)           :: scalarCanopyTempIter     ! trial temperature of the vegetation canopy at the current iteration (K)
 real(dp),intent(in)           :: scalarCanopyIceIter      ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiqIter      ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature of each snow/soil layer at the current iteration (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! trial volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! trial volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! trial matric head at the current iteration (m)
 ! output
 real(dp),intent(out)          :: mLayerTempDiff(:)        ! iteration increment for temperature (K)
 real(dp),intent(out)          :: scalarCanopyTempNew      ! new temperature of the vegetation canopy (K)
 real(dp),intent(out)          :: mLayerTempNew(:)         ! new temperature of each snow/soil layer (K)
 real(dp),intent(out)          :: mLayerVolFracIceNew(:)   ! new volumetric fraction of ice (-)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! new matric head (m)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! internal
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 ! initialize error control
 err=0; message="heatTransf/"

 ! define number of layers
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)  ! number of snow layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)                   ! total number of layers

 ! *****
 ! wrapper for the temperature change sub-routine...
 ! *************************************************
 ! used as a trick to both
 !  1) save typing; and
 !  2) "protect" variables through the use of the intent attribute
 call heatTransf_muster(&
                        ! input variables from heatTransf routine
                        dt,                                                        & ! intent(in): time step (seconds)
                        iter,                                                      & ! intent(in): current iteration count
                        firstSubstep,                                              & ! intent(in): flag to indicate if we are processing the first sub-step
                        scalarCanopyTempIter,                                      & ! intent(in): trial temperature of the vegetation canopy at the current iteration (K)
                        scalarCanopyIceIter,                                       & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                        scalarCanopyLiqIter,                                       & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                        mLayerTempIter,                                            & ! intent(in): trial temperature of each snow/soil layer at the current iteration (K)
                        mLayerVolFracIceIter,                                      & ! intent(in): trial volumetric fraction of ice at the current iteration (-)
                        mLayerVolFracLiqIter,                                      & ! intent(in): trial volumetric fraction of liquid water at the current iteration (-)
                        mLayerMatricHeadIter,                                      & ! intent(in): trial matric head at the current iteration (m)
                        ! model decisions
                        model_decisions(iLookDECISIONS%num_method)%iDecision,      & ! intent(in): choice of numerical method
                        model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,      & ! intent(in): method used to calculate flux derivatives
                        model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision,      & ! intent(in): type of lower boundary condition for thermodynamics
                        ! index variables
                        indx_data%var(iLookINDEX%nLayers)%dat(1),                  & ! intent(in): number of layers
                        indx_data%var(iLookINDEX%layerType)%dat,                   & ! intent(in): layer type (ix_soil or ix_snow)
                        ! general model parameters
                        mpar_data%var(iLookPARAM%wimplicit),                       & ! intent(in): weight assigned to start-of-step fluxes (-)
                        mpar_data%var(iLookPARAM%snowfrz_scale),                   & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                        mpar_data%var(iLookPARAM%lowerBoundTemp),                  & ! intent(in): temperature of the lower boundary (K)
                        ! vegetation parameters
                        mpar_data%var(iLookPARAM%heightCanopyTop),                 & ! intent(in): height of top of the vegetation canopy above ground surface (m)
                        mpar_data%var(iLookPARAM%heightCanopyBottom),              & ! intent(in): height of bottom of the vegetation canopy above ground surface (m)
                        mpar_data%var(iLookPARAM%specificHeatVeg),                 & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                        mpar_data%var(iLookPARAM%maxMassVegetation),               & ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
                        ! soil parameters
                        mpar_data%var(iLookPARAM%vGn_alpha),                       & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        mpar_data%var(iLookPARAM%vGn_n),                           & ! intent(in): van Genutchen "n" parameter (-)
                        mpar_data%var(iLookPARAM%theta_sat),                       & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),                       & ! intent(in): soil residual volumetric water content (-)
                        mpar_data%var(iLookPARAM%soil_dens_intr),                  & ! intent(in): intrinsic density of soil (kg m-3)
                        mpar_data%var(iLookPARAM%frac_sand),                       & ! intent(in): fraction of sand (-)
                        mpar_data%var(iLookPARAM%frac_silt),                       & ! intent(in): fraction of silt (-)
                        mpar_data%var(iLookPARAM%frac_clay),                       & ! intent(in): fraction of clay (-)
                        mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),               & ! intent(in): van Genutchen "m" parameter (-)
                        ! model state variables
                        ! NOTE: start-of-sub-step values -- protected by the intent(in) attribute
                        mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),          & ! intent(in): temperature of the vegetation canopy (K)
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),           & ! intent(in): mass of ice on the vegetation canopy (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),           & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                        mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1),              & ! intent(in): surface albedo (-)
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat,                   & ! intent(in): temperature of each layer (K)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,             & ! intent(in): volumetric fraction of ice in each layer (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,             & ! intent(in): volumetric fraction of liquid water in each layer (-)
                        ! model cooordinate variables (input)
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat,                  & ! intent(in): depth of each layer (m)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat,                 & ! intent(in): height at the mid-point of each layer (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat,                 & ! intent(in): height at the interface of each layer (m)
                        ! model diagnostic variables (intent inout) -- thermal properties constant over iterations
                        mvar_data%var(iLookMVAR%mLayerTcrit)%dat,                  & ! intent(in): critical soil temperature above which all water is unfrozen (K)
                        mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),   & ! intent(inout): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,           & ! intent(inout): volumetric heat capacity in each layer (J m-3 K-1) 
                        mvar_data%var(iLookMVAR%mLayerThermalC)%dat,               & ! intent(inout): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                        mvar_data%var(iLookMVAR%iLayerThermalC)%dat,               & ! intent(inout): thermal conductivity at the interface of each layer (W m-1 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat,             & ! intent(inout): volumetric fraction of air in each layer (-)
                        ! model diagnostic variables (output)
                        mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat,             & ! intent(out): derivative in the freezing curve (K-1)
                        mvar_data%var(iLookMVAR%iLayerInitNrgFlux)%dat,            & ! intent(out): energy flux at layer interfaces at the start of the time step (W m-2)
                        mvar_data%var(iLookMVAR%iLayerNrgFlux)%dat,                & ! intent(out): energy flux at layer interfaces at the end of the time step (W m-2)
                        ! output variables from heatTransf subroutine
                        mLayerTempDiff,                                            & ! intent(out): iteration increment for temperature (K) 
                        scalarCanopyTempNew,                                       & ! intent(out): new temperature of the vegetation canopy (K)
                        mLayerTempNew,                                             & ! intent(out): new temperature of each snow/soil layer (K)
                        mLayerVolFracIceNew,                                       & ! intent(out): new volumetric fraction of ice (-)
                        mLayerVolFracLiqNew,                                       & ! intent(out): new volumetric fraction of liquid water (-)
                        mLayerMatricHeadNew,                                       & ! intent(out): new matric head (m)
                        err,cmessage)                                                ! intent(out): error control
 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine heatTransf



 
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 ! ************************************************************************************************
 ! private subroutine: wrapper for the temperature change subroutine
 ! ************************************************************************************************
 subroutine heatTransf_muster(&
                              ! input variables from heatTransf routine
                              dt,                           & ! intent(in): time step (seconds)
                              iter,                         & ! intent(in): current iteration count
                              firstSubstep,                 & ! intent(in): flag to indicate if we are processing the first sub-step
                              scalarCanopyTempIter,         & ! intent(in): trial temperature of the vegetation canopy (K)
                              scalarCanopyIceIter,          & ! intent(in): trial mass of ice on the vegetation canopy (kg m-2)
                              scalarCanopyLiqIter,          & ! intent(in): trial mass of liquid water on the vegetation canopy (kg m-2)
                              mLayerTempIter,               & ! intent(in): trial temperature of each snow/soil layer at the current iteration (K)
                              mLayerVolFracIceIter,         & ! intent(in): trial volumetric fraction of ice at the current iteration (-)
                              mLayerVolFracLiqIter,         & ! intent(in): trial volumetric fraction of liquid water at the current iteration (-)
                              mLayerMatricHeadIter,         & ! intent(in): trial matric head at the current iteration (m)
                              ! model decisions
                              num_method,                   & ! intent(in): choice of numerical method
                              fDerivMeth,                   & ! intent(in): method used to calculate flux derivatives
                              bcLowrTdyn,                   & ! intent(in): type of lower boundary condition for thermodynamics
                              ! index variables
                              nLayers,                      & ! intent(in): number of layers
                              layerType,                    & ! intent(in): layer type (ix_soil or ix_snow)
                              ! general model parameters
                              wimplicit,                    & ! intent(in): weight assigned to start-of-step fluxes (-)
                              snowfrz_scale,                & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                              lowerBoundTemp,               & ! intent(in): temperature of the lower boundary (K)
                              ! vegetation parameters
                              heightCanopyTop,              & ! intent(in): height of top of the vegetation canopy above ground surface (m)
                              heightCanopyBottom,           & ! intent(in): height of bottom of the vegetation canopy above ground surface (m)
                              specificHeatVeg,              & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                              maxMassVegetation,            & ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
                              ! soil parameters
                              vGn_alpha,                    & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                        & ! intent(in): van Genutchen "n" parameter (-)
                              theta_sat,                    & ! intent(in): soil porosity (-)
                              theta_res,                    & ! intent(in): soil residual volumetric water content (-)
                              iden_soil,                    & ! intent(in): intrinsic density of soil (kg m-3)
                              frac_sand,                    & ! intent(in): fraction of sand (-)
                              frac_silt,                    & ! intent(in): fraction of silt (-)
                              frac_clay,                    & ! intent(in): fraction of clay (-)
                              vGn_m,                        & ! intent(in): van Genutchen "m" parameter (-)
                              ! model state variables
                              ! NOTE: start-of-sub-step values -- protected by the intent(in) attribute
                              scalarCanopyTemp,             & ! intent(in): temperature of the vegetation canopy (K)
                              scalarCanopyIce,              & ! intent(in): mass of ice on the vegetation canopy (kg m-2)
                              scalarCanopyLiq,              & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                              scalarAlbedo,                 & ! intent(in): surface albedo (-)
                              mLayerTemp,                   & ! intent(in): temperature of each snow/soil layer (K)
                              mLayerVolFracIce,             & ! intent(in): volumetric fraction of ice in each layer (-)
                              mLayerVolFracLiq,             & ! intent(in): volumetric fraction of liquid water in each layer (-)
                              ! model cooordinate variables (input)
                              mLayerDepth,                  & ! intent(in): depth of each layer (m)
                              mLayerHeight,                 & ! intent(in): height at the mid-point of each layer (m)
                              iLayerHeight,                 & ! intent(in): height at the interface of each layer (m)
                              ! model diagnostic variables (intent inout) -- thermal properties constant over iterations
                              mLayerTcrit,                  & ! intent(in): critical soil temperature above which all water is unfrozen (K)
                              scalarBulkVolHeatCapVeg,      & ! intent(inout): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                              mLayerVolHtCapBulk,           & ! intent(inout): volumetric heat capacity in each layer (J m-3 K-1) 
                              mLayerThermalC,               & ! intent(inout): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                              iLayerThermalC,               & ! intent(inout): thermal conductivity at the interface of each layer (W m-1 K-1)
                              mLayerVolFracAir,             & ! intent(inout): volumetric fraction of air in each layer (-)
                              ! model diagnostic variables (output)
                              mLayerdTheta_dTk,             & ! intent(out): derivative in the freezing curve (K-1)
                              iLayerInitNrgFlux,            & ! intent(out): energy flux at layer interfaces at the start of the time step (W m-2)
                              iLayerNrgFlux,                & ! intent(out): energy flux at layer interfaces at the end of the time step (W m-2)
                              ! output variables from heatTransf subroutine
                              mLayerTempDiff,               & ! intent(out): iteration increment for temperature (K)
                              scalarCanopyTempNew,          & ! intent(out): new vegetation temperature (K)
                              mLayerTempNew,                & ! intent(out): new temperature (K)
                              mLayerVolFracIceNew,          & ! intent(out): new volumetric fraction of ice (-)
                              mLayerVolFracLiqNew,          & ! intent(out): new volumetric fraction of liquid water (-)
                              mLayerMatricHeadNew,          & ! intent(out): new matric head (m)
                              err,message)                    ! intent(out): error control
 ! utility modules
 USE diagn_evar_module,only:diagn_evar                        ! compute diagnostic energy variables -- thermal conductivity and heat capacity
 USE vegNrgFlux_module,only:vegNrgFlux                        ! compute energy fluxes for vegetation and ground surface
 USE phseChange_module,only:phseChange                        ! compute change in phase over the time step
 USE snow_utils_module,only:dFracLiq_dTk                      ! differentiate the freezing curve w.r.t. temperature (snow)
 USE soil_utils_module,only:dTheta_dTk                        ! differentiate the freezing curve w.r.t. temperature (soil)
 USE conv_funcs_module,only:relhm2sphm                        ! compute specific humidity 
 USE tridagSolv_module,only:tridag                            ! solve tridiagonal system of equations
 implicit none
 ! input variables from the heatTransf subroutine
 real(dp),intent(in)            :: dt                         ! time step (seconds)
 integer(i4b),intent(in)        :: iter                       ! iteration count
 logical(i4b),intent(in)        :: firstSubStep               ! flag to indicate if we are processing the first sub-step
 real(dp),intent(in)            :: scalarCanopyTempIter       ! trial vegetation temperature (K)
 real(dp),intent(in)            :: scalarCanopyIceIter        ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiqIter        ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)            :: mLayerTempIter(:)          ! trial temperature at the current iteration (K)
 real(dp),intent(in)            :: mLayerVolFracIceIter(:)    ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)            :: mLayerVolFracLiqIter(:)    ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)            :: mLayerMatricHeadIter(:)    ! matric head at the current iteration (m)
 ! model decisions
 integer(i4b),intent(in)        :: num_method                 ! choice of numerical method
 integer(i4b),intent(in)        :: fDerivMeth                 ! method used to calculate flux derivatives
 integer(i4b),intent(in)        :: bcLowrTdyn                 ! type of upper boundary condition for thermodynamics
 ! model index variables
 integer(i4b),intent(in)        :: nLayers                    ! number of layers
 integer(i4b),intent(in)        :: layerType(:)               ! type of the layer (ix_soil or ix_snow)
 ! general model parameters
 real(dp),intent(in)            :: wimplicit                  ! weight assigned to start-of-step fluxes (-)
 real(dp),intent(in)            :: snowfrz_scale              ! scaling parameter for the snow freezing curve (K-1)
 real(dp),intent(in)            :: lowerBoundTemp             ! temperature of the lower boundary (K)
 ! vegetation parameters
 real(dp),intent(in)            :: heightCanopyTop            ! intent(in): height of top of the vegetation canopy above ground surface (m)
 real(dp),intent(in)            :: heightCanopyBottom         ! intent(in): height of bottom of the vegetation canopy above ground surface (m)
 real(dp),intent(in)            :: specificHeatVeg            ! intent(in): specific heat of vegetation (J kg-1 K-1)
 real(dp),intent(in)            :: maxMassVegetation          ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
 ! soil parameters
 real(dp),intent(in)            :: vGn_alpha                  ! van Genutchen "alpha" parameter
 real(dp),intent(in)            :: vGn_n                      ! van Genutchen "n" parameter
 real(dp),intent(in)            :: theta_sat                  ! soil porosity (-)
 real(dp),intent(in)            :: theta_res                  ! soil residual volumetric water content (-)
 real(dp),intent(in)            :: iden_soil                  ! intrinsic density of soil (kg m-3)
 real(dp),intent(in)            :: frac_sand                  ! fraction of sand (-)
 real(dp),intent(in)            :: frac_silt                  ! fraction of silt (-)
 real(dp),intent(in)            :: frac_clay                  ! fraction of clay (-)
 real(dp),intent(in)            :: vGn_m                      ! van Genutchen "m" parameter (-)
 ! model state variables
 ! NOTE: protected with the intent(in) attribute
 real(dp),intent(in)            :: scalarCanopyTemp           ! temperature of the vegetation canopy (K)
 real(dp),intent(in)            :: scalarCanopyIce            ! mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiq            ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: scalarAlbedo               ! surface albedo (-)
 real(dp),intent(in)            :: mLayerTemp(:)              ! temperature of each layer (K)
 real(dp),intent(in)            :: mLayerVolFracIce(:)        ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)        ! volumetric fraction of liquid water in each layer (-)
 ! model coordinate variables (intent in)
 real(dp),intent(in)            :: mLayerDepth(:)             ! depth of each layer (m)
 real(dp),intent(in)            :: mLayerHeight(:)            ! height at the mid-point of each layer (m)
 real(dp),intent(in)            :: iLayerHeight(0:)           ! height at the interface of each layer (m)
 ! model diagnostic variables (intent inout) -- thermal properties constant over iterations
 real(dp),intent(in)            :: mLayerTcrit(:)             ! critical soil temperature above which all water is unfrozen (K)
 real(dp),intent(inout)         :: scalarBulkVolHeatCapVeg    ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
 real(dp),intent(inout)         :: mLayerVolHtCapBulk(:)      ! volumetric heat capacity in each layer (J m-3 K-1) 
 real(dp),intent(inout)         :: mLayerThermalC(:)          ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),intent(inout)         :: iLayerThermalC(0:)         ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),intent(inout)         :: mLayerVolFracAir(:)        ! volumetric fraction of air in each layer (-)
 ! model diagnostic variables (intent out)
 real(dp),intent(out)           :: mLayerdTheta_dTk(:)        ! derivative in the freezing curve (K-1)
 real(dp),intent(out)           :: iLayerInitNrgFlux(0:)      ! energy flux at layer interfaces at the start of the time step (W m-2)
 real(dp),intent(out)           :: iLayerNrgFlux(0:)          ! energy flux at layer interfaces at the end of the time step (W m-2)
 ! output variables from the heatTransf subroutine
 real(dp),intent(out)           :: mLayerTempDiff(:)          ! iteration increment for temperature (K) 
 real(dp),intent(out)           :: scalarCanopyTempNew        ! new temperature of the vegetation canopy (K)
 real(dp),intent(out)           :: mLayerTempNew(:)           ! new temperature of each model layer (K)
 real(dp),intent(out)           :: mLayerVolFracIceNew(:)     ! new volumetric fraction of ice (-)
 real(dp),intent(out)           :: mLayerVolFracLiqNew(:)     ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)           :: mLayerMatricHeadNew(:)     ! new matric head (m)
 integer(i4b),intent(out)       :: err                        ! error code
 character(*),intent(out)       :: message                    ! error message
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! define local variables
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! define general local variables
 character(LEN=256)             :: cmessage                   ! error message of downwind routine
 integer(i4b)                   :: iLayer                     ! index of model layers
 logical(lgt)                   :: printflag                  ! .true. if print progress to the screen
 logical(lgt)                   :: fTranspire                 ! .true. if computing transpiration
 real(dp)                       :: exposedVAI                 ! exposed vegetation area index (m2 m-2)
 real(dp)                       :: canopyDepth                ! depth of the vegetation canopy (m)
 real(dp)                       :: theta                      ! total volumetric water content (liquid plus ice)
 real(dp)                       :: critDiff                   ! temperature difference from critical temperature (K)
 real(dp)                       :: maxdiffTemp(1)             ! maximum difference between temperature input and start-of-step temperature (K)
 real(dp),parameter             :: epsT=1.d-10                ! offset from Tcrit when re-setting iterations at the critical temperature (K)
 integer(i4b)                   :: nUnsat                     ! number of unsaturated layers
 ! define local variables for the fluxes at vegetation and ground surfaces
 real(dp)                       :: canopyNetFlux              ! net energy flux for the vegetation canopy (W m-2)
 real(dp)                       :: groundNetFlux              ! net energy flux for the ground surface (W m-2) 
 real(dp)                       :: dCanopyNetFlux_dCanopyTemp ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dGroundTemp ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dCanopyTemp ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dGroundTemp ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 ! define fluxes at the start of the sub-step
 real(dp)                       :: canopyNetFluxInit          ! net energy flux at the canopy at the start of the substep (W m-2)
 real(dp),dimension(0:nLayers)  :: iLayerNrgFluxInit          ! flux at layer interfaces of the snow-soil system at the start of the substep (W m-2)
 ! define the local variables for the solution
 real(dp)                       :: dTheta_dTkCanopy           ! derivative in fraction liquid water w.r.t. canopy temperature (K-1)
 real(dp),dimension(0:nLayers)  :: dFlux_dTempAbove           ! derivative in flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),dimension(0:nLayers)  :: dFlux_dTempBelow           ! derivative in flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 real(dp)                       :: nrg0,nrg1                  ! energy content at the start of the time step / current iteration (J m-3)
 real(dp)                       :: flx0,flx1                  ! fluxes at the start of the time step / current iteration (J m-3 s-1)
 real(dp)                       :: phse                       ! phase change term (J m-3)
 real(dp)                       :: wtim                       ! weighted time (s-1)
 real(dp)                       :: canopyResidual             ! residual for the canopy (J m-3)
 real(dp),dimension(nLayers)    :: mLayerResidual             ! residual for the snow/soil layers (J m-3)
 ! define tri-diagonal matrix elements for the canopy
 real(dp)                       :: diagCanopy                 ! diagonal element for the canopy
 real(dp)                       :: d_m1Canopy                 ! sub-diagonal element for the derivative in net ground flux w.r.t. canopy temperature (J m-3 K-1)
 real(dp)                       :: d_p1Canopy                 ! super-diagonal element for the derivative in net canopy flux w.r.t. ground temperature (J m-3 K-1)
 ! define tri-diagonal matrix elements for the snow-soil system
 real(dp),dimension(nLayers)    :: diagVector                 ! diagonal (J m-3 K-1)
 real(dp),dimension(nLayers-1)  :: d_m1Vector                 ! sub-diagonal (J m-3 K-1)
 real(dp),dimension(nLayers-1)  :: d_p1Vector                 ! super-diagonal (J m-3 K-1)
 ! define the tri-diagonal matrix
 real(dp),allocatable           :: d_m1(:)                    ! sub-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),allocatable           :: diag(:)                    ! diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),allocatable           :: d_p1(:)                    ! super-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),allocatable           :: rvec(:)                       ! residual vector (J m-3)
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="heatTransf_muster/"

 ! initialize print flag
 printflag=.false.

 ! iterate on the surface temperatire
 if(num_method==iterSurfEnergyBal)then
  err=20; message=trim(message)//'option "iterSurfEnergyBal" not implemented yet';return
 endif

 ! define the canopy depth
 canopyDepth = heightCanopyTop - heightCanopyBottom
 if(heightCanopyBottom > heightCanopyTop)then
  err=20; message=trim(message)//'height of the bottom of the canopy > top of the canopy'; return
 endif


 ! ***** compute diagnostic energy variables (thermal conductivity and volumetric heat capacity)
 if(iter==1)then   ! (constant over the iterations)
  call diagn_evar(&
                  ! input: state variables
                  ! NOTE: using start-of-substep variables
                  scalarCanopyIce,         & ! intent(in): mass of ice on the vegetation canopy (kg m-2)
                  scalarCanopyLiq,         & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                  mLayerVolFracIce,        & ! intent(in): volumetric fraction of ice in each layer (-)
                  mLayerVolFracLiq,        & ! intent(in): volumetric fraction of liquid water in each layer (-)
                  ! input: coordinate variables
                  layerType,               & ! intent(in): type of the layer (snow or soil)
                  mLayerHeight,            & ! intent(in): height of the layer mid-point (top of soil = 0)
                  iLayerHeight,            & ! intent(in): height of the layer interface (top of soil = 0)
                  ! input: model parameters
                  canopyDepth,             & ! intent(in): depth of the vegetation canopy (m)
                  specificHeatVeg,         & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                  maxMassVegetation,       & ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
                  iden_soil,               & ! intent(in): intrinsic density of soil (kg m-3)
                  theta_sat,               & ! intent(in): soil porosity (-)
                  frac_sand,               & ! intent(in): fraction of sand (-)
                  frac_silt,               & ! intent(in): fraction of silt (-)
                  frac_clay,               & ! intent(in): fraction of clay (-)
                  ! output: diagnostic variables
                  scalarBulkVolHeatCapVeg, & ! intent(out): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                  mLayerVolHtCapBulk,      & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1)
                  mLayerThermalC,          & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                  iLayerThermalC,          & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
                  mLayerVolFracAir,        & ! intent(out): volumetric fraction of air in each layer (-)
                  ! output: error control
                  err,cmessage)              ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif


 ! ***** compute energy fluxes at vegetation and ground surfaces
 call vegNrgFlux(&
                 ! input
                 dt,                            & ! intent(in): time step (seconds)
                 iter,                          & ! intent(in): iteration index
                 firstSubStep,                  & ! intent(in): flag to indicate if we are processing the first sub-step
                 scalarCanopyTempIter,          & ! intent(in): trial value of canopy temperature (K)
                 mLayerTempIter(1),             & ! intent(in): trial value of ground temperature (K)
                 scalarCanopyIceIter,           & ! intent(in): trial mass of ice on the vegetation canopy (kg m-2)
                 scalarCanopyLiqIter,           & ! intent(in): trial mass of liquid water on the vegetation canopy (kg m-2)
                 ! output
                 exposedVAI,                    & ! intent(out): exposed vegetation area index (m2 m-2)
                 canopyNetFlux,                 & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                 groundNetFlux,                 & ! intent(out): net energy flux for the ground surface (W m-2)
                 dCanopyNetFlux_dCanopyTemp,    & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                 dCanopyNetFlux_dGroundTemp,    & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                 dGroundNetFlux_dCanopyTemp,    & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                 dGroundNetFlux_dGroundTemp,    & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                 err,cmessage)                    ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! ***** compute fluxes at layer interfaces and their derivatives (J m-2 s-1)
 call iLayer_nrg(&
                 ! (model control variables)
                 fDerivMeth,                    & ! intent(in): method used to compute derivatives (numerical or analytical)
                 ! (input)
                 mLayerDepth,                   & ! intent(in): depth of each layer (m)
                 mLayerHeight,                  & ! intent(in): height of layer mid-points (m)
                 iLayerThermalC,                & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                 mLayerTempIter,                & ! intent(in): trial temperature at the current iteration (K)
                 groundNetFlux,                 & ! intent(in): total flux at the ground surface (W m-2)
                 dGroundNetFlux_dGroundTemp,    & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                 lowerBoundTemp,                & ! intent(in): temperature of the lower boundary (K)
                 ! (output)
                 iLayerNrgFlux,                 & ! intent(out): energy flux at the layer interfaces (W m-2)
                 dFlux_dTempAbove,              & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                 dFlux_dTempBelow,              & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                 err,cmessage)                    ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! ***** allocate space for the tri-diagonal matrix
 if(exposedVAI > 0._dp)then
  allocate(d_m1(nLayers),diag(nLayers+1),d_p1(nLayers),rvec(nLayers+1), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'problem allocating space for the tri-diag matrix'; return; endif
 else
  err=20; message=trim(message)//'non-vegetated surfaces not implemented yet'; return
 endif 


 ! ***** assign initial fluxes
 if(iter==1)then
  ! check that the temperature matches the temperature at the start of the step
  maxdiffTemp = maxval(abs(mLayerTempIter - mLayerTemp))
  if(maxdiffTemp(1) > 1.d-8)then; err=20; message=trim(message)//'first guess for temperature must match start-of-step value'; return; endif
  ! assign initial fluxes
  canopyNetFluxInit = canopyNetFlux
  iLayerNrgFluxInit = iLayerNrgFlux
 endif


 ! ***** compute the residual for the vegetation canopy
 ! (compute individual terms)
 nrg0 = scalarBulkVolHeatCapVeg*scalarCanopyTemp                      ! energy content at the start of the time step (J m-3)
 nrg1 = scalarBulkVolHeatCapVeg*scalarCanopyTempIter                  ! energy content at the current iteration (J m-3)
 flx0 = canopyNetFluxInit/canopyDepth                                 ! flux at the start of the time step (J m-3 s-1)
 flx1 = canopyNetFlux/canopyDepth                                     ! flux at the current iteration (J m-3 s-1)
 phse = LH_fus*(scalarCanopyIceIter - scalarCanopyIce)/canopyDepth    ! phase change term (J m-3)
 ! (compute residuals)
 canopyResidual = nrg1 - (nrg0 + (flx0*wimplicit + flx1*(1._dp - wimplicit))*dt + phse)


 ! ***** compute the residual vector for all snow/soil layers (J m-3)
 do iLayer=1,nLayers
  ! (compute individual terms)
  nrg0 = mLayerVolHtCapBulk(iLayer)*mLayerTemp(iLayer)                                        ! energy content at the start of the time step (J m-3)
  nrg1 = mLayerVolHtCapBulk(iLayer)*mLayerTempIter(iLayer)                                    ! energy content at the current iteration (J m-3)
  flx0 = -(iLayerInitNrgFlux(iLayer) - iLayerInitNrgFlux(iLayer-1))/mLayerDepth(iLayer)       ! flux at the start of the time step (J m-3 s-1)
  flx1 = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)               ! flux at the current iteration (J m-3 s-1)
  phse = LH_fus*iden_ice*(mLayerVolFracIceIter(iLayer) - mLayerVolFracIce(iLayer))            ! phase change term (J m-3)
  ! (compute residuals)
  mLayerResidual(iLayer) = nrg1 - (nrg0 + (flx0*wimplicit + flx1*(1._dp - wimplicit))*dt + phse)
 end do


 ! ***** compute the derivative in the freezing curve w.r.t. temperature (K-1)
 ! * vegetation canopy
 if(scalarCanopyIceIter > 0._dp)then
  dTheta_dTkCanopy = dFracLiq_dTk(scalarCanopyTempIter,snowfrz_scale)
 else
  dTheta_dTkCanopy = 0._dp
 endif
 ! * all snow-soil layers
 do iLayer=1,nLayers
  select case(layerType(iLayer))
   case(ix_snow) ! (snow layers)
    theta = mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iLayer)
    mLayerdTheta_dTk(iLayer) = dFracLiq_dTk(mLayerTempIter(iLayer),snowfrz_scale)*theta
   case(ix_soil) ! (soil layers)
    if(mLayerVolFracIceIter(iLayer)>0._dp)then
     mLayerdTheta_dTk(iLayer) = dTheta_dTk(mLayerTempIter(iLayer),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
    else
     mLayerdTheta_dTk(iLayer) = 0._dp
    endif
   case default; err=40; message=trim(message)//"cannot identify the layer as snow or soil"; return
  endselect
 end do

 print*, 'wimplicit = ', wimplicit


 ! compute the weighted time for end-of-step values
 wtim = (1._dp - wimplicit)*dt  ! weighted time

 ! ***** assemble the tri-diagonal matrix for the canopy
 diagCanopy = (wtim/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp + dGroundNetFlux_dCanopyTemp) + dTheta_dTkCanopy*LH_fus*iden_water + scalarBulkVolHeatCapVeg
 d_m1Canopy = (wtim/canopyDepth)*(-dGroundNetFlux_dCanopyTemp)
 d_p1Canopy = (wtim/canopyDepth)*( dCanopyNetFlux_dGroundTemp)

 ! ***** assemble the tri-diagonal matrix for the snow-soil layers
 diagVector = (wtim/mLayerDepth)*(-dFlux_dTempBelow(0:nLayers-1) + dFlux_dTempAbove(1:nLayers)) + mLayerdTheta_dTk*LH_fus*iden_water + mLayerVolHtCapBulk
 d_m1Vector = (wtim/mLayerDepth(1:nLayers-1))*(-dFlux_dTempAbove(1:nLayers-1) )
 d_p1Vector = (wtim/mLayerDepth(1:nLayers-1))*( dFlux_dTempBelow(1:nLayers-1) )

 ! ***** combine vectors
 d_m1 = (/d_m1Canopy,d_m1Vector/)
 diag = (/diagCanopy,diagVector/)
 d_p1 = (/d_p1Canopy,d_p1Vector/)
 rvec = (/canopyResidual,mLayerResidual/)

 write(*,'(a,100(f20.10,1x))') 'diag = ', diag
 write(*,'(a,100(f20.10,1x))') 'd_m1 = ', d_m1
 write(*,'(a,100(f20.10,1x))') 'd_p1 = ', d_p1

 pause ' assembled tri-diag system'

 ! ***** solve the tridiagonal system of equations -- returns mLayerTempDiff
 call tridag(d_m1,                    & ! intent(in): sub-diagonal elements of the tridiagonal system (J m-3 K-1)
             diag,                    & ! intent(in): diagonal elements of the tridiagonal system (J m-3 K-1)
             d_p1,                    & ! intent(in): super-diagonal elements of the tridiagonal system (J m-3 K-1)
             -rvec,                   & ! intent(in): residual vector (J m-3)
             mLayerTempDiff,          & ! intent(out): temperature increment (K)
             err,cmessage)              ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! adjust del temperature in cases where snow temperature exceeds Tfreeze -- use bi-section
 if(nSnow>0)then
  do iLayer=1,nSnow
   if(mLayerTempIter(iLayer)+mLayerTempDiff(iLayer) > Tfreeze)then
    mLayerTempDiff(iLayer) = (Tfreeze-mLayerTempIter(iLayer))*0.5_dp
   endif
  end do
 endif

 ! adjust del temperature in cases where soil temperature crosses the critical temperature
 do iLayer=nSnow+1,nLayers
  ! get the difference from the critical temperature (K)
  critDiff = mLayerTcrit(iLayer-nSnow) - mLayerTempIter(iLayer)
  ! set temperature to Tcrit in cases where temperatures cross Tcrit
  if(critDiff > 0._dp)then  ! (mLayerTempIter < Tcrit)
   if(mLayerTempDiff(iLayer) > critDiff) mLayerTempDiff(iLayer) = critDiff + epsT
  else                      ! (mLayerTempIter > Tcrit)
   if(mLayerTempDiff(iLayer) < critDiff) mLayerTempDiff(iLayer) = critDiff - epsT
  endif
 end do

 ! update temperature
 mLayerTempNew = mLayerTempIter + mLayerTempDiff

 ! compute phase change
 call phsechange(mLayerTempNew,       & ! intent(in): new temperature vector (K)
                 mLayerMatricHeadIter,& ! intent(in): matric head at the current iteration (m)
                 mLayerVolFracLiqIter,& ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                 mLayerVolFracIceIter,& ! intent(in): volumetric fraction of ice at the current iteration (-)
                 mLayerMatricHeadNew, & ! intent(out): new matric head (m)
                 mLayerVolFracLiqNew, & ! intent(out): new volumetric fraction of liquid water (-)
                 mLayerVolFracIceNew, & ! intent(out): new volumetric fraction of ice (-)
                 err,cmessage)          ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ***** update the fluxes at the layer interfaces
 do iLayer=0,nLayers
  if(iLayer==0)then;           iLayerNrgFlux(iLayer) = iLayerNrgFlux(iLayer) + dFlux_dTempBelow(iLayer)*mLayerTempDiff(iLayer+1)
  elseif(iLayer==nLayers)then; iLayerNrgFlux(iLayer) = iLayerNrgFlux(iLayer) + dFlux_dTempAbove(iLayer)*mLayerTempDiff(iLayer)
  else;                        iLayerNrgFlux(iLayer) = iLayerNrgFlux(iLayer) + dFlux_dTempAbove(iLayer)*mLayerTempDiff(iLayer) &
                                                                             + dFlux_dTempBelow(iLayer)*mLayerTempDiff(iLayer+1)
  endif
 end do ! (looping through layers)

 ! deallocate space for the tri-diagonal matrix
 deallocate(d_m1,diag,d_p1,rvec, stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem deallocating space for the tri-diag matrix'; return; endif

 ! ====================================================================================================================

 end subroutine heatTransf_muster


 ! ************************************************************************************************
 ! private subroutine: compute energy fluxes at layer interfaces, and their derivatives
 ! ************************************************************************************************
 subroutine iLayer_nrg(&
                       ! (model control variables)
                       fDerivMeth,                         & ! intent(in): method used to compute derivatives (numerical or analytical)
                       ! (input)
                       mLayerDepth,                        & ! intent(in): depth of each layer (m)
                       mLayerHeight,                       & ! intent(in): height of layer mid-points (m)
                       iLayerThermalC,                     & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                       mLayerTempTrial,                    & ! intent(in): trial temperature at the current iteration (K)
                       groundNetFlux,                      & ! intent(in): total flux at the ground surface (W m-2)
                       dGroundNetFlux_dGroundTemp,         & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                       lowerBoundTemp,                     & ! intent(in): temperature of the lower boundary (K)
                       ! (output)
                       iLayerNrgFlux,                      & ! intent(out): energy flux at the layer interfaces (W m-2)
                       dFlux_dTempAbove,                   & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                       dFlux_dTempBelow,                   & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                       err,message)                          ! intent(out): error control
 ! compute derivative in fluxes at layer interfaces w.r.t. temperature in the layer above and the layer below
 implicit none
 ! input
 integer(i4b),intent(in)       :: fDerivMeth                 ! intent(in): method used to calculate derivatives
 real(dp),intent(in)           :: mLayerDepth(:)             ! intent(in): depth of each layer (m)
 real(dp),intent(in)           :: mLayerHeight(:)            ! intent(in): height of layer mid-points (m)
 real(dp),intent(in)           :: iLayerThermalC(0:)         ! intent(in): thermal conductivity at layer interfaces (W m-1)
 real(dp),intent(in)           :: mLayerTempTrial(:)         ! intent(in): trial temperature at the current iteration (K)
 real(dp),intent(in)           :: groundNetFlux              ! intent(in): total flux at the ground surface (W m-2)
 real(dp),intent(in)           :: dGroundNetFlux_dGroundTemp ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(in)           :: lowerBoundTemp             ! intent(in): temperature of the lower boundary (K)
 ! output
 real(dp),intent(out)          :: iLayerNrgFlux(0:)          ! intent(out): energy flux at the layer interfaces (W m-2)
 real(dp),intent(out)          :: dFlux_dTempAbove(0:)       ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),intent(out)          :: dFlux_dTempBelow(0:)       ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 integer(i4b),intent(out)      :: err                        ! intent(out): error code
 character(*),intent(out)      :: message                    ! intent(out): error message
 ! local variables
 integer(i4b)                  :: iLayer                     ! index of model layers
 real(dp),parameter            :: dx=1.e-8_dp                ! finite difference increment (K)
 real(dp)                      :: dz                         ! height difference (m)
 real(dp)                      :: flux0,flux1,flux2          ! fluxes used to calculate derivatives (W m-2)
 ! initialize error control
 err=0; message='iLayer_nrg/'

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 ! compute fluxes within the domain -- positive downwards
 do iLayer=0,nLayers
  ! compute flux at the upper boundary -- positive downwards
  if(iLayer==0)then;           iLayerNrgFlux(0)       = groundNetFlux
  ! compute fluxes at the lower boundary -- positive downwards
  elseif(iLayer==nLayers)then; iLayerNrgFlux(nLayers) = -iLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempTrial(iLayer))/(mLayerDepth(iLayer)*0.5_dp)
  ! compute fluxes within the domain -- positive downwards
  else;                        iLayerNrgFlux(iLayer)  = -iLayerThermalC(iLayer)*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer)) / &
                                                                                (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
  endif ! (the type of layer)
 end do

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the derivative in fluxes at layer interfaces w.r.t temperature in the layer above and the layer below *****
 ! -------------------------------------------------------------------------------------------------------------------------

 ! initialize un-used elements
 dFlux_dTempBelow(nLayers) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

 ! loop through INTERFACES...
 do iLayer=0,nLayers

  ! ***** the upper boundary -- ** NOTE: dTotalSurfaceFlux_dTemp was computed previously using fDerivMeth
  if(iLayer==0)then  ! (upper boundary)
   dFlux_dTempBelow(iLayer) = dGroundNetFlux_dGroundTemp

  ! ***** the lower boundary
  elseif(iLayer==nLayers)then  ! (lower boundary)
   dz = mLayerDepth(iLayer)*0.5_dp
   if(fDerivMeth==analytical)then    ! ** analytical derivatives
    dFlux_dTempAbove(iLayer) = iLayerThermalC(iLayer)/dz
   else                              ! ** numerical derivatives
    flux0 = -iLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempTrial(iLayer)   ))/dz
    flux1 = -iLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempTrial(iLayer)+dx))/dz
    dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
   endif

  ! ***** internal layers
  else
   dz = (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
   if(fDerivMeth==analytical)then    ! ** analytical derivatives
    dFlux_dTempAbove(iLayer) =  iLayerThermalC(iLayer)/dz
    dFlux_dTempBelow(iLayer) = -iLayerThermalC(iLayer)/dz
   else                              ! ** numerical derivatives
    flux0 = -iLayerThermalC(iLayer)*( mLayerTempTrial(iLayer+1)     -  mLayerTempTrial(iLayer)    ) / dz
    flux1 = -iLayerThermalC(iLayer)*( mLayerTempTrial(iLayer+1)     - (mLayerTempTrial(iLayer)+dx)) / dz
    flux2 = -iLayerThermalC(iLayer)*((mLayerTempTrial(iLayer+1)+dx) -  mLayerTempTrial(iLayer)    ) / dz
    dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
    dFlux_dTempBelow(iLayer) = (flux2 - flux0)/dx
   endif

  endif  ! type of layer (upper, internal, or lower)

 end do  ! (looping through layers)
 
 end subroutine iLayer_nrg


 ! ************************************************************************************************
 ! private subroutine: compute the Jacobian matrix
 ! ************************************************************************************************
 subroutine cmpJacobian(&
                        ! input: model control
                        dt,                            & ! intent(in): time step (seconds)
                        iter,                          & ! intent(in): iteration index
                        firstSubStep,                  & ! intent(in): flag to indicate if we are processing the first sub-step
                        ! input: canopy variables
                        scalarCanopyTempIter,          & ! intent(in): trial value of canopy temperature (K)
                        scalarCanopyIceIter,           & ! intent(in): trial mass of ice on the vegetation canopy (kg m-2)
                        scalarCanopyLiqIter,           & ! intent(in): trial mass of liquid water on the vegetation canopy (kg m-2)
                        ! input: variables for the snow-soil system
                        mLayerDepth,                   & ! intent(in): depth of each layer (m)
                        mLayerHeight,                  & ! intent(in): height of layer mid-points (m)
                        iLayerThermalC,                & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                        mLayerTempIter,                & ! intent(in): trial temperature at the current iteration (K)
                        lowerBoundTemp,                & ! intent(in): temperature of the lower boundary (K)
                        ! output: error control
                        err,message)                     ! intent(out): error control
 implicit none
 ! input: model control
 real(dp),intent(in)           :: dt                     ! time step (seconds)
 integer(i4b),intent(in)       :: iter                   ! iteration index
 logical(i4b),intent(in)       :: firstSubStep           ! flag to indicate if we are processing the first sub-step
 ! input: canopy variables
 real(dp),intent(in)           :: canopyTempTrial        ! trial value of canopy temperature (K)
 real(dp),intent(in)           :: canopyIceTrial         ! trial value of mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)           :: canopyLiqTrial         ! trial value of mass of liquid water on the vegetation canopy (kg m-2)
 ! input: variables for the snow-soil system
 real(dp),intent(in)           :: mLayerDepth(:)         ! depth of each layer (m)
 real(dp),intent(in)           :: mLayerHeight(:)        ! height of layer mid-points (m)
 real(dp),intent(in)           :: iLayerThermalC(0:)     ! thermal conductivity at layer interfaces (W m-1)
 real(dp),intent(in)           :: mLayerTempTrial(:)     ! trial temperature at the current iteration (K)
 real(dp),intent(in)           :: lowerBoundTemp         ! temperature of the lower boundary (K)
 ! output: error control
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 ! ----------------------------------------------------------------------------------------------------------
 ! local variables


 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='cmpJacobian/'







 ! ***** compute energy fluxes at vegetation and ground surfaces
 call vegNrgFlux(&
                 ! input
                 dt,                            & ! intent(in): time step (seconds)
                 iter,                          & ! intent(in): iteration index
                 firstSubStep,                  & ! intent(in): flag to indicate if we are processing the first sub-step
                 scalarCanopyTempIter,          & ! intent(in): trial value of canopy temperature (K)
                 mLayerTempIter(1),             & ! intent(in): trial value of ground temperature (K)
                 scalarCanopyIceIter,           & ! intent(in): trial mass of ice on the vegetation canopy (kg m-2)
                 scalarCanopyLiqIter,           & ! intent(in): trial mass of liquid water on the vegetation canopy (kg m-2)
                 ! output
                 exposedVAI,                    & ! intent(out): exposed vegetation area index (m2 m-2)
                 canopyNetFlux,                 & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                 groundNetFlux,                 & ! intent(out): net energy flux for the ground surface (W m-2)
                 dCanopyNetFlux_dCanopyTemp,    & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                 dCanopyNetFlux_dGroundTemp,    & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                 dGroundNetFlux_dCanopyTemp,    & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                 dGroundNetFlux_dGroundTemp,    & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                 err,cmessage)                    ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! ***** compute fluxes at layer interfaces and their derivatives (J m-2 s-1)
 call iLayer_nrg(&
                 ! (model control variables)
                 fDerivMeth,                    & ! intent(in): method used to compute derivatives (numerical or analytical)
                 ! (input)
                 mLayerDepth,                   & ! intent(in): depth of each layer (m)
                 mLayerHeight,                  & ! intent(in): height of layer mid-points (m)
                 iLayerThermalC,                & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                 mLayerTempIter,                & ! intent(in): trial temperature at the current iteration (K)
                 groundNetFlux,                 & ! intent(in): total flux at the ground surface (W m-2)
                 dGroundNetFlux_dGroundTemp,    & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                 lowerBoundTemp,                & ! intent(in): temperature of the lower boundary (K)
                 ! (output)
                 iLayerNrgFlux,                 & ! intent(out): energy flux at the layer interfaces (W m-2)
                 dFlux_dTempAbove,              & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                 dFlux_dTempBelow,              & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                 err,cmessage)                    ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif





                        ! input: model control
                        ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                        bc_lower,                     & ! intent(in): index defining the type of boundary conditions
                        ! input: state variables
                        mLayerMatricHeadInput,        & ! intent(in): matric head in each layer (m)
                        mLayerVolFracLiqInput,        & ! intent(in): volumetric liquid water content (-)
                        mLayerVolFracIceInput,        & ! intent(in): volumetric ice content in each layer (-)
                        scalarAquiferStorageInput,    & ! intent(in): aquifer storage (m)
                        ! output: error control
                        err,message)                    ! intent(out): error control
 implicit none
 ! input: model control
 integer(i4b),intent(in)          :: ixRichards                  ! index defining the option for Richards' equation (moisture or mixdform)
 integer(i4b),intent(in)          :: bc_lower                    ! index defining the type of boundary conditions
 ! input: trial model state variables
 real(dp),intent(in)              :: mLayerMatricHeadInput(:)    ! input value of matric head in each layer  (m)
 real(dp),intent(in)              :: mLayerVolFracLiqInput(:)    ! input value of volumetric fraction of liquid water in each layer (-)
 real(dp),intent(in)              :: mLayerVolFracIceInput(:)    ! input value of volumetric fraction of ice in each layer (-)
 real(dp),intent(in)              :: scalarAquiferStorageInput   ! input value of aquifer storage (m)
 ! output: error control
 integer(i4b),intent(out)         :: err                         ! error code
 character(*),intent(out)         :: message                     ! error message
  ! -------------------------------------------------------------------------------------------------------------------------------------------------






end module heatTransf_module
