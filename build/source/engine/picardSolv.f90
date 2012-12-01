module picardSolv_module
USE nrtype
implicit none
private
public::picardSolv
contains

 ! ************************************************************************************************
 ! new subroutine: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine picardSolv(dt,maxiter,niter,err,message)
 ! provide access to subroutines
 USE diagn_evar_module,only:diagn_evar          ! compute diagnostic energy variables -- thermal conductivity and heat capacity
 USE heatTransf_module,only:heatTransf          ! compute change in temperature over the time step
 USE phseChange_module,only:phseChange          ! compute change in phase over the time step
 USE snowHydrol_module,only:snowHydrol          ! compute liquid water flow through the snowpack
 USE soilHydrol_module,only:soilHydrol          ! compute change in mass over the time step for the soil
 USE surfAlbedo_module,only:surfAlbedo          ! compute surface albedo
 USE snwDensify_module,only:snwDensify          ! compute densification of snow
 USE soil_utils_module,only:crit_soilT          ! compute the critical temperature above which all water is unfrozen
 ! provide access to data
 USE multiconst,only:&
                     gravity,      & ! acceleration of gravity              (m s-2)
                     Tfreeze,      & ! temperature at freezing              (K)
                     LH_fus,       & ! latent heat of fusion                (J kg-1)
                     iden_ice,     & ! intrinsic density of ice             (kg m-3)
                     iden_water      ! intrinsic density of liquid water    (kg m-3)
 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only:       &
  equilWaterTable,                 & ! equilibrium water table
  pseudoWaterTable,                & ! pseudo water table
  bigBucket,                       & ! a big bucket (lumped aquifer model)
  noExplicit                         ! no explicit groundwater parameterization
 ! model decision structures
 USE data_struc,only:model_decisions ! model decision structure
 USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
 ! general data structures
 USE data_struc,only:forcFileInfo                                  ! extract time step of forcing data
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX               ! named variables for structure elements
 implicit none
 ! define output
 real(dp),intent(in)                  :: dt                       ! time step (seconds)
 integer(i4b),intent(in)              :: maxiter                  ! maximum number of iterations
 integer(i4b),intent(out)             :: niter                    ! number of iterations
 integer(i4b),intent(out)             :: err                      ! error code
 character(*),intent(out)             :: message                  ! error message
 ! local pointers to snow parameters
 real(dp),pointer                     :: Fcapil                   ! capillary retention as a fraction of the total pore volume (-)
 real(dp),pointer                     :: snowfrz_scale            ! scaling parameter for the snow freezing curve (K-1)
 ! local pointers to soil/veg parameters
 real(dp),pointer                     :: soilAlbedo               ! soil albedo (-)
 real(dp),pointer                     :: vGn_alpha                ! van Genutchen "alpha" parameter
 real(dp),pointer                     :: vGn_n                    ! van Genutchen "n" parameter
 real(dp),pointer                     :: theta_sat                ! soil porosity (-)
 real(dp),pointer                     :: theta_res                ! soil residual volumetric water content (-)
 real(dp),pointer                     :: specificYield            ! specific yield (-)
 real(dp),pointer                     :: specificStorage          ! specific storage coefficient (m-1)
 real(dp),pointer                     :: f_impede                 ! ice impedence factor (-)
 real(dp),pointer                     :: rootingDepth             ! rooting depth (m)
 real(dp),pointer                     :: rootDistExp              ! exponent for the vertical distriution of root density (-)
 ! local pointers to algorithmic control parameters
 real(dp),pointer                     :: wimplicit                ! weight assigned to start-of-step fluxes (-)
 real(dp),pointer                     :: relConvTol_liquid        ! relative convergence tolerance for vol frac liq water (-)
 real(dp),pointer                     :: absConvTol_liquid        ! absolute convergence tolerance for vol frac liq water (-)
 real(dp),pointer                     :: relConvTol_matric        ! relative convergence tolerance for matric head (-)
 real(dp),pointer                     :: absConvTol_matric        ! absolute convergence tolerance for matric head (m)
 real(dp),pointer                     :: relConvTol_energy        ! relative convergence tolerance for energy (-)
 real(dp),pointer                     :: absConvTol_energy        ! absolute convergence tolerance for energy (J m-3)
 real(dp),pointer                     :: relConvTol_aquifr        ! relative convergence tolerance for aquifer storage (-)
 real(dp),pointer                     :: absConvTol_aquifr        ! absolute convergence tolerance for aquifer storage (J m-3)
 ! local pointers to derived model variables that are constant over the simulation period
 real(dp),pointer                     :: vGn_m                    ! van Genutchen "m" parameter (-)
 real(dp),pointer                     :: mLayerRootDensity(:)     ! fraction of roots in each soil layer (-)
 ! local pointers to scalar state variables
 real(dp),pointer                     :: scalarSWE                ! SWE (kg m-2)
 real(dp),pointer                     :: surfaceAlbedo            ! surface albedo (-) 
 real(dp),pointer                     :: scalarSfcMeltPond        ! ponded water caused by melt of the "snow without a layer" (kg m-2)
 real(dp),pointer                     :: scalarAquiferStorage     ! relative aquifer storage above the bottom of the soil profile (m)
 ! local pointers to diagnostic scalar variables
 real(dp),pointer                     :: scalarRainPlusMelt       ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 real(dp),pointer                     :: scalarSnowDepth          ! total snow depth (m)
 real(dp),pointer                     :: scalarWaterTableDepth    ! water table depth (m)
 real(dp),pointer                     :: scalarPotentialET        ! potential ET (kg m-2 s-1)
 real(dp),pointer                     :: scalarMassLiquid         ! evaporation/dew (kg m-2 s-1)
 real(dp),pointer                     :: scalarMassSolid          ! sublimation/frost (kg m-2 s-1)
 real(dp),pointer                     :: scalarSoilInflux         ! influx of water at the top of the soil profile (m s-1)
 real(dp),pointer                     :: scalarSoilBaseflow       ! total baseflow from the soil profile (m s-1)
 real(dp),pointer                     :: scalarSoilDrainage       ! drainage from the bottom of the soil profile (m s-1)
 real(dp),pointer                     :: scalarSoilEjection       ! total water ejected from soil layers (m s-1)
 real(dp),pointer                     :: scalarSoilWatBalError    ! error in the total soil water balance (kg m-2) 
 real(dp),pointer                     :: scalarTotalSoilLiq       ! total mass of liquid water in the soil (kg m-2)
 real(dp),pointer                     :: scalarTotalSoilIce       ! total mass of ice in the soil (kg m-2)
 ! local pointers to model coordinate variables
 real(dp),pointer                     :: iLayerHeight(:)          ! height of layer interfaces (m)
 ! local pointers to model state variables -- all layers
 real(dp),pointer                     :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer                     :: mLayerDepth(:)           ! depth of each layer (m)
 real(dp),pointer                     :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                     :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 ! local pointers to model state variables -- soil layers
 real(dp),pointer                     :: mLayerMatricHead(:)      ! matric head in each ***soil*** layer (m)
 ! local pointers to model diagnostic variables -- all layers
 real(dp),pointer                     :: mLayerVolFracAir(:)      ! volumetric fraction of air in each layer (-)
 real(dp),pointer                     :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer                     :: mLayerdTheta_dTk(:)      ! derivative in the freezing curve (K-1)
 real(dp),pointer                     :: mLayerMeltFreeze(:)      ! melt/freeze in each layer (kg m-3 s-1)
 real(dp),pointer                     :: mLayerInfilFreeze(:)     ! increase in volumetric ice content caused by freezing infiltrating flux (-)
 ! local pointers to model diagnostic variables -- soil layers
 real(dp),pointer                     :: mLayerTcrit(:)           ! critical soil temperature above which all water is unfrozen (K)
 real(dp),pointer                     :: mLayerdTheta_dPsi(:)     ! derivative in the soil water characteristic (m-1)
 real(dp),pointer                     :: iLayerInitLiqFluxSoil(:) ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
 real(dp),pointer                     :: iLayerLiqFluxSoil(:)     ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
 real(dp),pointer                     :: mLayerInitBaseflow(:)    ! baseflow from each soil layer at the start of the time step (m s-1)
 real(dp),pointer                     :: mLayerBaseflow(:)        ! baseflow from each soil layer (m s-1)
 real(dp),pointer                     :: mLayerInitTranspire(:)   ! transpiration loss from each soil layer at the start of the time step (m s-1)
 real(dp),pointer                     :: mLayerTranspire(:)       ! transpiration loss from each soil layer (m s-1)
 real(dp),pointer                     :: mLayerInitEjectWater(:)  ! water ejected from each soil layer at start-of-step (m s-1)
 real(dp),pointer                     :: mLayerEjectWater(:)      ! water ejected from each soil layer (m s-1)
 real(dp),pointer                     :: iLayerSatHydCond(:)      ! saturated hydraulic conductivity at layer interfaces (m s-1)
 ! local pointers to model diagnostic variables -- aquifer only
 real(dp),pointer                     :: scalarAquiferRootFrac    ! fraction of roots below the unsaturated zone (-)
 real(dp),pointer                     :: scalarInitAquiferTranspire ! initial aquifer transpiration rate (m s-1)
 real(dp),pointer                     :: scalarInitAquiferBaseflow  ! initial aquifer baseflow rate (m s-1)
 real(dp),pointer                     :: scalarInitAquiferRecharge  ! initial aquifer recharge rate (m s-1)
 real(dp),pointer                     :: scalarAquiferTranspire     ! aquifer tramspiration rate (m s-1)
 real(dp),pointer                     :: scalarAquiferBaseflow      ! aquifer baseflow rate (m s-1)
 real(dp),pointer                     :: scalarAquiferRecharge      ! aquifer recharge rate (m s-1)
 real(dp),pointer                     :: scalarAquiferBalError    ! error in the aquifer water balance (kg m-2) 
 ! local pointers to model index variables
 integer(i4b),pointer                 :: nLayers                  ! number of layers
 integer(i4b),pointer                 :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local model state variables
 real(dp),allocatable                 :: mLayerTempIter(:)        ! temperature vector at the current iteration (K)
 real(dp),allocatable                 :: mLayerTempNew(:)         ! temperature vector at the next iteration (K)
 real(dp),allocatable                 :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),allocatable                 :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice at the next iteration (-)
 real(dp),allocatable                 :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),allocatable                 :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the next iteration (-)
 real(dp),allocatable                 :: mLayerMatricHeadIter(:)  ! matric head at the current iteration (-)
 real(dp),allocatable                 :: mLayerMatricHeadNew(:)   ! matric head at the next iteration (-)
 real(dp)                             :: scalarAquiferStorageIter ! trial value of aquifer storage (m)
 real(dp)                             :: scalarAquiferStorageNew  ! aquifer storage at the end of the time step (m)
 ! define local model diagnostic variables
 real(dp),allocatable                 :: mLayerInfilFreezeNew(:)  ! increase in volumetric ice content caused by freezing infiltrating flux (-)
 ! define local error monitoring variables
 real(dp),allocatable                 :: mLayerTempIncr(:)        ! change in temperature from one iteration to the next (K)
 real(dp),allocatable                 :: mLayerTempIncrOld(:)     ! change in temperature from one iteration to the next -- previous iteration (K)
 real(dp),allocatable                 :: mLayerMatIncr(:)         ! change in matric head from one iteration to the next (m)
 real(dp),allocatable                 :: mLayerLiqIncr(:)         ! change in volumetric liquid water content from one iteration to the next (-)
 real(dp),allocatable                 :: mLayerNrgIncr(:)         ! change in energy from one iteration to the next (J m-3)
 real(dp)                             :: scalarAqiIncr            ! change in aquifer storage from one iteration to the next (m)
 real(dp),allocatable                 :: mLayerNrgError(:)        ! energy error in each layer (J m-3)
 real(dp),allocatable                 :: tempComponent(:)         ! temperature component of the energy increment (J m-3)
 real(dp),allocatable                 :: phseComponent(:)         ! phase component of the energy increment (J m-3)
 real(dp),allocatable                 :: inflComponent(:)         ! infiltration component of the energy increment (J m-3)
 ! define local variables
 character(len=256)                   :: cmessage                 ! error message of downwind routine
 integer(i4b)                         :: nLevels                  ! number of layers to consider in soil hydrology (e.g., # layers above water table)
 logical(lgt)                         :: freeze_infiltrate        ! .true to freeze infiltrating flux of liquid water
 logical(lgt)                         :: printflag                ! .true. if want to print
 integer(i4b)                         :: minLayer                 ! minimum layer to print
 integer(i4b)                         :: maxLayer                 ! maximum layer to print
 integer(i4b)                         :: nSnow                    ! number of snow layers
 integer(i4b)                         :: nSoil                    ! number of soil layers
 integer(i4b)                         :: iter                     ! iteration index
 integer(i4b),dimension(1)            :: liquid_pos               ! position of maximum error
 integer(i4b),dimension(1)            :: matric_pos               ! position of maximum error
 integer(i4b),dimension(1)            :: energy_pos               ! position of maximum error
 real(dp),dimension(1)                :: liquid_max               ! maximum absolute change in volumetric liquid water content for a given iteration (-)
 real(dp),dimension(1)                :: matric_max               ! maximum absolute change in matric head for a given iteration (m)
 real(dp),dimension(1)                :: energy_max               ! maximum absolute change in energy for a given iteration (J m-3)
 real(dp)                             :: aquifr_max               ! absolute change in aquifer storage for a given iteration (m)
 real(dp)                             :: theta                    ! volumetric fraction of total water, liquid plus ice (-)
 real(dp),parameter                   :: eps   = 1.d-10           ! small increment used to define ice content at the freezing point
 real(dp)                             :: checkCalcs               ! check the aquifer root density calculations
 real(dp)                             :: nrgRequired              ! case of "snow without a layer": energy required to melt all the snow (J m-2)
 real(dp)                             :: nrgAvailable             ! case of "snow without a layer": energy available to melt the snow (J m-2)
 real(dp)                             :: snwDensity               ! case of "snow without a layer": snow density (kg m-3)
 real(dp)                             :: volFrac_water            ! total volumetric fraction of water (liquid water plus ice) 
 real(dp)                             :: drainablePorosity        ! drainable porosity (-)
 integer(i4b)                         :: iLayer                   ! loop through model layers
 ! define balance check variables
 real(dp)                             :: totalChange              ! total change in volumetric liquid water content over a layer 
 real(dp)                             :: phaseChange              ! change in volumetric liquid water content associated with phase change
 real(dp)                             :: flux_Change              ! change in volumetric liquid water content associated with fluxes
 real(dp)                             :: evap_Change              ! change in volumetric liquid water content associated with transpiration
 real(dp)                             :: qbaseChange              ! change in volumetric liquid water content associated with baseflow
 real(dp)                             :: ejectChange              ! change in volumetric liquid water content associated with ejection of water
 real(dp)                             :: balanceSoilWater0        ! total soil storage at the start of the step (kg m-2)
 real(dp)                             :: balanceSoilWater1        ! total soil storage at the end of the step (kg m-2)
 real(dp)                             :: balanceSoilInflux        ! input to the soil zone
 real(dp)                             :: balanceSoilBaseflow      ! output from the soil zone
 real(dp)                             :: balanceSoilDrainage      ! output from the soil zone
 real(dp)                             :: balanceSoilEjection      ! output from the soil zone
 real(dp)                             :: balanceSoilTranspiration ! output from the soil zone
 real(dp)                             :: balanceAquifer0          ! total aquifer storage at the start of the step (kg m-2)
 real(dp)                             :: balanceAquifer1          ! total aquifer storage at the end of the step (kg m-2)
 ! initialize error control
 err=0; message="picardSolv/"

 ! assign local pointers to the model index structures
 nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType         => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! assign pointers to snow parameters
 snowfrz_scale     => mpar_data%var(iLookPARAM%snowfrz_scale)  ! scaling parameter for the snow freezing curve (K-1)
 Fcapil            => mpar_data%var(iLookPARAM%Fcapil)         ! capillary retention as a fraction of the total pore volume (-)

 ! assign pointers to soil/veg parameters
 soilAlbedo        => mpar_data%var(iLookPARAM%soilAlbedo)      ! soil albedo (-)
 vGn_alpha         => mpar_data%var(iLookPARAM%vGn_alpha)       ! van Genutchen "alpha" parameter (m-1)
 vGn_n             => mpar_data%var(iLookPARAM%vGn_n)           ! van Genutchen "n" parameter (-)
 theta_sat         => mpar_data%var(iLookPARAM%theta_sat)       ! soil porosity (-)
 theta_res         => mpar_data%var(iLookPARAM%theta_res)       ! soil residual volumetric water content (-)
 specificYield     => mpar_data%var(iLookPARAM%specificYield)   ! specific yield (-)
 specificStorage   => mpar_data%var(iLookPARAM%specificStorage) ! specific storage coefficient (m-1)
 f_impede          => mpar_data%var(iLookPARAM%f_impede)        ! ice impedence factor (-)
 rootingDepth      => mpar_data%var(iLookPARAM%rootingDepth)    ! rooting depth (m)
 rootDistExp       => mpar_data%var(iLookPARAM%rootDistExp)     ! root distribution exponent (-)

 ! assign pointers to algorithmic control parameters
 wimplicit         => mpar_data%var(iLookPARAM%wimplicit)            ! weight assigned to start-of-step fluxes (-)
 relConvTol_liquid => mpar_data%var(iLookPARAM%relConvTol_liquid)    ! relative convergence tolerance for vol frac liq water (-)
 absConvTol_liquid => mpar_data%var(iLookPARAM%absConvTol_liquid)    ! absolute convergence tolerance for vol frac liq water (-)
 relConvTol_matric => mpar_data%var(iLookPARAM%relConvTol_matric)    ! relative convergence tolerance for matric head (-)
 absConvTol_matric => mpar_data%var(iLookPARAM%absConvTol_matric)    ! absolute convergence tolerance for matric head (m)
 relConvTol_energy => mpar_data%var(iLookPARAM%relConvTol_energy)    ! relative convergence tolerance for energy (-)
 absConvTol_energy => mpar_data%var(iLookPARAM%absConvTol_energy)    ! absolute convergence tolerance for energy (J m-3)
 relConvTol_aquifr => mpar_data%var(iLookPARAM%relConvTol_aquifr)    ! relative convergence tolerance for aquifer storage (-)
 absConvTol_aquifr => mpar_data%var(iLookPARAM%absConvTol_aquifr)    ! absolute convergence tolerance for aquifer storage (m)

 ! assign pointers to model variables that are constant over the simulation period
 vGn_m             => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)           ! van Genutchen "m" parameter (-)
 mLayerRootDensity => mvar_data%var(iLookMVAR%mLayerRootDensity)%dat        ! fraction of roots in each soil layer (-)

 ! assign local pointers to scalar state variables
 scalarSWE            => mvar_data%var(iLookMVAR%scalarSWE)%dat(1)             ! SWE (kg m-2)
 surfaceAlbedo        => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)          ! surface albedo (-)
 scalarSfcMeltPond    => mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)     ! ponded water caused by melt of the "snow without a layer" (kg m-2)

 ! assign local pointers to diagnostic scalar variables
 scalarRainPlusMelt   => mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)    ! rain plus melt (m s-1)
 scalarSnowDepth      => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)       ! total snow depth (m)
 scalarPotentialET    => mvar_data%var(iLookMVAR%scalarPotentialET)%dat(1)     ! potential ET (kg m-2 s-1)
 scalarMassLiquid     => mvar_data%var(iLookMVAR%scalarMassLiquid)%dat(1)      ! transpiration (kg m-2 s-1)
 scalarMassSolid      => mvar_data%var(iLookMVAR%scalarMassSolid)%dat(1)       ! sublimation/frost (kg m-2 s-1)
 scalarSoilInflux     => mvar_data%var(iLookMVAR%scalarSoilInflux)%dat(1)      ! influx of water at the top of the soil profile (m s-1)
 scalarSoilBaseflow   => mvar_data%var(iLookMVAR%scalarSoilBaseflow)%dat(1)    ! baseflow from throughout the soil profile (m s-1)
 scalarSoilDrainage   => mvar_data%var(iLookMVAR%scalarSoilDrainage)%dat(1)    ! drainage from the bottom of the soil profile (m s-1)
 scalarSoilEjection   => mvar_data%var(iLookMVAR%scalarSoilEjection)%dat(1)    ! water ejected from soil layers (m s-1)
 scalarSoilWatBalError=> mvar_data%var(iLookMVAR%scalarSoilWatBalError)%dat(1) ! error in the total soil water balance (kg m-2)
 scalarTotalSoilLiq   => mvar_data%var(iLookMVAR%scalarTotalSoilLiq)%dat(1)    ! total mass of liquid water in the soil (kg m-2)
 scalarTotalSoilIce   => mvar_data%var(iLookMVAR%scalarTotalSoilIce)%dat(1)    ! total mass of ice in the soil (kg m-2)

 ! assign pointers to model coordinate variables
 iLayerHeight      => mvar_data%var(iLookMVAR%iLayerHeight)%dat             ! height of layer interfaces (m)

 ! assign pointers to model state variables -- all layers
 mLayerTemp        => mvar_data%var(iLookMVAR%mLayerTemp)%dat               ! temperature of each layer (K)
 mLayerDepth       => mvar_data%var(iLookMVAR%mLayerDepth)%dat              ! depth of each layer (m)
 mLayerVolFracIce  => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat         ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq  => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat         ! volumetric fraction of liquid water in each layer (-)

 ! assign pointers to model state variables -- soil layers
 mLayerMatricHead  => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat         ! matric head in each **soil** layer (m)

 ! assign pointers to model state variables -- aquifer
 scalarAquiferStorage  => mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1)  ! relative aquifer storage above the bottom of the soil profile (m)
 scalarWaterTableDepth => mvar_data%var(iLookMVAR%scalarWaterTableDepth)%dat(1) ! water table depth (m)
 
 ! assign pointers to model diagnostic variables -- all layers
 mLayerVolFracAir  => mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat         ! volumetric fraction of air in each layer (-)
 mLayerVolHtCapBulk=> mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat       ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerdTheta_dTk  => mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat         ! derivative in the freezing curve (K-1)
 mLayerMeltFreeze  => mvar_data%var(iLookMVAR%mLayerMeltFreeze)%dat         ! melt/freeze in each layer (kg m-3 s-1)
 mLayerInfilFreeze => mvar_data%var(iLookMVAR%mLayerInfilFreeze)%dat        ! increase in volumetric ice content caused by freezing infiltrating flux (-)

 ! assign pointers to model diagnostic variables -- soil only
 mLayerTcrit           => mvar_data%var(iLookMVAR%mLayerTcrit)%dat            ! critical soil temperature above which all water is unfrozen (K)
 mLayerdTheta_dPsi     => mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat      ! derivative in the soil water characteristic (m-1)
 iLayerInitLiqFluxSoil => mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat   ! (soil only) ! liquid flux at layer interfaces at the start of the time step (m s-1)
 iLayerLiqFluxSoil     => mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat       ! (soil only) ! liquid flux at layer interfaces at the end of the time step (m s-1)
 mLayerInitBaseflow    => mvar_data%var(iLookMVAR%mLayerInitBaseflow)%dat      ! (soil only) ! baseflow from each soil layer at start-of-step (m s-1)
 mLayerBaseflow        => mvar_data%var(iLookMVAR%mLayerBaseflow)%dat          ! (soil only) ! baseflow from each soil layer (m s-1)
 mLayerInitTranspire   => mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat     ! (soil only) ! transpiration loss from each soil layer at start-of-step (m s-1)
 mLayerTranspire       => mvar_data%var(iLookMVAR%mLayerTranspire)%dat         ! (soil only) ! transpiration loss from each soil layer (m s-1)
 mLayerInitEjectWater  => mvar_data%var(iLookMVAR%mLayerInitEjectWater)%dat    ! (soil only) ! water ejected from each soil layer at start-of-step (m s-1)
 mLayerEjectWater      => mvar_data%var(iLookMVAR%mLayerEjectWater)%dat        ! (soil only) ! water ejected from each soil layer (m s-1)
 iLayerSatHydCond      => mvar_data%var(iLookMVAR%iLayerSatHydCond)%dat        ! (soil only) ! saturated hydraulic conductivity at layer interfaces (m s-1)

 ! assign pointers to model diagnostic variables -- aquifer only 
 scalarAquiferRootFrac      => mvar_data%var(iLookMVAR%scalarAquiferRootFrac)%dat(1)   ! fraction of roots below the lowest unsatured layer (-)
 scalarInitAquiferTranspire => mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1)  ! initial transpiration from the aquifer (m s-1)
 scalarInitAquiferRecharge  => mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1)   ! initial aquifer recharge rate (m s-1)
 scalarInitAquiferBaseflow  => mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1)   ! initial aquifer baseflow rate (m s-1)
 scalarAquiferTranspire     => mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1)  ! transpiration from the aquifer (m s-1)
 scalarAquiferRecharge      => mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1)   ! aquifer recharge rate (m s-1)
 scalarAquiferBaseflow      => mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1)   ! aquifer baseflow rate (m s-1)
 scalarAquiferBalError      => mvar_data%var(iLookMVAR%scalarAquiferBalError)%dat(1)   ! error in the aquifer water balance (kg m-2)

 ! initialize print flag
 printflag=.false.

 ! set the freeze_infiltrate flag
 freeze_infiltrate = .true.

 ! define the maximum number of layers to print
 minLayer=40
 maxLayer=50

 ! allocate space for state variables at the start and end of the iteration
 allocate(mLayerTempIter(nLayers),      mLayerTempNew(nLayers),       &  ! all layers
          mLayerVolFracIceIter(nLayers),mLayerVolFracIceNew(nLayers), &  ! all layers
          mLayerVolFracLiqIter(nLayers),mLayerVolFracLiqNew(nLayers), &  ! all layers
          mLayerMatricHeadIter(nSoil),  mLayerMatricHeadNew(nSoil),   &  ! **soil layers only**
          stat=err)
 if(err/=0)then; err=20; message='problem allocating space for state variable vectors - 1'; return; endif

 ! allocate space for diagnostic variables
 allocate(mLayerInfilFreezeNew(nLayers),stat=err)
 if(err/=0)then; err=20; message='problem allocating space for diagnostic variable vectors - 1'; return; endif

 ! allocate space for error monitoring
 allocate(mLayerTempIncr(nLayers),    &  ! (all layers)
          mLayerTempIncrOld(nLayers), &  ! (all layers)
          tempComponent(nLayers),     &  ! (all layers)
          phseComponent(nLayers),     &  ! (all layers)
          inflComponent(nLayers),     &  ! (all layers)
          mLayerMatIncr(nSoil),       &  ! NOTE: nSoil
          mLayerLiqIncr(nLayers),     &  ! (all layers)
          mLayerNrgIncr(nLayers),     &  ! (all layers)
          mLayerNrgError(nLayers),    &  ! (all layers)
          stat=err)
 if(err/=0)then; err=20; message='problem allocating space for error monitoring'; return; endif

 ! initialize number of iterations
 niter=0

 if(printflag) print*, '************************************************************************'
 if(printflag) print*, '************************************************************************'
 if(printflag) print*, '************************************************************************'
 if(printflag) print*, 'layer type = ', layerType(minLayer:min(maxLayer,nLayers))
 !print*, 'initial temperature = ',          mLayerTemp(minLayer:min(maxLayer,nLayers))
 !print*, 'initial ice content = ',          mLayerVolFracIce(minLayer:min(maxLayer,nLayers))
 !print*, 'initial liquid water content = ', mLayerVolFracLiq(minLayer:min(maxLayer,nLayers))
 !print*, 'initial matric head =          ', mLayerMatricHead(minLayer:min(maxLayer,nSoil))
 !pause ' start of picardSolv'

 ! compute total soil moisture and ice (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
 scalarTotalSoilIce = sum(iden_ice  *mLayerVolFracIce(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))

 ! get the total water in the soil (liquid plus ice) at the start of the time step (kg m-2)
 balanceSoilWater0 = scalarTotalSoilLiq + scalarTotalSoilIce
 !print*, 'start of iteration, liquid, ice = ', scalarTotalSoilLiq, scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer0 = scalarAquiferStorage*iden_water

 ! initialize temperature
 mLayerTempIter = mLayerTemp

 ! initialize aquifer storage
 scalarAquiferStorageIter = scalarAquiferStorage

 ! identify the number of soil layers to use in the soil hydrology routine
 nLevels = nSoil  ! NOTE: always pass the full number of soil layers

 ! compute the fraction of roots below layers that are completely unsaturated (-)
 scalarAquiferRootFrac = 1._dp - sum(mLayerRootDensity(1:nLevels))

 ! check the aquifer root fraction is OK
 checkCalcs = 1._dp - ( min(iLayerHeight(nSnow+nLevels),rootingDepth) / rootingDepth)**rootDistExp
 if(abs(checkCalcs - scalarAquiferRootFrac) > epsilon(checkCalcs))then; err=20; message=trim(message)//'problem with the aquifer root density calculations'; return; endif

 ! compute the surface albedo (constant over the iterations)
 if(nSnow > 0)then
  call surfAlbedo(dt,&          ! input: time step (seconds)
                  err,cmessage) ! output: error control
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 else
  surfaceAlbedo = soilAlbedo
 endif

 ! compute diagnostic energy variables (thermal conductivity and volumetric heat capacity)
 ! (constant over the iterations)
 call diagn_evar(err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! calculate the critical soil temperature above which all water is unfrozen (K)
 do iLayer=nSnow+1,nSoil
  theta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)
  mLayerTcrit(iLayer-nSnow) = crit_soilT(theta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
 end do

 ! compute phase change "not accounted for" in the previous time step
 ! NOTE: input="start-of-step values" and output="Iter" to save extra copying
 call phsechange(mLayerTempIter,      & ! intent(in): new temperature vector (K)
                 mLayerMatricHead,    & ! intent(in): matric head at the current iteration (m)
                 mLayerVolFracLiq,    & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                 mLayerVolFracIce,    & ! intent(in): volumetric fraction of ice at the current iteration (-)
                 mLayerMatricHeadIter,& ! intent(out): new matric head (m)
                 mLayerVolFracLiqIter,& ! intent(out): new volumetric fraction of liquid water (-)
                 mLayerVolFracIceIter,& ! intent(out): new volumetric fraction of ice (-)
                 err,cmessage)          ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 
 ! initialize the melt/freeze vectors
 mLayerMeltFreeze  = 0._dp
 mLayerInfilFreeze = 0._dp

 !if(any(mLayerVolFracIce*iden_ice > 700._dp)) printflag=.true.

 ! ***** iterate
 do iter=1,maxiter

  ! increment number of iterations
  niter=niter+1

  !print*, '***********************************************************************'
  !print*, '***********************************************************************'

  !print*, 'iter = ', iter
  !print*, 'before heatTransf: mLayerTempIter(minLayer:min(maxLayer,nLayers)) = ',       mLayerTempIter(minLayer:min(maxLayer,nLayers))
  !print*, 'before heatTransf: mLayerVolFracLiqIter(minLayer:min(maxLayer,nLayers)) = ', mLayerVolFracLiqIter(minLayer:min(maxLayer,nLayers))
  !print*, 'before heatTransf: mLayerVolFracIceIter(minLayer:min(maxLayer,nLayers)) = ', mLayerVolFracIceIter(minLayer:min(maxLayer,nLayers))

  ! compute the temperature and ice content at the next iteration
  call heatTransf(dt,&                        ! time step (seconds)
                  iter,&                      ! current iteration count
                  mLayerTempIter,           & ! trial temperature at the current iteration (K)
                  mLayerVolFracIceIter,     & ! volumetric fraction of ice at the current iteration (-)
                  mLayerVolFracLiqIter,     & ! volumetric fraction of liquid water at the current iteration (-)
                  mLayerMatricHeadIter,     & ! matric head at the current iteration (m)
                  scalarAquiferStorageIter, & ! aquifer storage at the current iteration (m)
                  mLayerTempIncrOld,        & ! iteration increment for temperature from the previous iteration (K)
                  mLayerTempIncr,           & ! iteration increment for temperature (K)
                  mLayerTempNew,            & ! new temperature (K)
                  mLayerVolFracIceNew,      & ! new volumetric fraction of ice (-)
                  mLayerVolFracLiqNew,      & ! new volumetric fraction of liquid water (-)
                  mLayerMatricHeadNew,      & ! new matric head (m)
                  err,cmessage)               ! error control
  ! negative error code requires convergence check, so just check positive errors
  if(err>0)then; message=trim(message)//trim(cmessage); return; endif
  !print*, '*** after heatTransf'
  !write(*,'(a,10(f10.5,1x))') 'mLayerTempDiff =       ', mLayerTempNew(minLayer:min(maxLayer,nLayers)) - mLayerTempIter(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerTempIter =       ', mLayerTempIter(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerTempNew =        ', mLayerTempNew(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerVolFracIceIter = ', mLayerVolFracIceIter(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerVolFracIceNew =  ', mLayerVolFracIceNew(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerVolFracLiqIter = ', mLayerVolFracLiqIter(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerVolFracLiqNew =  ', mLayerVolFracLiqNew(minLayer:min(maxLayer,nLayers))

  if(printflag) print*, 'mLayerTempIter =       ', mLayerTempIter(minLayer:min(maxLayer,nLayers))
  if(printflag) print*, 'mLayerTempNew =        ', mLayerTempNew(minLayer:min(maxLayer,nLayers))
  if(printflag) print*, 'mLayerTempDiff =       ', mLayerTempNew(minLayer:min(maxLayer,nLayers)) - mLayerTempIter(minLayer:min(maxLayer,nLayers))
  if(printflag) print*, 'mLayerMatricHeadIter = ', mLayerMatricHeadIter(minLayer:min(maxLayer,nSoil))
  if(printflag) print*, 'mLayerMatricHeadNew = ',  mLayerMatricHeadNew(minLayer:min(maxLayer,nSoil))
  if(printflag) print*, 'mLayerVolFracLiqIter = ', mLayerVolFracLiqIter(minLayer:min(maxLayer,nLayers))
  if(printflag) print*, 'mLayerVolFracLiqNew =  ', mLayerVolFracLiqNew(minLayer:min(maxLayer,nLayers))
  if(printflag) print*, 'mLayerVolFracIceIter = ', mLayerVolFracIceIter(minLayer:min(maxLayer,nLayers))
  if(printflag) print*, 'mLayerVolFracIceNew =  ', mLayerVolFracIceNew(minLayer:min(maxLayer,nLayers))
  if(printflag) print*, 'after heatTransf'

  ! compute melt/freeze in each layer (kg m-3 s-1) -- melt is negative
  mLayerMeltFreeze = mLayerMeltFreeze + iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)/dt

  ! compute the temperature and phase components of the energy increment (J m-3)
  tempComponent = mLayerVolHtCapBulk*mLayerTempIncr
  phseComponent = iden_ice*LH_fus*(mLayerVolFracIceNew - mLayerVolFracIceIter)

  ! save the temperature increment
  mLayerTempIncrOld = mLayerTempIncr

  ! update the state variables
  mLayerTempIter       = mLayerTempNew

  ! only account for phase change after the first iteration
  !if(iter > 1)then
   mLayerMatricHeadIter = mLayerMatricHeadNew
   mLayerVolFracLiqIter = mLayerVolFracLiqNew
   !mLayerVolFracIceIter = mLayerVolFracIceNew
  !endif

  !print*, 'after heat transfer, nSnow, nSoil, nLayers, ice(nLayers) = ', nSnow, nSoil, nLayers, mLayerVolFracIceIter(nLayers)

  ! compute the volumetric liquid water content at the next iteration (note: only use snow vectors)
  ! NOTE: ice not modified in the snow hydrology routines, so can stay as "New"
  if(nSnow > 0)then
   call snowHydrol(dt,                               & ! time step (seconds)
                   iter,                             & ! iteration index
                   mLayerVolFracLiqIter(1:nSnow),    & ! volumetric fraction of liquid water at the current iteration (-)
                   mLayerVolFracIceNew(1:nSnow),     & ! volumetric fraction of ice at the current iteration (-)
                   mLayerVolFracLiqNew(1:nSnow),     & ! volumetric fraction of liquid water at the next iteration (-)
                   err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! compute the matric head at the next iteration (note liquid water and ice vectors are defined for all layers)
  ! NOTE: ice not modified in the soil hydrology routines, so can stay as "New"
  call soilHydrol(&
                  ! input
                  dt,&                                          ! time step (seconds)
                  iter,&                                        ! current iteration count
                  mLayerMatricHeadIter(1:nLevels),            & ! matric head in each layer at the current iteration (m)
                  mLayerVolFracLiqIter(nSnow+1:nSnow+nLevels),& ! volumetric fraction of liquid water at the current iteration (-)
                  mLayerVolFracIceNew(nSnow+1:nSnow+nLevels), & ! volumetric fraction of ice at the current iteration (-)
                  scalarAquiferStorageIter,                   & ! aquifer storage (m)
                  !output
                  mLayerMatricHeadNew(1:nLevels),             & ! matric head in each layer at the next iteration (m)
                  mLayerVolFracLiqNew(nSnow+1:nSnow+nLevels), & ! volumetric fraction of liquid water at the next iteration (-)
                  scalarAquiferStorageNew,                    & ! aquifer storage (m)
                  err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! copy across "saturated" layers
  if(nLevels < nSoil)then
   mLayerMatricHeadNew(nLevels+1:nSoil)         = mLayerMatricHeadIter(nLevels+1:nSoil)
   mLayerVolFracLiqNew(nSnow+nLevels+1:nLayers) = mLayerVolFracLiqIter(nSnow+nLevels+1:nLayers)
  endif
  !print*, '*** after hydrology'
  !write(*,'(a,50(e20.10,1x))') 'mLayerVolFracLiqIter = ', mLayerVolFracLiqIter(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,50(e20.10,1x))') 'mLayerVolFracLiqNew =  ', mLayerVolFracLiqNew(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,50(e20.10,1x))') 'mLayerMatricHeadIter = ', mLayerMatricHeadIter(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,50(e20.10,1x))') 'mLayerMatricHeadNew =  ', mLayerMatricHeadNew(minLayer:min(maxLayer,nLayers))
  !write(*,'(a, 1(e20.10,1x))') 'scalarAquiferStorageNew = ', scalarAquiferStorageNew

  ! compute the iteration increment for the matric head and volumetric fraction of liquid water
  mLayerMatIncr = mLayerMatricHeadNew - mLayerMatricHeadIter 
  mLayerLiqIncr = mLayerVolFracLiqNew - mLayerVolFracLiqIter 

  ! compute the iteration increment for aquifer storage
  scalarAqiIncr = scalarAquiferStorageNew - scalarAquiferStorageIter

  ! calculate the critical soil temperature above which all water is unfrozen (K)
  do iLayer=nSnow+1,nLayers
   theta = mLayerVolFracIceNew(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqNew(iLayer)
   mLayerTcrit(iLayer-nSnow) = crit_soilT(theta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
  end do

  ! option to compute phase change associated with infiltrating liquid water
  if(freeze_infiltrate)then
   ! compute phase change associated with infiltrating liquid water
   ! NOTE: input="New" and output="Iter" to save extra copying and get ready for the next iteration
   call phsechange(mLayerTempIter,      & ! intent(in): new temperature vector (K)
                   mLayerMatricHeadNew, & ! intent(in): matric head at the current iteration (m)
                   mLayerVolFracLiqNew, & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                   mLayerVolFracIceNew, & ! intent(in): volumetric fraction of ice at the current iteration (-)
                   mLayerMatricHeadIter,& ! intent(out): new matric head (m)
                   mLayerVolFracLiqIter,& ! intent(out): new volumetric fraction of liquid water (-)
                   mLayerVolFracIceIter,& ! intent(out): new volumetric fraction of ice (-)
                   err,cmessage)          ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  else 
   mLayerMatricHeadIter = mLayerMatricHeadNew
   mLayerVolFracLiqIter = mLayerVolFracLiqNew
   mLayerVolFracIceIter = mLayerVolFracIceNew
  endif
  !print*, '*** after phase change'
  !write(*,'(a,10(f10.5,1x))') 'mLayerVolFracLiqIter = ', mLayerVolFracLiqNew(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerVolFracLiqNew =  ', mLayerVolFracLiqIter(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerVolFracIceIter = ', mLayerVolFracIceNew(minLayer:min(maxLayer,nLayers))
  !write(*,'(a,10(f10.5,1x))') 'mLayerVolFracIceNew =  ', mLayerVolFracIceIter(minLayer:min(maxLayer,nLayers))

  ! compute melt/freeze of infiltrating liquid water in each layer (kg m-3 s-1) -- melt is negative
  mLayerInfilFreeze = mLayerInfilFreeze + iden_ice*(mLayerVolFracIceIter - mLayerVolFracIceNew)/dt

  ! compute increment in energy -- do here because of phase change (note: "iter" and "new" are reversed)
  inflComponent = iden_ice*LH_fus*(mLayerVolFracIceIter - mLayerVolFracIceNew)
  mLayerNrgIncr = tempComponent - phseComponent - inflComponent
  !print*, '*** nrg change'
  !write(*,'(a,1x,10(e20.10,1x))') 'tempComponent(1:9) = ', tempComponent(1:9)
  !write(*,'(a,1x,10(e20.10,1x))') 'phseComponent(1:9) = ', phseComponent(1:9)
  !write(*,'(a,1x,10(e20.10,1x))') 'inflComponent(1:9) = ', inflComponent(1:9)

  ! get ready for the next iteration
  scalarAquiferStorageIter = scalarAquiferStorageNew

  ! non-iterative check (do not expect convergence)
  if(maxiter==1) exit  ! NOTE: exit loop here to avoid return statement with error code

  ! compute maximum iteration increment
  liquid_max = maxval(abs(mLayerLiqIncr))
  matric_max = maxval(abs(mLayerMatIncr))
  energy_max = maxval(abs(mLayerNrgIncr))
  aquifr_max = abs(scalarAqiIncr)
  !print*, 'iter, maxiter = ', iter, maxiter
  !write(*,'(a,1x,2(e20.10,1x))') 'liquid_max, absConvTol_liquid = ', liquid_max, absConvTol_liquid
  !write(*,'(a,1x,2(e20.10,1x))') 'matric_max, absConvTol_matric = ', matric_max, absConvTol_matric
  !write(*,'(a,1x,2(e20.10,1x))') 'energy_max, absConvTol_energy = ', energy_max, absConvTol_energy
  !write(*,'(a,1x,1(e20.10,1x))') 'aquifr_max = ', aquifr_max
  !pause

  ! get position of maximum iteration increment
  liquid_pos = maxloc(abs(mLayerLiqIncr))
  matric_pos = maxloc(abs(mLayerMatIncr))
  energy_pos = maxloc(abs(mLayerNrgIncr))
  !print*, 'liquid_pos, matric_pos, energy_pos = ', liquid_pos, matric_pos, energy_pos
  !pause

  ! print current values
  !print*, 'mLayerTempIter       = ', mLayerTempIter(minLayer:min(maxLayer,nLayers))
  !print*, 'mLayerVolFracLiqIter = ', mLayerVolFracLiqIter(minLayer:min(maxLayer,nLayers))
  !print*, 'mLayerVolFracIceIter = ', mLayerVolFracIceIter(minLayer:min(maxLayer,nLayers))
  !print*, 'mLayerMatricHeadIter = ', mLayerMatricHeadIter(minLayer:min(maxLayer,nSoil))
  !print*, 'after phsechange'

  ! convergence check: 
  if(liquid_max(1) < absConvTol_liquid .and. &   ! volumetric fraction of liquid water (-)
     matric_max(1) < absConvTol_matric .and. &   ! matric head (m)
     energy_max(1) < absConvTol_energy .and. &   ! energy (J m-3)
     aquifr_max    < absConvTol_aquifr)      &   ! aquifer storage (m)
   exit

  ! check for lack of convergence
  if(niter==maxiter)then; err=-30; message=trim(message)//'failed to converge'; return; endif

 end do  ! (iterating)

 ! *****************************************************************************************************************************************
 ! *****************************************************************************************************************************************
 ! ***** END OF ITERATIONS *****************************************************************************************************************
 ! *****************************************************************************************************************************************
 ! *****************************************************************************************************************************************
 !pause 'after iterations'

 ! check matric head is not ridiculous
 do iLayer=1,nSoil
  if(mLayerMatricHeadIter(iLayer) > 100._dp)then
   write(message,'(a,i0,a,f9.1,a)')trim(message)//"matric head > 100 [iLayer=",iLayer,"; matricHead=",&
         mLayerMatricHeadIter(iLayer),"]"
   err=20; return
  endif
 end do

 !print*, 'check volumetric liquid water content...'
 !print*, mLayerVolFracLiqIter(nSnow+1)

 ! ** check that total volumetric water (liquid water plus ice) does not exceed porosity
 do iLayer=nSnow+1,nLayers
  ! compute total volumetric fraction filled with water (liquid plus ice)
  volFrac_water = mLayerVolFracIceIter(iLayer) + mLayerVolFracLiqIter(iLayer)
  ! check if the total volumetric water (liquid water plus ice) exceeds porosity
  if(volFrac_water > theta_sat)then
   write(*,'(a,i4,1x,2(f15.10,1x))') 'iLayer, mLayerVolFracIceIter(iLayer), mLayerVolFracLiqIter(iLayer) = ',&
                                      iLayer, mLayerVolFracIceIter(iLayer), mLayerVolFracLiqIter(iLayer)
   write(message,'(a,i0,a,i0,a)')trim(message)//"(liquid + ice) > porosity [iLayer=",iLayer,"; iSoil=",iLayer-nSnow,"]"
   err=20; return
  endif  ! (if the total volumetric water -- liquid water plus ice -- exceeds porosity)
 end do ! (looping through soil layers)

 !print*, 'check volumetric liquid water content again...'
 !print*, mLayerVolFracLiqIter(nSnow+1)

 ! compute sublimation
 if(nSnow>0)&
  mLayerVolFracIceIter(1) = mLayerVolFracIceIter(1) + scalarPotentialET*dt/(mLayerDepth(1)*iden_ice)

 ! compute change in snow density
 ! NOTE: input "iter" because of copying trick in phase change
 call snwDensify(dt,                  &  ! input:  time step (seconds)
                 mLayerVolFracLiqIter(1:nSnow),&  ! input:  volumetric fraction of liquid water after itertations (-)
                 mLayerVolFracIceIter(1:nSnow),&  ! input:  volumetric fraction of ice after itertations (-)
                 mLayerVolFracLiqNew(1:nSnow), &  ! output: volumetric fraction of liquid water after densification (-)
                 mLayerVolFracIceNew(1:nSnow), &  ! output: volumetric fraction of ice after densification (-)
                 err,cmessage)           ! output: error control
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! copy over state variables for soil layers
 mLayerMatricHeadNew                  = mLayerMatricHeadIter
 mLayerVolFracLiqNew(nSnow+1:nLayers) = mLayerVolFracLiqIter(nSnow+1:nLayers)
 mLayerVolFracIceNew(nSnow+1:nLayers) = mLayerVolFracIceIter(nSnow+1:nLayers)

 ! ***** compute melt for the case of "snow without a layer"
 if(nSnow==0 .and. scalarSWE > 0._dp)then
  ! only melt if temperature of the top soil layer is greater than Tfreeze
  if(mLayerTempNew(1) > Tfreeze)then
   ! compute the energy required to melt all the snow (J m-2)
   nrgRequired     = scalarSWE*LH_fus
   ! compute the energy available to melt the snow (J m-2)
   nrgAvailable    = mLayerVolHtCapBulk(1)*(mLayerTempNew(1) - Tfreeze)*mLayerDepth(1)
   ! compute the snow density (not saved)
   snwDensity      = scalarSWE/scalarSnowDepth
   ! compute the amount of melt, and update SWE (kg m-2)
   if(nrgAvailable > nrgRequired)then
    scalarSfcMeltPond  = scalarSWE
    scalarSWE          = 0._dp
   else
    scalarSfcMeltPond  = nrgAvailable/LH_fus
    scalarSWE          = scalarSWE - scalarSfcMeltPond
   endif
   ! update depth
   scalarSnowDepth = scalarSWE/snwDensity
   ! update temperature of the top soil layer (K)
   mLayerTempNew(1)= mLayerTempNew(1) - (scalarSfcMeltPond/mLayerDepth(1))/mLayerVolHtCapBulk(1)
  else  ! melt is zero if the temperature of the top soil layer is less than Tfreeze
   scalarSfcMeltPond = 0._dp  ! kg m-2
  endif ! (if the temperature of the top soil layer is greater than Tfreeze)
 else  ! melt is zero if the "snow without a layer" does not exist
  scalarSfcMeltPond = 0._dp  ! kg m-2
 endif ! (if the "snow without a layer" exists)

 ! compute total soil moisture and ice (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiqNew(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
 scalarTotalSoilIce = sum(iden_ice  *mLayerVolFracIceNew(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))

 ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
 balanceSoilWater1 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer1 = scalarAquiferStorageNew*iden_water

 ! compute the influx, drainage and ejection of water from the soil profile (m s-1)
 scalarSoilInflux       = (wimplicit*iLayerInitLiqFluxSoil(0)       + (1._dp - wimplicit)*iLayerLiqFluxSoil(0)      )
 scalarSoilBaseflow     = (wimplicit*sum(mLayerInitBaseflow)        + (1._dp - wimplicit)*sum(mLayerBaseflow)       )
 scalarSoilDrainage     = (wimplicit*iLayerInitLiqFluxSoil(nLevels) + (1._dp - wimplicit)*iLayerLiqFluxSoil(nLevels))
 scalarSoilEjection     = (wimplicit*sum(mLayerInitEjectWater)      + (1._dp - wimplicit)*sum(mLayerEjectWater)     )

 ! check ejected water
 if(scalarSoilEjection < 0._dp)then
  print*, 'mLayerInitEjectWater = ', mLayerInitEjectWater
  print*, 'mLayerEjectWater = ', mLayerEjectWater
  message=trim(message)//'ejected water < 0'; err=20; return
 endif

 ! get the input and output to/from the soil zone (kg m-2)
 balanceSoilInflux        = scalarSoilInflux*iden_water*dt
 balanceSoilBaseflow      = scalarSoilBaseflow*iden_water*dt
 balanceSoilDrainage      = scalarSoilDrainage*iden_water*dt
 balanceSoilEjection      = scalarSoilEjection*iden_water*dt
 balanceSoilTranspiration = scalarMassLiquid*dt - (wimplicit*scalarInitAquiferTranspire + (1._dp - wimplicit)*scalarAquiferTranspire)*iden_water*dt

 ! check the soil water balance
 scalarSoilWatBalError  = balanceSoilWater1 - (balanceSoilWater0 + (balanceSoilInflux + balanceSoilTranspiration - balanceSoilBaseflow - balanceSoilDrainage - balanceSoilEjection) )
 if(abs(scalarSoilWatBalError) > 1.d-3)then
  ! check the balance of each layer
  write(*,'(a)') 'water balance of each layer'
  write(*,'(a)') 'Temp0 (K), Temp1(K), Liq0 (-), Liq1 (-), Ice0 (-), Ice1 (-), totalChange (-), phaseChange (-), flux_Change, evap_Change, ejectChange, ',&
                 'phaseChange+flux_Change+evap_Change+ejectChange, totalChange - (phaseChange+flux_Change+evap_Change-ejectChange)'
  do iLayer=1,nSoil
   totalChange = mLayerVolFracLiqNew(iLayer+nSnow) - mLayerVolFracLiq(iLayer+nSnow) ! total change in volumetric liquid water content
   phaseChange = -(iden_ice/iden_water)*(mLayerVolFracIceNew(iLayer+nSnow) - mLayerVolFracIce(iLayer+nSnow))  ! change in liquid water content associated with freezing
   evap_Change = dt*mLayerTranspire(iLayer)/mLayerDepth(iLayer)
   qbaseChange = dt*mLayerBaseflow(iLayer)/mLayerDepth(iLayer)
   ejectChange = dt*mLayerEjectWater(iLayer)/mLayerDepth(iLayer)
   flux_Change = dt*(iLayerLiqFluxSoil(iLayer-1) - iLayerLiqFluxSoil(iLayer))/mLayerDepth(iLayer) ! change in volumetric liquid water content from the interface fluxes
   write(*,'(i4,1x,2(f15.8,1x),20(e15.5,1x))') iLayer, mLayerTemp(iLayer+nSnow), mLayerTempNew(iLayer+nSnow), &
                                                        mLayerVolFracLiq(iLayer+nSnow), mLayerVolFracLiqNew(iLayer+nSnow), &
                                                        mLayerVolFracIce(iLayer+nSnow), mLayerVolFracIceNew(iLayer+nSnow), &
                                                        totalChange, phaseChange, flux_Change, evap_Change, qbaseChange, ejectChange, &
                                                        phaseChange+flux_Change+evap_Change-qbaseChange-ejectChange, totalChange - (phaseChange+flux_Change+evap_Change-qbaseChange-ejectChange)
  end do
  ! print the total water balance
  print*, 'dt = ', dt
  print*, 'balanceSoilWater0 (kg m-2) = ',        balanceSoilWater0
  print*, 'balanceSoilWater1 (kg m-2) = ',        balanceSoilWater1
  print*, 'balanceSoilInflux (kg m-2) = ',        balanceSoilInflux
  print*, 'balanceSoilDrainage (kg m-2) = ',      balanceSoilDrainage
  print*, 'balanceSoilEjection (kg m-2) = ',      balanceSoilEjection
  print*, 'balanceSoilTranspiration (kg m-2) = ', balanceSoilTranspiration, sum(mLayerTranspire)*dt*iden_water
  write(message,'(a,e20.10,a)')trim(message)//"abs(scalarSoilWatBalError) > 1.d-3 [error = ",&
                               scalarSoilWatBalError," ]"
  err=10; return
 endif

 ! check the aquifer water balance
 select case(model_decisions(iLookDECISIONS%groundwatr)%iDecision)
  ! no explicit aquifer
  case(noExplicit,equilWaterTable)
   scalarAquiferBalError = 0._dp
  ! explicit aquifer
  case(bigBucket,pseudoWaterTable)
   ! check the aquifer water balance
   scalarAquiferBalError = balanceAquifer1 - (balanceAquifer0 + (iden_water*(wimplicit*scalarInitAquiferTranspire + (1._dp - wimplicit)*scalarAquiferTranspire)*dt) + &
                                                                (iden_water*(wimplicit*scalarInitAquiferRecharge  + (1._dp - wimplicit)*scalarAquiferRecharge) *dt) - &
                                                                (iden_water*(wimplicit*scalarInitAquiferBaseflow  + (1._dp - wimplicit)*scalarAquiferBaseflow) *dt) )
   ! print the terms in the aquifer balance if errors are sufficiently large
   if(abs(scalarAquiferBalError) > 1.d-3)then
    write(*,'(a,f20.10)') 'scalarAquiferBalError  = ', scalarAquiferBalError
    write(*,'(a,f20.10)') 'balanceAquifer1        = ', balanceAquifer1
    write(*,'(a,f20.10)') 'balanceAquifer0        = ', balanceAquifer0
    write(*,'(a,f20.10)') 'scalarAquiferTranspire = ', iden_water*(wimplicit*scalarInitAquiferTranspire + (1._dp - wimplicit)*scalarAquiferTranspire)*dt
    write(*,'(a,f20.10)') 'scalarAquiferRecharge  = ', iden_water*(wimplicit*scalarInitAquiferRecharge  + (1._dp - wimplicit)*scalarAquiferRecharge) *dt
    write(*,'(a,f20.10)') 'scalarAquiferBaseflow  = ', iden_water*(wimplicit*scalarInitAquiferBaseflow  + (1._dp - wimplicit)*scalarAquiferBaseflow) *dt
    write(message,'(a,e20.10,a)')trim(message)//"abs(scalarAquiferBalError) > 1.d-3 [error = ",scalarAquiferBalError," ]"
    err=10; return
   endif
  case default; err=20; message=trim(message)//'unknown groundwater parameterization'; return
 end select ! (selecting groundwater parameterization)

 ! update the state vectors
 mLayerTemp           = mLayerTempNew           ! New (after heatTransf)
 mLayerVolFracIce     = mLayerVolFracIceNew     ! New (after densification)
 mLayerVolFracLiq     = mLayerVolFracLiqNew     ! New (after densifcaction)
 mLayerMatricHead     = mLayerMatricHeadNew     ! New (after densification)
 scalarAquiferStorage = scalarAquiferStorageNew ! New (after groundwater)

 ! compuute total snow depth and SWE
 if(nSnow>0)then
  scalarSnowDepth = sum(mLayerDepth, mask = layerType==ix_snow)
  scalarSWE       = sum(iden_ice*  mLayerVolFracIce(1:nSnow)*mLayerDepth(1:nSnow)) + &  ! total ice
                    sum(iden_water*mLayerVolFracLiq(1:nSnow)*mLayerDepth(1:nSnow))
 endif

 ! deallocate space for state variables at the current iteration
 deallocate(mLayerTempIter,mLayerVolFracIceIter,mLayerMatricHeadIter,mLayerVolFracLiqIter,stat=err)
 if(err/=0)then; err=40; message='problem deallocating space for state variable vectors - 1'; return; endif

 ! deallocate space for state variables at the current iteration
 deallocate(mLayerTempNew,mLayerVolFracIceNew,mLayerMatricHeadNew,mLayerVolFracLiqNew,stat=err)
 if(err/=0)then; err=40; message='problem deallocating space for state variable vectors - 2'; return; endif

 ! deallocate space for diagnostic variables
 deallocate(mLayerInfilFreezeNew,stat=err)
 if(err/=0)then; err=40; message='problem deallocating space for diagnostic vectors'; return; endif

 ! deallocate space for error monitoring
 deallocate(tempComponent,phseComponent,inflComponent,mLayerTempIncr,mLayerTempIncrOld,mLayerLiqIncr,mLayerNrgIncr,mLayerNrgError,mLayerMatIncr,stat=err)
 if(err/=0)then; err=40; message='problem deallocating space for error monitoring vectors'; return; endif

 end subroutine picardSolv

end module picardSolv_module
