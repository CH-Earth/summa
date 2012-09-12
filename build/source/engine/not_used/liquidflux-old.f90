module liquidflux_module
USE nrtype
implicit none
private
public::liquidflux
contains

 ! ************************************************************************************************
 ! new subroutine: compute liquid water flux
 ! ************************************************************************************************
 subroutine liquidflux(err,message)
 USE multiconst,only:Cp_water,&   ! specific heat of water      (J kg-1 K-1)
                     Tfreeze, &   ! temperature of freezing     (K) 
                     LH_fus,  &   ! latent heat of fusion       (J kg-1)
                     LH_sub,  &   ! latent heat of sublimation  (J kg-1)
                     LH_vap       ! latent heat of vaporization (J kg-1)
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute energy fluxes within the domain
 implicit none
 ! dummy variables
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! define local pointers for model parameters
 real(dp),pointer              :: k_soil                   ! hydraulic condictivity of soil (kg m-2 s-1)
 real(dp),pointer              :: ClappH_bpar              ! Clapp-Hornberger "b" parameter (-)
 real(dp),pointer              :: matpot_sat               ! saturated soil suction (m)
 real(dp),pointer              :: soil_pore                ! soil porosity (-)
 real(dp),pointer              :: k_snow                   ! hydraulic conductivity of snow (kg m-2 s-1)
 real(dp),pointer              :: Fcapil                   ! capillary retention (fraction of total pore volume)
 real(dp),pointer              :: mw_exp                   ! exponent for meltwater flow through snow (-)
 ! define local pointers for model forcing variables
 real(dp),pointer              :: airtemp                  ! air temperature (K)
 ! define local pointers for model variables
 real(dp),pointer              :: rainfall                 ! computed rainfall rate (kg m-2 s-1)
 real(dp),pointer              :: latentHeat               ! latent heat flux at the surface (J m-2 s-1)
 real(dp),pointer              :: iLayerHeight(:)          ! height at the interface of each layer (m)
 real(dp),pointer              :: mLayerHeight(:)          ! height at the mid-point of each layer (m)
 real(dp),pointer              :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer              :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 real(dp),pointer              :: mLayerSolidPorosity(:)   ! solid porosity of each layer, voids between ice and dry solids (-)
 real(dp),pointer              :: iLayerDiffLiqFlux(:)     ! diffusive flux at layer interfaces (kg m-2 s-1)
 real(dp),pointer              :: iLayerGravLiqFlux(:)     ! gravitational flux at layer interfaces (kg m-2 s-1)
 real(dp),pointer              :: iLayerMassLiqFlux(:)     ! total liquid flux at layer interfaces (kg m-2 s-1)
 real(dp),pointer              :: iLayerEnthalpyLiqFlux(:) ! enthalpy associated with vertical flux of liquid water (J m-2 s-1)
 real(dp),pointer              :: mLayerLatDrainage(:)     ! lateral drainage (kg m-3 s-1)
 real(dp),pointer              :: mLayerLatEnthalpy(:)     ! enthalpy associated with lateral drainage (J m-3 s-1)
 real(dp),pointer              :: mLayerMassLiqFlux(:)     ! temporal derivative in liquid water storage (kg m-3 s-1)
 real(dp),pointer              :: mLayerEnthalpyLiqFlux(:) ! temporal derivative in enthalpy associated with vertical flux of liquid water (J m-3 s-1)
 ! define local pointers for model index variables
 integer(i4b),pointer          :: nLayers                 ! number of layers
 integer(i4b),pointer          :: layerType(:)            ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 real(dp),parameter            :: smoother=0.001_dp       ! amount of smoothing for relative saturation
 integer(i4b)                  :: iLayer                  ! index of model layer
 real(dp)                      :: iLayerVolFracLiq        ! volumetric fraction of liquid water at the layer interface (-)
 real(dp)                      :: iLayerSolidPorosity     ! soild porosity at the layer interface (-)
 real(dp)                      :: diff_coef               ! diffusion coefficient at a layer interface (m kg m-2 s-1)
 real(dp)                      :: RelSat                  ! relative saturation of snow layer (-)
 real(dp)                      :: sat_smooth              ! saturation smoothed using Kavetski-Kuczera logistic function (-)
 real(dp)                      :: frac_LatDrainage        ! fraction of lateral drainage in a given layer (-)
 real(dp)                      :: netMassLiqFlux          ! net flux of liquid water (kg m-3 s-1)
 real(dp)                      :: netEnthLiqFlux          ! enthalpy associated with the net liquid water flux (J m-3 s-1)
 
 ! initialize error control
 err=0; message="f-liquidflux/"

 ! assign pointers to model parameters
 k_soil                => mpar_data%var(iLookPARAM%k_soil)                   ! hydraulic conductivity of the soil (kg m-2 s-1)
 ClappH_bpar           => mpar_data%var(iLookPARAM%ClappH_bpar)              ! Clapp-Hornberger "b" parameter (-)
 matpot_sat            => mpar_data%var(iLookPARAM%matpot_sat)               ! saturated soil suction (m)
 soil_pore             => mpar_data%var(iLookPARAM%soil_pore)                ! soil porosity
 k_snow                => mpar_data%var(iLookPARAM%k_snow)                   ! hydraulic conductivity of snow (kg m-2 s-1)
 Fcapil                => mpar_data%var(iLookPARAM%Fcapil)                   ! capillary retention (fraction of total pore volume)
 mw_exp                => mpar_data%var(iLookPARAM%mw_exp)                   ! exponent for meltwater flow through snow (-)
 ! assign pointers to model forcing data
 airtemp               => forc_data%var(iLookFORCE%airtemp)                  ! air temperature (K)
 ! assign pointers to model variables
 rainfall              => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)     ! computed rainfall rate (kg m-2 s-1)
 LatentHeat            => mvar_data%var(iLookMVAR%scalarLatHeat)%dat(1)      ! latent heat flux at the surface (J m-2 s-1)
 iLayerHeight          => mvar_data%var(iLookMVAR%iLayerHeight)%dat          ! height at the interface of each layer (m)
 mLayerHeight          => mvar_data%var(iLookMVAR%mLayerHeight)%dat          ! height at the mid-point of each layer (m)
 mLayerTemp            => mvar_data%var(iLookMVAR%mLayerTemp)%dat            ! temperature of each layer (K)
 mLayerVolFracLiq      => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat      ! volumetric fraction of liquid water in each layer (-)
 mLayerSolidPorosity   => mvar_data%var(iLookMVAR%mLayerSolidPorosity)%dat   ! solid porosity of each layer, voids between ice and dry solids (-)
 iLayerDiffLiqFlux     => mvar_data%var(iLookMVAR%iLayerDiffLiqFlux)%dat     ! diffusive flux at layer interfaces (kg m-2 s-1)
 iLayerGravLiqFlux     => mvar_data%var(iLookMVAR%iLayerGravLiqFlux)%dat     ! gravitational flux at layer interfaces (kg m-2 s-1)
 iLayerMassLiqFlux     => mvar_data%var(iLookMVAR%iLayerMassLiqFlux)%dat     ! total liquid flux at layer interfaces (kg m-2 s-1)
 iLayerEnthalpyLiqFlux => mvar_data%var(iLookMVAR%iLayerEnthalpyLiqFlux)%dat ! enthalpy associated with vertical flux of liquid water (J m-2 s-1)
 mLayerLatDrainage     => mvar_data%var(iLookMVAR%mLayerLatDrainage)%dat     ! lateral drainage (kg m-3 s-1)
 mLayerLatEnthalpy     => mvar_data%var(iLookMVAR%mLayerLatEnthalpy)%dat     ! enthalpy associated with lateral drainage (J m-3 s-1)
 mLayerMassLiqFlux     => mvar_data%var(iLookMVAR%mLayerMassLiqFlux)%dat     ! temporal derivative in liquid water storage (kg m-3 s-1)
 mLayerEnthalpyLiqFlux => mvar_data%var(iLookMVAR%mLayerEnthalpyLiqFlux)%dat ! temporal derivative in enthalpy associated with vertical flux of liquid water (J m-3 s-1)
 ! assign pointers to index variables
 layerType             => indx_data%var(iLookINDEX%layerType)%dat            ! layer type (ix_soil or ix_snow)
 nLayers               => indx_data%var(iLookINDEX%nLayers)%dat(1)           ! number of layers

 ! NOTE: fluxes are positive downward

 !print*,'volFracLiq    = ',mLayerVolFracLiq
 !print*,'solidPorosity = ',mLayerSolidPorosity

 ! add rainfall to the upper boundary (kg m-2 s-1)
 iLayerDiffLiqFlux(nLayers) = 0._dp
 iLayerGravLiqFlux(nLayers) = rainfall
 iLayerMassLiqFlux(nLayers) = rainfall

 ! add the mass associated with the latent heat flux (kg m-2 s-1)
 if(mLayerTemp(nLayers)>Tfreeze) then ! unfrozen, use latent heat of vaporization
  iLayerMassLiqFlux(nLayers) = iLayerMassLiqFlux(nLayers) + LatentHeat/LH_vap
 else                                 ! frozen, use latent heat of sublimation
  iLayerMassLiqFlux(nLayers) = iLayerMassLiqFlux(nLayers) + LatentHeat/LH_sub
 endif

 ! compute the enthalpy associated with rainfall (J m-2 s-1)
 iLayerEnthalpyLiqFlux(nLayers) = (Cp_water*(airtemp-Tfreeze) + LH_fus) * rainfall

 ! iloop thru layers and compute the liquid water fluxes at the **bottom** of each layer
 do iLayer=1,nLayers
  select case(layerType(nLayers))
  ! **********************************************************************************************************
  ! (1) SOIL CASE: RICHARDS' EQUATION
  ! **********************************************************************************************************
  case(ix_soil)
   ! compute the volumetric water content and solid porosity at the bottom of the layer
   if(iLayer>1)then
    iLayerVolFracLiq    = 0.5_dp*(mLayerVolFracLiq(iLayer)    +mLayerVolFracLiq(iLayer-1))
    iLayerSolidPorosity = 0.5_dp*(mLayerSolidPorosity(iLayer) +mLayerSolidPorosity(iLayer-1))
   else
    iLayerVolFracLiq    = mlayerVolFracLiq(iLayer)
    iLayerSolidPorosity = mLayerSolidPorosity(iLayer)
   endif
   ! compute the diffusive flux at the bottom of the layer (kg m-2 s-1)
   if(iLayer>1)then
    ! check for super-saturation
    !if(iLayerVolFracLiq>iLayerSolidPorosity)stop 'FORTRAN STOP: super-saturated!'
    ! compute the diffusion coefficient at the bottom of the layer (m kg m-2 s-1)
    diff_coef = -k_soil*(ClappH_bpar*matpot_sat/iLayerSolidPorosity)*(iLayerVolFracLiq/iLayerSolidPorosity)**(ClappH_bpar+2._dp) 
    ! compute the diffusive flux at the bottom of the layer (kg m-2 s-1)
    iLayerDiffLiqFlux(ilayer-1) = diff_coef * ( (mLayerVolFracLiq(iLayer)-mLayerVolFracLiq(iLayer-1)) / (mLayerHeight(iLayer)-mLayerHeight(iLayer-1)) )
   else
    iLayerDiffLiqFlux(ilayer-1) = 0._dp
   endif
   ! compute the gravitational flux (kg m-2 s-1)
   iLayerGravLiqFlux(ilayer-1)     = k_soil * (iLayerVolFracLiq/iLayerSolidPorosity)**(2._dp*ClappH_bpar + 3._dp)
   ! compute the total liquid flux (kg m-2 s-1)
   iLayerMassLiqFlux(ilayer-1)     = iLayerGravLiqFlux(iLayer-1) + iLayerDiffLiqFlux(ilayer-1)
   ! compute the enthalpy associated with the liquid water flux (J m-2 s-1)
   iLayerEnthalpyLiqFlux(iLayer-1) = (Cp_water*(mLayerTemp(iLayer)-Tfreeze) + LH_fus) * iLayerMassLiqFlux(ilayer-1)
  ! **********************************************************************************************************
  ! (2) SNOW CASE: GRAVITY DRAINAGE
  ! **********************************************************************************************************
  case(ix_snow)
   ! compute relative saturation
   RelSat = (mLayerVolFracLiq(iLayer) - Fcapil) / (mLayerSolidPorosity(iLayer) - Fcapil)
   ! smooth relative saturation using logistic function from Kavetski and Kuczera (2007)
   sat_smooth = smoother * (RelSat/smoother + log(1._dp + exp(-(RelSat/smoother)))) 
   ! compute gravitational drainage at the BOTTOM of the layer (kg m-2 s-1)
   iLayerGravLiqFlux(iLayer-1)     = k_snow * sat_smooth**mw_exp
   ! set the diffusive flux to zero (kg m-2 s-1)
   iLayerDiffLiqFlux(ilayer-1)     = 0._dp
   ! compute the total liquid flux (kg m-2 s-1)
   iLayerMassLiqFlux(ilayer-1)     = iLayerGravLiqFlux(iLayer-1) + iLayerDiffLiqFlux(ilayer-1)
   ! compute the enthalpy associated with the liquid water flux (J m-2 s-1)
   iLayerEnthalpyLiqFlux(iLayer-1) = (Cp_water*(mLayerTemp(iLayer)-Tfreeze) + LH_fus) * iLayerMassLiqFlux(ilayer-1)
  ! **********************************************************************************************************
  case default; err=10; message=trim(message)//"unknownOption"; return
  end select
 end do  ! (looping through layers)
 ! ***********************************************************************************************************
 ! loop thru layers and compute the temporal derivative in model state variables
 ! ***********************************************************************************************************
 do iLayer=1,nLayers
  ! compute the fraction of lateral drainage
  if(mLayerVolFracLiq(iLayer)/mLayerSolidPorosity(iLayer) > 1._dp - smoother*20._dp*2._dp) then
   frac_LatDrainage = 1._dp / (1._dp + exp( (1._dp - mLayerVolFracLiq(iLayer)/mLayerSolidPorosity(iLayer) - smoother*20._dp) / smoother) )
  else
   frac_LatDrainage = 0._dp
  endif
  ! compute the net flux of liquid water (kg m-3 s-1)
  netMassLiqFlux = (iLayerMassLiqFlux(iLayer) - iLayerMassLiqFlux(ilayer-1)) / (iLayerHeight(iLayer) - iLayerHeight(iLayer-1))
  ! compute the enthalpy associated with the liquid water flux (J m-3 s-1)
  netEnthLiqFlux = (iLayerEnthalpyLiqFlux(iLayer) - iLayerEnthalpyLiqFlux(iLayer-1)) / (iLayerHeight(iLayer) - iLayerHeight(iLayer-1))
  ! compute the lateral drainage (kg m-3 s-1)
  mLayerLatDrainage(iLayer)     = netMassLiqFlux*frac_LatDrainage
  ! compute the enthalpy associated with lateral drainage (J m-3 s-1)
  mLayerLatEnthalpy(iLayer)     = netEnthLiqFlux*frac_LatDrainage
  ! compute the temporal derivative for the bulk density of liquid water (kg m-3 s-1)
  mLayerMassLiqFlux(iLayer)     = netMassLiqFlux - mLayerLatDrainage(iLayer)
  ! compute temporal derivative for the energy associated with the liquid water flux (J m-3 s-1)
  mLayerEnthalpyLiqFlux(iLayer) = netEnthLiqFlux !- mLayerLatEnthalpy(iLayer)
  if(mLayerVolFracLiq(iLayer)>soil_pore)then
   print*,iLayer
   print*,frac_LatDrainage
   print*,iLayerDiffLiqFlux(iLayer)
   print*,iLayerDiffLiqFlux(ilayer-1)
   print*,iLayerGravLiqFlux(iLayer)
   print*,iLayerGravLiqFlux(ilayer-1)
   print*,iLayerMassLiqFlux(iLayer)
   print*,iLayerMassLiqFlux(ilayer-1)
   print*,netMassLiqFlux
   print*,mLayerLatDrainage(iLayer)
   print*,mLayerMassLiqFlux(iLayer)
   stop ' in liquidflux'
  endif

 end do  ! (looping through layers)
 end subroutine liquidflux

end module liquidflux_module
