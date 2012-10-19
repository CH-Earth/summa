module var_derive_module
USE nrtype
implicit none
private
public::calcHeight
public::turbExchng
public::rootDensty
public::satHydCond
public::v_shortcut
contains

 ! **********************************************************************************************************
 ! new subroutine: compute snow height
 ! **********************************************************************************************************
 subroutine calcHeight(err,message)
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out) :: err               ! error code
 character(*),intent(out) :: message           ! error message
 ! declare pointers to data in model variable structures
 real(dp),pointer         :: mLayerDepth(:)    ! depth of the layer (m)
 real(dp),pointer         :: mLayerHeight(:)   ! height of the layer mid-point (m)
 real(dp),pointer         :: iLayerHeight(:)   ! height of the layer interface (m)
 ! declare pointers to model index variables
 integer(i4b),pointer     :: nLayers           ! number of layers
 integer(i4b),pointer     :: layerType(:)      ! type of the layer (ix_soil or ix_snow)
 ! declare local variables
 integer(i4b)             :: iLayer            ! loop through layers
 ! initialize error control
 err=0; message='calcHeight/'
 ! assign local pointers to the values in the model variable structures
 mLayerDepth    =>mvar_data%var(iLookMVAR%mLayerDepth)%dat             ! depth of the layer (m)
 mLayerHeight   =>mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height of the layer mid-point (m)
 iLayerHeight   =>mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height of the layer interface (m)
 ! assign local pointers to the model index structures
 nLayers        =>indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType      =>indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)
 ! ************************************************************************************************************************
 ! initialize layer height as the top of the snowpack -- positive downward
 iLayerHeight(0) = -sum(mLayerDepth, mask=layerType==ix_snow)
 ! loop through layers
 do iLayer=1,nLayers
  ! compute the height at the layer midpoint
  mLayerHeight(iLayer) = iLayerHeight(iLayer-1) + mLayerDepth(iLayer)/2._dp
  ! compute the height at layer interfaces
  iLayerHeight(iLayer) = iLayerHeight(iLayer-1) + mLayerDepth(iLayer)
 end do ! (looping through layers)
 !print*, 'layerType   = ',  layerType
 !print*, 'mLayerDepth = ',  mLayerDepth
 !print*, 'mLayerHeight = ', mLayerHeight
 !print*, 'iLayerHeight = ', iLayerHeight
 !print*, '************** '
 end subroutine calcHeight


 ! **********************************************************************************************************
 ! new subroutine: compute turbulent exchange coefficients
 ! **********************************************************************************************************
 subroutine turbExchng(err,message)
 ! used to compute derived model variables
 USE multiconst, only: vkc                     ! von Karman's constant
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out) :: err               ! error code
 character(*),intent(out) :: message           ! error message
 ! declare pointers to data in parameter structures
 real(dp),pointer         :: zon               ! roughness length (m)
 real(dp),pointer         :: mheight           ! measurement height (m)
 real(dp),pointer         :: bparam            ! parameter in Louis (1979) stability function
 real(dp),pointer         :: c_star            ! parameter in Louis (1979) stability function
 ! declare pointers to data in model variable structures
 real(dp),pointer         :: surfaceHeight     ! height of the layer interface (m)
 real(dp),pointer         :: ExNeut            ! exchange coefficient in neutral conditions
 real(dp),pointer         :: bprime            ! used in Louis (1979) stability function
 real(dp),pointer         :: cparam            ! used in Louis (1979) stability function
 ! initialize error control
 err=0; message='turbExchng/'
 ! assign local pointers to the values in the parameter structures
 zon            =>mpar_data%var(iLookPARAM%zon)                        ! roughness length (m)
 mheight        =>mpar_data%var(iLookPARAM%mheight)                    ! measurement height (m)
 bparam         =>mpar_data%var(iLookPARAM%bparam)                     ! parameter in Louis (1979) stability function
 c_star         =>mpar_data%var(iLookPARAM%c_star)                     ! parameter in Louis (1979) stability function
 ! assign local pointers to the values in the model variable structures
 surfaceHeight  =>mvar_data%var(iLookMVAR%iLayerHeight)%dat(0)         ! height of the upper-most layer interface (m)
 ExNeut         =>mvar_data%var(iLookMVAR%scalarExNeut)%dat(1)         ! exchange coefficient in neutral conditions
 bprime         =>mvar_data%var(iLookMVAR%scalarBprime)%dat(1)         ! used in Louis (1979) stability function
 cparam         =>mvar_data%var(iLookMVAR%scalarCparam)%dat(1)         ! used in Louis (1979) stability function
 ! check the measurement height is above the snow surface
 if(-surfaceHeight > mheight)then; err=20; message=trim(message)//'height of snow surface exceeds height of measurement'; return; endif
 ! ************************************************************************************************************************
 ! compute derived parameters for turbulent heat transfer -- note positive downwards, with 0 at soil surface, so snow height is negative
 ExNeut = (vkc**2._dp) / (log((mheight+surfaceHeight)/zon))**2._dp          ! exchange coefficient in neutral conditions
 bprime = bparam / 2._dp                                                    ! used in Louis (1979) stability function
 cparam = c_star * ExNeut * bparam * ((mheight+surfaceHeight)/zon)**0.5_dp  ! used in Louis (1979) stability function
 ! ************************************************************************************************************************
 end subroutine turbExchng

 ! **********************************************************************************************************
 ! new subroutine: compute vertical distribution of root density
 ! **********************************************************************************************************
 subroutine rootDensty(err,message)
 ! model decision structures
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only: &
 movingBoundary,             & ! moving lower boundary
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out) :: err                   ! error code
 character(*),intent(out) :: message               ! error message
 ! declare pointers to data in parameter structures
 real(dp),pointer         :: rootingDepth          ! rooting depth (m)
 real(dp),pointer         :: rootDistExp           ! exponent for the vertical distriution of root density (-)
 ! declare pointers to data in model variable structures
 real(dp),pointer         :: scalarAquiferRootFrac ! fraction of roots below the soil profile (-)
 real(dp),pointer         :: mLayerRootDensity(:)  ! fraction of roots in each soil layer (-)
 real(dp),pointer         :: iLayerHeight(:)       ! height of the layer interface (m)
 ! declare pointers to model index variables
 integer(i4b),pointer     :: nLayers               ! number of layers
 integer(i4b),pointer     :: layerType(:)          ! type of the layer (ix_soil or ix_snow)
 ! declare local variables
 integer(i4b)             :: nSoil,nSnow           ! number of soil and snow layers
 integer(i4b)             :: iLayer                ! loop through layers
 real(dp)                 :: fracRootLower         ! fraction of the rooting depth at the lower interface
 real(dp)                 :: fracRootUpper         ! fraction of the rooting depth at the upper interface
 ! initialize error control
 err=0; message='rootDensty/'

 ! assign local pointers to the values in the parameter structures
 rootingDepth          =>mpar_data%var(iLookPARAM%rootingDepth)                ! rooting depth (m)
 rootDistExp           =>mpar_data%var(iLookPARAM%rootDistExp)                 ! root distribution exponent (-)
 ! assign local pointers to the values in the model variable structures
 scalarAquiferRootFrac =>mvar_data%var(iLookMVAR%scalarAquiferRootFrac)%dat(1) ! fraction of roots below the soil profile (in the aquifer)
 mLayerRootDensity     =>mvar_data%var(iLookMVAR%mLayerRootDensity)%dat        ! fraction of roots in each soil layer (-)
 iLayerHeight          =>mvar_data%var(iLookMVAR%iLayerHeight)%dat             ! height of the layer interface (m)
 ! assign local pointers to the model index structures
 nLayers               =>indx_data%var(iLookINDEX%nLayers)%dat(1)              ! number of layers
 layerType             =>indx_data%var(iLookINDEX%layerType)%dat               ! layer type (ix_soil or ix_snow)
 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! for case of no explicit groundwater parameterization, check that the rooting depth is less than the soil depth
 if(model_decisions(iLookDECISIONS%groundwatr)%iDecision == noExplicit)then
  if(rootingDepth>iLayerHeight(nLayers))then; err=10; message=trim(message)//'rooting depth cannot exceed soil depth when there is no explicit gw representation'; return; endif
 endif

 ! compute the fraction of roots in each soil layer
 do iLayer=nSnow+1,nLayers
  if(iLayerHeight(iLayer-1)<rootingDepth)then
   ! compute the fraction of the rooting depth at the lower and upper interfaces
   fracRootLower = iLayerHeight(iLayer-1)/rootingDepth
   fracRootUpper = iLayerHeight(iLayer)/rootingDepth
   if(fracRootUpper>1._dp) fracRootUpper=1._dp
   ! compute the root density
   mLayerRootDensity(iLayer-nSnow) = fracRootUpper**rootDistExp - fracRootLower**rootDistExp
  else
   mLayerRootDensity(iLayer-nSnow) = 0._dp
  endif
  !print*, 'iLayerHeight(iLayer-1:iLayer) = ', iLayerHeight(iLayer-1:iLayer)
  !write(*,'(a,10(f11.5,1x))') 'mLayerRootDensity(iLayer-nSnow), fracRootUpper, fracRootLower, fracRootUpper**rootDistExp, fracRootLower**rootDistExp = ', &
  !                             mLayerRootDensity(iLayer-nSnow), fracRootUpper, fracRootLower, fracRootUpper**rootDistExp, fracRootLower**rootDistExp
 end do  ! (looping thru layers)
 ! check everything is OK
 if(rootingDepth < iLayerHeight(nLayers))then
  if(abs(sum(mLayerRootDensity) - 1._dp) > epsilon(rootingDepth))then
   message=trim(message)//'root density does not sum to one when rooting depth is within the soil profile'
   err=20; return
  endif
 endif

 end subroutine rootDensty

 ! **********************************************************************************************************
 ! new subroutine: compute vertical profile of saturated hydraulic conductivity
 ! **********************************************************************************************************
 subroutine satHydCond(err,message)
 ! model decision structures
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only: &
  constant,                  & ! constant hydraulic conductivity with depth
  exp_profile,               & ! exponential profile
  powerLaw_profile,          & ! power-law profile
  linear_profile               ! linear profile
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out) :: err                   ! error code
 character(*),intent(out) :: message               ! error message
 ! declare pointers to data in parameter structures
 real(dp),pointer         :: k_soil                ! saturated hydraulic conductivity at the compacted depth (m s-1)
 real(dp),pointer         :: compactedDepth        ! the depth at which k_soil reaches the compacted value given by CH78 (m)
 real(dp),pointer         :: zScale_TOPMODEL       ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! declare pointers to data in model variable structures
 real(dp),pointer         :: mLayerSatHydCond(:)   ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),pointer         :: iLayerSatHydCond(:)   ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 real(dp),pointer         :: mLayerHeight(:)       ! height at the mid-point of each layer (m)
 real(dp),pointer         :: iLayerHeight(:)       ! height at the interface of each layer (m)
 ! declare pointers to model index variables
 integer(i4b),pointer     :: nLayers               ! number of layers
 integer(i4b),pointer     :: layerType(:)          ! type of the layer (ix_soil or ix_snow)
 ! declare local variables
 integer(i4b)             :: nSoil,nSnow           ! number of soil and snow layers
 integer(i4b)             :: iLayer                ! loop through layers
 ! initialize error control
 err=0; message='satHydCond/'

 ! assign local pointers to the values in the parameter structures
 k_soil           => mpar_data%var(iLookPARAM%k_soil)              ! saturated hydraulic conductivity at the compacted depth (m s-1)
 compactedDepth   => mpar_data%var(iLookPARAM%compactedDepth)      ! the depth at which k_soil reaches the compacted value given by CH78 (m)
 zScale_TOPMODEL  => mpar_data%var(iLookPARAM%zScale_TOPMODEL)     ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! assign local pointers to the values in the model variable structures
 mLayerSatHydCond => mvar_data%var(iLookMVAR%mLayerSatHydCond)%dat ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 iLayerSatHydCond => mvar_data%var(iLookMVAR%iLayerSatHydCond)%dat ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 mLayerHeight     => mvar_data%var(iLookMVAR%mLayerHeight)%dat     ! height at the mid-point of each layer (m)
 iLayerHeight     => mvar_data%var(iLookMVAR%iLayerHeight)%dat     ! height at the interface of each layer (m)
 ! assign local pointers to the model index structures
 nLayers          => indx_data%var(iLookINDEX%nLayers)%dat(1)      ! number of layers
 layerType        => indx_data%var(iLookINDEX%layerType)%dat       ! layer type (ix_soil or ix_snow)
 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! loop through soil layers
 do iLayer=nSnow,nLayers
  select case(model_decisions(iLookDECISIONS%hc_profile)%iDecision)
   ! constant hydraulic conductivity with depth
   case(constant)
    iLayerSatHydCond(iLayer-nSnow) = k_soil
    if(iLayer > nSnow)& ! avoid layer 0
     mLayerSatHydCond(iLayer-nSnow) = k_soil
   ! exponential profile
   case(exp_profile)
    iLayerSatHydCond(iLayer-nSnow) = k_soil * exp(-(iLayerHeight(iLayer) - compactedDepth)/zScale_TOPMODEL)
    if(iLayer > nSnow)& ! avoid layer 0
     mLayerSatHydCond(iLayer-nSnow) = k_soil * exp(-(mLayerHeight(iLayer) - compactedDepth)/zScale_TOPMODEL)
   ! power-law and linear profile (not implemented yet)
   case(powerLaw_profile,linear_profile)
    message=trim(message)//"hydraulic conductivity profile not implemented yet [option="//trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)//"]"
    err=10; return
   ! error check (errors checked earlier also, so should not get here)
   case default
    message=trim(message)//"unknown hydraulic conductivity profile [option="//trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)//"]"
    err=10; return
  end select
  !if(iLayer > nSnow)& ! avoid layer 0
  ! write(*,'(i4,1x,2(f11.5,1x,e20.10,1x))') iLayer, mLayerHeight(iLayer), mLayerSatHydCond(iLayer-nSnow), iLayerHeight(iLayer), iLayerSatHydCond(iLayer-nSnow)
 end do  ! looping through soil layers

 end subroutine satHydCond


 ! **********************************************************************************************************
 ! new subroutine: compute "short-cut" variables
 ! **********************************************************************************************************
 subroutine v_shortcut(err,message)
 ! used to compute derived model variables
 USE multiconst, only:&
                       LH_fus,    &            ! latent heat of fusion                (J kg-1)
                       Cp_air,    &            ! specific heat of air                 (J kg-1 K-1)
                       Cp_ice,    &            ! specific heat of ice                 (J kg-1 K-1)
                       Cp_soil,   &            ! specific heat of soil                (J kg-1 K-1)
                       Cp_water,  &            ! specific heat of liquid water        (J kg-1 K-1)
                       iden_air,  &            ! intrinsic density of air             (kg m-3)
                       iden_ice,  &            ! intrinsic density of ice             (kg m-3)
                       iden_water,&            ! intrinsic density of liquid water    (kg m-3)
                       gravity,   &            ! gravitational acceleration           (m s-2)
                       Tfreeze                 ! freezing point of pure water         (K)
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out) :: err               ! error code
 character(*),intent(out) :: message           ! error message
 ! declare pointers to data in parameter structures
 real(dp),pointer         :: iden_soil         ! intrinsic density of soil (kg m-3)
 real(dp),pointer         :: frac_sand         ! fraction of sand (-)
 real(dp),pointer         :: frac_silt         ! fraction of silt (-)
 real(dp),pointer         :: frac_clay         ! fraction of clay (-)
 real(dp),pointer         :: theta_sat         ! soil porosity (-)
 real(dp),pointer         :: vGn_n             ! van Genutchen "n" parameter (-)
 real(dp),pointer         :: kappa             ! constant in the freezing curve function (m K-1)
 ! declare pointers to data in model variable structures
 real(dp),pointer         :: vGn_m             ! van Genutchen "m" parameter (-)
 real(dp),pointer         :: volHtCap_air      ! volumetric heat capacity of air (J m-3 K-1)
 real(dp),pointer         :: volHtCap_ice      ! volumetric heat capacity of ice (J m-3 K-1)
 real(dp),pointer         :: volHtCap_soil     ! volumetric heat capacity of soil (J m-3 K-1)
 real(dp),pointer         :: volHtCap_water    ! volumetric heat capacity of water (J m-3 K-1)
 real(dp),pointer         :: lambda_drysoil    ! thermal conductivity of dry soil (W m-1)
 real(dp),pointer         :: lambda_wetsoil    ! thermal conductivity of wet soil (W m-1)
 real(dp),pointer         :: volLatHt_fus      ! volumetric latent heat of fusion (J m-3)
 ! declare local variables
 real(dp)                 :: bulkden_soil      ! bulk density of soil (kg m-3)
 ! initialize error control
 err=0; message='v_shortcut/'
 ! assign local pointers to the values in the parameter structures
 iden_soil      =>mpar_data%var(iLookPARAM%soil_dens_intr)             ! intrinsic soil density (kg m-3)
 frac_sand      =>mpar_data%var(iLookPARAM%frac_sand)                  ! fraction of sand (-)
 frac_silt      =>mpar_data%var(iLookPARAM%frac_silt)                  ! fraction of silt (-)
 frac_clay      =>mpar_data%var(iLookPARAM%frac_clay)                  ! fraction of clay (-)
 theta_sat      =>mpar_data%var(iLookPARAM%theta_sat)                  ! soil porosity (-)
 vGn_n          =>mpar_data%var(iLookPARAM%vGn_n)                      ! van Genutchen "n" parameter (-)
 ! assign local pointers to the values in the model variable structures
 vGn_m          =>mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)          ! van Genutchen "m" parameter (-)
 kappa          =>mvar_data%var(iLookMVAR%scalarKappa)%dat(1)          ! constant in the freezing curve function (m K-1)
 volHtCap_air   =>mvar_data%var(iLookMVAR%scalarVolHtCap_air)%dat(1)   ! volumetric heat capacity of air (J m-3 K-1)
 volHtCap_ice   =>mvar_data%var(iLookMVAR%scalarVolHtCap_ice)%dat(1)   ! volumetric heat capacity of ice (J m-3 K-1)
 volHtCap_soil  =>mvar_data%var(iLookMVAR%scalarVolHtCap_soil)%dat(1)  ! volumetric heat capacity of soil (J m-3 K-1)
 volHtCap_water =>mvar_data%var(iLookMVAR%scalarVolHtCap_water)%dat(1) ! volumetric heat capacity of water (J m-3 K-1)
 lambda_drysoil =>mvar_data%var(iLookMVAR%scalarLambda_drysoil)%dat(1) ! thermal conductivity of dry soil (W m-1)
 lambda_wetsoil =>mvar_data%var(iLookMVAR%scalarLambda_wetsoil)%dat(1) ! thermal conductivity of wet soil (W m-1)
 volLatHt_fus   =>mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)
 ! ************************************************************************************************************************
 ! compute the van Genutchen "m" parameter
 vGn_m = 1._dp - 1._dp/vGn_n
 ! ************************************************************************************************************************
 ! compute the constant in the freezing curve function (m K-1)
 kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2
 ! ************************************************************************************************************************
 ! compute volumetric heat capacity (J m-3 K-1)
 volHtCap_air   = iden_air   * Cp_air
 volHtCap_ice   = iden_ice   * Cp_Ice
 volHtCap_soil  = iden_soil  * Cp_soil
 volHtCap_water = iden_water * Cp_water
 ! compute the thermal conductivity of dry and wet soils (W m-1)
 bulkden_soil   = iden_soil*(1._dp - theta_sat)
 lambda_drysoil = (0.135_dp*bulkden_soil + 64.7_dp) / (iden_soil - 0.947_dp*bulkden_soil)
 lambda_wetsoil = (8.80_dp*frac_sand + 2.92_dp*frac_clay) / (frac_sand + frac_clay)
 !print*, 'frac_sand, frac_silt, frac_clay = ', frac_sand, frac_silt, frac_clay
 !print*, 'lambda_drysoil, lambda_wetsoil = ', lambda_drysoil, lambda_wetsoil
 !print*, 'volHtCap_soil = ', volHtCap_soil
 ! compute the volumetric latent heat of fusion (J m-3)
 volLatHt_fus = iden_ice   * LH_fus
 ! ************************************************************************************************************************
 end subroutine v_shortcut

end module var_derive_module
