module flxMapping_module
implicit none
private
public::flxMapping
contains

 subroutine flxMapping(err,message)
 USE nrtype
 ! data types
 USE data_types, only: var_info        ! data type for metadata structure
 USE data_types, only: flux2state      ! data type for extended metadata structure, for flux-to-state mapping
 ! structures of named variables
 USE var_lookup, only: iLookFLUX       ! named variables for local flux variables
 ! metadata structures
 USE globalData, only: flux_meta       ! data structure for model fluxes
 USE globalData, only: flux2state_meta ! data structure for flux-to-state mapping
 ! named variables to describe the state variable type
 USE globalData,only:iname_nrgCanair   ! named variable defining the energy of the canopy air space
 USE globalData,only:iname_nrgCanopy   ! named variable defining the energy of the vegetation canopy
 USE globalData,only:iname_watCanopy   ! named variable defining the mass of water on the vegetation canopy
 USE globalData,only:iname_nrgLayer    ! named variable defining the energy state variable for snow+soil layers
 USE globalData,only:iname_watLayer    ! named variable defining the total water state variable for snow+soil layers
 USE globalData,only:iname_matLayer    ! named variable defining the matric head state variable for soil layers
 ! access missing values
 USE globalData,only:integerMissing    ! missing integer
 implicit none
 ! dummy variables
 integer(i4b),intent(out)       :: err                 ! error code
 character(*),intent(out)       :: message             ! error message
 ! local variables
 integer(i4b)                   :: iVar                ! variable index
 integer(i4b),parameter         :: integerUndefined=0  ! named variable to denote that the flux is undefined
 ! initialize error control
 err=0; message='flxMapping/'

 ! ** initialize flux-to-state mapping
 do iVar=1,size(flux_meta)
  flux2state_meta(iVar)%state1 = integerUndefined
  flux2state_meta(iVar)%state2 = integerUndefined
 end do

 ! ** define mapping between fluxes and states

 ! net energy and mass fluxes for the vegetation domain
 flux2state_meta(iLookFLUX%scalarCanopyNetLiqFlux)         = flux2state(state1=iname_watCanopy, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarCanairNetNrgFlux)         = flux2state(state1=iname_nrgCanair, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarCanopyNetNrgFlux)         = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer) 
 flux2state_meta(iLookFLUX%scalarGroundNetNrgFlux)         = flux2state(state1=iname_nrgLayer,  state2=iname_nrgCanopy)

 ! precipitation -- does not depend on state variables
 flux2state_meta(iLookFLUX%scalarRainfall)                  = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarSnowfall)                  = flux2state(state1=integerMissing, state2=integerMissing)
 
 ! shortwave radiation -- does not depend on state variables
 flux2state_meta(iLookFLUX%spectralIncomingDirect)          = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%spectralIncomingDiffuse)         = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarCanopySunlitPAR)           = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarCanopyShadedPAR)           = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%spectralBelowCanopyDirect)       = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%spectralBelowCanopyDiffuse)      = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarBelowCanopySolar)          = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarCanopyAbsorbedSolar)       = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarGroundAbsorbedSolar)       = flux2state(state1=integerMissing, state2=integerMissing)
 
 ! longwave radiation -- assume calculated when the canopy energy state variable is active OR when the ground energy state variable is active
 flux2state_meta(iLookFLUX%scalarLWRadCanopy)               = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadGround)               = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadUbound2Canopy)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadUbound2Ground)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadUbound2Ubound)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadCanopy2Ubound)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadCanopy2Ground)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadCanopy2Canopy)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadGround2Ubound)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWRadGround2Canopy)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWNetCanopy)               = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWNetGround)               = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLWNetUbound)               = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 
 ! turbulent heat transfer -- assume calculated when the canopy energy state variable is active OR when the ground energy state variable is active 
 flux2state_meta(iLookFLUX%scalarEddyDiffusCanopyTop)       = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarFrictionVelocity)          = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarWindspdCanopyTop)          = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarWindspdCanopyBottom)       = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarGroundResistance)          = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarCanopyResistance)          = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLeafResistance)            = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarSoilResistance)            = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarSenHeatTotal)              = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarSenHeatCanopy)             = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarSenHeatGround)             = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLatHeatTotal)              = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLatHeatCanopyEvap)         = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLatHeatCanopyTrans)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarLatHeatGround)             = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarCanopyAdvectiveHeatFlux)   = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarGroundAdvectiveHeatFlux)   = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarCanopySublimation)         = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarSnowSublimation)           = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 
 ! stomatal resistance and photosynthesis -- calculated when the canopy energy state variable is active
 flux2state_meta(iLookFLUX%scalarStomResistSunlit)          = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarStomResistShaded)          = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarPhotosynthesisSunlit)      = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarPhotosynthesisShaded)      = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)

 ! liquid water fluxes associated with evapotranspiration
 ! NOTE 1: calculated in the energy balance routines: energy balance must be calculated first in order for water to balance
 ! NOTE 2: if implement strang splitting, need to average fluxes from the start and end of the time step
 flux2state_meta(iLookFLUX%scalarCanopyTranspiration)       = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarCanopyEvaporation)         = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%scalarGroundEvaporation)         = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_meta(iLookFLUX%mLayerTranspire)                 = flux2state(state1=iname_matLayer,  state2=integerMissing)

 ! liquid and solid water fluxes through the canopy
 flux2state_meta(iLookFLUX%scalarThroughfallSnow)           = flux2state(state1=integerMissing,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarCanopySnowUnloading)       = flux2state(state1=integerMissing,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarThroughfallRain)           = flux2state(state1=iname_watCanopy, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarCanopyLiqDrainage)         = flux2state(state1=iname_watCanopy, state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarCanopyMeltFreeze)          = flux2state(state1=integerMissing,  state2=integerMissing)
 
 ! energy fluxes and for the snow and soil domains
 flux2state_meta(iLookFLUX%iLayerConductiveFlux)            = flux2state(state1=iname_nrgLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%iLayerAdvectiveFlux)             = flux2state(state1=iname_nrgLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%iLayerNrgFlux)                   = flux2state(state1=iname_nrgLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerNrgFlux)                   = flux2state(state1=iname_nrgLayer,  state2=integerMissing)
 
 ! liquid water fluxes for the snow domain
 flux2state_meta(iLookFLUX%iLayerLiqFluxSnow)               = flux2state(state1=iname_watLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerLiqFluxSnow)               = flux2state(state1=iname_watLayer,  state2=integerMissing)
 
 ! liquid water fluxes for the soil domain
 flux2state_meta(iLookFLUX%scalarRainPlusMelt)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarMaxInfilRate)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarInfiltration)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarExfiltration)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarSurfaceRunoff)             = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerSatHydCondMP)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerSatHydCond)                = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%iLayerSatHydCond)                = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerHydCond)                   = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%iLayerLiqFluxSoil)               = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerLiqFluxSoil)               = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerBaseflow)                  = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerColumnInflow)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%mLayerColumnOutflow)             = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarSoilBaseflow)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarSoilDrainage)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarAquiferRecharge)           = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarAquiferTranspire)          = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_meta(iLookFLUX%scalarAquiferBaseflow)           = flux2state(state1=iname_matLayer,  state2=integerMissing)

 ! ** copy across flux metadata
 do iVar=1,size(flux_meta)
  flux2state_meta(iVar)%var_info = flux_meta(iVar)
 end do

 ! ** check all variables are defined
 do iVar=1,size(flux_meta)
  if(flux2state_meta(iVar)%state1==integerUndefined .or. flux2state_meta(iVar)%state2==integerUndefined)then
   message=trim(message)//'flux-to-state mapping is undefined for variable "'//trim(flux_meta(iVar)%varname)//'"'
   err=20; return
  endif
 end do

 end subroutine flxMapping

end module flxMapping_module
