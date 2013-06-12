module newsnwfall_module
USE nrtype
implicit none
private
public::newsnwfall

contains

 ! ************************************************************************************************
 ! new subroutine: add new snowfall to the system, and increase number of snow layers if needed
 ! ************************************************************************************************
 subroutine newsnwfall(dt,             & ! time step (seconds)
                       err,message)      ! error control
 ! physical constants
 USE multiconst,only:&
                     Tfreeze,  & ! freezing point              (K)
                     LH_fus,   & ! latent heat of fusion       (J kg-1)
                     LH_vap,   & ! latent heat of vaporization (J kg-1)
                     LH_sub,   & ! latent heat of sublimation  (J kg-1)
                     iden_air, & ! intrinsic density of air    (kg m-3)
                     iden_ice, & ! intrinsic density of ice    (kg m-3)
                     iden_water  ! intrinsic density of water  (kg m-3)
 ! computational modules
 USE snow_utils_module,only:fracliquid,templiquid                  ! functions to compute temperature/liquid water
 ! data structures
 USE data_struc,only:mpar_meta,forc_meta,mvar_meta,indx_meta       ! metadata
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data       ! data structures
 USE data_struc,only:ix_soil,ix_snow,ix_mixd                       ! names of model layers
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX    ! named variables for structure elements
 ! add new snowfall to the system
 implicit none
 ! dummy variables
 real(dp),intent(in)                 :: dt                         ! time step (seconds)
 integer(i4b),intent(out)            :: err                        ! error code
 character(*),intent(out)            :: message                    ! error message
 ! assign pointers to model parameters
 real(dp),pointer                    :: fc_param                   ! freeezing curve parameter for snow (K-1)
 ! local pointers to model state variables
 real(dp),pointer                    :: surfaceLayerTemp           ! temperature of each layer (K)
 real(dp),pointer                    :: surfaceLayerDepth          ! depth of each layer (m)
 real(dp),pointer                    :: surfaceLayerVolFracIce     ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                    :: surfaceLayerVolFracLiq     ! volumetric fraction of liquid water in each layer (-)
 ! local pointers to diagnostic scalar variables
 real(dp),pointer                    :: scalarSWE                  ! SWE (kg m-2)
 real(dp),pointer                    :: scalarSnowDepth            ! total snow depth (m)
 real(dp),pointer                    :: scalarSnowfallTemp         ! computed temperature of fresh snow (K) 
 real(dp),pointer                    :: scalarNewSnowDensity       ! computed density of new snow (kg m-3)
 real(dp),pointer                    :: scalarThroughfallSnow      ! throughfall of snow through the canopy (kg m-2 s-1)
 real(dp),pointer                    :: scalarCanopySnowUnloading  ! unloading of snow from the canopy (kg m-2 s-1)
 ! local pointers to model index variables
 integer(i4b),pointer                :: nLayers                    ! number of layers
 integer(i4b),pointer                :: layerType(:)               ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 integer(i4b)                        :: nSoil                      ! number of soil layers
 integer(i4b)                        :: nSnow                      ! number of snow layers
 real(dp)                            :: newSnowfall                ! new snowfall -- throughfall and unloading (kg m-2 s-1)
 real(dp)                            :: newSnowDepth               ! new snow depth (m)
 real(dp),parameter                  :: densityCanopySnow=200._dp  ! density of snow on the vegetation canopy (kg m-3)
 real(dp)                            :: totalMassIceSurfLayer      ! total mass of ice in the surface layer (kg m-2)
 real(dp)                            :: totalDepthSurfLayer        ! total depth of the surface layer (m)
 real(dp)                            :: volFracWater               ! volumetric fraction of total water, liquid and ice (-)
 real(dp)                            :: fracLiq                    ! fraction of liquid water (-)
 ! initialize error control
 err=0; message="newsnwfall/"
 ! assign pointers to model parameters (compute layer temperature)
 fc_param                  => mpar_data%var(iLookPARAM%snowfrz_scale)                    ! freezing curve parameter for snow (K-1)
 ! assign local pointers to diagnostic scalar variables
 scalarSWE                 => mvar_data%var(iLookMVAR%scalarSWE)%dat(1)                  ! SWE (kg m-2)
 scalarSnowDepth           => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)            ! total snow depth (m)
 scalarSnowfallTemp        => mvar_data%var(iLookMVAR%scalarSnowfallTemp)%dat(1)         ! computed temperature of fresh snow (K)
 scalarNewSnowDensity      => mvar_data%var(iLookMVAR%scalarNewSnowDensity)%dat(1)       ! computed density of new snow (kg m-3) 
 scalarThroughfallSnow     => mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1)      ! throughfall of snow through the canopy (kg m-2 s-1)
 scalarCanopySnowUnloading => mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1)  ! unloading of snow from the canopy (kg m-2 s-1)

 ! assign local pointers to the model index structures
 nLayers                => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType              => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

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
  ! assign pointers to model state variables -- all surface layer only
  surfaceLayerTemp       => mvar_data%var(iLookMVAR%mLayerTemp)%dat(1)        ! temperature of the surface layer (K)
  surfaceLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat(1)       ! depth of the surface layer (m)
  surfaceLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1)  ! volumetric fraction of ice in the surface layer  (-)
  surfaceLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1)  ! volumetric fraction of liquid water in the surface layer (-)
  ! get the total mass of liquid water and ice (kg m-2)
  totalMassIceSurfLayer  = iden_ice*surfaceLayerVolFracIce*surfaceLayerDepth + newSnowfall*dt
  ! get the total snow depth
  totalDepthSurfLayer    = surfaceLayerDepth + newSnowDepth
  ! compute the new temperature
  surfaceLayerTemp       = (surfaceLayerTemp*surfaceLayerDepth + scalarSnowfallTemp*newSnowDepth) / totalDepthSurfLayer
  ! compute new volumetric fraction of liquid water and ice
  volFracWater = (totalMassIceSurfLayer/totalDepthSurfLayer)/iden_water + surfaceLayerVolFracLiq  ! volumetric fraction of total water (liquid and ice)
  fracLiq      = fracliquid(surfaceLayerTemp,fc_param)                           ! fraction of liquid water
  surfaceLayerVolFracIce = (1._dp - fracLiq)*volFracWater*(iden_water/iden_ice)  ! volumetric fraction of ice (-)
  surfaceLayerVolFracLiq =          fracLiq *volFracWater                        ! volumetric fraction of liquid water (-)
  ! update new layer depth (m)
  surfaceLayerDepth      = totalDepthSurfLayer
  ! re-compute snow depth and SWE
  scalarSnowDepth = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
  scalarSWE       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                          mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                          * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
 endif

 end subroutine newsnwfall

end module newsnwfall_module
