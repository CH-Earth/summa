module newsnwfall_module
USE nrtype
implicit none
private
public::newsnwfall
interface addOneLayer
 module procedure AddOneLayer_rv, AddOneLayer_iv
end interface AddOneLayer

contains

 ! ************************************************************************************************
 ! new subroutine: add new snowfall to the system, and increase number of snow layers if needed
 ! ************************************************************************************************
 subroutine newsnwfall(dt,            & ! time step (seconds)
                       err,message)     ! error control
 USE multiconst,only:&
                     Tfreeze,  & ! freezing point              (K)
                     LH_fus,   & ! latent heat of fusion       (J kg-1)
                     LH_vap,   & ! latent heat of vaporization (J kg-1)
                     LH_sub,   & ! latent heat of sublimation  (J kg-1)
                     iden_air, & ! intrinsic density of air    (kg m-3)
                     iden_ice, & ! intrinsic density of ice    (kg m-3)
                     iden_water  ! intrinsic density of water  (kg m-3)
 USE data_struc,only:mpar_meta,forc_meta,mvar_meta,indx_meta    ! metadata
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE data_struc,only:ix_soil,ix_snow,ix_mixd                    ! names of model layers
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 USE snow_utils_module,only:fracliquid,templiquid               ! functions to compute temperature/liquid water
 ! add new snowfall to the system
 implicit none
 ! dummy variables
 real(dp),intent(in)                 :: dt                       ! time step (seconds)
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! local pointers to model parameters (new snow density)
 real(dp),pointer                    :: newSnowDenMin            ! minimum new snow density (kg m-3)
 real(dp),pointer                    :: newSnowDenMult           ! multiplier for new snow density (kg m-3)
 real(dp),pointer                    :: newSnowDenScal           ! scaling factor for new snow density (K)
 ! local pointers to model parameters (control on the depth of snow layers)
 real(dp),pointer                    :: zmin                     ! minimum layer depth (m)
 real(dp),pointer                    :: zmax                     ! maximum layer depth (m)
 real(dp),pointer                    :: fc_param                 ! freeezing curve parameter for snow (K-1)
 ! local pointers to model forcing data
 real(dp),pointer                    :: airtemp                  ! air temperature at 2 meter height (K)
 ! local pointers to model state variables (surface layer only)
 real(dp),pointer                     :: surfaceLayerTemp        ! temperature of each layer (K)
 real(dp),pointer                     :: surfaceLayerDepth       ! depth of each layer (m)
 real(dp),pointer                     :: surfaceLayerVolFracIce  ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                     :: surfaceLayerVolFracLiq  ! volumetric fraction of liquid water in each layer (-)
 ! local pointers to model coordinate variables
 real(dp),pointer                     :: mLayerDepth(:)          ! depth of the layer (m)
 real(dp),pointer                     :: mLayerHeight(:)         ! height of the layer mid-point (m)
 real(dp),pointer                     :: iLayerHeight(:)         ! height of the layer interface (m)
 ! local pointers to diagnostic scalar variables
 real(dp),pointer                     :: scalarSnowDepth         ! total snow depth (m)
 real(dp),pointer                     :: scalarSWE               ! SWE (kg m-2)
 real(dp),pointer                     :: scalarSnowfall          ! snowfall flux (kg m-2 s-1)
 real(dp),pointer                     :: scalarSnowfallTemp      ! computed temperature of fresh snow (K) 
 ! local pointers to model index variables
 integer(i4b),pointer                :: nLayers                  ! number of layers
 integer(i4b),pointer                :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)                  :: cmessage                 ! error message of downwind routine
 integer(i4b)                        :: iLayer                   ! layer index
 integer(i4b)                        :: ix_lower                 ! lower bound of the vector
 integer(i4b)                        :: ix_upper                 ! upper bound of the vector
 integer(i4b)                        :: ivar                     ! variable index
 integer(i4b)                        :: nSoil                    ! number of soil layers
 integer(i4b)                        :: nSnow                    ! number of snow layers
 real(dp)                            :: newSnowDensity           ! new snow density (kg m-3)
 real(dp)                            :: newSnowDepth             ! newSnowDepth (m)
 real(dp)                            :: surfaceLayerSoilTemp     ! temperature of the top soil layer (K)
 real(dp)                            :: maxFrozenSnowTemp        ! maximum temperature when effectively all water is frozen (K)
 real(dp),parameter                  :: unfrozenLiq=0.01_dp      ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(dp)                            :: totalMassIceSurfLayer    ! total mass of ice in the surface layer (kg m-2)
 real(dp)                            :: totalDepthSurfLayer      ! total depth of the surface layer (m)
 real(dp)                            :: volFracWater             ! volumetric fraction of total water, liquid and ice (-)
 real(dp)                            :: fracLiq                  ! fraction of liquid water (-)
 real(dp),parameter                  :: fracTop=0.333333333_dp   ! fraction of old layer used for the top layer
 ! initialize error control
 err=0; message="newsnwfall/"
 ! assign pointers to model parameters (new snow density)
 newSnowDenMin          => mpar_data%var(iLookPARAM%newSnowDenMin)            ! minimum new snow density (kg m-3)
 newSnowDenMult         => mpar_data%var(iLookPARAM%newSnowDenMult)           ! multiplier for new snow density (kg m-3)
 newSnowDenScal         => mpar_data%var(iLookPARAM%newSnowDenScal)           ! scaling factor for new snow density (K)
 ! assign pointers to model parameters (control the depth of snow layers)
 zmin                   => mpar_data%var(iLookPARAM%zmin)                     ! minimum layer depth (m)
 zmax                   => mpar_data%var(iLookPARAM%zmax)                     ! maximum layer depth (m)
 fc_param               => mpar_data%var(iLookPARAM%snowfrz_scale)            ! freezing curve parameter for snow (K-1)
 ! assign pointers to model forcing variables
 airtemp                => forc_data%var(iLookFORCE%airtemp)                  ! air temperature at 2 meter height (K)
 ! assign local pointers to diagnostic scalar variables
 scalarSnowfall         => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1)     ! snowfall flux (kg m-2 s-1)
 scalarSnowfallTemp     => mvar_data%var(iLookMVAR%scalarSnowfallTemp)%dat(1) ! computed temperature of fresh snow (K)
 scalarSnowDepth        => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)    ! total snow depth (m)
 scalarSWE              => mvar_data%var(iLookMVAR%scalarSWE)%dat(1)          ! SWE (kg m-2)

 ! assign local pointers to the model index structures
 nLayers                => indx_data%var(iLookINDEX%nLayers)%dat(1)           ! number of layers
 layerType              => indx_data%var(iLookINDEX%layerType)%dat            ! layer type (ix_soil or ix_snow)

 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! compute snow density (kg m-3)
 newSnowDensity = newSnowDenMin + newSnowDenMult*exp((airtemp-Tfreeze)/newSnowDenScal)

 ! compute snow depth (m)
 newSnowDepth = dt*scalarSnowfall/newSnowDensity

 ! ***** special case of no snow layers
 if(nSnow==0)then
  ! increment depth and water equivalent
  scalarSnowDepth = scalarSnowDepth + newSnowDepth
  scalarSWE       = scalarSWE + dt*scalarSnowfall
  ! check if create the first snow layer
  if(scalarSnowDepth > 0.5_dp*(zmin + zmax))then
   ! add a layer to all model variables
   call addModelLayer(err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
   ! define index of the top (new) layer as ix_snow
   indx_data%var(iLookINDEX%layerType)%dat(1) = ix_snow
   ! define model coordinate variables -- use vectors directly, since just altered
   mvar_data%var(iLookMVAR%mLayerDepth)%dat(1)  = scalarSnowDepth               ! depth of the layer (m)
   mvar_data%var(iLookMVAR%mLayerHeight)%dat(1) = -scalarSnowDepth/2._dp        ! height of the layer mid-point (m)
   mvar_data%var(iLookMVAR%iLayerHeight)%dat(0) = -scalarSnowDepth              ! height of the layer interface (m)
   mvar_data%var(iLookMVAR%iLayerHeight)%dat(1) = 0._dp                         ! height of the layer interface (m)
   ! ***** define model state variables -- use vectors directly, since just altered
   ! compute surface layer temperature
   surfaceLayerSoilTemp = mvar_data%var(iLookMVAR%mLayerTemp)%dat(2)    ! temperature of the top soil layer (K)
   maxFrozenSnowTemp    = templiquid(unfrozenLiq,fc_param)              ! snow temperature at fraction "unfrozenLiq" (K)
   mvar_data%var(iLookMVAR%mLayerTemp)%dat(1)       = min(maxFrozenSnowTemp,surfaceLayerSoilTemp)    ! snow temperature  (K)
   ! compute the fraction of liquid water associated with the layer temperature
   fracLiq      = fracliquid(mvar_data%var(iLookMVAR%mLayerTemp)%dat(1),fc_param)
   ! compute volumeteric fraction of liquid water and ice
   volFracWater = (scalarSWE/scalarSnowDepth)/iden_ice  ! volumetric fraction of total water (liquid and ice)
   mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1) = (1._dp - fracLiq)*volFracWater*(iden_water/iden_ice)   ! volumetric fraction of ice (-)
   mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1) =          fracLiq *volFracWater                         ! volumetric fraction of liquid water (-)
   print*, 'created a new layer, nSnow = ', count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
  endif  ! if creating a new layer
  return
 endif

 ! ***** merge snowfall into the surface layer
 ! NOTE: may need to re-visit this based on enthalpy considerations...
 ! assign pointers to model state variables -- all surface layer only
 surfaceLayerTemp       => mvar_data%var(iLookMVAR%mLayerTemp)%dat(1)        ! temperature of the surface layer (K)
 surfaceLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat(1)       ! depth of the surface layer (m)
 surfaceLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1)  ! volumetric fraction of ice in the surface layer  (-)
 surfaceLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1)  ! volumetric fraction of liquid water in the surface layer (-)
 ! get the total mass of ice (kg m-2)
 totalMassIceSurfLayer  = iden_ice*surfaceLayerVolFracIce*surfaceLayerDepth + scalarSnowfall*dt
 ! get the total snow depth
 totalDepthSurfLayer    = surfaceLayerDepth + newSnowDepth
 ! compute the new temperature
 surfaceLayerTemp       = (surfaceLayerTemp*surfaceLayerDepth + scalarSnowfallTemp*newSnowDepth) / totalDepthSurfLayer
 ! compute new volumetric fraction of liquid water and ice
 volFracWater = (totalMassIceSurfLayer/totalDepthSurfLayer)/iden_water + surfaceLayerVolFracLiq  ! volumetric fraction of total water (liquid and ice)
 fracLiq      = fracliquid(surfaceLayerTemp,fc_param)                           ! fraction of liquid water
 surfaceLayerVolFracIce = (1._dp - fracLiq)*volFracWater*(iden_water/iden_ice)  ! volumetric fraction of ice (-)
 surfaceLayerVolFracLiq =          fracLiq *volFracWater                        ! volumetric fraction of liquid water (-)
 ! compute new layer depth (m)
 surfaceLayerDepth      = totalDepthSurfLayer

 ! ***** sub-divide the layer, if necessary
 if(surfaceLayerDepth > zmax)then
  ! add a layer to all model variables
  call addModelLayer(err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  ! assign pointers to model coordinate variables
  mLayerDepth  => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of the layer (m)
  mLayerHeight => mvar_data%var(iLookMVAR%mLayerHeight)%dat         ! height of the layer mid-point (m)
  iLayerHeight => mvar_data%var(iLookMVAR%iLayerHeight)%dat         ! height of the layer interface (m)
  ! modify the depth of each snow layer
  mLayerDepth(1) = fracTop*mLayerDepth(2)
  mLayerDepth(2) = (1._dp-fracTop)*mLayerDepth(2)
  ! loop over the top two layers, and modify the height at layer interfaces and at layer mid-points
  do iLayer=2,1,-1
   iLayerHeight(iLayer-1) = iLayerHeight(iLayer) - mLayerDepth(iLayer)
   mLayerHeight(iLayer)   = iLayerHeight(iLayer) - mLayerDepth(iLayer)*0.5_dp
  end do
  ! copy snow temperature and volumetric ice and liquid water content
  mvar_data%var(iLookMVAR%mLayerTemp)%dat(1)       = mvar_data%var(iLookMVAR%mLayerTemp)%dat(2)
  mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1) = mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(2)
  mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1) = mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(2)
  ! label the layers as snow
  indx_data%var(iLookINDEX%layerType)%dat(1:2) = ix_snow
  print*, 'created a new layer, nSnow = ', count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
 endif  ! (if sub-dividing layer)

 contains


  ! *********************************************************************************************
  ! new subroutine: add an additional layer to all model vectors
  ! *********************************************************************************************
  subroutine addModelLayer(err,message)
  implicit none
  ! dummy variables
  integer(i4b),intent(out)            :: err                      ! error code
  character(*),intent(out)            :: message                  ! error message
  ! initialize error control
  err=0; message='addModelLayer/'

  ! ***** add a layer to each model variable
  do ivar=1,size(mvar_data%var)
   ! define bounds
   select case(trim(mvar_meta(ivar)%vartype))
    case('midSnow'); ix_lower=1; ix_upper=nSnow
    case('midToto'); ix_lower=1; ix_upper=nLayers
    case('ifcSnow'); ix_lower=0; ix_upper=nSnow
    case('ifcToto'); ix_lower=0; ix_upper=nLayers
    case default; cycle
   end select
   ! add an additional layer -- only get to here if snow in the layer
   call AddOneLayer(mvar_data%var(ivar)%dat,ix_lower,ix_upper,err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  end do
  ! ***** modify the layer indices
  call AddOneLayer(indx_data%var(iLookINDEX%layerType)%dat,1,nLayers,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  indx_data%var(iLookINDEX%nLayers)%dat(1)  = indx_data%var(iLookINDEX%nLayers)%dat(1) + 1
  end subroutine addModelLayer


 end subroutine newsnwfall


 ! ************************************************************************************************
 ! new subroutine: add an additional snow layer
 ! ************************************************************************************************
 subroutine AddOneLayer_rv(datavec,ix_lower,ix_upper,err,message)
 ! Returns a new vector which has one more element than the input vector
 !  -- optionally copies data from the original vector to the new vector for elements (2:n)->(3:n+1),
 !      and copies element 1 into elements 1:2, and copies element 0 into element 0
 implicit none
 ! dummies
 real(dp),pointer,intent(inout)     :: datavec(:)  ! the original and the new vector
 integer(i4b),intent(in)            :: ix_lower    ! lower bound of the old vector
 integer(i4b),intent(in)            :: ix_upper    ! upper bound of the old vector
 integer(i4b),intent(out)           :: err         ! error code
 character(*),intent(out)           :: message     ! error message
 ! locals
 real(dp)                           :: tempvec(ix_lower:ix_upper)  ! temporary vector
 real(dp),parameter                 :: missingReal=-9999._dp
 ! initialize error control
 err=0; message='AddOneLayer_rv/'
 ! check the data vector is associated
 if(.not.associated(datavec))then; err=20; message='data vector is not associated'; return; endif 
 ! assign the data vector to the temporary vector
 tempvec=datavec
 ! reallocate space for the new vector
 deallocate(datavec,stat=err)
 if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; endif
 allocate(datavec(ix_lower:ix_upper+1),stat=err)
 if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; endif
 ! only copy data if the vector exists (can be a variable for snow, with no layers)
 if(ix_upper>0) datavec(2:ix_upper+1) = tempvec(1:ix_upper)
 ! set undefined elements to missing
 datavec(ix_lower:1) = missingReal
 end subroutine AddOneLayer_rv

 subroutine AddOneLayer_iv(datavec,ix_lower,ix_upper,err,message)
 ! Returns a new vector which has one more element than the input vector
 !  -- optionally copies data from the original vector to the new vector for elements (2:n)->(3:n+1),
 !      and copies element 1 into elements 1:2, and copies element 0 into element 0
 implicit none
 ! dummies
 integer(i4b),pointer,intent(inout) :: datavec(:)  ! the original and the new vector
 integer(i4b),intent(in)            :: ix_lower    ! lower bound of the old vector
 integer(i4b),intent(in)            :: ix_upper    ! upper bound of the old vector
 integer(i4b),intent(out)           :: err         ! error code
 character(*),intent(out)           :: message     ! error message
 ! locals
 integer(i4b)                       :: tempvec(ix_lower:ix_upper)  ! temporary vector
 integer(i4b),parameter             :: missingInteger=-9999
 ! initialize error control
 err=0; message='AddOneLayer_iv/'
 ! check the data vector is associated
 if(.not.associated(datavec))then; err=20; message='data vector is not associated'; return; endif
 ! assign the data vector to the temporary vector
 tempvec=datavec
 ! reallocate space for the new vector
 deallocate(datavec,stat=err)
 if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; endif
 allocate(datavec(ix_lower:ix_upper+1),stat=err)
 if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; endif
 ! only copy data if the vector exists (can be a variable for snow, with no layers)
 if(ix_upper>0) datavec(2:ix_upper+1) = tempvec(1:ix_upper)  ! copy elements to lower layers
 ! set undefined elements to missing
 datavec(ix_lower:1) = missingInteger
 end subroutine AddOneLayer_iv

end module newsnwfall_module
