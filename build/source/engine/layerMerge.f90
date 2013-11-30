module layerMerge_module
USE nrtype
! define look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
implicit none
private
public::layerMerge
interface removeOneLayer
 module procedure removeOneLayer_rv, removeOneLayer_iv
end interface removeOneLayer

contains

 ! ************************************************************************************************
 ! new subroutine: merge layers if the thickness is less than zmin
 ! ************************************************************************************************
 subroutine layerMerge(err,message)    ! error control
 ! physical constants
 USE multiconst,only:&
                     iden_ice,       & ! intrinsic density of ice             (kg m-3)
                     iden_water        ! intrinsic density of liquid water    (kg m-3)
 ! compute modules
 USE var_derive_module,only:calcHeight ! module to calculate height at layer interfaces and layer mid-point
 ! model decisions
 USE data_struc,only:model_decisions                            ! model decision structure
 USE var_lookup,only:iLookDECISIONS                             ! named variables for elements of the decision structure
 ! data structures
 USE data_struc,only:mpar_meta,forc_meta,mvar_meta,indx_meta                    ! metadata
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 ! dummy variables
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! local pointers to model parameters (control on the depth of snow layers)
 real(dp),pointer                    :: zmin                     ! minimum layer depth (m)
 real(dp),pointer                    :: zminLayer1               ! minimum layer depth for the 1st (top) layer(m) 
 real(dp),pointer                    :: zminLayer2               ! minimum layer depth for the 2nd layer(m) 
 real(dp),pointer                    :: zminLayer3               ! minimum layer depth for the 3rd layer(m) 
 real(dp),pointer                    :: zminLayer4               ! minimum layer depth for the 4th layer(m) 
 real(dp),pointer                    :: zminLayer5               ! minimum layer depth for the 5th (bottom) layer(m) 
 ! local pointers to model state variables
 real(dp),pointer                    :: mLayerDepth(:)           ! depth of each layer (m)
 real(dp),pointer                    :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                    :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 ! local pointers to diagnostic scalar variables
 real(dp),pointer                     :: scalarSnowDepth         ! total snow depth (m)
 real(dp),pointer                     :: scalarSWE               ! SWE (kg m-2)
 ! local pointers to model index variables
 integer(i4b),pointer                :: nLayers                  ! number of layers
 integer(i4b),pointer                :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)                  :: cmessage                 ! error message of downwind routine
 real(dp),dimension(5)               :: zminLayer                ! minimum layer depth in each layer (m)
 logical(lgt)                        :: removeLayer              ! flag to indicate need to remove a layer
 integer(i4b)                        :: nCheck                   ! number of layers to check for combination
 integer(i4b)                        :: iSnow                    ! index of snow layers (looping)
 integer(i4b)                        :: jSnow                    ! index of snow layer identified for combination with iSnow
 integer(i4b)                        :: kSnow                    ! index of the upper layer of the two layers identified for combination
 integer(i4b)                        :: nSnow                    ! number of snow layers
 integer(i4b)                        :: nSoil                    ! number of soil layers
 ! initialize error control
 err=0; message="layerMerge/"
 ! assign local pointers to model parameters (control the depth of snow layers)
 zmin             => mpar_data%var(iLookPARAM%zmin)              ! minimum layer depth (m)
 zminLayer1       => mpar_data%var(iLookPARAM%zminLayer1)        ! minimum layer depth for the 1st (top) layer (m)
 zminLayer2       => mpar_data%var(iLookPARAM%zminLayer2)        ! minimum layer depth for the 2nd layer (m)
 zminLayer3       => mpar_data%var(iLookPARAM%zminLayer3)        ! minimum layer depth for the 3rd layer (m)
 zminLayer4       => mpar_data%var(iLookPARAM%zminLayer4)        ! minimum layer depth for the 4th layer (m)
 zminLayer5       => mpar_data%var(iLookPARAM%zminLayer5)        ! minimum layer depth for the 5th (bottom) layer (m)

 ! assign local pointers to model state variables
 mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
 mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat     ! volumetric fraction of ice in each layer  (-)
 mLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat     ! volumetric fraction of liquid water in each layer (-)
 ! assign local pointers to diagnostic scalar variables
 scalarSnowDepth  => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)   ! total snow depth (m)
 scalarSWE        => mvar_data%var(iLookMVAR%scalarSWE)%dat(1)         ! SWE (kg m-2)
 ! assign local pointers to the model index structures
 nLayers          => indx_data%var(iLookINDEX%nLayers)%dat(1)          ! number of layers
 layerType        => indx_data%var(iLookINDEX%layerType)%dat           ! layer type (ix_soil or ix_snow)

 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! identify algorithmic control parameters to syb-divide and combine snow layers
 zminLayer = (/zminLayer1, zminLayer2, zminLayer3, zminLayer4, zminLayer5/)

 kSnow=0 ! initialize first layer to test (top layer)
 do ! attempt to remove multiple layers in a single time step (continuous do loop with exit clause)

  ! special case of >5 layers: add an offset to use maximum threshold from layer above
  if(model_decisions(iLookDECISIONS%snowLayers)%iDecision == rulesDependLayerIndex .and. nSnow > 5)then
   nCheck=5
  else
   nCheck=nSnow
  endif

  ! loop through snow layers
  do iSnow=kSnow+1,nCheck

   ! check if the layer depth is less than the depth threshold
   select case(model_decisions(iLookDECISIONS%snowLayers)%iDecision)
    case(sameRulesAllLayers);    removeLayer = (mLayerDepth(iSnow) < zmin)
    case(rulesDependLayerIndex); removeLayer = (mLayerDepth(iSnow) < zminLayer(iSnow))
    case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
   end select ! (option to combine/sub-divide snow layers)
   !print*, 'in layerMerge: iSnow, mLayerDepth(iSnow), zminLayer(iSnow), removeLayer = ', iSnow, mLayerDepth(iSnow), zminLayer(iSnow), removeLayer

   ! check if need to remove a layer
   if(removeLayer)then

    ! ***** handle special case of a single layer
    if(nSnow==1)then
     ! set the variables defining "snow without a layer"
     ! NOTE: ignoring cold content!!! Need to fix later...
     scalarSnowDepth = mLayerDepth(1)
     scalarSWE       = (mLayerVolFracIce(1)*iden_ice + mLayerVolFracLiq(1)*iden_water)*mLayerDepth(1)
     ! remove the top layer from all model variable vectors
     ! NOTE: nSnow-1 = 0, so routine removes layer #1
     call rmLyAllVars(nSnow-1,nSnow,nLayers,err,cmessage)
     if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
     ! update the total number of layers
     indx_data%var(iLookINDEX%nLayers)%dat(1) = nSoil
     ! adjust coordinate variables
     call calcHeight(err,cmessage)
     ! exit the do loop (no more snow layers to remove)
     return
    endif

    ! ***** identify the layer to combine
    if(iSnow==1)then
     jSnow = iSnow+1  ! upper-most layer, combine with its lower neighbor
    elseif(iSnow==nSnow)then
     jSnow = nSnow-1  ! lower-most layer, combine with its upper neighbor
    else
     if(mLayerDepth(iSnow-1)<mLayerDepth(iSnow+1))then; jSnow = iSnow-1; else; jSnow = iSnow+1; endif
    endif

    ! ***** combine layers
    ! identify the layer closest to the surface
    kSnow=min(iSnow,jSnow)
    ! combine layer with identified neighbor
    call layer_combine(kSnow,err,cmessage)
    if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
    ! re-assign pointers to the layer depth
    mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
    ! update the number of snow layers
    nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
    !print*, 'removed layer: mLayerDepth = ', mLayerDepth
    !pause

    ! exit the loop to try again
    exit

   endif  ! (if layer is below the mass threshold)

   kSnow=iSnow ! ksnow is used for completion test, so include here

  end do ! (looping through snow layers)

  !print*, 'ksnow = ', ksnow

  ! exit if finished
  if(kSnow==nCheck)exit

 end do ! continuous do

 ! handle special case of > 5 layers in the CLM option
 if(nSnow > 5 .and. model_decisions(iLookDECISIONS%snowLayers)%iDecision == rulesDependLayerIndex)then
  ! initial check to ensure everything is wonderful in the universe
  if(nSnow /= 6)then; err=5; message=trim(message)//'special case of > 5 layers: expect only six layers'; return; endif
  ! combine 5th layer with layer below
  call layer_combine(5,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  ! re-assign pointers to the layer depth
  mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
  ! update the number of snow layers
  nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
  !print*, 'removed layer: mLayerDepth = ', mLayerDepth
  !pause
  if(nSnow /= 5)then; err=5; message=trim(message)//'special case of > 5 layers: expect have fixed matters so there are now just 5 layers'; return; endif
 endif

 ! check that there are no more than 5 layers in the CLM option
 if(model_decisions(iLookDECISIONS%snowLayers)%iDecision == rulesDependLayerIndex)then
  if(nSnow > 5)then
   message=trim(message)//'expect no more than 5 layers when combination/sub-division rules depend on the layer index (CLM option)'
   err=20; return
  endif
 endif

 end subroutine layerMerge


 ! ***********************************************************************************************************
 ! new subroutine: combine snow layers and re-compute model state variables
 ! ***********************************************************************************************************
 subroutine layer_combine(iSnow,err,message)
 ! combines layer iSnow with iSnow+1
 USE multiconst,only:&
                     iden_ice,       & ! intrinsic density of ice             (kg m-3)
                     iden_water        ! intrinsic density of liquid water    (kg m-3)
 USE var_derive_module,only:calcHeight ! module to calculate height at layer interfaces and layer mid-point
 USE snow_utils_module,only:fracliquid                                          ! compute fraction of liquid water
 USE convE2Temp_module,only:E2T_nosoil,temp2ethpy                               ! convert temperature to enthalpy
 USE data_struc,only:mpar_meta,forc_meta,mvar_meta,indx_meta                    ! metadata
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 implicit none
 ! dummies
 integer(i4b),intent(in)             :: iSnow                    ! index of top layer to combine
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! local pointers to model parameters
 real(dp),pointer                    :: snowfrz_scale            ! scaling parameter for the freezing curve for snow (K-1)
 ! local pointers to model state variables
 real(dp),pointer                    :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer                    :: mLayerDepth(:)           ! depth of each layer (m)
 real(dp),pointer                    :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                    :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 ! local variables
 character(len=256)                  :: cmessage                 ! error message for downwind routine
 integer(i4b)                        :: nSnow                    ! number of snow layers
 integer(i4b)                        :: nSoil                    ! number of soil layers
 integer(i4b),pointer                :: nLayers                  ! total number of layers
 real(dp)                            :: massIce(2)               ! mass of ice in the two layers identified for combination (kg m-2)
 real(dp)                            :: massLiq(2)               ! mass of liquid water in the two layers identified for combination (kg m-2)
 real(dp)                            :: bulkDenWat(2)            ! bulk density if total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(dp)                            :: cBulkDenWat              ! combined bulk density of total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(dp)                            :: cTemp                    ! combined layer temperature
 real(dp)                            :: cDepth                   ! combined layer depth
 real(dp)                            :: cVolFracIce              ! combined layer volumetric fraction of ice
 real(dp)                            :: cVolFracLiq              ! combined layer volumetric fraction of liquid water
 real(dp)                            :: l1Enthalpy,l2Enthalpy    ! enthalpy in the two layers identified for combination (J m-3)
 real(dp)                            :: cEnthalpy                ! combined layer enthalpy (J m-3)
 real(dp)                            :: fLiq                     ! fraction of liquid water at the combined temperature cTemp
 real(dp),parameter                  :: eTol=1.e-4_dp            ! tolerance for the enthalpy-->temperature conversion (J m-3)
 ! initialize error control
 err=0; message="layer_combine/"

 ! identify the number of snow layers, and the total number of layers
 nLayers => indx_data%var(iLookINDEX%nLayers)%dat(1) 
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
 print*, '***** removing layer', iSnow

 ! ***** compute combined model state variables
 ! assign pointer to the freezing curve parameter
 snowfrz_scale   => mpar_data%var(iLookPARAM%snowfrz_scale)            ! scaling parameter for the freezing curve for snow (K-1)
 ! assign pointers to model state variables
 mLayerTemp       => mvar_data%var(iLookMVAR%mLayerTemp)%dat           ! temperature of each layer (K)
 mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
 mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat     ! volumetric fraction of ice in each layer  (-)
 mLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat     ! volumetric fraction of liquid water in each layer (-)
 !write(*,'(a,1x,6(i12,1x))')   'layer_combine, before merge: indx_data%var(iLookINDEX%layerType)%dat(1:6) = ', &
 !                                                            indx_data%var(iLookINDEX%layerType)%dat(1:6)
 !write(*,'(a,1x,6(f12.3,1x))') 'layer_combine, before merge: mvar_data%var(iLookMVAR%mLayerTemp)%dat(1:6) = ', &
 !                                                            mvar_data%var(iLookMVAR%mLayerTemp)%dat(1:6)
 ! compute combined depth
 cDepth       = mLayerDepth(isnow) + mLayerDepth(isnow+1)
 !write(*,'(a,10(f12.4,1x))') 'cDepth = ', cDepth
 ! compute mass of each layer (kg m-2)
 massIce(1:2) = iden_ice*mLayerVolFracIce(iSnow:iSnow+1)*mLayerDepth(iSnow:iSnow+1)
 massLiq(1:2) = iden_water*mLayerVolFracLiq(iSnow:iSnow+1)*mLayerDepth(iSnow:iSnow+1)
 !write(*,'(a,10(f12.4,1x))') 'massIce = ', massIce
 !write(*,'(a,10(f12.4,1x))') 'massLiq = ', massLiq
 ! compute bulk density of water (kg m-3)
 bulkDenWat(1:2) = (massIce(1:2) + massLiq(1:2))/mLayerDepth(iSnow:iSnow+1)
 cBulkDenWat     = (mLayerDepth(isnow)*bulkDenWat(1) + mLayerDepth(isnow+1)*bulkDenWat(2))/cDepth
 !write(*,'(a,10(f12.4,1x))') 'bulkDenWat  = ', bulkDenWat
 !write(*,'(a,10(f12.4,1x))') 'cBulkDenWat = ', cBulkDenWat
 ! compute enthalpy for each layer (J m-3)
 l1Enthalpy  = temp2ethpy(mLayerTemp(iSnow),  BulkDenWat(1),snowfrz_scale)
 l2Enthalpy  = temp2ethpy(mLayerTemp(iSnow+1),BulkDenWat(2),snowfrz_scale)
 !write(*,'(a,10(e20.9,1x))') 'l1Enthalpy = ', l1Enthalpy
 !write(*,'(a,10(e20.9,1x))') 'l2Enthalpy = ', l2Enthalpy
 ! compute combined enthalpy (J m-3)
 cEnthalpy   = (mLayerDepth(isnow)*l1Enthalpy + mLayerDepth(isnow+1)*l2Enthalpy)/cDepth
 ! convert enthalpy (J m-3) to temperature (K)
 call E2T_nosoil(cEnthalpy,cBulkDenWat,snowfrz_scale,cTemp,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! test enthalpy conversion
 if(abs(temp2ethpy(cTemp,cBulkDenWat,snowfrz_scale)/cBulkDenWat - cEnthalpy/cBulkDenWat) > eTol)then
  write(*,'(a,1x,f12.5,1x,2(e20.10,1x))') 'enthalpy test', cBulkDenWat, temp2ethpy(cTemp,cBulkDenWat,snowfrz_scale)/cBulkDenWat, cEnthalpy/cBulkDenWat
  message=trim(message)//'problem with enthalpy-->temperature conversion'
  err=20; return
 endif
 ! check temperature is within the two temperatures
 if(cTemp > max(mLayerTemp(iSnow),mLayerTemp(iSnow+1)))then; err=20; message=trim(message)//'merged temperature > max(temp1,temp2)'; return; endif
 if(cTemp < min(mLayerTemp(iSnow),mLayerTemp(iSnow+1)))then; err=20; message=trim(message)//'merged temperature < min(temp1,temp2)'; return; endif
 ! compute volumetric fraction of liquid water
 fLiq = fracLiquid(cTemp,snowfrz_scale)
 ! compute volumetric fraction of ice and liquid water
 cVolFracLiq =          fLiq *cBulkDenWat/iden_water
 cVolFracIce = (1._dp - fLiq)*cBulkDenWat/iden_ice
 !cVolFracIce = (massIce(1) + massIce(2))/(cDepth*iden_ice)
 !cVolFracLiq = (massLiq(1) + massLiq(2))/(cDepth*iden_water)
 ! remove a model layer from all model variable vectors
 call rmLyAllVars(iSnow,nSnow,nLayers,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! define the combined layer as snow
 indx_data%var(iLookINDEX%layerType)%dat(iSnow) = ix_snow
 ! update the total number of layers
 nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
 nSoil = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)
 indx_data%var(iLookINDEX%nLayers)%dat(1) = nSnow + nSoil

 ! ***** put state variables for the combined layer in the appropriate place
 mvar_data%var(iLookMVAR%mLayerTemp)%dat(iSnow)       = cTemp
 mvar_data%var(iLookMVAR%mLayerDepth)%dat(iSnow)      = cDepth
 mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(iSnow) = cVolFracIce
 mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(iSnow) = cVolFracLiq

 ! ***** adjust coordinate variables
 call calcHeight(err,cmessage) 
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
 !write(*,'(a,1x,6(i12,1x))')   'layer_combine, after merge:  indx_data%var(iLookINDEX%layerType)%dat(1:6) = ', &
 !                                                            indx_data%var(iLookINDEX%layerType)%dat(1:6)
 !write(*,'(a,1x,6(f12.3,1x))') 'layer_combine, after merge:  mvar_data%var(iLookMVAR%mLayerTemp)%dat(1:6) = ', &
 !                                                            mvar_data%var(iLookMVAR%mLayerTemp)%dat(1:6)
 end subroutine layer_combine

 ! ***********************************************************************************************************
 ! new subroutine: reduce the length of the vectors in data structures
 ! ***********************************************************************************************************
 subroutine rmLyAllVars(iSnow,nSnow,nLayers,err,message)
 ! removes layer "iSnow+1" and sets layer "iSnow" to a missing value
 ! (layer "iSnow" will be filled with a combined layer later)
 USE data_struc,only:mvar_meta,indx_meta   ! metadata
 USE data_struc,only:mvar_data,indx_data   ! data structures
 USE var_lookup,only:iLookMVAR,iLookINDEX  ! named variables for structure elements
 implicit none
 ! dummies
 integer(i4b),intent(in)            :: iSnow       ! new layer
 integer(i4b),intent(in)            :: nSnow       ! number of snow layers
 integer(i4b),intent(in)            :: nLayers     ! total number of layers
 integer(i4b),intent(out)           :: err         ! error code
 character(*),intent(out)           :: message     ! error message
 ! locals
 character(len=256)                 :: cmessage    ! error message for downwind routine
 integer(i4b)                       :: ix_lower    ! lower bound of the vector
 integer(i4b)                       :: ix_upper    ! upper bound of the vector
 integer(i4b)                       :: ivar
 ! initialize error control
 err=0; message="rmLyAllVars/"
 ! ***** loop through model variables and remove one layer
 do ivar=1,size(mvar_data%var)
  ! define bounds
  select case(trim(mvar_meta(ivar)%vartype))
   case('midSnow'); ix_lower=1; ix_upper=nSnow
   case('midToto'); ix_lower=1; ix_upper=nLayers
   case('ifcSnow'); ix_lower=0; ix_upper=nSnow
   case('ifcToto'); ix_lower=0; ix_upper=nLayers
   case default; cycle  ! no need to remove soil layers
  end select
  ! remove a layer for a model variable vector
  call removeOneLayer(mvar_data%var(ivar)%dat,ix_lower,ix_upper,iSnow,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 end do
 ! adjust the layer type (ix_soil or ix_snow) accordingly
 call removeOneLayer(indx_data%var(iLookINDEX%layerType)%dat,1,nLayers,iSnow,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 end subroutine rmLyAllVars


 ! ***********************************************************************************************************
 ! new subroutine: combine snow layers and reduce the length of the vectors in data structures
 ! ***********************************************************************************************************
 ! double precision
 subroutine removeOneLayer_rv(datavec,ix_lower,ix_upper,iSnow,err,message)
 ! Removes layer iSnow+1 from the input vector, set layer iSnow to a missing value,
 ! and reduce size of the vector by 1 element
 implicit none
 ! dummies
 real(dp),pointer,intent(inout)     :: datavec(:)  ! the original and the new vector
 integer(i4b),intent(in)            :: ix_lower    ! lower bound of the old vector
 integer(i4b),intent(in)            :: ix_upper    ! upper bound of the old vector
 integer(i4b),intent(in)            :: iSnow       ! new layer
 integer(i4b),intent(out)           :: err         ! error code
 character(*),intent(out)           :: message     ! error message
 ! locals
 real(dp)                           :: tempvec(ix_lower:ix_upper-1)  ! temporary vector
 real(dp),parameter                 :: missingDouble=-9999._dp
 ! initialize error control
 err=0; message='RemoveOneLayer_rv/'
 ! check the data vector is associated
 if(.not.associated(datavec))then; err=20; message='data vector is not associated'; return; endif
 ! copy elements across to the temporary vector
 if(iSnow>ix_lower) tempvec(ix_lower:iSnow-1) = datavec(ix_lower:iSnow-1)
 if(iSnow+1<ix_upper) tempvec(iSnow+1:ix_upper-1) = datavec(iSnow+2:ix_upper)  ! skip iSnow+1
 ! set the iSnow value to a stupid number -- fill in later
 if(iSnow>=ix_lower) tempvec(iSnow) = missingDouble
 ! adjust size of the data vector (one less element)
 deallocate(datavec,stat=err)
 if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; endif
 allocate(datavec(ix_lower:ix_upper-1),stat=err)
 if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; endif
 ! copy elements across to the data vector
 datavec=tempvec
 end subroutine RemoveOneLayer_rv

 ! integer
 subroutine RemoveOneLayer_iv(datavec,ix_lower,ix_upper,iSnow,err,message)
 ! Removes layer iSnow+1 from the input vector, set layer iSnow to a missing value,
 ! and reduce size of the vector by 1 element
 implicit none
 ! dummies
 integer(i4b),pointer,intent(inout) :: datavec(:)  ! the original and the new vector
 integer(i4b),intent(in)            :: ix_lower    ! lower bound of the old vector
 integer(i4b),intent(in)            :: ix_upper    ! upper bound of the old vector
 integer(i4b),intent(in)            :: iSnow       ! new layer
 integer(i4b),intent(out)           :: err         ! error code
 character(*),intent(out)           :: message     ! error message
 ! locals
 integer(i4b)                       :: tempvec(ix_lower:ix_upper-1)  ! temporary vector
 integer(i4b),parameter             :: missingInteger=-9999
 ! initialize error control
 err=0; message='RemoveOneLayer_iv/'
 ! check the data vector is associated
 if(.not.associated(datavec))then; err=20; message='data vector is not associated'; return; endif
 ! copy elements across to the temporary vector
 if(iSnow>ix_lower) tempvec(ix_lower:iSnow-1) = datavec(ix_lower:iSnow-1)
 if(iSnow+1<ix_upper) tempvec(iSnow+1:ix_upper-1) = datavec(iSnow+2:ix_upper)  ! skip iSnow+1
 ! set the iSnow value to a stupid number -- fill in later
 if(iSnow>=ix_lower) tempvec(iSnow) = missingInteger
 ! adjust size of the data vector (one less element)
 deallocate(datavec,stat=err)
 if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; endif
 allocate(datavec(ix_lower:ix_upper-1),stat=err)
 if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; endif
 ! copy elements across to the data vector
 datavec=tempvec
 end subroutine RemoveOneLayer_iv


end module layerMerge_module
