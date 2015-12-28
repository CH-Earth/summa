! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module layerMerge_module
! data types
USE nrtype
! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,   & ! number of snow layers
                    nSoil,   & ! number of soil layers
                    nLayers    ! total number of layers
! access named variables for snow and soil
USE data_struc,only:ix_soil,ix_snow            ! named variables for snow and soil
! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)
! look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
! provide access to external modules
USE var_derive_module,only:calcHeight ! module to calculate height at layer interfaces and layer mid-point
implicit none
private
public::layerMerge
interface removeOneLayer
 module procedure removeOneLayer_rv, removeOneLayer_iv
end interface removeOneLayer

contains


 ! *****************************************************************************************************************
 ! public subroutine layerMerge: merge layers if the thickness is less than zmin
 ! *****************************************************************************************************************
 subroutine layerMerge(&
                       ! input/output: model data structures
                       model_decisions,             & ! intent(in):    model decisions
                       mpar_data,                   & ! intent(in):    model parameters
                       indx_data,                   & ! intent(inout): type of each layer
                       mvar_data,                   & ! intent(inout): model variables for a local HRU
                       ! output: error control
                       err,message)                   ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! access the derived types to define the data structures
 USE data_struc,only:&
                     var_d,            & ! data vector (dp)
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! access named variables defining elements in the data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_d),intent(in)          :: mpar_data           ! model parameters
 type(var_ilength),intent(inout) :: indx_data           ! type of each layer
 type(var_dlength),intent(inout) :: mvar_data           ! model variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! model index variables
 ! NOTE: use pointers because the dimension length changes
 integer(i4b),pointer            :: layerType(:)        ! type of the layer (ix_soil or ix_snow)
 ! model state variables
 ! NOTE: use pointers because the dimension length changes
 real(dp),pointer                :: mLayerDepth(:)      ! depth of each layer (m)
 real(dp),pointer                :: mLayerVolFracIce(:) ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                :: mLayerVolFracLiq(:) ! volumetric fraction of liquid water in each layer (-)
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 real(dp),dimension(5)           :: zminLayer           ! minimum layer depth in each layer (m)
 logical(lgt)                    :: removeLayer         ! flag to indicate need to remove a layer
 integer(i4b)                    :: nCheck              ! number of layers to check for combination
 integer(i4b)                    :: iSnow               ! index of snow layers (looping)
 integer(i4b)                    :: jSnow               ! index of snow layer identified for combination with iSnow
 integer(i4b)                    :: kSnow               ! index of the upper layer of the two layers identified for combination
 ! initialize error control
 err=0; message="layerMerge/"
 ! --------------------------------------------------------------------------------------------------------
 ! associate variables to the data structures
 associate(&

 ! model decisions
 ix_snowLayers    => model_decisions(iLookDECISIONS%snowLayers)%iDecision, & ! decision for snow combination

 ! model parameters (control the depth of snow layers)
 zmin             => mpar_data%var(iLookPARAM%zmin),                       & ! minimum layer depth (m)
 zminLayer1       => mpar_data%var(iLookPARAM%zminLayer1),                 & ! minimum layer depth for the 1st (top) layer (m)
 zminLayer2       => mpar_data%var(iLookPARAM%zminLayer2),                 & ! minimum layer depth for the 2nd layer (m)
 zminLayer3       => mpar_data%var(iLookPARAM%zminLayer3),                 & ! minimum layer depth for the 3rd layer (m)
 zminLayer4       => mpar_data%var(iLookPARAM%zminLayer4),                 & ! minimum layer depth for the 4th layer (m)
 zminLayer5       => mpar_data%var(iLookPARAM%zminLayer5),                 & ! minimum layer depth for the 5th (bottom) layer (m)

 ! diagnostic scalar variables
 scalarSnowDepth  => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),      & ! total snow depth (m)
 scalarSWE        => mvar_data%var(iLookMVAR%scalarSWE)%dat(1)             & ! SWE (kg m-2)

 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------

 ! point to the model index structures
 layerType        => indx_data%var(iLookINDEX%layerType)%dat                 ! layer type (ix_soil or ix_snow)

 ! point to the model state variables
 mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat                ! depth of each layer (m)
 mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat           ! volumetric fraction of ice in each layer  (-)
 mLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat           ! volumetric fraction of liquid water in each layer (-)

 ! identify algorithmic control parameters to syb-divide and combine snow layers
 zminLayer = (/zminLayer1, zminLayer2, zminLayer3, zminLayer4, zminLayer5/)

 kSnow=0 ! initialize first layer to test (top layer)
 do ! attempt to remove multiple layers in a single time step (continuous do loop with exit clause)

  ! special case of >5 layers: add an offset to use maximum threshold from layer above
  if(ix_snowLayers == rulesDependLayerIndex .and. nSnow > 5)then
   nCheck=5
  else
   nCheck=nSnow
  endif

  ! loop through snow layers
  do iSnow=kSnow+1,nCheck

   ! check if the layer depth is less than the depth threshold
   select case(ix_snowLayers)
    case(sameRulesAllLayers);    removeLayer = (mLayerDepth(iSnow) < zmin)
    case(rulesDependLayerIndex); removeLayer = (mLayerDepth(iSnow) < zminLayer(iSnow))
    case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
   end select ! (option to combine/sub-divide snow layers)

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
     call rmLyAllVars(mvar_data,indx_data,nSnow-1,err,cmessage)
     if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
     ! update the total number of layers
     nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
     nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)
     nLayers = nSnow + nSoil
     ! save the number of layers
     indx_data%var(iLookINDEX%nSnow)%dat(1)   = nSnow
     indx_data%var(iLookINDEX%nSoil)%dat(1)   = nSoil
     indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers
     ! update coordinate variables
     call calcHeight(&
                     ! input/output: data structures
                     indx_data,   & ! intent(in): layer type
                     mvar_data,   & ! intent(inout): model variables for a local HRU
                     ! output: error control
                     err,cmessage)
     if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
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
    call layer_combine(mpar_data,mvar_data,indx_data,kSnow,err,cmessage)
    if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

    ! re-assign pointers to the layer depth
    mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
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
 if(nSnow > 5 .and. ix_snowLayers == rulesDependLayerIndex)then
  ! initial check to ensure everything is wonderful in the universe
  if(nSnow /= 6)then; err=5; message=trim(message)//'special case of > 5 layers: expect only six layers'; return; endif
  ! combine 5th layer with layer below
  call layer_combine(mpar_data,mvar_data,indx_data,5,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  ! re-assign pointers to the layer depth
  mLayerDepth => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
  !print*, 'removed layer: mLayerDepth = ', mLayerDepth
  !pause
  if(nSnow /= 5)then; err=5; message=trim(message)//'special case of > 5 layers: expect to reduced layers to exactly 5'; return; endif
 endif

 ! check that there are no more than 5 layers in the CLM option
 if(ix_snowLayers == rulesDependLayerIndex)then
  if(nSnow > 5)then
   message=trim(message)//'expect no more than 5 layers when combination/sub-division rules depend on the layer index (CLM option)'
   err=20; return
  endif
 endif

 ! end association to variables in the data structure
 end associate

 end subroutine layerMerge


 ! ***********************************************************************************************************
 ! private subroutine layer_combine: combine snow layers and re-compute model state variables
 ! ***********************************************************************************************************
 ! combines layer iSnow with iSnow+1
 ! ***********************************************************************************************************
 subroutine layer_combine(mpar_data,mvar_data,indx_data,iSnow,err,message)
 ! provide access to variables in the data structures
 USE data_struc,only:var_d                    ! data structures with fixed dimension
 USE data_struc,only:var_ilength,var_dlength  ! data vectors with variable length dimension
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX  ! named variables for structure elements
 ! provide access to external modules
 USE snow_utils_module,only:fracliquid                                          ! compute fraction of liquid water
 USE convE2Temp_module,only:E2T_nosoil,temp2ethpy                               ! convert temperature to enthalpy
 implicit none
 ! ------------------------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(var_d),intent(in)          :: mpar_data ! model parameters
 type(var_dlength),intent(inout) :: mvar_data ! model variables for a local HRU
 type(var_ilength),intent(inout) :: indx_data ! type of model layer
 ! input: snow layer indices
 integer(i4b),intent(in)         :: iSnow     ! index of top layer to combine
 ! output: error control
 integer(i4b),intent(out)        :: err       ! error code
 character(*),intent(out)        :: message   ! error message
 ! ------------------------------------------------------------------------------------------------------------
 ! model state variables
 ! NOTE: these are defined as pointers because the length of the data dimension changes
 real(dp),pointer                :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer                :: mLayerDepth(:)           ! depth of each layer (m)
 real(dp),pointer                :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 ! ------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=256)              :: cmessage                 ! error message for downwind routine
 real(dp)                        :: massIce(2)               ! mass of ice in the two layers identified for combination (kg m-2)
 real(dp)                        :: massLiq(2)               ! mass of liquid water in the two layers identified for combination (kg m-2)
 real(dp)                        :: bulkDenWat(2)            ! bulk density if total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(dp)                        :: cBulkDenWat              ! combined bulk density of total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(dp)                        :: cTemp                    ! combined layer temperature
 real(dp)                        :: cDepth                   ! combined layer depth
 real(dp)                        :: cVolFracIce              ! combined layer volumetric fraction of ice
 real(dp)                        :: cVolFracLiq              ! combined layer volumetric fraction of liquid water
 real(dp)                        :: l1Enthalpy,l2Enthalpy    ! enthalpy in the two layers identified for combination (J m-3)
 real(dp)                        :: cEnthalpy                ! combined layer enthalpy (J m-3)
 real(dp)                        :: fLiq                     ! fraction of liquid water at the combined temperature cTemp
 real(dp),parameter              :: eTol=1.e-4_dp            ! tolerance for the enthalpy-->temperature conversion (J m-3)
 ! initialize error control
 err=0; message="layer_combine/"

 ! associate variables with information in the data structures
 associate(&
 snowfrz_scale => mpar_data%var(iLookPARAM%snowfrz_scale)      & ! scaling parameter for the freezing curve for snow (K-1)
 ) ! end associate block

 ! ***** compute combined model state variables
 ! assign pointers to model state variables
 mLayerTemp       => mvar_data%var(iLookMVAR%mLayerTemp)%dat           ! temperature of each layer (K)
 mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
 mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat     ! volumetric fraction of ice in each layer  (-)
 mLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat     ! volumetric fraction of liquid water in each layer (-)

 ! compute combined depth
 cDepth       = mLayerDepth(isnow) + mLayerDepth(isnow+1)

 ! compute mass of each layer (kg m-2)
 massIce(1:2) = iden_ice*mLayerVolFracIce(iSnow:iSnow+1)*mLayerDepth(iSnow:iSnow+1)
 massLiq(1:2) = iden_water*mLayerVolFracLiq(iSnow:iSnow+1)*mLayerDepth(iSnow:iSnow+1)

 ! compute bulk density of water (kg m-3)
 bulkDenWat(1:2) = (massIce(1:2) + massLiq(1:2))/mLayerDepth(iSnow:iSnow+1)
 cBulkDenWat     = (mLayerDepth(isnow)*bulkDenWat(1) + mLayerDepth(isnow+1)*bulkDenWat(2))/cDepth

 ! compute enthalpy for each layer (J m-3)
 l1Enthalpy  = temp2ethpy(mLayerTemp(iSnow),  BulkDenWat(1),snowfrz_scale)
 l2Enthalpy  = temp2ethpy(mLayerTemp(iSnow+1),BulkDenWat(2),snowfrz_scale)

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

 ! remove a model layer from all model variable vectors
 call rmLyAllVars(mvar_data,indx_data,iSnow,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! define the combined layer as snow
 indx_data%var(iLookINDEX%layerType)%dat(iSnow) = ix_snow

 ! update the total number of layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)
 nLayers = nSnow + nSoil

 ! save the number of layers in the data structures
 indx_data%var(iLookINDEX%nSnow)%dat(1)   = nSnow
 indx_data%var(iLookINDEX%nSoil)%dat(1)   = nSoil
 indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers

 ! ***** put state variables for the combined layer in the appropriate place
 mvar_data%var(iLookMVAR%mLayerTemp)%dat(iSnow)       = cTemp
 mvar_data%var(iLookMVAR%mLayerDepth)%dat(iSnow)      = cDepth
 mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(iSnow) = cVolFracIce
 mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(iSnow) = cVolFracLiq

 ! ***** adjust coordinate variables
 call calcHeight(&
                 ! input/output: data structures
                 indx_data,   & ! intent(in): layer type
                 mvar_data,   & ! intent(inout): model variables for a local HRU
                 ! output: error control
                 err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! end association to data structures
 end associate

 end subroutine layer_combine


 ! ***********************************************************************************************************
 ! private subroutine rmLyAllVars: reduce the length of the vectors in data structures
 ! ***********************************************************************************************************
 ! removes layer "iSnow+1" and sets layer "iSnow" to a missing value
 ! (layer "iSnow" will be filled with a combined layer later)
 ! ***********************************************************************************************************
 subroutine rmLyAllVars(mvar_data,indx_data,iSnow,err,message)
 USE data_struc,only:mvar_meta                ! metadata
 USE data_struc,only:var_ilength,var_dlength  ! data vectors with variable length dimension
 USE var_lookup,only:iLookMVAR,iLookINDEX     ! named variables for structure elements
 implicit none
 ! input/output: data structures
 type(var_ilength),intent(inout) :: indx_data ! type of model layer
 type(var_dlength),intent(inout) :: mvar_data ! model variables for a local HRU
 ! input: snow layer indices
 integer(i4b),intent(in)         :: iSnow     ! new layer
 ! output: error control
 integer(i4b),intent(out)        :: err       ! error code
 character(*),intent(out)        :: message   ! error message
 ! locals
 character(len=256)              :: cmessage  ! error message for downwind routine
 integer(i4b)                    :: ix_lower  ! lower bound of the vector
 integer(i4b)                    :: ix_upper  ! upper bound of the vector
 integer(i4b)                    :: ivar
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
   case default; cycle  ! no need to remove soil layers or scalar variables
  end select
  ! remove a layer for a model variable vector
  call removeOneLayer(mvar_data%var(ivar)%dat,ix_lower,ix_upper,iSnow,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 end do
 ! adjust the layer type (ix_soil or ix_snow) accordingly
 call removeOneLayer(indx_data%var(iLookINDEX%layerType)%dat,1,nLayers,iSnow,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 end subroutine rmLyAllVars


 ! *****************************************************************************************************************
 ! private subroutine removeOneLayer: combine snow layers and reduce the length of the vectors in data structures
 ! *****************************************************************************************************************
 ! double precision
 ! Removes layer iSnow+1 from the input vector, set layer iSnow to a missing value,
 ! and reduce size of the vector by 1 element
 ! *****************************************************************************************************************
 subroutine removeOneLayer_rv(datavec,ix_lower,ix_upper,iSnow,err,message)
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


 ! *****************************************************************************************************************
 ! private subroutine RemoveOneLayer_iv: combine snow layers and reduce the length of the vectors in data structures
 ! *****************************************************************************************************************
 ! integer
 ! Removes layer iSnow+1 from the input vector, set layer iSnow to a missing value,
 ! and reduce size of the vector by 1 element
 ! *****************************************************************************************************************
 subroutine RemoveOneLayer_iv(datavec,ix_lower,ix_upper,iSnow,err,message)
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
