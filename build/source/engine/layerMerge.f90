! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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
USE nr_type

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

USE globalData,only:icefrz_mult     ! freezing curve scaling factor multipier of snow to ice, closer to a step function since ice does not hold water

! access named variables for layers
USE globalData,only:iname_snow        ! named variables for snow
USE globalData,only:iname_soil        ! named variables for soil
USE globalData,only:iname_glce        ! named variables for glacier ice
USE globalData,only:iname_lake        ! named variables for lake

! access metadata
USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta   ! metadata

! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water,     & ! intrinsic density of liquid water    (kg m-3)
                    Cp_ice,         & ! specific heat of ice                 (J kg-1 K-1)
                    Cp_water          ! specific heat of liquid water        (J kg-1 K-1)

! access the derived types to define the data structures
USE data_types,only:&
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! access named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookPROG,iLookINDEX  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                   ! named variables for elements of the decision structure

! look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index

! provide access to external modules
USE var_derive_module,only:calcHeight ! module to calculate height at layer interfaces and layer mid-point

! privacy
implicit none
private
public::layerMerge

contains


 ! *****************************************************************************************************************
 ! public subroutine layerMerge: merge layers if the thickness is less than zmin
 ! *****************************************************************************************************************
 subroutine layerMerge(&
                       ! input/output: model data structures
                       maxLayers,                   & ! intent(in):    maximum number of snow/firn/ice layers
                       tooMuchMelt,                 & ! intent(in):    flag to force merge of snow layers
                       model_decisions,             & ! intent(in):    model decisions
                       mpar_data,                   & ! intent(in):    model parameters
                       indx_data,                   & ! intent(inout): type of each layer
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,                   & ! intent(inout): model fluxes for a local HRU
                       ! output
                       mergedLayers,                & ! intent(out): flag to denote that layers were merged
                       err,message)                   ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 integer(i4b),intent(in)          :: maxLayers           ! maximum number of snow/firn/ice layers
 logical(lgt),intent(in)          :: tooMuchMelt         ! flag to denote that ice is insufficient to support melt
 type(model_options),intent(in)   :: model_decisions(:)  ! model decisions
 type(var_dlength),intent(in)     :: mpar_data           ! model parameters
 type(var_ilength),intent(inout)  :: indx_data           ! type of each layer
 type(var_dlength),intent(inout)  :: prog_data           ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout)  :: diag_data           ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout)  :: flux_data           ! model flux variables
 ! output
 logical(lgt),intent(out)         :: mergedLayers        ! flag to denote that layers were merged
 integer(i4b),intent(out)         :: err                 ! error code
 character(*),intent(out)         :: message             ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)               :: cmessage            ! error message of downwind routine
 real(rkind),dimension(maxLayers) :: zminLayer           ! minimum layer depth in each layer (m)
 real(rkind),dimension(5)         :: zminLayer_param     ! minimum layer depth in each layer (m) that has been set in the model parameters
 logical(lgt)                     :: removeLayer         ! flag to indicate need to remove a layer
 integer(i4b)                     :: nCheck              ! number of layers to check for combination
 integer(i4b)                     :: iLayer              ! layer index
 integer(i4b)                     :: jLayer              ! index of layer identified for combination with iLayer
 integer(i4b)                     :: kLayer              ! index of the upper layer of the two layers identified for combination
 integer(i4b)                     :: nSnow               ! number of snow layers
 integer(i4b)                     :: nLake               ! number of lake layers
 integer(i4b)                     :: nSoil               ! number of soil layers
 integer(i4b)                     :: nGlce               ! number of glacier ice layers
 integer(i4b)                     :: nLayers             ! total number of layers
 logical(lgt)                     :: doGlac              ! flag to denote that merging glacier ice
 integer(i4b)                     :: topLayer            ! index of the top layer of snow/ice
 integer(i4b)                     :: botLayer            ! index of the bottom layer of snow/ice
 ! --------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="layerMerge/"
 ! --------------------------------------------------------------------------------------------------------
 ! associate variables to the data structures
 associate(&

 ! model decisions
 ix_snowLayers    => model_decisions(iLookDECISIONS%snowLayers)%iDecision, & ! decision for snow combination

 ! model parameters (control the depth of snow layers)
 zmin             => mpar_data%var(iLookPARAM%zmin)%dat(1),                & ! minimum layer depth (m)
 zminLayer1       => mpar_data%var(iLookPARAM%zminLayer1)%dat(1),          & ! minimum layer depth for the 1st (top) layer (m)
 zminLayer2       => mpar_data%var(iLookPARAM%zminLayer2)%dat(1),          & ! minimum layer depth for the 2nd layer (m)
 zminLayer3       => mpar_data%var(iLookPARAM%zminLayer3)%dat(1),          & ! minimum layer depth for the 3rd layer (m)
 zminLayer4       => mpar_data%var(iLookPARAM%zminLayer4)%dat(1),          & ! minimum layer depth for the 4th layer (m)
 zminLayer5       => mpar_data%var(iLookPARAM%zminLayer5)%dat(1),          & ! minimum layer depth for the 5th (bottom) layer (m)
 noThetaChange    => indx_data%var(iLookINDEX%noThetaChange)%dat(1),       & ! number of layers with no change in total water content (bottom layers)

 ! diagnostic scalar variables
 scalarSnowDepth  => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),      & ! total snow depth (m)
 scalarSWE        => prog_data%var(iLookPROG%scalarSWE)%dat(1)             & ! SWE (kg m-2)

 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------

 ! identify algorithmic control parameters to sub-divide and combine layers
 zminLayer_param = (/zminLayer1, zminLayer2, zminLayer3, zminLayer4, zminLayer5/)
 if (maxLayers <= 5) then
    zminLayer = zminLayer_param(1:maxLayers)
 else
    zminLayer(1:5) = zminLayer_param
    do iLayer=6,maxLayers
      zminLayer(iLayer) = zminLayer(iLayer-1)*2._rkind
    end do
 end if

 ! intialize the modified layers flag
 mergedLayers=.false.

 ! initialize the number of layers
 nSnow    = indx_data%var(iLookINDEX%nSnow)%dat(1)
 nLake    = indx_data%var(iLookINDEX%nLake)%dat(1)
 nSoil    = indx_data%var(iLookINDEX%nSoil)%dat(1)
 nGlce    = indx_data%var(iLookINDEX%nGlce)%dat(1)
 nLayers  = indx_data%var(iLookINDEX%nLayers)%dat(1)
 ! NOTE: lake layers are never merged here: the ice cover of a lake is created, thickened and broken up by lakeIceCover
 doGlac=.false. ! initialize flag for glacier ice
 if (nSnow+nLake==0 .and. nGlce>0) then
   doGlac=.true.
   topLayer=nSnow+nLake+nSoil+1
   botLayer=nSnow+nLake+nSoil+nGlce-noThetaChange
 else
   topLayer=1
   botLayer=nSnow
 end if
 kLayer = topLayer-1 ! initialize first layer to test (top layer)
 do ! attempt to remove multiple layers in a single time step (continuous do loop with exit clause)

  ! set number of layers to check
  if(doGlac)then
    nCheck=botLayer
  elseif(ix_snowLayers == rulesDependLayerIndex .and. nSnow > maxLayers)then
    ! special case of >maxLayers layers: add an offset to use maximum threshold from layer above
    nCheck=maxLayers
  else
    nCheck=nSnow
  end if

  ! loop through snow/firn/ice layers
  do iLayer=kLayer+1,nCheck

    ! NOTE: the layer state is read through the accessors, since the layer vectors are reallocated as layers are removed
    ! check if the layer depth is less than the depth threshold
    if(doGlac)then
      removeLayer = (depth(iLayer) < zminLayer(iLayer-nSnow-nLake-nSoil))
    else
      select case(ix_snowLayers)
        case(sameRulesAllLayers);    removeLayer = (depth(iLayer) < zmin)
        case(rulesDependLayerIndex); removeLayer = (depth(iLayer) < zminLayer(iLayer))
        case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
      end select ! (option to combine/sub-divide snow layers)
    end if

    ! check if we have too much melt
    ! NOTE: assume that this is the top snow layer; need more trickery to relax this assumption
    if(tooMuchMelt .and. iLayer==topLayer) removeLayer = .true.

    ! check if need to remove a layer
    if(removeLayer)then

      ! flag that we modified a layer
      mergedLayers=.true.

      ! ***** handle special case of a single layer
      if(nSnow==1)then ! here assuming would not be merging glacier ice layers if had snow
        ! set the variables defining "snow without a layer"
        ! NOTE: ignoring cold content!!! Need to fix later...
        scalarSnowDepth = depth(1)
        scalarSWE       = (ice(1)*iden_ice + liq(1)*iden_water)*depth(1)
        ! remove the top layer from all model variable vectors
        ! NOTE: nSnow-1 = 0, so routine removes layer #1
        call rmLyAllVars(doGlac,prog_data,prog_meta,nSnow-1,nSnow,nGlce,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        call rmLyAllVars(doGlac,diag_data,diag_meta,nSnow-1,nSnow,nGlce,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        call rmLyAllVars(doGlac,flux_data,flux_meta,nSnow-1,nSnow,nGlce,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        call rmLyAllVars(doGlac,indx_data,indx_meta,nSnow-1,nSnow,nGlce,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if
        ! update the total number of layers
        nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
        nLayers = nSnow + nLake + nSoil + nGlce
        ! save the number of layers
        indx_data%var(iLookINDEX%nSnow)%dat(1)   = nSnow
        indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers
        ! update coordinate variables
        call calcHeight(&
                        ! input/output: data structures
                        indx_data,   & ! intent(in): layer type
                        prog_data,   & ! intent(inout): model variables for a local HRU
                        ! output: error control
                        err,cmessage)
        if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
        ! exit the do loop (no more snow layers to remove)
        return
      else if (doGlac .and. nGlce==1+noThetaChange)then
        err=20; message=trim(message)//'Melted entire water state of glacier, need to start with thicker top layers';return
      end if  ! (special case of 1 layer --> snow without a layer)
       
      ! ***** identify the layer to combine
      if(iLayer==topLayer)then
       jLayer = iLayer+1  ! upper-most layer, combine with its lower neighbor
      elseif(iLayer==botLayer)then
       jLayer = botLayer-1  ! lower-most layer, combine with its upper neighbor
      else
       if(depth(iLayer-1)<depth(iLayer+1))then; jLayer = iLayer-1; else; jLayer = iLayer+1; end if
      end if

      ! ***** combine layers
      ! identify the layer closest to the surface
      kLayer=min(iLayer,jLayer)
      ! combine layer with identified neighbor
      call layer_combine(doGlac,mpar_data,prog_data,diag_data,flux_data,indx_data,kLayer,err,cmessage)
      if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if

      ! update the number of snow layers
      nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
      nLake   = indx_data%var(iLookINDEX%nLake)%dat(1)
      nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
      nGlce   = indx_data%var(iLookINDEX%nGlce)%dat(1)
      nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)
      if(doGlac)then
        botLayer=nSnow+nLake+nSoil+nGlce-noThetaChange
      else
        botLayer=nSnow
      end if

      ! exit the loop to try again
      exit
    end if  ! (if layer is below the mass threshold)
    kLayer=iLayer ! ksnow is used for completion test, so include here
  end do ! (looping through snow layers)

  ! exit if finished
  if(kLayer==nCheck)exit
 end do ! continuous do

 ! handle special case of > maxLayers layers in the CLM option
 if(nSnow > maxLayers .and. ix_snowLayers == rulesDependLayerIndex)then
   ! flag that layers were merged
   mergedLayers=.true.
   ! initial check to ensure everything is wonderful in the universe
   if(nSnow /= maxLayers+1)then; err=5; message=trim(message)//'special case of >maxLayers layers: expect only one more'; return; end if
   ! combine maxLayers-th layer with layer below
   call layer_combine(doGlac,mpar_data,prog_data,diag_data,flux_data,indx_data,maxLayers,err,cmessage)
   ! update the number of snow layers
   nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
   nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if
   ! another check
   if(nSnow /= maxLayers)then; err=5; message=trim(message)//'special case of >maxLayers layers: expect to reduced layers to exactly maxLayers'; return; end if
 end if

 ! check that there are no more than maxLayers layers in the CLM option
 if(ix_snowLayers == rulesDependLayerIndex)then
   if(nSnow > maxLayers)then
     message=trim(message)//'expect no more than maxLayers layers when combination/sub-division rules depend on the layer index (CLM option)'
     err=20; return
   end if
 end if

 ! end association to variables in the data structure
 end associate

 contains

  ! accessors to the layer state, valid across the reallocation of the layer vectors (an associate name would be left
  ! pointing at the old allocation once a layer is added or removed)
  function depth(i); integer(i4b),intent(in) :: i; real(rkind) :: depth; depth = prog_data%var(iLookPROG%mLayerDepth)%dat(i);      end function depth
  function temp(i);  integer(i4b),intent(in) :: i; real(rkind) :: temp;  temp  = prog_data%var(iLookPROG%mLayerTemp)%dat(i);       end function temp
  function liq(i);   integer(i4b),intent(in) :: i; real(rkind) :: liq;   liq   = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(i); end function liq
  function ice(i);   integer(i4b),intent(in) :: i; real(rkind) :: ice;   ice   = prog_data%var(iLookPROG%mLayerVolFracIce)%dat(i); end function ice
  subroutine setDepth(i,x); integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerDepth)%dat(i)      = x; end subroutine setDepth
  subroutine setTemp(i,x);  integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerTemp)%dat(i)       = x; end subroutine setTemp
  subroutine setLiq(i,x);   integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(i) = x; end subroutine setLiq
  subroutine setIce(i,x);   integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerVolFracIce)%dat(i) = x; end subroutine setIce

 end subroutine layerMerge


 ! ***********************************************************************************************************
 ! private subroutine layer_combine: combine snow layers and re-compute model state variables
 ! ***********************************************************************************************************
 ! combines layer iLayer with iLayer+1
 ! ***********************************************************************************************************
 subroutine layer_combine(doGlac,mpar_data,prog_data,diag_data,flux_data,indx_data,iLayer,err,message)
 ! provide access to variables in the data structures
 USE var_lookup,only:iLookPARAM,iLookPROG,iLookINDEX              ! named variables for structure elements
 USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta      ! metadata
 USE data_types,only:var_ilength,var_dlength                      ! data vectors with variable length dimension
 ! provide access to external modules
 USE snow_utils_module,only:fracliquid                                         ! compute fraction of liquid water
 USE convertEnthalpyTemp_module,only:enthalpy2T_snLaGlWat,T2enthalpy_snLaGlWat ! convert temperature to liq+ice enthalpy for a snow/lake/glce layer
 implicit none
 ! ------------------------------------------------------------------------------------------------------------
 ! input/output: data structures
 logical(lgt),intent(in)         :: doGlac    ! flag to denote that merging glacier ice
 type(var_dlength),intent(in)    :: mpar_data ! model parameters
 type(var_dlength),intent(inout) :: prog_data ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data ! model flux variables
 type(var_ilength),intent(inout) :: indx_data ! type of model layer
 ! input: snow layer indices
 integer(i4b),intent(in)         :: iLayer     ! index of top layer to combine
 ! output: error control
 integer(i4b),intent(out)        :: err       ! error code
 character(*),intent(out)        :: message   ! error message
 ! ------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=256)              :: cmessage                 ! error message for downwind routine
 real(rkind)                     :: massIce(2)               ! mass of ice in the two layers identified for combination (kg m-2)
 real(rkind)                     :: massLiq(2)               ! mass of liquid water in the two layers identified for combination (kg m-2)
 real(rkind)                     :: bulkDenWat(2)            ! bulk density if total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(rkind)                     :: cBulkDenWat              ! combined bulk density of total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(rkind)                     :: cTemp                    ! combined layer temperature
 real(rkind)                     :: cDepth                   ! combined layer depth
 real(rkind)                     :: cVolFracIce              ! combined layer volumetric fraction of ice
 real(rkind)                     :: cVolFracLiq              ! combined layer volumetric fraction of liquid water
 real(rkind)                     :: l1Enthalpy,l2Enthalpy    ! enthalpy in the two layers identified for combination (J m-3)
 real(rkind)                     :: cEnthalpy                ! combined layer enthalpy (J m-3)
 real(rkind)                     :: fLiq                     ! fraction of liquid water at the combined temperature cTemp
 real(rkind)                     :: heatCap(2)               ! heat capacity of the ice and liquid water in the two layers (J m-2 K-1)
 real(rkind),parameter           :: eTol=1.e-1_rkind         ! tolerance for the enthalpy-->temperature conversion (J m-3)
 integer(i4b)                    :: nSnow                    ! number of snow layers
 integer(i4b)                    :: nLake                    ! number of lake layers
 integer(i4b)                    :: nSoil                    ! number of soil layers
 integer(i4b)                    :: nGlce                    ! number of glacier ice layers
 integer(i4b)                    :: nLayers                  ! total number of layers
 real(rkind)                     :: frz_scale_use            ! scaling parameter for the snow or glce freezing curve (K-1)
 integer(i4b)                    :: k                        ! index of the two layers to combine

 ! initialize error control
 err=0; message="layer_combine/"

 ! associate local variables with information in the data structures
 ! NOTE: the layer state is read and written through the accessors, since the layer vectors are reallocated when the layer is removed
 associate(&
 ! model parameters
 snowfrz_scale    => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)  & ! scaling parameter for the freezing curve for snow (K-1)
 ) ! (association of local variables with information in the data structures)

 ! initialize the number of layers
 nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
 nLake   = indx_data%var(iLookINDEX%nLake)%dat(1)
 nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
 nGlce   = indx_data%var(iLookINDEX%nGlce)%dat(1)
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)

 if(doGlac)then
  frz_scale_use = snowfrz_scale*icefrz_mult
 else
  frz_scale_use = snowfrz_scale
 end if

 ! compute combined depth
 cDepth = depth(iLayer) + depth(iLayer+1)

 ! compute mass of each layer (kg m-2), and the bulk density of water (kg m-3)
 do k=1,2
   massIce(k)    = iden_ice*ice(iLayer+k-1)*depth(iLayer+k-1)
   massLiq(k)    = iden_water*liq(iLayer+k-1)*depth(iLayer+k-1)
   bulkDenWat(k) = (massIce(k) + massLiq(k))/depth(iLayer+k-1)
 end do
 cBulkDenWat     = (depth(iLayer)*bulkDenWat(1) + depth(iLayer+1)*bulkDenWat(2))/cDepth

 ! compute enthalpy for each layer (J m-3)
 l1Enthalpy = T2enthalpy_snLaGlWat(temp(iLayer),  bulkDenWat(1),frz_scale_use)
 l2Enthalpy = T2enthalpy_snLaGlWat(temp(iLayer+1),bulkDenWat(2),frz_scale_use)

 ! compute combined enthalpy (J m-3)
 cEnthalpy = (depth(iLayer)*l1Enthalpy + depth(iLayer+1)*l2Enthalpy)/cDepth

 ! convert enthalpy (J m-3) to temperature (K)
 call enthalpy2T_snLaGlWat(cEnthalpy,cBulkDenWat,frz_scale_use,cTemp,.not.doGlac,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if

 ! test enthalpy conversion
 if(abs(T2enthalpy_snLaGlWat(cTemp,cBulkDenWat,frz_scale_use)/cBulkDenWat - cEnthalpy/cBulkDenWat) > eTol)then
  write(*,'(a,1x,f12.5,1x,2(e20.10,1x))') 'enthalpy test', cBulkDenWat, T2enthalpy_snLaGlWat(cTemp,cBulkDenWat,frz_scale_use)/cBulkDenWat, cEnthalpy/cBulkDenWat
  message=trim(message)//'problem with enthalpy-->temperature conversion'
  err=20; return
 end if

 ! check temperature is within the two temperatures
 ! NOTE: use tolerance, for cases of merging a layer that has just been split
 if(cTemp > max(temp(iLayer),temp(iLayer+1))+eTol)then; err=20; message=trim(message)//'merged temperature > max(temp1,temp2)'; return; end if
 if(cTemp < min(temp(iLayer),temp(iLayer+1))-eTol)then; err=20; message=trim(message)//'merged temperature < min(temp1,temp2)'; return; end if

 ! compute volumetric fraction of liquid water
 fLiq = fracliquid(cTemp,frz_scale_use)

 ! compute volumetric fraction of ice and liquid water
 cVolFracLiq =          fLiq *cBulkDenWat/iden_water
 cVolFracIce = (1._rkind - fLiq)*cBulkDenWat/iden_ice

 ! glacier ice: keep the ice and liquid water as they are, with the temperature that conserves the sensible heat (exact for
 ! a fixed composition). The liquid in melting ice is held at its residual content by drainage, and the ice freezing curve
 ! is so steep there that the liquid re-derived above from the merged temperature turns the temperature tolerance of the
 ! solver into a change of ice and liquid at every merge (the enthalpy conversion above is kept for its checks).
 if(doGlac)then
  do k=1,2
    heatCap(k) = Cp_ice*massIce(k) + Cp_water*massLiq(k)
  end do
  cTemp       = (heatCap(1)*temp(iLayer) + heatCap(2)*temp(iLayer+1))/(heatCap(1) + heatCap(2))
  cVolFracIce = (massIce(1) + massIce(2))/(iden_ice*cDepth)
  cVolFracLiq = (massLiq(1) + massLiq(2))/(iden_water*cDepth)
 end if

 ! end association of local variables with information in the data structures
 end associate

 ! remove a model layer from all model variable vectors
 call rmLyAllVars(doGlac,prog_data,prog_meta,iLayer,nSnow,nGlce,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 call rmLyAllVars(doGlac,diag_data,diag_meta,iLayer,nSnow,nGlce,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 call rmLyAllVars(doGlac,flux_data,flux_meta,iLayer,nSnow,nGlce,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 call rmLyAllVars(doGlac,indx_data,indx_meta,iLayer,nSnow,nGlce,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! define the combined layer as snow/glacier ice
 if (nSnow>0)then
   indx_data%var(iLookINDEX%layerType)%dat(iLayer) = iname_snow
 else
   indx_data%var(iLookINDEX%layerType)%dat(iLayer) = iname_glce
 end if

 ! save the number of layers in the data structures
 indx_data%var(iLookINDEX%nSnow)%dat(1)   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
 indx_data%var(iLookINDEX%nLake)%dat(1)   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_lake)
 indx_data%var(iLookINDEX%nSoil)%dat(1)   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)
 indx_data%var(iLookINDEX%nGlce)%dat(1)   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_glce)
 indx_data%var(iLookINDEX%nLayers)%dat(1) = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nSoil)%dat(1) &
                                          + indx_data%var(iLookINDEX%nGlce)%dat(1) + indx_data%var(iLookINDEX%nLake)%dat(1)

 ! update the number of \layers
 nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
 nLake   = indx_data%var(iLookINDEX%nLake)%dat(1)
 nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
 nGlce   = indx_data%var(iLookINDEX%nGlce)%dat(1)
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)

 ! ***** put state variables for the combined layer in the appropriate place
 call setTemp(iLayer, cTemp)
 call setDepth(iLayer, cDepth)
 call setIce(iLayer, cVolFracIce)
 call setLiq(iLayer, cVolFracLiq)

 ! ***** adjust coordinate variables
 call calcHeight(&
                 ! input/output: data structures
                 indx_data,   & ! intent(in): layer type
                 prog_data,   & ! intent(inout): model variables for a local HRU
                 ! output: error control
                 err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 contains

  ! accessors to the layer state, valid across the reallocation of the layer vectors (an associate name would be left
  ! pointing at the old allocation once a layer is added or removed)
  function depth(i); integer(i4b),intent(in) :: i; real(rkind) :: depth; depth = prog_data%var(iLookPROG%mLayerDepth)%dat(i);      end function depth
  function temp(i);  integer(i4b),intent(in) :: i; real(rkind) :: temp;  temp  = prog_data%var(iLookPROG%mLayerTemp)%dat(i);       end function temp
  function liq(i);   integer(i4b),intent(in) :: i; real(rkind) :: liq;   liq   = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(i); end function liq
  function ice(i);   integer(i4b),intent(in) :: i; real(rkind) :: ice;   ice   = prog_data%var(iLookPROG%mLayerVolFracIce)%dat(i); end function ice
  subroutine setDepth(i,x); integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerDepth)%dat(i)      = x; end subroutine setDepth
  subroutine setTemp(i,x);  integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerTemp)%dat(i)       = x; end subroutine setTemp
  subroutine setLiq(i,x);   integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(i) = x; end subroutine setLiq
  subroutine setIce(i,x);   integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerVolFracIce)%dat(i) = x; end subroutine setIce

 end subroutine layer_combine


 ! ***********************************************************************************************************
 ! private subroutine rmLyAllVars: reduce the length of the vectors in data structures
 ! ***********************************************************************************************************
 ! removes layer "iLayer+1" and sets layer "iLayer" to a missing value
 ! (layer "iLayer" will be filled with a combined layer later)
 ! iLayer indexes the whole column; the snow and glacier domain vectors start at the top of their domain
 ! ***********************************************************************************************************
 subroutine rmLyAllVars(doGlac,dataStruct,metaStruct,iLayer,nSnow,nGlce,nLayers,err,message)
 USE var_lookup,only:iLookVarType                 ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 USE f2008_funcs_module,only:cloneStruc           ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength      ! data vectors with variable length dimension
 USE data_types,only:var_info                     ! metadata structure
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 logical(lgt),intent(in)         :: doGlac         ! flag to denote that merging glacier ice
 class(*),intent(inout)          :: dataStruct     ! data structure
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in)         :: iLayer          ! new layer (index in the whole column)
 integer(i4b),intent(in)         :: nSnow,nGlce,nLayers ! number of snow, glacier ice layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! locals
 integer(i4b)                    :: iVar           ! variable index
 integer(i4b)                    :: ix_lower       ! lower bound of the vector
 integer(i4b)                    :: ix_upper       ! upper bound of the vector
 integer(i4b)                    :: ix_rm          ! index of the new layer within the vector
 real(rkind),allocatable         :: tempVec_rkind(:)  ! temporary vector (double precision)
 integer(i4b),allocatable        :: tempVec_i4b(:) ! temporary vector (integer)
 character(LEN=256)              :: cmessage       ! error message of downwind routine
 ! initialize error control
 err=0; message="rmLyAllVars/"

 ! check dimensions
 select type(dataStruct)
  type is (var_dlength); if(size(dataStruct%var) /= size(metaStruct)) err=20
  type is (var_ilength); if(size(dataStruct%var) /= size(metaStruct)) err=20
  class default; err=20; message=trim(message)//'unable to identify the data type'; return
 end select
 if(err/=0)then; message=trim(message)//'dimensions of data structure and metadata structures do not match'; return; end if

 ! ***** loop through model variables and remove one layer
 do iVar=1,size(metaStruct)

  ! define bounds, and the index of the new layer within the vector (the glacier domain vectors start below the layers above the glacier)
  if (doGlac)then
   select case(metaStruct(iVar)%varType)
    case(iLookVarType%midGlce); ix_lower=1; ix_upper=nGlce;   ix_rm=iLayer-(nLayers-nGlce)
    case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers; ix_rm=iLayer
    case(iLookVarType%ifcGlce); ix_lower=0; ix_upper=nGlce;   ix_rm=iLayer-(nLayers-nGlce)
    case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers; ix_rm=iLayer
    case default; cycle  ! no need to remove soil layers or scalar variables
   end select
  else
    select case(metaStruct(iVar)%varType)
     case(iLookVarType%midSnow); ix_lower=1; ix_upper=nSnow;   ix_rm=iLayer
     case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers; ix_rm=iLayer
     case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=nSnow;   ix_rm=iLayer
     case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers; ix_rm=iLayer
     case default; cycle  ! no need to remove soil layers or scalar variables
    end select
   end if

  ! remove layers
  select type(dataStruct)

   ! ** double precision
   type is (var_dlength)
    ! check allocated
    if(.not.allocated(dataStruct%var(iVar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! allocate the temporary vector
    allocate(tempVec_rkind(ix_lower:ix_upper-1), stat=err)
    if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
    ! copy elements across to the temporary vector
    if(ix_rm>=ix_lower)  tempVec_rkind(ix_rm)              = realMissing ! set merged layer to missing (fill in later)
    if(ix_rm>ix_lower)   tempVec_rkind(ix_lower:ix_rm-1)   = dataStruct%var(iVar)%dat(ix_lower:ix_rm-1)
    if(ix_rm+1<ix_upper) tempVec_rkind(ix_rm+1:ix_upper-1) = dataStruct%var(iVar)%dat(ix_rm+2:ix_upper)  ! skip ix_rm+1
    ! deallocate the data vector: strictly not necessary, but include to be safe
    deallocate(dataStruct%var(iVar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem deallocating data vector'; return; end if
    ! create the new data structure using the temporary vector as the source
    call cloneStruc(dataStruct%var(iVar)%dat, ix_lower, source=tempVec_rkind, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! deallocate the temporary data vector: strictly not necessary, but include to be safe
    deallocate(tempVec_rkind,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! ** integer
   type is (var_ilength)
    ! check allocated
    if(.not.allocated(dataStruct%var(iVar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! allocate the temporary vector
    allocate(tempVec_i4b(ix_lower:ix_upper-1), stat=err)
    if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
    ! copy elements across to the temporary vector
    if(ix_rm>=ix_lower)  tempVec_i4b(ix_rm)              = integerMissing ! set merged layer to missing (fill in later)
    if(ix_rm>ix_lower)   tempVec_i4b(ix_lower:ix_rm-1)   = dataStruct%var(iVar)%dat(ix_lower:ix_rm-1)
    if(ix_rm+1<ix_upper) tempVec_i4b(ix_rm+1:ix_upper-1) = dataStruct%var(iVar)%dat(ix_rm+2:ix_upper)  ! skip ix_rm+1
    ! deallocate the data vector: strictly not necessary, but include to be safe
    deallocate(dataStruct%var(iVar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem deallocating data vector'; return; end if
    ! create the new data structure using the temporary vector as the source
    call cloneStruc(dataStruct%var(iVar)%dat, ix_lower, source=tempVec_i4b, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! deallocate the temporary data vector: strictly not necessary, but include to be safe
    deallocate(tempVec_i4b,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! check that we found the data type
   class default; err=20; message=trim(message)//'unable to identify the data type'; return

  end select ! dependence on data types

 end do  ! looping through variables

 end subroutine rmLyAllVars

end module layerMerge_module
