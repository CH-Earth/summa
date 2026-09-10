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

module heat_Cp_Cm_module

! data types
USE nr_type

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength         ! data vector with variable length dimension (rkind)

! named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookDIAG,iLookINDEX,iLookPROG  ! named variables for structure elements

! physical constants
USE multiconst,only: gravity, &                          ! gravitational acceleration (m s-1)
                     Tfreeze, &                          ! freezing point of water (K)
                     Cp_soil,Cp_water,Cp_ice,Cp_air,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                     iden_water,iden_ice,iden_air,&      ! intrinsic density of water and ice (kg m-3)
                     LH_fus                              ! latent heat of fusion (J kg-1)

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair  ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy  ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy  ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy  ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer   ! named variable defining the energy state variable for layers
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for layers
USE globalData,only:iname_liqLayer   ! named variable defining the liquid  water state variable for layers
USE globalData,only:iname_matLayer   ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer   ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:iname_watAquifer ! named variable defining the water storage in the aquifer

! missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing real

USE globalData,only:icefrz_mult      ! freezing curve scaling factor multipier of snow to ice, closer to a step function since ice does not hold water

! domain types
USE globalData,only:iname_cas        ! named variables for canopy air space
USE globalData,only:iname_veg        ! named variables for vegetation canopy
USE globalData,only:iname_snow       ! named variables for snow
USE globalData,only:iname_soil       ! named variables for soil
USE globalData,only:iname_glce       ! named variables for glacier ice
USE globalData,only:iname_lake       ! named variables for lake
USE globalData,only:iname_aquifer    ! named variables for the aquifer

! privacy
implicit none
private
public::stateMultiplier
public::init_heatCapacity
public::heatCapacity
public::heatAdvectWat

contains


! **********************************************************************************************************
! public subroutine stateMultiplier: get scale factors for the temperature and water state vector
! **********************************************************************************************************
subroutine stateMultiplier(&
                      heatCapVeg,              & ! intent(in):  heat capacity for canopy
                      mLayerHeatCap,           & ! intent(in):  heat capacity for layers
                      ! input: data structures
                      indx_data,               & ! intent(in):  indices defining model states and layers
                      ! output
                      sMul,                    & ! intent(out): multiplier for state vector (used in the residual calculations)
                      err,message)               ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: data structures
  real(qp),intent(in)             :: heatCapVeg             ! volumetric heat capacity of vegetation (J m-3 K-1)
  real(qp),intent(in)             :: mLayerHeatCap(:)       ! volumetric heat capacity of layers (J m-3 K-1)
  type(var_ilength),intent(in)    :: indx_data              ! indices defining model states and layers
  ! output: state vectors
  real(qp),intent(inout)          :: sMul(:)    ! NOTE: qp  ! multiplier for state vector (used in the residual calculations)
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! state subsets
  integer(i4b)                    :: iLayer                 ! index of layer within the layer domains
  integer(i4b)                    :: ixStateSubset          ! index within the state subset
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  associate(&
    ! vector of energy and hydrology indices for the layer domains
    ixSnLaSoGlNrg       => indx_data%var(iLookINDEX%ixSnLaSoGlNrg)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the layer domains
    ixSnLaSoGlHyd       => indx_data%var(iLookINDEX%ixSnLaSoGlHyd)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the layer domains
    nSnLaSoGlNrg        => indx_data%var(iLookINDEX%nSnLaSoGlNrg )%dat(1)         ,& ! intent(in) : [i4b]    number of energy state variables in the layer domains
    nSnLaSoGlHyd        => indx_data%var(iLookINDEX%nSnLaSoGlHyd )%dat(1)         ,& ! intent(in) : [i4b]    number of hydrology state variables in the layer domains
    ! type of model state variabless
    ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
    ! number of layers
    nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in) : [i4b]    total number of layers
    )  ! end association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='stateMultiplier/'

    ! -----
    ! * define components of derivative matrices at start of time step (substep)...
    ! ------------------------------------------------------------------------------------------

    ! define the multiplier for the state vector for residual calculations (vegetation canopy)
    ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
    where(ixStateType_subset==iname_nrgCanair) sMul = Cp_air*iden_air ! volumetric heat capacity of air (J m-3 K-1)
    where(ixStateType_subset==iname_nrgCanopy) sMul = heatCapVeg      ! volumetric heat capacity of the vegetation (J m-3 K-1)
    where(ixStateType_subset==iname_watCanopy) sMul = 1._rkind        ! nothing else on the left hand side
    where(ixStateType_subset==iname_liqCanopy) sMul = 1._rkind        ! nothing else on the left hand side

    ! define the energy multiplier for the state vector for residual calculations (layer domains)
    if(nSnLaSoGlNrg>0)then
      do concurrent (iLayer=1:nLayers,ixSnLaSoGlNrg(iLayer)/=integerMissing) ! (loop through non-missing energy state variables in the layer domains)
        ixStateSubset        = ixSnLaSoGlNrg(iLayer) ! index within the state vector
        sMul(ixStateSubset)  = mLayerHeatCap(iLayer) ! transfer volumetric heat capacity to the state multiplier
      end do  ! looping through non-missing energy state variables in the layer domains
    endif

    ! define the hydrology multiplier and diagonal elements for the state vector for residual calculations (layer domains)
    if(nSnLaSoGlHyd>0)then
      do concurrent (iLayer=1:nLayers,ixSnLaSoGlHyd(iLayer)/=integerMissing) ! (loop through non-missing energy state variables in the layer domains)
        ixStateSubset        = ixSnLaSoGlHyd(iLayer) ! index within the state vector
        sMul(ixStateSubset)  = 1._rkind              ! state multiplier = 1 (nothing else on the left-hand-side)
      end do  ! looping through non-missing energy state variables in the layer domains
    endif

    ! define the scaling factor and diagonal elements for the aquifer
    where(ixStateType_subset==iname_watAquifer)  sMul = 1._rkind

  end associate
! end association to variables in the data structure where vector length does not change
end subroutine stateMultiplier

 ! **********************************************************************************************************
 ! public subroutine init_heatCapacity: compute start-of-step heat capacity
 ! **********************************************************************************************************
 subroutine init_heatCapacity(&
                       ! input: control variables
                       computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                       canopyDepth,             & ! intent(in):    canopy depth (m)
                       ! input/output: data structures
                       mpar_data,               & ! intent(in):    model parameters
                       indx_data,               & ! intent(in):    model layer indices
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       ! output: error control
                       err,message)               ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)         :: computeVegFlux         ! logical flag to denote if computing the vegetation flux
 real(rkind),intent(in)          :: canopyDepth            ! depth of the vegetation canopy (m)
 ! input/output: data structures
 type(var_dlength),intent(in)    :: mpar_data              ! model parameters
 type(var_ilength),intent(in)    :: indx_data              ! model layer indices
 type(var_dlength),intent(in)    :: prog_data              ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data              ! model diagnostic variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iLayer                 ! index of model layer
 integer(i4b)                    :: iSoil                  ! index of soil layer
  ! --------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! input: state variables
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),           & ! intent(in): canopy ice content (kg m-2)
 scalarCanopyLiquid      => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),           & ! intent(in): canopy liquid water content (kg m-2)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat,             & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat,             & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
 ! input: coordinate variables
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),                    & ! intent(in): number of snow layers
 nLake                   => indx_data%var(iLookINDEX%nLake)%dat(1),                    & ! intent(in): number of lake layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1),                  & ! intent(in): total number of layers
 layerType               => indx_data%var(iLookINDEX%layerType)%dat,                   & ! intent(in): layer type (iname_soil or iname_snow)
 ! input: heat capacity
 specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),          & ! intent(in): specific heat of vegetation (J kg-1 K-1)
 maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),        & ! intent(in): maximum mass of vegetation (kg m-2)
 ! input: depth varying soil parameters
 iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat,              & ! intent(in): intrinsic density of soil (kg m-3)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                   & ! intent(in): soil porosity (-)
 ! output: diagnostic variables
 scalarBulkVolHeatCapVeg => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),   & ! intent(out): volumetric heat capacity of the vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat            & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1)
 )  ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="init_heatCapacity/"

 ! initialize the soil layer
 iSoil=integerMissing

 ! compute the bulk volumetric heat capacity of vegetation (J m-3 K-1)
 if(computeVegFlux)then
  scalarBulkVolHeatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                            Cp_water*scalarCanopyLiquid/canopyDepth       + & ! liquid water component
                            Cp_ice*scalarCanopyIce/canopyDepth                ! ice component
 else
  scalarBulkVolHeatCapVeg = realMissing
 end if

 ! loop through layers
 do iLayer=1,nLayers
  ! get the soil layer
  if(iLayer>nSnow+nLake) iSoil = iLayer-nSnow-nLake

  select case(layerType(iLayer))
   ! * soil
   case(iname_soil)
    mLayerVolHtCapBulk(iLayer) = iden_soil(iSoil)  * Cp_soil  * ( 1._rkind - theta_sat(iSoil) ) + & ! soil component
                                 iden_ice          * Cp_ice   * mLayerVolFracIce(iLayer)        + & ! ice component
                                 iden_water        * Cp_water * mLayerVolFracLiq(iLayer)        + & ! liquid water component
                                 iden_air          * Cp_air   * ( theta_sat(iSoil) - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)) ) ! air component
   ! * snow, lake, and glacier ice
   case(iname_snow, iname_lake, iname_glce)
    mLayerVolHtCapBulk(iLayer) = iden_ice          * Cp_ice   * mLayerVolFracIce(iLayer)     + & ! ice component
                                 iden_water        * Cp_water * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                                 iden_air          * Cp_air   * ( 1._rkind - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)) ) ! air component
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute volumetric heat capacity'; return
  end select
 end do  ! looping through layers

 end associate
end subroutine init_heatCapacity

! **********************************************************************************************************
! public subroutine heatCapacity: compute diagnostic energy variable Cp (change in enthTemp with temperature)
! **********************************************************************************************************
subroutine heatCapacity(&
                      ! input: state variables
                      canopyDepth,             & ! intent(in):    canopy depth (m)
                      scalarCanopyIce,         & ! intent(in):    trial value for mass of ice on the vegetation canopy (kg m-2)
                      scalarCanopyLiquid,      & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                      scalarCanopyTemp,        & ! intent(in):    trial value of canopy temperature (K)
                      mLayerVolFracIce,        & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                      mLayerVolFracLiq,        & ! intent(in):    volumetric fraction of liquid water at the start of the sub-step (-)
                      mLayerTemp,              & ! intent(in):    trial value of layer temperature (K)
                      mLayerMatricHead,        & ! intent(in):    total water matric potential (m)
                      ! input: pre-computed derivatives
                      dTheta_dTkCanopy,        & ! intent(in):    derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
                      scalarFracLiqVeg,        & ! intent(in):    fraction of canopy liquid water (-)
                      mLayerdTheta_dTk,        & ! intent(in):    derivative of volumetric liquid water content w.r.t. temperature (K-1)
                      mLayerFracLiq,           & ! intent(in):    fraction of liquid water (-)
                      dVolTot_dPsi0,           & ! intent(in):    derivative in total water content w.r.t. total water matric potential (m-1)
                      ! input output data structures
                      mpar_data,               & ! intent(in):    model parameters
                      indx_data,               & ! intent(in):    model layer indices
                      ! output
                      heatCapVeg,              & ! intent(inout): heat capacity for canopy
                      mLayerHeatCap,           & ! intent(inout): heat capacity for snow and soil
                      dVolHtCapBulk_dPsi0,     & ! intent(inout): derivative in bulk heat capacity w.r.t. matric potential
                      dVolHtCapBulk_dTheta,    & ! intent(inout): derivative in bulk heat capacity w.r.t. volumetric water content
                      dVolHtCapBulk_dCanWat,   & ! intent(inout): derivative in bulk heat capacity w.r.t. volumetric water content
                      dVolHtCapBulk_dTk,       & ! intent(inout): derivative in bulk heat capacity w.r.t. temperature
                      dVolHtCapBulk_dTkCanopy, & ! intent(inout): derivative in bulk heat capacity w.r.t. temperature     
                      ! output: error control
                      err,message)               ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! provide access to external subroutines
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: canopyDepth             ! depth of the vegetation canopy (m)
  real(rkind),intent(in)          :: scalarCanopyIce         ! trial value of canopy ice content (kg m-2)
  real(rkind),intent(in)          :: scalarCanopyLiquid      ! trial value of canopy liquid content (kg m-2)
  real(rkind),intent(in)          :: scalarCanopyTemp        ! value of canopy temperature (kg m-2)
  real(rkind),intent(in)          :: mLayerVolFracLiq(:)     ! trial vector of volumetric liquid water content (-)
  real(rkind),intent(in)          :: mLayerVolFracIce(:)     ! trial vector of volumetric ice water content (-)
  real(rkind),intent(in)          :: mLayerTemp(:)           ! vector of temperature (-)
  real(rkind),intent(in)          :: mLayerMatricHead(:)     ! vector of total water matric potential (m)
  ! input: pre-computed derivatives 
  real(rkind),intent(in)          :: dTheta_dTkCanopy        ! derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)          :: scalarFracLiqVeg        ! fraction of canopy liquid water (-)
  real(rkind),intent(in)          :: mLayerdTheta_dTk(:)     ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)          :: mLayerFracLiq(:)        ! fraction of liquid water (-)
  real(rkind),intent(in)          :: dVolTot_dPsi0(:)        ! derivative in total water content w.r.t. total water matric potential (m-1)
  ! input/output: data structures 
  type(var_dlength),intent(in)    :: mpar_data               ! model parameters
  type(var_ilength),intent(in)    :: indx_data               ! model layer indices
  ! output 
  real(qp),intent(inout)          :: heatCapVeg              ! heat capacity for canopy
  real(qp),intent(inout)          :: mLayerHeatCap(:)        ! heat capacity for snow and soil
  real(rkind),intent(inout)       :: dVolHtCapBulk_dPsi0(:)  ! derivative in bulk heat capacity w.r.t. matric potential
  real(rkind),intent(inout)       :: dVolHtCapBulk_dTheta(:) ! derivative in bulk heat capacity w.r.t. volumetric water content
  real(rkind),intent(inout)       :: dVolHtCapBulk_dCanWat   ! derivative in bulk heat capacity w.r.t. volumetric water content
  real(rkind),intent(inout)       :: dVolHtCapBulk_dTk(:)    ! derivative in bulk heat capacity w.r.t. temperature
  real(rkind),intent(inout)       :: dVolHtCapBulk_dTkCanopy ! derivative in bulk heat capacity w.r.t. temperature
   ! output: error control 
  integer(i4b),intent(out)        :: err                     ! error code
  character(*),intent(out)        :: message                 ! error message
  ! -------------------------------------------------------- ------------------------------------------------------------------------
  ! local variables 
  integer(i4b)                    :: iState                  ! index of model state variable
  integer(i4b)                    :: iLayer                  ! index of model layer
  integer(i4b)                    :: ixFullVector            ! index within full state vector
  integer(i4b)                    :: ixDomainType            ! name of a given model domain
  integer(i4b)                    :: ixControlIndex          ! index within a given model domain
  real(rkind)                     :: fLiq                    ! fraction of liquid water
  real(rkind)                     :: Tcrit                   ! temperature where all water is unfrozen (K)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! associate variables in data structure
  associate(&
    ! input: coordinate variables
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)             ,& ! intent(in): number of snow layers
    nLake                   => indx_data%var(iLookINDEX%nLake)%dat(1)             ,& ! intent(in): number of lake layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)             ,& ! intent(in): number of soil layers
    ! mapping between the full state vector and the state subset
    ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat     ,& ! intent(in): [state subset] list of indices of the full state vector in the state subset
    ! type of domain, type of state variable, and index of control volume within domain
    ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat  ,& ! intent(in): [state subset] id of domain for desired model state variables
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat      ,& ! intent(in): index of the control volume for different domains (veg, snow, soil)
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat          ,& ! intent(in): indices defining the type of the state (iname_nrgLayer...)
    ! input: heat capacity
    specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1)   ,& ! intent(in): specific heat of vegetation (J kg-1 K-1)
    maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1) ,& ! intent(in): maximum mass of vegetation (kg m-2)
    ! input: depth varying soil parameters
    iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat       ,& ! intent(in): intrinsic density of soil (kg m-3)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat             & ! intent(in): soil porosity (-)
    )  ! end associate statement
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="heatCapacity/"

    ! loop through model state variables
    do iState=1,size(ixMapSubset2Full)

      ! -----
      ! - compute indices...
      ! --------------------

      ! get domain type, and index of the control volume within the domain
      ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
      ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
      ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

      ! check an energy state, since only need for energy state equations
      if(ixStateType(ixFullVector)==iname_nrgCanair .or. ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)then

        ! get the layer index
        select case(ixDomainType)
          case(iname_cas);     cycle ! canopy air space, do nothing (no water stored in canopy air space)
          case(iname_veg);     iLayer = integerMissing
          case(iname_snow);    iLayer = ixControlIndex
          case(iname_lake);    iLayer = ixControlIndex + nSnow
          case(iname_soil);    iLayer = ixControlIndex + nSnow + nLake
          case(iname_glce);    iLayer = ixControlIndex + nSnow + nLake + nSoil
          case(iname_aquifer); cycle ! aquifer: do nothing (no thermodynamics in the aquifer)
          case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_lake, iname_soil, iname_glce, iname_soil, or iname_aquifer'; return
        end select

        ! identify domain
        select case(ixDomainType)

          case(iname_veg)
            heatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                         Cp_water*scalarCanopyLiquid/canopyDepth       + & ! liquid water component
                         Cp_ice*scalarCanopyIce/canopyDepth                ! ice component

            ! derivatives
            fLiq = scalarFracLiqVeg
            dVolHtCapBulk_dCanWat = ( -Cp_ice*( fLiq-1._rkind ) + Cp_water*fLiq )/canopyDepth !this is iden_water/(iden_water*canopyDepth)
            if(scalarCanopyTemp < Tfreeze)then
              dVolHtCapBulk_dTkCanopy = iden_water * (-Cp_ice + Cp_water) * dTheta_dTkCanopy ! no derivative in air
            else
              dVolHtCapBulk_dTkCanopy = 0._rkind
            endif

          case(iname_snow, iname_lake, iname_glce)
            mLayerHeatCap(iLayer) =  iden_ice   * Cp_ice   * mLayerVolFracIce(iLayer) + & ! ice component
                                     iden_water * Cp_water * mLayerVolFracLiq(iLayer) + & ! liquid water component
                                     iden_air   * Cp_air   * ( 1._rkind - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)) ) ! air component
            ! derivatives
            fLiq = mLayerFracLiq(iLayer)
            dVolHtCapBulk_dTheta(iLayer) = iden_water * ( -Cp_ice*( fLiq-1._rkind ) + Cp_water*fLiq ) + iden_air * ( ( fLiq-1._rkind )*iden_water/iden_ice - fLiq ) * Cp_air
            if( mLayerTemp(iLayer) < Tfreeze)then
              dVolHtCapBulk_dTk(iLayer) = ( iden_water * (-Cp_ice + Cp_water) + iden_air * (iden_water/iden_ice - 1._rkind) * Cp_air ) * mLayerdTheta_dTk(iLayer)
            else
              dVolHtCapBulk_dTk(iLayer) = 0._rkind
            endif

          case(iname_soil)
            mLayerHeatCap(iLayer) =  iden_soil(ixControlIndex) * Cp_soil  * ( 1._rkind - theta_sat(ixControlIndex) ) + & ! soil component
                                     iden_ice                  * Cp_ice   * mLayerVolFracIce(iLayer)                 + & ! ice component
                                     iden_water                * Cp_water * mLayerVolFracLiq(iLayer)                 + & ! liquid water component
                                     iden_air                  * Cp_air   * ( theta_sat(ixControlIndex) - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)) )! air component
           ! derivatives
           dVolHtCapBulk_dTheta(iLayer) = iden_water* Cp_water ! shouldn't need, but compute anyway
           Tcrit = crit_soilT( mLayerMatricHead(ixControlIndex) )
           if( mLayerTemp(iLayer) < Tcrit)then
             dVolHtCapBulk_dPsi0(ixControlIndex) = (iden_ice * Cp_ice   - iden_air * Cp_air) * dVolTot_dPsi0(ixControlIndex)
             dVolHtCapBulk_dTk(iLayer) = (-iden_ice * Cp_ice + iden_water * Cp_water) * mLayerdTheta_dTk(iLayer)
           else
             dVolHtCapBulk_dPsi0(ixControlIndex) = (iden_water*Cp_water - iden_air * Cp_air) * dVolTot_dPsi0(ixControlIndex)
             dVolHtCapBulk_dTk(iLayer) = 0._rkind
           endif
        end select

      end if  ! if an energy layer
    end do  ! looping through state variables

  end associate

end subroutine heatCapacity

! **********************************************************************************************************
! public subroutine heatAdvectWat: compute diagnostic energy variable Cm (change in enthTemp with water)
! **********************************************************************************************************
subroutine heatAdvectWat(&
                      ! input: state variables
                      scalarCanopyTemp,        & ! intent(in):  value of canopy temperature (K)
                      mLayerTemp,              & ! intent(in):  vector of temperature (K)
                      mLayerMatricHead,        & ! intent(in):  vector of total water matric potential (-)
                      ! input data structures
                      mpar_data,               & ! intent(in):  model parameters
                      indx_data,               & ! intent(in):  model layer indices
                      ! output
                      scalarCanopyCm,          & ! intent(inout): Cm for vegetation (J kg K-1)
                      mLayerCm,                & ! intent(inout): Cm for layers (J kg K-1)
                      dCm_dPsi0,               & ! intent(inout): derivative in Cm w.r.t. matric potential (J kg)
                      dCm_dTk,                 & ! intent(inout): derivative in Cm w.r.t. temperature (J kg K-2)
                      dCm_dTkCanopy,           & ! intent(inout): derivative in Cm w.r.t. temperature (J kg K-2)
                      ! output: error control
                      err,message)               ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! provide access to external subroutines
  USE snow_utils_module,only:fracliquid     ! compute the fraction of liquid water (snow)
  USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists (soil)
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! input: state variables
  real(rkind),intent(in)               :: scalarCanopyTemp       ! value of canopy temperature (K)
  real(rkind),intent(in)               :: mLayerTemp(:)          ! vector of temperature (K)
  real(rkind),intent(in)               :: mLayerMatricHead(:)    ! vector of total water matric potential (-)
  ! input/output: data structures
  type(var_dlength),intent(in)         :: mpar_data              ! model parameters
  type(var_ilength),intent(in)         :: indx_data              ! model layer indices
  ! output: Cm and derivatives
  real(rkind),intent(inout)            :: scalarCanopyCm         ! Cm for vegetation (J kg K-1) use for LHS
  real(rkind),intent(inout)            :: mLayerCm(:)            ! Cm for layers (J kg K-1)
  real(rkind),intent(inout)            :: dCm_dPsi0(:)           ! derivative in Cm w.r.t. matric potential (J kg)
  real(rkind),intent(inout)            :: dCm_dTk(:)             ! derivative in Cm w.r.t. temperature (J kg K-2)
  real(rkind),intent(inout)            :: dCm_dTkCanopy          ! derivative in Cm w.r.t. temperature (J kg K-2)
  ! output: error control
  integer(i4b),intent(out)             :: err                    ! error code
  character(*),intent(out)             :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                         :: iState                 ! index of model state variable
  integer(i4b)                         :: iLayer                 ! index of model layer
  integer(i4b)                         :: ixFullVector           ! index within full state vector
  integer(i4b)                         :: ixDomainType           ! name of a given model domain
  integer(i4b)                         :: ixControlIndex         ! index within a given model domain
  real(rkind)                          :: diffT                  ! temperature difference from Tfreeze
  real(rkind)                          :: diff0                  ! temperature difference Tcrit from Tfreeze
  real(rkind)                          :: integral               ! integral of snow freezing curve
  real(rkind)                          :: fLiq                   ! fraction of liquid water
  real(rkind)                          :: dfLiq_dT               ! derivative of fraction of liquid water with temperature
  real(rkind)                          :: Tcrit                  ! temperature where all water is unfrozen (K)
  real(rkind)                          :: dTcrit_dPsi0           ! derivative of critical temperature with matric potential
  real(rkind)                          :: frz_scale_use          ! scaling parameter for the snow or glce freezing curve (K-1)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! associate variables in data structure
  associate(&
    ! input: coordinate variables
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)             ,& ! intent(in): [i4b] number of snow layers
    nLake                   => indx_data%var(iLookINDEX%nLake)%dat(1)             ,& ! intent(in): [i4b] number of lake layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)             ,& ! intent(in): [i4b] number of soil layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)           ,& ! intent(in): [i4b] total number of layers
    snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)     ,& ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
    noThetaChange           => indx_data%var(iLookINDEX%noThetaChange)%dat(1)     ,& ! intent(in): [i4b] number of layers with no change in total water content (bottom layers)
    ! mapping between the full state vector and the state subset
    ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat     ,& ! intent(in): [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! type of domain, type of state variable, and index of control volume within domain
    ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat  ,& ! intent(in): [i4b(:)] [state subset] id of domain for desired model state variables
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat      ,& ! intent(in): [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat           & ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    )  ! end associate statement
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="heatAdvectWat/"

    ! loop through model state variables
    do iState=1,size(ixMapSubset2Full)

      ! -----
      ! - compute indices...
      ! --------------------

      ! get domain type, and index of the control volume within the domain
      ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
      ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
      ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

      ! check an energy state, since only need for energy state equations
      if(ixStateType(ixFullVector)==iname_nrgCanair .or. ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)then

        ! get the layer index
        select case(ixDomainType)
          case(iname_cas);     cycle ! canopy air space, do nothing (no water stored in canopy air space)
          case(iname_veg);     iLayer = integerMissing
          case(iname_snow);    iLayer = ixControlIndex
          case(iname_lake);    iLayer = ixControlIndex + nSnow
          case(iname_soil);    iLayer = ixControlIndex + nSnow + nLake
          case(iname_glce);    iLayer = ixControlIndex + nSnow + nLake + nSoil
          case(iname_aquifer); cycle ! aquifer: do nothing (no thermodynamics in the aquifer)
          case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_lake, iname_soil, iname_glce, iname_soil, or iname_aquifer'; return
        end select

        ! identify domain
        select case(ixDomainType)

          case(iname_veg)
            ! Note that scalarCanopyCm/iden_water is computed
            diffT = scalarCanopyTemp - Tfreeze
            if(diffT>=0._rkind)then
              scalarCanopyCm =  Cp_water * diffT
              ! derivatives
              dCm_dTkCanopy  = Cp_water
            else
              integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
              fLiq = fracliquid(scalarCanopyTemp,snowfrz_scale)
              scalarCanopyCm = Cp_water * integral + Cp_ice * (diffT - integral) 
              ! derivatives
              dfLiq_dT = dFracLiq_dTk(scalarCanopyTemp,snowfrz_scale)
              dCm_dTkCanopy = Cp_water * fLiq + Cp_ice * (1._rkind - fLiq)
            end if

          case(iname_snow, iname_lake, iname_glce)
            diffT = mLayerTemp(iLayer) - Tfreeze
            if(diffT>=0._rkind)then
              mLayerCm(iLayer) = (iden_water * Cp_water - iden_air * Cp_air) * diffT
              ! derivatives
              dCm_dTk(iLayer) = iden_water * Cp_water - iden_air * Cp_air
            else
              frz_scale_use = snowfrz_scale
              if(ixDomainType==iname_lake .or. ixDomainType==iname_glce) frz_scale_use = snowfrz_scale*icefrz_mult

              fLiq = fracliquid(mLayerTemp(iLayer),frz_scale_use,iLayer>nLayers-noThetaChange)
              if(iLayer>nLayers-noThetaChange) then
                mLayerCm(iLayer) = 0._rkind ! no change in total water content, so no change enthalpy with water
                dCm_dTk(iLayer) = 0._rkind
              else
                integral = (1._rkind/frz_scale_use) * atan(frz_scale_use * diffT)
                mLayerCm(iLayer) = (iden_water * Cp_ice - iden_air * Cp_air * iden_water/iden_ice) * ( diffT - integral ) &
                                       + (iden_water * Cp_water - iden_air * Cp_air) * integral
                ! derivatives
                dfLiq_dT = dFracLiq_dTk(mLayerTemp(iLayer),frz_scale_use,iLayer>nLayers-noThetaChange)
                dCm_dTk(iLayer) = (iden_water * Cp_ice - iden_air * Cp_air * iden_water/iden_ice) * ( 1._rkind -fLiq ) &
                                 + (iden_water * Cp_water - iden_air * Cp_air) * fLiq
              endif
            end if

          case(iname_soil)
            diffT = mLayerTemp(iLayer) - Tfreeze
            Tcrit = crit_soilT( mLayerMatricHead(ixControlIndex) )
            diff0 = Tcrit - Tfreeze
            if( mLayerTemp(iLayer)>=Tcrit)then
              mLayerCm(iLayer) = (-iden_air * Cp_air + iden_water * Cp_water) * diffT
              ! derivatives
              dCm_dTk(iLayer) = -iden_air * Cp_air + iden_water * Cp_water
              dCm_dPsi0(ixControlIndex) = 0._rkind
            else        
              mLayerCm(iLayer) = -iden_air * Cp_air * diffT + iden_ice * Cp_ice * (mLayerTemp(iLayer)-Tcrit) &
                                     + iden_water * Cp_water * diff0
              ! derivatives
              dTcrit_dPsi0 = merge(gravity*Tfreeze/LH_fus,0._rkind,mLayerMatricHead(ixControlIndex)<=0._rkind)
              dCm_dTk(iLayer) = -iden_air * Cp_air + iden_ice * Cp_ice
              dCm_dPsi0(ixControlIndex) = (-iden_ice * Cp_ice + iden_water * Cp_water) * dTcrit_dPsi0
            endif

        end select

      end if  ! if an energy layer
    end do  ! looping through state variables

  end associate

end subroutine heatAdvectWat


end module heat_Cp_Cm_module
