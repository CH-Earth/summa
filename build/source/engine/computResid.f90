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

module computResid_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)

implicit none
private
public::computResid
contains

 ! **********************************************************************************************************
 ! public subroutine computResid: compute the residual vector and the Jacobian matrix
 ! **********************************************************************************************************
 subroutine computResid(&
                        ! input: model control
                        dt,                      & ! intent(in):    length of the time step (seconds)
                        nSnow,                   & ! intent(in):    number of snow layers
                        nSoil,                   & ! intent(in):    number of soil layers
                        nLayers,                 & ! intent(in):    total number of layers
                        nState,                  & ! intent(in):    total number of state variables
                        canopyDepth,             & ! intent(in):    depth of the vegetation canopy (m)
                        computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                        ! input: flux vectors
                        sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                        fVec,                    & ! intent(in):    flux vector
                        ! input: state variables (already disaggregated into scalars and vectors)
                        scalarCanairTempTrial,   & ! intent(in):    trial value for the temperature of the canopy air space (K)
                        scalarCanopyTempTrial,   & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                        mLayerTempTrial,         & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                        mLayerVolFracWatTrial,   & ! intent(in):    trial vector of total volumetric total water content (-)
                        ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                        scalarCanopyLiqTrial,    & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                        scalarCanopyIceTrial,    & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                        mLayerVolFracLiqTrial,   & ! intent(in):    trial value for the volumetric liquid water content in each snow and soil layer (-)
                        mLayerVolFracIceTrial,   & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                        ! input: data structures
                        prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                        flux_data,               & ! intent(in):    model fluxes for a local HRU
                        indx_data,               & ! intent(in):    index data
                        ! output
                        rVec,                    & ! intent(out):   residual vector
                        err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE var_lookup,only:iLookPROG                     ! named variables for structure elements
 USE var_lookup,only:iLookDIAG                     ! named variables for structure elements
 USE var_lookup,only:iLookFLUX                     ! named variables for structure elements
 USE var_lookup,only:iLookINDEX                    ! named variables for structure elements
 implicit none
 ! input: model control
 real(dp),intent(in)             :: dt                        ! length of the time step (seconds)
 integer(i4b),intent(in)         :: nSnow                     ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                   ! total number of layers in the snow+soil domain
 integer(i4b),intent(in)         :: nState                    ! total number of state variables
 real(dp),intent(in)             :: canopyDepth               ! depth of the vegetation canopy (m)
 logical(lgt),intent(in)         :: computeVegFlux            ! flag to indicate if computing fluxes over vegetation
 ! input: flux vectors
 real(qp),intent(in)             :: sMul(:)   ! NOTE: qp      ! state vector multiplier (used in the residual calculations)
 real(dp),intent(in)             :: fVec(:)                   ! flux vector
 ! input: state variables (already disaggregated into scalars and vectors)
 real(dp),intent(in)             :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
 real(dp),intent(in)             :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
 real(dp),intent(in)             :: mLayerTempTrial(:)        ! trial value for temperature of each snow/soil layer (K)
 real(dp),intent(in)             :: mLayerVolFracWatTrial(:)  ! trial vector of total volumetric total water content (-)
 ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
 real(dp),intent(in)             :: scalarCanopyLiqTrial      ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)             :: scalarCanopyIceTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)             :: mLayerVolFracLiqTrial(:)  ! trial value for volumetric fraction of liquid water (-)
 real(dp),intent(in)             :: mLayerVolFracIceTrial(:)  ! trial value for volumetric fraction of ice (-)
 ! input: data structures
 type(var_dlength),intent(in)    :: prog_data                 ! prognostic variables for a local HRU
 type(var_dlength),intent(in)    :: diag_data                 ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)    :: flux_data                 ! model fluxes for a local HRU
 type(var_ilength),intent(in)    :: indx_data                 ! indices defining model states and layers
 ! output
 real(qp),intent(out)            :: rVec(:)   ! NOTE: qp      ! residual vector
 integer(i4b),intent(out)        :: err                       ! error code
 character(*),intent(out)        :: message                   ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp)                        :: scalarCanopyWat           ! total water stored on the vegetation canopy at the start of the time step (kg m-2)
 real(dp)                        :: scalarCanopyWatTrial      ! total water stored on the vegetation canopy (kg m-2)
 real(dp)                        :: vThetaInit(nSoil)         ! liquid equivalent of total water at the start of the step (-)
 real(dp)                        :: vThetaTrial(nSoil)        ! liquid equivalent of total water at the current iteration (-)
 real(dp)                        :: rAdd(nState)              ! additional (sink) terms on the RHS of the state equation
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! link to the necessary variables for the residual computations
 associate(&
  ! model state variables (vegetation canopy)
  scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the canopy air space (K)
  scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the vegetation canopy (K)
  scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,&  ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)
  scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,&  ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
  ! model state variables (snow and soil domains)
  mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,&  ! intent(inout): [dp(:)] temperature of each snow/soil layer (K)
  mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of liquid water (-)
  mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of ice (-)
  ! layer depth
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)
  ! model fluxes (sink terms in the soil domain)
  mLayerTranspire         => flux_data%var(iLookFLUX%mLayerTranspire)%dat           ,&  ! intent(out): [dp] transpiration loss from each soil layer (m s-1)
  mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,&  ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
  mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,&  ! intent(out): [dp(:)] change in storage associated with compression of the soil matrix (-)
  ! model indices
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,&  ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
  ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat             &  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
 ) ! association to necessary variables for the residual computations

 ! ---
 ! * compute sink terms...
 ! -----------------------

 ! intialize additional terms on the RHS as zero
 rAdd(:) = 0._dp

 ! compute energy associated with melt freeze for the vegetation canopy
 if(computeVegFlux)then
  rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*(scalarCanopyIceTrial - scalarCanopyIce)/canopyDepth   ! energy associated with melt/freeze (J m-3)
 endif

 ! compute energy associated with melt/freeze for snow
 if(nSnow>0)&
 rAdd(ixSnowOnlyNrg) = rAdd(ixSnowOnlyNrg) + LH_fus*iden_ice*(mLayerVolFracIceTrial(1:nSnow) - mLayerVolFracIce(1:nSnow))       ! energy associated with melt/freeze (J m-3)

 ! compute energy associated with melt/freeze for soil
 rAdd(ixSoilOnlyNrg) = rAdd(ixSoilOnlyNrg) + LH_fus*iden_water*(mLayerVolFracIceTrial(nSnow+1:nLayers) - mLayerVolFracIce(nSnow+1:nLayers))     ! energy associated with melt/freeze (J m-3)

 ! sink terms for water (-)
 ! NOTE: state variable is volumetric water content, so melt-freeze is not included
 ! NOTE: ground evaporation was already included in the flux at the upper boundary
 ! NOTE: rAdd(ixSnowOnlyWat)=0, and is defined in the initialization above
 rAdd(ixSoilOnlyHyd)    = rAdd(ixSoilOnlyHyd) + dt*(mLayerTranspire(1:nSoil) - mLayerBaseflow(1:nSoil) )/mLayerDepth(nSnow+1:nLayers) - mLayerCompress(1:nSoil)

 ! ---
 ! * compute the residual vector...
 ! --------------------------------

 ! vegetation state variables (if they exist)
 if(computeVegFlux)then

  ! get the total water in the vegetation canopy (kg m-2)
  scalarCanopyWat      = scalarCanopyIce      + scalarCanopyIce      ! NOTE: No need to compute this at every iteration
  scalarCanopyWatTrial = scalarCanopyIceTrial + scalarCanopyIceTrial

  ! compute the residual vector for the vegetation canopy
  ! NOTE: sMul(ixVegWat) = 1, but include as it converts all variables to quadruple precision
  ! --> energy balance
  rVec(ixCasNrg) = sMul(ixCasNrg)*scalarCanairTempTrial - ( (sMul(ixCasNrg)*scalarCanairTemp + fVec(ixCasNrg)*dt) + rAdd(ixCasNrg) )
  rVec(ixVegNrg) = sMul(ixVegNrg)*scalarCanopyTempTrial - ( (sMul(ixVegNrg)*scalarCanopyTemp + fVec(ixVegNrg)*dt) + rAdd(ixVegNrg) )
  ! --> mass balance
  rVec(ixVegWat) = sMul(ixVegWat)*scalarCanopyWatTrial  - ( (sMul(ixVegWat)*scalarCanopyWat  + fVec(ixVegWat)*dt) + rAdd(ixVegWat) )

 endif

 ! compute the residual vector for the snow and soil sub-domains for energy
 rVec(ixSnowSoilNrg) = sMul(ixSnowSoilNrg)*mLayerTempTrial(1:nLayers) - ( (sMul(ixSnowSoilNrg)*mLayerTemp(1:nLayers)  + fVec(ixSnowSoilNrg)*dt) + rAdd(ixSnowSoilNrg) )

 ! compute the residual vector for the **snow** sub-domain for liquid water
 if(nSnow>0)&
 rVec(ixSnowOnlyWat) = mLayerVolFracWatTrial(1:nSnow) - ( (mLayerVolFracWatInit(1:nSnow)  + fVec(ixSnowOnlyWat)*dt) + rAdd(ixSnowOnlyWat) )

 ! compute the residual vector for the **soil** sub-domain for liquid water
 ! NOTE: No need to compute vThetaInit at every iteration
 vThetaInit(1:nSoil)  = mLayerVolFracLiq(nSnow+1:nLayers)      + mLayerVolFracIce(nSnow+1:nLayers)      ! liquid equivalent of total water at the start of the step
 vThetaTrial(1:nSoil) = mLayerVolFracLiqLocal(nSnow+1:nLayers) + mLayerVolFracIceLocal(nSnow+1:nLayers) ! liquid equivalent of total water at the current iteration
 rVec(ixSoilOnlyHyd)  = vThetaTrial(1:nSoil) - ( (vThetaInit(1:nSoil) + fVec(ixSoilOnlyHyd)*dt) + rAdd(ixSoilOnlyHyd) )

 ! print result
 if(globalPrintFlag) write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec

 ! end association with the necessary variabiles for the residual calculations
 end associate

 end subroutine computResid

end module computResid_module
