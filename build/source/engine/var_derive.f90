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

module var_derive_module
USE nrtype
implicit none
private
public::calcHeight
public::rootDensty
public::satHydCond
public::fracFuture
public::v_shortcut
contains


 ! **********************************************************************************************************
 ! public subroutine calcHeight: compute snow height
 ! **********************************************************************************************************
 subroutine calcHeight(&
                       ! input/output: data structures
                       indx_data,   & ! intent(in): layer type
                       prog_data,   & ! intent(inout): model variables for a local HRU
                       ! output: error control
                       err,message)
 ! access named variables for snow and soil
 USE globalData,only:ix_soil,ix_snow            ! named variables for snow and soil
 ! access to the derived types to define the data structures
 USE data_types,only:var_ilength    ! x%var(:)%dat (i4b)
 USE data_types,only:var_dlength    ! x%var(:)%dat (dp)
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookPROG,iLookINDEX  ! named variables for structure elements
 implicit none
 ! ----------------------------------------------------------------------------------
 ! dummy variables
 ! input/output: data structures
 type(var_ilength),intent(in)       :: indx_data      ! type of model layer
 type(var_dlength),intent(inout)    :: prog_data      ! model variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)           :: err            ! error code
 character(*),intent(out)           :: message        ! error message
 ! local variables
 integer(i4b)                       :: iLayer         ! loop through layers
 integer(i4b)                       :: ixLower(1)     ! index of the lower bound
 ! ----------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='calcHeight/'
 ! ----------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! associate the model index structures
 nLayers        => indx_data%var(iLookINDEX%nLayers)%dat(1),  &   ! total number of layers
 layerType      => indx_data%var(iLookINDEX%layerType)%dat,   &   ! layer type (ix_soil or ix_snow)
 ! associate the values in the model variable structures
 mLayerDepth    => prog_data%var(iLookPROG%mLayerDepth)%dat,  &   ! depth of the layer (m)
 mLayerHeight   => prog_data%var(iLookPROG%mLayerHeight)%dat, &   ! height of the layer mid-point (m)
 iLayerHeight   => prog_data%var(iLookPROG%iLayerHeight)%dat  &   ! height of the layer interface (m)
 ) ! end associate
 ! ----------------------------------------------------------------------------------

 ! initialize layer height as the top of the snowpack -- positive downward
 ixLower=lbound(iLayerHeight); if(ixLower(1) > 0)then; err=20; message=trim(message)//'unexpected lower bound for iLayerHeight'; return; end if
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

 ! end association to variables in the data structure
 end associate

 end subroutine calcHeight


 ! **********************************************************************************************************
 ! public subroutine rootDensty: compute vertical distribution of root density
 ! **********************************************************************************************************
 subroutine rootDensty(parData,indx_data,prog_data,diag_data,err,message)
 ! model decision structures
 USE globalData,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! look-up values for the choice of the rooting profile
 USE mDecisions_module,only: &
 powerLaw,                   & ! simple power-law rooting profile
 doubleExp                     ! the double exponential function of Xeng et al. (JHM 2001)
 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only: &
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization
 ! named variables
 USE globalData,only:ix_soil,ix_snow                                  ! named variables defining snow and soil layers
 USE var_lookup,only:iLookPARAM,iLookINDEX,iLookPROG,iLookDIAG        ! named variables for structure elements
 ! data types
 USE data_types,only:var_dlength    ! x%var(:)%dat (dp)
 USE data_types,only:var_ilength    ! x%var(:)%dat (i4b)
 implicit none
 ! declare input variables
 real(dp),intent(in)             :: parData(:)      ! vector of model parameters
 type(var_ilength),intent(in)    :: indx_data       ! data structure of model indices for a local HRU
 type(var_dlength),intent(in)    :: prog_data       ! data structure of model prognostic (state) variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data       ! data structure of model diagnostic variables for a local HRU
 ! declare output variables
 integer(i4b),intent(out)        :: err             ! error code
 character(*),intent(out)        :: message         ! error message
 ! declare local variables
 integer(i4b)                    :: iLayer          ! loop through layers
 real(dp)                        :: fracRootLower   ! fraction of the rooting depth at the lower interface
 real(dp)                        :: fracRootUpper   ! fraction of the rooting depth at the upper interface

 ! initialize error control
 err=0; message='rootDensty/'

 ! ----------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! associate the model decisions
 ixRootProfile         =>model_decisions(iLookDECISIONS%rootProfil)%iDecision,  & ! choice of the rooting profile
 ixGroundwater         =>model_decisions(iLookDECISIONS%groundwatr)%iDecision,  & ! choice of groundwater parameterization
 ! associate the values in the model parameter structures
 rootScaleFactor1      =>parData(iLookPARAM%rootScaleFactor1),                  & ! 1st scaling factor (m-1)
 rootScaleFactor2      =>parData(iLookPARAM%rootScaleFactor2),                  & ! 2nd scaling factor (m-1)
 rootingDepth          =>parData(iLookPARAM%rootingDepth),                      & ! rooting depth (m)
 rootDistExp           =>parData(iLookPARAM%rootDistExp),                       & ! root distribution exponent (-)
 ! associate the model index structures
 nSnow                 =>indx_data%var(iLookINDEX%nSnow)%dat(1),                & ! number of snow layers
 nLayers               =>indx_data%var(iLookINDEX%nLayers)%dat(1),              & ! total number of layers
 layerType             =>indx_data%var(iLookINDEX%layerType)%dat,               & ! layer type (ix_soil or ix_snow)
 iLayerHeight          =>prog_data%var(iLookPROG%iLayerHeight)%dat,             & ! height of the layer interface (m)
 ! associate the values in the model variable structures
 scalarAquiferRootFrac =>diag_data%var(iLookDIAG%scalarAquiferRootFrac)%dat(1), & ! fraction of roots below the soil profile (in the aquifer)
 mLayerRootDensity     =>diag_data%var(iLookDIAG%mLayerRootDensity)%dat         & ! fraction of roots in each soil layer (-)
 ) ! end associate
 ! ----------------------------------------------------------------------------------

 ! compute the fraction of roots in each soil layer
 do iLayer=nSnow+1,nLayers

  ! different options for the rooting profile
  select case(ixRootProfile)

   ! ** option 1: simple power-law profile
   case(powerLaw)
    if(iLayerHeight(iLayer-1)<rootingDepth)then
     ! compute the fraction of the rooting depth at the lower and upper interfaces
     if(iLayer==nSnow+1)then  ! height=0; avoid precision issues
      fracRootLower = 0._dp
     else
      fracRootLower = iLayerHeight(iLayer-1)/rootingDepth
     end if
     fracRootUpper = iLayerHeight(iLayer)/rootingDepth
     if(fracRootUpper>1._dp) fracRootUpper=1._dp
     ! compute the root density
     mLayerRootDensity(iLayer-nSnow) = fracRootUpper**rootDistExp - fracRootLower**rootDistExp
   else
    mLayerRootDensity(iLayer-nSnow) = 0._dp
   end if

   ! ** option 2: double expoential profile of Zeng et al. (JHM 2001)
   case(doubleExp)
    ! compute the cumulative fraction of roots at the top and bottom of the layer
    fracRootLower = 1._dp - 0.5_dp*(exp(-iLayerHeight(iLayer-1)*rootScaleFactor1) + exp(-iLayerHeight(iLayer-1)*rootScaleFactor2) )
    fracRootUpper = 1._dp - 0.5_dp*(exp(-iLayerHeight(iLayer  )*rootScaleFactor1) + exp(-iLayerHeight(iLayer  )*rootScaleFactor2) )
    ! compute the root density
    mLayerRootDensity(iLayer-nSnow) = fracRootUpper - fracRootLower
    write(*,'(a,10(f11.5,1x))') 'mLayerRootDensity(iLayer-nSnow), fracRootUpper, fracRootLower = ', &
                                 mLayerRootDensity(iLayer-nSnow), fracRootUpper, fracRootLower

   ! ** check
   case default; err=20; message=trim(message)//'unable to identify option for rooting profile'; return

  end select

 end do  ! (looping thru layers)

 ! check that root density is less than one
 if(sum(mLayerRootDensity) > 1._dp + epsilon(rootingDepth))then
  message=trim(message)//'problem with the root density calaculation'
  err=20; return
 end if

 ! compute fraction of roots in the aquifer
 if(sum(mLayerRootDensity) < 1._dp)then
  scalarAquiferRootFrac = 1._dp - sum(mLayerRootDensity)
 else
  scalarAquiferRootFrac = 0._dp
 end if
 
 ! check that roots in the aquifer are appropriate
 if(scalarAquiferRootFrac > epsilon(rootingDepth))then
  if(ixGroundwater /= bigBucket)then
   select case(ixRootProfile)
    case(powerLaw);  message=trim(message)//'roots in the aquifer only allowed for the big bucket gw parameterization: check that rooting depth < soil depth'
    case(doubleExp); message=trim(message)//'roots in the aquifer only allowed for the big bucket gw parameterization: increase soil depth to alow for exponential roots'
   end select
   err=10; return
  end if  ! if not the big bucket
 end if  ! if roots in the aquifer

 end associate

 end subroutine rootDensty


 ! **********************************************************************************************************
 ! public subroutine satHydCond: compute vertical profile of saturated hydraulic conductivity
 ! **********************************************************************************************************
 subroutine satHydCond(parData,indx_data,prog_data,flux_data,err,message)
 ! model decision structures
 USE globalData,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only: &
  constant,                  & ! constant hydraulic conductivity with depth
  powerLaw_profile             ! power-law profile
 ! named variables
 USE globalData,only:ix_soil,ix_snow                                  ! named variables defining snow and soil layers
 USE var_lookup,only:iLookPARAM,iLookINDEX,iLookPROG,iLookFLUX        ! named variables for structure elements
 ! data types
 USE data_types,only:var_dlength    ! x%var(:)%dat (dp)
 USE data_types,only:var_ilength    ! x%var(:)%dat (i4b)
 implicit none
 ! declare input variables
 real(dp),intent(in)             :: parData(:)      ! vector of model parameters
 type(var_ilength),intent(in)    :: indx_data       ! data structure of model indices for a local HRU
 type(var_dlength),intent(in)    :: prog_data       ! data structure of model prognostic (state) variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data       ! data structure of model fluxes for a local HRU
 ! declare output variables
 integer(i4b),intent(out)        :: err             ! error code
 character(*),intent(out)        :: message         ! error message
 ! declare local variables
 integer(i4b)                    :: iLayer          ! loop through layers
 ! initialize error control
 err=0; message='satHydCond/'
 ! ----------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! associate the values in the parameter structures
 k_soil             => parData(iLookPARAM%k_soil),                      & ! saturated hydraulic conductivity at the compacted depth (m s-1)
 k_macropore        => parData(iLookPARAM%k_macropore),                 & ! saturated hydraulic conductivity at the compacted depth for macropores (m s-1)
 compactedDepth     => parData(iLookPARAM%compactedDepth),              & ! the depth at which k_soil reaches the compacted value given by CH78 (m)
 zScale_TOPMODEL    => parData(iLookPARAM%zScale_TOPMODEL),             & ! exponent for the TOPMODEL-ish baseflow parameterization (-)
 ! associate the model index structures
 nSnow              => indx_data%var(iLookINDEX%nSnow)%dat(1),          & ! number of snow layers
 nLayers            => indx_data%var(iLookINDEX%nLayers)%dat(1),        & ! total number of layers
 layerType          => indx_data%var(iLookINDEX%layerType)%dat,         & ! layer type (ix_soil or ix_snow)
 ! associate the coordinate variables
 mLayerHeight       => prog_data%var(iLookPROG%mLayerHeight)%dat,       & ! height at the mid-point of each layer (m)
 iLayerHeight       => prog_data%var(iLookPROG%iLayerHeight)%dat,       & ! height at the interface of each layer (m)
 ! associate the values in the model variable structures
 mLayerSatHydCondMP => flux_data%var(ilookFLUX%mlayersathydcondmp)%dat, & ! saturated hydraulic conductivity for macropores at the mid-point of each layer (m s-1)
 mLayerSatHydCond   => flux_data%var(ilookFLUX%mlayersathydcond)%dat,   & ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 iLayerSatHydCond   => flux_data%var(ilookFLUX%ilayersathydcond)%dat    & ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 ) ! end associate
 ! ----------------------------------------------------------------------------------

 ! loop through soil layers
 ! NOTE: could do constant profile with the power-law profile with exponent=1, but keep constant profile decision for clarity
 do iLayer=nSnow,nLayers
  select case(model_decisions(iLookDECISIONS%hc_profile)%iDecision)
   ! constant hydraulic conductivity with depth
   case(constant)
    iLayerSatHydCond(iLayer-nSnow) = k_soil
    if(iLayer > nSnow)then ! avoid layer 0
     mLayerSatHydCond(iLayer-nSnow)   = k_soil
     mLayerSatHydCondMP(iLayer-nSnow) = k_macropore
    end if  ! if the mid-point of a layer
   ! power-law profile
   case(powerLaw_profile)
    ! (saturated hydraulic conductivity at layer interfaces)
    iLayerSatHydCond(iLayer-nSnow) = k_soil * ( (1._dp - iLayerHeight(iLayer)/iLayerHeight(nLayers))**(zScale_TOPMODEL - 1._dp) ) &
                                            / ( (1._dp -       compactedDepth/iLayerHeight(nLayers))**(zScale_TOPMODEL - 1._dp) )
    ! (saturated hydraulic conductivity at layer mid-points)
    if(iLayer > nSnow)then ! avoid layer 0
     ! (--> micropores)
     mLayerSatHydCond(iLayer-nSnow) = k_soil * ( (1._dp - mLayerHeight(iLayer)/iLayerHeight(nLayers))**(zScale_TOPMODEL - 1._dp) ) &
                                             / ( (1._dp -       compactedDepth/iLayerHeight(nLayers))**(zScale_TOPMODEL - 1._dp) )
     ! (--> macropores)
     mLayerSatHydCondMP(iLayer-nSnow) = k_macropore * ( (1._dp - mLayerHeight(iLayer)/iLayerHeight(nLayers))**(zScale_TOPMODEL - 1._dp) ) &
                                                    / ( (1._dp -       compactedDepth/iLayerHeight(nLayers))**(zScale_TOPMODEL - 1._dp) )
     !print*, 'compactedDepth = ', compactedDepth
     !print*, 'k_macropore    = ', k_macropore
     !print*, 'mLayerHeight(iLayer) = ', mLayerHeight(iLayer)
     !print*, 'iLayerHeight(nLayers) = ', iLayerHeight(nLayers)
     !print*, 'iLayer, mLayerSatHydCondMP(iLayer-nSnow) = ', mLayerSatHydCondMP(iLayer-nSnow)
    end if  ! if the mid-point of a layer
   ! error check (errors checked earlier also, so should not get here)
   case default
    message=trim(message)//"unknown hydraulic conductivity profile [option="//trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)//"]"
    err=10; return
  end select
  !if(iLayer > nSnow)& ! avoid layer 0
  ! write(*,'(i4,1x,2(f11.5,1x,e20.10,1x))') iLayer, mLayerHeight(iLayer), mLayerSatHydCond(iLayer-nSnow), iLayerHeight(iLayer), iLayerSatHydCond(iLayer-nSnow)
 end do  ! looping through soil layers
 !print*, trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)
 !print*, 'k_soil, k_macropore, zScale_TOPMODEL = ', k_soil, k_macropore, zScale_TOPMODEL
 !pause ' in satHydCond'
 end associate

 end subroutine satHydCond


 ! **********************************************************************************************************
 ! public subroutine fracFuture: compute the fraction of runoff in future time steps
 ! **********************************************************************************************************
 subroutine fracFuture(bpar_data,bvar_data,err,message)
 ! external functions
 USE soil_utils_module,only:gammp                     ! compute the cumulative probabilty based on the Gamma distribution
 ! model decision structures
 USE globalData,only:model_decisions                  ! model decision structure
 USE var_lookup,only:iLookDECISIONS                   ! named variables for elements of the decision structure
 ! look-up values for the sub-grid routing method
 USE mDecisions_module,only:      &
  timeDelay,&  ! time-delay histogram
  qInstant     ! instantaneous routing
 ! named variables 
 USE globalData,only:data_step                        ! time step of forcing data
 USE var_lookup,only:iLookBVAR,iLookBPAR              ! named variables for structure elements
 ! data types
 USE data_types,only:var_dlength    ! x%var(:)%dat (dp)
 implicit none
 ! input variables
 real(dp),intent(in)             :: bpar_data(:)           ! vector of basin-average model parameters
 ! output variables
 type(var_dlength),intent(inout) :: bvar_data              ! data structure of basin-average model variables
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! internal
 real(dp)                        :: dt                     ! data time step (s)
 integer(i4b)                    :: nTDH                   ! number of points in the time-delay histogram
 integer(i4b)                    :: iFuture                ! index in time delay histogram
 real(dp)                        :: aLambda                ! scale parameter in the Gamma distribution
 real(dp)                        :: tFuture                ! future time (end of step)
 real(dp)                        :: pSave                  ! cumulative probability at the start of the step
 real(dp)                        :: cumProb                ! cumulative probability at the end of the step
 real(dp)                        :: sumFrac                ! sum of runoff fractions in all steps
 real(dp),parameter              :: tolerFrac=0.01_dp      ! tolerance for fractional runoff
 ! initialize error control
 err=0; message='fracFuture/'
 ! ----------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ixRouting         => model_decisions(iLookDECISIONS%subRouting)%iDecision, & ! index for routing method
 routingGammaShape => bpar_data(iLookBPAR%routingGammaShape),               & ! shape parameter in Gamma distribution used for sub-grid routing (-)
 routingGammaScale => bpar_data(iLookBPAR%routingGammaScale),               & ! scale parameter in Gamma distribution used for sub-grid routing (s)
 runoffFuture      => bvar_data%var(iLookBVAR%routingRunoffFuture)%dat,     & ! runoff in future time steps (m s-1)
 fractionFuture    => bvar_data%var(iLookBVAR%routingFractionFuture)%dat    & ! fraction of runoff in future time steps (-)
 ) ! end associate
 ! ----------------------------------------------------------------------------------

 ! define time step
 dt =  data_step ! length of the data step (s)

 ! identify number of points in the time-delay histogram
 nTDH = size(runoffFuture)

 ! initialize runoffFuture
 runoffFuture(1:nTDH) = 0._dp

 !print*, 'nTDH = ', nTDH

 ! select option for sub-grid routing
 select case(ixRouting)

  ! ** instantaneous routing
  case(qInstant)
   fractionFuture(1)      = 1._dp
   fractionFuture(2:nTDH) = 0._dp

  ! ** time delay histogram
  case(timeDelay)
   ! initialize
   pSave   = 0._dp ! cumulative probability at the start of the step
   aLambda = routingGammaShape / routingGammaScale
   if(routingGammaShape <= 0._dp .or. aLambda < 0._dp)then
    message=trim(message)//'bad arguments for the Gamma distribution'
    err=20; return
   end if
   ! loop through time steps and compute fraction of runoff in future steps
   do iFuture = 1,nTDH
    tFuture = real(iFuture, kind(dt))*dt                  ! future time (end of step)
    cumProb = gammp(routingGammaShape,aLambda*tFuture)    ! cumulative probability at the end of the step
    fractionFuture(iFuture) = max(0._dp, cumProb - pSave) ! fraction of runoff in the current step
    pSave   = cumProb                                     ! save the cumulative probability for use in the next step
    if(fractionFuture(iFuture) < tiny(dt))then
     fractionFuture(iFuture:nTDH) = 0._dp
     exit
    end if
    !write(*,'(a,1x,i4,1x,3(f20.10,1x))') trim(message), iFuture, tFuture, cumProb, fractionFuture(iFuture)
   end do ! (looping through future time steps)
   ! check that we have enough bins
   sumFrac  = sum(fractionFuture)
   if(abs(1._dp - sumFrac) > tolerFrac)then
    message=trim(message)//'not enough bins for the time delay histogram -- fix hard-coded parameter in alloc_bvar'
    err=20; return
   end if
   ! ensure the fraction sums to one
   fractionFuture = fractionFuture/sumFrac

  ! ** error checking
  case default; err=20; message=trim(message)//'cannot find option for sub-grid routing'; return

 end select ! (select option for sub-grid routing)

 end associate

 end subroutine fracFuture


 ! **********************************************************************************************************
 ! public subroutine v_shortcut: compute "short-cut" variables
 ! **********************************************************************************************************
 subroutine v_shortcut(parData,diag_data,err,message)
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
 ! named variables
 USE globalData,only:ix_soil,ix_snow           ! named variables defining snow and soil layers
 USE var_lookup,only:iLookPARAM,iLookDIAG      ! named variables for structure elements
 ! data types
 USE data_types,only:var_dlength    ! x%var(:)%dat (dp)
 implicit none
 ! declare input variables
 real(dp),intent(in)             :: parData(:)      ! vector of model parameters
 type(var_dlength),intent(inout) :: diag_data         ! data structure of model variables for a local HRU
 ! declare output variables
 integer(i4b),intent(out)        :: err             ! error code
 character(*),intent(out)        :: message         ! error message
 ! declare local variables
 real(dp)                        :: bulkden_soil    ! bulk density of soil (kg m-3)
 ! initialize error control
 err=0; message='v_shortcut/'
 ! ----------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! associate values in the parameter structures
 iden_soil      =>parData(iLookPARAM%soil_dens_intr),                   & ! intrinsic soil density (kg m-3)
 frac_sand      =>parData(iLookPARAM%frac_sand),                        & ! fraction of sand (-)
 frac_silt      =>parData(iLookPARAM%frac_silt),                        & ! fraction of silt (-)
 frac_clay      =>parData(iLookPARAM%frac_clay),                        & ! fraction of clay (-)
 theta_sat      =>parData(iLookPARAM%theta_sat),                        & ! soil porosity (-)
 vGn_n          =>parData(iLookPARAM%vGn_n),                            & ! van Genutchen "n" parameter (-)
 ! associate values in the model variable structures
 vGn_m          =>diag_data%var(iLookDIAG%scalarVGn_m)%dat(1),          & ! van Genutchen "m" parameter (-)
 kappa          =>diag_data%var(iLookDIAG%scalarKappa)%dat(1),          & ! constant in the freezing curve function (m K-1)
 volHtCap_air   =>diag_data%var(iLookDIAG%scalarVolHtCap_air)%dat(1),   & ! volumetric heat capacity of air (J m-3 K-1)
 volHtCap_ice   =>diag_data%var(iLookDIAG%scalarVolHtCap_ice)%dat(1),   & ! volumetric heat capacity of ice (J m-3 K-1)
 volHtCap_soil  =>diag_data%var(iLookDIAG%scalarVolHtCap_soil)%dat(1),  & ! volumetric heat capacity of soil (J m-3 K-1)
 volHtCap_water =>diag_data%var(iLookDIAG%scalarVolHtCap_water)%dat(1), & ! volumetric heat capacity of water (J m-3 K-1)
 volLatHt_fus   =>diag_data%var(iLookDIAG%scalarvolLatHt_fus)%dat(1),   & ! volumetric latent heat of fusion (J m-3)
 lambda_drysoil =>diag_data%var(iLookDIAG%scalarLambda_drysoil)%dat(1), & ! thermal conductivity of dry soil (W m-1)
 lambda_wetsoil =>diag_data%var(iLookDIAG%scalarLambda_wetsoil)%dat(1)  & ! thermal conductivity of wet soil (W m-1)
 ) ! end associate
 ! ----------------------------------------------------------------------------------

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
 end associate

 end subroutine v_shortcut


end module var_derive_module
