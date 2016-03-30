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

module eval8summa_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! access the global print flag
USE globalData,only:globalPrintFlag

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! layer types
USE globalData,only:ix_soil,ix_snow ! named variables for snow and soil

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin                   ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation

implicit none
private
public::eval8summa

! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=tiny(1.0_dp)     ! a very small number
real(dp),parameter  :: veryBig=1.e+20_dp          ! a very big number
real(dp),parameter  :: dx = 1.e-8_dp              ! finite difference increment

contains

 ! **********************************************************************************************************
 ! public subroutine eval8summa: compute the residual vector and the Jacobian matrix
 ! **********************************************************************************************************
 subroutine eval8summa(&
                       ! input: model control
                       dt,                      & ! intent(in):    length of the time step (seconds)
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                       firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                       computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                       ! input: state vectors
                       stateVecTrial,           & ! intent(in):    model state vector
                       sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                       ! input: data structures
                       model_decisions,         & ! intent(in):    model decisions
                       type_data,               & ! intent(in):    type of vegetation and soil
                       attr_data,               & ! intent(in):    spatial attributes
                       mpar_data,               & ! intent(in):    model parameters
                       forc_data,               & ! intent(in):    model forcing data
                       bvar_data,               & ! intent(in):    average model variables for the entire basin
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       indx_data,               & ! intent(in):    index data
                       ! input-output: data structures
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,               & ! intent(inout): model fluxes for a local HRU
                       deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ! output
                       fluxVec,                 & ! intent(out):   flux vector
                       resVec,                  & ! intent(out):   residual vector
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! provide access to subroutines
 USE getVectorz_module,only:varExtract            ! extract variables from the state vector
 USE computFlux_module,only:computFlux            ! compute fluxes given a state vector
 ! provide access to indices that define elements of the data structures
 USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
 USE var_lookup,only:iLookTYPE                    ! named variables for structure elements
 USE var_lookup,only:iLookATTR                    ! named variables for structure elements
 USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
 USE var_lookup,only:iLookFORCE                   ! named variables for structure elements
 USE var_lookup,only:iLookBVAR                    ! named variables for structure elements
 USE var_lookup,only:iLookPROG                    ! named variables for structure elements
 USE var_lookup,only:iLookINDEX                   ! named variables for structure elements
 USE var_lookup,only:iLookDIAG                    ! named variables for structure elements
 USE var_lookup,only:iLookFLUX                    ! named variables for structure elements
 USE var_lookup,only:iLookDERIV                   ! named variables for structure elements
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                    ! length of the time step (seconds)
 integer(i4b),intent(in)         :: nSnow                 ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                 ! number of soil layers
 integer(i4b),intent(in)         :: nLayers               ! total number of layers
 logical(lgt),intent(in)         :: firstSubStep          ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall         ! flag to indicate if we are processing the first flux call
 logical(lgt),intent(in)         :: computeVegFlux        ! flag to indicate if computing fluxes over vegetation
 ! input: state vectors
 real(dp),intent(in)             :: stateVecTrial(:)      ! model state vector 
 real(qp),intent(in)             :: sMul(:)   ! NOTE: qp  ! state vector multiplier (used in the residual calculations)
 ! input: data structures
 type(model_options),intent(in)  :: model_decisions(:)    ! model decisions
 type(var_i),        intent(in)  :: type_data             ! type of vegetation and soil
 type(var_d),        intent(in)  :: attr_data             ! spatial attributes
 type(var_d),        intent(in)  :: mpar_data             ! model parameters
 type(var_d),        intent(in)  :: forc_data             ! model forcing data
 type(var_dlength),  intent(in)  :: bvar_data             ! model variables for the local basin
 type(var_dlength),  intent(in)  :: prog_data             ! prognostic variables for a local HRU
 type(var_ilength),  intent(in)  :: indx_data             ! indices defining model states and layers
 ! output: data structures
 type(var_dlength),intent(inout) :: diag_data             ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data             ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: deriv_data            ! derivatives in model fluxes w.r.t. relevant state variables
 ! output: flux and residual vectors
 real(dp),intent(out)            :: fluxVec(:)            ! flux vector
 real(qp),intent(out)            :: resVec(:) ! NOTE: qp  ! residual vector
 ! output: error control
 integer(i4b),intent(out)        :: err                   ! error code
 character(*),intent(out)        :: message               ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! state variables
 real(dp)                        :: scalarCanairTempTrial ! trial value for temperature of the canopy air space (K)
 real(dp)                        :: scalarCanopyTempTrial ! trial value for temperature of the vegetation canopy (K)
 real(dp)                        :: scalarCanopyWatTrial  ! trial value for liquid water storage in the canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerTempTrial       ! trial value for temperature of layers in the snow and soil domains (K)
 real(dp),dimension(nLayers)     :: mLayerVolFracWatTrial ! trial value for volumetric fraction of total water (-)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadTrial ! trial value for matric head (m)
 ! diagnostic variables
 real(dp)                        :: fracLiqVeg            ! fraction of liquid water on vegetation (-)
 real(dp),dimension(nSnow)       :: fracLiqSnow           ! fraction of liquid water in each snow layer (-)
 real(dp)                        :: scalarCanopyLiqTrial  ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyIceTrial  ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerVolFracLiqTrial ! trial value for volumetric fraction of liquid water (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracIceTrial ! trial value for volumetric fraction of ice (-)
 ! other local variables
 real(dp),dimension(nLayers)     :: mLayerVolFracWatInit  ! initial value for volumetric fraction of total water (-)
 real(dp)                        :: canopyDepth           ! depth of the vegetationb canopy 
 character(LEN=256)              :: cmessage              ! error message of downwind routine




 ! --------------------------------------------------------------------------------------------------------------------------------
 ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 associate(&

 ! vegetation parameters
 nVegState               => indx_data%var(iLookINDEX%nVegState)%dat(1)             ,&  ! intent(in): [i4b]    number of vegetation state variables
 heightCanopyTop         => mpar_data%var(iLookPARAM%heightCanopyTop)              ,&  ! intent(in): [dp] height of the top of the vegetation canopy (m)
 heightCanopyBottom      => mpar_data%var(iLookPARAM%heightCanopyBottom)           ,&  ! intent(in): [dp] height of the bottom of the vegetation canopy (m)

 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)                ,&  ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)

 ! soil parameters
 vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,&  ! intent(in): [dp] soil residual volumetric water content (-)
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)              ,&  ! intent(in): [dp] specific storage coefficient (m-1)
 fImpede                 => mpar_data%var(iLookPARAM%f_impede)                     ,&  ! intent(in): [dp] ice impedance parameter (-)

 ! diagnostic variables
 scalarBulkVolHeatCapVeg => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),&  ! intent(in): [dp   ] bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,&  ! intent(in): [dp(:)] bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)

 ! model state variables (vegetation canopy)
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the vegetation canopy (K)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,&  ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,&  ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)

 ! model state variables (ponded water)
 scalarSfcMeltPond       => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in): [dp] ponded water caused by melt of the "snow without a layer" (kg m-2)

 ! model state variables (snow and soil domains)
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,&  ! intent(inout): [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of ice (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of liquid water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,&  ! intent(inout): [dp(:)] matric head (m)
 scalarAquiferStorage    => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)   ,&  ! intent(inout): [dp   ] aquifer storage (m)

 ! indices for specific state variables
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
 ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat             &  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain

 ) ! association to variables in the data structures

 ! define canopy depth
 if(computeVegFlux)then
  canopyDepth = heightCanopyTop - heightCanopyBottom
 else
  canopyDepth = realMissing
 endif


 ! extract variables from the model state vector
 call varExtract(&
                 ! input
                 stateVecTrial,                             & ! intent(in):    model state vector (mixed units)
                 indx_data,                                 & ! intent(in):    indices defining model states and layers
                 snowfrz_scale,                             & ! intent(in):    scaling parameter for the snow freezing curve (K-1)
                 vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m, & ! intent(in):    van Genutchen soil parameters
                 ! output: variables for the vegetation canopy
                 fracLiqVeg,                                & ! intent(out):   fraction of liquid water on the vegetation canopy (-)
                 scalarCanairTempTrial,                     & ! intent(out):   trial value of canopy air temperature (K)
                 scalarCanopyTempTrial,                     & ! intent(out):   trial value of canopy temperature (K)
                 scalarCanopyWatTrial,                      & ! intent(out):   trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,                      & ! intent(out):   trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,                      & ! intent(out):   trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 fracLiqSnow,                               & ! intent(out):   volumetric fraction of water in each snow layer (-)
                 mLayerTempTrial,                           & ! intent(out):   trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,                     & ! intent(out):   trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,                     & ! intent(out):   trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,                     & ! intent(out):   trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,                     & ! intent(out):   trial vector of matric head (m)
                 ! output: error control
                 err,cmessage)                                ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! compute the fluxes for a given state vector
 call computFlux(&
                 ! input-output: model control
                 nSnow,                   & ! intent(in):    number of snow layers
                 nSoil,                   & ! intent(in):    number of soil layers
                 nLayers,                 & ! intent(in):    total number of layers
                 firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                 firstFluxCall,           & ! intent(inout): flag to denote the first flux call
                 computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                 canopyDepth,             & ! intent(in):    depth of the vegetation canopy (m)
                 scalarSfcMeltPond/dt,    & ! intent(in):    drainage from the surface melt pond (kg m-2 s-1)
                 ! input: state variables
                 scalarCanairTempTrial,   & ! intent(in):    trial value for the temperature of the canopy air space (K)
                 scalarCanopyTempTrial,   & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                 mLayerTempTrial,         & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                 mLayerMatricHeadTrial,   & ! intent(in):    trial value for the matric head in each soil layer (m)
                 ! input: diagnostic variables defining the liquid water and ice content
                 scalarCanopyLiqTrial,    & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                 scalarCanopyIceTrial,    & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                 mLayerVolFracLiqTrial,   & ! intent(in):    trial value for the volumetric liquid water content in each snow and soil layer (-)
                 mLayerVolFracIceTrial,   & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                 ! input: data structures
                 model_decisions,         & ! intent(in):    model decisions
                 type_data,               & ! intent(in):    type of vegetation and soil
                 attr_data,               & ! intent(in):    spatial attributes
                 mpar_data,               & ! intent(in):    model parameters
                 forc_data,               & ! intent(in):    model forcing data
                 bvar_data,               & ! intent(in):    average model variables for the entire basin
                 prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                 indx_data,               & ! intent(in):    index data
                 ! input-output: data structures
                 diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                 flux_data,               & ! intent(inout): model fluxes for a local HRU
                 deriv_data,              & ! intent(out):   derivatives in model fluxes w.r.t. relevant state variables
                 ! output: flux vector
                 fluxVec,                 & ! intent(out):   flux vector (mixed units)
                 ! output: error control
                 err,cmessage)              ! intent(out):   error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! end association with the information in the data structures
 end associate

 ! dummy specification of residuals
 resVec(:) = 0._dp


 end subroutine eval8summa
end module eval8summa_module
