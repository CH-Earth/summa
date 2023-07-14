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

module computJacob_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables to describe the state variable type
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ku             ! number of super-diagonal bands, assume ku>=3
USE globalData,only: kl             ! number of sub-diagonal bands, assume kl>=4
USE globalData,only: ixDiag         ! index for the diagonal band
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
 qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                       & ! a big bucket (lumped aquifer model)
 noExplicit                         ! no explicit groundwater parameterization

implicit none
! define constants
real(rkind),parameter     :: verySmall=tiny(1.0_rkind)     ! a very small number

private
public::computJacob
#ifdef SUNDIALS_ACTIVE
public::computJacob4kinsol
#endif

contains


! **********************************************************************************************************
! public subroutine computJacob: compute the Jacobian matrix
! **********************************************************************************************************
subroutine computJacob(&
             ! input: model control
					   dt,                         & ! intent(in):    length of the time step (seconds)
					   nSnow,                      & ! intent(in):    number of snow layers
					   nSoil,                      & ! intent(in):    number of soil layers
					   nLayers,                    & ! intent(in):    total number of layers
					   computeVegFlux,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
					   computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
					   ixMatrix,                   & ! intent(in):    form of the Jacobian matrix
					   ! input: data structures
					   indx_data,                  & ! intent(in):    index data
					   prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
					   diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
					   deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
					   dBaseflow_dMatric,          & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
					   ! input-output: Jacobian and its diagonal
					   dMat,                       & ! intent(inout): diagonal of the Jacobian matrix
					   aJac,                       & ! intent(out):   Jacobian matrix
					   ! output: error control
					   err,message)                  ! intent(out):   error code and error message
  ! -----------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  real(rkind),intent(in)            :: dt              ! length of the time step (seconds)
  integer(i4b),intent(in)           :: nSnow           ! number of snow layers
  integer(i4b),intent(in)           :: nSoil           ! number of soil layers
  integer(i4b),intent(in)           :: nLayers         ! total number of layers in the snow+soil domain
  logical(lgt),intent(in)           :: computeVegFlux  ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)           :: computeBaseflow ! flag to indicate if computing baseflow
  integer(i4b),intent(in)           :: ixMatrix        ! form of the Jacobian matrix
  ! input: data structures
  type(var_ilength),intent(in)      :: indx_data       ! indices defining model states and layers
  type(var_dlength),intent(in)      :: prog_data       ! prognostic variables for a local HRU
  type(var_dlength),intent(in)      :: diag_data       ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)      :: deriv_data      ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(in)            :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
  ! input-output: Jacobian and its diagonal
  real(rkind),intent(inout)         :: dMat(:)         ! diagonal of the Jacobian matrix
  real(rkind),intent(out)           :: aJac(:,:)       ! Jacobian matrix
  ! output variables
  integer(i4b),intent(out)          :: err             ! error code
  character(*),intent(out)          :: message         ! error message
  ! --------------------------------------------------------------
  ! * local variables
  ! --------------------------------------------------------------
  ! indices of model state variables
  integer(i4b)                      :: jState          ! index of state within the state subset
  integer(i4b)                      :: qState          ! index of cross-derivative state variable for baseflow
  integer(i4b)                      :: nrgState        ! energy state variable
  integer(i4b)                      :: watState        ! hydrology state variable
  integer(i4b)                      :: nState          ! number of state variables
  ! indices of model layers
  integer(i4b)                      :: iLayer          ! index of model layer
  integer(i4b)                      :: jLayer          ! index of model layer within the full state vector (hydrology)
  integer(i4b)                      :: pLayer          ! indices of soil layers (used for the baseflow derivatives)
  ! conversion factors
  real(rkind)                       :: convLiq2tot     ! factor to convert liquid water derivative to total water derivative
  ! --------------------------------------------------------------
  ! associate variables from data structures
  associate(&
    ! indices of model state variables
    ixCasNrg                     => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)                       ,& ! intent(in): [i4b] index of canopy air space energy state variable
    ixVegNrg                     => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)                       ,& ! intent(in): [i4b] index of canopy energy state variable
    ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                       ,& ! intent(in): [i4b] index of canopy hydrology state variable (mass)
    ixTopNrg                     => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)                       ,& ! intent(in): [i4b] index of upper-most energy state in the snow+soil subdomain
    ixTopHyd                     => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)                       ,& ! intent(in): [i4b] index of upper-most hydrology state in the snow+soil subdomain
    ixAqWat                      => indx_data%var(iLookINDEX%ixAqWat)%dat(1)                        ,& ! intent(in): [i4b] index of water storage in the aquifer
    ! vectors of indices for specfic state types within specific sub-domains IN THE FULL STATE VECTOR
    ixNrgLayer                   => indx_data%var(iLookINDEX%ixNrgLayer)%dat                        ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
    ! vector of energy indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilNrg                => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowOnlyNrg                => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow domain
    ixSoilOnlyNrg                => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
    ! vector of hydrology indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilHyd                => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    ixSnowOnlyHyd                => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow domain
    ixSoilOnlyHyd                => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! number of state variables of a specific type
    nSnowSoilNrg                 => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
    nSnowOnlyNrg                 => indx_data%var(iLookINDEX%nSnowOnlyNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
    nSoilOnlyNrg                 => indx_data%var(iLookINDEX%nSoilOnlyNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
    nSnowSoilHyd                 => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
    nSnowOnlyHyd                 => indx_data%var(iLookINDEX%nSnowOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    ! type and index of model control volume
    ixHydType                    => indx_data%var(iLookINDEX%ixHydType)%dat                         ,& ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
    ! mapping between states and model layers
    ixMapSubset2Full             => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat                  ,& ! intent(in): [i4b(:)] list of indices in the full state vector that are in the state subset
    ixMapFull2Subset             => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat                  ,& ! intent(in): [i4b(:)] list of indices in the state subset in each element of the full state vector
    ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
    dCanairNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanairTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. canopy air temperature
    dCanairNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanopyTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. canopy temperature
    dCanairNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dGroundTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. ground temperature
    dCanopyNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanairTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. canopy air temperature
    dCanopyNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanopyTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. canopy temperature
    dCanopyNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dGroundTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. ground temperature
    dCanopyNetFlux_dCanWat       => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanWat      )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy fluxes w.r.t. canopy total water content
    dGroundNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanairTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net ground flux w.r.t. canopy air temperature
    dGroundNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanopyTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net ground flux w.r.t. canopy temperature
    dGroundNetFlux_dCanWat       => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanWat      )%dat(1)  ,& ! intent(in): [dp]     derivative in net ground fluxes w.r.t. canopy total water content
    ! derivatives in evaporative fluxes w.r.t. relevant state variables
    dCanopyEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanair )%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy air temperature
    dCanopyEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanopy )%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy temperature
    dCanopyEvaporation_dTGround  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTGround )%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. ground temperature
    dCanopyEvaporation_dCanWat   => deriv_data%var(iLookDERIV%dCanopyEvaporation_dCanWat  )%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy total water content
    dGroundEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanair )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy air temperature
    dGroundEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanopy )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy temperature
    dGroundEvaporation_dTGround  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. ground temperature
    dGroundEvaporation_dCanWat   => deriv_data%var(iLookDERIV%dGroundEvaporation_dCanWat  )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy total water content
    ! derivatives in canopy water w.r.t canopy temperature
    dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy            )%dat(1)  ,& ! intent(in): [dp]     derivative of canopy liquid storage w.r.t. temperature
    dTheta_dTkCanopy             => deriv_data%var(iLookDERIV%dTheta_dTkCanopy            )%dat(1)  ,& ! intent(in): [dp]     derivative of volumetric liquid water content w.r.t. temperature
    ! derivatives in canopy liquid fluxes w.r.t. canopy water
    scalarCanopyLiqDeriv         => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv        )%dat(1)  ,& ! intent(in): [dp]     derivative in (throughfall + drainage) w.r.t. canopy liquid water
    ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
    dNrgFlux_dTempAbove          => deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove         )%dat     ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer above
    dNrgFlux_dTempBelow          => deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow         )%dat     ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer below
    ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. water state in layers above and below
    dNrgFlux_dWatAbove           => deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove          )%dat     ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer above
    dNrgFlux_dWatBelow           => deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow          )%dat     ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer below
    ! derivatives in soil transpiration w.r.t. canopy state variables
    mLayerdTrans_dTCanair        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanair       )%dat     ,& ! intent(in): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
    mLayerdTrans_dTCanopy        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanopy       )%dat     ,& ! intent(in): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
    mLayerdTrans_dTGround        => deriv_data%var(iLookDERIV%mLayerdTrans_dTGround       )%dat     ,& ! intent(in): derivatives in the soil layer transpiration flux w.r.t. ground temperature
    mLayerdTrans_dCanWat         => deriv_data%var(iLookDERIV%mLayerdTrans_dCanWat        )%dat     ,& ! intent(in): derivatives in the soil layer transpiration flux w.r.t. canopy total water
    ! derivatives in aquifer transpiration w.r.t. canopy state variables
    dAquiferTrans_dTCanair       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanair      )%dat(1)  ,&  ! intent(in): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
    dAquiferTrans_dTCanopy       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanopy      )%dat(1)  ,&  ! intent(in): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
    dAquiferTrans_dTGround       => deriv_data%var(iLookDERIV%dAquiferTrans_dTGround      )%dat(1)  ,&  ! intent(in): derivatives in the aquifer transpiration flux w.r.t. ground temperature
    dAquiferTrans_dCanWat        => deriv_data%var(iLookDERIV%dAquiferTrans_dCanWat       )%dat(1)  ,&  ! intent(in): derivatives in the aquifer transpiration flux w.r.t. canopy total water
    ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
    iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv      )%dat     ,& ! intent(in): [dp(:)]  derivative in vertical liquid water flux at layer interfaces
    ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
    dVolTot_dPsi0                => deriv_data%var(iLookDERIV%dVolTot_dPsi0               )%dat     ,& ! intent(in): [dp(:)]  derivative in total water content w.r.t. total water matric potential
    dCompress_dPsi               => deriv_data%var(iLookDERIV%dCompress_dPsi              )%dat     ,& ! intent(in): [dp(:)]  derivative in compressibility w.r.t matric head
    dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer above
    dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer below
    dq_dHydStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec    )%dat     ,& ! intent(in): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
    ! derivative in baseflow flux w.r.t. aquifer storage
    dBaseflow_dAquifer           => deriv_data%var(iLookDERIV%dBaseflow_dAquifer          )%dat(1)  ,&  ! intent(in): [dp(:)] erivative in baseflow flux w.r.t. aquifer storage (s-1)
    ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
    dq_dNrgStateAbove            => deriv_data%var(iLookDERIV%dq_dNrgStateAbove           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer above
    dq_dNrgStateBelow            => deriv_data%var(iLookDERIV%dq_dNrgStateBelow           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer below
    dq_dNrgStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec    )%dat     ,& ! intent(in): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
    ! derivative in liquid water fluxes for the soil and snow domain w.r.t temperature
    mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk            )%dat     ,& ! intent(in): [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature
    ! derivative in bulk heat capacity w.r.t. relevant state variables
    dVolHtCapBulk_dPsi0          => deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0         )%dat     ,& ! intent(in): [dp(:)]  derivative in bulk heat capacity w.r.t. matric potential
    dVolHtCapBulk_dTheta         => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTheta        )%dat     ,& ! intent(in): [dp(:)]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dCanWat        => deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWat       )%dat(1)  ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dTk            => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTk           )%dat     ,& ! intent(in): [dp(:)]  derivative in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dTkCanopy      => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanopy     )%dat(1)  ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. temperature
    ! derivatives in time
    mLayerdTemp_dt               => deriv_data%var(iLookDERIV%mLayerdTemp_dt              )%dat     ,& ! intent(in): [dp(:)] timestep change in layer temperature
    scalarCanopydTemp_dt         => deriv_data%var(iLookDERIV%scalarCanopydTemp_dt        )%dat(1)  ,& ! intent(in): [dp   ] timestep change in canopy temperature
    ! diagnostic variables
    scalarFracLiqVeg             => diag_data%var(iLookDIAG%scalarFracLiqVeg              )%dat(1)                ,& ! intent(in): [dp]     fraction of liquid water on vegetation (-)
    scalarBulkVolHeatCapVeg      => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg       )%dat(1)  ,& ! intent(in): [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
    mLayerFracLiqSnow            => diag_data%var(iLookDIAG%mLayerFracLiqSnow             )%dat     ,& ! intent(in): [dp(:)]  fraction of liquid water in each snow layer (-)
    mLayerVolHtCapBulk           => diag_data%var(iLookDIAG%mLayerVolHtCapBulk            )%dat     ,& ! intent(in): [dp(:)]  bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
    scalarSoilControl            => diag_data%var(iLookDIAG%scalarSoilControl             )%dat(1)  ,& ! intent(in): [dp]     soil control on infiltration, zero or one
    ! canopy and layer depth
    canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth             )%dat(1)  ,& ! intent(in): [dp   ]  canopy depth (m)
    mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth                   )%dat      & ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ) ! making association with data in structures
    ! --------------------------------------------------------------
    ! initialize error control
    err=0; message='computJacob/'

    ! *********************************************************************************************************************************************************
    ! * PART 0: PRELIMINARIES (INITIALIZE JACOBIAN AND COMPUTE TIME-VARIABLE DIAGONAL TERMS)
    ! *********************************************************************************************************************************************************
   
    ! get the number of state variables
    nState = size(dMat)

    ! initialize the Jacobian
    ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver
    aJac(:,:) = 0._rkind  ! analytical Jacobian matrix

    ! compute terms in the Jacobian for vegetation (excluding fluxes)
    ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
    if(ixVegNrg/=integerMissing)then
      dMat(ixVegNrg) = scalarBulkVolHeatCapVeg + LH_fus*iden_water*dTheta_dTkCanopy &
                       + dVolHtCapBulk_dTkCanopy * scalarCanopydTemp_dt
    endif

    ! compute additional terms for the Jacobian for the snow-soil domain (excluding fluxes)
    ! NOTE: energy for snow+soil is computed *within* the iteration loop as it includes phase change
    do iLayer=1,nLayers
      if(ixSnowSoilNrg(iLayer)/=integerMissing)then
       dMat(ixSnowSoilNrg(iLayer)) = mLayerVolHtCapBulk(iLayer) + LH_fus*iden_water*mLayerdTheta_dTk(iLayer) &
                                     + dVolHtCapBulk_dTk(iLayer) * mLayerdTemp_dt(iLayer)
      endif
    end do

    ! compute additional terms for the Jacobian for the soil domain (excluding fluxes)
    do iLayer=1,nSoil
      if(ixSoilOnlyHyd(iLayer)/=integerMissing)then
       dMat(ixSoilOnlyHyd(iLayer)) = dVolTot_dPsi0(iLayer) + dCompress_dPsi(iLayer)
      endif
    end do

    ! define the form of the matrix
    select case(ixMatrix)
      ! *********************************************************************************************************************************************************
      ! * PART 1: BAND MATRIX
      ! *********************************************************************************************************************************************************
      case(ixBandMatrix)

        ! check
        if(size(aJac,1)/=nBands .or. size(aJac,2)/=size(dMat))then
          message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nBands,nState)'
          err=20; return
        endif

        ! -----
        ! * energy and liquid fluxes over vegetation...
        ! ---------------------------------------------
        if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)

          ! * energy fluxes with the canopy water
          if(ixVegHyd/=integerMissing)then

            ! * cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
            if(ixCasNrg/=integerMissing) aJac(ixOffDiag(ixVegHyd,ixCasNrg),ixCasNrg) = -dCanopyEvaporation_dTCanair*dt
            ! dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy is the derivative in throughfall and canopy drainage with canopy temperature
            if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixVegHyd,ixVegNrg),ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy
            ! * liquid water fluxes for vegetation canopy (-), dt*scalarFracLiqVeg*scalarCanopyLiqDeriv is the derivative in throughfall and canopy drainage with canopy water
                                         aJac(ixDiag,                      ixVegHyd) = -scalarFracLiqVeg*(dCanopyEvaporation_dCanWat - scalarCanopyLiqDeriv)*dt + 1._rkind
            if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixVegHyd,ixTopNrg),ixTopNrg) = -dCanopyEvaporation_dTGround*dt

            ! * cross-derivative terms w.r.t. canopy water (kg-1 m2)
            if(ixTopHyd/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarFracLiqVeg*scalarCanopyLiqDeriv)/iden_water

            ! * cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
            ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fus/canopyDepth = J m-3; dLiq = kg m-2
            if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixVegNrg,ixVegHyd),ixVegHyd) = (-1._rkind + scalarFracLiqVeg)*LH_fus/canopyDepth &
                                                                                      + dVolHtCapBulk_dCanWat * scalarCanopydTemp_dt - (dt/canopyDepth) * dCanopyNetFlux_dCanWat ! dF/dLiq
            if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixTopNrg,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanWat)
          endif

          ! * cross-derivative terms w.r.t. canopy temperature (K-1)
          if(ixVegNrg/=integerMissing)then
            if(ixTopHyd/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water
          endif

          ! * energy fluxes with the canopy air space (J m-3 K-1)
          if(ixCasNrg/=integerMissing)then
                                         aJac(ixDiag,                      ixCasNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + dMat(ixCasNrg)
            if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixCasNrg,ixVegNrg),ixVegNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
            if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixCasNrg,ixTopNrg),ixTopNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)
          endif

          ! * energy fluxes with the vegetation canopy (J m-3 K-1)
          if(ixVegNrg/=integerMissing)then
            if(ixCasNrg/=integerMissing) aJac(ixOffDiag(ixVegNrg,ixCasNrg),ixCasNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
                                         aJac(ixDiag,                      ixVegNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + dMat(ixVegNrg)
            if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixVegNrg,ixTopNrg),ixTopNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)
          endif

          ! * energy fluxes with the surface (J m-3 K-1)
          if(ixTopNrg/=integerMissing)then
            if(ixCasNrg/=integerMissing) aJac(ixOffDiag(ixTopNrg,ixCasNrg),ixCasNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
            if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixTopNrg,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)
          endif

        endif  ! if there is a need to compute energy fluxes within vegetation

        ! -----
        ! * energy fluxes for the snow+soil domain...
        ! -------------------------------------------
        if(nSnowSoilNrg>0)then
          do iLayer=1,nLayers  ! loop through all layers in the snow+soil domain

            ! check if the state is in the subset
            if(ixSnowSoilNrg(iLayer)==integerMissing) cycle

            ! - define index within the state subset and the full state vector
            jState = ixSnowSoilNrg(iLayer)        ! index within the state subset

            ! - diagonal elements
            aJac(ixDiag,jState)   = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jState)

            ! - lower-diagonal elements
            if(iLayer>1)then
              if(ixSnowSoilNrg(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSnowSoilNrg(iLayer-1),jState),jState) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
            endif

            ! - upper diagonal elements
            if(iLayer<nLayers)then
              if(ixSnowSoilNrg(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSnowSoilNrg(iLayer+1),jState),jState) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
            endif

          end do  ! (looping through energy states in the snow+soil domain)
        endif   ! (if the subset includes energy state variables in the snow+soil domain)

        ! -----
        ! * liquid water fluxes for the snow domain...
        ! --------------------------------------------
        if(nSnowOnlyHyd>0)then
          do iLayer=1,nSnow  ! loop through layers in the snow domain

            ! - check that the snow layer is desired
            if(ixSnowOnlyHyd(iLayer)==integerMissing) cycle

            ! - define state indices for the current layer
            watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

            ! compute factor to convert liquid water derivative to total water derivative
            select case( ixHydType(iLayer) )
              case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow(iLayer)
              case default;         convLiq2tot = 1._rkind
            end select

            ! - diagonal elements
            aJac(ixDiag,watState) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot + dMat(watState)

            ! - lower-diagonal elements
            if(iLayer>1)then
              if(ixSnowOnlyHyd(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSnowOnlyHyd(iLayer-1),watState),watState) = 0._rkind  ! sub-diagonal: no dependence on other layers
            endif

            ! - upper diagonal elements
            if(iLayer<nSnow)then
              if(ixSnowOnlyHyd(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSnowOnlyHyd(iLayer+1),watState),watState) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot       ! dVol(below)/dLiq(above) -- (-)
            endif

          end do  ! (looping through liquid water states in the snow domain)
        endif   ! (if the subset includes hydrology state variables in the snow domain)

        ! -----
        ! * cross derivatives in the snow domain...
        ! ----------------------------------------
        if(nSnowOnlyHyd>0 .and. nSnowOnlyNrg>0)then
          do iLayer=1,nSnow  ! loop through layers in the snow domain

            ! - check that the snow layer is desired
            if(ixSnowOnlyNrg(iLayer)==integerMissing) cycle

            ! (define the energy state)
            nrgState = ixSnowOnlyNrg(iLayer)       ! index within the full state vector

            ! - define state indices for the current layer
            watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

            if(watstate/=integerMissing)then       ! (energy state for the current layer is within the state subset)

              ! - include derivatives of energy fluxes w.r.t water fluxes for current layer
              aJac(ixOffDiag(nrgState,watState),watState) = (-1._rkind + mLayerFracLiqSnow(iLayer))*LH_fus*iden_water  &
                                         + dVolHtCapBulk_dTheta(iLayer) * mLayerdTemp_dt(iLayer) &
                                         + (dt/mLayerDepth(iLayer))*(-dNrgFlux_dWatBelow(iLayer-1) + dNrgFlux_dWatAbove(iLayer)) ! (dF/dLiq)

              ! - include derivatives of water fluxes w.r.t energy fluxes for current layer
              aJac(ixOffDiag(watState,nrgState),nrgState) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT)

              ! (cross-derivative terms for the layer below)
              if(iLayer < nSnow)then
                aJac(ixOffDiag(ixSnowOnlyHyd(iLayer+1),nrgState),nrgState) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)        ! dVol(below)/dT(above) -- K-1
              endif ! (if there is a water state in the layer below the current layer in the given state subset)

              ! - include derivatives of heat capacity w.r.t water fluxes for surrounding layers starting with layer above
              if(iLayer>1)then
                if(ixSnowOnlyNrg(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSnowOnlyNrg(iLayer-1),watState),watState) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dWatBelow(iLayer-1) )
              endif

              ! (cross-derivative terms for the layer below)
              if(iLayer<nSnow)then
                if(ixSnowOnlyNrg(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSnowOnlyNrg(iLayer+1),watState),watState) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dWatAbove(iLayer  ) )
              elseif(iLayer==nSnow .and. nSoilOnlyNrg>0)then !bottom snow layer and there is soil below
                if(ixSoilOnlyNrg(1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyNrg(1),watState),watState) = (dt/mLayerDepth(nSnow+1))*(-dNrgFlux_dWatAbove(nSnow) )
              endif

            endif   ! (if the energy state for the current layer is within the state subset)

          end do  ! (looping through snow layers)
        endif   ! (if there are state variables for both water and energy in the snow domain)

        ! -----
        ! * liquid water fluxes for the soil domain...
        ! --------------------------------------------
        if(nSoilOnlyHyd>0)then
          do iLayer=1,nSoil

            ! - check that the soil layer is desired
            if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle

            ! - define state indices
            watState = ixSoilOnlyHyd(iLayer)         ! hydrology state index within the state subset

            ! - define indices of the soil layers
            jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

            ! - compute the diagonal elements
            ! all terms *excluding* baseflow
            aJac(ixDiag,watState) = (dt/mLayerDepth(jLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(watState)

            ! - compute the lower-diagonal elements
            if(iLayer>1)then
              if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(iLayer-1),watState),watState) = (dt/mLayerDepth(jLayer-1))*( dq_dHydStateBelow(iLayer-1))
            endif

            ! - compute the upper-diagonal elements
            if(iLayer<nSoil)then
              if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(iLayer+1),watState),watState) = (dt/mLayerDepth(jLayer+1))*(-dq_dHydStateAbove(iLayer))
            endif

            ! - only include banded terms for baseflow in banded structure
            if(computeBaseflow .and. nSoilOnlyHyd==nSoil)then
              do pLayer=1,nSoil
                qState = ixSoilOnlyHyd(pLayer)  ! hydrology state index within the state subset
                if(qState/=integerMissing)then
                  if((qstate - watState <= ku) .and. (watState - qstate <= kl)) &
                      aJac(ixOffDiag(watState,qState),qState) = (dt/mLayerDepth(jLayer))*dBaseflow_dMatric(iLayer,pLayer) + aJac(ixOffDiag(watState,qState),qState)
                endif
              end do
            endif

            ! - only include banded terms for surface infiltration below surface in banded structure; ixSoilOnlyHyd(1) - watState always <= kl
            if(ixSoilOnlyHyd(1)/=integerMissing)then
              if(watState - ixSoilOnlyHyd(1) <= ku) & 
                aJac(ixOffDiag(ixSoilOnlyHyd(1),watState),watState) = -(dt/mLayerDepth(1+nSnow))*dq_dHydStateLayerSurfVec(iLayer) + aJac(ixOffDiag(ixSoilOnlyHyd(1),watState),watState)
            endif
          end do  ! (looping through hydrology states in the soil domain)

          ! - include terms for surface infiltration above surface
          if(nSnowOnlyHyd>0)then !have snow above first soil layer
            if(ixSnowOnlyHyd(nSnow)/=integerMissing .and. ixSoilOnlyHyd(1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(1),ixSnowOnlyHyd(nSnow)),ixSnowOnlyHyd(nSnow)) = -(dt/mLayerDepth(1+nSnow))*dq_dHydStateLayerSurfVec(0)
          elseif(computeVegFlux)then !have vegetation above first soil layer, ixTopHyd = ixSoilOnlyHyd(1)
            if(ixVegHyd/=integerMissing .and. ixTopHyd/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixVegHyd),ixVegHyd) = -(dt/mLayerDepth(1+nSnow))*dq_dHydStateLayerSurfVec(0) + aJac(ixOffDiag(ixTopHyd,ixVegHyd),ixVegHyd)
          endif

        endif   ! (if the subset includes hydrology state variables in the soil domain)

        ! -----
        ! * liquid water fluxes for the aquifer...
        ! ----------------------------------------
        if(ixAqWat/=integerMissing) then
          aJac(ixDiag,ixAqWat) = -dBaseflow_dAquifer*dt + dMat(ixAqWat)
          if(ixSoilOnlyNrg(nSoil)/=integerMissing) aJac(ixOffDiag(ixAqWat,ixSoilOnlyNrg(nSoil)),ixSoilOnlyNrg(nSoil)) = -dq_dNrgStateAbove(nSoil)*dt ! dAquiferRecharge_dTk  = d_iLayerLiqFluxSoil(nSoil)_dTk
          if(ixSoilOnlyHyd(nSoil)/=integerMissing) aJac(ixOffDiag(ixAqWat,ixSoilOnlyHyd(nSoil)),ixSoilOnlyHyd(nSoil)) = -dq_dHydStateAbove(nSoil)*dt ! dAquiferRecharge_dWat = d_iLayerLiqFluxSoil(nSoil)_dWat
          ! - only include banded terms for derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration), would have to have few soil layers
          if(computeVegFlux)then
            if(ixCasNrg/=integerMissing)then
              if(ixAqWat-ixCasNrg <= kl) aJac(ixOffDiag(ixAqWat,ixCasNrg),ixCasNrg) = -dAquiferTrans_dTCanair*dt ! dVol/dT (K-1)
            endif
            if(ixVegNrg/=integerMissing)then
              if(ixAqWat-ixVegNrg <= kl) aJac(ixOffDiag(ixAqWat,ixVegNrg),ixVegNrg) = -dAquiferTrans_dTCanopy*dt ! dVol/dT (K-1)
            endif
            if(ixVegHyd/=integerMissing)then
              if(ixAqWat-ixVegHyd <= kl) aJac(ixOffDiag(ixAqWat,ixVegHyd),ixVegHyd) = -dAquiferTrans_dCanWat*dt  ! dVol/dLiq (kg m-2)-1
            endif
            if(ixTopNrg/=integerMissing)then
              if(ixAqWat-ixTopNrg <= kl) aJac(ixOffDiag(ixAqWat,ixTopNrg),ixTopNrg) = -dAquiferTrans_dTGround*dt ! dVol/dT (K-1)
            endif
          endif
        endif

        ! -----
        ! * cross derivatives in the soil domain...
        ! ----------------------------------------
        if(nSoilOnlyHyd>0 .and. nSoilOnlyNrg>0)then
          do iLayer=1,nSoilOnlyNrg

            ! - check that the soil layer is desired
            if(ixSoilOnlyNrg(iLayer)==integerMissing) cycle

            ! - define indices of the soil layers
            jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

            ! - define the energy state variable
            nrgState = ixNrgLayer(jLayer)       ! index within the full state vector

            ! - define index of hydrology state variable within the state subset
            watState = ixSoilOnlyHyd(iLayer)

            ! only compute derivatives if the water state for the current layer is within the state subset
            if(watstate/=integerMissing)then

              ! - include derivatives in liquid water fluxes w.r.t. temperature for current layer
              aJac(ixOffDiag(watState,nrgState),nrgState) = (dt/mLayerDepth(jLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance

              ! - compute lower diagonal elements
              if(iLayer>1)then
                if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(iLayer-1),nrgState),nrgState) = (dt/mLayerDepth(jLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
              endif

              ! compute upper-diagonal elements
              if(iLayer<nSoil)then
                if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(iLayer+1),nrgState),nrgState) = (dt/mLayerDepth(jLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1
              endif

              ! - include derivatives w.r.t. ground evaporation
              if(nSnow==0 .and. iLayer==1)then  ! upper-most soil layer, assume here that kl>=4
                if(computeVegFlux)then
                  if(ixCasNrg/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixCasNrg),ixCasNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanair/iden_water) ! dVol/dT (K-1)
                  if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanopy/iden_water) + aJac(ixOffDiag(ixTopHyd,ixVegNrg),ixVegNrg) ! dVol/dT (K-1)
                  if(ixVegHyd/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dCanWat/iden_water)  + aJac(ixOffDiag(ixTopHyd,ixVegHyd),ixVegHyd) ! dVol/dLiq (kg m-2)-1
                endif
                if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixTopNrg),ixTopNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac(ixOffDiag(ixTopHyd,ixTopNrg),ixTopNrg) ! dVol/dT (K-1)
              endif

              ! - only include banded terms for derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration) in banded structure
              if(computeVegFlux)then
                if(ixCasNrg/=integerMissing)then
                  if(watState-ixCasNrg <= kl) aJac(ixOffDiag(watState,ixCasNrg),ixCasNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTCanair(iLayer)) + aJac(ixOffDiag(watState,ixCasNrg),ixCasNrg) ! dVol/dT (K-1)
                endif
                if(ixVegNrg/=integerMissing)then
                  if(watState-ixVegNrg <= kl) aJac(ixOffDiag(watState,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTCanopy(iLayer)) + aJac(ixOffDiag(watState,ixVegNrg),ixVegNrg) ! dVol/dT (K-1)
                endif
                if(ixVegHyd/=integerMissing)then
                  if(watState-ixVegHyd <= kl) aJac(ixOffDiag(watState,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dCanWat(iLayer))  + aJac(ixOffDiag(watState,ixVegHyd),ixVegHyd) ! dVol/dLiq (kg m-2)-1
                endif
                if(ixTopNrg/=integerMissing)then
                  if(watState-ixTopNrg <= kl) aJac(ixOffDiag(watState,ixTopNrg),ixTopNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTGround(iLayer)) + aJac(ixOffDiag(watState,ixTopNrg),ixTopNrg) ! dVol/dT (K-1)
                endif
              endif

              ! - include derivatives in energy fluxes w.r.t. with respect to water for current layer
              aJac(ixOffDiag(nrgState,watState),watState) = dVolHtCapBulk_dPsi0(iLayer) * mLayerdTemp_dt(jLayer) &
                                                        + (dt/mLayerDepth(jLayer))*(-dNrgFlux_dWatBelow(jLayer-1) + dNrgFlux_dWatAbove(jLayer))
              if(mLayerdTheta_dTk(jLayer) > verySmall)then  ! ice is present
                aJac(ixOffDiag(nrgState,watState),watState) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water + aJac(ixOffDiag(nrgState,watState),watState)   ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
              endif

              ! - include derivatives of heat capacity w.r.t water fluxes for surrounding layers starting with layer above
              if(iLayer>1)then
                if(ixSoilOnlyNrg(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyNrg(iLayer-1),watState),watState) = (dt/mLayerDepth(jLayer-1))*( dNrgFlux_dWatBelow(jLayer-1) )
              elseif(iLayer==1 .and. nSnowOnlyNrg>0)then !top soil layer and there is snow above
                if(ixSnowOnlyNrg(nSnow)/=integerMissing) aJac(ixOffDiag(ixSnowOnlyNrg(nSnow),watState),watState) = (dt/mLayerDepth(nSnow))*( dNrgFlux_dWatBelow(nSnow) )
              endif

              ! (cross-derivative terms for the layer below)
              if(iLayer<nSoil)then
                if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyNrg(iLayer+1),watState),watState) = (dt/mLayerDepth(jLayer+1))*(-dNrgFlux_dWatAbove(jLayer  ) )
              endif

            endif   ! (if the water state for the current layer is within the state subset)

            ! - only include banded terms for surface infiltration below surface in banded structure; ixSoilOnlyHyd(1) - nrgState always <= kl
            if(ixSoilOnlyHyd(1)/=integerMissing)then
              if(nrgState - ixSoilOnlyHyd(1) <= ku) &
                  aJac(ixOffDiag(ixSoilOnlyHyd(1),nrgState),nrgState) = -(dt/mLayerDepth(1+nSnow))*dq_dNrgStateLayerSurfVec(iLayer) + aJac(ixOffDiag(ixSoilOnlyHyd(1),nrgState),nrgState)
            endif
          end do  ! (looping through energy states in the soil domain)

          ! - include terms for surface infiltration above surface
          if(nSnowOnlyNrg>0)then !have snow above first soil layer
            if(ixSnowOnlyNrg(nSnow)/=integerMissing .and. ixSoilOnlyHyd(1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(1),ixSnowOnlyNrg(nSnow)),ixSnowOnlyNrg(nSnow)) = -(dt/mLayerDepth(1+nSnow))*dq_dNrgStateLayerSurfVec(0)
          elseif(computeVegFlux)then !have vegetation above first soil layer, ixTopHyd = ixSoilOnlyHyd(1)
            if(ixVegNrg/=integerMissing .and. ixTopHyd/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixVegNrg),ixVegNrg) = -(dt/mLayerDepth(1+nSnow))*dq_dNrgStateLayerSurfVec(0) + aJac(ixOffDiag(ixTopHyd,ixVegNrg),ixVegNrg)
          endif

        endif   ! (if there are state variables for both water and energy in the soil domain)
 
      ! *********************************************************************************************************************************************************
      ! * PART 2: FULL MATRIX
      ! *********************************************************************************************************************************************************
      case(ixFullMatrix)

        ! check
        if(size(aJac,1)/=size(dMat) .or. size(aJac,2)/=size(dMat))then
          message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nState,nState)'
          err=20; return
        endif

        ! -----
        ! * energy and liquid fluxes over vegetation...
        ! ---------------------------------------------
        if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)

          ! * energy fluxes with the canopy water
          if(ixVegHyd/=integerMissing)then

            ! * cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
            if(ixCasNrg/=integerMissing) aJac(ixVegHyd,ixCasNrg) = -dCanopyEvaporation_dTCanair*dt
            ! dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy is the derivative in throughfall and canopy drainage with canopy temperature
            if(ixVegNrg/=integerMissing) aJac(ixVegHyd,ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy
            ! * liquid water fluxes for vegetation canopy (-)
                                         aJac(ixVegHyd,ixVegHyd) = -scalarFracLiqVeg*(dCanopyEvaporation_dCanWat - scalarCanopyLiqDeriv)*dt + 1._rkind
            if(ixTopNrg/=integerMissing) aJac(ixVegHyd,ixTopNrg) = -dCanopyEvaporation_dTGround*dt

            ! cross-derivative terms w.r.t. canopy water (kg-1 m2)
            if(ixTopHyd/=integerMissing) aJac(ixTopHyd,ixVegHyd) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarFracLiqVeg*scalarCanopyLiqDeriv)/iden_water

            ! cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
            ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fus/canopyDepth = J m-3; dLiq = kg m-2
            if(ixVegNrg/=integerMissing) aJac(ixVegNrg,ixVegHyd) = (-1._rkind + scalarFracLiqVeg)*LH_fus/canopyDepth &
                                                                   + dVolHtCapBulk_dCanWat * scalarCanopydTemp_dt - (dt/canopyDepth) * dCanopyNetFlux_dCanWat ! dF/dLiq
            if(ixTopNrg/=integerMissing) aJac(ixTopNrg,ixVegHyd) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanWat)
          endif

          ! * cross-derivative terms w.r.t. canopy temperature (K-1)
          if(ixVegNrg/=integerMissing)then
            if(ixTopHyd/=integerMissing) aJac(ixTopHyd,ixVegNrg) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water
          endif

          ! energy fluxes with the canopy air space (J m-3 K-1)
          if(ixCasNrg/=integerMissing)then
                                         aJac(ixCasNrg,ixCasNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + dMat(ixCasNrg)
            if(ixVegNrg/=integerMissing) aJac(ixCasNrg,ixVegNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
            if(ixTopNrg/=integerMissing) aJac(ixCasNrg,ixTopNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)
          endif

          ! * energy fluxes with the vegetation canopy (J m-3 K-1)
          if(ixVegNrg/=integerMissing)then
            if(ixCasNrg/=integerMissing) aJac(ixVegNrg,ixCasNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
                                         aJac(ixVegNrg,ixVegNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + dMat(ixVegNrg)
            if(ixTopNrg/=integerMissing) aJac(ixVegNrg,ixTopNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)
          endif

          ! energy fluxes with the surface (J m-3 K-1)
          if(ixTopNrg/=integerMissing)then
            if(ixCasNrg/=integerMissing) aJac(ixTopNrg,ixCasNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
            if(ixVegNrg/=integerMissing) aJac(ixTopNrg,ixVegNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)
          endif

         endif  ! if there is a need to compute energy fluxes within vegetation

        ! -----
        ! * energy fluxes for the snow+soil domain...
        ! -------------------------------------------
        if(nSnowSoilNrg>0)then
          do iLayer=1,nLayers  ! loop through all layers in the snow+soil domain

            ! check if the state is in the subset
            if(ixSnowSoilNrg(iLayer)==integerMissing) cycle

            ! - define index within the state subset and the full state vector
            jState = ixSnowSoilNrg(iLayer)        ! index within the state subset

            ! - diagonal elements
            aJac(jState,jState)   = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jState)

            ! - lower-diagonal elements
            if(iLayer>1)then
              if(ixSnowSoilNrg(iLayer-1)/=integerMissing) aJac(ixSnowSoilNrg(iLayer-1),jState) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
            endif

            ! - upper diagonal elements
            if(iLayer<nLayers)then
              if(ixSnowSoilNrg(iLayer+1)/=integerMissing) aJac(ixSnowSoilNrg(iLayer+1),jState) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
            endif

          end do  ! (looping through energy states in the snow+soil domain)
        endif   ! (if the subset includes energy state variables in the snow+soil domain)

        ! -----
        ! * liquid water fluxes for the snow domain...
        ! --------------------------------------------
        if(nSnowOnlyHyd>0)then
          do iLayer=1,nSnow  ! loop through layers in the snow domain

            ! - check that the snow layer is desired
            if(ixSnowOnlyHyd(iLayer)==integerMissing) cycle

            ! - define state indices for the current layer
            watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

            ! compute factor to convert liquid water derivative to total water derivative
            select case( ixHydType(iLayer) )
              case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow(iLayer)
              case default;         convLiq2tot = 1._rkind
            end select

            ! - diagonal elements
            aJac(watState,watState) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot + dMat(watState)

            ! - lower-diagonal elements
            if(iLayer>1)then
              if(ixSnowOnlyHyd(iLayer-1)/=integerMissing) aJac(ixSnowOnlyHyd(iLayer-1),watState) = 0._rkind  ! sub-diagonal: no dependence on other layers
            endif

            ! - upper diagonal elements
            if(iLayer<nSnow)then
              if(ixSnowOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSnowOnlyHyd(iLayer+1),watState) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot       ! dVol(below)/dLiq(above) -- (-)
            endif

           end do  ! (looping through liquid water states in the snow domain)
         endif   ! (if the subset includes hydrology state variables in the snow domain)

        ! -----
        ! * cross derivatives in the snow domain...
        ! ----------------------------------------
        if(nSnowOnlyHyd>0 .and. nSnowOnlyNrg>0)then
          do iLayer=1,nSnow  ! loop through layers in the snow domain

            ! - check that the snow layer is desired
            if(ixSnowOnlyNrg(iLayer)==integerMissing) cycle

            ! (define the energy state)
            nrgState = ixSnowOnlyNrg(iLayer)       ! index within the full state vector

            ! - define state indices for the current layer
            watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

            if(watstate/=integerMissing)then       ! (energy state for the current layer is within the state subset)

             ! - include derivatives of energy fluxes w.r.t water fluxes for current layer
             aJac(nrgState,watState) = (-1._rkind + mLayerFracLiqSnow(iLayer))*LH_fus*iden_water &
                                       + dVolHtCapBulk_dTheta(iLayer) * mLayerdTemp_dt(iLayer) &
                                       + (dt/mLayerDepth(iLayer))*(-dNrgFlux_dWatBelow(iLayer-1) + dNrgFlux_dWatAbove(iLayer)) ! (dF/dLiq)

             ! - include derivatives of water fluxes w.r.t energy fluxes for current layer
             aJac(watState,nrgState) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT)

             ! (cross-derivative terms for the layer below)
             if(iLayer < nSnow)then
               aJac(ixSnowOnlyHyd(iLayer+1),nrgState) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)        ! dVol(below)/dT(above) -- K-1
             endif ! (if there is a water state in the layer below the current layer in the given state subset)

             ! - include derivatives of heat capacity w.r.t water fluxes for surrounding layers starting with layer above
             if(iLayer>1)then
               if(ixSnowOnlyNrg(iLayer-1)/=integerMissing) aJac(ixSnowOnlyNrg(iLayer-1),watState) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dWatBelow(iLayer-1) )
             endif

             ! (cross-derivative terms for the layer below)
             if(iLayer<nSnow)then
               if(ixSnowOnlyNrg(iLayer+1)/=integerMissing) aJac(ixSnowOnlyNrg(iLayer+1),watState) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dWatAbove(iLayer  ) )
             elseif(iLayer==nSnow .and. nSoilOnlyNrg>0)then !bottom snow layer and there is soil below
               if(ixSoilOnlyNrg(1)/=integerMissing) aJac(ixSoilOnlyNrg(1),watState) = (dt/mLayerDepth(nSnow+1))*(-dNrgFlux_dWatAbove(nSnow) )
             endif

           endif   ! (if the energy state for the current layer is within the state subset)

          end do  ! (looping through snow layers)
        endif   ! (if there are state variables for both water and energy in the snow domain)

        ! -----
        ! * liquid water fluxes for the soil domain...
        ! --------------------------------------------
        if(nSoilOnlyHyd>0)then
          do iLayer=1,nSoil

            ! - check that the soil layer is desired
            if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle

            ! - define state indices
            watState = ixSoilOnlyHyd(iLayer)         ! hydrology state index within the state subset

            ! - define indices of the soil layers
            jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

            ! - compute the diagonal elements
            ! all terms *excluding* baseflow
            aJac(watState,watState) = (dt/mLayerDepth(jLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(watState)

            ! - compute the lower-diagonal elements
            if(iLayer>1)then
              if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer-1),watState) = (dt/mLayerDepth(jLayer-1))*( dq_dHydStateBelow(iLayer-1))
            endif

            ! - compute the upper-diagonal elements
            if(iLayer<nSoil)then
              if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer+1),watState) = (dt/mLayerDepth(jLayer+1))*(-dq_dHydStateAbove(iLayer))
            endif

            ! - include terms for baseflow
            if(computeBaseflow .and. nSoilOnlyHyd==nSoil)then
              do pLayer=1,nSoil
                qState = ixSoilOnlyHyd(pLayer)  ! hydrology state index within the state subset
                if(qState/=integerMissing) aJac(watState,qState) = (dt/mLayerDepth(jLayer))*dBaseflow_dMatric(iLayer,pLayer) + aJac(watState,qState)
              end do
            endif

            ! - include terms for surface infiltration below surface
            if(ixSoilOnlyHyd(1)/=integerMissing) aJac(ixSoilOnlyHyd(1),watState) = -(dt/mLayerDepth(1+nSnow))*dq_dHydStateLayerSurfVec(iLayer) + aJac(ixSoilOnlyHyd(1),watState)

          end do  ! (looping through hydrology states in the soil domain)

          ! - include terms for surface infiltration above surface
          if(nSnowOnlyNrg>0)then !have snow above first soil layer
            if(ixSnowOnlyHyd(nSnow)/=integerMissing .and. ixSoilOnlyHyd(1)/=integerMissing) aJac(ixSoilOnlyHyd(1),ixSnowOnlyHyd(nSnow)) = -(dt/mLayerDepth(1+nSnow))*dq_dHydStateLayerSurfVec(0)
          elseif(computeVegFlux)then !have vegetation above first soil layer, ixTopHyd = ixSoilOnlyHyd(1)
            if(ixVegHyd/=integerMissing .and. ixTopHyd/=integerMissing) aJac(ixTopHyd,ixVegHyd) = -(dt/mLayerDepth(1+nSnow))*dq_dHydStateLayerSurfVec(0) + aJac(ixTopHyd,ixVegHyd)
          endif

        endif   ! (if the subset includes hydrology state variables in the soil domain)

        ! -----
        ! * liquid water fluxes for the aquifer...
        ! ----------------------------------------
        if(ixAqWat/=integerMissing) then
          aJac(ixAqWat,ixAqWat) = -dBaseflow_dAquifer*dt + dMat(ixAqWat)
          if(ixSoilOnlyNrg(nSoil)/=integerMissing) aJac(ixAqWat,ixSoilOnlyNrg(nSoil)) = -dq_dNrgStateAbove(nSoil)*dt ! dAquiferRecharge_dTk  = d_iLayerLiqFluxSoil(nSoil)_dTk
          if(ixSoilOnlyHyd(nSoil)/=integerMissing) aJac(ixAqWat,ixSoilOnlyHyd(nSoil)) = -dq_dHydStateAbove(nSoil)*dt ! dAquiferRecharge_dWat = d_iLayerLiqFluxSoil(nSoil)_dWat
          ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
          if(computeVegFlux)then
            if(ixCasNrg/=integerMissing) aJac(ixAqWat,ixCasNrg) = -dAquiferTrans_dTCanair*dt ! dVol/dT (K-1)
            if(ixVegNrg/=integerMissing) aJac(ixAqWat,ixVegNrg) = -dAquiferTrans_dTCanopy*dt ! dVol/dT (K-1)
            if(ixVegHyd/=integerMissing) aJac(ixAqWat,ixVegHyd) = -dAquiferTrans_dCanWat*dt  ! dVol/dLiq (kg m-2)-1
            if(ixTopNrg/=integerMissing) aJac(ixAqWat,ixTopNrg) = -dAquiferTrans_dTGround*dt ! dVol/dT (K-1)
          endif
        endif

        ! -----
        ! * cross derivatives in the soil domain...
        ! ----------------------------------------
        if(nSoilOnlyHyd>0 .and. nSoilOnlyNrg>0)then
          do iLayer=1,nSoilOnlyNrg

            ! - check that the soil layer is desired
            if(ixSoilOnlyNrg(iLayer)==integerMissing) cycle

            ! - define indices of the soil layers
            jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

            ! - define the energy state variable
            nrgState = ixNrgLayer(jLayer)       ! index within the full state vector

            ! - define index of hydrology state variable within the state subset
            watState = ixSoilOnlyHyd(iLayer)

            ! only compute derivatives if the water state for the current layer is within the state subset
            if(watstate/=integerMissing)then

             ! - include derivatives in liquid water fluxes w.r.t. temperature for current layer
             aJac(watState,nrgState) = (dt/mLayerDepth(jLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance

             ! - compute lower diagonal elements
             if(iLayer>1)then
               if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer-1),nrgState) = (dt/mLayerDepth(jLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
             endif

             ! compute upper-diagonal elements
             if(iLayer<nSoil)then
               if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer+1),nrgState) = (dt/mLayerDepth(jLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1
             endif

              ! - include derivatives of energy w.r.t. ground evaporation
              if(nSnow==0 .and. iLayer==1)then  ! upper-most soil layer
                if(computeVegFlux)then
                  if(ixCasNrg/=integerMissing) aJac(ixTopHyd,ixCasNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanair/iden_water) ! dVol/dT (K-1)
                  if(ixVegNrg/=integerMissing) aJac(ixTopHyd,ixVegNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanopy/iden_water) + aJac(ixTopHyd,ixVegNrg) ! dVol/dT (K-1)
                  if(ixVegHyd/=integerMissing) aJac(ixTopHyd,ixVegHyd) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dCanWat/iden_water)  + aJac(ixTopHyd,ixVegHyd) ! dVol/dLiq (kg m-2)-1
                endif
                if(ixTopNrg/=integerMissing) aJac(ixTopHyd,ixTopNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac(ixTopHyd,ixTopNrg) ! dVol/dT (K-1)
              endif

              ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
              if(computeVegFlux)then
                if(ixCasNrg/=integerMissing) aJac(watState,ixCasNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTCanair(iLayer)) + aJac(watState,ixCasNrg) ! dVol/dT (K-1)
                if(ixVegNrg/=integerMissing) aJac(watState,ixVegNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTCanopy(iLayer)) + aJac(watState,ixVegNrg) ! dVol/dT (K-1)
                if(ixVegHyd/=integerMissing) aJac(watState,ixVegHyd) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dCanWat(iLayer))  + aJac(watState,ixVegHyd) ! dVol/dLiq (kg m-2)-1
                if(ixTopNrg/=integerMissing) aJac(watState,ixTopNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTGround(iLayer)) + aJac(watState,ixTopNrg) ! dVol/dT (K-1)
              endif

              ! - include derivatives in energy fluxes w.r.t. with respect to water for current layer
              aJac(nrgState,watState) = dVolHtCapBulk_dPsi0(iLayer) * mLayerdTemp_dt(jLayer) &
                                       + (dt/mLayerDepth(jLayer))*(-dNrgFlux_dWatBelow(jLayer-1) + dNrgFlux_dWatAbove(jLayer))
              if(mLayerdTheta_dTk(jLayer) > verySmall)then  ! ice is present
                aJac(nrgState,watState) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water + aJac(nrgState,watState)   ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
              endif

              ! - include derivatives of heat capacity w.r.t water fluxes for surrounding layers starting with layer above
              if(iLayer>1)then
                if(ixSoilOnlyNrg(iLayer-1)/=integerMissing) aJac(ixSoilOnlyNrg(iLayer-1),watState) = (dt/mLayerDepth(jLayer-1))*( dNrgFlux_dWatBelow(jLayer-1) )
              elseif(iLayer==1 .and. nSnowOnlyNrg>0)then !top soil layer and there is snow above
                if(ixSnowOnlyNrg(nSnow)/=integerMissing) aJac(ixSnowOnlyNrg(nSnow),watState) = (dt/mLayerDepth(nSnow))*( dNrgFlux_dWatBelow(nSnow) )
              endif

              ! (cross-derivative terms for the layer below)
              if(iLayer<nSoil)then
               if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSoilOnlyNrg(iLayer+1),watState) = (dt/mLayerDepth(jLayer+1))*(-dNrgFlux_dWatAbove(jLayer  ) )
              endif

            endif   ! (if the water state for the current layer is within the state subset)

            ! - include terms for surface infiltration below surface
            if(ixSoilOnlyHyd(1)/=integerMissing) aJac(ixSoilOnlyHyd(1),nrgState) = -(dt/mLayerDepth(1+nSnow))*dq_dNrgStateLayerSurfVec(iLayer) + aJac(ixSoilOnlyHyd(1),nrgState)

          end do  ! (looping through soil layers)

          ! - include terms for surface infiltration above surface
          if(nSnowOnlyNrg>0)then !have snow above first soil layer
            if(ixSnowOnlyNrg(nSnow)/=integerMissing .and. ixSoilOnlyHyd(1)/=integerMissing) aJac(ixSoilOnlyHyd(1),ixSnowOnlyNrg(nSnow)) = -(dt/mLayerDepth(1+nSnow))*dq_dNrgStateLayerSurfVec(0)
          elseif(computeVegFlux)then !have vegetation above first soil layer, ixTopHyd = ixSoilOnlyHyd(1)
            if(ixVegNrg/=integerMissing .and. ixTopHyd/=integerMissing) aJac(ixTopHyd,ixVegNrg) = -(dt/mLayerDepth(1+nSnow))*dq_dNrgStateLayerSurfVec(0) + aJac(ixTopHyd,ixVegNrg)
          endif

        endif   ! (if there are state variables for both water and energy in the soil domain)

      ! check
      case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'; return

    end select  ! type of matrix
    ! *********************************************************************************************************************************************************

    ! print the Jacobian
    if(globalPrintFlag)then
      select case(ixMatrix)
        case(ixBandMatrix)
          print*, '** banded analytical Jacobian:'
          write(*,'(a4,1x,100(i17,1x))') 'xCol', (iLayer, iLayer=min(iJac1,nState),min(iJac2,nState))
          do iLayer=kl+1,nBands
            write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJac(iLayer,jLayer),jLayer=min(iJac1,nState),min(iJac2,nState))
          end do
        case(ixFullMatrix)
          print*, '** full analytical Jacobian:'
          write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=min(iJac1,nState),min(iJac2,nState))
          do iLayer=min(iJac1,nState),min(iJac2,nState)
            write(*,'(i4,1x,100(e12.5,1x))') iLayer, aJac(min(iJac1,nState):min(iJac2,nState),iLayer)
          end do
      end select
    endif

    if(any(isNan(aJac)))then
      print *, '******************************* WE FOUND NAN IN JACOBIAN ************************************'
      stop 1
      message=trim(message)//'we found NaN'
      err=20; return
    endif

  ! end association to variables in the data structures
  end associate

end subroutine computJacob

#ifdef SUNDIALS_ACTIVE
! **********************************************************************************************************
! public function computJacob4kinsol: the interface to compute the Jacobian matrix dF/dy + c dF/dy' for IDA solver
! **********************************************************************************************************
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------
integer(c_int) function computJacob4kinsol(sunvec_y, sunvec_r, sunmat_J, &
                            user_data, sunvec_temp1, sunvec_temp2 &
                            ) result(ierr) bind(C, name='computJacob4kinsol')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fsundials_matrix_mod
  use fnvector_serial_mod
  use fsunmatrix_band_mod
  use fsunmatrix_dense_mod
  use type4kinsol 

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)                :: sunvec_y       ! solution N_Vector
  type(N_Vector)                :: sunvec_r       ! residual N_Vector
  type(SUNMatrix)               :: sunmat_J       ! Jacobian SUNMatrix
  type(c_ptr), value            :: user_data      ! user-defined data
  type(N_Vector)                :: sunvec_temp1   ! temporary N_Vector
  type(N_Vector)                :: sunvec_temp2   ! temporary N_Vector

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer       :: Jac(:,:)       ! Jacobian matrix
  type(data4kinsol), pointer    :: eqns_data      ! equations data
! ----------------------------------------------------------------

  ! get equations data from user-defined data
  call c_f_pointer(user_data,eqns_data)

  ! get data arrays from SUNDIALS vectors
  if (eqns_data%ixMatrix==ixBandMatrix) Jac(1:nBands, 1:eqns_data%nState) => FSUNBandMatrix_Data(sunmat_J)
  if (eqns_data%ixMatrix==ixFullMatrix) Jac(1:eqns_data%nState, 1:eqns_data%nState) => FSUNDenseMatrix_Data(sunmat_J)

  ! compute the analytical Jacobian matrix
  ! NOTE: The derivatives were computed in the previous call to computFlux
  !       This occurred either at the call to eval8summa at the start of systemSolv
  !        or in the call to eval8summa in the previous iteration
  call computJacob(&
                ! input: model control
	              eqns_data%dt_cur,                  & ! intent(in):    length of the time step (seconds)
                eqns_data%nSnow,                   & ! intent(in):    number of snow layers
                eqns_data%nSoil,                   & ! intent(in):    number of soil layers
                eqns_data%nLayers,                 & ! intent(in):    number of layers
                eqns_data%computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                (eqns_data%model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel),    & ! intent(in):    flag to indicate if we need to compute baseflow
                eqns_data%ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                ! input: data structures
                eqns_data%indx_data,               & ! intent(in):    index data
                eqns_data%prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                eqns_data%diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                eqns_data%deriv_data,              & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                eqns_data%dBaseflow_dMatric,       & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                ! input-output: Jacobian and its diagonal
                eqns_data%dMat,                    & ! intent(inout): diagonal of the Jacobian matrix
                Jac,                               & ! intent(out):   Jacobian matrix
                ! output: error control
                eqns_data%err,eqns_data%message)     ! intent(out):   error code and error message
  if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
  if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif                                  

  ! return success
  ierr = 0
  return
            
end function computJacob4kinsol
#endif

! **********************************************************************************************************
! private function: get the off-diagonal index in the band-diagonal matrix
! **********************************************************************************************************
function ixOffDiag(jState,iState)
  implicit none
  integer(i4b),intent(in)  :: jState    ! off-diagonal state
  integer(i4b),intent(in)  :: iState    ! diagonal state
  integer(i4b),parameter   :: ixDiag=kl+ku+1   ! index for the diagonal, the offset in the band Jacobian matrix
  integer(i4b)             :: ixOffDiag ! off-diagonal index in the band-diagonal matrix

  ixOffDiag = ixDiag + jState - iState
end function ixOffDiag

end module computJacob_module