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
USE nr_type

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,         & ! data vector with variable length dimension (i4b)
                    var_dlength,         & ! data vector with variable length dimension (rkind)
                    model_options,       & ! defines the model decisions
                    in_type_computJacob, & ! class for computJacob arguments
                    out_type_computJacob   ! class for computJacob arguments

! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables to describe the state variable type
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for layers
USE globalData,only:maxVolIceContent ! snow maximum volumetric ice content to store water (-)

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
private
public::computJacob
public::fluxJacAdd
public::ixInd
#ifdef SUNDIALS_ACTIVE
public::computJacob4kinsol
#endif
contains


! **********************************************************************************************************
! public subroutine computJacob: compute the Jacobian matrix
! **********************************************************************************************************
subroutine computJacob(&
                       ! input: model control
                       in_computJacob,             & ! intent(in):    model control 
                       ! input: data structures
                       indx_data,                  & ! intent(in):    index data
                       prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
                       deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                       dBaseflow_dWat,             & ! intent(in):    derivative in baseflow w.r.t. soil water characteristic
                       dBaseflow_dTk,              & ! intent(in):    derivative in baseflow w.r.t. temperature (m s-1 K-1)
                       ! input-output: Jacobian and its diagonal
                       dMat0,                      & ! intent(in):    diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
                       aJac,                       & ! intent(out):   Jacobian matrix
                       ! output: error control
                       out_computJacob)              ! intent(out):   error code and error message
  ! -----------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  type(in_type_computJacob),intent(in)    :: in_computJacob         ! model control 
  ! input: data structures
  type(var_ilength),intent(in)            :: indx_data              ! indices defining model states and layers
  type(var_dlength),intent(in)            :: prog_data              ! prognostic variables for a local HRU
  type(var_dlength),intent(in)            :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)            :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(in)                  :: dBaseflow_dWat(:,:)    ! derivative in baseflow w.r.t. soil water characteristic
  real(rkind),intent(in)                  :: dBaseflow_dTk(:,:)     ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  ! input-output: Jacobian and its diagonal
  real(rkind),intent(in)                  :: dMat0(:)               ! diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
  real(rkind),intent(out)                 :: aJac(:,:)              ! Jacobian matrix
  ! output variables
  type(out_type_computJacob),intent(out)  :: out_computJacob        ! error control
  ! -------------------------------------------------------------
  ! * local variables
  ! -------------------------------------------------------------
  real(rkind),allocatable                 :: dMat(:)                ! diagonal of the Jacobian matrix excluding fluxes, depending on the state vector
  ! indices of model state variables
  integer(i4b)                            :: nrgState               ! energy state variable
  integer(i4b)                            :: watState               ! hydrology state variable
  integer(i4b)                            :: nState                 ! number of state variables
  ! indices of model layers
  integer(i4b)                            :: iLayer                 ! index of model layer
  integer(i4b)                            :: jLayer                 ! index of model layer within the full state vector (hydrology)
  integer(i4b)                            :: qLayer                 ! indices of snow+lake+glce layers
  ! conversion factors
  real(rkind)                             :: maxVolIceContent_use   ! maximum volumetric ice content depending if snow or firn
  character(LEN=256)                      :: cmessage               ! error message of downwind routine
  logical(lgt)                            :: full                   ! flag to indicate if the matrix is full (true) or banded (false)
  ! --------------------------------------------------------------
  ! associate variables from data structures
  associate(&
    ! model control
    dt                           => in_computJacob % dt                                        ,& ! intent(in): length of the time step (seconds)
    nSnow                        => in_computJacob % nSnow                                     ,& ! intent(in): number of snow layers
    nLake                        => in_computJacob % nLake                                     ,& ! intent(in): number of lake layers
    nSoil                        => in_computJacob % nSoil                                     ,& ! intent(in): number of soil layers
    nGlce                        => in_computJacob % nGlce                                     ,& ! intent(in): number of glacier ice layers
    nLayers                      => in_computJacob % nLayers                                   ,& ! intent(in): total number of layers in the layer domains
    computeVegFlux               => in_computJacob % computeVegFlux                            ,& ! intent(in): flag to indicate if computing fluxes over vegetation
    computeBaseflow              => in_computJacob % computeBaseflow                           ,& ! intent(in): flag to indicate if computing baseflow
    ixMatrix                     => in_computJacob % ixMatrix                                  ,& ! intent(in): form of the Jacobian matrix
    ! indices of model state variables
    ixCasNrg                     => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)                  ,& ! intent(in): [i4b] index of canopy air space energy state variable
    ixVegNrg                     => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)                  ,& ! intent(in): [i4b] index of canopy energy state variable
    ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                  ,& ! intent(in): [i4b] index of canopy hydrology state variable (mass)
    ! vector of energy indices for the layer domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnLaSoGlNrg                => indx_data%var(iLookINDEX%ixSnLaSoGlNrg)%dat                ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the layer domains
    ixSoilOnlyNrg                => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat                ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
    ! vector of hydrology indices for the layer domains
    ! NOTE: states not in the subset are equal to integerMissing
    noThetaChange                => indx_data%var(iLookINDEX%noThetaChange)%dat(1)             ,& ! intent(in): [i4b]    number of layers with no change in total water content (bottom layers)  
    ixSnLaSoGlHyd                => indx_data%var(iLookINDEX%ixSnLaSoGlHyd)%dat                ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the layer domains
    ixSoilOnlyHyd                => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat                ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! number of state variables of a specific type
    nSnLaSoGlNrg                 => indx_data%var(iLookINDEX%nSnLaSoGlNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the layer domains
    nSnowOnlyNrg                 => indx_data%var(iLookINDEX%nSnowOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
    nSoilOnlyNrg                 => indx_data%var(iLookINDEX%nSoilOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
    nLakeOnlyNrg                 => indx_data%var(iLookINDEX%nLakeOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the lake domain
    nGlceOnlyNrg                 => indx_data%var(iLookINDEX%nGlceOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the glacier ice domain
    nSnLaSoGlHyd                 => indx_data%var(iLookINDEX%nSnLaSoGlHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the layer domains
    nSnowOnlyHyd                 => indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    nLakeOnlyHyd                 => indx_data%var(iLookINDEX%nLakeOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the lake domain
    nGlceOnlyHyd                 => indx_data%var(iLookINDEX%nGlceOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the glacier ice domain
    ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
    dCanopyNetFlux_dCanWat       => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanWat)%dat(1)   ,& ! intent(in): [dp]     derivative in net canopy fluxes w.r.t. canopy total water content
    ! derivatives in canopy water w.r.t canopy temperature
    dTheta_dTkCanopy             => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)         ,& ! intent(in): [dp]     derivative in volumetric liquid water content w.r.t. temperature
    dFracLiqVeg_dTkCanopy        => deriv_data%var(iLookDERIV%dFracLiqVeg_dTkCanopy)%dat(1)    ,& ! intent(in): [dp]     derivative in fraction of (throughfall + drainage)  w.r.t. temperature
    ! derivatives in energy fluxes at the interface of layers w.r.t. water state in layers above and below
    dNrgFlux_dWatAbove           => deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove)%dat          ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer above
    dNrgFlux_dWatBelow           => deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow)%dat          ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer below
    ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
    dVolTot_dPsi0                => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat               ,& ! intent(in): [dp(:)]  derivatives in total water content w.r.t. total water matric potential
    dCompress_dPsi               => deriv_data%var(iLookDERIV%dCompress_dPsi)%dat              ,& ! intent(in): [dp(:)]  derivatives in compressibility w.r.t matric head
    ! derivative in liquid water fluxes for the soil and snow domain w.r.t temperature
    dFracLiqWat_dTk              => deriv_data%var(iLookDERIV%dFracLiqWat_dTk)%dat             ,& ! intent(in): [dp(:)]  derivatives in fraction of liquid w.r.t. temperature
    mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat            ,& ! intent(in): [dp(:)]  derivatives in volumetric liquid water content w.r.t. temperature
    ! derivative in bulk heat capacity w.r.t. relevant state variables
    dVolHtCapBulk_dPsi0          => deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0)%dat         ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. matric potential
    dVolHtCapBulk_dTheta         => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTheta)%dat        ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dCanWat        => deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWat)%dat(1)    ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dTk            => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTk)%dat           ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dTkCanopy      => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanopy)%dat(1)  ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. temperature
    ! derivative in Cm w.r.t. relevant state variables
    dCm_dPsi0                    => deriv_data%var(iLookDERIV%dCm_dPsi0)%dat                   ,& ! intent(in): [dp(:)]  derivatives in heat capacity w.r.t. matric potential (J kg-1)
    dCm_dTk                      => deriv_data%var(iLookDERIV%dCm_dTk)%dat                     ,& ! intent(in): [dp(:)]  derivatives in heat capacity w.r.t. temperature (J kg-1 K-2)
    dCm_dTkCanopy                => deriv_data%var(iLookDERIV%dCm_dTkCanopy)%dat(1)            ,& ! intent(in): [dp   ]  derivative in heat capacity w.r.t. canopy temperature (J kg-1 K-2)
    ! derivatives in time
    mLayerdTemp_dt               => deriv_data%var(iLookDERIV%mLayerdTemp_dt)%dat             ,& ! intent(in):  [dp(:)] timestep change in layer temperature
    scalarCanopydTemp_dt         => deriv_data%var(iLookDERIV%scalarCanopydTemp_dt)%dat(1)    ,& ! intent(in):  [dp   ] timestep change in canopy temperature
    mLayerdWat_dt                => deriv_data%var(iLookDERIV%mLayerdWat_dt)%dat              ,& ! intent(in):  [dp(:)] timestep change in layer volumetric fraction of total water
    scalarCanopydWat_dt          => deriv_data%var(iLookDERIV%scalarCanopydWat_dt)%dat(1)     ,& ! intent(in):  [dp   ] timestep change in canopy total water
    ! diagnostic variables
    scalarFracLiqVeg             => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)          ,& ! intent(in): [dp]     fraction of liquid water on vegetation (-)
    scalarBulkVolHeatCapVeg      => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)   ,& ! intent(in): [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
    scalarCanopyCm               => diag_data%var(iLookDIAG%scalarCanopyCm)%dat(1)            ,& ! intent(in): [dp]     Cm for canopy vegetation (J kg-1)
    mLayerFracLiq                => diag_data%var(iLookDIAG%mLayerFracLiq)%dat                ,& ! intent(in): [dp(:)]  fraction of liquid water in each snow, lake, or glce layer (-)
    mLayerVolHtCapBulk           => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat           ,& ! intent(in): [dp(:)]  bulk volumetric heat capacity in each layer (J m-3 K-1)
    mLayerCm                     => diag_data%var(iLookDIAG%mLayerCm)%dat                     ,& ! intent(in): [dp(:)]  Cm for each layer (J m-3)
    ! canopy and layer depth
    canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)         ,& ! intent(in): [dp   ]  canopy depth (m)
    mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat                  ,& ! intent(in): [dp(:)]  depth of each layer in the sub-domain (m)
    ! output variables
    err                         => out_computJacob % err                                     ,& ! error code
    message                     => out_computJacob % cmessage                                 & ! error message
    ) ! making association with data in structures
    ! --------------------------------------------------------------
    ! initialize error control
    err=0; message='computJacob/'

    ! *********************************************************************************************************************************************************
    ! * PART 0: PRELIMINARIES (INITIALIZE JACOBIAN AND COMPUTE TIME-VARIABLE DIAGONAL TERMS)
    ! *********************************************************************************************************************************************************
    ! get the number of state variables
    nState = size(dMat0)

    ! initialize the Jacobian and diagonal
    ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver and dMat is modified below
    aJac(:,:) = 0._rkind  ! analytical Jacobian matrix
    allocate(dMat(nState)) 
    dMat = dMat0 ! dMat0(ixCasNrg) = Cp_air*iden_air and dMat0(Wat states) = 1.0

    if(computeVegFlux)then
      ! compute terms in the Jacobian for vegetation (excluding fluxes)
      if(ixVegNrg/=integerMissing)&
        dMat(ixVegNrg) = scalarBulkVolHeatCapVeg + LH_fus*iden_water*dTheta_dTkCanopy &
                         + dVolHtCapBulk_dTkCanopy * scalarCanopydTemp_dt &
                         + dCm_dTkCanopy * scalarCanopydWat_dt/canopyDepth &
                         + LH_fus * dFracLiqVeg_dTkCanopy * scalarCanopydWat_dt/canopyDepth
    endif

    ! compute terms for the Jacobian for the layer domain (excluding fluxes)
    do iLayer=1,nLayers
      if(ixSnLaSoGlNrg(iLayer)/=integerMissing)&
          dMat(ixSnLaSoGlNrg(iLayer)) = mLayerVolHtCapBulk(iLayer) + LH_fus*iden_water*mLayerdTheta_dTk(iLayer) &
                                        + dVolHtCapBulk_dTk(iLayer) * mLayerdTemp_dt(iLayer) &
                                        + dCm_dTk(iLayer) * mLayerdWat_dt(iLayer) &
                                        + LH_fus * iden_water * dFracLiqWat_dTk(iLayer) * mLayerdWat_dt(iLayer)
    end do

    ! compute terms for the Jacobian for the soil domain (excluding fluxes)
    do iLayer=1,nSoil
      if(ixSoilOnlyHyd(iLayer)/=integerMissing)& ! writes over dMat(ixSoilOnlyHyd(iLayer) = 1.0
          dMat(ixSoilOnlyHyd(iLayer)) = dVolTot_dPsi0(iLayer) + dCompress_dPsi(iLayer)
    end do

    ! compute the maximum volumetric ice content for the layer domains
    if(nGlce>0)then ! snow can be firn
      maxVolIceContent_use = min(maxVolIceContent+0.15,0.85_rkind) ! firn maximum volumetric ice content to store water (-)
    else ! snow
      maxVolIceContent_use = maxVolIceContent ! snow maximum volumetric ice content to store water (-)
    endif

    ! define the form of the matrix
    select case(ixMatrix)
     case(ixBandMatrix)
       ! check
       if(size(aJac,1)/=nBands .or. size(aJac,2)/=size(dMat))then
         message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nBands,nState)'
         err=20; return
       endif
       full = .false.
     case(ixFullMatrix)
       ! check
       if(size(aJac,1)/=size(dMat) .or. size(aJac,2)/=size(dMat))then
         message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nState,nState)'
         err=20; return
       endif
       full = .true.
     case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'; return
    end select

    ! *********************************************************************************************************************************************************
    ! * PART 1: COMPUTE CROSS-DERIVATIVE JACOBIAN TERMS 
    ! *********************************************************************************************************************************************************
    ! -----
    ! * cross derivatives in the vegetation...
    ! ---------------------------------------------
    if(computeVegFlux)then ! (derivatives only defined when vegetation protrudes over the surface)
      if(ixVegHyd/=integerMissing .and. ixVegNrg/=integerMissing)&
          ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fus/canopyDepth = J m-3; dLiq = kg m-2
          aJac(ixInd(full,ixVegNrg,ixVegHyd),ixVegHyd) = (-1._rkind + scalarFracLiqVeg)*LH_fus/canopyDepth &
                                                     + dVolHtCapBulk_dCanWat * scalarCanopydTemp_dt + scalarCanopyCm/canopyDepth &
                                                     - (dt/canopyDepth) * dCanopyNetFlux_dCanWat &
                                                     + LH_fus * scalarCanopydTemp_dt * dFracLiqVeg_dTkCanopy/canopyDepth
    endif  ! if there is a need to compute energy fluxes within vegetation

    ! -----
    ! * cross derivatives in the snow, lake, glce domains...
    ! ----------------------------------------
    if((nSnowOnlyHyd>0 .and. nSnowOnlyNrg>0) .or. (nLakeOnlyHyd>0 .and. nLakeOnlyNrg>0) .or. (nGlceOnlyHyd>0 .and. nGlceOnlyNrg>0))then
      do qLayer=1,nSnow+nLake+nGlce-noThetaChange ! loop through layers in the snow, lake, glce domains

        if(qLayer<=nSnow+nLake)then
          jLayer = qLayer
          iLayer = qLayer
        else
          jLayer = qLayer + nSoil
          iLayer = qLayer - nSnow - nLake
        endif
        ! - check that the layer is desired
        if(ixSnLaSoGlNrg(jLayer)==integerMissing) cycle
        ! (define the energy state)
        nrgState = ixSnLaSoGlNrg(jLayer)       ! index within the full state vector
        ! - define state indices for the current layer
        watState = ixSnLaSoGlHyd(jLayer)   ! hydrology state index within the state subset

        if(watState/=integerMissing)then       ! (water state for the current layer is within the state subset)
          ! - include derivatives of energy fluxes w.r.t water fluxes for current layer
          aJac(ixInd(full,nrgState,watState),watState) = (-1._rkind + mLayerFracLiq(jLayer))*LH_fus*iden_water  &
                                     + dVolHtCapBulk_dTheta(jLayer) * mLayerdTemp_dt(jLayer) + mLayerCm(jLayer) &
                                     + (dt/mLayerDepth(jLayer))*(-dNrgFlux_dWatBelow(jLayer-1) + dNrgFlux_dWatAbove(jLayer)) &
                                     + LH_fus*iden_water * mLayerdTemp_dt(jLayer) * dFracLiqWat_dTk(jLayer)    ! (dF/dLiq)
        endif ! (if the water state for the current layer is within the state subset)

      end do ! (looping through snow, lake, glce layers)
    endif ! (if there are state variables for both water and energy in the snow, lake, glce domains)

    ! -----
    ! * cross derivatives in the soil domain...
    ! ----------------------------------------
    if(nSoilOnlyHyd>0 .and. nSoilOnlyNrg>0)then
      do iLayer=1,nSoilOnlyNrg

        ! - check that the soil layer is desired
        if(ixSoilOnlyNrg(iLayer)==integerMissing) cycle
        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow+nLake            ! index of layer in the layer system
        ! - define the energy state variable
        nrgState = ixSoilOnlyNrg(iLayer)         ! index within the full state vector
        ! - define index of hydrology state variable within the state subset
        watState = ixSoilOnlyHyd(iLayer)

        ! only compute derivatives if the water state for the current layer is within the state subset
        if(watState/=integerMissing)then
          ! - include derivatives in energy fluxes w.r.t. with respect to water for current layer
          aJac(ixInd(full,nrgState,watState),watState) = (dt/mLayerDepth(jLayer))*(-dNrgFlux_dWatBelow(jLayer-1) + dNrgFlux_dWatAbove(jLayer))
          aJac(ixInd(full,nrgState,watState),watState) = mLayerCm(jLayer) * dVolTot_dPsi0(iLayer) + dVolHtCapBulk_dPsi0(iLayer) * mLayerdTemp_dt(jLayer) &
                                                        + dCm_dPsi0(iLayer) * mLayerdWat_dt(jLayer) + aJac(ixInd(full,nrgState,watState),watState)
          if(mLayerdTheta_dTk(jLayer) > tiny(1.0_rkind))& ! ice is present
              aJac(ixInd(full,nrgState,watState),watState) = -LH_fus*iden_water * dVolTot_dPsi0(iLayer) + aJac(ixInd(full,nrgState,watState),watState)   ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
        endif ! (if the water state for the current layer is within the state subset)

      end do ! (looping through energy states in the soil domain)
    endif ! (if there are state variables for both water and energy in the soil domain)

    ! *********************************************************************************************************************************************************
    ! * PART 2: COMPUTE FLUX JACOBIAN TERMS 
    ! *********************************************************************************************************************************************************
    call fluxJacAdd(full,maxVolIceContent_use,dt,nSnow,nLake,nSoil,nGlce,nLayers,computeVegFlux,computeBaseflow,&
                    indx_data,prog_data,diag_data,deriv_data,dBaseflow_dWat,dBaseflow_dTk,dMat,aJac,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    deallocate(dMat)

    ! *********************************************************************************************************************************************************
    ! * PART 3: JACOBIAN PRINT (IF DESIRED)
    ! *********************************************************************************************************************************************************
    if(globalPrintFlag .or. any(isNan(aJac)))then
      if(full) then
        print*, '** full analytical Jacobian:'
        write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=min(iJac1,nState),min(iJac2,nState))
        do iLayer=min(iJac1,nState),min(iJac2,nState)
          write(*,'(i4,1x,100(e12.5,1x))') iLayer, aJac(min(iJac1,nState):min(iJac2,nState),iLayer)
        end do
      else
        print*, '** banded analytical Jacobian:'
        write(*,'(a4,1x,100(i17,1x))') 'xCol', (iLayer, iLayer=min(iJac1,nState),min(iJac2,nState))
        do iLayer=kl+1,nBands
          write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJac(iLayer,jLayer),jLayer=min(iJac1,nState),min(iJac2,nState))
        end do
      endif
    endif
    if(any(isNan(aJac)))then; message=trim(message)//'NaN in Jacobian';err=20; return; endif

  end associate ! end association to variables in the data structures  

end subroutine computJacob

! ***********************************************************************************************************
! public subroutine to compute flux parts of the Jacobian that are shared between IDA and BE
! ***********************************************************************************************************
subroutine fluxJacAdd(&
                      ! input: model control
                      full,                       & ! intent(in):    flag to indicate if the matrix is full (true) or banded (false)
                      maxVolIceContent_use,       & ! intent(in):    maximum volumetric ice content depending if snow or firn
                      dt,                         & ! intent(in):    length of the time step (seconds)
                      nSnow,                      & ! intent(in):    number of snow layers
                      nLake,                      & ! intent(in):    number of lake layers
                      nSoil,                      & ! intent(in):    number of soil layers
                      nGlce,                      & ! intent(in):    number of glacier ice layers
                      nLayers,                    & ! intent(in):    total number of layers
                      computeVegFlux,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
                      ! input: data structures
                      indx_data,                  & ! intent(in):    index data
                      prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
                      deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                      dBaseflow_dWat,             & ! intent(in):    derivative in baseflow w.r.t. soil water characteristic
                      dBaseflow_dTk,              & ! intent(in):    derivative in baseflow w.r.t. temperature (m s-1 K-1)
                      ! input-output: Jacobian and its diagonal
                      dMat,                       & ! intent(in):    diagonal of the Jacobian matrix
                      aJac,                       & ! intent(inout): Jacobian matrix with flux terms added
                      ! output: error control
                      err,message)                  ! intent(out):   error code and error message
  ! -----------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  logical(lgt),intent(in)              :: full                       ! flag to indicate if the matrix is full (true) or banded (false)
  real(rkind),intent(in)               :: maxVolIceContent_use       ! maximum volumetric ice content depending if snow or firn
  real(rkind),intent(in)               :: dt                         ! length of the time step (seconds)
  integer(i4b),intent(in)              :: nSnow                      ! number of snow layers
  integer(i4b),intent(in)              :: nLake                      ! number of lake layers
  integer(i4b),intent(in)              :: nSoil                      ! number of soil layers
  integer(i4b),intent(in)              :: nGlce                      ! number of glacier ice layers
  integer(i4b),intent(in)              :: nLayers                    ! total number of layers in the layer domains
  logical(lgt),intent(in)              :: computeVegFlux             ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)              :: computeBaseflow            ! flag to indicate if computing baseflow
  ! input: data structures
  type(var_ilength),intent(in)         :: indx_data                  ! indices defining model states and layers
  type(var_dlength),intent(in)         :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(in)         :: diag_data                  ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)         :: deriv_data                 ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(in)               :: dBaseflow_dWat(:,:)        ! derivative in baseflow w.r.t. soil water characteristic
  real(rkind),intent(in)               :: dBaseflow_dTk(:,:)         ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  ! input-output: Jacobian and its diagonal
  real(rkind),intent(in)               :: dMat(:)                    ! diagonal of the Jacobian matrix
  real(rkind),intent(inout)            :: aJac(:,:)                  ! Jacobian matrix with flux terms added
  ! output variables
  integer(i4b),intent(out)             :: err                        ! error code
  character(*),intent(out)             :: message                    ! error message
  ! --------------------------------------------------------------
  ! * local variables
  ! --------------------------------------------------------------
  ! indices of model state variables
  integer(i4b)                         :: qState                     ! index of cross-derivative state variable for baseflow
  integer(i4b)                         :: nrgState                   ! energy state variable
  integer(i4b)                         :: watState                   ! hydrology state variable
  ! indices of model layers
  integer(i4b)                         :: iLayer,pLayer              ! index of model layer
  integer(i4b)                         :: jLayer                     ! index of model layer within the full state vector (hydrology)
  integer(i4b)                         :: qLayer                     ! indices of snow+lake+glce layers
  integer(i4b)                         :: endLayerWat                ! index of the last layer in a water domain
  integer(i4b)                         :: endLayerNrg                ! index of the last layer in an energy domain
  integer(i4b)                         :: denseLimit                 ! index of the limiting dense layer
  logical(i4b)                         :: solid                      ! flag to indicate if layer is solid ice (frozen lake or glacier ice)
  ! conversion factors
  real(rkind)                          :: convLiq2tot                ! factor to convert liquid water derivative to total water derivative
  ! --------------------------------------------------------------
  ! associate variables from data structures
  associate(&
    ! indices of model state variables
    ixCasNrg                     => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)                      ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                     => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)                      ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                      ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixTopNrg                     => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)                      ,& ! intent(in): [i4b]    index of upper-most energy state in the layer domains
    ixTopHyd                     => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)                      ,& ! intent(in): [i4b]    index of upper-most hydrology state in the layer domains
    ixAqWat                      => indx_data%var(iLookINDEX%ixAqWat)%dat(1)                       ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ! vector of energy indices for the layer domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnLaSoGlNrg                => indx_data%var(iLookINDEX%ixSnLaSoGlNrg)%dat                    ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the layer domains
    ixSoilOnlyNrg                => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat                    ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
    ! vector of hydrology indices for the layer domains
    ! NOTE: states not in the subset are equal to integerMissing
    noThetaChange                => indx_data%var(iLookINDEX%noThetaChange)%dat(1)                 ,& ! intent(in): [i4b]    number of layers with no change in total water content (bottom layers)  
    ixSnLaSoGlHyd                => indx_data%var(iLookINDEX%ixSnLaSoGlHyd)%dat                    ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the layer domains
    ixSoilOnlyHyd                => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat                    ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! number of state variables of a specific type
    nSnLaSoGlNrg                 => indx_data%var(iLookINDEX%nSnLaSoGlNrg)%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the layer domains
    nSnowOnlyNrg                 => indx_data%var(iLookINDEX%nSnowOnlyNrg)%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
    nSoilOnlyNrg                 => indx_data%var(iLookINDEX%nSoilOnlyNrg)%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
    nLakeOnlyNrg                 => indx_data%var(iLookINDEX%nLakeOnlyNrg)%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the lake domain
    nGlceOnlyNrg                 => indx_data%var(iLookINDEX%nGlceOnlyNrg)%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the glacier ice domain
    nSnLaSoGlHyd                 => indx_data%var(iLookINDEX%nSnLaSoGlHyd)%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the layer domains
    nSnowOnlyHyd                 => indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    nLakeOnlyHyd                 => indx_data%var(iLookINDEX%nLakeOnlyHyd)%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the lake domain
    nGlceOnlyHyd                 => indx_data%var(iLookINDEX%nGlceOnlyHyd)%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the glacier ice domain
    ! type and index of model control volume
    ixHydType                    => indx_data%var(iLookINDEX%ixHydType)%dat                        ,& ! intent(in): [i4b(:)] index of the type of hydrology states in layer domains
    ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
    dCanairNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanairTemp)%dat(1)   ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. canopy air temperature
    dCanairNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanopyTemp)%dat(1)   ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. canopy temperature
    dCanairNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dGroundTemp)%dat(1)   ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. ground temperature
    dCanopyNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanairTemp)%dat(1)   ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. canopy air temperature
    dCanopyNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanopyTemp)%dat(1)   ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. canopy temperature
    dCanopyNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dGroundTemp)%dat(1)   ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. ground temperature
    dGroundNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanairTemp)%dat(1)   ,& ! intent(in): [dp]     derivative in net ground flux w.r.t. canopy air temperature
    dGroundNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanopyTemp)%dat(1)   ,& ! intent(in): [dp]     derivative in net ground flux w.r.t. canopy temperature
    dGroundNetFlux_dCanWat       => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanWat)%dat(1)       ,& ! intent(in): [dp]     derivative in net ground fluxes w.r.t. canopy total water content
    ! derivatives in evaporative fluxes w.r.t. relevant state variables
    dCanopyEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanair)%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy air temperature
    dCanopyEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanopy)%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy temperature
    dCanopyEvaporation_dTGround  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTGround)%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. ground temperature
    dCanopyEvaporation_dCanWat   => deriv_data%var(iLookDERIV%dCanopyEvaporation_dCanWat)%dat(1)   ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy total water content
    dGroundEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanair)%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy air temperature
    dGroundEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanopy)%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy temperature
    dGroundEvaporation_dTGround  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround)%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. ground temperature
    dGroundEvaporation_dCanWat   => deriv_data%var(iLookDERIV%dGroundEvaporation_dCanWat)%dat(1)   ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy total water content
    ! derivatives in canopy water w.r.t canopy temperature
    dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy)%dat(1)             ,& ! intent(in): [dp]     derivative in canopy liquid storage w.r.t. temperature
    ! derivatives in canopy liquid fluxes w.r.t. canopy water
    scalarCanopyLiqDeriv         => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv)%dat(1)         ,& ! intent(in): [dp]     derivative in (throughfall + drainage) w.r.t. canopy liquid water
    ! derivatives in energy fluxes at the interface of layers w.r.t. temperature in layers above and below
    dNrgFlux_dTempAbove          => deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove)%dat             ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer above
    dNrgFlux_dTempBelow          => deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow)%dat             ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer below
    ! derivatives in energy fluxes at the interface of layers w.r.t. water state in layers above and below
    dNrgFlux_dWatAbove           => deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove)%dat              ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer above
    dNrgFlux_dWatBelow           => deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow)%dat              ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer below
    ! derivatives in soil transpiration w.r.t. canopy state variables
    mLayerdTrans_dTCanair        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanair)%dat           ,& ! intent(in): [dp(:)]  derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
    mLayerdTrans_dTCanopy        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanopy)%dat           ,& ! intent(in): [dp(:)]  derivatives in the soil layer transpiration flux w.r.t. canopy temperature
    mLayerdTrans_dTGround        => deriv_data%var(iLookDERIV%mLayerdTrans_dTGround)%dat           ,& ! intent(in): [dp(:)]  derivatives in the soil layer transpiration flux w.r.t. ground temperature
    mLayerdTrans_dCanWat         => deriv_data%var(iLookDERIV%mLayerdTrans_dCanWat)%dat            ,& ! intent(in): [dp(:)]  derivatives in the soil layer transpiration flux w.r.t. canopy total water
    ! derivatives in aquifer transpiration w.r.t. canopy state variables
    dAquiferTrans_dTCanair       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanair)%dat(1)       ,& ! intent(in): [dp]     derivative in the aquifer transpiration flux w.r.t. canopy air temperature
    dAquiferTrans_dTCanopy       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanopy)%dat(1)       ,& ! intent(in): [dp]     derivative in the aquifer transpiration flux w.r.t. canopy temperature
    dAquiferTrans_dTGround       => deriv_data%var(iLookDERIV%dAquiferTrans_dTGround)%dat(1)       ,& ! intent(in): [dp]     derivative in the aquifer transpiration flux w.r.t. ground temperature
    dAquiferTrans_dCanWat        => deriv_data%var(iLookDERIV%dAquiferTrans_dCanWat)%dat(1)        ,& ! intent(in): [dp]     derivative in the aquifer transpiration flux w.r.t. canopy total water
    ! derivative in liquid water fluxes at the interface of snow, lake, glce layers w.r.t. volumetric liquid water content in the layer above
    iLayerLiqFluxSnLaGlDeriv     => deriv_data%var(iLookDERIV%iLayerLiqFluxSnLaGlDeriv)%dat        ,& ! intent(in): [dp(:)]  derivative in vertical liquid water flux at layer interfaces
    scalarSurfaceIceMeltDeriv    => deriv_data%var(iLookDERIV%scalarSurfaceIceMeltDeriv)%dat(1)    ,& ! intent(in): [dp]     derivative in ice melt flux at top interface
    ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
    dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove)%dat               ,& ! intent(in): [dp(:)]  derivatives in  flux at layer interfaces w.r.t. states in the layer above
    dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow)%dat               ,& ! intent(in): [dp(:)]  derivatives in  flux at layer interfaces w.r.t. states in the layer below
    dq_dHydStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec)%dat        ,& ! intent(in): [dp(:)]  derivatives in  the flux in soil surface interface w.r.t. state variables in layers
     ! derivative in baseflow flux w.r.t. aquifer storage
    dBaseflow_dAquifer           => deriv_data%var(iLookDERIV%dBaseflow_dAquifer)%dat(1)           ,& ! intent(in): [dp(:)]  derivative in baseflow flux w.r.t. aquifer storage (s-1)
    ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
    dq_dNrgStateAbove            => deriv_data%var(iLookDERIV%dq_dNrgStateAbove)%dat               ,& ! intent(in): [dp(:)]  derivatives in  flux at layer interfaces w.r.t. states in the layer above
    dq_dNrgStateBelow            => deriv_data%var(iLookDERIV%dq_dNrgStateBelow)%dat               ,& ! intent(in): [dp(:)]  derivatives in  flux at layer interfaces w.r.t. states in the layer below
    dq_dNrgStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec)%dat        ,& ! intent(in): [dp(:)]  derivatives in  the flux in soil surface interface w.r.t. state variables in layers
      ! derivative in liquid water fluxes for the layer domains w.r.t temperature
    mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat                ,& ! intent(in): [dp(:)]  derivative in volumetric liquid water content w.r.t. temperature
    ! diagnostic variables
    scalarFracLiqVeg             => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)               ,& ! intent(in): [dp]     fraction of liquid water on vegetation (-)
    mLayerFracLiq                => diag_data%var(iLookDIAG%mLayerFracLiq)%dat                     ,& ! intent(in): [dp(:)]  fraction of liquid water in each snow, lake, or glce layer (-)
    scalarSoilControl            => diag_data%var(iLookDIAG%scalarSoilControl)%dat(1)              ,& ! intent(in): [dp]     soil control on infiltration for derivative
    scalarSoilControlBot         => diag_data%var(iLookDIAG%scalarSoilControlBot)%dat(1)           ,& ! intent(in): [dp]     soil control on bottom capillary fluxes for derivative
    mLayerVolFracIce             => prog_data%var(iLookPROG%mLayerVolFracIce)%dat                  ,& ! intent(in): [dp(:)]  volumetric fraction of ice in each layer start of step (-)
    ! canopy and layer depth
    canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)              ,& ! intent(in): [dp   ]  canopy depth (m)
    mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat                        & ! intent(in): [dp(:)]  depth of each layer in the sub-domain (m)
    ) ! making association with data in structures
    ! --------------------------------------------------------------
    ! initialize error control
    err=0; message='fluxJacAdd/'
    ! -----
    ! * energy and liquid fluxes over vegetation...
    ! ---------------------------------------------
    if(computeVegFlux)then ! (derivatives only defined when vegetation protrudes over the surface)

      ! * energy fluxes with the canopy water
      if(ixVegHyd/=integerMissing)then

        ! * cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
        if(ixCasNrg/=integerMissing) aJac(ixInd(full,ixVegHyd,ixCasNrg),ixCasNrg) = -dCanopyEvaporation_dTCanair*dt
        ! dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy is the derivative in throughfall and canopy drainage with canopy temperature
        if(ixVegNrg/=integerMissing) aJac(ixInd(full,ixVegHyd,ixVegNrg),ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy
        ! * liquid water fluxes for vegetation canopy (-), dt*scalarFracLiqVeg*scalarCanopyLiqDeriv is the derivative in throughfall and canopy drainage with canopy water
                                     aJac(ixInd(full,ixVegHyd,ixVegHyd),ixVegHyd) = -scalarFracLiqVeg*(dCanopyEvaporation_dCanWat - scalarCanopyLiqDeriv)*dt + dMat(ixVegHyd)
        if(ixTopNrg/=integerMissing) aJac(ixInd(full,ixVegHyd,ixTopNrg),ixTopNrg) = -dCanopyEvaporation_dTGround*dt

        ! * cross-derivative terms w.r.t. canopy water (kg-1 m2)
        if(nSnow>0)then
          if(ixTopHyd/=integerMissing) aJac(ixInd(full,ixTopHyd,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(1))*(-scalarFracLiqVeg*scalarCanopyLiqDeriv)/iden_water
        else
          if(ixTopHyd/=integerMissing) aJac(ixInd(full,ixTopHyd,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarFracLiqVeg*scalarCanopyLiqDeriv)/iden_water
        endif

        ! * cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
        if(ixTopNrg/=integerMissing) aJac(ixInd(full,ixTopNrg,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanWat)
      endif

      ! * cross-derivative terms w.r.t. canopy temperature (K-1)
      if(ixVegNrg/=integerMissing)then
        if(nSnow>0)then
          if(ixTopHyd/=integerMissing) aJac(ixInd(full,ixTopHyd,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(1))*(-scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water
        else
          if(ixTopHyd/=integerMissing) aJac(ixInd(full,ixTopHyd,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water
        endif
      endif

      ! * energy fluxes with the canopy air space (J m-3 K-1)
      if(ixCasNrg/=integerMissing)then
                                     aJac(ixInd(full,ixCasNrg,ixCasNrg),ixCasNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + dMat(ixCasNrg)
        if(ixVegNrg/=integerMissing) aJac(ixInd(full,ixCasNrg,ixVegNrg),ixVegNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
        if(ixTopNrg/=integerMissing) aJac(ixInd(full,ixCasNrg,ixTopNrg),ixTopNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)
      endif

      ! * energy fluxes with the vegetation canopy (J m-3 K-1)
      if(ixVegNrg/=integerMissing)then
        if(ixCasNrg/=integerMissing) aJac(ixInd(full,ixVegNrg,ixCasNrg),ixCasNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
                                     aJac(ixInd(full,ixVegNrg,ixVegNrg),ixVegNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + dMat(ixVegNrg)
        if(ixTopNrg/=integerMissing) aJac(ixInd(full,ixVegNrg,ixTopNrg),ixTopNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)
      endif

      ! * energy fluxes with the surface (J m-3 K-1)
      if(ixTopNrg/=integerMissing)then
        if(ixCasNrg/=integerMissing) aJac(ixInd(full,ixTopNrg,ixCasNrg),ixCasNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
        if(ixVegNrg/=integerMissing) aJac(ixInd(full,ixTopNrg,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)
      endif

    endif  ! if there is a need to compute energy fluxes within vegetation

    ! -----
    ! * energy fluxes for the layer domains...
    ! -------------------------------------------
    if(nSnLaSoGlNrg>0)then
      do iLayer=1,nLayers  ! loop through all layers in the layer domains

        ! check if the state is in the subset
        if(ixSnLaSoGlNrg(iLayer)==integerMissing) cycle
        ! - define index within the state subset and the full state vector
        nrgState = ixSnLaSoGlNrg(iLayer)        ! index within the state subset

        ! - diagonal elements
        aJac(ixInd(full,nrgState,nrgState),nrgState) = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(nrgState)

        ! - super-diagonal elements
        if(iLayer>1)then
          if(ixSnLaSoGlNrg(iLayer-1)/=integerMissing) aJac(ixInd(full,ixSnLaSoGlNrg(iLayer-1),nrgState),nrgState) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
        endif

        ! - sub-diagonal elements
        if(iLayer<nLayers)then
          if(ixSnLaSoGlNrg(iLayer+1)/=integerMissing) aJac(ixInd(full,ixSnLaSoGlNrg(iLayer+1),nrgState),nrgState) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
        endif

      end do ! (looping through energy states in the layer domains)
    endif ! (if the subset includes energy state variables in the layer domains)

    ! -----
    ! * liquid water fluxes for the snow, lake, glce domains...
    ! --------------------------------------------
    if(nSnowOnlyHyd+nLakeOnlyHyd+nGlceOnlyHyd>0)then
      do qLayer=1,nSnow+nLake+nGlce ! loop through layers in the snow, lake, glce domains

        if(qLayer<=nSnow+nLake)then
          jLayer = qLayer
          iLayer = qLayer
          endLayerWat = nSnow + nLake
          solid = .false.
          if(qLayer>nSnow .and. mLayerdTheta_dTk(jLayer) > tiny(1.0_rkind)) solid = .true. ! lake ice is solid
        else
          jLayer = qLayer + nSoil
          iLayer = qLayer - nSnow - nLake
          endLayerWat = nGlce - noThetaChange
          solid = .true.
        endif
        ! - check that the layer is desired
        if(ixSnLaSoGlHyd(jLayer)==integerMissing) cycle
        ! - define state indices for the current layer
        watState = ixSnLaSoGlHyd(jLayer)   ! hydrology state index within the state subset

        ! compute factor to convert liquid water derivative to total water derivative
        select case( ixHydType(jLayer) )
          case(iname_watLayer); convLiq2tot = mLayerFracLiq(jLayer)
          case default;         convLiq2tot = 1._rkind
        end select

        ! - diagonal elements, water does not move downwards in ice or upwards in snow
        ! NOTE: if unfrozen lake then dq_dHydStateBelow(0) should be some spillover from below, FIX with lakes  
        if(solid)then
          if(qLayer==nSnow+1 .or.(qLayer==nSnow+nLake+1 .and. nSoil==0))then
            aJac(ixInd(full,watState,watState),watState) = (dt/mLayerDepth(jLayer))*(-scalarSurfaceIceMeltDeriv*convLiq2tot) + dMat(watState)
          else
            aJac(ixInd(full,watState,watState),watState) = (dt/mLayerDepth(jLayer))*(-iLayerLiqFluxSnLaGlDeriv(jLayer-1)*convLiq2tot) + dMat(watState)
          endif
        else
          aJac(ixInd(full,watState,watState),watState) = (dt/mLayerDepth(jLayer))*iLayerLiqFluxSnLaGlDeriv(jLayer)*convLiq2tot + dMat(watState)
        endif

        ! - sub-diagonal elements for snow, sub-diagonal only (water does not move upwards in snow)
        if(qLayer<nSnow .and. mLayerVolFracIce(jLayer+1)<=maxVolIceContent_use)then
          if(ixSnLaSoGlHyd(jLayer+1)/=integerMissing) aJac(ixInd(full,ixSnLaSoGlHyd(jLayer+1),watState),watState) = -(dt/mLayerDepth(jLayer+1))*iLayerLiqFluxSnLaGlDeriv(jLayer)*convLiq2tot  ! dVol(below)/dLiq(above)
        endif
        ! - super-diagonal elements for ice, super-diagonal only (water does not move downwards in ice), but since water passes through immediately in ice, terms are 0

      end do ! (looping through liquid water states in the snow, lake, glce domains)
    endif ! (if the subset includes hydrology state variables in the snow, lake, glce domains)

    ! -----
    ! * cross derivatives in the snow, lake, glce domains...
    ! ----------------------------------------
    if((nSnowOnlyHyd>0 .and. nSnowOnlyNrg>0) .or. (nLakeOnlyHyd>0 .and. nLakeOnlyNrg>0) .or. (nGlceOnlyHyd>0 .and. nGlceOnlyNrg>0))then
      do qLayer=1,nSnow+nLake+nGlce-noThetaChange ! loop through layers in the snow, lake, glce domains

        if(qLayer<=nSnow+nLake)then
          jLayer = qLayer
          iLayer = qLayer
          endLayerWat = nSnow + nLake
          endLayerNrg = nSnow + nLake
          solid = .false.
          if(qLayer>nSnow .and. mLayerdTheta_dTk(jLayer) > tiny(1.0_rkind)) solid = .true. ! lake ice is solid
        else
          jLayer = qLayer + nSoil
          iLayer = qLayer - nSnow - nLake
          endLayerWat = nGlce - noThetaChange
          endLayerNrg = nGlce
          solid = .true.
        endif
        ! (define the energy state)
        nrgState = ixSnLaSoGlNrg(jLayer)       ! index within the full state vector
        ! - define state indices for the current layer
        watState = ixSnLaSoGlHyd(jLayer)   ! hydrology state index within the state subset

        if(nrgState/=integerMissing .and. watState/=integerMissing)then
          ! - include derivatives of water fluxes w.r.t energy fluxes for current layer, water does not move downwards in ice or upwards in snow
          ! NOTE: if unfrozen lake then dq_dNrgStateBelow(0) should be some spillover from below, FIX with lakes  
          if(solid)then
            if(qLayer==nSnow+1 .or. (qLayer==nSnow+nLake+1 .and. nSoil==0))then     
              aJac(ixInd(full,watState,nrgState),nrgState) = (dt/mLayerDepth(jLayer))*(-scalarSurfaceIceMeltDeriv*mLayerdTheta_dTk(jLayer))  ! (dVol/dT)
            else
              aJac(ixInd(full,watState,nrgState),nrgState) = (dt/mLayerDepth(jLayer))*(-iLayerLiqFluxSnLaGlDeriv(jLayer-1)*mLayerdTheta_dTk(jLayer))  ! (dVol/dT)
            endif
          else
            aJac(ixInd(full,watState,nrgState),nrgState) = (dt/mLayerDepth(jLayer))*iLayerLiqFluxSnLaGlDeriv(jLayer)*mLayerdTheta_dTk(jLayer)  ! (dVol/dT)
          endif
        endif ! (if the energy and water states for the current layer are within the state subset)

        if(watState/=integerMissing)then
          ! - include derivatives of heat capacity w.r.t water for layer above
          if(qLayer>1 .or. (qLayer==1 .and. nSnow==0 .and. nSoil>0))then ! have layer above
            if(ixSnLaSoGlNrg(jLayer-1)/=integerMissing) aJac(ixInd(full,ixSnLaSoGlNrg(jLayer-1),watState),watState) = (dt/mLayerDepth(jLayer-1))*( dNrgFlux_dWatBelow(jLayer-1) )
          endif

          ! - include derivatives of heat capacity w.r.t water for layer below
          if(iLayer<endLayerNrg .or. (iLayer==nSnow+nLake .and. nSoil+nGlce>0))then ! have layer below
            if(ixSnLaSoGlNrg(jLayer+1)/=integerMissing) aJac(ixInd(full,ixSnLaSoGlNrg(jLayer+1),watState),watState) = (dt/mLayerDepth(jLayer+1))*(-dNrgFlux_dWatAbove(jLayer) )
          endif
        endif ! (if the water state for the current layer is within the state subset)

        if(nrgState/=integerMissing)then
          ! - sub-diagonal elements for snow, sub-diagonal only (water does not move upwards in snow)
          if(qLayer<nSnow .and. mLayerVolFracIce(jLayer+1)<=maxVolIceContent_use)then
            if(ixSnLaSoGlHyd(jLayer+1)/=integerMissing) aJac(ixInd(full,ixSnLaSoGlHyd(jLayer+1),nrgState),nrgState) = -(dt/mLayerDepth(jLayer+1))*iLayerLiqFluxSnLaGlDeriv(jLayer)*mLayerdTheta_dTk(jLayer)  ! dVol(below)/dT(above)
          endif              
          ! - super-diagonal elements for ice, super-diagonal only (water does not move downwards in ice), but since water passes through immediately in ice, terms are 0
        endif ! (if the energy state for the current layer is within the state subset)

      end do ! (looping through snow, lake, glce layers)
    endif ! (if there are state variables for both water and energy in the snow, lake, glce domains)

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
        jLayer   = iLayer+nSnow+nLake            ! index of layer in the layer system

        ! - compute the diagonal elements
        ! all terms *excluding* baseflow
        aJac(ixInd(full,watState,watState),watState) = (dt/mLayerDepth(jLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(watState)

        ! - compute the super-diagonal elements
        if(iLayer>1)then
          if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixInd(full,ixSoilOnlyHyd(iLayer-1),watState),watState) = (dt/mLayerDepth(jLayer-1))*( dq_dHydStateBelow(iLayer-1))
        endif

        ! - compute the sub-diagonal elements
        if(iLayer<nSoil)then
          if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixInd(full,ixSoilOnlyHyd(iLayer+1),watState),watState) = (dt/mLayerDepth(jLayer+1))*(-dq_dHydStateAbove(iLayer))
        endif

        ! - include baseflow derivatives
        if(computeBaseflow .or. (nGlce>0 .and. nSoil>0))then ! include if glacier debris as it has lateral flow that behaves like baseflow lateral flow
          do pLayer=1,nSoil
            qState = ixSoilOnlyHyd(pLayer)  ! hydrology state index within the state subset
            if(qState/=integerMissing)then
              if((pLayer<=iLayer .and. watState - qstate <= kl) .or. (pLayer>iLayer .and. qstate - watState <= ku) .or. full) &
                aJac(ixInd(full,watState,qState),qState) = (dt/mLayerDepth(jLayer))*dBaseflow_dWat(iLayer,pLayer) + aJac(ixInd(full,watState,qState),qState)
            endif
          end do
        endif ! (if computed baseflow)

        ! - include derivatives for surface infiltration below surface
        if(ixSoilOnlyHyd(1)/=integerMissing .and. all(dq_dHydStateLayerSurfVec/=realMissing))then
           if(watState - ixSoilOnlyHyd(1) <= ku .or. full) &
               aJac(ixInd(full,ixSoilOnlyHyd(1),watState),watState) = -(dt/mLayerDepth(nSnow+nLake+1))*dq_dHydStateLayerSurfVec(iLayer) + aJac(ixInd(full,ixSoilOnlyHyd(1),watState),watState)
        endif
      end do ! (looping through hydrology states in the soil domain)

      ! - include derivatives for surface infiltration above surface if there is snow/lake (vegetation handled already)
      if(nSnow+nLake>0 .and. ixSoilOnlyHyd(1)/=integerMissing .and. all(dq_dHydStateLayerSurfVec/=realMissing))then 
        if(nSnow>0 .and. nLake==0)then ! have snow above first soil layer
          denseLimit = nSnow ! if passed through a too dense snowpack or lake, need to find top dense layer (bottom layer always included, dense or not)
          do pLayer=nSnow,1,-1
            if(mLayerVolFracIce(pLayer)<=maxVolIceContent_use) exit
            denseLimit = pLayer
          end do
          endLayerWat = nSnow
        else ! have lake above first soil layer
          denseLimit = nSnow+nLake
          endLayerWat = nSnow+nLake
          if(mLayerVolFracIce(denseLimit)>maxVolIceContent_use) endLayerWat = denseLimit-1_i4b ! lake frozen solid, no infiltration into soil so don't include
        endif
        do pLayer=denseLimit,endLayerWat
          if(ixSnLaSoGlHyd(pLayer)/=integerMissing)then
            ! compute factor to convert liquid water derivative to total water derivative
            select case( ixHydType(pLayer) )
              case(iname_watLayer); convLiq2tot = mLayerFracLiq(pLayer)
              case default;         convLiq2tot = 1._rkind
            end select
            if(ixSoilOnlyHyd(1) - ixSnLaSoGlHyd(pLayer) <= kl .or. full) &
                aJac(ixInd(full,ixSoilOnlyHyd(1),ixSnLaSoGlHyd(pLayer)),ixSnLaSoGlHyd(pLayer)) = -(dt/mLayerDepth(nSnow+nLake+1))*scalarSoilControl*iLayerLiqFluxSnLaGlDeriv(pLayer)*convLiq2tot + aJac(ixInd(full,ixSoilOnlyHyd(1),ixSnLaSoGlHyd(pLayer)),ixSnLaSoGlHyd(pLayer))
          endif
        end do ! (looping through snow/lake layers above soil until non-dense layer)
      endif ! (if snow or lake present above soil)

      ! - include derivatives for flux into bottom soil layer if there is glacier ice
      if(nSoil>0 .and. ixSoilOnlyHyd(nSoil)/=integerMissing .and. nGlce>0)then
        do pLayer=nSnow+nLake+nSoil+1,nLayers-noThetaChange
          if(ixSnLaSoGlHyd(pLayer)/=integerMissing)then
            ! compute factor to convert liquid water derivative to total water derivative
            select case( ixHydType(pLayer) )
              case(iname_watLayer); convLiq2tot = mLayerFracLiq(pLayer)
              case default;         convLiq2tot = 1._rkind
            end select
            if(ixSnLaSoGlHyd(pLayer) - ixSoilOnlyHyd(nSoil) <= ku .or. full) &
                aJac(ixInd(full,ixSoilOnlyHyd(nSoil),ixSnLaSoGlHyd(pLayer)),ixSnLaSoGlHyd(pLayer)) = (dt/mLayerDepth(nSnow+nLake+nSoil))*(-scalarSoilControlBot*iLayerLiqFluxSnLaGlDeriv(pLayer)*convLiq2tot) + aJac(ixInd(full,ixSoilOnlyHyd(nSoil),ixSnLaSoGlHyd(pLayer)),ixSnLaSoGlHyd(pLayer))
          endif
        end do ! (looping through glacier ice layers below soil)
      endif ! (if glacier ice present below soil)

    endif ! (if the subset includes hydrology state variables in the soil domain)

    ! -----
    ! * liquid water fluxes for the aquifer...
    ! ----------------------------------------
    if(ixAqWat/=integerMissing)then
      aJac(ixInd(full,ixAqWat,ixAqWat),ixAqWat) = -dBaseflow_dAquifer*dt + dMat(ixAqWat)
      if(nSoil>0)then
        if(ixSoilOnlyNrg(nSoil)/=integerMissing) aJac(ixInd(full,ixAqWat,ixSoilOnlyNrg(nSoil)),ixSoilOnlyNrg(nSoil)) = -dq_dNrgStateAbove(nSoil)*dt ! dAquiferRecharge_dTk  = d_iLayerLiqFluxSoil(nSoil)_dTk
        if(ixSoilOnlyHyd(nSoil)/=integerMissing) aJac(ixInd(full,ixAqWat,ixSoilOnlyHyd(nSoil)),ixSoilOnlyHyd(nSoil)) = -dq_dHydStateAbove(nSoil)*dt ! dAquiferRecharge_dWat = d_iLayerLiqFluxSoil(nSoil)_dWat
      endif
      ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
      if(computeVegFlux)then
        if(ixCasNrg/=integerMissing)then
          if(ixAqWat-ixCasNrg <= kl .or. full) aJac(ixInd(full,ixAqWat,ixCasNrg),ixCasNrg) = -dAquiferTrans_dTCanair*dt ! dVol/dT (K-1)
        endif
        if(ixVegNrg/=integerMissing)then
          if(ixAqWat-ixVegNrg <= kl .or. full) aJac(ixInd(full,ixAqWat,ixVegNrg),ixVegNrg) = -dAquiferTrans_dTCanopy*dt ! dVol/dT (K-1)
        endif
        if(ixVegHyd/=integerMissing)then
          if(ixAqWat-ixVegHyd <= kl .or. full) aJac(ixInd(full,ixAqWat,ixVegHyd),ixVegHyd) = -dAquiferTrans_dCanWat*dt  ! dVol/dLiq (kg m-2)-1
        endif
        if(ixTopNrg/=integerMissing)then
          if(ixAqWat-ixTopNrg <= kl .or. full) aJac(ixInd(full,ixAqWat,ixTopNrg),ixTopNrg) = -dAquiferTrans_dTGround*dt ! dVol/dT (K-1)
        endif
      endif
    endif ! (if aquifer water state is in the subset)å

    ! -----
    ! * cross derivatives in the soil domain...
    ! ----------------------------------------
    if(nSoilOnlyHyd>0 .and. (nSoilOnlyNrg>0 .or. computeVegFlux))then
      do iLayer=1,nSoilOnlyNrg

        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow+nLake            ! index of layer in the layer system
        ! - define the energy state variable
        nrgState = ixSoilOnlyNrg(iLayer)         ! index within the full state vector
        ! - define index of hydrology state variable within the state subset
        watState = ixSoilOnlyHyd(iLayer)

        if(watState/=integerMissing .and. nrgState/=integerMissing)then
          ! - include derivatives in liquid water fluxes w.r.t. temperature for current layer
          aJac(ixInd(full,watState,nrgState),nrgState) = (dt/mLayerDepth(jLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance

         ! - include derivatives w.r.t. ground evaporation
          if(nSnow==0 .and. iLayer==1)then 
            aJac(ixInd(full,ixTopHyd,ixTopNrg),ixTopNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac(ixInd(full,ixTopHyd,ixTopNrg),ixTopNrg) ! dVol/dT (K-1)
          endif
        endif   !(if both the energy and water states for the current layer are within the state subset)

        if(watState/=integerMissing)then
          ! - include derivatives w.r.t. ground evaporation
          if(nSnow+nLake==0 .and. iLayer==1)then 
            if(computeVegFlux)then ! surface soil layer, assume here that kl>=4
              if(ixCasNrg/=integerMissing) aJac(ixInd(full,ixTopHyd,ixCasNrg),ixCasNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanair/iden_water) ! dVol/dT (K-1)
              if(ixVegNrg/=integerMissing) aJac(ixInd(full,ixTopHyd,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanopy/iden_water) + aJac(ixInd(full,ixTopHyd,ixVegNrg),ixVegNrg) ! dVol/dT (K-1)
              if(ixVegHyd/=integerMissing) aJac(ixInd(full,ixTopHyd,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dCanWat/iden_water)  + aJac(ixInd(full,ixTopHyd,ixVegHyd),ixVegHyd) ! dVol/dLiq (kg m-2)-1
            endif
          endif

           ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
          if(computeVegFlux)then
            if(ixCasNrg/=integerMissing)then
              if(watState-ixCasNrg <= kl .or. full) aJac(ixInd(full,watState,ixCasNrg),ixCasNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTCanair(iLayer)) + aJac(ixInd(full,watState,ixCasNrg),ixCasNrg) ! dVol/dT (K-1)
            endif
            if(ixVegNrg/=integerMissing)then
              if(watState-ixVegNrg <= kl .or. full) aJac(ixInd(full,watState,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTCanopy(iLayer)) + aJac(ixInd(full,watState,ixVegNrg),ixVegNrg) ! dVol/dT (K-1)
            endif
            if(ixVegHyd/=integerMissing)then
              if(watState-ixVegHyd <= kl .or. full) aJac(ixInd(full,watState,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dCanWat(iLayer))  + aJac(ixInd(full,watState,ixVegHyd),ixVegHyd) ! dVol/dLiq (kg m-2)-1
            endif
            if(ixTopNrg/=integerMissing)then
              if(watState-ixTopNrg <= kl .or. full) aJac(ixInd(full,watState,ixTopNrg),ixTopNrg) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTGround(iLayer)) + aJac(ixInd(full,watState,ixTopNrg),ixTopNrg) ! dVol/dT (K-1)
            endif
          endif

           ! - include derivatives of heat capacity w.r.t water fluxes for layer above
          if(iLayer>1 .or. (iLayer==1 .and. nSnow+nLake>0))then ! have layer above
            if(ixSnLaSoGlNrg(jLayer-1)/=integerMissing) aJac(ixInd(full,ixSnLaSoGlNrg(jLayer-1),watState),watState) = (dt/mLayerDepth(jLayer-1))*( dNrgFlux_dWatBelow(jLayer-1) )
          endif

          ! - include derivatives of heat capacity w.r.t water fluxes for layer below
          if(iLayer<nSoil .or. (iLayer==nSoil .and. nGlce>0))then ! have layer below
            if(ixSnLaSoGlNrg(jLayer+1)/=integerMissing) aJac(ixInd(full,ixSnLaSoGlNrg(jLayer+1),watState),watState) = (dt/mLayerDepth(jLayer+1))*(-dNrgFlux_dWatAbove(jLayer) )
          endif

          ! - include baseflow derivatives
          if(computeBaseflow .or. (nGlce>0 .and. nSoil>0))then ! include if glacier debris as it has lateral flow that behaves like baseflow lateral flow
            do pLayer=1,nSoil
              qState = ixSoilOnlyNrg(pLayer)  ! hydrology state index within the state subset
              if(qState/=integerMissing)then
                if((pLayer<=iLayer .and. watState - qstate <= kl) .or. (pLayer>iLayer .and. qstate - watState <= ku) .or. full) &
                    aJac(ixInd(full,watState,qState),qState) = (dt/mLayerDepth(jLayer))*dBaseflow_dTk(iLayer,pLayer) + aJac(ixInd(full,watState,qState),qState)
              endif
            end do
          endif ! (if computed baseflow)
        endif ! (if the water state for the current layer is within the state subset)

        if(nrgState/=integerMissing)then
          ! - compute super-diagonal elements
          if(iLayer>1)then
            if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixInd(full,ixSoilOnlyHyd(iLayer-1),nrgState),nrgState) = (dt/mLayerDepth(jLayer-1))*( dq_dNrgStateBelow(iLayer-1)) + aJac(ixInd(full,ixSoilOnlyHyd(iLayer-1),nrgState),nrgState) ! K-1
          endif

          ! compute sub-diagonal elements
          if(iLayer<nSoil)then
            if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixInd(full,ixSoilOnlyHyd(iLayer+1),nrgState),nrgState) = (dt/mLayerDepth(jLayer+1))*(-dq_dNrgStateAbove(iLayer)) + aJac(ixInd(full,ixSoilOnlyHyd(iLayer+1),nrgState),nrgState) ! K-1
          endif

          ! - include derivatives for surface infiltration below surface
          if(ixSoilOnlyHyd(1)/=integerMissing .and. all(dq_dNrgStateLayerSurfVec/=realMissing))then
            if(nrgState - ixSoilOnlyHyd(1) <= ku .or. full) &
                aJac(ixInd(full,ixSoilOnlyHyd(1),nrgState),nrgState) = -(dt/mLayerDepth(nSnow+nLake+1))*dq_dNrgStateLayerSurfVec(iLayer) + aJac(ixInd(full,ixSoilOnlyHyd(1),nrgState),nrgState)
          endif
        endif ! (if the energy state for the current layer is within the state subset)
      end do ! (looping through energy states in the soil domain)

      ! - include derivatives for surface infiltration above surface if there is snow/lake (vegetation handled already)
      if(nSnow+nLake>0 .and. ixSoilOnlyHyd(1)/=integerMissing .and. all(dq_dNrgStateLayerSurfVec/=realMissing))then 
        if(nSnow>0 .and. nLake==0)then ! have snow above first soil layer
          denseLimit = nSnow ! if passed through a too dense snowpack or lake, need to find top dense layer (bottom layer always included, dense or not)
          do pLayer=nSnow,1,-1
            if(mLayerVolFracIce(pLayer)<=maxVolIceContent_use) exit
            denseLimit = pLayer
          end do
          endLayerNrg = nSnow
        else ! have lake above first soil layer
          denseLimit = nSnow+nLake ! bottom lake layer
          endLayerNrg = nSnow+nLake
        endif
        do pLayer=denseLimit,endLayerNrg
          if(ixSnLaSoGlNrg(pLayer)/=integerMissing)then
            if(ixSoilOnlyHyd(1) - ixSnLaSoGlNrg(pLayer) <= kl .or. full) &
                aJac(ixInd(full,ixSoilOnlyHyd(1),ixSnLaSoGlNrg(pLayer)),ixSnLaSoGlNrg(pLayer)) = -(dt/mLayerDepth(nSnow+nLake+1))*scalarSoilControl*iLayerLiqFluxSnLaGlDeriv(pLayer)*mLayerdTheta_dTk(pLayer) + aJac(ixInd(full,ixSoilOnlyHyd(1),ixSnLaSoGlNrg(pLayer)),ixSnLaSoGlNrg(pLayer))
          endif
        end do ! (looping through snow/lake layers above soil until non-dense layer)
      endif ! (if snow or lake present above soil)

      ! - include derivatives for flux into bottom soil layer if there is glacier ice
      if(nSoil>0 .and. ixSoilOnlyHyd(nSoil)/=integerMissing .and. nGlce>0)then
        do pLayer=nSnow+nLake+nSoil+1,nLayers-noThetaChange
          if(ixSnLaSoGlNrg(pLayer)/=integerMissing)then
            if(ixSnLaSoGlNrg(pLayer) - ixSoilOnlyHyd(nSoil) <= ku .or. full) &
                aJac(ixInd(full,ixSoilOnlyHyd(nSoil),ixSnLaSoGlNrg(pLayer)),ixSnLaSoGlNrg(pLayer)) = (dt/mLayerDepth(nSnow+nLake+nSoil))*(-scalarSoilControlBot*iLayerLiqFluxSnLaGlDeriv(pLayer)*mLayerdTheta_dTk(pLayer)) + aJac(ixInd(full,ixSoilOnlyHyd(nSoil),ixSnLaSoGlNrg(pLayer)),ixSnLaSoGlNrg(pLayer))
          endif
        end do ! (looping through glacier ice layers below soil)
      endif ! (if glacier ice present below soil)

    endif ! (if there are state variables for both water and energy in the soil domain)

  end associate ! end association to variables in the data structures

end subroutine fluxJacAdd

! **********************************************************************************************************
! public function: get the index in the band-diagonal matrix or full matrix
! **********************************************************************************************************
function ixInd(full,jState,iState)
  implicit none
  logical(lgt),intent(in)  :: full   ! true if using full matrix, false if using band-diagonal matrix
  integer(i4b),intent(in)  :: jState ! off-diagonal state
  integer(i4b),intent(in)  :: iState ! diagonal state
  integer(i4b)             :: ixInd  ! index in the band-diagonal matrix or full matrix

  if(full) then
    ixInd = jState
  else
    ixInd = ixDiag + jState - iState
  endif
end function ixInd

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
  use fsundials_core_mod
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

  ! class objects for subroutine arguments
  type(in_type_computJacob)     :: in_computJacob  ! intent(in)  computJacob arguments
  type(out_type_computJacob)    :: out_computJacob ! intent(out) computJacob arguments
! ----------------------------------------------------------------

  ! get equations data from user-defined data
  call c_f_pointer(user_data,eqns_data)

  ! get data arrays from SUNDIALS vectors
  if(eqns_data%ixMatrix==ixBandMatrix) Jac(1:nBands, 1:eqns_data%nState) => FSUNBandMatrix_Data(sunmat_J)
  if(eqns_data%ixMatrix==ixFullMatrix) Jac(1:eqns_data%nState, 1:eqns_data%nState) => FSUNDenseMatrix_Data(sunmat_J)

  ! compute the analytical Jacobian matrix
  ! NOTE: The derivatives were computed in the previous call to computFlux
  !       This occurred either at the call to eval8summa at the start of systemSolv
  !        or in the call to eval8summa in the previous iteration
  call initialize_computJacob ! pack in_computJacob object
  call computJacob(&
                ! input: model control
                in_computJacob, &
                ! input: data structures
                eqns_data%indx_data,               & ! intent(in):    index data
                eqns_data%prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                eqns_data%diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                eqns_data%deriv_data,              & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                eqns_data%dBaseflow_dWat,          & ! intent(in):    derivative in baseflow w.r.t. soil water characteristic
                eqns_data%dBaseflow_dTk,           & ! intent(in):    derivative in baseflow w.r.t. temperature (m s-1 K-1)
                ! input-output: Jacobian and its diagonal
                eqns_data%dMat,                    & ! intent(inout): diagonal of the Jacobian matrix
                Jac,                               & ! intent(out):   Jacobian matrix
                ! output: error control
                out_computJacob)                     ! intent(out):   error code and error message 
  call finalize_computJacob ! unpack out_computJacob object
  if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
  if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif                                  

  ! return success
  ierr = 0
  return

 contains

  subroutine initialize_computJacob
   ! *** Transfer data to in_computJacob class object from local variables ***
   call in_computJacob % initialize(eqns_data%dt_cur,eqns_data%nSnow,eqns_data%nLake,eqns_data%nSoil,eqns_data%nGlce,eqns_data%nLayers,eqns_data%computeVegFlux,&
      (eqns_data%model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel),eqns_data%ixMatrix)
  end subroutine initialize_computJacob

  subroutine finalize_computJacob
   ! *** Transfer data from out_computJacob class object to local variables ***
   call out_computJacob % finalize(eqns_data % err,eqns_data % message)
  end subroutine finalize_computJacob
            
end function computJacob4kinsol
#endif

end module computJacob_module
