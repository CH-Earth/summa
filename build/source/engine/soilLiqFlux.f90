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

module soilLiqFlux_module
! -----------------------------------------------------------------------------------------------------------

! data types
USE nr_type
USE data_types,only:&
                   var_ilength,           & ! x%var(:)%dat   (i4b)
                   var_dlength,           & ! x%var(:)%dat   (rkind)
                   in_type_soilLiqFlux,   & ! derived type for intent(in) arguments
                   io_type_soilLiqFlux,   & ! derived type for intent(inout) arguments
                   out_type_soilLiqFlux,  & ! derived type for intent(out) arguments
                   in_type_diagv_node,    & ! derived type for intent(in) arguments
                   out_type_diagv_node,   & ! derived type for intent(out) arguments
                   in_type_surfaceFlux,   & ! derived type for intent(in) arguments
                   io_type_surfaceFlux,   & ! derived type for intent(inout) arguments
                   out_type_surfaceFlux,  & ! derived type for intent(out) arguments
                   in_type_iLayerFlux,    & ! derived type for intent(in) arguments
                   out_type_iLayerFlux,   & ! derived type for intent(out) arguments
                   in_type_qDrainFlux,    & ! derived type for intent(in) arguments
                   io_type_qDrainFlux,    & ! derived type for intent(inout) arguments
                   out_type_qDrainFlux      ! derived type for intent(out) arguments

! missing values
USE globalData,only:integerMissing         ! missing integer
USE globalData,only:realMissing            ! missing real number

! constants
USE multiconst,only:iden_water             ! intrinsic density of water    (kg m-3)
USE globalData,only:veryBig                ! a very big number
USE globalData,only:verySmall              ! a small number
USE globalData,only:verySmaller            ! a smaller number than verySmall

! named variables
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements
USE var_lookup,only:iLookFLUX              ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookINDEX             ! named variables for structure elements

! model decisions
USE globalData,only:model_decisions        ! model decision structure
USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure

! provide access to look-up values for model decisions
USE mDecisions_module,only:   &
  ! look-up values for the choice of boundary conditions for hydrology
  prescribedHead,             & ! prescribed head
  funcBottomHead,             & ! function of matric head in the lower-most layer
  freeDrainage,               & ! free drainage
  liquidFlux,                 & ! liquid water flux
  zeroFlux,                   & ! zero flux
  ! look-up values for the choice of saturation excesssurface runoff parameterization
  zero_SE,                    & ! zero saturation excess surface runoff parameterization 
  homegrown_SE,               & ! homegrown saturation excess surface runoff parameterization 
  FUSEPRMS,                   & ! FUSE PRMS     surface runoff parameterization 
  FUSEAVIC,                   & ! FUSE ARNO/VIC surface runoff parameterization
  FUSETOPM,                   & ! FUSE TOPMODEL surface runoff parameterization
  ! look-up values for the maximum infiltration rate parameterization
  GreenAmpt,                  & ! Green-Ampt parameterization
  topmodel_GA,                & ! Green-Ampt parameterization with conductivity profile from TOPMODEL-ish parameterization
  noInfiltrationExcess,       & ! no infiltration excess runoff
  ! look-up values for the choice of hydraulic conductivity profile
  constant,                   & ! constant hydraulic conductivity with depth
  powerLaw_profile,           & ! power-law profile
  expLaw_profile,             & ! exponential profile
  ! look-up values for the choice of groundwater parameterization
  qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                  & ! a big bucket (lumped aquifer model)
  noExplicit                    ! no explicit groundwater parameterization

! -----------------------------------------------------------------------------------------------------------
implicit none
private
public::soilLiqFlux

! flag to denote if updating infiltration during iterations for testing purposes
logical(lgt),parameter :: updateInfil=.true.
contains

! ***************************************************************************************************************
! public subroutine soilLiqFlux: compute liquid water fluxes and their derivatives
! ***************************************************************************************************************
subroutine soilLiqFlux(&
                      ! input: model control, trial state variables, derivatives, and fluxes
                      in_soilLiqFlux,                & ! intent(in): model control, trial state variables, derivatives, and fluxes
                      ! input-output: data structures
                      mpar_data,                    & ! intent(in):    model parameters
                      indx_data,                    & ! intent(in):    model indices
                      prog_data,                    & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,                    & ! intent(inout): model fluxes for a local HRU
                      ! input-output: diagnostic variables, fluxes, and derivatives
                      io_soilLiqFlux,                & ! intent(inout): diagnostic variables, fluxes, and derivatives
                      ! output: error control
                      out_soilLiqFlux)                 ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control, trial state variables, derivatives, and fluxes
  type(in_type_soilLiqFlux),intent(in)    :: in_soilLiqFlux              ! model control, trial state variables, derivatives, and fluxes
  ! input-output: data structures
  type(var_dlength),intent(in)            :: mpar_data                   ! model parameters
  type(var_ilength),intent(in)            :: indx_data                   ! state vector geometry
  type(var_dlength),intent(in)            :: prog_data                   ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)         :: diag_data                   ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)         :: flux_data                   ! model fluxes for a local HRU
  ! input-output: diagnostic variables, fluxes, and derivatives
  type(io_type_soilLiqFlux),intent(inout) :: io_soilLiqFlux              ! diagnostic variables, fluxes, and derivatives
  ! output: error control
  type(out_type_soilLiqFlux),intent(out)  :: out_soilLiqFlux             ! error code and error message
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables: general
  character(LEN=256)                               :: cmessage            ! error message of downwind routine
  integer(i4b)                                     :: nSoil               ! number of soil layers
  integer(i4b)                                     :: nGlce               ! number of glacier ice layers
  integer(i4b)                                     :: ibeg,iend           ! start and end indices of the soil layers in concatanated snow-lake-soil-glce vector
  integer(i4b)                                     :: iLayer,iSoil        ! index of soil layer
  integer(i4b)                                     :: ixLayerDesired(1)   ! layer desired (scalar solution)
  integer(i4b)                                     :: ixTop               ! top layer in subroutine call
  integer(i4b)                                     :: ixBot               ! bottom layer in subroutine call
  ! transpiration sink term
  real(rkind),dimension(in_soilLiqFlux % nSoil)    :: mLayerTranspireFrac ! fraction of transpiration allocated to each soil layer (-)
  ! diagnostic variables
  real(rkind),dimension(in_soilLiqFlux % nSoil)    :: iceImpedeFac        ! ice impedence factor at layer mid-points (-)
  real(rkind),dimension(in_soilLiqFlux % nSoil)    :: dHydCond_dTemp      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  real(rkind),dimension(0:in_soilLiqFlux % nSoil)  :: iLayerHydCond       ! hydraulic conductivity at layer interface (m s-1)
  ! compute surface flux
  integer(i4b)                                     :: nRoots              ! number of soil layers with roots or layers that take infiltration
  integer(i4b)                                     :: ixIce               ! index of the lowest soil layer that contains ice
  real(rkind),dimension(0:in_soilLiqFlux % nSoil)  :: iLayerHeight        ! height of the layer interfaces (m)
  ! error control
  logical(lgt)                                     :: return_flag         ! flag for return statements
  ! -------------------------------------------------------------------------------------------------------------------------------------------------

  ! ** Initialize indices, error control, and get layer information **
  call initialize_soilLiqFlux; if (return_flag) return

  ! ** Compute transpiration, diagnostic variables, infiltration, and interface fluxes **
  call update_soilLiqFlux;     if (return_flag) return

  ! ** Final error control **
  call finalize_soilLiqFlux;   if (return_flag) return

contains

 subroutine initialize_soilLiqFlux
  ! **** Initial operations for soilLiqFlux module subroutine ****

  ! ** assign variables used in main associate block **
  nSoil = in_soilLiqFlux % nSoil ! get number of soil layers from input arguments
  nGlce = indx_data%var(iLookINDEX%nGlce)%dat(1) ! get number of glacier ice layers from index data structure

  ! get indices for the data structures
  ibeg = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nLake)%dat(1) + 1
  iend = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nLake)%dat(1) + nSoil

  ! get a copy of iLayerHeight (for soil layers only)
  ! NOTE: performance hit, though cannot define the shape (0:) with the associate construct
  iLayerHeight(0:nSoil) = prog_data%var(iLookPROG%iLayerHeight)%dat(ibeg-1:iend)  ! height of the layer interfaces (m)

  ! ** initialize error control **
  return_flag=.false.
  associate(&
    err                   => out_soilLiqFlux % err,                  & ! intent(out): error code
    message               => out_soilLiqFlux % cmessage              & ! intent(out): error message
  &)
   err=0; message='soilLiqFlux/' ! initialize error control
  end associate

  ! ** get the indices for the soil layers **
  associate(&
   scalarSolution => in_soilLiqFlux % scalarSolution,            & ! intent(in): flag to denote if implementing the scalar solution
   ixMatricHead   => indx_data%var(iLookINDEX%ixMatricHead)%dat, & ! intent(in): indices of soil layers where matric head is the state variable
   ixSoilOnlyHyd  => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat & ! intent(in): index in the state subset for hydrology state variables in the soil domain
  &)
   if (scalarSolution) then
     ixLayerDesired = pack(ixMatricHead, ixSoilOnlyHyd/=integerMissing)
     ixTop = ixLayerDesired(1)
     ixBot = ixLayerDesired(1)
   else
     ixTop = 1
     ixBot = nSoil
   end if
  end associate

  ! ** identify the number of layers that contain roots or take infiltration **
  associate(&
   rootingDepth => mpar_data%var(iLookPARAM%rootingDepth)%dat(1), & ! intent(in): rooting depth (m)
   err          => out_soilLiqFlux % err,                         & ! intent(out): error code
   message      => out_soilLiqFlux % cmessage                     & ! intent(out): error message
  &)
   nRoots = count(iLayerHeight(0:nSoil-1) < rootingDepth-verySmall)
   if(nGlce>0) nRoots = nSoil ! if glacier ice exists, then all soil layers are considered to have roots
   if(nRoots==0)then; message=trim(message)//'no layers with roots/infiltration'; err=20; return_flag=.true.; return; end if
  end associate

  ! ** identify lowest soil layer with ice **
  ! NOTE: cannot use count because there may be an unfrozen wedge
  associate(&
    mLayerVolFracIceTrial => in_soilLiqFlux % mLayerVolFracIceTrial & ! intent(in): volumetric fraction of ice at the current iteration (-)
  &)
   ixIce = 0  ! initialize the index of the ice layer (0 means no ice in the soil profile)
   do iLayer=1,nSoil ! (loop through soil layers)
     if (mLayerVolFracIceTrial(iLayer) > verySmaller) ixIce = iLayer
   end do
  end associate
 end subroutine initialize_soilLiqFlux

 subroutine update_soilLiqFlux
  ! **** Main computations for soilLiqFlux module subroutine ****

  if ( .not. (in_soilLiqFlux % scalarSolution .and. ixTop>1) ) then ! check the need to compute transpiration
   call compute_transpiration_sink; if (return_flag) return
  end if

  call compute_diagnostic_variables; if (return_flag) return

  call compute_surface_infiltration; if (return_flag) return

  call compute_interface_fluxes_derivatives; if (return_flag) return

  if ( .not. (in_soilLiqFlux % scalarSolution .and. ixTop<nSoil) ) then ! define the need to compute drainage
   call compute_drainage_flux; if (return_flag) return
  end if
 end subroutine update_soilLiqFlux

 subroutine finalize_soilLiqFlux
  ! **** Final operations for soilLiqFlux module subroutine ****

  ! final error control check for robustness
  associate(&
   err          => out_soilLiqFlux % err,                         & ! intent(out): error code
   message      => out_soilLiqFlux % cmessage                     & ! intent(out): error message
  &)
   if(err/=0)then; message=trim(message)//trim("finalize_soilLiqFlux: final error check failed"); return_flag=.true.; return; end if
  end associate
 end subroutine finalize_soilLiqFlux

 subroutine compute_transpiration_sink
  ! **** Compute the transpiration sink term ****

  call update_transpiration_loss_fraction
  call finalize_transpiration_loss_fraction; if (return_flag) return

  call update_transpiration_loss
 end subroutine compute_transpiration_sink

 subroutine update_transpiration_loss_fraction
  ! **** Update the fraction of transpiration loss from each soil layer *****
  associate(&
   scalarTranspireLim => diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1), & ! intent(in): weighted average of the transpiration limiting factor (-)
   mLayerRootDensity  => diag_data%var(iLookDIAG%mLayerRootDensity)%dat,     & ! intent(in): root density in each layer (-)
   mLayerTranspireLim => diag_data%var(iLookDIAG%mLayerTranspireLim)%dat     & ! intent(in): transpiration limiting factor in each layer (-)
  &)
   ! transpiration may be non-zero even if the soil moisture limiting factor is zero
   if (scalarTranspireLim > tiny(scalarTranspireLim)) then
    mLayerTranspireFrac(:) = mLayerRootDensity(:)*mLayerTranspireLim(:)/scalarTranspireLim
   else ! possibility of non-zero conductance and therefore transpiration in this case
    mLayerTranspireFrac(:) = mLayerRootDensity(:) / sum(mLayerRootDensity)
   end if
  end associate
 end subroutine update_transpiration_loss_fraction

 subroutine finalize_transpiration_loss_fraction
  ! **** Finalize operations for the fraction of transpiration loss from each soil layer *****
  associate(&
   err          => out_soilLiqFlux % err,     & ! intent(out): error code
   message      => out_soilLiqFlux % cmessage & ! intent(out): error message
  &)
   ! check fractions sum to one
   if (abs(sum(mLayerTranspireFrac) - 1._rkind) > verySmaller) then
     message=trim(message)//'fraction transpiration in soil layers does not sum to one'; err=20; return_flag=.true.; return
   end if
  end associate
 end subroutine finalize_transpiration_loss_fraction

 subroutine update_transpiration_loss
  ! **** Update transpiration loss from each soil layer (kg m-2 s-1 --> m s-1)*****
  associate(&
   scalarCanopyTranspiration => in_soilLiqFlux % scalarCanopyTranspiration, & ! canopy transpiration (kg m-2 s-1)
   mLayerTranspire           => io_soilLiqFlux % mLayerTranspire,   & ! transpiration loss from each soil layer (m s-1)
   ! intent(inout): derivatives in the soil layer transpiration flux ...
   mLayerdTrans_dCanWat  => io_soilLiqFlux % mLayerdTrans_dCanWat,  & ! ... w.r.t. canopy total water
   mLayerdTrans_dTCanair => io_soilLiqFlux % mLayerdTrans_dTCanair, & ! ... w.r.t. canopy air temperature
   mLayerdTrans_dTCanopy => io_soilLiqFlux % mLayerdTrans_dTCanopy, & ! ... w.r.t. canopy temperature
   mLayerdTrans_dTGround => io_soilLiqFlux % mLayerdTrans_dTGround, & ! ... w.r.t. ground temperature
   ! intent(in): derivative in canopy transpiration ...
   dCanopyTrans_dCanWat  => in_soilLiqFlux % dCanopyTrans_dCanWat,  & ! ... w.r.t. canopy total water content (s-1)
   dCanopyTrans_dTCanair => in_soilLiqFlux % dCanopyTrans_dTCanair, & ! ... w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTCanopy => in_soilLiqFlux % dCanopyTrans_dTCanopy, & ! ... w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTGround => in_soilLiqFlux % dCanopyTrans_dTGround, & ! ... w.r.t. ground temperature (kg m-2 s-1 K-1)
   ! intent(in): index of the upper boundary conditions for soil hydrology
   ixBcUpperSoilHydrology => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision &
  &)
   if (ixBcUpperSoilHydrology==prescribedHead) then ! special case of prescribed head -- no transpiration
    mLayerTranspire(:)      = 0._rkind
    ! derivatives in transpiration w.r.t. canopy state variables
    mLayerdTrans_dCanWat(:) = 0._rkind
    mLayerdTrans_dTCanair(:)= 0._rkind
    mLayerdTrans_dTCanopy(:)= 0._rkind
    mLayerdTrans_dTGround(:)= 0._rkind
   else
    mLayerTranspire(:) = mLayerTranspireFrac(:)*scalarCanopyTranspiration/iden_water
    ! * derivatives in transpiration w.r.t. canopy state variables *
    mLayerdTrans_dCanWat(:)  = mLayerTranspireFrac(:)*dCanopyTrans_dCanWat /iden_water
    mLayerdTrans_dTCanair(:) = mLayerTranspireFrac(:)*dCanopyTrans_dTCanair/iden_water
    mLayerdTrans_dTCanopy(:) = mLayerTranspireFrac(:)*dCanopyTrans_dTCanopy/iden_water
    mLayerdTrans_dTGround(:) = mLayerTranspireFrac(:)*dCanopyTrans_dTGround/iden_water
   end if
  end associate
 end subroutine update_transpiration_loss

 subroutine compute_diagnostic_variables
  ! **** compute diagnostic variables at the nodes throughout the soil profile ****
  type(in_type_diagv_node)  :: in_diagv_node  ! input data object for diagv_node
  type(out_type_diagv_node) :: out_diagv_node ! output data object for diagv_node

  do iSoil=ixTop,min(ixBot+1,nSoil) ! loop through soil layers

   call initialize_compute_diagnostic_variables(in_diagv_node)

   call update_compute_diagnostic_variables(in_diagv_node,out_diagv_node)

   call finalize_compute_diagnostic_variables(out_diagv_node); if (return_flag) return

  end do
 end subroutine compute_diagnostic_variables

 subroutine initialize_compute_diagnostic_variables(in_diagv_node)
  ! **** Initialize operations for the compute_diagnostic_variables subroutine ****
  type(in_type_diagv_node),intent(out) :: in_diagv_node  ! input data object for diagv_node
  ! interface local name space to input data object for diagv_node
  call in_diagv_node % initialize(iSoil,in_soilLiqFlux,diag_data,mpar_data,flux_data)
 end subroutine initialize_compute_diagnostic_variables

 subroutine update_compute_diagnostic_variables(in_diagv_node,out_diagv_node)
  ! **** Update operations for the compute_diagnostic_variables subroutine ****
  type(in_type_diagv_node) ,intent(in)  :: in_diagv_node  ! input data object for diagv_node
  type(out_type_diagv_node),intent(out) :: out_diagv_node ! output data object for diagv_node
  ! compute diagnostic variables
  call diagv_node(in_diagv_node,out_diagv_node)
 end subroutine update_compute_diagnostic_variables

 subroutine finalize_compute_diagnostic_variables(out_diagv_node)
  ! **** Finalize operations for the compute_diagnostic_variables subroutine ****
  type(out_type_diagv_node),intent(in) :: out_diagv_node ! output data object for diagv_node
  ! interface output data object for diagv_node to local name space
  associate(&
   err          => out_soilLiqFlux % err,     & ! error code
   message      => out_soilLiqFlux % cmessage & ! error message
  &)
   call out_diagv_node % finalize(iSoil,nSoil,io_soilLiqFlux,iceImpedeFac,&
                                  &dHydCond_dTemp,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end associate
 end subroutine finalize_compute_diagnostic_variables

 subroutine compute_surface_infiltration
  ! **** compute infiltration at the surface and its derivative w.r.t. mass in the upper soil layer ****
  ! NOTE: this needs to change if nLake>0
  type(in_type_surfaceFlux)  ::  in_surfaceFlux
  type(io_type_surfaceFlux)  ::  io_surfaceFlux
  type(out_type_surfaceFlux) :: out_surfaceFlux

  call initialize_compute_surface_infiltration(in_surfaceFlux,io_surfaceFlux)

  call update_compute_surface_infiltration(in_surfaceFlux,io_surfaceFlux,out_surfaceFlux)

  call finalize_compute_surface_infiltration(io_surfaceFlux,out_surfaceFlux); if (return_flag) return

 end subroutine compute_surface_infiltration

 subroutine initialize_compute_surface_infiltration(in_surfaceFlux,io_surfaceFlux)
  ! **** Initialize operations for compute_surface_infiltration ****
  type(in_type_surfaceFlux),intent(out) :: in_surfaceFlux
  type(io_type_surfaceFlux),intent(out) :: io_surfaceFlux
  ! set derivative w.r.t. state above to zero (does not exist)
  associate(&
   ! intent(inout): flux derivatives ...
   dq_dHydStateAbove => io_soilLiqFlux % dq_dHydStateAbove,& ! ... in layer interfaces w.r.t. state variables in the layer above
   dq_dNrgStateAbove => io_soilLiqFlux % dq_dNrgStateAbove & ! ... w.r.t. temperature in the layer above (m s-1 K-1)
  &)
   dq_dHydStateAbove(0) = 0._rkind
   dq_dNrgStateAbove(0) = 0._rkind
  end associate

  ! compute surface flux and its derivative...
  call in_surfaceFlux % initialize(nRoots,ixIce,nSoil,nGlce,ibeg,iend,in_soilLiqFlux,io_soilLiqFlux,&
                                 &model_decisions,prog_data,mpar_data,flux_data,diag_data,&
                                 &iLayerHeight,dHydCond_dTemp,iceImpedeFac)
  call io_surfaceFlux % initialize(nSoil,io_soilLiqFlux,iLayerHydCond)
 end subroutine initialize_compute_surface_infiltration

 subroutine update_compute_surface_infiltration(in_surfaceFlux,io_surfaceFlux,out_surfaceFlux)
  ! **** Update operations for compute_surface_infiltration ****
  type(in_type_surfaceFlux) ,intent(in)    :: in_surfaceFlux
  type(io_type_surfaceFlux) ,intent(inout) :: io_surfaceFlux
  type(out_type_surfaceFlux),intent(out)   :: out_surfaceFlux
  call surfaceFlux(io_soilLiqFlux,in_surfaceFlux,io_surfaceFlux,out_surfaceFlux)
 end subroutine update_compute_surface_infiltration

 subroutine finalize_compute_surface_infiltration(io_surfaceFlux,out_surfaceFlux)
  ! **** Finalize operations for compute_surface_infiltration ****
  type(io_type_surfaceFlux) ,intent(in) :: io_surfaceFlux
  type(out_type_surfaceFlux),intent(in) :: out_surfaceFlux

  ! interface object data components with local name space
  call io_surfaceFlux % finalize(nSoil,io_soilLiqFlux,iLayerHydCond)
  associate(&
   err     => out_soilLiqFlux % err,     & ! error code
   message => out_soilLiqFlux % cmessage & ! error message
  &)
   call out_surfaceFlux % finalize(io_soilLiqFlux,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end associate

  ! include base soil evaporation as the upper boundary flux
  associate(&
   iLayerLiqFluxSoil         => io_soilLiqFlux % iLayerLiqFluxSoil,      & ! liquid flux at soil layer interfaces (m s-1)
   scalarGroundEvaporation   => in_soilLiqFlux % scalarGroundEvaporation,& ! ground evaporation (kg m-2 s-1)
   scalarSurfaceInfiltration => io_soilLiqFlux % scalarInfiltration,     & ! surface infiltration rate (m s-1)
   dq_dHydStateBelow         => io_soilLiqFlux % dq_dHydStateBelow,      & ! derivative in the flux in layer interfaces w.r.t. state variables in the layer below
   dq_dNrgStateBelow         => io_soilLiqFlux % dq_dNrgStateBelow       & ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
  &)
   iLayerLiqFluxSoil(0) = scalarGroundEvaporation/iden_water + scalarSurfaceInfiltration

   dq_dHydStateBelow(0) = 0._rkind ! contribution will be in dq_dHydStateLayerSurfVec(1)
   dq_dNrgStateBelow(0) = 0._rkind ! contribution will be in dq_dNrgStateLayerSurfVec(1)
  end associate
 end subroutine finalize_compute_surface_infiltration

 subroutine compute_interface_fluxes_derivatives
  ! **** compute fluxes and derivatives at layer interfaces ****
  type(in_type_iLayerFlux)  :: in_iLayerFlux  ! input data object for iLayerFlux
  type(out_type_iLayerFlux) :: out_iLayerFlux ! output data object for iLayerFlux

  ! computing flux at the bottom of the layer
  do iLayer=ixTop,min(ixBot,nSoil-1)

   call initialize_compute_interface_fluxes_derivatives(in_iLayerFlux)

   call update_compute_interface_fluxes_derivatives(in_iLayerFlux,out_iLayerFlux)

   call finalize_compute_interface_fluxes_derivatives(out_iLayerFlux); if (return_flag) return

  end do
 end subroutine compute_interface_fluxes_derivatives

 subroutine initialize_compute_interface_fluxes_derivatives(in_iLayerFlux)
  ! **** Initialize operations for compute_interface_fluxes_derivatives subroutine ****
  type(in_type_iLayerFlux),intent(out) :: in_iLayerFlux  ! input data object for iLayerFlux
  ! interface local name space to iLayerFlux input object
  call in_iLayerFlux % initialize(iLayer,nSoil,ibeg,iend,in_soilLiqFlux,io_soilLiqFlux,&
                                 &prog_data,dHydCond_dTemp)
 end subroutine initialize_compute_interface_fluxes_derivatives

 subroutine update_compute_interface_fluxes_derivatives(in_iLayerFlux,out_iLayerFlux)
  ! **** Update operations for compute_interface_fluxes_derivatives subroutine ****
  type(in_type_iLayerFlux) ,intent(in)  :: in_iLayerFlux  ! input data object for iLayerFlux
  type(out_type_iLayerFlux),intent(out) :: out_iLayerFlux ! output data object for iLayerFlux
  ! compute fluxes at layer interface
  call iLayerFlux(in_iLayerFlux,out_iLayerFlux)
 end subroutine update_compute_interface_fluxes_derivatives

 subroutine finalize_compute_interface_fluxes_derivatives(out_iLayerFlux)
  ! **** Finalize operations for compute_interface_fluxes_derivatives subroutine
  type(out_type_iLayerFlux),intent(in) :: out_iLayerFlux ! output data object for iLayerFlux
  ! interface iLayerFlux output object to local name space
  associate(&
   err     => out_soilLiqFlux % err,                       & ! error code
   message => out_soilLiqFlux % cmessage                   & ! error message
  &)
   call out_iLayerFlux % finalize(iLayer,nSoil,io_soilLiqFlux,iLayerHydCond,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end associate
 end subroutine finalize_compute_interface_fluxes_derivatives

 subroutine compute_drainage_flux
  ! **** Compute the drainage flux from the bottom of the soil profile and its derivative ****
  type(in_type_qDrainFlux)  :: in_qDrainFlux
  type(io_type_qDrainFlux)  :: io_qDrainFlux
  type(out_type_qDrainFlux) :: out_qDrainFlux

  call initialize_compute_drainage_flux(in_qDrainFlux,io_qDrainFlux)

  call update_compute_drainage_flux(in_qDrainFlux,io_qDrainFlux,out_qDrainFlux)
 
  call finalize_compute_drainage_flux(io_qDrainFlux,out_qDrainFlux); if (return_flag) return

 end subroutine compute_drainage_flux

 subroutine initialize_compute_drainage_flux(in_qDrainFlux,io_qDrainFlux)
  ! **** Initialize operations for compute_drainage_flux ****
  type(in_type_qDrainFlux),intent(out) :: in_qDrainFlux
  type(io_type_qDrainFlux),intent(out) :: io_qDrainFlux
  call in_qDrainFlux % initialize(nSoil,nGlce,ibeg,iend,in_soilLiqFlux,io_soilLiqFlux,model_decisions,&
                                 &prog_data,mpar_data,flux_data,diag_data,iceImpedeFac,&
                                 &dHydCond_dTemp)
  call io_qDrainFlux % initialize(io_soilLiqFlux)
 end subroutine initialize_compute_drainage_flux

subroutine update_compute_drainage_flux(in_qDrainFlux,io_qDrainFlux,out_qDrainFlux)
  ! **** Update operations for compute_drainage_flux ****
  type(in_type_qDrainFlux) ,intent(in)   :: in_qDrainFlux
  type(io_type_qDrainFlux),intent(inout) :: io_qDrainFlux
  type(out_type_qDrainFlux),intent(out)  :: out_qDrainFlux
  call qDrainFlux(in_qDrainFlux,io_qDrainFlux,out_qDrainFlux)
 end subroutine update_compute_drainage_flux

 subroutine finalize_compute_drainage_flux(io_qDrainFlux,out_qDrainFlux)
  ! **** finalize operations for compute_drainage_flux ****
  type(io_type_qDrainFlux),intent(inout) :: io_qDrainFlux
  type(out_type_qDrainFlux),intent(in) :: out_qDrainFlux

  ! interface object data components with local name space
  call io_qDrainFlux % finalize(io_soilLiqFlux)
  associate(&
   err     => out_soilLiqFlux % err,                       & ! error code
   message => out_soilLiqFlux % cmessage                   & ! error message
  &)
   call out_qDrainFlux % finalize(nSoil,io_soilLiqFlux,iLayerHydCond,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end associate

  ! no dependence on the aquifer for drainage currently, but keep this here in case we want to couple some day
  ! NOTE: glacier coupling terms are computed in Jacobian directly
  associate(&
   ! derivatives in flux w.r.t. ...
   dq_dHydStateBelow => io_soilLiqFlux % dq_dHydStateBelow,& ! ... hydrology state variables in the layer below
   dq_dNrgStateBelow => io_soilLiqFlux % dq_dNrgStateBelow & ! ... temperature in the layer below (m s-1 K-1)
  &)
   dq_dHydStateBelow(nSoil) = 0._rkind  ! keep this here in case we want to couple some day....
   dq_dNrgStateBelow(nSoil) = 0._rkind  ! keep this here in case we want to couple some day....
  end associate
 end subroutine finalize_compute_drainage_flux
end subroutine soilLiqFlux

! ***************************************************************************************************************
! private subroutine diagv_node: compute transmittance and derivatives for model nodes
! ***************************************************************************************************************
subroutine diagv_node(in_diagv_node,out_diagv_node)
  USE soil_utils_module,only:iceImpede            ! compute the ice impedence factor
  USE soil_utils_module,only:volFracLiq           ! compute volumetric fraction of liquid water as a function of matric head
  USE soil_utils_module,only:hydCond_psi          ! compute hydraulic conductivity as a function of matric head
  USE soil_utils_module,only:hydCondMP_liq        ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
  USE soil_utils_module,only:dTheta_dPsi          ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
  USE soil_utils_module,only:dPsi_dTheta          ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  USE soil_utils_module,only:dHydCond_dPsi        ! compute derivative in hydraulic conductivity w.r.t. matric head (s-1)
  USE soil_utils_module,only:dIceImpede_dTemp     ! compute the derivative in the ice impedance factor w.r.t. temperature (K-1)
  ! compute hydraulic transmittance and derivatives for all layers
  implicit none
  ! input: model control, variables, derivatives, and parameters
  type(in_type_diagv_node),  intent(in)  :: in_diagv_node
  ! output: characteristic derivatives, transmittance variables, and error control
  type(out_type_diagv_node), intent(out) :: out_diagv_node
  ! local variables
  real(rkind)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
  real(rkind)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
  real(rkind)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
  real(rkind)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
  real(rkind)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
  real(rkind)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
  real(rkind)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
  real(rkind)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
  real(rkind)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
  real(rkind)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
  real(rkind)                      :: relSatMP                  ! relative saturation of macropores (-)
  logical(lgt)                     :: return_flag               ! flag for return statements

    call initialize_diagv_node

    call update_diagv_node;   if (return_flag) return

    call finalize_diagv_node; if (return_flag) return

contains

 subroutine initialize_diagv_node
  ! **** Initialize operations for diagv_node ****
  ! initialize error control
  return_flag=.false.
  associate(&
   err     => out_diagv_node % err    , & ! error code
   message => out_diagv_node % message  & ! error message
  &)
   err=0; message="diagv_node/"
  end associate
 end subroutine initialize_diagv_node

 subroutine update_diagv_node
  ! **** Update operations for diagv_node ****

   call update_diagv_node_characteristic_derivatives; if (return_flag) return

   call update_diagv_node_hydraulic_conductivity;     if (return_flag) return

 end subroutine update_diagv_node

 subroutine update_diagv_node_characteristic_derivatives
  ! **** Update operations for diagv_node: compute characteristic derivatives ****
  ! compute the derivative in the soil water characteristic
  associate(&
   ! input: state and diagnostic variables
   scalarMatricHeadLiqTrial => in_diagv_node % scalarMatricHeadLiqTrial, & ! liquid matric head in each layer (m)
   scalarVolFracLiqTrial    => in_diagv_node % scalarVolFracLiqTrial   , & ! volumetric fraction of liquid water in a given layer (-)
   ! input: soil parameters
   vGn_alpha => in_diagv_node % vGn_alpha, & ! van Genuchten "alpha" parameter (m-1)
   vGn_n     => in_diagv_node % vGn_n    , & ! van Genuchten "n" parameter (-)
   vGn_m     => in_diagv_node % vGn_m    , & ! van Genuchten "m" parameter (-)
   mpExp     => in_diagv_node % mpExp    , & ! empirical exponent in macropore flow equation (-)
   theta_sat => in_diagv_node % theta_sat, & ! soil porosity (-)
   theta_res => in_diagv_node % theta_res, & ! soil residual volumetric water content (-)
   ! output: derivative in the soil water characteristic
   scalardPsi_dTheta => out_diagv_node % scalardPsi_dTheta, & ! derivative in the soil water characteristic
   scalardTheta_dPsi => out_diagv_node % scalardTheta_dPsi, & ! derivative in the soil water characteristic
   ! output: error control
   err     => out_diagv_node % err    , & ! error code
   message => out_diagv_node % message  & ! error message
  &)

   scalardTheta_dPsi = dTheta_dPsi(scalarMatricHeadLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)

  end associate
 end subroutine update_diagv_node_characteristic_derivatives

 subroutine update_diagv_node_hydraulic_conductivity
  ! **** Update operations for diagv_node: compute hydraulic conductivity and derivatives ****
  ! compute hydraulic conductivity and its derivative in each soil layer
  associate(&
   scalarVolFracIceTrial    => in_diagv_node % scalarVolFracIceTrial, & ! volumetric fraction of ice in a given layer (-)
   f_impede  => in_diagv_node % f_impede,                             & ! ice impedence factor (-)
   iceImpedeFac  => out_diagv_node % iceImpedeFac                     & ! ice impedence factor in each layer (-)
  &)
   ! compute the ice impedence factor and its derivative w.r.t. volumetric liquid water content (-)
   call iceImpede(scalarVolFracIceTrial,f_impede, &  ! input
                   iceImpedeFac,dIceImpede_dLiq)     ! output
  end associate

  call update_diagv_node_hydraulic_conductivity_mixed_form; if (return_flag) return
 end subroutine update_diagv_node_hydraulic_conductivity

 subroutine update_diagv_node_hydraulic_conductivity_mixed_form
  ! **** Update operations for diagv_node: compute hydraulic conductivity and derivatives for mixed form of Richards' equation ****
  associate(&
   ! input: state and diagnostic variables
   scalarMatricHeadLiqTrial => in_diagv_node % scalarMatricHeadLiqTrial, & ! liquid matric head in each layer (m)
   scalarVolFracIceTrial    => in_diagv_node % scalarVolFracIceTrial   , & ! volumetric fraction of ice in a given layer (-)
   ! input: pre-computed derivatives
   dTheta_dTk    => in_diagv_node % dTheta_dTk   , & ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   dPsiLiq_dTemp => in_diagv_node % dPsiLiq_dTemp, & ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: soil parameters
   vGn_alpha => in_diagv_node % vGn_alpha, & ! van Genuchten "alpha" parameter (m-1)
   vGn_n     => in_diagv_node % vGn_n    , & ! van Genuchten "n" parameter (-)
   vGn_m     => in_diagv_node % vGn_m    , & ! van Genuchten "m" parameter (-)
   mpExp     => in_diagv_node % mpExp    , & ! empirical exponent in macropore flow equation (-)
   theta_sat => in_diagv_node % theta_sat, & ! soil porosity (-)
   theta_res => in_diagv_node % theta_res, & ! soil residual volumetric water content (-)
   theta_mp  => in_diagv_node % theta_mp , & ! volumetric liquid water content when macropore flow begins (-)
   f_impede  => in_diagv_node % f_impede , & ! ice impedence factor (-)
   ! input: saturated hydraulic conductivity ...
   scalarSatHydCond   => in_diagv_node % scalarSatHydCond,  & ! ... at the mid-point of a given layer (m s-1)
   scalarSatHydCondMP => in_diagv_node % scalarSatHydCondMP,& ! ... of macropores at the mid-point of a given layer (m s-1)
   ! output: derivative in the soil water characteristic
   scalardTheta_dPsi => out_diagv_node % scalardTheta_dPsi, & ! derivative in the soil water characteristic
   ! output: transmittance
   scalarHydCond => out_diagv_node % scalarHydCond, & ! hydraulic conductivity at layer mid-points (m s-1)
   scalarDiffuse => out_diagv_node % scalarDiffuse, & ! diffusivity at layer mid-points (m2 s-1)
   iceImpedeFac  => out_diagv_node % iceImpedeFac , & ! ice impedence factor in each layer (-)
   ! output: transmittance derivatives in ...
   dHydCond_dMatric => out_diagv_node % dHydCond_dMatric, & ! ... hydraulic conductivity w.r.t matric head (s-1)
   dHydCond_dTemp   => out_diagv_node % dHydCond_dTemp    & ! ... hydraulic conductivity w.r.t temperature (m s-1 K-1)
  &)

   ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
   hydCond_noIce = hydCond_psi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m)
   scalarDiffuse = realMissing ! not used, so cause problems
   ! compute the hydraulic conductivity of macropores (m s-1)
   localVolFracLiq = volFracLiq(scalarMatricHeadLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   scalarHydCondMP = hydCondMP_liq(localVolFracLiq,theta_sat,theta_mp,mpExp,scalarSatHydCondMP,scalarSatHydCond)
   scalarHydCond   = hydCond_noIce*iceImpedeFac + scalarHydCondMP

   ! compute derivative in hydraulic conductivity (m s-1)
   ! compute derivative for macropores
   if (localVolFracLiq > theta_mp) then
     relSatMP              = (localVolFracLiq - theta_mp)/(theta_sat - theta_mp)
     dHydCondMacro_dVolLiq = ((scalarSatHydCondMP - scalarSatHydCond)/(theta_sat - theta_mp))*mpExp*(relSatMP**(mpExp - 1._rkind))
     dHydCondMacro_dMatric = scalardTheta_dPsi*dHydCondMacro_dVolLiq
   else
     dHydCondMacro_dVolLiq = 0._rkind
     dHydCondMacro_dMatric = 0._rkind
   end if
   ! compute derivatives for micropores
   if (scalarVolFracIceTrial > verySmaller) then
     dK_dPsi__noIce        = dHydCond_dPsi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m)
     dHydCondMicro_dTemp   = dPsiLiq_dTemp*dK_dPsi__noIce  ! m s-1 K-1
     dHydCondMicro_dMatric = hydCond_noIce*dIceImpede_dLiq*scalardTheta_dPsi + dK_dPsi__noIce*iceImpedeFac
   else
     dHydCondMicro_dTemp   = 0._rkind
     dHydCondMicro_dMatric = dHydCond_dPsi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m)
   end if
   ! combine matric derivatives
   dHydCond_dMatric = dHydCondMicro_dMatric + dHydCondMacro_dMatric

   ! compute analytical derivative for change in ice impedance factor w.r.t. temperature
   call dIceImpede_dTemp(scalarVolFracIceTrial, & ! intent(in):  trial value of volumetric ice content (-)
                         dTheta_dTk,            & ! intent(in):  derivative in volumetric liquid water content w.r.t. temperature (K-1)
                         f_impede,              & ! intent(in):  ice impedance parameter (-)
                         dIceImpede_dT          ) ! intent(out): derivative in ice impedance factor w.r.t. temperature (K-1)
   ! compute derivative in hydraulic conductivity w.r.t. temperature
   dHydCond_dTemp = hydCond_noIce*dIceImpede_dT + dHydCondMicro_dTemp*iceImpedeFac

  end associate
 end subroutine update_diagv_node_hydraulic_conductivity_mixed_form

 subroutine finalize_diagv_node
  associate(&
   err     => out_diagv_node % err    , & ! error code
   message => out_diagv_node % message  & ! error message
  &)
   ! final error check
   if(err/=0)then; message=trim(message)//'unanticipated error in diagv_node'; return_flag=.true.; return; end if
  end associate
 end subroutine finalize_diagv_node

end subroutine diagv_node

! ***************************************************************************************************************
! private subroutine surfaceFlux: compute the surface flux and its derivative
! ***************************************************************************************************************
subroutine surfaceFlux(io_soilLiqFlux,in_surfaceFlux,io_surfaceFlux,out_surfaceFlux)
  USE soil_utils_module,only:volFracLiq            ! compute volumetric fraction of liquid water as a function of matric head (-)
  USE soil_utils_module,only:hydCond_psi           ! compute hydraulic conductivity as a function of matric head (m s-1)
  USE soil_utils_module,only:crit_soilT            ! compute critical temperature below which ice exists
  USE soil_utils_module,only:gammp,gammp_complex   ! compute the regularized lower incomplete Gamma function
  ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
  implicit none
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! input: use soilLiqFlux object for array dimensions
  type(io_type_soilLiqFlux) ,intent(in)    :: io_soilLiqFlux          ! input-output object for soilLiqFlux
  ! input: model control, variables, derivatives, soil layer depth, boundary conditions, fluxes, and transmittance and soil parameters
  type(in_type_surfaceFlux) ,intent(in)    :: in_surfaceFlux          ! input object for surfaceFlux
  ! input-output: hydraulic conductivity and diffusivity, and infiltration parameters
  type(io_type_surfaceFlux) ,intent(inout) :: io_surfaceFlux          ! input object for surfaceFlux
  ! output: runoff, infiltration, derivatives, and error control
  type(out_type_surfaceFlux),intent(out)   :: out_surfaceFlux         ! output object for surfaceFlux
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! general
  integer(i4b)                     :: iLayer                              ! index of soil layer
  real(rkind)                      :: Tcrit                               ! temperature where all water is unfrozen (K)
  real(rkind)                      :: fPart1,fPart2                       ! different parts of a function
  real(rkind)                      :: dPart1(1:in_surfaceFlux % nSoil)    ! derivatives for different parts of a function
  real(rkind)                      :: dPart2(1:in_surfaceFlux % nSoil)    ! derivatives for different parts of a function
  real(rkind)                      :: dfracCap(1:in_surfaceFlux % nSoil)  ! derivatives for different parts of a function
  real(rkind)                      :: dfInfRaw(1:in_surfaceFlux % nSoil)  ! derivatives for different parts of a function
  real(rkind)                      :: total_soil_depth                    ! total depth of soil (m)
  integer(i4b)                     :: ixInfRateMax_use                    ! topmodel_GA choice of the maximum infiltration rate for glacier domains
  integer(i4b)                     :: ix_hc_profile_use                   ! hydraulic conductivity profile used to shape the topmodel_GA infiltration rate
  ! head boundary condition
  real(rkind)                      :: cFlux                               ! capillary flux (m s-1)
  ! simplified Green-Ampt infiltration
  real(rkind)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
  real(rkind)                      :: rootZoneIce                         ! depth of ice in the root zone (m)
  real(rkind)                      :: availCapacity                       ! available storage capacity in the root zone (m)
  real(rkind)                      :: depthWettingFront                   ! depth to the wetting front (m)
  real(rkind)                      :: hydCondWettingFront                 ! hydraulic conductivity at the wetting front (m s-1)
  real(rkind)                      :: dHydCondWF_dDepth                   ! derivative in hydraulic conductivity at the wetting front w.r.t. its depth (s-1)
  real(rkind)                      :: refDepth_use                        ! reference depth for the power-law profile scaling (m)
  real(rkind)                      :: depthWF_use                         ! wetting-front depth, clamped off the base of the soil (m)
  ! saturated area associated with variable storage capacity
  real(rkind)                      :: fracCap                             ! fraction of pore space filled with liquid water and ice (-)
  real(rkind)                      :: fInfRaw                             ! infiltrating area before imposing solution constraints (-)
  real(rkind),parameter            :: maxFracCap=0.995_rkind              ! maximum fraction capacity -- used to avoid numerical problems associated with an enormous derivative
  real(rkind),parameter            :: scaleFactor=0.000001_rkind          ! scale factor for the smoothing function (-)
  real(rkind),parameter            :: qSurfScaleMax=1000._rkind           ! maximum surface runoff scaling factor (-)
  ! fraction of impermeable area associated with frozen ground
  real(rkind)                      :: alpha                               ! shape parameter in the Gamma distribution
  real(rkind)                      :: xLimg                               ! upper limit of the integral
  ! FUSE
  real(rkind),parameter            :: alpha_LSE=1.e3_rkind                ! smoothness parameter for LSE smoother function
  real(rkind),parameter            :: roundoff_tolerance = 1.e2_rkind * epsilon(1._rkind) ! tolerance for round-off error is near machine epsilon 
  real(rkind)                      :: S1                                  ! total water content in upper soil layer (m)
  real(rkind)                      :: S1_max                              ! Maximum storage in the upper layer (m)
  ! FUSE PRMS variables
  real(rkind)                      :: phi_tens                            ! fraction of total storage as tension storage (m)
  real(rkind)                      :: SatArea_max                         ! maximum saturated area (-)
  real(rkind)                      :: S1_T                                ! tension water content in upper soil layer (m)
  real(rkind)                      :: S1_T_max                            ! maximum tension water content in upper soil layer (m)
  ! FUSE ARNO/VIC variables
  logical(lgt),parameter :: smoother = .true.                             ! control for optional smoothing in base variable
  real(rkind)                      :: base                                ! base used in saturated area formula (-) ARNO/VIC 
  real(rkind)                      :: b_arnovic                           ! ARNO/VIC exponent (-) 
  real(rkind)                      :: S1_star                             ! total water content in upper FUSE layer computed with a smoothed min (m)
  ! FUSE TOPMODEL variables
  real(rkind)                      :: alpha_topmodel                      ! gamma shape
  real(rkind)                      :: chi_topmodel                        ! gamma scale
  real(rkind)                      :: lambda                              ! mean for alpha_topmodel
  real(rkind)                      :: mu                                  ! offset for alpha_topmodel
  real(rkind)                      :: x_crit                              ! critical x (random variable) value
  real(rkind),parameter            :: zeta_upper=1.e3_rkind               ! upper limit of integral (approaches infinity, but ~1000 provides an accurate result) 
  real(rkind)                      :: zeta_crit                           ! critical topographic index value (log space)
  real(rkind)                      :: zeta_crit_n                         ! critical topographic index value (power-transformed)
  real(rkind)                      :: n_topmodel                          ! TOPMODEL exponent exponent (must be sufficiently large to avoid divergence of lambda_n -- n>=3.5 or so)
  complex(rkind)                   :: lambda_n                            ! mean of the power-transformed topographic index
  ! derivatives
  real(rkind) :: dVolFracLiq_dWat(1:in_surfaceFlux % nSoil)  ! ... vol fraction of liquid w.r.t. water state variable in soil layers
  real(rkind) :: dVolFracIce_dWat(1:in_surfaceFlux % nSoil)  ! ... vol fraction of ice w.r.t. water state variable in soil layers
  real(rkind) :: dVolFracLiq_dTk(1:in_surfaceFlux % nSoil)   ! ... vol fraction of liquid w.r.t. temperature in soil layers
  real(rkind) :: dVolFracIce_dTk(1:in_surfaceFlux % nSoil)   ! ... vol fraction of ice w.r.t. temperature in soil layers
  real(rkind) :: dRootZoneLiq_dWat(1:in_surfaceFlux % nSoil) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
  real(rkind) :: dRootZoneIce_dWat(1:in_surfaceFlux % nSoil) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
  real(rkind) :: dRootZoneLiq_dTk(1:in_surfaceFlux % nSoil)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
  real(rkind) :: dRootZoneIce_dTk(1:in_surfaceFlux % nSoil)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
  real(rkind) :: dDepthWettingFront_dWat(1:in_surfaceFlux % nSoil) ! ... scalar depth of wetting front w.r.t. water state variable in root layers
  real(rkind) :: dDepthWettingFront_dTk(1:in_surfaceFlux % nSoil)  ! ... scalar depth of wetting front w.r.t. temperature in root layers
  real(rkind) :: dxMaxInfilRate_dWat(1:in_surfaceFlux % nSoil) ! ... scalar max infiltration rate w.r.t. water state variable in root layers
  real(rkind) :: dxMaxInfilRate_dTk(1:in_surfaceFlux % nSoil)  ! ... scalar max infiltration rate w.r.t. temperature in root layers
  real(rkind) :: dInfilArea_dWat(1:in_surfaceFlux % nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in soil layers
  real(rkind) :: dInfilArea_dTk(1:in_surfaceFlux % nSoil)   ! ... scalar infiltration rate w.r.t. temperature in soil layers
  real(rkind) :: dFrozenArea_dWat(1:in_surfaceFlux % nSoil) ! ... scalar frozen area w.r.t. water state variable in soil layers
  real(rkind) :: dFrozenArea_dTk(1:in_surfaceFlux % nSoil)  ! ... scalar frozen area w.r.t. temperature in soil layers
  real(rkind) :: dInfilRate_dWat(1:in_surfaceFlux % nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in soil layers
  real(rkind) :: dInfilRate_dTk(1:in_surfaceFlux % nSoil)   ! ... scalar infiltration rate w.r.t. temperature in soil layers
  ! error control
  logical(lgt) :: return_flag ! logical flag for return statements

  call initialize_surfaceFlux

  call update_surfaceFlux;   if (return_flag) return

  call finalize_surfaceFlux; if (return_flag) return

contains

 subroutine initialize_surfaceFlux
  ! **** Initialize operations for surfaceFlux ****
  ! allocate output object array components
  out_surfaceFlux % dq_dHydStateVec = io_soilLiqFlux % dq_dHydStateLayerSurfVec
  out_surfaceFlux % dq_dNrgStateVec = io_soilLiqFlux % dq_dNrgStateLayerSurfVec

  ! initialize error control
  return_flag=.false.
  associate(&
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)
   err=0; message="surfaceFlux/"
  end associate

  ! initialize derivatives
  associate(&
   ! output: derivatives in surface infiltration w.r.t. ...
   dq_dHydStateVec => out_surfaceFlux % dq_dHydStateVec, & ! ... hydrology state in every soil layer (m s-1 or s-1)
   dq_dNrgStateVec => out_surfaceFlux % dq_dNrgStateVec  & ! ... energy state in every soil layer (m s-1 K-1)
  &)
   dVolFracLiq_dWat(:)    = 0._rkind
   dVolFracIce_dWat(:)    = 0._rkind
   dVolFracLiq_dTk(:)     = 0._rkind
   dVolFracIce_dTk(:)     = 0._rkind
   dInfilArea_dWat(:)     = 0._rkind
   dInfilArea_dTk(:)      = 0._rkind
   dInfilRate_dWat(:)     = 0._rkind
   dInfilRate_dTk(:)      = 0._rkind
   dxMaxInfilRate_dWat(:) = 0._rkind
   dxMaxInfilRate_dTk(:)  = 0._rkind
   dFrozenArea_dWat(:)    = 0._rkind
   dFrozenArea_dTk(:)     = 0._rkind
   dq_dHydStateVec(:)     = 0._rkind
   dq_dNrgStateVec(:)     = 0._rkind ! energy state variable is temperature (transformed outside soilLiqFlux_module if needed)
  end associate

  ! initialize runoff values
  associate(&
   scalarSurfaceRunoff       => out_surfaceFlux % scalarSurfaceRunoff       , & ! surface runoff (m s-1)
   scalarSurfaceRunoff_IE    => out_surfaceFlux % scalarSurfaceRunoff_IE    , & ! infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE    => out_surfaceFlux % scalarSurfaceRunoff_SE      & ! saturation excess surface runoff (m s-1)
  &)
   scalarSurfaceRunoff       = 0._rkind 
   scalarSurfaceRunoff_IE    = 0._rkind  
   scalarSurfaceRunoff_SE    = 0._rkind  
  end associate

  ! total soil depth is needed for infiltration area methods, so compute here
  associate(&
   mLayerDepth  => in_surfaceFlux % mLayerDepth      & ! depth of soil layers (m) 
  &)
   total_soil_depth = sum(mLayerDepth)
  end associate

 end subroutine initialize_surfaceFlux

 subroutine update_surfaceFlux
  ! **** Update operations for surfaceFlux ****
  associate(&
   ! input: model control
   firstSplitOper => in_surfaceFlux % firstSplitOper, & ! flag indicating if desire to compute infiltration
   bc_upper       => in_surfaceFlux % bc_upper,       & ! index defining the type of boundary conditions
   ixInfRateMax   => in_surfaceFlux % ixInfRateMax,   & ! index defining the maximum infiltration rate method
   ix_hc_profile  => in_surfaceFlux % ix_hc_profile,  & ! index defining the hydraulic conductivity profile
   nGlce          => in_surfaceFlux % nGlce,          & ! number of glacier ice layers
   surfRun_SE     => in_surfaceFlux % surfRun_SE,     & ! index defining the saturation excess surface runoff method
   ! input to compute infiltration
   scalarRainPlusMelt => in_surfaceFlux % scalarRainPlusMelt, & ! rain plus melt  (m s-1)
   ! input output: infiltration
   scalarSurfaceInfiltration => io_surfaceFlux % scalarSurfaceInfiltration, & ! surface infiltration (m s-1)
   xMaxInfilRate             => io_surfaceFlux % xMaxInfilRate,             & ! maximum infiltration rate (m s-1)
   scalarSaturatedArea       => io_surfaceFlux % scalarSaturatedArea,       & ! fraction of area that is saturated (-)
   ! output: runoff 
   scalarSurfaceRunoff       => out_surfaceFlux % scalarSurfaceRunoff,      & ! surface runoff (m s-1)
   scalarSurfaceRunoff_IE    => out_surfaceFlux % scalarSurfaceRunoff_IE,   & ! infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE    => out_surfaceFlux % scalarSurfaceRunoff_SE,   & ! saturation excess surface runoff (m s-1)
   ! output: error control
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)
   ixInfRateMax_use = ixInfRateMax
   if(nGlce>0) ixInfRateMax_use = topmodel_GA ! glacier debris conductivity varies with depth, so the wetting-front conductivity must follow that profile
   ! the topmodel_GA infiltration rate re-derives the conductivity profile, so it has to use the same shape satHydCond used
   ix_hc_profile_use = ix_hc_profile
   if(nGlce>0) ix_hc_profile_use = expLaw_profile ! must match the override in satHydCond

   ! compute the surface flux and its derivative
   if (firstSplitOper .or. updateInfil) then
     select case(bc_upper)
       case(prescribedHead) ! head condition, no frozen area and all area infiltrates (for glacier debris no melt infiltrates)
         call update_surfaceFlux_prescribedHead; if (return_flag) return 
 
       case(liquidFlux)     ! flux condition
         ! compute volumetric fraction of liquid and ice water in each soil layer and their derivatives
         if(updateInfil) call update_volFracLiq_derivatives; if (return_flag) return

         ! Get infiltration area not considering frozen area, based on SE method
         select case(surfRun_SE) ! saturation excess surface runoff method, sets infiltration area (not considering frozen) and its derivatives
           case(zero_SE)         ! zero saturation excess surface runoff, all area infiltrates if not frozen
            io_surfaceFlux % scalarInfilArea = 1._rkind 
           case(homegrown_SE)    ! homegrown saturation excess surface runoff (original SUMMA method)
              call update_surfaceFlux_homegrown_infilArea;     if (return_flag) return
           case(FUSEPRMS)        ! FUSE PRMS surface runoff
             call update_surfaceFlux_FUSE_PRMS_infilArea;      if (return_flag) return
           case(FUSEAVIC)        ! FUSE ARNO/VIC surface runoff
             call update_surfaceFlux_FUSE_ARNO_VIC_infilArea;  if (return_flag) return
           case(FUSETOPM)        ! FUSE TOPMODEL surface runoff
             call update_surfaceFlux_FUSE_TOPMODEL_infilArea;  if (return_flag) return
           case default; err=20; message=trim(message)//'unknown saturation excess surface runoff method'; return_flag=.true.; return
         end select

         ! Calculate maximum infiltration rate and scalarFrozenArea (and their derivatives if needed)
         select case(ixInfRateMax_use)       ! maximum infiltration rate method (controls infiltration excess surface runoff)
           case(noInfiltrationExcess)    ! zero infiltration excess surface runoff
             call update_surfaceFlux_liquidFlux_noinfratemax
           case(GreenAmpt, topmodel_GA)  ! infiltration excess runoff possible
             call update_surfaceFlux_liquidFlux_calculate_infratemax;  if (return_flag) return
           case default; err=20; message=trim(message)//'unknown infiltration excess surface runoff method'; return_flag=.true.; return
         end select

         ! Compute total infiltration and derivatives
         call update_surfaceFlux_liquidFlux_infiltration;  if (return_flag) return
 
       case default; err=20; message=trim(message)//'unknown upper boundary condition for soil hydrology'; return_flag=.true.; return ! end of select of bc_upper
     end select
   endif

   ! compute surface runoff
   select case(bc_upper)
     case(prescribedHead)
       ! compute surface runoff, which is zero for the prescribed head case
       scalarSurfaceRunoff_IE = 0._rkind ! infiltration excess runoff 
       scalarSurfaceRunoff_SE = 0._rkind ! saturation excess runoff 
       scalarSurfaceRunoff    = 0._rkind ! total surface runoff

     case(liquidFlux)
       ! compute surface runoff (m s-1)
       scalarSurfaceRunoff = scalarRainPlusMelt - scalarSurfaceInfiltration
       scalarSurfaceRunoff_SE = scalarRainPlusMelt * scalarSaturatedArea
       if (scalarRainPlusMelt.gt.xMaxInfilRate) then ! infiltration excess surface runoff occurs
        ! saturation excess surface runoff computed by one of the saturation excess methods, remaining surface runoff is infiltration excess
        scalarSurfaceRunoff_IE = scalarSurfaceRunoff - scalarSurfaceRunoff_SE ! infiltration excess surface runoff     
       else ! infiltration excess runoff does not occur
        scalarSurfaceRunoff_SE = scalarSurfaceRunoff ! saturation excess surface runoff 
        scalarSurfaceRunoff_IE = 0._rkind            ! infiltration excess surface runoff 
       end if         
    end select  

  end associate
 end subroutine update_surfaceFlux

subroutine update_volFracLiq_derivatives
  ! **** updates the derivatives for volumetric fraction of liquid and ice water in each soil layer ****
  ! local variables
  integer(i4b)      :: nLayers         ! number of soil layers to process
  logical(lgt)      :: doIce           ! flag indicating whether ice derivatives are needed

  associate(&
   ! input: model control
   surfRun_SE          => in_surfaceFlux % surfRun_SE         , & ! index defining the saturation excess surface runoff method
   nRoots              => in_surfaceFlux % nRoots             , & ! number of layers that contain roots or take infiltration
   nSoil               => in_surfaceFlux % nSoil              , & ! total number of soil layers
   ! input: state and diagnostic variables
   mLayerTemp          => in_surfaceFlux % mLayerTemp         , & ! temperature (K)
   mLayerMatricHead    => in_surfaceFlux % mLayerMatricHead   , & ! matric head in each soil layer (m)
   ! input: pre-computed derivatives in ...
   dTheta_dTk          => in_surfaceFlux % dTheta_dTk         , & ! ... volumetric liquid water content w.r.t. temperature (K-1)
   dTheta_dPsi         => in_surfaceFlux % dTheta_dPsi        , & ! ... liquid water content w.r.t. liquid water matric potential (m-1)
   ! output: error control
   err      => out_surfaceFlux % err    , & ! error code
   message  => out_surfaceFlux % message  & ! error message
  &)

   ! determine number of layers to process and whether ice derivatives are needed
   if (surfRun_SE ==homegrown_SE) then ! need only root zone derivatives but need ice derivatives
     nLayers = nRoots
     doIce = .true.
   else ! might need entire soil column (FUSE methods), might need ice derivatives (infiltration excess method)
     if (ixInfRateMax_use == noInfiltrationExcess) then
       if (surfRun_SE ==zero_SE) then ! no derivatives needed
         nLayers = 0
         doIce = .false.
       else ! FUSE methods do not need ice derivatives
         nLayers = nSoil
         doIce = .false.
       end if
     else ! infiltration excess method needs ice derivatives in root zone
       if (surfRun_SE == zero_SE) then ! only need root zone derivatives
         nLayers = nRoots
         doIce = .true.
       else ! FUSE methods need soil column derivatives, will compute unused ice derivatives for layers beyond root zone 
         nLayers = nSoil
         doIce = .true.
       end if
     end if ! (if ixInfRateMax_use)
   end if ! (if homegrown_SE)

   if (nLayers > 0) then
     do iLayer=1,nLayers
       Tcrit = crit_soilT( mLayerMatricHead(iLayer) )
       if (mLayerTemp(iLayer) < Tcrit) then
         dVolFracLiq_dWat(iLayer) = 0._rkind
         if(doIce) dVolFracIce_dWat(iLayer) = dTheta_dPsi(iLayer)
       else
         dVolFracLiq_dWat(iLayer) = dTheta_dPsi(iLayer)
         if(doIce) dVolFracIce_dWat(iLayer) = 0._rkind
       end if
     end do
     dVolFracLiq_dTk(:) = dTheta_dTk(:) ! already zeroed out if not below critical temperature
     if(doIce) dVolFracIce_dTk(:) = -dVolFracLiq_dTk(:) ! often can and will simplify one of these terms out
   end if

  end associate
 end subroutine update_volFracLiq_derivatives

 subroutine update_surfaceFlux_FUSE_PRMS_infilArea
  ! **** Update operations for surfaceFlux: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- PRMS ****
  use soil_utils_module,only:LogSumExp  ! smooth max/min
  use soil_utils_module,only:SoftArgMax ! smooth arg max/min (for derivatives of LogSumExp)

  ! local variables
  real(rkind)           :: dS1_dLiq(1:in_surfaceFlux % nSoil)       ! derivative of S1 w.r.t. liquid water content
  real(rkind)           :: S1_T_derivatives(1:2)                    ! array of derivatives for S1_T
  real(rkind)           :: dS1_T_dS1                                ! derivative of S1_T w.r.t S1
  real(rkind)           :: dS1_T_dLiq(1:in_surfaceFlux % nSoil)     ! derivative of S1_T w.r.t liquid water content

  associate(&
   nSoil            => in_surfaceFlux % nSoil,            & ! number of soil layers
   mLayerVolFracLiq => in_surfaceFlux % mLayerVolFracLiq, & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth      => in_surfaceFlux % mLayerDepth,      & ! depth of soil layers (m) 
   theta_sat        => in_surfaceFlux % theta_sat,        & ! soil porosity (-)
   ! output: error control
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)

   ! validation of parameters
   SatArea_max = in_surfaceFlux % FUSE_Ac_max
   phi_tens    = in_surfaceFlux % FUSE_phi_tens
   ! validate input parameters 
   if ((SatArea_max<0._rkind).or.(SatArea_max>1._rkind)) then
    err=10; message=trim(message)//"FUSE PRMS surface runoff error: invalid SatArea_max (max saturated area) value"; return_flag=.true.; return
   end if
   if ((phi_tens<0._rkind).or.(phi_tens>1._rkind)) then
    err=10; message=trim(message)//"FUSE PRMS surface runoff error: invalid phi_tens (tension storage fraction) value"; return_flag=.true.; return
   end if

  ! compute water content in upper FUSE layer
   S1     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in upper FUSE layer (m)
   if (S1 <= 0._rkind) then; io_surfaceFlux % scalarInfilArea = 1._rkind; return; end if ! if no water, unsaturated and all area infiltrates
   S1_max = total_soil_depth * theta_sat                ! max water storage for upper FUSE layer (m)

  ! compute tension water content
   S1_T_max = phi_tens * S1_max
   S1_T     = LogSumExp(-alpha_LSE,[S1,S1_T_max],err) ! smooth approximation to S1_T=min(S1,S1_T_max)
   if(err/=0)then; err=10; message=trim(message)//"FUSE PRMS surface runoff: error in LogSumExp"; return_flag=.true.; return; end if
   if (S1_T < 0._rkind) then ! check for errors
    err=10; message=trim(message)//"FUSE PRMS surface runoff: S1_T is negative (may need to adjust magnitude of alpha_LSE)"; return_flag=.true.; return
   end if

   ! define the infiltrating area and derivatives for the non-frozen part of the cell/basin
   io_surfaceFlux % scalarInfilArea = 1._rkind - (S1_T/S1_T_max)*SatArea_max
   ! define the derivatives
   if(updateInfil)then
     dS1_dLiq           = mLayerDepth(:)
     S1_T_derivatives   = SoftArgMax(-alpha_LSE,[S1,S1_T_max])
     dS1_T_dS1          = S1_T_derivatives(1) 
     dS1_T_dLiq         = dS1_T_dS1 * dS1_dLiq(:)
     dInfilArea_dWat(:) = -(dS1_T_dLiq(:)/S1_T_max)*SatArea_max * dVolFracLiq_dWat(:)     
     dInfilArea_dTk(:)  = -(dS1_T_dLiq(:)/S1_T_max)*SatArea_max * dVolFracLiq_dTk(:)          
   endif ! else derivatives are zero
  end associate

 end subroutine update_surfaceFlux_FUSE_PRMS_infilArea

 subroutine update_surfaceFlux_FUSE_ARNO_VIC_infilArea
  ! **** Update operations for surfaceFlux: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- ARNO/VIC ****
  use soil_utils_module,only:LogSumExp  ! smooth max/min
  use soil_utils_module,only:SoftArgMax ! smooth arg max/min (for derivatives of LogSumExp)

  ! local variables
  real(rkind)            :: dS1_dLiq(1:in_surfaceFlux % nSoil)     ! derivative of S1 w.r.t. liquid water content
  real(rkind)            :: dS1_star_dS1                           ! derivative in S1_star w.r.t S1
  real(rkind)            :: dbase_dS1                              ! derivative of base w.r.t S1
  real(rkind)            :: S1_star_derivatives(1:2)               ! array of derivatives for S1_star from SoftArgMax function

  associate(&
   nSoil            => in_surfaceFlux % nSoil,            & ! number of soil layers
   mLayerVolFracLiq => in_surfaceFlux % mLayerVolFracLiq, & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth      => in_surfaceFlux % mLayerDepth,      & ! depth of soil layers (m) 
   theta_sat        => in_surfaceFlux % theta_sat,        & ! soil porosity (-)
   ! output: error control
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)

   ! validation of input parameters
   b_arnovic = in_surfaceFlux % FUSE_b ! interface ARNO/VIC exponent
   if ((b_arnovic < 0.001_rkind).or.(b_arnovic > 3._rkind)) then
    err=10; message=trim(message)//"FUSE ARNO/VIC exponent must be between 0.001 and 3"; return_flag=.true.; return
   end if

   ! compute water content in FUSE layers
   S1     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in FUSE layers (m)
   if (S1 <= 0._rkind) then; io_surfaceFlux % scalarInfilArea = 1._rkind; return; end if ! if no water, unsaturated and all area infiltrates
   S1_max = total_soil_depth * theta_sat                ! max water storage for FUSE layers (m)

   ! Original FUSE: SatArea = 1 - (1-S1/S1_max)**b_arnovic
   ! Optional: smoothed to prevent negative bases using a smooth approximation of S1_star = min(S1,S1_max)
   !           (Smoothed SatArea) = 1 - (1-S1_star/S1_max)**b_arnovic 
   if (smoother) then ! with smooth approximation of min(S1,S1_max)
    S1_star = LogSumExp(-alpha_LSE,[S1,S1_max],err) ! smooth approximation of min(S1,S1_max) to prevent negative bases
    if(err/=0)then; err=10; message=trim(message)//"FUSE ARNO/VIC surface runoff: error in LogSumExp"; return_flag=.true.; return; end if
   else               ! no smoothing
    S1_star = S1
   end if
   if (S1_star < 0._rkind) then ! check for errors
    err=10; message=trim(message)//&
    &"FUSE ARNO/VIC surface runoff: S1_star is negative (may need to apply smoothing or increase magnitude of alpha_LSE)";return_flag=.true.; return
   end if

   ! compute base value
   base = 1._rkind - S1_star/S1_max

   ! validate base value and add tolerance for round-off error
   if (base < -roundoff_tolerance) then ! if below zero outside of tolerance
    err=10; message=trim(message)//"FUSE ARNO/VIC base value is negative"; return_flag=.true.; return
   else if (base < 0._rkind) then       ! if below zero within tolerance
    base = 0._rkind
   end if

   ! define the infiltrating area and derivatives for the non-frozen part of the cell/basin
   io_surfaceFlux % scalarInfilArea = base**b_arnovic
   
   ! define the derivatives
   if(updateInfil)then
   ! compute derivatives needed for infiltration derivative
     dS1_dLiq = mLayerDepth(:)
   if (smoother) then ! with smooth approximation of min(S1,S1_max)
       S1_star_derivatives = SoftArgMax(-alpha_LSE,[S1,S1_max])
       dS1_star_dS1 = S1_star_derivatives(1)
   else               ! no smoothing
       dS1_star_dS1 = 1._rkind  ! S1_star = S1 if no smoothing
   end if
     dbase_dS1 = -1._rkind/S1_max * dS1_star_dS1
     dInfilArea_dWat(:) = b_arnovic*base**(b_arnovic-1._rkind)*dbase_dS1*dS1_dLiq(:) * dVolFracLiq_dWat(:)     
     dInfilArea_dTk(:)  = b_arnovic*base**(b_arnovic-1._rkind)*dbase_dS1*dS1_dLiq(:) * dVolFracLiq_dTk(:) 
    endif ! else derivatives are zero
  end associate

 end subroutine update_surfaceFlux_FUSE_ARNO_VIC_infilArea


 subroutine update_surfaceFlux_FUSE_TOPMODEL_infilArea
  ! **** Update operations for surfaceFlux: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- TOPMODEL ****
  ! local variables
  complex(rkind)                   :: F1,F2                              ! temporary storage for regularized lower incomplete gamma function values
  real(rkind)                      :: dS1_dLiq(1:in_surfaceFlux % nSoil) ! derivative in S1 w.r.t liquid water content 
  real(rkind)                      :: dzeta_crit_n_dS1                   ! derivative of zeta_crit_n w.r.t S1
  real(rkind)                      :: dzeta_crit_dzeta_crit_n            ! derivative of zeta_crit w.r.t zeta_crit_n
  real(rkind)                      :: dx_crit_dzeta_crit                 ! derivative of x_crit w.r.t zeta_crit
  real(rkind)                      :: dx_crit_dS1                        ! derivative of x_crit w.r.t S1
  real(rkind)                      :: dgammp_dx_crit                     ! derivative of gammp function in SatArea w.r.t x_crit

  associate(&
   nSoil            => in_surfaceFlux % nSoil,            & ! number of soil layers
   mLayerVolFracLiq => in_surfaceFlux % mLayerVolFracLiq, & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth      => in_surfaceFlux % mLayerDepth,      & ! depth of soil layers (m) 
   theta_sat        => in_surfaceFlux % theta_sat,        & ! soil porosity (-)
   ! output: error control
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)

   ! interface FUSE input parameters
   lambda       = in_surfaceFlux % FUSE_lambda
   chi_topmodel = in_surfaceFlux % FUSE_chi
   mu           = in_surfaceFlux % FUSE_mu
   n_topmodel   = in_surfaceFlux % FUSE_n

   ! compute water content in lower FUSE layer, here the entire soil column is used
   S1     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in lower FUSE layer (m)
   if (S1 <= 0._rkind) then; io_surfaceFlux % scalarInfilArea = 1._rkind; return; end if ! if no water, unsaturated and all area infiltrates
   S1_max = total_soil_depth * theta_sat                ! max water storage for lower FUSE layer (m)

   ! validate of parameters
   if ((lambda < 5._rkind ).or.(lambda > 10._rkind)) then
    err=10; message=trim(message)//"FUSE TOPMODEL lambda value must be between 5 and 10"; return_flag=.true.; return
   end if
   if (lambda <= mu) then
    err=10; message=trim(message)//"FUSE TOPMODEL lambda value must be greater than mu value"; return_flag=.true.; return
   end if
   if ((chi_topmodel < 2._rkind ).or.(chi_topmodel > 5._rkind)) then
    err=10; message=trim(message)//"FUSE TOPMODEL chi_topmodel value must be between 2 and 5"; return_flag=.true.; return
   end if
   if ((mu < 2.5_rkind ).or.(mu > 3.5_rkind)) then
    err=10; message=trim(message)//"FUSE TOPMODEL mu value must be between 2.5 and 3.5"; return_flag=.true.; return
   end if
   if ((n_topmodel < 3.5_rkind).or.(n_topmodel > 10._rkind)) then ! validate TOPMODEL exponent to avoid divergence of lambda_n
    err=10; message=trim(message)//"FUSE TOPMODEL exponent must be between 3.5 and 10"; return_flag=.true.; return
   end if
   ! validate water content values, these should be guaranteed by earlier checks but just in case
   if (S1 < 0._rkind) then; err=10; message=trim(message)//"negative water content value detected in lower FUSE layer"; return_flag=.true.; return; end if
   if (S1 > S1_max) then; err=10; message=trim(message)//"water content in lower FUSE layer exceeds max storage"; return_flag=.true.; return; end if

  ! check water content in lower FUSE layer 
  if (S1 > 0._rkind) then ! if some water is stored in lower FUSE layer
   ! set FUSE parameters - input parameters are lambda, chi_topmodel, and mu
   alpha_topmodel=(lambda-mu)/chi_topmodel

   ! * compute the mean power-transformed topographic index *
   ! compute regularized lower incomplete Gamma function values
   F1=gammp_complex(alpha_topmodel,(-(mu*n_topmodel - mu*chi_topmodel - (n_topmodel - chi_topmodel)*zeta_upper)/n_topmodel)/chi_topmodel)
   F2=gammp_complex(alpha_topmodel,(-(mu*n_topmodel - mu*chi_topmodel)/n_topmodel)/chi_topmodel)

   ! mean power-transformed topographic index (translated to Fortran from SageMath)
   lambda_n=(cmplx(-mu + zeta_upper,0._rkind,rkind)**alpha_topmodel*(F1 - 1)*exp(mu/n_topmodel)*gamma(alpha_topmodel)/cmplx(-(mu*n_topmodel - mu*chi_topmodel - &
           &(n_topmodel - chi_topmodel)*zeta_upper)/(n_topmodel*chi_topmodel),0._rkind,rkind)**alpha_topmodel - cmplx(-mu,0._rkind,rkind)**alpha_topmodel*(F2 - 1)*exp(mu/n_topmodel)*gamma(alpha_topmodel)/&
           &cmplx(-(mu*n_topmodel - mu*chi_topmodel)/(n_topmodel*chi_topmodel),0._rkind,rkind)**alpha_topmodel)/(cmplx(chi_topmodel,0._rkind,rkind)**alpha_topmodel*gamma(alpha_topmodel))

   ! compute critical zeta value
   ! note: to obtain physical topography values, only the real part of lambda_n is used 
   zeta_crit_n=lambda_n%re*S1_max/S1 ! power-transformed critical topographic index
   if (zeta_crit_n <= 0._rkind) then; err=10; message=trim(message)//"FUSE TOPMODEL zeta_crit_n <= 0"; return_flag=.true.; return; end if

   zeta_crit=n_topmodel*log(zeta_crit_n) ! critical topographic index in log space

   ! transform to x random variable and validate result
   x_crit=zeta_crit-mu
   if (x_crit < -roundoff_tolerance) then ! less than zero outside tolerance
     err=10; message=trim(message)//"FUSE TOPMODEL zeta_crit must be greater or equal to mu, try increasing lambda or decreasing mu";return_flag=.true.; return
   else if (x_crit < 0._rkind) then       ! less than zero but within tolerance
    x_crit = 0._rkind
   end if

   ! define the infiltrating area and derivatives for the non-frozen part of the cell/basin
   io_surfaceFlux % scalarInfilArea = gammp(alpha_topmodel,x_crit/chi_topmodel)

  else ! if (S1 == 0) no water is stored in lower FUSE layer (based on asymptotic behaviour of integral in eq. 9c of Clark et al. (2008))
   io_surfaceFlux % scalarInfilArea = 1._rkind
  end if

   ! define the derivatives
   if(updateInfil)then
     dS1_dLiq = mLayerDepth(:)    
     dzeta_crit_n_dS1 = -lambda_n%re*S1_max/S1**2_i4b  
     dzeta_crit_dzeta_crit_n = ( n_topmodel*zeta_crit_n**(n_topmodel-1._rkind) ) / zeta_crit_n**n_topmodel
     dx_crit_dzeta_crit = 1._rkind
     dx_crit_dS1 = dx_crit_dzeta_crit * dzeta_crit_dzeta_crit_n * dzeta_crit_n_dS1
     dgammp_dx_crit = ( (x_crit/chi_topmodel)**(alpha_topmodel-1._rkind) * exp(-x_crit/chi_topmodel) )/chi_topmodel/gamma(alpha_topmodel)
     dInfilArea_dWat(:) = dgammp_dx_crit * dx_crit_dS1 * dS1_dLiq(:) * dVolFracLiq_dWat(:)     
     dInfilArea_dTk(:)  = dgammp_dx_crit * dx_crit_dS1 * dS1_dLiq(:) * dVolFracLiq_dTk(:)
   endif ! else derivatives are zero
  end associate

 end subroutine update_surfaceFlux_FUSE_TOPMODEL_infilArea

 subroutine update_surfaceFlux_prescribedHead
  ! **** Update operations for surfaceFlux: prescribed pressure head condition ****
  associate(&
   ! input: state and diagnostic variables
   scalarMatricHeadLiq => in_surfaceFlux % scalarMatricHeadLiq , & ! liquid matric head in the upper-most soil layer (m)
   ! input: depth of each soil layer (m)
   mLayerDepth  => in_surfaceFlux % mLayerDepth  , & ! depth of each soil layer (m)
   ! input: diriclet boundary conditions
   upperBoundHead   => in_surfaceFlux % upperBoundHead  , & ! upper boundary condition for matric head (m)
   ! input: transmittance
   surfaceSatHydCond => in_surfaceFlux % surfaceSatHydCond , & ! saturated hydraulic conductivity at the surface (m s-1)
   dHydCond_dTemp    => in_surfaceFlux % dHydCond_dTemp    , & ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   iceImpedeFac      => in_surfaceFlux % iceImpedeFac      , & ! ice impedence factor in the upper-most soil layer (-)
   ! input: soil parameters
   vGn_alpha           => in_surfaceFlux % vGn_alpha           , & ! van Genuchten "alpha" parameter (m-1)
   vGn_n               => in_surfaceFlux % vGn_n               , & ! van Genuchten "n" parameter (-)
   vGn_m               => in_surfaceFlux % vGn_m               , & ! van Genuchten "m" parameter (-)
   ! input-output: hydraulic conductivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   surfaceHydCond => io_surfaceFlux % surfaceHydCond , & ! hydraulic conductivity (m s-1)
   ! output: infiltration
   scalarSurfaceInfiltration => io_surfaceFlux % scalarSurfaceInfiltration  , & ! surface infiltration (m s-1)
   ! output: derivatives in surface infiltration w.r.t. ...
   scalarSoilControl  => io_surfaceFlux % scalarSoilControl    , & ! soil control on infiltration for derivative
   dq_dHydStateVec    => out_surfaceFlux % dq_dHydStateVec     , & ! ... hydrology state in every soil layer (m s-1 or s-1)
   dq_dNrgStateVec    => out_surfaceFlux % dq_dNrgStateVec     , & ! ... energy state in every soil layer (m s-1 K-1)
   ! output: error control
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)

   ! compute transmission and the capillary flux
   surfaceHydCond = hydCond_psi(upperBoundHead,surfaceSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
   cflux = -surfaceHydCond*(scalarMatricHeadLiq - upperBoundHead) / (mLayerDepth(1)*0.5_rkind)

   ! compute the total flux (no glacier melt infiltration assumed for prescribed head condition)
   scalarSurfaceInfiltration = cflux + surfaceHydCond
   scalarSoilControl = 0._rkind

   ! compute the derivatives at the surface, only has a non-zero value for the upper-most soil layer
   if(updateInfil)then
     dq_dHydStateVec(1) = -surfaceHydCond/(mLayerDepth(1)/2._rkind)
     ! note: energy state variable is temperature (transformed outside soilLiqFlux_module if needed)
     dq_dNrgStateVec(1) = -(dHydCond_dTemp/2._rkind)*(scalarMatricHeadLiq - upperBoundHead)/(mLayerDepth(1)*0.5_rkind) + dHydCond_dTemp/2._rkind
   end if
   ! bottom melt infiltration is not considered for prescribed head condition, so derivatives are zero

   ! * additional assignment statements for surfaceFlux input-output object based on presribed head values *
   ! the infiltration is always constrained by the prescribed head so the maximum infiltration rate is set to missing
   io_surfaceFlux % xMaxInfilRate    = realMissing ! maximum infiltration rate (m s-1)
   ! no soil ice assumed for prescribed head condition
   io_surfaceFlux % scalarFrozenArea = 0._rkind      ! fraction of area that is considered impermeable due to soil ice (-)
   ! all area is available for infiltration, and to complement this saturated area (i.e., part where saturation excess runoff occurs) is set to zero
   io_surfaceFlux % scalarInfilArea     = 1._rkind ! fraction of area where water can infiltrate, may be frozen (-)
   io_surfaceFlux % scalarSaturatedArea = 0._rkind ! fraction of area that is considered saturated (-)

  end associate
 end subroutine update_surfaceFlux_prescribedHead

 subroutine update_surfaceFlux_homegrown_infilArea
  ! **** Update operations for surfaceFlux: homegrown saturation excess runoff condition ****
  call update_surfaceFlux_liquidFlux_computation_root_layers 
  call update_surfaceFlux_liquidFlux_computation_available_capacity; if (return_flag) return 
  call update_surfaceFlux_liquidFlux_computation_homegrown  ! this calculates infiltration area ignoring if frozen or not, depends on available capacity (depends on ice and root zone)
 end subroutine update_surfaceFlux_homegrown_infilArea

 subroutine update_surfaceFlux_liquidFlux_noinfratemax
  ! **** Update operations for surfaceFlux: no infiltration excess****
  associate(&
   ! input: model control
   surfRun_SE => in_surfaceFlux % surfRun_SE         & ! index defining the saturation excess surface runoff method
  &)
   io_surfaceFlux % xMaxInfilRate = veryBig ! set to a very large number so rainPlusMelt never exceeds this
   if (surfRun_SE /= homegrown_SE) then  ! frozen area (depends on ice and root zone)
    call update_surfaceFlux_liquidFlux_computation_root_layers
   end if
  end associate
  ! -- main computations - these always need to run
  call update_surfaceFlux_liquidFlux_computation_frozen_area
 end subroutine update_surfaceFlux_liquidFlux_noinfratemax

 subroutine update_surfaceFlux_liquidFlux_calculate_infratemax
  ! **** Update operations for surfaceFlux: infiltration excess possible - calculate max infiltration rate ****
  associate(&
   ! input: model control
   surfRun_SE => in_surfaceFlux % surfRun_SE         & ! index defining the saturation excess surface runoff method
  &)
   if (surfRun_SE /= homegrown_SE) then  ! infiltration rate max depends on available capacity (depends on ice and root zone) and frozen area (depends on ice and root zone)
     call update_surfaceFlux_liquidFlux_computation_root_layers 
     call update_surfaceFlux_liquidFlux_computation_available_capacity; if (return_flag) return
   end if
  end associate
  ! -- main computations - these always need to run
  call update_surfaceFlux_liquidFlux_computation_frozen_area
  call update_surfaceFlux_liquidFlux_computation_max_infiltration_rate
 end subroutine update_surfaceFlux_liquidFlux_calculate_infratemax

 subroutine update_surfaceFlux_liquidFlux_computation_root_layers 
  ! **** Update operations for surfaceFlux: root layer water computation ****
  associate(&
   ! input: model control
   nRoots              => in_surfaceFlux % nRoots            , & ! number of soil layers with roots (-)
   ! input: state and diagnostic variables
   mLayerVolFracLiq    => in_surfaceFlux % mLayerVolFracLiq  , & ! volumetric liquid water content in each soil layer (-)
   mLayerVolFracIce    => in_surfaceFlux % mLayerVolFracIce  , & ! volumetric ice content in each soil layer (-)
   ! input: depth of soil layers (m)
   mLayerDepth         => in_surfaceFlux % mLayerDepth      , & ! depth of each soil layer (m)
   iLayerHeight        => in_surfaceFlux % iLayerHeight     , & ! height at the interface of each layer for soil layers only (m)
   rootingDepth        => in_surfaceFlux % rootingDepth       & ! rooting depth (m)
  &)
 
   ! define the storage in the root zone (m) and derivatives, first initialize
   rootZoneLiq = 0._rkind
   rootZoneIce = 0._rkind
   dRootZoneLiq_dWat(:) = 0._rkind
   dRootZoneIce_dWat(:) = 0._rkind
   dRootZoneLiq_dTk(:)  = 0._rkind
   dRootZoneIce_dTk(:)  = 0._rkind

   ! process layers where the roots or infiltration extend to the bottom of the layer
   if (nRoots > 1) then
     do iLayer=1,nRoots-1
       rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(iLayer)*mLayerDepth(iLayer)
       rootZoneIce = rootZoneIce + mLayerVolFracIce(iLayer)*mLayerDepth(iLayer)
       if(updateInfil)then
         dRootZoneLiq_dWat(iLayer) = dVolFracLiq_dWat(iLayer)*mLayerDepth(iLayer)
         dRootZoneIce_dWat(iLayer) = dVolFracIce_dWat(iLayer)*mLayerDepth(iLayer)
         dRootZoneLiq_dTk(iLayer)  = dVolFracLiq_dTk(iLayer) *mLayerDepth(iLayer)
         dRootZoneIce_dTk(iLayer)  = dVolFracIce_dTk(iLayer) *mLayerDepth(iLayer)
       end if
     end do
   end if
   ! process layers where the roots or infiltration end in the current layer
   rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(nRoots)*min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
   rootZoneIce = rootZoneIce + mLayerVolFracIce(nRoots)*min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
   if(updateInfil)then
     dRootZoneLiq_dWat(nRoots) = dVolFracLiq_dWat(nRoots)*min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneIce_dWat(nRoots) = dVolFracIce_dWat(nRoots)*min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneLiq_dTk(nRoots)  = dVolFracLiq_dTk(nRoots)* min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneIce_dTk(nRoots)  = dVolFracIce_dTk(nRoots)* min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
   endif

  end associate
 end subroutine update_surfaceFlux_liquidFlux_computation_root_layers

 subroutine update_surfaceFlux_liquidFlux_computation_available_capacity 
  ! **** Update operations for surfaceFlux: compute and check available capacity to hold water ****
  associate(&
   ! input: soil parameters
   theta_sat           => in_surfaceFlux % theta_sat   , & ! soil porosity (-)
   rootingDepth        => in_surfaceFlux % rootingDepth, & ! rooting depth (m)
   ! output: error control
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)
   availCapacity = theta_sat*rootingDepth - rootZoneIce
   if (rootZoneLiq > availCapacity+verySmaller) then
     err=20; message=trim(message)//'liquid water in the root zone exceeds capacity'; return_flag=.true.; return
   end if

  end associate
 end subroutine update_surfaceFlux_liquidFlux_computation_available_capacity

 subroutine update_surfaceFlux_liquidFlux_computation_max_infiltration_rate
  ! **** Update operations for surfaceFlux: max infiltration rate and derivatives ****
  associate(&
   ! input: transmittance
   surfaceSatHydCond => in_surfaceFlux % surfaceSatHydCond , & ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: soil parameters
   zScale_TOPMODEL     => in_surfaceFlux % zScale_TOPMODEL     , & ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
   f_hydCond           => in_surfaceFlux % f_hydCond           , & ! decay rate of hydraulic conductivity with depth, exponential profile (m-1)
   compactedDepth      => in_surfaceFlux % compactedDepth      , & ! depth where k_soil reaches the compacted value, power-law profile (m)
   rootingDepth        => in_surfaceFlux % rootingDepth        , & ! rooting depth (m)
   wettingFrontSuction => in_surfaceFlux % wettingFrontSuction , & ! Green-Ampt wetting front suction (m)
   mLayerDepth         => in_surfaceFlux % mLayerDepth         , & ! depth of each soil layer (m)
   ! input-output: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate    => io_surfaceFlux % xMaxInfilRate , & ! maximum infiltration rate (m s-1)
   ! output: error control
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)
   ! define the depth to the wetting front (m) and derivatives
   depthWettingFront = (rootZoneLiq/availCapacity)*min(rootingDepth,total_soil_depth)
   if(updateInfil)then
     dDepthWettingFront_dWat(:)=( dRootZoneLiq_dWat(:)*min(rootingDepth, total_soil_depth) + dRootZoneIce_dWat(:)*depthWettingFront )/availCapacity
     dDepthWettingFront_dTk(:) =( dRootZoneLiq_dTk(:) *min(rootingDepth, total_soil_depth) + dRootZoneIce_dTk(:)*depthWettingFront  )/availCapacity
   end if

   ! process hydraulic conductivity-controlled infiltration rate
   select case(ixInfRateMax_use)  ! maximum infiltration rate parameterization (noInfExcess set in update_surfaceFlux)
    case(topmodel_GA)
     ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
     ! NOTE: this re-derives the conductivity profile rather than reading iLayerSatHydCond, so it must use the same shape as satHydCond
     select case(ix_hc_profile_use)
       case(expLaw_profile)   ! K(z) = K_0*exp(-f*z), decays with depth but stays finite at the base of the soil
         hydCondWettingFront = surfaceSatHydCond * exp(-f_hydCond*depthWettingFront)
         dHydCondWF_dDepth   = -f_hydCond*hydCondWettingFront
       ! power law, the same profile satHydCond builds: K(z) = K_0*(1 - z/R)**(zScale_TOPMODEL-1) for
       ! R = refDepth, decaying over the whole column and reaching zero at the base. It is not floored
       ! below compactedDepth -- see the NOTE in satHydCond.
       case(powerLaw_profile)
         refDepth_use = total_soil_depth
         if (total_soil_depth < compactedDepth) refDepth_use = compactedDepth + 1._rkind ! as in satHydCond
         ! clamp the wetting front off the base, where (1 - z/R)**(n-2) in the derivative is singular
         depthWF_use = min(depthWettingFront, 0.99_rkind*refDepth_use)
         hydCondWettingFront = surfaceSatHydCond * ( (1._rkind - depthWF_use/refDepth_use)**(zScale_TOPMODEL - 1._rkind) )
         dHydCondWF_dDepth   = -surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) &
                                * ( (1._rkind - depthWF_use/refDepth_use)**(zScale_TOPMODEL - 2._rkind) )/refDepth_use
         if (depthWettingFront > depthWF_use) dHydCondWF_dDepth = 0._rkind
       case(constant)         ! K uniform with depth, so the wetting front sees the surface conductivity
         hydCondWettingFront = surfaceSatHydCond
         dHydCondWF_dDepth   = 0._rkind
       case default; err=20; message=trim(message)//'unknown hydraulic conductivity profile for the topmodel_GA infiltration rate'; return_flag=.true.; return
     end select
     ! define the maximum infiltration rate (m s-1)
     xMaxInfilRate = hydCondWettingFront*( (wettingFrontSuction + depthWettingFront)/depthWettingFront )  ! maximum infiltration rate (m s-1)
     ! define the derivatives
     if(updateInfil)then
       fPart1    = hydCondWettingFront
       fPart2    = (wettingFrontSuction + depthWettingFront)/depthWettingFront
       dPart1(:) = dHydCondWF_dDepth*dDepthWettingFront_dWat(:)
       dPart2(:) = -dDepthWettingFront_dWat(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
       dxMaxInfilRate_dWat(:) = fPart1*dPart2(:) + fPart2*dPart1(:)
       dPart1(:) = dHydCondWF_dDepth*dDepthWettingFront_dTk(:)
       dPart2(:) = -dDepthWettingFront_dTk(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
       dxMaxInfilRate_dTk(:)  = fPart1*dPart2(:) + fPart2*dPart1(:)
     endif
    case(GreenAmpt)
      ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
      hydCondWettingFront = surfaceSatHydCond ! Green-Ampt assumes homogeneous soil, therefore the whole soil column has the same hydraulic conductivity
      ! define the maximum infiltration rate (m s-1)
      xMaxInfilRate = hydCondWettingFront * (1._rkind + (1._rkind - depthWettingFront/total_soil_depth) * wettingFrontSuction/depthWettingFront) ! Ks * (1 + (Md) * S/F)
      ! define the derivatives
      if(updateInfil)then
        dxMaxInfilRate_dWat(:) = -hydCondWettingFront*wettingFrontSuction*dDepthWettingFront_dWat(:)/depthWettingFront**2_i4b
        dxMaxInfilRate_dTk(:)  = -hydCondWettingFront*wettingFrontSuction*dDepthWettingFront_dTk(:)/depthWettingFront**2_i4b
      endif
   end select
  end associate
 end subroutine update_surfaceFlux_liquidFlux_computation_max_infiltration_rate

 subroutine update_surfaceFlux_liquidFlux_computation_homegrown
  ! **** Update operations for surfaceFlux: infiltrating area (ignoring frozen area) for homegrown saturation excess condition ****
  associate(&
   ! input: model control
   mLayerVolFracLiq => in_surfaceFlux % mLayerVolFracLiq    , & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth      => in_surfaceFlux % mLayerDepth         , & ! depth of each soil layer (m)
   ! input: soil parameters
   theta_sat        => in_surfaceFlux % theta_sat           , & ! soil porosity (-)
   qSurfScale       => in_surfaceFlux % qSurfScale          , & ! scaling factor in the surface runoff parameterization (-)
   ! input-output: surface runoff and infiltration flux (m s-1)
   scalarInfilArea  => io_surfaceFlux % scalarInfilArea       & ! fraction of area where water can infiltrate, may be frozen (-)
  &)
   ! define the infiltrating area and derivatives for the ignoring if frozen or not
   if (qSurfScale < qSurfScaleMax) then
     fracCap         = rootZoneLiq/(maxFracCap*availCapacity)                              ! fraction of available root zone filled with water
     fInfRaw         = 1._rkind - exp(-qSurfScale*(1._rkind - fracCap))                          ! infiltrating area -- allowed to violate solution constraints
     scalarInfilArea = min(0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor)), 1._rkind)   ! infiltrating area -- constrained
     ! define the derivatives
     if(updateInfil)then
       if (0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor))< 1._rkind) then
         dfracCap(:) = ( dRootZoneLiq_dWat(:)/maxFracCap + dRootZoneIce_dWat(:)*fracCap )/availCapacity
         dfInfRaw(:) = -qSurfScale*dfracCap(:) * exp(-qSurfScale*(1._rkind - fracCap))
         dInfilArea_dWat(:) = 0.5_rkind*dfInfRaw(:) * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
         dfracCap(:) = ( dRootZoneLiq_dTk(:)/maxFracCap + dRootZoneIce_dTk(:)*fracCap )/availCapacity
         dfInfRaw(:) = -qSurfScale*dfracCap(:) * exp(-qSurfScale*(1._rkind - fracCap))
         dInfilArea_dTk(:)  = 0.5_rkind*dfInfRaw(:) * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
       endif ! else derivatives are zero
     endif
   else
     scalarInfilArea = 1._rkind ! derivatives are zero
   end if
  end associate
 end subroutine update_surfaceFlux_liquidFlux_computation_homegrown

 subroutine update_surfaceFlux_liquidFlux_computation_frozen_area
  ! **** Update operations for surfaceFlux: get impermeable area due to soil freezing ****
  associate(&
   ! input: soil parameters
   soilIceScale        => in_surfaceFlux % soilIceScale       , & ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   soilIceCV           => in_surfaceFlux % soilIceCV          , & ! soil ice CV in Gamma distribution used to define frozen area (-)
   ! output: frozen area
   scalarFrozenArea    => io_surfaceFlux % scalarFrozenArea     & ! fraction of area that is considered impermeable due to soil ice (-)
  &)
   ! define the impermeable area and derivatives due to frozen ground
   if (rootZoneIce > tiny(rootZoneIce)) then  ! (avoid divide by zero)
      alpha = 1._rkind/(soilIceCV**2_i4b)     ! shape parameter in the Gamma distribution
      xLimg = alpha*soilIceScale/rootZoneIce  ! upper limit of the integral
     !if we use this, we will have a derivative of scalarFrozenArea w.r.t. water and temperature in each layer (through mLayerVolFracIce)
     ! Should fix to deal with frozen area in the root zone, calculations may be expensive
     !scalarFrozenArea = 1._rkind - gammp(alpha,xLimg)      ! fraction of frozen area
     !if(updateInfil)then
     !  dFrozenArea_dWat(:) = -dgammp_dx(alpha,xLimg)*(-alpha*soilIceScale/rootZoneIce**2_i4b)*dRootZoneIce_dWat(:)
     !  dFrozenArea_dTk(:)  = -dgammp_dx(alpha,xLimg)*(-alpha*soilIceScale/rootZoneIce**2_i4b)*dRootZoneIce_dTk(:)
     !end if
     scalarFrozenArea = 0._rkind
   else
     scalarFrozenArea = 0._rkind
   end if
  end associate
 end subroutine update_surfaceFlux_liquidFlux_computation_frozen_area

 subroutine update_surfaceFlux_liquidFlux_infiltration
  ! **** Update operations for surfaceFlux: final infiltration calculations ****
  ! local variables
  real(rkind) :: scalarInfilArea_unfrozen ! infiltration area that is not frozen
  real(rkind) :: rootZoneDepth            ! depth of active root-zone layers used in saturation limiter (m)
  real(rkind) :: compHeadRootZone         ! mean positive pressure head over active layers (m)
  real(rkind) :: compLimiter              ! smooth limiter that closes infiltration area under strong compression (-)
  real(rkind) :: compLimiter_prev         ! pre-limiter infiltration area used for chain rule (-)
  real(rkind) :: compArg                  ! argument of logistic compression limiter (-)
  real(rkind) :: dCompLimiter_dHead       ! derivative of compression limiter w.r.t. mean positive head (m-1)
  real(rkind),parameter :: compHeadCutoff=1._rkind   ! positive head where compression closure begins (m)
  real(rkind),parameter :: compHeadWidth =0.02_rkind ! smoothing width for compression closure (m)
  real(rkind),parameter :: headSmooth =1.e-4_rkind   ! smoothing for max(psi,0) approximation (m)
  real(rkind) :: posHead(1:in_surfaceFlux % nSoil)   ! smooth positive part of matric head (m)
  real(rkind) :: dPosHead_dPsi(1:in_surfaceFlux % nSoil) ! derivative of smooth positive head w.r.t. matric head (-)
  real(rkind) :: dCompHead_dWat(1:in_surfaceFlux % nSoil) ! derivative of mean positive head w.r.t. hyd state (m)
  integer(i4b) :: ixTop, ixBot             ! top and bottom layer indices for the active root zone

  ! compute infiltration
  associate(&
   ! input: model control and state
   ix_groundwatr     => in_surfaceFlux % ix_groundwatr,         & ! index defining the groundwater parameterization
   surfRun_SE        => in_surfaceFlux % surfRun_SE,            & ! index defining the saturation excess surface runoff method
   bc_lower          => in_surfaceFlux % bc_lower,              & ! index defining the lower boundary condition
   nSoil             => in_surfaceFlux % nSoil,                 & ! number of soil layers
   nGlce             => in_surfaceFlux % nGlce,                 & ! number of glacier debris layers
   nRoots            => in_surfaceFlux % nRoots,                & ! number of layers that contain roots or take infiltration (-)
   ixIce             => in_surfaceFlux % ixIce,                 & ! index of lowest ice layer
   mLayerMatricHead  => in_surfaceFlux % mLayerMatricHead,      & ! matric head in each soil layer (m)
   mLayerVolFracLiq  => in_surfaceFlux % mLayerVolFracLiq,      & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth       => in_surfaceFlux % mLayerDepth,           & ! depth of each soil layer (m)
   ! input: soil parameters
   theta_sat         => in_surfaceFlux % theta_sat,             & ! soil porosity (-)
   ! input: flux at the upper boundary
   scalarRainPlusMelt  => in_surfaceFlux % scalarRainPlusMelt,  & ! rain plus melt plus lake drainage, used as input to the soil zone before computing surface runoff (m s-1)
   ! input-output: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate       => io_surfaceFlux % xMaxInfilRate,       & ! maximum infiltration rate (m s-1)
   scalarSoilControl   => io_surfaceFlux % scalarSoilControl,   & ! soil control on infiltration for derivative
   scalarFrozenArea    => io_surfaceFlux % scalarFrozenArea,    & ! fraction of area that is considered impermeable due to soil ice (-)
   scalarInfilArea     => io_surfaceFlux % scalarInfilArea,     & ! fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea => io_surfaceFlux % scalarSaturatedArea, & ! saturated area fraction (-)
   scalarSurfaceInfiltration => io_surfaceFlux % scalarSurfaceInfiltration, & ! surface infiltration (m s-1)
   ! output: derivatives in surface infiltration w.r.t. ...
   dq_dHydStateVec => out_surfaceFlux % dq_dHydStateVec, & ! ... hydrology state in every soil layer (m s-1 or s-1)
   dq_dNrgStateVec => out_surfaceFlux % dq_dNrgStateVec, & ! ... energy state in every soil layer (m s-1 K-1)
   ! output: error control
   err     => out_surfaceFlux % err,     & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)

  ! check infiltration area and define saturated area
   if (scalarInfilArea < 0._rkind) then; err=20; message=trim(message)//'infiltration area less than zero'; return_flag=.true.; return; end if

   ! Saturation-area limiters
   if(surfRun_SE==homegrown_SE)then ! infiltration area based on all layers
     ixTop = ixIce + 1
     ixBot = nRoots
   else ! if not homegrown_SE, infiltration area based on all layers
     ixTop = 1
     ixBot = nSoil
   end if

  ! Close infiltration under saturation for blocked lower boundaries
  rootZoneDepth = sum(mLayerDepth(ixTop:ixBot))
  if (ixTop <= ixBot .and. (bc_lower/=freeDrainage .or. nGlce>0)) then ! glacier always has lower boundary zero flux (blocked boundary)
    ! drives infiltration area to zero once positive pressure becomes large
    posHead(:) = 0._rkind
    dPosHead_dPsi(:) = 0._rkind
    posHead(ixTop:ixBot) = 0.5_rkind*(mLayerMatricHead(ixTop:ixBot) + sqrt(mLayerMatricHead(ixTop:ixBot)**2_i4b + headSmooth**2_i4b)) ! smooth positive part of matric head (m)
    dPosHead_dPsi(ixTop:ixBot) = 0.5_rkind*(1._rkind + mLayerMatricHead(ixTop:ixBot)/sqrt(mLayerMatricHead(ixTop:ixBot)**2_i4b + headSmooth**2_i4b))

    ! compute derivatives of mean positive head w.r.t. water state variables
    dCompHead_dWat(:) = 0._rkind
    compHeadRootZone = sum(posHead(ixTop:ixBot)*mLayerDepth(ixTop:ixBot))/rootZoneDepth
    compArg = (compHeadRootZone - compHeadCutoff)/compHeadWidth
    compLimiter = 1._rkind/(1._rkind + exp(2._rkind*compArg))
    dCompLimiter_dHead = -(2._rkind/compHeadWidth)*compLimiter*(1._rkind - compLimiter)
    dCompHead_dWat(ixTop:ixBot) = (mLayerDepth(ixTop:ixBot)/rootZoneDepth) * dPosHead_dPsi(ixTop:ixBot)

    ! apply compression limiter to infiltration area and compute derivatives
    compLimiter_prev = scalarInfilArea
    scalarInfilArea = scalarInfilArea*compLimiter
    if(updateInfil)then
      dInfilArea_dWat(:) = dInfilArea_dWat(:)*compLimiter + compLimiter_prev*dCompLimiter_dHead*dCompHead_dWat(:)
      dInfilArea_dTk(:)  = dInfilArea_dTk(:)*compLimiter
    end if
   end if
   scalarSaturatedArea = 1._rkind - scalarInfilArea

   ! unfrozen infiltration area and infiltration (m s-1)
   scalarInfilArea_unfrozen=(1._rkind - scalarFrozenArea)*scalarInfilArea
   scalarSoilControl = 0._rkind
   scalarSurfaceInfiltration = scalarInfilArea_unfrozen * min(scalarRainPlusMelt,xMaxInfilRate)

   ! Compute total runoff derivatives, do w.r.t. infiltration only, scalarRainPlusMelt accounted for in computJacob* module
   if(updateInfil)then
     if (xMaxInfilRate > scalarRainPlusMelt) then
       scalarSoilControl = scalarInfilArea_unfrozen  ! derivative dependent on scalarRainPlusMelt (needed to compute scalarRainPlusMelt derivative inside computJacob*)
     elseif (xMaxInfilRate < scalarRainPlusMelt) then ! dInfilRate_d dependent on layers not at surface
       dInfilRate_dWat(:) = dxMaxInfilRate_dWat(:)
       dInfilRate_dTk(:)  = dxMaxInfilRate_dTk(:)
     end if
     ! Do not need to break into IE and SE components since they are never used separately in the Jacobian assembly
     dq_dHydStateVec(:) = (1._rkind - scalarFrozenArea)&
                         * ( dInfilArea_dWat(:)*min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dWat(:) )&
                         + (-dFrozenArea_dWat(:))*scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
     ! energy state variable is temperature (transformed outside soilLiqFlux_module if needed)
     dq_dNrgStateVec(:) = (1._rkind - scalarFrozenArea)&
                         * ( dInfilArea_dTk(:) *min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dTk(:)  )&
                         + (-dFrozenArea_dTk(:)) *scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
   end if
  end associate

  ! set surface hydraulic conductivity to missing (not used for flux condition)
  associate(&
   ! input-output: hydraulic conductivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   surfaceHydCond => io_surfaceFlux % surfaceHydCond   & ! hydraulic conductivity (m s-1)
  &)
   surfaceHydCond = realMissing
  end associate

 end subroutine update_surfaceFlux_liquidFlux_infiltration

 subroutine finalize_surfaceFlux
  ! **** Finalize operations for surfaceFlux ****
  ! final error check
  associate(&
   err     => out_surfaceFlux % err    , & ! error code
   message => out_surfaceFlux % message  & ! error message
  &)
   if(err/=0)then; message=trim(message)//'unanticipated error in surfaceFlux subroutine'; return_flag=.true.; return; end if
  end associate
 end subroutine finalize_surfaceFlux

end subroutine surfaceFlux

! ***************************************************************************************************************
! private subroutine iLayerFlux: compute the fluxes and derivatives at layer interfaces
! ***************************************************************************************************************
subroutine iLayerFlux(in_iLayerFlux,out_iLayerFlux)
  ! ---------------------------------------------------------------------------------------------------------------------------
  ! input: model control, state variables, coordinate variables, temperature derivatives, transmittance variables
  type(in_type_iLayerFlux),intent(in)   :: in_iLayerFlux   ! class object for input data
  ! output: transmittance variables and vertical flux at layer interface, derivatives, and error control
  type(out_type_iLayerFlux),intent(out) :: out_iLayerFlux  ! class object for output data
  ! ---------------------------------------------------------------------------------------------------------------------------
  ! local variables (named variables to provide index of 2-element vectors)
  integer(i4b),parameter           :: ixUpper=1            ! index of upper node in the 2-element vectors
  integer(i4b),parameter           :: ixLower=2            ! index of lower node in the 2-element vectors
  logical(lgt),parameter           :: useGeometric=.false. ! switch between the arithmetic and geometric mean
  ! local variables (Darcy flux)
  real(rkind)                      :: dPsi                 ! spatial difference in matric head (m)
  real(rkind)                      :: dz                   ! spatial difference in layer mid-points (m)
  real(rkind)                      :: cflux                ! capillary flux (m s-1)
  ! error control
  logical(lgt)                     :: return_flag          ! flag for return statements
  ! ---------------------------------------------------------------------------------------------------------------------------

  call initialize_iLayerFlux

  call update_iLayerFlux;   if (return_flag) return

  call finalize_iLayerFlux; if (return_flag) return

contains

 subroutine initialize_iLayerFlux
  ! **** Initialize operations for iLayerFlux ****
  return_flag=.false. ! initialize return flag
  associate(&
   err     => out_iLayerFlux % err    , & ! error code
   message => out_iLayerFlux % message  & ! error message
  &)
   ! initialize error control
   err=0; message="iLayerFlux/" ! initialize error control
  end associate
 end subroutine initialize_iLayerFlux

 subroutine update_iLayerFlux
  ! **** Update operations for iLayerFlux ****

  ! ** compute the fluxes
  call update_iLayerFlux_fluxes; if (return_flag) return

  ! ** compute the derivatives
  call update_iLayerFlux_derivatives; if (return_flag) return

 end subroutine update_iLayerFlux

 subroutine update_iLayerFlux_fluxes
  ! **** Update operations for iLayerFlux: compute fluxes ****
  associate(&
   ! input: state variables
   nodeMatricHeadLiqTrial => in_iLayerFlux % nodeMatricHeadLiqTrial, & ! liquid matric head at the soil nodes (m)
   ! input: model coordinate variables
   nodeHeight => in_iLayerFlux % nodeHeight, & ! height at the mid-point of the lower layer (m)
   ! input: transmittance
   nodeHydCondTrial => in_iLayerFlux % nodeHydCondTrial, & ! hydraulic conductivity at layer mid-points (m s-1)
   ! output: tranmsmittance at the layer interface (scalars)
   iLayerHydCond => out_iLayerFlux % iLayerHydCond, & ! hydraulic conductivity at the interface between layers (m s-1)
   ! output: vertical flux at the layer interface (scalars)
   iLayerLiqFluxSoil => out_iLayerFlux % iLayerLiqFluxSoil, & ! vertical flux of liquid water at the layer interface (m s-1)
   ! output: error control
   err     => out_iLayerFlux % err    , & ! error code
   message => out_iLayerFlux % message  & ! error message
  &)

   ! compute the vertical flux of liquid water
   ! compute the hydraulic conductivity at the interface
   if (useGeometric) then
     iLayerHydCond   = sqrt(nodeHydCondTrial(ixLower)   * nodeHydCondTrial(ixUpper))
   else
     iLayerHydCond   = (nodeHydCondTrial(ixLower)   + nodeHydCondTrial(ixUpper))*0.5_rkind
   end if

   dz = nodeHeight(ixLower) - nodeHeight(ixUpper)
   dPsi          = nodeMatricHeadLiqTrial(ixLower) - nodeMatricHeadLiqTrial(ixUpper)
   cflux         = -iLayerHydCond * dPsi/dz
   ! compute the total flux (add gravity flux, positive downwards)
   iLayerLiqFluxSoil = cflux + iLayerHydCond

  end associate
 end subroutine update_iLayerFlux_fluxes

 subroutine update_iLayerFlux_derivatives
  ! **** Update operations for iLayerFlux: compute derivatives ****
  ! * local variables (derivative in Darcy's flux) *
  ! derivatives at the layer interface
  real(rkind) :: dHydCondIface_dMatricAbove  ! hydraulic conductivity w.r.t. matric head in layer above
  real(rkind) :: dHydCondIface_dMatricBelow  ! hydraulic conductivity w.r.t. matric head in layer below
  associate(&
   ! input: temperature derivatives
   dPsiLiq_dTemp   => in_iLayerFlux % dPsiLiq_dTemp , & ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   dHydCond_dTemp  => in_iLayerFlux % dHydCond_dTemp, & ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: transmittance
   nodeHydCondTrial => in_iLayerFlux % nodeHydCondTrial, & ! hydraulic conductivity at layer mid-points (m s-1)
   ! input: transmittance derivatives
   dHydCond_dMatric => in_iLayerFlux % dHydCond_dMatric, & ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
   ! output: tranmsmittance at the layer interface (scalars)
   iLayerHydCond => out_iLayerFlux % iLayerHydCond, & ! hydraulic conductivity at the interface between layers (m s-1)
   ! output: derivatives in fluxes w.r.t. ...
   dq_dHydStateAbove => out_iLayerFlux % dq_dHydStateAbove, & ! ... matric head or volumetric liquid water in the layer above (m s-1 or s-1)
   dq_dHydStateBelow => out_iLayerFlux % dq_dHydStateBelow, & ! ... matric head or volumetric liquid water in the layer below (m s-1 or s-1)
   ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
   dq_dNrgStateAbove => out_iLayerFlux % dq_dNrgStateAbove, & ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
   dq_dNrgStateBelow => out_iLayerFlux % dq_dNrgStateBelow, & ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
   ! output: error control
   err     => out_iLayerFlux % err    , & ! error code
   message => out_iLayerFlux % message  & ! error message
  &)

   ! derivatives in hydraulic conductivity
   if (useGeometric) then
     dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_rkind/max(iLayerHydCond,verySmaller)
     dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_rkind/max(iLayerHydCond,verySmaller)
   else
     dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)/2._rkind
     dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)/2._rkind
   end if
   ! derivatives in the flux w.r.t. matric head
   dq_dHydStateAbove = -dHydCondIface_dMatricAbove*dPsi/dz + iLayerHydCond/dz + dHydCondIface_dMatricAbove
   dq_dHydStateBelow = -dHydCondIface_dMatricBelow*dPsi/dz - iLayerHydCond/dz + dHydCondIface_dMatricBelow
   ! derivative in the flux w.r.t. temperature
   dq_dNrgStateAbove = -(dHydCond_dTemp(ixUpper)/2._rkind)*dPsi/dz + iLayerHydCond*dPsiLiq_dTemp(ixUpper)/dz + dHydCond_dTemp(ixUpper)/2._rkind
   dq_dNrgStateBelow = -(dHydCond_dTemp(ixLower)/2._rkind)*dPsi/dz - iLayerHydCond*dPsiLiq_dTemp(ixLower)/dz + dHydCond_dTemp(ixLower)/2._rkind

  end associate
 end subroutine update_iLayerFlux_derivatives

 subroutine finalize_iLayerFlux
  ! **** Finalize operations for iLayerFlux ****
  associate(&
   err     => out_iLayerFlux % err    , & ! error code
   message => out_iLayerFlux % message  & ! error message
  &)
   ! final error check
   if(err/=0)then; message=trim(message)//'unanticipated error in iLayerFlux'; return_flag=.true.; return; end if
  end associate
 end subroutine finalize_iLayerFlux

end subroutine iLayerFlux

! ***************************************************************************************************************
! private subroutine qDrainFlux: compute the drainage flux from the bottom of the soil profile and its derivative
! ***************************************************************************************************************
subroutine qDrainFlux(in_qDrainFlux,io_qDrainFlux,out_qDrainFlux)
  USE soil_utils_module,only:hydCond_psi ! compute hydraulic conductivity as a function of matric head (m s-1)
  ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
  implicit none
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! input: model control, variables, boundary conditions, transmittance variables, and soil parameters
  type(in_type_qDrainFlux) ,intent(in)  :: in_qDrainFlux      ! object for qDrainFlux input data
  ! input-output: soil control for derivatives
  type(io_type_qDrainFlux),intent(inout):: io_qDrainFlux      ! object for qDrainFlux input-output data
  ! output: hydraulic conductivity and diffusivity, drainage fluxes and derivatives, and error control
  type(out_type_qDrainFlux),intent(out) :: out_qDrainFlux     ! object for qDrainFlux output data
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! local variables
   ! local variables (Darcy flux)
  real(rkind)                      :: zWater                  ! effective water table depth (m)
  real(rkind)                      :: dPsi                    ! spatial difference in matric head (m)
  real(rkind)                      :: dz                      ! spatial difference in layer mid-points (m)
  real(rkind)                      :: cflux                   ! capillary flux (m s-1)
  integer(i4b)                     :: bc_lower_use            ! mutable copy of lower boundary-condition index
  ! error control
  logical(lgt)                     :: return_flag             ! flag for return statements
  ! -----------------------------------------------------------------------------------------------------------------------------

   call initialize_qDrainFlux

   call update_qDrainFlux;   if (return_flag) return

   call finalize_qDrainFlux; if (return_flag) return

contains

 subroutine initialize_qDrainFlux
  ! ** Initialize operations for qDrainFlux **
  return_flag=.false. ! initialize return flag
  associate(&
   ! output: error control
   err     => out_qDrainFlux % err    , & ! error code
   message => out_qDrainFlux % message  & ! error message
  &)
   ! initialize error control
   err=0; message="qDrainFlux/"
  end associate
 end subroutine initialize_qDrainFlux

 subroutine update_qDrainFlux
  ! ** Update operations for qDrainFlux **
  associate(&
   ! input: model control
   bc_lower      => in_qDrainFlux % bc_lower, & ! index defining the type of boundary conditions
   nGlce         => in_qDrainFlux % nGlce,    & ! number of glacier ice layers
   ! output: error control
   err     => out_qDrainFlux % err    , &       ! error code
   message => out_qDrainFlux % message  &       ! error message
  &)
   bc_lower_use = bc_lower
   if(nGlce>0) bc_lower_use = zeroFlux ! if glacier debris, nothing can drain into impermeable glacier ice layer

   ! determine lower boundary condition
   select case(bc_lower_use)
     case(prescribedHead) ! specified matric head value
       call update_qDrainFlux_prescribedHead; if (return_flag) return
     case(funcBottomHead) ! specified matric head function
       call update_qDrainFlux_funcBottomHead; if (return_flag) return
     case(freeDrainage)   ! free drainage
       call update_qDrainFlux_freeDrainage;   if (return_flag) return
     case(zeroFlux)       ! zero flux from soil but allow for capillary flux upwards glacier ice melt
       if(nGlce>0)then
         call update_qDrainFlux_capillaryFlux; if (return_flag) return
       else
         call update_qDrainFlux_zeroFlux;      if (return_flag) return
       end if
     case default; err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return_flag=.true.; return
   end select 

  end associate
 end subroutine update_qDrainFlux

 subroutine update_qDrainFlux_prescribedHead
  ! ** Update operations for qDrainFlux: prescribed pressure head value at bottom boundary **
  associate(&
   ! input: state and diagnostic variables
   nodeMatricHeadLiq => in_qDrainFlux % nodeMatricHeadLiq, &  ! liquid matric head in the lowest unsaturated node (m)
   ! input: model coordinate variables
   nodeDepth  => in_qDrainFlux % nodeDepth , &                ! depth of the lowest unsaturated soil layer (m)
   ! input: diriclet boundary conditions
   lowerBoundHead  => in_qDrainFlux % lowerBoundHead , &      ! lower boundary condition for matric head (m)
   ! input: transmittance
   bottomSatHydCond  => in_qDrainFlux % bottomSatHydCond , &  ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   iceImpedeFac      => in_qDrainFlux % iceImpedeFac     , &  ! ice impedence factor in the upper-most soil layer (-)
   ! input: transmittance derivatives
   dHydCond_dTemp   => in_qDrainFlux % dHydCond_dTemp  , &    ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   vGn_alpha       => in_qDrainFlux % vGn_alpha      , &      ! van Genuchten "alpha" parameter (m-1)
   vGn_n           => in_qDrainFlux % vGn_n          , &      ! van Genuchten "n" parameter (-)
   vGn_m           => in_qDrainFlux % vGn_m          , &      ! van Genuchten "m" parameter (-)
   ! output: hydraulic conductivity at the bottom of the unsaturated zone
   bottomHydCond => out_qDrainFlux % bottomHydCond, &         ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   ! output: drainage flux from the bottom of the soil profile
   scalarDrainage => out_qDrainFlux % scalarDrainage, &       ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   dq_dHydStateUnsat => out_qDrainFlux % dq_dHydStateUnsat, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   dq_dNrgStateUnsat => out_qDrainFlux % dq_dNrgStateUnsat, & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
   ! output: error control
   err     => out_qDrainFlux % err    , &                     ! error code
   message => out_qDrainFlux % message  &                     ! error message
  &)

   ! compute flux
   bottomHydCond = hydCond_psi(lowerBoundHead,bottomSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
   cflux = -bottomHydCond*(lowerBoundHead  - nodeMatricHeadLiq) / (nodeDepth*0.5_rkind)
   scalarDrainage = cflux + bottomHydCond

   ! hydrology derivatives
   dq_dHydStateUnsat = bottomHydCond/(nodeDepth/2._rkind)
   ! energy derivatives
   dq_dNrgStateUnsat = -(dHydCond_dTemp/2._rkind)*(lowerBoundHead  - nodeMatricHeadLiq)/(nodeDepth*0.5_rkind)&
                     & + dHydCond_dTemp/2._rkind

  end associate
 end subroutine update_qDrainFlux_prescribedHead

 subroutine update_qDrainFlux_funcBottomHead
  ! ** Update operations for qDrainFlux: prescribed pressure head function at bottom boundary **
  associate(&
   ! input: state and diagnostic variables
   nodeMatricHeadLiq => in_qDrainFlux % nodeMatricHeadLiq, &  ! liquid matric head in the lowest unsaturated node (m)
   ! input: model coordinate variables
   nodeHeight => in_qDrainFlux % nodeHeight, &                ! height of the lowest unsaturated soil node (m)
   ! input: derivative in soil water characteristic
   node_dPsiLiq_dTemp  => in_qDrainFlux % node_dPsiLiq_dTemp , &  ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: transmittance
   surfaceSatHydCond => in_qDrainFlux % surfaceSatHydCond, &  ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: soil parameters
   kAnisotropic    => in_qDrainFlux % kAnisotropic   , &      ! anisotropy factor for lateral hydraulic conductivity (-)
   zScale_TOPMODEL => in_qDrainFlux % zScale_TOPMODEL, &      ! scale factor for TOPMODEL-ish baseflow parameterization (m)
   ! output: drainage flux from the bottom of the soil profile
   scalarDrainage => out_qDrainFlux % scalarDrainage, &       ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   dq_dHydStateUnsat => out_qDrainFlux % dq_dHydStateUnsat, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   dq_dNrgStateUnsat => out_qDrainFlux % dq_dNrgStateUnsat, & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
   ! output: error control
   err     => out_qDrainFlux % err    , &                     ! error code
   message => out_qDrainFlux % message  &                     ! error message
  &)

   ! compute flux
   zWater = nodeHeight - nodeMatricHeadLiq
   scalarDrainage = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)

   ! hydrology derivatives
   dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
   ! energy derivatives
   dq_dNrgStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)*node_dPsiLiq_dTemp/zScale_TOPMODEL

  end associate
 end subroutine update_qDrainFlux_funcBottomHead

 subroutine update_qDrainFlux_freeDrainage
  ! ** Update operations for qDrainFlux: free drainage at bottom boundary **
  associate(&
   ! input: transmittance
   nodeHydCond       => in_qDrainFlux % nodeHydCond    , &    ! hydraulic conductivity at the node itself (m s-1)
   ! input: transmittance derivatives
   dHydCond_dMatric => in_qDrainFlux % dHydCond_dMatric, &    ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
   dHydCond_dTemp   => in_qDrainFlux % dHydCond_dTemp  , &    ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   kAnisotropic    => in_qDrainFlux % kAnisotropic  , &       ! anisotropy factor for lateral hydraulic conductivity (-)
   ! output: drainage flux from the bottom of the soil profile
   scalarDrainage => out_qDrainFlux % scalarDrainage, &       ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   dq_dHydStateUnsat => out_qDrainFlux % dq_dHydStateUnsat, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   dq_dNrgStateUnsat => out_qDrainFlux % dq_dNrgStateUnsat, & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
   ! output: error control
   err     => out_qDrainFlux % err    , &                     ! error code
   message => out_qDrainFlux % message  &                     ! error message
  &)

   scalarDrainage = nodeHydCond*kAnisotropic ! compute flux

   ! hydrology derivatives
   dq_dHydStateUnsat = dHydCond_dMatric*kAnisotropic
   ! energy derivatives
   dq_dNrgStateUnsat = dHydCond_dTemp*kAnisotropic

  end associate
 end subroutine update_qDrainFlux_freeDrainage

 subroutine update_qDrainFlux_zeroFlux
  ! ** Update operations for qDrainFlux: zero flux condition at bottom boundary **
  associate(&
   ! output: drainage flux from the bottom of the soil profile
   scalarDrainage => out_qDrainFlux % scalarDrainage, &       ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   dq_dHydStateUnsat => out_qDrainFlux % dq_dHydStateUnsat, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   dq_dNrgStateUnsat => out_qDrainFlux % dq_dNrgStateUnsat  & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
  &)

   scalarDrainage = 0._rkind
   dq_dHydStateUnsat = 0._rkind
   dq_dNrgStateUnsat = 0._rkind

  end associate
 end subroutine update_qDrainFlux_zeroFlux


 subroutine update_qDrainFlux_capillaryFlux
  ! ** Compute the capillary flux at the bottom boundary **
  associate(&
   ! input: state and diagnostic variables
   scalarGlceMelt    => in_qDrainFlux % scalarGlceMelt,       &  ! glacier melt at the bottom of the soil zone (m s-1)
   nodeMatricHeadLiq => in_qDrainFlux % nodeMatricHeadLiq,    &  ! liquid matric head in the lowest soil layer (m)
   nodeDepth         => in_qDrainFlux % nodeDepth,            &  ! depth of the lowest soil layer (m)
   ! input: transmittance
   bottomSatHydCond  => in_qDrainFlux % bottomSatHydCond,     &  ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   ! input: transmittance derivatives
   dHydCond_dMatric => in_qDrainFlux % dHydCond_dMatric,      &  ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
   dHydCond_dTemp   => in_qDrainFlux % dHydCond_dTemp,        &  ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: derivatives in soil water characteristic
   node_dPsiLiq_dTemp  => in_qDrainFlux % node_dPsiLiq_dTemp, &  ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input-output: derivatives
   scalarSoilControlBot => io_qDrainFlux % scalarSoilControlBot, & ! soil control on bottom capillary fluxes for derivative
   ! output: drainage flux from the bottom of the soil profile
   scalarDrainage => out_qDrainFlux % scalarDrainage,         &  ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   dq_dHydStateUnsat => out_qDrainFlux % dq_dHydStateUnsat,   & ! ... state variable in lowest soil layer (m s-1 or s-1)
   dq_dNrgStateUnsat => out_qDrainFlux % dq_dNrgStateUnsat,   & ! ... energy state variable in lowest soil layer (m s-1 K-1)
   ! output: error control
   err     => out_qDrainFlux % err,     &  ! error code
   message => out_qDrainFlux % message  &  ! error message
  &)

   dz = nodeDepth*0.5_rkind
   ! compute the capillary flux
   dPsi  = -nodeMatricHeadLiq ! if not saturated, then matric head is negative
   cflux = -bottomSatHydCond * dPsi/dz
   scalarDrainage = cflux + bottomSatHydCond ! compute the total flux (add gravity flux, positive downwards)

   ! derivatives
   dq_dHydStateUnsat = -bottomSatHydCond/dz + dHydCond_dMatric
   dq_dNrgStateUnsat = -bottomSatHydCond*node_dPsiLiq_dTemp/dz + dHydCond_dTemp

   ! Constrain flux: can't drain into glacier and can't exceed melt supply
   scalarSoilControlBot = 0._rkind ! need for derivatives
   if (scalarDrainage < scalarGlceMelt .or. scalarDrainage > 0._rkind) then
     dq_dHydStateUnsat = 0._rkind
     dq_dNrgStateUnsat = 0._rkind
     if (scalarDrainage < scalarGlceMelt) scalarSoilControlBot = 1._rkind
   endif
   scalarDrainage = max(scalarGlceMelt, min(0._rkind, scalarDrainage))

  end associate
 end subroutine update_qDrainFlux_capillaryFlux


 subroutine finalize_qDrainFlux
  ! ** Finalize operations for qDrainFlux **
  associate(&
   ! output: error control
   err     => out_qDrainFlux % err    , & ! error code
   message => out_qDrainFlux % message  & ! error message
  &)
   ! final error check
   if(err/=0)then; message=trim(message)//'unanticipated error in qDrainFlux'; return_flag=.true.; return; end if
  end associate
 end subroutine finalize_qDrainFlux

end subroutine qDrainFlux

end module soilLiqFlux_module
