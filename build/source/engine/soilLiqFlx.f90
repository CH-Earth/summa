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

module soilLiqFlx_module
! -----------------------------------------------------------------------------------------------------------

! data types
USE nrtype
USE data_types,only:var_d                  ! x%var(:)       (rkind)
USE data_types,only:var_ilength            ! x%var(:)%dat   (i4b)
USE data_types,only:var_dlength            ! x%var(:)%dat   (rkind)
USE data_types,only:in_type_soilLiqFlx     ! derived type for intent(in) arguments
USE data_types,only:io_type_soilLiqFlx     ! derived type for intent(inout) arguments
USE data_types,only:out_type_soilLiqFlx    ! derived type for intent(out) arguments
USE data_types,only:in_type_diagv_node     ! derived type for intent(in) arguments 
USE data_types,only:out_type_diagv_node    ! derived type for intent(out) arguments 
USE data_types,only:in_type_surfaceFlx     ! derived type for intent(in) arguments
USE data_types,only:io_type_surfaceFlx     ! derived type for intent(inout) arguments
USE data_types,only:out_type_surfaceFlx    ! derived type for intent(out) arguments
USE data_types,only:in_type_iLayerFlux     ! derived type for intent(in) arguments 
USE data_types,only:out_type_iLayerFlux    ! derived type for intent(out) arguments 
USE data_types,only:in_type_qDrainFlux     ! derived type for intent(in) arguments 
USE data_types,only:out_type_qDrainFlux    ! derived type for intent(out) arguments 

! missing values
USE globalData,only:integerMissing         ! missing integer
USE globalData,only:realMissing            ! missing real number

! physical constants
USE multiconst,only:iden_water             ! intrinsic density of water    (kg m-3)

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
  ! look-up values for method used to compute derivative
  numerical,                  & ! numerical solution
  analytical,                 & ! analytical solution
  ! look-up values for the form of Richards' equation
  moisture,                   & ! moisture-based form of Richards' equation
  mixdform,                   & ! mixed form of Richards' equation
  ! look-up values for the type of hydraulic conductivity profile
  constant,                   & ! constant hydraulic conductivity with depth
  powerLaw_profile,           & ! power-law profile
  ! look-up values for the choice of groundwater parameterization
  qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                  & ! a big bucket (lumped aquifer model)
  noExplicit,                 & ! no explicit groundwater parameterization
  ! look-up values for the choice of boundary conditions for hydrology
  prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  funcBottomHead,             & ! function of matric head in the lower-most layer
  freeDrainage,               & ! free drainage
  liquidFlux,                 & ! liquid water flux
  zeroFlux                      ! zero flux

! -----------------------------------------------------------------------------------------------------------
implicit none
private
public :: soilLiqFlx
! constant parameters
real(rkind),parameter     :: verySmall=1.e-12_rkind       ! a very small number (used to avoid divide by zero)
real(rkind),parameter     :: dx=1.e-8_rkind               ! finite difference increment
contains
! ***************************************************************************************************************
! public subroutine soilLiqFlx: compute liquid water fluxes and their derivatives
! ***************************************************************************************************************
subroutine soilLiqFlx(&
                      ! input: model control, trial state variables, derivatives, and fluxes
                      in_soilLiqFlx,                & ! intent(in): model control, trial state variables, derivatives, and fluxes
                      ! input-output: data structures
                      mpar_data,                    & ! intent(in):    model parameters
                      indx_data,                    & ! intent(in):    model indices
                      prog_data,                    & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,                    & ! intent(inout): model fluxes for a local HRU
                      ! input-output: diagnostic variables, fluxes, and derivatives
                      io_soilLiqFlx,                & ! intent(inout): diagnostic variables, fluxes, and derivatives
                      ! output: error control
                      out_soilLiqFlx)                 ! intent(out): error control
  ! utility modules
  USE soil_utils_module,only:volFracLiq               ! compute volumetric fraction of liquid water
  USE soil_utils_module,only:matricHead               ! compute matric head (m)
  USE soil_utils_module,only:dTheta_dPsi              ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
  USE soil_utils_module,only:dPsi_dTheta              ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  USE soil_utils_module,only:hydCond_psi              ! compute hydraulic conductivity as a function of matric head
  USE soil_utils_module,only:hydCond_liq              ! compute hydraulic conductivity as a function of volumetric liquid water content
  USE soil_utils_module,only:hydCondMP_liq            ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control, trial state variables, derivatives, and fluxes
  type(in_type_soilLiqFlx),intent(in)    :: in_soilLiqFlx              ! model control, trial state variables, derivatives, and fluxes
  ! input-output: data structures
  type(var_dlength),intent(in)           :: mpar_data                  ! model parameters
  type(var_ilength),intent(in)           :: indx_data                  ! state vector geometry
  type(var_dlength),intent(in)           :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)        :: diag_data                  ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)        :: flux_data                  ! model fluxes for a local HRU
  ! input-output: diagnostic variables, fluxes, and derivatives
  type(io_type_soilLiqFlx),intent(inout) :: io_soilLiqFlx              ! diagnostic variables, fluxes, and derivatives
  ! output: error control
  type(out_type_soilLiqFlx),intent(out)  :: out_soilLiqFlx             ! error code and error message
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables: general
  character(LEN=256)                  :: cmessage                      ! error message of downwind routine
  integer(i4b)                        :: nSoil                         ! number of soil layers
  integer(i4b)                        :: ibeg,iend                     ! start and end indices of the soil layers in concatanated snow-soil vector
  integer(i4b)                        :: iLayer,iSoil                  ! index of soil layer
  integer(i4b)                        :: ixLayerDesired(1)             ! layer desired (scalar solution)
  integer(i4b)                        :: ixTop                         ! top layer in subroutine call
  integer(i4b)                        :: ixBot                         ! bottom layer in subroutine call
  ! transpiration sink term
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: mLayerTranspireFrac     ! fraction of transpiration allocated to each soil layer (-)
  ! diagnostic variables
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: iceImpedeFac            ! ice impedence factor at layer mid-points (-)
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: mLayerDiffuse           ! diffusivity at layer mid-point (m2 s-1)
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: dHydCond_dVolLiq        ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: dDiffuse_dVolLiq        ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: dHydCond_dTemp          ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  real(rkind),dimension(0:in_soilLiqFlx % nSoil)  :: iLayerHydCond           ! hydraulic conductivity at layer interface (m s-1)
  real(rkind),dimension(0:in_soilLiqFlx % nSoil)  :: iLayerDiffuse           ! diffusivity at layer interface (m2 s-1)
  ! compute surface flux
  integer(i4b)                                    :: nRoots                  ! number of soil layers with roots
  integer(i4b)                                    :: ixIce                   ! index of the lowest soil layer that contains ice
  real(rkind),dimension(0:in_soilLiqFlx % nSoil)  :: iLayerHeight            ! height of the layer interfaces (m)
  ! error control
  logical(lgt)                                    :: return_flag             ! flag for return statements
  ! -------------------------------------------------------------------------------------------------------------------------------------------------

  ! ** Initialize indices, error control, and get layer information ** 
  call initialize_soilLiqFlx; if (return_flag) return 

  ! ** Compute transpiration, diagnostic variables, infiltration, and interface fluxes **
  call update_soilLiqFlx;     if (return_flag) return

  ! ** Final error control **
  call finalize_soilLiqFlx;   if (return_flag) return
 
contains

 subroutine initialize_soilLiqFlx
  ! **** Initial operations for soilLiqFlx module subroutine ****

  ! ** assign variables used in main associate block **
  nSoil = in_soilLiqFlx % nSoil ! get number of soil layers from input arguments

  ! get indices for the data structures
  ibeg = indx_data%var(iLookINDEX%nSnow)%dat(1) + 1
  iend = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nSoil)%dat(1)

  ! get a copy of iLayerHeight
  ! NOTE: performance hit, though cannot define the shape (0:) with the associate construct
  iLayerHeight(0:nSoil) = prog_data%var(iLookPROG%iLayerHeight)%dat(ibeg-1:iend)  ! height of the layer interfaces (m)

  ! ** initialize error control **
  return_flag=.false.
  associate(&
    err                   => out_soilLiqFlx % err,                  & ! intent(out): error code
    message               => out_soilLiqFlx % cmessage              & ! intent(out): error message
  &)
   err=0; message='soilLiqFlx/' ! initialize error control
  end associate

  ! ** get the indices for the soil layers **
  associate(&
   scalarSolution => in_soilLiqFlx % scalarSolution,             & ! intent(in): flag to denote if implementing the scalar solution
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

  ! ** identify the number of layers that contain roots **
  associate(&
   rootingDepth => mpar_data%var(iLookPARAM%rootingDepth)%dat(1),& ! intent(in): rooting depth (m)
   err          => out_soilLiqFlx % err,                         & ! intent(out): error code
   message      => out_soilLiqFlx % cmessage                     & ! intent(out): error message
  &) 
   nRoots = count(iLayerHeight(0:nSoil-1) < rootingDepth-verySmall)
   if (nRoots==0) then
     message=trim(message)//'no layers with roots'
     err=20; return_flag=.true.; return
   end if
  end associate

  ! ** identify lowest soil layer with ice **
  ! NOTE: cannot use count because there may be an unfrozen wedge
  associate(&
    mLayerVolFracIceTrial => in_soilLiqFlx % mLayerVolFracIceTrial & ! intent(in): volumetric fraction of ice at the current iteration (-)
  &)
   ixIce = 0  ! initialize the index of the ice layer (0 means no ice in the soil profile)
   do iLayer=1,nSoil ! (loop through soil layers)
     if (mLayerVolFracIceTrial(iLayer) > verySmall) ixIce = iLayer
   end do
  end associate
 end subroutine initialize_soilLiqFlx

 subroutine update_soilLiqFlx
  ! **** Main computations for soilLiqFlx module subroutine ****

  if ( .not. (in_soilLiqFlx % scalarSolution .and. ixTop>1) ) then ! check the need to compute transpiration
   call compute_transpiration_sink; if (return_flag) return
  end if  

  call compute_diagnostic_variables; if (return_flag) return

  call compute_surface_infiltration; if (return_flag) return

  call compute_interface_fluxes_derivatives; if (return_flag) return

  if ( .not. (in_soilLiqFlx % scalarSolution .and. ixTop<nSoil) ) then ! define the need to compute drainage
   call compute_drainage_flux; if (return_flag) return
  end if
 end subroutine update_soilLiqFlx

 subroutine finalize_soilLiqFlx
  ! **** Final operations for soilLiqFlx module subroutine ****
 
  ! final error control check for robustness
  associate(&
   err          => out_soilLiqFlx % err,                         & ! intent(out): error code
   message      => out_soilLiqFlx % cmessage                     & ! intent(out): error message
  &)
   if (err/=0) then; message=trim(message)//trim("finalize_soilLiqFlx: final error check failed"); return_flag=.true.; return; end if
  end associate
 end subroutine finalize_soilLiqFlx

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
   err          => out_soilLiqFlx % err,     & ! intent(out): error code
   message      => out_soilLiqFlx % cmessage & ! intent(out): error message
  &)
   ! check fractions sum to one
   if (abs(sum(mLayerTranspireFrac) - 1._rkind) > verySmall) then
    message=trim(message)//'fraction transpiration in soil layers does not sum to one'
    err=20; return_flag=.true.; return
   end if
  end associate
 end subroutine finalize_transpiration_loss_fraction

 subroutine update_transpiration_loss
  ! **** Update transpiration loss from each soil layer (kg m-2 s-1 --> m s-1)*****
  associate(&
   scalarCanopyTranspiration => in_soilLiqFlx % scalarCanopyTranspiration, & ! canopy transpiration (kg m-2 s-1)
   mLayerTranspire           => io_soilLiqFlx % mLayerTranspire,   & ! transpiration loss from each soil layer (m s-1)
   ! intent(inout): derivatives in the soil layer transpiration flux ...
   mLayerdTrans_dCanWat  => io_soilLiqFlx % mLayerdTrans_dCanWat,  & ! ... w.r.t. canopy total water
   mLayerdTrans_dTCanair => io_soilLiqFlx % mLayerdTrans_dTCanair, & ! ... w.r.t. canopy air temperature
   mLayerdTrans_dTCanopy => io_soilLiqFlx % mLayerdTrans_dTCanopy, & ! ... w.r.t. canopy temperature
   mLayerdTrans_dTGround => io_soilLiqFlx % mLayerdTrans_dTGround, & ! ... w.r.t. ground temperature
   ! intent(in): derivative in canopy transpiration ...
   dCanopyTrans_dCanWat  => in_soilLiqFlx % dCanopyTrans_dCanWat,  & ! ... w.r.t. canopy total water content (s-1)
   dCanopyTrans_dTCanair => in_soilLiqFlx % dCanopyTrans_dTCanair, & ! ... w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTCanopy => in_soilLiqFlx % dCanopyTrans_dTCanopy, & ! ... w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTGround => in_soilLiqFlx % dCanopyTrans_dTGround, & ! ... w.r.t. ground temperature (kg m-2 s-1 K-1)
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
  call in_diagv_node % initialize(iSoil,in_soilLiqFlx,model_decisions,diag_data,mpar_data,flux_data)
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
   err          => out_soilLiqFlx % err,     & ! error code
   message      => out_soilLiqFlx % cmessage & ! error message
  &)
   call out_diagv_node % finalize(iSoil,nSoil,io_soilLiqFlx,mLayerDiffuse,iceImpedeFac,&
                                  &dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dTemp,err,cmessage)
   if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end associate
 end subroutine finalize_compute_diagnostic_variables

 subroutine compute_surface_infiltration
  ! **** compute infiltration at the surface and its derivative w.r.t. mass in the upper soil layer ****
  type(in_type_surfaceFlx)  ::  in_surfaceFlx
  type(io_type_surfaceFlx)  ::  io_surfaceFlx
  type(out_type_surfaceFlx) :: out_surfaceFlx

  call initialize_compute_surface_infiltration(in_surfaceFlx,io_surfaceFlx)

  call update_compute_surface_infiltration(in_surfaceFlx,io_surfaceFlx,out_surfaceFlx)

  call finalize_compute_surface_infiltration(io_surfaceFlx,out_surfaceFlx); if (return_flag) return

 end subroutine compute_surface_infiltration

 subroutine initialize_compute_surface_infiltration(in_surfaceFlx,io_surfaceFlx)
  ! **** Initialize operations for compute_surface_infiltration ****
  type(in_type_surfaceFlx),intent(out) :: in_surfaceFlx
  type(io_type_surfaceFlx),intent(out) :: io_surfaceFlx
  ! set derivative w.r.t. state above to zero (does not exist)
  associate(&
   ! intent(inout): flux derivatives ... 
   dq_dHydStateAbove => io_soilLiqFlx % dq_dHydStateAbove,& ! ... in layer interfaces w.r.t. state variables in the layer above
   dq_dNrgStateAbove => io_soilLiqFlx % dq_dNrgStateAbove & ! ... w.r.t. temperature in the layer above (m s-1 K-1)
  &)
   dq_dHydStateAbove(0) = 0._rkind
   dq_dNrgStateAbove(0) = 0._rkind
  end associate

  ! compute surface flux and its derivative...
  call in_surfaceFlx % initialize(nRoots,ixIce,nSoil,ibeg,iend,in_soilLiqFlx,io_soilLiqFlx,&
                                 &model_decisions,prog_data,mpar_data,flux_data,diag_data,&
                                 &iLayerHeight,dHydCond_dTemp,iceImpedeFac)
  call io_surfaceFlx % initialize(nSoil,io_soilLiqFlx,iLayerHydCond,iLayerDiffuse)
 end subroutine initialize_compute_surface_infiltration

 subroutine update_compute_surface_infiltration(in_surfaceFlx,io_surfaceFlx,out_surfaceFlx)
  ! **** Update operations for compute_surface_infiltration ****
  type(in_type_surfaceFlx) ,intent(in)    ::  in_surfaceFlx
  type(io_type_surfaceFlx) ,intent(inout) ::  io_surfaceFlx
  type(out_type_surfaceFlx),intent(out)   :: out_surfaceFlx
  call surfaceFlx(io_soilLiqFlx,in_surfaceFlx,io_surfaceFlx,out_surfaceFlx)
 end subroutine update_compute_surface_infiltration

 subroutine finalize_compute_surface_infiltration(io_surfaceFlx,out_surfaceFlx)
  ! **** Finalize operations for compute_surface_infiltration ****
  type(io_type_surfaceFlx) ,intent(in) :: io_surfaceFlx
  type(out_type_surfaceFlx),intent(in) :: out_surfaceFlx

  ! interface object data components with local name space
  call  io_surfaceFlx % finalize(nSoil,io_soilLiqFlx,iLayerHydCond,iLayerDiffuse)
  associate(&
   err     => out_soilLiqFlx % err,     & ! error code
   message => out_soilLiqFlx % cmessage & ! error message
  &)
   call out_surfaceFlx % finalize(io_soilLiqFlx,err,cmessage)
   if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end associate

  ! include base soil evaporation as the upper boundary flux
  associate(&
   iLayerLiqFluxSoil         => io_soilLiqFlx % iLayerLiqFluxSoil,      & ! liquid flux at soil layer interfaces (m s-1)
   scalarGroundEvaporation   => in_soilLiqFlx % scalarGroundEvaporation,& ! ground evaporation (kg m-2 s-1)
   scalarSurfaceInfiltration => io_soilLiqFlx % scalarInfiltration,     & ! surface infiltration rate (m s-1)
   dq_dHydStateBelow         => io_soilLiqFlx % dq_dHydStateBelow,      & ! derivative in the flux in layer interfaces w.r.t. state variables in the layer below
   dq_dNrgStateBelow         => io_soilLiqFlx % dq_dNrgStateBelow       & ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
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
  call in_iLayerFlux % initialize(iLayer,nSoil,ibeg,iend,in_soilLiqFlx,io_soilLiqFlx,model_decisions,&
                                 &prog_data,mLayerDiffuse,dHydCond_dTemp,dHydCond_dVolLiq,dDiffuse_dVolLiq)
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
   err     => out_soilLiqFlx % err,                       & ! error code
   message => out_soilLiqFlx % cmessage                   & ! error message
  &)
   call out_iLayerFlux % finalize(iLayer,nSoil,io_soilLiqFlx,iLayerHydCond,iLayerDiffuse,err,cmessage)
   if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end associate
 end subroutine finalize_compute_interface_fluxes_derivatives

 subroutine compute_drainage_flux
  ! **** Compute the drainage flux from the bottom of the soil profile and its derivative ****
  type(in_type_qDrainFlux)  :: in_qDrainFlux
  type(out_type_qDrainFlux) :: out_qDrainFlux

  call initialize_compute_drainage_flux(in_qDrainFlux)

  call update_compute_drainage_flux(in_qDrainFlux,out_qDrainFlux)

  call finalize_compute_drainage_flux(out_qDrainFlux); if (return_flag) return

 end subroutine compute_drainage_flux

 subroutine initialize_compute_drainage_flux(in_qDrainFlux)
  ! **** Initialize operations for compute_drainage_flux ****
  type(in_type_qDrainFlux),intent(out) :: in_qDrainFlux
  call in_qDrainFlux % initialize(nSoil,ibeg,iend,in_soilLiqFlx,io_soilLiqFlx,model_decisions,&
                                 &prog_data,mpar_data,flux_data,diag_data,iceImpedeFac,&
                                 &dHydCond_dVolLiq,dHydCond_dTemp)
 end subroutine initialize_compute_drainage_flux
 
 subroutine update_compute_drainage_flux(in_qDrainFlux,out_qDrainFlux)
  ! **** Update operations for compute_drainage_flux ****
  type(in_type_qDrainFlux) ,intent(in)  :: in_qDrainFlux
  type(out_type_qDrainFlux),intent(out) :: out_qDrainFlux
  call qDrainFlux(in_qDrainFlux,out_qDrainFlux)
 end subroutine update_compute_drainage_flux

 subroutine finalize_compute_drainage_flux(out_qDrainFlux)
  ! **** finalize operations for compute_drainage_flux ****
  type(out_type_qDrainFlux),intent(in) :: out_qDrainFlux
  associate(&
   err     => out_soilLiqFlx % err,                       & ! error code
   message => out_soilLiqFlx % cmessage                   & ! error message
  &)
   call out_qDrainFlux % finalize(nSoil,io_soilLiqFlx,iLayerHydCond,iLayerDiffuse,err,cmessage)
   if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end associate

  ! no dependence on the aquifer for drainage
  associate(&
   ! derivatives in flux w.r.t. ...
   dq_dHydStateBelow => io_soilLiqFlx % dq_dHydStateBelow,& ! ... hydrology state variables in the layer below
   dq_dNrgStateBelow => io_soilLiqFlx % dq_dNrgStateBelow & ! ... temperature in the layer below (m s-1 K-1)
  &)
   dq_dHydStateBelow(nSoil) = 0._rkind  ! keep this here in case we want to couple some day....
   dq_dNrgStateBelow(nSoil) = 0._rkind  ! keep this here in case we want to couple some day....
  end associate
 end subroutine finalize_compute_drainage_flux
end subroutine soilLiqFlx

! ***************************************************************************************************************
! private subroutine diagv_node: compute transmittance and derivatives for model nodes
! ***************************************************************************************************************
subroutine diagv_node(in_diagv_node,out_diagv_node) 
  USE soil_utils_module,only:iceImpede            ! compute the ice impedence factor
  USE soil_utils_module,only:volFracLiq           ! compute volumetric fraction of liquid water as a function of matric head
  USE soil_utils_module,only:matricHead           ! compute matric head (m)
  USE soil_utils_module,only:hydCond_psi          ! compute hydraulic conductivity as a function of matric head
  USE soil_utils_module,only:hydCond_liq          ! compute hydraulic conductivity as a function of volumetric liquid water content
  USE soil_utils_module,only:hydCondMP_liq        ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
  USE soil_utils_module,only:dTheta_dPsi          ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
  USE soil_utils_module,only:dPsi_dTheta          ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  USE soil_utils_module,only:dPsi_dTheta2         ! compute derivative in dPsi_dTheta (m)
  USE soil_utils_module,only:dHydCond_dLiq        ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
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
  real(rkind)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
  real(rkind)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
  real(rkind)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
  real(rkind)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
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
   ! input: model control
   ixRichards    => in_diagv_node % ixRichards, & ! index defining the option for Richards' equation (moisture or mixdform)
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

   select case(ixRichards)
     case(moisture)
       scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
       scalardTheta_dPsi = realMissing  ! deliberately cause problems if this is ever used
     case(mixdform)
       scalardTheta_dPsi = dTheta_dPsi(scalarMatricHeadLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
       scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select

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

  associate(&
   ! input: model control
   ixRichards    => in_diagv_node % ixRichards   , & ! index defining the option for Richards' equation (moisture or mixdform)
   ! output: error control
   err     => out_diagv_node % err    , & ! error code
   message => out_diagv_node % message  & ! error message
  &)
   select case(ixRichards)
     case(moisture) ! moisture-based form of Richards' equation
       call update_diagv_node_hydraulic_conductivity_moisture_form; if (return_flag) return
     case(mixdform) ! mixed form of Richards' equation -- just compute hydraulic condictivity
       call update_diagv_node_hydraulic_conductivity_mixed_form;    if (return_flag) return
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select 
  end associate
 end subroutine update_diagv_node_hydraulic_conductivity

 subroutine update_diagv_node_hydraulic_conductivity_moisture_form
  ! **** Update operations for diagv_node: compute hydraulic conductivity and derivatives for moisture form of Richards' equation ****

  ! validation
  associate(&
   ! output: error control
   err     => out_diagv_node % err    , & ! error code
   message => out_diagv_node % message  & ! error message
  &)
   ! haven't included macropores yet -- return with error for now
   err=20; message=trim(message)//'still need to include macropores for the moisture-based form of Richards eqn'
   return_flag=.true.; return
  end associate

  ! computation
  associate(&
   ! input: model control
   deriv_desired => in_diagv_node % deriv_desired, & ! flag indicating if derivatives are desired
   ! input: state and diagnostic variables
   scalarVolFracLiqTrial    => in_diagv_node % scalarVolFracLiqTrial   , & ! volumetric fraction of liquid water in a given layer (-)
   scalarVolFracIceTrial    => in_diagv_node % scalarVolFracIceTrial   , & ! volumetric fraction of ice in a given layer (-)
   ! input: soil parameters
   vGn_alpha => in_diagv_node % vGn_alpha, & ! van Genuchten "alpha" parameter (m-1)
   vGn_n     => in_diagv_node % vGn_n    , & ! van Genuchten "n" parameter (-)
   vGn_m     => in_diagv_node % vGn_m    , & ! van Genuchten "m" parameter (-)
   theta_sat => in_diagv_node % theta_sat, & ! soil porosity (-)
   theta_res => in_diagv_node % theta_res, & ! soil residual volumetric water content (-)
   ! input: saturated hydraulic conductivity ...
   scalarSatHydCond   => in_diagv_node % scalarSatHydCond,  & ! ... at the mid-point of a given layer (m s-1)
   ! output: derivative in the soil water characteristic
   scalardPsi_dTheta => out_diagv_node % scalardPsi_dTheta, & ! derivative in the soil water characteristic
   ! output: transmittance
   scalarHydCond => out_diagv_node % scalarHydCond, & ! hydraulic conductivity at layer mid-points (m s-1)
   scalarDiffuse => out_diagv_node % scalarDiffuse, & ! diffusivity at layer mid-points (m2 s-1)
   iceImpedeFac  => out_diagv_node % iceImpedeFac , & ! ice impedence factor in each layer (-)
   ! output: transmittance derivatives in ...
   dHydCond_dVolLiq => out_diagv_node % dHydCond_dVolLiq, & ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   dDiffuse_dVolLiq => out_diagv_node % dDiffuse_dVolLiq, & ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   dHydCond_dMatric => out_diagv_node % dHydCond_dMatric  & ! ... hydraulic conductivity w.r.t matric head (s-1)
  &)

   ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
   hydCond_noIce = hydCond_liq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m)
   scalarHydCond = hydCond_noIce*iceImpedeFac
   scalarDiffuse = scalardPsi_dTheta * scalarHydCond
   ! compute derivative in hydraulic conductivity (m s-1) and hydraulic diffusivity (m2 s-1)
   if (deriv_desired) then
     if (scalarVolFracIceTrial > epsilon(iceImpedeFac)) then
       dK_dLiq__noIce   = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)  ! [.true. = analytical]
       dHydCond_dVolLiq = hydCond_noIce*dIceImpede_dLiq + dK_dLiq__noIce*iceImpedeFac
     else
       dHydCond_dVolLiq = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)
     end if
       dPsi_dTheta2a    = dPsi_dTheta2(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m,.true.)   ! [.true. = analytical] compute derivative in dPsi_dTheta (m)
       dDiffuse_dVolLiq = dHydCond_dVolLiq*scalardPsi_dTheta + scalarHydCond*dPsi_dTheta2a
       dHydCond_dMatric = realMissing ! not used, so cause problems
   end if

  end associate
 end subroutine update_diagv_node_hydraulic_conductivity_moisture_form

 subroutine update_diagv_node_hydraulic_conductivity_mixed_form 
  ! **** Update operations for diagv_node: compute hydraulic conductivity and derivatives for mixed form of Richards' equation ****
  associate(&
   ! input: model control
   deriv_desired => in_diagv_node % deriv_desired, & ! flag indicating if derivatives are desired
   ! input: state and diagnostic variables
   scalarMatricHeadLiqTrial => in_diagv_node % scalarMatricHeadLiqTrial, & ! liquid matric head in each layer (m)
   scalarVolFracIceTrial    => in_diagv_node % scalarVolFracIceTrial   , & ! volumetric fraction of ice in a given layer (-)
   ! input: pre-computed deriavatives
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
   dHydCond_dVolLiq => out_diagv_node % dHydCond_dVolLiq, & ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   dDiffuse_dVolLiq => out_diagv_node % dDiffuse_dVolLiq, & ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
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
   if (deriv_desired) then 
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
     if (scalarVolFracIceTrial > verySmall) then
       dK_dPsi__noIce        = dHydCond_dPsi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)  ! analytical
       dHydCondMicro_dTemp   = dPsiLiq_dTemp*dK_dPsi__noIce  ! m s-1 K-1
       dHydCondMicro_dMatric = hydCond_noIce*dIceImpede_dLiq*scalardTheta_dPsi + dK_dPsi__noIce*iceImpedeFac
     else
       dHydCondMicro_dTemp   = 0._rkind
       dHydCondMicro_dMatric = dHydCond_dPsi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)
     end if
     ! combine derivatives
     dHydCond_dMatric = dHydCondMicro_dMatric + dHydCondMacro_dMatric

     ! compute analytical derivative for change in ice impedance factor w.r.t. temperature
     call dIceImpede_dTemp(scalarVolFracIceTrial, & ! intent(in):  trial value of volumetric ice content (-)
                           dTheta_dTk,            & ! intent(in):  derivative in volumetric liquid water content w.r.t. temperature (K-1)
                           f_impede,              & ! intent(in):  ice impedance parameter (-)
                           dIceImpede_dT          ) ! intent(out): derivative in ice impedance factor w.r.t. temperature (K-1)
     ! compute derivative in hydraulic conductivity w.r.t. temperature
     dHydCond_dTemp = hydCond_noIce*dIceImpede_dT + dHydCondMicro_dTemp*iceImpedeFac
     ! set values that are not used to missing
     dHydCond_dVolLiq = realMissing ! not used, so cause problems
     dDiffuse_dVolLiq = realMissing ! not used, so cause problems
   end if

  end associate
 end subroutine update_diagv_node_hydraulic_conductivity_mixed_form

 subroutine finalize_diagv_node
  ! **** Finalize operations for diagv_node ****
  associate(&
   deriv_desired => in_diagv_node % deriv_desired,        & ! flag indicating if derivatives are desired
   dHydCond_dVolLiq => out_diagv_node % dHydCond_dVolLiq, & ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   dDiffuse_dVolLiq => out_diagv_node % dDiffuse_dVolLiq, & ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   dHydCond_dMatric => out_diagv_node % dHydCond_dMatric  & ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  &)
   ! if derivatives are not desired, then set values to missing
   if (.not.deriv_desired) then
     dHydCond_dVolLiq   = realMissing ! not used, so cause problems
     dDiffuse_dVolLiq   = realMissing ! not used, so cause problems
     dHydCond_dMatric   = realMissing ! not used, so cause problems
   end if
  end associate

  associate(&
   err     => out_diagv_node % err    , & ! error code
   message => out_diagv_node % message  & ! error message
  &)
   ! final error check
   if (err /= 0_i4b) then
    message=trim(message)//'unanticipated error in diagv_node'
    return_flag=.true.; return
   end if
  end associate
 end subroutine finalize_diagv_node

end subroutine diagv_node

! ***************************************************************************************************************
! private subroutine surfaceFlx: compute the surface flux and its derivative
! ***************************************************************************************************************
subroutine surfaceFlx(io_soilLiqFlx,in_surfaceFlx,io_surfaceFlx,out_surfaceFlx)
  USE soil_utils_module,only:volFracLiq            ! compute volumetric fraction of liquid water as a function of matric head (-)
  USE soil_utils_module,only:hydCond_psi           ! compute hydraulic conductivity as a function of matric head (m s-1)
  USE soil_utils_module,only:hydCond_liq           ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
  USE soil_utils_module,only:dPsi_dTheta           ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  USE soil_utils_module,only:crit_soilT            ! compute critical temperature below which ice exists
  USE soil_utils_module,only:gammp                 ! compute the cumulative probabilty based on the Gamma distribution
  ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
  implicit none
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! input: use soilLiqFlx object for array dimensions
  type(io_type_soilLiqFlx) ,intent(in)    :: io_soilLiqFlx          ! input-output object for soilLiqFlx
  ! input: model control, variables, derivatives, soil layer depth, boundary conditions, fluxes, and transmittance and soil parameters
  type(in_type_surfaceFlx) ,intent(in)    :: in_surfaceFlx          ! input object for surfaceFlx
  ! input-output: hydraulic conductivity and diffusivity, and infiltration parameters
  type(io_type_surfaceFlx) ,intent(inout) :: io_surfaceFlx          ! input object for surfaceFlx
  ! output: runoff, infiltration, derivatives, and error control
  type(out_type_surfaceFlx),intent(out)   :: out_surfaceFlx         ! output object for surfaceFlx
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! general
  integer(i4b)                     :: iLayer                              ! index of soil layer
  real(rkind)                      :: Tcrit                               ! temperature where all water is unfrozen (K)
  real(rkind)                      :: fPart1,fPart2                       ! different parts of a function
  real(rkind)                      :: dPart1(1:in_surfaceFlx % nSoil)     ! derivatives for different parts of a function
  real(rkind)                      :: dPart2(1:in_surfaceFlx % nSoil)     ! derivatives for different parts of a function
  real(rkind)                      :: dfracCap(1:in_surfaceFlx % nSoil)   ! derivatives for different parts of a function
  real(rkind)                      :: dfInfRaw(1:in_surfaceFlx % nSoil)   ! derivatives for different parts of a function
  real(rkind)                      :: total_soil_depth                    ! total depth of soil (m)
  ! head boundary condition
  real(rkind)                      :: cFlux                               ! capillary flux (m s-1)
  real(rkind)                      :: dNum                                ! numerical derivative
  ! simplified Green-Ampt infiltration
  real(rkind)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
  real(rkind)                      :: rootZoneIce                         ! depth of ice in the root zone (m)
  real(rkind)                      :: availCapacity                       ! available storage capacity in the root zone (m)
  real(rkind)                      :: depthWettingFront                   ! depth to the wetting front (m)
  real(rkind)                      :: hydCondWettingFront                 ! hydraulic conductivity at the wetting front (m s-1)
  ! saturated area associated with variable storage capacity
  real(rkind)                      :: fracCap                             ! fraction of pore space filled with liquid water and ice (-)
  real(rkind)                      :: fInfRaw                             ! infiltrating area before imposing solution constraints (-)
  real(rkind),parameter            :: maxFracCap=0.995_rkind              ! maximum fraction capacity -- used to avoid numerical problems associated with an enormous derivative
  real(rkind),parameter            :: scaleFactor=0.000001_rkind          ! scale factor for the smoothing function (-)
  real(rkind),parameter            :: qSurfScaleMax=1000._rkind           ! maximum surface runoff scaling factor (-)
  ! fraction of impermeable area associated with frozen ground
  real(rkind)                      :: alpha                               ! shape parameter in the Gamma distribution
  real(rkind)                      :: xLimg                               ! upper limit of the integral
  ! derivatives in ...
  real(rkind) :: dVolFracLiq_dWat(1:in_surfaceFlx % nSoil)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
  real(rkind) :: dVolFracIce_dWat(1:in_surfaceFlx % nSoil)  ! ... vol fraction of ice w.r.t. water state variable in root layers
  real(rkind) :: dVolFracLiq_dTk(1:in_surfaceFlx % nSoil)   ! ... vol fraction of liquid w.r.t. temperature in root layers
  real(rkind) :: dVolFracIce_dTk(1:in_surfaceFlx % nSoil)   ! ... vol fraction of ice w.r.t. temperature in root layers
  real(rkind) :: dRootZoneLiq_dWat(1:in_surfaceFlx % nSoil) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
  real(rkind) :: dRootZoneIce_dWat(1:in_surfaceFlx % nSoil) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
  real(rkind) :: dRootZoneLiq_dTk(1:in_surfaceFlx % nSoil)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
  real(rkind) :: dRootZoneIce_dTk(1:in_surfaceFlx % nSoil)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
  real(rkind) :: dDepthWettingFront_dWat(1:in_surfaceFlx % nSoil) ! ... scalar depth of wetting front w.r.t. water state variable in root layers
  real(rkind) :: dDepthWettingFront_dTk(1:in_surfaceFlx % nSoil)  ! ... scalar depth of wetting front w.r.t. temperature in root layers
  real(rkind) :: dxMaxInfilRate_dWat(1:in_surfaceFlx % nSoil) ! ... scalar max infiltration rate w.r.t. water state variable in root layers
  real(rkind) :: dxMaxInfilRate_dTk(1:in_surfaceFlx % nSoil)  ! ... scalar max infiltration rate w.r.t. temperature in root layers
  real(rkind) :: dInfilArea_dWat(0:in_surfaceFlx % nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dInfilArea_dTk(0:in_surfaceFlx % nSoil)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind) :: dFrozenArea_dWat(0:in_surfaceFlx % nSoil) ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dFrozenArea_dTk(0:in_surfaceFlx % nSoil)  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
  real(rkind) :: dInfilRate_dWat(0:in_surfaceFlx % nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dInfilRate_dTk(0:in_surfaceFlx % nSoil)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  ! error control
  logical(lgt) :: return_flag ! logical flag for return statements

  call initialize_surfaceFlx

  call update_surfaceFlx;   if (return_flag) return

  call finalize_surfaceFlx; if (return_flag) return

contains

 subroutine initialize_surfaceFlx
  ! **** Initialize operations for surfaceFlx ****
 
  ! allocate output object array components
  out_surfaceFlx % dq_dHydStateVec = io_soilLiqFlx % dq_dHydStateLayerSurfVec
  out_surfaceFlx % dq_dNrgStateVec = io_soilLiqFlx % dq_dNrgStateLayerSurfVec

  ! initialize error control
  return_flag=.false.
  associate(&
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)
   err=0; message="surfaceFlx/"
  end associate
 
  ! initialize derivatives
  associate(&
   ! output: derivatives in surface infiltration w.r.t. ...
   dq_dHydStateVec => out_surfaceFlx % dq_dHydStateVec , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   dq_dNrgStateVec => out_surfaceFlx % dq_dNrgStateVec   & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
  &)
   dq_dHydStateVec(:) = 0._rkind
   dq_dNrgStateVec(:) = 0._rkind
  end associate
 end subroutine initialize_surfaceFlx

 subroutine update_surfaceFlx
  ! **** Update operations for surfaceFlx ****
  associate(&
   ! input: model control
   bc_upper => in_surfaceFlx % bc_upper, & ! index defining the type of boundary conditions
   ! output: error control
   err      => out_surfaceFlx % err    , & ! error code
   message  => out_surfaceFlx % message  & ! error message
  &)

   ! compute the surface flux and its derivative
   select case(bc_upper)
 
     case(prescribedHead) ! head condition
       call update_surfaceFlx_prescribedHead; if (return_flag) return 
 
     case(liquidFlux) ! flux condition
       call update_surfaceFlx_liquidFlux;     if (return_flag) return 
 
     case default; err=20; message=trim(message)//'unknown upper boundary condition for soil hydrology'; return_flag=.true.; return
 
   end select 

  end associate
 end subroutine update_surfaceFlx

 subroutine update_surfaceFlx_prescribedHead
  ! **** Update operations for surfaceFlx: prescribed pressure head condition ****
  associate(&
   ! input: model control
   deriv_desired  => in_surfaceFlx % deriv_desired  , & ! flag to indicate if derivatives are desired
   ixRichards     => in_surfaceFlx % ixRichards     , & ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   scalarMatricHeadLiq => in_surfaceFlx % scalarMatricHeadLiq , & ! liquid matric head in the upper-most soil layer (m)
   scalarVolFracLiq    => in_surfaceFlx % scalarVolFracLiq    , & ! volumetric liquid water content in the upper-most soil layer (-)
   ! input: depth of upper-most soil layer (m)
   mLayerDepth  => in_surfaceFlx % mLayerDepth  , & ! depth of upper-most soil layer (m)
   ! input: diriclet boundary conditions
   upperBoundHead   => in_surfaceFlx % upperBoundHead  , & ! upper boundary condition for matric head (m)
   upperBoundTheta  => in_surfaceFlx % upperBoundTheta , & ! upper boundary condition for volumetric liquid water content (-)
   ! input: transmittance
   surfaceSatHydCond => in_surfaceFlx % surfaceSatHydCond , & ! saturated hydraulic conductivity at the surface (m s-1)
   dHydCond_dTemp    => in_surfaceFlx % dHydCond_dTemp    , & ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   iceImpedeFac      => in_surfaceFlx % iceImpedeFac      , & ! ice impedence factor in the upper-most soil layer (-)
   ! input: soil parameters
   vGn_alpha           => in_surfaceFlx % vGn_alpha           , & ! van Genuchten "alpha" parameter (m-1)
   vGn_n               => in_surfaceFlx % vGn_n               , & ! van Genuchten "n" parameter (-)
   vGn_m               => in_surfaceFlx % vGn_m               , & ! van Genuchten "m" parameter (-)
   theta_sat           => in_surfaceFlx % theta_sat           , & ! soil porosity (-)
   theta_res           => in_surfaceFlx % theta_res           , & ! soil residual volumetric water content (-)
   ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   surfaceHydCond => io_surfaceFlx % surfaceHydCond , & ! hydraulic conductivity (m s-1)
   surfaceDiffuse => io_surfaceFlx % surfaceDiffuse , & ! hydraulic diffusivity at the surface (m
   ! output: runoff and infiltration
   scalarSurfaceRunoff       => out_surfaceFlx % scalarSurfaceRunoff       , & ! surface runoff (m s-1)
   scalarSurfaceInfiltration => out_surfaceFlx % scalarSurfaceInfiltration , & ! surface infiltration (m s-1)
   ! output: derivatives in surface infiltration w.r.t. ...
   dq_dHydStateVec => out_surfaceFlx % dq_dHydStateVec , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   dq_dNrgStateVec => out_surfaceFlx % dq_dNrgStateVec , & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
   ! output: error control
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)

   ! surface runoff iz zero for the head condition
   scalarSurfaceRunoff = 0._rkind

   ! compute transmission and the capillary flux
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture)
       ! compute the hydraulic conductivity and diffusivity at the boundary
       surfaceHydCond = hydCond_liq(upperBoundTheta,surfaceSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
       surfaceDiffuse = dPsi_dTheta(upperBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * surfaceHydCond
       ! compute the capillary flux
       cflux = -surfaceDiffuse*(scalarVolFracLiq - upperBoundTheta) / (mLayerDepth(1)*0.5_rkind)
     case(mixdform)
       ! compute the hydraulic conductivity and diffusivity at the boundary
       surfaceHydCond = hydCond_psi(upperBoundHead,surfaceSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
       surfaceDiffuse = realMissing
       ! compute the capillary flux
       cflux = -surfaceHydCond*(scalarMatricHeadLiq - upperBoundHead) / (mLayerDepth(1)*0.5_rkind)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select  ! end select form of Richards' eqn

   ! compute the total flux
   scalarSurfaceInfiltration = cflux + surfaceHydCond

   ! compute the derivative
   if (deriv_desired) then
     ! compute the hydrology derivative at the surface
     select case(ixRichards)  ! select form of Richards' equation
       case(moisture); dq_dHydStateVec(1) = -surfaceDiffuse/(mLayerDepth(1)/2._rkind)
       case(mixdform); dq_dHydStateVec(1) = -surfaceHydCond/(mLayerDepth(1)/2._rkind)
       case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
     end select
     ! compute the energy derivative at the surface
     dq_dNrgStateVec(1) = -(dHydCond_dTemp/2._rkind)*(scalarMatricHeadLiq - upperBoundHead)/(mLayerDepth(1)*0.5_rkind) + dHydCond_dTemp/2._rkind
   else
     dNum = 0._rkind
   end if

  end associate
 end subroutine update_surfaceFlx_prescribedHead

 subroutine update_surfaceFlx_liquidFlux 
  ! **** Update operations for surfaceFlx: flux condition ****

  ! force infiltration to be constant over the iterations
  associate(&
   ! input: model control
   firstSplitOper => in_surfaceFlx % firstSplitOper , & ! flag indicating if desire to compute infiltration
   ! output: derivatives in surface infiltration w.r.t. ...
   dq_dHydStateVec => out_surfaceFlx % dq_dHydStateVec , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   dq_dNrgStateVec => out_surfaceFlx % dq_dNrgStateVec   & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
  &)
   if (firstSplitOper) then
     call update_surfaceFlx_liquidFlux_computation; if (return_flag) return 
   else ! do not compute infiltration after first flux call in a splitting operation
     dq_dHydStateVec(:) = 0._rkind
     dq_dNrgStateVec(:) = 0._rkind
   end if 
  end associate

  call update_surfaceFlx_liquidFlux_infiltration ! final computations for infiltration and runoff
 end subroutine update_surfaceFlx_liquidFlux

 subroutine update_surfaceFlx_liquidFlux_computation 
  ! **** Update operations for surfaceFlx: flux condition -- main computations ****

   call update_surfaceFlx_liquidFlux_computation_root_layers 

   call update_surfaceFlx_liquidFlux_computation_available_capacity; if (return_flag) return 

   call update_surfaceFlx_liquidFlux_computation_wetting_front

   call update_surfaceFlx_liquidFlux_computation_infiltrating_area

   call update_surfaceFlx_liquidFlux_computation_validate_infiltration
 
   call update_surfaceFlx_liquidFlux_computation_impermeable_area

   call update_surfaceFlx_liquidFlux_computation_flux_derivatives

 end subroutine update_surfaceFlx_liquidFlux_computation 

 subroutine update_surfaceFlx_liquidFlux_computation_root_layers 
  ! **** Update operations for surfaceFlx: flux condition -- main computations (root layers) ****
  associate(&
   ! input: model control
   deriv_desired  => in_surfaceFlx % deriv_desired  , & ! flag to indicate if derivatives are desired
   ixRichards     => in_surfaceFlx % ixRichards     , & ! index defining the option for Richards' equation (moisture or mixdform)
   nRoots         => in_surfaceFlx % nRoots         , & ! number of layers that contain roots
   ! input: state and diagnostic variables
   mLayerTemp          => in_surfaceFlx % mLayerTemp          , & ! temperature (K)
   mLayerMatricHead    => in_surfaceFlx % mLayerMatricHead    , & ! matric head in each soil layer (m)
   mLayerVolFracLiq    => in_surfaceFlx % mLayerVolFracLiq    , & ! volumetric liquid water content in each soil layer (-)
   mLayerVolFracIce    => in_surfaceFlx % mLayerVolFracIce    , & ! volumetric ice content in each soil layer (-)
   ! input: pre-computed derivatives in ...
   ! note: all of these would need to be recomputed if wanted a numerical derivative
   dTheta_dTk             => in_surfaceFlx % dTheta_dTk             , & ! ... volumetric liquid water content w.r.t. temperature (K-1)
   dTheta_dPsi            => in_surfaceFlx % dTheta_dPsi            , & ! ... the soil water characteristic w.r.t. psi (m-1)
   mLayerdPsi_dTheta      => in_surfaceFlx % mLayerdPsi_dTheta      , & ! ... the soil water characteristic w.r.t. theta (m)
   ! input: depth of upper-most soil layer (m)
   mLayerDepth  => in_surfaceFlx % mLayerDepth  , & ! depth of upper-most soil layer (m)
   iLayerHeight => in_surfaceFlx % iLayerHeight , & ! height at the interface of each layer (m)
   ! input: soil parameters
   rootingDepth        => in_surfaceFlx % rootingDepth & ! rooting depth (m)
  &)

   ! process root layers only liquid and ice derivatives
   dVolFracLiq_dWat(:) = 0._rkind
   dVolFracIce_dWat(:) = 0._rkind
   dVolFracLiq_dTk(:)  = 0._rkind
   dVolFracIce_dTk(:)  = 0._rkind
   if (deriv_desired .and. nRoots > 0) then
     select case(ixRichards)  ! form of Richards' equation
       case(moisture)
         dVolFracLiq_dWat(:) = 1._rkind
         dVolFracIce_dWat(:) = mLayerdPsi_dTheta(:) - 1._rkind
       case(mixdform)
         do iLayer=1,nRoots
           Tcrit = crit_soilT( mLayerMatricHead(iLayer) )
           if (mLayerTemp(iLayer) < Tcrit) then
             dVolFracLiq_dWat(iLayer) = 0._rkind
             dVolFracIce_dWat(iLayer) = dTheta_dPsi(iLayer)
           else
             dVolFracLiq_dWat(iLayer) = dTheta_dPsi(iLayer)
             dVolFracIce_dWat(iLayer) = 0._rkind
           end if
         end do
     end select 
     dVolFracLiq_dTk(:) = dTheta_dTk(:) !already zeroed out if not below critical temperature
     dVolFracIce_dTk(:) = -dVolFracLiq_dTk(:) !often can and will simplify one of these terms out
   end if
 
   ! define the storage in the root zone (m) and derivatives
   rootZoneLiq = 0._rkind
   rootZoneIce = 0._rkind
   dRootZoneLiq_dWat(:) = 0._rkind
   dRootZoneIce_dWat(:) = 0._rkind
   dRootZoneLiq_dTk(:)  = 0._rkind
   dRootZoneIce_dTk(:)  = 0._rkind
 
   ! process layers where the roots extend to the bottom of the layer
   if (nRoots > 1) then
     do iLayer=1,nRoots-1
       rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(iLayer)*mLayerDepth(iLayer)
       rootZoneIce = rootZoneIce + mLayerVolFracIce(iLayer)*mLayerDepth(iLayer)
       dRootZoneLiq_dWat(iLayer) = dVolFracLiq_dWat(iLayer)*mLayerDepth(iLayer)
       dRootZoneIce_dWat(iLayer) = dVolFracIce_dWat(iLayer)*mLayerDepth(iLayer)
       dRootZoneLiq_dTk(iLayer)  = dVolFracLiq_dTk(iLayer) *mLayerDepth(iLayer)
       dRootZoneIce_dTk(iLayer)  = dVolFracIce_dTk(iLayer) *mLayerDepth(iLayer)
     end do
   end if
   ! process layers where the roots end in the current layer
   rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
   rootZoneIce = rootZoneIce + mLayerVolFracIce(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
   dRootZoneLiq_dWat(nRoots) = dVolFracLiq_dWat(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
   dRootZoneIce_dWat(nRoots) = dVolFracIce_dWat(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
   dRootZoneLiq_dTk(nRoots)  = dVolFracLiq_dTk(nRoots)* (rootingDepth - iLayerHeight(nRoots-1))
   dRootZoneIce_dTk(nRoots)  = dVolFracIce_dTk(nRoots)* (rootingDepth - iLayerHeight(nRoots-1))

  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_root_layers 

 subroutine update_surfaceFlx_liquidFlux_computation_available_capacity 
  ! **** Update operations for surfaceFlx: flux condition -- main computations (check available capacity) ****

  ! compute and check available capacity to hold water (m)
  associate(&
   ! input: soil parameters
   theta_sat           => in_surfaceFlx % theta_sat   , & ! soil porosity (-)
   rootingDepth        => in_surfaceFlx % rootingDepth, & ! rooting depth (m)
   ! output: error control
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)

   availCapacity = theta_sat*rootingDepth - rootZoneIce
   if (rootZoneLiq > availCapacity+verySmall) then
     message=trim(message)//'liquid water in the root zone exceeds capacity'
     err=20; return_flag=.true.; return
   end if

  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_available_capacity 

 subroutine update_surfaceFlx_liquidFlux_computation_wetting_front
  ! **** Update operations for surfaceFlx: flux condition -- main computations (wetting front and derivatives) ****
  associate(&
   ! input: depth of upper-most soil layer (m)
   mLayerDepth  => in_surfaceFlx % mLayerDepth  , & ! depth of upper-most soil layer (m)
   ! input: transmittance
   surfaceSatHydCond => in_surfaceFlx % surfaceSatHydCond , & ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: soil parameters
   zScale_TOPMODEL     => in_surfaceFlx % zScale_TOPMODEL     , & ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
   rootingDepth        => in_surfaceFlx % rootingDepth        , & ! rooting depth (m)
   wettingFrontSuction => in_surfaceFlx % wettingFrontSuction , & ! Green-Ampt wetting front suction (m)
   ! input-output: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate    => io_surfaceFlx % xMaxInfilRate  & ! maximum infiltration rate (m s-1)
  &)

   ! define the depth to the wetting front (m) and derivatives
   total_soil_depth = sum(mLayerDepth)
   depthWettingFront = (rootZoneLiq/availCapacity)*min(rootingDepth, total_soil_depth)
   dDepthWettingFront_dWat(:)=( dRootZoneLiq_dWat(:)*min(rootingDepth, total_soil_depth) + dRootZoneIce_dWat(:)*depthWettingFront )/availCapacity
   dDepthWettingFront_dTk(:) =( dRootZoneLiq_dTk(:) *min(rootingDepth, total_soil_depth) + dRootZoneIce_dTk(:)*depthWettingFront  )/availCapacity
 
   ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
   hydCondWettingFront =  surfaceSatHydCond * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 1._rkind) )
 
   ! define the maximum infiltration rate (m s-1) and derivatives
   xMaxInfilRate = hydCondWettingFront*( (wettingFrontSuction + depthWettingFront)/depthWettingFront )  ! maximum infiltration rate (m s-1)
   fPart1    = hydCondWettingFront
   fPart2    = (wettingFrontSuction + depthWettingFront)/depthWettingFront
   dPart1(:) = surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dWat(:))/total_soil_depth
   dPart2(:) = -dDepthWettingFront_dWat(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
   dxMaxInfilRate_dWat(:) = fPart1*dPart2(:) + fPart2*dPart1(:)
   dPart1(:) = surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dTk(:))/total_soil_depth
   dPart2(:) = -dDepthWettingFront_dTk(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
   dxMaxInfilRate_dTk(:)  = fPart1*dPart2(:) + fPart2*dPart1(:)

  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_wetting_front

 subroutine update_surfaceFlx_liquidFlux_computation_infiltrating_area
  ! **** Update operations for surfaceFlx: flux condition -- main computations (infiltrating area) ****
  associate(&
   ! input: model control
   nSoil          => in_surfaceFlx % nSoil           , & ! number of soil layers
   ! input: soil parameters
   qSurfScale       => in_surfaceFlx % qSurfScale    , & ! scaling factor in the surface runoff parameterization (-)
   ! input-output: surface runoff and infiltration flux (m s-1)
   scalarInfilArea  => io_surfaceFlx % scalarInfilArea & ! fraction of unfrozen area where water can infiltrate (-)
  &)
   ! define the infiltrating area and derivatives for the non-frozen part of the cell/basin
   if (qSurfScale < qSurfScaleMax) then
     fracCap         = rootZoneLiq/(maxFracCap*availCapacity)                              ! fraction of available root zone filled with water
     fInfRaw         = 1._rkind - exp(-qSurfScale*(1._rkind - fracCap))                          ! infiltrating area -- allowed to violate solution constraints
     scalarInfilArea = min(0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor)), 1._rkind)   ! infiltrating area -- constrained
     if (0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor))< 1._rkind) then
       dfracCap(:) = ( dRootZoneLiq_dWat(:)/maxFracCap + dRootZoneIce_dWat(:)*fracCap )/availCapacity
       dfInfRaw(:) = -qSurfScale*dfracCap(:) * exp(-qSurfScale*(1._rkind - fracCap))
       dInfilArea_dWat(1:nSoil) = 0.5_rkind*dfInfRaw(:) * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
       dfracCap(:) = ( dRootZoneLiq_dTk(:)/maxFracCap + dRootZoneIce_dTk(:)*fracCap )/availCapacity
       dfInfRaw(:) = -qSurfScale*dfracCap(:) * exp(-qSurfScale*(1._rkind - fracCap))
       dInfilArea_dTk(1:nSoil)  = 0.5_rkind*dfInfRaw(:) * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
     else ! scalarInfilArea = 1._rkind
       dInfilArea_dWat(1:nSoil) = 0._rkind
       dInfilArea_dTk(1:nSoil)  = 0._rkind
     end if
   else
     scalarInfilArea = 1._rkind
     dInfilArea_dWat(1:nSoil) = 0._rkind
     dInfilArea_dTk(1:nSoil)  = 0._rkind
   end if
   dInfilArea_dWat(0) = 0._rkind
   dInfilArea_dTk(0)  = 0._rkind
  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_infiltrating_area

 subroutine update_surfaceFlx_liquidFlux_computation_validate_infiltration
  ! **** Update operations for surfaceFlx: flux condition -- main computations (validate infiltration) ****
  associate(&
   ! input: model control
   nRoots         => in_surfaceFlx % nRoots, & ! number of layers that contain roots
   ixIce          => in_surfaceFlx % ixIce , & ! index of lowest ice layer
   ! input: state and diagnostic variables
   mLayerVolFracLiq    => in_surfaceFlx % mLayerVolFracLiq, & ! volumetric liquid water content in each soil layer (-)
   ! input: depth of upper-most soil layer (m)
   mLayerDepth  => in_surfaceFlx % mLayerDepth, & ! depth of upper-most soil layer (m)
   ! input: soil parameters
   theta_sat           => in_surfaceFlx % theta_sat, & ! soil porosity (-)
   ! input-output: surface runoff and infiltration flux (m s-1)
   scalarInfilArea  => io_surfaceFlx % scalarInfilArea & ! fraction of unfrozen area where water can infiltrate (-)
  &)
   ! check to ensure we are not infiltrating into a fully saturated column
   if (ixIce<nRoots) then
     if (sum(mLayerVolFracLiq(ixIce+1:nRoots)*mLayerDepth(ixIce+1:nRoots)) > 0.9999_rkind*theta_sat*sum(mLayerDepth(ixIce+1:nRoots))) scalarInfilArea=0._rkind
   end if
  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_validate_infiltration

 subroutine update_surfaceFlx_liquidFlux_computation_impermeable_area
  ! **** Update operations for surfaceFlx: flux condition -- main computations (impermeable area) ****
  associate(&
   ! input: model control
   nSoil          => in_surfaceFlx % nSoil , & ! number of soil layers
   ! input: pre-computed derivatives in ...
   ! note: all of these would need to be recomputed if wanted a numerical derivative
   above_soilLiqFluxDeriv => in_surfaceFlx % above_soilLiqFluxDeriv , & ! ... layer above soil (canopy or snow) liquid flux w.r.t. liquid water
   above_soildLiq_dTk     => in_surfaceFlx % above_soildLiq_dTk     , & ! ... layer above soil (canopy or snow) liquid flux w.r.t. temperature
   above_soilFracLiq      => in_surfaceFlx % above_soilFracLiq      , & ! ... liquid water layer above soil (canopy or snow) (-)
   ! input: flux at the upper boundary
   scalarRainPlusMelt => in_surfaceFlx % scalarRainPlusMelt , & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input: soil parameters
   soilIceScale        => in_surfaceFlx % soilIceScale        , & ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   soilIceCV           => in_surfaceFlx % soilIceCV           , & ! soil ice CV in Gamma distribution used to define frozen area (-)
   ! input-output: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate    => io_surfaceFlx % xMaxInfilRate    , & ! maximum infiltration rate (m s-1)
   scalarFrozenArea => io_surfaceFlx % scalarFrozenArea   & ! fraction of area that is considered impermeable due to soil ice (-)
  &)
   ! define the impermeable area and derivatives due to frozen ground
   if (rootZoneIce > tiny(rootZoneIce)) then  ! (avoid divide by zero)
     alpha            = 1._rkind/(soilIceCV**2_i4b)        ! shape parameter in the Gamma distribution
     xLimg            = alpha*soilIceScale/rootZoneIce  ! upper limit of the integral
 
     !if we use this, we will have a derivative of scalarFrozenArea w.r.t. water and temperature in each layer (through mLayerVolFracIce)
     scalarFrozenArea = 0._rkind
     dFrozenArea_dWat(1:nSoil) = 0._rkind
     dFrozenArea_dTk(1:nSoil)  = 0._rkind
   else
     scalarFrozenArea = 0._rkind
     dFrozenArea_dWat(1:nSoil) = 0._rkind
     dFrozenArea_dTk(1:nSoil)  = 0._rkind
   end if
   dFrozenArea_dWat(0) = 0._rkind
   dFrozenArea_dTk(0)  = 0._rkind
 
   if (xMaxInfilRate < scalarRainPlusMelt) then ! = dxMaxInfilRate_d, dependent on layers not at surface
     dInfilRate_dWat(0) = 0._rkind
     dInfilRate_dTk(0)  = 0._rkind
     dInfilRate_dWat(1:nSoil) = dxMaxInfilRate_dWat(:)
     dInfilRate_dTk(1:nSoil)  = dxMaxInfilRate_dTk(:)
   else ! = dRainPlusMelt_d, dependent on above layer (canopy or snow) water and temp
     dInfilRate_dWat(0) = above_soilLiqFluxDeriv*above_soilFracLiq
     dInfilRate_dTk(0)  = above_soilLiqFluxDeriv*above_soildLiq_dTk
     dInfilRate_dWat(1:nSoil) = 0._rkind
     dInfilRate_dTk(1:nSoil)  = 0._rkind
   end if
  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_impermeable_area

 subroutine update_surfaceFlx_liquidFlux_computation_flux_derivatives
  ! **** Update operations for surfaceFlx: flux condition -- main computations (flux derivatives) ****
  associate(&
   ! input: flux at the upper boundary
   scalarRainPlusMelt => in_surfaceFlx % scalarRainPlusMelt , & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input-output: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate    => io_surfaceFlx % xMaxInfilRate    , & ! maximum infiltration rate (m s-1)
   scalarInfilArea  => io_surfaceFlx % scalarInfilArea  , & ! fraction of unfrozen area where water can infiltrate (-)
   scalarFrozenArea => io_surfaceFlx % scalarFrozenArea , & ! fraction of area that is considered impermeable due to soil ice (-)
   ! output: derivatives in surface infiltration w.r.t. ...
   dq_dHydStateVec => out_surfaceFlx % dq_dHydStateVec  , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   dq_dNrgStateVec => out_surfaceFlx % dq_dNrgStateVec    & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
  &)
   ! dq w.r.t. infiltration only, scalarRainPlusMelt accounted for in computJacob module
   dq_dHydStateVec(:) = (1._rkind - scalarFrozenArea)&
                      & * ( dInfilArea_dWat(:)*min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dWat(:) )&
                      & + (-dFrozenArea_dWat(:))*scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
   dq_dNrgStateVec(:) = (1._rkind - scalarFrozenArea)&
                      & * ( dInfilArea_dTk(:) *min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dTk(:)  )&
                      & + (-dFrozenArea_dTk(:)) *scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_flux_derivatives

 subroutine update_surfaceFlx_liquidFlux_infiltration
  ! **** Update operations for surfaceFlx: flux condition -- final infiltration and runoff calculations ****
  associate(&
   ! input: flux at the upper boundary
   scalarRainPlusMelt => in_surfaceFlx % scalarRainPlusMelt , & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   surfaceHydCond => io_surfaceFlx % surfaceHydCond , & ! hydraulic conductivity (m s-1)
   surfaceDiffuse => io_surfaceFlx % surfaceDiffuse , & ! hydraulic diffusivity at the surface (m
   ! input-output: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate    => io_surfaceFlx % xMaxInfilRate    , & ! maximum infiltration rate (m s-1)
   scalarInfilArea  => io_surfaceFlx % scalarInfilArea  , & ! fraction of unfrozen area where water can infiltrate (-)
   scalarFrozenArea => io_surfaceFlx % scalarFrozenArea , & ! fraction of area that is considered impermeable due to soil ice (-)
   ! output: runoff and infiltration
   scalarSurfaceRunoff       => out_surfaceFlx % scalarSurfaceRunoff       , & ! surface runoff (m s-1)
   scalarSurfaceInfiltration => out_surfaceFlx % scalarSurfaceInfiltration   & ! surface infiltration (m s-1)
  &)
   ! compute infiltration (m s-1), if after first flux call in a splitting operation does not change
   scalarSurfaceInfiltration = (1._rkind - scalarFrozenArea)*scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
 
   ! compute surface runoff (m s-1)
   scalarSurfaceRunoff = scalarRainPlusMelt - scalarSurfaceInfiltration
 
   ! set surface hydraulic conductivity and diffusivity to missing (not used for flux condition)
   surfaceHydCond = realMissing
   surfaceDiffuse = realMissing
  end associate

 end subroutine update_surfaceFlx_liquidFlux_infiltration

 subroutine finalize_surfaceFlx
  ! **** Finalize operations for surfaceFlx ****
  ! final error check
  associate(&
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)
   if (err /= 0_i4b) then
    message=trim(message)//'unanticipated error in surfaceFlx subroutine'; return_flag=.true.; return
   end if
  end associate
 end subroutine finalize_surfaceFlx

end subroutine surfaceFlx

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
  real(rkind)                      :: dLiq                 ! spatial difference in volumetric liquid water (-)
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
  if (in_iLayerFlux % deriv_desired) then
    call update_iLayerFlux_derivatives; if (return_flag) return
  else
   ! output: derivatives in fluxes w.r.t. ...  
   out_iLayerFlux % dq_dHydStateAbove = realMissing ! ... matric head or volumetric lquid water in the layer above (m s-1 or s-1)
   out_iLayerFlux % dq_dHydStateBelow = realMissing ! ... matric head or volumetric lquid water in the layer below (m s-1 or s-1)
  end if
 end subroutine update_iLayerFlux
 
 subroutine update_iLayerFlux_fluxes
  ! **** Update operations for iLayerFlux: compute fluxes ****
  associate(&
   ! input: model control
   ixRichards    => in_iLayerFlux % ixRichards   , & ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state variables
   nodeMatricHeadLiqTrial => in_iLayerFlux % nodeMatricHeadLiqTrial, & ! liquid matric head at the soil nodes (m)
   nodeVolFracLiqTrial    => in_iLayerFlux % nodeVolFracLiqTrial   , & ! volumetric fraction of liquid water at the soil nodes (-)
   ! input: model coordinate variables
   nodeHeight => in_iLayerFlux % nodeHeight, & ! height at the mid-point of the lower layer (m)
   ! input: transmittance
   nodeHydCondTrial => in_iLayerFlux % nodeHydCondTrial, & ! hydraulic conductivity at layer mid-points (m s-1)
   nodeDiffuseTrial => in_iLayerFlux % nodeDiffuseTrial, & ! diffusivity at layer mid-points (m2 s-1)
   ! output: tranmsmittance at the layer interface (scalars)
   iLayerHydCond => out_iLayerFlux % iLayerHydCond, & ! hydraulic conductivity at the interface between layers (m s-1)
   iLayerDiffuse => out_iLayerFlux % iLayerDiffuse, & ! hydraulic diffusivity at the interface between layers (m2 s-1)
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
   ! compute the capillary flux
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture)
      iLayerDiffuse = sqrt(nodeDiffuseTrial(ixLower) * nodeDiffuseTrial(ixUpper))
      dLiq          = nodeVolFracLiqTrial(ixLower) - nodeVolFracLiqTrial(ixUpper)
      cflux         = -iLayerDiffuse * dLiq/dz
     case(mixdform)
      iLayerDiffuse = realMissing
      dPsi          = nodeMatricHeadLiqTrial(ixLower) - nodeMatricHeadLiqTrial(ixUpper)
      cflux         = -iLayerHydCond * dPsi/dz
     case default; err=10; message=trim(message)//"unable to identify option for Richards' equation"; return_flag=.true.; return
   end select
   ! compute the total flux (add gravity flux, positive downwards)
   iLayerLiqFluxSoil = cflux + iLayerHydCond

  end associate
 end subroutine update_iLayerFlux_fluxes

 subroutine update_iLayerFlux_derivatives
  ! **** Update operations for iLayerFlux: compute derivatives ****
  ! * local variables (derivative in Darcy's flux) *
  ! deriviatives at the layer interface
  real(rkind) :: dHydCondIface_dVolLiqAbove  ! hydraulic conductivity w.r.t. volumetric liquid water content in layer above
  real(rkind) :: dHydCondIface_dVolLiqBelow  ! hydraulic conductivity w.r.t. volumetric liquid water content in layer below
  real(rkind) :: dDiffuseIface_dVolLiqAbove  ! hydraulic diffusivity  w.r.t. volumetric liquid water content in layer above
  real(rkind) :: dDiffuseIface_dVolLiqBelow  ! hydraulic diffusivity  w.r.t. volumetric liquid water content in layer below
  real(rkind) :: dHydCondIface_dMatricAbove  ! hydraulic conductivity w.r.t. matric head in layer above
  real(rkind) :: dHydCondIface_dMatricBelow  ! hydraulic conductivity w.r.t. matric head in layer below
  associate(&
   ! input: model control
   ixRichards    => in_iLayerFlux % ixRichards   , & ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: temperature derivatives
   dPsiLiq_dTemp   => in_iLayerFlux % dPsiLiq_dTemp , & ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   dHydCond_dTemp  => in_iLayerFlux % dHydCond_dTemp, & ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: transmittance
   nodeHydCondTrial => in_iLayerFlux % nodeHydCondTrial, & ! hydraulic conductivity at layer mid-points (m s-1)
   nodeDiffuseTrial => in_iLayerFlux % nodeDiffuseTrial, & ! diffusivity at layer mid-points (m2 s-1)
   ! input: transmittance derivatives
   dHydCond_dVolLiq => in_iLayerFlux % dHydCond_dVolLiq, & ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   dDiffuse_dVolLiq => in_iLayerFlux % dDiffuse_dVolLiq, & ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   dHydCond_dMatric => in_iLayerFlux % dHydCond_dMatric, & ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
   ! output: tranmsmittance at the layer interface (scalars)
   iLayerHydCond => out_iLayerFlux % iLayerHydCond, & ! hydraulic conductivity at the interface between layers (m s-1)
   iLayerDiffuse => out_iLayerFlux % iLayerDiffuse, & ! hydraulic diffusivity at the interface between layers (m2 s-1)
   ! output: derivatives in fluxes w.r.t. ...  
   dq_dHydStateAbove => out_iLayerFlux % dq_dHydStateAbove, & ! ... matric head or volumetric lquid water in the layer above (m s-1 or s-1)
   dq_dHydStateBelow => out_iLayerFlux % dq_dHydStateBelow, & ! ... matric head or volumetric lquid water in the layer below (m s-1 or s-1)
   ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
   dq_dNrgStateAbove => out_iLayerFlux % dq_dNrgStateAbove, & ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
   dq_dNrgStateBelow => out_iLayerFlux % dq_dNrgStateBelow, & ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
   ! output: error control
   err     => out_iLayerFlux % err    , & ! error code
   message => out_iLayerFlux % message  & ! error message
  &)

   select case(ixRichards)  ! select form of Richards' equation
     case(moisture)
       ! still need to implement arithmetric mean for the moisture-based form
       if (.not.useGeometric) then
         message=trim(message)//'only currently implemented for geometric mean -- change local flag'
         err=20; return_flag=.true.; return
       end if
       ! derivatives in hydraulic conductivity at the layer interface (m s-1)
       dHydCondIface_dVolLiqAbove = dHydCond_dVolLiq(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_rkind/max(iLayerHydCond,verySmall)
       dHydCondIface_dVolLiqBelow = dHydCond_dVolLiq(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_rkind/max(iLayerHydCond,verySmall)
       ! derivatives in hydraulic diffusivity at the layer interface (m2 s-1)
       dDiffuseIface_dVolLiqAbove = dDiffuse_dVolLiq(ixUpper)*nodeDiffuseTrial(ixLower) * 0.5_rkind/max(iLayerDiffuse,verySmall)
       dDiffuseIface_dVolLiqBelow = dDiffuse_dVolLiq(ixLower)*nodeDiffuseTrial(ixUpper) * 0.5_rkind/max(iLayerDiffuse,verySmall)
       ! derivatives in the flux w.r.t. volumetric liquid water content
       dq_dHydStateAbove = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse/dz + dHydCondIface_dVolLiqAbove
       dq_dHydStateBelow = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse/dz + dHydCondIface_dVolLiqBelow
     case(mixdform)
       ! derivatives in hydraulic conductivity
       if (useGeometric) then
         dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_rkind/max(iLayerHydCond,verySmall)
         dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_rkind/max(iLayerHydCond,verySmall)
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
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select

  end associate
 end subroutine update_iLayerFlux_derivatives

 subroutine finalize_iLayerFlux
  ! **** Finalize operations for iLayerFlux ****
  associate(&
   err     => out_iLayerFlux % err    , & ! error code
   message => out_iLayerFlux % message  & ! error message
  &)
   ! final error check
   if (err /= 0_i4b) then
    message=trim(message)//'unanticipated error in iLayerFlux'
    return_flag=.true.; return
   end if
  end associate
 end subroutine finalize_iLayerFlux
 
end subroutine iLayerFlux

! ***************************************************************************************************************
! private subroutine qDrainFlux: compute the drainage flux from the bottom of the soil profile and its derivative
! ***************************************************************************************************************
subroutine qDrainFlux(in_qDrainFlux,out_qDrainFlux)
  USE soil_utils_module,only:volFracLiq  ! compute volumetric fraction of liquid water as a function of matric head (-)
  USE soil_utils_module,only:matricHead  ! compute matric head as a function of volumetric fraction of liquid water (m)
  USE soil_utils_module,only:hydCond_psi ! compute hydraulic conductivity as a function of matric head (m s-1)
  USE soil_utils_module,only:hydCond_liq ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
  USE soil_utils_module,only:dPsi_dTheta ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
  implicit none
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! input: model control, variables, boundary conditions, transmittance variables, and soil parameters
  type(in_type_qDrainFlux) ,intent(in)  :: in_qDrainFlux      ! object for qDrainFlux input data
  ! output: hydraulic conductivity and diffusivity, drainage fluxes and derivatives, and error control
  type(out_type_qDrainFlux),intent(out) :: out_qDrainFlux     ! object for qDrainFlux output data
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! local variables
  real(rkind)                      :: zWater                  ! effective water table depth (m)
  real(rkind)                      :: nodePsi                 ! matric head in the lowest unsaturated node (m)
  real(rkind)                      :: cflux                   ! capillary flux (m s-1)
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
   ! output: error control
   err     => out_qDrainFlux % err    , &       ! error code
   message => out_qDrainFlux % message  &       ! error message
  &)

   ! determine lower boundary condition
   select case(bc_lower)
     case(prescribedHead) ! specified matric head value
       call update_qDrainFlux_prescribedHead; if (return_flag) return
     case(funcBottomHead) ! specified matric head function
       call update_qDrainFlux_funcBottomHead; if (return_flag) return
     case(freeDrainage)   ! free drainage 
       call update_qDrainFlux_freeDrainage;   if (return_flag) return
     case(zeroFlux)       ! zero flux
       call update_qDrainFlux_zeroFlux;       if (return_flag) return
     case default;
        err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return_flag=.true.; return
   end select 

  end associate
 end subroutine update_qDrainFlux

 subroutine update_qDrainFlux_prescribedHead
  ! ** Update operations for qDrainFlux: prescribed pressure head value at bottom boundary **
  associate(&
   ! input: model control
   deriv_desired => in_qDrainFlux % deriv_desired, &          ! flag to indicate if derivatives are desired
   ixRichards    => in_qDrainFlux % ixRichards   , &          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   nodeMatricHeadLiq => in_qDrainFlux % nodeMatricHeadLiq, &  ! liquid matric head in the lowest unsaturated node (m)
   nodeVolFracLiq    => in_qDrainFlux % nodeVolFracLiq   , &  ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   nodeDepth  => in_qDrainFlux % nodeDepth , &                ! depth of the lowest unsaturated soil layer (m)
   ! input: diriclet boundary conditions
   lowerBoundHead  => in_qDrainFlux % lowerBoundHead , &      ! lower boundary condition for matric head (m)
   lowerBoundTheta => in_qDrainFlux % lowerBoundTheta, &      ! lower boundary condition for volumetric liquid water content (-)
   ! input: transmittance
   bottomSatHydCond  => in_qDrainFlux % bottomSatHydCond , &  ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   iceImpedeFac      => in_qDrainFlux % iceImpedeFac     , &  ! ice impedence factor in the upper-most soil layer (-)
   ! input: transmittance derivatives
   dHydCond_dTemp   => in_qDrainFlux % dHydCond_dTemp  , &    ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   vGn_alpha       => in_qDrainFlux % vGn_alpha      , &      ! van Genuchten "alpha" parameter (m-1)
   vGn_n           => in_qDrainFlux % vGn_n          , &      ! van Genuchten "n" parameter (-)
   vGn_m           => in_qDrainFlux % vGn_m          , &      ! van Genuchten "m" parameter (-)
   theta_sat       => in_qDrainFlux % theta_sat      , &      ! soil porosity (-)
   theta_res       => in_qDrainFlux % theta_res      , &      ! soil residual volumetric water content (-)
   ! output: hydraulic conductivity at the bottom of the unsaturated zone
   bottomHydCond => out_qDrainFlux % bottomHydCond, &         ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   bottomDiffuse => out_qDrainFlux % bottomDiffuse, &         ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
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
   select case(ixRichards)
     case(moisture) ! moisture-based form of Richards' equation
       ! compute the hydraulic conductivity and diffusivity at the boundary
       bottomHydCond = hydCond_liq(lowerBoundTheta,bottomSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
       bottomDiffuse = dPsi_dTheta(lowerBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * bottomHydCond
       ! compute the capillary flux
       cflux = -bottomDiffuse*(lowerBoundTheta - nodeVolFracLiq) / (nodeDepth*0.5_rkind)
     case(mixdform) ! mixed form of Richards' equation
       ! compute the hydraulic conductivity and diffusivity at the boundary
       bottomHydCond = hydCond_psi(lowerBoundHead,bottomSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
       bottomDiffuse = realMissing
       ! compute the capillary flux
       cflux = -bottomHydCond*(lowerBoundHead  - nodeMatricHeadLiq) / (nodeDepth*0.5_rkind)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select 
   scalarDrainage = cflux + bottomHydCond

   if (deriv_desired) then ! compute derivatives
     ! hydrology derivatives
     select case(ixRichards)  ! select form of Richards' equation
       case(moisture); dq_dHydStateUnsat = bottomDiffuse/(nodeDepth/2._rkind)
       case(mixdform); dq_dHydStateUnsat = bottomHydCond/(nodeDepth/2._rkind)
       case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
     end select
     ! energy derivatives
     dq_dNrgStateUnsat = -(dHydCond_dTemp/2._rkind)*(lowerBoundHead  - nodeMatricHeadLiq)/(nodeDepth*0.5_rkind)&
                       & + dHydCond_dTemp/2._rkind
   else     ! do not desire derivatives
     dq_dHydStateUnsat = realMissing
     dq_dNrgStateUnsat = realMissing
   end if
 
  end associate
 end subroutine update_qDrainFlux_prescribedHead

 subroutine update_qDrainFlux_funcBottomHead
  ! ** Update operations for qDrainFlux: prescribed pressure head function at bottom boundary **
  associate(&
   ! input: model control
   deriv_desired => in_qDrainFlux % deriv_desired, &          ! flag to indicate if derivatives are desired
   ixRichards    => in_qDrainFlux % ixRichards   , &          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   nodeMatricHeadLiq => in_qDrainFlux % nodeMatricHeadLiq, &  ! liquid matric head in the lowest unsaturated node (m)
   nodeVolFracLiq    => in_qDrainFlux % nodeVolFracLiq   , &  ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   nodeHeight => in_qDrainFlux % nodeHeight, &                ! height of the lowest unsaturated soil node (m)
   ! input: derivative in soil water characteristic
   node_dPsi_dTheta => in_qDrainFlux % node_dPsi_dTheta, &    ! derivative of the soil moisture characteristic w.r.t. theta (m)
   ! input: transmittance
   surfaceSatHydCond => in_qDrainFlux % surfaceSatHydCond, &  ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: soil parameters
   vGn_alpha       => in_qDrainFlux % vGn_alpha      , &      ! van Genuchten "alpha" parameter (m-1)
   vGn_n           => in_qDrainFlux % vGn_n          , &      ! van Genuchten "n" parameter (-)
   vGn_m           => in_qDrainFlux % vGn_m          , &      ! van Genuchten "m" parameter (-)
   theta_sat       => in_qDrainFlux % theta_sat      , &      ! soil porosity (-)
   theta_res       => in_qDrainFlux % theta_res      , &      ! soil residual volumetric water content (-)
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
   select case(ixRichards) ! select form of Richards' equation
     case(moisture); nodePsi = matricHead(nodeVolFracLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     case(mixdform); nodePsi = nodeMatricHeadLiq
   end select
   zWater = nodeHeight - nodePsi
   scalarDrainage = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)

   if (deriv_desired) then ! compute derivatives
     ! hydrology derivatives
     select case(ixRichards)  ! select form of Richards' equation
       case(moisture); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * node_dPsi_dTheta*exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
       case(mixdform); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
       case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
     end select
     ! energy derivatives
     err=20; message=trim(message)//"not yet implemented energy derivatives"; return_flag=.true.; return
   else     ! do not desire derivatives
     dq_dHydStateUnsat = realMissing
     dq_dNrgStateUnsat = realMissing
   end if

  end associate
 end subroutine update_qDrainFlux_funcBottomHead

 subroutine update_qDrainFlux_freeDrainage
  ! ** Update operations for qDrainFlux: free drainage at bottom boundary **
  associate(&
   ! input: model control
   deriv_desired => in_qDrainFlux % deriv_desired, &          ! flag to indicate if derivatives are desired
   ixRichards    => in_qDrainFlux % ixRichards   , &          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: transmittance
   nodeHydCond       => in_qDrainFlux % nodeHydCond    , &    ! hydraulic conductivity at the node itself (m s-1)
   ! input: transmittance derivatives
   dHydCond_dVolLiq => in_qDrainFlux % dHydCond_dVolLiq, &    ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
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

   if (deriv_desired) then ! compute derivatives
     ! hydrology derivatives
     select case(ixRichards)  ! select form of Richards' equation
       case(moisture); dq_dHydStateUnsat = dHydCond_dVolLiq*kAnisotropic
       case(mixdform); dq_dHydStateUnsat = dHydCond_dMatric*kAnisotropic
       case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
     end select
     ! energy derivatives
     dq_dNrgStateUnsat = dHydCond_dTemp*kAnisotropic
   else     ! do not desire derivatives
     dq_dHydStateUnsat = realMissing
     dq_dNrgStateUnsat = realMissing
   end if

  end associate
 end subroutine update_qDrainFlux_freeDrainage

 subroutine update_qDrainFlux_zeroFlux
  ! ** Update operations for qDrainFlux: zero flux condition at bottom boundary **
  associate(&
   ! input: model control
   deriv_desired => in_qDrainFlux % deriv_desired, &          ! flag to indicate if derivatives are desired
   ! output: drainage flux from the bottom of the soil profile
   scalarDrainage => out_qDrainFlux % scalarDrainage, &       ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   dq_dHydStateUnsat => out_qDrainFlux % dq_dHydStateUnsat, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   dq_dNrgStateUnsat => out_qDrainFlux % dq_dNrgStateUnsat  & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
  &)

   scalarDrainage = 0._rkind
   if (deriv_desired) then
     dq_dHydStateUnsat = 0._rkind
     dq_dNrgStateUnsat = 0._rkind
   else
     dq_dHydStateUnsat = realMissing
     dq_dNrgStateUnsat = realMissing
   end if

  end associate
 end subroutine update_qDrainFlux_zeroFlux

 subroutine finalize_qDrainFlux
  ! ** Finalize operations for qDrainFlux **
  associate(&
   ! output: error control
   err     => out_qDrainFlux % err    , & ! error code
   message => out_qDrainFlux % message  & ! error message
  &)
   ! final error check
   if (err /= 0_i4b) then
    message=trim(message)//'unanticipated error in qDrainFlux'
    return_flag=.true.; return
   end if
  end associate
 end subroutine finalize_qDrainFlux

end subroutine qDrainFlux

end module soilLiqFlx_module
