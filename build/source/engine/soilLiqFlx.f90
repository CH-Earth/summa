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
USE globalData,only:veryBig                ! a very big number
USE globalData,only:verySmall              ! a small number used as an additive constant to check if substantial difference among real numbers
USE globalData,only:verySmaller            ! a smaller number used as an additive constant to check if substantial difference among real numbers

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
  noInfiltrationExcess          ! no infiltration excess runoff

! -----------------------------------------------------------------------------------------------------------
implicit none
private
public::soilLiqFlx

! flag to denote if updating infiltration during iterations for testing purposes
logical(lgt),parameter :: updateInfil=.true. 
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
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: mLayerTranspireFrac ! fraction of transpiration allocated to each soil layer (-)
  ! diagnostic variables
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: iceImpedeFac        ! ice impedence factor at layer mid-points (-)
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: mLayerDiffuse       ! diffusivity at layer mid-point (m2 s-1)
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: dHydCond_dVolLiq    ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: dDiffuse_dVolLiq    ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
  real(rkind),dimension(in_soilLiqFlx % nSoil)    :: dHydCond_dTemp      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  real(rkind),dimension(0:in_soilLiqFlx % nSoil)  :: iLayerHydCond       ! hydraulic conductivity at layer interface (m s-1)
  real(rkind),dimension(0:in_soilLiqFlx % nSoil)  :: iLayerDiffuse       ! diffusivity at layer interface (m2 s-1)
  ! compute surface flux
  integer(i4b)                                    :: nRoots              ! number of soil layers with roots
  integer(i4b)                                    :: ixIce               ! index of the lowest soil layer that contains ice
  real(rkind),dimension(0:in_soilLiqFlx % nSoil)  :: iLayerHeight        ! height of the layer interfaces (m)
   ! error control
  logical(lgt)                                    :: return_flag         ! flag for return statements
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

  ! get a copy of iLayerHeight (for soil layers only)
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
   if(nRoots==0)then; message=trim(message)//'no layers with roots'; err=20; return_flag=.true.; return; end if
  end associate

  ! ** identify lowest soil layer with ice **
  ! NOTE: cannot use count because there may be an unfrozen wedge
  associate(&
    mLayerVolFracIceTrial => in_soilLiqFlx % mLayerVolFracIceTrial & ! intent(in): volumetric fraction of ice at the current iteration (-)
  &)
   ixIce = 0  ! initialize the index of the ice layer (0 means no ice in the soil profile)
   do iLayer=1,nSoil ! (loop through soil layers)
     if (mLayerVolFracIceTrial(iLayer) > verySmaller) ixIce = iLayer
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
   if(err/=0)then; message=trim(message)//trim("finalize_soilLiqFlx: final error check failed"); return_flag=.true.; return; end if
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
   if (abs(sum(mLayerTranspireFrac) - 1._rkind) > verySmaller) then
     message=trim(message)//'fraction transpiration in soil layers does not sum to one'; err=20; return_flag=.true.; return
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
   if(err/=0)then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
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
   if(err/=0)then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
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
   if(err/=0)then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
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
   if(err/=0)then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
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
   if (scalarVolFracIceTrial > epsilon(iceImpedeFac)) then
     dK_dLiq__noIce   = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)  ! [.true. = analytical]
     dHydCond_dVolLiq = hydCond_noIce*dIceImpede_dLiq + dK_dLiq__noIce*iceImpedeFac
   else
     dHydCond_dVolLiq = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)
   end if
   dPsi_dTheta2a    = dPsi_dTheta2(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m,.true.)   ! [.true. = analytical] compute derivative in dPsi_dTheta (m)
   dDiffuse_dVolLiq = dHydCond_dVolLiq*scalardPsi_dTheta + scalarHydCond*dPsi_dTheta2a
   dHydCond_dMatric = realMissing ! not used, so cause problems

  end associate
 end subroutine update_diagv_node_hydraulic_conductivity_moisture_form

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
! private subroutine surfaceFlx: compute the surface flux and its derivative
! ***************************************************************************************************************
subroutine surfaceFlx(io_soilLiqFlx,in_surfaceFlx,io_surfaceFlx,out_surfaceFlx)
  USE soil_utils_module,only:volFracLiq            ! compute volumetric fraction of liquid water as a function of matric head (-)
  USE soil_utils_module,only:hydCond_psi           ! compute hydraulic conductivity as a function of matric head (m s-1)
  USE soil_utils_module,only:hydCond_liq           ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
  USE soil_utils_module,only:dPsi_dTheta           ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  USE soil_utils_module,only:crit_soilT            ! compute critical temperature below which ice exists
  USE soil_utils_module,only:gammp,gammp_complex   ! compute the regularized lower incomplete Gamma function
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
  ! FUSE
  real(rkind),parameter            :: alpha_LSE=1.e3_rkind                ! smoothness parameter for LSE smoother function
  real(rkind),parameter            :: roundoff_tolerance = 1.e2_rkind * epsilon(1._rkind) ! tolerance for round-off error is near machine epsilon 
  real(rkind)                      :: S1                                  ! total water content in upper soil layer (m)
  real(rkind)                      :: S1_max                              ! Maximum storage in the upper layer (m)
  ! FUSE PRMS variables
  real(rkind)                      :: phi_tens                            ! fraction of total storage as tension storage (m)
  real(rkind)                      :: SatArea_max                              ! maximum saturated area (-)
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
  real(rkind)                      :: x_crit                              ! critical x (random variable) value
  real(rkind)                      :: zeta_crit_n                         ! critical topographic index value (power-transfomred)
  complex(rkind)                   :: lambda_n                            ! mean of the power-transformed topographic index
  real(rkind)                      :: n_topmodel                          ! TOPMODEL exponent exponent (must be sufficiently large to avoid divergence of lambda_n -- n>=3.5 or so)
  ! derivatives
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
  real(rkind) :: dInfilArea_dWat(1:in_surfaceFlx % nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dInfilArea_dTk(1:in_surfaceFlx % nSoil)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind) :: dFrozenArea_dWat(1:in_surfaceFlx % nSoil) ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dFrozenArea_dTk(1:in_surfaceFlx % nSoil)  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
  real(rkind) :: dInfilRate_dWat(1:in_surfaceFlx % nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dInfilRate_dTk(1:in_surfaceFlx % nSoil)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
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
   dVolFracLiq_dWat(:)    = 0._rkind
   dVolFracIce_dWat(:)    = 0._rkind
   dVolFracLiq_dTk(:)     = 0._rkind
   dVolFracIce_dTk(:)     = 0._rkind
   dInfilArea_dWat(:)     = 0._rkind
   dInfilArea_dTk(:)      = 0._rkind
   dxMaxInfilRate_dWat(:) = 0._rkind
   dxMaxInfilRate_dTk(:)  = 0._rkind
   dFrozenArea_dWat(:)    = 0._rkind
   dFrozenArea_dTk(:)     = 0._rkind
   dq_dHydStateVec(:)     = 0._rkind
   dq_dNrgStateVec(:)     = 0._rkind ! energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
  end associate

  ! initialize runoff and infiltration values
  associate(&
   scalarSurfaceRunoff       => out_surfaceFlx % scalarSurfaceRunoff       , & ! surface runoff (m s-1)
   scalarSurfaceRunoff_IE    => out_surfaceFlx % scalarSurfaceRunoff_IE    , & ! infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE    => out_surfaceFlx % scalarSurfaceRunoff_SE    , & ! saturation excess surface runoff (m s-1)
   scalarSurfaceInfiltration => out_surfaceFlx % scalarSurfaceInfiltration   & ! surface infiltration (m s-1)
  &)
   scalarSurfaceRunoff       = 0._rkind 
   scalarSurfaceRunoff_IE    = 0._rkind  
   scalarSurfaceRunoff_SE    = 0._rkind  
   scalarSurfaceInfiltration = 0._rkind 
  end associate

 end subroutine initialize_surfaceFlx

 subroutine update_surfaceFlx
  ! **** Update operations for surfaceFlx ****
  associate(&
   ! input: model control
   firstSplitOper => in_surfaceFlx % firstSplitOper, & ! flag indicating if desire to compute infiltration
   bc_upper   => in_surfaceFlx % bc_upper,           & ! index defining the type of boundary conditions
   ixInfRateMax => in_surfaceFlx % ixInfRateMax,     & ! index defining the maximum infiltration rate method
   surfRun_SE => in_surfaceFlx % surfRun_SE,         & ! index defining the saturation excess surface runoff method
   ! input to compute infiltration
   scalarRainPlusMelt => in_surfaceFlx % scalarRainPlusMelt, & ! rain plus melt  (m s-1)
   ! output: infiltration area and saturated area
   scalarInfilArea    => io_surfaceFlx % scalarInfilArea,      & ! fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea => io_surfaceFlx % scalarSaturatedArea, & ! saturated area fraction (-)
   ! output: runoff and infiltration 
   scalarSurfaceRunoff_SE    => out_surfaceFlx % scalarSurfaceRunoff_SE,    & ! saturation excess surface runoff (m s-1)
   scalarSurfaceRunoff       => out_surfaceFlx % scalarSurfaceRunoff,       & ! surface runoff (m s-1)
   ! output: derivatives in surface infiltration w.r.t. ...
   dq_dHydStateVec => out_surfaceFlx % dq_dHydStateVec, & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   dq_dNrgStateVec => out_surfaceFlx % dq_dNrgStateVec, & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
   ! output: error control
   err      => out_surfaceFlx % err    , & ! error code
   message  => out_surfaceFlx % message  & ! error message
  &)
 
   ! compute the surface flux and its derivative
   if (firstSplitOper .or. updateInfil) then
     select case(bc_upper)
       case(prescribedHead) ! head condition, no frozen area and all area infiltrates
         call update_surfaceFlx_prescribedHead; if (return_flag) return 
 
       case(liquidFlux)     ! flux condition
         ! compute volumetric fraction of liquid and ice water in each soil layer and their derivatives
         if(updateInfil) call update_volFracLiq_derivatives; if (return_flag) return

         ! Get infiltration area not considering frozen area, based on SE method
         select case(surfRun_SE) ! saturation excess surface runoff method, sets infiltration area (not considering frozen) and its derivatives
           case(zero_SE)         ! zero saturation excess surface runoff, all area infiltrates if not frozen
            io_surfaceFlx % scalarInfilArea = 1._rkind 
           case(homegrown_SE)    ! homegrown saturation excess surface runoff (original SUMMA method)
              call update_surfaceFlx_homegrown_infilArea;     if (return_flag) return
           case(FUSEPRMS)        ! FUSE PRMS surface runoff
             call update_surfaceFlx_FUSE_PRMS_infilArea;      if (return_flag) return
           case(FUSEAVIC)        ! FUSE ARNO/VIC surface runoff
             call update_surfaceFlx_FUSE_ARNO_VIC_infilArea;  if (return_flag) return
           case(FUSETOPM)        ! FUSE TOPMODEL surface runoff
             call update_surfaceFlx_FUSE_TOPMODEL_infilArea;  if (return_flag) return
           case default; err=20; message=trim(message)//'unknown saturation excess surface runoff method'; return_flag=.true.; return
         end select

         ! compute saturated area and saturation excess surface runoff
         if (scalarInfilArea < 0._rkind) then; err=20; message=trim(message)//'infiltration area less than zero'; return_flag=.true.; return; end if
         scalarSaturatedArea = 1._rkind - scalarInfilArea
         scalarSurfaceRunoff_SE = scalarRainPlusMelt * scalarSaturatedArea

         ! Calculate maximum infiltration rate and scalarFrozenArea (and their derivatives if needed)
         select case(ixInfRateMax)       ! maximum infiltration rate method (controls infiltration excess surface runoff)
           case(noInfiltrationExcess)    ! zero infiltration excess surface runoff
             call update_surfaceFlx_liquidFlux_noinfratemax
           case(GreenAmpt, topmodel_GA)  ! infiltration excess runoff possible
             call update_surfaceFlx_liquidFlux_calculate_infratemax;  if (return_flag) return
           case default; err=20; message=trim(message)//'unknown infiltration excess surface runoff method'; return_flag=.true.; return
         end select

         ! Compute total infiltration, gets infiltration excess surface runoff, modifies saturation excess surface runoff if more rain than can infiltrate
         call update_surfaceFlx_liquidFlux_infiltration;  if (return_flag) return
         
         ! update the derivatives for any combination of SE and IE parametrization options 
         if(updateInfil) call update_surfaceFlx_liquidFlux_derivatives

       case default; err=20; message=trim(message)//'unknown upper boundary condition for soil hydrology'; return_flag=.true.; return ! end of select of bc_upper
     end select 
   else ! do not compute infiltration after first flux call in a splitting operation unless updateInfil is true
     dq_dHydStateVec(:) = realMissing ! not used, so cause problems
     dq_dNrgStateVec(:) = realMissing ! not used, so cause problems
   end if 

  end associate
 end subroutine update_surfaceFlx

subroutine update_volFracLiq_derivatives
  ! **** Updates the derivatives for volumetric fraction of liquid and ice water in each soil layer ****
  ! local variables
  integer(i4b)      :: nLayers         ! number of soil layers to process
  logical(lgt)      :: doIce           ! flag indicating whether ice derivatives are needed

  associate(&
   ! input: model control
   ixInfRateMax        => in_surfaceFlx % ixInfRateMax       , & ! index defining the maximum infiltration rate method
   surfRun_SE          => in_surfaceFlx % surfRun_SE         , & ! index defining the saturation excess surface runoff method
   ixRichards          => in_surfaceFlx % ixRichards         , & ! index defining the option for Richards' equation (moisture or mixdform)
   nRoots              => in_surfaceFlx % nRoots             , & ! number of layers that contain roots
   nSoil               => in_surfaceFlx % nSoil              , & ! total number of soil layers
   ! input: state and diagnostic variables
   mLayerTemp          => in_surfaceFlx % mLayerTemp         , & ! temperature (K)
   mLayerMatricHead    => in_surfaceFlx % mLayerMatricHead   , & ! matric head in each soil layer (m)
   ! input: pre-computed derivatives in ...
   dTheta_dTk          => in_surfaceFlx % dTheta_dTk         , & ! ... volumetric liquid water content w.r.t. temperature (K-1)
   dTheta_dPsi         => in_surfaceFlx % dTheta_dPsi        , & ! ... the soil water characteristic w.r.t. psi (m-1)
   mLayerdPsi_dTheta   => in_surfaceFlx % mLayerdPsi_dTheta  , & ! ... the soil water characteristic w.r.t. theta (m)
   ! output: error control
   err      => out_surfaceFlx % err    , & ! error code
   message  => out_surfaceFlx % message  & ! error message
  &)

   ! determine number of layers to process and whether ice derivatives are needed
   if (surfRun_SE ==homegrown_SE) then ! need only root zone derivatives but need ice derivatives
     nLayers = nRoots
     doIce = .true.
   else ! might need entire soil column (FUSE methods), might need ice derivatives (infiltration excess method)
     if (ixInfRateMax == noInfiltrationExcess) then
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
     end if ! (if ixInfRateMax)
   end if ! (if homegrown_SE)

   if (nLayers > 0) then
     select case(ixRichards)  ! form of Richards' equation
       case(moisture)
         dVolFracLiq_dWat(:) = 1._rkind
         dVolFracIce_dWat(:) = mLayerdPsi_dTheta(:) - 1._rkind
       case(mixdform)
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
       case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
     end select 
     dVolFracLiq_dTk(:) = dTheta_dTk(:) !already zeroed out if not below critical temperature
     if(doIce) dVolFracIce_dTk(:) = -dVolFracLiq_dTk(:) !often can and will simplify one of these terms out
   end if

  end associate
 end subroutine update_volFracLiq_derivatives

 subroutine update_surfaceFlx_liquidFlux_derivatives
  ! **** Updates the derivatives for the liquid flux upper hydrology boundary condition if needed ****
  associate(&
   ! input: flux at the upper boundary
   scalarRainPlusMelt => in_surfaceFlx % scalarRainPlusMelt , & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate    => io_surfaceFlx % xMaxInfilRate    , & ! maximum infiltration rate (m s-1)
   scalarInfilArea  => io_surfaceFlx % scalarInfilArea  , & ! fraction of area where water can infiltrate, may be frozen (-)
   scalarFrozenArea => io_surfaceFlx % scalarFrozenArea , & ! fraction of area that is considered impermeable due to soil ice (-)
   ! output: derivatives in surface infiltration w.r.t. ...
   dq_dHydStateVec => out_surfaceFlx % dq_dHydStateVec, & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   dq_dNrgStateVec => out_surfaceFlx % dq_dNrgStateVec, & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
   ! output: error control
   err      => out_surfaceFlx % err    , & ! error code
   message  => out_surfaceFlx % message  & ! error message
  &)

   ! Compute total runoff derivatives, do w.r.t. infiltration only, scalarRainPlusMelt accounted for in computJacob* module
   ! Do not need to break into IE and SE components since they are never used separately in the Jacobian assembly
   dq_dHydStateVec(:) = (1._rkind - scalarFrozenArea)&
                       * ( dInfilArea_dWat(:)*min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dWat(:) )&
                       + (-dFrozenArea_dWat(:))*scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
   ! energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
   dq_dNrgStateVec(:) = (1._rkind - scalarFrozenArea)&
                       * ( dInfilArea_dTk(:) *min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dTk(:)  )&
                       + (-dFrozenArea_dTk(:)) *scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)

  end associate
 end subroutine update_surfaceFlx_liquidFlux_derivatives

 subroutine update_surfaceFlx_FUSE_PRMS_infilArea
  ! **** Update operations for surfaceFlx: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- PRMS ****
  use soil_utils_module,only:LogSumExp  ! smooth max/min
  use soil_utils_module,only:SoftArgMax ! smooth arg max/min (for derivatives of LogSumExp)

  ! local variables
  real(rkind)           :: dS1_dLiq(1:in_surfaceFlx % nSoil)        ! derivative of S1 w.r.t. liquid water content
  real(rkind)           :: S1_T_derivatives(1:2)                    ! array of derivatives for S1_T
  real(rkind)           :: dS1_T_dS1                                ! derivative of S1_T w.r.t S1
  real(rkind)           :: dS1_T_dLiq(1:in_surfaceFlx % nSoil)      ! derivative of S1_T w.r.t liquid water content

  associate(&
   nSoil            => in_surfaceFlx % nSoil,            & ! number of soil layers
   mLayerVolFracLiq => in_surfaceFlx % mLayerVolFracLiq, & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth      => in_surfaceFlx % mLayerDepth,      & ! depth of soil layers (m) 
   iLayerHeight     => in_surfaceFlx % iLayerHeight,     & ! height at the interface of each layer for soil layers only (m)
   theta_sat        => in_surfaceFlx % theta_sat,        & ! soil porosity (-)
   ! output: error control
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)

   ! validation of parameters
   SatArea_max = in_surfaceFlx % FUSE_Ac_max
   phi_tens    = in_surfaceFlx % FUSE_phi_tens
   ! validate input parameters 
   if ((SatArea_max<0._rkind).or.(SatArea_max>1._rkind)) then
    err=10; message=trim(message)//"FUSE PRMS surface runoff error: invalid SatArea_max (max saturated area) value"; return_flag=.true.; return
   end if
   if ((phi_tens<0._rkind).or.(phi_tens>1._rkind)) then
    err=10; message=trim(message)//"FUSE PRMS surface runoff error: invalid phi_tens (tension storage fraction) value"; return_flag=.true.; return
   end if

  ! compute water content in upper FUSE layer
   S1     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in upper FUSE layer (m)
   if (S1 <= 0._rkind) then; io_surfaceFlx % scalarInfilArea = 1._rkind; return; end if ! if no water, unsaturated and all area infiltrates
   S1_max = iLayerHeight(nSoil) * theta_sat             ! max water storage for upper FUSE layer (m)

  ! compute tension water content
   S1_T_max = phi_tens * S1_max
   S1_T     = LogSumExp(-alpha_LSE,[S1,S1_T_max],err) ! smooth approximation to S1_T=min(S1,S1_T_max)
   if(err/=0)then; err=10; message=trim(message)//"FUSE PRMS surface runoff: error in LogSumExp"; return_flag=.true.; return; end if
   if (S1_T < 0._rkind) then ! check for errors
    err=10; message=trim(message)//"FUSE PRMS surface runoff: S1_T is negative (may need to adjust magnitude of alpha_LSE)"; return_flag=.true.; return
   end if

   ! define the infiltrating area and derivatives for the non-frozen part of the cell/basin
   io_surfaceFlx % scalarInfilArea = 1._rkind - (S1_T/S1_T_max)*SatArea_max
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

 end subroutine update_surfaceFlx_FUSE_PRMS_infilArea

 subroutine update_surfaceFlx_FUSE_ARNO_VIC_infilArea
  ! **** Update operations for surfaceFlx: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- ARNO/VIC ****
  use soil_utils_module,only:LogSumExp  ! smooth max/min
  use soil_utils_module,only:SoftArgMax ! smooth arg max/min (for derivatives of LogSumExp)

  ! local variables
  real(rkind)            :: dS1_dLiq(1:in_surfaceFlx % nSoil)      ! derivative of S1 w.r.t. liquid water content
  real(rkind)            :: dS1_star_dS1                           ! derivative in S1_star w.r.t S1
  real(rkind)            :: dbase_dS1                              ! derivative of base w.r.t S1
  real(rkind)            :: S1_star_derivatives(1:2)               ! array of derivatives for S1_star from SoftArgMax function

  associate(&
   nSoil            => in_surfaceFlx % nSoil,            & ! number of soil layers
   mLayerVolFracLiq => in_surfaceFlx % mLayerVolFracLiq, & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth      => in_surfaceFlx % mLayerDepth,      & ! depth of soil layers (m) 
   iLayerHeight     => in_surfaceFlx % iLayerHeight,     & ! height at the interface of each layer for soil layers only (m)
   theta_sat        => in_surfaceFlx % theta_sat,        & ! soil porosity (-)
   ! output: error control
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)

   ! validation of input parameters
   b_arnovic = in_surfaceFlx % FUSE_b ! interface ARNO/VIC exponent
   if ((b_arnovic < 0.001_rkind).or.(b_arnovic > 3._rkind)) then
    err=10; message=trim(message)//"FUSE ARNO/VIC exponent must be between 0.001 and 3"; return_flag=.true.; return
   end if

   ! compute water content in FUSE layers
   S1     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in FUSE layers (m)
   if (S1 <= 0._rkind) then; io_surfaceFlx % scalarInfilArea = 1._rkind; return; end if ! if no water, unsaturated and all area infiltrates
   S1_max = iLayerHeight(nSoil) * theta_sat             ! max water storage for FUSE layers (m)

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
   io_surfaceFlx % scalarInfilArea = base**b_arnovic
   
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

 end subroutine update_surfaceFlx_FUSE_ARNO_VIC_infilArea


 subroutine update_surfaceFlx_FUSE_TOPMODEL_infilArea
  ! **** Update operations for surfaceFlx: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- TOPMODEL ****
  ! local variables
  real(rkind)                      :: lambda                            ! mean
  real(rkind)                      :: chi_topmodel                      ! scale
  real(rkind)                      :: mu                                ! offset
  real(rkind),parameter            :: zeta_upper=1.e3_rkind             ! upper limit of integral (approaches infinity, but ~1000 provides an accurate result) 
  real(rkind)                      :: zeta_crit                         ! critical topographic index value (log space)
  complex(rkind)                   :: F1,F2                             ! temporary storage for regularized lower incomplete gamma function values
  ! derivitive values
  real(rkind)                      :: dS1_dLiq(1:in_surfaceFlx % nSoil) ! derivative in S1 w.r.t liquid water content 
  real(rkind)                      :: dzeta_crit_n_dS1                  ! derivative of zeta_crit_n w.r.t S1
  real(rkind)                      :: dzeta_crit_dzeta_crit_n           ! derivative of zeta_crit w.r.t zeta_crit_n
  real(rkind)                      :: dx_crit_dzeta_crit                ! derivative of x_crit w.r.t zeta_crit
  real(rkind)                      :: dx_crit_dS1                       ! derivative of x_crit w.r.t S1
  real(rkind)                      :: dgammp_dx_crit                    ! derivative of gammp function in SatArea w.r.t x_crit

  associate(&
   nSoil            => in_surfaceFlx % nSoil,            & ! number of soil layers
   mLayerVolFracLiq => in_surfaceFlx % mLayerVolFracLiq, & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth      => in_surfaceFlx % mLayerDepth,      & ! depth of soil layers (m) 
   iLayerHeight     => in_surfaceFlx % iLayerHeight,     & ! height at the interface of each layer for soil layers only (m)
   theta_sat        => in_surfaceFlx % theta_sat,        & ! soil porosity (-)
   ! output: error control
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)

   ! interface FUSE input parameters
   lambda       = in_surfaceFlx % FUSE_lambda
   chi_topmodel = in_surfaceFlx % FUSE_chi
   mu           = in_surfaceFlx % FUSE_mu
   n_topmodel   = in_surfaceFlx % FUSE_n

   ! compute water content in lower FUSE layer, here the entire soil column is used
   S1     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in lower FUSE layer (m)
   if (S1 <= 0._rkind) then; io_surfaceFlx % scalarInfilArea = 1._rkind; return; end if ! if no water, unsaturated and all area infiltrates
   S1_max = iLayerHeight(nSoil) * theta_sat             ! max water storage for lower FUSE layer (m)

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
   io_surfaceFlx % scalarInfilArea = gammp(alpha_topmodel,x_crit/chi_topmodel)

  else ! if (S1 == 0) no water is stored in lower FUSE layer (based on asymptotic behaviour of integral in eq. 9c of Clark et al. (2008))
   io_surfaceFlx % scalarInfilArea = 1._rkind
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

 end subroutine update_surfaceFlx_FUSE_TOPMODEL_infilArea

 subroutine update_surfaceFlx_prescribedHead
  ! **** Update operations for surfaceFlx: prescribed pressure head condition ****
  associate(&
   ! input: model control
   ixRichards     => in_surfaceFlx % ixRichards     , & ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   scalarMatricHeadLiq => in_surfaceFlx % scalarMatricHeadLiq , & ! liquid matric head in the upper-most soil layer (m)
   scalarVolFracLiq    => in_surfaceFlx % scalarVolFracLiq    , & ! volumetric liquid water content in the upper-most soil layer (-)
   ! input: depth of each soil layer (m)
   mLayerDepth  => in_surfaceFlx % mLayerDepth  , & ! depth of each soil layer (m)
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
   surfaceDiffuse => io_surfaceFlx % surfaceDiffuse , & ! hydraulic diffusivity at the surface (m2 s-1)
   ! output: runoff and infiltration
   scalarSurfaceRunoff       => out_surfaceFlx % scalarSurfaceRunoff       , & ! surface runoff (m s-1)
   scalarSurfaceRunoff_IE    => out_surfaceFlx % scalarSurfaceRunoff_IE    , & ! infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE    => out_surfaceFlx % scalarSurfaceRunoff_SE    , & ! saturation excess surface runoff (m s-1)
   scalarSurfaceInfiltration => out_surfaceFlx % scalarSurfaceInfiltration , & ! surface infiltration (m s-1)
   ! output: derivatives in surface infiltration w.r.t. ...
   scalarSoilControl  => io_surfaceFlx % scalarSoilControl    , & ! soil control on infiltration for derivative
   dq_dHydStateVec    => out_surfaceFlx % dq_dHydStateVec     , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   dq_dNrgStateVec    => out_surfaceFlx % dq_dNrgStateVec     , & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
   ! output: error control
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)

   ! surface runoff iz zero for the head condition
   scalarSurfaceRunoff_IE = 0._rkind ! infiltration excess runoff 
   scalarSurfaceRunoff_SE = 0._rkind ! saturation excess runoff 
   scalarSurfaceRunoff    = 0._rkind 

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
   scalarSoilControl = 0._rkind 

   ! compute the derivatives at the surface, only has a non-zero value for the upper-most soil layer
   if(updateInfil)then
     select case(ixRichards)  ! select form of Richards' equation
       case(moisture); dq_dHydStateVec(1) = -surfaceDiffuse/(mLayerDepth(1)/2._rkind)
       case(mixdform); dq_dHydStateVec(1) = -surfaceHydCond/(mLayerDepth(1)/2._rkind)
       case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
     end select
     ! note: energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
     dq_dNrgStateVec(1) = -(dHydCond_dTemp/2._rkind)*(scalarMatricHeadLiq - upperBoundHead)/(mLayerDepth(1)*0.5_rkind) + dHydCond_dTemp/2._rkind
   end if

   ! * additional assignment statements for surfaceFlx input-output object based on presribed head values *
   ! the infiltration is always constrained by the prescribed head so the maximum infiltration rate is set to missing
   io_surfaceFlx % xMaxInfilRate    = realMissing ! maximum infiltration rate (m s-1)
   ! no soil ice assumed for prescribed head condition
   io_surfaceFlx % scalarFrozenArea = 0._rkind      ! fraction of area that is considered impermeable due to soil ice (-)
   ! all area is available for infiltration, and to complement this saturated area (i.e., part where saturation excess runoff occurs) is set to zero
   io_surfaceFlx % scalarInfilArea     = 1._rkind ! fraction of area where water can infiltrate, may be frozen (-)
   io_surfaceFlx % scalarSaturatedArea = 0._rkind ! fraction of area that is considered saturated (-)

  end associate
 end subroutine update_surfaceFlx_prescribedHead

 subroutine update_surfaceFlx_homegrown_infilArea
  ! **** Update operations for surfaceFlx: homegrown saturation excess runoff condition ****
  call update_surfaceFlx_liquidFlux_computation_root_layers 
  call update_surfaceFlx_liquidFlux_computation_available_capacity; if (return_flag) return 
  call update_surfaceFlx_liquidFlux_computation_homegrown  ! this calculates infiltration area ignoring if frozen or not, depends on available capacity (depends on ice and root zone)
 end subroutine update_surfaceFlx_homegrown_infilArea

 subroutine update_surfaceFlx_liquidFlux_noinfratemax
  ! **** Update operations for surfaceFlx: no infiltration excess****
  associate(&
   ! input: model control
   surfRun_SE => in_surfaceFlx % surfRun_SE         & ! index defining the saturation excess surface runoff method
  &)
   io_surfaceFlx % xMaxInfilRate = veryBig ! set to a very large number so rainPlusMelt never exceeds this
   if (surfRun_SE /= homegrown_SE) then  ! frozen area (depends on ice and root zone)
    call update_surfaceFlx_liquidFlux_computation_root_layers
   end if
  end associate
  ! -- main computations - these always need to run
  call update_surfaceFlx_liquidFlux_computation_frozen_area
 end subroutine update_surfaceFlx_liquidFlux_noinfratemax

 subroutine update_surfaceFlx_liquidFlux_calculate_infratemax
  ! **** Update operations for surfaceFlx: infiltration excess possible - calculate max infiltration rate ****
  associate(&
   ! input: model control
   surfRun_SE => in_surfaceFlx % surfRun_SE         & ! index defining the saturation excess surface runoff method
  &)
   if (surfRun_SE /= homegrown_SE) then  ! infiltration rate max depends on available capacity (depends on ice and root zone) and frozen area (depends on ice and root zone)
     call update_surfaceFlx_liquidFlux_computation_root_layers 
     call update_surfaceFlx_liquidFlux_computation_available_capacity; if (return_flag) return
   end if
  end associate
  ! -- main computations - these always need to run
  call update_surfaceFlx_liquidFlux_computation_frozen_area
  call update_surfaceFlx_liquidFlux_computation_max_infiltration_rate
 end subroutine update_surfaceFlx_liquidFlux_calculate_infratemax

 subroutine update_surfaceFlx_liquidFlux_computation_root_layers 
  ! **** Update operations for surfaceFlx: root layer water computation ****
  associate(&
   ! input: model control
   nRoots              => in_surfaceFlx % nRoots            , & ! number of soil layers with roots (-)
   ! input: state and diagnostic variables
   mLayerVolFracLiq    => in_surfaceFlx % mLayerVolFracLiq  , & ! volumetric liquid water content in each soil layer (-)
   mLayerVolFracIce    => in_surfaceFlx % mLayerVolFracIce  , & ! volumetric ice content in each soil layer (-)
   ! input: depth of soil layers (m)
   mLayerDepth         => in_surfaceFlx % mLayerDepth      , & ! depth of each soil layer (m)
   iLayerHeight        => in_surfaceFlx % iLayerHeight     , & ! height at the interface of each layer for soil layers only (m)
   rootingDepth        => in_surfaceFlx % rootingDepth       & ! rooting depth (m)
  &)
 
   ! define the storage in the root zone (m) and derivatives, first initialize
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
       if(updateInfil)then
         dRootZoneLiq_dWat(iLayer) = dVolFracLiq_dWat(iLayer)*mLayerDepth(iLayer)
         dRootZoneIce_dWat(iLayer) = dVolFracIce_dWat(iLayer)*mLayerDepth(iLayer)
         dRootZoneLiq_dTk(iLayer)  = dVolFracLiq_dTk(iLayer) *mLayerDepth(iLayer)
         dRootZoneIce_dTk(iLayer)  = dVolFracIce_dTk(iLayer) *mLayerDepth(iLayer)
       end if
     end do
   end if
   ! process layers where the roots end in the current layer
   rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(nRoots)*min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
   rootZoneIce = rootZoneIce + mLayerVolFracIce(nRoots)*min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
   if(updateInfil)then
     dRootZoneLiq_dWat(nRoots) = dVolFracLiq_dWat(nRoots)*min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneIce_dWat(nRoots) = dVolFracIce_dWat(nRoots)*min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneLiq_dTk(nRoots)  = dVolFracLiq_dTk(nRoots)* min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneIce_dTk(nRoots)  = dVolFracIce_dTk(nRoots)* min(mLayerDepth(nRoots),rootingDepth - iLayerHeight(nRoots-1))
   endif

  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_root_layers 

 subroutine update_surfaceFlx_liquidFlux_computation_available_capacity 
  ! **** Update operations for surfaceFlx: compute and check available capacity to hold water ****
  associate(&
   ! input: depth of each soil layer (m)
   mLayerDepth  => in_surfaceFlx % mLayerDepth  , & ! depth of each soil layer (m)
   ! input: soil parameters
   theta_sat           => in_surfaceFlx % theta_sat   , & ! soil porosity (-)
   rootingDepth        => in_surfaceFlx % rootingDepth, & ! rooting depth (m)
   ! output: error control
   err     => out_surfaceFlx % err    , & ! error code
   message => out_surfaceFlx % message  & ! error message
  &)
   availCapacity = theta_sat*rootingDepth - rootZoneIce
   if (rootZoneLiq > availCapacity+verySmaller) then
     err=20; message=trim(message)//'liquid water in the root zone exceeds capacity'; return_flag=.true.; return
   end if

  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_available_capacity 

 subroutine update_surfaceFlx_liquidFlux_computation_max_infiltration_rate
  ! **** Update operations for surfaceFlx: max infiltration rate and derivatives ****
  associate(&
   ! input: model control
   ixInfRateMax => in_surfaceFlx % ixInfRateMax , & ! index defining the maximum infiltration rate method (GreenAmpt, topmodel_GA, noInfiltrationExcess)
   ! input: transmittance
   surfaceSatHydCond => in_surfaceFlx % surfaceSatHydCond , & ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: soil parameters
   zScale_TOPMODEL     => in_surfaceFlx % zScale_TOPMODEL     , & ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
   rootingDepth        => in_surfaceFlx % rootingDepth        , & ! rooting depth (m)
   wettingFrontSuction => in_surfaceFlx % wettingFrontSuction , & ! Green-Ampt wetting front suction (m)
   ! input-output: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate    => io_surfaceFlx % xMaxInfilRate  & ! maximum infiltration rate (m s-1)
  &)
   ! compute total soil depth from layer thicknesses
   total_soil_depth = sum(in_surfaceFlx % mLayerDepth)

   ! define the depth to the wetting front (m) and derivatives
   depthWettingFront = (rootZoneLiq/availCapacity)*min(rootingDepth,total_soil_depth)
   if(updateInfil)then
     dDepthWettingFront_dWat(:)=( dRootZoneLiq_dWat(:)*min(rootingDepth,total_soil_depth) + dRootZoneIce_dWat(:)*depthWettingFront )/availCapacity
     dDepthWettingFront_dTk(:) =( dRootZoneLiq_dTk(:) *min(rootingDepth,total_soil_depth) + dRootZoneIce_dTk(:)*depthWettingFront  )/availCapacity
    end if

   ! process hydraulic conductivity-controlled infiltration rate
   select case(ixInfRateMax)  ! maximum infiltration rate parameterization (noInfExcess set in update_surfaceFlx)
    case(topmodel_GA)
     ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
     hydCondWettingFront = surfaceSatHydCond * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 1._rkind) )
     ! define the maximum infiltration rate (m s-1)
     xMaxInfilRate = hydCondWettingFront*( (wettingFrontSuction + depthWettingFront)/depthWettingFront )  ! maximum infiltration rate (m s-1)
     ! define the derivatives
     if(updateInfil)then
       fPart1    = hydCondWettingFront
       fPart2    = (wettingFrontSuction + depthWettingFront)/depthWettingFront
       dPart1(:) = surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dWat(:))/total_soil_depth
       dPart2(:) = -dDepthWettingFront_dWat(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
       dxMaxInfilRate_dWat(:) = fPart1*dPart2(:) + fPart2*dPart1(:)
       dPart1(:) = surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dTk(:))/total_soil_depth
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
 end subroutine update_surfaceFlx_liquidFlux_computation_max_infiltration_rate

 subroutine update_surfaceFlx_liquidFlux_computation_homegrown
  ! **** Update operations for surfaceFlx: infiltrating area (ignoring frozen area) for homegrown saturation excess condition ****
  associate(&
   ! input: model control
   nSoil            => in_surfaceFlx % nSoil               , & ! number of soil layers
   nRoots           => in_surfaceFlx % nRoots              , & ! number of layers that contain roots
   ixIce            => in_surfaceFlx % ixIce               , & ! index of lowest ice layer
   mLayerVolFracLiq => in_surfaceFlx % mLayerVolFracLiq    , & ! volumetric liquid water content in each soil layer (-)
   mLayerDepth      => in_surfaceFlx % mLayerDepth         , & ! depth of each soil layer (m)
   ! input: soil parameters
   theta_sat        => in_surfaceFlx % theta_sat           , & ! soil porosity (-)
   qSurfScale       => in_surfaceFlx % qSurfScale          , & ! scaling factor in the surface runoff parameterization (-)
   ! input-output: surface runoff and infiltration flux (m s-1)
   scalarInfilArea  => io_surfaceFlx % scalarInfilArea        & ! fraction of area where water can infiltrate, may be frozen (-)
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

   ! check to ensure we are not infiltrating into a fully saturated column
   if (ixIce<nRoots) then
     if (sum(mLayerVolFracLiq(ixIce+1:nRoots)*mLayerDepth(ixIce+1:nRoots)) > 0.9999_rkind*theta_sat*sum(mLayerDepth(ixIce+1:nRoots))) then 
      scalarInfilArea    = 0._rkind
      dInfilArea_dWat(:) = 0._rkind
      dInfilArea_dTk(:)  = 0._rkind
     end if
   end if
  end associate
 end subroutine update_surfaceFlx_liquidFlux_computation_homegrown

 subroutine update_surfaceFlx_liquidFlux_computation_frozen_area
  ! **** Update operations for surfaceFlx: get impermeable area due to soil freezing ****
  associate(&
   ! input: soil parameters
   soilIceScale        => in_surfaceFlx % soilIceScale       , & ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   soilIceCV           => in_surfaceFlx % soilIceCV          , & ! soil ice CV in Gamma distribution used to define frozen area (-)
   ! output: frozen area
   scalarFrozenArea    => io_surfaceFlx % scalarFrozenArea     & ! fraction of area that is considered impermeable due to soil ice (-)
  &)
   ! define the impermeable area and derivatives due to frozen ground
   if (rootZoneIce > tiny(rootZoneIce)) then  ! (avoid divide by zero)
      alpha = 1._rkind/(soilIceCV**2_i4b)     ! shape parameter in the Gamma distribution
      xLimg = alpha*soilIceScale/rootZoneIce  ! upper limit of the integral
     !if we use this, we will have a derivative of scalarFrozenArea w.r.t. water and temperature in each layer (through mLayerVolFracIce)
     ! Should fix to deal with frozen area in the root zone, calculations would be expensive
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
 end subroutine update_surfaceFlx_liquidFlux_computation_frozen_area

 subroutine update_surfaceFlx_liquidFlux_infiltration
  ! **** Update operations for surfaceFlx: final infiltration and runoff calculations ****
  ! local variables
  real(rkind) :: surfaceInfiltration      ! surface infiltration
  real(rkind) :: surfaceRunoff            ! surface runoff 
  real(rkind) :: surfaceRunoff_IE         ! infiltration excess component of surface runoff 
  real(rkind) :: surfaceRunoff_SE         ! saturation excess component of surface runoff
  real(rkind) :: scalarInfilArea_unfrozen ! infiltration area that is not frozen

  ! compute infiltration and runoff
  associate(&
   ! input: flux at the upper boundary
   scalarRainPlusMelt   => in_surfaceFlx % scalarRainPlusMelt, & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input-output: surface runoff and infiltration flux (m s-1)
   xMaxInfilRate       => io_surfaceFlx % xMaxInfilRate    ,   & ! maximum infiltration rate (m s-1)
   scalarInfilArea     => io_surfaceFlx % scalarInfilArea  ,   & ! fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea => io_surfaceFlx % scalarSaturatedArea, & ! fraction of area that is saturated (-)
   scalarSoilControl   => io_surfaceFlx % scalarSoilControl,   & ! soil control on infiltration for derivative
   scalarFrozenArea    => io_surfaceFlx % scalarFrozenArea,    & ! fraction of area that is considered impermeable due to soil ice (-)
   ! output: runoff and infiltration 
   scalarSurfaceRunoff_IE    => out_surfaceFlx % scalarSurfaceRunoff_IE,    & ! infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE    => out_surfaceFlx % scalarSurfaceRunoff_SE,    & ! saturation excess surface runoff (m s-1)
   scalarSurfaceRunoff       => out_surfaceFlx % scalarSurfaceRunoff,       & ! surface runoff (m s-1)
   scalarSurfaceInfiltration => out_surfaceFlx % scalarSurfaceInfiltration  & ! surface infiltration (m s-1)
  &)
   ! unfrozen infiltration area
   scalarInfilArea_unfrozen=(1._rkind - scalarFrozenArea)*scalarInfilArea
   ! soil control on infiltration for derivative if dependent on scalarRainPlusMelt (needed to compute scalarRainPlusMelt derivative inside computJacob*)
   scalarSoilControl = 0._rkind
   if (updateInfil .and. xMaxInfilRate > scalarRainPlusMelt) then
     scalarSoilControl = scalarInfilArea_unfrozen
   end if

   ! infiltration rate derivatives, will stay at zero if no infiltration excess or if infiltration not being updated
   if(updateInfil)then
     if (xMaxInfilRate < scalarRainPlusMelt) then ! = dxMaxInfilRate_d, dependent on layers not at surface
       dInfilRate_dWat(:) = dxMaxInfilRate_dWat(:)
       dInfilRate_dTk(:)  = dxMaxInfilRate_dTk(:)
     end if
   end if

   ! compute infiltration (m s-1)
   scalarSurfaceInfiltration = scalarInfilArea_unfrozen * min(scalarRainPlusMelt,xMaxInfilRate)  ! rain+melt falling on unfrozen unsaturated area, and smaller than max infiltration rate
 
   ! compute surface runoff (m s-1)
   scalarSurfaceRunoff = scalarRainPlusMelt - scalarSurfaceInfiltration
   if (scalarRainPlusMelt.gt.xMaxInfilRate) then ! infiltration excess surface runoff occurs
    ! saturation excess surface runoff computed by one of the saturation excess methods, remaining surface runoff is infiltration excess
    scalarSurfaceRunoff_IE = scalarSurfaceRunoff - scalarSurfaceRunoff_SE ! infiltration excess surface runoff     
   else ! infiltration excess runoff does not occur
    scalarSurfaceRunoff_SE = scalarSurfaceRunoff ! saturation excess surface runoff 
    scalarSurfaceRunoff_IE = 0._rkind            ! infiltration excess surface runoff 
   end if
  end associate

  ! set surface hydraulic conductivity and diffusivity to missing (not used for flux condition)
  associate(&
   ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   surfaceHydCond => io_surfaceFlx % surfaceHydCond , & ! hydraulic conductivity (m s-1)
   surfaceDiffuse => io_surfaceFlx % surfaceDiffuse   & ! hydraulic diffusivity at the surface (m2 s-1)
  &)
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
   if(err/=0)then; message=trim(message)//'unanticipated error in surfaceFlx subroutine'; return_flag=.true.; return; end if
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
  call update_iLayerFlux_derivatives; if (return_flag) return

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
   dq_dHydStateAbove => out_iLayerFlux % dq_dHydStateAbove, & ! ... matric head or volumetric liquid water in the layer above (m s-1 or s-1)
   dq_dHydStateBelow => out_iLayerFlux % dq_dHydStateBelow, & ! ... matric head or volumetric liquid water in the layer below (m s-1 or s-1)
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
         err=20; message=trim(message)//'only currently implemented for geometric mean -- change local flag'; return_flag=.true.; return
       end if
       ! derivatives in hydraulic conductivity at the layer interface (m s-1)
       dHydCondIface_dVolLiqAbove = dHydCond_dVolLiq(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_rkind/max(iLayerHydCond,verySmaller)
       dHydCondIface_dVolLiqBelow = dHydCond_dVolLiq(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_rkind/max(iLayerHydCond,verySmaller)
       ! derivatives in hydraulic diffusivity at the layer interface (m2 s-1)
       dDiffuseIface_dVolLiqAbove = dDiffuse_dVolLiq(ixUpper)*nodeDiffuseTrial(ixLower) * 0.5_rkind/max(iLayerDiffuse,verySmaller)
       dDiffuseIface_dVolLiqBelow = dDiffuse_dVolLiq(ixLower)*nodeDiffuseTrial(ixUpper) * 0.5_rkind/max(iLayerDiffuse,verySmaller)
       ! derivatives in the flux w.r.t. volumetric liquid water content
       dq_dHydStateAbove = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse/dz + dHydCondIface_dVolLiqAbove
       dq_dHydStateBelow = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse/dz + dHydCondIface_dVolLiqBelow
     case(mixdform)
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
   if(err/=0)then; message=trim(message)//'unanticipated error in iLayerFlux'; return_flag=.true.; return; end if
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
     case default; err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return_flag=.true.; return
   end select 

  end associate
 end subroutine update_qDrainFlux

 subroutine update_qDrainFlux_prescribedHead
  ! ** Update operations for qDrainFlux: prescribed pressure head value at bottom boundary **
  associate(&
   ! input: model control
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

   ! hydrology derivatives
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture); dq_dHydStateUnsat = bottomDiffuse/(nodeDepth/2._rkind)
     case(mixdform); dq_dHydStateUnsat = bottomHydCond/(nodeDepth/2._rkind)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select
   ! energy derivatives
   dq_dNrgStateUnsat = -(dHydCond_dTemp/2._rkind)*(lowerBoundHead  - nodeMatricHeadLiq)/(nodeDepth*0.5_rkind)&
                     & + dHydCond_dTemp/2._rkind
 
  end associate
 end subroutine update_qDrainFlux_prescribedHead

 subroutine update_qDrainFlux_funcBottomHead
  ! ** Update operations for qDrainFlux: prescribed pressure head function at bottom boundary **
  associate(&
   ! input: model control
   ixRichards    => in_qDrainFlux % ixRichards   , &          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   nodeMatricHeadLiq => in_qDrainFlux % nodeMatricHeadLiq, &  ! liquid matric head in the lowest unsaturated node (m)
   nodeVolFracLiq    => in_qDrainFlux % nodeVolFracLiq   , &  ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   nodeHeight => in_qDrainFlux % nodeHeight, &                ! height of the lowest unsaturated soil node (m)
   ! input: derivative in soil water characteristic
   node_dPsi_dTheta    => in_qDrainFlux % node_dPsi_dTheta   , &  ! derivative of the soil moisture characteristic w.r.t. theta (m)
   node_dPsiLiq_dTemp  => in_qDrainFlux % node_dPsiLiq_dTemp , &  ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
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

   ! hydrology derivatives
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * node_dPsi_dTheta*exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
     case(mixdform); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select
   ! energy derivatives
   dq_dNrgStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)*node_dPsiLiq_dTemp/zScale_TOPMODEL

  end associate
 end subroutine update_qDrainFlux_funcBottomHead

 subroutine update_qDrainFlux_freeDrainage
  ! ** Update operations for qDrainFlux: free drainage at bottom boundary **
  associate(&
   ! input: model control
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

   ! hydrology derivatives
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture); dq_dHydStateUnsat = dHydCond_dVolLiq*kAnisotropic
     case(mixdform); dq_dHydStateUnsat = dHydCond_dMatric*kAnisotropic
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select
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

end module soilLiqFlx_module
