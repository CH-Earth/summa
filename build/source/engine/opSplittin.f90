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

module opSplittin_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing double precision number
USE globalData,only:quadMissing      ! missing quadruple precision number

! domain types
USE globalData,only:iname_cas        ! named variables for the canopy air space
USE globalData,only:iname_veg        ! named variables for vegetation
USE globalData,only:iname_snow       ! named variables for snow
USE globalData,only:iname_soil       ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair  ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy  ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy  ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy  ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer   ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer   ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer   ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer   ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:iname_watAquifer ! named variable defining the water storage in the aquifer

! global metadata
USE globalData,only:flux_meta                        ! metadata on the model fluxes
USE globalData,only:diag_meta                        ! metadata on the model diagnostic variables
USE globalData,only:prog_meta                        ! metadata on the model prognostic variables
USE globalData,only:deriv_meta                       ! metadata on the model derivatives
USE globalData,only:flux2state_orig                  ! metadata on flux-to-state mapping (original state variables)
USE globalData,only:flux2state_liq                   ! metadata on flux-to-state mapping (liquid water state variables)

! constants
USE multiconst,only:&
                    gravity,       & ! acceleration of gravity              (m s-2)
                    Tfreeze,       & ! temperature at freezing              (K)
                    LH_fus,        & ! latent heat of fusion                (J kg-1)
                    LH_vap,        & ! latent heat of vaporization          (J kg-1)
                    LH_sub,        & ! latent heat of sublimation           (J kg-1)
                    Cp_air,        & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,      & ! intrinsic density of air             (kg m-3)
                    iden_ice,      & ! intrinsic density of ice             (kg m-3)
                    iden_water       ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookATTR        ! named variables for structure elements
USE var_lookup,only:iLookTYPE        ! named variables for structure elements
USE var_lookup,only:iLookPROG        ! named variables for structure elements
USE var_lookup,only:iLookDIAG        ! named variables for structure elements
USE var_lookup,only:iLookFLUX        ! named variables for structure elements
USE var_lookup,only:iLookFORCE       ! named variables for structure elements
USE var_lookup,only:iLookPARAM       ! named variables for structure elements
USE var_lookup,only:iLookINDEX       ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure

! look up structure for variable types
USE var_lookup,only:iLookVarType

! provide access to the number of flux variables
USE var_lookup,only:nFlux=>maxvarFlux ! number of model flux variables

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,                                                     & ! data vector (i4b)
                    var_d,                                                     & ! data vector (rkind)
                    var_flagVec,                                               & ! data vector with variable length dimension (i4b)
                    var_ilength,                                               & ! data vector with variable length dimension (i4b)
                    var_dlength,                                               & ! data vector with variable length dimension (rkind)
                    zLookup,                                                   & ! lookup tables
                    model_options,                                             & ! defines the model decisions
                    in_type_statefilter,out_type_statefilter,                  & ! classes for stateFilter objects
                    in_type_indexSplit,out_type_indexSplit,                    & ! classes for indexSplit objects
                    in_type_varSubstep,io_type_varSubstep,out_type_varSubstep    ! classes for varSubstep objects

! look-up values for the numerical method
USE mDecisions_module,only:       &
                    numrec       ,& ! home-grown backward Euler solution using free versions of Numerical recipes
                    kinsol       ,& ! SUNDIALS backward Euler solution using Kinsol
                    ida             ! SUNDIALS solution using IDA

! safety: set private unless specified otherwise
implicit none
private
public::opSplittin

! named variables for the coupling method
integer(i4b),parameter  :: fullyCoupled=1             ! 1st try: fully coupled solution
integer(i4b),parameter  :: stateTypeSplit=2           ! 2nd try: separate solutions for each state type

! named variables for the state variable split
integer(i4b),parameter  :: nrgSplit=1                 ! order in sequence for the energy operation
integer(i4b),parameter  :: massSplit=2                ! order in sequence for the mass operation

! named variables for the domain type split
integer(i4b),parameter  :: vegSplit=1                 ! order in sequence for the vegetation split
integer(i4b),parameter  :: snowSplit=2                ! order in sequence for the snow split
integer(i4b),parameter  :: soilSplit=3                ! order in sequence for the soil split
integer(i4b),parameter  :: aquiferSplit=4             ! order in sequence for the aquifer split

! named variables for the solution method
integer(i4b),parameter  :: vector=1                   ! vector solution method
integer(i4b),parameter  :: scalar=2                   ! scalar solution method
integer(i4b),parameter  :: nSolutions=2               ! number of solution methods

! named variables for the switch between states and domains
integer(i4b),parameter  :: fullDomain=1               ! full domain (veg+snow+soil)
integer(i4b),parameter  :: subDomain=2                ! sub domain (veg, snow, soil, and aquifer separately)

! maximum number of possible splits
integer(i4b),parameter  :: nStateTypes=2              ! number of state types (energy, water)
integer(i4b),parameter  :: nDomains=4                 ! number of domains (vegetation, snow, soil, and aquifer)

! control parameters
real(rkind),parameter   :: valueMissing=-9999._rkind     ! missing value
real(rkind),parameter   :: verySmall=1.e-12_rkind        ! a very small number (used to check consistency)
real(rkind),parameter   :: veryBig=1.e+20_rkind          ! a very big number
real(rkind),parameter   :: dx = 1.e-8_rkind              ! finite difference increment

contains


! **********************************************************************************************************
! public subroutine opSplittin: run the coupled energy-mass model for one timestep
!
! The logic of the solver is as follows:
! (1) Attempt different solutions in the following order: (a) fully coupled; (b) split by state type and by
!      domain type for a given energy and mass split (vegetation, snow, and soil); and (c) scalar solution
!      for a given state type and domain subset.
! (2) For a given split, compute a variable number of substeps (in varSubstep).
! **********************************************************************************************************
subroutine opSplittin(&
                      ! input: model control
                      nSnow,                & ! intent(in):    number of snow layers
                      nSoil,                & ! intent(in):    number of soil layers
                      nLayers,              & ! intent(in):    total number of layers
                      nState,               & ! intent(in):    total number of state variables
                      dt,                   & ! intent(in):    time step (s)
                      whole_step,           & ! intent(in):    length of whole step for surface drainage and average flux
                      firstSubStep,         & ! intent(in):    flag to denote first sub-step
                      firstInnerStep,       & ! intent(in):    flag to denote if the last time step in maxstep subStep
                      computeVegFlux,       & ! intent(in):    flag to denote if computing energy flux over vegetation
                       ! input/output: data structures
                      type_data,            & ! intent(in):    type of vegetation and soil
                      attr_data,            & ! intent(in):    spatial attributes
                      forc_data,            & ! intent(in):    model forcing data
                      mpar_data,            & ! intent(in):    model parameters
                      indx_data,            & ! intent(inout): index data
                      prog_data,            & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,            & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,            & ! intent(inout): model fluxes for a local HRU
                      bvar_data,            & ! intent(in):    model variables for the local basin
                      lookup_data,          & ! intent(in):    lookup tables
                      model_decisions,      & ! intent(in):    model decisions
                      ! output: model control
                      dtMultiplier,         & ! intent(out):   substep multiplier (-)
                      tooMuchMelt,          & ! intent(out):   flag to denote that ice is insufficient to support melt
                      stepFailure,          & ! intent(out):   flag to denote step failure
                      ixSolution,           & ! intent(out):   solution method used in this iteration
                      mean_step_dt,         & ! intent(out):   mean solution step for the time step
                      err,message)            ! intent(out):   error code and error message
  ! ---------------------------------------------------------------------------------------
  ! structure allocations
  USE allocspace_module,only:allocLocal                ! allocate local data structures
  ! population/extraction of state vectors
  USE indexState_module,only:indexSplit                ! get state indices
  USE varSubstep_module,only:varSubstep                ! complete substeps for a given split
  ! identify name of variable type (for error message)
  USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: model control
  integer(i4b),intent(in)         :: nSnow                          ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                          ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                        ! total number of layers
  integer(i4b),intent(in)         :: nState                         ! total number of state variables
  real(rkind),intent(in)          :: dt                             ! time step (seconds)
  real(rkind),intent(in)          :: whole_step                     ! length of whole step for surface drainage and average flux
  logical(lgt),intent(in)         :: firstSubStep                   ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(in)         :: computeVegFlux                 ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  ! input/output: data structures
  type(var_i),intent(in)          :: type_data                      ! type of vegetation and soil
  type(var_d),intent(in)          :: attr_data                      ! spatial attributes
  type(var_d),intent(in)          :: forc_data                      ! model forcing data
  type(var_dlength),intent(in)    :: mpar_data                      ! model parameters
  type(var_ilength),intent(inout) :: indx_data                      ! indices for a local HRU
  type(var_dlength),intent(inout) :: prog_data                      ! prognostic variables for a local HRU
  type(var_dlength),intent(inout) :: diag_data                      ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data                      ! model fluxes for a local HRU
  type(var_dlength),intent(in)    :: bvar_data                      ! model variables for the local basin
  type(zLookup),    intent(in)    :: lookup_data                    ! lookup tables
  type(model_options),intent(in)  :: model_decisions(:)             ! model decisions
  ! output: model control
  real(rkind),intent(out)         :: dtMultiplier                   ! substep multiplier (-)
  logical(lgt),intent(out)        :: tooMuchMelt                    ! flag to denote that ice is insufficient to support melt
  logical(lgt),intent(out)        :: stepFailure                    ! flag to denote step failure
  integer(i4b),intent(out)        :: ixSolution                     ! index of solution method (1,2)
  real(rkind),intent(out)         :: mean_step_dt                   ! mean solution step for the time step
  integer(i4b),intent(out)        :: err                            ! error code
  character(*),intent(out)        :: message                        ! error message
  ! ---------------------------------------------------------------------------------------
  ! * general local variables
  ! ---------------------------------------------------------------------------------------
  character(LEN=256)              :: cmessage                       ! error message of downwind routine
  integer(i4b)                    :: minLayer                       ! the minimum layer used in assigning flags for flux aggregations
  integer(i4b)                    :: iOffset                        ! offset to account for different indices in the soil domain
  integer(i4b)                    :: iMin(1),iMax(1)                ! bounds of a given vector
  integer(i4b)                    :: iLayer,jLayer                  ! index of model layer
  integer(i4b)                    :: iSoil                          ! index of soil layer
  integer(i4b)                    :: iVar                           ! index of variables in data structures
  logical(lgt)                    :: firstSuccess                   ! flag to define the first success
  logical(lgt)                    :: firstFluxCall                  ! flag to define the first flux call
  logical(lgt)                    :: reduceCoupledStep              ! flag to define the need to reduce the length of the coupled step
  logical(lgt)                    :: return_flag                    ! flag to indicate the execution of a return statement
  type(var_dlength)               :: prog_temp                      ! temporary model prognostic variables
  type(var_dlength)               :: diag_temp                      ! temporary model diagnostic variables
  type(var_dlength)               :: flux_temp                      ! temporary model fluxes
  type(var_dlength)               :: flux_mean                      ! mean model fluxes
  type(var_dlength)               :: flux_mntemp                    ! temporary mean model fluxes
  type(var_dlength)               :: deriv_data                     ! derivatives in model fluxes w.r.t. relevant state variables
  ! ------------------------------------------------------------------------------------------------------
  ! * operator splitting
  ! ------------------------------------------------------------------------------------------------------
  ! minimum timestep
  real(rkind),parameter           :: dtmin_coupled=1800._rkind      ! minimum time step for the fully coupled solution (seconds)
  real(rkind),parameter           :: dtmin_split=60._rkind          ! minimum time step for the fully split solution (seconds)
  real(rkind),parameter           :: dtmin_scalar=10._rkind         ! minimum time step for the scalar solution (seconds)
  real(rkind)                     :: dt_min                         ! minimum time step (seconds)
  real(rkind)                     :: dtInit                         ! initial time step (seconds)
  ! explicit error tolerance (depends on state type split, so defined here)
  real(rkind),parameter           :: errorTolLiqFlux=0.01_rkind     ! error tolerance in the explicit solution (liquid flux)
  real(rkind),parameter           :: errorTolNrgFlux=10._rkind      ! error tolerance in the explicit solution (energy flux)
  ! number of substeps taken for a given split
  integer(i4b)                    :: nSubsteps                      ! number of substeps taken for a given split
  ! named variables defining the coupling and solution method
  integer(i4b)                    :: ixCoupling                     ! index of coupling method (1,2)
  integer(i4b)                    :: ixStateThenDomain              ! switch between the state and domain (1,2)
  integer(i4b)                    :: tryDomainSplit                 ! (0,1) - flag to try the domain split
  ! actual number of splits
  integer(i4b)                    :: nStateTypeSplit                ! number of splits for the state type
  integer(i4b)                    :: nDomainSplit                   ! number of splits for the domain
  integer(i4b)                    :: nStateSplit                    ! number of splits for the states within a given domain
  ! indices for the state type and the domain split
  integer(i4b)                    :: iStateTypeSplit                ! index of the state type split
  integer(i4b)                    :: iDomainSplit                   ! index of the domain split
  integer(i4b)                    :: iStateSplit                    ! index of the state split
  ! flux masks
  logical(lgt)                    :: neededFlux(nFlux)              ! .true. if flux is needed at all
  logical(lgt)                    :: desiredFlux                    ! .true. if flux is desired for a given split
  type(var_ilength)               :: fluxCount                      ! number of times each flux is updated (should equal nSubsteps)
  type(var_flagVec)               :: fluxMask                       ! mask defining model fluxes
  ! state masks
  integer(i4b),dimension(nState)  :: stateCheck                     ! number of times each state variable is updated (should equal 1)
  logical(lgt),dimension(nState)  :: stateMask                      ! mask defining desired state variables
  integer(i4b)                    :: nSubset                        ! number of selected state variables for a given split
  ! flags
  logical(lgt)                    :: failure                        ! flag to denote failure of substepping
  logical(lgt)                    :: doAdjustTemp                   ! flag to adjust temperature after the mass split
  logical(lgt)                    :: failedMinimumStep              ! flag to denote failure of substepping for a given split
  integer(i4b)                    :: ixSaturation                   ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  integer(i4b)                    :: nCoupling
  logical(lgt)                    :: firstInnerStep                 ! flag to denote if the first time step in maxstep subStep
  ! mean steps 
  real(rkind)                     :: mean_step_state                ! mean step over the state (with or without domain splits)
  real(rkind)                     :: mean_step_solution             ! mean step for a solution (scalar or vector)
  logical(lgt)                    :: addFirstFlux                   ! flag to add the first flux to the mask
  ! loop control
  logical(lgt)                    :: exit_coupling,exit_stateTypeSplitting,exit_stateThenDomain,exit_domainSplit,exit_solution,exit_stateSplit
  logical(lgt)                    :: cycle_coupling,cycle_stateTypeSplitting,cycle_stateThenDomain,cycle_domainSplit,cycle_solution,cycle_stateSplit
  ! ------------------------ classes for subroutine arguments (classes defined in data_types module) ------------------------
  !      ** intent(in) arguments **         ||       ** intent(inout) arguments **        ||      ** intent(out) arguments **
  type(in_type_stateFilter) :: in_stateFilter;                                            type(out_type_stateFilter) :: out_stateFilter; ! stateFilter arguments
  type(in_type_indexSplit)  :: in_indexSplit;                                             type(out_type_indexSplit)  :: out_indexSplit;  ! indexSplit arguments
  type(in_type_varSubstep)  :: in_varSubstep;  type(io_type_varSubstep) :: io_varSubstep; type(out_type_varSubstep)  :: out_varSubstep;  ! varSubstep arguments
  ! ---------------------------------------------------------------------------------------

  call initialize_coupling; if (return_flag.eqv..true.) return ! select coupling options and allocate memory - return if error occurs
  coupling: do ixCoupling=1,nCoupling                          ! loop through different coupling strategies

    call initialize_stateTypeSplitting; if (return_flag.eqv..true.) return ! setup steps for stateTypeSplitting loop - return if error occurs
    stateTypeSplitting: do iStateTypeSplit=1,nStateTypeSplit               ! state splitting loop

      ! first try the state type split, then try the domain split within a given state type
      call initialize_stateThenDomain ! setup steps for stateThenDomain loop -- identify state-specific variables for a given state split
      stateThenDomain: do ixStateThenDomain=1,1+tryDomainSplit ! 1=state type split; 2=domain split within a given state type

        call initialize_domainSplit; if (return_flag.eqv..true.) return ! setup steps for domainSplit loop - return if error occurs
        domainSplit: do iDomainSplit=1,nDomainSplit                     ! domain splitting loop

          solution: do ixSolution=1,nSolutions ! trial with the vector then scalar solution

            call initialize_stateSplit; if (return_flag.eqv..true.) return ! setup steps for stateSplit loop - return if error occurs
            stateSplit: do iStateSplit=1,nStateSplit ! loop through layers (NOTE: nStateSplit=1 for the vector solution, hence no looping)

              ! define state subsets for a given split...
              call update_stateFilter; if (return_flag.eqv..true.) return ! get the mask for the state subset - return for a non-zero error code
              call validate_split ! verify that the split is valid
              if (cycle_domainSplit) cycle domainSplit
              if (cycle_solution) cycle solution
              if (return_flag.eqv..true.) return ! return for a non-zero error code
              call save_recover ! save/recover copies of variables and fluxes

              ! assemble vectors for a given split...
              call get_split_indices; if (return_flag.eqv..true.) return ! get indices for a given split - return for a non-zero error code
              call update_fluxMask; if (return_flag.eqv..true.) return ! define the mask of the fluxes used - return for a non-zero error code

              call solve_subset; if (return_flag.eqv..true.) return ! solve variable subset for one time step - return for a positive error code

              call assess_solution; if (return_flag.eqv..true.) return ! is solution a success or failure? - return for a recovering solution

              call try_other_solution_methods                  ! if solution failed to converge, try other splitting methods 
              if (cycle_coupling)        cycle coupling        ! exit loops if necessary
              if (cycle_stateThenDomain) cycle stateThenDomain
              if (cycle_solution)        cycle solution

              call confirm_variable_updates; if (return_flag.eqv..true.) return ! check that state variables updated - return if error 

              call success_check ! check for success
              if (exit_stateThenDomain) exit stateThenDomain ! exit loops if necessary
              if (exit_solution) exit solution
              if (return_flag.eqv..true.) return             ! return if error 

            end do stateSplit ! solution with split layers

          end do solution        ! trial with the full layer solution then the split layer solution
          call finalize_solution ! final steps following solution loop

        end do domainSplit ! domain type splitting loop

      end do stateThenDomain        ! switch between the state type and domain type splitting
      call finalize_stateThenDomain ! final steps following the stateThenDomain loop

    end do stateTypeSplitting                                          ! state type splitting loop
    call finalize_stateTypeSplitting; if (exit_coupling) exit coupling ! success = exit the coupling loop

  end do coupling        ! loop over coupling methods
  call finalize_coupling ! check variables and fluxes, and apply step halving if needed

 contains

  subroutine initialize_coupling
   ! *** initial steps for coupling loop ***
   ! initialize error control
   err=0; message="opSplittin/"

   associate(ixNumericalMethod => model_decisions(iLookDECISIONS%num_method)%iDecision) ! intent(in): [i4b] choice of numerical solver
    ! we just solve the fully coupled problem if IDA for now, splitting can happen otherwise
    select case(ixNumericalMethod)
     case(ida);            nCoupling = 1
     case(kinsol, numrec); nCoupling = 2
    end select
   end associate

   ! set the global print flag
   globalPrintFlag=.false.

   if (globalPrintFlag) print *, trim(message), dt

   ! initialize the first success call
   firstSuccess=.false.
   if (.not.firstInnerStep) firstSuccess=.true.

   ! initialize the flags
   tooMuchMelt=.false.  ! too much melt (merge snow layers)
   stepFailure=.false.  ! step failure

   ! initialize flag for the success of the substepping
   failure=.false.

   ! initialize the flux check
   neededFlux(:) = .false.

   ! initialize the state check
   stateCheck(:) = 0

   ! allocate local structures based on the number of snow and soil layers
   call allocate_memory
   if (return_flag.eqv..true.) return ! return if an error occurs during memory allocation 

   ! intialize the flux counter
   do iVar=1,size(flux_meta)  ! loop through fluxes
     fluxCount%var(iVar)%dat(:) = 0
   end do

   ! initialize the model fluxes
   do iVar=1,size(flux_meta)  ! loop through fluxes
    if (flux2state_orig(iVar)%state1==integerMissing .and. flux2state_orig(iVar)%state2==integerMissing) cycle ! flux does not depend on state (e.g., input)
    if (flux2state_orig(iVar)%state1==iname_watCanopy .and. .not.computeVegFlux) cycle ! use input fluxes in cases where there is no canopy
    if (firstInnerStep) flux_data%var(iVar)%dat(:) = 0._rkind
    flux_mean%var(iVar)%dat(:) = 0._rkind
   end do

   ! initialize derivatives
   do iVar=1,size(deriv_meta)
    deriv_data%var(iVar)%dat(:) = 0._rkind
   end do
  end subroutine initialize_coupling

  subroutine allocate_memory
   ! *** allocate memory for local structures ***
   return_flag=.false. ! initialize flag

   ! allocate space for the flux mask (used to define when fluxes are updated)
   call allocLocal(flux_meta(:),fluxMask,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the flux count (used to check that fluxes are only updated once)
   call allocLocal(flux_meta(:),fluxCount,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the temporary prognostic variable structure
   call allocLocal(prog_meta(:),prog_temp,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the temporary diagnostic variable structure
   call allocLocal(diag_meta(:),diag_temp,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the temporary flux variable structure
   call allocLocal(flux_meta(:),flux_temp,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the mean flux variable structure
   call allocLocal(flux_meta(:),flux_mean,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the temporary mean flux variable structure
   call allocLocal(flux_meta(:),flux_mntemp,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the derivative structure
   call allocLocal(deriv_meta(:),deriv_data,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end subroutine allocate_memory

  subroutine finalize_coupling
   ! *** final operations for coupling loop ***
   ! check that all state variables were updated
   if (any(stateCheck==0)) then
    message=trim(message)//'some state variables were not updated!'
    err=20; return
   endif

   ! check that the desired fluxes were computed
   do iVar=1,size(flux_meta)
    if (neededFlux(iVar) .and. any(fluxCount%var(iVar)%dat==0)) then
     print*, 'fluxCount%var(iVar)%dat = ', fluxCount%var(iVar)%dat
     message=trim(message)//'flux '//trim(flux_meta(iVar)%varname)//' was not computed'
     err=20; return
    end if
   end do

   ! use step halving if unable to complete the fully coupled solution in one substep
   if (ixCoupling/=fullyCoupled .or. nSubsteps>1) dtMultiplier=0.5_rkind
  end subroutine finalize_coupling

  subroutine initialize_stateTypeSplitting
   ! *** Initial steps to prepare for iterations of the stateTypeSplit loop ***
   return_flag=.false. ! initialize flag
   ! initialize the time step
   dtInit = min(merge(dt,            dtmin_coupled, ixCoupling==fullyCoupled), dt) ! initial time step
   dt_min = min(merge(dtmin_coupled, dtmin_split,   ixCoupling==fullyCoupled), dt) ! minimum time step

   ! keep track of the number of state splits
   associate(numberStateSplit => indx_data%var(iLookINDEX%numberStateSplit)%dat(1)) ! intent(inout): [i4b] number of state splitting solutions
    if (ixCoupling/=fullyCoupled) numberStateSplit = numberStateSplit + 1
   end associate

   ! define the number of operator splits for the state type
   select case(ixCoupling)
    case(fullyCoupled); nStateTypeSplit=1
    case(stateTypeSplit); nStateTypeSplit=nStateTypes
    case default; err=20; message=trim(message)//'coupling case not found'; return_flag=.true.; return
   end select  ! operator splitting option

   ! define if we wish to try the domain split
   select case(ixCoupling)
    case(fullyCoupled);   tryDomainSplit=0
    case(stateTypeSplit); tryDomainSplit=1
    case default; err=20; message=trim(message)//'coupling case not found'; return_flag=.true.; return
   end select  ! operator splitting option

   mean_step_dt = 0._rkind ! initialize mean step for the time step
   addFirstFlux = .true.     ! flag to add the first flux to the mask
  end subroutine initialize_stateTypeSplitting

  subroutine finalize_stateTypeSplitting
   ! *** Final operations subsequent to the stateTypeSplitting loop ***
   exit_coupling=.false. ! initialize flag for loop control 
   if (ixCoupling==fullyCoupled .and. .not.failure) then; exit_coupling=.true.; return; end if ! success = exit the coupling loop in opSplittin
  end subroutine finalize_stateTypeSplitting

  subroutine initialize_stateThenDomain
   ! *** Identify state-specific variables for a given state split ***
   doAdjustTemp = (ixCoupling/=fullyCoupled .and. iStateTypeSplit==massSplit) ! flag to adjust the temperature
   associate(&
    ixStateType => indx_data%var(iLookINDEX%ixStateType)%dat, & ! intent(in): [i4b(:)] indices defining the type of the state
    ixHydCanopy => indx_data%var(iLookINDEX%ixHydCanopy)%dat, & ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
    ixHydLayer  => indx_data%var(iLookINDEX%ixHydLayer)%dat   ) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain

    ! modify the state type names associated with the state vector
    if (ixCoupling/=fullyCoupled .and. iStateTypeSplit==massSplit) then ! if modifying state variables for the mass split
      if (computeVegFlux) then
        where(ixStateType(ixHydCanopy)==iname_watCanopy) ixStateType(ixHydCanopy)=iname_liqCanopy
      end if
      where(ixStateType(ixHydLayer)==iname_watLayer) ixStateType(ixHydLayer)=iname_liqLayer
      where(ixStateType(ixHydLayer)==iname_matLayer) ixStateType(ixHydLayer)=iname_lmpLayer
    end if
   end associate
  end subroutine initialize_stateThenDomain

  subroutine finalize_stateThenDomain
   ! *** Final steps following the stateThenDomain loop ***
   ! sum the mean steps for the time step over each state type split
   select case(ixStateThenDomain) 
     case(fullDomain); mean_step_dt = mean_step_dt + mean_step_solution/nStateTypeSplit
     case(subDomain);  mean_step_dt = mean_step_dt + mean_step_state/nStateTypeSplit
   end select
   associate(&
    ixStateType => indx_data%var(iLookINDEX%ixStateType)%dat, & ! intent(in): [i4b(:)] indices defining the type of the state
    ixHydCanopy => indx_data%var(iLookINDEX%ixHydCanopy)%dat, & ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
    ixHydLayer  => indx_data%var(iLookINDEX%ixHydLayer)%dat   ) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain

    ! * reset state variables for the mass split...
    ! modify the state type names associated with the state vector
    if (ixCoupling/=fullyCoupled .and. iStateTypeSplit==massSplit) then ! if modifying state variables for the mass split
      if (computeVegFlux) then
        where(ixStateType(ixHydCanopy)==iname_liqCanopy) ixStateType(ixHydCanopy)=iname_watCanopy
      end if
      where(ixStateType(ixHydLayer)==iname_liqLayer) ixStateType(ixHydLayer)=iname_watLayer
      where(ixStateType(ixHydLayer)==iname_lmpLayer) ixStateType(ixHydLayer)=iname_matLayer
    end if
   end associate  
  end subroutine finalize_stateThenDomain

  subroutine initialize_domainSplit
   ! *** initial operations to set up domainSplit loop ***
   return_flag=.false. ! initialize flag
   associate(numberDomainSplitNrg => indx_data%var(iLookINDEX%numberDomainSplitNrg )%dat(1),& ! intent(inout): [i4b] number of domain splitting solutions for energy (-)
             numberDomainSplitMass => indx_data%var(iLookINDEX%numberDomainSplitMass)%dat(1) )! intent(inout): [i4b] number of domain splitting solutions for mass   (-)
    ! keep track of the number of domain splits
    if (iStateTypeSplit==nrgSplit  .and. ixStateThenDomain==subDomain) numberDomainSplitNrg  = numberDomainSplitNrg  + 1
    if (iStateTypeSplit==massSplit .and. ixStateThenDomain==subDomain) numberDomainSplitMass = numberDomainSplitMass + 1
   end associate

   ! define the number of domain splits for the state type
   select case(ixStateThenDomain)
     case(fullDomain); nDomainSplit=1
     case(subDomain);  nDomainSplit=nDomains
     case default; err=20; message=trim(message)//'coupling case not found';
      return_flag=.true. ! return statement required in opSplittin
      return
   end select

   ! check that we haven't split the domain when we are fully coupled
   if (ixCoupling==fullyCoupled .and. nDomainSplit==nDomains) then
     message=trim(message)//'cannot split domains when fully coupled'
     return_flag=.true. ! return statement required in opSplittin
     err=20; return
   end if

   mean_step_state = 0._rkind ! initialize mean step for state
  end subroutine initialize_domainSplit

  subroutine finalize_solution
   ! *** final operations following solution loop ***
   ! sum the mean steps for the state over each domain split
   mean_step_state = mean_step_state + mean_step_solution/nDomainSplit
  end subroutine finalize_solution

  subroutine initialize_stateSplit
   ! *** initial operations to set up stateSplit loop ***
   return_flag=.false. ! initialize flag
   mean_step_solution = 0._rkind ! initialize mean step for a solution

   ! initialize error control
   err=0; message="opSplittin/"

   ! refine the time step
   if (ixSolution==scalar) then
    dtInit = min(dtmin_split, dt)    ! initial time step
    dt_min = min(dtmin_scalar, dt)   ! minimum time step
   end if

   ! initialize the first flux call
   firstFluxCall=.true.
   if (.not.firstInnerStep) firstFluxCall=.false.

   ! get the number of split layers
   select case(ixSolution)
    case(vector); nStateSplit=1
    case(scalar); nStateSplit=count(stateMask)
    case default; err=20; message=trim(message)//'unknown solution method'; 
     return_flag=.true. ! return statement required in opSplittin
     return
   end select
  end subroutine initialize_stateSplit

  ! **** stateFilter ****
  subroutine initialize_stateFilter
   call in_stateFilter % initialize(ixCoupling,ixSolution,ixStateThenDomain,iStateTypeSplit,iDomainSplit,iStateSplit)
  end subroutine initialize_stateFilter

  subroutine finalize_stateFilter
   call out_stateFilter % finalize(nSubset,err,cmessage)
  end subroutine finalize_stateFilter

  ! **** indexSplit ****
  subroutine initialize_indexSplit
   call in_indexSplit % initialize(nSnow,nSoil,nLayers,nSubset)
  end subroutine initialize_indexSplit

  subroutine finalize_indexSplit
   call out_indexSplit % finalize(err,cmessage)
  end subroutine finalize_indexSplit
  ! **** end indexSplit ****

  ! **** varSubstep ****
  subroutine initialize_varSubstep
   call in_varSubstep % initialize(dt,dtInit,dt_min,whole_step,nSubset,doAdjustTemp,firstSubStep,computeVegFlux,ixSolution,scalar,iStateSplit,fluxMask)
   call io_varSubstep % initialize(firstFluxCall,fluxCount,ixSaturation)
  end subroutine initialize_varSubstep

  subroutine finalize_varSubstep
   call io_varSubstep  % finalize(firstFluxCall,fluxCount,ixSaturation)
   call out_varSubstep % finalize(dtMultiplier,nSubsteps,failedMinimumStep,reduceCoupledStep,tooMuchMelt,err,cmessage)
  end subroutine finalize_varSubstep

  subroutine solve_subset 
   ! *** Solve variable subset for one time step ***
   return_flag=.false. ! initialize flag
   ! keep track of the number of scalar solutions
   associate(numberScalarSolutions => indx_data%var(iLookINDEX%numberScalarSolutions)%dat(1)) ! intent(inout): [i4b] number of scalar solutions
    if (ixSolution==scalar) numberScalarSolutions = numberScalarSolutions + 1
   end associate

   ! solve variable subset for one full time step
   call initialize_varSubstep
   call varSubstep(in_varSubstep,io_varSubstep,&                                            ! intent(inout): class objects for model control
                   model_decisions,lookup_data,type_data,attr_data,forc_data,mpar_data,&    ! intent(inout): data structures for model properties
                   indx_data,prog_data,diag_data,flux_data,flux_mean,deriv_data,bvar_data,&
                   out_varSubstep)                                                          ! intent(out): class object for model control
   call finalize_varSubstep
   if (err/=0) then 
    message=trim(message)//trim(cmessage) 
    if (err>0) then ! return for positive error codes
     return_flag=.true.; return
    end if 
   end if ! error control
  end subroutine solve_subset 

  subroutine assess_solution
   ! *** determine whether solution is a success or a failure ***
   return_flag=.false. ! initialize flag

   ! reduce coupled step if failed the minimum step for the scalar solution
   if (failedMinimumStep .and. ixSolution==scalar) reduceCoupledStep=.true.

   ! if too much melt (or some other need to reduce the coupled step) then return
   ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
   if (tooMuchMelt .or. reduceCoupledStep) then
     stepFailure=.true.
     err=0 ! recovering
     return_flag=.true. ! return statement required in opSplittin
     return
   end if

   ! define failure
   failure = (failedMinimumStep .or. err<0)
   if (.not.failure) firstSuccess=.true.

   ! if failed, need to reset the flux counter
   if (failure) then
     do iVar=1,size(flux_meta)
       iMin=lbound(flux_data%var(iVar)%dat)
       iMax=ubound(flux_data%var(iVar)%dat)
       do iLayer=iMin(1),iMax(1)
         if (fluxMask%var(iVar)%dat(iLayer)) fluxCount%var(iVar)%dat(iLayer) = fluxCount%var(iVar)%dat(iLayer) - nSubsteps
       end do
     end do
   end if
  end subroutine assess_solution

  subroutine try_other_solution_methods 
   ! *** if solution failed to converge, try other splitting methods *** 
   ! initialize flags
   cycle_coupling=.false.
   cycle_stateThenDomain=.false.
   cycle_solution=.false.

   ! try the fully split solution if failed to converge with a minimum time step in the coupled solution
   if (ixCoupling==fullyCoupled .and. failure) then; cycle_coupling=.true.; return; end if! return required to execute cycle statement in opSplittin

   ! try the scalar solution if failed to converge with a minimum time step in the split solution
   if (ixCoupling/=fullyCoupled) then
     select case(ixStateThenDomain)
       case(fullDomain); if (failure) cycle_stateThenDomain=.true.; return ! return required to execute cycle statement in opSplittin
       case(subDomain);  if (failure) cycle_solution=.true.; return
       case default; err=20; message=trim(message)//'unknown ixStateThenDomain case'
     end select
   end if
  end subroutine try_other_solution_methods 

  subroutine update_stateFilter
   ! *** Get the mask for the state subset ***
   return_flag=.false. ! initialize flag
   call initialize_stateFilter
   call stateFilter(in_stateFilter,indx_data,stateMask,out_stateFilter)
   call finalize_stateFilter
   if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! error control
  end subroutine update_stateFilter

  subroutine validate_split 
   ! *** Verify that the split is valid ***
   ! initialize flags
   cycle_domainSplit=.false.
   cycle_solution=.false.
   return_flag=.false.
   ! check that state variables exist
   if (nSubset==0) then; cycle_domainSplit=.true.; return; end if

   ! avoid redundant case where vector solution is of length 1
   if (ixSolution==vector .and. count(stateMask)==1) then; cycle_solution=.true.; return; end if

   ! check that we do not attempt the scalar solution for the fully coupled case
   if (ixCoupling==fullyCoupled .and. ixSolution==scalar) then
     message=trim(message)//'only apply the scalar solution to the fully split coupling strategy'
     err=20; return_flag=.true.; return
   end if

   ! reset the flag for the first flux call
   if (.not.firstSuccess) firstFluxCall=.true.
  end subroutine validate_split 

  subroutine save_recover
   ! save/recover copies of prognostic variables
   do iVar=1,size(prog_data%var)
     select case(failure)
       case(.false.); prog_temp%var(iVar)%dat(:) = prog_data%var(iVar)%dat(:)
       case(.true.);  prog_data%var(iVar)%dat(:) = prog_temp%var(iVar)%dat(:)
     end select
   end do 

   ! save/recover copies of diagnostic variables
   do iVar=1,size(diag_data%var)
     select case(failure)
       case(.false.); diag_temp%var(iVar)%dat(:) = diag_data%var(iVar)%dat(:)
       case(.true.);  diag_data%var(iVar)%dat(:) = diag_temp%var(iVar)%dat(:)
     end select
   end do 

   ! save/recover copies of model fluxes and mean fluxes
   do iVar=1,size(flux_data%var)
     select case(failure)
       case(.false.)
         flux_temp%var(iVar)%dat(:)   = flux_data%var(iVar)%dat(:)
         flux_mntemp%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:)
         addFirstFlux = .false.
       case(.true.)
         flux_data%var(iVar)%dat(:)   = flux_temp%var(iVar)%dat(:)
         flux_mean%var(iVar)%dat(:)   = flux_mntemp%var(iVar)%dat(:)
         if (addFirstFlux) addFirstFlux = .true.
     end select
   end do
  end subroutine save_recover

  subroutine get_split_indices
   ! *** Get indices for a given split ***
   return_flag=.false. ! initialize flag
   call initialize_indexSplit
   call indexSplit(in_indexSplit,stateMask,indx_data,out_indexSplit)
   call finalize_indexSplit
   if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end subroutine get_split_indices

  subroutine confirm_variable_updates
   ! *** check that state variables updated ***
   return_flag=.false. ! set flag
   ! check that state variables updated
   where(stateMask) stateCheck = stateCheck+1
   if (any(stateCheck>1)) then
     message=trim(message)//'state variable updated more than once!'
     err=20; return_flag=.true.; return
   end if
  end subroutine confirm_variable_updates

  subroutine success_check
   ! initialize flags
   return_flag=.false.
   exit_stateThenDomain=.false.
   exit_solution=.false.
   ! success = exit solution
   if (.not.failure) then
     ! sum the mean steps for the successful solution type
     mean_step_solution = mean_step_solution + (dt/nSubsteps)/nStateSplit
     select case(ixStateThenDomain)
       case(fullDomain); if (iStateSplit==nStateSplit) exit_stateThenDomain=.true. ! exit stateThenDomain
       case(subDomain);  if (iStateSplit==nStateSplit) exit_solution=.true. ! exit solution
       case default; err=20; message=trim(message)//'unknown ixStateThenDomain case'
     end select
   else ! failure
     call check_failure; return_flag=.true.; return ! check reason for failure and return
   end if  ! success check
  end subroutine success_check

  subroutine check_failure
   ! *** Analyze reason for failure ***
   if (ixSolution==scalar) then ! check that we did not fail for the scalar solution (last resort)
     message=trim(message)//'failed the minimum step for the scalar solution'
     err=20; return
   else ! check for an unexpected failure
     message=trim(message)//'unexpected failure'
     err=20; return
   end if
  end subroutine check_failure

  subroutine update_fluxMask
   ! *** update the fluxMask data structure ***
   return_flag=.false. ! initialize flag
 
   do iVar=1,size(flux_meta) ! loop through flux variables

    if (ixCoupling==fullyCoupled) then ! * identify flux mask for the fully coupled solution
     associate(ixStateType_subset => indx_data%var(iLookINDEX%ixStateType_subset)%dat) ! intent(in): [i4b(:)] indices of state types
      desiredFlux = any(ixStateType_subset==flux2state_orig(iVar)%state1) .or. any(ixStateType_subset==flux2state_orig(iVar)%state2)
     end associate

     ! make sure firstFluxCall fluxes are included in the mask
     if (firstFluxCall .and. addFirstFlux) then 
      if (iVar==iLookFlux%scalarSoilResistance) desiredFlux = .true.
      if (iVar==iLookFlux%scalarStomResistSunlit) desiredFlux = .true.
      if (iVar==iLookFlux%scalarStomResistShaded) desiredFlux = .true.
      if (iVar==iLookFlux%scalarPhotosynthesisSunlit) desiredFlux = .true.
      if (iVar==iLookFlux%scalarPhotosynthesisShaded) desiredFlux = .true.
     end if

     fluxMask%var(iVar)%dat = desiredFlux

    else ! * identify flux mask for the split solution

     associate(ixStateType_subset => indx_data%var(iLookINDEX%ixStateType_subset)%dat) ! intent(in): [i4b(:)] indices of state types
      select case(iStateTypeSplit) ! identify the flux mask for a given state split
       case(nrgSplit);  desiredFlux = any(ixStateType_subset==flux2state_orig(iVar)%state1) .or. any(ixStateType_subset==flux2state_orig(iVar)%state2)
       case(massSplit); desiredFlux = any(ixStateType_subset==flux2state_liq(iVar)%state1)  .or. any(ixStateType_subset==flux2state_liq(iVar)%state2)
       case default; err=20; message=trim(message)//'unable to identify split based on state type'; return_flag=.true.; return
      end select
     end associate

     ! make sure firstFluxCall fluxes are included in the mask
     if (firstFluxCall .and. addFirstFlux) then 
      if (iVar==iLookFlux%scalarSoilResistance) desiredFlux = .true.
      if (iVar==iLookFlux%scalarStomResistSunlit) desiredFlux = .true.
      if (iVar==iLookFlux%scalarStomResistShaded) desiredFlux = .true.
      if (iVar==iLookFlux%scalarPhotosynthesisSunlit) desiredFlux = .true.
      if (iVar==iLookFlux%scalarPhotosynthesisShaded) desiredFlux = .true.
     end if

     if (nDomains==1) then ! no domain splitting
      fluxMask%var(iVar)%dat = desiredFlux
     else ! domain splitting
      fluxMask%var(iVar)%dat = .false. ! initialize to .false.
      if (desiredFlux) then ! only need to proceed if the flux is desired
       select case(iDomainSplit) ! different domain splitting operations
        case(vegSplit) ! canopy fluxes -- (:1) gets the upper boundary(0) if it exists
         if (ixSolution==vector) then ! vector solution (should only be present for energy)
          fluxMask%var(iVar)%dat(:1) = desiredFlux
          if (ixStateThenDomain>1 .and. iStateTypeSplit/=nrgSplit) then
           message=trim(message)//'only expect a vector solution for the vegetation domain for energy'
           err=20; return_flag=.true.; return
          end if
         else                         ! scalar solution
          fluxMask%var(iVar)%dat(:1) = desiredFlux
         end if
        case(snowSplit,soilSplit) ! fluxes through snow and soil

         do iLayer=1,nLayers! loop through layers
          associate(ixLayerActive => indx_data%var(iLookINDEX%ixLayerActive)%dat) ! intent(in): [i4b(:)] indices for all active layers (inactive=integerMissing)
           if (ixLayerActive(iLayer)/=integerMissing) then

            ! get the offset (ixLayerActive=1,2,3,...nLayers, and soil vectors nSnow+1, nSnow+2, ..., nLayers)
            iOffset = merge(nSnow, 0, flux_meta(iVar)%vartype==iLookVarType%midSoil .or. flux_meta(iVar)%vartype==iLookVarType%ifcSoil)
            jLayer  = iLayer-iOffset

            ! identify the minimum layer
            select case(flux_meta(iVar)%vartype)
             case(iLookVarType%ifcToto, iLookVarType%ifcSnow, iLookVarType%ifcSoil); minLayer=merge(jLayer-1, jLayer, jLayer==1)
             case(iLookVarType%midToto, iLookVarType%midSnow, iLookVarType%midSoil); minLayer=jLayer
             case default; minLayer=integerMissing
            end select

            ! set desired layers
            select case(flux_meta(iVar)%vartype)
             case(iLookVarType%midToto,iLookVarType%ifcToto);                    fluxMask%var(iVar)%dat(minLayer:jLayer) = desiredFlux
             case(iLookVarType%midSnow,iLookVarType%ifcSnow); if (iLayer<=nSnow) fluxMask%var(iVar)%dat(minLayer:jLayer) = desiredFlux
             case(iLookVarType%midSoil,iLookVarType%ifcSoil); if (iLayer> nSnow) fluxMask%var(iVar)%dat(minLayer:jLayer) = desiredFlux
            end select

            ! add hydrology states for scalar variables
            if (iStateTypeSplit==massSplit .and. flux_meta(iVar)%vartype==iLookVarType%scalarv) then
             select case(iDomainSplit)
              case(snowSplit); if(iLayer==nSnow)   fluxMask%var(iVar)%dat = desiredFlux
              case(soilSplit); if(iLayer==nSnow+1) fluxMask%var(iVar)%dat = desiredFlux
             end select
            end if  ! if hydrology split and scalar

           end if    ! if the layer is active
          end associate
         end do   ! end looping through layers

        case(aquiferSplit) ! fluxes through aquifer
         fluxMask%var(iVar)%dat(:) = desiredFlux ! only would be firstFluxCall variables, no aquifer fluxes
        case default; err=20; message=trim(message)//'unable to identify split based on domain type'; return_flag=.true.; return ! check
       end select  ! domain split
      end if  ! end if flux is desired
     end if  ! end if domain splitting
    end if  ! end if not fully coupled

    ! define if the flux is desired
    if (desiredFlux) neededFlux(iVar)=.true.
    !if(desiredFlux) print*, flux_meta(iVar)%varname, fluxMask%var(iVar)%dat

    if ( globalPrintFlag .and. count(fluxMask%var(iVar)%dat)>0 ) print*, trim(flux_meta(iVar)%varname) ! * check

   end do  ! end looping through fluxes

  end subroutine update_fluxMask

end subroutine opSplittin


! **********************************************************************************************************
! private subroutine stateFilter: get a mask for the desired state variables
! **********************************************************************************************************
subroutine stateFilter(in_stateFilter,indx_data,stateMask,out_stateFilter)

  USE indexState_module,only:indxSubset                               ! get state indices
  implicit none
  ! input
  type(in_type_stateFilter),intent(in)   :: in_stateFilter            ! indices
  type(var_ilength),intent(inout)        :: indx_data                 ! indices for a local HRU
  ! output
  logical(lgt),intent(out)               :: stateMask(:)              ! mask defining desired state variables
  type(out_type_stateFilter),intent(out) :: out_stateFilter           ! number of selected state variables for a given split and error control
  ! local
  integer(i4b),allocatable               :: ixSubset(:)               ! list of indices in the state subset
  character(len=256)                     :: cmessage                  ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ! data structures
  associate(&
    ! indices for splitting methods
    ixCoupling        => in_stateFilter %  ixCoupling              ,& ! intent(in): [i4b] index of coupling method (1,2)
    ixSolution        => in_stateFilter %  ixSolution              ,& ! intent(in): [i4b] index of solution method (1,2)
    ixStateThenDomain => in_stateFilter %  ixStateThenDomain       ,& ! intent(in): [i4b] switch between full domain and sub domains
    iStateTypeSplit   => in_stateFilter %  iStateTypeSplit         ,& ! intent(in): [i4b] index of the state type split
    iDomainSplit      => in_stateFilter %  iDomainSplit            ,& ! intent(in): [i4b] index of the domain split
    iStateSplit       => in_stateFilter %  iStateSplit             ,& ! intent(in): [i4b] index of the layer split
    ! indices of model state variables
    ixStateType  => indx_data%var(iLookINDEX%ixStateType)%dat      ,& ! intent(in): [i4b(:)] indices defining the type of the state (ixNrgState...)
    ixNrgCanair  => indx_data%var(iLookINDEX%ixNrgCanair)%dat      ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
    ixNrgCanopy  => indx_data%var(iLookINDEX%ixNrgCanopy)%dat      ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
    ixHydCanopy  => indx_data%var(iLookINDEX%ixHydCanopy)%dat      ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
    ixNrgLayer   => indx_data%var(iLookINDEX%ixNrgLayer)%dat       ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
    ixHydLayer   => indx_data%var(iLookINDEX%ixHydLayer)%dat       ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
    ixWatAquifer => indx_data%var(iLookINDEX%ixWatAquifer)%dat     ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for water storage in the aquifer
    ixAllState   => indx_data%var(iLookINDEX%ixAllState)%dat       ,& ! intent(in): [i4b(:)] list of indices for all model state variables (1,2,3,...nState)
    ! number of layers
    nSnow        => indx_data%var(iLookINDEX%nSnow)%dat(1)         ,& ! intent(in): [i4b]    number of snow layers
    nSoil        => indx_data%var(iLookINDEX%nSoil)%dat(1)         ,& ! intent(in): [i4b]    number of soil layers
    nLayers      => indx_data%var(iLookINDEX%nLayers)%dat(1)       ,& ! intent(in): [i4b]    total number of layers
    ! output
    nSubset      => out_stateFilter % nSubset                      ,& ! intent(out): number of selected state variables for a given split
    err          => out_stateFilter % err                          ,& ! intent(out): error code
    message      => out_stateFilter % cmessage                       & ! intent(out): error message
    ) ! data structures
    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='stateFilter/'

    ! identify splitting option
    select case(ixCoupling)

      ! -----
      ! - fully coupled...
      ! ------------------

      ! use all state variables
      case(fullyCoupled); stateMask(:) = .true.

      ! -----
      ! - splitting by state type...
      ! ----------------------------

      ! initial split by state type
      case(stateTypeSplit)

        ! switch between full domain and sub domains
        select case(ixStateThenDomain)

          ! split into energy and mass
          case(fullDomain)
            select case(iStateTypeSplit)
              case(nrgSplit);  stateMask = (ixStateType==iname_nrgCanair .or. ixStateType==iname_nrgCanopy .or. ixStateType==iname_nrgLayer)
              case(massSplit); stateMask = (ixStateType==iname_liqCanopy .or. ixStateType==iname_liqLayer  .or. ixStateType==iname_lmpLayer .or. ixStateType==iname_watAquifer)
              case default; err=20; message=trim(message)//'unable to identify split based on state type'; return
            end select

          ! split into vegetation, snow, and soil
          case(subDomain)

            ! define state mask
            stateMask=.false. ! (initialize state mask)
            select case(iStateTypeSplit)

              ! define mask for energy
              case(nrgSplit)
                select case(iDomainSplit)
                  case(vegSplit)
                    if(ixNrgCanair(1)/=integerMissing) stateMask(ixNrgCanair) = .true.  ! energy of the canopy air space
                    if(ixNrgCanopy(1)/=integerMissing) stateMask(ixNrgCanopy) = .true.  ! energy of the vegetation canopy
                    stateMask(ixNrgLayer(1)) = .true.  ! energy of the upper-most layer in the snow+soil domain
                  case(snowSplit);   if(nSnow>1) stateMask(ixNrgLayer(2:nSnow)) = .true.    ! NOTE: (2:) because the top layer in the snow+soil domain included in vegSplit
                  case(soilSplit);   stateMask(ixNrgLayer(max(2,nSnow+1):nLayers)) = .true. ! NOTE: max(2,nSnow+1) gives second layer unless more than 2 snow layers
                  case(aquiferSplit) ! do nothing: no energy state variable for the aquifer domain
                  case default; err=20; message=trim(message)//'unable to identify model sub-domain'; return
                end select

              ! define mask for water
              case(massSplit)
                select case(iDomainSplit)
                  case(vegSplit);     if(ixHydCanopy(1)/=integerMissing) stateMask(ixHydCanopy) = .true.  ! hydrology of the vegetation canopy
                  case(snowSplit);    stateMask(ixHydLayer(1:nSnow)) = .true.  ! snow hydrology
                  case(soilSplit);    stateMask(ixHydLayer(nSnow+1:nLayers)) = .true.  ! soil hydrology
                  case(aquiferSplit); if(ixWatAquifer(1)/=integerMissing) stateMask(ixWatAquifer) = .true.  ! aquifer storage
                  case default; err=20; message=trim(message)//'unable to identify model sub-domain'; return
                end select

              ! check
              case default; err=20; message=trim(message)//'unable to identify the state type'; return
            end select  ! (split based on state type)

          ! check
          case default; err=20; message=trim(message)//'unable to identify the switch between full domains and sub domains'; return
        end select ! (switch between full domains and sub domains)

        ! check
      case default; err=20; message=trim(message)//'unable to identify coupling method'; return
    end select  ! (selecting solution method)

    ! identify scalar solutions
    if(ixSolution==scalar)then

      ! get the subset of indices
      call indxSubset(ixSubset, ixAllState, stateMask, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! get the mask
      stateMask(:) = .false.
      stateMask( ixSubset(iStateSplit) ) = .true.

      ! check
      if(count(stateMask)/=1)then
        message=trim(message)//'expect size=1 (scalar)'
        err=20; return
      endif

    endif

    ! get the number of selected state variables
    nSubset = count(stateMask)

    ! end associations
  end associate

end subroutine stateFilter

end module opSplittin_module
