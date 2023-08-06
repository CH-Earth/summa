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

module varSubstep_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas       ! named variables for the canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! global metadata
USE globalData,only:flux_meta       ! metadata on the model fluxes

! derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_flagVec,  & ! data vector with variable length dimension (i4b)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    zLookup,      & ! lookup tables
                    model_options   ! defines the model decisions

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! look up structure for variable types
USE var_lookup,only:iLookVarType

! constants
USE multiconst,only:&
                    Tfreeze,      & ! freezing temperature                 (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! look-up values for the numerical method
USE mDecisions_module,only:       &
                    numrec       ,&  ! home-grown backward Euler solution using free versions of Numerical recipes
                    kinsol       ,&  ! SUNDIALS backward Euler solution using Kinsol
                    ida              ! SUNDIALS solution using IDA

! safety: set private unless specified otherwise
implicit none
private
public::varSubstep

! algorithmic parameters
real(rkind),parameter     :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers

contains


! **********************************************************************************************************
! public subroutine varSubstep: run the model for a collection of substeps for a given state subset
! **********************************************************************************************************
subroutine varSubstep(&
                      ! input: model control
                      dt,                & ! intent(in)    : time step (s)
                      dtInit,            & ! intent(in)    : initial time step (seconds)
                      dt_min,            & ! intent(in)    : minimum time step (seconds)
                      whole_step,        & ! intent(in)    : length of whole step for surface drainage and average flux
                      nState,            & ! intent(in)    : total number of state variables
                      doAdjustTemp,      & ! intent(in)    : flag to indicate if we adjust the temperature
                      firstSubStep,      & ! intent(in)    : flag to denote first sub-step
                      firstFluxCall,     & ! intent(inout) : flag to indicate if we are processing the first flux call
                      computeVegFlux,    & ! intent(in)    : flag to denote if computing energy flux over vegetation
                      scalarSolution,    & ! intent(in)    : flag to denote implementing the scalar solution
                      iStateSplit,       & ! intent(in)    : index of the state in the splitting operation
                      fluxMask,          & ! intent(in)    : mask for the fluxes used in this given state subset
                      fluxCount,         & ! intent(inout) : number of times that fluxes are updated (should equal nSubsteps)
                       ! input/output: data structures
                      model_decisions,   & ! intent(in)    : model decisions
                      lookup_data,       & ! intent(in)    : lookup tables
                      type_data,         & ! intent(in)    : type of vegetation and soil
                      attr_data,         & ! intent(in)    : spatial attributes
                      forc_data,         & ! intent(in)    : model forcing data
                      mpar_data,         & ! intent(in)    : model parameters
                      indx_data,         & ! intent(inout) : index data
                      prog_data,         & ! intent(inout) : model prognostic variables for a local HRU
                      diag_data,         & ! intent(inout) : model diagnostic variables for a local HRU
                      flux_data,         & ! intent(inout) : model fluxes for a local HRU
                      flux_mean,         & ! intent(inout) : mean model fluxes for a local HRU
                      deriv_data,        & ! intent(inout) : derivatives in model fluxes w.r.t. relevant state variables
                      bvar_data,         & ! intent(in)    : model variables for the local basin
                      ! output: model control
                      ixSaturation,      & ! intent(inout) : index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      dtMultiplier,      & ! intent(out)   : substep multiplier (-)
                      nSubsteps,         & ! intent(out)   : number of substeps taken for a given split
                      failedMinimumStep, & ! intent(out)   : flag to denote success of substepping for a given split
                      reduceCoupledStep, & ! intent(out)   : flag to denote need to reduce the length of the coupled step
                      tooMuchMelt,       & ! intent(out)   : flag to denote that ice is insufficient to support melt
                      err,message)         ! intent(out)   : error code and error message
  ! ---------------------------------------------------------------------------------------
  ! structure allocations
  USE allocspace_module,only:allocLocal                ! allocate local data structures
  ! simulation of fluxes and residuals given a trial state vector
  USE getVectorz_module,only:popStateVec                ! populate the state vector
  USE getVectorz_module,only:varExtract                 ! extract variables from the state vector
  USE systemSolv_module,only:systemSolv                 ! solve the system of equations for one time step
  ! identify name of variable type (for error message)
  USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)             :: dt                            ! time step (seconds)
  real(rkind),intent(in)             :: dtInit                        ! initial time step (seconds)
  real(rkind),intent(in)             :: dt_min                        ! minimum time step (seconds)
  real(rkind),intent(in)             :: whole_step                    ! length of whole step for surface drainage and average flux
  integer(i4b),intent(in)            :: nState                        ! total number of state variables
  logical(lgt),intent(in)            :: doAdjustTemp                  ! flag to indicate if we adjust the temperature
  logical(lgt),intent(in)            :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)         :: firstFluxCall                 ! flag to define the first flux call
  logical(lgt),intent(in)            :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  logical(lgt),intent(in)            :: scalarSolution                ! flag to denote implementing the scalar solution
  integer(i4b),intent(in)            :: iStateSplit                   ! index of the state in the splitting operation
  type(var_flagVec),intent(in)       :: fluxMask                      ! flags to denote if the flux is calculated in the given state subset
  type(var_ilength),intent(inout)    :: fluxCount                     ! number of times that the flux is updated (should equal nSubsteps)
  ! input/output: data structures
  type(model_options),intent(in)     :: model_decisions(:)            ! model decisions
  type(zLookup),intent(in)           :: lookup_data                   ! lookup tables
  type(var_i),intent(in)             :: type_data                     ! type of vegetation and soil
  type(var_d),intent(in)             :: attr_data                     ! spatial attributes
  type(var_d),intent(in)             :: forc_data                     ! model forcing data
  type(var_dlength),intent(in)       :: mpar_data                     ! model parameters
  type(var_ilength),intent(inout)    :: indx_data                     ! indices for a local HRU
  type(var_dlength),intent(inout)    :: prog_data                     ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)    :: diag_data                     ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)    :: flux_data                     ! model fluxes for a local HRU
  type(var_dlength),intent(inout)    :: flux_mean                     ! mean model fluxes for a local HRU
  type(var_dlength),intent(inout)    :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  type(var_dlength),intent(in)       :: bvar_data                     ! model variables for the local basin
  ! output: model control
  integer(i4b),intent(inout)         :: ixSaturation                  ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  real(rkind),intent(out)            :: dtMultiplier                  ! substep multiplier (-)
  integer(i4b),intent(out)           :: nSubsteps                     ! number of substeps taken for a given split
  logical(lgt),intent(out)           :: failedMinimumStep             ! flag to denote success of substepping for a given split
  logical(lgt),intent(out)           :: reduceCoupledStep             ! flag to denote need to reduce the length of the coupled step
  logical(lgt),intent(out)           :: tooMuchMelt                   ! flag to denote that ice is insufficient to support melt
  integer(i4b),intent(out)           :: err                           ! error code
  character(*),intent(out)           :: message                       ! error message
  ! ---------------------------------------------------------------------------------------
  ! * general local variables
  ! ---------------------------------------------------------------------------------------
  ! error control
  character(LEN=256)                 :: cmessage                      ! error message of downwind routine
  ! general local variables
  integer(i4b)                       :: iVar                          ! index of variables in data structures
  integer(i4b)                       :: iSoil                         ! index of soil layers
  integer(i4b)                       :: ixLayer                       ! index in a given domain
  integer(i4b), dimension(1)         :: ixMin,ixMax                   ! bounds of a given flux vector
  ! time stepping
  real(rkind)                        :: dtSum                         ! sum of time from successful steps (seconds)
  real(rkind)                        :: dt_wght                       ! weight given to a given flux calculation
  real(rkind)                        :: dtSubstep                     ! length of a substep (s)
  real(rkind)                        :: maxstep                       ! maximum time step length (seconds)
  ! adaptive sub-stepping for the explicit solution
  logical(lgt)                       :: failedSubstep                 ! flag to denote success of substepping for a given split
  real(rkind),parameter              :: safety=0.85_rkind             ! safety factor in adaptive sub-stepping
  real(rkind),parameter              :: reduceMin=0.1_rkind           ! mimimum factor that time step is reduced
  real(rkind),parameter              :: increaseMax=4.0_rkind         ! maximum factor that time step is increased
  ! adaptive sub-stepping for the implicit solution
  integer(i4b)                       :: niter                         ! number of iterations taken
  integer(i4b),parameter             :: n_inc=5                       ! minimum number of iterations to increase time step
  integer(i4b),parameter             :: n_dec=15                      ! maximum number of iterations to decrease time step
  real(rkind),parameter              :: F_inc = 1.25_rkind            ! factor used to increase time step
  real(rkind),parameter              :: F_dec = 0.90_rkind            ! factor used to decrease time step
  ! state and flux vectors
  real(rkind)                        :: untappedMelt(nState)          ! un-tapped melt energy (J m-3 s-1)
  real(rkind)                        :: stateVecInit(nState)          ! initial state vector (mixed units)
  real(rkind)                        :: stateVecTrial(nState)         ! trial state vector (mixed units)
  real(rkind)                        :: stateVecPrime(nState)         ! trial state vector (mixed units)
  type(var_dlength)                  :: flux_temp                     ! temporary model fluxes
  ! flags
  logical(lgt)                       :: firstSplitOper                ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt)                       :: checkMassBalance              ! flag to check the mass balance
  logical(lgt)                       :: checkNrgBalance
  logical(lgt)                       :: waterBalanceError             ! flag to denote that there is a water balance error
  logical(lgt)                       :: nrgFluxModified               ! flag to denote that the energy fluxes were modified
  ! energy fluxes
  real(rkind)                        :: sumCanopyEvaporation          ! sum of canopy evaporation/condensation (kg m-2 s-1)
  real(rkind)                        :: sumLatHeatCanopyEvap          ! sum of latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
  real(rkind)                        :: sumSenHeatCanopy              ! sum of sensible heat flux from the canopy to the canopy air space (W m-2)
  real(rkind)                        :: sumSoilCompress               ! sum of total soil compression
  real(rkind),allocatable            :: sumLayerCompress(:)           ! sum of soil compression by layer
  ! ---------------------------------------------------------------------------------------
  ! point to variables in the data structures
  ! ---------------------------------------------------------------------------------------
  globalVars: associate(&
    ! model decisions
    ixNumericalMethod       => model_decisions(iLookDECISIONS%num_method)%iDecision   ,& ! intent(in):    [i4b]    choice of numerical solver
    ! number of layers
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):    [i4b]    total number of layers
    nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology variables in the soil domain
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! mapping between state vectors and control volumes
    ixLayerActive           => indx_data%var(iLookINDEX%ixLayerActive)%dat            ,& ! intent(in):    [i4b(:)] list of indices for all active layers (inactive=integerMissing)
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):    [i4b(:)] mapping of full state vector to the state subset
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):    [i4b(:)] index of control volume for different domains (veg, snow, soil)
    ! model state variables (vegetation canopy)
    scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature of the vegetation canopy (K)
    scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(inout): [dp]     mass of ice on the vegetation canopy (kg m-2)
    scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(inout): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(inout): [dp]     mass of total water on the vegetation canopy (kg m-2)
    ! model state variables (snow and soil domains)
    mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(inout): [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of ice (-)
    mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of liquid water (-)
    mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of total water (-)
    mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(inout): [dp(:)]  matric head (m)
    mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat        & ! intent(inout): [dp(:)]  matric potential of liquid water (m)
    )  ! end association with variables in the data structures
    ! *********************************************************************************************************************************************************
   
    ! initialize error control
    err=0; message='varSubstep/'

    ! initialize flag for the success of the substepping
    failedMinimumStep=.false.

    ! initialize the length of the substep
    dtSubstep = dtInit

    ! change maxstep with hard code here to make only the newton step loop in systemSolv* happen more frequently
    !   NOTE: this may just be amplifying the splitting error if maxstep is smaller than the full possible step
    maxstep = mpar_data%var(iLookPARAM%maxstep)%dat(1)  ! maximum time step (s)

    ! initalize flag for checking if energy fluxes had been modified
    nrgFluxModified = .false.

    ! allocate space for the temporary model flux structure
    call allocLocal(flux_meta(:),flux_temp,nSnow,nSoil,err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! initialize the model fluxes (some model fluxes are not computed in the iterations)
    do iVar=1,size(flux_data%var)
      flux_temp%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:)
    end do

    ! initialize the total energy fluxes (modified in updateProg)
    sumCanopyEvaporation = 0._rkind  ! canopy evaporation/condensation (kg m-2 s-1)
    sumLatHeatCanopyEvap = 0._rkind  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
    sumSenHeatCanopy     = 0._rkind  ! sensible heat flux from the canopy to the canopy air space (W m-2)
    sumSoilCompress      = 0._rkind  ! total soil compression
    allocate(sumLayerCompress(nSoil)); sumLayerCompress = 0._rkind ! soil compression by layer

    ! define the first flux call in a splitting operation
    firstSplitOper = (.not.scalarSolution .or. iStateSplit==1)

    ! initialize subStep
    dtSum     = 0._rkind  ! keep track of the portion of the time step that is completed
    nSubsteps = 0

    ! loop through substeps
    ! NOTE: continuous do statement with exit clause
    substeps: do
      dtSubstep = min(dtSubstep,maxstep)

      ! initialize error control
      err=0; message='varSubstep/'

      ! -----
      ! * populate state vectors...
      ! ---------------------------

      ! initialize state vectors
      call popStateVec(&
                      ! input
                      nState,                           & ! intent(in):    number of desired state variables
                      prog_data,                        & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                        & ! intent(in):    model diagnostic variables for a local HRU
                      indx_data,                        & ! intent(in):    indices defining model states and layers
                      ! output
                      stateVecInit,                     & ! intent(out):   initial model state vector (mixed units)
                      err,cmessage)                       ! intent(out):   error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

      ! -----
      ! * iterative solution...
      ! -----------------------
      ! solve the system of equations for a given state subset
      call systemSolv(&
                      ! input: model control
                      dtSubstep,         & ! intent(in):    time step (s)
                      whole_step,        & ! intent(in):    entire time step (s)
                      nState,            & ! intent(in):    total number of state variables
                      firstSubStep,      & ! intent(in):    flag to denote first sub-step
                      firstFluxCall,     & ! intent(inout): flag to indicate if we are processing the first flux call
                      firstSplitOper,    & ! intent(inout): flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,    & ! intent(in):    flag to denote if computing energy flux over vegetation
                      scalarSolution,    & ! intent(in):    flag to denote if implementing the scalar solution
                      ! input/output: data structures
                      lookup_data,       & ! intent(in):    lookup tables
                      type_data,         & ! intent(in):    type of vegetation and soil
                      attr_data,         & ! intent(in):    spatial attributes
                      forc_data,         & ! intent(in):    model forcing data
                      mpar_data,         & ! intent(in):    model parameters
                      indx_data,         & ! intent(inout): index data
                      prog_data,         & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,         & ! intent(inout): model diagnostic variables for a local HRU
                      flux_temp,         & ! intent(inout): model fluxes for a local HRU
                      bvar_data,         & ! intent(in):    model variables for the local basin
                      model_decisions,   & ! intent(in):    model decisions
                      stateVecInit,      & ! intent(in):    initial state vector
                      ! output: model control
                      deriv_data,        & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ixSaturation,      & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      stateVecTrial,     & ! intent(out):   updated state vector
                      stateVecPrime,     & ! intent(out):   updated state vector if need the prime space (ida)
                      untappedMelt,      & ! intent(out):   un-tapped melt energy (J m-3 s-1)
                      niter,             & ! intent(out):   number of iterations taken (numrec)
                      reduceCoupledStep, & ! intent(out):   flag to reduce the length of the coupled step
                      tooMuchMelt,       & ! intent(out):   flag to denote that ice is insufficient to support melt
                      err,cmessage)        ! intent(out):   error code and error message
      if(err/=0)then ! (check for errors, but do not fail yet)
        message=trim(message)//trim(cmessage)
        if(err>0) return
      endif
 
      ! if too much melt or need to reduce length of the coupled step then return
      ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
      if(tooMuchMelt .or. reduceCoupledStep) return

      ! identify failure, should not happen in SUNDIALS
      failedSubstep = (err<0)

      ! check
      if(globalPrintFlag)then
        print*, 'niter, failedSubstep, dtSubstep = ', niter, failedSubstep, dtSubstep
        print*, trim(cmessage)
      endif

      ! reduce step based on failure
      if(failedSubstep)then
        err=0; message='varSubstep/'  ! recover from failed convergence
        dtMultiplier  = 0.5_rkind        ! system failure: step halving
      else
        ! ** implicit Euler: adjust step length based on iteration count
        if(niter<n_inc)then
          dtMultiplier = F_inc
        elseif(niter>n_dec)then
          dtMultiplier = F_dec
        else
          dtMultiplier = 1._rkind
        endif
      endif  ! switch between failure and success

      ! check if we failed the substep
      if(failedSubstep)then

        ! check that the substep is greater than the minimum step
        if(dtSubstep*dtMultiplier<dt_min)then
          ! --> exit, and either (1) try another solution method; or (2) reduce coupled step
          failedMinimumStep=.true.
          exit subSteps

        else ! step is still OK
          dtSubstep = dtSubstep*dtMultiplier
          cycle subSteps
        endif  ! if step is less than the minimum

      endif  ! if failed the substep

      ! -----
      ! * update model fluxes...
      ! ------------------------

      ! NOTE: if we get to here then we are accepting the step of dtSubstep
      if(err/=0)then
        message=trim(message)//'expect err=0 if updating fluxes'
        return
      endif

      ! identify the need to check the mass balance
      select case(ixNumericalMethod)
        case(ida);            checkMassBalance = .false. ! IDA balance agreement levels are controlled by set tolerances
        case(kinsol, numrec); checkMassBalance = .true.  ! (.not.scalarSolution)
      end select

      ! identify the need to check the energy balance, DOES NOT WORK YET and only check if ixHowHeatCap == enthalpyFD
      checkNrgBalance = .false.

      ! update prognostic variables
      call updateProg(dtSubstep,nSnow,nSoil,nLayers,doAdjustTemp,computeVegFlux,untappedMelt,stateVecTrial,stateVecPrime,checkMassBalance, checkNrgBalance, & ! input: model control
                      model_decisions,lookup_data,mpar_data,indx_data,flux_temp,prog_data,diag_data,deriv_data,                                             & ! input-output: data structures
                      waterBalanceError,nrgFluxModified,err,cmessage)                                                                                         ! output: flags and error control
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        if(err>0) return
      endif

      ! if water balance error then reduce the length of the coupled step
      if(waterBalanceError)then
        message=trim(message)//'water balance error'
        reduceCoupledStep=.true.
        err=-20; return
      endif

      if(globalPrintFlag)&
      print*, trim(cmessage)//': dt = ', dtSubstep

      ! recover from errors in prognostic update
      if(err<0)then

        ! modify step
        err=0  ! error recovery
        dtSubstep = dtSubstep/2._rkind

        ! check minimum: fail minimum step if there is an error in the update
        if(dtSubstep<dt_min)then
          failedMinimumStep=.true.
          exit subSteps
        ! minimum OK -- try again
        else
          cycle substeps
        endif

      endif  ! if errors in prognostic update

      ! get the total energy fluxes (modified in updateProg), have to do differently
      if(nrgFluxModified .or. indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing)then
        sumCanopyEvaporation  = sumCanopyEvaporation  + dtSubstep*flux_temp%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)  ! canopy evaporation/condensation (kg m-2 s-1)
        sumLatHeatCanopyEvap  = sumLatHeatCanopyEvap  + dtSubstep*flux_temp%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
        sumSenHeatCanopy      = sumSenHeatCanopy      + dtSubstep*flux_temp%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)      ! sensible heat flux from the canopy to the canopy air space (W m-2)
      else
        sumCanopyEvaporation  = sumCanopyEvaporation  + dtSubstep*flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)  ! canopy evaporation/condensation (kg m-2 s-1)
        sumLatHeatCanopyEvap  = sumLatHeatCanopyEvap  + dtSubstep*flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
        sumSenHeatCanopy      = sumSenHeatCanopy      + dtSubstep*flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)      ! sensible heat flux from the canopy to the canopy air space (W m-2)
      endif  ! if energy fluxes were modified

      ! get the total soil compression
      if (count(indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat/=integerMissing)>0) then
        ! scalar compression
        if(.not.scalarSolution .or. iStateSplit==nSoil)&
        sumSoilCompress = sumSoilCompress + dtSubstep*diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) ! total soil compression
        ! vector compression
        do iSoil=1,nSoil
          if(indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(iSoil)/=integerMissing)&
          sumLayerCompress(iSoil) = sumLayerCompress(iSoil) + dtSubstep*diag_data%var(iLookDIAG%mLayerCompress)%dat(iSoil) ! soil compression in layers
        end do
      endif

      ! print progress
      if(globalPrintFlag)&
      write(*,'(a,1x,3(f13.2,1x))') 'updating: dtSubstep, dtSum, dt = ', dtSubstep, dtSum, dt

     ! increment fluxes
      dt_wght = dtSubstep/dt ! define weight applied to each sub-step
      do iVar=1,size(flux_meta)
        if(count(fluxMask%var(iVar)%dat)>0) then

          ! ** no domain splitting
          if(count(ixLayerActive/=integerMissing)==nLayers)then
            flux_mean%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
            fluxCount%var(iVar)%dat(:) = fluxCount%var(iVar)%dat(:) + 1

          ! ** domain splitting
          else
            ixMin=lbound(flux_data%var(iVar)%dat)
            ixMax=ubound(flux_data%var(iVar)%dat)
            do ixLayer=ixMin(1),ixMax(1)
              if(fluxMask%var(iVar)%dat(ixLayer)) then
                ! special case of the transpiration sink from soil layers: only computed for the top soil layer
                if(iVar==iLookFlux%mLayerTranspire)then
                  if(ixLayer==1) flux_mean%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
                ! standard case
                else
                  flux_mean%var(iVar)%dat(ixLayer) = flux_mean%var(iVar)%dat(ixLayer) + flux_temp%var(iVar)%dat(ixLayer)*dt_wght
                endif
                fluxCount%var(iVar)%dat(ixLayer) = fluxCount%var(iVar)%dat(ixLayer) + 1
              endif
            end do
          endif  ! (domain splitting)

        endif   ! (if the flux is desired)
      end do  ! (loop through fluxes)

      ! increment the number of substeps
      nSubsteps = nSubsteps+1

      ! increment the sub-step legth
      dtSum = dtSum + dtSubstep

      ! check that we have completed the sub-step
      if(dtSum >= dt-verySmall)then
        failedMinimumStep=.false.
        exit subSteps
      endif

      ! adjust length of the sub-step (make sure that we don't exceed the step)
      dtSubstep = min(dt - dtSum, max(dtSubstep*dtMultiplier, dt_min) )

    end do substeps  ! time steps for variable-dependent sub-stepping
    ! NOTE: if we get to here then we are accepting then dtSum should dt

    ! save the fluxes as averages
    do iVar=1,size(flux_meta)
      if(count(fluxMask%var(iVar)%dat)>0) flux_data%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:)
    enddo

    ! save the energy fluxes as averages
    flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) = sumCanopyEvaporation /dt      ! canopy evaporation/condensation (kg m-2 s-1)
    flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) = sumLatHeatCanopyEvap /dt      ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
    flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)     = sumSenHeatCanopy     /dt      ! sensible heat flux from the canopy to the canopy air space (W m-2)

    ! save the soil compression diagnostics as averages
    diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) = sumSoilCompress/dt
    do iSoil=1,nSoil
      if(indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(iSoil)/=integerMissing)&
      diag_data%var(iLookDIAG%mLayerCompress)%dat(iSoil) = sumLayerCompress(iSoil)/dt
    end do
    deallocate(sumLayerCompress)

  ! end associate statements
  end associate globalVars

  ! update error codes
  if(failedMinimumStep)then
    err=-20 ! negative = recoverable error
    message=trim(message)//'failed minimum step'
  endif
end subroutine varSubstep


! **********************************************************************************************************
! private subroutine updateProg: update prognostic variables
! **********************************************************************************************************
subroutine updateProg(dt,nSnow,nSoil,nLayers,doAdjustTemp,computeVegFlux,untappedMelt,stateVecTrial,stateVecPrime,checkMassBalance, checkNrgBalance, & ! input: model control
                      model_decisions,lookup_data,mpar_data,indx_data,flux_data,prog_data,diag_data,deriv_data,                                      & ! input-output: data structures
                      waterBalanceError,nrgFluxModified,err,message)                                                                                   ! output: flags and error control
  USE getVectorz_module,only:varExtract                             ! extract variables from the state vector
#ifdef SUNDIALS_ACTIVE
  USE updateVarsWithPrime_module,only:updateVarsWithPrime             ! update prognostic variables
#endif
  USE updateVars_module,only:updateVars                             ! update prognostic variables
  USE t2enthalpy_module, only:t2enthalpy                            ! compute enthalpy
  implicit none
  ! model control
  real(rkind)      ,intent(in)    :: dt                             ! time step (s)
  integer(i4b)     ,intent(in)    :: nSnow                          ! number of snow layers
  integer(i4b)     ,intent(in)    :: nSoil                          ! number of soil layers
  integer(i4b)     ,intent(in)    :: nLayers                        ! total number of layers
  logical(lgt)     ,intent(in)    :: doAdjustTemp                   ! flag to indicate if we adjust the temperature
  logical(lgt)     ,intent(in)    :: computeVegFlux                 ! flag to compute the vegetation flux
  real(rkind)      ,intent(in)    :: untappedMelt(:)                ! un-tapped melt energy (J m-3 s-1)
  real(rkind)      ,intent(in)    :: stateVecTrial(:)               ! trial state vector (mixed units)
  real(rkind)      ,intent(in)    :: stateVecPrime(:)               ! trial state vector (mixed units)
  logical(lgt)     ,intent(in)    :: checkMassBalance               ! flag to check the mass balance
  logical(lgt)     ,intent(in)    :: checkNrgBalance                ! flag to check the energy balance
  ! data structures
  type(model_options),intent(in)  :: model_decisions(:)             ! model decisions
  type(zLookup),intent(in)        :: lookup_data                    ! lookup tables
  type(var_dlength),intent(in)    :: mpar_data                      ! model parameters
  type(var_ilength),intent(in)    :: indx_data                      ! indices for a local HRU
  type(var_dlength),intent(inout) :: flux_data                      ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: prog_data                      ! prognostic variables for a local HRU
  type(var_dlength),intent(inout) :: diag_data                      ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: deriv_data                     ! derivatives in model fluxes w.r.t. relevant state variables
  ! flags and error control
  logical(lgt)     ,intent(out)   :: waterBalanceError              ! flag to denote that there is a water balance error
  logical(lgt)     ,intent(out)   :: nrgFluxModified                ! flag to denote that the energy fluxes were modified
  integer(i4b)     ,intent(out)   :: err                            ! error code
  character(*)     ,intent(out)   :: message                        ! error message
  ! ==================================================================================================================
  ! general
  integer(i4b)                    :: iState                         ! index of model state variable
  integer(i4b)                    :: ixSubset                       ! index within the state subset
  integer(i4b)                    :: ixFullVector                   ! index within full state vector
  integer(i4b)                    :: ixControlIndex                 ! index within a given domain
  real(rkind)                     :: volMelt                        ! volumetric melt (kg m-3)
  real(rkind),parameter           :: verySmall=epsilon(1._rkind)*2._rkind ! a very small number (deal with precision issues)
  ! mass balance
  real(rkind)                     :: canopyBalance0,canopyBalance1  ! canopy storage at start/end of time step
  real(rkind)                     :: soilBalance0,soilBalance1      ! soil storage at start/end of time step
  real(rkind)                     :: vertFlux                       ! change in storage due to vertical fluxes
  real(rkind)                     :: tranSink,baseSink,compSink     ! change in storage due to sink terms
  real(rkind)                     :: liqError                       ! water balance error
  real(rkind)                     :: fluxNet                        ! net water fluxes (kg m-2 s-1)
  real(rkind)                     :: superflousWat                  ! superflous water used for evaporation (kg m-2 s-1)
  real(rkind)                     :: superflousNrg                  ! superflous energy that cannot be used for evaporation (W m-2 [J m-2 s-1])
  character(LEN=256)              :: cmessage                       ! error message of downwind routine
  ! trial state variables
  real(rkind)                     :: scalarCanairTempTrial          ! trial value for temperature of the canopy air space (K)
  real(rkind)                     :: scalarCanopyTempTrial          ! trial value for temperature of the vegetation canopy (K)
  real(rkind)                     :: scalarCanopyWatTrial           ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),dimension(nLayers)  :: mLayerTempTrial                ! trial vector for temperature of layers in the snow and soil domains (K)
  real(rkind),dimension(nLayers)  :: mLayerVolFracWatTrial          ! trial vector for volumetric fraction of total water (-)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadTrial          ! trial vector for total water matric potential (m)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadLiqTrial       ! trial vector for liquid water matric potential (m)
  real(rkind)                     :: scalarAquiferStorageTrial      ! trial value for storage of water in the aquifer (m)
  ! diagnostic variables
  real(rkind)                     :: scalarCanopyLiqTrial           ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind)                     :: scalarCanopyIceTrial           ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),dimension(nLayers)  :: mLayerVolFracLiqTrial          ! trial vector for volumetric fraction of liquid water (-)
  real(rkind),dimension(nLayers)  :: mLayerVolFracIceTrial          ! trial vector for volumetric fraction of ice (-)
    ! derivative of state variables
  real(rkind)                     :: scalarCanairTempPrime          ! trial value for temperature of the canopy air space (K)
  real(rkind)                     :: scalarCanopyTempPrime          ! trial value for temperature of the vegetation canopy (K)
  real(rkind)                     :: scalarCanopyWatPrime           ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),dimension(nLayers)  :: mLayerTempPrime                ! trial vector for temperature of layers in the snow and soil domains (K)
  real(rkind),dimension(nLayers)  :: mLayerVolFracWatPrime          ! trial vector for volumetric fraction of total water (-)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadPrime          ! trial vector for total water matric potential (m)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadLiqPrime       ! trial vector for liquid water matric potential (m)
  real(rkind)                     :: scalarAquiferStoragePrime      ! trial value for storage of water in the aquifer (m)
  ! diagnostic variables
  real(rkind)                     :: scalarCanopyLiqPrime           ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind)                     :: scalarCanopyIcePrime           ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),dimension(nLayers)  :: mLayerVolFracLiqPrime          ! trial vector for volumetric fraction of liquid water (-)
  real(rkind),dimension(nLayers)  :: mLayerVolFracIcePrime          ! trial vector for volumetric fraction of ice (-)
  real(rkind)                     :: scalarCanairEnthalpyTrial      ! enthalpy of the canopy air space (J m-3)
  real(rkind)                     :: scalarCanopyEnthalpyTrial      ! enthalpy of the vegetation canopy (J m-3)
  real(rkind),dimension(nLayers)  :: mLayerEnthalpyTrial            ! enthalpy of snow + soil (J m-3)
  ! enthalpy derivatives
  real(rkind)                     :: dCanEnthalpy_dTk               ! derivatives in canopy enthalpy w.r.t. temperature
  real(rkind)                     :: dCanEnthalpy_dWat              ! derivatives in canopy enthalpy w.r.t. water state
  real(rkind),dimension(nLayers)  :: dEnthalpy_dTk                  ! derivatives in layer enthalpy w.r.t. temperature
  real(rkind),dimension(nLayers)  :: dEnthalpy_dWat                 ! derivatives in layer enthalpy w.r.t. water state
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! point to flux variables in the data structure
  associate(&
    ! model decisions
    ixNumericalMethod         => model_decisions(iLookDECISIONS%num_method)%iDecision       ,& ! intent(in):  [i4b] choice of numerical solver
    ! get indices for mass balance
    ixVegHyd                  => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                  ,& ! intent(in)   : [i4b]    index of canopy hydrology state variable (mass)
    ixSoilOnlyHyd             => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat                ,& ! intent(in)   : [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! get indices for the un-tapped melt
    ixNrgOnly                 => indx_data%var(iLookINDEX%ixNrgOnly)%dat                    ,& ! intent(in)   : [i4b(:)] list of indices for all energy states
    ixDomainType              => indx_data%var(iLookINDEX%ixDomainType)%dat                 ,& ! intent(in)   : [i4b(:)] indices defining the domain of the state (iname_veg, iname_snow, iname_soil)
    ixControlVolume           => indx_data%var(iLookINDEX%ixControlVolume)%dat              ,& ! intent(in)   : [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixMapSubset2Full          => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat             ,& ! intent(in)   : [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! water fluxes
    scalarRainfall            => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)             ,& ! intent(in)   : [dp]     rainfall rate (kg m-2 s-1)
    scalarThroughfallRain     => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)      ,& ! intent(in)   : [dp]     rain reaches ground without touching the canopy (kg m-2 s-1)
    scalarCanopyEvaporation   => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)    ,& ! intent(in)   : [dp]     canopy evaporation/condensation (kg m-2 s-1)
    scalarCanopyLiqDrainage   => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)    ,& ! intent(in)   : [dp]     drainage liquid water from vegetation canopy (kg m-2 s-1)
    iLayerLiqFluxSoil         => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat             ,& ! intent(in)   : [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
    iLayerNrgFlux             => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat                 ,& ! intent(in)   :
    mLayerNrgFlux             => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat                 ,& ! intent(out)  : [dp]     net energy flux for each layer within the snow+soil domain (J m-3 s-1)
    mLayerTranspire           => flux_data%var(iLookFLUX%mLayerTranspire)%dat               ,& ! intent(in)   : [dp(:)]  transpiration loss from each soil layer (m s-1)
    mLayerBaseflow            => flux_data%var(iLookFLUX%mLayerBaseflow)%dat                ,& ! intent(in)   : [dp(:)]  baseflow from each soil layer (m s-1)
    mLayerCompress            => diag_data%var(iLookDIAG%mLayerCompress)%dat                ,& ! intent(in)   : [dp(:)]  change in storage associated with compression of the soil matrix (-)
    scalarCanopySublimation   => flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1)    ,& ! intent(in)   : [dp]     sublimation of ice from the vegetation canopy (kg m-2 s-1)
    scalarSnowSublimation     => flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)      ,& ! intent(in)   : [dp]     sublimation of ice from the snow surface (kg m-2 s-1)
    ! energy fluxes
    scalarLatHeatCanopyEvap   => flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)    ,& ! intent(in)   : [dp]     latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
    scalarSenHeatCanopy       => flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)        ,& ! intent(in)   : [dp]     sensible heat flux from the canopy to the canopy air space (W m-2)
    ! domain depth
    canopyDepth               => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)          ,& ! intent(in)   : [dp   ]  canopy depth (m)
    mLayerDepth               => prog_data%var(iLookPROG%mLayerDepth)%dat                   ,& ! intent(in)   : [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model state variables (vegetation canopy)
    scalarCanairTemp          => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)           ,& ! intent(inout): [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp          => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)           ,& ! intent(inout): [dp]     temperature of the vegetation canopy (K)
    scalarCanopyIce           => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)            ,& ! intent(inout): [dp]     mass of ice on the vegetation canopy (kg m-2)
    scalarCanopyLiq           => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)            ,& ! intent(inout): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyWat           => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)            ,& ! intent(inout): [dp]     mass of total water on the vegetation canopy (kg m-2)
    ! model state variables (snow and soil domains)
    mLayerTemp                => prog_data%var(iLookPROG%mLayerTemp)%dat                    ,& ! intent(inout): [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracIce          => prog_data%var(iLookPROG%mLayerVolFracIce)%dat              ,& ! intent(inout): [dp(:)]  volumetric fraction of ice (-)
    mLayerVolFracLiq          => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat              ,& ! intent(inout): [dp(:)]  volumetric fraction of liquid water (-)
    mLayerVolFracWat          => prog_data%var(iLookPROG%mLayerVolFracWat)%dat              ,& ! intent(inout): [dp(:)]  volumetric fraction of total water (-)
    mLayerMatricHead          => prog_data%var(iLookPROG%mLayerMatricHead)%dat              ,& ! intent(inout): [dp(:)]  matric head (m)
    mLayerMatricHeadLiq       => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat           ,& ! intent(inout): [dp(:)]  matric potential of liquid water (m)
    ! enthalpy
    scalarCanairEnthalpy      => diag_data%var(iLookDIAG%scalarCanairEnthalpy)%dat(1)       ,& ! intent(inout): [dp]     enthalpy of the canopy air space (J m-3)
    scalarCanopyEnthalpy      => diag_data%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1)       ,& ! intent(inout): [dp]     enthalpy of the vegetation canopy (J m-3)
    mLayerEnthalpy            => diag_data%var(iLookDIAG%mLayerEnthalpy)%dat                ,& ! intent(inout): [dp(:)]  enthalpy of the snow+soil layers (J m-3)
    ! derivatives, diagnositic for enthalpy
    dTheta_dTkCanopy          => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)         ,& ! intent(in):    [dp]     derivative of volumetric liquid water content w.r.t. temperature
    dVolTot_dPsi0             => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat               ,& ! intent(in):    [dp(:)]  derivative in total water content w.r.t. total water matric potential
    mLayerdTheta_dTk          => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat            ,& ! intent(in):    [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature
    scalarFracLiqVeg          => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)           ,& ! intent(in):    [dp]     fraction of liquid water on vegetation (-)
    mLayerFracLiqSnow         => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat             ,& ! intent(in):    [dp(:)]  fraction of liquid water in each snow layer (-)
    ! model state variables (aquifer)
    scalarAquiferStorage      => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)       ,& ! intent(inout): [dp(:)]  storage of water in the aquifer (m)
    ! error tolerance
    absConvTol_liquid         => mpar_data%var(iLookPARAM%absConvTol_liquid)%dat(1)          & ! intent(in)   : [dp]     absolute convergence tolerance for vol frac liq water (-)
    ) ! associating flux variables in the data structure
    ! -------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='updateProg/'

    ! initialize water balancmLayerVolFracWatTrial error
    waterBalanceError=.false.

    ! get storage at the start of the step
    canopyBalance0 = merge(scalarCanopyLiq + scalarCanopyIce, realMissing, computeVegFlux)
    soilBalance0   = sum( (mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)  )*mLayerDepth(nSnow+1:nLayers) )

    ! -----
    ! * update states...
    ! ------------------

    ! initialize to state variable from the last update
    scalarCanairTempTrial     = scalarCanairTemp
    scalarCanopyTempTrial     = scalarCanopyTemp
    scalarCanopyWatTrial      = scalarCanopyWat
    scalarCanopyLiqTrial      = scalarCanopyLiq
    scalarCanopyIceTrial      = scalarCanopyIce
    mLayerTempTrial           = mLayerTemp
    mLayerVolFracWatTrial     = mLayerVolFracWat
    mLayerVolFracLiqTrial     = mLayerVolFracLiq
    mLayerVolFracIceTrial     = mLayerVolFracIce
    mLayerMatricHeadTrial     = mLayerMatricHead
    mLayerMatricHeadLiqTrial  = mLayerMatricHeadLiq
    scalarAquiferStorageTrial = scalarAquiferStorage

    ! extract states from the state vector
    call varExtract(&
                    ! input
                    stateVecTrial,            & ! intent(in):    model state vector (mixed units)
                    diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                    prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                    indx_data,                & ! intent(in):    indices defining model states and layers
                    ! output: variables for the vegetation canopy
                    scalarCanairTempTrial,    & ! intent(inout):   trial value of canopy air temperature (K)
                    scalarCanopyTempTrial,    & ! intent(inout):   trial value of canopy temperature (K)
                    scalarCanopyWatTrial,     & ! intent(inout):   trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,     & ! intent(inout):   trial value of canopy liquid water (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerTempTrial,          & ! intent(inout):   trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,    & ! intent(inout):   trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,    & ! intent(inout):   trial vector of volumetric liquid water content (-)
                    mLayerMatricHeadTrial,    & ! intent(inout):   trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial, & ! intent(inout):   trial vector of liquid water matric potential (m)
                    ! output: variables for the aquifer
                    scalarAquiferStorageTrial,& ! intent(inout):   trial value of storage of water in the aquifer (m)
                    ! output: error control
                    err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! initialize to state variable from the last update
    ! should all be set to previous values if splits, but for now operator splitting is not hooked up
    scalarCanairTempPrime     = realMissing
    scalarCanopyTempPrime     = realMissing
    scalarCanopyWatPrime      = realMissing
    scalarCanopyLiqPrime      = realMissing
    scalarCanopyIcePrime      = realMissing
    mLayerTempPrime           = realMissing
    mLayerVolFracWatPrime     = realMissing
    mLayerVolFracLiqPrime     = realMissing
    mLayerVolFracIcePrime     = realMissing
    mLayerMatricHeadPrime     = realMissing
    mLayerMatricHeadLiqPrime  = realMissing
    scalarAquiferStoragePrime = realMissing

    select case(ixNumericalMethod)
      case(ida)
#ifdef SUNDIALS_ACTIVE
        call varExtract(&
                  ! input
                  stateVecPrime,            & ! intent(in):    derivative of model state vector (mixed units)
                  diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                  prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                & ! intent(in):    indices defining model states and layers
                  ! output: variables for the vegetation canopy
                  scalarCanairTempPrime,    & ! intent(inout):   derivative of canopy air temperature (K)
                  scalarCanopyTempPrime,    & ! intent(inout):   derivative of canopy temperature (K)
                  scalarCanopyWatPrime,     & ! intent(inout):   derivative of canopy total water (kg m-2)
                  scalarCanopyLiqPrime,     & ! intent(inout):   derivative of canopy liquid water (kg m-2)
                  ! output: variables for the snow-soil domain
                  mLayerTempPrime,          & ! intent(inout):   derivative of layer temperature (K)
                  mLayerVolFracWatPrime,    & ! intent(inout):   derivative of volumetric total water content (-)
                  mLayerVolFracLiqPrime,    & ! intent(inout):   derivative of volumetric liquid water content (-)
                  mLayerMatricHeadPrime,    & ! intent(inout):   derivative of total water matric potential (m)
                  mLayerMatricHeadLiqPrime, & ! intent(inout):   derivative of liquid water matric potential (m)
                  ! output: variables for the aquifer
                  scalarAquiferStoragePrime,& ! intent(inout):   derivative of storage of water in the aquifer (m)
                  ! output: error control
                  err,cmessage)               ! intent(out):   error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

        ! update diagnostic variables
        call updateVarsWithPrime(&
                    ! input
                    .false.,                                   & ! intent(in):    logical flag if computing for Jacobian update
                    doAdjustTemp,                              & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                    mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                    indx_data,                                 & ! intent(in):    indices defining model states and layers
                    prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                    diag_data,                                 & ! intent(inout): model diagnostic variables for a local HRU
                    deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    ! output: variables for the vegetation canopy
                    scalarCanopyTempTrial,                     & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatTrial,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    scalarCanopyIceTrial,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                    scalarCanopyTempPrime,                     & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatPrime,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqPrime,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    scalarCanopyIcePrime,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                    mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                    mLayerTempPrime,                           & !
                    mLayerVolFracWatPrime,                     & ! intent(inout): Prime vector of volumetric total water content (-)
                    mLayerVolFracLiqPrime,                     & ! intent(inout): Prime vector of volumetric liquid water content (-)
                    mLayerVolFracIcePrime,                     & !
                    mLayerMatricHeadPrime,                     & ! intent(inout): Prime vector of total water matric potential (m)
                    mLayerMatricHeadLiqPrime,                  & ! intent(inout): Prime vector of liquid water matric potential (m)
                    ! output: error control
                    err,cmessage)                                ! intent(out):   error control
#endif
      case(kinsol, numrec)
         ! update diagnostic variables
        call updateVars(&
                 ! input
                 doAdjustTemp,             & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                 mpar_data,                & ! intent(in):    model parameters for a local HRU
                 indx_data,                & ! intent(in):    indices defining model states and layers
                 prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                 diag_data,                & ! intent(inout): model diagnostic variables for a local HRU
                 deriv_data,               & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! output: variables for the vegetation canopy
                 scalarCanopyTempTrial,    & ! intent(inout): trial value of canopy temperature (K)
                 scalarCanopyWatTrial,     & ! intent(inout): trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,     & ! intent(inout): trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,     & ! intent(inout): trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerTempTrial,          & ! intent(inout): trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,    & ! intent(inout): trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,    & ! intent(inout): trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,    & ! intent(inout): trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,    & ! intent(inout): trial vector of total water matric potential (m)
                 mLayerMatricHeadLiqTrial, & ! intent(inout): trial vector of liquid water matric potential (m)
                 ! output: error control
                 err,cmessage)               ! intent(out):   error control
    end select

    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! ----
    ! * check energy balance
    !------------------------
    if(checkNrgBalance)then
      ! compute enthalpy at t_{n+1}
      call t2enthalpy(&
                  .true.,                      & ! intent(in): logical flag to include phase change in enthalpy
                  ! input: data structures
                  diag_data,                   & ! intent(in):  model diagnostic variables for a local HRU
                  mpar_data,                   & ! intent(in):  parameter data structure
                  indx_data,                   & ! intent(in):  model indices
                  lookup_data,                 & ! intent(in):  lookup table data structure
                  ! input: state variables for the vegetation canopy
                  scalarCanairTempTrial,       & ! intent(in):  trial value of canopy air temperature (K)
                  scalarCanopyTempTrial,       & ! intent(in):  trial value of canopy temperature (K)
                  scalarCanopyWatTrial,        & ! intent(in):  trial value of canopy total water (kg m-2)
                  scalarCanopyIceTrial,        & ! intent(in):  trial value of canopy ice content (kg m-2)
                  ! input: variables for the snow-soil domain
                  mLayerTempTrial,             & ! intent(in):  trial vector of layer temperature (K)
                  mLayerVolFracWatTrial,       & ! intent(in):  trial vector of volumetric total water content (-)
                  mLayerMatricHeadTrial,       & ! intent(in):  trial vector of total water matric potential (m)
                  mLayerVolFracIceTrial,       & ! intent(in):  trial vector of volumetric fraction of ice (-)
                  ! input: pre-computed derivatives
                  dTheta_dTkCanopy,            & ! intent(in): derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
                  scalarFracLiqVeg,            & ! intent(in): fraction of canopy liquid water (-)
                  mLayerdTheta_dTk,            & ! intent(in): derivative of volumetric liquid water content w.r.t. temperature (K-1)
                  mLayerFracLiqSnow,           & ! intent(in): fraction of liquid water (-)
                  dVolTot_dPsi0,               & ! intent(in): derivative in total water content w.r.t. total water matric potential (m-1)
                  ! output: enthalpy
                  scalarCanairEnthalpyTrial,        & ! intent(out):  enthalpy of the canopy air space (J m-3)
                  scalarCanopyEnthalpyTrial,        & ! intent(out):  enthalpy of the vegetation canopy (J m-3)
                  mLayerEnthalpyTrial,             & ! intent(out):  enthalpy of each snow+soil layer (J m-3)
                  dCanEnthalpy_dTk,            & ! intent(out):  derivatives in canopy enthalpy w.r.t. temperature
                  dCanEnthalpy_dWat,           & ! intent(out):  derivatives in canopy enthalpy w.r.t. water state
                  dEnthalpy_dTk,               & ! intent(out):  derivatives in layer enthalpy w.r.t. temperature
                  dEnthalpy_dWat,              & ! intent(out):  derivatives in layer enthalpy w.r.t. water state
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif

    ! -----
    ! * check mass balance...
    ! -----------------------

    ! NOTE: should not need to do this, since mass balance is checked in the solver, and cannot do for IDA
    !   if do not check could cause problems if should modify nrgFlux
    if(checkMassBalance)then

      ! check mass balance for the canopy
      if(ixVegHyd/=integerMissing)then

        ! handle cases where fluxes empty the canopy
        fluxNet = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
        if(-fluxNet*dt > canopyBalance0)then

          ! --> first add water
          canopyBalance1 = canopyBalance0 + (scalarRainfall - scalarThroughfallRain)*dt

          ! --> next, remove canopy evaporation -- put the unsatisfied evap into sensible heat
          canopyBalance1 = canopyBalance1 + scalarCanopyEvaporation*dt
          if(canopyBalance1 < 0._rkind)then
            ! * get superfluous water and energy
            superflousWat = -canopyBalance1/dt     ! kg m-2 s-1
            superflousNrg = superflousWat*LH_vap   ! W m-2 (J m-2 s-1)
            ! * update fluxes and states
            canopyBalance1          = 0._rkind
            scalarCanopyEvaporation = scalarCanopyEvaporation + superflousWat
            scalarLatHeatCanopyEvap = scalarLatHeatCanopyEvap + superflousNrg
            scalarSenHeatCanopy     = scalarSenHeatCanopy - superflousNrg
          endif

          ! --> next, remove canopy drainage
          canopyBalance1 = canopyBalance1 -scalarCanopyLiqDrainage*dt
          if(canopyBalance1 < 0._rkind)then
            superflousWat           = -canopyBalance1/dt     ! kg m-2 s-1
            canopyBalance1          = 0._rkind
            scalarCanopyLiqDrainage = scalarCanopyLiqDrainage + superflousWat
          endif

          ! update the trial state
          scalarCanopyWatTrial = canopyBalance1

          ! set the modification flag
          nrgFluxModified = .true.

        else
          canopyBalance1  = canopyBalance0 + fluxNet*dt
          nrgFluxModified = .false.
        endif  ! cases where fluxes empty the canopy

        ! check the mass balance
        fluxNet  = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
        liqError = (canopyBalance0 + fluxNet*dt) - scalarCanopyWatTrial
        if(abs(liqError) > absConvTol_liquid*10._rkind)then  ! *10 because of precision issues
          !write(*,'(a,1x,f20.10)') 'dt = ', dt
          !write(*,'(a,1x,f20.10)') 'scalarCanopyWatTrial         = ', scalarCanopyWatTrial
          !write(*,'(a,1x,f20.10)') 'canopyBalance0               = ', canopyBalance0
          !write(*,'(a,1x,f20.10)') 'canopyBalance1               = ', canopyBalance1
          !write(*,'(a,1x,f20.10)') 'scalarRainfall*dt            = ', scalarRainfall*dt
          !write(*,'(a,1x,f20.10)') 'scalarCanopyLiqDrainage*dt   = ', scalarCanopyLiqDrainage*dt
          !write(*,'(a,1x,f20.10)') 'scalarCanopyEvaporation*dt   = ', scalarCanopyEvaporation*dt
          !write(*,'(a,1x,f20.10)') 'scalarThroughfallRain*dt     = ', scalarThroughfallRain*dt
          !write(*,'(a,1x,f20.10)') 'liqError                     = ', liqError
          waterBalanceError = .true.
          return
        endif  ! if there is a water balance error
      endif  ! if veg canopy

      ! check mass balance for soil, again not checked for IDA
      ! NOTE: fatal errors, though possible to recover using negative error codes
      if(count(ixSoilOnlyHyd/=integerMissing)==nSoil)then
        soilBalance1 = sum( (mLayerVolFracLiqTrial(nSnow+1:nLayers) + mLayerVolFracIceTrial(nSnow+1:nLayers) )*mLayerDepth(nSnow+1:nLayers) )
        vertFlux     = -(iLayerLiqFluxSoil(nSoil) - iLayerLiqFluxSoil(0))*dt           ! m s-1 --> m
        tranSink     = sum(mLayerTranspire)*dt                                         ! m s-1 --> m
        baseSink     = sum(mLayerBaseflow)*dt                                          ! m s-1 --> m
        compSink     = sum(mLayerCompress(1:nSoil) * mLayerDepth(nSnow+1:nLayers) )*dt ! m s-1 --> m
        liqError     = soilBalance1 - (soilBalance0 + vertFlux + tranSink - baseSink - compSink)
        if(abs(liqError) > absConvTol_liquid*10._rkind)then   ! *10 because of precision issues
          !write(*,'(a,1x,f20.10)') 'dt = ', dt
          !write(*,'(a,1x,f20.10)') 'soilBalance0      = ', soilBalance0
          !write(*,'(a,1x,f20.10)') 'soilBalance1      = ', soilBalance1
          !write(*,'(a,1x,f20.10)') 'vertFlux          = ', vertFlux
          !write(*,'(a,1x,f20.10)') 'tranSink          = ', tranSink
          !write(*,'(a,1x,f20.10)') 'baseSink          = ', baseSink
          !write(*,'(a,1x,f20.10)') 'compSink          = ', compSink
          !write(*,'(a,1x,f20.10)') 'liqError          = ', liqError
          waterBalanceError = .true.
          return
        endif  ! if there is a water balance error
      endif  ! if hydrology states exist in the soil domain
    endif  ! if checking the mass balance

    ! -----
    ! * remove untapped melt energy... always 0 at the moment but if use should be in solved as affects state
    ! --------------------------------

    ! only work with energy state variables
    if(size(ixNrgOnly)>0)then  ! energy state variables exist

      ! loop through energy state variables
      do iState=1,size(ixNrgOnly)

        ! get index of the control volume within the domain
        ixSubset       = ixNrgOnly(iState)             ! index within the state subset
        ixFullVector   = ixMapSubset2Full(ixSubset)    ! index within full state vector
        ixControlIndex = ixControlVolume(ixFullVector) ! index within a given domain

        ! compute volumetric melt (kg m-3)
        volMelt = dt*untappedMelt(ixSubset)/LH_fus  ! (kg m-3)

        ! update ice content
        select case( ixDomainType(ixFullVector) )
          case(iname_cas);  cycle ! do nothing, since there is no snow stored in the canopy air space
          case(iname_veg);  scalarCanopyIceTrial                        = scalarCanopyIceTrial                        - volMelt*canopyDepth  ! (kg m-2)
          case(iname_snow); mLayerVolFracIceTrial(ixControlIndex)       = mLayerVolFracIceTrial(ixControlIndex)       - volMelt/iden_ice     ! (-)
          case(iname_soil); mLayerVolFracIceTrial(ixControlIndex+nSnow) = mLayerVolFracIceTrial(ixControlIndex+nSnow) - volMelt/iden_water   ! (-)
          case default; err=20; message=trim(message)//'unable to identify domain type [remove untapped melt energy]'; return
        end select

        ! update liquid water content
        select case( ixDomainType(ixFullVector) )
          case(iname_cas);  cycle ! do nothing, since there is no snow stored in the canopy air space
          case(iname_veg);  scalarCanopyLiqTrial                        = scalarCanopyLiqTrial                        + volMelt*canopyDepth  ! (kg m-2)
          case(iname_snow); mLayerVolFracLiqTrial(ixControlIndex)       = mLayerVolFracLiqTrial(ixControlIndex)       + volMelt/iden_water   ! (-)
          case(iname_soil); mLayerVolFracLiqTrial(ixControlIndex+nSnow) = mLayerVolFracLiqTrial(ixControlIndex+nSnow) + volMelt/iden_water   ! (-)
          case default; err=20; message=trim(message)//'unable to identify domain type [remove untapped melt energy]'; return
        end select

      end do  ! looping through energy variables

      ! ========================================================================================================

      ! *** ice

      ! --> check if we removed too much water
      if(scalarCanopyIceTrial < 0._rkind  .or. any(mLayerVolFracIceTrial < 0._rkind) )then

        ! **
        ! canopy within numerical precision
        if(scalarCanopyIceTrial < 0._rkind)then

          if(scalarCanopyIceTrial > -verySmall)then
            scalarCanopyLiqTrial = scalarCanopyLiqTrial - scalarCanopyIceTrial
            scalarCanopyIceTrial = 0._rkind

          ! encountered an inconsistency: spit the dummy
          else
            print*, 'dt = ', dt
            print*, 'untappedMelt          = ', untappedMelt
            print*, 'untappedMelt*dt       = ', untappedMelt*dt
            print*, 'scalarCanopyiceTrial  = ', scalarCanopyIceTrial
            message=trim(message)//'melted more than the available water'
            err=20; return
          endif  ! (inconsistency)

        endif  ! if checking the canopy
        ! **
        ! snow+soil within numerical precision
        do iState=1,size(mLayerVolFracIceTrial)

          ! snow layer within numerical precision
          if(mLayerVolFracIceTrial(iState) < 0._rkind)then

            if(mLayerVolFracIceTrial(iState) > -verySmall)then
              mLayerVolFracLiqTrial(iState) = mLayerVolFracLiqTrial(iState) - mLayerVolFracIceTrial(iState)
              mLayerVolFracIceTrial(iState) = 0._rkind

            ! encountered an inconsistency: spit the dummy
            else
              print*, 'dt = ', dt
              print*, 'untappedMelt          = ', untappedMelt
              print*, 'untappedMelt*dt       = ', untappedMelt*dt
              print*, 'mLayerVolFracIceTrial = ', mLayerVolFracIceTrial
              message=trim(message)//'melted more than the available water'
              err=20; return
            endif  ! (inconsistency)

          endif  ! if checking a snow layer

        end do ! (looping through state variables)

      endif  ! (if we removed too much water)

      ! ========================================================================================================

      ! *** liquid water

      ! --> check if we removed too much water
      if(scalarCanopyLiqTrial < 0._rkind  .or. any(mLayerVolFracLiqTrial < 0._rkind) )then

        ! **
        ! canopy within numerical precision
        if(scalarCanopyLiqTrial < 0._rkind)then

          if(scalarCanopyLiqTrial > -verySmall)then
            scalarCanopyIceTrial = scalarCanopyIceTrial - scalarCanopyLiqTrial
            scalarCanopyLiqTrial = 0._rkind

          ! encountered an inconsistency: spit the dummy
          else
            print*, 'dt = ', dt
            print*, 'untappedMelt          = ', untappedMelt
            print*, 'untappedMelt*dt       = ', untappedMelt*dt
            print*, 'scalarCanopyLiqTrial  = ', scalarCanopyLiqTrial
            message=trim(message)//'frozen more than the available water'
            err=20; return
          endif  ! (inconsistency)
        endif  ! checking the canopy

        ! **
        ! snow+soil within numerical precision
        do iState=1,size(mLayerVolFracLiqTrial)

          ! snow layer within numerical precision
          if(mLayerVolFracLiqTrial(iState) < 0._rkind)then

            if(mLayerVolFracLiqTrial(iState) > -verySmall)then
              mLayerVolFracIceTrial(iState) = mLayerVolFracIceTrial(iState) - mLayerVolFracLiqTrial(iState)
              mLayerVolFracLiqTrial(iState) = 0._rkind

            ! encountered an inconsistency: spit the dummy
            else
              print*, 'dt = ', dt
              print*, 'untappedMelt          = ', untappedMelt
              print*, 'untappedMelt*dt       = ', untappedMelt*dt
              print*, 'mLayerVolFracLiqTrial = ', mLayerVolFracLiqTrial
              message=trim(message)//'frozen more than the available water'
              err=20; return
            endif  ! (inconsistency)

          endif  ! checking a snow layer

        end do ! (looping through state variables)

      endif  ! (if we removed too much water)

    endif  ! (if energy state variables exist)

    ! -----
    ! * update enthalpy as a diagnostic variable...
    ! --------------------------------
    scalarCanairEnthalpy = scalarCanairEnthalpyTrial
    scalarCanopyEnthalpy = scalarCanopyEnthalpyTrial
    mLayerEnthalpy = mLayerEnthalpyTrial

    ! -----
    ! * update prognostic variables...
    ! --------------------------------
    ! update state variables for the vegetation canopy
    scalarCanairTemp    = scalarCanairTempTrial    ! trial value of canopy air temperature (K)
    scalarCanopyTemp    = scalarCanopyTempTrial    ! trial value of canopy temperature (K)
    scalarCanopyWat     = scalarCanopyWatTrial     ! trial value of canopy total water (kg m-2)
    scalarCanopyLiq     = scalarCanopyLiqTrial     ! trial value of canopy liquid water (kg m-2)
    scalarCanopyIce     = scalarCanopyIceTrial     ! trial value of canopy ice content (kg m-2)

    ! update state variables for the snow+soil domain
    mLayerTemp          = mLayerTempTrial          ! trial vector of layer temperature (K)
    mLayerVolFracWat    = mLayerVolFracWatTrial    ! trial vector of volumetric total water content (-)
    mLayerVolFracLiq    = mLayerVolFracLiqTrial    ! trial vector of volumetric liquid water content (-)
    mLayerVolFracIce    = mLayerVolFracIceTrial    ! trial vector of volumetric ice water content (-)
    mLayerMatricHead    = mLayerMatricHeadTrial    ! trial vector of matric head (m)
    mLayerMatricHeadLiq = mLayerMatricHeadLiqTrial ! trial vector of matric head (m)

    ! update state variables for the aquifer
    scalarAquiferStorage = scalarAquiferStorageTrial

    ! end associations to info in the data structures
  end associate

end subroutine updateProg

end module varSubstep_module
