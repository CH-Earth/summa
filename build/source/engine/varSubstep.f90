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

module varSubstep_module

! data types
USE nrtype

! access missing values
USE multiconst,only:integerMissing  ! missing integer
USE multiconst,only:realMissing     ! missing double precision number
USE multiconst,only:quadMissing     ! missing quadruple precision number

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! safety: set private unless specified otherwise
implicit none
private
public::varSubstep

! algorithmic parameters
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers

contains


 ! **********************************************************************************************************
 ! public subroutine varSubstep: run the model for a collection of substeps for a given state subset
 ! **********************************************************************************************************
 subroutine varSubstep(&
                       ! input: model control
                       dt,                & ! intent(in)    : time step (s)
                       dt_min,            & ! intent(in)    : minimum time step (seconds)
                       errTol,            & ! intent(in)    : error tolerance for the explicit solution
                       nState,            & ! intent(in)    : total number of state variables
                       doAdjustTemp,      & ! intent(in)    : flag to indicate if we adjust the temperature
                       firstSubStep,      & ! intent(in)    : flag to denote first sub-step
                       firstFluxCall,     & ! intent(inout) : flag to indicate if we are processing the first flux call
                       explicitEuler,     & ! intent(in)    : flag to denote computing the explicit Euler solution
                       computeVegFlux,    & ! intent(in)    : flag to denote if computing energy flux over vegetation
                       fluxMask,          & ! intent(in)    : mask for the fluxes used in this given state subset
                       ! input/output: data structures
                       model_decisions,   & ! intent(in)    : model decisions
                       type_data,         & ! intent(in)    : type of vegetation and soil
                       attr_data,         & ! intent(in)    : spatial attributes
                       forc_data,         & ! intent(in)    : model forcing data
                       mpar_data,         & ! intent(in)    : model parameters
                       indx_data,         & ! intent(inout) : index data
                       prog_data,         & ! intent(inout) : model prognostic variables for a local HRU
                       diag_data,         & ! intent(inout) : model diagnostic variables for a local HRU
                       flux_data,         & ! intent(inout) : model fluxes for a local HRU
                       deriv_data,        & ! intent(inout) : derivatives in model fluxes w.r.t. relevant state variables
                       bvar_data,         & ! intent(in)    : model variables for the local basin
                       ! output: error control
                       dtMultiplier,      & ! intent(out)   : substep multiplier (-)
                       nSubsteps,         & ! intent(out)   : number of substeps taken for a given split
                       failedMinimumStep, & ! intent(out)   : flag to denote success of substepping for a given split
                       err,message)         ! intent(out)   : error code and error message
 ! ---------------------------------------------------------------------------------------
 ! structure allocations
 USE globalData,only:flux_meta                        ! metadata on the model fluxes
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE systemSolv_module,only:systemSolv                ! solve the system of equations for one time step
 USE getVectorz_module,only:popStateVec               ! populate the state vector
 USE getVectorz_module,only:varExtract                ! extract variables from the state vector
 USE updateVars_module,only:updateVars                ! update prognostic variables
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                            ! time step (seconds)
 real(dp),intent(in)             :: dt_min                        ! minimum time step (seconds)
 real(dp),intent(in)             :: errTol                        ! error tolerance in the explicit solution
 integer(i4b),intent(in)         :: nState                        ! total number of state variables
 logical(lgt),intent(in)         :: doAdjustTemp                  ! flag to indicate if we adjust the temperature
 logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall                 ! flag to define the first flux call
 logical(lgt),intent(in)         :: explicitEuler                 ! flag to denote computing the explicit Euler solution
 logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 logical(lgt),intent(in)         :: fluxMask(:)                   ! flags to denote if the flux is calculated in the given state subset
 ! input/output: data structures
 type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
 type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                     ! spatial attributes
 type(var_d),intent(in)          :: forc_data                     ! model forcing data
 type(var_d),intent(in)          :: mpar_data                     ! model parameters
 type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
 type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                     ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
 type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
 ! output: error control
 real(dp),intent(out)            :: dtMultiplier                  ! substep multiplier (-)
 integer(i4b),intent(out)        :: nSubsteps                     ! number of substeps taken for a given split
 logical(lgt),intent(out)        :: failedMinimumStep             ! flag to denote success of substepping for a given split
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 ! error control
 character(LEN=256)              :: cmessage                      ! error message of downwind routine
 ! time stepping
 real(dp)                        :: dtSum                         ! sum of time from successful steps (seconds)
 real(dp)                        :: dt_wght                       ! weight given to a given flux calculation
 real(dp)                        :: dtSubstep                     ! length of a substep (s) 
 integer(i4b)                    :: iVar                          ! index of variables in data structures
 ! adaptive sub-stepping for the explicit solution
 logical(lgt)                    :: failedSubstep                 ! flag to denote success of substepping for a given split
 real(dp)                        :: explicitError                 ! error in the explicit solution
 real(dp),parameter              :: safety=0.85_dp                ! safety factor in adaptive sub-stepping
 real(dp),parameter              :: reduceMin=0.1_dp              ! mimimum factor that time step is reduced
 real(dp),parameter              :: increaseMax=4.0_dp            ! maximum factor that time step is increased
 ! adaptive sub-stepping for the implicit solution
 integer(i4b)                    :: niter                         ! number of iterations taken
 integer(i4b),parameter          :: n_inc=5                       ! minimum number of iterations to increase time step
 integer(i4b),parameter          :: n_dec=15                      ! maximum number of iterations to decrease time step
 real(dp),parameter              :: F_inc = 1.25_dp               ! factor used to increase time step
 real(dp),parameter              :: F_dec = 0.90_dp               ! factor used to decrease time step
 ! state and flux vectors
 real(dp)                        :: stateVecInit(nState)          ! initial state vector (mixed units)
 real(dp)                        :: stateVecTrial(nState)         ! trial state vector (mixed units)
 type(var_dlength)               :: flux_temp                     ! temporary model fluxes
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 globalVars: associate(&
 ! number of layers
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):    [i4b]    total number of layers
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
 ! *********************************************************************************************************************************************************
 ! Procedure starts here

 ! initialize error control
 err=20; message='varSubstep/'

 ! initialize flag for the success of the substepping
 failedMinimumStep=.false.

 ! initialize the length of the substep
 dtSubstep = dt

 ! allocate space for the temporary model flux structure
 call allocLocal(flux_meta(:),flux_temp,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! initialize the model fluxes (some model fluxes are not computed in the iterations)
 do iVar=1,size(flux_data%var)
  flux_temp%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:)
 end do

 ! initialize subStep
 dtSum     = 0._dp  ! keep track of the portion of the time step that is completed
 nSubsteps = 0

 ! loop through substeps
 ! NOTE: continuous do statement with exit clause
 substeps: do

  ! increment substep
  nSubsteps = nSubsteps+1

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
                  dtSubstep,      & ! intent(in):    time step (s)
                  nState,         & ! intent(in):    total number of state variables
                  firstSubStep,   & ! intent(in):    flag to denote first sub-step
                  firstFluxCall,  & ! intent(inout): flag to indicate if we are processing the first flux call
                  explicitEuler,  & ! intent(in):    flag to denote computing the explicit Euler solution
                  computeVegFlux, & ! intent(in):    flag to denote if computing energy flux over vegetation
                  ! input/output: data structures
                  type_data,      & ! intent(in):    type of vegetation and soil
                  attr_data,      & ! intent(in):    spatial attributes
                  forc_data,      & ! intent(in):    model forcing data
                  mpar_data,      & ! intent(in):    model parameters
                  indx_data,      & ! intent(inout): index data
                  prog_data,      & ! intent(inout): model prognostic variables for a local HRU
                  diag_data,      & ! intent(inout): model diagnostic variables for a local HRU
                  flux_temp,      & ! intent(inout): model fluxes for a local HRU
                  bvar_data,      & ! intent(in):    model variables for the local basin
                  model_decisions,& ! intent(in):    model decisions
                  stateVecInit,   & ! intent(in):    initial state vector
                  ! output: model control
                  deriv_data,     & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  stateVecTrial,  & ! intent(out):   updated state vector
                  explicitError,  & ! intent(out):   error in the explicit solution
                  niter,          & ! intent(out):   number of iterations taken
                  err,cmessage)     ! intent(out):   error code and error message
  if(err>0)then; message=trim(message)//trim(cmessage); return; endif

  ! implicit Euler
  if(.not.explicitEuler)then

   ! failure: step halving
   failedSubStep = (err<0)
   if(failedSubstep)then
    err=0 ! recover from failed convergence
    dtMultiplier  = 0.5_dp   ! step halving

   ! success: adjust step length based on iteration count
   else    
    if(niter<n_inc)then
     dtMultiplier = F_inc
    elseif(niter>n_dec)then
     dtMultiplier = F_dec
    else
     dtMultiplier = 1._dp
    endif
   endif  ! successful step

  ! explicit Euler
  else
   failedSubstep  = (explicitError > errTol)
   if(failedSubstep)then
    dtMultiplier  = max(safety*sqrt(errTol/explicitError), reduceMin)
   else
    dtMultiplier  = min(safety*sqrt(errTol/explicitError), increaseMax)
   endif

  endif  ! switch between explicit and implicit Euler

  ! check if we failed the substep
  if(failedSubstep)then

   ! adjust substep length
   dtSubstep = dtSubstep*dtMultiplier  

   ! check that the substep is not tiny
   if(dtSubstep<dt_min)then
    failedMinimumStep=.true.
    exit subSteps
  
   ! if step size if ok then try again
   else
    cycle subSteps
   endif

  endif  ! if failed the substep

  ! -----
  ! * update model fluxes...
  ! ------------------------

  ! NOTE: we get to here if iterations are successful
  if(err/=0)then
   message=trim(message)//'expect err=0 if updating fluxes'
   return
  endif

  ! increment fluxes
  dt_wght = dtSubstep/dt ! (define weight applied to each splitting operation)
  do iVar=1,size(flux_meta)
   if(fluxMask(iVar)) flux_data%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
  end do

  ! -----
  ! * extract variables and update states...
  ! ----------------------------------------

  ! NOTE: updating prog_data structures

  ! new block (to use dimension information)
  updateProg : block

  ! trial state variables
  real(dp)                        :: scalarCanairTempTrial         ! trial value for temperature of the canopy air space (K)
  real(dp)                        :: scalarCanopyTempTrial         ! trial value for temperature of the vegetation canopy (K)
  real(dp)                        :: scalarCanopyWatTrial          ! trial value for liquid water storage in the canopy (kg m-2)
  real(dp),dimension(nLayers)     :: mLayerTempTrial               ! trial vector for temperature of layers in the snow and soil domains (K)
  real(dp),dimension(nLayers)     :: mLayerVolFracWatTrial         ! trial vector for volumetric fraction of total water (-)
  real(dp),dimension(nSoil)       :: mLayerMatricHeadTrial         ! trial vector for total water matric potential (m)
  real(dp),dimension(nSoil)       :: mLayerMatricHeadLiqTrial      ! trial vector for liquid water matric potential (m)

  ! diagnostic variables
  real(dp)                        :: scalarCanopyLiqTrial          ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(dp)                        :: scalarCanopyIceTrial          ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(dp),dimension(nLayers)     :: mLayerVolFracLiqTrial         ! trial vector for volumetric fraction of liquid water (-)
  real(dp),dimension(nLayers)     :: mLayerVolFracIceTrial         ! trial vector for volumetric fraction of ice (-)

  ! extract states from the state vector
  call varExtract(&
                  ! input
                  stateVecTrial,            & ! intent(in):    model state vector (mixed units)
                  diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                  prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                & ! intent(in):    indices defining model states and layers
                  ! output: variables for the vegetation canopy
                  scalarCanairTempTrial,    & ! intent(out):   trial value of canopy air temperature (K)
                  scalarCanopyTempTrial,    & ! intent(out):   trial value of canopy temperature (K)
                  scalarCanopyWatTrial,     & ! intent(out):   trial value of canopy total water (kg m-2)
                  scalarCanopyLiqTrial,     & ! intent(out):   trial value of canopy liquid water (kg m-2)
                  scalarCanopyIceTrial,     & ! intent(out):   trial value of canopy ice content (kg m-2)
                  ! output: variables for the snow-soil domain
                  mLayerTempTrial,          & ! intent(out):   trial vector of layer temperature (K)
                  mLayerVolFracWatTrial,    & ! intent(out):   trial vector of volumetric total water content (-)
                  mLayerVolFracLiqTrial,    & ! intent(out):   trial vector of volumetric liquid water content (-)
                  mLayerVolFracIceTrial,    & ! intent(out):   trial vector of volumetric ice water content (-)
                  mLayerMatricHeadTrial,    & ! intent(out):   trial vector of total water matric potential (m)
                  mLayerMatricHeadLiqTrial, & ! intent(out):   trial vector of liquid water matric potential (m)
                  ! output: error control
                  err,cmessage)               ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

  ! update diagnostic variables
  call updateVars(&
                  ! input
                  doAdjustTemp,                         & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                  mpar_data,                            & ! intent(in):    model parameters for a local HRU
                  indx_data,                            & ! intent(in):    indices defining model states and layers
                  prog_data,                            & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                            & ! intent(inout): model diagnostic variables for a local HRU
                  deriv_data,                           & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  ! output: variables for the vegetation canopy
                  scalarCanopyTempTrial,                & ! intent(inout): trial value of canopy temperature (K)
                  scalarCanopyWatTrial,                 & ! intent(inout): trial value of canopy total water (kg m-2)
                  scalarCanopyLiqTrial,                 & ! intent(inout): trial value of canopy liquid water (kg m-2)
                  scalarCanopyIceTrial,                 & ! intent(inout): trial value of canopy ice content (kg m-2)
                  ! output: variables for the snow-soil domain
                  mLayerTempTrial,                      & ! intent(inout): trial vector of layer temperature (K)
                  mLayerVolFracWatTrial,                & ! intent(inout): trial vector of volumetric total water content (-)
                  mLayerVolFracLiqTrial,                & ! intent(inout): trial vector of volumetric liquid water content (-)
                  mLayerVolFracIceTrial,                & ! intent(inout): trial vector of volumetric ice water content (-)
                  mLayerMatricHeadTrial,                & ! intent(inout): trial vector of total water matric potential (m)
                  mLayerMatricHeadLiqTrial,             & ! intent(inout): trial vector of liquid water matric potential (m)
                  ! output: error control
                  err,cmessage)                           ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

  ! update prog structures
  ! NOTE: need to do the copy rather than pass to avoid conflicting associate statements

  ! build elements of the state vector for the vegetation canopy
  scalarCanairTemp    = scalarCanairTempTrial    ! trial value of canopy air temperature (K)
  scalarCanopyTemp    = scalarCanopyTempTrial    ! trial value of canopy temperature (K)
  scalarCanopyWat     = scalarCanopyWatTrial     ! trial value of canopy total water (kg m-2) kg m-2
  scalarCanopyLiq     = scalarCanopyLiqTrial     ! trial value of canopy liquid water (kg m-2)kg m-2
  scalarCanopyIce     = scalarCanopyIceTrial     ! trial value of canopy ice content (kg m-2) kg m-2

  ! build elements of the state vector for the snow+soil domain
  mLayerTemp          = mLayerTempTrial          ! trial vector of layer temperature (K)
  mLayerVolFracWat    = mLayerVolFracWatTrial    ! trial vector of volumetric total water content (-)
  mLayerVolFracLiq    = mLayerVolFracLiqTrial    ! trial vector of volumetric liquid water content (-)
  mLayerVolFracIce    = mLayerVolFracIceTrial    ! trial vector of volumetric ice water content (-)
  mLayerMatricHead    = mLayerMatricHeadTrial    ! trial vector of matric head (m)
  mLayerMatricHeadLiq = mLayerMatricHeadLiqTrial ! trial vector of matric head (m)

  end block updateProg

  ! ------------------------------------------------------
  ! ------------------------------------------------------

  ! increment sub-step
  dtSum = dtSum + dtSubstep

  ! check that we have completed the sub-step
  if(dtSum >= dt-verySmall)then
   failedMinimumStep=.false.
   exit subSteps
  endif

  ! adjust length of the sub-step (make sure that we don't exceed the step)
  dtSubstep = min(dt - dtSum, dtSubstep*dtMultiplier)

 end do substeps  ! time steps for variable-dependent sub-stepping

 ! end associate statements
 end associate globalVars

 end subroutine varSubstep

end module varSubstep_module
