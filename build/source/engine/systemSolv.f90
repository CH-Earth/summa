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

module systemSolv_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE multiconst,only:integerMissing  ! missing integer
USE multiconst,only:realMissing     ! missing double precision number
USE multiconst,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:      &
 qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                      & ! a big bucket (lumped aquifer model)
 noExplicit                        ! no explicit groundwater parameterization

! safety: set private unless specified otherwise
implicit none
private
public::systemSolv

! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=1.e-12_dp        ! a very small number (used to check consistency)
real(dp),parameter  :: veryBig=1.e+20_dp          ! a very big number
real(dp),parameter  :: dx = 1.e-8_dp              ! finite difference increment

contains


 ! **********************************************************************************************************
 ! public subroutine systemSolv: run the coupled energy-mass model for one timestep
 ! **********************************************************************************************************
 subroutine systemSolv(&
                       ! input: model control
                       dt,             & ! intent(in):    time step (s)
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
                       ! output
                       deriv_data,     & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       stateVecTrial,  & ! intent(out):   updated state vector
                       explicitError,  & ! intent(out):   error in the explicit solution
                       niter,          & ! intent(out):   number of iterations taken
                       err,message)      ! intent(out):   error code and error message
 ! ---------------------------------------------------------------------------------------
 ! structure allocations
 USE globalData,only:flux_meta                        ! metadata on the model fluxes
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE eval8summa_module,only:eval8summa                ! simulation of fluxes and residuals given a trial state vector
 USE summaSolve_module,only:summaSolve                ! calculate the iteration increment, evaluate the new state, and refine if necessary
 USE getVectorz_module,only:getScaling                ! get the scaling vectors
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                            ! time step (seconds)
 integer(i4b),intent(in)         :: nState                        ! total number of state variables
 logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall                 ! flag to define the first flux call
 logical(lgt),intent(in)         :: explicitEuler                 ! flag to denote computing the explicit Euler solution
 logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                     ! spatial attributes
 type(var_d),intent(in)          :: forc_data                     ! model forcing data
 type(var_d),intent(in)          :: mpar_data                     ! model parameters
 type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
 type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_temp                     ! model fluxes for a local HRU
 type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
 type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
 real(dp),intent(in)             :: stateVecInit(:)               ! initial state vector (mixed units)
 ! output: model control
 type(var_dlength),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables 
 real(dp),intent(out)            :: stateVecTrial(:)              ! trial state vector (mixed units)
 real(dp),intent(out)            :: explicitError                 ! error in the explicit solution
 integer(i4b),intent(out)        :: niter                         ! number of iterations taken
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 character(LEN=256)              :: cmessage                      ! error message of downwind routine
 integer(i4b)                    :: iter                          ! iteration index
 integer(i4b)                    :: iVar                          ! index of variable
 integer(i4b)                    :: iLayer                        ! index of layer in the snow+soil domain
 integer(i4b)                    :: iState                        ! index of model state
 integer(i4b)                    :: nLeadDim                      ! length of the leading dimension of the Jacobian matrix (nBands or nState)
 integer(i4b)                    :: local_ixGroundwater           ! local index for groundwater representation
 real(dp),parameter              :: tempAccelerate=0.00_dp        ! factor to force initial canopy temperatures to be close to air temperature
 real(dp),parameter              :: xMinCanopyWater=0.0001_dp     ! minimum value to initialize canopy water (kg m-2)
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: forceFullMatrix=.false.       ! flag to force the use of the full Jacobian matrix
 integer(i4b)                    :: maxiter                       ! maximum number of iterations
 integer(i4b)                    :: ixMatrix                      ! form of matrix (band diagonal or full matrix)
 type(var_dlength)               :: flux_init                     ! model fluxes at the start of the time step 
 integer(i4b)                    :: ixSaturation                  ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),allocatable            :: dBaseflow_dMatric(:,:)        ! derivative in baseflow w.r.t. matric head (s-1)  ! NOTE: allocatable, since not always needed
 real(dp)                        :: stateVecNew(nState)           ! new state vector (mixed units)
 real(dp)                        :: fluxVec0(nState)              ! flux vector (mixed units)
 real(dp)                        :: fScale(nState)                ! characteristic scale of the function evaluations (mixed units)
 real(dp)                        :: xScale(nState)                ! characteristic scale of the state vector (mixed units)
 real(dp)                        :: dMat(nState)                  ! diagonal matrix (excludes flux derivatives)
 real(qp)                        :: sMul(nState)    ! NOTE: qp    ! multiplier for state vector for the residual calculations
 real(qp)                        :: rVec(nState)    ! NOTE: qp    ! residual vector
 real(dp)                        :: rAdd(nState)                  ! additional terms in the residual vector
 real(dp)                        :: fOld,fNew                     ! function values (-); NOTE: dimensionless because scaled
 logical(lgt)                    :: stateConstrained              ! flag to denote if the state was constrained in the explicit update
 logical(lgt)                    :: feasible                      ! flag to define the feasibility of the solution
 logical(lgt)                    :: converged                     ! convergence flag
 real(dp)                        :: resSinkNew(nState)            ! additional terms in the residual vector
 real(dp)                        :: fluxVecNew(nState)            ! new flux vector
 real(qp)                        :: resVecNew(nState)  ! NOTE: qp ! new residual vector
 real(dp)                        :: solutionError(nState)         ! vector of errors in the model solution
 real(dp),dimension(1)           :: errorTemp                     ! maximum error in explicit solution
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 globalVars: associate(&
 ! model decisions
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,& ! intent(in):    [i4b]    groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& ! intent(in):    [i4b]    spatial representation of groundwater (local-column or single-basin)
 ! accelerate solutuion for temperature
 airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,& ! intent(in):    [dp]     temperature of the upper boundary of the snow and soil domains (K)
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
 ! layer geometry
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in):    [i4b]    total number of layers
 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="systemSolv/"

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! -----
 ! * initialize...
 ! ---------------

 ! define maximum number of iterations
 maxiter = nint(mpar_data%var(iLookPARAM%maxiter))

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation)

 ! allocate space for the model fluxes at the start of the time step
 call allocLocal(flux_meta(:),flux_init,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the baseflow derivatives
 ! NOTE: needs allocation because only used when baseflow sinks are active
 if(ixGroundwater==qbaseTopmodel)then
  allocate(dBaseflow_dMatric(nSoil,nSoil),stat=err)  ! baseflow depends on total storage in the soil column, hence on matric head in every soil layer
 else
  allocate(dBaseflow_dMatric(0,0),stat=err)          ! allocate zero-length dimnensions to avoid passing around an unallocated matrix
 end if 
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the baseflow derivatives'; return; end if

 ! identify the matrix solution method
 ! (the type of matrix used to solve the linear system A.X=B)
 if(local_ixGroundwater==qbaseTopmodel .or. forceFullMatrix)then
  nLeadDim=nState         ! length of the leading dimension
  ixMatrix=ixFullMatrix   ! named variable to denote the full Jacobian matrix
 else
  nLeadDim=nBands         ! length of the leading dimension
  ixMatrix=ixBandMatrix   ! named variable to denote the band-diagonal matrix
 endif
   
 ! initialize the model fluxes (some model fluxes are not computed in the iterations)
 do iVar=1,size(flux_temp%var)
  flux_init%var(iVar)%dat(:) = flux_temp%var(iVar)%dat(:)
 end do
   
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! *** NUMERICAL SOLUTION FOR A GIVEN SUBSTEP AND SPLIT *********************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 
 ! -----
 ! * get scaling vectors...
 ! ------------------------
 
 ! initialize the global print flag
 globalPrintFlag=.false.
 
 ! initialize state vectors
 call getScaling(&
                 ! input
                 diag_data,                        & ! intent(in):    model diagnostic variables for a local HRU
                 indx_data,                        & ! intent(in):    indices defining model states and layers
                 ! output
                 fScale,                           & ! intent(out):   function scaling vector (mixed units)
                 xScale,                           & ! intent(out):   variable scaling vector (mixed units)
                 sMul,                             & ! intent(out):   multiplier for state vector (used in the residual calculations)
                 dMat,                             & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes) 
                 err,cmessage)                       ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
 
 ! -----
 ! * compute the initial function evaluation...
 ! --------------------------------------------
 
 ! initialize the trial state vectors
 stateVecTrial = stateVecInit
 
 ! need to intialize canopy water at a positive value
 if(ixVegHyd/=integerMissing)then
  if(stateVecTrial(ixVegHyd) < xMinCanopyWater) stateVecTrial(ixVegHyd) = stateVecTrial(ixVegHyd) + xMinCanopyWater
 endif
 
 ! try to accelerate solution for energy
 if(ixCasNrg/=integerMissing) stateVecTrial(ixCasNrg) = stateVecInit(ixCasNrg) + (airtemp - stateVecInit(ixCasNrg))*tempAccelerate
 if(ixVegNrg/=integerMissing) stateVecTrial(ixVegNrg) = stateVecInit(ixVegNrg) + (airtemp - stateVecInit(ixVegNrg))*tempAccelerate
 
 ! compute the flux and the residual vector for a given state vector
 ! NOTE 1: The derivatives computed in eval8summa are used to calculate the Jacobian matrix for the first iteration
 ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the first iteration increment
 call eval8summa(&
                 ! input: model control
                 dt,                      & ! intent(in):    length of the time step (seconds)
                 nSnow,                   & ! intent(in):    number of snow layers
                 nSoil,                   & ! intent(in):    number of soil layers
                 nLayers,                 & ! intent(in):    number of layers
                 nState,                  & ! intent(in):    number of state variables in the current subset
                 firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                 firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                 .true.,                  & ! intent(in):    flag to indicate if we are processing the first iteration in a splitting operation
                 computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                 ! input: state vectors
                 stateVecTrial,           & ! intent(in):    model state vector
                 fScale,                  & ! intent(in):    function scaling vector
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
                 flux_init,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                 deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! input-output: baseflow
                 ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                 dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                 ! output
                 feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                 fluxVec0,                & ! intent(out):   flux vector
                 rAdd,                    & ! intent(out):   additional (sink) terms on the RHS of the state equation
                 rVec,                    & ! intent(out):   residual vector
                 fOld,                    & ! intent(out):   function evaluation
                 err,cmessage)              ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
 
 ! check feasibility (state vector SHOULD be feasible at this point)
 if(.not.feasible)then
  message=trim(message)//'unfeasible state vector for the initial function evaluation'
  err=20; return
 endif
 
 ! copy over the initial flux structure since some model fluxes are not computed in the iterations
 do concurrent ( iVar=1:size(flux_meta) )
  flux_temp%var(iVar)%dat(:) = flux_init%var(iVar)%dat(:)
 end do
 
 ! ** if explicit Euler, then estimate state vector at the end of the time step
 if(explicitEuler)then
  call explicitUpdate(indx_data,mpar_data,prog_data,deriv_data,stateVecInit,sMul,   & ! input:  indices, parameters, derivatives, initial state vector, and state multiplier
                      fluxVec0*dt + rAdd,                                           & ! input:  right-hand-side of the state equation
                      stateVecTrial,stateConstrained,                               & ! output: trial state vector and flag to denote that it was constrained
                      err,cmessage)                                                   ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
 endif  ! if explicit Euler
 
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 
 ! **************************
 ! *** MAIN ITERATION LOOP...
 ! **************************
 
 ! iterate
 ! NOTE: this do loop is skipped in the explicitEuler solution (localMaxIter=0)
 do iter=1,maxIter
 
  ! print iteration count
  print*, '*** iter, maxiter, dt = ', iter, maxiter, dt

  ! keep track of the number of iterations
  niter = iter+1  ! +1 because xFluxResid was moved outside the iteration loop (for backwards compatibility)

  ! compute the next trial state vector
  !  1) Computes the Jacobian matrix based on derivatives from the last flux evaluation
  !  2) Computes the iteration increment based on Jacobian and residuals from the last flux evaluation
  !  3) Computes new fluxes and derivatives, new residuals, and (if necessary) refines the state vector
  ! NOTE: only returns the flux vector and function evaluation when the solution method is explicitEuler
  call summaSolve(&
                  ! input: model control
                  dt,                            & ! intent(in):    length of the time step (seconds)
                  explicitEuler,                 & ! intent(in):    logical flag to only return the flux and function evaluation
                  iter,                          & ! intent(in):    iteration index
                  nSnow,                         & ! intent(in):    number of snow layers
                  nSoil,                         & ! intent(in):    number of soil layers
                  nLayers,                       & ! intent(in):    total number of layers
                  nLeadDim,                      & ! intent(in):    length of the leading dimension of the Jacobian matrix (either nBands or nState) 
                  nState,                        & ! intent(in):    total number of state variables
                  ixMatrix,                      & ! intent(in):    type of matrix (full or band diagonal)
                  firstSubStep,                  & ! intent(in):    flag to indicate if we are processing the first sub-step
                  firstFluxCall,                 & ! intent(inout): flag to indicate if we are processing the first flux call
                  computeVegFlux,                & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  ! input: state vectors
                  stateVecTrial,                 & ! intent(in):    trial state vector
                  fScale,                        & ! intent(in):    function scaling vector
                  xScale,                        & ! intent(in):    "variable" scaling vector, i.e., for state variables
                  rVec,                          & ! intent(in):    residual vector
                  sMul,                          & ! intent(in):    state vector multiplier (used in the residual calculations)
                  dMat,                          & ! intent(inout): diagonal matrix (excludes flux derivatives)
                  fOld,                          & ! intent(in):    old function evaluation
                  ! input: data structures       
                  model_decisions,               & ! intent(in):    model decisions
                  type_data,                     & ! intent(in):    type of vegetation and soil
                  attr_data,                     & ! intent(in):    spatial attributes
                  mpar_data,                     & ! intent(in):    model parameters
                  forc_data,                     & ! intent(in):    model forcing data
                  bvar_data,                     & ! intent(in):    average model variables for the entire basin
                  prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                     & ! intent(in):    index data
                  ! input-output: data structures
                  diag_data,                     & ! intent(inout): model diagnostic variables for a local HRU
                  flux_temp,                     & ! intent(inout): model fluxes for a local HRU (temporary structure)
                  deriv_data,                    & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  ! input-output: baseflow       
                  ixSaturation,                  & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                  dBaseflow_dMatric,             & ! intent(inout): derivative in baseflow w.r.t. matric head (s-1)
                  ! output
                  stateVecNew,                   & ! intent(out):   new state vector
                  fluxVecNew,                    & ! intent(out):   new flux vector
                  resSinkNew,                    & ! intent(out):   additional (sink) terms on the RHS of the state equa
                  resVecNew,                     & ! intent(out):   new residual vector
                  fNew,                          & ! intent(out):   new function evaluation
                  converged,                     & ! intent(out):   convergence flag
                  err,cmessage)                    ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  
  ! update function evaluation, residual vector, and states
  ! NOTE 1: The derivatives computed in summaSolve are used to calculate the Jacobian matrix at the next iteration
  ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the new iteration increment
  if(.not.explicitEuler)then
   fOld          = fNew
   rVec          = resVecNew 
   stateVecTrial = stateVecNew
  endif
 
  ! print progress
  !write(*,'(a,10(f16.14,1x))') 'rVec                  = ', rVec(iJac1:iJac2)
  !write(*,'(a,10(f16.10,1x))') 'fluxVecNew            = ', fluxVecNew(iJac1:iJac2)*dt
  !write(*,'(a,10(f16.10,1x))') 'stateVecTrial         = ', stateVecTrial(iJac1:iJac2)
  !print*, 'PAUSE: check states and fluxes'; read(*,*) 
 
  ! exit iteration loop if converged
  if(converged .or. explicitEuler) exit
 
  ! check convergence
  if(iter==maxiter)then
   message=trim(message)//'failed to converge'
   err=-20; return
  endif
  !print*, 'PAUSE: iterating'; read(*,*)
 
 end do  ! iterating
 !print*, 'PAUSE: after iterations'; read(*,*)
  
 ! -----
 ! * update states...
 ! ------------------
 
 ! special case of explicit Euler
 if(explicitEuler)then

  ! estimate state vector at the end of the time step, based on the average of start-of-step and end-of step fluxes
  call explicitUpdate(indx_data,mpar_data,prog_data,deriv_data,stateVecInit,sMul,            & ! input:  indices, parameters, derivatives, initial state vector, and state multiplier
                      0.5_dp*(fluxVec0 + fluxVecNew)*dt + 0.5_dp*(rAdd + resSinkNew),        & ! input:  right-hand-side of the state equation
                      stateVecTrial,stateConstrained,                                        & ! output: trial state vector and flag to denote that it was constrained
                      err,cmessage)                                                            ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! check if the state is constrained
  if(stateConstrained)then
   message=trim(message)//'do not expect constrained state for the explicit Heun solution'
   err=20; return
  endif

  ! average start-of-step and end-of-step fluxes for explicit Euler
  do iVar=1,size(flux_temp%var)
   flux_temp%var(iVar)%dat(:) = 0.5_dp*(flux_init%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:) )
  end do
   
  ! estimate the solution error
  solutionError(:) = (fluxVec0(:)*dt + rAdd(:)) - (fluxVecNew(:)*dt + resSinkNew(:))
  errorTemp        = maxval( abs(solutionError) )
  explicitError    = errorTemp(1)

  ! print progress
  write(*,'(a,1x,10(f20.10,1x))') 'fluxVec0      = ', fluxVec0(iJac1:iJac2)
  write(*,'(a,1x,10(f20.10,1x))') 'fluxVecNew    = ', fluxVecNew(iJac1:iJac2)
  write(*,'(a,1x,10(f20.10,1x))') 'rAdd          = ', rAdd(iJac1:iJac2)
  write(*,'(a,1x,10(f20.10,1x))') 'resSinkNew    = ', resSinkNew(iJac1:iJac2)
  write(*,'(a,1x,10(f20.10,1x))') 'stateVecInit  = ', stateVecInit(iJac1:iJac2)
  write(*,'(a,1x,10(f20.10,1x))') 'stateVecTrial = ', stateVecTrial(iJac1:iJac2)
  write(*,'(a,1x,10(f20.10,1x))') 'stateVecNew   = ', stateVecNew(iJac1:iJac2)
  write(*,'(a,1x,10(f20.10,1x))') 'solutionError = ', solutionError(iJac1:iJac2)
  print*, 'dt = ', dt
  print*, 'PAUSE: checking state vector for the explicit Euler solution'; read(*,*)

 ! standard implicit solution
 else  ! switch between explicit and implicit Euler

  ! update temperatures (ensure new temperature is consistent with the fluxes)
  if(nSnowSoilNrg>0)then
   do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
    iState = ixSnowSoilNrg(iLayer)
    stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dt + resSinkNew(iState))/real(sMul(iState), dp)
   end do  ! looping through non-missing energy state variables in the snow+soil domain
  endif

  ! update volumetric water content in the snow (ensure change in state is consistent with the fluxes)
  ! NOTE: for soil water balance is constrained within the iteration loop
  if(nSnowSoilHyd>0)then
   do concurrent (iLayer=1:nSnow,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow domain)
    iState = ixSnowSoilHyd(iLayer)
    stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dt + resSinkNew(iState))
   end do  ! looping through non-missing energy state variables in the snow+soil domain
  endif

 endif  ! switch between explicit and implicit Euler
      
 ! end associate statements
 end associate globalVars

 end subroutine systemSolv

 ! **********************************************************************************************************
 ! private subroutine explicitUpdate: update the states using the explicit Euler method
 ! **********************************************************************************************************
 subroutine explicitUpdate(&
                           indx_data,             & ! intent(in)  : state indices
                           mpar_data,             & ! intent(in)  : model parameters
                           prog_data,             & ! intent(in)  : model prognostic variables
                           deriv_data,            & ! intent(in)  : state derivatives
                           stateVecInit,          & ! intent(in)  : initial state vector
                           stateVecMult,          & ! intent(in)  : state vector multiplier
                           riteHandSide,          & ! intent(in)  : right-hand-side of the state equation
                           stateVecNew,           & ! intent(out) : new state vector
                           constrained,           & ! intent(out) : flag to denote if the state was constrained
                           err,message)             ! intent(out) : error control
 USE var_lookup,only:iLookPROG                      ! named variables for structure elements
 USE var_lookup,only:iLookPARAM                     ! named variables for structure elements
 USE var_lookup,only:iLookINDEX                     ! named variables for structure elements
 USE var_lookup,only:iLookDERIV                     ! named variables for structure elements
 implicit none
 ! input
 type(var_ilength),intent(in)  :: indx_data         ! state indices
 type(var_d)      ,intent(in)  :: mpar_data         ! model parameters
 type(var_dlength),intent(in)  :: prog_data         ! model prognostic variables
 type(var_dlength),intent(in)  :: deriv_data        ! derivatives in model fluxes w.r.t. relevant state variables
 real(dp)         ,intent(in)  :: stateVecInit(:)   ! initial state vector   
 real(qp)         ,intent(in)  :: stateVecMult(:)   ! state vector multiplier
 real(dp)         ,intent(in)  :: riteHandSide(:)   ! right-hand-side of the state equation
 ! output
 real(dp)         ,intent(out) :: stateVecNew(:)    ! new state vector
 logical(lgt)     ,intent(out) :: constrained       ! flag to denote if the state was constrained
 integer(i4b)     ,intent(out) :: err               ! error code
 character(*)     ,intent(out) :: message           ! error message
 ! local variables
 integer(i4b)                  :: iState            ! state index
 integer(i4b)                  :: ixFullVector      ! index in the full state vector
 integer(i4b)                  :: ixControlIndex    ! index of the control volume for different domains (veg, snow, soil) 
 real(dp)                      :: convFactor        ! conversion factor 
 real(dp)                      :: valueMin,valueMax ! minimum and maximum state values    
 ! --------------------------------------------------------------------------------------------------------------

 ! make association with model indices defined in indexSplit
 associate(&
  theta_sat           => mpar_data%var(iLookPARAM%theta_sat),              & ! intent(in): [dp]     soil porosity (-)
  theta_res           => mpar_data%var(iLookPARAM%theta_res),              & ! intent(in): [dp]     soil residual volumetric water content (-)
  mLayerVolFracIce    => prog_data%var(iLookPROG%mLayerVolFracIce)%dat,    & ! intent(in): [dp(:)]  volumetric fraction of ice (-)
  dVolTot_dPsi0       => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat,     & ! intent(in): [dp(:)]  derivative in total water content w.r.t. total water matric potential
  ixControlVolume     => indx_data%var(iLookINDEX%ixControlVolume)%dat,    & ! intent(in): [i4b(:)] index of the control volume for different domains (veg, snow, soil)
  ixMapSubset2Full    => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat,   & ! intent(in): [i4b(:)] [state subset] list of indices of the full state vector in the state subset
  ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat, & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
  ixDomainType_subset => indx_data%var(iLookINDEX%ixDomainType_subset)%dat & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
 ) ! associations

 ! initialize error control
 err=0; message='explicitUpdate/'

 ! initialize the flag to denote if the state is constrained
 constrained=.false.

 ! loop through model states
 do iState=1,size(stateVecInit)

  ! get index of the control volume within the domain
  ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  ! get the flux-2-state conversion factor
  select case( ixStateType_subset(iState) )
   case(iname_nrgCanair,iname_nrgCanopy,iname_nrgLayer);                convFactor = 1._dp/real(stateVecMult(iState), dp)
   case(iname_watCanopy,iname_liqCanopy,iname_watLayer,iname_liqLayer); convFactor = 1._dp
   case(iname_matLayer,iname_lmpLayer);                                 convFactor = 1._dp/dVolTot_dPsi0(ixControlIndex)
   case default; err=20; message=trim(message)//'unable to identify the state type'; return
  end select

  ! update the state vector 
  stateVecNew(iState) = stateVecInit(iState) + convFactor*riteHandSide(iState)

  ! impose non-negativity constraints for the mass of water on the vegetation canopy
  if(ixStateType_subset(iState)==iname_watCanopy .or. ixStateType_subset(iState)==iname_liqCanopy)then
   if(stateVecNew(iState) < 0._dp)then
    stateVecNew(iState)=0._dp
    constrained=.true.
   endif
  endif

  ! impose minimum and maximum storage constraints for volumetric water
  if(ixStateType_subset(iState)==iname_watLayer .or. ixStateType_subset(iState)==iname_liqLayer)then
   select case( ixDomainType_subset(iState) )
    case(iname_snow)
     valueMin = 0._dp
     valueMax = merge(iden_ice/iden_water, 1._dp - mLayerVolFracIce(ixControlIndex), ixStateType_subset(iState)==iname_watLayer)
    case(iname_soil)
     valueMin = theta_res
     valueMax = theta_sat
    case default; err=20; message=trim(message)//'expect domain type to be iname_snow or iname_soil'; return
   end select
   if(stateVecNew(iState) < valueMin)then
    stateVecNew(iState)=valueMin
    constrained=.true.
   endif
   if(stateVecNew(iState) > valueMax)then
    stateVecNew(iState)=valueMax
    constrained=.true.
   endif
  endif

  ! impose below-freezing constraints for snow temperature
  if(ixDomainType_subset(iState)==iname_snow .and. ixStateType_subset(iState)==iname_nrgLayer)then
   if(stateVecNew(iState) > Tfreeze)then
    stateVecNew(iState)=Tfreeze
    constrained=.true.
   endif
  endif

 end do ! looping through states

 ! end association to the information in the data structures
 end associate

 end subroutine explicitUpdate





end module systemSolv_module
