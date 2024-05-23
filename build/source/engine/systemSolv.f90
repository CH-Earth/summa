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

module systemSolv_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

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

! global metadata
USE globalData,only:flux_meta       ! metadata on the model fluxes

! constants
USE multiconst,only:&
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,                        & ! data vector (i4b)
                    var_d,                        & ! data vector (rkind)
                    var_ilength,                  & ! data vector with variable length dimension (i4b)
                    var_dlength,                  & ! data vector with variable length dimension (rkind)
                    zLookup,                      & ! lookup tables
                    model_options,                & ! defines the model decisions
                    in_type_summaSolve4homegrown, & ! class for summaSolve4homegrown arguments
                    io_type_summaSolve4homegrown, & ! class for summaSolve4homegrown arguments
                    out_type_summaSolve4homegrown   ! class for summaSolve4homegrown arguments

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:&
                    localColumn,  & ! separate groundwater representation in each local soil column
                    singleBasin     ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
                    qbaseTopmodel,& ! TOPMODEL-ish baseflow parameterization
                    bigBucket,    & ! a big bucket (lumped aquifer model)
                    noExplicit      ! no explicit groundwater parameterization

 ! look-up values for the numerical method
USE mDecisions_module,only:&
                    homegrown    ,& ! homegrown backward Euler solution based on concepts from numerical recipes
                    kinsol       ,& ! SUNDIALS backward Euler solution using Kinsol
                    ida             ! SUNDIALS solution using IDA

! safety: set private unless specified otherwise
implicit none
private
public::systemSolv

contains


! **********************************************************************************************************
! public subroutine systemSolv: run the coupled energy-mass model for one timestep
! **********************************************************************************************************
subroutine systemSolv(&
                      ! input: model control
                      dt_cur,            & ! intent(in):    current stepsize
                      dt,                & ! intent(in):    entire time step (s)
                      nState,            & ! intent(in):    total number of state variables
                      nLayers,           & ! intent(in):    total number of layers
                      firstSubStep,      & ! intent(in):    flag to denote first sub-step
                      firstFluxCall,     & ! intent(inout): flag to indicate if we are processing the first flux call
                      firstSplitOper,    & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,    & ! intent(in):    flag to denote if computing energy flux over vegetation
                      scalarSolution,    & ! intent(in):    flag to denote if implementing the scalar solution
                      computMassBalance, & ! intent(in):    flag to compute mass balance
                      computNrgBalance,  & ! intent(in):    flag to compute energy balance
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
                      ! output
                      deriv_data,        & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ixSaturation,      & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      stateVecTrial,     & ! intent(out):   updated state vector
                      stateVecPrime,     & ! intent(out):   updated state vector if need the prime space (ida)
                      fluxVec,           & ! intent(out):   new flux vector
                      resSink,           & ! intent(out):   additional (sink) terms on the RHS of the state equa
                      resVec,            & ! intent(out):   new residual vector
                      untappedMelt,      & ! intent(out):   un-tapped melt energy (J m-3 s-1)
                      ! output: balances (only computed at this level for ida)
                      balance,           & ! intent(out):   balance of energy per state
                      ! output: model control
                      niter,             & ! intent(out):   number of iterations taken (homegrown)
                      nSteps,            & ! intent(out):   number of time steps taken in solver
                      reduceCoupledStep, & ! intent(out):   flag to reduce the length of the coupled step
                      tooMuchMelt,       & ! intent(out):   flag to denote that there was too much melt
                      err,message)         ! intent(out):   error code and error message
  ! ---------------------------------------------------------------------------------------
  ! structure allocations
  USE allocspace_module,only:allocLocal                     ! allocate local data structures
  ! state vector and solver
  USE getVectorz_module,only:getScaling                     ! get the scaling vectors
  USE enthalpyTemp_module,only:T2enthalpy_snwWat            ! convert temperature to liq+ice enthalpy for a snow layer
#ifdef SUNDIALS_ACTIVE
  USE tol4ida_module,only:popTol4ida                        ! populate tolerances
  USE eval8summaWithPrime_module,only:eval8summaWithPrime   ! get the fluxes and residuals
  USE summaSolve4ida_module,only:summaSolve4ida             ! solve DAE by IDA
  USE summaSolve4kinsol_module,only:summaSolve4kinsol       ! solve DAE by KINSOL
#endif
  USE eval8summa_module,only:eval8summa                     ! get the fluxes and residuals
  USE summaSolve4homegrown_module,only:summaSolve4homegrown ! solve DAE using homegrown solver

  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt_cur                        ! current stepsize
  real(rkind),intent(in)          :: dt                            ! entire time step for drainage pond rate
  integer(i4b),intent(in)         :: nState                        ! total number of state variables
  integer(i4b),intent(in)         :: nLayers                       ! total number of layers
  logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)      :: firstFluxCall                 ! flag to define the first flux call
  logical(lgt),intent(in)         :: firstSplitOper                ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  logical(lgt),intent(in)         :: scalarSolution                ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)         :: computMassBalance             ! flag to compute mass balance
  logical(lgt),intent(in)         :: computNrgBalance              ! flag to compute energy balance
  ! input/output: data structures
  type(zLookup),intent(in)        :: lookup_data                   ! lookup tables
  type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
  type(var_d),intent(in)          :: attr_data                     ! spatial attributes
  type(var_d),intent(in)          :: forc_data                     ! model forcing data
  type(var_dlength),intent(in)    :: mpar_data                     ! model parameters
  type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
  type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
  type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_temp                     ! model fluxes for a local HRU
  type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
  type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
  real(rkind),intent(in)          :: stateVecInit(:)               ! initial state vector (mixed units)
  ! output
  type(var_dlength),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  integer(i4b),intent(inout)      :: ixSaturation                  ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  real(rkind),intent(out)         :: stateVecTrial(:)              ! trial state vector (mixed units)
  real(rkind),intent(out)         :: stateVecPrime(:)              ! trial state vector (mixed units)
  real(rkind),intent(out)         :: fluxVec(nState)               ! flux vector (mixed units)
  real(rkind),intent(out)         :: resSink(nState)               ! additional terms in the residual vector homegrown solver
  real(qp),intent(out)            :: resVec(nState)    ! NOTE: qp  ! residual vector
  real(rkind),intent(out)         :: untappedMelt(:)               ! un-tapped melt energy (J m-3 s-1)
  ! output: balances (only computed at this level for ida)
  real(rkind),intent(out)         :: balance(nState)               ! balance per state
  ! output: model control
  integer(i4b),intent(out)        :: niter                         ! number of iterations taken
  integer(i4b),intent(out)        :: nSteps                        ! number of time steps taken in solver
  logical(lgt),intent(out)        :: reduceCoupledStep             ! flag to reduce the length of the coupled step
  logical(lgt),intent(out)        :: tooMuchMelt                   ! flag to denote that there was too much melt
  integer(i4b),intent(out)        :: err                           ! error code
  character(*),intent(out)        :: message                       ! error message
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
  real(rkind)                     :: bulkDensity                   ! bulk density of a given layer (kg m-3)
  real(rkind)                     :: volEnthalpy                   ! volumetric enthalpy of a given layer (J m-3)
  real(rkind),parameter           :: tinyStep=0.000001_rkind       ! stupidly small time step (s)
  ! ------------------------------------------------------------------------------------------------------
  ! * model solver
  ! ------------------------------------------------------------------------------------------------------
  logical(lgt),parameter          :: forceFullMatrix=.false.       ! flag to force the use of the full Jacobian matrix
  integer(i4b)                    :: ixMatrix                      ! form of matrix (band diagonal or full matrix)
  type(var_dlength)               :: flux_init                     ! model fluxes at the start of the time step
  real(rkind),allocatable         :: dBaseflow_dMatric(:,:)        ! derivative in baseflow w.r.t. matric head (s-1)  ! NOTE: allocatable, since not always needed
  real(rkind)                     :: stateVecNew(nState)           ! new state vector (mixed units)
  real(rkind)                     :: fluxVec0(nState)              ! flux vector (mixed units)
  real(rkind)                     :: dMat(nState)                  ! diagonal matrix (excludes flux derivatives)
  real(qp)                        :: sMul(nState)    ! NOTE: qp    ! multiplier for state vector for the residual calculations
  real(rkind)                     :: rAdd(nState)                  ! additional terms in the residual vector
  logical(lgt)                    :: feasible                      ! feasibility flag
  logical(lgt)                    :: sunSucceeds                   ! flag to indicate if SUNDIALS successfully solved the problem in current data step
  ! ida variables
  real(rkind)                     :: atol(nState)                  ! absolute tolerance ida
  real(rkind)                     :: rtol(nState)                  ! relative tolerance ida
  type(var_dlength)               :: flux_sum                      ! sum of fluxes model fluxes for a local HRU over a dt_cur
  real(rkind), allocatable        :: mLayerCmpress_sum(:)          ! sum of compression of the soil matrix
  ! ida solver variables outputted if use eval8summaWithPrime (not used here, just inside ida solver)
  logical(lgt)                    :: firstSplitOper0               ! flag to indicate if we are processing the first flux call in a splitting operation, changed inside eval8summaWithPrime
  real(rkind)                     :: scalarCanopyTempPrime         ! prime value for temperature of the vegetation canopy (K s-1)
  real(rkind)                     :: scalarCanopyWatPrime          ! prime value for total water content of the vegetation canopy (kg m-2 s-1)
  real(rkind)                     :: mLayerTempPrime(nLayers)      ! prime vector of temperature of each snow and soil layer (K s-1)
  real(rkind), allocatable        :: mLayerMatricHeadPrime(:)      ! prime vector of matric head of each snow and soil layer (m s-1)
  real(rkind)                     :: mLayerVolFracWatPrime(nLayers)! prime vector of volumetric total water content of each snow and soil layer (s-1)
  ! kinsol and homegrown solver variables
  real(rkind)                     :: fScale(nState)                ! characteristic scale of the function evaluations (mixed units)
  real(rkind)                     :: xScale(nState)                ! characteristic scale of the state vector (mixed units)
  real(qp)                        :: resVecNew(nState)  ! NOTE: qp ! new residual vector homegrown solver
  ! homegrown solver variables
  real(rkind)                     :: fOld,fNew                     ! function values (-); NOTE: dimensionless because scaled homegrown solver
  real(rkind)                     :: xMin,xMax                     ! state minimum and maximum (mixed units) homegrown solver
  integer(i4b)                    :: maxiter                       ! maximum number of iterations homegrown solver
  integer(i4b)                    :: localMaxIter                  ! maximum number of iterations (depends on solution type) homegrown solver
  integer(i4b), parameter         :: scalarMaxIter=100             ! maximum number of iterations for the scalar solution homegrown solver
  logical(lgt)                    :: converged                     ! convergence flag homegrown solver
  logical(lgt), parameter         :: post_massCons=.false.         ! “perfectly” conserve mass by pushing the errors into the states, turn off for now to agree with SUNDIALS
  ! class objects for call to summaSolve4homegrown
  type(in_type_summaSolve4homegrown)  :: in_SS4HG  ! object for intent(in)  summaSolve4homegrown arguments
  type(io_type_summaSolve4homegrown)  :: io_SS4HG  ! object for intent(io)  summaSolve4homegrown arguments
  type(out_type_summaSolve4homegrown) :: out_SS4HG ! object for intent(out) summaSolve4homegrown arguments
  ! flags
  logical(lgt) :: return_flag ! flag for handling systemSolv returns trigerred from internal subroutines 
  logical(lgt) :: exit_flag   ! flag for handling loop exit statements trigerred from internal subroutines 
  ! -----------------------------------------------------------------------------------------------------------

  call initialize_systemSolv; if (return_flag) return ! initialize variables and allocate arrays -- return if error

  call initial_function_evaluations; if (return_flag) return ! initial function evaluations -- return if error

  ! **************************
  ! * Solving the System
  ! **************************
  associate(ixNumericalMethod => model_decisions(iLookDECISIONS%num_method)%iDecision) ! intent(in): [i4b] choice of numerical solver
    select case(ixNumericalMethod)
      case(ida)    ! solve for general time step using IDA
        call solve_with_IDA; if (return_flag) return              ! solve using IDA -- return if error
      case(kinsol) ! solve for BE time step using KINSOL
        call solve_with_KINSOL; if (return_flag) return           ! solve using KINSOL -- return if error
      case(homegrown) ! solve for BE time step using Newton iterations
        call Newton_iterations_homegrown; if (return_flag) return ! Newton iterations using homegrown solver -- return if error
    end select
  end associate 
 
  call finalize_systemSolv ! set untapped melt to zero and deallocate arrays

contains

 subroutine initialize_systemSolv
  ! *** Initial setup operations for the systemSolv subroutine ***

  ! initialize error control
  err=0; message="systemSolv/"
  return_flag=.false. ! initialize return flag
  nSteps = 0 ! initialize number of time steps taken in solver

  ! check time step size
  if (dt_cur < tinyStep) then
    message=trim(message)//'dt is tiny'
    err=20; return_flag=.true.; return
  end if

  ! initialize the flags
  tooMuchMelt        = .false.   ! too much melt
  reduceCoupledStep  = .false.   ! need to reduce the length of the coupled step
  ! initialize balances
  balance(:) = realMissing

  associate(&
   ixSpatialGroundwater => model_decisions(iLookDECISIONS%spatial_gw)%iDecision,& ! intent(in): [i4b] spatial representation of groundwater (local-column or single-basin)
   ixGroundwater        => model_decisions(iLookDECISIONS%groundwatr)%iDecision & ! intent(in): [i4b] groundwater parameterization
   &)
  
   ! modify the groundwater representation for this single-column implementation
   select case(ixSpatialGroundwater)
    case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
    case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
    case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; 
     return_flag=.true.; return
   end select 

   call allocate_memory; if (return_flag) return 

   ! identify the matrix solution method, using the full matrix can be slow in many-layered systems
   ! (the type of matrix used to solve the linear system A.X=B)
   if (local_ixGroundwater==qbaseTopmodel .or. scalarSolution .or. forceFullMatrix .or. computeVegFlux) then
     nLeadDim=nState         ! length of the leading dimension
     ixMatrix=ixFullMatrix   ! named variable to denote the full Jacobian matrix
   else
     nLeadDim=nBands         ! length of the leading dimension
     ixMatrix=ixBandMatrix   ! named variable to denote the band-diagonal matrix
   end if
  end associate

  ! initialize the model fluxes (some model fluxes are not computed in the iterations)
  do iVar=1,size(flux_temp%var)
    flux_init%var(iVar)%dat(:) = flux_temp%var(iVar)%dat(:)
  end do

  ! initialize state vectors -- get scaling vectors
  call getScaling(diag_data,indx_data,fScale,xScale,sMul,dMat,err,cmessage)     
  if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! check for errors
 end subroutine initialize_systemSolv

 subroutine allocate_memory
  ! ** Allocate arrays used in systemSolv subroutine **
  associate(&
   nSnow             => indx_data%var(iLookINDEX%nSnow)%dat(1)              ,& ! intent(in): [i4b] number of snow layers
   nSoil             => indx_data%var(iLookINDEX%nSoil)%dat(1)              ,& ! intent(in): [i4b] number of soil layers
   ixNumericalMethod => model_decisions(iLookDECISIONS%num_method)%iDecision,& ! intent(in): [i4b] choice of numerical solver
   ixGroundwater     => model_decisions(iLookDECISIONS%groundwatr)%iDecision & ! intent(in): [i4b] groundwater parameterization
   &)
   ! allocate space for the model fluxes at the start of the time step
   call allocLocal(flux_meta(:),flux_init,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for mLayerCmpress_sum at the start of the time step
   if (ixNumericalMethod==ida) then
     allocate( mLayerCmpress_sum(nSoil) )
     allocate( mLayerMatricHeadPrime(nSoil) )
   else
     allocate( mLayerCmpress_sum(0) )     ! allocate zero-length dimensions to avoid passing around an unallocated matrix
     allocate( mLayerMatricHeadPrime(0) ) ! allocate zero-length dimensions to avoid passing around an unallocated matrix
   end if

   ! allocate space for the baseflow derivatives
   ! NOTE: needs allocation because only used when baseflow sinks are active
   if (ixGroundwater==qbaseTopmodel) then
    allocate(dBaseflow_dMatric(nSoil,nSoil),stat=err) ! baseflow depends on total storage in the soil column, hence on matric head in every soil layer
   else
    allocate(dBaseflow_dMatric(0,0),stat=err)         ! allocate zero-length dimensions to avoid passing around an unallocated matrix
   end if
   if (err/=0) then; err=20; message=trim(message)//'unable to allocate space for the baseflow derivatives'; return_flag=.true.; return; end if
  end associate

 end subroutine allocate_memory

 subroutine initial_function_evaluations
  ! ** Compute initial function evaluations **

  ! initialize the trial state vectors
  stateVecTrial = stateVecInit

  ! compute the initial flux and the residual vector, also gets values needed for the Jacobian matrix 
  associate(ixNumericalMethod => model_decisions(iLookDECISIONS%num_method)%iDecision) ! intent(in): [i4b] choice of numerical solver
   if (ixNumericalMethod==ida) then
     call initial_flux_and_residual_vectors_prime; if (return_flag) return
   else
     call initial_flux_and_residual_vectors; if (return_flag) return
   end if
  end associate

  if (.not.feasible) then; message=trim(message)//'state vector not feasible'; err=20; return_flag=.true.; return; end if

  ! copy over the initial flux structure since some model fluxes are not computed in the iterations
  do concurrent ( iVar=1:size(flux_meta) )
    flux_temp%var(iVar)%dat(:) = flux_init%var(iVar)%dat(:)
  end do

  ! check the need to merge snow layers
  associate(&
   nSnow            => indx_data%var(iLookINDEX%nSnow)%dat(1)        ,& ! intent(in): [i4b]   number of snow layers
   mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat ,& ! intent(in): [dp(:)] volumetric fraction of ice (-)
   mLayerVolFracLiq => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat ,& ! intent(in): [dp(:)] volumetric fraction of liquid water (-)
   mLayerTemp       => prog_data%var(iLookPROG%mLayerTemp)%dat       ,& ! intent(in): [dp(:)] temperature of each snow/soil layer (K)
   snowfrz_scale    => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1) & ! intent(in): [dp]    scaling parameter for the snow freezing curve (K-1)
   &)
   ! check the need to merge snow layers
   if (nSnow>0) then
     ! compute the energy required to melt the top snow layer (J m-2)
     bulkDensity = mLayerVolFracIce(1)*iden_ice + mLayerVolFracLiq(1)*iden_water
     volEnthalpy = T2enthalpy_snwWat(mLayerTemp(1),bulkDensity,snowfrz_scale)
     ! set flag and error codes for too much melt
     if (-volEnthalpy < flux_init%var(iLookFLUX%mLayerNrgFlux)%dat(1)*dt_cur) then
       tooMuchMelt = .true.
       message=trim(message)//'net flux in the top snow layer can melt all the snow in the top layer'
       err=-20; return ! negative error code to denote a warning
     end if
   end if
  end associate
 end subroutine initial_function_evaluations

 subroutine initial_flux_and_residual_vectors
  ! ** Compute initial flux and residual vectors ** 
  ! Note: prime initial values are 0 so it's fine to run the regular eval8summa with every solver choice
  associate(&
   nSnow => indx_data%var(iLookINDEX%nSnow)%dat(1),& ! intent(in): [i4b] number of snow layers
   nSoil => indx_data%var(iLookINDEX%nSoil)%dat(1) & ! intent(in): [i4b] number of soil layers
   &)
   call eval8summa(&
                    ! input: model control
                    dt_cur,                  & ! intent(in):    current stepsize
                    dt,                      & ! intent(in):    length of the entire time step (seconds) for drainage pond rate
                    nSnow,                   & ! intent(in):    number of snow layers
                    nSoil,                   & ! intent(in):    number of soil layers
                    nLayers,                 & ! intent(in):    number of layers
                    nState,                  & ! intent(in):    number of state variables in the current subset
                    .false.,                 & ! intent(in):    not inside Sundials solver
                    firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                    firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                    firstSplitOper,          & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                    computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                    scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                    ! input: state vectors
                    stateVecTrial,           & ! intent(in):    model state vector
                    fScale,                  & ! intent(in):    characteristic scale of the function evaluations
                    sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                    ! input: data structures
                    model_decisions,         & ! intent(in):    model decisions
                    lookup_data,             & ! intent(in):    lookup tables
                    type_data,               & ! intent(in):    type of vegetation and soil
                    attr_data,               & ! intent(in):    spatial attributes
                    mpar_data,               & ! intent(in):    model parameters
                    forc_data,               & ! intent(in):    model forcing data
                    bvar_data,               & ! intent(in):    average model variables for the entire basin
                    prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                    ! input-output: data structures
                    indx_data,               & ! intent(inout): index data
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
                    resVec,                  & ! intent(out):   residual vector
                    fOld,                    & ! intent(out):   function evaluation
                    err,cmessage)              ! intent(out):   error control
  end associate
  if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! check for errors
 end subroutine initial_flux_and_residual_vectors

 subroutine initial_flux_and_residual_vectors_prime
#ifdef SUNDIALS_ACTIVE
  ! ** Compute initial flux and residual vectors ** 
  ! Note: Need this extra subroutine to handle the case of enthalpy as a state variable, currently only implemented in the prime version
  !       If we implement it in the regular version, we can remove this subroutine
  associate(&
   nSnow            => indx_data%var(iLookINDEX%nSnow)%dat(1)                  , & ! intent(in):    [i4b]   number of snow layers
   nSoil            => indx_data%var(iLookINDEX%nSoil)%dat(1)                  , & ! intent(in):    [i4b]   number of soil layers
   scalarCanopyEnthalpy => diag_data%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1), & ! intent(inout): [dp]    enthalpy of the vegetation canopy (J m-2)
   scalarCanopyTemp => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)        , & ! intent(inout): [dp]    temperature of the vegetation canopy (K)
   scalarCanopyWat  => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)         , & ! intent(inout): [dp]    total water content of the vegetation canopy (kg m-2)
   mLayerTemp       => prog_data%var(iLookPROG%mLayerTemp)%dat                 , & ! intent(inout): [dp(:)] temperature of each snow/soil layer (K)
   mLayerMatricHead => prog_data%var(iLookPROG%mLayerMatricHead)%dat             & ! intent(out):   [dp(:)] matric head (m) 
   &)
   stateVecPrime(:) = 0._rkind ! prime initial values are 0
   firstSplitOper0 = firstSplitOper ! set the flag for the first split operation, do not want to reset it here
   call eval8summaWithPrime(&
                    ! input: model control
                    dt,                      & ! intent(in):    length of the entire time step (seconds) for drainage pond rate
                    nSnow,                   & ! intent(in):    number of snow layers
                    nSoil,                   & ! intent(in):    number of soil layers
                    nLayers,                 & ! intent(in):    total number of layers
                    nState,                  & ! intent(in):    total number of state variables in the current subset
                    .false.,                 & ! intent(in):    not inside Sundials solver                    
                    firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                    firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                    firstSplitOper0,         & ! intent(inout): flag to indicate if we are processing the first flux call in a splitting operation
                    computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                    scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                    ! input: state vectors   
                    stateVecTrial,           & ! intent(in):    model state vector
                    stateVecPrime,           & ! intent(in):    derivative of model state vector
                    sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                    ! input: data structures
                    model_decisions,         & ! intent(in):    model decisions
                    lookup_data,             & ! intent(in):    lookup table data structure
                    type_data,               & ! intent(in):    type of vegetation and soil
                    attr_data,               & ! intent(in):    spatial attributes
                    mpar_data,               & ! intent(in):    model parameters
                    forc_data,               & ! intent(in):    model forcing data
                    bvar_data,               & ! intent(in):    average model variables for the entire basin
                    prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                    ! input-output: data stuctures
                    indx_data,               & ! intent(inout): index data
                    diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                    flux_init,               & ! intent(inout): model fluxes for a local HRU
                    deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    ! input-output: values needed in case canopy gets buried
                    scalarCanopyEnthalpy,    & ! intent(inout): value for enthalpy of the vegetation canopy (J m-3)
                    scalarCanopyTemp,        & ! intent(inout): value for temperature of the vegetation canopy (K), also used to start enthalpy calculations
                    scalarCanopyWat,         & ! intent(inout): value for total water content of the vegetation canopy (kg m-2)
                    ! output: new values of variables needed in data window outside of internal IDA for rootfinding and to start enthalpy calculations
                    mLayerTemp,              & ! intent(inout): vector of layer temperature (K)
                    mLayerMatricHead,        & ! intent(out):   value for total water matric potential (m)
                  ! output: new prime values of variables needed in data window outside of internal IDA for Jacobian
                    scalarCanopyTempPrime,   & ! intent(out):   prime value for temperature of the vegetation canopy (K s-1)
                    scalarCanopyWatPrime,    & ! intent(out):   prime value for total water content of the vegetation canopy (kg m-2 s-1)
                    mLayerTempPrime,         & ! intent(out):   prime vector of temperature of each snow and soil layer (K s-1)
                    mLayerMatricHeadPrime,   & ! intent(out):   prime vector of matric head of each snow and soil layer (m s-1)
                    mLayerVolFracWatPrime,   & ! intent(out):   prime vector of volumetric total water content of each snow and soil layer (s-1)
                    ! input-output: baseflow    
                    ixSaturation,            & ! intent(inout): index of the lowest saturated layer
                    dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                    ! output: flux and residual vectors
                    feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                    fluxVec0,                & ! intent(out):   flux vector
                    rAdd,                    & ! intent(out):   sink terms on the RHS of the flux equation
                    resVec,                  & ! intent(out):   residual vector
                    err,cmessage)              ! intent(out):   error control
  end associate
  if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! check for errors
#endif
 end subroutine initial_flux_and_residual_vectors_prime

 subroutine Newton_step
  ! ** Compute the Newton step using concepts from numerical recipes ** 
  associate(&
   ! layer geometry
   nSnow => indx_data%var(iLookINDEX%nSnow)%dat(1),& ! intent(in): [i4b] number of snow layers
   nSoil => indx_data%var(iLookINDEX%nSoil)%dat(1) & ! intent(in): [i4b] number of soil layers
   )
   call in_SS4HG % initialize(dt_cur,dt,iter,nSnow,nSoil,nLayers,nLeadDim,nState,ixMatrix,firstSubStep,computeVegFlux,scalarSolution,fOld)
   call io_SS4HG % initialize(firstFluxCall,xMin,xMax,ixSaturation)
   call summaSolve4homegrown(in_SS4HG,&                                                                                ! input: model control
                            &stateVecTrial,fScale,xScale,resVec,sMul,dMat,&                                            ! input: state vectors
                            &model_decisions,lookup_data,type_data,attr_data,mpar_data,forc_data,bvar_data,prog_data,& ! input: data structures
                            &indx_data,diag_data,flux_temp,deriv_data,&                                                ! input-output: data structures
                            &dBaseflow_dMatric,io_SS4HG,&                                                              ! input-output: baseflow
                            &stateVecNew,fluxVec,resSink,resVecNew,out_SS4HG)                                          ! output
   call io_SS4HG % finalize(firstFluxCall,xMin,xMax,ixSaturation)
   call out_SS4HG % finalize(fNew,converged,err,cmessage)                
  end associate
  if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! check for errors
 
  ! save the computed functions, residuals, and solution
  fOld          = fNew
  resVec        = resVecNew
  stateVecTrial = stateVecNew
  stateVecPrime = stateVecTrial  !prime values not used here, dummy
  nSteps = 1 ! number of time steps taken in solver
 end subroutine Newton_step

 subroutine check_Newton_convergence
  ! ** Check for convergence of current Newton step **    

  ! exit iteration loop if converged
  if (converged) then; exit_flag=.true.; return; end if
      
  ! check convergence
  if (iter==localMaxiter) then
    message=trim(message)//'failed to converge'
    err=-20; return_flag=.true.; return
  end if
 end subroutine check_Newton_convergence    

 subroutine enforce_mass_conservation
  ! Post processing step to “perfectly” conserve mass by pushing the errors into the state variables
  ! NOTE: if the residual is large this will cause the state variables to be pushed outside of their bounds
  layerVars: associate(&
    nSnow        => indx_data%var(iLookINDEX%nSnow)%dat(1)         ,& ! intent(in): [i4b] number of snow layers
    ! vector of energy and hydrology indices for the snow and soil domains
    ixSnowSoilNrg => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat   ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowSoilHyd => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat   ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    nSnowSoilNrg  => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1),& ! intent(in): [i4b] number of energy state variables in the snow+soil domain
    nSnowSoilHyd  => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1) & ! intent(in): [i4b] number of hydrology state variables in the snow+soil domain
    )
  
    ! update temperatures (ensure new temperature is consistent with the fluxes)
    if (nSnowSoilNrg>0) then
      do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing) ! loop through non-missing energy state variables in the snow+soil domain
        iState = ixSnowSoilNrg(iLayer)
        stateVecTrial(iState) = stateVecInit(iState) + (fluxVec(iState)*dt_cur + resSink(iState))/real(sMul(iState), rkind)
        resVec(iState) = 0._qp
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    end if
    
    ! update volumetric water content in the snow (ensure change in state is consistent with the fluxes)
    ! NOTE: for soil water balance is constrained within the iteration loop
    if (nSnowSoilHyd>0) then
      do concurrent (iLayer=1:nSnow,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing water state variables in the snow domain)
        iState = ixSnowSoilHyd(iLayer)
        stateVecTrial(iState) = stateVecInit(iState) + (fluxVec(iState)*dt_cur + resSink(iState))
        resVec(iState) = 0._qp
      end do  ! looping through non-missing water state variables in the soil domain
    end if
  end associate layerVars
 end subroutine enforce_mass_conservation

 subroutine solve_with_IDA
#ifdef SUNDIALS_ACTIVE
  ! get tolerance vectors
  call popTol4ida(&
                  ! input
                  nState,        & ! intent(in):  number of desired state variables
                  prog_data,     & ! intent(in):  model prognostic variables for a local HRU
                  diag_data,     & ! intent(in):  model diagnostic variables for a local HRU
                  indx_data,     & ! intent(in):  indices defining model states and layers
                  mpar_data,     & ! intent(in):  model parameters
                  ! output
                  atol,          & ! intent(out): absolute tolerances vector (mixed units)
                  rtol,          & ! intent(out): relative tolerances vector (mixed units)
                  err,cmessage)    ! intent(out): error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! check for errors

  layerGeometry: associate(&
   nSnow => indx_data%var(iLookINDEX%nSnow)%dat(1),& ! intent(in): [i4b] number of snow layers
   nSoil => indx_data%var(iLookINDEX%nSoil)%dat(1) & ! intent(in): [i4b] number of soil layers
   )

   ! allocate space for the temporary flux_sum structure
   call allocLocal(flux_meta(:),flux_sum,nSnow,nSoil,err,cmessage)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return; end if

   ! initialize flux_sum
   do concurrent ( iVar=1:size(flux_meta) )
    flux_sum%var(iVar)%dat(:) = 0._rkind
   end do
   ! initialize sum of compression of the soil matrix
   mLayerCmpress_sum(:) = 0._rkind
   stateVecNew(:) = 0._rkind 
   stateVecPrime(:) = 0._rkind

   !---------------------------
   ! * solving F(y,y') = 0 by IDA, y is the state vector and y' is the time derivative vector dy/dt
   !---------------------------
   ! iterations and updates to trial state vector, fluxes, and derivatives are done inside IDA solver
   call summaSolve4ida(&
                       dt_cur,                  & ! intent(in):    current stepsize
                       dt,                      & ! intent(in):    entire time step for drainage pond rate
                       atol,                    & ! intent(in):    absolute tolerance
                       rtol,                    & ! intent(in):    relative tolerance
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    number of snow+soil layers
                       nState,                  & ! intent(in):    number of state variables in the current subset
                       ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                       firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                       computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                       scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                       computMassBalance,       & ! intent(in):    flag to compute mass balance
                       computNrgBalance,        & ! intent(in):    flag to compute energy balance
                       ! input: state vector
                       stateVecTrial,           & ! intent(in):    model state vector at the beginning of the data time step
                       sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                       dMat,                    & ! intent(inout): diagonal of the Jacobian matrix (excludes fluxes)
                       ! input: data structures
                       model_decisions,         & ! intent(in):    model decisions
                       lookup_data,             & ! intent(in):    lookup data
                       type_data,               & ! intent(in):    type of vegetation and soil
                       attr_data,               & ! intent(in):    spatial attributes
                       mpar_data,               & ! intent(in):    model parameters
                       forc_data,               & ! intent(in):    model forcing data
                       bvar_data,               & ! intent(in):    average model variables for the entire basin
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       ! input-output: data structures
                       indx_data,               & ! intent(inout): index data
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       flux_temp,               & ! intent(inout): model fluxes for a local HRU
                       flux_sum,                & ! intent(inout): sum of fluxes model fluxes for a local HRU over a data step
                       deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       mLayerCmpress_sum,       & ! intent(inout): sum of compression of the soil matrix
                       ! output
                       ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       sunSucceeds,             & ! intent(out):   flag to indicate if ida successfully solved the problem in current data step
                       tooMuchMelt,             & ! intent(inout): flag to denote that there was too much melt
                       nSteps,                  & ! intent(out):   number of time steps taken in solver
                       stateVecNew,             & ! intent(inout): model state vector (y) at the end of the data time step
                       stateVecPrime,           & ! intent(inout): derivative of model state vector (y') at the end of the data time step
                       balance,                 & ! intent(inout): balance per state
                       err,cmessage)              ! intent(out):   error control
   ! check if IDA is successful, only fail outright in the case of a non-recoverable error
   if ( .not.sunSucceeds ) then
    message=trim(message)//trim(cmessage)
    !if (err.ne.-20 .or. err=0) err = 20 ! 0 if infeasible solution, could happen since not using imposeConstraints 
    if (err.ne.-20) err = 20 ! -20 is a recoverable error
    return
   else
    if (tooMuchMelt) return !exit to start same step over after merge
   end if
   niter = 0  ! iterations are counted inside IDA solver

   ! save the computed solution
   stateVecTrial = stateVecNew

   ! compute average flux
   do iVar=1,size(flux_meta)
    flux_temp%var(iVar)%dat(:) = ( flux_sum%var(iVar)%dat(:) ) /  dt_cur
   end do

   ! compute the total change in storage associated with compression of the soil matrix (kg m-2)
   soilVars: associate(&
    ! layer geometry
    mLayerDepth        => prog_data%var(iLookPROG%mLayerDepth)%dat          ,& ! depth of each layer in the snow-soil sub-domain (m)
    mLayerCompress     => diag_data%var(iLookDIAG%mLayerCompress)%dat       ,& ! change in storage associated with compression of the soil matrix (-)
    scalarSoilCompress => diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) & ! total change in storage associated with compression of the soil matrix (kg m-2 s-1)
    )
    mLayerCompress = mLayerCmpress_sum /  dt_cur
    scalarSoilCompress = sum(mLayerCompress(1:nSoil)*mLayerDepth(nSnow+1:nLayers))*iden_water
   end associate soilVars
  end associate layerGeometry
#endif
 end subroutine solve_with_IDA

 subroutine solve_with_KINSOL
#ifdef SUNDIALS_ACTIVE
  associate(&
   nSnow => indx_data%var(iLookINDEX%nSnow)%dat(1),& ! intent(in): [i4b] number of snow layers
   nSoil => indx_data%var(iLookINDEX%nSoil)%dat(1) & ! intent(in): [i4b] number of soil layers
   )
   !---------------------------
   ! * solving F(y) = 0 from Backward Euler with KINSOL, y is the state vector 
   !---------------------------
   stateVecNew(:) = 0._rkind
   ! iterations and updates to trial state vector, fluxes, and derivatives are done inside IDA solver
   call summaSolve4kinsol(&
                          dt_cur,                  & ! intent(in):    data time step
                          dt,                      & ! intent(in):    length of the entire time step (seconds) for drainage pond rate
                          fScale,                  & ! intent(in):    characteristic scale of the function evaluations
                          xScale,                  & ! intent(in):    characteristic scale of the state vector
                          nSnow,                   & ! intent(in):    number of snow layers
                          nSoil,                   & ! intent(in):    number of soil layers
                          nLayers,                 & ! intent(in):    number of snow+soil layers
                          nState,                  & ! intent(in):    number of state variables in the current subset
                          ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                          firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                          computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                          scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                          ! input: state vector
                          stateVecTrial,           & ! intent(in):    model state vector at the beginning of the data time step
                          sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                          dMat,                    & ! intent(inout)  diagonal of the Jacobian matrix (excludes fluxes)
                          ! input: data structures
                          model_decisions,         & ! intent(in):    model decisions
                          lookup_data,             & ! intent(in):    lookup tables
                          type_data,               & ! intent(in):    type of vegetation and soil
                          attr_data,               & ! intent(in):    spatial attributes
                          mpar_data,               & ! intent(in):    model parameters
                          forc_data,               & ! intent(in):    model forcing data
                          bvar_data,               & ! intent(in):    average model variables for the entire basin
                          prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                          ! input-output: data structures
                          indx_data,               & ! intent(inout): index data
                          diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                          flux_temp,               & ! intent(inout): model fluxes for a local HRU
                          deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                          ! output
                          ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                          sunSucceeds,             & ! intent(out):   flag to indicate if ida successfully solved the problem in current data step
                          stateVecNew,             & ! intent(inout): model state vector (y) at the end of the data time step
                          fluxVec,                 & ! intent(out):   new flux vector
                          resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                          resVec,                  & ! intent(out):   new residual vector       
                          err,cmessage)              ! intent(out):   error control
  end associate

  ! check if KINSOL is successful, only fail outright in the case of a non-recoverable error
  if ( .not.sunSucceeds ) then
    message=trim(message)//trim(cmessage)
    !if(err.ne.-20 .or. err=0) err = 20 ! 0 if infeasible solution, should not happen with imposeConstraints 
    if (err.ne.-20) err = 20 ! -20 if hit maximum iterations
    return
  end if
  niter = 0  ! iterations are counted inside KINSOL solver
  nSteps = 1 ! number of time steps taken in solver

  ! save the computed solution
  stateVecTrial = stateVecNew
  stateVecPrime = stateVecTrial ! prime values not used here, dummy
#endif
 end subroutine solve_with_KINSOL

 subroutine Newton_iterations_homegrown
  ! ** Compute the backward Euler solution using Newton iterations from homegrown solver **

  ! define maximum number of iterations
  maxiter = nint(mpar_data%var(iLookPARAM%maxiter)%dat(1))

  ! correct the number of iterations
  localMaxIter = merge(scalarMaxIter, maxIter, scalarSolution)

  !---------------------------
  ! * solving F(y) = 0 from Backward Euler using concepts from numerical recipes, y is the state vector 
  !---------------------------
  ! iterate and update trial state vector, fluxes, and derivatives
  exit_flag=.false. ! initialize exit flag
  do iter=1,localMaxIter ! begin Newton iterations
    niter = iter+1                ! # of iterations -- +1 because xFluxResid was moved outside the iteration loop (for backwards compatibility)
    call Newton_step; if (return_flag) return ! compute Newton step -- return if error                
    call check_Newton_convergence ! check current Newton step for convergence
    if (exit_flag) exit           ! exit loop if convereged
    if (return_flag) return       ! return if error
  end do 

  if (post_massCons) call enforce_mass_conservation ! enforce mass conservation if desired
 end subroutine Newton_iterations_homegrown

 subroutine finalize_systemSolv
  ! set untapped melt energy to zero
  untappedMelt(:) = 0._rkind

  ! free memory
  deallocate(mLayerCmpress_sum)
  deallocate(mLayerMatricHeadPrime)
  deallocate(dBaseflow_dMatric)
 end subroutine finalize_systemSolv

end subroutine systemSolv

end module systemSolv_module
