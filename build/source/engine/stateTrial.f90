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

module stateTrial_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix

implicit none
private
public::stateTrial
contains

 ! *********************************************************************************************************
 ! public subroutine stateTrial: calculate the iteration increment, evaluate the new state, and refine if necessary
 ! *********************************************************************************************************
 subroutine stateTrial(&
                       ! input: model control
                       dt,                      & ! intent(in):    length of the time step (seconds)
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       nState,                  & ! intent(in):    total number of state variables
                       ixMatrix,                & ! intent(in):    type of matrix (full or band diagonal)
                       firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                       firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                       computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                       ! input: state vectors
                       stateVecTrial,           & ! intent(in):    trial state vector
                       rVec,                    & ! intent(in):    residual vector
                       aJac,                    & ! intent(in):    Jacobian matrix
                       fScale,                  & ! intent(in):    function scaling vector
                       xScale,                  & ! intent(in):    "variable" scaling vector, i.e., for state variables
                       sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                       dMat,                    & ! intent(in):    diagonal matrix (excludes flux derivatives)
                       fOld,                    & ! intent(in):    old function evaluation
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
                       flux_data,               & ! intent(inout): model fluxes for a local HRU
                       deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ! output
                       stateVecNew,             & ! intent(out):   new state vector
                       fluxVecNew,              & ! intent(out):   new flux vector
                       resVecNew,               & ! intent(out):   new residual vector
                       fNew,                    & ! intent(out):   new function evaluation
                       converged,               & ! intent(out):   convergence flag
                       err,message)               ! intent(out):   error control
 ! provide access to the matrix routines
 USE matrixOper_module, only: scaleMatrices
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                       ! length of the time step (seconds)
 integer(i4b),intent(in)         :: nSnow                    ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                    ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                  ! total number of layers
 integer(i4b),intent(in)         :: nState                   ! total number of state variables
 integer(i4b),intent(in)         :: ixMatrix                 ! type of matrix (full or band diagonal) 
 logical(lgt),intent(in)         :: firstSubStep             ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall            ! flag to indicate if we are processing the first flux call
 logical(lgt),intent(in)         :: computeVegFlux           ! flag to indicate if computing fluxes over vegetation
 ! input: state vectors
 real(dp),intent(in)             :: stateVecTrial(:)         ! trial state vector
 real(qp),intent(in)             :: rVec(:)   ! NOTE: qp     ! residual vector
 real(dp),intent(in)             :: aJac(:,:)                ! Jacobian matrix
 real(dp),intent(in)             :: fScale(:)                ! function scaling vector
 real(dp),intent(in)             :: xScale(:)                ! "variable" scaling vector, i.e., for state variables
 real(qp),intent(in)             :: sMul(:)   ! NOTE: qp     ! state vector multiplier (used in the residual calculations)
 real(dp),intent(in)             :: dMat(:)                  ! diagonal matrix (excludes flux derivatives) 
 real(dp),intent(in)             :: fOld                     ! old function evaluation
 ! input: data structures
 type(model_options),intent(in)  :: model_decisions(:)       ! model decisions
 type(var_i),        intent(in)  :: type_data                ! type of vegetation and soil
 type(var_d),        intent(in)  :: attr_data                ! spatial attributes
 type(var_d),        intent(in)  :: mpar_data                ! model parameters
 type(var_d),        intent(in)  :: forc_data                ! model forcing data
 type(var_dlength),  intent(in)  :: bvar_data                ! model variables for the local basin
 type(var_dlength),  intent(in)  :: prog_data                ! prognostic variables for a local HRU
 type(var_ilength),  intent(in)  :: indx_data                ! indices defining model states and layers
 ! output: data structures
 type(var_dlength),intent(inout) :: diag_data                ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: deriv_data               ! derivatives in model fluxes w.r.t. relevant state variables
 ! output: flux and residual vectors
 real(dp),intent(out)            :: stateVecNew(:)           ! new state vector
 real(dp),intent(out)            :: fluxVecNew(:)            ! new flux vector
 real(qp),intent(out)            :: resVecNew(:) ! NOTE: qp  ! new residual vector
 real(dp),intent(out)            :: fNew                     ! new function evaluation
 logical(lgt),intent(out)        :: converged                ! convergence flag
 ! output: error control
 integer(i4b),intent(out)        :: err                      ! error code
 character(*),intent(out)        :: message                  ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 real(dp),allocatable            :: aJacScaled(:,:)          ! scaled Jacobian matrix
 real(dp)                        :: rVecscaled(nState)       ! scaled residual vector
 integer(i4b),parameter          :: ixLineSearch=1001        ! step refinement = line search
 integer(i4b),parameter          :: ixTrustRegion=1002       ! step refinement = trust region
 integer(i4b),parameter          :: ixStepRefinement=ixLineSearch   ! decision for the numerical solution
 character(LEN=256)              :: cmessage                 ! error message of downwind routine
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='stateTrial/'

 ! allocate space
 select case(ixMatrix)
  case(ixFullMatrix); allocate(aJacScaled(nState,   nState), stat=err)
  case(ixBandMatrix); allocate(aJacScaled(nBands+kl,nState), stat=err)
  case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'
 end select
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the Jacobian matrix'; return; endif

 ! scale matrices
 call scaleMatrices(ixMatrix,nState,aJac,fScale,xScale,aJacScaled,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! scale the residual vector
 rVecScaled(1:nState) = fScale(:)*real(rVec(:), dp)   ! NOTE: residual vector is in quadruple precision

 ! compute the flux vector and the residual, and refine the iteration increment
 ! NOTE: aJacScaled is intent(inout) -- it is decomposed in lapackSolv
 select case(ixStepRefinement)
  case(ixLineSearch);  call lineSearchRefinement( stateVecTrial,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,cmessage)
  case(ixTrustRegion); call trustRegionRefinement(stateVecTrial,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,cmessage)
  case default; err=20; message=trim(message)//'unable to identify numerical solution'; return
 end select
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! deallocate space
 deallocate(aJacScaled, stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to de allocate space'; return; endif

 contains

  ! *********************************************************************************************************
  ! * internal subroutine lineSearchRefinement: refine the iteration increment using line searches
  ! *********************************************************************************************************
  subroutine lineSearchRefinement(stateVecTrial,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,message)
  ! provide access to the matrix routines
  USE matrixOper_module, only: lapackSolv
  USE matrixOper_module, only: computeGradient
  implicit none
  ! input
  real(dp),intent(in)            :: stateVecTrial(:)         ! trial state vector
  real(dp),intent(inout)         :: aJacScaled(:,:)          ! scaled jacobian matrix
  real(dp),intent(in)            :: rVecScaled(:)            ! scaled residual vector
  real(dp),intent(in)            :: fOld                     ! old function value
  ! output
  real(dp),intent(out)           :: stateVecNew(:)           ! new state vector
  real(dp),intent(out)           :: fluxVecNew(:)            ! new flux vector
  real(qp),intent(out)           :: resVecNew(:) ! NOTE: qp  ! new residual vector
  real(dp),intent(out)           :: fNew                     ! new function evaluation
  logical(lgt),intent(out)       :: converged                ! convergence flag
  integer(i4b),intent(out)       :: err                      ! error code
  character(*),intent(out)       :: message                  ! error message
  ! --------------------------------------------------------------------------------------------------------
  ! local
  character(len=256)             :: cmessage                 ! error message of downwind routine
  real(dp)                       :: newtStep(nState)         ! scaled newton step
  real(dp)                       :: xIncScaled(nState)       ! scaled newton step (restricted)
  real(dp)                       :: gradScaled(nState)       ! scaled gradient
  logical(lgt)                   :: feasible                 ! flag to denote the feasibility of the solution
  integer(i4b)                   :: iLine                    ! line search index
  integer(i4b),parameter         :: maxLineSearch=20         ! maximum number of backtracks
  real(dp),parameter             :: alpha=1.e-4_dp           ! check on gradient
  real(dp)                       :: xLambda                  ! backtrack magnitude
  real(dp)                       :: xLambdaTemp              ! temporary backtrack magnitude
  real(dp)                       :: slopeInit                ! initial slope
  real(dp)                       :: rVecTemp(nState)         ! temporary residual vector
  real(dp)                       :: rhs1,rhs2                ! rhs used to compute the cubic
  real(dp)                       :: aCoef,bCoef              ! coefficients in the cubic
  real(dp)                       :: disc                     ! temporary variable used in cubic
  real(dp)                       :: xLambdaPrev              ! previous lambda value (used in the cubic)
  real(dp)                       :: fPrev                    ! previous function evaluation (used in the cubic)
  ! --------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='lineSearchRefinement/'

  ! compute the newton step: use the lapack routines to solve the linear system A.X=B
  call lapackSolv(ixMatrix,nState,aJacScaled,-rVecScaled,newtStep,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! compute the gradient of the function vector
  call computeGradient(ixMatrix,nState,aJacScaled,rVecScaled,gradScaled,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! compute the initial slope
  slopeInit = dot_product(gradScaled,newtStep)

  ! initialize lambda
  xLambda=1._dp

  ! ***** LINE SEARCH LOOP...
  lineSearch: do iLine=1,maxLineSearch  ! try to refine the function by shrinking the step size

   ! compute the iteration increment
   xIncScaled = xLambda*newtStep

   ! compute the residual vector and function
   call eval8state(stateVecTrial,xIncScaled,stateVecNew,fluxVecNew,resVecNew,fNew,feasible,converged,err,message)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

   ! check feasibility
   if(.not.feasible) cycle

   ! check convergence
   if(converged) return

   ! check if the function is accepted
   if(fNew < fold + alpha*slopeInit*xLambda) return

   ! ***
   ! *** IF GET TO HERE WE BACKTRACK
   !      --> all remaining code simply computes the restricted step multiplier (xLambda)

   ! first backtrack: use quadratic
   if(iLine==1)then
    xLambdaTemp = -slopeInit / 2._dp*(fNew - fOld - slopeInit)

   ! subsequent backtracks: use cubic
   else

    ! check that we did not back-track all the way back to the original value
    if(iLine==maxLineSearch)then
     message=trim(message)//'backtracked all the way back to the original value'
     err=-20; return
    endif

    ! define rhs
    rhs1 = fNew - fOld - xLambda*slopeInit
    rhs2 = fPrev - fOld - xLambdaPrev*slopeInit

    ! define coefficients
    aCoef = (rhs1/(xLambda*xLambda) - rhs2*(xLambdaPrev*xLambdaPrev))/(xLambda - xLambdaPrev)
    bCoef = (-xLambdaPrev*rhs1/(xLambda*xLambda) + xLambda*rhs2/(xLambdaPrev*xLambdaPrev)) / (xLambda - xLambdaPrev)

    ! check if a quadratic
    if(aCoef==0._dp)then
     xLambdaTemp = -slopeInit/(2._dp*bCoef)

    ! calculate cubic
    else
     disc = bCoef*bCoef - 3._dp*aCoef*slopeInit
     if(disc < 0._dp)then
      xLambdaTemp = 0.5_dp*xLambda
     elseif(bCoef <= 0._dp)then
      xLambdaTemp = (-bCoef + sqrt(disc))/(3._dp*aCoef)
     else
      xLambdaTemp = -slopeInit/(bCoef + sqrt(disc))
     endif
    endif  ! calculating cubic

    ! constrain to <= 0.5*xLambda
    if(xLambdaTemp < 0.5_dp*xLambda) xLambdaTemp=0.5_dp*xLambda

   endif  ! subsequent backtracks

   ! save results
   xLambdaPrev = xLambda
   fPrev = fNew

   ! constrain lambda
   xLambda = max(xLambdaTemp, 0.1_dp*xLambda)

  end do lineSearch  ! backtrack loop

  end subroutine lineSearchRefinement


  ! *********************************************************************************************************
  ! * internal subroutine trustRegionRefinement: refine the iteration increment using trust regions
  ! *********************************************************************************************************
  subroutine trustRegionRefinement(stateVecTrial,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,message)
  ! provide access to the matrix routines
  USE matrixOper_module, only: lapackSolv
  USE matrixOper_module, only: computeGradient
  implicit none
  ! input
  real(dp),intent(in)            :: stateVecTrial(:)         ! trial state vector
  real(dp),intent(inout)         :: aJacScaled(:,:)          ! scaled jacobian matrix
  real(dp),intent(in)            :: rVecScaled(:)            ! scaled residual vector
  real(dp),intent(in)            :: fOld                     ! old function value
  ! output
  real(dp),intent(out)           :: stateVecNew(:)           ! new state vector
  real(dp),intent(out)           :: fluxVecNew(:)            ! new flux vector
  real(qp),intent(out)           :: resVecNew(:) ! NOTE: qp  ! new residual vector
  real(dp),intent(out)           :: fNew                     ! new function evaluation
  logical(lgt),intent(out)       :: converged                ! convergence flag
  integer(i4b),intent(out)       :: err                      ! error code
  character(*),intent(out)       :: message                  ! error message
  ! --------------------------------------------------------------------------------------------------------
  ! local variables

  ! .. needed ..


  ! --------------------------------------------------------------------------------------------------------
  err=0; message='trustRegionRefinement/'

  ! dummy
  message=trim(message)//'routine not implemented yet'
  err=20; return

  ! NOTE: may need a copy of aJacScaled: It is decomposed in lapackSolv


  end subroutine trustRegionRefinement


  ! *********************************************************************************************************
  ! * internal subroutine eval8state: compute the right-hand-side vector
  ! *********************************************************************************************************
  subroutine eval8state(stateVecTrial,stateVecIter,stateVecNew,fluxVecNew,resVecNew,fNew,feasible,converged,err,message)
  USE eval8summa_module,only:eval8summa                      ! simulation of fluxes and residuals given a trial state vector
  implicit none
  ! input
  real(dp),intent(in)            :: stateVecTrial(:)         ! trial state vector 
  real(dp),intent(in)            :: stateVecIter(:)          ! *scaled* iteration increment in model state vector
  ! output
  real(dp),intent(out)           :: stateVecNew(:)           ! updated state vector
  real(dp),intent(out)           :: fluxVecNew(:)            ! updated flux vector
  real(qp),intent(out)           :: resVecNew(:) ! NOTE: qp  ! updated residual vector
  real(dp),intent(out)           :: fNew                     ! new function value
  logical(lgt),intent(out)       :: feasible                 ! flag to denote the feasibility of the solution
  logical(lgt),intent(out)       :: converged                ! flag to denote the solution has converged
  integer(i4b),intent(out)       :: err                      ! error code
  character(*),intent(out)       :: message                  ! error message
  ! ----------------------------------------------------------------------------------------------------------
  ! local
  character(len=256)             :: cmessage                 ! error message of downwind routine
  real(dp)                       :: xInc(nState)             ! iteration increment (re-scaled to original units of the state vector)
  ! ----------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='eval8state/'

  ! re-scale the iteration increment
  xInc(:) = stateVecIter(:)*xScale(:)

  ! if enthalpy, then need to convert the iteration increment to temperature
  !if(nrgFormulation==ix_enthalpy) xInc(ixNrgOnly)/dMat(ixNrgOnly)

  ! impose solution constraints
  ! NOTE: we may not need to do this (or at least, do ALL of this), as we can probably rely on trust regions here
  !call imposeConstraints()

  ! update the state vector
  stateVecNew = stateVecTrial + xInc

  ! compute the flux and the residual vector for a given state vector
  call eval8summa(&
                  ! input: model control
                  dt,                      & ! intent(in):    length of the time step (seconds)
                  nSnow,                   & ! intent(in):    number of snow layers
                  nSoil,                   & ! intent(in):    number of soil layers
                  nLayers,                 & ! intent(in):    total number of layers
                  nState,                  & ! intent(in):    total number of state variables
                  firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                  firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                  computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  ! input: state vectors
                  stateVecNew,             & ! intent(in):    updated model state vector
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
                  flux_data,               & ! intent(inout): model fluxes for a local HRU
                  deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  ! output
                  feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                  fluxVecNew,              & ! intent(out):   new flux vector
                  resVecNew,               & ! intent(out):   new residual vector
                  fNew,                    & ! intent(out):   new function evaluation
                  err,cmessage)              ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! check convergence
  converged = checkConv(resVecNew,xInc,stateVecNew)

  end subroutine eval8state


  ! *********************************************************************************************************
  ! internal function checkConv: check convergence based on the residual vector
  ! *********************************************************************************************************
  function checkConv(rVec,xInc,xVec)
  ! provide access to named variables that define elements of the structure
  USE var_lookup,only:iLookPROG                       ! named variables for structure elements
  USE var_lookup,only:iLookPARAM                      ! named variables for structure elements
  USE var_lookup,only:iLookINDEX                      ! named variables for structure elements
  ! constants
  USE multiconst,only:iden_water                      ! intrinsic density of liquid water    (kg m-3)
  implicit none
  ! dummies
  real(qp),intent(in)       :: rVec(:)                ! residual vector (mixed units)
  real(dp),intent(in)       :: xInc(:)                ! iteration increment (mixed units)
  real(dp),intent(in)       :: xVec(:)                ! state vector (mixed units)
  logical(lgt)              :: checkConv              ! flag to denote convergence
  ! locals
  real(dp),dimension(nSoil) :: psiScale               ! scaling factor for matric head
  real(dp),parameter        :: xSmall=1.e-0_dp        ! a small offset
  real(dp)                  :: soilWatbalErr          ! error in the soil water balance
  real(dp)                  :: canopy_max             ! absolute value of the residual in canopy water (kg m-2)
  real(dp),dimension(1)     :: energy_max             ! maximum absolute value of the energy residual (J m-3)
  real(dp),dimension(1)     :: liquid_max             ! maximum absolute value of the volumetric liquid water content residual (-)
  real(dp),dimension(1)     :: matric_max             ! maximum absolute value of the matric head iteration increment (m)
  integer(i4b),dimension(1) :: energy_loc             ! location of maximum absolute value of the energy residual (-)
  integer(i4b),dimension(1) :: liquid_loc             ! location of maximum absolute value of the volumetric liquid water content residual (-)
  integer(i4b),dimension(1) :: matric_loc             ! location of maximum absolute value of the matric head increment (-)
  logical(lgt)              :: canopyConv             ! flag for canopy water balance convergence
  logical(lgt)              :: watbalConv             ! flag for soil water balance convergence
  logical(lgt)              :: liquidConv             ! flag for residual convergence
  logical(lgt)              :: matricConv             ! flag for matric head convergence
  logical(lgt)              :: energyConv             ! flag for energy convergence
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  associate(&
  ! convergence parameters
  absConvTol_liquid       => mpar_data%var(iLookPARAM%absConvTol_liquid)            ,&  ! intent(in): [dp] absolute convergence tolerance for vol frac liq water (-)
  absConvTol_matric       => mpar_data%var(iLookPARAM%absConvTol_matric)            ,&  ! intent(in): [dp] absolute convergence tolerance for matric head        (m)
  absConvTol_energy       => mpar_data%var(iLookPARAM%absConvTol_energy)            ,&  ! intent(in): [dp] absolute convergence tolerance for energy             (J m-3)
  ! layer depth
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)
  ! model indices
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,&  ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
  ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow+soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat             &  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
  ) ! making associations with variables in the data structures
  ! -------------------------------------------------------------------------------------------------------------------------------------------------

  ! check convergence based on the residuals for energy (J m-3)
  if(computeVegFlux)then
   canopy_max = real(abs(rVec(ixVegWat)), dp)*iden_water
   energy_max = real(maxval(abs( (/rVec(ixCasNrg), rVec(ixVegNrg), rVec(ixSnowSoilNrg)/) ) ), dp)
   energy_loc =      maxloc(abs( (/rVec(ixCasNrg), rVec(ixVegNrg), rVec(ixSnowSoilNrg)/) ) )
   canopyConv = (canopy_max    < absConvTol_liquid)  ! absolute error in canopy water balance (m)
  else
   energy_max = real(maxval(abs( rVec(ixSnowSoilNrg) ) ), dp)
   energy_loc =      maxloc(abs( rVec(ixSnowSoilNrg) ) )
   canopyConv = .true. ! disable check for canopy convergence when there is no canopy
  endif

  ! check convergence based on the residuals for volumetric liquid water content (-)
  liquid_max = real(maxval(abs( rVec(ixSnowSoilWat) ) ), dp)
  liquid_loc =      maxloc(abs( rVec(ixSnowSoilWat) ) )

  ! check convergence based on the iteration increment for matric head
  ! NOTE: scale by matric head to avoid unnecessairly tight convergence when there is no water
  psiScale   = abs(xVec(ixSoilOnlyHyd)) + xSmall ! avoid divide by zero
  matric_max = maxval(abs( xInc(ixSoilOnlyHyd)/psiScale ) )
  matric_loc = maxloc(abs( xInc(ixSoilOnlyHyd)/psiScale ) )

  ! compute the soil water balance error (m)
  soilWatbalErr = abs( sum(real(rVec(ixSoilOnlyHyd), dp)*mLayerDepth(nSnow+1:nLayers)) )

  ! convergence check
  watbalConv = (soilWatbalErr < absConvTol_liquid)  ! absolute error in total soil water balance (m)
  matricConv = (matric_max(1) < absConvTol_matric)  ! NOTE: based on iteration increment
  liquidConv = (liquid_max(1) < absConvTol_liquid)  ! (based on the residual)
  energyConv = (energy_max(1) < absConvTol_energy)  ! (based on the residual)

  ! final convergence check
  checkConv = (canopyConv .and. watbalConv .and. matricConv .and. liquidConv .and. energyConv)

  ! end associations with variables in the data structures
  end associate

  end function checkConv

 end subroutine stateTrial

end module stateTrial_module
