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

module summaSolve_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE multiconst,only:integerMissing  ! missing integer
USE multiconst,only:realMissing     ! missing double precision number
USE multiconst,only:quadMissing     ! missing quadruple precision number

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
                    qbaseTopmodel,& ! TOPMODEL-ish baseflow parameterization
                    bigBucket,    & ! a big bucket (lumped aquifer model)
                    noExplicit      ! no explicit groundwater parameterization

implicit none
private
public::summaSolve
contains

 ! *********************************************************************************************************
 ! public subroutine summaSolve: calculate the iteration increment, evaluate the new state, and refine if necessary
 ! *********************************************************************************************************
 subroutine summaSolve(&
                       ! input: model control
                       dt,                      & ! intent(in):    length of the time step (seconds)
                       iter,                    & ! intent(in):    iteration index
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       nLeadDim,                & ! intent(in):    length of the leading dimension of he Jacobian matrix (either nBands or nState)
                       nState,                  & ! intent(in):    total number of state variables
                       ixMatrix,                & ! intent(in):    type of matrix (full or band diagonal)
                       firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                       firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                       computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                       canopyDepth,             & ! intent(in):    depth of the vegetation canopy (m)
                       ! input: state vectors
                       stateVecTrial,           & ! intent(in):    trial state vector
                       fScale,                  & ! intent(in):    function scaling vector
                       xScale,                  & ! intent(in):    "variable" scaling vector, i.e., for state variables
                       rVec,                    & ! intent(in):    residual vector
                       sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                       dMat,                    & ! intent(inout): diagonal matrix (excludes flux derivatives)
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
                       ! input-output: baseflow
                       ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       dBaseflow_dMatric,       & ! intent(inout): derivative in baseflow w.r.t. matric head (s-1)
                       ! output
                       stateVecNew,             & ! intent(out):   new state vector
                       fluxVecNew,              & ! intent(out):   new flux vector
                       resSinkNew,              & ! intent(out):   additional (sink) terms on the RHS of the state equation
                       resVecNew,               & ! intent(out):   new residual vector
                       fNew,                    & ! intent(out):   new function evaluation
                       converged,               & ! intent(out):   convergence flag
                       err,message)               ! intent(out):   error control
 USE computJacob_module, only: computJacob
 USE matrixOper_module,  only: lapackSolv
 USE matrixOper_module,  only: scaleMatrices
 USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                       ! length of the time step (seconds)
 integer(i4b),intent(in)         :: iter                     ! interation index
 integer(i4b),intent(in)         :: nSnow                    ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                    ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                  ! total number of layers
 integer(i4b),intent(in)         :: nLeadDim                 ! length of the leading dimension of the Jacobian matrix (nBands or nState)
 integer(i4b),intent(in)         :: nState                   ! total number of state variables
 integer(i4b),intent(in)         :: ixMatrix                 ! type of matrix (full or band diagonal) 
 logical(lgt),intent(in)         :: firstSubStep             ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall            ! flag to indicate if we are processing the first flux call
 logical(lgt),intent(in)         :: computeVegFlux           ! flag to indicate if computing fluxes over vegetation
 real(dp),intent(in)             :: canopyDepth              ! depth of the vegetation canopy (m)
 ! input: state vectors
 real(dp),intent(in)             :: stateVecTrial(:)         ! trial state vector
 real(dp),intent(in)             :: fScale(:)                ! function scaling vector
 real(dp),intent(in)             :: xScale(:)                ! "variable" scaling vector, i.e., for state variables
 real(qp),intent(in)             :: rVec(:)   ! NOTE: qp     ! residual vector
 real(qp),intent(in)             :: sMul(:)   ! NOTE: qp     ! state vector multiplier (used in the residual calculations)
 real(dp),intent(inout)          :: dMat(:)                  ! diagonal matrix (excludes flux derivatives) 
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
 ! input-output: baseflow
 integer(i4b),intent(inout)      :: ixSaturation             ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),intent(inout)          :: dBaseflow_dMatric(:,:)   ! derivative in baseflow w.r.t. matric head (s-1)
 ! output: flux and residual vectors
 real(dp),intent(out)            :: stateVecNew(:)           ! new state vector
 real(dp),intent(out)            :: fluxVecNew(:)            ! new flux vector
 real(dp),intent(out)            :: resSinkNew(:)            ! sink terms on the RHS of the flux equation
 real(qp),intent(out)            :: resVecNew(:) ! NOTE: qp  ! new residual vector
 real(dp),intent(out)            :: fNew                     ! new function evaluation
 logical(lgt),intent(out)        :: converged                ! convergence flag
 ! output: error control
 integer(i4b),intent(out)        :: err                      ! error code
 character(*),intent(out)        :: message                  ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! Jacobian matrix
 real(dp)                        :: aJac(nLeadDim,nState)      ! Jacobian matrix
 real(dp)                        :: aJacScaled(nLeadDim,nState)      ! Jacobian matrix (scaled)
 real(dp)                        :: aJacScaledTemp(nLeadDim,nState)  ! Jacobian matrix (scaled) -- temporary copy since decomposed in lapack
 ! solution/step vectors
 real(dp),dimension(nState)      :: rVecScaled               ! residual vector (scaled)
 real(dp),dimension(nState)      :: newtStepScaled           ! full newton step (scaled)
 ! step size refinement
 logical(lgt)                    :: doRefine                 ! flag for step refinement
 integer(i4b),parameter          :: ixLineSearch=1001        ! step refinement = line search
 integer(i4b),parameter          :: ixTrustRegion=1002       ! step refinement = trust region
 integer(i4b),parameter          :: ixStepRefinement=ixLineSearch   ! decision for the numerical solution
 ! general
 integer(i4b)                    :: iLayer                   ! row index
 integer(i4b)                    :: jLayer                   ! column index
 character(LEN=256)              :: cmessage                 ! error message of downwind routine
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! associations to information in data structures
 associate(ixGroundwater => model_decisions(iLookDECISIONS%groundwatr)%iDecision)  ! intent(in): [i4b] groundwater parameterization 
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summaSolve/'

 ! -----
 ! * compute the Jacobian matrix...
 ! --------------------------------

 ! compute the analytical Jacobian matrix
 ! NOTE: The derivatives were computed in the previous call to computFlux
 !       This occurred either at the call to eval8summa at the start of systemSolv
 !        or in the call to eval8summa in the previous iteration (within lineSearchRefinement or trustRegionRefinement)
 call computJacob(&
                  ! input: model control
                  dt,                             & ! intent(in):    length of the time step (seconds)
                  nSnow,                          & ! intent(in):    number of snow layers
                  nSoil,                          & ! intent(in):    number of soil layers
                  nLayers,                        & ! intent(in):    total number of layers
                  canopyDepth,                    & ! intent(in):    depth of the vegetation canopy (m)
                  computeVegFlux,                 & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  (ixGroundwater==qbaseTopmodel), & ! intent(in):    flag to indicate if we need to compute baseflow
                  ixMatrix,                       & ! intent(in):    form of the Jacobian matrix
                  ! input: data structures        
                  indx_data,                      & ! intent(in):    index data
                  prog_data,                      & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                      & ! intent(in):    model diagnostic variables for a local HRU
                  deriv_data,                     & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                  dBaseflow_dMatric,              & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                  ! input-output: Jacobian and its diagonal
                  dMat,                           & ! intent(inout): diagonal of the Jacobian matrix
                  aJac,                           & ! intent(out):   Jacobian matrix
                  ! output: error control
                  err,cmessage)                     ! intent(out):   error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 ! -----
 ! * solve linear system...
 ! ------------------------

 ! scale the residual vector
 rVecScaled(1:nState) = fScale(:)*real(rVec(:), dp)   ! NOTE: residual vector is in quadruple precision

 ! scale matrices
 call scaleMatrices(ixMatrix,nState,aJac,fScale,xScale,aJacScaled,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 if(globalPrintFlag)then
  print*, '** SCALED banded analytical Jacobian:'
  write(*,'(a4,1x,100(i17,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
  do iLayer=kl+1,nBands
   write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJacScaled(iLayer,jLayer),jLayer=iJac1,iJac2)
  end do
 end if

 ! copy the scaled matrix, since it is decomposed in lapackSolv
 aJacScaledTemp = aJacScaled

 ! compute the newton step: use the lapack routines to solve the linear system A.X=B
 call lapackSolv(ixMatrix,nState,aJacScaledTemp,-rVecScaled,newtStepScaled,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 if(globalPrintFlag) write(*,'(a,1x,10(e17.10,1x))') 'newtStepScaled = ', newtStepScaled(iJac1:iJac2)

 ! -----
 ! * update, evaluate, and refine the state vector...
 ! --------------------------------------------------

 ! initialize the flag for step refinement
 doRefine=.true.

 ! compute the flux vector and the residual, and (if necessary) refine the iteration increment
 ! NOTE: in 99.9% of cases newtStep will be used (no refinement)
 select case(ixStepRefinement)
  case(ixLineSearch);  call lineSearchRefinement( doRefine,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,cmessage)
  case(ixTrustRegion); call trustRegionRefinement(doRefine,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,cmessage)
  case default; err=20; message=trim(message)//'unable to identify numerical solution'; return
 end select

 ! check warnings: negative error code = warning; in this case back-tracked to the original value
 ! NOTE: Accept the full newton step
 if(err<0)then
  doRefine=.false.;    call lineSearchRefinement( doRefine,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,cmessage)
 end if

 ! check errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 ! end association to info in data structures
 end associate

 contains

  ! *********************************************************************************************************
  ! * internal subroutine lineSearchRefinement: refine the iteration increment using line searches
  ! *********************************************************************************************************
  subroutine lineSearchRefinement(doLineSearch,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,message)
  ! provide access to the matrix routines
  USE matrixOper_module, only: computeGradient
  implicit none
  ! input
  logical(lgt),intent(in)        :: doLineSearch             ! flag to do the line search
  real(dp),intent(in)            :: stateVecTrial(:)         ! trial state vector
  real(dp),intent(in)            :: newtStepScaled(:)        ! scaled newton step
  real(dp),intent(in)            :: aJacScaled(:,:)          ! scaled jacobian matrix
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
  real(dp)                       :: gradScaled(nState)       ! scaled gradient
  real(dp)                       :: xInc(nState)             ! iteration increment (re-scaled to original units of the state vector)
  logical(lgt)                   :: feasible                 ! flag to denote the feasibility of the solution
  integer(i4b)                   :: iLine                    ! line search index
  integer(i4b),parameter         :: maxLineSearch=5          ! maximum number of backtracks
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

  ! check the need to compute the line search
  if(doLineSearch)then

   ! compute the gradient of the function vector
   call computeGradient(ixMatrix,nState,aJacScaled,rVecScaled,gradScaled,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
 
   ! compute the initial slope
   slopeInit = dot_product(gradScaled,newtStepScaled)

  end if  ! if computing the line search

  ! initialize lambda
  xLambda=1._dp

  ! ***** LINE SEARCH LOOP...
  lineSearch: do iLine=1,maxLineSearch  ! try to refine the function by shrinking the step size

   ! back-track along the search direction
   ! NOTE: start with back-tracking the scaled step
   xInc(:) = xLambda*newtStepScaled(:)

   ! re-scale the iteration increment
   xInc(:) = xInc(:)*xScale(:)

   ! if enthalpy, then need to convert the iteration increment to temperature
   !if(nrgFormulation==ix_enthalpy) xInc(ixNrgOnly) = xInc(ixNrgOnly)/dMat(ixNrgOnly)

   ! impose solution constraints
   ! NOTE: we may not need to do this (or at least, do ALL of this), as we can probably rely on the line search here
   !  (especially the feasibility check)
   call imposeConstraints(stateVecTrial,xInc,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

   ! compute the iteration increment
   stateVecNew = stateVecTrial + xInc

   ! compute the residual vector and function
   ! NOTE: This calls eval8summa in an internal subroutine
   !       The internal sub routine has access to all data
   !       Hence, we only need to include the variables of interest in lineSearch
   call eval8summa_wrapper(stateVecNew,fluxVecNew,resVecNew,fNew,feasible,err,message)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

   ! check line search
   if(globalPrintFlag)then
    write(*,'(a,1x,i4,1x,e17.10)' ) 'iLine, xLambda                 = ', iLine, xLambda
    write(*,'(a,1x,10(e17.10,1x))') 'fOld,fNew                      = ', fOld,fNew
    write(*,'(a,1x,10(e17.10,1x))') 'fold + alpha*slopeInit*xLambda = ', fold + alpha*slopeInit*xLambda
    write(*,'(a,1x,10(e17.10,1x))') 'resVecNew                      = ', resVecNew(iJac1:iJac2)
    write(*,'(a,1x,10(e17.10,1x))') 'xInc                           = ', xInc(iJac1:iJac2)
   end if

   ! check feasibility
   if(.not.feasible) cycle

   ! check convergence
   ! NOTE: some efficiency gains possible by scaling the full newton step outside the line search loop
   converged = checkConv(resVecNew,newtStepScaled*xScale,stateVecNew)
   if(converged) return

   ! early return if not computing the line search
   if(.not.doLineSearch) return

   ! check if the function is accepted
   if(fNew < fold + alpha*slopeInit*xLambda) return

   ! ***
   ! *** IF GET TO HERE WE BACKTRACK
   !      --> all remaining code simply computes the restricted step multiplier (xLambda)

   ! first backtrack: use quadratic
   if(iLine==1)then
    xLambdaTemp = -slopeInit / (2._dp*(fNew - fOld - slopeInit) )
    if(xLambdaTemp > 0.5_dp*xLambda) xLambdaTemp = 0.5_dp*xLambda

   ! subsequent backtracks: use cubic
   else

    ! check that we did not back-track all the way back to the original value
    if(iLine==maxLineSearch)then
     message=trim(message)//'backtracked all the way back to the original value'
     err=-20; return
    end if

    ! define rhs
    rhs1 = fNew - fOld - xLambda*slopeInit
    rhs2 = fPrev - fOld - xLambdaPrev*slopeInit

    ! define coefficients
    aCoef = (rhs1/(xLambda*xLambda) - rhs2/(xLambdaPrev*xLambdaPrev))/(xLambda - xLambdaPrev)
    bCoef = (-xLambdaPrev*rhs1/(xLambda*xLambda) + xLambda*rhs2/(xLambdaPrev*xLambdaPrev)) / (xLambda - xLambdaPrev)

    ! check if a quadratic
    if(aCoef==0._dp)then
     xLambdaTemp = -slopeInit/(2._dp*bCoef)

    ! calculate cubic
    else
     disc = bCoef*bCoef - 3._dp*aCoef*slopeInit
     if(disc < 0._dp)then
      xLambdaTemp = 0.5_dp*xLambda
     else
      xLambdaTemp = (-bCoef + sqrt(disc))/(3._dp*aCoef)
     end if
    end if  ! calculating cubic

    ! constrain to <= 0.5*xLambda
    if(xLambdaTemp > 0.5_dp*xLambda) xLambdaTemp=0.5_dp*xLambda

   end if  ! subsequent backtracks

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
  subroutine trustRegionRefinement(doTrustRefinement,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,fOld,stateVecNew,fluxVecNew,resVecNew,fNew,converged,err,message)
  ! provide access to the matrix routines
  USE matrixOper_module, only: lapackSolv
  USE matrixOper_module, only: computeGradient
  implicit none
  ! input
  logical(lgt),intent(in)        :: doTrustRefinement        ! flag to refine using trust regions
  real(dp),intent(in)            :: stateVecTrial(:)         ! trial state vector
  real(dp),intent(in)            :: newtStepScaled(:)        ! scaled newton step
  real(dp),intent(in)            :: aJacScaled(:,:)          ! scaled jacobian matrix
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
  stateVecNew = realMissing
  fluxVecNew  = realMissing
  resVecNew   = quadMissing
  fNew        = realMissing
  converged   = .true.


  message=trim(message)//'routine not implemented yet'
  err=20; return



  end subroutine trustRegionRefinement


  ! *********************************************************************************************************
  ! * internal subroutine eval8summa_wrapper: compute the right-hand-side vector
  ! *********************************************************************************************************
  ! NOTE: This is simply a wrapper routine for eval8summa, to reduce the number of calling arguments
  !       An internal subroutine, so have access to all data in the main subroutine
  subroutine eval8summa_wrapper(stateVecNew,fluxVecNew,resVecNew,fNew,feasible,err,message)
  USE eval8summa_module,only:eval8summa                      ! simulation of fluxes and residuals given a trial state vector
  implicit none
  ! input
  real(dp),intent(in)            :: stateVecNew(:)           ! updated state vector
  ! output
  real(dp),intent(out)           :: fluxVecNew(:)            ! updated flux vector
  real(qp),intent(out)           :: resVecNew(:) ! NOTE: qp  ! updated residual vector
  real(dp),intent(out)           :: fNew                     ! new function value
  logical(lgt),intent(out)       :: feasible                 ! flag to denote the feasibility of the solution
  integer(i4b),intent(out)       :: err                      ! error code
  character(*),intent(out)       :: message                  ! error message
  ! ----------------------------------------------------------------------------------------------------------
  ! local
  character(len=256)             :: cmessage                 ! error message of downwind routine
  ! ----------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='eval8summa_wrapper/'

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
                  canopyDepth,             & ! intent(in):    depth of the vegetation canopy (m)
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
                  ! input-output: baseflow
                  ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                  dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                  ! output
                  feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                  fluxVecNew,              & ! intent(out):   new flux vector
                  resSinkNew,              & ! intent(out):   additional (sink) terms on the RHS of the state equation
                  resVecNew,               & ! intent(out):   new residual vector
                  fNew,                    & ! intent(out):   new function evaluation
                  err,cmessage)              ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)


  end subroutine eval8summa_wrapper


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
  ! TEMPORARY CODE: use hard-coded convergence parameters
  real(dp),parameter        :: absConvTol_energy=1.e-0_dp   ! convergence tolerance for energy (J m-3)
  real(dp),parameter        :: absConvTol_liquid=1.e-8_dp   ! convergence tolerance for volumetric liquid water content (-)
  real(dp),parameter        :: absConvTol_matric=1.e-3_dp   ! convergence tolerance for matric head increment in soil layers (m)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  associate(&
  ! convergence parameters
  !absConvTol_liquid       => mpar_data%var(iLookPARAM%absConvTol_liquid)            ,&  ! intent(in): [dp] absolute convergence tolerance for vol frac liq water (-)
  !absConvTol_matric       => mpar_data%var(iLookPARAM%absConvTol_matric)            ,&  ! intent(in): [dp] absolute convergence tolerance for matric head        (m)
  !absConvTol_energy       => mpar_data%var(iLookPARAM%absConvTol_energy)            ,&  ! intent(in): [dp] absolute convergence tolerance for energy             (J m-3)
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
  end if

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

  ! print progress towards solution
  if(globalPrintFlag)then
   write(*,'(a,1x,i4,1x,4(e15.5,1x),3(i4,1x),5(L1,1x))') 'check convergence: ', iter, &
    fNew, matric_max(1), liquid_max(1), energy_max(1), matric_loc(1), liquid_loc(1), energy_loc(1), matricConv, liquidConv, energyConv, watbalConv, canopyConv
  end if

  ! end associations with variables in the data structures
  end associate

  end function checkConv


  ! *********************************************************************************************************
  ! internal subroutine imposeConstraints: impose solution constraints
  ! *********************************************************************************************************
  subroutine imposeConstraints(stateVecTrial,xInc,err,message)
  USE var_lookup,only:iLookINDEX                                  ! named variables for structure elements
  USE multiconst,only:Tfreeze                                     ! temperature at freezing (K)
  USE multiconst,only:gravity                                     ! acceleration of gravity (m s-2)
  USE multiconst,only:LH_fus                                      ! latent heat of fusion (J kg-1)
  implicit none
  ! dummies
  real(dp),intent(in)             :: stateVecTrial(:)             ! trial state vector
  real(dp),intent(inout)          :: xInc(:)                      ! iteration increment
  integer(i4b),intent(out)        :: err                          ! error code
  character(*),intent(out)        :: message                      ! error message
  ! -----------------------------------------------------------------------------------------------------
  ! local variables
  real(dp)                        :: cInc                         ! constrained temperature increment (K) -- simplified bi-section
  real(dp)                        :: xIncFactor                   ! scaling factor for the iteration increment (-)
  integer(i4b)                    :: iMax(1)                      ! index of maximum temperature
  real(dp),dimension(nSnow)       :: mLayerTempCheck              ! updated temperatures (K) -- used to check iteration increment for snow
  real(dp),dimension(nSnow)       :: mLayerVolFracLiqCheck        ! updated volumetric liquid water content (-) -- used to check iteration increment for snow
  logical(lgt),dimension(nSnow)   :: drainFlag                    ! flag to denote when drainage exceeds available capacity
  logical(lgt),dimension(nSoil)   :: crosFlag                     ! flag to denote temperature crossing from unfrozen to frozen (or vice-versa)
  logical(lgt)                    :: crosTempVeg                  ! flag to denoote where temperature crosses the freezing point
  integer(i4b)                    :: ixNrg,ixLiq                  ! index of energy and mass state variables in full state vector
  integer(i4b)                    :: iLayer                       ! index of model layer
  real(dp)                        :: xPsi00                       ! matric head after applying the iteration increment (m)
  real(dp)                        :: TcSoil                       ! critical point when soil begins to freeze (K)
  real(dp)                        :: critDiff                     ! temperature difference from critical (K)
  real(dp),parameter              :: epsT=1.e-7_dp                ! small interval above/below critical (K)
  real(dp),parameter              :: zMaxTempIncrement=1._dp      ! maximum temperature increment (K)
  ! -----------------------------------------------------------------------------------------------------
  ! associate variables with indices of model state variables
  associate(&
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,&  ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixTopNrg                => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most energy state in the snow-soil subdomain
  ixTopWat                => indx_data%var(iLookINDEX%ixTopWat)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most total water state in the snow-soil subdomain
  ixTopMat                => indx_data%var(iLookINDEX%ixTopMat)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most matric head state in the soil subdomain
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
  ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow-soil subdomain
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
  ixNrgOnly               => indx_data%var(iLookINDEX%ixNrgOnly)%dat                ,&  ! intent(in): [i4b(:)] list of indices for all energy states
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                 &  ! intent(in): [i4b(:)] type of each layer in the snow+soil domain (snow or soil)
  ) ! associating variables with indices of model state variables
  ! -----------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='imposeConstraints/'
  
  ! ** limit temperature increment to zMaxTempIncrement
  !if(any(abs(xInc(ixNrgOnly)) > zMaxTempIncrement))then
  ! iMax       = maxloc( abs(xInc(ixNrgOnly)) )                     ! index of maximum temperature increment
  ! xIncFactor = abs( zMaxTempIncrement/xInc(ixNrgOnly(iMax(1))) )  ! scaling factor for the iteration increment (-)
  ! xInc       = xIncFactor*xInc
  !end if
  
  ! vegetation
  if(computeVegFlux)then
   if(abs(xInc(ixVegNrg)) > 1._dp)then
    xIncFactor = abs(1._dp/xInc(ixVegNrg))  ! scaling factor for the iteration increment (-)
    xInc       = xIncFactor*xInc            ! scale iteration increments
   end if
  end if
  
  ! snow and soil
  if(any(abs(xInc(ixSnowSoilNrg)) > 1._dp))then
   iMax       = maxloc( abs(xInc(ixSnowSoilNrg)) )                   ! index of maximum temperature increment
   xIncFactor = abs( 1._dp/xInc(ixSnowSoilNrg(iMax(1))) )            ! scaling factor for the iteration increment (-)
   xInc       = xIncFactor*xInc
  end if
  
  ! ** impose solution constraints for vegetation
  ! (stop just above or just below the freezing point if crossing)
  if(computeVegFlux)then
  
   ! --------------------------------------------------------------------------------------------------------------------
   ! canopy temperatures
  
   ! initialize
   critDiff    = Tfreeze - stateVecTrial(ixVegNrg)
   crosTempVeg = .false.
  
   ! initially frozen (T < Tfreeze)
   if(critDiff > 0._dp)then
    if(xInc(ixVegNrg) > critDiff)then
     crosTempVeg = .true.
     cInc        = critDiff + epsT  ! constrained temperature increment (K)
    end if
  
   ! initially unfrozen (T > Tfreeze)
   else
    if(xInc(ixVegNrg) < critDiff)then
     crosTempVeg = .true.
     cInc        = critDiff - epsT  ! constrained temperature increment (K)
    end if
  
   end if  ! switch between frozen and unfrozen
  
   ! scale iterations
   if(crosTempVeg)then
    xIncFactor  = cInc/xInc(ixVegNrg)  ! scaling factor for the iteration increment (-)
    xInc        = xIncFactor*xInc      ! scale iteration increments
   end if
  
   ! --------------------------------------------------------------------------------------------------------------------
   ! canopy liquid water
  
   ! check if new value of storage will be negative
   if(stateVecTrial(ixVegWat)+xInc(ixVegWat) < 0._dp)then
    ! scale iteration increment
    cInc       = -0.5_dp*stateVecTrial(ixVegWat)                                  ! constrained iteration increment (K) -- simplified bi-section
    xIncFactor = cInc/xInc(ixVegWat)                                              ! scaling factor for the iteration increment (-)
    xInc       = xIncFactor*xInc                                                  ! new iteration increment
   end if
  
  end if  ! if computing fluxes through vegetation
  
  ! ** impose solution constraints for snow
  if(nSnow > 0)then
  
   ! --------------------------------------------------------------------------------------------------------------------
   ! loop through snow layers
   !checksnow: do iState=1,nSnow  ! necessary to ensure that NO layers rise above Tfreeze
    ! - get updated snow temperatures
    mLayerTempCheck = stateVecTrial(ixSnowOnlyNrg) + xInc(ixSnowOnlyNrg)
    ! - check sub-freezing temperatures for snow
    if(any(mLayerTempCheck > Tfreeze))then
     ! scale iteration increment
     iMax       = maxloc(mLayerTempCheck)                                          ! index of maximum temperature
     cInc       = 0.5_dp*(Tfreeze - stateVecTrial(ixSnowOnlyNrg(iMax(1))) )        ! constrained temperature increment (K) -- simplified bi-section
     xIncFactor = cInc/xInc(ixSnowOnlyNrg(iMax(1)))                                ! scaling factor for the iteration increment (-)
     xInc       = xIncFactor*xInc
   ! else    ! if snow temperature > freezing
   !  exit checkSnow
    end if   ! if snow temperature > freezing
   !end do checkSnow
  
   ! --------------------------------------------------------------------------------------------------------------------
   ! - check if drain more than what is available
   ! NOTE: change in total water is only due to liquid flux
  
   ! get new volumetric fraction of liquid water
   mLayerVolFracLiqCheck = stateVecTrial(ixSnowOnlyWat)+xInc(ixSnowOnlyWat)
   drainFlag(:) = .false.
  
   do iLayer=1,nSnow
    if(mLayerVolFracLiqCheck(iLayer) < 0._dp)then
     drainFlag(iLayer) = .true.
     xInc(ixSnowOnlyWat(iLayer)) = -0.5_dp*stateVecTrial(ixSnowOnlyWat(iLayer) )
    end if
   end do
  
  end if   ! if snow layers exist
  
  ! --------------------------------------------------------------------------------------------------------------------
  ! ** impose solution constraints for soil
  do iLayer=1,nSoil
  
   ! initialize crossing flag
   crosFlag(iLayer) = .false.
  
   ! identify indices for energy and mass state variables
   ixNrg = ixSnowSoilNrg(nSnow+iLayer)
   ixLiq = ixSnowSoilWat(nSnow+iLayer)
  
   ! identify the critical point when soil begins to freeze (TcSoil)
   xPsi00 = stateVecTrial(ixLiq) + xInc(ixLiq)
   TcSoil = Tfreeze + min(xPsi00,0._dp)*gravity*Tfreeze/LH_fus  ! (NOTE: J = kg m2 s-2, so LH_fus is in units of m2 s-2)
  
   ! get the difference from the current state and the crossing point (K)
   critDiff = TcSoil - stateVecTrial(ixNrg)
  
   ! * initially frozen (T < TcSoil)
   if(critDiff > 0._dp)then
  
    ! (check crossing above zero)
    if(xInc(ixNrg) > critDiff)then
     crosFlag(iLayer) = .true.
     xInc(ixNrg) = critDiff + epsT  ! set iteration increment to slightly above critical temperature
    end if
  
   ! * initially unfrozen (T > TcSoil)
   else
  
    ! (check crossing below zero)
    if(xInc(ixNrg) < critDiff)then
     crosFlag(iLayer) = .true.
     xInc(ixNrg) = critDiff - epsT  ! set iteration increment to slightly below critical temperature
    end if
  
   end if  ! (switch between initially frozen and initially unfrozen)
  
   ! place constraint for matric head
   if(xInc(ixLiq) > 1._dp .and. stateVecTrial(ixLiq) > 0._dp)then
    xInc(ixLiq) = 1._dp
   end if  ! if constraining matric head
  
  end do  ! (loop through soil layers)
  
  ! end association with variables with indices of model state variables
  end associate
  
  end subroutine imposeConstraints

 end subroutine summaSolve




end module summaSolve_module
