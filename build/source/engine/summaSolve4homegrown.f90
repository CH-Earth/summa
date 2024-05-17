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

module summaSolve4homegrown_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_snow       ! named variables for snow
USE globalData,only:iname_soil       ! named variables for soil

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! indices of elements of data structure
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

USE multiconst,only:&
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,                        & ! data vector (i4b)
                    var_d,                        & ! data vector (rkind)
                    var_ilength,                  & ! data vector with variable length dimension (i4b)
                    var_dlength,                  & ! data vector with variable length dimension (rkind)
                    zLookup,                      & ! lookup tables
                    model_options,                & ! defines the model decisions
                    in_type_computJacob,          & ! class for computJacob arguments
                    out_type_computJacob,         & ! class for computJacob arguments
                    in_type_lineSearchRefinement, & ! class for lineSearchRefinement arguments
                    out_type_lineSearchRefinement,& ! class for lineSearchRefinement arguments
                    in_type_summaSolve4homegrown, & ! class for summaSolve4homegrown arguments
                    io_type_summaSolve4homegrown, & ! class for summaSolve4homegrown arguments
                    out_type_summaSolve4homegrown   ! class for summaSolve4homegrown arguments


! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
  qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                      & ! a big bucket (lumped aquifer model)
  noExplicit                        ! no explicit groundwater parameterization

! look-up values for method used to compute derivative
USE mDecisions_module,only:       &
  numerical,                      & ! numerical solution
  analytical                        ! analytical solution

implicit none
private
public::summaSolve4homegrown
contains

 ! **************************************************************************************************************************
 ! public subroutine summaSolve4homegrown: calculate the iteration increment, evaluate the new state, and refine if necessary
 ! **************************************************************************************************************************
 subroutine summaSolve4homegrown(&
                        in_SS4HG,               & ! intent(in): model control and previous function value
                       ! input: model control
                       !firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                       ! input: state vectors
                       stateVecTrial,           & ! intent(in):    trial state vector
                       !xMin,xMax,               & ! intent(inout): brackets of the root
                       fScale,                  & ! intent(in):    characteristic scale of the function evaluations
                       xScale,                  & ! intent(in):    characteristic scale of the state vector
                       rVec,                    & ! intent(in):    residual vector
                       sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                       dMat,                    & ! intent(inout): diagonal matrix (excludes flux derivatives)
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
                       flux_data,               & ! intent(inout): model fluxes for a local HRU
                       deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ! input-output: baseflow
                       !ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       dBaseflow_dMatric,       & ! intent(inout): derivative in baseflow w.r.t. matric head (s-1)
                       io_SS4HG,                & ! intent(inout): first flux call flag, root brackets, index of lowest saturated layer
                       ! output
                       stateVecNew,             & ! intent(out):   new state vector
                       fluxVecNew,              & ! intent(out):   new flux vector
                       resSinkNew,              & ! intent(out):   additional (sink) terms on the RHS of the state equation
                       resVecNew,               & ! intent(out):   new residual vector
                       out_SS4HG)                 ! intent(out):   new function evaluation, convergence flag, and error control  
 USE computJacob_module, only: computJacob
 USE eval8summa_module,  only: imposeConstraints
 USE matrixOper_module,  only: lapackSolv
 USE matrixOper_module,  only: scaleMatrices
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 type(in_type_summaSolve4homegrown),intent(in)     :: in_SS4HG  ! model control variables and previous function evaluation
 type(io_type_summaSolve4homegrown),intent(inout)  :: io_SS4HG  ! first flux call flag and baseflow variables
 ! input: model control
! logical(lgt),intent(inout)      :: firstFluxCall            ! flag to indicate if we are processing the first flux call
 ! input: state vectors
 real(rkind),intent(in)          :: stateVecTrial(:)         ! trial state vector
! real(rkind),intent(inout)       :: xMin,xMax                ! brackets of the root
 real(rkind),intent(in)          :: fScale(:)                ! characteristic scale of the function evaluations
 real(rkind),intent(in)          :: xScale(:)                ! characteristic scale of the state vector
 real(qp),intent(in)             :: rVec(:)   ! NOTE: qp     ! residual vector
 real(qp),intent(inout)          :: sMul(:)   ! NOTE: qp     ! state vector multiplier (used in the residual calculations)
 real(rkind),intent(inout)       :: dMat(:)                  ! diagonal matrix (excludes flux derivatives)
 ! input: data structures
 type(model_options),intent(in)  :: model_decisions(:)       ! model decisions
 type(zLookup),      intent(in)  :: lookup_data              ! lookup tables
 type(var_i),        intent(in)  :: type_data                ! type of vegetation and soil
 type(var_d),        intent(in)  :: attr_data                ! spatial attributes
 type(var_dlength),  intent(in)  :: mpar_data                ! model parameters
 type(var_d),        intent(in)  :: forc_data                ! model forcing data
 type(var_dlength),  intent(in)  :: bvar_data                ! model variables for the local basin
 type(var_dlength),  intent(in)  :: prog_data                ! prognostic variables for a local HRU
 ! output: data structures
 type(var_ilength),intent(inout) :: indx_data                ! indices defining model states and layers
 type(var_dlength),intent(inout) :: diag_data                ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: deriv_data               ! derivatives in model fluxes w.r.t. relevant state variables
 ! input-output: baseflow
! integer(i4b),intent(inout)      :: ixSaturation             ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(rkind),intent(inout)       :: dBaseflow_dMatric(:,:)   ! derivative in baseflow w.r.t. matric head (s-1)
 ! output: flux and residual vectors
 real(rkind),intent(out)         :: stateVecNew(:)           ! new state vector
 real(rkind),intent(out)         :: fluxVecNew(:)            ! new flux vector
 real(rkind),intent(out)         :: resSinkNew(:)            ! sink terms on the RHS of the flux equation
 real(qp),intent(out)            :: resVecNew(:) ! NOTE: qp  ! new residual vector
 type(out_type_summaSolve4homegrown),intent(out)  :: out_SS4HG  ! new function evaluation, convergence flag, and error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! Jacobian matrix
 logical(lgt),parameter          :: doNumJacobian=.false.    ! flag to compute the numerical Jacobian matrix
 logical(lgt),parameter          :: testBandDiagonal=.false. ! flag to test the band diagonal Jacobian matrix
 real(rkind)                     :: nJac(in_SS4HG % nState,in_SS4HG % nState)      ! numerical Jacobian matrix
 real(rkind)                     :: aJac(in_SS4HG % nLeadDim,in_SS4HG % nState)      ! Jacobian matrix
 real(rkind)                     :: aJacScaled(in_SS4HG % nLeadDim,in_SS4HG % nState)      ! Jacobian matrix (scaled)
 real(rkind)                     :: aJacScaledTemp(in_SS4HG % nLeadDim,in_SS4HG % nState)  ! Jacobian matrix (scaled) -- temporary copy since decomposed in lapack
 ! solution/step vectors
 real(rkind),dimension(in_SS4HG % nState)   :: rVecScaled               ! residual vector (scaled)
 real(rkind),dimension(in_SS4HG % nState)   :: newtStepScaled           ! full newton step (scaled)
 ! step size refinement
 logical(lgt)                    :: doRefine                 ! flag for step refinement
 integer(i4b),parameter          :: ixLineSearch=1001        ! step refinement = line search
 integer(i4b),parameter          :: ixTrustRegion=1002       ! step refinement = trust region
 integer(i4b),parameter          :: ixStepRefinement=ixLineSearch   ! decision for the numerical solution
 ! general
 integer(i4b)                    :: mSoil                    ! number of soil layers in solution vector
 integer(i4b)                    :: iLayer                   ! row index
 integer(i4b)                    :: jLayer                   ! column index
 logical(lgt)                    :: globalPrintFlagInit      ! initial global print flag
 logical(lgt)                    :: return_flag              ! flag that controls execution of return statements
 character(LEN=256)              :: cmessage                 ! error message of downwind routine
 ! class objects for subroutine arguments
 type(in_type_computJacob)           :: in_computJacob  ! computJacob
 type(out_type_computJacob)          :: out_computJacob ! computJacob 
 type(in_type_lineSearchRefinement)  :: in_LSR  ! lineSearchRefinement
 type(out_type_lineSearchRefinement) :: out_LSR ! lineSearchRefinement 
 type(in_type_lineSearchRefinement)  :: in_TRR  ! trustRegionRefinement
 type(out_type_lineSearchRefinement) :: out_TRR ! trustRegionRefinement
 type(out_type_lineSearchRefinement) :: out_SRF ! safeRootFinder
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! ***** Compute the Newton Step *****

 call initialize_summaSolve4homegrown; if (return_flag) return ! initial setup -- return if error

 call update_Jacobian;                 if (return_flag) return ! compute Jacobian for Newton step -- return if error

 call solve_linear_system;             if (return_flag) return ! solve the linear system for the Newton step -- return if error

 call refine_Newton_step;              if (return_flag) return ! refine Newton step if needed -- return if error

 associate(err => out_SS4HG % err,message   => out_SS4HG % message)
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! check for errors
 end associate

 contains

  subroutine initialize_summaSolve4homegrown
   ! *** Initial steps for the summaSolve4homegrown algorithm (computing the Newton step) ***

   associate(&
    err       => out_SS4HG % err      ,& 
    message   => out_SS4HG % message   &     
   &)
    ! initialize error control
    err=0; message='summaSolve4homegrown/'
    return_flag=.false. ! initialize return flag
   
    ! choose Jacobian type
    select case(model_decisions(iLookDECISIONS%fDerivMeth)%iDecision) 
     case(numerical)
      err=20; message=trim(message)//'numerical derivatives are not implemented for BE homegrown solver';
      return_flag=.true.; return
     case(analytical); ! this is fine
     case default
      err=20; message=trim(message)//'expect choice numericl or analytic to calculate derivatives for Jacobian';
      return_flag=.true.; return
    end select
  
   end associate
 
   ! get the number of soil layers in the solution vector
   mSoil = size(indx_data%var(iLookINDEX%ixMatOnly)%dat)
  
   ! initialize the global print flag
   globalPrintFlagInit=globalPrintFlag
  end subroutine initialize_summaSolve4homegrown

  subroutine update_Jacobian
   ! *** Update Jacobian used for Newton step ***
  
   ! compute the analytical Jacobian matrix
   ! NOTE: The derivatives were computed in the previous call to computFlux
   !       This occurred either at the call to eval8summa at the start of systemSolv
   !        or in the call to eval8summa in the previous iteration (within lineSearchRefinement or trustRegionRefinement)
   associate(&
    err       => out_SS4HG % err      ,& 
    message   => out_SS4HG % message   &     
   &)
    call initialize_computJacob_summaSolve4homegrown
    call computJacob(in_computJacob,indx_data,prog_data,diag_data,deriv_data,dBaseflow_dMatric,dMat,aJac,out_computJacob)
    call finalize_computJacob_summaSolve4homegrown
    if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! (check for errors)
  
    ! compute the numerical Jacobian matrix
    if (doNumJacobian) then
     globalPrintFlag=.false.
     call numJacobian(stateVecTrial,dMat,nJac,err,cmessage)
     if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! (check for errors)
     globalPrintFlag=globalPrintFlagInit
    end if
   end associate
  end subroutine update_Jacobian

  subroutine solve_linear_system
   ! *** Solve the linear system for the Newton step using LAPACK routines ***

   associate(&
    nState         => in_SS4HG % nState         ,& ! intent(in): total number of state variables
    ixMatrix       => in_SS4HG % ixMatrix       ,& ! intent(in): type of matrix (full or band diagonal)
    err            => out_SS4HG % err           ,& ! intent(out): error code
    message        => out_SS4HG % message        & ! intent(out): error message    
   &)
    ! test the band diagonal matrix
    if (testBandDiagonal) then
     call testBandMat(check=.true.,err=err,message=cmessage)
     if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! (check for errors)
    end if
 
    ! scale the residual vector
    rVecScaled(1:nState) = fScale(:)*real(rVec(:), rkind)   ! NOTE: residual vector is in quadruple precision
   
    ! scale matrices
    call scaleMatrices(ixMatrix,nState,aJac,fScale,xScale,aJacScaled,err,cmessage)
    if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! check for errors
  
    ! debug statement for scaled Jacobian 
    if (globalPrintFlag .and. ixMatrix==ixBandMatrix) then
     print*, '** SCALED banded analytical Jacobian:'
     write(*,'(a4,1x,100(i17,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
     do iLayer=kl+1,nBands
      write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJacScaled(iLayer,jLayer),jLayer=min(iJac1,nState),min(iJac2,nState))
     end do
    end if
   
    ! copy the scaled matrix, since it is decomposed in lapackSolv
    aJacScaledTemp = aJacScaled
   
    ! compute the newton step: use the lapack routines to solve the linear system A.X=B
    call lapackSolv(ixMatrix,nState,aJacScaledTemp,-rVecScaled,newtStepScaled,err,cmessage)
    if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! check for errors
 
    ! debug statement for Newton step
    if (globalPrintFlag) write(*,'(a,1x,10(e17.10,1x))') 'newtStepScaled = ', newtStepScaled(min(iJac1,nState):min(iJac2,nState))

   end associate
  end subroutine solve_linear_system

  subroutine refine_Newton_step  
   ! *** Refine the Newton step if necessary ***  
  
   ! initialize the flag for step refinement
   doRefine=.true.
  
   ! * case 1: state vector
   ! compute the flux vector and the residual, and (if necessary) refine the iteration increment
   ! NOTE: in 99.9% of cases newtStep will be used (no refinement)

   associate(&
    fOld      => in_SS4HG  % fOld     ,&
    fnew      => out_SS4HG % fnew     ,& 
    converged => out_SS4HG % converged,&
    err       => out_SS4HG % err      ,& 
    message   => out_SS4HG % message   &     
   &)  
    if (size(stateVecTrial)>1) then

     ! try to backtrack
     select case(ixStepRefinement)
      case(ixLineSearch)  
       call in_LSR % initialize(doRefine,fOld)    
       call lineSearchRefinement(in_LSR,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,stateVecNew,fluxVecNew,resVecNew,out_LSR)
       call out_LSR % finalize(fNew,converged,err,cmessage)
      case(ixTrustRegion)
       call in_TRR % initialize(doRefine,fOld)
       call trustRegionRefinement(in_TRR,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,stateVecNew,fluxVecNew,resVecNew,out_TRR)
       call out_TRR % finalize(fNew,converged,err,cmessage)
      case default; err=20; message=trim(message)//'unable to identify numerical solution'; return_flag=.true.; return
     end select
  
     ! check warnings: negative error code = warning; in this case back-tracked to the original value
     ! NOTE: Accept the full newton step if back-tracked to the original value
     if (err<0) then
      doRefine=.false.;
      call in_LSR % initialize(doRefine,fOld)    
      call lineSearchRefinement(in_LSR,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,stateVecNew,fluxVecNew,resVecNew,out_LSR)
      call out_LSR % finalize(fNew,converged,err,cmessage)
     end if
  
    ! * case 2: scalar
    else
     call safeRootfinder(stateVecTrial,rVecScaled,newtStepScaled,stateVecNew,fluxVecNew,resVecNew,out_SRF)
     call out_SRF % finalize(fNew,converged,err,cmessage)
     if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! check for errors
    end if
   end associate
  end subroutine refine_Newton_step  

  ! *********************************************************************************************************
  ! * internal subroutine lineSearchRefinement: refine the iteration increment using line searches
  ! *********************************************************************************************************
  subroutine lineSearchRefinement(in_LSR,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,stateVecNew,fluxVecNew,resVecNew,out_LSR)
  ! provide access to the matrix routines
  USE matrixOper_module, only: computeGradient
  implicit none
  ! input
  type(in_type_lineSearchRefinement),intent(in)  :: in_LSR   ! class object for intent(in) arguments
  real(rkind),intent(in)         :: stateVecTrial(:)         ! trial state vector
  real(rkind),intent(in)         :: newtStepScaled(:)        ! scaled newton step
  real(rkind),intent(in)         :: aJacScaled(:,:)          ! scaled jacobian matrix
  real(rkind),intent(in)         :: rVecScaled(:)            ! scaled residual vector
  ! output
  real(rkind),intent(out)        :: stateVecNew(:)           ! new state vector
  real(rkind),intent(out)        :: fluxVecNew(:)            ! new flux vector
  real(qp),intent(out)           :: resVecNew(:) ! NOTE: qp  ! new residual vector
  type(out_type_lineSearchRefinement),intent(out) :: out_LSR ! class object for intent(out) arguments
  ! --------------------------------------------------------------------------------------------------------
  ! local
  character(len=256)             :: cmessage                      ! error message of downwind routine
  real(rkind)                    :: gradScaled(in_SS4HG % nState) ! scaled gradient
  real(rkind)                    :: xInc(in_SS4HG % nState)       ! iteration increment (re-scaled to original units of the state vector)
  logical(lgt)                   :: feasible                      ! flag to denote the feasibility of the solution
  integer(i4b)                   :: iLine                         ! line search index
  integer(i4b),parameter         :: maxLineSearch=5               ! maximum number of backtracks
  real(rkind),parameter          :: alpha=1.e-4_rkind             ! check on gradient
  real(rkind)                    :: xLambda                       ! backtrack magnitude
  real(rkind)                    :: xLambdaTemp                   ! temporary backtrack magnitude
  real(rkind)                    :: slopeInit                     ! initial slope
  real(rkind)                    :: rhs1,rhs2                     ! rhs used to compute the cubic
  real(rkind)                    :: aCoef,bCoef                   ! coefficients in the cubic
  real(rkind)                    :: disc                          ! temporary variable used in cubic
  real(rkind)                    :: xLambdaPrev                   ! previous lambda value (used in the cubic)
  real(rkind)                    :: fPrev                         ! previous function evaluation (used in the cubic)
  ! --------------------------------------------------------------------------------------------------------
  associate(&
   ! intent(in) variables
   doLineSearch => in_LSR % doSearch            ,& ! flag to do the line search
   fOld         => in_LSR % fOld                ,& ! old function value
   ! local variables
   nSnow        => in_SS4HG % nSnow             ,& ! number of snow layers
   nSoil        => in_SS4HG % nSoil             ,& ! number of soil layers
   nState       => in_SS4HG % nState            ,& ! total number of state variables
   ixMatrix     => in_SS4HG % ixMatrix          ,& ! type of matrix (full or band diagonal)
   ! intent(out) variables
   fNew      => out_LSR % fNew                  ,& ! new function evaluation
   converged => out_LSR % converged             ,& ! convergence flag
   err       => out_LSR % err                   ,& ! error code
   message   => out_LSR % message                & ! error message
  &)
   ! initialize error control
   err=0; message='lineSearchRefinement/'
   converged =.false.

   ! check the need to compute the line search
   if(doLineSearch)then

    ! compute the gradient of the function vector
    call computeGradient(ixMatrix,nState,aJacScaled,rVecScaled,gradScaled,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! compute the initial slope
    slopeInit = dot_product(gradScaled,newtStepScaled)

   end if  ! if computing the line search

   ! initialize lambda
   xLambda=1._rkind

   ! ***** LINE SEARCH LOOP...
   lineSearch: do iLine=1,maxLineSearch  ! try to refine the function by shrinking the step size

    ! back-track along the search direction
    ! NOTE: start with back-tracking the scaled step
    xInc(:) = xLambda*newtStepScaled(:)

    ! re-scale the iteration increment
    xInc(:) = xInc(:)*xScale(:)

    ! if enthalpy, then need to convert the iteration increment to temperature
    !if(nrgFormulation==ix_enthalpy) xInc(ixNrgOnly) = xInc(ixNrgOnly)/dMat(ixNrgOnly)

    ! state vector with proposed iteration increment
    stateVecNew = stateVecTrial + xInc

    ! impose solution constraints adjusting state vector and iteration increment
    ! NOTE: We may not need to do this (or at least, do ALL of this), as we can probably rely on the line search here
    call imposeConstraints(model_decisions,indx_data,prog_data,mpar_data,stateVecNew,stateVecTrial,nState,nSoil,nSnow,cmessage,err)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
    xInc = stateVecNew - stateVecTrial

    ! compute the residual vector and function
    ! NOTE: This calls eval8summa in an internal subroutine which has access to all data
    !       Hence, we only need to include the variables of interest in lineSearch
    call eval8summa_wrapper(stateVecNew,fluxVecNew,resVecNew,fNew,feasible,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! check line search
    if(globalPrintFlag)then
     write(*,'(a,1x,i4,1x,e17.10)' ) 'iLine, xLambda                 = ', iLine, xLambda
     write(*,'(a,1x,10(e17.10,1x))') 'fOld,fNew                      = ', fOld,fNew
     write(*,'(a,1x,10(e17.10,1x))') 'fold + alpha*slopeInit*xLambda = ', fold + alpha*slopeInit*xLambda
     write(*,'(a,1x,10(e17.10,1x))') 'resVecNew                      = ', resVecNew(min(iJac1,nState):min(iJac2,nState))
     write(*,'(a,1x,10(e17.10,1x))') 'xInc                           = ', xInc(min(iJac1,nState):min(iJac2,nState))
    end if

    ! check feasibility
    if(.not.feasible) cycle ! go back and impose constraints again

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
     xLambdaTemp = -slopeInit / (2._rkind*(fNew - fOld - slopeInit) )
     if(xLambdaTemp > 0.5_rkind*xLambda) xLambdaTemp = 0.5_rkind*xLambda

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
     if(aCoef==0._rkind)then
      xLambdaTemp = -slopeInit/(2._rkind*bCoef)

     ! calculate cubic
     else
      disc = bCoef*bCoef - 3._rkind*aCoef*slopeInit
      if(disc < 0._rkind)then
       xLambdaTemp = 0.5_rkind*xLambda
      else
       xLambdaTemp = (-bCoef + sqrt(disc))/(3._rkind*aCoef)
      end if
     end if  ! calculating cubic

     ! constrain to <= 0.5*xLambda
     if(xLambdaTemp > 0.5_rkind*xLambda) xLambdaTemp=0.5_rkind*xLambda

    end if  ! subsequent backtracks

    ! save results
    xLambdaPrev = xLambda
    fPrev = fNew

    ! constrain lambda
    xLambda = max(xLambdaTemp, 0.1_rkind*xLambda)

   end do lineSearch  ! backtrack loop
  end associate

  end subroutine lineSearchRefinement


  ! *********************************************************************************************************
  ! * internal subroutine trustRegionRefinement: refine the iteration increment using trust regions
  ! *********************************************************************************************************
  subroutine trustRegionRefinement(in_TRR,stateVecTrial,newtStepScaled,aJacScaled,rVecScaled,stateVecNew,fluxVecNew,resVecNew,out_TRR)
  ! provide access to the matrix routines
  USE matrixOper_module, only: lapackSolv
  USE matrixOper_module, only: computeGradient
  implicit none
  ! input
  type(in_type_lineSearchRefinement),intent(in) :: in_TRR            ! object for scalar intent(in) arguments -- reusing line search class
  real(rkind),intent(in)                        :: stateVecTrial(:)  ! trial state vector
  real(rkind),intent(in)                        :: newtStepScaled(:) ! scaled newton step
  real(rkind),intent(in)                        :: aJacScaled(:,:)   ! scaled jacobian matrix
  real(rkind),intent(in)                        :: rVecScaled(:)     ! scaled residual vector
  ! output
  real(rkind),intent(out)                       :: stateVecNew(:)    ! new state vector
  real(rkind),intent(out)                       :: fluxVecNew(:)     ! new flux vector
  real(qp),intent(out)                          :: resVecNew(:)      ! NOTE: qp  ! new residual vector
  type(out_type_lineSearchRefinement)           :: out_TRR           ! object for scalar intent(in) arguments -- reusing line search class
  ! --------------------------------------------------------------------------------------------------------
  ! local variables

  ! .. needed ..

  ! --------------------------------------------------------------------------------------------------------
  associate(&
   ! input
   doTrustRefinement => in_TRR % doSearch          ,&    ! flag to refine using trust regions
   fOld              => in_TRR % fOld              ,&    ! old function value
   nState            => in_SS4HG % nState          ,&    ! total number of state variables
   ! output
   fNew      => out_TRR % fNew                     ,&    ! new function evaluation
   converged => out_TRR % converged                ,&    ! convergence flag
   err       => out_TRR % err                      ,&    ! error code
   message   => out_TRR % message                   &    ! error message
  &)

   err=0; message='trustRegionRefinement/'
   converged =.false.

   ! check the need to refine the step
   if (doTrustRefinement) then

    ! (check vectors)
    if (size(stateVecTrial)/=nState .or. size(newtStepScaled)/=nState .or. size(rVecScaled)/=nState)then
     message=trim(message)//'unexpected size of input vectors'
     err=20; return
    end if

    ! (check matrix)
    if (size(aJacScaled,1)/=nState .or. size(aJacScaled,2)/=nState) then
     message=trim(message)//'unexpected size of Jacobian matrix'
     err=20; return
    end if

    ! dummy check for the function
    if (fold==realMissing) print*, 'missing'

    ! dummy
    stateVecNew = realMissing
    fluxVecNew  = realMissing
    resVecNew   = quadMissing
    fNew        = realMissing
    converged   = .true.


   end if  ! if doing the trust region refinement

   message=trim(message)//'routine not implemented yet'
   err=20; return

  end associate

  end subroutine trustRegionRefinement


  ! *********************************************************************************************************
  ! * internal subroutine safeRootfinder: refine the 1-d iteration increment using brackets
  ! *********************************************************************************************************
  subroutine safeRootfinder(stateVecTrial,rVecscaled,newtStepScaled,stateVecNew,fluxVecNew,resVecNew,out_SRF)
  USE,intrinsic :: ieee_arithmetic,only:ieee_is_nan          ! IEEE arithmetic (check NaN)
  USE globalData,only:dNaN                                   ! double precision NaN
  implicit none
  ! input
  real(rkind),intent(in)         :: stateVecTrial(:)         ! trial state vector
  real(rkind),intent(in)         :: rVecScaled(:)            ! scaled residual vector
  real(rkind),intent(in)         :: newtStepScaled(:)        ! scaled newton step
  ! output
  real(rkind),intent(out)        :: stateVecNew(:)           ! new state vector
  real(rkind),intent(out)        :: fluxVecNew(:)            ! new flux vector
  real(qp),intent(out)           :: resVecNew(:) ! NOTE: qp  ! new residual vector
  type(out_type_lineSearchRefinement),intent(out) :: out_SRF ! object for scalar intent(out) arguments (reusing lineSearchRefinement class)
  ! --------------------------------------------------------------------------------------------------------
  ! local variables
  character(len=256)             :: cmessage                 ! error message of downwind routine
  real(rkind),parameter          :: relTolerance=0.005_rkind ! force bi-section if trial is slightly larger than (smaller than) xmin (xmax)
  real(rkind)                    :: xTolerance               ! relTolerance*(xmax-xmin)
  real(rkind)                    :: xInc(in_SS4HG % nState)  ! iteration increment (re-scaled to original units of the state vector)
  real(rkind)                    :: rVec(in_SS4HG % nState)  ! residual vector (re-scaled to original units of the state equation)
  logical(lgt)                   :: feasible                 ! feasibility of the solution
  logical(lgt)                   :: doBisection              ! flag to do the bi-section
  logical(lgt)                   :: bracketsDefined          ! flag to define if the brackets are defined
  integer(i4b),parameter         :: nCheck=100               ! number of times to check the model state variables
  real(rkind),parameter          :: delX=1._rkind            ! trial increment
  ! --------------------------------------------------------------------------------------------------------
  associate(&
   iter           => in_SS4HG % iter           ,& ! intent(in): iteration index
   nSnow          => in_SS4HG % nSnow          ,& ! intent(in): number of snow layers
   nSoil          => in_SS4HG % nSoil          ,& ! intent(in): number of soil layers
   nState         => in_SS4HG % nState         ,& ! intent(in): total number of state
   xMin           => io_SS4HG % xMin           ,& ! intent(inout): bracket of the root   
   xMax           => io_SS4HG % xMax           ,& ! intent(inout): bracket of the root  
   fNew           => out_SRF % fNew            ,& ! intent(out): new function evaluation
   converged      => out_SRF % converged       ,& ! intent(out): convergence flag
   err            => out_SRF % err             ,& ! intent(out): error code
   message        => out_SRF % message          & ! intent(out): error message
  &)

   err=0; message='safeRootfinder/'
   converged = .false.

   ! check scalar
   if (size(stateVecTrial)/=1 .or. size(rVecScaled)/=1 .or. size(newtStepScaled)/=1) then
    message=trim(message)//'unexpected size of input vectors'
    err=20; return
   end if

   ! initialize brackets to rkind precision NaN
   if (iter==1) then
    xMax = dNaN
    xMin = dNaN
   end if

   ! get the residual vector
   rVec = real(rVecScaled, rkind)*real(fScale, rkind)

   ! update brackets
   if (rVec(1)<0._rkind) then
    xMin = stateVecTrial(1)
   else
    xMax = stateVecTrial(1)
   end if

   ! get the iteration increment
   xInc = newtStepScaled*xScale

   ! *****
   ! * case 1: the iteration increment is the same sign as the residual vector
   if (xInc(1)*rVec(1) > 0._rkind) then

    ! get brackets if they do not exist
    if ( ieee_is_nan(xMin) .or. ieee_is_nan(xMax) ) then
     call getBrackets(stateVecTrial,stateVecNew,xMin,xMax,err,cmessage)
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! check for errors
    end if

    ! use bi-section
    stateVecNew(1) = 0.5_rkind*(xMin + xMax)

   ! *****
   ! * case 2: the iteration increment is the correct sign
   else

    ! state vector with proposed iteration increment
    stateVecNew = stateVecTrial + xInc
     
    ! impose solution constraints adjusting state vector and iteration increment
    call imposeConstraints(model_decisions,indx_data,prog_data,mpar_data,stateVecNew,stateVecTrial,nState,nSoil,nSnow,cmessage,err)
    if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! check for errors
    xInc = stateVecNew - stateVecTrial

   end if  ! if the iteration increment is the same sign as the residual vector

   ! bi-section
   bracketsDefined = ( .not.ieee_is_nan(xMin) .and. .not.ieee_is_nan(xMax) )  ! check that the brackets are defined
   if (bracketsDefined) then
    xTolerance  = relTolerance*(xMax-xMin)
    doBisection = (stateVecNew(1)<xMin+xTolerance .or. stateVecNew(1)>xMax-xTolerance)
    if (doBisection) stateVecNew(1) = 0.5_rkind*(xMin+xMax)
   end if

   ! evaluate summa
   call eval8summa_wrapper(stateVecNew,fluxVecNew,resVecNew,fNew,feasible,err,cmessage)
   if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! check for errors

   ! check feasibility (should be feasible because of the call to imposeConstraints, except if canopyTemp>canopyTempMax (500._rkind)) 
   if (.not.feasible) then; err=20; message=trim(message)//'state vector not feasible'; return; end if

   ! check convergence
   converged = checkConv(resVecNew,xInc,stateVecNew)

  end associate

  end subroutine safeRootfinder

  ! *********************************************************************************************************
  ! * internal subroutine getBrackets: get the brackets
  ! *********************************************************************************************************
  subroutine getBrackets(stateVecTrial,stateVecNew,xMin,xMax,err,message)
  USE,intrinsic :: ieee_arithmetic,only:ieee_is_nan          ! IEEE arithmetic (check NaN)
  implicit none
  ! dummies
  real(rkind),intent(in)         :: stateVecTrial(:)         ! trial state vector
  real(rkind),intent(out)        :: stateVecNew(:)           ! new state vector
  real(rkind),intent(out)        :: xMin,xMax                ! constraints
  integer(i4b),intent(inout)     :: err                      ! error code
  character(*),intent(out)       :: message                  ! error message
  ! locals
  real(rkind)                    :: stateVecPrev(in_SS4HG % nState) ! iteration state vector
  integer(i4b)                   :: iCheck                          ! check the model state variables
  integer(i4b),parameter         :: nCheck=100                      ! number of times to check the model state variables
  logical(lgt)                   :: feasible                        ! feasibility of the solution
  real(rkind),parameter          :: delX=1._rkind                   ! trial increment
  real(rkind)                    :: xIncrement(in_SS4HG % nState)   ! trial increment
  ! initialize
  err=0; message='getBrackets/'

  ! initialize state vector
  stateVecNew = stateVecTrial
  stateVecPrev = stateVecNew

  ! get xIncrement
  xIncrement = -sign((/delX/),rVec)

  ! try the increment a few times
  do iCheck=1,nCheck

   ! state vector with proposed iteration increment
   stateVecNew = stateVecPrev + xIncrement

   ! impose solution constraints adjusting state vector and iteration increment
   associate(&
    nSnow          => in_SS4HG % nSnow          ,& ! intent(in): number of snow layers
    nSoil          => in_SS4HG % nSoil          ,& ! intent(in): number of soil layers
    nState         => in_SS4HG % nState          & ! intent(in): total number of state variables
   &)
    call imposeConstraints(model_decisions,indx_data,prog_data,mpar_data,stateVecNew,stateVecPrev,nState,nSoil,nSnow,cmessage,err)
   end associate
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
   xIncrement = stateVecNew - stateVecPrev

   ! evaluate summa
   associate(fnew => out_SS4HG % fnew)
    call eval8summa_wrapper(stateVecNew,fluxVecNew,resVecNew,fNew,feasible,err,cmessage)
   end associate
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

   ! check feasibility (should be feasible because of the call to imposeConstraints, except if canopyTemp>canopyTempMax (500._rkind)) 
   if(.not.feasible)then; message=trim(message)//'state vector not feasible'; err=20; return; endif

   ! update brackets
   if(real(resVecNew(1), rkind)<0._rkind)then
    xMin = stateVecNew(1)
   else
    xMax = stateVecNew(1)
   endif

   ! check that the brackets are defined
   if( .not.ieee_is_nan(xMin) .and. .not.ieee_is_nan(xMax) ) exit

   ! check that we found the brackets
   if(iCheck==nCheck)then
    message=trim(message)//'could not fix the problem where residual and iteration increment are of the same sign'
    err=20; return
   endif

   ! Save the state vector
    stateVecPrev = stateVecNew

  end do  ! multiple checks

  end subroutine getBrackets


  ! *********************************************************************************************************
  ! * internal subroutine numJacobian: compute the numerical Jacobian matrix
  ! *********************************************************************************************************
  subroutine numJacobian(stateVec,dMat,nJac,err,message)
  implicit none
  ! dummies
  real(rkind),intent(in)         :: stateVec(:)                ! trial state vector
  real(rkind),intent(in)         :: dMat(:)                    ! diagonal matrix
  ! output
  real(rkind),intent(out)        :: nJac(:,:)                  ! numerical Jacobian
  integer(i4b),intent(out)       :: err                        ! error code
  character(*),intent(out)       :: message                    ! error message
  ! ----------------------------------------------------------------------------------------------------------
  ! local
  character(len=256)             :: cmessage                   ! error message of downwind routine
  real(rkind),parameter          :: dx=1.e-8_rkind             ! finite difference increment
  real(rkind),dimension(in_SS4HG % nState)  :: stateVecPerturbed          ! perturbed state vector
  real(rkind),dimension(in_SS4HG % nState)  :: fluxVecInit,fluxVecJac     ! flux vector (mized units)
  real(qp),dimension(in_SS4HG % nState)     :: resVecInit,resVecJac ! qp  ! residual vector (mixed units)
  real(rkind)                    :: func                       ! function value
  logical(lgt)                   :: feasible                   ! flag to denote the feasibility of the solution
  integer(i4b)                   :: iJac                       ! index of row of the Jacobian matrix
  integer(i4b),parameter         :: ixNumFlux=1001             ! named variable for the flux-based form of the numerical Jacobian
  integer(i4b),parameter         :: ixNumRes=1002              ! named variable for the residual-based form of the numerical Jacobian
  integer(i4b)                   :: ixNumType=ixNumRes         ! method used to calculate the numerical Jacobian
  ! ----------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='numJacobian/'

  ! compute initial function evaluation
  call eval8summa_wrapper(stateVec,fluxVecInit,resVecInit,func,feasible,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
  if(.not.feasible)then; message=trim(message)//'initial state vector not feasible'; err=20; return; endif

  ! get a copy of the state vector to perturb
  stateVecPerturbed(:) = stateVec(:)

  ! loop through state variables
  do iJac=1,in_SS4HG % nState

   !print*, 'iJac = ', iJac
   !globalPrintFlag = merge(.true.,.false., iJac==1)

   ! perturb state vector
   stateVecPerturbed(iJac) = stateVec(iJac) + dx

   ! compute function evaluation
   call eval8summa_wrapper(stateVecPerturbed,fluxVecJac,resVecJac,func,feasible,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
   if(.not.feasible)then; message=trim(message)//'state vector not feasible'; err=20; return; endif
   !write(*,'(a,1x,2(f30.20,1x))') 'resVecJac(101:102)  = ', resVecJac(101:102)

   ! compute the row of the Jacobian matrix
   select case(ixNumType)
    case(ixNumRes);  nJac(:,iJac) = real(resVecJac - resVecInit, kind(rkind) )/dx  ! Jacobian based on residuals
    case(ixNumFlux); nJac(:,iJac) = -in_SS4HG % dt_cur*(fluxVecJac(:) - fluxVecInit(:))/dx     ! Jacobian based on fluxes
    case default; err=20; message=trim(message)//'Jacobian option not found'; return
   end select

   ! if flux option then add in the diagonal matrix
   if(ixNumType==ixNumFlux) nJac(iJac,iJac) = nJac(iJac,iJac) + dMat(iJac)

   ! set the state back to the input value
   stateVecPerturbed(iJac) = stateVec(iJac)

  end do  ! (looping through state variables)

  ! print the Jacobian
  print*, '** numerical Jacobian:', ixNumType==ixNumRes
  write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=min(iJac1,in_SS4HG % nState),min(iJac2,in_SS4HG % nState))
  do iLayer=min(iJac1,in_SS4HG % nState),min(iJac2,in_SS4HG % nState)
   write(*,'(i4,1x,100(e12.5,1x))') iLayer, nJac(min(iJac1,in_SS4HG % nState):min(iJac2,in_SS4HG % nState),iLayer)
  end do
  !print*, 'PAUSE: testing Jacobian'; read(*,*)

  end subroutine numJacobian

  ! *********************************************************************************************************
  ! * internal subroutine testBandMat: compute the full Jacobian matrix and decompose into a band matrix
  ! *********************************************************************************************************

  subroutine testBandMat(check,err,message)
  ! dummy variables
  logical(lgt),intent(in)         :: check                    ! flag to pause
  integer(i4b),intent(out)        :: err                      ! error code
  character(*),intent(out)        :: message                  ! error message
  ! local variables
  real(rkind)                     :: fullJac(in_SS4HG % nState,in_SS4HG % nState)   ! full Jacobian matrix
  real(rkind)                     :: bandJac(in_SS4HG % nLeadDim,in_SS4HG % nState) ! band Jacobian matrix
  integer(i4b)                    :: iState,jState            ! indices of the state vector
  character(LEN=256)              :: cmessage                 ! error message of downwind routine
  ! class objects for subroutine arguments
  type(in_type_computJacob)       :: in_computJacob
  type(out_type_computJacob)      :: out_computJacob
  ! initialize error control
  err=0; message='testBandMat/'

  ! check
  if (in_SS4HG % nLeadDim == in_SS4HG % nState) then
   message=trim(message)//'do not expect nLeadDim==nState: check that are computing the band diagonal matrix'//&
                          ' (is forceFullMatrix==.true.?)'
   err=20; return
  end if

  ! compute the full Jacobian matrix
  call initialize_computJacob_testBandMat
  call computJacob(in_computJacob,indx_data,prog_data,diag_data,deriv_data,dBaseflow_dMatric,dMat,fullJac,out_computJacob)
  call finalize_computJacob_testBandMat(err,cmessage)
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

  ! initialize band matrix
  bandJac(:,:) = 0._rkind

  ! transfer into the lapack band diagonal structure
  do iState=1,in_SS4HG % nState
   do jState=max(1,iState-ku),min(in_SS4HG % nState,iState+kl)
    bandJac(kl + ku + 1 + jState - iState, iState) = fullJac(jState,iState)
   end do
  end do

  ! print results
  print*, '** test banded analytical Jacobian:'
  write(*,'(a4,1x,100(i17,1x))') 'xCol', (iState, iState=iJac1,iJac2)
  do iState=kl+1,in_SS4HG % nLeadDim; write(*,'(i4,1x,100(e17.10,1x))') iState, bandJac(iState,iJac1:iJac2); end do

  ! check if the need to pause
  if(check)then
   print*, 'PAUSE: testing banded analytical Jacobian'
   read(*,*)
  endif

  end subroutine testBandMat

  subroutine initialize_computJacob_testBandMat
   ! *** Transfer data from out_computJacob class object to local variables in testBandMat ***
   associate(&
    dt_cur         => in_SS4HG % dt_cur         ,& ! intent(in): current stepsize
    nSnow          => in_SS4HG % nSnow          ,& ! intent(in): number of snow layers
    nSoil          => in_SS4HG % nSoil          ,& ! intent(in): number of soil layers
    nLayers        => in_SS4HG % nLayers        ,& ! intent(in): total number of layers
    computeVegFlux => in_SS4HG % computeVegFlux  & ! intent(in): flag to indicate if computing fluxes over vegetation
   &)
    call in_computJacob % initialize(dt_cur,nSnow,nSoil,nLayers,computeVegFlux,.false.,ixFullMatrix)
   end associate
  end subroutine initialize_computJacob_testBandMat

  subroutine finalize_computJacob_testBandMat(err,cmessage)
   ! *** Transfer data from out_computJacob class object to local variables in testBandMat ***
   ! Note: subroutine arguments are needed because testBandMat is an internal procedure 
   integer(i4b),intent(out)        :: err                      ! error code
   character(*),intent(out)        :: cmessage                 ! error message of downwind routine
   call out_computJacob % finalize(err,cmessage)
  end subroutine finalize_computJacob_testBandMat
 
  subroutine initialize_computJacob_summaSolve4homegrown
   ! *** Transfer data to in_computJacob class object from local variables in summaSolve4homegrown ***
   associate(&
    ixGroundwater  => model_decisions(iLookDECISIONS%groundwatr)%iDecision,&  ! intent(in): [i4b] groundwater parameterization
    dt_cur         => in_SS4HG % dt_cur         ,& ! intent(in): current stepsize
    nSnow          => in_SS4HG % nSnow          ,& ! intent(in): number of snow layers
    nSoil          => in_SS4HG % nSoil          ,& ! intent(in): number of soil layers
    nLayers        => in_SS4HG % nLayers        ,& ! intent(in): total number of layers
    ixMatrix       => in_SS4HG % ixMatrix       ,& ! intent(in): type of matrix (full or band diagonal)
    computeVegFlux => in_SS4HG % computeVegFlux  & ! intent(in): flag to indicate if computing fluxes over vegetation
   &)   
    call in_computJacob % initialize(dt_cur,nSnow,nSoil,nLayers,computeVegFlux,(ixGroundwater==qbaseTopmodel),ixMatrix)
   end associate
  end subroutine initialize_computJacob_summaSolve4homegrown

  subroutine finalize_computJacob_summaSolve4homegrown
   ! *** Transfer data from out_computJacob class object to local variables in summaSolve4homegrown ***
   associate(err => out_SS4HG % err)
    call out_computJacob % finalize(err,cmessage)
   end associate 
  end subroutine finalize_computJacob_summaSolve4homegrown

  ! *********************************************************************************************************
  ! * internal subroutine eval8summa_wrapper: compute the right-hand-side vector
  ! *********************************************************************************************************
  ! NOTE: This is simply a wrapper routine for eval8summa, to reduce the number of calling arguments
  !       An internal subroutine, so have access to all data in the main subroutine
  subroutine eval8summa_wrapper(stateVecNew,fluxVecNew,resVecNew,fNew,feasible,err,message)
  USE eval8summa_module,only:eval8summa                      ! simulation of fluxes and residuals given a trial state vector
  implicit none
  ! input
  real(rkind),intent(in)         :: stateVecNew(:)           ! updated state vector
  ! output
  real(rkind),intent(out)        :: fluxVecNew(:)            ! updated flux vector
  real(qp),intent(out)           :: resVecNew(:) ! NOTE: qp  ! updated residual vector
  real(rkind),intent(out)        :: fNew                     ! new function value
  logical(lgt),intent(out)       :: feasible                 ! flag to denote the feasibility of the solution
  integer(i4b),intent(out)       :: err                      ! error code
  character(*),intent(out)       :: message                  ! error message
  ! ----------------------------------------------------------------------------------------------------------
  ! local
  character(len=256)             :: cmessage                 ! error message of downwind routine
  ! ----------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='eval8summa_wrapper/'

  associate(&
   dt_cur         => in_SS4HG % dt_cur         ,& ! intent(in): current stepsize
   dt             => in_SS4HG % dt             ,& ! intent(in): entire time step for drainage pond rate
   nSnow          => in_SS4HG % nSnow          ,& ! intent(in): number of snow layers
   nSoil          => in_SS4HG % nSoil          ,& ! intent(in): number of soil layers
   nLayers        => in_SS4HG % nLayers        ,& ! intent(in): total number of layers
   nState         => in_SS4HG % nState         ,& ! intent(in): total number of state variables
   firstSubStep   => in_SS4HG % firstSubStep   ,& ! intent(in): flag to indicate if we are processing the first sub-step
   computeVegFlux => in_SS4HG % computeVegFlux ,& ! intent(in): flag to indicate if computing fluxes over vegetation
   scalarSolution => in_SS4HG % scalarSolution ,& ! intent(in): flag to denote if implementing the scalar solution
   firstFluxCall  => io_SS4HG % firstFluxCall  ,& ! intent(inout): flag to indicate if we are processing the first flux call  
   ixSaturation   => io_SS4HG % ixSaturation    & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)    
  &)
  ! compute the flux and the residual vector for a given state vector
   call eval8summa(&
                   ! input: model control
                   dt_cur,                  & ! intent(in):    current stepsize
                   dt,                      & ! intent(in):    length of the time step (seconds)
                   nSnow,                   & ! intent(in):    number of snow layers
                   nSoil,                   & ! intent(in):    number of soil layers
                   nLayers,                 & ! intent(in):    total number of layers
                   nState,                  & ! intent(in):    total number of state variables
                   .false.,                 & ! intent(in):    not inside Sundials solver
                   firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                   firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                   .false.,                 & ! intent(in):    not processing the first iteration in a splitting operation
                   computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                   scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                   ! input: state vectors
                   stateVecNew,             & ! intent(in):    updated model state vector
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
  end associate

  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
  
  end subroutine eval8summa_wrapper


  ! *********************************************************************************************************
  ! internal function checkConv: check convergence based on the residual vector
  ! *********************************************************************************************************
  function checkConv(rVec,xInc,xVec)
  implicit none
  ! dummies
  real(rkind),intent(in)       :: rVec(:)                 ! residual vector (mixed units)
  real(rkind),intent(in)       :: xInc(:)                 ! iteration increment (mixed units)
  real(rkind),intent(in)       :: xVec(:)                 ! state vector (mixed units)
  logical(lgt)                 :: checkConv               ! flag to denote convergence
  ! locals
  real(rkind),dimension(mSoil) :: psiScale                ! scaling factor for matric head
  real(rkind),parameter        :: xSmall=1.e-0_rkind      ! a small offset
  real(rkind),parameter        :: scalarTighten=0.1_rkind ! scaling factor for the scalar solution
  real(rkind)                  :: soilWatbalErr           ! error in the soil water balance
  real(rkind)                  :: canopy_max              ! absolute value of the residual in canopy water (kg m-2)
  real(rkind),dimension(1)     :: energy_max              ! maximum absolute value of the energy residual (J m-3)
  real(rkind),dimension(1)     :: liquid_max              ! maximum absolute value of the volumetric liquid water content residual (-)
  real(rkind),dimension(1)     :: matric_max              ! maximum absolute value of the matric head iteration increment (m)
  real(rkind)                  :: aquifer_max             ! absolute value of the residual in aquifer water (m)
  logical(lgt)                 :: canopyConv              ! flag for canopy water balance convergence
  logical(lgt)                 :: watbalConv              ! flag for soil water balance convergence
  logical(lgt)                 :: liquidConv              ! flag for residual convergence
  logical(lgt)                 :: matricConv              ! flag for matric head convergence
  logical(lgt)                 :: energyConv              ! flag for energy convergence
  logical(lgt)                 :: aquiferConv             ! flag for aquifer water balance convergence
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  associate(&
   ! model control
   iter                    => in_SS4HG % iter                                   ,& ! intent(in): iteration index
   nsnow                   => in_SS4HG % nsnow                                  ,& ! intent(in): number of snow layers
   scalarSolution          => in_SS4HG % scalarSolution                         ,& ! intent(in): flag to denote if implementing the scalar solution
   ! convergence parameters
   absConvTol_liquid       => mpar_data%var(iLookPARAM%absConvTol_liquid)%dat(1),&  ! intent(in): [dp] absolute convergence tolerance for vol frac liq water (-)
   absConvTol_matric       => mpar_data%var(iLookPARAM%absConvTol_matric)%dat(1),&  ! intent(in): [dp] absolute convergence tolerance for matric head        (m)
   absConvTol_energy       => mpar_data%var(iLookPARAM%absConvTol_energy)%dat(1),&  ! intent(in): [dp] absolute convergence tolerance for energy             (J m-3)
   ! layer depth
   mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat          ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)
   ! model indices
   ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)          ,&  ! intent(in): [i4b]    index of aquifer storage state variable
   ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)         ,&  ! intent(in): [i4b]    index of canopy air space energy state variable
   ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)         ,&  ! intent(in): [i4b]    index of canopy energy state variable
   ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)         ,&  ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
   ixNrgOnly               => indx_data%var(iLookINDEX%ixNrgOnly)%dat           ,&  ! intent(in): [i4b(:)] list of indices for all energy states
   ixHydOnly               => indx_data%var(iLookINDEX%ixHydOnly)%dat           ,&  ! intent(in): [i4b(:)] list of indices for all hydrology states
   ixMatOnly               => indx_data%var(iLookINDEX%ixMatOnly)%dat           ,&  ! intent(in): [i4b(:)] list of indices for matric head state variables in the state vector
   ixMatricHead            => indx_data%var(iLookINDEX%ixMatricHead)%dat        ,&  ! intent(in): [i4b(:)] list of indices for matric head in the soil vector
   fnew                    => out_SS4HG % fnew                                   &  ! intent(in): [dp] new function evaluations
  ) ! making associations with variables in the data structures
  ! -------------------------------------------------------------------------------------------------------------------------------------------------

  ! check convergence based on the canopy water balance
  if(ixVegHyd/=integerMissing)then
   canopy_max = real(abs(rVec(ixVegHyd)), rkind)*iden_water
   canopyConv = (canopy_max    < absConvTol_liquid)  ! absolute error in canopy water balance (mm)
  else
   canopy_max = realMissing
   canopyConv = .true.
  endif

  ! check convergence based on the residuals for energy (J m-3)
  if(size(ixNrgOnly)>0)then
   energy_max = real(maxval(abs( rVec(ixNrgOnly) )), rkind)
   energyConv = (energy_max(1) < absConvTol_energy)  ! (based on the residual)
  else
   energy_max = realMissing
   energyConv = .true.
  endif

  ! check convergence based on the residuals for volumetric liquid water content (-)
  if(size(ixHydOnly)>0)then
   liquid_max = real(maxval(abs( rVec(ixHydOnly) ) ), rkind)
   ! (tighter convergence for the scalar solution)
   if(scalarSolution)then
    liquidConv = (liquid_max(1) < absConvTol_liquid*scalarTighten)   ! (based on the residual)
   else
    liquidConv = (liquid_max(1) < absConvTol_liquid)                 ! (based on the residual)
   endif
  else
   liquid_max = realMissing
   liquidConv = .true.
  endif

  ! check convergence based on the iteration increment for matric head
  ! NOTE: scale by matric head to avoid unnecessairly tight convergence when there is no water
  if(size(ixMatOnly)>0)then
   psiScale   = abs( xVec(ixMatOnly) ) + xSmall ! avoid divide by zero
   matric_max = maxval(abs( xInc(ixMatOnly)/psiScale ) )
   matricConv = (matric_max(1) < absConvTol_matric)  ! NOTE: based on iteration increment
  else
   matric_max = realMissing
   matricConv = .true.
  endif

  ! check convergence based on the soil water balance error (m)
  if(size(ixMatOnly)>0)then
   soilWatBalErr = sum( real(rVec(ixMatOnly), rkind)*mLayerDepth(nSnow+ixMatricHead) )
   watbalConv    = (abs(soilWatbalErr) < absConvTol_liquid)  ! absolute error in total soil water balance (m)
  else
   soilWatbalErr = realMissing
   watbalConv    = .true.
  endif

  ! check convergence based on the aquifer storage
  if(ixAqWat/=integerMissing)then
   aquifer_max = real(abs(rVec(ixAqWat)), rkind)*iden_water
   aquiferConv = (aquifer_max    < absConvTol_liquid)  ! absolute error in aquifer water balance (mm)
  else
   aquifer_max = realMissing
   aquiferConv = .true.
  endif

  ! final convergence check
  checkConv = (canopyConv .and. watbalConv .and. matricConv .and. liquidConv .and. energyConv .and. aquiferConv)

  ! print progress towards solution
  if(globalPrintFlag)then
   write(*,'(a,1x,i4,1x,7(e15.5,1x),7(L1,1x))') 'check convergence: ', iter, &
    fNew, matric_max(1), liquid_max(1), energy_max(1), canopy_max, aquifer_max, soilWatBalErr, matricConv, liquidConv, energyConv, watbalConv, canopyConv, aquiferConv, watbalConv
  endif

  ! end associations with variables in the data structures
  end associate

  end function checkConv

 end subroutine summaSolve4homegrown


end module summaSolve4homegrown_module
