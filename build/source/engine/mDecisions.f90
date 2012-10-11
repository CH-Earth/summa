module mDecisions_module
USE nrtype
implicit none
private
public::mDecisions
! -----------------------------------------------------------------------------------------------------------
! ***** define look-up values for different model decisions *****
! -----------------------------------------------------------------------------------------------------------
! look-up values for the choice of numerical method
integer(i4b),parameter,public :: iterative         =11    ! iterative
integer(i4b),parameter,public :: nonIterative      =12    ! non-iterative
integer(i4b),parameter,public :: iterSurfEnergyBal =13    ! iterate only on the surface energy balance
! look-up values for method used to compute derivative
integer(i4b),parameter,public :: numerical         =21    ! numerical solution
integer(i4b),parameter,public :: analytical        =22    ! analytical solution
! look-up values for the form of Richards' equation
integer(i4b),parameter,public :: moisture          =31    ! moisture-based form of Richards' equation
integer(i4b),parameter,public :: mixdform          =32    ! mixed form of Richards' equation
! look-up values for the choice of groundwater parameterization
integer(i4b),parameter,public :: movingBoundary    =41    ! moving lower boundary
integer(i4b),parameter,public :: bigBucket         =42    ! a big bucket (lumped aquifer model)
integer(i4b),parameter,public :: noExplicit        =43    ! no explicit groundwater parameterization
! look-up values for the choice of boundary conditions for thermodynamics
integer(i4b),parameter,public :: prescribedTemp    =51    ! prescribed temperature
integer(i4b),parameter,public :: energyFlux        =52    ! energy flux
integer(i4b),parameter,public :: zeroFlux          =53    ! zero flux
! look-up values for the choice of boundary conditions for hydrology
integer(i4b),parameter,public :: liquidFlux        =61    ! liquid water flux
integer(i4b),parameter,public :: prescribedHead    =62    ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
integer(i4b),parameter,public :: funcBottomHead    =63    ! function of matric head in the lower-most layer
integer(i4b),parameter,public :: freeDrainage      =64    ! free drainage
integer(i4b),parameter,public :: groundwaterCouple =65    ! coupled to the groundwater sub-model (matric head=0 as a moving lower boundary)
! look-up values for the choice of stability function
integer(i4b),parameter,public :: standard          =71    ! standard MO similarity, a la Anderson (1976) 
integer(i4b),parameter,public :: louisInversePower =72    ! Louis (1979) inverse power function
integer(i4b),parameter,public :: mahrtExponential  =73    ! Mahrt (1987) exponential
! look-up values for the choice of compaction routine
integer(i4b),parameter,public :: constantSettlement=81    ! constant settlement rate
integer(i4b),parameter,public :: andersonEmpirical =82    ! semi-empirical method of Anderson (1976)
! look-up values for the choice of thermal conductivity
integer(i4b),parameter,public :: Yen1965           =91    ! Yen (1965)
integer(i4b),parameter,public :: Mellor1977        =92    ! Mellor (1977)
integer(i4b),parameter,public :: Jordan1991        =93    ! Jordan (1991)
integer(i4b),parameter,public :: Smirnova2000      =94    ! Smirnova et al. (2000)
! look-up values for the choice of albedo representation
integer(i4b),parameter,public :: funcSnowAge       =101   ! function of snow age
integer(i4b),parameter,public :: BATSlike          =102   ! BATS-like approach, with destructive metamorphism + soot content
! -----------------------------------------------------------------------------------------------------------
contains

 ! ************************************************************************************************
 ! new subroutine: save model decisions as named integers
 ! ************************************************************************************************
 subroutine mDecisions(err,message)
 ! model decision structures
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 ! initialize error control
 err=0; message='mDecisions/'

 ! -------------------------------------------------------------------------------------------------
 ! -------------------------------------------------------------------------------------------------

 ! read information from model decisions file, and populate model decisions structure
 call readoption(err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! (1) identify the numerical method
 select case(trim(model_decisions(iLookDECISIONS%num_method)%cDecision))
  case('itertive'); model_decisions(iLookDECISIONS%num_method)%iDecision = iterative           ! iterative
  case('non_iter'); model_decisions(iLookDECISIONS%num_method)%iDecision = nonIterative        ! non-iterative
  case('itersurf'); model_decisions(iLookDECISIONS%num_method)%iDecision = iterSurfEnergyBal   ! iterate only on the surface energy balance
  case default
   err=10; message=trim(message)//"unknown numerical [option="//trim(model_decisions(iLookDECISIONS%num_method)%cDecision)//"]"; return
 end select

 ! (2) identify the method used to calculate flux derivatives
 select case(trim(model_decisions(iLookDECISIONS%fDerivMeth)%cDecision))
  case('numericl'); model_decisions(iLookDECISIONS%fDerivMeth)%iDecision = numerical           ! numerical
  case('analytic'); model_decisions(iLookDECISIONS%fDerivMeth)%iDecision = analytical          ! analytical
  case default
   err=10; message=trim(message)//"unknown method used to calculate flux derivatives [option="//trim(model_decisions(iLookDECISIONS%fDerivMeth)%cDecision)//"]"; return
 end select

 ! (3) identify the form of Richards' equation
 select case(trim(model_decisions(iLookDECISIONS%f_Richards)%cDecision))
  case('moisture'); model_decisions(iLookDECISIONS%f_Richards)%iDecision = moisture            ! moisture-based form
  case('mixdform'); model_decisions(iLookDECISIONS%f_Richards)%iDecision = mixdform            ! mixed form
  case default
   err=10; message=trim(message)//"unknown form of Richards' equation [option="//trim(model_decisions(iLookDECISIONS%f_Richards)%cDecision)//"]"; return
 end select

 ! (4) identify the groundwater parameterization
 select case(trim(model_decisions(iLookDECISIONS%groundwatr)%cDecision))
  case('movBound'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = movingBoundary      ! moving lower boundary
  case('bigBuckt'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = bigBucket           ! a big bucket (lumped aquifer model)
  case('noXplict'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = noExplicit          ! no explicit groundwater parameterization
  case default
   err=10; message=trim(message)//"unknown groundwater parameterization [option="//trim(model_decisions(iLookDECISIONS%groundwatr)%cDecision)//"]"; return
 end select

 ! (5) identify the upper boundary conditions for thermodynamics
 select case(trim(model_decisions(iLookDECISIONS%bcUpprTdyn)%cDecision))
  case('presTemp'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = prescribedTemp      ! prescribed temperature
  case('nrg_flux'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = energyFlux          ! energy flux
  case default
   err=10; message=trim(message)//"unknown upper boundary conditions for thermodynamics [option="//trim(model_decisions(iLookDECISIONS%bcUpprTdyn)%cDecision)//"]"; return
 end select

 ! (6) identify the lower boundary conditions for thermodynamics
 select case(trim(model_decisions(iLookDECISIONS%bcLowrTdyn)%cDecision))
  case('presTemp'); model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision = prescribedTemp      ! prescribed temperature
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown lower boundary conditions for thermodynamics [option="//trim(model_decisions(iLookDECISIONS%bcLowrTdyn)%cDecision)//"]"; return
 end select

 ! (7) identify the upper boundary conditions for soil hydrology
 select case(trim(model_decisions(iLookDECISIONS%bcUpprSoiH)%cDecision))
  case('presHead'); model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision = prescribedHead      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  case('liq_flux'); model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision = liquidFlux          ! liquid water flux
  case default
   err=10; message=trim(message)//"unknown upper boundary conditions for soil hydrology [option="//trim(model_decisions(iLookDECISIONS%bcUpprSoiH)%cDecision)//"]"; return
 end select

 ! (8) identify the lower boundary conditions for soil hydrology
 select case(trim(model_decisions(iLookDECISIONS%bcLowrSoiH)%cDecision))
  case('presHead'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = prescribedHead      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  case('bottmPsi'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = funcBottomHead      ! function of matric head in the lower-most layer
  case('drainage'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = freeDrainage        ! free drainage
  case('gwCouple'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = groundwaterCouple   ! coupled to the groundwater sub-model (matric head=0 as a moving lower boundary)
  case default
   err=10; message=trim(message)//"unknown lower boundary conditions for soil hydrology [option="//trim(model_decisions(iLookDECISIONS%bcLowrSoiH)%cDecision)//"]"; return
 end select

 ! (9) identify the choice of atmospheric stability function
 select case(trim(model_decisions(iLookDECISIONS%astability)%cDecision))
  case('standard'); model_decisions(iLookDECISIONS%astability)%iDecision = standard            ! standard MO similarity, a la Anderson (1976)
  case('louisinv'); model_decisions(iLookDECISIONS%astability)%iDecision = louisInversePower   ! Louis (1979) inverse power function
  case('mahrtexp'); model_decisions(iLookDECISIONS%astability)%iDecision = mahrtExponential    ! Mahrt (1987) exponential
  case default
   err=10; message=trim(message)//"unknown stability function [option="//trim(model_decisions(iLookDECISIONS%astability)%cDecision)//"]"; return
 end select

 ! (10) choice of compaction routine
 select case(trim(model_decisions(iLookDECISIONS%compaction)%cDecision))
  case('consettl'); model_decisions(iLookDECISIONS%compaction)%iDecision = constantSettlement  ! constant settlement rate
  case('anderson'); model_decisions(iLookDECISIONS%compaction)%iDecision = andersonEmpirical   ! semi-empirical method of Anderson (1976)
  case default
   err=10; message=trim(message)//"unknown option for snow compaction [option="//trim(model_decisions(iLookDECISIONS%compaction)%cDecision)//"]"; return
 end select

 ! (11) choice of thermal conductivity
 select case(trim(model_decisions(iLookDECISIONS%thermlcond)%cDecision))
  case('tyen1965'); model_decisions(iLookDECISIONS%thermlcond)%iDecision = Yen1965             ! Yen (1965) 
  case('melr1977'); model_decisions(iLookDECISIONS%thermlcond)%iDecision = Mellor1977          ! Mellor (1977)
  case('jrdn1991'); model_decisions(iLookDECISIONS%thermlcond)%iDecision = Jordan1991          ! Jordan (1991)
  case('smnv2000'); model_decisions(iLookDECISIONS%thermlcond)%iDecision = Smirnova2000        ! Smirnova et al. (2000)
  case default
   err=10; message=trim(message)//"unknown option for thermal conductivity [option="//trim(model_decisions(iLookDECISIONS%thermlcond)%cDecision)//"]"; return
 end select

 ! (12) choice of albedo representation
 select case(trim(model_decisions(iLookDECISIONS%alb_method)%cDecision))
  case('fsnowage'); model_decisions(iLookDECISIONS%alb_method)%iDecision = funcSnowAge         ! function of snow age
  case('BATSlike'); model_decisions(iLookDECISIONS%alb_method)%iDecision = BATSlike            ! BATS-like approach, with destructive metamorphism + soot content
  case default
   err=10; message=trim(message)//"unknown option for snow albedo [option="//trim(model_decisions(iLookDECISIONS%alb_method)%cDecision)//"]"; return
 end select

 ! *****
 ! check for unrealistic options

 ! only allow the moving boundary gw parameterization with the moisture-based form of Richards' equation
 if(model_decisions(iLookDECISIONS%groundwatr)%iDecision == movingBoundary)then
  if(model_decisions(iLookDECISIONS%f_Richards)%iDecision /= moisture)then
   err=20; message=trim(message)//"moving boundary gw parameterization only allowed with the moisture-based from of Richards' equation"; return
  endif
 endif

 ! moving boundary gw parameterization must have "groundwaterCouple" as the lower boundary condition
 if(model_decisions(iLookDECISIONS%groundwatr)%iDecision == movingBoundary)then
  if(model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision /= groundwaterCouple)then
   err=20; message=trim(message)//'moving boundary gw parameterization must use "gwCouple" as the lower model boundary condition'; return
  endif
 endif

 end subroutine mDecisions

 ! ************************************************************************************************
 ! private subroutine: read information from model decisions file
 ! ************************************************************************************************
 subroutine readoption(err,message)
 ! used to read information from model decisions file
 USE ascii_util_module,only:file_open       ! open file
 USE ascii_util_module,only:get_vlines      ! get a vector of non-comment lines
 USE snow_fileManager,only:SETNGS_PATH      ! path for metadata files
 USE snow_fileManager,only:M_DECISIONS      ! definition of modeling options
 USE get_ixname_module,only:get_ixdecisions ! identify index of named variable
 USE data_struc,only:model_decisions        ! model decision structure
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b),parameter               :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 character(LEN=256),allocatable       :: charline(:)    ! vector of character strings
 integer(i4b)                         :: nDecisions     ! number of model decisions
 integer(i4b)                         :: iDecision      ! index of model decisions
 character(len=32)                    :: decision       ! name of model decision
 character(len=32)                    :: option         ! option for model decision
 integer(i4b)                         :: iVar           ! index of the decision in the data structure
 ! Start procedure here
 err=0; message='readoption/'
 ! build filename
 infile = trim(SETNGS_PATH)//trim(M_DECISIONS)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get a list of character strings from non-comment lines
 call get_vlines(unt,charline,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! close the file unit
 close(unt)
 ! get the number of model decisions
 nDecisions = size(charline)
 ! allocate space for the model decisions
 if(associated(model_decisions)) deallocate(model_decisions)
 allocate(model_decisions(nDecisions),stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemAllocateModelDecisions"; return; endif
 ! populate the model decisions structure
 do iDecision=1,nDecisions
  ! extract name of decision and the decision selected
  read(charline(iDecision),*,iostat=err) option, decision
  if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
  ! get the index of the decision in the data structure
  iVar = get_ixdecisions(trim(option))
  if(iVar<=0)then; err=40; message=trim(message)//"cannotFindDecisionIndex[name='"//trim(option)//"']"; return; endif
  ! populate the model decisions structure
  model_decisions(iVar)%cOption   = trim(option) 
  model_decisions(iVar)%cDecision = trim(decision)
 end do
 end subroutine readoption

end module mDecisions_module
