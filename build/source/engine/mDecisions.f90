module mDecisions_module
USE nrtype
implicit none
private
public::mDecisions
! -----------------------------------------------------------------------------------------------------------
! ***** define look-up values for different Noah-MP decisions *****
! -----------------------------------------------------------------------------------------------------------
! look-up values for the choice of function for the soil moisture control on stomatal resistance
integer(i4b),parameter,public :: NoahType          = 1    ! thresholded linear function of volumetric liquid water content
integer(i4b),parameter,public :: CLM_type          = 2    ! thresholded linear function of matric head
integer(i4b),parameter,public :: SiB_Type          = 3    ! exponential of the log of matric head
! look-up values for the choice of function for the soil moisture control on stomatal resistance
integer(i4b),parameter,public :: BallBerry         = 1    ! Ball-Berry
integer(i4b),parameter,public :: Jarvis            = 2    ! Jarvis
! -----------------------------------------------------------------------------------------------------------
! ***** define look-up values for different FUSE model decisions *****
! -----------------------------------------------------------------------------------------------------------
! look-up values for the choice of numerical method
integer(i4b),parameter,public :: iterative            =  11    ! iterative
integer(i4b),parameter,public :: nonIterative         =  12    ! non-iterative
integer(i4b),parameter,public :: iterSurfEnergyBal    =  13    ! iterate only on the surface energy balance
! look-up values for method used to compute derivative
integer(i4b),parameter,public :: numerical            =  21    ! numerical solution
integer(i4b),parameter,public :: analytical           =  22    ! analytical solution
! look-up values for the form of Richards' equation
integer(i4b),parameter,public :: moisture             =  31    ! moisture-based form of Richards' equation
integer(i4b),parameter,public :: mixdform             =  32    ! mixed form of Richards' equation
! look-up values for the choice of groundwater parameterization
integer(i4b),parameter,public :: equilWaterTable      =  41    ! equilibrium water table
integer(i4b),parameter,public :: pseudoWaterTable     =  42    ! pseudo water table
integer(i4b),parameter,public :: bigBucket            =  43    ! a big bucket (lumped aquifer model)
integer(i4b),parameter,public :: noExplicit           =  44    ! no explicit groundwater parameterization
! look-up values for the choice of hydraulic conductivity profile
integer(i4b),parameter,public :: constant             =  51    ! constant hydraulic conductivity with depth
integer(i4b),parameter,public :: exp_profile          =  52    ! exponential profile
integer(i4b),parameter,public :: powerLaw_profile     =  53    ! power-law profile
integer(i4b),parameter,public :: linear_profile       =  54    ! linear profile
! look-up values for the choice of boundary conditions for thermodynamics
integer(i4b),parameter,public :: prescribedTemp       =  61    ! prescribed temperature
integer(i4b),parameter,public :: energyFlux           =  62    ! energy flux
integer(i4b),parameter,public :: zeroFlux             =  63    ! zero flux
! look-up values for the choice of boundary conditions for hydrology
integer(i4b),parameter,public :: liquidFlux           =  71    ! liquid water flux
integer(i4b),parameter,public :: prescribedHead       =  72    ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
integer(i4b),parameter,public :: funcBottomHead       =  73    ! function of matric head in the lower-most layer
integer(i4b),parameter,public :: freeDrainage         =  74    ! free drainage
! look-up values for the choice of stability function
integer(i4b),parameter,public :: standard             =  81    ! standard MO similarity, a la Anderson (1976) 
integer(i4b),parameter,public :: louisInversePower    =  82    ! Louis (1979) inverse power function
integer(i4b),parameter,public :: mahrtExponential     =  83    ! Mahrt (1987) exponential
! look-up values for the choice of albedo representation
integer(i4b),parameter,public :: funcSnowAge          =  91    ! function of snow age
integer(i4b),parameter,public :: BATSlike             =  92    ! BATS-like approach, with destructive metamorphism + soot content
! look-up values for the choice of compaction routine
integer(i4b),parameter,public :: constantSettlement   = 101    ! constant settlement rate
integer(i4b),parameter,public :: andersonEmpirical    = 102    ! semi-empirical method of Anderson (1976)
! look-up values for the choice of method to combine and sub-divide snow layers
integer(i4b),parameter,public :: sameRulesAllLayers   = 111    ! same combination/sub-division rules applied to all layers
integer(i4b),parameter,public :: rulesDependLayerIndex= 112    ! combination/sub-dividion rules depend on layer index
! look-up values for the choice of thermal conductivity
integer(i4b),parameter,public :: Yen1965              = 121    ! Yen (1965)
integer(i4b),parameter,public :: Mellor1977           = 122    ! Mellor (1977)
integer(i4b),parameter,public :: Jordan1991           = 123    ! Jordan (1991)
integer(i4b),parameter,public :: Smirnova2000         = 124    ! Smirnova et al. (2000)
! look-up values for the choice of routing method
integer(i4b),parameter,public :: timeDelay            = 131    ! time-delay histogram
integer(i4b),parameter,public :: qInstant             = 132    ! instantaneous routing
! -----------------------------------------------------------------------------------------------------------
contains

 ! ************************************************************************************************
 ! new subroutine: save model decisions as named integers
 ! ************************************************************************************************
 subroutine mDecisions(err,message)
 ! model decision structures
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! Noah-MP decision structures
 USE noahmp_globals,only:DVEG               ! decision for dynamic vegetation
 USE noahmp_globals,only:OPT_RAD            ! decision for canopy radiation
 USE noahmp_globals,only:OPT_ALB            ! decision for snow albedo
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

 ! -------------------------------------------------------------------------------------------------

 ! (0) set Noah-MP options
 DVEG=3      ! option for dynamic vegetation
 OPT_RAD=1   ! option for canopy radiation
 OPT_ALB=1   ! option for snow albedo

 ! (N-03) identify the choice of function for the soil moisture control on stomatal resistance
 select case(trim(model_decisions(iLookDECISIONS%soilStress)%cDecision))
  case('NoahType'); model_decisions(iLookDECISIONS%soilStress)%iDecision = NoahType             ! thresholded linear function of volumetric liquid water content
  case('CLM_type'); model_decisions(iLookDECISIONS%soilStress)%iDecision = CLM_type             ! thresholded linear function of matric head
  case('SiB_Type'); model_decisions(iLookDECISIONS%soilStress)%iDecision = SiB_Type             ! exponential of the log of matric head
  case default
   err=10; message=trim(message)//"unknown numerical [option="//trim(model_decisions(iLookDECISIONS%soilStress)%cDecision)//"]"; return
 end select

 ! (N-04) identify the choice of function for stomatal resistance
 select case(trim(model_decisions(iLookDECISIONS%stomResist)%cDecision))
  case('BallBerry'); model_decisions(iLookDECISIONS%stomResist)%iDecision = BallBerry           ! Ball-Berry
  case('Jarvis'   ); model_decisions(iLookDECISIONS%stomResist)%iDecision = Jarvis              ! Jarvis
  case default
   err=10; message=trim(message)//"unknown numerical [option="//trim(model_decisions(iLookDECISIONS%stomResist)%cDecision)//"]"; return
 end select

 ! -------------------------------------------------------------------------------------------------

 ! (F-01) identify the numerical method
 select case(trim(model_decisions(iLookDECISIONS%num_method)%cDecision))
  case('itertive'); model_decisions(iLookDECISIONS%num_method)%iDecision = iterative           ! iterative
  case('non_iter'); model_decisions(iLookDECISIONS%num_method)%iDecision = nonIterative        ! non-iterative
  case('itersurf'); model_decisions(iLookDECISIONS%num_method)%iDecision = iterSurfEnergyBal   ! iterate only on the surface energy balance
  case default
   err=10; message=trim(message)//"unknown numerical [option="//trim(model_decisions(iLookDECISIONS%num_method)%cDecision)//"]"; return
 end select

 ! (F-02) identify the method used to calculate flux derivatives
 select case(trim(model_decisions(iLookDECISIONS%fDerivMeth)%cDecision))
  case('numericl'); model_decisions(iLookDECISIONS%fDerivMeth)%iDecision = numerical           ! numerical
  case('analytic'); model_decisions(iLookDECISIONS%fDerivMeth)%iDecision = analytical          ! analytical
  case default
   err=10; message=trim(message)//"unknown method used to calculate flux derivatives [option="//trim(model_decisions(iLookDECISIONS%fDerivMeth)%cDecision)//"]"; return
 end select

 ! (F-03) identify the form of Richards' equation
 select case(trim(model_decisions(iLookDECISIONS%f_Richards)%cDecision))
  case('moisture'); model_decisions(iLookDECISIONS%f_Richards)%iDecision = moisture            ! moisture-based form
  case('mixdform'); model_decisions(iLookDECISIONS%f_Richards)%iDecision = mixdform            ! mixed form
  case default
   err=10; message=trim(message)//"unknown form of Richards' equation [option="//trim(model_decisions(iLookDECISIONS%f_Richards)%cDecision)//"]"; return
 end select

 ! (F-04) identify the groundwater parameterization
 select case(trim(model_decisions(iLookDECISIONS%groundwatr)%cDecision))
  case('zEquilWT'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = equilWaterTable     ! equilibrium water table
  case('pseudoWT'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = pseudoWaterTable    ! pseudo water table
  case('bigBuckt'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = bigBucket           ! a big bucket (lumped aquifer model)
  case('noXplict'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = noExplicit          ! no explicit groundwater parameterization
  case default
   err=10; message=trim(message)//"unknown groundwater parameterization [option="//trim(model_decisions(iLookDECISIONS%groundwatr)%cDecision)//"]"; return
 end select

 ! (F-05) identify the hydraulic conductivity profile
 select case(trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision))
  case('constant'); model_decisions(iLookDECISIONS%hc_profile)%iDecision = constant            ! constant hydraulic conductivity with depth
  case('exp_prof'); model_decisions(iLookDECISIONS%hc_profile)%iDecision = exp_profile         ! exponential profile
  case('pow_prof'); model_decisions(iLookDECISIONS%hc_profile)%iDecision = powerLaw_profile    ! power-law profile
  case('lin_prof'); model_decisions(iLookDECISIONS%hc_profile)%iDecision = linear_profile      ! linear profile
  case default
   err=10; message=trim(message)//"unknown hydraulic conductivity profile [option="//trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)//"]"; return
 end select

 ! (F-06) identify the upper boundary conditions for thermodynamics
 select case(trim(model_decisions(iLookDECISIONS%bcUpprTdyn)%cDecision))
  case('presTemp'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = prescribedTemp      ! prescribed temperature
  case('nrg_flux'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = energyFlux          ! energy flux
  case default
   err=10; message=trim(message)//"unknown upper boundary conditions for thermodynamics [option="//trim(model_decisions(iLookDECISIONS%bcUpprTdyn)%cDecision)//"]"; return
 end select

 ! (F-07) identify the lower boundary conditions for thermodynamics
 select case(trim(model_decisions(iLookDECISIONS%bcLowrTdyn)%cDecision))
  case('presTemp'); model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision = prescribedTemp      ! prescribed temperature
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown lower boundary conditions for thermodynamics [option="//trim(model_decisions(iLookDECISIONS%bcLowrTdyn)%cDecision)//"]"; return
 end select

 ! (F-08) identify the upper boundary conditions for soil hydrology
 select case(trim(model_decisions(iLookDECISIONS%bcUpprSoiH)%cDecision))
  case('presHead'); model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision = prescribedHead      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  case('liq_flux'); model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision = liquidFlux          ! liquid water flux
  case default
   err=10; message=trim(message)//"unknown upper boundary conditions for soil hydrology [option="//trim(model_decisions(iLookDECISIONS%bcUpprSoiH)%cDecision)//"]"; return
 end select

 ! (F-09) identify the lower boundary conditions for soil hydrology
 select case(trim(model_decisions(iLookDECISIONS%bcLowrSoiH)%cDecision))
  case('presHead'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = prescribedHead      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  case('bottmPsi'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = funcBottomHead      ! function of matric head in the lower-most layer
  case('drainage'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = freeDrainage        ! free drainage
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown lower boundary conditions for soil hydrology [option="//trim(model_decisions(iLookDECISIONS%bcLowrSoiH)%cDecision)//"]"; return
 end select

 ! (F-10) identify the choice of atmospheric stability function
 select case(trim(model_decisions(iLookDECISIONS%astability)%cDecision))
  case('standard'); model_decisions(iLookDECISIONS%astability)%iDecision = standard            ! standard MO similarity, a la Anderson (1976)
  case('louisinv'); model_decisions(iLookDECISIONS%astability)%iDecision = louisInversePower   ! Louis (1979) inverse power function
  case('mahrtexp'); model_decisions(iLookDECISIONS%astability)%iDecision = mahrtExponential    ! Mahrt (1987) exponential
  case default
   err=10; message=trim(message)//"unknown stability function [option="//trim(model_decisions(iLookDECISIONS%astability)%cDecision)//"]"; return
 end select

 ! (F-11) choice of albedo representation
 select case(trim(model_decisions(iLookDECISIONS%alb_method)%cDecision))
  case('fsnowage'); model_decisions(iLookDECISIONS%alb_method)%iDecision = funcSnowAge         ! function of snow age
  case('BATSlike'); model_decisions(iLookDECISIONS%alb_method)%iDecision = BATSlike            ! BATS-like approach, with destructive metamorphism + soot content
  case default
   err=10; message=trim(message)//"unknown option for snow albedo [option="//trim(model_decisions(iLookDECISIONS%alb_method)%cDecision)//"]"; return
 end select

 ! (F-12) choice of snow compaction routine
 select case(trim(model_decisions(iLookDECISIONS%compaction)%cDecision))
  case('consettl'); model_decisions(iLookDECISIONS%compaction)%iDecision = constantSettlement  ! constant settlement rate
  case('anderson'); model_decisions(iLookDECISIONS%compaction)%iDecision = andersonEmpirical   ! semi-empirical method of Anderson (1976)
  case default
   err=10; message=trim(message)//"unknown option for snow compaction [option="//trim(model_decisions(iLookDECISIONS%compaction)%cDecision)//"]"; return
 end select

 ! (F-13) choice of method to combine and sub-divide snow layers
 select case(trim(model_decisions(iLookDECISIONS%snowLayers)%cDecision))
  case('jrdn1991'); model_decisions(iLookDECISIONS%snowLayers)%iDecision = sameRulesAllLayers    ! SNTHERM option: same combination/sub-dividion rules applied to all layers
  case('CLM_2010'); model_decisions(iLookDECISIONS%snowLayers)%iDecision = rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
  case default
   err=10; message=trim(message)//"unknown option for combination/sub-division of snow layers [option="//trim(model_decisions(iLookDECISIONS%snowLayers)%cDecision)//"]"; return
 end select

 ! (F-14) choice of thermal conductivity
 select case(trim(model_decisions(iLookDECISIONS%thermlcond)%cDecision))
  case('tyen1965'); model_decisions(iLookDECISIONS%thermlcond)%iDecision = Yen1965             ! Yen (1965) 
  case('melr1977'); model_decisions(iLookDECISIONS%thermlcond)%iDecision = Mellor1977          ! Mellor (1977)
  case('jrdn1991'); model_decisions(iLookDECISIONS%thermlcond)%iDecision = Jordan1991          ! Jordan (1991)
  case('smnv2000'); model_decisions(iLookDECISIONS%thermlcond)%iDecision = Smirnova2000        ! Smirnova et al. (2000)
  case default
   err=10; message=trim(message)//"unknown option for thermal conductivity [option="//trim(model_decisions(iLookDECISIONS%thermlcond)%cDecision)//"]"; return
 end select

 ! (F-15) choice of routing method
 select case(trim(model_decisions(iLookDECISIONS%subRouting)%cDecision))
  case('timeDlay'); model_decisions(iLookDECISIONS%subRouting)%iDecision = timeDelay           ! time-delay histogram
  case('qInstant'); model_decisions(iLookDECISIONS%subRouting)%iDecision = qInstant            ! instantaneous routing
  case default
   err=10; message=trim(message)//"unknown option for sub-grid routing [option="//trim(model_decisions(iLookDECISIONS%subRouting)%cDecision)//"]"; return
 end select

 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 ! check for consistency among options
 select case(model_decisions(iLookDECISIONS%groundwatr)%iDecision)
  case(equilWaterTable,pseudoWaterTable)
   if(model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision /= zeroFlux)then
    message=trim(message)//'lower boundary condition for soil hydology must be zeroFlux with (zEquilWT or pseudoWT) options for groundwater'
    err=20; return 
   endif
 end select

 print*, 'check decisions', model_decisions(iLookDECISIONS%snowLayers)%iDecision, sameRulesAllLayers, rulesDependLayerIndex
 pause

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
