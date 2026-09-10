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

module read_pinit_module
USE nr_type
! check for when model decisions are undefined
USE mDecisions_module,only: unDefined
USE globalData,only:model_decisions
USE globalData,only:realMissing
USE multiconst,only:secprhour  ! number of seconds in an hour
USE var_lookup,only:iLookDECISIONS,iLookPARAM,iLookBPAR
implicit none
private
public::read_pinit
contains


 ! ************************************************************************************************
 ! public subroutine read_pinit: read default model parameter values and constraints
 ! ************************************************************************************************
 subroutine read_pinit(filenm,isLocal,absEnergyFac,mpar_meta,parFallback,err,message)
 ! used to read metadata on the forcing data file
 USE summaFileManager,only:SETTINGS_PATH   ! path for input parameter and other configuration files
 USE ascii_util_module,only:file_open      ! open ascii file
 USE ascii_util_module,only:split_line     ! extract the list of variable names from the character string
 USE data_types,only:var_info              ! data type for metadata
 USE data_types,only:par_info              ! data type for parameter constraints
 USE get_ixname_module,only:get_ixParam    ! identify index of named variable for local column model parameters
 USE get_ixname_module,only:get_ixBpar     ! identify index of named variable for basin-average model parameters
 implicit none
 ! define input
 character(*),intent(in)                :: filenm         ! name of file containing default values and constraints of model parameters
 logical(lgt),intent(in)                :: isLocal        ! .true. if the file describes local column parameters
 real(rkind),intent(in)                 :: absEnergyFac   ! multiplier for absolute value of energy state variable (for enthalpy or temperature)
 type(var_info),intent(in)              :: mpar_meta(:)   ! metadata for model parameters
 ! define output
 type(par_info),intent(out)             :: parFallback(:) ! default values and constraints of model parameters
 integer(i4b),intent(out)               :: err            ! error code
 character(*),intent(out)               :: message        ! error message
 ! define general variables
 logical(lgt),parameter                 :: backwardsCompatible=.false. ! .true. if skip check that all parameters are populated
 character(len=256)                     :: cmessage       ! error message for downwind routine
 character(LEN=256)                     :: infile         ! input filename
 integer(i4b)                           :: unt            ! file unit (free unit output from file_open)
 integer(i4b)                           :: iline          ! loop through lines in the file
 integer(i4b),parameter                 :: maxLines=1000  ! maximum lines in the file
 character(LEN=256)                     :: temp           ! single line of information
 ! define local variables for the default model parameters
 integer(i4b)                           :: iend           ! check for the end of the file
 character(LEN=256)                     :: ffmt           ! file format
 character(LEN=32)                      :: varName        ! name of variable
 type(par_info)                         :: parTemp        ! temporary parameter structure
 character(LEN=2)                       :: dLim           ! column delimiter
 integer(i4b)                           :: iVar           ! index of model variable
 ! Start procedure here
 err=0; message="read_pinit/"
 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename and update error message
 infile = trim(SETTINGS_PATH)//trim(filenm)
 message=trim(message)//'file='//trim(infile)//' - '
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! **********************************************************************************************
 ! (2) read default model parameter values and constraints
 ! **********************************************************************************************
 ! fill parameter vector with missing data
 parFallback(:)%default_val = realMissing
 parFallback(:)%lower_limit = realMissing
 parFallback(:)%upper_limit = realMissing
 ! ---------------------------------------------------------------------------------------------
 ! read format code
 ! ---------------------------------------------------------------------------------------------
 do iline=1,maxLines
  ! (read through comment lines)
  read(unt,'(a)',iostat=iend) temp  ! read a line of data
  if(iend/=0)then; err=20; message=trim(message)//'got to end of file before found the format code'; return; end if
  if (temp(1:1)=='!')cycle
  ! (read in format string -- assume that the first non-comment line is the format code)
  read(temp,*)ffmt  ! read in format string
  exit
  if(iLine==maxLines)then; err=20; message=trim(message)//'problem finding format code -- no non-comment line after start of parameter definitions'; return; end if
 end do ! looping through lines
 ! ---------------------------------------------------------------------------------------------
 ! read in default values of model parameters, and parameter constraints
 ! ---------------------------------------------------------------------------------------------
 do iline=1,maxLines
  ! (read through comment lines)
  read(unt,'(a)',iostat=iend) temp  ! read a line of data
  if(iend/=0)exit !end of file
  if (temp(1:1)=='!')cycle
  ! (save data into a temporary variables)
  read(temp,trim(ffmt),iostat=err) varName, dLim, parTemp%default_val, dLim, parTemp%lower_limit, dLim, parTemp%upper_limit
  if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; end if
  ! skip deprecated parameters (for backwards compatibility with existing parameter input files)
  if (trim(varName) == 'upperBoundTheta' .or. trim(varName) == 'lowerBoundTheta')then
    write(*,'(a)') "WARNING: deprecated parameter '"//trim(varName)//"' found in parameter input file -- ignoring this parameter"
    cycle
  end if
  ! (identify the index of the variable in the data structure)
  if(isLocal)then
   iVar = get_ixParam(trim(varName))
  else
   iVar = get_ixBpar(trim(varName))
  end if
  ! (check that we have successfully found the parameter)
  if(iVar>0)then
   if(iVar>size(parFallback))then
    err=35; message=trim(message)//"indexOutOfRange[var="//trim(varName)//"]"; return
   end if
   ! (put data in the structure)
   parFallback(iVar)=parTemp
  else
   err=40; message=trim(message)//"variable in parameter file not present in data structure [var="//trim(varName)//"]"; return
  end if
 end do  ! (looping through lines in the file)

 ! add these defaults for backwards compatibility pre Sundials, FUSE, and glacier and lake domains
 if (isLocal) then ! dealing with parameters for local column
  ! BE solver parameters
  if (parFallback(iLookPARAM%be_steps)%default_val < 0.99_rkind*realMissing) then
   parFallback(iLookPARAM%be_steps)%default_val = 1._rkind
  end if
  ! IDA solver parameters
  call set_ida_defaults(absEnergyFac, parFallback, err, cmessage)
  if (err /= 0) then; message = trim(message)//trim(cmessage); return; end if
  ! set FUSE parameter defaults
  call set_FUSE_defaults(parFallback, err, cmessage)
  if (err /= 0) then; message = trim(message)//trim(cmessage); return; end if
  ! glacier and lake parameters
  if (parFallback(iLookPARAM%albedoFrznWatVisible)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%albedoFrznWatVisible)%default_val = 0.6_rkind
  end if
  if (parFallback(iLookPARAM%albedoFrznWatNearIR)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%albedoFrznWatNearIR)%default_val = 0.4_rkind
  end if
  if (parFallback(iLookPARAM%albedoOpenWatVisible)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%albedoOpenWatVisible)%default_val = 0.06_rkind
  end if
  if (parFallback(iLookPARAM%albedoOpenWatNearIR)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%albedoOpenWatNearIR)%default_val = 0.06_rkind
  end if
  if (parFallback(iLookPARAM%z0Water)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%z0Water)%default_val = 0.0005_rkind
  end if
  if (parFallback(iLookPARAM%z0Ice)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%z0Ice)%default_val = 0.0010_rkind
  end if
  if (parFallback(iLookPARAM%glacierWindFactor)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%glacierWindFactor)%default_val = 1._rkind ! 
  end if
  if (parFallback(iLookPARAM%glacierTempReduction)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%glacierTempReduction)%default_val = 0._rkind
  end if
 else
  ! glacier parameters
  if (parFallback(iLookBPAR%glacStor_kIce)%default_val < 0.99_rkind*realMissing) then ! 5-29
   parFallback(iLookBPAR%glacStor_kIce)%default_val = 15._rkind*secprhour! convert from hours to seconds
  end if
  if (parFallback(iLookBPAR%glacStor_kSnow)%default_val < 0.99_rkind*realMissing) then ! 30-149
   parFallback(iLookBPAR%glacStor_kSnow)%default_val = 90._rkind*secprhour ! convert from hours to seconds
  end if
  if (parFallback(iLookBPAR%glacStor_kFirn)%default_val < 0.99_rkind*realMissing) then ! 150-1000
    parFallback(iLookBPAR%glacStor_kFirn)%default_val = 575._rkind*secprhour ! convert from hours to seconds
  endif
  if (parFallback(iLookBPAR%debrisConc)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookBPAR%debrisConc)%default_val = 5.0_rkind ! 0.1 to 6.4 kg/m3 following Anderson and Anderson (2018)
  endif
  if (parFallback(iLookBPAR%wallErosionRate)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookBPAR%wallErosionRate)%default_val = 8.0_rkind ! 1 to 15 mm yr-1 following Anderson and Anderson (2016)
  endif
  if (parFallback(iLookBPAR%debrisCritStress)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookBPAR%debrisCritStress)%default_val = 80000 ! 20000-100000 Pa follow Mayer and Licciulli (2021)
  endif
  if (parFallback(iLookBPAR%latMoraineWidth)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookBPAR%latMoraineWidth)%default_val = 200._rkind ! from looking at Alaska glaciers (m)
  endif
  if (parFallback(iLookPARAM%f_hydCond)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%f_hydCond)%default_val = 3._rkind ! 1-5 m-1 for supraglacial debris and weathered shallow till
  end if
 end if

 ! check we have populated all variables
 ! NOTE: ultimately need a need a parameter dictionary to ensure that the parameters used are populated
 if(.not.backwardsCompatible)then  ! if we add new variables in future versions of the code, then some may be missing in the input file
  if(any(parFallback(:)%default_val < 0.99_rkind*realMissing))then
   do iVar=1,size(parFallback)
    if(parFallback(iVar)%default_val < 0.99_rkind*realMissing)then
     err=40; message=trim(message)//"variableNonexistent[var="//trim(mpar_meta(iVar)%varName)//"]"; return
    end if
   end do
  end if
 ! populate parameters that were not included in the original control files
 else ! (need backwards compatibility)
  if(isLocal)then
   if(model_decisions(iLookDECISIONS%cIntercept)%iDecision == unDefined)then
    parFallback(iLookPARAM%canopyWettingFactor)%default_val = 1._rkind             ! maximum wetted fraction of the canopy (-)
    parFallback(iLookPARAM%canopyWettingExp)%default_val    = 0.666666667_rkind    ! exponent in canopy wetting function (-)
   end if
  end if
 end if
 
 ! close file unit
 close(unt)
 end subroutine read_pinit

 ! ************************************************************************************************
 ! Subroutine to separate the default settings of the IDA solver from the rest of the model parameters
 ! ************************************************************************************************
 subroutine set_ida_defaults(absEnergyFac, parFallback, err, message)
 USE data_types,only:par_info              ! data type for parameter constraints
 implicit none
 ! define input
 real(rkind),intent(in)                 :: absEnergyFac   ! multiplier for absolute value of energy state variable (for enthalpy or temperature)
 type(par_info),intent(out)             :: parFallback(:) ! default values and constraints of model parameters
 integer(i4b),intent(out)               :: err            ! error code
 character(*),intent(out)               :: message        ! error message
 ! local varaibles
 integer(i4b)                           :: i   
 real(rkind)                            :: default_relTol = 1.e-5_rkind 
 real(rkind)                            :: default_absTol = 1.e-5_rkind
 integer(i4b), dimension(7)             :: relTol_paramIndx = [iLookPARAM%relTolTempCas, iLookPARAM%relTolTempVeg, iLookPARAM%relTolWatVeg, &
                                                               iLookPARAM%relTolTempSoilSnow, iLookPARAM%relTolWatSnow, iLookPARAM%relTolMatric, &
                                                               iLookPARAM%relTolAquifr]
 integer(i4b), dimension(3)             :: absTolTemp_paramIndx = [iLookPARAM%absTolTempCas, iLookPARAM%absTolTempVeg, iLookPARAM%absTolTempSoilSnow]
 integer(i4b), dimension(4)             :: absTolWat_paramIndx =  [iLookPARAM%absTolWatVeg, iLookPARAM%absTolWatSnow, iLookPARAM%absTolMatric, &
                                                                   iLookPARAM%absTolAquifr]
 err=0 ! initialize error code
 message="set_ida_defaults/"
 
  ! Relative Tolerances
  do i = 1, size(relTol_paramIndx)
    if (parFallback(relTol_paramIndx(i))%default_val < 0.99_rkind*realMissing) then
      parFallback(relTol_paramIndx(i))%default_val = default_relTol
    end if
  end do

  ! Absolute Tolerances
  do i = 1, size(absTolTemp_paramIndx)
    if (parFallback(absTolTemp_paramIndx(i))%default_val < 0.99_rkind*realMissing) then
      parFallback(absTolTemp_paramIndx(i))%default_val = default_absTol*absEnergyFac ! scale by absolute energy multiplier
    end if
  end do
  do i = 1, size(absTolWat_paramIndx)
    if (parFallback(absTolWat_paramIndx(i))%default_val < 0.99_rkind*realMissing) then
      parFallback(absTolWat_paramIndx(i))%default_val = default_absTol
    end if
  end do

  ! IDA Solver Parameters
  if (parFallback(iLookPARAM%idaMaxOrder)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%idaMaxOrder)%default_val = 5
  end if
  if (parFallback(iLookPARAM%idaMaxInternalSteps)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%idaMaxInternalSteps)%default_val = 999999 ! IDA default is 500, this is often too small for us
  end if
  if (parFallback(iLookPARAM%idaInitStepSize)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%idaInitStepSize)%default_val = 0
  end if
  if (parFallback(iLookPARAM%idaMinStepSize)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%idaMinStepSize)%default_val = 0 ! IDA default is 0
  end if
  if (parFallback(iLookPARAM%idaMaxStepSize)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%idaMaxStepSize)%default_val = 0 ! 0 means IDA's default of infinity
  end if
  if (parFallback(iLookPARAM%idaMaxErrTestFail)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%idaMaxErrTestFail)%default_val = 50 ! IDA default is 10
  end if
  if (parFallback(iLookPARAM%idaMaxDataWindowSteps)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%idaMaxDataWindowSteps)%default_val = 1.e6_rkind ! default is infinity, if 1e10 or larger then treat as infinity
  end if
  if (parFallback(iLookPARAM%idaDetectEvents)%default_val < 0.99_rkind*realMissing) then
    parFallback(iLookPARAM%idaDetectEvents)%default_val = 1._rkind ! default is to detect events (0 means do not detect events)
  end if
 end subroutine set_ida_defaults

 
 ! ************************************************************************************************
 ! Subroutine to set the FUSE default values if they are not already set
 ! ************************************************************************************************
 subroutine set_FUSE_defaults(parFallback, err, message)
  USE data_types       ,only:par_info                   ! data type for parameter constraints
  USE mDecisions_module,only:FUSEPRMS,FUSEAVIC,FUSETOPM ! model decision parameters 
  implicit none
  ! define dummy arguments
  type(par_info),intent(inout)           :: parFallback(:) ! default values and constraints of model parameters
  integer(i4b),intent(out)               :: err            ! error code
  character(*),intent(out)               :: message        ! error message
  !local variables
  logical(lgt)                           :: warning_flag   ! flag for warnings to standard output 

  ! initialize error control
  err=0
  message="set_FUSE_defaults/"
  warning_flag=.false.

  ! set FUSE parameter defaults for backwards compatibility
  if (parFallback(iLookPARAM%FUSE_Ac_max)%default_val == realMissing) then   ! FUSE PRMS max saturated area
   parFallback(iLookPARAM%FUSE_Ac_max)%default_val=0.95_rkind; warning_flag=.true.
  end if
  if (parFallback(iLookPARAM%FUSE_phi_tens)%default_val == realMissing) then ! FUSE PRMS tension storage fraction
   parFallback(iLookPARAM%FUSE_phi_tens)%default_val=0.5_rkind; warning_flag=.true.
  end if
  if (parFallback(iLookPARAM%FUSE_b)%default_val == realMissing) then        ! FUSE ARNO/VIC exponent
   parFallback(iLookPARAM%FUSE_b)%default_val=2._rkind; warning_flag=.true.
  end if
  if (parFallback(iLookPARAM%FUSE_lambda)%default_val == realMissing) then   ! FUSE TOPMODEL gamma distribution lambda parameter
   parFallback(iLookPARAM%FUSE_lambda)%default_val=7._rkind; warning_flag=.true.
  end if
  if (parFallback(iLookPARAM%FUSE_chi)%default_val == realMissing) then      ! FUSE TOPMODEL gamma distribution chi    parameter
   parFallback(iLookPARAM%FUSE_chi)%default_val=3._rkind; warning_flag=.true.
  end if
  if (parFallback(iLookPARAM%FUSE_mu)%default_val == realMissing) then       ! FUSE TOPMODEL gamma distribution mu     parameter
   parFallback(iLookPARAM%FUSE_mu)%default_val=3._rkind; warning_flag=.true.
  end if
  if (parFallback(iLookPARAM%FUSE_n)%default_val == realMissing) then        ! FUSE TOPMODEL exponent
   parFallback(iLookPARAM%FUSE_n)%default_val=4._rkind; warning_flag=.true.
  end if

  ! issue a warning if FUSE model decision choices used but default parameters not found in local parameters file
  if ((model_decisions(iLookDECISIONS%surfRun_SE)%iDecision == FUSEPRMS).or.&
     &(model_decisions(iLookDECISIONS%surfRun_SE)%iDecision == FUSEAVIC).or.&
     &(model_decisions(iLookDECISIONS%surfRun_SE)%iDecision == FUSETOPM)) then
     if (warning_flag) then
      print '(a136)', " WARNING: some FUSE parameters required by model decisions but are not in the local parameters file&
                      & -- default values have been assumed."
      print '(a1)',   " "
     end if
  end if

 end subroutine set_FUSE_defaults

end module read_pinit_module
