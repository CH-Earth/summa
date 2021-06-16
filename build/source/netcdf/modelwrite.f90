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

module modelwrite_module

! NetCDF types
USE netcdf
USE netcdf_util_module,only:netcdf_err                    ! netcdf error handling function

! top-level data types
USE nrtype

! missing values
USE globalData,only: integerMissing, realMissing

! provide access to global data
USE globalData,only:gru_struc                             ! gru->hru mapping structure

! provide access to the derived types to define the data structures
USE data_types,only:&
                    ! final data vectors
                    dlength,             & ! var%dat
                    ilength,             & ! var%dat
                    ! no spatial dimension
                    var_i,               & ! x%var(:)            (i4b)
                    var_i8,              & ! x%var(:)            integer(8)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    ! no variable dimension
                    hru_i,               & ! x%hru(:)            (i4b)
                    hru_d,               & ! x%hru(:)            (dp)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)     (i4b)
                    gru_double,          & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_int8,        & ! x%gru(:)%hru(:)%var(:)     integer(8)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec      ! x%gru(:)%hru(:)%var(:)%dat (dp)

! vector lengths
USE var_lookup, only: maxvarFreq ! number of output frequencies
USE var_lookup, only: maxvarStat ! number of statistics

implicit none
private
public::writeParm
public::writeData
public::writeBasin
public::writeTime
public::writeRestart
! define dimension lengths
integer(i4b),parameter      :: maxSpectral=2              ! maximum number of spectral bands
contains

 ! **********************************************************************************************************
 ! public subroutine writeParm: write model parameters
 ! **********************************************************************************************************
 subroutine writeParm(ispatial,struct,meta,err,message)
 USE globalData,only:ncid                        ! netcdf file ids
 USE data_types,only:var_info                    ! metadata info
 USE var_lookup,only:iLookStat                   ! index in statistics vector
 USE var_lookup,only:iLookFreq                   ! index in vector of model output frequencies
 implicit none

 ! declare input variables
 integer(i4b)  ,intent(in)   :: iSpatial         ! hydrologic response unit
 class(*)      ,intent(in)   :: struct           ! data structure
 type(var_info),intent(in)   :: meta(:)          ! metadata structure
 integer(i4b)  ,intent(out)  :: err              ! error code
 character(*)  ,intent(out)  :: message          ! error message
 ! local variables
 integer(i4b)                :: iVar             ! loop through variables

 ! initialize error control
 err=0;message="writeParm/"

 ! loop through local column model parameters
 do iVar = 1,size(meta)

  ! check that the variable is desired
  if (meta(iVar)%statIndex(iLookFREQ%timestep)==integerMissing) cycle

  ! initialize message
  message=trim(message)//trim(meta(iVar)%varName)//'/'

  ! HRU data
  if (iSpatial/=integerMissing) then
   select type (struct)
    class is (var_i)
     err = nf90_put_var(ncid(iLookFreq%timestep),meta(iVar)%ncVarID(iLookFreq%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
    class is (var_i8)
     err = nf90_put_var(ncid(iLookFreq%timestep),meta(iVar)%ncVarID(iLookFreq%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
    class is (var_d)
     err = nf90_put_var(ncid(iLookFreq%timestep),meta(iVar)%ncVarID(iLookFreq%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
    class is (var_dlength)
     err = nf90_put_var(ncid(iLookFreq%timestep),meta(iVar)%ncVarID(iLookFreq%timestep),(/struct%var(iVar)%dat/),start=(/iSpatial,1/),count=(/1,size(struct%var(iVar)%dat)/))
    class default; err=20; message=trim(message)//'unknown variable type (with HRU)'; return
   end select
   call netcdf_err(err,message); if (err/=0) return

  ! GRU data
  else
   select type (struct)
    class is (var_d)
     err = nf90_put_var(ncid(iLookFreq%timestep),meta(iVar)%ncVarID(iLookFreq%timestep),(/struct%var(iVar)/),start=(/1/),count=(/1/))
    class is (var_i8)
     err = nf90_put_var(ncid(iLookFreq%timestep),meta(iVar)%ncVarID(iLookFreq%timestep),(/struct%var(iVar)/),start=(/1/),count=(/1/))
    class default; err=20; message=trim(message)//'unknown variable type (no HRU)'; return
   end select
  end if
  call netcdf_err(err,message); if (err/=0) return

  ! re-initialize message
  message="writeParm/"
 end do  ! looping through local column model parameters

 end subroutine writeParm

 ! **************************************************************************************
 ! public subroutine writeData: write model time-dependent data
 ! **************************************************************************************
 subroutine writeData(finalizeStats,outputTimestep,nHRUrun,maxLayers,meta,stat,dat,map,indx,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE var_lookup,only:iLookIndex                     ! index into index structure
 USE var_lookup,only:iLookStat                      ! index into stat structure
 USE globalData,only:outFreq,ncid                   ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none
 ! declare dummy variables
 logical(lgt)  ,intent(in)        :: finalizeStats(:)  ! flags to finalize statistics
 integer(i4b)  ,intent(in)        :: outputTimestep(:) ! output time step
 integer(i4b)  ,intent(in)        :: nHRUrun           ! number of HRUs in the run domain
 integer(i4b)  ,intent(in)        :: maxLayers         ! maximum number of layers
 type(var_info),intent(in)        :: meta(:)           ! meta data
 class(*)      ,intent(in)        :: stat              ! stats data
 class(*)      ,intent(in)        :: dat               ! timestep data
 integer(i4b)  ,intent(in)        :: map(:)            ! map into stats child struct
 type(gru_hru_intVec) ,intent(in) :: indx              ! index data
 integer(i4b)  ,intent(out)       :: err               ! error code
 character(*)  ,intent(out)       :: message           ! error message
 ! local variables
 integer(i4b)                     :: iGRU              ! grouped response unit counter
 integer(i4b)                     :: iHRU              ! hydrologic response unit counter
 integer(i4b)                     :: iVar              ! variable index
 integer(i4b)                     :: iStat             ! statistics index
 integer(i4b)                     :: iFreq             ! frequency index
 integer(i4b)                     :: ncVarID           ! used only for time
 integer(i4b)                     :: nSnow             ! number of snow layers
 integer(i4b)                     :: nSoil             ! number of soil layers
 integer(i4b)                     :: nLayers           ! total number of layers
 ! output arrays
 integer(i4b)                     :: datLength         ! length of each data vector
 integer(i4b)                     :: maxLength         ! maximum length of each data vector
 real(rkind)                         :: realVec(nHRUrun)  ! real vector for all HRUs in the run domain
 real(rkind)                         :: realArray(nHRUrun,maxLayers+1)  ! real array for all HRUs in the run domain
 integer(i4b)                     :: intArray(nHRUrun,maxLayers+1)   ! integer array for all HRUs in the run domain
 integer(i4b)                     :: dataType          ! type of data
 integer(i4b),parameter           :: ixInteger=1001    ! named variable for integer
 integer(i4b),parameter           :: ixReal=1002       ! named variable for real
 ! initialize error control
 err=0;message="writeData/"

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! skip frequencies that are not needed
  if(.not.outFreq(iFreq)) cycle

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  do iVar = 1,size(meta)

   ! handle time first
   if (meta(iVar)%varName=='time')then
    ! get variable index
    err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID)
    call netcdf_err(err,message); if (err/=0) return
    ! define HRUs and GRUs (only write once)
    iGRU=1; iHRU=1
    ! data bound write
    select type(dat) ! forcStruc
     class is (gru_hru_double)   ! x%gru(:)%hru(:)%var(:)
      err = nf90_put_var(ncid(iFreq),ncVarID,(/dat%gru(iGRU)%hru(iHRU)%var(iVar)/),start=(/outputTimestep(iFreq)/),count=(/1/))
      call netcdf_err(err,message); if (err/=0) return
      cycle ! move onto the next variable
     class default; err=20; message=trim(message)//'time variable must be of type gru_hru_double (forcing data structure)'; return
    end select
   end if  ! id time

   ! define the statistics index
   iStat = meta(iVar)%statIndex(iFreq)

   ! check that the variable is desired
   if (iStat==integerMissing.or.trim(meta(iVar)%varName)=='unknown') cycle

   ! stats output: only scalar variable type
   if(meta(iVar)%varType==iLookVarType%scalarv) then
    select type(stat)
     class is (gru_hru_doubleVec)

      ! loop through HRUs and GRUs, and place data in the single vector
      do iGRU=1,size(gru_struc)
       do iHRU=1,gru_struc(iGRU)%hruCount
        realVec(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix) = stat%gru(iGRU)%hru(iHRU)%var(map(iVar))%dat(iFreq)
       end do
      end do

      ! write data
      err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realVec,start=(/1,outputTimestep(iFreq)/),count=(/nHRUrun,1/))

     class default; err=20; message=trim(message)//'stats must be scalarv and of type gru_hru_doubleVec'; return
    end select  ! stat

   ! non-scalar variables: regular data structures
   else

    ! initialize the data vectors
    select type (dat)
     class is (gru_hru_doubleVec); realArray(:,:) = realMissing;    dataType=ixReal
     class is (gru_hru_intVec);     intArray(:,:) = integerMissing; dataType=ixInteger
     class default; err=20; message=trim(message)//'data must not be scalarv and either of type gru_hru_doubleVec or gru_hru_intVec'; return
    end select

    ! loop thru GRUs and HRUs
    do iGRU=1,size(gru_struc)
     do iHRU=1,gru_struc(iGRU)%hruCount

      ! get the model layers
      nSoil   = indx%gru(iGRU)%hru(iHRU)%var(iLookIndex%nSoil)%dat(1)
      nSnow   = indx%gru(iGRU)%hru(iHRU)%var(iLookIndex%nSnow)%dat(1)
      nLayers = indx%gru(iGRU)%hru(iHRU)%var(iLookIndex%nLayers)%dat(1)

      ! get the length of each data vector
      select case (meta(iVar)%varType)
       case(iLookVarType%wLength); datLength = maxSpectral
       case(iLookVarType%midToto); datLength = nLayers
       case(iLookVarType%midSnow); datLength = nSnow
       case(iLookVarType%midSoil); datLength = nSoil
       case(iLookVarType%ifcToto); datLength = nLayers+1
       case(iLookVarType%ifcSnow); datLength = nSnow+1
       case(iLookVarType%ifcSoil); datLength = nSoil+1
       case default; cycle
      end select ! vartype

      ! get the data vectors
      select type (dat)
       class is (gru_hru_doubleVec); realArray(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = dat%gru(iGRU)%hru(iHRU)%var(iVar)%dat(:)
       class is (gru_hru_intVec);     intArray(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = dat%gru(iGRU)%hru(iHRU)%var(iVar)%dat(:)
       class default; err=20; message=trim(message)//'data must not be scalarv and either of type gru_hru_doubleVec or gru_hru_intVec'; return
      end select

     end do  ! HRU loop
    end do  ! GRU loop

    ! get the maximum length of each data vector
    select case (meta(iVar)%varType)
     case(iLookVarType%wLength); maxLength = maxSpectral
     case(iLookVarType%midToto); maxLength = maxLayers
     case(iLookVarType%midSnow); maxLength = maxLayers-nSoil
     case(iLookVarType%midSoil); maxLength = nSoil
     case(iLookVarType%ifcToto); maxLength = maxLayers+1
     case(iLookVarType%ifcSnow); maxLength = (maxLayers-nSoil)+1
     case(iLookVarType%ifcSoil); maxLength = nSoil+1
     case default; cycle
    end select ! vartype

    ! write the data vectors
    select case(dataType)
     case(ixReal);    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realArray(1:nHRUrun,1:maxLength),start=(/1,1,outputTimestep(iFreq)/),count=(/nHRUrun,maxLength,1/))
     case(ixInteger); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),intArray(1:nHRUrun,1:maxLength),start=(/1,1,outputTimestep(iFreq)/),count=(/nHRUrun,maxLength,1/))
     case default; err=20; message=trim(message)//'data must be of type integer or real'; return
    end select ! data type

   end if ! not scalarv

   ! process error code
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)//'_'//trim(get_statName(iStat))
   call netcdf_err(err,message); if (err/=0) return

  end do ! iVar
 end do ! iFreq

 end subroutine writeData

 ! **************************************************************************************
 ! public subroutine writeBasin: write basin-average variables
 ! **************************************************************************************
 subroutine writeBasin(iGRU,finalizeStats,outputTimestep,meta,stat,dat,map,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE globalData,only:outFreq,ncid                   ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none

 ! declare dummy variables
 integer(i4b)  ,intent(in)     :: iGRU              ! GRU index
 logical(lgt)  ,intent(in)     :: finalizeStats(:)  ! flags to finalize statistics
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 type(var_info),intent(in)     :: meta(:)           ! meta data
 type(dlength) ,intent(in)     :: stat(:)           ! stats data
 type(dlength) ,intent(in)     :: dat(:)            ! timestep data
 integer(i4b)  ,intent(in)     :: map(:)            ! map into stats child struct
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iStat             ! statistics index
 integer(i4b)                  :: iFreq             ! frequency index
 ! initialize error control
 err=0;message="f-writeBasin/"

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! skip frequencies that are not needed
  if(.not.outFreq(iFreq)) cycle

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  do iVar = 1,size(meta)

   ! define the statistics index
   iStat = meta(iVar)%statIndex(iFreq)

   ! check that the variable is desired
   if (iStat==integerMissing.or.trim(meta(iVar)%varName)=='unknown') cycle

   ! stats/data output - select data type
   select case (meta(iVar)%varType)

    case (iLookVarType%scalarv)
     err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),(/stat(map(iVar))%dat(iFreq)/),start=(/iGRU,outputTimestep(iFreq)/),count=(/1,1/))

    case (iLookVarType%routing)
     if (iFreq==1 .and. outputTimestep(iFreq)==1) then
      err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),(/dat(iVar)%dat/),start=(/1/),count=(/1000/))
     end if

    case default
     err=40; message=trim(message)//"unknownVariableType[name='"//trim(meta(iVar)%varName)//"';type='"//trim(get_varTypeName(meta(iVar)%varType))//    "']"; return
   end select ! variable type

   ! process error code
   if (err.ne.0) message=trim(message)//trim(meta(iVar)%varName)//'_'//trim(get_statName(iStat))
   call netcdf_err(err,message); if (err/=0) return

  end do ! iVar
 end do ! iFreq

 end subroutine writeBasin

 ! **************************************************************************************
 ! public subroutine writeTime: write current time to all files
 ! **************************************************************************************
 subroutine writeTime(finalizeStats,outputTimestep,meta,dat,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE globalData,only:ncid                           ! output file IDs
 USE var_lookup,only:iLookStat                      ! index into stat structure
 implicit none

 ! declare dummy variables
 logical(lgt)  ,intent(in)     :: finalizeStats(:)  ! flags to finalize statistics
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 type(var_info),intent(in)     :: meta(:)           ! meta data
 integer       ,intent(in)     :: dat(:)            ! timestep data
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iFreq             ! frequency index
 integer(i4b)                  :: ncVarID           ! used only for time
 ! initialize error control
 err=0;message="f-writeTime/"

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  do iVar = 1,size(meta)

   ! check instantaneous
   if (meta(iVar)%statIndex(iFreq)/=iLookStat%inst) cycle

   ! get variable id in file
   err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID)
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
   call netcdf_err(err,message)
   if (err/=0) then; err=20; return; end if

   ! add to file
   err = nf90_put_var(ncid(iFreq),ncVarID,(/dat(iVar)/),start=(/outputTimestep(iFreq)/),count=(/1/))
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
   call netcdf_err(err,message)
   if (err/=0) then; err=20; return; end if

  end do ! iVar
 end do ! iFreq

 end subroutine writeTime

 ! *********************************************************************************************************
 ! public subroutine printRestartFile: print a re-start file
 ! *********************************************************************************************************
 subroutine writeRestart(filename,         & ! intent(in): name of restart file
                         nGRU,             & ! intent(in): number of GRUs
                         nHRU,             & ! intent(in): number of HRUs
                         prog_meta,        & ! intent(in): prognostics metadata
                         prog_data,        & ! intent(in): prognostics data
                         bvar_meta,        & ! intent(in): basin (gru) variable metadata
                         bvar_data,        & ! intent(in): basin (gru) variable data
                         maxLayers,        & ! intent(in): maximum number of layers
                         maxSnowLayers,    & ! intent(in): maximum number of snow layers
                         indx_meta,        & ! intent(in): index metadata
                         indx_data,        & ! intent(in): index data
                         err,message)        ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! access the derived types to define the data structures
 USE data_types,only:var_info               ! metadata
 ! access named variables defining elements in the data structures
 USE var_lookup,only:iLookINDEX             ! named variables for structure elements
 USE var_lookup,only:iLookVarType           ! named variables for structure elements
 USE var_lookup,only:iLookBVAR              ! named variables for structure elements
 ! constants
 USE globalData,only:gru_struc              ! gru-hru mapping structures
 ! external routines
 USE netcdf_util_module,only:nc_file_close  ! close netcdf file
 USE netcdf_util_module,only:nc_file_open   ! open netcdf file
 USE globalData,only:nTimeDelay             ! number of timesteps in the time delay histogram
 
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input
 character(len=256),intent(in)      :: filename      ! name of the restart file
 integer(i4b),intent(in)            :: nGRU          ! number of GRUs
 integer(i4b),intent(in)            :: nHRU          ! number of HRUs
 type(var_info),intent(in)          :: prog_meta(:)  ! prognostic variable metadata
 type(gru_hru_doubleVec),intent(in) :: prog_data     ! prognostic vars
 type(var_info),intent(in)          :: bvar_meta(:)  ! basin variable metadata
 type(gru_doubleVec),intent(in)     :: bvar_data     ! basin variables
 type(var_info),intent(in)          :: indx_meta(:)  ! metadata
 type(gru_hru_intVec),intent(in)    :: indx_data     ! indexing vars
 ! output: error control
 integer(i4b),intent(out)           :: err           ! error code
 character(*),intent(out)           :: message       ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! dummy variables
 integer(i4b), intent(in)           :: maxLayers     ! maximum number of total layers
 integer(i4b), intent(in)           :: maxSnowLayers ! maximum number of snow layers

 ! local variables
 integer(i4b)                       :: ncid          ! netcdf file id
 integer(i4b),allocatable           :: ncVarID(:)    ! netcdf variable id
 integer(i4b)                       :: ncSnowID      ! index variable id
 integer(i4b)                       :: ncSoilID      ! index variable id

 integer(i4b)                       :: nSoil         ! number of soil layers
 integer(i4b)                       :: nSnow         ! number of snow layers
 integer(i4b)                       :: maxSnow       ! maximum number of snow layers
 integer(i4b)                       :: maxSoil       ! maximum number of soil layers
 integer(i4b)                       :: nLayers       ! number of total layers
 integer(i4b),parameter             :: nSpectral=2   ! number of spectal bands
 integer(i4b),parameter             :: nScalar=1     ! size of a scalar
 integer(i4b)                       :: nProgVars     ! number of prognostic variables written to state file

 integer(i4b)                       :: hruDimID      ! variable dimension ID
 integer(i4b)                       :: gruDimID      ! variable dimension ID
 integer(i4b)                       :: tdhDimID      ! variable dimension ID
 integer(i4b)                       :: scalDimID     ! variable dimension ID
 integer(i4b)                       :: specDimID     ! variable dimension ID
 integer(i4b)                       :: midSnowDimID  ! variable dimension ID
 integer(i4b)                       :: midSoilDimID  ! variable dimension ID
 integer(i4b)                       :: midTotoDimID  ! variable dimension ID
 integer(i4b)                       :: ifcSnowDimID  ! variable dimension ID
 integer(i4b)                       :: ifcSoilDimID  ! variable dimension ID
 integer(i4b)                       :: ifcTotoDimID  ! variable dimension ID

 character(len=32),parameter        :: hruDimName    ='hru'      ! dimension name for HRUs
 character(len=32),parameter        :: gruDimName    ='gru'      ! dimension name for GRUs
 character(len=32),parameter        :: tdhDimName    ='tdh'      ! dimension name for time-delay basin variables
 character(len=32),parameter        :: scalDimName   ='scalarv'  ! dimension name for scalar data
 character(len=32),parameter        :: specDimName   ='spectral' ! dimension name for spectral bands
 character(len=32),parameter        :: midSnowDimName='midSnow'  ! dimension name for snow-only layers
 character(len=32),parameter        :: midSoilDimName='midSoil'  ! dimension name for soil-only layers
 character(len=32),parameter        :: midTotoDimName='midToto'  ! dimension name for layered varaiables
 character(len=32),parameter        :: ifcSnowDimName='ifcSnow'  ! dimension name for snow-only layers
 character(len=32),parameter        :: ifcSoilDimName='ifcSoil'  ! dimension name for soil-only layers
 character(len=32),parameter        :: ifcTotoDimName='ifcToto'  ! dimension name for layered variables

 integer(i4b)                       :: cHRU          ! count of HRUs
 integer(i4b)                       :: iHRU          ! index of HRUs
 integer(i4b)                       :: iGRU          ! index of GRUs
 integer(i4b)                       :: iVar          ! variable index
 logical(lgt)                       :: okLength      ! flag to check if the vector length is OK
 character(len=256)                 :: cmessage      ! downstream error message
 ! --------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='writeRestart/'

 ! size of prognostic variable vector
 nProgVars = size(prog_meta)
 allocate(ncVarID(nProgVars+1))     ! include 1 additional basin variable in ID array (possibly more later)

 ! maximum number of soil layers
 maxSoil = gru_struc(1)%hruInfo(1)%nSoil

 ! maximum number of snow layers
 maxSnow = maxSnowLayers
 
 ! create file
 err = nf90_create(trim(filename),nf90_classic_model,ncid)
 message='iCreate[create]'; call netcdf_err(err,message); if(err/=0)return

 ! define dimensions
                err = nf90_def_dim(ncid,trim(hruDimName)    ,nHRU       ,    hruDimID); message='iCreate[hru]'     ; call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(gruDimName)    ,nGRU       ,    gruDimID); message='iCreate[gru]'     ; call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(tdhDimName)    ,nTimeDelay ,    tdhDimID); message='iCreate[tdh]'     ; call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(scalDimName)   ,nScalar    ,   scalDimID); message='iCreate[scalar]'  ; call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(specDimName)   ,nSpectral  ,   specDimID); message='iCreate[spectral]'; call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(midSoilDimName),maxSoil    ,midSoilDimID); message='iCreate[ifcSoil]' ; call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(midTotoDimName),maxLayers  ,midTotoDimID); message='iCreate[midToto]' ; call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(ifcSoilDimName),maxSoil+1  ,ifcSoilDimID); message='iCreate[ifcSoil]' ; call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(ifcTotoDimName),maxLayers+1,ifcTotoDimID); message='iCreate[ifcToto]' ; call netcdf_err(err,message); if(err/=0)return
 if (maxSnow>0) err = nf90_def_dim(ncid,trim(midSnowDimName),maxSnow    ,midSnowDimID); message='iCreate[ifcSnow]' ; call netcdf_err(err,message); if(err/=0)return
 if (maxSnow>0) err = nf90_def_dim(ncid,trim(ifcSnowDimName),maxSnow+1  ,ifcSnowDimID); message='iCreate[ifcSnow]' ; call netcdf_err(err,message); if(err/=0)return
 ! re-initialize error control
 err=0; message='writeRestart/'

 ! define prognostic variables
 do iVar = 1,nProgVars
  if (prog_meta(iVar)%varType==iLookvarType%unknown) cycle

  ! define variable
  select case(prog_meta(iVar)%varType)
   case(iLookvarType%scalarv);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,  scalDimID /),ncVarID(iVar))
   case(iLookvarType%wLength);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,  specDimID /),ncVarID(iVar))
   case(iLookvarType%midSoil);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,midSoilDimID/),ncVarID(iVar))
   case(iLookvarType%midToto);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,midTotoDimID/),ncVarID(iVar))
   case(iLookvarType%ifcSoil);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,ifcSoilDimID/),ncVarID(iVar))
   case(iLookvarType%ifcToto);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,ifcTotoDimID/),ncVarID(iVar))
   case(iLookvarType%midSnow); if (maxSnow>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,midSnowDimID/),ncVarID(iVar))
   case(iLookvarType%ifcSnow); if (maxSnow>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,ifcSnowDimID/),ncVarID(iVar))
  end select

  ! check errors
  if(err/=0)then
   message=trim(message)//trim(cmessage)//' [variable '//trim(prog_meta(iVar)%varName)//']'
   return
  end if

  ! add parameter description
  err = nf90_put_att(ncid,ncVarID(iVar),'long_name',trim(prog_meta(iVar)%vardesc))
  call netcdf_err(err,message)

  ! add parameter units
  err = nf90_put_att(ncid,ncVarID(iVar),'units',trim(prog_meta(iVar)%varunit))
  call netcdf_err(err,message)

 end do ! iVar
 
 ! define selected basin variables (derived) -- e.g., hillslope routing
 err = nf90_def_var(ncid, trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varName), nf90_double, (/gruDimID, tdhDimID /), ncVarID(nProgVars+1))
 err = nf90_put_att(ncid,ncVarID(nProgVars+1),'long_name',trim(bvar_meta(iLookBVAR%routingRunoffFuture)%vardesc));   call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncVarID(nProgVars+1),'units'    ,trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varunit));   call netcdf_err(err,message)

 ! define index variables - snow
 err = nf90_def_var(ncid,trim(indx_meta(iLookIndex%nSnow)%varName),nf90_int,(/hruDimID/),ncSnowID); call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSnowID,'long_name',trim(indx_meta(iLookIndex%nSnow)%vardesc));           call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSnowID,'units'    ,trim(indx_meta(iLookIndex%nSnow)%varunit));           call netcdf_err(err,message)

 ! define index variables - soil
 err = nf90_def_var(ncid,trim(indx_meta(iLookIndex%nSoil)%varName),nf90_int,(/hruDimID/),ncSoilID); call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSoilID,'long_name',trim(indx_meta(iLookIndex%nSoil)%vardesc));           call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSoilID,'units'    ,trim(indx_meta(iLookIndex%nSoil)%varunit));           call netcdf_err(err,message)

 ! end definition phase
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return

 ! write variables
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   cHRU = gru_struc(iGRU)%hruInfo(iHRU)%hru_ix
   do iVar = 1,size(prog_meta)

    ! excape if this variable is not used
    if (prog_meta(iVar)%varType==iLookvarType%unknown) cycle

    ! actual number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
    nLayers = nSoil + nSnow

    ! check size
    ! NOTE: this may take time that we do not wish to use
    okLength=.true.
    select case (prog_meta(iVar)%varType)
     case(iLookVarType%scalarv);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nScalar  )
     case(iLookVarType%wlength);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSpectral)
     case(iLookVarType%midSoil);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSoil    )
     case(iLookVarType%midToto);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nLayers  )
     case(iLookVarType%ifcSoil);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSoil+1  )
     case(iLookVarType%ifcToto);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nLayers+1)
     case(iLookVarType%midSnow); if (nSnow>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSnow    )
     case(iLookVarType%ifcSnow); if (nSnow>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSnow+1  )
     case default; err=20; message=trim(message)//'unknown var type'; return
    end select

    ! error check
    if(.not.okLength)then
     message=trim(message)//'bad vector length for variable '//trim(prog_meta(iVar)%varname)
     err=20; return
    endif

    ! write data
    select case (prog_meta(iVar)%varType)
     case(iLookVarType%scalarv);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nScalar  /))
     case(iLookVarType%wlength);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSpectral/))
     case(iLookVarType%midSoil);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSoil    /))
     case(iLookVarType%midToto);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nLayers  /))
     case(iLookVarType%ifcSoil);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSoil+1  /))
     case(iLookVarType%ifcToto);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nLayers+1/))
     case(iLookVarType%midSnow); if (nSnow>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSnow    /))
     case(iLookVarType%ifcSnow); if (nSnow>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSnow+1  /))
     case default; err=20; message=trim(message)//'unknown var type'; return
    end select

    ! error check
    if (err.ne.0) message=trim(message)//'writing variable:'//trim(prog_meta(iVar)%varName)
    call netcdf_err(err,message); if (err/=0) return
    err=0; message='writeRestart/'

   end do ! iVar loop

   ! write index variables
   err=nf90_put_var(ncid,ncSnowID,(/indx_data%gru(iGRU)%hru(iHRU)%var(iLookIndex%nSnow)%dat/),start=(/cHRU/),count=(/1/))
   err=nf90_put_var(ncid,ncSoilID,(/indx_data%gru(iGRU)%hru(iHRU)%var(iLookIndex%nSoil)%dat/),start=(/cHRU/),count=(/1/))

  end do ! iHRU loop
  
  ! write selected basin variables
  err=nf90_put_var(ncid,ncVarID(nProgVars+1),(/bvar_data%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat/),  start=(/iGRU/),count=(/1,nTimeDelay/))
  
 end do  ! iGRU loop

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(ncVarID)

 end subroutine writeRestart

end module modelwrite_module
