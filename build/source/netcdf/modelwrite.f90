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
USE netcdf_util_module,only:netcdf_err  ! netcdf error handling function

! top-level data types
USE nr_type

! missing values
USE globalData,only: integerMissing, realMissing

! output constraints
USE globalData,only:maxSnowLayers       ! maximum number of snow layers
USE globalData,only:maxSoilLayers       ! maximum number of soil layers
USE globalData,only:maxLayers           ! maximum number of layers
USE globalData,only:nTimeDelay          ! number of timesteps in the time delay histogram
USE globalData,only:nSpecBand           ! maximum number of spectral bands
USE globalData,only:allowRoutingOutput  ! flag to allow routing variable output

! provide access to global data
USE globalData,only:nGRUrun             ! number of GRUs in the run (rank)
USE globalData,only:nHRUrun             ! number of HRUs in the run (rank)
USE globalData,only:gru_struc           ! gru->hru mapping structure

! provide access to the derived types to define the data structures
USE data_types,only:&
                    ! final data vectors
                    dlength,             & ! var%dat
                    ilength,             & ! var%dat
                    ! no spatial dimension
                    var_i,               & ! x%var(:)                   (i4b)
                    var_i8,              & ! x%var(:)                   (i8b)
                    var_d,               & ! x%var(:)                   (rkind)
                    var_ilength,         & ! x%var(:)%dat               (i4b)
                    var_dlength,         & ! x%var(:)%dat               (rkind)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)            (i4b)
                    gru_double,          & ! x%gru(:)%var(:)            (rkind)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat        (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat        (rkind)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (rkind)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec      ! x%gru(:)%hru(:)%var(:)%dat (rkind)

! vector lengths
USE var_lookup, only: maxvarFreq ! number of output frequencies
USE var_lookup, only: maxvarStat ! number of statistics   

implicit none
private
public::writeParam
public::writeData_fullSeries
public::writeData_perStep
public::writeTime
public::writeRestart

contains

 ! **********************************************************************************************************
 ! public subroutine writeParam: write model parameters
 ! **********************************************************************************************************
 subroutine writeParam(iSpatial,struct,meta,err,message)
 USE globalData,only:ncid                        ! netcdf file ids
 USE data_types,only:var_info                    ! metadata info
 USE var_lookup,only:iLookSTAT                   ! index in statistics vector
 USE var_lookup,only:iLookFREQ                   ! index in vector of model output frequencies
 implicit none

 ! declare input variables
 integer(i4b)  ,intent(in)   :: iSpatial         ! HRU index or GRU index
 class(*)      ,intent(in)   :: struct           ! data structure
 type(var_info),intent(in)   :: meta(:)          ! metadata structure
 integer(i4b)  ,intent(out)  :: err              ! error code
 character(*)  ,intent(out)  :: message          ! error message
 ! local variables
 integer(i4b)                :: iVar             ! loop through variables

 ! initialize error control
 err=0;message="writeParam/"

 ! loop through local column model parameters
 do iVar = 1,size(meta)

  ! check that the variable is desired
  if (meta(iVar)%statIndex(iLookFREQ%timestep)==integerMissing) cycle

  ! initialize message
  message=trim(message)//trim(meta(iVar)%varName)//':'

  select type (struct)
   class is (var_i)
    err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
   class is (var_i8)
    err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
   class is (var_d)
    err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
   class is (var_dlength)
    err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)%dat/),start=(/iSpatial,1/),count=(/1,size(struct%var(iVar)%dat)/))
   class default; err=20; message=trim(message)//'parameter type must be var_i, var_i8, var_d, or var_dlength'; return
  end select
  call netcdf_err(err,message); if (err/=0) return

  ! re-initialize message
  message="writeParam/"
 end do  ! looping through local column model parameters

 end subroutine writeParam

 ! **************************************************************************************
 ! public subroutine writeData_fullSeries: write buffered model time-dependent data for each HRU
 ! **************************************************************************************
 subroutine writeData_fullSeries(finalizeStats,maxWrite,meta,datt,map,indx,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE var_lookup,only:iLookFREQ                      ! index into freq structure
 USE globalData,only:outFreq,ncid                   ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none
 ! declare dummy variables
 logical(lgt)  ,intent(in)      :: finalizeStats(:)               ! flags to finalize statistics
 integer(i4b)  ,intent(in)      :: maxWrite                       ! maximum number of steps written
 type(var_info),intent(in)      :: meta(:)                        ! meta data
 class(*)      ,intent(in)      :: datt(:)                        ! timestep or buffer data
 integer(i4b)  ,intent(in)      :: map(:)                         ! map into stats child struct
 type(gru_hru_intVec),intent(in):: indx                           ! index data
 integer(i4b)  ,intent(out)     :: err                            ! error code
 character(*)  ,intent(out)     :: message                        ! error message
 ! local variables
 integer(i4b)                   :: iGRU                           ! grouped response unit counter
 integer(i4b)                   :: iHRU                           ! hydrologic response unit counter
 integer(i4b)                   :: iVar                           ! variable index
 integer(i4b)                   :: iStat                          ! statistics index
 integer(i4b)                   :: iFreq                          ! frequency index
 integer(i4b)                   :: iTime                          ! time index
 integer(i4b)                   :: ncVarID                        ! used only for time
 integer(i4b)                   :: nSpace                         ! number of spatial data elements
 ! output arrays
 real(rkind)                    :: timeBuffer(maxWrite)           ! buffer for all time steps
 real(rkind)                    :: realBuffer(nHRUrun,maxWrite)   ! buffer for all HRUs in the run domain + time steps

 ! initialize error control
 err=0

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! skip frequencies that are not needed
  if(.not.outFreq(iFreq)) cycle

  ! restrict attention to the timestep data if buffered write
  if(iFreq/=iLookFREQ%timestep) cycle

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  iVarLoop: do iVar = 1,size(meta)

   ! initialize message
   message="writeData_fullSeries/"//trim(meta(iVar)%varName)

   ! ****************************************************************************
   ! *** write time information -- instantaneous
   ! ****************************************************************************

   ! handle time first
   if(trim(meta(iVar)%varName)=='time')then
    message=trim(message)//':' ! add statistic (none) to message 

    ! get variable index
    err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID)
    call netcdf_err(err,message); if (err/=0) return

    ! define HRUs and GRUs (only write once)
    iGRU=1; iHRU=1

    ! data bound array access
    select type (datt) ! forcStruc
     class is (gru_hru_double) ! x%gru(:)%hru(:)%var(:)
      do iTime=1,maxWrite
       timeBuffer(iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%var(iVar)
      end do
     class default; err=20; message=trim(message)//'time variable must be of type gru_hru_double (forcing data structure)'; return
    end select  ! type of data structure

    ! write time
    err = nf90_put_var(ncid(iFreq),ncVarID,(/timeBuffer/),start=(/1/),count=(/maxWrite/))
    call netcdf_err(err,message); if (err/=0) return
    cycle ! move onto the next variable

   end if  ! if time

   ! ****************************************************************************
   ! *** write scalar variables
   ! ****************************************************************************

   ! define the statistics index
   iStat = meta(iVar)%statIndex(iFreq)
   message=trim(message)//'_'//trim(get_statName(iStat))//':' ! add statistic to message

   ! check that the variable is desired, currently do not write large variables (unknown and routing) as they are large and slow things down a lot
   if (iStat==integerMissing .or. meta(iVar)%varType==iLookVarType%unknown .or. meta(iVar)%varType==integerMissing) cycle
   if (meta(iVar)%varType==iLookVarType%routing .and. .not.allowRoutingOutput) cycle ! routing variable write can be turned on with the allowRoutingOutput flag

   ! buffered output: only scalar variable type
   if(meta(iVar)%varType==iLookVarType%scalarv) then

    ! initialize the data vectors
    select type (datt)
     class is (gru_hru_double); nSpace = nHRUrun; realBuffer(:,:) = realMissing
     class is (gru_hru_int);    nSpace = nHRUrun; realBuffer(:,:) = realMissing
     class is (gru_double);     nSpace = nGRUrun; realBuffer(:,:) = realMissing
     class is (gru_int);        nSpace = nGRUrun; realBuffer(:,:) = realMissing
     class default; err=20; message=trim(message)//'data is not scalarv so should be either of type gru_hru_[double or int] or gru_[double or int]'; return
    end select

    ! loop through time, HRUs and GRU
    do iTime=1,maxWrite
     do iGRU=1,size(gru_struc)
      do iHRU=1,gru_struc(iGRU)%hruCount

       ! get the data vectors
       select type (datt)
        class is (gru_hru_double); realBuffer(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%var(map(iVar))
        class is (gru_hru_int);    realBuffer(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%var(map(iVar))
        class is (gru_double); realBuffer(iGRU,iTime) = datt(iTime)%gru(iGRU)%var(map(iVar)); exit ! only need to get the GRU-level data once
        class is (gru_int);    realBuffer(iGRU,iTime) = datt(iTime)%gru(iGRU)%var(map(iVar)); exit ! only need to get the GRU-level data once
       end select ! time step data structure

      end do  ! HRU loop
     end do  ! GRU loop
    end do  ! time loop

    ! write the data vectors
    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nSpace,1:maxWrite),start=(/1,1/),count=(/nSpace,maxWrite/))
    call netcdf_err(err,message); if (err/=0) return

   else ! cannot write non-scalar variables in buffered write -- too complicated and slow, so not currently supported
    write(*,*)'WARNING: cannot output non-scalar type data when using the buffered write option (writeFullSeries), skipping variable '//trim(meta(iVar)%varName); cycle

   end if ! not scalarv

  end do iVarLoop ! iVar
 end do ! iFreq

 end subroutine writeData_fullSeries

 ! **************************************************************************************
 ! public subroutine writeData_perStep: write per-step model time-dependent data for each HRU
 ! **************************************************************************************
 subroutine writeData_perStep(finalizeStats,outputTimestep,maxLengthAll,meta,stat,datt,map,indx,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE var_lookup,only:maxvarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE var_lookup,only:iLookINDEX                     ! index into index structure
 USE var_lookup,only:iLookFREQ                      ! index into freq structure
 USE globalData,only:outFreq,ncid                   ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none
 ! declare dummy variables
 logical(lgt)  ,intent(in)      :: finalizeStats(:)               ! flags to finalize statistics
 integer(i4b)  ,intent(in)      :: outputTimestep(:)              ! output time step
 integer(i4b)  ,intent(in)      :: maxLengthAll                   ! maxLength all data
 type(var_info),intent(in)      :: meta(:)                        ! meta data
 class(*)      ,intent(in)      :: stat                           ! stats data
 class(*)      ,intent(in)      :: datt                           ! timestep data
 integer(i4b)  ,intent(in)      :: map(:)                         ! map into stats child struct
 type(gru_hru_intVec),intent(in):: indx                           ! index data
 integer(i4b)  ,intent(out)     :: err                            ! error code
 character(*)  ,intent(out)     :: message                        ! error message
 ! local variables
 integer(i4b)                   :: iGRU                           ! grouped response unit counter
 integer(i4b)                   :: iHRU                           ! hydrologic response unit counter
 integer(i4b)                   :: iVar                           ! variable index
 integer(i4b)                   :: iStat                          ! statistics index
 integer(i4b)                   :: iFreq                          ! frequency index
 integer(i4b)                   :: ncVarID                        ! used only for time
 integer(i4b)                   :: nSnow                          ! number of snow layers
 integer(i4b)                   :: nSoil                          ! number of soil layers
 integer(i4b)                   :: nLayers                        ! total number of layers
 integer(i4b)                   :: nSpace                         ! number of spatial data elements
 ! output arrays
 integer(i4b)                   :: datLength                      ! length of each data vector
 integer(i4b)                   :: maxLength                      ! maximum length of each data vector
 real(rkind)                    :: timeStep                       ! timestep value written to file
 real(rkind)                    :: realBuffer(nHRUrun)            ! buffer for all HRUs in the run domain + time steps
 real(rkind)                    :: realArray(nHRUrun,maxLengthAll)! real array for all HRUs in the run domain
 integer(i4b)                   :: intArray(nHRUrun,maxLengthAll) ! integer array for all HRUs in the run domain
 integer(i4b)                   :: dataType                       ! type of data
 integer(i4b),parameter         :: ixInteger=1001                 ! named variable for integer
 integer(i4b),parameter         :: ixReal=1002                    ! named variable for real

 ! initialize error control
 err=0

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! skip frequencies that are not needed
  if(.not.outFreq(iFreq)) cycle

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  iVarLoop: do iVar = 1,size(meta)

   ! initialize message
   message="writeData_perStep/"//trim(meta(iVar)%varName)

   ! ****************************************************************************
   ! *** write time information -- instantaneous
   ! ****************************************************************************

   ! handle time first
   if(trim(meta(iVar)%varName)=='time')then
    message=trim(message)//':' ! add statistic (none) to message

    ! get variable index
    err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID)
    call netcdf_err(err,message); if (err/=0) return

    ! define HRUs and GRUs (only write once)
    iGRU=1; iHRU=1

    ! data bound array access
    select type (datt) ! forcStruc
     class is (gru_hru_double) ! x%gru(:)%hru(:)%var(:)
      timeStep = datt%gru(iGRU)%hru(iHRU)%var(iVar)
     class default; err=20; message=trim(message)//'time variable must be of type gru_hru_double (forcing data structure)'; return
    end select  ! type of data structure

    ! write time
    err = nf90_put_var(ncid(iFreq),ncVarID,(/timeStep/),start=(/outputTimestep(iFreq)/),count=(/1/))
    call netcdf_err(err,message); if (err/=0) return
    cycle ! move onto the next variable

   end if  ! if time

   ! ****************************************************************************
   ! *** write scalar variables
   ! ****************************************************************************

   ! define the statistics index
   iStat = meta(iVar)%statIndex(iFreq)
   message=trim(message)//'_'//trim(get_statName(iStat))//':' ! add statistic to message

   ! check that the variable is desired, currently do not write large variables (unknown and routing) as they are large and slow things down a lot
   if (iStat==integerMissing .or. meta(iVar)%varType==iLookVarType%unknown .or. meta(iVar)%varType==integerMissing) cycle
   if (meta(iVar)%varType==iLookVarType%routing .and. .not.allowRoutingOutput) cycle ! routing variable write can be turned on with the allowRoutingOutput flag

   ! stats output: only scalar variable type
   if(meta(iVar)%varType==iLookVarType%scalarv) then

    ! initialize the data vectors
    select type (stat)
     class is (gru_hru_doubleVec); nSpace = nHRUrun; realBuffer(:) = realMissing
     class is (gru_doubleVec);     nSpace = nGRUrun; realBuffer(:) = realMissing
      class default; message=trim(message)//'stats must be scalarv and of type gru_hru_doubleVec or gru_doubleVec'; err=20; return
    end select

    ! loop thru GRUs and HRUs
    do iGRU=1,size(gru_struc)
     do iHRU=1,gru_struc(iGRU)%hruCount

      ! get the data vectors
      select type (stat)
       class is (gru_hru_doubleVec); realBuffer(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix) = stat%gru(iGRU)%hru(iHRU)%var(map(iVar))%dat(iFreq)
       class is (gru_doubleVec); realBuffer(iGRU) = stat%gru(iGRU)%var(map(iVar))%dat(iFreq); exit ! only need to get the GRU-level data once
      end select  ! stat data structure

     end do  ! HRU loop
    end do  ! GRU loop

    ! write the data vectors
    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nSpace),start=(/1,outputTimestep(iFreq)/),count=(/nSpace,1/))
    call netcdf_err(err,message); if (err/=0) return

   ! ****************************************************************************
   ! *** write non-scalar variables (regular data structures -- instantaneous)
   ! ****************************************************************************

   ! non-scalar variables: regular data structures
   else

    ! initialize the data vectors
    select type (datt)
     class is (gru_hru_doubleVec); nSpace = nHRUrun; realArray(:,:) = realMissing;   dataType=ixReal
     class is (gru_hru_intVec);    nSpace = nHRUrun; intArray(:,:) = integerMissing; dataType=ixInteger
     class is (gru_doubleVec);     nSpace = nGRUrun; realArray(:,:) = realMissing;   dataType=ixReal
     class is (gru_intVec);        nSpace = nGRUrun; intArray(:,:) = integerMissing; dataType=ixInteger
     class default; message=trim(message)//'data is not scalarv so should be either of type gru_hru_[double or int]Vec or gru_[double or int]Vec';err=20; return
    end select

    ! loop thru GRUs and HRUs
    do iGRU=1,size(gru_struc)
     do iHRU=1,gru_struc(iGRU)%hruCount

      ! get the model layers
      nSoil   = indx%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)
      nSnow   = indx%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
      nLayers = indx%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nLayers)%dat(1)

      ! get the length of each data vector
      select case (meta(iVar)%varType)
       case (iLookVarType%wLength); datLength = nSpecBand
       case (iLookVarType%midToto); datLength = nLayers
       case (iLookVarType%midSnow); datLength = nSnow
       case (iLookVarType%midSoil); datLength = nSoil
       case (iLookVarType%ifcToto); datLength = nLayers+1
       case (iLookVarType%ifcSnow); datLength = nSnow+1
       case (iLookVarType%ifcSoil); datLength = nSoil+1
       case (iLookVarType%routing); datLength = nTimeDelay
       case default; cycle iVarLoop
       ! case parSoil only in parameters (mpar, not written here)
       ! case unknown skipped above
      end select ! varType

      ! get the data vectors
      select type (datt)
       class is (gru_hru_doubleVec); realArray(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = datt%gru(iGRU)%hru(iHRU)%var(iVar)%dat(:)
       class is (gru_hru_intVec);     intArray(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = datt%gru(iGRU)%hru(iHRU)%var(iVar)%dat(:)
       class is (gru_doubleVec); realArray(iGRU,1:datLength) = datt%gru(iGRU)%var(iVar)%dat(:); exit ! only need to get the GRU-level data once
       class is (gru_intVec);     intArray(iGRU,1:datLength) = datt%gru(iGRU)%var(iVar)%dat(:); exit ! only need to get the GRU-level data once
      end select

     end do  ! HRU loop
    end do  ! GRU loop

    ! get the maximum length of each data vector
    select case (meta(iVar)%varType)
     case (iLookVarType%wLength); maxLength = nSpecBand
     case (iLookVarType%midToto); maxLength = maxLayers
     case (iLookVarType%midSnow); maxLength = maxSnowLayers
     case (iLookVarType%midSoil); maxLength = maxSoilLayers
     case (iLookVarType%ifcToto); maxLength = maxLayers+1
     case (iLookVarType%ifcSnow); maxLength = maxSnowLayers+1
     case (iLookVarType%ifcSoil); maxLength = maxSoilLayers+1
     case (iLookVarType%routing); maxLength = nTimeDelay
     case default; cycle iVarLoop ! move onto the next variable
    end select ! varType

    ! write the data vectors
    if(maxLength==0) cycle iVarLoop ! skip if there is no length
    select case (dataType)
     case (ixReal);    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realArray(1:nSpace,1:maxLength),start=(/1,1,outputTimestep(iFreq)/),count=(/nSpace,maxLength,1/))
     case (ixInteger); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),intArray(1:nSpace,1:maxLength),start=(/1,1,outputTimestep(iFreq)/),count=(/nSpace,maxLength,1/))
    end select ! data type
    call netcdf_err(err,message); if (err/=0) return

   end if ! not scalarv

  end do iVarLoop ! iVar
 end do ! iFreq

 end subroutine writeData_perStep

 ! **************************************************************************************
 ! public subroutine writeTime: write current time to all files
 ! **************************************************************************************
 subroutine writeTime(finalizeStats,outputTimestep,meta,datt,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE globalData,only:ncid                           ! output file IDs
 USE var_lookup,only:iLookSTAT                      ! index into stat structure
 implicit none

 ! declare dummy variables
 logical(lgt)  ,intent(in)     :: finalizeStats(:)  ! flags to finalize statistics
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 type(var_info),intent(in)     :: meta(:)           ! meta data
 integer       ,intent(in)     :: datt(:)           ! timestep data
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iFreq             ! frequency index
 integer(i4b)                  :: ncVarID           ! used only for time
 ! initialize error control
 err=0;message="writeTime/"

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  do iVar = 1,size(meta)

   ! check instantaneous
   if (meta(iVar)%statIndex(iFreq)/=iLookSTAT%inst) cycle

   ! get variable id in file
   err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID)
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
   call netcdf_err(err,message)
   if (err/=0) then; err=20; return; end if

   ! add to file
   err = nf90_put_var(ncid(iFreq),ncVarID,(/datt(iVar)/),start=(/outputTimestep(iFreq)/),count=(/1/))
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
   call netcdf_err(err,message)
   if (err/=0) then; err=20; return; end if

  end do ! iVar
 end do ! iFreq

 end subroutine writeTime

 ! *********************************************************************************************************
 ! public subroutine writeRestart: print a re-start file
 ! *********************************************************************************************************
 subroutine writeRestart(filename,         & ! intent(in): name of restart file
                         nGRU_local,       & ! intent(in): number of GRUs assigned to this rank
                         nHRU_local,       & ! intent(in): number of HRUs assigned to this rank
                         prog_meta,        & ! intent(in): prognostics metadata
                         prog_data,        & ! intent(in): prognostics data
                         bvar_meta,        & ! intent(in): basin (gru) variable metadata
                         bvar_data,        & ! intent(in): basin (gru) variable data
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
 USE def_output_module,only: write_hru_info ! write HRU information to netcdf file
 
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input
 character(len=256),intent(in)      :: filename      ! name of the restart file
 integer(i4b),intent(in)            :: nGRU_local    ! number of GRUs assigned to this rank
 integer(i4b),intent(in)            :: nHRU_local    ! number of HRUs assigned to this rank
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
 ! local variables
 integer(i4b)                       :: ncid          ! netcdf file id
 integer(i4b),allocatable           :: ncVarID(:)    ! netcdf variable id
 integer(i4b)                       :: ncSnowID      ! index variable id
 integer(i4b)                       :: ncSoilID      ! index variable id
 integer(i4b)                       :: nSoil         ! number of soil layers
 integer(i4b)                       :: nSnow         ! number of snow layers
 integer(i4b)                       :: nLayers       ! number of total layers
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

 ! create file
 err = nf90_create(trim(filename),NF90_NETCDF4,ncid)
 message='iCreate[create]'; call netcdf_err(err,message); if(err/=0)return

 ! define dimensions
 err = nf90_def_dim(ncid,trim(gruDimName)    ,nGRU_local       ,    gruDimID); message='iCreate[gru]'     ; call netcdf_err(err,message); if(err/=0)return
 err = nf90_def_dim(ncid,trim(hruDimName)    ,nHRU_local       ,    hruDimID); message='iCreate[hru]'     ; call netcdf_err(err,message); if(err/=0)return
 err = nf90_def_dim(ncid,trim(tdhDimName)    ,nTimeDelay       ,    tdhDimID); message='iCreate[tdh]'     ; call netcdf_err(err,message); if(err/=0)return
 err = nf90_def_dim(ncid,trim(scalDimName)   ,nScalar          ,   scalDimID); message='iCreate[scalar]'  ; call netcdf_err(err,message); if(err/=0)return
 err = nf90_def_dim(ncid,trim(specDimName)   ,nSpecBand        ,   specDimID); message='iCreate[spectral]'; call netcdf_err(err,message); if(err/=0)return
 err = nf90_def_dim(ncid,trim(midTotoDimName),maxLayers        ,midTotoDimID); message='iCreate[midToto]' ; call netcdf_err(err,message); if(err/=0)return
 err = nf90_def_dim(ncid,trim(ifcTotoDimName),maxLayers+1      ,ifcTotoDimID); message='iCreate[ifcToto]' ; call netcdf_err(err,message); if(err/=0)return
 
 if(maxSoilLayers>0)then
   err = nf90_def_dim(ncid,trim(midSoilDimName),maxSoilLayers    ,midSoilDimID); message='iCreate[midSoil]' ; call netcdf_err(err,message); if(err/=0)return
   err = nf90_def_dim(ncid,trim(ifcSoilDimName),maxSoilLayers+1  ,ifcSoilDimID); message='iCreate[ifcSoil]' ; call netcdf_err(err,message); if(err/=0)return
 endif
 
 if(maxSnowLayers>0)then
   err = nf90_def_dim(ncid,trim(midSnowDimName),maxSnowLayers    ,midSnowDimID); message='iCreate[midSnow]' ; call netcdf_err(err,message); if(err/=0)return
   err = nf90_def_dim(ncid,trim(ifcSnowDimName),maxSnowLayers+1  ,ifcSnowDimID); message='iCreate[ifcSnow]' ; call netcdf_err(err,message); if(err/=0)return
 endif

 ! re-initialize error control
 err=0; message='writeRestart/'

 ! define prognostic variables

 ! NOTE: soil/snow prognostic variable types are only present when the
 ! corresponding maximum layer count is greater than zero.

 do iVar = 1,nProgVars
  if (prog_meta(iVar)%varType==iLookVarType%unknown) cycle

  ! define variable
  select case(prog_meta(iVar)%varType)
   case (iLookVarType%scalarv);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/hruDimID,  scalDimID /),ncVarID(iVar))
   case (iLookVarType%wLength);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/hruDimID,  specDimID /),ncVarID(iVar))
   case (iLookVarType%midToto);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/hruDimID,midTotoDimID/),ncVarID(iVar))
   case (iLookVarType%ifcToto);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/hruDimID,ifcTotoDimID/),ncVarID(iVar))
   case (iLookVarType%midSoil); if (maxSoilLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/hruDimID,midSoilDimID/),ncVarID(iVar))
   case (iLookVarType%ifcSoil); if (maxSoilLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/hruDimID,ifcSoilDimID/),ncVarID(iVar))
   case (iLookVarType%midSnow); if (maxSnowLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/hruDimID,midSnowDimID/),ncVarID(iVar))
   case (iLookVarType%ifcSnow); if (maxSnowLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/hruDimID,ifcSnowDimID/),ncVarID(iVar))
  end select

  ! check errors
  if(err/=0)then
   message=trim(message)//trim(nf90_strerror(err))//' [variable '//trim(prog_meta(iVar)%varName)//']'
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
 call netcdf_err(err,message); if(err/=0)return

 ! define attributes
 err = nf90_put_att(ncid,ncVarID(nProgVars+1),'long_name',trim(bvar_meta(iLookBVAR%routingRunoffFuture)%vardesc));   call netcdf_err(err,message); if(err/=0)return
 err = nf90_put_att(ncid,ncVarID(nProgVars+1),'units'    ,trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varunit));   call netcdf_err(err,message); if(err/=0)return

 ! define index variables - snow
 err = nf90_def_var(ncid,trim(indx_meta(iLookINDEX%nSnow)%varName),nf90_int,(/hruDimID/),ncSnowID); call netcdf_err(err,message); if(err/=0)return
 err = nf90_put_att(ncid,ncSnowID,'long_name',trim(indx_meta(iLookINDEX%nSnow)%vardesc));           call netcdf_err(err,message); if(err/=0)return
 err = nf90_put_att(ncid,ncSnowID,'units'    ,trim(indx_meta(iLookINDEX%nSnow)%varunit));           call netcdf_err(err,message); if(err/=0)return

 ! define index variables - soil
 err = nf90_def_var(ncid,trim(indx_meta(iLookINDEX%nSoil)%varName),nf90_int,(/hruDimID/),ncSoilID); call netcdf_err(err,message); if(err/=0)return
 err = nf90_put_att(ncid,ncSoilID,'long_name',trim(indx_meta(iLookINDEX%nSoil)%vardesc));           call netcdf_err(err,message); if(err/=0)return
 err = nf90_put_att(ncid,ncSoilID,'units'    ,trim(indx_meta(iLookINDEX%nSoil)%varunit));           call netcdf_err(err,message); if(err/=0)return

 ! end definition phase
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return

 ! write variables
 do iGRU = 1,nGRU_local
  do iHRU = 1,gru_struc(iGRU)%hruCount
   cHRU = gru_struc(iGRU)%hruInfo(iHRU)%hru_ix
   do iVar = 1,size(prog_meta)

    ! excape if this variable is not used
    if (prog_meta(iVar)%varType==iLookVarType%unknown) cycle

    ! actual number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
    nLayers = nSnow + nSoil

    ! check size
    ! NOTE: this may take time that we do not wish to use
    okLength=.true.
    select case (prog_meta(iVar)%varType)
     case (iLookVarType%scalarv);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nScalar  )
     case (iLookVarType%wlength);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSpecBand)
     case (iLookVarType%midSoil);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSoil    )
     case (iLookVarType%midToto);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nLayers  )
     case (iLookVarType%ifcSoil); if (nSoil>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSoil+1  )
     case (iLookVarType%ifcToto); if (nSoil>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nLayers+1)
     case (iLookVarType%midSnow); if (nSnow>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSnow    )
     case (iLookVarType%ifcSnow); if (nSnow>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSnow+1  )
     case default; err=20; message=trim(message)//'unknown var type'; return
    end select

    ! error check
    if(.not.okLength)then
     message=trim(message)//'bad vector length for variable '//trim(prog_meta(iVar)%varName)
     err=20; return
    endif

    ! write data
    select case (prog_meta(iVar)%varType)
     case (iLookVarType%scalarv);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nScalar  /))
     case (iLookVarType%wlength);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSpecBand/))
     case (iLookVarType%midSoil);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSoil    /))
     case (iLookVarType%midToto);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nLayers  /))
     case (iLookVarType%ifcSoil); if (nSoil>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSoil+1  /))
     case (iLookVarType%ifcToto); if (nSoil>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nLayers+1/))
     case (iLookVarType%midSnow); if (nSnow>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSnow    /))
     case (iLookVarType%ifcSnow); if (nSnow>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSnow+1  /))
     case default; err=20; message=trim(message)//'unknown var type'; return
    end select

    ! error check
    if (err/=0) message=trim(message)//'writing variable:'//trim(prog_meta(iVar)%varName)
    call netcdf_err(err,message); if (err/=0) return
    err=0; message='writeRestart/'

   end do ! iVar loop

   ! write index variables
   err=nf90_put_var(ncid,ncSnowID,(/indx_data%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat/),start=(/cHRU/),count=(/1/)); call netcdf_err(err,message); if (err/=0) return
   err=nf90_put_var(ncid,ncSoilID,(/indx_data%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat/),start=(/cHRU/),count=(/1/)); call netcdf_err(err,message); if (err/=0) return

  end do ! iHRU loop
  
  ! write selected basin variables
  err=nf90_put_var(ncid,ncVarID(nProgVars+1),(/bvar_data%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat/),  start=(/iGRU/),count=(/1,nTimeDelay/))
  call netcdf_err(err,message)
  if (err/=0) return 

 end do  ! iGRU loop

 ! write HRU dimension and ID for file
 call write_hru_info(ncid, gruDimID, hruDimID, err, cmessage); if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(ncVarID)

 end subroutine writeRestart

end module modelwrite_module
