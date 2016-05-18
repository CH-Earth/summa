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

module modelwrite_module
USE nrtype
USE netcdf
USE netcdf_util_module,only:netcdf_err                    ! netcdf error handling function
implicit none
private
public::writeParm
public::writeData
public::writeBasin
public::writeTime
! define dimension lengths
integer(i4b),parameter      :: maxSpectral=2              ! maximum number of spectral bands
contains

 ! **********************************************************************************************************
 ! public subroutine writeParm: write model parameters
 ! **********************************************************************************************************
 subroutine writeParm(iHRU,dat,meta,err,message)
 USE data_types,only:var_info,var_d,var_i        ! metadata structure type
 USE var_lookup,only:iLookStat                   ! to index into write flag
 USE multiconst,only:integerMissing
 USE globalData,only:ncid                        ! netcdf file ids
 implicit none

 ! declare input variables
 integer(i4b)  ,intent(in)   :: iHRU             ! hydrologic response unit
 class(*)      ,intent(in)   :: dat(:)           ! local attributes
 type(var_info),intent(in)   :: meta(:)          ! metadata structure
 integer(i4b)  ,intent(out)  :: err              ! error code
 character(*)  ,intent(out)  :: message          ! error message
 ! local variables
 integer(i4b)                :: iVar             ! loop through variables
 integer(i4b)  ,parameter    :: modelTime=1      ! these particular data are only output in the timestep file

 ! initialize error control
 err=0;message="f-writeParm/"

 ! loop through local column model parameters
 do iVar = 1,size(meta)

  ! check that the variable is desired
  if (.not.meta(iVar)%statFlag(iLookStat%inst)) cycle

  ! initialize message
  message=trim(message)//trim(meta(iVar)%varName)//'/'

  ! write data
  if (iHRU.ne.integerMissing) then
   select type (dat)
    type is (integer)
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/dat(iVar)/),start=(/iHRU/),count=(/1/))
    type is (real(dp))
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/dat(iVar)/),start=(/iHRU/),count=(/1/))
    class default; err=20; message=trim(message)//'unkonwn dat type (with HRU)'; return
   endselect
  else
   select type (dat)
    type is (real(dp))
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/dat(iVar)/),start=(/1/),count=(/1/))
    class default; err=20; message=trim(message)//'unkonwn dat type (no HRU)'; return
   endselect
  endif
  call netcdf_err(err,message); if (err/=0) return

  ! re-initialize message
  message="f-writeParm/"
 end do  ! looping through local column model parameters

 end subroutine writeParm

 ! **************************************************************************************
 ! public subroutine writeData: write model time-dependent data
 ! **************************************************************************************
 subroutine writeData(modelTimestep,outputTimestep,meta,stat,dat,map,indx,iHRU,err,message)
 USE data_types,only:var_info,dlength,ilength       ! type structures for passing
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE var_lookup,only:iLookIndex                     ! index into index structure
 USE var_lookup,only:iLookStat                      ! index into stat structure
 USE globalData,only:outFreq,nFreq,ncid             ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none

 ! declare dummy variables
 type(var_info),intent(in)     :: meta(:)           ! meta data
 class(*)      ,intent(in)     :: stat(:)           ! stats data
 class(*)      ,intent(in)     :: dat(:)            ! timestep data
 type(ilength) ,intent(in)     :: indx(:)           ! index data
 integer(i4b)  ,intent(in)     :: map(:)            ! map into stats child struct
 integer(i4b)  ,intent(in)     :: iHRU              ! hydrologic response unit
 integer(i4b)  ,intent(in)     :: modelTimestep     ! model time step
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iStat             ! statistics index
 integer(i4b)                  :: iFreq             ! frequency index
 integer(i4b)                  :: ncVarID           ! used only for time
 integer(i4b)                  :: nSnow             ! number of snow layers
 integer(i4b)                  :: nSoil             ! number of soil layers
 integer(i4b)                  :: nLayers           ! total number of layers
 integer(i4b)                  :: midSnowStartIndex ! start index of the midSnow vector for a given timestep
 integer(i4b)                  :: midSoilStartIndex ! start index of the midSoil vector for a given timestep
 integer(i4b)                  :: midTotoStartIndex ! start index of the midToto vector for a given timestep
 integer(i4b)                  :: ifcSnowStartIndex ! start index of the ifcSnow vector for a given timestep
 integer(i4b)                  :: ifcSoilStartIndex ! start index of the ifcSoil vector for a given timestep
 integer(i4b)                  :: ifcTotoStartIndex ! start index of the ifcToto vector for a given timestep

 ! initialize error control
 err=0;message="writeData/"

 ! model layers
 nSoil             = indx(iLookIndex%nSoil)%dat(1)
 nSnow             = indx(iLookIndex%nSnow)%dat(1)
 nLayers           = indx(iLookIndex%nLayers)%dat(1)
 ! model indices
 midSnowStartIndex = indx(iLookIndex%midSnowStartIndex)%dat(1)
 midSoilStartIndex = indx(iLookIndex%midSoilStartIndex)%dat(1)
 midTotoStartIndex = indx(iLookIndex%midTotoStartIndex)%dat(1)
 ifcSnowStartIndex = indx(iLookIndex%ifcSnowStartIndex)%dat(1)
 ifcSoilStartIndex = indx(iLookIndex%ifcSoilStartIndex)%dat(1)
 ifcTotoStartIndex = indx(iLookIndex%ifcTotoStartIndex)%dat(1)

 do iFreq = 1,nFreq
  ! check that the timestep is desired
  if (mod(modelTimestep,outFreq(iFreq)).ne.0) cycle

   ! loop through model variables
   do iVar = 1,size(meta)

    ! handle time first
    if (meta(iVar)%varName=='time') then    
     select type(stat)
      type is (dlength)
       err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID) 
       call netcdf_err(err,message); if (err/=0) return
       err = nf90_put_var(ncid(iFreq),ncVarID,(/stat(iVar)%dat(iLookStat%inst)/),start=(/outputTimestep(iFreq)/),count=(/1,1/))
       call netcdf_err(err,message); if (err/=0) return
       cycle
     class default; err=20; message=trim(message)//'time variable must be of type dlength'; return; 
     endselect
    endif

    ! check that the variable is desired
    if (meta(iVar)%outFreq.ne.iFreq) cycle

    ! loop through output stats
    do iStat = 1,maxVarStat
     ! check that the variable is desired
     if ((.not.meta(iVar)%statFlag(iStat)).or.(trim(meta(iVar)%varName)=='unknown')) cycle

     ! stats/dats output - select data type
     if (meta(iVar)%varType==iLookVarType%scalarv) then
       select type(stat)
        type is (ilength)
         err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/stat(map(iVar))%dat(iStat)/),start=(/iHRU,outputTimestep(iFreq)/),count=(/1,1/))
        type is (dlength)
         err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/stat(map(iVar))%dat(iStat)/),start=(/iHRU,outputTimestep(iFreq)/),count=(/1,1/))
        class default; err=20; message=trim(message)//'stats must be scalarv and either ilength of dlength'; return
       endselect  ! stat 
     else
      select type (dat)
       type is (dlength)
        selectcase (meta(iVar)%varType)
         case(iLookVarType%wLength); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,1,outputTimestep(iFreq)/),count=(/1,maxSpectral,1/))
         case(iLookVarType%midToto); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midTotoStartIndex/),count=(/1,nLayers/))
         case(iLookVarType%midSnow); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midSnowStartIndex/),count=(/1,nSnow/))
         case(iLookVarType%midSoil); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midSoilStartIndex/),count=(/1,nSoil/))
         case(iLookVarType%ifcToto); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcTotoStartIndex/),count=(/1,nLayers/))
         case(iLookVarType%ifcSnow); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcSnowStartIndex/),count=(/1,nSnow/))
         case(iLookVarType%ifcSoil); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcSoilStartIndex/),count=(/1,nSoil/))
        endselect ! vartype
       type is (ilength)
        selectcase (meta(iVar)%varType)
         case(iLookVarType%wLength); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,1,outputTimestep(iFreq)/),count=(/1,maxSpectral,1/))
         case(iLookVarType%midToto); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midTotoStartIndex/),count=(/1,nLayers/))
         case(iLookVarType%midSnow); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midSnowStartIndex/),count=(/1,nSnow/))
         case(iLookVarType%midSoil); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midSoilStartIndex/),count=(/1,nSoil/))
         case(iLookVarType%ifcToto); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcTotoStartIndex/),count=(/1,nLayers/))
         case(iLookVarType%ifcSnow); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcSnowStartIndex/),count=(/1,nSnow/))
         case(iLookVarType%ifcSoil); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcSoilStartIndex/),count=(/1,nSoil/))
        endselect ! vartype
      endselect ! dat
     endif ! sacalarv

     ! process error code
     if (err.ne.0) message=trim(message)//trim(meta(iVar)%varName)//'_'//trim(get_statName(iStat))
     call netcdf_err(err,message); if (err/=0) return

    enddo ! iStat
   enddo ! iVar
  enddo ! iFreq

 end subroutine writeData

 ! **************************************************************************************
 ! public subroutine writeBasin: write basin-average variables
 ! **************************************************************************************
 subroutine writeBasin(modelTimestep,outputTimestep,meta,stat,dat,map,err,message)
 USE data_types,only:var_info,dlength,ilength       ! type structures for passing
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE var_lookup,only:iLookIndex                     ! index into index structure
 USE var_lookup,only:iLookStat                      ! index into stat structure
 USE globalData,only:outFreq,nFreq,ncid             ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none

 ! declare dummy variables
 type(var_info),intent(in)     :: meta(:)           ! meta data
 type(dlength) ,intent(in)     :: stat(:)           ! stats data
 type(dlength) ,intent(in)     :: dat(:)            ! timestep data
 integer(i4b)  ,intent(in)     :: map(:)            ! map into stats child struct
 integer(i4b)  ,intent(in)     :: modelTimestep     ! model time step
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iStat             ! statistics index
 integer(i4b)                  :: iFreq             ! frequency index
 ! initialize error control
 err=0;message="f-writeBasin/"

 do iFreq = 1,nFreq
  ! check that the timestep is desired
  if (mod(modelTimestep,outFreq(iFreq)).ne.0) cycle

   ! loop through model variables
   do iVar = 1,size(meta)

    ! check that the variable is desired
    if (meta(iVar)%outFreq.ne.iFreq) cycle

    ! loop through output stats
    do iStat = 1,maxVarStat
     ! check that the variable is desired
     if ((.not.meta(iVar)%statFlag(iStat)).or.(trim(meta(iVar)%varName)=='unknown')) cycle

     ! stats/dats output - select data type
     selectcase (meta(iVar)%varType)

      case (iLookVarType%scalarv)
       err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/stat(map(iVar))%dat(iStat)/),start=(/outputTimestep(iFreq)/),count=(/1/))

      case (iLookVarType%routing)
       if (modelTimestep==1) then
        err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(map(iVar))%dat/),start=(/1/),count=(/1000/))
       endif

      case default
       err=40; message=trim(message)//"unknownVariableType[name='"//trim(meta(iVar)%varName)//"';type='"//trim(get_varTypeName(meta(iVar)%varType))//    "']"; return
      endselect ! variable type

     ! process error code
     if (err.ne.0) message=trim(message)//trim(meta(iVar)%varName)//'_'//trim(get_statName    (iStat))
     call netcdf_err(err,message); if (err/=0) return

    enddo ! iStat
   enddo ! iVar
  enddo ! iFreq

 end subroutine writeBasin

 ! **************************************************************************************
 ! public subroutine writeTime: write current time to all files 
 ! **************************************************************************************
 subroutine writeTime(modelTimestep,outputTimestep,meta,dat,err,message)
 USE data_types,only:var_info,dlength,ilength       ! type structures for passing
 USE var_lookup,only:maxVarStat,iLookStat           ! index into stats structure
 USE globalData,only:outFreq,nFreq,ncid             ! output file information
 implicit none

 ! declare dummy variables
 type(var_info),intent(in)     :: meta(:)           ! meta data
 integer       ,intent(in)     :: dat(:)            ! timestep data
 integer(i4b)  ,intent(in)     :: modelTimestep     ! model time step
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iFreq             ! frequency index
 integer(i4b)                  :: ncVarID           ! used only for time
 ! initialize error control
 err=0;message="f-writeTime/"

 do iFreq = 1,nFreq
  ! check that the timestep is desired
  if (mod(modelTimestep,outFreq(iFreq)).ne.0) cycle

   ! loop through model variables
   do iVar = 1,size(meta)

    ! if variable is desired
    if (.not.meta(iVar)%statFlag(iLookStat%inst)) cycle

    ! get variable id in file
    err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID) 
    if (err.gt.0) message=trim(message)//trim(meta(iVar)%varName)
    call netcdf_err(err,message); if (err/=0) then; err=20; return; endif

    ! add to file
    err = nf90_put_var(ncid(iFreq),ncVarID,(/dat(iVar)/),start=(/outputTimestep(iFreq)/),count=(/1/))
    if (err.gt.0) message=trim(message)//trim(meta(iVar)%varName)
    call netcdf_err(err,message); if (err/=0) then; err=20; return; endif

   enddo ! iVar
  enddo ! iFreq

 end subroutine writeTime 

end module modelwrite_module
