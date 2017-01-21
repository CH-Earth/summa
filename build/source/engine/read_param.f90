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

module read_param_module
USE nrtype
implicit none
private
public::read_param
contains


 ! ************************************************************************************************
 ! public subroutine read_param: read trial model parameter values
 ! ************************************************************************************************
 subroutine read_param(nHRU,typeStruct,mparStruct,bparStruct,err,message)
 ! used to read model initial conditions
 USE summaFileManager,only:SETNGS_PATH               ! path for metadata files
 USE summaFileManager,only:PARAMETER_TRIAL           ! file with parameter trial values
 USE ascii_util_module,only:file_open                ! open file
 USE ascii_util_module,only:split_line               ! extract the list of variable names from the character string
 USE ascii_util_module,only:get_vlines               ! get a list of character strings from non-comment lines
 USE get_ixname_module,only:get_ixparam,get_ixbpar   ! access function to find index of elements in structure
 USE data_types,only:gru_hru_int                     ! spatial integer data type: x%hru(:)%var(:)
 USE data_types,only:gru_hru_double                  ! spatial double data type: x%hru(:)%var(:)
 USE data_types,only:gru_double                      ! spatial double data type: x%gru(:)%var(:)
 USE globalData,only:index_map                       ! mapping from global HRUs to the elements in the data structures
 USE globalData,only:gru_struc                       ! gru-hru mapping structures
 USE var_lookup,only:iLookPARAM,iLookTYPE            ! named variables to index elements of the data vectors
 implicit none
 ! define input
 integer(i4b),        intent(in)    :: nHRU             ! number of global HRUs
 type(gru_hru_int),   intent(in)    :: typeStruct       ! local classification of soil veg etc. for each HRU
 ! define output
 type(gru_hru_double),intent(inout) :: mparStruct       ! model parameters
 type(gru_double)    ,intent(inout) :: bparStruct       ! basin parameters
 integer(i4b),        intent(out)   :: err              ! error code
 character(*),        intent(out)   :: message          ! error message
 ! define local variables
 character(len=1024)                :: cmessage         ! error message for downwind routine
 character(LEN=1024)                :: infile           ! input filename
 integer(i4b)                       :: unt              ! file unit (free unit output from file_open)
 integer(i4b)                       :: iline            ! loop through lines in the file
 integer(i4b),parameter             :: maxLines=1000    ! maximum lines in the file
 integer(i4b)                       :: iend             ! check for the end of the file
 integer(i4b),parameter             :: sLen=4096        ! string length for line of parameter data
 character(LEN=sLen)                :: temp             ! single line of information
 character(LEN=sLen),allocatable    :: charline(:)      ! vector of character strings
 character(LEN=64),allocatable      :: varnames(:)      ! vector of variable names
 character(LEN=64),allocatable      :: chardata(:)      ! vector of character data
 logical(lgt)                       :: foundHRU(nHRU)   ! vector of flags to check that an HRU has been found in parameter data
 integer(i4b)                       :: hruIndex         ! HRU identifier in the file
 integer(i4b)                       :: iHRU             ! index of HRU within data vector
 integer(i4b)                       :: localHRU,iGRU    ! index of HRU and GRU within data structure
 integer(i4b)                       :: ipar,jpar        ! index of model parameter
 integer(i4b)                       :: nPars            ! number of model parameters
 integer(i4b)                       :: nDataLine        ! number of data lines in the file
 ! Start procedure here
 err=0; message="read_param/"

 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(PARAMETER_TRIAL)
 ! open file

 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! **********************************************************************************************
 ! (2) read the parameter names
 ! **********************************************************************************************
 ! loop through file until reach the first non-comment line (list of variable names)
 do iline=1,maxLines
  read(unt,'(a)',iostat=iend)temp; if(iend/=0)exit    ! read line of data
  if (temp(1:1)=='!')cycle
  ! extract the list of variable names from the character string
  call split_line(temp,varnames,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  exit
 end do
 ! save the number of parameters
 nPars = size(varnames)
 ! check that there are at least 2 "words" -- must modify at least one parameter
 if (nPars < 2) return
! if(nPars < 2)then
!  message=trim(message)//'expect need to modify at least one parameter [file = '//trim(infile)//']'
!  err=20; return
! end if
 ! check that the first parameter is the HRU index
 if(varnames(1) /= 'hruIndex' .and. varnames(1) /= 'hruId')then
  message=trim(message)//'expect first parameter name to be '//"'hruIndex' or 'hruId' [file = "//trim(infile)//']'
  err=20; return
 end if

 ! **********************************************************************************************
 ! (3) read parameter data (continue reading from previous point in the file)
 ! **********************************************************************************************
 ! get a list of character strings from non-comment lines
 call get_vlines(unt,charline,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 nDataLine = size(charline)

 ! **********************************************************************************************
 ! (4) populate the model parameter vectors
 ! **********************************************************************************************
 ! initialize the check HRU vector
 foundHRU(:) = .false. ! logical array to ensure that all HRUs are populated

 ! allocate space for the character data
 allocate(chardata(nPars),stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problem: allocate chardata"; return; end if

 ! loop through the HRUs
 dataLineLoop: do iline=1,nDataLine

  ! get the HRU index
  read(charline(iline),*,iostat=err) hruIndex
  if(err/=0)then;err=41;message=trim(message)//"problem: internal read [data='"//trim(charline(iline))//"']"; return; endif

  ! identify the HRU index
  hruLoop: do iHRU=1,nHRU
   iGRU=index_map(iHRU)%gru_ix
   localHRU=index_map(iHRU)%localHRU  
   if(hruIndex == typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex))then    
    if (foundHRU(iHRU)) then
     err=51;message=trim(message)//"duplicate HRU found in the parameter file at line: "//new_line(' ')//charline(iline); return
    else
     foundHRU(iHRU) = .true. 
     ! get the vector of parameters for a given layer, and the HRU index                                                            
     read(charline(iline),*,iostat=err) chardata
     if(err/=0)then;err=40;message=trim(message)//"problemInternalRead [data='"//trim(charline(iline))//"']"; return; endif                                                                                                                                  
     ! loop through the model parameters
     do ipar=2,nPars  ! start at #2 because the first "word" is the HRU index

      ! get the variable index
      jpar = get_ixparam(trim(varnames(ipar)))

      ! if variable does not exist in mpar struc, try bpar struct
      if(jpar<=0)then
       jpar = get_ixbpar(trim(varnames(ipar)))

       ! if not in either structure, then spit error
       if (jpar<=0)then; err=40; message=trim(message)//"cannotFindVariableIndex[name='"//trim(varnames(ipar))//"']"; return; endif

       ! populate the appropriate element of the basin par vector
       read(chardata(ipar),*,iostat=err) bparStruct%gru(iGRU)%var(jpar)
       if(err/=0)then;err=42;message=trim(message)//"problemInternalRead[data='"//trim(chardata(ipar))//"']"; return; endif
      end if

      ! populate the appropriate element of the local parameter vector
      read(chardata(ipar),*,iostat=err) mparStruct%gru(iGRU)%hru(localHRU)%var(jpar)
      if(err/=0)then;err=42;message=trim(message)//"problemInternalRead[data='"//trim(chardata(ipar))//"']"; return; endif

     end do    ! (looping through model parameters)

    end if 
    exit hruLoop
   endif
  end do hruLoop
 end do dataLineLoop    ! (looping through HRUs)

 ! check that all HRUs are populated
 if(count(foundHRU) /= nHRU)then
  do iHRU=1,nHRU
   if(.not.foundHRU(iHRU))then
    iGRU=index_map(iHRU)%gru_ix
    localHRU=index_map(iHRU)%localHRU   
    write(message,'(a,i0,a)') trim(message)//'unable to identify HRU in parameter file [index = ',&
                               typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex),'; file='//trim(infile)//']'
    err=20; return
   end if
  end do  ! looping through HRUs
 end if   ! if some HRUs are not populated

 ! **********************************************************************************************
 deallocate(varnames,charline,chardata,stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemDeallocate"; return; end if
 ! **********************************************************************************************

 ! check parameters
 do iGRU=1,size(gru_struc)
  do iHRU=1,gru_struc(iGRU)%hruCount
   call checkParam(mparStruct%gru(iGRU)%hru(iHRU),err,cmessage)
   if(err/=0)then
    write(message,'(a,i0,a,i0,a)') trim(message)//trim(cmessage)//' (iHRU = ', iHRU, '; iGRU = ', iGRU, ')'
    return
   endif
  end do  ! HRU
 end do  ! GRU

 end subroutine read_param

 ! ************************************************************************************************
 ! private subroutine checkParam: check trial model parameter values
 ! ************************************************************************************************
 subroutine checkParam(mpar_data,err,message)
 ! used to check model parameters
 USE data_types,only:var_d                      ! spatial double data type: x%hru(:)%var(:)
 USE var_lookup,only:iLookPARAM                 ! named variables to index elements of the data vectors
 implicit none
 ! dummy 
 type(var_d), intent(inout) :: mpar_data        ! model parameters
 integer(i4b),intent(out)   :: err              ! error code
 character(*),intent(out)   :: message          ! error message
 ! associations
 associate(&
 ! canopy geometry
 heightCanopyTop        => mpar_data%var(iLookPARAM%heightCanopyTop),   & ! intent(in): [dp] height at the top of the vegetation canopy (m)
 heightCanopyBottom     => mpar_data%var(iLookPARAM%heightCanopyBottom),& ! intent(in): [dp] height at the bottom of the vegetation canopy (m)
 ! transpiration
 critSoilWilting        => mpar_data%var(iLookPARAM%critSoilWilting),   & ! intent(in): [dp] critical vol. liq. water content when plants are wilting (-)
 critSoilTranspire      => mpar_data%var(iLookPARAM%critSoilTranspire), & ! intent(in): [dp] critical vol. liq. water content when transpiration is limited (-)
 ! soil properties
 theta_sat              => mpar_data%var(iLookPARAM%theta_sat),         & ! intent(in): soil porosity (-)
 theta_res              => mpar_data%var(iLookPARAM%theta_res)          & ! intent(in): soil residual volumetric water content (-)
 ) ! associations to parameters
 ! Start procedure here
 err=0; message="checkParam/"

 ! check canopy geometry
 if(heightCanopyTop < heightCanopyBottom)then
  write(message,'(a,i0,a)') trim(message)//'height of canopy top is less than the height of the canopy bottom' 
  err=20; return
 endif

 ! check transpiration
 if(critSoilTranspire < critSoilWilting)then
  write(message,'(a,i0,a)') trim(message)//'critical point for transpiration is less than the wilting point' 
  err=20; return
 endif

 ! check porosity
 if(theta_sat < theta_res)then
  write(message,'(a,i0,a)') trim(message)//'porosity is less than the residual liquid water content'
  err=20; return
 endif

 ! end associations
 end associate

 end subroutine checkParam

end module read_param_module
