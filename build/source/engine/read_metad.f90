module read_metad_module
USE nrtype
implicit none
private
public::read_metad
! define named variables to define the type of data structures used
integer(i4b),parameter :: ix_time =1001
integer(i4b),parameter :: ix_force=1002
integer(i4b),parameter :: ix_attr =1003
integer(i4b),parameter :: ix_type =1004
integer(i4b),parameter :: ix_param=1005
integer(i4b),parameter :: ix_mvar =1006
integer(i4b),parameter :: ix_index=1007
integer(i4b),parameter :: ix_bpar =1008
integer(i4b),parameter :: ix_bvar =1009
contains

 ! ************************************************************************************************
 ! (1) new subroutine: populate metadata structures
 ! ************************************************************************************************
 subroutine read_metad(err,message)
 ! used to populate metadata structures with metadata
 USE snow_fileManager,only:SETNGS_PATH                ! path for metadata files
 USE snow_fileManager,only:META_TIME,META_ATTR,META_TYPE,META_FORCE         ! name of metadata files (global)
 USE snow_fileManager,only:META_LOCALPARAM,META_LOCALMVAR,META_LOCALINDEX   ! name of metadata files (local column)
 USE snow_fileManager,only:META_BASINPARAM,META_BASINMVAR                   ! name of metadata files (basin-average)
 USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta  ! metadata structures
 USE data_struc,only:mpar_meta,mvar_meta,indx_meta            ! metadata structures
 USE data_struc,only:bpar_meta,bvar_meta                      ! metadata structures
 implicit none
 ! declare variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 character(len=1024)                  :: cmessage    ! error message for downstream routine
 ! initialize errors
 err=0; message="read_metad/"
 ! populate time structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_TIME),ix_time,time_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'time/'//trim(cmessage); return; endif
 ! populate local attributes structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_ATTR),ix_attr,attr_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'attr/'//trim(cmessage); return; endif
 ! populate local category structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_TYPE),ix_type,type_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'type/'//trim(cmessage); return; endif
 ! populate forcing structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_FORCE),ix_force,forc_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'forc/'//trim(cmessage); return; endif
 ! populate local parameter structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_LOCALPARAM),ix_param,mpar_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'param/'//trim(cmessage); return; endif
 ! populate local model variable structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_LOCALMVAR),ix_mvar,mvar_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'mvar/'//trim(cmessage); return; endif
 ! populate local model variable structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_LOCALINDEX),ix_index,indx_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'indx/'//trim(cmessage); return; endif
 ! populate basin parameter structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_BASINPARAM),ix_bpar,bpar_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'bpar/'//trim(cmessage); return; endif
 ! populate basin model variable structure with metadata
 call v_metadata(trim(SETNGS_PATH)//trim(META_BASINMVAR),ix_bvar,bvar_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'bvar/'//trim(cmessage); return; endif
 end subroutine read_metad


 ! ************************************************************************************************
 ! (1) new subroutine: read metadata from a file
 ! ************************************************************************************************
 subroutine v_metadata(infile,ivar_lookup,meta_vec,err,message)
 ! used to read metadata from an input file and populate the appropriate metadata structure
 USE data_struc,only:var_info                                   ! metadata structure
 USE ascii_util_module,only:file_open
 ! identify indices of named variables
 USE get_ixname_module,only:get_ixTime
 USE get_ixname_module,only:get_ixAttr
 USE get_ixname_module,only:get_ixType
 USE get_ixname_module,only:get_ixForce
 USE get_ixname_module,only:get_ixParam
 USE get_ixname_module,only:get_ixMvar
 USE get_ixname_module,only:get_ixIndex
 USE get_ixname_module,only:get_ixBpar
 USE get_ixname_module,only:get_ixBvar
 implicit none
 ! define input
 character(*),intent(in)              :: infile         ! input filename
 integer(i4b),intent(in)              :: ivar_lookup    ! used to identify variable type
 ! define input/output
 type(var_info),intent(inout),pointer :: meta_vec(:)    ! vector of metadata
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 integer(i4b),parameter               :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                         :: iline          ! loop through lines in the file 
 integer(i4b),parameter               :: maxLines=1000  ! maximum lines in the file 
 character(LEN=256)                   :: temp           ! single lime of information
 integer(i4b)                         :: iend           ! check for the end of the file
 character(LEN=256)                   :: ffmt           ! file format
 type(var_info)                       :: metaTemp       ! temporary metadata structure
 character(len=1)                     :: dLim(4)        ! column delimiter
 integer(i4b)                         :: ivar           ! index of model variable
 ! Start procedure here
 err=0; message="v_metadata/"
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get to the start of the variable descriptions 
 do iline=1,maxLines
  read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit    ! read line of data
  if (temp(1:1)/='!') exit  ! assume first line not comment is format code
 end do ! looping through file to find the format code
 ! read in format string
 read(temp,*)ffmt
 ! loop through the lines in the file
 do iline=1,maxLines
  ! read a line of data and exit iif an error code (character read, so only possible error is end of file)
  read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit
  ! check that the line is not a comment
  if (temp(1:1)=='!')cycle
  ! save data into a temporary structure
  read(temp,trim(ffmt),iostat=err) metaTemp%varname,dLim(1),metaTemp%vardesc,dLim(2),metaTemp%varunit,dLim(3),&
                                   metaTemp%vartype,dLim(4),metaTemp%v_write
  if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
  ! check that the delimiters are in the correct place
  if(any(dLim /= '|'))then
   message=trim(message)//'delimiter is not in the correct place; line = ['//trim(temp)//']; filename = '//trim(infile)
   err=32; return
  endif
  ! identify the index of the named variable
  select case(ivar_lookup)
   case(ix_time);  ivar = get_ixTime(metaTemp%varname)
   case(ix_attr);  ivar = get_ixAttr(metaTemp%varname)
   case(ix_type);  ivar = get_ixType(metaTemp%varname)
   case(ix_force); ivar = get_ixForce(metaTemp%varname)
   case(ix_param); ivar = get_ixParam(metaTemp%varname)
   case(ix_mvar);  print*, trim(metaTemp%varname); ivar = get_ixMvar(metaTemp%varname)
   case(ix_index); ivar = get_ixIndex(metaTemp%varname)
   case(ix_bpar);  ivar = get_ixBpar(metaTemp%varname)
   case(ix_bvar);  ivar = get_ixBvar(metaTemp%varname)
   case default; err=35; message=trim(message)//"caseNotFound"; return
  end select
  if(ivar<=0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(metaTemp%varname)//"]"; return; endif
  ! check if index is within range
  if(ivar>size(meta_vec))then; err=50; message=trim(message)//"variableExceedsVectorSize[var="//trim(metaTemp%varname)//"]"; return; endif
  ! put data into the metadata vector
  meta_vec(ivar) = metaTemp 
 enddo  ! looping through lines in the file
 ! check that all elements are populated
 if(any(meta_vec(:)%varname==''))then
  do iline=1,size(meta_vec)
   print*,iline,' -> ',trim(meta_vec(iline)%varname)
  end do
  err=40; message=trim(message)//"someVariablesNotPopulated"; return
 endif
 ! close file unit
 close(unt)
 end subroutine v_metadata

end module read_metad_module
