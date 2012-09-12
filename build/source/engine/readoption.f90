module readoption_module
USE nrtype
implicit none
private
public::readoption
contains

 ! ************************************************************************************************
 ! new subroutine: read information on model forcing fils
 ! ************************************************************************************************
 subroutine readoption(err,message)
 ! used to read metadata on the forcing data file
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
 err=0; message="f-fuse/readoption/"
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
  model_decisions(iVar)%option   = trim(option) 
  model_decisions(iVar)%decision = trim(decision)
 end do
 end subroutine readoption

end module readoption_module
