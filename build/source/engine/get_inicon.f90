module get_inicon_module
implicit none
contains

 subroutine get_inicon(err,message)
 ! used to read in initial conditions
 USE nrtype                                                ! variable types, etc.
 USE snow_fileManager,only:SETNGS_PATH,M_INITIALCOND       ! path/name of initial conditions file
 USE multiparam                                            ! use the model parameter structures
 USE multistate                                            ! state structures
 implicit none
 ! define dummy variables
 integer(i4b),intent(out)::err
 character(*),intent(out)::message
 ! define local variables
 logical(lgt)::xist                        ! .TRUE. if the file exists
 integer(i4b),parameter  :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 character(LEN=256)      :: ffmt           ! file format
 character(LEN=256)      :: temp           ! single lime of information
 integer(i4b)            :: iend           ! check for the end of the file
 integer(i4b)            :: iline          ! loop through lines in the file 
 integer(i4b),parameter  :: maxLines=1000  ! maximum lines in the file
 integer(i4b),parameter  :: maxState=100   ! maximum number of model state variables
 type(state_info)        :: stateTemp      ! temporary state structure
 character(len=2)        :: dLim           ! column delimiter
 integer(i4b),parameter  :: ipart1=0       ! look-up value for part1 of the file read
 integer(i4b),parameter  :: ipart2=1       ! look-up value for part1 of the file read
 integer(i4b)            :: ipart          ! switch between part1 and part2
 integer(i4b)            :: istate         ! loop through model state variables
 integer(i4b)            :: ilayer         ! loop through model layers
 integer(i4b)            :: numCols        ! maximum number of columns desired
 real(dp),dimension(:),allocatable :: tempIni ! temporary vector to hold initial conditions
 ! Start procedure here
 err=0; message="f-fuseGetIniCon/a-OK"
 ! check if the initial conditions file exists
 inquire(file=trim(SETNGS_PATH)//trim(M_INITIALCOND),exist=xist) ! Check for existence of masterfile
 if(.not.xist)then
   message="f-fuseGetIniCon/FileNotFound[file='"//trim(SETNGS_PATH)//trim(M_INITIALCOND)//"']"
   err=10; return
 endif
 ! open initial conditions file
 open(unt,file=trim(SETNGS_PATH)//trim(M_INITIALCOND),status="old",action="read",iostat=err)
 if(err/=0)then
   message="f-fuseGetIniCon/OpenError['"//trim(SETNGS_PATH)//trim(M_INITIALCOND)//"']"
   err=20; return
 endif
 ! initialize to read part1
 ipart=ipart1
 ! loop through the lines in the file
 do iline=1,maxLines
  ! read a line of data and exit iif an error code (character read, so only possible error is end of file)
  read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit
  ! check that the line is not a comment
  if (temp(1:1)=='!')cycle
  ! now read from the file
  select case(ipart)
   ! *************************
   ! ***** process part1 *****
   case(ipart1)
    read(temp,*,iostat=err)ffmt      ! read the format string
    if (err/=0) then; err=20; message="f-fuseGetIniCon/problemWithFormatOfFormatString"; return; endif
    do istate=1,maxState
     ! check that we do not have a comment line
     read(unt,'(a)',iostat=iend)temp; if (iend/=0 .or. temp(1:1)=="!")exit
     ! read information for a given state variable
     read(temp,trim(ffmt),iostat=err) stateTemp%sname,dLim,stateTemp%sdesc,dLim,stateTemp%sunit,dLim,stateTemp%icol
     if (err/=0) then; err=20; message="f-fuseGetIniCon/problemReadingStateMetadata"; return; endif
     ! put the state information in the structure
     select case(trim(stateTemp%sname))
      case('layerTemp'      ); StateMeta%layerTemp       = stateTemp   ! temperature of each layer (K)
      case('layerBulkDenWat'); StateMeta%layerBulkDenWat = stateTemp   ! bulk density of water in each layer -- both liquid and ice (kg m-3)
      case('layerDepth'     ); StateMeta%layerDepth      = stateTemp   ! depth of each layer (m)
      case('layerEnthalpy'  ); StateMeta%layerEnthalpy   = stateTemp   ! enthalpy of each layer (J m-2)
      case('surfaceAlbedo'  ); StateMeta%surfaceAlbedo   = stateTemp   ! albedo of the surface, soil or snow (-)
      case default
       err=10; message="f-fuseGetIniCon/stateNotFound["//trim(stateTemp%sname)//"]"; return
     end select
    end do ! (looping through state variables)
    ipart=ipart2
    cycle
   ! *************************
   ! ***** process part2 *****
   case(ipart2)
    ! declare the surface albedo
    if (associated(mState%surfaceAlbedo)) deallocate(mState%surfaceAlbedo)
    allocate(mState%surfaceAlbedo,stat=err)
    if(err/=0)then; err=20; message="f-fuseGetIniCon/problemAllocate/surfaceAlbedo"; return; endif
    ! read the number of snow layers, the number of soil layers, and the surface albedo
    read(temp,*,iostat=err) nSnow, nSoil, mState%surfaceAlbedo
    if (err/=0) then; err=20; message="f-fuseGetIniCon/problemWithStringFormat[string="//trim(temp)//"]"; return; endif
    ! get the number of state variables
    nLayer = nSoil + nSnow
    ! define the surface albedo for the soil
    if (nSnow==0) mState%surfaceAlbedo = MPARAM%soilAlbedo%VALUE
    ! deallocate arrays
    if (associated(mState%layerTemp))       deallocate(mState%layerTemp)
    if (associated(mState%layerBulkDenWat)) deallocate(mState%layerBulkDenWat)
    if (associated(mState%layerDepth))      deallocate(mState%layerDepth)
    if (associated(mState%layerEnthalpy))   deallocate(mState%layerEnthalpy)
    ! allocate arrays
    allocate(mState%layerTemp(nLayer),mState%layerBulkDenWat(nLayer),mState%layerDepth(nLayer), &
             mState%layerEnthalpy(nLayer),mState%surfaceAlbedo,stat=err)
    if(err/=0)then; err=20; message="f-fuseGetIniCon/problemAllocate/DataStructures"; return; endif
    ! get the number of columns desired, and allocate space for the temporary vector
    numCols = maxval((/stateMeta%layerTemp%icol,stateMeta%layerBulkDenWat%icol,stateMeta%layerDepth%icol,stateMeta%layerEnthalpy%icol/))
    if(allocated(tempIni)) deallocate(tempIni); allocate(tempIni(numCols),stat=err)
    if(err/=0)then; err=20; message="f-fuseGetIniCon/problemAllocate/temporaryVector"; return; endif
    ! read initial states
    do ilayer=1,nLayer
     read(unt,*,iostat=err) tempIni
     if (err/=0) then; err=20; write(message,'(a,i0,a)')"f-fuseGetIniCon/problemReadingSoilIni[layer=",ilayer,"]"; return; endif
     if (stateMeta%layerTemp%icol > 0)       mState%layerTemp(ilayer)       = tempIni(stateMeta%layerTemp%icol)
     if (stateMeta%layerBulkDenWat%icol > 0) mState%layerBulkDenWat(ilayer) = tempIni(stateMeta%layerBulkDenWat%icol)
     if (stateMeta%layerDepth%icol > 0)      mState%layerDepth(ilayer)      = tempIni(stateMeta%layerDepth%icol)
     if (stateMeta%layerEnthalpy%icol > 0)   mState%layerEnthalpy(ilayer)   = tempIni(stateMeta%layerEnthalpy%icol)
    end do   ! looping through state variables
    ! deallocate the temporary vector
    deallocate(tempIni)
   ! ********************************
   ! ***** process default case *****
   case default
    err=20; message="f-fuseGetIniCon/cannotIdentifyDesiredAction"; return 
  end select
 end do ! (looping through maxLines)
 ! close file unit
 close(unt)
 end subroutine get_inicon

end module get_inicon_module
