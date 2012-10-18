program multi_driver
! used to evaluate different methods for simulating snow processes
! *****************************************************************************
! use desired modules
! *****************************************************************************
USE nrtype                                                  ! variable types, etc.
USE snow_fileManager,only:fuse_SetDirsUndPhiles             ! sets directories and filenames
USE snow_fileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
USE allocspace_module,only:init_metad                       ! module to allocate space for metadata structures
USE mDecisions_module,only:mDecisions                       ! module to read model decisions
USE read_metad_module,only:read_metad                       ! module to populate metadata structures
USE def_output_module,only:def_output                       ! module to define model output
USE ffile_info_module,only:ffile_info                       ! module to read information on forcing datafile
USE read_pinit_module,only:read_pinit                       ! module to read initial model parameter values
USE read_icond_module,only:read_icond                       ! module to read initial conditions
USE read_param_module,only:read_param                       ! module to read model parameter sets
USE ConvE2Temp_module,only:E2T_lookup                       ! module to calculate a look-up table for the temperature-enthalpy conversion
USE var_derive_module,only:calcHeight                       ! module to calculate height at layer interfaces and layer mid-point
USE var_derive_module,only:turbExchng                       ! module to calculate turbulaent exchange coefficients for neutral conditons
USE var_derive_module,only:v_shortcut                       ! module to calculate "short-cut" variables
USE var_derive_module,only:rootDensty                       ! module to calculate the vertical distribution of roots
USE var_derive_module,only:satHydCond                       ! module to calculate the saturated hydraulic conductivity in each soil layer
USE read_force_module,only:read_force                       ! module to read model forcing data
USE derivforce_module,only:derivforce                       ! module to compute derived forcing data
USE modelwrite_module,only:writeParam,writeForce,writeModel ! module to write model output
USE coupled_em_module,only:coupled_em                       ! module to run the coupled energy and mass model
USE data_struc,only:forcFileInfo                            ! information on forcing data file
USE data_struc,only:time_data,forc_data                     ! time and forcing data structures
USE data_struc,only:mpar_meta,mpar_data,mpar_sets           ! model parameter information
USE data_struc,only:indx_data,indx_meta                     ! index data structures
USE var_lookup,only:iLookTIME,iLookINDEX,iLookFORCE         ! identifies element of the index and forcing structures
implicit none

! *****************************************************************************
! (0) variable definitions
! *****************************************************************************
integer(i4b)              :: iParSet=0                      ! loop through parameter sets
integer(i4b)              :: nParSets=0                     ! number of parameter sets
integer(i4b)              :: iStep=0                        ! index of model time step
integer(i4b)              :: jStep=0                        ! index of model output
character(len=8)          :: cdate1=''                      ! initial date
character(len=10)         :: ctime1=''                      ! initial time
character(len=32)         :: output_fileSuffix=''           ! suffix for the output file 
character(len=256)        :: fuseFileManager=''             ! path/name of file defining directories and files
character(len=256)        :: fileout=''                     ! output filename
integer(i4b)              :: err=0                          ! error code
character(len=512)        :: message=''                     ! error message
integer(i4b),pointer      :: nSnow=>null()                  ! number of snow layers
integer(i4b),pointer      :: nSoil=>null()                  ! number of soil layers
integer(i4b),pointer      :: nLayers=>null()                ! total number of layers
integer(i4b),pointer      :: midSnowStartIndex=>null()      ! start index of the midSnow vector for a given timestep
integer(i4b),pointer      :: midSoilStartIndex=>null()      ! start index of the midSoil vector for a given timestep
integer(i4b),pointer      :: midTotoStartIndex=>null()      ! start index of the midToto vector for a given timestep
integer(i4b),pointer      :: ifcSnowStartIndex=>null()      ! start index of the ifcSnow vector for a given timestep
integer(i4b),pointer      :: ifcSoilStartIndex=>null()      ! start index of the ifcSoil vector for a given timestep
integer(i4b),pointer      :: ifcTotoStartIndex=>null()      ! start index of the ifcToto vector for a given timestep
real(dp)                  :: dt_init=0._dp                  ! used to initialize the length of the sub-step

! *****************************************************************************
! (1) inital priming -- get command line arguments, identify files, etc.
! *****************************************************************************
! get the initial time
call date_and_time(cdate1,ctime1)
print*,ctime1
! get command-line arguments for the output file suffix
call getarg(1,output_fileSuffix)
if (len_trim(output_fileSuffix) == 0) then
 print*,'1st command-line argument missing, expect text string defining the output file suffix'; stop
endif
! get command-line argument for the muster file
call getarg(2,fuseFileManager) ! path/name of file defining directories and files
if (len_trim(fuseFileManager) == 0) then
 print*,'2nd command-line argument missing, expect path/name of muster file'; stop
endif
! set directories and files -- fuseFileManager used as command-line argument
call fuse_SetDirsUndPhiles(fuseFileManager,err,message); call handle_err(err,message)
! read model decisions
call mDecisions(err,message); call handle_err(err,message)

! *****************************************************************************
! (2) read model metadata
! *****************************************************************************
! initialize model metadata structures
call init_metad(err,message); call handle_err(err,message)
! read metadata on all model variables
call read_metad(err,message); call handle_err(err,message)
! read description of model forcing datafile
call ffile_info(err,message); call handle_err(err,message)
! read default values and constraints for model parameters
call read_pinit(err,message); call handle_err(err,message) 

! *****************************************************************************
! (3) read trial model parameter values -- and allocate space for parameter structures
! *****************************************************************************
! read trial model parameter values -- and allocate space for parameter structures
call read_param(nParSets,err,message); call handle_err(err,message)

! *****************************************************************************
! (4) loop through the model parameter sets
! *****************************************************************************
do iParSet=1,nParSets

 ! assign the parameter structure to the appropriate parameter set
 mpar_data => mpar_sets(iParSet)
 ! read description of model initial conditions -- also initializes model structure components
 call read_icond(err,message); call handle_err(err,message)
 ! compute derived model variables that are pretty much constant
 call E2T_lookup(err,message); call handle_err(err,message) ! calculate a look-up table for the temperature-enthalpy conversion
 call turbExchng(err,message); call handle_err(err,message) ! calculate turbulent exchange coefficients under neutral conditions
 call rootDensty(err,message); call handle_err(err,message) ! calculate vertical distribution of root density
 call calcHeight(err,message); call handle_err(err,message) ! calculate height at layer interfaces and layer mid-point
 call satHydCond(err,message); call handle_err(err,message) ! calculate saturated hydraulic conductivity in each soil layer
 call v_shortcut(err,message); call handle_err(err,message) ! calculate "short-cut" variables such as volumetric heat capacity
 ! define the filename for model spinup
 write(fileout,'(a,i0,a,i0,a)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_spinup'//trim(output_fileSuffix)//'.nc'
 ! define the file if the first parameter set
 if(iParSet==1) then
  call def_output(fileout,err,message); call handle_err(err,message)
 endif
 ! write model parameters to the model output file
 call writeParam(fileout,iParSet,err,message); call handle_err(err,message)


 ! initialize time step length
 dt_init = 10._dp ! seconds

 ! initialize time step index
 jstep=1

 ! ****************************************************************************
 ! (5) loop through time
 ! ****************************************************************************
 do istep=1,forcFileInfo%numtim

  ! assign pointers to model layers
  nSnow   => indx_data%var(iLookINDEX%nSnow)%dat(1)
  nSoil   => indx_data%var(iLookINDEX%nSoil)%dat(1)
  nLayers => indx_data%var(iLookINDEX%nLayers)%dat(1)

  ! assign pointers to model indices
  midSnowStartIndex => indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1)
  midSoilStartIndex => indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1)
  midTotoStartIndex => indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1)
  ifcSnowStartIndex => indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1)
  ifcSoilStartIndex => indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1)
  ifcTotoStartIndex => indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1)

  ! re-compute turbulent exchange coefficients for neutral conditions (the height of the surface may change)
  call turbExchng(err,message); call handle_err(err,message)


  ! ***************************************************************************
  ! (6) read/write forcing data
  ! ***************************************************************************
  ! read a line of forcing data (if not already opened, open file, and get to the correct place)
  call read_force(istep,err,message); call handle_err(err,message)
  ! compute derived forcing variables
  call derivforce(err,message); call handle_err(err,message)

  ! *****************************************************************************
  ! (7) create a new NetCDF output file, and write parameters and forcing data
  ! *****************************************************************************
  ! check the start of a new water year
  if(time_data%var(iLookTIME%im)  ==10 .and. &   ! month = October
     time_data%var(iLookTIME%id)  ==1  .and. &   ! day = 1
     time_data%var(iLookTIME%ih)  ==1  .and. &   ! hour = 1
     time_data%var(iLookTIME%imin)==0)then       ! minute = 0
   ! define the filename
   write(fileout,'(a,i0,a,i0,a)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_',&
                                  time_data%var(iLookTIME%iyyy),'-',time_data%var(iLookTIME%iyyy)+1,&
                                  trim(output_fileSuffix)//'.nc'
   ! define the file if the first parameter set
   if(iParSet==1) then
    call def_output(fileout,err,message); call handle_err(err,message)
   endif
   ! write model parameters to the model output file
   call writeParam(fileout,iParSet,err,message); call handle_err(err,message)
    ! re-initalize the indices for midSnow, midSoil, midToto, and ifcToto
    jStep=1
    indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = 1
  endif
  ! write the forcing data to the model output file
  if(iParSet==1) call writeForce(fileout,jstep,err,message); call handle_err(err,message)
  !stop 'FORTRAN STOP: after call to writeForce'

  ! ****************************************************************************
  ! (8) run the model
  ! ****************************************************************************
  print*, time_data%var, nSnow
  ! run the model for a single parameter set and time step
  call coupled_em(dt_init,err,message); call handle_err(err,message) 
  !if(istep>8332) stop 'FORTRAN STOP: after call to coupled_em'

  ! write the model output to the NetCDF file
  call writeModel(fileout,iParSet,jstep,err,message); call handle_err(err,message)
  !if(istep>7900) call handle_err(20,'stopping on a specified step: after call to writeModel')
  
  ! increment the model indices
  midSnowStartIndex = midSnowStartIndex + nSnow
  midSoilStartIndex = midSoilStartIndex + nSoil
  midTotoStartIndex = midTotoStartIndex + nLayers
  ifcSnowStartIndex = ifcSnowStartIndex + nSnow+1 
  ifcSoilStartIndex = ifcSoilStartIndex + nSoil+1 
  ifcTotoStartIndex = ifcTotoStartIndex + nLayers+1 

  ! increment the time index
  jstep = jstep+1

 end do  ! (looping through time)
 call stop_program('end of first parameter set')
end do  ! (looping through model parameter sets)
call stop_program('finished looping through parameter sets')

contains

 subroutine handle_err(err,message)
 ! used to handle error codes
 USE data_struc,only:mvar_data            ! variable data structure
 USE var_lookup,only:iLookMVAR            ! named variables defining elements in data structure
 implicit none
 ! define dummy variables
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 ! return if A-OK
 if(err==0) return
 ! process error messages
 if (err>0) then
  write(*,'(a)') 'FATAL ERROR: '//trim(message)
 else
  write(*,'(a)') 'WARNING: '//trim(message); print*,'(can keep going, but stopping anyway)'
 endif
 ! dump variables
 print*, 'error, variable dump:'
 print*, 'dt = ', dt_init
 print*, 'istep = ', istep
 if(associated(forc_data))then
  print*, 'pptrate            = ', forc_data%var(iLookFORCE%pptrate)
  print*, 'airtemp            = ', forc_data%var(iLookFORCE%airtemp)
 endif
 if(associated(mvar_data))then
  print*, 'scalarRainPlusMelt = ', mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)
  write(*,'(a,100(f11.5,1x))') 'mLayerDepth        = ', mvar_data%var(iLookMVAR%mLayerDepth)%dat
  write(*,'(a,100(f11.5,1x))') 'mLayerTemp         = ', mvar_data%var(iLookMVAR%mLayerTemp)%dat
  write(*,'(a,100(f11.5,1x))') 'mLayerVolFracIce   = ', mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat
  write(*,'(a,100(f11.5,1x))') 'mLayerVolFracLiq   = ', mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat
  print*, 'mLayerMatricHead   = ', mvar_data%var(iLookMVAR%mLayerMatricHead)%dat
 endif
 print*,'error code = ', err
 if(associated(time_data)) print*, time_data%var, nSnow
 write(*,'(a)') trim(message)
 stop
 end subroutine handle_err

 subroutine stop_program(message)
 ! used to stop program execution
 implicit none
 ! define dummy variables
 character(*),intent(in)::message
 ! define the local variables
 integer(i4b),parameter :: outunit=6               ! write to screen
 character(len=8)       :: cdate2                  ! final date
 character(len=10)      :: ctime2                  ! final time
 ! get the final date and time
 call date_and_time(cdate2,ctime2)
 ! print initial and final date and time
 write(outunit,*) 'initial date/time = '//'ccyy='//cdate1(1:4)//' - mm='//cdate1(5:6)//' - dd='//cdate1(7:8), &
                                         ' - hh='//ctime1(1:2)//' - mi='//ctime1(3:4)//' - ss='//ctime1(5:10)
 write(outunit,*) 'final date/time   = '//'ccyy='//cdate2(1:4)//' - mm='//cdate2(5:6)//' - dd='//cdate2(7:8), &
                                         ' - hh='//ctime2(1:2)//' - mi='//ctime2(3:4)//' - ss='//ctime2(5:10)
 ! stop with message
 print*,'FORTRAN STOP: '//trim(message)
 stop
 end subroutine

end program multi_driver
