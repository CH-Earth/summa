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

module summa_util
! utilities to manage summa simulation

! data types
USE nrtype                              ! high-level data types

! global data
USE globalData,only:integerMissing      ! missing integer value
USE globalData,only:realMissing         ! missing double precision value

! provide access to file IDs
USE globalData,only:ncid               ! file id of netcdf output file

! privacy
implicit none
private

! routines to make public
public::getCommandArguments
public::stop_program
public::handle_err
contains

 ! **************************************************************************************************
 ! * obtain the command line arguments
 ! **************************************************************************************************
 subroutine getCommandArguments(summa1_struc,err,message)
 ! data types
 USE summa_type, only:summa1_type_dec                         ! master summa data type
 ! provide access to named parameters
 USE globalData,only:iRunModeFull,iRunModeGRU,iRunModeHRU
 USE globalData,only:ixProgress_it,ixProgress_im,ixProgress_id,ixProgress_ih,ixProgress_never
 USE globalData,only:ixRestart_iy,ixRestart_im,ixRestart_id,ixRestart_end,ixRestart_never
 USE globalData,only:noNewFiles,newFileEveryOct1
 ! provide access to runtime options
 USE globalData,only: startGRU          ! index of the starting GRU for parallelization run
 USE globalData,only: checkHRU          ! index of the HRU for a single HRU run
 USE globalData,only: iRunMode          ! define the current running mode
 USE globalData,only: newOutputFile     ! define option for new output file
 USE globalData,only: ixProgress        ! define frequency to write progress
 USE globalData,only: ixRestart         ! define frequency to write restart files
 USE globalData,only: output_fileSuffix ! suffix for the output file
 implicit none
 ! dummy variables
 type(summa1_type_dec),intent(inout)   :: summa1_struc        ! master summa data structure
 integer(i4b),intent(out)              :: err                 ! error code
 character(*),intent(out)              :: message             ! error message
 ! local variables
 integer(i4b)                          :: iArgument           ! index of command line argument
 integer(i4b)                          :: nArgument           ! number of command line arguments
 character(len=256),allocatable        :: argString(:)        ! string to store command line arguments
 integer(i4b)                          :: nLocalArgument      ! number of command line arguments to read for a switch
 character(len=70), parameter          :: spaces = ''         ! setting a blank string
 ! version information generated during compiling
 INCLUDE 'summaversion.inc'
 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&
  nGRU                 => summa1_struc%nGRU                ,& ! number of grouped response units
  nHRU                 => summa1_struc%nHRU                ,& ! number of global hydrologic response units
  summaFileManagerFile => summa1_struc%summaFileManagerFile & ! path/name of file defining directories and files
 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='getCommandArguments/'

 ! check number of command-line arguments
 nArgument = command_argument_count()
 if (nArgument < 1) then
  call printCommandHelp()
 end if

 ! read command line arguments
 allocate(argString(nArgument))
 do iArgument = 1,nArgument
  call get_command_argument(iArgument,argString(iArgument))

  ! print versions if needed
  if (trim(argString(iArgument)) == '-v' .or. trim(argString(iArgument)) == '--version') then
   print "(A)", '----------------------------------------------------------------------'
   print "(A)", '     SUMMA - Structure for Unifying Multiple Modeling Alternatives    '
   print "(A)", spaces(1:int(real(70 - len_trim(summaVersion) - 9) / 2))//'Version: '   //trim(summaVersion)
   print "(A)", spaces(1:int(real(70 - len_trim(buildTime) - 12) / 2))  //'Build Time: '//trim(buildTime)
   print "(A)", spaces(1:int(real(70 - len_trim(gitBranch) - 12) / 2))  //'Git Branch: '//trim(gitBranch)
   print "(A)", spaces(1:int(real(70 - len_trim(gitHash) - 10) / 2))    //'Git Hash: '  //trim(gitHash)
   print "(A)", '----------------------------------------------------------------------'
   if (nArgument == 1) stop
  end if

 end do ! reading command-line arguments

 ! initialize command line argument variables
 startGRU = integerMissing; checkHRU = integerMissing
 nGRU = integerMissing; nHRU = integerMissing
 newOutputFile = noNewFiles
 iRunMode = iRunModeFull

 ! loop through all command arguments
 nLocalArgument = 0
 do iArgument = 1,nArgument
  if (nLocalArgument>0) then; nLocalArgument = nLocalArgument -1; cycle; end if ! skip the arguments have been read
  select case (trim(argString(iArgument)))

   case ('-m', '--master')
    ! update arguments
    nLocalArgument = 1
    if (iArgument+nLocalArgument>nArgument)then
     message="missing argument file_suffix; type 'summa.exe --help' for correct usage"
     err=1; return
    endif
    ! get name of master control file
    summaFileManagerFile=trim(argString(iArgument+1))
    print "(A)", "file_master is '"//trim(summaFileManagerFile)//"'."

   ! define the formation of new output files
   case ('-n', '--newFile')
    ! check that the number of command line arguments is correct
    nLocalArgument = 1  ! expect just one argument for new output files
    if (iArgument+nLocalArgument>nArgument)then
     message="missing argument file_suffix; type 'summa.exe --help' for correct usage"
     err=1; return
    endif
    ! get the decision for the formation of new output files
    select case( trim(argString(iArgument+1)) )
     case('noNewFiles');       newOutputFile = noNewFiles
     case('newFileEveryOct1'); newOutputFile = newFileEveryOct1
     case default
      message='unknown option for new output file: expect "noNewFiles" or "newFileEveryOct1"'
      err=1; return
    end select

   case ('-s', '--suffix')
    ! define file suffix
    nLocalArgument = 1
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument) then
     message="missing argument file_suffix; type 'summa.exe --help' for correct usage"
     err=1; return
    endif
    output_fileSuffix=trim(argString(iArgument+1))
    print "(A)", "file_suffix is '"//trim(output_fileSuffix)//"'."

   case ('-h', '--hru')
    ! define a single HRU run
    if (iRunMode == iRunModeGRU)then
     message="single-HRU run and GRU-parallelization run cannot be both selected."
     err=1; return
    endif
    iRunMode=iRunModeHRU
    nLocalArgument = 1
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument) call handle_err(1,"missing argument checkHRU; type 'summa.exe --help' for correct usage")
    read(argString(iArgument+1),*) checkHRU ! read the index of the HRU for a single HRU run
    nHRU=1; nGRU=1                          ! nHRU and nGRU are both one in this case
    ! examines the checkHRU is correct
    if (checkHRU<1) then
     message="illegal iHRU specification; type 'summa.exe --help' for correct usage"
     err=1; return
    else
     print '(A)',' Single-HRU run activated. HRU '//trim(argString(iArgument+1))//' is selected for simulation.'
    end if

   case ('-g','--gru')
    ! define a GRU parallelization run; get the starting GRU and countGRU
    if (iRunMode == iRunModeHRU)then
     message="single-HRU run and GRU-parallelization run cannot be both selected."
     err=1; return
    endif
    iRunMode=iRunModeGRU
    nLocalArgument = 2
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument)then
     message="missing argument startGRU or countGRU; type 'summa.exe --help' for correct usage"
     err=1; return
    endif
    read(argString(iArgument+1),*) startGRU ! read the argument of startGRU
    read(argString(iArgument+2),*) nGRU     ! read the argument of countGRU
    if (startGRU<1 .or. nGRU<1) then
     message='startGRU and countGRU must be larger than 1.'
     err=1; return
    else
     print '(A)', ' GRU-Parallelization run activated. '//trim(argString(iArgument+2))//' GRUs are selected for simulation.'
    end if

   case ('-p', '--progress')
    ! define the frequency to print progress
    nLocalArgument = 1
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument)then
     message="missing argument freqProgress; type 'summa.exe --help' for correct usage"
     err=1; return
    endif
    select case (trim(argString(iArgument+1)))
     case ('t' , 'timestep');  ixProgress = ixProgress_it
     case ('h' , 'hour');      ixProgress = ixProgress_ih
     case ('d' , 'day');       ixProgress = ixProgress_id  ! default
     case ('m' , 'month');     ixProgress = ixProgress_im
     case ('n' , 'never');     ixProgress = ixProgress_never
     case default
      message='unknown frequency to print progress'
      err=1; return
    end select

   case ('-r', '--restart')
    ! define the frequency to write restart files
    nLocalArgument = 1
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument)then
     message="missing argument freqRestart; type 'summa.exe --help' for correct usage"
     err=1; return
    endif
    select case (trim(argString(iArgument+1)))
     case ('y' , 'year');  ixRestart = ixRestart_iy
     case ('m' , 'month'); ixRestart = ixRestart_im
     case ('d' , 'day');   ixRestart = ixRestart_id
     case ('e' , 'end');   ixRestart = ixRestart_end
     case ('n' , 'never'); ixRestart = ixRestart_never
     case default
      message='unknown frequency to write restart files'
      err=1; return
    end select

   ! do nothing
   case ('-v','--version')

   ! print help message
   case ('--help')
    call printCommandHelp

   case default
    call printCommandHelp
    message='unknown command line option'
    err=1; return

  end select
 end do  ! looping through command line arguments

 ! check if master_file has been received.
 if (len(trim(summaFileManagerFile))==0)then
  message="master_file is not received; type 'summa.exe --help' for correct usage"
  err=1; return
 endif

 ! set startGRU for full run
 if (iRunMode==iRunModeFull) startGRU=1

 ! end associate statements
 end associate summaVars

 end subroutine getCommandArguments

 ! **************************************************************************************************
 ! print the correct command line usage of SUMMA
 ! **************************************************************************************************
 subroutine printCommandHelp()
 implicit none
 ! command line usage
 print "(//A)",'Usage: summa.exe -m master_file [-s fileSuffix] [-g startGRU countGRU] [-h iHRU] [-r freqRestart] [-p freqProgress] [-c]'
 print "(A,/)",  ' summa.exe          summa executable'
 print "(A)",  'Running options:'
 print "(A)",  ' -m --master        Define path/name of master file (required)'
 print "(A)",  ' -n --newFile       Define frequency [noNewFiles,newFileEveryOct1] of new output files'
 print "(A)",  ' -s --suffix        Add fileSuffix to the output files'
 print "(A)",  ' -g --gru           Run a subset of countGRU GRUs starting from index startGRU'
 print "(A)",  ' -h --hru           Run a single HRU with index of iHRU'
 print "(A)",  ' -r --restart       Define frequency [y,m,d,e,never] to write restart files'
 print "(A)",  ' -p --progress      Define frequency [m,d,h,never] to print progress'
 print "(A)",  ' -v --version       Display version information of the current build'
 stop
 end subroutine printCommandHelp

 ! **************************************************************************************************
 ! error handler
 ! **************************************************************************************************
 subroutine handle_err(err,message)
 USE netcdf_util_module,only:nc_file_close             ! module to handle netcdf stuff for inputs and outputs
 implicit none
 ! dummy variables
 integer(i4b),intent(in)            :: err             ! error code
 character(*),intent(in)            :: message         ! error message
 ! local variables
 integer(i4b)                       :: iFreq           ! loop through output frequencies
 integer(i4b)                       :: nc_err          ! error code of nc_close
 character(len=256)                 :: cmessage        ! error message of the downwind routine
 ! ---------------------------------------------------------------------------------------
 ! return if A-OK
 if(err==0) return

 ! process error messages
 if (err>0) then
  write(*,'(//a/)') 'FATAL ERROR: '//trim(message)
 else
  write(*,'(//a/)') 'WARNING: '//trim(message); print*,'(can keep going, but stopping anyway)'
 endif

 ! close any remaining output files
 do iFreq = 1,size(ncid)
  if (ncid(iFreq)/=integerMissing) then
   call nc_file_close(ncid(iFreq),nc_err,cmessage)
   if(nc_err/=0) print*, trim(cmessage)
  end if
 end do

 stop 1
 end subroutine handle_err

 ! **************************************************************************************************
 ! stop_program: stop program execution
 ! **************************************************************************************************
 subroutine stop_program(err,message)
 ! used to stop program execution
 ! desired modules
 USE netcdf                                            ! netcdf libraries
 USE time_utils_module,only:elapsedSec                 ! calculate the elapsed time
 ! global data
 USE globalData,only: nThreads                         ! number of threads
 USE globalData,only: startInit                        ! date/time for the start of the initialization
 USE globalData,only: elapsedInit                      ! elapsed time for the initialization
 USE globalData,only: elapsedSetup                     ! elapsed time for the parameter setup
 USE globalData,only: elapsedRestart                   ! elapsed time to read the restart data
 USE globalData,only: elapsedRead                      ! elapsed time for the data read
 USE globalData,only: elapsedWrite                     ! elapsed time for the stats/write
 USE globalData,only: elapsedPhysics                   ! elapsed time for the physics
 implicit none
 ! define dummy variables
 integer(i4b),intent(in)            :: err             ! error code
 character(*),intent(in)            :: message         ! error messgage
 ! define the local variables
 integer(i4b),parameter             :: outunit=6       ! write to screen
 integer(i4b)                       :: endModelRun(8)  ! final time
 integer(i4b)                       :: localErr        ! local error code
 integer(i4b)                       :: iFreq           ! loop through output frequencies
 real(rkind)                           :: elpSec          ! elapsed seconds

 ! close any remaining output files
 ! NOTE: use the direct NetCDF call with no error checking since the file may already be closed
 do iFreq = 1,size(ncid)
  if (ncid(iFreq)/=integerMissing) localErr = nf90_close(ncid(iFreq))
 end do

 ! get the final date and time
 call date_and_time(values=endModelRun)
 elpSec = elapsedSec(startInit,endModelRun)

 ! print initial and final date and time
 write(outunit,"(/,A,I4,'-',I2.2,'-',I2.2,2x,I2,':',I2.2,':',I2.2,'.',I3.3)") 'initial date/time = ',startInit(1:3),  startInit(5:8)
 write(outunit,"(A,I4,'-',I2.2,'-',I2.2,2x,I2,':',I2.2,':',I2.2,'.',I3.3)")   '  final date/time = ',endModelRun(1:3),endModelRun(5:8)

 ! print elapsed time for the initialization
 write(outunit,"(/,A,1PG15.7,A)")                                             '     elapsed init = ', elapsedInit,           ' s'
 write(outunit,"(A,1PG15.7)")                                                 '    fraction init = ', elapsedInit/elpSec

 ! print elapsed time for the parameter setup
 write(outunit,"(/,A,1PG15.7,A)")                                             '    elapsed setup = ', elapsedSetup,          ' s'
 write(outunit,"(A,1PG15.7)")                                                 '   fraction setup = ', elapsedSetup/elpSec

 ! print elapsed time to read the restart data
 write(outunit,"(/,A,1PG15.7,A)")                                             '  elapsed restart = ', elapsedRestart,        ' s'
 write(outunit,"(A,1PG15.7)")                                                 ' fraction restart = ', elapsedRestart/elpSec

 ! print elapsed time for the data read
 write(outunit,"(/,A,1PG15.7,A)")                                             '     elapsed read = ', elapsedRead,           ' s'
 write(outunit,"(A,1PG15.7)")                                                 '    fraction read = ', elapsedRead/elpSec

 ! print elapsed time for the data write
 write(outunit,"(/,A,1PG15.7,A)")                                             '    elapsed write = ', elapsedWrite,          ' s'
 write(outunit,"(A,1PG15.7)")                                                 '   fraction write = ', elapsedWrite/elpSec

 ! print elapsed time for the physics
 write(outunit,"(/,A,1PG15.7,A)")                                             '  elapsed physics = ', elapsedPhysics,        ' s'
 write(outunit,"(A,1PG15.7)")                                                 ' fraction physics = ', elapsedPhysics/elpSec

 ! print total elapsed time
 write(outunit,"(/,A,1PG15.7,A)")                                             '     elapsed time = ', elpSec,                ' s'
 write(outunit,"(A,1PG15.7,A)")                                               '       or           ', elpSec/60_rkind,          ' m'
 write(outunit,"(A,1PG15.7,A)")                                               '       or           ', elpSec/3600_rkind,        ' h'
 write(outunit,"(A,1PG15.7,A/)")                                              '       or           ', elpSec/86400_rkind,       ' d'

 ! print the number of threads
 write(outunit,"(A,i10,/)")                                                   '   number threads = ', nThreads

 ! stop with message
 if(err==0)then
  print*,'FORTRAN STOP: '//trim(message)
  stop
 else
  print*,'FATAL ERROR: '//trim(message)
  stop 1
 endif

 end subroutine

end module summa_util
