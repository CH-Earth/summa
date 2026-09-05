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
USE nr_type                             ! high-level data types
USE data_types, only: cli_options       ! command-line-interface options
USE summa_type, only: summa1_type_dec   ! master summa data type

! named parameters

USE globalData,only:iRunModeFull,iRunModeGRU,iRunModeHRU

USE globalData,only:ixProgress_it,ixProgress_im,ixProgress_id,ixProgress_ih,ixProgress_never
USE globalData,only:ixRestart_iy,ixRestart_im,ixRestart_id,ixRestart_end,ixRestart_never

USE globalData,only:noNewFiles,newFileEveryOct1

! global data
USE globalData, only: iulog              ! I/O unit for logging messages
USE globalData, only: integerMissing     ! missing integer value
USE globalData, only: realMissing        ! missing double precision value

! provide access to file IDs
USE globalData,only:ncid                 ! file id of netcdf output file

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
 
 implicit none
 
 ! dummy variables
 type(summa1_type_dec),intent(inout)   :: summa1_struc        ! master summa data structure
 integer(i4b),intent(out)              :: err                 ! error code
 character(*),intent(out)              :: message             ! error message

 type(cli_options)                     :: cli_opts            ! command line interface options
 character(len=256)                    :: cmessage            ! error message of downwind routine

 err=0
 message='getCommandArguments/'

 ! parse the command-line arguments
 call parse_command_args(cli_opts, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! apply the command line arguments
 call apply_command_args(cli_opts,summa1_struc,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine getCommandArguments

 ! --------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------
 
 ! **************************************************************************************************
 ! parse the command argyments
 ! **************************************************************************************************
 subroutine parse_command_args(opts,err,message)

 ! dummy arguments
 type(cli_options), intent(out) :: opts
 integer(i4b),      intent(out) :: err
 character(*),      intent(out) :: message

 ! locals
 integer(i4b)                   :: n_arg         ! number of command line arguments
 integer(i4b)                   :: i             ! looping
 character(len=:) , allocatable :: a, v, vn      ! command line arguments
 character(len=:) , allocatable :: program_name  ! name of executable program
 character(len=256)             :: cmessage      ! error message of downwind routine

 err = 0
 message = 'parse_command_args/'

 ! name of executable program
 call get_arg(0, program_name)

 ! set defaults
 opts%suffix   = ''
 opts%new_file = noNewFiles
 opts%run_mode = iRunModeFull
 opts%progress = ixProgress_id
 opts%restart  = ixRestart_never

 ! number of command-line arguments
 n_arg = command_argument_count()
 
 ! check number of command-line arguments
 if(n_arg < 1)then
   call printCommandHelp()
   err=20; return
 endif

 ! parse command-line arguments

 i = 1
 do while (i <= n_arg)
   call get_arg(i,a)

   select case (trim(a))

     case ('--help')
       opts%show_help = .true.
       i = i + 1

     case ('-v','--version')
       opts%show_version = .true.
       i = i + 1

     case ('-m','--master')
       call require_next(i, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

       opts%master_file = trim(v)
       write(iulog,*) "master_file is '"//trim(opts%master_file)//"'."
       i = i + 2

     case ('-c','--config')
       call require_next(i, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
        opts%config_file = trim(v)
        write(iulog,*) "config_file is '"//trim(opts%config_file)//"'."
        i = i + 2

     case ('-s','--suffix')
       call require_next(i, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
       opts%suffix = trim(v)
       write(iulog,*) "file_suffix is '"//trim(opts%suffix)//"'." 
       i = i + 2

     case ('-n','--newFile')
       call require_next(i, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
       select case(trim(v))
         case ('noNewFiles');       opts%new_file = noNewFiles
         case ('newFileEveryOct1'); opts%new_file = newFileEveryOct1
         case default 
           message = trim(message)//'unknown option for new output file: expect "noNewFiles" or "newFileEveryOct1"'
           err = 1; return
       end select
       i = i + 2

     case ('-h','--hru')
     
       opts%run_mode = iRunModeHRU  
       
       call require_next(i, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

       call parse_integer(v,'iHRU',opts%hru_index,err,cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
       i = i + 2

     case ('-g','--gru')
      
       opts%run_mode = iRunModeGRU 
       
       call require_next(i, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
       call parse_integer(v,'startGRU', opts%start_gru, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

       call require_next(i+1, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

       call parse_integer(v, 'countGRU', opts%count_gru, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
       i = i + 3

     case ('-p','--progress')
       call require_next(i, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
       select case(trim(v))
         case ('t','timestep'); opts%progress = ixProgress_it
         case ('h','hour');     opts%progress = ixProgress_ih
         case ('d','day');      opts%progress = ixProgress_id
         case ('m','month');    opts%progress = ixProgress_im
         case ('n','never');    opts%progress = ixProgress_never
         case default
           message = trim(message)//'unknown frequency to print progress: "'//trim(v)//'"'
           err = 1; return
       end select
       i = i + 2
      
     case ('-r','--restart')
       call require_next(i, n_arg, a, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        
       select case(trim(v))
         case ('y','year');  opts%restart = ixRestart_iy
         case ('m','month'); opts%restart = ixRestart_im
         case ('d','day');   opts%restart = ixRestart_id
         case ('e','end');   opts%restart = ixRestart_end
         case ('n','never'); opts%restart = ixRestart_never
         case default
           message = trim(message)//'unknown frequency to write restart files: "'//trim(v)//'"'
           err = 1; return
       end select
       i = i + 2

     case ('--param')

       call require_next(i, n_arg, a, vn, err, cmessage)  ! param name
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
       call require_next(i+1, n_arg, a, v, err, cmessage)   ! param value
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
       call append_param(opts, vn, v, err, cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
       i = i + 3


     case default
        if (len_trim(a) > 0 .and. a(1:1) == '-') then
          err = 1; cmessage = 'unknown option: '//trim(a)//'; type "'//trim(program_name)//' --help" for usage'
        else
          err = 1; cmessage = "unexpected positional argument: "//trim(a)//'; type "'//trim(program_name)//' --help" for usage'
        end if
   
   end select

   ! process error code
   if(err/=0)then
    message=trim(message)//trim(cmessage)
    err=20; return
   endif

 end do  ! looping through arguments

 ! Early exits
 if (opts%show_help) then
   call printCommandHelp()
   stop 0
 end if
 if (opts%show_version) then
   call printVersionInfo(n_arg)
   stop 0
 end if

 ! validate command-line options

 if(opts%hru_index /= integerMissing .and. &
    opts%start_gru /= integerMissing)then
      message = trim(message)// &
                'single-HRU run and GRU-parallelization run cannot both be selected'
      err = 1; return
 endif

 if(opts%run_mode == iRunModeGRU)then
   if(opts%start_gru < 1 .or. opts%count_gru < 1)then
     message = trim(message)//'startGRU and countGRU must be at least 1'
     err = 1;return
    endif
  endif

 ! list parameters supplied by the CLI

 if(allocated(opts%param_name))then
   write(iulog,*) 'Parameters adjusted:'
   do i=1,size(opts%param_name)
     write(iulog,*) trim(opts%param_name(i)), opts%param_value(i)
   enddo
 endif

 end subroutine parse_command_args

 ! --------------------------------------------------------------------------------------------------
 ! Helpers
 ! --------------------------------------------------------------------------------------------------

 subroutine get_arg(i, arg)
   integer(i4b), intent(in) :: i
   character(len=:), allocatable, intent(out) :: arg
   integer(i4b)             :: L
   call get_command_argument(i, length=L)
   allocate(character(len=L) :: arg)
   call get_command_argument(i, arg)
 end subroutine get_arg

 ! --------------------------------------------------------------------------------------------------

 subroutine require_next(i, narg, opt, val, err, message)
   
   integer, intent(in) :: i, narg
   character(len=*), intent(in)               :: opt
   character(len=:), allocatable, intent(out) :: val
   integer(i4b),      intent(out) :: err 
   character(len=*),  intent(out) :: message

   character(len=:) , allocatable :: program_name  ! name of executable program
   
   err = 0
   message = 'require_next/'
   
   if (i+1 > narg) then
     call get_arg(0, program_name)
     message = trim(message)//'missing value after '//trim(opt)//'; type "'//trim(program_name)//' --help" for usage'
     err = 1; return
   end if
   call get_arg(i+1, val)
 
 end subroutine require_next

 ! --------------------------------------------------------------------------------------------------

 subroutine append_param(opts, name, value_string, err, message)

   type(cli_options), intent(inout) :: opts
   character(len=*),  intent(in)    :: name
   character(len=*),  intent(in)    :: value_string
   integer(i4b),      intent(out)   :: err
   character(len=*),  intent(out)   :: message

   real(rkind)        :: value
   character(len=256) :: cmessage

   err = 0
   message = 'append_param/'

   ! parse parameter value
   call parse_real(value_string, name, value, err, cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! append parameter name and value
   if(.not.allocated(opts%param_name))then
     opts%param_name  = [trim(name)]
     opts%param_value = [value]
   else
     opts%param_name  = [opts%param_name, trim(name)]
     opts%param_value = [opts%param_value, value]
   endif

 end subroutine append_param

 ! --------------------------------------------------------------------------------------------------

 subroutine parse_integer(value, name, result, err, message)
   character(len=*), intent(in)  :: value
   character(len=*), intent(in)  :: name
   integer(i4b),     intent(out) :: result
   integer(i4b),     intent(out) :: err
   character(len=*), intent(out) :: message

   integer :: ios

   err = 0
   message = 'parse_integer/'

   read(value,*,iostat=ios) result
   if(ios/=0)then
     message = trim(message)//'invalid '//trim(name)//' specification: "'//trim(value)//'"'
     err = 1; return
   endif

 end subroutine parse_integer

 ! --------------------------------------------------------------------------------------------------

 subroutine parse_real(value, name, result, err, message)

   character(len=*), intent(in)  :: value
   character(len=*), intent(in)  :: name
   real(rkind),      intent(out) :: result
   integer(i4b),     intent(out) :: err
   character(len=*), intent(out) :: message

   integer :: ios

   err = 0
   message = 'parse_real/'

   read(value,*,iostat=ios) result
   if(ios/=0)then
     message = trim(message)//'invalid '//trim(name)// &
               ' specification: "'//trim(value)//'"'
     err = 1
     return
   endif

 end subroutine parse_real


 ! **************************************************************************************************
 ! apply the command argyments
 ! **************************************************************************************************
 subroutine apply_command_args(opts, summa1_struc, err, message)

   ! global run controls
   USE globalData, only: iRunMode
   USE globalData, only: startGRU
   USE globalData, only: checkHRU
   USE globalData, only: newOutputFile
   USE globalData, only: output_fileSuffix
   USE globalData, only: ixProgress
   USE globalData, only: ixRestart

   ! build options
   USE build_options, only: ngen_active

   implicit none

   ! dummy variables
   type(cli_options),     intent(in)    :: opts
   type(summa1_type_dec), intent(inout) :: summa1_struc
   integer(i4b),          intent(out)   :: err
   character(*),          intent(out)   :: message

   err = 0
   message = 'apply_command_args/'

   ! *** NextGen runtime configuration
   
   if(ngen_active)then
   
     checkHRU      = integerMissing
     startGRU      = integerMissing
     newOutputFile = noNewFiles
     ixProgress    = ixProgress_never
     iRunMode      = iRunModeGRU
   
     summa1_struc%nGRU_user  = 1
     summa1_struc%nHRU_check = integerMissing
   
     return
   
   endif

   ! *** file names and output controls

   if(allocated(opts%master_file)) &
     summa1_struc%summaFileManagerFile = opts%master_file

   if(allocated(opts%config_file)) &
     summa1_struc%summaConfigFile = opts%config_file

   if(allocated(opts%suffix)) &
     output_fileSuffix = opts%suffix

   newOutputFile = opts%new_file
   ixProgress    = opts%progress
   ixRestart     = opts%restart


   ! *** run mode
   
   iRunMode = opts%run_mode
   
   select case(iRunMode)
   
     case (iRunModeFull)
   
       startGRU = 1
       checkHRU = integerMissing
   
       summa1_struc%nGRU_user  = integerMissing
       summa1_struc%nHRU_check = integerMissing
   
   
     case (iRunModeHRU)
   
       checkHRU = opts%hru_index
   
       summa1_struc%nHRU_check = 1
       summa1_struc%nGRU_user  = 1
   
       startGRU = integerMissing
   
   
     case (iRunModeGRU)
   
       startGRU = opts%start_gru
   
       summa1_struc%nGRU_user  = opts%count_gru
       summa1_struc%nHRU_check = integerMissing
   
       checkHRU = integerMissing
   
   
     case default
   
       message = trim(message)//'unknown run mode'
       err = 1
       return
   
   end select

   ! *** parameter overrides passed through the CLI

   if(allocated(opts%param_name))then
     summa1_struc%param_name  = opts%param_name
     summa1_struc%param_value = opts%param_value
   endif

   ! *** informational output

   select case(iRunMode)

     case (iRunModeHRU)
       write(iulog,'(A,I0,A)') &
         ' Single-HRU run activated. HRU ',checkHRU,' is selected for simulation.'

     case (iRunModeGRU)
       write(iulog,'(A,I0,A)') &
         ' GRU-parallelization run activated. ', summa1_struc%nGRU_user,' GRUs are selected for simulation.'

   end select

 end subroutine apply_command_args

 ! **************************************************************************************************
 ! print the SUMMA version information
 ! **************************************************************************************************
 subroutine printVersionInfo(n_arg)
 implicit none
 integer(i4b), intent(in) :: n_arg
 
 INCLUDE 'summaversion.inc' ! version information generated during compiling

 print "(A)", '----------------------------------------------------------------------'
 print "(A)", '     SUMMA - Structure for Unifying Multiple Modeling Alternatives    '
 print "(A)", repeat(' ', max(0, (70 - len('Version: ')    - len_trim(summaVersion)) / 2))//'Version: '//trim(summaVersion)
 print "(A)", repeat(' ', max(0, (70 - len('Build Time: ') - len_trim(buildTime))    / 2))//'Build Time: '//trim(buildTime)
 print "(A)", repeat(' ', max(0, (70 - len('Git Branch: ') - len_trim(gitBranch))    / 2))//'Git Branch: '//trim(gitBranch)
 print "(A)", repeat(' ', max(0, (70 - len('Git Hash: ')   - len_trim(gitHash))      / 2))//'Git Hash: '//trim(gitHash)
 print "(A)", '----------------------------------------------------------------------'
 
 if(n_arg == 1) stop 0
 
 end subroutine printVersionInfo

 ! **************************************************************************************************
 ! print the correct command line usage of SUMMA
 ! **************************************************************************************************
 subroutine printCommandHelp()
 implicit none
 
 character(len=:), allocatable :: exe
 call get_arg(0, exe)
 
 ! command line usage
 print "(//A)",'Usage: '//trim(exe)//' -m master_file [-c config_file] [-s fileSuffix] [-g startGRU countGRU] [-h iHRU] [-r freqRestart] [-p freqProgress]'
 print "(A,/)", 'Running executable: '//trim(exe)
 print "(A)",  'Running options:'
 print "(A)",  ' -m --master        Define path/name of master file (required)'
 print "(A)",  ' -c --config        Define path/name of TOML configuration file'
 print "(A)",  ' -n --newFile       Define frequency [noNewFiles,newFileEveryOct1] of new output files'
 print "(A)",  ' -s --suffix        Add fileSuffix to the output files'
 print "(A)",  ' -g --gru           Run a subset of countGRU GRUs starting from index startGRU'
 print "(A)",  ' -h --hru           Run a single HRU with index of iHRU'
 print "(A)",  ' -r --restart       Define frequency [y,m,d,e,never] to write restart files'
 print "(A)",  ' -p --progress      Define frequency [m,d,h,never] to print progress'
 print "(A)",  ' -v --version       Display version information of the current build'
 print "(A)",  ' --help             Display command-line usage'
 stop 0
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
 USE globalData,only: iulog                            ! I/O unit for logging messages

 implicit none
 ! define dummy variables
 integer(i4b),intent(in)            :: err             ! error code
 character(*),intent(in)            :: message         ! error messgage
 ! define the local variables
 integer(i4b)                       :: endModelRun(8)  ! final time
 integer(i4b)                       :: localErr        ! local error code
 integer(i4b)                       :: iFreq           ! loop through output frequencies
 real(rkind)                        :: elpSec          ! elapsed seconds

 ! close any remaining output files
 ! NOTE: use the direct NetCDF call with no error checking since the file may already be closed
 do iFreq = 1,size(ncid)
  if (ncid(iFreq)/=integerMissing) localErr = nf90_close(ncid(iFreq))
 end do
#ifndef NGEN_ACTIVE
 ! get the final date and time
 call date_and_time(values=endModelRun)
 elpSec = elapsedSec(startInit,endModelRun)

 ! print initial and final date and time
 write(iulog,"(/,A,I4,'-',I2.2,'-',I2.2,2x,I2,':',I2.2,':',I2.2,'.',I3.3)") 'initial date/time = ',startInit(1:3),  startInit(5:8)
 write(iulog,"(A,I4,'-',I2.2,'-',I2.2,2x,I2,':',I2.2,':',I2.2,'.',I3.3)")   '  final date/time = ',endModelRun(1:3),endModelRun(5:8)

 ! print elapsed time for the initialization
 write(iulog,"(/,A,1PG15.7,A)")                                             '     elapsed init = ', elapsedInit,           ' s'
 write(iulog,"(A,1PG15.7)")                                                 '    fraction init = ', elapsedInit/elpSec

 ! print elapsed time for the parameter setup
 write(iulog,"(/,A,1PG15.7,A)")                                             '    elapsed setup = ', elapsedSetup,          ' s'
 write(iulog,"(A,1PG15.7)")                                                 '   fraction setup = ', elapsedSetup/elpSec

 ! print elapsed time to read the restart data
 write(iulog,"(/,A,1PG15.7,A)")                                             '  elapsed restart = ', elapsedRestart,        ' s'
 write(iulog,"(A,1PG15.7)")                                                 ' fraction restart = ', elapsedRestart/elpSec

 ! print elapsed time for the data read
 write(iulog,"(/,A,1PG15.7,A)")                                             '     elapsed read = ', elapsedRead,           ' s'
 write(iulog,"(A,1PG15.7)")                                                 '    fraction read = ', elapsedRead/elpSec

 ! print elapsed time for the data write
 write(iulog,"(/,A,1PG15.7,A)")                                             '    elapsed write = ', elapsedWrite,          ' s'
 write(iulog,"(A,1PG15.7)")                                                 '   fraction write = ', elapsedWrite/elpSec

 ! print elapsed time for the physics
 write(iulog,"(/,A,1PG15.7,A)")                                             '  elapsed physics = ', elapsedPhysics,        ' s'
 write(iulog,"(A,1PG15.7)")                                                 ' fraction physics = ', elapsedPhysics/elpSec

 ! print total elapsed time
 write(iulog,"(/,A,1PG15.7,A)")                                             '     elapsed time = ', elpSec,                ' s'
 write(iulog,"(A,1PG15.7,A)")                                               '       or           ', elpSec/60_rkind,          ' m'
 write(iulog,"(A,1PG15.7,A)")                                               '       or           ', elpSec/3600_rkind,        ' h'
 write(iulog,"(A,1PG15.7,A/)")                                              '       or           ', elpSec/86400_rkind,       ' d'

 ! print the number of threads
 write(iulog,"(A,i10,/)")                                                   '   number threads = ', nThreads
#endif
 ! stop with message
 if(err==0)then
  write(iulog,*) 'FORTRAN STOP: '//trim(message)
  stop
 else
  write(iulog,*) 'FATAL ERROR: '//trim(message)
  stop 1
 endif

 end subroutine

end module summa_util
