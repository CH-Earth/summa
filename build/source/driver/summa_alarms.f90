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

module summa_alarms
! used to set alarms to write model output

! named variables to define new output files
USE globalData, only: noNewFiles              ! no new output files
USE globalData, only: newFileEveryOct1        ! create a new file on Oct 1 every year (start of the USA water year)

! named variables to write restart files
USE globalData, only: ixRestart_iy            ! named variable to print a re-start file once per year
USE globalData, only: ixRestart_im            ! named variable to print a re-start file once per month
USE globalData, only: ixRestart_id            ! named variable to print a re-start file once per day
USE globalData, only: ixRestart_end           ! named variable to print a re-start file at the end of a run
USE globalData, only: ixRestart_never         ! named variable to print a re-start file never

! named variables to print progress
USE globalData, only: ixProgress_im           ! named variable to print progress once per month
USE globalData, only: ixProgress_id           ! named variable to print progress once per day
USE globalData, only: ixProgress_ih           ! named variable to print progress once per hour
USE globalData, only: ixProgress_it           ! named variable to print progress once per timestep
USE globalData, only: ixProgress_never        ! named variable to print progress never

! named variable for time structures
USE var_lookup,only:iLookTIME                 ! named variables for time data structure
USE var_lookup,only:iLookFreq                 ! named variables for the frequency structure

! structure dimensions
USE var_lookup,only:maxvarFreq                ! maximum number of output files

! safety: set private unless specified otherwise
implicit none
private
public::summa_setWriteAlarms
contains

 ! used to set alarms to write model output
 subroutine summa_setWriteAlarms(oldTime, newTime, endTime,       &   ! time vectors
                                 newOutputFile, defNewOutputFile, &   ! flag to define new output file
                                 ixRestart,     printRestart,     &   ! flag to print the restart file
                                 ixProgress,    printProgress,    &   ! flag to print simulation progress
                                 resetStats,    finalizeStats,    &   ! flags to reset and finalize stats
                                 statCounter,                     &   ! statistics counter
                                 err,message)                         ! error control
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nrtype                                                  ! variable types, etc.
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables: time vectors
 integer(i4b),intent(in)               :: oldTime(:)         ! time vector from the previous time step
 integer(i4b),intent(in)               :: newTime(:)         ! time vector from the current time step
 integer(i4b),intent(in)               :: endTime(:)         ! time vector at the end of the simulation
 ! dummy variables: model decisions
 integer(i4b),intent(in)               :: newOutputFile      ! option for the new output file
 integer(i4b),intent(in)               :: ixRestart          ! option to write the restart file
 integer(i4b),intent(in)               :: ixProgress         ! option to print simulation progress
 logical(lgt),intent(in)               :: resetStats(:)      ! flags to reset statistics
 ! dummy variables: alarms
 logical(lgt),intent(out)              :: defNewOutputFile   ! flag to define new output file
 logical(lgt),intent(out)              :: printRestart       ! flag to write the restart file
 logical(lgt),intent(out)              :: printProgress      ! flag to print simulation progress
 ! dummy variables: controls on statistics output
 logical(lgt),intent(out)              :: finalizeStats(:)   ! flags to finalize statistics
 integer(i4b),intent(out)              :: statCounter(:)     ! index in model output for different output frequencies
 ! dummy variables: error control
 integer(i4b),intent(out)              :: err                ! error code
 character(*),intent(out)              :: message            ! error message
 ! ---------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                          :: iFreq              ! loop through frequencies
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_setWriteAlarms/'

 ! *****************************************************************************
 ! *** define the need to create the model output file
 ! *****************************************************************************

 ! define the need to create a new output file
 select case(newOutputFile)

  ! (don't create a new output files)
  case(noNewFiles); defNewOutputFile=.false.

  ! (check for the start of the USA water year)
  case(newFileEveryOct1)
   defNewOutputFile = (newTime(iLookTIME%im) == 10 .and. &             ! month = October
                       newTime(iLookTIME%im) /= oldTime(iLookTIME%im)) ! first timestep in October

  ! (check that we found the option)
  case default; err=20; message=trim(message)//'unable to identify the option to define new output files'; return

 end select

 ! *****************************************************************************
 ! *** define the need to create a restart file
 ! *****************************************************************************
 select case(ixRestart)
  case(ixRestart_iy);    printRestart = (newTime(iLookTIME%im) == 1 .and. newTime(iLookTIME%id) == 1 .and. &
                                         newTime(iLookTIME%ih) == 0 .and. newTime(iLookTIME%imin) == 0)
  case(ixRestart_im);    printRestart = (newTime(iLookTIME%id) == 1 .and. newTime(iLookTIME%ih) == 0 .and. &
                                         newTime(iLookTIME%imin) == 0)
  case(ixRestart_id);    printRestart = (newTime(iLookTIME%ih) == 0 .and. newTime(iLookTIME%imin) == 0)
  case(ixRestart_end);   printRestart = (newTime(iLookTIME%im)   == endTime(iLookTIME%im) .and. &
                                         newTime(iLookTIME%id)   == endTime(iLookTIME%id) .and. &
                                         newTime(iLookTIME%ih)   == endTime(iLookTIME%ih) .and. &
                                         newTime(iLookTIME%imin) == endTime(iLookTIME%imin))    ! newTime does not have a '24h', won't write ending state if end_h=24
  case(ixRestart_never); printRestart = .false.
  case default; err=20; message=trim(message)//'unable to identify option for the restart file'; return
 end select

 ! *****************************************************************************
 ! *** define the need to print progress
 ! *****************************************************************************
 select case(ixProgress)
  case(ixProgress_im);    printProgress = (newTime(iLookTIME%im) /= oldTime(iLookTIME%im))  ! start month missed
  case(ixProgress_id);    printProgress = (newTime(iLookTIME%id) /= oldTime(iLookTIME%id))  ! start day missed
  case(ixProgress_ih);    printProgress = (newTime(iLookTIME%imin) == 0)
  case(ixProgress_it);    printProgress = .true.
  case(ixProgress_never); printProgress = .false.
  case default; err=20; message=trim(message)//'unable to identify option to print progress'; return
 end select

 ! *****************************************************************************
 ! *** reset counters/flags for model statistics
 ! *****************************************************************************

 ! reset output counters/flags
 do iFreq=1,maxVarFreq  ! loop through output frequencies

   ! define the need to finalize statistics
   ! NOTE: time vector is configured so that ih=0 at the start of the day, hence day in oldTime and timeStruct%var differ
   select case(iFreq)
    case(iLookFreq%day     ); finalizeStats(iFreq)=(oldTime(iLookTIME%id  )/=newTime(iLookTIME%id  ))  ! daily aggregation
    case(iLookFreq%month   ); finalizeStats(iFreq)=(oldTime(iLookTIME%im  )/=newTime(iLookTIME%im  ))  ! monthly aggregation
    case(iLookFreq%annual  ); finalizeStats(iFreq)=(oldTime(iLookTIME%iyyy)/=newTime(iLookTIME%iyyy))  ! yearly (annual) aggregation
    case(iLookFreq%timestep); finalizeStats(iFreq)=.true.          ! timestep-level output (no temporal aggregation)
    case default; err=20; message=trim(message)//'unable to identify output frequency'; return
   end select

   ! reset ouput timestep
   if(resetStats(iFreq)) statCounter(iFreq)=1

 end do ! looping through output frequencies

 end subroutine summa_setWriteAlarms

end module summa_alarms
