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

module qTimeDelay_module

! data types
USE nr_type

! constants, time information
USE globalData,only:data_step   ! length of the data step (s)
USE globalData,only:realMissing ! missing real number

! look-up values for the sub-grid routing method
USE mDecisions_module,only:      &
 timeDelay,&  ! time-delay histogram
 qInstant     ! instantaneous routing

implicit none
private
public::qOverland
public::qGlacier
contains


 ! *************************************************************************************************************
 ! public subroutine qOverland: compute the time delay in runoff in a basin (places runoff in future time steps)
 ! *************************************************************************************************************
 subroutine qOverland(&
                      ! input
                      ixRouting,             & ! index for routing method
                      averageTotalRunoff,    & ! total runoff to the channel from all active components (m s-1)
                      fracFuture,            & ! fraction of runoff in future time steps (m s-1)
                      qFuture,               & ! runoff in future time steps (m s-1)
                      ! output
                      averageInstantRunoff,  & ! instantaneous runoff (m s-1)
                      averageRoutedRunoff,   & ! routed runoff (m s-1)
                      err,message)             ! error control
 implicit none
 ! input
 integer(i4b),intent(in)    :: ixRouting                 ! index for routing method
 real(rkind),intent(in)     :: averageTotalRunoff        ! total runoff to the channel from all active components (m s-1)
 real(rkind),intent(in)     :: fracFuture(:)             ! fraction of runoff in future time steps (m s-1)
 real(rkind),intent(inout)  :: qFuture(:)                ! runoff in future time steps (m s-1)
 ! output
 real(rkind),intent(out)    :: averageInstantRunoff      ! instantaneous runoff (m s-1)
 real(rkind),intent(out)    :: averageRoutedRunoff       ! routed runoff (m s-1)
 integer(i4b),intent(out)   :: err                       ! error code
 character(*),intent(out)   :: message                   ! error message
 ! internal
 integer(i4b)               :: nTDH                      ! number of points in the time-delay histogram
 integer(i4b)               :: iFuture                   ! index in time delay histogram
 ! initialize error control
 err=0; message='qOverland/'

 ! assign instantaneous runoff (m s-1)  (Note: this variable is redundant with averageTotalRunoff, could remove)
 averageInstantRunoff = averageTotalRunoff

 ! compute routed runoff (m s-1)
 select case(ixRouting)  ! (select option for sub-grid routing)
  ! ** instantaneous routing
  case(qInstant)
   averageRoutedRunoff = averageInstantRunoff

  ! ** time delay histogram
  case(timeDelay)
   ! identify number of points in the time-delay histogram
   nTDH = size(qFuture)
   ! place a fraction of runoff in future steps
   qFuture(1:nTDH) = qFuture(1:nTDH) + averageInstantRunoff*fracFuture(1:nTDH)
   ! save the routed runoff
   averageRoutedRunoff = qFuture(1)
   ! move array back
   do iFuture=2,nTDH
    qFuture(iFuture-1) = qFuture(iFuture)
   end do
   qFuture(nTDH) = 0._rkind

  ! ** error checking
  case default; err=20; message=trim(message)//'cannot find option for sub-grid routing'; return

 end select ! (select option for sub-grid routing)
 ! For open water SUMMA doesn't run any calculations
 !  the values for any output variables in the netCDF will stay at the value at which they were initialized, which may be a large negative
 ! Coast may be similarly large and negative
 !if (averageRoutedRunoff < 0._rkind) averageRoutedRunoff = realMissing

 end subroutine qOverland


 ! *************************************************************************************************************
 ! public subroutine qGlacier: compute the time delay in glacier melt to route to stream in a basin (places runoff in future time steps)
 ! *************************************************************************************************************
 subroutine qGlacier(&
                     ! input
                     glacStor_kIce,         & ! storage coefficient ice reservoir (hours)
                     glacStor_kSnow,        & ! storage coefficient snow reservoir (hours)
                     glacStor_kFirn,        & ! storage coefficient firn reservoir (hours)
                     glacIceMelt,           & ! total melt into ice reservoir (m s-1)
                     glacSnowMelt,          & ! total melt into snow reservoir (m s-1)
                     glacFirnMelt,          & ! total melt into firn reservoir (m s-1)
                     glacierAblArea,        & ! per glacier ablation area (m2)
                     glacierAccArea,        & ! per glacier acumulation area (m2)   
                     nGlac,                 & ! number of glaciers in the basin
                     ! output        
                     qIceFuture,            & ! per glacier ice reservoir runoff in future time steps (m3 s-1)
                     qSnowFuture,           & ! per glacier snow reservoir runoff in future time steps (m3 s-1)
                     qFirnFuture,           & ! per glacier firn reservoir runoff in future time steps (m3 s-1)
                     glacierRoutedRunoff,   & ! routed glacier runoff (m s-1)
                     err,message)             ! error control
 implicit none
 ! input
 real(rkind),intent(in)     :: glacStor_kIce             ! storage coefficient ice reservoir (hours)
 real(rkind),intent(in)     :: glacStor_kSnow            ! storage coefficient snow reservoir (hours)
 real(rkind),intent(in)     :: glacStor_kFirn            ! storage coefficient firn reservoir (hours)
 real(rkind),intent(in)     :: glacIceMelt               ! total melt into ice reservoir (m s-1)
 real(rkind),intent(in)     :: glacSnowMelt              ! total melt into snow reservoir (m s-1)
 real(rkind),intent(in)     :: glacFirnMelt              ! total melt into firn reservoir (m s-1)
 real(rkind),intent(in)     :: glacierAblArea(:)         ! per glacier ablation area (m2)
 real(rkind),intent(in)     :: glacierAccArea(:)         ! per glacier acumulation area (m2)
integer(i4b),intent(in)     :: nGlac                     ! number of glaciers in the basin
 real(rkind),intent(inout)  :: qIceFuture(:)             ! per glacier ice reservoir runoff in future time steps (m s-1)
 real(rkind),intent(inout)  :: qSnowFuture(:)            ! per glacier snow reservoir runoff in future time steps (m s-1)
 real(rkind),intent(inout)  :: qFirnFuture(:)            ! per glacier firn reservoir runoff in future time steps (m s-1)
 ! output
 real(rkind),intent(out)    :: glacierRoutedRunoff       ! routed glacier runoff (m s-1)
 integer(i4b),intent(out)   :: err                       ! error code
 character(*),intent(out)   :: message                   ! error message
 ! internal
 real(rkind)                :: qIce                      ! hourly ice reservoir runoff (m s-1)
 real(rkind)                :: qSnow                     ! hourly snow reservoir runoff (m s-1) 
 real(rkind)                :: qFirn                     ! hourly firn reservoir runoff (m s-1)
 real(rkind)                :: frac                      ! fraction of glacier area
 real(rkind)                :: glacierAblTotal           ! total ablation area (m2)
 real(rkind)                :: glacierAccTotal           ! total accumulation area (m2)
 integer(i4b)               :: iGlacier                  ! index for glaciers
 ! initialize error control
 err=0; message='qGlacier/' 

 glacierRoutedRunoff = 0._rkind

 glacierAblTotal = sum(glacierAblArea(1:nGlac))
 glacierAccTotal = sum(glacierAccArea(1:nGlac))

 do iGlacier=1,nGlac
   ! ice reservoir runoff (m s-1)
   frac = 0._rkind
   if (glacierAblTotal>0._rkind) frac = glacierAblArea(iGlacier)/glacierAblTotal
   qIce = qIceFuture(iGlacier) + glacIceMelt*frac - glacIceMelt*frac*exp(-data_step/glacStor_kIce)
   qIceFuture(iGlacier) = qIce*exp(-data_step/glacStor_kIce) ! place runoff in future time steps 

   ! snow reservoir runoff (m s-1)
   qSnow = qSnowFuture(iGlacier) + glacSnowMelt*frac - glacSnowMelt*frac*exp(-data_step/glacStor_kSnow)
   qSnowFuture(iGlacier) = qSnow*exp(-data_step/glacStor_kSnow) ! place runoff in future time steps

   ! firn reservoir runoff (m s-1)
   frac = 0._rkind ! reset fraction
   if (glacierAccTotal>0._rkind) frac = glacierAccArea(iGlacier)/glacierAccTotal
   qFirn = qFirnFuture(iGlacier) + glacFirnMelt*frac - glacFirnMelt*frac*exp(-data_step/glacStor_kFirn)
   qFirnFuture(iGlacier) = qFirn*exp(-data_step/glacStor_kFirn) ! place runoff in future time steps 

   ! routed glacier runoff (m s-1)
   glacierRoutedRunoff = glacierRoutedRunoff + qIce + qSnow + qFirn 
 end do

 end subroutine qGlacier


end module qTimeDelay_module
