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

module allocspace_progStuct_module

  ! data types
  USE nrtype
  
  ! provide access to the derived types to define the data structures
  USE data_types,only:&
                      ! final data vectors
                      dlength,             & ! var%dat
                      ilength,             & ! var%dat
                      ! no spatial dimension
                      var_i,               & ! x%var(:)            (i4b)
                      var_i8,              & ! x%var(:)            integer(8)
                      var_d,               & ! x%var(:)            (dp)
                      var_flagVec,         & ! x%var(:)%dat        (logical)
                      var_ilength,         & ! x%var(:)%dat        (i4b)
                      var_dlength,         & ! x%var(:)%dat        (dp)
                      ! gru dimension
                      gru_int,             & ! x%gru(:)%var(:)     (i4b)
                      gru_int8,            & ! x%gru(:)%var(:)     integer(8)
                      gru_double,          & ! x%gru(:)%var(:)     (dp)
                      gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                      gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                      ! gru+hru dimension
                      gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                      gru_hru_int8,        & ! x%gru(:)%hru(:)%var(:)     integer(8)
                      gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                      gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                      gru_hru_doubleVec      ! x%gru(:)%hru(:)%var(:)%dat (dp)
  
  ! metadata structure
  USE data_types,only:var_info               ! data type for metadata
  
  ! access missing values
  USE globalData,only:integerMissing         ! missing integer
  USE globalData,only:realMissing            ! missing double precision number
  
  USE globalData,only: nTimeDelay            ! number of timesteps in the time delay histogram
  USE globalData,only: nBand                 ! number of spectral bands
  
  ! access variable types
  USE var_lookup,only:iLookVarType           ! look up structure for variable typed
  USE var_lookup,only:maxvarFreq             ! allocation dimension (output frequency)
  
  ! privacy
  implicit none
  private
  public::allocGlobal_porgStruct

  
  ! -----------------------------------------------------------------------------------------------------------------------------------
  contains
 ! ************************************************************************************************
 ! public subroutine allocGlobal_progStruct: allocate space for progStruct_timestep_start
 ! Modified copy of the subroutine allocGlobal() from allocspace.f90 specificly for allocating
 ! the array progStruct_timestep_start 
 ! ************************************************************************************************
  subroutine allocGlobal_porgStruct(metaStruct,dataStruct,nSnow,err,message)
    ! NOTE: safety -- ensure only used in allocGlobal
    USE globalData,only: gru_struc     ! gru-hru mapping structures
    USE allocspace_module, only:allocLocal
    implicit none
    ! input
    type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
    integer(i4b),intent(in)         :: nSnow
    ! output
    class(*),intent(out)            :: dataStruct     ! data structure
    integer(i4b),intent(out)        :: err            ! error code
    character(*),intent(out)        :: message        ! error message
    ! local variables
    logical(lgt)                    :: check          ! .true. if structure is already allocated
    integer(i4b)                    :: iHRU           ! loop index through HRUs
    integer(i4b)                    :: iGRU           ! loop index through GRUs
    integer(i4b)                    :: nGRU           ! number of GRUs
    logical(lgt)                    :: spatial        ! spatial flag
    character(len=256)              :: cmessage       ! error message of the downwind routine
    ! initialize error control
    err=0; message='allocGlobal_porgStruct/'
    ! initialize allocation check
    check=.false.
   
    ! get the number of GRUs
    nGRU = size(gru_struc)
   
    ! * allocate GRU dimension
    select type(dataStruct)
     ! gru dimension only
     class is (gru_int);           if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     class is (gru_int8);          if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     class is (gru_intVec);        if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     class is (gru_double);        if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     class is (gru_doubleVec);     if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     ! gru+hru dimensions
     class is (gru_hru_int);       if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     class is (gru_hru_int8);      if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     class is (gru_hru_intVec);    if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     class is (gru_hru_double);    if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
     class is (gru_hru_doubleVec); if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
    end select
   
    ! check errors
    if(check) then; err=20; message=trim(message)//'GRU structure was unexpectedly allocated already'; return; end if
    if(err/=0)then; err=20; message=trim(message)//'problem allocating GRU dimension'; return; end if
   
    ! * allocate HRU dimension
    do iGRU=1,nGRU
     ! allocate the HRU dimension
     select type(dataStruct)
      class is (gru_hru_int);       if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
      class is (gru_hru_int8);      if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
      class is (gru_hru_intVec);    if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
      class is (gru_hru_double);    if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
      class is (gru_hru_doubleVec); if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
      class default  ! do nothing: It is acceptable to not be any of these specified cases
     end select
     ! check errors
     if(check) then; err=20; message=trim(message)//'HRU structure was unexpectedly allocated already'; return; end if
     if(err/=0)then; err=20; message=trim(message)//'problem allocating HRU dimension'; return; end if
    end do
   
    ! * allocate local data structures where there is a spatial dimension
    gruLoop: do iGRU=1,nGRU
   
     ! initialize the spatial flag
     spatial=.false.
   
     ! loop through HRUs
     hruLoop: do iHRU=1,gru_struc(iGRU)%hruCount
   
      ! get the number of snow and soil layers
      associate(&
      ! nSnow => gru_struc(iGRU)%hruInfo(iHRU)%nSnow, & ! number of snow layers for each HRU
      nSoil => gru_struc(iGRU)%hruInfo(iHRU)%nSoil  ) ! number of soil layers for each HRU
   
      ! allocate space for structures WITH an HRU dimension
      select type(dataStruct)
       class is (gru_hru_int);       call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
       class is (gru_hru_int8);      call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
       class is (gru_hru_intVec);    call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
       class is (gru_hru_double);    call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
       class is (gru_hru_doubleVec); call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
       class default; exit hruLoop
      end select
   
      ! error check
      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
   
      ! end association to info in data structures
      end associate
   
     end do hruLoop ! loop through HRUs
   
     ! allocate space for structures *WITHOUT* an HRU dimension
     select type(dataStruct)
      class is (gru_double);    call allocLocal(metaStruct,dataStruct%gru(iGRU),nSnow=0,nSoil=0,err=err,message=cmessage); spatial=.true.
      class is (gru_doubleVec); call allocLocal(metaStruct,dataStruct%gru(iGRU),nSnow=0,nSoil=0,err=err,message=cmessage); spatial=.true.
      class default
       if(.not.spatial) exit gruLoop  ! no need to allocate spatial dimensions if none exist for a given variable
       cycle gruLoop  ! can have an HRU dimension if we get to here
     end select
   
     ! error check
     if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
   
    end do gruLoop ! loop through GRUs
   
    ! * allocate local data structures where there is no spatial dimension
    select type(dataStruct)
     class is (var_i);         call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
     class is (var_i8);        call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
     class is (var_d);         call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
     class is (var_ilength);   call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
     class is (var_dlength);   call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
     ! check identified the data type
     class default; if(.not.spatial)then; err=20; message=trim(message)//'unable to identify derived data type'; return; end if
    end select
   
    ! error check
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
   
  end subroutine allocGlobal_porgStruct
  
end module allocspace_progStuct_module
  