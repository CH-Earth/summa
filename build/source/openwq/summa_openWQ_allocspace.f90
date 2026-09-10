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

module allocspace_progStruct_module

  ! data types
  USE nr_type

  ! provide access to the derived type used for the openWQ start-of-timestep snapshot
  USE data_types,only:gru_hru_dom_doubleVec   ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (dp)

  ! metadata structure
  USE data_types,only:var_info                ! data type for metadata

  ! privacy
  implicit none
  private
  public::allocGlobal_progStruct

  ! -----------------------------------------------------------------------------------------------------------------------------------
  contains
 ! ************************************************************************************************
 ! public subroutine allocGlobal_progStruct: allocate space for progStruct_timestep_start
 ! Modified copy of the subroutine allocGlobal() from allocspace.f90, specialized for allocating
 ! the array progStruct_timestep_start used by the openWQ coupling.
 !
 ! It differs from allocGlobal in two ways:
 !   1) it only ever handles the gru_hru_dom_doubleVec data structure (progStruct), and
 !   2) the number of snow layers is forced to nSnow (= maxSnowLayers) for every domain so that
 !      the snapshot buffer is large enough to hold the state at the start of the timestep even
 !      after the physics adds snow layers during the step.
 !
 ! The spatial layout mirrors the "spatial domain" data structures used throughout SummaSundials:
 ! every HRU is subdivided into one or more domains (upland, glacier, wetland, ...), each with its
 ! own layer stack, hence the gru -> hru -> dom -> var -> dat nesting.
 ! ************************************************************************************************
  subroutine allocGlobal_progStruct(metaStruct,dataStruct,nSnow,err,message)
    ! NOTE: safety -- ensure only used for the openWQ progStruct snapshot
    USE globalData,only: gru_struc            ! gru-hru-dom mapping structures
    USE allocspace_module,only:allocLocal
    implicit none
    ! input
    type(var_info),intent(in)             :: metaStruct(:)  ! metadata structure
    integer(i4b),intent(in)               :: nSnow          ! forced (maximum) number of snow layers for the snapshot buffer
    ! output
    type(gru_hru_dom_doubleVec),intent(inout) :: dataStruct ! data structure
    integer(i4b),intent(out)              :: err            ! error code
    character(*),intent(out)              :: message        ! error message
    ! local variables
    integer(i4b)                          :: iHRU           ! loop index through HRUs
    integer(i4b)                          :: iGRU           ! loop index through GRUs
    integer(i4b)                          :: iDOM           ! loop index through domains
    integer(i4b)                          :: nGRU           ! number of GRUs
    character(len=256)                    :: cmessage       ! error message of the downwind routine
    ! initialize error control
    err=0; message='allocGlobal_progStruct/'

    ! get the number of GRUs
    nGRU = size(gru_struc)

    ! * allocate GRU dimension
    if(allocated(dataStruct%gru))then
      err=20; message=trim(message)//'GRU dimension was unexpectedly allocated already'; return
    end if
    allocate(dataStruct%gru(nGRU),stat=err)
    if(err/=0)then; err=20; message=trim(message)//'problem allocating GRU dimension'; return; end if

    ! * allocate HRU, DOM and local (variable) dimensions
    do iGRU=1,nGRU

      allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err)
      if(err/=0)then; err=20; message=trim(message)//'problem allocating HRU dimension'; return; end if

      do iHRU=1,gru_struc(iGRU)%hruCount

        allocate(dataStruct%gru(iGRU)%hru(iHRU)%dom(gru_struc(iGRU)%hruInfo(iHRU)%domCount),stat=err)
        if(err/=0)then; err=20; message=trim(message)//'problem allocating DOM dimension'; return; end if

        do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount

          ! get the number of lake, soil and glacier ice layers for this domain
          ! (snow layers are forced to nSnow because they vary through the timestep)
          associate(&
          nLake => gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake, & ! number of lake layers for this domain
          nSoil => gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil, & ! number of soil layers for this domain
          nGlce => gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce  ) ! number of glacier ice layers for this domain

          call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU)%dom(iDOM), &
                          nSnow=nSnow,nLake=nLake,nSoil=nSoil,nGlce=nGlce,nGlac=0,err=err,message=cmessage)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

          end associate

        end do ! loop through domains
      end do ! loop through HRUs
    end do ! loop through GRUs

  end subroutine allocGlobal_progStruct

end module allocspace_progStruct_module
