! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
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

module allocspace_module
! data types
USE nrtype
! metadata structures
USE globalData,only:time_meta,forc_meta,attr_meta,type_meta        ! metadata structures
USE globalData,only:prog_meta,diag_meta,flux_meta,deriv_meta       ! metadata structures
USE globalData,only:mpar_meta,indx_meta                            ! metadata structures
USE globalData,only:bpar_meta,bvar_meta                            ! metadata structures
USE globalData,only:model_decisions                                ! model decision structure
! time structures
USE globalData,only:refTime,startTime,finshTime                    ! reference time, start time, and end time for the model simulation
! structure information
USE globalData,only:structInfo             ! information on the data structures                  
! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,               & ! x%var(:)            (i4b)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    spatial_int,         & ! x%hru(:)%var(:)     (i4b)
                    spatial_double,      & ! x%hru(:)%var(:)     (dp)
                    spatial_intVec,      & ! x%hru(:)%var(:)%dat (i4b)
                    spatial_doubleVec      ! x%hru(:)%var(:)%dat (dp)
! metadata structure
USE data_types,only:var_info
implicit none
private
public::initStruct
public::allocLocal
! define missing values
integer(i4b),parameter :: missingInteger=-9999
real(dp),parameter     :: missingDouble=-9999._dp
! define fixed dimensions
integer(i4b),parameter :: nBand=2         ! number of spectral bands
integer(i4b),parameter :: nTimeDelay=2000 ! number of elements in the time delay histogram
! -----------------------------------------------------------------------------------------------------------------------------------
contains

 ! ************************************************************************************************
 ! public subroutine initStruct: initialize data structures
 ! ************************************************************************************************
 subroutine initStruct(&
                       ! input: model control
                       nHRU,       &    ! number of HRUs
                       nSnow,      &    ! number of snow layers for each HRU
                       nSoil,      &    ! number of soil layers for each HRU
                       ! input: data structures
                       timeStruct, &    ! model time data
                       forcStruct, &    ! model forcing data
                       attrStruct, &    ! local attributes for each HRU
                       typeStruct, &    ! local classification of soil veg etc. for each HRU
                       mparStruct, &    ! model parameters
                       indxStruct, &    ! model indices
                       progStruct, &    ! model prognostic (state) variables
                       diagStruct, &    ! model diagnostic variables
                       fluxStruct, &    ! model fluxes
                       derivStruct,&    ! model derivatives
                       bparStruct, &    ! basin-average parameters
                       bvarStruct, &    ! basin-average variables
                       ! output: error control
                       err,message)
 implicit none
 ! input: model control
 integer(i4b),            intent(in)   :: nHRU          ! number of HRUs
 integer(i4b),            intent(in)   :: nSnow(:)      ! number of snow layers for each HRU
 integer(i4b),            intent(in)   :: nSoil(:)      ! number of soil layers for each HRU
 ! output: data structures
 type(var_i),             intent(out)  :: timeStruct    ! model time data
 type(spatial_double),    intent(out)  :: forcStruct    ! model forcing data
 type(spatial_double),    intent(out)  :: attrStruct    ! local attributes for each HRU
 type(spatial_int),       intent(out)  :: typeStruct    ! local classification of soil veg etc. for each HRU
 type(spatial_double),    intent(out)  :: mparStruct    ! model parameters
 type(spatial_intVec),    intent(out)  :: indxStruct    ! model indices
 type(spatial_doubleVec), intent(out)  :: progStruct    ! model prognostic (state) variables
 type(spatial_doubleVec), intent(out)  :: diagStruct    ! model diagnostic variables
 type(spatial_doubleVec), intent(out)  :: fluxStruct    ! model fluxes
 type(var_dlength),       intent(out)  :: derivStruct   ! model derivatives
 type(var_d),             intent(out)  :: bparStruct    ! basin-average parameters
 type(var_dlength),       intent(out)  :: bvarStruct    ! basin-average variables
 ! output: error control
 integer(i4b),intent(out)              :: err           ! error code
 character(*),intent(out)              :: message       ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                          :: iStruct       ! index of data structure 
 character(len=256)                    :: cmessage      ! error message of downwind routine
 ! initialize errors
 err=0; message="initStruct/"

 ! loop through data structures
 do iStruct=1,size(structInfo)
  ! allocate space
  select case(trim(structInfo(iStruct)%structName))
   case('time'); call allocGlobal(nHRU, time_meta,  timeStruct,  err, cmessage)   ! model forcing data
   case('forc'); call allocGlobal(nHRU, forc_meta,  forcStruct,  err, cmessage)   ! model forcing data
   case('attr'); call allocGlobal(nHRU, attr_meta,  attrStruct,  err, cmessage)   ! local attributes for each HRU
   case('type'); call allocGlobal(nHRU, type_meta,  typeStruct,  err, cmessage)   ! local classification of soil veg etc. for each HRU  
   case('mpar'); call allocGlobal(nHRU, mpar_meta,  mparStruct,  err, cmessage)   ! model parameters
   case('indx'); call allocGlobal(nHRU, indx_meta,  indxStruct,  err, cmessage)   ! model variables
   case('prog'); call allocGlobal(nHRU, prog_meta,  progStruct,  err, cmessage)   ! model prognostic (state) variables
   case('diag'); call allocGlobal(nHRU, diag_meta,  diagStruct,  err, cmessage)   ! model diagnostic variables
   case('flux'); call allocGlobal(nHRU, flux_meta,  fluxStruct,  err, cmessage)   ! model fluxes
   case('deriv');call allocGlobal(nHRU, deriv_meta, derivStruct, err, cmessage)   ! model derivatives
   case('bpar'); call allocGlobal(nHRU, bpar_meta,  bparStruct,  err, cmessage)   ! basin-average parameters
   case('bvar'); call allocGlobal(nHRU, bvar_meta,  bvarStruct,  err, cmessage)   ! basin-average variables
   case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
  end select
  ! check errors
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage)//'[structure =  '//trim(structInfo(iStruct)%structName)//']'; return; endif
 end do  ! looping through data structures

 contains

  ! ************************************************************************************************
  ! internal subroutine allocGlobal: allocate space for global data structures 
  ! ************************************************************************************************
  subroutine allocGlobal(nHRU,metaStruct,dataStruct,err,message)
  implicit none
  ! input
  integer(i4b),intent(in)         :: nHRU           ! number of HRUs
  type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
  ! output
  class(*),intent(out)            :: dataStruct     ! data structure
  integer(i4b),intent(out)        :: err            ! error code
  character(*),intent(out)        :: message        ! error message
  ! local variables
  logical(lgt)                    :: check          ! .true. if structure is already allocated
  integer(i4b)                    :: iHRU           ! loop htrough HRUs
  character(len=256)              :: cmessage       ! error message of the downwind routine
  ! initialize error control
  err=0; message='allocGlobal/'
  
  ! initialize allocation check
  check=.false.
  
  ! allocate HRU dimension
  select type(dataStruct)
   type is (spatial_int);       if(allocated(dataStruct%hru))then; check=.true.; else; allocate(dataStruct%hru(nHRU),stat=err); endif 
   type is (spatial_intVec);    if(allocated(dataStruct%hru))then; check=.true.; else; allocate(dataStruct%hru(nHRU),stat=err); endif
   type is (spatial_double);    if(allocated(dataStruct%hru))then; check=.true.; else; allocate(dataStruct%hru(nHRU),stat=err); endif
   type is (spatial_doubleVec); if(allocated(dataStruct%hru))then; check=.true.; else; allocate(dataStruct%hru(nHRU),stat=err); endif
   class default  ! do nothing: It is acceptable to not be any of these specified cases
  end select

  ! check errors
  if(check) then; err=20; message=trim(message)//'structure was unexpectedly allocated already'; return; endif
  if(err/=0)then; err=20; message=trim(message)//'problem allocating'; return; endif

  ! allocate local data structures
  select type(dataStruct)
   ! structures with an HRU dimension
   type is (spatial_int);       do iHRU=1,nHRU; call allocLocal(metaStruct,dataStruct%hru(iHRU),nSnow(iHRU),nSoil(iHRU),err=err,message=cmessage); end do
   type is (spatial_intVec);    do iHRU=1,nHRU; call allocLocal(metaStruct,dataStruct%hru(iHRU),nSnow(iHRU),nSoil(iHRU),err=err,message=cmessage); end do
   type is (spatial_double);    do iHRU=1,nHRU; call allocLocal(metaStruct,dataStruct%hru(iHRU),nSnow(iHRU),nSoil(iHRU),err=err,message=cmessage); end do
   type is (spatial_doubleVec); do iHRU=1,nHRU; call allocLocal(metaStruct,dataStruct%hru(iHRU),nSnow(iHRU),nSoil(iHRU),err=err,message=cmessage); end do
   ! structures without an HRU dimension
   type is (var_i);        call allocLocal(metaStruct,dataStruct,err=err,message=cmessage) 
   type is (var_d);        call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
   type is (var_ilength);  call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
   type is (var_dlength);  call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
   ! check identified the data type
   class default; err=20; message=trim(message)//'unable to identify derived data type for the variable dimension'; return
  end select
  ! error check
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  end subroutine allocGlobal

 end subroutine initStruct

 ! ************************************************************************************************
 ! public subroutine allocLocal: allocate space for local data structures
 ! ************************************************************************************************
 subroutine allocLocal(metaStruct,dataStruct,nSnow,nSoil,err,message)
 implicit none
 ! input-output
 type(var_info),intent(in)        :: metaStruct(:)  ! metadata structure
 class(*),intent(inout)           :: dataStruct     ! data structure
 ! optional input
 integer(i4b),intent(in),optional :: nSnow          ! number of snow layers
 integer(i4b),intent(in),optional :: nSoil          ! number of soil layers
 ! output
 integer(i4b),intent(out)         :: err            ! error code
 character(*),intent(out)         :: message        ! error message
 ! local
 integer(i4b)                     :: iVar           ! loop through variables in the metadata structure
 logical(lgt)                     :: check          ! .true. if the variables are allocated
 integer(i4b)                     :: nVars          ! number of variables in the metadata structure
 integer(i4b)                     :: nLayers        ! total number of layers
 character(len=256)               :: cmessage       ! error message of the downwind routine
 ! initialize error control
 err=0; message='allocLocal/'

 ! get the number of variables in the metadata structure
 nVars = size(metaStruct)

 ! check if nSnow and nSoil are present
 if(present(nSnow) .or. present(nSoil))then
  ! check both are present
  if(.not.present(nSoil))then; err=20; message=trim(message)//'expect nSoil to be present when nSnow is present'; return; endif
  if(.not.present(nSnow))then; err=20; message=trim(message)//'expect nSnow to be present when nSoil is present'; return; endif
  nLayers = nSnow+nSoil

 ! It is possible that nSnow and nSoil are actually needed here, so we return an error if the optional arguments are missing when needed
 else
  do iVar=1,size(metaStruct)  ! loop through variables in the metadata structure
   select case(trim(metaStruct(iVar)%vartype))
    case('midSnow','ifcSnow'); err=20; message=trim(message)//'nSnow is missing'; return
    case('midSoil','ifcSoil'); err=20; message=trim(message)//'nSoil is missing'; return
    case('midToto','ifcToto'); err=20; message=trim(message)//'nLayers is missing'; return
   end select
  end do ! loop through variables in the metadata structure
 endif

 ! initialize allocation check
 check=.false.

 ! allocate the dimension for model variables
 select type(dataStruct)
  type is (var_i);       if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); endif; return
  type is (var_d);       if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); endif; return
  type is (var_ilength); if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); endif
  type is (var_dlength); if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); endif
  class default; err=20; message=trim(message)//'unable to identify derived data type for the variable dimension'; return
 end select
 ! check errors
 if(check) then; err=20; message=trim(message)//'structure was unexpectedly allocated already'; return; endif
 if(err/=0)then; err=20; message=trim(message)//'problem allocating'; return; endif

 ! allocate the dimension for model data
 select type(dataStruct)
  type is (var_ilength); call allocateDat_int(metaStruct,nSnow,nSoil,nLayers,dataStruct,err,cmessage) 
  type is (var_dlength); call allocateDat_dp( metaStruct,nSnow,nSoil,nLayers,dataStruct,err,cmessage) 
  class default; err=20; message=trim(message)//'unable to identify derived data type for the data dimension'; return
 end select
 ! check errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine allocLocal

 ! ************************************************************************************************
 ! private subroutine allocateDat_dp: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_dp(metadata,nSnow,nSoil,nLayers, & ! input
                           varData,err,message)            ! output
 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_dlength),intent(inout)   :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 ! initialize error control
 err=0; message='allocateDat_dp/'

 ! loop through variables in the data structure
 do iVar=1,size(metadata)

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate structures
  else
   select case(trim(metadata(iVar)%vartype))
    case('scalarv'); allocate(varData%var(iVar)%dat(1),stat=err)
    case('wLength'); allocate(varData%var(iVar)%dat(nBand),stat=err)
    case('midSnow'); allocate(varData%var(iVar)%dat(nSnow),stat=err)
    case('midSoil'); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case('midToto'); allocate(varData%var(iVar)%dat(nLayers),stat=err)
    case('ifcSnow'); allocate(varData%var(iVar)%dat(0:nSnow),stat=err)
    case('ifcSoil'); allocate(varData%var(iVar)%dat(0:nSoil),stat=err)
    case('ifcToto'); allocate(varData%var(iVar)%dat(0:nLayers),stat=err)
    case('routing'); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err)
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=special (and valid) case that is allocated later (initialize with zero-length vector)
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   ! check error
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; endif
   ! set to missing
   varData%var(iVar)%dat(:) = missingDouble
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_dp

 ! ************************************************************************************************
 ! private subroutine allocateDat_int: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_int(metadata,nSnow,nSoil,nLayers, & ! input
                            varData,err,message)            ! output
 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_ilength),intent(inout)   :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 ! initialize error control
 err=0; message='allocateDat_int/'

 ! loop through variables in the data structure
 do iVar=1,size(metadata)

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate structures
  else
   select case(trim(metadata(iVar)%vartype))
    case('scalarv'); allocate(varData%var(iVar)%dat(1),stat=err)
    case('wLength'); allocate(varData%var(iVar)%dat(nBand),stat=err)
    case('midSnow'); allocate(varData%var(iVar)%dat(nSnow),stat=err)
    case('midSoil'); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case('midToto'); allocate(varData%var(iVar)%dat(nLayers),stat=err)
    case('ifcSnow'); allocate(varData%var(iVar)%dat(0:nSnow),stat=err)
    case('ifcSoil'); allocate(varData%var(iVar)%dat(0:nSoil),stat=err)
    case('ifcToto'); allocate(varData%var(iVar)%dat(0:nLayers),stat=err)
    case('routing'); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err)
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=special (and valid) case that is allocated later (initialize with zero-length vector)
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   ! check error
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; endif
   ! set to missing
   varData%var(iVar)%dat(:) = missingInteger
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_int

end module allocspace_module
