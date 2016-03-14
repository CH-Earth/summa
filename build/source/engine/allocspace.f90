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
! number of variables in each data structure
USE var_lookup,only:maxvarTime,maxvarForc,maxvarAttr,maxvarType    ! maximum number variables in each data structure
USE var_lookup,only:maxvarProg,maxvarDiag,maxvarFlux,maxvarDeriv   ! maximum number variables in each data structure
USE var_lookup,only:maxvarMpar,maxvarMvar,maxvarIndx               ! maximum number variables in each data structure
USE var_lookup,only:maxvarBpar,maxvarBvar                          ! maximum number variables in each data structure
USE var_lookup,only:maxvarDecisions                                ! maximum number of decisions
! metadata structures
USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta        ! metadata structures
USE data_struc,only:prog_meta,diag_meta,flux_meta,deriv_meta       ! metadata structures
USE data_struc,only:mpar_meta,mvar_meta,indx_meta                  ! metadata structures
USE data_struc,only:bpar_meta,bvar_meta                            ! metadata structures
USE data_struc,only:model_decisions                                ! model decision structure
! time structures
USE data_struc,only:refTime,startTime,finshTime                    ! reference time, start time, and end time for the model simulation
implicit none
private
public::init_metad
public::initStruct
! define missing values
integer(i4b),parameter :: missingInteger=-9999
real(dp),parameter     :: missingDouble=-9999._dp
! define fixed dimensions
integer(i4b),parameter :: nBand=2         ! number of spectral bands
integer(i4b),parameter :: nTimeDelay=2000 ! number of elements in the time delay histogram
! -----------------------------------------------------------------------------------------------------------------------------------
! data structure information
integer(i4b),parameter               :: nStruct=14  ! number of data structures
type info ! data structure information
 character(len=32)                   :: structName  ! name of the data structure
 character(len=32)                   :: lookName    ! name of the look-up variables
 integer(i4b)                        :: nVar        ! number of variables in each data structure
end type info
! populate structure information
type(info),parameter,dimension(nStruct) :: structInfo=(/&
                                            info('time',  'TIME' , maxvarTime ), & ! the time data structure
                                            info('forc',  'FORCE', maxvarForc ), & ! the forcing data structure
                                            info('attr',  'ATTR' , maxvarAttr ), & ! the attribute data structure
                                            info('type',  'TYPE' , maxvarType ), & ! the type data structure
                                            info('dpar',  'PARAM', maxvarMpar ), & ! the model parameter data structure
                                            info('mpar',  'PARAM', maxvarMpar ), & ! the model parameter data structure
                                            info('mvar',  'MVAR' , maxvarMvar ), & ! the model variable data structure
                                            info('bpar',  'BPAR' , maxvarBpar ), & ! the basin parameter data structure
                                            info('bvar',  'BVAR' , maxvarBvar ), & ! the basin variable data structure
                                            info('indx',  'INDEX', maxvarIndx ), & ! the model index data structure
                                            info('prog',  'PROG',  maxvarProg),  & ! the prognostic (state) variable data structure
                                            info('diag',  'DIAG' , maxvarDiag ), & ! the diagnostic variable data structure
                                            info('flux',  'FLUX' , maxvarFlux ), & ! the flux data structure
                                            info('deriv', 'DERIV', maxvarDeriv) /) ! the model derivative data structure
! -----------------------------------------------------------------------------------------------------------------------------------
contains


 ! ************************************************************************************************
 ! public subroutine init_metad: initialize metadata structures
 ! ************************************************************************************************
 subroutine init_metad(err,message)
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! initialize errors
 err=0; message="init_model/"

 ! ensure metadata structures are deallocated
 if (allocated(time_meta))  deallocate(time_meta)   ! 1
 if (allocated(forc_meta))  deallocate(forc_meta)   ! 2
 if (allocated(attr_meta))  deallocate(attr_meta)   ! 3
 if (allocated(type_meta))  deallocate(type_meta)   ! 4
 if (allocated(mpar_meta))  deallocate(mpar_meta)   ! 5
 if (allocated(mvar_meta))  deallocate(mvar_meta)   ! 6
 if (allocated(bpar_meta))  deallocate(bpar_meta)   ! 7
 if (allocated(bvar_meta))  deallocate(bvar_meta)   ! 8
 if (allocated(indx_meta))  deallocate(indx_meta)   ! 9
 if (allocated(prog_meta))  deallocate(prog_meta)  ! 10
 if (allocated(diag_meta))  deallocate(diag_meta)   ! 11
 if (allocated(flux_meta))  deallocate(flux_meta)   ! 12
 if (allocated(deriv_meta)) deallocate(deriv_meta)  ! 13

 ! allocate metadata structures
 allocate(time_meta(maxvarTime),forc_meta(maxvarForc),attr_meta(maxvarAttr),type_meta(maxvarType),&
          prog_meta(maxvarProg),diag_meta(maxvarDiag),flux_meta(maxvarFlux),deriv_meta(maxvarDeriv),&
          mpar_meta(maxvarMpar),mvar_meta(maxvarMvar),indx_meta(maxvarIndx),&
          bpar_meta(maxvarBpar),bvar_meta(maxvarBvar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateMetadata"; return; endif

 ! ensure decisions structure is deallocated
 if(allocated(model_decisions)) deallocate(model_decisions)

 ! allocate space for the model decisions
 allocate(model_decisions(maxvarDecisions),stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemAllocateModelDecisions"; return; endif

 ! ensure model time structures are deallolcated
 if (allocated(refTime))   deallocate(refTime)
 if (allocated(startTime)) deallocate(startTime)
 if (allocated(finshTime)) deallocate(finshTime)

 ! allocate space for the model time structures (1st level)
 allocate(refTime,startTime,finshTime,stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemAllocateTimeStructures(1st level)"; return; endif

 ! allocate space for the model time structures (2nd level)
 allocate(refTime%var(maxvarTime),startTime%var(maxvarTime),finshTime%var(maxvarTime),stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemAllocateTimeStructures(2nd level)"; return; endif

 ! set time vectors to missing
 refTime%var(:)   = missingInteger
 startTime%var(:) = missingInteger
 finshTime%var(:) = missingInteger

 end subroutine init_metad


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
                       dparStruct, &    ! default model parameters
                       mparStruct, &    ! model parameters
                       mvarStruct, &    ! model variables
                       indxStruct, &    ! model indices
                       progStruct, &    ! model prognostic (state) variables
                       diagStruct, &    ! model diagnostic variables
                       fluxStruct, &    ! model fluxes
                       derivStruct,&    ! model derivatives
                       bparStruct, &    ! basin-average parameters
                       bvarStruct, &    ! basin-average variables
                       ! output: error control
                       err,message)
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_i,               & ! x%var(:)            (i4b)
                     var_d,               & ! x%var(:)            (dp)
                     var_ilength,         & ! x%var(:)%dat        (i4b)
                     var_dlength,         & ! x%var(:)%dat        (dp)
                     spatial_int,         & ! x%hru(:)%var(:)     (i4b)
                     spatial_double,      & ! x%hru(:)%var(:)     (dp)
                     spatial_intVec,      & ! x%hru(:)%var(:)%dat (i4b)
                     spatial_doubleVec      ! x%hru(:)%var(:)%dat (dp)
 implicit none
 ! input: model control
 integer(i4b),            intent(in)   :: nHRU          ! number of HRUs
 integer(i4b),            intent(in)   :: nSnow(:)      ! number of snow layers for each HRU
 integer(i4b),            intent(in)   :: nSoil(:)      ! number of soil layers for each HRU
 ! input: data structures
 type(var_i),             intent(out)  :: timeStruct    ! model time data
 type(spatial_double),    intent(out)  :: forcStruct    ! model forcing data
 type(spatial_double),    intent(out)  :: attrStruct    ! local attributes for each HRU
 type(spatial_int),       intent(out)  :: typeStruct    ! local classification of soil veg etc. for each HRU
 type(spatial_double),    intent(out)  :: dparStruct    ! default model parameters
 type(spatial_double),    intent(out)  :: mparStruct    ! model parameters
 type(spatial_doubleVec), intent(out)  :: mvarStruct    ! model variables
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
 integer(i4b)                          :: iHRU          ! index of HRUs 
 integer(i4b)                          :: iStruct       ! index of data structure 
 logical(lgt)                          :: check         ! flag to check if the structure is already allocated
 integer(i4b)                          :: nLayers       ! total number of layers in the snow+soil domian
 character(len=256)                    :: cMessage      ! error message of downwind routine
 ! initialize errors
 err=0; message="initStruct/"

 ! -----
 ! * allocate spatial dimensions...
 ! --------------------------------

 ! loop through data structures
 do iStruct=1,nStruct
  check=.false. ! initialize check
  ! allocate spatial dimension
  select case(trim(structInfo(iStruct)%structName))
   case('time'); cycle   ! do nothing: no spatial dimension 
   case('forc'); if(allocated(forcStruct%hru))then; check=.true.; else; allocate(forcStruct%hru(nHRU),stat=err); endif
   case('attr'); if(allocated(attrStruct%hru))then; check=.true.; else; allocate(attrStruct%hru(nHRU),stat=err); endif 
   case('type'); if(allocated(typeStruct%hru))then; check=.true.; else; allocate(typeStruct%hru(nHRU),stat=err); endif 
   case('dpar'); if(allocated(dparStruct%hru))then; check=.true.; else; allocate(dparStruct%hru(nHRU),stat=err); endif
   case('mpar'); if(allocated(mparStruct%hru))then; check=.true.; else; allocate(mparStruct%hru(nHRU),stat=err); endif
   case('mvar'); if(allocated(mvarStruct%hru))then; check=.true.; else; allocate(mvarStruct%hru(nHRU),stat=err); endif
   case('indx'); if(allocated(indxStruct%hru))then; check=.true.; else; allocate(indxStruct%hru(nHRU),stat=err); endif
   case('prog'); if(allocated(progStruct%hru))then; check=.true.; else; allocate(progStruct%hru(nHRU),stat=err); endif
   case('diag'); if(allocated(diagStruct%hru))then; check=.true.; else; allocate(diagStruct%hru(nHRU),stat=err); endif 
   case('flux'); if(allocated(fluxStruct%hru))then; check=.true.; else; allocate(fluxStruct%hru(nHRU),stat=err); endif
   case('deriv');cycle   ! do nothing: no spatial dimension
   case('bpar'); cycle   ! do nothing: no spatial dimension 
   case('bvar'); cycle   ! do nothing: no spatial dimension
   case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
  end select
  ! check errors
  if(check) then; err=20; message=trim(message)//'structure '//trim(structInfo(iStruct)%structName)//' was unexpectedly allocated already'; return; endif
  if(err/=0)then; err=20; message=trim(message)//'problem allocating the '//trim(structInfo(iStruct)%structName)//' structure'; return; endif
 end do  ! looping through data structures

 ! -----
 ! * allocate variable dimensions...
 ! ---------------------------------

 do iStruct=1,nStruct  ! loop through data structures
  do iHRU=1,nHRU  ! loop through HRUs
   check=.false. ! initialize check
   ! cycle for variables without a HRU dimension
   select case(trim(structInfo(iStruct)%structName))
    case('time','deriv','bpar','bvar'); if(iHRU>1) cycle
   end select
   ! allocate variable dimension
   select case(trim(structInfo(iStruct)%structName))
    case('time'); if(allocated(timeStruct%var)          )then; check=.true.; else; allocate(timeStruct%var(structInfo(iStruct)%nVar),stat=err); endif
    case('forc'); if(allocated(forcStruct%hru(iHRU)%var))then; check=.true.; else; allocate(forcStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('attr'); if(allocated(attrStruct%hru(iHRU)%var))then; check=.true.; else; allocate(attrStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('type'); if(allocated(typeStruct%hru(iHRU)%var))then; check=.true.; else; allocate(typeStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('dpar'); if(allocated(dparStruct%hru(iHRU)%var))then; check=.true.; else; allocate(dparStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('mpar'); if(allocated(mparStruct%hru(iHRU)%var))then; check=.true.; else; allocate(mparStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('mvar'); if(allocated(mvarStruct%hru(iHRU)%var))then; check=.true.; else; allocate(mvarStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('indx'); if(allocated(indxStruct%hru(iHRU)%var))then; check=.true.; else; allocate(indxStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('prog'); if(allocated(progStruct%hru(iHRU)%var))then; check=.true.; else; allocate(progStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif 
    case('diag'); if(allocated(diagStruct%hru(iHRU)%var))then; check=.true.; else; allocate(diagStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('flux'); if(allocated(fluxStruct%hru(iHRU)%var))then; check=.true.; else; allocate(fluxStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('deriv');if(allocated(derivStruct%var)         )then; check=.true.; else; allocate(derivStruct%var(structInfo(iStruct)%nVar),stat=err); endif
    case('bpar'); if(allocated(bparStruct%var)          )then; check=.true.; else; allocate(bparStruct%var(structInfo(iStruct)%nVar),stat=err); endif 
    case('bvar'); if(allocated(bvarStruct%var)          )then; check=.true.; else; allocate(bvarStruct%var(structInfo(iStruct)%nVar),stat=err); endif
    case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
   end select
   ! check errors
   if(check) then; err=20;write(message,'(a,i0)') trim(message)//'structure '//trim(structInfo(iStruct)%structName)//' was unexpectedly allocated already for HRU ', iHRU; return; endif
   if(err/=0)then; err=20;write(message,'(a,i0)') trim(message)//'problem allocating the '//trim(structInfo(iStruct)%structName)//' structure for HRU ', iHRU; return; endif
  end do  ! looping through data structures
 end do  ! loop through HRUs

 ! -----
 ! * allocate specific variables...
 ! --------------------------------

 do iStruct=1,nStruct  ! loop through data structures
  do iHRU=1,nHRU  ! loop through HRUs
   ! get the number of layers in the snow+soil sub-domain
   nLayers = nSnow(iHRU) + nSoil(iHRU)
   if(trim(structInfo(iStruct)%structName)=='bvar' .and. iHRU>1) cycle
   ! cycle for variables without a HRU dimension
   select case(trim(structInfo(iStruct)%structName))
    case('deriv','bvar'); if(iHRU>1) cycle
   end select
   ! allocate data dimension
   select case(trim(structInfo(iStruct)%structName))
    case('time'); cycle  ! do nothing: no data dimension 
    case('forc'); cycle  ! do nothing: no data dimension
    case('attr'); cycle  ! do nothing: no data dimension
    case('type'); cycle  ! do nothing: no data dimension
    case('dpar'); cycle  ! do nothing: no data dimension
    case('mpar'); cycle  ! do nothing: no data dimension
    case('mvar'); call allocateDat_dp( mvar_meta, nSnow(iHRU),nSoil(iHRU),nLayers,mvarStruct%hru(iHRU),err,cmessage)
    case('indx'); call allocateDat_int(indx_meta, nSnow(iHRU),nSoil(iHRU),nLayers,indxStruct%hru(iHRU),err,cmessage)
    case('prog'); call allocateDat_dp( prog_meta, nSnow(iHRU),nSoil(iHRU),nLayers,progStruct%hru(iHRU),err,cmessage) 
    case('diag'); call allocateDat_dp( diag_meta, nSnow(iHRU),nSoil(iHRU),nLayers,diagStruct%hru(iHRU),err,cmessage)  ! NOTE: no HRU dimension
    case('flux'); call allocateDat_dp( flux_meta, nSnow(iHRU),nSoil(iHRU),nLayers,fluxStruct%hru(iHRU),err,cmessage)  ! NOTE: no HRU dimension
    case('deriv');call allocateDat_dp( deriv_meta,nSnow(iHRU),nSoil(iHRU),nLayers,derivStruct,err,cmessage) ! NOTE: no HRU dimension
    case('bpar'); cycle  ! do nothing: no data dimension 
    case('bvar'); call allocateDat_dp( bvar_meta,nSnow(iHRU),nSoil(iHRU),nLayers,bvarStruct,err,cmessage)   ! NOTE: no HRU dimension
    case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
   end select
   ! check errors
   if(err/=0)then
    write(message,'(a,i0)') trim(message)//trim(cmessage)//'; attempting to allocate the '//trim(structInfo(iStruct)%structName)//' structure for HRU ', iHRU
    err=20; return
   endif
  end do  ! looping through data structures
 end do  ! loop through HRUs

 end subroutine initStruct

 ! ************************************************************************************************
 ! private subroutine allocateDat_dp: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_dp(metadata,nSnow,nSoil,nLayers, & ! input
                           varData,err,message)            ! output
 ! provide access to derived data types
 USE data_struc,only:var_info       ! metadata structure
 USE data_struc,only:var_dlength    ! x%var(:)%dat (dp)
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

  ! allocate dimension
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
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=initialize with zero-length vector
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; endif
   varData%var(iVar)%dat(:) = missingDouble
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_dp

 ! ************************************************************************************************
 ! private subroutine allocateDat_int: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_int(metadata,nSnow,nSoil,nLayers, & ! input
                            varData,err,message)            ! output
 ! provide access to derived data types
 USE data_struc,only:var_info       ! metadata structure
 USE data_struc,only:var_ilength    ! x%var(:)%dat (i4b)
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

  ! allocate dimension
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
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=initialize with zero-length vector
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; endif
   varData%var(iVar)%dat(:) = missingInteger
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_int

end module allocspace_module
