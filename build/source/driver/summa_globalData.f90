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

module summa_globalData
! used to declare and allocate global summa data structures

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing double precision number

! size of data structures
USE var_lookup,only:maxvarForc      ! forcing data:             maximum number variables
USE var_lookup,only:maxvarProg      ! prognostic variables:     maximum number variables
USE var_lookup,only:maxvarDiag      ! diagnostic variables:     maximum number variables
USE var_lookup,only:maxvarFlux      ! model fluxes:             maximum number variables
USE var_lookup,only:maxvarIndx      ! model indices:            maximum number variables
USE var_lookup,only:maxvarBvar      ! basin-average variables:  maximum number variables

! metadata structures
USE globalData,only:time_meta,forc_meta,attr_meta,type_meta ! metadata structures
USE globalData,only:prog_meta,diag_meta,flux_meta           ! metadata structures
USE globalData,only:mpar_meta,indx_meta                     ! metadata structures
USE globalData,only:bpar_meta,bvar_meta                     ! metadata structures
USE globalData,only:averageFlux_meta                        ! metadata for time-step average fluxes

! statistics metadata structures
USE globalData,only:statForc_meta                           ! child metadata for stats
USE globalData,only:statProg_meta                           ! child metadata for stats
USE globalData,only:statDiag_meta                           ! child metadata for stats
USE globalData,only:statFlux_meta                           ! child metadata for stats
USE globalData,only:statIndx_meta                           ! child metadata for stats
USE globalData,only:statBvar_meta                           ! child metadata for stats

! mapping from original to child structures
USE globalData,only:forcChild_map                           ! index of the child data structure: stats forc
USE globalData,only:progChild_map                           ! index of the child data structure: stats prog
USE globalData,only:diagChild_map                           ! index of the child data structure: stats diag
USE globalData,only:fluxChild_map                           ! index of the child data structure: stats flux
USE globalData,only:indxChild_map                           ! index of the child data structure: stats indx
USE globalData,only:bvarChild_map                           ! index of the child data structure: stats bvar

! safety: set private unless specified otherwise
implicit none
private
public::summa_defineGlobalData
contains

 subroutine summa_defineGlobalData(err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nrtype                                                  ! variable types, etc.
 ! subroutines and functions: initial priming
 USE,intrinsic :: ieee_arithmetic                            ! IEEE arithmetic (obviously)
 ! subroutines and functions: define metadata structures
 USE popMetadat_module,only:popMetadat                       ! module to populate metadata structures
 USE flxMapping_module,only:flxMapping                       ! module to map fluxes to states
 USE checkStruc_module,only:checkStruc                       ! module to check metadata structures
 USE childStruc_module,only:childStruc                       ! module to create a child data structure
 ! miscellaneous global data
 USE globalData,only:dNaN                                    ! double precision NaN
 USE globalData,only:doJacobian                              ! flag to compute the Jacobian
 USE globalData,only:structInfo                              ! information on the data structures
 ! named variables that describe elements of child  model structures
 USE var_lookup,only:iLookVarType                            ! look-up values for variable type structure
 USE var_lookup,only:childFLUX_MEAN                          ! look-up values for timestep-average model fluxes
 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 integer(i4b),intent(out)              :: err                ! error code
 character(*),intent(out)              :: message            ! error message
 ! local variables
 character(LEN=256)                    :: cmessage           ! error message of downwind routine
 logical(lgt), dimension(maxvarFlux)   :: flux_mask          ! mask defining desired flux variables
 logical(lgt), dimension(maxvarForc)   :: statForc_mask      ! mask defining forc stats
 logical(lgt), dimension(maxvarProg)   :: statProg_mask      ! mask defining prog stats
 logical(lgt), dimension(maxvarDiag)   :: statDiag_mask      ! mask defining diag stats
 logical(lgt), dimension(maxvarFlux)   :: statFlux_mask      ! mask defining flux stats
 logical(lgt), dimension(maxvarIndx)   :: statIndx_mask      ! mask defining indx stats
 logical(lgt), dimension(maxvarBvar)   :: statBvar_mask      ! mask defining bvar stats
 integer(i4b)                          :: iStruct            ! index of data structure
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_defineGlobalData/'

 ! initialize the Jacobian flag
 doJacobian=.false.        ! initialize the Jacobian flag

 ! define double precision NaNs (shared in globalData)
 dNaN = ieee_value(1._rkind, ieee_quiet_nan)

 ! populate metadata for all model variables
 call popMetadat(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! define mapping between fluxes and states
 call flxMapping(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! check data structures
 call checkStruc(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! define the mask to identify the subset of variables in the "child" data structure (just scalar variables)
 flux_mask = (flux_meta(:)%vartype==iLookVarType%scalarv)

 ! create the averageFlux metadata structure
 call childStruc(flux_meta, flux_mask, averageFlux_meta, childFLUX_MEAN, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! child metadata structures - so that we do not carry full stats structures around everywhere
 ! only carry stats for variables with output frequency > model time step
 statForc_mask = (forc_meta(:)%vartype==iLookVarType%scalarv.and.forc_meta(:)%varDesire)
 statProg_mask = (prog_meta(:)%vartype==iLookVarType%scalarv.and.prog_meta(:)%varDesire)
 statDiag_mask = (diag_meta(:)%vartype==iLookVarType%scalarv.and.diag_meta(:)%varDesire)
 statFlux_mask = (flux_meta(:)%vartype==iLookVarType%scalarv.and.flux_meta(:)%varDesire)
 statIndx_mask = (indx_meta(:)%vartype==iLookVarType%scalarv.and.indx_meta(:)%varDesire)
 statBvar_mask = (bvar_meta(:)%vartype==iLookVarType%scalarv.and.bvar_meta(:)%varDesire)

 ! create the stats metadata structures
 do iStruct=1,size(structInfo)
  select case (trim(structInfo(iStruct)%structName))
   case('forc'); call childStruc(forc_meta,statForc_mask,statForc_meta,forcChild_map,err,cmessage)
   case('prog'); call childStruc(prog_meta,statProg_mask,statProg_meta,progChild_map,err,cmessage)
   case('diag'); call childStruc(diag_meta,statDiag_mask,statDiag_meta,diagChild_map,err,cmessage)
   case('flux'); call childStruc(flux_meta,statFlux_mask,statFlux_meta,fluxChild_map,err,cmessage)
   case('indx'); call childStruc(indx_meta,statIndx_mask,statIndx_meta,indxChild_map,err,cmessage)
   case('bvar'); call childStruc(bvar_meta,statBvar_mask,statBvar_meta,bvarChild_map,err,cmessage)
  end select
  ! check errors
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[statistics for =  '//trim(structInfo(iStruct)%structName)//']'; return; endif
 end do ! iStruct

 ! set all stats metadata to correct var types
 statForc_meta(:)%vartype = iLookVarType%outstat
 statProg_meta(:)%vartype = iLookVarType%outstat
 statDiag_meta(:)%vartype = iLookVarType%outstat
 statFlux_meta(:)%vartype = iLookVarType%outstat
 statIndx_meta(:)%vartype = iLookVarType%outstat
 statBvar_meta(:)%vartype = iLookVarType%outstat

 end subroutine summa_defineGlobalData

end module summa_globalData
