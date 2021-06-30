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

module volicePackFida_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! named variables for snow and soil
USE globalData,only:iname_snow          ! named variables for snow
USE globalData,only:iname_soil          ! named variables for soil

! named variables for parent structures
USE var_lookup,only:iLookINDEX          ! named variables for structure elements

! physical constants
USE multiconst,only:&
                    Tfreeze,  & ! freezing point              (K)
                    LH_fus,   & ! latent heat of fusion       (J kg-1)
                    LH_vap,   & ! latent heat of vaporization (J kg-1)
                    LH_sub,   & ! latent heat of sublimation  (J kg-1)
                    iden_air, & ! intrinsic density of air    (kg m-3)
                    iden_ice, & ! intrinsic density of ice    (kg m-3)
                    iden_water  ! intrinsic density of water  (kg m-3)

! privacy
implicit none
private
public::volicePackFida
public::newsnwfall

contains


 ! ************************************************************************************************
 ! public subroutine volicePackFida: combine and sub-divide layers if necessary)
 ! ************************************************************************************************
 subroutine volicePackFida(&
                       ! input/output: model data structures
                       tooMuchMelt,                 & ! intent(in):    flag to force merge of snow layers
                       model_decisions,             & ! intent(in):    model decisions
                       mpar_data,                   & ! intent(in):    model parameters
                       indx_data,                   & ! intent(inout): type of each layer
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,                   & ! intent(inout): model fluxes for a local HRU
                       ! output
                       modifiedLayers,              & ! intent(out): flag to denote that layers were modified
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 ! external subroutine
 USE layerMerge_module,only:layerMerge   ! merge snow layers if they are too thin
 USE layerDivide_module,only:layerDivide ! sub-divide layers if they are too thick
 implicit none
 ! ------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 logical(lgt),intent(in)         :: tooMuchMelt         ! flag to denote that ice is insufficient to support melt
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_ilength),intent(inout) :: indx_data           ! type of each layer
 type(var_dlength),intent(inout) :: prog_data           ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data           ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data           ! model flux variables
 ! output
 logical(lgt),intent(out)        :: modifiedLayers      ! flag to denote that we modified the layers
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! ------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 logical(lgt)                    :: mergedLayers        ! flag to denote that layers were merged
 logical(lgt)                    :: divideLayer         ! flag to denote that a layer was divided
 ! initialize error control
 err=0; message='volicePackFida/'

 ! divide snow layers if too thick
 call layerDivide(&
                  ! input/output: model data structures
                  model_decisions,             & ! intent(in):    model decisions
                  mpar_data,                   & ! intent(in):    model parameters
                  indx_data,                   & ! intent(inout): type of each layer
                  prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                  diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,                   & ! intent(inout): model fluxes for a local HRU
                  ! output
                  divideLayer,                 & ! intent(out): flag to denote that layers were modified
                  err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=65; message=trim(message)//trim(cmessage); return; end if
 
 if(divideLayer) print *, 'divideLayer'

 ! merge snow layers if they are too thin
 call layerMerge(&
                 ! input/output: model data structures
                 tooMuchMelt,                 & ! intent(in):    flag to force merge of snow layers
                 model_decisions,             & ! intent(in):    model decisions
                 mpar_data,                   & ! intent(in):    model parameters
                 indx_data,                   & ! intent(inout): type of each layer
                 prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                 diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                 flux_data,                   & ! intent(inout): model fluxes for a local HRU
                 ! output
                 mergedLayers,                & ! intent(out): flag to denote that layers were modified
                 err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=65; message=trim(message)//trim(cmessage); return; end if
 
 if(mergedLayers) print *, 'mergedLayers'

 ! update the number of layers
 indx_data%var(iLookINDEX%nSnow)%dat(1)   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
 indx_data%var(iLookINDEX%nSoil)%dat(1)   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)
 indx_data%var(iLookINDEX%nLayers)%dat(1) = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nSoil)%dat(1)

 ! flag if layers were modified
 modifiedLayers = (mergedLayers .or. divideLayer)

 end subroutine volicePackFida


end module volicePackFida_module
