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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

! **************************************************************************************************
! SUMMA parameter sampling
!
! Provides SUMMA-specific routines for generating parameter samples and constructing complete
! parameter override vectors for model evaluation. This module separates parameter-sampling
! strategy from the parameter specification, model evaluation, and calibration workflow.
! **************************************************************************************************

module summa_parameter_sampling

  ! data types
  USE nr_type,    only: i4b,rkind

  ! parameter-search data types
  USE parameter_search, only: parameter_spec,parameter_search_info

  implicit none
  private

  public :: generate_parameter_sample 

contains

  ! **************************************************************************************************
  ! Generate a SUMMA parameter sample.
  !
  ! Generates one feasible parameter vector using the configured parameter-search strategy and
  ! constructs the complete SUMMA parameter override vector for model evaluation. The sampled
  ! parameter vector contains only parameters included in the search, whereas the override vector
  ! also includes non-sampled parameters required by calibration constraints.
  ! **************************************************************************************************
  
  subroutine generate_parameter_sample(param_spec,search,          &
                                       param_value,param_override, &
                                       err,message)
  
    ! parameter sampling
    USE parameter_search, only: sample_parameters
  
    ! SUMMA parameter overrides
    USE summa_parameter_spec, only: build_summa_parameter_overrides
  
    implicit none
  
    ! dummy variables
    type(parameter_spec),        intent(in)  :: param_spec        ! SUMMA parameter specification
    type(parameter_search_info), intent(in)  :: search            ! parameter-search information
    real(rkind),                 intent(out) :: param_value(:)    ! sampled parameter values
    real(rkind),                 intent(out) :: param_override(:) ! complete SUMMA parameter overrides
    integer(i4b),                intent(out) :: err               ! error code
    character(*),                intent(out) :: message           ! error message
  
    ! local variables
    character(len=256) :: cmessage
  
    err=0
    message='generate_parameter_sample/'
  
    ! generate a feasible parameter vector in the parameter-search space
    call sample_parameters(search,param_value,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif
  
    ! construct the complete SUMMA parameter override vector
    call build_summa_parameter_overrides(param_spec,         &
                                         search%param_names, &
                                         param_value,        &
                                         param_override,     &
                                         err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif
  
  end subroutine generate_parameter_sample

end module summa_parameter_sampling
