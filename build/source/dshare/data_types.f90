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

MODULE data_types
 ! used to define model data structures
 USE nrtype
 USE multiconst,only:integerMissing
 implicit none
 private

 ! ***********************************************************************************************************
 ! Define the model decisions
 ! ***********************************************************************************************************
 ! the model decision structure
 type,public  :: model_options
  character(len=64)                      :: cOption='notPopulatedYet'
  character(len=64)                      :: cDecision='notPopulatedYet'
  integer(i4b)                           :: iDecision=integerMissing
 end type model_options

 ! ***********************************************************************************************************
 ! Define metadata for model forcing datafile
 ! ***********************************************************************************************************
 ! define a derived type for the data in the file
 type,public  :: file_info
  character(len=256)                     :: filenmDesc='notPopulatedYet' ! name of file that describes the data
  character(len=256)                     :: filenmData='notPopulatedYet' ! name of data file
  integer(i4b)                           :: ncols                    ! number of columns in the file
  integer(i4b)                           :: ixFirstHRU               ! index of the first HRU to share the same data
  integer(i4b),allocatable               :: time_ix(:)               ! column index for each time variable
  integer(i4b),allocatable               :: data_ix(:)               ! column index for each forcing data variable
 end type file_info

 ! ***********************************************************************************************************
 ! Define metadata on model parameters
 ! ***********************************************************************************************************
 ! define a data type to store model parameter information
 type,public  :: par_info
  real(dp)                               :: default_val              ! default parameter value
  real(dp)                               :: lower_limit              ! lower bound
  real(dp)                               :: upper_limit              ! upper bound
 endtype par_info

 ! ***********************************************************************************************************
 ! Define variable metadata
 ! ***********************************************************************************************************
 ! define derived type for model variables, including name, decription, and units
 type,public :: var_info
  character(len=64)                      :: varname=''       ! variable name
  character(len=128)                     :: vardesc=''       ! variable description
  character(len=64)                      :: varunit=''       ! variable units
  character(len=32)                      :: vartype=''       ! variable type (scalar, model layers, etc.)
  logical(lgt)                           :: v_write=.FALSE.  ! flag to write variable to the output file
 endtype var_info

 ! ***********************************************************************************************************
 ! Define summary of data structures
 ! ***********************************************************************************************************
 ! data structure information
 type,public :: struct_info 
  character(len=32)                      :: structName  ! name of the data structure
  character(len=32)                      :: lookName    ! name of the look-up variables
  integer(i4b)                           :: nVar        ! number of variables in each data structure
 end type struct_info

 ! ***********************************************************************************************************
 ! Define hierarchal derived data types
 ! ***********************************************************************************************************
 ! define derived types to hold multivariate data for a single variable (different variables have different length)
 ! NOTE: use derived types here to facilitate adding the "variable" dimension
 ! ** double precision type
 type, public :: dlength
  real(dp),allocatable                       :: dat(:) 
 endtype dlength
 ! ** integer type
 type, public :: ilength
  integer(i4b),allocatable                   :: dat(:) 
 endtype ilength

 ! define derived types to hold data for multiple variables
 ! NOTE: use derived types here to facilitate adding extra dimensions (e.g., spatial)
 ! ** double precision type of variable length
 type, public :: var_dlength
  type(dlength),allocatable                  :: var(:) 
 endtype var_dlength
 ! ** integer type of variable length
 type, public :: var_ilength
  type(ilength),allocatable                  :: var(:) 
 endtype var_ilength
 ! ** double precision type of fixed length
 type, public :: var_d
  real(dp),allocatable                       :: var(:) 
 endtype var_d
 ! ** integer type of variable length
 type, public :: var_i
  integer(i4b),allocatable                   :: var(:) 
 endtype var_i

 ! define derived types to hold spatial dimensions
 ! ** double precision type of variable length
 type, public :: spatial_doubleVec
  type(var_dlength),allocatable :: hru(:)
 endtype spatial_doubleVec
 ! ** integer type of variable length
 type, public :: spatial_intVec
  type(var_ilength),allocatable :: hru(:)
 endtype spatial_intVec
 ! ** double precision type of fixed length
 type, public :: spatial_double
  type(var_d),allocatable       :: hru(:)
 endtype spatial_double
 ! ** integer type of variable length
 type, public :: spatial_int
  type(var_i),allocatable       :: hru(:)
 endtype spatial_int

END MODULE data_types

