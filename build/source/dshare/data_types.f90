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

MODULE data_types
 ! used to define model data structures
 USE nr_type, integerMissing=>nr_integerMissing
 USE var_lookup,only:maxvarFreq
 USE var_lookup,only:maxvarStat
 USE var_lookup,only:maxvarDecisions  ! maximum number of decisions
 USE var_lookup,only:iLookPARAM       ! lookup indices for model parameters
 USE var_lookup,only:iLookFLUX        ! lookup indices for flux data
 USE var_lookup,only:iLookDERIV       ! lookup indices for derivative data
 USE var_lookup,only:iLookFORCE       ! lookup indices for forcing data 
 USE var_lookup,only:iLookDIAG        ! lookup indices for diagnostic variable data
 USE var_lookup,only:iLookDECISIONS   ! lookup indices for elements of the decision structure
 USE var_lookup,only:iLookPROG        ! lookup indices for prognostic variables
 implicit none
 private

 ! ***********************************************************************************************************
 ! Define the model decisions
 ! ***********************************************************************************************************
 ! the model decision structure
 type,public  :: model_options
  character(len=64)                      :: cOption   = 'notPopulatedYet'
  character(len=64)                      :: cDecision = 'notPopulatedYet'
  integer(i4b)                           :: iDecision = integerMissing
 end type model_options

 ! ***********************************************************************************************************
 ! Define metadata for model forcing datafile
 ! ***********************************************************************************************************
 ! define a derived type for the data in the file
 type,public  :: file_info
  character(len=256)                     :: filenmData='notPopulatedYet'  ! name of data file
  integer(i4b)                           :: nVars                         ! number of variables in the file
  integer(i4b)                           :: nTimeSteps                    ! number of variables in the file
  integer(i4b),allocatable               :: var_ix(:)                     ! index of each forcing data variable in the data structure
  integer(i4b),allocatable               :: data_id(:)                    ! netcdf variable id for each forcing data variable
  character(len=256),allocatable         :: varName(:)                    ! netcdf variable name for each forcing data variable
  real(rkind)                            :: firstJulDay                   ! first julian day in forcing file
  real(rkind)                            :: convTime2Days                 ! factor to convert time to days
 end type file_info

 ! ***********************************************************************************************************
 ! Define metadata on model parameters
 ! ***********************************************************************************************************
 ! define a data type to store model parameter information
 type,public  :: par_info
  real(rkind)                            :: default_val                   ! default parameter value
  real(rkind)                            :: lower_limit                   ! lower bound
  real(rkind)                            :: upper_limit                   ! upper bound
 endtype par_info

 ! ***********************************************************************************************************
 ! Define variable metadata
 ! ***********************************************************************************************************
 ! define derived type for model variables, including name, description, and units
 type,public :: var_info
  character(len=64)                      :: varName   = 'empty'           ! variable name
  character(len=128)                     :: vardesc   = 'empty'           ! variable description
  character(len=64)                      :: varunit   = 'empty'           ! variable units
  integer(i4b)                           :: varType   = integerMissing    ! variable type
  integer(i4b),dimension(maxvarFreq)     :: ncVarID   = integerMissing    ! netcdf variable id (missing if frequency is not desired)
  integer(i4b),dimension(maxvarFreq)     :: statIndex = integerMissing    ! index of desired statistic for temporal aggregation
  logical(lgt)                           :: varDesire = .false.           ! flag to denote if the variable is desired for model output
 endtype var_info

 ! define extended data type (include indices to map onto parent data type)
 type,extends(var_info),public :: extended_info
  integer(i4b)                           :: ixParent                      ! index in the parent data structure
 endtype extended_info

 ! define extended data type (includes named variables for the states affected by each flux)
 type,extends(var_info),public :: flux2state
  integer(i4b)                           :: state1                        ! named variable of the 1st state affected by the flux
  integer(i4b)                           :: state2                        ! named variable of the 2nd state affected by the flux
 endtype flux2state

 ! ***********************************************************************************************************
 ! Define summary of data structures
 ! ***********************************************************************************************************
 ! data structure information
 type,public :: struct_info
  character(len=32)                      :: structName                    ! name of the data structure
  character(len=32)                      :: lookName                      ! name of the look-up variables
  integer(i4b)                           :: nVar                          ! number of variables in each data structure
 end type struct_info

 ! ***********************************************************************************************************
 ! Define data types to map between GRUs and HRUs
 ! ***********************************************************************************************************

 ! hru info data structure
 type, public :: hru_info
  integer(i4b)                           :: hru_nc                        ! index of the hru in the netcdf file
  integer(i4b)                           :: hru_ix                        ! index of the hru in the run domain
  integer(i8b)                           :: hru_id                        ! id (non-sequential number) of the hru
  integer(i4b)                           :: nSnow                         ! number of snow layers
  integer(i4b)                           :: nSoil                         ! number of soil layers
 endtype hru_info

 ! define mapping from GRUs to the HRUs
 type, public :: gru2hru_map
  integer(i8b)                           :: gru_id                        ! id of the gru
  integer(i4b)                           :: hruCount                      ! total number of hrus in the gru
  type(hru_info), allocatable            :: hruInfo(:)                    ! basic information of HRUs within the gru
  integer(i4b)                           :: gru_nc                        ! index of gru in the netcdf file
 endtype gru2hru_map

 ! define the mapping from the HRUs to the GRUs
 type, public :: hru2gru_map
  integer(i4b)                           :: gru_ix                        ! index of gru which the hru belongs to
  integer(i4b)                           :: localHRU_ix                   ! index of a hru within a gru (start from 1 per gru)
 endtype hru2gru_map

 ! ***********************************************************************************************************
 ! Define model coupling structure
 ! ***********************************************************************************************************

 type, public :: q_coupling
   integer(i8b)                          :: id                            ! identifier of the runoff element
   real(rkind)                           :: qsim                          ! simulated runoff for this element (m s-1)
 end type q_coupling

 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------

 ! ***********************************************************************************************************
 ! Define hierarchal derived data types
 ! ***********************************************************************************************************
  ! define derived types to hold look-up tables for each soil layer
 ! ** double precision type
 type, public :: dLookup
  real(rkind),allocatable                :: lookup(:)                     ! lookup(:)
 endtype dLookup
 ! ** double precision type for a variable number of soil layers; variable length
 type, public :: vLookup
  type(dLookup),allocatable              :: var(:)                        ! var(:)%lookup(:)
 endtype vLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: zLookup
  type(vLookup),allocatable              :: z(:)                          ! z(:)%var(:)%lookup(:)
 endtype zLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: hru_z_vLookup
  type(zLookup),allocatable              :: hru(:)                        ! hru(:)%z(:)%var(:)%lookup(:)
 endtype hru_z_vLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: gru_hru_z_vLookup
  type(hru_z_vLookup),allocatable        :: gru(:)                        ! gru(:)%hru(:)%z(:)%var(:)%lookup(:)
 endtype gru_hru_z_vLookup
 ! define derived types to hold multiVariate data for a single variable (different variables have different length)
 ! NOTE: use derived types here to facilitate adding the "variable" dimension
 ! ** double precision type
 type, public :: dlength
  real(rkind),allocatable                :: dat(:)                        ! dat(:)
 endtype dlength
 ! ** integer type (4 byte)
 type, public :: ilength
  integer(i4b),allocatable               :: dat(:)                        ! dat(:)
 endtype ilength
 ! ** integer type (8 byte)
 type, public :: i8length
  integer(i8b),allocatable               :: dat(:)                        ! dat(:)
 endtype i8length
 ! ** logical type
 type, public :: flagVec
  logical(lgt),allocatable               :: dat(:)                        ! dat(:)
 endtype flagVec

 ! define derived types to hold data for multiple variables
 ! NOTE: use derived types here to facilitate adding extra dimensions (e.g., spatial)

 ! ** double precision type of variable length
 type, public :: var_dlength
  type(dlength),allocatable              :: var(:)                        ! var(:)%dat
 endtype var_dlength
 ! ** integer type of variable length (4 byte)
 type, public :: var_ilength
  type(ilength),allocatable              :: var(:)                        ! var(:)%dat
 endtype var_ilength
 ! ** integer type of variable length (8 byte)
 type, public :: var_i8length
  type(i8length),allocatable             :: var(:)                        ! var(:)%dat
 endtype var_i8length
 ! ** logical type of variable length
 type, public :: var_flagVec
  type(flagVec),allocatable              :: var(:)                        ! var(:)%dat
 endtype var_flagVec

 ! ** double precision type of fixed length
 type, public :: var_d
  real(rkind),allocatable                :: var(:)                        ! var(:)
 endtype var_d
 ! ** integer type of fixed length (4 byte)
 type, public :: var_i
  integer(i4b),allocatable               :: var(:)                        ! var(:)
 endtype var_i
 ! ** integer type of fixed length (8 byte)
 type, public :: var_i8
  integer(i8b),allocatable               :: var(:)                        ! var(:)
 endtype var_i8

 ! ** double precision type of fixed length
 type, public :: hru_d
  real(rkind),allocatable                :: hru(:)                        ! hru(:)
 endtype hru_d
 ! ** integer type of fixed length (4 byte)
 type, public :: hru_i
  integer(i4b),allocatable               :: hru(:)                        ! hru(:)
 endtype hru_i
 ! ** integer type of fixed length (8 byte)
 type, public :: hru_i8
  integer(i8b),allocatable               :: hru(:)                        ! hru(:)
 endtype hru_i8

 ! define derived types to hold JUST the HRU dimension
 ! ** double precision type of variable length
 type, public :: hru_doubleVec
  type(var_dlength),allocatable          :: hru(:)                        ! hru(:)%var(:)%dat
 endtype hru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: hru_intVec
  type(var_ilength),allocatable          :: hru(:)                        ! hru(:)%var(:)%dat
 endtype hru_intVec
 ! ** double precision type of fixed length
 type, public :: hru_double
  type(var_d),allocatable                :: hru(:)                        ! hru(:)%var(:)
 endtype hru_double
 ! ** integer type of fixed length (4 byte)
 type, public :: hru_int
  type(var_i),allocatable                :: hru(:)                        ! hru(:)%var(:)
 endtype hru_int
 ! ** integer type of fixed length (8 byte)
 type, public :: hru_int8
  type(var_i8),allocatable               :: hru(:)                        ! hru(:)%var(:)
 endtype hru_int8

 ! define derived types to hold JUST the HRU dimension
 ! ** double precision type of variable length
 type, public :: gru_doubleVec
  type(var_dlength),allocatable          :: gru(:)                        ! gru(:)%var(:)%dat
 endtype gru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: gru_intVec
  type(var_ilength),allocatable          :: gru(:)                        ! gru(:)%var(:)%dat
 endtype gru_intVec
 ! ** double precision type of fixed length
 type, public :: gru_double
  type(var_d),allocatable                :: gru(:)                        ! gru(:)%var(:)
 endtype gru_double
 ! ** integer type of variable length (4 byte)
 type, public :: gru_int
  type(var_i),allocatable                :: gru(:)                        ! gru(:)%var(:)
 endtype gru_int
 ! ** integer type of variable length (8 byte)
 type, public :: gru_int8
  type(var_i8),allocatable               :: gru(:)                        ! gru(:)%var(:)
 endtype gru_int8

 ! define derived types to hold BOTH the GRU and HRU dimension
 ! ** double precision type of variable length
 type, public :: gru_hru_doubleVec
  type(hru_doubleVec),allocatable        :: gru(:)                        ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: gru_hru_intVec
  type(hru_intVec),allocatable           :: gru(:)                        ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_intVec
 ! ** double precision type of fixed length
 type, public :: gru_hru_double
  type(hru_double),allocatable           :: gru(:)                        ! gru(:)%hru(:)%var(:)
 endtype gru_hru_double
 ! ** integer type of variable length (4 byte)
 type, public :: gru_hru_int
  type(hru_int),allocatable              :: gru(:)                        ! gru(:)%hru(:)%var(:)
 endtype gru_hru_int
 ! ** integer type of variable length (8 byte)
 type, public :: gru_hru_int8
  type(hru_int8),allocatable             :: gru(:)                        ! gru(:)%hru(:)%var(:)
 endtype gru_hru_int8
 ! ** double precision type of fixed length
 type, public :: gru_d
  type(hru_d),allocatable                :: gru(:)                        ! gru(:)%hru(:)
 endtype gru_d
 ! ** integer type of fixed length
 type, public :: gru_i
  type(hru_i),allocatable                :: gru(:)                        ! gru(:)%hru(:)
 endtype gru_i

 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------------------------------------------
 !
 integer(i4b),parameter :: len_msg=256 ! length of character string used in class definitions

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subroutines in computFlux
 ! ***********************************************************************************************************
 ! Note: class procedures are located in the contains block of this (data_types) module
 ! ** vegNrgFlux
 type, public :: in_type_vegNrgFlux ! class for intent(in) arguments in vegNrgFlux call
   logical(lgt)             :: firstSubStep                      ! intent(in): flag to indicate if we are processing the first sub-step
   logical(lgt)             :: firstFluxCall                     ! intent(in): flag to indicate if we are processing the first flux call
   logical(lgt)             :: computeVegFlux                    ! intent(in): flag to indicate if we need to compute fluxes over vegetation
   logical(lgt)             :: checkLWBalance                    ! intent(in): flag to check longwave balance
   real(rkind)              :: upperBoundTemp                    ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
   real(rkind)              :: scalarCanairTempTrial             ! intent(in): trial value of the canopy air space temperature (K)
   real(rkind)              :: scalarCanopyTempTrial             ! intent(in): trial value of canopy temperature (K)
   real(rkind)              :: mLayerTempTrial_1                 ! intent(in): trial value of ground temperature (K)
   real(rkind)              :: scalarCanopyIceTrial              ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
   real(rkind)              :: scalarCanopyLiqTrial              ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
   real(rkind)              :: dCanLiq_dTcanopy                  ! intent(in): derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)
  contains
   procedure :: initialize => initialize_in_vegNrgFlux
 end type in_type_vegNrgFlux

 type, public :: out_type_vegNrgFlux ! class for intent(out) arguments in vegNrgFlux call
   real(rkind)              :: scalarCanopyTranspiration               ! intent(out): canopy transpiration (kg m-2 s-1)
   real(rkind)              :: scalarCanopyEvaporation                 ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
   real(rkind)              :: scalarGroundEvaporation                 ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
   real(rkind)              :: scalarCanairNetNrgFlux                  ! intent(out): net energy flux for the canopy air space (W m-2)
   real(rkind)              :: scalarCanopyNetNrgFlux                  ! intent(out): net energy flux for the vegetation canopy (W m-2)
   real(rkind)              :: scalarGroundNetNrgFlux                  ! intent(out): net energy flux for the ground surface (W m-2)
   real(rkind)              :: dCanairNetFlux_dCanairTemp              ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
   real(rkind)              :: dCanairNetFlux_dCanopyTemp              ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
   real(rkind)              :: dCanairNetFlux_dGroundTemp              ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
   real(rkind)              :: dCanopyNetFlux_dCanairTemp              ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
   real(rkind)              :: dCanopyNetFlux_dCanopyTemp              ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
   real(rkind)              :: dCanopyNetFlux_dGroundTemp              ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
   real(rkind)              :: dGroundNetFlux_dCanairTemp              ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
   real(rkind)              :: dGroundNetFlux_dCanopyTemp              ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
   real(rkind)              :: dGroundNetFlux_dGroundTemp              ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
   real(rkind)              :: dCanopyEvaporation_dCanWat              ! intent(out): derivative in canopy evaporation w.r.t. canopy total water content (s-1)
   real(rkind)              :: dCanopyEvaporation_dTCanair             ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyEvaporation_dTCanopy             ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyEvaporation_dTGround             ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dGroundEvaporation_dCanWat              ! intent(out): derivative in ground evaporation w.r.t. canopy total water content (s-1)
   real(rkind)              :: dGroundEvaporation_dTCanair             ! intent(out): derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dGroundEvaporation_dTCanopy             ! intent(out): derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dGroundEvaporation_dTGround             ! intent(out): derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyTrans_dCanWat                    ! intent(out): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
   real(rkind)              :: dCanopyTrans_dTCanair                   ! intent(out): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyTrans_dTCanopy                   ! intent(out): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyTrans_dTGround                   ! intent(out): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyNetFlux_dCanWat                  ! intent(out): derivative in net canopy fluxes w.r.t. canopy total water content (J kg-1 s-1)
   real(rkind)              :: dGroundNetFlux_dCanWat                  ! intent(out): derivative in net ground fluxes w.r.t. canopy total water content (J kg-1 s-1)
   integer(i4b)             :: err                                     ! intent(out): error code
   character(len=len_msg)   :: cmessage                                ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_vegNrgFlux
 end type out_type_vegNrgFlux
 ! ** end vegNrgFlux

 ! ** snowSoilNrgFlux
 type, public :: in_type_snowSoilNrgFlux ! class for intent(in) arguments in snowSoilNrgFlux call
   logical(lgt)             :: scalarSolution                    ! intent(in): flag to denote if implementing the scalar solution
   real(rkind)              :: scalarGroundNetNrgFlux            ! intent(in): net energy flux for the ground surface (W m-2)
   real(rkind), allocatable :: iLayerLiqFluxSnow(:)              ! intent(in): liquid flux at the interface of each snow layer (m s-1)
   real(rkind), allocatable :: iLayerLiqFluxSoil(:)              ! intent(in): liquid flux at the interface of each soil layer (m s-1)
   real(rkind), allocatable :: mLayerTempTrial(:)                ! intent(in): temperature in each layer at the current iteration (m)
   real(rkind), allocatable :: dThermalC_dWatAbove(:)            ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
   real(rkind), allocatable :: dThermalC_dWatBelow(:)            ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
   real(rkind), allocatable :: dThermalC_dTempAbove(:)           ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
   real(rkind), allocatable :: dThermalC_dTempBelow(:)           ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
  contains
   procedure :: initialize => initialize_in_snowSoilNrgFlux
 end type in_type_snowSoilNrgFlux

 type, public :: io_type_snowSoilNrgFlux ! class for intent(inout) arguments in snowSoilNrgFlux call
   real(rkind)              :: dGroundNetFlux_dGroundTemp        ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  contains
   procedure :: initialize => initialize_io_snowSoilNrgFlux
   procedure :: finalize   => finalize_io_snowSoilNrgFlux
 end type io_type_snowSoilNrgFlux

 type, public :: out_type_snowSoilNrgFlux ! class for intent(inout) arguments in snowSoilNrgFlux call
   real(rkind), allocatable :: iLayerNrgFlux(:)                  ! intent(out): energy flux at the layer interfaces (W m-2)
   real(rkind), allocatable :: dNrgFlux_dTempAbove(:)            ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
   real(rkind), allocatable :: dNrgFlux_dTempBelow(:)            ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
   real(rkind), allocatable :: dNrgFlux_dWatAbove(:)             ! intent(out): derivatives in the flux w.r.t. water state in the layer above (J m-2 s-1 K-1)
   real(rkind), allocatable :: dNrgFlux_dWatBelow(:)             ! intent(out): derivatives in the flux w.r.t. water state in the layer below (J m-2 s-1 K-1)
   integer(i4b)             :: err                               ! intent(out): error code
   character(len=len_msg)   :: cmessage                          ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_snowSoilNrgFlux
 end type out_type_snowSoilNrgFlux
 ! ** end snowSoilNrgFlux

 ! ** vegLiqFlux
 type, public :: in_type_vegLiqFlux ! class for intent(in) arguments in vegLiqFlux call
   logical(lgt)             :: computeVegFlux                    ! intent(in): flag to denote if computing energy flux over vegetation
   real(rkind)              :: scalarCanopyLiqTrial              ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
   real(rkind)              :: scalarRainfall                    ! intent(in): rainfall rate (kg m-2 s-1)
  contains
   procedure :: initialize => initialize_in_vegLiqFlux
 end type in_type_vegLiqFlux

 type, public :: out_type_vegLiqFlux ! class for intent(out) arguments in vegLiqFlux call
   real(rkind)              :: scalarThroughfallRain             ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   real(rkind)              :: scalarCanopyLiqDrainage           ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   real(rkind)              :: scalarThroughfallRainDeriv        ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
   real(rkind)              :: scalarCanopyLiqDrainageDeriv      ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
   integer(i4b)             :: err                               ! intent(out): error code
   character(len=len_msg)   :: cmessage                          ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_vegLiqFlux
 end type out_type_vegLiqFlux
 ! ** end vegLiqFlux

 ! ** snowLiqFlux
 type, public :: in_type_snowLiqFlux ! class for intent(in) arguments in snowLiqFlux call
   integer(i4b)             :: nSnow                             ! intent(in):    number of snow layers
   logical(lgt)             :: firstFluxCall                     ! intent(in):    the first flux call (compute variables that are constant over the iterations)
   logical(lgt)             :: scalarSolution                    ! intent(in):    flag to indicate the scalar solution
   real(rkind)              :: scalarThroughfallRain             ! intent(in):    rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
   real(rkind)              :: scalarCanopyLiqDrainage           ! intent(in):    liquid drainage from the vegetation canopy (kg m-2 s-1)
   real(rkind), allocatable :: mLayerVolFracLiqTrial(:)          ! intent(in):    trial value of volumetric fraction of liquid water at the current iteration (-)
  contains
   procedure :: initialize => initialize_in_snowLiqFlux
 end type in_type_snowLiqFlux

 type, public :: io_type_snowLiqFlux ! class for intent(inout) arguments in snowLiqFlux call
   real(rkind), allocatable :: iLayerLiqFluxSnow(:)              ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
   real(rkind), allocatable :: iLayerLiqFluxSnowDeriv(:)         ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
  contains
   procedure :: initialize => initialize_io_snowLiqFlux
   procedure :: finalize   => finalize_io_snowLiqFlux
 end type io_type_snowLiqFlux

 type, public :: out_type_snowLiqFlux ! class for intent(out) arguments in snowLiqFlux call
   integer(i4b)             :: err                               ! intent(out):   error code
   character(len=len_msg)   :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize => finalize_out_snowLiqFlux
 end type out_type_snowLiqFlux
 ! ** end snowLiqFlux

 ! ** soilLiqFlux
 type, public :: in_type_soilLiqFlux ! class for intent(in) arguments in soilLiqFlux call
   integer(i4b)             :: nSoil                             ! intent(in):    number of soil layers
   logical(lgt)             :: firstSplitOper                    ! intent(in):    flag indicating first flux call in a splitting operation
   logical(lgt)             :: scalarSolution                    ! intent(in):    flag to indicate the scalar solution
   real(rkind)              :: scalarAquiferStorageTrial         ! intent(in):    trial value of aquifer storage (m)
   real(rkind), allocatable :: mLayerTempTrial(:)                ! intent(in):    trial temperature at the current iteration (K)
   real(rkind), allocatable :: mLayerMatricHeadTrial(:)          ! intent(in):    matric potential (m)
   real(rkind), allocatable :: mLayerMatricHeadLiqTrial(:)       ! intent(in):    liquid water matric potential (m)
   real(rkind), allocatable :: mLayerVolFracLiqTrial(:)          ! intent(in):    volumetric fraction of liquid water (-)
   real(rkind), allocatable :: mLayerVolFracIceTrial(:)          ! intent(in):    volumetric fraction of ice (-)
   real(rkind), allocatable :: mLayerdTheta_dTk(:)               ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind), allocatable :: dPsiLiq_dTemp(:)                  ! intent(in):    derivative in liquid water matric potential w.r.t. temperature (m K-1)
   real(rkind)              :: dCanopyTrans_dCanWat              ! intent(in):    derivative in canopy transpiration w.r.t. canopy total water content (s-1)
   real(rkind)              :: dCanopyTrans_dTCanair             ! intent(in):    derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyTrans_dTCanopy             ! intent(in):    derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyTrans_dTGround             ! intent(in):    derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
   real(rkind)              :: scalarCanopyTranspiration         ! intent(in):    canopy transpiration (kg m-2 s-1)
   real(rkind)              :: scalarGroundEvaporation           ! intent(in):    ground evaporation (kg m-2 s-1)
   real(rkind)              :: scalarRainPlusMelt                ! intent(in):    rain plus melt (m s-1)
  contains
   procedure :: initialize => initialize_in_soilLiqFlux
 end type in_type_soilLiqFlux

 type, public :: io_type_soilLiqFlux ! class for intent(inout) arguments in soilLiqFlux call
   real(rkind)              :: scalarMaxInfilRate                ! intent(inout): maximum infiltration rate (m s-1)
   real(rkind)              :: scalarInfilArea                   ! intent(inout): fraction of area where water can infiltrate, may be frozen (-)
   real(rkind)              :: scalarSaturatedArea               ! intent(inout): fraction of area that is considered saturated (-)
   real(rkind)              :: scalarFrozenArea                  ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
   real(rkind)              :: scalarSoilControl                 ! intent(inout): soil control on infiltration for derivative
   real(rkind)              :: scalarSurfaceRunoff               ! intent(inout): surface runoff (m s-1)
   real(rkind)              :: scalarSurfaceRunoff_IE            ! intent(inout): infiltration excess surface runoff (m s-1)
   real(rkind)              :: scalarSurfaceRunoff_SE            ! intent(inout): saturation excess surface runoff (m s-1)
   real(rkind), allocatable :: mLayerdTheta_dPsi(:)              ! intent(inout): derivative in liquid water content w.r.t. matric potential (m-1)
   real(rkind), allocatable :: dHydCond_dMatric(:)               ! intent(inout): derivative in hydraulic conductivity w.r.t matric head (s-1)
   real(rkind)              :: scalarInfiltration                ! intent(inout): surface infiltration rate (m s-1)
   real(rkind), allocatable :: iLayerLiqFluxSoil(:)              ! intent(inout): liquid fluxes at layer interfaces (m s-1)
   real(rkind), allocatable :: mLayerTranspire(:)                ! intent(inout): transpiration loss from each soil layer (m s-1)
   real(rkind), allocatable :: mLayerHydCond(:)                  ! intent(inout): hydraulic conductivity in each layer (m s-1)
   real(rkind), allocatable :: dq_dHydStateAbove(:)              ! intent(inout): derivatives in the flux w.r.t. matric head in the layer above (s-1)
   real(rkind), allocatable :: dq_dHydStateBelow(:)              ! intent(inout): derivatives in the flux w.r.t. matric head in the layer below (s-1)
   real(rkind), allocatable :: dq_dHydStateLayerSurfVec(:)       ! intent(inout): derivative in surface infiltration w.r.t. hydrology state in every soil layer (m s-1 or s-1)
   real(rkind), allocatable :: dq_dNrgStateAbove(:)              ! intent(inout): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
   real(rkind), allocatable :: dq_dNrgStateBelow(:)              ! intent(inout): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
   real(rkind), allocatable :: dq_dNrgStateLayerSurfVec(:)       ! intent(inout): derivative in surface infiltration w.r.t. energy state in every soil layer (m s-1 K-1)
   real(rkind), allocatable :: mLayerdTrans_dTCanair(:)          ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
   real(rkind), allocatable :: mLayerdTrans_dTCanopy(:)          ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
   real(rkind), allocatable :: mLayerdTrans_dTGround(:)          ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. ground temperature
   real(rkind), allocatable :: mLayerdTrans_dCanWat(:)           ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy total water  
  contains
   procedure :: initialize => initialize_io_soilLiqFlux
   procedure :: finalize   => finalize_io_soilLiqFlux
 end type io_type_soilLiqFlux

 type, public :: out_type_soilLiqFlux ! class for intent(out) arguments in soilLiqFlux call
   integer(i4b)             :: err                               ! intent(out):   error code
   character(len=len_msg)   :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize   => finalize_out_soilLiqFlux
 end type out_type_soilLiqFlux
 ! ** end soilLiqFlux

 ! ** groundwatr
 type, public :: in_type_groundwatr  ! class for intent(in) arguments in groundwatr call
   integer(i4b)             :: nSnow                             ! intent(in):    number of snow layers
   integer(i4b)             :: nSoil                             ! intent(in):    number of soil layers
   integer(i4b)             :: nLayers                           ! intent(in):    total number of layers
   logical(lgt)             :: firstFluxCall                     ! intent(in):    logical flag to compute index of the lowest saturated layer
   real(rkind), allocatable :: dVolTot_dPsi0(:)                  ! intent(in):    derivative in total volumetric water content w.r.t. matric head (m-1)
   real(rkind), allocatable :: mLayerdTheta_dPsi(:)              ! intent(in):    derivative in liquid water content w.r.t. matric potential (m-1)
   real(rkind), allocatable :: mLayerdTheta_dTk(:)               ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind), allocatable :: mLayerVolFracLiqTrial(:)          ! intent(in):    volumetric fraction of liquid water (-)
   real(rkind), allocatable :: mLayerVolFracIceTrial(:)          ! intent(in):    volumetric fraction of ice (-)
  contains
   procedure :: initialize => initialize_in_groundwatr
 end type in_type_groundwatr

 type, public :: io_type_groundwatr  ! class for intent(io) arguments in groundwatr call
   integer(i4b)             :: ixSaturation                      ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
  contains
   procedure :: initialize => initialize_io_groundwatr
   procedure :: finalize   => finalize_io_groundwatr 
 end type io_type_groundwatr

 type, public :: out_type_groundwatr ! class for intent(out) arguments in groundwatr call
   real(rkind), allocatable :: mLayerBaseflow(:)                 ! intent(out):   baseflow from each soil layer (m s-1)
   real(rkind), allocatable :: dBaseflow_dWat(:,:)               ! intent(out):   derivative in baseflow w.r.t. soil water characteristic
   real(rkind), allocatable :: dBaseflow_dTk(:,:)                ! intent(out):   derivative in baseflow w.r.t. temperature (m s-1 K-1)
   integer(i4b)             :: err                               ! intent(out):   error code
   character(len=len_msg)   :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize => finalize_out_groundwatr
 end type out_type_groundwatr
 ! ** end groundwatr

 ! ** bigAquifer
 type, public :: in_type_bigAquifer  ! class for intent(in) arguments in bigAquifer call
   real(rkind)              :: scalarAquiferStorageTrial         ! intent(in):    trial value of aquifer storage (m)
   real(rkind)              :: scalarCanopyTranspiration         ! intent(in):    canopy transpiration (kg m-2 s-1)
   real(rkind)              :: scalarSoilDrainage                ! intent(in):    soil drainage (m s-1)
   real(rkind)              :: dCanopyTrans_dCanWat              ! intent(in):    derivative in canopy transpiration w.r.t. canopy total water content (s-1)
   real(rkind)              :: dCanopyTrans_dTCanair             ! intent(in):    derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyTrans_dTCanopy             ! intent(in):    derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
   real(rkind)              :: dCanopyTrans_dTGround             ! intent(in):    derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
  contains
   procedure :: initialize => initialize_in_bigAquifer
 end type in_type_bigAquifer

 type, public :: io_type_bigAquifer  ! class for intent(inout) arguments in bigAquifer call
   real(rkind)              :: dAquiferTrans_dTCanair            ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
   real(rkind)              :: dAquiferTrans_dTCanopy            ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
   real(rkind)              :: dAquiferTrans_dTGround            ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
   real(rkind)              :: dAquiferTrans_dCanWat             ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
  contains
   procedure :: initialize => initialize_io_bigAquifer
   procedure :: finalize   => finalize_io_bigAquifer
 end type io_type_bigAquifer

 type, public :: out_type_bigAquifer  ! class for intent(out) arguments in bigAquifer call
   real(rkind)              :: scalarAquiferTranspire            ! intent(out):   transpiration loss from the aquifer (m s-1)
   real(rkind)              :: scalarAquiferRecharge             ! intent(out):   recharge to the aquifer (m s-1)
   real(rkind)              :: scalarAquiferBaseflow             ! intent(out):   total baseflow from the aquifer (m s-1)
   real(rkind)              :: dBaseflow_dAquifer                ! intent(out):   change in baseflow flux w.r.t. aquifer storage (s-1)
   integer(i4b)             :: err                               ! intent(out):   error code
   character(len=len_msg)   :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize   => finalize_out_bigAquifer
 end type out_type_bigAquifer
 ! ** end bigAquifer

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subroutines in soilLiqFlux
 ! ***********************************************************************************************************

 ! ** diagv_node
 type, public :: in_type_diagv_node ! intent(in) data
   ! input: state and diagnostic variables
   real(rkind)            :: scalarMatricHeadLiqTrial  ! liquid matric head in each layer (m)
   real(rkind)            :: scalarVolFracLiqTrial     ! volumetric fraction of liquid water in a given layer (-)
   real(rkind)            :: scalarVolFracIceTrial     ! volumetric fraction of ice in a given layer (-)
   ! input: pre-computed derivatives
   real(rkind)            :: dTheta_dTk                ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind)            :: dPsiLiq_dTemp             ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: soil parameters
   real(rkind)            :: vGn_alpha                 ! van Genuchten "alpha" parameter (m-1)
   real(rkind)            :: vGn_n                     ! van Genuchten "n" parameter (-)
   real(rkind)            :: vGn_m                     ! van Genuchten "m" parameter (-)
   real(rkind)            :: mpExp                     ! empirical exponent in macropore flow equation (-)
   real(rkind)            :: theta_sat                 ! soil porosity (-)
   real(rkind)            :: theta_res                 ! soil residual volumetric water content (-)
   real(rkind)            :: theta_mp                  ! volumetric liquid water content when macropore flow begins (-)
   real(rkind)            :: f_impede                  ! ice impedence factor (-)
   ! input: saturated hydraulic conductivity
   real(rkind)            :: scalarSatHydCond          ! saturated hydraulic conductivity at the mid-point of a given layer (m s-1)
   real(rkind)            :: scalarSatHydCondMP        ! saturated hydraulic conductivity of macropores at the mid-point of a given layer (m s-1)
  contains
   procedure :: initialize => initialize_in_diagv_node
 end type in_type_diagv_node

 type, public :: out_type_diagv_node ! intent(out) data
   ! output: derivative in the soil water characteristic
   real(rkind)            :: scalardPsi_dTheta         ! derivative in the soil water characteristic
   real(rkind)            :: scalardTheta_dPsi         ! derivative in the soil water characteristic
   ! output: transmittance
   real(rkind)            :: scalarHydCond             ! hydraulic conductivity at layer mid-points (m s-1)
   real(rkind)            :: scalarDiffuse             ! diffusivity at layer mid-points (m2 s-1)
   real(rkind)            :: iceImpedeFac              ! ice impedence factor in each layer (-)
   ! output: transmittance derivatives
   real(rkind)            :: dHydCond_dMatric          ! derivative in hydraulic conductivity w.r.t matric head (s-1)
   real(rkind)            :: dHydCond_dTemp            ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! output: error control
   integer(i4b)           :: err                       ! error code
   character(len=len_msg) :: message                   ! error message
  contains
   procedure :: finalize   => finalize_out_diagv_node
 end type out_type_diagv_node
 ! ** end diagv_node

 ! ** surfaceFlux
 type, public :: in_type_surfaceFlux ! intent(in) data
   ! input: model control
   logical(lgt) :: firstSplitOper   ! flag indicating if desire to compute infiltration
   integer(i4b) :: bc_lower         ! index defining the type of lower boundary conditions
   integer(i4b) :: ixInfRateMax     ! index defining the maximum infiltration rate method (GreenAmpt or topmodel_GA)
   integer(i4b) :: surfRun_SE       ! index defining the saturation excess surface runoff method
   integer(i4b) :: ix_groundwatr    ! index defining the groundwater parameterization
   integer(i4b) :: bc_upper         ! index defining the type of boundary conditions
   integer(i4b) :: nRoots           ! number of layers that contain roots
   integer(i4b) :: ixIce            ! index of lowest ice layer
   integer(i4b) :: nSoil            ! number of soil layers
   ! input: state and diagnostic variables
   real(rkind),allocatable :: mLayerTemp(:)        ! temperature (K)
   real(rkind)             :: scalarMatricHeadLiq  ! liquid matric head in the upper-most soil layer (m)
   real(rkind),allocatable :: mLayerMatricHead(:)  ! matric head in each soil layer (m)
   real(rkind)             :: scalarVolFracLiq     ! volumetric liquid water content in the upper-most soil layer (-)
   real(rkind)             :: scalarTotalSoilLiq   ! total liquid water in the soil column (kg m-2)
   real(rkind),allocatable :: mLayerVolFracLiq(:)  ! volumetric liquid water content in each soil layer (-)
   real(rkind),allocatable :: mLayerVolFracIce(:)  ! volumetric ice content in each soil layer (-)
   ! input: pre-computed derivatives (all of these would need to be recomputed if wanted a numerical derivative)
   real(rkind),allocatable :: dTheta_dTk(:)        ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind),allocatable :: dTheta_dPsi(:)       ! derivative in liquid water content w.r.t. matric potential (m-1)
   ! input: depth of each soil layer (m)
   real(rkind),allocatable :: mLayerDepth(:)       ! depth of each soil layer (m)
   real(rkind),allocatable :: iLayerHeight(:)      ! height at the interface of each layer (m)
   ! input: diriclet boundary conditions
   real(rkind) :: upperBoundHead      ! upper boundary condition for matric head (m)
   ! input: flux at the upper boundary
   real(rkind) :: scalarRainPlusMelt  ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input: transmittance
   real(rkind) :: surfaceSatHydCond   ! saturated hydraulic conductivity at the surface (m s-1)
   real(rkind) :: dHydCond_dTemp      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   real(rkind) :: iceImpedeFac        ! ice impedence factor in the upper-most soil layer (-)
   ! input: soil parameters
   real(rkind) :: vGn_alpha           ! van Genuchten "alpha" parameter (m-1)
   real(rkind) :: vGn_n               ! van Genuchten "n" parameter (-)
   real(rkind) :: vGn_m               ! van Genuchten "m" parameter (-)
   real(rkind) :: theta_sat           ! soil porosity (-)
   real(rkind) :: theta_res           ! soil residual volumetric water content (-)
   real(rkind) :: qSurfScale          ! scaling factor in the surface runoff parameterization (-)
   real(rkind) :: zScale_TOPMODEL     ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
   real(rkind) :: rootingDepth        ! rooting depth (m)
   real(rkind) :: wettingFrontSuction ! Green-Ampt wetting front suction (m)
   real(rkind) :: soilIceScale        ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   real(rkind) :: soilIceCV           ! soil ice CV in Gamma distribution used to define frozen area (-)
   ! input: FUSE parameters
   real(rkind) :: FUSE_Ac_max   ! FUSE PRMS max saturated area
   real(rkind) :: FUSE_phi_tens ! FUSE PRMS tension fraction
   real(rkind) :: FUSE_b        ! FUSE ARNO/VIC exponent
   real(rkind) :: FUSE_lambda   ! FUSE TOPMODEL gamma distribution lambda parameter
   real(rkind) :: FUSE_chi      ! FUSE TOPMODEL chi   distribution lambda parameter
   real(rkind) :: FUSE_mu       ! FUSE TOPMODEL mu    distribution lambda parameter
   real(rkind) :: FUSE_n        ! FUSE TOPMODEL exponent
  contains
   procedure :: initialize => initialize_in_surfaceFlux
 end type in_type_surfaceFlux 

 type, public :: io_type_surfaceFlux    ! intent(inout) data
   ! input-output: hydraulic conductivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   real(rkind) :: surfaceHydCond            ! hydraulic conductivity (m s-1)
   ! input-output: surface runoff and infiltration flux (m s-1)
   real(rkind) :: xMaxInfilRate             ! maximum infiltration rate (m s-1)
   real(rkind) :: scalarInfilArea           ! fraction of area where water can infiltrate, may be frozen (-)
   real(rkind) :: scalarSaturatedArea       ! fraction of area that is considered saturated (-)
   real(rkind) :: scalarFrozenArea          ! fraction of area that is considered impermeable due to soil ice (-)
   real(rkind) :: scalarSoilControl         ! soil control on infiltration for derivative
   real(rkind) :: scalarSurfaceInfiltration ! surface infiltration rate (m s-1)
  contains
   procedure :: initialize => initialize_io_surfaceFlux
   procedure :: finalize   => finalize_io_surfaceFlux
 end type io_type_surfaceFlux
 
 type, public :: out_type_surfaceFlux ! intent(out) data
   ! output: runoff and infiltration
   real(rkind) :: scalarSurfaceRunoff       ! surface runoff (m s-1)
   real(rkind) :: scalarSurfaceRunoff_IE    ! infiltration excess surface runoff (m s-1)
   real(rkind) :: scalarSurfaceRunoff_SE    ! saturation excess surface runoff (m s-1)
   ! output: derivatives in surface infiltration w.r.t. ...
   real(rkind),allocatable :: dq_dHydStateVec(:) ! ... hydrology state in every soil layer (m s-1 or s-1)
   real(rkind),allocatable :: dq_dNrgStateVec(:) ! ... energy state in every soil layer (m s-1 K-1)
   ! output: error control
   integer(i4b)            :: err     ! error code
   character(len=len_msg)  :: message ! error message
  contains
   procedure :: finalize   => finalize_out_surfaceFlux
 end type out_type_surfaceFlux 
 ! ** end surfaceFlux

 ! ** iLayerFlux
 type, public :: in_type_iLayerFlux ! intent(in) data
   ! input: state variables
   real(rkind),allocatable :: nodeMatricHeadLiqTrial(:) ! liquid matric head at the soil nodes (m)
   ! input: model coordinate variables
   real(rkind),allocatable :: nodeHeight(:)             ! height at the mid-point of the lower layer (m)
   ! input: temperature derivatives
   real(rkind),allocatable :: dPsiLiq_dTemp(:)          ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   real(rkind),allocatable :: dHydCond_dTemp(:)         ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: transmittance
   real(rkind),allocatable :: nodeHydCondTrial(:)       ! hydraulic conductivity at layer mid-points (m s-1)
   ! input: transmittance derivatives
   real(rkind),allocatable :: dHydCond_dMatric(:)       ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
  contains
   procedure :: initialize => initialize_in_iLayerFlux
 end type in_type_iLayerFlux

 type, public :: out_type_iLayerFlux ! intent(out) data
   ! output: tranmsmittance at the layer interface (scalars)
   real(rkind) :: iLayerHydCond      ! hydraulic conductivity at the interface between layers (m s-1)
   ! output: vertical flux at the layer interface (scalars)
   real(rkind) :: iLayerLiqFluxSoil  ! vertical flux of liquid water at the layer interface (m s-1)
   ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
   real(rkind) :: dq_dHydStateAbove  ! derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
   real(rkind) :: dq_dHydStateBelow  ! derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1)
   ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
   real(rkind) :: dq_dNrgStateAbove  ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
   real(rkind) :: dq_dNrgStateBelow  ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
   ! output: error control
   integer(i4b)           :: err     ! error code
   character(len=len_msg) :: message ! error message
  contains
   procedure :: finalize => finalize_out_iLayerFlux
 end type out_type_iLayerFlux
 ! ** end iLayerFlux

 ! ** qDrainFlux
 type, public :: in_type_qDrainFlux ! intent(in) data
   ! input: model control
   integer(i4b) :: bc_lower                  ! index defining the type of boundary conditions
   ! input: state and diagnostic variables
   real(rkind)  :: nodeMatricHeadLiq         ! liquid matric head in the lowest unsaturated node (m)
   ! input: model coordinate variables
   real(rkind)  :: nodeDepth                 ! depth of the lowest unsaturated soil layer (m)
   real(rkind)  :: nodeHeight                ! height of the lowest unsaturated soil node (m)
   ! input: diriclet boundary conditions
   real(rkind)  :: lowerBoundHead            ! lower boundary condition for matric head (m)
   ! input: derivative in soil water characteristic
   real(rkind)  :: node_dPsiLiq_dTemp        ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: transmittance
   real(rkind)  :: surfaceSatHydCond         ! saturated hydraulic conductivity at the surface (m s-1)
   real(rkind)  :: bottomSatHydCond          ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind)  :: nodeHydCond               ! hydraulic conductivity at the node itself (m s-1)
   real(rkind)  :: iceImpedeFac              ! ice impedence factor in the upper-most soil layer (-)
   ! input: transmittance derivatives
   real(rkind)  :: dHydCond_dMatric          ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
   real(rkind)  :: dHydCond_dTemp            ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   real(rkind)  :: vGn_alpha                 ! van Genuchten "alpha" parameter (m-1)
   real(rkind)  :: vGn_n                     ! van Genuchten "n" parameter (-)
   real(rkind)  :: vGn_m                     ! van Genuchten "m" parameter (-)
   real(rkind)  :: theta_sat                 ! soil porosity (-)
   real(rkind)  :: theta_res                 ! soil residual volumetric water content (-)
   real(rkind)  :: kAnisotropic              ! anisotropy factor for lateral hydraulic conductivity (-)
   real(rkind)  :: zScale_TOPMODEL           ! scale factor for TOPMODEL-ish baseflow parameterization (m)
  contains
   procedure :: initialize => initialize_in_qDrainFlux
 end type in_type_qDrainFlux

 type, public :: out_type_qDrainFlux ! intent(out) data
   ! output: hydraulic conductivity at the bottom of the unsaturated zone
   real(rkind) :: bottomHydCond      ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   ! output: drainage flux from the bottom of the soil profile
   real(rkind) :: scalarDrainage     ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux
   real(rkind) :: dq_dHydStateUnsat  ! change in drainage flux w.r.t. change in state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind) :: dq_dNrgStateUnsat  ! change in drainage flux w.r.t. change in energy state variable in lowest unsaturated node (m s-1 K-1)
   ! output: error control
   integer(i4b)           :: err     ! error code
   character(len=len_msg) :: message ! error message
  contains
   procedure :: finalize => finalize_out_qDrainFlux
 end type out_type_qDrainFlux
 ! ** end qDrainFlux

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subroutines in opSplittin
 ! ***********************************************************************************************************
 ! ** stateFilter
 type, public :: out_type_stateFilter ! class for intent(out) arguments in stateFilter call
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: cmessage                    ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_stateFilter
 end type out_type_stateFilter
 ! ** end stateFilter

 ! ** indexSplit
 type, public :: in_type_indexSplit  ! class for intent(in) arguments in indexSplit call
   integer(i4b)             :: nSnow                       ! intent(in): number of snow layers
   integer(i4b)             :: nSoil                       ! intent(in): number of soil layers
   integer(i4b)             :: nLayers                     ! intent(in): total number of layers
   integer(i4b)             :: nSubset                     ! intent(in): number of states in the subset
  contains
   procedure :: initialize => initialize_in_indexSplit
 end type in_type_indexSplit

 type, public :: out_type_indexSplit ! class for intent(out) arguments in indexSplit call
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: cmessage                    ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_indexSplit
 end type out_type_indexSplit
 ! ** end indexSplit

 ! ** varSubstep
 type, public :: in_type_varSubstep  ! class for intent(in) arguments in varSubstep call
   real(rkind)              :: dt                          ! intent(in): time step (s)
   real(rkind)              :: dtInit                      ! intent(in): initial time step (seconds)
   real(rkind)              :: dt_min                      ! intent(in): minimum time step (seconds)
   real(rkind)              :: whole_step                  ! intent(in): length of whole step for surface drainage and average flux
   integer(i4b)             :: nSubset                     ! intent(in): total number of variables in the state subset
   logical(lgt)             :: doAdjustTemp                ! intent(in): flag to indicate if we adjust the temperature
   logical(lgt)             :: firstSubStep                ! intent(in): flag to denote first sub-step
   logical(lgt)             :: computeVegFlux              ! intent(in): flag to denote if computing energy flux over vegetation
   logical(lgt)             :: scalarSolution              ! intent(in): flag to denote computing the scalar solution
   integer(i4b)             :: iStateSplit                 ! intent(in): index of the layer in the splitting operation
   type(var_flagVec)        :: fluxMask                    ! intent(in): mask for the fluxes used in this given state subset
  contains
   procedure :: initialize => initialize_in_varSubstep
 end type in_type_varSubstep

 type, public :: io_type_varSubstep  ! class for intent(inout) arguments in varSubstep call
   logical(lgt)             :: firstFluxCall               ! intent(inout): flag to indicate if we are processing the first flux call
   type(var_ilength)        :: fluxCount                   ! intent(inout): number of times fluxes are updated (should equal nsubstep)
   integer(i4b)             :: ixSaturation                ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
  contains
   procedure :: initialize => initialize_io_varSubstep
   procedure :: finalize   => finalize_io_varSubstep
 end type io_type_varSubstep

 type, public :: out_type_varSubstep  ! class for intent(out) arguments in varSubstep call
   real(rkind)              :: dtMultiplier                ! intent(out): substep multiplier (-)
   integer(i4b)             :: nSubsteps                   ! intent(out): number of substeps taken for a given split
   logical(lgt)             :: failedMinimumStep           ! intent(out): flag for failed substeps
   logical(lgt)             :: reduceCoupledStep           ! intent(out): flag to reduce the length of the coupled step
   logical(lgt)             :: tooMuchMelt                 ! intent(out): flag to denote that ice is insufficient to support melt
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: cmessage                    ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_varSubstep
 end type out_type_varSubstep
 ! ** end varSubstep

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subroutines in summaSolv4homegrown
 ! ***********************************************************************************************************

 type, public :: in_type_computJacob  ! class for intent(in) arguments in computJacob call
   ! input: model control
   real(rkind)              :: dt                          ! intent(in): length of the time step (seconds)
   integer(i4b)             :: nSnow                       ! intent(in): number of snow layers
   integer(i4b)             :: nSoil                       ! intent(in): number of soil layers
   integer(i4b)             :: nLayers                     ! intent(in): total number of layers in the snow+soil domain
   logical(lgt)             :: computeVegFlux              ! intent(in): flag to indicate if computing fluxes over vegetation
   logical(lgt)             :: computeBaseflow             ! intent(in): flag to indicate if computing baseflow
   integer(i4b)             :: ixMatrix                    ! intent(in): form of the Jacobian matrix
  contains
   procedure :: initialize => initialize_in_computJacob
 end type in_type_computJacob

 type, public :: out_type_computJacob  ! class for intent(out) arguments in computJacob call
   ! output: error control
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: cmessage                    ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_computJacob
 end type out_type_computJacob

 type, public :: in_type_lineSearchRefinement  ! class for intent(in) arguments in lineSearchRefinement call
   logical(lgt)             :: doSearch                    ! intent(in): flag to do the line search
   real(rkind)              :: fOld                        ! intent(in): old function value
  contains
   procedure :: initialize => initialize_in_lineSearchRefinement
 end type in_type_lineSearchRefinement

 type, public :: out_type_lineSearchRefinement  ! class for intent(out) arguments in lineSearchRefinement call
   real(rkind)              :: fNew                        ! intent(out): new function evaluation
   logical(lgt)             :: converged                   ! intent(out): convergence flag
   ! output: error control
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: message                     ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_lineSearchRefinement
 end type out_type_lineSearchRefinement

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subroutines in systemSolv
 ! ***********************************************************************************************************

 type, public :: in_type_summaSolv4homegrown  ! class for intent(in) arguments in summaSolv4homegrown call
   real(rkind)              :: dt_cur                   ! intent(in): current stepsize
   real(rkind)              :: dt                       ! intent(in): entire time step for drainage pond rate
   integer(i4b)             :: iter                     ! intent(in): iteration index
   integer(i4b)             :: nSnow                    ! intent(in): number of snow layers
   integer(i4b)             :: nSoil                    ! intent(in): number of soil layers
   integer(i4b)             :: nLayers                  ! intent(in): total number of layers
   integer(i4b)             :: nLeadDim                 ! intent(in): length of the leading dimension of the Jacobian matrix (nBands or nState)
   integer(i4b)             :: nState                   ! intent(in): total number of state variables
   integer(i4b)             :: ixMatrix                 ! intent(in): type of matrix (full or band diagonal)
   logical(lgt)             :: firstSubStep             ! intent(in): flag to indicate if we are processing the first sub-step
   logical(lgt)             :: computeVegFlux           ! intent(in): flag to indicate if computing fluxes over vegetation
   logical(lgt)             :: scalarSolution           ! intent(in): flag to denote if implementing the scalar solution
   real(rkind)              :: fOld                     ! intent(in): old function evaluation
  contains
   procedure :: initialize => initialize_in_summaSolv4homegrown
 end type in_type_summaSolv4homegrown

 type, public :: io_type_summaSolv4homegrown  ! class for intent(inout) arguments in summaSolv4homegrown call
   logical(lgt)             :: firstFluxCall            ! intent(inout): flag to indicate if we are processing the first flux call
   real(rkind)              :: xMin,xMax                ! intent(inout): brackets of the root
   integer(i4b)             :: ixSaturation             ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
  contains
   procedure :: initialize => initialize_io_summaSolv4homegrown
   procedure :: finalize   => finalize_io_summaSolv4homegrown
 end type io_type_summaSolv4homegrown

 type, public :: out_type_summaSolv4homegrown  ! class for intent(out) arguments in summaSolv4homegrown call
   real(rkind)              :: fNew                     ! intent(out): new function evaluation
   logical(lgt)             :: converged                ! intent(out): convergence flag
   integer(i4b)             :: err                      ! intent(out): error code
   character(len=len_msg)   :: message                  ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_summaSolv4homegrown
 end type out_type_summaSolv4homegrown

contains
 
 ! **** vegNrgFlux ****
 subroutine initialize_in_vegNrgFlux(in_vegNrgFlux,firstSubStep,firstFluxCall,computeVegFlux,checkLWBalance,&
                                     scalarCanairTempTrial,scalarCanopyTempTrial,mLayerTempTrial,scalarCanopyIceTrial,&
                                     scalarCanopyLiqTrial,forc_data,deriv_data)
  class(in_type_vegNrgFlux),intent(out) :: in_vegNrgFlux               ! class object for intent(in) vegNrgFlux arguments
  logical(lgt),intent(in)               :: firstSubStep                ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(in)               :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(in)               :: computeVegFlux              ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)               :: checkLWBalance              ! flag to check longwave balance
  real(rkind),intent(in)                :: scalarCanairTempTrial       ! trial value for temperature of the canopy air space (K)
  real(rkind),intent(in)                :: scalarCanopyTempTrial       ! trial value for temperature of the vegetation canopy (K)
  real(rkind),intent(in)                :: mLayerTempTrial(:)          ! trial value for temperature of each snow/soil layer (K)
  real(rkind),intent(in)                :: scalarCanopyIceTrial        ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),intent(in)                :: scalarCanopyLiqTrial        ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  type(var_d),intent(in)                :: forc_data                   ! model forcing data
  type(var_dlength),intent(in)          :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   upperBoundTemp               => forc_data%var(iLookFORCE%airtemp),                 & ! intent(in): [dp]     temperature of the upper boundary of the snow and soil domains (K)
   dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy)%dat(1) ) ! intent(out): [dp] derivative of canopy liquid storage w.r.t. temperature
   ! intent(in) arguments
   in_vegNrgFlux % firstSubStep=firstSubStep                      ! intent(in): flag to indicate if we are processing the first sub-step
   in_vegNrgFlux % firstFluxCall=firstFluxCall                    ! intent(in): flag to indicate if we are processing the first flux call
   in_vegNrgFlux % computeVegFlux=computeVegFlux                  ! intent(in): flag to indicate if we need to compute fluxes over vegetation
   in_vegNrgFlux % checkLWBalance=checkLWBalance                  ! intent(in): flag to check longwave balance
   in_vegNrgFlux % upperBoundTemp=upperBoundTemp                  ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
   in_vegNrgFlux % scalarCanairTempTrial=scalarCanairTempTrial    ! intent(in): trial value of the canopy air space temperature (K)
   in_vegNrgFlux % scalarCanopyTempTrial=scalarCanopyTempTrial    ! intent(in): trial value of canopy temperature (K)
   in_vegNrgFlux % mLayerTempTrial_1=mLayerTempTrial(1)           ! intent(in): trial value of ground temperature (K)
   in_vegNrgFlux % scalarCanopyIceTrial=scalarCanopyIceTrial      ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
   in_vegNrgFlux % scalarCanopyLiqTrial=scalarCanopyLiqTrial      ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
   in_vegNrgFlux % dCanLiq_dTcanopy=dCanLiq_dTcanopy              ! intent(in): derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)        
  end associate
 end subroutine initialize_in_vegNrgFlux

 subroutine finalize_out_vegNrgFlux(out_vegNrgFlux,flux_data,deriv_data,err,cmessage)
  class(out_type_vegNrgFlux),intent(in) :: out_vegNrgFlux              ! class object for intent(out) vegNrgFlux arguments
  type(var_dlength),intent(inout)       :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(inout)       :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from snowSoilNrgFlux

  ! intent(out) arguments: evapotranspiration values and net energy fluxes
  associate(&
    scalarCanopyTranspiration    => flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1), & ! intent(out): [dp] canopy transpiration (kg m-2 s-1)
    scalarCanopyEvaporation      => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1),   & ! intent(out): [dp] canopy evaporation/condensation (kg m-2 s-1)
    scalarGroundEvaporation      => flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1),   & ! intent(out): [dp] ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
    scalarCanairNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1),    & ! intent(out): [dp] net energy flux for the canopy air space        (W m-2)
    scalarCanopyNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1),    & ! intent(out): [dp] net energy flux for the vegetation canopy       (W m-2)
    scalarGroundNetNrgFlux       => flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1)     ) ! intent(out): [dp] net energy flux for the ground surface          (W m-2)
   scalarCanopyTranspiration  =out_vegNrgFlux % scalarCanopyTranspiration   ! intent(out): canopy transpiration (kg m-2 s-1)
   scalarCanopyEvaporation    =out_vegNrgFlux % scalarCanopyEvaporation     ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
   scalarGroundEvaporation    =out_vegNrgFlux % scalarGroundEvaporation     ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
   scalarCanairNetNrgFlux     =out_vegNrgFlux % scalarCanairNetNrgFlux      ! intent(out): net energy flux for the canopy air space (W m-2)
   scalarCanopyNetNrgFlux     =out_vegNrgFlux % scalarCanopyNetNrgFlux      ! intent(out): net energy flux for the vegetation canopy (W m-2)
   scalarGroundNetNrgFlux     =out_vegNrgFlux % scalarGroundNetNrgFlux      ! intent(out): net energy flux for the ground surface (W m-2)
  end associate

  ! intent(out) arguments: net canopy flux derivatives
  associate(&
    dCanairNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanairTemp)%dat(1), & ! intent(out): [dp] derivative in net canopy air space flux w.r.t. canopy air temperature
    dCanairNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanopyTemp)%dat(1), & ! intent(out): [dp] derivative in net canopy air space flux w.r.t. canopy temperature
    dCanairNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dGroundTemp)%dat(1), & ! intent(out): [dp] derivative in net canopy air space flux w.r.t. ground temperature
    dCanopyNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanairTemp)%dat(1), & ! intent(out): [dp] derivative in net canopy flux w.r.t. canopy air temperature
    dCanopyNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanopyTemp)%dat(1), & ! intent(out): [dp] derivative in net canopy flux w.r.t. canopy temperature
    dCanopyNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dGroundTemp)%dat(1), & ! intent(out): [dp] derivative in net canopy flux w.r.t. ground temperature
    dGroundNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanairTemp)%dat(1), & ! intent(out): [dp] derivative in net ground flux w.r.t. canopy air temperature
    dGroundNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanopyTemp)%dat(1), & ! intent(out): [dp] derivative in net ground flux w.r.t. canopy temperature
    dGroundNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp)%dat(1)  ) ! intent(out): [dp] derivative in net ground flux w.r.t. ground temperature
   dCanairNetFlux_dCanairTemp =out_vegNrgFlux % dCanairNetFlux_dCanairTemp  ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
   dCanairNetFlux_dCanopyTemp =out_vegNrgFlux % dCanairNetFlux_dCanopyTemp  ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
   dCanairNetFlux_dGroundTemp =out_vegNrgFlux % dCanairNetFlux_dGroundTemp  ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
   dCanopyNetFlux_dCanairTemp =out_vegNrgFlux % dCanopyNetFlux_dCanairTemp  ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
   dCanopyNetFlux_dCanopyTemp =out_vegNrgFlux % dCanopyNetFlux_dCanopyTemp  ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
   dCanopyNetFlux_dGroundTemp =out_vegNrgFlux % dCanopyNetFlux_dGroundTemp  ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
   dGroundNetFlux_dCanairTemp =out_vegNrgFlux % dGroundNetFlux_dCanairTemp  ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
   dGroundNetFlux_dCanopyTemp =out_vegNrgFlux % dGroundNetFlux_dCanopyTemp  ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
   dGroundNetFlux_dGroundTemp =out_vegNrgFlux % dGroundNetFlux_dGroundTemp  ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  end associate

  ! intent(out) arguments: canopy evaporation derivatives
  associate(&
    dCanopyEvaporation_dCanWat   => deriv_data%var(iLookDERIV%dCanopyEvaporation_dCanWat)%dat(1),  & ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy total water content
    dCanopyEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanair)%dat(1), & ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy air temperature
    dCanopyEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanopy)%dat(1), & ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy temperature
    dCanopyEvaporation_dTGround  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTGround)%dat(1), & ! intent(out): [dp] derivative in canopy evaporation w.r.t. ground temperature
    dGroundEvaporation_dCanWat   => deriv_data%var(iLookDERIV%dGroundEvaporation_dCanWat)%dat(1),  & ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy total water content
    dGroundEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanair)%dat(1), & ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy air temperature
    dGroundEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanopy)%dat(1), & ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy temperature
    dGroundEvaporation_dTGround  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround)%dat(1)  ) ! intent(out): [dp] derivative in ground evaporation w.r.t. ground temperature
   dCanopyEvaporation_dCanWat  =out_vegNrgFlux % dCanopyEvaporation_dCanWat  ! intent(out): derivative in canopy evaporation w.r.t. canopy total water content (s-1)
   dCanopyEvaporation_dTCanair =out_vegNrgFlux % dCanopyEvaporation_dTCanair ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dCanopyEvaporation_dTCanopy =out_vegNrgFlux % dCanopyEvaporation_dTCanopy ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dCanopyEvaporation_dTGround =out_vegNrgFlux % dCanopyEvaporation_dTGround ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
   dGroundEvaporation_dCanWat  =out_vegNrgFlux % dGroundEvaporation_dCanWat  ! intent(out): derivative in ground evaporation w.r.t. canopy total water content (s-1)
   dGroundEvaporation_dTCanair =out_vegNrgFlux % dGroundEvaporation_dTCanair ! intent(out): derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dGroundEvaporation_dTCanopy =out_vegNrgFlux % dGroundEvaporation_dTCanopy ! intent(out): derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dGroundEvaporation_dTGround =out_vegNrgFlux % dGroundEvaporation_dTGround ! intent(out): derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  end associate

  ! intent(out) arguments: canopy transpiration and net flux derivatives
  associate(& 
    dCanopyTrans_dCanWat     => deriv_data%var(iLookDERIV%dCanopyTrans_dCanWat)%dat(1),   & ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
    dCanopyTrans_dTCanair    => deriv_data%var(iLookDERIV%dCanopyTrans_dTCanair)%dat(1),  & ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dTCanopy    => deriv_data%var(iLookDERIV%dCanopyTrans_dTCanopy)%dat(1),  & ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dTGround    => deriv_data%var(iLookDERIV%dCanopyTrans_dTGround)%dat(1),  & ! intent(out): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
    dCanopyNetFlux_dCanWat   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanWat)%dat(1), & ! intent(out): [dp] derivative in net canopy fluxes w.r.t. canopy total water content
    dGroundNetFlux_dCanWat   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanWat)%dat(1)  ) ! intent(out): [dp] derivative in net ground fluxes w.r.t. canopy total water content
   dCanopyTrans_dCanWat       =out_vegNrgFlux % dCanopyTrans_dCanWat  ! intent(out): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
   dCanopyTrans_dTCanair      =out_vegNrgFlux % dCanopyTrans_dTCanair ! intent(out): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTCanopy      =out_vegNrgFlux % dCanopyTrans_dTCanopy ! intent(out): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTGround      =out_vegNrgFlux % dCanopyTrans_dTGround ! intent(out): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
   dCanopyNetFlux_dCanWat     =out_vegNrgFlux % dCanopyNetFlux_dCanWat! intent(out): derivative in net canopy fluxes w.r.t. canopy total water content (J kg-1 s-1)
   dGroundNetFlux_dCanWat     =out_vegNrgFlux % dGroundNetFlux_dCanWat! intent(out): derivative in net ground fluxes w.r.t. canopy total water content (J kg-1 s-1)
  end associate

   ! intent(out) arguments: error control
   err                        =out_vegNrgFlux % err                   ! intent(out): error code
   cmessage                   =out_vegNrgFlux % cmessage              ! intent(out): error message
 end subroutine finalize_out_vegNrgFlux
 ! **** end vegNrgFlux ****

 ! **** snowSoilNrgFlux ****
 subroutine initialize_in_snowSoilNrgFlux(in_snowSoilNrgFlux,scalarSolution,firstFluxCall,mLayerTempTrial,flux_data,deriv_data)
  class(in_type_snowSoilNrgFlux),intent(out) :: in_snowSoilNrgFlux               ! class object for intent(in) snowSoilNrgFlux arguments
  logical(lgt),intent(in)               :: scalarSolution              ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)               :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  real(rkind),intent(in)                :: mLayerTempTrial(:)          ! trial value for temperature of each snow/soil layer (K)
  type(var_dlength),intent(in)          :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(in)          :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   scalarGroundNetNrgFlux   => flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1), & ! intent(out): [dp] net energy flux for the ground surface (W m-2)
   iLayerLiqFluxSnow        => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat,         & ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   iLayerLiqFluxSoil        => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat,         & ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
   dThermalC_dWatAbove      => deriv_data%var(iLookDERIV%dThermalC_dWatAbove)%dat,     & ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. water state in the layer above
   dThermalC_dWatBelow      => deriv_data%var(iLookDERIV%dThermalC_dWatBelow)%dat,     & ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. water state in the layer above
   dThermalC_dTempAbove     => deriv_data%var(iLookDERIV%dThermalC_dTempAbove)%dat,    & ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. energy state in the layer above
   dThermalC_dTempBelow     => deriv_data%var(iLookDERIV%dThermalC_dTempBelow)%dat     ) ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. energy state in the layer above
   ! intent(in) arguments
   in_snowSoilNrgFlux % scalarSolution = scalarSolution .and. .not.firstFluxCall ! intent(in): flag to denote if implementing the scalar solution
   in_snowSoilNrgFlux % scalarGroundNetNrgFlux=scalarGroundNetNrgFlux            ! intent(in): net energy flux for the ground surface (W m-2)
   in_snowSoilNrgFlux % iLayerLiqFluxSnow     =iLayerLiqFluxSnow                 ! intent(in): liquid flux at the interface of each snow layer (m s-1)
   in_snowSoilNrgFlux % iLayerLiqFluxSoil     =iLayerLiqFluxSoil                 ! intent(in): liquid flux at the interface of each soil layer (m s-1)
   in_snowSoilNrgFlux % mLayerTempTrial       =mLayerTempTrial                   ! intent(in): temperature in each layer at the current iteration (m)
   in_snowSoilNrgFlux % dThermalC_dWatAbove   =dThermalC_dWatAbove               ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
   in_snowSoilNrgFlux % dThermalC_dWatBelow   =dThermalC_dWatBelow               ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
   in_snowSoilNrgFlux % dThermalC_dTempAbove  =dThermalC_dTempAbove              ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
   in_snowSoilNrgFlux % dThermalC_dTempBelow  =dThermalC_dTempBelow              ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
  end associate
 end subroutine initialize_in_snowSoilNrgFlux

 subroutine initialize_io_snowSoilNrgFlux(io_snowSoilNrgFlux,deriv_data)
  class(io_type_snowSoilNrgFlux),intent(out) :: io_snowSoilNrgFlux                 ! class object for intent(inout) snowSoilNrgFlux arguments
  type(var_dlength),intent(in)          :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   dGroundNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp)%dat(1) ) ! intent(out): [dp] derivative in net ground flux w.r.t. ground temperature
   ! intent(inout) arguments
   io_snowSoilNrgFlux % dGroundNetFlux_dGroundTemp=dGroundNetFlux_dGroundTemp ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  end associate
 end subroutine initialize_io_snowSoilNrgFlux

 subroutine finalize_io_snowSoilNrgFlux(io_snowSoilNrgFlux,deriv_data)
  class(io_type_snowSoilNrgFlux),intent(in)  :: io_snowSoilNrgFlux                 ! class object for intent(inout) snowSoilNrgFlux arguments
  type(var_dlength),intent(inout)       :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   dGroundNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp)%dat(1) ) ! intent(out): [dp] derivative in net ground flux w.r.t. ground temperature
   ! intent(inout) arguments
   dGroundNetFlux_dGroundTemp=io_snowSoilNrgFlux % dGroundNetFlux_dGroundTemp ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  end associate
 end subroutine finalize_io_snowSoilNrgFlux

 subroutine finalize_out_snowSoilNrgFlux(out_snowSoilNrgFlux,flux_data,deriv_data,err,cmessage)
  class(out_type_snowSoilNrgFlux),intent(in) :: out_snowSoilNrgFlux              ! class object for intent(out) snowSoilNrgFlux arguments
  type(var_dlength),intent(inout)       :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(inout)       :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from snowSoilNrgFlux
  associate(&
   iLayerNrgFlux         => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat,         & ! intent(out): [dp(0:)] vertical energy flux at the interface of snow and soil layers
   dNrgFlux_dTempAbove   => deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove)%dat, & ! intent(out): [dp(:)] derivatives in the flux w.r.t. temperature in the layer above
   dNrgFlux_dTempBelow   => deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow)%dat, & ! intent(out): [dp(:)] derivatives in the flux w.r.t. temperature in the layer below
   dNrgFlux_dWatAbove    => deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove)%dat,  & ! intent(out):  [dp(:)] derivatives in the flux w.r.t. water state in the layer above
   dNrgFlux_dWatBelow    => deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow)%dat   ) ! intent(out): [dp(:)] derivatives in the flux w.r.t. water state in the layer below
   ! intent(out) arguments
   iLayerNrgFlux       =out_snowSoilNrgFlux % iLayerNrgFlux          ! intent(out): energy flux at the layer interfaces (W m-2)
   dNrgFlux_dTempAbove =out_snowSoilNrgFlux % dNrgFlux_dTempAbove    ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
   dNrgFlux_dTempBelow =out_snowSoilNrgFlux % dNrgFlux_dTempBelow    ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
   dNrgFlux_dWatAbove  =out_snowSoilNrgFlux % dNrgFlux_dWatAbove     ! intent(out): derivatives in the flux w.r.t. water state in the layer above (J m-2 s-1 K-1)
   dNrgFlux_dWatBelow  =out_snowSoilNrgFlux % dNrgFlux_dWatBelow     ! intent(out): derivatives in the flux w.r.t. water state in the layer below (J m-2 s-1 K-1)
   err                 =out_snowSoilNrgFlux % err                    ! intent(out): error code
   cmessage            =out_snowSoilNrgFlux % cmessage               ! intent(out): error message
  end associate
 end subroutine finalize_out_snowSoilNrgFlux
 ! **** end snowSoilNrgFlux ****
 
 ! **** vegLiqFlux ****
 subroutine initialize_in_vegLiqFlux(in_vegLiqFlux,computeVegFlux,scalarCanopyLiqTrial,flux_data)
  class(in_type_vegLiqFlux),intent(out)   :: in_vegLiqFlux               ! class object for intent(in) vegLiqFlux arguments
  logical(lgt),intent(in)                 :: computeVegFlux              ! flag to indicate if computing fluxes over vegetation
  real(rkind),intent(in)                  :: scalarCanopyLiqTrial        ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  type(var_dlength),intent(in)            :: flux_data                   ! model fluxes for a local HRU
  associate(scalarRainfall => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)) ! intent(in): [dp] rainfall rate (kg m-2 s-1)
  ! intent(in) arguments
  in_vegLiqFlux % computeVegFlux      =computeVegFlux        ! intent(in): flag to denote if computing energy flux over vegetation
  in_vegLiqFlux % scalarCanopyLiqTrial=scalarCanopyLiqTrial  ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
  in_vegLiqFlux % scalarRainfall      =scalarRainfall        ! intent(in): rainfall rate (kg m-2 s-1)
  end associate
 end subroutine initialize_in_vegLiqFlux

 subroutine finalize_out_vegLiqFlux(out_vegLiqFlux,flux_data,deriv_data,err,cmessage)
  class(out_type_vegLiqFlux),intent(in)   :: out_vegLiqFlux              ! class object for intent(out) vegLiqFlux arguments
  type(var_dlength),intent(inout)         :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(inout)         :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  integer(i4b),intent(out)                :: err                         ! error code
  character(*),intent(out)                :: cmessage                    ! error message from vegLiqFlux
  associate( &
   scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),         & ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1),       & ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarThroughfallRainDeriv   => deriv_data%var(iLookDERIV%scalarThroughfallRainDeriv  )%dat(1),& ! intent(out): [dp] derivative in throughfall w.r.t. canopy liquid water
   scalarCanopyLiqDrainageDeriv => deriv_data%var(iLookDERIV%scalarCanopyLiqDrainageDeriv)%dat(1) ) ! intent(out): [dp] derivative in canopy drainage w.r.t. canopy liquid water
   ! intent(out) arguments
   scalarThroughfallRain       =out_vegLiqFlux % scalarThroughfallRain       ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage     =out_vegLiqFlux % scalarCanopyLiqDrainage     ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarThroughfallRainDeriv  =out_vegLiqFlux % scalarThroughfallRainDeriv  ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
   scalarCanopyLiqDrainageDeriv=out_vegLiqFlux % scalarCanopyLiqDrainageDeriv! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
   err                         =out_vegLiqFlux % err                         ! intent(out): error code
   cmessage                    =out_vegLiqFlux % cmessage                    ! intent(out): error control
  end associate
 end subroutine finalize_out_vegLiqFlux
 ! **** end vegLiqFlux ****

 ! **** snowLiqFlux ****
 subroutine initialize_in_snowLiqFlux(in_snowLiqFlux,nSnow,firstFluxCall,scalarSolution,mLayerVolFracLiqTrial,flux_data)
  class(in_type_snowLiqFlux),intent(out)  :: in_snowLiqFlux              ! class object for intent(in) snowLiqFlux arguments            
  integer(i4b),intent(in)                 :: nSnow                       ! number of snow layers
  logical(lgt),intent(in)                 :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(in)                 :: scalarSolution              ! flag to denote if implementing the scalar solution
  real(rkind),intent(in)                  :: mLayerVolFracLiqTrial(:)    ! trial value for volumetric fraction of liquid water (-)
  type(var_dlength),intent(in)            :: flux_data                   ! model fluxes for a local HRU
  associate(&
   scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),  & ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1))  ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  ! intent(in) arguments
  in_snowLiqFlux % nSnow                  =nSnow                          ! intent(in): number of snow layers
  in_snowLiqFlux % firstFluxCall          =firstFluxCall                  ! intent(in): the first flux call (compute variables that are constant over the iterations)
  in_snowLiqFlux % scalarSolution         =(scalarSolution .and. .not.firstFluxCall) ! intent(in): flag to indicate the scalar solution
  in_snowLiqFlux % scalarThroughfallRain  =scalarThroughfallRain          ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
  in_snowLiqFlux % scalarCanopyLiqDrainage=scalarCanopyLiqDrainage        ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
  in_snowLiqFlux % mLayerVolFracLiqTrial  =mLayerVolFracLiqTrial(1:nSnow) ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
  end associate
 end subroutine initialize_in_snowLiqFlux 

 subroutine initialize_io_snowLiqFlux(io_snowLiqFlux,flux_data,deriv_data)
  class(io_type_snowLiqFlux),intent(out)  :: io_snowLiqFlux               ! class object for intent(inout) snowLiqFlux arguments
  type(var_dlength),intent(in)            :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(in)            :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
    iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat,       & ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
    iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv)%dat ) ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
  io_snowLiqFlux % iLayerLiqFluxSnow      =iLayerLiqFluxSnow       ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
  io_snowLiqFlux % iLayerLiqFluxSnowDeriv =iLayerLiqFluxSnowDeriv  ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
  end associate
 end subroutine initialize_io_snowLiqFlux

 subroutine finalize_io_snowLiqFlux(io_snowLiqFlux,flux_data,deriv_data)
  class(io_type_snowLiqFlux),intent(in)   :: io_snowLiqFlux              ! class object for intent(inout) snowLiqFlux arguments
  type(var_dlength),intent(inout)         :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(inout)         :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
    iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat,       & ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
    iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv)%dat ) ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
  ! intent(inout) arguments
  iLayerLiqFluxSnow     =io_snowLiqFlux % iLayerLiqFluxSnow               ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
  iLayerLiqFluxSnowDeriv=io_snowLiqFlux % iLayerLiqFluxSnowDeriv          ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
  end associate
 end subroutine finalize_io_snowLiqFlux

 subroutine finalize_out_snowLiqFlux(out_snowLiqFlux,err,cmessage)
  class(out_type_snowLiqFlux),intent(in)  :: out_snowLiqFlux             ! class object for intent(out) snowLiqFlux arguments
  integer(i4b),intent(out)                :: err                         ! error code
  character(*),intent(out)                :: cmessage                    ! error message from snowLiqFlux
  ! intent(out) arguments
  err     =out_snowLiqFlux % err        ! intent(out):   error code
  cmessage=out_snowLiqFlux % cmessage   ! intent(out):   error message
 end subroutine finalize_out_snowLiqFlux
 ! **** end snowLiqFlux ****

 ! **** soilLiqFlux ****
 subroutine initialize_in_soilLiqFlux(in_soilLiqFlux,nSnow,nSoil,nlayers,firstSplitOper,scalarSolution,firstFluxCall,scalarAquiferStorageTrial,&
                                     mLayerTempTrial,mLayerMatricHeadTrial,mLayerMatricHeadLiqTrial,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,&
                                     flux_data,deriv_data)
  class(in_type_soilLiqFlux),intent(out) :: in_soilLiqFlux              ! class object for intent(in) soilLiqFlux arguments
  integer(i4b),intent(in)                :: nSnow                       ! number of snow layers
  integer(i4b),intent(in)                :: nSoil                       ! number of soil layers
  integer(i4b),intent(in)                :: nLayers                     ! total number of layers
  logical(lgt),intent(in)                :: firstSplitOper              ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)                :: scalarSolution              ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)                :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  real(rkind),intent(in)                 :: scalarAquiferStorageTrial   ! trial value of aquifer storage (m)
  real(rkind),intent(in)                 :: mLayerTempTrial(:)          ! trial value for temperature of each snow/soil layer (K)
  real(rkind),intent(in)                 :: mLayerMatricHeadTrial(:)    ! trial value for the total water matric potential (m)
  real(rkind),intent(in)                 :: mLayerMatricHeadLiqTrial(:) ! trial value for the liquid water matric potential (m)
  real(rkind),intent(in)                 :: mLayerVolFracLiqTrial(:)    ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in)                 :: mLayerVolFracIceTrial(:)    ! trial value for volumetric fraction of ice (-)
  type(var_dlength),intent(in)           :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(in)           :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables

  ! intent(in) arguments: model control
  in_soilLiqFlux % nSoil         =nSoil                                         ! intent(in): number of soil layers
  in_soilLiqFlux % firstSplitOper=firstSplitOper                                ! intent(in): flag indicating first flux call in a splitting operation
  in_soilLiqFlux % scalarSolution=(scalarSolution .and. .not.firstFluxCall)     ! intent(in): flag to indicate the scalar solution

  ! intent(in) arguments: aquifer variables needed for FUSE parameterizations
  in_soilLiqFlux % scalarAquiferStorageTrial = scalarAquiferStorageTrial        ! intent(in): trial value of aquifer storage (m)

  ! intent(in) arguments: trial temperature, matric potential, and volumetric fractions
  in_soilLiqFlux % mLayerTempTrial=mLayerTempTrial(nSnow+1:nLayers)             ! intent(in): trial temperature at the current iteration (K)
  in_soilLiqFlux % mLayerMatricHeadTrial   =mLayerMatricHeadTrial(1:nSoil)      ! intent(in): matric potential (m)
  in_soilLiqFlux % mLayerMatricHeadLiqTrial=mLayerMatricHeadLiqTrial(1:nSoil)   ! intent(in): liquid water matric potential (m)
  in_soilLiqFlux % mLayerVolFracLiqTrial=mLayerVolFracLiqTrial(nSnow+1:nLayers) ! intent(in): volumetric fraction of liquid water (-)
  in_soilLiqFlux % mLayerVolFracIceTrial=mLayerVolFracIceTrial(nSnow+1:nLayers) ! intent(in): volumetric fraction of ice (-)

  ! intent(in) arguments: derivatives for liquid water
  associate(&
   mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat, & ! intent(in): [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
   dPsiLiq_dTemp                => deriv_data%var(iLookDERIV%dPsiLiq_dTemp)%dat     ) ! intent(in): [dp(:)] derivative in the liquid water matric potential w.r.t. temperature
   in_soilLiqFlux % mLayerdTheta_dTk=mLayerdTheta_dTk(nSnow+1:nLayers)           ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
   in_soilLiqFlux % dPsiLiq_dTemp=dPsiLiq_dTemp(1:nSoil)                         ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
  end associate

   ! intent(in) arguments: canopy transpiration derivatives
  associate(&
   dCanopyTrans_dCanWat         => deriv_data%var(iLookDERIV%dCanopyTrans_dCanWat)%dat(1),  & ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
   dCanopyTrans_dTCanair        => deriv_data%var(iLookDERIV%dCanopyTrans_dTCanair)%dat(1), & ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTCanopy        => deriv_data%var(iLookDERIV%dCanopyTrans_dTCanopy)%dat(1), & ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTGround        => deriv_data%var(iLookDERIV%dCanopyTrans_dTGround)%dat(1)  ) ! intent(out): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
   in_soilLiqFlux % dCanopyTrans_dCanWat  =dCanopyTrans_dCanWat     ! intent(in): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
   in_soilLiqFlux % dCanopyTrans_dTCanair =dCanopyTrans_dTCanair    ! intent(in): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   in_soilLiqFlux % dCanopyTrans_dTCanopy =dCanopyTrans_dTCanopy    ! intent(in): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
   in_soilLiqFlux % dCanopyTrans_dTGround =dCanopyTrans_dTGround    ! intent(in): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
  end associate

  ! intent(in) arguments: evaporative fluxes and rain plus melt
  associate(&
   scalarCanopyTranspiration    => flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1), & ! intent(out): [dp] canopy transpiration (kg m-2 s-1)
   scalarGroundEvaporation      => flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1),   & ! intent(out): [dp] ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
   scalarRainPlusMelt           => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1)         ) ! intent(out): [dp] rain plus melt (m s-1)
   in_soilLiqFlux % scalarCanopyTranspiration=scalarCanopyTranspiration                          ! intent(in): canopy transpiration (kg m-2 s-1)
   in_soilLiqFlux % scalarGroundEvaporation  =scalarGroundEvaporation                            ! intent(in): ground evaporation (kg m-2 s-1)
   in_soilLiqFlux % scalarRainPlusMelt       =scalarRainPlusMelt                                 ! intent(in): rain plus melt (m s-1)
  end associate
 end subroutine initialize_in_soilLiqFlux

 subroutine initialize_io_soilLiqFlux(io_soilLiqFlux,nSoil,dHydCond_dMatric,flux_data,diag_data,deriv_data)
  class(io_type_soilLiqFlux),intent(out) :: io_soilLiqFlux              ! class object for intent(inout) soilLiqFlux arguments
  integer(i4b),intent(in)                :: nSoil                       ! number of soil layers
  real(rkind),intent(in)                 :: dHydCond_dMatric(nSoil)     ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  type(var_dlength),intent(in)           :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(in)           :: diag_data                   ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)           :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables

  ! intent(inout) arguments: max infiltration rate, frozen area, and surface runoff
  associate(&
   scalarMaxInfilRate     => flux_data%var(iLookFLUX%scalarMaxInfilRate)%dat(1),     & ! intent(out): [dp] maximum infiltration rate (m s-1)
   scalarInfilArea        => diag_data%var(iLookDIAG%scalarInfilArea)%dat(1),        & ! intent(out): [dp] fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea    => diag_data%var(iLookDIAG%scalarSaturatedArea)%dat(1),    & ! intent(out): [dp] fraction of area that is considered saturated (-)
   scalarFrozenArea       => diag_data%var(iLookDIAG%scalarFrozenArea)%dat(1),       & ! intent(out): [dp] fraction of area that is considered impermeable due to soil ice (-)
   scalarSoilControl      => diag_data%var(iLookDIAG%scalarSoilControl)%dat(1),      & ! intent(out): [dp] soil control on infiltration for derivative
   scalarSurfaceRunoff    => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1),    & ! intent(out): [dp] surface runoff (m s-1)
   scalarSurfaceRunoff_IE => flux_data%var(iLookFLUX%scalarSurfaceRunoff_IE)%dat(1), & ! intent(out): [dp] infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE => flux_data%var(iLookFLUX%scalarSurfaceRunoff_SE)%dat(1)  ) ! intent(out): [dp] saturation excess surface runoff (m s-1)
   io_soilLiqFlux % scalarMaxInfilRate      =scalarMaxInfilRate       ! intent(inout): maximum infiltration rate (m s-1)
   io_soilLiqFlux % scalarInfilArea         =scalarInfilArea          ! intent(inout): fraction of area where water can infiltrate, may be frozen (-)
   io_soilLiqFlux % scalarSaturatedArea     =scalarSaturatedArea      ! intent(inout): fraction of area that is considered saturated (-)
   io_soilLiqFlux % scalarFrozenArea        =scalarFrozenArea         ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
   io_soilLiqFlux % scalarSoilControl       =scalarSoilControl        ! intent(inout): soil control on infiltration for derivative
   io_soilLiqFlux % scalarSurfaceRunoff     =scalarSurfaceRunoff      ! intent(inout): surface runoff (m s-1)
   io_soilLiqFlux % scalarSurfaceRunoff_IE  =scalarSurfaceRunoff_IE   ! intent(inout): infiltration excess surface runoff (m s-1)
   io_soilLiqFlux % scalarSurfaceRunoff_SE  =scalarSurfaceRunoff_SE   ! intent(inout): saturation excess surface runoff (m s-1)
  end associate

  ! intent(inout) arguments: derivatives, fluxes, and layer properties
  associate(& 
   mLayerdTheta_dPsi            => deriv_data%var(iLookDERIV%mLayerdTheta_dPsi)%dat,   & ! intent(out): [dp(:)] derivative in liquid water content w.r.t. matric potential (m-1)
   scalarInfiltration           => flux_data%var(iLookFLUX%scalarInfiltration)%dat(1), & ! intent(out): [dp] infiltration of water into the soil profile (m s-1)
   iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat,     & ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
   mLayerTranspire              => flux_data%var(iLookFLUX%mLayerTranspire)%dat,       & ! intent(out): [dp(:)] transpiration loss from each soil layer (m s-1)
   mLayerHydCond                => flux_data%var(iLookFLUX%mLayerHydCond)%dat          ) ! intent(out): [dp(:)]  hydraulic conductivity in each soil layer (m s-1)
   io_soilLiqFlux % mLayerdTheta_dPsi       =mLayerdTheta_dPsi        ! intent(inout): derivative in liquid water content w.r.t. matric potential (m-1)
   io_soilLiqFlux % dHydCond_dMatric        =dHydCond_dMatric         ! intent(inout): derivative in hydraulic conductivity w.r.t matric head (s-1)
   io_soilLiqFlux % scalarInfiltration      =scalarInfiltration       ! intent(inout): surface infiltration rate (m s-1)
   io_soilLiqFlux % iLayerLiqFluxSoil       =iLayerLiqFluxSoil        ! intent(inout): liquid fluxes at layer interfaces (m s-1)
   io_soilLiqFlux % mLayerTranspire         =mLayerTranspire          ! intent(inout): transpiration loss from each soil layer (m s-1)
   io_soilLiqFlux % mLayerHydCond           =mLayerHydCond            ! intent(inout): hydraulic conductivity in each layer (m s-1)
  end associate

  ! intent(inout) arguments: flux and surface infiltration derivatives
  associate(&
   dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
   dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
   dq_dHydStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec)%dat, & ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
   dq_dNrgStateAbove            => deriv_data%var(iLookDERIV%dq_dNrgStateAbove)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
   dq_dNrgStateBelow            => deriv_data%var(iLookDERIV%dq_dNrgStateBelow)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
   dq_dNrgStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec)%dat  ) ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
   io_soilLiqFlux % dq_dHydStateAbove       =dq_dHydStateAbove        ! intent(inout): derivatives in the flux w.r.t. matric head in the layer above (s-1)
   io_soilLiqFlux % dq_dHydStateBelow       =dq_dHydStateBelow        ! intent(inout): derivatives in the flux w.r.t. matric head in the layer below (s-1)
   io_soilLiqFlux % dq_dHydStateLayerSurfVec=dq_dHydStateLayerSurfVec ! intent(inout): derivative in surface infiltration w.r.t. hydrology state in every soil layer (m s-1 or s-1)
   io_soilLiqFlux % dq_dNrgStateAbove       =dq_dNrgStateAbove        ! intent(inout): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
   io_soilLiqFlux % dq_dNrgStateBelow       =dq_dNrgStateBelow        ! intent(inout): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
   io_soilLiqFlux % dq_dNrgStateLayerSurfVec=dq_dNrgStateLayerSurfVec ! intent(inout): derivative in surface infiltration w.r.t. energy state in every soil layer (m s-1 K-1)
  end associate

  ! intent(inout) arguments: transpiration flux derivatives
  associate(&
   mLayerdTrans_dTCanair        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanair)%dat,  & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
   mLayerdTrans_dTCanopy        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanopy)%dat, & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
   mLayerdTrans_dTGround        => deriv_data%var(iLookDERIV%mLayerdTrans_dTGround)%dat, & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. ground temperature
   mLayerdTrans_dCanWat         => deriv_data%var(iLookDERIV%mLayerdTrans_dCanWat)%dat   ) ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy total water
   io_soilLiqFlux % mLayerdTrans_dTCanair   =mLayerdTrans_dTCanair    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
   io_soilLiqFlux % mLayerdTrans_dTCanopy   =mLayerdTrans_dTCanopy    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
   io_soilLiqFlux % mLayerdTrans_dTGround   =mLayerdTrans_dTGround    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. ground temperature
   io_soilLiqFlux % mLayerdTrans_dCanWat    =mLayerdTrans_dCanWat     ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy total water 
  end associate
 end subroutine initialize_io_soilLiqFlux

 subroutine finalize_io_soilLiqFlux(io_soilLiqFlux,nSoil,dHydCond_dMatric,flux_data,diag_data,deriv_data)
  class(io_type_soilLiqFlux),intent(in) :: io_soilLiqFlux             ! class object for intent(inout) soilLiqFlux arguments
  integer(i4b),intent(in)               :: nSoil                      ! number of soil layers
  real(rkind),intent(out)               :: dHydCond_dMatric(nSoil)    ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  type(var_dlength),intent(inout)       :: flux_data                  ! model fluxes for a local HRU
  type(var_dlength),intent(inout)       :: diag_data                  ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)       :: deriv_data                 ! derivatives in model fluxes w.r.t. relevant state variables

  ! intent(inout) arguments: max infiltration rate, frozen area, and surface runoff
  associate(&
   scalarMaxInfilRate     => flux_data%var(iLookFLUX%scalarMaxInfilRate)%dat(1),     & ! intent(out): [dp] maximum infiltration rate (m s-1)
   scalarInfilArea        => diag_data%var(iLookDIAG%scalarInfilArea)%dat(1),        & ! intent(out): [dp] fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea    => diag_data%var(iLookDIAG%scalarSaturatedArea)%dat(1),    & ! intent(out): [dp] fraction of area that is considered saturated (-)
   scalarFrozenArea       => diag_data%var(iLookDIAG%scalarFrozenArea)%dat(1),       & ! intent(out): [dp] fraction of area that is considered impermeable due to soil ice (-)
   scalarSoilControl      => diag_data%var(iLookDIAG%scalarSoilControl)%dat(1),      & ! intent(out): [dp] soil control on infiltration for derivative
   scalarSurfaceRunoff    => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1),    & ! intent(out): [dp] surface runoff (m s-1)
   scalarSurfaceRunoff_IE => flux_data%var(iLookFLUX%scalarSurfaceRunoff_IE)%dat(1), & ! intent(out): [dp] infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE => flux_data%var(iLookFLUX%scalarSurfaceRunoff_SE)%dat(1)  ) ! intent(out): [dp] saturation excess surface runoff (m s-1)
   scalarMaxInfilRate      =io_soilLiqFlux % scalarMaxInfilRate       ! intent(inout): maximum infiltration rate (m s-1)
   scalarInfilArea         =io_soilLiqFlux % scalarInfilArea          ! intent(inout): fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea     =io_soilLiqFlux % scalarSaturatedArea      ! intent(inout): fraction of area that is considered saturated (-)
   scalarFrozenArea        =io_soilLiqFlux % scalarFrozenArea         ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
   scalarSoilControl       =io_soilLiqFlux % scalarSoilControl        ! intent(inout): soil control on infiltration for derivative
   scalarSurfaceRunoff     =io_soilLiqFlux % scalarSurfaceRunoff      ! intent(inout): surface runoff (m s-1)
   scalarSurfaceRunoff_IE  =io_soilLiqFlux % scalarSurfaceRunoff_IE   ! intent(inout): infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE  =io_soilLiqFlux % scalarSurfaceRunoff_SE   ! intent(inout): saturation excess surface runoff (m s-1)
  end associate

  ! intent(inout) arguments: derivatives, fluxes, and layer properties
  associate(& 
   mLayerdTheta_dPsi      => deriv_data%var(iLookDERIV%mLayerdTheta_dPsi)%dat,   & ! intent(out): [dp(:)] derivative in liquid water content w.r.t. matric potential (m-1)
   scalarInfiltration     => flux_data%var(iLookFLUX%scalarInfiltration)%dat(1), & ! intent(out): [dp] infiltration of water into the soil profile (m s-1)
   iLayerLiqFluxSoil      => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat,     & ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
   mLayerTranspire        => flux_data%var(iLookFLUX%mLayerTranspire)%dat,       & ! intent(out): [dp(:)] transpiration loss from each soil layer (m s-1)
   mLayerHydCond          => flux_data%var(iLookFLUX%mLayerHydCond)%dat          ) ! intent(out): [dp(:)]  hydraulic conductivity in each soil layer (m s-1)
   mLayerdTheta_dPsi       =io_soilLiqFlux % mLayerdTheta_dPsi        ! intent(inout): derivative in liquid water content w.r.t. matric potential (m-1)
   dHydCond_dMatric        =io_soilLiqFlux % dHydCond_dMatric         ! intent(inout): derivative in hydraulic conductivity w.r.t matric head (s-1)
   scalarInfiltration      =io_soilLiqFlux % scalarInfiltration       ! intent(inout): surface infiltration rate (m s-1)
   iLayerLiqFluxSoil       =io_soilLiqFlux % iLayerLiqFluxSoil        ! intent(inout): liquid fluxes at layer interfaces (m s-1)
   mLayerTranspire         =io_soilLiqFlux % mLayerTranspire          ! intent(inout): transpiration loss from each soil layer (m s-1)
   mLayerHydCond           =io_soilLiqFlux % mLayerHydCond            ! intent(inout): hydraulic conductivity in each layer (m s-1)
  end associate

  ! intent(inout) arguments: flux and surface infiltration derivatives
  associate(&
   dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
   dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
   dq_dHydStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec)%dat, & ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
   dq_dNrgStateAbove            => deriv_data%var(iLookDERIV%dq_dNrgStateAbove)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
   dq_dNrgStateBelow            => deriv_data%var(iLookDERIV%dq_dNrgStateBelow)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
   dq_dNrgStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec)%dat  ) ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
   dq_dHydStateAbove       =io_soilLiqFlux % dq_dHydStateAbove        ! intent(inout): derivatives in the flux w.r.t. matric head in the layer above (s-1)
   dq_dHydStateBelow       =io_soilLiqFlux % dq_dHydStateBelow        ! intent(inout): derivatives in the flux w.r.t. matric head in the layer below (s-1)
   dq_dHydStateLayerSurfVec=io_soilLiqFlux % dq_dHydStateLayerSurfVec ! intent(inout): derivative in surface infiltration w.r.t. hydrology state in every soil layer (m s-1 or s-1)
   dq_dNrgStateAbove       =io_soilLiqFlux % dq_dNrgStateAbove        ! intent(inout): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
   dq_dNrgStateBelow       =io_soilLiqFlux % dq_dNrgStateBelow        ! intent(inout): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
   dq_dNrgStateLayerSurfVec=io_soilLiqFlux % dq_dNrgStateLayerSurfVec ! intent(inout): derivative in surface infiltration w.r.t. energy state in every soil layer (m s-1 K-1)
  end associate

  ! intent(inout) arguments: transpiration flux derivatives
  associate(&
   mLayerdTrans_dTCanair        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanair)%dat, & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
   mLayerdTrans_dTCanopy        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanopy)%dat, & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
   mLayerdTrans_dTGround        => deriv_data%var(iLookDERIV%mLayerdTrans_dTGround)%dat, & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. ground temperature
   mLayerdTrans_dCanWat         => deriv_data%var(iLookDERIV%mLayerdTrans_dCanWat)%dat   ) ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy total water
   mLayerdTrans_dTCanair   =io_soilLiqFlux % mLayerdTrans_dTCanair      ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
   mLayerdTrans_dTCanopy   =io_soilLiqFlux % mLayerdTrans_dTCanopy      ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
   mLayerdTrans_dTGround   =io_soilLiqFlux % mLayerdTrans_dTGround      ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. ground temperature
   mLayerdTrans_dCanWat    =io_soilLiqFlux % mLayerdTrans_dCanWat       ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy total water 
  end associate
 end subroutine finalize_io_soilLiqFlux

 subroutine finalize_out_soilLiqFlux(out_soilLiqFlux,err,cmessage)
  class(out_type_soilLiqFlux),intent(in) :: out_soilLiqFlux             ! class object for intent(out) soilLiqFlux arguments
  integer(i4b),intent(out)               :: err                         ! error code
  character(*),intent(out)               :: cmessage                    ! error message from groundwatr
  ! intent(out) arguments
  err                     =out_soilLiqFlux % err       ! intent(out):   error code
  cmessage                =out_soilLiqFlux % cmessage  ! intent(out):   error message
 end subroutine finalize_out_soilLiqFlux
 ! **** end soilLiqFlux ****

 ! **** groundwatr ****
 subroutine initialize_in_groundwatr(in_groundwatr,nSnow,nSoil,nLayers,firstFluxCall,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,deriv_data)
  class(in_type_groundwatr),intent(out) :: in_groundwatr               ! class object for intent(in) groundwatr arguments
  integer(i4b),intent(in)               :: nSnow                       ! number of snow layers
  integer(i4b),intent(in)               :: nSoil                       ! number of soil layers
  integer(i4b),intent(in)               :: nLayers                     ! total number of layers
  logical(lgt),intent(in)               :: firstFluxCall               ! logical flag to compute index of the lowest saturated layer
  real(rkind),intent(in)                :: mLayerVolFracLiqTrial(:)    ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in)                :: mLayerVolFracIceTrial(:)    ! trial value for volumetric fraction of ice (-)
  type(var_dlength),intent(in)          :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
 
  associate(&
   dVolTot_dPsi0     => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat,         & ! intent(in): [dp(:)] derivative in total volumetric water content w.r.t. matric head (m-1)
   mLayerdTheta_dPsi => deriv_data%var(iLookDERIV%mLayerdTheta_dPsi)%dat,     & ! intent(in): [dp(:)] derivative in liquid water content w.r.t. matric potential (m-1)
   mLayerdTheta_dTk  => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat(nSnow+1:nSnow+nSoil) )! intent(in): [dp(:)] derivative in volumetric liquid water content w.r.t. temperature (K-1)
   ! intent(in) arguments
   in_groundwatr % nSnow                    = nSnow                                  ! intent(in):    number of snow layers
   in_groundwatr % nSoil                    = nSoil                                  ! intent(in):    number of soil layers
   in_groundwatr % nLayers                  = nLayers                                ! intent(in):    total number of layers
   in_groundwatr % firstFluxCall            = firstFluxCall                          ! intent(in):    logical flag to compute index of the lowest saturated layer
   in_groundwatr % dVolTot_dPsi0            = dVolTot_dPsi0                          ! intent(in):    derivative in total volumetric water content w.r.t. matric head (m-1)
   in_groundwatr % mLayerdTheta_dPsi        = mLayerdTheta_dPsi                      ! intent(in):    derivative in liquid water content w.r.t. matric potential (m-1)
   in_groundwatr % mLayerdTheta_dTk         = mLayerdTheta_dTk                       ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
   in_groundwatr % mLayerVolFracLiqTrial    = mLayerVolFracLiqTrial(nSnow+1:nLayers) ! intent(in):    volumetric fraction of liquid water (-)
   in_groundwatr % mLayerVolFracIceTrial    = mLayerVolFracIceTrial(nSnow+1:nLayers) ! intent(in):    volumetric fraction of ice (-)
  end associate
 end subroutine initialize_in_groundwatr

 subroutine initialize_io_groundwatr(io_groundwatr,ixSaturation)
  class(io_type_groundwatr),intent(out) :: io_groundwatr ! class object for intent(inout) groundwatr arguments
  integer(i4b),intent(in)               :: ixSaturation  ! index of lowest saturated layer (NOTE: only computed on the first iteration)
  ! intent(inout) arguments
  io_groundwatr % ixSaturation = ixSaturation ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
 end subroutine initialize_io_groundwatr
 
 subroutine finalize_io_groundwatr(io_groundwatr,ixSaturation)
  class(io_type_groundwatr),intent(in)  :: io_groundwatr ! class object for intent(inout) groundwatr arguments
  integer(i4b),intent(out)              :: ixSaturation  ! index of lowest saturated layer (NOTE: only computed on the first iteration)
  ! intent(inout) arguments
  ixSaturation = io_groundwatr % ixSaturation ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
 end subroutine finalize_io_groundwatr

 subroutine finalize_out_groundwatr(out_groundwatr,dBaseflow_dWat,dBaseflow_dTk,flux_data,err,cmessage)
  class(out_type_groundwatr),intent(in) :: out_groundwatr          ! class object for intent(out) groundwatr arguments
  real(rkind),intent(out)               :: dBaseflow_dWat(:,:)     ! derivative in baseflow w.r.t. soil water characteristic
  real(rkind),intent(out)               :: dBaseflow_dTk(:,:)      ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  type(var_dlength),intent(inout)       :: flux_data               ! model fluxes for a local HRU
  integer(i4b),intent(out)              :: err                     ! error code
  character(*),intent(out)              :: cmessage                ! error message from groundwatr
  associate(&
   mLayerBaseflow => flux_data%var(iLookFLUX%mLayerBaseflow)%dat ) ! intent(out): [dp(:)]  baseflow from each soil layer (m s-1)
   ! intent(out) arguments
   mLayerBaseflow = out_groundwatr % mLayerBaseflow                ! intent(out):   baseflow from each soil layer (m s-1)
   dBaseflow_dWat = out_groundwatr % dBaseflow_dWat                ! intent(out):   derivative in baseflow w.r.t. soil water characteristic
   dBaseflow_dTk  = out_groundwatr % dBaseflow_dTk                 ! intent(out):   derivative in baseflow w.r.t. temperature (m s-1 K-1)
   err            = out_groundwatr % err                           ! intent(out):   error code
   cmessage       = out_groundwatr % cmessage                      ! intent(out):   error message
  end associate
 end subroutine finalize_out_groundwatr
 ! **** end groundwatr ****

 ! **** bigAquifer ****
 subroutine initialize_in_bigAquifer(in_bigAquifer,scalarAquiferStorageTrial,flux_data,deriv_data)
  class(in_type_bigAquifer),intent(out) :: in_bigAquifer             ! class object for intent(in) bigAquifer arguments
  real(rkind),intent(in)                :: scalarAquiferStorageTrial ! trial value of aquifer storage (m)
  type(var_dlength),intent(in)          :: flux_data                 ! model fluxes for a local HRU
  type(var_dlength),intent(in)          :: deriv_data                ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   scalarCanopyTranspiration    => flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1), &  ! intent(out): [dp]    canopy transpiration (kg m-2 s-1)
   scalarSoilDrainage           => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1),        &  ! intent(out): [dp]     drainage from the soil profile (m s-1)
   dCanopyTrans_dCanWat         => deriv_data%var(iLookDERIV%dCanopyTrans_dCanWat)%dat(1),    &  ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
   dCanopyTrans_dTCanair        => deriv_data%var(iLookDERIV%dCanopyTrans_dTCanair)%dat(1),   &  ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTCanopy        => deriv_data%var(iLookDERIV%dCanopyTrans_dTCanopy)%dat(1),   &  ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dCanopyTrans_dTGround        => deriv_data%var(iLookDERIV%dCanopyTrans_dTGround)%dat(1) )     ! intent(out): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
   ! intent(in) arguments
   in_bigAquifer % scalarAquiferStorageTrial = scalarAquiferStorageTrial ! intent(in): trial value of aquifer storage (m)
   in_bigAquifer % scalarCanopyTranspiration = scalarCanopyTranspiration ! intent(in): canopy transpiration (kg m-2 s-1)
   in_bigAquifer % scalarSoilDrainage        = scalarSoilDrainage        ! intent(in): soil drainage (m s-1)
   in_bigAquifer % dCanopyTrans_dCanWat      = dCanopyTrans_dCanWat      ! intent(in): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
   in_bigAquifer % dCanopyTrans_dTCanair     = dCanopyTrans_dTCanair     ! intent(in): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   in_bigAquifer % dCanopyTrans_dTCanopy     = dCanopyTrans_dTCanopy     ! intent(in): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
   in_bigAquifer % dCanopyTrans_dTGround     = dCanopyTrans_dTGround     ! intent(in): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
  end associate
 end subroutine initialize_in_bigAquifer
 
 subroutine initialize_io_bigAquifer(io_bigAquifer,deriv_data)
  class(io_type_bigAquifer),intent(out) :: io_bigAquifer  ! class object for intent(inout) bigAquifer arguments
  type(var_dlength),intent(in)          :: deriv_data     ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   dAquiferTrans_dTCanair       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanair)%dat(1), & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
   dAquiferTrans_dTCanopy       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanopy)%dat(1), & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
   dAquiferTrans_dTGround       => deriv_data%var(iLookDERIV%dAquiferTrans_dTGround)%dat(1), & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. ground temperature
   dAquiferTrans_dCanWat        => deriv_data%var(iLookDERIV%dAquiferTrans_dCanWat)%dat(1) )   ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy total water
   ! intent(inout) arguments
   io_bigAquifer % dAquiferTrans_dTCanair = dAquiferTrans_dTCanair       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
   io_bigAquifer % dAquiferTrans_dTCanopy = dAquiferTrans_dTCanopy       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
   io_bigAquifer % dAquiferTrans_dTGround = dAquiferTrans_dTGround       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
   io_bigAquifer % dAquiferTrans_dCanWat  = dAquiferTrans_dCanWat        ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
  end associate
 end subroutine initialize_io_bigAquifer
 
 subroutine finalize_io_bigAquifer(io_bigAquifer,deriv_data)
  class(io_type_bigAquifer),intent(in)  :: io_bigAquifer  ! class object for intent(inout) bigAquifer arguments
  type(var_dlength),intent(inout)       :: deriv_data     ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   dAquiferTrans_dTCanair       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanair)%dat(1), & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
   dAquiferTrans_dTCanopy       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanopy)%dat(1), & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
   dAquiferTrans_dTGround       => deriv_data%var(iLookDERIV%dAquiferTrans_dTGround)%dat(1), & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. ground temperature
   dAquiferTrans_dCanWat        => deriv_data%var(iLookDERIV%dAquiferTrans_dCanWat)%dat(1) )   ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy total water
   ! intent(inout) arguments
   dAquiferTrans_dTCanair = io_bigAquifer % dAquiferTrans_dTCanair       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
   dAquiferTrans_dTCanopy = io_bigAquifer % dAquiferTrans_dTCanopy       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
   dAquiferTrans_dTGround = io_bigAquifer % dAquiferTrans_dTGround       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
   dAquiferTrans_dCanWat  = io_bigAquifer % dAquiferTrans_dCanWat        ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
  end associate
 end subroutine finalize_io_bigAquifer

 subroutine finalize_out_bigAquifer(out_bigAquifer,flux_data,deriv_data,err,cmessage)
  class(out_type_bigAquifer),intent(in) :: out_bigAquifer ! class object for intent(out) bigAquifer arguments
  type(var_dlength),intent(inout)       :: flux_data      ! model fluxes for a local HRU
  type(var_dlength),intent(inout)       :: deriv_data     ! derivatives in model fluxes w.r.t. relevant state variables
  integer(i4b),intent(out)              :: err            ! error code
  character(*),intent(out)              :: cmessage       ! error message from bigAquifer
  associate(&
   scalarAquiferTranspire       => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1), & ! intent(out): [dp] transpiration loss from the aquifer (m s-1)
   scalarAquiferRecharge        => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1),  & ! intent(out): [dp] recharge to the aquifer (m s-1)
   scalarAquiferBaseflow        => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1),  & ! intent(out): [dp] total baseflow from the aquifer (m s-1)
   dBaseflow_dAquifer           => deriv_data%var(iLookDERIV%dBaseflow_dAquifer)%dat(1) )    ! intent(out): [dp(:)] derivative in baseflow flux w.r.t. aquifer storage (s-1)
   ! intent(out) arguments
   scalarAquiferTranspire = out_bigAquifer % scalarAquiferTranspire      ! intent(out):   transpiration loss from the aquifer (m s-1)
   scalarAquiferRecharge  = out_bigAquifer % scalarAquiferRecharge       ! intent(out):   recharge to the aquifer (m s-1)
   scalarAquiferBaseflow  = out_bigAquifer % scalarAquiferBaseflow       ! intent(out):   total baseflow from the aquifer (m s-1)
   dBaseflow_dAquifer     = out_bigAquifer % dBaseflow_dAquifer          ! intent(out):   change in baseflow flux w.r.t. aquifer storage (s-1)
   err                    = out_bigAquifer % err                         ! intent(out):   error code
   cmessage               = out_bigAquifer % cmessage                    ! intent(out):   error message
  end associate
 end subroutine finalize_out_bigAquifer
 ! **** end bigAquifer ****

 ! **** diagv_node ****
 subroutine initialize_in_diagv_node(in_diagv_node,iSoil,in_soilLiqFlux,diag_data,mpar_data,flux_data)
  class(in_type_diagv_node),intent(out) :: in_diagv_node                    ! class object for input diagv_node variables
  integer(i4b),intent(in)               :: iSoil                            ! index of soil layer
  type(in_type_soilLiqFlux),intent(in)  :: in_soilLiqFlux                   ! input data for soilLiqFlux
  type(var_dlength),intent(in)          :: diag_data                        ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)          :: mpar_data                        ! model parameters
  type(var_dlength),intent(in)          :: flux_data                        ! model fluxes for a local HRU

  associate(&
   ! intent(in): state variables
   mLayerMatricHeadLiqTrial => in_soilLiqFlux % mLayerMatricHeadLiqTrial, & ! liquid matric head in each layer at the current iteration (m)
   mLayerVolFracLiqTrial    => in_soilLiqFlux % mLayerVolFracLiqTrial,    & ! volumetric fraction of liquid water at the current iteration (-)
   mLayerVolFracIceTrial    => in_soilLiqFlux % mLayerVolFracIceTrial,    & ! volumetric fraction of ice at the current iteration (-)
   mLayerdTheta_dTk         => in_soilLiqFlux % mLayerdTheta_dTk,         & ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   dPsiLiq_dTemp            => in_soilLiqFlux % dPsiLiq_dTemp,            & ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! intent(in): van Genuchten and other soil parameters..
   vGn_alpha          => mpar_data%var(iLookPARAM%vGn_alpha)%dat,        & ! "alpha" parameter (m-1)
   vGn_n              => mpar_data%var(iLookPARAM%vGn_n)%dat,            & ! "n" parameter (-)
   vGn_m              => diag_data%var(iLookDIAG%scalarVGn_m)%dat,       & ! "m" parameter (-)
   mpExp              => mpar_data%var(iLookPARAM%mpExp)%dat(1),         & ! empirical exponent in macropore flow equation (-)
   theta_sat          => mpar_data%var(iLookPARAM%theta_sat)%dat,        & ! soil porosity (-)
   theta_res          => mpar_data%var(iLookPARAM%theta_res)%dat,        & ! soil residual volumetric water content (-)
   theta_mp           => mpar_data%var(iLookPARAM%theta_mp)%dat(1),      & ! volumetric liquid water content when macropore flow begins (-)
   f_impede           => mpar_data%var(iLookPARAM%f_impede)%dat(1),      & ! ice impedence factor (-)
   mLayerSatHydCond   => flux_data%var(iLookFLUX%mLayerSatHydCond)%dat,  & ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
   mLayerSatHydCondMP => flux_data%var(iLookFLUX%mLayerSatHydCondMP)%dat & ! saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
  &)
   ! input: state variables
   in_diagv_node % scalarMatricHeadLiqTrial = mLayerMatricHeadLiqTrial(iSoil) ! liquid matric head in each layer (m)
   in_diagv_node % scalarVolFracLiqTrial    = mLayerVolFracLiqTrial(iSoil)    ! volumetric fraction of liquid water in a given layer (-)
   in_diagv_node % scalarVolFracIceTrial    = mLayerVolFracIceTrial(iSoil)    ! volumetric fraction of ice in a given layer (-)
   in_diagv_node % dTheta_dTk               = mLayerdTheta_dTk(iSoil)         ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   in_diagv_node % dPsiLiq_dTemp            = dPsiLiq_dTemp(iSoil)            ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: soil parameters
   in_diagv_node % vGn_alpha          = vGn_alpha(iSoil)          ! van Genuchten "alpha" parameter (m-1)
   in_diagv_node % vGn_n              = vGn_n(iSoil)              ! van Genuchten "n" parameter (-)
   in_diagv_node % vGn_m              = vGn_m(iSoil)              ! van Genuchten "m" parameter (-)
   in_diagv_node % mpExp              = mpExp                     ! empirical exponent in macropore flow equation (-)
   in_diagv_node % theta_sat          = theta_sat(iSoil)          ! soil porosity (-)
   in_diagv_node % theta_res          = theta_res(iSoil)          ! soil residual volumetric water content (-)
   in_diagv_node % theta_mp           = theta_mp                  ! volumetric liquid water content when macropore flow begins (-)
   in_diagv_node % f_impede           = f_impede                  ! ice impedence factor (-)
   in_diagv_node % scalarSatHydCond   = mLayerSatHydCond(iSoil)   ! saturated hydraulic conductivity at the mid-point of a given layer (m s-1)
   in_diagv_node % scalarSatHydCondMP = mLayerSatHydCondMP(iSoil) ! saturated hydraulic conductivity of macropores at the mid-point of a given layer (m s-1)
  end associate
 end subroutine initialize_in_diagv_node

 subroutine finalize_out_diagv_node(out_diagv_node,iSoil,nSoil,io_soilLiqFlux,iceImpedeFac,&
                                   &dHydCond_dTemp,err,cmessage)
  class(out_type_diagv_node),intent(in)   :: out_diagv_node ! class object for output diagv_node variables
  integer(i4b),intent(in)                 :: nSoil,iSoil    ! number of soil layers and index
  type(io_type_soilLiqFlux),intent(inout) :: io_soilLiqFlux ! input-output class object for soilLiqFlux
  real(rkind),intent(inout) :: iceImpedeFac(1:nSoil)        ! ice impedence factor at layer mid-points (-)
  real(rkind),intent(inout) :: dHydCond_dTemp(1:nSoil)      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  integer(i4b),intent(out)  :: err                          ! error code
  character(*),intent(out)  :: cmessage                     ! error message

  associate(&
   ! hydraulic conductivity and derivatives
   mLayerdTheta_dPsi => io_soilLiqFlux % mLayerdTheta_dPsi,     & ! derivative in liquid water content w.r.t. matric potential (m-1)
   mLayerHydCond     => io_soilLiqFlux % mLayerHydCond,         & ! hydraulic conductivity in each soil layer (m s-1)
   dHydCond_dMatric  => io_soilLiqFlux % dHydCond_dMatric       & ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  &)
   ! output: derivative in the soil water characteristic
   mLayerdTheta_dPsi(iSoil) = out_diagv_node % scalardTheta_dPsi ! derivative in liquid water content w.r.t. matric potential (m-1)
   ! output: transmittance
   mLayerHydCond(iSoil) = out_diagv_node % scalarHydCond ! hydraulic conductivity at layer mid-points (m s-1)
   iceImpedeFac(iSoil)  = out_diagv_node % iceImpedeFac  ! ice impedence factor in each layer (-)
   dHydCond_dMatric(iSoil) = out_diagv_node % dHydCond_dMatric ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
   dHydCond_dTemp(iSoil)   = out_diagv_node % dHydCond_dTemp   ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! output: error control
   err      = out_diagv_node % err     ! error code
   cmessage = out_diagv_node % message ! error message
  end associate
 end subroutine finalize_out_diagv_node
 ! **** end diagv_node ****

 ! **** surfaceFlux ****
 subroutine initialize_in_surfaceFlux(in_surfaceFlux,nRoots,ixIce,nSoil,ibeg,iend,in_soilLiqFlux,io_soilLiqFlux,&
                                    &model_decisions,prog_data,mpar_data,flux_data,diag_data,&
                                    &iLayerHeight,dHydCond_dTemp,iceImpedeFac)
  class(in_type_surfaceFlux),intent(out) :: in_surfaceFlux ! input object for surfaceFlux
  integer(i4b),intent(in)                :: nRoots         ! number of soil layers with roots
  integer(i4b),intent(in)                :: ixIce          ! index of the lowest soil layer that contains ice
  integer(i4b),intent(in)                :: nSoil          ! number of soil layers
  integer(i4b),intent(in)                :: ibeg,iend      ! start and end indices of the soil layers in concatanated snow-soil vector
  type(in_type_soilLiqFlux),intent(in)   :: in_soilLiqFlux ! input data for soilLiqFlux
  type(io_type_soilLiqFlux),intent(in)   :: io_soilLiqFlux ! input-output class object for soilLiqFlux
  type(model_options),intent(in)         :: model_decisions(maxvarDecisions) ! the model decision structure
  type(var_dlength),intent(in)           :: prog_data      ! prognostic variables for a local HRU
  type(var_dlength),intent(in)           :: mpar_data      ! model parameters
  type(var_dlength),intent(in)           :: flux_data      ! model fluxes for a local HRU
  type(var_dlength),intent(in)           :: diag_data      ! diagnostic variables for a local HRU
  real(rkind),intent(in) :: iLayerHeight(0:nSoil)          ! height of the layer interfaces (m)
  real(rkind),intent(in) :: dHydCond_dTemp(1:nSoil)        ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  real(rkind),intent(in) :: iceImpedeFac(1:nSoil)          ! ice impedence factor at layer mid-points (-)

  associate(&
   ! model control
   firstSplitOper         => in_soilLiqFlux % firstSplitOper,                     & ! flag to compute infiltration
   ixBcLowerSoilHydrology => model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision,& ! index of the lower boundary conditions for soil hydrology
   ixBcUpperSoilHydrology => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,& ! index defining the type of boundary conditions
   ixInfRateMax           => model_decisions(iLookDECISIONS%infRateMax)%iDecision,& ! index of the maximum infiltration rate parameterization
   surfRun_SE             => model_decisions(iLookDECISIONS%surfRun_SE)%iDecision,& ! index defining the saturation excess surface runoff method
   ix_groundwatr          => model_decisions(iLookDECISIONS%groundwatr)%iDecision & ! index defining the groundwater parameterization
  &)
   ! intent(in): model control
   in_surfaceFlux % firstSplitOper = firstSplitOper          ! flag indicating if desire to compute infiltration
   in_surfaceFlux % bc_lower       = ixBcLowerSoilHydrology  ! index defining the type of boundary conditions at the bottom of the soil
   in_surfaceFlux % bc_upper       = ixBcUpperSoilHydrology  ! index defining the type of boundary conditions (Neumann or Dirichlet)
   in_surfaceFlux % ixInfRateMax   = ixInfRateMax            ! index defining the maximum infiltration rate parameterization (GreenAmpt or topmodel_GA)
   in_surfaceFlux % surfRun_SE     = surfRun_SE              ! index defining the saturation excess surface runoff method
   in_surfaceFlux % ix_groundwatr  = ix_groundwatr           ! index defining the groundwater parameterization
   in_surfaceFlux % nRoots         = nRoots                  ! number of layers that contain roots
   in_surfaceFlux % ixIce          = ixIce                   ! index of lowest ice layer
   in_surfaceFlux % nSoil          = nSoil                   ! number of soil layers
  end associate

  associate(&
   ! state variables
   mLayerTempTrial          => in_soilLiqFlux % mLayerTempTrial,          & ! intent(in): temperature in each layer at the current iteration (m)
   mLayerMatricHeadLiqTrial => in_soilLiqFlux % mLayerMatricHeadLiqTrial, & ! liquid matric head in each layer at the current iteration (m)
   mLayerMatricHeadTrial    => in_soilLiqFlux % mLayerMatricHeadTrial,    & ! intent(in): matric head in each layer at the current iteration (m)
   mLayerVolFracLiqTrial    => in_soilLiqFlux % mLayerVolFracLiqTrial,    & ! volumetric fraction of liquid water at the current iteration (-)
   mLayerVolFracIceTrial    => in_soilLiqFlux % mLayerVolFracIceTrial,    & ! volumetric fraction of ice at the current iteration (-)
   scalarTotalSoilLiq       => diag_data%var(iLookDIAG%scalarTotalSoilLiq)%dat(1) & ! total liquid water in the soil column (kg m-2)
  &)
   ! intent(in): state variables
   in_surfaceFlux % mLayerTemp          = mLayerTempTrial             ! temperature (K)
   in_surfaceFlux % scalarMatricHeadLiq = mLayerMatricHeadLiqTrial(1) ! liquid matric head in the upper-most soil layer (m)
   in_surfaceFlux % mLayerMatricHead    = mLayerMatricHeadTrial       ! matric head in each soil layer (m)
   in_surfaceFlux % scalarVolFracLiq    = mLayerVolFracLiqTrial(1)    ! volumetric liquid water content the upper-most soil layer (-)
   in_surfaceFlux % scalarTotalSoilLiq  = scalarTotalSoilLiq          ! total liquid water in the soil column (kg m-2)
   in_surfaceFlux % mLayerVolFracLiq    = mLayerVolFracLiqTrial       ! volumetric liquid water content in each soil layer (-)
   in_surfaceFlux % mLayerVolFracIce    = mLayerVolFracIceTrial       ! volumetric ice content in each soil layer (-)
  end associate

  associate(&
   ! pre-computed derivatives
   mLayerdTheta_dTk       => in_soilLiqFlux % mLayerdTheta_dTk,      & ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   mLayerdTheta_dPsi      => io_soilLiqFlux % mLayerdTheta_dPsi      & ! derivative in liquid water content w.r.t. matric potential (m-1)
  &)
   ! intent(in): pre-computed derivatives
   in_surfaceFlux % dTheta_dTk             = mLayerdTheta_dTk       ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   in_surfaceFlux % dTheta_dPsi            = mLayerdTheta_dPsi      ! derivative in liquid water content w.r.t. matric potential (m-1)
  end associate

  associate(&
   ! depth of each soil layer (m)
   mLayerDepth         => prog_data%var(iLookPROG%mLayerDepth)%dat(ibeg:iend) & ! depth of the layer (m)
  &)
   ! intent(in): depth of each soil layer (m)
   in_surfaceFlux % mLayerDepth     = mLayerDepth  ! depth of each soil layer (m)
   in_surfaceFlux % iLayerHeight    = iLayerHeight ! height at the interface of each layer (m)
  end associate

  associate(&
   ! boundary conditions
   upperBoundHead      => mpar_data%var(iLookPARAM%upperBoundHead)%dat(1) & ! upper boundary condition for matric head (m)
  &)
   ! intent(in): boundary conditions
   in_surfaceFlux % upperBoundHead  = upperBoundHead  ! upper boundary condition (m)
  end associate

  associate(&
   ! flux at the upper boundary
   scalarRainPlusMelt  => in_soilLiqFlux % scalarRainPlusMelt & ! rain plus melt (m s-1)
  &)
   ! intent(in): flux at the upper boundary
   in_surfaceFlux % scalarRainPlusMelt = scalarRainPlusMelt ! rain plus melt (m s-1)
  end associate

  associate(&
   ! transmittance
   iLayerSatHydCond    => flux_data%var(iLookFLUX%iLayerSatHydCond)%dat & ! saturated hydraulic conductivity at the interface of each layer (m s-1)
  &)
   ! intent(in): transmittance
   in_surfaceFlux % surfaceSatHydCond = iLayerSatHydCond(0) ! saturated hydraulic conductivity at the surface (m s-1)
   in_surfaceFlux % dHydCond_dTemp    = dHydCond_dTemp(1)   ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   in_surfaceFlux % iceImpedeFac      = iceImpedeFac(1)     ! ice impedence factor in the upper-most soil layer (-)
  end associate

  associate(&
   ! soil parameters
   vGn_alpha           => mpar_data%var(iLookPARAM%vGn_alpha)%dat,             & ! "alpha" parameter (m-1)
   vGn_n               => mpar_data%var(iLookPARAM%vGn_n)%dat,                 & ! "n" parameter (-)
   vGn_m               => diag_data%var(iLookDIAG%scalarVGn_m)%dat,            & ! "m" parameter (-)
   theta_sat           => mpar_data%var(iLookPARAM%theta_sat)%dat,             & ! soil porosity (-)
   theta_res           => mpar_data%var(iLookPARAM%theta_res)%dat,             & ! soil residual volumetric water content (-)
   qSurfScale          => mpar_data%var(iLookPARAM%qSurfScale)%dat(1),         & ! scaling factor in the surface runoff parameterization (-)
   zScale_TOPMODEL     => mpar_data%var(iLookPARAM%zScale_TOPMODEL)%dat(1),    & ! TOPMODEL scaling factor (m)
   rootingDepth        => mpar_data%var(iLookPARAM%rootingDepth)%dat(1),       & ! rooting depth (m)
   wettingFrontSuction => mpar_data%var(iLookPARAM%wettingFrontSuction)%dat(1),& ! Green-Ampt wetting front suction (m)
   soilIceScale        => mpar_data%var(iLookPARAM%soilIceScale)%dat(1),       & ! scaling factor for depth of soil ice, used to get frozen fraction (m)
   soilIceCV           => mpar_data%var(iLookPARAM%soilIceCV)%dat(1)           & ! CV of depth of soil ice, used to get frozen fraction (-)
  &)
   ! intent(in): soil parameters
   in_surfaceFlux % vGn_alpha           = vGn_alpha(1)        ! van Genuchten "alpha" parameter (m-1)
   in_surfaceFlux % vGn_n               = vGn_n(1)            ! van Genuchten "n" parameter (-)
   in_surfaceFlux % vGn_m               = vGn_m(1)            ! van Genuchten "m" parameter (-)
   in_surfaceFlux % theta_sat           = theta_sat(1)        ! soil porosity (-)
   in_surfaceFlux % theta_res           = theta_res(1)        ! soil residual volumetric water content (-)
   in_surfaceFlux % qSurfScale          = qSurfScale          ! scaling factor in the surface runoff parameterization (-)
   in_surfaceFlux % zScale_TOPMODEL     = zScale_TOPMODEL     ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
   in_surfaceFlux % rootingDepth        = rootingDepth        ! rooting depth (m)
   in_surfaceFlux % wettingFrontSuction = wettingFrontSuction ! Green-Ampt wetting front suction (m)
   in_surfaceFlux % soilIceScale        = soilIceScale        ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   in_surfaceFlux % soilIceCV           = soilIceCV           ! soil ice CV in Gamma distribution used to define frozen area (-)
  end associate

  ! intent(in): FUSE parameters
  associate(&
   FUSE_Ac_max   => mpar_data%var(iLookPARAM%FUSE_Ac_max  )%dat(1), & ! FUSE PRMS max saturated area
   FUSE_phi_tens => mpar_data%var(iLookPARAM%FUSE_phi_tens)%dat(1), & ! FUSE PRMS tension fraction
   FUSE_b        => mpar_data%var(iLookPARAM%FUSE_b       )%dat(1), & ! FUSE ARNO/VIC exponent
   FUSE_lambda   => mpar_data%var(iLookPARAM%FUSE_lambda  )%dat(1), & ! FUSE TOPMODEL gamma distribution lambda parameter
   FUSE_chi      => mpar_data%var(iLookPARAM%FUSE_chi     )%dat(1), & ! FUSE TOPMODEL chi   distribution lambda parameter
   FUSE_mu       => mpar_data%var(iLookPARAM%FUSE_mu      )%dat(1), & ! FUSE TOPMODEL mu    distribution lambda parameter
   FUSE_n        => mpar_data%var(iLookPARAM%FUSE_n       )%dat(1)  & ! FUSE TOPMODEL exponent
  &)
   in_surfaceFlux % FUSE_Ac_max   = FUSE_Ac_max   ! FUSE PRMS max saturated area
   in_surfaceFlux % FUSE_phi_tens = FUSE_phi_tens ! FUSE PRMS tension fraction
   in_surfaceFlux % FUSE_b        = FUSE_b        ! FUSE ARNO/VIC exponent
   in_surfaceFlux % FUSE_lambda   = FUSE_lambda   ! FUSE TOPMODEL gamma distribution lambda parameter
   in_surfaceFlux % FUSE_chi      = FUSE_chi      ! FUSE TOPMODEL chi   distribution lambda parameter
   in_surfaceFlux % FUSE_mu       = FUSE_mu       ! FUSE TOPMODEL mu    distribution lambda parameter
   in_surfaceFlux % FUSE_n        = FUSE_n        ! FUSE TOPMODEL exponent
  end associate
 end subroutine initialize_in_surfaceFlux

 subroutine initialize_io_surfaceFlux(io_surfaceFlux,nSoil,io_soilLiqFlux,iLayerHydCond)
  class(io_type_surfaceFlux),intent(out) :: io_surfaceFlux ! input-output object for surfaceFlux
  integer(i4b),intent(in)                :: nSoil          ! number of soil layers
  type(io_type_soilLiqFlux),intent(in)   :: io_soilLiqFlux ! input-output class object for soilLiqFlux
  real(rkind),intent(in) :: iLayerHydCond(0:nSoil)         ! hydraulic conductivity at layer interface (m s-1)

  associate(&
   ! fluxes at layer interfaces and surface runoff
   xMaxInfilRate       => io_soilLiqFlux % scalarMaxInfilRate,      & ! maximum infiltration rate (m s-1)
   scalarInfilArea     => io_soilLiqFlux % scalarInfilArea,         & ! fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea => io_soilLiqFlux % scalarSaturatedArea,     & ! fraction of area that is considered saturated (-)
   scalarFrozenArea    => io_soilLiqFlux % scalarFrozenArea,        & ! fraction of area that is considered impermeable due to soil ice (-)
   scalarSoilControl   => io_soilLiqFlux % scalarSoilControl,       & ! soil control on infiltration for derivative
   scalarSurfaceInfiltration => io_soilLiqFlux % scalarInfiltration & ! surface infiltration (m s-1)
  &)
   ! intent(inout): hydraulic conductivity at the surface
   io_surfaceFlux % surfaceHydCond = iLayerHydCond(0)         ! hydraulic conductivity at the surface (m s-1)
   ! intent(inout): fluxes at layer interfaces and surface runoff
   io_surfaceFlux % xMaxInfilRate       = xMaxInfilRate       ! maximum infiltration rate (m s-1)
   io_surfaceFlux % scalarInfilArea     = scalarInfilArea     ! fraction of area where water can infiltrate, may be frozen (-)
   io_surfaceFlux % scalarSaturatedArea = scalarSaturatedArea ! fraction of area that is considered saturated (-)
   io_surfaceFlux % scalarFrozenArea    = scalarFrozenArea    ! fraction of area that is considered impermeable due to soil ice (-)
   io_surfaceFlux % scalarSoilControl   = scalarSoilControl   ! soil control on infiltration for derivative
   io_surfaceFlux % scalarSurfaceInfiltration = scalarSurfaceInfiltration ! surface infiltration (m s-1)
  end associate
 end subroutine initialize_io_surfaceFlux

 subroutine finalize_io_surfaceFlux(io_surfaceFlux,nSoil,io_soilLiqFlux,iLayerHydCond)
  class(io_type_surfaceFlux),intent(in)   :: io_surfaceFlux ! input-output object for surfaceFlux
  integer(i4b),intent(in)                 :: nSoil          ! number of soil layers
  type(io_type_soilLiqFlux),intent(inout) :: io_soilLiqFlux ! input-output class object for soilLiqFlux
  real(rkind),intent(inout) :: iLayerHydCond(0:nSoil)       ! hydraulic conductivity at layer interface (m s-1)

  associate(&
   ! fluxes at layer interfaces and surface runoff
   xMaxInfilRate       => io_soilLiqFlux % scalarMaxInfilRate,  & ! maximum infiltration rate (m s-1)
   scalarInfilArea     => io_soilLiqFlux % scalarInfilArea,     & ! fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea => io_soilLiqFlux % scalarSaturatedArea, & ! fraction of area that is considered saturated (-)
   scalarFrozenArea    => io_soilLiqFlux % scalarFrozenArea,    & ! fraction of area that is considered impermeable due to soil ice (-)
   scalarSoilControl   => io_soilLiqFlux % scalarSoilControl,   & ! soil control on infiltration for derivative
   scalarSurfaceInfiltration => io_soilLiqFlux % scalarInfiltration & ! surface infiltration rate (m s-1)
  &)
   ! intent(inout): hydraulic conductivity at the surface
   iLayerHydCond(0) = io_surfaceFlux % surfaceHydCond         ! hydraulic conductivity at the surface (m s-1) 
   ! intent(inout): fluxes at layer interfaces and surface runoff
   xMaxInfilRate       = io_surfaceFlux % xMaxInfilRate       ! maximum infiltration rate (m s-1)                                   
   scalarInfilArea     = io_surfaceFlux % scalarInfilArea     ! fraction of area where water can infiltrate, may be frozen (-)
   scalarSaturatedArea = io_surfaceFlux % scalarSaturatedArea ! fraction of area that is considered saturated (-)
   scalarFrozenArea    = io_surfaceFlux % scalarFrozenArea    ! fraction of area that is considered impermeable due to soil ice (-)
   scalarSoilControl = io_surfaceFlux % scalarSoilControl     ! soil control on infiltration for derivative
   scalarSurfaceInfiltration = io_surfaceFlux % scalarSurfaceInfiltration ! surface infiltration (m s-1)
  end associate
 end subroutine finalize_io_surfaceFlux

 subroutine finalize_out_surfaceFlux(out_surfaceFlux,io_soilLiqFlux,err,message)
  class(out_type_surfaceFlux),intent(in)  :: out_surfaceFlux ! output object for surfaceFlux
  type(io_type_soilLiqFlux),intent(inout) :: io_soilLiqFlux  ! input-output class object for soilLiqFlux
  integer(i4b),intent(out)  :: err       ! error code
  character(*),intent(out)  :: message   ! error message

  associate(&
   ! intent(out): surface runoff and infiltration
   scalarSurfaceRunoff       => io_soilLiqFlux % scalarSurfaceRunoff,    & ! surface runoff (m s-1)
   scalarSurfaceRunoff_IE    => io_soilLiqFlux % scalarSurfaceRunoff_IE, & ! infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE    => io_soilLiqFlux % scalarSurfaceRunoff_SE, & ! saturation excess surface runoff (m s-1)
    ! intent(inout): derivatives in surface infiltration in the upper-most soil layer w.r.t ... 
   dq_dHydStateLayerSurfVec => io_soilLiqFlux % dq_dHydStateLayerSurfVec, & ! ... hydrology state above soil snow or canopy and every soil layer (m s-1 or s-1)
   dq_dNrgStateLayerSurfVec => io_soilLiqFlux % dq_dNrgStateLayerSurfVec  & ! ... temperature above soil snow or canopy and every soil layer (m s-1 or s-1)
  &)
   ! intent(out): surface runoff
   scalarSurfaceRunoff       = out_surfaceFlux % scalarSurfaceRunoff    ! surface runoff (m s-1)
   scalarSurfaceRunoff_IE    = out_surfaceFlux % scalarSurfaceRunoff_IE ! infiltration excess surface runoff (m s-1)
   scalarSurfaceRunoff_SE    = out_surfaceFlux % scalarSurfaceRunoff_SE ! saturation excess surface runoff (m s-1)
   ! intent(inout): derivatives in surface infiltration in the upper-most soil layer w.r.t. ...
   dq_dHydStateLayerSurfVec  = out_surfaceFlux % dq_dHydStateVec ! ... hydrology state in every soil layer (m s-1 or s-1)
   dq_dNrgStateLayerSurfVec  = out_surfaceFlux % dq_dNrgStateVec ! ... energy state in every soil layer (m s-1 K-1)
  end associate
  ! intent(out): error control
  err      = out_surfaceFlux % err     ! error code
  message  = out_surfaceFlux % message ! error message
 end subroutine finalize_out_surfaceFlux
 ! **** end surfaceFlux ****

 ! **** iLayerFlux ****
 subroutine initialize_in_iLayerFlux(in_iLayerFlux,iLayer,nSoil,ibeg,iend,in_soilLiqFlux,io_soilLiqFlux,&
                                    &prog_data,dHydCond_dTemp)
  class(in_type_iLayerFlux),intent(out) :: in_iLayerFlux  ! class object for input iLayerFlux variables
  integer(i4b),intent(in)               :: nSoil,iLayer   ! number of soil layers and index
  integer(i4b),intent(in)               :: ibeg,iend     ! start and end indices of the soil layers in concatanated snow-soil vector
  type(in_type_soilLiqFlux),intent(in)  :: in_soilLiqFlux ! input class object for soilLiqFlux
  type(io_type_soilLiqFlux),intent(in)  :: io_soilLiqFlux ! input-output class object for soilLiqFlux
  type(var_dlength),intent(in)          :: prog_data      ! prognostic variables for a local HRU
  real(rkind),intent(in) :: dHydCond_dTemp(1:nSoil)   ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)

  associate(&
   ! intent(in): state variables (adjacent layers)
   mLayerMatricHeadLiqTrial => in_soilLiqFlux % mLayerMatricHeadLiqTrial, & ! liquid matric head in each layer at the current iteration (m)
   mLayerVolFracLiqTrial    => in_soilLiqFlux % mLayerVolFracLiqTrial,    & ! volumetric fraction of liquid water at the current iteration (-)
   ! intent(in): model coordinate variables (adjacent layers)
   mLayerHeight => prog_data%var(iLookPROG%mLayerHeight)%dat(ibeg:iend), & ! height of the layer mid-point (m)
   ! intent(in): temperature derivatives
   dPsiLiq_dTemp => in_soilLiqFlux % dPsiLiq_dTemp, & ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! intent(in): transmittance (adjacent layers)
   mLayerHydCond => io_soilLiqFlux % mLayerHydCond, & ! hydraulic conductivity in each soil layer (m s-1)
   ! intent(in): transmittance derivatives (adjacent layers)
   dHydCond_dMatric => io_soilLiqFlux % dHydCond_dMatric & ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  &)
   ! intent(in): state variables (adjacent layers)
   in_iLayerFlux % nodeMatricHeadLiqTrial = mLayerMatricHeadLiqTrial(iLayer:iLayer+1) ! liquid matric head at the soil nodes (m)
   ! intent(in): model coordinate variables (adjacent layers)
   in_iLayerFlux % nodeHeight = mLayerHeight(iLayer:iLayer+1) ! height of the soil nodes (m)
   ! intent(in): temperature derivatives
   in_iLayerFlux % dPsiLiq_dTemp  = dPsiLiq_dTemp(iLayer:iLayer+1)  ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   in_iLayerFlux % dHydCond_dTemp = dHydCond_dTemp(iLayer:iLayer+1) ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! intent(in): transmittance (adjacent layers)
   in_iLayerFlux % nodeHydCondTrial = mLayerHydCond(iLayer:iLayer+1) ! hydraulic conductivity at the soil nodes (m s-1)
   ! intent(in): transmittance derivatives (adjacent layers) for hydraulic ...
   in_iLayerFlux % dHydCond_dMatric = dHydCond_dMatric(iLayer:iLayer+1) ! ... conductivity w.r.t. change in matric head (s-1)
  end associate
 end subroutine initialize_in_iLayerFlux

 subroutine finalize_out_iLayerFlux(out_iLayerFlux,iLayer,nSoil,io_soilLiqFlux,iLayerHydCond,err,cmessage)
  class(out_type_iLayerFlux),intent(in)  :: out_iLayerFlux ! class object for output iLayerFlux variables
  integer(i4b),intent(in)                :: nSoil,iLayer   ! number of soil layers and index
  type(io_type_soilLiqFlux),intent(inout) :: io_soilLiqFlux  ! input-output class object for soilLiqFlux
  real(rkind),intent(inout) :: iLayerHydCond(0:nSoil) ! hydraulic conductivity at layer interface (m s-1)
  integer(i4b),intent(out)  :: err       ! error code
  character(*),intent(out)  :: cmessage  ! error message

  associate(&
   ! intent(out): vertical flux at the layer interface (scalars)
   iLayerLiqFluxSoil => io_soilLiqFlux % iLayerLiqFluxSoil,& ! liquid flux at soil layer interfaces (m s-1)
   ! intent(out): derivatives in fluxes in the layer above and layer below w.r.t ...
   dq_dHydStateAbove => io_soilLiqFlux % dq_dHydStateAbove,& ! ... state variables in the layer above
   dq_dHydStateBelow => io_soilLiqFlux % dq_dHydStateBelow,& ! ... state variables in the layer below
   dq_dNrgStateAbove => io_soilLiqFlux % dq_dNrgStateAbove,& ! ... temperature in the layer above (m s-1 K-1)
   dq_dNrgStateBelow => io_soilLiqFlux % dq_dNrgStateBelow & ! ... temperature in the layer below (m s-1 K-1)
  &)
   ! intent(out): tranmsmittance at the layer interface (scalars)
   iLayerHydCond(iLayer) = out_iLayerFlux % iLayerHydCond         ! hydraulic conductivity at the interface between layers (m s-1)
   ! intent(out): vertical flux at the layer interface (scalars)
   iLayerLiqFluxSoil(iLayer) = out_iLayerFlux % iLayerLiqFluxSoil ! vertical flux of liquid water at the layer interface (m s-1)
   ! intent(out): derivatives in fluxes in the layer above and layer below w.r.t. ... 
   dq_dHydStateAbove(iLayer) = out_iLayerFlux % dq_dHydStateAbove ! ... matric head or volumetric lquid water in the layer above (m s-1 or s-1)
   dq_dHydStateBelow(iLayer) = out_iLayerFlux % dq_dHydStateBelow ! ... matric head or volumetric lquid water in the layer below (m s-1 or s-1)
   dq_dNrgStateAbove(iLayer) = out_iLayerFlux % dq_dNrgStateAbove ! ... temperature in the layer above (m s-1 K-1)
   dq_dNrgStateBelow(iLayer) = out_iLayerFlux % dq_dNrgStateBelow ! ... temperature in the layer below (m s-1 K-1)
   ! intent(out): error control
   err      = out_iLayerFlux % err     ! error code
   cmessage = out_iLayerFlux % message ! error message
  end associate
 end subroutine finalize_out_iLayerFlux
 ! **** end iLayerFlux ****

 ! **** qDrainFlux ****
 subroutine initialize_in_qDrainFlux(in_qDrainFlux,nSoil,ibeg,iend,in_soilLiqFlux,io_soilLiqFlux,model_decisions,&
                                    &prog_data,mpar_data,flux_data,diag_data,iceImpedeFac,&
                                    &dHydCond_dTemp)
  class(in_type_qDrainFlux),intent(out) :: in_qDrainFlux  ! class object for input qDrainFlux variables
  integer(i4b),intent(in)               :: nSoil          ! number of soil layers
  integer(i4b),intent(in)               :: ibeg,iend      ! start and end indices of the soil layers in concatanated snow-soil vector
  type(in_type_soilLiqFlux),intent(in)  :: in_soilLiqFlux ! input class object for soilLiqFlux
  type(io_type_soilLiqFlux),intent(in)  :: io_soilLiqFlux ! input-output class object for soilLiqFlux
  type(model_options),intent(in)        :: model_decisions(maxvarDecisions) ! the model decision structure
  type(var_dlength),intent(in)          :: prog_data      ! prognostic variables for a local HRU
  type(var_dlength),intent(in)          :: mpar_data      ! model parameters
  type(var_dlength),intent(in)          :: flux_data      ! model fluxes for a local HRU
  type(var_dlength),intent(in)          :: diag_data      ! diagnostic variables for a local HRU
  real(rkind),intent(in) :: iceImpedeFac(1:nSoil)     ! ice impedence factor at layer mid-points (-)
  real(rkind),intent(in) :: dHydCond_dTemp(1:nSoil)   ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)

  associate(&
   ! intent(in): model control
   ixBcLowerSoilHydrology => model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision,& ! index of the lower boundary conditions for soil hydrology
   ! intent(in): state variables
   mLayerMatricHeadLiqTrial => in_soilLiqFlux % mLayerMatricHeadLiqTrial, & ! liquid matric head in each layer at the current iteration (m)
   mLayerVolFracLiqTrial    => in_soilLiqFlux % mLayerVolFracLiqTrial,    & ! volumetric fraction of liquid water at the current iteration (-)
   ! intent(in): model coordinate variables
   mLayerDepth  => prog_data%var(iLookPROG%mLayerDepth)%dat(ibeg:iend), & ! depth of the layer (m)
   mLayerHeight => prog_data%var(iLookPROG%mLayerHeight)%dat(ibeg:iend),& ! height of the layer mid-point (m)
   ! intent(in): boundary conditions
   lowerBoundHead  => mpar_data%var(iLookPARAM%lowerBoundHead)%dat(1), & ! lower boundary condition for matric head (m)
   ! intent(in): derivative in the soil water characteristic
   dPsiLiq_dTemp     => in_soilLiqFlux % dPsiLiq_dTemp,         & ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! intent(in): transmittance
   iLayerSatHydCond => flux_data%var(iLookFLUX%iLayerSatHydCond)%dat,& ! saturated hydraulic conductivity at the interface of each layer (m s-1)
   mLayerHydCond    => io_soilLiqFlux % mLayerHydCond,                & ! hydraulic conductivity in each soil layer (m s-1)
   ! intent(in): transmittance derivatives
   dHydCond_dMatric => io_soilLiqFlux % dHydCond_dMatric,& ! derivative in hydraulic conductivity w.r.t matric head (s-1)
   ! intent(in): soil parameters
   vGn_alpha       => mpar_data%var(iLookPARAM%vGn_alpha)%dat,         & ! "alpha" parameter (m-1)
   vGn_n           => mpar_data%var(iLookPARAM%vGn_n)%dat,             & ! "n" parameter (-)
   vGn_m           => diag_data%var(iLookDIAG%scalarVGn_m)%dat,        & ! "m" parameter (-)
   theta_sat       => mpar_data%var(iLookPARAM%theta_sat)%dat,         & ! soil porosity (-)
   theta_res       => mpar_data%var(iLookPARAM%theta_res)%dat,         & ! soil residual volumetric water content (-)
   kAnisotropic    => mpar_data%var(iLookPARAM%kAnisotropic)%dat(1),   & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
   zScale_TOPMODEL => mpar_data%var(iLookPARAM%zScale_TOPMODEL)%dat(1) & ! TOPMODEL scaling factor (m)
  &)
   ! intent(in): model control
   in_qDrainFlux % bc_lower      = ixBcLowerSoilHydrology ! index defining the type of boundary conditions
   ! intent(in): state variables
   in_qDrainFlux % nodeMatricHeadLiq = mLayerMatricHeadLiqTrial(nSoil) ! liquid matric head in the lowest unsaturated node (m)
   ! intent(in): model coordinate variables
   in_qDrainFlux % nodeDepth  = mLayerDepth(nSoil)  ! depth of the lowest unsaturated soil layer (m)
   in_qDrainFlux % nodeHeight = mLayerHeight(nSoil) ! height of the lowest unsaturated soil node (m)
   ! intent(in): boundary conditions
   in_qDrainFlux % lowerBoundHead  = lowerBoundHead  ! lower boundary condition (m)
   ! intent(in): derivative in the soil water characteristic
   in_qDrainFlux % node_dPsiLiq_dTemp = dPsiLiq_dTemp(nSoil)     ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! intent(in): transmittance
   in_qDrainFlux % surfaceSatHydCond = iLayerSatHydCond(0)     ! saturated hydraulic conductivity at the surface (m s-1)
   in_qDrainFlux % bottomSatHydCond  = iLayerSatHydCond(nSoil) ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   in_qDrainFlux % nodeHydCond       = mLayerHydCond(nSoil)    ! hydraulic conductivity at the node itself (m s-1)
   in_qDrainFlux % iceImpedeFac      = iceImpedeFac(nSoil)     ! ice impedence factor in the lower-most soil layer (-)
   in_qDrainFlux % dHydCond_dMatric = dHydCond_dMatric(nSoil) ! ... matric head (s-1)
   in_qDrainFlux % dHydCond_dTemp   = dHydCond_dTemp(nSoil)   ! ... temperature (m s-1 K-1)
   ! intent(in): soil parameters
   in_qDrainFlux % vGn_alpha       = vGn_alpha(nSoil) ! van Genuchten "alpha" parameter (m-1)
   in_qDrainFlux % vGn_n           = vGn_n(nSoil)     ! van Genuchten "n" parameter (-)
   in_qDrainFlux % vGn_m           = vGn_m(nSoil)     ! van Genuchten "m" parameter (-)
   in_qDrainFlux % theta_sat       = theta_sat(nSoil) ! soil porosity (-)
   in_qDrainFlux % theta_res       = theta_res(nSoil) ! soil residual volumetric water content (-)
   in_qDrainFlux % kAnisotropic    = kAnisotropic     ! anisotropy factor for lateral hydraulic conductivity (-)
   in_qDrainFlux % zScale_TOPMODEL = zScale_TOPMODEL  ! TOPMODEL scaling factor (m)
  end associate
 end subroutine initialize_in_qDrainFlux

 subroutine finalize_out_qDrainFlux(out_qDrainFlux,nSoil,io_soilLiqFlux,iLayerHydCond,err,cmessage)
  class(out_type_qDrainFlux),intent(in) :: out_qDrainFlux    ! class object for output qDrainFlux variables
  integer(i4b),intent(in)   :: nSoil                         ! number of soil layers
  type(io_type_soilLiqFlux),intent(inout) :: io_soilLiqFlux  ! input-output class object for soilLiqFlux
  real(rkind),intent(inout) :: iLayerHydCond(0:nSoil)        ! hydraulic conductivity at layer interface (m s-1)
  integer(i4b),intent(out)  :: err                           ! error code
  character(*),intent(out)  :: cmessage                      ! error message

  associate(&
   ! intent(out): drainage flux
   iLayerLiqFluxSoil => io_soilLiqFlux % iLayerLiqFluxSoil,& ! liquid flux at soil layer interfaces (m s-1)
   ! intent(out): derivatives in drainage flux w.r.t. ...
   dq_dHydStateAbove => io_soilLiqFlux % dq_dHydStateAbove,& ! ... state variables in the layer above
   dq_dNrgStateAbove => io_soilLiqFlux % dq_dNrgStateAbove & ! ... temperature in the layer above (m s-1 K-1)
  &)
   ! intent(out): hydraulic conductivity at the bottom
   iLayerHydCond(nSoil) = out_qDrainFlux % bottomHydCond ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   ! intent(out): drainage flux
   iLayerLiqFluxSoil(nSoil) = out_qDrainFlux % scalarDrainage    ! drainage flux (m s-1)
   ! intent(out): derivatives in drainage flux w.r.t. ...
   dq_dHydStateAbove(nSoil) = out_qDrainFlux % dq_dHydStateUnsat ! ... change in hydrology state in lowest unsaturated node (m s-1 or s-1)
   dq_dNrgStateAbove(nSoil) = out_qDrainFlux % dq_dNrgStateUnsat ! ... change in energy state in lowest unsaturated node (m s-1 or s-1)
  end associate
  ! intent(out): error control
  err      = out_qDrainFlux % err     ! error code 
  cmessage = out_qDrainFlux % message ! error message
 end subroutine finalize_out_qDrainFlux
 ! **** end qDrainFlux ****

 ! **** stateFilter ****
 subroutine finalize_out_stateFilter(out_stateFilter,err,cmessage)
  class(out_type_stateFilter),intent(in) :: out_stateFilter   ! class object for intent(out) stateFilter arguments
  integer(i4b),intent(out)               :: err               ! intent(out): error code
  character(*),intent(out)               :: cmessage          ! intent(out): error message
  err      = out_stateFilter % err                            ! intent(out): error code
  cmessage = out_stateFilter % cmessage                       ! intent(out): error message
 end subroutine finalize_out_stateFilter
 ! **** end stateFilter ****

 ! **** indexSplit ****
 subroutine initialize_in_indexSplit(in_indexSplit,nSnow,nSoil,nLayers,nSubset)
  class(in_type_indexSplit),intent(out) :: in_indexSplit    ! class object for intent(in) indexSplit arguments
  integer(i4b),intent(in)               :: nSnow            ! intent(in): number of snow layers
  integer(i4b),intent(in)               :: nSoil            ! intent(in): number of soil layers
  integer(i4b),intent(in)               :: nLayers          ! intent(in): total number of layers
  integer(i4b),intent(in)               :: nSubset          ! intent(in): number of states in the subset
  in_indexSplit % nSnow   = nSnow                           ! intent(in): number of snow layers          
  in_indexSplit % nSoil   = nSoil                           ! intent(in): number of soil layers
  in_indexSplit % nLayers = nLayers                         ! intent(in): total number of layers
  in_indexSplit % nSubset = nSubset                         ! intent(in): number of states in the subset
 end subroutine initialize_in_indexSplit

 subroutine finalize_out_indexSplit(out_indexSplit,err,cmessage)
  class(out_type_indexSplit),intent(in) :: out_indexSplit   ! class object for intent(out) indexSplit arguments
  integer(i4b),intent(out)              :: err              ! intent(out): error code
  character(*),intent(out)              :: cmessage         ! intent(out): error message
  err      = out_indexSplit % err                           ! intent(out): error code    
  cmessage = out_indexSplit % cmessage                      ! intent(out): error message
 end subroutine finalize_out_indexSplit
 ! **** end indexSplit ****

 ! **** varSubstep ****
 subroutine initialize_in_varSubstep(in_varSubstep,dt,dtInit,dt_min,whole_step,nSubset,&
                                     doAdjustTemp,firstSubStep,computeVegFlux,ixSolution,scalar,iStateSplit,fluxMask)
  class(in_type_varSubstep),intent(out) :: in_varSubstep  ! class object for intent(in) varSubstep arguments
  real(rkind),intent(in)                :: dt             ! time step (s)
  real(rkind),intent(in)                :: dtInit         ! initial time step (s)
  real(rkind),intent(in)                :: dt_min         ! minimum time step (s) 
  real(rkind),intent(in)                :: whole_step     ! length of whole step for surface drainage and average flux
  integer(i4b),intent(in)               :: nSubset        ! total number of variables in the state subset
  logical(lgt),intent(in)               :: doAdjustTemp   ! flag to indicate if we adjust the temperature
  logical(lgt),intent(in)               :: firstSubStep   ! flag to denote first sub-step
  logical(lgt),intent(in)               :: computeVegFlux ! flag to denote if computing energy flux over vegetation
  integer(i4b),intent(in)               :: ixSolution     ! index of solution method
  integer(i4b),intent(in)               :: scalar         ! scalar solution method
  integer(i4b),intent(in)               :: iStateSplit    ! index of the layer in the splitting operation
  type(var_flagVec),intent(in)          :: fluxMask       ! mask for the fluxes used in this given state subset
 
  ! intent(in) arguments
  in_varSubstep % dt             = dt                     ! intent(in): time step (s)
  in_varSubstep % dtInit         = dtInit                 ! intent(in): initial time step (s)
  in_varSubstep % dt_min         = dt_min                 ! intent(in): minimum time step (s)
  in_varSubstep % whole_step     = whole_step             ! intent(in): length of whole step for surface drainage and average flux
  in_varSubstep % nSubset        = nSubset                ! intent(in): total number of variables in the state subset
  in_varSubstep % doAdjustTemp   = doAdjustTemp           ! intent(in): flag to indicate if we adjust the temperature
  in_varSubstep % firstSubStep   = firstSubStep           ! intent(in): flag to denote first sub-step
  in_varSubstep % computeVegFlux = computeVegFlux         ! intent(in): flag to denote if computing energy flux over vegetation
  in_varSubstep % scalarSolution = (ixSolution==scalar)   ! intent(in): flag to denote computing the scalar solution
  in_varSubstep % iStateSplit    = iStateSplit            ! intent(in): index of the layer in the splitting operation
  in_varSubstep % fluxMask       = fluxMask               ! intent(in): mask for the fluxes used in this given state subset
 end subroutine initialize_in_varSubstep

 subroutine initialize_io_varSubstep(io_varSubstep,firstFluxCall,fluxCount,ixSaturation)
  class(io_type_varSubstep),intent(out) :: io_varSubstep  ! class object for intent(inout) varSubstep arguments
  logical(lgt),intent(in)               :: firstFluxCall  ! flag to indicate if we are processing the first flux call
  type(var_ilength),intent(in)          :: fluxCount      ! number of times fluxes are updated (should equal nsubstep)
  integer(i4b),intent(in)               :: ixSaturation   ! index of the lowest saturated layer (NOTE: only computed on the first iteration)

  ! intent(inout) arguments
  io_varSubstep % firstFluxCall = firstFluxCall           ! intent(inout): flag to indicate if we are processing the first flux call
  io_varSubstep % fluxCount     = fluxCount               ! intent(inout): number of times fluxes are updated (should equal nsubstep)
  io_varSubstep % ixSaturation  = ixSaturation            ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
 end subroutine initialize_io_varSubstep

 subroutine finalize_io_varSubstep(io_varSubstep,firstFluxCall,fluxCount,ixSaturation)
  class(io_type_varSubstep),intent(in)  :: io_varSubstep  ! class object for intent(inout) varSubstep arguments
  logical(lgt),intent(out)              :: firstFluxCall  ! flag to indicate if we are processing the first flux call
  type(var_ilength),intent(out)         :: fluxCount      ! number of times fluxes are updated (should equal nsubstep)
  integer(i4b),intent(out)              :: ixSaturation   ! index of the lowest saturated layer (NOTE: only computed on the first iteration)

  ! intent(inout) arguments
  firstFluxCall = io_varSubstep % firstFluxCall           ! intent(inout): flag to indicate if we are processing the first flux call
  fluxCount     = io_varSubstep % fluxCount               ! intent(inout): number of times fluxes are updated (should equal nsubstep)
  ixSaturation  = io_varSubstep % ixSaturation            ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
 end subroutine finalize_io_varSubstep

 subroutine finalize_out_varSubstep(out_varSubstep,dtMultiplier,nSubsteps,failedMinimumStep,reduceCoupledStep,tooMuchMelt,err,cmessage)
  class(out_type_varSubstep),intent(in) :: out_varSubstep    ! class object for intent(out) varSubstep arguments
  real(rkind),intent(out)               :: dtMultiplier      ! substep multiplier (-)
  integer(i4b),intent(out)              :: nSubsteps         ! number of substeps taken for a given split
  logical(lgt),intent(out)              :: failedMinimumStep ! flag for failed substeps
  logical(lgt),intent(out)              :: reduceCoupledStep ! flag to reduce the length of the coupled step
  logical(lgt),intent(out)              :: tooMuchMelt       ! flag to denote that ice is insufficient to support melt
  integer(i4b),intent(out)              :: err               ! error code
  character(*),intent(out)              :: cmessage          ! error message                                          

  ! intent(out) arguments
  dtMultiplier      = out_varSubstep % dtMultiplier       ! intent(out): substep multiplier (-)
  nSubsteps         = out_varSubstep % nSubsteps          ! intent(out): number of substeps taken for a given split
  failedMinimumStep = out_varSubstep % failedMinimumStep  ! intent(out): flag for failed substeps
  reduceCoupledStep = out_varSubstep % reduceCoupledStep  ! intent(out): flag to reduce the length of the coupled step
  tooMuchMelt       = out_varSubstep % tooMuchMelt        ! intent(out): flag to denote that ice is insufficient to support melt
  err               = out_varSubstep % err                ! intent(out): error code
  cmessage          = out_varSubstep % cmessage           ! intent(out): error message                                          
 end subroutine finalize_out_varSubstep
 ! **** end varSubstep ****

 ! **** computJacob ****
subroutine initialize_in_computJacob(in_computJacob,dt,nSnow,nSoil,nLayers,computeVegFlux,computeBaseflow,ixMatrix)
  class(in_type_computJacob),intent(out) :: in_computJacob           ! class object for intent(in) computJacob arguments
  real(rkind),intent(in)              :: dt                          ! intent(in): length of the time step (seconds)
  integer(i4b),intent(in)             :: nSnow                       ! intent(in): number of snow layers
  integer(i4b),intent(in)             :: nSoil                       ! intent(in): number of soil layers
  integer(i4b),intent(in)             :: nLayers                     ! intent(in): total number of layers in the snow+soil domain
  logical(lgt),intent(in)             :: computeVegFlux              ! intent(in): flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)             :: computeBaseflow             ! intent(in): flag to indicate if computing baseflow
  integer(i4b),intent(in)             :: ixMatrix                    ! intent(in): form of the Jacobian matrix                         
 
  ! intent(in) arguments
  in_computJacob % dt               =  dt                            ! intent(in): length of the time step (seconds)                    
  in_computJacob % nSnow            =  nSnow                         ! intent(in): number of snow layers
  in_computJacob % nSoil            =  nSoil                         ! intent(in): number of soil layers
  in_computJacob % nLayers          =  nLayers                       ! intent(in): total number of layers in the snow+soil domain
  in_computJacob % computeVegFlux   =  computeVegFlux                ! intent(in): flag to indicate if computing fluxes over vegetation
  in_computJacob % computeBaseflow  =  computeBaseflow               ! intent(in): flag to indicate if computing baseflow
  in_computJacob % ixMatrix         =  ixMatrix                      ! intent(in): form of the Jacobian matrix                         
 end subroutine initialize_in_computJacob

 subroutine finalize_out_computJacob(out_computJacob,err,cmessage)
  class(out_type_computJacob),intent(in) :: out_computJacob           ! class object for intent(out) computJacob arguments
  integer(i4b),intent(out)               :: err                       ! intent(out): error code
  character(*),intent(out)               :: cmessage                  ! intent(out): error message
  ! intent(out) arguments
  err               = out_computJacob % err                           ! intent(out): error code
  cmessage          = out_computJacob % cmessage                      ! intent(out): error message                                          
 end subroutine finalize_out_computJacob

 ! **** lineSearchRefinement ****
 subroutine initialize_in_lineSearchRefinement(in_lineSearchRefinement,doSearch,fOld)
  class(in_type_lineSearchRefinement),intent(out) :: in_lineSearchRefinement   ! class object for intent(out) arguments
  logical(lgt),intent(in)                         :: doSearch                  ! intent(in): flag to do the line search
  real(rkind) ,intent(in)                         :: fOld                      ! intent(in): old function value
  in_lineSearchRefinement % doSearch     = doSearch                  ! intent(in): flag to do the line search
  in_lineSearchRefinement % fOld         = fOld                      ! intent(in): old function value
 end subroutine initialize_in_lineSearchRefinement
 
 subroutine finalize_out_lineSearchRefinement(out_lineSearchRefinement,fNew,converged,err,message)
  class(out_type_lineSearchRefinement),intent(in) :: out_lineSearchRefinement   ! class object for intent(out) arguments
  real(rkind) ,intent(out)   :: fNew                                  ! intent(out): new function evaluation
  logical(lgt),intent(out)   :: converged                             ! intent(out): convergence flag
  integer(i4b),intent(out)   :: err                                   ! intent(out): error code
  character(*),intent(out)   :: message                               ! intent(out): error message
  fNew      = out_lineSearchRefinement % fNew                         ! intent(out): new function evaluation
  converged = out_lineSearchRefinement % converged                    ! intent(out): convergence flag
  err       = out_lineSearchRefinement % err                          ! intent(out): error code
  message   = out_lineSearchRefinement % message                      ! intent(out): error message
 end subroutine finalize_out_lineSearchRefinement

 ! **** summaSolv4homegrown ****

 subroutine initialize_in_summaSolv4homegrown(in_SS4NR,dt_cur,dt,iter,nSnow,nSoil,nLayers,nLeadDim,nState,ixMatrix,firstSubStep,computeVegFlux,scalarSolution,fOld)
  class(in_type_summaSolv4homegrown),intent(out)    :: in_SS4NR   ! class object for intent(out) arguments
  real(rkind) ,intent(in) :: dt_cur                   ! intent(in): current stepsize
  real(rkind) ,intent(in) :: dt                       ! intent(in): entire time step for drainage pond rate
  integer(i4b),intent(in) :: iter                     ! intent(in): iteration index
  integer(i4b),intent(in) :: nSnow                    ! intent(in): number of snow layers
  integer(i4b),intent(in) :: nSoil                    ! intent(in): number of soil layers
  integer(i4b),intent(in) :: nLayers                  ! intent(in): total number of layers
  integer(i4b),intent(in) :: nLeadDim                 ! intent(in): length of the leading dimension of the Jacobian matrix (nBands or nState)
  integer(i4b),intent(in) :: nState                   ! intent(in): total number of state variables
  integer(i4b),intent(in) :: ixMatrix                 ! intent(in): type of matrix (full or band diagonal)
  logical(lgt),intent(in) :: firstSubStep             ! intent(in): flag to indicate if we are processing the first sub-step
  logical(lgt),intent(in) :: computeVegFlux           ! intent(in): flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in) :: scalarSolution           ! intent(in): flag to denote if implementing the scalar solution
  real(rkind) ,intent(in) :: fOld                     ! intent(in): old function evaluation

  in_SS4NR % dt_cur         = dt_cur        
  in_SS4NR % dt             = dt            
  in_SS4NR % iter           = iter          
  in_SS4NR % nSnow          = nSnow         
  in_SS4NR % nSoil          = nSoil         
  in_SS4NR % nLayers        = nLayers       
  in_SS4NR % nLeadDim       = nLeadDim      
  in_SS4NR % nState         = nState        
  in_SS4NR % ixMatrix       = ixMatrix      
  in_SS4NR % firstSubStep   = firstSubStep  
  in_SS4NR % computeVegFlux = computeVegFlux 
  in_SS4NR % scalarSolution = scalarSolution
  in_SS4NR % fOld           = fOld            
 end subroutine initialize_in_summaSolv4homegrown

 subroutine initialize_io_summaSolv4homegrown(io_SS4NR,firstFluxCall,xMin,xMax,ixSaturation)
  class(io_type_summaSolv4homegrown),intent(out)    :: io_SS4NR   ! class object for intent(inout) arguments
  logical(lgt),intent(in) :: firstFluxCall ! intent(inout): flag to indicate if we are processing the first flux call
  real(rkind) ,intent(in) :: xMin,xMax     ! intent(inout): brackets of the root
  integer(i4b),intent(in) :: ixSaturation  ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)

  io_SS4NR % firstFluxCall = firstFluxCall
  io_SS4NR % xMin          = xMin    
  io_SS4NR % xMax          = xMax    
  io_SS4NR % ixSaturation  = ixSaturation 
 end subroutine initialize_io_summaSolv4homegrown

 subroutine finalize_io_summaSolv4homegrown(io_SS4NR,firstFluxCall,xMin,xMax,ixSaturation)
  class(io_type_summaSolv4homegrown),intent(in)    :: io_SS4NR   ! class object for intent(inout) arguments
  logical(lgt),intent(out) :: firstFluxCall ! intent(inout): flag to indicate if we are processing the first flux call
  real(rkind) ,intent(out) :: xMin,xMax     ! intent(inout): brackets of the root
  integer(i4b),intent(out) :: ixSaturation  ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)

  firstFluxCall = io_SS4NR % firstFluxCall
  xMin          = io_SS4NR % xMin    
  xMax          = io_SS4NR % xMax    
  ixSaturation  = io_SS4NR % ixSaturation 
 end subroutine finalize_io_summaSolv4homegrown

 subroutine finalize_out_summaSolv4homegrown(out_SS4NR,fNew,converged,err,message)
  class(out_type_summaSolv4homegrown),intent(in)    :: out_SS4NR   ! class object for intent(out) arguments
  real(rkind) ,intent(out) :: fNew      ! intent(out): new function evaluation
  logical(lgt),intent(out) :: converged ! intent(out): convergence flag
  integer(i4b),intent(out) :: err       ! intent(out): error code
  character(*),intent(out) :: message   ! intent(out): error message

  fNew      = out_SS4NR % fNew       
  converged = out_SS4NR % converged
  err       = out_SS4NR % err      
  message   = out_SS4NR % message  
 end subroutine finalize_out_summaSolv4homegrown

END MODULE data_types
