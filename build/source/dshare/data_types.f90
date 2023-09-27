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
 USE nrtype, integerMissing=>nr_integerMissing
 USE var_lookup,only:maxvarFreq
 USE var_lookup,only:maxvarStat
 USE var_lookup,only:iLookFLUX        ! lookup indices for flux data
 USE var_lookup,only:iLookDERIV       ! lookup indices for derivative data
 implicit none
 ! constants necessary for variable defs
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
  character(len=64)                      :: varname   = 'empty'           ! variable name
  character(len=128)                     :: vardesc   = 'empty'           ! variable description
  character(len=64)                      :: varunit   = 'empty'           ! variable units
  integer(i4b)                           :: vartype   = integerMissing    ! variable type
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
  integer(8)                             :: hru_id                        ! id (non-sequential number) of the hru
  integer(i4b)                           :: nSnow                         ! number of snow layers
  integer(i4b)                           :: nSoil                         ! number of soil layers
 endtype hru_info

 ! define mapping from GRUs to the HRUs
 type, public :: gru2hru_map
  integer(8)                             :: gru_id                        ! id of the gru
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
 ! define derived types to hold multivariate data for a single variable (different variables have different length)
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
  integer(8),allocatable                 :: dat(:)                        ! dat(:)
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
  integer(8),allocatable                 :: var(:)                        ! var(:)
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
  integer(8),allocatable                 :: hru(:)                        ! hru(:)
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
 ! ** integer type of variable length (8 byte)
 type, public :: hru_int8Vec
  type(var_i8length),allocatable         :: hru(:)                        ! hru(:)%var(:)%dat
 endtype hru_int8Vec
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
 ! ** integer type of variable length (8 byte)
 type, public :: gru_int8Vec
  type(var_i8length),allocatable         :: gru(:)                        ! gru(:)%var(:)%dat
 endtype gru_int8Vec
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
 ! ** integer type of variable length (8 byte)
 type, public :: gru_hru_int8Vec
  type(hru_int8Vec),allocatable          :: gru(:)                        ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_int8Vec
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

 ! define derived types used to simplify passing subroutine arguments
 ! ** vegNrgFlux
 type, public :: in_type_vegNrgFlux ! derived type for intent(in) arguments in vegNrgFlux call
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
 end type in_type_vegNrgFlux

 type, public :: out_type_vegNrgFlux ! derived type for intent(out) arguments in vegNrgFlux call
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
  character(:),allocatable :: cmessage                                ! intent(out): error message
 end type out_type_vegNrgFlux
 ! ** end vegNrgFlux

 ! ** ssdNrgFlux
 type, public :: in_type_ssdNrgFlux ! derived type for intent(in) arguments in ssdNrgFlux call
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
   procedure :: initialize => initialize_in_ssdNrgFlux
 end type in_type_ssdNrgFlux

 type, public :: io_type_ssdNrgFlux ! derived type for intent(inout) arguments in ssdNrgFlux call
   real(rkind)              :: dGroundNetFlux_dGroundTemp        ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  contains
   procedure :: initialize => initialize_io_ssdNrgFlux
   procedure :: finalize   => finalize_io_ssdNrgFlux
 end type io_type_ssdNrgFlux

 type, public :: out_type_ssdNrgFlux ! derived type for intent(inout) arguments in ssdNrgFlux call
   real(rkind), allocatable :: iLayerNrgFlux(:)                  ! intent(out): energy flux at the layer interfaces (W m-2)
   real(rkind), allocatable :: dNrgFlux_dTempAbove(:)            ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
   real(rkind), allocatable :: dNrgFlux_dTempBelow(:)            ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
   real(rkind), allocatable :: dNrgFlux_dWatAbove(:)             ! intent(out): derivatives in the flux w.r.t. water state in the layer above (J m-2 s-1 K-1)
   real(rkind), allocatable :: dNrgFlux_dWatBelow(:)             ! intent(out): derivatives in the flux w.r.t. water state in the layer below (J m-2 s-1 K-1)
   integer(i4b)             :: err                               ! intent(out): error code
   character(:),allocatable :: cmessage                          ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_ssdNrgFlux
 end type out_type_ssdNrgFlux
 ! ** end ssdNrgFlux

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
   character(:),allocatable :: cmessage                          ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_vegLiqFlux
 end type out_type_vegLiqFlux
 ! ** end vegLiqFlux

 ! ** snowLiqFlx
 type, public :: in_type_snowLiqFlx ! class for intent(in) arguments in snowLiqFlx call
   integer(i4b)             :: nSnow                             ! intent(in):    number of snow layers
   logical(lgt)             :: firstFluxCall                     ! intent(in):    the first flux call (compute variables that are constant over the iterations)
   logical(lgt)             :: scalarSolution                    ! intent(in):    flag to indicate the scalar solution
   real(rkind)              :: scalarThroughfallRain             ! intent(in):    rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
   real(rkind)              :: scalarCanopyLiqDrainage           ! intent(in):    liquid drainage from the vegetation canopy (kg m-2 s-1)
   real(rkind), allocatable :: mLayerVolFracLiqTrial(:)          ! intent(in):    trial value of volumetric fraction of liquid water at the current iteration (-)
  contains
   procedure :: initialize => initialize_in_snowLiqFlx
 end type in_type_snowLiqFlx

 type, public :: io_type_snowLiqFlx ! class for intent(inout) arguments in snowLiqFlx call
   real(rkind), allocatable :: iLayerLiqFluxSnow(:)              ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
   real(rkind), allocatable :: iLayerLiqFluxSnowDeriv(:)         ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
  contains
   procedure :: initialize => initialize_io_snowLiqFlx
   procedure :: finalize   => finalize_io_snowLiqFlx
 end type io_type_snowLiqFlx

 type, public :: out_type_snowLiqFlx ! derived type for intent(out) arguments in snowLiqFlx call
   integer(i4b)             :: err                               ! intent(out):   error code
   character(:),allocatable :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize => finalize_out_snowLiqFlx
 end type out_type_snowLiqFlx
 ! ** end snowLiqFlx

 ! ** soilLiqFlx
 type, public :: in_type_soilLiqFlx ! derived type for intent(in) arguments in soilLiqFlx call
  integer(i4b)             :: nSoil                             ! intent(in):    number of soil layers
  logical(lgt)             :: firstSplitOper                    ! intent(in):    flag indicating first flux call in a splitting operation
  logical(lgt)             :: scalarSolution                    ! intent(in):    flag to indicate the scalar solution
  logical(lgt)             :: deriv_desired                     ! intent(in):    flag indicating if derivatives are desired
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
  real(rkind)              :: above_soilLiqFluxDeriv            ! intent(in):    derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
  real(rkind)              :: above_soildLiq_dTk                ! intent(in):    derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
  real(rkind)              :: above_soilFracLiq                 ! intent(in):    fraction of liquid water layer above soil (canopy or snow) (-)
  real(rkind)              :: scalarCanopyTranspiration         ! intent(in):    canopy transpiration (kg m-2 s-1)
  real(rkind)              :: scalarGroundEvaporation           ! intent(in):    ground evaporation (kg m-2 s-1)
  real(rkind)              :: scalarRainPlusMelt                ! intent(in):    rain plus melt (m s-1)
 end type in_type_soilLiqFlx

 type, public :: io_type_soilLiqFlx ! derived type for intent(inout) arguments in soilLiqFlx call
  real(rkind)              :: scalarMaxInfilRate                ! intent(inout): maximum infiltration rate (m s-1)
  real(rkind)              :: scalarInfilArea                   ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
  real(rkind)              :: scalarFrozenArea                  ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
  real(rkind)              :: scalarSurfaceRunoff               ! intent(inout): surface runoff (m s-1)
  real(rkind), allocatable :: mLayerdTheta_dPsi(:)              ! intent(inout): derivative in the soil water characteristic w.r.t. psi (m-1)
  real(rkind), allocatable :: mLayerdPsi_dTheta(:)              ! intent(inout): derivative in the soil water characteristic w.r.t. theta (m)
  real(rkind), allocatable :: dHydCond_dMatric(:)               ! intent(inout): derivative in hydraulic conductivity w.r.t matric head (s-1)
  real(rkind)              :: scalarInfiltration                ! intent(inout): surface infiltration rate (m s-1) -- controls on infiltration only computed for iter==1
  real(rkind), allocatable :: iLayerLiqFluxSoil(:)              ! intent(inout): liquid fluxes at layer interfaces (m s-1)
  real(rkind), allocatable :: mLayerTranspire(:)                ! intent(inout): transpiration loss from each soil layer (m s-1)
  real(rkind), allocatable :: mLayerHydCond(:)                  ! intent(inout): hydraulic conductivity in each layer (m s-1)
  real(rkind), allocatable :: dq_dHydStateAbove(:)              ! intent(inout): derivatives in the flux w.r.t. matric head in the layer above (s-1)
  real(rkind), allocatable :: dq_dHydStateBelow(:)              ! intent(inout): derivatives in the flux w.r.t. matric head in the layer below (s-1)
  real(rkind), allocatable :: dq_dHydStateLayerSurfVec(:)       ! intent(inout): derivative in surface infiltration w.r.t. hydrology state in above soil snow or canopy and every soil layer  (m s-1 or s-1)
  real(rkind), allocatable :: dq_dNrgStateAbove(:)              ! intent(inout): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
  real(rkind), allocatable :: dq_dNrgStateBelow(:)              ! intent(inout): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
  real(rkind), allocatable :: dq_dNrgStateLayerSurfVec(:)       ! intent(inout): derivative in surface infiltration w.r.t. energy state in above soil snow or canopy and every soil layer (m s-1 K-1)
  real(rkind), allocatable :: mLayerdTrans_dTCanair(:)          ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
  real(rkind), allocatable :: mLayerdTrans_dTCanopy(:)          ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
  real(rkind), allocatable :: mLayerdTrans_dTGround(:)          ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. ground temperature
  real(rkind), allocatable :: mLayerdTrans_dCanWat(:)           ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy total water  
 end type io_type_soilLiqFlx

 type, public :: out_type_soilLiqFlx ! derived type for intent(out) arguments in soilLiqFlx call
  integer(i4b)             :: err                               ! intent(out):   error code
  character(:),allocatable :: cmessage                          ! intent(out):   error message
 end type out_type_soilLiqFlx
 ! ** end soilLiqFlx

 ! ** groundwatr
 type, public :: in_type_groundwatr  ! derived type for intent(in) arguments in groundwatr call
   integer(i4b)             :: nSnow                             ! intent(in):    number of snow layers
   integer(i4b)             :: nSoil                             ! intent(in):    number of soil layers
   integer(i4b)             :: nLayers                           ! intent(in):    total number of layers
   logical(lgt)             :: firstFluxCall                     ! intent(in):    logical flag to compute index of the lowest saturated layer
   real(rkind), allocatable :: mLayerdTheta_dPsi(:)              ! intent(in):    derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
   real(rkind), allocatable :: mLayerMatricHeadLiqTrial(:)       ! intent(in):    liquid water matric potential (m)
   real(rkind), allocatable :: mLayerVolFracLiqTrial(:)          ! intent(in):    volumetric fraction of liquid water (-)
   real(rkind), allocatable :: mLayerVolFracIceTrial(:)          ! intent(in):    volumetric fraction of ice (-)
  contains
   procedure :: initialize => initialize_in_groundwatr
 end type in_type_groundwatr

 type, public :: io_type_groundwatr  ! derived type for intent(io) arguments in groundwatr call
   integer(i4b)             :: ixSaturation                      ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
  contains
   procedure :: initialize => initialize_io_groundwatr
   procedure :: finalize   => finalize_io_groundwatr 
 end type io_type_groundwatr

 type, public :: out_type_groundwatr ! derived type for intent(out) arguments in groundwatr call
   real(rkind), allocatable :: mLayerBaseflow(:)                 ! intent(out):   baseflow from each soil layer (m s-1)
   real(rkind), allocatable :: dBaseflow_dMatric(:,:)            ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
   integer(i4b)             :: err                               ! intent(out):   error code
   character(:),allocatable :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize => finalize_out_groundwatr
 end type out_type_groundwatr
 ! ** end groundwatr

 ! ** bigAquifer
 type, public :: in_type_bigAquifer  ! derived type for intent(in) arguments in bigAquifer call
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

 type, public :: io_type_bigAquifer  ! derived type for intent(inout) arguments in bigAquifer call
   real(rkind)              :: dAquiferTrans_dTCanair            ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
   real(rkind)              :: dAquiferTrans_dTCanopy            ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
   real(rkind)              :: dAquiferTrans_dTGround            ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
   real(rkind)              :: dAquiferTrans_dCanWat             ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
  contains
   procedure :: initialize => initialize_io_bigAquifer
   procedure :: finalize   => finalize_io_bigAquifer
 end type io_type_bigAquifer

 type, public :: out_type_bigAquifer  ! derived type for intent(out) arguments in bigAquifer call
   real(rkind)              :: scalarAquiferTranspire            ! intent(out):   transpiration loss from the aquifer (m s-1)
   real(rkind)              :: scalarAquiferRecharge             ! intent(out):   recharge to the aquifer (m s-1)
   real(rkind)              :: scalarAquiferBaseflow             ! intent(out):   total baseflow from the aquifer (m s-1)
   real(rkind)              :: dBaseflow_dAquifer                ! intent(out):   change in baseflow flux w.r.t. aquifer storage (s-1)
   integer(i4b)             :: err                               ! intent(out):   error code
   character(:),allocatable :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize   => finalize_out_bigAquifer
 end type out_type_bigAquifer
 ! ** end bigAquifer

contains

 ! ** ssdNrgFlux
 subroutine initialize_in_ssdNrgFlux(in_ssdNrgFlux,scalarSolution,firstFluxCall,mLayerTempTrial,flux_data,deriv_data)
  class(in_type_ssdNrgFlux),intent(out) :: in_ssdNrgFlux               ! class object for intent(in) ssdNrgFlux arguments
  logical(lgt),intent(in)               :: scalarSolution              ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)               :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  real(rkind),intent(in)                :: mLayerTempTrial(:)          ! trial value for temperature of each snow/soil layer (K)
  type(var_dlength),intent(in)          :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(in)          :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   scalarGroundNetNrgFlux       => flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1), & ! intent(out): [dp] net energy flux for the ground surface (W m-2)
   iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat,         & ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat,         & ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
   dThermalC_dWatAbove          => deriv_data%var(iLookDERIV%dThermalC_dWatAbove)%dat,     & ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. water state in the layer above
   dThermalC_dWatBelow          => deriv_data%var(iLookDERIV%dThermalC_dWatBelow)%dat,     & ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. water state in the layer above
   dThermalC_dTempAbove         => deriv_data%var(iLookDERIV%dThermalC_dTempAbove)%dat,    & ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. energy state in the layer above
   dThermalC_dTempBelow         => deriv_data%var(iLookDERIV%dThermalC_dTempBelow)%dat     ) ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. energy state in the layer above
   ! intent(in) arguments
   in_ssdNrgFlux % scalarSolution=scalarSolution .and. .not.firstFluxCall ! intent(in): flag to denote if implementing the scalar solution
   in_ssdNrgFlux % scalarGroundNetNrgFlux=scalarGroundNetNrgFlux          ! intent(in): net energy flux for the ground surface (W m-2)
   in_ssdNrgFlux % iLayerLiqFluxSnow=iLayerLiqFluxSnow                    ! intent(in): liquid flux at the interface of each snow layer (m s-1)
   in_ssdNrgFlux % iLayerLiqFluxSoil=iLayerLiqFluxSoil                    ! intent(in): liquid flux at the interface of each soil layer (m s-1)
   in_ssdNrgFlux % mLayerTempTrial=mLayerTempTrial                        ! intent(in): temperature in each layer at the current iteration (m)
   in_ssdNrgFlux % dThermalC_dWatAbove=dThermalC_dWatAbove                ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
   in_ssdNrgFlux % dThermalC_dWatBelow=dThermalC_dWatBelow                ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
   in_ssdNrgFlux % dThermalC_dTempAbove=dThermalC_dTempAbove              ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
   in_ssdNrgFlux % dThermalC_dTempBelow=dThermalC_dTempBelow              ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
  end associate
 end subroutine initialize_in_ssdNrgFlux

 subroutine initialize_io_ssdNrgFlux(io_ssdNrgFlux,deriv_data)
  class(io_type_ssdNrgFlux),intent(out) :: io_ssdNrgFlux                 ! class object for intent(inout) ssdNrgFlux arguments
  type(var_dlength),intent(in)          :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   dGroundNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp)%dat(1) ) ! intent(out): [dp] derivative in net ground flux w.r.t. ground temperature
   ! intent(inout) arguments
   io_ssdNrgFlux % dGroundNetFlux_dGroundTemp=dGroundNetFlux_dGroundTemp ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  end associate
 end subroutine initialize_io_ssdNrgFlux

 subroutine finalize_io_ssdNrgFlux(io_ssdNrgFlux,deriv_data)
  class(io_type_ssdNrgFlux),intent(in)  :: io_ssdNrgFlux                 ! class object for intent(inout) ssdNrgFlux arguments
  type(var_dlength),intent(inout)       :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
   dGroundNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp)%dat(1) ) ! intent(out): [dp] derivative in net ground flux w.r.t. ground temperature
   ! intent(inout) arguments
   dGroundNetFlux_dGroundTemp=io_ssdNrgFlux % dGroundNetFlux_dGroundTemp ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  end associate
 end subroutine finalize_io_ssdNrgFlux

 subroutine finalize_out_ssdNrgFlux(out_ssdNrgFlux,flux_data,deriv_data,err,cmessage)
  class(out_type_ssdNrgFlux),intent(in) :: out_ssdNrgFlux              ! class object for intent(out) ssdNrgFlux arguments
  type(var_dlength),intent(inout)       :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(inout)       :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from vegLiqFlux
  associate(&
   iLayerNrgFlux                => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat,         & ! intent(out): [dp(0:)] vertical energy flux at the interface of snow and soil layers
   dNrgFlux_dTempAbove          => deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove)%dat, & ! intent(out): [dp(:)] derivatives in the flux w.r.t. temperature in the layer above
   dNrgFlux_dTempBelow          => deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow)%dat, & ! intent(out): [dp(:)] derivatives in the flux w.r.t. temperature in the layer below
   dNrgFlux_dWatAbove           => deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove)%dat,  & ! intent(out):  [dp(:)] derivatives in the flux w.r.t. water state in the layer above
   dNrgFlux_dWatBelow           => deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow)%dat   ) ! intent(out): [dp(:)] derivatives in the flux w.r.t. water state in the layer below
   ! intent(out) arguments
   iLayerNrgFlux      =out_ssdNrgFlux % iLayerNrgFlux                     ! intent(out): energy flux at the layer interfaces (W m-2)
   dNrgFlux_dTempAbove=out_ssdNrgFlux % dNrgFlux_dTempAbove               ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
   dNrgFlux_dTempBelow=out_ssdNrgFlux % dNrgFlux_dTempBelow               ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
   dNrgFlux_dWatAbove =out_ssdNrgFlux % dNrgFlux_dWatAbove                ! intent(out): derivatives in the flux w.r.t. water state in the layer above (J m-2 s-1 K-1)
   dNrgFlux_dWatBelow =out_ssdNrgFlux % dNrgFlux_dWatBelow                ! intent(out): derivatives in the flux w.r.t. water state in the layer below (J m-2 s-1 K-1)
   err                =out_ssdNrgFlux % err                               ! intent(out): error code
   cmessage           =out_ssdNrgFlux % cmessage                          ! intent(out): error message
  end associate
 end subroutine finalize_out_ssdNrgFlux
 ! ** end ssdNrgFlux
 
 ! ** vegLiqFlux
 subroutine initialize_in_vegLiqFlux(in_vegLiqFlux,computeVegFlux,scalarCanopyLiqTrial,flux_data) ! SJT
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
 ! ** end vegLiqFlux

 ! ** snowLiqFlx
 subroutine initialize_in_snowLiqFlx(in_snowLiqFlx,nSnow,firstFluxCall,scalarSolution,mLayerVolFracLiqTrial,flux_data)
  class(in_type_snowLiqFlx),intent(out)   :: in_snowLiqFlx               ! class object of intent(in) arguments            
  integer(i4b),intent(in)                 :: nSnow                       ! number of snow layers
  logical(lgt),intent(in)                 :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(in)                 :: scalarSolution              ! flag to denote if implementing the scalar solution
  real(rkind),intent(in)                  :: mLayerVolFracLiqTrial(:)    ! trial value for volumetric fraction of liquid water (-)
  type(var_dlength),intent(in)            :: flux_data                   ! model fluxes for a local HRU
  associate(&
   scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),         & ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)) ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  ! intent(in) arguments
  in_snowLiqFlx % nSnow                  =nSnow                          ! intent(in): number of snow layers
  in_snowLiqFlx % firstFluxCall          =firstFluxCall                  ! intent(in): the first flux call (compute variables that are constant over the iterations)
  in_snowLiqFlx % scalarSolution         =(scalarSolution .and. .not.firstFluxCall) ! intent(in): flag to indicate the scalar solution
  in_snowLiqFlx % scalarThroughfallRain  =scalarThroughfallRain          ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
  in_snowLiqFlx % scalarCanopyLiqDrainage=scalarCanopyLiqDrainage        ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
  in_snowLiqFlx % mLayerVolFracLiqTrial  =mLayerVolFracLiqTrial(1:nSnow) ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
  end associate
 end subroutine initialize_in_snowLiqFlx 

 subroutine initialize_io_snowLiqFlx(io_snowLiqFlx,flux_data,deriv_data)
  class(io_type_snowLiqFlx),intent(out)   :: io_snowLiqFlx               ! class object for intent(inout) arguments
  type(var_dlength),intent(in)            :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(in)            :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
    iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat                  ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
    iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv      )%dat)        ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
  io_snowLiqFlx % iLayerLiqFluxSnow      =iLayerLiqFluxSnow       ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
  io_snowLiqFlx % iLayerLiqFluxSnowDeriv =iLayerLiqFluxSnowDeriv  ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
  end associate
 end subroutine initialize_io_snowLiqFlx

 subroutine finalize_io_snowLiqFlx(io_snowLiqFlx,flux_data,deriv_data)
  class(io_type_snowLiqFlx),intent(in)    :: io_snowLiqFlx               ! class object for intent(inout) arguments
  type(var_dlength),intent(inout)         :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(inout)         :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  associate(&
    iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat                  ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
    iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv      )%dat)        ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
  ! intent(inout) arguments
  iLayerLiqFluxSnow     =io_snowLiqFlx % iLayerLiqFluxSnow        ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
  iLayerLiqFluxSnowDeriv=io_snowLiqFlx % iLayerLiqFluxSnowDeriv   ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
  end associate
 end subroutine finalize_io_snowLiqFlx

 subroutine finalize_out_snowLiqFlx(out_snowLiqFlx,err,cmessage)
  class(out_type_snowLiqFlx),intent(in)   :: out_snowLiqFlx              ! class object for intent(out) arguments
  integer(i4b),intent(out)                :: err                         ! error code
  character(*),intent(out)                :: cmessage                    ! error message from vegLiqFlux
  ! intent(out) arguments
  err     =out_snowLiqFlx % err                                   ! intent(out):   error code
  cmessage=out_snowLiqFlx % cmessage                              ! intent(out):   error message
 end subroutine finalize_out_snowLiqFlx
 ! ** end snowLiqFlx

 ! ** groundwatr
 subroutine initialize_in_groundwatr(in_groundwatr,nSnow,nSoil,nLayers,firstFluxCall,mLayerMatricHeadLiqTrial,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,deriv_data)
  class(in_type_groundwatr),intent(out) :: in_groundwatr               ! class object for intent(in) arguments
  integer(i4b),intent(in)               :: nSnow                       ! number of snow layers
  integer(i4b),intent(in)               :: nSoil                       ! number of soil layers
  integer(i4b),intent(in)               :: nLayers                     ! total number of layers
  logical(lgt),intent(in)               :: firstFluxCall               ! logical flag to compute index of the lowest saturated layer
  real(rkind),intent(in)                :: mLayerMatricHeadLiqTrial(:) ! trial value for the liquid water matric potential (m)
  real(rkind),intent(in)                :: mLayerVolFracLiqTrial(:)    ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in)                :: mLayerVolFracIceTrial(:)    ! trial value for volumetric fraction of ice (-)
  type(var_dlength),intent(in)          :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
 
  associate(&
   mLayerdTheta_dPsi            => deriv_data%var(iLookDERIV%mLayerdTheta_dPsi)%dat )! intent(out): [dp(:)] derivative in the soil water characteristic w.r.t. psi
   ! intent(in) arguments
   in_groundwatr % nSnow                    = nSnow                                  ! intent(in):    number of snow layers
   in_groundwatr % nSoil                    = nSoil                                  ! intent(in):    number of soil layers
   in_groundwatr % nLayers                  = nLayers                                ! intent(in):    total number of layers
   in_groundwatr % firstFluxCall            = firstFluxCall                          ! intent(in):    logical flag to compute index of the lowest saturated layer
   in_groundwatr % mLayerdTheta_dPsi        = mLayerdTheta_dPsi                      ! intent(in):    derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
   in_groundwatr % mLayerMatricHeadLiqTrial = mLayerMatricHeadLiqTrial               ! intent(in):    liquid water matric potential (m)
   in_groundwatr % mLayerVolFracLiqTrial    = mLayerVolFracLiqTrial(nSnow+1:nLayers) ! intent(in):    volumetric fraction of liquid water (-)
   in_groundwatr % mLayerVolFracIceTrial    = mLayerVolFracIceTrial(nSnow+1:nLayers) ! intent(in):    volumetric fraction of ice (-)
  end associate
 end subroutine initialize_in_groundwatr

 subroutine initialize_io_groundwatr(io_groundwatr,ixSaturation)
  class(io_type_groundwatr),intent(out) :: io_groundwatr ! class object for intent(inout) arguments
  integer(i4b),intent(in)               :: ixSaturation  ! index of lowest saturated layer (NOTE: only computed on the first iteration)
  ! intent(inout) arguments
  io_groundwatr % ixSaturation = ixSaturation ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
 end subroutine initialize_io_groundwatr
 
 subroutine finalize_io_groundwatr(io_groundwatr,ixSaturation)
  class(io_type_groundwatr),intent(in)  :: io_groundwatr ! class object for intent(inout) arguments
  integer(i4b),intent(out)              :: ixSaturation  ! index of lowest saturated layer (NOTE: only computed on the first iteration)
  ! intent(inout) arguments
  ixSaturation = io_groundwatr % ixSaturation ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
 end subroutine finalize_io_groundwatr

 subroutine finalize_out_groundwatr(out_groundwatr,dBaseflow_dMatric,flux_data,err,cmessage)
  class(out_type_groundwatr),intent(in) :: out_groundwatr              ! class object for intent(out) arguments
  real(rkind),intent(out)               :: dBaseflow_dMatric(:,:)      ! derivative in baseflow w.r.t. matric head (s-1)
  type(var_dlength),intent(inout)       :: flux_data                   ! model fluxes for a local HRU
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from vegLiqFlux
  associate(&
   mLayerBaseflow               => flux_data%var(iLookFLUX%mLayerBaseflow)%dat )     ! intent(out): [dp(:)]  baseflow from each soil layer (m s-1)
   ! intent(out) arguments
   mLayerBaseflow    = out_groundwatr % mLayerBaseflow                               ! intent(out):   baseflow from each soil layer (m s-1)
   dBaseflow_dMatric = out_groundwatr % dBaseflow_dMatric                            ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
   err               = out_groundwatr % err                                          ! intent(out):   error code
   cmessage          = out_groundwatr % cmessage                                     ! intent(out):   error message
  end associate
 end subroutine finalize_out_groundwatr
 ! ** end groundwatr

 ! ** bigAquifer
 subroutine initialize_in_bigAquifer(in_bigAquifer,scalarAquiferStorageTrial,flux_data,deriv_data)
  class(in_type_bigAquifer),intent(out) :: in_bigAquifer             ! class object for intent(in) arguments
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
  class(io_type_bigAquifer),intent(out) :: io_bigAquifer  ! class object for intent(inout) arguments
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
  class(io_type_bigAquifer),intent(in)  :: io_bigAquifer  ! class object for intent(inout) arguments
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
  class(out_type_bigAquifer),intent(in) :: out_bigAquifer ! class object for intent(out) arguments
  type(var_dlength),intent(inout)       :: flux_data      ! model fluxes for a local HRU
  type(var_dlength),intent(inout)       :: deriv_data     ! derivatives in model fluxes w.r.t. relevant state variables
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from vegLiqFlux
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
 ! ** end bigAquifer
END MODULE data_types
