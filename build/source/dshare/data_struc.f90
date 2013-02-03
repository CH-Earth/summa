MODULE data_struc
 ! used to define model data structures
 USE nrtype
 implicit none
 private
 ! ***********************************************************************************************************
 ! Define the model decisions
 ! ***********************************************************************************************************
 ! the model decision structure
 type model_options
  character(len=64)                      :: cOption
  character(len=64)                      :: cDecision
  integer(i4b)                           :: iDecision
 end type model_options
 type(model_options),pointer,save,public :: model_decisions(:)      ! the decision structure
 ! ***********************************************************************************************************
 ! Define metadata for model forcing datafile 
 ! ***********************************************************************************************************
 ! define a derived type for the data in the file
 type,public  :: file_info
  character(len=256)                     :: filenm=''
  integer(i4b)                           :: ncols                    ! number of columns in the file
  integer(i4b),pointer                   :: time_ix(:) => null()     ! column index for each time variable
  integer(i4b),pointer                   :: data_ix(:) => null()     ! column index for each forcing data variable
  real(dp)                               :: data_step                ! time step of the data
  integer(i4b)                           :: istart                   ! start index of the simulation
  integer(i4b)                           :: numtim                   ! number of time steps in the simulation
 end type file_info
 ! and save all the data in a single data structure
 type(file_info),pointer,save,public     :: forcFileInfo => null()   ! file info for model forcing data
 ! ***********************************************************************************************************
 ! Define metadata on model parameters
 ! ***********************************************************************************************************
 ! define a data type to store model parameter information
 type,public  :: par_info
  real(dp)                               :: default_val              ! default parameter value
  real(dp)                               :: lower_limit              ! lower bound
  real(dp)                               :: upper_limit              ! upper bound
 endtype par_info
 ! define a vector, with a separate element for each parameter (variable)
 type(par_info),pointer,save,public      :: parFallback(:) => null() ! fixed parameter structure 
 ! ***********************************************************************************************************
 ! Define variable metadata
 ! ***********************************************************************************************************
 ! define derived type for model variables, including name, decription, and units
 type,public :: var_info
  character(len=64)                      :: varname=''               ! variable name
  CHARACTER(len=128)                     :: vardesc=''               ! variable description
  character(len=64)                      :: varunit=''               ! variable units
  character(len=32)                      :: vartype=''               ! variable type (scalar, model layers, etc.)
  logical(lgt)                           :: v_write=.FALSE.          ! flag to write variable to the output file
 endtype var_info
 ! define arrays of metadata
 type(var_info),pointer,save,public      :: time_meta(:) => null()   ! model time information
 type(var_info),pointer,save,public      :: forc_meta(:) => null()   ! model forcing data
 type(var_info),pointer,save,public      :: attr_meta(:) => null()   ! local attributes
 type(var_info),pointer,save,public      :: type_meta(:) => null()   ! local classification of veg, soil, etc.
 type(var_info),pointer,save,public      :: mpar_meta(:) => null()   ! model parameters
 type(var_info),pointer,save,public      :: mvar_meta(:) => null()   ! model variables
 type(var_info),pointer,save,public      :: indx_meta(:) => null()   ! model indices
 ! ***********************************************************************************************************
 ! Define hierarchal derived data types
 ! ***********************************************************************************************************
 ! define named variables to describe the layer type
 integer(i4b),parameter,public      :: ix_soil=1001             ! named variable to denote a soil layer
 integer(i4b),parameter,public      :: ix_snow=1002             ! named variable to denote a snow layer
 integer(i4b),parameter,public      :: ix_mixd=1003             ! named variable to denote a mixed layer
 ! define derived types to hold multivariate data for a single variable (different variables have different length)
 ! NOTE: use derived types here to facilitate adding the "variable" dimension
 ! ** double precision type
 type, public :: dlength
  real(dp),pointer                       :: dat(:) => null()
 endtype dlength
 ! ** integer type
 type, public :: ilength
  integer(i4b),pointer                   :: dat(:) => null()
 endtype ilength
 ! define derived types to hold data for multiple variables
 ! NOTE: use derived types here to facilitate adding extra dimensions (e.g., spatial)
 ! ** double precision type of variable length
 type, public :: var_dlength 
  type(dlength),pointer                  :: var(:) => null()
 endtype var_dlength
 ! ** integer type of variable length
 type, public :: var_ilength
  type(ilength),pointer                  :: var(:) => null()
 endtype var_ilength
 ! ** double precision type of fixed length
 type, public :: var_d
  real(dp),pointer                       :: var(:) => null()
 endtype var_d
 ! ** integer type of variable length
 type, public :: var_i
  integer(i4b),pointer                   :: var(:) => null()
 endtype var_i
 ! define top-level derived types
 ! NOTE: either allocate directly, or use to point to higher dimensional structures
 type(var_i),pointer,save,public         :: time_data => null()      ! model time data
 type(var_d),pointer,save,public         :: forc_data => null()      ! model forcing data
 type(var_d),pointer,save,public         :: attr_data => null()      ! local attributes
 type(var_i),pointer,save,public         :: type_data => null()      ! local classification of veg, soil, etc.
 type(var_d),pointer,save,public         :: mpar_data => null()      ! model parameters
 type(var_dlength),pointer,save,public   :: mvar_data => null()      ! model variables
 type(var_ilength),pointer,save,public   :: indx_data => null()      ! model indices
 ! ***********************************************************************************************************
 ! Define common variables
 ! ***********************************************************************************************************
 real(dp),save,public                    :: refJulday                ! reference time in fractional julian days
 real(dp),save,public                    :: fracJulday               ! fractional julian days since the start of year
 integer(i4b),save,public                :: yearLength               ! number of days in the current year
 integer(i4b),save,public                :: urbanVegCategory=1       ! vegetation category for urban areas
 ! ***********************************************************************************************************
 ! Define ancillary data structures
 ! ***********************************************************************************************************
 type(var_d),pointer,save,public         :: mpar_sets(:) => null()   ! structure to hold multiple parameter sets
 type(var_i),pointer,save,public         :: refTime      => null()   ! reference time for the model simulation
 ! ***********************************************************************************************************
END MODULE data_struc
! ***********************************************************************************************************************

