module calibration_types

  USE nr_type, only: i4b,i8b,rkind,lgt
  USE nr_type, only: integerMissing=>nr_integerMissing

  implicit none
  private

  public :: calib_info
  public :: calib_param_info
  public :: target_info
  public :: nsga2_info
  public :: param_bounds_info
  public :: ordered_constraint
  public :: multiplier_dependency

  ! Derived type definitions go here.

 ! -----------------------------------------------------------------------------------------------------------
 ! Configuration information for one calibration parameter.
 ! Each calibration parameter has one scalar search variable representing
 ! a physical value, spatial multiplier, or parameter ratio.
 type, public :: calib_param_info
   character(len=64) :: name = ''                    ! Parameter name
   logical(lgt)      :: has_lower = .false.          ! Lower bound explicitly configured
   logical(lgt)      :: has_upper = .false.          ! Upper bound explicitly configured
   real(rkind)       :: lower_bound = 0._rkind       ! Lower calibration bound
   real(rkind)       :: upper_bound = 0._rkind       ! Upper calibration bound
   character(len=16) :: transformation = 'none'      ! Search transformation
   integer(i4b)      :: param_type = 1               ! Physical=1, spatial multiplier=2, ratio=3
   character(len=64) :: reference_name = ''          ! Reference parameter for a ratio
   character(len=16) :: reference_reduction = 'none' ! none: elementwise; min: local reference minimum
 end type calib_param_info

 ! -----------------------------------------------------------------------------------------------------------
 type :: param_bounds_info
  character(len=64)              :: name               ! parameter name, optionally name@domainType
  real(rkind)                    :: lower              ! lower bound replacing the parameter file's
  real(rkind)                    :: upper              ! upper bound replacing the parameter file's
 end type param_bounds_info

 ! -----------------------------------------------------------------------------------------------------------
 type :: ordered_constraint
  character(len=64), allocatable :: parameters(:)      ! Ordered parameter names
  real(rkind)                    :: gap_fraction       ! Minimum adjacent gap as fraction of total range
 end type ordered_constraint

 ! -----------------------------------------------------------------------------------------------------------
 type :: multiplier_dependency
   character(len=64), allocatable :: parameters(:)    ! Multiplicative parameter dependency.
 end type multiplier_dependency

 ! -----------------------------------------------------------------------------------------------------------
 ! One calibration target: an observed series, the simulated variable it is compared against, and the
 ! metric that scores the two together.  A calibration owns a list of these and evaluates every one of
 ! them for each parameter trial, so a single-objective calibration is simply a list of length one.
 type,public  :: target_info
  character(len=64)              :: name = ''               ! target name, used to label the output
  character(len=64)              :: variable = 'streamflow' ! simulated variable compared against the observations
  ! Which spatial unit of the domain the simulated variable is scored over.  The whole domain is right
  ! for a single basin and wrong for anything else: two basins in one domain have their own storage,
  ! and a reach outlet temperature belongs to one reach rather than to an average over the basin.
  character(len=16)              :: spatial_unit = 'domain'     ! domain, gru, hru, reach
  integer(i8b)                   :: spatial_id = integerMissing ! id of that unit, as the domain numbers it
  character(len=:),  allocatable :: obs_path                ! path to the observation file
  character(len=:),  allocatable :: obs_file                ! observation file
  character(len=:),  allocatable :: vname_obs               ! name of the observed variable within that file
  character(len=16)              :: metric = 'kge'          ! kge, kgep, nse, mae, rmse
  character(len=16)              :: obs_transform = 'none'  ! none, log, power, boxcox
  real(rkind)                    :: boxcox_lambda = 0.25_rkind ! Box-Cox exponent (used only for boxcox)
  real(rkind)                    :: weight = 1.0_rkind      ! DDS weight: positive = active, zero = diagnostic
                                                            ! NSGA-II currently requires every weight > 0
  character(len=64)              :: obs_units = ''          ! units of the observations, overriding the file's
  ! Putting the simulated and observed series on the same footing.  A satellite storage product, for
  ! instance, reports a departure in millimetres from a multi-year mean, while the model carries the
  ! rate storage is changing at, so that target accumulates and references a baseline.
  logical(lgt)                   :: accumulate = .false.    ! integrate the simulated rate into its quantity
  character(len=:), allocatable  :: baseline_start          ! start of the anomaly baseline (YYYY-MM-DD)
  character(len=:), allocatable  :: baseline_end            ! end of the anomaly baseline (YYYY-MM-DD)
  character(len=16)              :: balance = ''            ! seasonal or annual balance between extremes
  real(rkind)                    :: balance_window = 60._rkind ! half-width of the search for an extreme (days)
  real(rkind)                    :: low_flow = 1._rkind     ! fraction of the lowest observed values scored
  logical(lgt)                   :: flow_duration = .false. ! score the sorted series (flow-duration curves)
 end type target_info

 ! -----------------------------------------------------------------------------------------------------------
 ! NSGA-II settings (Deb et al., 2002).  The sample budget is population_size times the number of
 ! generations, the first of which is the random initial population.
 type,public  :: nsga2_info
  integer(i4b)                   :: population_size = 100          ! members kept each generation
  real(rkind)                    :: crossover_probability = 0.9_rkind ! chance a pair of parents is crossed
  real(rkind)                    :: crossover_eta = 20._rkind      ! SBX distribution index
  real(rkind)                    :: mutation_probability = -1._rkind ! chance each variable mutates; negative is 1/nParam
  real(rkind)                    :: mutation_eta = 20._rkind       ! polynomial-mutation distribution index
 end type nsga2_info

 ! -----------------------------------------------------------------------------------------------------------
 ! calibration configuration
 type,public  :: calib_info
  character(len=:),  allocatable :: metric             ! KGE, KGEp, NSE, RMSE, MAE
  character(len=:),  allocatable :: obs_transform      ! none, log, power, boxcox
  real(rkind)                    :: boxcox_lambda = 0.25_rkind ! Default exponent inherited by targets
  character(len=:),  allocatable :: start_date         ! start of the calibration time period
  character(len=:),  allocatable :: end_date           ! end of the calibration time period
  type(calib_param_info),      allocatable :: param_list(:)   ! Parameters included in optimization
  type(param_bounds_info),     allocatable :: param_bounds(:) ! bounds that replace the parameter files'
  type(ordered_constraint),    allocatable :: ordered(:)      ! ordered parameter constraints
  type(multiplier_dependency), allocatable :: multiplier(:)   ! multiplicative parameter dependency
  type(target_info),           allocatable :: targets(:)      ! calibration targets, one per objective
  integer(i4b)                   :: n_samples      = 1000        ! number of parameter samples to evaluate
  real(rkind)                    :: default_metric = 9999._rkind ! default objective function
  character(len=16)              :: algorithm      = 'dds'       ! dds (on the scalarized targets) or nsga2 (on each)
  type(nsga2_info)               :: nsga2                        ! NSGA-II settings
  logical(lgt)                   :: write_aligned  = .false.     ! flag to write the aligned sim/obs time series
 end type calib_info


end module calibration_types
