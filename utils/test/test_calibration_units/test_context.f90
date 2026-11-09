! Review-only holders: actual parameter metadata types, synthetic defaults.
module globalData
  use nr_type, only: i4b,rkind,integerMissing=>nr_integerMissing
  use data_types, only: gru2hru_map,var_info,par_info,struct_info
  use var_lookup, only: maxvarMpar,maxvarBpar
  implicit none
  integer :: iulog=6
  integer,parameter :: upland=1,glacCln1=2,glacCln2=3,glacDbr=4,wetland=5
  type(gru2hru_map),allocatable :: gru_struc(:)
  type(var_info) :: mpar_meta(maxvarMpar),bpar_meta(maxvarBpar)
  type(par_info) :: localParFallback(maxvarMpar),basinParFallback(maxvarBpar)
  type(struct_info) :: structInfo(0)
end module
module summa_type
  use calibration_types, only: calib_info
  implicit none
  type :: config_info
    type(calib_info) :: calib
    logical :: is_case_root=.true.
  end type
end module
