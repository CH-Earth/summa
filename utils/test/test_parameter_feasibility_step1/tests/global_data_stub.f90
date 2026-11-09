! Test-only globalData holder. Uses real SUMMA mapping/metadata types.
! Never add this substitute module to a production SUMMA build.
module globalData
  use data_types, only: gru2hru_map,var_info
  use var_lookup, only: maxvarMpar
  implicit none
  type(gru2hru_map), allocatable :: gru_struc(:)
  type(var_info) :: mpar_meta(maxvarMpar)
end module globalData
