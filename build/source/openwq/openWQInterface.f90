! Copyright 2020, Diogo Costa (diogo.pinhodacosta@canada.ca)
! This file is part of OpenWQ model.

! This program, openWQ, is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! C bindings of the OpenWQ hydrolink (included in openWQ.f90).
! The arguments are described in OpenWQ_hydrolink.h.

interface

   function create_openwq_c() bind(C, name="create_openwq")
      use iso_c_binding
      implicit none
      type(c_ptr) :: create_openwq_c
   end function

   function openwq_decl_c(    &
      openWQ,                 &
      num_col,                &
      nCanopy_2openwq,        &
      nSnow_2openwq,          &
      nSoil_2openwq,          &
      nRunoff_2openwq,        &
      nAquifer_2openwq,       &
      nLake_2openwq,          &
      nYdirec_2openwq,        &
      colId,                  &
      colDom,                 &
      has_glacier,            &
      num_reach,              &
      reachId) bind(C, name="openwq_decl")
      use iso_c_binding
      implicit none
      integer(c_int) :: openwq_decl_c
      type(c_ptr),          intent(in), value :: openWQ
      integer(c_int),       intent(in), value :: num_col
      integer(c_int),       intent(in), value :: nCanopy_2openwq
      integer(c_int),       intent(in), value :: nSnow_2openwq
      integer(c_int),       intent(in), value :: nSoil_2openwq
      integer(c_int),       intent(in), value :: nRunoff_2openwq
      integer(c_int),       intent(in), value :: nAquifer_2openwq
      integer(c_int),       intent(in), value :: nLake_2openwq
      integer(c_int),       intent(in), value :: nYdirec_2openwq
      integer(c_long_long), intent(in)        :: colId(*)
      integer(c_int),       intent(in)        :: colDom(*)
      integer(c_int),       intent(in), value :: has_glacier
      integer(c_int),       intent(in), value :: num_reach
      integer(c_long_long), intent(in)        :: reachId(*)
   end function

   function openwq_set_reach_state_c( &
      openWQ,                         &
      n_reach,                        &
      reachVol_m3,                    &
      airTemp_K,                      &
      SWrad_Wm2,                      &
      area_m2) bind(C, name="openwq_set_reach_state")
      use iso_c_binding
      implicit none
      integer(c_int) :: openwq_set_reach_state_c
      type(c_ptr),    intent(in), value :: openWQ
      integer(c_int), intent(in), value :: n_reach
      real(c_double), intent(in)        :: reachVol_m3(*)
      real(c_double), intent(in)        :: airTemp_K(*)
      real(c_double), intent(in)        :: SWrad_Wm2(*)
      real(c_double), intent(in)        :: area_m2(*)
   end function

   function openwq_run_time_start_c(   &
      openWQ,                          &
      index_col,                       &
      nSnow_2openwq,                   &
      nLake_2openwq,                   &
      nSoil_2openwq,                   &
      simtime_summa,                   &
      soilMoist_depVar_summa_frac,     &
      soilTemp_depVar_summa_K,         &
      airTemp_depVar_summa_K,          &
      SWrad_depVar_summa_Wm2,          &
      sweWatVol_stateVar_summa_m3,     &
      lakeWatVol_stateVar_summa_m3,    &
      canopyWatVol_stateVar_summa_m3,  &
      soilWatVol_stateVar_summa_m3,    &
      aquiferWatVol_stateVar_summa_m3, &
      col_area_m2) bind(C, name="openwq_run_time_start")
      use iso_c_binding
      implicit none
      integer(c_int) :: openwq_run_time_start_c
      type(c_ptr),    intent(in), value :: openWQ
      integer(c_int), intent(in), value :: index_col
      integer(c_int), intent(in), value :: nSnow_2openwq
      integer(c_int), intent(in), value :: nLake_2openwq
      integer(c_int), intent(in), value :: nSoil_2openwq
      integer(c_int), intent(in)        :: simtime_summa(5)
      real(c_double), intent(in)        :: soilMoist_depVar_summa_frac(*)
      real(c_double), intent(in)        :: soilTemp_depVar_summa_K(*)
      real(c_double), intent(in), value :: airTemp_depVar_summa_K
      real(c_double), intent(in), value :: SWrad_depVar_summa_Wm2
      real(c_double), intent(in)        :: sweWatVol_stateVar_summa_m3(*)
      real(c_double), intent(in)        :: lakeWatVol_stateVar_summa_m3(*)
      real(c_double), intent(in), value :: canopyWatVol_stateVar_summa_m3
      real(c_double), intent(in)        :: soilWatVol_stateVar_summa_m3(*)
      real(c_double), intent(in), value :: aquiferWatVol_stateVar_summa_m3
      real(c_double), intent(in), value :: col_area_m2
   end function

   function openwq_run_space_c(   &
      openWQ,                     &
      simtime,                    &
      source, ix_s, iy_s, iz_s,   &
      recipient, ix_r, iy_r, iz_r, &
      wflux_s2r,                  &
      wmass_source,               &
      out_to_stream) bind(C, name="openwq_run_space")
      use iso_c_binding
      implicit none
      integer(c_int) :: openwq_run_space_c
      type(c_ptr),    intent(in), value :: openWQ
      integer(c_int), intent(in)        :: simtime(5)
      integer(c_int), intent(in), value :: source
      integer(c_int), intent(in), value :: ix_s
      integer(c_int), intent(in), value :: iy_s
      integer(c_int), intent(in), value :: iz_s
      integer(c_int), intent(in), value :: recipient
      integer(c_int), intent(in), value :: ix_r
      integer(c_int), intent(in), value :: iy_r
      integer(c_int), intent(in), value :: iz_r
      real(c_double), intent(in), value :: wflux_s2r
      real(c_double), intent(in), value :: wmass_source
      integer(c_int), intent(in), value :: out_to_stream
   end function

   function openwq_run_space_in_c( &
      openWQ,                      &
      simtime,                     &
      source_EWF_name,             &
      recipient, ix_r, iy_r, iz_r, &
      wflux_s2r) bind(C, name="openwq_run_space_in")
      use iso_c_binding
      implicit none
      integer(c_int) :: openwq_run_space_in_c
      type(c_ptr),       intent(in), value :: openWQ
      integer(c_int),    intent(in)        :: simtime(5)
      character(c_char), intent(in)        :: source_EWF_name(*)
      integer(c_int),    intent(in), value :: recipient
      integer(c_int),    intent(in), value :: ix_r
      integer(c_int),    intent(in), value :: iy_r
      integer(c_int),    intent(in), value :: iz_r
      real(c_double),    intent(in), value :: wflux_s2r
   end function

   function openwq_set_watervol_c( &
      openWQ, icmp, ix, iy, iz, vol_m3) bind(C, name="openwq_set_watervol")
      use iso_c_binding
      implicit none
      integer(c_int) :: openwq_set_watervol_c
      type(c_ptr),    intent(in), value :: openWQ
      integer(c_int), intent(in), value :: icmp
      integer(c_int), intent(in), value :: ix
      integer(c_int), intent(in), value :: iy
      integer(c_int), intent(in), value :: iz
      real(c_double), intent(in), value :: vol_m3
   end function

   function openwq_set_fluxvol_c( &
      openWQ, iflux, ix, iy, iz, flux_vol_m3) bind(C, name="openwq_set_fluxvol")
      use iso_c_binding
      implicit none
      integer(c_int) :: openwq_set_fluxvol_c
      type(c_ptr),    intent(in), value :: openWQ
      integer(c_int), intent(in), value :: iflux
      integer(c_int), intent(in), value :: ix
      integer(c_int), intent(in), value :: iy
      integer(c_int), intent(in), value :: iz
      real(c_double), intent(in), value :: flux_vol_m3
   end function

   function openwq_run_time_end_c( &
      openWQ,                      &
      simtime) bind(C, name="openwq_run_time_end")
      use iso_c_binding
      implicit none
      integer(c_int) :: openwq_run_time_end_c
      type(c_ptr),    intent(in), value :: openWQ
      integer(c_int), intent(in)        :: simtime(5)
   end function

end interface
