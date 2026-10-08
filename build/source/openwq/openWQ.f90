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

! Fortran wrapper of the OpenWQ hydrolink class (see OpenWQ_hydrolink.h for the arguments).

module openWQ

   USE, intrinsic :: iso_c_binding
   USE nr_type
   implicit none
   private
   public :: CLASSWQ_openwq

   include "openWQInterface.f90"

   type CLASSWQ_openwq
      private
      type(c_ptr) :: ptr  ! C++ CLASSWQ_openwq object
   contains
      procedure :: decl                   => openwq_decl
      procedure :: openwq_set_reach_state => openwq_set_reach_state
      procedure :: openwq_run_time_start  => openwq_run_time_start
      procedure :: openwq_run_space       => openwq_run_space
      procedure :: openwq_run_space_in    => openwq_run_space_in
      procedure :: openwq_set_watervol    => openwq_set_watervol
      procedure :: openwq_set_fluxvol     => openwq_set_fluxvol
      procedure :: openwq_run_time_end    => openwq_run_time_end
   end type

   interface CLASSWQ_openwq
      procedure create_openwq
   end interface

contains

   function create_openwq()
      type(CLASSWQ_openwq) :: create_openwq
      create_openwq%ptr = create_openwq_c()
   end function

   integer function openwq_decl( &
      this,                      &
      num_col,                   &
      nCanopy_2openwq,           &
      nSnow_2openwq,             &
      nSoil_2openwq,             &
      nRunoff_2openwq,           &
      nAquifer_2openwq,          &
      nLake_2openwq,             &
      nYdirec_2openwq,           &
      colId,                     &
      colDom,                    &
      has_glacier,               &
      num_reach,                 &
      reachId)
      class(CLASSWQ_openwq)            :: this
      integer(i4b), intent(in)         :: num_col
      integer(i4b), intent(in)         :: nCanopy_2openwq
      integer(i4b), intent(in)         :: nSnow_2openwq
      integer(i4b), intent(in)         :: nSoil_2openwq
      integer(i4b), intent(in)         :: nRunoff_2openwq
      integer(i4b), intent(in)         :: nAquifer_2openwq
      integer(i4b), intent(in)         :: nLake_2openwq
      integer(i4b), intent(in)         :: nYdirec_2openwq
      integer(c_long_long), intent(in) :: colId(:)
      integer(c_int), intent(in)       :: colDom(:)
      integer(i4b), intent(in)         :: has_glacier
      integer(i4b), intent(in)         :: num_reach
      integer(c_long_long), intent(in) :: reachId(:)
      openwq_decl = openwq_decl_c( &
         this%ptr,                 &
         num_col,                  &
         nCanopy_2openwq,          &
         nSnow_2openwq,            &
         nSoil_2openwq,            &
         nRunoff_2openwq,          &
         nAquifer_2openwq,         &
         nLake_2openwq,            &
         nYdirec_2openwq,          &
         colId,                    &
         colDom,                   &
         has_glacier,              &
         num_reach,                &
         reachId)
   end function

   integer function openwq_set_reach_state( &
      this, n_reach, reachVol_m3, airTemp_K, SWrad_Wm2, area_m2)
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: n_reach
      real(rkind),  intent(in)   :: reachVol_m3(:)
      real(rkind),  intent(in)   :: airTemp_K(:)
      real(rkind),  intent(in)   :: SWrad_Wm2(:)
      real(rkind),  intent(in)   :: area_m2(:)
      openwq_set_reach_state = openwq_set_reach_state_c( &
         this%ptr, n_reach, reachVol_m3, airTemp_K, SWrad_Wm2, area_m2)
   end function

   integer function openwq_run_time_start(   &
      this,                                  &
      index_col,                             &
      nSnow_2openwq,                         &
      nLake_2openwq,                         &
      nSoil_2openwq,                         &
      simtime,                               &
      soilMoist_depVar_summa_frac,           &
      soilTemp_depVar_summa_K,               &
      airTemp_depVar_summa_K,                &
      SWrad_depVar_summa_Wm2,                &
      sweWatVol_stateVar_summa_m3,           &
      lakeWatVol_stateVar_summa_m3,          &
      canopyWatVol_stateVar_summa_m3,        &
      soilWatVol_stateVar_summa_m3,          &
      aquiferWatVol_stateVar_summa_m3,       &
      col_area_m2)
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: index_col     ! 0-based column index
      integer(i4b), intent(in)   :: nSnow_2openwq
      integer(i4b), intent(in)   :: nLake_2openwq
      integer(i4b), intent(in)   :: nSoil_2openwq
      integer(i4b), intent(in)   :: simtime(5)
      real(rkind),  intent(in)   :: soilMoist_depVar_summa_frac(:)
      real(rkind),  intent(in)   :: soilTemp_depVar_summa_K(:)
      real(rkind),  intent(in)   :: airTemp_depVar_summa_K
      real(rkind),  intent(in)   :: SWrad_depVar_summa_Wm2
      real(rkind),  intent(in)   :: sweWatVol_stateVar_summa_m3(:)
      real(rkind),  intent(in)   :: lakeWatVol_stateVar_summa_m3(:)
      real(rkind),  intent(in)   :: canopyWatVol_stateVar_summa_m3
      real(rkind),  intent(in)   :: soilWatVol_stateVar_summa_m3(:)
      real(rkind),  intent(in)   :: aquiferWatVol_stateVar_summa_m3
      real(rkind),  intent(in)   :: col_area_m2
      openwq_run_time_start = openwq_run_time_start_c( &
         this%ptr,                              &
         index_col,                             &
         nSnow_2openwq,                         &
         nLake_2openwq,                         &
         nSoil_2openwq,                         &
         simtime,                               &
         soilMoist_depVar_summa_frac,           &
         soilTemp_depVar_summa_K,               &
         airTemp_depVar_summa_K,                &
         SWrad_depVar_summa_Wm2,                &
         sweWatVol_stateVar_summa_m3,           &
         lakeWatVol_stateVar_summa_m3,          &
         canopyWatVol_stateVar_summa_m3,        &
         soilWatVol_stateVar_summa_m3,          &
         aquiferWatVol_stateVar_summa_m3,       &
         col_area_m2)
   end function

   integer function openwq_run_space(  &
      this,                            &
      simtime,                         &
      source, ix_s, iy_s, iz_s,        &
      recipient, ix_r, iy_r, iz_r,     &
      wflux_s2r, wmass_source,         &
      out_to_stream)
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5)
      integer(i4b), intent(in)   :: source
      integer(i4b), intent(in)   :: ix_s
      integer(i4b), intent(in)   :: iy_s
      integer(i4b), intent(in)   :: iz_s
      integer(i4b), intent(in)   :: recipient
      integer(i4b), intent(in)   :: ix_r
      integer(i4b), intent(in)   :: iy_r
      integer(i4b), intent(in)   :: iz_r
      real(rkind),  intent(in)   :: wflux_s2r
      real(rkind),  intent(in)   :: wmass_source
      integer(i4b), intent(in)   :: out_to_stream
      openwq_run_space = openwq_run_space_c( &
         this%ptr,                           &
         simtime,                            &
         source, ix_s, iy_s, iz_s,           &
         recipient, ix_r, iy_r, iz_r,        &
         wflux_s2r, wmass_source,            &
         out_to_stream)
   end function

   integer function openwq_run_space_in(  &
      this,                               &
      simtime,                            &
      source_EWF_name,                    &
      recipient, ix_r, iy_r, iz_r,        &
      wflux_s2r)
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5)
      character(*), intent(in)   :: source_EWF_name
      integer(i4b), intent(in)   :: recipient
      integer(i4b), intent(in)   :: ix_r
      integer(i4b), intent(in)   :: iy_r
      integer(i4b), intent(in)   :: iz_r
      real(rkind),  intent(in)   :: wflux_s2r
      openwq_run_space_in = openwq_run_space_in_c( &
         this%ptr,                                 &
         simtime,                                  &
         trim(source_EWF_name)//c_null_char,       &
         recipient, ix_r, iy_r, iz_r,              &
         wflux_s2r)
   end function

   integer function openwq_set_watervol( &
      this, icmp, ix, iy, iz, vol_m3)
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: icmp   ! 0-based compartment index
      integer(i4b), intent(in)   :: ix     ! 1-based cell indices
      integer(i4b), intent(in)   :: iy
      integer(i4b), intent(in)   :: iz
      real(rkind),  intent(in)   :: vol_m3
      openwq_set_watervol = openwq_set_watervol_c( &
         this%ptr, icmp, ix, iy, iz, vol_m3)
   end function

   integer function openwq_set_fluxvol( &
      this, iflux, ix, iy, iz, flux_vol_m3)
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: iflux  ! 0-based export index
      integer(i4b), intent(in)   :: ix     ! 1-based cell indices
      integer(i4b), intent(in)   :: iy
      integer(i4b), intent(in)   :: iz
      real(rkind),  intent(in)   :: flux_vol_m3
      openwq_set_fluxvol = openwq_set_fluxvol_c( &
         this%ptr, iflux, ix, iy, iz, flux_vol_m3)
   end function

   integer function openwq_run_time_end(  &
      this,                               &
      simtime)
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5)
      openwq_run_time_end = openwq_run_time_end_c( &
         this%ptr,                                 &
         simtime)
   end function

end module openWQ
