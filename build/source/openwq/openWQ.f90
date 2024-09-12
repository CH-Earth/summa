module openwq
   
 USE, intrinsic :: iso_c_binding
 USE nrtype
 private
 public :: CLASSWQ_openwq

 include "openWQInterface.f90"

 type CLASSWQ_openwq
    private
    type(c_ptr) :: ptr ! pointer to openWQ class

 contains
   !  procedure :: get_num => openWQ_get_num
    procedure :: decl => openWQ_init
    procedure :: openwq_run_time_start => openwq_run_time_start
    procedure :: openwq_run_space => openwq_run_space
    procedure :: openwq_run_space_in => openwq_run_space_in
    procedure :: openwq_run_time_end => openwq_run_time_end

 end type

 interface CLASSWQ_openwq
    procedure create_openwq
 end interface
 contains
    function create_openwq()
        implicit none
        type(CLASSWQ_openwq) :: create_openwq
        create_openwq%ptr = create_openwq_c()
    end function

    ! supposed to be decl but needed to openWQ_decl in the interface file
    ! returns integer of either a failure(-1) or success(0)
   integer function openWQ_init( &
      this,                      & ! openwq object
      num_hru,                   & ! num HRU
      nCanopy_2openwq,           & ! num layers of canopy (fixed to 1)
      nSnow_2openwq,             & ! num layers of snow (fixed to max of 5 because it varies)
      nSoil_2openwq,             & ! num layers of snoil (variable)
      nRunoff_2openwq,           & ! num layers of runoff (fixed to 1)
      nAquifer_2openwq,          & ! num layers of aquifer (fixed to 1)
      nYdirec_2openwq)                 ! num of layers in y-dir (set to 1 because not used in summa)
      
      implicit none
      class(CLASSWQ_openwq) :: this
      integer(i4b), intent(in) :: num_hru
      integer(i4b), intent(in) :: nCanopy_2openwq
      integer(i4b), intent(in) :: nSnow_2openwq
      integer(i4b), intent(in) :: nSoil_2openwq
      integer(i4b), intent(in) :: nRunoff_2openwq
      integer(i4b), intent(in) :: nAquifer_2openwq
      
      integer(i4b), intent(in) :: nYdirec_2openwq

      openWQ_init = openwq_decl_c(  &
         this%ptr,                  & ! openwq object
         num_hru,                   & ! num HRU
         nCanopy_2openwq,           & ! num layers of canopy (fixed to 1)
         nSnow_2openwq,             & ! num layers of snow (fixed to max of 5 because it varies)
         nSoil_2openwq,             & ! num layers of snoil (variable)
         nRunoff_2openwq,           & ! num layers of runoff (fixed to 1)
         nAquifer_2openwq,          & ! num layers of aquifer (fixed to 1)
         nYdirec_2openwq)                 ! num of layers in y-dir (set to 1 because not used in summa)

   end function


   integer function openwq_run_time_start(   &
      this,                                  &
      last_hru_flag,                         &
      hru_index,                             &
      nSnow_2openwq,                         &
      nSoil_2openwq,                         &
      simtime,                               &
      soilMoist_depVar_summa_frac,           &                    
      soilTemp_depVar_summa_K,               &
      airTemp_depVar_summa_K,                &
      sweWatVol_stateVar_summa_m3,           &
      canopyWatVol_stateVar_summa_m3,        &
      soilWatVol_stateVar_summa_m3,          &
      aquiferWatVol_stateVar_summa_m3)
      
      implicit none
      class(CLASSWQ_openwq)      :: this
      logical(1), intent(in)     :: last_hru_flag
      integer(i4b), intent(in)   :: hru_index
      integer(i4b), intent(in)   :: nSnow_2openwq
      integer(i4b), intent(in)   :: nSoil_2openwq
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars
      real(rkind),  intent(in)   :: airTemp_depVar_summa_K
      real(rkind),  intent(in)   :: soilTemp_depVar_summa_K(nSoil_2openwq)
      real(rkind),  intent(in)   :: soilMoist_depVar_summa_frac(nSoil_2openwq)
      real(rkind),  intent(in)   :: canopyWatVol_stateVar_summa_m3
      real(rkind),  intent(in)   :: sweWatVol_stateVar_summa_m3(nSnow_2openwq)
      real(rkind),  intent(in)   :: soilWatVol_stateVar_summa_m3(nSoil_2openwq)
      real(rkind),  intent(in)   :: aquiferWatVol_stateVar_summa_m3

      openwq_run_time_start = openwq_run_time_start_c( &
         this%ptr,                              & 
         last_hru_flag,                         &
         hru_index,                             &
         nSnow_2openwq,                         &
         nSoil_2openwq,                         &
         simtime,                               &
         soilMoist_depVar_summa_frac,           &                    
         soilTemp_depVar_summa_K,               &
         airTemp_depVar_summa_K,                &
         sweWatVol_stateVar_summa_m3,           &
         canopyWatVol_stateVar_summa_m3,        &
         soilWatVol_stateVar_summa_m3,          &
         aquiferWatVol_stateVar_summa_m3)
   
   end function

   integer function openwq_run_space(  &
      this,                            &
      simtime,                         &
      source,ix_s,iy_s,iz_s,           &
      recipient,ix_r,iy_r,iz_r,        &
      wflux_s2r,wmass_source)

      implicit none
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars
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

      openwq_run_space = openwq_run_space_c( &
         this%ptr,                           &
         simtime,                            &
         source,ix_s,iy_s,iz_s,              &
         recipient,ix_r,iy_r,iz_r,           &
         wflux_s2r,wmass_source)
   
   end function

   integer function openwq_run_space_in(  &
      this,                               &
      simtime,                            &
      source_EWF_name,                    &
      recipient,ix_r,iy_r,iz_r,           &
      wflux_s2r)

      implicit none
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars
      integer(i4b), intent(in)   :: recipient
      integer(i4b), intent(in)   :: ix_r
      integer(i4b), intent(in)   :: iy_r
      integer(i4b), intent(in)   :: iz_r
      real(rkind),  intent(in)   :: wflux_s2r
      character(*), intent(in)   :: source_EWF_name

      openwq_run_space_in = openwq_run_space_in_c( &
         this%ptr,                                 &
         simtime,                                  &
         source_EWF_name,                          &
         recipient,ix_r,iy_r,iz_r,                 &
         wflux_s2r)

   end function


   integer function openwq_run_time_end(  &
      this,                               &
      simtime)

      implicit none
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars

      openwq_run_time_end = openwq_run_time_end_c( &
         this%ptr,                                 &
         simtime)

   end function

end module openwq