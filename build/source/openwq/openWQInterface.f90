! OpenWQ C Interface
! This file contains the Fortran functions that are callable from C.
! These function are mapped to the C functions in defined in OpenWQ_interface.h
! and implmeneted in OpenWQ_interface.c


interface
    function create_openwq_c() bind(C, name="create_openwq")

        use iso_c_binding
        implicit none
        type(c_ptr) :: create_openwq_c

    end function

    function openwq_decl_c(     &
        openWQ,                 &
        num_hru,                &
        nCanopy_2openwq,        &
        nSnow_2openwq,          &
        nSoil_2openwq,          &
        nRunoff_2openwq,        &
        nAquifer_2openwq,       &
        y_direction) bind(C, name="openwq_decl")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_decl_c ! returns a return value of 0 (success) or -1 (failure)
        type(c_ptr), intent(in), value :: openWQ
        integer(c_int), intent(in), value  :: num_hru
        integer(c_int), intent(in), value  :: nCanopy_2openwq
        integer(c_int), intent(in), value  :: nSnow_2openwq
        integer(c_int), intent(in), value  :: nSoil_2openwq
        integer(c_int), intent(in), value  :: nAquifer_2openwq
        integer(c_int), intent(in), value  :: nRunoff_2openwq
        integer(c_int), intent(in), value  :: y_direction

    end function

    function openwq_run_time_start_c(&
        openWQ,                             &
        last_hru_flag,                      &
        hru_index,                          &
        nSnow_2openwq,                      &
        nSoil_2openwq,                      &
        simtime_summa,                      &
        soilMoist_depVar_summa_frac,        &                
        soilTemp_depVar_summa_K,            &
        airTemp_depVar_summa_K,             &
        sweWatVol_stateVar_summa_m3,        &
        canopyWatVol_stateVar_summa_m3,     &
        soilWatVol_stateVar_summa_m3,       &
        aquiferWatVol_stateVar_summa_m3) bind(C, name="openwq_run_time_start")

        use iso_c_binding
        implicit none
        integer(c_int)                       :: openwq_run_time_start_c ! returns 0 (success) or -1 (failure)
        type(c_ptr),    intent(in), value    :: openWQ
        logical(c_bool),   intent(in)        :: last_hru_flag
        integer(c_int), intent(in), value    :: hru_index
        integer(c_int), intent(in), value    :: nSnow_2openwq
        integer(c_int), intent(in), value    :: nSoil_2openwq
        integer(c_int), intent(in)           :: simtime_summa(5)
        real(c_double), intent(in)           :: soilMoist_depVar_summa_frac(nSoil_2openwq)
        real(c_double), intent(in)           :: soilTemp_depVar_summa_K(nSoil_2openwq)
        real(c_double), intent(in), value    :: airTemp_depVar_summa_K
        real(c_double), intent(in)           :: sweWatVol_stateVar_summa_m3(nSnow_2openwq)
        real(c_double), intent(in), value    :: canopyWatVol_stateVar_summa_m3
        real(c_double), intent(in)           :: soilWatVol_stateVar_summa_m3(nSoil_2openwq)
        real(c_double), intent(in), value    :: aquiferWatVol_stateVar_summa_m3

    end function

    function openwq_run_space_c(&
        openWQ, &
        simtime, &
        source,ix_s,iy_s,iz_s, &
        recipient,ix_r,iy_r,iz_r, &
        wflux_s2r, &
        wmass_source) bind(C, name="openwq_run_space")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_space_c ! returns 0 (success) or -1 (failure)
        type(c_ptr),    intent(in), value      :: openWQ
        integer(c_int), intent(in)             :: simtime(5)
        integer(c_int), intent(in), value      :: source
        integer(c_int), intent(in), value      :: ix_s
        integer(c_int), intent(in), value      :: iy_s 
        integer(c_int), intent(in), value      :: iz_s
        integer(c_int), intent(in), value      :: recipient
        integer(c_int), intent(in), value      :: ix_r
        integer(c_int), intent(in), value      :: iy_r
        integer(c_int), intent(in), value      :: iz_r
        real(c_double), intent(in), value      :: wflux_s2r
        real(c_double), intent(in), value      :: wmass_source

    end function

    function openwq_run_space_in_c( &
        openWQ, &
        simtime, &
        source_EWF_name, &
        recipient,ix_r,iy_r,iz_r, &
        wflux_s2r) bind(C, name="openwq_run_space_in")

        USE iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_space_in_c
        type(c_ptr), intent(in), value         :: openWQ
        integer(c_int), intent(in)             :: simtime(5)
        integer(c_int), intent(in), value      :: recipient
        integer(c_int), intent(in), value      :: ix_r
        integer(c_int), intent(in), value      :: iy_r
        integer(c_int), intent(in), value      :: iz_r
        real(c_double), intent(in), value      :: wflux_s2r
        character(c_char), intent(in)          :: source_EWF_name

    end function

    function openwq_run_time_end_c( &
        openWQ, &
        simtime) bind(C, name="openwq_run_time_end")

        USE iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_time_end_c ! returns 0 (success) or -1 (failure)
        type(c_ptr),    intent(in), value   :: openWQ
        integer(c_int), intent(in)          :: simtime(5)

    end function

end interface