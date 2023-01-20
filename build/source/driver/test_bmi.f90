program test_bmi

    use summa_bmi, only: initialize, update, finalize, &
                        get_output_name, get_output_units, &
                        get_num_basins, get_num_output_fields, &
                        get_latlons, get_basin_field, modelTimeStep
    use globalData,only:numtim

    implicit none

    character(64), allocatable  :: dir
    character(64)               :: fldid, fldun
    real, allocatable           :: lats(:), lons(:), Q(:)
    integer                     :: istat, iseq, idt, count, ifld, nflds, ndt

    write(0,*) ""
    write(0,*) "INITIALIZING "
    write(0,*) "-------------------------------------------------------"
    istat = initialize(dir)!, iseq)

    count = get_num_basins()
    write(0,*) "The number of basins is", count
    allocate(lats(count), lons(count), Q(count))
    if(istat/=0) stop istat

    call get_latlons(lats, lons)
    write(0,*) "The latitudes are",  lats
    write(0,*) "The longitudes are", lons
    nflds = get_num_output_fields()

    write(0,*) ""
    write(0,*) "BEGIN TIME LOOP"
    write(0,*) "-------------------------------------------------------"
    DO modelTimeStep = 1, numtim
        istat = update()
        if(istat/=0) stop istat
        write(0,*) "------------- TS ", modelTimeStep, "-------------"
        call get_basin_field( 1, Q)
        call get_output_name( 1, fldid)
        call get_output_units(1, fldun)
        write(0,*) fldid, " max = ", maxval(Q), ", min = ", minval(Q), " " , fldun
    enddo

    write(0,*) ""
    write(0,*) "FINALIZING"
    write(0,*) "-------------------------------------------------------"
    istat = finalize()
    if(istat/=0) stop istat
    deallocate(lats, lons, Q)
    stop 0

end program
