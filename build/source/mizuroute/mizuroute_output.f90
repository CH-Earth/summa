module mizuroute_output_module

  use nrtype
  use netcdf

  use mizuroute_types, only: mizuroute_info
  use mizuroute_types, only: mizuroute_domain

  implicit none
  private

  public :: define_mizuroute_output
  public :: write_mizuroute_output

contains

  !-----------------------------------------------------------------------
  ! Add mizuRoute output variables to an existing SUMMA NetCDF file
  !-----------------------------------------------------------------------
  subroutine define_mizuroute_output(ncid, info, domain, ierr, message)

    use globaldata,     only: routeMethods
    use init_mizuRoute, only: route_method_name

    integer(i4b),           intent(in)  :: ncid
    type(mizuroute_info),   intent(in)  :: info
    type(mizuroute_domain), intent(in)  :: domain
    integer(i4b),           intent(out) :: ierr
    character(*),           intent(out) :: message

    integer(i4b) :: dim_time, dim_hru, dim_seg, dim_method
    integer(i4b) :: varid_hru, varid_seg, varid_method
    integer(i4b) :: varid_uparea, varid_qbasin, varid_Qreach
    integer(i4b), dimension(2) :: dimids_basin
    integer(i4b), dimension(3) :: dimids_reach
    integer(i4b) :: iRoute

    logical(lgt) :: in_define
    integer(i4b) :: ierr_enddef

    character(len=32) :: attName

    ierr = 0
    message = 'define_mizuroute_output/'

    in_define = .false.

    netcdf_block: block

      ! enter (re)-define mode
      ierr = nf90_redef(ncid); if(ierr/=nf90_noerr) exit netcdf_block
      in_define = .true.

      ! existing SUMMA time dimension
      ierr = nf90_inq_dimid(ncid, 'time', dim_time); if(ierr/=nf90_noerr) exit netcdf_block

      ! dimensions
      ierr = nf90_def_dim(ncid, 'mizu_hru', info%n_hru,         dim_hru);    if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_def_dim(ncid, 'seg',      info%n_seg,         dim_seg);    if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_def_dim(ncid, 'method',   size(routeMethods), dim_method); if(ierr/=nf90_noerr) exit netcdf_block

      dimids_basin = (/ dim_hru, dim_time /)
      dimids_reach = (/ dim_method, dim_seg, dim_time /)

      ! upstream area
      ierr = nf90_def_var(ncid, 'upArea', NF90_DOUBLE, (/dim_seg/), varid_uparea);    if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_uparea, 'long_name', &
                      'drainage area above the downstream end of each river reach');  if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_uparea, 'units', 'm2');                         if(ierr/=nf90_noerr) exit netcdf_block

      ! basin runoff on river-network HRUs
      ierr = nf90_def_var(ncid, 'q_basin', NF90_DOUBLE, dimids_basin, varid_qbasin);  if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_qbasin, 'long_name', &
                                'runoff on mizuRoute river-network HRUs');            if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_qbasin, 'units', 'm s-1');                      if(ierr/=nf90_noerr) exit netcdf_block

      ! routed streamflow
      ierr = nf90_def_var(ncid, 'Q_reach', NF90_DOUBLE, dimids_reach, varid_Qreach);  if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_Qreach, 'long_name', &
                          'streamflow at the downstream end of each river reach');    if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_Qreach, 'units', 'm3 s-1');                     if(ierr/=nf90_noerr) exit netcdf_block

      ! coordinate variables: mizu_hru
      ierr = nf90_def_var(ncid, 'mizu_hru', NF90_INT, (/dim_hru/), varid_hru);        if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_hru, 'units', '-');                             if(ierr/=nf90_noerr) exit netcdf_block

      ! coordinate variables: seg
      ierr = nf90_def_var(ncid, 'seg', NF90_INT, (/dim_seg/), varid_seg);             if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_seg, 'units', '-');                             if(ierr/=nf90_noerr) exit netcdf_block

      ! coordinate variables: method
      ierr = nf90_def_var(ncid, 'method', NF90_INT, (/dim_method/), varid_method);    if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_method, 'long_name', 'routing method');         if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_att(ncid, varid_method, 'source', 'mizuRoute');                 if(ierr/=nf90_noerr) exit netcdf_block

      ! routing method descriptions
      do iRoute = 1,size(routeMethods)
        write(attName,'("method_",I0)') routeMethods(iRoute)
        ierr = nf90_put_att(ncid, varid_method, trim(attName), &
                            route_method_name(routeMethods(iRoute))); if(ierr/=nf90_noerr) exit netcdf_block
      enddo

      ! leave define mode
      ierr = nf90_enddef(ncid); if(ierr/=nf90_noerr) exit netcdf_block
      in_define = .false.

      ! upstream area
      ierr = nf90_put_var(ncid, varid_uparea, domain%river_network%driver%totArea); if(ierr/=nf90_noerr) exit netcdf_block

      ! coordinate data
      ierr = nf90_put_var(ncid, varid_hru,    domain%river_network%driver%hru_id);  if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_var(ncid, varid_seg,    domain%river_network%driver%seg_id);  if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_put_var(ncid, varid_method, routeMethods);                        if(ierr/=nf90_noerr) exit netcdf_block
    
    end block netcdf_block

    ! process NetCDF errors
    if(ierr/=nf90_noerr)then
      message = trim(message)//trim(nf90_strerror(ierr))
      if(in_define) ierr_enddef = nf90_enddef(ncid)
      return
    endif

    ierr = 0

  end subroutine define_mizuroute_output

  !-----------------------------------------------------------------------
  ! Write mizuRoute streamflow to an existing host-model NetCDF file
  !-----------------------------------------------------------------------
  subroutine write_mizuroute_output(ncid, istart, numtim, info, domain, ierr, message)
  
    integer(i4b),            intent(in)  :: ncid
    integer(i4b),            intent(in)  :: istart
    integer(i4b),            intent(in)  :: numtim
    type(mizuroute_info),    intent(in)  :: info
    type(mizuroute_domain),  intent(in)  :: domain
    integer(i4b),            intent(out) :: ierr
    character(*),            intent(out) :: message
  
    integer(i4b) :: varid_qbasin, varid_Qreach
    integer(i4b) :: iRoute
    integer(i4b), dimension(2) :: start2_basin, count2_basin
    integer(i4b), dimension(3) :: start3_reach, count3_reach
  
    ierr = 0
    message = 'write_mizuroute_output/'
  
    netcdf_block: block
  
      ! get variable ID
      ierr = nf90_inq_varid(ncid, 'q_basin', varid_qbasin); if(ierr/=nf90_noerr) exit netcdf_block
      ierr = nf90_inq_varid(ncid, 'Q_reach', varid_Qreach); if(ierr/=nf90_noerr) exit netcdf_block
  
      ! basin runoff
      start2_basin = (/1,          istart/)
      count2_basin = (/info%n_hru, numtim/)
      
      ierr = nf90_put_var(ncid, varid_qbasin, &
                          domain%river_network%driver%basin_runoff(:,1:numtim), &
                          start=start2_basin, count=count2_basin)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! reach stream for each active routing method
      do iRoute = 1,size(domain%river_network%driver%method)

        start3_reach = (/iRoute,          1, istart/)
        count3_reach = (/     1, info%n_seg, numtim/)
 
        ierr = nf90_put_var(ncid, varid_Qreach, &
                            domain%river_network%driver%method(iRoute)%streamflow(:,1:numtim), &
                            start=start3_reach, count=count3_reach)
        if(ierr/=nf90_noerr) exit netcdf_block

      enddo
  
    end block netcdf_block
  
    if(ierr/=nf90_noerr)then
      message = trim(message)//trim(nf90_strerror(ierr))
      return
    endif
  
  end subroutine write_mizuroute_output


end module mizuroute_output_module
