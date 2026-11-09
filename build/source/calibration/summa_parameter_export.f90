
! Export calibrated physical SUMMA parameters to a trial-parameter NetCDF file.
module summa_parameter_export

  USE nr_type, only: i4b,rkind
  USE netcdf
  USE data_types, only: gru_hru_dom_doubleVec,gru_double
  USE globalData, only: gru_struc,mpar_meta,bpar_meta,maxSoilLayers,integerMissing
  USE get_ixname_module, only: get_ixVarType,get_ixParam,get_ixBpar
  USE parameter_search, only: parameter_spec
  USE summa_parameter_fields, only: summa_parameter_map

  implicit none
  private

  type :: export_parameter
    character(len=64) :: name = ''
    integer(i4b)      :: ixM = 0   ! Model parameter index
    integer(i4b)      :: ixB = 0   ! Basin parameter index
  end type export_parameter

  public :: write_physical_trial_parameters

contains

  ! **************************************************************************************************
  ! Write the effective physical values of sampled calibration parameters.
  !
  ! The calibration adapter has already applied the selected parameter sample
  ! to mpar and bpar. Parameter names and metadata come from the specification.
  !
  ! HRUs with several domains get a dom dimension, as read_trial_parameters reads it.
  ! **************************************************************************************************
  subroutine write_physical_trial_parameters(filename,trial_param_names, &
                                             spec,param_map,mpar,bpar,sample_id,err,message)

    implicit none

    ! Dummy arguments
    character(*),                intent(in)  :: filename             ! parameter trial filename
    character(len=64),           intent(in)  :: trial_param_names(:) ! names of trial model parameters
    type(parameter_spec),        intent(in)  :: spec                 ! Complete calibration parameter specification
    type(summa_parameter_map),   intent(in)  :: param_map(:)         ! Mapping from calibration parameters to SUMMA parameter indices
    type(gru_hru_dom_doubleVec), intent(in)  :: mpar                 ! Physical HRU parameter fields for the selected sample
    type(gru_double),            intent(in)  :: bpar                 ! Physical basin parameter fields for the selected sample
    integer(i4b),                intent(in)  :: sample_id            ! Calibration sample index (zero for default parameters)
    integer(i4b),                intent(out) :: err                  ! Error code
    character(*),                intent(out) :: message              ! Error message

    ! NetCDF identifiers
    integer(i4b) :: ncid,dim_hru,dim_gru,dim_depth,dim_dom
    integer(i4b) :: var_hru,var_gru
    integer(i4b), allocatable :: var_param(:)

    ! Spatial indices
    integer(i4b) :: nGRU,nHRU,nDom
    integer(i4b) :: g,h,i,k,n,ix,d,iAt
    integer(i4b) :: soil_type

    ! Parameter values and identifiers
    integer(i4b), allocatable :: hru_id(:),gru_id(:)
    real(rkind), allocatable :: hru_value(:,:),gru_value(:),soil_value(:,:,:)
    character(len=64) :: base_name

    ! List of trial model parameters
    type(export_parameter), allocatable :: export_list(:)
    integer(i4b) :: nExport,ixM,ixB
    
    ! Error handling
    integer(i4b) :: ierr,ierr_close
    logical :: file_open,has_soil,is_soil

    err=0
    message='write_physical_trial_parameters/'
    file_open=.false.

    ! ---------------------------------------------------------------------------------------
    ! Determine the number of GRUs and HRUs.
    ! ---------------------------------------------------------------------------------------
    nGRU=size(bpar%gru)
    nHRU=0
    nDom=1

    do g=1,nGRU
      nHRU=nHRU+size(mpar%gru(g)%hru)
      do h=1,size(mpar%gru(g)%hru)
        nDom=max(nDom,size(mpar%gru(g)%hru(h)%dom))
      enddo
    enddo

    ! ---------------------------------------------------------------------------------------
    ! Combine the original trial parameters with the sampled calibration parameters.
    ! ---------------------------------------------------------------------------------------

    allocate(export_list(size(trial_param_names)+size(spec%params)))
    nExport=0
    
    ! Parameters present in the original trial parameter file.
    do i=1,size(trial_param_names)
    
      ixM=get_ixParam(trim(trial_param_names(i)))
      ixB=get_ixBpar(trim(trial_param_names(i)))
    
      if(ixM==integerMissing .and. ixB==integerMissing)then
        message=trim(message)//'unknown trial parameter: '//trim(trial_param_names(i))
        err=20; return
      endif
    
      nExport=nExport+1
      export_list(nExport)%name=trim(trial_param_names(i))
    
      if(ixM/=integerMissing) export_list(nExport)%ixM=ixM
      if(ixB/=integerMissing) export_list(nExport)%ixB=ixB
    
    enddo
    
    ! Add sampled parameters that are not already in the original file.
    do i=1,size(spec%params)
    
      if(.not.spec%params(i)%sampled) cycle
    
      ! name@type is written under its parameter name: the file holds every domain's value
      base_name=spec%params(i)%name
      iAt=index(base_name,'@')
      if(iAt>0) base_name=base_name(1:iAt-1)

      if(nExport>0)then
        if(any(export_list(1:nExport)%name==trim(base_name))) cycle
      endif
    
      nExport=nExport+1
      export_list(nExport)%name=trim(base_name)
      export_list(nExport)%ixM=param_map(i)%ixM
      export_list(nExport)%ixB=param_map(i)%ixB
    
    enddo

    ! Determine whether any exported parameter requires a soil-depth dimension.
    soil_type=get_ixVarType('parSoil')
    has_soil=.false.

    do i=1,nExport
      if(export_list(i)%ixM>0)then
        if(mpar_meta(export_list(i)%ixM)%varType==soil_type) has_soil=.true.
      endif
    enddo

    ! ---------------------------------------------------------------------------------------
    ! Allocate arrays and construct HRU/GRU identifiers.
    ! ---------------------------------------------------------------------------------------
    allocate(var_param(nExport))
    allocate(hru_id(nHRU),gru_id(nGRU))
    allocate(hru_value(nHRU,nDom),gru_value(nGRU))
    if(has_soil) allocate(soil_value(nHRU,maxSoilLayers,nDom))

    k=0
    do g=1,nGRU
      gru_id(g)=int(gru_struc(g)%gru_id,i4b)

      do h=1,size(mpar%gru(g)%hru)
        k=k+1
        hru_id(k)=int(gru_struc(g)%hruInfo(h)%hru_id,i4b)
      enddo
    enddo

    ! ---------------------------------------------------------------------------------------
    ! Create the NetCDF file.
    ! ---------------------------------------------------------------------------------------
    netcdf_block: block

      ierr=nf90_create(trim(filename),nf90_clobber,ncid)
      if(ierr/=nf90_noerr) exit netcdf_block
      file_open=.true.

      ! Define spatial dimensions.
      ierr=nf90_def_dim(ncid,'hru',nHRU,dim_hru)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr=nf90_def_dim(ncid,'gru',nGRU,dim_gru)
      if(ierr/=nf90_noerr) exit netcdf_block

      if(has_soil)then
        ierr=nf90_def_dim(ncid,'depth',maxSoilLayers,dim_depth)
        if(ierr/=nf90_noerr) exit netcdf_block
      endif

      if(nDom>1)then
        ierr=nf90_def_dim(ncid,'dom',nDom,dim_dom)
        if(ierr/=nf90_noerr) exit netcdf_block
      endif

      ! Define HRU and GRU identifiers.
      ierr=nf90_def_var(ncid,'hruId',nf90_int,[dim_hru],var_hru)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr=nf90_def_var(ncid,'gruId',nf90_int,[dim_gru],var_gru)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -------------------------------------------------------------------------------------
      ! Define sampled parameter variables.
      ! -------------------------------------------------------------------------------------
      do i=1,nExport

        if(export_list(i)%ixB>0)then

          ! Basin parameter: one value per GRU.
          ierr=nf90_def_var(ncid,trim(export_list(i)%name), &
                            nf90_double,[dim_gru],var_param(i))
          if(ierr/=nf90_noerr) exit netcdf_block

          ! Parameter metadata.
          ierr=nf90_put_att(ncid,var_param(i),'long_name',trim(bpar_meta(export_list(i)%ixB)%varDesc))
          if(ierr/=nf90_noerr) exit netcdf_block
          
          ierr=nf90_put_att(ncid,var_param(i),'units',trim(bpar_meta(export_list(i)%ixB)%varUnit))
          if(ierr/=nf90_noerr) exit netcdf_block

        else

          ! Model parameter: one value per HRU or soil layer.
          ix=export_list(i)%ixM
          is_soil=mpar_meta(ix)%varType==soil_type

          if(is_soil .and. nDom>1)then
            ierr=nf90_def_var(ncid,trim(export_list(i)%name), &
                              nf90_double,[dim_hru,dim_depth,dim_dom],var_param(i))
          elseif(is_soil)then
            ierr=nf90_def_var(ncid,trim(export_list(i)%name), &
                              nf90_double,[dim_hru,dim_depth],var_param(i))
          elseif(nDom>1)then
            ierr=nf90_def_var(ncid,trim(export_list(i)%name), &
                              nf90_double,[dim_hru,dim_dom],var_param(i))
          else
            ierr=nf90_def_var(ncid,trim(export_list(i)%name), &
                              nf90_double,[dim_hru],var_param(i))
          endif
          if(ierr/=nf90_noerr) exit netcdf_block

          ! Parameter metadata.
          ierr=nf90_put_att(ncid,var_param(i),'long_name',trim(mpar_meta(export_list(i)%ixM)%varDesc))
          if(ierr/=nf90_noerr) exit netcdf_block
          
          ierr=nf90_put_att(ncid,var_param(i),'units',trim(mpar_meta(export_list(i)%ixM)%varUnit))
          if(ierr/=nf90_noerr) exit netcdf_block

        endif

      enddo

      ! Sample metadata.
      ierr=nf90_put_att(ncid,nf90_global,'calibration_sample',sample_id)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr=nf90_put_att(ncid,nf90_global,'parameter_space','physical')
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr=nf90_enddef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -------------------------------------------------------------------------------------
      ! Write HRU and GRU identifiers.
      ! -------------------------------------------------------------------------------------
      ierr=nf90_put_var(ncid,var_hru,hru_id)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr=nf90_put_var(ncid,var_gru,gru_id)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -------------------------------------------------------------------------------------
      ! Write physical parameter values.
      ! -------------------------------------------------------------------------------------
      do i=1,nExport

        if(export_list(i)%ixB>0)then

          ! Basin parameter values.
          ix=export_list(i)%ixB

          do g=1,nGRU
            gru_value(g)=bpar%gru(g)%var(ix)
          enddo

          ierr=nf90_put_var(ncid,var_param(i),gru_value)

        else

          ! HRU parameter values.
          ix=export_list(i)%ixM
          is_soil=mpar_meta(ix)%varType==soil_type
          k=0

          ! Domains an HRU does not have are never read, and left as fill values.
          if(is_soil)then
            soil_value=nf90_fill_double
          else
            hru_value=nf90_fill_double
          endif

          do g=1,nGRU
            do h=1,size(mpar%gru(g)%hru)

              k=k+1

              do d=1,size(mpar%gru(g)%hru(h)%dom)
                if(is_soil)then
                  ! Soil parameters retain their vertical distribution.
                  n=size(mpar%gru(g)%hru(h)%dom(d)%var(ix)%dat)
                  soil_value(k,1:n,d)=mpar%gru(g)%hru(h)%dom(d)%var(ix)%dat
                else
                  ! Scalar parameters have one physical value per HRU domain.
                  hru_value(k,d)=mpar%gru(g)%hru(h)%dom(d)%var(ix)%dat(1)
                endif
              enddo

            enddo
          enddo

          if(is_soil .and. nDom>1)then
            ierr=nf90_put_var(ncid,var_param(i),soil_value)
          elseif(is_soil)then
            ierr=nf90_put_var(ncid,var_param(i),soil_value(:,:,1))
          elseif(nDom>1)then
            ierr=nf90_put_var(ncid,var_param(i),hru_value)
          else
            ierr=nf90_put_var(ncid,var_param(i),hru_value(:,1))
          endif

        endif

        if(ierr/=nf90_noerr) exit netcdf_block

      enddo

    end block netcdf_block

    ! ---------------------------------------------------------------------------------------
    ! Handle NetCDF errors and close the file.
    ! ---------------------------------------------------------------------------------------
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      err=20
    endif

    if(file_open)then
      ierr_close=nf90_close(ncid)

      if(ierr_close/=nf90_noerr .and. err==0)then
        message=trim(message)//trim(nf90_strerror(ierr_close))
        err=20
      endif
    endif

  end subroutine write_physical_trial_parameters

end module summa_parameter_export

