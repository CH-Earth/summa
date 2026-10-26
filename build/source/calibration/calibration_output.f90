! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
module calibration_output_module

  USE nr_type, only: i4b,rkind,lgt
  USE netcdf

  USE data_types,       only: target_info
  USE parameter_search, only: parameter_spec

  implicit none
  private

  public :: create_calibration_output
  public :: write_calibration_output
  public :: write_generation_output
  public :: write_pareto_front
  public :: write_failed_trial
  public :: close_calibration_output

contains

  ! **************************************************************************************************
  ! Create calibration output file.
  !
  ! Creates a calibration NetCDF file containing one variable for each parameter represented in the
  ! parameter specification and an objective-function variable. Parameter metadata are stored as
  ! variable attributes, including native bounds, parameter transformation, sampling status, and
  ! ordered-constraint information.
  !
  ! Parameter and objective values are written by global trial index along the fixed sample dimension.
  ! An NSGA-II run also records the population kept after each generation and the final Pareto front.
  ! **************************************************************************************************
  subroutine create_calibration_output(filename,spec,nSamples,nWorkers,node_index,global_case_group, &
                                       case_name,targets,algorithm,population_size,ncid,ierr,message)
    implicit none
    character(*),         intent(in)  :: filename
    type(parameter_spec), intent(in)  :: spec
    integer(i4b),         intent(in)  :: nSamples
    integer(i4b),         intent(in)  :: nWorkers
    integer(i4b),         intent(in)  :: node_index
    integer(i4b),         intent(in)  :: global_case_group
    character(*),         intent(in)  :: case_name
    type(target_info),    intent(in)  :: targets(:)
    character(*),         intent(in)  :: algorithm        ! dds or nsga2
    integer(i4b),         intent(in)  :: population_size  ! NSGA-II population; unused by DDS
    integer(i4b),         intent(out) :: ncid
    integer(i4b),         intent(out) :: ierr
    character(*),         intent(out) :: message
    integer(i4b) :: dim_sample,dim_time,dim_target,dim_name
    integer(i4b) :: varid_sample,varid_objective,varid_param
    integer(i4b) :: varid_worker_rank,varid_target_name
    integer(i4b) :: varid_start_time,varid_end_time
    integer(i4b) :: dim_member,dim_generation,varid
    integer(i4b), dimension(2) :: time_dims
    integer(i4b), dimension(2) :: objective_dims
    integer(i4b), dimension(2) :: name_dims
    integer(i4b) :: iTarget
    integer(i4b) :: nTarget
    character(len=len(targets(1)%name))  :: target_name
    character(len=:), allocatable :: constraint_names
    character(len=64) :: attribute_name
    integer(i4b) :: iParam
    integer(i4b) :: iConstraint
    integer(i4b) :: iOrdered
    integer(i4b) :: sampled
    logical(lgt) :: file_open
    integer(i4b) :: ierr_close

    ierr=0
    message='create_calibration_output/'
    file_open=.false.
    nTarget=size(targets)
    if(nTarget < 1)then
      message=trim(message)//'a calibration output file needs at least one calibration target'
      ierr=20; return
    endif
    netcdf_block: block

      ! create calibration output file
      ierr=nf90_create(trim(filename),NF90_CLOBBER,ncid)
      if(ierr/=nf90_noerr) exit netcdf_block
      file_open=.true.

      ! fixed sample dimension
      ierr=nf90_def_dim(ncid,'sample',nSamples,dim_sample)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! sample coordinate
      ierr=nf90_def_var(ncid,'sample',NF90_INT,(/dim_sample/),varid_sample)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_sample,'long_name','parameter trial')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_sample,'units','-')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! time-component dimension
      ierr=nf90_def_dim(ncid,'time_component',8,dim_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      time_dims=(/dim_time,dim_sample/)

      ! calibration-target dimension: one objective value per target, per trial
      ierr=nf90_def_dim(ncid,'target',nTarget,dim_target)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_def_dim(ncid,'name_length',len(target_name),dim_name)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -----------------------------------------------------------------------------------------------
      ! Parameter variables and metadata
      ! -----------------------------------------------------------------------------------------------
      do iParam=1,size(spec%params)
        ierr=nf90_def_var(ncid,trim(spec%params(iParam)%name),NF90_DOUBLE, (/dim_sample/),varid_param)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid_param,'long_name', trim(spec%params(iParam)%long_name))
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid_param,'units', trim(spec%params(iParam)%units))
        if(ierr/=nf90_noerr) exit netcdf_block

        ! native physical parameter bounds
        ierr=nf90_put_att(ncid,varid_param,'lower_bound',spec%params(iParam)%lower)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid_param,'upper_bound',spec%params(iParam)%upper)
        if(ierr/=nf90_noerr) exit netcdf_block

        ! parameter transformation used during search
        ierr=nf90_put_att(ncid,varid_param,'transformation', trim(spec%params(iParam)%transformation))
        if(ierr/=nf90_noerr) exit netcdf_block

        ! parameter sampling status
        if(spec%params(iParam)%sampled)then
          sampled=1
        else
          sampled=0
        endif
        ierr=nf90_put_att(ncid,varid_param,'sampled',sampled)
        if(ierr/=nf90_noerr) exit netcdf_block

        ! scalar trial value supplied by the model-specific parameter specification
        ierr=nf90_put_att(ncid,varid_param,'trial_value',spec%params(iParam)%trial_value)
        if(ierr/=nf90_noerr) exit netcdf_block

      enddo

      ! -----------------------------------------------------------------------------------------------
      ! Rank/timing information
      ! -----------------------------------------------------------------------------------------------

      ! worker rank
      ierr=nf90_def_var(ncid,'worker_rank',NF90_INT,(/dim_sample/),varid_worker_rank)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_worker_rank,'long_name', 'MPI worker rank responsible for parameter trial')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_worker_rank,'units','-')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! parameter-trial start time
      ierr=nf90_def_var(ncid,'dispatch_time',NF90_INT,time_dims,varid_start_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_start_time,'long_name','parameter trial dispatch time')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_start_time,'components', &
                        'year month day UTC_offset hour minute second millisecond')
      if(ierr/=nf90_noerr) exit netcdf_block
      
      ! parameter-trial end time
      ierr=nf90_def_var(ncid,'completion_time',NF90_INT,time_dims,varid_end_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_end_time,'long_name','parameter trial completion time')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_end_time,'components', &
                        'year month day UTC_offset hour minute second millisecond')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -----------------------------------------------------------------------------------------------
      ! Objective function: one value per calibration target, for every trial
      ! -----------------------------------------------------------------------------------------------
      objective_dims=(/dim_target,dim_sample/)
      ierr=nf90_def_var(ncid,'objective',NF90_DOUBLE,objective_dims,varid_objective)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'long_name','calibration objective function')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'units','-')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! the metric and transformation behind each target's value, in target order
      ierr=nf90_put_att(ncid,varid_objective,'metric',trim(joined_field(targets,'metric')))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'obs_transform',trim(joined_field(targets,'obs_transform')))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'target',trim(joined_field(targets,'name')))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'variable',trim(joined_field(targets,'variable')))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'sense',trim(joined_field(targets,'sense')))
      if(ierr/=nf90_noerr) exit netcdf_block

      ! a failed trial has no objectives and ranks last on every target
      ierr=nf90_def_var(ncid,'trial_failed',NF90_INT,(/dim_sample/),varid)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid,'long_name','1 if the model refused or failed on the trial; its objectives are fill values')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! target names, so a reader can label the objectives without parsing an attribute
      name_dims=(/dim_name,dim_target/)
      ierr=nf90_def_var(ncid,'target_name',NF90_CHAR,name_dims,varid_target_name)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_target_name,'long_name','name of each calibration target')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -----------------------------------------------------------------------------------------------
      ! NSGA-II: the population kept after each generation, and the non-dominated trials
      ! -----------------------------------------------------------------------------------------------
      if(trim(algorithm) == 'nsga2')then
        ierr=nf90_def_dim(ncid,'member',population_size,dim_member)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_def_dim(ncid,'generation',nSamples/population_size,dim_generation)
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(ncid,'birth_generation',NF90_INT,(/dim_sample/),varid)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid,'long_name','generation that proposed the trial, 1 the random initial population')
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(ncid,'population',NF90_INT,(/dim_member,dim_generation/),varid)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid,'long_name','sample index of each member kept after the generation')
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(ncid,'population_rank',NF90_INT,(/dim_member,dim_generation/),varid)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid,'long_name','non-domination front of each member, 1 is non-dominated')
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(ncid,'population_crowding',NF90_DOUBLE,(/dim_member,dim_generation/),varid)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid,'long_name','crowding distance of each member in its front, infinite at its extremes')
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(ncid,'pareto_front',NF90_INT,(/dim_sample/),varid)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid,'long_name','1 if no other trial is at least as good on every target and better on one')
        if(ierr/=nf90_noerr) exit netcdf_block
      endif

      ! -----------------------------------------------------------------------------------------------
      ! Global metadata
      ! -----------------------------------------------------------------------------------------------
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'title','SUMMA calibration parameter trials')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'case_name',trim(case_name))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'mpi_workers',nWorkers)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'node_index',node_index)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'global_case_group',global_case_group)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'parameter_space','physical')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'algorithm',trim(algorithm))
      if(ierr/=nf90_noerr) exit netcdf_block

      ! ordered parameter constraints
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'n_ordered_constraints',size(spec%ordered))
      if(ierr/=nf90_noerr) exit netcdf_block
      
      do iConstraint=1,size(spec%ordered)
        ! construct ordered parameter list
        constraint_names=''
        do iOrdered=1,size(spec%ordered(iConstraint)%param_index)
          iParam=spec%ordered(iConstraint)%param_index(iOrdered)
          if(iOrdered > 1) constraint_names=constraint_names//','
          constraint_names=constraint_names//trim(spec%params(iParam)%name)
        enddo
      
        ! write parameter names
        write(attribute_name,'("ordered_constraint_",I0)') iConstraint
        ierr=nf90_put_att(ncid,NF90_GLOBAL,trim(attribute_name),constraint_names)
        if(ierr/=nf90_noerr) exit netcdf_block
      
        ! write gap fraction
        ierr=nf90_put_att(ncid,NF90_GLOBAL,trim(attribute_name)//'_gap_fraction', &
                          spec%ordered(iConstraint)%gap_fraction)
        if(ierr/=nf90_noerr) exit netcdf_block
      enddo

      ! leave define mode
      ierr=nf90_enddef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! target names, written once now that the file holds them
      do iTarget=1,nTarget
        target_name=targets(iTarget)%name
        ierr=nf90_put_var(ncid,varid_target_name,target_name, &
                          start=(/1,iTarget/),count=(/len(target_name),1/))
        if(ierr/=nf90_noerr) exit netcdf_block
      enddo

    end block netcdf_block

    ! process NetCDF errors
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      if(file_open) ierr_close=nf90_close(ncid)
      return
    endif
    ierr=0

  end subroutine create_calibration_output

  ! **************************************************************************************************
  ! Join one field of every calibration target into a comma-separated list, in target order.
  !
  ! Written as a variable attribute so a reader can see what each objective is without opening the
  ! configuration that produced it.
  ! **************************************************************************************************
  function joined_field(targets,field) result(joined)
    USE metrics, only: metric_is_maximized
    implicit none
    type(target_info), intent(in) :: targets(:)
    character(*),      intent(in) :: field
    character(len=:), allocatable :: joined
    integer(i4b) :: iTarget

    joined=''
    do iTarget=1,size(targets)
      if(iTarget > 1) joined=joined//', '
      select case(trim(field))
        case ('name');          joined=joined//trim(targets(iTarget)%name)
        case ('variable');      joined=joined//trim(targets(iTarget)%variable)
        case ('metric');        joined=joined//trim(targets(iTarget)%metric)
        case ('obs_transform'); joined=joined//trim(targets(iTarget)%obs_transform)
        case ('sense')
          if(metric_is_maximized(targets(iTarget)%metric))then
            joined=joined//'maximize'
          else
            joined=joined//'minimize'
          endif
        case default;           joined=joined//'unknown'
      end select
    enddo

  end function joined_field

  ! **************************************************************************************************
  ! Write one calibration trial.
  !
  ! Writes the complete physical-space parameter override vector and the objective value of every
  ! calibration target for one parameter trial. Parameter names must correspond to variables defined
  ! when the calibration output file was created.
  ! **************************************************************************************************
  subroutine write_calibration_output(ncid,isample,worker_rank,param_names,param_values, &
                                      objective,failed,start_time,end_time,ierr,message)
    implicit none
    integer(i4b), intent(in) :: ncid
    integer(i4b), intent(in) :: isample
    integer(i4b), intent(in) :: worker_rank
    character(*), intent(in) :: param_names(:)
    real(rkind),  intent(in) :: param_values(:)
    real(rkind),  intent(in) :: objective(:)
    logical(lgt), intent(in) :: failed           ! the trial failed, so its objectives are written as fill values
    integer(i4b), intent(in) :: start_time(8)
    integer(i4b), intent(in) :: end_time(8)
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message
    integer(i4b) :: varid_sample, varid_worker_rank
    integer(i4b) :: varid_start_time,varid_end_time
    integer(i4b) :: varid_param,varid_objective,varid_failed
    integer(i4b) :: iParam
    integer(i4b), dimension(1) :: start1,count1
    integer(i4b), dimension(2) :: start2,count2

    ierr=0
    message='write_calibration_output/'
    if(size(param_names) /= size(param_values))then
      message=trim(message)//'parameter name and value vectors have different sizes'
      ierr=20; return
    endif
    start1=(/isample/)
    count1=(/1/)
    start2=(/1,isample/)
    count2=(/8,1/)
    netcdf_block: block

      ! sample coordinate
      ierr=nf90_inq_varid(ncid,'sample',varid_sample)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_sample,(/isample/),start=start1,count=count1)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! parameter values
      do iParam=1,size(param_names)
        ierr=nf90_inq_varid(ncid,trim(param_names(iParam)),varid_param)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_var(ncid,varid_param,(/param_values(iParam)/), start=start1,count=count1)
        if(ierr/=nf90_noerr) exit netcdf_block
      enddo

      ! worker rank
      ierr=nf90_inq_varid(ncid,'worker_rank',varid_worker_rank)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_worker_rank,(/worker_rank/),start=start1,count=count1)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! parameter-trial start times
      ierr=nf90_inq_varid(ncid,'dispatch_time',varid_start_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_start_time,start_time, start=start2,count=count2)
      if(ierr/=nf90_noerr) exit netcdf_block
     
      ! parameter-trial end times
      ierr=nf90_inq_varid(ncid,'completion_time',varid_end_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_end_time,end_time, start=start2,count=count2)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! objective function: every target's value for this trial
      ierr=nf90_inq_varid(ncid,'objective',varid_objective)
      if(ierr/=nf90_noerr) exit netcdf_block
      if(failed)then
        ierr=nf90_put_var(ncid,varid_objective,spread(NF90_FILL_DOUBLE,1,size(objective)), &
                          start=(/1,isample/),count=(/size(objective),1/))
      else
        ierr=nf90_put_var(ncid,varid_objective,objective, &
                          start=(/1,isample/),count=(/size(objective),1/))
      endif
      if(ierr/=nf90_noerr) exit netcdf_block

      ! failure flag
      ierr=nf90_inq_varid(ncid,'trial_failed',varid_failed)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_failed,(/merge(1,0,failed)/),start=start1,count=count1)
      if(ierr/=nf90_noerr) exit netcdf_block

    end block netcdf_block
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      return
    endif
    ierr=0

  end subroutine write_calibration_output

  ! **************************************************************************************************
  ! Write one NSGA-II generation: the generation that proposed trials first..last, and the members
  ! kept after its selection with their fronts and crowding distances.
  ! **************************************************************************************************
  subroutine write_generation_output(ncid,iGen,first,last,member,rank,distance,ierr,message)
    implicit none
    integer(i4b), intent(in)  :: ncid
    integer(i4b), intent(in)  :: iGen           ! generation
    integer(i4b), intent(in)  :: first,last     ! trials it proposed
    integer(i4b), intent(in)  :: member(:)      ! sample index of each member kept
    integer(i4b), intent(in)  :: rank(:)        ! their fronts
    real(rkind),  intent(in)  :: distance(:)    ! their crowding distances
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message
    integer(i4b) :: varid

    ierr=0
    message='write_generation_output/'
    netcdf_block: block
      ierr=nf90_inq_varid(ncid,'birth_generation',varid)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid,spread(iGen,1,last-first+1),start=(/first/),count=(/last-first+1/))
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr=nf90_inq_varid(ncid,'population',varid)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid,member,start=(/1,iGen/),count=(/size(member),1/))
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr=nf90_inq_varid(ncid,'population_rank',varid)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid,rank,start=(/1,iGen/),count=(/size(rank),1/))
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr=nf90_inq_varid(ncid,'population_crowding',varid)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid,distance,start=(/1,iGen/),count=(/size(distance),1/))
      if(ierr/=nf90_noerr) exit netcdf_block
    end block netcdf_block
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      return
    endif
    ierr=0

  end subroutine write_generation_output

  ! **************************************************************************************************
  ! Write which trials are non-dominated among every trial evaluated.
  ! **************************************************************************************************
  subroutine write_pareto_front(ncid,front,ierr,message)
    implicit none
    integer(i4b), intent(in)  :: ncid
    logical(lgt), intent(in)  :: front(:)       ! one flag per trial
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message
    integer(i4b) :: varid

    ierr=0
    message='write_pareto_front/'
    ierr=nf90_inq_varid(ncid,'pareto_front',varid)
    if(ierr==nf90_noerr) ierr=nf90_put_var(ncid,varid,merge(1,0,front))
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      return
    endif
    ierr=0

  end subroutine write_pareto_front

  ! **************************************************************************************************
  ! Save a failed trial as a standalone repro case, in <output_path>/failed_trials/<algorithm>_g<generation>_s<sample>/
  ! (DDS has no generations, so dds_s<sample>/).
  !
  ! The directory holds trial.json (parameters, error, worker, source revision), the configuration,
  ! the spun-up initial state every trial starts from,
  ! and rerun.sh, which reruns the trial with the serial driver built beside this one. Returns the
  ! directory in dir.
  ! **************************************************************************************************
  subroutine write_failed_trial(config,sample_id,worker_rank,param_names,param_values,error_message,dir)
    USE summa_type,       only: config_info
    USE summaFileManager, only: OUTPUT_PATH,STATE_PATH,SETTINGS_PATH,MODEL_INITCOND
    implicit none
    type(config_info), intent(in)  :: config
    integer(i4b),      intent(in)  :: sample_id
    integer(i4b),      intent(in)  :: worker_rank
    character(*),      intent(in)  :: param_names(:)
    real(rkind),       intent(in)  :: param_values(:)
    character(*),      intent(in)  :: error_message
    character(len=:), allocatable, intent(out) :: dir
    character(len=:), allocatable :: exe,serial_exe,init_state,source_rev,launch_dir
    character(len=32)   :: val
    character(len=4096) :: line
    integer(i4b) :: unit,iParam,ix,ios
    logical(lgt) :: rerunnable

    ! one directory per trial
    write(val,'(i6.6)') sample_id
    if(trim(config%calib%algorithm) == 'nsga2')then
      write(line,'(a,i0,a)') 'nsga2_g',(sample_id-1)/config%calib%nsga2%population_size+1,'_s'//trim(val)
    else
      line=trim(config%calib%algorithm)//'_s'//trim(val)
    endif
    dir=trim(OUTPUT_PATH)//'failed_trials/'//trim(line)//'/'
    call execute_command_line('mkdir -p "'//dir//'"')

    ! the executable that ran the trial, and its serial driver (_opt dropped)
    call getcwd(line)
    launch_dir=trim(line)
    call get_command_argument(0,line)
    exe=trim(line)
    if(exe(1:1) /= '/') exe=launch_dir//'/'//exe
    ix=index(exe,'_opt',back=.true.)
    serial_exe=exe
    if(ix > index(exe,'/',back=.true.)) serial_exe=exe(1:ix-1)//exe(ix+4:)

    ! the source revision of the checkout the executable lives in
    call execute_command_line('git -C "'//exe(1:index(exe,'/',back=.true.))//'" describe --always --dirty > "'// &
                              dir//'git_describe.txt" 2>/dev/null')
    source_rev='unknown'
    open(newunit=unit,file=dir//'git_describe.txt',status='old',action='read',iostat=ios)
    if(ios==0)then
      read(unit,'(a)',iostat=ios) line
      if(ios==0 .and. len_trim(line) > 0) source_rev=trim(line)
      close(unit,status='delete')
    endif

    ! the spun-up state the trial started from, as SUMMA resolves it
    if(STATE_PATH == '')then
      init_state=trim(SETTINGS_PATH)//trim(MODEL_INITCOND)
    else
      init_state=trim(STATE_PATH)//trim(MODEL_INITCOND)
    endif
    call execute_command_line('cp "'//init_state//'" "'//dir//'initial_state.nc"')

    ! a manifest case is built from a template, so only a single configuration file can be rerun as it stands
    rerunnable=allocated(config%config_file) .and. .not.allocated(config%manifest_file)
    if(allocated(config%config_file)) call execute_command_line('cp "'//config%config_file//'" "'//dir//'config.toml"')

    ! trial.json
    open(newunit=unit,file=dir//'trial.json',status='replace',action='write')
    write(unit,'(a)') '{'
    write(unit,'(a,i0,a)') '  "sample": ',sample_id,','
    write(unit,'(a)')      '  "algorithm": "'//trim(config%calib%algorithm)//'",'
    if(trim(config%calib%algorithm) == 'nsga2') &
      write(unit,'(a,i0,a)') '  "generation": ',(sample_id-1)/config%calib%nsga2%population_size+1,','
    write(unit,'(a,i0,a)') '  "worker_rank": ',worker_rank,','
    if(allocated(config%case_name))   write(unit,'(a)') '  "case": "'//json_escape(config%case_name)//'",'
    if(allocated(config%config_file)) write(unit,'(a)') '  "config_file": "'//json_escape(config%config_file)//'",'
    write(unit,'(a)') '  "initial_state": "'//json_escape(init_state)//'",'
    write(unit,'(a)') '  "executable": "'//json_escape(exe)//'",'
    write(unit,'(a)') '  "source_revision": "'//json_escape(source_rev)//'",'
    write(unit,'(a)') '  "error": "'//json_escape(trim(error_message))//'",'
    write(unit,'(a)') '  "parameters": {'
    do iParam=1,size(param_names)
      write(val,'(es24.16)') param_values(iParam)
      if(iParam < size(param_names)) val=trim(adjustl(val))//','
      write(unit,'(a)') '    "'//trim(param_names(iParam))//'": '//trim(adjustl(val))
    enddo
    write(unit,'(a)') '  }'
    write(unit,'(a)') '}'
    close(unit)

    ! rerun.sh
    open(newunit=unit,file=dir//'rerun.sh',status='replace',action='write')
    write(unit,'(a)') '#!/bin/bash'
    write(unit,'(a)') '# Reruns this failed calibration trial with the serial executable; trial.json says how it failed.'
    write(unit,'(a)') '# Set SUMMA_EXE to use another build. Output goes to output/ and rerun.log beside this script.'
    write(unit,'(a)') 'set -euo pipefail'
    if(.not.rerunnable)then
      write(unit,'(a)') 'echo "this trial came from a run manifest or file manager, not a single TOML configuration;"'
      write(unit,'(a)') 'echo "rerun it by hand from trial.json"; exit 1'
    else
      write(unit,'(a)') 'here="$(cd "$(dirname "$0")" && pwd)"'
      write(unit,'(a)') 'SUMMA_EXE="${SUMMA_EXE:-'//serial_exe//'}"'
      write(unit,'(a)') 'mkdir -p "${here}/output"'
      write(unit,'(a)') 'sed -e "s|^\([[:space:]]*state_path[[:space:]]*=\).*|\1 \"${here}/\"|" \'
      write(unit,'(a)') '    -e "s|^\([[:space:]]*init_condition[[:space:]]*=\).*|\1 \"initial_state.nc\"|" \'
      write(unit,'(a)') '    -e "s|^\([[:space:]]*output_path[[:space:]]*=\).*|\1 \"${here}/output/\"|" \'
      write(unit,'(a)') '    -e "s|^\([[:space:]]*work_path[[:space:]]*=\).*|\1 \"${here}/output\"|" \'
      write(unit,'(a)') '    "${here}/config.toml" > "${here}/rerun.toml"'
      write(unit,'(a)') '# relative paths in the configuration resolve from where the calibration was launched'
      write(unit,'(a)') 'cd "'//launch_dir//'"'
      write(unit,'(a)') '"${SUMMA_EXE}" -c "${here}/rerun.toml" \'
      do iParam=1,size(param_names)
        write(val,'(es24.16)') param_values(iParam)
        write(unit,'(a)') '  --param '//trim(param_names(iParam))//' '//trim(adjustl(val))//' \'
      enddo
      write(unit,'(a)') '  2>&1 | tee "${here}/rerun.log"'
    endif
    close(unit)
    call execute_command_line('chmod +x "'//dir//'rerun.sh"')

  contains

    ! escape a string for a JSON value
    function json_escape(s) result(e)
      character(*), intent(in)      :: s
      character(len=:), allocatable :: e
      integer(i4b) :: i
      e=''
      do i=1,len(s)
        select case(s(i:i))
          case('"','\'); e=e//'\'//s(i:i)
          case default;  e=e//s(i:i)
        end select
      enddo
    end function json_escape

  end subroutine write_failed_trial

  ! **************************************************************************************************
  ! Close calibration output file.
  ! **************************************************************************************************
  subroutine close_calibration_output(ncid,ierr,message)
    implicit none
    integer(i4b), intent(in)  :: ncid
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    ierr=0
    message='close_calibration_output/'
    ierr=nf90_close(ncid)
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      return
    endif
    ierr=0

  end subroutine close_calibration_output

end module calibration_output_module
