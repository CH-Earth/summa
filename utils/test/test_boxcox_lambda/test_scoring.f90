program test_boxcox_scoring
 use nr_type,only:rkind,i4b,lgt
 use metrics,only:compute_metric,metric_undefined
 use summa_type,only:config_info
 use simulated_series,only:sim_series_type
 use scorer_under_test,only:score_target
 use read_observations_module,only:input_obs
 use write_evaluation_module,only:write_count,written_objective,written_boxcox_lambda
 use summa_simulation,only:scalarize_objectives
 use globalData,only:iulog
 use, intrinsic :: ieee_arithmetic
 implicit none
 type(config_info) :: config
 type(sim_series_type) :: series(1)
 real(rkind) :: t(3),value,expected,qn,values(3)
 integer :: err,checks,j,referr
 logical(lgt) :: failed
 character(1024) :: message,refmessage
 checks=0
 open(newunit=iulog,status='scratch',action='readwrite')
 allocate(config%calib%targets(3))
 config%calib%start_date='1983-10-01'
 config%calib%end_date='1990-09-30'
 do j=1,3
   config%calib%targets(j)%variable='streamflow'
   config%calib%targets(j)%metric='nse'
   config%calib%targets(j)%obs_transform='boxcox'
   config%calib%targets(j)%weight=0._rkind
   config%calib%targets(j)%obs_path='test/'
   config%calib%targets(j)%obs_file='test.nc'
   config%calib%targets(j)%vname_obs='q_obs'
 enddo
 config%calib%targets(1)%name='kge'
 config%calib%targets(1)%metric='kge'
 config%calib%targets(1)%obs_transform='none'
 config%calib%targets(1)%weight=1._rkind
 config%calib%targets(2)%name='boxcox_half'
 config%calib%targets(2)%boxcox_lambda=0.5_rkind
 config%calib%targets(3)%name='boxcox_log'
 config%calib%targets(3)%boxcox_lambda=0._rkind
 t=[1._rkind,2._rkind,3._rkind]
 input_obs=[0.1_rkind,1._rkind,4._rkind]
 series(1)%values=[0.2_rkind,1.1_rkind,3.7_rkind]
 do j=1,3
   call score(j)
   call check(err==0 .and. .not.failed,'target-specific score is valid')
   call compute_metric(input_obs,series(1)%values,config%calib%targets(j)%metric, &
                       config%calib%targets(j)%obs_transform,expected, &
                       boxcox_lambda=config%calib%targets(j)%boxcox_lambda,err=referr,message=refmessage)
   call check(referr==0 .and. value==expected,'scorer forwards the configured exponent')
   values(j)=value
 enddo
 call check(scalarize_objectives(config,values)==values(1),'DDS still uses only the KGE target')
 call check(write_count==0,'calibration metrics do not require simulation output')
 config%write_timeseries=.true.
 call score(2)
 call check(err==0 .and. write_count==1,'evaluation writer is reached')
 call check(written_boxcox_lambda==0.5_rkind,'evaluation writer receives target exponent')
 call check(written_objective==value,'evaluation writer receives target score')
 call score(3)
 call check(err==0 .and. written_boxcox_lambda==0._rkind,'zero exponent is not replaced by a default')
 input_obs(1)=0._rkind
 call score(3)
 call check(err==0 .and. .not.failed .and. ieee_is_nan(value),'undefined zero-lambda diagnostic is nonfatal')
 config%calib%targets(3)%weight=1._rkind
 call score(3)
 call check(err==metric_undefined .and. failed,'undefined active Box-Cox target fails trial')
 config%calib%targets(3)%weight=0._rkind
 config%calib%targets(3)%boxcox_lambda=ieee_value(0._rkind,ieee_quiet_nan)
 call score(3)
 call check(err==20 .and. .not.failed,'invalid diagnostic lambda remains a configuration error')
 close(iulog)
 write(*,'(A,I0,A)') 'PASS: production target scorer, per-target lambda and diagnostic behavior (',checks,' assertions).'
contains
 subroutine score(j)
 integer,intent(in) :: j
 failed=.false.
 call score_target(config,j,1,series,t,'days since 1980-01-01',value,failed,err,message)
 end subroutine
 subroutine check(ok,label)
 logical,intent(in) :: ok
 character(*),intent(in) :: label
 checks=checks+1
 if(.not.ok)then
   write(*,'(A)') 'FAIL: '//label//' '//trim(message)
   stop 1
 endif
 end subroutine
end program
