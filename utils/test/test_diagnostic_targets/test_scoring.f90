program test_scoring
 use nr_type,only:rkind,i4b,lgt
 use metrics,only:metric_undefined
 use summa_type,only:config_info
 use simulated_series,only:sim_series_type
 use scorer_under_test,only:score_target
 use read_observations_module,only:input_obs,input_error
 use timeseries_alignment,only:alignment_error
 use write_evaluation_module,only:write_count,written_objective,write_error
 use globalData,only:iulog
 use, intrinsic :: ieee_arithmetic
 implicit none
 type(config_info) :: config
 type(sim_series_type) :: series(1)
 real(rkind) :: t(3),value,qn
 integer :: err,checks,ioerr
 logical(lgt) :: failed
 character(1024) :: message,line
 checks=0
 open(newunit=iulog,status='scratch',action='readwrite')
 allocate(config%calib%targets(1))
 config%calib%start_date='1983-10-01'
 config%calib%end_date='1990-09-30'
 config%calib%targets(1)%name='nse_log'
 config%calib%targets(1)%variable='streamflow'
 config%calib%targets(1)%metric='nse'
 config%calib%targets(1)%obs_transform='log'
 config%calib%targets(1)%weight=0._rkind
 config%calib%targets(1)%obs_path='test/'
 config%calib%targets(1)%obs_file='test.nc'
 config%calib%targets(1)%vname_obs='q_obs'
 t=[1._rkind,2._rkind,3._rkind]
 input_obs=[0._rkind,1._rkind,2._rkind]
 series(1)%values=input_obs
 qn=ieee_value(0._rkind,ieee_quiet_nan)
 call score()
 call check(err==0 .and. .not.failed,'valid diagnostic does not fail trial')
 call check(abs(value-1._rkind)<1.e-14_rkind,'valid diagnostic score')
 call check(write_count==0,'calibration metrics do not need full time-series output')
 input_obs=0._rkind; series(1)%values=0._rkind
 config%write_timeseries=.true.
 call score()
 call check(err==0 .and. .not.failed,'undefined diagnostic is nonfatal')
 call check(ieee_is_nan(value),'undefined diagnostic returns NaN')
 call check(write_count==1 .and. ieee_is_nan(written_objective),'undefined diagnostic reaches evaluation writer')
 rewind(iulog)
 read(iulog,'(A)',iostat=ioerr) line
 call check(ioerr==0 .and. index(line,'undefined diagnostic:')>0,'undefined diagnostic reason logged')
 config%calib%targets(1)%weight=1._rkind
 call score()
 call check(err==metric_undefined .and. failed,'undefined active target fails trial')
 config%calib%targets(1)%weight=0._rkind
 config%calib%targets(1)%metric='typo'
 call score()
 call check(err==20 .and. .not.failed,'diagnostic metric typo remains fatal')
 config%calib%targets(1)%metric='nse'
 config%calib%targets(1)%obs_transform='typo'
 call score()
 call check(err==20 .and. .not.failed,'diagnostic transform typo remains fatal')
 config%calib%targets(1)%obs_transform='log'
 input_error=99
 call score()
 call check(err==99 .and. .not.failed,'diagnostic input error remains fatal')
 input_error=0; alignment_error=98
 call score()
 call check(err==98 .and. .not.failed,'diagnostic alignment error remains fatal')
 alignment_error=0
 config%calib%targets(1)%variable='typo'
 call score()
 call check(err==20 .and. .not.failed,'missing simulated series remains fatal')
 config%calib%targets(1)%variable='streamflow'
 config%calib%targets(1)%metric='mae'
 config%calib%targets(1)%obs_transform='none'
 input_obs=qn; series(1)%values=[1._rkind,2._rkind,3._rkind]
 call score()
 call check(err==0 .and. .not.failed .and. ieee_is_nan(value),'no finite pairs: diagnostic missing')
 config%calib%targets(1)%weight=1._rkind
 call score()
 call check(err==metric_undefined .and. failed,'no finite pairs: active trial failure')
 config%calib%targets(1)%weight=0._rkind
 write_error=97
 call score()
 call check(err==97 .and. .not.failed,'output error is not swallowed for undefined diagnostic')
 close(iulog)
 write(*,'(A,I0,A)') 'PASS: production target scorer with test I/O/alignment substitutes (',checks,' assertions).'
contains
 subroutine score()
 failed=.false.
 call score_target(config,1,1,series,t,'days since 1980-01-01',value,failed,err,message)
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
