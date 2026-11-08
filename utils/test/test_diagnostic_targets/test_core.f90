program test_core
 use nr_type, only: rkind,i4b
 use metrics, only: compute_metric,metric_undefined
 use summa_type, only: config_info
 use summa_simulation, only: scalarize_objectives
 use settings_under_test, only: check_search_settings
 use, intrinsic :: ieee_arithmetic
 implicit none
 real(rkind) :: q(5),s(5),v,expected,eps,qn,inf,values(2),zero(0)
 type(config_info) :: config
 integer :: err,checks
 character(1024) :: message
 q=[0._rkind,0.1_rkind,1._rkind,2._rkind,4._rkind]
 s=[0.01_rkind,0.2_rkind,0.8_rkind,2.2_rkind,3.8_rkind]
 qn=ieee_value(0._rkind,ieee_quiet_nan)
 inf=ieee_value(0._rkind,ieee_positive_inf)
 checks=0
 call compute_metric(q,q,'kge','none',v,err,message)
 call check(err==0 .and. abs(v-1._rkind)<1.e-14_rkind,'perfect KGE')
 call compute_metric(q,s,'nse','log',v,err,message)
 eps=sum(q)/(100*size(q))
 expected=1._rkind-sum((log(q+eps)-log(s+eps))**2)/ &
         sum((log(q+eps)-sum(log(q+eps))/size(q))**2)
 call check(err==0 .and. abs(v-expected)<1.e-14_rkind,'NSE-log formula/epsilon')
 call compute_metric([2._rkind],[1._rkind],'mae','none',v,err,message)
 call check(err==0 .and. v==1._rkind,'single-pair MAE is valid')
 call compute_metric([2._rkind],[1._rkind],'rmse','none',v,err,message)
 call check(err==0 .and. v==1._rkind,'single-pair RMSE is valid')
 call compute_metric([1._rkind,qn,2._rkind],[1._rkind,5._rkind,3._rkind], &
                     'mae','none',v,err,message)
 call check(err==0 .and. v==0.5_rkind,'finite-pair mask is preserved')
 call undefined([1._rkind],[1._rkind],'nse','none','one-pair efficiency')
 call undefined([0._rkind,0._rkind],[0._rkind,0._rkind],'nse','log','zero-flow log')
 call undefined(q,[-1._rkind,0.1_rkind,1._rkind,2._rkind,4._rkind], &
                'nse','log','nonpositive log argument')
 call undefined([2._rkind,2._rkind],[1._rkind,2._rkind],'nse','none','constant observed NSE')
 call undefined(spread(0.1_rkind,1,7),spread(0.2_rkind,1,7),'kge','none','constant rounding case')
 call undefined(q,spread(1._rkind,1,5),'kge','none','constant simulated KGE')
 call undefined([-1._rkind,1._rkind],[0._rkind,2._rkind],'kge','none','zero observed mean')
 call undefined([1._rkind,2._rkind],[-1._rkind,1._rkind],'kgep','none','zero simulated mean')
 call undefined([qn,qn],[1._rkind,2._rkind],'mae','none','no finite pairs')
 call undefined(zero,zero,'nse','log','empty vectors')
 call undefined([-1._rkind,1._rkind],[1._rkind,2._rkind],'mae','boxcox','negative Box-Cox')
 call undefined([-1._rkind,1._rkind],[1._rkind,2._rkind],'mae','0.5','fractional negative base')
 call undefined([0._rkind,0._rkind],[0._rkind,0._rkind],'mae','-1','reciprocal zero')
 call compute_metric(q,q,'mae','typo',v,err,message)
 call check(err==20,'bad transformation is a configuration error')
 call compute_metric(zero,zero,'nse','typo',v,err,message)
 call check(err==20,'bad transformation not hidden by missing pairs')
 call compute_metric(zero,zero,'typo','none',v,err,message)
 call check(err==20,'bad metric not hidden by missing pairs')
 call compute_metric(q,q(1:4),'nse','none',v,err,message)
 call check(err==20,'dimension mismatch is a configuration error')
 call compute_metric(q,q,'nse','NaN',v,err,message)
 call check(err==20,'NaN power exponent rejected')
 call compute_metric([-1._rkind,1._rkind],[-2._rkind,2._rkind], &
                     'mae','2',v,err,message)
 call check(err==0 .and. v==3._rkind,'integer power accepts negative arguments')
 allocate(config%calib%targets(2))
 config%calib%algorithm='dds'
 config%calib%n_samples=20
 config%calib%targets(1)%name='kge'
 config%calib%targets(1)%metric='kge'
 config%calib%targets(1)%weight=1._rkind
 config%calib%targets(2)%name='nse_log'
 config%calib%targets(2)%metric='nse'
 config%calib%targets(2)%weight=0._rkind
 call check_search_settings(config,err,message)
 call check(err==0,'DDS with diagnostic accepted')
 values=[0.8_rkind,qn]
 v=scalarize_objectives(config,values)
 call check(v==0.8_rkind,'zero-weight NaN skipped')
 values(2)=inf
 v=scalarize_objectives(config,values)
 call check(v==0.8_rkind,'zero-weight infinity skipped')
 values(2)=-9999._rkind
 v=scalarize_objectives(config,values)
 call check(v==0.8_rkind,'zero-weight sentinel skipped')
 config%calib%targets(2)%metric='rmse'
 config%calib%targets(2)%weight=0.25_rkind
 values=[0.8_rkind,2._rkind]
 v=scalarize_objectives(config,values)
 call check(abs(v-0.3_rkind)<1.e-14_rkind,'active weighted error metric oriented correctly')
 config%calib%targets%weight=0._rkind
 call check_search_settings(config,err,message)
 call check(err/=0,'all-zero optimization refused')
 config%calib%targets(1)%weight=1._rkind
 config%calib%targets(2)%weight=-1._rkind
 call check_search_settings(config,err,message)
 call check(err/=0,'negative weight refused')
 config%calib%targets(2)%weight=qn
 call check_search_settings(config,err,message)
 call check(err/=0,'NaN weight refused')
 config%calib%targets(2)%weight=inf
 call check_search_settings(config,err,message)
 call check(err/=0,'infinite weight refused')
 config%calib%targets(2)%weight=0._rkind
 config%calib%algorithm='nsga2'
 config%calib%nsga2%population_size=10
 call check_search_settings(config,err,message)
 call check(err/=0,'NSGA-II diagnostic guard')
 config%calib%targets(2)%weight=1._rkind
 call check_search_settings(config,err,message)
 call check(err==0,'positive-weight NSGA-II unaffected')
 write(*,'(A,I0,A)') 'PASS: metrics, scalarizer and optimization settings (',checks,' assertions).'
contains
 subroutine check(ok,label)
 logical,intent(in) :: ok
 character(*),intent(in) :: label
 checks=checks+1
 if(.not.ok)then
   write(*,'(A)') 'FAIL: '//label//' '//trim(message)
   stop 1
 endif
 end subroutine
 subroutine undefined(obs,sim,metric,transform,label)
 real(rkind),intent(in) :: obs(:),sim(:)
 character(*),intent(in) :: metric,transform,label
 call compute_metric(obs,sim,metric,transform,v,err,message)
 call check(err==metric_undefined,trim(label)//' status')
 call check(ieee_is_nan(v),trim(label)//' NaN')
 end subroutine
end program
