program test_boxcox_metrics
 use nr_type,only:rkind,i4b
 use metrics,only:compute_metric,metric_undefined
 use calibration_types,only:calib_info,target_info
 use, intrinsic :: ieee_arithmetic
 implicit none
 type(calib_info) :: calib
 type(target_info) :: target
 real(rkind) :: q(5),s(5),qo(5),so(5),v,ref,old,lambda,qn,inf,empty(0)
 real(rkind),parameter :: lambdas(7)=[0.2_rkind,0.5_rkind,1._rkind,2._rkind, &
                                     -0.5_rkind,-1._rkind,0._rkind]
 character(4),parameter :: names(5)=[character(4)::'kge','kgep','nse','mae','rmse']
 integer :: i,j,err,referr,checks
 character(1024) :: message,refmessage
 checks=0
 q=[0.01_rkind,0.2_rkind,1._rkind,2._rkind,10._rkind]
 s=[0.02_rkind,0.15_rkind,0.9_rkind,2.2_rkind,9._rkind]
 qn=ieee_value(0._rkind,ieee_quiet_nan)
 inf=ieee_value(0._rkind,ieee_positive_inf)
 call check(calib%boxcox_lambda==0.25_rkind,'calibration default is 0.25')
 call check(target%boxcox_lambda==0.25_rkind,'target default is 0.25')
 do i=1,size(names)
   ! Exact regression against the old arithmetic, including a zero-flow case.
   do j=1,2
     qo=q; so=s
     if(j==2)then
       qo(1)=0._rkind; so(1)=0._rkind
     endif
     call compute_metric(qo,so,trim(names(i)),'boxcox',old,err=err,message=message)
     call check(err==0,'lambda omitted with keyword errors is valid')
     call compute_metric(qo,so,trim(names(i)),'boxcox',v,boxcox_lambda=0.25_rkind,err=err,message=message)
     call check(err==0,'explicit default call remains valid')
     call check(v==old,'omitted and explicit default agree bit for bit')
     qo=(qo**0.25_rkind-1._rkind)/0.25_rkind
     so=(so**0.25_rkind-1._rkind)/0.25_rkind
     call compute_metric(qo,so,trim(names(i)),'none',ref,err=referr,message=refmessage)
     call check(referr==0 .and. v==ref,'default equals historical Box-Cox expression exactly')
   enddo
   do j=1,size(lambdas)
     lambda=lambdas(j)
     if(lambda==0._rkind)then
       qo=log(q); so=log(s)
     else
       qo=(q**lambda-1._rkind)/lambda
       so=(s**lambda-1._rkind)/lambda
     endif
     call compute_metric(qo,so,trim(names(i)),'none',ref,err=referr,message=refmessage)
     call compute_metric(q,s,trim(names(i)),'boxcox',v,boxcox_lambda=lambda,err=err,message=message)
     call check(err==0 .and. referr==0,'custom exponent gives valid metric')
     call check(abs(v-ref)<2.e-13_rkind*max(1._rkind,abs(ref)),'custom exponent matches direct formula')
   enddo
 enddo
 ! Near-zero exponents must approach unshifted log, without cancellation.
 do j=1,2
   lambda=1.e-12_rkind
   if(j==2) lambda=-lambda
   call compute_metric(log(q),log(s),'nse','none',ref,err=referr,message=refmessage)
   call compute_metric(q,s,'nse','boxcox',v,boxcox_lambda=lambda,err=err,message=message)
   call check(err==0 .and. abs(v-ref)<1.e-11_rkind,'stable near-zero exponent')
 enddo
 call compute_metric([0._rkind,1._rkind],[1._rkind,4._rkind],'mae','boxcox',v,boxcox_lambda=0.5_rkind,err=err, &
                     message=message)
 call check(err==0 .and. v==2._rkind,'zero values allowed for positive lambda')
 call compute_metric(q,s,'mae','boxcox',v,boxcox_lambda=0._rkind,err=err,message=message)
 call compute_metric(q,s,'mae','log',ref,err=referr,message=refmessage)
 call check(err==0 .and. referr==0 .and. abs(v-ref)>1.e-4_rkind,'unshifted lambda=0 differs from legacy shifted log')
 do j=1,2
   lambda=0._rkind
   if(j==2) lambda=-0.5_rkind
   call compute_metric([0._rkind,1._rkind],[1._rkind,2._rkind],'mae','boxcox',v,boxcox_lambda=lambda,err=err, &
                       message=message)
   call check(err==metric_undefined .and. ieee_is_nan(v),'zero invalid for nonpositive lambda')
 enddo
 call compute_metric([-1._rkind,1._rkind],[1._rkind,2._rkind],'mae','boxcox',v,boxcox_lambda=2._rkind,err=err, &
                     message=message)
 call check(err==metric_undefined,'negative inputs rejected even for integer lambda')
 call compute_metric(q,s,'mae','boxcox',v,boxcox_lambda=qn,err=err,message=message)
 call check(err==20 .and. ieee_is_nan(v),'NaN lambda is a configuration error')
 call compute_metric(q,s,'mae','boxcox',v,boxcox_lambda=inf,err=err,message=message)
 call check(err==20,'positive infinite lambda rejected')
 call compute_metric(q,s,'mae','boxcox',v,boxcox_lambda=-inf,err=err,message=message)
 call check(err==20,'negative infinite lambda rejected')
 call compute_metric(empty,empty,'mae','boxcox',v,boxcox_lambda=qn,err=err,message=message)
 call check(err==20,'bad lambda not hidden by empty data')
 call compute_metric(empty,empty,'mae','boxcox',v,boxcox_lambda=0.5_rkind,err=err,message=message)
 call check(err==metric_undefined,'empty valid configuration remains an undefined metric')
 call compute_metric(q,s,'nse','none',v,boxcox_lambda=-0.7_rkind,err=err,message=message)
 call compute_metric(q,s,'nse','none',ref,err=referr,message=refmessage)
 call check(err==0 .and. v==ref,'lambda does not affect non-Box-Cox transformations')
 call compute_metric([1.e300_rkind],[1._rkind],'mae','boxcox',v,boxcox_lambda=2._rkind,err=err,message=message)
 call check(err==metric_undefined,'power overflow is reported without an FPE')
 call compute_metric([1.e300_rkind],[1._rkind],'mae','boxcox',v,boxcox_lambda=huge(1._rkind),err=err, &
                     message=message)
 call check(err==metric_undefined,'lambda*log(Q) overflow is guarded')
 call compute_metric([0._rkind],[1._rkind],'mae','boxcox',v,boxcox_lambda=1.e-310_rkind,err=err,message=message)
 call check(err==metric_undefined,'zero-flow reciprocal overflow is guarded')
 call compute_metric([1.e300_rkind],[1._rkind],'mae','boxcox',v,boxcox_lambda=-huge(1._rkind),err=err, &
                     message=message)
 call check(err==0 .and. ieee_is_finite(v),'extreme negative exponent avoids spurious overflow')
 call compute_metric([exp(7.1_rkind)],[1._rkind],'mae','boxcox',v,boxcox_lambda=100._rkind,err=err, &
                     message=message)
 ref=exp(710._rkind-log(100._rkind))-0.01_rkind
 call check(err==0 .and. abs(v/ref-1._rkind)<1.e-12_rkind,'scaled result survives intermediate power overflow')
 write(*,'(A,I0,A)') 'PASS: Box-Cox defaults, formulas, limits and domain checks (',checks,' assertions).'
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
end program
