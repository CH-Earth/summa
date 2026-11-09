program test_boxcox_output
 use nr_type,only:rkind
 use netcdf
 use calibration_types,only:target_info
 use parameter_search,only:parameter_spec
 use calibration_output_under_test,only:create_calibration_output
 use write_evaluation_module,only:write_evaluation
 use, intrinsic :: ieee_arithmetic
 implicit none
 type(target_info) :: targets(4)
 type(parameter_spec) :: spec
 integer :: err,ncid,varid,checks,before
 character(1024) :: message
 real(rkind) :: empty(0),bounds(2,0),qn
 integer :: counts(0)
 checks=0
 allocate(spec%params(0),spec%ordered(0))
 targets(1)%name='kge'; targets(1)%weight=1._rkind
 targets(2)%name='bc_default'; targets(2)%obs_transform='boxcox'; targets(2)%weight=0._rkind
 targets(3)%name='bc_zero'; targets(3)%obs_transform='boxcox'; targets(3)%boxcox_lambda=0._rkind
 targets(4)%name='bc_negative'; targets(4)%obs_transform='boxcox'; targets(4)%boxcox_lambda=-0.5_rkind
 targets(3:4)%weight=0._rkind
 call create_calibration_output('test.nc',spec,4,1,0,0,'test',targets,'dds',0,ncid,err,message)
 call check(err==0,'calibration file definition succeeds')
 call check(allocated(lambda_vector),'calibration writer emits lambda values')
 call check(size(lambda_vector)==4,'one lambda slot per target')
 call check(lambda_vector(1)==NF90_FILL_DOUBLE,'non-Box-Cox target has missing exponent')
 call check(lambda_vector(2)==0.25_rkind,'default exponent recorded')
 call check(lambda_vector(3)==0._rkind,'zero exponent not confused with missing value')
 call check(lambda_vector(4)==-0.5_rkind,'negative exponent recorded')
 call check(all(weight_vector==[1._rkind,0._rkind,0._rkind,0._rkind]),'diagnostic weights unchanged')
 err=nf90_inq_varid(ncid,'target_boxcox_lambda',varid)
 call check(err==0 .and. has_fill(varid),'lambda variable declares _FillValue')
 call check(fill_attrs(varid)==NF90_FILL_DOUBLE,'lambda fill agrees with stored missing values')
 call write_evaluation(ncid,'bc_half',.false.,empty,bounds,empty,empty,counts,'days','m3 s-1','nse','boxcox', &
                       0.8_rkind,boxcox_lambda=0.5_rkind,err=err,message=message)
 call check(err==0,'production evaluation writer accepts exponent')
 call check(has_lambda(nvars) .and. lambda_attrs(nvars)==0.5_rkind,'evaluation objective records custom exponent')
 call write_evaluation(ncid,'bc_default',.false.,empty,bounds,empty,empty,counts,'days','m3 s-1','nse','boxcox', &
                       0.8_rkind,err=err,message=message)
 call check(err==0 .and. lambda_attrs(nvars)==0.25_rkind,'writer with omitted lambda records historical default')
 call write_evaluation(ncid,'bc_zero',.false.,empty,bounds,empty,empty,counts,'days','m3 s-1','nse','boxcox', &
                       0.8_rkind,boxcox_lambda=0._rkind,err=err,message=message)
 call check(err==0 .and. has_lambda(nvars) .and. lambda_attrs(nvars)==0._rkind,'evaluation writer records lambda zero')
 call write_evaluation(ncid,'plain',.false.,empty,bounds,empty,empty,counts,'days','m3 s-1','nse','none', &
                       0.8_rkind,boxcox_lambda=0.5_rkind,err=err,message=message)
 call check(err==0 .and. .not.has_lambda(nvars),'no misleading lambda attribute for plain target')
 qn=ieee_value(0._rkind,ieee_quiet_nan)
 before=nvars
 call write_evaluation(ncid,'invalid',.false.,empty,bounds,empty,empty,counts,'days','m3 s-1','nse','boxcox', &
                       0.8_rkind,boxcox_lambda=qn,err=err,message=message)
 call check(err==20 .and. nvars==before,'bad evaluation exponent rejected before defining metadata')
 write(*,'(A,I0,A)') 'PASS: production metadata writers with test NetCDF interface (',checks,' assertions).'
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
