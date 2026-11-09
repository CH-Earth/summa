program test_boxcox_config
 use nr_type,only:rkind
 use calibration_types,only:calib_info
 use tomlf_all,only:toml_table,toml_array,put_real,put_string
 use parsers_under_test,only:parse_calibration_settings,parse_calibration_targets
 use legacy_under_test,only:config_info,set_default_calibration_target
 use, intrinsic :: ieee_arithmetic
 implicit none
 type(toml_table),target :: parent
 type(toml_table),pointer :: table
 type(toml_array),target :: targets
 type(calib_info) :: calib
 type(config_info) :: config
 character(:),allocatable :: path,variable
 character(1024) :: message
 real(rkind) :: qn,inf
 integer :: err,checks,j
 checks=0
 qn=ieee_value(0._rkind,ieee_quiet_nan)
 inf=ieee_value(0._rkind,ieee_positive_inf)
 table=>parent
 path='obs/'; variable='q_obs'
 allocate(targets%tables(3))
 parent%targets=>targets
 do j=1,3
   call put_string(targets%tables(j),'obs_file','test.nc')
   call put_string(targets%tables(j),'obs_transform','boxcox')
 enddo
 call parse_calibration_settings(table,calib,err,message)
 call check(err==0 .and. calib%boxcox_lambda==0.25_rkind,'no global key keeps historical default')
 call parse_calibration_targets(table,calib,path,variable,err,message)
 call check(err==0 .and. all(calib%targets%boxcox_lambda==0.25_rkind),'all omitted targets use default')
 call put_real(parent,'boxcox_lambda',0.4_rkind)
 call parse_calibration_settings(table,calib,err,message)
 call check(err==0 .and. calib%boxcox_lambda==0.4_rkind,'global exponent parsed')
 call put_real(targets%tables(2),'boxcox_lambda',0._rkind)
 call put_real(targets%tables(3),'boxcox_lambda',-0.5_rkind)
 call parse_calibration_targets(table,calib,path,variable,err,message)
 call check(err==0,'explicit target list parsed')
 call check(calib%targets(1)%boxcox_lambda==0.4_rkind,'target inherits global value')
 call check(calib%targets(2)%boxcox_lambda==0._rkind,'zero override retained exactly')
 call check(calib%targets(3)%boxcox_lambda==-0.5_rkind,'negative override retained')
 call check(all(calib%targets%obs_transform=='boxcox'),'target transformation unchanged')
 ! Reparse a reused calibration object: overrides cannot leak to another target.
 call put_real(parent,'boxcox_lambda',0.6_rkind)
 call parse_calibration_settings(table,calib,err,message)
 call parse_calibration_targets(table,calib,path,variable,err,message)
 call check(err==0 .and. calib%targets(1)%boxcox_lambda==0.6_rkind,'reparse inherits the new general default')
 call check(calib%targets(2)%boxcox_lambda==0._rkind,'reparse preserves explicit zero')
 call put_string(parent,'boxcox_lambda','bad')
 call parse_calibration_settings(table,calib,err,message)
 call check(err==20,'global wrong type rejected')
 call put_real(parent,'boxcox_lambda',qn)
 call parse_calibration_settings(table,calib,err,message)
 call check(err==20,'global NaN rejected')
 call put_real(parent,'boxcox_lambda',inf)
 call parse_calibration_settings(table,calib,err,message)
 call check(err==20,'global infinity rejected')
 call put_real(parent,'boxcox_lambda',0.4_rkind)
 call parse_calibration_settings(table,calib,err,message)
 call put_string(targets%tables(2),'boxcox_lambda','bad')
 call parse_calibration_targets(table,calib,path,variable,err,message)
 call check(err==20,'target wrong type rejected')
 call put_real(targets%tables(2),'boxcox_lambda',qn)
 call parse_calibration_targets(table,calib,path,variable,err,message)
 call check(err==20,'target NaN rejected')
 call put_real(targets%tables(2),'boxcox_lambda',inf)
 call parse_calibration_targets(table,calib,path,variable,err,message)
 call check(err==20,'target infinity rejected')
 call put_real(targets%tables(2),'boxcox_lambda',0._rkind)
 call parse_calibration_targets(table,calib,path,variable,err,message)
 call check(err==0,'valid configuration recovers after parse error')
 config%calib%metric='nse'; config%calib%obs_transform='boxcox'
 config%calib%boxcox_lambda=0.7_rkind
 config%obs%obs_file='test.nc'; config%obs%obs_path=path
 call set_default_calibration_target(config,err,message)
 call check(err==0 .and. size(config%calib%targets)==1,'legacy implicit target synthesized')
 call check(config%calib%targets(1)%boxcox_lambda==0.7_rkind,'legacy implicit target copies general exponent')
 config%calib%targets=calib%targets
 config%calib%boxcox_lambda=0.9_rkind
 call set_default_calibration_target(config,err,message)
 call check(err==0 .and. size(config%calib%targets)==3,'explicit targets still take precedence')
 call check(config%calib%targets(2)%boxcox_lambda==0._rkind,'explicit exponent not overwritten by legacy adapter')
 deallocate(targets%tables)
 write(*,'(A,I0,A)') 'PASS: production parsers with test TOML accessors and legacy target adapter (',checks,' assertions).'
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
