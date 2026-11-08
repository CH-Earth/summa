! Test-only substitutes for observation I/O, time alignment and output.
! The production score_target routine itself is extracted without edits.
module globalData
 use nr_type, only: rkind,i4b,nr_integerMissing,nr_realMissing
 integer(i4b),parameter :: integerMissing=nr_integerMissing
 real(rkind),parameter :: realMissing=nr_realMissing
 integer :: iulog
 integer :: ncid(1)=[1]
end module
module var_lookup
 type :: freq_type
   integer :: timestep=1
 end type
 type(freq_type) :: iLookFREQ
end module
module simulated_series
 use nr_type,only:rkind,i4b
 type :: sim_series_type
   real(rkind),allocatable :: values(:)
   character(32) :: units='m3 s-1'
 end type
 contains
 function find_simulated_series(series,variable,unit,id) result(ix)
 use nr_type,only:i8b,nr_integerMissing
 type(sim_series_type),intent(in) :: series(:)
 character(*),intent(in) :: variable
 integer(i4b),intent(in) :: unit
 integer(i8b),intent(in) :: id
 integer(i4b) :: ix
 ix=nr_integerMissing
 if(size(series)>0 .and. trim(variable)=='streamflow') ix=1
 end function
end module
module read_observations_module
 use nr_type,only:rkind,i4b
 real(rkind),allocatable :: input_obs(:)
 integer :: input_error=0
 contains
 subroutine read_observations(path,file,variable,t,bounds,v,tu,vu,err,message)
 character(*),intent(in) :: path,file,variable
 real(rkind),allocatable,intent(out) :: t(:),bounds(:,:),v(:)
 character(:),allocatable,intent(out) :: tu,vu
 integer(i4b),intent(out) :: err
 character(*),intent(out) :: message
 integer :: j,n
 err=input_error; message='observation input error'
 if(err/=0) return
 n=size(input_obs)
 allocate(t(n),bounds(2,n))
 do j=1,n
  t(j)=real(j,rkind); bounds(:,j)=[real(j-1,rkind),real(j,rkind)]
 enddo
 v=input_obs; tu='days since 1980-01-01'; vu='m3 s-1'; message=''
 end subroutine
end module
module timeseries_alignment
 use nr_type,only:rkind,i4b
 integer :: alignment_error=0
 contains
 subroutine align_timeseries(ts,vs,tsu,vsu,to,vo,tou,vou,first,last,bounds, &
                             ta,ba,ns,sa,oa,err,message)
 real(rkind),intent(in) :: ts(:),vs(:),to(:),vo(:),bounds(:,:)
 character(*),intent(in) :: tsu,vsu,tou,vou,first,last
 real(rkind),allocatable,intent(out) :: ta(:),ba(:,:),sa(:),oa(:)
 integer(i4b),allocatable,intent(out) :: ns(:)
 integer(i4b),intent(out) :: err
 character(*),intent(out) :: message
 err=alignment_error; message='alignment error'
 if(err/=0) return
 ta=to; ba=bounds; sa=vs; oa=vo
 allocate(ns(size(vo))); ns=1; message=''
 end subroutine
 subroutine select_aligned_window(t,b,n,s,o,units,first,last,err,message)
 real(rkind),allocatable,intent(inout) :: t(:),b(:,:),s(:),o(:)
 integer(i4b),allocatable,intent(inout) :: n(:)
 character(*),intent(in) :: units,first,last
 integer(i4b),intent(out) :: err
 character(*),intent(out) :: message
 err=20; message='unexpected baseline/window call in harness'
 end subroutine
end module
module series_transform
 use nr_type,only:rkind,i4b
 contains
 subroutine accumulate_series(t,v,u,out,err,message)
 real(rkind),intent(in) :: t(:),v(:)
 character(*),intent(in) :: u
 real(rkind),allocatable,intent(out) :: out(:)
 integer(i4b),intent(out) :: err
 character(*),intent(out) :: message
 err=20; message='unexpected accumulate call in harness'
 end subroutine
 subroutine remove_baseline_mean(t,v,u,first,last,out,err,message)
 real(rkind),intent(in) :: t(:),v(:)
 character(*),intent(in) :: u,first,last
 real(rkind),allocatable,intent(out) :: out(:)
 integer(i4b),intent(out) :: err
 character(*),intent(out) :: message
 err=20; message='unexpected baseline call in harness'
 end subroutine
 subroutine rebase_matched_series(t,s,o,ns,u,first,last,sb,ob,nb,err,message)
 real(rkind),intent(in) :: t(:)
 real(rkind),intent(inout) :: s(:),o(:)
 integer(i4b),intent(in) :: ns(:)
 character(*),intent(in) :: u,first,last
 real(rkind),intent(out) :: sb,ob
 integer(i4b),intent(out) :: nb,err
 character(*),intent(out) :: message
 err=20; message='unexpected rebase call in harness'
 end subroutine
 subroutine balance_between_extremes(t,s,u,to,tou,balance,window,out,err,message)
 real(rkind),intent(in) :: t(:),s(:),to(:),window
 character(*),intent(in) :: u,tou,balance
 real(rkind),allocatable,intent(out) :: out(:)
 integer(i4b),intent(out) :: err
 character(*),intent(out) :: message
 err=20; message='unexpected balance call in harness'
 end subroutine
end module
module write_evaluation_module
 use nr_type,only:rkind,i4b
 integer :: write_count=0,write_error=0
 real(rkind) :: written_objective
 contains
 subroutine write_evaluation(ncid,name,aligned,t,b,o,s,ns,tu,vu,metric,transform,obj,err,message)
 integer(i4b),intent(in) :: ncid,ns(:)
 character(*),intent(in) :: name,tu,vu,metric,transform
 logical,intent(in) :: aligned
 real(rkind),intent(in) :: t(:),b(:,:),o(:),s(:),obj
 integer(i4b),intent(out) :: err
 character(*),intent(out) :: message
 write_count=write_count+1; written_objective=obj
 err=write_error; message='test output error'
 end subroutine
end module
