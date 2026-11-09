! TEST ONLY: capture the metadata written through the NetCDF Fortran interface.
! This module does NOT create or read a NetCDF file.
module netcdf
 use nr_type,only:rkind
 implicit none
 integer,parameter :: NF90_NOERR=0,NF90_CLOBBER=0,NF90_INT=4,NF90_DOUBLE=6,NF90_CHAR=2,NF90_GLOBAL=-1
 real(rkind),parameter :: NF90_FILL_DOUBLE=9.9692099683868690e36_rkind
 integer :: nvars=0,ndims=0,ngroups=1000
 character(64) :: names(100)=''
 real(rkind) :: lambda_attrs(100)=0._rkind,fill_attrs(100)=0._rkind
 logical :: has_lambda(100)=.false.,has_fill(100)=.false.
 real(rkind),allocatable :: lambda_vector(:),weight_vector(:)
 interface nf90_put_att
   module procedure att_text,att_real,att_integer
 end interface
 interface nf90_put_var
   module procedure var_real_scalar,var_real_vector,var_real_matrix,var_int_vector,var_text
 end interface
contains
 integer function nf90_create(path,mode,ncid) result(err)
 character(*),intent(in) :: path
 integer,intent(in) :: mode
 integer,intent(out) :: ncid
 ncid=1; err=0
 end function
 integer function nf90_def_dim(ncid,name,n,dimid) result(err)
 integer,intent(in) :: ncid,n
 character(*),intent(in) :: name
 integer,intent(out) :: dimid
 ndims=ndims+1; dimid=ndims; err=0
 end function
 integer function nf90_def_var(ncid,name,xtype,dimids,varid) result(err)
 integer,intent(in) :: ncid,xtype
 character(*),intent(in) :: name
 integer,optional,intent(in) :: dimids(:)
 integer,intent(out) :: varid
 nvars=nvars+1; varid=nvars; names(varid)=name; err=0
 end function
 integer function nf90_enddef(ncid) result(err)
 integer,intent(in) :: ncid
 err=0
 end function
 integer function nf90_redef(ncid) result(err)
 integer,intent(in) :: ncid
 err=0
 end function
 integer function nf90_close(ncid) result(err)
 integer,intent(in) :: ncid
 err=0
 end function
 function nf90_strerror(code) result(message)
 integer,intent(in) :: code
 character(32) :: message
 write(message,'(A,I0)') 'test NetCDF status ',code
 end function
 integer function nf90_inq_ncid(ncid,name,gid) result(err)
 integer,intent(in) :: ncid
 character(*),intent(in) :: name
 integer,intent(out) :: gid
 gid=-1; err=1
 end function
 integer function nf90_def_grp(ncid,name,gid) result(err)
 integer,intent(in) :: ncid
 character(*),intent(in) :: name
 integer,intent(out) :: gid
 ngroups=ngroups+1; gid=ngroups; err=0
 end function
 integer function nf90_inq_varid(ncid,name,varid) result(err)
 integer,intent(in) :: ncid
 character(*),intent(in) :: name
 integer,intent(out) :: varid
 integer :: j
 err=1; varid=-1
 do j=nvars,1,-1
   if(trim(names(j))==name)then
     varid=j; err=0; return
   endif
 enddo
 end function
 integer function att_text(ncid,varid,name,value) result(err)
 integer,intent(in) :: ncid,varid
 character(*),intent(in) :: name,value
 err=0
 end function
 integer function att_integer(ncid,varid,name,value) result(err)
 integer,intent(in) :: ncid,varid,value
 character(*),intent(in) :: name
 err=0
 end function
 integer function att_real(ncid,varid,name,value) result(err)
 integer,intent(in) :: ncid,varid
 character(*),intent(in) :: name
 real(rkind),intent(in) :: value
 err=0
 if(varid<1) return
 if(name=='boxcox_lambda')then
   has_lambda(varid)=.true.; lambda_attrs(varid)=value
 endif
 if(name=='_FillValue')then
   has_fill(varid)=.true.; fill_attrs(varid)=value
 endif
 end function
 integer function var_real_scalar(ncid,varid,value,start,count) result(err)
 integer,intent(in) :: ncid,varid
 real(rkind),intent(in) :: value
 integer,optional,intent(in) :: start(:),count(:)
 err=0
 end function
 integer function var_real_vector(ncid,varid,values,start,count) result(err)
 integer,intent(in) :: ncid,varid
 real(rkind),intent(in) :: values(:)
 integer,optional,intent(in) :: start(:),count(:)
 err=0
 if(names(varid)=='target_boxcox_lambda') lambda_vector=values
 if(names(varid)=='target_weight') weight_vector=values
 end function
 integer function var_real_matrix(ncid,varid,values,start,count) result(err)
 integer,intent(in) :: ncid,varid
 real(rkind),intent(in) :: values(:,:)
 integer,optional,intent(in) :: start(:),count(:)
 err=0
 end function
 integer function var_int_vector(ncid,varid,values,start,count) result(err)
 integer,intent(in) :: ncid,varid,values(:)
 integer,optional,intent(in) :: start(:),count(:)
 err=0
 end function
 integer function var_text(ncid,varid,value,start,count) result(err)
 integer,intent(in) :: ncid,varid
 character(*),intent(in) :: value
 integer,optional,intent(in) :: start(:),count(:)
 err=0
 end function
end module
