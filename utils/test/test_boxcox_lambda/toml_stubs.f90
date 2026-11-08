! TEST ONLY: in-memory substitutes for the tomlf_all get_value/table interface.
! Production parser routines are extracted unchanged; no TOML text is parsed here.
module tomlf_all
 use nr_type,only:i4b,i8b,rkind,lgt
 implicit none
 type :: toml_key
   character(:),allocatable :: key
 end type
 type :: entry
   character(64) :: key=''
   integer :: kind=0
   character(256) :: text=''
   real(rkind) :: number=0._rkind
   integer(i8b) :: integer_value=0_i8b
   logical :: boolean=.false.
 end type
 type :: toml_table
   type(entry),allocatable :: entries(:)
   type(toml_array),pointer :: targets=>null()
 contains
   procedure :: get_keys
   procedure :: has_key
 end type
 type :: toml_array
   type(toml_table),pointer :: tables(:)=>null()
 end type
 interface get_value
   module procedure get_string,get_real,get_integer,get_integer8,get_logical,get_array,get_table
 end interface
 interface len
   module procedure array_len
 end interface
contains
 integer function index_of(table,key) result(k)
 class(toml_table),intent(in) :: table
 character(*),intent(in) :: key
 integer :: j
 k=0
 if(.not.allocated(table%entries)) return
 do j=1,size(table%entries)
   if(trim(table%entries(j)%key)==key)then
     k=j; return
   endif
 enddo
 end function
 logical function has_key(table,key) result(found)
 class(toml_table),intent(in) :: table
 character(*),intent(in) :: key
 found=index_of(table,key)>0
 if(key=='target') found=associated(table%targets)
 end function
 subroutine get_keys(table,keys)
 class(toml_table),intent(in) :: table
 type(toml_key),allocatable,intent(out) :: keys(:)
 integer :: j,n,m
 n=0
 if(allocated(table%entries)) n=size(table%entries)
 m=n
 if(associated(table%targets)) m=m+1
 allocate(keys(m))
 do j=1,n
   keys(j)%key=trim(table%entries(j)%key)
 enddo
 if(m>n) keys(m)%key='target'
 end subroutine
 subroutine put(table,e)
 type(toml_table),intent(inout) :: table
 type(entry),intent(in) :: e
 integer :: k
 k=index_of(table,trim(e%key))
 if(k>0)then
   table%entries(k)=e
 else if(allocated(table%entries))then
   table%entries=[table%entries,e]
 else
   table%entries=[e]
 endif
 end subroutine
 subroutine put_real(table,key,value)
 type(toml_table),intent(inout) :: table
 character(*),intent(in) :: key
 real(rkind),intent(in) :: value
 type(entry) :: e
 e%key=key; e%kind=2; e%number=value
 call put(table,e)
 end subroutine
 subroutine put_string(table,key,value)
 type(toml_table),intent(inout) :: table
 character(*),intent(in) :: key,value
 type(entry) :: e
 e%key=key; e%kind=1; e%text=value
 call put(table,e)
 end subroutine
 subroutine get_string(table,key,value,stat)
 type(toml_table),pointer,intent(in) :: table
 character(*),intent(in) :: key
 character(:),allocatable,intent(out) :: value
 integer,intent(out) :: stat
 integer :: k
 stat=1; k=index_of(table,key)
 if(k==0) return
 if(table%entries(k)%kind/=1) return
 value=trim(table%entries(k)%text); stat=0
 end subroutine
 subroutine get_real(table,key,value,stat)
 type(toml_table),pointer,intent(in) :: table
 character(*),intent(in) :: key
 real(rkind),intent(inout) :: value
 integer,intent(out) :: stat
 integer :: k
 stat=1; k=index_of(table,key)
 if(k==0) return
 if(table%entries(k)%kind/=2) return
 value=table%entries(k)%number; stat=0
 end subroutine
 subroutine get_integer(table,key,value,stat)
 type(toml_table),pointer,intent(in) :: table
 character(*),intent(in) :: key
 integer(i4b),intent(inout) :: value
 integer,intent(out) :: stat
 integer :: k
 stat=1; k=index_of(table,key)
 if(k==0) return
 if(table%entries(k)%kind/=3) return
 value=int(table%entries(k)%integer_value,i4b); stat=0
 end subroutine
 subroutine get_integer8(table,key,value,stat)
 type(toml_table),pointer,intent(in) :: table
 character(*),intent(in) :: key
 integer(i8b),intent(inout) :: value
 integer,intent(out) :: stat
 integer :: k
 stat=1; k=index_of(table,key)
 if(k==0) return
 if(table%entries(k)%kind/=3) return
 value=table%entries(k)%integer_value; stat=0
 end subroutine
 subroutine get_logical(table,key,value,stat)
 type(toml_table),pointer,intent(in) :: table
 character(*),intent(in) :: key
 logical(lgt),intent(inout) :: value
 integer,intent(out) :: stat
 integer :: k
 stat=1; k=index_of(table,key)
 if(k==0) return
 if(table%entries(k)%kind/=4) return
 value=table%entries(k)%boolean; stat=0
 end subroutine
 subroutine get_array(table,key,value,requested,stat)
 type(toml_table),pointer,intent(in) :: table
 character(*),intent(in) :: key
 type(toml_array),pointer,intent(out) :: value
 logical,optional,intent(in) :: requested
 integer,intent(out) :: stat
 nullify(value); stat=1
 if(key/='target') return
 if(associated(table%targets))then
   value=>table%targets; stat=0
 endif
 end subroutine
 subroutine get_table(array,k,value,stat)
 type(toml_array),pointer,intent(in) :: array
 integer,intent(in) :: k
 type(toml_table),pointer,intent(out) :: value
 integer,intent(out) :: stat
 value=>array%tables(k); stat=0
 end subroutine
 integer function array_len(array) result(n)
 type(toml_array),intent(in) :: array
 n=size(array%tables)
 end function
end module
