program test
USE nrtype
implicit none
type var_lookup
 integer(i4b)            :: i1=1001
 integer(i4b)            :: i2=1002
endtype var_lookup
type(var_lookup),parameter :: ilookup
print*,ilookup%i1
print*,ilookup%i2








stop
end
