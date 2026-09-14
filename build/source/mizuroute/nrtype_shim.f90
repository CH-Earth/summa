module nrtype

  use nr_type, only: &
      I8B, I4B, I2B, I1B, &
      SP, DP, SPC, DPC, LGT

  implicit none

  integer,      parameter :: wp = DP

  integer(I4B), parameter :: strLen     = 256
  integer(I4B), parameter :: FileStrLen = 300
  integer(I4B), parameter :: gageStrLen = 30

end module nrtype
