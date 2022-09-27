! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module bigAquifer_module
! -----------------------------------------------------------------------------------------------------------

! numerical recipes data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! physical constants
USE multiconst,only:&
                    LH_vap,  & ! latent heat of vaporization   (J kg-1)
                    iden_water ! intrinsic density of water    (kg m-3)

! -----------------------------------------------------------------------------------------------------------
implicit none
private
public::bigAquifer
contains


! ***************************************************************************************************************
! public subroutine bigAquifer: compute aquifer water fluxes and their derivatives
! ***************************************************************************************************************
subroutine bigAquifer(&
                    ! input: state variables and fluxes
                    scalarAquiferStorageTrial,    & ! intent(in):  trial value of aquifer storage (m)
                    scalarCanopyTranspiration,    & ! intent(in):  canopy transpiration (kg m-2 s-1)
                    scalarSoilDrainage,           & ! intent(in):  soil drainage (m s-1)
                    ! input: pre-computed derivatives
                    dCanopyTrans_dCanWat,         & ! intent(in): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
                    dCanopyTrans_dTCanair,        & ! intent(in): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                    dCanopyTrans_dTCanopy,        & ! intent(in): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
                    dCanopyTrans_dTGround,        & ! intent(in): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
                    ! input: diagnostic variables and parameters
                    mpar_data,                    & ! intent(in):  model parameter structure
                    diag_data,                    & ! intent(in):  diagnostic variable structure
                    ! output: fluxes
                    scalarAquiferTranspire,       & ! intent(out): transpiration loss from the aquifer (m s-1)
                    scalarAquiferRecharge,        & ! intent(out): recharge to the aquifer (m s-1)
                    scalarAquiferBaseflow,        & ! intent(out): total baseflow from the aquifer (m s-1)
                    dBaseflow_dAquifer,           & ! intent(out): change in baseflow flux w.r.t. aquifer storage (s-1)
                    ! output: derivatives in transpiration w.r.t. canopy state variables
                    dAquiferTrans_dTCanair,        & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
                    dAquiferTrans_dTCanopy,        & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
                    dAquiferTrans_dTGround,        & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. ground temperature
                    dAquiferTrans_dCanWat,         & ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy total water
                    ! output: error control
                    err,message)                    ! intent(out): error control
  ! named variables
  USE var_lookup,only:iLookDIAG              ! named variables for structure elements
  USE var_lookup,only:iLookPARAM             ! named variables for structure elements
  ! data types
  USE data_types,only:var_dlength            ! x%var(:)%dat   (rkind)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: state variables, fluxes, and parameters
  real(rkind),intent(in)              :: scalarAquiferStorageTrial    ! trial value of aquifer storage (m)
  real(rkind),intent(in)              :: scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
  real(rkind),intent(in)              :: scalarSoilDrainage           ! soil drainage (m s-1)
  ! input: pre-computed derivatves
  real(rkind),intent(in)              :: dCanopyTrans_dCanWat          ! derivative in canopy transpiration w.r.t. canopy total water content (s-1)
  real(rkind),intent(in)              :: dCanopyTrans_dTCanair         ! derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  real(rkind),intent(in)              :: dCanopyTrans_dTCanopy         ! derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
  real(rkind),intent(in)              :: dCanopyTrans_dTGround         ! derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
  ! input: diagnostic variables and parameters
  type(var_dlength),intent(in)     :: mpar_data                    ! model parameters
  type(var_dlength),intent(in)     :: diag_data                    ! diagnostic variables for a local HRU
  ! output: fluxes
  real(rkind),intent(out)             :: scalarAquiferTranspire       ! transpiration loss from the aquifer (m s-1)
  real(rkind),intent(out)             :: scalarAquiferRecharge        ! recharge to the aquifer (m s-1)
  real(rkind),intent(out)             :: scalarAquiferBaseflow        ! total baseflow from the aquifer (m s-1)
  real(rkind),intent(out)             :: dBaseflow_dAquifer           ! change in baseflow flux w.r.t. aquifer storage (s-1)
  ! output: derivatives in transpiration w.r.t. canopy state variables
  real(rkind),intent(inout)           :: dAquiferTrans_dTCanair      ! derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
  real(rkind),intent(inout)           :: dAquiferTrans_dTCanopy      ! derivatives in the aquifer transpiration flux w.r.t. canopy temperature
  real(rkind),intent(inout)           :: dAquiferTrans_dTGround      ! derivatives in the aquifer transpiration flux w.r.t. ground temperature
  real(rkind),intent(inout)           :: dAquiferTrans_dCanWat       ! derivatives in the aquifer transpiration flux w.r.t. canopy total water
  ! output: error control
  integer(i4b),intent(out)         :: err                          ! error code
  character(*),intent(out)         :: message                      ! error message
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  real(rkind)                         :: aquiferTranspireFrac         ! fraction of total transpiration that comes from the aquifer (-)
  real(rkind)                         :: xTemp                        ! temporary variable (-)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='bigAquifer/'

  ! make association between local variables and the information in the data structures
  associate(&
    ! model diagnostic variables: contribution of the aquifer to transpiration
    scalarTranspireLim     => diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1),     & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
    scalarAquiferRootFrac  => diag_data%var(iLookDIAG%scalarAquiferRootFrac)%dat(1),  & ! intent(in): [dp] fraction of roots below the lowest soil layer (-)
    scalarTranspireLimAqfr => diag_data%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1), & ! intent(in): [dp] transpiration limiting factor for the aquifer (-)
    ! model parameters: baseflow flux
    aquiferBaseflowRate    => mpar_data%var(iLookPARAM%aquiferBaseflowRate)%dat(1),   & ! intent(in): [dp] tbaseflow rate when aquiferStorage = aquiferScaleFactor (m s-1)
    aquiferScaleFactor     => mpar_data%var(iLookPARAM%aquiferScaleFactor)%dat(1),    & ! intent(in): [dp] scaling factor for aquifer storage in the big bucket (m)
    aquiferBaseflowExp     => mpar_data%var(iLookPARAM%aquiferBaseflowExp)%dat(1)     & ! intent(in): [dp] baseflow exponent (-)
    )  ! associating local variables with the information in the data structures

    ! compute aquifer transpiration (m s-1)
    aquiferTranspireFrac   = scalarAquiferRootFrac*scalarTranspireLimAqfr/scalarTranspireLim   ! fraction of total transpiration that comes from the aquifer (-)
    scalarAquiferTranspire = aquiferTranspireFrac*scalarCanopyTranspiration/iden_water         ! aquifer transpiration (kg m-2 s-1 --> m s-1)
    ! derivatives in transpiration w.r.t. canopy state variables
    dAquiferTrans_dCanWat  = aquiferTranspireFrac*dCanopyTrans_dCanWat /iden_water
    dAquiferTrans_dTCanair = aquiferTranspireFrac*dCanopyTrans_dTCanair/iden_water
    dAquiferTrans_dTCanopy = aquiferTranspireFrac*dCanopyTrans_dTCanopy/iden_water
    dAquiferTrans_dTGround = aquiferTranspireFrac*dCanopyTrans_dTGround/iden_water

    ! compute aquifer recharge (transfer variables -- included for generality for basin-wide aquifer)
    scalarAquiferRecharge = scalarSoilDrainage ! m s-1

    ! compute the aquifer baseflow (m s-1)
    xTemp                 = scalarAquiferStorageTrial/aquiferScaleFactor
    scalarAquiferBaseflow = aquiferBaseflowRate*(xTemp**aquiferBaseflowExp)

    ! compute the derivative in the net aquifer flux
    dBaseflow_dAquifer    = -(aquiferBaseflowExp*aquiferBaseflowRate*(xTemp**(aquiferBaseflowExp - 1._rkind)))/aquiferScaleFactor

    ! end association to data in structures
  end associate

end subroutine bigAquifer

end module bigAquifer_module
