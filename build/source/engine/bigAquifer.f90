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
                    LH_vap,  &      ! latent heat of vaporization   (J kg-1)
                    iden_water      ! intrinsic density of water    (kg m-3)
! -----------------------------------------------------------------------------------------------------------
implicit none
private
public :: bigAquifer
contains
! ***************************************************************************************************************
! public subroutine bigAquifer: compute aquifer water fluxes and their derivatives
! ***************************************************************************************************************
subroutine bigAquifer(&
                      ! input: state variables, fluxes, and pre-computed derivatives
                      in_bigAquifer,                & ! intent(in):    state variables, fluxes, and pre-computed derivatives
                      ! input: diagnostic variables and parameters
                      mpar_data,                    & ! intent(in):    model parameter structure
                      diag_data,                    & ! intent(in):    diagnostic variable structure
                      ! input-output: derivatives in transpiration w.r.t. canopy state variables
                      io_bigAquifer,                & ! intent(inout): derivatives in transpiration w.r.t. canopy state variables
                      ! output: fluxes and error control
                      out_bigAquifer)                 ! intent(out):   fluxes and error control
  ! named variables
  USE var_lookup,only:iLookDIAG                     ! named variables for structure elements
  USE var_lookup,only:iLookPARAM                    ! named variables for structure elements
  ! data types
  USE data_types,only:var_dlength                   ! x%var(:)%dat [rkind]
  USE data_types,only:in_type_bigAquifer            ! derived typ for intent(in) arguments
  USE data_types,only:io_type_bigAquifer            ! derived typ for intent(inout) arguments
  USE data_types,only:out_type_bigAquifer           ! derived typ for intent(out) arguments
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: state variables, fluxes, and pre-computed derivatives
  type(in_type_bigAquifer),intent(in)    :: in_bigAquifer                ! state variables, fluxes, and pre-computed derivatives
  ! input: diagnostic variables and parameters
  type(var_dlength),intent(in)           :: mpar_data                    ! model parameters
  type(var_dlength),intent(in)           :: diag_data                    ! diagnostic variables for a local HRU
  ! input-output: derivatives in transpiration w.r.t. canopy state variables
  type(io_type_bigAquifer),intent(inout) :: io_bigAquifer                ! derivatives in transpiration w.r.t. canopy state variables
  ! output: fluxes and error control
  type(out_type_bigAquifer),intent(out)  :: out_bigAquifer               ! fluxes and error control
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  real(rkind)                            :: aquiferTranspireFrac         ! fraction of total transpiration that comes from the aquifer (-)
  real(rkind)                            :: xTemp                        ! temporary variable (-)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! make association between local variables and the information in the data structures
  associate(&
    ! input: state variables, fluxes, and parameters
    scalarAquiferStorageTrial => in_bigAquifer % scalarAquiferStorageTrial,   & ! intent(in): [dp] trial value of aquifer storage (m)
    scalarCanopyTranspiration => in_bigAquifer % scalarCanopyTranspiration,   & ! intent(in): [dp] canopy transpiration (kg m-2 s-1)
    scalarSoilDrainage        => in_bigAquifer % scalarSoilDrainage,          & ! intent(in): [dp] soil drainage (m s-1)
    ! input: pre-computed derivatves
    dCanopyTrans_dCanWat  => in_bigAquifer % dCanopyTrans_dCanWat,    & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
    dCanopyTrans_dTCanair => in_bigAquifer % dCanopyTrans_dTCanair,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dTCanopy => in_bigAquifer % dCanopyTrans_dTCanopy,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dTGround => in_bigAquifer % dCanopyTrans_dTGround,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
    ! input: model diagnostic variables: contribution of the aquifer to transpiration
    scalarTranspireLim     => diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1),     & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
    scalarAquiferRootFrac  => diag_data%var(iLookDIAG%scalarAquiferRootFrac)%dat(1),  & ! intent(in): [dp] fraction of roots below the lowest soil layer (-)
    scalarTranspireLimAqfr => diag_data%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1), & ! intent(in): [dp] transpiration limiting factor for the aquifer (-)
    ! input: model parameters: baseflow flux
    aquiferBaseflowRate    => mpar_data%var(iLookPARAM%aquiferBaseflowRate)%dat(1),   & ! intent(in): [dp] tbaseflow rate when aquiferStorage = aquiferScaleFactor (m s-1)
    aquiferScaleFactor     => mpar_data%var(iLookPARAM%aquiferScaleFactor)%dat(1),    & ! intent(in): [dp] scaling factor for aquifer storage in the big bucket (m)
    aquiferBaseflowExp     => mpar_data%var(iLookPARAM%aquiferBaseflowExp)%dat(1),    & ! intent(in): [dp] baseflow exponent (-)
    ! input-output: derivatives in transpiration w.r.t. canopy state variables
    dAquiferTrans_dTCanair => io_bigAquifer % dAquiferTrans_dTCanair, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
    dAquiferTrans_dTCanopy => io_bigAquifer % dAquiferTrans_dTCanopy, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
    dAquiferTrans_dTGround => io_bigAquifer % dAquiferTrans_dTGround, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
    dAquiferTrans_dCanWat  => io_bigAquifer % dAquiferTrans_dCanWat,  & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
    ! output: fluxes
    scalarAquiferTranspire => out_bigAquifer % scalarAquiferTranspire,& ! intent(out): transpiration loss from the aquifer (m s-1)
    scalarAquiferRecharge  => out_bigAquifer % scalarAquiferRecharge, & ! intent(out): recharge to the aquifer (m s-1)
    scalarAquiferBaseflow  => out_bigAquifer % scalarAquiferBaseflow, & ! intent(out): total baseflow from the aquifer (m s-1)
    dBaseflow_dAquifer     => out_bigAquifer % dBaseflow_dAquifer,    & ! intent(out): change in baseflow flux w.r.t. aquifer storage (s-1)
    ! output: error control
    err                    => out_bigAquifer % err,                   & ! intent(out): error code
    message                => out_bigAquifer % cmessage               & ! intent(out): error message
    )  ! end associating local variables with the information in the data structures
    err=0; message='bigAquifer/' ! initialize error control

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
    if (xTemp<0._rkind) xTemp = 0._rkind ! otherwise will give NaN in next line
    scalarAquiferBaseflow = aquiferBaseflowRate*(xTemp**aquiferBaseflowExp)

    ! compute the derivative in the net aquifer flux
    dBaseflow_dAquifer    = -(aquiferBaseflowExp*aquiferBaseflowRate*(xTemp**(aquiferBaseflowExp - 1._rkind)))/aquiferScaleFactor

  end associate ! end association to data in structure

end subroutine bigAquifer

end module bigAquifer_module
