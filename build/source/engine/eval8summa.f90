! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
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

module eval8summa_module
! data types
USE nrtype
! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
! layer types
USE data_struc,only:ix_soil,ix_snow ! named variables for snow and soil
! provide access to the derived types to define the data structures
USE data_struc,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)
! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookMVAR       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
implicit none
private
! common variables
real(dp),parameter :: valueMissing=-9999._dp ! missing value
contains

 ! **********************************************************************************************************
 ! public subroutine eval8summa: compute the residual vector and the Jacobian matrix
 ! **********************************************************************************************************
 subroutine eval8summa(&
                       ! input
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       nState,                  & ! intent(in):    total number of state variables
                       stateVec,                & ! intent(in):    model state vector
                       mpar_data,               & ! intent(in):    model parameters for a local HRU
                       mvar_data,               & ! intent(in):    model variables for a local HRU
                       indx_data,               & ! intent(in):    indices defining model states and layers
                       ! output
                       rVec,                    & ! intent(out):   residual vector
                       aJac,                    & ! intent(out):   analytical Jacobian matrix (either band or full)
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE getVectorz_module,only:varExtract            ! extract variables from the state vector
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 integer(i4b),intent(in)         :: nSnow         ! number of snow layers
 integer(i4b),intent(in)         :: nSoil         ! number of soil layers
 integer(i4b),intent(in)         :: nLayers       ! total number of layers
 integer(i4b),intent(in)         :: nState        ! total number of state variables
 real(dp),intent(in)             :: stateVec(:)   ! model state vector (mixed units)
 type(var_d),intent(in)          :: mpar_data     ! model parameters for a local HRU
 type(var_dlength),intent(in)    :: mvar_data     ! model variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data     ! indices defining model states and layers
 ! output
 real(dp),intent(out)            :: rVec(:)       ! residual vector
 real(dp),intent(out)            :: aJac(:,:)     ! analytical Jacobian matrix (either band or full)
 integer(i4b),intent(out)        :: err           ! error code
 character(*),intent(out)        :: message       ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! state variables
 real(dp)                        :: scalarCanairTempTrial        ! trial value for temperature of the canopy air space (K)
 real(dp)                        :: scalarCanopyTempTrial        ! trial value for temperature of the vegetation canopy (K)
 real(dp)                        :: scalarCanopyWatTrial         ! trial value for liquid water storage in the canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerTempTrial              ! trial value for temperature of layers in the snow and soil domains (K)
 real(dp),dimension(nLayers)     :: mLayerVolFracWatTrial        ! trial value for volumetric fraction of total water (-)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadTrial        ! trial value for matric head (m)
 ! diagnostic variables
 real(dp)                        :: scalarCanopyLiqTrial         ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyIceTrial         ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerVolFracLiqTrial        ! trial value for volumetric fraction of liquid water (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracIceTrial        ! trial value for volumetric fraction of ice (-)
 ! fraction of liquid water
 real(dp)                        :: fracLiqVeg                   ! fraction of liquid water on the vegetation canopy (-)
 real(dp),dimension(nSnow)       :: fracLiqSnow                  ! volumetric fraction of water in each snow layer (-)
 ! error control
 character(len=256)              :: cMessage                     ! error message of downwind routine
 





 ! extract summa variables from the model state vector
 call varExtract(&
                 ! input
                 stateVec,                      & ! intent(in):    model state vector (mixed units)
                 mpar_data,                     & ! intent(in):    model parameters
                 mvar_data,                     & ! intent(in):    model variables for a local HRU
                 indx_data,                     & ! intent(in):    indices defining model states and layers
                 ! output: variables for the vegetation canopy
                 fracLiqVeg,                    & ! intent(out):   fraction of liquid water on the vegetation canopy (-)
                 scalarCanairTempTrial,         & ! intent(out):   trial value of canopy air temperature (K)
                 scalarCanopyTempTrial,         & ! intent(out):   trial value of canopy temperature (K)
                 scalarCanopyWatTrial,          & ! intent(out):   trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,          & ! intent(out):   trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,          & ! intent(out):   trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 fracLiqSnow,                   & ! intent(out):   volumetric fraction of water in each snow layer (-)
                 mLayerTempTrial,               & ! intent(out):   trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,         & ! intent(out):   trial vector of volumetric total water content (-) 
                 mLayerVolFracLiqTrial,         & ! intent(out):   trial vector of volumetric liquid water content (-) 
                 mLayerVolFracIceTrial,         & ! intent(out):   trial vector of volumetric ice water content (-) 
                 mLayerMatricHeadTrial,         & ! intent(out):   trial vector of matric head (m)
                 ! output: error control 
                 err,cmessage)                    ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------


 ! set values to missing
 rVec(:) = valueMissing
 aJac(:,:) = valueMissing






 ! start of internal subroutines 
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 !contains









 end subroutine eval8summa
end module eval8summa_module
