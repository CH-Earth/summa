module noahMP_veg_module
USE nrtype
! named variables for snow and soil
USE data_struc,only:ix_soil,ix_snow        ! named variables for snow and soil
! -------------------------------------------------------------------------------------------------
implicit none
private
public::noahMP_veg
! number of soil and snow layers
integer(i4b)                  :: nSoil     ! number of soil layers
integer(i4b)                  :: nSnow     ! number of snow layers
integer(i4b)                  :: nLayers   ! total number of layers
contains


 ! ************************************************************************************************
 ! new subroutine: compute energy and mass fluxes for vegetation
 ! ************************************************************************************************
 subroutine noahMP_veg(err,message)                             ! intent(out): error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 ! output
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! internal
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 ! initialize error control
 err=0; message="noahMP_veg/"

 ! define number of layers
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)  ! number of snow layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)                   ! total number of layers

 ! compute vegetation phenology





 end subroutine noahMP_veg


end module noahMP_veg_module
