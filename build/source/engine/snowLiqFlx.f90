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

module snowLiqFlx_module
USE nrtype                                    ! numerical recipes data types
USE multiconst,only:iden_ice,iden_water       ! intrinsic density of ice and water (kg m-3)
implicit none
private
public::snowLiqFlx
contains


 ! ************************************************************************************************
 ! public subroutine snowLiqFlx: compute liquid water flux through the snowpack
 ! ************************************************************************************************
 subroutine snowLiqFlx(&
                       ! input: model control
                       nSnow,                   & ! intent(in): number of snow layers
                       firstFluxCall,           & ! intent(in): the first flux call
                       ! input: forcing for the snow domain
                       scalarThroughfallRain,   & ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
                       scalarCanopyLiqDrainage, & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                       ! input: model state vector
                       mLayerVolFracLiqTrial,   & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
                       ! output: fluxes and derivatives
                       iLayerLiqFluxSnow,       & ! intent(out): vertical liquid water flux at layer interfaces (m s-1)
                       iLayerLiqFluxSnowDeriv,  & ! intent(out): derivative in vertical liquid water flux at layer interfaces (m s-1)
                       ! output: error control
                       err,message)               ! intent(out): error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:mpar_data,mvar_data                                  ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: nSnow                      ! number of snow layers
 logical(lgt),intent(in)       :: firstFluxCall              ! the first flux call
 ! input: forcing for the snow domain
 real(dp),intent(in)           :: scalarThroughfallRain      ! computed throughfall rate (kg m-2 s-1)
 real(dp),intent(in)           :: scalarCanopyLiqDrainage    ! computed drainage of liquid water (kg m-2 s-1)
 ! input: model state vector
 real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)   ! trial value of volumetric fraction of liquid water at the current iteration (-)
 ! output: fluxes and derivatives
 real(dp),intent(out)          :: iLayerLiqFluxSnow(0:)      ! vertical liquid water flux at layer interfaces (m s-1)
 real(dp),intent(out)          :: iLayerLiqFluxSnowDeriv(0:) ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)            :: cmessage                   ! error message of downwind routine
 ! ------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='snowLiqFlx/'

 ! ** calculate fluxes and derivatives for liquid water flow through snow
 call snowLiqFlx_muster(&
                        ! input: model control
                        nSnow,                                                      & ! intent(in): number of snow layers
                        firstFluxCall,                                              & ! intent(in): the first flux call
                        ! input: forcing for the snow domain
                        scalarThroughfallRain,                                      & ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
                        scalarCanopyLiqDrainage,                                    & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                        ! input: model state vector
                        mLayerVolFracLiqTrial,                                      & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
                        ! input: snow properties and parameters
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow),     & ! intent(in): volumetric ice content at the start of the time step (-)
                        mpar_data%var(iLookPARAM%Fcapil),                           & ! intent(in): capillary retention as a fraction of the total pore volume (-)
                        mpar_data%var(iLookPARAM%k_snow),                           & ! intent(in): hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
                        mpar_data%var(iLookPARAM%mw_exp),                           & ! intent(in): exponent for meltwater flow (-)
                        ! input/output: diagnostic variables -- only computed for the first iteration
                        mvar_data%var(iLookMVAR%mLayerPoreSpace)%dat,               & ! intent(inout): pore space in each snow layer (-)
                        mvar_data%var(iLookMVAR%mLayerThetaResid)%dat,              & ! intent(inout): esidual volumetric liquid water content in each snow layer (-)
                        ! output: fluxes and derivatives
                        iLayerLiqFluxSnow,                                          & ! intent(out): vertical liquid water flux at layer interfaces (m s-1)
                        iLayerLiqFluxSnowDeriv,                                     & ! intent(out): derivative in vertical liquid water flux at layer interfaces (m s-1)
                        ! output: error control
                        err,cmessage)                                                 ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine snowLiqFlx


 ! ************************************************************************************************
 ! * private subroutine: calculate fluxes and derivatives for liquid water flow through snow
 ! ************************************************************************************************
 subroutine snowLiqFlx_muster(&
                              ! input: model control
                              nSnow,                                                      & ! intent(in): number of snow layers
                              firstFluxCall,                                              & ! intent(in): the first flux call
                              ! input: forcing for the snow domain
                              scalarThroughfallRain,                                      & ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
                              scalarCanopyLiqDrainage,                                    & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                              ! input: model state vector
                              mLayerVolFracLiqTrial,                                      & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
                              ! input: snow properties and parameters
                              mLayerVolFracIce,                                           & ! intent(in): volumetric ice content at the start of the time step (-)
                              Fcapil,                                                     & ! intent(in): capillary retention as a fraction of the total pore volume (-)
                              k_snow,                                                     & ! intent(in): hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
                              mw_exp,                                                     & ! intent(in): exponent for meltwater flow (-)
                              ! input/output: diagnostic variables -- only computed for the first iteration
                              mLayerPoreSpace,                                            & ! intent(inout): pore space in each snow layer (-)
                              mLayerThetaResid,                                           & ! intent(inout): residual volumetric liquid water content in each snow layer (-)
                              ! output: fluxes and derivatives
                              iLayerLiqFluxSnow,                                          & ! intent(out): vertical liquid water flux at layer interfaces (m s-1)
                              iLayerLiqFluxSnowDeriv,                                     & ! intent(out): derivative in vertical liquid water flux at layer interfaces (m s-1)
                              ! output: error control
                              err,message)                                                  ! intent(out): error control
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: nSnow                      ! number of snow layers
 logical(lgt),intent(in)       :: firstFluxCall              ! the first flux call
 ! input: forcing for the snow domain
 real(dp),intent(in)           :: scalarThroughfallRain      ! computed throughfall rate (kg m-2 s-1)
 real(dp),intent(in)           :: scalarCanopyLiqDrainage    ! computed drainage of liquid water (kg m-2 s-1)
 ! input: model state vector
 real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)   ! trial value of volumetric fraction of liquid water at the current iteration (-)
 ! input: snow properties and parameters
 real(dp),intent(in)           :: mLayerVolFracIce(:)        ! volumetric ice content at the start of the time step (-)
 real(dp),intent(in)           :: Fcapil                     ! capillary retention as a fraction of the total pore volume (-)
 real(dp),intent(in)           :: k_snow                     ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
 real(dp),intent(in)           :: mw_exp                     ! exponent for meltwater flow (-)
 ! input/output: diagnostic variables -- only computed for the first iteration
 real(dp),intent(inout)        :: mLayerPoreSpace(:)         ! pore space in each snow layer (-)
 real(dp),intent(inout)        :: mLayerThetaResid(:)        ! residual volumetric liquid water content in each snow layer (-)
 ! output: fluxes and derivatives
 real(dp),intent(out)          :: iLayerLiqFluxSnow(0:)      ! vertical liquid water flux at layer interfaces (m s-1)
 real(dp),intent(out)          :: iLayerLiqFluxSnowDeriv(0:) ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! ---------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                  :: iLayer                     ! layer index
 real(dp)                      :: multResid                  ! multiplier for the residual water content (-)
 real(dp),parameter            :: residThrs=550._dp          ! ice density threshold to reduce residual liquid water content (kg m-3)
 real(dp),parameter            :: residScal=10._dp           ! scaling factor for residual liquid water content reduction factor (kg m-3)
 real(dp),parameter            :: maxVolIceContent=0.7_dp    ! maximum volumetric ice content to store water (-)
 real(dp)                      :: availCap                   ! available storage capacity [0,1] (-)
 real(dp)                      :: relSaturn                  ! relative saturation [0,1] (-)
 real(dp),parameter            :: dx = 1.e-8_dp              ! finite difference increment
 ! ---------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='snowLiqFlx_muster/'

 ! check that the input vectors match nSnow
 if(size(mLayerVolFracLiqTrial)/=nSnow .or. size(mLayerVolFracIce)/=nSnow .or. &
    size(iLayerLiqFluxSnow)/=nSnow+1 .or. size(iLayerLiqFluxSnowDeriv)/=nSnow+1) then
  err=20; message=trim(message)//'size mismatch of input/output vectors'; return
 endif

 ! check the meltwater exponent is >=1
 if(mw_exp<1._dp)then; err=20; message=trim(message)//'meltwater exponent < 1'; return; endif

 ! define the liquid flux at the upper boundary (m s-1)
 iLayerLiqFluxSnow(0)      = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water
 iLayerLiqFluxSnowDeriv(0) = 0._dp

 ! compute properties fixed over the time step
 if(firstFluxCall)then
  ! loop through snow layers
  do iLayer=1,nSnow
   ! compute the reduction in liquid water holding capacity at high snow density (-)
   multResid = 1._dp / ( 1._dp + exp( (mLayerVolFracIce(iLayer)*iden_ice - residThrs) / residScal) )
   ! compute the pore space (-)
   mLayerPoreSpace(iLayer)  = 1._dp - mLayerVolFracIce(iLayer)
   ! compute the residual volumetric liquid water content (-)
   mLayerThetaResid(iLayer) = Fcapil*mLayerPoreSpace(iLayer) * multResid
  end do  ! (looping through snow layers)
 endif  ! (if the first flux call)

 ! compute fluxes
 do iLayer=1,nSnow  ! (loop through snow layers)
  ! ** allow liquid water to pass through under very high density
  if(mLayerVolFracIce(iLayer) > maxVolIceContent)then ! NOTE: use start-of-step ice content, to avoid convergence problems
   iLayerLiqFluxSnow(iLayer)      = iLayerLiqFluxSnow(iLayer-1)
   iLayerLiqFluxSnowDeriv(iLayer) = 0._dp
  ! ** typical flux computations
  else
   ! check that flow occurs
   if(mLayerVolFracLiqTrial(iLayer) > mLayerThetaResid(iLayer))then
    ! compute the relative saturation (-)
    availCap  = mLayerPoreSpace(iLayer) - mLayerThetaResid(iLayer)                 ! available capacity
    relSaturn = (mLayerVolFracLiqTrial(iLayer) - mLayerThetaResid(iLayer)) / availCap    ! relative saturation
    !print*, 'mLayerVolFracLiqTrial(iLayer) = ', mLayerVolFracLiqTrial(iLayer)
    !print*, 'mLayerPoreSpace(iLayer), mLayerThetaResid(iLayer) = ', mLayerThetaResid(iLayer)
    !print*, 'iLayer, availCap, relSaturn, k_snow = ', iLayer, availCap, relSaturn, k_snow
    ! compute the flux and derivative (m s-1)
    iLayerLiqFluxSnow(iLayer)      = k_snow*relSaturn**mw_exp
    iLayerLiqFluxSnowDeriv(iLayer) = ( (k_snow*mw_exp)/availCap ) * relSaturn**(mw_exp - 1._dp)
    ! check the derivative
    !relSaturn1 = (mLayerVolFracLiqTrial(iLayer)+dx - mLayerThetaResid(iLayer)) / availCap    ! relative saturation
    !testFlux   =  k_snow*relSaturn1**mw_exp
    !write(*,'(a,1x,10(e25.10,1x))') 'iLayerLiqFluxSnow(iLayer), testFlux, iLayerLiqFluxSnowDeriv(iLayer), (testFlux - iLayerLiqFluxSnow(iLayer))/dx = ', &
    !                                 iLayerLiqFluxSnow(iLayer), testFlux, iLayerLiqFluxSnowDeriv(iLayer), (testFlux - iLayerLiqFluxSnow(iLayer))/dx
   else  ! flow does not ocur
    iLayerLiqFluxSnow(iLayer)      = 0._dp
    iLayerLiqFluxSnowDeriv(iLayer) = 0._dp
   endif  ! storage above residual content
  endif  ! check for very high density
 end do  ! loop through snow layers

 end subroutine snowLiqFlx_muster


end module snowLiqFlx_module
