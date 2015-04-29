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

module radTransfr_module
USE nrtype
implicit none
private
public::radTransfr
contains

 ! ************************************************************************************************
 ! public subroutine: compute radiation absorbed by vegetation and ground
 ! ************************************************************************************************
 subroutine radTransfr(dt,&          ! input: time step (seconds)
                       err,message)  ! output: error control
 USE data_struc,only:mpar_data,mvar_data,model_decisions     ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookDECISIONS     ! named variables for structure elements
 ! compute radiation absorbed by vegetation and ground
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! local pointers to model forcing data
 real(dp),pointer              :: sw_down                    ! downwelling shortwave radiation (W m-2)
 ! local pointers to model parameters
 real(dp),pointer              :: Frad_vis                   ! fraction of radiation in the visible part of the spectrum (-)
 ! local pointers to model variables


 ! compute the surface albedo (constant over the iterations)
 if(nSnow > 0)then
  call snowAlbedo(dt,&          ! input: time step (seconds)
                  err,cmessage) ! output: error control
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 else
  surfaceAlbedoVis = soilAlbedoVis
  surfaceAlbedoNir = soilAlbedoNir
  surfaceAlbedo    = Frad_vis*(soilAlbedoVis) + (1._dp - Frad_vis)*soilAlbedoNir
 endif
 ! compute radiation in the visible and near-infra-red part of the spectrum
 swradVis = Frad_vis*sw_down
 swradNir = (1._dp - Frad_vis)*sw_down
 ! compute the radiation absorbed at the ground surface
 swradAbsGroundVis = swradVis*(1._dp - surfaceAlbedoVis)
 swradAbsGroundNir = swradNir*(1._dp - surfaceAlbedoNir)
 swradAbsGround    = swradAbsGroundVis + swradAbsGroundNir




 ! ************************************************************************************************
 ! private subroutine: compute snow albedo
 ! ************************************************************************************************
 subroutine snowAlbedo(dt,&          ! input: time step (seconds)
                       err,message)  ! output: error control
 USE data_struc,only:mpar_data,mvar_data,model_decisions     ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookDECISIONS     ! named variables for structure elements
 USE mDecisions_module,only:funcSnowAge,BATSlike             ! named variables for albedo options
 ! compute the snow albedo
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: snw_crit                   ! critical mass necessary for albedo refreshment (kg m-2)
 real(dp),pointer              :: alb_fresh                  ! fresh snow albedo (-)
 real(dp),pointer              :: alb_dry                    ! minimum snow albedo during winter (-)
 real(dp),pointer              :: alb_wet                    ! minimum snow albedo during spring (-)
 real(dp),pointer              :: alb_decay                  ! temporal decay factor for snow albedo (s-1)
 real(dp),pointer              :: alb_scale                  ! temporal albedo scaling factor (s)
 real(dp),pointer              :: soot_load                  ! albedo decay associated with soot load (-)
 ! local pointers to model forcing data
 real(dp),pointer              :: snowfall                   ! snowfall (kg m-2 s-1)
 ! local pointers to model state variables
 real(dp),pointer              :: surfaceAlbedo              ! surface albedo (-)
 real(dp),pointer              :: surfaceTemp                ! temperature of the top layer (K)
 real(dp),pointer              :: surfaceVolFracLiq          ! volumetric fraction of liquid water the top snow layer (-)
 ! local variables
 real(dp),parameter            :: liqWinter=0.001            ! volumetric fraction of liquid water that defines "winter" or "spring" conditions
 real(dp)                      :: alb_min                    ! minimum surface albedo (-)
 real(dp)                      :: dAlb_dt                    ! change in albedo with time (s-1)
 ! initialize error control
 err=0; message="snowAlbedo/"

 ! assign pointers to the surface albedo
 surfaceAlbedo => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)

 ! assign pointers to the snowfall
 snowfall => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1)

 ! ***** compute increase in snow albedo
 if(snowfall > 0._dp)then
  ! assign pointers to albedo refreshment parameters
  snw_crit   => mpar_data%var(iLookPARAM%snw_crit)           ! critical mass necessary for albedo refreshment (kg m-2)
  alb_fresh  => mpar_data%var(iLookPARAM%alb_fresh)          ! fresh snow albedo (-)
  ! compute increase in snow albedo with time (s-1)
  dAlb_dt = min(snowfall/snw_crit, (alb_fresh - surfaceAlbedo)/dt)

 ! ***** compute decrease in snow albedo
 else

  ! identify the albedo decay method
  select case(model_decisions(iLookDECISIONS%alb_method)%iDecision)

   ! --------------------------------------------------------------------------------------------------
   ! method 1: albedo decay a function of snow age
   ! --------------------------------------------------------------------------------------------------
   case(funcSnowAge)
    ! assign pointers to model parameters
    alb_dry   => mpar_data%var(iLookPARAM%alb_dry)           ! minimum snow albedo during winter (-)
    alb_wet   => mpar_data%var(iLookPARAM%alb_wet)           ! minimum snow albedo during spring (-)
    alb_decay => mpar_data%var(iLookPARAM%alb_decay)         ! temporal decay factor for snow albedo (s-1)
    ! assign pointers to the volumetric fraction of liquid water in the upper snow layer
    surfaceVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1)
    ! compute the minimum snow albedo (-)
    if(surfaceVolFracLiq < liqWinter)then; alb_min = alb_dry ! winter conditions
    else; alb_min = alb_wet ! spring conditions
    endif
    ! compute decrease in snow albedo with time (s-1)
    dAlb_dt = -alb_decay*(surfaceAlbedo - alb_min)
    ! update albedo (use explicit solution)
    surfaceAlbedo    = surfaceAlbedo + dAlb_dt*dt
    ! assume albedo is constant in the visible and near infra-red part of the spectrum
    surfaceAlbedoVis = surfaceAlbedo
    surfaceAlbedoNir = surfaceAlbedo

   ! --------------------------------------------------------------------------------------------------
   ! method 2: albedo decay based on a BATS-like approach, with destructive metamorphism + soot content
   ! --------------------------------------------------------------------------------------------------
   case(BATSlike)






    err=10; message=trim(message)//'"batslike" method not implemented yet'; return

   ! check for unknown albedo decay method
   case default
    err=10; message=trim(message)//"unknownOption"; return

  end select  ! (identifying the albedo decay method)

 endif ! (if snowing: switch between albedo increase or decay)


 end subroutine snowAlbedo


end module radTransfr_module
