module stateLimit_module
USE nrtype
implicit none
private
! define bounds for enthalpy
real(dp),save,dimension(2)   :: boundsEnthalpySoil=huge(dp)
real(dp),save,dimension(2)   :: boundsEnthalpySnow=huge(dp)
! define bounds for the bulk density of water
real(dp),save,dimension(2)   :: boundsBdenWatrSoil=huge(dp)
real(dp),save,dimension(2)   :: boundsBdenWatrSnow=huge(dp)
! define bounds for snow depth
real(dp),save,dimension(2)   :: boundsDepthSnow=huge(dp)
! define bounds for snow albedo
real(dp),save,dimension(2)   :: boundsAlbedoSnow=huge(dp)
! make routines public, but keep data private
public::get_limits
public::stateLimit
contains

 ! ************************************************************************************************
 ! new subroutine: identify reasonable limits for the model state variables
 ! ************************************************************************************************
 subroutine get_limits(err,message)
 ! store reasonable limits for the model state variables as derived model parameters
 USE multiconst, only: Tfreeze, &                   ! freezing point of water (K)
                       wat_den,ice_den              ! density of water and ice (kg m-3)
 USE data_struc,only:mpar_data                      ! parameter data structures
 USE var_lookup,only:iLookPARAM                     ! named variables for structure elements
 USE ConvE2Temp_module,only:temp2ethpy              ! convert temperature to enthalpy
 implicit none
 ! dummy variables
 integer(i4b),intent(out)                  :: err             ! error code
 character(*),intent(out)                  :: message         ! error message
 ! local pointers to data in structures
 real(dp),pointer                          :: soil_dens       ! soil density (kg m-3)
 real(dp),pointer                          :: soil_pore       ! soil porosity (-)
 real(dp),pointer                          :: soil_frz        ! freezing curve parameter for soil (K-1)
 real(dp),pointer                          :: snow_frz        ! freezing curve parameter for snow(K-1)
 real(dp),pointer                          :: zmax            ! maximum layer depth (m)
 real(dp),pointer                          :: alb_fresh       ! fresh snow albedo (-)
 ! local variables
 real(dp),parameter                        :: Tmin=100._dp    ! minimum allowable temperature
 real(dp),parameter                        :: Tmax=500._dp    ! maximum allowable temperature
 real(dp),parameter                        :: eps=epsilon(Tmin) ! very small number
 ! initialize error control
 err=0; message="get_limits/"
 ! assign pointers to model parameters
 soil_dens => mpar_data%var(iLookPARAM%soil_dens_bulk)        ! hydraulic conductivity of the soil (kg m-2 s-1)
 soil_pore => mpar_data%var(iLookPARAM%soil_pore)             ! soil porosity (-)
 soil_frz  => mpar_data%var(iLookPARAM%soil_frz)              ! freezing curve parameter for soil(K-1)
 snow_frz  => mpar_data%var(iLookPARAM%snow_frz)              ! freezing curve parameter for snow(K-1)
 zmax      => mpar_data%var(iLookPARAM%zmax)                  ! maximum layer depth (m)
 alb_fresh => mpar_data%var(iLookPARAM%alb_fresh)             ! fresh snow albedo (-)
 ! define bounds for enthalpy
 boundsEnthalpySoil(1) = temp2ethpy(Tmin,soil_dens,soil_pore*wat_den,soil_frz)
 boundsEnthalpySoil(2) = temp2ethpy(Tmax,soil_dens,soil_pore*wat_den,soil_frz)
 boundsEnthalpySnow(1) = temp2ethpy(Tmin,0._dp,wat_den,snow_frz)
 boundsEnthalpySnow(2) = temp2ethpy(Tfreeze,0._dp,wat_den,snow_frz)
 ! define bounds for the bulk density of water
 boundsBdenWatrSoil(1) = eps
 boundsBdenWatrSoil(2) = soil_pore*wat_den
 boundsBdenWatrSnow(1) = 10._dp
 boundsBdenWatrSnow(2) = ice_den
 ! define bounds for snow depth
 boundsDepthSnow(1) = eps
 boundsDepthSnow(2) = zmax*10._dp
 ! define bounds for snow albedo
 boundsAlbedoSnow(1) = 0.1_dp
 boundsAlbedoSnow(2) = alb_fresh
 end subroutine get_limits

 ! ************************************************************************************************
 ! new subroutine: populate state bound vectors
 ! ************************************************************************************************
 subroutine stateLimit(err,message)
 USE data_struc,only:indx_data,ix_soil,ix_snow                ! model indices
 USE var_lookup,only:iLookINDEX,iLookMVAR                     ! named variables to define elements in the structure
 USE data_struc,only:stateLookup                              ! look-up table for identifying model states
 USE data_struc,only:xboundLower,xboundUpper                  ! state bounds
 implicit none
 ! dummy variables
 integer(i4b),intent(out)                  :: err             ! error code
 character(*),intent(out)                  :: message         ! error message
 ! local pointers to data in structures
 integer(i4b),pointer                      :: layerType(:)    ! type of the layer (ix_soil or ix_snow)
 integer(i4b),pointer                      :: nLayTotl        ! total number of layers
 ! local variables
 integer(i4b)                              :: nState          ! number of state variables
 real(dp),dimension(:),allocatable         :: tmpvec          ! temporary vector for all layers
 real(dp),dimension(:),allocatable         :: xlimit          ! temporary vector for all states
 integer(i4b)                              :: ibound          ! loop thru bounds
 ! initialize error control
 err=0; message="stateLimit/" 
 ! perform some basic checks
 if(.not.associated(stateLookup))then; err=10; message='stateLookup undefined'; return; endif
 nState=size(stateLookup)
 ! get the total number of layers and the layer type
 nLayTotl  => indx_data%var(iLookINDEX%nLayers)%dat(1)
 layerType => indx_data%var(iLookINDEX%layerType)%dat
 ! allocate space for the bound vectors
 allocate(xboundLower(nState),xboundUpper(nState),xlimit(nState),tmpvec(nLayTotl),stat=err)
 if(err/=0)then; err=20; message='problem allocating space for bound vectors'; return; endif
 xboundLower(:)=-huge(tmpvec); xboundUpper(:)=-huge(tmpvec); xlimit=-huge(tmpvec)
 ! loop through the lower and upper bounds
 do ibound=1,2
  ! populate the bounds vectors for enthalpy
  where(layerType==ix_soil) tmpvec=boundsEnthalpySoil(ibound)
  where(layerType==ix_snow) tmpvec=boundsEnthalpySnow(ibound)
  xlimit = unpack(tmpvec,stateLookup==iLookMVAR%mLayerEnthalpy,xlimit)
  ! populate the bounds vectors for the bulk density of water
  where(layerType==ix_soil) tmpvec=boundsBdenWatrSoil(ibound)
  where(layerType==ix_snow) tmpvec=boundsBdenWatrSnow(ibound)
  xlimit = unpack(tmpvec,stateLookup==iLookMVAR%mLayerBulkDenWat,xlimit)
  ! populate the bounds vectors for the snow depth
  where(stateLookup==iLookMVAR%mLayerDepth) xlimit=boundsDepthSnow(ibound)
  ! populate the bounds vectors for snow albedo
  where(stateLookup==iLookMVAR%scalarAlbedo) xlimit=boundsAlbedoSnow(ibound)
  ! ...and, save in the final bound vectors
  if(ibound==1) xboundLower=xlimit
  if(ibound==2) xboundUpper=xlimit
 end do  ! (looping thru bound vectors)
 ! deallocate space
 deallocate(xlimit,tmpvec,stat=err)
 if(err/=0)then; err=20; message='problem deallocating space for bound vectors'; return; endif
 end subroutine stateLimit

end module stateLimit_module
