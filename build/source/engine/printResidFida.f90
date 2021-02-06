

module printResidFida_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! define access to state variables to print
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
! privacy
implicit none
private
public::printResidFida
contains

 ! **********************************************************************************************************
 ! public subroutine printResidFida: compute the residual vector
 ! **********************************************************************************************************
 subroutine printResidFida(&
                        ! input: model control
                        nSnow,                     & ! intent(in):    number of snow layers
                        nSoil,                     & ! intent(in):    number of soil layers
                        nLayers,                   & ! intent(in):    total number of layers
                        ! input: data structures
                        indx_data,                 & ! intent(in):    index data
                        ! output
                        rAdd,                      & ! intent(out):   additional (sink) terms on the RHS of the state equation
                        rVec)                       ! intent(out):   residual vector

 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 integer(i4b),intent(in)         :: nSnow                     ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                   ! total number of layers in the snow+soil domain
 type(var_ilength),intent(in)    :: indx_data                 ! indices defining model states and layers
 ! output
 real(dp),intent(in)            :: rAdd(:)                   ! additional (sink) terms on the RHS of the state equation
 real(qp),intent(in)            :: rVec(:)   ! NOTE: qp      ! residual vector
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 integer(i4b)                    :: iLayer                    ! index of layer within the snow+soil domain
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! link to the necessary variables for the residual computations
 associate(&
  ! number of state variables of a specific type
  nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
  nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
  nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
  ! model indices
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy energy state variable
  ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
  ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in): [i4b]    index of water storage in the aquifer
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
  ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
  ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
  ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
  ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in): [i4b(:)] named variables defining the type of hydrology states in snow+soil domain
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                 & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
 ) ! association to necessary variables for the residual computations
 ! --------------------------------------------------------------------------------------------------------------------------------


 if(ixVegNrg/=integerMissing) print *, 'rAdd(ixVegNrg) = ', rAdd(ixVegNrg) 

 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   
   select case( layerType(iLayer) )
    case(iname_snow)
     print *, 'rAdd( ixSnowSoilNrg(iLayer) ) = ', rAdd( ixSnowSoilNrg(iLayer) ) 
    case(iname_soil); print *, 'rAdd( ixSnowSoilNrg(iLayer) ) = ', rAdd( ixSnowSoilNrg(iLayer) )  
   end select
  end do  
 endif
 
 if(nSoilOnlyHyd>0)then   
  do concurrent (iLayer=1:nSoil,ixSoilOnlyHyd(iLayer)/=integerMissing)   
   print *, 'rAdd( ixSoilOnlyHyd(iLayer) ) = ', rAdd( ixSoilOnlyHyd(iLayer) ) 
  end do  
 endif

 if(ixCasNrg/=integerMissing) print *, 'rVec(ixCasNrg) = ', rVec(ixCasNrg) 
 if(ixVegNrg/=integerMissing) print *, 'rVec(ixVegNrg) = ', rVec(ixVegNrg) 
 if(ixVegHyd/=integerMissing)then     
   print *, 'rVec(ixVegHyd) = ', rVec(ixVegHyd) 
 endif

 if(nSnowSoilNrg>0)then 
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)  
   print *, 'rVec( ixSnowSoilNrg(iLayer) ) = ', rVec( ixSnowSoilNrg(iLayer) ) 
  end do  
 endif

 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)  
 print *, 'rVec( ixSnowSoilHyd(iLayer) ) = ', rVec( ixSnowSoilHyd(iLayer) ) 
  end do  
 endif

  if(ixAqWat/=integerMissing)  print *, ' rVec(ixAqWat) = ', rVec(ixAqWat)

 end associate

 end subroutine printResidFida

end module printResidFida_module
