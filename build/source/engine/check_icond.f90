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

module check_icond_module
USE nr_type

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing real number

! constants
USE globalData,only:maxVolIceContent ! snow maximum volumetric ice content to store water (-)
USE globalData,only:verySmall        ! a small number

! access domain types
USE globalData,only:upland           ! horizontal domain type for upland areas
USE globalData,only:glacCln1         ! first horizontal domain type for glacier clean areas
USE globalData,only:glacCln2         ! second horizontal domain type for glacier clean areas
USE globalData,only:glacDbr          ! horizontal domain type for glacier debris areas
USE globalData,only:wetland          ! horizontal domain type for wetland areas

USE globalData,only:icefrz_mult      ! freezing curve scaling factor multipier of snow to ice, closer to a step function since ice does not hold water

implicit none
private
public::check_icond
contains

 ! ************************************************************************************************
 ! public subroutine check_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine check_icond(nGRU,                          & ! number of GRUs
                        bvarData,                      & ! model basin variables
                        progData,                      & ! intent(inout): model prognostic (state) variables
                        diagData,                      & ! intent(inout): model diagnostic variables
                        mparData,                      & ! intent(in):    model parameters
                        indxData,                      & ! intent(in):    layer index data
                        lookupData,                    & ! intent(in):    lookup table data
                        attrData,                      & ! intent(in):    model attributes
                        checkEnthalpy,                 & ! intent(in):    flag if to check enthalpy for consistency
                        no_dom_vars,                   & ! intent(out):   flag that domain variables are not in initial conditions
                        no_ice_vars,                   & ! intent(out):   flag that glacier ice variables are not in initial conditions
                        no_ablfrac,                    & ! intent(out):   flag that glacier ablation fraction variable is not in initial conditions
                        no_icond_enth,                 & ! intent(in):    flag that enthalpy not in initial conditions
                        use_lookup,                    & ! intent(in):    flag to use the lookup table for soil enthalpy                             
                        err,message)                     ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nr_type
 USE var_lookup,only:iLookBVAR                           ! variable lookup structure
 USE var_lookup,only:iLookPARAM                          ! variable lookup structure
 USE var_lookup,only:iLookPROG                           ! variable lookup structure
 USE var_lookup,only:iLookDIAG                           ! variable lookup structure
 USE var_lookup,only:iLookINDEX                          ! variable lookup structure
 USE var_lookup,only:iLookATTR                           ! variable lookup structure
 USE globalData,only:gru_struc                           ! gru-hru mapping structures
 USE data_types,only:gru_doubleVec                       ! gru double precision structure no hru no domain
 USE data_types,only:gru_hru_doubleVec                   ! hru double precision structure no domain
 USE data_types,only:gru_hru_dom_doubleVec               ! full double precision structure
 USE data_types,only:gru_hru_dom_intVec                  ! full integer structure
 USE data_types,only:gru_hru_dom_z_vLookup               ! full lookup structure
 USE data_types,only:gru_hru_double                      ! hru double precision structure no domain no depth
 USE globalData,only:iname_soil,iname_snow,iname_glce,iname_lake ! named variables to describe the type of layer
 USE multiconst,only:&
                       LH_fus,    &                      ! latent heat of fusion                (J kg-1)
                       iden_ice,  &                      ! intrinsic density of ice             (kg m-3)
                       iden_water,&                      ! intrinsic density of liquid water    (kg m-3)
                       gravity,   &                      ! gravitational acceleration           (m s-2)
                       Tfreeze                           ! freezing point of pure water         (K)
 USE updatState_module,only:updatSnLaGl                  ! update snow lake glce states
 USE updatState_module,only:updatSoil                    ! update soil states
 USE convertEnthalpyTemp_module,only:T2enthTemp_cas      ! convert temperature to enthalpy for canopy air space
 USE convertEnthalpyTemp_module,only:T2enthTemp_veg      ! convert temperature to enthalpy for vegetation
 USE convertEnthalpyTemp_module,only:T2enthTemp_snLaGl   ! convert temperature to enthalpy for snow, lake, and ice
 USE convertEnthalpyTemp_module,only:T2enthTemp_soil     ! convert temperature to enthalpy for soil
 
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 integer(i4b),intent(in)                   :: nGRU                       ! number of grouped response units
 type(gru_doubleVec),intent(inout)         :: bvarData                   ! basin variables
 type(gru_hru_dom_doubleVec),intent(inout) :: progData                   ! prognostic vars
 type(gru_hru_dom_doubleVec),intent(inout) :: diagData                   ! diagnostic vars
 type(gru_hru_dom_doubleVec),intent(in)    :: mparData                   ! parameters
 type(gru_hru_dom_intVec),intent(in)       :: indxData                   ! layer indexes
 type(gru_hru_dom_z_vLookup),intent(in)    :: lookupData                 ! lookup table data
 type(gru_hru_double),intent(in)           :: attrData                   ! attributes
 logical(lgt),intent(in)                   :: checkEnthalpy              ! if true either need enthTemp as starting residual value, or for state variable initialization
 logical(lgt),intent(out)                  :: no_dom_vars                ! if true, no domain area/elevation in initial conditions
 logical(lgt),intent(out)                  :: no_ice_vars                ! if true, no glacier ice variables in initial conditions
 logical(lgt),intent(out)                  :: no_ablfrac                 ! if true, no glacier ablation fraction in initial conditions
 logical(lgt),intent(in)                   :: no_icond_enth              ! if true, no enthalpy in initial conditions
 logical(lgt),intent(in)                   :: use_lookup                 ! flag to use the lookup table for soil enthalpy, otherwise use hypergeometric function
 integer(i4b),intent(out)                  :: err                        ! error code
 character(*),intent(out)                  :: message                    ! returned error message
 ! locals
 character(len=256)                        :: cmessage                   ! downstream error message
 integer(i4b)                              :: iGRU,iHRU,iDOM,iGlac       ! loop index
 ! temporary variables for realism checks
 integer(i4b)                              :: iLayer                     ! index of model layer
 integer(i4b)                              :: iSoil                      ! index of soil layer
 real(rkind)                               :: fLiq                       ! fraction of liquid water on the vegetation canopy (-)
 real(rkind)                               :: vGn_m                      ! van Genutchen "m" parameter (-)
 real(rkind)                               :: scalarTheta                ! liquid water equivalent of total water [liquid water + ice] (-)
 real(rkind)                               :: h1,h2                      ! used to check depth and height are consistent
 real(rkind)                               :: kappa                      ! constant in the freezing curve function (m K-1)
 integer(i4b)                              :: nSnow                      ! number of snow layers
 integer(i4b)                              :: nLake                      ! number of lake layers
 integer(i4b)                              :: nSoil                      ! number of soil layers
 integer(i4b)                              :: nGlce                      ! number of glacier ice layers
 integer(i4b)                              :: nLayers                    ! total number of layers
 real(rkind),parameter                     :: xTol=1.e-10_rkind          ! small tolerance to address precision issues
 real(rkind),parameter                     :: areaTol=1.e-1_rkind        ! tolerance to address precision issues in glacier area summation
 real(rkind),parameter                     :: canIceTol=1.e-3_rkind      ! small tolerance to allow existence of canopy ice for above-freezing temperatures (kg m-2)
 real(rkind)                               :: remaining_area             ! remaining area to be distributed
 real(rkind)                               :: remaining_elev             ! remaining elevation to be distributed
 real(rkind)                               :: remaining_tan_slope        ! remaining tan slope to be distributed (area-weighted)
 real(rkind)                               :: remaining_aspect_sin       ! remaining sine component for circular aspect mean
 real(rkind)                               :: remaining_aspect_cos       ! remaining cosine component for circular aspect mean
 real(rkind),parameter                     :: deg2rad=PI_D/180._rkind    ! convert degrees to radians
 real(rkind),parameter                     :: rad2deg=180._rkind/PI_D    ! convert radians to degrees
 real(rkind),parameter                     :: aspect_tol=1.e-12_rkind    ! tolerance for undefined circular mean
 real(rkind)                               :: glacierAblAreaTot          ! total basin glacier ablation area from bvarData (m2)
 real(rkind)                               :: glacierAccAreaTot          ! total basin glacier accumulation area from bvarData (m2)
 real(rkind)                               :: area                       ! glacier area for a single glacier (m2)
 real(rkind)                               :: ratio                      ! ratio of glacier area to basin area
 real(rkind)                               :: frz_scale_use              ! scaling parameter for the snow or glce freezing curve (K-1)
 real(rkind)                               :: maxVolIceContent_use       ! maximum volumetric ice content depending if snow or firn
 ! --------------------------------------------------------------------------------------------------------

 ! Start procedure here
 err=0; message="check_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! Check that the initial conditions do not conflict with parameters, structure, etc.
 ! --------------------------------------------------------------------------------------------------------

 ! check and correct domain area and elevation, and ensure that the area is positive, and make backwards compatible
 do iGRU = 1,nGRU
   glacierAblAreaTot = 0._rkind
   glacierAccAreaTot = 0._rkind
   do iHRU=1,gru_struc(iGRU)%hruCount
     ! update the HRU area and elevation
     remaining_area = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)
     remaining_elev = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)*attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%elevation)
     remaining_tan_slope = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)*attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%tan_slope)
     remaining_aspect_sin = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)*sin(attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%aspect)*deg2rad)
     remaining_aspect_cos = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)*cos(attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%aspect)*deg2rad)
     do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
       associate(typeDOM => gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                 DOMarea => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1), &
                 DOMelev => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1), &
                 DOMtan_slope => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMtan_slope)%dat(1), &
                 DOMaspect => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMaspect)%dat(1), &
                 DOMcontourLength => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMcontourLength)%dat(1), &
                 glacMass4AreaChange => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%glacMass4AreaChange)%dat(1), &
                 scalarGlceWE => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarGlceWE)%dat(1), &
                 glacierAblFrac => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarAblFrac)%dat(1) )
         if(no_ice_vars)then ! set glacier ice variables to zero if they do not exist in the initial conditions file
           glacMass4AreaChange = 0._rkind
           scalarGlceWE = 0._rkind
         endif
         if (typeDOM.ne.upland .and. DOMarea>0._rkind) then
           if(no_dom_vars)then; err=20; message=trim(message)//'problem with getting variable id, var= DOMarea or DOMelev'; return; endif
           remaining_area = remaining_area - DOMarea
           remaining_elev = remaining_elev - DOMarea * DOMelev
           remaining_tan_slope = remaining_tan_slope - DOMarea * DOMtan_slope
           remaining_aspect_sin = remaining_aspect_sin - DOMarea*sin(DOMaspect*deg2rad)
           remaining_aspect_cos = remaining_aspect_cos - DOMarea*cos(DOMaspect*deg2rad)
         endif
         if (typeDOM==glacCln1 .or. typeDOM==glacCln2 .or. typeDOM==glacDbr) then
           if(no_ice_vars) write(*,'(A)') 'WARNING: glacier domain found but glacMass4AreaChange or scalarGlceWE missing, setting both to 0.0.'
           if(no_ablfrac)then; err=20; message=trim(message)//': problem with getting variable id, var= scalarAblFrac'; return; endif ! scalarAblFrac must exist if glacier domain
           glacierAblAreaTot = glacierAblAreaTot + glacierAblFrac*DOMarea
           glacierAccAreaTot = glacierAccAreaTot + (1.0_rkind-glacierAblFrac)*DOMarea
         else ! not glacier domain, ensure glacier variables are zero
           glacMass4AreaChange = 0._rkind
           scalarGlceWE = 0._rkind
           glacierAblFrac = 0._rkind
         end if
       end associate
     end do
     do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
       associate(typeDOM => gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                 DOMarea => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1), &
                 DOMelev => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1), &
                 DOMtan_slope => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMtan_slope)%dat(1), &
                 DOMaspect => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMaspect)%dat(1), &
                 DOMcontourLength => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMcontourLength)%dat(1) )
         if (typeDOM==upland) then
           DOMarea = remaining_area
           if(remaining_area>0._rkind)then
             DOMelev = remaining_elev/remaining_area
             DOMtan_slope = remaining_tan_slope/remaining_area
             if(remaining_aspect_sin**2 + remaining_aspect_cos**2 > aspect_tol)then
               DOMaspect = modulo(atan2(remaining_aspect_sin,remaining_aspect_cos)*rad2deg,360._rkind)
             else
               DOMaspect = 0._rkind
             endif
             DOMcontourLength = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%contourLength) ! for now, just set to the HRU contour length, but could be improved in the future
           else
             if (remaining_area<-xTol) write(*,'(A,E22.16,A)') 'WARNING: area of upland HRU (=', remaining_area, ') < 0. Resetting to 0.0'
             DOMelev = realMissing
             DOMarea = 0._rkind
             DOMtan_slope = realMissing
             DOMaspect = realMissing
             DOMcontourLength = 0._rkind
           end if
         end if
       end associate
     end do
   end do
   ! if they exist, check that the glacier areas are consistent with the basin areas, correct for grid tolerance (only warn if out of tolerance)
   if (gru_struc(iGRU)%nGlac>0) then
     if(sum(bvarData%gru(iGRU)%var(iLookBVAR%glacierAblArea)%dat)>0._rkind)then
       ratio = glacierAblAreaTot/sum(bvarData%gru(iGRU)%var(iLookBVAR%glacierAblArea)%dat)
       if ( abs( ratio - 1._rkind) >areaTol ) write(*,*) 'WARNING: glacier ablation domain area in GRU ',iGRU,' starts at ',ratio, &
       ' times the basin variable value but both will be calculated from the grid data after a year.'
     else ! = 0.0
        if (glacierAblAreaTot>0._rkind) then 
          write(*,'(A)') 'WARNING: basin glacier ablation area is zero, but there are glacier domains in the basin and scalarAblFrac>0. Resetting basin glacier ablation area to domain areas.'
          bvarData%gru(iGRU)%var(iLookBVAR%glacierAblArea)%dat(:) = glacierAblAreaTot/gru_struc(iGRU)%nGlac
        end if
     end if
     if(sum(bvarData%gru(iGRU)%var(iLookBVAR%glacierAccArea)%dat)>0._rkind )then
       ratio = glacierAccAreaTot/sum(bvarData%gru(iGRU)%var(iLookBVAR%glacierAccArea)%dat)
       if ( abs( ratio - 1._rkind) >areaTol ) write(*,*) 'WARNING: glacier accumulation domain area in GRU ',iGRU,' starts at ',ratio, &
       ' times the basin variable value but both will be calculated from the grid data after a year.'
     else ! = 0.0
       if (glacierAccAreaTot>0._rkind) then
         write(*,'(A)') 'WARNING: basin glacier accumulation area is zero, but there are glacier domains in the basin and scalarAblFrac<1. Resetting basin glacier accumulation area to domain areas.'
         bvarData%gru(iGRU)%var(iLookBVAR%glacierAccArea)%dat(:) = glacierAccAreaTot/gru_struc(iGRU)%nGlac
       end if
     end if
     ! check that the glacier areas are positive
     do iGlac = 1, gru_struc(iGRU)%nGlac
       area = bvarData%gru(iGRU)%var(iLookBVAR%glacierAblArea)%dat(iGlac) + bvarData%gru(iGRU)%var(iLookBVAR%glacierAccArea)%dat(iGlac)
       if (area<=0._rkind) write(*,*) 'WARNING: zero area for glacier ',iGlac,' in GRU ',iGRU,', need to set area to a postive value if want to build a glacier.'
     end do
   end if ! if glaciers

 enddo ! iGRU loop

 ! check for realistic values of albedo
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   do iDOM = 1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
    ! ensure the spectral average albedo is realistic
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) > mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMax)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMax)%dat(1)
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinWinter)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinWinter)%dat(1)
    ! ensure the visible albedo is realistic
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) > mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMaxVisible)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMaxVisible)%dat(1)
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) < mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinVisible)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinVisible)%dat(1)
    ! ensure the nearIR albedo is realistic
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) > mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMaxNearIR)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMaxNearIR)%dat(1)
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) < mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinNearIR)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinNearIR)%dat(1)
   end do
  end do
 end do

 ! ensure the initial conditions are consistent with the constitutive functions
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   do iDOM = 1,gru_struc(iGRU)%hruInfo(iHRU)%domCount

    ! associate local variables with variables in the data structures
    associate(&
    ! state variables in the canopy air space
    scalarCanairTemp     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanairTemp)%dat(1)     ,& ! canopy air temperature (K)
    scalarCanairEnthalpy => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanairEnthalpy)%dat(1) ,& ! canopy air enthalpy (J m-3)
    ! state variables in the vegetation canopy
    scalarCanopyTemp     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyTemp)%dat(1)     ,& ! canopy temperature (K)
    scalarCanopyEnthTemp => diagData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1) ,& ! canopy temperature component of enthalpy (J m-3)
    scalarCanopyEnthalpy => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyEnthalpy)%dat(1) ,& ! canopy enthalpy (J m-3)
    scalarCanopyLiq      => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyLiq)%dat(1)      ,& ! mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyIce      => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyIce)%dat(1)      ,& ! mass of ice on the vegetation canopy (kg m-2)
    heightCanopyTop      => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%heightCanopyTop)%dat(1)     ,& ! height of the top of the canopy layer (m)
    heightCanopyBottom   => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%heightCanopyBottom)%dat(1)  ,& ! height of the bottom of the canopy layer (m)
    specificHeatVeg      => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%specificHeatVeg)%dat(1)     ,& ! specific heat of vegetation (J kg-1 K-1)
    maxMassVegetation    => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%maxMassVegetation)%dat(1)   ,& ! maximum mass of vegetation (kg m-2)
    ! state variables in the layer domains
    mLayerTemp           => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerTemp)%dat              ,& ! temperature (K)
    mLayerEnthTemp       => diagData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookDIAG%mLayerEnthTemp)%dat          ,& ! temperature component of enthalpy (J m-3)
    mLayerEnthalpy       => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerEnthalpy)%dat          ,& ! enthalpy (J m-3)
    mLayerVolFracLiq     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracLiq)%dat        ,& ! volumetric fraction of liquid water in each snow layer (-)
    mLayerVolFracIce     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracIce)%dat        ,& ! volumetric fraction of ice in each snow layer (-)
    mLayerMatricHead     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerMatricHead)%dat        ,& ! matric head (m)
    layerType            => indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat              ,& ! type of layer (ix_soil or ix_snow)
    noThetaChange        => indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%noThetaChange)%dat(1)       ,& ! number of layers with no change in total water content (bottom layers)
    ! depth varying soil properties
    soil_dens_intr       => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%soil_dens_intr)%dat         ,& ! intrinsic soil density             (kg m-3)
    vGn_alpha            => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_alpha)%dat              ,& ! van Genutchen "alpha" parameter (m-1)
    vGn_n                => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_n)%dat                  ,& ! van Genutchen "n" parameter (-)
    theta_sat            => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_sat)%dat              ,& ! soil porosity (-)
    theta_res            => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_res)%dat              ,& ! soil residual volumetric water content (-)
    ! snow parameters
    snowfrz_scale        => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%snowfrz_scale)%dat(1)       ,& ! scaling parameter for the snow freezing curve (K-1)
    FCapil               => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%FCapil)%dat(1)               & ! fraction of pore space in tension storage (-)
    )  ! (associate local variables with model parameters)

     ! compute the constant in the freezing curve function (m K-1)
     kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2

     ! check canopy ice content for unrealistic situations
     if(scalarCanopyIce > canIceTol .and. scalarCanopyTemp > Tfreeze)then
      ! ice content > threshold, terminate run
      write(message,'(A,E22.16,A,E9.3,A,F7.3,A,F7.3,A)') trim(message)//'canopy ice (=',scalarCanopyIce,') > canIceTol (=',canIceTol,') when canopy temperature (=',scalarCanopyTemp,') > Tfreeze (=',Tfreeze,')'
      err=20; return
     else if(scalarCanopyIce > 0._rkind .and. scalarCanopyTemp > Tfreeze)then
      ! if here, ice content < threshold. Could be sublimation on previous timestep or simply wrong input. Print a warning
      write(*,'(A,E22.16,A,F7.3,A,F7.3,A)') 'WARNING: canopy ice content in restart file (=',scalarCanopyIce,') > 0 when canopy temperature (=',scalarCanopyTemp,') > Tfreeze (=',Tfreeze,'). Continuing.',NEW_LINE('a')
     end if
     scalarTheta = scalarCanopyIce + scalarCanopyLiq

     if(checkEnthalpy)then ! enthalpy as state variable or in residual
       if(no_icond_enth)then ! no enthalpy in icond file
         call T2enthTemp_cas(&
                    scalarCanairTemp,       & ! intent(in): canopy air temperature (K)
                    scalarCanairEnthalpy)     ! intent(out): enthalpy of the canopy air space (J m-3)
 
         call T2enthTemp_veg(&
                    (heightCanopyTop-heightCanopyBottom), & ! intent(in): canopy depth (m)
                    specificHeatVeg,        & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                    maxMassVegetation,      & ! intent(in): maximum mass of vegetation (kg m-2)
                    snowfrz_scale,          & ! intent(in): scaling parameter for the snow freezing curve  (K-1)
                    scalarCanopyTemp,       & ! intent(in): canopy temperature (K)
                    scalarTheta,            & ! intent(in): canopy water content (kg m-2)
                    scalarCanopyEnthTemp)     ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
         scalarCanopyEnthalpy = scalarCanopyEnthTemp - LH_fus * scalarCanopyIce/ (heightCanopyTop-heightCanopyBottom)
       else ! enthalpy is in the icond file
         scalarCanopyEnthTemp = scalarCanopyEnthalpy + LH_fus * scalarCanopyIce/ (heightCanopyTop-heightCanopyBottom)
       end if    
     end if

     ! number of layers
     nSnow   = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
     nLake   = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
     nSoil   = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
     nGlce   = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce
     nLayers = nSnow + nLake + nSoil + nGlce

     ! compute the maximum volumetric ice content for the layer domains
     if(nGlce>0)then ! snow can be firn
       maxVolIceContent_use = min(maxVolIceContent+0.15,0.85_rkind) ! firn maximum volumetric ice content to store water (-)
     else ! snow
       maxVolIceContent_use = maxVolIceContent ! snow maximum volumetric ice content to store water (-)
     endif

     ! loop through all layers
     do iLayer=1,nLayers

      ! *****
      ! * check that the initial volumetric fraction of liquid water and ice is reasonable...
      ! *************************************************************************************
      select case(layerType(iLayer))

       ! ***** snow, ice, lake, volume expansion allowed
       case(iname_snow, iname_lake, iname_glce)
        scalarTheta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)
        ! (check liquid water)
        if(mLayerVolFracLiq(iLayer) < 0._rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water < 0: layer = ',iLayer; err=20; return; end if
        if(mLayerVolFracLiq(iLayer) > 1._rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water > 1: layer = ',iLayer; err=20; return; end if
        ! (check ice)
        if (layerType(iLayer)==iname_snow) then
          if(mLayerVolFracIce(iLayer) > maxVolIceContent_use+0.1_rkind)then; write(message,'(a,f4.2,a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice > ',maxVolIceContent_use,': layer = ',iLayer; err=20; return; end if
          if(scalarTheta > maxVolIceContent_use+0.1_rkind)then; write(message,'(a,f4.2,a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] > ',maxVolIceContent_use,': layer = '    ,iLayer; err=20; return; end if
        else ! glacier ice or lake (could be all ice)
          if(mLayerVolFracIce(iLayer) > 1._rkind  )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice > 1: layer = ',iLayer; err=20; return; end if
          if(scalarTheta > 1._rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] > 1: layer = '      ,iLayer; err=20; return; end if
        end if
        if (layerType(iLayer)==iname_lake) then ! lake could be all liquid
          if(mLayerVolFracIce(iLayer) < 0._rkind  )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0: layer = '   ,iLayer; err=20; return; end if
        else if (layerType(iLayer)==iname_glce) then ! glacier ice should be mostly ice
          if(mLayerVolFracIce(iLayer) < 0.80_rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0.80: layer = ',iLayer; err=20; return; end if
        else if (layerType(iLayer)==iname_snow) then ! 
          if(mLayerVolFracIce(iLayer) < 0.05_rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0.05: layer = ',iLayer; err=20; return; end if
        end if
        ! check total water
        if(scalarTheta < 0.05_rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] < 0.05: layer = ',iLayer; err=20; return; end if

       ! ***** soil, no volume expansion
       case(iname_soil)
        iSoil       = iLayer - nSnow - nLake
        if(vGn_n(iSoil) <= 1._rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot have van Genutchen n <= 1: soil layer = ',iSoil; err=20; return; end if
        if(vGn_alpha(iSoil) >= 0._rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot have van Genutchen alpha >= 0: soil layer = ',iSoil; err=20; return; end if
        vGn_m       = 1._rkind - 1._rkind/vGn_n(iSoil)
        scalarTheta = mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)
        ! (check liquid water)
        if(mLayerVolFracLiq(iLayer) < theta_res(iSoil)-xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water < theta_res: layer = ',iLayer; err=20; return; end if
        if(mLayerVolFracLiq(iLayer) > theta_sat(iSoil)+xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water > theta_sat: layer = ',iLayer; err=20; return; end if
        ! (check ice)
        if(mLayerVolFracIce(iLayer) < 0._rkind             )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0: layer = '        ,iLayer; err=20; return; end if
        if(mLayerVolFracIce(iLayer) > theta_sat(iSoil)+xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice > theta_sat: layer = ',iLayer; err=20; return; end if
        ! check total water
        if(scalarTheta < theta_res(iSoil)-xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] < theta_res: layer = ',iLayer; err=20; return; end if
        if(scalarTheta > theta_sat(iSoil)+xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] > theta_sat: layer = ',iLayer; err=20; return; end if

       case default
        write(*,*) 'Cannot recognize case in initial vol water/ice check: type=', layerType(iLayer)
        err=20; message=trim(message)//'cannot identify layer type'; return
      end select

      ! *****
      ! * check that the initial conditions are consistent with the constitutive functions...
      ! *************************************************************************************
      select case(layerType(iLayer))

       ! ** snow, lake, ice
       case(iname_snow, iname_lake, iname_glce)
         
        if (layerType(iLayer)==iname_snow .or. layerType(iLayer)==iname_glce)then 
          ! check that temperature is less than freezing
          if(mLayerTemp(iLayer) > Tfreeze)then
           message=trim(message)//'initial snow or glacier ice temperature is greater than freezing'
           err=20; return
          end if
        endif
        if (layerType(iLayer)==iname_snow) frz_scale_use = snowfrz_scale
        if (layerType(iLayer)==iname_glce .or. layerType(iLayer)==iname_lake) frz_scale_use = snowfrz_scale*icefrz_mult
        
        ! ensure consistency among state variables
        call updatSnLaGl(&
                        iLayer>nLayers-noThetaChange,   & ! intent(in):  flag that no liquid water in layer
                        mLayerTemp(iLayer),             & ! intent(in):  temperature (K)
                        scalarTheta,                    & ! intent(in):  volumetric fraction of total water (-)
                        frz_scale_use,                  & ! intent(in):  scaling parameter for the freezing curve (K-1)
                        mLayerVolFracLiq(iLayer),       & ! intent(out): volumetric fraction of liquid water (-)
                        mLayerVolFracIce(iLayer),       & ! intent(out): volumetric fraction of ice (-)
                        fLiq,                           & ! intent(out): fraction of liquid water (-)
                        err,cmessage)                     ! intent(out): error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

        if(checkEnthalpy)then ! enthalpy as state variable or in residual
          if(no_icond_enth)then ! no enthalpy in icond file
            call T2enthTemp_snLaGl(&
                        iLayer>nLayers-noThetaChange,   & ! intent(in):  flag that no liquid water in layer
                        frz_scale_use,                  & ! intent(in):  scaling parameter for the freezing curve  (K-1)
                        mLayerTemp(iLayer),             & ! intent(in):  layer temperature (K)
                        scalarTheta,                    & ! intent(in):  volumetric total water content (-)
                        mLayerEnthTemp(iLayer))           ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
            mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_ice * LH_fus * mLayerVolFracIce(iLayer)
          else
            mLayerEnthTemp(iLayer) = mLayerEnthalpy(iLayer) + iden_ice * LH_fus * mLayerVolFracIce(iLayer)
          end if
        endif

       ! ** soil
       case(iname_soil)

       ! ensure consistency among state variables
       call updatSoil(&
                      mLayerTemp(iLayer),              & ! intent(in): layer temperature (K)
                      mLayerMatricHead(iSoil),         & ! intent(in): matric head (m)
                      vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m, & ! intent(in): van Genutchen soil parameters
                      scalarTheta,                     & ! intent(out): volumetric fraction of total water (-)
                      mLayerVolFracLiq(iLayer),        & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer),        & ! intent(out): volumetric fraction of ice (-)
                      err,cmessage)                      ! intent(out): error control
       if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

       if(checkEnthalpy)then ! enthalpy as state variable or in residual
         if(no_icond_enth)then ! no enthalpy in icond file
           call T2enthTemp_soil(&
                       use_lookup,                      & ! intent(in):  flag to use the lookup table for soil enthalpy
                       soil_dens_intr(iSoil),           & ! intent(in):  intrinsic soil density (kg m-3)
                       vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m, & ! intent(in): van Genutchen soil parameters
                       iSoil,                           & ! intent(in):  index of the control volume within the domain
                       lookupData%gru(iGRU)%hru(iHRU)%dom(iDOM),  & ! intent(in):  lookup table data structure
                       realMissing,                     & ! intent(in):  lower value of integral (not computed)
                       mLayerTemp(iLayer),              & ! intent(in):  layer temperature (K)
                       mLayerMatricHead(iSoil),         & ! intent(in):  matric head (m)
                       mLayerEnthTemp(iLayer),          & ! intent(out): temperature component of enthalpy soil layer (J m-3)
                       err,cmessage)                      ! intent(out): error control
           if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
           mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_water * LH_fus * mLayerVolFracIce(iLayer)
         else
           mLayerEnthTemp(iLayer) = mLayerEnthalpy(iLayer) + iden_water * LH_fus * mLayerVolFracIce(iLayer)
         end if
       endif

       case default; err=10; message=trim(message)//'unknown case for model layer'; return
      end select

     end do  ! (looping through layers)
    
    ! end association to variables in the data structures
    end associate

    ! if snow layers exist, compute snow depth and SWE
    if(nSnow>0)then
     progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSWE)%dat(1) = sum( (progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                                                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice)  * &
                                                                                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
    end if  ! if snow layers exist

    ! check that the layering is consistent
    do iLayer=1,nLayers
     h1 = sum(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerDepth)%dat(1:iLayer)) ! sum of the depths up to the current layer
     h2 = progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%iLayerHeight)%dat(iLayer) - progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%iLayerHeight)%dat(0)  ! difference between snow-atm interface and bottom of layer
     if(abs(h1 - h2) > verySmall)then
      write(message,'(a,1x,i0,a,f6.3,a,f6.3)') trim(message)//'mis-match between layer depth and layer height; layer = ', iLayer, '; sum depths = ',h1,'; height = ',h2
      err=20; return
     end if
    end do

   end do ! iDOM
  end do ! iHRU
 end do ! iGRU

 end subroutine check_icond

end module check_icond_module
