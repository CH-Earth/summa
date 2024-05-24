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

module enthalpyTemp_module

! constants
USE multiconst, only: gravity, &                          ! gravitational acceleration (m s-1)
                      Tfreeze, &                          ! freezing point of water (K)
                      Cp_soil,Cp_water,Cp_ice,Cp_air,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                      iden_water,iden_ice,iden_air,&      ! intrinsic density of water and ice (kg m-3)
                      LH_fus                              ! latent heat of fusion (J kg-1)

! data types
USE nrtype
USE data_types,only:var_iLength                    ! var(:)%dat(:)
USE data_types,only:var_dLength                    ! var(:)%dat(:)
USE data_types,only:zLookup                        ! z(:)%var(:)%lookup(:)

! indices within parameter structure
USE var_lookup,only:iLookPARAM                     ! named variables to define structure element
USE var_lookup,only:iLookINDEX                     ! named variables to define structure element
USE var_lookup,only:iLookLOOKUP                    ! named variables to define structure element
USE var_lookup,only:iLookDIAG                      ! named variables for structure elements

! data dimensions
USE var_lookup,only:maxvarLookup                   ! maximum number of variables in the lookup tables

! domain types
USE globalData,only:iname_cas                      ! named variables for canopy air space
USE globalData,only:iname_veg                      ! named variables for vegetation canopy
USE globalData,only:iname_snow                     ! named variables for snow
USE globalData,only:iname_soil                     ! named variables for soil
USE globalData,only:iname_aquifer                  ! named variables for the aquifer

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair                ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy                ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_nrgLayer                 ! named variable defining the energy state variable for snow+soil layers

! missing values
USE globalData,only:integerMissing                 ! missing integer
USE globalData,only:realMissing                    ! missing real number

implicit none
public::T2H_lookup_snWat
public::T2L_lookup_soil
public::enthalpy2T_snwWat
public::T2enthalpy_snwWat
public::T2enthTemp_cas
public::T2enthTemp_veg
public::T2enthTemp_snow
public::T2enthTemp_soil
public::enthTemp_or_enthalpy
public::enthalpy2T_cas
public::enthalpy2T_veg
public::enthalpy2T_snow
public::enthalpy2T_soil
private::hyp_2F1_real
private::brent, brent0, diff_H_veg, diff_H_snow, diff_H_soil

! define the snow look-up table used to compute temperature based on enthalpy
integer(i4b),parameter               :: nlook=10001       ! number of elements in the lookup table
real(rkind),dimension(nlook),public  :: H_lookup          ! enthalpy values (J kg-1)
real(rkind),dimension(nlook),public  :: T_lookup          ! temperature values (K)
contains


! ************************************************************************************************************************
! public subroutine T2H_lookup_snWat: define a look-up table to liquid + ice enthalpy based on temperature
!                                     appropriate when no dry mass, as in snow
! ************************************************************************************************************************
subroutine T2H_lookup_snWat(mpar_data,                     &  ! intent(in):    parameter data structure
                           err,message)
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines 
  USE nr_utility_module,only:arth                       ! use to build vectors with regular increments
  USE spline_int_module,only:spline,splint              ! use for cubic spline interpolation
  implicit none
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare dummy variables
  type(var_dlength),intent(in)  :: mpar_data            ! model parameters
  integer(i4b),intent(out)      :: err                  ! error code
  character(*),intent(out)      :: message              ! error message
  ! declare local variables
  character(len=128)            :: cmessage             ! error message in downwind routine
  real(rkind),parameter         :: T_start=260.0_rkind  ! start temperature value where all liquid water is assumed frozen (K)
  real(rkind)                   :: T_incr,H_incr        ! temperature/enthalpy increments
  real(rkind),dimension(nlook)  :: Tk                   ! initial temperature vector
  real(rkind),dimension(nlook)  :: Hy                   ! initial enthalpy vector
  real(rkind),parameter         :: waterWght=1._rkind   ! weight applied to total water (kg m-3) --- cancels out
  real(rkind),dimension(nlook)  :: H2                   ! 2nd derivatives of the interpolating function at tabulated points
  real(rkind)                   :: dT                   ! derivative of temperature with enthalpy at H_lookup
  integer(i4b)                  :: ilook                ! loop through lookup table
  ! -------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="T2H_lookup_snWat/"

  ! associate
  associate( snowfrz_scale => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1) )

    ! define initial temperature vector
    T_incr = (Tfreeze - T_start) / real(nlook-1, kind(rkind))  ! temperature increment
    Tk     = arth(T_start,T_incr,nlook)
    ! ***** compute specific enthalpy (NOTE: J m-3 --> J kg-1) *****

    do ilook=1,nlook
      Hy(ilook) = T2enthalpy_snwWat(Tk(ilook),waterWght,snowfrz_scale)/waterWght  ! (J m-3 --> J kg-1)
    end do

    ! define the final enthalpy vector
    H_incr   = (-Hy(1)) / real(nlook-1, kind(rkind))  ! enthalpy increment
    H_lookup = arth(Hy(1),H_incr,nlook)

    ! use cubic spline interpolation to obtain temperature values at the desired values of enthalpy
    call spline(Hy,Tk,1.e30_rkind,1.e30_rkind,H2,err,cmessage)  ! get the second derivatives
    if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

    do ilook=1,nlook
      call splint(Hy,Tk,H2,H_lookup(ilook),T_lookup(ilook),dT,err,cmessage)
      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
    end do

  end associate

 end subroutine T2H_lookup_snWat

! ************************************************************************************************************************
! public subroutine T2L_lookup_soil: define a look-up table to compute integral of soil Clapeyron equation liquid water
!                                    matric potential from temperature
! ************************************************************************************************************************
subroutine T2L_lookup_soil(nSoil,                         &  ! intent(in):    number of soil layers
                           mpar_data,                     &  ! intent(in):    parameter data structure
                           lookup_data,                   &  ! intent(inout): lookup table data structure
                           err,message)
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines                    
  USE nr_utility_module,only:arth                       ! use to build vectors with regular increments
  USE spline_int_module,only:spline,splint              ! use for cubic spline interpolation
  USE soil_utils_module,only:volFracLiq                 ! use to compute the volumetric fraction of liquid water
  implicit none
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare dummy variables
  integer(i4b),intent(in)       :: nSoil
  type(var_dlength),intent(in)  :: mpar_data            ! model parameters
  type(zLookup),intent(inout)   :: lookup_data          ! lookup tables
  integer(i4b),intent(out)      :: err                  ! error code
  character(*),intent(out)      :: message              ! error message
  ! declare local variables
  character(len=128)            :: cmessage             ! error message in downwind routine
  integer(i4b),parameter        :: nLook=500            ! number of elements in the lookup table
  integer(i4b),parameter        :: nIntegr8=10000       ! number of points used in the numerical integration
  real(rkind),parameter         :: T_lower=260.0_rkind  ! lowest temperature value where all liquid water is assumed frozen (K)
  real(rkind),dimension(nLook)  :: xTemp                ! temporary vector
  real(rkind)                   :: xIncr                ! temporary increment
  real(rkind)                   :: T_incr               ! temperature increment
  real(rkind)                   :: dL                   ! derivative of integral with temperature at T_test
  integer(i4b)                  :: iVar                 ! loop through variables
  integer(i4b)                  :: iSoil                ! loop through soil layers
  integer(i4b)                  :: iLook                ! loop through lookup table
  integer(i4b)                  :: jIntegr8             ! index for numerical integration
  logical(lgt)                  :: check                ! flag to check allocation
  real(rkind)                   :: vGn_m                ! van Genuchten "m" parameter (-)
  real(rkind)                   :: vFracLiq             ! volumetric fraction of liquid water (-)
  real(rkind)                   :: matricHead           ! matric head (m)
  ! -------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="T2L_lookup_soil/"

  ! get the values of temperature for the lookup table
  xIncr = 1._rkind/real(nLook-1, kind(rkind))
  xTemp = T_lower + (Tfreeze - T_lower)*sqrt(sqrt(arth(0._rkind,xIncr,nLook))) ! use sqrt(sqrt()) to give more values near freezing

  ! -----
  ! * allocate space for the lookup table...
  ! ----------------------------------------

  ! initialize checks
  check=.false.

  ! allocate space for soil layers
  if(allocated(lookup_data%z))then; check=.true.; else; allocate(lookup_data%z(nSoil), stat=err); endif
  if(check) then; err=20; message=trim(message)//'lookup table z dimension was unexpectedly allocated already'; return; end if
  if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table z dimension dimension'; return; end if

  ! allocate space for the variables in the lookup table
  do iSoil=1,nSoil
    if(allocated(lookup_data%z(iSoil)%var))then; check=.true.; else; allocate(lookup_data%z(iSoil)%var(maxvarLookup), stat=err); endif
    if(check) then; err=20; message=trim(message)//'lookup table var dimension was unexpectedly allocated already'; return; end if
    if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table var dimension dimension'; return; end if

    ! allocate space for the values in the lookup table
    do iVar=1,maxvarLookup
      if(allocated(lookup_data%z(iSoil)%var(iVar)%lookup))then; check=.true.; else; allocate(lookup_data%z(iSoil)%var(iVar)%lookup(nLook), stat=err); endif
      if(check) then; err=20; message=trim(message)//'lookup table value dimension was unexpectedly allocated already'; return; end if
      if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table vaule dimension dimension'; return; end if

    end do ! (looping through variables)
  end do ! (looping through soil layers)

  ! loop through soil layers
  do iSoil=1,nSoil

    ! -----
    ! * make association to variables in the data structures...
    ! ---------------------------------------------------------

    associate(&

      ! associate model parameters
      snowfrz_scale  => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)           , & ! scaling parameter for freezing     (K-1)
      soil_dens_intr => mpar_data%var(iLookPARAM%soil_dens_intr)%dat(iSoil)      , & ! intrinsic soil density             (kg m-3)
      theta_sat      => mpar_data%var(iLookPARAM%theta_sat)%dat(iSoil)           , & ! soil porosity                      (-)
      theta_res      => mpar_data%var(iLookPARAM%theta_res)%dat(iSoil)           , & ! volumetric residual water content  (-)
      vGn_alpha      => mpar_data%var(iLookPARAM%vGn_alpha)%dat(iSoil)           , & ! van Genuchten "alpha" parameter    (m-1)
      vGn_n          => mpar_data%var(iLookPARAM%vGn_n)%dat(iSoil)               , & ! van Genuchten "n" parameter        (-)

      ! associate values in the lookup table
      Tk            => lookup_data%z(iSoil)%var(iLookLOOKUP%temperature)%lookup  , & ! temperature (K)
      Ly            => lookup_data%z(iSoil)%var(iLookLOOKUP%psiLiq_int)%lookup   , & ! integral of mLayerPsiLiq from Tfreeze to Tk (K)
      L2            => lookup_data%z(iSoil)%var(iLookLOOKUP%deriv2)%lookup         & ! second derivative of the interpolating function

      ) ! end associate statement

      ! compute vGn_m
      vGn_m = 1._rkind - 1._rkind/vGn_n

      ! -----
      ! * populate the lookup table...
      ! ------------------------------

      ! initialize temperature and integral
      Tk(nLook) = Tfreeze
      Ly(nLook) = 0._rkind

      ! loop through lookup table
      do iLook=(nLook-1),1,-1

        ! update temperature and integral
        Tk(iLook) = Tk(iLook+1)
        Ly(iLook) = Ly(iLook+1)

        ! get the temperature increment for the numerical integration
        T_incr = (xTemp(iLook)-xTemp(iLook+1))/real(nIntegr8, kind(rkind))

        ! numerical integration between different values of the lookup table
        do jIntegr8=1,nIntegr8

          ! update temperature
          Tk(iLook)  = Tk(iLook) + T_incr

          ! compute the volumetric liquid water and ice content at the mid point of the temperature increment
          matricHead = (LH_fus/gravity)*(Tk(iLook) - Tfreeze - T_incr/2._rkind)/Tfreeze
          vFracLiq   = volFracLiq(matricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)

          ! compute integral
          Ly(iLook)  = Ly(iLook) + vFracLiq*T_incr
  
        end do  ! numerical integration

      end do  ! loop through lookup table

      ! use cubic spline interpolation to obtain integral values at the desired values of temperature
      call spline(Tk,Ly,1.e30_rkind,1.e30_rkind,L2,err,cmessage)  ! get the second derivatives
      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

    ! end asssociation to variables in the data structures
    end associate

  end do  ! (looping through soil layers)
end subroutine T2L_lookup_soil


! ************************************************************************************************************************
! public subroutine enthalpy2T_snwWat: compute temperature based on specific temperature component of liquid + ice enthalpy 
!                                      appropriate when no dry mass, as in snow. Uses look-up table for enthalpy
! ************************************************************************************************************************
subroutine enthalpy2T_snwWat(Hy,BulkDenWater,fc_param,Tk,err,message)
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare dummy variables
  real(rkind),intent(in)      :: Hy            ! total enthalpy (J m-3)
  real(rkind),intent(in)      :: BulkDenWater  ! bulk density of water (kg m-3)
  real(rkind),intent(in)      :: fc_param      ! freezing curve parameter (K-1)
  real(rkind),intent(out)     :: Tk            ! initial temperature guess / final temperature value (K)
  integer(i4b),intent(out)    :: err           ! error code
  character(*),intent(out)    :: message       ! error message
  ! declare local variables
  real(rkind),parameter       :: dx=1.d-8      ! finite difference increment (J kg-1)
  real(rkind),parameter       :: atol=1.d-12   ! convergence criteria (J kg-1)
  real(rkind)                 :: H_spec        ! specific enthalpy (J kg-1)
  real(rkind)                 :: H_incr        ! enthalpy increment
  integer(i4b)                :: niter=15      ! maximum number of iterations
  integer(i4b)                :: iter          ! iteration index
  integer(i4b)                :: i0            ! position in lookup table
  real(rkind)                 :: Tg0,Tg1       ! trial temperatures (K)
  real(rkind)                 :: Ht0,Ht1       ! specific enthalpy, based on the trial temperatures (J kg-1)
  real(rkind)                 :: f0,f1         ! function evaluations (difference between enthalpy guesses)
  real(rkind)                 :: dh            ! enthalpy derivative
  real(rkind)                 :: dT            ! temperature increment
  ! -------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="enthalpy2T_snwWat/"
  ! convert input of total enthalpy (J m-3) to total specific enthalpy (J kg-1)
  H_spec = Hy/BulkDenWater ! (NOTE: no soil)
 
  ! ***** get initial guess and derivative assuming all water is frozen
  if(H_spec<H_lookup(1))then ! process cases below the limit of the look-up table
    ! get temperature guess
    Tg0 = (H_spec - H_lookup(1))/Cp_ice + T_lookup(1)
    Tg1 = Tg0+dx
    ! compute enthalpy
    Ht0 = T2enthalpy_snwWat(Tg0,1._rkind,fc_param)
    Ht1 = T2enthalpy_snwWat(Tg1,1._rkind,fc_param)
    ! compute function evaluations
    f0  = Ht0 - H_spec
    f1  = Ht1 - H_spec

  ! ***** get initial guess and derivative from the look-up table
  else
    ! get enthalpy increment
    H_incr = H_lookup(2) - H_lookup(1)
    ! get position in lookup table
    i0 = ceiling( (H_spec - H_lookup(1)) / H_incr, kind(i4b) )
    ! check found the appropriate value in the look-up table
    if(H_spec < H_lookup(i0) .or. H_spec > H_lookup(i0+1) .or. &
       i0 < 1 .or. i0+1 > nlook)then
     err=10; message=trim(message)//'problem finding appropriate value in lookup table'; return
    end if
    ! get temperature guess
    Tg0 = T_lookup(i0)
    Tg1 = T_lookup(i0+1)
    ! compute function evaluations
    f0  = H_lookup(i0) - H_spec
    f1  = H_lookup(i0+1) - H_spec
  end if

  ! compute initial derivative
  dh  = (f1 - f0) / (Tg1 - Tg0)
  ! compute initial change in T
  dT  = -f0/dh
  ! exit if already close enough
  if(abs(dT)<atol)then
    Tk = Tg0+dT
    return
  end if

  ! **** iterate a little
  do iter=1,niter
    ! comute new value of Tg
    Tg1 = Tg0+dT
    ! get new function evaluation
    Ht1 = T2enthalpy_snwWat(Tg1,1._rkind,fc_param)
    f1  = Ht1 - H_spec
    ! compute derivative of dT
    dh  = (f1 - f0)/dT
    ! compute change in T
    dT  = -f1/dh
    ! exit if converged
    if(abs(dT)<atol)then
      Tk = Tg1+dT
      return
    end if
    ! get ready for next iteration -- save old function evaluation and temperature
    f0  = f1
    Tg0 = Tg1
    ! and check for convergence
    if(iter==niter)then; err=20; message=trim(message)//"failedToConverge"; return; end if
  end do  ! (iteration loop)
end subroutine enthalpy2T_snwWat


! ************************************************************************************************************************
! public function T2enthalpy_snwWat: compute liquid and ice mixture enthalpy based on temperature and mass (J m-3) for a
!                                    layer only where the layer has no dry mass, as in snow.
!                                    NOTE: enthalpy is a relative value, defined as zero at Tfreeze where all water is liquid
! ************************************************************************************************************************
function T2enthalpy_snwWat(Tk,BulkDenWater,fc_param)
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! declare dummy variables
  real(rkind),intent(in)  :: Tk                ! layer temperature (K)
  real(rkind),intent(in)  :: BulkDenWater      ! bulk density of water (kg m-3)
  real(rkind),intent(in)  :: fc_param          ! freezing curve parameter (K-1)
  real(rkind)             :: T2enthalpy_snwWat ! return value of the function, total specific enthalpy (J m-3)
  ! declare local variables
  real(rkind)             :: frac_liq          ! fraction of liquid water
  real(rkind)             :: enthTempWater     ! temperature component of specific enthalpy for total water (liquid and ice) (J kg-1)
  real(rkind)             :: enthMass          ! mass component of specific enthalpy (J kg-1)
  ! -------------------------------------------------------------------------------------------------------------------------
  ! compute the fraction of liquid water in the given layer
  frac_liq     = 1._rkind / ( 1._rkind + ( fc_param*( Tfreeze - min(Tk,Tfreeze) ) )**2_i4b )

  ! compute the temperature component of enthalpy for total water (J kg-1)
  ! NOTE: negative enthalpy means require energy to bring to Tfreeze
  if(Tk< Tfreeze) enthTempWater = Cp_ice*(Tk - Tfreeze) - (Cp_water - Cp_ice)*(atan(fc_param*(Tfreeze - Tk))/fc_param)
  if(Tk>=Tfreeze) enthTempWater = Cp_water*(Tk - Tfreeze)

  ! compute the mass component of enthalpy -- energy required to melt ice (J kg-1)
  ! NOTE: negative enthalpy means require energy to bring to Tfreeze
  enthMass = -LH_fus*(1._rkind - frac_liq)

  ! finally, compute the total enthalpy (J m-3)
  T2enthalpy_snwWat = BulkDenWater*(enthTempWater + enthMass) !+ BulkDenSoil*enthTempSoil
end function T2enthalpy_snwWat


! ************************************************************************************************************************
! public subroutine T2enthTemp_cas: compute temperature component of enthalpy from temperature and total water content, canopy air space
! ************************************************************************************************************************
subroutine T2enthTemp_cas(&
                      scalarCanairTemp,       & ! intent(in):  canopy air temperature (K)
                      scalarCanairEnthalpy,   & ! intent(out): enthalpy of the canopy air space (J m-3)
                      err,message)              ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  ! input: variables for the canopy air space
  real(rkind),intent(in)           :: scalarCanairTemp      ! canopy air temperature (K)
  ! output: enthalpy
  real(rkind),intent(out)          :: scalarCanairEnthalpy  ! enthalpy of the canopy air space (J m-3)
  ! output: error control
  integer(i4b),intent(out)         :: err                   ! error code
  character(*),intent(out)         :: message               ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="T2enthTemp_cas/"

  scalarCanairEnthalpy = Cp_air * iden_air * (scalarCanairTemp - Tfreeze)

end subroutine T2enthTemp_cas

! ************************************************************************************************************************
! public subroutine T2enthTemp_veg: compute temperature component of enthalpy from temperature and total water content, canopy
! ************************************************************************************************************************
subroutine T2enthTemp_veg(&
                      canopyDepth,            & ! intent(in):  canopy depth (m)
                      specificHeatVeg,        & ! intent(in):  specific heat of vegetation (J kg-1 K-1)
                      maxMassVegetation,      & ! intent(in):  maximum mass of vegetation (kg m-2)
                      snowfrz_scale,          & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                      scalarCanopyTemp,       & ! intent(in):  canopy temperature (K)
                      scalarCanopyWat,        & ! intent(in):  canopy total water (kg m-2)
                      scalarCanopyEnthTemp,   & ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
                      err,message)              ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  real(rkind),intent(in)           :: canopyDepth           ! canopy depth (m)
  real(rkind),intent(in)           :: specificHeatVeg       ! specific heat of vegetation (J kg-1 K-1)
  real(rkind),intent(in)           :: maxMassVegetation     ! maximum mass of vegetation (kg m-2)
  real(rkind),intent(in)           :: snowfrz_scale         ! scaling parameter for the snow freezing curve  (K-1)
  ! input: variables for the vegetation canopy
  real(rkind),intent(in)           :: scalarCanopyTemp      ! canopy temperature (K)
  real(rkind),intent(in)           :: scalarCanopyWat       ! canopy total water (kg m-2)
  ! output: enthalpy
  real(rkind),intent(out)          :: scalarCanopyEnthTemp  ! temperature component of enthalpy of the vegetation canopy (J m-3)
  ! output: error control
  integer(i4b),intent(out)         :: err                   ! error code
  character(*),intent(out)         :: message               ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  real(rkind)                      :: diffT                 ! temperature difference of temp from Tfreeze
  real(rkind)                      :: integral              ! integral of snow freezing curve
  ! enthalpy
  real(rkind)                      :: enthVeg               ! enthalpy of the vegetation (J m-3)
  real(rkind)                      :: enthLiq               ! enthalpy of the liquid region (J m-3)
  real(rkind)                      :: enthIce               ! enthalpy of the ice region (J m-3)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="T2enthTemp_veg/"

  diffT   = scalarCanopyTemp - Tfreeze
  enthVeg = specificHeatVeg * maxMassVegetation * diffT / canopyDepth

  if(diffT>=0._rkind)then
    enthLiq = Cp_water * scalarCanopyWat * diffT / canopyDepth
    enthIce = 0._rkind
  else
    integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
    enthLiq  = Cp_water * scalarCanopyWat * integral / canopyDepth
    enthIce  = Cp_ice * scalarCanopyWat * ( diffT - integral ) / canopyDepth
  endif

  scalarCanopyEnthTemp = enthVeg + enthLiq + enthIce

end subroutine T2enthTemp_veg

! ************************************************************************************************************************
! public subroutine T2enthTemp_snow: compute temperature component of enthalpy from temperature and total water content, snow layer
! ************************************************************************************************************************
subroutine T2enthTemp_snow(&
                      snowfrz_scale,          & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                      mLayerTemp,             & ! intent(in):  layer temperature (K)
                      mLayerVolFracWat,       & ! intent(in):  volumetric total water content (-)
                      mLayerEnthTemp,         & ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
                      err,message)              ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  real(rkind),intent(in)           :: snowfrz_scale         ! scaling parameter for the snow freezing curve  (K-1)
  ! input: variables for the snow domain
  real(rkind),intent(in)           :: mLayerTemp            ! layer temperature (K)
  real(rkind),intent(in)           :: mLayerVolFracWat      ! volumetric total water content (-)
  ! output: enthalpy
  real(rkind),intent(out)          :: mLayerEnthTemp        ! temperature component of enthalpy of each snow layer (J m-3)
  ! output: error control
  integer(i4b),intent(out)         :: err                   ! error code
  character(*),intent(out)         :: message               ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  real(rkind)                      :: diffT                 ! temperature difference of temp from Tfreeze
  real(rkind)                      :: integral              ! integral of snow freezing curve
  ! enthalpy
  real(rkind)                      :: enthLiq               ! enthalpy of the liquid region (J m-3)
  real(rkind)                      :: enthIce               ! enthalpy of the ice region (J m-3)
  real(rkind)                      :: enthAir               ! enthalpy of air (J m-3)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="T2enthTemp_snow/"

  diffT    = mLayerTemp - Tfreeze  ! diffT<0._rkind because snow is frozen

  if(diffT==0._rkind)then ! only need for upper bound
    enthLiq = 0._rkind
    enthIce = 0._rkind
    enthAir = 0._rkind
  else
    integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
    enthLiq  = iden_water * Cp_water * mLayerVolFracWat * integral
    enthIce  = iden_water * Cp_ice * mLayerVolFracWat * ( diffT - integral )
    enthAir  = iden_air * Cp_air * ( diffT - mLayerVolFracWat * ( (iden_water/iden_ice)*(diffT-integral) + integral ) )
  endif

  mLayerEnthTemp = enthLiq + enthIce + enthAir

end subroutine T2enthTemp_snow


! ************************************************************************************************************************
! public subroutine T2enthTemp_soil: compute temperature component of enthalpy from temperature and total water content, soil layer
! ************************************************************************************************************************
subroutine T2enthTemp_soil(&
                      use_lookup,               & ! intent(in):  flag to use the lookup table for soil enthalpy
                      soil_dens_intr,           & ! intent(in):  intrinsic soil density (kg m-3)
                      vGn_alpha,                & ! intent(in):  van Genutchen "alpha" parameter
                      vGn_n,                    & ! intent(in):  van Genutchen "n" parameter
                      theta_sat,                & ! intent(in):  soil porosity (-)
                      theta_res,                & ! intent(in):  soil residual volumetric water content (-)
                      vGn_m,                    & ! intent(in):  van Genutchen "m" parameter (-)
                      ixControlIndex,           & ! intent(in):  index of the control volume within the domain
                      lookup_data,              & ! intent(in):  lookup table data structure
                      integral_frz_low0,        & ! intent(in):  integral_frz_low if computed outside, else realMissing
                      mLayerTemp,               & ! intent(in):  layer temperature (K)
                      mLayerMatricHead,         & ! intent(in):  total water matric potential (m)
                      mLayerEnthTemp,           & ! intent(out): temperature component of enthalpy soil layer (J m-3)
                      err,message)                ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
  USE spline_int_module,only:splint         ! use for cubic spline interpolation
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  logical(lgt),intent(in)          :: use_lookup             ! flag to use the lookup table for soil enthalpy, otherwise use hypergeometric function
 ! input: data structures
  real(rkind),intent(in)           :: soil_dens_intr         ! intrinsic soil density (kg m-3)
  real(rkind),intent(in)           :: vGn_alpha              ! van Genutchen "alpha" parameter
  real(rkind),intent(in)           :: vGn_n                  ! van Genutchen "n" parameter
  real(rkind),intent(in)           :: theta_sat              ! soil porosity (-)
  real(rkind),intent(in)           :: theta_res              ! soil residual volumetric water content (-)
  real(rkind),intent(in)           :: vGn_m                  ! van Genutchen "m" parameter (-)
  integer(i4b),intent(in)          :: ixControlIndex         ! index within a given model domain
  type(zLookup),intent(in)         :: lookup_data            ! lookup tables
  real(rkind),intent(in)           :: integral_frz_low0      ! integral_frz_low if computed outside, else realMissing
  ! input: variables for the soil domain
  real(rkind),intent(in)           :: mLayerTemp             ! layer temperature (K)
  real(rkind),intent(in)           :: mLayerMatricHead       ! total water matric potential (m)
  ! output: enthalpy
  real(rkind),intent(out)          :: mLayerEnthTemp         ! temperature component of enthalpy of soil layer (J m-3)
  ! output: error control
  integer(i4b),intent(out)         :: err                    ! error code
  character(*),intent(out)         :: message                ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  character(len=128)               :: cmessage               ! error message in downwind routine
  real(rkind)                      :: Tcrit                  ! temperature where all water is unfrozen (K)
  real(rkind)                      :: volFracWat             ! volumetric fraction of total water, liquid+ice (-)
  real(rkind)                      :: diff0                  ! temperature difference of Tcrit from Tfreeze
  real(rkind)                      :: diffT                  ! temperature difference of temp from Tfreeze
  real(rkind)                      :: dL                     ! derivative of soil lookup table with temperature at layer temperature
  real(rkind)                      :: arg                    ! argument of soil hypergeometric function
  real(rkind)                      :: gauss_hg_T             ! soil hypergeometric function result
  real(rkind)                      :: integral_unf           ! integral of unfrozen soil water content (from Tfreeze to Tcrit)
  real(rkind)                      :: integral_frz_low       ! lower limit of integral of frozen soil water content (from Tfreeze to Tcrit)
  real(rkind)                      :: integral_frz_upp       ! upper limit of integral of frozen soil water content (from Tfreeze to soil temperature)
  real(rkind)                      :: xConst                 ! constant in the freezing curve function (m K-1)
  real(rkind)                      :: mLayerPsiLiq           ! liquid water matric potential (m)
  ! enthalpy
  real(rkind)                      :: enthSoil               ! enthalpy of soil particles (J m-3)
  real(rkind)                      :: enthLiq                ! enthalpy of the liquid region (J m-3)
  real(rkind)                      :: enthIce                ! enthalpy of the ice region (J m-3)
  real(rkind)                      :: enthAir                ! enthalpy of air (J m-3)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="T2enthTemp_soil/"

  Tcrit      = crit_soilT( mLayerMatricHead )
  volFracWat = volFracLiq(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
  diffT      = mLayerTemp - Tfreeze
  diff0      = Tcrit - Tfreeze

  ! *** compute enthalpy of water for unfrozen conditions
  if(mlayerTemp>=Tcrit)then
    enthLiq= iden_water * Cp_water * volFracWat * diffT
    enthIce= 0._rkind

  ! *** compute enthalpy of water for frozen conditions
  else
    ! *** compute integral of mLayerPsiLiq from Tfreeze to layer temperature
    ! get the unfrozen water content
    integral_unf = ( Tcrit - Tfreeze ) * volFracWat

    ! get the frozen water content
    if(use_lookup)then ! cubic spline interpolation for integral of mLayerPsiLiq from Tfreeze to layer temperature
      ! make associate to the the lookup table
      lookVars: associate(&
        Tk => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%temperature)%lookup,  & ! temperature (K)
        Ly => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%psiLiq_int)%lookup,   & ! integral of mLayerPsiLiq from Tfreeze to Tk (K)
        L2 => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%deriv2)%lookup        & ! second derivative of the interpolating function
        ) ! end associate statement

        ! get the lower limit of the integral
        if(diff0<0._rkind)then
          if(integral_frz_low0>=0)then ! = realMissing if non-compute
            integral_frz_low = integral_frz_low0
          else
            call splint(Tk,Ly,L2,Tcrit,integral_frz_low,dL,err,cmessage)
            if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
          endif
        else ! Tcrit=Tfreeze, i.e. mLayerMatricHeadTrial(ixControlIndex)>0
          integral_frz_low = 0._rkind
        end if
        ! get the upper limit of the integral
        call splint(Tk,Ly,L2,mlayerTemp,integral_frz_upp,dL,err,cmessage)
        if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

      end associate lookVars

    else ! hypergeometric function for integral of mLayerPsiLiq from Tfreeze to layer temperature
      ! get the lower limit of the integral
      if(diff0<0._rkind)then
        if(integral_frz_low0>=0)then ! = realMissing if non-compute
          integral_frz_low = integral_frz_low0
        else
          arg              = (vGn_alpha * mLayerMatricHead)**vGn_n
          gauss_hg_T       = hyp_2F1_real(vGn_m,1._rkind/vGn_n,1._rkind + 1._rkind/vGn_n,-arg)
          integral_frz_low = diff0 * ( (theta_sat - theta_res)*gauss_hg_T + theta_res )
        endif
      else ! Tcrit=Tfreeze, i.e. mLayerMatricHeadTrial(ixControlIndex)>0
        integral_frz_low = 0._rkind
      end if
      ! get the upper limit of the integral
      xConst           = LH_fus/(gravity*Tfreeze)        ! m K-1 (NOTE: J = kg m2 s-2)
      mLayerPsiLiq     = xConst*diffT   ! liquid water matric potential from the Clapeyron eqution, DIFFERENT from the liquid water matric potential used in the flux calculations
      arg              = (vGn_alpha * mLayerPsiLiq)**vGn_n
      gauss_hg_T       = hyp_2F1_real(vGn_m,1._rkind/vGn_n,1._rkind + 1._rkind/vGn_n,-arg)
      integral_frz_upp = diffT * ( (theta_sat - theta_res)*gauss_hg_T + theta_res )
    endif

    enthLiq   = iden_water * Cp_water * (integral_unf + integral_frz_upp - integral_frz_low)
    enthIce   = iden_ice * Cp_ice * ( volFracWat * diffT - (integral_unf + integral_frz_upp - integral_frz_low) )

  endif ! (if frozen conditions)

  enthSoil = soil_dens_intr * Cp_soil * ( 1._rkind - theta_sat ) * diffT
  enthAir  = iden_air * Cp_air * ( 1._rkind - theta_sat - volFracWat ) * diffT

  mLayerEnthTemp = enthLiq + enthIce + enthSoil + enthAir

end subroutine T2enthTemp_soil

! ************************************************************************************************************************
! public subroutine enthTemp_or_enthalpy: add energy associated with thaw/freeze to temperature component of enthalpy to get total enthalpy, H, or vice versa
! ************************************************************************************************************************
subroutine enthTemp_or_enthalpy(&
                      ! input: data structures
                      do_enthTemp2enthalpy,    & ! intent(in):    flag if enthalpy is to be computed from temperature component of enthalpy, or vice versa if false
                      diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                      indx_data,               & ! intent(in):    model indices
                      ! input: ice content change
                      scalarCanopyIce,         & ! intent(in):    value of canopy ice content (kg m-2) or prime ice content (kg m-2 s-1)
                      mLayerVolFracIce,        & ! intent(in):    vector of volumetric fraction of ice (-) or prime volumetric fraction of ice (s-1)
                      ! input/output: enthalpy
                      scalarCanopyH,           & ! intent(inout): enthTemp to enthalpy of the vegetation canopy (J m-3), or vice versa if do_enthTemp2enthalpy false
                      mLayerH,                 & ! intent(inout): enthTemp to enthalpy of each snow+soil layer (J m-3), or vice versa if do_enthTemp2enthalpy false
                      ! output: error control
                      err,message)               ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  ! input: data structures
  logical(lgt),intent(in)          :: do_enthTemp2enthalpy       ! flag if enthalpy is to be computed from temperature component of enthalpy, or vice versa if false
  type(var_dlength),intent(in)     :: diag_data                  ! diagnostic variables for a local HRU
  type(var_ilength),intent(in)     :: indx_data                  ! model indices
  ! input: ice content change
  real(rkind),intent(in)           :: scalarCanopyIce            ! value for canopy ice content (kg m-2) or prime ice content (kg m-2 s-1)
  real(rkind),intent(in)           :: mLayerVolFracIce(:)        ! vector of volumetric fraction of ice (-) or prime volumetric fraction of ice (s-1)
  ! input output: enthalpy
  real(rkind),intent(inout)        :: scalarCanopyH              ! enthTemp to enthalpy of the vegetation canopy (J m-3), or vice versa if do_enthTemp2enthalpy false
  real(rkind),intent(inout)        :: mLayerH(:)                 ! enthTemp to enthalpy of each snow+soil layer (J m-3), or vice versa if do_enthTemp2enthalpy false
  ! output: error control
  integer(i4b),intent(out)         :: err                        ! error code
  character(*),intent(out)         :: message                    ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  integer(i4b)                     :: iState                     ! index of model state variable
  integer(i4b)                     :: iLayer                     ! index of model layer
  integer(i4b)                     :: ixFullVector               ! index within full state vector
  integer(i4b)                     :: ixDomainType               ! name of a given model domain
  integer(i4b)                     :: ixControlIndex             ! index within a given model domain
   ! ------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  associate(&
    ! number of model layers, and layer type
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)            ,& ! intent(in): [i4b]    total number of snow layers
    ! mapping between the full state vector and the state subset
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat    ,& ! intent(in): [i4b(:)] list of indices in the state subset for each state in the full state vector
    ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat    ,& ! intent(in): [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! type of domain, type of state variable, and index of control volume within domain
    ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat ,& ! intent(in): [i4b(:)] [state subset] id of domain for desired model state variables
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat     ,& ! intent(in): [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat         ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
   ! canopy depth
    canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)  & ! intent(in): [dp]     canopy depth (m)
    ) ! end associate statement
    ! -----------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message="enthTemp_or_enthalpy/"

    ! loop through model state variables
    do iState=1,size(ixMapSubset2Full)

      ! -----
      ! - compute indices...
      ! --------------------

      ! get domain type, and index of the control volume within the domain
      ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
      ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
      ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

      ! check an energy state
      if(ixStateType(ixFullVector)==iname_nrgCanair .or. ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)then

        ! get the layer index
        select case(ixDomainType)
          case(iname_cas); cycle ! canopy air space: do nothing (no water stored in canopy air space)
          case(iname_veg)
            if (do_enthTemp2enthalpy)then
              scalarCanopyH = scalarCanopyH - LH_fus * scalarCanopyIce/ canopyDepth
            else 
              scalarCanopyH = scalarCanopyH + LH_fus * scalarCanopyIce/ canopyDepth
            end if
          case(iname_snow)
            iLayer = ixControlIndex
            if (do_enthTemp2enthalpy)then
              mLayerH(iLayer) = mLayerH(iLayer) - iden_ice * LH_fus * mLayerVolFracIce(iLayer)
            else
              mLayerH(iLayer) = mLayerH(iLayer) + iden_ice * LH_fus * mLayerVolFracIce(iLayer)
            end if
          case(iname_soil)
            iLayer = ixControlIndex + nSnow
            if (do_enthTemp2enthalpy)then
              mLayerH(iLayer) = mLayerH(iLayer) - iden_water * LH_fus * mLayerVolFracIce(iLayer)
            else
              mLayerH(iLayer) = mLayerH(iLayer) + iden_water * LH_fus * mLayerVolFracIce(iLayer)
            end if
          case(iname_aquifer); cycle ! aquifer: do nothing (no thermodynamics in the aquifer)
          case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
        end select

      end if  ! if an energy layer
    end do  ! looping through state variables

  end associate

end subroutine enthTemp_or_enthalpy

! ************************************************************************************************************************
! public subroutine enthalpy2T_cas: compute temperature from enthalpy, canopy air space
! ************************************************************************************************************************
subroutine enthalpy2T_cas(&
                      computJac,              & ! intent(in):    flag if computing for Jacobian update
                      scalarCanairEnthalpy,   & ! intent(in):    enthalpy of the canopy air space (J m-3)
                      scalarCanairTemp,       & ! intent(out):   canopy air temperature (K)
                      dCanairTemp_dEnthalpy,  & ! intent(inout): derivative of canopy air temperature with enthalpy
                      err,message)              ! intent(out):   error control
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  logical(lgt),intent(in)          :: computJac             ! flag if computing for Jacobian update
  ! input: enthalpy state variables
  real(rkind),intent(in)           :: scalarCanairEnthalpy  ! enthalpy of the canopy air space (J m-3)
  ! output: temperature diagnostic variables
  real(rkind),intent(out)          :: scalarCanairTemp      ! canopy air temperature (K)
   ! output: derivatives
  real(rkind),intent(inout)        :: dCanairTemp_dEnthalpy ! derivative of canopy air temperature with enthalpy
  ! output: error control
  integer(i4b),intent(out)         :: err                   ! error code
  character(*),intent(out)         :: message               ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="enthalpy2T_cas/"

  scalarCanairTemp = scalarCanairEnthalpy / ( Cp_air*iden_air ) + Tfreeze
  if(computJac) dCanairTemp_dEnthalpy = 1._rkind / ( Cp_air*iden_air )

end subroutine enthalpy2T_cas


! ************************************************************************************************************************
! public subroutine enthalpy2T_veg: compute temperature from enthalpy and total water content, canopy
! ************************************************************************************************************************
subroutine enthalpy2T_veg(&
                      computJac,              & ! intent(in):    flag if computing for Jacobian update
                      canopyDepth,            & ! intent(in):    canopy depth (m)
                      specificHeatVeg,        & ! intent(in):    specific heat of vegetation (J kg-1 K-1)
                      maxMassVegetation,      & ! intent(in):    maximum mass of vegetation (kg m-2)
                      snowfrz_scale,          & ! intent(in):    scaling parameter for the snow freezing curve  (K-1)
                      scalarCanopyEnthalpy,   & ! intent(in):    enthalpy of the vegetation canopy (J m-3)
                      scalarCanopyWat,        & ! intent(in):    canopy total water (kg m-2)
                      scalarCanopyTemp,       & ! intent(inout): canopy temperature (K)
                      dCanopyTemp_dEnthalpy,  & ! intent(inout): derivative of canopy temperature with enthalpy
                      dCanopyTemp_dCanWat,    & ! intent(inout): derivative of canopy temperature with canopy water
                      err,message)              ! intent(out):   error control
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines
  USE snow_utils_module,only:fracliquid         ! compute the fraction of liquid water (snow)
  USE snow_utils_module,only:dFracLiq_dTk       ! differentiate the freezing curve w.r.t. temperature (snow)
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  logical(lgt),intent(in)          :: computJac             ! flag if computing for Jacobian update
  ! input: data structures
  real(rkind),intent(in)           :: canopyDepth           ! canopy depth (m)
  real(rkind),intent(in)           :: specificHeatVeg       ! specific heat of vegetation (J kg-1 K-1)
  real(rkind),intent(in)           :: maxMassVegetation     ! maximum mass of vegetation (kg m-2)
  real(rkind),intent(in)           :: snowfrz_scale         ! scaling parameter for the snow freezing curve  (K-1)
   ! input: enthalpy state variables
  real(rkind),intent(in)           :: scalarCanopyEnthalpy  ! enthalpy of the vegetation canopy (J m-3)
  ! input: water state variables
  real(rkind),intent(in)           :: scalarCanopyWat       ! trial value for canopy total water (kg m-2)
   ! output: temperature diagnostic variables
  real(rkind),intent(inout)        :: scalarCanopyTemp      ! trial value for canopy temperature (K)
  ! output: derivatives
  real(rkind),intent(inout)        :: dCanopyTemp_dEnthalpy ! derivative of canopy temperature with enthalpy
  real(rkind),intent(inout)        :: dCanopyTemp_dCanWat   ! derivative of canopy temperature with canopy water
  ! output: error control
  integer(i4b),intent(out)         :: err                   ! error code
  character(*),intent(out)         :: message               ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
   ! declare local variables
  real(rkind)                      :: T                  ! temperature (K)
  real(rkind)                      :: H                  ! enthalpy (J m-3)
  real(rkind)                      :: diffT              ! temperature difference of temp from Tfreeze
  real(rkind)                      :: integral           ! integral of snow freezing curve
  real(rkind)                      :: fLiq               ! fraction liquid 
  real(rkind)                      :: vec(9)             ! vector of parameters for the enthalpy function
   ! variable derivatives
  real(rkind)                      :: dT_dEnthalpy       ! derivative of temperature with enthalpy state variable
  real(rkind)                      :: dT_dWat            ! derivative of temperature with water state variable
  real(rkind)                      :: dH_dT              ! derivative of enthalpy with temperature
  real(rkind)                      :: dH_dWat            ! derivative of enthalpy with water state variable
  real(rkind)                      :: dfLiq_dT           ! derivative of fraction liquid water with temperature
  real(rkind)                      :: denthIce_dT        ! derivative of enthalpy of ice with temperature
  real(rkind)                      :: denthLiq_dT        ! derivative of enthalpy of liquid water with temperature
  real(rkind)                      :: denthVeg_dT        ! derivative of enthalpy of vegetation with temperature
  real(rkind)                      :: denthIce_dWat      ! derivative of enthalpy of ice with water state variable
  real(rkind)                      :: denthLiq_dWat      ! derivative of enthalpy of liquid water with water state variable
  real(rkind)                      :: denthVeg_dWat      ! derivative of enthalpy of vegetation with water state variable 
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="enthalpy2T_veg/"
 
  ! ***** get temperature if unfrozen vegetation
  T            = scalarCanopyEnthalpy * canopyDepth / ( specificHeatVeg * maxMassVegetation + Cp_water * scalarCanopyWat ) + Tfreeze
  if(computJac)then  
    dT_dEnthalpy = canopyDepth / ( specificHeatVeg * maxMassVegetation + Cp_water * scalarCanopyWat )
    dT_dWat      = -Cp_water * scalarCanopyEnthalpy * canopyDepth / ( specificHeatVeg * maxMassVegetation + Cp_water * scalarCanopyWat )**2_i4b
  endif

  ! ***** iterate to find temperature if ice exists
  if( T<Tfreeze )then
    T = min(scalarCanopyTemp,Tfreeze) ! initial guess

    ! find the root of the function
    ! inputs = function, lower bound, upper bound, initial point, tolerance, integer flag if want detail
    ! and the vector of parameters, not.snow_layers
    vec      = 0._rkind
    vec(1:6) = (/scalarCanopyEnthalpy, canopyDepth, specificHeatVeg, maxMassVegetation, snowfrz_scale, scalarCanopyWat/)
    T = brent(diff_H_veg, T, 200._rkind, Tfreeze, vec)

    ! compute Jacobian terms
    if(computJac)then
    ! NOTE: dintegral_dT = fLiq
      diffT    = T - Tfreeze
      integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
      fLiq     = fracLiquid(T, snowfrz_scale)

      ! w.r.t. temperature, NOTE: dintegral_dT = fLiq
      dfLiq_dT    = dFracLiq_dTk(T,snowfrz_scale)
      denthLiq_dT = Cp_water * scalarCanopyWat * fLiq / canopyDepth
      denthIce_dT = Cp_ice * scalarCanopyWat * (1._rkind - fLiq) / canopyDepth
      denthVeg_dT = specificHeatVeg * maxMassVegetation / canopyDepth
      dH_dT       = denthVeg_dT + denthLiq_dT + denthIce_dT + LH_fus * dfLiq_dT * scalarCanopyWat / canopyDepth

      ! w.r.t. layer water content
      denthLiq_dWat = Cp_water * diffT / canopyDepth
      denthIce_dWat = 0._rkind
      denthVeg_dWat = 0._rkind
      dH_dWat       = denthVeg_dWat + denthLiq_dWat + denthIce_dWat - LH_fus * (1._rkind - fLiq) / canopyDepth

      dT_dEnthalpy = 1._rkind / dH_dT
      dT_dWat      = dH_dWat / dH_dT
    endif
  endif ! (if ice exists)

  ! update temperature and derivatives
  scalarCanopyTemp = T
  if(computJac)then
    dCanopyTemp_dEnthalpy = dT_dEnthalpy
    dCanopyTemp_dCanWat   = dT_dWat
  endif

end subroutine enthalpy2T_veg

! ************************************************************************************************************************
! public subroutine enthalpy2T_snow: compute temperature from enthalpy and total water content, snow layer
! ************************************************************************************************************************
subroutine enthalpy2T_snow(&
                      computJac,         & ! intent(in):    flag if computing for Jacobian update
                      snowfrz_scale,     & ! intent(in):    scaling parameter for the snow freezing curve (K-1)
                      mLayerEnthalpy,    & ! intent(in):    enthalpy of snow+soil layer (J m-3)
                      mLayerVolFracWat,  & ! intent(in):    volumetric total water content (-)
                      mLayerTemp,        & ! intent(inout): layer temperature (K)
                      dTemp_dEnthalpy,   & ! intent(inout): derivative of layer temperature with enthalpy
                      dTemp_dTheta,      & ! intent(inout): derivative of layer temperature with volumetric total water content
                      err,message)         ! intent(out):   error control
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines
  USE snow_utils_module,only:fracliquid    ! compute the fraction of liquid water (snow)
  USE snow_utils_module,only:dFracLiq_dTk  ! differentiate the freezing curve w.r.t. temperature (snow)
 
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  logical(lgt),intent(in)          :: computJac          ! flag if computing for Jacobian update
  ! input: data structures
  real(rkind),intent(in)           :: snowfrz_scale      ! scaling parameter for the snow freezing curve  (K-1)
  ! input: enthalpy state variables
  real(rkind),intent(in)           :: mLayerEnthalpy     ! enthalpy of each snow+soil layer (J m-3)
  ! input: water state variables
  real(rkind),intent(in)           :: mLayerVolFracWat   ! volumetric total water content (-)
  ! output: temperature diagnostic variables
  real(rkind),intent(inout)        :: mLayerTemp         ! layer temperature (K)
  ! output: derivatives
  real(rkind),intent(inout)        :: dTemp_dEnthalpy    ! derivative of layer temperature with enthalpy
  real(rkind),intent(inout)        :: dTemp_dTheta       ! derivative of layer temperature with volumetric total water content
   ! output: error control
  integer(i4b),intent(out)         :: err                ! error code
  character(*),intent(out)         :: message            ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  real(rkind)                      :: T                  ! temperature (K)
  real(rkind)                      :: H                  ! enthalpy (J m-3)
  real(rkind)                      :: diffT              ! temperature difference of temp from Tfreeze
  real(rkind)                      :: integral           ! integral of snow freezing curve
  real(rkind)                      :: fLiq               ! fraction liquid 
  real(rkind)                      :: vec(9)             ! vector of parameters for the enthalpy function
  real(rkind)                      :: l_bound            ! lower bound for the enthalpy function
   ! variable derivatives
  real(rkind)                      :: dT_dEnthalpy       ! derivative of temperature with enthalpy state variable
  real(rkind)                      :: dT_dWat            ! derivative of temperature with water state variable
  real(rkind)                      :: dH_dT              ! derivative of enthalpy with temperature
  real(rkind)                      :: dH_dWat            ! derivative of enthalpy with water state variable
  real(rkind)                      :: dfLiq_dT           ! derivative of fraction liquid water with temperature
  real(rkind)                      :: denthIce_dT        ! derivative of enthalpy of ice with temperature
  real(rkind)                      :: denthLiq_dT        ! derivative of enthalpy of liquid water with temperature
  real(rkind)                      :: denthAir_dT        ! derivative of enthalpy of air with temperature
  real(rkind)                      :: denthIce_dWat      ! derivative of enthalpy of ice with water state variable
  real(rkind)                      :: denthLiq_dWat      ! derivative of enthalpy of liquid water with water state variable
  real(rkind)                      :: denthAir_dWat      ! derivative of enthalpy of air with water state variable
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="enthalpy2T_snow/"

  ! ***** iterate to find temperature, ice always exists
  T  = mLayerTemp ! initial guess, will be less than Tfreeze since was a solution

  ! find the root of the function
  ! inputs = function, lower bound, upper bound, initial point, tolerance, integer flag if want detail
  ! and the vector of parameters, snow_layer
  vec = 0._rkind
  vec(1:3) = (/mLayerEnthalpy, snowfrz_scale, mLayerVolFracWat/)
  if(mLayerEnthalpy>0._rkind)then
    T = Tfreeze+ 0.1_rkind ! need to merge layers, trigger the merge
  else
    l_bound = diff_H_snow(200._rkind, vec)
    if (l_bound > 0._rkind) then
      T = Tfreeze + 0.1_rkind ! need to merge layers, trigger the merge
    else
      T = brent(diff_H_snow, T, 200._rkind, Tfreeze, vec)
    end if
  endif

  ! compute Jacobian terms
  if(computJac)then
    ! NOTE: dintegral_dT = fLiq
    diffT    = T - Tfreeze
    integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
    fLiq     = fracLiquid(T, snowfrz_scale)
 
    ! w.r.t. temperature, NOTE: dintegral_dT = fLiq
    dfLiq_dT    = dFracLiq_dTk(T,snowfrz_scale)
    denthLiq_dT = iden_water * Cp_water * mLayerVolFracWat * fLiq
    denthIce_dT = iden_water * Cp_ice * mLayerVolFracWat * (1._rkind - fLiq)
    denthAir_dT = iden_air * Cp_air * (1._rkind - mLayerVolFracWat * ( (iden_water/iden_ice)*(1._rkind-fLiq) + fLiq ) )
    dH_dT       = denthLiq_dT + denthIce_dT + denthAir_dT + iden_water * LH_fus * dfLiq_dT * mLayerVolFracWat

    ! w.r.t. layer water content
    denthLiq_dWat = iden_water * Cp_water * integral
    denthIce_dWat = iden_water * Cp_ice * ( diffT - integral )
    denthAir_dWat = -iden_air * Cp_air * ( (iden_water/iden_ice)*(diffT-integral) + integral )
    dH_dWat       = denthLiq_dWat + denthIce_dWat + denthAir_dWat - iden_water * LH_fus * (1._rkind - fLiq)

    dT_dEnthalpy = 1._rkind / dH_dT
    dT_dWat      = dH_dWat / dH_dT
  endif

  ! update temperature and derivatives
  mLayerTemp = T
  if(computJac)then  
    dTemp_dEnthalpy = dT_dEnthalpy
    dTemp_dTheta    = dT_dWat
  endif

end subroutine enthalpy2T_snow

! ************************************************************************************************************************
! public subroutine enthalpy2T_soil: compute temperature from enthalpy and total water content, soil layer
! ************************************************************************************************************************
subroutine enthalpy2T_soil(&
                      computJac,                & ! intent(in):    flag if computing for Jacobian update
                      use_lookup,               & ! intent(in):    flag to use the lookup table for soil enthalpy
                      soil_dens_intr,           & ! intent(in):    intrinsic soil density (kg m-3)
                      vGn_alpha,                & ! intent(in):    van Genutchen "alpha" parameter
                      vGn_n,                    & ! intent(in):    van Genutchen "n" parameter
                      theta_sat,                & ! intent(in):    soil porosity (-)
                      theta_res,                & ! intent(in):    soil residual volumetric water content (-)
                      vGn_m,                    & ! intent(in):    van Genutchen "m" parameter (-)
                      ixControlIndex,           & ! intent(in):    index of the control volume within the domain
                      lookup_data,              & ! intent(in):    lookup table data structure
                      mLayerEnthalpy,           & ! intent(in):    enthalpy of each snow+soil layer (J m-3)
                      mLayerMatricHead,         & ! intent(in):    total water matric potential (m)
                      mLayerTemp,               & ! intent(inout): layer temperature (K)
                      dTemp_dEnthalpy,          & ! intent(inout): derivative of layer temperature with enthalpy
                      dTemp_dTheta,             & ! intent(inout): derivative of layer temperature with volumetric total water content
                      dTemp_dPsi0,              & ! intent(inout): derivative of layer temperature with total water matric potential
                      err,message)                ! intent(out):   error control         
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines
  USE spline_int_module,only:splint               ! use for cubic spline interpolation               
  USE soil_utils_module,only:dTheta_dTk           ! differentiate the freezing curve w.r.t. temperature (soil)
  USE soil_utils_module,only:crit_soilT           ! compute critical temperature below which ice exists
  USE soil_utils_module,only:volFracLiq           ! compute volumetric fraction of liquid water
  USE soil_utils_module,only:dTheta_dPsi          ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)

  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  logical(lgt),intent(in)          :: computJac              ! flag if computing for Jacobian update
  logical(lgt),intent(in)          :: use_lookup             ! flag to use the lookup table for soil enthalpy, otherwise use hypergeometric function
  ! input: data structures
  real(rkind),intent(in)           :: soil_dens_intr         ! intrinsic soil density (kg m-3)
  real(rkind),intent(in)           :: vGn_alpha              ! van Genutchen "alpha" parameter
  real(rkind),intent(in)           :: vGn_n                  ! van Genutchen "n" parameter
  real(rkind),intent(in)           :: theta_sat              ! soil porosity (-)
  real(rkind),intent(in)           :: theta_res              ! soil residual volumetric water content (-)
  real(rkind),intent(in)           :: vGn_m                  ! van Genutchen "m" parameter (-)
  integer(i4b),intent(in)          :: ixControlIndex         ! index within a given model domain
  type(zLookup),intent(in)         :: lookup_data            ! lookup tables
  ! input: enthalpy state variables
  real(rkind),intent(in)           :: mLayerEnthalpy         ! enthalpy of each snow+soil layer (J m-3)
  ! input: water state variables
  real(rkind),intent(in)           :: mLayerMatricHead       ! total water matric potential (m)
  ! output: temperature diagnostic variables
  real(rkind),intent(inout)        :: mLayerTemp             ! layer temperature (K)
  ! output: derivatives
  real(rkind),intent(inout)        :: dTemp_dEnthalpy        ! derivative of layer temperature with enthalpy
  real(rkind),intent(inout)        :: dTemp_dTheta           ! derivative of layer temperature with volumetric total water content
  real(rkind),intent(inout)        :: dTemp_dPsi0            ! derivative of layer temperature with total water matric potential
  ! output: error control
  integer(i4b),intent(out)         :: err                    ! error code
  character(*),intent(out)         :: message                ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  character(len=128)               :: cmessage               ! error message in downwind routine
  real(rkind)                      :: Tcrit                  ! temperature where all water is unfrozen (K)
  real(rkind)                      :: volFracWat             ! volumetric fraction of total water, liquid+ice (-)
  real(rkind)                      :: diff0                  ! temperature difference of Tcrit from Tfreeze
  real(rkind)                      :: dTcrit_dPsi0           ! derivative of temperature where all water is unfrozen (K) with matric head
  real(rkind)                      :: dL                     ! derivative of soil lookup table with temperature at layer temperature
  real(rkind)                      :: integral_unf           ! integral of unfrozen soil water content (from Tfreeze to Tcrit)
  real(rkind)                      :: integral_frz_low       ! lower limit of integral of frozen soil water content (from Tfreeze to Tcrit)
  real(rkind)                      :: xConst                 ! constant in the freezing curve function (m K-1)
  real(rkind)                      :: mLayerPsiLiq           ! liquid water matric potential (m)
  real(rkind)                      :: T                      ! temperature (K)
  real(rkind)                      :: H                      ! enthalpy (J m-3)
  real(rkind)                      :: diffT                  ! temperature difference of temp from Tfreeze
  real(rkind)                      :: fLiq                   ! fraction liquid water
  real(rkind)                      :: integral_frz_upp       ! upper limit of integral of frozen soil water content (from Tfreeze to soil temperature)
  real(rkind)                      :: arg                    ! argument of soil hypergeometric function
  real(rkind)                      :: gauss_hg_T             ! soil hypergeometric function result
  real(rkind)                      :: vec(9)                 ! vector of parameters for the enthalpy function
  ! variable derivatives
  real(rkind)                      :: dvolFracWat_dPsi0      ! derivative of the soil water content w.r.t. matric head
  real(rkind)                      :: dintegral_unf_dWat     ! derivative of integral of unfrozen soil water content with water content
  real(rkind)                      :: dintegral_frz_low_dWat ! derivative of integral of frozen soil water content with water content
  real(rkind)                      :: dT_dEnthalpy           ! derivative of temperature with enthalpy state variable
  real(rkind)                      :: dT_dWat                ! derivative of temperature with water state variable
  real(rkind)                      :: dH_dT                  ! derivative of enthalpy with temperature
  real(rkind)                      :: dH_dWat                ! derivative of enthalpy with water state variable
  real(rkind)                      :: dfLiq_dT               ! derivative of fraction liquid water with temperature
  real(rkind)                      :: dintegral_frz_upp_dT   ! derivative of integral of frozen soil water content with temperature
  real(rkind)                      :: denthSoil_dT           ! derivative of enthalpy of soil with temperature
  real(rkind)                      :: denthIce_dT            ! derivative of enthalpy of ice with temperature
  real(rkind)                      :: denthLiq_dT            ! derivative of enthalpy of liquid water with temperature
  real(rkind)                      :: denthAir_dT            ! derivative of enthalpy of air with temperature
  real(rkind)                      :: denthSoil_dWat         ! derivative of enthalpy of soil with water state variable
  real(rkind)                      :: denthIce_dWat          ! derivative of enthalpy of ice with water state variable
  real(rkind)                      :: denthLiq_dWat          ! derivative of enthalpy of liquid water with water state variable
  real(rkind)                      :: denthAir_dWat          ! derivative of enthalpy of air with water state variable
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="enthalpy2T_soil/"

  Tcrit             = crit_soilT(mLayerMatricHead)
  volFracWat        = volFracLiq(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
  dTcrit_dPsi0      = merge(gravity*Tfreeze/LH_fus,0._rkind,mLayerMatricHead<=0._rkind)
  dvolFracWat_dPsi0 = dTheta_dPsi(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) 

  ! ***** get temperature if unfrozen soil
  T            = mLayerEnthalpy / ( iden_water * Cp_water * volFracWat + soil_dens_intr * Cp_soil * (1._rkind - theta_sat) &
                                   + iden_air * Cp_air * (1._rkind - theta_sat - volFracWat) ) + Tfreeze
  if(computJac)then  
    dT_dEnthalpy = 1._rkind / ( iden_water * Cp_water * volFracWat + soil_dens_intr*Cp_soil*(1._rkind - theta_sat) &
                               + iden_air*Cp_air*(1._rkind - theta_sat - volFracWat) )
    dT_dWat      = -iden_water * Cp_water * dvolFracWat_dPsi0 * mLayerEnthalpy / ( iden_water * Cp_water * volFracWat &
                                     + soil_dens_intr * Cp_soil * (1._rkind - theta_sat) + iden_air * Cp_air * (1._rkind - theta_sat - volFracWat) )**2_i4b
  endif

  ! ***** iterate to find temperature if ice exists
  if( T<Tcrit )then
    T = min(mLayerTemp,Tcrit) ! initial guess

    ! *** compute integral of mLayerPsiLiq from Tfreeze to layer temperature
    ! get the unfrozen water content of enthalpy
    integral_unf       = ( Tcrit - Tfreeze ) * volFracWat ! unfrozen water content
    if(computJac) dintegral_unf_dWat = dTcrit_dPsi0 * volFracWat + ( Tcrit - Tfreeze ) * dTheta_dPsi(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)

    ! get the frozen water content of enthalpy, statrt with lower limit of the integral
    diff0 = Tcrit - Tfreeze
    if (diff0<0._rkind)then

      if(use_lookup)then ! cubic spline interpolation for integral of mLayerPsiLiq from Tfreeze to layer temperature
        ! make associate to the the lookup table
        lookVars: associate(&
          Tk => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%temperature)%lookup,  & ! temperature (K)
          Ly => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%psiLiq_int)%lookup,   & ! integral of mLayerPsiLiq from Tfreeze to Tk (K)
          L2 => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%deriv2)%lookup        & ! second derivative of the interpolating function
          ) ! end associate statement

          call splint(Tk,Ly,L2,Tcrit,integral_frz_low,dL,err,cmessage)
          if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
          if(computJac) dintegral_frz_low_dWat = dL * dTcrit_dPsi0

        end associate lookVars

      else ! hypergeometric function for integral of mLayerPsiLiq from Tfreeze to layer temperature
        arg              = (vGn_alpha * mLayerMatricHead)**vGn_n
        gauss_hg_T       = hyp_2F1_real(vGn_m,1._rkind/vGn_n,1._rkind + 1._rkind/vGn_n,-arg)
        integral_frz_low = diff0 * ( (theta_sat - theta_res)*gauss_hg_T + theta_res )
        if(computJac) dintegral_frz_low_dWat = volFracWat * dTcrit_dPsi0
      endif
    else ! Tcrit=Tfreeze, i.e. mLayerMatricHead>0
      integral_frz_low       = 0._rkind 
      dintegral_frz_low_dWat = 0._rkind
    end if

    ! find the root of the function
    ! inputs = function, lower bound, upper bound, initial point, tolerance, integer flag if want detail
    ! and the vector of parameters, not.snow_layer, lookup data
    vec(1:9) = (/mLayerEnthalpy, soil_dens_intr, vGn_alpha, vGn_n, theta_sat, theta_res, vGn_m, integral_frz_low, mLayerMatricHead/)
    T = brent(diff_H_soil, T, 200._rkind, Tcrit, vec, use_lookup, lookup_data, ixControlIndex)

  ! compute Jacobian terms
    if(computJac)then
      ! NOTE: here fLiq is the total liquid fraction, not fraction of water fraction that is liquid
      xConst       = LH_fus/(gravity*Tfreeze)        ! m K-1 (NOTE: J = kg m2 s-2)
      diffT        = T - Tfreeze
      mLayerPsiLiq = xConst*diffT   ! liquid water matric potential from the Clapeyron eqution, DIFFERENT from the liquid water matric potential used in the flux calculations
      arg          = (vGn_alpha * mLayerPsiLiq)**vGn_n
      fLiq         = volFracLiq(mLayerPsiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) 

      ! get the upper limit of the integral
      if(use_lookup)then ! cubic spline interpolation for integral of mLayerPsiLiq from Tfreeze to layer temperature
        ! make associate to the the lookup table
        lookVars2: associate(&
          Tk => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%temperature)%lookup,  & ! temperature (K)
          Ly => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%psiLiq_int)%lookup,   & ! integral of mLayerPsiLiq from Tfreeze to Tk (K)
          L2 => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%deriv2)%lookup        & ! second derivative of the interpolating function
          ) ! end associate statement

          ! integral of mLayerPsiLiq from Tfreeze to layer temperature
          call splint(Tk,Ly,L2,T,integral_frz_upp,dL,err,cmessage)
          if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
          dintegral_frz_upp_dT = dL

        end associate lookVars2

      else ! hypergeometric function for integral of mLayerPsiLiq from Tfreeze to layer temperature
        gauss_hg_T             = hyp_2F1_real(vGn_m,1._rkind/vGn_n,1._rkind + 1._rkind/vGn_n,-arg)
        integral_frz_upp       = diffT * ( (theta_sat - theta_res)*gauss_hg_T + theta_res )
        dintegral_frz_upp_dT   = fLiq
      end if

      ! w.r.t. temperature
      dfLiq_dT     = dTheta_dTk(T,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
      denthLiq_dT  = iden_water * Cp_water * dintegral_frz_upp_dT
      denthIce_dT  = iden_ice * Cp_ice * ( volFracWat  - dintegral_frz_upp_dT ) 
      denthSoil_dT = soil_dens_intr * Cp_soil * ( 1._rkind - theta_sat )
      denthAir_dT  = iden_air * Cp_air * ( 1._rkind - theta_sat - volFracWat ) 
      dH_dT        = denthLiq_dT + denthIce_dT + denthSoil_dT + denthAir_dT + iden_water * LH_fus * dfLiq_dT

      ! w.r.t. layer water content
      denthLiq_dWat  = iden_water * Cp_water * (dintegral_unf_dWat - dintegral_frz_low_dWat)
      denthIce_dWat  = iden_ice * Cp_ice * ( dvolFracWat_dPsi0 * diffT - (dintegral_unf_dWat - dintegral_frz_low_dWat) )
      denthSoil_dWat = 0._rkind
      denthAir_dWat  = -iden_air * Cp_air * dvolFracWat_dPsi0 * diffT      
      dH_dWat        = denthLiq_dWat + denthIce_dWat + denthAir_dWat - iden_water * LH_fus * dvolFracWat_dPsi0  

      dT_dEnthalpy = 1._rkind / dH_dT
      dT_dWat      = dH_dWat / dH_dT
    endif
  end if ! (if ice exists)

  ! update temperature and derivatives
  mLayerTemp = T
  if(computJac)then  
    dTemp_dEnthalpy = dT_dEnthalpy
    dTemp_dTheta    = realMissing ! do not use
    dTemp_dPsi0     = dT_dWat
  endif

end subroutine enthalpy2T_soil

!----------------------------------------------------------------------
! private function: compute hypergeometric function with real arguments into real result
!----------------------------------------------------------------------
 function hyp_2F1_real(a_real, b_real, c_real, z_real)
  !--------------------------------------------------------------------
  USE hyp_2F1_module,only:HYP_2F1 ! use for hypergeometric function
  implicit none
  real(rkind),intent(in) :: a_real, b_real, c_real, z_real
  complex(rkind)         :: a_complex, b_complex, c_complex, z_complex, result
  real(rkind)            :: hyp_2F1_real
  
  a_complex = CMPLX(a_real, 0._rkind, rkind)
  b_complex = CMPLX(b_real, 0._rkind, rkind)
  c_complex = CMPLX(c_real, 0._rkind, rkind)
  z_complex = CMPLX(z_real, 0._rkind, rkind)
  result = HYP_2F1(a_complex, b_complex, c_complex, z_complex)
  hyp_2F1_real = REAL(result, rkind)
   
end function hyp_2F1_real

!----------------------------------------------------------------------
! private function: Brent's method to find a root of a function
!----------------------------------------------------------------------
function brent0 (fun, x1, x2, fx1, fx2, tol_x, tol_f, detail, vec, use_lookup, lookup_data, ixControlIndex)
  !
  ! Description of algorithm: 
  ! Find a root of function f(x) given intial bracketing interval [a,b]
  ! where f(a) and f(b) must have opposite signs. At a typical step we have
  ! three points a, b, and c such that f(b)f(c)<0, and a may coincide with
  ! c. The points a, b, and c change during the algorithm, and the root
  ! always lies in either [b,c] or [c, b]. The value b is the best
  ! approximation to the root and a is the previous value of b. 
  ! 
  ! The iteration uses following selection of algorithms 
  ! when bracket shrinks reasonablly fast,
  !  - Linear interporation if a == b
  !  - Quadratic interporation if a != b and the point is in the bracket.
  ! othrwise 
  ! - Bisection.
  ! 
  ! Inputs:
  !   fun: function to be solved
  !   tol_x: tolerance for x
  !   tol_f: tolerance for f(x)
  !   x1, x2: Upper bound and lower bound for a function
  !   detail: output result of iteration if detail is >0 
  !  
  ! Based on zeroin.f in netlib
  ! modified from fzero.f90 by Yoki Okawa, Jan 30, 2009

  implicit none
  real(rkind) :: brent0
  integer, parameter :: d = rkind
  real(rkind), intent(IN) :: x1, x2, fx1, fx2, vec(9), tol_x, tol_f
  real(rkind), external :: fun
  integer, intent(IN) :: detail
  logical(lgt), intent(in), optional :: use_lookup
  type(zLookup),intent(in), optional :: lookup_data
  integer(i4b), intent(in), optional :: ixControlIndex
  
  integer :: i, exitflag, disp
  real(rkind) :: a, b, c, diff,e, fa, fb, fc, p, q, r, s, tol1, xm, tmp
  real(rkind), parameter :: EPS = epsilon(a)
  integer, parameter :: imax = 100  ! maximum number of iteration
  
  exitflag = 0
  if (detail /= 0) then
    disp = 1
  else
    disp = 0
  end if
    
  ! intialize values
  a = x1
  b = x2
  c = x2
  fa = fx1
  fb = fx2
  fc = fx2
    
  ! check sign
  if ( (fa>0. .and. fb>0. )  .or.  (fa>0. .and. fb>0. )) then
    write(*,*)  'Error (brent0.f90): Root must be bracketed by two inputs'
    write(*, "(' x1 = ', 1F8.4, ' x2 = ', 1F8.4, ' f(x1) = ', 1F15.4, ' f(x2) = ', 1F15.4)") a,b,fa,fb
    write(*,*) 'press any key to halt the program'
    read(*,*)
    stop
  end if
  
  if (disp == 1 ) then 
    write(*,*) 'Brents method to find a root of f(x)'
    write(*,*) ' '
    write(*,*) '  i           x          bracketsize       f(x)'
  end if
  
  ! main iteration
  do i = 1, imax
    ! rename c and adjust bounding interval if both a(=b) and c are same sign
    if ((fb > 0.  .and. fc > 0) .or. (fb <0. .and. fc < 0. ) ) then 
      c = a
      fc = fa
      e = b-a
      diff = e
    end if
    
    ! if c is better guess than b, use it. 
    if (abs(fc) < abs(fb) ) then 
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    end if

    if (disp == 1) then 
      tmp = c-b
      write(*,"('  ', 1I2, 3F16.10)") i, b, abs(b-c), fb
    end if
    
    ! convergence check
    tol1=2.0_rkind* EPS * abs(b) + 0.5_rkind*tol_x
    xm = 0.5_rkind * (c - b)
    if (abs(xm) < tol1 .or. abs(fb) <= tol_f )  then
      exitflag = 1
      exit
    end if
     
    ! try inverse quadratic interpolation
    if (abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then 
      s = fb/fa
        if (abs(a - c) < EPS) then 
        p = 2.0_rkind *xm * s
        q = 1.0_rkind  - s
      else
        q = fa/fc
        r = fb/fc
        p = s * (2.0_rkind * xm * q * (q -r ) - (b - a) * (r - 1.0_rkind))
        q = (q - 1.0_rkind ) * (r - 1.0_rkind) * (s - 1.0_rkind) 
      end if
      
      ! accept if q is not too small to stay in bound
      if (p > 0.0_rkind) q = -q
      p = abs(p)                
      if (2.0 * p < min(3.0 * xm * q - abs(tol1* q), abs(e *q))) then 
        e = d
        diff = p / q
      else   ! interpolation failed. use bisection
        diff= xm 
        e = d
      end if
    else  ! quadratic interpolation bounds moves too slowly, use bisection
      diff = xm
      e = d
    end if
    
    ! update last bound
    a = b
    fa = fb
    
    ! move the best guess
    if (abs(d) > tol1) then 
      b = b + diff
    else
      b = b + sign(tol1, xm)
    end if
    
    ! evaluate new trial root
    if(present(use_lookup))then
      fb = fun(b, vec, use_lookup, lookup_data, ixControlIndex)
    else
      fb = fun(b, vec)
    end if

  end do 
  
  ! case for non convergence
  if (exitflag /= 1 ) then 
    write(*,*) 'Error (brent0.f90) :  convergence was not attained'
    write(*,*) 'Initial value:'
    write(*,"(4F10.5)" )   x1, x2, fx1, fx2
    write(*,*) ' '
    write(*,*) 'final value:'
    write(*,"('x = '  ,1F6.4, ':                f(x1) = ' ,  1F6.4  )" )  b,  fb  
  else if( disp == 1) then
    write(*,*) 'Brents method was converged.'
    write(*,*) ''
  end if
  brent0 = b
  return
  
  end function brent0

!----------------------------------------------------------------------
! private function: Find an initial guess of bracket and call brent0
!----------------------------------------------------------------------
  function brent (fun, x0, LowerBound, UpperBound, vec, use_lookup, lookup_data, ixControlIndex)
    ! 
    ! Inputs
    !   fun: function to evaluate
    !   x0: Initial guess
    !   LowerBound, UpperBound : Lower and upper bound of the function 
    
    implicit none
    real(rkind) :: brent
    integer, parameter :: d = rkind
    real(rkind), intent(IN) :: x0, vec(9)
    real(rkind), external :: fun
    real(rkind), intent(IN) :: LowerBound, UpperBound
    logical(lgt), intent(in), optional :: use_lookup
    type(zLookup),intent(in), optional :: lookup_data
    integer(i4b), intent(in), optional :: ixControlIndex
    
    real(rkind) :: a , b , olda, oldb, fa, fb, folda, foldb
    real(rkind), parameter :: sqrt2 = sqrt(2.0_d)! change in dx
    integer, parameter :: maxiter = 40, detail = 0
    real(rkind) :: dx  ! change in bracket
    integer :: iter, exitflag, disp
    real(rkind) :: sgn
    real(rkind), parameter :: tol_x = 1.e-5_rkind, tol_f = 1.e0_rkind

    a  = x0 ! lower bracket
    b =  x0 ! upper bracket
    exitflag = 0  ! flag to see we found the bracket

    if(present(use_lookup))then
      sgn = fun(x0, vec, use_lookup, lookup_data, ixControlIndex) ! sign of initial guess
    else  
      sgn = fun(x0, vec)
    endif
    ! set disp variable
    if (detail /= 0) then
      disp = 1
    else
      disp = 0
    end if
    fa = sgn
    fb = sgn

    if(abs(sgn) <= tol_f ) then ! if solution didn't change, initial guess is the solution
      brent = x0
      return
    end if  
    
    ! set initial change dx
    if (abs(x0)<0.00000002_rkind) then 
      dx = 1.0_rkind/50.0_rkind
    else
      dx = 1.0_rkind/50.0_rkind * x0
    end if
    
    if (disp == 1) then 
      write(*,*) 'Search for initial guess for Brents method'
      write(*,*) 'find two points whose sign for f(x) is different '
      write(*,*) 'x1 searches downwards, x2 searches upwards with increasing increment'
      write(*,*) ' '
      write(*,*) '  i           x1               x2            f(x1)            f(x2)'
      write(*,"(1I4,4F17.6)") 0, a, b, fa, fb
    end if
    
    ! main loop to extend a and b
    do iter = 1, maxiter
      ! update boundary
      olda = a 
      oldb = b
      folda = fa
      foldb = fb
      a = a - dx
      b = b + dx
      dx = dx * sqrt2
      
      ! boundary check
      if (a < LowerBound ) a = LowerBound
      if (b > UpperBound ) b = UpperBound

      if(present(use_lookup))then
        fa = fun(a, vec, use_lookup, lookup_data, ixControlIndex)
        fb = fun(b, vec, use_lookup, lookup_data, ixControlIndex)
      else  
        fa = fun(a, vec)
        fb = fun(b, vec)
      end if
      
      if (disp == 1) write(*,"(1I4,4F17.6)") iter, a, b, fa, fb
      
      ! check if sign of functions changed or not
      if (( (sgn >= 0 ) .and.  (fa <= 0) ) .or. & 
          ( (sgn <= 0 ) .and.  (fa >= 0  ) ))then  ! sign of a changed 
        ! use a and olda as bracket
        b = olda
        fb = folda
        exitflag = 1
        exit
      else if  (( (sgn >= 0 ) .and.  (fb <= 0  ) ) .or. & 
                ( (sgn <= 0 ) .and.  (fb >= 0  ) )) then ! sign of b changed
        a = oldb
        fa = foldb
        exitflag = 1
        exit
      end if
      
    end do
    
    ! case for non convergence
    if (exitflag /=  1 ) then   
      write(*,*) ' Error (brent2) : Proper initial value for Brents method could not be found'
      write(*,*) ' Change initial guess and try again. '
      write(*,*) ' You might want to try disp = 1 option too'
      write(*,*) '  i           x1               x2            f(x1)            f(x2)'
      write(*,"(1I4,4F17.6)") iter, a, b, fa, fb
      write(*,*) '  press any key to abort the program'
      read(*,*) 
      stop
    else if (disp == 1) then
      write(*,*) '  Initial guess was found.'
      write(*,*) ''
    end if
    
    ! call brent0
    if(present(use_lookup))then
      brent = brent0(fun, a, b, fa, fb, tol_x, tol_f, detail, vec, use_lookup, lookup_data, ixControlIndex)
    else
      brent = brent0(fun, a, b, fa, fb, tol_x, tol_f, detail, vec)
    end if
    
    end function brent  

  !----------------------------------------------------------------------
  ! private functions for temperature to enthalpy conversion for Brent's method
  !----------------------------------------------------------------------
  function diff_H_veg ( scalarCanopyTemp, vec)
    USE snow_utils_module,only:fracliquid     ! compute volumetric fraction of liquid water
    implicit none
    real(rkind) :: diff_H_veg
    real(rkind) , intent(IN) :: scalarCanopyTemp, vec(8) 
    real(rkind) :: scalarCanopyEnthalpy, scalarCanopyEnthTemp, scalarCanopyWat, scalarCanopyIce
    real(rkind) :: canopyDepth, specificHeatVeg, maxMassVegetation, snowfrz_scale, fLiq
    integer(i4b) :: err
    character(256) :: cmessage
  
    scalarCanopyEnthalpy = vec(1)
    canopyDepth          = vec(2)
    specificHeatVeg      = vec(3)
    maxMassVegetation    = vec(4)
    snowfrz_scale        = vec(5)
    scalarCanopyWat      = vec(6)
  
    call T2enthTemp_veg(canopyDepth, specificHeatVeg, maxMassVegetation, snowfrz_scale, scalarCanopyTemp, &
                        scalarCanopyWat, scalarCanopyEnthTemp, err, cmessage)
    fLiq  = fracliquid(scalarCanopyTemp, snowfrz_scale)
    diff_H_veg = scalarCanopyEnthTemp - LH_fus * scalarCanopyWat* (1._rkind - fLiq)/ canopyDepth - scalarCanopyEnthalpy
  
  end function diff_H_veg
  !----------------------------------------------------------------------
  function diff_H_snow ( mLayerTemp, vec)
    USE snow_utils_module,only:fracliquid     ! compute volumetric fraction of liquid water
    implicit none
    real(rkind) :: diff_H_snow
    real(rkind) , intent(IN) :: mLayerTemp, vec(9) 
    real(rkind) :: mLayerEnthalpy, mLayerEnthTemp, mLayerVolFracWat, mLayerVolFracIce, snowfrz_scale, fLiq
    integer(i4b) :: err
    character(256) :: cmessage
  
    mLayerEnthalpy   = vec(1)
    snowfrz_scale    = vec(2)
    mLayerVolFracWat = vec(3)
  
    call T2enthTemp_snow(snowfrz_scale, mLayerTemp, mLayerVolFracWat, mLayerEnthTemp, err, cmessage)
    fLiq   = fracliquid(mLayerTemp, snowfrz_scale)
    diff_H_snow = mLayerEnthTemp - iden_water * LH_fus * mLayerVolFracWat * (1._rkind - fLiq) - mLayerEnthalpy
  
  end function diff_H_snow
  !----------------------------------------------------------------------
  function diff_H_soil ( mLayerTemp, vec, use_lookup, lookup_data, ixControlIndex)
    USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water based on matric head
    implicit none
    real(rkind) :: diff_H_soil
    real(rkind) , intent(in) :: mLayerTemp, vec(9) 
    logical(lgt), intent(in) :: use_lookup
    type(zLookup),intent(in) :: lookup_data
    integer(i4b), intent(in) :: ixControlIndex
    real(rkind) :: mLayerEnthalpy, mLayerEnthTemp, mLayerMatricHead, volFracWat, xConst, mLayerPsiLiq, fLiq
    real(rkind) :: soil_dens_intr, vGn_alpha, vGn_n, theta_sat, theta_res, vGn_m, integral_frz_low
    integer(i4b) :: err
    character(256) :: cmessage
  
    mLayerEnthalpy   = vec(1)
    soil_dens_intr   = vec(2)
    vGn_alpha        = vec(3)
    vGn_n            = vec(4)
    theta_sat        = vec(5)
    theta_res        = vec(6)
    vGn_m            = vec(7)
    integral_frz_low = vec(8)
    mLayerMatricHead = vec(9)
  
    call T2enthTemp_soil(use_lookup, soil_dens_intr, vGn_alpha, vGn_n, theta_sat, theta_res, vGn_m, &
                         ixControlIndex, lookup_data, integral_frz_low, mLayerTemp, mLayerMatricHead, &
                         mLayerEnthTemp, err, cmessage)

    volFracWat   = volFracLiq(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    xConst       = LH_fus/(gravity*Tfreeze)        ! m K-1 (NOTE: J = kg m2 s-2)
    mLayerPsiLiq = xConst*(mLayerTemp - Tfreeze)   ! liquid water matric potential from the Clapeyron eqution
    fLiq         = volFracLiq(mLayerPsiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) ! here fLiq is the total liquid fraction, not fraction of water fraction that is liquid
    diff_H_soil = mLayerEnthTemp - iden_water * LH_fus * (volFracWat - fLiq) - mLayerEnthalpy
  
  end function diff_H_soil


end module enthalpyTemp_module