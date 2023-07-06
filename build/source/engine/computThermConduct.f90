
module computThermConduct_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength         ! data vector with variable length dimension (rkind)

! physical constants
USE multiconst,only:&
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    ! thermal conductivity
                    lambda_air,  & ! thermal conductivity of air   (J s-1 m-1)
                    lambda_ice,  & ! thermal conductivity of ice   (J s-1 m-1)
                    lambda_water   ! thermal conductivity of water (J s-1 m-1)

! missing values
USE globalData,only:integerMissing ! missing integer
USE globalData,only:realMissing    ! missing real number

! named variables that define the layer type
USE globalData,only:iname_snow     ! snow
USE globalData,only:iname_soil     ! soil

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements


! provide access to named variables for thermal conductivity of soil
USE globalData,only:model_decisions  ! model decision structure
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure

! provide access to look-up values for model decisions
USE mDecisions_module,only:      &
 ! look-up values for choice of thermal conductivity representation for snow
 Yen1965,                        & ! Yen (1965)
 Mellor1977,                     & ! Mellor (1977)
 Jordan1991,                     & ! Jordan (1991)
 Smirnova2000,                   & ! Smirnova et al. (2000)
 ! look-up values for choice of thermal conductivity representation for soil
 funcSoilWet,                    & ! function of soil wetness
 mixConstit,                     & ! mixture of constituents
 hanssonVZJ,                     & ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
 ! look-up values for the form of Richards' equation
 moisture,                       & ! moisture-based form of Richards' equation
 mixdform                          ! mixed form of Richards' equation

! privacy
implicit none
private
public::computThermConduct

! algorithmic parameters
real(rkind),parameter     :: valueMissing=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
real(rkind),parameter     :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers
real(rkind),parameter     :: mpe=1.e-6_rkind         ! prevents overflow error if division by zero
real(rkind),parameter     :: dx=1.e-6_rkind          ! finite difference increment
contains

! **********************************************************************************************************
! public subroutine computThermConduct: recompute diagnostic energy variables (thermal conductivity and heat capacity)
! **********************************************************************************************************
subroutine computThermConduct(&
                    ! input: control variables
                    computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                    nLayers,                 & ! intent(in):    total number of layers
                    canopyDepth,             & ! intent(in):    canopy depth (m)
                    ! input: state variables
                    scalarCanopyIce,         & ! intent(in):    canopy ice content (kg m-2)
                    scalarCanopyLiquid,      & ! intent(in):    canopy liquid water content (kg m-2)
                    mLayerTemp,              & ! intent(in):    temperature at the current iteration (K)
                    mLayerMatricHead,        & ! intent(in):    matric head at the current iteration(m)
                    mLayerVolFracIce,        & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                    mLayerVolFracLiq,        & ! intent(in):    volumetric fraction of liquid water at the start of the sub-step (-)
                    ! input: pre-computed derivatives
                    mLayerdTheta_dTk,        & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                    mLayerFracLiqSnow,       & ! intent(in):    fraction of liquid water (-)
                    ! input/output: data structures
                    mpar_data,               & ! intent(in):    model parameters
                    indx_data,               & ! intent(in):    model layer indices
                    prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                    diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                    ! output: derivatives
                    dThermalC_dWatAbove,     & ! intent(out):   derivative in the thermal conductivity w.r.t. water state in the layer above
                    dThermalC_dWatBelow,     & ! intent(out):   derivative in the thermal conductivity w.r.t. water state in the layer above
                    dThermalC_dTempAbove,    & ! intent(out):   derivative in the thermal conductivity w.r.t. energy state in the layer above
                    dThermalC_dTempBelow,    & ! intent(out):   derivative in the thermal conductivity w.r.t. energy state in the layer above
                    ! output: error control
                    err,message)               ! intent(out):   error control

  ! utility modules
  USE snow_utils_module,only:tcond_snow     ! compute thermal conductivity of snow
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  USE soil_utils_module,only:dTheta_dPsi    ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
  USE soil_utils_module,only:dPsi_dTheta    ! compute derivative of the soil moisture characteristic w.r.t. theta (m)

  implicit none
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  logical(lgt),intent(in)              :: computeVegFlux           ! logical flag to denote if computing the vegetation flux
  integer(i4b),intent(in)              :: nLayers                  ! total number of layers in the snow+soil domain
  real(rkind),intent(in)               :: canopyDepth              ! depth of the vegetation canopy (m)
  ! input: trial model state variables
  real(rkind),intent(in)               :: scalarCanopyIce          ! trial value of canopy ice content (kg m-2)
  real(rkind),intent(in)               :: scalarCanopyLiquid
  real(rkind),intent(in)               :: mLayerTemp(:)            ! temperature in each layer at the current iteration (m)
  real(rkind),intent(in)               :: mLayerMatricHead(:)      ! matric head in each layer at the current iteration (m)
  real(rkind),intent(in)               :: mLayerVolFracIce(:)      ! volumetric fraction of ice at the current iteration (-)
  real(rkind),intent(in)               :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid at the current iteration (-)
  ! input: pre-computed derivatives
  real(rkind),intent(in)               :: mLayerdTheta_dTk(:)      ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)               :: mLayerFracLiqSnow(:)     ! fraction of liquid water (-)
  ! input/output: data structures
  type(var_dlength),intent(in)         :: mpar_data                ! model parameters
  type(var_ilength),intent(in)         :: indx_data                ! model layer indices
  type(var_dlength),intent(in)         :: prog_data                ! model prognostic variables for a local HRU
  type(var_dlength),intent(inout)      :: diag_data                ! model diagnostic variables for a local HRU
  ! output: derivatives
  real(rkind),intent(out)              :: dThermalC_dWatAbove(0:)  ! derivative in the thermal conductivity w.r.t. water state in the layer above
  real(rkind),intent(out)              :: dThermalC_dWatBelow(0:)  ! derivative in the thermal conductivity w.r.t. water state in the layer above
  real(rkind),intent(out)              :: dThermalC_dTempAbove(0:) ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  real(rkind),intent(out)              :: dThermalC_dTempBelow(0:) ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  ! output: error control
  integer(i4b),intent(out)             :: err                      ! error code
  character(*),intent(out)             :: message                  ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  character(LEN=256)                   :: cmessage                 ! error message of downwind routine
  integer(i4b)                         :: iLayer                   ! index of model layer
  integer(i4b)                         :: iSoil                    ! index of soil layer
  real(rkind)                          :: TCn                      ! thermal conductivity below the layer interface (W m-1 K-1)
  real(rkind)                          :: TCp                      ! thermal conductivity above the layer interface (W m-1 K-1)
  real(rkind)                          :: zdn                      ! height difference between interface and lower value (m)
  real(rkind)                          :: zdp                      ! height difference between interface and upper value (m)
  real(rkind)                          :: bulkden_soil             ! bulk density of soil (kg m-3)
  real(rkind)                          :: lambda_drysoil           ! thermal conductivity of dry soil (W m-1)
  real(rkind)                          :: lambda_wetsoil           ! thermal conductivity of wet soil (W m-1)
  real(rkind)                          :: lambda_wet               ! thermal conductivity of the wet material
  real(rkind)                          :: relativeSat              ! relative saturation (-)
  real(rkind)                          :: kerstenNum               ! the Kersten number (-), defining weight applied to conductivity of the wet medium
  real(rkind)                          :: den                      ! denominator in the thermal conductivity calculations
  real(rkind)                          :: Tcrit                    ! temperature where all water is unfrozen (K)
  real(rkind),dimension(nLayers)       :: dThermalC_dWat           ! derivative in thermal conductivity w.r.t. matric head or volumetric liquid water
  real(rkind),dimension(nLayers)       :: dThermalC_dNrg           ! derivative in thermal conductivity w.r.t. temperature
  real(rkind)                          :: dlambda_wet_dWat         ! derivative in thermal conductivity of wet material w.r.t.soil water state variable
  real(rkind)                          :: dlambda_wet_dTk          ! derivative in thermal conductivity of wet material w.r.t. temperature
  real(rkind)                          :: dkerstenNum_dWat         ! derivative in Kersten number w.r.t. soil water state variable
  real(rkind)                          :: dVolFracLiq_dWat         ! derivative in vol fraction of liquid w.r.t. water state variable
  real(rkind)                          :: dVolFracIce_dWat         ! derivative in vol fraction of ice w.r.t. water state variable
  real(rkind)                          :: dVolFracLiq_dTk          ! derivative in vol fraction of liquid w.r.t. temperature
  real(rkind)                          :: dVolFracIce_dTk          ! derivative in vol fraction of ice w.r.t. temperature
  ! local variables to reproduce the thermal conductivity of Hansson et al. VZJ 2005
  real(rkind),parameter                :: c1=0.55_rkind            ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
  real(rkind),parameter                :: c2=0.8_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
  real(rkind),parameter                :: c3=3.07_rkind            ! optimized parameter from Hansson et al. VZJ 2005 (-)
  real(rkind),parameter                :: c4=0.13_rkind            ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
  real(rkind),parameter                :: c5=4._rkind              ! optimized parameter from Hansson et al. VZJ 2005 (-)
  real(rkind),parameter                :: f1=13.05_rkind           ! optimized parameter from Hansson et al. VZJ 2005 (-)
  real(rkind),parameter                :: f2=1.06_rkind            ! optimized parameter from Hansson et al. VZJ 2005 (-)
  real(rkind)                          :: fArg,xArg                ! temporary variables (see Hansson et al. VZJ 2005 for details)
  real(rkind)                          :: dxArg_dWat,dxArg_dTk     ! derivates of the temporary variables with respect to soil water state variable and temperature

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! associate variables in data structure
  associate(&
    ! input: model decisions
    ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision,      & ! intent(in):  [i4b]   index of the form of Richards' equation
    ixThCondSnow            => model_decisions(iLookDECISIONS%thCondSnow)%iDecision,      & ! intent(in):  [i4b]   choice of method for thermal conductivity of snow
    ixThCondSoil            => model_decisions(iLookDECISIONS%thCondSoil)%iDecision,      & ! intent(in):  [i4b]   choice of method for thermal conductivity of soil
    ! input: coordinate variables
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),                    & ! intent(in):  [dp]    number of snow layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1),                    & ! intent(in):  [dp]    number of soil layers
    layerType               => indx_data%var(iLookINDEX%layerType)%dat,                   & ! intent(in):  [dp(:)] layer type (iname_soil or iname_snow)
    mLayerHeight            => prog_data%var(iLookPROG%mLayerHeight)%dat,                 & ! intent(in):  [dp(:)] height at the mid-point of each layer (m)
    iLayerHeight            => prog_data%var(iLookPROG%iLayerHeight)%dat,                 & ! intent(in):  [dp(:)] height at the interface of each layer (m)
    ! input: heat capacity and thermal conductivity
    specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),          & ! intent(in):  [dp]    specific heat of vegetation (J kg-1 K-1)
    maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),        & ! intent(in):  [dp]    maximum mass of vegetation (kg m-2)
    fixedThermalCond_snow   => mpar_data%var(iLookPARAM%fixedThermalCond_snow)%dat(1),    & ! intent(in):  [dp]    temporally constant thermal conductivity of snow (W m-1 K-1)
    ! input: depth varying soil parameters
    iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat,              & ! intent(in):  [dp(:)] intrinsic density of soil (kg m-3)
    thCond_soil             => mpar_data%var(iLookPARAM%thCond_soil)%dat,                 & ! intent(in):  [dp(:)] thermal conductivity of soil (W m-1 K-1)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                   & ! intent(in):  [dp(:)] soil porosity (-)
    frac_sand               => mpar_data%var(iLookPARAM%frac_sand)%dat,                   & ! intent(in):  [dp(:)] fraction of sand (-)
    frac_silt               => mpar_data%var(iLookPARAM%frac_silt)%dat,                   & ! intent(in):  [dp(:)] fraction of silt (-)
    frac_clay               => mpar_data%var(iLookPARAM%frac_clay)%dat,                   & ! intent(in):  [dp(:)] fraction of clay (-)
    vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat,                  & ! intent(in):  [dp(:)] van Genutchen "m" parameter (-)
    vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat,                       & ! intent(in):  [dp(:)] van Genutchen "n" parameter (-)
    vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat,                   & ! intent(in):  [dp(:)] van Genutchen "alpha" parameter (m-1)
    theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat,                   & ! intent(in):  [dp(:)] soil residual volumetric water content (-)
    ! input: snow parameters
    snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),            & ! intent(in):  [dp]    scaling parameter for the snow freezing curve (K-1)
    ! output: diagnostic variables and derivatives (diagnostic as may be treated as constant)
    mLayerThermalC          => diag_data%var(iLookDIAG%mLayerThermalC)%dat,               & ! intent(out): [dp(:)] thermal conductivity at the mid-point of each layer (W m-1 K-1)
    iLayerThermalC          => diag_data%var(iLookDIAG%iLayerThermalC)%dat,               & ! intent(out): [dp(:)] thermal conductivity at the interface of each layer (W m-1 K-1)
    mLayerVolFracAir        => diag_data%var(iLookDIAG%mLayerVolFracAir)%dat              & ! intent(out): [dp(:)] volumetric fraction of air in each layer (-)
    )  ! association of local variables with information in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="computThermConduct/"

    ! initialize the soil layer
    iSoil=integerMissing

    ! loop through layers
    do iLayer=1,nLayers

      ! get the soil layer
      if(iLayer>nSnow) iSoil = iLayer-nSnow

      ! compute the thermal conductivity of dry and wet soils (W m-1)
      ! NOTE: this is actually constant over the simulation, and included here for clarity
      if(ixThCondSoil == funcSoilWet .and. layerType(iLayer)==iname_soil)then
        bulkden_soil   = iden_soil(iSoil)*( 1._rkind - theta_sat(iSoil) )
        lambda_drysoil = (0.135_rkind*bulkden_soil + 64.7_rkind) / (iden_soil(iSoil) - 0.947_rkind*bulkden_soil)
        lambda_wetsoil = (8.80_rkind*frac_sand(iSoil) + 2.92_rkind*frac_clay(iSoil)) / (frac_sand(iSoil) + frac_clay(iSoil))
      end if

      ! *****
      ! * compute the volumetric fraction of air in each layer...
      ! *********************************************************
      select case(layerType(iLayer))
        case(iname_soil); mLayerVolFracAir(iLayer) = theta_sat(iSoil) - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
        case(iname_snow); mLayerVolFracAir(iLayer) = 1._rkind - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
        case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute volumetric fraction of air'; return
      end select

      ! *****
      ! * compute the thermal conductivity of snow and soil and derivatives at the mid-point of each layer...
      ! ***************************************************************************************************
      dThermalC_dWat(iLayer) = 0._rkind
      dThermalC_dNrg(iLayer) = 0._rkind

      select case(layerType(iLayer))

        ! ***** soil
        case(iname_soil)

          select case(ixRichards)  ! (form of Richards' equation)
            case(moisture)
              dVolFracLiq_dWat = 1._rkind
              dVolFracIce_dWat = dPsi_dTheta(mLayerVolFracLiq(iLayer),vGn_alpha(iSoil),theta_res(iSoil),theta_sat(iSoil),vGn_n(iSoil),vGn_m(iSoil)) - 1._rkind
            case(mixdform)
              Tcrit = crit_soilT( mLayerMatricHead(iSoil) )
              if(mLayerTemp(iLayer) < Tcrit) then
                dVolFracLiq_dWat = 0._rkind
                dVolFracIce_dWat = dTheta_dPsi(mLayerMatricHead(iSoil),vGn_alpha(iSoil),theta_res(iSoil),theta_sat(iSoil),vGn_n(iSoil),vGn_m(iSoil))
              else
                dVolFracLiq_dWat = dTheta_dPsi(mLayerMatricHead(iSoil),vGn_alpha(iSoil),theta_res(iSoil),theta_sat(iSoil),vGn_n(iSoil),vGn_m(iSoil))
                dVolFracIce_dWat = 0._rkind
              endif
          end select
          dVolFracLiq_dTk = mLayerdTheta_dTk(iLayer) !already zeroed out if not below critical temperature
          dVolFracIce_dTk = -dVolFracLiq_dTk !often can and will simplify one of these terms out

          ! select option for thermal conductivity of soil
          select case(ixThCondSoil)

            ! ** function of soil wetness
            case(funcSoilWet)

              ! compute the thermal conductivity of the wet material (W m-1)
              lambda_wet  = lambda_wetsoil**( 1._rkind - theta_sat(iSoil) ) * lambda_water**theta_sat(iSoil) * lambda_ice**(theta_sat(iSoil) - mLayerVolFracLiq(iLayer))
              dlambda_wet_dWat = -lambda_wet * log(lambda_ice) * dVolFracLiq_dWat
              dlambda_wet_dTk  = -lambda_wet * log(lambda_ice) * dVolFracLiq_dTk

              relativeSat = (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))/theta_sat(iSoil)  ! relative saturation
              ! drelativeSat_dWat = dPsi0_dWat/theta_sat(iLayer), and drelativeSat_dTk = 0 (so dkerstenNum_dTk = 0)
              ! compute the Kersten number (-)
              if(relativeSat > 0.1_rkind)then ! log10(0.1) = -1
                kerstenNum = log10(relativeSat) + 1._rkind
                dkerstenNum_dWat = (dVolFracIce_dWat + dVolFracLiq_dWat) / ( theta_sat(iSoil) * relativeSat * log(10._rkind) )
              else
                kerstenNum = 0._rkind  ! dry thermal conductivity
                dkerstenNum_dWat = 0._rkind
              endif
              ! ...and, compute the thermal conductivity
              mLayerThermalC(iLayer) = kerstenNum*lambda_wet + (1._rkind - kerstenNum)*lambda_drysoil

              ! compute derivatives
              dThermalC_dWat(iLayer) = dkerstenNum_dWat * ( lambda_wet - lambda_drysoil ) + kerstenNum*dlambda_wet_dWat
              dThermalC_dNrg(iLayer) = kerstenNum*dlambda_wet_dTk

            ! ** mixture of constituents
            case(mixConstit)
                mLayerThermalC(iLayer) = thCond_soil(iSoil) * ( 1._rkind - theta_sat(iSoil) ) + & ! soil component
                                         lambda_ice         * mLayerVolFracIce(iLayer)     + & ! ice component
                                         lambda_water       * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                                         lambda_air         * mLayerVolFracAir(iLayer)         ! air component
              ! compute derivatives
              dThermalC_dWat(iLayer) = lambda_ice*dVolFracIce_dWat + lambda_water*dVolFracLiq_dWat + lambda_air*(-dVolFracIce_dWat - dVolFracLiq_dWat)
              dThermalC_dNrg(iLayer) = (lambda_ice - lambda_water) * dVolFracIce_dTk

            ! ** test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
            case(hanssonVZJ)
              fArg  = 1._rkind + f1*mLayerVolFracIce(iLayer)**f2
              xArg  = mLayerVolFracLiq(iLayer) + fArg*mLayerVolFracIce(iLayer)
              dxArg_dWat = dVolFracLiq_dWat + dVolFracIce_dWat * (1._rkind + f1*(f2+1)*mLayerVolFracIce(iLayer)**f2)
              dxArg_dTk  = dVolFracIce_dTk * f1*(f2+1)*mLayerVolFracIce(iLayer)**f2
              ! ...and, compute the thermal conductivity
              mLayerThermalC(iLayer) = c1 + c2*xArg + (c1 - c4)*exp(-(c3*xArg)**c5)

              ! compute derivatives
              dThermalC_dWat(iLayer) = ( c2 - c5*c3*(c3*xArg)**(c5-1)*(c1 - c4)*exp(-(c3*xArg)**c5) ) * dxArg_dWat
              dThermalC_dNrg(iLayer) = ( c2 - c5*c3*(c3*xArg)**(c5-1)*(c1 - c4)*exp(-(c3*xArg)**c5) ) * dxArg_dTk

            ! ** check
            case default; err=20; message=trim(message)//'unable to identify option for thermal conductivity of soil'; return

          end select  ! option for the thermal conductivity of soil

        ! ***** snow
        case(iname_snow)
          dVolFracIce_dWat = ( 1._rkind - mLayerFracLiqSnow(iLayer) )*(iden_water/iden_ice)
          dVolFracIce_dTk = -mLayerdTheta_dTk(iLayer)*(iden_water/iden_ice)

          ! temporally constant thermal conductivity
          if(ixThCondSnow==Smirnova2000)then
            mLayerThermalC(iLayer) = fixedThermalCond_snow
            dThermalC_dWat(iLayer) = 0._rkind
            dThermalC_dNrg(iLayer) = 0._rkind
            ! thermal conductivity as a function of snow density
          else
            call tcond_snow(mLayerVolFracIce(iLayer)*iden_ice,  & ! input: snow density (kg m-3)
                            mLayerThermalC(iLayer),             & ! output: thermal conductivity (W m-1 K-1)
                            err,cmessage)                         ! output: error control
            if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
            select case(ixThCondSnow)
              case(Yen1965)
                dThermalC_dWat(iLayer) = 2._rkind * 3.217d-6 * mLayerVolFracIce(iLayer) * iden_ice**2_i4b * dVolFracIce_dWat
                dThermalC_dNrg(iLayer) = 2._rkind * 3.217d-6 * mLayerVolFracIce(iLayer) * iden_ice**2_i4b * dVolFracIce_dTk
              case(Mellor1977)
                dThermalC_dWat(iLayer) = 2._rkind * 2.576d-6 * mLayerVolFracIce(iLayer) * iden_ice**2_i4b * dVolFracIce_dWat
                dThermalC_dNrg(iLayer) = 2._rkind * 2.576d-6 * mLayerVolFracIce(iLayer) * iden_ice**2_i4b * dVolFracIce_dTk
              case(Jordan1991)
                dThermalC_dWat(iLayer) = ( 7.75d-5 + 2._rkind * 1.105d-6 * mLayerVolFracIce(iLayer) * iden_ice ) * (lambda_ice-lambda_air) * iden_ice * dVolFracIce_dWat
                dThermalC_dNrg(iLayer) = ( 7.75d-5 + 2._rkind * 1.105d-6 * mLayerVolFracIce(iLayer) * iden_ice ) * (lambda_ice-lambda_air) * iden_ice * dVolFracIce_dTk
            end select  ! option for the thermal conductivity of snow
          end if

        ! * error check
        case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute thermal conductivity'; return

      end select

    end do  ! looping through layers

    ! *****
    ! * compute the thermal conductivity of snow at the interface of each layer...
    ! ****************************************************************************
    do iLayer=1,nLayers-1  ! (loop through layers)
    ! get temporary variables
      TCn = mLayerThermalC(iLayer)    ! thermal conductivity below the layer interface (W m-1 K-1)
      TCp = mLayerThermalC(iLayer+1)  ! thermal conductivity above the layer interface (W m-1 K-1)
      zdn = iLayerHeight(iLayer)   - mLayerHeight(iLayer) ! height difference between interface and lower value (m)
      zdp = mLayerHeight(iLayer+1) - iLayerHeight(iLayer) ! height difference between interface and upper value (m)
      den = TCn*zdp + TCp*zdn  ! denominator
      ! compute thermal conductivity
      if(TCn+TCp > epsilon(TCn))then
        iLayerThermalC(iLayer) = (TCn*TCp*(zdn + zdp)) / den
        dThermalC_dWatBelow(iLayer) = ( TCn*(zdn + zdp) - iLayerThermalC(iLayer)*zdn ) / den * dThermalC_dWat(iLayer+1)
        dThermalC_dWatAbove(iLayer) = ( TCp*(zdn + zdp) - iLayerThermalC(iLayer)*zdp ) / den * dThermalC_dWat(iLayer)
        dThermalC_dTempBelow(iLayer) = ( TCn*(zdn + zdp) - iLayerThermalC(iLayer)*zdn ) / den * dThermalC_dNrg(iLayer+1)
        dThermalC_dTempAbove(iLayer) = ( TCp*(zdn + zdp) - iLayerThermalC(iLayer)*zdp ) / den * dThermalC_dNrg(iLayer)
      else
        iLayerThermalC(iLayer) = (TCn*zdn +  TCp*zdp) / (zdn + zdp)
        dThermalC_dWatBelow(iLayer) = zdp / (zdn + zdp) * dThermalC_dWat(iLayer+1)
        dThermalC_dWatAbove(iLayer) = zdn / (zdn + zdp) * dThermalC_dWat(iLayer)
        dThermalC_dTempBelow(iLayer) = zdp / (zdn + zdp) * dThermalC_dNrg(iLayer+1)
        dThermalC_dTempAbove(iLayer) = zdn / (zdn + zdp) * dThermalC_dNrg(iLayer)
      endif
    end do  ! looping through layers

    ! special case of hansson
    if(ixThCondSoil==hanssonVZJ)then
      iLayerThermalC(0) = 28._rkind*(0.5_rkind*(iLayerHeight(1) - iLayerHeight(0)))
      dThermalC_dWatBelow(0) = 0._rkind
      dThermalC_dTempBelow(0) = 0._rkind
    else
      iLayerThermalC(0) = mLayerThermalC(1)
      dThermalC_dWatBelow(0) = dThermalC_dWat(1)
      dThermalC_dTempBelow(0) = dThermalC_dNrg(1)
    end if
    dThermalC_dWatAbove(0) = realMissing
    dThermalC_dTempAbove(0) = realMissing

    ! assume the thermal conductivity at the domain boundaries is equal to the thermal conductivity of the layer
    iLayerThermalC(nLayers) = mLayerThermalC(nLayers)
    dThermalC_dWatAbove(nLayers) = dThermalC_dWat(nLayers)
    dThermalC_dTempAbove(nLayers) = dThermalC_dNrg(nLayers)
    dThermalC_dWatBelow(nLayers) = realMissing
    dThermalC_dTempBelow(nLayers) = realMissing

  ! end association to variables in the data structure
  end associate

end subroutine computThermConduct


end module computThermConduct_module
