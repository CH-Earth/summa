module noahmp_globals

  ! Maybe most of these can be moved to a REDPRM use statement?
  ! MPC -- yes, all of these variables can be local to REDPRM (see additional comments)
  use nrtype
  use module_sf_noahlsm, only: &
       &                       SLCATS,     &   !                               -- only used in REDPRM: number of soil categories
       &                       LUCATS,     &   !                               -- only used in REDPRM: number of land use categories
       &                       CSOIL_DATA, &   ! CSOIL  = CSOIL_DATA           -- only used in REDPRM: soil heat capacity [J M-3 K-1]
       &                       BB,         &   ! BEXP   = BB(SOILTYP)          -- only used in REDPRM: Clapp-Hornberger exponent
       &                       SATDK,      &   ! DKSAT  = SATDK(SOILTYP)       -- only used in REDPRM: saturated soil hydraulic conductivity
       &                       SATDW,      &   ! DWSAT  = SATDW(SOILTYP)       -- only used in REDPRM: saturated soil hydraulic diffusivity
       &                       F11,        &   ! F1     = F11(SOILTYP)         -- only used in REDPRM: soil thermal diffusivity/conductivity coef
       &                       SATPSI,     &   ! PSISAT = SATPSI(SOILTYP)      -- only used in REDPRM: saturated soil matric potential
       &                       QTZ,        &   ! QUARTZ = QTZ(SOILTYP)         -- only used in REDPRM: soil quartz content
       &                       DRYSMC,     &   ! SMCDRY = DRYSMC(SOILTYP)      -- only used in REDPRM: dry soil moisture threshold where direct evap from top layer ends
       &                       MAXSMC,     &   ! SMCMAX = MAXSMC(SOILTYP)      -- only used in REDPRM: porosity, saturated value of soil moisture (volumetric)
       &                       REFSMC,     &   ! SMCREF = REFSMC(SOILTYP)      -- only used in REDPRM: reference soil moisture (field capacity)
       &                       WLTSMC,     &   ! SMCWLT = WLTSMC(SOILTYP)      -- only used in REDPRM: wilting point soil moisture (volumetric)
       &                       RSTBL,      &   ! RSMIN  = RSTBL(VEGTYP)        -- only used in REDPRM: minimum Canopy Resistance [s/m]
       &                       RGLTBL,     &   ! RGL    = RGLTBL(VEGTYP)       -- only used in REDPRM: parameter used in radiation stress function
       &                       HSTBL,      &   ! HS     = HSTBL(VEGTYP)        -- only used in REDPRM: parameter used in vapor pressure deficit function
       &                       NROTBL,     &   ! NROOT  = NROTBL(VEGTYP)       -- only used in REDPRM: rooting depth [as the number of layers]
       &                       TOPT_DATA,  &   ! TOPT   = TOPT_DATA            -- only used in REDPRM: optimum transpiration air temperature
       &                       RSMAX_DATA, &   ! RSMAX  = RSMAX_DATA           -- only used in REDPRM: maximum stomatal resistance
       &                       ZBOT_DATA,  &   ! ZBOT   = ZBOT_DATA            -- only used in REDPRM: Depth (m) of lower boundary soil temperature
       &                       CZIL_DATA,  &   ! CZIL   = CZIL_DATA            -- only used in REDPRM: Calculate roughness length of heat
       &                       FRZK_DATA,  &   ! FRZK   = FRZK_DATA            -- only used in REDPRM: ** Used to get FRZX, to compute maximum infiltration rate (in INFIL)
       &                       SLOPE_DATA, &   ! SLOPE  = SLOPE_DATA(SLOPETYP) -- only used in REDPRM: slope index (0 - 1)
       &                       REFDK_DATA, &   ! REFDK  = REFDK_DATA           -- only used in REDPRM: ** Used to get KDT, to compute maximum infiltration rate (in INFIL)
       &                       REFKDT_DATA     ! REFKDT = REFKDT_DATA          -- only used in REDPRM: ** Used to get KDT, to compute maximum infiltration rate (in INFIL)

  implicit none

! ==================================================================================================
!------------------------------------------------------------------------------------------!
! Physical Constants:                                                                      !
!------------------------------------------------------------------------------------------!

  REAL(DP), PARAMETER :: GRAV   = 9.80616d0   !acceleration due to gravity (m/s2)
  REAL(DP), PARAMETER :: SB     = 5.67E-08  !Stefan-Boltzmann constant (w/m2/k4)
  REAL(DP), PARAMETER :: VKC    = 0.40d0      !von Karman constant
  REAL(DP), PARAMETER :: TFRZ   = 273.16d0    !freezing/melting point (k)
  REAL(DP), PARAMETER :: HSUB   = 2.8440E06 !latent heat of sublimation (j/kg)
  REAL(DP), PARAMETER :: HVAP   = 2.5104E06 !latent heat of vaporization (j/kg)
  REAL(DP), PARAMETER :: HFUS   = 0.3336E06 !latent heat of fusion (j/kg)
  REAL(DP), PARAMETER :: CWAT   = 4.188E06  !specific heat capacity of water (j/m3/k)
  REAL(DP), PARAMETER :: CICE   = 2.094E06  !specific heat capacity of ice (j/m3/k)
  REAL(DP), PARAMETER :: CPAIR  = 1004.64d0   !heat capacity dry air at const pres (j/kg/k)
  REAL(DP), PARAMETER :: TKWAT  = 0.6d0       !thermal conductivity of water (w/m/k)
  REAL(DP), PARAMETER :: TKICE  = 2.2d0       !thermal conductivity of ice (w/m/k)
  REAL(DP), PARAMETER :: TKAIR  = 0.023d0     !thermal conductivity of air (w/m/k)
  REAL(DP), PARAMETER :: RAIR   = 287.04d0    !gas constant for dry air (j/kg/k)
  REAL(DP), PARAMETER :: RW     = 461.269d0   !gas constant for  water vapor (j/kg/k)
  REAL(DP), PARAMETER :: DENH2O = 1000.     !density of water (kg/m3)
  REAL(DP), PARAMETER :: DENICE = 917.      !density of ice (kg/m3)

!------------------------------------------------------------------------------------------!
! From the VEGPARM.TBL tables, as functions of vegetation category.
!------------------------------------------------------------------------------------------!
  INTEGER :: NROOT        !rooting depth [as the number of layers] ( Assigned in REDPRM )
  REAL(DP)    :: RGL          !parameter used in radiation stress function ( Assigned in REDPRM )
  REAL(DP)    :: RSMIN        !minimum Canopy Resistance [s/m] ( Assigned in REDPRM )
  REAL(DP)    :: HS           !parameter used in vapor pressure deficit function ( Assigned in REDPRM )
  REAL(DP)    :: RSMAX        !maximum stomatal resistance ( Assigned in REDPRM )
  REAL(DP)    :: TOPT         !optimum transpiration air temperature.

  ! MPC change: make variables private for a given thread
  !$omp threadprivate(NROOT, RGL, RSMIN, HS, RSMAX, TOPT)

!------------------------------------------------------------------------------------------!
! From the SOILPARM.TBL tables, as functions of soil category.
!------------------------------------------------------------------------------------------!
  REAL(DP)    :: BEXP         !B parameter ( Assigned in REDPRM )
  REAL(DP)    :: SMCDRY       !dry soil moisture threshold where direct evap from top
                          !layer ends (volumetric) ( Assigned in REDPRM )
  REAL(DP)    :: F1           !soil thermal diffusivity/conductivity coef ( Assigned in REDPRM )
  REAL(DP)    :: SMCMAX       !porosity, saturated value of soil moisture (volumetric)
  REAL(DP)    :: SMCREF       !reference soil moisture (field capacity) (volumetric) ( Assigned in REDPRM )
  REAL(DP)    :: PSISAT       !saturated soil matric potential ( Assigned in REDPRM )
  REAL(DP)    :: DKSAT        !saturated soil hydraulic conductivity ( Assigned in REDPRM )
  REAL(DP)    :: DWSAT        !saturated soil hydraulic diffusivity ( Assigned in REDPRM )
  REAL(DP)    :: SMCWLT       !wilting point soil moisture (volumetric) ( Assigned in REDPRM )
  REAL(DP)    :: QUARTZ       !soil quartz content ( Assigned in REDPRM )

  ! MPC change: make variables private for a given thread
  !$omp threadprivate(BEXP, SMCDRY, F1, SMCMAX, SMCREF, PSISAT, DKSAT, DWSAT, SMCWLT, QUARTZ)

!------------------------------------------------------------------------------------------!
! From the GENPARM.TBL file
!------------------------------------------------------------------------------------------!
  REAL(DP)    :: SLOPE       !slope index (0 - 1) ( Assigned in REDPRM )
  REAL(DP)    :: CSOIL       !vol. soil heat capacity [j/m3/K] ( Assigned in REDPRM )
  REAL(DP)    :: ZBOT        !Depth (m) of lower boundary soil temperature ( Assigned in REDPRM )
  REAL(DP)    :: CZIL        !Calculate roughness length of heat ( Assigned in REDPRM )

  ! MPC note: FRZK_DATA, REFDK_DATA, and REFKDT_DATA are used in REDPRM to compute KDT and FRZX
  !           (FRZK, REFDK, and REFKDT are local variables within REDPRM and do not need to be thread private)

  REAL(DP)    :: KDT         !used in compute maximum infiltration rate (in INFIL) ( Assigned in REDPRM )
  REAL(DP)    :: FRZX        !used in compute maximum infiltration rate (in INFIL) ( Assigned in REDPRM )

  ! MPC change: make variables private for a given thread
  !$omp threadprivate(SLOPE, CSOIL, ZBOT, CZIL, KDT, FRZX)


! LSM GENERAL PARAMETERS
! =====================================options for different schemes================================
! options for dynamic vegetation:
! 1 -> off (use table LAI; use FVEG = SHDFAC from input)
! 2 -> on (together with OPT_CRS = 1)
! 3 -> off (use table LAI; calculate FVEG)
! 4 -> off (use table LAI; use maximum vegetation fraction)

  INTEGER :: DVEG    != 2   !

! options for canopy stomatal resistance
! 1-> Ball-Berry; 2->Jarvis

  INTEGER :: OPT_CRS != 1    !(must 1 when DVEG = 2)

! options for soil moisture factor for stomatal resistance
! 1-> Noah (soil moisture)
! 2-> CLM  (matric potential)
! 3-> SSiB (matric potential)

  INTEGER :: OPT_BTR != 1    !(suggested 1)

! options for runoff and groundwater
! 1 -> TOPMODEL with groundwater (Niu et al. 2007 JGR) ;
! 2 -> TOPMODEL with an equilibrium water table (Niu et al. 2005 JGR) ;
! 3 -> original surface and subsurface runoff (free drainage)
! 4 -> BATS surface and subsurface runoff (free drainage)

  INTEGER :: OPT_RUN != 1    !(suggested 1)

! options for surface layer drag coeff (CH & CM)
! 1->M-O ; 2->original Noah (Chen97); 3->MYJ consistent; 4->YSU consistent.

  INTEGER :: OPT_SFC != 1    !(1 or 2 or 3 or 4)

! options for supercooled liquid water (or ice fraction)
! 1-> no iteration (Niu and Yang, 2006 JHM); 2: Koren's iteration

  INTEGER :: OPT_FRZ != 1    !(1 or 2)

! options for frozen soil permeability
! 1 -> linear effects, more permeable (Niu and Yang, 2006, JHM)
! 2 -> nonlinear effects, less permeable (old)

  INTEGER :: OPT_INF != 1    !(suggested 1)

! options for radiation transfer
! 1 -> modified two-stream (gap = F(solar angle, 3D structure ...)<1-FVEG)
! 2 -> two-stream applied to grid-cell (gap = 0)
! 3 -> two-stream applied to vegetated fraction (gap=1-FVEG)

  INTEGER :: OPT_RAD != 1    !(suggested 1)

! options for ground snow surface albedo
! 1-> BATS; 2 -> CLASS

  INTEGER :: OPT_ALB != 2    !(suggested 2)

! options for partitioning  precipitation into rainfall & snowfall
! 1 -> Jordan (1991); 2 -> BATS: when SFCTMP<TFRZ+2.2d0 ; 3-> SFCTMP<TFRZ

  INTEGER :: OPT_SNF != 1    !(suggested 1)

! options for lower boundary condition of soil temperature
! 1 -> zero heat flux from bottom (ZBOT and TBOT not used)
! 2 -> TBOT at ZBOT (8m) read from a file (original Noah)

  INTEGER :: OPT_TBOT != 2   !(suggested 2)

! options for snow/soil temperature time scheme (only layer 1)
! 1 -> semi-implicit; 2 -> full implicit (original Noah)

  INTEGER :: OPT_STC != 1    !(suggested 1)
! ==================================================================================================
! runoff parameters used for SIMTOP and SIMGM:
  REAL(DP), PARAMETER :: TIMEAN = 10.5d0   !gridcell mean topgraphic index (global mean)
  REAL(DP), PARAMETER :: FSATMX = 0.38d0   !maximum surface saturated fraction (global mean)

! adjustable parameters for snow processes

  REAL(DP), PARAMETER :: M      = 1.0d0 ! 2.50d0   !melting factor (-)
  REAL(DP), PARAMETER :: Z0SNO  = 0.002d0  !snow surface roughness length (m) (0.002d0)
  REAL(DP), PARAMETER :: SSI    = 0.03d0   !liquid water holding capacity for snowpack (m3/m3) (0.03d0)
  REAL(DP), PARAMETER :: SWEMX  = 1.00d0   !new snow mass to fully cover old snow (mm)
                                     !equivalent to 10mm depth (density = 100 kg/m3)

! NOTES: things to add or improve
! 1. lake model: explicit representation of lake water storage, sunlight through lake
!    with different purity, turbulent mixing of surface laker water, snow on frozen lake, etc.
! 2. shallow snow wihtout a layer: melting energy
! 3. urban model to be added.
! 4. irrigation
!------------------------------------------------------------------------------------------!
END MODULE NOAHMP_GLOBALS
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
MODULE NOAHMP_VEG_PARAMETERS
    use nrtype

    IMPLICIT NONE

    INTEGER, PARAMETER :: MAX_VEG_PARAMS = 33
    INTEGER, PARAMETER :: MVT   = 27
    INTEGER, PARAMETER :: MBAND = 2

    INTEGER, PRIVATE :: ISURBAN
    INTEGER :: ISWATER
    INTEGER :: ISBARREN
    INTEGER :: ISSNOW
    INTEGER :: EBLFOREST

    REAL(DP) :: CH2OP(MVT)       !maximum intercepted h2o per unit lai+sai (mm)
    REAL(DP) :: DLEAF(MVT)       !characteristic leaf dimension (m)
    REAL(DP) :: Z0MVT(MVT)       !momentum roughness length (m)
    REAL(DP) :: HVT(MVT)         !top of canopy (m)
    REAL(DP) :: HVB(MVT)         !bottom of canopy (m)
    REAL(DP) :: DEN(MVT)         !tree density (no. of trunks per m2)
    REAL(DP) :: RC(MVT)          !tree crown radius (m)
    REAL(DP) :: SAIM(MVT,12)     !monthly stem area index, one-sided
    REAL(DP) :: LAIM(MVT,12)     !monthly leaf area index, one-sided
    REAL(DP) :: SLA(MVT)         !single-side leaf area per Kg [m2/kg]
    REAL(DP) :: DILEFC(MVT)      !coeficient for leaf stress death [1/s]
    REAL(DP) :: DILEFW(MVT)      !coeficient for leaf stress death [1/s]
    REAL(DP) :: FRAGR(MVT)       !fraction of growth respiration  !original was 0.3d0
    REAL(DP) :: LTOVRC(MVT)      !leaf turnover [1/s]

    REAL(DP) :: C3PSN(MVT)       !photosynthetic pathway: 0. = c4, 1. = c3
    REAL(DP) :: KC25(MVT)        !co2 michaelis-menten constant at 25c (pa)
    REAL(DP) :: AKC(MVT)         !q10 for kc25
    REAL(DP) :: KO25(MVT)        !o2 michaelis-menten constant at 25c (pa)
    REAL(DP) :: AKO(MVT)         !q10 for ko25
    REAL(DP) :: VCMX25(MVT)      !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    REAL(DP) :: AVCMX(MVT)       !q10 for vcmx25
    REAL(DP) :: BP(MVT)          !minimum leaf conductance (umol/m**2/s)
    REAL(DP) :: MP(MVT)          !slope of conductance-to-photosynthesis relationship
    REAL(DP) :: QE25(MVT)        !quantum efficiency at 25c (umol co2 / umol photon)
    REAL(DP) :: AQE(MVT)         !q10 for qe25
    REAL(DP) :: RMF25(MVT)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
    REAL(DP) :: RMS25(MVT)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
    REAL(DP) :: RMR25(MVT)       !root maintenance respiration at 25c (umol co2/kg bio/s)
    REAL(DP) :: ARM(MVT)         !q10 for maintenance respiration
    REAL(DP) :: FOLNMX(MVT)      !foliage nitrogen concentration when f(n)=1 (%)
    REAL(DP) :: TMIN(MVT)        !minimum temperature for photosynthesis (k)

    REAL(DP) :: XL(MVT)          !leaf/stem orientation index
    REAL(DP) :: RHOL(MVT,MBAND)  !leaf reflectance: 1=vis, 2=nir
    REAL(DP) :: RHOS(MVT,MBAND)  !stem reflectance: 1=vis, 2=nir
    REAL(DP) :: TAUL(MVT,MBAND)  !leaf transmittance: 1=vis, 2=nir
    REAL(DP) :: TAUS(MVT,MBAND)  !stem transmittance: 1=vis, 2=nir

    REAL(DP) :: MRP(MVT)         !microbial respiration parameter (umol co2 /kg c/ s)
    REAL(DP) :: CWPVT(MVT)       !empirical canopy wind parameter

    REAL(DP) :: WRRAT(MVT)       !wood to non-wood ratio
    REAL(DP) :: WDPOOL(MVT)      !wood pool (switch 1 or 0) depending on woody or not [-]
    REAL(DP) :: TDLEF(MVT)       !characteristic T for leaf freezing [K]

    INTEGER :: IK,IM
    REAL(DP) :: TMP10(MVT*MBAND)
    REAL(DP) :: TMP11(MVT*MBAND)
    REAL(DP) :: TMP12(MVT*MBAND)
    REAL(DP) :: TMP13(MVT*MBAND)
    REAL(DP) :: TMP14(MVT*12)
    REAL(DP) :: TMP15(MVT*12)
    REAL(DP) :: TMP16(MVT*5)

    real(DP) slarea(MVT)
    real(DP) eps(MVT,5)

CONTAINS
  subroutine read_mp_veg_parameters(FILENAME_VEGTABLE,DATASET_IDENTIFIER)
    implicit none
    character(len=*), intent(in) :: FILENAME_VEGTABLE
    character(len=*), intent(in) :: DATASET_IDENTIFIER
    integer :: ierr

    ! Temporary arrays used in reshaping namelist arrays
    REAL(DP) :: TMP10(MVT*MBAND)
    REAL(DP) :: TMP11(MVT*MBAND)
    REAL(DP) :: TMP12(MVT*MBAND)
    REAL(DP) :: TMP13(MVT*MBAND)
    REAL(DP) :: TMP14(MVT*12)
    REAL(DP) :: TMP15(MVT*12)
    REAL(DP) :: TMP16(MVT*5)

    integer :: NVEG
    character(len=256) :: VEG_DATASET_DESCRIPTION

    NAMELIST / noah_mp_usgs_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_usgs_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_modis_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_modis_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    ! MPC change: enable use of alternative veg tables
    !  - in this case, tables using attributes from other models used in the PLUMBER experiment

    NAMELIST / noah_mp_plumberCABLE_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberCABLE_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_plumberCHTESSEL_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberCHTESSEL_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_plumberSUMMA_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberSUMMA_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
    CH2OP  = -1.E36
    DLEAF  = -1.E36
    Z0MVT  = -1.E36
    HVT    = -1.E36
    HVB    = -1.E36
    DEN    = -1.E36
    RC     = -1.E36
    RHOL   = -1.E36
    RHOS   = -1.E36
    TAUL   = -1.E36
    TAUS   = -1.E36
    XL     = -1.E36
    CWPVT  = -1.E36
    C3PSN  = -1.E36
    KC25   = -1.E36
    AKC    = -1.E36
    KO25   = -1.E36
    AKO    = -1.E36
    AVCMX  = -1.E36
    AQE    = -1.E36
    LTOVRC = -1.E36
    DILEFC = -1.E36
    DILEFW = -1.E36
    RMF25  = -1.E36
    SLA    = -1.E36
    FRAGR  = -1.E36
    TMIN   = -1.E36
    VCMX25 = -1.E36
    TDLEF  = -1.E36
    BP     = -1.E36
    MP     = -1.E36
    QE25   = -1.E36
    RMS25  = -1.E36
    RMR25  = -1.E36
    ARM    = -1.E36
    FOLNMX = -1.E36
    WDPOOL = -1.E36
    WRRAT  = -1.E36
    MRP    = -1.E36
    SAIM   = -1.E36
    LAIM   = -1.E36
    SLAREA = -1.E36
    EPS    = -1.E36

    open(15, file=trim(FILENAME_VEGTABLE), status='old', form='formatted', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,'("****** Error ******************************************************")')
       write(*,'("Cannot find file MPTABLE.TBL")')
       write(*,'("STOP")')
       write(*,'("*******************************************************************")')
       call wrf_error_fatal("STOP in Noah-MP read_mp_veg_parameters")
    endif

    if ( trim(DATASET_IDENTIFIER) == "USGS" ) then
       read(15,noah_mp_usgs_veg_categories)
       read(15,noah_mp_usgs_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "MODIFIED_IGBP_MODIS_NOAH" ) then
       read(15,noah_mp_modis_veg_categories)
       read(15,noah_mp_modis_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberCABLE" ) then
       read(15,noah_mp_plumberCABLE_veg_categories)
       read(15,noah_mp_plumberCABLE_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberCHTESSEL" ) then
       read(15,noah_mp_plumberCHTESSEL_veg_categories)
       read(15,noah_mp_plumberCHTESSEL_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberSUMMA" ) then
       read(15,noah_mp_plumberSUMMA_veg_categories)
       read(15,noah_mp_plumberSUMMA_parameters)
    else
       write(*,'("Unrecognized DATASET_IDENTIFIER in subroutine READ_MP_VEG_PARAMETERS")')
       write(*,'("DATASET_IDENTIFIER = ''", A, "''")') trim(DATASET_IDENTIFIER)
       call wrf_error_fatal("STOP in Noah-MP read_mp_veg_parameters")
    endif
    close(15)

    ! Problem.  Namelist reading of 2-d arrays doesn't work well when the arrays are declared with larger dimension than the
    ! variables in the provided namelist.  So we need to reshape the 2-d arrays after we've read them.

    if ( MVT > NVEG ) then

       !
       ! Reshape the 2-d arrays:
       !

       TMP10 = reshape( RHOL, (/ MVT*size(RHOL,2) /))
       TMP11 = reshape( RHOS, (/ MVT*size(RHOS,2) /))
       TMP12 = reshape( TAUL, (/ MVT*size(TAUL,2) /))
       TMP13 = reshape( TAUS, (/ MVT*size(TAUS,2) /))
       TMP14 = reshape( SAIM, (/ MVT*size(SAIM,2) /))
       TMP15 = reshape( LAIM, (/ MVT*size(LAIM,2) /))
       TMP16 = reshape( EPS,  (/ MVT*size(EPS ,2) /))

       RHOL(1:NVEG,:) = reshape( TMP10, (/ NVEG, size(RHOL,2) /))
       RHOS(1:NVEG,:) = reshape( TMP11, (/ NVEG, size(RHOS,2) /))
       TAUL(1:NVEG,:) = reshape( TMP12, (/ NVEG, size(TAUL,2) /))
       TAUS(1:NVEG,:) = reshape( TMP13, (/ NVEG, size(TAUS,2) /))
       SAIM(1:NVEG,:) = reshape( TMP14, (/ NVEG, size(SAIM,2) /))
       LAIM(1:NVEG,:) = reshape( TMP15, (/ NVEG, size(LAIM,2) /))
       EPS(1:NVEG,:)  = reshape( TMP16, (/ NVEG, size(EPS,2)  /))

       RHOL(NVEG+1:MVT,:) = -1.E36
       RHOS(NVEG+1:MVT,:) = -1.E36
       TAUL(NVEG+1:MVT,:) = -1.E36
       TAUS(NVEG+1:MVT,:) = -1.E36
       SAIM(NVEG+1:MVT,:) = -1.E36
       LAIM(NVEG+1:MVT,:) = -1.E36
       EPS( NVEG+1:MVT,:) = -1.E36
    endif

  end subroutine read_mp_veg_parameters

END MODULE NOAHMP_VEG_PARAMETERS
! ==================================================================================================
! ==================================================================================================
MODULE NOAHMP_RAD_PARAMETERS
    use nrtype

    IMPLICIT NONE

    INTEGER I                ! loop index
    INTEGER, PARAMETER :: MSC   = 9
    INTEGER, PARAMETER :: MBAND = 2

    REAL(DP) :: ALBSAT(MSC,MBAND)   !saturated soil albedos: 1=vis, 2=nir
    REAL(DP) :: ALBDRY(MSC,MBAND)   !dry soil albedos: 1=vis, 2=nir
    REAL(DP) :: ALBICE(MBAND)       !albedo land ice: 1=vis, 2=nir
    REAL(DP) :: ALBLAK(MBAND)       !albedo frozen lakes: 1=vis, 2=nir
    REAL(DP) :: OMEGAS(MBAND)       !two-stream parameter omega for snow
    REAL(DP) :: BETADS              !two-stream parameter betad for snow
    REAL(DP) :: BETAIS              !two-stream parameter betad for snow
    REAL(DP) :: EG(2)               !emissivity

! saturated soil albedos: 1=vis, 2=nir
    DATA(ALBSAT(I,1),I=1,8)/0.15d0,0.11d0,0.10d0,0.09d0,0.08d0,0.07d0,0.06d0,0.05d0/
    DATA(ALBSAT(I,2),I=1,8)/0.30d0,0.22d0,0.20d0,0.18d0,0.16d0,0.14d0,0.12d0,0.10d0/

! dry soil albedos: 1=vis, 2=nir
    DATA(ALBDRY(I,1),I=1,8)/0.27d0,0.22d0,0.20d0,0.18d0,0.16d0,0.14d0,0.12d0,0.10d0/
    DATA(ALBDRY(I,2),I=1,8)/0.54d0,0.44d0,0.40d0,0.36d0,0.32d0,0.28d0,0.24d0,0.20d0/

! albedo land ice: 1=vis, 2=nir
    DATA (ALBICE(I),I=1,MBAND) /0.80d0, 0.55d0/

! albedo frozen lakes: 1=vis, 2=nir
    DATA (ALBLAK(I),I=1,MBAND) /0.60d0, 0.40d0/

! omega,betad,betai for snow
    DATA (OMEGAS(I),I=1,MBAND) /0.8d0, 0.4d0/
    DATA BETADS, BETAIS /0.5d0, 0.5d0/

! emissivity ground surface
      DATA EG /0.97d0, 0.98d0/ ! 1-soil;2-lake

END MODULE NOAHMP_RAD_PARAMETERS
! ==================================================================================================

MODULE NOAHMP_ROUTINES
  USE NOAHMP_GLOBALS
  use nrtype
  IMPLICIT NONE

  public  :: noahmp_options
  public  :: REDPRM

  public :: PHENOLOGY
  public ::  RADIATION
  private :: ALBEDO
  private :: SNOW_AGE
  private :: SNOWALB_BATS
  private :: SNOWALB_CLASS
  private :: GROUNDALB
  public ::  TWOSTREAM
  public ::  STOMATA
  public ::  CANRES

contains
!
! ==================================================================================================
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  SUBROUTINE PHENOLOGY (VEGTYP , ISURBAN, SNOWH  , TV     , LAT   , YEARLEN , JULIAN , & !in
                        LAI    , SAI    , TROOT  , HTOP  , ELAI    , ESAI   , IGS)

! --------------------------------------------------------------------------------------------------
! vegetation phenology considering vegeation canopy being buries by snow and evolution in time
! --------------------------------------------------------------------------------------------------
  USE NOAHMP_VEG_PARAMETERS
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs
  INTEGER                , INTENT(IN   ) :: VEGTYP !vegetation type
  INTEGER                , INTENT(IN   ) :: ISURBAN!urban category
  REAL(DP)                   , INTENT(IN   ) :: SNOWH  !snow height [m]
  REAL(DP)                   , INTENT(IN   ) :: TV     !vegetation temperature (k)
  REAL(DP)                   , INTENT(IN   ) :: LAT    !latitude (radians)
  INTEGER                , INTENT(IN   ) :: YEARLEN!Number of days in the particular year
  REAL(DP)                   , INTENT(IN   ) :: JULIAN !Julian day of year (fractional) ( 0 <= JULIAN < YEARLEN )
  REAL(DP)                   , INTENT(IN   ) :: TROOT  !root-zone averaged temperature (k)
  REAL(DP)                   , INTENT(INOUT) :: LAI    !LAI, unadjusted for burying by snow
  REAL(DP)                   , INTENT(INOUT) :: SAI    !SAI, unadjusted for burying by snow

! outputs
  REAL(DP)                   , INTENT(OUT  ) :: HTOP   !top of canopy layer (m)
  REAL(DP)                   , INTENT(OUT  ) :: ELAI   !leaf area index, after burying by snow
  REAL(DP)                   , INTENT(OUT  ) :: ESAI   !stem area index, after burying by snow
  REAL(DP)                   , INTENT(OUT  ) :: IGS    !growing season index (0=off, 1=on)

! locals

  REAL(DP)                                   :: DB     !thickness of canopy buried by snow (m)
  REAL(DP)                                   :: FB     !fraction of canopy buried by snow
  REAL(DP)                                   :: SNOWHC !critical snow depth at which short vege
                                                   !is fully covered by snow

  INTEGER                                :: K       !index
  INTEGER                                :: IT1,IT2 !interpolation months
  REAL(DP)                                   :: DAY     !current day of year ( 0 <= DAY < YEARLEN )
  REAL(DP)                                   :: WT1,WT2 !interpolation weights
  REAL(DP)                                   :: T       !current month (1.00d0, ..., 12.00d0)
! --------------------------------------------------------------------------------------------------

  IF ( DVEG == 1 .or. DVEG == 3 .or. DVEG == 4 ) THEN

     IF (LAT >= 0.) THEN
        ! Northern Hemisphere
        DAY = JULIAN
     ELSE
        ! Southern Hemisphere.  DAY is shifted by 1/2 year.
        DAY = MOD ( JULIAN + ( 0.5d0 * YEARLEN ) , DBLE(YEARLEN) )
     ENDIF

     T = 12. * DAY / DBLE(YEARLEN)
     IT1 = T + 0.5d0
     IT2 = IT1 + 1
     WT1 = (IT1+0.5d0) - T
     WT2 = 1.-WT1
     IF (IT1 .LT.  1) IT1 = 12
     IF (IT2 .GT. 12) IT2 = 1

     LAI = WT1*LAIM(VEGTYP,IT1) + WT2*LAIM(VEGTYP,IT2)
     SAI = WT1*SAIM(VEGTYP,IT1) + WT2*SAIM(VEGTYP,IT2)
  ENDIF
  IF (SAI < 0.05d0) SAI = 0.0d0  ! MB: SAI CHECK
  IF (LAI < 0.05d0 .OR. SAI == 0.0d0) LAI = 0.0d0  ! MB: LAI CHECK

  IF ( ( VEGTYP == ISWATER ) .OR. ( VEGTYP == ISBARREN ) .OR. ( VEGTYP == ISSNOW ) .or. ( VEGTYP == ISURBAN) ) THEN
     LAI  = 0.
     SAI  = 0.
  ENDIF

!buried by snow

     DB = MIN( MAX(SNOWH - HVB(VEGTYP),0.), HVT(VEGTYP)-HVB(VEGTYP) )
     FB = DB / MAX(1.E-06,HVT(VEGTYP)-HVB(VEGTYP))
     !print*, 'HVB(VEGTYP), HVT(VEGTYP), DB, FB = ', HVB(VEGTYP), HVT(VEGTYP), DB, FB

     IF(HVT(VEGTYP)> 0. .AND. HVT(VEGTYP) <= 0.5d0) THEN
       SNOWHC = HVT(VEGTYP)*EXP(-SNOWH/0.1d0)
       IF (SNOWH>SNOWHC) THEN
        FB = 1.
       ELSE
        FB = SNOWH/SNOWHC
       ENDIF
     ENDIF

     ELAI =  LAI*(1.-FB)
     ESAI =  SAI*(1.-FB)
     IF (ESAI < 0.05d0) ESAI = 0.0d0  ! MB: ESAI CHECK
     IF (ELAI < 0.05d0 .OR. ESAI == 0.0d0) ELAI = 0.0d0  ! MB: LAI CHECK

     IF (TV .GT. TMIN(VEGTYP)) THEN
         IGS = 1.
     ELSE
         IGS = 0.
     ENDIF

     HTOP = HVT(VEGTYP)

  END SUBROUTINE PHENOLOGY
! ==================================================================================================
! ==================================================================================================
  SUBROUTINE RADIATION (VEGTYP  ,IST     ,ISC     ,ICE     ,NSOIL   , & !in
                        SNEQVO  ,SNEQV   ,DT      ,COSZ    ,SNOWH   , & !in
                        TG      ,TV      ,FSNO    ,QSNOW   ,FWET    , & !in
                        ELAI    ,ESAI    ,SMC     ,SOLAD   ,SOLAI   , & !in
                        FVEG    ,ILOC    ,JLOC    ,                   & !in
                        ALBOLD  ,TAUSS   ,                            & !inout
                        FSUN    ,LAISUN  ,LAISHA  ,PARSUN  ,PARSHA  , & !out
                        SAV     ,SAG     ,FSR     ,FSA     ,FSRV    , &
                        FSRG    ,BGAP    ,WGAP)            !out
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
  INTEGER, INTENT(IN)                  :: ILOC
  INTEGER, INTENT(IN)                  :: JLOC
  INTEGER, INTENT(IN)                  :: VEGTYP !vegetation type
  INTEGER, INTENT(IN)                  :: IST    !surface type
  INTEGER, INTENT(IN)                  :: ISC    !soil color type (1-lighest; 8-darkest)
  INTEGER, INTENT(IN)                  :: ICE    !ice (ice = 1)
  INTEGER, INTENT(IN)                  :: NSOIL  !number of soil layers

  REAL(DP), INTENT(IN)                     :: DT     !time step [s]
  REAL(DP), INTENT(IN)                     :: QSNOW  !snowfall (mm/s)
  REAL(DP), INTENT(IN)                     :: SNEQVO !snow mass at last time step(mm)
  REAL(DP), INTENT(IN)                     :: SNEQV  !snow mass (mm)
  REAL(DP), INTENT(IN)                     :: SNOWH  !snow height (mm)
  REAL(DP), INTENT(IN)                     :: COSZ   !cosine solar zenith angle (0-1)
  REAL(DP), INTENT(IN)                     :: TG     !ground temperature (k)
  REAL(DP), INTENT(IN)                     :: TV     !vegetation temperature (k)
  REAL(DP), INTENT(IN)                     :: ELAI   !LAI, one-sided, adjusted for burying by snow
  REAL(DP), INTENT(IN)                     :: ESAI   !SAI, one-sided, adjusted for burying by snow
  REAL(DP), INTENT(IN)                     :: FWET   !fraction of canopy that is wet
  REAL(DP), DIMENSION(1:NSOIL), INTENT(IN) :: SMC    !volumetric soil water [m3/m3]
  REAL(DP), DIMENSION(1:2)    , INTENT(IN) :: SOLAD  !incoming direct solar radiation (w/m2)
  REAL(DP), DIMENSION(1:2)    , INTENT(IN) :: SOLAI  !incoming diffuse solar radiation (w/m2)
  REAL(DP), INTENT(IN)                     :: FSNO   !snow cover fraction (-)
  REAL(DP), INTENT(IN)                     :: FVEG   !green vegetation fraction [0.0d0-1.0d0]

! inout
  REAL(DP),                  INTENT(INOUT) :: ALBOLD !snow albedo at last time step (CLASS type)
  REAL(DP),                  INTENT(INOUT) :: TAUSS  !non-dimensional snow age.

! output
  REAL(DP), INTENT(OUT)                    :: FSUN   !sunlit fraction of canopy (-)
  REAL(DP), INTENT(OUT)                    :: LAISUN !sunlit leaf area (-)
  REAL(DP), INTENT(OUT)                    :: LAISHA !shaded leaf area (-)
  REAL(DP), INTENT(OUT)                    :: PARSUN !average absorbed par for sunlit leaves (w/m2)
  REAL(DP), INTENT(OUT)                    :: PARSHA !average absorbed par for shaded leaves (w/m2)
  REAL(DP), INTENT(OUT)                    :: SAV    !solar radiation absorbed by vegetation (w/m2)
  REAL(DP), INTENT(OUT)                    :: SAG    !solar radiation absorbed by ground (w/m2)
  REAL(DP), INTENT(OUT)                    :: FSA    !total absorbed solar radiation (w/m2)
  REAL(DP), INTENT(OUT)                    :: FSR    !total reflected solar radiation (w/m2)

!jref:start
  REAL(DP), INTENT(OUT)                    :: FSRV    !veg. reflected solar radiation (w/m2)
  REAL(DP), INTENT(OUT)                    :: FSRG    !ground reflected solar radiation (w/m2)
  REAL(DP), INTENT(OUT)                    :: BGAP
  REAL(DP), INTENT(OUT)                    :: WGAP
!jref:end

! local
  REAL(DP)                                 :: FAGE   !snow age function (0 - new snow)
  REAL(DP), DIMENSION(1:2)                 :: ALBGRD !ground albedo (direct)
  REAL(DP), DIMENSION(1:2)                 :: ALBGRI !ground albedo (diffuse)
  REAL(DP), DIMENSION(1:2)                 :: ALBD   !surface albedo (direct)
  REAL(DP), DIMENSION(1:2)                 :: ALBI   !surface albedo (diffuse)
  REAL(DP), DIMENSION(1:2)                 :: FABD   !flux abs by veg (per unit direct flux)
  REAL(DP), DIMENSION(1:2)                 :: FABI   !flux abs by veg (per unit diffuse flux)
  REAL(DP), DIMENSION(1:2)                 :: FTDD   !down direct flux below veg (per unit dir flux)
  REAL(DP), DIMENSION(1:2)                 :: FTID   !down diffuse flux below veg (per unit dir flux)
  REAL(DP), DIMENSION(1:2)                 :: FTII   !down diffuse flux below veg (per unit dif flux)
!jref:start
  REAL(DP), DIMENSION(1:2)                 :: FREVI
  REAL(DP), DIMENSION(1:2)                 :: FREVD
  REAL(DP), DIMENSION(1:2)                 :: FREGI
  REAL(DP), DIMENSION(1:2)                 :: FREGD
!jref:end

  REAL(DP)                                 :: FSHA   !shaded fraction of canopy
  REAL(DP)                                 :: VAI    !total LAI + stem area index, one sided

  REAL(DP),PARAMETER :: MPE = 1.E-6
  LOGICAL VEG  !true: vegetated for surface temperature calculation

! --------------------------------------------------------------------------------------------------

! surface abeldo

   CALL ALBEDO (VEGTYP ,IST    ,ISC    ,ICE    ,NSOIL  , & !in
                DT     ,COSZ   ,FAGE   ,ELAI   ,ESAI   , & !in
                TG     ,TV     ,SNOWH  ,FSNO   ,FWET   , & !in
                SMC    ,SNEQVO ,SNEQV  ,QSNOW  ,FVEG   , & !in
                ILOC   ,JLOC   ,                         & !in
                ALBOLD ,TAUSS                          , & !inout
                ALBGRD ,ALBGRI ,ALBD   ,ALBI   ,FABD   , & !out
                FABI   ,FTDD   ,FTID   ,FTII   ,FSUN   , & !)   !out
                FREVI  ,FREVD   ,FREGD ,FREGI  ,BGAP   , & !inout
                WGAP)

! surface radiation

     FSHA = 1.-FSUN
     LAISUN = ELAI*FSUN
     LAISHA = ELAI*FSHA
     VAI = ELAI+ ESAI
     IF (VAI .GT. 0.) THEN
        VEG = .TRUE.
     ELSE
        VEG = .FALSE.
     END IF

   CALL SURRAD (MPE    ,FSUN   ,FSHA   ,ELAI   ,VAI    , & !in
                LAISUN ,LAISHA ,SOLAD  ,SOLAI  ,FABD   , & !in
                FABI   ,FTDD   ,FTID   ,FTII   ,ALBGRD , & !in
                ALBGRI ,ALBD   ,ALBI   ,ILOC   ,JLOC   , & !in
                PARSUN ,PARSHA ,SAV    ,SAG    ,FSA    , & !out
                FSR    ,                                 & !out
                FREVI  ,FREVD  ,FREGD  ,FREGI  ,FSRV   , & !inout
                FSRG)

  END SUBROUTINE RADIATION
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  SUBROUTINE ALBEDO (VEGTYP ,IST    ,ISC    ,ICE    ,NSOIL  , & !in
                     DT     ,COSZ   ,FAGE   ,ELAI   ,ESAI   , & !in
                     TG     ,TV     ,SNOWH  ,FSNO   ,FWET   , & !in
                     SMC    ,SNEQVO ,SNEQV  ,QSNOW  ,FVEG   , & !in
                     ILOC   ,JLOC   ,                         & !in
                     ALBOLD ,TAUSS                          , & !inout
                     ALBGRD ,ALBGRI ,ALBD   ,ALBI   ,FABD   , & !out
                     FABI   ,FTDD   ,FTID   ,FTII   ,FSUN   , & !out
                     FREVI  ,FREVD  ,FREGD  ,FREGI  ,BGAP   , & !out
                     WGAP)

! --------------------------------------------------------------------------------------------------
! surface albedos. also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation.
! also sunlit fraction of the canopy.
! --------------------------------------------------------------------------------------------------
  USE NOAHMP_VEG_PARAMETERS
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
  INTEGER,                  INTENT(IN)  :: ILOC
  INTEGER,                  INTENT(IN)  :: JLOC
  INTEGER,                  INTENT(IN)  :: NSOIL  !number of soil layers
  INTEGER,                  INTENT(IN)  :: VEGTYP !vegetation type
  INTEGER,                  INTENT(IN)  :: IST    !surface type
  INTEGER,                  INTENT(IN)  :: ISC    !soil color type (1-lighest; 8-darkest)
  INTEGER,                  INTENT(IN)  :: ICE    !ice (ice = 1)

  REAL(DP),                     INTENT(IN)  :: DT     !time step [sec]
  REAL(DP),                     INTENT(IN)  :: QSNOW  !snowfall
  REAL(DP),                     INTENT(IN)  :: COSZ   !cosine solar zenith angle for next time step
  REAL(DP),                     INTENT(IN)  :: SNOWH  !snow height (mm)
  REAL(DP),                     INTENT(IN)  :: TG     !ground temperature (k)
  REAL(DP),                     INTENT(IN)  :: TV     !vegetation temperature (k)
  REAL(DP),                     INTENT(IN)  :: ELAI   !LAI, one-sided, adjusted for burying by snow
  REAL(DP),                     INTENT(IN)  :: ESAI   !SAI, one-sided, adjusted for burying by snow
  REAL(DP),                     INTENT(IN)  :: FSNO   !fraction of grid covered by snow
  REAL(DP),                     INTENT(IN)  :: FWET   !fraction of canopy that is wet
  REAL(DP),                     INTENT(IN)  :: SNEQVO !snow mass at last time step(mm)
  REAL(DP),                     INTENT(IN)  :: SNEQV  !snow mass (mm)
  REAL(DP),                     INTENT(IN)  :: FVEG   !green vegetation fraction [0.0d0-1.0d0]
  REAL(DP), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    !volumetric soil water (m3/m3)

! inout
  REAL(DP),                  INTENT(INOUT)  :: ALBOLD !snow albedo at last time step (CLASS type)
  REAL(DP),                  INTENT(INOUT)  :: TAUSS  !non-dimensional snow age

! output
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: ALBGRD !ground albedo (direct)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: ALBGRI !ground albedo (diffuse)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: ALBD   !surface albedo (direct)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: ALBI   !surface albedo (diffuse)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FABD   !flux abs by veg (per unit direct flux)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FABI   !flux abs by veg (per unit diffuse flux)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FTDD   !down direct flux below veg (per unit dir flux)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FTID   !down diffuse flux below veg (per unit dir flux)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FTII   !down diffuse flux below veg (per unit dif flux)
  REAL(DP),                     INTENT(OUT) :: FSUN   !sunlit fraction of canopy (-)
!jref:start
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FREVD
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FREVI
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FREGD
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: FREGI
  REAL(DP), INTENT(OUT) :: BGAP
  REAL(DP), INTENT(OUT) :: WGAP
!jref:end

! ------------------------------------------------------------------------
! ------------------------ local variables -------------------------------
! local
  REAL(DP)                 :: FAGE     !snow age function
  REAL(DP)                 :: ALB
  INTEGER              :: IB       !indices
  INTEGER              :: NBAND    !number of solar radiation wave bands
  INTEGER              :: IC       !direct beam: ic=0; diffuse: ic=1

  REAL(DP)                 :: WL       !fraction of LAI+SAI that is LAI
  REAL(DP)                 :: WS       !fraction of LAI+SAI that is SAI
  REAL(DP)                 :: MPE      !prevents overflow for division by zero

  REAL(DP), DIMENSION(1:2) :: RHO      !leaf/stem reflectance weighted by fraction LAI and SAI
  REAL(DP), DIMENSION(1:2) :: TAU      !leaf/stem transmittance weighted by fraction LAI and SAI
  REAL(DP), DIMENSION(1:2) :: FTDI     !down direct flux below veg per unit dif flux = 0
  REAL(DP), DIMENSION(1:2) :: ALBSND   !snow albedo (direct)
  REAL(DP), DIMENSION(1:2) :: ALBSNI   !snow albedo (diffuse)

  REAL(DP)                 :: VAI      !ELAI+ESAI
  REAL(DP)                 :: GDIR     !average projected leaf/stem area in solar direction
  REAL(DP)                 :: EXT      !optical depth direct beam per unit leaf + stem area

! --------------------------------------------------------------------------------------------------

  NBAND = 2
  MPE = 1.E-06
  BGAP = 0.
  WGAP = 0.

! initialize output because solar radiation only done if COSZ > 0

  DO IB = 1, NBAND
    ALBD(IB) = 0.
    ALBI(IB) = 0.
    ALBGRD(IB) = 0.
    ALBGRI(IB) = 0.
    FABD(IB) = 0.
    FABI(IB) = 0.
    FTDD(IB) = 0.
    FTID(IB) = 0.
    FTII(IB) = 0.
    IF (IB.EQ.1) FSUN = 0.
  END DO

  IF(COSZ <= 0) GOTO 100

! weight reflectance/transmittance by LAI and SAI

  DO IB = 1, NBAND
    VAI = ELAI + ESAI
    WL  = ELAI / MAX(VAI,MPE)
    WS  = ESAI / MAX(VAI,MPE)
    RHO(IB) = MAX(RHOL(VEGTYP,IB)*WL+RHOS(VEGTYP,IB)*WS, MPE)
    TAU(IB) = MAX(TAUL(VEGTYP,IB)*WL+TAUS(VEGTYP,IB)*WS, MPE)
  END DO

! snow age

   CALL SNOW_AGE (DT,TG,SNEQVO,SNEQV,TAUSS,FAGE)

! snow albedos: only if COSZ > 0 and FSNO > 0

  IF(OPT_ALB == 1) &
     CALL SNOWALB_BATS (NBAND, FSNO,COSZ,FAGE,ALBSND,ALBSNI)
  IF(OPT_ALB == 2) THEN
     CALL SNOWALB_CLASS (NBAND,QSNOW,DT,ALB,ALBOLD,ALBSND,ALBSNI,ILOC,JLOC)
     ALBOLD = ALB
  END IF

! ground surface albedo

  CALL GROUNDALB (NSOIL   ,NBAND   ,ICE     ,IST     ,ISC     , & !in
                  FSNO    ,SMC     ,ALBSND  ,ALBSNI  ,COSZ    , & !in
                  TG      ,ILOC    ,JLOC    ,                   & !in
                  ALBGRD  ,ALBGRI  )                              !out


! loop over NBAND wavebands to calculate surface albedos and solar
! fluxes for unit incoming direct (IC=0) and diffuse flux (IC=1)

  DO IB = 1, NBAND
      IC = 0      ! direct
      CALL TWOSTREAM (IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
                      FWET   ,TV      ,ALBGRD  ,ALBGRI  ,RHO    , & !in
                      TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
                      FABD   ,ALBD    ,FTDD    ,FTID    ,GDIR   , &!)   !out
                      FREVD  ,FREGD   ,BGAP    ,WGAP)

      IC = 1      ! diffuse
      CALL TWOSTREAM (IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
                      FWET   ,TV      ,ALBGRD  ,ALBGRI  ,RHO    , & !in
                      TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
                      FABI   ,ALBI    ,FTDI    ,FTII    ,GDIR   , & !)   !out
                      FREVI  ,FREGI   ,BGAP    ,WGAP)

  END DO

! sunlit fraction of canopy. set FSUN = 0 if FSUN < 0.01d0.

  EXT = GDIR/COSZ * SQRT(1.-RHO(1)-TAU(1))
  FSUN = (1.-EXP(-EXT*VAI)) / MAX(EXT*VAI,MPE)
  EXT = FSUN

  IF (EXT .LT. 0.01d0) THEN
     WL = 0.
  ELSE
     WL = EXT
  END IF
  FSUN = WL

100 CONTINUE

  END SUBROUTINE ALBEDO
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  SUBROUTINE SURRAD (MPE     ,FSUN    ,FSHA    ,ELAI    ,VAI     , & !in
                     LAISUN  ,LAISHA  ,SOLAD   ,SOLAI   ,FABD    , & !in
                     FABI    ,FTDD    ,FTID    ,FTII    ,ALBGRD  , & !in
                     ALBGRI  ,ALBD    ,ALBI    ,ILOC    ,JLOC    , & !in
                     PARSUN  ,PARSHA  ,SAV     ,SAG     ,FSA     , & !out
                     FSR     , & !)                                       !out
                     FREVI   ,FREVD   ,FREGD   ,FREGI   ,FSRV    , &
                     FSRG) !inout

! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  INTEGER, INTENT(IN)              :: ILOC
  INTEGER, INTENT(IN)              :: JLOC
  REAL(DP), INTENT(IN)                 :: MPE     !prevents underflow errors if division by zero

  REAL(DP), INTENT(IN)                 :: FSUN    !sunlit fraction of canopy
  REAL(DP), INTENT(IN)                 :: FSHA    !shaded fraction of canopy
  REAL(DP), INTENT(IN)                 :: ELAI    !leaf area, one-sided
  REAL(DP), INTENT(IN)                 :: VAI     !leaf + stem area, one-sided
  REAL(DP), INTENT(IN)                 :: LAISUN  !sunlit leaf area index, one-sided
  REAL(DP), INTENT(IN)                 :: LAISHA  !shaded leaf area index, one-sided

  REAL(DP), DIMENSION(1:2), INTENT(IN) :: SOLAD   !incoming direct solar radiation (w/m2)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: SOLAI   !incoming diffuse solar radiation (w/m2)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FABD    !flux abs by veg (per unit incoming direct flux)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FABI    !flux abs by veg (per unit incoming diffuse flux)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FTDD    !down dir flux below veg (per incoming dir flux)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FTID    !down dif flux below veg (per incoming dir flux)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FTII    !down dif flux below veg (per incoming dif flux)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: ALBGRD  !ground albedo (direct)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: ALBGRI  !ground albedo (diffuse)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: ALBD    !overall surface albedo (direct)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: ALBI    !overall surface albedo (diffuse)

  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FREVD    !overall surface albedo veg (direct)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FREVI    !overall surface albedo veg (diffuse)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FREGD    !overall surface albedo grd (direct)
  REAL(DP), DIMENSION(1:2), INTENT(IN) :: FREGI    !overall surface albedo grd (diffuse)

! output

  REAL(DP), INTENT(OUT)                :: PARSUN  !average absorbed par for sunlit leaves (w/m2)
  REAL(DP), INTENT(OUT)                :: PARSHA  !average absorbed par for shaded leaves (w/m2)
  REAL(DP), INTENT(OUT)                :: SAV     !solar radiation absorbed by vegetation (w/m2)
  REAL(DP), INTENT(OUT)                :: SAG     !solar radiation absorbed by ground (w/m2)
  REAL(DP), INTENT(OUT)                :: FSA     !total absorbed solar radiation (w/m2)
  REAL(DP), INTENT(OUT)                :: FSR     !total reflected solar radiation (w/m2)
  REAL(DP), INTENT(OUT)                :: FSRV    !reflected solar radiation by vegetation
  REAL(DP), INTENT(OUT)                :: FSRG    !reflected solar radiation by ground

! ------------------------ local variables ----------------------------------------------------
  INTEGER                          :: IB      !waveband number (1=vis, 2=nir)
  INTEGER                          :: NBAND   !number of solar radiation waveband classes

  REAL(DP)                             :: ABS     !absorbed solar radiation (w/m2)
  REAL(DP)                             :: RNIR    !reflected solar radiation [nir] (w/m2)
  REAL(DP)                             :: RVIS    !reflected solar radiation [vis] (w/m2)
  REAL(DP)                             :: LAIFRA  !leaf area fraction of canopy
  REAL(DP)                             :: TRD     !transmitted solar radiation: direct (w/m2)
  REAL(DP)                             :: TRI     !transmitted solar radiation: diffuse (w/m2)
  REAL(DP), DIMENSION(1:2)             :: CAD     !direct beam absorbed by canopy (w/m2)
  REAL(DP), DIMENSION(1:2)             :: CAI     !diffuse radiation absorbed by canopy (w/m2)
! ---------------------------------------------------------------------------------------------
   NBAND = 2

! zero summed solar fluxes

    SAG = 0.
    SAV = 0.
    FSA = 0.

! loop over nband wavebands

  DO IB = 1, NBAND

! absorbed by canopy

    CAD(IB) = SOLAD(IB)*FABD(IB)
    CAI(IB) = SOLAI(IB)*FABI(IB)
    SAV     = SAV + CAD(IB) + CAI(IB)
    FSA     = FSA + CAD(IB) + CAI(IB)

! transmitted solar fluxes incident on ground

    TRD = SOLAD(IB)*FTDD(IB)
    TRI = SOLAD(IB)*FTID(IB) + SOLAI(IB)*FTII(IB)

! solar radiation absorbed by ground surface

    ABS = TRD*(1.-ALBGRD(IB)) + TRI*(1.-ALBGRI(IB))
    SAG = SAG + ABS
    FSA = FSA + ABS
  END DO

! partition visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves

     LAIFRA = ELAI / MAX(VAI,MPE)
     IF (FSUN .GT. 0.) THEN
        PARSUN = (CAD(1)+FSUN*CAI(1)) * LAIFRA / MAX(LAISUN,MPE)
        PARSHA = (FSHA*CAI(1))*LAIFRA / MAX(LAISHA,MPE)
     ELSE
        PARSUN = 0.
        PARSHA = (CAD(1)+CAI(1))*LAIFRA /MAX(LAISHA,MPE)
     ENDIF

! reflected solar radiation

     RVIS = ALBD(1)*SOLAD(1) + ALBI(1)*SOLAI(1)
     RNIR = ALBD(2)*SOLAD(2) + ALBI(2)*SOLAI(2)
     FSR  = RVIS + RNIR

! reflected solar radiation of veg. and ground (combined ground)
     FSRV = FREVD(1)*SOLAD(1)+FREVI(1)*SOLAI(1)+FREVD(2)*SOLAD(2)+FREVI(2)*SOLAI(2)
     FSRG = FREGD(1)*SOLAD(1)+FREGI(1)*SOLAI(1)+FREGD(2)*SOLAD(2)+FREGI(2)*SOLAI(2)

  END SUBROUTINE SURRAD
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  SUBROUTINE SNOW_AGE (DT,TG,SNEQVO,SNEQV,TAUSS,FAGE)
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! ------------------------ code history ------------------------------------------------------------
! from BATS
! ------------------------ input/output variables --------------------------------------------------
!input
   REAL(DP), INTENT(IN) :: DT        !main time step (s)
   REAL(DP), INTENT(IN) :: TG        !ground temperature (k)
   REAL(DP), INTENT(IN) :: SNEQVO    !snow mass at last time step(mm)
   REAL(DP), INTENT(IN) :: SNEQV     !snow water per unit ground area (mm)

!output
   REAL(DP), INTENT(OUT) :: FAGE     !snow age

!input/output
   REAL(DP), INTENT(INOUT) :: TAUSS      !non-dimensional snow age
!local
   REAL(DP)            :: TAGE       !total aging effects
   REAL(DP)            :: AGE1       !effects of grain growth due to vapor diffusion
   REAL(DP)            :: AGE2       !effects of grain growth at freezing of melt water
   REAL(DP)            :: AGE3       !effects of soot
   REAL(DP)            :: DELA       !temporary variable
   REAL(DP)            :: SGE        !temporary variable
   REAL(DP)            :: DELS       !temporary variable
   REAL(DP)            :: DELA0      !temporary variable
   REAL(DP)            :: ARG        !temporary variable
! See Yang et al. (1997) J.of Climate for detail.
!---------------------------------------------------------------------------------------------------

   IF(SNEQV.LE.0.0d0) THEN
          TAUSS = 0.
   ELSE IF (SNEQV.GT.800.) THEN
          TAUSS = 0.
   ELSE
          DELA0 = 1.E-6*DT
          ARG   = 5.E3*(1./TFRZ-1./TG)
          AGE1  = EXP(ARG)
          AGE2  = EXP(AMIN1(0.,10.*ARG))
          AGE3  = 0.3d0
          TAGE  = AGE1+AGE2+AGE3

          DELA  = DELA0*TAGE
          DELS  = AMAX1(0.0d0,SNEQV-SNEQVO) / SWEMX
          SGE   = (TAUSS+DELA)*(1.0d0-DELS)
          TAUSS = AMAX1(0.,SGE)
   ENDIF

   FAGE= TAUSS/(TAUSS+1.)

  END SUBROUTINE SNOW_AGE
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  SUBROUTINE SNOWALB_BATS (NBAND,FSNO,COSZ,FAGE,ALBSND,ALBSNI)
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  INTEGER,INTENT(IN) :: NBAND  !number of waveband classes

  REAL(DP),INTENT(IN) :: COSZ    !cosine solar zenith angle
  REAL(DP),INTENT(IN) :: FSNO    !snow cover fraction (-)
  REAL(DP),INTENT(IN) :: FAGE    !snow age correction

! output

  REAL(DP), DIMENSION(1:2),INTENT(OUT) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
  REAL(DP), DIMENSION(1:2),INTENT(OUT) :: ALBSNI !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  INTEGER :: IB          !waveband class

  REAL(DP) :: FZEN                 !zenith angle correction
  REAL(DP) :: CF1                  !temperary variable
  REAL(DP) :: SL2                  !2.*SL
  REAL(DP) :: SL1                  !1/SL
  REAL(DP) :: SL                   !adjustable parameter
  REAL(DP), PARAMETER :: C1 = 0.2d0  !default in BATS
  REAL(DP), PARAMETER :: C2 = 0.5d0  !default in BATS
!  REAL(DP), PARAMETER :: C1 = 0.2d0 * 2. ! double the default to match Sleepers River's
!  REAL(DP), PARAMETER :: C2 = 0.5d0 * 2. ! snow surface albedo (double aging effects)
! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        ALBSND(1: NBAND) = 0.
        ALBSNI(1: NBAND) = 0.

! when cosz > 0

        SL=2.0d0
        SL1=1./SL
        SL2=2.*SL
        CF1=((1.+SL1)/(1.+SL2*COSZ)-SL1)
        FZEN=AMAX1(CF1,0.)

        ALBSNI(1)=0.95d0*(1.-C1*FAGE)
        ALBSNI(2)=0.65d0*(1.-C2*FAGE)

        ALBSND(1)=ALBSNI(1)+0.4d0*FZEN*(1.-ALBSNI(1))    !  vis direct
        ALBSND(2)=ALBSNI(2)+0.4d0*FZEN*(1.-ALBSNI(2))    !  nir direct

  END SUBROUTINE SNOWALB_BATS
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  SUBROUTINE SNOWALB_CLASS (NBAND,QSNOW,DT,ALB,ALBOLD,ALBSND,ALBSNI,ILOC,JLOC)
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  INTEGER,INTENT(IN) :: ILOC !grid index
  INTEGER,INTENT(IN) :: JLOC !grid index
  INTEGER,INTENT(IN) :: NBAND  !number of waveband classes

  REAL(DP),INTENT(IN) :: QSNOW     !snowfall (mm/s)
  REAL(DP),INTENT(IN) :: DT        !time step (sec)
  REAL(DP),INTENT(IN) :: ALBOLD    !snow albedo at last time step

! in & out

  REAL(DP),                INTENT(INOUT) :: ALB        !
! output

  REAL(DP), DIMENSION(1:2),INTENT(OUT) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
  REAL(DP), DIMENSION(1:2),INTENT(OUT) :: ALBSNI !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  INTEGER :: IB          !waveband class

! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        ALBSND(1: NBAND) = 0.
        ALBSNI(1: NBAND) = 0.

! when cosz > 0

         ALB = 0.55d0 + (ALBOLD-0.55d0) * EXP(-0.01d0*DT/3600.)


! 1 mm fresh snow(SWE) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
! here assume 1cm snow depth will fully cover the old snow

         IF (QSNOW > 0.) then
           ALB = ALB + MIN(QSNOW*DT,SWEMX) * (0.84d0-ALB)/(SWEMX)
         ENDIF

         ALBSNI(1)= ALB         ! vis diffuse
         ALBSNI(2)= ALB         ! nir diffuse
         ALBSND(1)= ALB         ! vis direct
         ALBSND(2)= ALB         ! nir direct

  END SUBROUTINE SNOWALB_CLASS
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  SUBROUTINE GROUNDALB (NSOIL   ,NBAND   ,ICE     ,IST     ,ISC     , & !in
                        FSNO    ,SMC     ,ALBSND  ,ALBSNI  ,COSZ    , & !in
                        TG      ,ILOC    ,JLOC    ,                   & !in
                        ALBGRD  ,ALBGRI  )                              !out
! --------------------------------------------------------------------------------------------------
  USE NOAHMP_RAD_PARAMETERS
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
!input

  INTEGER,                  INTENT(IN)  :: ILOC   !grid index
  INTEGER,                  INTENT(IN)  :: JLOC   !grid index
  INTEGER,                  INTENT(IN)  :: NSOIL  !number of soil layers
  INTEGER,                  INTENT(IN)  :: NBAND  !number of solar radiation waveband classes
  INTEGER,                  INTENT(IN)  :: ICE    !value of ist for land ice
  INTEGER,                  INTENT(IN)  :: IST    !surface type
  INTEGER,                  INTENT(IN)  :: ISC    !soil color class (1-lighest; 8-darkest)
  REAL(DP),                     INTENT(IN)  :: FSNO   !fraction of surface covered with snow (-)
  REAL(DP),                     INTENT(IN)  :: TG     !ground temperature (k)
  REAL(DP),                     INTENT(IN)  :: COSZ   !cosine solar zenith angle (0-1)
  REAL(DP), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    !volumetric soil water content (m3/m3)
  REAL(DP), DIMENSION(1:    2), INTENT(IN)  :: ALBSND !direct beam snow albedo (vis, nir)
  REAL(DP), DIMENSION(1:    2), INTENT(IN)  :: ALBSNI !diffuse snow albedo (vis, nir)

!output

  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: ALBGRD !ground albedo (direct beam: vis, nir)
  REAL(DP), DIMENSION(1:    2), INTENT(OUT) :: ALBGRI !ground albedo (diffuse: vis, nir)

!local

  INTEGER                               :: IB     !waveband number (1=vis, 2=nir)
  REAL(DP)                                  :: INC    !soil water correction factor for soil albedo
  REAL(DP)                                  :: ALBSOD !soil albedo (direct)
  REAL(DP)                                  :: ALBSOI !soil albedo (diffuse)
! --------------------------------------------------------------------------------------------------

  DO IB = 1, NBAND
        INC = MAX(0.11d0-0.40d0*SMC(1), 0.)
        IF (IST .EQ. 1)  THEN                     !soil
           ALBSOD = MIN(ALBSAT(ISC,IB)+INC,ALBDRY(ISC,IB))
           ALBSOI = ALBSOD
        ELSE IF (TG .GT. TFRZ) THEN               !unfrozen lake, wetland
           ALBSOD = 0.06d0/(MAX(0.01d0,COSZ)**1.7d0 + 0.15d0)
           ALBSOI = 0.06d0
        ELSE                                      !frozen lake, wetland
           ALBSOD = ALBLAK(IB)
           ALBSOI = ALBSOD
        END IF

! increase desert and semi-desert albedos

        IF (IST .EQ. 1 .AND. ISC .EQ. 9) THEN
           ALBSOD = ALBSOD + 0.10d0
           ALBSOI = ALBSOI + 0.10d0
        end if

        ALBGRD(IB) = ALBSOD*(1.-FSNO) + ALBSND(IB)*FSNO
        ALBGRI(IB) = ALBSOI*(1.-FSNO) + ALBSNI(IB)*FSNO
  END DO

  END SUBROUTINE GROUNDALB
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  SUBROUTINE TWOSTREAM (IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
                        FWET   ,T       ,ALBGRD  ,ALBGRI  ,RHO    , & !in
                        TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
                        FAB    ,FRE     ,FTD     ,FTI     ,GDIR   , & !)   !out
                        FREV   ,FREG    ,BGAP    ,WGAP)

! --------------------------------------------------------------------------------------------------
! use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse
! flux given an underlying surface with known albedo.
! --------------------------------------------------------------------------------------------------
  USE NOAHMP_VEG_PARAMETERS
  USE NOAHMP_RAD_PARAMETERS
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

   INTEGER,              INTENT(IN)  :: ILOC    !grid index
   INTEGER,              INTENT(IN)  :: JLOC    !grid index
   INTEGER,              INTENT(IN)  :: IST     !surface type
   INTEGER,              INTENT(IN)  :: IB      !waveband number
   INTEGER,              INTENT(IN)  :: IC      !0=unit incoming direct; 1=unit incoming diffuse
   INTEGER,              INTENT(IN)  :: VEGTYP  !vegetation type

   REAL(DP),                 INTENT(IN)  :: COSZ    !cosine of direct zenith angle (0-1)
   REAL(DP),                 INTENT(IN)  :: VAI     !one-sided leaf+stem area index (m2/m2)
   REAL(DP),                 INTENT(IN)  :: FWET    !fraction of lai, sai that is wetted (-)
   REAL(DP),                 INTENT(IN)  :: T       !surface temperature (k)

   REAL(DP), DIMENSION(1:2), INTENT(IN)  :: ALBGRD  !direct  albedo of underlying surface (-)
   REAL(DP), DIMENSION(1:2), INTENT(IN)  :: ALBGRI  !diffuse albedo of underlying surface (-)
   REAL(DP), DIMENSION(1:2), INTENT(IN)  :: RHO     !leaf+stem reflectance
   REAL(DP), DIMENSION(1:2), INTENT(IN)  :: TAU     !leaf+stem transmittance
   REAL(DP),                 INTENT(IN)  :: FVEG    !green vegetation fraction [0.0d0-1.0d0]

! output

   REAL(DP), DIMENSION(1:2), INTENT(OUT) :: FAB     !flux abs by veg layer (per unit incoming flux)
   REAL(DP), DIMENSION(1:2), INTENT(OUT) :: FRE     !flux refl above veg layer (per unit incoming flux)
   REAL(DP), DIMENSION(1:2), INTENT(OUT) :: FTD     !down dir flux below veg layer (per unit in flux)
   REAL(DP), DIMENSION(1:2), INTENT(OUT) :: FTI     !down dif flux below veg layer (per unit in flux)
   REAL(DP),                 INTENT(OUT) :: GDIR    !projected leaf+stem area in solar direction
   REAL(DP), DIMENSION(1:2), INTENT(OUT) :: FREV    !flux reflected by veg layer   (per unit incoming flux)
   REAL(DP), DIMENSION(1:2), INTENT(OUT) :: FREG    !flux reflected by ground (per unit incoming flux)

! local
   REAL(DP)                              :: OMEGA   !fraction of intercepted radiation that is scattered
   REAL(DP)                              :: OMEGAL  !omega for leaves
   REAL(DP)                              :: BETAI   !upscatter parameter for diffuse radiation
   REAL(DP)                              :: BETAIL  !betai for leaves
   REAL(DP)                              :: BETAD   !upscatter parameter for direct beam radiation
   REAL(DP)                              :: BETADL  !betad for leaves
   REAL(DP)                              :: EXT     !optical depth of direct beam per unit leaf area
   REAL(DP)                              :: AVMU    !average diffuse optical depth

   REAL(DP)                              :: COSZI   !0.001d0 <= cosz <= 1.000d0
   REAL(DP)                              :: ASU     !single scattering albedo
   REAL(DP)                              :: CHIL    ! -0.4d0 <= xl <= 0.6d0

   REAL(DP)                              :: TMP0,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP7,TMP8,TMP9
   REAL(DP)                              :: P1,P2,P3,P4,S1,S2,U1,U2,U3
   REAL(DP)                              :: B,C,D,D1,D2,F,H,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10
   REAL(DP)                              :: PHI1,PHI2,SIGMA
   REAL(DP)                              :: FTDS,FTIS,FRES
   REAL(DP)                              :: DENFVEG
   REAL(DP)                              :: VAI_SPREAD
!jref:start
   REAL(DP)                              :: FREVEG,FREBAR,FTDVEG,FTIVEG,FTDBAR,FTIBAR
   REAL(DP)                              :: THETAZ
!jref:end

!  variables for the modified two-stream scheme
!  Niu and Yang (2004), JGR

   REAL(DP), PARAMETER :: PAI = 3.14159265d0
   REAL(DP) :: HD       !crown depth (m)
   REAL(DP) :: BB       !vertical crown radius (m)
   REAL(DP) :: THETAP   !angle conversion from SZA
   REAL(DP) :: FA       !foliage volume density (m-1)
   REAL(DP) :: NEWVAI   !effective LSAI (-)

   REAL(DP),INTENT(INOUT) :: BGAP     !between canopy gap fraction for beam (-)
   REAL(DP),INTENT(INOUT) :: WGAP     !within canopy gap fraction for beam (-)

   REAL(DP) :: KOPEN    !gap fraction for diffue light (-)
   REAL(DP) :: GAP      !total gap fraction for beam ( <=1-shafac )

! -----------------------------------------------------------------
! compute within and between gaps
     VAI_SPREAD = VAI
     if(VAI == 0.0d0) THEN
         GAP     = 1.0d0
         KOPEN   = 1.0d0
     ELSE
         IF(OPT_RAD == 1) THEN
           DENFVEG = -LOG(MAX(1.0d0-FVEG,0.01d0))/(PAI*RC(VEGTYP)**2)
           HD      = HVT(VEGTYP) - HVB(VEGTYP)
           BB      = 0.5d0 * HD
           THETAP  = ATAN(BB/RC(VEGTYP) * TAN(ACOS(MAX(0.01d0,COSZ))) )
           ! BGAP    = EXP(-DEN(VEGTYP) * PAI * RC(VEGTYP)**2/COS(THETAP) )
           BGAP    = EXP(-DENFVEG * PAI * RC(VEGTYP)**2/COS(THETAP) )
           FA      = VAI/(1.33d0 * PAI * RC(VEGTYP)**3.0d0 *(BB/RC(VEGTYP))*DENFVEG)
           NEWVAI  = HD*FA
           WGAP    = (1.0d0-BGAP) * EXP(-0.5d0*NEWVAI/COSZ)
           GAP     = MIN(1.0d0-FVEG, BGAP+WGAP)

           KOPEN   = 0.05d0
         END IF

         IF(OPT_RAD == 2) THEN
           GAP     = 0.0d0
           KOPEN   = 0.0d0
         END IF

         IF(OPT_RAD == 3) THEN
           GAP     = 1.0d0-FVEG
           KOPEN   = 1.0d0-FVEG
         END IF
     end if

! calculate two-stream parameters OMEGA, BETAD, BETAI, AVMU, GDIR, EXT.
! OMEGA, BETAD, BETAI are adjusted for snow. values for OMEGA*BETAD
! and OMEGA*BETAI are calculated and then divided by the new OMEGA
! because the product OMEGA*BETAI, OMEGA*BETAD is used in solution.
! also, the transmittances and reflectances (TAU, RHO) are linear
! weights of leaf and stem values.

     COSZI  = MAX(0.001d0, COSZ)
     CHIL   = MIN( MAX(XL(VEGTYP), -0.4d0), 0.6d0)
     IF (ABS(CHIL) .LE. 0.01d0) CHIL = 0.01d0
     PHI1   = 0.5d0 - 0.633d0*CHIL - 0.330d0*CHIL*CHIL
     PHI2   = 0.877d0 * (1.-2.*PHI1)
     GDIR   = PHI1 + PHI2*COSZI
     EXT    = GDIR/COSZI
     AVMU   = ( 1. - PHI1/PHI2 * LOG((PHI1+PHI2)/PHI1) ) / PHI2
     OMEGAL = RHO(IB) + TAU(IB)
     TMP0   = GDIR + PHI2*COSZI
     TMP1   = PHI1*COSZI
     ASU    = 0.5d0*OMEGAL*GDIR/TMP0 * ( 1.-TMP1/TMP0*LOG((TMP1+TMP0)/TMP1) )
     BETADL = (1.+AVMU*EXT)/(OMEGAL*AVMU*EXT)*ASU
     BETAIL = 0.5d0 * ( RHO(IB)+TAU(IB) + (RHO(IB)-TAU(IB))   &
            * ((1.+CHIL)/2.)**2 ) / OMEGAL

! adjust omega, betad, and betai for intercepted snow

     IF (T .GT. TFRZ) THEN                                !no snow
        TMP0 = OMEGAL
        TMP1 = BETADL
        TMP2 = BETAIL
     ELSE
        TMP0 =   (1.-FWET)*OMEGAL        + FWET*OMEGAS(IB)
        TMP1 = ( (1.-FWET)*OMEGAL*BETADL + FWET*OMEGAS(IB)*BETADS ) / TMP0
        TMP2 = ( (1.-FWET)*OMEGAL*BETAIL + FWET*OMEGAS(IB)*BETAIS ) / TMP0
     END IF

     OMEGA = TMP0
     BETAD = TMP1
     BETAI = TMP2

! absorbed, reflected, transmitted fluxes per unit incoming radiation

     B = 1. - OMEGA + OMEGA*BETAI
     C = OMEGA*BETAI
     TMP0 = AVMU*EXT
     D = TMP0 * OMEGA*BETAD
     F = TMP0 * OMEGA*(1.-BETAD)
     TMP1 = B*B - C*C
     H = SQRT(TMP1) / AVMU
     SIGMA = TMP0*TMP0 - TMP1
     if ( ABS (SIGMA) < 1.e-6_dp) SIGMA = SIGN(1.e-6_dp,SIGMA)
     P1 = B + AVMU*H
     P2 = B - AVMU*H
     P3 = B + TMP0
     P4 = B - TMP0
     S1 = EXP(-H*VAI)
     ! MPC change: precision issues
     S2 = max(epsilon(VAI),EXP(-EXT*VAI))
     IF (IC .EQ. 0) THEN
        U1 = B - C/ALBGRD(IB)
        U2 = B - C*ALBGRD(IB)
        U3 = F + C*ALBGRD(IB)
     ELSE
        U1 = B - C/ALBGRI(IB)
        U2 = B - C*ALBGRI(IB)
        U3 = F + C*ALBGRI(IB)
     END IF
     TMP2 = U1 - AVMU*H
     TMP3 = U1 + AVMU*H
     D1 = P1*TMP2/S1 - P2*TMP3*S1
     TMP4 = U2 + AVMU*H
     TMP5 = U2 - AVMU*H
     D2 = TMP4/S1 - TMP5*S1
     H1 = -D*P4 - C*F
     TMP6 = D - H1*P3/SIGMA
     TMP7 = ( D - C - H1/SIGMA*(U1+TMP0) ) * S2
     H2 = ( TMP6*TMP2/S1 - P2*TMP7 ) / D1
     H3 = - ( TMP6*TMP3*S1 - P1*TMP7 ) / D1
     H4 = -F*P3 - C*D
     TMP8 = H4/SIGMA
     TMP9 = ( U3 - TMP8*(U2-TMP0) ) * S2
     H5 = - ( TMP8*TMP4/S1 + TMP9 ) / D2
     H6 = ( TMP8*TMP5*S1 + TMP9 ) / D2
     H7 = (C*TMP2) / (D1*S1)
     H8 = (-C*TMP3*S1) / D1
     H9 = TMP4 / (D2*S1)
     H10 = (-TMP5*S1) / D2

! downward direct and diffuse fluxes below vegetation
! Niu and Yang (2004), JGR.

     IF (IC .EQ. 0) THEN
        FTDS = S2                           *(1.0d0-GAP) + GAP
        FTIS = (H4*S2/SIGMA + H5*S1 + H6/S1)*(1.0d0-GAP)
     ELSE
        FTDS = 0.
        FTIS = (H9*S1 + H10/S1)*(1.0d0-KOPEN) + KOPEN
     END IF
     FTD(IB) = FTDS
     FTI(IB) = FTIS

! flux reflected by the surface (veg. and ground)

     IF (IC .EQ. 0) THEN
        FRES   = (H1/SIGMA + H2 + H3)*(1.0d0-GAP  ) + ALBGRD(IB)*GAP
        FREVEG = (H1/SIGMA + H2 + H3)*(1.0d0-GAP  )
        FREBAR = ALBGRD(IB)*GAP                   !jref - separate veg. and ground reflection
     ELSE
        FRES   = (H7 + H8) *(1.0d0-KOPEN) + ALBGRI(IB)*KOPEN
        FREVEG = (H7 + H8) *(1.0d0-KOPEN) + ALBGRI(IB)*KOPEN
        FREBAR = 0                                !jref - separate veg. and ground reflection
     END IF
     FRE(IB) = FRES

     FREV(IB) = FREVEG
     FREG(IB) = FREBAR

! flux absorbed by vegetation

     FAB(IB) = 1. - FRE(IB) - (1.-ALBGRD(IB))*FTD(IB) &
                            - (1.-ALBGRI(IB))*FTI(IB)

!if(iloc == 1.and.jloc ==  2) then
  !write(*,'(a7,2i2,5(a6,f8.4d0),2(a9,f8.4d0))') "ib,ic: ",ib,ic," GAP: ",GAP," FTD: ",FTD(IB)," FTI: ",FTI(IB)," FRE: ", &
  !       FRE(IB)," FAB: ",FAB(IB)," ALBGRD: ",ALBGRD(IB)," ALBGRI: ",ALBGRI(IB)
!end if

  END SUBROUTINE TWOSTREAM
! ==================================================================================================
  ! ==================================================================================================
! ----------------------------------------------------------------------
  SUBROUTINE STOMATA (VEGTYP  ,MPE     ,APAR    ,FOLN    ,ILOC    , JLOC, & !in
                      TV      ,EI      ,EA      ,SFCTMP  ,SFCPRS  , & !in
                      O2      ,CO2     ,IGS     ,BTRAN   ,RB      , & !in
                      RS      ,PSN     )                              !out
! --------------------------------------------------------------------------------------------------
  USE NOAHMP_VEG_PARAMETERS
  USE MODULE_MODEL_CONSTANTS
  use nrtype
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
      INTEGER,INTENT(IN)  :: ILOC   !grid index
      INTEGER,INTENT(IN)  :: JLOC   !grid index
      INTEGER,INTENT(IN)  :: VEGTYP !vegetation physiology type

      REAL(DP), INTENT(IN)    :: IGS    !growing season index (0=off, 1=on)
      REAL(DP), INTENT(IN)    :: MPE    !prevents division by zero errors
      REAL(DP), INTENT(IN)    :: TV     !foliage temperature (k)
      REAL(DP), INTENT(IN)    :: EI     !vapor pressure inside leaf (sat vapor press at tv) (pa)
      REAL(DP), INTENT(IN)    :: EA     !vapor pressure of canopy air (pa)
      REAL(DP), INTENT(IN)    :: APAR   !par absorbed per unit lai (w/m2)
      REAL(DP), INTENT(IN)    :: O2     !atmospheric o2 concentration (pa)
      REAL(DP), INTENT(IN)    :: CO2    !atmospheric co2 concentration (pa)
      REAL(DP), INTENT(IN)    :: SFCPRS !air pressure at reference height (pa)
      REAL(DP), INTENT(IN)    :: SFCTMP !air temperature at reference height (k)
      REAL(DP), INTENT(IN)    :: BTRAN  !soil water transpiration factor (0 to 1)
      REAL(DP), INTENT(IN)    :: FOLN   !foliage nitrogen concentration (%)
      REAL(DP), INTENT(IN)    :: RB     !boundary layer resistance (s/m)

! output
      REAL(DP), INTENT(OUT)   :: RS     !leaf stomatal resistance (s/m)
      REAL(DP), INTENT(OUT)   :: PSN    !foliage photosynthesis (umol co2 /m2/ s) [always +]

! in&out
      REAL(DP)               :: RLB    !boundary layer resistance (s m2 / umol)
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
      INTEGER :: ITER     !iteration index
      INTEGER :: NITER    !number of iterations

      DATA NITER /3/
      SAVE NITER

      REAL(DP) :: AB          !used in statement functions
      REAL(DP) :: BC          !used in statement functions
      REAL(DP) :: F1          !generic temperature response (statement function)
      REAL(DP) :: F2          !generic temperature inhibition (statement function)
      REAL(DP) :: TC          !foliage temperature (degree Celsius)
      REAL(DP) :: CS          !co2 concentration at leaf surface (pa)
      REAL(DP) :: KC          !co2 Michaelis-Menten constant (pa)
      REAL(DP) :: KO          !o2 Michaelis-Menten constant (pa)
      REAL(DP) :: A,B,C,Q     !intermediate calculations for RS
      REAL(DP) :: R1,R2       !roots for RS
      REAL(DP) :: FNF         !foliage nitrogen adjustment factor (0 to 1)
      REAL(DP) :: PPF         !absorb photosynthetic photon flux (umol photons/m2/s)
      REAL(DP) :: WC          !Rubisco limited photosynthesis (umol co2/m2/s)
      REAL(DP) :: WJ          !light limited photosynthesis (umol co2/m2/s)
      REAL(DP) :: WE          !export limited photosynthesis (umol co2/m2/s)
      REAL(DP) :: CCP          !co2 compensation point (pa)
      REAL(DP) :: CI          !internal co2 (pa)
      REAL(DP) :: AWC         !intermediate calculation for wc
      REAL(DP) :: VCMX        !maximum rate of carbonylation (umol co2/m2/s)
      REAL(DP) :: J           !electron transport (umol co2/m2/s)
      REAL(DP) :: CEA         !constrain ea or else model blows up
      REAL(DP) :: CF          !s m2/umol -> s/m

      F1(AB,BC) = AB**((BC-25.)/10.)
      F2(AB) = 1. + EXP((-2.2E05+710.*(AB+273.16))/(8.314*(AB+273.16)))
      REAL(DP) :: T
! ---------------------------------------------------------------------------------------------

! MPC change
! Original code uses gs = g0 + g1*hs*An/cs (g0 is not limited by soil moisture stress)
! Following Bonan et al. (JGR 2011) Improving canopy processes in the CLM... we replace the original code
!  with gs = g0*bTran + g1*hs*An/cs, where bTran is the soil moisture stress function

! initialize RS=RSMAX and PSN=0 because will only do calculations
! for APAR > 0, in which case RS <= RSMAX and PSN >= 0

         CF = SFCPRS/(8.314*SFCTMP)*1.e06
         RS = 1./(BP(VEGTYP)*BTRAN) * CF    ! MPC change: include BTRAN multiplier
         PSN = 0.
         !write(*,'(a,1x,20(f16.6,1x))') 'TV-TFRZ, RS, CF = ', TV-TFRZ, RS, CF

         IF (APAR .LE. 0.) RETURN

         FNF = MIN( FOLN/MAX(MPE,FOLNMX(VEGTYP)), 1.0 )
         TC  = TV-TFRZ
         PPF = 4.6*APAR
         J   = PPF*QE25(VEGTYP)
         !write(*,'(a,1x,10(f20.10,1x))') 'J, APAR, QE25(VEGTYP) = ', J, APAR, QE25(VEGTYP)
         KC  = KC25(VEGTYP) * F1(AKC(VEGTYP),TC)
         KO  = KO25(VEGTYP) * F1(AKO(VEGTYP),TC)
         AWC = KC * (1.+O2/KO)
         CCP  = 0.5_dp*KC/KO*O2*0.21_dp
         VCMX = VCMX25(VEGTYP) / F2(TC) * FNF * BTRAN * F1(AVCMX(VEGTYP),TC)
         !write(*,'(a,1x,20(f14.8,1x))') 'F1, F2, VCMX25(VEGTYP), AVCMX(VEGTYP), TC = ', F1(AVCMX(VEGTYP),TC), F2(TC), VCMX25(VEGTYP), AVCMX(VEGTYP), TC

! first guess ci
         CI = 0.7*CO2*C3PSN(VEGTYP) + 0.4*CO2*(1.-C3PSN(VEGTYP))
! rb: s/m -> s m**2 / umol
         RLB = RB/CF
! constrain ea
         CEA = MAX(0.25*EI*C3PSN(VEGTYP)+0.40*EI*(1.-C3PSN(VEGTYP)), MIN(EA,EI) )


! ci iteration
!jref: C3PSN is equal to 1 for all veg types.
       DO ITER = 1, NITER

            WJ = MAX(CI-CCP,0.)*J/(CI+2.*CCP)*C3PSN(VEGTYP)  + J*(1.-C3PSN(VEGTYP))
            WC = MAX(CI-CCP,0.)*VCMX/(CI+AWC)*C3PSN(VEGTYP) + VCMX*(1.-C3PSN(VEGTYP))
            WE = 0.5*VCMX*C3PSN(VEGTYP) + 4000.*VCMX*CI/SFCPRS*(1.-C3PSN(VEGTYP))
            PSN = MIN(WJ,WC,WE) * IGS

            CS = MAX( CO2-1.37*RLB*SFCPRS*PSN, MPE )

            A = MP(VEGTYP)*PSN*SFCPRS*CEA / (CS*EI) + BTRAN*BP(VEGTYP)      ! MPC change: include BTRAN multiplier for 2nd term
            B = ( MP(VEGTYP)*PSN*SFCPRS/CS + BTRAN*BP(VEGTYP) ) * RLB - 1.  ! MPC change: include BTRAN multiplier for 2nd term in brackets
            C = -RLB
            IF (B .GE. 0.) THEN
               Q = -0.5*( B + SQRT(B*B-4.*A*C) )
            ELSE
               Q = -0.5*( B - SQRT(B*B-4.*A*C) )
            END IF
            R1 = Q/A
            R2 = C/Q
            RS = MAX(R1,R2)
            CI = MAX( CS-PSN*SFCPRS*1.65*RS, 0. )
       END DO

! rs, rb:  s m**2 / umol -> s/m
         RS = RS*CF

  END SUBROUTINE STOMATA
! ==================================================================================================
  SUBROUTINE CANRES (PAR   ,SFCTMP,RCSOIL ,EAH   ,SFCPRS , & !in
                     RC    ,PSN   ,ILOC   ,JLOC  )           !out

! --------------------------------------------------------------------------------------------------
! calculate canopy resistance which depends on incoming solar radiation,
! air temperature, atmospheric water vapor pressure deficit at the
! lowest model level, and soil moisture (preferably unfrozen soil
! moisture rather than total)
! --------------------------------------------------------------------------------------------------
! source:  Jarvis (1976), Noilhan and Planton (1989, MWR), Jacquemin and
! Noilhan (1990, BLM). Chen et al (1996, JGR, Vol 101(D3), 7251-7268),
! eqns 12-14 and table 2 of sec. 3.1.2
! --------------------------------------------------------------------------------------------------
!niu    USE module_Noahlsm_utility
! --------------------------------------------------------------------------------------------------
    IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs

    INTEGER,                  INTENT(IN)  :: ILOC   !grid index
    INTEGER,                  INTENT(IN)  :: JLOC   !grid index
    REAL(DP),                     INTENT(IN)  :: PAR    !par absorbed per unit sunlit lai (w/m2)
    REAL(DP),                     INTENT(IN)  :: SFCTMP !canopy air temperature
    REAL(DP),                     INTENT(IN)  :: SFCPRS !surface pressure (pa)
    REAL(DP),                     INTENT(IN)  :: EAH    !water vapor pressure (pa)
    REAL(DP),                     INTENT(IN)  :: RCSOIL !soil moisture stress factor

!outputs

    REAL(DP),                     INTENT(OUT) :: RC     !canopy resistance per unit LAI
    REAL(DP),                     INTENT(OUT) :: PSN    !foliage photosynthesis (umolco2/m2/s)

!local

    REAL(DP)                                  :: RCQ
    REAL(DP)                                  :: RCS
    REAL(DP)                                  :: RCT
    REAL(DP)                                  :: FF
    REAL(DP)                                  :: Q2     !water vapor mixing ratio (kg/kg)
    REAL(DP)                                  :: Q2SAT  !saturation Q2
    REAL(DP)                                  :: DQSDT2 !d(Q2SAT)/d(T)

! RSMIN, RSMAX, TOPT, RGL, HS are canopy stress parameters set in REDPRM
! ----------------------------------------------------------------------
! initialize canopy resistance multiplier terms.
! ----------------------------------------------------------------------
    RC     = 0.0_dp
    RCS    = 0.0_dp
    RCT    = 0.0_dp
    RCQ    = 0.0_dp

!  compute Q2 and Q2SAT

    Q2 = 0.622_dp *  EAH  / (SFCPRS - 0.378_dp * EAH) !specific humidity [kg/kg]
    Q2 = Q2 / (1.0_dp + Q2)                        !mixing ratio [kg/kg]
    CALL CALHUM(SFCTMP, SFCPRS, Q2SAT, DQSDT2)

! contribution due to incoming solar radiation

    FF  = 2.0 * PAR / RGL
    RCS = (FF + RSMIN / RSMAX) / (1.0+ FF)
    RCS = MAX (RCS,0.0001)

! contribution due to air temperature

    RCT = 1.0- 0.0016* ( (TOPT - SFCTMP)**2.0)
    RCT = MAX (RCT,0.0001)

! contribution due to vapor pressure deficit

    RCQ = 1.0/ (1.0+ HS * MAX(0.,Q2SAT-Q2))
    RCQ = MAX (RCQ,0.01)

! determine canopy resistance due to all factors

    RC  = RSMIN / (RCS * RCT * RCQ * RCSOIL)
    PSN = -999.99       ! PSN not applied for dynamic carbon

  END SUBROUTINE CANRES
! ==================================================================================================
  SUBROUTINE CALHUM(SFCTMP, SFCPRS, Q2SAT, DQSDT2)
        IMPLICIT NONE
        REAL(DP), INTENT(IN)       :: SFCTMP, SFCPRS
        REAL(DP), INTENT(OUT)      :: Q2SAT, DQSDT2
        REAL(DP), PARAMETER        :: A2=17.67_dp,A3=273.15_dp,A4=29.65_dp, ELWV=2.501D6,         &
                                  A23M4=A2*(A3-A4), E0=0.611_dp, RV=461.0_dp,             &
                                  EPSILON=0.622_dp
        REAL(DP)                   :: ES, SFCPRSX
! Q2SAT: saturated mixing ratio
        ES = E0 * EXP ( ELWV/RV*(1./A3 - 1./SFCTMP) )
! convert SFCPRS from Pa to KPa
        SFCPRSX = SFCPRS*1.E-3
        Q2SAT = EPSILON * ES / (SFCPRSX-ES)
! convert from  g/g to g/kg
        Q2SAT = Q2SAT * 1.E3
! Q2SAT is currently a 'mixing ratio'
! DQSDT2 is calculated assuming Q2SAT is a specific humidity
        DQSDT2=(Q2SAT/(1+Q2SAT))*A23M4/(SFCTMP-A4)**2
! DG Q2SAT needs to be in g/g when returned for SFLX
        Q2SAT = Q2SAT / 1.E3
  END SUBROUTINE CALHUM
! ==================================================================================================
! ********************************* end of carbon subroutines *****************************
! ==================================================================================================
  SUBROUTINE REDPRM (VEGTYP,SOILTYP,SLOPETYP,ZSOIL,NSOIL,ISURBAN)

!niu  use module_sf_noahlsm_param_init

    IMPLICIT NONE
! ----------------------------------------------------------------------
! Internally set (default valuess)
! all soil and vegetation parameters required for the execusion oF
! the Noah lsm are defined in VEGPARM.TBL, SOILPARM.TB, and GENPARM.TBL.
! ----------------------------------------------------------------------
!     Vegetation parameters:
!             CMXTBL: MAX CNPY Capacity
!              NROOT: Rooting depth
!
! ----------------------------------------------------------------------
!      Soil parameters:
!       SSATPSI: SAT (saturation) soil potential
!        SSATDW: SAT soil diffusivity
!           F1: Soil thermal diffusivity/conductivity coef.
!        QUARTZ: Soil quartz content
!  Modified by F. Chen (12/22/97)  to use the STATSGO soil map
!  Modified By F. Chen (01/22/00)  to include PLaya, Lava, and White San
!  Modified By F. Chen (08/05/02)  to include additional parameters for the Noah
! NOTE: SATDW = BB*SATDK*(SATPSI/MAXSMC)
!         F11 = log10(SATPSI) + BB*log10(MAXSMC) + 2.0d0
!       REFSMC1=MAXSMC*(5.79E-9/SATDK)**(1/(2*BB+3)) 5.79E-9 m/s= 0.5d0 mm
!       REFSMC=REFSMC1+1./3.(MAXSMC-REFSMC1)
!       WLTSMC1=MAXSMC*(200./SATPSI)**(-1./BB)    (Wetzel and Chang, 198
!       WLTSMC=WLTSMC1-0.5d0*WLTSMC1
! Note: the values for playa is set for it to have a thermal conductivit
! as sand and to have a hydrulic conductivity as clay
!
! ----------------------------------------------------------------------
! BLANK        OCEAN/SEA
!       CSOIL_DATA: soil heat capacity [J M-3 K-1]
!       ZBOT_DATA: depth[M] of lower boundary soil temperature
!       CZIL_DATA: calculate roughness length of heat
!       SMLOW_DATA and MHIGH_DATA: two soil moisture wilt, soil moisture referen
!                 parameters
! Set maximum number of soil- and veg- in data statement.
! ----------------------------------------------------------------------
    INTEGER, PARAMETER     :: MAX_SOILTYP=30,MAX_VEGTYP=30

! Veg parameters
    INTEGER, INTENT(IN)    :: VEGTYP
    INTEGER, INTENT(IN)    :: ISURBAN
! Soil parameters
    INTEGER, INTENT(IN)    :: SOILTYP
! General parameters
    INTEGER, INTENT(IN)    :: SLOPETYP
! General parameters
    INTEGER, INTENT(IN)    :: NSOIL
! Layer parameters
    REAL(DP),DIMENSION(NSOIL),INTENT(IN) :: ZSOIL

! Locals
    REAL(DP)    :: REFDK
    REAL(DP)    :: REFKDT
    REAL(DP)    :: FRZK
    REAL(DP)                   :: FRZFACT
    INTEGER                :: I
    CHARACTER(len=256)     :: message
! ----------------------------------------------------------------------
!
    IF (SOILTYP .gt. SLCATS) THEN
       call wrf_message('SOILTYP must be less than SLCATS:')
       write(message, '("SOILTYP = ", I6, ";    SLCATS = ", I6)') SOILTYP, SLCATS
       call wrf_message(trim(message))
       call wrf_error_fatal ('REDPRM: Error: too many input soil types')
    END IF
    IF (VEGTYP .gt. LUCATS) THEN
       call wrf_message('VEGTYP must be less than LUCATS:')
       write(message, '("VEGTYP = ", I6, ";    LUCATS = ", I6)') VEGTYP, LUCATS
       call wrf_message(trim(message))
       call wrf_error_fatal ('Error: too many input landuse types')
    END IF

! ----------------------------------------------------------------------
!  SET-UP SOIL PARAMETERS
! ----------------------------------------------------------------------
    CSOIL  = CSOIL_DATA
    BEXP   = BB (SOILTYP)
    DKSAT  = SATDK (SOILTYP)
    DWSAT  = SATDW (SOILTYP)
    F1     = F11 (SOILTYP)
    PSISAT = SATPSI (SOILTYP)
    QUARTZ = QTZ (SOILTYP)
    SMCDRY = DRYSMC (SOILTYP)
    SMCMAX = MAXSMC (SOILTYP)
    SMCREF = REFSMC (SOILTYP)
    SMCWLT = WLTSMC (SOILTYP)

    IF(VEGTYP==ISURBAN)THEN
       SMCMAX = 0.45d0
       SMCREF = 0.42d0
       SMCWLT = 0.40d0
       SMCDRY = 0.40d0
       CSOIL  = 3.E6
    ENDIF

! ----------------------------------------------------------------------
! Set-up universal parameters (not dependent on SOILTYP, VEGTYP)
! ----------------------------------------------------------------------
    ZBOT   = ZBOT_DATA
    CZIL   = CZIL_DATA

    FRZK   = FRZK_DATA
    REFDK  = REFDK_DATA
    REFKDT = REFKDT_DATA
    KDT    = REFKDT * DKSAT / REFDK
    SLOPE  = SLOPE_DATA (SLOPETYP)

! adjust FRZK parameter to actual soil type: FRZK * FRZFACT

    if(SOILTYP /= 14) then
      FRZFACT = (SMCMAX / SMCREF) * (0.412d0 / 0.468d0)
      FRZX = FRZK * FRZFACT
    end if

!    write(*,*) FRZK, FRZX, KDT, SLOPE, SLOPETYP
! ----------------------------------------------------------------------
! SET-UP VEGETATION PARAMETERS
! ----------------------------------------------------------------------
    ! Six redprm_canres variables:
    TOPT = TOPT_DATA
    RGL = RGLTBL (VEGTYP)
    RSMAX = RSMAX_DATA
    RSMIN = RSTBL (VEGTYP)
    HS = HSTBL (VEGTYP)
    NROOT = NROTBL (VEGTYP)

    IF(VEGTYP==ISURBAN)THEN
       RSMIN=400.0d0
    ENDIF

!    SHDFAC = SHDTBL(VEGTYP)
!    IF (VEGTYP .eq. BARE) SHDFAC = 0.0d0

    IF (NROOT .gt. NSOIL) THEN
       WRITE (*,*) 'Warning: too many root layers'
       write (*,*) 'NROOT = ', nroot
       write (*,*) 'NSOIL = ', nsoil
       call wrf_error_fatal("STOP in Noah-MP")
    END IF

! ----------------------------------------------------------------------

  END  SUBROUTINE REDPRM
! ==================================================================================================

  subroutine noahmp_options(idveg     ,iopt_crs  ,iopt_btr  ,iopt_run  ,iopt_sfc  ,iopt_frz , &
                             iopt_inf  ,iopt_rad  ,iopt_alb  ,iopt_snf  ,iopt_tbot, iopt_stc )

  implicit none

  INTEGER,  INTENT(IN) :: idveg     !dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
  INTEGER,  INTENT(IN) :: iopt_crs  !canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  INTEGER,  INTENT(IN) :: iopt_btr  !soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  INTEGER,  INTENT(IN) :: iopt_run  !runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  INTEGER,  INTENT(IN) :: iopt_sfc  !surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  INTEGER,  INTENT(IN) :: iopt_frz  !supercooled liquid water (1-> NY06; 2->Koren99)
  INTEGER,  INTENT(IN) :: iopt_inf  !frozen soil permeability (1-> NY06; 2->Koren99)
  INTEGER,  INTENT(IN) :: iopt_rad  !radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
  INTEGER,  INTENT(IN) :: iopt_alb  !snow surface albedo (1->BATS; 2->CLASS)
  INTEGER,  INTENT(IN) :: iopt_snf  !rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
  INTEGER,  INTENT(IN) :: iopt_tbot !lower boundary of soil temperature (1->zero-flux; 2->Noah)

  INTEGER,  INTENT(IN) :: iopt_stc  !snow/soil temperature time scheme (only layer 1)
                                    ! 1 -> semi-implicit; 2 -> full implicit (original Noah)

! -------------------------------------------------------------------------------------------------

  dveg = idveg

  opt_crs  = iopt_crs
  opt_btr  = iopt_btr
  opt_run  = iopt_run
  opt_sfc  = iopt_sfc
  opt_frz  = iopt_frz
  opt_inf  = iopt_inf
  opt_rad  = iopt_rad
  opt_alb  = iopt_alb
  opt_snf  = iopt_snf
  opt_tbot = iopt_tbot
  opt_stc  = iopt_stc

  end subroutine noahmp_options

END MODULE NOAHMP_ROUTINES
! ==================================================================================================

MODULE MODULE_SF_NOAHMPLSM

  USE NOAHMP_ROUTINES
  USE NOAHMP_GLOBALS
  USE NOAHMP_VEG_PARAMETERS

END MODULE MODULE_SF_NOAHMPLSM
