module noahmp_globals
  USE nrtype

  ! Maybe most of these can be moved to a REDPRM use statement?
  ! MPC -- yes, all of these variables can be local to REDPRM (see additional comments)
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

  REAL(rkind), PARAMETER :: GRAV   = 9.80616   !acceleration due to gravity (m/s2)
  REAL(rkind), PARAMETER :: SB     = 5.67E-08  !Stefan-Boltzmann constant (w/m2/k4)
  REAL(rkind), PARAMETER :: VKC    = 0.40      !von Karman constant
  REAL(rkind), PARAMETER :: TFRZ   = 273.16    !freezing/melting point (k)
  REAL(rkind), PARAMETER :: HSUB   = 2.8440E06 !latent heat of sublimation (j/kg)
  REAL(rkind), PARAMETER :: HVAP   = 2.5104E06 !latent heat of vaporization (j/kg)
  REAL(rkind), PARAMETER :: HFUS   = 0.3336E06 !latent heat of fusion (j/kg)
  REAL(rkind), PARAMETER :: CWAT   = 4.188E06  !specific heat capacity of water (j/m3/k)
  REAL(rkind), PARAMETER :: CICE   = 2.094E06  !specific heat capacity of ice (j/m3/k)
  REAL(rkind), PARAMETER :: CPAIR  = 1004.64   !heat capacity dry air at const pres (j/kg/k)
  REAL(rkind), PARAMETER :: TKWAT  = 0.6       !thermal conductivity of water (w/m/k)
  REAL(rkind), PARAMETER :: TKICE  = 2.2       !thermal conductivity of ice (w/m/k)
  REAL(rkind), PARAMETER :: TKAIR  = 0.023     !thermal conductivity of air (w/m/k)
  REAL(rkind), PARAMETER :: RAIR   = 287.04    !gas constant for dry air (j/kg/k)
  REAL(rkind), PARAMETER :: RW     = 461.269   !gas constant for  water vapor (j/kg/k)
  REAL(rkind), PARAMETER :: DENH2O = 1000.     !density of water (kg/m3)
  REAL(rkind), PARAMETER :: DENICE = 917.      !density of ice (kg/m3)

!------------------------------------------------------------------------------------------!
! From the VEGPARM.TBL tables, as functions of vegetation category.
!------------------------------------------------------------------------------------------!
  INTEGER :: NROOT        !rooting depth [as the number of layers] ( Assigned in REDPRM )
  REAL(rkind)    :: RGL          !parameter used in radiation stress function ( Assigned in REDPRM )
  REAL(rkind)    :: RSMIN        !minimum Canopy Resistance [s/m] ( Assigned in REDPRM )
  REAL(rkind)    :: HS           !parameter used in vapor pressure deficit function ( Assigned in REDPRM )
  REAL(rkind)    :: RSMAX        !maximum stomatal resistance ( Assigned in REDPRM )
  REAL(rkind)    :: TOPT         !optimum transpiration air temperature.

  ! MPC change: make variables private for a given thread
  !$omp threadprivate(NROOT, RGL, RSMIN, HS, RSMAX, TOPT)

!------------------------------------------------------------------------------------------!
! From the SOILPARM.TBL tables, as functions of soil category.
!------------------------------------------------------------------------------------------!
  REAL(rkind)    :: BEXP         !B parameter ( Assigned in REDPRM )
  REAL(rkind)    :: SMCDRY       !dry soil moisture threshold where direct evap from top
                          !layer ends (volumetric) ( Assigned in REDPRM )
  REAL(rkind)    :: F1           !soil thermal diffusivity/conductivity coef ( Assigned in REDPRM )
  REAL(rkind)    :: SMCMAX       !porosity, saturated value of soil moisture (volumetric)
  REAL(rkind)    :: SMCREF       !reference soil moisture (field capacity) (volumetric) ( Assigned in REDPRM )
  REAL(rkind)    :: PSISAT       !saturated soil matric potential ( Assigned in REDPRM )
  REAL(rkind)    :: DKSAT        !saturated soil hydraulic conductivity ( Assigned in REDPRM )
  REAL(rkind)    :: DWSAT        !saturated soil hydraulic diffusivity ( Assigned in REDPRM )
  REAL(rkind)    :: SMCWLT       !wilting point soil moisture (volumetric) ( Assigned in REDPRM )
  REAL(rkind)    :: QUARTZ       !soil quartz content ( Assigned in REDPRM )

  ! MPC change: make variables private for a given thread
  !$omp threadprivate(BEXP, SMCDRY, F1, SMCMAX, SMCREF, PSISAT, DKSAT, DWSAT, SMCWLT, QUARTZ)

!------------------------------------------------------------------------------------------!
! From the GENPARM.TBL file
!------------------------------------------------------------------------------------------!
  REAL(rkind)    :: SLOPE       !slope index (0 - 1) ( Assigned in REDPRM )
  REAL(rkind)    :: CSOIL       !vol. soil heat capacity [j/m3/K] ( Assigned in REDPRM )
  REAL(rkind)    :: ZBOT        !Depth (m) of lower boundary soil temperature ( Assigned in REDPRM )
  REAL(rkind)    :: CZIL        !Calculate roughness length of heat ( Assigned in REDPRM )

  ! MPC note: FRZK_DATA, REFDK_DATA, and REFKDT_DATA are used in REDPRM to compute KDT and FRZX
  !           (FRZK, REFDK, and REFKDT are local variables within REDPRM and do not need to be thread private)

  REAL(rkind)    :: KDT         !used in compute maximum infiltration rate (in INFIL) ( Assigned in REDPRM )
  REAL(rkind)    :: FRZX        !used in compute maximum infiltration rate (in INFIL) ( Assigned in REDPRM )

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
! 1 -> Jordan (1991); 2 -> BATS: when SFCTMP<TFRZ+2.2 ; 3-> SFCTMP<TFRZ

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
  REAL(rkind), PARAMETER :: TIMEAN = 10.5   !gridcell mean topgraphic index (global mean)
  REAL(rkind), PARAMETER :: FSATMX = 0.38   !maximum surface saturated fraction (global mean)

! adjustable parameters for snow processes

  REAL(rkind), PARAMETER :: M      = 1.0 ! 2.50   !melting factor (-)
  REAL(rkind), PARAMETER :: Z0SNO  = 0.002  !snow surface roughness length (m) (0.002)
  REAL(rkind), PARAMETER :: SSI    = 0.03   !liquid water holding capacity for snowpack (m3/m3) (0.03)
  REAL(rkind), PARAMETER :: SWEMX  = 1.00   !new snow mass to fully cover old snow (mm)
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

    REAL(rkind) :: CH2OP(MVT)       !maximum intercepted h2o per unit lai+sai (mm)
    REAL(rkind) :: DLEAF(MVT)       !characteristic leaf dimension (m)
    REAL(rkind) :: Z0MVT(MVT)       !momentum roughness length (m)
    REAL(rkind) :: HVT(MVT)         !top of canopy (m)
    REAL(rkind) :: HVB(MVT)         !bottom of canopy (m)
    REAL(rkind) :: DEN(MVT)         !tree density (no. of trunks per m2)
    REAL(rkind) :: RC(MVT)          !tree crown radius (m)
    REAL(rkind) :: SAIM(MVT,12)     !monthly stem area index, one-sided
    REAL(rkind) :: LAIM(MVT,12)     !monthly leaf area index, one-sided
    REAL(rkind) :: SLA(MVT)         !single-side leaf area per Kg [m2/kg]
    REAL(rkind) :: DILEFC(MVT)      !coeficient for leaf stress death [1/s]
    REAL(rkind) :: DILEFW(MVT)      !coeficient for leaf stress death [1/s]
    REAL(rkind) :: FRAGR(MVT)       !fraction of growth respiration  !original was 0.3
    REAL(rkind) :: LTOVRC(MVT)      !leaf turnover [1/s]

    REAL(rkind) :: C3PSN(MVT)       !photosynthetic pathway: 0. = c4, 1. = c3
    REAL(rkind) :: KC25(MVT)        !co2 michaelis-menten constant at 25c (pa)
    REAL(rkind) :: AKC(MVT)         !q10 for kc25
    REAL(rkind) :: KO25(MVT)        !o2 michaelis-menten constant at 25c (pa)
    REAL(rkind) :: AKO(MVT)         !q10 for ko25
    REAL(rkind) :: VCMX25(MVT)      !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    REAL(rkind) :: AVCMX(MVT)       !q10 for vcmx25
    REAL(rkind) :: BP(MVT)          !minimum leaf conductance (umol/m**2/s)
    REAL(rkind) :: MP(MVT)          !slope of conductance-to-photosynthesis relationship
    REAL(rkind) :: QE25(MVT)        !quantum efficiency at 25c (umol co2 / umol photon)
    REAL(rkind) :: AQE(MVT)         !q10 for qe25
    REAL(rkind) :: RMF25(MVT)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
    REAL(rkind) :: RMS25(MVT)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
    REAL(rkind) :: RMR25(MVT)       !root maintenance respiration at 25c (umol co2/kg bio/s)
    REAL(rkind) :: ARM(MVT)         !q10 for maintenance respiration
    REAL(rkind) :: FOLNMX(MVT)      !foliage nitrogen concentration when f(n)=1 (%)
    REAL(rkind) :: TMIN(MVT)        !minimum temperature for photosynthesis (k)

    REAL(rkind) :: XL(MVT)          !leaf/stem orientation index
    REAL(rkind) :: RHOL(MVT,MBAND)  !leaf reflectance: 1=vis, 2=nir
    REAL(rkind) :: RHOS(MVT,MBAND)  !stem reflectance: 1=vis, 2=nir
    REAL(rkind) :: TAUL(MVT,MBAND)  !leaf transmittance: 1=vis, 2=nir
    REAL(rkind) :: TAUS(MVT,MBAND)  !stem transmittance: 1=vis, 2=nir

    REAL(rkind) :: MRP(MVT)         !microbial respiration parameter (umol co2 /kg c/ s)
    REAL(rkind) :: CWPVT(MVT)       !empirical canopy wind parameter

    REAL(rkind) :: WRRAT(MVT)       !wood to non-wood ratio
    REAL(rkind) :: WDPOOL(MVT)      !wood pool (switch 1 or 0) depending on woody or not [-]
    REAL(rkind) :: TDLEF(MVT)       !characteristic T for leaf freezing [K]

    INTEGER :: IK,IM
    REAL(rkind) :: TMP10(MVT*MBAND)
    REAL(rkind) :: TMP11(MVT*MBAND)
    REAL(rkind) :: TMP12(MVT*MBAND)
    REAL(rkind) :: TMP13(MVT*MBAND)
    REAL(rkind) :: TMP14(MVT*12)
    REAL(rkind) :: TMP15(MVT*12)
    REAL(rkind) :: TMP16(MVT*5)

    real(rkind) slarea(MVT)
    real(rkind) eps(MVT,5)

CONTAINS
  subroutine read_mp_veg_parameters(FILENAME_VEGTABLE,DATASET_IDENTIFIER)
    implicit none
    character(len=*), intent(in) :: FILENAME_VEGTABLE
    character(len=*), intent(in) :: DATASET_IDENTIFIER
    integer :: ierr

    ! Temporary arrays used in reshaping namelist arrays
    REAL(rkind) :: TMP10(MVT*MBAND)
    REAL(rkind) :: TMP11(MVT*MBAND)
    REAL(rkind) :: TMP12(MVT*MBAND)
    REAL(rkind) :: TMP13(MVT*MBAND)
    REAL(rkind) :: TMP14(MVT*12)
    REAL(rkind) :: TMP15(MVT*12)
    REAL(rkind) :: TMP16(MVT*5)

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

    REAL(rkind) :: ALBSAT(MSC,MBAND)   !saturated soil albedos: 1=vis, 2=nir
    REAL(rkind) :: ALBDRY(MSC,MBAND)   !dry soil albedos: 1=vis, 2=nir
    REAL(rkind) :: ALBICE(MBAND)       !albedo land ice: 1=vis, 2=nir
    REAL(rkind) :: ALBLAK(MBAND)       !albedo frozen lakes: 1=vis, 2=nir
    REAL(rkind) :: OMEGAS(MBAND)       !two-stream parameter omega for snow
    REAL(rkind) :: BETADS              !two-stream parameter betad for snow
    REAL(rkind) :: BETAIS              !two-stream parameter betad for snow
    REAL(rkind) :: EG(2)               !emissivity

! saturated soil albedos: 1=vis, 2=nir
    DATA(ALBSAT(I,1),I=1,8)/0.15,0.11,0.10,0.09,0.08,0.07,0.06,0.05/
    DATA(ALBSAT(I,2),I=1,8)/0.30,0.22,0.20,0.18,0.16,0.14,0.12,0.10/

! dry soil albedos: 1=vis, 2=nir
    DATA(ALBDRY(I,1),I=1,8)/0.27,0.22,0.20,0.18,0.16,0.14,0.12,0.10/
    DATA(ALBDRY(I,2),I=1,8)/0.54,0.44,0.40,0.36,0.32,0.28,0.24,0.20/

! albedo land ice: 1=vis, 2=nir
    DATA (ALBICE(I),I=1,MBAND) /0.80, 0.55/

! albedo frozen lakes: 1=vis, 2=nir
    DATA (ALBLAK(I),I=1,MBAND) /0.60, 0.40/

! omega,betad,betai for snow
    DATA (OMEGAS(I),I=1,MBAND) /0.8, 0.4/
    DATA BETADS, BETAIS /0.5, 0.5/

! emissivity ground surface
      DATA EG /0.97, 0.98/ ! 1-soil;2-lake

END MODULE NOAHMP_RAD_PARAMETERS
! ==================================================================================================

MODULE NOAHMP_ROUTINES
    use nrtype
  USE NOAHMP_GLOBALS
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
  REAL(rkind)                   , INTENT(IN   ) :: SNOWH  !snow height [m]
  REAL(rkind)                   , INTENT(IN   ) :: TV     !vegetation temperature (k)
  REAL(rkind)                   , INTENT(IN   ) :: LAT    !latitude (radians)
  INTEGER                , INTENT(IN   ) :: YEARLEN!Number of days in the particular year
  REAL(rkind)                   , INTENT(IN   ) :: JULIAN !Julian day of year (fractional) ( 0 <= JULIAN < YEARLEN )
  real(rkind)                   , INTENT(IN   ) :: TROOT  !root-zone averaged temperature (k)
  REAL(rkind)                   , INTENT(INOUT) :: LAI    !LAI, unadjusted for burying by snow
  REAL(rkind)                   , INTENT(INOUT) :: SAI    !SAI, unadjusted for burying by snow

! outputs
  REAL(rkind)                   , INTENT(OUT  ) :: HTOP   !top of canopy layer (m)
  REAL(rkind)                   , INTENT(OUT  ) :: ELAI   !leaf area index, after burying by snow
  REAL(rkind)                   , INTENT(OUT  ) :: ESAI   !stem area index, after burying by snow
  REAL(rkind)                   , INTENT(OUT  ) :: IGS    !growing season index (0=off, 1=on)

! locals

  REAL(rkind)                                   :: DB     !thickness of canopy buried by snow (m)
  REAL(rkind)                                   :: FB     !fraction of canopy buried by snow
  REAL(rkind)                                   :: SNOWHC !critical snow depth at which short vege
                                                   !is fully covered by snow

  INTEGER                                :: K       !index
  INTEGER                                :: IT1,IT2 !interpolation months
  REAL(rkind)                                   :: DAY     !current day of year ( 0 <= DAY < YEARLEN )
  REAL(rkind)                                   :: WT1,WT2 !interpolation weights
  REAL(rkind)                                   :: T       !current month (1.00, ..., 12.00)
! --------------------------------------------------------------------------------------------------

  IF ( DVEG == 1 .or. DVEG == 3 .or. DVEG == 4 ) THEN

     IF (LAT >= 0.) THEN
        ! Northern Hemisphere
        DAY = JULIAN
     ELSE
        ! Southern Hemisphere.  DAY is shifted by 1/2 year.
        DAY = MOD ( JULIAN + ( 0.5 * YEARLEN ) , REAL(YEARLEN) )
     ENDIF

     T = 12. * DAY / REAL(YEARLEN)
     IT1 = T + 0.5
     IT2 = IT1 + 1
     WT1 = (IT1+0.5) - T
     WT2 = 1.-WT1
     IF (IT1 .LT.  1) IT1 = 12
     IF (IT2 .GT. 12) IT2 = 1

     LAI = WT1*LAIM(VEGTYP,IT1) + WT2*LAIM(VEGTYP,IT2)
     SAI = WT1*SAIM(VEGTYP,IT1) + WT2*SAIM(VEGTYP,IT2)
  ENDIF
  IF (SAI < 0.05) SAI = 0.0  ! MB: SAI CHECK
  IF (LAI < 0.05 .OR. SAI == 0.0) LAI = 0.0  ! MB: LAI CHECK

  IF ( ( VEGTYP == ISWATER ) .OR. ( VEGTYP == ISBARREN ) .OR. ( VEGTYP == ISSNOW ) .or. ( VEGTYP == ISURBAN) ) THEN
     LAI  = 0.
     SAI  = 0.
  ENDIF

!buried by snow

     DB = MIN( MAX(SNOWH - HVB(VEGTYP),0.), HVT(VEGTYP)-HVB(VEGTYP) )
     FB = DB / MAX(1.E-06,HVT(VEGTYP)-HVB(VEGTYP))
     !print*, 'HVB(VEGTYP), HVT(VEGTYP), DB, FB = ', HVB(VEGTYP), HVT(VEGTYP), DB, FB

     IF(HVT(VEGTYP)> 0. .AND. HVT(VEGTYP) <= 0.5) THEN
       SNOWHC = HVT(VEGTYP)*EXP(-SNOWH/0.1)
       IF (SNOWH>SNOWHC) THEN
        FB = 1.
       ELSE
        FB = SNOWH/SNOWHC
       ENDIF
     ENDIF

     ELAI =  LAI*(1.-FB)
     ESAI =  SAI*(1.-FB)
     IF (ESAI < 0.05) ESAI = 0.0  ! MB: ESAI CHECK
     IF (ELAI < 0.05 .OR. ESAI == 0.0) ELAI = 0.0  ! MB: LAI CHECK

     IF (TV .GT. TMIN(VEGTYP)) THEN
         IGS = 1.
     ELSE
         IGS = 0.
     ENDIF

     HTOP = HVT(VEGTYP)

  END SUBROUTINE PHENOLOGY
! ==================================================================================================
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

  REAL(rkind), INTENT(IN)                     :: DT     !time step [s]
  REAL(rkind), INTENT(IN)                     :: QSNOW  !snowfall (mm/s)
  REAL(rkind), INTENT(IN)                     :: SNEQVO !snow mass at last time step(mm)
  REAL(rkind), INTENT(IN)                     :: SNEQV  !snow mass (mm)
  REAL(rkind), INTENT(IN)                     :: SNOWH  !snow height (mm)
  REAL(rkind), INTENT(IN)                     :: COSZ   !cosine solar zenith angle (0-1)
  REAL(rkind), INTENT(IN)                     :: TG     !ground temperature (k)
  REAL(rkind), INTENT(IN)                     :: TV     !vegetation temperature (k)
  REAL(rkind), INTENT(IN)                     :: ELAI   !LAI, one-sided, adjusted for burying by snow
  REAL(rkind), INTENT(IN)                     :: ESAI   !SAI, one-sided, adjusted for burying by snow
  REAL(rkind), INTENT(IN)                     :: FWET   !fraction of canopy that is wet
  REAL(rkind), DIMENSION(1:NSOIL), INTENT(IN) :: SMC    !volumetric soil water [m3/m3]
  REAL(rkind), DIMENSION(1:2)    , INTENT(IN) :: SOLAD  !incoming direct solar radiation (w/m2)
  REAL(rkind), DIMENSION(1:2)    , INTENT(IN) :: SOLAI  !incoming diffuse solar radiation (w/m2)
  REAL(rkind), INTENT(IN)                     :: FSNO   !snow cover fraction (-)
  REAL(rkind), INTENT(IN)                     :: FVEG   !green vegetation fraction [0.0-1.0]

! inout
  REAL(rkind),                  INTENT(INOUT) :: ALBOLD !snow albedo at last time step (CLASS type)
  REAL(rkind),                  INTENT(INOUT) :: TAUSS  !non-dimensional snow age.

! output
  REAL(rkind), INTENT(OUT)                    :: FSUN   !sunlit fraction of canopy (-)
  REAL(rkind), INTENT(OUT)                    :: LAISUN !sunlit leaf area (-)
  REAL(rkind), INTENT(OUT)                    :: LAISHA !shaded leaf area (-)
  REAL(rkind), INTENT(OUT)                    :: PARSUN !average absorbed par for sunlit leaves (w/m2)
  REAL(rkind), INTENT(OUT)                    :: PARSHA !average absorbed par for shaded leaves (w/m2)
  REAL(rkind), INTENT(OUT)                    :: SAV    !solar radiation absorbed by vegetation (w/m2)
  REAL(rkind), INTENT(OUT)                    :: SAG    !solar radiation absorbed by ground (w/m2)
  REAL(rkind), INTENT(OUT)                    :: FSA    !total absorbed solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                    :: FSR    !total reflected solar radiation (w/m2)

!jref:start
  REAL(rkind), INTENT(OUT)                    :: FSRV    !veg. reflected solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                    :: FSRG    !ground reflected solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                    :: BGAP
  REAL(rkind), INTENT(OUT)                    :: WGAP
!jref:end

! local
  REAL(rkind)                                 :: FAGE   !snow age function (0 - new snow)
  REAL(rkind), DIMENSION(1:2)                 :: ALBGRD !ground albedo (direct)
  REAL(rkind), DIMENSION(1:2)                 :: ALBGRI !ground albedo (diffuse)
  REAL(rkind), DIMENSION(1:2)                 :: ALBD   !surface albedo (direct)
  REAL(rkind), DIMENSION(1:2)                 :: ALBI   !surface albedo (diffuse)
  REAL(rkind), DIMENSION(1:2)                 :: FABD   !flux abs by veg (per unit direct flux)
  REAL(rkind), DIMENSION(1:2)                 :: FABI   !flux abs by veg (per unit diffuse flux)
  REAL(rkind), DIMENSION(1:2)                 :: FTDD   !down direct flux below veg (per unit dir flux)
  REAL(rkind), DIMENSION(1:2)                 :: FTID   !down diffuse flux below veg (per unit dir flux)
  REAL(rkind), DIMENSION(1:2)                 :: FTII   !down diffuse flux below veg (per unit dif flux)
!jref:start
  REAL(rkind), DIMENSION(1:2)                 :: FREVI
  REAL(rkind), DIMENSION(1:2)                 :: FREVD
  REAL(rkind), DIMENSION(1:2)                 :: FREGI
  REAL(rkind), DIMENSION(1:2)                 :: FREGD
!jref:end

  REAL(rkind)                                 :: FSHA   !shaded fraction of canopy
  REAL(rkind)                                 :: VAI    !total LAI + stem area index, one sided

  REAL(rkind),PARAMETER :: MPE = 1.E-6
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

  REAL(rkind),                     INTENT(IN)  :: DT     !time step [sec]
  REAL(rkind),                     INTENT(IN)  :: QSNOW  !snowfall
  REAL(rkind),                     INTENT(IN)  :: COSZ   !cosine solar zenith angle for next time step
  REAL(rkind),                     INTENT(IN)  :: SNOWH  !snow height (mm)
  REAL(rkind),                     INTENT(IN)  :: TG     !ground temperature (k)
  REAL(rkind),                     INTENT(IN)  :: TV     !vegetation temperature (k)
  REAL(rkind),                     INTENT(IN)  :: ELAI   !LAI, one-sided, adjusted for burying by snow
  REAL(rkind),                     INTENT(IN)  :: ESAI   !SAI, one-sided, adjusted for burying by snow
  REAL(rkind),                     INTENT(IN)  :: FSNO   !fraction of grid covered by snow
  REAL(rkind),                     INTENT(IN)  :: FWET   !fraction of canopy that is wet
  REAL(rkind),                     INTENT(IN)  :: SNEQVO !snow mass at last time step(mm)
  REAL(rkind),                     INTENT(IN)  :: SNEQV  !snow mass (mm)
  REAL(rkind),                     INTENT(IN)  :: FVEG   !green vegetation fraction [0.0-1.0]
  REAL(rkind), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    !volumetric soil water (m3/m3)

! inout
  REAL(rkind),                  INTENT(INOUT)  :: ALBOLD !snow albedo at last time step (CLASS type)
  REAL(rkind),                  INTENT(INOUT)  :: TAUSS  !non-dimensional snow age

! output
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBGRD !ground albedo (direct)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBGRI !ground albedo (diffuse)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBD   !surface albedo (direct)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBI   !surface albedo (diffuse)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FABD   !flux abs by veg (per unit direct flux)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FABI   !flux abs by veg (per unit diffuse flux)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FTDD   !down direct flux below veg (per unit dir flux)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FTID   !down diffuse flux below veg (per unit dir flux)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FTII   !down diffuse flux below veg (per unit dif flux)
  REAL(rkind),                     INTENT(OUT) :: FSUN   !sunlit fraction of canopy (-)
!jref:start
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FREVD
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FREVI
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FREGD
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FREGI
  REAL(rkind), INTENT(OUT) :: BGAP
  REAL(rkind), INTENT(OUT) :: WGAP
!jref:end

! ------------------------------------------------------------------------
! ------------------------ local variables -------------------------------
! local
  REAL(rkind)                 :: FAGE     !snow age function
  REAL(rkind)                 :: ALB
  INTEGER              :: IB       !indices
  INTEGER              :: NBAND    !number of solar radiation wave bands
  INTEGER              :: IC       !direct beam: ic=0; diffuse: ic=1

  REAL(rkind)                 :: WL       !fraction of LAI+SAI that is LAI
  REAL(rkind)                 :: WS       !fraction of LAI+SAI that is SAI
  REAL(rkind)                 :: MPE      !prevents overflow for division by zero

  REAL(rkind), DIMENSION(1:2) :: RHO      !leaf/stem reflectance weighted by fraction LAI and SAI
  REAL(rkind), DIMENSION(1:2) :: TAU      !leaf/stem transmittance weighted by fraction LAI and SAI
  REAL(rkind), DIMENSION(1:2) :: FTDI     !down direct flux below veg per unit dif flux = 0
  REAL(rkind), DIMENSION(1:2) :: ALBSND   !snow albedo (direct)
  REAL(rkind), DIMENSION(1:2) :: ALBSNI   !snow albedo (diffuse)

  REAL(rkind)                 :: VAI      !ELAI+ESAI
  REAL(rkind)                 :: GDIR     !average projected leaf/stem area in solar direction
  REAL(rkind)                 :: EXT      !optical depth direct beam per unit leaf + stem area

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

! sunlit fraction of canopy. set FSUN = 0 if FSUN < 0.01.

  EXT = GDIR/COSZ * SQRT(1.-RHO(1)-TAU(1))
  FSUN = (1.-EXP(-EXT*VAI)) / MAX(EXT*VAI,MPE)
  EXT = FSUN

  IF (EXT .LT. 0.01) THEN
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
  REAL(rkind), INTENT(IN)                 :: MPE     !prevents underflow errors if division by zero

  REAL(rkind), INTENT(IN)                 :: FSUN    !sunlit fraction of canopy
  REAL(rkind), INTENT(IN)                 :: FSHA    !shaded fraction of canopy
  REAL(rkind), INTENT(IN)                 :: ELAI    !leaf area, one-sided
  REAL(rkind), INTENT(IN)                 :: VAI     !leaf + stem area, one-sided
  REAL(rkind), INTENT(IN)                 :: LAISUN  !sunlit leaf area index, one-sided
  REAL(rkind), INTENT(IN)                 :: LAISHA  !shaded leaf area index, one-sided

  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: SOLAD   !incoming direct solar radiation (w/m2)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: SOLAI   !incoming diffuse solar radiation (w/m2)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FABD    !flux abs by veg (per unit incoming direct flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FABI    !flux abs by veg (per unit incoming diffuse flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FTDD    !down dir flux below veg (per incoming dir flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FTID    !down dif flux below veg (per incoming dir flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FTII    !down dif flux below veg (per incoming dif flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: ALBGRD  !ground albedo (direct)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: ALBGRI  !ground albedo (diffuse)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: ALBD    !overall surface albedo (direct)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: ALBI    !overall surface albedo (diffuse)

  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FREVD    !overall surface albedo veg (direct)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FREVI    !overall surface albedo veg (diffuse)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FREGD    !overall surface albedo grd (direct)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FREGI    !overall surface albedo grd (diffuse)

! output

  REAL(rkind), INTENT(OUT)                :: PARSUN  !average absorbed par for sunlit leaves (w/m2)
  REAL(rkind), INTENT(OUT)                :: PARSHA  !average absorbed par for shaded leaves (w/m2)
  REAL(rkind), INTENT(OUT)                :: SAV     !solar radiation absorbed by vegetation (w/m2)
  REAL(rkind), INTENT(OUT)                :: SAG     !solar radiation absorbed by ground (w/m2)
  REAL(rkind), INTENT(OUT)                :: FSA     !total absorbed solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                :: FSR     !total reflected solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                :: FSRV    !reflected solar radiation by vegetation
  REAL(rkind), INTENT(OUT)                :: FSRG    !reflected solar radiation by ground

! ------------------------ local variables ----------------------------------------------------
  INTEGER                          :: IB      !waveband number (1=vis, 2=nir)
  INTEGER                          :: NBAND   !number of solar radiation waveband classes

  REAL(rkind)                             :: ABS     !absorbed solar radiation (w/m2)
  REAL(rkind)                             :: RNIR    !reflected solar radiation [nir] (w/m2)
  REAL(rkind)                             :: RVIS    !reflected solar radiation [vis] (w/m2)
  REAL(rkind)                             :: LAIFRA  !leaf area fraction of canopy
  REAL(rkind)                             :: TRD     !transmitted solar radiation: direct (w/m2)
  REAL(rkind)                             :: TRI     !transmitted solar radiation: diffuse (w/m2)
  REAL(rkind), DIMENSION(1:2)             :: CAD     !direct beam absorbed by canopy (w/m2)
  REAL(rkind), DIMENSION(1:2)             :: CAI     !diffuse radiation absorbed by canopy (w/m2)
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
   REAL(rkind), INTENT(IN) :: DT        !main time step (s)
   REAL(rkind), INTENT(IN) :: TG        !ground temperature (k)
   REAL(rkind), INTENT(IN) :: SNEQVO    !snow mass at last time step(mm)
   REAL(rkind), INTENT(IN) :: SNEQV     !snow water per unit ground area (mm)

!output
   REAL(rkind), INTENT(OUT) :: FAGE     !snow age

!input/output
   REAL(rkind), INTENT(INOUT) :: TAUSS      !non-dimensional snow age
!local
   REAL(rkind)            :: TAGE       !total aging effects
   REAL(rkind)            :: AGE1       !effects of grain growth due to vapor diffusion
   REAL(rkind)            :: AGE2       !effects of grain growth at freezing of melt water
   REAL(rkind)            :: AGE3       !effects of soot
   REAL(rkind)            :: DELA       !temporary variable
   REAL(rkind)            :: SGE        !temporary variable
   REAL(rkind)            :: DELS       !temporary variable
   REAL(rkind)            :: DELA0      !temporary variable
   REAL(rkind)            :: ARG        !temporary variable
! See Yang et al. (1997) J.of Climate for detail.
!---------------------------------------------------------------------------------------------------

   IF(SNEQV.LE.0.0) THEN
          TAUSS = 0.
   ELSE IF (SNEQV.GT.800.) THEN
          TAUSS = 0.
   ELSE
          DELA0 = 1.E-6*DT
          ARG   = 5.E3*(1./TFRZ-1./TG)
          AGE1  = EXP(ARG)
          AGE2  = EXP(AMIN1(0.,10.*ARG))
          AGE3  = 0.3
          TAGE  = AGE1+AGE2+AGE3

          DELA  = DELA0*TAGE
          DELS  = AMAX1(0.0,SNEQV-SNEQVO) / SWEMX
          SGE   = (TAUSS+DELA)*(1.0-DELS)
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

  REAL(rkind),INTENT(IN) :: COSZ    !cosine solar zenith angle
  REAL(rkind),INTENT(IN) :: FSNO    !snow cover fraction (-)
  REAL(rkind),INTENT(IN) :: FAGE    !snow age correction

! output

  REAL(rkind), DIMENSION(1:2),INTENT(OUT) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
  REAL(rkind), DIMENSION(1:2),INTENT(OUT) :: ALBSNI !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  INTEGER :: IB          !waveband class

  REAL(rkind) :: FZEN                 !zenith angle correction
  REAL(rkind) :: CF1                  !temperary variable
  REAL(rkind) :: SL2                  !2.*SL
  REAL(rkind) :: SL1                  !1/SL
  REAL(rkind) :: SL                   !adjustable parameter
  REAL(rkind), PARAMETER :: C1 = 0.2  !default in BATS
  REAL(rkind), PARAMETER :: C2 = 0.5  !default in BATS
!  REAL(rkind), PARAMETER :: C1 = 0.2 * 2. ! double the default to match Sleepers River's
!  REAL(rkind), PARAMETER :: C2 = 0.5 * 2. ! snow surface albedo (double aging effects)
! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        ALBSND(1: NBAND) = 0.
        ALBSNI(1: NBAND) = 0.

! when cosz > 0

        SL=2.0
        SL1=1./SL
        SL2=2.*SL
        CF1=((1.+SL1)/(1.+SL2*COSZ)-SL1)
        FZEN=AMAX1(CF1,0.)

        ALBSNI(1)=0.95*(1.-C1*FAGE)
        ALBSNI(2)=0.65*(1.-C2*FAGE)

        ALBSND(1)=ALBSNI(1)+0.4*FZEN*(1.-ALBSNI(1))    !  vis direct
        ALBSND(2)=ALBSNI(2)+0.4*FZEN*(1.-ALBSNI(2))    !  nir direct

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

  REAL(rkind),INTENT(IN) :: QSNOW     !snowfall (mm/s)
  REAL(rkind),INTENT(IN) :: DT        !time step (sec)
  REAL(rkind),INTENT(IN) :: ALBOLD    !snow albedo at last time step

! in & out

  REAL(rkind),                INTENT(INOUT) :: ALB        !
! output

  REAL(rkind), DIMENSION(1:2),INTENT(OUT) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
  REAL(rkind), DIMENSION(1:2),INTENT(OUT) :: ALBSNI !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  INTEGER :: IB          !waveband class

! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        ALBSND(1: NBAND) = 0.
        ALBSNI(1: NBAND) = 0.

! when cosz > 0

         ALB = 0.55 + (ALBOLD-0.55) * EXP(-0.01*DT/3600.)


! 1 mm fresh snow(SWE) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
! here assume 1cm snow depth will fully cover the old snow

         IF (QSNOW > 0.) then
           ALB = ALB + MIN(QSNOW*DT,SWEMX) * (0.84-ALB)/(SWEMX)
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
  REAL(rkind),                     INTENT(IN)  :: FSNO   !fraction of surface covered with snow (-)
  REAL(rkind),                     INTENT(IN)  :: TG     !ground temperature (k)
  REAL(rkind),                     INTENT(IN)  :: COSZ   !cosine solar zenith angle (0-1)
  REAL(rkind), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    !volumetric soil water content (m3/m3)
  REAL(rkind), DIMENSION(1:    2), INTENT(IN)  :: ALBSND !direct beam snow albedo (vis, nir)
  REAL(rkind), DIMENSION(1:    2), INTENT(IN)  :: ALBSNI !diffuse snow albedo (vis, nir)

!output

  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBGRD !ground albedo (direct beam: vis, nir)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBGRI !ground albedo (diffuse: vis, nir)

!local

  INTEGER                               :: IB     !waveband number (1=vis, 2=nir)
  REAL(rkind)                                  :: INC    !soil water correction factor for soil albedo
  REAL(rkind)                                  :: ALBSOD !soil albedo (direct)
  REAL(rkind)                                  :: ALBSOI !soil albedo (diffuse)
! --------------------------------------------------------------------------------------------------

  DO IB = 1, NBAND
        INC = MAX(0.11-0.40*SMC(1), 0.)
        IF (IST .EQ. 1)  THEN                     !soil
           ALBSOD = MIN(ALBSAT(ISC,IB)+INC,ALBDRY(ISC,IB))
           ALBSOI = ALBSOD
        ELSE IF (TG .GT. TFRZ) THEN               !unfrozen lake, wetland
           ALBSOD = 0.06/(MAX(0.01,COSZ)**1.7 + 0.15)
           ALBSOI = 0.06
        ELSE                                      !frozen lake, wetland
           ALBSOD = ALBLAK(IB)
           ALBSOI = ALBSOD
        END IF

! increase desert and semi-desert albedos

        IF (IST .EQ. 1 .AND. ISC .EQ. 9) THEN
           ALBSOD = ALBSOD + 0.10
           ALBSOI = ALBSOI + 0.10
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

   REAL(rkind),                 INTENT(IN)  :: COSZ    !cosine of direct zenith angle (0-1)
   REAL(rkind),                 INTENT(IN)  :: VAI     !one-sided leaf+stem area index (m2/m2)
   REAL(rkind),                 INTENT(IN)  :: FWET    !fraction of lai, sai that is wetted (-)
   REAL(rkind),                 INTENT(IN)  :: T       !surface temperature (k)

   REAL(rkind), DIMENSION(1:2), INTENT(IN)  :: ALBGRD  !direct  albedo of underlying surface (-)
   REAL(rkind), DIMENSION(1:2), INTENT(IN)  :: ALBGRI  !diffuse albedo of underlying surface (-)
   REAL(rkind), DIMENSION(1:2), INTENT(IN)  :: RHO     !leaf+stem reflectance
   REAL(rkind), DIMENSION(1:2), INTENT(IN)  :: TAU     !leaf+stem transmittance
   REAL(rkind),                 INTENT(IN)  :: FVEG    !green vegetation fraction [0.0-1.0]

! output

   REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FAB     !flux abs by veg layer (per unit incoming flux)
   REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FRE     !flux refl above veg layer (per unit incoming flux)
   REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FTD     !down dir flux below veg layer (per unit in flux)
   REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FTI     !down dif flux below veg layer (per unit in flux)
   REAL(rkind),                 INTENT(OUT) :: GDIR    !projected leaf+stem area in solar direction
   REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FREV    !flux reflected by veg layer   (per unit incoming flux)
   REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FREG    !flux reflected by ground (per unit incoming flux)

! local
   REAL(rkind)                              :: OMEGA   !fraction of intercepted radiation that is scattered
   REAL(rkind)                              :: OMEGAL  !omega for leaves
   REAL(rkind)                              :: BETAI   !upscatter parameter for diffuse radiation
   REAL(rkind)                              :: BETAIL  !betai for leaves
   REAL(rkind)                              :: BETAD   !upscatter parameter for direct beam radiation
   REAL(rkind)                              :: BETADL  !betad for leaves
   REAL(rkind)                              :: EXT     !optical depth of direct beam per unit leaf area
   REAL(rkind)                              :: AVMU    !average diffuse optical depth

   REAL(rkind)                              :: COSZI   !0.001 <= cosz <= 1.000
   REAL(rkind)                              :: ASU     !single scattering albedo
   REAL(rkind)                              :: CHIL    ! -0.4 <= xl <= 0.6

   REAL(rkind)                              :: TMP0,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP7,TMP8,TMP9
   REAL(rkind)                              :: P1,P2,P3,P4,S1,S2,U1,U2,U3
   REAL(rkind)                              :: B,C,D,D1,D2,F,H,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10
   REAL(rkind)                              :: PHI1,PHI2,SIGMA
   REAL(rkind)                              :: FTDS,FTIS,FRES
   REAL(rkind)                              :: DENFVEG
   REAL(rkind)                              :: VAI_SPREAD
!jref:start
   REAL(rkind)                              :: FREVEG,FREBAR,FTDVEG,FTIVEG,FTDBAR,FTIBAR
   REAL(rkind)                              :: THETAZ
!jref:end

!  variables for the modified two-stream scheme
!  Niu and Yang (2004), JGR

   REAL(rkind), PARAMETER :: PAI = 3.14159265
   REAL(rkind) :: HD       !crown depth (m)
   REAL(rkind) :: BB       !vertical crown radius (m)
   REAL(rkind) :: THETAP   !angle conversion from SZA
   REAL(rkind) :: FA       !foliage volume density (m-1)
   REAL(rkind) :: NEWVAI   !effective LSAI (-)

   REAL(rkind),INTENT(INOUT) :: BGAP     !between canopy gap fraction for beam (-)
   REAL(rkind),INTENT(INOUT) :: WGAP     !within canopy gap fraction for beam (-)

   REAL(rkind) :: KOPEN    !gap fraction for diffue light (-)
   REAL(rkind) :: GAP      !total gap fraction for beam ( <=1-shafac )

! -----------------------------------------------------------------
! compute within and between gaps
     VAI_SPREAD = VAI
     if(VAI == 0.0) THEN
         GAP     = 1.0
         KOPEN   = 1.0
     ELSE
         IF(OPT_RAD == 1) THEN
           DENFVEG = -LOG(MAX(1.0-FVEG,0.01))/(PAI*RC(VEGTYP)**2)
           HD      = HVT(VEGTYP) - HVB(VEGTYP)
           BB      = 0.5 * HD
           THETAP  = ATAN(BB/RC(VEGTYP) * TAN(ACOS(MAX(0.01,COSZ))) )
           ! BGAP    = EXP(-DEN(VEGTYP) * PAI * RC(VEGTYP)**2/COS(THETAP) )
           BGAP    = EXP(-DENFVEG * PAI * RC(VEGTYP)**2/COS(THETAP) )
           FA      = VAI/(1.33 * PAI * RC(VEGTYP)**3.0 *(BB/RC(VEGTYP))*DENFVEG)
           NEWVAI  = HD*FA
           WGAP    = (1.0-BGAP) * EXP(-0.5*NEWVAI/COSZ)
           GAP     = MIN(1.0-FVEG, BGAP+WGAP)

           KOPEN   = 0.05
         END IF

         IF(OPT_RAD == 2) THEN
           GAP     = 0.0
           KOPEN   = 0.0
         END IF

         IF(OPT_RAD == 3) THEN
           GAP     = 1.0-FVEG
           KOPEN   = 1.0-FVEG
         END IF
     end if

! calculate two-stream parameters OMEGA, BETAD, BETAI, AVMU, GDIR, EXT.
! OMEGA, BETAD, BETAI are adjusted for snow. values for OMEGA*BETAD
! and OMEGA*BETAI are calculated and then divided by the new OMEGA
! because the product OMEGA*BETAI, OMEGA*BETAD is used in solution.
! also, the transmittances and reflectances (TAU, RHO) are linear
! weights of leaf and stem values.

     COSZI  = MAX(0.001, COSZ)
     CHIL   = MIN( MAX(XL(VEGTYP), -0.4), 0.6)
     IF (ABS(CHIL) .LE. 0.01) CHIL = 0.01
     PHI1   = 0.5 - 0.633*CHIL - 0.330*CHIL*CHIL
     PHI2   = 0.877 * (1.-2.*PHI1)
     GDIR   = PHI1 + PHI2*COSZI
     EXT    = GDIR/COSZI
     AVMU   = ( 1. - PHI1/PHI2 * LOG((PHI1+PHI2)/PHI1) ) / PHI2
     OMEGAL = RHO(IB) + TAU(IB)
     TMP0   = GDIR + PHI2*COSZI
     TMP1   = PHI1*COSZI
     ASU    = 0.5*OMEGAL*GDIR/TMP0 * ( 1.-TMP1/TMP0*LOG((TMP1+TMP0)/TMP1) )
     BETADL = (1.+AVMU*EXT)/(OMEGAL*AVMU*EXT)*ASU
     BETAIL = 0.5 * ( RHO(IB)+TAU(IB) + (RHO(IB)-TAU(IB))   &
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
     if ( ABS (SIGMA) < 1.e-6 ) SIGMA = SIGN(1.e-6,REAL(SIGMA))
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
        FTDS = S2                           *(1.0-GAP) + GAP
        FTIS = (H4*S2/SIGMA + H5*S1 + H6/S1)*(1.0-GAP)
     ELSE
        FTDS = 0.
        FTIS = (H9*S1 + H10/S1)*(1.0-KOPEN) + KOPEN
     END IF
     FTD(IB) = FTDS
     FTI(IB) = FTIS

! flux reflected by the surface (veg. and ground)

     IF (IC .EQ. 0) THEN
        FRES   = (H1/SIGMA + H2 + H3)*(1.0-GAP  ) + ALBGRD(IB)*GAP
        FREVEG = (H1/SIGMA + H2 + H3)*(1.0-GAP  )
        FREBAR = ALBGRD(IB)*GAP                   !jref - separate veg. and ground reflection
     ELSE
        FRES   = (H7 + H8) *(1.0-KOPEN) + ALBGRI(IB)*KOPEN
        FREVEG = (H7 + H8) *(1.0-KOPEN) + ALBGRI(IB)*KOPEN
        FREBAR = 0                                !jref - separate veg. and ground reflection
     END IF
     FRE(IB) = FRES

     FREV(IB) = FREVEG
     FREG(IB) = FREBAR

! flux absorbed by vegetation

     FAB(IB) = 1. - FRE(IB) - (1.-ALBGRD(IB))*FTD(IB) &
                            - (1.-ALBGRI(IB))*FTI(IB)

!if(iloc == 1.and.jloc ==  2) then
  !write(*,'(a7,2i2,5(a6,f8.4),2(a9,f8.4))') "ib,ic: ",ib,ic," GAP: ",GAP," FTD: ",FTD(IB)," FTI: ",FTI(IB)," FRE: ", &
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
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
      INTEGER,INTENT(IN)  :: ILOC   !grid index
      INTEGER,INTENT(IN)  :: JLOC   !grid index
      INTEGER,INTENT(IN)  :: VEGTYP !vegetation physiology type

      REAL(rkind), INTENT(IN)    :: IGS    !growing season index (0=off, 1=on)
      REAL(rkind), INTENT(IN)    :: MPE    !prevents division by zero errors

      REAL(rkind), INTENT(IN)    :: TV     !foliage temperature (k)
      REAL(rkind), INTENT(IN)    :: EI     !vapor pressure inside leaf (sat vapor press at tv) (pa)
      REAL(rkind), INTENT(IN)    :: EA     !vapor pressure of canopy air (pa)
      REAL(rkind), INTENT(IN)    :: APAR   !par absorbed per unit lai (w/m2)
      REAL(rkind), INTENT(IN)    :: O2     !atmospheric o2 concentration (pa)
      REAL(rkind), INTENT(IN)    :: CO2    !atmospheric co2 concentration (pa)
      REAL(rkind), INTENT(IN)    :: SFCPRS !air pressure at reference height (pa)
      REAL(rkind), INTENT(IN)    :: SFCTMP !air temperature at reference height (k)
      REAL(rkind), INTENT(IN)    :: BTRAN  !soil water transpiration factor (0 to 1)
      REAL(rkind), INTENT(IN)    :: FOLN   !foliage nitrogen concentration (%)
      REAL(rkind), INTENT(IN)    :: RB     !boundary layer resistance (s/m)

! output
      REAL(rkind), INTENT(OUT)   :: RS     !leaf stomatal resistance (s/m)
      REAL(rkind), INTENT(OUT)   :: PSN    !foliage photosynthesis (umol co2 /m2/ s) [always +]

! in&out
      REAL(rkind)                :: RLB    !boundary layer resistance (s m2 / umol)
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
      INTEGER :: ITER     !iteration index
      INTEGER :: NITER    !number of iterations

      DATA NITER /3/
      SAVE NITER

      REAL(rkind) :: AB          !used in statement functions
      REAL(rkind) :: BC          !used in statement functions
      REAL(rkind) :: F1          !generic temperature response (statement function)
      REAL(rkind) :: F2          !generic temperature inhibition (statement function)
      REAL(rkind) :: TC          !foliage temperature (degree Celsius)
      REAL(rkind) :: CS          !co2 concentration at leaf surface (pa)
      REAL(rkind) :: KC          !co2 Michaelis-Menten constant (pa)
      REAL(rkind) :: KO          !o2 Michaelis-Menten constant (pa)
      REAL(rkind) :: A,B,C,Q     !intermediate calculations for RS
      REAL(rkind) :: R1,R2       !roots for RS
      REAL(rkind) :: FNF         !foliage nitrogen adjustment factor (0 to 1)
      REAL(rkind) :: PPF         !absorb photosynthetic photon flux (umol photons/m2/s)
      REAL(rkind) :: WC          !Rubisco limited photosynthesis (umol co2/m2/s)
      REAL(rkind) :: WJ          !light limited photosynthesis (umol co2/m2/s)
      REAL(rkind) :: WE          !export limited photosynthesis (umol co2/m2/s)
      REAL(rkind) :: CP          !co2 compensation point (pa)
      REAL(rkind) :: CI          !internal co2 (pa)
      REAL(rkind) :: AWC         !intermediate calculation for wc
      REAL(rkind) :: VCMX        !maximum rate of carbonylation (umol co2/m2/s)
      REAL(rkind) :: J           !electron transport (umol co2/m2/s)
      REAL(rkind) :: CEA         !constrain ea or else model blows up
      REAL(rkind) :: CF          !s m2/umol -> s/m

      F1(AB,BC) = AB**((BC-25.)/10.)
      F2(AB) = 1. + EXP((-2.2E05+710.*(AB+273.16))/(8.314*(AB+273.16)))
      REAL(rkind) :: T
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
         CP  = 0.5*KC/KO*O2*0.21
         VCMX = VCMX25(VEGTYP) / F2(TC) * FNF * BTRAN * F1(AVCMX(VEGTYP),TC)
         !write(*,'(a,1x,20(f14.8,1x))') 'F1, F2, VCMX25(VEGTYP), AVCMX(VEGTYP), TC = ', F1(AVCMX(VEGTYP),TC), F2(TC), VCMX25(VEGTYP), AVCMX(VEGTYP), TC

! first guess ci

         CI = 0.7*CO2*C3PSN(VEGTYP) + 0.4*CO2*(1.-C3PSN(VEGTYP))

         !write(*,'(a,1x,10(f20.10,1x))') 'KC25(VEGTYP), AKC(VEGTYP), KO25(VEGTYP), AKO(VEGTYP) = ', KC25(VEGTYP), AKC(VEGTYP), KO25(VEGTYP), AKO(VEGTYP)
         !write(*,'(a,1x,10(f20.10,1x))') 'CO2, CI, CP, KC, KO = ', CO2, CI, CP, KC, KO

! rb: s/m -> s m**2 / umol

         RLB = RB/CF

! constrain ea

         CEA = MAX(0.25*EI*C3PSN(VEGTYP)+0.40*EI*(1.-C3PSN(VEGTYP)), MIN(EA,EI) )

         !print*, '**'
         !write(*,'(a,1x,20(f14.8,1x))') 'BTRAN, VCMX, MP(VEGTYP), EA, EI, CEA/EI, O2, CO2, KC, KO, J, TC, RLB, SFCPRS = ', &
         !                                BTRAN, VCMX, MP(VEGTYP), EA, EI, CEA/EI, O2, CO2, KC, KO, J, TC, RLB, SFCPRS

! ci iteration
!jref: C3PSN is equal to 1 for all veg types.
       DO ITER = 1, NITER

            WJ = MAX(CI-CP,0.)*J/(CI+2.*CP)*C3PSN(VEGTYP)  + J*(1.-C3PSN(VEGTYP))
            WC = MAX(CI-CP,0.)*VCMX/(CI+AWC)*C3PSN(VEGTYP) + VCMX*(1.-C3PSN(VEGTYP))
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

            !write(*,'(a,1x10(f14.10,1x))') 'WJ, WC, WE, PSN, CI, CS, RS, TV, VCMX, J = ', &
            !                                WJ, WC, WE, PSN, CI, CS, RS, TV, VCMX, J

       END DO
       !pause

! rs, rb:  s m**2 / umol -> s/m

         RS = RS*CF
         !write(*,'(a,1x,10(f20.10,1x))') 'RS, G, CEA, EA, EI = ', RS, 1./RS, CEA, EA, EI

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
    REAL(rkind),                     INTENT(IN)  :: PAR    !par absorbed per unit sunlit lai (w/m2)
    REAL(rkind),                     INTENT(IN)  :: SFCTMP !canopy air temperature
    REAL(rkind),                     INTENT(IN)  :: SFCPRS !surface pressure (pa)
    REAL(rkind),                     INTENT(IN)  :: EAH    !water vapor pressure (pa)
    REAL(rkind),                     INTENT(IN)  :: RCSOIL !soil moisture stress factor

!outputs

    REAL(rkind),                     INTENT(OUT) :: RC     !canopy resistance per unit LAI
    REAL(rkind),                     INTENT(OUT) :: PSN    !foliage photosynthesis (umolco2/m2/s)

!local

    REAL(rkind)                                  :: RCQ
    REAL(rkind)                                  :: RCS
    REAL(rkind)                                  :: RCT
    REAL(rkind)                                  :: FF
    REAL(rkind)                                  :: Q2     !water vapor mixing ratio (kg/kg)
    REAL(rkind)                                  :: Q2SAT  !saturation Q2
    REAL(rkind)                                  :: DQSDT2 !d(Q2SAT)/d(T)

! RSMIN, RSMAX, TOPT, RGL, HS are canopy stress parameters set in REDPRM
! ----------------------------------------------------------------------
! initialize canopy resistance multiplier terms.
! ----------------------------------------------------------------------
    RC     = 0.0
    RCS    = 0.0
    RCT    = 0.0
    RCQ    = 0.0

!  compute Q2 and Q2SAT

    Q2 = 0.622 *  EAH  / (SFCPRS - 0.378 * EAH) !specific humidity [kg/kg]
    Q2 = Q2 / (1.0 + Q2)                        !mixing ratio [kg/kg]

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
        REAL(rkind), INTENT(IN)       :: SFCTMP, SFCPRS
        REAL(rkind), INTENT(OUT)      :: Q2SAT, DQSDT2
        REAL(rkind), PARAMETER        :: A2=17.67,A3=273.15,A4=29.65, ELWV=2.501E6,         &
                                  A23M4=A2*(A3-A4), E0=0.611, RV=461.0,             &
                                  EPSILON=0.622
        REAL(rkind)                   :: ES, SFCPRSX
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
!         F11 = ALOG10(SATPSI) + BB*ALOG10(MAXSMC) + 2.0
!       REFSMC1=MAXSMC*(5.79E-9/SATDK)**(1/(2*BB+3)) 5.79E-9 m/s= 0.5 mm
!       REFSMC=REFSMC1+1./3.(MAXSMC-REFSMC1)
!       WLTSMC1=MAXSMC*(200./SATPSI)**(-1./BB)    (Wetzel and Chang, 198
!       WLTSMC=WLTSMC1-0.5*WLTSMC1
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
    REAL(rkind),DIMENSION(NSOIL),INTENT(IN) :: ZSOIL

! Locals
    REAL(rkind)    :: REFDK
    REAL(rkind)    :: REFKDT
    REAL(rkind)    :: FRZK
    REAL(rkind)                   :: FRZFACT
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
       SMCMAX = 0.45
       SMCREF = 0.42
       SMCWLT = 0.40
       SMCDRY = 0.40
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
      FRZFACT = (SMCMAX / SMCREF) * (0.412 / 0.468)
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
       RSMIN=400.0
    ENDIF

!    SHDFAC = SHDTBL(VEGTYP)
!    IF (VEGTYP .eq. BARE) SHDFAC = 0.0

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
    use nrtype

  USE NOAHMP_ROUTINES
  USE NOAHMP_GLOBALS
  USE NOAHMP_VEG_PARAMETERS

END MODULE MODULE_SF_NOAHMPLSM
