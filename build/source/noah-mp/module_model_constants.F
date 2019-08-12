!WRF:MODEL_LAYER:CONSTANTS
!

 MODULE module_model_constants
   use nrtype
   !  2. Following are constants for use in defining real number bounds.

   !  A really small number.

   REAL(DP)    , PARAMETER :: epsilon         = 1.E-15

   !  4. Following is information related to the physical constants.

   !  These are the physical constants used within the model.

! JM NOTE -- can we name this grav instead?
   REAL(DP)    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)

#if ( NMM_CORE == 1 )
   REAL(DP)    , PARAMETER :: r_d          = 287.04
   REAL(DP)    , PARAMETER :: cp           = 1004.6
#else
   REAL(DP)    , PARAMETER :: r_d          = 287.
   REAL(DP)    , PARAMETER :: cp           = 7.*r_d/2.
#endif

   REAL(DP)    , PARAMETER :: r_v          = 461.6
   REAL(DP)    , PARAMETER :: cv           = cp-r_d
   REAL(DP)    , PARAMETER :: cpv          = 4.*r_v
   REAL(DP)    , PARAMETER :: cvv          = cpv-r_v
   REAL(DP)    , PARAMETER :: cvpm         = -cv/cp
   REAL(DP)    , PARAMETER :: cliq         = 4190.
   REAL(DP)    , PARAMETER :: cice         = 2106.
   REAL(DP)    , PARAMETER :: psat         = 610.78
   REAL(DP)    , PARAMETER :: rcv          = r_d/cv
   REAL(DP)    , PARAMETER :: rcp          = r_d/cp
   REAL(DP)    , PARAMETER :: rovg         = r_d/g
   REAL(DP)    , PARAMETER :: c2           = cp * rcv
   real    , parameter :: mwdry        = 28.966 ! molecular weight of dry air (g/mole)

   REAL(DP)    , PARAMETER :: p1000mb      = 100000.
   REAL(DP)    , PARAMETER :: t0           = 300.
   REAL(DP)    , PARAMETER :: p0           = p1000mb
   REAL(DP)    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL(DP)    , PARAMETER :: cvovcp       = 1./cpovcv
   REAL(DP)    , PARAMETER :: rvovrd       = r_v/r_d

   REAL(DP)    , PARAMETER :: reradius     = 1./6370.0e03

   REAL(DP)    , PARAMETER :: asselin      = .025
!   REAL(DP)    , PARAMETER :: asselin      = .0
   REAL(DP)    , PARAMETER :: cb           = 25.

   REAL(DP)    , PARAMETER :: XLV0         = 3.15E6
   REAL(DP)    , PARAMETER :: XLV1         = 2370.
   REAL(DP)    , PARAMETER :: XLS0         = 2.905E6
   REAL(DP)    , PARAMETER :: XLS1         = 259.532

   REAL(DP)    , PARAMETER :: XLS          = 2.85E6
   REAL(DP)    , PARAMETER :: XLV          = 2.5E6
   REAL(DP)    , PARAMETER :: XLF          = 3.50E5

   REAL(DP)    , PARAMETER :: rhowater     = 1000.
   REAL(DP)    , PARAMETER :: rhosnow      = 100.
   REAL(DP)    , PARAMETER :: rhoair0      = 1.28
!
   REAL(DP)    , PARAMETER :: n_ccn0       = 1.0E8
!
   REAL(DP)    , PARAMETER :: DEGRAD       = 3.1415926/180.
   REAL(DP)    , PARAMETER :: DPD          = 360./365.

   REAL(DP)    , PARAMETER ::  SVP1=0.6112
   REAL(DP)    , PARAMETER ::  SVP2=17.67
   REAL(DP)    , PARAMETER ::  SVP3=29.65
   REAL(DP)    , PARAMETER ::  SVPT0=273.15
   REAL(DP)    , PARAMETER ::  EP_1=R_v/R_d-1.
   REAL(DP)    , PARAMETER ::  EP_2=R_d/R_v
   REAL(DP)    , PARAMETER ::  KARMAN=0.4
   REAL(DP)    , PARAMETER ::  EOMEG=7.2921E-5
   REAL(DP)    , PARAMETER ::  STBOLT=5.67051E-8

   REAL(DP)    , PARAMETER ::  prandtl = 1./3.0
                                         ! constants for w-damping option
   REAL(DP)    , PARAMETER ::  w_alpha = 0.3 ! strength m/s/s
   REAL(DP)    , PARAMETER ::  w_beta  = 1.0 ! activation cfl number

       REAL(DP) , PARAMETER ::  pq0=379.90516
       REAL(DP) , PARAMETER ::  epsq2=0.2
       REAL(DP) , PARAMETER ::  a2=17.2693882
       REAL(DP) , PARAMETER ::  a3=273.16
       REAL(DP) , PARAMETER ::  a4=35.86
       REAL(DP) , PARAMETER ::  epsq=1.e-12
       REAL(DP) , PARAMETER ::  p608=rvovrd-1.
!#if ( NMM_CORE == 1 )
       REAL(DP) , PARAMETER ::  climit=1.e-20
       REAL(DP) , PARAMETER ::  cm1=2937.4
       REAL(DP) , PARAMETER ::  cm2=4.9283
       REAL(DP) , PARAMETER ::  cm3=23.5518
!       REAL(DP) , PARAMETER ::  defc=8.0
!       REAL(DP) , PARAMETER ::  defm=32.0
       REAL(DP) , PARAMETER ::  defc=0.0
       REAL(DP) , PARAMETER ::  defm=99999.0
       REAL(DP) , PARAMETER ::  epsfc=1./1.05
       REAL(DP) , PARAMETER ::  epswet=0.0
       REAL(DP) , PARAMETER ::  fcdif=1./3.
#ifdef HWRF
       REAL(DP) , PARAMETER ::  fcm=0.0
#else
       REAL(DP) , PARAMETER ::  fcm=0.00003
#endif
       REAL(DP) , PARAMETER ::  gma=-r_d*(1.-rcp)*0.5
       REAL(DP) , PARAMETER ::  p400=40000.0
       REAL(DP) , PARAMETER ::  phitp=15000.0
       REAL(DP) , PARAMETER ::  pi2=2.*3.1415926
       REAL(DP) , PARAMETER ::  plbtm=105000.0
       REAL(DP) , PARAMETER ::  plomd=64200.0
       REAL(DP) , PARAMETER ::  pmdhi=35000.0
       REAL(DP) , PARAMETER ::  q2ini=0.50
       REAL(DP) , PARAMETER ::  rfcp=0.25/cp
       REAL(DP) , PARAMETER ::  rhcrit_land=0.75
       REAL(DP) , PARAMETER ::  rhcrit_sea=0.80
       REAL(DP) , PARAMETER ::  rlag=14.8125
       REAL(DP) , PARAMETER ::  rlx=0.90
       REAL(DP) , PARAMETER ::  scq2=50.0
       REAL(DP) , PARAMETER ::  slopht=0.001
       REAL(DP) , PARAMETER ::  tlc=2.*0.703972477
       REAL(DP) , PARAMETER ::  wa=0.15
       REAL(DP) , PARAMETER ::  wght=0.35
       REAL(DP) , PARAMETER ::  wpc=0.075
       REAL(DP) , PARAMETER ::  z0land=0.10
#ifdef HWRF
       REAL(DP) , PARAMETER ::  z0max=0.01
#else
       REAL(DP) , PARAMETER ::  z0max=0.008
#endif
       REAL(DP) , PARAMETER ::  z0sea=0.001
!#endif


   !  Earth

   !  The value for P2SI *must* be set to 1.0 for Earth
   !  Although, now we may not need this declaration here (see above)
   !REAL(DP)    , PARAMETER :: P2SI         = 1.0

   !  Orbital constants:

   INTEGER , PARAMETER :: PLANET_YEAR = 365
   REAL(DP) , PARAMETER :: OBLIQUITY = 23.5
   REAL(DP) , PARAMETER :: ECCENTRICITY = 0.014
   REAL(DP) , PARAMETER :: SEMIMAJORAXIS = 1.0 ! In AU
   ! Don't know the following values, so we'll fake them for now
   REAL(DP) , PARAMETER :: zero_date = 0.0   ! Time of perihelion passage
   !  Fraction into the year (from perhelion) of the
   !  occurrence of the Northern Spring Equinox
   REAL(DP) , PARAMETER :: EQUINOX_FRACTION= 0.0

 CONTAINS
   SUBROUTINE init_module_model_constants
   END SUBROUTINE init_module_model_constants
 END MODULE module_model_constants
