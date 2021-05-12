!WRF:MODEL_LAYER:CONSTANTS
!

 MODULE module_model_constants
  USE nrtype

   !  2. Following are constants for use in defining real number bounds.

   !  A really small number.

   REAL(SUMMA_PREC)    , PARAMETER :: epsilon         = 1.E-15

   !  4. Following is information related to the physical constants.

   !  These are the physical constants used within the model.

! JM NOTE -- can we name this grav instead?
   REAL(SUMMA_PREC)    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)

#if ( NMM_CORE == 1 )
   REAL(SUMMA_PREC)    , PARAMETER :: r_d          = 287.04
   REAL(SUMMA_PREC)    , PARAMETER :: cp           = 1004.6
#else
   REAL(SUMMA_PREC)    , PARAMETER :: r_d          = 287.
   REAL(SUMMA_PREC)    , PARAMETER :: cp           = 7.*r_d/2.
#endif

   REAL(SUMMA_PREC)    , PARAMETER :: r_v          = 461.6
   REAL(SUMMA_PREC)    , PARAMETER :: cv           = cp-r_d
   REAL(SUMMA_PREC)    , PARAMETER :: cpv          = 4.*r_v
   REAL(SUMMA_PREC)    , PARAMETER :: cvv          = cpv-r_v
   REAL(SUMMA_PREC)    , PARAMETER :: cvpm         = -cv/cp
   REAL(SUMMA_PREC)    , PARAMETER :: cliq         = 4190.
   REAL(SUMMA_PREC)    , PARAMETER :: cice         = 2106.
   REAL(SUMMA_PREC)    , PARAMETER :: psat         = 610.78
   REAL(SUMMA_PREC)    , PARAMETER :: rcv          = r_d/cv
   REAL(SUMMA_PREC)    , PARAMETER :: rcp          = r_d/cp
   REAL(SUMMA_PREC)    , PARAMETER :: rovg         = r_d/g
   REAL(SUMMA_PREC)    , PARAMETER :: c2           = cp * rcv
   real    , parameter :: mwdry        = 28.966 ! molecular weight of dry air (g/mole)

   REAL(SUMMA_PREC)    , PARAMETER :: p1000mb      = 100000.
   REAL(SUMMA_PREC)    , PARAMETER :: t0           = 300.
   REAL(SUMMA_PREC)    , PARAMETER :: p0           = p1000mb
   REAL(SUMMA_PREC)    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL(SUMMA_PREC)    , PARAMETER :: cvovcp       = 1./cpovcv
   REAL(SUMMA_PREC)    , PARAMETER :: rvovrd       = r_v/r_d

   REAL(SUMMA_PREC)    , PARAMETER :: reradius     = 1./6370.0e03

   REAL(SUMMA_PREC)    , PARAMETER :: asselin      = .025
!   REAL(SUMMA_PREC)    , PARAMETER :: asselin      = .0
   REAL(SUMMA_PREC)    , PARAMETER :: cb           = 25.

   REAL(SUMMA_PREC)    , PARAMETER :: XLV0         = 3.15E6
   REAL(SUMMA_PREC)    , PARAMETER :: XLV1         = 2370.
   REAL(SUMMA_PREC)    , PARAMETER :: XLS0         = 2.905E6
   REAL(SUMMA_PREC)    , PARAMETER :: XLS1         = 259.532

   REAL(SUMMA_PREC)    , PARAMETER :: XLS          = 2.85E6
   REAL(SUMMA_PREC)    , PARAMETER :: XLV          = 2.5E6
   REAL(SUMMA_PREC)    , PARAMETER :: XLF          = 3.50E5

   REAL(SUMMA_PREC)    , PARAMETER :: rhowater     = 1000.
   REAL(SUMMA_PREC)    , PARAMETER :: rhosnow      = 100.
   REAL(SUMMA_PREC)    , PARAMETER :: rhoair0      = 1.28
!
   REAL(SUMMA_PREC)    , PARAMETER :: n_ccn0       = 1.0E8
!
   REAL(SUMMA_PREC)    , PARAMETER :: DEGRAD       = 3.1415926/180.
   REAL(SUMMA_PREC)    , PARAMETER :: DPD          = 360./365.

   REAL(SUMMA_PREC)    , PARAMETER ::  SVP1=0.6112
   REAL(SUMMA_PREC)    , PARAMETER ::  SVP2=17.67
   REAL(SUMMA_PREC)    , PARAMETER ::  SVP3=29.65
   REAL(SUMMA_PREC)    , PARAMETER ::  SVPT0=273.15
   REAL(SUMMA_PREC)    , PARAMETER ::  EP_1=R_v/R_d-1.
   REAL(SUMMA_PREC)    , PARAMETER ::  EP_2=R_d/R_v
   REAL(SUMMA_PREC)    , PARAMETER ::  KARMAN=0.4
   REAL(SUMMA_PREC)    , PARAMETER ::  EOMEG=7.2921E-5
   REAL(SUMMA_PREC)    , PARAMETER ::  STBOLT=5.67051E-8

   REAL(SUMMA_PREC)    , PARAMETER ::  prandtl = 1./3.0
                                         ! constants for w-damping option
   REAL(SUMMA_PREC)    , PARAMETER ::  w_alpha = 0.3 ! strength m/s/s
   REAL(SUMMA_PREC)    , PARAMETER ::  w_beta  = 1.0 ! activation cfl number

       REAL(SUMMA_PREC) , PARAMETER ::  pq0=379.90516
       REAL(SUMMA_PREC) , PARAMETER ::  epsq2=0.2
       REAL(SUMMA_PREC) , PARAMETER ::  a2=17.2693882
       REAL(SUMMA_PREC) , PARAMETER ::  a3=273.16
       REAL(SUMMA_PREC) , PARAMETER ::  a4=35.86
       REAL(SUMMA_PREC) , PARAMETER ::  epsq=1.e-12
       REAL(SUMMA_PREC) , PARAMETER ::  p608=rvovrd-1.
!#if ( NMM_CORE == 1 )
       REAL(SUMMA_PREC) , PARAMETER ::  climit=1.e-20
       REAL(SUMMA_PREC) , PARAMETER ::  cm1=2937.4
       REAL(SUMMA_PREC) , PARAMETER ::  cm2=4.9283
       REAL(SUMMA_PREC) , PARAMETER ::  cm3=23.5518
!       REAL(SUMMA_PREC) , PARAMETER ::  defc=8.0
!       REAL(SUMMA_PREC) , PARAMETER ::  defm=32.0
       REAL(SUMMA_PREC) , PARAMETER ::  defc=0.0
       REAL(SUMMA_PREC) , PARAMETER ::  defm=99999.0
       REAL(SUMMA_PREC) , PARAMETER ::  epsfc=1./1.05
       REAL(SUMMA_PREC) , PARAMETER ::  epswet=0.0
       REAL(SUMMA_PREC) , PARAMETER ::  fcdif=1./3.
#ifdef HWRF
       REAL(SUMMA_PREC) , PARAMETER ::  fcm=0.0
#else
       REAL(SUMMA_PREC) , PARAMETER ::  fcm=0.00003
#endif
       REAL(SUMMA_PREC) , PARAMETER ::  gma=-r_d*(1.-rcp)*0.5
       REAL(SUMMA_PREC) , PARAMETER ::  p400=40000.0
       REAL(SUMMA_PREC) , PARAMETER ::  phitp=15000.0
       REAL(SUMMA_PREC) , PARAMETER ::  pi2=2.*3.1415926
       REAL(SUMMA_PREC) , PARAMETER ::  plbtm=105000.0
       REAL(SUMMA_PREC) , PARAMETER ::  plomd=64200.0
       REAL(SUMMA_PREC) , PARAMETER ::  pmdhi=35000.0
       REAL(SUMMA_PREC) , PARAMETER ::  q2ini=0.50
       REAL(SUMMA_PREC) , PARAMETER ::  rfcp=0.25/cp
       REAL(SUMMA_PREC) , PARAMETER ::  rhcrit_land=0.75
       REAL(SUMMA_PREC) , PARAMETER ::  rhcrit_sea=0.80
       REAL(SUMMA_PREC) , PARAMETER ::  rlag=14.8125
       REAL(SUMMA_PREC) , PARAMETER ::  rlx=0.90
       REAL(SUMMA_PREC) , PARAMETER ::  scq2=50.0
       REAL(SUMMA_PREC) , PARAMETER ::  slopht=0.001
       REAL(SUMMA_PREC) , PARAMETER ::  tlc=2.*0.703972477
       REAL(SUMMA_PREC) , PARAMETER ::  wa=0.15
       REAL(SUMMA_PREC) , PARAMETER ::  wght=0.35
       REAL(SUMMA_PREC) , PARAMETER ::  wpc=0.075
       REAL(SUMMA_PREC) , PARAMETER ::  z0land=0.10
#ifdef HWRF
       REAL(SUMMA_PREC) , PARAMETER ::  z0max=0.01
#else
       REAL(SUMMA_PREC) , PARAMETER ::  z0max=0.008
#endif
       REAL(SUMMA_PREC) , PARAMETER ::  z0sea=0.001
!#endif


   !  Earth

   !  The value for P2SI *must* be set to 1.0 for Earth
   !  Although, now we may not need this declaration here (see above)
   !REAL(SUMMA_PREC)    , PARAMETER :: P2SI         = 1.0

   !  Orbital constants:

   INTEGER , PARAMETER :: PLANET_YEAR = 365
   REAL(SUMMA_PREC) , PARAMETER :: OBLIQUITY = 23.5
   REAL(SUMMA_PREC) , PARAMETER :: ECCENTRICITY = 0.014
   REAL(SUMMA_PREC) , PARAMETER :: SEMIMAJORAXIS = 1.0 ! In AU
   ! Don't know the following values, so we'll fake them for now
   REAL(SUMMA_PREC) , PARAMETER :: zero_date = 0.0   ! Time of perihelion passage
   !  Fraction into the year (from perhelion) of the
   !  occurrence of the Northern Spring Equinox
   REAL(SUMMA_PREC) , PARAMETER :: EQUINOX_FRACTION= 0.0

 CONTAINS
   SUBROUTINE init_module_model_constants
   END SUBROUTINE init_module_model_constants
 END MODULE module_model_constants
