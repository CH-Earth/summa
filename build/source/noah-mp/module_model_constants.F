!WRF:MODEL_LAYER:CONSTANTS
!

 MODULE module_model_constants
  USE nrtype

   !  2. Following are constants for use in defining real number bounds.

   !  A really small number.

   REAL(rkind)    , PARAMETER :: epsilon         = 1.E-15

   !  4. Following is information related to the physical constants.

   !  These are the physical constants used within the model.

! JM NOTE -- can we name this grav instead?
   REAL(rkind)    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)

#if ( NMM_CORE == 1 )
   REAL(rkind)    , PARAMETER :: r_d          = 287.04
   REAL(rkind)    , PARAMETER :: cp           = 1004.6
#else
   REAL(rkind)    , PARAMETER :: r_d          = 287.
   REAL(rkind)    , PARAMETER :: cp           = 7.*r_d/2.
#endif

   REAL(rkind)    , PARAMETER :: r_v          = 461.6
   REAL(rkind)    , PARAMETER :: cv           = cp-r_d
   REAL(rkind)    , PARAMETER :: cpv          = 4.*r_v
   REAL(rkind)    , PARAMETER :: cvv          = cpv-r_v
   REAL(rkind)    , PARAMETER :: cvpm         = -cv/cp
   REAL(rkind)    , PARAMETER :: cliq         = 4190.
   REAL(rkind)    , PARAMETER :: cice         = 2106.
   REAL(rkind)    , PARAMETER :: psat         = 610.78
   REAL(rkind)    , PARAMETER :: rcv          = r_d/cv
   REAL(rkind)    , PARAMETER :: rcp          = r_d/cp
   REAL(rkind)    , PARAMETER :: rovg         = r_d/g
   REAL(rkind)    , PARAMETER :: c2           = cp * rcv
   real    , parameter :: mwdry        = 28.966 ! molecular weight of dry air (g/mole)

   REAL(rkind)    , PARAMETER :: p1000mb      = 100000.
   REAL(rkind)    , PARAMETER :: t0           = 300.
   REAL(rkind)    , PARAMETER :: p0           = p1000mb
   REAL(rkind)    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL(rkind)    , PARAMETER :: cvovcp       = 1./cpovcv
   REAL(rkind)    , PARAMETER :: rvovrd       = r_v/r_d

   REAL(rkind)    , PARAMETER :: reradius     = 1./6370.0e03

   REAL(rkind)    , PARAMETER :: asselin      = .025
!   REAL(rkind)    , PARAMETER :: asselin      = .0
   REAL(rkind)    , PARAMETER :: cb           = 25.

   REAL(rkind)    , PARAMETER :: XLV0         = 3.15E6
   REAL(rkind)    , PARAMETER :: XLV1         = 2370.
   REAL(rkind)    , PARAMETER :: XLS0         = 2.905E6
   REAL(rkind)    , PARAMETER :: XLS1         = 259.532

   REAL(rkind)    , PARAMETER :: XLS          = 2.85E6
   REAL(rkind)    , PARAMETER :: XLV          = 2.5E6
   REAL(rkind)    , PARAMETER :: XLF          = 3.50E5

   REAL(rkind)    , PARAMETER :: rhowater     = 1000.
   REAL(rkind)    , PARAMETER :: rhosnow      = 100.
   REAL(rkind)    , PARAMETER :: rhoair0      = 1.28
!
   REAL(rkind)    , PARAMETER :: n_ccn0       = 1.0E8
!
   REAL(rkind)    , PARAMETER :: DEGRAD       = 3.1415926/180.
   REAL(rkind)    , PARAMETER :: DPD          = 360./365.

   REAL(rkind)    , PARAMETER ::  SVP1=0.6112
   REAL(rkind)    , PARAMETER ::  SVP2=17.67
   REAL(rkind)    , PARAMETER ::  SVP3=29.65
   REAL(rkind)    , PARAMETER ::  SVPT0=273.15
   REAL(rkind)    , PARAMETER ::  EP_1=R_v/R_d-1.
   REAL(rkind)    , PARAMETER ::  EP_2=R_d/R_v
   REAL(rkind)    , PARAMETER ::  KARMAN=0.4
   REAL(rkind)    , PARAMETER ::  EOMEG=7.2921E-5
   REAL(rkind)    , PARAMETER ::  STBOLT=5.67051E-8

   REAL(rkind)    , PARAMETER ::  prandtl = 1./3.0
                                         ! constants for w-damping option
   REAL(rkind)    , PARAMETER ::  w_alpha = 0.3 ! strength m/s/s
   REAL(rkind)    , PARAMETER ::  w_beta  = 1.0 ! activation cfl number

       REAL(rkind) , PARAMETER ::  pq0=379.90516
       REAL(rkind) , PARAMETER ::  epsq2=0.2
       REAL(rkind) , PARAMETER ::  a2=17.2693882
       REAL(rkind) , PARAMETER ::  a3=273.16
       REAL(rkind) , PARAMETER ::  a4=35.86
       REAL(rkind) , PARAMETER ::  epsq=1.e-12
       REAL(rkind) , PARAMETER ::  p608=rvovrd-1.
!#if ( NMM_CORE == 1 )
       REAL(rkind) , PARAMETER ::  climit=1.e-20
       REAL(rkind) , PARAMETER ::  cm1=2937.4
       REAL(rkind) , PARAMETER ::  cm2=4.9283
       REAL(rkind) , PARAMETER ::  cm3=23.5518
!       REAL(rkind) , PARAMETER ::  defc=8.0
!       REAL(rkind) , PARAMETER ::  defm=32.0
       REAL(rkind) , PARAMETER ::  defc=0.0
       REAL(rkind) , PARAMETER ::  defm=99999.0
       REAL(rkind) , PARAMETER ::  epsfc=1./1.05
       REAL(rkind) , PARAMETER ::  epswet=0.0
       REAL(rkind) , PARAMETER ::  fcdif=1./3.
#ifdef HWRF
       REAL(rkind) , PARAMETER ::  fcm=0.0
#else
       REAL(rkind) , PARAMETER ::  fcm=0.00003
#endif
       REAL(rkind) , PARAMETER ::  gma=-r_d*(1.-rcp)*0.5
       REAL(rkind) , PARAMETER ::  p400=40000.0
       REAL(rkind) , PARAMETER ::  phitp=15000.0
       REAL(rkind) , PARAMETER ::  pi2=2.*3.1415926
       REAL(rkind) , PARAMETER ::  plbtm=105000.0
       REAL(rkind) , PARAMETER ::  plomd=64200.0
       REAL(rkind) , PARAMETER ::  pmdhi=35000.0
       REAL(rkind) , PARAMETER ::  q2ini=0.50
       REAL(rkind) , PARAMETER ::  rfcp=0.25/cp
       REAL(rkind) , PARAMETER ::  rhcrit_land=0.75
       REAL(rkind) , PARAMETER ::  rhcrit_sea=0.80
       REAL(rkind) , PARAMETER ::  rlag=14.8125
       REAL(rkind) , PARAMETER ::  rlx=0.90
       REAL(rkind) , PARAMETER ::  scq2=50.0
       REAL(rkind) , PARAMETER ::  slopht=0.001
       REAL(rkind) , PARAMETER ::  tlc=2.*0.703972477
       REAL(rkind) , PARAMETER ::  wa=0.15
       REAL(rkind) , PARAMETER ::  wght=0.35
       REAL(rkind) , PARAMETER ::  wpc=0.075
       REAL(rkind) , PARAMETER ::  z0land=0.10
#ifdef HWRF
       REAL(rkind) , PARAMETER ::  z0max=0.01
#else
       REAL(rkind) , PARAMETER ::  z0max=0.008
#endif
       REAL(rkind) , PARAMETER ::  z0sea=0.001
!#endif


   !  Earth

   !  The value for P2SI *must* be set to 1.0 for Earth
   !  Although, now we may not need this declaration here (see above)
   !REAL(rkind)    , PARAMETER :: P2SI         = 1.0

   !  Orbital constants:

   INTEGER , PARAMETER :: PLANET_YEAR = 365
   REAL(rkind) , PARAMETER :: OBLIQUITY = 23.5
   REAL(rkind) , PARAMETER :: ECCENTRICITY = 0.014
   REAL(rkind) , PARAMETER :: SEMIMAJORAXIS = 1.0 ! In AU
   ! Don't know the following values, so we'll fake them for now
   REAL(rkind) , PARAMETER :: zero_date = 0.0   ! Time of perihelion passage
   !  Fraction into the year (from perhelion) of the
   !  occurrence of the Northern Spring Equinox
   REAL(rkind) , PARAMETER :: EQUINOX_FRACTION= 0.0

 CONTAINS
   SUBROUTINE init_module_model_constants
   END SUBROUTINE init_module_model_constants
 END MODULE module_model_constants
