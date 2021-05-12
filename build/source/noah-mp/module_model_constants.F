!WRF:MODEL_LAYER:CONSTANTS
!

 MODULE module_model_constants
  USE nrtype

   !  2. Following are constants for use in defining real number bounds.

   !  A really small number.

   REAL(rk)    , PARAMETER :: epsilon         = 1.E-15

   !  4. Following is information related to the physical constants.

   !  These are the physical constants used within the model.

! JM NOTE -- can we name this grav instead?
   REAL(rk)    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)

#if ( NMM_CORE == 1 )
   REAL(rk)    , PARAMETER :: r_d          = 287.04
   REAL(rk)    , PARAMETER :: cp           = 1004.6
#else
   REAL(rk)    , PARAMETER :: r_d          = 287.
   REAL(rk)    , PARAMETER :: cp           = 7.*r_d/2.
#endif

   REAL(rk)    , PARAMETER :: r_v          = 461.6
   REAL(rk)    , PARAMETER :: cv           = cp-r_d
   REAL(rk)    , PARAMETER :: cpv          = 4.*r_v
   REAL(rk)    , PARAMETER :: cvv          = cpv-r_v
   REAL(rk)    , PARAMETER :: cvpm         = -cv/cp
   REAL(rk)    , PARAMETER :: cliq         = 4190.
   REAL(rk)    , PARAMETER :: cice         = 2106.
   REAL(rk)    , PARAMETER :: psat         = 610.78
   REAL(rk)    , PARAMETER :: rcv          = r_d/cv
   REAL(rk)    , PARAMETER :: rcp          = r_d/cp
   REAL(rk)    , PARAMETER :: rovg         = r_d/g
   REAL(rk)    , PARAMETER :: c2           = cp * rcv
   real    , parameter :: mwdry        = 28.966 ! molecular weight of dry air (g/mole)

   REAL(rk)    , PARAMETER :: p1000mb      = 100000.
   REAL(rk)    , PARAMETER :: t0           = 300.
   REAL(rk)    , PARAMETER :: p0           = p1000mb
   REAL(rk)    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL(rk)    , PARAMETER :: cvovcp       = 1./cpovcv
   REAL(rk)    , PARAMETER :: rvovrd       = r_v/r_d

   REAL(rk)    , PARAMETER :: reradius     = 1./6370.0e03

   REAL(rk)    , PARAMETER :: asselin      = .025
!   REAL(rk)    , PARAMETER :: asselin      = .0
   REAL(rk)    , PARAMETER :: cb           = 25.

   REAL(rk)    , PARAMETER :: XLV0         = 3.15E6
   REAL(rk)    , PARAMETER :: XLV1         = 2370.
   REAL(rk)    , PARAMETER :: XLS0         = 2.905E6
   REAL(rk)    , PARAMETER :: XLS1         = 259.532

   REAL(rk)    , PARAMETER :: XLS          = 2.85E6
   REAL(rk)    , PARAMETER :: XLV          = 2.5E6
   REAL(rk)    , PARAMETER :: XLF          = 3.50E5

   REAL(rk)    , PARAMETER :: rhowater     = 1000.
   REAL(rk)    , PARAMETER :: rhosnow      = 100.
   REAL(rk)    , PARAMETER :: rhoair0      = 1.28
!
   REAL(rk)    , PARAMETER :: n_ccn0       = 1.0E8
!
   REAL(rk)    , PARAMETER :: DEGRAD       = 3.1415926/180.
   REAL(rk)    , PARAMETER :: DPD          = 360./365.

   REAL(rk)    , PARAMETER ::  SVP1=0.6112
   REAL(rk)    , PARAMETER ::  SVP2=17.67
   REAL(rk)    , PARAMETER ::  SVP3=29.65
   REAL(rk)    , PARAMETER ::  SVPT0=273.15
   REAL(rk)    , PARAMETER ::  EP_1=R_v/R_d-1.
   REAL(rk)    , PARAMETER ::  EP_2=R_d/R_v
   REAL(rk)    , PARAMETER ::  KARMAN=0.4
   REAL(rk)    , PARAMETER ::  EOMEG=7.2921E-5
   REAL(rk)    , PARAMETER ::  STBOLT=5.67051E-8

   REAL(rk)    , PARAMETER ::  prandtl = 1./3.0
                                         ! constants for w-damping option
   REAL(rk)    , PARAMETER ::  w_alpha = 0.3 ! strength m/s/s
   REAL(rk)    , PARAMETER ::  w_beta  = 1.0 ! activation cfl number

       REAL(rk) , PARAMETER ::  pq0=379.90516
       REAL(rk) , PARAMETER ::  epsq2=0.2
       REAL(rk) , PARAMETER ::  a2=17.2693882
       REAL(rk) , PARAMETER ::  a3=273.16
       REAL(rk) , PARAMETER ::  a4=35.86
       REAL(rk) , PARAMETER ::  epsq=1.e-12
       REAL(rk) , PARAMETER ::  p608=rvovrd-1.
!#if ( NMM_CORE == 1 )
       REAL(rk) , PARAMETER ::  climit=1.e-20
       REAL(rk) , PARAMETER ::  cm1=2937.4
       REAL(rk) , PARAMETER ::  cm2=4.9283
       REAL(rk) , PARAMETER ::  cm3=23.5518
!       REAL(rk) , PARAMETER ::  defc=8.0
!       REAL(rk) , PARAMETER ::  defm=32.0
       REAL(rk) , PARAMETER ::  defc=0.0
       REAL(rk) , PARAMETER ::  defm=99999.0
       REAL(rk) , PARAMETER ::  epsfc=1./1.05
       REAL(rk) , PARAMETER ::  epswet=0.0
       REAL(rk) , PARAMETER ::  fcdif=1./3.
#ifdef HWRF
       REAL(rk) , PARAMETER ::  fcm=0.0
#else
       REAL(rk) , PARAMETER ::  fcm=0.00003
#endif
       REAL(rk) , PARAMETER ::  gma=-r_d*(1.-rcp)*0.5
       REAL(rk) , PARAMETER ::  p400=40000.0
       REAL(rk) , PARAMETER ::  phitp=15000.0
       REAL(rk) , PARAMETER ::  pi2=2.*3.1415926
       REAL(rk) , PARAMETER ::  plbtm=105000.0
       REAL(rk) , PARAMETER ::  plomd=64200.0
       REAL(rk) , PARAMETER ::  pmdhi=35000.0
       REAL(rk) , PARAMETER ::  q2ini=0.50
       REAL(rk) , PARAMETER ::  rfcp=0.25/cp
       REAL(rk) , PARAMETER ::  rhcrit_land=0.75
       REAL(rk) , PARAMETER ::  rhcrit_sea=0.80
       REAL(rk) , PARAMETER ::  rlag=14.8125
       REAL(rk) , PARAMETER ::  rlx=0.90
       REAL(rk) , PARAMETER ::  scq2=50.0
       REAL(rk) , PARAMETER ::  slopht=0.001
       REAL(rk) , PARAMETER ::  tlc=2.*0.703972477
       REAL(rk) , PARAMETER ::  wa=0.15
       REAL(rk) , PARAMETER ::  wght=0.35
       REAL(rk) , PARAMETER ::  wpc=0.075
       REAL(rk) , PARAMETER ::  z0land=0.10
#ifdef HWRF
       REAL(rk) , PARAMETER ::  z0max=0.01
#else
       REAL(rk) , PARAMETER ::  z0max=0.008
#endif
       REAL(rk) , PARAMETER ::  z0sea=0.001
!#endif


   !  Earth

   !  The value for P2SI *must* be set to 1.0 for Earth
   !  Although, now we may not need this declaration here (see above)
   !REAL(rk)    , PARAMETER :: P2SI         = 1.0

   !  Orbital constants:

   INTEGER , PARAMETER :: PLANET_YEAR = 365
   REAL(rk) , PARAMETER :: OBLIQUITY = 23.5
   REAL(rk) , PARAMETER :: ECCENTRICITY = 0.014
   REAL(rk) , PARAMETER :: SEMIMAJORAXIS = 1.0 ! In AU
   ! Don't know the following values, so we'll fake them for now
   REAL(rk) , PARAMETER :: zero_date = 0.0   ! Time of perihelion passage
   !  Fraction into the year (from perhelion) of the
   !  occurrence of the Northern Spring Equinox
   REAL(rk) , PARAMETER :: EQUINOX_FRACTION= 0.0

 CONTAINS
   SUBROUTINE init_module_model_constants
   END SUBROUTINE init_module_model_constants
 END MODULE module_model_constants
