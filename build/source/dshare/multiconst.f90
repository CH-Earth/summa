MODULE multiconst
 USE nrtype
 ! define physical constants
 REAL(DP), PARAMETER         :: ave_slp      =  101325.0_dp      ! mean sea level pressure              (Pa)
 REAL(DP), PARAMETER         :: vkc          =       0.4_dp      ! von Karman constant                  (-)
 REAL(DP), PARAMETER         :: satvpfrz     =     610.8_dp      ! sat vapour pressure at 273.16K       (Pa)
 REAL(DP), PARAMETER         :: w_ratio      =       0.622_dp    ! molecular ratio water to dry air     (-)
 REAL(DP), PARAMETER         :: R_da         =     287.053_dp    ! gas constant for dry air             (Pa K-1 m3 kg-1; J kg-1 K-1)
 REAL(DP), PARAMETER         :: R_wv         =     461.285_dp    ! gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
 REAL(DP), PARAMETER         :: gravity      =       9.80616_dp  ! acceleration of gravity              (m s-2)
 REAL(DP), PARAMETER         :: Cp_air       =    1005._dp       ! specific heat of air                 (J kg-1 K-1)
 REAL(DP), PARAMETER         :: Cp_ice       =    2114._dp       ! specific heat of ice                 (J kg-1 K-1)
 REAL(DP), PARAMETER         :: Cp_soil      =     850._dp       ! specific heat of soil                (J kg-1 K-1)
 REAL(DP), PARAMETER         :: Cp_water     =    4181._dp       ! specific heat of liquid water        (J kg-1 K-1)
 REAL(DP), PARAMETER         :: Tfreeze      =     273.16_dp     ! temperature at freezing              (K)
 REAL(DP), PARAMETER         :: TriplPt      =     273.16_dp     ! triple point of water                (K)
 REAL(DP), PARAMETER         :: LH_fus       =  333700.0_dp      ! latent heat of fusion                (J kg-1)
 REAL(DP), PARAMETER         :: LH_vap       = 2501000.0_dp      ! latent heat of vaporization          (J kg-1)
 REAL(DP), PARAMETER         :: LH_sub       = 2834700.0_dp      ! latent heat of sublimation           (J kg-1)
 REAL(DP), PARAMETER         :: sigma        =       5.6705d-8   ! Stefan Boltzman constant             (W m-2 K-4)
 REAL(DP), PARAMETER         :: em_sno       =       0.99_dp     ! emissivity of snow                   (-)
 REAL(DP), PARAMETER         :: lambda_air   =       0.023_dp    ! thermal conductivity of air          (W m-1 K-1)
 REAL(DP), PARAMETER         :: lambda_ice   =       2.29_dp     ! thermal conductivity of ice          (W m-1 K-1)
 REAL(DP), PARAMETER         :: lambda_soil  =       3.21_dp     ! thermal conductivity of soil         (W m-1 K-1)
 REAL(DP), PARAMETER         :: lambda_water =       0.60_dp     ! thermal conductivity of liquid water (W m-1 K-1)
 REAL(DP), PARAMETER         :: iden_air     =       1.293_dp    ! intrinsic density of air             (kg m-3)
 REAL(DP), PARAMETER         :: iden_ice     =     917.0_dp      ! intrinsic density of ice             (kg m-3)
 REAL(DP), PARAMETER         :: iden_water   =    1000.0_dp      ! intrinsic density of liquid water    (kg m-3)
 REAL(DP), PARAMETER         :: secprday     =   86400._dp       ! number of seconds in a day
 REAL(DP), PARAMETER         :: secprhour    =    3600._dp       ! number of seconds in an hour
 REAL(DP), PARAMETER         :: secprmin     =      60._dp       ! number of seconds in a minute
END MODULE multiconst
