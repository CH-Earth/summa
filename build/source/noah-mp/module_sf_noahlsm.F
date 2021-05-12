MODULE module_sf_noahlsm
  USE nrtype
  USE module_model_constants

!   REAL(SUMMA_PREC), PARAMETER    :: CP = 1004.5
  REAL(SUMMA_PREC), PARAMETER      :: RD = 287.04, SIGMA = 5.67E-8,                 &
                          CPH2O = 4.218E+3,CPICE = 2.106E+3,            &
                          LSUBF = 3.335E+5,                             &
                          EMISSI_S = 0.95

! VEGETATION PARAMETERS
        INTEGER :: LUCATS , BARE
        INTEGER :: NATURAL
        integer, PARAMETER :: NLUS=50
        CHARACTER(LEN=256) LUTYPE
        INTEGER, DIMENSION(1:NLUS) :: NROTBL
        real, dimension(1:NLUS) ::  SNUPTBL, RSTBL, RGLTBL, HSTBL,                &
                                    SHDTBL, MAXALB,                               &
                                    EMISSMINTBL, EMISSMAXTBL,                     &
                                    LAIMINTBL, LAIMAXTBL,                         &
                                    Z0MINTBL, Z0MAXTBL,                           &
                                    ALBEDOMINTBL, ALBEDOMAXTBL
        REAL(SUMMA_PREC) ::   TOPT_DATA,CMCMAX_DATA,CFACTR_DATA,RSMAX_DATA

! SOIL PARAMETERS
        INTEGER :: SLCATS
        INTEGER, PARAMETER :: NSLTYPE=30
        CHARACTER(LEN=256) SLTYPE
        REAL(SUMMA_PREC), DIMENSION (1:NSLTYPE) :: BB,DRYSMC,F11,                           &
        MAXSMC, REFSMC,SATPSI,SATDK,SATDW, WLTSMC,QTZ
        ! MPC add van Genutchen parameters
        REAL(SUMMA_PREC), DIMENSION (1:NSLTYPE) :: theta_res, theta_sat,                    &
        vGn_alpha, vGn_n, k_soil

! LSM GENERAL PARAMETERS
        INTEGER :: SLPCATS
        INTEGER, PARAMETER :: NSLOPE=30
        REAL(SUMMA_PREC), DIMENSION (1:NSLOPE) :: SLOPE_DATA
        REAL(SUMMA_PREC) ::  SBETA_DATA,FXEXP_DATA,CSOIL_DATA,SALP_DATA,REFDK_DATA,           &
                 REFKDT_DATA,FRZK_DATA,ZBOT_DATA,  SMLOW_DATA,SMHIGH_DATA,        &
                        CZIL_DATA
        REAL(SUMMA_PREC) ::  LVCOEF_DATA

        CHARACTER*256  :: err_message
        integer, private :: iloc, jloc

END MODULE module_sf_noahlsm
