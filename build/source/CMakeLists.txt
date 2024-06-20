#===============================================================================
# List of Source Files for SUMMA
# All Variables are reachable from the parent scope
#===============================================================================

# NOAHMP modules
set(NOAHMP
    ${NOAHMP_DIR}/module_model_constants.F
    ${NOAHMP_DIR}/module_sf_noahutl.F
    ${NOAHMP_DIR}/module_sf_noahlsm.F
    ${NOAHMP_DIR}/module_sf_noahmplsm.F
    CACHE INTERNAL "NOAHMP")

# Free versions of numerical recipes utilities for NOAH-MP modules
set(NRUTIL
    ${ENGINE_DIR}/f2008funcs.f90
    ${ENGINE_DIR}/nr_utility.f90
    ${ENGINE_DIR}/nrtype.f90
    CACHE INTERNAL "NRUTIL")

# Free versions of numerical recipes procedures for SUMMA modules
set(NRPROC
    ${ENGINE_DIR}/expIntegral.f90
    ${ENGINE_DIR}/spline_int.f90
    ${ENGINE_DIR}/hyp_2F1.f90
    CACHE INTERNAL "NRPROC")

# Hook-up modules
set(HOOKUP
    ${HOOKUP_DIR}/ascii_util.f90
    ${HOOKUP_DIR}/summaFileManager.f90
    CACHE INTERNAL "HOOKUP")

# Data modules
set(DATAMS
    ${DSHARE_DIR}/data_types.f90
    ${DSHARE_DIR}/flxMapping.f90
    ${DSHARE_DIR}/get_ixname.f90
    ${DSHARE_DIR}/globalData.f90
    ${DSHARE_DIR}/multiconst.f90
    ${DSHARE_DIR}/outpt_stat.f90
    ${DSHARE_DIR}/popMetadat.f90
    ${DSHARE_DIR}/var_lookup.f90
    CACHE INTERNAL "DATAMS")
set(DATAMS_SUNDIALS
    ${DSHARE_DIR}/type4ida.f90
    ${DSHARE_DIR}/type4kinsol.f90
    CACHE INTERNAL "DATAMS_SUNDIALS")
set(DATAMS_ACTORS
    ${ACTORS_DIR}/global/actor_data_types.f90
    CACHE INTERNAL "DATAMS_ACTORS")

# Utility modules
set(UTILMS
    ${ENGINE_DIR}/matrixOper.f90
    ${ENGINE_DIR}/mDecisions.f90
    ${ENGINE_DIR}/snow_utils.f90
    ${ENGINE_DIR}/soil_utils.f90
    ${ENGINE_DIR}/time_utils.f90
    ${ENGINE_DIR}/updatState.f90
    CACHE INTERNAL "UTILMS")
set(UTILMS_SUNDIALS
    ${ENGINE_DIR}/soil_utilsAddPrime.f90
    ${ENGINE_DIR}/updatStateWithPrime.f90
    CACHE INTERNAL "UTILMS_SUNDIALS")

# NetCDF routines
set(NETCDF
    ${NETCDF_DIR}/def_output.f90
    ${NETCDF_DIR}/modelwrite.f90
    ${NETCDF_DIR}/netcdf_util.f90
    ${NETCDF_DIR}/read_icond.f90
    CACHE INTERNAL "NETCDF")

# Preliminary modules
set(PRELIM
    ${ENGINE_DIR}/allocspace.f90
    ${ENGINE_DIR}/check_icond.f90
    ${ENGINE_DIR}/checkStruc.f90
    ${ENGINE_DIR}/childStruc.f90
    ${ENGINE_DIR}/conv_funcs.f90
    ${ENGINE_DIR}/ffile_info.f90
    ${ENGINE_DIR}/read_pinit.f90
    ${ENGINE_DIR}/read_attrb.f90
    ${ENGINE_DIR}/paramCheck.f90
    ${ENGINE_DIR}/pOverwrite.f90
    ${ENGINE_DIR}/sunGeomtry.f90
    ${ENGINE_DIR}/read_param.f90
    CACHE INTERNAL "PRELIM")

# Model run support modules
set(MODRUN
    ${ENGINE_DIR}/canopySnow.f90
    ${ENGINE_DIR}/derivforce.f90
    ${ENGINE_DIR}/enthalpyTemp.f90
    ${ENGINE_DIR}/getVectorz.f90
    ${ENGINE_DIR}/indexState.f90
    ${ENGINE_DIR}/layerMerge.f90
    ${ENGINE_DIR}/layerDivide.f90
    ${ENGINE_DIR}/qTimeDelay.f90
    ${ENGINE_DIR}/snowAlbedo.f90
    ${ENGINE_DIR}/snwCompact.f90
    ${ENGINE_DIR}/tempAdjust.f90
    ${ENGINE_DIR}/updateVars.f90
    ${ENGINE_DIR}/var_derive.f90
    ${ENGINE_DIR}/volicePack.f90
    CACHE INTERNAL "MODRUN")
set(MODRUN_NOT_ACTORS
    ${ENGINE_DIR}/read_force.f90
    CACHE INTERNAL "MODRUN_NOT_ACTORS")
set(MODRUN_SUNDIALS
    ${ENGINE_DIR}/tol4ida.f90
    ${ENGINE_DIR}/updateVarsWithPrime.f90
    CACHE INTERNAL "MODRUN_SUNDIALS")

# Solver main modules
set(SOLVER
    ${ENGINE_DIR}/bigAquifer.f90
    ${ENGINE_DIR}/computFlux.f90
    ${ENGINE_DIR}/computHeatCap.f90
    ${ENGINE_DIR}/computJacob.f90
    ${ENGINE_DIR}/computResid.f90
    ${ENGINE_DIR}/computSnowDepth.f90
    ${ENGINE_DIR}/computThermConduct.f90
    ${ENGINE_DIR}/coupled_em.f90
    ${ENGINE_DIR}/diagn_evar.f90
    ${ENGINE_DIR}/eval8summa.f90
    ${ENGINE_DIR}/groundwatr.f90
    ${ENGINE_DIR}/opSplittin.f90
    ${ENGINE_DIR}/snowLiqFlx.f90
    ${ENGINE_DIR}/soilLiqFlx.f90
    ${ENGINE_DIR}/ssdNrgFlux.f90
    ${ENGINE_DIR}/stomResist.f90
    ${ENGINE_DIR}/summaSolve4homegrown.f90
    ${ENGINE_DIR}/systemSolv.f90
    ${ENGINE_DIR}/varSubstep.f90
    ${ENGINE_DIR}/vegLiqFlux.f90
    ${ENGINE_DIR}/vegNrgFlux.f90
    ${ENGINE_DIR}/vegPhenlgy.f90
    ${ENGINE_DIR}/vegSWavRad.f90
    CACHE INTERNAL "SOLVER")
set(SOLVER_NOT_ACTORS
    ${ENGINE_DIR}/run_oneGRU.f90
    ${ENGINE_DIR}/run_oneHRU.f90
    CACHE INTERNAL "SOLVER_NOT_ACTORS")
set(SOLVER_SUNDIALS
    ${ENGINE_DIR}/computJacobWithPrime.f90
    ${ENGINE_DIR}/computResidWithPrime.f90
    ${ENGINE_DIR}/eval8summaWithPrime.f90
    ${ENGINE_DIR}/summaSolve4ida.f90
    ${ENGINE_DIR}/summaSolve4kinsol.f90
    CACHE INTERNAL "SOLVER_SUNDIALS")

# Driver support modules
set(DRIVER
    ${DRIVER_DIR}/summa_type.f90
    ${DRIVER_DIR}/summa_setup.f90
    ${DRIVER_DIR}/summa_restart.f90
    ${DRIVER_DIR}/summa_alarms.f90
    ${DRIVER_DIR}/summa_globalData.f90
    CACHE INTERNAL "DRIVER")
set(DRIVER_NOT_ACTORS
    ${DRIVER_DIR}/summa_util.f90
    ${DRIVER_DIR}/summa_defineOutput.f90
    ${DRIVER_DIR}/summa_init.f90
    ${DRIVER_DIR}/summa_forcing.f90
    ${DRIVER_DIR}/summa_modelRun.f90
    ${DRIVER_DIR}/summa_writeOutput.f90
    CACHE INTERNAL "DRIVER_NOT_ACTORS")
set(DRIVER_NEXTGEN
    ${DRIVER_DIR}/summa_bmi.f90
    CACHE INTERNAL "DRIVER_NEXTGEN")
