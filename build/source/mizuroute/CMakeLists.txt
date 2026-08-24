set(MIZU_TYPE
    ${SUMMA_MIZU_DIR}/nrtype_shim.f90
    ${MIZU_DIR}/public_var.f90
    ${MIZU_DIR}/datetime_data.f90
    ${MIZU_DIR}/dataTypes.f90
    ${SUMMA_MIZU_DIR}/mizuRoute_types.f90
)

set(MIZU_BASE
    ${MIZU_DIR}/var_lookup.f90
    ${MIZU_DIR}/nr_utils.f90
    ${MIZU_DIR}/gamma_func.f90
    ${MIZU_DIR}/ncio_utils.f90
    ${MIZU_DIR}/ascii_utils.f90
    ${MIZU_DIR}/base_route.f90
)

set(MIZU_DATA
    ${SUMMA_MIZU_DIR}/mizuRoute_types.f90
    ${SUMMA_MIZU_DIR}/water_balance_shim.f90
    ${SUMMA_MIZU_DIR}/globalData_shim.f90
    ${SUMMA_MIZU_DIR}/pio_utils_shim.f90
)

set(MIZU_TOPO
    ${MIZU_DIR}/allocation.f90
    ${MIZU_DIR}/popMetadat.f90
    ${MIZU_DIR}/hydraulic.f90
    ${MIZU_DIR}/network_topo.f90
    ${MIZU_DIR}/read_param.f90
    ${MIZU_DIR}/process_param.f90
    ${MIZU_DIR}/process_ntopo.f90
    ${MIZU_DIR}/read_streamSeg.f90
    ${MIZU_DIR}/write_streamSeg.f90
)

set(MIZU_INIT
    ${SUMMA_MIZU_DIR}/init_model_data_shim.f90
)

set(MIZU_RMAP
    ${MIZU_DIR}/standalone/read_remap.f90
    ${MIZU_DIR}/process_remap.f90
)

set(MIZU_ROUTE
    ${MIZU_DIR}/advection_diffusion.f90
    ${MIZU_DIR}/data_assimilation.f90
    ${MIZU_DIR}/accum_runoff.f90
    ${MIZU_DIR}/irf_route.f90
    ${MIZU_DIR}/kwt_route.f90
    ${MIZU_DIR}/kwe_route.f90
    ${MIZU_DIR}/mc_route.f90
    ${MIZU_DIR}/dfw_route.f90
)

set(MIZU_COUPLING
    ${SUMMA_MIZU_DIR}/init_mizuRoute.f90
    ${SUMMA_MIZU_DIR}/network_routing.f90
)
