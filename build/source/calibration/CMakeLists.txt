# Generic calibration routines
set(CALIBR_GENERIC
    ${CALIBR_DIR}/parameter_search.f90
    ${CALIBR_DIR}/parameter_fields.f90
    ${CALIBR_DIR}/parameter_constraints.f90
    ${CALIBR_DIR}/model_calibration.f90
    ${CALIBR_DIR}/nsga2.f90
)

# TOML-dependent calibration configuration
if(USE_TOML)
    list(APPEND CALIBR_GENERIC
        ${CALIBR_DIR}/calibration_config.f90
    )
endif()

# SUMMA-specific calibration routines
set(CALIBR_SUMMA
    ${CALIBR_DIR}/summa_parameter_fields.f90
    ${CALIBR_DIR}/summa_parameter_spec.f90
    ${CALIBR_DIR}/summa_parameter_export.f90
    ${CALIBR_DIR}/calibration_output.f90
)

# Existing SUMMA calibration sources
set(CALIBR
    ${CALIBR_GENERIC}
    ${CALIBR_SUMMA}
)

# MPI-dependent calibration
set(CALIBR_MPI
    ${CALIBR_DIR}/calibration_parallel.f90
    ${CALIBR_DIR}/dynamic_case_dispatch.f90
    ${CALIBR_DIR}/summa_parameter_sampling.f90
)
