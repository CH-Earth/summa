function(compile_with_ida PARENT_DIR, DIR_SUNDIALS)
    find_package(LAPACK REQUIRED)
    set(EXEC_NAME summa_sundials)

    link_directories(${DIR_SUNDIALS}/lib64)
    set(CMAKE_BUILD_RPATH "${DIR_SUNDIALS}/lib64")
    set(SUMMA_INCLUDES 
        "$ENV{EBROOTNETCDFMINFORTRAN}/include"
        "${DIR_SUNDIALS}/include"
        "${DIR_SUNDIALS}/fortran"
        ${netCDF_INCLUDES}
        ${LAPACK_INCLUDES})
    
    set(SUMMA_LIBS
        -lsundials_fnvecmanyvector_mod 
        -lsundials_fida_mod 
        -lsundials_fnvecserial_mod 
        -lsundials_fsunlinsoldense_mod 
        -lsundials_fsunmatrixdense_mod 
        -lnetcdff
        -lopenblas
        ${netCDF_LIBRARIES}
        ${LAPACK_LIBRARIES}
        SUMMA_NOAHMP)
    
    set(SUMMA_ACTORS_INCLUDES 
        ${CAF_INCLUDES}
        "$ENV{EBROOTNETCDFMINFORTRAN}/include"
        ${LAPACK_INCLUDES}
        "${DIR_SUNDIALS}/include"
        "${PARENT_DIR}/build/includes/global"
        "${PARENT_DIR}/build/includes/summa_actor"
        "${PARENT_DIR}/build/includes/gru_actor"
        "${PARENT_DIR}/build/includes/job_actor"
        "${PARENT_DIR}/build/includes/file_access_actor"
        "${PARENT_DIR}/build/includes/hru_actor")

    set(SUMMA_ACTORS_LIBS   
        ${CAF_LIBRARIES}
        ${netCDF_LIBRARIES}
        ${LAPACK_LIBRARIES}
        -lopenblas
        -lcaf_core
        -lcaf_io
        summa
        -lnetcdff
        -lsundials_fnvecmanyvector_mod 
        -lsundials_fida_mod 
        -lsundials_fnvecserial_mod 
        -lsundials_fsunlinsoldense_mod 
        -lsundials_fsunmatrixdense_mod)


    set(ACTORS_DIR ${PARENT_DIR}/build/source/actors)
    set(DRIVER_DIR ${PARENT_DIR}/build/source/driver)
    set(DSHARE_DIR ${PARENT_DIR}/build/source/dshare)
    set(ENGINE_DIR ${PARENT_DIR}/build/source/engine)
    set(HOOKUP_DIR ${PARENT_DIR}/build/source/hookup)
    set(NETCDF_DIR ${PARENT_DIR}/build/source/netcdf)
    set(NOAHMP_DIR ${PARENT_DIR}/build/source/noah-mp)
    set(FILE_ACCESS_DIR ${ACTORS_DIR}/file_access_actor)
    set(JOB_ACTOR_DIR ${ACTORS_DIR}/job_actor)
    set(HRU_ACTOR_DIR ${ACTORS_DIR}/hru_actor)
    set(GRU_ACTOR_DIR ${ACTORS_DIR}/gru_actor)
    set(SUMMA_DSHARE_DIR ${PARENT_DIR}/build/summa-sundials/build/source/dshare)
    set(SUMMA_ENGINE_DIR ${PARENT_DIR}/build/summa-sundials/build/source/engine)
    set(SUMMA_NOAHMP_DIR ${PARENT_DIR}/build/summa-sundials/build/source/noah-mp)

    set(NRUTIL
        ${SUMMA_ENGINE_DIR}/nrtype.f90
        ${SUMMA_ENGINE_DIR}/f2008funcs.f90
        ${SUMMA_ENGINE_DIR}/nr_utility.f90)
    
    set(NRPROC
        ${SUMMA_ENGINE_DIR}/expIntegral.f90
        ${SUMMA_ENGINE_DIR}/spline_int.f90)
    
    SET(HOOKUP
        ${HOOKUP_DIR}/ascii_util.f90
        ${HOOKUP_DIR}/summaActors_FileManager.f90)
    
    SET(DATAMS
        ${SUMMA_DSHARE_DIR}/multiconst.f90
        ${SUMMA_DSHARE_DIR}/var_lookup.f90
        ${DSHARE_DIR}/data_types.f90
        ${DSHARE_DIR}/globalData.f90
        ${SUMMA_DSHARE_DIR}/flxMapping.f90)
    
    SET(DEPENDS_ON_FILEMANAGER
        ${SUMMA_DSHARE_DIR}/get_ixname.f90
        ${SUMMA_DSHARE_DIR}/popMetadat.f90
        ${SUMMA_DSHARE_DIR}/outpt_stat.f90)
    
    SET(UTILMS
        ${SUMMA_ENGINE_DIR}/time_utils.f90
        ${ENGINE_DIR}/sundials/mDecisions.f90
        ${SUMMA_ENGINE_DIR}/snow_utils.f90
        ${SUMMA_ENGINE_DIR}/soil_utils.f90
        ${SUMMA_ENGINE_DIR}/soil_utilsAddSundials.f90
        ${SUMMA_ENGINE_DIR}/updatState.f90
        ${SUMMA_ENGINE_DIR}/updatStateSundials.f90
        ${SUMMA_ENGINE_DIR}/matrixOper.f90)
    
    set(SOLVER
        ${ENGINE_DIR}/vegPhenlgy.f90
        ${SUMMA_ENGINE_DIR}/diagn_evar.f90
        ${SUMMA_ENGINE_DIR}/stomResist.f90
        ${SUMMA_ENGINE_DIR}/groundwatr.f90
        ${SUMMA_ENGINE_DIR}/vegSWavRad.f90
        ${SUMMA_ENGINE_DIR}/vegNrgFlux.f90
        ${SUMMA_ENGINE_DIR}/ssdNrgFlux.f90
        ${SUMMA_ENGINE_DIR}/vegLiqFlux.f90
        ${SUMMA_ENGINE_DIR}/snowLiqFlx.f90
        ${SUMMA_ENGINE_DIR}/soilLiqFlx.f90
        ${SUMMA_ENGINE_DIR}/bigAquifer.f90
        ${SUMMA_ENGINE_DIR}/computFlux.f90
        ${SUMMA_ENGINE_DIR}/type4IDA.f90
        ${SUMMA_ENGINE_DIR}/tol4IDA.f90 
        ${SUMMA_ENGINE_DIR}/computEnthalpy.f90
        ${SUMMA_ENGINE_DIR}/computHeatCap.f90
        ${SUMMA_ENGINE_DIR}/computThermConduct.f90
        ${SUMMA_ENGINE_DIR}/computResid.f90
        ${SUMMA_ENGINE_DIR}/computJacob.f90
        ${SUMMA_ENGINE_DIR}/eval8summa.f90
        ${SUMMA_ENGINE_DIR}/summaSolve.f90
        ${SUMMA_ENGINE_DIR}/systemSolv.f90
        ${SUMMA_ENGINE_DIR}/computResidSundials.f90
        ${SUMMA_ENGINE_DIR}/eval8summaSundials.f90
        ${SUMMA_ENGINE_DIR}/computJacobSundials.f90
        ${SUMMA_ENGINE_DIR}/computSnowDepth.f90
        ${SUMMA_ENGINE_DIR}/summaSolveSundialsIDA.f90
        ${SUMMA_ENGINE_DIR}/systemSolvSundials.f90
        ${SUMMA_ENGINE_DIR}/varSubstep.f90
        ${SUMMA_ENGINE_DIR}/varSubstepSundials.f90
        ${SUMMA_ENGINE_DIR}/opSplittin.f90
        ${ENGINE_DIR}/sundials/coupled_em.f90)
    
    set(INTERFACE
        ${ACTORS_DIR}/global/cppwrap_datatypes.f90
        ${ACTORS_DIR}/global/cppwrap_auxiliary.f90
        ${ACTORS_DIR}/global/cppwrap_metadata.f90)
    
    set(FILE_ACCESS_INTERFACE
        ${FILE_ACCESS_DIR}/fortran_code/output_structure.f90
        ${FILE_ACCESS_DIR}/fortran_code/cppwrap_fileAccess.f90
        ${FILE_ACCESS_DIR}/fortran_code/read_attribute.f90
        ${FILE_ACCESS_DIR}/fortran_code/read_forcing.f90
        ${FILE_ACCESS_DIR}/fortran_code/read_param.f90
        ${FILE_ACCESS_DIR}/fortran_code/read_initcond.f90
        ${FILE_ACCESS_DIR}/fortran_code/writeOutputFromOutputStructure.f90
        ${FILE_ACCESS_DIR}/fortran_code/write_to_netcdf.f90)
    
    set(JOB_INTERFACE
        ${JOB_ACTOR_DIR}/job_actor.f90)

    set(HRU_INTERFACE
        ${HRU_ACTOR_DIR}/fortran_code/model_run.f90
        ${HRU_ACTOR_DIR}/fortran_code/setup_hru.f90
        ${HRU_ACTOR_DIR}/fortran_code/restart.f90
        ${HRU_ACTOR_DIR}/fortran_code/hru_actor.f90
        ${HRU_ACTOR_DIR}/fortran_code/init_hru_actor.f90
        ${HRU_ACTOR_DIR}/fortran_code/outputStrucWrite.f90
        ${HRU_ACTOR_DIR}/fortran_code/hru_writeOutput.f90)

    set(PRELIM
        ${SUMMA_ENGINE_DIR}/conv_funcs.f90
        ${SUMMA_ENGINE_DIR}/sunGeomtry.f90
        ${SUMMA_ENGINE_DIR}/convE2Temp.f90
        ${ENGINE_DIR}/allocspaceActors.f90
        ${ENGINE_DIR}/alloc_fileAccess.f90
        ${SUMMA_ENGINE_DIR}/checkStruc.f90
        ${SUMMA_ENGINE_DIR}/childStruc.f90
        ${ENGINE_DIR}/ffile_info.f90
        ${ENGINE_DIR}/read_dimension.f90
        ${ENGINE_DIR}/read_pinit.f90
        ${SUMMA_ENGINE_DIR}/pOverwrite.f90
        ${SUMMA_ENGINE_DIR}/paramCheck.f90
        ${ENGINE_DIR}/check_icondActors.f90)
    
    set(NOAHMP
        ${SUMMA_NOAHMP_DIR}/module_model_constants.F
        ${SUMMA_NOAHMP_DIR}/module_sf_noahutl.F
        ${SUMMA_NOAHMP_DIR}/module_sf_noahlsm.F
        ${SUMMA_NOAHMP_DIR}/module_sf_noahmplsm.F)
    
    set(MODRUN
        ${SUMMA_ENGINE_DIR}/indexState.f90
        ${SUMMA_ENGINE_DIR}/getVectorz.f90
        ${SUMMA_ENGINE_DIR}/t2enthalpy.f90
        ${SUMMA_ENGINE_DIR}/updateVars.f90
        ${SUMMA_ENGINE_DIR}/updateVarsSundials.f90
        ${SUMMA_ENGINE_DIR}/var_derive.f90
        ${ENGINE_DIR}/derivforce.f90
        ${SUMMA_ENGINE_DIR}/snowAlbedo.f90
        ${SUMMA_ENGINE_DIR}/canopySnow.f90
        ${SUMMA_ENGINE_DIR}/tempAdjust.f90
        ${SUMMA_ENGINE_DIR}/snwCompact.f90
        ${SUMMA_ENGINE_DIR}/layerMerge.f90
        ${SUMMA_ENGINE_DIR}/layerDivide.f90
        ${SUMMA_ENGINE_DIR}/volicePack.f90
        ${SUMMA_ENGINE_DIR}/qTimeDelay.f90)

    set(NETCDF
        ${NETCDF_DIR}/netcdf_util.f90
        ${NETCDF_DIR}/def_output.f90
        ${NETCDF_DIR}/writeOutput.f90
        ${NETCDF_DIR}/read_icondActors.f90)

    set(DRIVER
        ${DRIVER_DIR}/summaActors_type.f90
        ${DRIVER_DIR}/summaActors_util.f90
        ${DRIVER_DIR}/summaActors_globalData.f90
        ${DRIVER_DIR}/summaActors_alarms.f90)

    set(COMM_ALL
        ${NRPROC}
        ${DATAMS}
        ${INTERFACE}
        ${HOOKUP}
        ${DEPENDS_ON_FILEMANAGER}
        ${UTILMS})

    set(SUMMA_ALL
        ${NETCDF}
        ${PRELIM}
        ${MODRUN}
        ${SOLVER}
        ${DRIVER}
        ${JOB_INTERFACE}
        ${FILE_ACCESS_INTERFACE}
        ${HRU_INTERFACE}
        ${GRU_INTERFACE})

    set(ACTORS_GLOBAL
        ${ACTORS_DIR}/global/global.cpp
        ${ACTORS_DIR}/global/timing_info.cpp
        ${ACTORS_DIR}/global/message_atoms.cpp
        ${ACTORS_DIR}/global/settings_functions.cpp
        ${ACTORS_DIR}/global/auxiliary.cpp)

    set(SUMMA_ACTOR
        ${ACTORS_DIR}/summa_actor/summa_actor.cpp
        ${ACTORS_DIR}/summa_actor/summa_client.cpp
        ${ACTORS_DIR}/summa_actor/summa_server.cpp
        ${ACTORS_DIR}/summa_actor/summa_backup_server.cpp
        ${ACTORS_DIR}/summa_actor/batch/batch.cpp
        ${ACTORS_DIR}/summa_actor/batch/batch_container.cpp
        ${ACTORS_DIR}/summa_actor/client/client.cpp
        ${ACTORS_DIR}/summa_actor/client/client_container.cpp)

    set(HRU_ACTOR
        ${ACTORS_DIR}/hru_actor/cpp_code/hru_actor.cpp)

    set(FILE_ACCESS_ACTOR
        ${ACTORS_DIR}/file_access_actor/cpp_code/file_access_actor.cpp
        ${ACTORS_DIR}/file_access_actor/cpp_code/forcing_file_info.cpp
        ${ACTORS_DIR}/file_access_actor/cpp_code/output_container.cpp)

    set(JOB_ACTOR
        ${ACTORS_DIR}/job_actor/job_actor.cpp
        ${ACTORS_DIR}/job_actor/GRUinfo.cpp)

    set(MAIN
        ${ACTORS_DIR}/main.cpp)
    
    add_library(SUMMA_NOAHMP OBJECT
        ${NOAHMP}
        ${NRUTIL})
        target_compile_options(SUMMA_NOAHMP PRIVATE ${SUMMA_NOAHMP_OPTIONS})
    add_library(SUMMA_COMM OBJECT
        ${COMM_ALL})
        target_compile_options(SUMMA_COMM PRIVATE ${SUMMA_ALL_OPTIONS})
        target_include_directories(SUMMA_COMM PRIVATE ${SUMMA_INCLUDES})
        target_link_libraries(SUMMA_COMM PUBLIC ${SUMMA_LIBS})

    # Build SUMMA Shared Library
    add_library(summa SHARED
        ${SUMMA_ALL})
    target_compile_options(summa PRIVATE ${SUMMA_ALL_OPTIONS})
    target_include_directories(summa PUBLIC ${SUMMA_INCLUDES})
    target_link_libraries(summa PUBLIC ${SUMMA_LIBS} SUMMA_COMM)

    # Build Summa-Actors executable
    add_executable(${EXEC_NAME}
        ${ACTORS_GLOBAL}
        ${HRU_ACTOR}
        ${FILE_ACCESS_ACTOR}
        ${JOB_ACTOR}
        ${SUMMA_ACTOR}
        ${SUMMA_CLIENT}
        ${SUMMA_SERVER}
        ${MAIN})
        set_property(TARGET ${EXEC_NAME} PROPERTY LINKER_LANGUAGE Fortran)
        target_include_directories(${EXEC_NAME} PUBLIC ${SUMMA_ACTORS_INCLUDES})
        target_link_libraries( ${EXEC_NAME} ${SUMMA_ACTORS_LIBS})
endfunction()