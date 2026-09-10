#===============================================================================
# Compilation directives for the (optional) OpenWQ framework
#===============================================================================
if(CMAKE_BUILD_TYPE MATCHES Debug)
    SET(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic -fPIC")
else()    
    SET(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic ${OpenMP_CXX_FLAGS} -fPIC")
endif()

# Check if the openwq directory exists relative to the current CMakeLists.txt
message(STATUS "Checking for openwq directory in ${CMAKE_CURRENT_SOURCE_DIR}")
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/openwq")
    message(FATAL_ERROR "OpenWQ files not found. Please check that the openwq \
        repository was cloned into build/source/openwq.")
endif()

SET(OPENWQ_DIR ${CMAKE_CURRENT_SOURCE_DIR}/openwq CACHE INTERNAL "OPENWQ_DIR")
# Source layout mirrors the list in the OpenWQ repo's own CMakeLists.txt
# (build/source/openwq/openwq/CMakeLists.txt). Keep in sync when OpenWQ is updated.
FILE(GLOB OPENWQ_SOURCES
          ${OPENWQ_DIR}/src/*.cpp
          ${OPENWQ_DIR}/src/global/*.cpp
          ${OPENWQ_DIR}/src/extwatflux_ss/*.cpp
          ${OPENWQ_DIR}/src/extwatflux_ss/json_SetConfigInfo/*.cpp
          ${OPENWQ_DIR}/src/initiate/*.cpp
          ${OPENWQ_DIR}/src/readjson/*.cpp
          ${OPENWQ_DIR}/src/readjson/json_utils/*.cpp
          ${OPENWQ_DIR}/src/couplercalls/*.cpp
          ${OPENWQ_DIR}/src/output/*.cpp
          ${OPENWQ_DIR}/src/output/json_SetConfigInfo/*.cpp
          ${OPENWQ_DIR}/src/compute/*.cpp
          ${OPENWQ_DIR}/src/compute/json_SetConfigInfo/*.cpp
          ${OPENWQ_DIR}/src/compute/solver_backeuler/*.cpp
          ${OPENWQ_DIR}/src/compute/solver_sundials/*.cpp
          ${OPENWQ_DIR}/src/units/*.cpp
          ${OPENWQ_DIR}/src/utils/*.cpp
          ${OPENWQ_DIR}/src/utils/models_toolbox/*.cpp
          ${OPENWQ_DIR}/src/models_CH/*.cpp
          ${OPENWQ_DIR}/src/models_CH/bgc_flex/*.cpp
          ${OPENWQ_DIR}/src/models_CH/PHREEQC/*.cpp
          ${OPENWQ_DIR}/src/models_TD/*.cpp
          ${OPENWQ_DIR}/src/models_TD/native_adv/*.cpp
          ${OPENWQ_DIR}/src/models_TD/native_advdisp/*.cpp
          ${OPENWQ_DIR}/src/models_LE/*.cpp
          ${OPENWQ_DIR}/src/models_LE/native_boundmix/*.cpp
          ${OPENWQ_DIR}/src/models_SI/*.cpp
          ${OPENWQ_DIR}/src/models_SI/langmuir/*.cpp
          ${OPENWQ_DIR}/src/models_SI/freundlich/*.cpp
          ${OPENWQ_DIR}/src/models_TS/*.cpp
          ${OPENWQ_DIR}/src/models_TS/hype_hvb/*.cpp
          ${OPENWQ_DIR}/src/models_TS/hype_mmf/*.cpp)

SET(OPENWQ_INCLUDES ${OPENWQ_DIR}/src/ CACHE INTERNAL "OPENWQ_INCLUDES")

set (OPENWQ_COUPLER
    ${CMAKE_CURRENT_SOURCE_DIR}/openWQ.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/summa_openWQ_allocspace.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/summa_openWQ.f90
    CACHE INTERNAL "OPENWQ_COUPLER")

FIND_PACKAGE(Armadillo REQUIRED)
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C)
FIND_PACKAGE(OpenMP)

# Current OpenWQ (develop) hard-depends on PhreeqcRM (geochemistry engine) and on
# SUNDIALS/CVODE for its chemistry/sediment ODE solver (src/compute/solver_sundials).
# These are separate from SUMMA's own optional SUNDIALS (IDA/KINSOL) use.
FIND_PACKAGE(PhreeqcRM REQUIRED)
FIND_PACKAGE(SUNDIALS REQUIRED)

SET(HYDROLINK_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/OpenWQ_hydrolink.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OpenWQ_interface.cpp
    CACHE INTERNAL "HYDROLINK_SOURCES")

ADD_LIBRARY(openWQ OBJECT ${OPENWQ_SOURCES} ${HYDROLINK_SOURCES})
# OpenWQ's headers use the raw HDF5 C API (hid_t etc.), pulled in transitively via
# <armadillo> when Armadillo is built with ARMA_USE_HDF5. H5_USE_110_API keeps the
# 1.10-era signatures Armadillo's HDF5 wrapper expects when building against HDF5 >= 1.12.
TARGET_COMPILE_DEFINITIONS(openWQ PRIVATE H5_USE_110_API)
TARGET_INCLUDE_DIRECTORIES(openWQ PRIVATE
    ${ARMADILLO_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIRS}
    ${PhreeqcRM_INCLUDE_DIRS}
    ${OPENWQ_INCLUDES})
# HDF5: link the actual C library. ${HDF5_C_LIBRARY_hdf5} (used by older revisions of this
# file) is not populated by modern FindHDF5; ${HDF5_C_LIBRARIES} is. The undefined _H5* symbols
# come from Armadillo's HDF5 wrapper (compiled into these objects via ARMA_USE_HDF5).
TARGET_LINK_LIBRARIES(openWQ PUBLIC
    ${ARMADILLO_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${HDF5_C_LIBRARY_hdf5}
    PhreeqcRM::PhreeqcRM
    SUNDIALS::cvode
    SUNDIALS::nvecserial
    SUNDIALS::core)
IF(OpenMP_CXX_FOUND)
    MESSAGE("** OpenWQ **: OpenMP enabled!")
    TARGET_LINK_LIBRARIES(openWQ PUBLIC OpenMP::OpenMP_CXX)
ENDIF()

