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
FILE(GLOB OPENWQ_SOURCES 
          ${OPENWQ_DIR}/src/*.cpp 
          ${OPENWQ_DIR}/src/global/*.cpp
          ${OPENWQ_DIR}/src/extwatflux_ss/*.cpp
          ${OPENWQ_DIR}/src/initiate/*.cpp
          ${OPENWQ_DIR}/src/chem/*.cpp
          ${OPENWQ_DIR}/src/watertransp/*.cpp
          ${OPENWQ_DIR}/src/readjson/*.cpp
          ${OPENWQ_DIR}/src/couplercalls/*.cpp
          ${OPENWQ_DIR}/src/output/*.cpp
          ${OPENWQ_DIR}/src/solver/*.cpp
          ${OPENWQ_DIR}/src/units/*.cpp
          ${OPENWQ_DIR}/src/utils/*.cpp)

SET(OPENWQ_INCLUDES ${OPENWQ_DIR}/src/ CACHE INTERNAL "OPENWQ_INCLUDES")

set (OPENWQ_COUPLER
    ${CMAKE_CURRENT_SOURCE_DIR}/openWQ.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/summa_openWQ_allocspace.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/summa_openWQ.f90
    CACHE INTERNAL "OPENWQ_COUPLER")

FIND_PACKAGE(Armadillo REQUIRED) 
FIND_PACKAGE(HDF5 REQUIRED)
FIND_PACKAGE(OpenMP)


SET(HYDROLINK_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/OpenWQ_hydrolink.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OpenWQ_interface.cpp 
    CACHE INTERNAL "HYDROLINK_SOURCES")

ADD_LIBRARY(openWQ OBJECT ${OPENWQ_SOURCES} ${HYDROLINK_SOURCES})
TARGET_INCLUDE_DIRECTORIES(openWQ PRIVATE
    ${ARMADILLO_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIRS}
    ${OPENWQ_INCLUDES})
TARGET_LINK_LIBRARIES(openWQ PUBLIC
    ${ARMADILLO_LIBRARIES}
    ${HDF5_C_LIBRARY_hdf5})
IF(OpenMP_CXX_FOUND)
    MESSAGE("** OpenWQ **: OpenMP enabled!")
    TARGET_LINK_LIBRARIES(openWQ PUBLIC OpenMP::OpenMP_CXX)
ENDIF()

