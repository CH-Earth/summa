cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(summa VERSION 3.2.0 LANGUAGES Fortran)

# Add options for build type
set(CMAKE_CONFIGURATION_TYPES Release Debug)

# Set Default Executable Name
set(EXEC_NAME summa.exe)
set(F_MASTER "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(EXEC_DIR "${F_MASTER}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXEC_DIR})

# Use Find files for NetCDF and OpenBLAS
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/cmake/")

# NetCDF is found with a custom FindNetCDF.cmake file
find_package(NetCDF REQUIRED)
list(APPEND EXT_TARGETS NetCDF::NetCDF)

# Attempt to find LAPACK first
find_package(LAPACK QUIET)
if(LAPACK_FOUND)
  message(STATUS "Using LAPACK")
  list(APPEND EXT_TARGETS ${LAPACK_LIBRARIES})
else()
  # If LAPACK not found, fallback to OpenBLAS
  message(STATUS "LAPACK not found, trying OpenBLAS")
  set(BLA_VENDOR OpenBLAS)
  find_package(OpenBLAS REQUIRED)
  list(APPEND EXT_TARGETS OpenBLAS::OpenBLAS)
endif()


# Set Compiler Flags
set(FLAGS_OPT $ENV{FLAGS_OPT}) # get optional user-specified flags from environment variables
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("\nSetting SUMMA Debug Options")
    add_compile_definitions(DEBUG)
    set(FLAGS_NOAH -g -O0 -fbacktrace -fbounds-check -ffree-form -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors ${FLAGS_OPT})
    set(FLAGS_ALL -g -O0 -fbacktrace -fbounds-check -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors -cpp ${FLAGS_OPT})
    set(FLAGS_CXX -g -O0 -fbounds-check -Wfatal-errors -std=c++17 ${FLAGS_OPT})
else()
    message("\nSetting SUMMA Release Options")
    set(FLAGS_NOAH -O3 -ffree-form -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors ${FLAGS_OPT})
    set(FLAGS_ALL  -O3 -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors -cpp ${FLAGS_OPT})
    set(FLAGS_CXX  -O3 -Wfatal-errors -std=c++17 ${FLAGS_OPT})
endif()

#=========================================================================================
# COMPILE PART 1: Define directory paths
#=========================================================================================

# Define directories that contains source code
set(DRIVER_DIR ${F_MASTER}/build/source/driver)
set(DSHARE_DIR ${F_MASTER}/build/source/dshare)
set(ENGINE_DIR ${F_MASTER}/build/source/engine)
set(HOOKUP_DIR ${F_MASTER}/build/source/hookup)
set(NETCDF_DIR ${F_MASTER}/build/source/netcdf)
set(NOAHMP_DIR ${F_MASTER}/build/source/noah-mp)

#=========================================================================================
# COMPILE PART 2: Assemble all of the SUMMA sub-routines
#=========================================================================================

# SUMMA Source Files are defined in the CMakeLists.txt file in the subdirectory
add_subdirectory(${F_MASTER}/build/source/)

#=========================================================================================
# COMPILE PART 3: Collect the subroutines into build groups depending on build type
#=========================================================================================

set(COMM_ALL ${NRPROC} ${HOOKUP} ${DATAMS} ${UTILMS})
set(SUMMA_ALL ${NETCDF} ${PRELIM} ${MODRUN} ${SOLVER} ${DRIVER})
# Add non-actor files
set(SUMMA_ALL ${SUMMA_ALL} ${PRELIM_NOT_ACTORS} ${MODRUN_NOT_ACTORS}
              ${SOLVER_NOT_ACTORS} ${DRIVER_NOT_ACTORS})
set(MAIN_SUMMA ${DRIVER_DIR}/summa_driver.f90)

# Define version number, not working correctly
set(VERSIONFILE     ${DRIVER_DIR}/summaversion.inc)
execute_process(COMMAND "    ${GIT_EXECUTABLE} tag | tail -n 1" OUTPUT_VARIABLE VERSION)
execute_process(COMMAND "date" OUTPUT_VARIABLE BULTTIM)
execute_process(COMMAND "    ${GIT_EXECUTABLE} describe --long --all --always | sed -e's/heads\///'" OUTPUT_VARIABLE GITBRCH)
execute_process(COMMAND "    ${GIT_EXECUTABLE} rev-parse HEAD" OUTPUT_VARIABLE GITHASH)

#=========================================================================================
# COMPILE PART 4: Do the compilation
#=========================================================================================

# update version information, not working correctly
file(WRITE  ${VERSIONFILE} "character(len=64), parameter     :: summaVersion = '${VERSION}'\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: buildTime = ''\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: gitBranch = '${GITBRCH}'\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: gitHash = '${GITHASH}'")


# Build SUMMA_NOAHMP Object
add_library(SUMMA_NOAHMP OBJECT ${NOAHMP} ${NRUTIL})
target_compile_options(SUMMA_NOAHMP PRIVATE ${FLAGS_NOAH})

# Build SUMMA_COMM Object
add_library(SUMMA_COMM OBJECT ${COMM_ALL})
target_compile_options(SUMMA_COMM PRIVATE ${FLAGS_ALL})
target_link_libraries(SUMMA_COMM PUBLIC SUMMA_NOAHMP ${EXT_TARGETS}) # added flags to the link step

add_library(summa SHARED ${SUMMA_ALL})
target_compile_options(summa PRIVATE ${FLAGS_ALL})
target_link_libraries(summa PUBLIC ${EXT_TARGETS} SUMMA_NOAHMP SUMMA_COMM)

add_executable(${EXEC_NAME} ${MAIN_SUMMA})
set_property(TARGET ${EXEC_NAME} PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(${EXEC_NAME} summa ${EXT_TARGETS})