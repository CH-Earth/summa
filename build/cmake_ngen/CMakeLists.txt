cmake_minimum_required(VERSION 3.10)
enable_language( Fortran )

#Get the iso_c_fortran binding module to build as part of this build
add_subdirectory(../iso_c_fortran_bmi ${CMAKE_BINARY_DIR}/iso_c_bmi)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../../cmake/")

project(summabmi VERSION 1.0.0 DESCRIPTION "Summa-Sundials BMI Module Shared Library")

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)

#### Add variables for individual libraries that are used within the *.pc.in files
set(SUMMA_LIB_NAME_CMAKE summabmi)
set(SUMMA_LIB_DESC_CMAKE "Summa-Sundials BMI Module Shared Library")

# Make sure these are compiled with this directive UNCLEAR IF WE NEED THESE
add_compile_definitions(BMI_ACTIVE NGEN_FORCING_ACTIVE NGEN_OUTPUT_ACTIVE NGEN_ACTIVE)


set(ROOT_DIR ../../../)
set(F_MASTER $(ROOT_DIR)/summa)

SET (CMAKE_Fortran_COMPILER gfortran)

include(bmi.cmake)
# --------------------------------------------------------------------------------------------
# When compiling Summa-Sundials it can be compiled in a Cluster or Debug mode, which is set by
# the CMAKE_BUILD_TYPES variables. The default is Release, non-cluster mode.

######### SET THE PATHS TO THE SUNDIALS LIBRARIES AND INCLUDE FILES #########
#############################################################################

# Add options for build type
set(CMAKE_BUILD_TYPES Debug Release Cluster_Debug Cluster)

# Set default configuration type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set compiler flags
if(CMAKE_BUILD_TYPE MATCHES Debug OR CMAKE_BUILD_TYPE MATCHES Cluster_Debug)
    message("Debug build.")
    add_definitions(-DDEBUG)
    set(FLAGS_NOAH -g -O0 -ffree-form -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors)
    set(FLAGS_ALL -g -O0 -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors)
else()
    message("Release build.")
    set(FLAGS_NOAH -O3 -ffree-form -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors)
    set(FLAGS_ALL -O3 -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors)
endif()


# Packages that are required
find_package(netCDF REQUIRED)
find_package(LAPACK REQUIRED)

set(DIR_SUNDIALS= $(ROOT_DIR)/sundials/instdir)
set(INC_SUNDIALS=-I$(DIR_SUNDIALS)/include -I$(DIR_SUNDIALS)/fortran)

# Set include directories
if(CMAKE_BUILD_TYPE MATCHES Cluster OR CMAKE_BUILD_TYPE MATCHES Cluster_Debug)
    set(INCLUDES $(EBROOTNETCDFMINFORTRAN)/include)
    set(LDFLAGS $(EBROOTNETCDFMINFORTRAN)/lib64 -L$(EBROOTOPENBLAS)/lib)
    set(LIBRARIES $(LDFLAGS) -lnetcdff -lopenblas)

    message("\nBuilding IDA\n")

    link_directories(${DIR_SUNDIALS}/lib64)
    set(CMAKE_BUILD_RPATH "${DIR_SUNDIALS}/lib64")
    set(LIB_SUNDIALS $(DIR_SUNDIALS)/lib64 -lsundials_fnvecmanyvector_mod -lsundials_fida_mod -lsundials_fnvecserial_mod -lsundials_fsunlinsoldense_mod -lsundials_fsunmatrixdense_mod)

else()
    set(INCLUDES /opt/local/include -I/opt/local/lib)
    set(LDFLAGS /opt/local/lib)
    set(LIBRARIES $(LDFLAGS) -llapack -lgfortran -lnetcdff -lnetcdf)

    message("\nBuilding IDA\n")

    link_directories(${DIR_SUNDIALS}/lib)
    set(CMAKE_BUILD_RPATH "${DIR_SUNDIALS}/lib")
    set(LIB_SUNDIALS $(DIR_SUNDIALS)/lib -lsundials_fnvecmanyvector_mod -lsundials_fida_mod -lsundials_fnvecserial_mod -lsundials_fsunlinsoldense_mod -lsundials_fsunmatrixdense_mod)


endif()


message("\nBuilding SUMMA\n")
compile_with_bmi($(ROOT_DIR), $(F_MASTER), $(DIR_SUNDIALS), $(INCLUDES), $(LIBRARIES), $(INC_SUNDIALS), $(LIB_SUNDIALS), $(FLAGS_ALL), $(FLAGS_NOAH))


target_link_libraries(summabmi PUBLIC iso_c_bmi)
target_compile_options(summabmi PUBLIC -cpp -DNGEN_ACTIVE)

set_target_properties(summabmi PROPERTIES VERSION ${PROJECT_VERSION})

include(GNUInstallDirs)

install(TARGETS summabmi
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(summabmi.pc.in summabmi.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/summabmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
