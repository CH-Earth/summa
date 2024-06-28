include(FindPackageHandleStandardArgs)

# Find the NetCDF C library and include directory
find_library(NetCDF_C_LIBRARY NAMES netcdf)
find_path(NetCDF_C_INCLUDE_DIR NAMES netcdf.h)

# Find the NetCDF Fortran library
find_library(NetCDF_F90_LIBRARY NAMES netcdff)
find_path(NetCDF_F90_INCLUDE_DIR NAMES netcdf.mod)

set (NetCDF_LIBRARIES ${NetCDF_C_LIBRARY} ${NetCDF_F90_LIBRARY})
set (NetCDF_INCLUDE_DIRS ${NetCDF_C_INCLUDE_DIR} ${NetCDF_F90_INCLUDE_DIR})

find_package_handle_standard_args(NetCDF DEFAULT_MSG
        NetCDF_LIBRARIES NetCDF_INCLUDE_DIRS)

if(NetCDF_FOUND)
    mark_as_advanced (NetCDF_C_LIBRARY NetCDF_C_INCLUDE_DIR NetCDF_F90_LIBRARY NetCDF_DIR)

    add_library(NetCDF::NetCDF INTERFACE IMPORTED)

    set_target_properties(NetCDF::NetCDF  PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${NetCDF_LIBRARIES}")

    message(STATUS "NetCDF incl for all components -- ${NetCDF_INCLUDE_DIRS}")
    message(STATUS "NetCDF lib for all components -- ${NetCDF_LIBRARIES}")
endif()