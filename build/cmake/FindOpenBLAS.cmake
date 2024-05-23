# Find the OpenBLAS library
find_library(OpenBLAS_LIBRARY NAMES openblas)

# Find the OpenBLAS include directory
find_path(OpenBLAS_INCLUDE_DIR NAMES cblas.h)

if(OpenBLAS_LIBRARY AND OpenBLAS_INCLUDE_DIR)
  add_library(OpenBLAS::OpenBLAS INTERFACE IMPORTED)
  set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${OpenBLAS_LIBRARY}")
else()
  message(FATAL_ERROR "OpenBLAS not found")
endif()