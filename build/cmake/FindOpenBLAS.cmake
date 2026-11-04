# Find the OpenBLAS library
find_library(OpenBLAS_LIBRARY NAMES openblas)

# A cached include directory can vanish, as a versioned macOS SDK does when Xcode updates
if(OpenBLAS_INCLUDE_DIR AND NOT EXISTS "${OpenBLAS_INCLUDE_DIR}")
  unset(OpenBLAS_INCLUDE_DIR CACHE)
endif()

# Find the OpenBLAS include directory, beside the library first so the macOS vecLib headers are not taken for it
if(OpenBLAS_LIBRARY)
  get_filename_component(OpenBLAS_LIBRARY_DIR "${OpenBLAS_LIBRARY}" DIRECTORY)
  find_path(OpenBLAS_INCLUDE_DIR NAMES cblas.h PATHS "${OpenBLAS_LIBRARY_DIR}/../include" PATH_SUFFIXES openblas NO_DEFAULT_PATH)
endif()
find_path(OpenBLAS_INCLUDE_DIR NAMES cblas.h)

if(OpenBLAS_LIBRARY AND OpenBLAS_INCLUDE_DIR)
  add_library(OpenBLAS::OpenBLAS INTERFACE IMPORTED)
  set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${OpenBLAS_LIBRARY}")
else()
  message(FATAL_ERROR "OpenBLAS not found")
endif()
