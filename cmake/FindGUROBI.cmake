#
# Try to find GUROBI
# Once done this will define
#
# GUROBI_FOUND           - system has GUROBI
# GUROBI_INCLUDE_DIRS    - the GUROBI include directories
# GUROBI_LIBRARIES       - Link these to use GUROBI
#

# Hardcoded search paths
set(SEARCH_PATHS_FOR_HEADERS
    "$ENV{GUROBI_HOME}/include"
    "/Library/gurobi810/mac64/include"
    "C:\\dev\\gurobi810\\win64\\include"
)

set(SEARCH_PATHS_FOR_LIBRARIES
    "$ENV{GUROBI_HOME}/lib"
    "/Library/gurobi810/mac64/lib"
    "C:\\dev\\gurobi810\\win64\\lib"
)

find_path(GUROBI_INCLUDE_DIR gurobi_c++.h
  PATHS ${SEARCH_PATHS_FOR_HEADERS}
)


find_library( GUROBI_C_LIBRARY
              NAMES gurobi81
              PATHS ${SEARCH_PATHS_FOR_LIBRARIES}
              )


find_library( GUROBI_CXX_LIBRARY
              NAMES gurobi_c++ gurobi_c++md2017
              PATHS ${SEARCH_PATHS_FOR_LIBRARIES}
              )

# Check that Gurobi was successfully found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_INCLUDE_DIR)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_C_LIBRARY)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_CXX_LIBRARY)

set(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})
set(GUROBI_LIBRARIES ${GUROBI_C_LIBRARY} ${GUROBI_CXX_LIBRARY} )

# Hide variables from CMake-Gui options
mark_as_advanced(GUROBI_LIBRARIES GUROBI_INCLUDE_DIRS GUROBI_INCLUDE_DIR)
