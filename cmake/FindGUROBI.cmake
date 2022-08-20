# ------------------------------------------------------------------------------
#      Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
#      https://3d.bk.tudelft.nl/liangliang/
#
#      This file is part of Easy3D. If it is useful in your research/work,
#      I would be grateful if you show your appreciation by citing it:
#       ------------------------------------------------------------------
#           Liangliang Nan.
#           Easy3D: a lightweight, easy-to-use, and efficient C++
#           library for processing and rendering 3D data. 2018.
#       ------------------------------------------------------------------
#
#      Easy3D is free software; you can redistribute it and/or modify
#      it under the terms of the GNU General Public License Version 3
#      as published by the Free Software Foundation.
#
#      Easy3D is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#      GNU General Public License for more details.
#
#      You should have received a copy of the GNU General Public License
#      along with this program. If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# This file sets up Gurobi for CMake. Once done this will define
#
#   GUROBI_FOUND           - system has GUROBI
#   GUROBI_INCLUDE_DIRS    - the GUROBI include directories
#   GUROBI_LIBRARIES       - Link these to use GUROBI
#
#  In your CMakeLists file, you need to add, e.g. (modify it if necessary):
#        if (GUROBI_FOUND)
#            message(STATUS "Gurobi include dir: " ${GUROBI_INCLUDE_DIRS})
#            message(STATUS "Gurobi libraries: " ${GUROBI_LIBRARIES})
#            target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_GUROBI)
#            target_include_directories(${PROJECT_NAME} PRIVATE ${GUROBI_INCLUDE_DIRS})
#            target_link_libraries(${PROJECT_NAME} PRIVATE ${GUROBI_LIBRARIES})
#        endif()
# ------------------------------------------------------------------------------


# Is it already configured?
if (NOT GUROBI_FOUND)

    # Hardcoded search paths
    set(SEARCH_PATHS_FOR_HEADERS
            "$ENV{GUROBI_HOME}/include"
            "/Library/gurobi951/macos_universal2/include"
            "/Library/gurobi950/mac64/include"
            "/Library/gurobi950/macos_universal2/include"
            "C:\\dev\\gurobi950\\win64\\include"
            "/Library/gurobi901/mac64/include"
            "/home/liangliang/dev/gurobi901/include"
            "C:\\dev\\gurobi901\\win64\\include"
            )

    set(SEARCH_PATHS_FOR_LIBRARIES
            "$ENV{GUROBI_HOME}/lib"
            "/Library/gurobi951/macos_universal2/lib"
            "/Library/gurobi950/mac64/lib"
            "/Library/gurobi950/macos_universal2/lib"
            "C:\\dev\\gurobi950\\win64\\lib"
            "/Library/gurobi901/mac64/lib"
            "/home/liangliang/dev/gurobi901/lib"
            "C:\\dev\\gurobi901\\win64\\lib"
            )

    find_path(GUROBI_INCLUDE_DIR gurobi_c++.h
            PATHS ${SEARCH_PATHS_FOR_HEADERS}
            )


    find_library(GUROBI_C_LIBRARY
            NAMES gurobi95 gurobi90 libgurobi
            PATHS ${SEARCH_PATHS_FOR_LIBRARIES}
            )

    find_library(GUROBI_CXX_LIBRARY_DEBUG
            NAMES gurobi_c++ gurobi_c++mdd2017
            PATHS ${SEARCH_PATHS_FOR_LIBRARIES}
            )

    find_library(GUROBI_CXX_LIBRARY_RELEASE
            NAMES gurobi_c++ gurobi_c++md2017
            PATHS ${SEARCH_PATHS_FOR_LIBRARIES}
            )

    # setup header file directories
    set(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})

    # setup libraries files
    set(GUROBI_LIBRARIES
            debug ${GUROBI_CXX_LIBRARY_DEBUG}
            optimized ${GUROBI_CXX_LIBRARY_RELEASE}
            ${GUROBI_C_LIBRARY}
            )

endif ()

# Check that Gurobi was successfully found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_INCLUDE_DIRS)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARIES)

# Hide variables from CMake-Gui options
mark_as_advanced(
        GUROBI_INCLUDE_DIRS
        GUROBI_INCLUDE_DIR
        GUROBI_LIBRARIES
        GUROBI_CXX_LIBRARY_DEBUG
        GUROBI_CXX_LIBRARY_RELEASE
        GUROBI_C_LIBRARY
)
