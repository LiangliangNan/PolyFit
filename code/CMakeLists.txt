cmake_minimum_required(VERSION 3.12)

if (APPLE)
        set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum OS X deployment version" FORCE)
        message(STATUS "Minimum OS X deployment version: ${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif ()


get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(PolyFit)

# One shouldn't generate the BUILD project directly in the SOURCES folder!
if ( ${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
        if ( NOT SAME_BUILD_AND_SOURCE_FOLDER_WARNING_ALREADY_ISSUED )
                message(FATAL_ERROR "It is not advised to BUILD the binaries directly in the SOURCE folder!\n If you want to proceed with this option, just CONFIGURE the project once again" )
                set( SAME_BUILD_AND_SOURCE_FOLDER_WARNING_ALREADY_ISSUED TRUE )
        endif()
endif()

################################################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
        # The /MP flag is available and works as intended with the MSVC compiler
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
elseif(CMAKE_BUILD_PARALLEL_LEVEL)
        # Use build system parallelism (e.g., `make -j` or `ninja`)
        set(CMAKE_BUILD_PARALLEL_LEVEL 8) # Example: Use 8 threads
endif()

################################################################################

# Detects whether this is a top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
        set(POLYFIT_TOPLEVEL_PROJECT ON)
else()
        set(POLYFIT_TOPLEVEL_PROJECT OFF)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

################################################################################
### Where to put the libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

### Configuration
set(POLYFIT_CODE_DIR       ${CMAKE_CURRENT_LIST_DIR})
set(POLYFIT_INCLUDE_DIR    ${POLYFIT_CODE_DIR}  ${CMAKE_CURRENT_BINARY_DIR})
set(POLYFIT_3RD_PARTY_DIR  ${POLYFIT_CODE_DIR}/3rd_party)

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

if (POLICY CMP0167)
        cmake_policy(SET CMP0167 NEW)
endif ()
find_package(Boost REQUIRED)

################################################################################

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${CMAKE_BINARY_DIR}/lib CACHE
    PATH "Directory where all the .lib files are dumped." FORCE)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_BINARY_DIR}/bin CACHE
    PATH "Directory where .exe and .dll files are dumped." FORCE)

################################################################################

add_subdirectory(3rd_party/QGLViewer)
add_subdirectory(3rd_party/scip)
add_subdirectory(3rd_party/soplex)
add_subdirectory(3rd_party/glew)
add_subdirectory(3rd_party/pybind11)

add_subdirectory(basic)
add_subdirectory(math)
add_subdirectory(method)
add_subdirectory(model)
add_subdirectory(renderer)


add_subdirectory(Examples)
add_subdirectory(PolyFit)

# Python bindings and examples
add_subdirectory(python)

