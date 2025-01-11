cmake_minimum_required(VERSION 3.12)
project(PyPolyFit)

set(CMAKE_CXX_STANDARD 11)  # If you use pybind11::overload_cast<...>, then it requires compiling in C++14 mode

add_subdirectory(3rd_pybind11)
pybind11_add_module(${PROJECT_NAME} MODULE bindings/polyfit.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE method)

target_compile_definitions(${PROJECT_NAME} PRIVATE PYBIND11_SIMPLE_GIL_SAFE_ABI) # to use the Python Stable ABI

set(POLYFIT_MODULE_NAME "PyPolyFit") # Sets the base name of the output file
target_compile_definitions(${PROJECT_NAME} PUBLIC "PyPolyFit_NAME=${POLYFIT_MODULE_NAME}")

# Platform-specific suffix for Python extensions
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(POLYFIT_MODULE_SUFFIX ".pyd")
else()
    set(POLYFIT_MODULE_SUFFIX ".so")
endif()
# Set the target properties with the computed suffix
set_target_properties(${PROJECT_NAME} PROPERTIES
        OUTPUT_NAME ${POLYFIT_MODULE_NAME}  # Sets the base name of the output file
        SUFFIX "${POLYFIT_MODULE_SUFFIX}"   # Use the computed suffix
        FOLDER "python"
        )
