# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project(PolyFit-Example)

cmake_minimum_required(VERSION 3.1...3.22)

#create_single_source_cgal_program("polyfit_example_without_input_planes.cpp")
#create_single_source_cgal_program("polyfit_example_user_provided_planes.cpp")
#create_single_source_cgal_program("polyfit_example_model_complexty_control.cpp")
#create_single_source_cgal_program("polyfit_example_with_region_growing.cpp")

add_executable(${PROJECT_NAME}
        polyfit_example_without_input_planes.cpp
        )

target_link_libraries(${PROJECT_NAME} PRIVATE easy3d_core easy3d_fileio polyfit)