cmake_minimum_required(VERSION 3.12)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

set(${PROJECT_NAME}_HEADERS
        polyfit.h
        )

set(${PROJECT_NAME}_SOURCES
        polyfit.cpp
        )

add_library(${PROJECT_NAME} STATIC
        ${${PROJECT_NAME}_SOURCES}
        ${${PROJECT_NAME}_HEADERS}
        )

target_link_libraries(${PROJECT_NAME} PUBLIC easy3d::core easy3d::util linear_program_solver)

include(../cmake/UseCGAL.cmake)