cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(math_HEADERS
    box.h
    line.h
    math_common.h
    math_types.h
    matrix.h
    plane.h
    polygon2d.h
    principal_axes.h
    quaternion.h
    semi_definite_symmetric_eigen.h
    vecg.h
    linear_program.h
    linear_program_solver.h
    )

set(math_SOURCES
    math_types.cpp
    polygon2d.cpp
    principal_axes.cpp
    quaternion.cpp
    semi_definite_symmetric_eigen.cpp
    linear_program.cpp
    linear_program_io.cpp
    linear_program_solver.cpp
    linear_program_solver_GLPK.cpp
    linear_program_solver_LPSOLVE.cpp
    linear_program_solver_SCIP.cpp
    linear_program_solver_GUROBI.cpp
    )



add_library(${PROJECT_NAME} SHARED ${math_SOURCES} ${math_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
                      FOLDER "PolyFit")


target_include_directories(${PROJECT_NAME} PRIVATE
  ${POLYFIT_scip_DIR}
  ${POLYFIT_scip_DIR}/scip
  )

target_compile_definitions(${PROJECT_NAME} PRIVATE MATH_EXPORTS)


if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()


include(../cmake/FindGUROBI.cmake)
if (GUROBI_FOUND)
    message(STATUS "Gurobi include dir: " ${GUROBI_INCLUDE_DIRS})
    message(STATUS "Gurobi libraries: " ${GUROBI_LIBRARIES})

    target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_GUROBI)

    target_include_directories(${PROJECT_NAME} PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${GUROBI_LIBRARIES})
endif()


target_link_libraries( ${PROJECT_NAME} PRIVATE basic 3rd_scip 3rd_lpsolve 3rd_glpk 3rd_soplex )
