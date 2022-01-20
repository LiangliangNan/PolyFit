cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(solver_HEADERS
        linear_program.h
        linear_program_solver.h
        linear_program_solver_SCIP.h
        linear_program_solver_common.h
        stop_watch.h
        )

set(solver_SOURCES
        linear_program.cpp
        linear_program_io.cpp
        linear_program_solver.cpp
        linear_program_solver_SCIP.cpp
        stop_watch.cpp
        )


add_library(${PROJECT_NAME} SHARED ${solver_SOURCES} ${solver_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "solvers")

target_compile_definitions(${PROJECT_NAME} PRIVATE SOLVER_EXPORTS)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${POLYFIT_ROOT}
        ${POLYFIT_3RD_PARTY_DIR}/3rd_scip
        ${POLYFIT_3RD_PARTY_DIR}/3rd_scip/scip
        )

target_link_libraries(${PROJECT_NAME} PRIVATE 3rd_scip 3rd_soplex)