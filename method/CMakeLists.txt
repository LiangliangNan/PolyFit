cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(method_HEADERS
    alpha_shape_CGAL4.10_and_earlier.h
    alpha_shape_CGAL4.11_and_later.h
    alpha_shape_mesh.h
    alpha_shape.h
    cgal_types.h
    face_selection.h
    hypothesis_generator.h
    method_common.h
    method_global.h
    )

set(method_SOURCES
    alpha_shape_mesh.cpp
    face_selection.cpp
    hypothesis_generator.cpp
    method_global.cpp
    )


#target_include_directories(${PROJECT_NAME} PRIVATE
#    ${POLYFIT_scip_DIR}
#    ${POLYFIT_scip_DIR}/scip
#    )




add_library(${PROJECT_NAME} SHARED ${method_SOURCES} ${method_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
                      FOLDER "PolyFit")


find_package(CGAL REQUIRED)
if(CGAL_FOUND)
    message("Found CGAL.")
    include(${CGAL_USE_FILE})
    message("CGAL_CXX_FLAGS_INIT: ${CGAL_CXX_FLAGS_INIT}")
    add_definitions(${CGAL_CXX_GLAFS_INIT})
    message("CGAL_LIBRARIES: ${CGAL_LIBRARIES}")
    message("CGAL_3RD_PARTY_LIBRARIES: ${CGAL_3RD_PARTY_LIBRARIES}")
    target_link_libraries(${PROJECT_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
else()
    message(FATAL_ERROR "CGAL is required.")
endif()

find_package(GMP REQUIRED)
    if(GMP_FOUND)
    message("GMP_INCLUDE_DIR: ${GMP_INCLUDE_DIR}")
    message("GMP_LIBRARIES: ${GMP_LIBRARIES}")
    target_link_libraries(${PROJECT_NAME} ${GMP_LIBRARIES})
endif()


target_compile_definitions(${PROJECT_NAME} PRIVATE METHOD_EXPORTS)


if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()

target_link_libraries( ${PROJECT_NAME} basic math model)
