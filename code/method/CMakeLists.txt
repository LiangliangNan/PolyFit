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
        reconstruction.h
        )

set(method_SOURCES
        alpha_shape_mesh.cpp
        face_selection.cpp
        hypothesis_generator.cpp
        method_global.cpp
        reconstruction.cpp
        )


add_library(${PROJECT_NAME} SHARED ${method_SOURCES} ${method_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
        FOLDER "PolyFit")

target_include_directories(${PROJECT_NAME} PUBLIC ${POLYFIT_INCLUDE_DIR})
find_package(CGAL REQUIRED)
if (CGAL_FOUND)
    message(STATUS "Found CGAL ${CGAL_VERSION}")
    target_link_libraries(${PROJECT_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
else ()
    message(FATAL_ERROR "CGAL is required.")
endif ()

find_package(GMP REQUIRED)
if (GMP_FOUND)
    message(STATUS "GMP_INCLUDE_DIR: ${GMP_INCLUDE_DIR}")
    message(STATUS "GMP_LIBRARIES: ${GMP_LIBRARIES}")
    target_link_libraries(${PROJECT_NAME} ${GMP_LIBRARIES})
endif ()


target_compile_definitions(${PROJECT_NAME} PRIVATE METHOD_EXPORTS)


if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
            _CRT_SECURE_NO_WARNINGS
            _CRT_SECURE_NO_DEPRECATE
            )
endif ()

target_link_libraries(${PROJECT_NAME} basic math model)

# RPATH settings for macOS and Linux
if (APPLE)
    # For macOS, @loader_path ensures the library looks for dependencies relative to the location of the module.
    set_target_properties(${PROJECT_NAME} PROPERTIES
            BUILD_RPATH "@loader_path"
            INSTALL_RPATH "@loader_path"
            )
elseif (UNIX)
    # For Linux, $ORIGIN ensures the library looks for dependencies relative to the location of the module.
    set_target_properties(${PROJECT_NAME} PROPERTIES
            BUILD_RPATH "$ORIGIN"
            INSTALL_RPATH "$ORIGIN"
            )
endif ()
