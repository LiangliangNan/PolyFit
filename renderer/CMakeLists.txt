cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(renderer_HEADERS
    opengl_info.h
    point_set_render.h
    renderer_common.h
    rendering_styles.h
    surface_render.h
    )

set(renderer_SOURCES
    opengl_info.cpp
    point_set_render.cpp
    surface_render.cpp
    )


add_library(${PROJECT_NAME} SHARED ${renderer_SOURCES} ${renderer_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
                      FOLDER "PolyFit")

target_include_directories(${PROJECT_NAME} PRIVATE ${POLYFIT_glew_DIR}/include)

find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE ${OPENGL_LIBRARIES} basic math model 3rd_glew)


target_compile_definitions(${PROJECT_NAME} PRIVATE RENDERER_EXPORTS)


if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()
