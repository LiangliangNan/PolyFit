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

target_include_directories(${PROJECT_NAME} PRIVATE ${POLYFIT_INCLUDE_DIR})

set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE ${OPENGL_LIBRARIES} basic math model 3rd_glew)

target_compile_definitions(${PROJECT_NAME} PRIVATE RENDERER_EXPORTS)

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()

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
