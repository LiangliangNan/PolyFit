set(PROJECT_NAME 3rd_glew)

add_library(${PROJECT_NAME} STATIC src/glew.c)

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "3rd_party")
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(GLEW_INCLUDE_DIR  ${CMAKE_CURRENT_LIST_DIR}/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${GLEW_INCLUDE_DIR})

target_compile_definitions(${PROJECT_NAME} PUBLIC GLEW_STATIC GLEW_NO_GLU)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
target_compile_definitions(${PROJECT_NAME} PRIVATE GLEW_BUILD PUBLIC GLEW_NO_GLU)

# prefer GLVND
set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)
if (OPENGL_opengl_LIBRARY)
    set(FOUND_GL_LIBRARIES ${OPENGL_opengl_LIBRARY})
else ()
    set(FOUND_GL_LIBRARIES ${OPENGL_gl_LIBRARY})
endif ()
if (OPENGL_glx_LIBRARY)
    list(APPEND FOUND_GL_LIBRARIES ${OPENGL_glx_LIBRARY})
endif ()
message(STATUS "OpenGL libraries: ${FOUND_GL_LIBRARIES}")
target_link_libraries(${PROJECT_NAME} PUBLIC ${FOUND_GL_LIBRARIES})

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()
