set(PROJECT_NAME 3rd_glew)

add_library(${PROJECT_NAME} STATIC src/glew.c)

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "3rd_party")
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(GLEW_INCLUDE_DIR  ${CMAKE_CURRENT_LIST_DIR}/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${GLEW_INCLUDE_DIR})

target_compile_definitions(${PROJECT_NAME} PUBLIC GLEW_STATIC GLEW_NO_GLU)

set (OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)
message(STATUS "OpenGL libraries: ${OPENGL_LIBRARIES}")
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES})


if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()
