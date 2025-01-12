set(PROJECT_NAME 3rd_QGLViewer)

include(../../cmake/UseQt.cmake)

find_package(OpenGL REQUIRED)

# QGLViewer target.
set(QGLViewer_SRC
    "QGLViewer/VRender/BackFaceCullingOptimizer.cpp"
    "QGLViewer/VRender/BSPSortMethod.cpp"
    "QGLViewer/VRender/EPSExporter.cpp"
    "QGLViewer/VRender/Exporter.cpp"
    "QGLViewer/VRender/FIGExporter.cpp"
    "QGLViewer/VRender/gpc.cpp"
    "QGLViewer/VRender/NVector3.cpp"
    "QGLViewer/VRender/ParserGL.cpp"
    "QGLViewer/VRender/Primitive.cpp"
    "QGLViewer/VRender/PrimitivePositioning.cpp"
    "QGLViewer/VRender/TopologicalSortMethod.cpp"
    "QGLViewer/VRender/Vector2.cpp"
    "QGLViewer/VRender/Vector3.cpp"
    "QGLViewer/VRender/VisibilityOptimizer.cpp"
    "QGLViewer/VRender/VRender.cpp"
    "QGLViewer/camera.cpp"
    "QGLViewer/constraint.cpp"
    "QGLViewer/frame.cpp"
    "QGLViewer/keyFrameInterpolator.cpp"
    "QGLViewer/manipulatedCameraFrame.cpp"
    "QGLViewer/manipulatedFrame.cpp"
    "QGLViewer/mouseGrabber.cpp"
    "QGLViewer/qglviewer.cpp"
    "QGLViewer/quaternion.cpp"
    "QGLViewer/saveSnapshot.cpp"
    "QGLViewer/vec.cpp")

add_library(${PROJECT_NAME} STATIC ${QGLViewer_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "3rd_party")

target_include_directories(${PROJECT_NAME} PRIVATE ${POLYFIT_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC ${QtLibs} OpenGL::GL OpenGL::GLU)

target_compile_definitions(${PROJECT_NAME} PUBLIC QGLVIEWER_STATIC)

# get a clean windows
if (WIN32 OR MSVC)
  target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX)
endif()