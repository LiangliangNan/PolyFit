if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

project(libQGLViewer LANGUAGES CXX VERSION 2.9.1)
cmake_minimum_required(VERSION 3.16)

include(../cmake/UseQt.cmake)

find_package(OpenGL REQUIRED)

# QGLViewer target.
set(QGLViewer_SRC
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/BackFaceCullingOptimizer.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/BSPSortMethod.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/EPSExporter.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/Exporter.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/FIGExporter.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/gpc.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/NVector3.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/ParserGL.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/Primitive.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/PrimitivePositioning.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/TopologicalSortMethod.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/Vector2.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/Vector3.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/VisibilityOptimizer.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/VRender/VRender.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/camera.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/constraint.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/frame.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/keyFrameInterpolator.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/manipulatedCameraFrame.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/manipulatedFrame.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/mouseGrabber.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/qglviewer.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/quaternion.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/saveSnapshot.cpp"
    "${PROJECT_SOURCE_DIR}/QGLViewer/vec.cpp")
add_library(3rd_QGLViewer SHARED ${QGLViewer_SRC})
target_include_directories(3rd_QGLViewer INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(3rd_QGLViewer PUBLIC ${QtLibs} OpenGL::GL OpenGL::GLU)