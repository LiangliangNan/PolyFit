cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

# ------------------------------------------------------------------------------
# Qt
# ------------------------------------------------------------------------------

## we will use cmake autouic/automoc/autorcc features
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set( QT5_ROOT_PATH CACHE PATH "Qt5 root directory (i.e. where the 'bin' folder lies)" )
if ( QT5_ROOT_PATH )
        list( APPEND CMAKE_PREFIX_PATH ${QT5_ROOT_PATH} )
endif()

# find qt5 components
find_package(Qt5 COMPONENTS Core Gui OpenGL Widgets Xml REQUIRED)

# in the case no Qt5Config.cmake file could be found, cmake will explicitly ask the user for the QT5_DIR containing it!
# thus no need to keep additional variables and checks

# Starting with the QtCore lib, find the bin and root directories
get_target_property(QT5_LIB_LOCATION Qt5::Core LOCATION_${CMAKE_BUILD_TYPE})
get_filename_component(QT_BINARY_DIR ${QT5_LIB_LOCATION} DIRECTORY)

# Apple uses frameworks - move up until we get to the base directory to set the bin directory properly
if ( APPLE )
        get_filename_component(QT_BINARY_DIR ${QT_BINARY_DIR} DIRECTORY)
        set(QT_BINARY_DIR "${QT_BINARY_DIR}/bin")

        set( MACDEPLOYQT "${QT_BINARY_DIR}/macdeployqt" )
endif()

# set QT5_ROOT_PATH if it wasn't set by the user
if ( NOT QT5_ROOT_PATH )
        get_filename_component(QT5_ROOT_PATH ${QT_BINARY_DIR} DIRECTORY)
endif()

include_directories(${Qt5OpenGL_INCLUDE_DIRS}
                    ${Qt5Widgets_INCLUDE_DIRS}
                    ${Qt5Core_INCLUDE_DIRS}
                    ${Qt5Gui_INCLUDE_DIRS}
                    ${Qt5Xml_INCLUDE_DIRS}
                   )

# turn on QStringBuilder for more efficient string construction
#	see https://doc.qt.io/qt-5/qstring.html#more-efficient-string-construction
add_definitions( -DQT_USE_QSTRINGBUILDER )


set(qglviewer_HEADERS
    camera.h
    config.h
    constraint.h
    domUtils.h
    frame.h
    ImageInterface.h
    keyFrameInterpolator.h
    manipulatedCameraFrame.h
    manipulatedFrame.h
    mouseGrabber.h
    qglviewer.h
    quaternion.h
    vec.h
    VRender/AxisAlignedBox.h
    VRender/Exporter.h
    VRender/gpc.h
    VRender/NVector3.h
    VRender/Optimizer.h
    VRender/ParserGL.h
    VRender/Primitive.h
    VRender/PrimitivePositioning.h
    VRender/SortMethod.h
    VRender/Types.h
    VRender/Vector2.h
    VRender/Vector3.h
    VRender/VRender.h
    )

set(qglviewer_SOURCES
    camera.cpp
    constraint.cpp
    frame.cpp
    ImageInterface.cpp
    keyFrameInterpolator.cpp
    manipulatedCameraFrame.cpp
    manipulatedFrame.cpp
    mouseGrabber.cpp
    qglviewer.cpp
    quaternion.cpp
    saveSnapshot.cpp
    vec.cpp
    VRender/BackFaceCullingOptimizer.cpp
    VRender/BSPSortMethod.cpp
    VRender/EPSExporter.cpp
    VRender/Exporter.cpp
    VRender/FIGExporter.cpp
    VRender/gpc.cpp
    VRender/NVector3.cpp
    VRender/ParserGL.cpp
    VRender/Primitive.cpp
    VRender/PrimitivePositioning.cpp
    VRender/TopologicalSortMethod.cpp
    VRender/Vector2.cpp
    VRender/Vector3.cpp
    VRender/VisibilityOptimizer.cpp
    VRender/VRender.cpp
    )


add_library(${PROJECT_NAME} SHARED ${qglviewer_SOURCES} ${qglviewer_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
                      FOLDER "3rd_party")


target_include_directories(${PROJECT_NAME} PRIVATE
    ${POLYFIT_qglviewer_DIR}/shared
    ${POLYFIT_qglviewer_DIR}/bfp
    ${POLYFIT_qglviewer_DIR}/bfp/bfp_LUSOL
    ${POLYFIT_qglviewer_DIR}/bfp/bfp_LUSOL/LUSOL
    ${POLYFIT_qglviewer_DIR}/colamd
    )


target_compile_definitions(${PROJECT_NAME} PRIVATE CREATE_QGLVIEWER_DLL)

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()


target_link_libraries( ${PROJECT_NAME} Qt5::Core )
target_link_libraries( ${PROJECT_NAME} Qt5::Gui )
target_link_libraries( ${PROJECT_NAME} Qt5::Widgets )
target_link_libraries( ${PROJECT_NAME} Qt5::OpenGL )
target_link_libraries( ${PROJECT_NAME} Qt5::Xml )


set (OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)
message(STATUS "OpenGL libraries: ${OPENGL_LIBRARIES}")
target_link_libraries( ${PROJECT_NAME} ${OPENGL_LIBRARIES} )
