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
# find_package(Qt5 COMPONENTS Core Gui OpenGL Widgets REQUIRED)
find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Xml REQUIRED)

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


# turn on QStringBuilder for more efficient string construction
#	see https://doc.qt.io/qt-5/qstring.html#more-efficient-string-construction
add_definitions( -DQT_USE_QSTRINGBUILDER )

set(${PROJECT_NAME}_HEADERS
    main_window.h
    paint_canvas.h
    dlg/weight_panel_click.h
    dlg/weight_panel_manual.h
    dlg/wgt_render.h
    )

set(${PROJECT_NAME}_SOURCES
    main_window.cpp
    main.cpp
    paint_canvas.cpp
    dlg/weight_panel_click.cpp
    dlg/weight_panel_manual.cpp
    dlg/wgt_render.cpp
    )

set(${PROJECT_NAME}_FORMS
    main_window.ui
    dlg/weight_panel_click.ui
    dlg/weight_panel_manual.ui
    dlg/wgt_render.ui
    )

set(${PROJECT_NAME}_RESOURCES
    Resources/PolyFit.qrc
    Resources/PolyFit.rc
    Resources/PolyFit.icns
    )

if( APPLE )
    # set also the icon
    # This tells CMake where to find and install the ICON file
    set(${PROJECT_NAME}_ICON ${CMAKE_CURRENT_SOURCE_DIR}/Resources/PolyFit.icns)
    set_source_files_properties(${${PROJECT_NAME}_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

    add_executable( ${PROJECT_NAME} MACOSX_BUNDLE ${${PROJECT_NAME}_ICON} ${${PROJECT_NAME}_HEADERS} ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_FORMS} ${${PROJECT_NAME}_RESOURCES} )

    set_source_files_properties(${${PROJECT_NAME}_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

    ## Note: On macOS, high-DPI support is enabled by settings in the Info.plist file. Make sure they are present.
    ## For more information, see: http://doc.qt.io/qt-5/highdpi.html
    ## Enable high DPI support by using my Info.plist file in which I added NSHighResolutionCapable to be true.
    set_target_properties(${PROJECT_NAME} PROPERTIES
            INSTALL_RPATH "@executable_path/../Frameworks"
            MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME}
            MACOSX_BUNDLE_INFO_PLIST "${CMAKE_CURRENT_SOURCE_DIR}/Mac/PolyFit.plist"
            MACOSX_BUNDLE_ICON_FILE PolyFit.icns
            )
else()
    add_executable( ${PROJECT_NAME} ${${PROJECT_NAME}_HEADERS} ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_FORMS} ${${PROJECT_NAME}_RESOURCES} )
endif()


set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "PolyFit")


# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
target_compile_definitions(${PROJECT_NAME} PRIVATE QT_DEPRECATED_WARNINGS)


target_link_libraries( ${PROJECT_NAME} Qt5::Core )
target_link_libraries( ${PROJECT_NAME} Qt5::Gui )
target_link_libraries( ${PROJECT_NAME} Qt5::Widgets )
target_link_libraries( ${PROJECT_NAME} Qt5::OpenGL )
target_link_libraries( ${PROJECT_NAME} Qt5::Xml )

target_link_libraries( ${PROJECT_NAME} basic math model method renderer 3rd_glew 3rd_QGLViewer-2.6.3)