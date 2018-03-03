#-------------------------------------------------
#
# Project created by QtCreator 2014-11-28T00:30:46
#
#-------------------------------------------------

CONFIG -= qt

TARGET = math
TEMPLATE = lib

DEFINES += MATH_EXPORTS
win32 { DEFINES += WIN32 WIN64 }

CONFIG(debug, debug|release) { DEFINES += _DEBUG }
CONFIG(release, debug|release) { DEFINES += NDEBUG }


SOURCES += \
    math_types.cpp \
    polygon2d.cpp \
    principal_axes.cpp \
    quaternion.cpp \
    semi_definite_symmetric_eigen.cpp

HEADERS += \
    box.h \
    line.h \
    math_common.h \
    math_types.h \
    matrix.h \
    plane.h \
    polygon2d.h \
    principal_axes.h \
    quaternion.h \
    semi_definite_symmetric_eigen.h \
    vecg.h

INCLUDEPATH += . \
#    $$quote($(CGAL_DIR)/include) \
#    $$quote($(CGAL_DIR)/include) \
     $$PWD/../3rd_party/numeric_stuff/SuiteSparse-4.4.5/SuiteSparse_config \
     $$PWD/../3rd_party/numeric_stuff/SuiteSparse-4.4.5/AMD/Include

DEPENDPATH += .



unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}

## Liangliang: there must be a better way to do this. Please let me know. liangliang.nan@gmail.com
#macx {
#    MAC_SDK  = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk
#    if( !exists( $$MAC_SDK) ) {
#        MAC_SDK  = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk
#        if( !exists( $$MAC_SDK) ) {
#            MAC_SDK  = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk
#            if( !exists( $$MAC_SDK) ) {
#                MAC_SDK  = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk
#                if( !exists( $$MAC_SDK) ) {
#                    error("The selected Mac OSX SDK does not exist at $$MAC_SDK!")
#                }
#            }
#        }
#    }
#    macx:QMAKE_MAC_SDK = $$MAC_SDK
#}



win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../basic/release/ -lbasic
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../basic/debug/ -lbasic
else:unix: LIBS += -L$$OUT_PWD/../basic/ -lbasic

INCLUDEPATH += $$PWD/../basic
DEPENDPATH += $$PWD/../basic
