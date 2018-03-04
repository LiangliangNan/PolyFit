#-------------------------------------------------
#
# Project created by QtCreator 2014-11-28T00:19:41
#
#-------------------------------------------------

CONFIG -= qt

TARGET = basic
TEMPLATE = lib

DEFINES += BASIC_EXPORTS


win32 { DEFINES += WIN32 WIN64}


#Liangliang： I don't understand why sometimes this is not needed.
win32 { LIBS += -lshell32 }

CONFIG(debug, debug|release) { DEFINES += _DEBUG }
CONFIG(release, debug|release) { DEFINES += NDEBUG }


SOURCES += \
    assertions.cpp \
    attribute_adapter.cpp \
    attribute_life_cycle.cpp \
    attribute_manager.cpp \
    attribute_serializer.cpp \
    attribute_store.cpp \
    basic_types.cpp \
    counted.cpp \
    file_utils.cpp \
    logger.cpp \
    progress.cpp \
    rat.cpp \
    raw_attribute_store.cpp \
    stop_watch.cpp

HEADERS += \
    assertions.h \
    attribute_adapter.h \
    attribute_copier.h \
    attribute_life_cycle.h \
    attribute_manager.h \
    attribute_serializer.h \
    attribute_store.h \
    attribute.h \
    basic_common.h \
    basic_types.h \
    canvas.h \
    color.h \
    counted.h \
    dlist.h \
    file_utils.h \
    generic_attributes_io.h \
    line_stream.h \
    logger.h \
    pointer_iterator.h \
    progress.h \
    rat.h \
    raw_attribute_store.h \
    record_id.h \
    smart_pointer.h \
    stop_watch.h

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

