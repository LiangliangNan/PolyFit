#-------------------------------------------------
#
# Project created by QtCreator 2014-11-28T00:30:46
#
#-------------------------------------------------

CONFIG -= qt

TARGET = renderer
TEMPLATE = lib

DEFINES += RENDERER_EXPORTS GLEW_BUILD
win32: { DEFINES += WIN32 WIN64 }
win32: { LIBS += -lopengl32 -lglu32 }
else: macx { LIBS += -lGL -lGLU -L/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries }


INCLUDEPATH += . \
    $$quote($(OPENGL_DIR)/include)


CONFIG(debug, debug|release) { DEFINES += _DEBUG }
CONFIG(release, debug|release) { DEFINES += NDEBUG }


SOURCES += \
    opengl_info.cpp \
    point_set_render.cpp \
    surface_render.cpp \
    glew.c

HEADERS += \
    glew.h \
    opengl_info.h \
    point_set_render.h \
    renderer_common.h \
    rendering_styles.h \
    surface_render.h


INCLUDEPATH += . \
#    $$quote($(CGAL_DIR)/include) \
#    $$quote($(CGAL_DIR)/include) \
#    $$PWD/../3rd_party/numeric_stuff/SuiteSparse-4.4.5/SuiteSparse_config \
#    $$PWD/../3rd_party/numeric_stuff/SuiteSparse-4.4.5/AMD/Include

DEPENDPATH += .



unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}


win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../basic/release/ -lbasic
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../basic/debug/ -lbasic
else:unix: LIBS += -L$$OUT_PWD/../basic/ -lbasic

INCLUDEPATH += $$PWD/../basic
DEPENDPATH += $$PWD/../basic


win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../math/release/ -lmath
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../math/debug/ -lmath
else:unix: LIBS += -L$$OUT_PWD/../math/ -lmath

INCLUDEPATH += $$PWD/../math
DEPENDPATH += $$PWD/../math

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../model/release/ -lmodel
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../model/debug/ -lmodel
else:unix: LIBS += -L$$OUT_PWD/../model/ -lmodel

INCLUDEPATH += $$PWD/../model
DEPENDPATH += $$PWD/../model
