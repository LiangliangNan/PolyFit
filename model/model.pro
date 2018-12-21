#-------------------------------------------------
#
# Project created by QtCreator 2014-11-28T00:30:46
#
#-------------------------------------------------

CONFIG -= qt

TARGET = model
TEMPLATE = lib

DEFINES += MODEL_EXPORTS
win32 { DEFINES += WIN32 WIN64 }

CONFIG(debug, debug|release) { DEFINES += _DEBUG }
CONFIG(release, debug|release) { DEFINES += NDEBUG }


SOURCES += \
    kdtree_search.cpp \
    map_builder.cpp \
    map_cells.cpp \
    map_copier.cpp \
    map_editor.cpp \
    map_enumerator.cpp \
    map_geometry.cpp \
    map_io.cpp \
    map_serializer_obj.cpp \
    map_serializer.cpp \
    map.cpp \
    point_set_io.cpp \
    point_set_serializer_vg.cpp \
    point_set.cpp \
    kdtree/kdTree.cpp

HEADERS += \
    iterators.h \
    kdtree_search.h \
    map_attributes.h \
    map_builder.h \
    map_cells.h \
    map_circulators.h \
    map_copier.h \
    map_editor.h \
    map_enumerator.h \
    map_geometry.h \
    map_io.h \
    map_serializer_obj.h \
    map_serializer.h \
    map.h \
    model_common.h \
    point_set_io.h \
    point_set_serializer_vg.h \
    point_set.h \
    vertex_group.h \
    kdtree/kdTree.h \
    kdtree/PriorityQueue.h \
    kdtree/QueryGrid.h \
    kdtree/vector2D.h \
    kdtree/vector3D.h

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
