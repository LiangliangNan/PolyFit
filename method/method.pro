#-------------------------------------------------
#
# Project created by QtCreator 2014-11-28T00:30:46
#
#-------------------------------------------------

CONFIG -= qt

TARGET = method
TEMPLATE = lib

DEFINES += METHOD_EXPORTS
win32 { DEFINES += WIN32 WIN64 }

CONFIG(debug, debug|release) { DEFINES += _DEBUG }
CONFIG(release, debug|release) { DEFINES += NDEBUG }



SOURCES += \
    alpha_shape_boundary.cpp \
    alpha_shape_mesh.cpp \
    face_selection_gurobi.cpp \
    face_selection_lpsolve.cpp \
    face_selection.cpp \
    hypothesis_generator.cpp \
    method_global.cpp \
    polyfit_info.cpp
HEADERS += \
    alpha_shape_boundary.h \
    alpha_shape_CGAL4.10_and_earlier.h \
    alpha_shape_CGAL4.11_and_later.h \
    alpha_shape_mesh.h \
    alpha_shape.h \
    cgal_types.h \
    face_selection.h \
    hypothesis_generator.h \
    method_common.h \
    method_global.h \
    polyfit_info.h

INCLUDEPATH += . \
#    $$quote($(CGAL_DIR)/include) \
#    $$quote($(CGAL_DIR)/include) \
#     $$PWD/../3rd_party/numeric_stuff/SuiteSparse-4.4.5/SuiteSparse_config \
#     $$PWD/../3rd_party/numeric_stuff/SuiteSparse-4.4.5/AMD/Include

DEPENDPATH += .



unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}

#Specify which sdk to use in your .pro file as follows:
#macx {
#    QMAKE_MAC_SDK = macosx10.9
#}


win32: {
    INCLUDEPATH += $$quote($(CGAL_DIR)/include) \
                   $$quote($(BOOST_DIR)) \
                   $$quote($(CGAL_DIR)/auxiliary/gmp/include) \
                   $$quote($(GUROBI_DIR)/include)

    DEPENDPATH +=  $$quote($(CGAL_DIR)/include) \
                   $$quote($(BOOST_DIR)) \
                   $$quote($(CGAL_DIR)/auxiliary/gmp/include) \
                   $$quote($(GUROBI_DIR)/include)

#    CONFIG(release, debug|release): LIBS += $(CGAL_DIR)/lib/libCGAL-vc140-mt-4.11.1.lib
#    else:CONFIG(debug, debug|release): LIBS += $(CGAL_DIR)/lib/libCGAL-vc140-mt-gd-4.11.1.lib
    LIBS += -L$(CGAL_DIR)/lib  # CGAL has the auto_link feature
    LIBS += -L$(CGAL_DIR)/auxiliary/gmp/lib/ -llibgmp-10

#    LIBS += $(GUROBI_DIR)/lib/gurobi75.lib
#    CONFIG(release, debug|release): LIBS += $(GUROBI_DIR)/lib/gurobi_c++md2017.lib
#    else:CONFIG(debug, debug|release): LIBS += $(GUROBI_DIR)/lib/gurobi_c++mdd2017.lib
    LIBS += -L$(GUROBI_DIR)/lib   # the actual library files names will be specified in "face_selection_gurobi.cpp"
}
else: {
    INCLUDEPATH += $$PWD/../../../../../../opt/local/include
    DEPENDPATH += $$PWD/../../../../../../opt/local/include

    LIBS += -L$$PWD/../../../../../../opt/local/lib/ -lgmp.10
    LIBS += -L$$PWD/../../../../../../opt/local/lib/ -lCGAL.13.0.1


    # Gurobi
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.11  # seems it does not work with higher version macOS?

    INCLUDEPATH += $$PWD/../../../../../../Library/gurobi752/mac64/include
    DEPENDPATH += $$PWD/../../../../../../Library/gurobi752/mac64/include

#    QMAKE_CXXFLAGS += -std=c++11 -stdlib=libstdc++
#    LIBS += -L$$PWD/../../../../../../Library/gurobi752/mac64/lib/ -lgurobi_stdc++
#    PRE_TARGETDEPS += $$PWD/../../../../../../Library/gurobi752/mac64/lib/libgurobi_stdc++.a

    LIBS += -L$$PWD/../../../../../../Library/gurobi752/mac64/lib/ -lgurobi_g++4.2
    PRE_TARGETDEPS += $$PWD/../../../../../../Library/gurobi752/mac64/lib/libgurobi_g++4.2.a
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


win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../3rd_lpsolve/release/ -l3rd_lpsolve
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../3rd_lpsolve/debug/ -l3rd_lpsolve
else:unix: LIBS += -L$$OUT_PWD/../3rd_lpsolve/ -l3rd_lpsolve

INCLUDEPATH += $$PWD/../3rd_lpsolve
DEPENDPATH += $$PWD/../3rd_lpsolve

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../3rd_lpsolve/release/lib3rd_lpsolve.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../3rd_lpsolve/debug/lib3rd_lpsolve.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../3rd_lpsolve/release/3rd_lpsolve.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../3rd_lpsolve/debug/3rd_lpsolve.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../3rd_lpsolve/lib3rd_lpsolve.a
