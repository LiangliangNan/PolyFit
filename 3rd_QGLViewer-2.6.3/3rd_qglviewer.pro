
TEMPLATE = lib
TARGET = 3rd_qglviewer

CONFIG *= qt opengl
QT *= opengl gui core widgets xml

DEFINES += CREATE_QGLVIEWER_DLL

win32 { DEFINES += WIN32 WIN64 }
win32: LIBS += -lopengl32 -lglu32


CONFIG(debug, debug|release) { DEFINES += _DEBUG }
CONFIG(release, debug|release) { DEFINES += NDEBUG }


HEADERS = \
    camera.h \
    config.h \
    constraint.h \
    domUtils.h \
    frame.h \
    ImageInterface.h \
    keyFrameInterpolator.h \
    manipulatedCameraFrame.h \
    manipulatedFrame.h \
    mouseGrabber.h \
    qglviewer.h \
    quaternion.h \
    vec.h \
    VRender/AxisAlignedBox.h \
    VRender/Exporter.h \
    VRender/gpc.h \
    VRender/NVector3.h \
    VRender/Optimizer.h \
    VRender/ParserGL.h \
    VRender/Primitive.h \
    VRender/PrimitivePositioning.h \
    VRender/SortMethod.h \
    VRender/Types.h \
    VRender/Vector2.h \
    VRender/Vector3.h \
    VRender/VRender.h

SOURCES = \
    camera.cpp \
    constraint.cpp \
    frame.cpp \
    ImageInterface.cpp \
    keyFrameInterpolator.cpp \
    manipulatedCameraFrame.cpp \
    manipulatedFrame.cpp \
    mouseGrabber.cpp \
    qglviewer.cpp \
    quaternion.cpp \
    saveSnapshot.cpp \
    vec.cpp \
    VRender/BackFaceCullingOptimizer.cpp \
    VRender/BSPSortMethod.cpp \
    VRender/EPSExporter.cpp \
    VRender/Exporter.cpp \
    VRender/FIGExporter.cpp \
    VRender/gpc.cpp \
    VRender/NVector3.cpp \
    VRender/ParserGL.cpp \
    VRender/Primitive.cpp \
    VRender/PrimitivePositioning.cpp \
    VRender/TopologicalSortMethod.cpp \
    VRender/Vector2.cpp \
    VRender/Vector3.cpp \
    VRender/VisibilityOptimizer.cpp \
    VRender/VRender.cpp

# -----------------------------------
# --  I m a g e I n t e r f a c e  --
# -----------------------------------
FORMS *= \
    ImageInterface.ui \
    VRenderInterface.ui


# ---------------
# --  U n i x  --
# ---------------
unix {
  # Make much smaller libraries (and packages) by removing debugging informations
  QMAKE_CFLAGS_RELEASE -= -g
  QMAKE_CXXFLAGS_RELEASE -= -g
}


# -----------------------
# --  S G I   I r i x  --
# -----------------------
irix-cc|irix-n32 {
  QMAKE_CFLAGS_RELEASE   -= -O3 -O2 -OPT:Olimit=30000
  QMAKE_LFLAGS_RELEASE   -= -O3 -O2 -OPT:Olimit=30000
  QMAKE_CXXFLAGS_RELEASE -= -O3 -O2 -OPT:Olimit=30000
  QMAKE_CFLAGS_RELEASE   *= -IPA -Ofast=IP35
  QMAKE_LFLAGS_RELEASE   *= -IPA -Ofast=IP35
  QMAKE_CXXFLAGS_RELEASE *= -IPA -Ofast=IP35
  QMAKE_CFLAGS           *= -LANG:std
  QMAKE_LFLAGS           *= -LANG:std
  QMAKE_CXXFLAGS         *= -LANG:std
  QMAKE_CFLAGS           *= -woff 1424,3201,1110,1188
  QMAKE_CXXFLAGS         *= -woff 1424,3201,1110,1188
  QMAKE_LIBS_OPENGL      -= -lXi
}


# -------------------
# --  M a c O S X  --
# -------------------
macx|darwin-g++ {
#  # This setting creates a Mac framework. Comment out this line to create a dylib instead.
#  !staticlib: CONFIG *= lib_bundle

#  lib_bundle {
#	FRAMEWORK_HEADERS.version = Versions
#	# Should be $$replace(TRANSLATIONS, .ts, .qm), but 'replace' is only available in Qt 4.3
#	FRAMEWORK_HEADERS.files = $${QGL_HEADERS} qglviewer.icns qglviewer_fr.qm
#	FRAMEWORK_HEADERS.path = Headers
#	QMAKE_BUNDLE_DATA += FRAMEWORK_HEADERS

#	DESTDIR = $${HOME_DIR}/Library/Frameworks/

#	# For a Framework, 'include' and 'lib' do no make sense.
#	# These and prefix will all define the DESTDIR, in that order in case several are defined
#	!isEmpty( INCLUDE_DIR ) {
#	  DESTDIR = $${INCLUDE_DIR}
#	}

#	!isEmpty( LIB_DIR ) {
#	  DESTDIR = $${LIB_DIR}
#	}

#	!isEmpty( PREFIX ) {
#	  DESTDIR = $${PREFIX}
#	}

#	QMAKE_POST_LINK=cd $$DESTDIR/QGLViewer.framework/Headers && (test -L QGLViewer || ln -s . QGLViewer)

#	#QMAKE_LFLAGS_SONAME  = -Wl,-install_name,@executable_path/../Frameworks/
#	#QMAKE_LFLAGS_SONAME  = -Wl,-install_name,

#	# Framework already installed, with includes
#	INSTALLS -= include target
#  } else {
#	#QMAKE_LFLAGS_SONAME  = -Wl,-install_name,libQGLViewer.dylib
#  }
}


# ---------------------
# --  W i n d o w s  --
# ---------------------
win32 {

  # Needed by Intel C++, (icl.exe) so that WINGDIAPI is a defined symbol in gl.h.
  DEFINES *= WIN32

  # Make sure to have C++ files, PentiumPro code, few warnings, add
  # support to RTTI and Exceptions, and generate debug info "program database".
  # Any feedback on these flags is welcome.
  !win32-g++ {
        QMAKE_CXXFLAGS = -TP -GR -Zi
        DEFINES += NOMINMAX
        win32-msvc {
          QMAKE_CXXFLAGS *= -GX
        } else {
          QMAKE_CXXFLAGS *= -EHs
        }
  }
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

SUBDIRS += \
    3rd_qglviewer.pro

DISTFILES += \
    qglviewer_fr.qm \
    3rd_qglviewer.vcproj \
    qglviewer.icns \
    qglviewer-icon.xpm \
    3rd_qglviewer.vcxproj \
    3rd_qglviewer.vcxproj.filters \
    qglviewer_fr.ts
