#		    l i b Q G L V i e w e r
#	C o m p i l a t i o n    c o n f i g u r a t i o n

# Run "qmake; make; make install" to compile and install the library on Unix systems.
# Optional arguments can tune install paths (as in "qmake PREFIX=$HOME"). See doc/download.html for details.

TEMPLATE = lib
TARGET = QGLViewer
VERSION = 2.9.1
CONFIG *= qt opengl warn_on shared thread create_prl rtti no_keywords

QGL_HEADERS = \
	  qglviewer.h \
	  camera.h \
	  manipulatedFrame.h \
	  manipulatedCameraFrame.h \
	  frame.h \
	  constraint.h \
	  keyFrameInterpolator.h \
	  mouseGrabber.h \
	  quaternion.h \
	  vec.h \
	  domUtils.h \
	  config.h

SOURCES = \
	  qglviewer.cpp \
	  camera.cpp \
	  manipulatedFrame.cpp \
	  manipulatedCameraFrame.cpp \
	  frame.cpp \
	  saveSnapshot.cpp \
	  constraint.cpp \
	  keyFrameInterpolator.cpp \
	  mouseGrabber.cpp \
	  quaternion.cpp \
	  vec.cpp

HEADERS *= $${QGL_HEADERS}
DISTFILES *= qglviewer-icon.xpm
DESTDIR = $${PWD}

TRANSLATIONS = qglviewer_fr.ts

QT *= xml opengl

equals (QT_MAJOR_VERSION, 5) {
	QT *= gui widgets
}
equals (QT_MAJOR_VERSION, 6) {
	QT *= gui widgets openglwidgets
}

!isEmpty( QGLVIEWER_STATIC ) {
  CONFIG *= staticlib
}

# -----------------------------------
# --  I m a g e I n t e r f a c e  --
# -----------------------------------
FORMS *= ImageInterface.ui

# ---------------------------------------------
# --  V e c t o r i a l   R e n d e r i n g  --
# ---------------------------------------------
# In case of compilation troubles with vectorial rendering, uncomment this line
# DEFINES *= NO_VECTORIAL_RENDER

contains( DEFINES, NO_VECTORIAL_RENDER ) {
  message( Vectorial rendering disabled )
} else {
  FORMS *= VRenderInterface.ui

  SOURCES *= \
	VRender/BackFaceCullingOptimizer.cpp \
	VRender/BSPSortMethod.cpp \
	VRender/EPSExporter.cpp \
	VRender/Exporter.cpp \
	VRender/FIGExporter.cpp \
	VRender/gpc.cpp \
	VRender/ParserGL.cpp \
	VRender/Primitive.cpp \
	VRender/PrimitivePositioning.cpp \
	VRender/TopologicalSortMethod.cpp \
	VRender/VisibilityOptimizer.cpp \
	VRender/Vector2.cpp \
	VRender/Vector3.cpp \
	VRender/NVector3.cpp \
	VRender/VRender.cpp

  VRENDER_HEADERS = \
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

  HEADERS *= $${VRENDER_HEADERS}
}


# ---------------
# --  U n i x  --
# ---------------
unix {
	CONFIG -= debug debug_and_release
	CONFIG *= release

	# INCLUDE_DIR and LIB_DIR specify where to install the include files and the library.
	# Use qmake INCLUDE_DIR=... LIB_DIR=... , or qmake PREFIX=... to customize your installation.
	isEmpty( PREFIX ) {
		PREFIX_=/usr/local
	} else {
		PREFIX_=$${PREFIX}
	}
	isEmpty( LIB_DIR ) {
		macx|darwin-g++ {
			LIB_DIR_ = /Library/Frameworks
		} else {
			LIB_DIR_ = $${PREFIX_}/lib
		}
	} else {
		LIB_DIR_ = $${LIB_DIR}
	}
	isEmpty( INCLUDE_DIR ) {
		macx|darwin-g++ {
			isEmpty( PREFIX ) {
				INCLUDE_DIR_ = $${PWD}/Library/Developer/Headers
			} else {
				INCLUDE_DIR_ = $${PREFIX}/Headers
			}
		} else {
			INCLUDE_DIR_ = $${PREFIX_}/include
		}
	} else {
		INCLUDE_DIR_ = $${INCLUDE_DIR}
	}
	isEmpty( DOC_DIR ) {
		macx|darwin-g++ {
			isEmpty( PREFIX ) {
				DOC_DIR = $${PWD}/Library/Developer/Shared/Documentation/QGLViewer
			} else {
				DOC_DIR = $${PREFIX}/Shared/Documentation/QGLViewer
			}
		} else {
			DOC_DIR = $${PREFIX_}/share/doc/QGLViewer
		}
	}

	# GLUT for Unix architecture
	!isEmpty( USE_GLUT ) {
		QMAKE_LIBS_OPENGL *= -lglut
	}

	macx|darwin-g++ {
		# GLU is part of the OpenGL framework
	} else {
		QMAKE_LIBS_OPENGL *= -lGLU

		isEmpty( NO_QT_VERSION_SUFFIX ) {
			equals (QT_MAJOR_VERSION, 4) {
				TARGET = $$join(TARGET,,,-qt4)
			}
			equals (QT_MAJOR_VERSION, 5) {
				TARGET = $$join(TARGET,,,-qt5)
			}
			equals (QT_MAJOR_VERSION, 6) {
				TARGET = $$join(TARGET,,,-qt6)
			}
		}
	}

	MOC_DIR = .moc
	OBJECTS_DIR = .obj

	# Adds a -P option so that "make install" as root creates files owned by root and links are preserved.
	# This is not a standard option, and it may have to be removed on old Unix flavors.
	!hpux {
		QMAKE_COPY_FILE = $${QMAKE_COPY_FILE} -P
	}

	# Make much smaller libraries (and packages) by removing debugging informations
	QMAKE_CFLAGS_RELEASE -= -g
	QMAKE_CXXFLAGS_RELEASE -= -g

	# install header
	include.path = $${INCLUDE_DIR_}/QGLViewer
	# Should be $$replace(TRANSLATIONS, .ts, .qm), but 'replace' only appeared in Qt 4.3
	include.files = $${QGL_HEADERS} qglviewer_fr.qm

	# install documentation html
	documentation.path = $${DOC_DIR}
	documentation.files = ../doc/*.html ../doc/*.css  ../doc/*.qch

	# install documentation images
	docImages.path = $${DOC_DIR}/images
	docImages.files = ../doc/images/*

	# install documentation examples
	#docExamples.path = $${DOC_DIR}/examples
	#docExamples.files = ../examples/*../examples/*/*

	# install documentation refManual
	docRefManual.path = $${DOC_DIR}/refManual
	docRefManual.files = ../doc/refManual/*

	# install static library
	#staticlib.extra = make -f Makefile.Release staticlib
	#staticlib.path = $${LIB_DIR_}
	#staticlib.files = lib$${TARGET}.a

	# install library
	target.path = $${LIB_DIR_}

	# "make install" configuration options
	INSTALLS *= target include documentation docImages docRefManual

	# "make uninstall" for all targets
	target.uninstall = @echo "uninstall"
	include.uninstall = @echo "uninstall"  
	documentation.uninstall = @echo "uninstall"  
	docImages.uninstall = @echo "uninstall"  
	docRefManual.uninstall = @echo "uninstall"  
}


# -------------------
# --  M a c O S X  --
# -------------------
macx|darwin-g++ {
	# Default setting creates a Mac framework. Comment out this line to create a dylib instead.
	!staticlib: CONFIG *= lib_bundle

	include.files *= qglviewer.icns

    # Or whatever exists in /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/
    #QMAKE_MAC_SDK = macosx10.15

	lib_bundle {
		FRAMEWORK_HEADERS.version = Versions
		# Should be $$replace(TRANSLATIONS, .ts, .qm), but 'replace' is only available in Qt 4.3
		FRAMEWORK_HEADERS.files = $${QGL_HEADERS} qglviewer.icns qglviewer_fr.qm
		FRAMEWORK_HEADERS.path = Headers
		QMAKE_BUNDLE_DATA += FRAMEWORK_HEADERS

		# So that the path QGLViewer/*.h exists
		QMAKE_POST_LINK=cd $$DESTDIR/QGLViewer.framework/Headers && (test -L QGLViewer || ln -s . QGLViewer)

		# Specific paths for the installation of the framework.
		!isEmpty( LIB_DIR ) {
			target.path = $${LIB_DIR}
		}

		# Framework already contains includes
		INSTALLS -= include
	}

	# GLUT for Mac architecture
	!isEmpty( USE_GLUT ) {
		QMAKE_LIBS_OPENGL -= -lglut
		QMAKE_LIBS_OPENGL *= -framework GLUT -lobjc
	}
}


# ---------------------
# --  W i n d o w s  --
# ---------------------
win32 {
	# Windows requires a debug lib version to link against debug applications
	CONFIG *= debug_and_release build_all

	# Needed by Intel C++, (icl.exe) so that WINGDIAPI is a defined symbol in gl.h.
	DEFINES *= WIN32

	staticlib {
		DEFINES *= QGLVIEWER_STATIC
	} else {
		DEFINES *= CREATE_QGLVIEWER_DLL
	}

	CONFIG *= embed_manifest_dll

	# Use native OpenGL drivers with Qt5.5
	# No longer implicit since the ANGLE driver is now an alternative
	LIBS += -lopengl32 -lglu32

	# TP : C++ source code
	# GR : Enables run-time type information (RTTI).
	# Zi : Generates complete debugging information (removed)
	# EHs : The exception-handling model that catches C++ exceptions only and tells the
	#       compiler to assume that functions declared as extern "C" may throw an exception.
	# FS : Enable parallel compilation
	# Any feedback on these flags is welcome.
	!win32-g++ {
		QMAKE_CXXFLAGS *= -TP -GR
		DEFINES += NOMINMAX
		win32-msvc {
			QMAKE_CXXFLAGS *= -EH -FS
		} else {
			QMAKE_CXXFLAGS *= -EHs
		}
	}
}


build_pass:CONFIG(debug, debug|release) {
  unix: TARGET = $$join(TARGET,,,_debug)
  else: TARGET = $$join(TARGET,,,d)
}
