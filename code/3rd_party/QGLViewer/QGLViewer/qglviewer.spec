%define version_major 2
%define version_minor 9
%define version_revision 1

Name:		libQGLViewer
Version:	%{version_major}.%{version_minor}.%{version_revision}
Release:	1	

Summary:	Qt based OpenGL generic 3D viewer library.
License:	GPL
Group:		Development/C++
Source:		%{name}-%{version}.tar.gz
URL:		http://www.libqglviewer.com
Buildroot:      %{_tmppath}/%{name}-%{version}-buildroot
Vendor:		Alma
Requires: 	libqt >= 4.0.0

%description
libQGLViewer is a C++ library based on Qt that eases the creation of OpenGL 3D viewers. It provides
some of the typical 3D viewer functionalities, such as the possibility to move the camera using the
mouse, which lacks in most of the other APIs. Other features include mouse manipulated frames,
interpolated keyFrames, object selection, stereo display, screenshot saving and much more. It is
used by OpenGL beginners as well as to create complex applications, being fully customizable and
easy to extend.

%package devel
Summary: The libQGLViewer header files, documentation and examples.
Group: Development/Libraries
Requires: %{name} = 0:%{version}


%description devel
This package contains the header files for libQGLViewer. Install this package to develop programs
that uses libQGLViewer. A reference documentation and pedagogical examples are included.

%prep
%define docDir %{_defaultdocdir}/QGLViewer
%define includeDir %{_includedir}/QGLViewer
%define libDir %{_libdir}
%setup -q -n %{name}-%{version}

%build
# if [[ -z "${QTDIR}" ]]
# then
  # echo "QTDIR undefined - Trying to autodetect..."
  # autoDetect=$(locate lib/libqt | head -1 | sed s:"/lib/libqt.*":"":)
  # if [[ -d $autoDetect ]]
  # then
    # export QTDIR=$autoDetect
  # else
    # echo "Compilation error - QTDIR is undefined, unable to run qmake"
    # echo "Use export QTDIR=... ([ba]sh) or setenv QTDIR ... ([t]csh) and re-run"
    # exit 1
  # fi
# fi

# export PATH=${PATH}:${QTDIR}/bin
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${QTDIR}/lib

cd QGLViewer
qmake-qt4
make %{?_smp_mflags}
make staticlib
make clean
#cd ../designerPlugin
#qmake
#make %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT
%{__install} -d $RPM_BUILD_ROOT%{includeDir}
%{__install} --mode=644 QGLViewer/*.h $RPM_BUILD_ROOT%{includeDir}
#%{__install} --mode=644 QGLViewer/qglviewer.cw $RPM_BUILD_ROOT%{includeDir}

%{__install} -d $RPM_BUILD_ROOT%{libDir}
%{__install} --mode=644 QGLViewer/libQGLViewer.so.%{version} $RPM_BUILD_ROOT%{libDir}
ln -s libQGLViewer.so.%{version} $RPM_BUILD_ROOT%{libDir}/libQGLViewer.so.%{version_major}.%{version_minor}
ln -s libQGLViewer.so.%{version} $RPM_BUILD_ROOT%{libDir}/libQGLViewer.so.%{version_major}
ln -s libQGLViewer.so.%{version} $RPM_BUILD_ROOT%{libDir}/libQGLViewer.so
%{__install} --mode=644 QGLViewer/libQGLViewer.a $RPM_BUILD_ROOT%{libDir}

%{__install} -d $RPM_BUILD_ROOT%{includeDir}/designerPlugin
%{__install} --mode=644 designerPlugin/* $RPM_BUILD_ROOT%{includeDir}/designerPlugin

# %{__install} -d $RPM_BUILD_ROOT%{_mandir}/man3
%{__install} -d $RPM_BUILD_ROOT%{docDir}
%{__install} -d $RPM_BUILD_ROOT%{docDir}/refManual
%{__install} -d $RPM_BUILD_ROOT%{docDir}/images
# %{__install} -d $RPM_BUILD_ROOT%{docDir}/examples
# %{__install} -d $RPM_BUILD_ROOT%{docDir}/examples/contribs
# %{__install} doc/man/man3/QGLViewer.3 $RPM_BUILD_ROOT%{_mandir}/man3/
# %{__install} doc/man/man3/qglviewer_* $RPM_BUILD_ROOT%{_mandir}/man3/
%{__install} --mode=644 doc/*.html doc/*.css $RPM_BUILD_ROOT%{docDir}
%{__install} --mode=644 INSTALL README LICENCE CHANGELOG $RPM_BUILD_ROOT%{docDir}
%{__install} --mode=644 doc/refManual/* $RPM_BUILD_ROOT%{docDir}/refManual
%{__install} --mode=644 doc/images/* $RPM_BUILD_ROOT%{docDir}/images
# %{__install} --mode=644 doc/examples/*.html $RPM_BUILD_ROOT%{docDir}/examples
# %{__install} --mode=644 examples/examples.pro $RPM_BUILD_ROOT%{docDir}/examples
# %{__install} --mode=644 examples/contribs/contribs.pro $RPM_BUILD_ROOT%{docDir}/examples/contribs
for dir in $(find examples -type d)
do
  %{__install} -d $RPM_BUILD_ROOT%{docDir}/$dir
  for file in $(find $dir -maxdepth 1 -type f)
  do
    %{__install} --mode=644 $file $RPM_BUILD_ROOT%{docDir}/$dir
  done
done

# Which ROBUST path should be used ?
#%{__install} -d /usr/lib/qt4/plugins/designer/
#%{__install} --mode=644 designerPlugin/*.so /usr/lib/qt4/plugins/designer/

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{libDir}/%{name}.so
%{libDir}/%{name}.so.%{version_major}.%{version_minor}
%{libDir}/%{name}.so.%{version_major}
%{libDir}/%{name}.so.%{version}
%{libDir}/%{name}.a

%files devel
%defattr(-,root,root)
%dir %{includeDir}
%{includeDir}/*.h
#%{includeDir}/qglviewer.cw

%dir %{includeDir}/designerPlugin
%{includeDir}/designerPlugin/*

# %doc %{_mandir}/man3/QGLViewer.3.bz2
# %doc %{_mandir}/man3/qglviewer_*

%dir %{docDir}
%doc %{docDir}/*.html
%doc %{docDir}/*.css
%doc %{docDir}/README
%doc %{docDir}/LICENCE
%doc %{docDir}/INSTALL
%doc %{docDir}/CHANGELOG
%docdir %{docDir}/refManual
%doc %{docDir}/refManual/*
%dir %{docDir}/images
%{docDir}/images/*
%dir %{docDir}/examples
%dir %{docDir}/examples/*
%{docDir}/examples/*/*

%changelog
* Sat Dec 31 2022 Gilles Debunne <contact@libQGLViewer.com> 2.9.1
- Qt5 compilation error fix.

* Sat Dec 31 2022 Gilles Debunne <contact@libQGLViewer.com> 2.9.0
- CMake compilation, Qt6 fixes, high dpi fixes.

* Sun Mar 13 2022 Gilles Debunne <contact@libQGLViewer.com> 2.8.0
- Updates for Qt6 compatibility.

* Fri Nov 17 2019 Gilles Debunne <contact@libQGLViewer.com> 2.7.2
- Update include.

* Fri Nov 17 2017 Gilles Debunne <contact@libQGLViewer.com> 2.7.1
- Fix deprecated message in constructor.

* Wed Jun 14 2017 Gilles Debunne <contact@libQGLViewer.com> 2.7.0
- QGLViewer extends QOpenGLWidget instead of the deprecated QGLWidget.

* Tue Sep 27 2016 Gilles Debunne <contact@libQGLViewer.com> 2.6.4
- Fix ARM build, remove contribs examples' compilation.

* Thu Jul 16 2015 Gilles Debunne <contact@libQGLViewer.com> 2.6.3
- Build library in QGLViewer directory. Simplified examples.pri.

* Fri Jan 23 2015 Gilles Debunne <contact@libQGLViewer.com> 2.6.2
- Fix qreal issue in Quaternion::FromRotationMatrix().

* Thu Jan 22 2015 Gilles Debunne <contact@libQGLViewer.com> 2.6.1
- Bug fix in Camera::pointUnderPixel introduced by the switch to qreal.

* Mon Nov 17 2014 Gilles Debunne <contact@libQGLViewer.com> 2.6.0
- The entire API now consistently uses the qreal (i.e double) type for floating point numbers.

* Thu Nov 6 2014 Gilles Debunne <contact@libQGLViewer.com> 2.5.4
- More minor compilation warning and error fixes.

* Mon Sep 1 2014 Gilles Debunne <contact@libQGLViewer.com> 2.5.3
- Minor compilation warning and error fixes for recent compilers.

* Wed Mar 5 2014 Gilles Debunne <contact@libQGLViewer.com> 2.5.2
- New sceneUpVector and camera's rotation around it.

* Wed Jan 22 2014 Gilles Debunne <contact@libQGLViewer.com> 2.5.1
- Revolve around point reset to scene center.

* Thu Dec 19 2013 Gilles Debunne <contact@libQGLViewer.com> 2.5.0
- Refactor of all mouse bindings methods. New clear bindings methods.

* Sun Nov 17 2013 Gilles Debunne <contact@libQGLViewer.com> 2.4.1
- Using update() instead of updateGL(). Matrix computation is cached. Full reindent.

* Sat May 25 2013 Gilles Debunne <contact@libQGLViewer.com> 2.4.0
- Supports Qt5 (as well as Qt3 and Qt4), source available on github.

* Thu May 24 2012 Gilles Debunne <contact@libQGLViewer.com> 2.3.17
- Compilation of the examples using a framework fixed on Mac. 

* Tue Apr 10 2012 Gilles Debunne <contact@libQGLViewer.com> 2.3.16
- Added -lGLU library when linking on Linux.

* Thu Mar 8 2012 Gilles Debunne <contact@libQGLViewer.com> 2.3.15
- Added a QGLViewer prefix so that headers are correctly when using the framework on mac.

* Fri Mar 2 2012 Gilles Debunne <contact@libQGLViewer.com> 2.3.14
- Added symbol export for designer on windows. Framework installed in user's library folder on mac.

* Fri Feb 24 2012 Gilles Debunne <contact@libQGLViewer.com> 2.3.13
- Handle lower case compilation folder names on windows/Qt 4.8.

* Tue Feb 14 2012 Gilles Debunne <contact@libQGLViewer.com> 2.3.12
- Fixed compilation for windows/Qt 4.8, better designer installation, removed setPhysicalDistanceToScreen.

* Thu Nov 17 2011 Gilles Debunne <contact@libQGLViewer.com> 2.3.11
- Unneeded Qt3 reference removed.

* Wed Jun 14 2011 Gilles Debunne <contact@libQGLViewer.com> 2.3.10
- Patches in designerPlugin and Qt3 support with Visual 2008.

* Sat Dec 4 2010 Gilles Debunne <contact@libQGLViewer.com> 2.3.9
- Vec switched to double.

* Sun Nov 7 2010 Gilles Debunne <contact@libQGLViewer.com> 2.3.8
- Vec operations return double, Support for Qt3 signal/slot cw file removed.

* Thu Oct 28 2010 Gilles Debunne <contact@libQGLViewer.com> 2.3.7
- Fixes for compilation on windows.

* Sat May 29 2010 Gilles Debunne <contact@libQGLViewer.com> 2.3.6
- Minor fixes for Qt3 and gcc 4.3 compilation. No more LD_LIBRARY_PATH tuning.

* Sun Mar 1 2010 Gilles Debunne <contact@libQGLViewer.com> 2.3.5
- Mac compilation improvements, no_keywords CONFIG option.

* Tue Sep 1 2009 Gilles Debunne <contact@libQGLViewer.com> 2.3.4
- Patches for Qt3 compilation.

* Tue Jul 14 2009 Gilles Debunne <contact@libQGLViewer.com> 2.3.3
- New snapshotToClipboard method.

* Tue Jul 7 2009 Gilles Debunne <contact@libQGLViewer.com> 2.3.2
- Fixed tiled snapshot rendering when using screen coordinates.

* Wed Oct 1 2008 Gilles Debunne <contact@libQGLViewer.com> 2.3.1
- Exceptions added to the GPL Open Source license, Qt2 no longer supported.

* Mon Sep 1 2008 Gilles Debunne <contact@libQGLViewer.com> 2.3.0
- New examples, default package compilation using Qt4, new version numbering.

* Tue Aug 28 2007 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.6-3
- make install problem with static lib compilation fixed (and bug submitted to Trolltech).

* Tue Jul 25 2007 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.6-2
- Missing event includes in examples (qt4).

* Tue Jul 4 2007 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.6-1
- Misspelling in ui files (qt3) and missing include in keyboardAndMouse example (qt4).

* Tue May 29 2007 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.6-0
- Camera saves and restores scene center and radius. Camera::getModelViewProjectionMatrix

* Sun Apr 1 2007 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.5-1
- A new clippingPlane example.

* Wed Jan 31 2007 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.5-0
- 1px shift in pointUnderPixel fixed. Minor bugs fixed.

* Wed Jan 23 2007 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.4-2
- Qt4 include naming convention.

* Tue Dec 12 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.4-1
- Black screen bug introduced by DRIVE mode fixed.

* Tue Nov 28 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.4-0
- New DRIVE mouse mode, standardCamera example, minor bug fixes.

* Wed Jul 26 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.3-1
- Missing brace in .pro added.

* Wed Jul 12 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.3-0
- Camera::interpolateTo(). Constraint bug fix. 

* Mon Jun 7 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.2-3
- Path in examples' .pro for minGW debug compilation.

* Mon May 29 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.2-2
- Documentation problem fixed. Compiles with Qt 2.3.

* Wed May 15 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.2-1
- Missing .ui file added. 

* Wed May 14 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.2-0
- Snapshot can be created at an arbitrary size with optional oversampling anti-aliassing.

* Wed Mar 29 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.1-1
- libQGLViewer version added in documentation pages footers.

* Wed Mar 8 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.1-0
- Minor bug fixes in saveSnapshot. Designer plugin installation directory can now be changed.

* Thu Feb 23 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.0-2
- Truncation warnings on windows fixed, patch for a moc bug with Qt 4.1.1 on VC 6 (thanks Juergen). 

* Wed Feb 22 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.0-1
- Warning with Qt 2.3 fixed. Selection problem with certain compilers fixed on BlobWar example. 

* Wed Feb 8 2006 Gilles Debunne <Gilles.Debunne@imag.fr> 2.2.0-0
- New Camera methods, new examples, many major improvements and bug fixes.

* Mon Nov 14 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-8
- VRender dialog bug with Qt4 fixed.

* Fri Oct 21 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-7
- API documentation restored.

* Fri Oct 7 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-6
- "weak definition" error message fixed with gcc4/Qt4.

* Thu Oct 6 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-5
- Spurious QTDIR in .pro caused error messages with Qt4.

* Tue Oct 4 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-4
- Agora example fixed. Minor typos fixed.

* Wed Sep 21 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-3
- Compilation options fixed for Visual 6. Qt2.3 patches in help window.

* Tue Sep 20 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-2
- Compilation options fixed, especially for MinGW.

* Mon Sep 19 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-1
- shared compilation option added for MinGW.

* Sun Sep 18 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.1-0
- saveSnapshot() file dialog problem fixed.

* Wed Sep 14 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.1.0-0
- Now runs with Qt version 2, 3 and 4. Dual licensing.

* Wed Aug 10 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.0.1-0
- Bounding box warning when viewing EPS fixed.

* Thu Jul 25 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.0.0-5
- Unclosed parenthesis fixed in saveSnapshot with Qt 2.3.

* Thu Jul 21 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.0.0-4
- Division by 0 in manipulatedFrame fixed. terrain and agora examples fixed.

* Thu Jul 7 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.0.0-3
- Patches for Qt 3.1, 64 bits architectures and gcc 4.0 compilation.

* Wed Jun 29 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.0.0-2
- float/double warnings fixed for .NET windows compilation.

* Wed Jun 28 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.0.0-1
- VRender code patched for windows compilation.

* Wed Jun 22 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 2.0.0-0
- New major release. API cleaned up and documentation entirely rewritten.

* Wed Apr 28 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.9-6
- zClippingCoef method, .NET compilation, better help window.

* Wed Mar 2 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.9-5
- saveToFile and restoreFromFile renamed. New setStateFileName method.

* Wed Feb 9 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.9-4
- Minor bug fixes in restoreFromFile.

* Wed Jan 13 2005 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.9-3
- Minor bug fixes. displaymessage(), Vec class. New startScreenCoordinatesSystem() z range.

* Wed Nov 23 2004 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.9-2
- Minor bug fixes. Better MouseGrabber and MultiView examples. New sizeHint() method.

* Wed Nov 17 2004 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.9-1
- Many new utility functions. Many new examples. Binary distribution. Fully tested on Linux, Windows and Mac.

* Wed Jul 7 2004 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.8-3
- QAccel abandonned for key bindings. Library installed in /usr instead of /usr/local.

* Mon Jun 14 2004 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.8-2
- Patches for compilation on HP UX architecture. Anaglyph example.

* Wed Jun 9 2004 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.8-1
- Key customization modified. New key and mouse descriptions can be added in the help window.

* Wed May 5 2004 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.7-1
- New select procedure based on GL_SELECT. Frame inverse() method.

* Wed Mar 17 2004 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.6-2
- Doc installed in /usr/local/share/doc instead of /usr/share/doc. Patch for gcc 2.96.

* Wed Feb 3 2004 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.6-1
- Signal mechanism modified, MouseMotion renamed MouseAction, new project URL.

* Wed Dec 24 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.5-8
- New dom syntax, startScreenCoordinatesSystem orientation, bug fixes.

* Wed Nov 26 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.5-7
- ZOOM_ON_REGION added.

* Fri Nov 17 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.5-6
- Minor bug fixes with MouseGrabber.

* Fri Nov 7 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.5-5
- keyboardAndMouse example, minor bug fixes.

* Fri Oct 31 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.5-4
- QTextEdit patch with Qt 2.3.

* Wed Oct 24 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.5-3
- Patch for nVidia bug with anti-aliased fonts.

* Wed Oct 22 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.5-2
- Bug fixes, better help window.

* Thu Oct 2 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.5-1
- GLUT dependency removed. drawText uses Qt.

* Wed Sep 24 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.4-2
- Seg fault on exit fixed. Minor fixes.

* Wed Jul 16 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.4-1
- Mouse bindings configuration

* Wed Jun 25 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.3-1
- Doxygen search engine, FAQ page, Z-buffer display, constraints in KFI.

* Mon May 5 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.2-2
- /usr/lib changed to /usr/local/lib

* Thu Apr 17 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.2-1
- help() uses popup windows. DLL created for windows.

* Wed Apr 10 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.1-1
- A ManipulatedCameraFrame class. double in Quaternion.

* Wed Mar 26 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.3.0-1
- Many changes in the API. Documentation updated. A new MouseGrabber class.

* Wed Mar 19 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.9-3
- ORTHO camera improvements, better default help().

* Wed Mar 5 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.9-2
- Slerp interpolation fixed. Tiny Camera matrix improvements.

* Wed Feb 26 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.9-1
- No more camera referenceFrame, slerp interpolation and new install paths.

* Wed Jan 29 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.8-3
- pixelGLRatio function, minor changes.

* Wed Jan 22 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.8-2
- Minor bug fixes. GL state saving optimized.

* Wed Jan 15 2003 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.8-1
- SPECIAL key disappears. New trackball features.

* Thu Dec 12 2002 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.7-2
- Minor improvements, draw3DText. Mac and Windows compatible release.
- Documentation and examples added to the distribution.

* Wed Dec 4 2002 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.7-1
- KeyFrameInterpolator and EPSRender. Bug fixes and new trackball

* Thu Sep 10 2002 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.6-1
- New features and bug fixes. See CHANGELOG for details. Cleaner spec.

* Thu Jul 25 2002 Xavier Decoret <Xavier.Decoret@imag.fr> 1.2.5-2
- Links with qt-mt (multithread) so it works fine with Mandrake libqt3-devel
- fix the spec file: files were copied in /usr directories during rebuild!

* Tue Jul 16 2002 Gilles Debunne <Gilles.Debunne@imag.fr> 1.2.5-1
- First rpm release
