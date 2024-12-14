/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/ - liangliang.nan@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <QApplication>
#include <QSurfaceFormat>

#include <iostream>
#include <QLocale>
#include <QTranslator>
#include <QDebug>

#include "main_window.h"

#include <time.h>


int main(int argc, char **argv)
{
    srand(time(nullptr));

	//Locale management
	//Force 'English' locale to get a consistent behavior everywhere
	QLocale locale = QLocale(QLocale::English);
	locale.setNumberOptions(QLocale::c().numberOptions());
	QLocale::setDefault(locale);

#ifdef Q_OS_UNIX
	//We reset the numeric locale for POSIX functions
	//See http://qt-project.org/doc/qt-5/qcoreapplication.html#locale-settings
	setlocale(LC_NUMERIC, "C");
#endif

#if (QT_VERSION >= QT_VERSION_CHECK(5, 6, 0) && (QT_VERSION < QT_VERSION_CHECK(6, 0, 0)))
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif

	// Note: Calling QSurfaceFormat::setDefaultFormat() before constructing the
	//       QApplication instance is mandatory on some platforms(for example, macOS)
	//       when an OpenGL core profile context is requested. This is to ensure
	//       that resource sharing between contexts stays functional as all internal
	//       contexts are created using the correct version and profile.
	QSurfaceFormat format = QSurfaceFormat::defaultFormat();
	format.setVersion(4, 3);
	format.setProfile(QSurfaceFormat::CompatibilityProfile);
	format.setDepthBufferSize(24);
	format.setStencilBufferSize(8);
	format.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
	format.setSamples(4);
#ifndef NDEBUG
	format.setOption(QSurfaceFormat::DebugContext);
#endif
	QSurfaceFormat::setDefaultFormat(format);

	QApplication app(argc, argv);

	MainWindow window;	
	window.show();

	return app.exec();
};
