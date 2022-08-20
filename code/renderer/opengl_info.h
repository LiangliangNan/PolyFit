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


#ifndef _OPENGL_GLINFORMATION_H_
#define _OPENGL_GLINFORMATION_H_

#include "renderer_common.h"

#ifdef _WIN32
#include <Windows.h>
#endif

#include <GL/glew.h>
#include <string>



#ifndef NDEBUG
#define ogf_check_gl {\
	GLInfo::check_gl(__FILE__, __LINE__) ;\
}
#else
#define ogf_check_gl
#endif

class RENDERER_API GLInfo
{
public:
	/**
	* Prints the last GL error to the Logger.
	*/
	static void check_gl(const std::string& file, int line) ;

	static std::string gl_vendor() ;
	static std::string gl_renderer() ;
	static std::string gl_version() ;
	static std::string gl_extensions() ;

	static std::string glew_version() ;
	static std::string glsl_version() ;
} ;


#endif


