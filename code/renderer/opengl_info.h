/* ---------------------------------------------------------------------------
 * Copyright (C) 2017 Liangliang Nan <liangliang.nan@gmail.com>
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of PolyFit. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 *
 *     Liangliang Nan and Peter Wonka.
 *     PolyFit: Polygonal Surface Reconstruction from Point Clouds.
 *     ICCV 2017.
 *
 *  For more information:
 *  https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html
 * ---------------------------------------------------------------------------
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


