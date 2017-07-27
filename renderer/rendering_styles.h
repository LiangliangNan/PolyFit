#ifndef _RENDERER_RENDERING_STYLES_H_
#define _RENDERER_RENDERING_STYLES_H_
/*
*  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
*  Copyright (C) 2000-2005 INRIA - Project ALICE
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy - levy@loria.fr
*
*     Project ALICE
*     LORIA, INRIA Lorraine,
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX
*     FRANCE
*
*  Note that the GNU General Public License does not permit incorporating
*  the Software into proprietary programs.
*
* As an exception to the GPL, Graphite can be linked with the following (non-GPL) libraries:
*     Qt, SuperLU, WildMagic and CGAL
*/


#include "renderer_common.h"
#include "../basic/color.h"



struct PointStyle {
	PointStyle() : visible(false), color(0.0, 0.0, 0.0, 1.0), size(4) {}
	bool  visible ;
	Color color ;
	float size ;
} ;

std::ostream& operator<<(std::ostream& out, const PointStyle& ps) ;
std::istream& operator>>(std::istream& in, PointStyle& ps) ;

//________________________________________________________

struct EdgeStyle {
	EdgeStyle() : visible(false), color(1.0, 0.0, 0.0, 1.0), width(1) {}
	bool  visible ;
	Color color ;
	float width ;
} ;

std::ostream& operator<<(std::ostream& out, const EdgeStyle& es) ;
std::istream& operator>>(std::istream& in, EdgeStyle& es) ;

//________________________________________________________

struct SurfaceStyle {
	SurfaceStyle() : visible(false), color(0.33f, 0.67f, 1.0f, 0.5f) {}
	bool  visible ;
	Color color ;
} ;

std::ostream& operator<<(std::ostream& out, const SurfaceStyle& ss) ;
std::istream& operator>>(std::istream& in, SurfaceStyle& ss) ;


#endif