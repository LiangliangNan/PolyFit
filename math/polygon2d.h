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



#ifndef _MATH_GEOMETRY_POLYGON2D_H_
#define _MATH_GEOMETRY_POLYGON2D_H_

#include "math_common.h"
#include "math_types.h"
#include <vector>





namespace Geom {

	void MATH_API save_polygon(const Polygon2d& P, const std::string& file_name) ;

	double MATH_API signed_area(const Polygon2d& P) ;

	inline double area(const Polygon2d& P) { 
		return ::fabs(signed_area(P)) ; 
	}

	void MATH_API convex_hull(const Polygon2d& P, Polygon2d& result) ;

	/**
	* V1 and V2 are the normalized axes of the
	* minimum area enclosing rectangle.
	*/
	void MATH_API minimum_area_enclosing_rectangle(
		const Polygon2d& P, vec2& V1, vec2& V2
		) ;

	// NOTE: works for both convex and non-convex polygons.
	bool MATH_API point_is_in_polygon(const Polygon2d& P, const vec2& p);

	/** 
	* Note: the barycenter of a polygon is not 
	* the barycenter of its vertices 
	*/
	vec2 MATH_API barycenter(const Polygon2d& P) ;

	vec2 MATH_API vertices_barycenter(const Polygon2d& P) ;

	bool MATH_API polygon_is_convex(const Polygon2d& P) ;

}



#endif
