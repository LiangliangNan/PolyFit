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


#include "math_types.h"
#include "../basic/basic_types.h"
#include "../basic/assertions.h"



namespace Geom {

	double triangle_beauty(
		const vec2& p1, const vec2& p2, const vec2& p3
		) {
			vec2 p1p2 = p2 - p1 ;
			vec2 p2p3 = p3 - p2 ;
			vec2 p3p1 = p1 - p3 ;
			double l12 = ::sqrt(dot(p1p2, p1p2)) ;
			double l23 = ::sqrt(dot(p2p3, p2p3)) ;
			double l31 = ::sqrt(dot(p3p1, p3p1)) ;
			double d = 0.5 * (l12 + l23 + l31) ;
			return 8.0 * (d - l12) * (d - l23) * (d - l31) / (l12 * l23 * l31) ;
	}

	double triangle_beauty(
		const vec3& p1, const vec3& p2, const vec3& p3
		) {
			vec3 p1p2 = p2 - p1 ;
			vec3 p2p3 = p3 - p2 ;
			vec3 p3p1 = p1 - p3 ;
			double l12 = ::sqrt(dot(p1p2, p1p2)) ;
			double l23 = ::sqrt(dot(p2p3, p2p3)) ;
			double l31 = ::sqrt(dot(p3p1, p3p1)) ;
			double d = 0.5 * (l12 + l23 + l31) ;
			return 
				8.0 * (d - l12) * (d - l23) * (d - l31) / (l12 * l23 * l31) ;
	}


	static inline void perp(const vec2& p1, const vec2& p2, vec2& p, vec2& v) {
		v = p2 - p1 ;
		v = vec2(-v.y, v.x) ;
		p = barycenter(p1,p2) ;
	}


	vec3 perpendicular(const vec3& V) {
		int min_index = 0 ;
		double c = ::fabs(V[0]) ;
		double cur = ::fabs(V[1]) ;
		if(cur < c) {
			min_index = 1 ;
			c = cur ;
		}
		cur = ::fabs(V[2]) ;
		if(cur < c) {
			min_index = 2 ;
		}
		vec3 result ;
		switch(min_index) {
			case 0:
				result = vec3(0, -V.z, V.y) ;
				break ;
			case 1:
				result = vec3(V.z, 0, -V.x) ;
				break ;
			case 2:
				result = vec3(-V.y, V.x, 0) ;
				break ;
		}
		return result ;
	}

	inline double det3x3(
		double a00,  double a01,  double a02,
		double a10,  double a11,  double a12,
		double a20,  double a21,  double a22
		) {
			double m01 = a00*a11 - a10*a01;
			double m02 = a00*a21 - a20*a01;
			double m12 = a10*a21 - a20*a11;
			return m01*a22 - m02*a12 + m12*a02;
	}

	vec3 tetra_circum_center(
		const vec3& p, const vec3& q, 
		const vec3& r, const vec3& s
		) {
			vec3 qp = q-p ;
			double qp2 = length2(qp) ;
			vec3 rp = r-p ;
			double rp2 = length2(rp) ;
			vec3 sp = s-p ;
			double sp2 = length2(sp) ;

			double num_x = det3x3(
				qp.y,qp.z,qp2,
				rp.y,rp.z,rp2,
				sp.y,sp.z,sp2
				);

			double num_y = det3x3(
				qp.x,qp.z,qp2,
				rp.x,rp.z,rp2,
				sp.x,sp.z,sp2
				) ;

			double num_z = det3x3(
				qp.x,qp.y,qp2,
				rp.x,rp.y,rp2,
				sp.x,sp.y,sp2
				) ;

			double den = det3x3(
				qp.x,qp.y,qp.z,
				rp.x,rp.y,rp.z,
				sp.x,sp.y,sp.z
				) ;

			ogf_assert(::fabs(den) > 1e-30) ;

			den *= 2.0 ;

			return vec3(
				p.x + num_x / den,
				p.y - num_y / den,
				p.z + num_z / den 
				) ;

	}


}
