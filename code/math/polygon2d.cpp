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



#include "polygon2d.h"
#include "math_types.h"

#include <fstream>
#include <cstdlib> // for qsort



namespace Geom {


	//------------------------ Helpers for convex hull ---------------------
	// Note: this implementation uses qsort, it would be possible to rewrite
	// it using STL's algorithms (would be cleaner and faster)

	inline bool ccw(
		vec2* P, int i, int j, int k
		) {
			const vec2& pi = P[i] ; 
			const vec2& pj = P[j] ;
			const vec2& pk = P[k] ;
			double a = pi.x - pj.x ;
			double b = pi.y - pj.y ;
			double c = pk.x - pj.x ;
			double d = pk.y - pj.y ;
			return (a*d - b*c <= 0) ;
	}

	inline int cmp_sgn(
		const vec2& p1, const vec2& p2
		) {
			double b1 = p1.x - p2.x ;
			if(b1 > 0) {
				return 1 ;
			}
			if(b1 < 0) {
				return -1 ;
			}
			double b2 = p1.y - p2.y ;
			if(b2 > 0) {
				return 1 ;
			}
			if(b2 < 0) {
				return -1 ;
			}
			return 0 ;
	} 

	// The comparison functions used by qsort takes pointer
	// to the things to sort as arguments.
	static int cmpl(vec2* p1, vec2* p2) {
		return cmp_sgn(*p1, *p2) ;
	}

	static int cmph(vec2* p1, vec2* p2) {
		return cmp_sgn(*p2, *p1) ;
	}

	typedef int (*cmpfunc)(const void*, const void*) ;

	static int make_chain(
		vec2* V, int n,
		int (*cmp)(vec2*, vec2*)
		) {
			int i, j, s = 1;
			vec2 t;

			::qsort(V, n, sizeof(vec2), cmpfunc(cmp));

			for (i=2; i<n; i++) {
				for (j=s; j>=1 && ccw(V, i, j, j-1); j--){}
				s = j+1;
				t = V[s]; V[s] = V[i]; V[i] = t;
			}
			return s;
	}

	//------------------------------------------------------------------------------


	void convex_hull(const Polygon2d& PP, Polygon2d& result) {
		result.clear() ;
		int n = PP.size() ;
		vec2* P = new vec2[n+1] ;
		{ for(int i=0; i<n; i++) {
			P[i] = PP[i] ;
		}}
		int u = make_chain(P, n, cmpl);  
		P[n] = P[0];
		int ch = u+make_chain(P+u, n-u+1, cmph);  
		{for(int i=0; i<ch; i++) {
			result.push_back(P[i]) ;
		}}
		delete[] P ;
	}

	void minimum_area_enclosing_rectangle(
		const Polygon2d& PP, 
		vec2& S, vec2& T
		) {

			// Note: this implementation has O(n2) complexity :-(
			// (where n is the number of vertices in the convex hull)
			// If this appears to be a bottleneck, use a smarter
			// implementation with better complexity.

			Polygon2d P ;
			convex_hull(PP, P) ;

			int N = P.size() ;

			// Add the first vertex at the end of P
			P.push_back(P[0]) ;

			double min_area = Numeric::big_double ;

			for(int i=1; i<=N; i++) {
				vec2 Si = P[i] - P[i-1] ;

				if(length2(Si) < 1e-20) {
					continue ;
				}

				vec2 Ti(-Si.y, Si.x) ;
				Si = normalize(Si) ;
				Ti = normalize(Ti) ;
				double s0 =  Numeric::big_double ;
				double s1 = -Numeric::big_double ;
				double t0 =  Numeric::big_double ;
				double t1 = -Numeric::big_double ; 
				for(int j=1; j<N; j++) {
					vec2 D = P[j] - P[0] ;
					double s = dot(Si, D) ;
					s0 = ogf_min(s0, s) ;
					s1 = ogf_max(s1, s) ;
					double t = dot(Ti, D) ;
					t0 = ogf_min(t0, t) ;
					t1 = ogf_max(t1, t) ;
				}
				double area = (s1 - s0) * (t1 - t0) ;
				if(area < min_area) {
					min_area = area ;
					if((s1 - s0) < (t1 - t0)) {
						S = Si ;
						T = Ti ;
					} else {
						S = Ti ;
						T = Si ;
					}
				}
			}
	}


	bool point_is_in_polygon(const Polygon2d& polygon, const vec2& p) {
		bool inside = false;
		std::size_t n = polygon.size();
		for (std::size_t i = 0, j = n - 1; i < n; j = i, ++i) {
			const vec2& u0 = polygon[i];
			const vec2& u1 = polygon[j];  // current edge

			if (((u0.y <= p.y) && (p.y < u1.y)) ||  // U1 is above the ray, U0 is on or below the ray
				((u1.y <= p.y) && (p.y < u0.y)))    // U0 is above the ray, U1 is on or below the ray
			{
				// find x-intersection of current edge with the ray. 
				// Only consider edge crossings on the ray to the right of P.
				double x = u0.x + (p.y - u0.y) * (u1.x - u0.x) / (u1.y - u0.y);
				if (x > p.x)
					inside = !inside;
			}
		}

		return inside;
	}


	// http://astronomy.swin.edu.au/~pbourke/geometry/polyarea/
	vec2 barycenter(const Polygon2d& P) {
		ogf_assert(P.size() > 0) ;

		double A = signed_area(P) ;

		if(::fabs(A) < 1e-30) {
			return P[0] ;
		}

		double x = 0.0 ;
		double y = 0.0 ;
		for(unsigned int i=0; i<P.size(); i++) {
			unsigned int j = (i+1) % P.size() ;
			const vec2& t1 = P[i] ;
			const vec2& t2 = P[j] ;
			double d = (t1.x * t2.y - t2.x * t1.y) ;
			x += (t1.x + t2.x) * d ;
			y += (t1.y + t2.y) * d ;
		}

		return vec2(
			x / (6.0 * A),
			y / (6.0 * A)
			) ;
	}


	vec2 vertices_barycenter(const Polygon2d& P) {
		ogf_assert(P.size() != 0) ;
		double x = 0 ;
		double y = 0 ;
		for(unsigned int i=0; i<P.size(); i++) {
			x += P[i].x ;
			y += P[i].y ;
		}
		x /= double(P.size()) ;
		y /= double(P.size()) ;
		return vec2(x,y) ;
	}

	static inline Sign opposite(Sign sign_in) {
		return Sign(-int(sign_in)) ;
	}

	static inline bool intersect_segments(
		const vec2& p1, const vec2& p2,
		const vec2& q1, const vec2& q2,
		vec2& result
		) {

			vec2 Vp = p2 - p1 ;
			vec2 Vq = q2 - q1 ;
			vec2 pq = q1 - p1 ;

			double a =  Vp.x ;
			double b = -Vq.x ;
			double c =  Vp.y ;
			double d = -Vq.y ;

			double delta = a*d-b*c ;
			if(delta == 0.0) {
				return false ;
			}

			double tp = (d * pq.x -b * pq.y) / delta ;

			result = vec2(
				(1.0 - tp) * p1.x + tp * p2.x,
				(1.0 - tp) * p1.y + tp * p2.y
				) ;

			return true ;
	}

	void save_polygon(const Polygon2d& P, const std::string& file_name) {
		std::ofstream out(file_name.c_str()) ;
		{for(unsigned int i=0; i<P.size(); i++) {
			out << "v " << P[i].x << " " << P[i].y << std::endl ;
			out << "vt " << P[i].x << " " << P[i].y << std::endl ;
		}}
		out << "f " ;
		{for(unsigned int i=0; i<P.size(); i++) {
			out << i+1 << "/" << i+1 << " " ;
		}}
		out << std::endl ;
	}

	// http://astronomy.swin.edu.au/~pbourke/geometry/polyarea/
	double signed_area(const Polygon2d& P) {
		double result = 0 ;
		for(unsigned int i=0; i<P.size(); i++) {
			unsigned int j = (i+1) % P.size() ;
			const vec2& t1 = P[i] ;
			const vec2& t2 = P[j] ;
			result += t1.x * t2.y - t2.x * t1.y ;
		}
		result /= 2.0 ;
		return result ;
	}

}