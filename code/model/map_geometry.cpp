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


#include "map_geometry.h"
#include "map_attributes.h"
#include "map_copier.h"
#include "map_builder.h"
#include "../math/polygon2d.h"


namespace Geom {

	vec3 facet_normal(const Map::Facet* f) {
		vec3 result(0,0,0) ;
		Map::Halfedge* cir = f->halfedge();
		do {
			vec3 v0 = vector(cir) ;
			vec3 v1 = vector(cir->prev()->opposite()) ;
			vec3 n = cross(v0, v1) ;
			result = result + n ;
			cir = cir->next() ;
		} while(cir != f->halfedge()) ;
		result = normalize(result) ;
		return result ;
	}


	vec3 vertex_normal(const Map::Vertex* v) {
		vec3 result(0,0,0) ;
		Map::Halfedge* cir = v->halfedge();
		unsigned int count = 0 ;
		do {
			if (!cir->is_border()) {
				count++ ;
				vec3 v0 = vector(cir->next()) ;
				vec3 v1 = vector(cir->opposite());
				vec3 n = cross(v0, v1) ;
				result = result + n ;
			}
			cir = cir->next_around_vertex() ;
		} while (cir != v->halfedge());
		result = normalize(result);
		return result;
	}


	vec3 triangle_normal(const Map::Facet* f){
		ogf_assert(f->is_triangle()) ;
		vec3 result = (
			cross(vector(f->halfedge()->next()),
			vector(f->halfedge()->opposite()))
			) + (
			cross(vector(f->halfedge()),
			vector(f->halfedge()->prev()->opposite()))
			) + ( 
			cross(vector(f->halfedge()->next()->next()),
			vector(f->halfedge()->next()->opposite()))
			) ;
		result = normalize(result);
		return result;
	}

	Plane3d facet_plane(const Map::Facet* f) {
		return Plane3d(
			f->halfedge()->vertex()->point(), 
			f->halfedge()->next()->vertex()->point(), 
			f->halfedge()->next()->next()->vertex()->point()
			);
	}

	Polygon3d facet_polygon(const Map::Facet* f) {
		Polygon3d plg ;
		Map::Halfedge* cir = f->halfedge();
		do {
			plg.push_back(cir->vertex()->point());
			cir = cir->next() ;
		} while(cir != f->halfedge()) ;
		return plg;
	}

	/*
	// I do not trust this one for the moment ...
	double facet_area(const Map::Facet* f) {
	vec3 n = facet_normal(f) ;
	vec3 w(0,0,0) ;
	Map::Halfedge* it = f->halfedge() ;
	do {
	vec3 v1(
	it-> vertex()-> point().x,
	it-> vertex()-> point().y,
	it-> vertex()-> point().z
	) ;
	vec3 v2(
	it-> next()-> vertex()-> point().x,
	it-> next()-> vertex()-> point().y,
	it-> next()-> vertex()-> point().z
	) ;
	w = w + (v1 ^ v2) ;
	it = it->next() ;
	} while(it != f->halfedge()) ;
	return 0.5 * ::fabs(w * n) ;
	}
	*/

	double facet_area(const Map::Facet* f) {
		double result = 0 ;
		Map::Halfedge* h = f->halfedge() ;
		const vec3& p = h->vertex()->point() ;
		h = h->next() ;
		do {
			result += triangle_area(
				p,
				h->vertex()->point(),
				h->next()->vertex()->point() 
				) ;
			h = h->next() ;
		} while(h != f->halfedge()) ;
		return result ;
	}


	double border_length(Map::Halfedge* start) {
		ogf_assert(start->is_border()) ;
		double result = 0 ;
		Map::Halfedge* cur = start ;
		do {
			result += edge_length(cur) ;
			cur = cur->next() ;
		} while(cur != start) ;
		return result ;
	}

	double map_area(const Map* map) {
		double result = 0 ;
		FOR_EACH_FACET_CONST(Map, map, it) {
			result += facet_area(it) ;
		}
		return result ;
	}

	Box3d bounding_box(const Map* map) {
		ogf_debug_assert(map->size_of_vertices() > 0);
		Box3d result ;
		FOR_EACH_VERTEX_CONST(Map, map, it) {
			result.add_point(it->point()) ;
		}
		return result ;
	}

	Map* duplicate(const Map* map) {
		if (!map)
			return nil;

		Map* result = new Map;

		if(result != nil) {
			MapBuilder builder(result) ;
			MapCopier copier ;
			copier.set_copy_all_attributes(true) ;
			builder.begin_surface() ;
			copier.copy(builder, const_cast<Map*>(map)) ;
			builder.end_surface() ;
		}
		return result ;
	}
}
