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


#ifndef _GEOM_MAP_GEOMETRY_H_
#define _GEOM_MAP_GEOMETRY_H_

#include "model_common.h"
#include "map.h"
#include "../math/polygon2d.h"


// Adds some functions related to Map to the Geom namespace.
namespace Geom {

	inline vec3 vector(const Map::Halfedge* h) {
		return h->vertex()->point() - h->prev()->vertex()->point();
	}

	MODEL_API vec3 facet_normal(const Map::Facet* f) ; 
	MODEL_API vec3 vertex_normal(const Map::Vertex* v) ;
	MODEL_API vec3 triangle_normal(const Map::Facet* f) ;

	// I assume the facet is planar
	MODEL_API Plane3d	facet_plane(const Map::Facet* f) ;
	MODEL_API Polygon3d	facet_polygon(const Map::Facet* f) ;

	double MODEL_API facet_area(const Map::Facet* f) ;

	inline double edge_length(const Map::Halfedge* h) {
		return length(vector(h)) ;
	}

	// average edge length around $v$
	inline double average_edge_length(const Map::Vertex* v) {
		Map::Halfedge* cir = v->halfedge();
		unsigned int count = 0 ;
		double total_len = 0;
		do {
			total_len += edge_length(cir);
			count++ ;
			cir = cir->next_around_vertex() ;
		} while (cir != v->halfedge());
		return total_len / count;
	}

	inline double average_edge_length(const Map::Facet* f) {
		Map::Halfedge* cir = f->halfedge();
		unsigned int count = 0 ;
		double total_len = 0;
		do {
			total_len += edge_length(cir);
			count++ ;
			cir = cir->next() ;
		} while (cir != f->halfedge());
		return total_len / count;
	}

	inline double triangle_area(const Map::Facet* f) {
		ogf_assert(f->is_triangle()) ;
		Map::Halfedge* h1 = f->halfedge() ;
		Map::Halfedge* h2 = h1->next() ;
		Map::Halfedge* h3 = h2->next() ;            
		return triangle_area(
			h1->vertex()->point(),
			h2->vertex()->point(),
			h3->vertex()->point()
			) ;
	}

	double MODEL_API map_area(const Map* map) ;

	Box3d MODEL_API bounding_box(const Map* map) ;

	MODEL_API Map*  duplicate(const Map* map);
}


#endif

