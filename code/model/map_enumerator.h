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



#ifndef _GEOM_MAP_ENUMERATOR_H_
#define _GEOM_MAP_ENUMERATOR_H_

#include "model_common.h"
#include "map_attributes.h"


class Map;


class MODEL_API MapEnumerator {
public:
	/**
	* returns the number of vertices.
	*/
	static int enumerate_vertices(
		Map* map, Attribute<Map::Vertex, int>& id,
		int start = 0, int step = 1
		) ;

	/**
	* returns the number of halfedges.
	*/
	static int enumerate_halfedges(
		Map* map, Attribute<Map::Halfedge, int>& id,
		int start = 0, int step = 1
		) ;

	/**
	* returns the number of facets.
	*/
	static int enumerate_facets(
		Map* map, Attribute<Map::Facet, int>& id,
		int start = 0, int step = 1
		) ;

	/**
	* returns the number of planar components.
	* Note: assume the model is watertight (or you need glue first)
	*/
	static int enumerate_planar_components(
		Map* map, Attribute<Map::Facet, int>& id
	);
};


#endif

