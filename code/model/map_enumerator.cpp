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


#include "map_enumerator.h"
#include "map.h"
#include "map_geometry.h"
#include "../basic/assertions.h"

#include <stack>

int MapEnumerator::enumerate_vertices(Map* map, Attribute<Map::Vertex, int>& id, int start, int step)
{
	ogf_assert(start >= 0) ;
	ogf_assert(step > 0) ;
	ogf_assert(
		id.attribute_manager() == map->vertex_attribute_manager()
		) ;
	int cur_id = start ;
	FOR_EACH_VERTEX(Map,map, it) {
		id[it] = cur_id ;
		cur_id += step ;
	}
	return (cur_id - start) / step ;
}

int MapEnumerator::enumerate_halfedges(Map* map, Attribute<Map::Halfedge, int>& id, int start, int step) 
{
	ogf_assert(start >= 0) ;
	ogf_assert(step > 0) ;
	ogf_assert(
		id.attribute_manager() == map->halfedge_attribute_manager()
		) ;
	int cur_id = start ;
	FOR_EACH_HALFEDGE(Map,map, it) {
		id[it] = cur_id ;
		cur_id += step ;
	}
	return (cur_id - start) / step ;
}

int MapEnumerator::enumerate_facets(Map* map, Attribute<Map::Facet, int>& id, int start, int step) 
{
	ogf_assert(start >= 0) ;
	ogf_assert(step > 0) ;
	ogf_assert(
		id.attribute_manager() == map->facet_attribute_manager()
		) ;
	int cur_id = start ;
	FOR_EACH_FACET(Map, map, it) {
		id[it] = cur_id ;
		cur_id += step ;
	}
	return (cur_id - start) / step ;
}



static void propagate_planar_component(
	Attribute<Map::Facet, int>& id,
	Map::Facet* f, int cur_id
)
{
	std::stack<Map::Facet*> stack;
	stack.push(f);

	while (!stack.empty()) {
		Map::Facet* top = stack.top();
		stack.pop();
		if (id[top] == -1) {
			id[top] = cur_id;
			const vec3& n_top = Geom::facet_normal(top);
			Map::Halfedge* it = top->halfedge();
			do {
				Map::Facet* cur = it->opposite()->facet();
				if (cur != nil && id[cur] == -1) {
					const vec3& n_cur = Geom::facet_normal(cur);
					float error = (std::abs(dot(n_top, n_cur)) - 1.0f);
					if (std::abs(error) < 1e-5)
						stack.push(cur);
				}
				it = it->next();
			} while (it != top->halfedge());
		}
	}
}


int MapEnumerator::enumerate_planar_components(Map* map, Attribute<Map::Facet, int>& id)
{
	ogf_assert(id.attribute_manager() == map->facet_attribute_manager());
	{ FOR_EACH_FACET(Map, map, it) {
		id[it] = -1;
	}}
	int cur_id = 0;
	{ FOR_EACH_FACET(Map, map, it) {
		if (id[it] == -1) {
			propagate_planar_component(id, it, cur_id);
			cur_id++;
		}
	}}
	return cur_id;
}
