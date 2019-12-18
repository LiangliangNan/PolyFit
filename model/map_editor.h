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



#ifndef _GEOM_MAP_EDITOR_H_
#define _GEOM_MAP_EDITOR_H_

#include "model_common.h"
#include "map.h"
#include "map_attributes.h"


class MapComponent;

class MODEL_API MapEditor : public MapMutator {
public:
	MapEditor(Map* target = nil);
	virtual void set_target(Map* target);

	// _________________ CGAL interface __________

	// returns the new edge
	bool can_split_facet(Halfedge* h, Halfedge* g);
	Halfedge* split_facet(Halfedge* h, Halfedge* g);

	/**
	* Insert a new vertex in the edge referred to by h.
	* If triangulate is set to true, triangulate the affected facets.
	*/
	Vertex* split_edge(Halfedge* h, double alpha = 0.5);
	bool can_join_edges(Vertex* v);

	/**
	* It is not allowed to collapse an edge of a
	*  triangle that has its two other edges on the border.
	*/
	// TODO: check if it works for polygons.
	bool can_collapse_edge(Halfedge* h);
	bool collapse_edge(Halfedge* h);

	void erase_facet(Halfedge* h);

    void reorient_facet(Halfedge* h);

    /**
    * Checks wether the two specified half-edges can be glued,
    * from a topological point of view. h0 and h1 should point
    * in reverse direction, and should be on the border.
    */
    bool can_glue(Halfedge* h0, Halfedge* h1);

    /**
    * h0 and h1 should point in reversed direction, and
    * should be on the border.
    */
    bool glue(Halfedge* h0, Halfedge* h1);

	//_____________ copy attributes _______________

	void copy_attributes(Vertex* to, Vertex* from);
	void copy_attributes(Halfedge* to, Halfedge* from);
	void copy_attributes(Facet* to, Facet* from);

protected:

	//_________________ utilities ____________________


	/**
	* Removes faces having only two edges.
	* @param f0 is an interior halfedge of the face to be removed.
	*/
	void remove_null_face(Halfedge* f0);

	/**
	* Checks wheter the vertices pointed by h0 and h1 can be
	* merged. It is called twice by can_glue(), once per
	* orientation of the edges.
	*/
	bool can_merge_vertices(Halfedge* h0, Halfedge* h1);

    /**
    * To be explained by Nico.
    * Note: should be called with both (h0,h1) and with (h1,h0)
    */
    bool orbits_are_compatible(Halfedge* h0, Halfedge* h1);

	/**
	* Checks the existence of an half_edge e such that
	* e->vertex() = v1 and e->opposite()->vertex() = v2
	*/
	bool halfedge_exists_between_vertices(Vertex* v1, Vertex* v2);

	bool halfedges_on_same_vertex(Halfedge* h1, Halfedge* h2);

	bool halfedges_on_same_facet(Halfedge* h1, Halfedge* h2);

protected:
	MapFacetAttribute<bool> border_facet_;
	MapVertexNormal			vertex_normal_;
	MapVertexLock			is_locked_;
};



#endif

