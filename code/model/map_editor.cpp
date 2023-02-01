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


#include "map_editor.h"
#include "map_geometry.h"
#include "../basic/assertions.h"

#include <set>
#include <cfloat>



MapEditor::MapEditor(Map* target) : MapMutator(target) {
	if (target != nil) {
		set_target(target);
	}
}

void MapEditor::set_target(Map* target) {
	MapMutator::set_target(target);

	if (vertex_normal_.is_bound())
		vertex_normal_.unbind();
	if (MapVertexNormal::is_defined(target))
		vertex_normal_.bind(target);

	if (is_locked_.is_bound())
		is_locked_.unbind();
	if (MapVertexLock::is_defined(target))
		is_locked_.bind(target);
}

//__________ CGAL interface __________________________________

bool MapEditor::can_split_facet(
	Halfedge* h, Halfedge* g
	)
{
	if (h == g) {
		return false;
	}

	if (h->facet() == nil)
		h = h->next()->opposite();
	if (g->facet() == nil)
		g = g->next()->opposite();

	if (!halfedges_on_same_facet(h, g)) {
		return false;
	}
	if (h->next() == g || g->next() == h) {
		return false;
	}
	return true;
}

Map::Halfedge* MapEditor::split_facet(Halfedge* h, Halfedge* g)
{
	if (!can_split_facet(h, g)) {
		return nil;
	}

	if (h->facet() == nil)
		h = h->next()->opposite();
	if (g->facet() == nil)
		g = g->next()->opposite();

	Map::Facet* f = h->facet();

	Halfedge* result = new_edge();

	make_sequence(result->opposite(), g->next());
	make_sequence(result, h->next());
	make_sequence(g, result);
	make_sequence(h, result->opposite());

	set_halfedge_vertex(result, h->vertex());
	set_halfedge_vertex(result->opposite(), g->vertex());

	make_facet_key(result->opposite(), h->facet());

	set_facet_on_orbit(result, new_facet(f));
	make_facet_key(result);

	return result;
}


MapEditor::Vertex* MapEditor::split_edge(Halfedge* h, double alpha) {
	Vertex* ov1 = h->vertex();
	Vertex* ov2 = h->opposite()->vertex();

	vec3 p = h->vertex()->point()*(1. - alpha) + h->opposite()->vertex()->point()*alpha;

	Vertex* v = new_vertex();

	Halfedge* z1 = h->next();
	Halfedge* z2 = h->opposite()->prev();
	Halfedge* r = new_edge(h);
	make_sequence(h, r);
	make_sequence(r->opposite(), h->opposite());
	make_sequence(r, z1);
	make_sequence(z2, r->opposite());

	make_vertex_key(r, h->vertex());
	make_vertex_key(h, v);
	set_halfedge_vertex(r->opposite(), v);

	set_halfedge_facet(r, h->facet());
	set_halfedge_facet(r->opposite(), h->opposite()->facet());

	v->set_point(p);

	return v;
}


void MapEditor::erase_facet(Halfedge* h) {

	ogf_assert(!h->is_border());

	// Note: I have the feeling that this code should be
	//  much simpler ...

	std::vector<Halfedge*> edges_to_delete;

	delete_facet(h->facet());
	Halfedge* end = h;
	do {
		set_halfedge_facet(h, nil);
		Halfedge* g = h->next();
		bool h_opp_on_border = h->opposite()->is_border();
		bool g_opp_on_border = g->opposite()->is_border();
		if (h_opp_on_border) {
			// remove vertex X if it looks like that : 
			//
			//   ............
			//   ======X=====
			//
			// but not like that :
			//
			//     .....|
			//     .....|
			//     -----X-----
			//          |..... 
			//          |.....

			if (g_opp_on_border) {
				if (g->opposite()->next() == h->opposite()) {
					delete_vertex(h->vertex());
				}
				else {
					make_vertex_key(h->opposite()->prev());
					Halfedge* z1 = h->opposite()->prev();
					Halfedge* z2 = g->opposite()->next();
					make_sequence(z1, z2);
				}
			}
			edges_to_delete.push_back(h);
		}
		else {
			if (h->next()->opposite()->is_border()) {
				Halfedge* next =
					h->next()->opposite()->next();
				make_sequence(h, next);
				make_vertex_key(h);
			}
			if (h->prev()->opposite()->is_border()) {
				Halfedge* prev = h->prev()->opposite()->prev();
				make_sequence(prev, h);
				make_vertex_key(prev);
			}
		}
		h = g;
	} while (h != end);

	for (
		std::vector<Halfedge*>::iterator
		it = edges_to_delete.begin();
	it != edges_to_delete.end(); it++
		) {
		delete_edge(*it);
	}
}

//_____________________ utilities _____________________________

void MapEditor::remove_null_face(Halfedge* f0) {
	Halfedge* border = f0->next()->opposite();

	//remove facet
	delete_facet(f0->facet());

	// set link (fake set_opposite) to others
	make_sequence(f0, border->next());
	make_sequence(border->prev(), f0);
	set_halfedge_facet(f0, border->facet());

	// set links from others
	if (!f0->is_border()){
		make_facet_key(f0);
	}
	make_vertex_key(f0);
	make_vertex_key(f0->opposite());
	delete_edge(border);
}

bool MapEditor::can_merge_vertices(
	Halfedge* h0, Halfedge* h1
	) {
	// It's OK it they are already the same !
	if (h0->vertex() == h1->vertex()) {
		return true;
	}
	return (
		orbits_are_compatible(h0, h1) &&
		orbits_are_compatible(h1, h0)
		);
};

bool MapEditor::orbits_are_compatible(
	Halfedge* h0, Halfedge* h1
	) {

	Halfedge* cir_h0 = h0;
	do {
		// Number of potential opposites half_edges
		// (should not be greater than 1)
		int nb_common = 0;
		Halfedge* hh0 = cir_h0->opposite();
		Halfedge* cir_h1 = h1;
		do {
			Halfedge* hh1 = cir_h1->opposite();
			if (
				hh0->vertex() == hh1->vertex() ||
				(
				hh0->vertex() == h0->opposite()->vertex() &&
				hh1->vertex() == h1->opposite()->vertex()
				) || (
				hh0->vertex() == h1->opposite()->vertex() &&
				hh1->vertex() == h0->opposite()->vertex()
				)
				) {
				if (
                    (hh0->opposite()->is_border() && hh1->is_border()) ||
                    (hh0->is_border() && hh1->opposite()->is_border())
                    )
                {
					// Found a potential opposite edge.
					nb_common++;
				}
				else {
					// Potential opposite edge not on the border.
					return false;
				}
			}
			cir_h1 = cir_h1->next_around_vertex();
		} while (cir_h1 != h1);
		if (nb_common > 1) {
			return false;
		}
		cir_h0 = cir_h0->next_around_vertex();
	} while (cir_h0 != h0);
	return true;
}

bool MapEditor::halfedge_exists_between_vertices(Vertex* v1, Vertex* v2) {
	Halfedge* cir = v1->halfedge();
	do {
		if (cir->opposite()->vertex() == v2) {
			return true;
		}
		cir = cir->next_around_vertex();
	} while (cir != v1->halfedge());
	return false;
}


bool MapEditor::halfedges_on_same_vertex(Halfedge* h1, Halfedge* h2) {
	Halfedge* it = h1;
	do {
		if (it == h2) {
			return true;
		}
		it = it->next_around_vertex();
	} while (it != h1);
	return false;
}

bool MapEditor::halfedges_on_same_facet(Halfedge* h1, Halfedge* h2) {
	if (h1->facet() == nil || h2->facet() == nil)
		return false;

	Halfedge* it = h1;
	do {
		if (it == h2) {
			return true;
		}
		it = it->next();
	} while (it != h1);
	return false;
}

bool MapEditor::can_collapse_edge(Halfedge* h) {

	{// disable glueing problems

		if ((!h->is_border() &&
			h->next()->opposite()->facet()
			== h->prev()->opposite()->facet())
			|| (!h->opposite()->is_border() &&
			h->opposite()->next()->opposite()->facet()
			== h->opposite()->prev()->opposite()->facet())
			)
			return false;
	}
	{// don't remove more than one vertex
		if (// it's a triangle
			h->next()->next()->next() == h
			// the vertex is alone on border
			&& h->next()->opposite()->is_border()
			&& h->prev()->opposite()->is_border()
			)
			return false;

		// the same on the other face
		if (// it's a triangle
			h->opposite()->next()->next()->next() == h
			// the vertex is alone on border
			&& h->opposite()->next()->opposite()->is_border()
			&& h->opposite()->prev()->opposite()->is_border()
			)
			return false;
	}

	// don't do stupid holes
	{
		if (
			(h->is_border() && h->next()->next()->next() == h) ||
			(
			h->opposite()->is_border() &&
			h->opposite()->next()->next()->next() == h->opposite()
			)
			) {
			return false;
		}
	}

	// don't merge holes (i.e. don't split surface)
	{
		if (
			!h->is_border() &&
			!h->opposite()->is_border() &&
			h->vertex()->is_on_border() &&
			h->opposite()->vertex()->is_on_border()
			) {
			return false;
		}
	}

	// be carefull of the toblerone case (don't remove volumes)
	{
		Halfedge* cir = h;
		int nb_twice = 0;

		std::set<Vertex*> marked;

		// do { 
		//    cir->opposite()->vertex()->set_id(0);
		//    cir = cir->next_around_vertex();
		// } while ( cir != h);

		cir = h->opposite();
		do {
			marked.insert(cir->opposite()->vertex());
			cir = cir->next_around_vertex();
		} while (cir != h->opposite());

		cir = h;
		do {
			if (
				marked.find(cir->opposite()->vertex()) != marked.end()
				) {
				nb_twice++;
			}
			marked.insert(cir->opposite()->vertex());
			cir = cir->next_around_vertex();
		} while (cir != h);

		if (h->next()->next()->next() == h)  {
			nb_twice--;
		}
		if (h->opposite()->next()->next()->next() == h->opposite()) {
			nb_twice--;
		}

		if (nb_twice > 0) {
			//std::cerr<<" \nbe carefull of the toblerone case";
			return false;
		}
	}
	return true;


	/*

	if(
	!h->next()->opposite()->is_border() ||
	!h->prev()->opposite()->is_border()
	) {
	return false ;
	}
	if(
	!h->opposite()->next()->opposite()->is_border() ||
	!h->opposite()->prev()->opposite()->is_border()
	) {
	return false ;
	}
	return true ;
	*/
}

bool MapEditor::collapse_edge(Halfedge* h) {

	if (!can_collapse_edge(h)) {
		return false;
	}

	Vertex* dest = h->vertex();

	// everyone has the same vertex
	{
		Vertex* v = h->opposite()->vertex();
		Halfedge* i;
		Halfedge* ref;
		i = ref = h->opposite();
		do {
			Halfedge* local = i;
			set_halfedge_vertex(local, dest);
			i = i->next_around_vertex();
		} while (i != ref);
		delete_vertex(v);
	}


	// remove links to current edge (facet & vertex)
	Halfedge* hDle = h;
	if (!h->is_border()) {
		set_facet_halfedge(hDle->facet(), hDle->next());
	}

	if (!h->opposite()->is_border()) {
		set_facet_halfedge(
			hDle->opposite()->facet(), hDle->opposite()->next()
			);
	}
	set_vertex_halfedge(dest, hDle->next()->opposite());

	Halfedge* f0 = h->next();
	Halfedge* f1 = h->opposite()->prev();

	// update prev/next links
	{
		Halfedge* e = h->next();
		make_sequence(hDle->prev(), e);
		e = hDle->opposite()->next();
		make_sequence(hDle->opposite()->prev(), e);
	}

	// remove null faces
	if (f0->next()->next() == f0) {
		remove_null_face(f0);
	}

	if (f1->next()->next() == f1) {
		remove_null_face(f1);
	}

	delete_edge(hDle);

	return true;
}

void MapEditor::reorient_facet(Map::Halfedge* first) {
    if (first == nil) {
        return;
    }
    Map::Halfedge* last = first;
    Map::Halfedge* prev = first;
    Map::Halfedge* start = first;
    first = first->next();
    Map::Vertex* new_v = start->vertex();
    while (first != last) {
        Map::Vertex*  tmp_v = first->vertex();
        set_halfedge_vertex(first, new_v);
        set_vertex_halfedge(first->vertex(), first);
        new_v = tmp_v;
        Map::Halfedge* next = first->next();
        set_halfedge_next(first, prev);
        set_halfedge_prev(first, next);
        prev = first;
        first = next;
    }
    set_halfedge_vertex(start, new_v);
    set_vertex_halfedge(start->vertex(), start);
    Map::Halfedge* next = start->next();
    set_halfedge_next(start, prev);
    set_halfedge_prev(start, next);
}



bool MapEditor::glue(
    Halfedge* h0, Halfedge* h1
    )
{
    if (!can_glue(h0, h1)) {
        return false;
    }

    vec3 new_p0 = Geom::barycenter(
        h0->vertex()->point(), h1->opposite()->vertex()->point()
        );

    vec3 new_p1 = Geom::barycenter(
        h1->vertex()->point(), h0->opposite()->vertex()->point()
        );


    // merge vertices if necessary

    Vertex* dest0 = h0->vertex();
    Vertex* dest1 = h1->vertex();

    Vertex* org0 = h0->opposite()->vertex();
    Vertex* org1 = h1->opposite()->vertex();

    if (is_locked_[dest1]) {
        is_locked_[org0] = true;
    }

    if (is_locked_[dest0]) {
        is_locked_[org1] = true;
    }

    if (org0 != dest1) {
        set_vertex_on_orbit(h1, org0);
        delete_vertex(dest1);
    }

    if (org1 != dest0) {
        set_vertex_on_orbit(h0, org1);
        delete_vertex(dest0);
    }

    // set halfedge connections

    make_sequence(h1->prev(), h0->next());
    make_sequence(h0->prev(), h1->next());
    make_opposite(h0->opposite(), h1->opposite());
    make_vertex_key(h0->opposite());
    make_vertex_key(h1->opposite());

    org1->set_point(new_p0);
    org0->set_point(new_p1);

    delete_halfedge(h0);
    delete_halfedge(h1);

    return true;
}

bool MapEditor::can_glue(Halfedge* h0, Halfedge* h1) {

    // Checks that both Halfedges are on the border.
    if (!h0->is_border() || !h1->is_border()) {
        return false;
    }

    // don't glue two halfedges on a same face
    if (
        h0->opposite()->facet() == h1->opposite()->facet()
        ) {
        return false;
    }

    // don't merge two vertices on a same halfedge
    if (
        halfedge_exists_between_vertices(
        h0->vertex(), h1->opposite()->vertex()
        ) ||
        halfedge_exists_between_vertices(
        h1->vertex(), h0->opposite()->vertex()
        )
        ) {
        return false;
    }

    if (
        !can_merge_vertices(h0, h1->opposite()) ||
        !can_merge_vertices(h1, h0->opposite())
        ) {
        return false;
    }

    return true;
}


void MapEditor::copy_attributes(Vertex* to, Vertex* from) {
	target()->vertex_attribute_manager()->copy_record(to, from);
}

void MapEditor::copy_attributes(Halfedge* to, Halfedge* from) {
	target()->halfedge_attribute_manager()->copy_record(to, from);
}

void MapEditor::copy_attributes(Facet* to, Facet* from) {
	target()->facet_attribute_manager()->copy_record(to, from);
}
