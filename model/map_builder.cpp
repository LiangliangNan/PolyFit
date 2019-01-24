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


#include "map_builder.h"
#include "map.h"
#include "../basic/logger.h"

#include <algorithm>


struct FacetKey {
	FacetKey(Map::Vertex* a1, Map::Vertex* a2, Map::Vertex* a3) {
		v1 = a1 ; v2 = a2 ; v3 = a3 ;
	}
	bool operator<(const FacetKey& rhs) const {
		if(v1 < rhs.v1) { return true ; }
		if(v1 > rhs.v1) { return false ; }
		if(v2 < rhs.v2) { return true ; }
		if(v2 > rhs.v2) { return false ; }
		return (v3 < rhs.v3) ; 
	}
	Map::Vertex* v1 ;
	Map::Vertex* v2 ;
	Map::Vertex* v3 ;
} ;

static std::set<FacetKey>* all_facet_keys ;

//_________________________________________________________

AbstractMapBuilder::~AbstractMapBuilder() { }

void AbstractMapBuilder::create_vertices(unsigned int nb_vertices, bool with_colors) {
	for(unsigned int i=0; i<nb_vertices; i++) {
		add_vertex(vec3(0.0, 0.0, 0.0)) ;
	}
}

void AbstractMapBuilder::reserve_vertices(unsigned int nb_vertices) {
	// Do nothing by default
}

void AbstractMapBuilder::reserve_facets(unsigned int nb_facets) {
	// Do nothing by default
}


//_________________________________________________________

/*
Map* MapBuilder::map() {
return target() ;
}
*/
void MapBuilder::begin_surface() {
	transition(initial, surface) ;
	star_.bind(target()) ;
	is_locked_.bind(target()) ;
	nb_non_manifold_v_ = 0 ;
	nb_duplicate_e_    = 0 ;
	nb_isolated_v_     = 0 ;
	all_facet_keys = new std::set<FacetKey> ;
}

void MapBuilder::end_surface() {
	delete all_facet_keys ; 
	all_facet_keys = nil ;

	transition(surface, final) ;
	terminate_surface() ;
	vertex_.clear() ;
	star_.unbind() ;
	if(color_.is_bound()) {
		color_.unbind() ;
	}
	if(!quiet_ && (nb_non_manifold_v_ != 0)) {
		Logger::warn("MapBuilder") 
			<< "Encountered " << nb_non_manifold_v_
			<< " non-manifold vertices (fixed)" << std::endl ;
	}
	if(!quiet_ && (nb_duplicate_e_ != 0)) {
		Logger::warn("MapBuilder") 
			<< "Encountered " << nb_duplicate_e_
			<< " duplicated edges (fixed)" << std::endl ;
	}
	if(!quiet_ && (nb_isolated_v_ != 0)) {
		Logger::warn("MapBuilder") 
			<< "Encountered " << nb_isolated_v_
			<< " isolated vertices (removed)" << std::endl ;
	}
	if (!quiet_ && target()->size_of_facets() == 0) {
		Logger::warn("MapBuilder") 
			<< "0 facets. Is it a point cloud?" << std::endl ;
	}
}

void MapBuilder::reset() {
	transition(final, initial) ;
}

void MapBuilder::add_vertex(unsigned int id, const vec3& p) {
	transition(surface, surface) ;
	add_vertex_internal(id, p) ;
}

void MapBuilder::add_vertex(const vec3& p) {
	transition(surface, surface) ;
	add_vertex_internal(p) ;
}

void MapBuilder::add_vertex_internal(unsigned int id, const vec3& p) {
	Map::Vertex* new_v = new_vertex() ;
	new_v->set_point(p) ;
	while(vertex_.size() <= id) {
		vertex_.push_back(nil) ;
	}
	vertex_[id] = new_v ;
}

void MapBuilder::add_vertex_internal(const vec3& p) {
	add_vertex_internal(vertex_.size(), p) ;
}

void MapBuilder::begin_facet() {
	facet_vertex_.clear() ;
}

void MapBuilder::end_facet() {
	int nb_vertices = facet_vertex_.size() ;

	if(nb_vertices < 3) {
		if(!quiet_) {
			Logger::err("MapBuilder")
				<< "Facet with less than 3 vertices, ignored"
				<< std::endl ;
		}
		return ;
	}

	if(false && nb_vertices == 3) {
		std::vector<Map::Vertex*> W = facet_vertex_ ;
		std::sort(W.begin(), W.end()) ;
		FacetKey k(W[0], W[1],W[2]) ;
		if(all_facet_keys->find(k) != all_facet_keys->end()) {
			std::cerr << '.' << std::flush ;
			return ;
		}
		all_facet_keys->insert(k) ;
	}


	// Detect duplicated vertices
	{for(int i=0; i<nb_vertices; i++) {
		for(int j=i+1; j<nb_vertices; j++) {
			if(facet_vertex_[i] == facet_vertex_[j]) {
				if(!quiet_) {
					Logger::err("MapBuilder") 
						<< "Facet with duplicated vertices, ignored"
						<< std::endl ;
				}
				return ;
			}
		}
	}}

	// Detect and remove non-manifold edges by duplicating
	//   the corresponding vertices.
	{for(int from=0; from<nb_vertices; from++) {
		int to=((from+1) % nb_vertices) ;
		if(find_halfedge_between(facet_vertex_[from], facet_vertex_[to]) != nil) 
		{
			nb_duplicate_e_++ ;
			facet_vertex_[from] = copy_vertex(facet_vertex_[from]) ;
			facet_vertex_[to] = copy_vertex(facet_vertex_[to]) ;
		}
	}}


	begin_facet_internal() ;
	for(int i=0; i<nb_vertices; i++) {
		add_vertex_to_facet_internal(facet_vertex_[i]) ;
	}
	end_facet_internal() ;
}


Map::Vertex* MapBuilder::copy_vertex(Vertex* from) {
	// Note: tex coords are not copied, since 
	//  an halfedge does not exist in the copy to
	//  store the texvertex pointer.

	Vertex* result = new_vertex() ;
	result->set_point(from->point()) ;
	if(is_locked_[from]) {
		is_locked_[result] = true ;
	}
	return result ;
}

void MapBuilder::add_vertex_to_facet(int i) {
	if(i < 0 || i >= int(vertex_.size())) {
		if(!quiet_) {
			Logger::err("MapBuilder") << "vertex index " 
				<< i << " out of range"
				<< std::endl ;
		}
		return ;
	}
	facet_vertex_.push_back(vertex_[i]) ;
}

void MapBuilder::lock_vertex(int i) {
	if(i < 0 || i >= int(vertex_.size())) {
		if(!quiet_) {
			Logger::err("MapBuilder") << "vertex index " 
				<< i << " out of range"
				<< std::endl ;
		}
		return ;
	}
	is_locked_[vertex(i)] = true ;
}

void MapBuilder::unlock_vertex(int i) {
	if(i < 0 || i >= int(vertex_.size())) {
		if(!quiet_) {
			Logger::err("MapBuilder") << "vertex index " 
				<< i << " out of range"
				<< std::endl ;
		}
		return ;
	}
	is_locked_[vertex(i)] = false ;
}

void MapBuilder::begin_facet_internal() {
	transition(surface, facet) ;
	current_facet_ = new_facet() ;
	first_vertex_in_facet_ = nil ;
	current_vertex_ = nil ;
	first_halfedge_in_facet_ = nil ;
	current_halfedge_ = nil ;
}

void MapBuilder::end_facet_internal() {
	transition(facet, surface) ;
	Map::Halfedge* h = 
		new_halfedge_between(current_vertex_, first_vertex_in_facet_) ;

	make_sequence(current_halfedge_, h) ;
	make_sequence(h, first_halfedge_in_facet_) ;
}

void MapBuilder::add_vertex_to_facet_internal(Vertex* v) {
	transition(facet, facet) ;

	if(first_vertex_in_facet_ == nil) {
		first_vertex_in_facet_ = v ;
	} else {
		Map::Halfedge* new_halfedge = 
			new_halfedge_between(current_vertex_, v) ;

		if(first_halfedge_in_facet_ == nil) {
			first_halfedge_in_facet_ = new_halfedge ;
			make_facet_key(first_halfedge_in_facet_) ;
		} else {
			make_sequence(current_halfedge_, new_halfedge) ;
		}
		current_halfedge_ = new_halfedge ;
	}
	current_vertex_ = v ;
}

Map::Vertex* MapBuilder::current_vertex() {
	return vertex_[vertex_.size() - 1] ;
}

Map::Vertex* MapBuilder::vertex(int i) {
	return vertex_[i] ;
}

Map::Facet* MapBuilder::current_facet() {
	return current_facet_ ;
}

Map::Halfedge* MapBuilder::new_halfedge_between(
	Vertex* from, Vertex* to
	) {

		ogf_assert(from != to) ;

		// Non-manifold edges have been removed
		// by the higher level public functions.
		// Let us do a sanity check anyway ...
		ogf_assert(find_halfedge_between(from, to) == nil) ;

		Map::Halfedge* result = new_halfedge() ;
		set_halfedge_facet(result, current_facet_) ;
		set_halfedge_vertex(result, to) ;

		Map::Halfedge* opposite = find_halfedge_between(to, from) ;
		if(opposite != nil) {
			make_opposite(result, opposite) ;
		}

		star_[from].push_back(result) ;
		set_vertex_halfedge(to, result) ;

		return result ;
}

Map::Halfedge* MapBuilder::find_halfedge_between(
	Vertex* from, Vertex* to
	) {
		Star& star = star_[from] ;
		for(Star::iterator it = star.begin(); it != star.end(); it++) {
			Map::Halfedge* cur = *it ;
			if(cur->vertex() == to) {
				return cur ;
			}
		}
		return nil ;
}

bool MapBuilder::vertex_is_manifold(Vertex* v) {

	if(v->halfedge() == nil) {
		// Isolated vertex, removed in subsequent steps
		return true ;
	}

	// A vertex is manifold if the stored and the 
	// computed stars match (the number of halfedges
	// are tested).
	// Note: this test is valid only if the borders
	// have been constructed.

	return (int(star_[v].size()) == v->degree()) ;
}

bool MapBuilder::split_non_manifold_vertex(Vertex* v) {
	if(vertex_is_manifold(v)) {
		return false ;
	}

	std::set<Map::Halfedge*> star ;
	{ for(unsigned int i=0; i<star_[v].size(); i++) {
		star.insert(star_[v][i]) ;
	}}

	// For the first wedge, reuse the vertex
	disconnect_vertex(
		v->halfedge()->opposite(), v, star
		) ;

	// There should be other wedges (else the vertex
	// would have been manifold)
	ogf_assert(!star.empty()) ;

	{
		MapVertexLock is_locked(target()) ;
		is_locked[v] = true ;
		// Create the vertices for the remaining wedges.
		while(!star.empty()) {
			Vertex* new_v = copy_vertex(v) ;
			is_locked[new_v] = true ;
			Map::Halfedge* h = *(star.begin()) ;
			disconnect_vertex(h, new_v, star) ;
		}
	}

	return true ;
}


void MapBuilder::disconnect_vertex(
								   Map::Halfedge* start_in, Map::Vertex* v, std::set<Map::Halfedge*>& star
								   ) 
{

	Map::Halfedge* start = start_in ;

	star_[v].clear() ;

	//   Important note: in this class, all the Stars correspond to the
	// set of halfedges radiating FROM a vertex (i.e. h->vertex() != v
	// if h belongs to Star(v) ).

	ogf_assert(star.find(start) != star.end()) ;


	//   Note: the following code needs a special handling
	// of borders, since borders are not correctly connected
	// in a non-manifold vertex (the standard next_around_vertex
	// iteration does not traverse the whole star)

	while(!start->is_border()) {
		start = start->prev()->opposite() ;
		if(start == start_in) {
			break ;
		}
	}
	set_vertex_halfedge(v,start->opposite()) ;

	Map::Halfedge* cur = start ;
	set_halfedge_vertex(cur->opposite(), v) ;
	star_[v].push_back(cur) ;
	std::set<Map::Halfedge*>::iterator it = star.find(cur) ;
	ogf_assert(it != star.end()) ;
	star.erase(it) ;

	while(!cur->opposite()->is_border()) {
		cur = cur->opposite()->next() ;
		if(cur == start) {
			break ;
		}
		set_halfedge_vertex(cur->opposite(), v) ;
		std::set<Map::Halfedge*>::iterator it = star.find(cur) ;
		ogf_assert(it != star.end()) ;
		star.erase(it) ;
		star_[v].push_back(cur) ;
	} 

	if(start->is_border()) {
		make_sequence(cur->opposite(),start) ;
	}
}


void MapBuilder::terminate_surface() {

	// Step 1 : create the border halfedges, and setup the 'opposite'
	//   and 'vertex' pointers.
	FOR_EACH_HALFEDGE(Map, target(), it) {
		if(it->opposite() == nil) {
			Map::Halfedge* h = new_halfedge() ;
			make_opposite(h, it) ;
			set_halfedge_vertex(h, it->prev()->vertex()) ;

			// Used later to fix non-manifold vertices.
			star_[ it->vertex() ].push_back(h) ;
		}
	}

	// Step 2 : setup the 'next' and 'prev' pointers of the border.
	//
	FOR_EACH_HALFEDGE(Map, target(), it) {
		if(it->facet() == nil) {

			Map::Halfedge* next = it->opposite() ;
			while(next->facet() != nil) {
				next = next->prev()->opposite() ;
			}
			set_halfedge_next(it, next) ;

			Map::Halfedge* prev = it->opposite() ;
			while(prev->facet() != nil) {
				prev = prev->next()->opposite() ;
			}
			set_halfedge_prev(it, prev) ;
		}
	}

	// Step 3 : fix non-manifold vertices

	for(unsigned int i=0; i<vertex_.size(); i++) {
		if(vertex_[i] == nil) {
			continue ;
		}
		if(split_non_manifold_vertex(vertex_[i])) {
			nb_non_manifold_v_++ ;
		}
	}

	// Step 4 : check for isolated vertices
	for(unsigned int i=0; i<vertex_.size(); i++) {
		if(vertex_[i] == nil) {
			continue ;
		}
		if(star_[vertex_[i]].size() == 0) {
			nb_isolated_v_++ ;
			delete_vertex(vertex_[i]) ;
		}
	}

	target()->assert_is_valid() ;
} 

void MapBuilder::transition(state from, state to) {
	check_state(from) ;
	state_ = to ;
}

void MapBuilder::check_state(state s) {
	if(state_ != s) {
		Logger::err("MapBuilder") << "invalid state: \'" 
			<< state_to_string(state_)
			<< "\', expected \'"
			<< state_to_string(s) 
			<< "\'"
			<< std::endl ;
		bool ok = false ;
		ogf_assert(ok) ;
	}
}

std::string MapBuilder::state_to_string(state s) {
	switch(s) {
		case initial:
			return "initial" ;
		case surface:
			return "surface" ;
		case facet: 
			return "facet" ;
		case final:
			return "final" ;
		default:
			return "unknown" ;
	}
	return "unknown" ;
}

void MapBuilder::set_vertex(unsigned int v, const vec3& p) {
	vertex(v)->set_point(p) ;
}

void MapBuilder::set_vertex_color(unsigned int v, const Color& c) {
	if(!color_.is_bound()) {
		color_.bind(target(), "color") ;
	}
	color_[vertex(v)] = c ;
}

void MapBuilder::create_vertices(unsigned int nb_v, bool with_colors) {
	AbstractMapBuilder::create_vertices(nb_v) ;
	if(with_colors && !color_.is_bound()) {
		color_.bind(target(), "color") ;
	}
}
