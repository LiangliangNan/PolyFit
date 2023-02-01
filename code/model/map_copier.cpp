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

#include "map_copier.h"
#include "map_builder.h"



MapCopier::MapCopier() { 
	copy_all_attributes_ = false ;
	current_source_ = nil ;
	current_destination_ = nil ;
	// By default: copy vertex locks
	declare_vertex_attribute_to_copy("lock") ;
}

MapCopier::~MapCopier() { 
}        

void MapCopier::copy(
					 MapBuilder& builder, Map* source,
					 MapVertexAttribute<int>& vertex_id, 
					 int& cur_vertex_id
					 )
{
	bind_attribute_copiers(builder.target(), source) ;

	// Step 1 : clear vertex ids
	FOR_EACH_VERTEX(Map, source, it) {
		vertex_id[it] = -1 ;
	}

	// Step 2: enumerate vertices
	FOR_EACH_VERTEX(Map, source, it) {
		vertex_id[it] = cur_vertex_id ;
		builder.add_vertex(it->point()) ;
		copy_vertex_attributes(builder.current_vertex(), it) ;
		cur_vertex_id++ ;
	}

	// Step 3: create facets
	FOR_EACH_FACET(Map, source, it) {
		Map::Halfedge* h = it->halfedge() ;
		builder.begin_facet() ;
		do {
			builder.add_vertex_to_facet(vertex_id[h->vertex()]) ;
			h = h->next() ;
		} while(h != it->halfedge()) ;
		builder.end_facet() ;
		copy_facet_attributes(builder.current_facet(), it) ;
		// TODO: copy halfedge attributes
	}
}

template <class RECORD> inline void bind_attribute_copiers(
	std::vector< AttributeCopier<RECORD> >& copiers,
	AttributeManager* to, AttributeManager* from,
	const std::set<std::string>& attributes_to_copy,
	bool copy_all_attributes
	) {
		copiers.clear() ;
		std::vector<std::string> names ;
		from->list_named_attributes(names) ;
		for(unsigned int i=0; i<names.size(); i++) {
			if(copy_all_attributes || (attributes_to_copy.find(names[i]) != attributes_to_copy.end())) {
				bind_source(copiers, from, names[i]) ;
			}
		}
		bind_destinations(copiers, to) ;
}

void MapCopier::bind_attribute_copiers(Map* destination, Map* source) {
	::bind_attribute_copiers(
		vertex_attribute_copiers_,
		destination->vertex_attribute_manager(), source->vertex_attribute_manager(), 
		vertex_attributes_to_copy_,
		copy_all_attributes_
		) ;

	::bind_attribute_copiers(
		halfedge_attribute_copiers_,
		destination->halfedge_attribute_manager(), source->halfedge_attribute_manager(), 
		halfedge_attributes_to_copy_,
		copy_all_attributes_
		) ;

	::bind_attribute_copiers(
		facet_attribute_copiers_,
		destination->facet_attribute_manager(), source->facet_attribute_manager(), 
		facet_attributes_to_copy_,
		copy_all_attributes_
		) ;

}

