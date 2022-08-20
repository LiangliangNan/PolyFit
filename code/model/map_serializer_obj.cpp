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


#include "map_serializer_obj.h"
#include "../basic/logger.h"
#include "../basic/file_utils.h"
#include "../basic/line_stream.h"
#include "../model/map_builder.h"
#include "../model/map_enumerator.h"
#include "../basic/generic_attributes_io.h"


MapSerializer_obj::MapSerializer_obj() {
	read_supported_ = true ;
	write_supported_ = true ;
}

bool MapSerializer_obj::serialize_read(
									   const std::string& file_name, Map* mesh
									   ) 
{
	current_directory_ = FileUtils::dir_name(file_name) ;
	bool flag = MapSerializer::serialize_read(file_name, mesh) ;
	current_directory_ = "" ;
	return flag;
}


void MapSerializer_obj::read_mtl_lib(std::istream& input) {
	LineInputStream in(input) ;
	std::string keyword ;
	std::string cur_mtl = "default" ;
	while(!in.eof()) {
		in.get_line() ;
		in >> keyword ;
		if(keyword == "newmtl") {
			in >> cur_mtl ;
		} else if(keyword == "Kd") {
			Color c ;
			in >> c ;
			material_lib_[cur_mtl] = c ;
		}
	}        
}

bool MapSerializer_obj::do_read(std::istream& input, AbstractMapBuilder& builder) {
	MapFacetAttribute<Color> color_ ;
	LineInputStream in(input) ;
	MapBuilder* concrete_builder = dynamic_cast<MapBuilder*>(&builder) ;

	builder.begin_surface() ;
	while(!in.eof()) {
		in.get_line() ;

		std::string keyword ;
		in >> keyword ;

		if(keyword == "v") {
			vec3 p ;
			in >> p ;
			builder.add_vertex(p) ;
		} else if(keyword == "vt") {
			vec2 q ;
			in >> q ;
			//builder.add_tex_vertex(q) ;
		} else if(keyword == "f") {
			builder.begin_facet() ;
			while(!in.eol()) {
				std::string s ;
				in >> s ;
				if(s.length() > 0) {
					std::istringstream v_input(s) ;
					int index ;
					v_input >> index ;
					builder.add_vertex_to_facet(index - 1) ;
					char c ;
					v_input >> c ;
					if(c == '/') {
						v_input >> index ;
						//builder.set_corner_tex_vertex(index - 1) ;
					}
				}
			}
			builder.end_facet() ;
			if(color_.is_bound() && concrete_builder != nil) {
				color_[concrete_builder->current_facet()] = current_material_ ;
			}
		} else if(keyword == "#") {
			std::string second_keyword ;
			in >> second_keyword ;
			if(second_keyword == "anchor") {
				int index ;
				in >> index ;
				builder.lock_vertex(index - 1) ;
			} 
		} else if(keyword == "mtllib") {
			std::string mtl_lib_filename ;
			in >> mtl_lib_filename ;
			mtl_lib_filename = current_directory_ + "/" + mtl_lib_filename ;
			std::ifstream mtl_lib_in(mtl_lib_filename.c_str()) ;
			if(mtl_lib_in) {
				Logger::out("MapSerializer_obj") << "using material lib: " << mtl_lib_filename << std::endl ;
				read_mtl_lib(mtl_lib_in) ;
			} else {
				Logger::err("MapSerializer_obj") << mtl_lib_filename << ": no such file" << std::endl ;
			}
		} else if(keyword == "usemtl") {
			std::string material ;
			in >> material ;
			MaterialLib::iterator it = material_lib_.find(material) ;
			if(it == material_lib_.end()) {
				current_material_ = Color(0.7f, 0.7f, 0.7f, 1.0f) ;
			} else {
				current_material_ = it->second ;
			}
			if(!color_.is_bound() && concrete_builder != nil) {
				color_.bind(concrete_builder->target(), "color") ;
			}
		}
	}

	builder.end_surface() ;
	color_.unbind() ;
	material_lib_.clear() ;
	return true ;
}

bool MapSerializer_obj::do_write(std::ostream& out, const Map* mesh) const {
	// Obj files numbering starts with 1
	Attribute<Vertex, int>	vertex_id(mesh->vertex_attribute_manager());
	MapEnumerator::enumerate_vertices(const_cast<Map*>(mesh), vertex_id, 1);

	// Output Vertices
	FOR_EACH_VERTEX_CONST(Map, mesh, it) {
		out << "v " << it->point() << std::endl ;
	} 

	// Output facets
	FOR_EACH_FACET_CONST(Map, mesh, it) {
		Map::Halfedge* jt = it->halfedge() ;
		out << "f " ;
		do {
			out << vertex_id[jt->vertex()] << " ";
			jt = jt->next() ;
		} while(jt != it->halfedge()) ;
		out << std::endl ;
	}

	MapVertexLock is_locked(const_cast<Map*>(mesh));
	FOR_EACH_VERTEX_CONST(Map, mesh, it) {
		if(is_locked[it]) {
			out << "# anchor " << vertex_id[it] << std::endl ;
		}
	}

	return true ;
}

//_________________________________________________________

MapSerializer_eobj::MapSerializer_eobj() : MapSerializer_obj() {
}

static Map::Halfedge* find_halfedge_between(Map::Vertex* v1, Map::Vertex* v2) {
	Map::Halfedge* h = v2->halfedge() ;
	do {
		if(h->opposite()->vertex() == v1) {
			return h ;
		}
		h = h->next_around_vertex() ;
	} while(h != v2->halfedge()) ;
	return nil ;
}


bool MapSerializer_eobj::do_read(std::istream& input, AbstractMapBuilder& out) 
{
	MapBuilder& builder = dynamic_cast<MapBuilder&>(out) ;

	Map* map = builder.target() ;

	std::vector< SerializedAttribute<Map::Vertex>   > v_attributes ;
	std::vector< SerializedAttribute<Map::Halfedge> > h_attributes ;
	std::vector< SerializedAttribute<Map::Facet>    > f_attributes ;

	std::vector<Map::Facet*> facets ;

	bool surface_terminated = false ;

	LineInputStream in(input) ;

	builder.begin_surface() ;
	while(!in.eof()) {
		in.get_line() ;
		std::string keyword ;

		in >> keyword ;

		if(keyword == "]]>") {
			break ; // end of CDATA in gsg file.
		} else if(keyword == "v") {
			vec3 p ;
			in >> p ;
			builder.add_vertex(p) ;
		} else if(keyword == "vt") {
			vec2 q ;
			in >> q ;
			//builder.add_tex_vertex(q) ;
		} else if(keyword == "f") {
			builder.begin_facet() ;
			while(!in.eol()) {
				std::string s ;
				in >> s ;
				if(s.length() > 0) {
					std::istringstream v_input(s) ;
					int index ;
					v_input >> index ;
					builder.add_vertex_to_facet(index - 1) ;
					char c ;
					v_input >> c ;
					if(c == '/') {
						v_input >> index ;
						//builder.set_corner_tex_vertex(index - 1) ;
					}
				}
			}
			builder.end_facet() ;
			facets.push_back(builder.current_facet()) ;
		} else if(keyword == "#") {
			std::string second_keyword ;
			in >> second_keyword ;

			if(second_keyword == "attribute") {
				if(!surface_terminated) {
					// Quick and dirty fix to ensure that
					// border edges exist before we put
					// attributes on them !!!
					builder.terminate_surface() ;
					surface_terminated = true ;
				}
				std::string name ;
				std::string location ;
				std::string type ;
				in >> name >> location >> type ;
				std::cerr << "Attribute " << name << " on " << location << " : " << type << std::endl ;
				if(location == "vertex") {
					v_attributes.push_back(SerializedAttribute<Map::Vertex>()) ;
					v_attributes.rbegin()->bind(map->vertex_attribute_manager(), name, type) ;
				} else if(location == "halfedge") {
					h_attributes.push_back(SerializedAttribute<Map::Halfedge>()) ;
					h_attributes.rbegin()->bind(map->halfedge_attribute_manager(), name, type) ;
				} else if(location == "facet") {
					f_attributes.push_back(SerializedAttribute<Map::Facet>()) ;
					f_attributes.rbegin()->bind(map->facet_attribute_manager(), name, type) ;
				} else {
					Logger::warn("MapSerializer_eobj") << "Invalid attribute location:" 
						<< location << std::endl ; 
				}
			} else if(second_keyword == "attrs") {

				std::string location ;
				in >> location ;
				if(location == "v") {
					int id ;
					in >> id ;
					id-- ;
					Map::Vertex* v = builder.vertex(id) ;
					serialize_read_attributes(in.line(), v, v_attributes) ;
				} else if(location == "h") {
					int id1, id2 ;
					in >> id1 >> id2 ;
					id1-- ; id2-- ;
					Map::Vertex* v1 = builder.vertex(id1) ;
					Map::Vertex* v2 = builder.vertex(id2) ;
					Map::Halfedge* h = find_halfedge_between(v1,v2) ;
					if(h == nil) {
						Logger::warn("MapSerializer_eobj") << "Halfedge does not exist" << std::endl ;
					} else {
						serialize_read_attributes(in.line(), h, h_attributes) ;
					}
				} else if(location == "f") {
					int id ;
					in >> id ;
					id-- ;
					Map::Facet* f = facets[id] ;
					serialize_read_attributes(in.line(), f, f_attributes) ;
				} else {
					Logger::warn("MapSerializer_eobj") << "Invalid attribute location:" 
						<< location << std::endl ; 
				}
			} else if(second_keyword == "anchor") {
				int index ;
				in >> index ;
				builder.lock_vertex(index - 1) ;
			} else if(second_keyword == "END") {
				break ;
			}
		} 
	}

	builder.end_surface() ;

	return true ;
}



// Note: we cannot use MapSerializer_obj::serialize_write,
// since it saves anchors using special keywords. In our
// case, anchors are handled by the generic mechanism.

bool MapSerializer_eobj::do_write(std::ostream& output, const Map* mesh) const {
	Attribute<Vertex,int> vertex_id(
		mesh->vertex_attribute_manager()
		) ;

	// Obj files numbering starts with 1 (instead of 0)
	MapEnumerator::enumerate_vertices(const_cast<Map*>(mesh), vertex_id, 1) ;

	// Output Vertices
	{ FOR_EACH_VERTEX_CONST(Map, mesh, it) {
		output << "v " << it->point() << std::endl ;
	}} 

	// Output facets

	{ FOR_EACH_FACET_CONST(Map, mesh, it) {
		Halfedge* jt = it->halfedge() ;
		output << "f " ;
		do {
			output << vertex_id[jt-> vertex()] << " " ;
			jt = jt->next() ;
		} while(jt != it->halfedge()) ;
		output << std::endl ;
	}}


	{
		std::vector<SerializedAttribute<Map::Vertex> > attributes ;
		if(get_serializable_attributes(mesh->vertex_attribute_manager(), attributes, output, "vertex")) {
			int vid = 1 ;
			FOR_EACH_VERTEX_CONST(Map, mesh, it) {
				const Map::Vertex* v = it ;
				output << "# attrs v " << vid << " " ;
				serialize_write_attributes(output, v, attributes) ;
				output << std::endl ;
				vid++ ;
			}
		}
	}

	{
		std::vector<SerializedAttribute<Map::Halfedge> > attributes ;
		if(get_serializable_attributes(mesh->halfedge_attribute_manager(), attributes, output, "halfedge")) {
			FOR_EACH_HALFEDGE_CONST(Map, mesh, it) {
				const Map::Halfedge* h = it ;
				output << "# attrs h " 
					<< vertex_id[h->opposite()->vertex()] << " " << vertex_id[h->vertex()] << " " ;
				serialize_write_attributes(output, h, attributes) ;
				output << std::endl ;
			}
		}
	}

	{
		std::vector<SerializedAttribute<Map::Facet> > attributes ;
		if(get_serializable_attributes(mesh->facet_attribute_manager(), attributes, output, "facet")) {
			int fid = 1 ;
			FOR_EACH_FACET_CONST(Map, mesh, it) {
				const Map::Facet* f = it ;
				output << "# attrs f " << fid << " " ;
				serialize_write_attributes(output, f, attributes) ;
				output << std::endl ;
				fid++ ;
			}
		}
	}

	output << "# END" << std::endl ;

	return true ;
}

