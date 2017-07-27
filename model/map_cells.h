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



#ifndef _GEOM_MAP_TYPES_H_
#define _GEOM_MAP_TYPES_H_

#include "model_common.h"
#include "../basic/counted.h"
#include "../basic/smart_pointer.h"
#include "../math/math_types.h"
#include "../basic/record_id.h"
#include "../basic/attribute_manager.h"


class Map ;
class MapMutator ;

namespace MapTypes {

	class Vertex ;
	class Halfedge ;
	class Facet ;
	class MapMutator ;

	//___________________________________________________

	/**
	* Combinatorial Element. Combel is the base class for 
	* vertices, half-edges and polygons.
	*/
	class MODEL_API Combel : public Record {
	public:
		Combel()  { }
		~Combel() { }
	} ;

	//___________________________________________________

	/**
	* A vertex of a Map. Each Vertex has a geometry (i.e. a vec3)
	*/
	class MODEL_API Vertex : public Combel {
	public:
		Vertex() : halfedge_(nil) {  }
		Vertex(const vec3& p) : halfedge_(nil), point_(p) {  }
		~Vertex() { halfedge_ = nil ; }

		const vec3& point() const     { return point_ ; }
		vec3& point()                 { return point_ ; }
		void set_point(const vec3& p) { point_ = p ;    }

		// the halfedge points to this vertex
		Halfedge* halfedge() const { return halfedge_ ; }

		bool is_valid() const ;
		void assert_is_valid() const ;

		unsigned int degree() const ;
		bool is_on_border() const ;

		/* Returns true if "this" and v are in each other's one-rings */
		bool is_connected(const Vertex* v) const;

	protected:
		void set_halfedge(Halfedge* h) { halfedge_ = h ; }
		friend class ::Map ;
		friend class ::MapMutator ;

	private:
		Halfedge* halfedge_ ;
		vec3 point_ ;
	} ;

	//______________________________________________

	/**
	* Each edge of a Map is composed of two Halfedges.
	*/

	class MODEL_API Halfedge : public Combel {
	public:
		Halfedge() : 
		  opposite_(nil), next_(nil), 
			  prev_(nil), facet_(nil), vertex_(nil) {
		  }
		  ~Halfedge() { 
			  opposite_ = nil ; next_ = nil ;
			  prev_ = nil ; facet_ = nil ; vertex_ = nil ;
		  }

		  Halfedge* opposite() const { return opposite_ ; }
		  Halfedge* next() const { return next_ ; }
		  Halfedge* prev() const { return prev_ ; }

		  Halfedge* next_around_vertex() const {
			  return opposite()->prev() ;
		  }
		  Halfedge* prev_around_vertex() const {
			  return next()->opposite() ;
		  }

		  Facet*  facet() const { return facet_ ; }
		  Vertex* vertex() const { return vertex_ ; }

		  bool is_border() const { return facet_ == nil ; }
		  bool is_border_edge() const { 
			  return is_border() || opposite()->is_border() ; 
		  }

		  /** One halfedge per facet exactly is the facet key. */
		  bool is_facet_key() const ;

		  /** One halfedge per vertex exactly is the vertex key. */
		  bool is_vertex_key() const ;

		  /** 
		  * One halfedge per edge exactly is the edge key. 
		  * Note: this can be used for loops, to traverse one halfedge 
		  * per edge exactly (for instance, to draw the mesh).
		  */
		  bool is_edge_key() const ;
		  Halfedge* edge_key() const { 
			  return is_edge_key() ? const_cast<Halfedge*>(this) : opposite() ; 
		  }

		  bool is_valid() const ;
		  void assert_is_valid() const ;

	protected:
		void set_opposite(Halfedge* h) { opposite_ = h; }
		void set_next(Halfedge* h) { next_ = h; }
		void set_prev(Halfedge* h) { prev_ = h; }
		void set_facet(Facet* f) { facet_ = f ; }
		void set_vertex(Vertex* v) { vertex_ = v ; }

		friend class ::Map ;
		friend class ::MapMutator ;

	private:
		Halfedge* opposite_ ;
		Halfedge* next_ ;
		Halfedge* prev_ ;
		Facet* facet_ ;
		Vertex* vertex_ ;
	} ;

	//______________________________________________


	/**
	* A Facet of a Map.
	*/

	class MODEL_API Facet : public Combel {
	public:
		Facet() : halfedge_(nil) { }
		~Facet() { halfedge_ = nil ; }

		Halfedge* halfedge() const { return halfedge_ ; }

		int degree() const ;
		int nb_edges() const { return degree() ; }
		int nb_vertices() const { return degree() ; }

		bool is_on_border() const ;
		bool is_triangle() const ;
		bool is_valid() const ;

		void assert_is_valid() const ;

	protected:
		void set_halfedge(Halfedge* h) { halfedge_ = h ; }
		friend class ::Map ;
		friend class ::MapMutator ;

	private:
		Halfedge* halfedge_ ;
	} ;

	//_________________________________________________________

	inline bool Halfedge::is_facet_key() const {
		return (facet_->halfedge() == this) ;
	}

	inline bool Halfedge::is_vertex_key() const {
		return (vertex_->halfedge() == this) ;
	}

	inline bool Halfedge::is_edge_key() const {
		// TODO: if the GarbageCollector is added, 
		// watch out, the edge keys can change...
		return (this < opposite_) ;
	}

	//_________________________________________________________

}

#endif

