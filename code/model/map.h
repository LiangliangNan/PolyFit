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




#ifndef _GEOM_MAP_H_
#define _GEOM_MAP_H_

#include "model_common.h"
#include "iterators.h"
#include "map_cells.h"
#include "../basic/dlist.h"
#include "../basic/attribute.h"
#include "../basic/counted.h"
#include "../basic/smart_pointer.h"

#include <vector>
#include <list>



class Map;

/**
* Has its members add() and remove() called
* whenever a vertex is added (respectively
* removed) from the structure. This class
* is meant to be subclassed by client code.
*/
template <class CELL> 
class MapCombelObserver {
public:
	MapCombelObserver(Map* m) : map_(m) { }
	virtual ~MapCombelObserver() { }
	virtual void add(CELL* c) = 0;
	virtual void remove(CELL* c) = 0; 

protected:
	Map* map() { return map_; }

private:
	Map* map_ ;
} ;

/** 
* DEFINITION
* 
* The boundary representation of a 3d-polyhedron P consists of vertices, edges and 
* facets. The vertices are points in space. The edges are straight line segments. 
* The facets are planar polygons. We restrict here the facets to be simple planar 
* polygons without holes and the boundary of the polyhedron to be an oriented 
* 2-manifold. Thus facets are consistently oriented and an edge is incident to 
* exactly two facets. We restrict the representation further that an edge has two 
* distinct incident endpoints and following duality that an edge has two distinct 
* incident facets. 
*
* Map implements a topological Map. It is highly inspired by CGAL's Polyhedron class, 
* with the difference that each halfedge is allocated separately. All the elements are
* allocated in a high-performance chunk-based allocator. Map is nearly source-level 
* compatible with CGAL polyhedron. The operations modifying the Map are implemented in
* the MapEditor class.
*/

class MapMutator ;


class MODEL_API Map : public Counted
{
public:
	typedef SmartPointer<Map>					Ptr;

	// __________________ types ___________________

	typedef MapTypes::Vertex					Vertex ;
	typedef MapTypes::Halfedge					Halfedge ;
	typedef MapTypes::Facet						Facet ;

	typedef DList<Vertex>::iterator				Vertex_iterator ;
	typedef DList<Halfedge>::iterator			Halfedge_iterator ;
	typedef DList<Facet>::iterator				Facet_iterator ;

	typedef DList<Vertex>::const_iterator		Vertex_const_iterator ;
	typedef DList<Halfedge>::const_iterator		Halfedge_const_iterator ;
	typedef DList<Facet>::const_iterator		Facet_const_iterator ;

	typedef GenericAttributeManager<Vertex>		VertexAttributeManager ;
	typedef GenericAttributeManager<Halfedge>	HalfedgeAttributeManager ;
	typedef GenericAttributeManager<Facet>		FacetAttributeManager ;

public:

	// ______________ constructor and destructor ___________________

	Map() : bbox_is_valid_(false) {}

	virtual ~Map();

	// ________________________ access _____________________________

	Vertex_iterator vertices_begin()    { return vertices_.begin() ;  }
	Vertex_iterator vertices_end()      { return vertices_.end() ;    }
	Halfedge_iterator halfedges_begin() { return halfedges_.begin() ; }
	Halfedge_iterator halfedges_end()   { return halfedges_.end() ;   }
	Facet_iterator facets_begin()       { return facets_.begin() ;    }
	Facet_iterator facets_end()         { return facets_.end() ;      } 

	Vertex_const_iterator vertices_begin() const	{ return vertices_.begin() ;	}
	Vertex_const_iterator vertices_end() const		{ return vertices_.end() ;		}
	Halfedge_const_iterator halfedges_begin() const { return halfedges_.begin() ;	}
	Halfedge_const_iterator halfedges_end() const	{ return halfedges_.end() ;		}
	Facet_const_iterator facets_begin() const		{ return facets_.begin() ;		}
	Facet_const_iterator facets_end() const			{ return facets_.end() ;		} 

	int size_of_vertices() const  { return vertices_.size() ;  }
	int size_of_halfedges() const { return halfedges_.size() ; }
	int size_of_facets() const    { return facets_.size() ;    }

	const Box3d& bbox() const;
	void invalidate_bbox() { bbox_is_valid_ = false; }

	// ___________________ attributes _______________________

	VertexAttributeManager* vertex_attribute_manager() const {
		return const_cast<VertexAttributeManager*>(
			&vertex_attribute_manager_ 
			) ;
	}

	HalfedgeAttributeManager* halfedge_attribute_manager() const {
		return const_cast<HalfedgeAttributeManager*>(
			&halfedge_attribute_manager_ 
			) ;
	}

	FacetAttributeManager* facet_attribute_manager() const {
		return const_cast<FacetAttributeManager*>(
			&facet_attribute_manager_ 
			) ;
	}         

	// ___________________ observers ________________________

	void add_vertex_observer(MapCombelObserver<Vertex>* obs) ;
	void add_halfedge_observer(MapCombelObserver<Halfedge>* obs) ;
	void add_facet_observer(MapCombelObserver<Facet>* obs) ;

	void remove_vertex_observer(MapCombelObserver<Vertex>* obs) ;
	void remove_halfedge_observer(MapCombelObserver<Halfedge>* obs) ;
	void remove_facet_observer(MapCombelObserver<Facet>* obs) ;

	// __________________ modification ______________________

	void clear() ;
	void clear_inactive_items() ;

	// __________________ stored normals ____________________

	void compute_vertex_normals();
	void compute_facet_normals();

	// __________________ predicated ________________________

	bool is_triangulated() const ;

	// __________________ low level (for experts only) ______

	/**
	* checks the validity of the combinatorial structure.
	*/
	bool is_valid() const ;

	/**
	* checks the validity of the combinatorial structure,
	* and stops in an assertion failure if an error is detected.
	* This function can be used to debug low-level operations,
	* such as those in the MeshBuilder and the MapEditor classes.
	*/
	void assert_is_valid() const ;

protected:

	/**
	* allocates a pair of halfedges connected with next(), prev(),
	* and opposite() links.
	*/
	Halfedge* new_edge() ;

	/**
	* deallocates a pair of halfedges, connected with an opposite() link.
	*/
	void delete_edge(Halfedge* h) ;

	Vertex* new_vertex() ;
	Halfedge* new_halfedge() ;
	Facet* new_facet() ;

	/** copies geometry and attributes from rhs. */
	Vertex* new_vertex(const Vertex* rhs) ;
	/** copies attributes from rhs. */
	Halfedge* new_halfedge(const Halfedge* rhs) ;
	/** copies attributes from rhs. */
	Facet* new_facet(const Facet* rhs) ;

	void delete_vertex(Vertex* v) ;
	void delete_halfedge(Halfedge* h) ;
	void delete_facet(Facet* f) ;

	void activate_vertex(Vertex* v) ;
	void activate_halfedge(Halfedge* h) ;
	void activate_facet(Facet* f) ;

	void deactivate_vertex(Vertex* v) ;
	void deactivate_halfedge(Halfedge* h) ;
	void deactivate_facet(Facet* f) ;

	friend class ::MapMutator ;

protected:
	// MeshCombelObservers notification
	void notify_add_vertex(Vertex* v) ;
	void notify_add_halfedge(Halfedge* h) ;
	void notify_add_facet(Facet* f) ;

	void notify_remove_vertex(Vertex* v) ;
	void notify_remove_halfedge(Halfedge* h) ;
	void notify_remove_facet(Facet* f) ;

private:
	// elements
	DList<Vertex>	vertices_ ;
	DList<Halfedge> halfedges_ ;
	DList<Facet>	facets_ ;

	VertexAttributeManager		vertex_attribute_manager_ ;
	HalfedgeAttributeManager	halfedge_attribute_manager_ ;
	FacetAttributeManager		facet_attribute_manager_ ;

	std::vector< MapCombelObserver<Vertex>* >	vertex_observers_ ;
	std::vector< MapCombelObserver<Halfedge>* >	halfedge_observers_ ;
	std::vector< MapCombelObserver<Facet>* >	facet_observers_ ;

	mutable bool	bbox_is_valid_;
	mutable Box3d	bbox_;
} ;


//______________________________________________________________


/*
* MapMutator is the base class for the classes that can modify the topology of a mesh.
*/
class MODEL_API MapMutator
{
public:

	// _________________ Types ____________________

	typedef Map::Vertex                   Vertex ;
	typedef Map::Halfedge                 Halfedge ;
	typedef Map::Facet                    Facet ;

	typedef Map::Vertex_iterator          Vertex_iterator ;
	typedef Map::Halfedge_iterator        Halfedge_iterator ;
	typedef Map::Facet_iterator           Facet_iterator ;

	typedef Map::Vertex_const_iterator    Vertex_const_iterator ;
	typedef Map::Halfedge_const_iterator  Halfedge_const_iterator ;
	typedef Map::Facet_const_iterator     Facet_const_iterator ;

protected:
	MapMutator(Map* target = nil) : target_(target) { }
	virtual ~MapMutator() { target_ = nil;  } 

public:
	Map* target() const { return target_ ; }
	virtual void set_target(Map* target) { target_ = target ; }

protected:
	void make_sequence(Halfedge* h1, Halfedge* h2) {
		h1->set_next(h2) ;
		h2->set_prev(h1) ;
	}

	void make_opposite(Halfedge* h1, Halfedge* h2) {
		h1->set_opposite(h2) ;
		h2->set_opposite(h1) ;
	}

	void set_vertex_on_orbit(Halfedge* h, Vertex* v) ;
	void set_facet_on_orbit(Halfedge* h, Facet* f) ;

	void make_vertex_key(Halfedge* h) { 
		h->vertex()->set_halfedge(h) ; 
	}

	void make_vertex_key(Halfedge* h, Vertex* v) { 
		v->set_halfedge(h) ;
		h->set_vertex(v) ;
	}

	void make_facet_key(Halfedge* h) {
		h->facet()->set_halfedge(h) ;
	}

	void make_facet_key(Halfedge* h, Facet* f) {
		f->set_halfedge(h) ;
		h->set_facet(f) ;
	}

	// _________ create / destroy (see Map) _________

	Halfedge* new_edge() { return target_->new_edge() ; }
	void delete_edge(Halfedge* h) { target_->delete_edge(h) ; }

	Vertex*    new_vertex()     { return target_->new_vertex() ;     }
	Halfedge*  new_halfedge()   { return target_->new_halfedge() ;   }
	Facet*     new_facet()      { return target_->new_facet() ;      }

	Vertex* new_vertex(const Vertex* rhs)			{ return target_->new_vertex(rhs) ; }
	Halfedge* new_halfedge(const Halfedge* rhs)		{ return target_->new_halfedge(rhs) ; }
	Halfedge* new_edge(Halfedge* rhs) ;
	Facet* new_facet(const Facet* rhs)				{ return target_->new_facet(rhs) ; }


	void delete_vertex(Vertex* v)			{ target_->delete_vertex(v) ;		}
	void delete_halfedge(Halfedge* h)		{ target_->delete_halfedge(h) ;	}
	void delete_facet(Facet* f)				{ target_->delete_facet(f) ;		}

	// _________ activate / deactivate (see Map) _________

	void activate_vertex(Vertex* v)		{ target_->activate_vertex(v);   }
	void activate_halfedge(Halfedge* h) { target_->activate_halfedge(h); }
	void activate_facet(Facet* f)		{ target_->activate_facet(f);	 }

	void deactivate_vertex(Vertex* v)	  { target_->deactivate_vertex(v);  }
	void deactivate_halfedge(Halfedge* h) { target_->deactivate_halfedge(h); }
	void deactivate_facet(Facet* f)		  { target_->deactivate_facet(f);    }

	// _________ basic functions ____________________

	void set_vertex_halfedge(Vertex* v, Halfedge* h) { v->set_halfedge(h) ; }
	void set_halfedge_opposite(Halfedge* h1, Halfedge* h2) { h1->set_opposite(h2) ;}
	void set_halfedge_next(Halfedge* h1, Halfedge* h2) { h1->set_next(h2) ;}
	void set_halfedge_prev(Halfedge* h1, Halfedge* h2) { h1->set_prev(h2) ;}
	void set_halfedge_facet(Halfedge* h, Facet* f) { h->set_facet(f) ;}
	void set_halfedge_vertex(Halfedge* h, Vertex* v) { h->set_vertex(v) ;}
	void set_facet_halfedge(Facet* f, Halfedge* h) { f->set_halfedge(h) ;}

private:
	Map* target_ ;
} ;



//_________________________________________________________

/**
* The constructor and destructor of this
* specialization call mesh->add_vertex_observer(this)
* and mesh->remove_vertex_observer(this) respectively.
*/
template<> 
class MODEL_API MapCombelObserver<Map::Vertex> {
public:
	MapCombelObserver(Map* m) ;
	virtual ~MapCombelObserver() ;
	virtual void add(Map::Vertex* c) { }
	virtual void remove(Map::Vertex* c) { }

protected:
	Map* map() { return map_; }

private:
	Map* map_ ;
} ;

/**
* The constructor and destructor of this
* specialization call mesh->add_halfedge_observer(this)
* and mesh->remove_halfedge_observer(this) respectively.
*/
template<> 
class MODEL_API MapCombelObserver<Map::Halfedge> {
public:
	MapCombelObserver(Map* m) ;
	virtual ~MapCombelObserver() ;
	virtual void add(Map::Halfedge* c) { }
	virtual void remove(Map::Halfedge* c) { }

protected:
	Map* map() { return map_; }

private:
	Map* map_ ;
} ;


/**
* The constructor and destructor of this
* specialization call map->add_facet_observer(this)
* and map->remove_facet_observer(this) respectively.
*/
template<> 
class MODEL_API MapCombelObserver<Map::Facet> {
public:
	MapCombelObserver(Map* m) ;
	virtual ~MapCombelObserver() ;
	virtual void add(Map::Facet* c) { }
	virtual void remove(Map::Facet* c) { }

protected:
	Map* map() { return map_; }

private:
	Map* map_ ;
} ;


#endif

