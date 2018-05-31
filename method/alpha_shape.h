/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/ - liangliang.nan@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef _ALPHA_SHAPE_H_
#define _ALPHA_SHAPE_H_

#include "cgal_types.h"

/*
Breaking change since CGAL 4.11: The dangerous implicit conversions between weighted Points
and points in the concept Kernel have been disabled.
Constructors offering to build a weighted point from a point(and reversely)
are still requested by the concept Kernel but must now be marked with the
explicit specifier.
*/
#if CGAL_VERSION_NR >= 1041100000 
#include "alpha_shape_CGAL4.11_and_later.h"
#else
#include "alpha_shape_CGAL4.10_and_earlier.h"
#endif

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>

#include "../model/point_set.h"
#include "../model/vertex_group.h"

#include <string>



/* A vertex class with an additional member */
template < class Gt, class myVb = CGAL::Triangulation_hierarchy_vertex_base_2<Gt> >
class AS_vertex_base : public  myVb
{
	typedef myVb                              Base;
public:
	typedef typename myVb::Vertex_handle      Vertex_handle;
	typedef typename myVb::Face_handle        Face_handle;
	typedef typename myVb::Point              Point;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename myVb::template Rebind_TDS<TDS2>::Other	Vb2;
		typedef AS_vertex_base<Gt, Vb2>                         Other;
	};

public:
	AS_vertex_base() : Base(), index_(-1) {}
	AS_vertex_base(const Point & p) : Base(p), index_(-1) {}
	AS_vertex_base(const Point & p, Face_handle f) : Base(f, p), index_(-1) {}
	AS_vertex_base(Face_handle f) : Base(f), index_(-1) {}

	void set_index(int idx) { index_ = idx; }
	int  index() const { return index_; }

private:
	int index_;
};



typedef CGAL::Alpha_shape_vertex_base_2<K>					Avb;
typedef AS_vertex_base<Avb>									Av;
typedef CGAL::Triangulation_face_base_2<K>					Tf;
typedef CGAL::Alpha_shape_face_base_2<K, Tf>				Af;

typedef CGAL::Triangulation_default_data_structure_2<K, Av, Af> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>				Dt;
typedef CGAL::Triangulation_hierarchy_2<Dt>					Ht;
typedef CGAL::Alpha_shape_CGAL<Ht>							CGAL_AlphaShape;

typedef CGAL_AlphaShape::Finite_vertices_iterator			Vertices_iterator;
typedef CGAL_AlphaShape::Finite_faces_iterator				Faces_iterator;
typedef CGAL_AlphaShape::Face_handle						Face_handle;
typedef CGAL_AlphaShape::Vertex_handle						Vertex_handle;


class AlphaShape : public CGAL_AlphaShape  
{
public:
	// the input is a list of points.
	// NOTE: an 'id' (same as the order in the list) will be given 
	//       during the construction of the alpha shape (this allow 
	//       me to retrieve the original points)
	template <typename InputIterator>
	AlphaShape(InputIterator first, InputIterator beyond);  

	std::vector<Face_handle> get_all_finite_facets();

	// Determine the number of connected solid components. 
	// ids returns the component id for each finite facet.
	int enumerate_solid_components(std::map<Face_handle, int>& ids);

	// get the solid_component that has the largest area.
	std::vector<Face_handle> extract_largest_solid_component();

protected:
	// Liangliang: ref to CGAL4.2: Alpha_shape_2.h : line 1341
	void traverse(const Face_handle& pFace, std::map<Face_handle, int>& ids, int cur_id);

private:
	int num_duplicated_;
};


template <typename InputIterator>
AlphaShape::AlphaShape(InputIterator first, InputIterator beyond) {
	num_duplicated_ = 0;
	InputIterator it = first;
	for (int id = 0; it != beyond; ++it, ++id) {
		const Point2& p = *it;
		Vertex_handle vh = Ht::insert(p);
		if (vh->index() == -1) {
			vh->set_index(id);
		}
		else
			++num_duplicated_;
	}

	//if (num_duplicated_ > 0)
	//	Logger::out("-") << num_duplicated_ << " points ignored in triangulation" << std::endl;

	if (dimension() == 2) {
		// Compute the associated _interval_face_map
		initialize_interval_face_map();

		// Compute the associated _interval_edge_map
		initialize_interval_edge_map();

		// Compute the associated _interval_vertex_map
		initialize_interval_vertex_map();

		// merge the two maps
		initialize_alpha_spectrum();
	}

	//Alpha_iterator alpha = find_optimal_alpha(1);
	//if (alpha == alpha_end())
	//	set_alpha(0.0f);
	//else 
	//	set_alpha(*alpha);
}


inline std::vector<Face_handle> AlphaShape::get_all_finite_facets()
{
	std::vector<Face_handle> facets;

	for (Finite_faces_iterator fit = finite_faces_begin(); fit != finite_faces_end(); ++fit) {
		Face_handle pFace = fit;
		CGAL_triangulation_postcondition(pFace != NULL);
		if (classify(pFace) == AlphaShape::INTERIOR)
			facets.push_back(pFace);
	}

	return facets;
}


// Liangliang: ref to CGAL4.2: Alpha_shape_2.h : line 1341
inline void AlphaShape::traverse(const Face_handle& pFace, std::map<Face_handle, int>& ids, int cur_id)
{
	std::list<Face_handle> faces;
	faces.push_back(pFace);
	Face_handle pNeighbor, fh;

	while (!faces.empty()) {
		fh = faces.front();
		faces.pop_front();
		for (int i = 0; i < 3; i++)
		{
			pNeighbor = fh->neighbor(i);
			CGAL_triangulation_assertion(pNeighbor != NULL);
			if (classify(pNeighbor) == AlphaShape::INTERIOR){
				int& id = ids[pNeighbor];
				if (id == -1){
					id = cur_id;
					faces.push_back(pNeighbor);
				}
			}
		}
	}
}

// Determine the number of connected solid components 
inline int AlphaShape::enumerate_solid_components(std::map<Face_handle, int>& components)
{
	if (number_of_vertices() == 0)
		return 0;

	std::map<Face_handle, int> status;
	Finite_faces_iterator face_it = finite_faces_begin();
	for (; face_it != finite_faces_end(); ++face_it)
		status[face_it] = -1; 

	unsigned int cur_id = 0;
	// only finite faces
	for (face_it = finite_faces_begin(); face_it != finite_faces_end(); ++face_it)
	{
		Face_handle pFace = face_it;
		CGAL_triangulation_postcondition(pFace != NULL);

		if (classify(pFace) == AlphaShape::INTERIOR) {
			int& id = status[pFace];
			if (id == -1) {
				// we traverse only interior faces
				traverse(pFace, status, cur_id);
				++cur_id;
			}
		}
	}

	// discard the ones with id==-1
	components.clear();
	std::map<Face_handle, int>::const_iterator it = status.begin(); 
	for (; it != status.end(); ++it) {
		int id = it->second;
		if (id > -1) {
			components[it->first] = id;
		}
	}
	return cur_id;
}



inline std::vector<Face_handle> AlphaShape::extract_largest_solid_component() {
	std::map<Face_handle, int> ids;
	int num = enumerate_solid_components(ids);
	std::vector<double> area(num, 0.0);

	std::map<Face_handle, int>::const_iterator it = ids.begin();
	for (; it != ids.end(); ++it) {
		Face_handle f = it->first;
		
		// convert the facet into a 2D polygon
		Polygon2 plg;
		for (unsigned int j = 0; j < 3; ++j) {
			Vertex_handle vh = f->vertex(j);
			const Point2& p = vh->point();
			plg.push_back(p);
		}

		int id = it->second;
		area[id] += std::fabs(plg.area());
	}

	//////////////////////////////////////////////////////////////////////////

	int largest_id = 0;
	double largest_area = area[0];

	for (int i = 1; i<num; ++i) {
		double ar = area[i];
		if (ar > largest_area) {
			largest_area = ar;
			largest_id = i;
		}
	}

	//////////////////////////////////////////////////////////////////////////

	std::vector<Face_handle> result;
	for (it = ids.begin(); it != ids.end(); ++it) {
		Face_handle f = it->first;
		int id = it->second;
		if (id == largest_id)
			result.push_back(f);
	}
	return result;
}


#endif 
