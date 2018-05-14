/*
Copyright (C) 2017  Liangliang Nan
http://web.siat.ac.cn/~liangliang/ - liangliang.nan@gmail.com

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

#ifndef _FACE_SELECTION_H_
#define _FACE_SELECTION_H_

#include "method_common.h"
#include "polyfit_info.h"
#include "../math/math_types.h"
#include "../math/linear_program.h"
#include "../math/linear_program_solver.h"
#include "../model/map_attributes.h"


#include <vector>
#include <string>
#include <map>


class Map;
class PointSet;
class VertexGroup;
class PolyFitInfo;

namespace MapTypes {
	class Vertex;
	class Facet;
	class Halfedge;
}


// to determine if a face should be selected or not
class METHOD_API FaceSelection
{
public:
	FaceSelection(PointSet* pset, Map* model);
	~FaceSelection() {}

	virtual void optimize(PolyFitInfo* polyfit_info, LinearProgramSolver::SolverName solver_name);

private:
	PointSet* pset_;
	Map*      model_;

	LinearProgram	program_;

	MapFacetAttribute<VertexGroup*> facet_attrib_supporting_vertex_group_;
	MapFacetAttribute<double>		facet_attrib_supporting_point_num_;
	MapFacetAttribute<double>		facet_attrib_facet_area_;
	MapFacetAttribute<double>		facet_attrib_covered_area_;

	MapFacetAttribute<Plane3d*>					facet_attrib_supporting_plane_;
	MapVertexAttribute< std::set<Plane3d*> >	vertex_source_planes_;
	MapHalfedgeAttribute< std::set<Plane3d*> >	edge_source_planes_;

protected:
	// faces (each face is represented by its halfedge) intersecting at 
	// a common edge (i.e., line segments having the same end points). 
	struct FaceStar : public std::vector < MapTypes::Halfedge* > {};

	class IntersectionAdjacency {
	public:
		std::vector<FaceStar> extract(Map* model, const std::vector<Plane3d*>& supporting_planes);
	private:
		void init(Map* model, const std::vector<Plane3d*>& supporting_planes);
	private:
		// each end point of an edge is denoted by its index.
		// 	How to use: edge_container_[min_end_point_index][max_end_point_index]
		std::map< std::size_t, std::map< std::size_t, std::set<MapTypes::Halfedge*> > > edge_container_;
		MapVertexAttribute< std::set<Plane3d*> >	vertex_source_planes_;
		std::map< Plane3d*, int > plane_index_;
	};

	IntersectionAdjacency adjacency_;
};

#endif