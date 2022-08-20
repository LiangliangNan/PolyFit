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

#ifndef _FACE_SELECTION_H_
#define _FACE_SELECTION_H_

#include "method_common.h"
#include "hypothesis_generator.h"
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

	virtual void optimize(const HypothesisGenerator::Adjacency& adjacency, LinearProgramSolver::SolverName solver_name);

    // NOTE: the adjacency is the one extracted after the face optimization step
    virtual void re_orient(const HypothesisGenerator::Adjacency& adjacency, LinearProgramSolver::SolverName solver_name);

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
};

#endif
