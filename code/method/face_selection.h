/* ---------------------------------------------------------------------------
 * Copyright (C) 2017 Liangliang Nan <liangliang.nan@gmail.com>
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of PolyFit. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 *
 *     Liangliang Nan and Peter Wonka.
 *     PolyFit: Polygonal Surface Reconstruction from Point Clouds.
 *     ICCV 2017.
 *
 *  For more information:
 *  https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html
 * ---------------------------------------------------------------------------
 */

#ifndef _FACE_SELECTION_H_
#define _FACE_SELECTION_H_

#include <method/method_common.h>
#include <math/math_types.h>
#include <math/linear_program.h>
#include <math/linear_program_solver.h>
#include <model/map_attributes.h>


#include <vector>
#include <string>
#include <map>


class Map;
class PointSet;
class VertexGroup;
class HypothesisGenerator;

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

    void optimize(HypothesisGenerator *generator,
                  LinearProgramSolver::SolverName solver_name,  // solver name
                  double data_fitting,       // weight for data fitting term
                  double model_coverage,     // weight for model coverage term
                  double model_complexity    // weight for model complexity term)
    );

protected:
    // NOTE: the adjacency is the one extracted after the face optimization step
    void re_orient(HypothesisGenerator* generator, LinearProgramSolver::SolverName solver_name);

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
