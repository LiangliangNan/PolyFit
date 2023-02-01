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


#include "face_selection.h"
#include "method_global.h"
#include "../basic/stop_watch.h"
#include "../model/point_set.h"
#include "../model/map_geometry.h"
#include "../basic/logger.h"
#include "../model/map_editor.h"

#include <algorithm>


#define ALLOW_OPEN_SURFACE


FaceSelection::FaceSelection(PointSet* pset, Map* model)
	: pset_(pset)
	, model_(model)
{
}


void FaceSelection::optimize(const HypothesisGenerator::Adjacency& adjacency, LinearProgramSolver::SolverName solver_name) {
	if (pset_ == 0 || model_ == 0)
		return;

	facet_attrib_supporting_vertex_group_.bind_if_defined(model_, Method::facet_attrib_supporting_vertex_group);
	if (!facet_attrib_supporting_vertex_group_.is_bound()) {
		Logger::err("-") << "attribute " << Method::facet_attrib_supporting_vertex_group << " doesn't exist" << std::endl;
		return;
	}
	facet_attrib_supporting_point_num_.bind_if_defined(model_, Method::facet_attrib_supporting_point_num);
	if (!facet_attrib_supporting_point_num_.is_bound()) {
		Logger::err("-") << "attribute " << Method::facet_attrib_supporting_point_num << " doesn't exist" << std::endl;
		return;
	}
	facet_attrib_facet_area_.bind_if_defined(model_, Method::facet_attrib_facet_area);
	if (!facet_attrib_facet_area_.is_bound()) {
		Logger::err("-") << "attribute " << Method::facet_attrib_facet_area << " doesn't exist" << std::endl;
		return;
	}
	facet_attrib_covered_area_.bind_if_defined(model_, Method::facet_attrib_covered_area);
	if (!facet_attrib_covered_area_.is_bound()) {
		Logger::err("-") << "attribute " << Method::facet_attrib_covered_area << " doesn't exist" << std::endl;
		return;
	}

	edge_source_planes_.bind_if_defined(model_, "EdgeSourcePlanes");
	vertex_source_planes_.bind_if_defined(model_, "VertexSourcePlanes");
	facet_attrib_supporting_plane_.bind_if_defined(model_, "FacetSupportingPlane");

	//////////////////////////////////////////////////////////////////////////

	double total_points = double(pset_->points().size());
	std::size_t idx = 0;
	MapFacetAttribute<std::size_t>	facet_indices(model_);
	FOR_EACH_FACET(Map, model_, it) {
		Map::Facet* f = it;
		facet_indices[f] = idx;
		++idx;
	}

	//-------------------------------------

	StopWatch w;

	//-------------------------------------

	// binary variables:
	// x[0] ... x[num_faces - 1] : binary labels of all the input faces
	// x[num_faces] ... x[num_faces + num_edges] : binary labels of all the intersecting edges (remain or not)
    // x[num_faces + num_edges] ... x[num_faces + num_edges + num_edges] : binary labels of corner edges (sharp edge or not)
    // extra variables for open surfaces
    // x[num_faces + num_edges + num_edges] ... x[num_faces + num_edges + num_edges + num_edges] : binary labels of boundary edges

	Logger::out("-") << "formulating binary program...." << std::endl;
	w.start();

	std::size_t num_faces = model_->size_of_facets();
	std::size_t num_edges = 0;

	typedef typename HypothesisGenerator::SuperEdge SuperEdge;
	std::map<const SuperEdge*, std::size_t> edge_usage_status;	// keep or remove an intersecting edges
	for (std::size_t i = 0; i < adjacency.size(); ++i) {
		const SuperEdge& fan = adjacency[i];
		if (fan.size() == 4) {
			std::size_t var_idx = num_faces + num_edges;
			edge_usage_status[&fan] = var_idx;
			++num_edges;
		}
	}

	//double coeff_data_fitting = Method::lambda_data_fitting / total_points;
	//double coeff_coverage = Method::lambda_model_coverage / model_->bbox().area();
    //double coeff_complexity = Method::lambda_model_complexity / double(fans.size());
    // choose a better scale
	double coeff_data_fitting = Method::lambda_data_fitting;
	double coeff_coverage = total_points * Method::lambda_model_coverage / model_->bbox().area();
	double coeff_complexity = total_points * Method::lambda_model_complexity / double(adjacency.size());

	program_.clear();
	LinearObjective* objective = program_.create_objective(LinearObjective::MINIMIZE);

	std::map<const SuperEdge*, std::size_t> edge_sharp_status;	// the edge is sharp or not
	std::size_t num_sharp_edges = 0;
	for (std::size_t i = 0; i < adjacency.size(); ++i) {
		const SuperEdge& fan = adjacency[i];
		if (fan.size() == 4) {
			std::size_t var_idx = num_faces + num_edges + num_sharp_edges;
			edge_sharp_status[&fan] = var_idx;

			// accumulate model complexity term
			objective->add_coefficient(var_idx, coeff_complexity);
			++num_sharp_edges;
		}
	}
	assert(num_edges == num_sharp_edges);

#ifdef ALLOW_OPEN_SURFACE
    double coeff_boundary = total_points * Method::lambda_boundary / double(num_edges);

    std::map<const SuperEdge*, std::size_t> edge_boundary_status;	// the edge is sharp or not
    std::size_t num_boundary_edges = 0;
    for (std::size_t i = 0; i < adjacency.size(); ++i) {
        const SuperEdge& fan = adjacency[i];
        if (fan.size() == 4) {
            std::size_t var_idx = num_faces + num_edges + num_sharp_edges + num_boundary_edges;
            edge_boundary_status[&fan] = var_idx;

            // accumulate boundary term
            objective->add_coefficient(var_idx, coeff_boundary);
            ++num_boundary_edges;
        }
    }
    assert(num_edges == num_boundary_edges);
#endif

	FOR_EACH_FACET(Map, model_, it) {
		Map::Facet* f = it;
		std::size_t var_idx = facet_indices[f];

		// accumulate data fitting term
		double num = facet_attrib_supporting_point_num_[f];
		objective->add_coefficient(var_idx, -coeff_data_fitting * num);

		// accumulate model coverage term
		double uncovered_area = (facet_attrib_facet_area_[f] - facet_attrib_covered_area_[f]);
		objective->add_coefficient(var_idx, coeff_coverage * uncovered_area);
	}

#ifdef ALLOW_OPEN_SURFACE
    std::size_t total_variables = num_faces + num_edges + num_sharp_edges + num_boundary_edges;
#else
    std::size_t total_variables = num_faces + num_edges + num_sharp_edges;
#endif

	Logger::out("-") << "#total variables: " << total_variables << std::endl;
	Logger::out(" ") << "    - face is selected: " << num_faces << std::endl;
	Logger::out(" ") << "    - edge is used: " << num_edges << std::endl;
    Logger::out(" ") << "    - edge is sharp: " << num_sharp_edges << std::endl;
#ifdef ALLOW_OPEN_SURFACE
    Logger::out(" ") << "    - edge is boundary: " << num_boundary_edges << std::endl;
#endif

#if 1
	const std::vector<Variable*>& variables = program_.create_n_variables(total_variables);
	for (std::size_t i = 0; i < total_variables; ++i) {
		Variable* v = variables[i];
		v->set_variable_type(Variable::BINARY);
	}
#else // Liangliang: I was just curious about how the results look like if all variables 
	//             are relaxed to be continuous.
	const std::vector<Variable*>& variables = program_.create_n_variables(total_variables);
	for (std::size_t i = 0; i < total_variables; ++i) {
		Variable* v = variables[i];
		v->set_variable_type(Variable::CONTINUOUS);
		v->set_bounds(Variable::DOUBLE, 0, 1);
	}
#endif

	//////////////////////////////////////////////////////////////////////////

	// Add constraints: the number of faces associated with an edge must be either 2 or 0
	for (std::size_t i = 0; i < adjacency.size(); ++i) {
		LinearConstraint* c = program_.create_constraint(LinearConstraint::FIXED, 0.0, 0.0);
		const SuperEdge& fan = adjacency[i];
		for (std::size_t j = 0; j < fan.size(); ++j) {
			MapTypes::Facet* f = fan[j]->facet();
			std::size_t var_idx = facet_indices[f];
			c->add_coefficient(var_idx, 1.0);
		}

		if (fan.size() == 4) {
#ifdef ALLOW_OPEN_SURFACE
            std::size_t var_edge_boundary_idx = edge_boundary_status[&fan];
            c->add_coefficient(var_edge_boundary_idx, 1.0);
#endif

            std::size_t var_edge_usage_idx = edge_usage_status[&fan];
            c->add_coefficient(var_edge_usage_idx, -2.0);  //
		}
		else { // boundary edge
			   // will be set to 0 (i.e., we don't allow open surface)
		}
	}

	// Add constraints: for the sharp edges. The explanation of posing this constraint can be found here:
	// https://user-images.githubusercontent.com/15526536/30185644-12085a9c-942b-11e7-831d-290dd2a4d50c.png
    const double M = 1.0;
	for (std::size_t i = 0; i < adjacency.size(); ++i) {
		const SuperEdge& fan = adjacency[i];
		if (fan.size() != 4)
			continue;

		// if an edge is sharp, the edge must be selected first:
		// X[var_edge_usage_idx] >= X[var_edge_sharp_idx]	
        LinearConstraint* c = program_.create_constraint();
		std::size_t var_edge_usage_idx = edge_usage_status[&fan];
        c->add_coefficient(var_edge_usage_idx, 1.0);
		std::size_t var_edge_sharp_idx = edge_sharp_status[&fan];
        c->add_coefficient(var_edge_sharp_idx, -1.0);
        c->set_bound(LinearConstraint::LOWER, 0.0);

		for (std::size_t j = 0; j < fan.size(); ++j) {
			MapTypes::Facet* f1 = fan[j]->facet();
			Plane3d* plane1 = facet_attrib_supporting_plane_[f1];
			std::size_t fid1 = facet_indices[f1];
			for (std::size_t k = j + 1; k < fan.size(); ++k) {
				MapTypes::Facet* f2 = fan[k]->facet();
				Plane3d* plane2 = facet_attrib_supporting_plane_[f2];
				std::size_t fid2 = facet_indices[f2];

				if (plane1 != plane2) {
					// the constraint is:
					//X[var_edge_sharp_idx] + M * (3 - (X[fid1] + X[fid2] + X[var_edge_usage_idx])) >= 1
					// which equals to  
					//X[var_edge_sharp_idx] - M * X[fid1] - M * X[fid2] - M * X[var_edge_usage_idx] >= 1 - 3M
					c = program_.create_constraint();
					c->add_coefficient(var_edge_sharp_idx, 1.0);
					c->add_coefficient(fid1, -M);
					c->add_coefficient(fid2, -M);
					c->add_coefficient(var_edge_usage_idx, -M);
					c->set_bound(LinearConstraint::LOWER, 1.0 - 3.0 * M);
				}
			}
		}

#ifdef ALLOW_OPEN_SURFACE
        // if an edge is boundary, the edge must be selected first:
        // X[var_edge_usage_idx] >= X[var_boundary_idx]
        LinearConstraint* boundary_constraint = program_.create_constraint();
        boundary_constraint->add_coefficient(var_edge_usage_idx, 1.0);
        std::size_t var_edge_boundary_idx = edge_boundary_status[&fan];
        boundary_constraint->add_coefficient(var_edge_boundary_idx, -1.0);
        boundary_constraint->set_bound(LinearConstraint::LOWER, 0.0);
#endif
	}

	Logger::out("-") << "#total constraints: " << program_.constraints().size() << std::endl;
	Logger::out("-") << "formulating binary program done. " << w.elapsed() << " sec" << std::endl;

	//////////////////////////////////////////////////////////////////////////

	// Optimize model
	Logger::out("-") << "solving the binary program. Please wait..." << std::endl;
	w.start();

	LinearProgramSolver solver;
	if (solver.solve(&program_, solver_name)) {
		Logger::out("-") << "solving the binary program done. " << w.elapsed() << " sec" << std::endl;

		// mark results
		const std::vector<double>& X = solver.solution();
		std::vector<Map::Facet*> to_delete;
		FOR_EACH_FACET(Map, model_, it) {
			Map::Facet* f = it;
			std::size_t idx = facet_indices[f];
			//if (static_cast<int>(X[idx]) == 0) { // Liangliang: be careful, floating point!!!
			//if (static_cast<int>(X[idx]) != 1) { // Liangliang: be careful, floating point!!!
			if (static_cast<int>(std::round(X[idx])) == 0) {
				to_delete.push_back(f);
			}
		}

		MapEditor editor(model_);
		for (std::size_t i = 0; i < to_delete.size(); ++i) {
			Map::Facet* f = to_delete[i];
			editor.erase_facet(f->halfedge());
		}

		//////////////////////////////////////////////////////////////////////////

		// mark the sharp edges
		MapHalfedgeAttribute<bool> edge_is_sharp(model_, "SharpEdge");
		FOR_EACH_EDGE(Map, model_, it)
			edge_is_sharp[it] = false;

		for (std::size_t i = 0; i < adjacency.size(); ++i) {
			const SuperEdge& fan = adjacency[i];
			if (fan.size() != 4)
				continue;

			std::size_t idx_sharp_var = edge_sharp_status[&fan];
			if (static_cast<int>(X[idx_sharp_var]) == 1) {
				for (std::size_t j = 0; j < fan.size(); ++j) {
					Map::Halfedge* e = fan[j];
					Map::Facet* f = e->facet();
					if (f) { // some faces may be deleted
						std::size_t fid = facet_indices[f];
						// if (static_cast<int>(X[fid]) == 1) { // Liangliang: be careful, floating point!!!
						if (static_cast<int>(std::round(X[fid])) == 1) {
							edge_is_sharp[e] = true;
							break;
						}
					}
				}
			}
		}
	}
	else {
		Logger::out("-") << "solving the binary program failed." << w.elapsed() << std::endl;
	}

	facet_attrib_supporting_vertex_group_.unbind();
	facet_attrib_supporting_point_num_.unbind();
	facet_attrib_facet_area_.unbind();
	facet_attrib_covered_area_.unbind();

	vertex_source_planes_.unbind();
	edge_source_planes_.unbind();
	facet_attrib_supporting_plane_.unbind();
}


void FaceSelection::re_orient(const HypothesisGenerator::Adjacency& adjacency, LinearProgramSolver::SolverName solver_name) {
	if (model_ == nullptr)
		return;

#if 1
	// check if input is legal
	for (std::size_t i = 0; i < adjacency.size(); ++i) {
		const auto& SuperEdge = adjacency[i];
		if (SuperEdge.size() != 2) {
			Logger::err("-") << "number of faces associated with an edge should be 2." << std::endl;
			return;
		}
	}
#endif

	std::size_t idx = 0;
	MapFacetAttribute<std::size_t>	facet_indices(model_);
	FOR_EACH_FACET(Map, model_, it) {
		Map::Facet* f = it;
		facet_indices[f] = idx;
		++idx;
	}

	//-------------------------------------

	StopWatch w;
	Logger::out("-") << "re-orientation..." << std::endl;

	//-------------------------------------

	// binary variables:
	// x[0] ... x[num_faces - 1] : binary labels of all the input faces

	Logger::out("-") << "formulating binary program...." << std::endl;
	w.start();

	program_.clear();

	const std::vector<Variable*>& variables = program_.create_n_variables(model_->size_of_facets());
	for (std::size_t i = 0; i < variables.size(); ++i) {
		Variable* v = variables[i];
		v->set_variable_type(Variable::BINARY);
	}

	LinearObjective* objective = program_.create_objective(LinearObjective::MINIMIZE);
	FOR_EACH_FACET(Map, model_, it) {
		Map::Facet* f = it;
		std::size_t var_idx = facet_indices[f];
		objective->add_coefficient(var_idx, 1.0);
	}

	//////////////////////////////////////////////////////////////////////////

	typedef typename HypothesisGenerator::SuperEdge SuperEdge;

	// Add constraints:
	for (std::size_t i = 0; i < adjacency.size(); ++i) {
		const SuperEdge& fan = adjacency[i];
		assert(fan.size() == 2);

		MapTypes::Halfedge* h0 = fan[0];
		MapTypes::Halfedge* h1 = fan[1];

		MapTypes::Facet* f0 = h0->facet();
		MapTypes::Facet* f1 = h1->facet();
		std::size_t var_idx0 = facet_indices[f0];
		std::size_t var_idx1 = facet_indices[f1];

		if (dot(Geom::vector(h0), Geom::vector(h1)) > 0) { // one must flip: x_i + x_j = 1
			LinearConstraint* c = program_.create_constraint(LinearConstraint::FIXED, 1.0, 1.0);
			c->add_coefficient(var_idx0, 1.0);
			c->add_coefficient(var_idx1, 1.0);
		}
		else { // both flip, or both not: x_i - x_j = 0
			LinearConstraint* c = program_.create_constraint(LinearConstraint::FIXED, 0.0, 0.0);
			c->add_coefficient(var_idx0, 1.0);
			c->add_coefficient(var_idx1, -1.0);
		}
	}

	Logger::out("-") << "#total variables: " << program_.variables().size() << std::endl;
	Logger::out("-") << "#total constraints: " << program_.constraints().size() << std::endl;
	Logger::out("-") << "formulating binary program done. " << w.elapsed() << " sec" << std::endl;

	//////////////////////////////////////////////////////////////////////////

	// Optimize model
	Logger::out("-") << "solving the binary program. Please wait..." << std::endl;
	w.start();

	LinearProgramSolver solver;
	if (solver.solve(&program_, solver_name)) {
		Logger::out("-") << "solving the binary program done. " << w.elapsed() << " sec" << std::endl;

		MapFacetAttribute<bool> visited(model_);
		FOR_EACH_FACET(Map, model_, it) {
			Map::Facet* f = it;
			visited[f] = false;
		}

		const std::vector<double>& X = solver.solution();
		MapEditor editor(model_);
		FOR_EACH_FACET(Map, model_, it) {
			Map::Facet* f = it;
			std::size_t fid = facet_indices[f];
			if (static_cast<int>(std::round(X[fid])) == 1 && !visited[f]) {
				editor.reorient_facet(f->halfedge());
				visited[f] = true;
			}
		}
		// Note: A border edge is now parallel to its opposite edge.
		// We scan all border edges for this property. If it holds, we
		// reorient the associated hole and search again until no border
		// edge with that property exists any longer. Then, all holes are
		// reoriented.
		FOR_EACH_HALFEDGE(Map, model_, it) {
			if (it->is_border() && it->vertex() == it->opposite()->vertex()) {
				editor.reorient_facet(it);
			}
		}
		model_->compute_facet_normals();
	}
	else {
		Logger::out("-") << "solving the binary program failed. " << w.elapsed() << " sec." << std::endl;
	}
}
