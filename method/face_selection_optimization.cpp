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

#include "face_selection.h"
#include "polyfit_info.h"
#include "method_global.h"
#include "../model/point_set.h"
#include "../model/map_editor.h"
#include "../model/map_geometry.h"
#include "../basic/stop_watch.h"
#include "../basic/logger.h"
#include "../math/linear_program_solver.h"


void FaceSelection::optimize(PolyFitInfo* polyfit_info) {
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

	Logger::out("-") << "extracting face adjacencies...." << std::endl;
	StopWatch w;
	const std::vector< FaceStar >& fans = adjacency_.extract(model_, polyfit_info->planes);
	Logger::out("-") << "done. " << w.elapsed() << " seconds." << std::endl;

	//-------------------------------------

	// binary variables:
	// x[0] ... x[num_faces - 1] : binary labels of all the input faces
	// x[num_faces] ... x[num_faces + num_edges] : binary labels of all the intersecting edges (remain or not)
	// x[num_faces + num_edges] ... x[num_faces + num_edges + num_edges] : binary labels of corner edges (sharp edge of not)

	Logger::out("-") << "formulating binary program...." << std::endl;
	w.start();

	std::size_t num_faces = model_->size_of_facets();
	std::size_t num_edges = 0;
	std::map<const FaceStar*, std::size_t> edge_usage_status;	// keep or remove an intersecting edges
	for (std::size_t i = 0; i < fans.size(); ++i) {
		const FaceStar& fan = fans[i];
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
	double coeff_complexity = total_points * Method::lambda_model_complexity / double(fans.size());

	typedef Variable<double>			Variable;
	typedef LinearExpression<double>	Objective;
	typedef LinearConstraint<double>	Constraint;
	typedef LinearProgram<double>		LinearProgram;

	Objective obj;

	std::map<const FaceStar*, std::size_t> edge_sharp_status;	// the edge is sharp or not
	std::size_t num_sharp_edges = 0;
	for (std::size_t i = 0; i < fans.size(); ++i) {
		const FaceStar& fan = fans[i];
		if (fan.size() == 4) {
			std::size_t var_idx = num_faces + num_edges + num_sharp_edges;
			edge_sharp_status[&fan] = var_idx;

			// accumulate model complexity term
			obj.add_coefficient(var_idx, coeff_complexity);
			++num_sharp_edges;
		}
	}
	assert(num_edges == num_sharp_edges);

	FOR_EACH_FACET(Map, model_, it) {
		Map::Facet* f = it;
		std::size_t var_idx = facet_indices[f];

		// accumulate data fitting term
		double num = facet_attrib_supporting_point_num_[f];
		obj.add_coefficient(var_idx, -coeff_data_fitting * num);

		// accumulate model coverage term
		double uncovered_area = (facet_attrib_facet_area_[f] - facet_attrib_covered_area_[f]);
		obj.add_coefficient(var_idx, coeff_coverage * uncovered_area);
	}
	program_.set_objective(obj, LinearProgram::MINIMIZE);

	std::size_t total_variables = num_faces + num_edges + num_sharp_edges;
	Logger::out("-") << "#total variables: " << total_variables << std::endl;
	Logger::out(" ") << "    - face is selected: " << num_faces << std::endl;
	Logger::out(" ") << "    - edge is used: " << num_edges << std::endl;
	Logger::out(" ") << "    - edge is sharp: " << num_sharp_edges << std::endl;

	typedef LinearProgram::Variable Variable;
	for (std::size_t i = 0; i < total_variables; ++i) {
		program_.add_variable(Variable(Variable::BINARY));
	}

	//////////////////////////////////////////////////////////////////////////

	typedef LinearProgram::Constraint Constraint;

	// Add constraints: the number of faces associated with an edge must be either 2 or 0
	std::size_t var_edge_used_idx = 0;
	for (std::size_t i = 0; i < fans.size(); ++i) {
		Constraint constraint;

		const FaceStar& fan = fans[i];
		for (std::size_t j = 0; j < fan.size(); ++j) {
			MapTypes::Facet* f = fan[j]->facet();
			std::size_t var_idx = facet_indices[f];
			constraint.add_coefficient(var_idx, 1.0);
		}

		if (fan.size() == 4) {
			std::size_t var_idx = num_faces + var_edge_used_idx;
			constraint.add_coefficient(var_idx, -2.0);  // 
			++var_edge_used_idx;
		}
		else { // boundary edge
			   // will be set to 0 (i.e., we don't allow open surface)
		}

		constraint.set_bounds(Constraint::FIXED, 0.0, 0.0);
		program_.add_constraint(constraint);
	}

	// Add constraints: for the sharp edges
	double M = 1.0;
	for (std::size_t i = 0; i < fans.size(); ++i) {
		const FaceStar& fan = fans[i];
		if (fan.size() != 4)
			continue;

		Constraint edge_used_constraint;
		std::size_t var_edge_usage_idx = edge_usage_status[&fan];
		edge_used_constraint.add_coefficient(var_edge_usage_idx, 1.0);
		std::size_t var_edge_sharp_idx = edge_sharp_status[&fan];
		edge_used_constraint.add_coefficient(var_edge_sharp_idx, -1.0);
		edge_used_constraint.set_bounds(Constraint::LOWER, 0.0, 0.0);
		program_.add_constraint(edge_used_constraint);

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

					Constraint edge_sharp_constraint;
					edge_sharp_constraint.add_coefficient(var_edge_sharp_idx, 1.0);
					edge_sharp_constraint.add_coefficient(fid1, -M);
					edge_sharp_constraint.add_coefficient(fid2, -M);
					edge_sharp_constraint.add_coefficient(var_edge_usage_idx, -M);
					edge_sharp_constraint.set_bounds(Constraint::LOWER, 1.0 - 3.0 * M, 0.0);
					program_.add_constraint(edge_sharp_constraint);
				}
			}
		}
	}

	Logger::out("-") << "#total constraints: " << program_.constraints().size() << std::endl;
	Logger::out("-") << "formulating binary program done. " << w.elapsed() << " sec" << std::endl;

	//////////////////////////////////////////////////////////////////////////

	// Optimize model
	Logger::out("-") << "solving the binary program...." << std::endl;
	w.start();

	LinearProgramSolver solver;
	if (solver.solve(&program_, Method::solver_name)) {
		Logger::out("-") << "solving the binary program done. " << w.elapsed() << " sec" << std::endl;

		// mark results
		const std::vector<double>& X = solver.get_result();
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

		for (std::size_t i = 0; i < fans.size(); ++i) {
			const FaceStar& fan = fans[i];
			if (fan.size() != 4)
				continue;

			std::size_t idx_sharp_var = edge_sharp_status[&fan];
			if (static_cast<int>(X[idx_sharp_var]) == 1) {
				for (std::size_t j = 0; j < fan.size(); ++j) {
					Map::Halfedge* e = fan[j];
					Map::Facet* f = e->facet();
					if (f) { // some faces may be deleted
						std::size_t fid = facet_indices[f];
						if (static_cast<int>(X[fid]) == 1) {
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
