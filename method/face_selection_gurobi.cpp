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

#include <gurobi_c++.h>


#if (_MSC_VER == 1800)  // vs2013
#pragma comment(lib, "gurobi70.lib")
#ifdef _DEBUG
#pragma comment(lib, "gurobi_c++mdd2013.lib")
#else
#pragma comment(lib, "gurobi_c++md2013.lib")
#endif
#elif (_MSC_VER == 1900) // vs 2015
#pragma comment(lib, "gurobi75.lib")
#ifdef _DEBUG
#pragma comment(lib, "gurobi_c++mdd2015.lib")
#else
#pragma comment(lib, "gurobi_c++md2015.lib")
#endif
#elif (_MSC_VER == 1910 || _MSC_VER == 1911) // vs 2017
#pragma comment(lib, "gurobi75.lib")
#ifdef _DEBUG
#pragma comment(lib, "gurobi_c++mdd2017.lib")
#else
#pragma comment(lib, "gurobi_c++md2017.lib")
#endif
#endif



void FaceSelection::optimize_Gurobi(PolyFitInfo* polyfit_info, bool prune_faces) {
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
	face_attrib_supporting_plane_.bind_if_defined(model_, "FacetSupportingPlane");

	//////////////////////////////////////////////////////////////////////////

	Logger::out("-") << "optimization using GUTOBI solver..." << std::endl;

	try {
		StopWatch w;
		Logger::out("-") << "extracting adjacencies...." << std::endl;
		if (prune_faces)
			prune_boundary_faces(model_);
		const std::vector< FaceStar >& fans = adjacency_.extract(model_, polyfit_info->planes);
		Logger::out("-") << "done. " << w.elapsed() << " seconds." << std::endl;

		double total_points = double(pset_->points().size());
		int idx = 0;
		MapFacetAttribute<int>	facet_indices(model_);
		FOR_EACH_FACET(Map, model_, it) {
			Map::Facet* f = it;
			facet_indices[f] = idx;
			++idx;
		}

		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);

		GRBModel model = GRBModel(env);

		//////////////////////////////////////////////////////////////////////////
		// create variables

		int num_faces = model_->size_of_facets();
		std::vector<GRBVar> X;

		// the variables representing the binary labels of all the input faces
		for (int i = 0; i < num_faces; ++i) {						//
			GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			X.push_back(var);
		}

		std::map<const FaceStar*, int> edge_usage_status, edge_sharp_status;

		// the variables representing the binary labels of intersecting edges (remain or not)
		int num_edge_is_used = 0;
		for (std::size_t i = 0; i < fans.size(); ++i) {
			const FaceStar& fan = fans[i];
			if (fan.size() == 4) {
				GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				edge_usage_status[&fan] = static_cast<int>(X.size());
				X.push_back(var);
				++num_edge_is_used;
			}
		}

		GRBLinExpr term_data_fitting = 0.0;
		GRBLinExpr term_coverage = 0.0;
		GRBLinExpr term_complexity = 0.0;

		// the variables representing the binary labels of corner edges (sharp edge of not)
		int num_sharp_edge = 0;
		for (std::size_t i = 0; i < fans.size(); ++i) {
			const FaceStar& fan = fans[i];
			if (fan.size() == 4) {
				GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				edge_sharp_status[&fan] = static_cast<int>(X.size());	
				X.push_back(var);
				++num_sharp_edge;
				term_complexity += var;

			}
		}

		// Integrate new variables
		model.update();

		//////////////////////////////////////////////////////////////////////////
		// Set objective
		FOR_EACH_FACET(Map, model_, it) {
			Map::Facet* f = it;
			int idx = facet_indices[f];
			double num = facet_attrib_supporting_point_num_[f];
			term_data_fitting += X[idx] * num;
			double uncovered_area = (facet_attrib_facet_area_[f] - facet_attrib_covered_area_[f]);
			term_coverage += X[idx] * uncovered_area;	
		}

		term_data_fitting = 1.0 - term_data_fitting / total_points;
		term_coverage /= model_->bbox().area();
		term_complexity /= double(fans.size());

		GRBLinExpr obj =
			Method::lambda_data_fitting * term_data_fitting	+ 
			Method::lambda_model_coverage * term_coverage +
			Method::lambda_model_complexity * term_complexity;
		model.setObjective(obj, GRB_MINIMIZE);

		Logger::out("-") << "   #variables for face selection: " << num_faces << std::endl;
		Logger::out("-") << "   #variables for edge usage: " << edge_usage_status.size() << std::endl;
		Logger::out("-") << "   #variables for sharp edges: " << edge_sharp_status.size() << std::endl;
		Logger::out("-") << "   #total binary variables: " << X.size() << std::endl;

		//////////////////////////////////////////////////////////////////////////

		// Add constraints: the number of faces associated with an edge must be either 2 or 0
		std::size_t num = 0;
		int constraint_idx = 0;
		for (std::size_t i = 0; i < fans.size(); ++i) {
			const FaceStar& fan = fans[i];
			std::size_t count = fan.size();
			GRBLinExpr num_facet = 0;
			for (std::size_t j = 0; j < count; ++j) {
				MapTypes::Facet* f = fan[j]->facet();
				int fid = facet_indices[f];
				num_facet += X[fid];
			}
			if (count >= 2) {
				model.addConstr(num_facet == 2 * X[num_faces + constraint_idx]);
				++constraint_idx;
			}
			else {
				model.addConstr(num_facet == 0);	// make sure the surface is closed (no open edge)
			}

			++num;
		}

		// Add constraints: for the sharp edges
		int M = 3;
		for (std::size_t i = 0; i < fans.size(); ++i) {
			const FaceStar& fan = fans[i];
			std::size_t count = fan.size();
			if (count != 4) 
				continue;

			int idx_edge_var = edge_usage_status[&fan];
			int idx_sharp_var = edge_sharp_status[&fan];
			model.addConstr(X[idx_edge_var] >= X[idx_sharp_var]);
			++num;

			for (int j = 0; j < count; ++j) {
				MapTypes::Facet* f1 = fan[j]->facet();
				Plane3d* plane1 = face_attrib_supporting_plane_[f1];
				int fid1 = facet_indices[f1];
				for (int k = j + 1; k < count; ++k) {
					MapTypes::Facet* f2 = fan[k]->facet();
					Plane3d* plane2 = face_attrib_supporting_plane_[f2];
					int fid2 = facet_indices[f2];

					if (plane1 != plane2) {
						model.addConstr(X[idx_sharp_var] + M * (3 - (X[fid1] + X[fid2] + X[idx_edge_var])) >= 1);
						++num;
					}
				}
			}
		}

		Logger::out("-") << "   #total constraints: " << num << std::endl;

		//////////////////////////////////////////////////////////////////////////
		// Optimize model
		w.start();
		Logger::out("-") << "minimizing energy..." << std::endl;
		model.optimize();

		int optimstatus = model.get(GRB_IntAttr_Status);
		if (optimstatus == GRB_OPTIMAL) {
			double objval = model.get(GRB_DoubleAttr_ObjVal);
			Logger::out("-") << "done. objective: " << truncate_digits(objval, 3) << ". " << w.elapsed() << " sec." << std::endl;
		}
		else if (optimstatus == GRB_INF_OR_UNBD) {
			Logger::err("-") << "model is infeasible or unbounded" << std::endl;
			return;
		}
		else if (optimstatus == GRB_INFEASIBLE) {
			Logger::err("-") << "model is infeasible" << std::endl;
			return;
		}
		else if (optimstatus == GRB_UNBOUNDED) {
			Logger::err("-") << "model is unbounded" << std::endl;
			return;
		}
		else {
			Logger::err("-") << "optimization was stopped with status = " << optimstatus << std::endl;
			return;
		}

		//////////////////////////////////////////////////////////////////////////

		// mark results
		std::vector<Map::Facet*> to_delete;
		FOR_EACH_FACET(Map, model_, it) {
			Map::Facet* f = it;
			int idx = facet_indices[f];
			int value = static_cast<int>(X[idx].get(GRB_DoubleAttr_X));
			if (value == 0)
				to_delete.push_back(f);
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
			int status = static_cast<int>(X[idx_sharp_var].get(GRB_DoubleAttr_X));
			if (status == 1) {
				for (int j = 0; j < fan.size(); ++j) {
					Map::Halfedge* e = fan[j];
					Map::Facet* f = e->facet();
					if (f) { // some faces may be deleted
						int fid = facet_indices[f];
						int face_status = static_cast<int>(X[fid].get(GRB_DoubleAttr_X));
						if (face_status == 1) {
							edge_is_sharp[e] = true;
							break;
						}
					}
				}
			}
		}
	}
	catch (GRBException e) {
		Logger::err("-") << "Error code = " << e.getErrorCode() << std::endl;
		Logger::err("-") << e.getMessage() << std::endl;
	}
	catch (...) {
		Logger::err("-") << "Exception during optimization" << std::endl;
	}

	facet_attrib_supporting_vertex_group_.unbind();
	facet_attrib_supporting_point_num_.unbind();
	facet_attrib_facet_area_.unbind();
	facet_attrib_covered_area_.unbind();

	vertex_source_planes_.unbind();
	edge_source_planes_.unbind();
	face_attrib_supporting_plane_.unbind();
}