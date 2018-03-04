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

#include "../3rd_lpsolve/lp_lib.h"


// PolyFit is always linked against both lpsolve and Gurobi.
// Let the user chose the solver.

//#ifdef _DEBUG
//#pragma comment(lib, "../x64/Debug/3rd_lpsolve.lib")
//#else
//#pragma comment(lib, "../x64/Release/3rd_lpsolve.lib")
//#endif


int FaceSelection::num_of_constriants(const std::vector<FaceStar>& fans) {
	int num = static_cast<int>(fans.size());

	for (std::size_t i = 0; i < fans.size(); ++i) {
		const FaceStar& fan = fans[i];
		std::size_t count = fan.size();
		if (count == 4) {
			++num;

			for (int j = 0; j < count; ++j) {
				MapTypes::Facet* f1 = fan[j]->facet();
				Plane3d* plane1 = face_attrib_supporting_plane_[f1];
				for (int k = j + 1; k < count; ++k) {
					MapTypes::Facet* f2 = fan[k]->facet();
					Plane3d* plane2 = face_attrib_supporting_plane_[f2];

					if (plane1 != plane2) {
						++num;
					}
				}
			}
		}
	}
	return num;
}


void FaceSelection::optimize_lp_solve(PolyFitInfo* polyfit_info, bool prune_faces) {
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

	StopWatch w;
	Logger::out("-") << "extracting adjacencies...." << std::endl;
	if (prune_faces)
		prune_boundary_faces(model_);
	const std::vector< FaceStar >& fans = adjacency_.extract(model_, polyfit_info->planes);
	Logger::out("-") << "done. " << w.elapsed() << " seconds." << std::endl;

	double total_points = double(pset_->points().size());
	int idx = 0;
	MapFacetAttribute<int>	facet_indices(model_);
	std::vector<Map::Facet*> all_faces;
	FOR_EACH_FACET(Map, model_, it) {
		Map::Facet* f = it;
		all_faces.push_back(f);
		facet_indices[f] = idx;
		++idx;
	}

	//////////////////////////////////////////////////////////////////////////

	Logger::out("-") << "optimization using LP_SOLVE solver..." << std::endl;

	int num_faces = model_->size_of_facets();

	std::map<const FaceStar*, int> edge_usage_status, edge_sharp_status;

	// the variables representing the binary labels of intersecting edges (remain or not)
	int index_remained_edge = 0;
	for (std::size_t i = 0; i < fans.size(); ++i) {
		const FaceStar& fan = fans[i];
		if (fan.size() == 4) {
			edge_usage_status[&fan] = num_faces + index_remained_edge;
			++index_remained_edge;
		}
	}

	int index_sharp_edge = 0;
	for (std::size_t i = 0; i < fans.size(); ++i) {
		const FaceStar& fan = fans[i];
		if (fan.size() == 4) {
			edge_sharp_status[&fan] = num_faces + index_remained_edge + index_sharp_edge;
			++index_sharp_edge;
		}
	}

	int num_variables = num_faces + index_remained_edge + index_sharp_edge;

	/* Create a new LP model */
	lprec* lp = make_lp(0, num_variables);
	if (!lp) {
		Logger::err("-") << "error in creating a LP model" << std::endl;
		return;
	}

	set_verbose(lp, SEVERE);
	for (int i = 1; i <= num_variables; ++i) {
		set_binary(lp, i, TRUE); /* sets variable i to binary */
	}

	//////////////////////////////////////////////////////////////////////////
	// Set objective

	double coeff_data_fitting = Method::lambda_data_fitting / total_points;
	double coeff_coverage = Method::lambda_model_coverage / model_->bbox().area();
	double coeff_complexity = Method::lambda_model_complexity / double(fans.size());

	// the variables representing the binary labels of corner edges
	std::vector<double> row(num_variables + 1, 0);	// lp_solve uses 1-based arrays
	// The LP_SOLVE manual says the first element is ignored, so any value doesn't effect.
	//row[0] = Method::lambda_data_fitting;

	int var_complexity_offset = num_faces + index_remained_edge;
	for (std::size_t i = 0; i < fans.size(); ++i) {
		const FaceStar& fan = fans[i];
		if (fan.size() == 4) {
			row[var_complexity_offset + 1] += coeff_complexity;	// lp_solve uses 1-based arrays
			++var_complexity_offset;
		}
	}

	FOR_EACH_FACET(Map, model_, it) {
		Map::Facet* f = it;
		int idx = facet_indices[f];
		double num = facet_attrib_supporting_point_num_[f];
		row[idx + 1] -= coeff_data_fitting * num;		// lp_solve uses 1-based arrays

		double uncovered_area = facet_attrib_facet_area_[f] - facet_attrib_covered_area_[f];
		row[idx + 1] += coeff_coverage * uncovered_area;// lp_solve uses 1-based arrays
	}

	set_obj_fn(lp, row.data());

	Logger::out("-") << "   #variables for face selection: " << num_faces << std::endl;
	Logger::out("-") << "   #variables for edge usage: " << index_remained_edge << std::endl;
	Logger::out("-") << "   #variables for sharp edges: " << index_sharp_edge << std::endl;
	Logger::out("-") << "   #total binary variables: " << num_variables << std::endl;

	//////////////////////////////////////////////////////////////////////////

	int num = num_of_constriants(fans);
	Logger::out("-") << "   #total constraints: " << num << std::endl;

	/* num constraints will be added, so allocate memory for it in advance to make things faster */
	resize_lp(lp, num, get_Ncolumns(lp));

	set_add_rowmode(lp, TRUE);

	// Add constraints: the number of faces associated with an edge must be either 2 or 0
	int constraint_idx = 0;
	for (std::size_t i = 0; i < fans.size(); ++i) {
		const FaceStar& fan = fans[i];
		std::size_t count = fan.size();

		std::vector<int> colno;
		std::vector<double> sparserow;
		for (std::size_t j = 0; j < count; ++j) {
			MapTypes::Facet* f = fan[j]->facet();
			int fid = facet_indices[f];
			colno.push_back(fid + 1);	 // lp_solve uses 1-based arrays
			sparserow.push_back(1.0);
		}

		if (count >= 2) { //model.addConstr(sum == 2 * X[num_faces + constraint_idx]);
			colno.push_back(num_faces + constraint_idx + 1); // lp_solve uses 1-based arrays
			sparserow.push_back(-2.0);
			++constraint_idx;
		}
		else { //model.addConstr(num_facet == 0);	// make sure the surface is closed (no open edge)
			// no need
		}
		add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), EQ, 0.0);
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

		if (count == 4) {
			int colno[2];
			REAL sparserow[2]; 
			colno[0] = idx_edge_var + 1;	sparserow[0] = 1.0;	 // lp_solve uses 1-based arrays
			colno[1] = idx_sharp_var + 1;	sparserow[1] = -1.0; // lp_solve uses 1-based arrays
			//model.addConstr(X[idx_edge_var] >= X[idx_corner_var]);
			add_constraintex(lp, 2, sparserow, colno, GE, 0.0); 

			for (int j = 0; j < count; ++j) {
				MapTypes::Facet* f1 = fan[j]->facet();
				Plane3d* plane1 = face_attrib_supporting_plane_[f1];
				int fid1 = facet_indices[f1];
				for (int k = j + 1; k < count; ++k) {
					MapTypes::Facet* f2 = fan[k]->facet();
					Plane3d* plane2 = face_attrib_supporting_plane_[f2];
					int fid2 = facet_indices[f2];

					if (plane1 != plane2) {
						//model.addConstr(X[idx_corner_var] + M * (3 - (X[fid1] + X[fid2] + X[idx_edge_var])) >= 1);
						//// equals to  
						//model.addConstr(X[idx_corner_var] - M * X[fid1] - M * X[fid2] - M * X[idx_edge_var] >= 1 - 3M);
						int col_no[4];
						REAL sparse_row[4];
						col_no[0] = idx_sharp_var + 1;	sparse_row[0] = 1.0; // lp_solve uses 1-based arrays
						col_no[1] = fid1 + 1;			sparse_row[1] = -M;	 // lp_solve uses 1-based arrays
						col_no[2] = fid2 + 1;			sparse_row[2] = -M;	 // lp_solve uses 1-based arrays
						col_no[3] = idx_edge_var + 1;	sparse_row[3] = -M;	 // lp_solve uses 1-based arrays
						add_constraintex(lp, 4, sparse_row, col_no, GE, 1 - 3 * M);
					}
				}
			}
		}
	}
	set_add_rowmode(lp, FALSE);

	////////////////////////////////////////////////////////////////////////////
	// Optimize model
	w.start();
	Logger::out("-") << "minimizing energy..." << std::endl;
	int ret = solve(lp);

	switch (ret) {
	case 0:
	{
		double objval = get_objective(lp);
		Logger::out("-") << "done. objective: " << Method::lambda_data_fitting + truncate_digits(objval, 3) << ". " << w.elapsed() << " sec." << std::endl;

		std::vector<double> results(num_variables);
		get_variables(lp, results.data());

		std::vector<Map::Facet*> to_delete;
		for (std::size_t i = 0; i < num_faces; ++i) {
			int value = static_cast<int>(results[i]);
			if (value == 0) {
				to_delete.push_back(all_faces[i]);
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

			int idx_sharp_var = edge_sharp_status[&fan];
			int status = static_cast<int>(results[idx_sharp_var]);
			if (status == 1) {
				for (int j = 0; j < fan.size(); ++j) {
					Map::Halfedge* e = fan[j];
					Map::Facet* f = e->facet();
					if (f) { // some faces may be deleted
						int fid = facet_indices[f];
						int face_status = static_cast<int>(results[fid]);
						if (face_status == 1) {
							edge_is_sharp[e] = true;
							break;
						}
					}
				}
			}
		}
	}
	break;

	case -2:
		Logger::err("-") << "Out of memory" << std::endl;
		break;
	case 1:
		Logger::err("-") << "The model is sub-optimal. Only happens if there are integer variables and there is already an integer solution found. The solution is not guaranteed the most optimal one." << std::endl;
		break;
	case 2:
		Logger::err("-") << "The model is infeasible" << std::endl;
		break;
	case 3:
		Logger::err("-") << "The model is unbounded" << std::endl;
		break;
	case 4:
		Logger::err("-") << "The model is degenerative" << std::endl;
		break;
	case 5:
		Logger::err("-") << "Numerical failure encountered" << std::endl;
		break;
	case 6:
		Logger::err("-") << "The abort() routine was called" << std::endl;
		break;
	case 7:
		Logger::err("-") << "A timeout occurred" << std::endl;
		break;
	case 9:
		Logger::err("-") << "The model could be solved by presolve. This can only happen if presolve is active via set_presolve()" << std::endl;
		break;
	case 25:
		Logger::err("-") << "Accuracy error encountered" << std::endl;
		break;
	default:
		break;
	}

	delete_lp(lp);

	facet_attrib_supporting_vertex_group_.unbind();
	facet_attrib_supporting_point_num_.unbind();
	facet_attrib_facet_area_.unbind();
	facet_attrib_covered_area_.unbind();

	vertex_source_planes_.unbind();
	edge_source_planes_.unbind();
	face_attrib_supporting_plane_.unbind();
}
