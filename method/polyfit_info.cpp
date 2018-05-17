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

#include "polyfit_info.h"
#include "method_global.h"
#include "alpha_shape_mesh.h"
#include "../basic/logger.h"
#include "../basic/progress.h"
#include "../basic/stop_watch.h"
#include "../model/map_geometry.h"
#include "../model/point_set.h"
#include "../model/map.h"
#include "../model/map_circulators.h"
#include "../model/kdtree_search.h"

#include <vector>
#include <algorithm>


void PolyFitInfo::clear() {
	for (std::size_t i = 0; i < planes.size(); ++i)
		delete planes[i];
	planes.clear();
}


bool PolyFitInfo::ready_for_optimization(Map* mesh) const {
	return (
		planes.size() > 4 &&
		MapFacetAttribute<double>::is_defined(mesh, Method::facet_attrib_supporting_point_num) &&
		MapFacetAttribute<double>::is_defined(mesh, Method::facet_attrib_facet_area) &&
		MapFacetAttribute<double>::is_defined(mesh, Method::facet_attrib_covered_area)
		);
}


float PolyFitInfo::compute_point_confidences(PointSet* pset, int s1 /* = 6 */, int s2 /* = 16 */, int s3 /* = 32 */) {
	std::vector<vec3>& points = pset->points();
	std::vector<float>& planar_qualities = pset->planar_qualities();

	if (planar_qualities.size() != points.size())
		planar_qualities.resize(points.size());

	KdTreeSearch_var kdtree = new KdTreeSearch;
	kdtree->begin();
	kdtree->add_vertex_set(pset);
	kdtree->end();

	std::vector<int> neighbor_size;
	neighbor_size.push_back(s1);
	neighbor_size.push_back(s2);
	neighbor_size.push_back(s3);

	std::vector< std::vector<double> > eigen_values(3);
	for (int i = 0; i < 3; ++i)
		eigen_values[i].resize(3);

	double total = 0;
	ProgressLogger progress(points.size());
	for (std::size_t i = 0; i < points.size(); ++i) {
		const vec3& p = points[i];
		for (int j = 0; j < 3; ++j) {
			std::vector<unsigned int> neighbors;
			std::vector<double> sqr_distances;
			kdtree->find_closest_K_points(p, neighbor_size[j], neighbors, sqr_distances);

			PrincipalAxes3d pca;
			pca.begin();
			double avg = 0;
			for (unsigned int k = 0; k < neighbors.size(); ++k) {
				pca.add_point(points[neighbors[k]]);
				if (j == 0) {
					avg += std::sqrt(sqr_distances[k]);
				}
			}
			pca.end();
			total += (avg / neighbors.size());

			for (int k = 0; k < 3; ++k)
				eigen_values[j][k] = pca.eigen_value(3 - k - 1); // eigen values are sorted in descending order

			assert(eigen_values[j][0] <= eigen_values[j][1] && eigen_values[j][1] <= eigen_values[j][2]);
		}

		double conf = 0.0;
		for (int j = 0; j < 3; ++j) {
			conf += (1 - 3.0 * eigen_values[j][0] / (eigen_values[j][0] + eigen_values[j][1] + eigen_values[j][2])) * (eigen_values[j][1] / eigen_values[j][2]);
		}
		conf /= 3.0;
		planar_qualities[i] = static_cast<float>(conf);
		progress.next();
	}
	return static_cast<float>(total / points.size());
}


void PolyFitInfo::generate(PointSet* pset, Map* mesh, bool use_conficence /* = false */) {
	facet_attrib_supporting_vertex_group.bind(mesh, Method::facet_attrib_supporting_vertex_group);

	StopWatch w;
	Logger::out("-") << "computing point confidences..." << std::endl;
	double avg_spacing = compute_point_confidences(pset, 6, 16, 25);
	float radius = static_cast<float>(avg_spacing)* 5.0f;
	Logger::out("-") << "done. avg spacing: " << avg_spacing << ". " << w.elapsed() << " sec." << std::endl;

	std::vector<VertexGroup::Ptr>& groups = pset->groups();
	const std::vector<vec3>& pts = pset->points();

	float max_dist = 0;
	for (std::size_t i = 0; i < groups.size(); ++i) {
		VertexGroup* g = groups[i];
		const Plane3d& plane = g->plane();
		for (std::size_t j = 0; j < g->size(); ++j) {
			int idx = g->at(j);
			const vec3& p = pts[idx];
			float sdist = plane.squared_ditance(p);
			max_dist = std::max(max_dist, sdist);
		}
	}
	max_dist = std::sqrt(max_dist);

	//////////////////////////////////////////////////////////////////////////

	MapFacetAttribute<double>	facet_attrib_supporting_point_num(mesh, Method::facet_attrib_supporting_point_num);
	MapFacetAttribute<double>	facet_attrib_facet_area(mesh, Method::facet_attrib_facet_area);
	MapFacetAttribute<double>	facet_attrib_covered_area(mesh, Method::facet_attrib_covered_area);

	Logger::out("-") << "computing face confidences..." << std::endl;
	w.start();
	ProgressLogger progress(mesh->size_of_facets());
	FOR_EACH_FACET(Map, mesh, it) {
		Map::Facet* f = it;

		double face_area = Geom::facet_area(f);
		if (face_area < 1e-16) {
			Logger::err("-") << "degenerate facet with area: " << face_area << std::endl;
			FacetHalfedgeCirculator cir(f);
			for (; !cir->end(); ++cir) {
				Logger::err("-") << cir->vertex()->point() << std::endl;
			}
			Logger::err("-") << std::endl;
			continue;
		}

		VertexGroup* g = facet_attrib_supporting_vertex_group[f];

		std::vector<unsigned int> points;
		double num = facet_points_projected_in(pset, g, f, max_dist, points);
		if (use_conficence)
			facet_attrib_supporting_point_num[f] = num;
		else
			facet_attrib_supporting_point_num[f] = static_cast<double>(points.size());

		facet_attrib_facet_area[f] = face_area;

		Map::Ptr alpha_mesh = AlphaShapeMesh::apply(pset, points, g->plane(), radius);
		double covered_area = 0;
		if (alpha_mesh) {
			FOR_EACH_FACET(Map, alpha_mesh, it) 
				covered_area += Geom::triangle_area(it);
		}
		facet_attrib_covered_area[f] = covered_area;

		if (covered_area > face_area) {
			// this may not be an error (floating point precision limit)
			facet_attrib_covered_area[f] = face_area;
		}
		progress.next();
	}

	facet_attrib_supporting_vertex_group.unbind();

	Logger::out("-") << "done. " << w.elapsed() << " sec." << std::endl;
}


float PolyFitInfo::facet_points_projected_in(PointSet* pset, VertexGroup* g, MapTypes::Facet* f, float max_dist, std::vector<unsigned int>& points) {
	const Plane3d& plane = g->plane();
	const vec3& orig = plane.point();
	const vec3& base1 = plane.base1();
	const vec3& base2 = plane.base2();

	const Polygon3d& plg3d = Geom::facet_polygon(f);
	const Polygon2d& plg2d = Geom::to_2d(orig, base1, base2, plg3d);
	const std::vector<vec3>& pts = pset->points();
	const std::vector<float>& confidences = pset->planar_qualities();

	points.clear();
	float epsilon = max_dist * 0.5f;// considering noise and outliers
	float count = 0.0f;
	for (int i = 0; i < g->size(); ++i) {
		unsigned int idx = g->at(i);
		const vec3& p = pts[idx];
		const vec2& q = Geom::to_2d(orig, base1, base2, p);
		if (Geom::point_is_in_polygon(plg2d, q)) {
			points.push_back(idx);
			float dist = std::sqrt(plane.squared_ditance(p));
			if (dist < epsilon) { // in case of numerical issues (floating point precision)
				count += (1 - dist / epsilon) * confidences[idx];
			}
		}
	}

	return count;
}

