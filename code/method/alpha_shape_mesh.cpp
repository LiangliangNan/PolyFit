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

#include <method/alpha_shape_mesh.h>
#include <math/math_types.h>
#include <model/vertex_group.h>
#include <basic/stop_watch.h>
#include <basic/logger.h>
#include <model/map.h>
#include <model/map_builder.h>
#include <method/alpha_shape.h>


// return the indices of boundary points
Map* AlphaShapeMesh::apply(AlphaShape* as, const Plane3d& plane, float radius) {
    double alpha_value = radius * radius;
	as->set_alpha(alpha_value);

	std::map<int, int> new_index;
	std::vector<vec3>  points;
	typedef std::vector<int> Triangle;
	std::vector<Triangle> faces;

	int current_idx = 0;
	for (AlphaShape::Finite_faces_iterator fit = as->finite_faces_begin(); fit != as->finite_faces_end(); ++fit) {
		if (as->classify(fit) == AlphaShape::INTERIOR) {
			Triangle tri;
			for (int i = 0; i < 3; ++i) {
				AlphaShape::Vertex_handle vh = fit->vertex(i);
				int old_idx = vh->index();
				if (new_index.find(old_idx) == new_index.end()) {
					const vec3& p = plane.to_3d(to_my_point(vh->point()));
					points.push_back(p);
					new_index[old_idx] = current_idx;
					tri.push_back(current_idx);
					++current_idx;
				}
				else {
					int idx = new_index[old_idx];
					tri.push_back(idx);
				}
			}

			faces.push_back(tri);
		}
	}

	if (faces.empty())
        return nullptr;

	Map* mesh = new Map;
	MapBuilder builder(mesh);
	builder.set_quiet(true);
	builder.begin_surface();

	for (std::size_t i = 0; i < points.size(); ++i)
		builder.add_vertex(points[i]);

	for (std::size_t i = 0; i < faces.size(); ++i) {
		const Triangle& tri = faces[i];
		builder.begin_facet();
        for (std::size_t j = 0; j < 3; ++j) {
            int idx = tri[j];
			builder.add_vertex_to_facet(idx);
		}
		builder.end_facet();
	}

	builder.end_surface();

	return mesh;
}



Map* AlphaShapeMesh::apply(const VertexGroup* g, float radius) {
	const PointSet* pset = g->point_set();
	if (!pset)
        return nullptr;

	std::size_t num_input = g->size();
	const std::vector<vec3>& points = pset->points();

	std::list<Point2> pts;
	const Plane3d& plane = g->plane();
	for (std::size_t i = 0; i < num_input; ++i) {
        unsigned int idx = g->at(i);
		const vec3& p = points[idx];
		const vec2& q = plane.to_2d(p);
		const Point2& qq = to_cgal_point(q);
		pts.push_back(qq);
	}

	AlphaShape as(pts.begin(), pts.end());

	return apply(&as, plane, radius);
}


Map* AlphaShapeMesh::apply(const PointSet* pset, const std::vector<unsigned int>& point_indices, const Plane3d& plane, float radius) {
	if (point_indices.size() < 10) {
		//Logger::out("-") << "very few points - no need to compute AlphaShapeMesh" << std::endl;
        return nullptr;
	}

	const std::vector<vec3>& points = pset->points();

	std::list<Point2> pts;
	for (std::size_t i = 0; i < point_indices.size(); ++i) {
        unsigned int idx = point_indices[i];
		const vec3& p = points[idx];
		const vec2& q = plane.to_2d(p);
		const Point2& qq = to_cgal_point(q);
		pts.push_back(qq);
	}

	AlphaShape as(pts.begin(), pts.end());
	return apply(&as, plane, radius);
}
