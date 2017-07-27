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


#include "alpha_shape_boundary.h"
#include "../math/math_types.h"
#include "../model/vertex_group.h"
#include "../basic/stop_watch.h"
#include "../basic/logger.h"

#include "alpha_shape.h"




std::vector<unsigned int> AlphaShapeBoundary::apply(VertexGroup* g, int percent/* = 99*/) {
	std::vector<unsigned int> result;

	PointSet* pset = g->point_set();
	if (!pset)
		return result;

	std::size_t num_input = g->size();
	const std::vector<vec3>& points = pset->points();

	std::list<Point2> pts;
	Plane3d plane = g->plane();
	for (std::size_t i = 0; i < num_input; ++i) {
		int idx = g->at(i);
		const vec3& p = points[idx];
		const vec2& q = plane.to_2d(p); 
		const Point2& qq = to_cgal_point(q);
		pts.push_back(qq);
	}

 	AlphaShape as(pts.begin(), pts.end());  

	//////////////////////////////////////////////////////////////////////////

	std::set<unsigned int> boundaries;
	std::vector<Face_handle> faces = as.extract_largest_solid_component();
	for (unsigned int i=0; i<faces.size(); ++i) {
		Face_handle f = faces[i];
		for (unsigned int j=0; j<3; ++j) {
			Vertex_handle vh = f->vertex(j);
			const Point2& p = vh->point();
			if (as.classify(p) == AlphaShape::REGULAR) {
				int idx = vh->index();	// this is the index in the current vertex group
				idx = g->at(idx);		// this is the index in the original point cloud
				boundaries.insert(idx);
			}
		}
	}

	result.insert(result.end(), boundaries.begin(), boundaries.end());
	return result;
}

