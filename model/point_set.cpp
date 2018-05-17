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

#include "point_set.h"
#include "vertex_group.h"


PointSet::PointSet() : bbox_is_valid_(false)
{
}

PointSet::~PointSet() {
}

const Box3d& PointSet::bbox() const {
	if (!bbox_is_valid_) {
		Box3d result;
		for (int i = 0; i < points_.size(); ++i) {
			result.add_point(points_[i]);
		}
		bbox_ = result;
		bbox_is_valid_ = true;
	}
	return bbox_;
}

void PointSet::delete_points(const std::vector<unsigned int>& indices) {
	std::vector<int> keep(num_points(), 1);
	for (std::size_t i = 0; i < indices.size(); ++i) {
		unsigned int id = indices[i];
		keep[id] = 0;
	}
	
	std::vector<vec3>  new_points;
	std::vector<vec3>  new_colors;
	std::vector<vec3>  new_normals;
	std::vector<float>  new_planar_qualities;
	for (std::size_t i = 0; i < num_points(); ++i) {
		if (keep[i] == 1) 
			new_points.push_back(points_[i]);
	}
	if (has_colors()) {
		for (std::size_t i = 0; i < num_points(); ++i) {
			if (keep[i] == 1)
				new_colors.push_back(colors_[i]);
		}
	}
	if (has_normals()) {
		for (std::size_t i = 0; i < num_points(); ++i) {
			if (keep[i] == 1)
				new_normals.push_back(normals_[i]);
		}
	}
	if (has_planar_qualities()) {
		for (std::size_t i = 0; i < num_points(); ++i) {
			if (keep[i] == 1)
				new_planar_qualities.push_back(planar_qualities_[i]);
		}
	}

	points_ = new_points;
	normals_ = new_normals;
	colors_ = new_colors;
	planar_qualities_ = new_planar_qualities;
}

std::vector<unsigned int> PointSet::idle_points() const {
	std::vector<int> remained(num_points(), 1);
	for (std::size_t i = 0; i < groups_.size(); ++i) {
		VertexGroup* g = groups_[i];
		for (std::size_t j = 0; j < g->size(); ++j) {
			unsigned int id = g->at(j);
			remained[id] = 0;
		}
	}

	std::vector<unsigned int> results;
	for (std::size_t i = 0; i < num_points(); ++i) {
		if (remained[i])
            results.push_back(static_cast<unsigned int>(i));
	}
	return results;
}


void PointSet::fit_plane(VertexGroup::Ptr g) {
	PrincipalAxes3d pca;
	pca.begin();
	for (std::size_t j = 0; j < g->size(); ++j) {
		pca.add_point(points_[g->at(j)]);
	}
	pca.end();

	Plane3d plane(pca.center(), pca.axis(2)); // eigen vectors have been normalized
	g->set_plane(plane);
}
