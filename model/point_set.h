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

#ifndef _POINT_SET_H_
#define _POINT_SET_H_


#include "model_common.h"
#include "../basic/basic_types.h"
#include "../math/math_types.h"
#include "../basic/counted.h"
#include "../basic/smart_pointer.h"

#include "vertex_group.h"
#include <list>


class VertexGroup;

class MODEL_API PointSet : public Counted
{
public:
	typedef SmartPointer<PointSet>	Ptr;

public:
	PointSet();
	~PointSet();

    unsigned int  num_points() const { return static_cast<unsigned int>(points_.size()); }

	std::vector<vec3>& points() { return points_; }
	std::vector<vec3>& colors() { return colors_; }
	std::vector<vec3>& normals() { return normals_; }
	std::vector<float>& planar_qualities() { return planar_qualities_; }
	const std::vector<vec3>& points() const { return points_; }
	const std::vector<vec3>& colors() const { return colors_; }
	const std::vector<vec3>& normals() const { return normals_; }
	const std::vector<float>& planar_qualities() const { return planar_qualities_; }

	bool    has_normals() const { return normals_.size() > 0 && normals_.size() == points_.size(); }
	bool	has_colors() const  { return colors_.size() > 0 && colors_.size() == points_.size(); }
	bool    has_planar_qualities() const { return planar_qualities_.size() > 0 && planar_qualities_.size() == points_.size(); }
	
	void	delete_points(const std::vector<unsigned int>& indices);

	//////////////////////////////////////////////////////////////////////////

	std::vector<VertexGroup::Ptr>& groups() { return groups_; }
	const std::vector<VertexGroup::Ptr>& groups() const { return groups_; }

	// the points that don't belong to any vertex groups
	std::vector<unsigned int> idle_points() const;

	void fit_plane(VertexGroup::Ptr g);

	const Box3d& bbox() const;
	void invalidate_bbox() { bbox_is_valid_ = false; }

private:
	std::vector<vec3>  points_;
	std::vector<vec3>  colors_;
	std::vector<vec3>  normals_;
	std::vector<float> planar_qualities_;

	mutable bool	bbox_is_valid_;
	mutable Box3d	bbox_;

	std::vector<VertexGroup::Ptr>	groups_;

};


#endif // _POINT_SET_H_
