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

#ifndef _POINT_SET_H_
#define _POINT_SET_H_


#include <model/model_common.h>
#include <basic/basic_types.h>
#include <math/math_types.h>
#include <basic/counted.h>
#include <basic/smart_pointer.h>

#include <model/vertex_group.h>
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
