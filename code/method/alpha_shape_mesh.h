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

#pragma once

#include <method/method_common.h>
#include <math/math_types.h>
#include <string>
#include <vector>

class Map;
class PointSet;
class VertexGroup;
class AlphaShape;

// convert the AlphaShape into a mesh representation
class METHOD_API AlphaShapeMesh
{
public:
	// plane is used to lift the 2D points into 3D space
	static Map* apply(AlphaShape* as, const Plane3d& plane, float radius);

	// the input is a vertex group
	static Map* apply(const VertexGroup* g, float radius);

	// the input is a subset of a point cloud and the points lie on plane
	static Map* apply(const PointSet* pset, const std::vector<unsigned int>& point_indices, const Plane3d& plane, float radius);

};