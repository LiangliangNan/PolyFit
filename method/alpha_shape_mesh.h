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

#pragma once

#include "method_common.h"
#include "../math/math_types.h"
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