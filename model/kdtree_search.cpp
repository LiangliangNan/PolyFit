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

#include "kdtree_search.h"
#include "kdtree/kdTree.h"
#include "../model/point_set.h"



#define get_tree(x) ((kdtree::KdTree*)(x))


KdTreeSearch::KdTreeSearch()  {
	points_num_ = 0;
	tree_ = nil;
}


KdTreeSearch::~KdTreeSearch() {
    delete get_tree(tree_);
}


void KdTreeSearch::begin()  {
	vertices_.clear();

    delete get_tree(tree_);
	tree_ = nil;
}


void KdTreeSearch::end()  {
	points_num_ = vertices_.size();

	kdtree::Vector3D* points = new kdtree::Vector3D[points_num_];
	std::list<vec3*>::iterator it = vertices_.begin();
	unsigned int idx = 0;
	for (; it != vertices_.end(); ++it) {
		vec3* v = *it;
		points[idx].x = v->x;
		points[idx].y = v->y;
		points[idx].z = v->z;
		++idx;
	}

	unsigned int maxBucketSize = 16 ;	// number of points per bucket
	tree_ = new kdtree::KdTree(points, points_num_, maxBucketSize );
	delete [] points;
}


void KdTreeSearch::add_point(vec3* v)  {
	vertices_.push_back(v);
}


void KdTreeSearch::add_vertex_set(PointSet* vs)  {
	std::vector<vec3>& points = vs->points();
	for (std::size_t i = 0; i < points.size(); ++i)
		vertices_.push_back(&(points[i]));
}


int KdTreeSearch::find_closest_point(const vec3& p) const {
	kdtree::Vector3D v3d( p.x, p.y, p.z );
	get_tree(tree_)->setNOfNeighbours( 1 );
	get_tree(tree_)->queryPosition( v3d );

	unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
	if (num == 1) 
		return get_tree(tree_)->getNeighbourPositionIndex(0);
	else
		return -1;
}

int KdTreeSearch::find_closest_point(const vec3& p, double& squared_distance) const {
	kdtree::Vector3D v3d( p.x, p.y, p.z );
	get_tree(tree_)->setNOfNeighbours( 1 );
	get_tree(tree_)->queryPosition( v3d );

	unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
	if (num == 1) {
		squared_distance = get_tree(tree_)->getSquaredDistance(0);
		return get_tree(tree_)->getNeighbourPositionIndex(0);
	} else 
		return -1;
}

void KdTreeSearch::find_closest_K_points(
	const vec3& p, unsigned int k, std::vector<unsigned int>& neighbors
	)  const {
		kdtree::Vector3D v3d( p.x, p.y, p.z );
		get_tree(tree_)->setNOfNeighbours( k );
		get_tree(tree_)->queryPosition( v3d );

		unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
		if (num == k) {
			neighbors.resize(k);
			for (unsigned int i=0; i<k; ++i) {
				neighbors[i] = get_tree(tree_)->getNeighbourPositionIndex(i);
			}		
		} else
			std::cerr << "less than " << k << " points found" << std::endl;
}

void KdTreeSearch::find_closest_K_points(
	const vec3& p, unsigned int k, std::vector<unsigned int>& neighbors, std::vector<double>& squared_distances
	)  const {
		kdtree::Vector3D v3d( p.x, p.y, p.z );
		get_tree(tree_)->setNOfNeighbours( k );
		get_tree(tree_)->queryPosition( v3d );

		unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
		if (num == k) {
			neighbors.resize(k);
			squared_distances.resize(k);
			for (unsigned int i=0; i<k; ++i) {
				neighbors[i] = get_tree(tree_)->getNeighbourPositionIndex(i);
				squared_distances[i] = get_tree(tree_)->getSquaredDistance(i);
			}		
		} else
			std::cerr << "less than " << k << " points found" << std::endl;
}



void KdTreeSearch::find_points_in_radius(
	const vec3& p, double squared_radius, std::vector<unsigned int>& neighbors
	)  const {
		kdtree::Vector3D v3d( p.x, p.y, p.z );
		get_tree(tree_)->queryRange( v3d, squared_radius, true );

		unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
		neighbors.resize(num);
		for (unsigned int i=0; i<num; ++i) {
			neighbors[i] = get_tree(tree_)->getNeighbourPositionIndex(i);
		}	
}


void KdTreeSearch::find_points_in_radius(
	const vec3& p, double squared_radius, std::vector<unsigned int>& neighbors, std::vector<double>& squared_distances
	)  const {
		kdtree::Vector3D v3d( p.x, p.y, p.z );
		get_tree(tree_)->queryRange( v3d, squared_radius, true );

		unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
		neighbors.resize(num);
		squared_distances.resize(num);
		for (unsigned int i=0; i<num; ++i) {
			neighbors[i] = get_tree(tree_)->getNeighbourPositionIndex(i);
			squared_distances[i] = get_tree(tree_)->getSquaredDistance(i);
		}	
}


unsigned int KdTreeSearch::find_points_in_cylinder(
	const vec3& p1, const vec3& p2, double radius, 
	std::vector<unsigned int>& neighbors, std::vector<double>& squared_distances,
	bool bToLine
	) const {
		kdtree::Vector3D s( p1.x, p1.y, p1.z );
		kdtree::Vector3D t( p2.x, p2.y, p2.z );
		get_tree(tree_)->queryLineIntersection( s, t, radius, bToLine, true );

		unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();

		neighbors.resize(num);
		squared_distances.resize(num);
		for (unsigned int i=0; i<num; ++i) {
			neighbors[i] = get_tree(tree_)->getNeighbourPositionIndex(i);
			squared_distances[i] = get_tree(tree_)->getSquaredDistance(i);
		}	

		return num;
}

unsigned int KdTreeSearch::find_points_in_cylinder(
	const vec3& p1, const vec3& p2, double radius, 
	std::vector<unsigned int>& neighbors,
	bool bToLine
	) const {
		kdtree::Vector3D s( p1.x, p1.y, p1.z );
		kdtree::Vector3D t( p2.x, p2.y, p2.z );
		get_tree(tree_)->queryLineIntersection( s, t, radius, bToLine, true );

		unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
		neighbors.resize(num);
		for (unsigned int i=0; i<num; ++i) {
			neighbors[i] = get_tree(tree_)->getNeighbourPositionIndex(i);
		}

		return num;
}


unsigned int KdTreeSearch::find_points_in_cone(
	const vec3& eye, const vec3& p1, const vec3& p2, double angle_range, 
	std::vector<unsigned int>& neighbors, std::vector<double>& squared_distances,
	bool bToLine
	) const {
		kdtree::Vector3D eye3d( eye.x, eye.y, eye.z );
		kdtree::Vector3D s( p1.x, p1.y, p1.z );
		kdtree::Vector3D t( p2.x, p2.y, p2.z ); 
		get_tree(tree_)->queryConeIntersection( eye3d, s, t, angle_range, bToLine, true );

		unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
		neighbors.resize(num);
		squared_distances.resize(num);
		for (unsigned int i=0; i<num; ++i) {
			neighbors[i] = get_tree(tree_)->getNeighbourPositionIndex(i);
			squared_distances[i] = get_tree(tree_)->getSquaredDistance(i);
		}

		return num;
}

unsigned int KdTreeSearch::find_points_in_cone(
	const vec3& eye, const vec3& p1, const vec3& p2, double angle_range, 
	std::vector<unsigned int>& neighbors,
	bool bToLine
	) const {
		kdtree::Vector3D eye3d( eye.x, eye.y, eye.z );
		kdtree::Vector3D s( p1.x, p1.y, p1.z );
		kdtree::Vector3D t( p2.x, p2.y, p2.z );
		get_tree(tree_)->queryConeIntersection( eye3d, s, t, angle_range, bToLine, true);

		unsigned int num = get_tree(tree_)->getNOfFoundNeighbours();
		neighbors.resize(num);
		for (unsigned int i=0; i<num; ++i) {
			neighbors[i] = get_tree(tree_)->getNeighbourPositionIndex(i);
		}

		return num;
}
