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

#ifndef __KDTREE_KDTREE_SEARCH__
#define __KDTREE_KDTREE_SEARCH__

#include "model_common.h"
#include "../math/math_types.h"
#include "../basic/counted.h"
#include "../basic/smart_pointer.h"

#include <list>



class PointSet;

class MODEL_API KdTreeSearch : public Counted {
public:
	KdTreeSearch();
	virtual ~KdTreeSearch();

	//______________ tree construction __________________________

	virtual void begin() ;
	virtual void add_point(vec3* v) ;
	virtual void add_vertex_set(PointSet* vs) ;
	virtual void end() ;

	//________________ closest point ____________________________

	// return the index of the closest point, -1 if not found
	// NOTE: *squared* distance is returned
	virtual int find_closest_point(const vec3& p, double& squared_distance) const ;
	virtual int find_closest_point(const vec3& p) const ;

	//_________________ K-nearest neighbors ____________________

	// NOTE: *squared* distances are returned
	virtual void find_closest_K_points(
		const vec3& p, unsigned int k, 
		std::vector<unsigned int>& neighbors, std::vector<double>& squared_distances
		) const ;

	virtual void find_closest_K_points(
		const vec3& p, unsigned int k, 
		std::vector<unsigned int>& neighbors
		) const;

	//___________________ radius search __________________________

	// fixed-radius kNN	search. Search for all points in the range.
	// NOTE: *squared* radius of query ball
	virtual void find_points_in_radius(const vec3& p, double squared_radius, 
		std::vector<unsigned int>& neighbors
		) const ;

	virtual void find_points_in_radius(const vec3& p, double squared_radius, 
		std::vector<unsigned int>& neighbors, std::vector<double>& squared_distances
		) const ;

	//____________________ cylinder range search _________________

	// Search for the nearest points whose distances to line segment $v1$-$v2$ are smaller 
	// than $radius$. If $bToLine$ is true, the points found are ordered by their distances 
	// to the line segment. Otherwise, they are ordered by their distances to $v1$.
	// NOTE: it is radius (instead of *squared* radius).
	unsigned int find_points_in_cylinder(
		const vec3& p1, const vec3& p2, double radius, 
		std::vector<unsigned int>& neighbors, std::vector<double>& squared_distances,
		bool bToLine = true
		) const ;

	unsigned int find_points_in_cylinder(
		const vec3& p1, const vec3& p2, double radius, 
		std::vector<unsigned int>& neighbors,
		bool bToLine = true
		) const ;
	
	//_______________________ cone range search __________________

	// Search for the nearest points $P_i$ with an cone from $v1$ to $v2$ defined by v1 and v2. 
	// As a result, the angle between $v1$$P_i$ and $v1$$v2$ is smaller than $angle_range$.
	// Search for the nearest points P_i where the angle between $v1$-P_i and $v1$-$v2$ is 
	// smaller than $angle$.
	// NOTE: angle is in radian.
	unsigned int find_points_in_cone(
		const vec3& eye, const vec3& p1, const vec3& p2, double angle_range, 
		std::vector<unsigned int>& neighbors, std::vector<double>& squared_distances,
		bool bToLine = true
		) const ;

	unsigned int find_points_in_cone(
		const vec3& eye, const vec3& p1, const vec3& p2, double angle_range, 
		std::vector<unsigned int>& neighbors,
		bool bToLine = true
		) const ;

protected:
	std::list<vec3*>	vertices_;
	unsigned int		points_num_;
	void*				tree_;
} ;



typedef	SmartPointer<KdTreeSearch>	KdTreeSearch_var;
#endif


