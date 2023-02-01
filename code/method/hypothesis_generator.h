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

#ifndef _HYPOTHESIS_GENERATOR_
#define _HYPOTHESIS_GENERATOR_

#include "method_common.h"
#include "../math/polygon2d.h"
#include "../model/vertex_group.h"
#include "../model/map_attributes.h"

#include <string>
#include <vector>


class Map;
class PointSet;
class VertexGroup;
class MapEditor;
class ProgressLogger;

namespace MapTypes {
	class Vertex;
	class Facet;
	class Halfedge;
}

class METHOD_API HypothesisGenerator
{
public:
	HypothesisGenerator(PointSet* pset);
	~HypothesisGenerator();

	void refine_planes();

	Map* generate();

	void compute_confidences(Map* mesh, bool use_conficence = false);

	// Intersection: a set of 'faces' intersecting at a common edge
	struct SuperEdge : public std::vector<MapTypes::Halfedge*> {
		const vec3* s;
		const vec3* t;
	};
	typedef typename std::vector<SuperEdge>		Adjacency;

	// the adjacency information will be used to formulate the hard constraints.
	Adjacency extract_adjacency(Map* mesh);

	bool ready_for_optimization(Map* mesh) const;

private:
	// construct mesh for the bbox of the point set
	Map* construct_bbox_mesh();

	Map* compute_proxy_mesh(Map* bbox_mesh);

	// pairwise cut
	void pairwise_cut(Map* mesh);

private:
	void collect_valid_planes();

	void merge(VertexGroup* g1, VertexGroup* g2);

	// test if face 'f' insects plane 'plane'
	bool do_intersect(MapTypes::Facet* f, Plane3d* plane);

    struct Intersection {
        enum Type { EXISTING_VERTEX, NEW_VERTEX };

        Intersection(Type t) : type(t), vtx(0), edge(0) {}
        Type type;

        // for EXISTING_VERTEX
        Map::Vertex* vtx;

        // for NEW_VERTEX
        MapTypes::Halfedge* edge;
        vec3				pos;
    };

	// compute the intersecting points of a face 'f' and a 'plane'. The intersecting points are returned 
	// by 'existing_vts' (if the plane intersects the face at its vertices) and 'new_vts' (if the plane 
	// intersects the face at its edges).
    void compute_intersections(
            MapTypes::Facet* f,
            Plane3d* plane,
            std::vector<Intersection>& intersections
    );

	std::vector<MapTypes::Facet*> cut(MapTypes::Facet* f, Plane3d* cutter, Map* mesh);

	// split an existing edge, meanwhile, assign the new edges the original source faces (the old edge 
	// lies in the intersection of the two faces)
	MapTypes::Vertex* split_edge(const Intersection& ep, MapEditor* editor, Plane3d* cutting_plane);

	// collect all faces in 'mesh' that intersect 'face'
	std::set<Plane3d *> collect_cutting_planes(MapTypes::Facet* face, Map* mesh);

	void triplet_intersection();

	// query the intersecting point for existing data base, i.e., triplet_intersection_
	vec3* query_intersection(Plane3d* plane1, Plane3d* plane2, Plane3d* plane3);

	// compute the intersection of a plane triplet
	// returns true if the intersection exists (p returns the point)
	bool intersection_plane_triplet(const Plane3d* plane1, const Plane3d* plane2, const Plane3d* plane3, vec3& p);

	// the pairwise intersection may result in tiny faces and we may have numerical problems when computing the 
	// face confidences where face area is the denominator. To handle this, we simply remove these degenerate 
	// faces by collapsing the edges.
	void remove_degenerated_facets(Map* mesh);

	// std::vector<unsigned int>& points returns the point indices projected in f.
	// returns the 'number' of points projected in f (accounts for a notion of confidence)
	float facet_points_projected_in(PointSet* pset, VertexGroup* g, MapTypes::Facet* f, float max_dist, std::vector<unsigned int>& points);

	// returns average spacing
	float compute_point_confidences(PointSet* pset, int s1 = 6, int s2 = 16, int s3 = 32, ProgressLogger* progress = nullptr);

	// clear cached intermediate results
	void clear();

private:
	PointSet* pset_;

	MapFacetAttribute<VertexGroup*> facet_attrib_supporting_vertex_group_;
	MapFacetAttribute<Plane3d*>		facet_attrib_supporting_plane_;

	std::vector<VertexGroup::Ptr>		plane_segments_;
	std::map<VertexGroup*, Plane3d*>	vertex_group_plane_;

	std::vector<Plane3d*>  supporting_planes_;		// including the bbox face planes
	float				   max_dist_;				// maximum distance to the supporting plane
	
	// How to use: triplet_intersection_[plane_min][plane_mid][plane_max]
	std::map<Plane3d*, std::map<Plane3d*, std::map<Plane3d*, vec3*> > >  triplet_intersection_;

	// precomputed intersecting points of all plane triplets
	std::vector< vec3*>		intersecting_points_;

	// to avoid numerical issues (there are always small differences when computing the intersecting 
	// point of a plane triplet), I store how a edge is computed (from two planes). Then, I just need 
	// to query the intersecting point of a plane triplet 
	// from a precomputed table. By doing so, I can avoid this numerical issues.  
	MapHalfedgeAttribute< std::set<Plane3d*> > edge_source_planes_;

	// to avoid numerical issues (due to floating point precision, there are always small differences 
	// when computing the intersection of a plane triplet), I store how a vertex is computed (from
	// three planes). Then, I just need to compare the three plane to identify if two points are the same. 
	MapVertexAttribute< std::set<Plane3d*> > vertex_source_planes_;
};

#endif
