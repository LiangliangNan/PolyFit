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
class PolyFitInfo;

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

	Map* generate(PolyFitInfo* polyfit_info);

private:
	// construct mesh for the bbox of the point set
	Map* construct_bbox_mesh(std::vector<Plane3d*>& supporting_planes);

	Map* compute_proxy_mesh(Map* bbox_mesh);

	// pairwise cut
	void pairwise_cut(Map* mesh);

private:
	void collect_valid_planes(std::vector<Plane3d*>& supporting_planes);

	void merge(VertexGroup* g1, VertexGroup* g2, float max_dist);

	// test if face 'f' insects plane 'plane'
	bool do_intersect(MapTypes::Facet* f, MapTypes::Facet* plane);


	struct EdgePos {
		EdgePos(MapTypes::Halfedge* e, const vec3& p) : edge(e), pos(p) {}
		MapTypes::Halfedge* edge;
		vec3				pos;
	};

	// compute the intersecting points of a face 'f' and a 'plane'. The intersecting points are returned 
	// by 'existing_vts' (if the plane intersects the face at its vertices) and 'new_vts' (if the plane 
	// intersects the face at its edges).
	void compute_intersections(
		MapTypes::Facet* f,
		MapTypes::Facet* plane,
		std::vector<MapTypes::Vertex*>& existing_vts,
		std::vector<EdgePos>& new_vts
		);

	std::vector<MapTypes::Facet*> cut(MapTypes::Facet* f, MapTypes::Facet* cutter, Map* mesh);

	// split an existing edge, meanwhile, assign the new edges the original source faces (the old edge 
	// lies in the intersection of the two faces)
	MapTypes::Vertex* split_edge(const EdgePos& ep, MapEditor* editor, MapTypes::Facet* cutter);

	// collect all faces in 'mesh' that intersect 'face'
	std::set<MapTypes::Facet*> collect_intersecting_faces(MapTypes::Facet* face, Map* mesh);

	void triplet_intersection(const std::vector<Plane3d*>& supporting_planes);

	// query the intersecting point for existing data base, i.e., triplet_intersection_
	bool query_intersection(Plane3d* plane1, Plane3d* plane2, Plane3d* plane3, vec3& p);

	// compute the intersection of a plane triplet
	// returns true if the intersection exists (p returns the point)
	bool intersection_plane_triplet(const Plane3d* plane1, const Plane3d* plane2, const Plane3d* plane3, vec3& p);

private:
	PointSet* pset_;

	MapFacetAttribute<VertexGroup*> facet_attrib_supporting_vertex_group_;

	MapFacetAttribute<Plane3d*>		facet_attrib_supporting_plane_;

	std::vector<VertexGroup::Ptr>		plane_segments_;
	std::map<VertexGroup*, Plane3d*>	vertex_group_plane_;
	
	// How to use: triplet_intersection_[plane_min][plane_mid][plane_max]
	std::map<Plane3d*, std::map<Plane3d*, std::map<Plane3d*, vec3> > >  triplet_intersection_;

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
