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

#include "hypothesis_generator.h"
#include "cgal_types.h"
#include "alpha_shape.h"
#include "polyfit_info.h"
#include "method_global.h"
#include "../basic/progress.h"
#include "../basic/logger.h"
#include "../basic/assertions.h"
#include "../basic/stop_watch.h"
#include "../model/vertex_group.h"
#include "../model/point_set.h"
#include "../model/iterators.h"
#include "../model/map.h"
#include "../model/map_builder.h"
#include "../model/map_editor.h"
#include "../model/map_circulators.h"
#include "../model/map_geometry.h"

#include <CGAL/convex_hull_2.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <algorithm>



#define	REMOVE_DEGENERATE_FACES


#ifdef REMOVE_DEGENERATE_FACES

static void remove_degenerated_facets(Map* mesh) {
	// You can't collect all the edges and then collapse them one by one, 
	// because collapsing one edge affects other neighboring edges.
	// std::vector<Map::Halfedge*> to_collapse;

	MapEditor editor(mesh);
	bool degenerated_facet_found = false;
	int count = 0;
	do {
		degenerated_facet_found = false;
		FOR_EACH_EDGE(Map, mesh, it) {
			double len = Geom::edge_length(it);
			if (len <= Method::coincident_threshold) {
				degenerated_facet_found = true;
				if (editor.collapse_edge(it)) {
					++count;
					break;  // Liangliang: only one edge can be collapsed (one may affect others)
				}
			}
		}
	} while (degenerated_facet_found);

	if (count > 0) {
		Logger::out("-") << count << " degenerate edges collapsed" << std::endl;
	}

	FOR_EACH_EDGE(Map, mesh, it) {
		double len = Geom::edge_length(it);
		if (len <= Method::coincident_threshold) {
			Logger::out("-") << count << "very short edge detected. length: " << len << std::endl;
		}
	}
}

#endif


HypothesisGenerator::HypothesisGenerator(PointSet* pset)
	: pset_(pset)
{
}


HypothesisGenerator::~HypothesisGenerator()
{
}


static std::list<unsigned int> points_on_plane(VertexGroup* g, const Plane3d& plane, float dist_threshold) {
	std::list<unsigned int> result;

	const std::vector<vec3>& points = g->point_set()->points();
	for (std::size_t i = 0; i < g->size(); ++i) {
		int idx = g->at(i);
		const vec3& p = points[idx];

		float sdist = plane.squared_ditance(p);
		float dist = std::sqrt(sdist);
		if (dist < dist_threshold)
			result.push_back(idx);
	}
	return result;
}

void HypothesisGenerator::merge(VertexGroup* g1, VertexGroup* g2, float max_dist) {
	std::vector<VertexGroup::Ptr>& groups = pset_->groups();
	const std::vector<vec3>& points = pset_->points();

	std::vector<unsigned int> points_indices;
	points_indices.insert(points_indices.end(), g1->begin(), g1->end());
	points_indices.insert(points_indices.end(), g2->begin(), g2->end());

	VertexGroup::Ptr g = new VertexGroup;
	g->insert(g->end(), points_indices.begin(), points_indices.end());
	g->set_point_set(pset_);
	g->set_color(fused_color(g1->color(), static_cast<float>(g1->size()), g2->color(), static_cast<float>(g2->size())));
	pset_->fit_plane(g);
	groups.push_back(g);

	std::vector<VertexGroup::Ptr>::iterator pos = std::find(groups.begin(), groups.end(), g1);
	if (pos != groups.end()) {
		groups.erase(pos);
	}
	else {
		std::cerr << "fatal error: vertex group doesn't exist" << std::endl;
	}

	pos = std::find(groups.begin(), groups.end(), g2);
	if (pos != groups.end()) {
		groups.erase(pos);
	}
	else {
		std::cerr << "fatal error: vertex group doesn't exist" << std::endl;
	}
}


void HypothesisGenerator::refine_planes() {
	std::vector<VertexGroup::Ptr>& groups = pset_->groups();
	const std::vector<vec3>& points = pset_->points();

	std::size_t num = groups.size();

	float avg_max_dist = 0;
	for (std::size_t i = 0; i < groups.size(); ++i) {
		VertexGroup* g = groups[i];
		pset_->fit_plane(g);

		float g_max_dist = -FLT_MAX;
		const Plane3d& plane = g->plane();
		for (std::size_t j = 0; j < g->size(); ++j) {
			int idx = g->at(j);
			const vec3& p = points[idx];
			float sdist = plane.squared_ditance(p);
			g_max_dist = std::max(g_max_dist, std::sqrt(sdist));
		}

		avg_max_dist += g_max_dist;
	}
	avg_max_dist /= groups.size();
	avg_max_dist /= 2.0f;

	float theta = 10.0f;				// in degree
	theta = static_cast<float>(M_PI * theta / 180.0f);	// in radian
	bool merged = false;
	do
	{
		merged = false;
		std::sort(groups.begin(), groups.end(), VertexGroupCmpIncreasing());

		for (std::size_t i = 0; i < groups.size(); ++i) {
			VertexGroup* g1 = groups[i];
			const Plane3d& plane1 = g1->plane();
			const vec3& n1 = plane1.normal();
			float num_threshold = g1->size() / 5.0f;
			for (std::size_t j = i + 1; j < groups.size(); ++j) {
				VertexGroup* g2 = groups[j];
				const Plane3d& plane2 = g2->plane();
				const vec3& n2 = plane2.normal();
				if (std::abs(dot(n1, n2)) > std::cos(theta)) {
					const std::list<unsigned int>& set1on2 = points_on_plane(g1, plane2, avg_max_dist);
					const std::list<unsigned int>& set2on1 = points_on_plane(g2, plane1, avg_max_dist);
					if (set1on2.size() > num_threshold || set2on1.size() > num_threshold) {
						merge(g1, g2, avg_max_dist);
						merged = true;
						break;
					}
				}
			}
			if (merged)
				break;
		}
	} while (merged);

	std::sort(groups.begin(), groups.end(), VertexGroupCmpDecreasing());

	if (num - groups.size() > 0) {
		Logger::out("-") << num - groups.size() << " planar segments merged" << std::endl;
	}
}


void HypothesisGenerator::collect_valid_planes(std::vector<Plane3d*>& supporting_planes) {
	supporting_planes.clear();
	plane_segments_.clear();
	vertex_group_plane_.clear();

	std::vector<VertexGroup::Ptr>& groups = pset_->groups();
	for (std::size_t i = 0; i < groups.size(); ++i) {
		VertexGroup* g = groups[i];
		pset_->fit_plane(g);

		plane_segments_.push_back(g);
		Plane3d* plane = new Plane3d(g->plane());
		supporting_planes.push_back(plane);
		vertex_group_plane_[g] = plane;
	}
}



static void check_source_planes(Map* mesh) {
	MapFacetAttribute<Plane3d*>					face_supporting_plane(mesh, "FacetSupportingPlane");
	MapHalfedgeAttribute< std::set<Plane3d*> >	edge_source_planes(mesh, "EdgeSourcePlanes");
	MapVertexAttribute< std::set<Plane3d*> >	vertex_source_planes(mesh, "VertexSourcePlanes");

	FOR_EACH_FACET(Map, mesh, it) {
		if (face_supporting_plane[it] == nil)
			std::cerr << "fatal error: face_supporting_plane[it] == nil" << std::endl;
	}

	FOR_EACH_HALFEDGE(Map, mesh, it) {
		const std::set<Plane3d*>& tmp = edge_source_planes[it];
		if (tmp.size() != 2)
			std::cerr << "fatal error: edge_source_planes[it].size() != 2. Size = " << tmp.size() << std::endl;
	}

	FOR_EACH_VERTEX(Map, mesh, it) {
		const std::set<Plane3d*>& tmp = vertex_source_planes[it];
		if (tmp.size() != 3)
			std::cerr << "vertex_source_planes[it].size() != 3. Size = " << tmp.size() << std::endl;
	}
}


Map* HypothesisGenerator::construct_bbox_mesh(std::vector<Plane3d*>& supporting_planes) {
	Box3d box = pset_->bbox();
	float delta = box.radius() * 0.01f;

	Map* mesh = new Map;
	MapBuilder builder(mesh);

	MapFacetAttribute<Plane3d*> face_supporting_plane(mesh, "FacetSupportingPlane");
	MapHalfedgeAttribute< std::set<Plane3d*> >	edge_source_planes(mesh, "EdgeSourcePlanes");
	MapVertexAttribute< std::set<Plane3d*> >	vertex_source_planes(mesh, "VertexSourcePlanes");

	float xmin = box.x_min() - delta, xmax = box.x_max() + delta;
	float ymin = box.y_min() - delta, ymax = box.y_max() + delta;
	float zmin = box.z_min() - delta, zmax = box.z_max() + delta;

	builder.begin_surface();

	builder.add_vertex(vec3(xmin, ymin, zmin));  // 0
	builder.add_vertex(vec3(xmax, ymin, zmin));  // 1
	builder.add_vertex(vec3(xmax, ymin, zmax));  // 2
	builder.add_vertex(vec3(xmin, ymin, zmax));  // 3
	builder.add_vertex(vec3(xmax, ymax, zmax));  // 4
	builder.add_vertex(vec3(xmax, ymax, zmin));  // 5
	builder.add_vertex(vec3(xmin, ymax, zmin));  // 6
	builder.add_vertex(vec3(xmin, ymax, zmax));  // 7


	builder.begin_facet();
	builder.add_vertex_to_facet(0);
	builder.add_vertex_to_facet(1);
	builder.add_vertex_to_facet(2);
	builder.add_vertex_to_facet(3);
	builder.end_facet();
	MapTypes::Facet* f = builder.current_facet();
	Plane3d* plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(1);
	builder.add_vertex_to_facet(5);
	builder.add_vertex_to_facet(4);
	builder.add_vertex_to_facet(2);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(1);
	builder.add_vertex_to_facet(0);
	builder.add_vertex_to_facet(6);
	builder.add_vertex_to_facet(5);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(4);
	builder.add_vertex_to_facet(5);
	builder.add_vertex_to_facet(6);
	builder.add_vertex_to_facet(7);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(0);
	builder.add_vertex_to_facet(3);
	builder.add_vertex_to_facet(7);
	builder.add_vertex_to_facet(6);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(2);
	builder.add_vertex_to_facet(4);
	builder.add_vertex_to_facet(7);
	builder.add_vertex_to_facet(3);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.end_surface();

	// assign the original planes for each edge
	FOR_EACH_HALFEDGE(Map, mesh, it) {
		Plane3d* plane1 = face_supporting_plane[it->facet()];
		Plane3d* plane2 = face_supporting_plane[it->opposite()->facet()];
		edge_source_planes[it].insert(plane1);
		edge_source_planes[it].insert(plane2);
	}

	// assign the original planes for each vertex
	FOR_EACH_VERTEX(Map, mesh, it) {
		VertexInHalfedgeCirculator cit(it);
		for (; !cit->end(); ++cit) {
			Plane3d* plane = face_supporting_plane[cit->facet()];
			vertex_source_planes[it].insert(plane);
		}
		if (cit->size() != 3)
			std::cout << "fatal_error. A bbox mesh corner does not relate to 3 planes" << std::endl;
	}

	return mesh;
}


Map* HypothesisGenerator::compute_proxy_mesh(Map* bbox_mesh) {
	MapFacetAttribute<Plane3d*>					bbox_mesh_face_supporting_plane(bbox_mesh, "FacetSupportingPlane");
	MapHalfedgeAttribute< std::set<Plane3d*> >	bbox_mesh_edge_source_planes(bbox_mesh, "EdgeSourcePlanes");
	MapVertexAttribute< std::set<Plane3d*> >	bbox_mesh_vertex_source_planes(bbox_mesh, "VertexSourcePlanes");

	Map* mesh = new Map;
	MapBuilder builder(mesh);

	MapFacetAttribute<Color> color(mesh, "color");
	MapFacetAttribute<VertexGroup*> facet_supporting_vertex_group(mesh, Method::facet_attrib_supporting_vertex_group);
	MapFacetAttribute<Plane3d*>		face_supporting_plane(mesh, "FacetSupportingPlane");
	MapHalfedgeAttribute< std::set<Plane3d*> >	edge_source_planes(mesh, "EdgeSourcePlanes");
	MapVertexAttribute< std::set<Plane3d*> >	vertex_source_planes(mesh, "VertexSourcePlanes");

	builder.begin_surface();
	int idx = 0;
	for (std::size_t i = 0; i < plane_segments_.size(); ++i) {
		VertexGroup* g = plane_segments_[i];
		Plane3d* plane = vertex_group_plane_[g];

		std::vector<vec3> points;
		std::vector< std::set<Plane3d*> > point_source_planes;
		FOR_EACH_EDGE_CONST(Map, bbox_mesh, it) {
			const vec3& s = it->prev()->vertex()->point();
			const vec3& t = it->vertex()->point();
			Sign ss = plane->orient(s);
			Sign st = plane->orient(t);
			if ((ss == POSITIVE && st == NEGATIVE) || (ss == NEGATIVE && st == POSITIVE)) {
				vec3 p;
				if (plane->intersection(Line3d::from_two_points(s, t), p)) {
					points.push_back(p);

					std::set<Plane3d*> planes = bbox_mesh_edge_source_planes[it];
					planes.insert(plane);
					point_source_planes.push_back(planes);
				}
				else
					Logger::err("-") << "fatal error. Should have intersection" << std::endl;
			}
			else {
				if (ss == ZERO) {
					points.push_back(s);

					std::set<Plane3d*> planes = bbox_mesh_vertex_source_planes[it->prev()->vertex()];
					point_source_planes.push_back(planes);
				}
				else if (st == ZERO) {
					points.push_back(t);

					std::set<Plane3d*> planes = bbox_mesh_vertex_source_planes[it->vertex()];
					point_source_planes.push_back(planes);
				}
				else {
					// no intersection with the plane
				}
			}
		}

		if (points.size() >= 3) {
			std::list<Point3> pts;
			for (std::size_t i = 0; i < points.size(); ++i) {
				const vec3& p = points[i];
				vec2 q = plane->to_2d(p);
				pts.push_back(Point3(q.x, q.y, double(i))); // trick: put the point index as the 'z' component
			}

			typedef CGAL::Projection_traits_xy_3<K>  Projection;

			std::list<Point3> hull;
			CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull), Projection());

			std::vector<vec3> ch;
			std::vector< std::set<Plane3d*> > ch_source_planes;

			for (std::list<Point3>::iterator it = hull.begin(); it != hull.end(); ++it) {
				int idx = int(it->z());
				ch.push_back(points[idx]);
				ch_source_planes.push_back(point_source_planes[idx]);
			}

			if (ch.size() >= 3) {
				for (std::size_t j = 0; j < ch.size(); ++j) {
					builder.add_vertex(ch[j]);
					MapTypes::Vertex* v = builder.current_vertex();
					vertex_source_planes[v] = ch_source_planes[j];
				}
				builder.begin_facet();
				for (int k = idx; k < idx + ch.size(); ++k)
					builder.add_vertex_to_facet(k);
				builder.end_facet();

				Map::Facet* f = builder.current_facet();
				color[f] = g->color();
				facet_supporting_vertex_group[f] = g;
				face_supporting_plane[f] = plane;

				// assign each edge the source planes
				{
					FacetHalfedgeCirculator cir(f);
					for (; !cir->end(); ++cir) {
						MapTypes::Halfedge* h = cir->halfedge();
						edge_source_planes[h].insert(plane);

						FOR_EACH_FACET(Map, bbox_mesh, fc) {
							Plane3d* bbox_plane = bbox_mesh_face_supporting_plane[fc];
							if ((bbox_plane->squared_ditance(h->vertex()->point()) < 1e-6) && (bbox_plane->squared_ditance(h->prev()->vertex()->point()) < 1e-6)) {
								edge_source_planes[h].insert(bbox_plane);
								break;
							}
						}
						std::set<Plane3d*> tmp = edge_source_planes[h];
						if (tmp.size() != 2) {
							std::cout << "fatal error: edge_source_planes[h].size() != 2. Size = " << tmp.size() << std::endl;
						}
					}
				}

				idx += int(ch.size());
			}
			else {
				Logger::err("-") << "fatal error. Check if this is a degenerate case" << std::endl;
				std::cout << "output CH points: " << ch.size() << std::endl;
				for (std::size_t i = 0; i < ch.size(); ++i) {
					const vec3& p = ch[i];
					std::cout << "\t" << p << std::endl;
				}
			}
		}
		else {
			Logger::err("-") << "fatal error. Check if this is a degenerate case" << std::endl;
		}
	}

	builder.end_surface();

	FOR_EACH_HALFEDGE(Map, mesh, it) {
		const std::set<Plane3d*>& tmp = edge_source_planes[it];
		if (tmp.size() == 2 && edge_source_planes[it->opposite()].size() != 2)
			edge_source_planes[it->opposite()] = tmp;
	}

	check_source_planes(bbox_mesh);
	check_source_planes(mesh);

	return mesh;
}


static bool halfedge_exists_between_vertices(Map::Vertex* v1, Map::Vertex* v2) {
	Map::Halfedge* cir = v1->halfedge();
	do {
		if (cir->opposite()->vertex() == v2) {
			return true;
		}
		cir = cir->next_around_vertex();
	} while (cir != v1->halfedge());
	return false;
}


// test if face f and plane intersect
bool HypothesisGenerator::do_intersect(MapTypes::Facet* f, MapTypes::Facet* plane)
{
	std::vector<MapTypes::Vertex*> existing_vts;
	std::vector<EdgePos> new_vts;
	compute_intersections(f, plane, existing_vts, new_vts);

	if (existing_vts.size() == 2) {
		if (!halfedge_exists_between_vertices(existing_vts[0], existing_vts[1]))
			return true;
	}
	else if (existing_vts.size() + new_vts.size() == 2)
		return true;

	return false;
}


std::set<Map::Facet*> HypothesisGenerator::collect_intersecting_faces(Map::Facet* face, Map* mesh) {
	std::set<Map::Facet*> intersecting_faces;
	FOR_EACH_FACET(Map, mesh, it) {
		Map::Facet* f = it;
		if (f != face && facet_attrib_supporting_vertex_group_[f] != facet_attrib_supporting_vertex_group_[face]) {
			if (do_intersect(f, face)) {
				intersecting_faces.insert(f);
			}
		}
	}
	return intersecting_faces;
}


void HypothesisGenerator::compute_intersections(
	MapTypes::Facet* f,
	MapTypes::Facet* cutter,
	std::vector<MapTypes::Vertex*>& existing_vts,
	std::vector<EdgePos>& new_vts)
{
	existing_vts.clear();
	new_vts.clear();

	double on_plane_threshold = Method::coincident_threshold * Method::coincident_threshold;
	Plane3d* plane = facet_attrib_supporting_plane_[cutter];

	Map::Halfedge* h = f->halfedge();
	do {
		const vec3& s = h->opposite()->vertex()->point();
		const vec3& t = h->vertex()->point();
		if (plane->squared_ditance(t) <= on_plane_threshold) {		// plane cuts at vertex 't'
			existing_vts.push_back(h->vertex());
		}
		else if (plane->squared_ditance(s) > on_plane_threshold) {	// cut at the edge
			const std::set<Plane3d*>& source_planes = edge_source_planes_[h];
			if (source_planes.size() == 2) {  // if the edge was computed from two faces, I use the source faces for computing the intersecting point
				if (plane->intersection(s, t)) {
					Plane3d* plane1 = *(source_planes.begin());
					Plane3d* plane2 = *(source_planes.rbegin());
					Plane3d* plane3 = plane;

					if (plane3 != plane1 && plane3 != plane2) {
						vec3 p;
						if (query_intersection(plane1, plane2, plane3, p))
							new_vts.push_back(EdgePos(h, p));
						else {
							if (intersection_plane_triplet(plane1, plane2, plane3, p)) {
								triplet_intersection_[plane1][plane2][plane3] = p; // store the intersection in our data base
								new_vts.push_back(EdgePos(h, p));
							}
							else
								Logger::warn("-") << "fatal error. should have intersection. " << new_vts.size() << std::endl;
						}
					}
				}
			}
			else {
				vec3 p;
				if (plane->intersection(s, t, p)) {
					new_vts.push_back(EdgePos(h, p));
				}
			}
		}
		h = h->next();
	} while (h != f->halfedge());
}


// split an existing edge, meanwhile, assign the new edges the original source faces (the old edge lies in the intersection of the two faces)
MapTypes::Vertex* HypothesisGenerator::split_edge(const EdgePos& ep, MapEditor* editor, MapTypes::Facet* cutter) {
	const std::set<Plane3d*>& sfs = edge_source_planes_[ep.edge];

	MapTypes::Vertex* v = editor->split_edge(ep.edge);
	v->set_point(ep.pos);

	if (sfs.size() == 2) {
		MapTypes::Halfedge* h = v->halfedge();
		edge_source_planes_[h] = sfs;
		edge_source_planes_[h->opposite()] = sfs;

		h = h->next();
		edge_source_planes_[h] = sfs;
		edge_source_planes_[h->opposite()] = sfs;
	}

	vertex_source_planes_[v] = sfs;
	vertex_source_planes_[v].insert(facet_attrib_supporting_plane_[cutter]);

	return v;
}


std::vector<Map::Facet*> HypothesisGenerator::cut(MapTypes::Facet* f, MapTypes::Facet* cutter, Map* mesh) {
	std::vector<Map::Facet*> new_faces;

	std::vector<Map::Vertex*> existing_vts;
	std::vector<EdgePos> new_vts;
	compute_intersections(f, cutter, existing_vts, new_vts);

	{ // this seem redundant, because in previous step it has been confirmed that cut is needed
		if (existing_vts.size() + new_vts.size() != 2)
			return new_faces;
		else if (existing_vts.size() == 2) {
			if (halfedge_exists_between_vertices(existing_vts[0], existing_vts[1]))
				return new_faces;
		}
	}

	//////////////////////////////////////////////////////////////////////////

	VertexGroup* g = facet_attrib_supporting_vertex_group_[f];
	MapEditor editor(mesh);

	if (existing_vts.size() + new_vts.size() == 2) {
		Map::Vertex* v0 = nil;
		Map::Vertex* v1 = nil;

		if (existing_vts.size() == 2) {
			v0 = existing_vts[0];
			v1 = existing_vts[1];

			if (halfedge_exists_between_vertices(v0, v1)) {
				std::cout << "---------fatal error: halfedge_exists_between_vertices(v0, v1)" << std::endl;
			}
		}
		else if (existing_vts.size() == 1) {
			v0 = existing_vts[0];
			v1 = split_edge(new_vts[0], &editor, cutter);
		}
		else if (new_vts.size() == 2) {
			v0 = split_edge(new_vts[0], &editor, cutter);
			v1 = split_edge(new_vts[1], &editor, cutter);
		}

		Map::Halfedge* h0 = v0->halfedge();
		if (h0->facet() != f) {
			do {
				h0 = h0->next()->opposite();
				if (h0->facet() == f)
					break;
			} while (h0 != v0->halfedge());
		}
		assert(h0->facet() == f);

		Map::Halfedge* h1 = v1->halfedge();
		if (h1->facet() != f) {
			do {
				h1 = h1->next()->opposite();
				if (h1->facet() == f)
					break;
			} while (h1 != v1->halfedge());
		}
		assert(h1->facet() == f);

		if (editor.can_split_facet(h0, h1)) {
			Map::Halfedge* h = editor.split_facet(h0, h1);
			if (h) {
				edge_source_planes_[h].insert(facet_attrib_supporting_plane_[f]);
				edge_source_planes_[h].insert(facet_attrib_supporting_plane_[cutter]);
				edge_source_planes_[h->opposite()].insert(facet_attrib_supporting_plane_[f]);
				edge_source_planes_[h->opposite()].insert(facet_attrib_supporting_plane_[cutter]);

				Map::Facet* f1 = h->facet();
				facet_attrib_supporting_vertex_group_[f1] = g;
				new_faces.push_back(f1);
				Map::Facet* f2 = h->opposite()->facet();
				facet_attrib_supporting_vertex_group_[f2] = g;
				new_faces.push_back(f2);
				return new_faces;
			}
			else
				Logger::warn("-") << "fatal error. should have intersection. " << new_vts.size() << std::endl;
		}
		else
			Logger::warn("-") << "fatal error. should have intersection. " << new_vts.size() << std::endl;
	}

	return new_faces;
}


void HypothesisGenerator::pairwise_cut(Map* mesh)
{
	std::vector<MapTypes::Facet*> all_faces;
	FOR_EACH_FACET(Map, mesh, it) {
		MapTypes::Facet* f = it;
		all_faces.push_back(f);
	}

	for (std::size_t i = 0; i < all_faces.size(); ++i) {
		MapTypes::Facet* f = all_faces[i];

		std::set<MapTypes::Facet*> intersecting_faces = collect_intersecting_faces(f, mesh);
		if (intersecting_faces.empty())
			continue;
		std::vector<MapTypes::Facet*> cutting_faces;
		cutting_faces.insert(cutting_faces.end(), intersecting_faces.begin(), intersecting_faces.end());

		// 1. f will be cut by all the intersecting_faces
		//    note: after each cut, the original face doesn't exist any more and it is replaced by multiple pieces.
		//          then each piece will be cut by another face.
		std::vector<MapTypes::Facet*> pieces;
		pieces.push_back(f);
		while (!intersecting_faces.empty()) {
			std::vector<MapTypes::Facet*> new_faces;
			MapTypes::Facet* cutter = *(intersecting_faces.begin());
			for (std::size_t j = 0; j < pieces.size(); ++j) {
				std::vector<MapTypes::Facet*> tmp = cut(pieces[j], cutter, mesh);
				new_faces.insert(new_faces.end(), tmp.begin(), tmp.end());
			}
			pieces = new_faces;
			intersecting_faces.erase(cutter);
		}

		// 2. all the cutting_faces will be cut by f.
		for (std::size_t j = 0; j < cutting_faces.size(); ++j) {
			Map::Facet* fa = cutting_faces[j];
			cut(fa, f, mesh);
		}
	}
}


// compute the intersection of a plane triplet
// returns true if the intersection exists (p returns the point)
bool HypothesisGenerator::intersection_plane_triplet(const Plane3d* plane1, const Plane3d* plane2, const Plane3d* plane3, vec3& p) {
	if (plane1 == nil || plane2 == nil || plane3 == nil) {
		Logger::err("-") << "null planes" << std::endl;
		return false;
	}
	if (plane1 == plane2 || plane2 == plane3) {
		Logger::err("-") << "identical planes" << std::endl;
		return false;
	}

	CGAL::Object obj = CGAL::intersection(to_cgal_plane(*plane1), to_cgal_plane(*plane2), to_cgal_plane(*plane3));

	// pt is the intersection point of the 3 planes 
	if (const Point3* pt = CGAL::object_cast<Point3>(&obj)) {
		p = to_my_point(*pt);
		return true;
	}
	else if (const Plane3* plane = CGAL::object_cast<Plane3>(&obj)) {
		Logger::warn("-") << "3 faces lie on the same supporting plane" << std::endl;
		return false;
	}
	else if (const Line3* line = CGAL::object_cast<Line3>(&obj)) {
		Logger::warn("-") << "3 faces intersect at the same line" << std::endl;
		return false;
	}

	return false;
}


void HypothesisGenerator::triplet_intersection(const std::vector<Plane3d*>& supporting_planes) {
	triplet_intersection_.clear();

	std::vector<Plane3d*> all_planes = supporting_planes;
	std::sort(all_planes.begin(), all_planes.end());

	for (std::size_t i = 0; i < all_planes.size(); ++i) {
		Plane3d* plane1 = all_planes[i];
		for (std::size_t j = i + 1; j < all_planes.size(); ++j) {
			Plane3d* plane2 = all_planes[j];
			for (std::size_t k = j + 1; k < all_planes.size(); ++k) {
				Plane3d* plane3 = all_planes[k];

				assert(plane1 < plane2);
				assert(plane2 < plane3);

				vec3 p;
				if (intersection_plane_triplet(plane1, plane2, plane3, p))
					triplet_intersection_[plane1][plane2][plane3] = p; // store the intersection in our data base
			}
		}
	}
}


bool HypothesisGenerator::query_intersection(Plane3d* plane1, Plane3d* plane2, Plane3d* plane3, vec3& p) {
	Plane3d* min_plane = ogf_min(plane1, plane2, plane3);
	Plane3d* max_plane = ogf_max(plane1, plane2, plane3);

	Plane3d* mid_plane = 0;
	if (plane1 != min_plane && plane1 != max_plane)
		mid_plane = plane1;
	else if (plane2 != min_plane && plane2 != max_plane)
		mid_plane = plane2;
	else if (plane3 != min_plane && plane3 != max_plane)
		mid_plane = plane3;

	if (triplet_intersection_.find(min_plane) == triplet_intersection_.end())
		return false;

	std::map<Plane3d*, std::map<Plane3d*, vec3> >& tmp2 = triplet_intersection_[min_plane];
	if (tmp2.find(mid_plane) == tmp2.end())
		return false;

	std::map<Plane3d*, vec3>& tmp3 = tmp2[mid_plane];
	if (tmp3.find(max_plane) == tmp3.end())
		return false;

	p = tmp3[max_plane];
	return true;
}



Map* HypothesisGenerator::generate(PolyFitInfo* polyfit_info) {
	if (!pset_)
		return nil;

	if (pset_->groups().empty()) {
		Logger::warn("-") << "planar segments do not exist" << std::endl;
		return nil;
	}

	collect_valid_planes(polyfit_info->planes);

	Map* bbox_mesh = construct_bbox_mesh(polyfit_info->planes);
	Map* mesh = compute_proxy_mesh(bbox_mesh);
	if (!mesh)
		return nil;

	check_source_planes(mesh);

	facet_attrib_supporting_vertex_group_.bind(mesh, Method::facet_attrib_supporting_vertex_group);
	facet_attrib_supporting_plane_.bind(mesh, "FacetSupportingPlane");
	edge_source_planes_.bind(mesh, "EdgeSourcePlanes");
	vertex_source_planes_.bind(mesh, "VertexSourcePlanes");

	triplet_intersection(polyfit_info->planes);
	pairwise_cut(mesh);
	check_source_planes(mesh);

#ifdef REMOVE_DEGENERATE_FACES
	remove_degenerated_facets(mesh);
	check_source_planes(mesh);
	FOR_EACH_FACET(Map, mesh, it) {
		if (Geom::facet_area(it) < 1e-10) 
			std::cerr << "degenerate face detected" << std::endl;
	}
#endif

	facet_attrib_supporting_vertex_group_.unbind();
	facet_attrib_supporting_plane_.unbind();
	edge_source_planes_.unbind();
	vertex_source_planes_.unbind();

	return mesh;
}