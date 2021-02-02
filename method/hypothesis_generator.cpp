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

#include "hypothesis_generator.h"
#include "cgal_types.h"
#include "alpha_shape.h"
#include "method_global.h"
#include "alpha_shape_mesh.h"
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
#include "../model/kdtree_search.h"

#include <CGAL/convex_hull_2.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <algorithm>




// For debugging (identifying topological issues)
//#define DISPLAY_ADJACENCY_STATISTICS


HypothesisGenerator::HypothesisGenerator(PointSet* pset)
	: pset_(pset)
{
}


HypothesisGenerator::~HypothesisGenerator()
{
	clear();
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

void HypothesisGenerator::merge(VertexGroup* g1, VertexGroup* g2) {
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
		Logger::err("-") << "fatal error: vertex group doesn't exist" << std::endl;
	}

	pos = std::find(groups.begin(), groups.end(), g2);
	if (pos != groups.end()) {
		groups.erase(pos);
	}
	else {
		Logger::err("-") << "fatal error: vertex group doesn't exist" << std::endl;
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
						merge(g1, g2);
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


void HypothesisGenerator::collect_valid_planes() {
	supporting_planes_.clear();
	plane_segments_.clear();
	vertex_group_plane_.clear();

	std::vector<VertexGroup::Ptr>& groups = pset_->groups();
	for (std::size_t i = 0; i < groups.size(); ++i) {
		VertexGroup* g = groups[i];
		pset_->fit_plane(g);

		plane_segments_.push_back(g);
		Plane3d* plane = new Plane3d(g->plane());
		supporting_planes_.push_back(plane);
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


Map* HypothesisGenerator::construct_bbox_mesh() {
	Box3d box = pset_->bbox();
    float delta = box.radius() * 0.05f;

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
	supporting_planes_.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(1);
	builder.add_vertex_to_facet(5);
	builder.add_vertex_to_facet(4);
	builder.add_vertex_to_facet(2);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes_.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(1);
	builder.add_vertex_to_facet(0);
	builder.add_vertex_to_facet(6);
	builder.add_vertex_to_facet(5);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes_.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(4);
	builder.add_vertex_to_facet(5);
	builder.add_vertex_to_facet(6);
	builder.add_vertex_to_facet(7);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes_.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(0);
	builder.add_vertex_to_facet(3);
	builder.add_vertex_to_facet(7);
	builder.add_vertex_to_facet(6);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes_.push_back(plane);
	face_supporting_plane[f] = plane;

	builder.begin_facet();
	builder.add_vertex_to_facet(2);
	builder.add_vertex_to_facet(4);
	builder.add_vertex_to_facet(7);
	builder.add_vertex_to_facet(3);
	builder.end_facet();
	f = builder.current_facet();
	plane = new Plane3d(Geom::facet_plane(f));
	supporting_planes_.push_back(plane);
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
			Logger::err("-") << "fatal_error. A bbox mesh corner does not relate to 3 planes" << std::endl;
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
							Logger::err("-") << "fatal error: edge_source_planes[h].size() != 2. Size = " << tmp.size() << std::endl;
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
bool HypothesisGenerator::do_intersect(MapTypes::Facet* f, Plane3d* plane)
{
    std::vector<Intersection> vts;
    compute_intersections(f, plane, vts);

    return (vts.size() > 1);
}


std::set<Plane3d *> HypothesisGenerator::collect_cutting_planes(MapTypes::Facet *face, Map *mesh) {
	std::set<Plane3d*> cutting_planes;
	FOR_EACH_FACET(Map, mesh, it) {
		Map::Facet* f = it;
		if (f != face) {
		    Plane3d* plane = facet_attrib_supporting_plane_[f];
		    if (plane != facet_attrib_supporting_plane_[face]) {
			    if (do_intersect(f, facet_attrib_supporting_plane_[face]))
                    cutting_planes.insert(plane);
			}
		}
	}
	return cutting_planes;
}


void HypothesisGenerator::compute_intersections(
        MapTypes::Facet* f,
        Plane3d* plane,
        std::vector<Intersection>& vts)
{
    vts.clear();

    Map::Halfedge* h = f->halfedge();
    do {
        const vec3& s = h->opposite()->vertex()->point();
        const vec3& t = h->vertex()->point();
        if (plane->squared_ditance(t) <= Method::snap_sqr_distance_threshold) {		// plane cuts at vertex 't'
            Intersection it(Intersection::EXISTING_VERTEX);
            it.vtx = h->vertex();
            vts.push_back(it);
        }
        else if (plane->squared_ditance(s) > Method::snap_sqr_distance_threshold) {	// cut at the edge
            const std::set<Plane3d*>& source_planes = edge_source_planes_[h];
            if (source_planes.size() == 2) {  // if the edge was computed from two faces, I use the source faces for computing the intersecting point
                if (plane->intersection(s, t)) {
                    Plane3d* plane1 = *(source_planes.begin());
                    Plane3d* plane2 = *(source_planes.rbegin());
                    Plane3d* plane3 = plane;
                    if (plane3 != plane1 && plane3 != plane2) {
                        vec3* p = query_intersection(plane1, plane2, plane3);
                        if (p) {
                            Intersection it(Intersection::NEW_VERTEX);
                            it.edge = h;
                            it.pos = *p;
                            vts.push_back(it);
                        }
                        else {
                            vec3 q;
                            if (intersection_plane_triplet(plane1, plane2, plane3, q)) {
                                triplet_intersection_[plane1][plane2][plane3] = p; // store the intersection in our data base
                                Intersection it(Intersection::NEW_VERTEX);
                                it.edge = h;
                                it.pos = q;
                                vts.push_back(it);
                            }
                            else
                                Logger::warn("-") << "fatal error. should have intersection. " << std::endl;
                        }
                    }
                    else {
                        Logger::warn("-") << "fatal error. should have 3 different planes. " << std::endl;
                    }
                }
            }
            else {
                vec3 p;
                if (plane->intersection(s, t, p)) {
                    Intersection it(Intersection::NEW_VERTEX);
                    it.edge = h;
                    it.pos = p;
                    vts.push_back(it);
                }
            }
        }
        h = h->next();
    } while (h != f->halfedge());
}


// split an existing edge, meanwhile, assign the new edges the original source faces (the old edge lies in the intersection of the two faces)
MapTypes::Vertex* HypothesisGenerator::split_edge(const Intersection& ep, MapEditor* editor, Plane3d* cutter) {
	const std::set<Plane3d*>& sfs = edge_source_planes_[ep.edge];
	assert(sfs.size() == 2);

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
    else {
        Logger::warn("-") << "edge_source_planes.size != 2" << std::endl;
    }

    vertex_source_planes_[v] = sfs;
    vertex_source_planes_[v].insert(cutter);

    return v;
}


std::vector<Map::Facet*> HypothesisGenerator::cut(MapTypes::Facet* f, Plane3d* cutter, Map* mesh) {
	std::vector<Map::Facet*> new_faces;

    std::vector<Intersection> vts;
    compute_intersections(f, cutter, vts);
    if (vts.size() < 2) { // no actual intersection
        return new_faces;
    }

    MapEditor editor(mesh);

    Map::Vertex* v0 = nil;
    Map::Vertex* v1 = nil;
    if (vts.size() == 2) {
		// test if the two intersecting points are both very close to an existing vertex.
		// Since we allow snapping, we test if the two intersecting points are the same.
		if (vts[0].type == Intersection::NEW_VERTEX && vts[1].type == Intersection::NEW_VERTEX) {
            v0 = split_edge(vts[0], &editor, cutter);
            v1 = split_edge(vts[1], &editor, cutter);
        }
        else if (vts[0].type == Intersection::NEW_VERTEX && vts[1].type == Intersection::EXISTING_VERTEX) {
            v0 = split_edge(vts[0], &editor, cutter);
            v1 = vts[1].vtx;
        }
        else if (vts[0].type == Intersection::EXISTING_VERTEX && vts[1].type == Intersection::NEW_VERTEX) {
            v0 = vts[0].vtx;
            v1 = split_edge(vts[1], &editor, cutter);
        }
        else if (vts[0].type == Intersection::EXISTING_VERTEX && vts[1].type == Intersection::EXISTING_VERTEX) {
            if (halfedge_exists_between_vertices(vts[0].vtx, vts[1].vtx))
			    return new_faces;
            else {
                v0 = vts[0].vtx;
                v1 = vts[1].vtx;
            }
        }
	}
	else {
        Logger::warn("-") << "This might be an error: number of intersecting points is " << vts.size() << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////

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

    VertexGroup* g = facet_attrib_supporting_vertex_group_[f];
	if (editor.can_split_facet(h0, h1)) {
		Map::Halfedge* h = editor.split_facet(h0, h1);
		if (h) {
			edge_source_planes_[h].insert(facet_attrib_supporting_plane_[f]);
			edge_source_planes_[h].insert(cutter);
			edge_source_planes_[h->opposite()].insert(facet_attrib_supporting_plane_[f]);
			edge_source_planes_[h->opposite()].insert(cutter);

			Map::Facet* f1 = h->facet();
			facet_attrib_supporting_vertex_group_[f1] = g;
			new_faces.push_back(f1);
			Map::Facet* f2 = h->opposite()->facet();
			facet_attrib_supporting_vertex_group_[f2] = g;
			new_faces.push_back(f2);
		}
		else
			Logger::warn("-") << "fatal error. should have intersection." << std::endl;
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

	std::map< MapTypes::Facet*, std::set<Plane3d *> > face_cutters;
    for (std::size_t i = 0; i < all_faces.size(); ++i) {
        MapTypes::Facet *f = all_faces[i];
        face_cutters[f] = collect_cutting_planes(f, mesh);
    }

	ProgressLogger progress(all_faces.size());
	for (std::size_t i = 0; i < all_faces.size(); ++i) {
		MapTypes::Facet* f = all_faces[i];

		std::set<Plane3d*>& cutting_planes = face_cutters[f];
		if (cutting_planes.empty())
			continue;

		// f will be cut by all the intersecting_faces
		// note: after each cut, the original face doesn't exist any more and it is replaced by multiple pieces.
		//       then each piece will be cut by another face.
		std::vector<MapTypes::Facet*> faces_to_be_cut;
		faces_to_be_cut.push_back(f);
		while (!cutting_planes.empty()) {
			std::set<MapTypes::Facet*> new_faces;		// stores the new faces
			std::set<MapTypes::Facet*> remained_faces;	// faces that will be cut later
            Plane3d* cutter = *(cutting_planes.begin());
			for (std::size_t j = 0; j < faces_to_be_cut.size(); ++j) {
				MapTypes::Facet* current_face = faces_to_be_cut[j];
				std::vector<MapTypes::Facet*> tmp = cut(current_face, cutter, mesh);
				new_faces.insert(tmp.begin(), tmp.end());
				if (tmp.empty()) {
					remained_faces.insert(current_face);
				}
			}
			faces_to_be_cut = std::vector<MapTypes::Facet*>(new_faces.begin(), new_faces.end());
			faces_to_be_cut.insert(faces_to_be_cut.end(), remained_faces.begin(), remained_faces.end());
            cutting_planes.erase(cutter);
		}

		progress.next();
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


void HypothesisGenerator::triplet_intersection() {
	triplet_intersection_.clear();

	std::sort(supporting_planes_.begin(), supporting_planes_.end());

	for (std::size_t i = 0; i < supporting_planes_.size(); ++i) {
		Plane3d* plane1 = supporting_planes_[i];
		for (std::size_t j = i + 1; j < supporting_planes_.size(); ++j) {
			Plane3d* plane2 = supporting_planes_[j];
			for (std::size_t k = j + 1; k < supporting_planes_.size(); ++k) {
				Plane3d* plane3 = supporting_planes_[k];

				assert(plane1 < plane2);
				assert(plane2 < plane3);

				vec3 p;
				if (intersection_plane_triplet(plane1, plane2, plane3, p)) {
					vec3* new_point = new vec3(p);
					triplet_intersection_[plane1][plane2][plane3] = new_point; // store the intersection in our data base
					intersecting_points_.push_back(new_point);
				}
			}
		}
	}
}


void HypothesisGenerator::remove_degenerated_facets(Map* mesh) {
	// You can't collect all the edges and then collapse them one by one, 
	// because collapsing one edge affects other neighboring edges.
	// std::vector<Map::Halfedge*> to_collapse;

	MapEditor editor(mesh);
	bool degenerated_facet_found = false;
	int count = 0;
	do {
		degenerated_facet_found = false;
		FOR_EACH_EDGE(Map, mesh, it) {
			double sd = distance2(it->vertex()->point(), it->prev()->vertex()->point());
			if (sd < Method::snap_sqr_distance_threshold) {
				degenerated_facet_found = true;
				if (editor.collapse_edge(it)) {
					++count;
					break;  // Liangliang: only one edge can be collapsed (one may affect others)
				}
			}
		}
	} while (degenerated_facet_found);

	if (count > 0)
		Logger::out("-") << count << " degenerate edges collapsed" << std::endl;
}


vec3* HypothesisGenerator::query_intersection(Plane3d* plane1, Plane3d* plane2, Plane3d* plane3) {
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
		return 0;

	std::map<Plane3d*, std::map<Plane3d*, vec3*> >& tmp2 = triplet_intersection_[min_plane];
	if (tmp2.find(mid_plane) == tmp2.end())
		return 0;

	std::map<Plane3d*, vec3*>& tmp3 = tmp2[mid_plane];
	if (tmp3.find(max_plane) == tmp3.end())
		return 0;

	return tmp3[max_plane];
}



Map* HypothesisGenerator::generate() {
	if (!pset_)
		return nil;

	if (pset_->groups().empty()) {
		Logger::warn("-") << "planar segments do not exist" << std::endl;
		return nil;
	}

	collect_valid_planes();

	Map* bbox_mesh = construct_bbox_mesh();
	Map* mesh = compute_proxy_mesh(bbox_mesh);
	if (!mesh)
		return nil;

	check_source_planes(mesh);

	facet_attrib_supporting_vertex_group_.bind(mesh, Method::facet_attrib_supporting_vertex_group);
	facet_attrib_supporting_plane_.bind(mesh, "FacetSupportingPlane");
	edge_source_planes_.bind(mesh, "EdgeSourcePlanes");
	vertex_source_planes_.bind(mesh, "VertexSourcePlanes");

	triplet_intersection();
	pairwise_cut(mesh);
	check_source_planes(mesh);

	remove_degenerated_facets(mesh);
	check_source_planes(mesh);

	facet_attrib_supporting_vertex_group_.unbind();
	facet_attrib_supporting_plane_.unbind();
	edge_source_planes_.unbind();
	vertex_source_planes_.unbind();

	return mesh;
}


void HypothesisGenerator::clear() {
	for (std::size_t i = 0; i < supporting_planes_.size(); ++i)
		delete supporting_planes_[i];
	supporting_planes_.clear();

	for (std::size_t i = 0; i < intersecting_points_.size(); ++i)
		delete intersecting_points_[i];
	intersecting_points_.clear();
}




float HypothesisGenerator::compute_point_confidences(PointSet* pset, int s1 /* = 6 */, int s2 /* = 16 */, int s3 /* = 32 */, ProgressLogger* progress) {
	std::vector<vec3>& points = pset->points();
	std::vector<float>& planar_qualities = pset->planar_qualities();

	if (planar_qualities.size() != points.size())
		planar_qualities.resize(points.size());

	KdTreeSearch_var kdtree = new KdTreeSearch;
	kdtree->begin();
	kdtree->add_vertex_set(pset);
	kdtree->end();

	std::vector<int> neighbor_size;
	neighbor_size.push_back(s1);
	neighbor_size.push_back(s2);
	neighbor_size.push_back(s3);

	std::vector< std::vector<double> > eigen_values(3);
	for (int i = 0; i < 3; ++i)
		eigen_values[i].resize(3);

	double total = 0;
	for (std::size_t i = 0; i < points.size(); ++i) {
		const vec3& p = points[i];
		for (int j = 0; j < 3; ++j) {
			std::vector<unsigned int> neighbors;
			std::vector<double> sqr_distances;
			kdtree->find_closest_K_points(p, neighbor_size[j], neighbors, sqr_distances);

			PrincipalAxes3d pca;
			pca.begin();
			double avg = 0;
			for (unsigned int k = 0; k < neighbors.size(); ++k) {
				pca.add_point(points[neighbors[k]]);
				if (j == 0) {
					avg += std::sqrt(sqr_distances[k]);
				}
			}
			pca.end();
			total += (avg / neighbors.size());

			for (int k = 0; k < 3; ++k)
				eigen_values[j][k] = pca.eigen_value(3 - k - 1); // eigen values are sorted in descending order

			assert(eigen_values[j][0] <= eigen_values[j][1] && eigen_values[j][1] <= eigen_values[j][2]);
		}

		double conf = 0.0;
		for (int j = 0; j < 3; ++j) {
			conf += (1 - 3.0 * eigen_values[j][0] / (eigen_values[j][0] + eigen_values[j][1] + eigen_values[j][2])) * (eigen_values[j][1] / eigen_values[j][2]);
		}
		conf /= 3.0;
		planar_qualities[i] = static_cast<float>(conf);

		if (progress)
		    progress->next();
	}
	return static_cast<float>(total / points.size());
}


void HypothesisGenerator::compute_confidences(Map* mesh, bool use_conficence /* = false */) {
	facet_attrib_supporting_vertex_group_.bind(mesh, Method::facet_attrib_supporting_vertex_group);

	StopWatch w;
	Logger::out("-") << "computing point confidences..." << std::endl;
    ProgressLogger progress(pset_->num_points() + mesh->size_of_facets());
	double avg_spacing = compute_point_confidences(pset_, 6, 16, 25, &progress);
	float radius = static_cast<float>(avg_spacing)* 5.0f;
	Logger::out("-") << "done. avg spacing: " << avg_spacing << ". " << w.elapsed() << " sec." << std::endl;

	std::vector<VertexGroup::Ptr>& groups = pset_->groups();
	const std::vector<vec3>& pts = pset_->points();

	float max_dist = 0;
	for (std::size_t i = 0; i < groups.size(); ++i) {
		VertexGroup* g = groups[i];
		const Plane3d& plane = g->plane();
		for (std::size_t j = 0; j < g->size(); ++j) {
			int idx = g->at(j);
			const vec3& p = pts[idx];
			float sdist = plane.squared_ditance(p);
			max_dist = std::max(max_dist, sdist);
		}
	}
	max_dist = std::sqrt(max_dist);

	//////////////////////////////////////////////////////////////////////////

	MapFacetAttribute<double>	facet_attrib_supporting_point_num(mesh, Method::facet_attrib_supporting_point_num);
	MapFacetAttribute<double>	facet_attrib_facet_area(mesh, Method::facet_attrib_facet_area);
	MapFacetAttribute<double>	facet_attrib_covered_area(mesh, Method::facet_attrib_covered_area);

	Logger::out("-") << "computing face confidences..." << std::endl;
	w.start();
	FOR_EACH_FACET(Map, mesh, it) {
		Map::Facet* f = it;

		double face_area = Geom::facet_area(f);
		if (face_area < 1e-16) {
			Logger::err("-") << "degenerate facet with area: " << face_area << std::endl;
			FacetHalfedgeCirculator cir(f);
			for (; !cir->end(); ++cir) {
				Logger::err("-") << cir->vertex()->point() << std::endl;
			}
			Logger::err("-") << std::endl;
			continue;
		}

		VertexGroup* g = facet_attrib_supporting_vertex_group_[f];

		std::vector<unsigned int> points;
		double num = facet_points_projected_in(pset_, g, f, max_dist, points);
		if (use_conficence)
			facet_attrib_supporting_point_num[f] = num;
		else
			facet_attrib_supporting_point_num[f] = static_cast<double>(points.size());

		facet_attrib_facet_area[f] = face_area;

		Map::Ptr alpha_mesh = AlphaShapeMesh::apply(pset_, points, g->plane(), radius);
		double covered_area = 0;
		if (alpha_mesh) {
			FOR_EACH_FACET(Map, alpha_mesh, it)
				covered_area += Geom::triangle_area(it);
		}
		facet_attrib_covered_area[f] = covered_area;

		if (covered_area > face_area) {
			// this may not be an error (floating point precision limit)
			facet_attrib_covered_area[f] = face_area;
		}
		progress.next();
	}

	facet_attrib_supporting_vertex_group_.unbind();

	Logger::out("-") << "done. " << w.elapsed() << " sec." << std::endl;
}


float HypothesisGenerator::facet_points_projected_in(PointSet* pset, VertexGroup* g, MapTypes::Facet* f, float max_dist, std::vector<unsigned int>& points) {
	const Plane3d& plane = g->plane();
	const vec3& orig = plane.point();
	const vec3& base1 = plane.base1();
	const vec3& base2 = plane.base2();

	const Polygon3d& plg3d = Geom::facet_polygon(f);
	const Polygon2d& plg2d = Geom::to_2d(orig, base1, base2, plg3d);
	const std::vector<vec3>& pts = pset->points();
	const std::vector<float>& confidences = pset->planar_qualities();

	points.clear();
	float epsilon = max_dist * 0.5f;// considering noise and outliers
	float count = 0.0f;
	for (int i = 0; i < g->size(); ++i) {
		unsigned int idx = g->at(i);
		const vec3& p = pts[idx];
		const vec2& q = Geom::to_2d(orig, base1, base2, p);
		if (Geom::point_is_in_polygon(plg2d, q)) {
			points.push_back(idx);
			float dist = std::sqrt(plane.squared_ditance(p));
			if (dist < epsilon) { // in case of numerical issues (floating point precision)
				count += (1 - dist / epsilon) * confidences[idx];
			}
		}
	}

	return count;
}



HypothesisGenerator::Adjacency HypothesisGenerator::extract_adjacency(Map* mesh) {
	vertex_source_planes_.bind(mesh, "VertexSourcePlanes");

	// an edge is denoted by its two end points
	typedef typename std::map< vec3*, std::set<MapTypes::Halfedge*> >	Edge_map;
	typedef typename std::map< vec3*, Edge_map >						Face_pool;
	Face_pool face_pool;

	FOR_EACH_HALFEDGE(Map, mesh, h) {
		if (h->facet() == 0)
			continue;

		Map::Vertex* sd = h->opposite()->vertex();
		Map::Vertex* td = h->vertex();

		std::set<Plane3d*>& set_s = vertex_source_planes_[sd];
		std::set<Plane3d*>& set_t = vertex_source_planes_[td];
		CGAL_assertion(set_s.size() == 3);
		CGAL_assertion(set_t.size() == 3);

		std::vector<Plane3d*> s_planes(set_s.begin(), set_s.end());
		CGAL_assertion(s_planes[0] < s_planes[1]);
		CGAL_assertion(s_planes[1] < s_planes[2]);
		vec3* s = triplet_intersection_[s_planes[0]][s_planes[1]][s_planes[2]];

		std::vector<Plane3d*> t_planes(set_t.begin(), set_t.end());
		CGAL_assertion(t_planes[0] < t_planes[1]);
		CGAL_assertion(t_planes[1] < t_planes[2]);
		vec3* t = triplet_intersection_[t_planes[0]][t_planes[1]][t_planes[2]];

		if (s > t)
			std::swap(s, t);
		face_pool[s][t].insert(h);
	}

#ifdef DISPLAY_ADJACENCY_STATISTICS
	std::map<std::size_t, std::size_t>   num_each_sized_fans;
	for (std::size_t i = 1; i < 20; ++i)
		num_each_sized_fans[i] = 0;
#endif

	Adjacency fans;
	Face_pool::const_iterator it = face_pool.begin();
	for (; it != face_pool.end(); ++it) {
		const vec3* s = it->first;
		const Edge_map& tmp = it->second;
		Edge_map::const_iterator cur = tmp.begin();
		for (; cur != tmp.end(); ++cur) {
			const vec3* t = cur->first;
			const std::set<Map::Halfedge*>& faces = cur->second;
			SuperEdge fan;
			fan.s = s;
			fan.t = t;
			fan.insert(fan.end(), faces.begin(), faces.end());
			fans.push_back(fan);

#ifdef DISPLAY_ADJACENCY_STATISTICS
			++num_each_sized_fans[fan.size()];
#endif
		}
	}

#ifdef DISPLAY_ADJACENCY_STATISTICS
	std::map<std::size_t, std::size_t>::iterator pos = num_each_sized_fans.begin();
	for (; pos != num_each_sized_fans.end(); ++pos) {
		if (pos->second > 0)
			std::cout << "\t" << pos->first << " - sized fans: " << pos->second << std::endl;
	}
#endif

	vertex_source_planes_.unbind();

	return fans;
}


bool HypothesisGenerator::ready_for_optimization(Map* mesh) const {
	return (
		supporting_planes_.size() > 4 &&
		MapFacetAttribute<double>::is_defined(mesh, Method::facet_attrib_supporting_point_num) &&
		MapFacetAttribute<double>::is_defined(mesh, Method::facet_attrib_facet_area) &&
		MapFacetAttribute<double>::is_defined(mesh, Method::facet_attrib_covered_area)
		);
}
