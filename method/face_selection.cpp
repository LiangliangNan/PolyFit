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


#include "face_selection.h"
#include "method_global.h"
#include "../model/point_set.h"
#include "../model/map_geometry.h"
#include "../basic/logger.h"
#include "../model/map_editor.h"

#include <algorithm>


// For debugging (identifying topological issues)
//#define DISPLAY_ADJACENCY_STATISTICS


FaceSelection::FaceSelection(PointSet* pset, Map* model)
	: pset_(pset)
	, model_(model)
{
	Logger::out("-") << "tolerance for coincident vertices: " << Method::coincident_threshold << std::endl;

	FOR_EACH_EDGE(Map, model, it) {
		double len = Geom::edge_length(it);
		if (len <= Method::coincident_threshold) {
			std::cout << "---[Purification] very short edge detected. length: " << len << std::endl;
		}
	}

}


FaceSelection::~FaceSelection()
{
}


void FaceSelection::prune_boundary_faces(Map* model) {
	std::vector<MapTypes::Facet*> boundary_faces;

	FOR_EACH_FACET(Map, model, it) {
		MapTypes::Halfedge* h = it->halfedge();
		do {
			if (h->is_border_edge()) {
				boundary_faces.push_back(it);
				break;
			}
			h = h->next();
		} while (h != it->halfedge());
	}

	MapEditor editor(model);
	for (std::size_t i = 0; i < boundary_faces.size(); ++i) {
		Map::Facet* f = boundary_faces[i];
		editor.erase_facet(f->halfedge());
	}
	Logger::out("-") << boundary_faces.size() << " boundary faces pruned" << std::endl;
}


void FaceSelection::IntersectionAdjacency::init(Map* model, const std::vector<Plane3d*>& supporting_planes) {
	edge_container_.clear();
	plane_index_.clear();
	vertex_source_planes_.bind_if_defined(model, "VertexSourcePlanes");

	int num = int(supporting_planes.size());
	for (int i = 0; i < num; ++i) {
		Plane3d* plane = supporting_planes[i];
		plane_index_[plane] = i;
	}

	FOR_EACH_HALFEDGE(Map, model, it) {
		if (it->facet() == nil)
			continue;

		Map::Vertex* s = it->vertex();
		Map::Vertex* t = it->opposite()->vertex();

		const std::set<Plane3d*>& set_s = vertex_source_planes_[s];
		const std::set<Plane3d*>& set_t = vertex_source_planes_[t];

		std::vector<int> s_indices, t_indices;
		std::set<Plane3d*>::const_iterator pos = set_s.begin();
		for (; pos != set_s.end(); ++pos) {
			int idx = plane_index_[*pos];
			s_indices.push_back(idx);
		}
		pos = set_t.begin();
		for (; pos != set_t.end(); ++pos) {
			int idx = plane_index_[*pos];
			t_indices.push_back(idx);
		}
		assert(s_indices.size() == 3);
		assert(t_indices.size() == 3);
		std::sort(s_indices.begin(), s_indices.end());
		std::sort(t_indices.begin(), t_indices.end());

		int idx_s = ((s_indices[0] * num) + s_indices[1]) * num + s_indices[2];
		int idx_t = ((t_indices[0] * num) + t_indices[1]) * num + t_indices[2];

		if (idx_s > idx_t)
			ogf_swap(idx_s, idx_t);

		edge_container_[idx_s][idx_t].insert(it);
	}

	vertex_source_planes_.unbind();
}


std::vector<FaceStar> FaceSelection::IntersectionAdjacency::extract(Map* model, const std::vector<Plane3d*>& supporting_planes) {
	init(model, supporting_planes);

#ifdef DISPLAY_ADJACENCY_STATISTICS
	std::map<std::size_t, std::size_t>   num_each_sized_fans;
	for (std::size_t i = 1; i < 20; ++i)
		num_each_sized_fans[i] = 0;
#endif

	std::vector<FaceStar> fans;
	std::map< std::size_t, std::map< std::size_t, std::set<MapTypes::Halfedge*> > >::const_iterator it = edge_container_.begin();
	for (; it != edge_container_.end(); ++it) {
		const std::map < std::size_t, std::set<MapTypes::Halfedge*> >& tmp = it->second;
		std::map< std::size_t, std::set<MapTypes::Halfedge*> >::const_iterator cur = tmp.begin();
		for (; cur != tmp.end(); ++cur) {
			const std::set<MapTypes::Halfedge*>& edges = cur->second;
			FaceStar fan;
			fan.insert(fan.end(), edges.begin(), edges.end());
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

	return fans;
}
