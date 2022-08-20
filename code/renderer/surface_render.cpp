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

#include "surface_render.h"
#include "opengl_info.h"
#include "../basic/logger.h"
#include "../basic/canvas.h"
// these two are here for cylinders
#include "../model/map_geometry.h" 
#include "../math/quaternion.h"

/* this is how we can safely include GLU */
#if defined(__APPLE__) && defined(__MACH__)
#   include <OpenGL/glu.h>
#else
#    include <GL/glu.h>
#endif

static GLUquadric* g_quadric = 0;

SurfaceRender::SurfaceRender(Canvas* cvs) 
	: canvas_(cvs)
{
	per_face_color_ = true;

	surface_style_.visible = true;
	surface_style_.color = Color(1.0f, 1.0f, 0.0f, 1.0f);

	mesh_style_.visible = false;
	mesh_style_.color = Color(0.8f, 0.35f, 0.0f, 1.0f);
	mesh_style_.width = 1;

	sharp_edge_style_.visible = true;
	sharp_edge_style_.color = Color(0.8f, 0.35f, 0.0f, 1.0f);
	sharp_edge_style_.width = 1;

	if (g_quadric == 0)
		g_quadric = gluNewQuadric();
}


SurfaceRender::~SurfaceRender() {
	if (g_quadric)
		gluDeleteQuadric(g_quadric);
}


void SurfaceRender::draw(Map* mesh, bool interacting) {
	if (surface_style_.visible) {
		if (mesh_style_.visible) {
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(0.5f, -0.0001f);
		}

		draw_surface(mesh);

		if (mesh_style_.visible)
			glDisable(GL_POLYGON_OFFSET_FILL);
	}

	if (mesh_style_.visible)
		draw_mesh(mesh);

	if (sharp_edge_style_.visible)
		draw_corner_edges(mesh, interacting);
}

const SurfaceStyle& SurfaceRender::surface_style() const {
	return surface_style_;
}

void SurfaceRender::set_surface_style(const SurfaceStyle& x) {
	surface_style_ = x;
}

const EdgeStyle& SurfaceRender::mesh_style() const {
	return mesh_style_;
}

void SurfaceRender::set_mesh_style(const EdgeStyle& x) {
	mesh_style_ = x;
}

const EdgeStyle& SurfaceRender::sharp_edge_style() const {
	return sharp_edge_style_;
}

void SurfaceRender::set_sharp_edge_style(const EdgeStyle& x) {
	sharp_edge_style_ = x;
}

void SurfaceRender::draw_surface(Map* mesh) {
	if (per_face_color_)
		facet_color_.bind_if_defined(mesh, "color");

	glEnable(GL_MULTISAMPLE);

	glEnable(GL_LIGHTING);
	glColor4fv(surface_style_.color.data());

	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	if (!MapFacetNormal::is_defined(mesh))
		mesh->compute_facet_normals();
	MapFacetNormal normal(mesh);

	FOR_EACH_FACET_CONST(Map, mesh, it) {
		const vec3& n = normal[it];

		glNormal3fv(n.data());
		if (facet_color_.is_bound())
			glColor4fv(facet_color_[it].data());

		glBegin(GL_POLYGON);
		Map::Halfedge* jt = it->halfedge();
		do {
			const vec3& p = jt->vertex()->point();
			glVertex3fv(p.data());
			jt = jt->next();
		} while (jt != it->halfedge());
		glEnd();
	}

	if (facet_color_.is_bound()) {
		facet_color_.unbind();
	}
}


void SurfaceRender::draw_mesh(Map* mesh) {
	glColor4fv(mesh_style_.color.data());

	glDisable(GL_LIGHTING);
	glLineWidth(mesh_style_.width);
	glBegin(GL_LINES);
	FOR_EACH_EDGE_CONST(Map, mesh, it) {
//		if (it->is_border_edge())
//			continue;
		const vec3& p = it->vertex()->point();
		const vec3& q = it->opposite()->vertex()->point();
		glVertex3fv(p.data());
		glVertex3fv(q.data());
	}
	glEnd();
}


void SurfaceRender::draw_corner_edges(Map* mesh, bool interacting) {
	glColor4fv(sharp_edge_style_.color.data());

    MapHalfedgeAttribute<bool> edge_is_sharp(mesh, "SharpEdge");

	if (interacting) {
		glDisable(GL_LIGHTING);
		glLineWidth(sharp_edge_style_.width + 1.0f);
		glBegin(GL_LINES);
		FOR_EACH_EDGE_CONST(Map, mesh, it) {
            if (edge_is_sharp[it] || it->is_border_edge()) {
                const vec3& t = it->vertex()->point();
                glVertex3fv(t.data());
                const vec3& s = it->opposite()->vertex()->point();
                glVertex3fv(s.data());
            }
		}
		glEnd();
	}
	else {
		const vec3& c = mesh->bbox().center();
		float ratio = canvas_->get_camera()->pixelGLRatio(c.x, c.y, c.z);
		float r = (sharp_edge_style_.width + 1) * ratio;

		int slices = 20;

		// current implementation uses gluCylinder()
		// TODO: add "segment as cylinder" support (using GLSL).

		glEnable(GL_LIGHTING);

		glShadeModel(GL_SMOOTH);

		FOR_EACH_EDGE_CONST(Map, mesh, it) {
            if (edge_is_sharp[it] || it->is_border_edge()) {
                glPushMatrix();
                const vec3& t = it->vertex()->point();
                glTranslated(t.x, t.y, t.z);
                gluSphere(g_quadric, r, slices, slices);
                glPopMatrix();

                glPushMatrix();
                const vec3& s = it->opposite()->vertex()->point();
                glTranslated(s.x, s.y, s.z);
                gluSphere(g_quadric, r, slices, slices);
                glMultMatrixd(Quaternion(vec3(0, 0, 1), Geom::vector(it)).matrix());
                gluCylinder(g_quadric, r, r, Geom::edge_length(it), slices, 1);
                glPopMatrix();
            }
		}
	}
}


void SurfaceRender::set_per_face_color(bool x) {
	per_face_color_ = x;
}

