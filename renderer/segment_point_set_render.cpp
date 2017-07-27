#include "segment_point_set_render.h"
#include "../geom/point_set.h"
#include "../geom/vertex_group.h"


#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif




SegmentPointSetRender::SegmentPointSetRender(PointSet* obj)
: PlainPointSetRender(obj)
{
}


SegmentPointSetRender::~SegmentPointSetRender(void)
{
}


void SegmentPointSetRender::draw() {
	if (vertices_style_.visible) {
		// Liangliang: my experience is that it hurts a lot the rendering performance.
		// If not needed, you can set to "GL_FALSE" (e.g., for large scale scenes).
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

		if (points_as_spheres_ && has_points_shaders_) {
			activate_points_shaders();
			glEnable(GL_MULTISAMPLE);
		}
		else
			glDisable(GL_MULTISAMPLE);

		draw_boundaries();

		draw_segments();

		if (points_as_spheres_ && has_points_shaders_)
			deactivate_points_shaders();
		
		// Liangliang: my experience is that it hurts a lot the rendering performance.
		// If not needed, you can set to "GL_FALSE" (e.g., for large scale scenes).
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

		glEnable(GL_MULTISAMPLE);
		draw_faces();
	}
}


void SegmentPointSetRender::draw_segments() {
	if (!target())
		return;

	int num = target()->num_points();
	if (num < 1)
		return;

	const std::vector<vec3>& points = target()->points();
	const std::vector<vec3>& normals = target()->normals();
	const std::vector<vec3>& colors = target()->colors();
	const std::vector<VertexGroup*>& groups = target()->groups();

	glPointSize(vertices_style_.size);
	if (target()->has_normals()) {
		glEnable(GL_LIGHTING);

		glBegin(GL_POINTS);
		for (std::size_t i = 0; i < groups.size(); ++i) {
			VertexGroup* g = groups[i];
			glColor3fv(g->color().data());

			for (std::size_t j = 0; j < g->size(); ++j) {
				unsigned int idx = g->at(j);
				const vec3& p = points[idx];
				const vec3& n = normals[idx];
				glNormal3fv(n.data());
				glVertex3fv(p.data());
			}
		}
		glEnd();
	}
	else {
		glDisable(GL_LIGHTING); // always off for points without normals

		glBegin(GL_POINTS);
		for (std::size_t i = 0; i < groups.size(); ++i) {
			VertexGroup* g = groups[i];
			glColor3fv(g->color().data());

			for (std::size_t j = 0; j < g->size(); ++j) {
				unsigned int idx = g->at(j);
				const vec3& p = points[idx];
				glVertex3fv(p.data());
			}
		}
		glEnd();
	}
}



void SegmentPointSetRender::draw_faces() {
	if (!target())
		return;

	int num = target()->num_points();
	if (num < 1)
		return;

	const std::vector<vec3>& points = target()->points();
	const std::vector<VertexGroup*>& groups = target()->groups();

	glEnable(GL_LIGHTING);
	for (std::size_t i = 0; i < groups.size(); ++i) {
		VertexGroup* g = groups[i];
		const Plane3d& plane = g->plane();

		glColor3fv(g->color().data());
		glNormal3fv(plane.normal().data());

		const std::vector<unsigned int>& f = g->facet();
		glBegin(GL_POLYGON);
		for (std::size_t j = 0; j < f.size(); ++j) {
			unsigned int idx = f[j];
			const vec3& p = points[idx];
			vec3 q = plane.projection(p);
			glVertex3fv(q.data());
		}
		glEnd();
	}
}