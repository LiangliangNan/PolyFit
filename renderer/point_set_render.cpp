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



#include "point_set_render.h"
#include "../model/point_set.h"
#include "../model/vertex_group.h"

#include <GL/glew.h>



PointSetRender::PointSetRender(Canvas* cvs)
	: canvas_(cvs)
{
	per_point_color_ = false;

	point_set_style_.visible = true;
	point_set_style_.color = Color(85 / 255.0f, 170 / 255.0f, 1.0f);
	point_set_style_.size = 3;

	vertex_group_style_.visible = true;
	vertex_group_style_.color = Color(85 / 255.0f, 170 / 255.0f, 1.0f);
	vertex_group_style_.size = 3;
}


PointSetRender::~PointSetRender(void)
{
}

const PointStyle& PointSetRender::point_set_style() const {
	return point_set_style_;
}

void PointSetRender::set_point_set_style(const PointStyle& x) {
	point_set_style_ = x;
}

const PointStyle& PointSetRender::vertex_group_style() const {
	return vertex_group_style_;
}

void PointSetRender::set_vertex_group_style(const PointStyle& x) {
	vertex_group_style_ = x;
}

void PointSetRender::set_per_point_color(bool x) {
	per_point_color_ = x;
}


void PointSetRender::draw(PointSet*	pset) {
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	glDisable(GL_MULTISAMPLE);

	if (point_set_style_.visible) {
		if (per_point_color_ && pset->has_colors())
			draw_point_set_per_point_color(pset);
		else
			draw_point_set_uniform_color(pset);
	}

	if (vertex_group_style_.visible)
		draw_vertex_groups(pset);

	glEnable(GL_MULTISAMPLE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}


void PointSetRender::draw_point_set_per_point_color(PointSet* pset) {
	if (!pset)
		return;

	if (!pset->has_colors())
		return;

	int num = pset->num_points();
	if (num < 1)
		return;

	float* points = &(pset->points()[0].x);
	float* colors = &(pset->colors()[0].x);

	glPointSize(point_set_style_.size);

	if (pset->has_normals()) {
		float* normals = &(pset->normals()[0].x);
		glEnable(GL_LIGHTING);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, points);
		glNormalPointer(GL_FLOAT, 0, normals);
		glColorPointer(3, GL_FLOAT, 0, colors);
		glDrawArrays(GL_POINTS, 0, num);
		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
	}
	else {
		glDisable(GL_LIGHTING); // always off for points without normals

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, points);
		glColorPointer(3, GL_FLOAT, 0, colors);
		glDrawArrays(GL_POINTS, 0, num);
		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
		glEnable(GL_LIGHTING);
	}
}


void PointSetRender::draw_point_set_uniform_color(PointSet* pset) {
	if (!pset)
		return;

	int num = pset->num_points();
	if (num < 1)
		return;

	float* points = &(pset->points()[0].x);

	glPointSize(point_set_style_.size);
	glColor3fv(point_set_style_.color.data());
	if (pset->has_normals()) {
		float* normals = &(pset->normals()[0].x);
		glEnable(GL_LIGHTING);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, points);
		glNormalPointer(GL_FLOAT, 0, normals);
		glDrawArrays(GL_POINTS, 0, num);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
	}
	else {
		glDisable(GL_LIGHTING); // always off for points without normals

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, points);
		glDrawArrays(GL_POINTS, 0, num);
		glDisableClientState(GL_VERTEX_ARRAY);
		glEnable(GL_LIGHTING);
	}
}

void PointSetRender::draw_vertex_groups(PointSet* pset) {
	if (!pset)
		return;

	int num = pset->num_points();
	if (num < 1)
		return;

	const std::vector<vec3>& points = pset->points();
	const std::vector<vec3>& normals = pset->normals();
	const std::vector<vec3>& colors = pset->colors();
	const std::vector<VertexGroup::Ptr>& groups = pset->groups();

	glPointSize(vertex_group_style_.size);
	if (pset->has_normals()) {
		glEnable(GL_LIGHTING);

		glBegin(GL_POINTS);
		for (std::size_t i = 0; i < groups.size(); ++i) {
			VertexGroup* g = groups[i];
			if (!g->is_visible())
				continue;

			if (g->is_highlighted())
				glColor3f(0.0f, 1.0f, 1.0f);
			else
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
			if (!g->is_visible())
				continue;

			if (g->is_highlighted())
				glColor3f(0.0f, 1.0f, 1.0f);
			else
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

