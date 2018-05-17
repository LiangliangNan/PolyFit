
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


#ifndef _POINT_SET_RENDERER_H_
#define _POINT_SET_RENDERER_H_

#include "renderer_common.h"
#include "rendering_styles.h"


class Canvas;
class PointSet;

class RENDERER_API PointSetRender
{
public:
	PointSetRender(Canvas* cvs);
	~PointSetRender(void);

	bool per_point_color() const { return per_point_color_; }
	void set_per_point_color(bool x);

	const PointStyle& point_set_style() const;
	void set_point_set_style(const PointStyle& x);

	const PointStyle& vertex_group_style() const;
	void set_vertex_group_style(const PointStyle& x);

	virtual void draw(PointSet*	pset);

protected:
	// whole point set
	virtual void draw_point_set_per_point_color(PointSet* pset);
	virtual void draw_point_set_uniform_color(PointSet* pset);

	// segments (vertex groups)
	virtual void draw_vertex_groups(PointSet* pset);

protected:
	Canvas*			canvas_;

	PointStyle		point_set_style_;
	PointStyle		vertex_group_style_;

	bool			per_point_color_;
};



#endif