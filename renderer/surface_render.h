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

#ifndef _RENDERER_SURFACE_RENDERER_H_
#define _RENDERER_SURFACE_RENDERER_H_

#include "renderer_common.h"
#include "rendering_styles.h"
#include "../model/map_attributes.h"


class Map;
class Canvas;

class RENDERER_API SurfaceRender
{
public:
	SurfaceRender(Canvas* cvs) ;
	~SurfaceRender();

	virtual void draw(Map* mesh, bool interacting);

	//___________________________________________________________

	bool per_face_color() const { return per_face_color_; }
	void set_per_face_color(bool x);

	bool smooth_shading() const ;
	void set_smooth_shading(bool x) ;

	const SurfaceStyle& surface_style() const ;
	void set_surface_style(const SurfaceStyle& x) ;

	const EdgeStyle& mesh_style() const ;
	void set_mesh_style(const EdgeStyle& x) ;

	const EdgeStyle& sharp_edge_style() const;
	void set_sharp_edge_style(const EdgeStyle& x);

protected:
	virtual void draw_surface(Map* mesh);
	virtual void draw_mesh(Map* mesh);
	virtual void draw_corner_edges(Map* mesh, bool interacting);

protected:
	Canvas* canvas_;

	bool	per_face_color_;
	MapFacetAttribute<Color>	facet_color_;

	SurfaceStyle surface_style_ ;

	EdgeStyle    mesh_style_;
	EdgeStyle	 sharp_edge_style_;
} ;



#endif
