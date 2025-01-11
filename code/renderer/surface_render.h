/* ---------------------------------------------------------------------------
 * Copyright (C) 2017 Liangliang Nan <liangliang.nan@gmail.com>
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of PolyFit. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 *
 *     Liangliang Nan and Peter Wonka.
 *     PolyFit: Polygonal Surface Reconstruction from Point Clouds.
 *     ICCV 2017.
 *
 *  For more information:
 *  https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html
 * ---------------------------------------------------------------------------
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
