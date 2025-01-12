
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


#ifndef _POINT_SET_RENDERER_H_
#define _POINT_SET_RENDERER_H_

#include <renderer/renderer_common.h>
#include <renderer/rendering_styles.h>


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