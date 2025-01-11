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

#ifndef WGT_RENDER_H
#define WGT_RENDER_H

#include <QWidget>

#include "ui_wgt_render.h"


class MainWindow;
class SurfaceRender;
class PointSetRender;

class WgtRender : public QWidget, public Ui::WidgetRender
{
	Q_OBJECT

public:
	WgtRender(QWidget *parent = 0);
	~WgtRender() {}

	virtual void updatePanel() ;

private Q_SLOTS:

	//---------------------------------------------------------------------

	void setPerFaceColor(bool b);
	void setShowSurface(bool b);
	void setShowSharpEdges(bool b);

	void setSurfaceColor();
	void setSharpEdgeColor();

	void setSharpEdgeWidth(double v);

	//---------------------------------------------------------------------

	void setPerPointColor(bool b);
	void setShowPointSet(bool b);
	void setShowSegments(bool);

	// color
	void setPointSetColor();

	// size
	void setPointsSize(double);
	void setSegmentVerticesSize(double);

private:
	MainWindow*		mainWindow_;

	SurfaceRender*	mesh_render_;
	PointSetRender* point_set_render_;
};

#endif 
