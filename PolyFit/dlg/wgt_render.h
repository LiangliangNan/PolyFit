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

#ifndef WGT_RENDER_H
#define WGT_RENDER_H

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
