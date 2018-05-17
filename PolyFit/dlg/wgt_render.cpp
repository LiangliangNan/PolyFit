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

#include "wgt_render.h"
#include "main_window.h"
#include "paint_canvas.h"
#include "../basic/file_utils.h"
#include "../renderer/surface_render.h"
#include "../renderer/point_set_render.h"
#include <QColorDialog>
#include <QFileDialog>



WgtRender::WgtRender(QWidget *parent) 
{
	setupUi(this);

	mainWindow_ = dynamic_cast<MainWindow*>(parent);
	mesh_render_ = mainWindow_->canvas()->mesh_render();
	point_set_render_ = mainWindow_->canvas()->point_set_render();

	//////////////////////////////////////////////////////////////////////////
	
	updatePanel();

	connect(checkBoxPerFaceColor, SIGNAL(toggled(bool)), this, SLOT(setPerFaceColor(bool)));

	// show / hide
	connect(checkBoxSurface, SIGNAL(toggled(bool)), this, SLOT(setShowSurface(bool)));
	connect(checkBoxSharpEdges, SIGNAL(toggled(bool)), this, SLOT(setShowSharpEdges(bool)));
	// color
	connect(toolButtonSurfaceColor, SIGNAL(clicked()), this, SLOT(setSurfaceColor()));
	connect(toolButtonSharpEdgesColor, SIGNAL(clicked()), this, SLOT(setSharpEdgeColor()));

	// width
	connect(doubleSpinBoxSharpEdgesWidth, SIGNAL(valueChanged(double)), this, SLOT(setSharpEdgeWidth(double)));

	//////////////////////////////////////////////////////////////////////////

	connect(checkBoxPerPointColor, SIGNAL(toggled(bool)), this, SLOT(setPerPointColor(bool)));

	// show / hide
	connect(checkBoxPointSet, SIGNAL(toggled(bool)), this, SLOT(setShowPointSet(bool)));
	connect(checkBoxSegments, SIGNAL(toggled(bool)), this, SLOT(setShowSegments(bool)));

	// color
	connect(toolButtonPointSetColor, SIGNAL(clicked()), this, SLOT(setPointSetColor()));

	// size
	connect(doubleSpinBoxVerticesSize, SIGNAL(valueChanged(double)), this, SLOT(setPointsSize(double)));
	connect(doubleSpinBoxSegmentVerticesSize, SIGNAL(valueChanged(double)), this, SLOT(setSegmentVerticesSize(double)));
}


void WgtRender::setPerFaceColor(bool b) {
	mesh_render_->set_per_face_color(b);
	mainWindow_->canvas()->update();
}


void WgtRender::setShowSurface(bool b) {
	SurfaceStyle s = mesh_render_->surface_style();
	s.visible = b;
	mesh_render_->set_surface_style(s);
	mainWindow_->canvas()->update();
}	


void WgtRender::setShowSharpEdges(bool b) {
	EdgeStyle s = mesh_render_->sharp_edge_style();
	s.visible = b;
	mesh_render_->set_sharp_edge_style(s);
	mainWindow_->canvas()->update();
}


void WgtRender::setSurfaceColor() {
	SurfaceStyle s = mesh_render_->surface_style();
	QColor orig(s.color.r() * 255, s.color.g() * 255, s.color.b() * 255, s.color.a() * 255);
	QColor c = QColorDialog::getColor(orig, this);	
	if (c.isValid()) {
		//s.color = Color(c.redF(), c.greenF(), c.blueF(), c.alphaF());
		s.color = Color(c.redF(), c.greenF(), c.blueF(), 0.5f);
		mesh_render_->set_surface_style(s);
		mainWindow_->canvas()->update();

		QPixmap pixmap(25, 19);
		pixmap.fill(c);
		toolButtonSurfaceColor->setIcon(QIcon(pixmap));
	}
}


void WgtRender::setSharpEdgeColor() {
	EdgeStyle s = mesh_render_->sharp_edge_style();
	QColor orig(s.color.r() * 255, s.color.g() * 255, s.color.b() * 255, s.color.a() * 255);
	QColor c = QColorDialog::getColor(orig, this);
	if (c.isValid()) {
		s.color = Color(c.redF(), c.greenF(), c.blueF(), c.alphaF());
		mesh_render_->set_sharp_edge_style(s);
		mainWindow_->canvas()->update();

		QPixmap pixmap(25, 19);
		pixmap.fill(c);
		toolButtonSharpEdgesColor->setIcon(QIcon(pixmap));
	}
}


void WgtRender::setSharpEdgeWidth(double v) {
	EdgeStyle s = mesh_render_->sharp_edge_style();
	s.width = v;
	mesh_render_->set_sharp_edge_style(s);
	mainWindow_->canvas()->update();
}


void WgtRender::updatePanel() {
	checkBoxPerFaceColor->setChecked(mesh_render_->per_face_color());

	const SurfaceStyle& ss = mesh_render_->surface_style();
	checkBoxSurface->setChecked(ss.visible);
	const EdgeStyle& cs = mesh_render_->sharp_edge_style();
	checkBoxSharpEdges->setChecked(cs.visible);

	QPixmap pixmap(25, 19);
	pixmap.fill(QColor(ss.color.r() * 255, ss.color.g() * 255, ss.color.b() * 255, ss.color.a() * 255));
	toolButtonSurfaceColor->setIcon(QIcon(pixmap));
	pixmap.fill(QColor(cs.color.r() * 255, cs.color.g() * 255, cs.color.b() * 255, cs.color.a() * 255));
	toolButtonSharpEdgesColor->setIcon(QIcon(pixmap));

	doubleSpinBoxSharpEdgesWidth->setValue(cs.width);

	//////////////////////////////////////////////////////////////////////////

	checkBoxPerPointColor->setChecked(point_set_render_->per_point_color());

	const PointStyle& ps = point_set_render_->point_set_style();
	checkBoxPointSet->setChecked(ps.visible);
	const PointStyle& vs = point_set_render_->vertex_group_style();
	checkBoxSegments->setChecked(vs.visible);

	pixmap.fill(QColor(ps.color.r() * 255, ps.color.g() * 255, ps.color.b() * 255, ps.color.a() * 255));
	toolButtonPointSetColor->setIcon(QIcon(pixmap));

	doubleSpinBoxVerticesSize->setValue(ps.size);
	doubleSpinBoxSegmentVerticesSize->setValue(vs.size);
}

void WgtRender::setPerPointColor(bool b){
	point_set_render_->set_per_point_color(b);
	mainWindow_->canvas()->update();
}


void WgtRender::setShowPointSet(bool b) {
	PointStyle s = point_set_render_->point_set_style();
	s.visible = b;
	point_set_render_->set_point_set_style(s);
	mainWindow_->canvas()->update();
}


void WgtRender::setShowSegments(bool b){
	PointStyle s = point_set_render_->vertex_group_style();
	s.visible = b;
	point_set_render_->set_vertex_group_style(s);
	mainWindow_->canvas()->update_all();
}



// color
void WgtRender::setPointSetColor(){
	PointStyle s = point_set_render_->point_set_style();
	QColor orig(s.color.r() * 255, s.color.g() * 255, s.color.b() * 255, s.color.a() * 255);
	QColor c = QColorDialog::getColor(orig, this);
	if (c.isValid()) {
		//s.color = Color(c.redF(), c.greenF(), c.blueF(), c.alphaF());
		s.color = Color(c.redF(), c.greenF(), c.blueF(), 1.0f);
		point_set_render_->set_point_set_style(s);
		mainWindow_->canvas()->update();

		QPixmap pixmap(25, 19);
		pixmap.fill(c);
		toolButtonPointSetColor->setIcon(QIcon(pixmap));
	}
}


// size
void WgtRender::setPointsSize(double v){
	PointStyle s = point_set_render_->point_set_style();
	s.size = v;
	point_set_render_->set_point_set_style(s);
	mainWindow_->canvas()->update();
}


void WgtRender::setSegmentVerticesSize(double v) {
	PointStyle s = point_set_render_->vertex_group_style();
	s.size = v;
	point_set_render_->set_vertex_group_style(s);
	mainWindow_->canvas()->update_all();
}

