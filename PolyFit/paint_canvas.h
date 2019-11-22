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

#ifndef PAINTCANVAS_H
#define PAINTCANVAS_H

#include <GL/glew.h>

#include "../3rd_QGLViewer-2.6.3/qglviewer.h"
#include "../basic/color.h"
#include "../math/math_types.h"
#include "../basic/canvas.h"
#include "../model/point_set.h"
#include "../model/map.h"
#include "../math/linear_program_solver.h"


class MainWindow;
class SurfaceRender;
class PointSetRender;
class HypothesisGenerator;

class PaintCanvas : public QGLViewer
{
	Q_OBJECT

public:
	PaintCanvas(QWidget *parent, QGLFormat format);
	~PaintCanvas();

public:

	void update_graphics();
	void update_all();

	//////////////////////////////////////////////////////////////////////////

	MainWindow* mainWindow() const { return main_window_; }

	vec2 projectionOf(const vec3& p);          // point to screen
	vec3 unProjectionOf(double winx, double winy, double winz);  // screen to point

	//////////////////////////////////////////////////////////////////////////
	
	void setMesh(Map* mesh);
	void setPointSet(PointSet* pset);

	// the active object
	Map*		hypothesisMesh() const ;
	Map*		optimizedMesh() const;
	PointSet*	pointSet() const;
	
	SurfaceRender*	mesh_render() const { return mesh_render_; }
	PointSetRender* point_set_render() const { return point_set_render_; }

	void snapshotScreen(const QString& fileName);

	void clear();

	//////////////////////////////////////////////////////////////////////////

protected:
	virtual void draw();
	virtual void init();

	// Keyboard events functions
	virtual void keyPressEvent(QKeyEvent *e);

public Q_SLOTS:
	void fitScreen() ;

	void showCoordinateSystem(bool);

	//////////////////////////////////////////////////////////////////////////

	void refinePlanes();
	void generateFacetHypothesis();
	void generateQualityMeasures();
	void optimization();

	void setShowInput(bool);
	void setShowCandidates(bool);
	void setShowResult(bool);

	void saveStateAsMappleFormat();

private :
	void drawCornerAxis();

protected:
	MainWindow*	main_window_;
	vec3		light_pos_;

	int     coord_system_region_size_;
	bool	show_coord_sys_;

	Map::Ptr		hypothesis_mesh_;
	Map::Ptr		optimized_mesh_;
	PointSet::Ptr	point_set_;

	bool	show_input_;
	bool	show_candidates_;
	bool	show_result_;

	SurfaceRender*	mesh_render_;
	PointSetRender* point_set_render_;

	HypothesisGenerator* hypothesis_;

	bool		show_hint_text_;
	QString     hint_text_;
	QString     hint_text2nd_;

	bool        show_mouse_hint_;
};


#endif // PAINTCANVAS_H
