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

#ifndef PAINTCANVAS_H
#define PAINTCANVAS_H

#include <GL/glew.h>

#include "../3rd_party/QGLViewer/QGLViewer/qglviewer.h"
#include "../basic/color.h"
#include "../math/math_types.h"
#include "../model/point_set.h"
#include "../model/map.h"

class MainWindow;
class SurfaceRender;
class PointSetRender;
class HypothesisGenerator;

class PaintCanvas : public QGLViewer
{
	Q_OBJECT

public:
	PaintCanvas(QWidget *parent);
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
