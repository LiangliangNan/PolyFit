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


#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>

#include <math/math_types.h>
#include <basic/logger.h>
#include <basic/progress.h>
#include <math/linear_program_solver.h>

#include "ui_main_window.h"

class QLabel;
class QComboBox;
class PaintCanvas;
class QProgressBar;
class WeightPanelClick;
class WeightPanelManual;
class WgtRender;

class MainWindow 
	: public QMainWindow
	, public LoggerClient
	, public ProgressClient
	, public Ui::PolyFitClass
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();

	PaintCanvas* canvas() { return mainCanvas_; }

	virtual void out_message(const std::string& msg);
	virtual void warn_message(const std::string& msg);
	virtual void err_message(const std::string& msg);
	virtual void status_message(const std::string& msg, int timeout);
	virtual void notify_progress(std::size_t value);

	void updateWeights();
	void disableActions(bool b);

	void defaultRenderingForCandidates();
	void defaultRenderingForResult();

	LinearProgramSolver::SolverName active_solver() const;

public Q_SLOTS:
	bool open();
	bool saveReconstruction();
	bool saveCandidateFaces();

	void updateStatusBar();

	void snapshotScreen();

	void resetWeights();
	void setManualInputWeights(bool);

	void about();

private:
	void createActions(); 
	void createStatusBar();
	void createToolBar();

	void createRenderingPanel();

	void readSettings();
	void writeSettings();
	
	bool doOpen(const QString &fileName);
	bool doSave(const QString &fileName);

	void setCurrentFile(const QString &fileName);
	
	QString strippedName(const QString &fullFileName);

protected:
	void closeEvent(QCloseEvent *e);

private:
	PaintCanvas*	mainCanvas_;

	QString			pointCloudFileName_;
	QString			hypothesisMeshFileName_;
	QString			optimizedMeshFileName_;
	QString			curDataDirectory_;
	QString			curCameraConfigFileDirectory_;

	QProgressBar*	progress_bar_;
	QComboBox*		solverBox_;

	QLabel *statusLabel_,
		*numPointsLabel_,
		*numHypoFacesLabel_,
		*numOptimizedFacesLabel_;

	//////////////////////////////////////////////////////////////////////////

	WgtRender*	wgtRender_;

	WeightPanelClick*	panelClick_;
	WeightPanelManual*	panelManual_;

	float default_fitting_;
	float default_coverage_;
	float default_complexity_;
};

#endif // TESTQGLVIEWER_H
