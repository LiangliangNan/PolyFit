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


#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>

#include "ui_main_window.h"
#include "../math/math_types.h"
#include "../basic/logger.h"
#include "../basic/progress.h"
#include "../math/linear_program_solver.h"

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
	MainWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
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
	bool save();
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
	void dragEnterEvent(QDragEnterEvent *e);
	void dropEvent(QDropEvent *e);
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
