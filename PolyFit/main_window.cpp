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

#include <QMessageBox>
#include <QFileDialog>
#include <QLabel>
#include <QStatusBar>
#include <QSettings>
#include <QCloseEvent>
#include <QPlainTextEdit>
#include <QGroupBox>
#include <QColorDialog>
#include <QProgressBar>
#include <QMimeData>
#include <QComboBox>

#include "main_window.h"
#include "paint_canvas.h"

#include "dlg/wgt_render.h"
#include "dlg/weight_panel_click.h"
#include "dlg/weight_panel_manual.h"

#include "../basic/logger.h"
#include "../basic/file_utils.h"
#include "../model/map.h"
#include "../model/point_set.h"
#include "../model/map_attributes.h"
#include "../model/map_copier.h"
#include "../model/map_builder.h"
#include "../model/map_io.h"
#include "../model/map_enumerator.h"
#include "../model/point_set_io.h"
#include "../method/method_global.h"
#include "../basic/attribute_serializer.h"


MainWindow::MainWindow(QWidget *parent, Qt::WindowFlags flags)
: QMainWindow(parent, flags)
, curDataDirectory_(".")
{	
	setupUi(this);

	//////////////////////////////////////////////////////////////////////////

	Logger::initialize();
	Logger::instance()->register_client(this);
	Logger::instance()->set_value(Logger::LOG_REGISTER_FEATURES, "*"); // log everything
	Logger::instance()->set_value(Logger::LOG_FILE_NAME, "PolyFit.log");
	//Logger::instance()->set_value("log_features",
	//	"EigenSolver;MapBuilder;MapParameterizer\
	//	LinearSolver");

	// Liangliang: added the time stamp in the log file
	std::string tstr = String::from_current_time();
	Logger::out("") << "--- started at: " << tstr << " ---" << std::endl;
	
	Progress::instance()->set_client(this) ;

	AttributeSerializer::initialize();

	register_attribute_type<int>("int");
	register_attribute_type<float>("float");
	register_attribute_type<double>("double");
	register_attribute_type<bool>("bool") ;
	register_attribute_type<std::string>("string") ;
	register_attribute_type<vec2>("vec2") ;
	register_attribute_type<vec3>("vec3") ;
	register_attribute_type<vec4>("vec4") ;
	register_attribute_type<mat2>("mat2") ;
	register_attribute_type<mat3>("mat3") ;
	register_attribute_type<mat4>("mat4") ;
	register_attribute_type<Color>("Color");

	// ensure backward compatibility with .eobj files generated before.
	// PointXd/VectorXd do not exist anymore.
	register_attribute_type_alias("Vector2d", "vec2") ;
	register_attribute_type_alias("Vector3d", "vec3") ;
	register_attribute_type_alias("Point2d", "vec2") ;
	register_attribute_type_alias("Point3d", "vec3") ;

	//////////////////////////////////////////////////////////////////////////

	// Setup the format to allow anti-aliasing if the graphic driver allows this.
	QGLFormat format = QGLFormat::defaultFormat();
	format.setProfile(QGLFormat::CompatibilityProfile);
	format.setSampleBuffers(true); // you can also call setOption(QGL::SampleBuffers)
	format.setSamples(8);  // 8 is enough

	mainCanvas_ = new PaintCanvas(this, format);
	mainCanvas_->setAttribute(Qt::WA_MouseTracking);
	mainCanvas_->setMouseTracking(true);
	layoutCanvas->addWidget(mainCanvas_);

	//////////////////////////////////////////////////////////////////////////

	setWindowState(Qt::WindowMaximized);
	setFocusPolicy(Qt::StrongFocus);
	showMaximized();

	createRenderingPanel();

	createActions();
	createStatusBar();
	createToolBar();

	readSettings();
	setWindowTitle("PolyFit");

	setContextMenuPolicy(Qt::CustomContextMenu);
	setWindowIcon(QIcon(":/Resources/PolyFit.png"));

	setAcceptDrops(true);
	disableActions(true);
}


MainWindow::~MainWindow()
{
	if (wgtRender_)	delete wgtRender_;

	//////////////////////////////////////////////////////////////////////////

	AttributeSerializer::terminate();
	Progress::instance()->set_client(nil) ;
	Logger::instance()->unregister_client(this);
	Logger::terminate();
}


void MainWindow::out_message(const std::string& msg) {
	plainTextEditOutput->moveCursor(QTextCursor::End);
	plainTextEditOutput->insertPlainText(QString::fromStdString(msg));
	plainTextEditOutput->repaint();
	plainTextEditOutput->update();
}


void MainWindow::warn_message(const std::string& msg) {
	plainTextEditOutput->moveCursor(QTextCursor::End);
	plainTextEditOutput->insertPlainText(QString::fromStdString(msg));
	plainTextEditOutput->repaint();
	plainTextEditOutput->update();
}


void MainWindow::err_message(const std::string& msg) {
	plainTextEditOutput->moveCursor(QTextCursor::End);
	plainTextEditOutput->insertPlainText(QString::fromStdString(msg));
	plainTextEditOutput->repaint();
	plainTextEditOutput->update();
}


void MainWindow::status_message(const std::string& msg, int timeout) {
	statusBar()->showMessage(QString::fromStdString(msg), timeout);
}


void MainWindow::notify_progress(std::size_t value) {
	progress_bar_->setValue(value);
	progress_bar_->setTextVisible(value != 0);
	mainCanvas_->update_all();
}


void MainWindow::dragEnterEvent(QDragEnterEvent *e) {
	if (e->mimeData()->hasUrls()) {
		e->acceptProposedAction();
	}
}

void MainWindow::dropEvent(QDropEvent *e) {
    foreach (const QUrl &url, e->mimeData()->urls()) {
        const QString &fileName = url.toLocalFile();
        doOpen(fileName);
    }
}


void MainWindow::createActions() {
	connect(actionOpen, SIGNAL(triggered()), this, SLOT(open()));
	connect(actionSave, SIGNAL(triggered()), this, SLOT(save()));

	connect(actionSnapshot, SIGNAL(triggered()), this, SLOT(snapshotScreen()));

	connect(actionRefinePlanes, SIGNAL(triggered()), mainCanvas_, SLOT(refinePlanes()));
	connect(actionGenerateFacetHypothesis, SIGNAL(triggered()), mainCanvas_, SLOT(generateFacetHypothesis()));
	connect(actionGenerateQualityMeasures, SIGNAL(triggered()), mainCanvas_, SLOT(generateQualityMeasures()));
	connect(actionOptimization, SIGNAL(triggered()), mainCanvas_, SLOT(optimization()));

	connect(checkBoxShowInput, SIGNAL(toggled(bool)), mainCanvas_, SLOT(setShowInput(bool)));
	connect(checkBoxShowCandidates, SIGNAL(toggled(bool)), mainCanvas_, SLOT(setShowCandidates(bool)));
	connect(checkBoxShowResult, SIGNAL(toggled(bool)), mainCanvas_, SLOT(setShowResult(bool)));

	wgtRender_ = new WgtRender(this);
	layoutRenderer->addWidget(wgtRender_);

	// about menu
	connect(actionAbout, SIGNAL(triggered()), this, SLOT(about()));
}


void MainWindow::createRenderingPanel() {
	default_fitting_ = truncate_digits(Method::lambda_data_fitting, 3);
	default_coverage_ = truncate_digits(Method::lambda_model_coverage, 3);
	default_complexity_ = truncate_digits(Method::lambda_model_complexity, 3);

	panelClick_ = new WeightPanelClick(this);
	panelManual_ = new WeightPanelManual(this);

	verticalLayoutWeights->addWidget(panelClick_);
	verticalLayoutWeights->addWidget(panelManual_);
	panelManual_->setVisible(false);

	connect(panelClick_, SIGNAL(weights_changed()), panelManual_, SLOT(updateUI()));

	connect(pushButtonDefaultWeight, SIGNAL(pressed()), this, SLOT(resetWeights()));
	connect(checkBoxManualInputWeights, SIGNAL(toggled(bool)), this, SLOT(setManualInputWeights(bool)));
}


void MainWindow::updateWeights() {
	if (panelManual_->isEnabled()) {
		panelManual_->updateWeights();
		panelClick_->updateUI();
	}
}


void MainWindow::disableActions(bool b) {
	actionRefinePlanes->setDisabled(b);
	actionGenerateFacetHypothesis->setDisabled(b);
	actionGenerateQualityMeasures->setDisabled(b);
	actionOptimization->setDisabled(b);
}


void MainWindow::defaultRenderingForCandidates() {
	wgtRender_->checkBoxPerFaceColor->setChecked(true);
}


void MainWindow::defaultRenderingForResult() {
	wgtRender_->checkBoxPerFaceColor->setChecked(false);
}


void MainWindow::resetWeights() {
	Method::lambda_data_fitting = default_fitting_;
	Method::lambda_model_coverage = default_coverage_;
	Method::lambda_model_complexity = default_complexity_;

	panelClick_->updateUI();
	panelManual_->updateUI();
}


void MainWindow::setManualInputWeights(bool b) {
	if (b) {
		panelClick_->setVisible(false);
		panelManual_->setVisible(true);
	}
	else {
		panelClick_->setVisible(true);
		panelManual_->setVisible(false);
	}
}


void MainWindow::closeEvent(QCloseEvent *event)
{
	writeSettings();
	event->accept();
}

void MainWindow::createStatusBar()
{	
	statusLabel_ = new QLabel("Ready");
	statusLabel_->setFixedWidth(scrollArea->width());
	statusLabel_->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	statusBar()->addWidget(statusLabel_, 1);

	QLabel* space1 = new QLabel;
	statusBar()->addWidget(space1, 1);

	int length = 200;
	numPointsLabel_ = new QLabel;
	numPointsLabel_->setFixedWidth(length);
	numPointsLabel_->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	statusBar()->addPermanentWidget(numPointsLabel_, 1);

	numHypoFacesLabel_ = new QLabel;
	numHypoFacesLabel_->setFixedWidth(length);
	numHypoFacesLabel_->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	statusBar()->addPermanentWidget(numHypoFacesLabel_, 1);

	numOptimizedFacesLabel_ = new QLabel;
	numOptimizedFacesLabel_->setFixedWidth(length);
	numOptimizedFacesLabel_->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	statusBar()->addPermanentWidget(numOptimizedFacesLabel_, 1);

	QLabel* space2 = new QLabel;
	statusBar()->addWidget(space2, 1);

	//////////////////////////////////////////////////////////////////////////

	progress_bar_ = new QProgressBar;
	progress_bar_->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
	progress_bar_->setFixedWidth(400);
	statusBar()->addPermanentWidget(progress_bar_, 1);

	//////////////////////////////////////////////////////////////////////////

	updateStatusBar();
}


void MainWindow::createToolBar()
{
	solverBox_ = new QComboBox(this);
	solverBox_->setFixedHeight(23);
	solverBox_->setEditable(false);
#ifdef HAS_GUROBI
    solverBox_->addItem("GUROBI");
#endif
    solverBox_->addItem("SCIP");
	solverBox_->addItem("GLPK");
	solverBox_->addItem("LPSOLVE");

	QLabel* label = new QLabel(this);
	label->setText("    Solver");
	label->setAlignment(Qt::AlignLeft);

	QVBoxLayout* layout = new QVBoxLayout;
	layout->addWidget(solverBox_);
	layout->addWidget(label);

	QWidget* widget = new QWidget(this);
	widget->setLayout(layout);

	toolBar->insertWidget(actionRefinePlanes, widget);

	toolBar->insertSeparator(actionRefinePlanes);
}


void MainWindow::updateStatusBar()
{
	QString points = "#points: 0";
	QString hypoFaces = "#faces(hypo): 0";
	QString optimizedFaces = "#faces(result): 0";

	if (mainCanvas_->pointSet()) {
		points = QString("#points: %1").arg(mainCanvas_->pointSet()->num_points());
	}
	if (mainCanvas_->hypothesisMesh()) {
		hypoFaces = QString("#faces(candidates): %1").arg(mainCanvas_->hypothesisMesh()->size_of_facets());
	}
	if (mainCanvas_->optimizedMesh()) {
		// I need to report the number of planar faces, instead of the original face candidates
		//optimizedFaces = QString("#faces(final): %1").arg(mainCanvas_->optimizedMesh()->size_of_facets());
		MapFacetAttribute<int> attrib(mainCanvas_->optimizedMesh());
		int num = MapEnumerator::enumerate_planar_components(mainCanvas_->optimizedMesh(), attrib);
		optimizedFaces = QString("#faces(final): %1").arg(num);
	}

	numPointsLabel_->setText(points);
	numHypoFacesLabel_->setText(hypoFaces);
	numOptimizedFacesLabel_->setText(optimizedFaces);
}


void MainWindow::about()
{
#if defined (ENV_32_BIT)
	QString title = QMessageBox::tr("<h3>PolyFit (32-bit)</h3>");
#elif defined (ENV_64_BIT)
	QString title = QMessageBox::tr("<h3>PolyFit (64-bit)</h3>");
#else
	QString title = QMessageBox::tr("<h3>PolyFit"</h3>);
#endif

#ifndef NDEBUG
	title += QMessageBox::tr(" (Debug Version)");
#endif

	QString text = QMessageBox::tr(
		"<p>PolyFit implements our <span style=\"font-style:italic;\">hypothesis and selection</span> based surface reconstruction method described in the following paper:</p>"
		
		"<p>--------------------------------------------------------------------------<br>"
		"<span style=\"font-style:italic;\">Liangliang Nan</span> and <span style=\"font-style:italic;\">Peter Wonka</span>.<br>"
		"<a href=\"https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html\">PolyFit: Polygonal Surface Reconstruction from Point Clouds.</a><br>"
		"ICCV 2017.<br>"
		"--------------------------------------------------------------------------</p>"

		"<p>Extract planes? You can use my <a href=\"https://3d.bk.tudelft.nl/liangliang/software.html\">Mapple</a> program for plane extraction. Please refer to the ReadMe files for more details.</p>"

		"<p>For comments, suggestions, or any issues, please contact me at<br>"
		"<a href=\"mailto:liangliang.nan@gmail.com\">liangliang.nan@gmail.com</a>.</p>"
		"<p>Liangliang Nan<br>"
		"<a href=\"https://3d.bk.tudelft.nl/liangliang/\">https://3d.bk.tudelft.nl/liangliang/</a><br>"
		"@July.18, 2017</p>"
	);

	QMessageBox::about(this, "About PolyFit", title + text);
}

void MainWindow::readSettings()
{
	QSettings settings("LiangliangNan", "PolyFit");
	curDataDirectory_ = settings.value("currentDirectory").toString();	
}

void MainWindow::writeSettings()
{
	QSettings settings("LiangliangNan", "PolyFit");
	settings.setValue("currentDirectory", curDataDirectory_);
}


void MainWindow::setCurrentFile(const QString &fileName)
{
	curDataDirectory_ = fileName.left(fileName.lastIndexOf("/") + 1); // path includes "/"

	setWindowModified(false);

	QString shownName = "Untitled";
	if (!fileName.isEmpty())
		shownName = strippedName(fileName);

	setWindowTitle(tr("%1[*] - %2").arg(shownName).arg(tr("PolyFit")));
}


bool MainWindow::open()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open file"), curDataDirectory_,
		tr("Supported Format (*.vg *.bvg *.obj)")
		);

	if (fileName.isEmpty())
		return false;

	return doOpen(fileName);
}


bool MainWindow::save()
{
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save file"), optimizedMeshFileName_,
		tr("Mesh (*.obj)")
		);

	if (fileName.isEmpty())
		return false;

	bool success = false;
	std::string ext = FileUtils::extension(fileName.toStdString());
	String::to_lowercase(ext);

	if (ext == "obj")
		success = MapIO::save(fileName.toStdString(), canvas()->optimizedMesh());
	else 
		success = PointSetIO::save(fileName.toStdString(), canvas()->pointSet());

	if (success) {
		setCurrentFile(fileName);
		status_message("File saved", 500);
		return true;
	}
	else {
		status_message("Saving failed", 500);
		return false;
	}
}


bool MainWindow::doOpen(const QString &fileName)
{
	std::string name = fileName.toStdString();
	std::string ext = FileUtils::extension(name);
	String::to_lowercase(ext);

	Map* mesh = nil;
	PointSet* pset = nil;
	if (ext == "obj") {
		mesh = MapIO::read(name);
		if (mesh)
			optimizedMeshFileName_ = fileName;
	}
	else {
		pset = PointSetIO::read(name);
		if (pset)
			pointCloudFileName_ = fileName;
	}

	if (mesh)
		canvas()->setMesh(mesh);

	if (pset) {
		canvas()->clear();
		canvas()->setPointSet(pset);

		hypothesisMeshFileName_ = pointCloudFileName_;
		int idx = fileName.lastIndexOf(".");
		hypothesisMeshFileName_.truncate(idx);
		hypothesisMeshFileName_.append(".eobj");

		optimizedMeshFileName_ = pointCloudFileName_;
		idx = fileName.lastIndexOf(".");
		optimizedMeshFileName_.truncate(idx);
		optimizedMeshFileName_.append(".obj");
	}

	if (pset || mesh) {
		updateStatusBar();
		setCurrentFile(fileName);
		setWindowTitle(tr("%1[*] - %2").arg(strippedName(fileName)).arg(tr("PolyFit")));
		status_message("File loaded", 500);

		if (pset) {
			checkBoxShowInput->setChecked(true);
			checkBoxShowCandidates->setChecked(true);
			checkBoxShowResult->setChecked(true);
			actionRefinePlanes->setDisabled(false);
			actionGenerateFacetHypothesis->setDisabled(true);
			actionGenerateQualityMeasures->setDisabled(true);
			actionOptimization->setDisabled(true);
		}
		return true;
	} 	
	else {	
		status_message("Open failed", 500);
		return false;
	} 
}


QString MainWindow::strippedName(const QString &fullFileName)
{
	return QFileInfo(fullFileName).fileName();
}


void MainWindow::snapshotScreen() {
	const std::string& fileName = optimizedMeshFileName_.toStdString();
	const std::string& snapshot = FileUtils::replace_extension(fileName, "png");

	QString snapshotFileName = QFileDialog::getSaveFileName(this,
		tr("Save Snapshot"), QString::fromStdString(snapshot),
		tr("PNG Image (*.png)")
	);

	canvas()->snapshotScreen(snapshotFileName);
}


LinearProgramSolver::SolverName MainWindow::active_solver() const {
	const QString& solverString = solverBox_->currentText();

	if (solverString == "GLPK")
		return LinearProgramSolver::GLPK;
#ifdef HAS_GUROBI
	else if (solverString == "GUROBI")
		return LinearProgramSolver::GUROBI;
#endif
	else if (solverString == "LPSOLVE")
		return LinearProgramSolver::LPSOLVE;
    else // default to SCIP
		return LinearProgramSolver::SCIP;
}
