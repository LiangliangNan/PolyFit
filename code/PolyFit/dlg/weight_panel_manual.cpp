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

#include "weight_panel_manual.h"

#include <math/math_types.h>
#include <method/method_global.h>

#include "main_window.h"
#include "paint_canvas.h"

#include "ui_weight_panel_manual.h"


WeightPanelManual::WeightPanelManual(QWidget *parent)
	: QDialog(parent)
{
	setupUi(this);
	mainWindow_ = dynamic_cast<MainWindow*>(parent);
	updateUI();
}


void WeightPanelManual::updateUI() {
	float fitting = truncate_digits(Method::weight_data_fitting, 3);
	float coverage = truncate_digits(Method::weight_model_coverage, 3);
	float complexity = truncate_digits(Method::weight_model_complexity, 3);

	QString text_fitting = QString("%1").arg(fitting);
	QString text_coverage = QString("%1").arg(coverage);
	QString text_complexity = QString("%1").arg(complexity);

	lineEditFitting->setText(text_fitting);
	lineEditCoverage->setText(text_coverage);
	lineEditComplexity->setText(text_complexity);
}


void WeightPanelManual::updateWeights() {
	Method::weight_data_fitting = lineEditFitting->text().toFloat();
	Method::weight_model_coverage = lineEditCoverage->text().toFloat();
	Method::weight_model_complexity = lineEditComplexity->text().toFloat();
}