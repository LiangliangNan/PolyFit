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


#ifndef WEIGHT_PANEL_MANUAL_H
#define WEIGHT_PANEL_MANUAL_H

#include <QDialog>
#include <ui_weight_panel_manual.h>

class MainWindow;

class WeightPanelManual : public QDialog, public Ui::WeightPanelManualClass
{
	Q_OBJECT

public:
	WeightPanelManual(QWidget *parent);
	~WeightPanelManual() {}

	void updateWeights();

public Q_SLOTS:
	void updateUI();

private:
	MainWindow*  mainWindow_;
};

#endif // WEIGHT_PANEL_MANUAL_H
