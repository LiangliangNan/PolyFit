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


#ifndef WEIGHT_PANEL_CLICK_H
#define WEIGHT_PANEL_CLICK_H

#include <QDialog>
#include <ui_weight_panel_click.h>


class MainWindow;

class WeightPanelClick : public QDialog, public Ui::WeightPanelClickClass
{
	Q_OBJECT

public:
	WeightPanelClick(QWidget *parent);
	~WeightPanelClick() {}

	virtual void paintEvent(QPaintEvent *);

	virtual void mousePressEvent(QMouseEvent *);
	virtual void mouseMoveEvent(QMouseEvent *);

	void computeWeight();

public Q_SLOTS:
	void updateUI();

Q_SIGNALS:
	void weights_changed();

private:
	MainWindow*  mainWindow_;

	QPointF  pos_fitting_;
	QPointF  pos_coverage_;
	QPointF  pos_complexity_;

	float fitting_;
	float coverage_;
	float complexity_;

	QPointF  pos_;
	QPolygonF triangle_;
};

#endif // WEIGHT_PANEL_CLICK_H
