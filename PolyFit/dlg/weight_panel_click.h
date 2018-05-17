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


#ifndef WEIGHT_PANEL_CLICK_H
#define WEIGHT_PANEL_CLICK_H

#include <QWidget>
#include "ui_weight_panel_click.h"
#include <QGradient>


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
