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


#ifndef WEIGHT_PANEL_MANUAL_H
#define WEIGHT_PANEL_MANUAL_H

#include <QWidget>
#include "ui_weight_panel_manual.h"


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
