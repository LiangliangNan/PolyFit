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

#include "weight_panel_click.h"
#include <QPainter>
#include <QMouseEvent>
#include <QMessageBox>
#include "main_window.h"
#include "../math/math_types.h"
#include "../method/method_global.h"


static QPointF pointFromWeights(
	const QPointF& n, const QPointF& f, const QPointF& c,
	float wn, float wf, float wc)
{
	float x = (n.x() * wn + f.x() * wf + c.x() * wc) / (wn + wf + wc);
	float y = (n.y() * wn + f.y() * wf + c.y() * wc) / (wn + wf + wc);
	return QPointF(x, y);
}


static float area_of_triangle(QPointF a, QPointF b, QPointF c) {
	float area = Geom::triangle_area(vec2(a.x(), a.y()), vec2(b.x(), b.y()), vec2(c.x(), c.y()));
	return area;
}


WeightPanelClick::WeightPanelClick(QWidget *parent)
	: QDialog(parent)
	, triangle_(3)
{
	setupUi(this);

	mainWindow_ = dynamic_cast<MainWindow*>(parent);

	int x_offset = -20;
	pos_fitting_ = QPointF(160 + x_offset, 24);
	pos_coverage_ = QPointF(80 + x_offset, 180);
	pos_complexity_ = QPointF(250 + x_offset, 180);

	triangle_[0] = pos_fitting_;
	triangle_[1] = pos_coverage_;
	triangle_[2] = pos_complexity_;

	pos_ = pointFromWeights(pos_fitting_, pos_coverage_, pos_complexity_, Method::lambda_data_fitting, Method::lambda_model_coverage, Method::lambda_model_complexity);

	fitting_ = truncate_digits(Method::lambda_data_fitting, 3);
	coverage_ = truncate_digits(Method::lambda_model_coverage, 3);
	complexity_ = truncate_digits(Method::lambda_model_complexity, 3);

	setMouseTracking(true);
}


static QGradient gradient(const QColor &color, const QRectF &rect) {
	QColor c = color;
	c.setAlpha(160);
	QLinearGradient result(rect.topLeft(), rect.bottomRight());
	result.setColorAt(0, c.dark(150));
	result.setColorAt(0.5, c.light(200));
	result.setColorAt(1, c.dark(150));
	return result;
}


void WeightPanelClick::paintEvent(QPaintEvent* event) {
	QPainter painter(this);	
	painter.setRenderHint(QPainter::Antialiasing);

	QPalette pal = palette();
	QPen pen = pal.text().color();
	pen.setWidth(1);
	painter.setPen(pen);
	painter.setBrush(gradient(Qt::green, triangle_.boundingRect()));
	painter.drawPolygon(triangle_);

	//////////////////////////////////////////////////////////////////////////

	// paint the shape name
	painter.setBrush(pal.text());
	QString text_fitting = QString("Fitting (%1)").arg(fitting_);			painter.drawText(pos_fitting_.x() - 30, pos_fitting_.y() - 5, text_fitting);
	QString text_coverage = QString("Coverage (%1)").arg(coverage_);		painter.drawText(pos_coverage_.x() - 40, pos_coverage_.y() + 20, text_coverage);
	QString text_complexity = QString("Complexity (%1)").arg(complexity_);	painter.drawText(pos_complexity_.x() - 60, pos_complexity_.y() + 20, text_complexity);

	painter.setPen(Qt::red);
	painter.setBrush(Qt::red);
	painter.drawEllipse(pos_, 5, 5);

	painter.setRenderHint(QPainter::Antialiasing, false);
	painter.end();
}


void WeightPanelClick::mousePressEvent(QMouseEvent* event)
{
	QPointF p = event->pos();

	if (triangle_.containsPoint(p, Qt::OddEvenFill)) {
		pos_ = p;
		computeWeight();
		update();
	}
	else {
		//QMessageBox::warning(this, "Warning", "Click outside of triangle doesn't make sense");
	}

	QWidget::mousePressEvent(event);
}


void WeightPanelClick::mouseMoveEvent(QMouseEvent* event)
{
	QPointF p = event->pos();

	if (triangle_.containsPoint(p, Qt::OddEvenFill)) {
		float area_sum = area_of_triangle(pos_fitting_, pos_coverage_, pos_complexity_);

		float area_data_fitting = area_of_triangle(pos_complexity_, pos_coverage_, p);
		float area_model_coverage = area_of_triangle(pos_complexity_, pos_fitting_, p);
		float area_model_complexity = area_of_triangle(pos_fitting_, pos_coverage_, p);

		float fitting = truncate_digits(area_data_fitting / area_sum, 3);
		float coverage = truncate_digits(area_model_coverage / area_sum, 3);
		float complexity = truncate_digits(area_model_complexity / area_sum, 3);
		std::string str_fitting = std::to_string(fitting);		
		std::string str_coverage = std::to_string(coverage);	
		std::string str_complexity = std::to_string(complexity);

		std::string weights = 
			"Fitting: " + str_fitting.substr(0, str_fitting.find(".") + 4) +
			"; Coverage: " + str_coverage.substr(0, str_coverage.find(".") + 4) +
			"; Complexity: " + str_complexity.substr(0, str_complexity.find(".") + 4);
		mainWindow_->status_message(weights, 2000);
	}

	QWidget::mousePressEvent(event);
}


void WeightPanelClick::computeWeight() {
	float area_sum = area_of_triangle(pos_fitting_, pos_coverage_, pos_complexity_);

	float area_data_fitting = area_of_triangle(pos_complexity_, pos_coverage_, pos_);
	float area_model_coverage = area_of_triangle(pos_complexity_, pos_fitting_, pos_);
	float area_model_complexity = area_of_triangle(pos_fitting_, pos_coverage_, pos_);

	Method::lambda_data_fitting = area_data_fitting / area_sum;
	Method::lambda_model_coverage = area_model_coverage / area_sum;
	Method::lambda_model_complexity = area_model_complexity / area_sum;

	fitting_ = truncate_digits(Method::lambda_data_fitting, 3);
	coverage_ = truncate_digits(Method::lambda_model_coverage, 3);
	complexity_ = truncate_digits(Method::lambda_model_complexity, 3);

	emit weights_changed();

	update();
}


void WeightPanelClick::updateUI() {
	float sum = Method::lambda_data_fitting + Method::lambda_model_coverage + Method::lambda_model_complexity;
	fitting_ = truncate_digits(Method::lambda_data_fitting / sum, 3);
	coverage_ = truncate_digits(Method::lambda_model_coverage / sum, 3);
	complexity_ = truncate_digits(Method::lambda_model_complexity / sum, 3);

	pos_ = pointFromWeights(pos_fitting_, pos_coverage_, pos_complexity_, Method::lambda_data_fitting, Method::lambda_model_coverage, Method::lambda_model_complexity);
	update();
}
