/****************************************************************************

 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of the QGLViewer library version 2.6.3.

 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 versions 2.0 or 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.
 In addition, as a special exception, Gilles Debunne gives you certain 
 additional rights, described in the file GPL_EXCEPTION in this package.

 libQGLViewer uses dual licensing. Commercial/proprietary software must
 purchase a libQGLViewer Commercial License.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/

#include "domUtils.h"
#include "manipulatedCameraFrame.h"
#include "qglviewer.h"

#include <QMouseEvent>

using namespace qglviewer;
using namespace std;

/*! Default constructor.

 flySpeed() is set to 0.0 and sceneUpVector() is (0,1,0). The pivotPoint() is set to (0,0,0).

  \attention Created object is removeFromMouseGrabberPool(). */
ManipulatedCameraFrame::ManipulatedCameraFrame()
    : driveSpeed_(0.0), sceneUpVector_(0.0, 1.0, 0.0), rotatesAroundUpVector_(false), zoomsOnPivotPoint_(true)
{
	setFlySpeed(0.0);
	removeFromMouseGrabberPool();
	connect(&flyTimer_, SIGNAL(timeout()), SLOT(flyUpdate()));
}

/*! Equal operator. Calls ManipulatedFrame::operator=() and then copy attributes. */
ManipulatedCameraFrame& ManipulatedCameraFrame::operator=(const ManipulatedCameraFrame& mcf)
{
	ManipulatedFrame::operator=(mcf);

	setFlySpeed(mcf.flySpeed());
	setSceneUpVector(mcf.sceneUpVector());
	setRotatesAroundUpVector(mcf.rotatesAroundUpVector_);
	setZoomsOnPivotPoint(mcf.zoomsOnPivotPoint_);

	return *this;
}

/*! Copy constructor. Performs a deep copy of all members using operator=(). */
ManipulatedCameraFrame::ManipulatedCameraFrame(const ManipulatedCameraFrame& mcf)
	: ManipulatedFrame(mcf)
{
	removeFromMouseGrabberPool();
	connect(&flyTimer_, SIGNAL(timeout()), SLOT(flyUpdate()));
	(*this)=(mcf);
}

////////////////////////////////////////////////////////////////////////////////

/*! Overloading of ManipulatedFrame::spin().

Rotates the ManipulatedCameraFrame around its pivotPoint() instead of its origin. */
void ManipulatedCameraFrame::spin()
{
	rotateAroundPoint(spinningQuaternion(), pivotPoint());
}

#ifndef DOXYGEN
/*! Called for continuous frame motion in fly mode (see QGLViewer::MOVE_FORWARD). Emits
  manipulated(). */
void ManipulatedCameraFrame::flyUpdate()
{
	static Vec flyDisp(0.0, 0.0, 0.0);
	switch (action_)
	{
		case QGLViewer::MOVE_FORWARD:
			flyDisp.z = -flySpeed();
			translate(localInverseTransformOf(flyDisp));
			break;
		case QGLViewer::MOVE_BACKWARD:
			flyDisp.z = flySpeed();
			translate(localInverseTransformOf(flyDisp));
			break;
		case QGLViewer::DRIVE:
			flyDisp.z = flySpeed() * driveSpeed_;
			translate(localInverseTransformOf(flyDisp));
			break;
		default:
			break;
	}

	// Needs to be out of the switch since ZOOM/fastDraw()/wheelEvent use this callback to trigger a final draw().
	// #CONNECTION# wheelEvent.
	Q_EMIT manipulated();
}

Vec ManipulatedCameraFrame::flyUpVector() const {
	qWarning("flyUpVector() is deprecated. Use sceneUpVector() instead.");
	return sceneUpVector();
}

void ManipulatedCameraFrame::setFlyUpVector(const Vec& up) {
	qWarning("setFlyUpVector() is deprecated. Use setSceneUpVector() instead.");
	setSceneUpVector(up);
}

#endif

/*! This method will be called by the Camera when its orientation is changed, so that the
sceneUpVector (private) is changed accordingly. You should not need to call this method. */
void ManipulatedCameraFrame::updateSceneUpVector()
{
	sceneUpVector_ = inverseTransformOf(Vec(0.0, 1.0, 0.0));
}

////////////////////////////////////////////////////////////////////////////////
//          S t a t e   s a v i n g   a n d   r e s t o r i n g               //
////////////////////////////////////////////////////////////////////////////////

/*! Returns an XML \c QDomElement that represents the ManipulatedCameraFrame.

 Adds to the ManipulatedFrame::domElement() the ManipulatedCameraFrame specific informations in a \c
 ManipulatedCameraParameters child QDomElement.

 \p name is the name of the QDomElement tag. \p doc is the \c QDomDocument factory used to create
 QDomElement.

 Use initFromDOMElement() to restore the ManipulatedCameraFrame state from the resulting
 \c QDomElement.

 See Vec::domElement() for a complete example. See also Quaternion::domElement(),
 Frame::domElement(), Camera::domElement()... */
QDomElement ManipulatedCameraFrame::domElement(const QString& name, QDomDocument& document) const
{
	QDomElement e = ManipulatedFrame::domElement(name, document);
	QDomElement mcp = document.createElement("ManipulatedCameraParameters");
	mcp.setAttribute("flySpeed", QString::number(flySpeed()));
	DomUtils::setBoolAttribute(mcp, "rotatesAroundUpVector", rotatesAroundUpVector());
	DomUtils::setBoolAttribute(mcp, "zoomsOnPivotPoint", zoomsOnPivotPoint());
	mcp.appendChild(sceneUpVector().domElement("sceneUpVector", document));
	e.appendChild(mcp);
	return e;
}

/*! Restores the ManipulatedCameraFrame state from a \c QDomElement created by domElement().

First calls ManipulatedFrame::initFromDOMElement() and then initializes ManipulatedCameraFrame
specific parameters. */
void ManipulatedCameraFrame::initFromDOMElement(const QDomElement& element)
{
	// No need to initialize, since default sceneUpVector and flySpeed are not meaningful.
	// It's better to keep current ones. And it would destroy constraint() and referenceFrame().
	// *this = ManipulatedCameraFrame();
	ManipulatedFrame::initFromDOMElement(element);

	QDomElement child=element.firstChild().toElement();
	while (!child.isNull())
	{
		if (child.tagName() == "ManipulatedCameraParameters")
		{
			setFlySpeed(DomUtils::qrealFromDom(child, "flySpeed", flySpeed()));
			setRotatesAroundUpVector(DomUtils::boolFromDom(child, "rotatesAroundUpVector", false));
			setZoomsOnPivotPoint(DomUtils::boolFromDom(child, "zoomsOnPivotPoint", false));

			QDomElement schild=child.firstChild().toElement();
			while (!schild.isNull())
			{
				if (schild.tagName() == "sceneUpVector")
					setSceneUpVector(Vec(schild));

				schild = schild.nextSibling().toElement();
			}
		}
		child = child.nextSibling().toElement();
	}
}


////////////////////////////////////////////////////////////////////////////////
//                 M o u s e    h a n d l i n g                               //
////////////////////////////////////////////////////////////////////////////////

#ifndef DOXYGEN
/*! Protected internal method used to handle mouse events. */
void ManipulatedCameraFrame::startAction(int ma, bool withConstraint)
{
	ManipulatedFrame::startAction(ma, withConstraint);

	switch (action_)
	{
		case QGLViewer::MOVE_FORWARD:
		case QGLViewer::MOVE_BACKWARD:
		case QGLViewer::DRIVE:
			flyTimer_.setSingleShot(false);
			flyTimer_.start(10);
			break;
		case QGLViewer::ROTATE:
			constrainedRotationIsReversed_ = transformOf(sceneUpVector_).y < 0.0;
			break;
		default:
			break;
	}
}

void ManipulatedCameraFrame::zoom(qreal delta, const Camera * const camera) {
	const qreal sceneRadius = camera->sceneRadius();
	if (zoomsOnPivotPoint_) {
		Vec direction = position() - camera->pivotPoint();
		if (direction.norm() > 0.02 * sceneRadius || delta > 0.0)
			translate(delta * direction);
	} else {
		const qreal coef = qMax(fabs((camera->frame()->coordinatesOf(camera->pivotPoint())).z), 0.2 * sceneRadius);
		Vec trans(0.0, 0.0, -coef * delta);
		translate(inverseTransformOf(trans));
	}
}

#endif

/*! Overloading of ManipulatedFrame::mouseMoveEvent().

Motion depends on mouse binding (see <a href="../mouse.html">mouse page</a> for details). The
resulting displacements are basically inverted from those of a ManipulatedFrame. */
void ManipulatedCameraFrame::mouseMoveEvent(QMouseEvent* const event, Camera* const camera)
{
	// #CONNECTION# QGLViewer::mouseMoveEvent does the update().
	switch (action_)
	{
		case QGLViewer::TRANSLATE:
		{
			const QPoint delta = prevPos_ - event->pos();
			Vec trans(delta.x(), -delta.y(), 0.0);
			// Scale to fit the screen mouse displacement
			switch (camera->type())
			{
				case Camera::PERSPECTIVE :
					trans *= 2.0 * tan(camera->fieldOfView()/2.0) *
							fabs((camera->frame()->coordinatesOf(pivotPoint())).z) / camera->screenHeight();
					break;
				case Camera::ORTHOGRAPHIC :
				{
					GLdouble w,h;
					camera->getOrthoWidthHeight(w, h);
					trans[0] *= 2.0 * w / camera->screenWidth();
					trans[1] *= 2.0 * h / camera->screenHeight();
					break;
				}
			}
			translate(inverseTransformOf(translationSensitivity()*trans));
			break;
		}

		case QGLViewer::MOVE_FORWARD:
		{
			Quaternion rot = pitchYawQuaternion(event->x(), event->y(), camera);
			rotate(rot);
			//#CONNECTION# wheelEvent MOVE_FORWARD case
			// actual translation is made in flyUpdate().
			//translate(inverseTransformOf(Vec(0.0, 0.0, -flySpeed())));
			break;
		}

		case QGLViewer::MOVE_BACKWARD:
		{
			Quaternion rot = pitchYawQuaternion(event->x(), event->y(), camera);
			rotate(rot);
			// actual translation is made in flyUpdate().
			//translate(inverseTransformOf(Vec(0.0, 0.0, flySpeed())));
			break;
		}

		case QGLViewer::DRIVE:
		{
			Quaternion rot = turnQuaternion(event->x(), camera);
			rotate(rot);
			// actual translation is made in flyUpdate().
			driveSpeed_ = 0.01 * (event->y() - pressPos_.y());
			break;
		}

		case QGLViewer::ZOOM:
		{
			zoom(deltaWithPrevPos(event, camera), camera);
			break;
		}

		case QGLViewer::LOOK_AROUND:
		{
			Quaternion rot = pitchYawQuaternion(event->x(), event->y(), camera);
			rotate(rot);
			break;
		}

		case QGLViewer::ROTATE:
		{
			Quaternion rot;
			if (rotatesAroundUpVector_) {
				// Multiply by 2.0 to get on average about the same speed as with the deformed ball
				qreal dx = 2.0 * rotationSensitivity() * (prevPos_.x() - event->x()) / camera->screenWidth();
				qreal dy = 2.0 * rotationSensitivity() * (prevPos_.y() - event->y()) / camera->screenHeight();
				if (constrainedRotationIsReversed_) dx = -dx;
				Vec verticalAxis = transformOf(sceneUpVector_);
				rot = Quaternion(verticalAxis, dx) * Quaternion(Vec(1.0, 0.0, 0.0), dy);
			} else {
				Vec trans = camera->projectedCoordinatesOf(pivotPoint());
				rot = deformedBallQuaternion(event->x(), event->y(), trans[0], trans[1], camera);
			}
			//#CONNECTION# These two methods should go together (spinning detection and activation)
			computeMouseSpeed(event);
			setSpinningQuaternion(rot);
			spin();
			break;
		}

		case QGLViewer::SCREEN_ROTATE:
		{
			Vec trans = camera->projectedCoordinatesOf(pivotPoint());

			const qreal angle = atan2(event->y() - trans[1], event->x() - trans[0]) - atan2(prevPos_.y()-trans[1], prevPos_.x()-trans[0]);

			Quaternion rot(Vec(0.0, 0.0, 1.0), angle);
			//#CONNECTION# These two methods should go together (spinning detection and activation)
			computeMouseSpeed(event);
			setSpinningQuaternion(rot);
			spin();
			updateSceneUpVector();
			break;
		}

		case QGLViewer::ROLL:
		{
			const qreal angle = M_PI * (event->x() - prevPos_.x()) / camera->screenWidth();
			Quaternion rot(Vec(0.0, 0.0, 1.0), angle);
			rotate(rot);
			setSpinningQuaternion(rot);
			updateSceneUpVector();
			break;
		}

		case QGLViewer::SCREEN_TRANSLATE:
		{
			Vec trans;
			int dir = mouseOriginalDirection(event);
			if (dir == 1)
				trans.setValue(prevPos_.x() - event->x(), 0.0, 0.0);
			else if (dir == -1)
				trans.setValue(0.0, event->y() - prevPos_.y(), 0.0);

			switch (camera->type())
			{
				case Camera::PERSPECTIVE :
					trans *= 2.0 * tan(camera->fieldOfView()/2.0) *
							fabs((camera->frame()->coordinatesOf(pivotPoint())).z) / camera->screenHeight();
					break;
				case Camera::ORTHOGRAPHIC :
				{
					GLdouble w,h;
					camera->getOrthoWidthHeight(w, h);
					trans[0] *= 2.0 * w / camera->screenWidth();
					trans[1] *= 2.0 * h / camera->screenHeight();
					break;
				}
			}

			translate(inverseTransformOf(translationSensitivity()*trans));
			break;
		}

		case QGLViewer::ZOOM_ON_REGION:
		case QGLViewer::NO_MOUSE_ACTION:
			break;
	}

	if (action_ != QGLViewer::NO_MOUSE_ACTION)
	{
		prevPos_ = event->pos();
		if (action_ != QGLViewer::ZOOM_ON_REGION)
			// ZOOM_ON_REGION should not emit manipulated().
			// prevPos_ is used to draw rectangle feedback.
			Q_EMIT manipulated();
	}
}


/*! This is an overload of ManipulatedFrame::mouseReleaseEvent(). The QGLViewer::MouseAction is
  terminated. */
void ManipulatedCameraFrame::mouseReleaseEvent(QMouseEvent* const event, Camera* const camera)
{
	if ((action_ == QGLViewer::MOVE_FORWARD) || (action_ == QGLViewer::MOVE_BACKWARD) || (action_ == QGLViewer::DRIVE))
		flyTimer_.stop();

	if (action_ == QGLViewer::ZOOM_ON_REGION)
		camera->fitScreenRegion(QRect(pressPos_, event->pos()));

	ManipulatedFrame::mouseReleaseEvent(event, camera);
}

/*! This is an overload of ManipulatedFrame::wheelEvent().

The wheel behavior depends on the wheel binded action. Current possible actions are QGLViewer::ZOOM,
QGLViewer::MOVE_FORWARD, QGLViewer::MOVE_BACKWARD. QGLViewer::ZOOM speed depends on
wheelSensitivity() while QGLViewer::MOVE_FORWARD and QGLViewer::MOVE_BACKWARD depend on flySpeed().
See QGLViewer::setWheelBinding() to customize the binding. */
void ManipulatedCameraFrame::wheelEvent(QWheelEvent* const event, Camera* const camera)
{
	//#CONNECTION# QGLViewer::setWheelBinding, ManipulatedFrame::wheelEvent.
	switch (action_)
	{
		case QGLViewer::ZOOM:
		{
			zoom(wheelDelta(event), camera);
			Q_EMIT manipulated();
			break;
		}
		case QGLViewer::MOVE_FORWARD:
		case QGLViewer::MOVE_BACKWARD:
			//#CONNECTION# mouseMoveEvent() MOVE_FORWARD case
			translate(inverseTransformOf(Vec(0.0, 0.0, 0.2*flySpeed()*event->delta())));
			Q_EMIT manipulated();
			break;
		default:
			break;
	}

	// #CONNECTION# startAction should always be called before
	if (previousConstraint_)
		setConstraint(previousConstraint_);

	// The wheel triggers a fastDraw. A final update() is needed after the last wheel event to
	// polish the rendering using draw(). Since the last wheel event does not say its name, we use
	// the flyTimer_ to trigger flyUpdate(), which emits manipulated. Two wheel events
	// separated by more than this delay millisec will trigger a draw().
	const int finalDrawAfterWheelEventDelay = 400;

	// Starts (or prolungates) the timer.
	flyTimer_.setSingleShot(true);
	flyTimer_.start(finalDrawAfterWheelEventDelay);

	// This could also be done *before* manipulated is emitted, so that isManipulated() returns false.
	// But then fastDraw would not be used with wheel.
	// Detecting the last wheel event and forcing a final draw() is done using the timer_.
	action_ = QGLViewer::NO_MOUSE_ACTION;
}

////////////////////////////////////////////////////////////////////////////////

/*! Returns a Quaternion that is a rotation around current camera Y, proportionnal to the horizontal mouse position. */
Quaternion ManipulatedCameraFrame::turnQuaternion(int x, const Camera* const camera)
{
	return Quaternion(Vec(0.0, 1.0, 0.0), rotationSensitivity()*(prevPos_.x()-x)/camera->screenWidth());
}

/*! Returns a Quaternion that is the composition of two rotations, inferred from the
  mouse pitch (X axis) and yaw (sceneUpVector() axis). */
Quaternion ManipulatedCameraFrame::pitchYawQuaternion(int x, int y, const Camera* const camera)
{
	const Quaternion rotX(Vec(1.0, 0.0, 0.0), rotationSensitivity()*(prevPos_.y()-y)/camera->screenHeight());
	const Quaternion rotY(transformOf(sceneUpVector()), rotationSensitivity()*(prevPos_.x()-x)/camera->screenWidth());
	return rotY * rotX;
}
