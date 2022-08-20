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

#ifndef QGLVIEWER_MANIPULATED_CAMERA_FRAME_H
#define QGLVIEWER_MANIPULATED_CAMERA_FRAME_H

#include "manipulatedFrame.h"

namespace qglviewer {
/*! \brief The ManipulatedCameraFrame class represents a ManipulatedFrame with Camera specific mouse bindings.
  \class ManipulatedCameraFrame manipulatedCameraFrame.h QGLViewer/manipulatedCameraFrame.h

  A ManipulatedCameraFrame is a specialization of a ManipulatedFrame, designed to be set as the
  Camera::frame(). Mouse motions are basically interpreted in a negated way: when the mouse goes to
  the right, the ManipulatedFrame translation goes to the right, while the ManipulatedCameraFrame
  has to go to the \e left, so that the \e scene seems to move to the right.

  A ManipulatedCameraFrame rotates around its pivotPoint(), which corresponds to the
  associated Camera::pivotPoint().

  A ManipulatedCameraFrame can also "fly" in the scene. It basically moves forward, and turns
  according to the mouse motion. See flySpeed(), sceneUpVector() and the QGLViewer::MOVE_FORWARD and
  QGLViewer::MOVE_BACKWARD QGLViewer::MouseAction.

  See the <a href="../mouse.html">mouse page</a> for a description of the possible actions that can
  be performed using the mouse and their bindings.
  \nosubgrouping */
class QGLVIEWER_EXPORT ManipulatedCameraFrame : public qglviewer::ManipulatedFrame
{
#ifndef DOXYGEN
	friend class Camera;
	friend class ::QGLViewer;
#endif

	Q_OBJECT

public:
	ManipulatedCameraFrame();
	/*! Virtual destructor. Empty. */
	virtual ~ManipulatedCameraFrame() {}

	ManipulatedCameraFrame(const ManipulatedCameraFrame& mcf);
	ManipulatedCameraFrame& operator=(const ManipulatedCameraFrame& mcf);

	/*! @name Pivot point */
	//@{
public:
	/*! Returns the point the ManipulatedCameraFrame pivot point, around which the camera rotates.

	It is defined in the world coordinate system. Default value is (0,0,0).

	When the ManipulatedCameraFrame is associated to a Camera, Camera::pivotPoint() also
	returns this value. This point can interactively be changed using the mouse (see
	Camera::setPivotPointFromPixel() and QGLViewer::RAP_FROM_PIXEL and QGLViewer::RAP_IS_CENTER
	in the <a href="../mouse.html">mouse page</a>). */
	Vec pivotPoint() const { return pivotPoint_; }
	/*! Sets the pivotPoint(), defined in the world coordinate system. */
	void setPivotPoint(const Vec& point) { pivotPoint_ = point; }

#ifndef DOXYGEN
	Vec revolveAroundPoint() const { qWarning("revolveAroundPoint() is deprecated, use pivotPoint() instead"); return  pivotPoint(); }
	void setRevolveArountPoint(const Vec& point) { qWarning("setRevolveAroundPoint() is deprecated, use setPivotPoint() instead"); setPivotPoint(point); }
#endif
	//@}

	/*! @name Camera manipulation */
	//@{
public:
	/*! Returns \c true when the frame's rotation is constrained around the sceneUpVector(),
		and \c false otherwise, when the rotation is completely free (default).

		In free mode, the associated camera can be arbitrarily rotated in the scene, along its
		three axis, thus possibly leading to any arbitrary orientation.

		When you setRotatesAroundUpVector() to \c true, the sceneUpVector() defines a
		'vertical' direction around which the camera rotates. The camera can rotate left
		or right, around this axis. It can also be moved up or down to show the 'top' and
		'bottom' views of the scene. As a result, the sceneUpVector() will always appear vertical
		in the scene, and the horizon is preserved and stays projected along the camera's
		horizontal axis.

		Note that setting this value to \c true when the sceneUpVector() is not already
		vertically projected will break these invariants. It will also limit the possible movement
		of the camera, possibly up to a lock when the sceneUpVector() is projected horizontally.
		Use Camera::setUpVector() to define the sceneUpVector() and align the camera before calling
		this method to ensure this does not happen. */
	bool rotatesAroundUpVector() const { return rotatesAroundUpVector_; }
	/*! Sets the value of rotatesAroundUpVector().

	   Default value is false (free rotation). */
	void setRotatesAroundUpVector(bool constrained) { rotatesAroundUpVector_ = constrained; }

	/*! Returns whether or not the QGLViewer::ZOOM action zooms on the pivot point.

	  When set to \c false (default), a zoom action will move the camera along its Camera::viewDirection(),
	  i.e. back and forth along a direction perpendicular to the projection screen.

	  setZoomsOnPivotPoint() to \c true will move the camera along an axis defined by the
	  Camera::pivotPoint() and its current position instead. As a result, the projected position of the
	  pivot point on screen will stay the same during a zoom. */
	bool zoomsOnPivotPoint() const { return zoomsOnPivotPoint_; }
	/*! Sets the value of zoomsOnPivotPoint().

	   Default value is false. */
	void setZoomsOnPivotPoint(bool enabled) { zoomsOnPivotPoint_ = enabled; }

private:
#ifndef DOXYGEN
	void zoom(qreal delta, const Camera * const camera);
#endif
	//@}

	/*! @name Fly parameters */
	//@{
public Q_SLOTS:
	/*! Sets the flySpeed(), defined in OpenGL units.

	Default value is 0.0, but it is modified according to the QGLViewer::sceneRadius() when the
	ManipulatedCameraFrame is set as the Camera::frame(). */
	void setFlySpeed(qreal speed) { flySpeed_ = speed; }

	/*! Sets the sceneUpVector(), defined in the world coordinate system.

	Default value is (0,1,0), but it is updated by the Camera when this object is set as its Camera::frame().
	Using Camera::setUpVector() instead is probably a better solution. */
	void setSceneUpVector(const Vec& up) { sceneUpVector_ = up; }

public:
	/*! Returns the fly speed, expressed in OpenGL units.

	It corresponds to the incremental displacement that is periodically applied to the
	ManipulatedCameraFrame position when a QGLViewer::MOVE_FORWARD or QGLViewer::MOVE_BACKWARD
	QGLViewer::MouseAction is proceeded.

	\attention When the ManipulatedCameraFrame is set as the Camera::frame(), this value is set
	according to the QGLViewer::sceneRadius() by QGLViewer::setSceneRadius(). */
	qreal flySpeed() const { return flySpeed_; }

	/*! Returns the up vector of the scene, expressed in the world coordinate system.

	In 'fly mode' (corresponding to the QGLViewer::MOVE_FORWARD and QGLViewer::MOVE_BACKWARD
	QGLViewer::MouseAction bindings), horizontal displacements of the mouse rotate
	the ManipulatedCameraFrame around this vector. Vertical displacements rotate always around the
	Camera \c X axis.

	This value is also used when setRotationIsConstrained() is set to \c true to define the up vector
	(and incidentally the 'horizon' plane) around which the camera will rotate.

	Default value is (0,1,0), but it is updated by the Camera when this object is set as its Camera::frame().
	Camera::setOrientation() and Camera::setUpVector()) direclty modify this value and should be used
	instead. */
	Vec sceneUpVector() const { return sceneUpVector_; }

#ifndef DOXYGEN
	Vec flyUpVector() const;
	void setFlyUpVector(const Vec& up);
#endif
	//@}

	/*! @name Mouse event handlers */
	//@{
protected:
	virtual void mouseReleaseEvent(QMouseEvent* const event, Camera* const camera);
	virtual void mouseMoveEvent   (QMouseEvent* const event, Camera* const camera);
	virtual void wheelEvent       (QWheelEvent* const event, Camera* const camera);
	//@}

	/*! @name Spinning */
	//@{
protected Q_SLOTS:
	virtual void spin();
	//@}

	/*! @name XML representation */
	//@{
public:
	virtual QDomElement domElement(const QString& name, QDomDocument& document) const;
public Q_SLOTS:
	virtual void initFromDOMElement(const QDomElement& element);
	//@}

#ifndef DOXYGEN
protected:
	virtual void startAction(int ma, bool withConstraint=true); // int is really a QGLViewer::MouseAction
#endif

private Q_SLOTS:
	virtual void flyUpdate();

private:
	void updateSceneUpVector();
	Quaternion turnQuaternion(int x, const Camera* const camera);
	Quaternion pitchYawQuaternion(int x, int y, const Camera* const camera);

private:
	// Fly mode data
	qreal flySpeed_;
	qreal driveSpeed_;
	Vec sceneUpVector_;
	QTimer flyTimer_;

	bool rotatesAroundUpVector_;
	// Inverse the direction of an horizontal mouse motion. Depends on the projected
	// screen orientation of the vertical axis when the mouse button is pressed.
	bool constrainedRotationIsReversed_;

	bool zoomsOnPivotPoint_;

	Vec pivotPoint_;
};

} // namespace qglviewer

#endif // QGLVIEWER_MANIPULATED_CAMERA_FRAME_H
