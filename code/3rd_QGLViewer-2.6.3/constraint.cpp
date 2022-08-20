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

#include "constraint.h"
#include "frame.h"
#include "camera.h"
#include "manipulatedCameraFrame.h"

using namespace qglviewer;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
//                                  Constraint                                //
////////////////////////////////////////////////////////////////////////////////

/*! Default constructor.

translationConstraintType() and rotationConstraintType() are set to AxisPlaneConstraint::FREE.
translationConstraintDirection() and rotationConstraintDirection() are set to (0,0,0). */
AxisPlaneConstraint::AxisPlaneConstraint()
	: translationConstraintType_(FREE), rotationConstraintType_(FREE)
{
	// Do not use set since setRotationConstraintType needs a read.
}

/*! Simply calls setTranslationConstraintType() and setTranslationConstraintDirection(). */
void AxisPlaneConstraint::setTranslationConstraint(Type type, const Vec& direction)
{
	setTranslationConstraintType(type);
	setTranslationConstraintDirection(direction);
}

/*! Defines the translationConstraintDirection(). The coordinate system where \p direction is expressed depends on your class implementation. */
void AxisPlaneConstraint::setTranslationConstraintDirection(const Vec& direction)
{
	if ((translationConstraintType()!=AxisPlaneConstraint::FREE) && (translationConstraintType()!=AxisPlaneConstraint::FORBIDDEN))
	{
		const qreal norm = direction.norm();
		if (norm < 1E-8)
		{
			qWarning("AxisPlaneConstraint::setTranslationConstraintDir: null vector for translation constraint");
			translationConstraintType_ = AxisPlaneConstraint::FREE;
		}
		else
			translationConstraintDir_ = direction/norm;
	}
}

/*! Simply calls setRotationConstraintType() and setRotationConstraintDirection(). */
void AxisPlaneConstraint::setRotationConstraint(Type type, const Vec& direction)
{
	setRotationConstraintType(type);
	setRotationConstraintDirection(direction);
}

/*! Defines the rotationConstraintDirection(). The coordinate system where \p direction is expressed depends on your class implementation. */
void AxisPlaneConstraint::setRotationConstraintDirection(const Vec& direction)
{
	if ((rotationConstraintType()!=AxisPlaneConstraint::FREE) && (rotationConstraintType()!=AxisPlaneConstraint::FORBIDDEN))
	{
		const qreal norm = direction.norm();
		if (norm < 1E-8)
		{
			qWarning("AxisPlaneConstraint::setRotationConstraintDir: null vector for rotation constraint");
			rotationConstraintType_ = AxisPlaneConstraint::FREE;
		}
		else
			rotationConstraintDir_ = direction/norm;
	}
}

/*! Set the Type() of the rotationConstraintType(). Default is AxisPlaneConstraint::FREE.

 Depending on this value, the Frame will freely rotate (AxisPlaneConstraint::FREE), will only be able
 to rotate around an axis (AxisPlaneConstraint::AXIS), or will not able to rotate at all
 (AxisPlaneConstraint::FORBIDDEN).

 Use Frame::setOrientation() to define the orientation of the constrained Frame before it gets
 constrained.

 \attention An AxisPlaneConstraint::PLANE Type() is not meaningful for rotational constraints and
 will be ignored. */
void AxisPlaneConstraint::setRotationConstraintType(Type type)
{
	if (rotationConstraintType() == AxisPlaneConstraint::PLANE)
	{
		qWarning("AxisPlaneConstraint::setRotationConstraintType: the PLANE type cannot be used for a rotation constraints");
		return;
	}

	rotationConstraintType_ = type;
}


////////////////////////////////////////////////////////////////////////////////
//                               LocalConstraint                              //
////////////////////////////////////////////////////////////////////////////////

/*! Depending on translationConstraintType(), constrain \p translation to be along an axis or
  limited to a plane defined in the Frame local coordinate system by
  translationConstraintDirection(). */
void LocalConstraint::constrainTranslation(Vec& translation, Frame* const frame)
{
	Vec proj;
	switch (translationConstraintType())
	{
	case AxisPlaneConstraint::FREE:
		break;
	case AxisPlaneConstraint::PLANE:
		proj = frame->rotation().rotate(translationConstraintDirection());
		translation.projectOnPlane(proj);
		break;
	case AxisPlaneConstraint::AXIS:
		proj = frame->rotation().rotate(translationConstraintDirection());
		translation.projectOnAxis(proj);
		break;
	case AxisPlaneConstraint::FORBIDDEN:
		translation = Vec(0.0, 0.0, 0.0);
		break;
	}
}

/*! When rotationConstraintType() is AxisPlaneConstraint::AXIS, constrain \p rotation to be a rotation
  around an axis whose direction is defined in the Frame local coordinate system by
  rotationConstraintDirection(). */
void LocalConstraint::constrainRotation(Quaternion& rotation, Frame* const)
{
	switch (rotationConstraintType())
	{
	case AxisPlaneConstraint::FREE:
		break;
	case AxisPlaneConstraint::PLANE:
		break;
	case AxisPlaneConstraint::AXIS:
	{
		Vec axis = rotationConstraintDirection();
		Vec quat = Vec(rotation[0], rotation[1], rotation[2]);
		quat.projectOnAxis(axis);
		rotation = Quaternion(quat, 2.0*acos(rotation[3]));
	}
		break;
	case AxisPlaneConstraint::FORBIDDEN:
		rotation = Quaternion(); // identity
		break;
	}
}

////////////////////////////////////////////////////////////////////////////////
//                               WorldConstraint                              //
////////////////////////////////////////////////////////////////////////////////

/*! Depending on translationConstraintType(), constrain \p translation to be along an axis or
  limited to a plane defined in the world coordinate system by
  translationConstraintDirection(). */
void WorldConstraint::constrainTranslation(Vec& translation, Frame* const frame)
{
	Vec proj;
	switch (translationConstraintType())
	{
	case AxisPlaneConstraint::FREE:
		break;
	case AxisPlaneConstraint::PLANE:
		if (frame->referenceFrame())
		{
			proj = frame->referenceFrame()->transformOf(translationConstraintDirection());
			translation.projectOnPlane(proj);
		}
		else
			translation.projectOnPlane(translationConstraintDirection());
		break;
	case AxisPlaneConstraint::AXIS:
		if (frame->referenceFrame())
		{
			proj = frame->referenceFrame()->transformOf(translationConstraintDirection());
			translation.projectOnAxis(proj);
		}
		else
			translation.projectOnAxis(translationConstraintDirection());
		break;
	case AxisPlaneConstraint::FORBIDDEN:
		translation = Vec(0.0, 0.0, 0.0);
		break;
	}
}

/*! When rotationConstraintType() is AxisPlaneConstraint::AXIS, constrain \p rotation to be a rotation
  around an axis whose direction is defined in the world coordinate system by
  rotationConstraintDirection(). */
void WorldConstraint::constrainRotation(Quaternion& rotation, Frame* const frame)
{
	switch (rotationConstraintType())
	{
	case AxisPlaneConstraint::FREE:
		break;
	case AxisPlaneConstraint::PLANE:
		break;
	case AxisPlaneConstraint::AXIS:
	{
		Vec quat(rotation[0], rotation[1], rotation[2]);
		Vec axis = frame->transformOf(rotationConstraintDirection());
		quat.projectOnAxis(axis);
		rotation = Quaternion(quat, 2.0*acos(rotation[3]));
		break;
	}
	case AxisPlaneConstraint::FORBIDDEN:
		rotation = Quaternion(); // identity
		break;
	}
}

////////////////////////////////////////////////////////////////////////////////
//                               CameraConstraint                              //
////////////////////////////////////////////////////////////////////////////////

/*! Creates a CameraConstraint, whose constrained directions are defined in the \p camera coordinate
  system. */
CameraConstraint::CameraConstraint(const Camera* const camera)
	: AxisPlaneConstraint(), camera_(camera)
{}

/*! Depending on translationConstraintType(), constrain \p translation to be along an axis or
  limited to a plane defined in the camera() coordinate system by
  translationConstraintDirection(). */
void CameraConstraint::constrainTranslation(Vec& translation, Frame* const frame)
{
	Vec proj;
	switch (translationConstraintType())
	{
	case AxisPlaneConstraint::FREE:
		break;
	case AxisPlaneConstraint::PLANE:
		proj = camera()->frame()->inverseTransformOf(translationConstraintDirection());
		if (frame->referenceFrame())
			proj = frame->referenceFrame()->transformOf(proj);
		translation.projectOnPlane(proj);
		break;
	case AxisPlaneConstraint::AXIS:
		proj = camera()->frame()->inverseTransformOf(translationConstraintDirection());
		if (frame->referenceFrame())
			proj = frame->referenceFrame()->transformOf(proj);
		translation.projectOnAxis(proj);
		break;
	case AxisPlaneConstraint::FORBIDDEN:
		translation = Vec(0.0, 0.0, 0.0);
		break;
	}
}

/*! When rotationConstraintType() is AxisPlaneConstraint::AXIS, constrain \p rotation to be a rotation
  around an axis whose direction is defined in the camera() coordinate system by
  rotationConstraintDirection(). */
void CameraConstraint::constrainRotation(Quaternion& rotation, Frame* const frame)
{
	switch (rotationConstraintType())
	{
	case AxisPlaneConstraint::FREE:
		break;
	case AxisPlaneConstraint::PLANE:
		break;
	case AxisPlaneConstraint::AXIS:
	{
		Vec axis = frame->transformOf(camera()->frame()->inverseTransformOf(rotationConstraintDirection()));
		Vec quat = Vec(rotation[0], rotation[1], rotation[2]);
		quat.projectOnAxis(axis);
		rotation = Quaternion(quat, 2.0*acos(rotation[3]));
	}
		break;
	case AxisPlaneConstraint::FORBIDDEN:
		rotation = Quaternion(); // identity
		break;
	}
}
