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
#include "quaternion.h"
#include <stdlib.h> // RAND_MAX

// All the methods are declared inline in Quaternion.h
using namespace qglviewer;
using namespace std;

/*! Constructs a Quaternion that will rotate from the \p from direction to the \p to direction.

Note that this rotation is not uniquely defined. The selected axis is usually orthogonal to \p from
and \p to, minimizing the rotation angle. This method is robust and can handle small or almost identical vectors. */
Quaternion::Quaternion(const Vec& from, const Vec& to)
{
	const qreal epsilon = 1E-10;

	const qreal fromSqNorm = from.squaredNorm();
	const qreal toSqNorm   = to.squaredNorm();
	// Identity Quaternion when one vector is null
	if ((fromSqNorm < epsilon) || (toSqNorm < epsilon))
	{
		q[0]=q[1]=q[2]=0.0;
		q[3]=1.0;
	}
	else
	{
		Vec axis = cross(from, to);
		const qreal axisSqNorm = axis.squaredNorm();

		// Aligned vectors, pick any axis, not aligned with from or to
		if (axisSqNorm < epsilon)
			axis = from.orthogonalVec();

		qreal angle = asin(sqrt(axisSqNorm / (fromSqNorm * toSqNorm)));

		if (from*to < 0.0)
			angle = M_PI-angle;

		setAxisAngle(axis, angle);
	}
}

/*! Returns the image of \p v by the Quaternion inverse() rotation.

rotate() performs an inverse transformation. Same as inverse().rotate(v). */
Vec Quaternion::inverseRotate(const Vec& v) const
{
	return inverse().rotate(v);
}

/*! Returns the image of \p v by the Quaternion rotation.

See also inverseRotate() and operator*(const Quaternion&, const Vec&). */
Vec Quaternion::rotate(const Vec& v) const
{
	const qreal q00 = 2.0 * q[0] * q[0];
	const qreal q11 = 2.0 * q[1] * q[1];
	const qreal q22 = 2.0 * q[2] * q[2];

	const qreal q01 = 2.0 * q[0] * q[1];
	const qreal q02 = 2.0 * q[0] * q[2];
	const qreal q03 = 2.0 * q[0] * q[3];

	const qreal q12 = 2.0 * q[1] * q[2];
	const qreal q13 = 2.0 * q[1] * q[3];

	const qreal q23 = 2.0 * q[2] * q[3];

	return Vec((1.0 - q11 - q22)*v[0] + (      q01 - q23)*v[1] + (      q02 + q13)*v[2],
			(      q01 + q23)*v[0] + (1.0 - q22 - q00)*v[1] + (      q12 - q03)*v[2],
			(      q02 - q13)*v[0] + (      q12 + q03)*v[1] + (1.0 - q11 - q00)*v[2] );
}

/*! Set the Quaternion from a (supposedly correct) 3x3 rotation matrix.

  The matrix is expressed in European format: its three \e columns are the images by the rotation of
  the three vectors of an orthogonal basis. Note that OpenGL uses a symmetric representation for its
  matrices.

  setFromRotatedBasis() sets a Quaternion from the three axis of a rotated frame. It actually fills
  the three columns of a matrix with these rotated basis vectors and calls this method. */
void Quaternion::setFromRotationMatrix(const qreal m[3][3])
{
	// Compute one plus the trace of the matrix
	const qreal onePlusTrace = 1.0 + m[0][0] + m[1][1] + m[2][2];

	if (onePlusTrace > 1E-5)
	{
		// Direct computation
		const qreal s = sqrt(onePlusTrace) * 2.0;
		q[0] = (m[2][1] - m[1][2]) / s;
		q[1] = (m[0][2] - m[2][0]) / s;
		q[2] = (m[1][0] - m[0][1]) / s;
		q[3] = 0.25 * s;
	}
	else
	{
		// Computation depends on major diagonal term
		if ((m[0][0] > m[1][1])&(m[0][0] > m[2][2]))
		{
			const qreal s = sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2.0;
			q[0] = 0.25 * s;
			q[1] = (m[0][1] + m[1][0]) / s;
			q[2] = (m[0][2] + m[2][0]) / s;
			q[3] = (m[1][2] - m[2][1]) / s;
		}
		else
			if (m[1][1] > m[2][2])
			{
				const qreal s = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2.0;
				q[0] = (m[0][1] + m[1][0]) / s;
				q[1] = 0.25 * s;
				q[2] = (m[1][2] + m[2][1]) / s;
				q[3] = (m[0][2] - m[2][0]) / s;
			}
			else
			{
				const qreal s = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2.0;
				q[0] = (m[0][2] + m[2][0]) / s;
				q[1] = (m[1][2] + m[2][1]) / s;
				q[2] = 0.25 * s;
				q[3] = (m[0][1] - m[1][0]) / s;
			}
	}
	normalize();
}

#ifndef DOXYGEN
void Quaternion::setFromRotationMatrix(const float m[3][3])
{
	qWarning("setFromRotationMatrix now expects a double[3][3] parameter");

	qreal mat[3][3];
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			mat[i][j] = qreal(m[i][j]);

	setFromRotationMatrix(mat);
}

void Quaternion::setFromRotatedBase(const Vec& X, const Vec& Y, const Vec& Z)
{
	qWarning("setFromRotatedBase is deprecated, use setFromRotatedBasis instead");
	setFromRotatedBasis(X,Y,Z);
}
#endif

/*! Sets the Quaternion from the three rotated vectors of an orthogonal basis.

  The three vectors do not have to be normalized but must be orthogonal and direct (X^Y=k*Z, with k>0).

  \code
  Quaternion q;
  q.setFromRotatedBasis(X, Y, Z);
  // Now q.rotate(Vec(1,0,0)) == X and q.inverseRotate(X) == Vec(1,0,0)
  // Same goes for Y and Z with Vec(0,1,0) and Vec(0,0,1).
  \endcode

  See also setFromRotationMatrix() and Quaternion(const Vec&, const Vec&). */
void Quaternion::setFromRotatedBasis(const Vec& X, const Vec& Y, const Vec& Z)
{
	qreal m[3][3];
	qreal normX = X.norm();
	qreal normY = Y.norm();
	qreal normZ = Z.norm();

	for (int i=0; i<3; ++i)
	{
		m[i][0] = X[i] / normX;
		m[i][1] = Y[i] / normY;
		m[i][2] = Z[i] / normZ;
	}

	setFromRotationMatrix(m);
}

/*! Returns the axis vector and the angle (in radians) of the rotation represented by the Quaternion.
 See the axis() and angle() documentations. */
void Quaternion::getAxisAngle(Vec& axis, qreal& angle) const
{
	angle = 2.0 * acos(q[3]);
	axis = Vec(q[0], q[1], q[2]);
	const qreal sinus = axis.norm();
	if (sinus > 1E-8)
		axis /= sinus;

	if (angle > M_PI)
	{
		angle = 2.0 * qreal(M_PI) - angle;
		axis = -axis;
	}
}

/*! Returns the normalized axis direction of the rotation represented by the Quaternion.

It is null for an identity Quaternion. See also angle() and getAxisAngle(). */
Vec Quaternion::axis() const
{
	Vec res = Vec(q[0], q[1], q[2]);
	const qreal sinus = res.norm();
	if (sinus > 1E-8)
		res /= sinus;
	return (acos(q[3]) <= M_PI/2.0) ? res : -res;
}

/*! Returns the angle (in radians) of the rotation represented by the Quaternion.

 This value is always in the range [0-pi]. Larger rotational angles are obtained by inverting the
 axis() direction.

 See also axis() and getAxisAngle(). */
qreal Quaternion::angle() const
{
	const qreal angle = 2.0 * acos(q[3]);
	return (angle <= M_PI) ? angle : 2.0*M_PI - angle;
}

/*! Returns an XML \c QDomElement that represents the Quaternion.

 \p name is the name of the QDomElement tag. \p doc is the \c QDomDocument factory used to create
 QDomElement.

 When output to a file, the resulting QDomElement will look like:
 \code
 <name q0=".." q1=".." q2=".." q3=".." />
 \endcode

 Use initFromDOMElement() to restore the Quaternion state from the resulting \c QDomElement. See
 also the Quaternion(const QDomElement&) constructor.

 See the Vec::domElement() documentation for a complete QDomDocument creation and saving example.

 See also Frame::domElement(), Camera::domElement(), KeyFrameInterpolator::domElement()... */
QDomElement Quaternion::domElement(const QString& name, QDomDocument& document) const
{
	QDomElement de = document.createElement(name);
	de.setAttribute("q0", QString::number(q[0]));
	de.setAttribute("q1", QString::number(q[1]));
	de.setAttribute("q2", QString::number(q[2]));
	de.setAttribute("q3", QString::number(q[3]));
	return de;
}

/*! Restores the Quaternion state from a \c QDomElement created by domElement().

 The \c QDomElement should contain the \c q0, \c q1 , \c q2 and \c q3 attributes. If one of these
 attributes is missing or is not a number, a warning is displayed and these fields are respectively
 set to 0.0, 0.0, 0.0 and 1.0 (identity Quaternion).

 See also the Quaternion(const QDomElement&) constructor. */
void Quaternion::initFromDOMElement(const QDomElement& element)
{
	Quaternion q(element);
	*this = q;
}

/*! Constructs a Quaternion from a \c QDomElement representing an XML code of the form
  \code< anyTagName q0=".." q1=".." q2=".." q3=".." />\endcode

  If one of these attributes is missing or is not a number, a warning is displayed and the associated
  value is respectively set to 0, 0, 0 and 1 (identity Quaternion).

  See also domElement() and initFromDOMElement(). */
Quaternion::Quaternion(const QDomElement& element)
{
	QStringList attribute;
	attribute << "q0" << "q1" << "q2" << "q3";
	for (int i=0; i<attribute.size(); ++i)
		q[i] = DomUtils::qrealFromDom(element, attribute[i], ((i<3)?0.0:1.0));
}

/*! Returns the Quaternion associated 4x4 OpenGL rotation matrix.

 Use \c glMultMatrixd(q.matrix()) to apply the rotation represented by Quaternion \c q to the
 current OpenGL matrix.

 See also getMatrix(), getRotationMatrix() and inverseMatrix().

 \attention The result is only valid until the next call to matrix(). Use it immediately (as shown
 above) or consider using getMatrix() instead.

 \attention The matrix is given in OpenGL format (row-major order) and is the transpose of the
 actual mathematical European representation. Consider using getRotationMatrix() instead. */
const GLdouble* Quaternion::matrix() const
{
	static GLdouble m[4][4];
	getMatrix(m);
	return (const GLdouble*)(m);
}

/*! Fills \p m with the OpenGL representation of the Quaternion rotation.

Use matrix() if you do not need to store this matrix and simply want to alter the current OpenGL
matrix. See also getInverseMatrix() and Frame::getMatrix(). */
void Quaternion::getMatrix(GLdouble m[4][4]) const
{
	const qreal q00 = 2.0 * q[0] * q[0];
	const qreal q11 = 2.0 * q[1] * q[1];
	const qreal q22 = 2.0 * q[2] * q[2];

	const qreal q01 = 2.0 * q[0] * q[1];
	const qreal q02 = 2.0 * q[0] * q[2];
	const qreal q03 = 2.0 * q[0] * q[3];

	const qreal q12 = 2.0 * q[1] * q[2];
	const qreal q13 = 2.0 * q[1] * q[3];

	const qreal q23 = 2.0 * q[2] * q[3];

	m[0][0] = 1.0 - q11 - q22;
	m[1][0] =       q01 - q23;
	m[2][0] =       q02 + q13;

	m[0][1] =       q01 + q23;
	m[1][1] = 1.0 - q22 - q00;
	m[2][1] =       q12 - q03;

	m[0][2] =       q02 - q13;
	m[1][2] =       q12 + q03;
	m[2][2] = 1.0 - q11 - q00;

	m[0][3] = 0.0;
	m[1][3] = 0.0;
	m[2][3] = 0.0;

	m[3][0] = 0.0;
	m[3][1] = 0.0;
	m[3][2] = 0.0;
	m[3][3] = 1.0;
}

/*! Same as getMatrix(), but with a \c GLdouble[16] parameter. See also getInverseMatrix() and Frame::getMatrix(). */
void Quaternion::getMatrix(GLdouble m[16]) const
{
	static GLdouble mat[4][4];
	getMatrix(mat);
	int count = 0;
	for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			m[count++] = mat[i][j];
}

/*! Fills \p m with the 3x3 rotation matrix associated with the Quaternion.

  See also getInverseRotationMatrix().

  \attention \p m uses the European mathematical representation of the rotation matrix. Use matrix()
  and getMatrix() to retrieve the OpenGL transposed version. */
void Quaternion::getRotationMatrix(qreal m[3][3]) const
{
	static GLdouble mat[4][4];
	getMatrix(mat);
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			// Beware of transposition
			m[i][j] = qreal(mat[j][i]);
}

/*! Returns the associated 4x4 OpenGL \e inverse rotation matrix. This is simply the matrix() of the
  inverse().

  \attention The result is only valid until the next call to inverseMatrix(). Use it immediately (as
  in \c glMultMatrixd(q.inverseMatrix())) or use getInverseMatrix() instead.

  \attention The matrix is given in OpenGL format (row-major order) and is the transpose of the
  actual mathematical European representation. Consider using getInverseRotationMatrix() instead. */
const GLdouble* Quaternion::inverseMatrix() const
{
	static GLdouble m[4][4];
	getInverseMatrix(m);
	return (const GLdouble*)(m);
}

/*! Fills \p m with the OpenGL matrix corresponding to the inverse() rotation.

Use inverseMatrix() if you do not need to store this matrix and simply want to alter the current
OpenGL matrix. See also getMatrix(). */
void Quaternion::getInverseMatrix(GLdouble m[4][4]) const
{
	inverse().getMatrix(m);
}

/*! Same as getInverseMatrix(), but with a \c GLdouble[16] parameter. See also getMatrix(). */
void Quaternion::getInverseMatrix(GLdouble m[16]) const
{
	inverse().getMatrix(m);
}

/*! \p m is set to the 3x3 \e inverse rotation matrix associated with the Quaternion.

 \attention This is the classical mathematical rotation matrix. The OpenGL format uses its
 transposed version. See inverseMatrix() and getInverseMatrix(). */
void Quaternion::getInverseRotationMatrix(qreal m[3][3]) const
{
	static GLdouble mat[4][4];
	getInverseMatrix(mat);
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			// Beware of transposition
			m[i][j] = qreal(mat[j][i]);
}


/*! Returns the slerp interpolation of Quaternions \p a and \p b, at time \p t.

 \p t should range in [0,1]. Result is \p a when \p t=0 and \p b when \p t=1.

 When \p allowFlip is \c true (default) the slerp interpolation will always use the "shortest path"
 between the Quaternions' orientations, by "flipping" the source Quaternion if needed (see
 negate()). */
Quaternion Quaternion::slerp(const Quaternion& a, const Quaternion& b, qreal t, bool allowFlip)
{
	qreal cosAngle = Quaternion::dot(a, b);

	qreal c1, c2;
	// Linear interpolation for close orientations
	if ((1.0 - fabs(cosAngle)) < 0.01)
	{
		c1 = 1.0 - t;
		c2 = t;
	}
	else
	{
		// Spherical interpolation
		qreal angle    = acos(fabs(cosAngle));
		qreal sinAngle = sin(angle);
		c1 = sin(angle * (1.0 - t)) / sinAngle;
		c2 = sin(angle * t) / sinAngle;
	}

	// Use the shortest path
	if (allowFlip && (cosAngle < 0.0))
		c1 = -c1;

	return Quaternion(c1*a[0] + c2*b[0], c1*a[1] + c2*b[1], c1*a[2] + c2*b[2], c1*a[3] + c2*b[3]);
}

/*! Returns the slerp interpolation of the two Quaternions \p a and \p b, at time \p t, using
  tangents \p tgA and \p tgB.

  The resulting Quaternion is "between" \p a and \p b (result is \p a when \p t=0 and \p b for \p
  t=1).

  Use squadTangent() to define the Quaternion tangents \p tgA and \p tgB. */
Quaternion Quaternion::squad(const Quaternion& a, const Quaternion& tgA, const Quaternion& tgB, const Quaternion& b, qreal t)
{
	Quaternion ab = Quaternion::slerp(a, b, t);
	Quaternion tg = Quaternion::slerp(tgA, tgB, t, false);
	return Quaternion::slerp(ab, tg, 2.0*t*(1.0-t), false);
}

/*! Returns the logarithm of the Quaternion. See also exp(). */
Quaternion Quaternion::log()
{
	qreal len = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

	if (len < 1E-6)
		return Quaternion(q[0], q[1], q[2], 0.0);
	else
	{
		qreal coef = acos(q[3]) / len;
		return Quaternion(q[0]*coef, q[1]*coef, q[2]*coef, 0.0);
	}
}

/*! Returns the exponential of the Quaternion. See also log(). */
Quaternion Quaternion::exp()
{
	qreal theta = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

	if (theta < 1E-6)
		return Quaternion(q[0], q[1], q[2], cos(theta));
	else
	{
		qreal coef = sin(theta) / theta;
		return Quaternion(q[0]*coef, q[1]*coef, q[2]*coef, cos(theta));
	}
}

/*! Returns log(a. inverse() * b). Useful for squadTangent(). */
Quaternion Quaternion::lnDif(const Quaternion& a, const Quaternion& b)
{
	Quaternion dif = a.inverse()*b;
	dif.normalize();
	return dif.log();
}

/*! Returns a tangent Quaternion for \p center, defined by \p before and \p after Quaternions.

 Useful for smooth spline interpolation of Quaternion with squad() and slerp(). */
Quaternion Quaternion::squadTangent(const Quaternion& before, const Quaternion& center, const Quaternion& after)
{
	Quaternion l1 = Quaternion::lnDif(center,before);
	Quaternion l2 = Quaternion::lnDif(center,after);
	Quaternion e;
	for (int i=0; i<4; ++i)
		e.q[i] = -0.25 * (l1.q[i] + l2.q[i]);
	e = center*(e.exp());

	// if (Quaternion::dot(e,b) < 0.0)
	// e.negate();

	return e;
}

ostream& operator<<(ostream& o, const Quaternion& Q)
{
	return o << Q[0] << '\t' << Q[1] << '\t' << Q[2] << '\t' << Q[3];
}

/*! Returns a random unit Quaternion.

You can create a randomly directed unit vector using:
\code
Vec randomDir = Quaternion::randomQuaternion() * Vec(1.0, 0.0, 0.0); // or any other Vec
\endcode

\note This function uses rand() to create pseudo-random numbers and the random number generator can
be initialized using srand().*/
Quaternion Quaternion::randomQuaternion()
{
	// The rand() function is not very portable and may not be available on your system.
	// Add the appropriate include or replace by an other random function in case of problem.
	qreal seed = rand()/(qreal)RAND_MAX;
	qreal r1 = sqrt(1.0 - seed);
	qreal r2 = sqrt(seed);
	qreal t1 = 2.0 * M_PI * (rand()/(qreal)RAND_MAX);
	qreal t2 = 2.0 * M_PI * (rand()/(qreal)RAND_MAX);
	return Quaternion(sin(t1)*r1, cos(t1)*r1, sin(t2)*r2, cos(t2)*r2);
}
