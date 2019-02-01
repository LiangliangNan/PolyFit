#ifndef _MATH_QUATERNION_H_
#define _MATH_QUATERNION_H_

#include "math_common.h"
#include "math_types.h"

/*The Quaternion class represents 3D rotations and orientations.

The Quaternion is an appropriate (although not very intuitive) representation for 3D rotations and
orientations. Many tools are provided to ease the definition of a Quaternion: see constructors,
set_axis_angle(), set_from_rotation_matrix(), set_from_rotated_basis().

You can apply the rotation represented by the Quaternion to 3D points using rotate() and
inverse_rotate().

You can apply the Quaternion q rotation to the OpenGL matrices using: glMultMatrixd(q.matrix());
equvalent to glRotate(q.angle()*180.0/M_PI, q.axis().x, q.axis().y, q.axis().z);

Internal representation
The internal representation of a Quaternion corresponding to a rotation around axis, with an angle
alpha is made of four doubles q[i]:
{q[0],q[1],q[2]} = sin(alpha/2) * {axis[0],axis[1],axis[2]}
q[3] = cos(alpha/2)
Note that certain implementations place the cosine term in first position (instead of last here).

The Quaternion is always normalized, so that its inverse() is actually its conjugate.
*/

class MATH_API Quaternion
{
public:
	/* Defining a Quaternion */
	/* Default constructor, builds an identity rotation. */
	Quaternion()
	{ q[0]=q[1]=q[2]=0.0;  q[3]=1.0; }

	/* Constructor from rotation axis (non null) and angle (in radians). See also set_axis_angle(). */
	Quaternion(const vec3& axis, double angle)
	{
		set_axis_angle(axis, angle);
	}

	/* Constructs a Quaternion that will rotate from the \p from direction to the \p to direction.
	Note that this rotation is not uniquely defined. The selected axis is usually orthogonal to \p from
	and \p to, minimizing the rotation angle. This method is robust and can handle small or almost identical vectors. */
	Quaternion(const vec3& from, const vec3& to);

	/* Constructor from the four values of a Quaternion. First three values are axis*sin(angle/2) and
	last one is cos(angle/2).

	\attention The identity Quaternion is Quaternion(0,0,0,1) and not Quaternion(0,0,0,0) (which is
	not unitary). The default Quaternion() creates such identity Quaternion. */
	Quaternion(double q0, double q1, double q2, double q3)
	{ q[0]=q0;    q[1]=q1;    q[2]=q2;    q[3]=q3; }

	/* Copy constructor. */
	Quaternion(const Quaternion& Q)
	{ for (int i=0; i<4; ++i) q[i] = Q.q[i]; }

	/* Equal operator. */
	Quaternion& operator=(const Quaternion& Q) {
		for (int i=0; i<4; ++i)
			q[i] = Q.q[i];
		return (*this);
	}

	/* Sets the Quaternion as a rotation of axis and angle (in radians).
	\p axis does not need to be normalized. A null axis will result in an identity Quaternion. */
	void set_axis_angle(const vec3& axis, double angle);

	/* Sets the Quaternion value. */
	void set_value(double q0, double q1, double q2, double q3)
	{ q[0]=q0;    q[1]=q1;    q[2]=q2;    q[3]=q3; }

	/* Set the Quaternion from a (supposedly correct) 3x3 rotation matrix.
	The matrix is expressed in European format: its three columns are the images by the rotation of
	the three vectors of an orthogonal basis. Note that OpenGL uses a symmetric representation for its
	matrices.*/
	void set_from_rotation_matrix(const double m[3][3]);

	/* set_from_rotated_basis() sets a Quaternion from the three axis of a rotated frame. It actually fills
	the three columns of a matrix with these rotated basis vectors and calls this method. */
	void set_from_rotated_basis(const vec3& X, const vec3& Y, const vec3& Z);

	/* Returns the normalized axis direction of the rotation represented by the Quaternion.
	It is null for an identity Quaternion. See also angle() and get_axis_angle(). */
	vec3 axis() const;

	/* Returns the angle (in radians) of the rotation represented by the Quaternion.
	This value is always in the range [0-pi]. Larger rotational angles are obtained by inverting the
	axis() direction. See also axis() and get_axis_angle(). */
	double angle() const;

	/* Returns the axis vector and the angle (in radians) of the rotation represented by the Quaternion.*/
	void get_axis_angle(vec3& axis, double& angle) const;

	/* Bracket operator, with a constant return value. i must range in [0..3]. */
	double operator[](int i) const { return q[i]; }

	/* Bracket operator returning an l-value. i must range in [0..3]. */
	double& operator[](int i) { return q[i]; }

	/* Rotation computations */

	/* Returns the composition of the a and b rotations.
	The order is important. When applied to a Vec v (see operator*(const Quaternion&, const Vec&)
	and rotate()) the resulting Quaternion acts as if b was applied first and then a was applied. 
	This is obvious since the image v' of v by the composited rotation satisfies: 
	v'= (a*b) * v = a * (b*v) 
	Note that a*b usually differs from b*a.
	\attention For efficiency reasons, the resulting Quaternion is not normalized. Use normalize() in
	case of numerical drift with small rotation composition. */
	friend Quaternion operator*(const Quaternion& a, const Quaternion& b)
	{
		return Quaternion(
			a.q[3]*b.q[0] + b.q[3]*a.q[0] + a.q[1]*b.q[2] - a.q[2]*b.q[1],
			a.q[3]*b.q[1] + b.q[3]*a.q[1] + a.q[2]*b.q[0] - a.q[0]*b.q[2],
			a.q[3]*b.q[2] + b.q[3]*a.q[2] + a.q[0]*b.q[1] - a.q[1]*b.q[0],
			a.q[3]*b.q[3] - b.q[0]*a.q[0] - a.q[1]*b.q[1] - a.q[2]*b.q[2]
			);
	}

	/* Quaternion rotation is composed with q.
	See operator*(), since this is equivalent to this = this * q.
	\note For efficiency reasons, the resulting Quaternion is not normalized.
	You may normalize() it after each application in case of numerical drift. */
	Quaternion& operator*=(const Quaternion &q) {
		*this = (*this)*q;
		return *this;
	}

	/* Returns the image of v by the rotation q.
	Same as q.rotate(v). See rotate() and inverse_rotate(). */
	friend vec3 operator*(const Quaternion& q, const vec3& v) { return q.rotate(v); }
	
	/* Returns the image of v by the Quaternion rotation.
	See also inverse_rotate() and operator*(const Quaternion&, const Vec&). */
	vec3 rotate(const vec3& v) const;

	/* Returns the image of v by the Quaternion inverse() rotation.
	rotate() performs an inverse transformation. Same as inverse().rotate(v). */
	vec3 inverse_rotate(const vec3& v) const;

    /* Inversion
    Returns the inverse Quaternion (inverse rotation).
	Result has a negated axis() direction and the same angle(). A composition (see operator*()) of a
	Quaternion and its inverse() results in an identity function.
	Use invert() to actually modify the Quaternion. */
	Quaternion inverse() const { return Quaternion(-q[0], -q[1], -q[2], q[3]); }

	/* Inverses the Quaternion (same rotation angle(), but negated axis()).
	See also inverse(). */
	void invert() { q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; }

	/* Negates all the coefficients of the Quaternion.
	This results in an other representation of the same rotation (opposite rotation angle, but with
	a negated axis direction: the two cancel out). However, note that the results of axis() and
	angle() are unchanged after a call to this method since angle() always returns a value in [0,pi].
	This method is mainly useful for Quaternion interpolation, so that the spherical interpolation 
	takes the shortest path on the unit sphere. See slerp() for details. */
	void negate() { invert(); q[3] = -q[3]; }

	/* Normalizes the Quaternion coefficients.
	This method should not need to be called since we only deal with unit Quaternions. This is however
	useful to prevent numerical drifts, especially with small rotational increments. See also
	normalized(). */
	double normalize() {
		const double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
		for (int i=0; i<4; ++i)
			q[i] /= norm;
		return norm;
	}

	/* Returns a normalized version of the Quaternion.
	See also normalize(). */
	Quaternion normalized() const {
		double Q[4];
		const double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
		for (int i=0; i<4; ++i)
			Q[i] = q[i] / norm;
		return Quaternion(Q[0], Q[1], Q[2], Q[3]);
	}

	/* Returns the Quaternion associated 4x4 OpenGL rotation matrix.
	Use glMultMatrixd(q.matrix()) to apply the rotation represented by Quaternion q to the
	current OpenGL matrix. See also get_matrix(), get_rotation_matrix() and inverse_matrix().
	\attention The result is only valid until the next call to matrix(). Use it immediately (as shown
	above) or consider using get_matrix() instead.
	\attention The matrix is given in OpenGL format (row-major order) and is the transpose of the
	actual mathematical European representation. Consider using get_rotation_matrix() instead. */
	const double* matrix() const;
	/* Fills m with the OpenGL representation of the Quaternion rotation.
	Use matrix() if you do not need to store this matrix and simply want to alter the current OpenGL
	matrix. See also get_inverse_matrix(). */
	void get_matrix(double m[4][4]) const;
	/* Same as get_matrix(), but with a double[16] parameter. See also get_inverse_matrix(). */
	void get_matrix(double m[16]) const;

	/* Fills m with the 3x3 rotation matrix associated with the Quaternion.
	See also get_inverse_rotation_matrix().
	\attention m uses the European mathematical representation of the rotation matrix. Use matrix()
	and get_matrix() to retrieve the OpenGL transposed version. */
	void get_rotation_matrix(double m[3][3]) const;

	/* Returns the associated 4x4 OpenGL inverse rotation matrix. This is simply the matrix() of the inverse().
	\attention The result is only valid until the next call to inverse_matrix(). Use it immediately (as
	in glMultMatrixd(q.inverse_matrix())) or use get_inverse_matrix() instead.
	\attention The matrix is given in OpenGL format (row-major order) and is the transpose of the
	actual mathematical European representation. Consider using get_inverse_rotation_matrix() instead. */
	const double* inverse_matrix() const;
	/* Fills m with the OpenGL matrix corresponding to the inverse() rotation.
	Use inverse_matrix() if you do not need to store this matrix and simply want to alter the current
	OpenGL matrix. See also get_matrix(). */
	void get_inverse_matrix(double m[4][4]) const;
	/* Same as get_inverse_matrix(), but with a double[16] parameter. See also get_matrix(). */
	void get_inverse_matrix(double m[16]) const;

	/* m is set to the 3x3 inverse rotation matrix associated with the Quaternion.
	\attention This is the classical mathematical rotation matrix. The OpenGL format uses its
	transposed version. See inverse_matrix() and get_inverse_matrix(). */
	void get_inverse_rotation_matrix(double m[3][3]) const;

	/* Slerp(Spherical Linear intERPolation) interpolation */
	/* Returns the slerp interpolation of Quaternions a and b, at time t. t should range in [0,1]. 
	Result is a when t=0 and b when t=1.
	When allowFlip is true (default) the slerp interpolation will always use the "shortest path" between 
	the Quaternions' orientations, by "flipping" the source Quaternion if needed (see negate()). */
	static Quaternion slerp(const Quaternion& a, const Quaternion& b, double t, bool allowFlip=true);

	/* Returns the slerp interpolation of the two Quaternions a and b, at time t, using tangents tgA and tgB.
	The resulting Quaternion is "between" a and b (result is a when t=0 and b for t=1).
	Use squad_tangent() to define the Quaternion tangents tgA and tgB. */
	static Quaternion squad(const Quaternion& a, const Quaternion& tgA, const Quaternion& tgB, const Quaternion& b, float t);

	/* Returns the "dot" product of a and b: a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]. */
	static double dot(const Quaternion& a, const Quaternion& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]; }

	/* Returns the logarithm of the Quaternion. See also exp(). */
	Quaternion log();
	/* Returns the exponential of the Quaternion. See also log(). */
	Quaternion exp();

	/* Returns log(a. inverse() * b). Useful for squad_tangent(). */
	static Quaternion ln_dif(const Quaternion& a, const Quaternion& b);
	/* Returns a tangent Quaternion for \p center, defined by \p before and \p after Quaternions.
	Useful for smooth spline interpolation of Quaternion with squad() and slerp(). */
	static Quaternion squad_tangent(const Quaternion& before, const Quaternion& center, const Quaternion& after);

	/* Returns a random unit Quaternion.
	You can create a randomly directed unit vector using:
	Vec randomDir = Quaternion::random_quaternion() * Vec(1.0, 0.0, 0.0); // or any other Vec
	\note This function uses rand() to create pseudo-random numbers and the random number generator can
	be initialized using srand().*/
	static Quaternion random_quaternion();

private:
	/* The internal data representation is private, use operator[] to access values. */
	double q[4];
};


std::ostream& operator<<(std::ostream& o, const Quaternion&);

#endif // _MATH_QUATERNION_H_
