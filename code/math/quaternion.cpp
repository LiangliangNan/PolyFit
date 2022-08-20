
#include "quaternion.h"
#include <cstdlib> // RAND_MAX
#include <iostream>


Quaternion::Quaternion(const vec3& from, const vec3& to)
{
	const double epsilon = 1E-10f;

	const double fromSqNorm = from.length2();
	const double toSqNorm   = to.length2();
	// Identity Quaternion when one vector is null
	if ((fromSqNorm < epsilon) || (toSqNorm < epsilon))
	{
		q[0]=q[1]=q[2]=0.0;
		q[3]=1.0;
	}
	else
	{
		vec3 axis = cross(from, to);
		const double axisSqNorm = axis.length2();

		// Aligned vectors, pick any axis, not aligned with from or to
		if (axisSqNorm < epsilon) 
			axis = Geom::perpendicular(from);

		double angle = std::asin(std::sqrt(axisSqNorm / (fromSqNorm * toSqNorm)));

		if (::dot(from, to) < 0.0)
			angle = M_PI-angle;

		set_axis_angle(axis, angle);
	}
}


void Quaternion::set_axis_angle(const vec3& axis, double angle)
{
	const double norm = axis.length();
	if (norm < 1E-8)
	{
		// Null rotation
		q[0] = 0.0;      q[1] = 0.0;      q[2] = 0.0;      q[3] = 1.0;
	}
	else
	{
		const double sin_half_angle = std::sin(angle / 2.0);
		q[0] = sin_half_angle*axis[0]/norm;
		q[1] = sin_half_angle*axis[1]/norm;
		q[2] = sin_half_angle*axis[2]/norm;
		q[3] = cos(angle / 2.0);
	}
}


vec3 Quaternion::inverse_rotate(const vec3& v) const
{
	return inverse().rotate(v);
}


vec3 Quaternion::rotate(const vec3& v) const
{
	const double q00 = 2.0l * q[0] * q[0];
	const double q11 = 2.0l * q[1] * q[1];
	const double q22 = 2.0l * q[2] * q[2];

	const double q01 = 2.0l * q[0] * q[1];
	const double q02 = 2.0l * q[0] * q[2];
	const double q03 = 2.0l * q[0] * q[3];

	const double q12 = 2.0l * q[1] * q[2];
	const double q13 = 2.0l * q[1] * q[3];

	const double q23 = 2.0l * q[2] * q[3];

	return vec3(
		(1.0 - q11 - q22)*v[0] + (      q01 - q23)*v[1] + (      q02 + q13)*v[2],
		(      q01 + q23)*v[0] + (1.0 - q22 - q00)*v[1] + (      q12 - q03)*v[2],
		(      q02 - q13)*v[0] + (      q12 + q03)*v[1] + (1.0 - q11 - q00)*v[2] );
}


void Quaternion::set_from_rotation_matrix(const double m[3][3])
{
	// Compute one plus the trace of the matrix
	const double onePlusTrace = 1.0 + m[0][0] + m[1][1] + m[2][2];

	if (onePlusTrace > 1E-5)
	{
		// Direct computation
		const double s = std::sqrt(onePlusTrace) * 2.0;
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
			const double s = std::sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2.0; 
			q[0] = 0.25 * s;
			q[1] = (m[0][1] + m[1][0]) / s; 
			q[2] = (m[0][2] + m[2][0]) / s; 
			q[3] = (m[1][2] - m[2][1]) / s;
		}
		else
			if (m[1][1] > m[2][2])
			{ 
				const double s = std::sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2.0; 
				q[0] = (m[0][1] + m[1][0]) / s; 
				q[1] = 0.25 * s;
				q[2] = (m[1][2] + m[2][1]) / s; 
				q[3] = (m[0][2] - m[2][0]) / s;
			}
			else
			{ 
				const double s = std::sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2.0; 
				q[0] = (m[0][2] + m[2][0]) / s; 
				q[1] = (m[1][2] + m[2][1]) / s; 
				q[2] = 0.25 * s;
				q[3] = (m[0][1] - m[1][0]) / s;
			}
	}
	normalize();
}


void Quaternion::set_from_rotated_basis(const vec3& X, const vec3& Y, const vec3& Z)
{
	double m[3][3];
	double normX = X.length();
	double normY = Y.length();
	double normZ = Z.length();

	for (int i=0; i<3; ++i)
	{
		m[i][0] = X[i] / normX;
		m[i][1] = Y[i] / normY;
		m[i][2] = Z[i] / normZ;
	}

	set_from_rotation_matrix(m);
}


void Quaternion::get_axis_angle(vec3& axis, double& angle) const
{
	angle = 2.0 * std::acos(q[3]);
	axis = vec3(q[0], q[1], q[2]);
	const double sinus = axis.length();
	if (sinus > 1E-8)
		axis = axis/sinus;

	if (angle > M_PI)
	{
		angle = 2.0*M_PI - angle;
		axis = -axis;
	}
}


vec3 Quaternion::axis() const
{
	vec3 res = vec3(q[0], q[1], q[2]);
	const double sinus = res.length();
	if (sinus > 1E-8)
		res = res/sinus;
	return (std::acos(q[3]) <= M_PI/2.0) ? res : -res;
}


double Quaternion::angle() const
{
	const double angle = 2.0 * std::acos(q[3]);
	return (angle <= M_PI) ? angle : 2.0*M_PI - angle;
}


const double* Quaternion::matrix() const
{
	static double m[4][4];
	get_matrix(m);
	return (const double*)(m);
}


void Quaternion::get_matrix(double m[4][4]) const
{
	const double q00 = 2.0l * q[0] * q[0];
	const double q11 = 2.0l * q[1] * q[1];
	const double q22 = 2.0l * q[2] * q[2];

	const double q01 = 2.0l * q[0] * q[1];
	const double q02 = 2.0l * q[0] * q[2];
	const double q03 = 2.0l * q[0] * q[3];

	const double q12 = 2.0l * q[1] * q[2];
	const double q13 = 2.0l * q[1] * q[3];

	const double q23 = 2.0l * q[2] * q[3];

	m[0][0] = 1.0l - q11 - q22;
	m[1][0] =        q01 - q23;
	m[2][0] =        q02 + q13;

	m[0][1] =        q01 + q23;
	m[1][1] = 1.0l - q22 - q00;
	m[2][1] =        q12 - q03;

	m[0][2] =        q02 - q13;
	m[1][2] =        q12 + q03;
	m[2][2] = 1.0l - q11 - q00;

	m[0][3] = 0.0l;
	m[1][3] = 0.0l;
	m[2][3] = 0.0l;

	m[3][0] = 0.0l;
	m[3][1] = 0.0l;
	m[3][2] = 0.0l;
	m[3][3] = 1.0l;
}


void Quaternion::get_matrix(double m[16]) const
{
	static double mat[4][4];
	get_matrix(mat);
	int count = 0;
	for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			m[count++] = mat[i][j];
}


void Quaternion::get_rotation_matrix(double m[3][3]) const
{
	static double mat[4][4];
	get_matrix(mat);
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			// Beware of transposition
			m[i][j] = mat[j][i];
}


const double* Quaternion::inverse_matrix() const
{
	static double m[4][4];
	get_inverse_matrix(m);
	return (const double*)(m);
}


void Quaternion::get_inverse_matrix(double m[4][4]) const
{
	inverse().get_matrix(m);
}


void Quaternion::get_inverse_matrix(double m[16]) const
{
	inverse().get_matrix(m);
}


void Quaternion::get_inverse_rotation_matrix(double m[3][3]) const
{
	static double mat[4][4];
	get_inverse_matrix(mat);
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			// Beware of transposition
			m[i][j] = mat[j][i];
}


Quaternion Quaternion::slerp(const Quaternion& a, const Quaternion& b, double t, bool allowFlip)
{
	double cosAngle = Quaternion::dot(a, b);

	double c1, c2;
	// Linear interpolation for close orientations
	if ((1.0 - fabs(cosAngle)) < 0.01)
	{
		c1 = 1.0 - t;
		c2 = t;
	}
	else
	{
		// Spherical interpolation
		double angle    = std::acos(std::fabs(cosAngle));
		double sinAngle = std::sin(angle);
		c1 = std::sin(angle * (1.0 - t)) / sinAngle;
		c2 = std::sin(angle * t) / sinAngle;
	}

	// Use the shortest path
	if (allowFlip && (cosAngle < 0.0))
		c1 = -c1;

	return Quaternion(c1*a[0] + c2*b[0], c1*a[1] + c2*b[1], c1*a[2] + c2*b[2], c1*a[3] + c2*b[3]);
}


Quaternion Quaternion::squad(const Quaternion& a, const Quaternion& tgA, const Quaternion& tgB, const Quaternion& b, float t)
{
	Quaternion ab = Quaternion::slerp(a, b, t);
	Quaternion tg = Quaternion::slerp(tgA, tgB, t, false);
	return Quaternion::slerp(ab, tg, 2.0*t*(1.0-t), false);
}


Quaternion Quaternion::log()
{
	double len = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

	if (len < 1E-6)
		return Quaternion(q[0], q[1], q[2], 0.0);
	else
	{
		double coef = std::acos(q[3]) / len;
		return Quaternion(q[0]*coef, q[1]*coef, q[2]*coef, 0.0);
	}
}


Quaternion Quaternion::exp()
{
	double theta = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

	if (theta < 1E-6)
		return Quaternion(q[0], q[1], q[2], std::cos(theta));
	else
	{
		double coef = std::sin(theta) / theta;
		return Quaternion(q[0]*coef, q[1]*coef, q[2]*coef, std::cos(theta));
	}
}


Quaternion Quaternion::ln_dif(const Quaternion& a, const Quaternion& b)
{
	Quaternion dif = a.inverse()*b;
	dif.normalize();
	return dif.log();
}


Quaternion Quaternion::squad_tangent(const Quaternion& before, const Quaternion& center, const Quaternion& after)
{
	Quaternion l1 = Quaternion::ln_dif(center,before);
	Quaternion l2 = Quaternion::ln_dif(center,after);
	Quaternion e;
	for (int i=0; i<4; ++i)
		e.q[i] = -0.25 * (l1.q[i] + l2.q[i]);
	e = center*(e.exp());

	// if (Quaternion::dot(e,b) < 0.0)
	// e.negate();

	return e;
}


Quaternion Quaternion::random_quaternion()
{
	// The rand() function is not very portable and may not be available on your system.
	// Add the appropriate include or replace by an other random function in case of problem.
	double seed = rand()/(double)RAND_MAX;
	double r1 = std::sqrt(1.0 - seed);
	double r2 = std::sqrt(seed);
	double t1 = 2.0 * M_PI * (rand()/(double)RAND_MAX);
	double t2 = 2.0 * M_PI * (rand()/(double)RAND_MAX);
	return Quaternion(std::sin(t1)*r1, std::cos(t1)*r1, std::sin(t2)*r2, std::cos(t2)*r2);
}


std::ostream& operator<<(std::ostream& o, const Quaternion& Q)
{
	return o << Q[0] << '\t' << Q[1] << '\t' << Q[2] << '\t' << Q[3];
}
