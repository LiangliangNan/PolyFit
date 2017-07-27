#ifndef _CGAL_TYPES_H_
#define _CGAL_TYPES_H_


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include "../math/math_types.h"



typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;

typedef K::Point_2					Point2;
typedef K::Vector_2					Vector2;

typedef K::Point_3					Point3;
typedef K::Vector_3					Vector3;

typedef K::Line_3					Line3;

typedef K::Plane_3					Plane3;

typedef CGAL::Polygon_2<K>			Polygon2;


inline Point2	to_cgal_point(const vec2& p)	{ return Point2(p.x, p.y); }
inline Vector2  to_cgal_vector(const vec2& v)	{ return Vector2(v.x, v.y); }

inline Point3	to_cgal_point(const vec3& p)	{ return Point3(p.x, p.y, p.z); }
inline Vector3  to_cgal_vector(const vec3& v)	{ return Vector3(v.x, v.y, v.z); }

inline vec2		to_my_point(const Point2& p)	{ return vec2(p.x(), p.y()); }
inline vec2		to_my_vector(const Vector2& v)	{ return vec2(v.x(), v.y()); }

inline vec3		to_my_point(const Point3& p)	{ return vec3(p.x(), p.y(), p.z()); }
inline vec3		to_my_vector(const Vector3& v)	{ return vec3(v.x(), v.y(), v.z()); }

inline Plane3	to_cgal_plane(const Plane3d& plane) {
	return Plane3(to_cgal_point(plane.point()), to_cgal_vector(plane.normal()));
}

// the coordinate system is defined by (orig, base1, base2)
inline Point2  convert_to_2d(const Point3& orig, const Vector3& base1, const Vector3& base2, const Point3& p) {
	Vector3 vec = p - orig;
	double x = vec * base1;
	double y = vec * base2;
	return Point2(x, y);
}


inline Point3  convert_to_3d(const Point3& orig, const Vector3& base1, const Vector3& base2, const Point2& p) {
	return orig + base1 * p.x() + base2 * p.y();
}


#endif