/*
*  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
*  Copyright (C) 2000-2005 INRIA - Project ALICE
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy - levy@loria.fr
*
*     Project ALICE
*     LORIA, INRIA Lorraine,
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX
*     FRANCE
*
*  Note that the GNU General Public License does not permit incorporating
*  the Software into proprietary programs.
*
* As an exception to the GPL, Graphite can be linked with the following (non-GPL) libraries:
*     Qt, SuperLU, WildMagic and CGAL
*/


#ifndef _MATH_TYPES_VECG_H_
#define _MATH_TYPES_VECG_H_

#include "math_common.h"
#include "../basic/basic_types.h"
#include "../basic/assertions.h"
#include <iostream>
#include <cfloat>


template <int DIM, class T> 
class vecng {
public:
	typedef vecng<DIM,T> thisclass ;

	vecng() { for(unsigned int i=0; i<DIM; i++) { data_[i] = T(0); } }

	// This one should never be called : a template constructor cannot be a copy constructor
	template<class T2> explicit vecng(const vecng<DIM,T2>& rhs) {  
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] = T(rhs[i]) ;
		}
	}

	// to avoid compilation problems
	template<class T2, int DIM2> explicit vecng(const vecng<DIM2,T2>& rhs) {  
		ogf_debug_assert(DIM2==DIM);
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] = T(rhs[i]) ;
		}
	}
	template<class T2> explicit vecng(const T2* rhs) {
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] = T(rhs[i]) ;
		}
	}

	thisclass& operator=(const thisclass& rhs) {
		memcpy(data_, rhs.data(), DIM*sizeof(T));
		return *this ;
	}

	unsigned int dimension() const { return (unsigned int)DIM; } 

	T* data()             { return data_ ; }

	const T* data() const { return data_ ; }

	inline T& operator[](unsigned int idx) {
		ogf_debug_assert(idx < DIM) ;
		return data()[idx] ;
	}

	inline const T& operator[](unsigned int idx) const {
		ogf_debug_assert(idx < DIM) ;
		return data()[idx] ;
	}

	inline T length2() const { 
		T result = T(0) ;
		for(unsigned int i=0; i<DIM; i++) {
			result += data_[i]*data_[i] ;
		}
		return result ;
	}

	inline T length() const {
		return sqrt(length2()) ;
	}

	inline T distance2(const thisclass &rhs) const { 
		T result = T(0) ;
		for(unsigned int i=0; i<DIM; i++) {
			T val = rhs.data_[i]-data_[i];
			result += val*val ;
		}
		return result ;
	}

	// operators
	inline thisclass& operator+=(const thisclass& v) { 
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] += v.data_[i] ;
		}
		return *this ; 
	}

	inline thisclass& operator-=(const thisclass& v) { 
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] -= v.data_[i] ;
		}
		return *this ; 
	}

	inline thisclass& operator*=(const thisclass& v) {
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] *= v.data_[i] ;
		}
		return *this ; 
	}

	inline thisclass& operator/=(const thisclass& v) { 
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] /= v.data_[i] ;
		}
		return *this ; 
	}

	template <class T2> inline thisclass& operator*=(T2 s) { 
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] *= T(s) ;
		}
		return *this ; 
	}

	template <class T2> inline thisclass& operator/=(T2 s) { 
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] /= T(s) ;
		}
		return *this ; 
	}

	inline thisclass operator+ (const thisclass& v) const {
		thisclass result(*this) ;
		for(unsigned int i=0; i<DIM; i++) {
			result.data_[i] += v.data_[i] ;
		}
		return result ;
	}

	inline thisclass operator- (const thisclass& v) const {
		thisclass result(*this) ;
		for(unsigned int i=0; i<DIM; i++) {
			result.data_[i] -= v.data_[i] ;
		}
		return result ;
	}

	template <class T2> inline thisclass operator* (T2 s) const {
		thisclass result(*this) ;
		for(unsigned int i=0; i<DIM; i++) {
			result.data_[i] *= T(s) ;
		}
		return result ;
	}

	template <class T2> inline thisclass operator/ (T2 s) const {
		thisclass result(*this) ;
		for(unsigned int i=0; i<DIM; i++) {
			result.data_[i] /= T(s) ;
		}
		return result ;
	}

	inline thisclass operator- () const {
		thisclass result ;
		for(unsigned int i=0; i<DIM; i++) {
			result.data_[i] = -data_[i] ;
		}
		return result ;
	}


private:
	T data_[DIM] ;
} ;


template <int DIM, class T> inline T dot(const vecng<DIM,T>& v1, const vecng<DIM,T>& v2) {  
	T result = 0 ;
	for(unsigned int i=0; i<DIM; i++) {
		result += v1[i]*v2[i] ;
	}
	return result ;
}

template <int DIM, class T> inline vecng<DIM,T> operator-(const vecng<DIM,T>& v1) { 
	vecng<DIM,T> result ;
	for(unsigned int i=0; i<DIM; i++) {
		result[i] = -v1[i] ;
	}
	return result ;
}

template <class T2, int DIM, class T> inline vecng<DIM,T> operator*(T2 s, const vecng<DIM,T>& v) { 
	vecng<DIM,T> result ;
	for(unsigned int i=0; i<DIM; i++) {
		result[i] = T(s)*v[i] ;
	}
	return result ;
}

template <int DIM, class T> inline vecng<DIM,T> operator+(const vecng<DIM,T>& v1, const vecng<DIM,T>& v2) { 
	vecng<DIM,T> result ;
	for(unsigned int i=0; i<DIM; i++) {
		result[i] = v1[i]+v2[i] ;
	}
	return result ;
}

template <int DIM, class T> inline vecng<DIM,T> operator-(const vecng<DIM,T>& v1, const vecng<DIM,T>& v2) { 
	vecng<DIM,T> result ;
	for(unsigned int i=0; i<DIM; i++) {
		result[i] = v1[i]-v2[i] ;
	}
	return result ;
}

// Compatibility with GLSL
template <int DIM, class T> inline T length(const vecng<DIM,T>& v)	{ return v.length() ;  }
template <int DIM, class T> inline T length2(const vecng<DIM,T>& v) { return v.length2() ; }
template <int DIM, class T> inline T distance2(const vecng<DIM,T>& v1, const vecng<DIM,T>& v2)	{ return v2.distance2(v1) ; }
template <int DIM, class T> inline T distance(const vecng<DIM,T>& v1, const vecng<DIM,T>& v2)	{ return length(v2 - v1) ;  }
template <int DIM, class T> inline vecng<DIM,T> normalize(const vecng<DIM,T>& v) { 
	T s = length(v) ;
	if(s > 1e-30) { s = T(1)/s ; }
	return s * v ; 
}
template <int DIM, class T> inline vecng<DIM,T> mix(const vecng<DIM,T>& v1, const vecng<DIM,T>& v2, T s) { 
	return (T(1) - s) * v1 + s * v2 ; 
}


//-------------------- vec2 -------------------------------------------------------------------

template <class T> 
class vecng<2,T> {
public:
	typedef vecng<2,T> thisclass ;

	vecng() : x(0), y(0) { }
	vecng(T x_in, T y_in) : x(x_in), y(y_in) {  }
	template<class T2> explicit vecng(const vecng<2,T2> & v)
		: x(v.x), y(v.y) {}

	template<class T2> explicit vecng(const T2* v)
		: x(v[0]), y(v[1]) {}

	inline T length2() const { return x*x+y*y ; }
	inline T length() const { return sqrt(x*x+y*y) ; }
	inline T distance2(const thisclass& rhs) const {
		T dx = rhs.x-x;
		T dy = rhs.y-y;			
		return dx*dx+dy*dy;
	}

	// operators
	inline thisclass& operator+=(const thisclass& v) { x += v.x ; y += v.y ; return *this ; }
	inline thisclass& operator-=(const thisclass& v) { x -= v.x ; y -= v.y ; return *this ; }
	inline thisclass& operator*=(const thisclass& v) { x *= v.x ; y *= v.y ; return *this ; }
	inline thisclass& operator/=(const thisclass& v) { x /= v.x ; y /= v.y ; return *this ; }
	template <class T2> inline thisclass& operator*=(T2 s) { x *= T(s) ; y *= T(s) ; return *this ; }
	template <class T2> inline thisclass& operator/=(T2 s) { x /= T(s) ; y /= T(s) ; return *this ; }

	inline thisclass operator+ (const thisclass& v) const {return thisclass(x+v.x, y+v.y); }
	inline thisclass operator- (const thisclass& v) const {return thisclass(x-v.x, y-v.y); }
	template <class T2> inline thisclass operator* (T2 s) const {return thisclass(x*T(s), y*T(s)); }
	template <class T2> inline thisclass operator/ (T2 s) const {return thisclass(x/T(s), y/T(s)); }
	inline thisclass operator- () const {return thisclass(-x, -y);}

	unsigned int dimension() const { return (unsigned int)2; }

	T* data() { return &x ; }
	const T* data() const { return &x ; }

	inline T& operator[](unsigned int idx) {
		ogf_debug_assert(idx < 2) ;
		return data()[idx] ;
	}

	inline const T& operator[](unsigned int idx) const {
		ogf_debug_assert(idx < 2) ;
		return data()[idx] ;
	}

	T x ;
	T y ;
} ;

template <class T> inline T dot(const vecng<2,T>& v1, const vecng<2,T>& v2) {  return v1.x*v2.x + v1.y*v2.y ;  }

template <class T> inline  T det(const vecng<2,T>& v1, const vecng<2,T>& v2) {
	return v1.x*v2.y - v1.y*v2.x ;
}


template <class T> inline vecng<2,T> operator-(const vecng<2,T>& v1) { 
	return vecng<2,T>(-v1.x, -v1.y) ; 
}
template <class T2, class T> inline vecng<2,T> operator*(T2 s, const vecng<2,T>& v) { 
	return vecng<2,T>(T(s)*v.x, T(s)*v.y) ;   
}

template <class T> inline vecng<2,T> operator+(const vecng<2,T>& v1, const vecng<2,T>& v2) { 
	return vecng<2,T>(v1.x+v2.x, v1.y+v2.y) ; 
}

template <class T> inline vecng<2,T> operator-(const vecng<2,T>& v1, const vecng<2,T>& v2) { 
	return vecng<2,T>(v1.x - v2.x, v1.y - v2.y) ; 
}


//---------------- vec3 ------------------------------------------------------------------------

template <class T> 
class vecng<3,T> {
public:
	typedef vecng<3,T> thisclass ;

	vecng() : x(0), y(0), z(0) { }
	vecng(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {  }
	template<class T2> explicit vecng(const vecng<3,T2> & v) : x(v.x), y(v.y), z(v.z) {}
	template<class T2> explicit vecng(const T2* v)
		: x(v[0]), y(v[1]), z(v[2]) {}

	inline T length2() const { return x*x+y*y+z*z ; }
	inline T length() const { return sqrt(x*x+y*y+z*z) ; }
	inline T distance2(const thisclass& rhs) const {
		T dx = rhs.x-x;
		T dy = rhs.y-y;			
		T dz = rhs.z-z;			
		return dx*dx+dy*dy+dz*dz;
	}

	// operators
	inline thisclass& operator+=(const thisclass& v) { x += v.x ; y += v.y ; z += v.z ; return *this ; }
	inline thisclass& operator-=(const thisclass& v) { x -= v.x ; y -= v.y ; z -= v.z ; return *this ; }
	inline thisclass& operator*=(const thisclass& v) { x *= v.x ; y *= v.y ; z *= v.z ; return *this ; }
	inline thisclass& operator/=(const thisclass& v) { x /= v.x ; y /= v.y ; z /= v.z ; return *this ; }
	template <class T2> inline thisclass& operator*=(T2 s) { x *= T(s) ; y *= T(s) ; z *= T(s) ; return *this ; }
	template <class T2> inline thisclass& operator/=(T2 s) { x /= T(s) ; y /= T(s) ; z /= T(s) ; return *this ; }


	inline thisclass operator+ (const thisclass& v) const {return thisclass(x+v.x, y+v.y, z+v.z); }
	inline thisclass operator- (const thisclass& v) const {return thisclass(x-v.x, y-v.y, z-v.z); }
	template <class T2> inline thisclass operator* (T2 s) const {return thisclass(x*T(s), y*T(s), z*T(s)); }
	template <class T2> inline thisclass operator/ (T2 s) const {return thisclass(x/T(s), y/T(s), z/T(s)); }

	inline thisclass operator- () const {return thisclass(-x, -y, -z);}
	unsigned int dimension() const { return (unsigned int)3; }

	T* data() { return &x ; }
	const T* data() const { return &x ; }

	inline T& operator[](unsigned int idx) {
		ogf_debug_assert(idx < 3) ;
		return data()[idx] ;
	}

	inline const T& operator[](unsigned int idx) const {
		ogf_debug_assert(idx < 3) ;
		return data()[idx] ;
	}

	T x ;
	T y ;
	T z ;
} ;

template <class T> inline T dot(const vecng<3,T>& v1, const vecng<3,T>& v2) {  
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;  
}

template <class T> inline  vecng<3,T> cross(const vecng<3,T>& v1, const vecng<3,T>& v2) {
	return vecng<3,T>(
		v1.y*v2.z - v1.z*v2.y,
		v1.z*v2.x - v1.x*v2.z,
		v1.x*v2.y - v1.y*v2.x
		) ;
}

template <class T> inline vecng<3,T> operator-(const vecng<3,T>& v1) { return vecng<3,T>(-v1.x, -v1.y, -v1.z) ; }
template <class T2, class T> inline vecng<3,T> operator*(T2 s, const vecng<3,T>& v) { 
	return vecng<3,T>(T(s)*v.x, T(s)*v.y, T(s)*v.z) ;   
}

template <class T> inline vecng<3,T> operator+(const vecng<3,T>& v1, const vecng<3,T>& v2) { 
	return vecng<3,T>(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z) ; 
}

template <class T> inline vecng<3,T> operator-(const vecng<3,T>& v1, const vecng<3,T>& v2) { 
	return vecng<3,T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z) ; 
}


// ----------------- vec4 ----------------------------------------------------------------------------------

template <class T> 
class vecng<4,T> {
public:
	typedef vecng<4,T> thisclass ;

	vecng() : x(0), y(0), z(0), w(0) { }
	vecng(T x_in, T y_in, T z_in, T w_in) : x(x_in), y(y_in), z(z_in), w(w_in) {  }
	template<class T2> explicit vecng(const vecng<4,T2> & v)
		: x(v.x), y(v.y), z(v.z), w(v.w) {}
	template<class T2> explicit vecng(const T2* v)
		: x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}

	inline T length2() const { return x*x+y*y+z*z+w*w ; }
	inline T length() const { return sqrt(x*x+y*y+z*z+w*w) ; }
	inline T distance2(const thisclass& rhs) const {
		T dx = rhs.x-x;
		T dy = rhs.y-y;			
		T dz = rhs.z-z;			
		T dw = rhs.w-w;
		return dx*dx+dy*dy+dz*dz+dw*dw;
	}

	unsigned int dimension() const { return (unsigned int)4; }

	// operators
	inline thisclass& operator+=(const thisclass& v) { x += v.x ; y += v.y ; z += v.z ; w += v.w ; return *this ; }
	inline thisclass& operator-=(const thisclass& v) { x -= v.x ; y -= v.y ; z -= v.z ; w -= v.w ; return *this ; }
	inline thisclass& operator*=(const thisclass& v) { x *= v.x ; y *= v.y ; z *= v.z ; w *= v.w ; return *this ; }
	inline thisclass& operator/=(const thisclass& v) { x /= v.x ; y /= v.y ; z /= v.z ; w /= v.w ; return *this ; }
	template <class T2> inline thisclass& operator*=(T2 s) { 
		x *= T(s) ; y *= T(s) ; z *= T(s) ; w *= T(s) ; return *this ; 
	}
	template <class T2> inline thisclass& operator/=(T2 s) { 
		x /= T(s) ; y /= T(s) ; z /= T(s) ; w /= T(s) ; return *this ; 
	}
	inline thisclass operator+ (const thisclass& v) const {return thisclass(x+v.x, y+v.y, z+v.z, w+v.w); }
	inline thisclass operator- (const thisclass& v) const {return thisclass(x-v.x, y-v.y, z-v.z, w-v.w); }
	template <class T2> inline thisclass operator* (T2 s) const {return thisclass(x*T(s), y*T(s), z*T(s), w*T(s)); }
	template <class T2> inline thisclass operator/ (T2 s) const {return thisclass(x/T(s), y/T(s), z/T(s), w/T(s)); }
	inline thisclass operator- () const {return thisclass(-x, -y, -z, -w);}


	T* data() { return &x ; }
	const T* data() const { return &x ; }

	inline T& operator[](unsigned int idx) {
		ogf_debug_assert(idx < 4) ;
		return data()[idx] ;
	}

	inline const T& operator[](unsigned int idx) const {
		ogf_debug_assert(idx < 4) ;
		return data()[idx] ;
	}

	T x ;
	T y ;
	T z ;
	T w ; 
} ;

template <class T> inline T dot(const vecng<4,T>& v1, const vecng<4,T>& v2) {  
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;  
}

template <class T> inline vecng<4,T> operator-(const vecng<4,T>& v1) { return vecng<4,T>(-v1.x, -v1.y, -v1.z, -v1.w) ; }
template <class T2, class T> inline vecng<4,T> operator*(T2 s, const vecng<4,T>& v) { 
	return vecng<4,T>(T(s)*v.x, T(s)*v.y, T(s)*v.z, T(s)*v.w) ;   
}

template <class T> inline vecng<4,T> operator+(const vecng<4,T>& v1, const vecng<4,T>& v2) { 
	return vecng<4,T>(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z, v1.w+v2.w) ; 
}

template <class T> inline vecng<4,T> operator-(const vecng<4,T>& v1, const vecng<4,T>& v2) { 
	return vecng<4,T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w) ; 
}





template <int DIM, class T> inline std::ostream& operator<<(std::ostream& out, const ::vecng<DIM,T>& v) {
	for(unsigned int i=0; i<DIM; i++) {
		out << v[i] << " " ;
	}
	return out ;
}

template <int DIM, class T> inline std::istream& operator>>(std::istream& in, ::vecng<DIM,T>& v) {
	for(unsigned int i=0; i<DIM; i++) {
		in >> v[i] ;
	}
	return in ;
}

template <class T> inline std::ostream& operator<<(std::ostream& out, const ::vecng<2,T>& v) {
	return out << v.x << " " << v.y ;
}

template <class T> inline std::istream& operator>>(std::istream& in, ::vecng<2,T>& v) {
	return in >> v.x >> v.y ;
}

template <class T> inline std::ostream& operator<<(std::ostream& out, const ::vecng<3,T>& v) {
	return out << v.x << " " << v.y << " " << v.z  ;
}

template <class T> inline std::istream& operator>>(std::istream& in, ::vecng<3,T>& v) {
	return in >> v.x >> v.y >> v.z ;
}

template <class T> inline std::ostream& operator<<(std::ostream& out, const ::vecng<4,T>& v) {
	return out << v.x << " " << v.y << " " << v.z  << " " << v.w ;
}

template <class T> inline std::istream& operator>>(std::istream& in, ::vecng<4,T>& v) {
	return in >> v.x >> v.y >> v.z >> v.w ;
}

#endif

