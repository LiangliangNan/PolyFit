
#ifndef _COLOR_H_
#define _COLOR_H_
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

#include "../basic/basic_types.h"
#include "../basic/assertions.h"

#include <cassert>
#include <iostream>


template <class T> 
class GenericColor
{
public:
	GenericColor(T r=0, T g=0, T b=0, T a=1.0f) ;
	template<class T2> explicit GenericColor(const T2* v, T a=1.0f) {
		components_[0] = v[0] ;
		components_[1] = v[1] ;
		components_[2] = v[2] ;
		components_[3] = a ;
	}

	T r() const ;
	T g() const ;
	T b() const ;
	T a() const ;

	void set_r(T r) ;
	void set_g(T g) ;
	void set_b(T b) ;
	void set_a(T a) ;

	void set(T r, T g, T b, T a=1.0f) ; 	
	
	T& operator[](int i) ;
	const T& operator[](int i) const ;

	// Low-level access
	const T* data() const ;
	T*       data() ;

	// color key 
	bool operator<(const GenericColor& rhs) const {
		if(components_[0] < rhs.components_[0]) { return true ; }
		if(components_[0] > rhs.components_[0]) { return false ; }
		if(components_[1] < rhs.components_[1]) { return true ; }
		if(components_[1] > rhs.components_[1]) { return false ; }
		if(components_[2] < rhs.components_[2]) { return true ; }
		if(components_[2] > rhs.components_[2]) { return false ; }
		return (components_[3] < rhs.components_[3]) ; 
	}

private:
	T components_[4] ;
} ;



template <class T> inline
GenericColor<T>::GenericColor(T r, T g, T b, T a) {
	components_[0] = r ;
	components_[1] = g ;
	components_[2] = b ;
	components_[3] = a ;
}

template <class T> inline
T GenericColor<T>::r() const {
	return components_[0] ;
}

template <class T> inline
T GenericColor<T>::g() const {
	return components_[1] ;
}

template <class T> inline
T GenericColor<T>::b() const {
	return components_[2] ;
}

template <class T> inline
T GenericColor<T>::a() const {
	return components_[3] ;
}

// Low-level access
template <class T> inline
const T* GenericColor<T>::data() const { 
	return components_; 
}

template <class T> inline
T* GenericColor<T>::data() {
	return components_; 
}

template <class T> inline
void GenericColor<T>::set(T r, T g, T b, T a) {
	components_[0] = r ;
	components_[1] = g ;
	components_[2] = b ;
	components_[3] = a ;
}

template <class T> inline
void GenericColor<T>::set_r(T r) {
	components_[0] = r ;
}

template <class T> inline
void GenericColor<T>::set_g(T g) {
	components_[1] = g ;
}

template <class T> inline
void GenericColor<T>::set_b(T b) {
	components_[2] = b ;
}

template <class T> inline
void GenericColor<T>::set_a(T a) {
	components_[3] = a ;
}

template <class T> inline
T& GenericColor<T>::operator[](int i) {
	ogf_assert(i >= 0 && i <= 3) ;
	return components_[i] ;
}

template <class T> inline
const T& GenericColor<T>::operator[](int i) const {
	ogf_assert(i >= 0 && i <= 3) ;
	return components_[i] ;
}


template <class T> inline
std::ostream& operator<<(std::ostream& output, const GenericColor<T>& color) {
	return output << 
		color[0] << " " << color[1] << " " << color[2] << " " << color[3] ;
}

template <class T> inline
std::istream& operator>>(std::istream& input, GenericColor<T>& color) {
	return input >> color[0] >> color[1] >> color[2] >> color[3] ;
}



//_______________________ Colors

typedef GenericColor<Numeric::int8>    Color_int8 ;
typedef GenericColor<Numeric::uint8>   Color_uint8 ;
typedef GenericColor<Numeric::int16>   Color_int16 ;
typedef GenericColor<Numeric::uint16>  Color_uint16 ;
typedef GenericColor<Numeric::int32>   Color_int32 ;
typedef GenericColor<Numeric::uint32>  Color_uint32 ;
typedef GenericColor<Numeric::float32> Color_float32 ;
typedef GenericColor<Numeric::float64> Color_float64 ;

typedef		Color_float32	Color ;


inline Color random_color(bool allow_dark = false) {
	float min_rgb = 0.3f;
	if (allow_dark)
		min_rgb = 0.0f;
	
	return Color(Numeric::random_float32(min_rgb, 1.0f), Numeric::random_float32(min_rgb, 1.0f), Numeric::random_float32(min_rgb, 1.0f));
}

inline Color fused_color(const Color& c1, float w1, const Color& c2, float w2) {
	float r = (c1.r() * w1 + c2.r() * w2) / (w1 + w2);
	float g = (c1.g() * w1 + c2.g() * w2) / (w1 + w2);
	float b = (c1.b() * w1 + c2.b() * w2) / (w1 + w2);
	float a = (c1.a() * w1 + c2.a() * w2) / (w1 + w2);
	return Color(r, g, b, a);
}

#endif
