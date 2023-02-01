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

#ifndef _BASIC_TYPES_H_
#define _BASIC_TYPES_H_

#include "basic_common.h"
#include <string>
#include <vector>

#include <cmath>
#include <cstring>
#include <sstream>

#if (defined _WIN32) || (defined _WIN64)

#include <windows.h>
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#else

#include <unistd.h>

#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



//____________________________________________________________________________


/*
* A Namespace gathering typedefs for memory management.
*/

namespace Memory {

	typedef unsigned char  byte ;
	typedef unsigned char  word8 ;
	typedef unsigned short word16 ;
	typedef unsigned int   word32 ;

	typedef byte* pointer ;

	inline void clear(void* addr, size_t size) {
		::memset(addr, 0, size) ;
	}

	inline void copy(void* to, const void* from, size_t size) {
		::memcpy(to, from, size) ;
	}

} 

//_______________________________________________________________________

/**
* A namespace gathering typedefs
* corresponding to numbers. These types
* names have the form (u)int<size> or float<size>,
* where the (optional) u denotes an unsigned type,
* and the size is in bits.
*/

namespace Numeric {

	typedef char					int8 ;
	typedef short					int16 ;
	typedef int						int32 ;
#ifdef WIN32
	typedef __int64					int64 ;
#else
	typedef long long int			int64 ;	
#endif

	typedef unsigned char			uint8 ;
	typedef unsigned short			uint16 ;
	typedef unsigned int			uint32 ;

#ifdef WIN32
	typedef unsigned __int64		uint64 ;
#else
	typedef unsigned long long int	uint64 ;
#endif    

	typedef float					float32 ;
	typedef double					float64 ;

	typedef void*					pointer;

	// -----------------------------------------------------------------

	extern BASIC_API const float big_float ;
	extern BASIC_API const float small_float ;
	extern BASIC_API const double big_double ;
	extern BASIC_API const double small_double ;

	bool BASIC_API is_nan(float32 x) ;
	bool BASIC_API is_nan(float64 x) ;

	int32	BASIC_API random_int32() ;
	float32	BASIC_API random_float32() ;
	float64	BASIC_API random_float64() ;

	/* Random real number in the range [min, max] */
	inline float32 random_float32(float32 min, float32 max) { return min + random_float32() * (max - min); }
	inline float64 random_float64(float64 min, float64 max) { return min + random_float64() * (max - min); }

	/** Index of maximum of 2 values. */
	template < typename T >
	size_t index_of_max(T a, T b) {
		return a > b ? 0 : 1;
	}

	/** Index of maximum of 2 values by magnitude. */
	template < typename T >
	size_t index_of_max_abs(T a, T b) {
		return index_of_max(std::fabs(a), std::fabs(b));
	}

	/** Index of minimum of 2 values. */
	template < typename T >
	size_t index_of_min(T a, T b) {
		return a < b ? 0 : 1;
	}

	/** Index of minimum of 2 values by magnitude. */
	template < typename T >
	size_t index_of_min_abs(T a, T b) {
		return index_of_min(std::fabs(a), std::fabs(b));
	}

	/** Index of maximum of 3 values. */
	template < typename T >
	size_t index_of_max(T a, T b, T c) {
		return a > b ? (c > a ? 2 : 0) : (b > c ? 1 : 2);
	}

	/** Index of maximum of 3 values by magnitude. */
	template < typename T >
	size_t index_of_max_abs(T a, T b, T c) {
		return index_of_max(std::fabs(a), std::fabs(b), std::fabs(c));
	}

	/** Index of minimum of 3 values. */
	template < typename T >
	size_t index_of_min(T a, T b, T c) {
		return a < b ? (c < a ? 2 : 0) : (b < c ? 1 : 2);
	}

	/** Index of minimum of 3 values by magnitude. */
	template < typename T >
	size_t index_of_min_abs(T a, T b, T c) {
		return index_of_min(std::fabs(a), std::fabs(b), std::fabs(c));
	}

} 

//_______________________________________________________________________

namespace String {

	void BASIC_API split_string(
		const std::string& in, 
		char separator,
		std::vector<std::string>& out,
		bool skip_empty_fields = true ) ;

	void BASIC_API join_strings(
		const std::vector<std::string>& in,
		char separator,
		std::string& out
		) ;

	void BASIC_API join_strings(
		const std::vector<std::string>& in,
		const std::string& separator,
		std::string& out
		) ;


	std::string BASIC_API join_strings(
		const std::vector<std::string>& in,
		char separator
		) ;

	std::string BASIC_API join_strings(
		const std::vector<std::string>& in,
		const std::string& separator
		) ;

	// return the number of substrings that have been replaced.
	// NOTE: if 'isolated' set to true, only isolated (by spaces) 
	//       substrings will be replaced.
	int BASIC_API replace_substring(
		std::string& in, 
		const std::string& original_substring, 
		const std::string& new_substring,
		bool isolated = true
		) ;

	// -----------------------------------------------------------------

	template <typename T>
	std::string from_value(T v)	{
		std::ostringstream ss;
		ss << v << '\0';
		//ss << v;	
		return ss.str();
	}

	template <typename T>
	T to_value(const std::string& str) {
		std::istringstream ss(str);
		T result;
		return (ss >> result) ? result : 0;
	}

	// -----------------------------------------------------------------

	void BASIC_API to_lowercase(std::string& in) ;
	void BASIC_API to_uppercase(std::string& in) ;

	inline std::string BASIC_API char_to_string(char c) {
		char s[2] ;
		s[0] = c ;
		s[1] = '\0' ;
		return std::string(s) ;
	}

	std::string BASIC_API quote(const std::string& s, char quotes = '\"') ;

	// format example: "Fri Jan 09 11:39:32 2015"
	std::string BASIC_API from_current_time() ;
} 

//_______________________________________________________________________

#define nil 0

//_______________________________________________________________________

template <class T> 
inline T ogf_max(T x1, T x2) {
	return x1 > x2 ? x1 : x2;
}

template <class T> 
inline T ogf_min(T x1, T x2) {
	return x1 < x2 ? x1 : x2;
}

template <class T>
inline T ogf_max(T x1, T x2, T x3) {
	return ogf_max( ogf_max(x1, x2), x3 );
}

template <class T>
inline T ogf_min(T x1, T x2, T x3) {
	return ogf_min( ogf_min(x1, x2), x3 );
}

enum Sign { NEGATIVE=-1, ZERO=0, POSITIVE=1 } ;

template <class T> 
inline Sign ogf_sgn(T x) {
	return (x > 0) ? POSITIVE : (
		(x < 0) ? NEGATIVE : ZERO
		);
}

template <class T> 
inline T ogf_abs(T x)  {
	return (x > 0) ? x : -x;
}

template <class T> 
inline T ogf_sqr(T x)  {
	return x*x;
}

template <class T> 
inline void ogf_clamp(T& x, T min, T max) {
	if(x < min) {
		x = min ;
	} else if(x > max) {
		x = max ;
	}
}


template <class T> 
inline void ogf_swap(T& x, T& y) {
	T z = x ;
	x = y ;
	y = z ;
}


// round the given floating point number v to be num_digits.
template <class T>
inline T truncate_digits(const T& v, int num_digits) {
	T tmp = pow(10.0f, num_digits);
	long long des = static_cast<long long>((v < 0) ? (v * tmp - 0.5f) : (v * tmp + 0.5f));
	T result = T(des) / tmp;
	return result;
}


#endif
