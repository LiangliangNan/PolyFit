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


#ifndef _BASIC_ASSERTION_H_
#define _BASIC_ASSERTION_H_

#include "basic_common.h"
#include <string>



void BASIC_API ogf_assertion_failed(
									const std::string& condition_string,
									const std::string& file, int line ) ;

void BASIC_API ogf_range_assertion_failed(
	double value, double min_value, double max_value, 
	const std::string& file, int line ) ;

void BASIC_API ogf_should_not_have_reached(
	const std::string& file, int line ) ;


//________________________________________________________________________________


// Three levels of assert:
// use ogf_assert() and ogf_range_assert()               for non-expensive asserts
// use ogf_debug_assert() and ogf_debug_range_assert()   for expensive asserts
// use ogf_parano_assert() and ogf_parano_range_assert() for very exensive asserts

#define ogf_assert(x) {									\
	if(!(x)) {											\
	::ogf_assertion_failed(#x,__FILE__, __LINE__) ;		\
	}													\
} 

#define ogf_range_assert(x,min_val,max_val) {			\
	if(((x) < (min_val)) || ((x) > (max_val))) {		\
	::ogf_range_assertion_failed(x, min_val, max_val,	\
	__FILE__, __LINE__									\
	) ;													\
	}													\
}

#define ogf_assert_not_reached {						\
	::ogf_should_not_have_reached(__FILE__, __LINE__) ;	\
}

#ifndef NDEBUG
#define ogf_debug_assert(x) ogf_assert(x)
#define ogf_debug_range_assert(x,min_val,max_val) ogf_range_assert(x,min_val,max_val)
#else
#define ogf_debug_assert(x) 
#define ogf_debug_range_assert(x,min_val,max_val) 
#endif

#ifdef	PARANOID_DEBUG
#define ogf_parano_assert(x) ogf_assert(x)
#define ogf_parano_range_assert(x,min_val,max_val) ogf_range_assert(x,min_val,max_val)
#else
#define ogf_parano_assert(x) 
#define ogf_parano_range_assert(x,min_val,max_val) 
#endif

#endif
