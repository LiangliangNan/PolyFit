
#ifndef _COUNTED_H_
#define _COUNTED_H_
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

#include "basic_common.h"
#include "basic_types.h"
#include "smart_pointer.h"
#include "assertions.h"



//____________________________________________________________________________

/**
* This is the base class to be used for objects having
* "reference count" memory management. They can be 
* referred to by using SmartPointer<T>, calling ref()
* and unref() when necessary.
* @see SmartPointer
*/

class BASIC_API Counted {

public:
	Counted() ;
	virtual ~Counted() ;

	void ref() const ;
	void unref() const ;
	bool is_shared() const ;

	static void ref(const Counted* counted) ;
	static void unref(const Counted* counted) ;

protected:
private:
	int nb_refs_ ;
} ;

//____________________________________________________________________________

inline Counted::Counted() : nb_refs_(0) {
}

inline void Counted::ref() const {
	Counted* non_const_this = (Counted *)this ;
	non_const_this->nb_refs_++ ;
}

inline void Counted::unref() const {
	Counted* non_const_this = (Counted *)this ;    
	non_const_this->nb_refs_-- ;

	ogf_assert(nb_refs_ >= 0) ;

	if(nb_refs_ == 0) {
		delete this ;
	}
}

inline bool Counted::is_shared() const {
	return (nb_refs_ > 1) ;
}



#endif
