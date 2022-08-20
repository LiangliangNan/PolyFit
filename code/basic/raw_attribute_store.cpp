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


#include "raw_attribute_store.h"
#include <typeinfo>
#include <string.h>


void RawAttributeStore::clear() {
	for(
		std::vector<Memory::pointer>::iterator 
		it=data_.begin(); it!=data_.end(); it++
		) {
			// Paranoid stuff: reset freed memory to zero
#ifdef OGF_ATTRIBUTE_CHECK
			Memory::clear(*it, CHUNK_SIZE * item_size_) ;
#endif
			delete[] *it ;
	}
	data_.clear() ;
}

RawAttributeStore::~RawAttributeStore() {
	clear() ;
}

void RawAttributeStore::grow() {
	Memory::pointer chunk = new Memory::byte[
		CHUNK_SIZE * item_size_
	] ;
	// Paranoid stuff: initialize allocated memory to zero
#ifdef OGF_ATTRIBUTE_CHECK
	Memory::clear(chunk, CHUNK_SIZE * item_size_) ;
#endif
	data_.push_back(chunk) ;
	//        std::cerr << "RawAttributeStore (" << typeid(*this).name() 
	//                  << ") grow" << std::endl ;
}