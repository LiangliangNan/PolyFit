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


#include "rat.h"


void RAT::clear() {
	RawAttributeStore::clear() ;
	free_list_.forget() ;
}

RecordId RAT::new_record_id() {
	ogf_attribute_assert(!is_full()) ;
	RecordId result = free_list_ ;
	free_list_ = cell(free_list_) ;
	cell(result).unfree() ;
	return result ;
}

void RAT::delete_record_id(RecordId record) {
	RecordId& ref = cell(record) ;
	ogf_attribute_assert(!ref.is_free()) ;
	ref = free_list_ ;
	ref.free() ;
	free_list_ = record ;
}


void RAT::grow() {
	RawAttributeStore::grow() ;
	unsigned int chunk = nb_chunks() - 1 ;
	for(unsigned int i=0; i<CHUNK_SIZE-1; i++) {
		cell(chunk,i)=RecordId(chunk,i+1,true) ;
	}
	cell(chunk,CHUNK_SIZE-1) = free_list_ ;
	free_list_ = RecordId(chunk,0) ;
}

