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


#ifndef __RAT_H_
#define __RAT_H_

#include "basic_common.h" 
#include "raw_attribute_store.h"

/**
* RAT (Record Allocation Table) is used internally by
* the AttributeManager. It manages the free list.
*/

class BASIC_API RAT : public RawAttributeStore {
public:
	RAT() : RawAttributeStore(sizeof(RecordId)) { }

	bool is_full() const { return free_list_.is_nil() ; }

	virtual void clear() ; 

	/**
	* returns a new unique RecordId. If the RAT is full,
	* fails in an assertion check (grow() should be called
	* before).
	*/
	RecordId new_record_id() ;

	/**
	* adds a RecordId to the free list. It can
	* then be returned by a subsequent call to
	* new_record_id()
	*/
	void delete_record_id(RecordId record) ;

	/**
	* adds all the items of the new chunk to 
	* the free list.
	*/
	virtual void grow() ;

protected:
	RecordId& cell(unsigned int chunk, unsigned int offset) {
		return *reinterpret_cast<RecordId*>(
			RawAttributeStore::data(chunk,offset)
			) ;
	}

	RecordId& cell(RecordId index) {
		return cell(index.chunk(),index.offset()) ;
	}

	const RecordId& cell(unsigned int chunk, unsigned int offset) const {
		return *reinterpret_cast<RecordId*>(
			RawAttributeStore::data(chunk,offset)
			) ;
	}

	const RecordId& cell(RecordId index) const {
		return cell(index.chunk(),index.offset()) ;
	}

private:
	RecordId free_list_ ;

	friend class AttributeManager ;
} ;


#endif

