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


#ifndef _RAW_ATTRIBUTE_STORE__
#define _RAW_ATTRIBUTE_STORE__

#include "basic_common.h"
#include "record_id.h"

#include <vector>
#include <set>


/**
* Base class for the low-level storage for attributes. Client
* code should not need to use this directly. The storage space
* is allocated in chunks, to make dynamic growing more efficient
* (i.e. without needing to copy the data).
*/
class BASIC_API RawAttributeStore {
public:
	enum { CHUNK_SIZE = RecordId::MAX_OFFSET + 1 } ;

	RawAttributeStore(unsigned int item_size) : item_size_(item_size) { }
	virtual void clear() ; 
	virtual ~RawAttributeStore() ;
	unsigned int item_size() const { return item_size_ ; }
	unsigned int nb_chunks() const { return (unsigned int)data_.size() ; }
	unsigned int capacity() const { 
		return (unsigned int) data_.size() * CHUNK_SIZE ;
	}

	Memory::pointer data(
		unsigned int chunk, unsigned int offset
		) const {
			ogf_attribute_assert(chunk < data_.size()) ;
			ogf_attribute_assert(offset < CHUNK_SIZE) ;
			return &(data_[chunk][offset * item_size_]) ;
	}

	Memory::pointer data(const Record& r) const {
		return data(r.record_id().chunk(), r.record_id().offset()) ;
	}

	virtual void grow() ;

private:
	unsigned int item_size_ ;
	std::vector<Memory::pointer> data_ ;
} ;




#endif

