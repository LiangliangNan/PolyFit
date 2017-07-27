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

#ifndef _ATTRIBUTE_STORE_H_
#define _ATTRIBUTE_STORE_H_


#include "basic_common.h"
#include "raw_attribute_store.h"
#include "attribute_life_cycle.h"
#include "../basic/smart_pointer.h"
#include "../basic/counted.h"

#include <typeinfo>



class AttributeManager ;
class Record ;

/**
* stores an attribute, and knows how to construct,copy,destroy
* instances of the attribute. This class should not be used
* directly by client code.
*/
class BASIC_API AttributeStore : public Counted, public RawAttributeStore {
public:

	AttributeStore(
		AttributeLifeCycle* life_cycle,
		AttributeManager* manager = nil
		) : RawAttributeStore( life_cycle->item_size() ), 
		life_cycle_(life_cycle), manager_(nil) {
			bind(manager) ;
	}

	virtual ~AttributeStore() ;

	void construct(
		Memory::pointer addr, Record* record = 0
		) {
			life_cycle_->construct(addr,record) ;
	}

	void destroy(
		Memory::pointer addr, Record* record = 0
		) {
			life_cycle_->destroy(addr,record) ;
	}

	void copy(
		Memory::pointer lhs, Record* record_lhs,
		Memory::pointer rhs, const Record* record_rhs
		) {
			life_cycle_->copy(lhs,record_lhs,rhs,record_rhs) ;
	}

	void copy_construct(
		Memory::pointer lhs, Record* record_lhs,
		Memory::pointer rhs, const Record* record_rhs
		) {
			life_cycle_->copy_construct(lhs,record_lhs,rhs,record_rhs) ;
	}

	void bind(AttributeManager* manager) ;
	AttributeManager* attribute_manager() const { return manager_; }

	virtual const std::type_info& attribute_type_id() const = 0 ;

	/** returns an empty AttributeStore() of the same type. */
	virtual AttributeStore* clone() = 0 ;

protected:
	AttributeLifeCycle_var life_cycle_ ;
	AttributeManager* manager_ ;
} ;

typedef SmartPointer<AttributeStore> AttributeStore_var ;

//_________________________________________________________

/**
* A typed AttributeStore, templated by the 
* Record class and the Attribute class. This
* is used for static and dynamic type checking
* in the AttributeManager.
*/
template <class ATTRIBUTE> 
class GenericAttributeStore : public AttributeStore {
public:
	GenericAttributeStore(
		AttributeLifeCycle* life_cycle,
		AttributeManager* manager = nil
		) : AttributeStore(life_cycle, manager) { 
	}        
	virtual ~GenericAttributeStore() { }
	virtual const std::type_info& attribute_type_id() const {
		return typeid(ATTRIBUTE) ;
	}
	virtual AttributeStore* clone() {
		return new GenericAttributeStore<ATTRIBUTE>(life_cycle_) ;
	}
} ;



#endif

