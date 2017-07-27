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

#ifndef _ATTRIBUTE_LIFE_CYCLE_H_
#define _ATTRIBUTE_LIFE_CYCLE_H_

#include "basic_common.h"
#include "../basic/counted.h"
#include "../basic/smart_pointer.h"

#include <string.h>



class Record ;

/**
* manages the life cycle (creation, copy, destruction) of
* the attributes. Note that AttributeLifeCycle is meant to
* be subclassed. GenericAttributeLifeCycle provides a default
* implementation calling the constructors and the destructor
* of the objects.
*/

class BASIC_API AttributeLifeCycle : public Counted {
public:

	/**
	* @param item_size size of the attribute
	* @param notify if set, the virtual functions
	*    notify_xxx are called when the corresponding
	*    action is performed. All these functions should
	*    be overloaded.
	* @param pod (plain ordinary datatype) if set, bitwise
	*    copy is used rather than calling the virtual_xxx
	*    functions.
	*/
	AttributeLifeCycle(
		unsigned int item_size, bool notify, bool pod
		) : item_size_(item_size), notify_(notify), pod_(pod) {
	}

	void construct(
		Memory::pointer addr, Record* record = 0
		) ;
	void destroy(
		Memory::pointer addr, Record* record = 0
		) ;
	void copy(
		Memory::pointer lhs, Record* record_lhs,
		Memory::pointer rhs, const Record* record_rhs
		) ;
	void copy_construct(
		Memory::pointer lhs, Record* record_lhs,
		Memory::pointer rhs, const Record* record_rhs
		) ;

	void pod_construct(Memory::pointer addr) { 
		Memory::clear(addr, item_size()) ;
	}
	void pod_destroy(Memory::pointer addr) { 
		Memory::clear(addr, item_size()) ;
	}
	void pod_copy(
		Memory::pointer lhs, Memory::pointer rhs
		) {
			Memory::copy(lhs, rhs, item_size()) ;
	}
	void pod_copy_construct(
		Memory::pointer lhs, Memory::pointer rhs
		) {
			Memory::copy(lhs, rhs, item_size()) ;
	}

	virtual void virtual_construct(Memory::pointer addr) ;
	virtual void virtual_destroy(Memory::pointer addr) ;
	virtual void virtual_copy(
		Memory::pointer lhs, Memory::pointer rhs 
		) ;
	virtual void virtual_copy_construct(
		Memory::pointer lhs, Memory::pointer rhs 
		) ;

	/**
	* if notify is set, this function is called after each record
	* creation. This function is meant to be overloaded.
	*/
	virtual void notify_construct(
		Memory::pointer addr, Record* record
		) ;

	/**
	* if notify is set, this function is called before each record
	* destruction. This function is meant to be overloaded.
	*/
	virtual void notify_destroy(
		Memory::pointer addr, Record* record
		) ;

	/**
	* if notify is set, this function is called before each record
	* copy. This function is meant to be overloaded.
	*/
	virtual void notify_copy(
		Memory::pointer lhs, Record* record_lhs,
		Memory::pointer rhs, const Record* record_rhs
		) ;

	/**
	* if notify is set, this function is called after each record
	* creation using a copy constructor. This function is meant to 
	* be overloaded.
	*/
	virtual void notify_copy_construct(
		Memory::pointer lhs, Record* record_lhs,
		Memory::pointer rhs, const Record* record_rhs
		) ;

	unsigned int item_size() const { return item_size_ ; }

	virtual AttributeLifeCycle* clone() = 0 ;

protected:
	unsigned int item_size_ ;
	bool notify_ ;
	bool pod_ ;
} ;

typedef SmartPointer<AttributeLifeCycle> AttributeLifeCycle_var;

//_________________________________________________________

inline void AttributeLifeCycle::construct(
	Memory::pointer addr, Record* record
	) {
		if(pod_) {
			pod_construct(addr) ;
		} else {
			virtual_construct(addr) ;
		}
		if(notify_) {
			notify_construct(addr, record) ;
		}
}

inline void AttributeLifeCycle::destroy(
										Memory::pointer addr, Record* record
										) 
{
	if(notify_) {
		notify_destroy(addr, record) ;
	}
	if(pod_) {
		pod_destroy(addr) ;
	} else {
		virtual_destroy(addr) ;
	}
}

inline void AttributeLifeCycle::copy(
									 Memory::pointer lhs, Record* record_lhs,
									 Memory::pointer rhs, const Record* record_rhs
									 ) 
{
	if(notify_) {
		notify_copy(lhs, record_lhs, rhs, record_rhs) ;
	}
	if(pod_) {
		pod_copy(lhs, rhs) ;
	} else {
		virtual_copy(lhs, rhs) ;
	}
}

inline void AttributeLifeCycle::copy_construct(
	Memory::pointer lhs, Record* record_lhs,
	Memory::pointer rhs, const Record* record_rhs
	) {
		if(pod_) {
			pod_copy_construct(lhs, rhs) ;
		} else {
			virtual_copy_construct(lhs,rhs) ;
		}
		if(notify_) {
			notify_copy_construct(lhs, record_lhs, rhs, record_rhs) ;
		}
}


//____________________________________________________________________________

/**
* defines default life cycle operations using the constructors,
* operator= and destructor of the T type. The RECORD type is
* used for static type checking.
*/

template <class ATTRIBUTE> 
class GenericAttributeLifeCycle : public AttributeLifeCycle {
public:
	GenericAttributeLifeCycle(
		bool notify = false, bool pod = false
		) : AttributeLifeCycle(sizeof(ATTRIBUTE), notify, pod) {
	}

	virtual void virtual_construct(Memory::pointer addr) ;
	virtual void virtual_destroy(Memory::pointer addr) ;
	virtual void virtual_copy(
		Memory::pointer lhs, Memory::pointer rhs 
		) ;
	virtual void virtual_copy_construct(
		Memory::pointer lhs, Memory::pointer rhs 
		) ;

	virtual AttributeLifeCycle* clone()  ;
} ;



//_____________________________________________________________________________________

// If your compiler barks somewhere in this file, it is
// likely that it lacks the new.h header, in that case,
// uncomment the following line:
// #define DOES_NOT_HAVE_NEW_H

#ifdef DOES_NOT_HAVE_NEW_H
inline void *operator new(size_t, void* p) { return p; }
#else
#include <new>
#endif

#include <string.h>
#include <typeinfo>

// The following code makes use of C++ oddities,
// enabling to separate memory management from
// object initialization and termination. 
// The correct way to call a constructor for an
// object of class T located at address addr is:
//    new(addr) T()
// which is the placement syntax for operator new.
// This means that new() uses the address addr instead
// of allocating space for the object.
// The correct way to call a destructor for an 
// object of class T located at address addr is:
//    addr-> ~T()
// Note that the placement syntax of operator new and
// direct calls to the destructor DO NOT allocate NOR
// deallocate memory.


//_________________________________________________________

template <class T>     
void GenericAttributeLifeCycle<T>::virtual_construct(
	Memory::pointer addr
	) {
		//  Initialize memory to zero, so that uninitialized 
		// pointers can be more easily detected.
		Memory::clear(addr, item_size()) ;
		// Placement syntax of operator new (see comments at
		// the beginning of this file).
		new (addr) T() ;
}

template <class T>     
void GenericAttributeLifeCycle<T>::virtual_destroy(
	Memory::pointer addr
	) {
		T* item = reinterpret_cast<T*>(addr) ;
		// Direct call of the destructor (see comments at the
		// beginning of this file).
		item-> ~T() ;
		//  Reset memory to zero, so that access to
		// pointers in destroyed object can be more easily detected.
		Memory::clear(addr, item_size()) ;
}

template <class T>     
void GenericAttributeLifeCycle<T>::virtual_copy(
	Memory::pointer to, Memory::pointer from 
	) {
		T* to_item = reinterpret_cast<T*>(to) ;
		T* from_item = reinterpret_cast<T*>(from) ;
		// calls operator= for type T.
		*to_item = *from_item ;
}

template <class T>     
void GenericAttributeLifeCycle<T>::virtual_copy_construct(
	Memory::pointer to, Memory::pointer from 
	) {
		//  Initialize memory to zero, so that uninitialized 
		// pointers can be more easily detected.
		Memory::clear(to, item_size()) ;
		T* from_item = reinterpret_cast<T*>(from) ;
		// Placement syntax of operator new (see comments at the
		// beginning of the file).
		new (to) T(*from_item) ;
}

template <class T>
AttributeLifeCycle* GenericAttributeLifeCycle<T>::clone() {
	return new GenericAttributeLifeCycle<T>(notify_, pod_) ;
}



#endif

