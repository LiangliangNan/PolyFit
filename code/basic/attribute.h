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


#ifndef ___ATTRIBUTE__
#define ___ATTRIBUTE__


#include "attribute_manager.h"
#include "attribute_store.h"
#include "attribute_life_cycle.h"


class /*BASIC_API*/ AttributeBase {
protected:
    static AttributeStore* resolve_named_attribute_store(
        AttributeManager* manager, const std::string& name
    ) {
        return manager->resolve_named_attribute_store(name) ;
    }

    static void bind_named_attribute_store(
        AttributeManager* manager, 
        const std::string& name, AttributeStore* as
    ) {
        manager->bind_named_attribute_store(name,as) ;
    }
} ;

/**
 * enables access and (creation or retreival) to the attributes
 * stored in an AttributeManager. Note that distinct Attributes
 * can refer to the same attribute (they can share the same
 * AttributeStore).
 */

template <class RECORD, class ATTRIBUTE> 
class Attribute : public AttributeBase {
public:

    typedef RECORD Record;
    typedef GenericAttributeManager<RECORD> AttributeManager;
    typedef GenericAttributeStore<ATTRIBUTE> AttributeStore;
    typedef GenericAttributeLifeCycle<ATTRIBUTE> DefaultAttributeLifeCycle;
    typedef Attribute<RECORD,ATTRIBUTE> thisclass;
    
    typedef ::AttributeManager::Mode Mode ;

    /**
     * constructs an unbound Attribute.
     */
    Attribute() { } 

    /**
     * constructs a transient Attribute. If specified, the 
     * optional AttributeLifeCycle is used to manage attribute 
     * creation and destruction.
     */
    Attribute(
        AttributeManager* manager, AttributeLifeCycle* life_cycle=nil
    ) ;

    /**
     * constructs or retreives a persistent Attribute (according
     * to the specified mode). mode is one of FIND, CREATE, 
     * FIND_OR_CREATE (default). Note that two attributes of
     * different types cannot have the same name (the type
     * is dynamically detected by the AttributeManager, using
     * RTTI). If specified, the optional AttributeLifeCycle is
     * used to manage attribute creation and destruction.
     */
    Attribute(
        AttributeManager* manager, const std::string& name, 
        Mode mode=AttributeManager::FIND_OR_CREATE, 
        AttributeLifeCycle* life_cycle=nil
    ) ;
    
    /**
     * constructs an Attribute refering to the same attribute
     * as rhs.
     */
    Attribute(const thisclass& rhs):store_(rhs.store_) { }

    /**
     * makes this Attribute refer to the same attribute as rhs.
     */
    thisclass& operator=(const thisclass& rhs) { 
        store_=rhs.store_; return *this;
    }

    bool is_bound() const { return !store_.is_nil(); }

    AttributeManager* attribute_manager() const {
        ogf_assert(is_bound()) ;
        AttributeManager* result = dynamic_cast<AttributeManager*>(
            store_->attribute_manager()
        ) ;
        ogf_assert(result != nil) ;
        return result ;
    }
    
    void bind(
        AttributeManager* manager, AttributeLifeCycle* life_cycle=nil
    ) ;

    void bind(
        AttributeManager* manager, const std::string& name, 
        Mode mode=AttributeManager::FIND_OR_CREATE, 
        AttributeLifeCycle* life_cycle=nil
    ) ;
    
    bool bind_if_defined(
        AttributeManager* manager, const std::string& name
    ) ;

    void unbind() ;

    ATTRIBUTE& operator[](const RECORD& record) {
        return *data(record) ;
    }

    const ATTRIBUTE& operator[](const RECORD& record) const {
        return *data(record) ;
    }
    
    ATTRIBUTE& operator[](const RECORD* record) {
        return *data(*record) ;
    }

    const ATTRIBUTE& operator[](const RECORD* record) const {
        return *data(*record) ;
    }

    /**
     * Checks whether manager has an attribute of this type
     * bound with the specified name.
     */
    static bool is_defined(
        AttributeManager* manager, const std::string& name
    ) ; 

	// Liangliang: for each access of type id
	virtual const std::type_info& attribute_type_id() const {
		return typeid(ATTRIBUTE) ;
	}

protected:
    ATTRIBUTE* data(const RECORD& record) const {
        return reinterpret_cast<ATTRIBUTE*>(store_->data(record)) ;
    }

private:
    AttributeStore_var store_ ;
} ;


//__________________________________________________________

template <class RECORD, class ATTRIBUTE> inline 
void Attribute<RECORD,ATTRIBUTE>::bind(
    AttributeManager* manager, AttributeLifeCycle* life_cycle
) {
    if(life_cycle == nil) {
        life_cycle = new DefaultAttributeLifeCycle ;
    }
    store_ = new AttributeStore(life_cycle, manager);
}

template <class RECORD, class ATTRIBUTE> inline 
bool Attribute<RECORD,ATTRIBUTE>::is_defined(
    AttributeManager* manager, const std::string& name
) {
    return (
        manager->named_attribute_is_bound(name) &&
        dynamic_cast<AttributeStore*>(
            resolve_named_attribute_store(manager,name)
        ) != nil
    ) ;
} 

template <class RECORD, class ATTRIBUTE> inline 
void Attribute<RECORD,ATTRIBUTE>::bind(
    AttributeManager* manager, const std::string& name, 
    Mode mode, AttributeLifeCycle* life_cycle
) {
    switch(mode) {
    case AttributeManager::FIND:
    {
        ogf_assert(manager->named_attribute_is_bound(name)) ;
        ogf_assert(life_cycle == nil) ;
        store_ = resolve_named_attribute_store(manager,name) ;
        ogf_assert(store_ != nil) ;
        // Sanity check, checks the attribute type.
        ::AttributeStore* check_type = store_ ;
        ogf_assert(dynamic_cast<AttributeStore*>(check_type) != nil) ;
        break ;
    }
    case AttributeManager::CREATE:
    {
        ogf_assert(!manager->named_attribute_is_bound(name)) ;
        if(life_cycle == nil) {
            life_cycle = new DefaultAttributeLifeCycle ;
        }
        store_ = new AttributeStore(life_cycle,manager);
        bind_named_attribute_store(manager,name,store_) ;
        break ;
    }
    case AttributeManager::FIND_OR_CREATE:
    {
        if(manager->named_attribute_is_bound(name)) {
            store_ = resolve_named_attribute_store(manager,name) ;
            ogf_assert(store_ != nil) ;
            // Sanity check, checks the attribute type.
            ::AttributeStore* check_type = store_ ;
            ogf_assert(dynamic_cast<AttributeStore*>(check_type) != nil) ;
        } else {
            if(life_cycle == nil) {
                life_cycle = new DefaultAttributeLifeCycle ;
            }
            store_ = new AttributeStore(life_cycle,manager) ;
            bind_named_attribute_store(manager,name,store_) ;
        }
        break ;
    }
    }
}

template <class RECORD, class ATTRIBUTE> inline 
void Attribute<RECORD,ATTRIBUTE>::unbind() {
    store_.forget() ;
}

template <class RECORD, class ATTRIBUTE> inline 
Attribute<RECORD,ATTRIBUTE>::Attribute(
    AttributeManager* manager, AttributeLifeCycle* life_cycle
) {
    bind(manager,life_cycle) ;
}

template <class RECORD, class ATTRIBUTE> inline 
Attribute<RECORD,ATTRIBUTE>::Attribute(
    AttributeManager* manager, const std::string& name, 
    Mode mode, AttributeLifeCycle* life_cycle
) {
    bind(manager,name,mode,life_cycle) ;
}

template <class RECORD, class ATTRIBUTE> inline 
bool Attribute<RECORD,ATTRIBUTE>::bind_if_defined(
    AttributeManager* manager,
    const std::string& name
) {
    if(is_defined(manager,name)) {
        bind(manager, name) ;
        return true ;
    }
    return false ;
}


#endif
