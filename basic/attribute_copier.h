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

#ifndef ___ATTRIBUTE_COPIER__
#define ___ATTRIBUTE_COPIER__


#include "attribute.h"
#include "../basic/logger.h"

#include <string>
#include <typeinfo>
#include <vector>


template <class RECORD> class AttributeCopier {
public:
	AttributeCopier() { 
		source_manager_ = nil ; source_store_ = nil ;
		dest_manager_ = nil ; dest_store_ = nil ;
	}
	bool source_is_bound() const { return source_store_ != nil ; }
	void bind_source(AttributeManager* manager, const std::string& name) {
		ogf_parano_assert(manager->record_type_id() == typeid(RECORD)) ;
		source_manager_ = manager ;
		if(source_manager_->named_attribute_is_bound(name)) {
			source_store_ = source_manager_->resolve_named_attribute_store(name) ;
		} else {
			source_store_ = nil ;
		}
		source_name_ = name ;
	}
	bool destination_is_bound() const { return dest_store_ != nil ; }
	void bind_destination(AttributeManager* manager, const std::string& name) {
		ogf_assert(source_is_bound()) ;
		dest_manager_ = manager ;
		if(dest_manager_->named_attribute_is_bound(name)) {
			dest_store_ = dest_manager_->resolve_named_attribute_store(name) ;
		} else {
			dest_store_ = source_store_->clone() ;
			dest_store_->bind(dest_manager_) ;
			dest_manager_->bind_named_attribute_store(name, dest_store_) ;
		}
		if(source_store_->attribute_type_id() != dest_store_->attribute_type_id()) {
			Logger::err("AttributeCopier") << "Attribute " << name 
				<< " does not have same type as source attribute"
				<< std::endl ;
			unbind_destination() ;
		}
		dest_name_ = name ;
	}
	void bind_destination(AttributeManager* manager) {
		bind_destination(manager, source_name_) ;
	}
	void unbind_source() {
		source_manager_ = nil ;
		source_store_ = nil ;
		source_name_ = "" ;
	}
	void unbind_destination() {
		dest_manager_ = nil ;
		dest_store_ = nil ;
		dest_name_ = "" ; 
	}
	AttributeManager* source_attribute_manager() const { return source_manager_ ; }
	const std::string& source_name() const { return source_name_ ; }
	AttributeManager* destination_attribute_manager() const { return dest_manager_ ; }
	const std::string& destination_name() const { return dest_name_ ; }
	void copy_attribute(RECORD* to, RECORD* from) {
		ogf_assert(source_is_bound()) ;
		ogf_assert(destination_is_bound()) ;
		Memory::pointer to_addr = dest_store_->data(*to) ;
		Memory::pointer from_addr = source_store_->data(*from) ;
		source_store_->copy(to_addr, to, from_addr, from) ;
	}
private:
	AttributeManager* source_manager_ ;
	AttributeStore* source_store_ ;
	AttributeManager* dest_manager_ ;
	AttributeStore* dest_store_ ;
	std::string source_name_ ;
	std::string dest_name_ ;
} ;

//_________________________________________________________________________________________________________

template <class RECORD> inline void copy_attributes(
	std::vector<AttributeCopier<RECORD> >& copiers,
	RECORD* to, RECORD* from
	) {
		for(unsigned int i=0; i<copiers.size(); i++) {
			if(copiers[i].source_is_bound() && copiers[i].destination_is_bound()) {
				copiers[i].copy_attribute(to, from) ;
			}
		}
}

template <class RECORD> inline void bind_source(
	std::vector<AttributeCopier<RECORD> >& copiers,
	AttributeManager* source_attribute_manager,
	const std::string& name
	) {
		for(unsigned int i=0; i<copiers.size(); i++) {
			if(copiers[i].source_name() == name) {
				if(copiers[i].source_attribute_manager() == source_attribute_manager) {
					return ;
				}
				copiers[i].bind_source(source_attribute_manager, name) ;
				return ;
			}
		}
		copiers.push_back(AttributeCopier<RECORD>()) ;
		copiers.rbegin()->bind_source(source_attribute_manager, name) ;
}

template <class RECORD> inline void bind_destinations(
	std::vector<AttributeCopier<RECORD> >& copiers,
	AttributeManager* destination_attribute_manager
	) {
		for(unsigned int i=0; i<copiers.size(); i++) {
			if(copiers[i].source_is_bound()) {
				if(!copiers[i].destination_is_bound() ||
					copiers[i].destination_name() != copiers[i].source_name() ||
					copiers[i].destination_attribute_manager() != destination_attribute_manager
					) {
						copiers[i].bind_destination(destination_attribute_manager, copiers[i].source_name()) ;
				}
			}
		}
}


#endif
