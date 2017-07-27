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


#include "attribute_manager.h"
#include "../basic/logger.h"
#include <typeinfo>


AttributeManager::~AttributeManager() {
	// Makes sure that some attributes are not still connected to
	// this AttributeManager. 
	// Displays an error message to help debugging.

	if(attributes_.size() != named_attributes_.size()) {
		Logger::err("AttributeManager") 
			<< "Fatal error: AttributeManager destroyed "
			<< "before one of its attributes"
			<< std::endl ;

		int nb_dangling_attrs = 
			attributes_.size() - named_attributes_.size() ;

		Logger::err("AttributeManager")
			<< "found " << nb_dangling_attrs << " dangling attribute(s)"
			<< std::endl ;
		Logger::err("AttributeManager")
			<< "Could come from a Map destroyed "
			<< "before an algorithm class plugged on it"
			<< std::endl ;
		Logger::err("AttributeManager")
			<< "To fix the problem: put the algorithm class in a scope { }"
			<< std::endl ;
		Logger::err("AttributeManager")
			<< "     or make the algorithm class unbind its attributes"
			<< " at the end of apply()"
			<< std::endl ;
		Logger::err("AttributeManager")
			<< "exiting..."
			<< std::endl ;
		ogf_assert(false) ;
	}
}

void AttributeManager::new_record(Record* record) {
	if(rat_.is_full()) {
		rat_.grow() ;
		for(std::set<AttributeStore*>::iterator 
			it=attributes_.begin(); it!=attributes_.end(); it++
			) {
				(*it)->grow() ;
				ogf_attribute_assert((*it)->capacity() == capacity()) ;
		}
	}
	record->set_record_id(rat_.new_record_id()) ;
	unsigned int chunk  = record->record_id().chunk() ;
	unsigned int offset = record->record_id().offset() ;
	for(std::set<AttributeStore*>::iterator 
		it=attributes_.begin(); it!=attributes_.end(); it++
		) {
			(*it)->construct(  (*it)->data(chunk,offset), record ) ;
	} 
	size_++ ;
}

void AttributeManager::new_record(Record* record, const Record* from) {
	if(rat_.is_full()) {
		rat_.grow() ;
		for(std::set<AttributeStore*>::iterator 
			it=attributes_.begin(); it!=attributes_.end(); it++
			) {
				(*it)->grow() ;
				ogf_attribute_assert((*it)->capacity() == capacity()) ;
		}
	}
	record->set_record_id(rat_.new_record_id()) ;
	unsigned int chunk  = record->record_id().chunk() ;
	unsigned int offset = record->record_id().offset() ;

	unsigned int chunk_from  = from->record_id().chunk() ;
	unsigned int offset_from = from->record_id().offset() ; 

	for(std::set<AttributeStore*>::iterator 
		it=attributes_.begin(); it!=attributes_.end(); it++
		) {
			(*it)->copy_construct(  
				(*it)->data(chunk,offset), record,
				(*it)->data(chunk_from,offset_from), from
				) ;
	} 
	size_++ ;
}

void AttributeManager::copy_record(Record* record, const Record* from) {
	unsigned int chunk  = record->record_id().chunk() ;
	unsigned int offset = record->record_id().offset() ;

	unsigned int chunk_from  = from->record_id().chunk() ;
	unsigned int offset_from = from->record_id().offset() ; 

	for(std::set<AttributeStore*>::iterator 
		it=attributes_.begin(); it!=attributes_.end(); it++
		) {
			(*it)->copy(  
				(*it)->data(chunk,offset), record,
				(*it)->data(chunk_from,offset_from), from
				) ;
	} 
}

void AttributeManager::delete_record(Record* record) {
	unsigned int chunk  = record->record_id().chunk() ;
	unsigned int offset = record->record_id().offset() ;
	for(std::set<AttributeStore*>::iterator 
		it=attributes_.begin(); it!=attributes_.end(); it++
		) {
			(*it)->destroy(  (*it)->data(chunk,offset), record ) ;
	}    
	rat_.delete_record_id(record->record_id()) ;
	record->record_id().forget() ;
	size_-- ;
}

void AttributeManager::clear() {
	for(std::set<AttributeStore*>::iterator 
		it=attributes_.begin(); it!=attributes_.end(); it++
		) {
			for(unsigned int chunk=0; chunk<rat_.nb_chunks(); chunk++) {
				for(unsigned int offset=0; offset<RAT::CHUNK_SIZE; offset++) {
					if(!rat_.cell(chunk,offset).is_free()) {
						(*it)->destroy( (*it)->data(chunk,offset) ) ;
					} 
				}
			}
			(*it)->clear() ;
	} 
	rat_.clear() ;
	size_ = 0 ;
}

void AttributeManager::list_named_attributes(
	std::vector<std::string>& names
	) {
		names.clear() ;
		for(
			std::map<std::string, AttributeStore_var>::iterator
			it=named_attributes_.begin();
		it!=named_attributes_.end(); it++
			) {
				names.push_back(it->first) ;
		}
}

bool AttributeManager::named_attribute_is_bound(
	const std::string& name
	) {
		return (
			named_attributes_.find(name) != 
			named_attributes_.end()
			) ;
}

void AttributeManager::register_attribute_store(AttributeStore* as) {
	ogf_assert(
		attributes_.find(as) == attributes_.end() 
		) ;
	attributes_.insert(as) ;

	for(unsigned int chunk=0; chunk<rat_.nb_chunks(); chunk++) {
		for(unsigned int offset=0; offset<RAT::CHUNK_SIZE; offset++) {
			if(!rat_.cell(chunk,offset).is_free()) {
				as->construct( as->data(chunk,offset) ) ;
			}
		}
	}
}

void AttributeManager::unregister_attribute_store(AttributeStore* as) {

	for(unsigned int chunk=0; chunk<rat_.nb_chunks(); chunk++) {
		for(unsigned int offset=0; offset<RAT::CHUNK_SIZE; offset++) {
			if(!rat_.cell(chunk,offset).is_free()) {
				as->destroy( as->data(chunk,offset) ) ;
			} 
		}
	}

	std::set<AttributeStore*>::iterator it = attributes_.find(as) ;
	ogf_assert(it != attributes_.end()) ;
	attributes_.erase(it) ;
}

void AttributeManager::bind_named_attribute_store(
	const std::string& name, AttributeStore* as
	) {
		ogf_assert( !named_attribute_is_bound(name) ) ;
		named_attributes_[name] = as ;
}

AttributeStore* AttributeManager::resolve_named_attribute_store(
	const std::string& name
	) {
		std::map<std::string, AttributeStore_var>::iterator 
			it=named_attributes_.find(name) ;
		ogf_assert(it != named_attributes_.end()) ;
		return it->second ;
}

void AttributeManager::delete_named_attribute(
	const std::string& name
	) {
		std::map<std::string, AttributeStore_var>::iterator 
			it=named_attributes_.find(name) ;
		ogf_assert(it != named_attributes_.end()) ;
		ogf_assert(!it->second->is_shared()) ;
		named_attributes_.erase(it) ;
}

const std::type_info& AttributeManager::resolve_named_attribute_type_id( const std::string& name ){
	AttributeStore* store = resolve_named_attribute_store( name );
	return store->attribute_type_id();
}

