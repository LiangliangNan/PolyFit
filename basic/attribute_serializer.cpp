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



#include "attribute_serializer.h"



AttributeSerializer::SerializerMap* AttributeSerializer::type_to_serializer_ = nil ;
AttributeSerializer::SerializerMap* AttributeSerializer::name_to_serializer_  = nil ;
AttributeSerializer::StringMap*     AttributeSerializer::alias_to_name_       = nil ;    
AttributeSerializer::StringMap*     AttributeSerializer::type_to_name_        = nil ;    

void AttributeSerializer::initialize() {
	ogf_assert(type_to_serializer_ == nil) ;
	type_to_serializer_ = new SerializerMap ;
	ogf_assert(name_to_serializer_ == nil) ;
	name_to_serializer_ = new SerializerMap ;
	ogf_assert(alias_to_name_ == nil) ;
	alias_to_name_ = new StringMap ;
	ogf_assert(type_to_name_ == nil) ;
	type_to_name_ = new StringMap ;
}

void AttributeSerializer::terminate() {
	delete type_to_serializer_ ;
	type_to_serializer_ = nil ;
	delete name_to_serializer_ ;
	name_to_serializer_ = nil ;
	delete alias_to_name_ ;
	alias_to_name_ = nil ;
	delete type_to_name_ ;
	type_to_name_ = nil ;
}

AttributeSerializer* AttributeSerializer::resolve_by_type(const std::type_info& attribute_type) {
	ogf_assert(type_to_serializer_ != nil) ;
	SerializerMap::iterator it = type_to_serializer_->find(attribute_type.name()) ;
	if(it == type_to_serializer_->end()) {
		return nil ;
	}
	return it->second ;
}

AttributeSerializer* AttributeSerializer::resolve_by_name(const std::string& type_name_in) {
	ogf_assert(alias_to_name_ != nil) ;
	std::string type_name = type_name_in ;
	{
		StringMap::const_iterator it = alias_to_name_->find(type_name_in) ;
		if(it != alias_to_name_->end()) {
			type_name = it->second ;
		}
	}
	ogf_assert(name_to_serializer_ != nil) ;       
	SerializerMap::iterator it = name_to_serializer_->find(type_name) ;
	if(it == name_to_serializer_->end()) {
		return nil ;
	}
	return it->second ;
}


std::string AttributeSerializer::find_name_by_type(const std::type_info& attribute_type) {
	ogf_assert(type_to_name_ != nil) ; 
	StringMap::iterator it = type_to_name_->find(attribute_type.name()) ;
	if(it == type_to_name_->end()) {
		return "unknown" ;
	}
	return it->second ;
}

void AttributeSerializer::bind(
							   const std::type_info& attribute_type, const std::string& attribute_type_name, 
							   AttributeSerializer* serializer
							   ) 
{
	ogf_assert(resolve_by_type(attribute_type) == nil) ;
	ogf_assert(resolve_by_name(attribute_type_name) == nil) ;
	(*type_to_serializer_)[attribute_type.name()] = serializer ;
	(*name_to_serializer_)[attribute_type_name] = serializer ;
	(*type_to_name_)[attribute_type.name()] = attribute_type_name ;
}
