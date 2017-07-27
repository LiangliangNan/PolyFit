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

#ifndef ___ATTRIBUTE_SERIALIZER__
#define ___ATTRIBUTE_SERIALIZER__

#include "basic_common.h"
#include "attribute.h"

#include <map>
#include <string>
#include <typeinfo>


/**
* AttributeSerializer is used to save and load attributes attached
* to an object. This is the base class for serializing the value of 
* an attribute and creating an AttributeStore given its type name.
*/
class BASIC_API AttributeSerializer : public Counted {
public:
	static void initialize() ;
	static void terminate() ;

	static AttributeSerializer* resolve_by_type(const std::type_info& attribute_type) ;
	static AttributeSerializer* resolve_by_name(const std::string& attribute_type_name) ;
	static std::string find_name_by_type(const std::type_info& attribute_type) ;
	static void register_type_alias(const std::string& alias, const std::string& name) {
		(*alias_to_name_)[alias] = name ;
	}

	static void bind(
		const std::type_info& attribute_type,
		const std::string& attribute_type_name,
		AttributeSerializer* serializer
		) ;

	virtual AttributeStore* create_attribute_store(AttributeManager* manager) = 0 ;
	virtual bool serialize_read(std::istream& in,   Memory::pointer addr) = 0 ;
	virtual bool serialize_write(std::ostream& out, Memory::pointer addr) = 0 ;

private:
	typedef std::map<std::string, SmartPointer<AttributeSerializer> > SerializerMap ;
	typedef std::map<std::string, std::string> StringMap ;
	static SerializerMap* type_to_serializer_ ;
	static SerializerMap* name_to_serializer_ ;
	static StringMap*     alias_to_name_ ;
	static StringMap*     type_to_name_ ;
} ;

typedef SmartPointer<AttributeSerializer> AttributeSerializer_var ;

//_________________________________________________________________________________________________________

/**
* Default implementation of AttributeSerializer.
*/
template <class ATTRIBUTE> class GenericAttributeSerializer : public AttributeSerializer {
public:
	typedef GenericAttributeStore<ATTRIBUTE> AttributeStore;
	typedef GenericAttributeLifeCycle<ATTRIBUTE> DefaultAttributeLifeCycle;

	virtual AttributeStore* create_attribute_store(AttributeManager* manager) {
		DefaultAttributeLifeCycle* life_cycle = new DefaultAttributeLifeCycle ;
		return new AttributeStore(life_cycle, manager) ;
	}
	virtual bool serialize_read(std::istream& in, Memory::pointer addr) {
		ATTRIBUTE& attr = *reinterpret_cast<ATTRIBUTE*>(addr) ;
		in >> attr ;
		return true ;
	}
	virtual bool serialize_write(std::ostream& out, Memory::pointer addr) {
		ATTRIBUTE& attr = *reinterpret_cast<ATTRIBUTE*>(addr) ;
		out << attr ;
		return true ;
	}

} ;

//_________________________________________________________________________________________________________

/**
* Use this class to declare a new serializable attribute type.
* In the common.cpp file of the library, add:
* ogf_register_attribute_type<MyAttributeType>("MyAttributeType") ;
*/
template <class T> 
class register_attribute_type {
public:
	register_attribute_type(const std::string& type_name) {
		AttributeSerializer::bind(typeid(T), type_name, new GenericAttributeSerializer<T>()) ;
	}
} ;

//_________________________________________________________________________________________________________

inline void register_attribute_type_alias(const std::string& alias, const std::string& type_name) {
	AttributeSerializer::register_type_alias(alias, type_name) ;
}

//_________________________________________________________________________________________________________


/**
* SerializedAttributeRef is what SerializedAttribute::operator[] returns.
* It is just meant to overload operator<< and operator>>.
*/
class SerializedAttributeRef {
public:
	SerializedAttributeRef(
		AttributeSerializer* ser, Memory::pointer addr
		) : serializer_(ser), addr_(addr) {
	}
	AttributeSerializer* serializer() const { return serializer_ ; }
	Memory::pointer addr() const { return addr_ ; }
private:
	AttributeSerializer* serializer_ ;
	Memory::pointer addr_ ;
} ;

inline std::istream& operator>>(std::istream& in, const SerializedAttributeRef& r) {
	r.serializer()->serialize_read(in, r.addr()) ;
	return in ;
}

inline std::ostream& operator<<(std::ostream& out, const SerializedAttributeRef& r) {
	r.serializer()->serialize_write(out, r.addr()) ;
	return out ;
}

//_________________________________________________________________________________________________________

/**
* SerializedAttribute allows writing attribute values into a stream,
* reading attribute values from a stream, and creating an attribute
* from its type name. 
*/
template <class RECORD> class SerializedAttribute : public AttributeBase {
public:
	typedef GenericAttributeManager<RECORD> AttributeManager;

	void bind(AttributeManager* manager, const std::string& name) {
		attribute_manager_ = manager ;
		attribute_store_ = resolve_named_attribute_store(manager, name) ;
		if(attribute_store_ != nil) {
			serializer_ = AttributeSerializer::resolve_by_type(attribute_store_->attribute_type_id()) ;
		}
		name_ = name ;
	}

	void bind(AttributeManager* manager, const std::string& name, const std::string& attribute_type_name) {
		attribute_manager_ = manager ;
		serializer_ = AttributeSerializer::resolve_by_name(attribute_type_name) ;
		if(serializer_ != nil) {
			if(attribute_manager_->named_attribute_is_bound(name)) {
				attribute_store_ = resolve_named_attribute_store(attribute_manager_,name) ;
				ogf_assert(
					AttributeSerializer::find_name_by_type(
					attribute_store_->attribute_type_id()
					) == attribute_type_name
					) ;
			} else {
				attribute_store_ = serializer_->create_attribute_store(attribute_manager_) ;
				bind_named_attribute_store(attribute_manager_,name,attribute_store_) ;
			}
		}
		name_ = name ;
	}

	void unbind() {
		attribute_manager_ = nil ;
		attribute_store_ = nil ;
		serializer_ = nil ;
	}

	SerializedAttribute() {
		attribute_manager_ = nil ;
		attribute_store_ = nil ;
		serializer_ = nil ;
	}

	SerializedAttribute(AttributeManager* manager, const std::string& name) {
		bind(manager, name) ;
	}

	SerializedAttribute(
		AttributeManager* manager, const std::string& name, const std::string& attribute_type_name
		) {
			bind(manager, name, attribute_type_name) ;
	}

	SerializedAttribute(const SerializedAttribute& rhs) {
		attribute_manager_ = rhs.attribute_manager_ ;
		attribute_store_   = rhs.attribute_store_ ;
		serializer_        = rhs.serializer_ ;
		name_              = rhs.name_ ;
	}

	bool is_bound() const {
		return (attribute_manager_ != nil) && (attribute_store_ != nil) && (serializer_ != nil) ;
	}

	const std::string& name() const { return name_ ; }

	std::string type_name() const {
		ogf_assert(attribute_store_ != nil) ;
		return AttributeSerializer::find_name_by_type(attribute_store_->attribute_type_id()) ;
	}

	SerializedAttributeRef operator[](const RECORD* record) {
		return SerializedAttributeRef(serializer_, attribute_store_->data(*record)) ;
	}

private:
	AttributeManager* attribute_manager_ ;
	AttributeStore* attribute_store_ ;
	AttributeSerializer* serializer_ ;
	std::string name_ ;
} ;



#endif
