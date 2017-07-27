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


#ifndef _BASIC_GENERIC_ATTRIBUTES_IO_H_
#define _BASIC_GENERIC_ATTRIBUTES_IO_H_

#include "basic_common.h"
#include "attribute_manager.h"
#include "attribute_serializer.h"
#include <iostream>


template <class T> 
inline bool get_serializable_attributes(
	GenericAttributeManager<T>* manager, std::vector<SerializedAttribute<T> >& attributes,
	std::ostream& out, const std::string& location, const std::string& attribute_kw = "# attribute"
	) {
		bool result = false ;
		std::vector<std::string> names ;
		manager->list_named_attributes(names) ;
		for(unsigned int i=0; i<names.size(); i++) {
			attributes.push_back(SerializedAttribute<T>()) ;
			attributes.rbegin()->bind(manager, names[i]) ;
			if(attributes.rbegin()->is_bound()) {
				std::cerr << "Attribute " << names[i] << " on " << location << " : " 
					<< attributes.rbegin()->type_name() << std::endl ;
				out << attribute_kw << " " << names[i] << " " << location << " " 
					<< attributes.rbegin()->type_name() << std::endl ;
				result = true ;
			} else {
				std::cerr << "Attribute " << names[i] << " on " << location 
					<< " is not serializable" << std::endl ;
				attributes.pop_back() ;
			}
		}
		return result ;
}

template <class T> 
inline void serialize_read_attributes(
	std::istream& in, const T* item, std::vector<SerializedAttribute<T> >& attributes
	) {
		for(unsigned int i=0; i<attributes.size(); i++) {
			in >> attributes[i][item] ;
		}
}

template <class T> 
inline void serialize_write_attributes(
	std::ostream& out, const T* item, std::vector<SerializedAttribute<T> >& attributes
	) {
		for(unsigned int i=0; i<attributes.size(); i++) {
			out << attributes[i][item] << " " ;
		}
}


#endif

