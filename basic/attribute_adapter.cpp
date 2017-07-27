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


#include "attribute_adapter.h"


AttributeAdapterBase::SecondaryType AttributeAdapterBase::parse_name(std::string& name) {
	int point_pos = name.find('.')  ;
	if(point_pos < 0) {
		return ATTR_VALUE ;
	}
	std::string extension = name.substr(point_pos+1, name.length() - point_pos - 1) ;
	name = name.substr(0, point_pos) ;
	if(extension == "x" || extension == "real") {
		return ATTR_X ;
	} else if(extension == "y" || extension == "imag" || extension == "imaginary") {
		return ATTR_Y ;
	} else if(extension == "z") {
		return ATTR_Z ;
	} else if(extension == "norm" || extension == "modulus" || extension == "length") {
		return ATTR_NORM ;
	} 
	return ATTR_ERROR ;
}
