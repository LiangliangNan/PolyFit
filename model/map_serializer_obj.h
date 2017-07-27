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

#ifndef _MAP_SERIALIZER_OBJ_H_
#define _MAP_SERIALIZER_OBJ_H_

#include "model_common.h"
#include "map_serializer.h"
#include "../basic/color.h"


class MODEL_API MapSerializer_obj : public MapSerializer 
{
public:
	MapSerializer_obj();

	virtual bool serialize_read(
		const std::string& file_name, Map* mesh
		) ;

protected:
	virtual bool do_read(std::istream& input, AbstractMapBuilder& builder) ;        
	virtual bool do_write(std::ostream& output, const Map* mesh) const; 
	void read_mtl_lib(std::istream& input) ;

protected:
	std::string	current_directory_ ;

	typedef Color Material ;
	typedef std::map<std::string, Material> MaterialLib ;

	MaterialLib	material_lib_ ;
	Material	current_material_ ;
} ;


//_________________________________________________________

/**
* Extended obj file format, adds attributes.
*/
class MODEL_API MapSerializer_eobj : public MapSerializer_obj {
public:
	MapSerializer_eobj() ;

	virtual bool do_read(std::istream& input, AbstractMapBuilder& builder) ;        
	virtual bool do_write(std::ostream& output, const Map* mesh) const; 

} ;


#endif

