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



#ifndef _MAP_SERIALIZER_H_
#define _MAP_SERIALIZER_H_

#include "model_common.h"
#include "../model/map.h"
#include "../basic/counted.h"
#include "../basic/smart_pointer.h"
#include <iostream>


class AbstractMapBuilder;

class MODEL_API MapSerializer : public Counted
{
public:
	typedef Map::Vertex		Vertex ;
	typedef Map::Halfedge		Halfedge ;
	typedef Map::Facet			Facet ;

	typedef Map::Vertex_iterator	Vertex_iterator ;
	typedef Map::Halfedge_iterator	Halfedge_iterator ;
	typedef Map::Facet_iterator		Facet_iterator ;

public:
    MapSerializer()
		: read_supported_(false)
		, write_supported_(false)
	{} 

	virtual ~MapSerializer() {}

	virtual bool serialize_read(const std::string& file_name, Map* mesh) ;
	virtual bool serialize_write(const std::string& file_name, const Map* mesh) const ;

public:
	/**
	* checks whether the stream should be opened
	* in text or binary mode. Default returns false.
	*/
	virtual bool binary() const ;

	/**
	* checks whether reading and writing to streams is supported.
	*/
	virtual bool streams_supported() const ;

	/**
	* checks whether reading is implemented.
	*/
	virtual bool read_supported() const ;

	/**
	* checks whether writing is implemented.
	*/
	virtual bool write_supported() const ;

protected:
	virtual bool do_read(std::istream& in, AbstractMapBuilder& builder) ; 
	virtual bool do_write(std::ostream& out, const Map* mesh) const;

protected:
	bool read_supported_ ;
	bool write_supported_ ;

} ;

typedef SmartPointer<MapSerializer> MapSerializer_var ;


#endif

