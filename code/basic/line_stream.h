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


#ifndef _LINE_STREAM_H_
#define _LINE_STREAM_H_

#include "basic_types.h"
#include "assertions.h"

#include <iostream>
#include <sstream>



namespace IO {

	//____________________________________________________________________________

	class LineInputStream {
	public:
		LineInputStream(std::istream& in) : in_(in), line_in_(nil) {   }
		~LineInputStream() {
			delete line_in_ ; 
			line_in_ = nil ;
		}

		bool eof() const { return in_.eof() ; }
		bool eol() const { return line_in_ == nil || line_in_->eof() ; }
		//bool ok() const { return in_ != nil; } // Liangliang: changed for vs2013
		bool ok() const { return !in_.fail(); }

		void get_line() {
			getline(in_, buffer_);
			//in_.getline(buffer_, 65536) ;
			delete line_in_ ; 
			line_in_ = new std::istringstream(buffer_) ;
		}

		std::istream& line() { 
			ogf_assert(line_in_ != nil) ;
			return *line_in_ ; 
		}

		const std::string& current_line() const {
			return buffer_;
		}

		template <class T> LineInputStream& operator>>(T& param) {
			*line_in_ >> param;
			return *this;
		}

	private:
		std::istream& in_ ;
		std::istringstream* line_in_ ;
		std::string buffer_;
		//char buffer_[65536] ;
	} ;
}


typedef IO::LineInputStream     LineInputStream;



#endif
