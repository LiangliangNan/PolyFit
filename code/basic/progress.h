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

#ifndef __BASIC_PROGRESS__
#define __BASIC_PROGRESS__

#include "basic_common.h"
#include "basic_types.h"


class ProgressClient ;

/**
* For internal use, client code should use UserProgress.
*/
class BASIC_API Progress {
public:
	Progress() ;
	virtual ~Progress() ;

	static Progress* instance() ;

	virtual void notify(std::size_t new_val) ;

	void set_client(ProgressClient* c) { client_ = c ; }

	void push() ;
	void pop() ;

	void cancel()            { canceled_ = true ;  }
	void clear_canceled()    { canceled_ = false ; }
	bool is_canceled() const { return canceled_ ;  }
private:
	static Progress* instance_ ;
	ProgressClient* client_ ;
	int  level_ ;
	bool canceled_ ;
} ;

//_________________________________________________________

/**
* For internal use, client code do not need to use this one.
*/
class BASIC_API ProgressClient {
public:
	virtual void notify_progress(std::size_t new_val) = 0;
	virtual ~ProgressClient() ;
} ;

//_________________________________________________________

class BASIC_API ProgressLogger {
public:
	ProgressLogger(std::size_t max_val = 100, const std::string& task_name = "", bool quiet = false) ;
	virtual ~ProgressLogger() ;

	virtual void notify(std::size_t new_val) ;
	virtual void next() ;
	bool is_canceled() const {
		return Progress::instance()->is_canceled() ;
	}
	void reset() { notify(0) ; }
	void reset(std::size_t max_val) ;

protected:
	virtual void update() ;

private:
	std::size_t max_val_ ;
	std::string task_name_ ;
	std::size_t cur_val_ ;
	std::size_t cur_percent_ ;
	bool quiet_ ;
} ;


#endif

