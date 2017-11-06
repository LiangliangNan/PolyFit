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


#include "progress.h"
#include "assertions.h"


Progress* Progress::instance_ = nil ;

Progress::Progress() : client_(nil), level_(0), canceled_(false) {
}

Progress::~Progress() {
}

void Progress::push() { 
	level_++ ; 
	if(level_ == 1) {
		clear_canceled() ;
	}
}

void Progress::pop() {
	ogf_assert(level_ > 0) ;
	level_-- ;
}

void Progress::notify(std::size_t new_val) {
	if(client_ != nil && level_ < 2) {
		client_->notify_progress(new_val) ;
	}
}

Progress* Progress::instance() {
	if(instance_ == nil) {
		instance_ = new Progress ;
	}
	return instance_ ;
}


//_________________________________________________________


ProgressClient::~ProgressClient() {
	if (Progress::instance())
		delete Progress::instance();
}

//_________________________________________________________

ProgressLogger::ProgressLogger(std::size_t max_val, const std::string& task_name, bool quiet)
: max_val_(max_val), task_name_(task_name), quiet_(quiet) 
{
	cur_val_ = 0 ; 
	cur_percent_ = 0 ;
	Progress::instance()->push() ;
	if(!quiet_) {
		Progress::instance()->notify(0) ;
	}
}

void ProgressLogger::reset(std::size_t max_val)  {
	max_val_ = max_val ;
	reset() ;
}

ProgressLogger::~ProgressLogger() {
    Progress::instance()->notify(100) ;
    Progress::instance()->notify(0) ;
    Progress::instance()->pop() ;
}

void ProgressLogger::next() {
	cur_val_++ ;
	update() ;
}

void ProgressLogger::notify(std::size_t new_val) {
	cur_val_ = new_val ;
	update() ;
}


void ProgressLogger::update() {
	std::size_t percent = cur_val_ * 100 / ogf_max<std::size_t>(1, max_val_-1) ;
	
	if(percent != cur_percent_) {
		cur_percent_ = percent ;

		if(!quiet_) {
			Progress::instance()->notify(ogf_min<std::size_t>(cur_percent_, 100)) ;
		}
	}
}

