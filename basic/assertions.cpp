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


#include "assertions.h"
#include "logger.h"

#include <cstdlib>


#ifndef WIN32
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#endif



static void ogf_abort() {
#ifdef WIN32
	std::abort() ;
#else
	int ogf_pid = getpid() ;
	Logger::err("Assert") << "Current pid is: " << ogf_pid << std::endl ;
	Logger::err("Assert") << "Going to bed" << std::endl ;
	Logger::err("Assert") << "Use: \'gdb MeshStudio " << ogf_pid << "\' and then \'where\' to see the stack trace" << std::endl ;

	kill(ogf_pid, SIGSTOP) ;

	/*
	if(!fork()) {
	char spid[1024] ;
	sprintf("%d", spid, ogf_pid) ;
	execl("gdb", "gdb", "MeshStudio", spid, nil) ;
	} else {
	kill(ogf_pid, SIGSTOP) ;
	}
	*/
#endif        
}

void ogf_assertion_failed(
						  const std::string& condition_string,
						  const std::string& file, int line
						  ) 
{
	Logger::err("Assert") << "Assertion failed: " << condition_string << std::endl ;
	Logger::err("Assert") << "File: " << file << std::endl ;
	Logger::err("Assert") << "Line: " << line << std::endl ;
	ogf_abort() ;
}

void ogf_range_assertion_failed(
								double value, double min_value, double max_value, 
								const std::string& file, int line
								) 
{
	Logger::err("Assert") << "Range assertion failed: " 
		<< value << " in " 
		<< "[ " << min_value << " ... " << max_value << " ]"
		<< std::endl ;
	Logger::err("Assert") << "File: " << file << std::endl ;
	Logger::err("Assert") << "Line: " << line << std::endl ;
	ogf_abort() ;
}

void ogf_should_not_have_reached(
								 const std::string& file, int line
								 ) 
{
	Logger::err("Assert") << "Control should not have reached this point:" << std::endl ;
	Logger::err("Assert") << "File: " << file << std::endl ;
	Logger::err("Assert") << "Line: " << line << std::endl ;
	ogf_abort() ;
}


