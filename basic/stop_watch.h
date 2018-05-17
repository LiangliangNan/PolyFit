/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/ - liangliang.nan@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef _STOP_WATCH_H_
#define _STOP_WATCH_H_

#include "basic_common.h"


#ifdef WIN32
#	include <windows.h>
#else 
#	include <sys/time.h>
#endif // WIN32

//______________________________________________________________________


/**
* 
* The purpose of this file is to make a timer function
* that is as precise as possible on any given platform.
*
* usage example:
*   {
*      StopWatch w ;
*      // do task_1 ...
*      std::cout << "task_1 done. time: " << w.elapsed() << " sec.";
*	   w.start();
*      // do task_2 ...
*      std::cout << "task_2 done. time: " << w.elapsed() << " sec.";
*   } 
*/

class BASIC_API StopWatch 
{
public :
	StopWatch() ; // the watch will automatically start in construction
	~StopWatch() ;

	void  start() ;

	// returns user elapsed time since the construction / start in sec.
	double elapsed() const ;

private:

#ifdef WIN32
	LONGLONG  freq_;
	LONGLONG  start_count_;
#else
	timeval start_time_;
#endif

} ;


#endif

