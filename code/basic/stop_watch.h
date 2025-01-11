/* ---------------------------------------------------------------------------
 * Copyright (C) 2017 Liangliang Nan <liangliang.nan@gmail.com>
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of PolyFit. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 *
 *     Liangliang Nan and Peter Wonka.
 *     PolyFit: Polygonal Surface Reconstruction from Point Clouds.
 *     ICCV 2017.
 *
 *  For more information:
 *  https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html
 * ---------------------------------------------------------------------------
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

