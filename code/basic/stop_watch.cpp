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

#include "stop_watch.h"
#include "basic_types.h" // for "round()"
#include <iostream>



//_________________________________________________________

StopWatch::StopWatch() {
	start();
}

StopWatch::~StopWatch() {}


void StopWatch::start() {
#ifdef WIN32
	LARGE_INTEGER  largeInteger;
	QueryPerformanceFrequency(&largeInteger);
	freq_ = largeInteger.QuadPart;
	QueryPerformanceCounter(&largeInteger);
	start_count_ = largeInteger.QuadPart;
#else
	gettimeofday(&start_time_, 0);
#endif // WIN32
}

double StopWatch::elapsed() const {
#ifdef WIN32
	LARGE_INTEGER  largeInteger;
	QueryPerformanceCounter(&largeInteger);
	LONGLONG now_count = largeInteger.QuadPart;
	double time = (double)( (now_count - start_count_) / static_cast<double>(freq_) );
#else
	timeval now;
	gettimeofday(&now, 0);
    double time = (now.tv_sec - start_time_.tv_sec) + (now.tv_usec - start_time_.tv_usec) / 1.0e6;
#endif  // WIN32
    return truncate_digits(time, 2);
}
