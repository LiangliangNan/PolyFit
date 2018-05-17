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

#ifndef _MATH_GEOMETRY_LINE_H_
#define _MATH_GEOMETRY_LINE_H_

#include "math_common.h"
#include "vecg.h"



template <int DIM, class FT> 
class GenericLine {
public:
	typedef vecng<DIM, FT>		 Point ;
	typedef vecng<DIM, FT>		 Vector ;
	typedef GenericLine<DIM, FT> thisclass ;

public:
	static GenericLine from_point_and_direction(const Point& p, const Vector& dir) { return GenericLine(p, dir) ; }
	static GenericLine from_two_points(const Point& p, const Point& q)   { return GenericLine(p, q - p) ; }

	GenericLine() {}
	void set_point(const Point& p) { p_ = p; }
	void set_dir(const Vector& dir) { dir_ = normalize(dir); }

	const Vector& direction() const { return dir_; }

	const Point&  point() const { return p_; }

	GenericLine opposite() const { return GenericLine(p_, -dir_); }

	// the projection of p on this line
	Point  projection(const Point &p) const { return p_ + dir_ * dot(p - p_, dir_) ; }

	FT	   squared_ditance(const Point &p) const { length2(projection(p) - p); }

private:  // Ambiguities exist for this one.
    GenericLine(const Point & p, const Vector & dir) ;

private:
	Point	p_;
	Vector  dir_;
};



template <int DIM, class FT> inline
GenericLine<DIM, FT>::GenericLine(const Point & p, const Vector & dir) : p_(p) {
    dir_ = normalize(dir);

#ifndef NDEBUG // degenerate case
    if (length(dir_) < 1e-15) {
        std::cerr << "degenerate line constructed from point (" 
			<< p << ") and direction (" << dir << ")" << std::endl;
    }
#endif
}


#endif

