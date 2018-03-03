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


#ifndef _MATH_GEOMETRY_BOX_H_
#define _MATH_GEOMETRY_BOX_H_

#include "math_common.h"
#include "vecg.h"
#include "../basic/assertions.h"


template <class T>
class GenericBox2d {
public:
	GenericBox2d() : initialized_(false), x_min_(1e30), y_min_(1e30), x_max_(-1e30), y_max_(-1e30) { }

	bool initialized() const { return initialized_; }
	void clear() { initialized_ = false; }

	T x_min() const { ogf_debug_assert(initialized_); return x_min_; }
	T y_min() const { ogf_debug_assert(initialized_); return y_min_; }
	T x_max() const { ogf_debug_assert(initialized_); return x_max_; }
	T y_max() const { ogf_debug_assert(initialized_); return y_max_; }

	T min(unsigned axis) const { return (axis == 0) ? x_min_ : y_min_; }
	T max(unsigned axis) const { return (axis == 0) ? x_max_ : y_max_; }

	T width()  const { return x_max() - x_min(); }
	T height() const { return y_max() - y_min(); }

	T area() const { return width() * height(); }

	vecng<2, T> center() const {
		return vecng<2, T>((x_max() + x_min()) / 2, (y_max() + y_min()) / 2);
	}
	T radius() const {
        return T(0.5) * ::sqrt(ogf_sqr(x_max() - x_min()) + ogf_sqr(y_max() - y_min()));
	}


	void add_point(const vecng<2, T>& p) {
		if (!initialized_) {
			x_min_ = p.x;
			y_min_ = p.y;
			x_max_ = p.x;
			y_max_ = p.y;
			initialized_ = true;
		}
		else {
			x_min_ = ogf_min(x_min_, p.x);
			y_min_ = ogf_min(y_min_, p.y);
			x_max_ = ogf_max(x_max_, p.x);
			y_max_ = ogf_max(y_max_, p.y);
		}
	}

	void add_box(const GenericBox2d<T>& b) {
		if (b.initialized()) {
			add_point(vecng<2, T>(b.x_min(), b.y_min()));
			add_point(vecng<2, T>(b.x_max(), b.y_max()));
		}
	}

private:
	bool initialized_;
	T x_min_;
	T y_min_;
	T x_max_;
	T y_max_;
};

//_________________________________________________________________________

template <class T>
class GenericBox3d {
public:
	GenericBox3d() : initialized_(false),
		x_min_(1e30f), y_min_(1e30f), z_min_(1e30f),
		x_max_(-1e30f), y_max_(-1e30f), z_max_(-1e30f) {
	}

	bool initialized() const { return initialized_; }
	void clear() { initialized_ = false; }

	T x_min() const { ogf_debug_assert(initialized_); return x_min_; }
	T y_min() const { ogf_debug_assert(initialized_); return y_min_; }
	T z_min() const { ogf_debug_assert(initialized_); return z_min_; }
	T x_max() const { ogf_debug_assert(initialized_); return x_max_; }
	T y_max() const { ogf_debug_assert(initialized_); return y_max_; }
	T z_max() const { ogf_debug_assert(initialized_); return z_max_; }

	T min(unsigned axis) const { return (axis == 0) ? x_min_ : ((axis == 1) ? y_min_ : z_min_); }
	T max(unsigned axis) const { return (axis == 0) ? x_max_ : ((axis == 1) ? y_max_ : z_max_); }

	T width()  const { return x_max() - x_min(); }
	T height() const { return y_max() - y_min(); }
	T depth()  const { return z_max() - z_min(); }

	T area() const { return T(2.0) * (width() * height() + height() * depth() + depth() * width()); }

	vecng<3, T> center() const {
		return vecng<3, T>(
			T(0.5)*(x_max() + x_min()),
			T(0.5)*(y_max() + y_min()),
			T(0.5)*(z_max() + z_min())
			);
	}

	T radius() const {
		return T(0.5) * ::sqrt(
			ogf_sqr(x_max() - x_min()) + ogf_sqr(y_max() - y_min()) + ogf_sqr(z_max() - z_min())
			);
	}

	void add_point(const vecng<3, T>& p) {
		if (!initialized_) {
			x_min_ = p.x;
			y_min_ = p.y;
			z_min_ = p.z;
			x_max_ = p.x;
			y_max_ = p.y;
			z_max_ = p.z;
			initialized_ = true;
		}
		else {
			x_min_ = ogf_min(x_min_, p.x);
			y_min_ = ogf_min(y_min_, p.y);
			z_min_ = ogf_min(z_min_, p.z);
			x_max_ = ogf_max(x_max_, p.x);
			y_max_ = ogf_max(y_max_, p.y);
			z_max_ = ogf_max(z_max_, p.z);
		}
	}

	void add_box(const GenericBox3d<T>& b) {
		add_point(vecng<3, T>(b.x_min(), b.y_min(), b.z_min()));
		add_point(vecng<3, T>(b.x_max(), b.y_max(), b.z_max()));
	}

private:
	bool initialized_;
	T x_min_;
	T y_min_;
	T z_min_;
	T x_max_;
	T y_max_;
	T z_max_;
};



#endif
