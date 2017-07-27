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


#ifndef __MESH_TOOLS_MATH_PPAL_AXIS__
#define __MESH_TOOLS_MATH_PPAL_AXIS__

#include "math_common.h"
#include "math_types.h"
#include "../basic/assertions.h"

// eigen values are sorted in descending order, 
// eigen vectors are sorted in accordance.

/**
* PrincipalAxes3d enables the center and inertia axes of
* a cloud of 3d points to be computed.
*/
class MATH_API PrincipalAxes3d {
public:
	PrincipalAxes3d() ;
	void begin() ;
	void add_point(const vec3& p, double weight = 1.0) ;
	void end() ;

	vec3 center() const ;
	const vec3& axis(int i) const ;
	double eigen_value(int i) const ; 


private:
	double	center_[3] ;
	vec3	axis_[3] ;
	double	eigen_value_[3] ;

	double	M_[6] ;
	int		nb_points_ ;
	double	sum_weights_ ;
} ;

//_________________________________________________________

/**
* PrincipalAxes2d enables the center and inertia axes of
* a cloud of 2d points to be computed.
*/
class MATH_API PrincipalAxes2d {
public:
	PrincipalAxes2d() ;
	void begin() ;
	void add_point(const vec2& p, double weight = 1.0) ;
	void end() ;

	vec2 center() const ;
	const vec2& axis(int i) const ;
	double eigen_value(int i) const ; 


private:
	double	center_[2] ;
	vec2	axis_[2] ;
	double	eigen_value_[2] ;

	double	M_[3] ;
	int		nb_points_ ;
	double	sum_weights_ ;
} ;

//_________________________________________________________

inline vec3 PrincipalAxes3d::center() const {
	return vec3(float(center_[0]), float(center_[1]), float(center_[2]));
}

inline const vec3& PrincipalAxes3d::axis(int i) const {
	ogf_assert(i >= 0 && i < 3) ;
	return axis_[i] ;
}

inline double PrincipalAxes3d::eigen_value(int i) const {
	ogf_assert(i >= 0 && i < 3) ;
	return eigen_value_[i] ;
}

//______________________________________________________

inline vec2 PrincipalAxes2d::center() const {
	return vec2(float(center_[0]), float(center_[1]));
}

inline const vec2& PrincipalAxes2d::axis(int i) const {
	ogf_assert(i >= 0 && i < 2) ;
	return axis_[i] ;
}

inline double PrincipalAxes2d::eigen_value(int i) const {
	ogf_assert(i >= 0 && i < 2) ;
	return eigen_value_[i] ;
}


#endif

