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


#include "principal_axes.h"
#include "semi_definite_symmetric_eigen.h"


PrincipalAxes3d::PrincipalAxes3d() {

}

void PrincipalAxes3d::begin() {
	nb_points_ = 0 ;
	sum_weights_ = 0 ;
	center_[0] = center_[1] = center_[2] = 0 ;
	M_[0] = M_[1] = M_[2] = M_[3] = M_[4] = M_[5] = 0 ;
}

void PrincipalAxes3d::end() {
	center_[0] /= sum_weights_ ;
	center_[1] /= sum_weights_ ;
	center_[2] /= sum_weights_ ;

	// If the system is under-determined, 
	//   return the trivial basis.
	if(nb_points_ < 4) {
		axis_[0] = vec3(1,0,0) ;
		axis_[1] = vec3(0,1,0) ;
		axis_[2] = vec3(0,0,1) ;
		eigen_value_[0] = 1.0 ;
		eigen_value_[1] = 1.0 ;
		eigen_value_[2] = 1.0 ;
	} else {
		double x = center_[0] ;
		double y = center_[1] ;
		double z = center_[2] ;

		M_[0] = M_[0]/sum_weights_ - x*x ;
		M_[1] = M_[1]/sum_weights_ - x*y ;
		M_[2] = M_[2]/sum_weights_ - y*y ;
		M_[3] = M_[3]/sum_weights_ - x*z ;
		M_[4] = M_[4]/sum_weights_ - y*z ;
		M_[5] = M_[5]/sum_weights_ - z*z ;

		if( M_[0] <= 0 ) {
			M_[0] = 1.e-30 ; 
		}
		if( M_[2] <= 0 ) {
			M_[2] = 1.e-30 ; 
		}
		if( M_[5] <= 0 ) {
			M_[5] = 1.e-30 ; 
		}

		double eigen_vectors[9] ;
		MatrixUtil::eigen_symmetric(M_, 3, eigen_vectors, eigen_value_) ;

		axis_[0] = vec3(
			float(eigen_vectors[0]), float(eigen_vectors[1]), float(eigen_vectors[2])
		) ;

		axis_[1] = vec3(
			float(eigen_vectors[3]), float(eigen_vectors[4]), float(eigen_vectors[5])
		) ;

		axis_[2] = vec3(
			float(eigen_vectors[6]), float(eigen_vectors[7]), float(eigen_vectors[8])
		) ;

		// Normalize the eigen vectors
		for(int i=0; i<3; i++) {
			axis_[i] = normalize(axis_[i]) ;
		}
	}

}

void PrincipalAxes3d::add_point(const vec3& p, double weight) {
	center_[0] += p.x * weight ;
	center_[1] += p.y * weight ;
	center_[2] += p.z * weight ;

	double x = p.x ;
	double y = p.y ; 
	double z = p.z ;

	M_[0] += weight * x*x ;
	M_[1] += weight * x*y ;
	M_[2] += weight * y*y ;
	M_[3] += weight * x*z ;
	M_[4] += weight * y*z ;
	M_[5] += weight * z*z ;

	nb_points_++ ;
	sum_weights_ += weight ;
}

//_________________________________________________________

PrincipalAxes2d::PrincipalAxes2d() {

}

void PrincipalAxes2d::begin() {
	nb_points_ = 0 ;
	sum_weights_ = 0 ;
	center_[0] = center_[1] = 0 ;
	M_[0] = M_[1] = M_[2] = 0 ;
}

void PrincipalAxes2d::end() {

	center_[0] /= sum_weights_ ;
	center_[1] /= sum_weights_ ;

	// If the system is under-determined, 
	//  return the trivial basis.
	if(nb_points_ < 3) {
		axis_[0] = vec2(1,0) ;
		axis_[1] = vec2(0,1) ;
		eigen_value_[0] = 1.0 ;
		eigen_value_[1] = 1.0 ;
	} else {
		double x = center_[0] ;
		double y = center_[1] ;

		M_[0] = M_[0]/sum_weights_ - x*x ;
		M_[1] = M_[1]/sum_weights_ - x*y ;
		M_[2] = M_[2]/sum_weights_ - y*y ;

		if( M_[0] <= 0 ) {
			M_[0] = 1.e-30 ; 
		}

		if( M_[2] <= 0 ) {
			M_[2] = 1.e-30 ; 
		}

		double eigen_vectors[4] ;
		MatrixUtil::eigen_symmetric(M_, 2, eigen_vectors, eigen_value_) ;

		axis_[0] = vec2(float(eigen_vectors[0]), float(eigen_vectors[1]));
		axis_[1] = vec2(float(eigen_vectors[2]), float(eigen_vectors[3]));

		// Normalize the eigen vectors
		for(int i=0; i<2; i++) {
			axis_[i] = normalize(axis_[i]) ;
		}
	}

}

void PrincipalAxes2d::add_point(const vec2& p, double weight) {

	double x = p.x ;
	double y = p.y ; 

	center_[0] += x * weight ;
	center_[1] += y * weight ;

	M_[0] += weight * x*x ;
	M_[1] += weight * x*y ;
	M_[2] += weight * y*y ;

	nb_points_++ ;
	sum_weights_ += weight ;
}
