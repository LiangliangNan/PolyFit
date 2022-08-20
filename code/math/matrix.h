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



#ifndef __MATH_LINEAR_ALGEBRA_MATRIX__
#define __MATH_LINEAR_ALGEBRA_MATRIX__

#include "../basic/assertions.h"
#include <iostream>



/**
* A class for representing matrices where the coefficients
* may be of arbitrary types.
*/

template <int DIM, class FT> class Matrix ;
template <int DIM, class FT> Matrix<DIM, FT> operator*(
	FT op1, const Matrix<DIM, FT>& op2
	) ;


template <int DIM, class FT> 
class Matrix 
{
public:
	Matrix() ;
	void load_zero() ;
	void load_identity() ;

	FT& operator()(int i, int j) ;
	const FT& operator()(int i, int j) const ;

	Matrix<DIM, FT>& operator+=(const Matrix<DIM, FT>& rhs) ;
	Matrix<DIM, FT>& operator-=(const Matrix<DIM, FT>& rhs) ;
	Matrix<DIM, FT>& operator*=(FT rhs) ;
	Matrix<DIM, FT>& operator/=(FT rhs) ;

	Matrix<DIM, FT> operator+(const Matrix<DIM, FT>& op2) const ;
	Matrix<DIM, FT> operator-(const Matrix<DIM, FT>& op2) const ;
	Matrix<DIM, FT> operator*(const Matrix<DIM, FT>& op2) const ;

	Matrix<DIM, FT> inverse() const ;
	Matrix<DIM, FT> transpose() const ;

	// cast to Scalar array
	operator FT*() { return &(coeff_[0][0]); }

	// cast to const Scalar array
	operator const FT*() const { return &(coeff_[0][0]); }

	// Routines for interfacing with Fortran, OpenGL etc...
	const FT* data() const { return &(coeff_[0][0]) ; }
	FT* data() { return &(coeff_[0][0]) ; }

	void get_lower_triangle(FT* store) {
		for(unsigned int i=0; i<DIM; i++) {
			for(unsigned int j=0; j<=i; j++) {
				*store++ = coeff_[i][j] ;
			}
		}
	}

private:
	FT coeff_[DIM][DIM] ;
} ;

//_______________________________________________________________________

template <int DIM, class FT> inline 
Matrix<DIM, FT>::Matrix() {
	load_identity() ;
}


template <int DIM, class FT> inline 
FT& Matrix<DIM, FT>::operator()(int i, int j) {
	ogf_assert(i >= 0 && i < DIM) ;
	ogf_assert(j >= 0 && j < DIM) ;
	return coeff_[i][j] ;
}

template <int DIM, class FT> inline 
const FT& Matrix<DIM, FT>::operator()(int i, int j) const {
	ogf_assert(i >= 0 && i < DIM) ;
	ogf_assert(j >= 0 && j < DIM) ;
	return coeff_[i][j] ;
}

template <int DIM, class FT> inline 
Matrix<DIM, FT>& Matrix<DIM, FT>::operator+=(const Matrix<DIM, FT>& rhs) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] += rhs.coeff_[i][j] ;
		}
	}
	return *this ;
}

template <int DIM, class FT> inline 
Matrix<DIM, FT>& Matrix<DIM, FT>::operator-=(const Matrix<DIM, FT>& rhs) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] -= rhs.coeff_[i][j] ;
		}
	}
	return *this ;
}

template <int DIM, class FT> inline 
Matrix<DIM, FT>& Matrix<DIM, FT>::operator*=(FT rhs) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] *= rhs ;
		}
	}
	return *this ;
}

template <int DIM, class FT> inline 
Matrix<DIM, FT>& Matrix<DIM, FT>::operator/=(FT rhs) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] /= rhs ;
		}
	}
	return *this ;
}

template <int DIM, class FT> inline 
Matrix<DIM, FT> Matrix<DIM, FT>::operator+(const Matrix<DIM, FT>& op2) const {
	Matrix<DIM, FT> result = *this ;
	result += op2 ;
	return result ;
}

template <int DIM, class FT> inline 
Matrix<DIM, FT> Matrix<DIM, FT>::operator-(const Matrix<DIM, FT>& op2) const {
	Matrix<DIM, FT> result = *this ;
	result -= op2 ;
	return result ;

}

template <int DIM, class FT> inline 
Matrix<DIM, FT> Matrix<DIM, FT>::operator*(const Matrix<DIM, FT>& op2) const {
	Matrix<DIM, FT> result ;
	result.load_zero() ;
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			for(int k=0; k<DIM; k++) {
				result.coeff_[i][j] += coeff_[i][k] * op2.coeff_[k][j] ;
			}
		}
	}
	return result ;
}

template <int DIM, class FT> inline 
std::ostream& operator << (std::ostream& output, const Matrix<DIM, FT>& m) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			output << m(i,j) << " " ;
		}
	}    
	return output ;
}

template <int DIM, class FT> inline 
std::istream& operator >> (std::istream& input, Matrix<DIM, FT>& m) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			input >> m(i,j) ;
		}
	}    
	return input ;
}



//_______________________________________________________________________

template <int DIM, class FT> void
Matrix<DIM, FT>::load_zero() {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] = FT(0) ;
		}
	}
}

template <int DIM, class FT> void
Matrix<DIM, FT>::load_identity() {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] = (i==j) ? FT(1) : FT(0) ;
		}
	}
}

template <int DIM, class FT> Matrix<DIM, FT> 
Matrix<DIM, FT>::inverse() const {
	FT val, val2;
	int i, j, k, ind;
	Matrix<DIM, FT> tmp = (*this);
	Matrix<DIM, FT> result;

	result.load_identity();

	for (i = 0; i != DIM; i++) {
		val = tmp(i,i);			/* find pivot */
		ind = i;
		for (j = i + 1; j != DIM; j++) {
			if (std::fabs(tmp(j,i)) > std::fabs(val)) {
				ind = j;
				val = tmp(j,i);
			}
		}

		if (ind != i) {			
			for (j = 0; j != DIM; j++) {
				val2 = result(i,j);
				result(i,j) = result(ind,j);
				result(ind,j) = val2;           /* swap columns */
				val2 = tmp(i,j);
				tmp(i,j) = tmp(ind,j);
				tmp(ind,j) = val2;
			}
		}

		ogf_assert(val != 0.0);

		for (j = 0; j != DIM; j++) {
			tmp(i,j)    /= val;
			result(i,j) /= val;
		}

		for (j = 0; j != DIM; j++) {		
			if (j == i)
				continue;                       /* eliminate column */
			val = tmp(j,i);
			for (k = 0; k != DIM; k++) {
				tmp(j,k)     -= tmp(i,k)     * val;
				result(j,k)  -= result(i,k)  * val;
			}
		}
	}

	return result;
}

template <int DIM, class FT> Matrix<DIM, FT> 
Matrix<DIM, FT>::transpose() const {
	Matrix<DIM, FT> result ;
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			result(i,j) = (*this)(j,i) ;
		}
	}
	return result ;
}

template <int N, class FT> inline
void mult(const Matrix<N, FT>& M, const FT* x, FT* y) {
	for(int i=0; i<N; i++) {
		y[i] = 0 ;
		for(int j=0; j<N; j++) {
			y[i] += M(i,j)*x[j] ;
		}
	}
}


#endif
