/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  slinsolver.h
 * @brief Sparse Linear Solver virtual base class.
 */
#ifndef _SLINSOLVER_H_
#define _SLINSOLVER_H_


#include <assert.h>
#include <string.h>

#include "spxdefines.h"
#include "svector.h"
#include "ssvector.h"
#include "dsvector.h"
#include "dvector.h"
#include "didxset.h"

namespace soplex
{
/**@brief   Sparse Linear Solver virtual base class.
   @ingroup Algo
   
   Class SLinSolver provides a class for solving sparse linear systems with
   a matrix \f$A\f$ and arbitrary right-hand side vectors. For doing so, the
   matrix must be first #load%ed to an #SLinSolver object as an array of
   pointers to the \em column \ref SVector "SVectors" of this matrix.
*/
class SLinSolver
{
public:

   //---------------------------------------
   /**@name Types */
   //@{
   /// status flags of the SLinSolver class.
   enum Status
   {
      /** The SLinSolver is ready for solving linear systems with the
          loaded matrix */
      OK       = 0,
      /** The loaded matrix allows only for instable solutions to be
          computed */
      INSTABLE = 1,
      /// The loaded matrix is singular.
      SINGULAR = 2,
      /// No matrix has yet been loaded.
      UNLOADED = 4,
      /// An error has occurred.
      ERROR    = 8
   };
   //@}

   //---------------------------------------
   /**@name Miscellaneous */
   //@{
   /// returns the name of the SLinSolver.
   virtual const char* getName() const = 0;

   /// returns the Status of the SLinSolver.
   virtual Status status() const = 0;

   /// unloads any matrix.
   virtual void clear() = 0;

   /// returns current memory consumption.
   virtual int memory() const = 0;

   /// returns dimension of loaded matrix.
   virtual int dim() const = 0;

   /// loads \p dim column vectors \p vec into the solver.
   /** Initializes SLinSolver for the solution of linear systems
       with the matrix consisting of \p dim column vectors given in \p vec.
   */
   virtual Status load(const SVector* vec[], int dim) = 0;
   
   /// returns a stability number (0: singularity, 1: perfect stability).
   /** Returns a stability parameter between 0 and 1, where 0 indicates
       singularity, while 1 indicates perfect stability.
   */
   virtual Real stability() const = 0;

   /// return estimate for the condition number based on the diagonal of U
   virtual Real conditionEstimate(int type = 0) const = 0;

   /// returns statistical information in form of a string.
   virtual std::string statistics() const = 0;

   /// Substitute column \p idx with \p subst.
   /** The change method is used to modify the loaded matrix by substituting
       column \p idx with the new vector \p subst. One may also pass the
       optional parameter \p eta to the solution of #solveRight() if
       readily  availabble. This may improve on the performance of the update.
   */
   virtual Status change(int idx, const SVector& subst, const SSVector* eta = 0) = 0;

   /// consistency check.
   virtual bool isConsistent() const = 0;

   /// get number of factorizations
   virtual int getFactorCount() const = 0;
   //@}


   /**@name Solving linear systems
      For solving linear systems with an SLinSolver object, it must
      have previously been loaded with the matrix to use.
      
      Two types of systems can be solved \f$A x = b\f$ and \f$x^T A = b^T\f$.
      Method names related to the first and second type are solveRight() and
      solveLeft(), respectively.
      
      The methods receive their right hand-side vector \f$b\f$ as a
      \c const parameter, that will hence be unchanged after termination.
      
      Some methods are available with two parameters for right hand-side
      vectors. Then two system are solved in one method invocation. This
      should generally be faster than solving two systems seperately.
      
      The result vector(s) are allways given as the first parameter(s). Two
      types of result vectors are supported, Vector and SSVector.
   */
   //@{
   /// Solves \f$Ax=b\f$.
   virtual void solveRight (Vector& x, const Vector& b) /* const */ = 0;
   /// Solves \f$Ax=b\f$.
   virtual void solveRight (SSVector& x, const SSVector& b) /* const */ = 0;
   virtual void solveRight (SSVector& x, const SVector& b) /* const */ = 0;

   /** @brief Solves \f$Ax=b\f$.
       Possibly sets up internal data structures suitable for an optimized
       subsequent change() call with \f$b\f$ as entering column.
   */
   virtual void solveRight4update(SSVector& x, const SVector& b) = 0;
   
   /// Solves \f$Ax=b\f$ and \f$Ay=d\f$.
   virtual void solve2right4update(SSVector& x,
                                   Vector& y,
                                   const SVector& b,
                                   SSVector& d ) = 0;
   /// sparse version of solving two systems of equations
   virtual void solve2right4update(SSVector& x,
                                   SSVector& y,
                                   const SVector& b,
                                   SSVector& d ) = 0;
   /// Solves \f$Ax=b\f$, \f$Ay=d\f$ and \f$Az=e\f$.
   virtual void solve3right4update(SSVector& x,
                                   Vector& y,
                                   Vector& z,
                                   const SVector& b,
                                   SSVector& d,
                                   SSVector& e) = 0;
   /// sparse version of solving three systems of equations
   virtual void solve3right4update(SSVector& x,
                                   SSVector& y,
                                   SSVector& z,
                                   const SVector& b,
                                   SSVector& d,
                                   SSVector& e) = 0;
   /// solves \f$x^TA=b^T\f$.
   virtual void solveLeft (Vector& x, const Vector& b) /* const */ = 0;
   virtual void solveLeft (SSVector& x, const SSVector& b) /* const */ = 0;
   /// sparse version of solving one system of equations with transposed basis matrix
   virtual void solveLeft (SSVector& x, const SVector& b) /* const */ = 0;
   /// solves \f$x^TA=b^T\f$ and \f$x^TA=rhs2^T\f$ internally using \f$rhs2\f$.
   virtual void solveLeft (SSVector& x,
                           Vector& two,
                           const SVector& b,
                           SSVector& rhs2) /* const */ = 0;
   /// sparse version of solving two systems of equations with transposed basis matrix
   virtual void solveLeft (SSVector& x,
                           SSVector& two,
                           const SVector& b,
                           SSVector& rhs2) /* const */ = 0;
   /// solves \f$x^TA=b^T\f$, \f$y^TA=d^T\f$ and \f$z^TA=e^T\f$
   virtual void solveLeft (SSVector& x, Vector& y, Vector& z,
                           const SVector& b, SSVector& d, SSVector& e) = 0;
   /// sparse version of solving three systems of equations with transposed basis matrix
   virtual void solveLeft (SSVector& x, SSVector& y, SSVector& z,
                           const SVector& b, SSVector& d, SSVector& e) = 0;
   //@}


   //---------------------------------------
   /**@name Constructors / Destructors */
   //@{
   /// default constructor
   SLinSolver()
      : spxout(0)
   {}
   /// destructor
   virtual ~SLinSolver()
   {}
   /// clone function for polymorphism
   virtual SLinSolver* clone() const = 0;
   //@}

   /// message handler
   SPxOut* spxout;


};

} // namespace soplex
#endif // _SLINSOLVER_H_
