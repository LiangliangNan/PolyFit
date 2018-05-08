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

/**@file  slufactor.h
 * @brief Implementation of Sparse Linear Solver.
 */
#ifndef _SLUFACTOR_H_
#define _SLUFACTOR_H_

#include <assert.h>

#include "spxdefines.h"
#include "timerfactory.h"
#include "dvector.h"
#include "slinsolver.h"
#include "clufactor.h"

namespace soplex
{
/// maximum nr. of factorization updates allowed before refactorization.
#define MAXUPDATES      1000

/**@brief   Implementation of Sparse Linear Solver.
 * @ingroup Algo
 * 
 * This class implements a SLinSolver interface by using the sparse LU
 * factorization implemented in CLUFactor.
 */
class SLUFactor : public SLinSolver, protected CLUFactor
{
public:

   //--------------------------------
   /**@name Types */
   //@{
   /// Specifies how to perform \ref soplex::SLUFactor::change "change" method.
   enum UpdateType
   {
      ETA = 0,       ///<
      FOREST_TOMLIN  ///<
   };
   /// for convenience
   typedef SLinSolver::Status Status;
   //@}

private:

   //--------------------------------
   /**@name Private data */
   //@{
   DVector    vec;           ///< Temporary vector
   SSVector   ssvec;         ///< Temporary semi-sparse vector
   //@}

protected:

   //--------------------------------
   /**@name Protected data */
   //@{
   bool       usetup;        ///< TRUE iff update vector has been setup
   UpdateType uptype;        ///< the current \ref soplex::SLUFactor::UpdateType "UpdateType".
   SSVector   eta;           ///< 
   SSVector   forest;        ///< ? Update vector set up by solveRight4update() and solve2right4update()
   Real       lastThreshold; ///< pivoting threshold of last factorization
   //@}

   //--------------------------------
   /**@name Control Parameters */
   //@{
   /// minimum threshold to use.
   Real minThreshold;
   /// minimum stability to achieve by setting threshold.
   Real minStability;
   /// |x| < epsililon is considered to be 0.
   Real epsilon;
   /// Time spent in solves
   Timer* solveTime;
   Timer::TYPE timerType;
   /// Number of solves
   int     solveCount;
   //@}

protected:

   //--------------------------------
   /**@name Protected helpers */
   //@{
   ///
   void freeAll();
   ///
   void changeEta(int idx, SSVector& eta);
   //@}


public:

   //--------------------------------
   /**@name Update type */
   //@{
   /// returns the current update type uptype.
   UpdateType utype() const
   {
      return uptype;
   }

   /// sets update type.
   /** The new UpdateType becomes valid only after the next call to
       method load().
   */
   void setUtype(UpdateType tp)
   {
      uptype = tp;
   }

   /// sets minimum Markowitz threshold.
   void setMarkowitz(Real m)
   {
      if( m < 0.01 )
         m = 0.01;

      if( m > 0.99 )
         m = 0.99;

      minThreshold = m;
      lastThreshold = m;
   }

   /// returns Markowitz threshold.
   Real markowitz()
   {
      return lastThreshold;
   }
   //@}

   //--------------------------------
   /**@name Derived from SLinSolver
      See documentation of \ref soplex::SLinSolver "SLinSolver" for a 
      documentation of these methods.
   */
   //@{
   ///
   void clear();
   ///
   int dim() const
   {
      return thedim;
   }
   ///
   int memory() const
   {
      return nzCnt + l.start[l.firstUnused];
   }
   ///
   const char* getName() const
   {
      return (uptype == SLUFactor::ETA) ? "SLU-Eta" : "SLU-Forest-Tomlin";
   }
   ///
   Status status() const
   {
      return Status(stat);
   }
   ///
   Real stability() const;
   /// return condition number estimate based on the diagonal of U
   Real conditionEstimate(int type = 0) const;
   ///
   std::string statistics() const;
   ///
   Status load(const SVector* vec[], int dim);
   //@}

public:

   //--------------------------------
   /**@name Solve */
   //@{
   /// Solves \f$Ax=b\f$.
   void solveRight (Vector& x, const Vector& b);
   void solveRight(SSVector& x, const SSVector& b)
   {
      x.unSetup();
      solveRight((Vector&) x, (const Vector&) b);
   }
   /// Solves \f$Ax=b\f$.
   void solveRight (SSVector& x, const SVector& b);
   /// Solves \f$Ax=b\f$.
   void solveRight4update(SSVector& x, const SVector& b);
   /// Solves \f$Ax=b\f$ and \f$Ay=d\f$.
   void solve2right4update(SSVector& x, Vector& y, const SVector& b, SSVector& d);
   /// Sparse version of solving two systems of equations
   void solve2right4update(SSVector& x, SSVector& y, const SVector& b, SSVector& d);
   /// Solves \f$Ax=b\f$, \f$Ay=d\f$ and \f$Az=e\f$.
   void solve3right4update(SSVector& x, Vector& y, Vector& z,
                           const SVector& b, SSVector& d, SSVector& e);
   /// sparse version of solving three systems of equations
   void solve3right4update(SSVector& x, SSVector& y, SSVector& z,
                           const SVector& b, SSVector& d, SSVector& e);
   /// sparse version of solving one system of equations with transposed basis matrix
   void solveLeft(Vector& x, const Vector& b);
   void solveLeft(SSVector& x, const SSVector& b)
   {
      x.unSetup();
      solveLeft((Vector&) x, (const Vector&) b);
   }
   /// Solves \f$Ax=b\f$.
   void solveLeft(SSVector& x, const SVector& b);
   /// Solves \f$Ax=b\f$ and \f$Ay=d\f$.
   void solveLeft(SSVector& x, Vector& y, const SVector& b, SSVector& d);
   /// sparse version of solving two systems of equations with transposed basis matrix
   void solveLeft(SSVector& x, SSVector& two, const SVector& b, SSVector& rhs2);
   /// Solves \f$Ax=b\f$, \f$Ay=d\f$ and \f$Az=e\f$.
   void solveLeft(SSVector& x, Vector& y, Vector& z,
                  const SVector& b, SSVector& d, SSVector& e);
   /// sparse version of solving three systems of equations with transposed basis matrix
   void solveLeft(SSVector& x, SSVector& y, SSVector& z,
                  const SVector& b, SSVector& d, SSVector& e);
   ///
   Status change(int idx, const SVector& subst, const SSVector* eta = 0);
   //@}

   //--------------------------------
   /**@name Miscellaneous */
   //@{
   /// time spent in factorizations
   Real getFactorTime() const
   {
      return factorTime->time();
   }
   /// reset FactorTime
   void resetFactorTime()
   {
      factorTime->reset();
   }
   /// number of factorizations performed
   int getFactorCount() const
   {
      return factorCount;
   }
   /// time spent in solves
   Real getSolveTime() const
   {
      return solveTime->time();
   }
   /// reset SolveTime
   void resetSolveTime()
   {
      solveTime->reset();
   }
   /// number of solves performed
   int getSolveCount() const
   {
      return solveCount;
   }
   /// reset timers and counters
   void resetCounters()
   {
      factorTime->reset();
      solveTime->reset();
      factorCount = 0;
      solveCount = 0;
   }
   /// prints the LU factorization to stdout.
   void dump() const;

   /// consistency check.
   bool isConsistent() const;
   //@}

   //------------------------------------
   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   SLUFactor();
   /// assignment operator.
   SLUFactor& operator=(const SLUFactor& old);
   /// copy constructor.
   SLUFactor(const SLUFactor& old);
   /// destructor.
   virtual ~SLUFactor();
   /// clone function for polymorphism
   inline virtual SLinSolver* clone() const
   {
      return new SLUFactor(*this);
   }
   //@}

private:

   //------------------------------------
   /**@name Private helpers */
   //@{
   /// used to implement the assignment operator
   void assign(const SLUFactor& old);
   //@}
};

} // namespace soplex
#endif // _SLUFACTOR_H_
