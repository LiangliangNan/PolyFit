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

/**@file  slufactor_rational.h
 * @brief Implementation of Sparse Linear Solver with Rational precision.
 */
#ifndef _SLUFACTOR_RATIONAL_H_
#define _SLUFACTOR_RATIONAL_H_

#include <assert.h>

#include "spxdefines.h"
#include "timerfactory.h"
#include "dvector.h"
#include "slinsolver_rational.h"
#include "clufactor_rational.h"
#include "rational.h"

namespace soplex
{
/// maximum nr. of factorization updates allowed before refactorization.
#define MAXUPDATES      1000

/**@brief   Implementation of Sparse Linear Solver with Rational precision.
 * @ingroup Algo
 *
 * This class implements a SLinSolverRational interface by using the sparse LU
 * factorization implemented in CLUFactorRational.
 */
class SLUFactorRational : public SLinSolverRational, protected CLUFactorRational
{
public:

   //--------------------------------
   /**@name Types */
   //@{
   /// Specifies how to perform \ref soplex::SLUFactorRational::change "change" method.
   enum UpdateType
   {
      ETA = 0,       ///<
      FOREST_TOMLIN  ///<
   };
   /// for convenience
   typedef SLinSolverRational::Status Status;
   //@}

private:

   //--------------------------------
   /**@name Private data */
   //@{
   DVectorRational    vec;           ///< Temporary vector
   SSVectorRational   ssvec;         ///< Temporary semi-sparse vector
   //@}

protected:

   //--------------------------------
   /**@name Protected data */
   //@{
   bool       usetup;             ///< TRUE iff update vector has been setup
   UpdateType uptype;             ///< the current \ref soplex::SLUFactor::UpdateType "UpdateType".
   SSVectorRational   eta;        ///<
   SSVectorRational   forest;     ///< ? Update vector set up by solveRight4update() and solve2right4update()
   Rational       lastThreshold;  ///< pivoting threshold of last factorization
   //@}

   //--------------------------------
   /**@name Control Parameters */
   //@{
   /// minimum threshold to use.
   Rational minThreshold;
   /// minimum stability to achieve by setting threshold.
   Rational minStability;
   /// Time spent in solves
   Timer*  solveTime;
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
   void changeEta(int idx, SSVectorRational& eta);
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
   void setMarkowitz(const Rational& m)
   {
      if( m < 0.01 )
      {
         minThreshold = 0.01;
         lastThreshold = 0.01;
      }
      else if( m > 0.99 )
      {
         minThreshold = 0.99;
         lastThreshold = 0.99;
      }
      else
      {
         minThreshold = m;
         lastThreshold = m;
      }
   }

   /// returns Markowitz threshold.
   Rational markowitz()
   {
      return lastThreshold;
   }
   //@}

   //--------------------------------
   /**@name Derived from SLinSolverRational
      See documentation of \ref soplex::SLinSolverRational "SLinSolverRational" for a
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
      return (uptype == SLUFactorRational::ETA) ? "SLU-Eta" : "SLU-Forest-Tomlin";
   }
   ///
   Status status() const
   {
      return Status(stat);
   }
   ///
   Rational stability() const;
   ///
   std::string statistics() const;
   ///
   Status load(const SVectorRational* vec[], int dim);
   //@}

public:

   //--------------------------------
   /**@name Solve */
   //@{
   /// Solves \f$Ax=b\f$.
   void solveRight (VectorRational& x, const VectorRational& b);
   /// Solves \f$Ax=b\f$.
   void solveRight (SSVectorRational& x, const SVectorRational& b);
   /// Solves \f$Ax=b\f$.
   void solveRight4update(SSVectorRational& x, const SVectorRational& b);
   /// Solves \f$Ax=b\f$ and \f$Ay=d\f$.
   void solve2right4update(SSVectorRational& x, VectorRational& y, const SVectorRational& b, SSVectorRational& d);
   /// Solves \f$Ax=b\f$, \f$Ay=d\f$ and \f$Az=e\f$.
   void solve3right4update(SSVectorRational& x, VectorRational& y, VectorRational& z,
                           const SVectorRational& b, SSVectorRational& d, SSVectorRational& e);
   /// Solves \f$Ax=b\f$.
   void solveLeft(VectorRational& x, const VectorRational& b);
   /// Solves \f$Ax=b\f$.
   void solveLeft(SSVectorRational& x, const SVectorRational& b);
   /// Solves \f$Ax=b\f$ and \f$Ay=d\f$.
   void solveLeft(SSVectorRational& x, VectorRational& y, const SVectorRational& b, SSVectorRational& d);
   /// Solves \f$Ax=b\f$, \f$Ay=d\f$ and \f$Az=e\f$.
   void solveLeft(SSVectorRational& x, VectorRational& y, VectorRational& z,
                  const SVectorRational& b, SSVectorRational& d, SSVectorRational& e);
   ///
   Status change(int idx, const SVectorRational& subst, const SSVectorRational* eta = 0);
   //@}

   //--------------------------------
   /**@name Miscellaneous */
   //@{
   /// time spent in factorizations
   Real getFactorTime() const
   {
      return factorTime->time();
   }
   /// set time limit on factorization
   void setTimeLimit(const Real limit)
   {
      timeLimit = limit;
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
   SLUFactorRational();
   /// assignment operator.
   SLUFactorRational& operator=(const SLUFactorRational& old);
   /// copy constructor.
   SLUFactorRational(const SLUFactorRational& old);
   /// destructor.
   virtual ~SLUFactorRational();
   /// clone function for polymorphism
   inline virtual SLinSolverRational* clone() const
   {
      return new SLUFactorRational(*this);
   }
   //@}

private:

   //------------------------------------
   /**@name Private helpers */
   //@{
   /// used to implement the assignment operator
   void assign(const SLUFactorRational& old);
   //@}
};

} // namespace soplex
#endif // _SLUFACTOR_RATIONAL_H_
