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

/**@file  spxdevexpr.h
 * @brief Devex pricer.
 */
#ifndef _SPXDEVEXPR_H_
#define _SPXDEVEXPR_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxpricer.h"

namespace soplex
{

/**@brief   Devex pricer.
   @ingroup Algo

   The Devex Pricer for SoPlex implements an approximate steepest edge pricing,
   that does without solving an extra linear system and computing the scalar
   products.

   See SPxPricer for a class documentation.

   @todo There seem to be problems with this pricer especially on the 
         greenbe[ab] problems with the entering algorithm 
         (row representation?).
*/
class SPxDevexPR : public SPxPricer
{
private:

   //-------------------------------------
   /**@name Data */
   //@{
   Real  last;           ///< penalty, selected at last iteration.
   DataArray<IdxElement> prices;   ///< temporary array of precomputed pricing values
   DataArray<IdxElement> pricesCo; ///< temporary array of precomputed pricing values
   DIdxSet bestPrices;   ///< set of best pricing candidates
   DIdxSet bestPricesCo; ///< set of best pricing candidates
   bool refined;         ///< has a refinement step already been tried?
   ///@}

   //-------------------------------------
   /**@name Private helpers */
   //@{
   /// set entering/leaving algorithm
   void setupWeights(SPxSolver::Type);
   /// build up vector of pricing values for later use
   int buildBestPriceVectorLeave(Real feastol);
   /// internal implementation of SPxPricer::selectLeave()
   int selectLeaveX(Real feastol, int start = 0, int incr = 1);
   /// implementation of sparse pricing in the leaving Simplex
   int selectLeaveSparse(Real feastol);
   /// implementation of hyper sparse pricing in the leaving Simplex
   int selectLeaveHyper(Real feastol);
   /// build up vector of pricing values for later use
   SPxId buildBestPriceVectorEnterDim(Real& best, Real feastol);
   SPxId buildBestPriceVectorEnterCoDim(Real& best, Real feastol);
   /// choose the best entering index among columns and rows but prefer sparsity
   SPxId selectEnterX(Real tol);
   /// implementation of sparse pricing in the entering Simplex (slack variables)
   SPxId selectEnterSparseDim(Real& best, Real feastol);
   /// implementation of sparse pricing in the entering Simplex
   SPxId selectEnterSparseCoDim(Real& best, Real feastol);
   /// SPxPricer::selectEnter() in dense case (slack variabels)
   SPxId selectEnterDenseDim(Real& best, Real feastol, int start = 0, int incr = 1);
   /// SPxPricer::selectEnter() in dense case
   SPxId selectEnterDenseCoDim(Real& best, Real feastol, int start = 0, int incr = 1);
   /// implementation of hyper sparse pricing in the entering Simplex
   SPxId selectEnterHyperDim(Real& best, Real feastol);
   /// implementation of hyper sparse pricing in the entering Simplex
   SPxId selectEnterHyperCoDim(Real& best, Real feastol);
   //@}

public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   SPxDevexPR() 
      : SPxPricer("Devex")
      , last(0)
      , refined(false)
   {}
   /// copy constructor
   SPxDevexPR( const SPxDevexPR& old)
      : SPxPricer(old)
      , last(old.last)
      , refined(false)
   {}
   /// assignment operator
   SPxDevexPR& operator=( const SPxDevexPR& rhs)
   {
      if(this != &rhs)
      {
         SPxPricer::operator=(rhs);
         last = rhs.last;
      }

      return *this;
   }  
   /// destructor
   virtual ~SPxDevexPR()
   {}
   /// clone function for polymorphism
   inline virtual SPxPricer* clone()  const 
   {
      return new SPxDevexPR(*this);
   }
   //@}

   //-------------------------------------
   /**@name Access / modification */
   //@{
   /// sets the solver
   virtual void load(SPxSolver* base);
   /// set entering/leaving algorithm
   virtual void setType(SPxSolver::Type);
   /// set row/column representation
   virtual void setRep(SPxSolver::Representation);
   ///
   virtual int selectLeave();
   ///
   virtual SPxId selectEnter();
   ///
   virtual void left4(int n, SPxId id);
   ///
   virtual void entered4(SPxId id, int n);
   /// \p n vectors have been added to loaded LP.
   virtual void addedVecs (int n);
   /// \p n covectors have been added to loaded LP.
   virtual void addedCoVecs(int n);
   //@}

   //-------------------------------------
   /**@name Consistency check */
   //@{
   /// consistency check
   virtual bool isConsistent() const;
   //@}
};

} // namespace soplex
#endif // _SPXDEVEXPR_H_
