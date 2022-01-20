/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996      Roland Wunderling                              */
/*                  1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  solbase.h
 * @brief Class for storing a primal-dual solution with basis information
 */
#ifndef _SOLBASE_H_
#define _SOLBASE_H_

/* undefine SOPLEX_DEBUG flag from including files; if SOPLEX_DEBUG should be defined in this file, do so below */
#ifdef SOPLEX_DEBUG
#define SOPLEX_DEBUG_SOLBASE
#undef SOPLEX_DEBUG
#endif

#include <assert.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "basevectors.h"
#include "spxsolver.h" // needed for basis information

namespace soplex
{
/**@class   SolBase
 * @brief   Class for storing a primal-dual solution with basis information
 * @ingroup Algo
 */
template< class R >
class SolBase
{
   friend class SoPlex;
   template < class S > friend class SolBase;

public:
   /// is the stored solution primal feasible?
   bool isPrimalFeasible() const
   {
      return _isPrimalFeasible;
   }

   /// gets the primal solution vector; returns true on success
   bool getPrimal(VectorBase<R>& vector) const
   {
      vector = _primal;

      return _isPrimalFeasible;
   }

   /// gets the vector of slack values; returns true on success
   bool getSlacks(VectorBase<R>& vector) const
   {
      vector = _slacks;

      return _isPrimalFeasible;
   }

   /// is a primal unbounded ray available?
   bool hasPrimalRay() const
   {
      return _hasPrimalRay;
   }

   /// gets the primal unbounded ray if available; returns true on success
   bool getPrimalRay(VectorBase<R>& vector) const
   {
      if( _hasPrimalRay )
         vector = _primalRay;

      return _hasPrimalRay;
   }

   /// is a dual solution available?
   bool isDualFeasible() const
   {
      return _isDualFeasible;
   }

   /// gets the dual solution vector; returns true on success
   bool getDual(VectorBase<R>& vector) const
   {
      vector = _dual;

      return _isDualFeasible;
   }

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedCost(VectorBase<R>& vector) const
   {
      vector = _redCost;

      return _isDualFeasible;
   }

   /// is a dual farkas ray available?
   bool hasDualFarkas() const
   {
      return _hasDualFarkas;
   }

   /// gets the Farkas proof if available; returns true on success
   bool getDualFarkas(VectorBase<R>& vector) const
   {
      if( _hasDualFarkas )
         vector = _dualFarkas;

      return _hasDualFarkas;
   }

   /// returns total size of primal solution
   int totalSizePrimal(const int base = 2) const
   {
      int size = 0;

      if( _isPrimalFeasible )
         size += totalSizeRational(_primal.get_const_ptr(), _primal.dim(), base);

      if( _hasPrimalRay )
         size += totalSizeRational(_primalRay.get_const_ptr(), _primalRay.dim(), base);

      return size;
   }

   /// returns total size of dual solution
   int totalSizeDual(const int base = 2) const
   {
      int size = 0;

      if( _isDualFeasible )
         size += totalSizeRational(_dual.get_const_ptr(), _dual.dim(), base);

      if( _hasDualFarkas )
         size += totalSizeRational(_dualFarkas.get_const_ptr(), _dualFarkas.dim(), base);

      return size;
   }

   /// returns size of least common multiple of denominators in primal solution
   int dlcmSizePrimal(const int base = 2) const
   {
      int size = 0;

      if( _isPrimalFeasible )
         size += dlcmSizeRational(_primal.get_const_ptr(), _primal.dim(), base);

      if( _hasPrimalRay )
         size += dlcmSizeRational(_primalRay.get_const_ptr(), _primalRay.dim(), base);

      return size;
   }

   /// returns  size of least common multiple of denominators in dual solution
   int dlcmSizeDual(const int base = 2) const
   {
      int size = 0;

      if( _isDualFeasible )
         size += dlcmSizeRational(_dual.get_const_ptr(), _dual.dim(), base);

      if( _hasDualFarkas )
         size += dlcmSizeRational(_dualFarkas.get_const_ptr(), _dualFarkas.dim(), base);

      return size;
   }

   /// returns size of largest denominator in primal solution
   int dmaxSizePrimal(const int base = 2) const
   {
      int size = 0;

      if( _isPrimalFeasible )
         size += dmaxSizeRational(_primal.get_const_ptr(), _primal.dim(), base);

      if( _hasPrimalRay )
         size += dmaxSizeRational(_primalRay.get_const_ptr(), _primalRay.dim(), base);

      return size;
   }

   /// returns size of largest denominator in dual solution
   int dmaxSizeDual(const int base = 2) const
   {
      int size = 0;

      if( _isDualFeasible )
         size += dmaxSizeRational(_dual.get_const_ptr(), _dual.dim(), base);

      if( _hasDualFarkas )
         size += dmaxSizeRational(_dualFarkas.get_const_ptr(), _dualFarkas.dim(), base);

      return size;
   }

   /// invalidate solution
   void invalidate()
   {
      _isPrimalFeasible = false;
      _hasPrimalRay = false;
      _isDualFeasible = false;
      _hasDualFarkas = false;
   }

private:
   DVectorBase<R> _primal;
   DVectorBase<R> _slacks;
   DVectorBase<R> _primalRay;
   DVectorBase<R> _dual;
   DVectorBase<R> _redCost;
   DVectorBase<R> _dualFarkas;

   R _objVal;

   unsigned int _isPrimalFeasible:1;
   unsigned int _hasPrimalRay:1;
   unsigned int _isDualFeasible:1;
   unsigned int _hasDualFarkas:1;

   /// default constructor only for friends
   SolBase<R>()
      : _objVal(0)
   {
      invalidate();
   }

   /// assignment operator only for friends
   SolBase<R>& operator=(const SolBase<R>& sol)
   {
      if( this != &sol )
      {

         _isPrimalFeasible = sol._isPrimalFeasible;
         _primal = sol._primal;
         _slacks = sol._slacks;
         _objVal = sol._objVal;

         _hasPrimalRay = sol._hasPrimalRay;
         if( _hasPrimalRay )
            _primalRay = sol._primalRay;

         _isDualFeasible = sol._isDualFeasible;
         _dual = sol._dual;
         _redCost = sol._redCost;

         _hasDualFarkas = sol._hasDualFarkas;
         if( _hasDualFarkas )
            _dualFarkas = sol._dualFarkas;
      }

      return *this;
   }

   /// assignment operator only for friends
   template < class S >
   SolBase<R>& operator=(const SolBase<S>& sol)
   {
      if( (SolBase<S>*)this != &sol )
      {

         _isPrimalFeasible = sol._isPrimalFeasible;
         _primal = sol._primal;
         _slacks = sol._slacks;
         _objVal = R(sol._objVal);

         _hasPrimalRay = sol._hasPrimalRay;
         if( _hasPrimalRay )
            _primalRay = sol._primalRay;

         _isDualFeasible = sol._isDualFeasible;
         _dual = sol._dual;
         _redCost = sol._redCost;

         _hasDualFarkas = sol._hasDualFarkas;
         if( _hasDualFarkas )
            _dualFarkas = sol._dualFarkas;
      }

      return *this;
   }
};
} // namespace soplex

/* reset the SOPLEX_DEBUG flag to its original value */
#undef SOPLEX_DEBUG
#ifdef SOPLEX_DEBUG_SOLBASE
#define SOPLEX_DEBUG
#undef SOPLEX_DEBUG_SOLBASE
#endif

#endif // _SOLBASE_H_
