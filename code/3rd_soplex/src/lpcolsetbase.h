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

/**@file  lpcolsetbase.h
 * @brief Set of LP columns.
 */
#ifndef _LPCOLSETBASE_H_
#define _LPCOLSETBASE_H_

#include <assert.h>

#include "spxdefines.h"
#include "basevectors.h"
#include "datakey.h"
#include "lpcolbase.h"

namespace soplex
{
/**@brief   Set of LP columns.
 * @ingroup Algebra
 *
 *  Class LPColSetBase implements a set of \ref LPColBase "LPColBase%s". Unless for memory limitations, any number of LPColBase%s may be
 *  #add%ed to an LPColSetBase. Single or multiple LPColBase%s may be #add%ed to an LPColSetBase, where each method add() comes with
 *  two different signatures. One with and one without a parameter, used for returning the \ref DataKey "DataKeys"
 *  assigned to the new LPColBase%s by the set. See DataKey for a more detailed description of the concept of keys. For the
 *  concept of renumbering LPColBase%s within an LPColSetBase after removal of some LPColBase%s, see DataSet.
 *
 * @see        DataSet, DataKey
 */
template < class R >
class LPColSetBase : protected SVSetBase<R>
{
   template < class S > friend class LPColSetBase;

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   DVectorBase<R> low;     ///< vector of lower bounds.
   DVectorBase<R> up;      ///< vector of upper bounds.
   DVectorBase<R> object;  ///< vector of objective coefficients.

   //@}

protected:

   DataArray < int > scaleExp;   ///< column scaling factors (stored as bitshift)

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Protected helpers */
   //@{

   /// Returns the complete SVSetBase.
   const SVSetBase<R>* colSet() const
   {
      return this;
   }

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Inquiry */
   //@{

   /// Returns the number of LPColBase%s currently in LPColSetBase.
   int num() const
   {
      return SVSetBase<R>::num();
   }

   /// Returns maximum number of LPColBase%s currently fitting into LPColSetBase.
   int max() const
   {
      return SVSetBase<R>::max();
   }

   ///
   const VectorBase<R>& maxObj() const
   {
      return object;
   }

   /// Returns vector of objective values w.r.t. maximization.
   VectorBase<R>& maxObj_w()
   {
      return object;
   }

   ///
   const R& maxObj(int i) const
   {
      return object[i];
   }

   /// Returns objective value (w.r.t. maximization) of \p i 'th LPColBase in LPColSetBase.
   R& maxObj_w(int i)
   {
      return object[i];
   }

   ///
   const R& maxObj(const DataKey& k) const
   {
      return object[number(k)];
   }

   /// Returns objective value (w.r.t. maximization) of LPColBase with DataKey \p k in LPColSetBase.
   R& maxObj_w(const DataKey& k)
   {
      return object[number(k)];
   }

   ///
   const VectorBase<R>& lower() const
   {
      return low;
   }

   /// Returns vector of lower bound values.
   VectorBase<R>& lower_w()
   {
      return low;
   }

   ///
   const R& lower(int i) const
   {
      return low[i];
   }

   /// Returns lower bound of \p i 'th LPColBase in LPColSetBase.
   R& lower_w(int i)
   {
      return low[i];
   }

   ///
   const R& lower(const DataKey& k) const
   {
      return low[number(k)];
   }

   /// Returns lower bound of LPColBase with DataKey \p k in LPColSetBase.
   R& lower_w(const DataKey& k)
   {
      return low[number(k)];
   }

   ///
   const VectorBase<R>& upper() const
   {
      return up;
   }

   /// Returns vector of upper bound values.
   VectorBase<R>& upper_w()
   {
      return up;
   }

   ///
   const R& upper(int i) const
   {
      return up[i];
   }

   /// Returns upper bound of \p i 'th LPColBase in LPColSetBase.
   R& upper_w(int i)
   {
      return up[i];
   }

   ///
   const R& upper(const DataKey& k) const
   {
      return up[number(k)];
   }

   /// Returns upper bound of LPColBase with DataKey \p k in LPColSetBase.
   R& upper_w(const DataKey& k)
   {
      return up[number(k)];
   }

   ///
   SVectorBase<R>& colVector_w(int i)
   {
      return SVSetBase<R>::operator[](i);
   }

   /// Returns colVector of \p i 'th LPColBase in LPColSetBase.
   const SVectorBase<R>& colVector(int i) const
   {
      return SVSetBase<R>::operator[](i);
   }

   /// Returns writeable colVector of LPColBase with DataKey \p k in LPColSetBase.
   SVectorBase<R>& colVector_w(const DataKey& k)
   {
      return SVSetBase<R>::operator[](k);
   }

   /// Returns colVector of LPColBase with DataKey \p k in LPColSetBase.
   const SVectorBase<R>& colVector(const DataKey& k) const
   {
      return SVSetBase<R>::operator[](k);
   }

   /// Returns DataKey of \p i 'th LPColBase in LPColSetBase.
   DataKey key(int i) const
   {
      return SVSetBase<R>::key(i);
   }

   /// Returns number of LPColBase with DataKey \p k in LPColSetBase.
   int number(const DataKey& k) const
   {
      return SVSetBase<R>::number(k);
   }

   /// Does DataKey \p k belong to LPColSetBase ?
   bool has(const DataKey& k) const
   {
      return SVSetBase<R>::has(k);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Extension
    *
    *  All extension methods come with two signatures, one of which providing a parameter to return the assigned
    *  DataKey(s). See DataSet for a more detailed description. All extension methods are designed to automatically
    *  reallocate memory if required.
    */
   //@{

   ///
   void add(const LPColBase<R>& pcol)
   {
      DataKey k;
      add(k, pcol);
   }

   /// Adds p pcol to LPColSetBase.
   void add(DataKey& pkey, const LPColBase<R>& pcol)
   {
      add(pkey, pcol.obj(), pcol.lower(), pcol.colVector(), pcol.upper());
   }

   ///
   void add(const R& pobj, const R& plower, const SVectorBase<R>& pcolVector, const R& pupper, const int& pscaleExp = 0)
   {
      DataKey k;
      add(k, pobj, plower, pcolVector, pupper, pscaleExp);
   }

   /// Adds LPColBase consisting of objective value \p obj, lower bound \p lower, column vector \p colVector and upper bound \p upper to LPColSetBase.
   void add(DataKey& newkey, const R& obj, const R& newlower, const SVectorBase<R>& newcolVector, const R& newupper, const int& newscaleExp = 0)
   {
      SVSetBase<R>::add(newkey, newcolVector);

      if( num() > low.dim() )
      {
         low.reDim(num());
         up.reDim(num());
         object.reDim(num());
         scaleExp.reSize(num());
      }

      low[num() - 1] = newlower;
      up[num() - 1] = newupper;
      object[num() - 1] = obj;
      scaleExp[num() - 1] = newscaleExp;
   }

   /// Adds LPColBase consisting of left hand side \p lhs, column vector \p colVector, and right hand side \p rhs to LPColSetBase.
   template < class S >
   void add(const S* obj, const S* lowerValue, const S* colValues, const int* colIndices, int colSize, const S* upperValue)
   {
      DataKey k;
      add(k, obj, lowerValue, colValues, colIndices, colSize, upperValue);
   }

   /// Adds LPColBase consisting of left hand side \p lhs, column vector \p colVector, and right hand side \p rhs to
   /// LPColSetBase, with DataKey \p key.
   template < class S >
   void add(DataKey& newkey, const S* objValue, const S* lowerValue, const S* colValues, const int* colIndices, int colSize, const S* upperValue)
   {
      SVSetBase<R>::add(newkey, colValues, colIndices, colSize);

      if( num() > low.dim() )
      {
         low.reDim(num());
         up.reDim(num());
         object.reDim(num());
      }

      low[num() - 1] = *lowerValue;
      up[num() - 1] = *upperValue;
      object[num() - 1] = *objValue;
   }

   ///
   void add(const LPColSetBase<R>& newset)
   {
      int i = num();

      SVSetBase<R>::add(newset);

      if( num() > low.dim() )
      {
         low.reDim(num());
         up.reDim(num());
         object.reDim(num());
         scaleExp.reSize(num());
      }

      for( int j = 0; i < num(); ++i, ++j )
      {
         low[i] = newset.lower(j);
         up[i] = newset.upper(j);
         object[i] = newset.maxObj(j);
         scaleExp[i] = newset.scaleExp[j];
      }
   }

   /// Adds all LPColBase%s of \p set to LPColSetBase.
   void add(DataKey keys[], const LPColSetBase<R>& newset)
   {
      int i = num();

      add(newset);

      for( int j = 0; i < num(); ++i, ++j )
         keys[j] = key(i);
   }

   /// Extends column \p n to fit \p newmax nonzeros.
   void xtend(int n, int newmax)
   {
      SVSetBase<R>::xtend(colVector_w(n), newmax);
   }

   /// Extends column with DataKey \p key to fit \p newmax nonzeros.
   void xtend(const DataKey& pkey, int pnewmax)
   {
      SVSetBase<R>::xtend(colVector_w(pkey), pnewmax);
   }

   ///
   void add2(const DataKey& k, int n, const int idx[], const R val[])
   {
      SVSetBase<R>::add2(colVector_w(k), n, idx, val);
   }

   /// Adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th colVector.
   void add2(int i, int n, const int idx[], const R val[])
   {
      SVSetBase<R>::add2(colVector_w(i), n, idx, val);
   }

   /// Adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th colVector.
   template < class S >
   void add2(int i, int n, const int idx[], const S val[])
   {
      SVSetBase<R>::add2(colVector_w(i), n, idx, val);
   }

   ///
   SVectorBase<R>& create(int pnonzeros = 0, const R& pobj = 1, const R& plw = 0, const R& pupp = 1, const int& pscaleExp = 0)
   {
      DataKey k;
      return create(k, pnonzeros, pobj, plw, pupp, pscaleExp);
   }

   /// Creates new LPColBase with specified arguments and returns a reference to its column vector.
   SVectorBase<R>& create(DataKey& newkey, int nonzeros = 0, const R& obj = 1, const R& newlow = 0, const R& newup = 1, const int& newscaleExp = 0)
   {
      if( num() + 1 > low.dim() )
      {
         low.reDim(num() + 1);
         up.reDim(num() + 1);
         object.reDim(num() + 1);
         scaleExp.reSize(num() + 1);
      }

      low[num()] = newlow;
      up[num()] = newup;
      object[num()] = obj;
      scaleExp[num()] = newscaleExp;

      return *SVSetBase<R>::create(newkey, nonzeros);
   }

   //@}


   // ------------------------------------------------------------------------------------------------------------------
   /**@name Shrinking
    *
    *  See DataSet for a description of the renumbering of the remaining LPColBase%s in a LPColSetBase after the call of
    *  a removal method.
    */
   //@{

   /// Removes \p i 'th LPColBase.
   void remove(int i)
   {
      SVSetBase<R>::remove(i);
      low[i] = low[num()];
      up[i] = up[num()];
      object[i] = object[num()];
      scaleExp[i] = scaleExp[num()];
      low.reDim(num());
      up.reDim(num());
      object.reDim(num());
      scaleExp.reSize(num());
   }

   /// Removes LPColBase with DataKey \p k.
   void remove(const DataKey& k)
   {
      remove(number(k));
   }

   /// Removes multiple elements.
   void remove(int perm[])
   {
      int n = num();

      SVSetBase<R>::remove(perm);

      for( int i = 0; i < n; ++i )
      {
         if( perm[i] >= 0 && perm[i] != i )
         {
            low[perm[i]] = low[i];
            up[perm[i]] = up[i];
            object[perm[i]] = object[i];
            scaleExp[perm[i]] = scaleExp[i];
         }
      }

      low.reDim(num());
      up.reDim(num());
      object.reDim(num());
      scaleExp.reSize(num());
   }

   /// Removes LPColBase%s with numbers \p nums, where \p n is the length of the array \p nums
   void remove(const int nums[], int n)
   {
      DataArray < int > perm(num());
      remove(nums, n, perm.get_ptr());
   }

   /// Removes LPColBase%s with numbers \p nums, where \p n is the length of the array \p nums, and stores the index permutation in array \p perm.
   void remove(const int nums[], int n, int* perm)
   {
      SVSetBase<R>::remove(nums, n, perm);

      int j = num();

      for( int i = 0; i < j; ++i )
      {
         if( perm[i] >= 0 && perm[i] != i )
         {
            low[perm[i]] = low[i];
            up[perm[i]] = up[i];
            object[perm[i]] = object[i];
            scaleExp[perm[i]] = scaleExp[i];
         }
      }

      low.reDim(num());
      up.reDim(num());
      object.reDim(num());
      scaleExp.reSize(num());
   }

   /// Removes all LPColBase%s from the set.
   void clear()
   {
      SVSetBase<R>::clear();
      low.reDim(num());
      up.reDim(num());
      object.reDim(num());
      scaleExp.clear();
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Memory Management
    *  See SVSet for a description of the memory management methods.
   */
   //@{

   /// Reallocates memory to be able to store \p newmax LPColBase%s.
   void reMax(int newmax = 0)
   {
      SVSetBase<R>::reMax(newmax);
      up.reSize(max());
      low.reSize(max());
      object.reSize(max());
      scaleExp.reSize(max());
   }

   /// Returns used nonzero memory.
   int memSize() const
   {
      return SVSetBase<R>::memSize();
   }

   /// Returns length of nonzero memory.
   int memMax() const
   {
      return SVSetBase<R>::memMax();
   }

   /// Resets length of nonzero memory.
   void memRemax(int newmax)
   {
      SVSetBase<R>::memRemax(newmax);
   }

   /// Garbage collection in nonzero memory.
   void memPack()
   {
      SVSetBase<R>::memPack();
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Miscellaneous */
   //@{

   /// Checks consistency.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if( low.dim() != object.dim() )
         return MSGinconsistent("LPColSetBase");
      if( low.dim() != up.dim() )
         return MSGinconsistent("LPColSetBase");
      if( low.dim() != num() )
         return MSGinconsistent("LPColSetBase");

      return low.isConsistent() && up.isConsistent() && SVSetBase<R>::isConsistent();
#else
      return true;
#endif
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Constructors / Destructors */
   //@{

   /// Default constructor.
   /** The user can specify the initial maximum number of columns \p max and the initial maximum number of nonzero
    *  entries \p memmax. If these parameters are omitted, a default size is used. However, one can add an arbitrary
    *  number of columns to the LPColSetBase, which may result in automated memory realllocation.
   */
   explicit
   LPColSetBase<R>(int pmax = -1, int pmemmax = -1)
      : SVSetBase<R>(pmax, pmemmax), low(0), up(0), object(0), scaleExp(0)
   {
      assert(isConsistent());
   }

   /// Assignment operator.
   LPColSetBase<R>& operator=(const LPColSetBase<R>& rs)
   {
      if( this != &rs )
      {
         SVSetBase<R>::operator=(rs);
         low = rs.low;
         up = rs.up;
         object = rs.object;
         scaleExp = rs.scaleExp;

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   LPColSetBase<R>& operator=(const LPColSetBase<S>& rs)
   {
      if( this != (const LPColSetBase<R>*)(&rs) )
      {
         SVSetBase<R>::operator=(rs);
         low = rs.low;
         up = rs.up;
         object = rs.object;
         scaleExp = rs.scaleExp;

         assert(isConsistent());
      }

      return *this;
   }

   /// Copy constructor.
   LPColSetBase<R>(const LPColSetBase<R>& rs)
      : SVSetBase<R>(rs)
      , low(rs.low)
      , up(rs.up)
      , object(rs.object)
      , scaleExp(rs.scaleExp)
   {
      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   LPColSetBase<R>(const LPColSetBase<S>& rs)
      : SVSetBase<R>(rs)
      , low(rs.low)
      , up(rs.up)
      , object(rs.object)
      , scaleExp(rs.scaleExp)
   {
      assert(isConsistent());
   }

   /// Destructor.
   virtual ~LPColSetBase<R>()
   {}

   //@}
};

} // namespace soplex
#endif // _LPCOLSETBASE_H_
