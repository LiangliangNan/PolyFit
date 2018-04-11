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

/**@file  lprowsetbase.h
 * @brief Set of LP columns.
 */
#ifndef _LPROWSETBASE_H_
#define _LPROWSETBASE_H_


#include <assert.h>

#include "spxdefines.h"
#include "basevectors.h"
#include "datakey.h"
#include "lprowbase.h"

namespace soplex
{
/**@brief   Set of LP rows.
 * @ingroup Algebra
 *
 *  Class LPRowSetBase implements a set of \ref LPRowBase "LPRowBase%s". Unless for memory limitations, any number of
 *  LPRowBase%s may be #add%ed to an LPRowSetBase. Single or multiple LPRowBase%s may be added to an LPRowSetBase, where
 *  each method add() comes with two different signatures. One with and one without a parameter, used for returning the
 *  Keys assigned to the new LPRowBase%s by the set. See DataKey for a more detailed description of the concept of
 *  keys. For the concept of renumbering LPRowBase%s within an LPRowSetBase after removal of some LPRows see DataSet.
 *
 * @see        DataSet, DataKey
*/
template < class R >
class LPRowSetBase : protected SVSetBase<R>
{
   template < class S > friend class LPRowSetBase;

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   DVectorBase<R> left;    ///< vector of left hand sides (lower bounds) of LPRowBase%s.
   DVectorBase<R> right;   ///< vector of right hand sides (upper bounds) of LPRowBase%s.
   DVectorBase<R> object;  ///< vector of objective coefficients.

   //@}

protected:

   DataArray < int > scaleExp;   ///< row scaling factors (stored as bitshift)

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Helpers */
   //@{

   /// Returns the complete SVSet.
   const SVSetBase<R>* rowSet() const
   {
      return this;
   }

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Access / modification */
   //@{

   /// Returns the number of LPRowBase%s in LPRowSetBase.
   int num() const
   {
      return SVSetBase<R>::num();
   }

   /// Returns the maximum number of LPRowBase%s that fit.
   int max() const
   {
      return SVSetBase<R>::max();
   }

   /// Returns the vector of lhs values.
   const VectorBase<R>& lhs() const
   {
      return left;
   }

   /// Returns the vector of lhs values.
   VectorBase<R>& lhs_w()
   {
      return left;
   }

   /// Returns the lhs of the \p i 'th LPRowBase.
   const R& lhs(int i) const
   {
      return left[i];
   }

   /// Returns the lhs of the \p i 'th LPRowBase.
   R& lhs_w(int i)
   {
      return left[i];
   }

   /// Returns the lhs of the LPRowBase with DataKey \p k in LPRowSetBase.
   const R& lhs(const DataKey& k) const
   {
      return left[number(k)];
   }

   /// Returns the lhs of the LPRowBase with DataKey \p k in LPRowSetBase.
   R& lhs_w(const DataKey& k)
   {
      return left[number(k)];
   }

   /// Returns the vector of rhs values.
   const VectorBase<R>& rhs() const
   {
      return right;
   }

   /// Returns the vector of rhs values (writeable).
   VectorBase<R>& rhs_w()
   {
      return right;
   }

   /// Returns the rhs of the \p i 'th LPRowBase.
   const R& rhs(int i) const
   {
      return right[i];
   }

   /// Returns the rhs of the \p i 'th LPRowBase (writeable).
   R& rhs_w(int i)
   {
      return right[i];
   }

   /// Returns the rhs of the LPRowBase with DataKey \p k in LPRowSetBase.
   const R& rhs(const DataKey& k) const
   {
      return right[number(k)];
   }

   /// Returns the rhs of the LPRowBase with DataKey \p k in LPRowSetBase (writeable).
   R& rhs_w(const DataKey& k)
   {
      return right[number(k)];
   }

   /// Returns the vector of objective coefficients.
   const VectorBase<R>& obj() const
   {
      return object;
   }

   /// Returns the vector of objective coefficients (writeable).
   VectorBase<R>& obj_w()
   {
      return object;
   }

   /// Returns the objective coefficient of the \p i 'th LPRowBase.
   const R& obj(int i) const
   {
      return object[i];
   }

   /// Returns the objective coefficient of the \p i 'th LPRowBase (writeable).
   R& obj_w(int i)
   {
      return object[i];
   }

   /// Returns the objective coefficient of the LPRowBase with DataKey \p k in LPRowSetBase.
   const R& obj(const DataKey& k) const
   {
      return object[number(k)];
   }

   /// Returns the objective coefficient of the LPRowBase with DataKey \p k in LPRowSetBase (writeable).
   R& obj_w(const DataKey& k)
   {
      return object[number(k)];
   }

   /// Returns a writable rowVector of the \p i 'th LPRowBase.
   SVectorBase<R>& rowVector_w(int i)
   {
      return SVSetBase<R>::operator[](i);
   }

   /// Returns the rowVector of the \p i 'th LPRowBase.
   const SVectorBase<R>& rowVector(int i) const
   {
      return SVSetBase<R>::operator[](i);
   }

   /// Returns a writable rowVector of the LPRowBase with DataKey \p k.
   SVectorBase<R>& rowVector_w(const DataKey& k)
   {
      return SVSetBase<R>::operator[](k);
   }

   /// Returns the rowVector of the LPRowBase with DataKey \p k.
   const SVectorBase<R>& rowVector(const DataKey& k) const
   {
      return SVSetBase<R>::operator[](k);
   }

   /// Returns the inequalitiy type of the \p i 'th LPRowBase.
   typename LPRowBase<R>::Type type(int i) const
   {
      if( rhs(i) >= double(infinity) )
         return LPRowBase<R>::GREATER_EQUAL;
      if( lhs(i) <= double(-infinity) )
         return LPRowBase<R>::LESS_EQUAL;
      if( lhs(i) == rhs(i) )
         return LPRowBase<R>::EQUAL;

      return LPRowBase<R>::RANGE;
   }

   /// Returns the inequality type of the LPRowBase with DataKey \p k.
   typename LPRowBase<R>::Type type(const DataKey& k) const
   {
      return type(number(k));
   }

   /// Changes the inequality type of row \p i to \p type.
   void setType(int i, typename LPRowBase<R>::Type t)
   {
      switch( t )
      {
      case LPRowBase<R>::LESS_EQUAL:
         lhs_w(i) = -infinity;
         break;
      case LPRowBase<R>::EQUAL:
         if( lhs_w(i) > -infinity )
            rhs_w(i) = lhs(i);
         else
            lhs_w(i) = rhs(i);
         break;
      case LPRowBase<R>::GREATER_EQUAL:
         rhs_w(i) = infinity;
         break;
      case LPRowBase<R>::RANGE:
         MSG_ERROR( std::cerr << "EROWST01 RANGE not supported in LPRowSet::setType()" << std::endl );
         throw SPxInternalCodeException("XROWST01 This should never happen.");
      default:
         throw SPxInternalCodeException("XROWST02 This should never happen.");
      }
   }

   /// Returns the value of the \p i'th LPRowBase.
   const R& value(int i) const
   {
      if( rhs(i) < infinity )
         return rhs(i);
      else
      {
         assert(lhs(i) > -infinity);
         return lhs(i);
      }
   }

   /// Returns the value of the LPRowBase with DataKey \p k.
   /** The \em value of a row depends on its type: if the inequality is of type "greater or equal", the value is the lhs
    *  of the row. Otherwise, the value is the rhs.
   */
   const R& value(const DataKey& k) const
   {
      return value(number(k));
   }

   /// Returns the DataKey of the \p i 'th LPRowBase in LPRowSetBase.
   DataKey key(int i) const
   {
      return SVSetBase<R>::key(i);
   }

   /// Returns the number of the LPRowBase with DataKey \p k in LPRowSetBase.
   int number(const DataKey& k) const
   {
      return SVSetBase<R>::number(k);
   }

   /// does DataKey \p k belong to LPRowSetBase ?
   bool has(const DataKey& k) const
   {
      return SVSetBase<R>::has(k);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Extension
    *
    *  Extension methods come with two signatures, one of them providing a parameter to return the assigned
    *  DataKey(s). See DataSet for a more detailed description. All extension methods will automatically rearrange or
    *  allocate more memory if required.
   */
   //@{

   ///
   void add(const LPRowBase<R>& row)
   {
      DataKey k;
      add(k, row);
   }

   /// Adds \p row to LPRowSetBase.
   void add(DataKey& pkey, const LPRowBase<R>& prow)
   {
      add(pkey, prow.lhs(), prow.rowVector(), prow.rhs(), prow.obj());
   }

   /// Adds LPRowBase consisting of left hand side \p lhs, row vector \p rowVector, and right hand side \p rhs to LPRowSetBase.
   void add(const R& plhs, const SVectorBase<R>& prowVector, const R& prhs, const R& pobj = 0, const int& pscaleExp = 0)
   {
      DataKey k;
      add(k, plhs, prowVector, prhs, pobj, pscaleExp);
   }

   /// Adds LPRowBase consisting of left hand side \p lhs, row vector \p rowVector, and right hand side \p rhs to LPRowSetBase.
   template < class S >
   void add(const S* lhsValue, const S* rowValues, const int* rowIndices, int rowSize, const S* rhsValue, const S* objValue = 0)
   {
      assert(lhsValue != 0);
      assert(rowSize <= 0 || rowValues != 0);
      assert(rowSize <= 0 || rowIndices != 0);
      assert(rhsValue != 0);

      DataKey k;
      add(k, lhsValue, rowValues, rowIndices, rowSize, rhsValue, objValue);
   }

   /// Adds LPRowBase consisting of left hand side \p lhs, row vector \p rowVector, and right hand side \p rhs to
   /// LPRowSetBase, with DataKey \p key.
   template < class S >
   void add(DataKey& newkey, const S* lhsValue, const S* rowValues, const int* rowIndices, int rowSize, const S* rhsValue, const S* objValue = 0)
   {
      assert(lhsValue != 0);
      assert(rowSize <= 0 || rowValues != 0);
      assert(rowSize <= 0 || rowIndices != 0);
      assert(rhsValue != 0);

      SVSetBase<R>::add(newkey, rowValues, rowIndices, rowSize);

      if( num() > left.dim() )
      {
         left.reDim(num());
         right.reDim(num());
         object.reDim(num());
      }

      left[num() - 1] = *lhsValue;
      right[num() - 1] = *rhsValue;
      if( objValue != 0 )
         object[num() - 1] = *objValue;
      else
         object[num() - 1] = 0;
   }

   /// Adds LPRowBase consisting of left hand side \p lhs, row vector \p rowVector, and right hand side \p rhs to
   /// LPRowSetBase, with DataKey \p key.
   void add(DataKey& newkey, const R& newlhs, const SVectorBase<R>& newrowVector, const R& newrhs, const R& newobj = 0, const int& newscaleExp = 0)
   {
      SVSetBase<R>::add(newkey, newrowVector);

      if( num() > left.dim() )
      {
         left.reDim(num());
         right.reDim(num());
         object.reDim(num());
         scaleExp.reSize(num());
      }

      left[num() - 1] = newlhs;
      right[num() - 1] = newrhs;
      object[num() - 1] = newobj;
      scaleExp[num() - 1] = newscaleExp;
   }

   ///
   void add(const LPRowSetBase<R>& newset)
   {
      int i = num();

      SVSetBase<R>::add(newset);

      if( num() > left.dim() )
      {
         left.reDim(num());
         right.reDim(num());
         object.reDim(num());
         scaleExp.reSize(num());
      }

      for( int j = 0; i < num(); ++i, ++j )
      {
         left[i] = newset.lhs(j);
         right[i] = newset.rhs(j);
         object[i] = newset.obj(j);
         scaleExp[i] = newset.scaleExp[j];
      }
   }

   /// Adds all LPRowBase%s of \p set to LPRowSetBase.
   void add(DataKey keys[], const LPRowSetBase<R>& set)
   {
      int i = num();

      add(set);

      for( int j = 0; i < num(); ++i, ++j )
         keys[j] = key(i);
   }

   /// Extends row \p n to fit \p newmax nonzeros.
   void xtend(int n, int newmax)
   {
      SVSetBase<R>::xtend(rowVector_w(n), newmax);
   }

   /// Extends row with DataKey \p key to fit \p newmax nonzeros.
   void xtend(const DataKey& pkey, int pnewmax)
   {
      SVSetBase<R>::xtend(rowVector_w(pkey), pnewmax);
   }

   /// Adds \p n nonzero (\p idx, \p val)-pairs to rowVector with DataKey \p k.
   void add2(const DataKey& k, int n, const int idx[], const R val[])
   {
      SVSetBase<R>::add2(rowVector_w(k), n, idx, val);
   }

   /// Adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th rowVector.
   void add2(int i, int n, const int idx[], const R val[])
   {
      SVSetBase<R>::add2(rowVector_w(i), n, idx, val);
   }

   /// Adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th rowVector.
   template < class S >
   void add2(int i, int n, const int idx[], const S val[])
   {
      SVSetBase<R>::add2(rowVector_w(i), n, idx, val);
   }

   /// Creates new LPRowBase with specified parameters and returns a reference to its row vector.
   SVectorBase<R>& create(int pnonzeros = 0, const R& plhs = 0, const R& prhs = 1, const R& pobj = 0, const int& pscaleExp = 0)
   {
      DataKey k;
      return create(k, pnonzeros, plhs, prhs, pobj, pscaleExp);
   }

   /// Creates new LPRowBase with specified parameters and returns a reference to its row vector.
   SVectorBase<R>& create(DataKey& newkey, int nonzeros = 0, const R& newlhs = 0, const R& newrhs = 1, const R& newobj = 0, const int& newscaleExp = 0)
   {
      if( num() + 1 > left.dim() )
      {
         left.reDim(num() + 1);
         right.reDim(num() + 1);
         object.reDim(num() + 1);
         scaleExp.reSize(num() + 1);
      }

      left[num()] = newlhs;
      right[num()] = newrhs;
      object[num()] = newobj;
      scaleExp[num()] = newscaleExp;

      return *SVSetBase<R>::create(newkey, nonzeros);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Shrinking
    *
    *  See DataSet for a description of the renumbering of the remaining LPRowBase%s in a LPRowSetBase after the call of
    *  a removal method.
    */
   //@{

   /// Removes \p i 'th LPRowBase.
   void remove(int i)
   {
      SVSetBase<R>::remove(i);
      left[i] = left[num()];
      right[i] = right[num()];
      object[i] = object[num()];
      scaleExp[i] = scaleExp[num()];
      left.reDim(num());
      right.reDim(num());
      object.reDim(num());
      scaleExp.reSize(num());
   }

   /// Removes LPRowBase with DataKey \p k.
   void remove(const DataKey& k)
   {
      remove(number(k));
   }

   /// Removes multiple LPRowBase%s.
   void remove(int perm[])
   {
      int j = num();

      SVSetBase<R>::remove(perm);

      for( int i = 0; i < j; ++i )
      {
         if( perm[i] >= 0 && perm[i] != i )
         {
            left[perm[i]] = left[i];
            right[perm[i]] = right[i];
            object[perm[i]] = object[i];
            scaleExp[perm[i]] = scaleExp[i];
         }
      }

      left.reDim (num());
      right.reDim(num());
      object.reDim(num());
      scaleExp.reSize(num());
   }

   /// Removes \p n LPRowBase%s with row numbers given by \p nums.
   void remove(const int nums[], int n)
   {
      DataArray<int> perm(num());
      remove(nums, n, perm.get_ptr());
   }

   /// Removes \p n LPRowBase%s with row numbers given by \p nums,
   /// Stores permutation of row indices in \p perm.
   void remove(const int nums[], int n, int* perm)
   {
      SVSetBase<R>::remove(nums, n, perm);

      int j = num();

      for( int i = 0; i < j; ++i )
      {
         if( perm[i] >= 0 && perm[i] != i )
         {
            left[perm[i]] = left[i];
            right[perm[i]] = right[i];
            object[perm[i]] = object[i];
            scaleExp[perm[i]] = scaleExp[i];
         }
      }

      left.reDim (num());
      right.reDim(num());
      object.reDim(num());
      scaleExp.reSize(num());
   }

   /// Removes all LPRowBase%s.
   void clear()
   {
      SVSetBase<R>::clear();
      left.reDim(num());
      right.reDim(num());
      object.reDim(num());
      scaleExp.clear();
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Memory Management
    *
    *  For a description of the memory management methods, see the documentation of SVSet, which has been used for
    *  implementating LPRowSetBase.
    */
   //@{

   /// Reallocates memory to be able to store \p newmax LPRowBase%s.
   void reMax(int newmax = 0)
   {
      SVSetBase<R>::reMax(newmax);
      left.reSize (max());
      right.reSize(max());
      object.reSize(max());
      scaleExp.reSize(max());
   }

   /// Returns number of used nonzero entries.
   int memSize() const
   {
      return SVSetBase<R>::memSize();
   }

   /// Returns length of nonzero memory.
   int memMax() const
   {
      return SVSetBase<R>::memMax();
   }

   /// Reallocates memory to be able to store \p newmax nonzeros.
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
   /**@name Consistency check */

   /// Checks consistency.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      const int ldim = left.dim();

      if( ldim != right.dim() )
         return MSGinconsistent("LPRowSetBase");
      if( ldim != object.dim() )
         return MSGinconsistent("LPRowSetBase");
      if( ldim != num() )
         return MSGinconsistent("LPRowSetBase");

      return left.isConsistent() && right.isConsistent() && object.isConsistent() && SVSetBase<R>::isConsistent();
#else
      return true;
#endif
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Construction / Destruction */
   //@{

   /// Default constructor.
   /** The user can specify the initial maximum number of rows \p max and the initial maximum number of nonzero entries
    *  \p memmax. If these parameters are omitted, a default size is used. However, one can add an arbitrary number of
    *  rows to the LPRowSetBase, which may result in automated memory realllocation.
    */
   explicit
   LPRowSetBase<R>(int pmax = -1, int pmemmax = -1)
      : SVSetBase<R>(pmax, pmemmax), left(0), right(0), object(0), scaleExp(0)
   {
      assert(isConsistent());
   }

   /// Assignment operator.
   LPRowSetBase<R>& operator=(const LPRowSetBase<R>& rs)
   {
      if( this != &rs )
      {
         SVSetBase<R>::operator=(rs);
         left = rs.left;
         right = rs.right;
         object = rs.object;
         scaleExp = rs.scaleExp;

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   LPRowSetBase<R>& operator=(const LPRowSetBase<S>& rs)
   {
      if( this != (const LPRowSetBase<R>*)(&rs) )
      {
         SVSetBase<R>::operator=(rs);
         left = rs.left;
         right = rs.right;
         object = rs.object;
         scaleExp = rs.scaleExp;

         assert(isConsistent());
      }

      return *this;
   }

   /// Copy constructor.
   LPRowSetBase<R>(const LPRowSetBase<R>& rs)
      : SVSetBase<R>(rs)
      , left(rs.left)
      , right(rs.right)
      , object(rs.object)
      , scaleExp(rs.scaleExp)
   {
      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   LPRowSetBase<R>(const LPRowSetBase<S>& rs)
      : SVSetBase<R>(rs)
      , left(rs.left)
      , right(rs.right)
      , object(rs.object)
      , scaleExp(rs.scaleExp)
   {
      assert(isConsistent());
   }

   /// Destructor.
   virtual ~LPRowSetBase<R>()
   {}

   //@}
};
} // namespace soplex
#endif // _LPROWSETBASE_H_
