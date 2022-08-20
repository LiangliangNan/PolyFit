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

/**@file  idxset.h
 * @brief Set of indices.
 */
#ifndef _IDXSET_H_
#define _IDXSET_H_

#include "spxdefines.h"
#include "spxalloc.h"
#include <assert.h>

namespace soplex
{
/**@brief   Set of indices.
   @ingroup Elementary

   Class IdxSet provides a set of indices. At construction it must be given
   an array of int where to store the indice and its length. The array will
   from then on be managed by the IdxSet.
   
   Indices are implicitely numbered from 0 thru size()-1. They can be
   accessed (and altered) via method index() with the desired index number as
   argument.  Range checking is performed in the debug version.
   
   Indices may be added or removed from the set, by calling add() or
   remove() methods, respectively. However, no IdxSet can hold more then
   max() indices, i.e. the number given at the constructor.
   
   When removing indices, the remaining ones are renumbered. However, all
   indices before the first removed index keep their number unchanged.

   The internal structure of an IdxSet consists of an array #idx storing the
   indices, its length len, and the actually used number of indices #num.
   The class IdxSet doesn't allocate memory for the #idx array. Instead, the
   user has to provide an adequate buffer to the constructor.

   An IdxSet cannot be extended to fit more than max() elements. If
   necessary, the user must explicitely provide the IdxSet with a
   suitable memory. Alternatively, one can use \ref DIdxSet "DIdxSets" 
   which provide the required memory managemant.
*/
class IdxSet
{
protected:

   //---------------------------------------
   /**@name Data */
   //@{
   int  num;           ///< number of used indices
   int  len;           ///< length of array \ref soplex::IdxSet::idx "idx"
   int* idx;           ///< array of indices
   bool freeArray;     ///< true iff \ref soplex::IdxSet::idx "idx" should be freed inside of this object
   //@}

public:

   //---------------------------------------
   /**@name Construction / destruction */
   //@{
   /// constructor.
   /** The constructur receives the index memory \p imem to use for saving
       its indices. This must be large enough to fit \p n indices. \p l can
       be given to construct an #IdxSet initialized to the \p l first
       indices in \p imem.
    */
   IdxSet(int n, int imem[], int l = 0)
      : num(l), len(n), idx(imem), freeArray(false)
   {
      assert(isConsistent());
   }

   /// default constructor.
   /** The default constructor creates an index set with an empty index
       space. You cannot store any indices in an #IdxSet created with
       the default constructor.
   */
   IdxSet()
      : num(0), len(0), idx(0), freeArray(false)
   {
      assert(isConsistent());
   }
   
   /// destructor.
   virtual ~IdxSet()
   {
      if(freeArray)
         spx_free(idx);
   }

   /// assignment operator.
   /** The assignment operator copies all nonzeros of the right handside
       #IdxSet to the left one. This implies, that the latter must have
       enough index memory.
    */
   IdxSet& operator=(const IdxSet& set);
   /// copy constructor.
   IdxSet(const IdxSet&);
   //@}

   //---------------------------------------
   /**@name Access */
   //@{
   /// access \p n 'th index.
   int index(int n) const
   {
      assert(n >= 0 && n < size() && idx != 0);
      return idx[n];
   }
   /// returns the number of used indices.
   int size() const
   {
      return num;
   }
   /// returns the maximal number of indices which can be stored in IdxSet.
   int max() const
   {
      return len;
   }

   /// returns the maximal index.
   int dim() const;

   /// returns the position of index \p i.
   /** Finds the position of the first index \p i in the #IdxSet. If no index \p i is
       available in the #IdxSet, -1 is returned. Otherwise,
       index(pos(\p i)) == \p i holds.
    */
   int pos(int i) const;
   //@}

   //---------------------------------------
   /**@name Modification */
   //@{
   /// appends \p n uninitialized indices.
   void add(int n)
   {
      assert(n >= 0 && n + size() <= max());
      num += n;
   }

   /// appends all indices of \p set.
   void add(const IdxSet& set)
   {
      add(set.size(), set.idx);
   }

   /// appends \p n indices in \p i.
   void add(int n, const int i[]);

   /// appends index \p i.
   void addIdx(int i)
   {
      assert(size() < max());
      idx[num++] = i;
   }
   /// removes indices at position numbers \p n through \p m.
   void remove(int n, int m);

   /// removes \p n 'th index.
   void remove(int n)
   {
//      /**@todo Shouldn't this be an assert instead of an if (see add()) */
//      if (n < size() && n >= 0)
//         idx[n] = idx[--num];
      assert(n >= 0 && n < size()); 
      idx[n] = idx[--num];
   }

   /// removes all indices.
   void clear()
   {
      num = 0;
   }
   //@}

   //---------------------------------------
   /**@name Consistency check */
   //@{
   /// consistency check.
   bool isConsistent() const;
   //@}
};

} // namespace soplex
#endif // _IDXSET_H_
