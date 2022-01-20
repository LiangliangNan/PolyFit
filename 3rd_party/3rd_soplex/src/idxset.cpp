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

#include "idxset.h"

namespace soplex
{

int IdxSet::dim() const
{
   int ddim = -1;

   for (int i = 0; i < size(); i++)
      if (ddim < idx[i])
         ddim = idx[i];

   return ddim; 
}

int IdxSet::pos(int i) const
{
   for(int n = 0; n < size(); n++)
      if (idx[n] == i)
         return n;

   return -1;
}

void IdxSet::add(int n, const int i[])
{
   assert(n >= 0 && size() + n <= max());
   for (int j = 0; j < n; j++)
      idx[size() + j] = i[j];
   add(n);
}

void IdxSet::remove(int n, int m)
{
   assert(n <= m && m < size() && n >= 0);
   ++m;

   int cpy = m - n;
   int newnum = num - cpy;
   cpy = (size() - m >= cpy) ? cpy : size() - m;

   do
   {
      --num;
      --cpy;
      idx[n + cpy] = idx[num];
   }
   while (cpy > 0);
   num = newnum;
}

IdxSet& IdxSet::operator=(const IdxSet& rhs)
{
   if (this != &rhs)
   {
      if (idx != 0 && max() < rhs.size())
      {
         if (freeArray)
            spx_free(idx);
         idx = 0;
      }

      if (idx == 0)
      {
         len = rhs.size();
         spx_alloc(idx, len);
         freeArray = true;
      }

      for (num = 0; num < rhs.size(); ++num)
         idx[num] = rhs.idx[num];
   }

   assert(size() == rhs.size());
   assert(size() <= max());
   assert(isConsistent());

   return *this;
}

IdxSet::IdxSet(const IdxSet& old)
   : len(old.len)
   , idx(0)
{
   spx_alloc(idx, len);

   for (num = 0; num < old.num; num++)
      idx[num] = old.idx[num];

   freeArray = true;
   
   assert(size() == old.size());
   assert(size() <= max());
   assert(isConsistent());
}

bool IdxSet::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   int i, j;

   if (len > 0 && idx == 0)
      return MSGinconsistent("IdxSet");

   for (i = 0; i < size(); ++i)
   {
      if (index(i) < 0)
         return MSGinconsistent("IdxSet");

      for (j = 0; j < i; j++)
         if (index(i) == index(j))
            return MSGinconsistent("IdxSet");
   }
#endif

   return true;
}
} // namespace soplex
