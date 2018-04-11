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

#include "didxset.h"
#include "spxalloc.h"

namespace soplex
{

void DIdxSet::setMax(int newmax)
{ 
   assert(idx   != 0);
   assert(max() >  0);

   len = (newmax < size()) ? size() : newmax;
   len = (len < 1) ? 1 : len;

   assert(len > 0);

   spx_realloc(idx, len);
}

DIdxSet::DIdxSet(const IdxSet& old)
   : IdxSet()
{
   len = old.size();
   len = (len < 1) ? 1 : len;
   spx_alloc(idx, len);

   IdxSet::operator= ( old );
}

DIdxSet::DIdxSet(const DIdxSet& old)
   : IdxSet()
{
   len = old.size();
   len = (len < 1) ? 1 : len;
   spx_alloc(idx, len);

   IdxSet::operator= ( old );
}

DIdxSet::DIdxSet(int n)
   : IdxSet()
{
   len = (n < 1) ? 1 : n;
   spx_alloc(idx, len);
}

DIdxSet::~DIdxSet()
{
   if(idx)   
      spx_free(idx);
}
} // namespace soplex
