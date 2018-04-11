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

#include "updatevector.h"

namespace soplex
{

UpdateVector& UpdateVector::operator=(const UpdateVector& rhs)
{
   if (this != &rhs)
   {
      theval   = rhs.theval;
      thedelta = rhs.thedelta;
      DVector::operator=(rhs);

      assert(UpdateVector::isConsistent());
   }
   return *this;
}

UpdateVector::UpdateVector( const UpdateVector& base)
   : DVector( base )
   , theval(base.theval)
   , thedelta(base.thedelta)
{
   assert(UpdateVector::isConsistent());
}

bool UpdateVector::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (dim() != thedelta.dim())
      return MSGinconsistent("UpdateVector");

   return DVector::isConsistent() && thedelta.isConsistent();
#else
   return true;
#endif
}
} // namespace soplex
