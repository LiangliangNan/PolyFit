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

#include <iostream>

#include "spxdefines.h"
#include "spxbasis.h"
#include "spxsolver.h"

namespace soplex
{

SPxBasis::Desc::Desc(const SPxSolver& base)
{
   reSize(base.nRows(), base.nCols());

   if (base.rep() == SPxSolver::ROW)
   {
      stat   = &rowstat;
      costat = &colstat;
   }
   else
   {
      assert(base.rep() == SPxSolver::COLUMN);

      stat   = &colstat;
      costat = &rowstat;
   }

   assert(Desc::isConsistent());
}

SPxBasis::Desc::Desc(const Desc& old)
   : rowstat(old.rowstat)
   , colstat(old.colstat)
{
   if (old.stat == &old.rowstat)
   {
      assert(old.costat == &old.colstat);
      
      stat   = &rowstat;
      costat = &colstat;
   }
   else
   {
      assert(old.costat == &old.rowstat);
      
      stat   = &colstat;
      costat = &rowstat;
   }

   assert(Desc::isConsistent());
}

SPxBasis::Desc& SPxBasis::Desc::operator=(const SPxBasis::Desc& rhs)
{
   if (this != &rhs)
   {
      rowstat = rhs.rowstat;
      colstat = rhs.colstat;
      
      if (rhs.stat == &rhs.rowstat)
      {
         assert(rhs.costat == &rhs.colstat);
         
         stat   = &rowstat;
         costat = &colstat;
      }
      else
      {
         assert(rhs.costat == &rhs.rowstat);
         
         stat   = &colstat;
         costat = &rowstat;
      }

      assert(Desc::isConsistent());
   }
   return *this;
}

void SPxBasis::Desc::reSize(int rowDim, int colDim)
{

   assert(rowDim >= 0);
   assert(colDim >= 0);

   int noldrows = rowstat.size();
   int noldcols = colstat.size();

   rowstat.reSize(rowDim);
   colstat.reSize(colDim);

   for( int i = rowDim - 1; i >= noldrows; i-- )
      rowstat[i] = D_UNDEFINED;

   for( int i = colDim - 1; i >= noldcols; i-- )
      colstat[i] = D_UNDEFINED;
}

void SPxBasis::Desc::dump() const
{
   int i;

   // Dump regardless of the verbosity level if this method is called.

   std::cout << "DBDESC01 column status: ";
   for(i = 0; i < nCols(); i++)
      std::cout << colStatus(i);
   std::cout << std::endl;

   std::cout << "DBDESC02 row status:    ";
   for(i = 0; i < nRows(); i++)
      std::cout << rowStatus(i);
   std::cout << std::endl;
}

bool SPxBasis::Desc::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   return rowstat.isConsistent() && colstat.isConsistent();
#else
   return true;
#endif
}

std::ostream& operator<<(std::ostream& os, const SPxBasis::Desc::Status& stat)
{
   char text;
   
   switch(stat)
   {
   case SPxBasis::Desc::P_ON_LOWER :
      text = 'L';
      break;
   case SPxBasis::Desc::P_ON_UPPER :
      text = 'U';
      break;
   case SPxBasis::Desc::P_FREE :
      text = 'F';
      break;
   case SPxBasis::Desc::P_FIXED :
      text = 'X';
      break;
   case SPxBasis::Desc::D_FREE :
      text = 'f';
      break;
   case SPxBasis::Desc::D_ON_UPPER :
      text = 'u';
      break;
   case SPxBasis::Desc::D_ON_LOWER :
      text = 'l';
      break;
   case SPxBasis::Desc::D_ON_BOTH :
      text = 'x';
      break;
   case SPxBasis::Desc::D_UNDEFINED :
      text = '.';
      break;
   default :
      os << std::endl << "Invalid status <" << int(stat) << ">" << std::endl;
      throw SPxInternalCodeException("XSPXDE01 This should never happen.");
   }
   os << text;

   return os;
}

} // namespace soplex
