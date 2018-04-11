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

#include <assert.h>

#include "spxdefines.h"
#include "spxweightpr.h"
#include "exceptions.h"

namespace soplex
{

void SPxWeightPR::setRep(SPxSolver::Representation rep)
{
   if (rep == SPxSolver::ROW)
   {
      penalty = rPenalty.get_const_ptr();
      coPenalty = cPenalty.get_const_ptr();
   }
   else
   {
      penalty = cPenalty.get_const_ptr();
      coPenalty = rPenalty.get_const_ptr();
   }
}

void SPxWeightPR::setType(SPxSolver::Type tp)
{
   if (thesolver && tp == SPxSolver::LEAVE)
   {
      leavePenalty.reDim( thesolver->dim() );
      computeLeavePenalty( 0, thesolver->dim() );
   }
}

void SPxWeightPR::computeLeavePenalty(int start, int end)
{
   const SPxBasis& basis = solver()->basis();

   for (int i = start; i < end; ++i)
   {
      SPxId id = basis.baseId(i);
      if (id.type() == SPxId::ROW_ID)
         leavePenalty[i] = rPenalty[ thesolver->number(id) ];
      else
         leavePenalty[i] = cPenalty[ thesolver->number(id) ];
   }
}

void SPxWeightPR::computeRP(int start, int end)
{
   for (int i = start; i < end; ++i)
   {
      /**@todo TK04NOV98 here is a bug.
       *       solver()->rowVector(i).length() could be zero, so
       *       solver()->rowVector(i).length2() is also zero and we
       *       get an arithmetic exception.
       */
      assert(solver()->rowVector(i).length() > 0);

      rPenalty[i] = (solver()->rowVector(i) * solver()->maxObj()) * objlength
                    / solver()->rowVector(i).length2();
      ASSERT_WARN( "WWGTPR01", rPenalty[i] > -1 - solver()->epsilon() );
   }
}

void SPxWeightPR::computeCP(int start, int end)
{
   for (int i = start; i < end; ++i)
   {
      cPenalty[i] = solver()->maxObj(i) * objlength;
      ASSERT_WARN( "WWGTPR02", cPenalty[i] > -1 - solver()->epsilon() );
   }
}

void SPxWeightPR::load(SPxSolver* base)
{
   thesolver = base;

   rPenalty.reDim(base->nRows());
   cPenalty.reDim(base->nCols());

   objlength = 1 / solver()->maxObj().length();
   computeCP(0, base->nCols());
   computeRP(0, base->nRows());
}

int SPxWeightPR::selectLeave()
{
   const Real* test = thesolver->fTest().get_const_ptr();
   Real type = 1 - 2 * (thesolver->rep() == SPxSolver::COLUMN ? 1 : 0);
   Real best = type * infinity;
   int lastIdx = -1;
   Real x;
   int i;

   for (i = solver()->dim() - 1; i >= 0; --i)
   {
      x = test[i];
      if (x < -theeps)
      {
         x *= leavePenalty[i];
         if (type * (x-best) < 0.0)
         {
            best = x;
            lastIdx = i;
         }
      }
   }
   assert(isConsistent());
   return lastIdx;
}

#if 0
/**@todo remove this code */
// ??? This is the old (buggy) version
int SPxWeightPR::selectLeave()
{
   const Real* test = thesolver->fTest().get_const_ptr();
   Real type = 1 - 2 * (thesolver->rep() == SPxSolver::COLUMN);
   Real best = type * infinity;
   int lastIdx = -1;
   Real x;
   int i;

   for (i = solver()->dim() - 1; i >= 0; --i)
   {
      x = test[i];
      if (x < -theeps)
      {
         x *= leavePenalty[i];
         if (x < best)
         {
            best = x;
            lastIdx = i;
         }
      }
   }
   assert(isConsistent());
   return lastIdx;
}
#endif

SPxId SPxWeightPR::selectEnter()
{
   const Vector& rTest = (solver()->rep() == SPxSolver::ROW)
                         ? solver()->test() : solver()->coTest();
   const Vector& cTest = (solver()->rep() == SPxSolver::ROW)
                         ? solver()->coTest() : solver()->test();
   const SPxBasis::Desc& ds = solver()->basis().desc();
   Real best = infinity;
   SPxId lastId;
   Real x;
   int i;

   for (i = solver()->nRows() - 1; i >= 0; --i)
   {
      x = rTest[i];
      if (x < -theeps)
      {
         x *= -x;
         switch (ds.rowStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER :
         case SPxBasis::Desc::D_ON_LOWER :
            x *= 1 + rPenalty[i];
            break;
         case SPxBasis::Desc::P_ON_UPPER :
         case SPxBasis::Desc::D_ON_UPPER :
            x *= 1 - rPenalty[i];
            break;
         case SPxBasis::Desc::P_FREE :
         case SPxBasis::Desc::D_FREE :
            return SPxId(solver()->rId(i));
         case SPxBasis::Desc::D_ON_BOTH :
            if (solver()->pVec()[i] > solver()->upBound()[i])
               x *= 1 + rPenalty[i];
            else
               x *= 1 - rPenalty[i];
            break;
         case SPxBasis::Desc::D_UNDEFINED :
         case SPxBasis::Desc::P_FIXED :
         default:
            throw SPxInternalCodeException("XWGTPR01 This should never happen.");
         }
         if (x < best)
         {
            best = x;
            lastId = solver()->rId(i);
         }
      }
   }

   for (i = solver()->nCols() - 1; i >= 0; --i)
   {
      x = cTest[i];
      if (x < -theeps)
      {
         x *= -x;
         switch (ds.colStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER :
         case SPxBasis::Desc::D_ON_LOWER :
            x *= 1 + cPenalty[i];
            break;
         case SPxBasis::Desc::P_ON_UPPER :
         case SPxBasis::Desc::D_ON_UPPER :
            x *= 1 - cPenalty[i];
            break;
         case SPxBasis::Desc::P_FREE :
         case SPxBasis::Desc::D_FREE :
            return SPxId(solver()->cId(i));
         case SPxBasis::Desc::D_ON_BOTH :
            if (solver()->coPvec()[i] > solver()->ucBound()[i])
               x *= 1 + cPenalty[i];
            else
               x *= 1 - cPenalty[i];
            break;
         case SPxBasis::Desc::P_FIXED :
         case SPxBasis::Desc::D_UNDEFINED :
         default:
            throw SPxInternalCodeException("XWGTPR02 This should never happen.");
         }
         if (x < best)
         {
            best = x;
            lastId = solver()->cId(i);
         }
      }
   }
   assert(isConsistent());
   return lastId;
}

void SPxWeightPR::addedVecs(int)
{
   if (solver()->rep() == SPxSolver::ROW)
   {
      int start = rPenalty.dim();
      rPenalty.reDim(solver()->nRows());
      computeRP(start, solver()->nRows());
   }
   else
   {
      int start = cPenalty.dim();
      cPenalty.reDim(solver()->nCols());
      computeCP(start, solver()->nCols());
   }
   if (solver()->type() == SPxSolver::LEAVE)
   {
      int start = leavePenalty.dim();
      leavePenalty.reDim( solver()->dim() );
      computeLeavePenalty( start, solver()->dim() );
   }
}

void SPxWeightPR::addedCoVecs(int)
{
   if (solver()->rep() == SPxSolver::COLUMN)
   {
      int start = rPenalty.dim();
      rPenalty.reDim(solver()->nRows());
      computeRP(start, solver()->nRows());
   }
   else
   {
      int start = cPenalty.dim();
      cPenalty.reDim(solver()->nCols());
      computeCP(start, solver()->nCols());
   }
   if (solver()->type() == SPxSolver::LEAVE)
   {
      int start = leavePenalty.dim();
      leavePenalty.reDim( solver()->dim() );
      computeLeavePenalty( start, solver()->dim() );
   }
}

void SPxWeightPR::removedVec(int i)
{
   assert(solver() != 0);

   if (solver()->rep() == SPxSolver::ROW)
   {
      rPenalty[i] = rPenalty[rPenalty.dim()];
      rPenalty.reDim(solver()->nRows());
   }
   else
   {
      cPenalty[i] = cPenalty[cPenalty.dim()];
      cPenalty.reDim(solver()->nCols());
   }
}

void SPxWeightPR::removedVecs(const int perm[])
{
   assert(solver() != 0);

   if (solver()->rep() == SPxSolver::ROW)
   {
      int j = rPenalty.dim();
      for (int i = 0; i < j; ++i)
      {
         if (perm[i] >= 0)
            rPenalty[perm[i]] = rPenalty[i];
      }
      rPenalty.reDim(solver()->nRows());
   }
   else
   {
      int j = cPenalty.dim();
      for (int i = 0; i < j; ++i)
      {
         if (perm[i] >= 0)
            cPenalty[perm[i]] = cPenalty[i];
      }
      cPenalty.reDim(solver()->nCols());
   }
}

void SPxWeightPR::removedCoVec(int i)
{
   assert(solver() != 0);

   if (solver()->rep() == SPxSolver::COLUMN)
   {
      rPenalty[i] = rPenalty[rPenalty.dim()];
      rPenalty.reDim(solver()->nRows());
   }
   else
   {
      cPenalty[i] = cPenalty[cPenalty.dim()];
      cPenalty.reDim(solver()->nCols());
   }
}

void SPxWeightPR::removedCoVecs(const int perm[])
{
   assert(solver() != 0);

   if (solver()->rep() == SPxSolver::COLUMN)
   {
      int j = rPenalty.dim();
      for (int i = 0; i < j; ++i)
      {
         if (perm[i] >= 0)
            rPenalty[perm[i]] = rPenalty[i];
      }
      rPenalty.reDim(solver()->nRows());
   }
   else
   {
      int j = cPenalty.dim();
      for (int i = 0; i < j; ++i)
      {
         if (perm[i] >= 0)
            cPenalty[perm[i]] = cPenalty[i];
      }
      cPenalty.reDim(solver()->nCols());
   }
}

bool SPxWeightPR::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (solver() != 0)
   {
      if (rPenalty.dim() != solver()->nRows())
         return MSGinconsistent("SPxWeightPR");
      if (cPenalty.dim() != solver()->nCols())
         return MSGinconsistent("SPxWeightPR");
   }
#endif

   return true;
}
} // namespace soplex
