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
#include "spxout.h"

namespace soplex
{

void SPxBasis::reDim()
{

   assert(theLP != 0);

   MSG_DEBUG( std::cout << "DCHBAS01 SPxBasis::reDim():"
                     << " matrixIsSetup=" << matrixIsSetup
                     << " fatorized=" << factorized
                     << std::endl; )

   thedesc.reSize (theLP->nRows(), theLP->nCols());

   if (theLP->dim() != matrix.size())
   {
      MSG_INFO3( (*spxout), (*spxout) << "ICHBAS02 basis redimensioning invalidates factorization"
                           << std::endl; )

      matrix.reSize (theLP->dim());
      theBaseId.reSize(theLP->dim());
      matrixIsSetup = false;
      factorized    = false;
   }

   MSG_DEBUG( std::cout << "DCHBAS03 SPxBasis::reDim(): -->"
                     << " matrixIsSetup=" << matrixIsSetup
                     << " fatorized=" << factorized
                     << std::endl; )

   assert( matrix.size()    >= theLP->dim() );
   assert( theBaseId.size() >= theLP->dim() );
}

/* adapt basis and basis descriptor to added rows */
void SPxBasis::addedRows(int n)
{
   assert(theLP != 0);

   if( n > 0 )
   {
      reDim();

      if (theLP->rep() == SPxSolver::COLUMN)
      {
         /* after adding rows in column representation, reDim() should set these bools to false. */
         assert( !matrixIsSetup && !factorized );

         for (int i = theLP->nRows() - n; i < theLP->nRows(); ++i)
         {
            thedesc.rowStatus(i) = dualRowStatus(i);
            baseId(i) = theLP->SPxLP::rId(i);
         }
      }
      else
      {
         assert(theLP->rep() == SPxSolver::ROW);

         for (int i = theLP->nRows() - n; i < theLP->nRows(); ++i)
            thedesc.rowStatus(i) = dualRowStatus(i);
      }

      /* If matrix was set up, load new basis vectors to the matrix.
       * In the row case, the basis is not effected by adding rows. However,
       * since @c matrix stores references to the rows in the LP (SPxLP), a realloc
       * in SPxLP (e.g. due to space requirements) might invalidate these references.
       * We therefore have to "reload" the matrix if it is set up. Note that reDim()
       * leaves @c matrixIsSetup untouched if only row have been added, since the basis
       * matrix already has the correct size. */
      if (status() > NO_PROBLEM && matrixIsSetup)
         loadMatrixVecs();

      /* update basis status */
      switch (status())
      {
      case PRIMAL:
      case UNBOUNDED:
         setStatus(REGULAR);
         break;
      case OPTIMAL:
      case INFEASIBLE:
         setStatus(DUAL);
         break;
      case NO_PROBLEM:
      case SINGULAR:
      case REGULAR:
      case DUAL:
         break;
      default:
         MSG_ERROR( std::cerr << "ECHBAS04 Unknown basis status!" << std::endl; )
         throw SPxInternalCodeException("XCHBAS01 This should never happen.");
      }
   }
}

void SPxBasis::removedRow(int i)
{

   assert(status() >  NO_PROBLEM);
   assert(theLP    != 0);

   if (theLP->rep() == SPxSolver::ROW)
   {
      if (theLP->isBasic(thedesc.rowStatus(i)))
      {
         setStatus(NO_PROBLEM);
         factorized = false;

         MSG_DEBUG( std::cout << "DCHBAS05 Warning: deleting basic row!\n"; )
      }
   }
   else
   {
      assert(theLP->rep() == SPxSolver::COLUMN);
      factorized = false;
      if (!theLP->isBasic(thedesc.rowStatus(i)))
      {
         setStatus(NO_PROBLEM);
         MSG_DEBUG( std::cout << "DCHBAS06 Warning: deleting nonbasic row!\n"; )
      }
      else if (status() > NO_PROBLEM && matrixIsSetup)
      {
         for (int j = theLP->dim(); j >= 0; --j)
         {
            SPxId id = baseId(j);

            if (id.isSPxRowId() && !theLP->has(SPxRowId(id)))
            {
               baseId(j) = baseId(theLP->dim());

               if (j < theLP->dim())
                  matrix[j] = &theLP->vector(baseId(j));
               break;
            }
         }
      }
   }
   thedesc.rowStatus(i) = thedesc.rowStatus(theLP->nRows());
   reDim();
}

void SPxBasis::removedRows(const int perm[])
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   int i;
   int n = thedesc.nRows();

   if (theLP->rep() == SPxSolver::ROW)
   {
      for (i = 0; i < n; ++i)
      {
         if (perm[i] != i)
         {
            if (perm[i] < 0)               // row got removed
            {
               if (theLP->isBasic(thedesc.rowStatus(i)))
               {
                  setStatus(NO_PROBLEM);
                  factorized = matrixIsSetup = false;
                  MSG_DEBUG( std::cout << "DCHBAS07 Warning: deleting basic row!\n"; )
               }
            }
            else                            // row was moved
               thedesc.rowStatus(perm[i]) = thedesc.rowStatus(i);
         }
      }
   }
   else
   {
      assert(theLP->rep() == SPxSolver::COLUMN);

      factorized    = false;
      matrixIsSetup = false;

      for (i = 0; i < n; ++i)
      {
         if (perm[i] != i)
         {
            if (perm[i] < 0)               // row got removed
            {
               if (!theLP->isBasic(thedesc.rowStatus(i)))
                  setStatus(NO_PROBLEM);
            }
            else                            // row was moved
               thedesc.rowStatus(perm[i]) = thedesc.rowStatus(i);
         }
      }
   }
   reDim();
}


static SPxBasis::Desc::Status
primalColStatus(int i, const SPxLP* theLP)
{
   assert(theLP != 0);

   if (theLP->upper(i) < infinity)
   {
      if (theLP->lower(i) > -infinity)
      {
         if (theLP->lower(i) == theLP->SPxLP::upper(i))
            return SPxBasis::Desc::P_FIXED;
         /*
             else
                 return (-theLP->lower(i) < theLP->upper(i))
                             ? SPxBasis::Desc::P_ON_LOWER
                          : SPxBasis::Desc::P_ON_UPPER;
         */
         else if (theLP->maxObj(i) == 0)
            return (-theLP->lower(i) < theLP->upper(i))
               ? SPxBasis::Desc::P_ON_LOWER
               : SPxBasis::Desc::P_ON_UPPER;
         else
            return (theLP->maxObj(i) < 0)
               ? SPxBasis::Desc::P_ON_LOWER
               : SPxBasis::Desc::P_ON_UPPER;
      }
      else
         return SPxBasis::Desc::P_ON_UPPER;
   }
   else if (theLP->lower(i) > -infinity)
      return SPxBasis::Desc::P_ON_LOWER;
   else
      return SPxBasis::Desc::P_FREE;
}


/* adapt basis and basis descriptor to added columns */
void SPxBasis::addedCols(int n)
{
   assert(theLP != 0);

   if( n > 0 )
   {
      reDim();

      if (theLP->rep() == SPxSolver::ROW)
      {
         /* after adding columns in row representation, reDim() should set these bools to false. */
         assert( !matrixIsSetup && !factorized );

         for (int i = theLP->nCols() - n; i < theLP->nCols(); ++i)
         {
            thedesc.colStatus(i) = primalColStatus(i, theLP);
            baseId(i) = theLP->SPxLP::cId(i);
         }
      }
      else
      {
         assert(theLP->rep() == SPxSolver::COLUMN);

         for (int i = theLP->nCols() - n; i < theLP->nCols(); ++i)
            thedesc.colStatus(i) = primalColStatus(i, theLP);
      }

      /* If matrix was set up, load new basis vectors to the matrix
       * In the column case, the basis is not effected by adding columns. However,
       * since @c matrix stores references to the columns in the LP (SPxLP), a realloc
       * in SPxLP (e.g. due to space requirements) might invalidate these references.
       * We therefore have to "reload" the matrix if it is set up. Note that reDim()
       * leaves @c matrixIsSetup untouched if only columns have been added, since the
       * basis matrix already has the correct size. */
      if (status() > NO_PROBLEM && matrixIsSetup)
         loadMatrixVecs();

      switch (status())
      {
      case DUAL:
      case INFEASIBLE:
         setStatus(REGULAR);
         break;
      case OPTIMAL:
      case UNBOUNDED:
         setStatus(PRIMAL);
         break;
      case NO_PROBLEM:
      case SINGULAR:
      case REGULAR:
      case PRIMAL:
         break;
      default:
         MSG_ERROR( std::cerr << "ECHBAS08 Unknown basis status!" << std::endl; )
         throw SPxInternalCodeException("XCHBAS02 This should never happen.");
      }
   }
}

void SPxBasis::removedCol(int i)
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   if (theLP->rep() == SPxSolver::COLUMN)
   {
      if (theLP->isBasic(thedesc.colStatus(i)))
         setStatus(NO_PROBLEM);
   }
   else
   {
      assert(theLP->rep() == SPxSolver::ROW);
      factorized = false;
      if (!theLP->isBasic(thedesc.colStatus(i)))
         setStatus(NO_PROBLEM);
      else if (status() > NO_PROBLEM)
      {
         for (int j = theLP->dim(); j >= 0; --j)
         {
            SPxId id = baseId(j);
            if (id.isSPxColId() && !theLP->has(SPxColId(id)))
            {
               baseId(j) = baseId(theLP->dim());
               if ( matrixIsSetup &&
                    j < theLP->dim() )
                  matrix[j] = &theLP->vector(baseId(j));
               break;
            }
         }
      }
   }

   thedesc.colStatus(i) = thedesc.colStatus(theLP->nCols());
   reDim();
}

void SPxBasis::removedCols(const int perm[])
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   int i;
   int n = thedesc.nCols();

   if (theLP->rep() == SPxSolver::COLUMN)
   {
      for (i = 0; i < n; ++i)
      {
         if (perm[i] < 0)           // column got removed
         {
            if (theLP->isBasic(thedesc.colStatus(i)))
               setStatus(NO_PROBLEM);
         }
         else                        // column was potentially moved
            thedesc.colStatus(perm[i]) = thedesc.colStatus(i);
      }
   }
   else
   {
      assert(theLP->rep() == SPxSolver::ROW);
      factorized = matrixIsSetup = false;
      for (i = 0; i < n; ++i)
      {
         if (perm[i] != i)
         {
            if (perm[i] < 0)               // column got removed
            {
               if (!theLP->isBasic(thedesc.colStatus(i)))
                  setStatus(NO_PROBLEM);
            }
            else                            // column was moved
               thedesc.colStatus(perm[i]) = thedesc.colStatus(i);
         }
      }
   }

   reDim();
}


/**
 * mark the basis as not factorized
 */
void SPxBasis::invalidate()
{
   if( factorized || matrixIsSetup )
   {
      MSG_INFO3( (*spxout), (*spxout) << "ICHBAS09 explicit invalidation of factorization" << std::endl; )
   }

   factorized    = false;
   matrixIsSetup = false;
}


/**
 * Create the initial slack basis descriptor and set up the basis matrix accordingly.
 * This code has been adapted from SPxBasis::addedRows() and SPxBasis::addedCols().
 */
void SPxBasis::restoreInitialBasis()
{
   assert(!factorized);

   MSG_INFO3( (*spxout), (*spxout) << "ICHBAS10 setup slack basis" << std::endl; )

   if (theLP->rep() == SPxSolver::COLUMN)
   {
      for (int i = 0; i < theLP->nRows(); ++i)
      {
         thedesc.rowStatus(i) = dualRowStatus(i);
         baseId(i) = theLP->SPxLP::rId(i);
      }

      for (int i = 0; i < theLP->nCols(); ++i)
         thedesc.colStatus(i) = primalColStatus(i, theLP);
   }
   else
   {
      assert(theLP->rep() == SPxSolver::ROW);

      for (int i = 0; i < theLP->nRows(); ++i)
         thedesc.rowStatus(i) = dualRowStatus(i);

      for (int i = 0; i < theLP->nCols(); ++i)
      {
         thedesc.colStatus(i) = primalColStatus(i, theLP);
         baseId(i) = theLP->SPxLP::cId(i);
      }
   }

   /* if matrix was set up, load new basis vectors to the matrix */
   if (status() > NO_PROBLEM && matrixIsSetup)
      loadMatrixVecs();

   /* update basis status */
   setStatus(REGULAR);
}

/**
 * @todo Implement changedRow(), changedCol(), changedElement() in a more clever
 * way. For instance, the basis won't be singular (but maybe infeasible) if the 
 * change doesn't affect the basis rows/columns.
 *
 * The following methods (changedRow(), changedCol(), changedElement()) radically
 * change the current basis to the original (slack) basis also present after 
 * loading the LP. The reason is that through the changes, the current basis may
 * become singular. Going back to the initial basis is quite inefficient, but 
 * correct.
 */

/**@todo is this correctly implemented?
 */
void SPxBasis::changedRow(int /*row*/)
{
   invalidate();
   restoreInitialBasis();
}

/**@todo is this correctly implemented?
 */
void SPxBasis::changedCol(int /*col*/)
{
   invalidate();
   restoreInitialBasis();
}

/**@todo is this correctly implemented?
 */
void SPxBasis::changedElement(int /*row*/, int /*col*/)
{
   invalidate();
   restoreInitialBasis();
}
} // namespace soplex
