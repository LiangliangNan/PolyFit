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
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "spxdefines.h"
#include "spxbasis.h"
#include "didxset.h"
#include "dvector.h"
#include "spxsolver.h"
#include "mpsinput.h"
#include "spxout.h"
#include "exceptions.h"

namespace soplex
{
SPxBasis::Desc::Status
SPxBasis::dualStatus(const SPxColId& id) const
{
   return dualColStatus(static_cast<SPxLP*>(theLP)->number(id));
}

SPxBasis::Desc::Status
SPxBasis::dualStatus(const SPxRowId& id) const
{
   return dualRowStatus((static_cast<SPxLP*>(theLP))->number(id));
}

SPxBasis::Desc::Status
SPxBasis::dualRowStatus(int i) const
{
   assert(theLP != 0);

   if (theLP->rhs(i) < infinity)
   {
      if (theLP->lhs(i) > -infinity)
      {
         if (theLP->lhs(i) == theLP->rhs(i))
            return Desc::D_FREE;
         else
            return Desc::D_ON_BOTH;
      }
      else
         return Desc::D_ON_LOWER;
   }
   else if (theLP->lhs(i) > -infinity)
      return Desc::D_ON_UPPER;
   else
      return Desc::D_UNDEFINED;
}

SPxBasis::Desc::Status
SPxBasis::dualColStatus(int i) const
{
   assert(theLP != 0);

   if (theLP->SPxLP::upper(i) < infinity)
   {
      if (theLP->SPxLP::lower(i) > -infinity)
      {
         if (theLP->SPxLP::lower(i) == theLP->SPxLP::upper(i))
            return Desc::D_FREE;
         else
            return Desc::D_ON_BOTH;
      }
      else
         return Desc::D_ON_LOWER;
   }
   else if (theLP->SPxLP::lower(i) > -infinity)
      return Desc::D_ON_UPPER;
   else
      return Desc::D_UNDEFINED;
}

void SPxBasis::loadMatrixVecs()
{
   assert(theLP != 0);
   assert(theLP->dim() == matrix.size());

   MSG_INFO3( (*spxout), (*spxout) << "IBASIS01 loadMatrixVecs() invalidates factorization" << std::endl; )

   int i;
   nzCount = 0;
   for (i = theLP->dim() - 1; i >= 0; --i)
   {
      matrix[i] = &theLP->vector(baseId(i));
      nzCount += matrix[i]->size();
   }
   matrixIsSetup = true;
   factorized = false;
   if (factor != 0)
      factor->clear();
}

bool SPxBasis::isDescValid(const Desc& ds)
{

   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   int basisdim;

   if ( ds.nRows() != theLP->nRows() || ds.nCols() != theLP->nCols() )
   {
      MSG_DEBUG( std::cout << "IBASIS20 Dimension mismatch\n" );
      return false;
   }

   basisdim = 0;
   for ( int row = ds.nRows()-1; row >= 0; --row )
   {
      if ( ds.rowstat[row] >= 0 )
      {
         if ( ds.rowstat[row] != dualRowStatus(row) )
         {
            MSG_DEBUG( std::cout << "IBASIS21 Basic row " << row << " with incorrect dual status " << dualRowStatus(row) << "\n");
            return false;
         }
      }
      else
      {
         basisdim++;
         if ( (ds.rowstat[row] == Desc::P_FIXED && theLP->SPxLP::lhs(row) != theLP->SPxLP::rhs(row))
              || (ds.rowstat[row] == Desc::P_ON_UPPER && theLP->SPxLP::rhs(row) >= infinity)
              || (ds.rowstat[row] == Desc::P_ON_LOWER && theLP->SPxLP::lhs(row) <= -infinity) )
         {
            MSG_DEBUG( std::cout << "IBASIS22 Nonbasic row with incorrect status: lhs=" << theLP->SPxLP::lhs(row) << ", rhs=" << theLP->SPxLP::rhs(row) << ", stat=" << ds.rowstat[row] << "\n");
            return false;
         }
      }
   }

   for ( int col = ds.nCols()-1; col >= 0; --col )
   {
      if ( ds.colstat[col] >= 0 )
      {
         if ( ds.colstat[col] !=  dualColStatus(col) )
         {
            MSG_DEBUG( std::cout << "IBASIS23 Basic column " << col << " with incorrect dual status " << ds.colstat[col] << " != " << dualColStatus(col) << "\n");
            return false;
         }
      }
      else
      {
         basisdim++;
         if ( (ds.colstat[col] == Desc::P_FIXED && theLP->SPxLP::lower(col) != theLP->SPxLP::upper(col))
              || (ds.colstat[col] == Desc::P_ON_UPPER && theLP->SPxLP::upper(col) >= infinity)
              || (ds.colstat[col] == Desc::P_ON_LOWER && theLP->SPxLP::lower(col) <= -infinity) )
         {
            MSG_DEBUG( std::cout << "IBASIS24 Nonbasic column with incorrect status: lower=" << theLP->SPxLP::lower(col) << ", upper=" << theLP->SPxLP::upper(col) << ", stat=" << ds.colstat[col] << "\n");
            return false;
         }
      }
   }

   if ( basisdim != theLP->nCols() )
   {
      MSG_DEBUG( std::cout << "IBASIS25 Incorrect basis dimension " << basisdim << " != " << theLP->nCols() << "\n" );
      return false;
   }

   // basis descriptor valid
   return true;
}

/*
    Loading a #Desc# into the basis can be done more efficiently, by
    explicitely programming both cases, for the rowwise and for the columnwise
    representation. This implementation hides this distinction in the use of
    methods #isBasic()# and #vector()#.
 */
void SPxBasis::loadDesc(const Desc& ds)
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);
   assert(ds.nRows() == theLP->nRows());
   assert(ds.nCols() == theLP->nCols());

   SPxId none;
   int   i;
   int   j;
   bool consistent = true;

   MSG_INFO3( (*spxout), (*spxout) << "IBASIS02 loading of Basis invalidates factorization" << std::endl; )

   lastin      = none;
   lastout     = none;
   lastidx     = -1;
   iterCount   = 0;
   updateCount = 0;

   if (&ds != &thedesc)
   {
      thedesc = ds;
      setRep();
   }

   assert(theLP->dim() == matrix.size());

   nzCount = 0;
   for (j = i = 0; i < theLP->nRows(); ++i)
   {
      /* for columns and rows with D_... status, the correct D_... status depends on bounds and sides; if a basis
       * descriptor is loaded after changing bounds or sides, e.g. in the refine() method, we have to correct them
       */
      if (thedesc.rowStatus(i) >= 0)
         thedesc.rowStatus(i) = dualRowStatus(i);
      else if (thedesc.rowStatus(i) == SPxBasis::Desc::P_FIXED && theLP->SPxLP::lhs(i) != theLP->SPxLP::rhs(i))
      {
         if (theLP->SPxLP::lhs(i) > -infinity && theLP->SPxLP::maxRowObj(i) < 0.0)
            thedesc.rowStatus(i) = SPxBasis::Desc::P_ON_LOWER;
         else if (theLP->SPxLP::rhs(i) < infinity)
            thedesc.rowStatus(i) = SPxBasis::Desc::P_ON_UPPER;
         else
            thedesc.rowStatus(i) = SPxBasis::Desc::P_FREE;
      }

      if (theLP->isBasic(thedesc.rowStatus(i)))
      {
         assert(theLP->dim() == matrix.size());
         assert(j <= matrix.size());

         if( j == matrix.size() )
         {
            // too many basic variables
            consistent = false;
            break;
         }

         SPxRowId id = theLP->SPxLP::rId(i);
         theBaseId[j] = id;
         matrix[j] = &theLP->vector(id);
         nzCount += matrix[j++]->size();
      }
   }

   for (i = 0; i < theLP->nCols(); ++i)
   {
      /* for columns and rows with D_... status, the correct D_... status depends on bounds and sides; if a basis
       * descriptor is loaded after changing bounds or sides, e.g. in the refine() method, we have to correct them
       */
      if (thedesc.colStatus(i) >= 0)
         thedesc.colStatus(i) = dualColStatus(i);
      else if (thedesc.colStatus(i) == SPxBasis::Desc::P_FIXED && theLP->SPxLP::lower(i) != theLP->SPxLP::upper(i))
      {
         if (theLP->SPxLP::lower(i) <= -infinity && theLP->SPxLP::upper(i) >= infinity)
            thedesc.colStatus(i) = SPxBasis::Desc::P_FREE;
         else if (theLP->SPxLP::upper(i) >= infinity || (theLP->SPxLP::lower(i) > -infinity && theLP->SPxLP::maxObj(i) < 0.0))
            thedesc.colStatus(i) = SPxBasis::Desc::P_ON_LOWER;
         else
            thedesc.colStatus(i) = SPxBasis::Desc::P_ON_UPPER;
      }

      if (theLP->isBasic(thedesc.colStatus(i)))
      {
         assert(theLP->dim() == matrix.size());
         assert(j <= matrix.size());

         if( j == matrix.size() )
         {
            // too many basic variables
            consistent = false;
            break;
         }

         SPxColId id = theLP->SPxLP::cId(i);
         theBaseId[j] = id;
         matrix[j] = &theLP->vector(id);
         nzCount += matrix[j++]->size();
      }
   }

   if( j < matrix.size() )
      consistent = false;  // not enough basic variables

   /* if dimensions are inconsistent, restore slack basis
    * if dimensions are consistent, then we have setup a correct matrix
    */
   if( !consistent )
      restoreInitialBasis();
   else
      matrixIsSetup = true;

   assert(isDescValid(thedesc));

   factorized = false;
   if (factor != 0)
      factor->clear();
}

void SPxBasis::setRep()
{
   assert(theLP != 0);

   reDim();
   minStab = 0.0;

   if (theLP->rep() == SPxSolver::ROW)
   {
      thedesc.stat   = &thedesc.rowstat;
      thedesc.costat = &thedesc.colstat;
   }
   else
   {
      thedesc.stat   = &thedesc.colstat;
      thedesc.costat = &thedesc.rowstat;
   }
}

void SPxBasis::load(SPxSolver* lp, bool initSlackBasis)
{
   assert(lp != 0);
   theLP = lp;

   setOutstream(*theLP->spxout);

   setRep();

   if( initSlackBasis )
   {
      restoreInitialBasis();
      loadDesc(thedesc);
   }
}

void SPxBasis::loadBasisSolver(SLinSolver* p_solver, const bool destroy)
{
   assert(!freeSlinSolver || factor != 0);

   setOutstream(*p_solver->spxout);

   MSG_INFO3( (*spxout), (*spxout) << "IBASIS03 loading of Solver invalidates factorization"
                        << std::endl; )

   if(freeSlinSolver)
   {
      delete factor;
      factor = 0;
   }

   factor = p_solver;
   factorized = false;
   factor->clear();
   freeSlinSolver = destroy;
}




/** 
 *  The specification is taken from the
 *
 *  ILOG CPLEX 7.0 Reference Manual, Appendix E, Page 543.
 *
 *  This routine should read valid BAS format files. 
 *
 *  @return true if the file was read correctly.
 *
 *  Here is a very brief outline of the format:
 *
 *  The format is in a form similar to an MPS file. The basic assumption is that all (column)
 *  variables are nonbasic at their lower bound and all row (variables) are basic; only the
 *  differences to this rule are given. Each data line contains an indicator, a variable name and
 *  possibly a row/constraint name. The following meaning applies with respect to the indicators:
 *
 *  - XU: the variable is basic, the row is nonbasic at its upper bound 
 *  - XL: the variable is basic, the row is nonbasic at its lower bound 
 *  - UL: the variable is nonbasic and at its upper bound
 *  - LL: the variable is nonbasic and at its lower bound
 *
 *  The CPLEX format contains an additional indicator 'BS', but this is unsupported here.
 *
 *  Nonbasic variables without lower bound have the following default status for SoPlex:
 *  - at their upper bound if finite,
 *  - at zero if free.
 */
bool SPxBasis::readBasis(
   std::istream&  is, 
   const NameSet* rowNames, 
   const NameSet* colNames)
{
   assert(theLP != 0);

   /* prepare names */
   const NameSet* rNames = rowNames;
   const NameSet* cNames = colNames;

   NameSet* p_colNames = 0;
   NameSet* p_rowNames = 0;

   if ( colNames == 0 )
   {
      int nCols = theLP->nCols();
      std::stringstream name;

      spx_alloc(p_colNames);
      p_colNames = new (p_colNames) NameSet();
      p_colNames->reMax(nCols);
      for (int j = 0; j < nCols; ++j)
      {
         name << "x" << j;
         DataKey key = theLP->colId(j);
         p_colNames->add(key, name.str().c_str());
      }
      cNames = p_colNames;
   }

   if ( rNames == 0 )
   {
      int nRows = theLP->nRows();
      std::stringstream name;

      spx_alloc(p_rowNames);
      p_rowNames = new (p_rowNames) NameSet();
      p_rowNames->reMax(nRows);
      for (int i = 0; i < nRows; ++i)
      {
         name << "C" << i;
         DataKey key = theLP->rowId(i);
         p_rowNames->add(key, name.str().c_str());
      }
      rNames = p_rowNames;
   }

   /* load default basis if necessary */
   if (status() == NO_PROBLEM)
      load(theLP, false);

   /* initialize with standard settings */
   Desc l_desc(thedesc);

   for (int i = 0; i < theLP->nRows(); i++)
      l_desc.rowstat[i] = dualRowStatus(i);

   for (int i = 0; i < theLP->nCols(); i++)
   {
      if (theLP->SPxLP::lower(i) == theLP->SPxLP::upper(i))
         l_desc.colstat[i] = Desc::P_FIXED;
      else if (theLP->SPxLP::lower(i) <= -infinity && theLP->SPxLP::upper(i) >= infinity)
         l_desc.colstat[i] = Desc::P_FREE;
      else if (theLP->SPxLP::lower(i) <= -infinity)
         l_desc.colstat[i] = Desc::P_ON_UPPER;
      else
         l_desc.colstat[i] = Desc::P_ON_LOWER;
   }

   MPSInput mps(is);

   if (mps.readLine() && (mps.field0() != 0) && !strcmp(mps.field0(), "NAME"))
   {
      while (mps.readLine())
      {
         int c = -1;
         int r = -1;

         if ((mps.field0() != 0) && !strcmp(mps.field0(), "ENDATA"))
         {
            mps.setSection(MPSInput::ENDATA);
            break;
         }
         if ((mps.field1() == 0) || (mps.field2() == 0))
            break;

         if ((c = cNames->number(mps.field2())) < 0)
            break;

         if (*mps.field1() == 'X')
            if (mps.field3() == 0 || (r = rNames->number(mps.field3())) < 0)
               break;

         if (!strcmp(mps.field1(), "XU"))
         {
            l_desc.colstat[c] = dualColStatus(c);
            if( theLP->LPRowSet::type(r) == LPRow::GREATER_EQUAL )
               l_desc.rowstat[r] = Desc::P_ON_LOWER;
            else if( theLP->LPRowSet::type(r) == LPRow::EQUAL )
               l_desc.rowstat[r] = Desc::P_FIXED;
            else
               l_desc.rowstat[r] = Desc::P_ON_UPPER;
         }
         else if (!strcmp(mps.field1(), "XL"))
         {
            l_desc.colstat[c] = dualColStatus(c);
            if( theLP->LPRowSet::type(r) == LPRow::LESS_EQUAL )
               l_desc.rowstat[r] = Desc::P_ON_UPPER;
            else if( theLP->LPRowSet::type(r) == LPRow::EQUAL )
               l_desc.rowstat[r] = Desc::P_FIXED;
            else
               l_desc.rowstat[r] = Desc::P_ON_LOWER;
         }
         else if (!strcmp(mps.field1(), "UL"))
         {
            l_desc.colstat[c] = Desc::P_ON_UPPER;
         }
         else if (!strcmp(mps.field1(), "LL"))
         {
            l_desc.colstat[c] = Desc::P_ON_LOWER;
         }
         else
         {
            mps.syntaxError();
            break;
         }
      }
   }
   if (!mps.hasError())
   {
      if (mps.section() == MPSInput::ENDATA)
      {
         // force basis to be different from NO_PROBLEM
         // otherwise the basis will be overwritten at later stages.
         setStatus(REGULAR);
         loadDesc(l_desc);
      }
      else
         mps.syntaxError();
   }

   if ( rowNames == 0 )
   {
      p_rowNames->~NameSet();
      spx_free(p_rowNames);
   }
   if ( colNames == 0 )
   {
      p_colNames->~NameSet();
      spx_free(p_colNames);
   }

#ifndef NDEBUG
   MSG_DEBUG( thedesc.dump() );
#endif

   return !mps.hasError();
}


/* Get row name - copied from spxmpswrite.cpp 
 *
 * @todo put this in a common file and unify accross different formats (mps, lp, basis).
 */
static const char* getRowName(
   const SPxLP*   lp,
   int            idx,
   const NameSet* rnames, 
   char*          buf)
{
   assert(buf != 0);
   assert(idx >= 0);
   assert(idx < lp->nRows());

   if (rnames != 0) 
   {
      DataKey key = lp->rId(idx);

      if (rnames->has(key))
         return (*rnames)[key];
   }
   spxSnprintf(buf, 16, "C%d", idx);
   
   return buf;
}
   
/* Get column name - copied from spxmpswrite.cpp 
 *
 * @todo put this in a common file and unify accross different formats (mps, lp, basis).
 */
static const char* getColName(
   const SPxLP*   lp,
   int            idx,
   const NameSet* cnames, 
   char*          buf)
{
   assert(buf != 0);
   assert(idx >= 0);
   assert(idx < lp->nCols());

   if (cnames != 0) 
   {
      DataKey key = lp->cId(idx);

      if (cnames->has(key))
         return (*cnames)[key];
   }
   spxSnprintf(buf, 16, "x%d", idx);
   
   return buf;
}

/* writes a file in MPS basis format to \p os.
 *
 * See SPxBasis::readBasis() for a short description of the format.
 */
void SPxBasis::writeBasis(
   std::ostream&  os, 
   const NameSet* rowNames, 
   const NameSet* colNames,
   const bool cpxFormat
   ) const
{
   assert(theLP != 0);

   os.setf(std::ios::left);
   os << "NAME  soplex.bas\n";

   /* do not write basis if there is none */
   if (status() == NO_PROBLEM)
   {
      os << "ENDATA" << std::endl;
      return;
   }

   /* start writing */
   char buf[255];
   int row = 0;
   for (int col = 0; col < theLP->nCols(); col++)
   {
      if ( thedesc.colStatus(col) > 0 )
      {
         /* Find non basic row */
         for (; row < theLP->nRows(); row++)
         {
            if ( thedesc.rowStatus(row) < 0 )
               break;
         }

         assert( row != theLP->nRows() );

         if( thedesc.rowStatus( row ) == Desc::P_ON_UPPER && (!cpxFormat || theLP->LPRowSet::type(row) == LPRow::RANGE) )
            os << " XU ";
         else
            os << " XL ";

         os << std::setw(8) << getColName(theLP, col, colNames, buf);

         /* break in two parts since buf is reused */
         os << "       " 
            << getRowName(theLP, row, rowNames, buf)
            << std::endl;

         row++;
      }
      else
      {
         if ( thedesc.colStatus( col ) == Desc::P_ON_UPPER )
         {
            os << " UL "
               << getColName(theLP, col, colNames, buf)
               << std::endl;
         }
         else
         {
            /* Default is all non-basic variables on lower bound (if finite) or at zero (if free).
             * nothing to do in this case.
             */
            assert(theLP->lower(col) <= -infinity || thedesc.colStatus(col) == Desc::P_ON_LOWER || thedesc.colStatus(col) == Desc::P_FIXED);
            assert(theLP->lower(col) > -infinity || theLP->upper(col) < infinity || thedesc.colStatus(col) == Desc::P_FREE);
         }
      }
   }

#ifndef NDEBUG
   MSG_DEBUG( thedesc.dump() );

   // Check that we covered all nonbasic rows - the remaining should be basic.
   for (; row < theLP->nRows(); row++)
   {
      if ( thedesc.rowStatus(row) < 0 )
         break;
   }
   assert( row == theLP->nRows() );

#endif // NDEBUG

   os << "ENDATA" << std::endl;
}

void SPxBasis::printMatrix() const
{

   assert(matrixIsSetup);

   for( int i = 0; i < matrix.size(); i++ )
   {
      std::cout << "C" << i << "=" << *matrix[i] << std::endl;
   }
}

void SPxBasis::printMatrixMTX(int number)
{
   int dim;
   int nnz;
   char filename[SPX_MAXSTRLEN];

   dim = matrix.size();
   nnz = nzCount;
   spxSnprintf(filename, SPX_MAXSTRLEN, "basis/basis%d.mtx", number);
   std::cout << "printing basis matrix to file " << filename << "\n";
   FILE * basisfile;
   basisfile = fopen (filename,"w");
   // print marker necessary for reading the file in Matlab
   fprintf(basisfile, "%%%%MatrixMarket matrix coordinate real general\n");
   // print matrix information
   fprintf(basisfile, "%d %d %d\n", dim, dim, nnz );
   // print matrix data
   for( int i = 0; i < matrix.size(); ++i )
   {
      for( int j = 0; j < baseVec(i).size(); ++j )
      {
         int idx = baseVec(i).index(j);
         Real val = baseVec(i).value(j);
         fprintf(basisfile, "%d %d %.13" REAL_FORMAT "\n",i+1,idx+1,val);
      }
   }
   fclose (basisfile);

   return;
}

void SPxBasis::change(
   int i,
   SPxId& id,
   const SVector* enterVec,
   const SSVector* eta)
{

   assert(matrixIsSetup);
   assert(!id.isValid() || (enterVec != 0));
   assert(factor != 0);

   lastidx = i;
   lastin  = id;

   if (id.isValid() && i >= 0)
   {
      assert(enterVec != 0);

      // update the counter for nonzeros in the basis matrix
      nzCount      = nzCount - matrix[i]->size() + enterVec->size();
      // let the new id enter the basis
      matrix[i]    = enterVec;
      lastout      = theBaseId[i];
      theBaseId[i] = id;

      ++iterCount;
      ++updateCount;

      MSG_DEBUG( std::cout << "factor_stats: iteration= " << iteration()
         << " update= " << updateCount
         << " total_update= " << totalUpdateCount
         << " nonzero_B= " << nzCount
         << " nonzero_LU= " << factor->memory()
         << " factor_fill= " << lastFill
         << " time= " << theLP->time()
         << std::endl; )

      // never factorize? Just do it !
      if (!factorized)
         factorize();

      // too much memory growth ?
      else if (Real(factor->memory()) > 1000 + factor->dim() + lastMem * memFactor)
      {
         MSG_INFO3( (*spxout), (*spxout) << "IBASIS04 memory growth factor triggers refactorization"
                              << " memory= " << factor->memory()
                              << " lastMem= " << lastMem
                              << " memFactor= " << memFactor
                              << std::endl; )
         factorize();
      }

      // relative fill too high ?
      else if (Real(factor->memory()) > lastFill * Real(nzCount))
      {
         MSG_INFO3( (*spxout), (*spxout) << "IBASIS04 fill factor triggers refactorization"
                              << " memory= " << factor->memory()
                              << " nzCount= " << nzCount
                              << " lastFill= " << lastFill
                              << std::endl; )

         factorize();
      }
      // absolute fill in basis matrix too high ?
      else if (nzCount > lastNzCount)
      {
         MSG_INFO3( (*spxout), (*spxout) << "IBASIS05 nonzero factor triggers refactorization"
                              << " nzCount= " << nzCount
                              << " lastNzCount= " << lastNzCount
                              << " nonzeroFactor= " << nonzeroFactor
                              << std::endl; )
         factorize();
      }
      // too many updates ?
      else if (updateCount >= maxUpdates)
      {
         MSG_INFO3( (*spxout), (*spxout) << "IBASIS06 update count triggers refactorization"
                              << " updateCount= " << updateCount
                              << " maxUpdates= " << maxUpdates
                              << std::endl; )
         factorize();
      }
      else
      {
         try
         {
#ifdef MEASUREUPDATETIME
            theTime.start();
#endif
            factor->change(i, *enterVec, eta);
            totalUpdateCount++;
#ifdef MEASUREUPDATETIME
            theTime.stop();
#endif
         }
         catch( ... )
         {
            MSG_INFO3( (*spxout), (*spxout) << "IBASIS13 problems updating factorization; refactorizing basis"
               << std::endl; )

#ifdef MEASUREUPDATETIME
            theTime.stop();
#endif

            // singularity was detected in update; we refactorize
            factorize();

            // if factorize() detects singularity, an exception is thrown, hence at this point we have a regular basis
            // and can try the update again
            assert(status() >= SPxBasis::REGULAR);
            try
            {
#ifdef MEASUREUPDATETIME
               theTime.start();
#endif
               factor->change(i, *enterVec, eta);
               totalUpdateCount++;
#ifdef MEASUREUPDATETIME
               theTime.stop();
#endif
            }
            // with a freshly factorized, regular basis, the update is unlikely to fail; if this happens nevertheless,
            // we have to invalidate the basis to have the statuses correct
            catch( const SPxException& F )
            {
               MSG_INFO3( (*spxout), (*spxout) << "IBASIS14 problems updating factorization; invalidating factorization"
                  << std::endl; )

#ifdef MEASUREUPDATETIME
               theTime.stop();
#endif

               factorized = false;
               throw F;
            }
         }

         assert(minStab > 0.0);

         if (factor->status() != SLinSolver::OK || factor->stability() < minStab)
         {
            MSG_INFO3( (*spxout), (*spxout) << "IBASIS07 stability triggers refactorization"
                                 << " stability= " << factor->stability()
                                 << " minStab= " << minStab
                                 << std::endl; )
            factorize();
         }
      }
   }
   else
      lastout = id;
}

void SPxBasis::factorize()
{

   assert(factor != 0);

   if (!matrixIsSetup)
      loadDesc(thedesc);

   assert(matrixIsSetup);

   updateCount = 0;

   switch(factor->load(matrix.get_ptr(), matrix.size()))
   {
   case SLinSolver::OK :
      if (status() == SINGULAR)
         setStatus(REGULAR);

      factorized = true;
      minStab = factor->stability();

      // This seems always to be about 1e-7
      if (minStab > 1e-4)
         minStab *= 0.001;
      if (minStab > 1e-5)
         minStab *= 0.01;
      if (minStab > 1e-6)
         minStab *= 0.1;
      break;
   case SLinSolver::SINGULAR :
      setStatus(SINGULAR);
      factorized = false;
      break;
   default :
      MSG_ERROR( std::cerr << "EBASIS08 error: unknown status of factorization.\n"; )
      factorized = false;
      throw SPxInternalCodeException("XBASIS01 This should never happen.");
   }

   // get nonzero count of factorization
   lastMem    = factor->memory();
   // compute fill ratio between factorization and basis matrix (multiplied with tolerance)
   lastFill   = fillFactor * Real(lastMem) / Real(nzCount > 0 ? nzCount : 1);
   lastNzCount = int(nonzeroFactor * Real(nzCount > 0 ? nzCount : 1));

   if (status() == SINGULAR)
   {
      throw SPxStatusException("Cannot factorize singular matrix");
   }
}

Vector& SPxBasis::multWithBase(Vector& x) const
{
   assert(status() > SINGULAR);
   assert(theLP->dim() == x.dim());

   int i;
   DVector tmp(x);

   if (!matrixIsSetup)
      (const_cast<SPxBasis*>(this))->loadDesc(thedesc);

   assert( matrixIsSetup );

   for (i = x.dim() - 1; i >= 0; --i)
      x[i] = *(matrix[i]) * tmp;

   return x;
}

void SPxBasis::multWithBase(SSVector& x, SSVector& result) const
{
   assert(status() > SINGULAR);
   assert(theLP->dim() == x.dim());
   assert(x.dim() == result.dim());

   if (!matrixIsSetup)
      (const_cast<SPxBasis*>(this))->loadDesc(thedesc);

   result.clear();

   assert(matrixIsSetup);

   for( int i = 0; i < x.dim(); ++i )
      result.add(i, (*matrix[i]) * x);

   return;
}


Vector& SPxBasis::multBaseWith(Vector& x) const
{
   assert(status() > SINGULAR);
   assert(theLP->dim() == x.dim());

   int i;
   DVector tmp(x);

   if (!matrixIsSetup)
      (const_cast<SPxBasis*>(this))->loadDesc(thedesc);

   assert( matrixIsSetup );

   x.clear();
   for (i = x.dim() - 1; i >= 0; --i)
   {
      if (tmp[i] != 0.0)
         x.multAdd(tmp[i], *(matrix[i]));
   }

   return x;
}

void SPxBasis::multBaseWith(SSVector& x, SSVector& result) const
{
   assert(status() > SINGULAR);
   assert(theLP->dim() == x.dim());
   assert(x.dim() == result.dim());

   if (!matrixIsSetup)
      (const_cast<SPxBasis*>(this))->loadDesc(thedesc);

   assert(matrixIsSetup);

   result.clear();
   if( x.isSetup() )
   {
      for( int i = 0; i < x.size(); ++i )
      {
         int idx = x.index(i);
         result.multAdd(x[idx], (*matrix[idx]));
      }
   }
   else
   {
      for( int i = 0; i < x.dim(); ++i )
         result.multAdd(x[i], (*matrix[i]));
   }

   return;
}

/* compute an estimated condition number for the current basis matrix
 * by computing estimates of the norms of B and B^-1 using the power method.
 * maxiters and tolerance control the accuracy of the estimate.
 */
Real SPxBasis::condition(int maxiters, Real tolerance)
{
   int dimension = matrix.size();
   int miniters = 3;    // minimum number of power method iterations
   int i;
   int c;
   Real norm;
   Real norminv;
   Real norm1;
   Real norm2;

   // catch corner case of empty matrix
   if( dimension == 0 )
      return 1.0;

   SSVector x(dimension);
   SSVector y(dimension);

   // check whether a regular basis matrix is available
   if( status() < REGULAR )
      return 0;

   if( !matrixIsSetup )
      (const_cast<SPxBasis*>(this))->loadDesc(thedesc);

   if( !factorized )
      factorize();

   // initialize vectors
   norm1 = 1.0 / (Real) dimension;
   for( i = 0; i < dimension; i++ )
      x.add(i, norm1);
   y = x;

   // compute norm of B
   for( c = 0; c < maxiters; ++c )
   {
      norm2 = norm1;

      // y = B*x
      multBaseWith(x, y);
      norm1 = y.length();

      // stop if converged
      if( c >= miniters && spxAbs(norm1 - norm2) < tolerance * norm1 )
         break;

      // x = B^T*y and normalize
      multWithBase(y, x);
      norm2 = 1.0 / x.length();
      x *= norm2;
   }
   norm = norm1;

   // reinitialize vectors
   x.clear();
   y.clear();
   norm1 = 1.0 / (Real) dimension;
   for( i = 0; i < dimension; i++ )
      x.add(i, norm1);
   y = x;

   // compute norm of B^-1
   for( c = 0; c < maxiters; ++c )
   {
      norm2 = norm1;

      // x = B^-1*y
      factor->solveRight(x, y);
      x.setup();
      norm1 = x.length();

      // stop if converged
      if( c >= miniters && spxAbs(norm1 - norm2) < tolerance * norm1 )
         break;

      // y = B^-T*x and normalize
      factor->solveLeft(y, x);
      y.setup();
      norm2 = 1.0 / y.length();
      y *= norm2;
   }
   norminv = norm1;

   return norm * norminv;
}

/* compute condition number estimation based on the diagonal of the LU factorization */
Real SPxBasis::getFastCondition(int type)
{
   Real cond = infinity;

   if( factorized )
      cond = factor->conditionEstimate(type);

   return cond;
}

void SPxBasis::dump()
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);
   assert(thedesc.nRows() == theLP->nRows());
   assert(thedesc.nCols() == theLP->nCols());
   assert(theLP->dim() == matrix.size());

   int i, basesize;

   // Dump regardless of the verbosity level if this method is called.

   std::cout << "DBASIS09 Basis entries:";
   basesize = 0;
   for (i = 0; i < theLP->nRows(); ++i)
   {
      if (theLP->isBasic(thedesc.rowStatus(i)))
      {
         if(basesize % 10 == 0)
            std::cout << std::endl << "DBASIS10 ";
         SPxRowId id = theLP->SPxLP::rId(i);
         std::cout << "\tR" << theLP->number(id);
         basesize++;
      }
   }

   for (i = 0; i < theLP->nCols(); ++i)
   {
      if (theLP->isBasic(thedesc.colStatus(i)))
      {
         if(basesize % 10 == 0)
            std::cout << std::endl << "DBASIS11 ";
         SPxColId id = theLP->SPxLP::cId(i);
         std::cout << "\tC" << theLP->number(id);
         basesize++;
      }
   }
   std::cout << std::endl;

   assert(basesize == matrix.size());
}


bool SPxBasis::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   int primals = 0;
   int i;

   if (status() > NO_PROBLEM)
   {
      if (theLP == 0)
         return MSGinconsistent("SPxBasis");

      if (theBaseId.size() != theLP->dim() || matrix.size() != theLP->dim())
         return MSGinconsistent("SPxBasis");

      if (thedesc.nCols() != theLP->nCols() || thedesc.nRows() != theLP->nRows())
         return MSGinconsistent("SPxBasis");

      for (i = 0; i < thedesc.nRows(); ++i)
      {
         if (thedesc.rowStatus(i) >= 0)
         {
            if (thedesc.rowStatus(i) != dualRowStatus(i))
               return MSGinconsistent("SPxBasis");
         }
         else
            ++primals;
      }
      
      for (i = 0; i < thedesc.nCols(); ++i)
      {
         if (thedesc.colStatus(i) >= 0)
         {
            if (thedesc.colStatus(i) != dualColStatus(i))
               return MSGinconsistent("SPxBasis");
         }
         else
            ++primals;
      }
      if (primals != thedesc.nCols())
         return MSGinconsistent("SPxBasis");
   }
   return thedesc.isConsistent() && theBaseId.isConsistent() 
      && matrix.isConsistent() && factor->isConsistent();
#else
   return true;
#endif // CONSISTENCY_CHECKS
}

SPxBasis::SPxBasis(Timer::TYPE ttype)
   : theLP (0)
   , matrixIsSetup (false)
   , factor (0)
   , factorized (false)
   , maxUpdates (200)
   , nonzeroFactor(10.0)
   , fillFactor(5.0)
   , memFactor(1.5)
   , iterCount (0)
   , lastIterCount(0)
   , iterDegenCheck(0)
   , updateCount(0)
   , totalUpdateCount(0)
   , nzCount (1)
   , lastMem(0)
   , lastFill(0)
   , lastNzCount(0)
   , theTime(0)
   , timerType(ttype)
   , lastidx(0)
   , minStab(0.0)
   , thestatus (NO_PROBLEM)
   , freeSlinSolver(false)
   , spxout(0)
{
   // info: is not consistent at this moment, e.g. because theLP == 0

   theTime = TimerFactory::createTimer(timerType);
}

/**@warning Do not change the LP object.
 *  Only pointer to that object is copied.
 *  Hint: no problem, we use this function for copy
 *   constructor of SPxSolver
 */
SPxBasis::SPxBasis(const SPxBasis& old)
   : theLP(old.theLP)
   , theBaseId(old.theBaseId)
   , matrix(old.matrix)
   , matrixIsSetup(old.matrixIsSetup)
   , factor(old.factor)
   , factorized(old.factorized)
   , maxUpdates(old.maxUpdates)
   , nonzeroFactor(old.nonzeroFactor)
   , fillFactor(old.fillFactor)
   , memFactor(old.memFactor)
   , iterCount(old.iterCount)
   , lastIterCount(old.lastIterCount)
   , iterDegenCheck(old.iterDegenCheck)
   , updateCount(old.updateCount)
   , totalUpdateCount(old.totalUpdateCount)
   , nzCount(old.nzCount)
   , lastMem(old.lastMem)
   , lastFill(old.lastFill)
   , lastNzCount(old.lastNzCount)
   , theTime(old.theTime)
   , timerType(old.timerType)
   , lastin(old.lastin)
   , lastout(old.lastout)
   , lastidx(old.lastidx)
   , minStab(old.minStab)
   , thestatus(old.thestatus)
   , thedesc(old.thedesc)
   , spxout(old.spxout)
{

   this->factor = old.factor->clone();
   freeSlinSolver = true;

   assert(SPxBasis::isConsistent());
}

SPxBasis::~SPxBasis()
{

   assert(!freeSlinSolver || factor != 0);

   if(freeSlinSolver)
   {
      delete factor;
      factor = 0;
   }

   spx_free(theTime);
}


/**@warning  Note that we do not create a deep copy of the corresponding SPxSolver object.
 *  Only the reference to that object is copied.
 */
SPxBasis& SPxBasis::operator=(const SPxBasis& rhs)
{

   assert(!freeSlinSolver || factor != 0);

   if (this != &rhs)
   {
      theLP         = rhs.theLP;
      theBaseId     = rhs.theBaseId;
      matrix        = rhs.matrix;
      matrixIsSetup = rhs.matrixIsSetup;

      if(freeSlinSolver)
      {
         delete factor;
         factor = 0;
      }
      factor = rhs.factor->clone();
      freeSlinSolver = true;
      
      factorized    = rhs.factorized;
      maxUpdates    = rhs.maxUpdates;
      nonzeroFactor = rhs.nonzeroFactor;
      fillFactor    = rhs.fillFactor;
      memFactor     = rhs.memFactor;
      iterCount     = rhs.iterCount;
      nzCount       = rhs.nzCount;
      lastFill      = rhs.lastFill;
      lastNzCount   = rhs.lastNzCount;
      lastin        = rhs.lastin;
      lastout       = rhs.lastout;
      lastidx       = rhs.lastidx;
      minStab       = rhs.minStab;
      thestatus     = rhs.thestatus;
      thedesc       = rhs.thedesc;

      assert(SPxBasis::isConsistent());
   }

   return *this;
}


//
// Auxiliary functions.
//

// Pretty-printing of basis status.
std::ostream& operator<<( std::ostream& os,
                          const SPxBasis::SPxStatus& status )
{
   switch ( status )
      {
      case SPxBasis::NO_PROBLEM:
         os << "NO_PROBLEM";
         break;
      case SPxBasis::SINGULAR:
         os << "SINGULAR";
         break;
      case SPxBasis::REGULAR:
         os << "REGULAR";
         break;
      case SPxBasis::DUAL:
         os << "DUAL";
         break;
      case SPxBasis::PRIMAL:
         os << "PRIMAL";
         break;
      case SPxBasis::OPTIMAL:
         os << "OPTIMAL";
         break;
      case SPxBasis::UNBOUNDED:
         os << "UNBOUNDED";
         break;
      case SPxBasis::INFEASIBLE:
         os << "INFEASIBLE";
         break;
      default:
         os << "?other?";
         break;
      }
   return os;
}


} // namespace soplex
