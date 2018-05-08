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

#ifndef SOPLEX_LEGACY
#include <iostream>
#include <assert.h>

#include "soplex.h"
#include "statistics.h"

namespace soplex
{
   /// check scaling of LP
   void SoPlex::_checkScaling(SPxLPReal* origLP) const
   {
      MSG_INFO1( spxout, spxout << "DEBUG: checking correctness of scaled LP" << std::endl; )
      assert(_realLP->nCols() == origLP->nCols());
      assert(_realLP->nRows() == origLP->nRows());
      assert(_realLP->isScaled() && !origLP->isScaled());
      bool correct = true;

      MSG_INFO1( spxout, spxout << "DEBUG: checking rows..." << std::endl; )
      for( int i = 0; i < origLP->nRows(); ++i )
      {
         assert(EQ(origLP->lhs(i), _realLP->lhsUnscaled(i)));
         assert(EQ(origLP->rhs(i), _realLP->rhsUnscaled(i)));

         DSVectorReal row;
         _realLP->getRowVectorUnscaled(i, row);

         assert(origLP->rowVector(i).size() == row.size());

         for( int j = 0; j < row.size(); ++j)
         {
            if( NE(row.value(j), origLP->rowVector(i).value(j)) )
            {
               MSG_INFO1( spxout, spxout << "DEBUG: scaling error in row " << i << ", col " << j
                          << ": orig " << origLP->rowVector(i).value(j)
                          << ", unscaled: " << row.value(j) << std::endl; )
               correct = false;
            }
         }
      }

      MSG_INFO1( spxout, spxout << "DEBUG: checking cols..." << std::endl; )
      for( int i = 0; i < origLP->nCols(); ++i )
      {
         assert(EQ(origLP->lower(i), _realLP->lowerUnscaled(i)));
         assert(EQ(origLP->upper(i), _realLP->upperUnscaled(i)));
         assert(EQ(origLP->obj(i), _realLP->objUnscaled(i)));

         DSVectorReal col;
         _realLP->getColVectorUnscaled(i, col);

         assert(origLP->colVector(i).size() == col.size());

         for( int j = 0; j < col.size(); ++j)
         {
            if( NE(col.value(j), origLP->colVector(i).value(j), _solver.feastol()) )
            {
               MSG_INFO1( spxout, spxout << "DEBUG: scaling error in col " << i << ", row " << j
                          << ": orig " << origLP->colVector(i).value(j)
                          << ", unscaled: " << col.value(j) << std::endl; )
               correct = false;
            }
         }
      }
      if( !correct )
      {
         MSG_INFO1( spxout, spxout << "DEBUG: scaling check failed" << std::endl; )
      }
      assert(correct);
   }



   void SoPlex::_checkBasisScaling()
   {
      if( _status != SPxSolver::OPTIMAL )
      {
         MSG_INFO1( spxout, spxout << "DEBUG: skipping test on non optimal bases\n" );
         return;
      }

      assert(&_solver == _realLP);
      DVector** binvcol = 0;
      DVector** binvrow = 0;
      int* inds = 0;
      int basisdim = _solver.nRows(); // do all operations with regard to the column basis
      bool colrep = (_solver.rep() == SPxSolver::COLUMN);
      spx_alloc(binvcol, basisdim);
      spx_alloc(binvrow, basisdim);
      spx_alloc(inds, basisdim);

      if( colrep )
      {
         MSG_INFO1( spxout, spxout << "DEBUG: computing columns of inverted basis matrix\n";)
         // collect columns of the basis inverse
         for( int i = 0; i < basisdim; ++i)
         {
            binvcol[i] = new DVector(basisdim);
            binvcol[i]->clear();
            assert(getBasisInverseColReal(i, binvcol[i]->get_ptr(), 0, 0, true));
         }
      }
      MSG_INFO1( spxout, spxout << "DEBUG: computing rows of inverted basis matrix\n";)
      // collect rows of the basis inverse
      for( int i = 0; i < basisdim; ++i)
      {
         binvrow[i] = new DVector(basisdim);
         binvrow[i]->clear();
         assert(getBasisInverseRowReal(i, binvrow[i]->get_ptr(), 0, 0, true));
      }

      if( colrep )
      {
         MSG_INFO1( spxout, spxout << "DEBUG: checking columns for identity after multiplying with basis matrix\n";)
         // multiply with (unscaled) basis matrix and check result (should be unitvecs)
         for( int i = 0; i < basisdim; ++i)
         {
            DVector result(*binvcol[i]);
            assert(multBasis(result.get_ptr(), true));
            Real sumerror = 0.0;
            for( int j = 0; j < basisdim; ++j)
            {
               Real error = 0.0;
               if( j != i )
                  error = spxAbs(result[j]);
               else
                  error = spxAbs(result[j] - 1.0);
               if( error > _solver.feastol() )
                  MSG_INFO1( spxout, spxout << "ERROR: col " << i << " " << j << ", " << result[j] << std::endl );
               sumerror += error;
            }
            assert(_solver.rep() == SPxSolver::ROW || sumerror < _solver.feastol());
         }
      }

      MSG_INFO1( spxout, spxout << "DEBUG: checking rows for identity after multiplying with basis matrix\n";)
      for( int i = 0; i < basisdim; ++i)
      {
         DVector result(*binvrow[i]);
         assert(multBasisTranspose(result.get_ptr(), true));
         Real sumerror = 0.0;
         for( int j = 0; j < basisdim; ++j)
         {
            Real error = 0.0;
            if( j != i )
               error = spxAbs(result[j]);
            else
               error = spxAbs(result[j] - 1.0);
            if( error > _solver.feastol() )
               MSG_INFO1( spxout, spxout << "ERROR: row " << i << " " << j << ", " << result[j] << std::endl );
            sumerror += error;
         }
         assert(sumerror < _solver.feastol());
      }

      if( _solver.isScaled() )
      {
         MSG_INFO1( spxout, spxout << "DEBUG: unscaling LP\n"; )
//         _solver.setRep(SPxSolver::COLUMN);
         _solver.unscaleLPandReloadBasis();
//         _solver.setBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
//         _solver.solve();

         DVector** binvcol2 = 0;
         DVector** binvrow2 = 0;
         spx_alloc(binvcol2, basisdim);
         spx_alloc(binvrow2, basisdim);

         if( colrep )
         {
            MSG_INFO1( spxout, spxout << "DEBUG: computing columns of inverted basis matrix again\n";)
            // collect columns of the basis inverse
            for( int i = 0; i < basisdim; ++i)
            {
               binvcol2[i] = new DVector(basisdim);
               binvcol2[i]->clear();
               assert(getBasisInverseColReal(i, binvcol2[i]->get_ptr(), 0, 0, false));
            }
         }

         MSG_INFO1( spxout, spxout << "DEBUG: computing rows of inverted basis matrix again\n";)
         // collect rows of the basis inverse
         for( int i = 0; i < basisdim; ++i)
         {
            binvrow2[i] = new DVector(basisdim);
            binvrow2[i]->clear();
            assert(getBasisInverseRowReal(i, binvrow2[i]->get_ptr(), 0, 0, false));
         }

         MSG_INFO1( spxout, spxout << "DEBUG: checking rows and columns of scaled/unscaled inverted of basis matrix\n";)
         for( int i = 0; i < basisdim; ++i)
         {
            Real sumerror = 0.0;
            for( int j = 0; j < basisdim; ++j)
            {
               if( colrep )
               {
                  if( NE((*binvcol[i])[j], (*binvcol2[i])[j], _solver.feastol()) )
                  {
                     MSG_INFO1( spxout, spxout << "ERROR: col " << i << " " << j << ", " << (*binvcol[i])[j] << " " << (*binvcol2[i])[j] << std::endl );
                     sumerror += spxAbs((*binvcol[i])[j] - (*binvcol2[i])[j]);
                  }
               }
               if( NE((*binvrow[i])[j], (*binvrow2[i])[j], _solver.feastol()) )
               {
                  MSG_INFO1( spxout, spxout << "ERROR: row " << i << " " << j << ", " << (*binvrow[i])[j] /  (*binvrow2[i])[j] << std::endl );
                  sumerror += spxAbs((*binvrow[i])[j] - (*binvrow2[i])[j]);
               }
            }
            assert(sumerror < _solver.feastol());
         }

         spx_free(binvcol2);
         spx_free(binvrow2);
      }

      spx_free(inds);
      spx_free(binvrow);
      spx_free(binvcol);
   }

} // namespace soplex
#endif
