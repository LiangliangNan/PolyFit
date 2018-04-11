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

#define ALLOWED_UNSCALE_PERCENTAGE    0.1
#define MIN_OPT_CALLS_WITH_SCALING     10

namespace soplex
{
   /// solves real LP
   void SoPlex::_optimizeReal()
   {
      assert(_realLP != 0);
      assert(_realLP == &_solver);

      _solReal.invalidate();
      ++_optimizeCalls;

      // start timing
      _statistics->solvingTime->start();

      if( boolParam(SoPlex::PERSISTENTSCALING) )
      {
         // scale original problem; overwriting _realLP
         if( _scaler && !_realLP->isScaled() && _reapplyPersistentScaling() )
         {
#ifdef SOPLEX_DEBUG
            SPxLPReal* origLP = 0;
            spx_alloc(origLP);
            origLP = new (origLP) SPxLPReal(*_realLP);
#endif
            _scaler->scale(*_realLP, true);
            _isRealLPScaled = _realLP->isScaled(); // a scaler might decide not to apply scaling
#ifdef SOPLEX_DEBUG
            _checkScalingReal(origLP);
#endif
         }
         // unscale previously scaled problem, overwriting _realLP
         else if( !_scaler && _realLP->isScaled() )
         {
            _solver.unscaleLPandReloadBasis();
            _isRealLPScaled = false;
            ++_unscaleCalls;
         }
      }

      // remember that last solve was in floating-point
      _lastSolveMode = SOLVEMODE_REAL;

      // solve and store solution; if we have a starting basis, do not apply preprocessing; if we are solving from
      // scratch, apply preprocessing according to parameter settings
      if( !_hasBasis && realParam(SoPlex::OBJLIMIT_LOWER) == -realParam(SoPlex::INFTY) && realParam(SoPlex::OBJLIMIT_UPPER) == realParam(SoPlex::INFTY) )
         _preprocessAndSolveReal(true);
      else
         _preprocessAndSolveReal(false);

      _statistics->finalBasisCondition = _solver.getFastCondition();

      // stop timing
      _statistics->solvingTime->stop();
   }



   /// check whether persistent scaling is supposed to be reapplied again after unscaling
   bool SoPlex::_reapplyPersistentScaling() const
   {
      if( (_unscaleCalls > _optimizeCalls * ALLOWED_UNSCALE_PERCENTAGE) && _optimizeCalls > MIN_OPT_CALLS_WITH_SCALING )
         return false;
      else
         return true;
   }



   /// checks result of the solving process and solves again without preprocessing if necessary
   void SoPlex::_evaluateSolutionReal(SPxSimplifier::Result simplificationStatus)
   {
      // if the simplifier detected infeasibility or unboundedness we optimize again
      // just to get the proof (primal or dual ray)
      // todo get infeasibility proof from simplifier
      switch( simplificationStatus )
      {
      case SPxSimplifier::INFEASIBLE:
      case SPxSimplifier::DUAL_INFEASIBLE:
      case SPxSimplifier::UNBOUNDED:
         MSG_INFO1( spxout, spxout << "simplifier detected infeasibility or unboundedness - solve again without simplifying" << std::endl; )
         _hasBasis = false;
         _preprocessAndSolveReal(false);
         return;

      case SPxSimplifier::VANISHED:
         _status = SPxSolver::OPTIMAL;
         _storeSolutionRealFromPresol();
         return;

      case SPxSimplifier::OKAY:
         _status = _solver.status();
      }

      // process result
      switch( _status )
      {
      case SPxSolver::OPTIMAL:
         _storeSolutionReal(!_isRealLPLoaded || _isRealLPScaled);
         // apply polishing on original problem
         if( _applyPolishing )
         {
            int polishing = intParam(SoPlex::SOLUTION_POLISHING);
            setIntParam(SoPlex::SOLUTION_POLISHING, polishing);
            _preprocessAndSolveReal(false);
         }
         break;

      case SPxSolver::UNBOUNDED:
      case SPxSolver::INFEASIBLE:
      case SPxSolver::INForUNBD:
         // in case of infeasibility or unboundedness, we currently can not unsimplify, but have to solve the original LP again
         if( !_isRealLPLoaded )
         {
            MSG_INFO1( spxout, spxout << " --- loading original problem" << std::endl; )
            _solver.changeObjOffset(realParam(SoPlex::OBJ_OFFSET));
            // we cannot do more to remove violations
            _resolveWithoutPreprocessing(simplificationStatus);
         }
         else
         {
            _storeSolutionReal(false);
         }
         break;

      case SPxSolver::SINGULAR:
         // if preprocessing was applied, try to run again without to avoid singularity
         if( !_isRealLPLoaded )
         {
            MSG_INFO3( spxout, spxout << "encountered singularity - trying to solve again without simplifying" << std::endl; )
            _preprocessAndSolveReal(false);
            return;
         }
         _hasBasis = false;
         break;

      case SPxSolver::ABORT_CYCLING:
         // if preprocessing was applied, try to run again without to avoid cycling
         if( !_isRealLPLoaded )
         {
            MSG_INFO3( spxout, spxout << "encountered cycling - trying to solve again without simplifying" << std::endl; )
            _preprocessAndSolveReal(false);
            return;
         }
         // FALLTHROUGH
      case SPxSolver::ABORT_TIME:
      case SPxSolver::ABORT_ITER:
      case SPxSolver::ABORT_VALUE:
      case SPxSolver::REGULAR:
      case SPxSolver::RUNNING:
         _storeSolutionReal(false);
         break;

      default:
         _hasBasis = false;
         break;
      }
   }



   /// solves real LP with/without preprocessing
   void SoPlex::_preprocessAndSolveReal(bool applySimplifier)
   {
      _solver.changeObjOffset(realParam(SoPlex::OBJ_OFFSET));
      _statistics->preprocessingTime->start();

      _applyPolishing = false;

      if( applySimplifier )
         _enableSimplifierAndScaler();
      else
         _disableSimplifierAndScaler();

      // create a copy of the LP when simplifying or when using internal scaling, i.e. w/o persistent scaling
      bool copyLP = (_simplifier != 0 || (_scaler && !_isRealLPScaled));

      _solver.setTerminationValue(intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE
         ? realParam(SoPlex::OBJLIMIT_UPPER) : realParam(SoPlex::OBJLIMIT_LOWER));

      if( _isRealLPLoaded )
      {
         assert(_realLP == &_solver);

         // preprocessing is always applied to the LP in the solver; hence we have to create a copy of the original LP
         // if simplifier is turned on
         if( copyLP )
         {
            _realLP = 0;
            spx_alloc(_realLP);
            _realLP = new (_realLP) SPxLPReal(_solver);
            _isRealLPLoaded = false;
         }
      }
      else
      {
         assert(_realLP != &_solver);

         // load real LP and basis if available
         if( _hasBasis )
         {
            assert(_basisStatusRows.size() == numRowsReal());
            assert(_basisStatusCols.size() == numColsReal());

            _solver.loadLP(*_realLP, false);
            _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
         }
         // load real LP and set up slack basis
         else
            _solver.loadLP(*_realLP, true);

         // if there is no simplifier, then the original and the transformed problem are identical and it is more
         // memory-efficient to keep only the problem in the solver
         if( !copyLP )
         {
            _realLP->~SPxLPReal();
            spx_free(_realLP);
            _realLP = &_solver;
            _isRealLPLoaded = true;
         }
      }

      // assert that we have two problems if and only if we apply the simplifier
      assert(_realLP == &_solver || copyLP);
      assert(_realLP != &_solver || !copyLP);

      // apply problem simplification
      SPxSimplifier::Result simplificationStatus = SPxSimplifier::OKAY;
      if( _simplifier )
      {
         assert(!_isRealLPLoaded);
         // do not remove bounds of boxed variables or sides of ranged rows if bound flipping is used; also respect row-boundflip parameter
         bool keepbounds = intParam(SoPlex::RATIOTESTER) == SoPlex::RATIOTESTER_BOUNDFLIPPING;
         if( intParam(SoPlex::REPRESENTATION) == SoPlex::REPRESENTATION_ROW
             || (intParam(SoPlex::REPRESENTATION) == SoPlex::REPRESENTATION_AUTO
                 && (_solver.nCols() + 1) * realParam(SoPlex::REPRESENTATION_SWITCH) < (_solver.nRows() + 1)) )
            keepbounds &= boolParam(SoPlex::ROWBOUNDFLIPS);
         simplificationStatus = _simplifier->simplify(_solver, realParam(SoPlex::EPSILON_ZERO), realParam(SoPlex::FEASTOL), realParam(SoPlex::OPTTOL), keepbounds);
         _solver.changeObjOffset(_simplifier->getObjoffset() + realParam(SoPlex::OBJ_OFFSET));
         _solver.setScalingInfo(false);
         _applyPolishing = true;
         _solver.setSolutionPolishing(SPxSolver::POLISH_OFF);
      }

      _statistics->preprocessingTime->stop();

      // run the simplex method if problem has not been solved by the simplifier
      if( simplificationStatus == SPxSimplifier::OKAY )
      {
         if( _scaler && !_solver.isScaled() )
         {
            _scaler->scale(_solver, false);
         }

         _solveRealLPAndRecordStatistics();
      }

      _evaluateSolutionReal(simplificationStatus);
   }



   /// loads original problem into solver and solves again after it has been solved to infeasibility or unboundedness with preprocessing
   void SoPlex::_resolveWithoutPreprocessing(SPxSimplifier::Result simplificationStatus)
   {
      assert(!_isRealLPLoaded || _scaler != 0);
      assert(_simplifier != 0 || _scaler != 0);
      assert(_status == SPxSolver::INFEASIBLE || _status == SPxSolver::INForUNBD || _status == SPxSolver::UNBOUNDED);

      // if simplifier was active, then we unsimplify to get the basis
      if( _simplifier )
      {
         assert(!_simplifier->isUnsimplified());
         assert(simplificationStatus == SPxSimplifier::OKAY);

         // get temporary solution vectors for transformed problem
         DVectorReal primal(_solver.nCols());
         DVectorReal slacks(_solver.nRows());
         DVectorReal dual(_solver.nRows());
         DVectorReal redCost(_solver.nCols());

         _basisStatusRows.reSize(numRowsReal());
         _basisStatusCols.reSize(numColsReal());
         assert(_basisStatusRows.size() >= _solver.nRows());
         assert(_basisStatusCols.size() >= _solver.nCols());

         // get solution data from transformed problem
         _solver.getPrimal(primal);
         _solver.getSlacks(slacks);
         _solver.getDual(dual);
         _solver.getRedCost(redCost);

         // unscale vectors
         if( _scaler && _solver.isScaled())
         {
            _scaler->unscalePrimal(_solver, primal);
            _scaler->unscaleSlacks(_solver, slacks);
            _scaler->unscaleDual(_solver, dual);
            _scaler->unscaleRedCost(_solver, redCost);
         }

         // get basis of transformed problem
         _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(), _basisStatusRows.size(), _basisStatusCols.size());

         try
         {
            _simplifier->unsimplify(primal, dual, slacks, redCost, _basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(), false);
            _simplifier->getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(), _basisStatusRows.size(), _basisStatusCols.size());
            _hasBasis = true;
         }
         catch( const SPxException& E )
         {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> during unsimplification. Resolving without simplifier and scaler.\n" );
            _hasBasis = false;
         }
      }
      // if the original problem is not in the solver because of scaling, we also need to store the basis
      else if( _scaler != 0 )
      {
         _basisStatusRows.reSize(numRowsReal());
         _basisStatusCols.reSize(numColsReal());
         assert(_basisStatusRows.size() == _solver.nRows());
         assert(_basisStatusCols.size() == _solver.nCols());

         _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(), _basisStatusRows.size(), _basisStatusCols.size());
         _hasBasis = true;
      }

      // resolve the original problem
      _preprocessAndSolveReal(false);
      return;
   }



   /// verify computed solution based on status and resolve if claimed primal or dual feasibility is not fulfilled
   void SoPlex::_verifySolutionReal()
   {
      assert(_hasSolReal);
      if( !_solReal._isPrimalFeasible && !_solReal._isDualFeasible )
      {
         _hasSolReal = false;
         return;
      }

      MSG_INFO1( spxout, spxout << " --- verifying computed solution" << std::endl; )

      Real sumviol = 0;
      Real boundviol = 0;
      Real rowviol = 0;
      Real dualviol = 0;
      Real redcostviol = 0;

      if( _solReal._isPrimalFeasible )
      {
         (void) getBoundViolationReal(boundviol, sumviol);
         (void) getRowViolationReal(rowviol, sumviol);
      }
      if( _solReal._isDualFeasible )
      {
         (void) getDualViolationReal(dualviol, sumviol);
         (void) getRedCostViolationReal(redcostviol, sumviol);
      }

      if( boundviol >= _solver.feastol() || rowviol >= _solver.feastol() || dualviol >= _solver.opttol() || redcostviol >= _solver.opttol())
      {
         assert(&_solver == _realLP);
         assert(_isRealLPLoaded);
         MSG_INFO3( spxout, spxout << "bound violation: " << boundviol
                                   << ", row violation: " << rowviol
                                   << ", dual violation: " << dualviol
                                   << ", redcost violation: " << redcostviol << std::endl; )
         MSG_INFO1( spxout, spxout << " --- detected violations in original problem space -- solve again without presolving" << std::endl; )

         if( _isRealLPScaled )
         {
            _solver.unscaleLPandReloadBasis();
            _isRealLPScaled = false;
            ++_unscaleCalls;
         }

         _preprocessAndSolveReal(false);
      }
   }


   /// stores solution data from the solver, possibly after applying unscaling and unsimplifying
   void SoPlex::_storeSolutionReal(bool verify)
   {
      // prepare storage for basis (enough to fit the original basis)
      _basisStatusRows.reSize(numRowsReal());
      _basisStatusCols.reSize(numColsReal());

      // prepare storage for the solution data (only in transformed space due to unscaling), w/o setting it to zero
      _solReal._primal.reDim(_solver.nCols(), false);
      _solReal._slacks.reDim(_solver.nRows(), false);
      _solReal._dual.reDim(_solver.nRows(), false);
      _solReal._redCost.reDim(_solver.nCols(), false);

      // check primal status consistency and query solution status
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::ERROR);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NO_RATIOTESTER);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NO_PRICER);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NO_SOLVER);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NOT_INIT);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::SINGULAR);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NO_PROBLEM);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::UNBOUNDED);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::INFEASIBLE);
      assert(_solver.basis().status() != SPxBasis::UNBOUNDED || status() == SPxSolver::UNBOUNDED);
      assert(_solver.basis().status() == SPxBasis::UNBOUNDED || _solver.basis().status() == SPxBasis::NO_PROBLEM || status() != SPxSolver::UNBOUNDED);

      _solReal._isPrimalFeasible = (status() == SPxSolver::OPTIMAL
         || ((_solver.basis().status() == SPxBasis::PRIMAL || _solver.basis().status() == SPxBasis::UNBOUNDED)
            && _solver.shift() < 10.0 * realParam(SoPlex::EPSILON_ZERO)));

      _solReal._hasPrimalRay = (status() == SPxSolver::UNBOUNDED && _isRealLPLoaded);

      // check dual status consistency and query solution status
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::ERROR);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NO_RATIOTESTER);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NO_PRICER);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NO_SOLVER);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NOT_INIT);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::SINGULAR);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NO_PROBLEM);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::UNBOUNDED);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::INFEASIBLE);
      assert(_solver.basis().status() != SPxBasis::INFEASIBLE || status() == SPxSolver::INFEASIBLE);
      assert(_solver.basis().status() == SPxBasis::INFEASIBLE || _solver.basis().status() == SPxBasis::NO_PROBLEM || status() != SPxSolver::INFEASIBLE);

      _solReal._isDualFeasible = (status() == SPxSolver::OPTIMAL
         || ((_solver.basis().status() == SPxBasis::DUAL || _solver.basis().status() == SPxBasis::INFEASIBLE)
            && _solver.shift() < 10.0 * realParam(SoPlex::EPSILON_ZERO)));

      _solReal._hasDualFarkas = (status() == SPxSolver::INFEASIBLE && _isRealLPLoaded);

      // get infeasibility or unboundedness proof if available
      if( _solReal._hasPrimalRay )
      {
         _solReal._primalRay.reDim(_solver.nCols(), false);
         _solver.getPrimalray(_solReal._primalRay);
      }

      if( _solReal._hasDualFarkas )
      {
         _solReal._dualFarkas.reDim(_solver.nRows(), false);
         _solver.getDualfarkas(_solReal._dualFarkas);
      }

      // get solution data from the solver; independent of solution status
      _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(),
                       _basisStatusRows.size(), _basisStatusCols.size());

      _solver.getPrimal(_solReal._primal);
      _solver.getSlacks(_solReal._slacks);
      _solver.getDual(_solReal._dual);
      _solver.getRedCost(_solReal._redCost);

      _hasBasis = true;

      // get primal and/or dual objective function value depending on status
      _solver.forceRecompNonbasicValue();
      _solReal._objVal = _solver.objValue();

      // infeasible solutions shall also be stored and be accessible
      _hasSolReal = true;

      // unscale vectors
      if( _solver.isScaled() && !_isRealLPLoaded )
         _unscaleSolutionReal(_solver, false);

      // get unsimplified solution data from simplifier
      if( _simplifier )
      {
         assert(!_simplifier->isUnsimplified());
         assert(_simplifier->result() == SPxSimplifier::OKAY);
         assert(_realLP != &_solver);

         try
         {
            // pass solution data of transformed problem to simplifier
            _simplifier->unsimplify(_solReal._primal, _solReal._dual, _solReal._slacks, _solReal._redCost,
                                    _basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(), status() == SPxSolver::OPTIMAL);
         }
         catch( const SPxException& E )
         {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> during unsimplification. Resolving without simplifier and scaler.\n" );
            _hasBasis = false;
            _preprocessAndSolveReal(false);
            return;
         }

         // copy unsimplified solution data from simplifier (size and dimension is adapted during copy)
         _solReal._primal  = _simplifier->unsimplifiedPrimal();
         _solReal._slacks  = _simplifier->unsimplifiedSlacks();
         _solReal._dual    = _simplifier->unsimplifiedDual();
         _solReal._redCost = _simplifier->unsimplifiedRedCost();

         // overwrite the transformed basis with the unsimplified one
         _simplifier->getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(), _basisStatusRows.size(), _basisStatusCols.size());

         // load original problem but don't setup a slack basis
         _loadRealLP(false);

         assert(_realLP == &_solver);

         // load unsimplified basis into solver
         assert(_basisStatusRows.size() == numRowsReal());
         assert(_basisStatusCols.size() == numColsReal());
         _solver.setBasisStatus(SPxBasis::REGULAR);
         _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
         _hasBasis = true;
      }
      // load realLP into the solver again (internal scaling was applied)
      else if( _realLP != &_solver )
      {
         assert(_solver.isScaled());
         _loadRealLP(false);
      }

      // unscale stored solution (removes persistent scaling)
      if( _isRealLPScaled )
         _unscaleSolutionReal(*_realLP, true);

      // check solution for violations and solve again if necessary
      if( verify )
         _verifySolutionReal();

      assert(_solver.nCols() == numColsReal());
      assert(_solver.nRows() == numRowsReal());
   }



   void SoPlex::_storeSolutionRealFromPresol()
   {
      assert(_simplifier);
      assert(_simplifier->result() == SPxSimplifier::VANISHED);

      // prepare storage for basis (enough to fit the original basis)
      _basisStatusRows.reSize(numRowsReal());
      _basisStatusCols.reSize(numColsReal());

      // prepare storage for the solution data and initialize it to zero
      _solReal._primal.reDim(numColsReal(), true);
      _solReal._slacks.reDim(numRowsReal(), true);
      _solReal._dual.reDim(numRowsReal(), true);
      _solReal._redCost.reDim(numColsReal(), true);

      // load original LP and setup slack basis for unsimplifying
      _loadRealLP(true);

      // store slack basis
      _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(),
                       _basisStatusRows.size(), _basisStatusCols.size());

      assert(!_simplifier->isUnsimplified());

      try
      {
         // unsimplify basis and solution data
         _simplifier->unsimplify(_solReal._primal, _solReal._dual, _solReal._slacks, _solReal._redCost,
                                 _basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());

      }
      catch( const SPxException& E )
      {
         MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> during unsimplification. Resolving without simplifier and scaler.\n" );
         _preprocessAndSolveReal(false);
         return;
      }

      // copy unsimplified solution data from simplifier
      _solReal._primal  = _simplifier->unsimplifiedPrimal();
      _solReal._slacks  = _simplifier->unsimplifiedSlacks();
      _solReal._dual    = _simplifier->unsimplifiedDual();
      _solReal._redCost = _simplifier->unsimplifiedRedCost();

      // unscale stored solution (removes persistent scaling)
      if( _isRealLPScaled )
         _unscaleSolutionReal(*_realLP, true);

      // compute the original objective function value
      _solReal._objVal = realParam(SoPlex::OBJ_OFFSET);
      for( int i = 0; i < numColsReal(); ++i )
         _solReal._objVal += _solReal._primal[i] * objReal(i);

      // store the unsimplified basis
      _simplifier->getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(), _basisStatusRows.size(), _basisStatusCols.size());
      _hasBasis = true;
      _hasSolReal = true;
      _solReal._isPrimalFeasible = true;
      _solReal._isDualFeasible = true;

      // check solution for violations and solve again if necessary
      _verifySolutionReal();
   }



   /// load original LP and possibly setup a slack basis
   void SoPlex::_loadRealLP(bool initBasis)
   {
      _solver.loadLP(*_realLP, initBasis);
      _isRealLPLoaded = true;
      _realLP->~SPxLPReal();
      spx_free(_realLP);
      _realLP = &_solver;
      if( initBasis )
         _solver.init();
   }



   /// unscales stored solution to remove internal or external scaling of LP
   void SoPlex::_unscaleSolutionReal(SPxLPReal& LP, bool persistent)
   {
      MSG_INFO1( spxout, spxout << " --- unscaling " << (persistent ? "external" : "internal") <<" solution" << std::endl; )
      assert(_scaler);
      assert(!persistent || (boolParam(SoPlex::PERSISTENTSCALING) && _isRealLPScaled));
      _scaler->unscalePrimal(LP, _solReal._primal);
      _scaler->unscaleSlacks(LP, _solReal._slacks);
      _scaler->unscaleDual(LP, _solReal._dual);
      _scaler->unscaleRedCost(LP, _solReal._redCost);
      if( _solReal.hasPrimalRay() )
         _scaler->unscalePrimalray(LP, _solReal._primalRay);
      if( _solReal.hasDualFarkas() )
         _scaler->unscaleDualray(LP, _solReal._dualFarkas);
   }
} // namespace soplex
#endif
