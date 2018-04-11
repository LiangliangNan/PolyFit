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

/**@file  soplex.h
 * @brief Preconfigured SoPlex LP solver
 */

#ifndef _SOPLEX_H_
#define _SOPLEX_H_

#ifndef SOPLEX_LEGACY
#include <string.h>

///@todo SoPlex should also have an spxout object to avoid using a global one
#include "spxdefines.h"
#include "basevectors.h"
#include "spxsolver.h"
#include "slufactor.h"
#include "slufactor_rational.h"

///@todo try to move to cpp file by forward declaration
#include "spxsimplifier.h"
#include "spxmainsm.h"

#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxleastsqsc.h"
#include "spxgeometsc.h"

#include "spxstarter.h"
#include "spxweightst.h"
#include "spxsumst.h"
#include "spxvectorst.h"

#include "spxpricer.h"
#include "spxautopr.h"
#include "spxdantzigpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxsteeppr.h"
#include "spxsteepexpr.h"
#include "spxhybridpr.h"

#include "spxratiotester.h"
#include "spxdefaultrt.h"
#include "spxharrisrt.h"
#include "spxfastrt.h"
#include "spxboundflippingrt.h"

#include "sol.h"

#define DEFAULT_RANDOM_SEED   0   // used to suppress output when the seed was not changed

///@todo implement automatic rep switch, based on row/col dim
///@todo introduce status codes for SoPlex, especially for rational solving

///@todo record and return "best" solutions found during IR (Ambros)
///@todo implement main IR loop for primal and dual feasible case with fail otherwise (Ambros)
///@todo implement statistical info (time, factor time, iters, ...) since last call to solveReal() or solveRational() (Ambros?)
///@todo implement performInfeasibilityIR (Ambros?)
///@todo extend IR loop to infeasible case (Dan?)
///@todo extend IR loop to unbounded case (Dan?)

///@todo interface rational reconstruction code for rational vectors
///@todo integrate rational reconstruction into IR loop
///@todo templatize SPxSolver and necessary components (SLUFactor, pricer, ratiotester)
///@todo integrate rational SPxSolver and distinguish between original and transformed rational LP
///@todo rational scalers
///@todo rational simplifier

namespace soplex
{

/**@class SoPlex
 * @brief   Preconfigured SoPlex LP-solver.
 * @ingroup Algo
 */
class SoPlex
{
public:

   //**@name Construction and destruction */
   //@{

   /// default constructor
   SoPlex();

   /// assignment operator
   SoPlex& operator=(const SoPlex& rhs);

   /// copy constructor
   SoPlex(const SoPlex& rhs);

   /// destructor
   virtual ~SoPlex();

   //@}


   //**@name Access to the real LP */
   //@{

   /// returns number of rows
   int numRowsReal() const;

   /// returns number of columns
   int numColsReal() const;

   /// returns number of nonzeros
   int numNonzerosReal() const;

   /// returns smallest non-zero element in absolute value
   Real minAbsNonzeroReal() const;

   /// returns biggest non-zero element in absolute value
   Real maxAbsNonzeroReal() const;

   /// returns (unscaled) coefficient
   Real coefReal(int row, int col) const;

   /// returns vector of row \p i, ignoring scaling
   const SVectorReal& rowVectorRealInternal(int i) const;

   /// gets vector of row \p i
   void getRowVectorReal(int i, DSVectorReal& row) const;

   /// returns right-hand side vector, ignoring scaling
   const VectorReal& rhsRealInternal() const;

   /// gets right-hand side vector
   void getRhsReal(DVectorReal& rhs) const;

   /// returns right-hand side of row \p i
   Real rhsReal(int i) const;

   /// returns left-hand side vector, ignoring scaling
   const VectorReal& lhsRealInternal() const;

   /// gets left-hand side vector
   void getLhsReal(DVectorReal& lhs) const;

   /// returns left-hand side of row \p i
   Real lhsReal(int i) const;

   /// returns inequality type of row \p i
   LPRowReal::Type rowTypeReal(int i) const;

   /// returns vector of col \p i, ignoring scaling
   const SVectorReal& colVectorRealInternal(int i) const;

   /// gets vector of col \p i
   void getColVectorReal(int i, DSVectorReal& col) const;

   /// returns upper bound vector
   const VectorReal& upperRealInternal() const;

   /// returns upper bound of column \p i
   Real upperReal(int i) const;

   /// gets upper bound vector
   void getUpperReal(DVectorReal& upper) const;

   /// returns lower bound vector
   const VectorReal& lowerRealInternal() const;

   /// returns lower bound of column \p i
   Real lowerReal(int i) const;

   /// gets lower bound vector
   void getLowerReal(DVectorReal& lower) const;

   /// gets objective function vector
   void getObjReal(VectorReal& obj) const;

   /// returns objective value of column \p i
   Real objReal(int i) const;

   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorReal& maxObjRealInternal() const;

   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   Real maxObjReal(int i) const;

   /// gets number of available dual norms
   void getNdualNorms(int& nnormsRow, int& nnormsCol) const;

   /// gets steepest edge norms and returns false if they are not available
   bool getDualNorms(int& nnormsRow, int& nnormsCol, Real* norms) const;

   /// sets steepest edge norms and returns false if that's not possible
   bool setDualNorms(int nnormsRow, int nnormsCol, Real* norms);

   /// pass integrality information about the variables to the solver
   void setIntegralityInformation(int ncols, int* intInfo);

   //@}


   //**@name Access to the rational LP */
   //@{

   /// returns number of rows
   int numRowsRational() const;

   /// returns number of columns
   int numColsRational() const;

   /// returns number of nonzeros
   int numNonzerosRational() const;

   /// returns smallest non-zero element in absolute value
   Rational minAbsNonzeroRational() const;

   /// returns biggest non-zero element in absolute value
   Rational maxAbsNonzeroRational() const;

   /// gets row \p i
   void getRowRational(int i, LPRowRational& lprow) const;

   /// gets rows \p start, ..., \p end.
   void getRowsRational(int start, int end, LPRowSetRational& lprowset) const;

   /// returns vector of row \p i
   const SVectorRational& rowVectorRational(int i) const;

   /// returns right-hand side vector
   const VectorRational& rhsRational() const;

   /// returns right-hand side of row \p i
   const Rational& rhsRational(int i) const;

   /// returns left-hand side vector
   const VectorRational& lhsRational() const;

   /// returns left-hand side of row \p i
   const Rational& lhsRational(int i) const;

   /// returns inequality type of row \p i
   LPRowRational::Type rowTypeRational(int i) const;

   /// gets column \p i
   void getColRational(int i, LPColRational& lpcol) const;

   /// gets columns \p start, ..., \p end
   void getColsRational(int start, int end, LPColSetRational& lpcolset) const;

   /// returns vector of column \p i
   const SVectorRational& colVectorRational(int i) const;

   /// returns upper bound vector
   const VectorRational& upperRational() const;

   /// returns upper bound of column \p i
   const Rational& upperRational(int i) const;

   /// returns lower bound vector
   const VectorRational& lowerRational() const;

   /// returns lower bound of column \p i
   const Rational& lowerRational(int i) const;

   /// gets objective function vector
   void getObjRational(VectorRational& obj) const;

   /// gets objective value of column \p i
   void getObjRational(int i, Rational& obj) const;

   /// returns objective value of column \p i
   Rational objRational(int i) const;

   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorRational& maxObjRational() const;

   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   const Rational& maxObjRational(int i) const;

   //@}


   //**@name Modification of the real LP */
   //@{

   /// adds a single row
   void addRowReal(const LPRowReal& lprow);

   /// adds multiple rows
   void addRowsReal(const LPRowSetReal& lprowset);

   /// adds a single column
   void addColReal(const LPCol& lpcol);

   /// adds multiple columns
   void addColsReal(const LPColSetReal& lpcolset);

   /// replaces row \p i with \p lprow
   void changeRowReal(int i, const LPRowReal& lprow);

   /// changes left-hand side vector for constraints to \p lhs
   void changeLhsReal(const VectorReal& lhs);

   /// changes left-hand side of row \p i to \p lhs
   void changeLhsReal(int i, const Real& lhs);

   /// changes right-hand side vector to \p rhs
   void changeRhsReal(const VectorReal& rhs);

   /// changes right-hand side of row \p i to \p rhs
   void changeRhsReal(int i, const Real& rhs);

   /// changes left- and right-hand side vectors
   void changeRangeReal(const VectorReal& lhs, const VectorReal& rhs);

   /// changes left- and right-hand side of row \p i
   void changeRangeReal(int i, const Real& lhs, const Real& rhs);

   /// replaces column \p i with \p lpcol
   void changeColReal(int i, const LPColReal& lpcol);

   /// changes vector of lower bounds to \p lower
   void changeLowerReal(const VectorReal& lower);

   /// changes lower bound of column i to \p lower
   void changeLowerReal(int i, const Real& lower);

   /// changes vector of upper bounds to \p upper
   void changeUpperReal(const VectorReal& upper);

   /// changes \p i 'th upper bound to \p upper
   void changeUpperReal(int i, const Real& upper);

   /// changes vectors of column bounds to \p lower and \p upper
   void changeBoundsReal(const VectorReal& lower, const VectorReal& upper);

   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsReal(int i, const Real& lower, const Real& upper);

   /// changes objective function vector to \p obj
   void changeObjReal(const VectorReal& obj);

   /// changes objective coefficient of column i to \p obj
   void changeObjReal(int i, const Real& obj);

   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementReal(int i, int j, const Real& val);

   /// removes row \p i
   void removeRowReal(int i);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void removeRowsReal(int perm[]);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsReal() may be passed
   /// as buffer memory
   void removeRowsReal(int idx[], int n, int perm[] = 0);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsReal() may be passed as buffer
   /// memory
   void removeRowRangeReal(int start, int end, int perm[] = 0);

   /// removes column i
   void removeColReal(int i);

   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void removeColsReal(int perm[]);

   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void removeColsReal(int idx[], int n, int perm[] = 0);

   /// removes columns \p start to \p end including both; an array \p perm of size #numColsReal() may be passed as
   /// buffer memory
   void removeColRangeReal(int start, int end, int perm[] = 0);

   /// clears the LP
   void clearLPReal();

   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, if sync mode is manual
   void syncLPReal();

   //@}


   //**@name Modification of the rational LP */
   //@{

   /// adds a single row
   void addRowRational(const LPRowRational& lprow);

#ifdef SOPLEX_WITH_GMP
   /// adds a single row
   void addRowRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int rowSize, const mpq_t* rhs);

   /// adds a set of rows
   void addRowsRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int* rowStarts, const int* rowLengths, const int numRows, const int numValues, const mpq_t* rhs);
#endif

   /// adds multiple rows
   void addRowsRational(const LPRowSetRational& lprowset);

   /// adds a single column
   void addColRational(const LPColRational& lpcol);

#ifdef SOPLEX_WITH_GMP
   /// adds a single column
   void addColRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int colSize, const mpq_t* upper);

   /// adds a set of columns
   void addColsRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int* colStarts, const int* colLengths, const int numCols, const int numValues, const mpq_t* upper);
#endif

   /// adds multiple columns
   void addColsRational(const LPColSetRational& lpcolset);

   /// replaces row \p i with \p lprow
   void changeRowRational(int i, const LPRowRational& lprow);

   /// changes left-hand side vector for constraints to \p lhs
   void changeLhsRational(const VectorRational& lhs);

   /// changes left-hand side of row \p i to \p lhs
   void changeLhsRational(int i, const Rational& lhs);

#ifdef SOPLEX_WITH_GMP
   /// changes left-hand side of row \p i to \p lhs
   void changeLhsRational(int i, const mpq_t* lhs);
#endif

   /// changes right-hand side vector to \p rhs
   void changeRhsRational(const VectorRational& rhs);

#ifdef SOPLEX_WITH_GMP
   /// changes right-hand side vector to \p rhs
   void changeRhsRational(const mpq_t* rhs, int rhsSize);
#endif

   /// changes right-hand side of row \p i to \p rhs
   void changeRhsRational(int i, const Rational& rhs);

   /// changes left- and right-hand side vectors
   void changeRangeRational(const VectorRational& lhs, const VectorRational& rhs);

   /// changes left- and right-hand side of row \p i
   void changeRangeRational(int i, const Rational& lhs, const Rational& rhs);

#ifdef SOPLEX_WITH_GMP
   /// changes left- and right-hand side of row \p i
   void changeRangeRational(int i, const mpq_t* lhs, const mpq_t* rhs);
#endif

   /// replaces column \p i with \p lpcol
   void changeColRational(int i, const LPColRational& lpcol);

   /// changes vector of lower bounds to \p lower
   void changeLowerRational(const VectorRational& lower);

   /// changes lower bound of column i to \p lower
   void changeLowerRational(int i, const Rational& lower);

#ifdef SOPLEX_WITH_GMP
   /// changes lower bound of column i to \p lower
   void changeLowerRational(int i, const mpq_t* lower);
#endif

   /// changes vector of upper bounds to \p upper
   void changeUpperRational(const VectorRational& upper);

   /// changes \p i 'th upper bound to \p upper
   void changeUpperRational(int i, const Rational& upper);

#ifdef SOPLEX_WITH_GMP
   /// changes upper bound of column i to \p upper
   void changeUpperRational(int i, const mpq_t* upper);
#endif

   /// changes vectors of column bounds to \p lower and \p upper
   void changeBoundsRational(const VectorRational& lower, const VectorRational& upper);

   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsRational(int i, const Rational& lower, const Rational& upper);

#ifdef SOPLEX_WITH_GMP
   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsRational(int i, const mpq_t* lower, const mpq_t* upper);
#endif

   /// changes objective function vector to \p obj
   void changeObjRational(const VectorRational& obj);

   /// changes objective coefficient of column i to \p obj
   void changeObjRational(int i, const Rational& obj);

#ifdef SOPLEX_WITH_GMP
   /// changes objective coefficient of column i to \p obj
   void changeObjRational(int i, const mpq_t* obj);
#endif

   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementRational(int i, int j, const Rational& val);

#ifdef SOPLEX_WITH_GMP
   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementRational(int i, int j, const mpq_t* val);
#endif

   /// removes row \p i
   void removeRowRational(int i);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the new
   /// index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsRational()
   void removeRowsRational(int perm[]);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsRational() may be
   /// passed as buffer memory
   void removeRowsRational(int idx[], int n, int perm[] = 0);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsRational() may be passed as
   /// buffer memory
   void removeRowRangeRational(int start, int end, int perm[] = 0);

   /// removes column i
   void removeColRational(int i);

   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsRational()
   void removeColsRational(int perm[]);

   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsRational() may be
   /// passed as buffer memory
   void removeColsRational(int idx[], int n, int perm[] = 0);

   /// removes columns \p start to \p end including both; an array \p perm of size #numColsRational() may be passed as
   /// buffer memory
   void removeColRangeRational(int start, int end, int perm[] = 0);

   /// clears the LP
   void clearLPRational();

   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, if sync mode is manual
   void syncLPRational();

   //@}


   //**@name Solving and general solution query */
   //@{

   /// optimize the given LP
   SPxSolver::Status optimize();

   // old name for backwards compatibility
   SPxSolver::Status solve()
   {
      return optimize();
   }

   /// returns the current solver status
   SPxSolver::Status status() const;

   /// is stored primal solution feasible?
   bool isPrimalFeasible() const;

   /// is a primal feasible solution available?
   bool hasPrimal() const;

   /// is a primal unbounded ray available?
   bool hasPrimalRay() const;

   /// is stored dual solution feasible?
   bool isDualFeasible() const;

   /// is a dual feasible solution available?
   bool hasDual() const;

   /// is Farkas proof of infeasibility available?
   bool hasDualFarkas() const;

   //@}


   //**@name Query for the real solution data */
   //@{

   /// returns the objective value if a primal solution is available
   Real objValueReal();

   /// gets the primal solution vector if available; returns true on success
   bool getPrimalReal(VectorReal& vector);

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksReal(VectorReal& vector);

   /// gets the primal ray if available; returns true on success
   bool getPrimalRayReal(VectorReal& vector);

   /// gets the dual solution vector if available; returns true on success
   bool getDualReal(VectorReal& vector);

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedCostReal(VectorReal& vector);

   /// gets the Farkas proof if available; returns true on success
   bool getDualFarkasReal(VectorReal& vector);

   /// gets violation of bounds; returns true on success
   bool getBoundViolationReal(Real& maxviol, Real& sumviol);

   /// gets violation of constraints; returns true on success
   bool getRowViolationReal(Real& maxviol, Real& sumviol);

   /// gets violation of reduced costs; returns true on success
   bool getRedCostViolationReal(Real& maxviol, Real& sumviol);

   /// gets violation of dual multipliers; returns true on success
   bool getDualViolationReal(Real& maxviol, Real& sumviol);

   //@}


   //**@name Query for the rational solution data */
   //@{

   /// returns the objective value if a primal solution is available
   Rational objValueRational();

   /// gets the primal solution vector if available; returns true on success
   bool getPrimalRational(VectorRational& vector);

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksRational(VectorRational& vector);

   /// gets the primal ray if LP is unbounded; returns true on success
   bool getPrimalRayRational(VectorRational& vector);

   /// gets the dual solution vector if available; returns true on success
   bool getDualRational(VectorRational& vector);

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedCostRational(VectorRational& vector);

   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool getDualFarkasRational(VectorRational& vector);

   /// gets violation of bounds; returns true on success
   bool getBoundViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of constraints; returns true on success
   bool getRowViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of reduced costs; returns true on success
   bool getRedCostViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of dual multipliers; returns true on success
   bool getDualViolationRational(Rational& maxviol, Rational& sumviol);

#ifdef SOPLEX_WITH_GMP
   /// gets the primal solution vector if available; returns true on success
   bool getPrimalRational(mpq_t* vector, const int size);

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksRational(mpq_t* vector, const int size);

   /// gets the primal ray if LP is unbounded; returns true on success
   bool getPrimalRayRational(mpq_t* vector, const int size);

   /// gets the dual solution vector if available; returns true on success
   bool getDualRational(mpq_t* vector, const int size);

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedCostRational(mpq_t* vector, const int size);

   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool getDualFarkasRational(mpq_t* vector, const int size);
#endif

   /// get size of primal solution
   int totalSizePrimalRational(const int base = 2);

   /// get size of dual solution
   int totalSizeDualRational(const int base = 2);

   /// get size of least common multiple of denominators in primal solution
   int dlcmSizePrimalRational(const int base = 2);

   /// get size of least common multiple of denominators in dual solution
   int dlcmSizeDualRational(const int base = 2);

   /// get size of largest denominator in primal solution
   int dmaxSizePrimalRational(const int base = 2);

   /// get size of largest denominator in dual solution
   int dmaxSizeDualRational(const int base = 2);

   //@}


   //**@name Access and modification of basis information */
   //@{

   /// is an advanced starting basis available?
   bool hasBasis() const;

   /// returns the current basis status
   SPxBasis::SPxStatus basisStatus() const;

   /// returns basis status for a single row
   SPxSolver::VarStatus basisRowStatus(int row) const;

   /// returns basis status for a single column
   SPxSolver::VarStatus basisColStatus(int col) const;

   /// gets current basis via arrays of statuses
   void getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const;

   /// gets the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
   void getBasisInd(int* bind) const;

   /// compute condition number estimate based on the diagonal of the LU factorization; returns true on success
   /// type = 0: max/min ratio
   /// type = 1: trace of U (sum of diagonal elements)
   /// type = 2: product of diagonal elements
   bool getFastCondition(Real& condition, int type = 0);

   /// computes an estimated condition number for the current basis matrix using the power method; returns true on success
   bool getEstimatedCondition(Real& condition);

   /// computes the exact condition number for the current basis matrix using the power method; returns true on success
   bool getExactCondition(Real& condition);

   /// computes row \p r of basis inverse; returns true on success
   /// @param r which row of the basis inverse is computed
   /// @param coef values of result vector (not packed but scattered)
   /// @param inds indices of result vector (NULL if not to be used)
   /// @param ninds number of nonzeros in result vector
   /// @param unscale determines whether the result should be unscaled according to the original LP data
   bool getBasisInverseRowReal(int r, Real* coef, int* inds = NULL, int* ninds = NULL, bool unscale = true);

   /// computes column \p c of basis inverse; returns true on success
   /// @param c which column of the basis inverse is computed
   /// @param coef values of result vector (not packed but scattered)
   /// @param inds indices of result vector (NULL if not to be used)
   /// @param ninds number of nonzeros in result vector
   /// @param unscale determines whether the result should be unscaled according to the original LP data
   bool getBasisInverseColReal(int c, Real* coef, int* inds = NULL, int* ninds = NULL, bool unscale = true);

   /// computes dense solution of basis matrix B * \p sol = \p rhs; returns true on success
   bool getBasisInverseTimesVecReal(Real* rhs, Real* sol, bool unscale = true);

   /// multiply with basis matrix; B * \p vec (inplace)
   /// @param vec (dense) vector to be multiplied with
   /// @param unscale determines whether the result should be unscaled according to the original LP data
   bool multBasis(Real* vec, bool unscale = true);

   /// multiply with transpose of basis matrix; \p vec * B^T (inplace)
   /// @param vec (dense) vector to be multiplied with
   /// @param unscale determines whether the result should be unscaled according to the original LP data
   bool multBasisTranspose(Real* vec, bool unscale = true);

   /// compute rational basis inverse; returns true on success
   bool computeBasisInverseRational();

   /// gets an array of indices for the columns of the rational basis matrix; bind[i] >= 0 means that the i-th column of
   /// the basis matrix contains variable bind[i]; bind[i] < 0 means that the i-th column of the basis matrix contains
   /// the slack variable for row -bind[i]-1; performs rational factorization if not available; returns true on success
   bool getBasisIndRational(DataArray<int>& bind);

   /// computes row r of basis inverse; performs rational factorization if not available; returns true on success
   bool getBasisInverseRowRational(const int r, SSVectorRational& vec);

   /// computes column c of basis inverse; performs rational factorization if not available; returns true on success
   bool getBasisInverseColRational(const int c, SSVectorRational& vec);

   /// computes solution of basis matrix B * sol = rhs; performs rational factorization if not available; returns true
   /// on success
   bool getBasisInverseTimesVecRational(const SVectorRational& rhs, SSVectorRational& sol);

   /// sets starting basis via arrays of statuses
   void setBasis(const SPxSolver::VarStatus rows[], const SPxSolver::VarStatus cols[]);

   /// clears starting basis
   void clearBasis();

   //@}


   //**@name Statistical information */
   //@{

   /// number of iterations since last call to solve
   int numIterations() const;

   /// time spent in last call to solve
   Real solveTime() const;

   /// statistical information in form of a string
   std::string statisticString() const;

   /// name of starter
   const char* getStarterName();

   /// name of simplifier
   const char* getSimplifierName();

   /// name of scaling method
   const char* getScalerName();

   /// name of currently loaded pricer
   const char* getPricerName();

   /// name of currently loaded ratiotester
   const char* getRatiotesterName();

   //@}


   //**@name File I/O */
   //@{

   /// reads LP file in LP or MPS format according to READMODE parameter; gets row names, column names, and
   /// integer variables if desired; returns true on success
   bool readFile(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool writeFileReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* intvars = 0, const bool unscale = true) const;

   /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool writeFileRational(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* intvars = 0) const;

   /// writes the dual of the real LP to file; LP or MPS format is chosen from the extension in \p filename;
   /// if \p rowNames and \p colNames are \c NULL, default names are used; if \p intVars is not \c NULL,
   /// the variables contained in it are marked as integer; returns true on success
   bool writeDualFileReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* intvars = 0) const;

   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed; returns true on success
   bool readBasisFile(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0);

   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used;
   /// returns true on success
   bool writeBasisFile(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const bool cpxFormat = false) const;

   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void writeStateReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const bool cpxFormat = false) const;

   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void writeStateRational(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const bool cpxFormat = false) const;

   //@}


   //**@name Parameters */
   //@{

   /// boolean parameters
   typedef enum
   {
      /// should lifting be used to reduce range of nonzero matrix coefficients?
      LIFTING = 0,

      /// should LP be transformed to equality form before a rational solve?
      EQTRANS = 1,

      /// should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?
      TESTDUALINF = 2,

      /// should a rational factorization be performed after iterative refinement?
      RATFAC = 3,

      /// should the decomposition based dual simplex be used to solve the LP? Setting this to true forces the solve mode to
      /// SOLVEMODE_REAL and the basis representation to REPRESENTATION_ROW
      USEDECOMPDUALSIMPLEX = 4,

      /// should the degeneracy be computed for each basis?
      COMPUTEDEGEN = 5,

      /// should the dual of the complementary problem be used in the decomposition simplex?
      USECOMPDUAL = 6,

      /// should row and bound violations be computed explicitly in the update of reduced problem in the decomposition simplex
      EXPLICITVIOL = 7,

      /// should cycling solutions be accepted during iterative refinement?
      ACCEPTCYCLING = 8,

      /// apply rational reconstruction after each iterative refinement?
      RATREC = 9,

      /// round scaling factors for iterative refinement to powers of two?
      POWERSCALING = 10,

      /// continue iterative refinement with exact basic solution if not optimal?
      RATFACJUMP = 11,

      /// use bound flipping also for row representation?
      ROWBOUNDFLIPS = 12,

      /// use persistent scaling?
      PERSISTENTSCALING = 13,

      /// perturb the entire problem or only the relevant bounds of s single pivot?
      FULLPERTURBATION = 14,

      /// number of boolean parameters
      BOOLPARAM_COUNT = 15
   } BoolParam;

   /// integer parameters
   typedef enum
   {
      /// objective sense
      OBJSENSE = 0,

      /// type of computational form, i.e., column or row representation
      REPRESENTATION = 1,

      /// type of algorithm, i.e., primal or dual
      ALGORITHM = 2,

      /// type of LU update
      FACTOR_UPDATE_TYPE = 3,

      /// maximum number of updates without fresh factorization
      FACTOR_UPDATE_MAX = 4,

      /// iteration limit (-1 if unlimited)
      ITERLIMIT = 5,

      /// refinement limit (-1 if unlimited)
      REFLIMIT = 6,

      /// stalling refinement limit (-1 if unlimited)
      STALLREFLIMIT = 7,

      /// display frequency
      DISPLAYFREQ = 8,

      /// verbosity level
      VERBOSITY = 9,

      /// type of simplifier
      SIMPLIFIER = 10,

      /// type of scaler
      SCALER = 11,

      /// type of starter used to create crash basis
      STARTER = 12,

      /// type of pricer
      PRICER = 13,

      /// type of ratio test
      RATIOTESTER = 14,

      /// mode for synchronizing real and rational LP
      SYNCMODE = 15,

      /// mode for reading LP files
      READMODE = 16,

      /// mode for iterative refinement strategy
      SOLVEMODE = 17,

      /// mode for a posteriori feasibility checks
      CHECKMODE = 18,

      /// type of timer
      TIMER = 19,

      /// mode for hyper sparse pricing
      HYPER_PRICING = 20,

      /// minimum number of stalling refinements since last pivot to trigger rational factorization
      RATFAC_MINSTALLS = 21,

      /// maximum number of conjugate gradient iterations in least square scaling
      LEASTSQ_MAXROUNDS = 22,

      /// mode for solution polishing
      SOLUTION_POLISHING = 23,

      /// the number of iterations before the decomposition simplex initialisation is terminated.
      DECOMP_ITERLIMIT = 24,

      /// the maximum number of rows that are added in each iteration of the decomposition based simplex
      DECOMP_MAXADDEDROWS = 25,

      /// the iteration frequency at which the decomposition solve output is displayed.
      DECOMP_DISPLAYFREQ = 26,

      /// the verbosity of the decomposition based simplex
      DECOMP_VERBOSITY = 27,

      /// print condition number during the solve
      PRINTCONDITION = 28,

      /// number of integer parameters
      INTPARAM_COUNT = 29
   } IntParam;

   /// values for parameter OBJSENSE
   enum
   {
      /// minimization
      OBJSENSE_MINIMIZE = -1,

      /// maximization
      OBJSENSE_MAXIMIZE = 1
   };

   /// values for parameter REPRESENTATION
   enum
   {
      /// automatic choice according to number of rows and columns
      REPRESENTATION_AUTO = 0,

      /// column representation Ax - s = 0, lower <= x <= upper, lhs <= s <= rhs
      REPRESENTATION_COLUMN = 1,

      /// row representation (lower,lhs) <= (x,Ax) <= (upper,rhs)
      REPRESENTATION_ROW = 2
   };

   /// values for parameter ALGORITHM
   enum
   {
      /// primal simplex algorithm, i.e., entering for column and leaving for row representation
      ALGORITHM_PRIMAL = 0,

      /// dual simplex algorithm, i.e., leaving for column and entering for row representation
      ALGORITHM_DUAL = 1
   };

   /// values for parameter FACTOR_UPDATE_TYPE
   enum
   {
      /// product form update
      FACTOR_UPDATE_TYPE_ETA = 0,

      /// Forrest-Tomlin type update
      FACTOR_UPDATE_TYPE_FT = 1
   };

   /// values for parameter VERBOSITY
   enum
   {
      /// only error output
      VERBOSITY_ERROR = 0,

      /// only error and warning output
      VERBOSITY_WARNING = 1,

      /// only error, warning, and debug output
      VERBOSITY_DEBUG = 2,

      /// standard verbosity level
      VERBOSITY_NORMAL = 3,

      /// high verbosity level
      VERBOSITY_HIGH = 4,

      /// full verbosity level
      VERBOSITY_FULL = 5
   };

   /// values for parameter SIMPLIFIER
   enum
   {
      /// no simplifier
      SIMPLIFIER_OFF = 0,

      /// automatic choice
      SIMPLIFIER_AUTO = 1
   };

   /// values for parameter SCALER
   enum
   {
      /// no scaler
      SCALER_OFF = 0,

      /// equilibrium scaling on rows or columns
      SCALER_UNIEQUI = 1,

      /// equilibrium scaling on rows and columns
      SCALER_BIEQUI = 2,

      /// geometric mean scaling on rows and columns, max 1 round
      SCALER_GEO1 = 3,

      /// geometric mean scaling on rows and columns, max 8 rounds
      SCALER_GEO8 = 4,

       /// least square scaling
      SCALER_LEASTSQ = 5,

      /// geometric mean scaling (max 8 rounds) followed by equilibrium scaling (rows and columns)
      SCALER_GEOEQUI = 6
   };

   /// values for parameter STARTER
   enum
   {
      /// slack basis
      STARTER_OFF = 0,

      /// greedy crash basis weighted by objective, bounds, and sides
      STARTER_WEIGHT = 1,

      /// crash basis from a greedy solution
      STARTER_SUM = 2,

      /// generic solution-based crash basis
      STARTER_VECTOR = 3
   };

   /// values for parameter PRICER
   enum
   {
      /// automatic pricer
      PRICER_AUTO = 0,

      /// Dantzig pricer
      PRICER_DANTZIG = 1,

      /// partial multiple pricer based on Dantzig pricing
      PRICER_PARMULT = 2,

      /// devex pricer
      PRICER_DEVEX = 3,

      /// steepest edge pricer with initialization to unit norms
      PRICER_QUICKSTEEP = 4,

      /// steepest edge pricer with exact initialization of norms
      PRICER_STEEP = 5
   };

   /// values for parameter RATIOTESTER
   enum
   {
      /// textbook ratio test without stabilization
      RATIOTESTER_TEXTBOOK = 0,

      /// standard Harris ratio test
      RATIOTESTER_HARRIS = 1,

      /// modified Harris ratio test
      RATIOTESTER_FAST = 2,

      /// bound flipping ratio test for long steps in the dual simplex
      RATIOTESTER_BOUNDFLIPPING = 3
   };

   /// values for parameter SYNCMODE
   enum
   {
      /// store only real LP
      SYNCMODE_ONLYREAL = 0,

      /// automatic sync of real and rational LP
      SYNCMODE_AUTO = 1,

      /// user sync of real and rational LP
      SYNCMODE_MANUAL = 2
   };

   /// values for parameter READMODE
   enum
   {
      /// standard floating-point parsing
      READMODE_REAL = 0,

      /// rational parsing
      READMODE_RATIONAL = 1
   };

   /// values for parameter SOLVEMODE
   enum
   {
      /// apply standard floating-point algorithm
      SOLVEMODE_REAL = 0,

      /// decide depending on tolerances whether to apply iterative refinement
      SOLVEMODE_AUTO = 1,

      /// force iterative refinement
      SOLVEMODE_RATIONAL = 2
   };

   /// values for parameter CHECKMODE
   enum
   {
      /// floating-point check
      CHECKMODE_REAL = 0,

      /// decide according to READMODE
      CHECKMODE_AUTO = 1,

      /// rational check
      CHECKMODE_RATIONAL = 2
   };

   /// values for parameter TIMER
   enum
   {
      /// disable timing
      TIMER_OFF = 0,

      /// cpu or user time
      TIMER_CPU = 1,

      /// wallclock time
      TIMER_WALLCLOCK = 2
   };

   /// values for parameter HYPER_PRICING
   enum
   {
      /// never
      HYPER_PRICING_OFF = 0,

      /// decide according to problem size
      HYPER_PRICING_AUTO = 1,

      /// always
      HYPER_PRICING_ON = 2
   };

   /// values for parameter SOLUTION_POLISHING
   enum
   {
      /// no solution polishing
      POLISHING_OFF = 0,

      /// maximize number of basic slack variables, i.e. more variables on bounds
      POLISHING_INTEGRALITY = 1,

      /// minimize number of basic slack variables, i.e. more variables between bounds
      POLISHING_FRACTIONALITY = 2
   };

   /// real parameters
   typedef enum
   {
      /// primal feasibility tolerance
      FEASTOL = 0,

      /// dual feasibility tolerance
      OPTTOL = 1,

      /// general zero tolerance
      EPSILON_ZERO = 2,

      /// zero tolerance used in factorization
      EPSILON_FACTORIZATION = 3,

      /// zero tolerance used in update of the factorization
      EPSILON_UPDATE = 4,

      /// pivot zero tolerance used in factorization
      EPSILON_PIVOT = 5,

      /// infinity threshold
      INFTY = 6,

      /// time limit in seconds (INFTY if unlimited)
      TIMELIMIT = 7,

      /// lower limit on objective value
      OBJLIMIT_LOWER = 8,

      /// upper limit on objective value
      OBJLIMIT_UPPER = 9,

      /// working tolerance for feasibility in floating-point solver during iterative refinement
      FPFEASTOL = 10,

      /// working tolerance for optimality in floating-point solver during iterative refinement
      FPOPTTOL = 11,

      /// maximum increase of scaling factors between refinements
      MAXSCALEINCR = 12,

      /// lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
      LIFTMINVAL = 13,

      /// upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
      LIFTMAXVAL = 14,

      /// sparse pricing threshold (\#violations < dimension * SPARSITY_THRESHOLD activates sparse pricing)
      SPARSITY_THRESHOLD = 15,

      /// threshold on number of rows vs. number of columns for switching from column to row representations in auto mode
      REPRESENTATION_SWITCH = 16,

      /// geometric frequency at which to apply rational reconstruction
      RATREC_FREQ = 17,

      /// minimal reduction (sum of removed rows/cols) to continue simplification
      MINRED = 18,

      /// refactor threshold for nonzeros in last factorized basis matrix compared to updated basis matrix
      REFAC_BASIS_NNZ = 19,

      /// refactor threshold for fill-in in current factor update compared to fill-in in last factorization
      REFAC_UPDATE_FILL = 20,

      /// refactor threshold for memory growth in factorization since last refactorization
      REFAC_MEM_FACTOR = 21,

      /// accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)
      LEASTSQ_ACRCY = 22,

      /// objective offset
      OBJ_OFFSET = 23,

      /// number of real parameters
      REALPARAM_COUNT = 24
   } RealParam;

#ifdef SOPLEX_WITH_RATIONALPARAM
   /// rational parameters
   typedef enum
   {
      /// number of rational parameters
      RATIONALPARAM_COUNT = 0
   } RationalParam;
#endif

   /// class of parameter settings
   class Settings
   {
   public:
      static struct BoolParam {
         /// constructor
         BoolParam();
         /// array of names for boolean parameters
         std::string name[SoPlex::BOOLPARAM_COUNT];
         /// array of descriptions for boolean parameters
         std::string description[SoPlex::BOOLPARAM_COUNT];
         /// array of default values for boolean parameters
         bool defaultValue[SoPlex::BOOLPARAM_COUNT];
      } boolParam;

      static struct IntParam {
         /// constructor
         IntParam();
          /// array of names for integer parameters
         std::string name[SoPlex::INTPARAM_COUNT];
         /// array of descriptions for integer parameters
         std::string description[SoPlex::INTPARAM_COUNT];
         /// array of default values for integer parameters
         int defaultValue[SoPlex::INTPARAM_COUNT];
         /// array of lower bounds for int parameter values
         int lower[SoPlex::INTPARAM_COUNT];
         /// array of upper bounds for int parameter values
         int upper[SoPlex::INTPARAM_COUNT];
      } intParam;

      static struct RealParam {
         /// constructor
         RealParam();
         /// array of names for real parameters
         std::string name[SoPlex::REALPARAM_COUNT];
         /// array of descriptions for real parameters
         std::string description[SoPlex::REALPARAM_COUNT];
         /// array of default values for real parameters
         Real defaultValue[SoPlex::REALPARAM_COUNT];
         /// array of lower bounds for real parameter values
         Real lower[SoPlex::REALPARAM_COUNT];
         /// array of upper bounds for real parameter values
         Real upper[SoPlex::REALPARAM_COUNT];
      } realParam;

#ifdef SOPLEX_WITH_RATIONALPARAM
      static struct RationalParam {
         /// constructor
         RationalParam();
         /// array of names for rational parameters
         std::string name[SoPlex::RATIONALPARAM_COUNT];
         /// array of descriptions for rational parameters
         std::string description[SoPlex::RATIONALPARAM_COUNT];
         /// array of default values for rational parameters
         Rational defaultValue[SoPlex::RATIONALPARAM_COUNT];
         /// array of lower bounds for rational parameter values
         Rational lower[SoPlex::RATIONALPARAM_COUNT];
         /// array of upper bounds for rational parameter values
         Rational upper[SoPlex::RATIONALPARAM_COUNT];
      } rationalParam;
#endif

      /// array of current boolean parameter values
      bool _boolParamValues[SoPlex::BOOLPARAM_COUNT];

      /// array of current integer parameter values
      int _intParamValues[SoPlex::INTPARAM_COUNT];

      /// array of current real parameter values
      Real _realParamValues[SoPlex::REALPARAM_COUNT];

#ifdef SOPLEX_WITH_RATIONALPARAM
      /// array of current rational parameter values
      Rational _rationalParamValues[SoPlex::RATIONALPARAM_COUNT];
#endif

      /// default constructor initializing default settings
      Settings();

      /// copy constructor
      Settings(const Settings& settings);

      /// assignment operator
      Settings& operator=(const Settings& settings);
   };

   mutable SPxOut spxout;

   /// returns boolean parameter value
   bool boolParam(const BoolParam param) const;

   /// returns integer parameter value
   int intParam(const IntParam param) const;

   /// returns real parameter value
   Real realParam(const RealParam param) const;

#ifdef SOPLEX_WITH_RATIONALPARAM
   /// returns rational parameter value
   Rational rationalParam(const RationalParam param) const;
#endif

   /// returns current parameter settings
   const Settings& settings() const;

   /// sets boolean parameter value; returns true on success
   bool setBoolParam(const BoolParam param, const bool value, const bool init = true);

   /// sets integer parameter value; returns true on success
   bool setIntParam(const IntParam param, const int value, const bool init = true);

   /// sets real parameter value; returns true on success
   bool setRealParam(const RealParam param, const Real value, const bool init = true);

#ifdef SOPLEX_WITH_RATIONALPARAM
   /// sets rational parameter value; returns true on success
   bool setRationalParam(const RationalParam param, const Rational value, const bool init = true);
#endif

   /// sets parameter settings; returns true on success
   bool setSettings(const Settings& newSettings, const bool init = true);

   /// resets default parameter settings
   void resetSettings(const bool quiet = false, const bool init = true);

   /// print non-default parameter values
   void printUserSettings();

   /// writes settings file; returns true on success
   bool saveSettingsFile(const char* filename, const bool onlyChanged = false) const;

   /// reads settings file; returns true on success
   bool loadSettingsFile(const char* filename);

   /// parses one setting string and returns true on success; note that string is modified
   bool parseSettingsString(char* line);

   //@}


   //**@name Statistics */
   //@{

   /// prints solution statistics
   void printSolutionStatistics(std::ostream& os);

   /// prints statistics on solving process
   void printSolvingStatistics(std::ostream& os);

   /// prints short statistics
   void printShortStatistics(std::ostream& os);

   /// prints complete statistics
   void printStatistics(std::ostream& os);

   /// prints status
   void printStatus(std::ostream& os, SPxSolver::Status status);

   //@}


   //**@name Miscellaneous */
   //@{

   /// prints version and compilation options
   void printVersion() const;

   /// checks if real LP and rational LP are in sync; dimensions will always be compared,
   /// vector and matrix values only if the respective parameter is set to true.
   /// If quiet is set to true the function will only display which vectors are different.
   bool areLPsInSync(const bool checkVecVals = true, const bool checkMatVals = false, const bool quiet = false) const;

   /// set the random seeds of the solver instance
   void setRandomSeed(unsigned int seed);

   /// returns the current random seed of the solver instance
   unsigned int randomSeed() const;

   //@}

private:

   //**@name Statistics on solving process */
   //@{

   /// class of statistics
   class Statistics;

   /// statistics since last call to solveReal() or solveRational()
   Statistics* _statistics;

   //@}


   //**@name Parameter settings */
   //@{

   Settings* _currentSettings;

   Rational _rationalPosInfty;
   Rational _rationalNegInfty;
   Rational _rationalFeastol;
   Rational _rationalOpttol;
   Rational _rationalMaxscaleincr;

   //@}


   //**@name Data for the real LP */
   //@{

   SPxSolver _solver;
   SLUFactor _slufactor;
   SPxMainSM _simplifierMainSM;
   SPxEquiliSC _scalerUniequi;
   SPxEquiliSC _scalerBiequi;
   SPxGeometSC _scalerGeo1;
   SPxGeometSC _scalerGeo8;
   SPxGeometSC _scalerGeoequi;
   SPxLeastSqSC _scalerLeastsq;
   SPxWeightST _starterWeight;
   SPxSumST _starterSum;
   SPxVectorST _starterVector;
   SPxAutoPR _pricerAuto;
   SPxDantzigPR _pricerDantzig;
   SPxParMultPR _pricerParMult;
   SPxDevexPR _pricerDevex;
   SPxSteepPR _pricerQuickSteep;
   SPxSteepExPR _pricerSteep;
   SPxDefaultRT _ratiotesterTextbook;
   SPxHarrisRT _ratiotesterHarris;
   SPxFastRT _ratiotesterFast;
   SPxBoundFlippingRT _ratiotesterBoundFlipping;

   SPxLPReal* _realLP; // the real LP is also used as the original LP for the decomposition dual simplex
   SPxLPReal* _decompLP; // used to store the original LP for the decomposition dual simplex
   SPxSimplifier* _simplifier;
   SPxScaler* _scaler;
   SPxStarter* _starter;

   bool _isRealLPLoaded; // true indicates that the original LP is loaded in the _solver variable, hence all actions 
                         // are performed on the original LP.
   bool _isRealLPScaled;
   bool _applyPolishing;

   DVectorReal _manualLower;
   DVectorReal _manualUpper;
   DVectorReal _manualLhs;
   DVectorReal _manualRhs;
   DVectorReal _manualObj;
   SPxLPReal _manualRealLP;

   //@}


   //**@name Data for the rational LP */
   //@{

   SPxLPRational* _rationalLP;
   SLUFactorRational _rationalLUSolver;
   DataArray<int> _rationalLUSolverBind;

   LPColSetRational _slackCols;
   DVectorRational _unboundedLower;
   DVectorRational _unboundedUpper;
   DVectorRational _unboundedLhs;
   DVectorRational _unboundedRhs;
   DSVectorRational _tauColVector;
   DVectorRational _feasObj;
   DVectorRational _feasLhs;
   DVectorRational _feasRhs;
   DVectorRational _feasLower;
   DVectorRational _feasUpper;
   DVectorRational _modLower;
   DVectorRational _modUpper;
   DVectorRational _modLhs;
   DVectorRational _modRhs;
   DVectorRational _modObj;
   DSVectorRational _primalDualDiff;
   DataArray< SPxSolver::VarStatus > _storedBasisStatusRows;
   DataArray< SPxSolver::VarStatus > _storedBasisStatusCols;
   DataArray< UnitVectorRational* > _unitMatrixRational;
   bool _storedBasis;
   int _beforeLiftRows;
   int _beforeLiftCols;

   /// type of bounds and sides
   typedef enum
   {
      /// both bounds are infinite
      RANGETYPE_FREE = 0,

      /// lower bound is finite, upper bound is infinite
      RANGETYPE_LOWER = 1,

      /// upper bound is finite, lower bound is infinite
      RANGETYPE_UPPER = 2,

      /// lower and upper bound finite, but different
      RANGETYPE_BOXED = 3,

      /// lower bound equals upper bound
      RANGETYPE_FIXED = 4
   } RangeType;

   DataArray< RangeType > _colTypes;
   DataArray< RangeType > _rowTypes;

   //@}


   //**@name Data for the Decomposition Based Dual Simplex */
   //@{

   /** row violation structure
    */
   struct RowViolation
   {
      Real               violation;          /**< the violation of the row */
      int                idx;                /**< index of corresponding row */
   };

   /** Compare class for row violations
    */
   struct RowViolationCompare
   {
   public:
      /** constructor
       */
      RowViolationCompare()
         : entry(0)
      {
      }

      const RowViolation*  entry;

      Real operator() (
         RowViolation      i,
         RowViolation      j
         ) const
      {
         return i.violation - j.violation;
      }
   };


   typedef enum
   {
      // is the original problem currently being solved.
      DECOMP_ORIG = 0,

      // is the reduced problem currently being solved.
      DECOMP_RED = 1,

      // is the complementary problem currently being solved.
      DECOMP_COMP = 2
   } decompStatus;

   // the expected sign of the dual variables.
   enum DualSign
   {
      IS_POS = 0,
      IS_NEG = 1,
      IS_FREE = 2
   };

   SPxSolver _compSolver; // adding a solver to contain the complementary problem. It is too confusing to switch
                          // the LP for the reduced and complementary problem in the one solver variable. The reduced
                          // problem will be stored in _solver and the complementary problem will be stored in
                          // _compSolver.
   SLUFactor _compSlufactor; // I don't know whether this is necessary, but it is a test for now.

   SPxBasis _decompTransBasis;   // the basis required for the transformation to form the reduced problem

   DVector _transformedObj;       // the objective coefficients of the transformed problem
   DVector _decompFeasVector;       // feasibility vector calculated using unshifted bounds.
   LPRowSet _transformedRows;    // a set of the original rows that have been transformed using the original basis.
   SPxColId _compSlackColId;     // column id of the primal complementary problem slack column.
   SPxRowId _compSlackDualRowId; // row id in the dual of complementary problem related to the slack column.
   bool* _decompReducedProbRows;    // flag to indicate the inclusion of a row in the reduced problem.
   bool* _decompReducedProbCols;    // flag to indicate the inclusion of a col in the reduced problem.
   int* _decompRowStatus;
   int* _decompColStatus;
   int* _decompCompProbColIDsIdx;   // the index to _decompPrimalColIDs for a given original col.
   DataArray < SPxRowId > _decompReducedProbRowIDs;   // the row IDs for the related rows in the reduced problem
   DataArray < SPxRowId > _decompReducedProbColRowIDs;// the row IDs for the related cols in the reduced problem
   DataArray < SPxColId > _decompReducedProbColIDs;   // the col IDs for the related cols in the reduced problem
   DataArray < SPxRowId > _decompPrimalRowIDs;        // the primal row IDs from the original problem
   DataArray < SPxColId > _decompPrimalColIDs;        // the primal col IDs from the original problem
   DataArray < SPxRowId > _decompElimPrimalRowIDs;    // the primal row IDs eliminated in the complementary problem
   DataArray < SPxRowId > _decompDualRowIDs;          // the dual row IDs from the complementary problem
   DataArray < SPxColId > _decompDualColIDs;          // the dual col IDs from the complementary problem
   DataArray < SPxColId > _decompFixedVarDualIDs;     // the column ids related to the fixed variables.
   DataArray < SPxColId > _decompVarBoundDualIDs;     // the column ids related to the variable bound constraints.

   DataArray < SPxColId > _decompCompPrimalFixedVarIDs;  // the column ids related to the fixed variables in the complementary primal.
   DataArray < SPxColId > _decompCompPrimalVarBoundIDs;  // the column ids related to the variable bound constraints in the complementary primal.

   DataArray < SPxRowId > _decompCompPrimalRowIDs;        // the primal row IDs from the complementary problem
   DataArray < SPxColId > _decompCompPrimalColIDs;        // the primal col IDs from the complementary problem

   int _nDecompViolBounds;       // the number of violated bound constraints
   int _nDecompViolRows;         // the number of violated rows
   int* _decompViolatedBounds;   // the violated bounds given by the solution to the IDS reduced problem
   int* _decompViolatedRows;     // the violated rows given by the solution to the IDS reduced problem


   int* _fixedOrigVars;    // the original variables that are at their bounds in the reduced problem.
                           // 1: fixed to upper, -1: fixed to lower, 0: unfixed.
   int _nPrimalRows;       // the number of original problem rows included in the complementary problem
   int _nPrimalCols;       // the number of original problem columns included in the complementary problem
   int _nElimPrimalRows;   // the number of primal rows from the original problem eliminated from the complementary prob
   int _nDualRows;         // the number of dual rows in the complementary problem. NOTE: _nPrimalRows = _nDualCols
   int _nDualCols;         // the number of dual columns in the complementary problem. NOTE: _nPrimalRows = _nDualCols
   int _nCompPrimalRows;   // the number of rows in the complementary primal problem. NOTE: _nPrimalRows = _nCompPrimalRows
   int _nCompPrimalCols;   // the number of dual columns in the complementary problem. NOTE: _nPrimalCols = _nCompPrimalCols

   int _decompDisplayLine;     // the count for the display line

   NameSet* _rowNames;      // the row names from the input file
   NameSet* _colNames;      // the col names from the input file

   // Statistic information
   int numIncludedRows;    // the number of rows currently included in the reduced problem.
   int numDecompIter;         // the number of iterations of the decomposition dual simplex algorithm.
   int numRedProbIter;     // the number of simplex iterations performed in the reduced problem.
   int numCompProbIter;    // the number of iterations of the complementary problem.

   // problem statistics
   int numProbRows;
   int numProbCols;
   int numNonzeros;
   Real minAbsNonzero;
   Real maxAbsNonzero;

   int origCountLower;
   int origCountUpper;
   int origCountBoxed;
   int origCountFreeCol;

   int origCountLhs;
   int origCountRhs;
   int origCountRanged;
   int origCountFreeRow;


   decompStatus _currentProb;

   //@}


   //**@name Solution data */
   //@{

   SPxSolver::Status _status;
   int _lastSolveMode;

   DataArray< SPxSolver::VarStatus > _basisStatusRows;
   DataArray< SPxSolver::VarStatus > _basisStatusCols;

   SolReal _solReal;
   SolRational _solRational;
   SolRational _workSol;

   bool _hasBasis;
   bool _hasSolReal;
   bool _hasSolRational;

   //@}

   //**@name Miscellaneous */
   //@{

   int  _optimizeCalls;
   int  _unscaleCalls;

   Rational _rationalPosone;
   Rational _rationalNegone;
   Rational _rationalZero;

   //@}

   //**@name Constant helper methods */
   //@{

   /// extends sparse vector to hold newmax entries if and only if it holds no more free entries
   void _ensureDSVectorRationalMemory(DSVectorRational& vec, const int newmax) const;

   /// creates a permutation for removing rows/columns from an array of indices
   void _idxToPerm(int* idx, int idxSize, int* perm, int permSize) const;

   /// creates a permutation for removing rows/columns from a range of indices
   void _rangeToPerm(int start, int end, int* perm, int permSize) const;

   /// checks consistency
   bool _isConsistent() const;

   /// should solving process be stopped?
   bool _isSolveStopped(bool& stoppedTime, bool& stoppedIter) const;

   /// determines RangeType from real bounds
   RangeType _rangeTypeReal(const Real& lower, const Real& upper) const;

   /// determines RangeType from rational bounds
   RangeType _rangeTypeRational(const Rational& lower, const Rational& upper) const;

   /// switches RANGETYPE_LOWER to RANGETYPE_UPPER and vice versa
   RangeType _switchRangeType(const RangeType& rangeType) const;

   /// checks whether RangeType corresponds to finite lower bound
   bool _lowerFinite(const RangeType& rangeType) const;

   /// checks whether RangeType corresponds to finite upper bound
   bool _upperFinite(const RangeType& rangeType) const;

   //@}


   //**@name Non-constant helper methods */
   //@{

   /// adds a single row to the real LP and adjusts basis
   void _addRowReal(const LPRowReal& lprow);

   /// adds a single row to the real LP and adjusts basis
   void _addRowReal(Real lhs, const SVectorReal& lprow, Real rhs);

   /// adds multiple rows to the real LP and adjusts basis
   void _addRowsReal(const LPRowSetReal& lprowset);

   /// adds a single column to the real LP and adjusts basis
   void _addColReal(const LPColReal& lpcol);

   /// adds a single column to the real LP and adjusts basis
   void _addColReal(Real obj, Real lower, const SVectorReal& lpcol, Real upper);

   /// adds multiple columns to the real LP and adjusts basis
   void _addColsReal(const LPColSetReal& lpcolset);

   /// replaces row \p i with \p lprow and adjusts basis
   void _changeRowReal(int i, const LPRowReal& lprow);

   /// changes left-hand side vector for constraints to \p lhs and adjusts basis
   void _changeLhsReal(const VectorReal& lhs);

   /// changes left-hand side of row \p i to \p lhs and adjusts basis
   void _changeLhsReal(int i, const Real& lhs);

   /// changes right-hand side vector to \p rhs and adjusts basis
   void _changeRhsReal(const VectorReal& rhs);

   /// changes right-hand side of row \p i to \p rhs and adjusts basis
   void _changeRhsReal(int i, const Real& rhs);

   /// changes left- and right-hand side vectors and adjusts basis
   void _changeRangeReal(const VectorReal& lhs, const VectorReal& rhs);

   /// changes left- and right-hand side of row \p i and adjusts basis
   void _changeRangeReal(int i, const Real& lhs, const Real& rhs);

   /// replaces column \p i with \p lpcol and adjusts basis
   void _changeColReal(int i, const LPColReal& lpcol);

   /// changes vector of lower bounds to \p lower and adjusts basis
   void _changeLowerReal(const VectorReal& lower);

   /// changes lower bound of column i to \p lower and adjusts basis
   void _changeLowerReal(int i, const Real& lower);

   /// changes vector of upper bounds to \p upper and adjusts basis
   void _changeUpperReal(const VectorReal& upper);

   /// changes \p i 'th upper bound to \p upper and adjusts basis
   void _changeUpperReal(int i, const Real& upper);

   /// changes vectors of column bounds to \p lower and \p upper and adjusts basis
   void _changeBoundsReal(const VectorReal& lower, const VectorReal& upper);

   /// changes bounds of column \p i to \p lower and \p upper and adjusts basis
   void _changeBoundsReal(int i, const Real& lower, const Real& upper);

   /// changes matrix entry in row \p i and column \p j to \p val and adjusts basis
   void _changeElementReal(int i, int j, const Real& val);

   /// removes row \p i and adjusts basis
   void _removeRowReal(int i);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void _removeRowsReal(int perm[]);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsReal() may be passed
   /// as buffer memory
   void _removeRowsReal(int idx[], int n, int perm[]);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsReal() may be passed as buffer
   /// memory
   void _removeRowRangeReal(int start, int end, int perm[]);

   /// removes column i
   void _removeColReal(int i);

   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void _removeColsReal(int perm[]);

   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void _removeColsReal(int idx[], int n, int perm[]);

   /// removes columns \p start to \p end including both; an array \p perm of size #numColsReal() may be passed as
   /// buffer memory
   void _removeColRangeReal(int start, int end, int perm[]);

   /// invalidates solution
   void _invalidateSolution();

   /// enables simplifier and scaler according to current parameters
   void _enableSimplifierAndScaler();

   /// disables simplifier and scaler
   void _disableSimplifierAndScaler();

   /// ensures that the rational LP is available; performs no sync
   void _ensureRationalLP();

   /// ensures that the real LP and the basis are loaded in the solver; performs no sync
   void _ensureRealLPLoaded();

   /// call floating-point solver and update statistics on iterations etc.
   void _solveRealLPAndRecordStatistics();

   /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool _readFileReal(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// reads rational LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool _readFileRational(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// completes range type arrays after adding columns and/or rows
   void _completeRangeTypesRational();

   /// recomputes range types from scratch using real LP
   void _recomputeRangeTypesReal();

   /// recomputes range types from scratch using rational LP
   void _recomputeRangeTypesRational();

   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, without looking at the sync mode
   void _syncLPReal(bool time = true);

   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, without looking at the sync mode
   void _syncLPRational(bool time = true);

   /// synchronizes rational solution with real solution, i.e., copies (rounded) rational solution to real solution
   void _syncRealSolution();

   /// synchronizes real solution with rational solution, i.e., copies real solution to rational solution
   void _syncRationalSolution();

   /// returns pointer to a constant unit vector available until destruction of the SoPlex class
   const UnitVectorRational* _unitVectorRational(const int i);

   /// parses one line in a settings file and returns true on success; note that the string is modified
   bool _parseSettingsLine(char* line, const int lineNumber);

   //@}


   //**@name Private solving methods implemented in solverational.cpp */
   //@{

   /// solves rational LP
   void _optimizeRational();

   /// solves current problem with iterative refinement and recovery mechanism
   void _performOptIRStable(SolRational& sol,
      bool acceptUnbounded,
      bool acceptInfeasible,
      int minRounds,
      bool& primalFeasible,
      bool& dualFeasible,
      bool& infeasible,
      bool& unbounded,
      bool& stopped,
      bool& stoppedIter,
      bool& error);

   /// performs iterative refinement on the auxiliary problem for testing unboundedness
   void _performUnboundedIRStable(SolRational& sol, bool& hasUnboundedRay, bool& stopped, bool& stoppedIter, bool& error);

   /// performs iterative refinement on the auxiliary problem for testing feasibility
   void _performFeasIRStable(SolRational& sol, bool& withDualFarkas, bool& stopped, bool& stoppedIter, bool& error);

   /// reduces matrix coefficient in absolute value by the lifting procedure of Thiele et al. 2013
   void _lift();

   /// undoes lifting
   void _project(SolRational& sol);

   /// store basis
   void _storeBasis();

   /// restore basis
   void _restoreBasis();

   /// stores objective, bounds, and sides of real LP
   void _storeLPReal();

   /// restores objective, bounds, and sides of real LP
   void _restoreLPReal();

   /// introduces slack variables to transform inequality constraints into equations for both rational and real LP,
   /// which should be in sync
   void _transformEquality();

   /// undoes transformation to equality form
   void _untransformEquality(SolRational& sol);

   /// transforms LP to unboundedness problem by moving the objective function to the constraints, changing right-hand
   /// side and bounds to zero, and adding an auxiliary variable for the decrease in the objective function
   void _transformUnbounded();

   /// undoes transformation to unboundedness problem
   void _untransformUnbounded(SolRational& sol, bool unbounded);

   /// transforms LP to feasibility problem by removing the objective function, shifting variables, and homogenizing the
   /// right-hand side
   void _transformFeasibility();

   /// undoes transformation to feasibility problem
   void _untransformFeasibility(SolRational& sol, bool infeasible);

   /** computes radius of infeasibility box implied by an approximate Farkas' proof

     Given constraints of the form \f$ lhs <= Ax <= rhs \f$, a farkas proof y should satisfy \f$ y^T A = 0 \f$ and
     \f$ y_+^T lhs - y_-^T rhs > 0 \f$, where \f$ y_+, y_- \f$ denote the positive and negative parts of \f$ y \f$.
     If \f$ y \f$ is approximate, it may not satisfy \f$ y^T A = 0 \f$ exactly, but the proof is still valid as long
     as the following holds for all potentially feasible \f$ x \f$:

     \f[
         y^T Ax < (y_+^T lhs - y_-^T rhs)              (*)
     \f]

     we may therefore calculate \f$ y^T A \f$ and \f$ y_+^T lhs - y_-^T rhs \f$ exactly and check if the upper and lower
     bounds on \f$ x \f$ imply that all feasible \f$ x \f$ satisfy (*), and if not then compute bounds on \f$ x \f$ to
     guarantee (*).  The simplest way to do this is to compute

     \f[
     B = (y_+^T lhs - y_-^T rhs) / \sum_i(|(y^T A)_i|)
     \f]

     noting that if every component of \f$ x \f$ has \f$ |x_i| < B \f$, then (*) holds.

     \f$ B \f$ can be increased by iteratively including variable bounds smaller than \f$ B \f$.  The speed of this
     method can be further improved by using interval arithmetic for all computations.  For related information see
     Sec. 4 of Neumaier and Shcherbina, Mathematical Programming A, 2004.

     Set transformed to true if this method is called after _transformFeasibility().
     */
   void _computeInfeasBox(SolRational& sol, bool transformed);

   /// solves real LP during iterative refinement
   SPxSolver::Status _solveRealForRational(bool fromscratch, VectorReal& primal, VectorReal& dual,
                                           DataArray< SPxSolver::VarStatus >& basisStatusRows,
                                           DataArray< SPxSolver::VarStatus >& basisStatusCols, bool& returnedBasis);

   /// solves real LP with recovery mechanism
   SPxSolver::Status _solveRealStable(bool acceptUnbounded, bool acceptInfeasible, VectorReal& primal, VectorReal& dual,
                                      DataArray< SPxSolver::VarStatus >& basisStatusRows,
                                      DataArray< SPxSolver::VarStatus >& basisStatusCols, bool& returnedBasis, const bool forceNoSimplifier = false);

   /// computes rational inverse of basis matrix as defined by _rationalLUSolverBind
   void _computeBasisInverseRational();

   /// factorizes rational basis matrix in column representation
   void _factorizeColumnRational(SolRational& sol, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols, bool& stoppedTime, bool& stoppedIter, bool& error, bool& optimal);

   /// attempts rational reconstruction of primal-dual solution
   bool _reconstructSolutionRational(SolRational& sol, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols, const Rational& denomBoundSquared);
   //@}


   //**@name Private solving methods implemented in solvereal.cpp */
   //@{

   /// solves real LP
   void _optimizeReal();

   /// checks result of the solving process and solves again without preprocessing if necessary
   void _evaluateSolutionReal(SPxSimplifier::Result simplificationStatus);

   /// solves real LP with/without preprocessing
   void _preprocessAndSolveReal(bool applyPreprocessing);

   /// loads original problem into solver and solves again after it has been solved to optimality with preprocessing
   void _resolveWithoutPreprocessing(SPxSimplifier::Result simplificationStatus);

   /// verify computed solution and resolve if necessary
   void _verifySolutionReal();

   /// stores solution of the real LP; before calling this, the real LP must be loaded in the solver and solved (again)
   void _storeSolutionReal(bool verify = true);

   /// stores solution from the simplifier because problem vanished in presolving step
   void _storeSolutionRealFromPresol();

   /// unscales stored solution to remove internal or external scaling of LP
   void _unscaleSolutionReal(SPxLPReal& LP, bool persistent = true);

   /// load original LP and possibly setup a slack basis
   void _loadRealLP(bool initBasis);

   /// check scaling of LP
   void _checkScaling(SPxLPReal* origLP) const;

   /// check correctness of (un)scaled basis matrix operations
   void _checkBasisScaling();

   /// check whether persistent scaling is supposed to be reapplied again after unscaling
   bool _reapplyPersistentScaling() const;

   /// solves LP using the decomposition based dual simplex
   void _solveDecompositionDualSimplex();

   /// creating copies of the original problem that will be manipulated to form the reduced and complementary problems
   void _createDecompReducedAndComplementaryProblems();

   /// forms the reduced problem
   void _formDecompReducedProblem(bool& stop);

   /// solves the reduced problem
   void _solveDecompReducedProblem();

   /// forms the complementary problem
   void _formDecompComplementaryProblem();

   /// simplifies the problem and solves
   void _decompSimplifyAndSolve(SPxSolver& solver, SLUFactor& sluFactor, bool fromScratch, bool applyPreprocessing);

   /// loads original problem into solver and solves again after it has been solved to optimality with preprocessing
   void _decompResolveWithoutPreprocessing(SPxSolver& solver, SLUFactor& sluFactor, SPxSimplifier::Result result);

   /// identifies the columns of the row-form basis that correspond to rows with zero dual multipliers.
   void _getZeroDualMultiplierIndices(Vector feasVector, int* nonposind, int* colsforremoval,
         int* nnonposind, bool& stop);

   /// retrieves the compatible columns from the constraint matrix
   void _getCompatibleColumns(Vector feasVector, int* nonposind, int* compatind, int* rowsforremoval, int* colsforremoval,
      int nnonposind, int* ncompatind, bool formRedProb, bool& stop);

   /// computes the reduced problem objective coefficients
   void _computeReducedProbObjCoeff(bool& stop);

   /// computes the compatible bound constraints and adds them to the reduced problem
   void _getCompatibleBoundCons(LPRowSet& boundcons, int* compatboundcons, int* nonposind, int* ncompatboundcons,
         int nnonposind, bool& stop);

   /// computes the rows to remove from the complementary problem
   void _getRowsForRemovalComplementaryProblem(int* nonposind, int* bind, int* rowsforremoval, int* nrowsforremoval,
         int nnonposind);

   /// removing rows from the complementary problem.
   void _deleteAndUpdateRowsComplementaryProblem(SPxRowId rangedRowIds[], int& naddedrows);

   /// evaluates the solution of the reduced problem for the DBDS
   void _evaluateSolutionDecomp(SPxSolver& solver, SLUFactor& sluFactor, SPxSimplifier::Result result);

   /// update the reduced problem with additional columns and rows
   void _updateDecompReducedProblem(Real objVal, DVector dualVector, DVector redcostVector, DVector compPrimalVector,
      DVector compDualVector);

   /// update the reduced problem with additional columns and rows based upon the violated original bounds and rows
   void _updateDecompReducedProblemViol(bool allrows);

   /// builds the update rows with those violated in the complmentary problem
   void _findViolatedRows(Real compObjValue, DataArray<RowViolation>& violatedrows, int& nviolatedrows);

   /// update the dual complementary problem with additional columns and rows
   void _updateDecompComplementaryDualProblem(bool origObj);

   /// update the primal complementary problem with additional columns and rows
   void _updateDecompComplementaryPrimalProblem(bool origObj);

   /// checking the optimality of the original problem.
   void _checkOriginalProblemOptimality(Vector primalVector, bool printViol);

   /// updating the slack column coefficients to adjust for equality constraints
   void _updateComplementaryDualSlackColCoeff();

   /// updating the slack column coefficients to adjust for equality constraints
   void _updateComplementaryPrimalSlackColCoeff();

   /// removing the dual columns related to the fixed variables
   void _removeComplementaryDualFixedPrimalVars(int* currFixedVars);

   /// removing the dual columns related to the fixed variables
   void _identifyComplementaryDualFixedPrimalVars(int* currFixedVars);

   /// updating the dual columns related to the fixed primal variables.
   void _updateComplementaryDualFixedPrimalVars(int* currFixedVars);

   /// removing the dual columns related to the fixed variables
   void _identifyComplementaryPrimalFixedPrimalVars(int* currFixedVars);

   /// updating the dual columns related to the fixed primal variables.
   void _updateComplementaryPrimalFixedPrimalVars(int* currFixedVars);

   /// updating the complementary dual problem with the original objective function
   void _setComplementaryDualOriginalObjective();

   /// updating the complementary primal problem with the original objective function
   void _setComplementaryPrimalOriginalObjective();

   /// determining which bound the primal variables will be fixed to.
   int getOrigVarFixedDirection(int colNum);

   /// checks the dual feasibility of the current basis
   bool checkBasisDualFeasibility(Vector feasVec);

   /// returns the expected sign of the dual variables for the reduced problem
   DualSign getExpectedDualVariableSign(int rowNumber);

   /// returns the expected sign of the dual variables for the original problem
   DualSign getOrigProbDualVariableSign(int rowNumber);

   /// prints a display line of the flying table for the DBDS 
   void printDecompDisplayLine(SPxSolver& solver, const SPxOut::Verbosity origVerb, bool force, bool forceHead);

   /// stores the problem statistics of the original problem
   void getOriginalProblemStatistics();

   /// stores the problem statistics of the original problem
   void printOriginalProblemStatistics(std::ostream& os);

   /// gets the coefficient of the slack variable in the primal complementary problem
   Real getCompSlackVarCoeff(int primalRowNum);

   /// gets violation of bounds; returns true on success
   bool getDecompBoundViolation(Real& maxviol, Real& sumviol);

   /// gets violation of constraints; returns true on success
   bool getDecompRowViolation(Real& maxviol, Real& sumviol);

   /// function call to terminate the decomposition simplex
   bool decompTerminate(Real timeLimit);

   /// function to build a basis for the original problem as given by the solution to the reduced problem
   void _writeOriginalProblemBasis(const char* filename, NameSet* rowNames, NameSet* colNames, bool cpxFormat);

   /// function to retrieve the original problem row basis status from the reduced and complementary problems
   void getOriginalProblemBasisRowStatus(DataArray< int >& degenerateRowNums,
      DataArray< SPxSolver::VarStatus >& degenerateRowStatus, int& nDegenerateRows, int& nNonBasicRows);

   /// function to retrieve the column status for the original problem basis from the reduced and complementary problems
   void getOriginalProblemBasisColStatus(int& nNonBasicCols);

   //@}
};
}
#else
#include "soplexlegacy.h"

namespace soplex
{
   typedef SoPlexLegacy SoPlex;
}
#endif
#endif // _SOPLEX_H_
