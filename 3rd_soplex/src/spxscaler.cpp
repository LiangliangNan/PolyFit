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

/**@file  spxscaler.cpp
 * @brief LP scaling base class.
 */

#include <cmath>

#include <iostream>
#include <assert.h>
#include "spxscaler.h"
#include "spxlp.h"
#include "dsvector.h"
#include "dvector.h"
#include <limits>

namespace soplex
{

std::ostream& operator<<(std::ostream& s, const SPxScaler& sc)
{
   const DataArray < int >& colscaleExp = *(sc.m_activeColscaleExp);
   DataArray < int > rowccaleExp = *(sc.m_activeRowscaleExp);

   s << sc.getName() << " scaler:" << std::endl;
   s << "colscale = [ ";
   for(int ci = 0; ci < colscaleExp.size(); ++ci )
      s << colscaleExp[ci] << " ";
   s << "]" << std::endl;

   s << "rowscale = [ ";
   for(int ri = 0; ri < rowccaleExp.size(); ++ri )
      s << rowccaleExp[ri] << " ";
   s << "]" << std::endl;

   return s;
}


SPxScaler::SPxScaler(
   const char* name, 
   bool        colFirst, 
   bool        doBoth,
   SPxOut*     outstream)
   : m_name(name)
   , m_activeColscaleExp(0)
   , m_activeRowscaleExp(0)
   , m_colFirst(colFirst)
   , m_doBoth(doBoth)
   , spxout(outstream)
{
   assert(SPxScaler::isConsistent());
}

SPxScaler::SPxScaler(const SPxScaler& old)
   : m_name(old.m_name)
   , m_activeColscaleExp(old.m_activeColscaleExp)
   , m_activeRowscaleExp(old.m_activeRowscaleExp)
   , m_colFirst(old.m_colFirst)
   , m_doBoth(old.m_doBoth)
   , spxout(old.spxout)
{
   assert(SPxScaler::isConsistent());
}

SPxScaler::~SPxScaler()
{
   m_name = 0;
}

SPxScaler& SPxScaler::operator=(const SPxScaler& rhs)
{
   if (this != &rhs)
   {
      m_name     = rhs.m_name;
      m_activeColscaleExp = rhs.m_activeColscaleExp;
      m_activeRowscaleExp = rhs.m_activeRowscaleExp;
      m_colFirst = rhs.m_colFirst;
      m_doBoth   = rhs.m_doBoth;
      spxout     = rhs.spxout;

      assert(SPxScaler::isConsistent());
   }
   return *this;
}

const char* SPxScaler::getName() const
{

   return m_name;
}

void SPxScaler::setOrder(bool colFirst)
{

   m_colFirst = colFirst;
}

void SPxScaler::setBoth(bool both)
{

   m_doBoth = both;
}

void SPxScaler::setRealParam(Real param, const char* name)
{}

void SPxScaler::setIntParam(int param, const char* name)
{}

void SPxScaler::setup(SPxLP& lp)
{
   assert(lp.isConsistent());
   m_activeColscaleExp = &lp.LPColSetBase<Real>::scaleExp;
   m_activeRowscaleExp = &lp.LPRowSetBase<Real>::scaleExp;
   m_activeColscaleExp->reSize(lp.nCols());
   m_activeRowscaleExp->reSize(lp.nRows());

   for( int i = 0; i < lp.nCols(); ++i)
      (*m_activeColscaleExp)[i] = 0;
   for( int i = 0; i < lp.nRows(); ++i)
      (*m_activeRowscaleExp)[i] = 0;

   lp.lp_scaler = this;
}

int SPxScaler::computeScaleExp(const SVector& vec, const DataArray<int>& oldScaleExp) const
{
   Real maxi = 0.0;

   // find largest absolute value after applying existing scaling factors
   for( int i = 0; i < vec.size(); ++i )
   {
      Real x = spxAbs(spxLdexp(vec.value(i), oldScaleExp[vec.index(i)]));

      if( GT(x, maxi) )
         maxi = x;
   }
   // empty rows/cols are possible
   if( maxi == 0.0 )
      return 0;
   // get exponent corresponding to new scaling factor
   else
   {
      int scaleExp;
      spxFrexp(1.0 / maxi, &(scaleExp));
      return scaleExp - 1;
   }
}

#ifndef SOPLEX_LEGACY
int SPxScaler::computeScaleExp(const SVectorBase<Rational>& vec, const DataArray<int>& oldScaleExp) const
{
   return 0;
}
#endif

void SPxScaler::applyScaling(SPxLPBase<Real>& lp)
{
   assert(lp.nCols() == m_activeColscaleExp->size());
   assert(lp.nRows() == m_activeRowscaleExp->size());

   DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;
   DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   for( int i = 0; i < lp.nRows(); ++i )
   {
      SVector& vec = lp.rowVector_w(i);
      int exp1;
      int exp2 = rowscaleExp[i];

      for( int j = 0; j < vec.size(); ++j)
      {
         exp1 = colscaleExp[vec.index(j)];
         vec.value(j) = spxLdexp(vec.value(j), exp1 + exp2);
      }

      lp.maxRowObj_w(i) = spxLdexp(lp.maxRowObj(i), exp2);

      if( lp.rhs(i) < infinity )
         lp.rhs_w(i) = spxLdexp(lp.rhs_w(i), exp2);

      if( lp.lhs(i) > -infinity )
         lp.lhs_w(i) = spxLdexp(lp.lhs_w(i), exp2);

      MSG_DEBUG( std::cout << "DEBUG: rowscaleExp(" << i << "): " << exp2 << std::endl; )
   }

   for( int i = 0; i < lp.nCols(); ++i )
   {
      SVector& vec = lp.colVector_w(i);
      int exp1;
      int exp2 = colscaleExp[i];

      for( int j = 0; j < vec.size(); ++j)
      {
         exp1 = rowscaleExp[vec.index(j)];
         vec.value(j) = spxLdexp(vec.value(j), exp1 + exp2);
      }

      lp.maxObj_w(i) = spxLdexp(lp.maxObj_w(i), exp2);

      if( lp.upper(i) < infinity )
         lp.upper_w(i) = spxLdexp(lp.upper_w(i), -exp2);

      if( lp.lower(i) > -infinity )
         lp.lower_w(i) = spxLdexp(lp.lower_w(i), -exp2);

      MSG_DEBUG( std::cout << "DEBUG: colscaleExp(" << i << "): " << exp2 << std::endl; )
   }

   lp.setScalingInfo(true);
   assert(lp.isConsistent());
}

/// unscale SPxLP
void SPxScaler::unscale(SPxLPBase<Real>& lp)
{
   assert(lp.isScaled());

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;
   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   for( int i = 0; i < lp.nRows(); ++i )
   {
      SVector& vec = lp.rowVector_w(i);

      int exp1;
      int exp2 = rowscaleExp[i];

      for( int j = 0; j < vec.size(); ++j)
      {
         exp1 = colscaleExp[vec.index(j)];
         vec.value(j) = spxLdexp(vec.value(j), -exp1 - exp2);
      }

      lp.maxRowObj_w(i) = spxLdexp(lp.maxRowObj(i), -exp2);

      if( lp.rhs(i) < infinity )
         lp.rhs_w(i) = spxLdexp(lp.rhs_w(i), -exp2);

      if( lp.lhs(i) > -infinity )
         lp.lhs_w(i) = spxLdexp(lp.lhs_w(i), -exp2);
   }

   for( int i = 0; i < lp.nCols(); ++i )
   {
      SVector& vec = lp.colVector_w(i);

      int exp1;
      int exp2 = colscaleExp[i];

      for( int j = 0; j < vec.size(); ++j)
      {
         exp1 = rowscaleExp[vec.index(j)];
         vec.value(j) = spxLdexp(vec.value(j), -exp1 - exp2);
      }

      lp.maxObj_w(i) = spxLdexp(lp.maxObj_w(i), -exp2);

      if( lp.upper(i) < infinity )
         lp.upper_w(i) = spxLdexp(lp.upper_w(i), exp2);

      if( lp.lower(i) > -infinity )
         lp.lower_w(i) = spxLdexp(lp.lower_w(i), exp2);
   }

   lp._isScaled = false;
   assert(lp.isConsistent());
}

/// returns scaling factor for column \p i
/// todo pass the LP?!
int SPxScaler::getColScaleExp(int i) const
{
   return (*m_activeColscaleExp)[i];
}

/// returns scaling factor for row \p i
/// todo pass the LP?!
int SPxScaler::getRowScaleExp(int i) const
{
   return (*m_activeRowscaleExp)[i];
}

/// gets unscaled column \p i
void SPxScaler::getColUnscaled(const SPxLP& lp, int i, DSVector& vec) const
{
   assert(lp.isScaled());
   assert(i < lp.nCols());
   assert(i >= 0);
   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;
   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   vec = lp.LPColSet::colVector(i);

   int exp1;
   int exp2 = colscaleExp[i];

   const SVectorReal& col = lp.colVector(i);
   vec.setMax(col.size());
   vec.clear();

   for( int j = 0; j < col.size(); j++ )
   {
      exp1 = rowscaleExp[col.index(j)];
      vec.add(col.index(j), spxLdexp(col.value(j), -exp1 - exp2));
   }
}

/// returns maximum absolute value of unscaled column \p i
Real SPxScaler::getColMaxAbsUnscaled(const SPxLP& lp, int i) const
{
   assert(i < lp.nCols());
   assert(i >= 0);

   DataArray < int >& colscaleExp = *m_activeColscaleExp;
   DataArray < int >& rowscaleExp = *m_activeRowscaleExp;
   const SVector& colVec = lp.LPColSet::colVector(i);

   Real max = 0.0;
   int exp1;
   int exp2 = colscaleExp[i];

   for( int j = 0; j < colVec.size(); j++ )
   {
      exp1 = rowscaleExp[colVec.index(j)];
      Real abs = spxAbs(spxLdexp(colVec.value(j), -exp1 - exp2));
      if( abs > max )
         max = abs;
   }

   return max;
}

/// returns minimum absolute value of unscaled column \p i
Real SPxScaler::getColMinAbsUnscaled(const SPxLP& lp, int i) const
{
   assert(i < lp.nCols());
   assert(i >= 0);

   DataArray < int >& colscaleExp = *m_activeColscaleExp;
   DataArray < int >& rowscaleExp = *m_activeRowscaleExp;
   const SVector& colVec = lp.LPColSet::colVector(i);

   Real min = infinity;
   int exp1;
   int exp2 = colscaleExp[i];

   for( int j = 0; j < colVec.size(); j++ )
   {
      exp1 = rowscaleExp[colVec.index(j)];
      Real abs = spxAbs(spxLdexp(colVec.value(j), -exp1 - exp2));
      if( abs < min )
         min = abs;
   }

   return min;
}


/// returns unscaled upper bound \p i
Real SPxScaler::upperUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(lp.isScaled());
   assert(i < lp.nCols());
   assert(i >= 0);

   if( lp.LPColSet::upper(i) < infinity )
   {
      const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;
      return spxLdexp(lp.LPColSet::upper(i) , colscaleExp[i]);
   }
   else
      return lp.LPColSet::upper(i);
}


/// gets unscaled upper bound vector
void SPxScaler::getUpperUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   assert(lp.isScaled());
   assert(lp.LPColSet::upper().dim() == vec.dim());

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   for( int i = 0; i < lp.LPColSet::upper().dim(); i++)
      vec[i] = spxLdexp(lp.LPColSet::upper()[i], colscaleExp[i]);
}


/// returns unscaled upper bound vector of \p lp
Real SPxScaler::lowerUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(lp.isScaled());
   assert(i < lp.nCols());
   assert(i >= 0);

   if( lp.LPColSet::lower(i) > -infinity )
   {
      const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;
      return spxLdexp(lp.LPColSet::lower(i), colscaleExp[i]);
   }
   else
      return lp.LPColSet::lower(i);
}


/// returns unscaled lower bound vector of \p lp
void SPxScaler::getLowerUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   assert(lp.isScaled());
   assert(lp.LPColSet::lower().dim() == vec.dim());

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   for( int i = 0; i < lp.LPColSet::lower().dim(); i++)
      vec[i] = spxLdexp(lp.LPColSet::lower()[i], colscaleExp[i]);
}

/// returns unscaled objective function coefficient of \p i
Real SPxScaler::maxObjUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(lp.isScaled());
   assert(i < lp.nCols());
   assert(i >= 0);

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   return spxLdexp(lp.LPColSet::maxObj(i) , -colscaleExp[i]);
}


/// gets unscaled objective function coefficient of \p i
void SPxScaler::getMaxObjUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   assert(lp.isScaled());
   assert(lp.LPColSet::maxObj().dim() == vec.dim());

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   for( int i = 0; i < lp.LPColSet::maxObj().dim(); i++)
      vec[i] = spxLdexp(lp.LPColSet::maxObj()[i], -colscaleExp[i]);
}

/// gets unscaled row \p i
void SPxScaler::getRowUnscaled(const SPxLP& lp, int i, DSVector& vec) const
{
   assert(lp.isScaled());
   assert(i < lp.nRows());
   assert(i >= 0);

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;
   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;
   int exp1;
   int exp2 = rowscaleExp[i];

   const SVectorReal& row = lp.rowVector(i);
   vec.setMax(row.size());
   vec.clear();

   for( int j = 0; j < row.size(); j++ )
   {
      exp1 = colscaleExp[row.index(j)];
      vec.add(row.index(j), spxLdexp(row.value(j), -exp1 - exp2));
   }
}

/// returns maximum absolute value of unscaled row \p i
Real SPxScaler::getRowMaxAbsUnscaled(const SPxLP& lp, int i) const
{
   assert(i < lp.nRows());
   assert(i >= 0);
   DataArray < int >& colscaleExp = *m_activeColscaleExp;
   DataArray < int >& rowscaleExp = *m_activeRowscaleExp;
   const SVector& rowVec = lp.LPRowSet::rowVector(i);

   Real max = 0.0;

   int exp1;
   int exp2 = rowscaleExp[i];

   for( int j = 0; j < rowVec.size(); j++ )
   {
      exp1 = colscaleExp[rowVec.index(j)];
      Real abs = spxAbs(spxLdexp(rowVec.value(j), -exp1 - exp2));

      if( GT(abs, max) )
         max = abs;
   }

   return max;
}

/// returns minimum absolute value of unscaled row \p i
Real SPxScaler::getRowMinAbsUnscaled(const SPxLP& lp, int i) const
{
   assert(i < lp.nRows());
   assert(i >= 0);
   DataArray < int >& colscaleExp = *m_activeColscaleExp;
   DataArray < int >& rowscaleExp = *m_activeRowscaleExp;
   const SVector& rowVec = lp.LPRowSet::rowVector(i);

   Real min = infinity;

   int exp1;
   int exp2 = rowscaleExp[i];

   for( int j = 0; j < rowVec.size(); j++ )
   {
      exp1 = colscaleExp[rowVec.index(j)];
      Real abs = spxAbs(spxLdexp(rowVec.value(j), -exp1 - exp2));

      if( LT(abs, min) )
         min = abs;
   }

   return min;
}

/// returns unscaled right hand side \p i
Real SPxScaler::rhsUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(lp.isScaled());
   assert(i < lp.nRows());
   assert(i >= 0);

   if( lp.LPRowSet::rhs(i) < infinity )
   {
      const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;
      return spxLdexp(lp.LPRowSet::rhs(i) , -rowscaleExp[i]);
   }
   else
      return lp.LPRowSet::rhs(i);
}


/// gets unscaled right hand side vector
void SPxScaler::getRhsUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   assert(lp.isScaled());
   assert(lp.LPRowSet::rhs().dim() == vec.dim());

   for( int i = 0; i < lp.LPRowSet::rhs().dim(); i++)
   {
      const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;
      vec[i] = spxLdexp(lp.LPRowSet::rhs()[i], -rowscaleExp[i]);
   }
}


/// returns unscaled left hand side \p i of \p lp
Real SPxScaler::lhsUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(lp.isScaled());
   assert(i < lp.nRows());
   assert(i >= 0);

   if( lp.LPRowSet::lhs(i) > -infinity )
   {
      const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;
      return spxLdexp(lp.LPRowSet::lhs(i) , -rowscaleExp[i]);
   }
   else
      return lp.LPRowSet::lhs(i);
}

/// returns unscaled left hand side vector of \p lp
void SPxScaler::getLhsUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   assert(lp.isScaled());
   assert(lp.LPRowSet::lhs().dim() == vec.dim());

   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   for( int i = 0; i < lp.LPRowSet::lhs().dim(); i++)
      vec[i] = spxLdexp(lp.LPRowSet::lhs()[i], -rowscaleExp[i]);
}

/// returns unscaled coefficient of \p lp
Real SPxScaler::getCoefUnscaled(const SPxLPBase<Real>& lp, int row, int col) const
{
   assert(lp.isScaled());
   assert(row < lp.nRows());
   assert(col < lp.nCols());

   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;
   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   return spxLdexp(lp.colVector(col)[row], - rowscaleExp[row] - colscaleExp[col]);
}

void SPxScaler::unscalePrimal(const SPxLPBase<Real>& lp, Vector& x) const
{
   assert(lp.isScaled());

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   assert(x.dim() == colscaleExp.size());

   for( int j = 0; j < x.dim(); ++j )
      x[j] = spxLdexp(x[j], colscaleExp[j]);
}

void SPxScaler::unscaleSlacks(const SPxLPBase<Real>& lp, Vector& s) const
{
   assert(lp.isScaled());

   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   assert(s.dim() == rowscaleExp.size());

   for( int i = 0; i < s.dim(); ++i )
      s[i] = spxLdexp(s[i], -rowscaleExp[i]);
}

void SPxScaler::unscaleDual(const SPxLPBase<Real>& lp, Vector& pi) const
{
   assert(lp.isScaled());

   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   assert(pi.dim() == rowscaleExp.size());

   for( int i = 0; i < pi.dim(); ++i )
      pi[i] = spxLdexp(pi[i], rowscaleExp[i]);
}

void SPxScaler::unscaleRedCost(const SPxLPBase<Real>& lp, Vector& r) const
{
   assert(lp.isScaled());

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   assert(r.dim() == colscaleExp.size());

   for( int j = 0; j < r.dim(); ++j )
      r[j] = spxLdexp(r[j], -colscaleExp[j]);
}

void SPxScaler::unscalePrimalray(const SPxLPBase<Real>& lp, Vector& ray) const
{
   assert(lp.isScaled());

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   assert(ray.dim() == colscaleExp.size());

   for( int j = 0; j < ray.dim(); ++j )
      ray[j] = spxLdexp(ray[j], colscaleExp[j]);
}

void SPxScaler::unscaleDualray(const SPxLPBase<Real>& lp, Vector& ray) const
{
   assert(lp.isScaled());

   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   assert(ray.dim() == rowscaleExp.size());

   for( int i = 0; i < ray.dim(); ++i )
      ray[i] = spxLdexp(ray[i], rowscaleExp[i]);
}

void SPxScaler::scaleObj(const SPxLPBase<Real>& lp, VectorReal& origObj) const
{
   assert(lp.isScaled());

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   for( int i = 0; i < origObj.dim(); ++i )
   {
      origObj[i] = spxLdexp(origObj[i], colscaleExp[i]);
   }
}

Real SPxScaler::scaleObj(const SPxLPBase<Real>& lp, int i, Real origObj) const
{
   assert(lp.isScaled());
   assert(i < lp.nCols());
   assert(i >= 0);

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;
   int exp = colscaleExp[i];

   return spxLdexp(origObj, exp);
}

Real SPxScaler::scaleElement(const SPxLPBase<Real>& lp, int row, int col, Real val) const
{
   assert(lp.isScaled());
   assert(col < lp.nCols());
   assert(col >= 0);
   assert(row < lp.nRows());
   assert(row >= 0);

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;
   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   return spxLdexp(val, colscaleExp[col] + rowscaleExp[row]);
}

Real SPxScaler::scaleLower(const SPxLPBase<Real>& lp, int col, Real lower) const
{
   assert(lp.isScaled());
   assert(col < lp.nCols());
   assert(col >= 0);

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   return spxLdexp(lower, -colscaleExp[col]);
}

Real SPxScaler::scaleUpper(const SPxLPBase<Real>& lp, int col, Real upper) const
{
   assert(lp.isScaled());
   assert(col < lp.nCols());
   assert(col >= 0);

   const DataArray < int >& colscaleExp = lp.LPColSetBase<Real>::scaleExp;

   return spxLdexp(upper, -colscaleExp[col]);
}

Real SPxScaler::scaleLhs(const SPxLPBase<Real>& lp, int row, Real lhs) const
{
   assert(lp.isScaled());
   assert(row < lp.nRows());
   assert(row >= 0);

   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   return spxLdexp(lhs, rowscaleExp[row]);
}

Real SPxScaler::scaleRhs(const SPxLPBase<Real>& lp, int row, Real rhs) const
{
   assert(lp.isScaled());
   assert(row < lp.nRows());
   assert(row >= 0);

   const DataArray < int >& rowscaleExp = lp.LPRowSetBase<Real>::scaleExp;

   return spxLdexp(rhs, rowscaleExp[row]);
}

Real SPxScaler::minAbsColscale() const
{
   const DataArray < int >& colscaleExp = *m_activeColscaleExp;

   Real mini = infinity;

   for( int i = 0; i < colscaleExp.size(); ++i )
      if( spxAbs(spxLdexp(1.0, colscaleExp[i])) < mini )
         mini = spxAbs(spxLdexp(1.0, colscaleExp[i]));

   return mini;
}

Real SPxScaler::maxAbsColscale() const
{
   const DataArray < int >& colscaleExp = *m_activeColscaleExp;

   Real maxi = 0.0;

   for( int i = 0; i < colscaleExp.size(); ++i )
      if( spxAbs(spxLdexp(1.0, colscaleExp[i])) > maxi )
         maxi = spxAbs(spxLdexp(1.0, colscaleExp[i]));


   return maxi;
}

Real SPxScaler::minAbsRowscale() const
{
   const DataArray < int >& rowscaleExp = *m_activeRowscaleExp;

   int mini = std::numeric_limits<int>::max();

   for( int i = 0; i < rowscaleExp.size(); ++i )
      if( rowscaleExp[i] < mini )
         mini = rowscaleExp[i];

   return spxLdexp(1.0, mini);
}

Real SPxScaler::maxAbsRowscale() const
{
   const DataArray < int >& rowscaleExp = *m_activeRowscaleExp;

   int maxi = std::numeric_limits<int>::min();

   for( int i = 0; i < rowscaleExp.size(); ++i )
      if( rowscaleExp[i] > maxi )
         maxi = rowscaleExp[i];

   return spxLdexp(1.0, maxi);
}

/** \f$\max_{j\in\mbox{ cols}}
 *   \left(\frac{\max_{i\in\mbox{ rows}}|a_ij|}
 *              {\min_{i\in\mbox{ rows}}|a_ij|}\right)\f$
 */
Real SPxScaler::maxColRatio(const SPxLP& lp) const
{

   Real pmax = 0.0;

   for(int i = 0; i < lp.nCols(); ++i )
   {
      const SVector& vec  = lp.colVector(i);
      Real           mini = infinity;
      Real           maxi = 0.0;

      for(int j = 0; j < vec.size(); ++j)
      {
         Real x = spxAbs(vec.value(j));

         if( isZero(x) )
            continue;
         if( x < mini )
            mini = x;
         if( x > maxi )
            maxi = x;
      }

      if( mini == infinity )
         continue;

      Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }
   return pmax;
}

/** \f$\max_{i\in\mbox{ rows}}
 *   \left(\frac{\max_{j\in\mbox{ cols}}|a_ij|}
 *              {\min_{j\in\mbox{ cols}}|a_ij|}\right)\f$
 */
Real SPxScaler::maxRowRatio(const SPxLP& lp) const
{

   Real pmax = 0.0;

   for(int i = 0; i < lp.nRows(); ++i )
   {
      const SVector& vec  = lp.rowVector(i);
      Real           mini = infinity;
      Real           maxi = 0.0;

      for(int j = 0; j < vec.size(); ++j)
      {
         Real x = spxAbs(vec.value(j));

         if( isZero(x) )
            continue;
         if( x < mini )
            mini = x;
         if( x > maxi )
            maxi = x;
      }

      if( mini == infinity )
         continue;

      Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }
   return pmax;
}

void SPxScaler::computeExpVec(const std::vector<Real>& vec, DataArray<int>& vecExp)
{
   assert(vec.size() == unsigned(vecExp.size()));

   for( unsigned i = 0; i < vec.size(); ++i )
   {
       frexp(vec[i], &(vecExp[int(i)]));
       vecExp[int(i)] -= 1;
   }
}

bool SPxScaler::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   return m_activeColscaleExp->isConsistent() && m_activeRowscaleExp->isConsistent();
#else
   return true;
#endif
}


} // namespace soplex
