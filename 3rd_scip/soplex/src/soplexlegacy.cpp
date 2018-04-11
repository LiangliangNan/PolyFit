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

#include "soplexlegacy.h"
#include "exceptions.h"

namespace soplex
{
SoPlexLegacy::SoPlexLegacy(SPxOut& outstream, SPxSolver::Type p_type, SPxSolver::Representation p_rep)
   : m_solver(p_type, p_rep)
   , m_postScaler(0)
   , m_simplifier(0)
   , m_vanished(false)
   , m_freePostScaler(false)
   , m_freeSimplifier(false)
{
   m_solver.setOutstream(outstream);
   m_solver.setBasisSolver(&m_slu);
   m_solver.setTester(new SPxBoundFlippingRT(), true);
   m_solver.setPricer(new SPxSteepPR(), true);
   m_solver.setStarter(0);

   assert(SoPlexLegacy::isConsistent());
}

SoPlexLegacy::~SoPlexLegacy()
{
   assert(!m_freePostScaler || m_postScaler != 0);
   assert(!m_freeSimplifier || m_simplifier != 0);

   if(m_freePostScaler)
   {
      delete m_postScaler;
      m_postScaler = 0;
   }

   if(m_freeSimplifier)
   {
      delete m_simplifier;
      m_simplifier = 0;
   }
}

SoPlexLegacy& SoPlexLegacy::operator=(const SoPlexLegacy& base)
{
   assert(!m_freePostScaler || m_postScaler != 0);
   assert(!m_freeSimplifier || m_simplifier != 0);

   if(this != &base)
   {
      SPxLP::operator=(base);
      m_slu = base.m_slu;  // call of SLinSolver::clone() SPxBasis assignment operator not necessary (done by m_solver.setSolver(&m_slu) below)
      m_solver = base.m_solver;
      m_vanished = base.m_vanished;
      m_solver.setBasisSolver(&m_slu);

      // m_postScaler
      if(m_freePostScaler)
      {
         delete m_postScaler;
         m_postScaler = 0;
      }
      if(base.m_postScaler == 0)
      {
         m_postScaler = 0;
         m_freePostScaler = false;
      }
      else
      {
         m_postScaler = base.m_postScaler->clone();
         m_freePostScaler = true;
      }

      // m_simplifier
      if(m_freeSimplifier)
      {
         delete m_simplifier;
         m_simplifier = 0;
      }
      if(base.m_simplifier == 0)
      {
         m_simplifier = 0;
         m_freeSimplifier = false;
      }
      else
      {
         m_simplifier = base.m_simplifier->clone();
         m_freeSimplifier = true;
      }


   }

   return *this;
}


SoPlexLegacy::SoPlexLegacy(const SoPlexLegacy& old)
   : SPxLP(old)
   , m_slu(old.m_slu)  // call of SLinSolver::clone() SPxBasis copy constructor not necessary (done by m_solver.setSolver(&m_slu) below)
   , m_solver(old.m_solver)
   , m_vanished(old.m_vanished)
{
   m_solver.setBasisSolver(&m_slu);

   // m_postScaler
   if(old.m_postScaler == 0)
   {
      m_postScaler = 0;
      m_freePostScaler = false;
   }
   else
   {
      m_postScaler = old.m_postScaler->clone();
      m_freePostScaler = true;
   }

   // m_simplifier
   if(old.m_simplifier == 0)
   {
      m_simplifier = 0;
      m_freeSimplifier = false;
   }
   else
   {
      m_simplifier = old.m_simplifier->clone();
      m_freeSimplifier = true;
   }
}



void SoPlexLegacy::setPostScaler(SPxScaler* x, const bool destroy)
{

   assert(!m_freePostScaler || m_postScaler != 0);

   if(m_freePostScaler)
   {
      delete m_postScaler;
      m_postScaler = 0;
   }
   m_postScaler = x;
   if( m_postScaler )
      m_postScaler->setOutstream(*spxout);
   m_freePostScaler = destroy;
}

void SoPlexLegacy::setSimplifier(SPxSimplifier* x, const bool destroy)
{

   assert(!m_freeSimplifier || m_simplifier != 0);

   if(m_freeSimplifier)
   {
      delete m_simplifier;
      m_simplifier = 0;
   }
   m_simplifier = x;
   m_freeSimplifier = destroy;
}

Real SoPlexLegacy::objValue() const
{

   DVector x(nCols());

   getPrimal(x);

   return x * maxObj() * Real(spxSense());
}

SPxSolver::Status SoPlexLegacy::solve()
{

   if (nRows() <= 0 && nCols() <= 0) // no problem loaded
      throw SPxStatusException("XSOLVR01 No Problem loaded");

   // assume presolver did NOT solve problem
   m_vanished = false;

   // working LP
   SPxLP work(*this);
   work.setOutstream(*spxout);

   // keep useful bounds for bfrt
   bool keepbounds = (!strcmp(m_solver.ratiotester()->getName(), "Bound Flipping"));

   // should the LP be simplified ?
   if (m_simplifier != 0)
   {
      switch(m_simplifier->simplify(work, m_solver.epsilon(), m_solver.feastol(), m_solver.opttol(), keepbounds))
      {
      case SPxSimplifier::UNBOUNDED :
         m_solver.setBasisStatus(SPxBasis::UNBOUNDED);
         return SPxSolver::UNBOUNDED;
      case SPxSimplifier::INFEASIBLE :
         m_solver.setBasisStatus(SPxBasis::INFEASIBLE);
         return SPxSolver::INFEASIBLE;
      case SPxSimplifier::VANISHED :
         m_vanished = true;
         return SPxSolver::OPTIMAL;
      case SPxSimplifier::OKAY:
         break;
      default:
         throw SPxInternalCodeException("XRSOLVR01 This should never happen.");
      }
   }

   // should the LP be scaled after simplifing?
   if (m_postScaler != 0)
      m_postScaler->scale(work, false);

   /* if a basis of correct size was set externally, try to load it into the transformed LP */
   if ( m_colsbasisstatus.size() == work.nCols() && m_rowsbasisstatus.size() == work.nRows() )
   {
      m_solver.loadLP(work);
      if ( m_solver.isBasisValid(m_rowsbasisstatus, m_colsbasisstatus) )
      {
         m_solver.setBasis(m_rowsbasisstatus.get_const_ptr(), m_colsbasisstatus.get_const_ptr());
      }
      m_colsbasisstatus.clear();
      m_rowsbasisstatus.clear();
   }
   /* if a basis with status at least REGULAR exists (loaded via readBasisFile()
    * or available from previous simplex run), we check whether it can be (re)used
    * for the newly loaded LP.
    */
   else if ( m_solver.basis().status() <= SPxBasis::SINGULAR )
   {
      m_solver.loadLP(work);
   }
   else
   {
      SPxBasis::Desc oldbasisdesc(m_solver.basis().desc());
      m_solver.loadLP(work);
      if(m_solver.basis().isDescValid(oldbasisdesc))
         m_solver.loadBasis(oldbasisdesc);
   }

   return m_solver.solve();
}

SPxSolver::Status SoPlexLegacy::getPrimal(Vector& x) const
{

   if (has_simplifier())
   {
      if (!m_simplifier->isUnsimplified())
         unsimplify();

      x = m_simplifier->unsimplifiedPrimal();

      if (m_vanished)
         return SPxSolver::OPTIMAL;
      else
         return m_solver.status();
   }

   // else
   SPxSolver::Status stat = m_solver.getPrimal(x);

   // unscale postscaling
   if (m_postScaler != 0)
      m_postScaler->unscalePrimal(m_solver, x);

   return stat;
}

SPxSolver::Status SoPlexLegacy::getSlacks(Vector& s) const
{

   if (has_simplifier())
   {
      if (!m_simplifier->isUnsimplified())
         unsimplify();

      s = m_simplifier->unsimplifiedSlacks();

      if (m_vanished)
         return SPxSolver::OPTIMAL;
      else
         return m_solver.status();
   }

   // else
   SPxSolver::Status stat = m_solver.getSlacks(s);

   // unscale postscaling
   if (m_postScaler != 0)
      m_postScaler->unscaleSlacks(m_solver, s);

   return stat;
}

SPxSolver::Status SoPlexLegacy::getDual(Vector& pi) const
{

   if (has_simplifier())
   {
      if (!m_simplifier->isUnsimplified())
         unsimplify();

      pi = m_simplifier->unsimplifiedDual();

      if (m_vanished)
         return SPxSolver::OPTIMAL;
      else
         return m_solver.status();
   }

   // else
   SPxSolver::Status stat = m_solver.getDual(pi);

   // unscale postscaling
   if (m_postScaler != 0)
      m_postScaler->unscaleDual(m_solver, pi);

   return stat;
}

SPxSolver::Status SoPlexLegacy::getRedCost(Vector& rdcost) const
{

   if (has_simplifier())
   {
      if (!m_simplifier->isUnsimplified())
         unsimplify();

      rdcost = m_simplifier->unsimplifiedRedCost();

      if (m_vanished)
         return SPxSolver::OPTIMAL;
      else
         return m_solver.status();
   }

   // else
   SPxSolver::Status stat = m_solver.getRedCost(rdcost);

   // unscale postscaling
   if (m_postScaler != 0)
      m_postScaler->unscaleRedCost(m_solver, rdcost);

   return stat;
}

SPxSolver::VarStatus SoPlexLegacy::getBasisRowStatus(int i) const
{
   SPxBasis::SPxStatus b_status = m_solver.getBasisStatus();
   if((b_status == SPxBasis::NO_PROBLEM || (has_simplifier() && b_status == SPxBasis::SINGULAR)) && !m_vanished)
      return SPxSolver::UNDEFINED;

   if (has_simplifier())
   {
      if (!m_simplifier->isUnsimplified())
         unsimplify();

      return m_simplifier->getBasisRowStatus(i);
   }
   else
      return m_solver.getBasisRowStatus(i);
}

SPxSolver::VarStatus SoPlexLegacy::getBasisColStatus(int j) const
{
   SPxBasis::SPxStatus b_status = m_solver.getBasisStatus();
   if((b_status == SPxBasis::NO_PROBLEM || (has_simplifier() && b_status == SPxBasis::SINGULAR)) && !m_vanished)
      return SPxSolver::UNDEFINED;

   if (has_simplifier())
   {
      if (!m_simplifier->isUnsimplified())
         unsimplify();

      return m_simplifier->getBasisColStatus(j);
   }
   else
      return m_solver.getBasisColStatus(j);
}

SPxSolver::Status SoPlexLegacy::getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const
{
   SPxBasis::SPxStatus b_status = m_solver.getBasisStatus();
   if((b_status == SPxBasis::NO_PROBLEM || (has_simplifier() && b_status == SPxBasis::SINGULAR)) && !m_vanished)
   {
      int i;

      if (cols)
         for (i = nCols() - 1; i >= 0; --i)
            cols[i] = SPxSolver::UNDEFINED;

      if (rows)
         for (i = nRows() - 1; i >= 0; --i)
            rows[i] = SPxSolver::UNDEFINED;

      return m_solver.status();
   }

   if (has_simplifier())
   {
      if (!m_simplifier->isUnsimplified())
         unsimplify();

      m_simplifier->getBasis(rows, cols);
      return m_solver.status();
   }

   else
      return m_solver.getBasis(rows, cols);
}

SPxSolver::Status SoPlexLegacy::getPrimalray(Vector& primalray) const
{
   /// Does not work yet with presolve
   if (has_simplifier())
   {
      MSG_ERROR( std::cerr << "ESOLVR02 Primal ray with presolving not yet implemented" << std::endl; )
      throw SPxStatusException("XSOLVR02 Primal ray with presolving not yet implemented");
   }
   SPxSolver::Status stat = m_solver.getPrimalray(primalray);

   if (m_postScaler != 0)
      m_postScaler->unscalePrimal(m_solver, primalray);

   return stat;
}

SPxSolver::Status SoPlexLegacy::getDualfarkas(Vector& dualfarkas) const
{
   /// Does not work yet with presolve
   if (has_simplifier())
   {
      MSG_ERROR( std::cerr << "ESOLVR02 Dual farkas with presolving not yet implemented" << std::endl; )
      throw SPxStatusException("XSOLVR03 Dual farkas with presolving not yet implemented");
      //      return SPxSolver::ERROR;
   }
   SPxSolver::Status stat = m_solver.getDualfarkas(dualfarkas);

   if (m_postScaler != 0)
      m_postScaler->unscaleDual(m_solver, dualfarkas);

   return stat;
}

void SoPlexLegacy::qualConstraintViolation(
   Real& maxviol,
   Real& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   DVector solu( nCols() );

   getPrimal( solu );

   for( int row = 0; row < nRows(); ++row )
   {
      const SVector& rowvec = rowVector( row );

      Real val = 0.0;

      for( int col = 0; col < rowvec.size(); ++col )
         val += rowvec.value( col ) * solu[rowvec.index( col )];

      Real viol = 0.0;

      assert(lhs( row ) <= rhs( row ));

      if (val < lhs( row ))
         viol = spxAbs(val - lhs( row ));
      else
         if (val > rhs( row ))
            viol = spxAbs(val - rhs( row ));

      if (viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

void SoPlexLegacy::qualBoundViolation(
   Real& maxviol,
   Real& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   DVector solu( nCols() );

   getPrimal( solu );

   for( int col = 0; col < nCols(); ++col )
   {
      assert( lower( col ) <= upper( col ));

      Real viol = 0.0;

      if (solu[col] < lower( col ))
         viol = spxAbs( solu[col] - lower( col ));
      else
         if (solu[col] > upper( col ))
            viol = spxAbs( solu[col] - upper( col ));

      if (viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

bool SoPlexLegacy::readBasisFile(
   const char*    filename,
   const NameSet* rowNames,
   const NameSet* colNames
   )
{
   // init solver using original LP
   m_solver.loadLP(*this);
   return m_solver.readBasisFile(filename, rowNames, colNames);
}

bool SoPlexLegacy::writeBasisFile(
   const char*    filename,
   const NameSet* rowNames,
   const NameSet* colNames
   )
{
   assert(rep() == SPxSolver::COLUMN);

   /* make sure the basis of the original problem is written */
   unsimplify();

   std::ofstream file(filename);
   std::ostream& os = file;

   os.setf(std::ios::left);
   os << "NAME  soplex.bas\n";

   /* do not write basis if there is none */
   if( status() == SPxSolver::NO_PROBLEM )
   {
      os << "ENDATA" << std::endl;
      return true;
   }

   /* start writing */
   char buf[SPX_MAXSTRLEN];
   int row = 0;
   for( int col = 0; col < nCols(); col++ )
   {
      if( getBasisColStatus(col) == SPxSolver::BASIC )
      {
         /* Find non basic row */
         for( ; row < nRows(); row++ )
         {
            if( getBasisRowStatus(row) != SPxSolver::BASIC )
               break;
         }

         assert(row != nRows());

         os << ( getBasisRowStatus(row) == SPxSolver::ON_UPPER ? " XU " : " XL " )
         << std::setw(8) << getColName(col, colNames, buf);

         /* break in two parts since buf is reused */
         os << "       "
         << getRowName(row, rowNames, buf)
         << std::endl;

         row++;
      }
      else
      {
         if ( getBasisColStatus(col) == SPxSolver::ON_UPPER )
         {
            os << " UL "
            << getColName(col, colNames, buf)
            << std::endl;
         }
      }
   }

   #ifndef NDEBUG
   // Check that we covered all nonbasic rows - the remaining should be basic.
   for( ; row < nRows(); row++ )
   {
      if( getBasisRowStatus(row) != SPxSolver::BASIC )
         break;
   }
   assert(row == nRows());
   #endif // NDEBUG

   os << "ENDATA" << std::endl;
   return true;
}

bool SoPlexLegacy::writeState(
   const char*    filename,
   const NameSet* rowNames,
   const NameSet* colNames ) const
{

   return m_solver.writeState(filename, rowNames, colNames);
}

void SoPlexLegacy::unsimplify() const
{
   if (!has_simplifier() || m_simplifier->isUnsimplified())
      return;

   DVector psp_x(m_solver.nCols());  // primal solution (prescaled simplified postscaled)
   DVector psp_y(m_solver.nRows());  // dual solution   (prescaled simplified postscaled)
   DVector psp_s(m_solver.nRows());  // slacks          (prescaled simplified postscaled)
   DVector psp_r(m_solver.nCols());  // reduced costs   (prescaled simplified postscaled)

   if (! m_vanished) {
      // If solver status is not regular or optimal, do nothing.
      const SPxSolver::Status stat = status();
      if( stat != SPxSolver::OPTIMAL && stat != SPxSolver::REGULAR )
         return;

      m_solver.getPrimal(psp_x);
      m_solver.getDual(psp_y);
      m_solver.getSlacks(psp_s);
      m_solver.getRedCost(psp_r);

      // unscale postscaling
      if (m_postScaler != 0)
      {
         m_postScaler->unscalePrimal(m_solver, psp_x);
         m_postScaler->unscaleDual(m_solver, psp_y);
         m_postScaler->unscaleSlacks(m_solver, psp_s);
         m_postScaler->unscaleRedCost(m_solver, psp_r);
      }
   }
   else {
      // If there is no sensible solution, do nothing.
      const SPxSolver::Status  stat = status();
      if (stat != SPxSolver::OPTIMAL)
         return;

      psp_x.reDim(0);
      psp_y.reDim(0);
      psp_s.reDim(0);
      psp_r.reDim(0);
   }

   // unsimplify
   if(m_vanished)
   {
      m_simplifier->unsimplify(psp_x, psp_y, psp_s, psp_r, NULL, NULL);
   }
   else
   {
      SPxSolver::VarStatus *rows, *cols;
      rows = NULL;
      cols = NULL;
      try
      {
         spx_alloc(rows, m_solver.nRows());
         spx_alloc(cols, m_solver.nCols());

         m_solver.getBasis(rows, cols, m_solver.nRows(), m_solver.nCols());
         m_simplifier->unsimplify(psp_x, psp_y, psp_s, psp_r, rows, cols);
      }
      catch( const SPxMemoryException& x)
      {
         spx_free(rows);
         spx_free(cols);
         throw x;
      }

      spx_free(rows);
      spx_free(cols);
   }
}
} // namespace soplex
