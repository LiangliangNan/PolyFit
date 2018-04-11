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
#include "spxautopr.h"
#include "spxout.h"

namespace soplex
{

void SPxAutoPR::load(SPxSolver* p_solver)
{
   steep.load(p_solver);
   devex.load(p_solver);
   thesolver = p_solver;
   setType(p_solver->type());
}

void SPxAutoPR::clear()
{
   steep.clear();
   devex.clear();
   thesolver = 0;
}

void SPxAutoPR::setEpsilon(Real eps)
{
   steep.setEpsilon(eps);
   devex.setEpsilon(eps);
   theeps = eps;
}

void SPxAutoPR::setType(SPxSolver::Type tp)
{
   activepricer->setType(tp);
}

void SPxAutoPR::setRep(SPxSolver::Representation rep)
{
   steep.setRep(rep);
   devex.setRep(rep);
}

bool SPxAutoPR::setActivePricer(SPxSolver::Type type)
{
   // switch to steep as soon as switchIters is reached
   if( activepricer == &devex && thesolver->iterations() >= switchIters )
   {
      activepricer = &steep;
      activepricer->setType(type);
      return true;
   }

   // use devex for the iterations < switchIters
   else if( activepricer == &steep && thesolver->iterations() < switchIters  )
   {
      activepricer = &devex;
      activepricer->setType(type);
      return true;
   }

   return false;
}

int SPxAutoPR::selectLeave()
{
   if( setActivePricer(SPxSolver::LEAVE) )
      MSG_INFO1( (*thesolver->spxout), (*thesolver->spxout) << " --- active pricer: " << activepricer->getName() << std::endl; )

   return activepricer->selectLeave();
}

void SPxAutoPR::left4(int n, SPxId id)
{
   activepricer->left4(n, id);
}

SPxId SPxAutoPR::selectEnter()
{
   if( setActivePricer(SPxSolver::ENTER) )
      MSG_INFO1( (*thesolver->spxout), (*thesolver->spxout) << " --- active pricer: " << activepricer->getName() << std::endl; )

   return activepricer->selectEnter();
}

void SPxAutoPR::entered4(SPxId id, int n)
{
   activepricer->entered4(id, n);
}

} // namespace soplex
