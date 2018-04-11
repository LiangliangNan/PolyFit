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
#include "spxhybridpr.h"
#include "spxout.h"

namespace soplex
{

bool SPxHybridPR::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (thesolver != 0 &&
      (thesolver != steep.solver() ||
         thesolver != devex.solver() ||
         thesolver != parmult.solver()))
      return MSGinconsistent("SPxHybridPR");

   return steep.isConsistent()
          && devex.isConsistent()
          && parmult.isConsistent();
#else
   return true;
#endif
}

void SPxHybridPR::load(SPxSolver* p_solver)
{
   steep.load(p_solver);
   devex.load(p_solver);
   parmult.load(p_solver);
   thesolver = p_solver;
   setType(p_solver->type());
}

void SPxHybridPR::clear()
{
   steep.clear();
   devex.clear();
   parmult.clear();
   thesolver = 0;
}

void SPxHybridPR::setEpsilon(Real eps)
{
   steep.setEpsilon(eps);
   devex.setEpsilon(eps);
   parmult.setEpsilon(eps);
}

void SPxHybridPR::setType(SPxSolver::Type tp)
{
   if (tp == SPxSolver::LEAVE)
   {
      thepricer = &steep;
      thesolver->setPricing(SPxSolver::FULL);
   }
   else
   {
      if (thesolver->dim() > hybridFactor * thesolver->coDim())
      {
         /**@todo I changed from devex to steepest edge pricing here 
          *       because of numerical difficulties, this should be 
          *       investigated.
          */
         // thepricer = &devex;
         thepricer = &steep;
         thesolver->setPricing(SPxSolver::FULL);
      }
      else
      {
         thepricer = &parmult;
         thesolver->setPricing(SPxSolver::PARTIAL);
      }
   }
   
   MSG_INFO1( (*thesolver->spxout), (*thesolver->spxout) << "IPRHYB01 switching to "
                        << thepricer->getName() << std::endl; )

   thepricer->setType(tp);
}

void SPxHybridPR::setRep(SPxSolver::Representation rep)
{
   steep.setRep(rep);
   devex.setRep(rep);
   parmult.setRep(rep);
}

int SPxHybridPR::selectLeave()
{
   return thepricer->selectLeave();
}

void SPxHybridPR::left4(int n, SPxId id)
{
   thepricer->left4(n, id);
}

SPxId SPxHybridPR::selectEnter()
{
   return thepricer->selectEnter();
}

void SPxHybridPR::entered4(SPxId id, int n)
{
   thepricer->entered4(id, n);
}

void SPxHybridPR::addedVecs (int n)
{
   steep.addedVecs(n);
   devex.addedVecs(n);
   parmult.addedVecs(n);
}

void SPxHybridPR::addedCoVecs(int n)
{
   steep.addedCoVecs(n);
   devex.addedCoVecs(n);
   parmult.addedCoVecs(n);
}

} // namespace soplex
