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

/**@file  validation.cpp
 * @brief Validation object for soplex solutions
 */

#include "validation.h"

namespace soplex {

/// updates the external solution used for validation
bool Validation::updateExternalSolution(char* solution)
{
   validate = true;
   validatesolution = solution;

   if( strncmp(solution, "+infinity", 9 ) == 0 )
      return true;
   else if ( strncmp(solution, "-infinity", 9) == 0 )
      return true;
   else
   {
      char* tailptr;
      strtod(solution, &tailptr);
      if (*tailptr) {
         //conversion failed because the input wasn't a number
         return false;
      }
   }
   return true;
}



/// updates the tolerance used for validation
bool Validation::updateValidationTolerance(char* tolerance)
{
   char* tailptr;
   validatetolerance = strtod(tolerance, &tailptr);
   if (*tailptr) {
      //conversion failed because the input wasn't a number
      return false;
   }
   return true;
}



/// validates the soplex solution using the external solution
void Validation::validateSolveReal(SoPlex& soplex)
{
#ifndef SOPLEX_LEGACY
   bool passedValidation = true;
   std::string reason = "";
   Real objViolation = 0.0;
   Real maxBoundViolation = 0.0;
   Real maxRowViolation = 0.0;
   Real maxRedCostViolation = 0.0;
   Real maxDualViolation = 0.0;
   Real sumBoundViolation = 0.0;
   Real sumRowViolation = 0.0;
   Real sumRedCostViolation = 0.0;
   Real sumDualViolation = 0.0;
   Real sol;

   std::ostream& os = soplex.spxout.getStream(SPxOut::INFO1);

   if( strncmp(validatesolution, "+infinity", 9 ) == 0 )
      sol =  soplex.realParam(SoPlex::INFTY);
   else if ( strncmp(validatesolution, "-infinity", 9) == 0 )
      sol =  -soplex.realParam(SoPlex::INFTY);
   else
   {
      sol = atof(validatesolution);
   }

   objViolation = spxAbs(sol - soplex.objValueReal());
   if( ! EQ(objViolation, 0.0, validatetolerance) )
   {
      passedValidation = false;
      reason += "Objective Violation; ";
   }
   if( SPxSolver::OPTIMAL == soplex.status() )
   {
      soplex.getBoundViolationReal(maxBoundViolation, sumBoundViolation);
      soplex.getRowViolationReal(maxRowViolation, sumRowViolation);
      soplex.getRedCostViolationReal(maxRedCostViolation, sumRedCostViolation);
      soplex.getDualViolationReal(maxDualViolation, sumDualViolation);
      if( ! LE(maxBoundViolation, validatetolerance) )
      {
         passedValidation = false;
         reason += "Bound Violation; ";
      }
      if( ! LE(maxRowViolation, validatetolerance) )
      {
         passedValidation = false;
         reason += "Row Violation; ";
      }
      if( ! LE(maxRedCostViolation, validatetolerance) )
      {
         passedValidation = false;
         reason += "Reduced Cost Violation; ";
      }
      if( ! LE(maxDualViolation, validatetolerance) )
      {
         passedValidation = false;
         reason += "Dual Violation; ";
      }
   }

   os << "\n";
   os << "Validation          :";
   if(passedValidation)
      os << " Success\n";
   else
   {
      reason[reason.length()-2] = ']';
      os << " Fail [" + reason + "\n";
   }
   os << "   Objective        : " << std::scientific << std::setprecision(8) << objViolation << std::fixed << "\n";
   os << "   Bound            : " << std::scientific << std::setprecision(8) << maxBoundViolation << std::fixed << "\n";
   os << "   Row              : " << std::scientific << std::setprecision(8) << maxRowViolation << std::fixed << "\n";
   os << "   Reduced Cost     : " << std::scientific << std::setprecision(8) << maxRedCostViolation << std::fixed << "\n";
   os << "   Dual             : " << std::scientific << std::setprecision(8) << maxDualViolation << std::fixed << "\n";
#endif
}

} /* namespace soplex */
