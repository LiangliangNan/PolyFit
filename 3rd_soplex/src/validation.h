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

/**@file  validation.h
 * @brief Validation object for soplex solutions
 */

#ifndef SRC_VALIDATION_H_
#define SRC_VALIDATION_H_

#include "soplex.h"

namespace soplex {

class Validation
{
public:

   /// should the soplex solution be validated?
   bool           validate;

   /// external solution used for validation
   char*          validatesolution;

   /// tolerance used for validation
   double         validatetolerance;

   /// default constructor
   Validation()
   {
      validate = false;
      validatetolerance = 1e-5;
      validatesolution = 0;
   }

   /// default destructor
   ~Validation()
   {
      ;
   }

   /// updates the external solution used for validation
   bool updateExternalSolution(char* solution);

   /// updates the tolerance used for validation
   bool updateValidationTolerance(char* tolerance);

   /// validates the soplex solution using the external solution
   void validateSolveReal(SoPlex& soplex);
};

} /* namespace soplex */

#endif /* SRC_VALIDATION_H_ */
