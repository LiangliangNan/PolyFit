/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_validation.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for validation
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_validation.h"

/** validate the result of the solve
 *
 *  the validation includes
 *
 *  - checking the feasibility of the incumbent solution in the original problem (using SCIPcheckSolOrig())
 *
 *  - checking if the objective bounds computed by SCIP agree with external primal and dual reference bounds.
 *
 *  All external reference bounds the original problem space and the original objective sense.
 *
 *  For infeasible problems, +/-SCIPinfinity() should be passed as reference bounds depending on the objective sense
 *  of the original problem.
 */
SCIP_RETCODE SCIPvalidateSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             primalreference,    /**< external primal reference value for the problem, or SCIP_UNKNOWN */
   SCIP_Real             dualreference,      /**< external dual reference value for the problem, or SCIP_UNKNOWN */
   SCIP_Real             reftol,             /**< relative tolerance for acceptable violation of reference values */
   SCIP_Bool             quiet,              /**< TRUE if no status line should be printed */
   SCIP_Bool*            feasible,           /**< pointer to store if the best solution is feasible in the original problem,
                                              *   or NULL */
   SCIP_Bool*            primalboundcheck,   /**< pointer to store if the primal bound respects the given dual reference
                                              *   value, or NULL */
   SCIP_Bool*            dualboundcheck      /**< pointer to store if the dual bound respects the given primal reference
                                              *   value, or NULL */
   )
{
   SCIP_Bool localfeasible;
   SCIP_Bool localprimalboundcheck;
   SCIP_Bool localdualboundcheck;
   SCIP_Real primviol;
   SCIP_Real dualviol;
   assert(scip != NULL);

   /* if no problem exists, there is no need for validation */
   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      if( feasible != NULL )
         *feasible = TRUE;
      if( primalboundcheck != NULL )
         *primalboundcheck = TRUE;
      if( dualboundcheck != NULL )
         *dualboundcheck = TRUE;

      return SCIP_OKAY;
   }

   localfeasible = TRUE;
   localdualboundcheck = TRUE;

   /* check the best solution for feasibility in the original problem */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_SOL* bestsol = SCIPgetBestSol(scip);
      SCIP_Real checkfeastolfac;
      SCIP_Real oldfeastol;

      assert(bestsol != NULL);

      /* scale feasibility tolerance by set->num_checkfeastolfac */
      oldfeastol = SCIPfeastol(scip);
      SCIP_CALL( SCIPgetRealParam(scip, "numerics/checkfeastolfac", &checkfeastolfac) );
      if( !SCIPisEQ(scip, checkfeastolfac, 1.0) )
      {
         SCIP_CALL( SCIPchgFeastol(scip, oldfeastol * checkfeastolfac) );
      }

      SCIP_CALL( SCIPcheckSolOrig(scip, bestsol, &localfeasible, !quiet, TRUE) );

      /* restore old feasibilty tolerance */
      if( !SCIPisEQ(scip, checkfeastolfac, 1.0) )
      {
         SCIP_CALL( SCIPchgFeastol(scip, oldfeastol) );
      }
   }
   else
   {
      localfeasible = TRUE;
   }

   primviol = 0.0;
   dualviol = 0.0;
   /* check the primal and dual bounds computed by SCIP against the external reference values within reference tolerance */
   /* solution for an infeasible problem */
   if( SCIPgetNSols(scip) > 0 && ((SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE && SCIPisInfinity(scip, dualreference))
            || (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE && SCIPisInfinity(scip, -dualreference))) )
      localprimalboundcheck = FALSE;
   else
   {
      /* check if reference primal bound is not better than the proven dual bound and, if SCIP claims to be optimal,
       * if the
       */
      SCIP_Real pb = SCIPgetPrimalbound(scip);
      SCIP_Real db = SCIPgetDualbound(scip);

      /* compute the relative violation between the primal bound and dual reference value, and vice versa */
      if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
      {
         if( dualreference != SCIP_UNKNOWN ) /*lint !e777 */
            primviol = SCIPrelDiff(dualreference, pb);
         if( primalreference != SCIP_UNKNOWN ) /*lint !e777 */
            dualviol = SCIPrelDiff(db, primalreference);
      }
      else
      {
         if( dualreference != SCIP_UNKNOWN ) /*lint !e777 */
            primviol = SCIPrelDiff(pb, dualreference);

         if( primalreference != SCIP_UNKNOWN ) /*lint !e777 */
            dualviol = SCIPrelDiff(primalreference, db);
      }
      primviol = MAX(primviol, 0.0);
      dualviol = MAX(dualviol, 0.0);

      localprimalboundcheck = EPSP(reftol, primviol);
      localdualboundcheck = EPSP(reftol, dualviol);
   }

   if( !quiet )
   {
      SCIPinfoMessage(scip, NULL, "Validation         : ");
      if( ! localfeasible )
         SCIPinfoMessage(scip, NULL, "Fail (infeasible)");
      else if( ! localprimalboundcheck )
         SCIPinfoMessage(scip, NULL, "Fail (primal bound)");
      else if( ! localdualboundcheck )
         SCIPinfoMessage(scip, NULL, "Fail (dual bound)");
      else
         SCIPinfoMessage(scip, NULL, "Success");
      SCIPinfoMessage(scip, NULL, "\n");
      SCIPinfoMessage(scip, NULL, "  %-17s: %10u\n", "cons violation", !localfeasible); /*lint !e705*/
      SCIPinfoMessage(scip, NULL, "  %-17s: %10.8g (reference: %16.9e)\n", "primal violation", primviol, dualreference);
      SCIPinfoMessage(scip, NULL, "  %-17s: %10.8g (reference: %16.9e)\n", "dual violation", dualviol, primalreference);
   }

   if( feasible != NULL )
      *feasible = localfeasible;
   if( primalboundcheck != NULL )
      *primalboundcheck = localprimalboundcheck;
   if( dualboundcheck != NULL )
      *dualboundcheck = localdualboundcheck;

   return SCIP_OKAY;
}
