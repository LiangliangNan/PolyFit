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

/**@file   scip_general.c
 * @ingroup OTHER_CFILES
 * @brief  general public methods
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

#include "blockmemshell/memory.h"
#include "lpi/lpi.h"
#include "scip/exprinterpret.h"
#include "scip/clock.h"
#include "scip/debug.h"
#include "scip/dialog.h"
#include "scip/interrupt.h"
#include "scip/mem.h"
#include "scip/message_default.h"
#include "scip/nlp.h"
#include "scip/pub_message.h"
#include "scip/retcode.h"
#include "scip/scipbuildflags.h"
#include "scip/scipcoreplugins.h"
#include "scip/scip_general.h"
#include "scip/scipgithash.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_solvingstats.h"
#include "scip/set.h"
#include "scip/solve.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/syncstore.h"

#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

#ifdef SCIP_WITH_ZLIB
#include <zlib.h>
#endif

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPgetStage
#undef SCIPhasPerformedPresolve
#undef SCIPisStopped

/** returns complete SCIP version number in the format "major . minor tech"
 *
 *  @return complete SCIP version
 */
SCIP_Real SCIPversion(
   void
   )
{
   return (SCIP_Real)(SCIP_VERSION)/100.0;
}

/** returns SCIP major version
 *
 *  @return major SCIP version
 */
int SCIPmajorVersion(
   void
   )
{
   return SCIP_VERSION/100;
}

/** returns SCIP minor version
 *
 *  @return minor SCIP version
 */
int SCIPminorVersion(
   void
   )
{
   return (SCIP_VERSION/10) % 10; /*lint !e778*/
}

/** returns SCIP technical version
 *
 *  @return technical SCIP version
 */
int SCIPtechVersion(
   void
   )
{
   return SCIP_VERSION % 10; /*lint !e778*/
}

/** returns SCIP sub version number
 *
 *  @return subversion SCIP version
 */
int SCIPsubversion(
   void
   )
{
   return SCIP_SUBVERSION;
}

/** prints a version information line to a file stream via the message handler system
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
void SCIPprintVersion(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert( scip != NULL );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "SCIP version %d.%d.%d",
      SCIPmajorVersion(), SCIPminorVersion(), SCIPtechVersion());
#if SCIP_SUBVERSION > 0
   SCIPmessageFPrintInfo(scip->messagehdlr, file, ".%d", SCIPsubversion());
#endif

   SCIPmessageFPrintInfo(scip->messagehdlr, file, " [precision: %d byte]", (int)sizeof(SCIP_Real));

#ifndef BMS_NOBLOCKMEM
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " [memory: block]");
#else
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " [memory: standard]");
#endif
#ifndef NDEBUG
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " [mode: debug]");
#else
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " [mode: optimized]");
#endif
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " [LP solver: %s]", SCIPlpiGetSolverName());
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " [GitHash: %s]", SCIPgetGitHash());
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "%s\n", SCIP_COPYRIGHT);
}

/** prints detailed information on the compile-time flags
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
void SCIPprintBuildOptions(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert( scip != NULL );

   /* compiler */
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Compiler: ");
#if defined(__INTEL_COMPILER)
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Intel %d\n", __INTEL_COMPILER);
#elif defined(__clang__)
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "clang %d.%d.%d\n", __clang_major__, __clang_minor__, __clang_patchlevel__);
#elif defined(_MSC_VER)
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "microsoft visual c %d\n", _MSC_FULL_VER);
#elif defined(__GNUC__)
#if defined(__GNUC_PATCHLEVEL__)
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "gcc %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "gcc %d.%d\n", __GNUC__, __GNUC_MINOR__);
#endif
#else
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "unknown\n");
#endif

   /* build flags */
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "\nBuild options:\n%s", SCIPgetBuildFlags());
}

/** prints error message for the given SCIP_RETCODE via the error prints method */
void SCIPprintError(
   SCIP_RETCODE          retcode             /**< SCIP return code causing the error */
   )
{
   SCIPmessagePrintError("SCIP Error (%d): ", retcode);
   SCIPretcodePrintError(retcode);
   SCIPmessagePrintError("\n");
}

/*
 * general SCIP methods
 */

/** internal method to create SCIP */
static
SCIP_RETCODE doScipCreate(
   SCIP**                scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_ALLOC( BMSallocMemory(scip) );

   /* all members are initialized to NULL */
   BMSclearMemory(*scip);

   /* create a default message handler */
   SCIP_CALL( SCIPcreateMessagehdlrDefault(&(*scip)->messagehdlr, TRUE, NULL, FALSE) );

   SCIP_CALL( SCIPmemCreate(&(*scip)->mem) );
   SCIP_CALL( SCIPsetCreate(&(*scip)->set, (*scip)->messagehdlr, (*scip)->mem->setmem, *scip) );
   SCIP_CALL( SCIPinterruptCreate(&(*scip)->interrupt) );
   SCIP_CALL( SCIPdialoghdlrCreate((*scip)->set, &(*scip)->dialoghdlr) );
   SCIP_CALL( SCIPclockCreate(&(*scip)->totaltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPsyncstoreCreate( &(*scip)->syncstore ) );

   /* include additional core functionality */
   SCIP_CALL( SCIPincludeCorePlugins(*scip) );

   SCIPclockStart((*scip)->totaltime, (*scip)->set);

   SCIP_CALL( SCIPnlpInclude((*scip)->set, SCIPblkmem(*scip)) );

   if( strcmp(SCIPlpiGetSolverName(), "NONE") != 0 )
   {
      SCIP_CALL( SCIPsetIncludeExternalCode((*scip)->set, SCIPlpiGetSolverName(), SCIPlpiGetSolverDesc()) );
   }
   if( strcmp(SCIPexprintGetName(), "NONE") != 0 )
   {
      SCIP_CALL( SCIPsetIncludeExternalCode((*scip)->set, SCIPexprintGetName(), SCIPexprintGetDesc()) );
   }

#ifdef SCIP_WITH_ZLIB
   SCIP_CALL( SCIPsetIncludeExternalCode((*scip)->set, "ZLIB " ZLIB_VERSION, "General purpose compression library by J. Gailly and M. Adler (zlib.net)") );
#endif

   return SCIP_OKAY;
}

/** creates and initializes SCIP data structures
 *
 *  @note The SCIP default message handler is installed. Use the method SCIPsetMessagehdlr() to install your own
 *        message handler or SCIPsetMessagehdlrLogfile() and SCIPsetMessagehdlrQuiet() to write into a log
 *        file and turn off/on the display output, respectively.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @post After calling this method \SCIP reached the solving stage \ref SCIP_STAGE_INIT
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcreate(
   SCIP**                scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_FINALLY( doScipCreate(scip), (void)SCIPfree(scip) );

   return SCIP_OKAY;
}

/** frees SCIP data structures
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method \SCIP reached the solving stage \ref SCIP_STAGE_FREE
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPfree(
   SCIP**                scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);
   if( *scip == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcheckStage(*scip, "SCIPfree", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPfreeProb(*scip) );
   assert((*scip)->set->stage == SCIP_STAGE_INIT);

   /* switch stage to FREE */
   (*scip)->set->stage = SCIP_STAGE_FREE;

   SCIP_CALL( SCIPsyncstoreRelease(&(*scip)->syncstore) );
   SCIP_CALL( SCIPsetFree(&(*scip)->set, (*scip)->mem->setmem) );
   SCIP_CALL( SCIPdialoghdlrFree(*scip, &(*scip)->dialoghdlr) );
   SCIPclockFree(&(*scip)->totaltime);
   SCIPinterruptFree(&(*scip)->interrupt);
   SCIP_CALL( SCIPmemFree(&(*scip)->mem) );

   /* release message handler */
   SCIP_CALL( SCIPmessagehdlrRelease(&(*scip)->messagehdlr) );

   BMSfreeMemory(scip);

   return SCIP_OKAY;
}

#undef SCIPgetStage
#undef SCIPhasPerformedPresolve
#undef SCIPisStopped

/** returns current stage of SCIP
 *
 *  @return the current SCIP stage
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_STAGE SCIPgetStage(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->stage;
}

/** outputs SCIP stage and solution status if applicable via the message handler
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 *
 *  @note If limits have been changed between the solution and the call to this function, the status is recomputed and
 *        thus may to correspond to the original status.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPprintStage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintStage", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "initialization");
      break;
   case SCIP_STAGE_PROBLEM:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "problem creation / modification");
      break;
   case SCIP_STAGE_TRANSFORMING:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "problem transformation");
      break;
   case SCIP_STAGE_TRANSFORMED:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "problem transformed");
      break;
   case SCIP_STAGE_INITPRESOLVE:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "presolving is being initialized");
      break;
   case SCIP_STAGE_PRESOLVING:
      if( SCIPsolveIsStopped(scip->set, scip->stat, TRUE) )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "solving was interrupted [");
         SCIP_CALL( SCIPprintStatus(scip, file) );
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "]");
      }
      else
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "presolving process is running");
      break;
   case SCIP_STAGE_EXITPRESOLVE:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "presolving is being exited");
      break;
   case SCIP_STAGE_PRESOLVED:
      if( SCIPsolveIsStopped(scip->set, scip->stat, TRUE) )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "solving was interrupted [");
         SCIP_CALL( SCIPprintStatus(scip, file) );
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "]");
      }
      else
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "problem is presolved");
      break;
   case SCIP_STAGE_INITSOLVE:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "solving process initialization");
      break;
   case SCIP_STAGE_SOLVING:
      if( SCIPsolveIsStopped(scip->set, scip->stat, TRUE) )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "solving was interrupted [");
         SCIP_CALL( SCIPprintStatus(scip, file) );
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "]");
      }
      else
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "solving process is running");
      break;
   case SCIP_STAGE_SOLVED:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "problem is solved [");
      SCIP_CALL( SCIPprintStatus(scip, file) );
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "]");

      /* We output that the objective limit has been reached if no solution respecting the objective limit has been
       * found (nlimsolsfound == 0) and the primal bound is finite. Note that it still might be that the original
       * problem is infeasible, even without the objective limit, i.e., we cannot be sure that we actually reached the
       * objective limit. */
      if( scip->primal->nlimsolsfound == 0 && !SCIPisInfinity(scip, (SCIP_Real)SCIPgetObjsense(scip) * SCIPgetPrimalbound(scip))  )
         SCIPmessageFPrintInfo(scip->messagehdlr, file, " (objective limit reached)");

      break;
   case SCIP_STAGE_EXITSOLVE:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "solving process deinitialization");
      break;
   case SCIP_STAGE_FREETRANS:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "freeing transformed problem");
      break;
   case SCIP_STAGE_FREE:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "freeing SCIP");
      break;
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** gets solution status
 *
 *  @return SCIP solution status
 *
 *  See \ref SCIP_Status "SCIP_STATUS" for a complete list of all possible solving status.
 */
SCIP_STATUS SCIPgetStatus(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetStatus", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( scip->set->stage == SCIP_STAGE_INIT || scip->set->stage == SCIP_STAGE_FREE )
      return SCIP_STATUS_UNKNOWN;
   else
   {
      assert(scip->stat != NULL);

      return scip->stat->status;
   }
}

/** outputs solution status
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  See \ref SCIP_Status "SCIP_STATUS" for a complete list of all possible solving status.
 */
SCIP_RETCODE SCIPprintStatus(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintStatus", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( SCIPgetStatus(scip) )
   {
   case SCIP_STATUS_UNKNOWN:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "unknown");
      break;
   case SCIP_STATUS_USERINTERRUPT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "user interrupt");
      break;
   case SCIP_STATUS_NODELIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "node limit reached");
      break;
   case SCIP_STATUS_TOTALNODELIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "total node limit reached");
      break;
   case SCIP_STATUS_STALLNODELIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "stall node limit reached");
      break;
   case SCIP_STATUS_TIMELIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "time limit reached");
      break;
   case SCIP_STATUS_MEMLIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "memory limit reached");
      break;
   case SCIP_STATUS_GAPLIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "gap limit reached");
      break;
   case SCIP_STATUS_SOLLIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "solution limit reached");
      break;
   case SCIP_STATUS_BESTSOLLIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "solution improvement limit reached");
      break;
   case SCIP_STATUS_RESTARTLIMIT:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "restart limit reached");
      break;
   case SCIP_STATUS_OPTIMAL:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "optimal solution found");
      break;
   case SCIP_STATUS_INFEASIBLE:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "infeasible");
      break;
   case SCIP_STATUS_UNBOUNDED:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "unbounded");
      break;
   case SCIP_STATUS_INFORUNBD:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "infeasible or unbounded");
      break;
   case SCIP_STATUS_TERMINATE:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "termination signal received");
      break;
   default:
      SCIPerrorMessage("invalid status code <%d>\n", SCIPgetStatus(scip));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** returns whether the current stage belongs to the transformed problem space
 *
 *  @return Returns TRUE if the \SCIP instance is transformed, otherwise FALSE
 */
SCIP_Bool SCIPisTransformed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return ((int)scip->set->stage >= (int)SCIP_STAGE_TRANSFORMING);
}

/** returns whether the solution process should be probably correct
 *
 *  @note This feature is not supported yet!
 *
 *  @return Returns TRUE if \SCIP is exact solving mode, otherwise FALSE
 */
SCIP_Bool SCIPisExactSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return (scip->set->misc_exactsolve);
}

/** returns whether the presolving process would be finished given no more presolving reductions are found in this
 *  presolving round
 *
 *  Checks whether the number of presolving rounds is not exceeded and the presolving reductions found in the current
 *  presolving round suffice to trigger another presolving round.
 *
 *  @note if subsequent presolvers find more reductions, presolving might continue even if the method returns FALSE
 *  @note does not check whether infeasibility or unboundedness was already detected in presolving (which would result
 *        in presolving being stopped although the method returns TRUE)
 *
 *  @return Returns TRUE if presolving is finished if no further reductions are detected
 */
SCIP_Bool SCIPisPresolveFinished(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   int maxnrounds;
   SCIP_Bool finished;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisPresolveFinished", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* get maximum number of presolving rounds */
   maxnrounds = scip->set->presol_maxrounds;
   if( maxnrounds == -1 )
      maxnrounds = INT_MAX;

   /* don't abort, if enough changes were applied to the variables */
   finished = (scip->transprob->nvars == 0
      || (scip->stat->npresolfixedvars - scip->stat->lastnpresolfixedvars
         + scip->stat->npresolaggrvars - scip->stat->lastnpresolaggrvars
         + scip->stat->npresolchgvartypes - scip->stat->lastnpresolchgvartypes
         + (scip->stat->npresolchgbds - scip->stat->lastnpresolchgbds)/10.0
         + (scip->stat->npresoladdholes - scip->stat->lastnpresoladdholes)/10.0
         <= scip->set->presol_abortfac * scip->transprob->nvars)); /*lint !e653*/

   /* don't abort, if enough changes were applied to the constraints */
   finished = finished
      && (scip->transprob->nconss == 0
         || (scip->stat->npresoldelconss - scip->stat->lastnpresoldelconss
            + scip->stat->npresoladdconss - scip->stat->lastnpresoladdconss
            + scip->stat->npresolupgdconss - scip->stat->lastnpresolupgdconss
            + scip->stat->npresolchgsides - scip->stat->lastnpresolchgsides
            <= scip->set->presol_abortfac * scip->transprob->nconss));

   /* don't abort, if enough changes were applied to the coefficients (assume a 1% density of non-zero elements) */
   finished = finished
      && (scip->transprob->nvars == 0 || scip->transprob->nconss == 0
         || (scip->stat->npresolchgcoefs - scip->stat->lastnpresolchgcoefs
            <= scip->set->presol_abortfac * 0.01 * scip->transprob->nvars * scip->transprob->nconss));

#ifdef SCIP_DISABLED_CODE
   /* since 2005, we do not take cliques and implications into account when deciding whether to stop presolving */
   /* don't abort, if enough new implications or cliques were found (assume 100 implications per variable) */
   finished = finished
      && (scip->stat->nimplications - scip->stat->lastnpresolimplications
         <= scip->set->presol_abortfac * 100 * scip->transprob->nbinvars)
      && (SCIPcliquetableGetNCliques(scip->cliquetable) - scip->stat->lastnpresolcliques
         <= scip->set->presol_abortfac * scip->transprob->nbinvars);
#endif

   /* abort if maximal number of presolving rounds is reached */
   finished = finished || (scip->stat->npresolrounds + 1 >= maxnrounds);

   return finished;
}

/** returns whether SCIP has performed presolving during the last solve
 *
 *  @return Returns TRUE if presolving was performed during the last solve
 */
SCIP_Bool SCIPhasPerformedPresolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPhasPerformedPresolve", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->performpresol;
}

/** returns whether the user pressed CTRL-C to interrupt the solving process
 *
 *  @return Returns TRUE if Ctrl-C was pressed, otherwise FALSE.
 */ /*lint -e715*/
SCIP_Bool SCIPpressedCtrlC(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return SCIPinterrupted();
}

/** returns whether the solving process should be / was stopped before proving optimality;
 *  if the solving process should be / was stopped, the status returned by SCIPgetStatus() yields
 *  the reason for the premature abort
 *
 *  @return Returns TRUE if solving process is stopped/interrupted, otherwise FALSE.
 */
SCIP_Bool SCIPisStopped(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisStopped", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsolveIsStopped(scip->set, scip->stat, FALSE);
}

/** includes information about an external code linked into the SCIP library */
SCIP_RETCODE SCIPincludeExternalCodeInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of external code */
   const char*           description         /**< description of external code, or NULL */
   )
{
   assert(scip != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeExternalCodeInformation", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsetIncludeExternalCode(scip->set, name, description) );

   return SCIP_OKAY;
}

/** returns an array of names of currently included external codes */
char** SCIPgetExternalCodeNames(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->extcodenames;
}

/** returns an array of the descriptions of currently included external codes
 *
 *  @note some descriptions may be NULL
 */
char** SCIPgetExternalCodeDescriptions(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->extcodedescs;
}

/** returns the number of currently included information on external codes */
int SCIPgetNExternalCodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nextcodes;
}

/** prints information on external libraries to a file stream via the message handler system
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
void SCIPprintExternalCodes(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "External libraries: ");
   if( scip->set->nextcodes == 0 )
   {
      SCIPinfoMessage(scip, file, "none\n");
      return;
   }
   SCIPinfoMessage(scip, file, "\n");

   for( i = 0; i < scip->set->nextcodes; ++i )
   {
      SCIPinfoMessage(scip, file, "  %-20s %s\n", scip->set->extcodenames[i], scip->set->extcodedescs[i] != NULL ? scip->set->extcodedescs[i] : "");
   }
}
