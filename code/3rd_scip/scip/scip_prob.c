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

/**@file   scip_prob.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for global and local (sub)problems
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
#include "scip/benders.h"
#include "scip/clock.h"
#include "scip/concurrent.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/dcmp.h"
#include "scip/debug.h"
#include "scip/lp.h"
#include "scip/pricer.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/pub_cons.h"
#include "scip/pub_event.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/pub_sol.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/reader.h"
#include "scip/reopt.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/struct_cons.h"
#include "scip/struct_lp.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/syncstore.h"
#include "scip/tree.h"
#include <stdio.h>
#include <string.h>

/** creates empty problem and initializes all solving data structures (the objective sense is set to MINIMIZE)
 *  If the problem type requires the use of variable pricers, these pricers should be added to the problem with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method, \SCIP reaches the following stage:
 *        - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPcreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   SCIP_DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   SCIP_DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   SCIP_DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   SCIP_DECL_PROBINITSOL ((*probinitsol)),   /**< solving process initialization method of transformed data */
   SCIP_DECL_PROBEXITSOL ((*probexitsol)),   /**< solving process deinitialization method of transformed data */
   SCIP_DECL_PROBCOPY    ((*probcopy)),      /**< copies user data if you want to copy it to a subscip, or NULL */
   SCIP_PROBDATA*        probdata            /**< user problem data set by the reader */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateProb", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE) );

   /* free old problem */
   SCIP_CALL( SCIPfreeProb(scip) );
   assert(scip->set->stage == SCIP_STAGE_INIT);

   /* switch stage to PROBLEM */
   scip->set->stage = SCIP_STAGE_PROBLEM;

   SCIP_CALL( SCIPstatCreate(&scip->stat, scip->mem->probmem, scip->set, NULL, NULL, scip->messagehdlr) );

   SCIP_CALL( SCIPprobCreate(&scip->origprob, scip->mem->probmem, scip->set, name,
         probdelorig, probtrans, probdeltrans, probinitsol, probexitsol, probcopy, probdata, FALSE) );

   /* create solution pool for original solution candidates */
   SCIP_CALL( SCIPprimalCreate(&scip->origprimal) );

   /* create conflict pool for storing conflict constraints */
   SCIP_CALL( SCIPconflictstoreCreate(&scip->conflictstore, scip->set) );

   /* initialize reoptimization structure, if needed */
   SCIP_CALL( SCIPenableReoptimization(scip, scip->set->reopt_enable) );

   SCIP_CALL( SCIPdecompstoreCreate(&scip->decompstore, SCIPblkmem(scip), SCIP_DECOMPSTORE_CAPA) );

   return SCIP_OKAY;
}

/** creates empty problem and initializes all solving data structures (the objective sense is set to MINIMIZE)
 *  all callback methods will be set to NULL and can be set afterwards, if needed, via SCIPsetProbDelorig(),
 *  SCIPsetProbTrans(), SCIPsetProbDeltrans(), SCIPsetProbInitsol(), SCIPsetProbExitsol(), and
 *  SCIPsetProbCopy()
 *  If the problem type requires the use of variable pricers, these pricers should be added to the problem with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method, \SCIP reaches the following stage:
 *        - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPcreateProbBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< problem name */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateProbBasic", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPcreateProb(scip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** sets callback to free user data of original problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetProbDelorig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_PROBDELORIG ((*probdelorig))    /**< frees user data of original problem */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbDelorig", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPprobSetDelorig(scip->origprob, probdelorig);

   return SCIP_OKAY;
}

/** sets callback to create user data of transformed problem by transforming original user data
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetProbTrans(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_PROBTRANS   ((*probtrans))      /**< creates user data of transformed problem by transforming original user data */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbTrans", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPprobSetTrans(scip->origprob, probtrans);

   return SCIP_OKAY;
}

/** sets callback to free user data of transformed problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetProbDeltrans(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_PROBDELTRANS((*probdeltrans))   /**< frees user data of transformed problem */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbDeltrans", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPprobSetDeltrans(scip->origprob, probdeltrans);

   return SCIP_OKAY;
}

/** sets solving process initialization callback of transformed data
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetProbInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_PROBINITSOL ((*probinitsol))    /**< solving process initialization method of transformed data */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbInitsol", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPprobSetInitsol(scip->origprob, probinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization callback of transformed data
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetProbExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_PROBEXITSOL ((*probexitsol))    /**< solving process deinitialization method of transformed data */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbExitsol", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPprobSetExitsol(scip->origprob, probexitsol);

   return SCIP_OKAY;
}

/** sets callback to copy user data to a subscip
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetProbCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_PROBCOPY    ((*probcopy))       /**< copies user data if you want to copy it to a subscip, or NULL */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbCopy", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPprobSetCopy(scip->origprob, probcopy);

   return SCIP_OKAY;
}

/** reads problem from file and initializes all solving data structures
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @post After the method was called, \SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT if reading failed (usually, when a SCIP_READERROR occurs)
 *       - \ref SCIP_STAGE_PROBLEM if the problem file was successfully read
 */
SCIP_RETCODE SCIPreadProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< problem file name */
   const char*           extension           /**< extension of the desired file reader,
                                              *   or NULL if file extension should be used */
   )
{
   SCIP_RETCODE retcode;
   SCIP_RESULT result;
   SCIP_Bool usevartable;
   SCIP_Bool useconstable;
   int i;
   char* tmpfilename;
   char* fileextension;

   assert(scip != NULL);
   assert(filename != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPreadProb", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPgetBoolParam(scip, "misc/usevartable", &usevartable) );
   SCIP_CALL( SCIPgetBoolParam(scip, "misc/useconstable", &useconstable) );

   if( !usevartable || !useconstable )
   {
      SCIPerrorMessage("Cannot read problem if vartable or constable is disabled. Make sure parameters 'misc/usevartable' and 'misc/useconstable' are set to TRUE.\n");
      return SCIP_READERROR;
   }

   /* try all readers until one could read the file */
   result = SCIP_DIDNOTRUN;

   /* copy filename */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpfilename, filename, (int)strlen(filename)+1) );

   fileextension = NULL;
   if( extension == NULL )
   {
      /* get extension from filename */
      SCIPsplitFilename(tmpfilename, NULL, NULL, &fileextension, NULL);
   }

   for( i = 0; i < scip->set->nreaders && result == SCIP_DIDNOTRUN; ++i )
   {
      retcode = SCIPreaderRead(scip->set->readers[i], scip->set, filename,
            extension != NULL ? extension : fileextension, &result);

      /* check for reader errors */
      if( retcode == SCIP_NOFILE || retcode == SCIP_READERROR )
         goto TERMINATE;
      SCIP_CALL( retcode );
   }

   switch( result )
   {
   case SCIP_DIDNOTRUN:
      retcode = SCIP_PLUGINNOTFOUND;
      break;
   case SCIP_SUCCESS:
      if( scip->origprob != NULL )
      {
         SCIP_Real readingtime;

         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
            "original problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
            scip->origprob->nvars, scip->origprob->nbinvars, scip->origprob->nintvars,
            scip->origprob->nimplvars, scip->origprob->ncontvars,
            scip->origprob->nconss);

         /* in full verbose mode we will also print the number of constraints per constraint handler */
         if( scip->set->disp_verblevel == SCIP_VERBLEVEL_FULL )
         {
            int* nconss;
            int c;
            int h;

            SCIP_CALL( SCIPallocClearBufferArray(scip, &nconss, scip->set->nconshdlrs) );

            /* loop over all constraints and constraint-handlers to count for each type the amount of original
             * constraints
             */
            for( c = scip->origprob->nconss - 1; c >= 0; --c )
            {
               for( h = scip->set->nconshdlrs - 1; h >= 0; --h )
               {
                  if( scip->origprob->conss[c]->conshdlr == scip->set->conshdlrs[h] )
                  {
                     ++(nconss[h]);
                     break;
                  }
               }
               /* constraint handler should be found */
               assert(h >= 0);
            }

            /* loop over all constraints handlers for printing the number of original constraints */
            for( h = 0; h < scip->set->nconshdlrs; ++h )
            {
               if( nconss[h] > 0 )
               {
                  SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                     "%7d constraints of type <%s>\n", nconss[h], SCIPconshdlrGetName(scip->set->conshdlrs[h]));
               }
            }

            SCIPfreeBufferArray(scip, &nconss);
         }

         /* in case the permutation seed is different to 0, permute the original problem */
         if( scip->set->random_permutationseed > 0 )
         {
            SCIP_Bool permuteconss;
            SCIP_Bool permutevars;
            int permutationseed;

            permuteconss = scip->set->random_permuteconss;
            permutevars = scip->set->random_permutevars;
            permutationseed = scip->set->random_permutationseed;

            SCIP_CALL( SCIPpermuteProb(scip, (unsigned int)permutationseed, permuteconss, permutevars, permutevars, permutevars, permutevars) );
         }

         /* get reading time */
         readingtime = SCIPgetReadingTime(scip);

         /* display timing statistics */
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "Reading Time: %.2f\n", readingtime);
      }
      retcode = SCIP_OKAY;
      break;
   default:
      assert(i < scip->set->nreaders);
      SCIPerrorMessage("invalid result code <%d> from reader <%s> reading file <%s>\n",
         result, SCIPreaderGetName(scip->set->readers[i]), filename);
      retcode = SCIP_READERROR;
   }  /*lint !e788*/

 TERMINATE:
   /* free buffer array */
   SCIPfreeBufferArray(scip, &tmpfilename);

   /* check if reading time should belong to solving time */
   if( scip->set->time_reading )
   {
      SCIP_Real readingtime;

      /* get reading time */
      readingtime = SCIPgetReadingTime(scip);

      /* add reading time to solving time */
      SCIPclockSetTime(scip->stat->solvingtime, readingtime);
   }

   return retcode;
}

/** write original or transformed problem */
static
SCIP_RETCODE writeProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< output file (or NULL for standard output) */
   const char*           extension,          /**< extension of the desired file reader,
                                              *   or NULL if file extension should be used */
   SCIP_Bool             transformed,        /**< output the transformed problem? */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;
   char* tmpfilename;
   char* fileextension;
   char* compression;
   FILE* file;

   assert(scip != NULL );

   fileextension = NULL;
   compression = NULL;
   file = NULL;
   tmpfilename = NULL;

   if( filename != NULL &&  filename[0] != '\0' )
   {
      int success;

      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPerrorMessage("cannot create file <%s> for writing\n", filename);
         SCIPprintSysError(filename);
         return SCIP_FILECREATEERROR;
      }

      /* get extension from filename,
       * if an error occurred, close the file before returning */
      if( BMSduplicateMemoryArray(&tmpfilename, filename, strlen(filename)+1) == NULL )
      {
         (void) fclose(file);
         SCIPerrorMessage("Error <%d> in function call\n", SCIP_NOMEMORY);
         return SCIP_NOMEMORY;
      }

      SCIPsplitFilename(tmpfilename, NULL, NULL, &fileextension, &compression);

      if( compression != NULL )
      {
         SCIPmessagePrintWarning(scip->messagehdlr, "currently it is not possible to write files with any compression\n");
         BMSfreeMemoryArray(&tmpfilename);
         (void) fclose(file);
         return SCIP_FILECREATEERROR;
      }

      if( extension == NULL && fileextension == NULL )
      {
         SCIPmessagePrintWarning(scip->messagehdlr, "filename <%s> has no file extension, select default <cip> format for writing\n", filename);
      }

      if( transformed )
         retcode = SCIPprintTransProblem(scip, file, extension != NULL ? extension : fileextension, genericnames);
      else
         retcode = SCIPprintOrigProblem(scip, file, extension != NULL ? extension : fileextension, genericnames);

      BMSfreeMemoryArray(&tmpfilename);

      success = fclose(file);
      if( success != 0 )
      {
         SCIPerrorMessage("An error occurred while closing file <%s>\n", filename);
         return SCIP_FILECREATEERROR;
      }
   }
   else
   {
      /* print to stdout */
      if( transformed )
         retcode = SCIPprintTransProblem(scip, NULL, extension, genericnames);
      else
         retcode = SCIPprintOrigProblem(scip, NULL, extension, genericnames);
   }

   /* check for write errors */
   if( retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** writes original problem to file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPwriteOrigProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< output file (or NULL for standard output) */
   const char*           extension,          /**< extension of the desired file reader,
                                              *   or NULL if file extension should be used */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteOrigProblem", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert( scip != NULL );
   assert( scip->origprob != NULL );

   retcode = writeProblem(scip, filename, extension, FALSE, genericnames);

   /* check for write errors */
   if( retcode == SCIP_FILECREATEERROR || retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** writes transformed problem which are valid in the current node to file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note If you want the write all constraints (including the once which are redundant for example), you need to set
 *        the parameter <write/allconss> to TRUE
 */
SCIP_RETCODE SCIPwriteTransProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< output file (or NULL for standard output) */
   const char*           extension,          /**< extension of the desired file reader,
                                              *   or NULL if file extension should be used */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteTransProblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert( scip != NULL );
   assert( scip->transprob != NULL );

   retcode = writeProblem(scip, filename, extension, TRUE, genericnames);

   /* check for write errors */
   if( retcode == SCIP_FILECREATEERROR || retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** frees problem and solution process data
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After this method was called, SCIP is in the following stage:
 *       - \ref SCIP_STAGE_INIT
 */
SCIP_RETCODE SCIPfreeProb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool transsolorig;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeProb", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE) );

   /* if we free the problem, we do not have to transfer transformed solutions to the original space, so temporarily disable it */
   transsolorig = scip->set->misc_transsolsorig;
   scip->set->misc_transsolsorig = FALSE;

   SCIP_CALL( SCIPfreeTransform(scip) );
   /* for some reason the free transform can generate events caught in the globalbnd eventhander
    * which requires the concurrent so it must be freed afterwards this happened o instance fiber
    */
   SCIP_CALL( SCIPfreeConcurrent(scip) );

   assert(scip->set->stage == SCIP_STAGE_INIT || scip->set->stage == SCIP_STAGE_PROBLEM);
   scip->set->misc_transsolsorig = transsolorig;

   if( scip->set->stage == SCIP_STAGE_PROBLEM )
   {
      int i;

      /* free concsolvers and deinitialize the syncstore */
      if( scip->set->nconcsolvers > 0 )
      {
         assert(SCIPsyncstoreIsInitialized(scip->syncstore));

         SCIP_CALL( SCIPsetFreeConcsolvers(scip->set) );
         SCIP_CALL( SCIPsyncstoreExit(scip->syncstore) );
      }

      /* deactivate all pricers */
      for( i = scip->set->nactivepricers-1; i >= 0; --i )
      {
         SCIP_CALL( SCIPpricerDeactivate(scip->set->pricers[i], scip->set) );
      }
      assert(scip->set->nactivepricers == 0);

      /* deactivate all Benders' decomposition */
      for( i = scip->set->nactivebenders-1; i >= 0; --i )
      {
         SCIP_CALL( SCIPbendersDeactivate(scip->set->benders[i], scip->set) );
      }
      assert(scip->set->nactivebenders == 0);

      /* free all debug data */
      SCIP_CALL( SCIPdebugFreeDebugData(scip->set) );

      /* free original primal solution candidate pool, original problem and problem statistics data structures */
      if( scip->reopt != NULL )
      {
         SCIP_CALL( SCIPreoptFree(&scip->reopt, scip->set, scip->origprimal, scip->mem->probmem) );
      }
      SCIPdecompstoreFree(&scip->decompstore, SCIPblkmem(scip));
      SCIP_CALL( SCIPconflictstoreFree(&scip->conflictstore, scip->mem->probmem, scip->set, scip->stat, scip->reopt) );
      SCIP_CALL( SCIPprimalFree(&scip->origprimal, scip->mem->probmem) );
      SCIP_CALL( SCIPprobFree(&scip->origprob, scip->messagehdlr, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
      SCIP_CALL( SCIPstatFree(&scip->stat, scip->mem->probmem) );

      /* readers */
      for( i = 0; i < scip->set->nreaders; ++i )
      {
         SCIP_CALL( SCIPreaderResetReadingTime(scip->set->readers[i]) );
      }

      /* switch stage to INIT */
      scip->set->stage = SCIP_STAGE_INIT;
   }
   assert(scip->set->stage == SCIP_STAGE_INIT);

   return SCIP_OKAY;
}

/** permutes parts of the problem data structure
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *
 *  @todo This need to be changed to use the new random number generator implemented in random.c
 */
SCIP_RETCODE SCIPpermuteProb(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          randseed,           /**< seed value for random generator */
   SCIP_Bool             permuteconss,       /**< should the list of constraints in each constraint handler be permuted? */
   SCIP_Bool             permutebinvars,     /**< should the list of binary variables be permuted? */
   SCIP_Bool             permuteintvars,     /**< should the list of integer variables be permuted? */
   SCIP_Bool             permuteimplvars,    /**< should the list of implicit integer variables be permuted? */
   SCIP_Bool             permutecontvars     /**< should the list of continuous integer variables be permuted? */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSHDLR** conshdlrs;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Bool permuted;
   int nconshdlrs;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int nvars;
   int j;

   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPpermuteProb", FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, NULL) );

   assert(nvars == 0 || vars != NULL);
   assert(nvars == nbinvars+nintvars+nimplvars+SCIPgetNContVars(scip));

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);
   assert(nconshdlrs == 0 || conshdlrs != NULL);

   /* create a random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, randseed, TRUE) );

   /* The constraint handler should not be permuted since they are called w.r.t. to certain properties; besides
    * that the "conshdlrs" array should stay in the order as it is since this array is used to copy the plugins for
    * sub-SCIPs and contains the dependencies between the constraint handlers; for example the linear constraint
    * handler stays in front of all constraint handler which can upgrade a linear constraint (such as logicor,
    * setppc, and knapsack).
    */

   permuted = FALSE;

   /* for each constraint handler, permute its constraints */
   if( permuteconss )
   {
      int i;

      /* we must only permute active constraints */
      if( SCIPisTransformed(scip) && !SCIPprobIsPermuted(scip->transprob) )
      {
         /* loop over all constraint handlers */
         for( i = 0; i < nconshdlrs; ++i )
         {
            SCIP_CONS** conss;
            int nconss;

            conss = SCIPconshdlrGetConss(conshdlrs[i]);
            nconss = SCIPconshdlrGetNActiveConss(conshdlrs[i]);

            assert(nconss == 0 || conss != NULL);

            SCIPrandomPermuteArray(randnumgen, (void**)conss, 0, nconss);

            /* readjust the mapping of constraints to array positions */
            for( j = 0; j < nconss; ++j )
               conss[j]->consspos = j;

            permuted = TRUE;
         }
      }
      else if( !SCIPisTransformed(scip) && !SCIPprobIsPermuted(scip->origprob) )
      {
         SCIP_CONS** conss = scip->origprob->conss;
         int nconss = scip->origprob->nconss;

         SCIPrandomPermuteArray(randnumgen, (void**)conss, 0, nconss);

         for( j = 0; j < nconss; ++j )
         {
            assert(conss[j]->consspos == -1);
            conss[j]->addarraypos = j;
         }

         permuted = TRUE;
      }
   }

   /* permute binary variables */
   if( permutebinvars && !SCIPprobIsPermuted(scip->origprob) )
   {
      SCIPrandomPermuteArray(randnumgen, (void**)vars, 0, nbinvars);

      /* readjust the mapping of variables to array positions */
      for( j = 0; j < nbinvars; ++j )
         vars[j]->probindex = j;

      permuted = TRUE;
   }

   /* permute general integer variables */
   if( permuteintvars && !SCIPprobIsPermuted(scip->origprob) )
   {
      SCIPrandomPermuteArray(randnumgen, (void**)vars, nbinvars, nbinvars+nintvars);

      /* readjust the mapping of variables to array positions */
      for( j = nbinvars; j < nbinvars+nintvars; ++j )
         vars[j]->probindex = j;

      permuted = TRUE;
   }

   /* permute general integer variables */
   if( permuteimplvars && !SCIPprobIsPermuted(scip->origprob) )
   {
      SCIPrandomPermuteArray(randnumgen, (void**)vars, nbinvars+nintvars, nbinvars+nintvars+nimplvars);

      /* readjust the mapping of variables to array positions */
      for( j = nbinvars+nintvars; j < nbinvars+nintvars+nimplvars; ++j )
         vars[j]->probindex = j;

      permuted = TRUE;
   }

   /* permute general integer variables */
   if( permutecontvars && !SCIPprobIsPermuted(scip->origprob) )
   {
      SCIPrandomPermuteArray(randnumgen, (void**)vars, nbinvars+nintvars+nimplvars, nvars);

      /* readjust the mapping of variables to array positions */
      for( j = nbinvars+nintvars+nimplvars; j < nvars; ++j )
         vars[j]->probindex = j;

      permuted = TRUE;
   }

   if( permuted && SCIPisTransformed(scip) )
   {
      assert(!SCIPprobIsPermuted(scip->transprob));

      /* mark tranformed problem as permuted */
      SCIPprobMarkPermuted(scip->transprob);

      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "permute transformed problem using random seed %u\n", randseed);
   }
   else if( permuted && !SCIPisTransformed(scip) )
   {
      assert(!SCIPprobIsPermuted(scip->origprob));

      /* mark original problem as permuted */
      SCIPprobMarkPermuted(scip->origprob);

      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "permute original problem using random seed %u\n", randseed);
   }

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}

/** gets user problem data
 *
 *  @return a SCIP_PROBDATA pointer, or NULL if no problem data was allocated
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_PROBDATA* SCIPgetProbData(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetProbData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobGetData(scip->origprob);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      return SCIPprobGetData(scip->transprob);

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** sets user problem data
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPsetProbData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< user problem data to use */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetData(scip->origprob, probdata);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      SCIPprobSetData(scip->transprob, probdata);
      return SCIP_OKAY;

   case SCIP_STAGE_INIT:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }
}

/** returns name of the current problem instance
 *
 *  @return name of the current problem instance
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
const char* SCIPgetProbName(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetProbName", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprobGetName(scip->origprob);
}

/** sets name of the current problem instance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPsetProbName(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name to be set */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbName", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprobSetName(scip->origprob, name);
}

/** changes the objective function of the original problem.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVED
 *
 *  @note This method should be only used to change the objective function during two reoptimization runs and is only
 *        recommended to an experienced user.
 *
 *  @note All variables not given in \p vars array are assumed to have an objective coefficient of zero.
 */
SCIP_RETCODE SCIPchgReoptObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJSENSE         objsense,           /**< new objective function */
   SCIP_VAR**            vars,               /**< original problem variables */
   SCIP_Real*            coefs,              /**< objective coefficients */
   int                   nvars               /**< variables in vars array */
   )
{
   SCIP_VAR** origvars;
   int norigvars;
   int i;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgReoptObjective", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || coefs != NULL);

   origvars = scip->origprob->vars;
   norigvars = scip->origprob->nvars;

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(scip, "objective function need to be set:\n");
   for( i = 0; i < nvars; i++ )
   {
      if( !SCIPisZero(scip, coefs[i]) )
         SCIPdebugMsg(scip, "%s%g <%s> ", SCIPisPositive(scip, coefs[i]) ? "+" : "", coefs[i], SCIPvarGetName(vars[i]));
   }
   SCIPdebugMsg(scip, "\n");
#endif

   /* Set all coefficients of original variables to 0, since we will add the new objective coefficients later. */
   for( i = 0; i < norigvars; i++ )
   {
      SCIP_CALL( SCIPchgVarObj(scip, origvars[i], 0.0) );
   }

   if( scip->set->stage >= SCIP_STAGE_TRANSFORMED )
   {
      /* In order to avoid numerical troubles, also explicitly set all transformed objective coefficients to 0. */
      for( i = 0; i < scip->transprob->nvars; i++ )
      {
         SCIP_CALL( SCIPchgVarObj(scip, scip->transprob->vars[i], 0.0) );
      }
   }

   /* reset objective data of original problem */
   scip->origprob->objscale = 1.0;
   scip->origprob->objsense = objsense;
   scip->origprob->objoffset = 0.0;
   scip->origprob->objisintegral = FALSE;

   if( scip->set->stage >= SCIP_STAGE_TRANSFORMED )
   {
      /* reset objective data of transformed problem */
      scip->transprob->objscale = 1.0;
      scip->transprob->objsense = objsense;
      scip->transprob->objoffset = 0.0;
      scip->transprob->objisintegral = FALSE;
   }

   /* set new objective values */
   for( i = 0; i < nvars; ++i )
   {
      if( !SCIPvarIsOriginal(vars[i]) )
      {
         SCIPerrorMessage("Cannot handle variable <%s> (status: %d) in SCIPchgReoptObjective().\n",
               SCIPvarGetName(vars[i]), SCIPvarGetStatus(vars[i]));
         return SCIP_INVALIDDATA;
      }

      /* Add coefficients because this gets transferred to the transformed problem (the coefficients were set to 0 above). */
      SCIP_CALL( SCIPaddVarObj(scip, vars[i], coefs[i]) );
   }

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(scip, "new objective function:\n");
   for( i = 0; i < norigvars; i++ )
   {
      SCIP_Real objval = SCIPvarGetObj(origvars[i]);
      if( !SCIPisZero(scip, objval) )
         SCIPdebugMsg(scip, "%s%g <%s> ", SCIPisPositive(scip, objval) ? "+" : "", objval, SCIPvarGetName(origvars[i]));
   }
   SCIPdebugMsg(scip, "\n");
#endif

   return SCIP_OKAY;
}

/** returns objective sense of original problem
 *
 *  @return objective sense of original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_OBJSENSE SCIPgetObjsense(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetObjsense", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->objsense;
}

/** sets objective sense of problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetObjsense(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJSENSE         objsense            /**< new objective sense */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetObjsense", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   if( objsense != SCIP_OBJSENSE_MAXIMIZE && objsense != SCIP_OBJSENSE_MINIMIZE )
   {
      SCIPerrorMessage("invalid objective sense\n");
      return SCIP_INVALIDDATA;
   }

   SCIPprobSetObjsense(scip->origprob, objsense);

   return SCIP_OKAY;
}

/** adds offset of objective function
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 */
SCIP_RETCODE SCIPaddObjoffset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             addval              /**< value to add to objective offset */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddObjoffset", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPprobAddObjoffset(scip->transprob, addval);
   SCIP_CALL( SCIPprimalUpdateObjoffset(scip->primal, SCIPblkmem(scip), scip->set, scip->stat, scip->eventfilter,
         scip->eventqueue, scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp) );

   return SCIP_OKAY;
}

/** adds offset of objective function to original problem and to all existing solution in original space
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPaddOrigObjoffset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             addval              /**< value to add to objective offset */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddOrigObjoffset", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   scip->origprob->objoffset += addval;
   SCIPprimalAddOrigObjoffset(scip->origprimal, scip->set, addval);

   return SCIP_OKAY;
}

/** returns the objective offset of the original problem
 *
 *  @return the objective offset of the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetOrigObjoffset(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetOrigObjoffset", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->origprob->objoffset;
}

/** returns the objective scale of the original problem
 *
 *  @return the objective scale of the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetOrigObjscale(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetOrigObjscale", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->origprob->objscale;
}

/** returns the objective offset of the transformed problem
 *
 *  @return the objective offset of the transformed problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetTransObjoffset(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetTransObjoffset", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->transprob->objoffset;
}

/** returns the objective scale of the transformed problem
 *
 *  @return the objective scale of the transformed problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetTransObjscale(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetTransObjscale", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->transprob->objscale;
}

/** sets limit on objective function, such that only solutions better than this limit are accepted
 *
 *  @note SCIP will only look for solutions with a strictly better objective value, thus, e.g., prune
 *        all branch-and-bound nodes with dual bound equal or worse to the objective limit.
 *        However, SCIP will also collect solutions with objective value worse than the objective limit and
 *        use them to run improvement heuristics on them.
 *  @note If SCIP can prove that there exists no solution with a strictly better objective value, the solving status
 *        will normally be infeasible (the objective limit is interpreted as part of the problem).
 *        The only exception is that by chance, SCIP found a solution with the same objective value and thus
 *        proved the optimality of this solution, resulting in solution status optimal.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetObjlimit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             objlimit            /**< new primal objective limit */
   )
{
   SCIP_Real oldobjlimit;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetObjlimit", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetObjlim(scip->origprob, objlimit);
      break;
   case SCIP_STAGE_PRESOLVED:
      oldobjlimit = SCIPprobGetObjlim(scip->origprob, scip->set);
      assert(oldobjlimit == SCIPprobGetObjlim(scip->transprob, scip->set)); /*lint !e777*/
      if( SCIPtransformObj(scip, objlimit) > SCIPprobInternObjval(scip->transprob, scip->origprob, scip->set, oldobjlimit) && ! scip->set->reopt_enable)
      {
         SCIPerrorMessage("cannot relax objective limit from %.15g to %.15g in presolved stage.\n", oldobjlimit, objlimit);
         return SCIP_INVALIDDATA;
      }
      SCIPprobSetObjlim(scip->origprob, objlimit);
      SCIPprobSetObjlim(scip->transprob, objlimit);
      SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
            scip->eventqueue, scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp) );
      break;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_SOLVING:
      oldobjlimit = SCIPprobGetObjlim(scip->origprob, scip->set);
      assert(oldobjlimit == SCIPprobGetObjlim(scip->transprob, scip->set)); /*lint !e777*/
      if( SCIPtransformObj(scip, objlimit) > SCIPprobInternObjval(scip->transprob, scip->origprob, scip->set, oldobjlimit) )
      {
         SCIPerrorMessage("cannot relax objective limit from %.15g to %.15g after problem was transformed.\n", oldobjlimit, objlimit);
         return SCIP_INVALIDDATA;
      }
      SCIPprobSetObjlim(scip->origprob, objlimit);
      SCIPprobSetObjlim(scip->transprob, objlimit);
      SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
            scip->eventqueue, scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp) );
      break;

   default:
      SCIPerrorMessage("method is not callable in SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   } /*lint !e788*/

   return SCIP_OKAY;
}

/** returns current limit on objective function
 *
 *  @return the current objective limit of the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetObjlimit(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetObjlimit", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobGetObjlim(scip->origprob, scip->set);
}

/** informs SCIP, that the objective value is always integral in every feasible solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This function should be used to inform SCIP that the objective function is integral, helping to improve the
 *        performance. This is useful when using column generation. If no column generation (pricing) is used, SCIP
 *        automatically detects whether the objective function is integral or can be scaled to be integral. However, in
 *        any case, the user has to make sure that no variable is added during the solving process that destroys this
 *        property.
 */
SCIP_RETCODE SCIPsetObjIntegral(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetObjIntegral", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetObjIntegral(scip->origprob);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIPprobSetObjIntegral(scip->transprob);
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("method is not callable in SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   } /*lint !e788*/
}

/** returns whether the objective value is known to be integral in every feasible solution
 *
 *  @return TRUE, if objective value is known to be always integral, otherwise FALSE
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note If no pricing is performed, SCIP automatically detects whether the objective function is integral or can be
 *        scaled to be integral, helping to improve performance. This function returns the result. Otherwise
 *        SCIPsetObjIntegral() can be used to inform SCIP. However, in any case, the user has to make sure that no
 *        variable is added during the solving process that destroys this property.
 */
SCIP_Bool SCIPisObjIntegral(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   int v;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisObjIntegral", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* if the user explicitly added the information that there is an integral objective, return TRUE */
      if( SCIPprobIsObjIntegral(scip->origprob) )
         return TRUE;

      /* if there exist unknown variables, we cannot conclude that the objective value is always integral */
      if ( scip->set->nactivepricers != 0 )
         return FALSE;

      /* if the objective value offset is fractional, the value itself is possibly fractional */
      if ( ! SCIPisIntegral(scip, scip->origprob->objoffset) )
         return FALSE;

      /* scan through the variables */
      for (v = 0; v < scip->origprob->nvars; ++v)
      {
         SCIP_Real obj;

         /* get objective value of variable */
         obj = SCIPvarGetObj(scip->origprob->vars[v]);

         /* check, if objective value is non-zero */
         if ( ! SCIPisZero(scip, obj) )
         {
            /* if variable's objective value is fractional, the problem's objective value may also be fractional */
            if ( ! SCIPisIntegral(scip, obj) )
               break;

            /* if variable with non-zero objective value is continuous, the problem's objective value may be fractional */
            if ( SCIPvarGetType(scip->origprob->vars[v]) == SCIP_VARTYPE_CONTINUOUS )
               break;
         }
      }

      /* we do not store the result, since we assume that the original problem might be changed */
      if ( v == scip->origprob->nvars )
         return TRUE;
      return FALSE;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      return SCIPprobIsObjIntegral(scip->transprob);

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return FALSE; /*lint !e527*/
   } /*lint !e788*/
}

/** returns the Euclidean norm of the objective function vector (available only for transformed problem)
 *
 *  @return the Euclidean norm of the transformed objective function vector
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Real SCIPgetObjNorm(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetObjNorm", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( scip->lp->objsqrnormunreliable )
      SCIPlpRecalculateObjSqrNorm(scip->set, scip->lp);
   assert(!scip->lp->objsqrnormunreliable);

   return SCIPlpGetObjNorm(scip->lp);
}

/** adds variable to the problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to add */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* avoid inserting the same variable twice */
   if( SCIPvarGetProbindex(var) != -1 )
      return SCIP_OKAY;

   /* insert the negation variable x instead of the negated variable x' in x' = offset - x */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert(SCIPvarGetNegationVar(var) != NULL);
      SCIP_CALL( SCIPaddVar(scip, SCIPvarGetNegationVar(var)) );
      return SCIP_OKAY;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot add transformed variables to original problem\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPprobAddVar(scip->origprob, scip->mem->probmem, scip->set, scip->lp, scip->branchcand,
            scip->eventfilter, scip->eventqueue, var) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      /* check variable's status */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot add original variables to transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      else if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot add fixed or aggregated variables to transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPprobAddVar(scip->transprob, scip->mem->probmem, scip->set, scip->lp,
            scip->branchcand, scip->eventfilter, scip->eventqueue, var) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** adds variable to the problem and uses it as pricing candidate to enter the LP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddPricedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             score               /**< pricing score of variable (the larger, the better the variable) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddPricedVar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* insert the negation variable x instead of the negated variable x' in x' = offset - x */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert(SCIPvarGetNegationVar(var) != NULL);
      SCIP_CALL( SCIPaddPricedVar(scip, SCIPvarGetNegationVar(var), score) );
      return SCIP_OKAY;
   }

   /* add variable to problem if not yet inserted */
   if( SCIPvarGetProbindex(var) == -1 )
   {
      /* check variable's status */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot add original variables to transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      else if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot add fixed or aggregated variables to transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPprobAddVar(scip->transprob, scip->mem->probmem, scip->set, scip->lp,
            scip->branchcand, scip->eventfilter, scip->eventqueue, var) );
   }

   /* add variable to pricing storage */
   SCIP_CALL( SCIPpricestoreAddVar(scip->pricestore, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, var, score,
         (SCIPtreeGetCurrentDepth(scip->tree) == 0)) );

   return SCIP_OKAY;
}

/** removes variable from the problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @warning The variable is not deleted from the constraints when in SCIP_STAGE_PROBLEM.  In this stage, it is the
 *           user's responsibility to ensure the variable has been removed from all constraints or the constraints
 *           deleted.
 */
SCIP_RETCODE SCIPdelVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to delete */
   SCIP_Bool*            deleted             /**< pointer to store whether marking variable to be deleted was successful */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(deleted != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPdelVar", FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot remove transformed variables from original problem\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPprobDelVar(scip->origprob, scip->mem->probmem, scip->set, scip->eventqueue, var, deleted) );

      /* delete the variables from the problems that were marked to be deleted */
      SCIP_CALL( SCIPprobPerformVarDeletions(scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->cliquetable, scip->lp, scip->branchcand) );

      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      /* check variable's status */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot remove original variables from transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      else if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot remove fixed or aggregated variables from transformed problem\n");
         return SCIP_INVALIDDATA;
      }

      SCIP_CALL( SCIPprobDelVar(scip->transprob, scip->mem->probmem, scip->set, scip->eventqueue, var, deleted) );

      return SCIP_OKAY;
   case SCIP_STAGE_FREETRANS:
      /* in FREETRANS stage, we don't need to remove the variable, because the transformed problem is freed anyways */
      *deleted = FALSE;

      return SCIP_OKAY;
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** gets variables of the problem along with the numbers of different variable types; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note Variables in the vars array are ordered: binaries first, then integers, implicit integers and continuous last.
 */
SCIP_RETCODE SCIPgetVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarsData", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( vars != NULL )
         *vars = scip->origprob->vars;
      if( nvars != NULL )
         *nvars = scip->origprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->origprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->origprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->origprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->origprob->ncontvars;
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      if( vars != NULL )
         *vars = scip->transprob->vars;
      if( nvars != NULL )
         *nvars = scip->transprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->transprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->transprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->transprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->transprob->ncontvars;
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** gets array with active problem variables
 *
 *  @return array with active problem variables
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note Variables in the array are ordered: binaries first, then integers, implicit integers and continuous last.
 *
 *  @warning If your are using the methods which add or change bound of variables (e.g., SCIPchgVarType(), SCIPfixVar(),
 *           SCIPaggregateVars(), and SCIPmultiaggregateVar()), it can happen that the internal variable array (which is
 *           accessed via this method) gets resized and/or resorted. This can invalid the data pointer which is returned
 *           by this method.
 */
SCIP_VAR** SCIPgetVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->vars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      return scip->transprob->vars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of active problem variables
 *
 *  @return the number of active problem variables
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      return scip->transprob->nvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of binary active problem variables
 *
 *  @return the number of binary active problem variables
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNBinVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNBinVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nbinvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      return scip->transprob->nbinvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of integer active problem variables
 *
 *  @return the number of integer active problem variables
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNIntVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNIntVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nintvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      return scip->transprob->nintvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of implicit integer active problem variables
 *
 *  @return the number of implicit integer active problem variables
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNImplVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNImplVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nimplvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      return scip->transprob->nimplvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of continuous active problem variables
 *
 *  @return the number of continuous active problem variables
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNContVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNContVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->ncontvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      return scip->transprob->ncontvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}


/** gets number of active problem variables with a non-zero objective coefficient
 *
 *  @note In case of the original problem the number of variables is counted. In case of the transformed problem the
 *        number of variables is just returned since it is stored internally
 *
 *  @return the number of active problem variables with a non-zero objective coefficient
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNObjVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNObjVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobGetNObjVars(scip->origprob, scip->set);

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return SCIPprobGetNObjVars(scip->transprob, scip->set);

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}


/** gets array with fixed and aggregated problem variables; data may become invalid after
 *  calls to SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 *
 *  @return an array with fixed and aggregated problem variables; data may become invalid after
 *          calls to SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_VAR** SCIPgetFixedVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetFixedVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return NULL;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->fixedvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of fixed or aggregated problem variables
 *
 *  @return the number of fixed or aggregated problem variables
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNFixedVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNFixedVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return 0;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nfixedvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets variables of the original problem along with the numbers of different variable types; data may become invalid
 *  after a call to SCIPchgVarType()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPgetOrigVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetOrigVarsData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( vars != NULL )
      *vars = scip->origprob->vars;
   if( nvars != NULL )
      *nvars = scip->origprob->nvars;
   if( nbinvars != NULL )
      *nbinvars = scip->origprob->nbinvars;
   if( nintvars != NULL )
      *nintvars = scip->origprob->nintvars;
   if( nimplvars != NULL )
      *nimplvars = scip->origprob->nimplvars;
   if( ncontvars != NULL )
      *ncontvars = scip->origprob->ncontvars;

   return SCIP_OKAY;
}

/** gets array with original problem variables; data may become invalid after
 *  a call to SCIPchgVarType()
 *
 *  @return an array with original problem variables; data may become invalid after
 *          a call to SCIPchgVarType()
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_VAR** SCIPgetOrigVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetOrigVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->vars;
}

/** gets number of original problem variables
 *
 *  @return the number of original problem variables
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
int SCIPgetNOrigVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNOrigVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->nvars;
}

/** gets number of binary variables in the original problem
 *
 *  @return the number of binary variables in the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
int SCIPgetNOrigBinVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNOrigBinVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->nbinvars;
}

/** gets the number of integer variables in the original problem
 *
 *  @return the number of integer variables in the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
int SCIPgetNOrigIntVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNOrigIntVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->nintvars;
}

/** gets number of implicit integer variables in the original problem
 *
 *  @return the number of implicit integer variables in the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
int SCIPgetNOrigImplVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNOrigImplVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->nimplvars;
}

/** gets number of continuous variables in the original problem
 *
 *  @return the number of continuous variables in the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
int SCIPgetNOrigContVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNOrigContVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->ncontvars;
}

/** gets number of all problem variables created during creation and solving of problem;
 *  this includes also variables that were deleted in the meantime
 *
 *  @return the number of all problem variables created during creation and solving of problem;
 *          this includes also variables that were deleted in the meantime
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
int SCIPgetNTotalVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNTotalVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert(scip->stat != NULL);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      return scip->stat->nvaridx;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}


/** gets variables of the original or transformed problem along with the numbers of different variable types;
 *  the returned problem space (original or transformed) corresponds to the given solution;
 *  data may become invalid after calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and
 *  SCIPmultiaggregateVar()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetSolVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that selects the problem space, NULL for current solution */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetSolVarsData", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->set->stage == SCIP_STAGE_PROBLEM || (sol != NULL && SCIPsolIsOriginal(sol)) )
   {
      if( vars != NULL )
         *vars = scip->origprob->vars;
      if( nvars != NULL )
         *nvars = scip->origprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->origprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->origprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->origprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->origprob->ncontvars;
   }
   else
   {
      if( vars != NULL )
         *vars = scip->transprob->vars;
      if( nvars != NULL )
         *nvars = scip->transprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->transprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->transprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->transprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->transprob->ncontvars;
   }

   return SCIP_OKAY;
}

/** returns variable of given name in the problem, or NULL if not existing
 *
 *  @return variable of given name in the problem, or NULL if not existing
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_VAR* SCIPfindVar(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of variable to find */
   )
{
   SCIP_VAR* var;

   assert(name != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPfindVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobFindVar(scip->origprob, name);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      var = SCIPprobFindVar(scip->transprob, name);
      if( var == NULL )
         return SCIPprobFindVar(scip->origprob, name);
      else
         return var;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** returns TRUE iff all potential variables exist in the problem, and FALSE, if there may be additional variables,
 *  that will be added in pricing and improve the objective value
 *
 *  @return TRUE, if all potential variables exist in the problem; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Bool SCIPallVarsInProb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPallVarsInProb", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return (scip->set->nactivepricers == 0);
}

/** adds constraint to the problem; if constraint is only valid locally, it is added to the local subproblem of the
 *  current node (and all of its subnodes); otherwise it is added to the global problem;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPaddCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddCons", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   {
      SCIP_CALL( SCIPprobAddCons(scip->origprob, scip->set, scip->stat, cons) );

      if( scip->set->reopt_enable )
      {
         SCIP_CALL( SCIPreoptAddCons(scip->reopt, scip->set, scip->mem->probmem, cons) );
      }
   }
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
      SCIP_CALL( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
      assert( SCIPtreeGetCurrentDepth(scip->tree) >= 0 ||  scip->set->stage == SCIP_STAGE_PRESOLVED
         || scip->set->stage == SCIP_STAGE_INITSOLVE );
      if( SCIPtreeGetCurrentDepth(scip->tree) <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
         SCIPconsSetLocal(cons, FALSE);
      if( SCIPconsIsGlobal(cons) )
      {
         SCIP_CALL( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
      }
      else
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) > SCIPtreeGetEffectiveRootDepth(scip->tree));
         SCIP_CALL( SCIPnodeAddCons(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->tree, cons) );
      }
      return SCIP_OKAY;

   case SCIP_STAGE_EXITSOLVE:
      SCIP_CALL( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was added, or from the problem, if it was a problem constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPdelCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to delete */
   )
{
   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPdelCons", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(cons->addconssetchg == NULL);
      SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->reopt) );
      return SCIP_OKAY;

      /* only added constraints can be removed in (de-)initialization process of presolving, otherwise the reduction
       * might be wrong
       */
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
      assert(SCIPconsIsAdded(cons));
      /*lint -fallthrough*/

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_EXITSOLVE:
      SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->reopt) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** returns original constraint of given name in the problem, or NULL if not existing
 *
 *  @return original constraint of given name in the problem, or NULL if not existing
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_CONS* SCIPfindOrigCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint to find */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPfindOrigCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      return SCIPprobFindCons(scip->origprob, name);

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** returns constraint of given name in the problem, or NULL if not existing
 *
 *  @return constraint of given name in the problem, or NULL if not existing
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_CONS* SCIPfindCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint to find */
   )
{
   SCIP_CONS* cons;

   assert(name != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPfindCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobFindCons(scip->origprob, name);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      cons = SCIPprobFindCons(scip->transprob, name);
      if( cons == NULL )
         return SCIPprobFindCons(scip->origprob, name);
      else
         return cons;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of upgraded constraints
 *
 *  @return number of upgraded constraints
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNUpgrConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNUpgrConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return 0;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->stat->npresolupgdconss;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets total number of globally valid constraints currently in the problem
 *
 *  @return total number of globally valid constraints currently in the problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nconss;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nconss;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets array of globally valid constraints currently in the problem
 *
 *  @return array of globally valid constraints currently in the problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @warning If your are using the method SCIPaddCons(), it can happen that the internal constraint array (which is
 *           accessed via this method) gets resized. This can invalid the pointer which is returned by this method.
 */
SCIP_CONS** SCIPgetConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->conss;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->conss;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets total number of constraints in the original problem
 *
 *  @return total number of constraints in the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
int SCIPgetNOrigConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNOrigConss", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->nconss;
}

/** gets array of constraints in the original problem
 *
 *  @return array of constraints in the original problem
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_CONS** SCIPgetOrigConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetOrigConss", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->origprob->conss;
}

/** computes the number of check constraint in the current node (loop over all constraint handler and cumulates the
 *  number of check constraints)
 *
 *  @return returns the number of check constraints
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetNCheckConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;
   int ncheckconss;
   int c;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNCheckConss", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);
   assert(conshdlrs != NULL);

   ncheckconss = 0;

   /* loop over all constraint handler and collect the number of constraints which need to be checked */
   for( c = 0; c < nconshdlrs; ++c )
   {
      assert(conshdlrs[c] != NULL);
      ncheckconss += SCIPconshdlrGetNCheckConss(conshdlrs[c]);
   }

   return ncheckconss;
}

/*
 * local subproblem methods
 */

/** adds a conflict to a given node or globally to the problem if @p node == NULL.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note this method will release the constraint
 */
SCIP_RETCODE SCIPaddConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to add conflict (or NULL if global) */
   SCIP_CONS*            cons,               /**< constraint representing the conflict */
   SCIP_NODE*            validnode,          /**< node at whichaddConf the constraint is valid (or NULL) */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             iscutoffinvolved    /**< is a cutoff bound involved in this conflict */
   )
{
   SCIP_Real primalbound;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(scip->conflictstore != NULL);
   assert(conftype != SCIP_CONFTYPE_BNDEXCEEDING || iscutoffinvolved);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConflict", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( iscutoffinvolved )
      primalbound = SCIPgetCutoffbound(scip);
   else
      primalbound = -SCIPinfinity(scip);

   /* add a global conflict */
   if( node == NULL )
   {
      SCIP_CALL( SCIPaddCons(scip, cons) );
   }
   /* add a local conflict */
   else
   {
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, validnode) );
   }

   if( node == NULL || SCIPnodeGetType(node) != SCIP_NODETYPE_PROBINGNODE )
   {
      /* add the conflict to the conflict store */
      SCIP_CALL( SCIPconflictstoreAddConflict(scip->conflictstore, scip->mem->probmem, scip->set, scip->stat, scip->tree,
            scip->transprob, scip->reopt, cons, conftype, iscutoffinvolved, primalbound) );
   }

   /* mark constraint to be a conflict */
   SCIPconsMarkConflict(cons);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** tries to remove conflicts depending on an old cutoff bound if the improvement of the new incumbent is good enough
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPclearConflictStore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENT*           event               /**< event data */
   )
{
   assert(scip != NULL);
   assert(event != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND);
   assert(SCIPeventGetSol(event) != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPclearConflictStore", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconflictstoreCleanNewIncumbent(scip->conflictstore, scip->set, scip->stat, scip->mem->probmem,
         scip->transprob, scip->reopt, scip->primal->cutoffbound) );

   return SCIP_OKAY;
}

/** adds constraint to the given node (and all of its subnodes), even if it is a global constraint;
 *  It is sometimes desirable to add the constraint to a more local node (i.e., a node of larger depth) even if
 *  the constraint is also valid higher in the tree, for example, if one wants to produce a constraint which is
 *  only active in a small part of the tree although it is valid in a larger part.
 *  In this case, one should pass the more global node where the constraint is valid as "validnode".
 *  Note that the same constraint cannot be added twice to the branching tree with different "validnode" parameters.
 *  If the constraint is valid at the same node as it is inserted (the usual case), one should pass NULL as "validnode".
 *  If the "validnode" is the root node, it is automatically upgraded into a global constraint, but still only added to
 *  the given node. If a local constraint is added to the root node, it is added to the global problem instead.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddConsNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to add constraint to */
   SCIP_CONS*            cons,               /**< constraint to add */
   SCIP_NODE*            validnode           /**< node at which the constraint is valid, or NULL */
   )
{
   assert(cons != NULL);
   assert(node != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConsNode", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( validnode != NULL )
   {
      int validdepth;

      validdepth = SCIPnodeGetDepth(validnode);
      if( validdepth > SCIPnodeGetDepth(node) )
      {
         SCIPerrorMessage("cannot add constraint <%s> valid in depth %d to a node of depth %d\n",
            SCIPconsGetName(cons), validdepth, SCIPnodeGetDepth(node));
         return SCIP_INVALIDDATA;
      }
      if( cons->validdepth != -1 && cons->validdepth != validdepth )
      {
         SCIPerrorMessage("constraint <%s> is already marked to be valid in depth %d - cannot mark it to be valid in depth %d\n",
            SCIPconsGetName(cons), cons->validdepth, validdepth);
         return SCIP_INVALIDDATA;
      }
      if( validdepth <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
         SCIPconsSetLocal(cons, FALSE);
      else
         cons->validdepth = validdepth;
   }

   if( SCIPnodeGetDepth(node) <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
   {
      SCIPconsSetLocal(cons, FALSE);
      SCIP_CALL( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
   }
   else
   {
      SCIP_CALL( SCIPnodeAddCons(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, cons) );
   }

   return SCIP_OKAY;
}

/** adds constraint locally to the current node (and all of its subnodes), even if it is a global constraint;
 *  It is sometimes desirable to add the constraint to a more local node (i.e., a node of larger depth) even if
 *  the constraint is also valid higher in the tree, for example, if one wants to produce a constraint which is
 *  only active in a small part of the tree although it is valid in a larger part.
 *
 *  If the constraint is valid at the same node as it is inserted (the usual case), one should pass NULL as "validnode".
 *  If the "validnode" is the root node, it is automatically upgraded into a global constraint, but still only added to
 *  the given node. If a local constraint is added to the root node, it is added to the global problem instead.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note The same constraint cannot be added twice to the branching tree with different "validnode" parameters. This is
 *        the case due to internal data structures and performance issues. In such a case you should try to realize your
 *        issue using the method SCIPdisableCons() and SCIPenableCons() and control these via the event system of SCIP.
 */
SCIP_RETCODE SCIPaddConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to add */
   SCIP_NODE*            validnode           /**< node at which the constraint is valid, or NULL */
   )
{
   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConsLocal", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPaddConsNode(scip, SCIPtreeGetCurrentNode(scip->tree), cons, validnode) );

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the given node (and all subnodes);
 *  if the method is called at the root node, the constraint is globally deleted from the problem;
 *  the constraint deletion is being remembered at the given node, s.t. after leaving the node's subtree, the constraint
 *  is automatically enabled again, and after entering the node's subtree, it is automatically disabled;
 *  this may improve performance because redundant checks on this constraint are avoided, but it consumes memory;
 *  alternatively, use SCIPdisableCons()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPdelConsNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to disable constraint in */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   )
{
   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPdelConsNode", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* only added constraints can be removed in (de-)initialization process of presolving, otherwise the reduction
    * might be wrong
    */
   if( scip->set->stage == SCIP_STAGE_INITPRESOLVE || scip->set->stage == SCIP_STAGE_EXITPRESOLVE )
      assert(SCIPconsIsAdded(cons));

   if( SCIPnodeGetDepth(node) <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
   {
      SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->reopt) );
   }
   else
   {
      SCIP_CALL( SCIPnodeDelCons(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, cons) );
   }

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the current node (and all subnodes);
 *  if the method is called during problem modification or at the root node, the constraint is globally deleted from
 *  the problem;
 *  the constraint deletion is being remembered at the current node, s.t. after leaving the current subtree, the
 *  constraint is automatically enabled again, and after reentering the current node's subtree, it is automatically
 *  disabled again;
 *  this may improve performance because redundant checks on this constraint are avoided, but it consumes memory;
 *  alternatively, use SCIPdisableCons()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 *
 */
SCIP_RETCODE SCIPdelConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   )
{
   SCIP_NODE* node;

   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPdelConsLocal", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(cons->addconssetchg == NULL);
      SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->reopt) );
      return SCIP_OKAY;

      /* only added constraints can be removed in (de-)initialization process of presolving, otherwise the reduction
       * might be wrong
       */
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
      assert(SCIPconsIsAdded(cons));
      /*lint -fallthrough*/

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      node = SCIPtreeGetCurrentNode(scip->tree);

      if( SCIPnodeGetDepth(node) <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
      {
         SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->reopt) );
      }
      else
      {
         SCIP_CALL( SCIPnodeDelCons(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, cons) );
      }
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** gets estimate of best primal solution w.r.t. original problem contained in current subtree
 *
 *  @return estimate of best primal solution w.r.t. original problem contained in current subtree
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetLocalOrigEstimate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODE* node;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLocalOrigEstimate", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);
   return node != NULL ? SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPnodeGetEstimate(node)) : SCIP_INVALID;
}

/** gets estimate of best primal solution w.r.t. transformed problem contained in current subtree
 *
 *  @return estimate of best primal solution w.r.t. transformed problem contained in current subtree
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetLocalTransEstimate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODE* node;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLocalTransEstimate", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);

   return node != NULL ? SCIPnodeGetEstimate(node) : SCIP_INVALID;
}

/** gets dual bound of current node
 *
 *  @return dual bound of current node
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetLocalDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODE* node;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLocalDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);
   return node != NULL ? SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPnodeGetLowerbound(node)) : SCIP_INVALID;
}

/** gets lower bound of current node in transformed problem
 *
 *  @return lower bound  of current node in transformed problem
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetLocalLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODE* node;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLocalLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);

   return node != NULL ? SCIPnodeGetLowerbound(node) : SCIP_INVALID;
}

/** gets dual bound of given node
 *
 *  @return dual bound of a given node
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetNodeDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node to get dual bound for */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNodeDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPnodeGetLowerbound(node));
}

/** gets lower bound of given node in transformed problem
 *
 *  @return lower bound  of given node in transformed problem
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetNodeLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node to get dual bound for */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNodeLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPnodeGetLowerbound(node);
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the current node's dual bound (in
 *  original problem space), sets the current node's dual bound to the new value
 *
 *  @note the given new bound has to be a dual bound, i.e., it has to be valid for the original problem.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPupdateLocalDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdateLocalDualbound", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* since no root node, for which we could update the dual bound, has been create yet, update the dual bound stored in
       * the problem data
       */
      SCIPprobUpdateDualbound(scip->origprob, newbound);
      break;

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
      /* since no root node, for which we could update the dual bound, has been create yet, update the dual bound stored in
       * the problem data
       */
      SCIPprobUpdateDualbound(scip->transprob, SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, newbound));
      break;

   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, SCIPtreeGetCurrentNode(scip->tree), SCIPprobInternObjval(scip->transprob, scip->origprob, scip->set, newbound)) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return SCIP_INVALIDCALL; /*lint !e527*/
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** if given value is larger than the current node's lower bound (in transformed problem), sets the current node's
 *  lower bound to the new value
 *
 *  @note the given new bound has to be a lower bound, i.e., it has to be valid for the transformed problem.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPupdateLocalLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdateLocalLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
      /* since no root node, for which we could update the lower bound, has been created yet, update the dual bound stored
       * in the problem data
       */
      SCIPprobUpdateDualbound(scip->transprob, SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, newbound));
      break;

   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, SCIPtreeGetCurrentNode(scip->tree), newbound) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return SCIP_INVALIDCALL; /*lint !e527*/
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the node's dual bound,
 *  sets the node's dual bound to the new value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPupdateNodeDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to update dual bound for */
   SCIP_Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdateNodeDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPupdateNodeLowerbound(scip, node, SCIPprobInternObjval(scip->transprob, scip->origprob, scip->set, newbound)) );

   return SCIP_OKAY;
}

/** if given value is larger than the node's lower bound (in transformed problem), sets the node's lower bound
 *  to the new value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPupdateNodeLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to update lower bound for */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdateNodeLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(node, scip->stat, scip->set, scip->tree, scip->transprob, scip->origprob, newbound);

   /* if lowerbound exceeds the cutoffbound the node will be marked to be cutoff
    *
    * If the node is an inner node (,not a child node,) we need to cutoff the node manually if we exceed the
    * cutoffbound. This is only relevant if a user updates the lower bound; in the main solving process of SCIP the
    * lowerbound is only changed before branching and the given node is always a child node. Therefore, we only check
    * for a cutoff here in the user function instead of in SCIPnodeUpdateLowerbound().
    */
   if( SCIPisGE(scip, newbound, scip->primal->cutoffbound) )
   {
      SCIP_CALL( SCIPnodeCutoff(node, scip->set, scip->stat, scip->tree, scip->transprob, scip->origprob, scip->reopt,
            scip->lp, scip->mem->probmem) );
   }

   return SCIP_OKAY;
}

/** change the node selection priority of the given child
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgChildPrio(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            child,              /**< child to update the node selection priority */
   SCIP_Real             priority            /**< node selection priority value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgChildPrio", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPnodeGetType(child) != SCIP_NODETYPE_CHILD )
      return SCIP_INVALIDDATA;

   SCIPchildChgNodeselPrio(scip->tree, child, priority);

   return SCIP_OKAY;
}
