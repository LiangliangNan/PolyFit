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

/**@file   benders_default.c
 * @ingroup OTHER_CFILES
 * @brief  default Benders' decomposition plugin
 * @author Stephen J. Maher
 *
 * The default Benders' decomposition plugin is provided to simplify the interaction the Benders' decomposition
 * framework within SCIP. This plugin is included in the SCIP package by default. Using the default Benders'
 * decomposition plugin, a problem can be solved by Benders' decomposition by calling
 *
 * SCIPcreateBendersDefault(master problem, array of subproblems, number of subproblems)
 *
 * where "master problem" is a SCIP instance of the master problem, "array of subproblems" is an array of SCIP instances
 * that are the Benders' decomposition subproblems and "number of subproblems" is an integer indicating the number of
 * subproblems for this decomposition.
 *
 * A key feature of the default Benders' decomposition plugin is the automatic generation of the variable mapping
 * between the variables of the master problem and the subproblems.
 *
 * In the traditional application of Benders' decomposition, master problem variables are fixed to a solution value and
 * modify the RHS of the second stage constraints. The implementation within SCIP requires that a variable is created
 * in the subproblem for every master problem variable that appears in the subproblem constraints. This variable MUST
 * have the same name as the corresponding variable in the master problem. This name is used to automatically generate
 * the mapping between the master problem and the corresponding subproblem variables.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/benders_default.h"
#include "scip/bendersdefcuts.h"
#include "scip/pub_benders.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip.h"

#define BENDERS_NAME                "default"
#define BENDERS_DESC                "default implementation of Benders' decomposition"
#define BENDERS_PRIORITY            0
#define BENDERS_CUTLP            TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO        TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX         TRUE   /**< should Benders' cut be generated for relaxation solutions */
#define BENDERS_SHAREAUXVARS    FALSE   /**< should this Benders' share the highest priority Benders' aux vars */


/*
 * Data structures
 */

/** Benders' decomposition data */
struct SCIP_BendersData
{
   SCIP**                subproblems;        /**< the Benders' decomposition subproblems */
   SCIP_HASHMAP*         mastervartosubindex;/**< hash map from the master variable to an index for the subproblemn variables */
   SCIP_HASHMAP*         subvartomastervar;  /**< hashmap from the subproblem variable to the master variable */
   SCIP_VAR***           subproblemvars;     /**< the subproblem variables corresponding to master problem variables */
   int                   nmastervars;        /**< the number of variables in the master problem */
   int                   nsubproblems;       /**< the number of subproblems */
   SCIP_Bool             created;            /**< flag to indicate that the Benders' decomposition Data was created */
   SCIP_Bool             subprobscopied;     /**< were the subproblems copied during the SCIP copy */
   SCIP_Bool             mappingcreated;     /**< flag to indicate whether the variable mapping has been created */
};




/*
 * Local methods
 */

/** creates the Benders' decomposition data */
static
SCIP_RETCODE createBendersData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   SCIP_BENDERSDATA**    bendersdata,        /**< the Benders' decomposition data */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   int i;

   assert(scip != NULL);
   assert(subproblems != NULL);

   (*bendersdata)->nsubproblems = nsubproblems;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*bendersdata)->subproblems, nsubproblems) );

   /* Copying the subproblem to the Benders' decomposition data. */
   for( i = 0; i < nsubproblems; i++ )
      (*bendersdata)->subproblems[i] = subproblems[i];

   (*bendersdata)->created = TRUE;

   return SCIP_OKAY;
}


/** Creates the variable mappings between the master problem and the corresponding variable in the subproblem.
 *
 *  TODO: add the functionality to allow the user to provide an array of hashmaps for mapping between the master problem
 *  variables and the corresponding subproblem variables.
 *  TODO: check for uniqueness of names in this function.
 */
static
SCIP_RETCODE createVariableMappings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_VAR** vars;
   int nsubproblems;
   int nvars;
   char varname[SCIP_MAXSTRLEN];
   int i;
   int j;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   nsubproblems = bendersdata->nsubproblems;

   /* getting the master problem variable data */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   bendersdata->nmastervars = nvars;

   /* creating the hashmaps for the mapping between the master variables and the sub variables */
   SCIP_CALL( SCIPhashmapCreate(&bendersdata->mastervartosubindex, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&bendersdata->subvartomastervar, SCIPblkmem(scip), nvars*nsubproblems) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &bendersdata->subproblemvars, nsubproblems) );
   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &bendersdata->subproblemvars[i], nvars) );
   }

   /* this loop stores a mapping between the master problem variables and their counterpart in the subproblems. For each
    * master problem variable, the variable name is used to search for any corresponding variables in each of the
    * subproblems. If a corresponding variable exists, then a mapping is inserted into subvartomastervar and
    * mastervartosubvar hashmaps
    */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_VAR* origvar;
      SCIP_VAR* subvar;
      SCIP_Real scalar;
      SCIP_Real constant;
      const char* origvarname;
      int charcount = SCIPgetSubscipDepth(scip)*2;

      /* getting the original variable for the master variable
       * NOTE: This retrieved variable is the original variable. It may be a bug in regards to other parts of the code.
       * The process maps the subproblem variable to the original master variable. It was original supposed to be a
       * mapping between the subproblem variables and the transformed master variable.
       */
      origvar = vars[i];
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      /* retrieving the var name */
      origvarname = SCIPvarGetName(origvar);
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s", &origvarname[charcount]);

      /* retrieving the subproblem variable for the given master variable */
      for( j = 0; j < nsubproblems; j++ )
      {
         /* find the corresponding subproblem variable for a given master problem variable using the variable name. */
         subvar = SCIPfindVar(bendersdata->subproblems[j], varname);

         /* adding the subvariable to master variable mapping into the hash map */
         if( subvar != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(bendersdata->subvartomastervar, subvar, origvar) );
         }

         /* storing the subproblem variable */
         bendersdata->subproblemvars[j][i] = subvar;

         if( subvar != NULL )
         {
            SCIP_CALL( SCIPcaptureVar(bendersdata->subproblems[j], bendersdata->subproblemvars[j][i]) );
         }
      }

      /* storing the mapping of the master variable to the variable index */
      SCIP_CALL( SCIPhashmapInsertInt(bendersdata->mastervartosubindex, vars[i], i) );
   }

   bendersdata->mappingcreated = TRUE;

   return SCIP_OKAY;
}



/*
 * Callback methods for Benders' decomposition
 */

/** copy method for Benders' decomposition plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BENDERSCOPY(bendersCopyDefault)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;   /* the source Benders' decomposition data */

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   /* including the Benders' decomposition in the target SCIP.
    * NOTE: this method uses the same subproblems as the main SCIP. In a parallel setting, this will not be thread safe.
    * It would be cleaner to copy the subproblems also.
    */
   SCIP_CALL( SCIPincludeBendersDefault(scip) );

   /* if the Benders' decomposition is active, then it must be created in the copy */
   if( SCIPbendersIsActive(benders) )
   {
      SCIP** subproblems;
      int i;

      /* copying the subproblems if the threadsafe flag was set to TRUE */
      if( threadsafe )
      {
         /* allocating memory for the subproblems array */
         SCIP_CALL( SCIPallocBufferArray(scip, &subproblems, bendersdata->nsubproblems) );

         for( i = 0; i < bendersdata->nsubproblems; i++ )
         {
            SCIP_Bool valid;

            /* creating the SCIP instance for the subproblem */
            SCIP_CALL( SCIPcreate(&subproblems[i]) );

            /* the original problem is copied so that the variable mappings are created correctly.
             * TODO: use a varmap to create the mappings for the copies
             */
            SCIP_CALL( SCIPcopyOrig(bendersdata->subproblems[i], subproblems[i], NULL, NULL, "", TRUE, FALSE, FALSE,
                  &valid) );
            assert(valid);
         }
      }
      else
         subproblems = bendersdata->subproblems;

      SCIP_CALL( SCIPcreateBendersDefault(scip, subproblems, bendersdata->nsubproblems) );

      /* freeing the buffer memory for the subproblems */
      if( threadsafe )
      {
         SCIP_BENDERS* targetbenders;
         SCIP_BENDERSDATA* targetbendersdata;

         targetbenders = SCIPfindBenders(scip, BENDERS_NAME);
         assert(targetbenders != NULL);

         targetbendersdata = SCIPbendersGetData(targetbenders);

         /* indicating that the subproblems have been copied */
         targetbendersdata->subprobscopied = TRUE;

         SCIPfreeBufferArray(scip, &subproblems);
      }
   }

   return SCIP_OKAY;
}

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
/**! [SnippetBendersFreeDefault] */
static
SCIP_DECL_BENDERSFREE(bendersFreeDefault)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   assert(bendersdata != NULL);

   /* should have been freed in bendersExitDefault (if mappingcreated), or not been created at the first place */
   assert(bendersdata->subproblemvars == NULL);
   assert(bendersdata->subvartomastervar == NULL);
   assert(bendersdata->mastervartosubindex == NULL);
   if( bendersdata->created )
   {
      /* if the subproblems were copied, then the copy needs to be freed */
      if( bendersdata->subprobscopied )
      {
         for( i = bendersdata->nsubproblems - 1; i >= 0; i-- )
         {
            SCIP_CALL( SCIPfree(&bendersdata->subproblems[i]) );
         }
      }

      SCIPfreeBlockMemoryArray(scip, &bendersdata->subproblems, bendersdata->nsubproblems);
   }

   SCIPfreeBlockMemory(scip, &bendersdata);

   return SCIP_OKAY;
}
/**! [SnippetBendersFreeDefault] */


/** initialization method of Benders' decomposition (called after problem was transformed) */
static
SCIP_DECL_BENDERSINIT(bendersInitDefault)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(benders != NULL);

   /* creating the variable mappings */
   SCIP_CALL( createVariableMappings(scip, benders) );

   return SCIP_OKAY;
}

/** deinitialization method of Benders' decomposition (called before transformed problem is freed and the Benders'
 * decomposition is active)
 */
static
SCIP_DECL_BENDERSEXIT(bendersExitDefault)
{
   SCIP_BENDERSDATA* bendersdata;
   int i;
   int j;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   assert(bendersdata != NULL);

   if( bendersdata->mappingcreated )
   {
      for( i = bendersdata->nsubproblems - 1; i >= 0; i-- )
      {
         for( j = 0; j < bendersdata->nmastervars; j++ )
         {
            if( bendersdata->subproblemvars[i][j] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(bendersdata->subproblems[i], &bendersdata->subproblemvars[i][j]) );
            }
         }
         SCIPfreeBlockMemoryArray(scip, &bendersdata->subproblemvars[i], bendersdata->nmastervars);
      }
      SCIPfreeBlockMemoryArray(scip, &bendersdata->subproblemvars, bendersdata->nsubproblems);

      /* free hash map */
      SCIPhashmapFree(&bendersdata->subvartomastervar);
      SCIPhashmapFree(&bendersdata->mastervartosubindex);
   }

   return SCIP_OKAY;
}

/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
/**! [SnippetBendersGetvarDefault] */
static
SCIP_DECL_BENDERSGETVAR(bendersGetvarDefault)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   SCIP_VAR* origvar;
   SCIP_Real scalar;
   SCIP_Real constant;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);

   bendersdata = SCIPbendersGetData(benders);

   if( probnumber == -1 )
   {
      origvar = var;
      /* The variable needs to be transformed back into an original variable. If the variable is already original, then
       * this function just returns the same variable
       */
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      /* using the original variable, the master variable can be retrieved from the hash map */
      (*mappedvar) = (SCIP_VAR*) SCIPhashmapGetImage(bendersdata->subvartomastervar, origvar);

      if( (*mappedvar) == NULL )
         (*mappedvar) = (SCIP_VAR*) SCIPhashmapGetImage(bendersdata->subvartomastervar, var);
   }
   else
   {
      int masterindex;
      /* The variable needs to be transformed back into an original variable. If the variable is already original, then
       * this function just returns the same variable
       */

      /* we are requesting the subproblem variable for a master problem variable
       * The master problem variable is a transformed variable. The original variable is not required.
       * NOTE: Currently the original variable is being used. This may not be correct and should be the transformed
       * variable.
       */
      masterindex = SCIPhashmapGetImageInt(bendersdata->mastervartosubindex, var);
      (*mappedvar) = bendersdata->subproblemvars[probnumber][masterindex];
   }

   return SCIP_OKAY;
}
/**! [SnippetBendersGetvarDefault] */

/** the method for creating the Benders' decomposition subproblem. This method is called during the initialisation stage
 *  (after the master problem was transformed)
 *
 *  This method must create the SCIP instance for the subproblem and add the required variables and constraints. In
 *  addition, the settings required for the solving the problem must be set here. However, some settings will be
 *  overridden by the standard solving method included in the Benders' decomposition framework. If a special solving
 *  method is desired, the user can implement the bendersSolvesubDefault callback.
 */
static
SCIP_DECL_BENDERSCREATESUB(bendersCreatesubDefault)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   /* adding the subproblem to the Benders' decomposition structure */
   SCIP_CALL( SCIPaddBendersSubproblem(scip, benders, bendersdata->subproblems[probnumber]) );

   return SCIP_OKAY;
}


/*
 * Benders' decomposition specific interface methods
 */

/** Creates a default Benders' decomposition algorithm and activates it in SCIP
 *
 * @note Every variable that appears in the subproblem constraints must be created in the corresponding subproblem with
 * the same name as in the master problem.
 *
 * @note The default Benders' decomposition implementation relies on unique variable names in the master problem and in
 * each of the subproblems. This is required because a variable mapping is made between the master problem variables and
 * the counterparts in the subproblems. This mapping is created using the variable names.
 */
SCIP_RETCODE SCIPcreateBendersDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_BENDERS* benders;
   SCIP_BENDERSDATA* bendersdata;
   int maxrestarts;

   assert(scip != NULL);
   assert(subproblems != NULL);
   assert(nsubproblems > 0);

   benders = SCIPfindBenders(scip, BENDERS_NAME);
   bendersdata = SCIPbendersGetData(benders);

   /* turning restarts off */
   SCIP_CALL( SCIPgetIntParam(scip, "presolving/maxrestarts", &maxrestarts) );
   if( SCIPisParamFixed(scip, "presolving/maxrestarts") && maxrestarts != 0)
   {
      SCIPerrorMessage("The number of restarts is fixed to %d. The default Benders' decomposition requires the number"
         " of restarts to be 0.", maxrestarts);
      return SCIP_ERROR;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
      SCIP_CALL( SCIPfixParam(scip, "presolving/maxrestarts") );
   }

   SCIP_CALL( createBendersData(scip, subproblems, &bendersdata, nsubproblems) );

   SCIP_CALL( SCIPactivateBenders(scip, benders, nsubproblems) );

   return SCIP_OKAY;
}

/** creates the default Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeBendersDefault(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_BENDERS* benders;

   /* create default Benders' decomposition data */
   bendersdata = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &bendersdata) );
   BMSclearMemory(bendersdata);

   benders = NULL;

   /* include Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, BENDERS_CUTLP,
         BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, BENDERS_SHAREAUXVARS, bendersGetvarDefault, bendersCreatesubDefault,
         bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBendersCopy(scip, benders, bendersCopyDefault) );
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeDefault) );
   SCIP_CALL( SCIPsetBendersInit(scip, benders, bendersInitDefault) );
   SCIP_CALL( SCIPsetBendersExit(scip, benders, bendersExitDefault) );

   /* OPTIONAL: including the default cuts for Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );

   return SCIP_OKAY;
}
