/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   debug.c
 * @brief  methods for debugging
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/scip.h"
#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/struct_scip.h"

#ifdef WITH_DEBUG_SOLUTION

#define SCIP_HASHSIZE_DEBUG        500    /**< minimum size of hash map for storing whether a solution is valid for the node */

struct SCIP_DebugSolData
{
   char**                solnames;           /**< variable names in the solution */
   SCIP_Real*            solvals;            /**< solution value array (only nonzero entries) */
   int                   nsolvals;           /**< number of entries in the debug solution */
   int                   solsize;            /**< size of the array entries */
   SCIP_SOL*             debugsol;           /**< a debug solution */
   SCIP_STAGE            debugsolstage;      /**< solving stage of debug solution */
   SCIP_HASHMAP*         solinnode;          /**< maps nodes to bools, storing whether the solution is valid for the node */
   SCIP_Bool             falseptr;           /**< pointer to value FALSE used for hashmap */
   SCIP_Bool             trueptr;            /**< pointer to value TRUE used for hashmap */
   SCIP_Bool             solisachieved;      /**< means if current best solution is better than the given debug solution */
   SCIP_Real             debugsolval;        /**< objective value for debug solution */
   SCIP_Bool             debugsoldisabled;   /**< flag indicating if debugging of solution was disabled or not */
   SCIP_Bool             warningprinted;     /**< flag indicating if a warning was already printed */
};


/** creates debug solution data */
SCIP_RETCODE SCIPdebugSolDataCreate(
   SCIP_DEBUGSOLDATA**   debugsoldata        /**< pointer to debug solution data */
   )
{
   assert(debugsoldata != NULL);

   SCIP_ALLOC( BMSallocMemory(debugsoldata) );

   (*debugsoldata)->solnames = NULL;
   (*debugsoldata)->solvals = NULL;
   (*debugsoldata)->nsolvals = 0;
   (*debugsoldata)->solsize = 0;
   (*debugsoldata)->debugsol = NULL;
   (*debugsoldata)->debugsolstage = SCIP_STAGE_INIT;
   (*debugsoldata)->solinnode = NULL;
   (*debugsoldata)->falseptr = FALSE;
   (*debugsoldata)->trueptr = TRUE;
   (*debugsoldata)->solisachieved = FALSE;
   (*debugsoldata)->debugsolval = 0.0;
   (*debugsoldata)->debugsoldisabled = TRUE;
   (*debugsoldata)->warningprinted = FALSE;

   return SCIP_OKAY;
}

#ifdef SCIP_MORE_DEBUG
/** comparison method for sorting variables w.r.t. to their name */
static
SCIP_DECL_SORTPTRCOMP(sortVarsAfterNames)
{
   return strcmp(SCIPvarGetName((SCIP_VAR*)elem1), SCIPvarGetName((SCIP_VAR*)elem2));
}
#endif

/* checks whether the parameter is specified */
static
SCIP_Bool debugSolutionAvailable(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;

   assert(set != NULL);

   debugsoldata = SCIPsetGetDebugSolData(set);

   /* check whether a debug solution is specified */
    if( strcmp(set->misc_debugsol, "-") == 0 )
    {
       if( !debugsoldata->warningprinted )
       {
          SCIPmessagePrintWarning(SCIPgetMessagehdlr(set->scip), "SCIP is compiled with 'DEBUGSOL=true' but no debug solution is given:\n ");
          SCIPmessagePrintWarning(SCIPgetMessagehdlr(set->scip), "*** Please set the parameter 'misc/debugsol' and reload the problem again to use the debugging-mechanism ***\n\n");
          debugsoldata->warningprinted = TRUE;
       }
       return FALSE;
    }
    else
    {
       debugsoldata->warningprinted = FALSE;
       return TRUE;
    }
}

/** reads solution from given file into given arrays */
static
SCIP_RETCODE readSolfile(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           solfilename,        /**< solution filename to read */
   SCIP_SOL**            debugsolptr,
   SCIP_Real*            debugsolvalptr,
   SCIP_STAGE*           debugsolstageptr,
   char***               names,              /**< pointer to store the array of variable names */
   SCIP_Real**           vals,               /**< pointer to store the array of solution values */
   int*                  nvals,              /**< pointer to store the number of non-zero elements */
   int*                  valssize            /**< pointer to store the length of the variable names and solution values arrays */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* solvalues;
   SCIP_FILE* file;
   SCIP_SOL* debugsol;
   SCIP_Real debugsolval;
   int nonvalues;
   int nfound;
   int i;
   SCIP_Bool unknownvariablemessage;

   assert(set != NULL);
   assert(solfilename != NULL);
   assert(names != NULL);
   assert(*names == NULL);
   assert(vals != NULL);
   assert(*vals == NULL);
   assert(nvals != NULL);
   assert(valssize != NULL);

   printf("***** debug: reading solution file <%s>\n", solfilename);

   /* open solution file */
   file = SCIPfopen(solfilename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open solution file <%s> specified in scip/debug.h\n", solfilename);
      SCIPprintSysError(solfilename);
      return SCIP_NOFILE;
   }

   /* read data */
   nonvalues = 0;
   *valssize = 0;
   unknownvariablemessage = FALSE;

   while( !SCIPfeof(file) )
   {
      char buf[SCIP_MAXSTRLEN];
      char name[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      char valuestring[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real val;
      int nread;

      if( SCIPfgets(buf, SCIP_MAXSTRLEN, file) == NULL )
      {
         if( SCIPfeof(file) )
            break;
         else
            return SCIP_READERROR;
      }

      /* there are some lines which may preceed the solution information */
      if( strncasecmp(buf, "solution status:", 16) == 0 || strncasecmp(buf, "objective value:", 16) == 0 ||
         strncasecmp(buf, "Log started", 11) == 0 || strncasecmp(buf, "Variable Name", 13) == 0 ||
         strncasecmp(buf, "All other variables", 19) == 0 || strncasecmp(buf, "\n", 1) == 0 ||
         strncasecmp(buf, "NAME", 4) == 0 || strncasecmp(buf, "ENDATA", 6) == 0 )    /* allow parsing of SOL-format on the MIPLIB 2003 pages */
      {
         ++nonvalues;
         continue;
      }

      /* cppcheck-suppress invalidscanf */
      nread = sscanf(buf, "%s %s %s\n", name, valuestring, objstring);
      if( nread < 2 )
      {
         printf("invalid input line %d in solution file <%s>: <%s>\n", *nvals + nonvalues, solfilename, name);
         SCIPfclose(file);
         return SCIP_READERROR;
      }

      /* find the variable */
      var = SCIPfindVar(set->scip, name);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPverbMessage(set->scip, SCIP_VERBLEVEL_NORMAL, NULL, "unknown variable <%s> in line %d of solution file <%s>\n",
               name, *nvals + nonvalues, solfilename);
            SCIPverbMessage(set->scip, SCIP_VERBLEVEL_NORMAL, NULL, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* cast the value, check first for inv(alid) or inf(inite) ones that need special treatment */
      if( strncasecmp(valuestring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(valuestring, "+inf", 4) == 0 || strncasecmp(valuestring, "inf", 3) == 0 )
         val = SCIPsetInfinity(set);
      else if( strncasecmp(valuestring, "-inf", 4) == 0 )
         val = -SCIPsetInfinity(set);
      else
      {
         /* cppcheck-suppress invalidscanf */
         nread = sscanf(valuestring, "%lf", &val);
         if( nread != 1 )
         {
            SCIPerrorMessage("Invalid solution value <%s> for variable <%s> in line %d of solution file <%s>.\n",
                             valuestring, name, *nvals + nonvalues, solfilename);
            SCIPfclose(file);
            return SCIP_READERROR;
         }
      }

      /* allocate memory */
      if( *nvals >= *valssize )
      {
         *valssize = MAX(2 * *valssize, (*nvals)+1);
         SCIP_ALLOC( BMSreallocMemoryArray(names, *valssize) );
         SCIP_ALLOC( BMSreallocMemoryArray(vals, *valssize) );
      }
      assert(*nvals < *valssize);

      /* store solution value in sorted list */
      for( i = *nvals; i > 0 && strcmp(name, (*names)[i-1]) < 0; --i )
      {
         (*names)[i] = (*names)[i-1];
         (*vals)[i] = (*vals)[i-1];
      }
      SCIP_ALLOC( BMSduplicateMemoryArray(&(*names)[i], name, strlen(name)+1) );
      SCIPdebugMsg(set->scip, "found variable <%s>: value <%g>\n", (*names)[i], val);
      (*vals)[i] = val;
      (*nvals)++;
   }

   /* get memory for SCIP solution */
   SCIP_ALLOC( BMSallocMemoryArray(&vars, *valssize) );
   SCIP_ALLOC( BMSallocMemoryArray(&solvalues, *valssize) );

   debugsolval = 0.0;
   nfound = 0;

   /* get solution value */
   for( i = 0; i < *nvals; ++i)
   {
      SCIP_VAR* var;
      var = SCIPfindVar(set->scip, (*names)[i]);
      if( var != NULL )
      {
         vars[nfound] = var;
         solvalues[nfound] = (*vals)[i];
         ++nfound;
         debugsolval += (*vals)[i] * SCIPvarGetObj(var);
      }
   }
   SCIPdebugMsg(set->scip, "Debug Solution value is %g.\n", debugsolval);

#ifdef SCIP_MORE_DEBUG
   SCIPsortPtrReal((void**)vars, solvalues, sortVarsAfterNames, nfound);

   for( i = 0; i < nfound - 1; ++i)
   {
      assert(strcmp(SCIPvarGetName(vars[i]), SCIPvarGetName(vars[i + 1])) != 0);
   }
#endif

   if( debugsolptr != NULL )
   {
      /* create SCIP solution */
      SCIP_CALL( SCIPcreateOrigSol(set->scip, &debugsol, NULL) );
      *debugsolstageptr = SCIPgetStage(set->scip);

      /* set SCIP solution values */
      SCIP_CALL( SCIPsetSolVals(set->scip, debugsol, nfound, vars, solvalues ) );
   }

   BMSfreeMemoryArray(&vars);
   BMSfreeMemoryArray(&solvalues);

   if( debugsolptr != NULL )
      *debugsolptr = debugsol;

   if( debugsolvalptr != NULL )
      *debugsolvalptr = debugsolval;

   /* close file */
   SCIPfclose(file);

   printf("***** debug: read %d non-zero entries (%d variables found)\n", *nvals, nfound);

   return SCIP_OKAY;
}

/** reads feasible solution to check from file */
static
SCIP_RETCODE readSolution(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;

   assert(set != NULL);

   debugsoldata = SCIPsetGetDebugSolData(set);

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   if( debugsoldata == NULL || debugsoldata->nsolvals > 0 )
      return SCIP_OKAY;

   SCIP_CALL( readSolfile(set, set->misc_debugsol, &debugsoldata->debugsol, &debugsoldata->debugsolval,
         &debugsoldata->debugsolstage, &(debugsoldata->solnames), &(debugsoldata->solvals), &(debugsoldata->nsolvals),
         &(debugsoldata->solsize)) );

   return SCIP_OKAY;
}

/** gets value of given variable in debugging solution */
static
SCIP_RETCODE getSolutionValue(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to get solution value for */
   SCIP_Real*            val                 /**< pointer to store solution value */
   )
{
   SCIP_VAR* solvar;
   SCIP_DEBUGSOLDATA* debugsoldata;
   SCIP_Real scalar;
   SCIP_Real constant;
   const char* name;
   int left;
   int right;
   int middle;
   int cmp;

   assert(set != NULL);
   assert(var != NULL);
   assert(val != NULL);

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   debugsoldata = SCIPsetGetDebugSolData(set);
   assert(debugsoldata != NULL);

   /* allow retrieving solution values only if referring to the SCIP instance that is debugged */
   if( !SCIPdebugSolIsEnabled(set->scip) )
   {
      *val = SCIP_UNKNOWN;
      return SCIP_OKAY;
   }

   SCIP_CALL( readSolution(set) );
   SCIPsetDebugMsg(set, "Now handling variable <%s>, which has status %d, is of type %d, and was deleted: %d, negated: %d, transformed: %d\n",
      SCIPvarGetName(var), SCIPvarGetStatus(var), SCIPvarGetType(var), SCIPvarIsDeleted(var), SCIPvarIsNegated(var),SCIPvarIsTransformedOrigvar(var));

   /* ignore deleted variables */
   if( SCIPvarIsDeleted(var) )
   {
      SCIPsetDebugMsg(set, "**** unknown solution value for deleted variable <%s>\n", SCIPvarGetName(var));
      *val = SCIP_UNKNOWN;
      return SCIP_OKAY;
   }

   /* retransform variable onto original variable space */
   solvar = var;
   scalar = 1.0;
   constant = 0.0;
   if( SCIPvarIsNegated(solvar) )
   {
      scalar = -1.0;
      constant = SCIPvarGetNegationConstant(solvar);
      solvar = SCIPvarGetNegationVar(solvar);
   }

   if( SCIPvarIsTransformed(solvar) )
   {
      SCIP_CALL( SCIPvarGetOrigvarSum(&solvar, &scalar, &constant) );
      if( solvar == NULL )
      {
         /* if no original counterpart, then maybe someone added a value for the transformed variable, so search for var (or its negation) */
         SCIPsetDebugMsg(set, "variable <%s> has no original counterpart\n", SCIPvarGetName(var));
         solvar = var;
         scalar = 1.0;
         constant = 0.0;
         if( SCIPvarIsNegated(solvar) )
         {
            scalar = -1.0;
            constant = SCIPvarGetNegationConstant(solvar);
            solvar = SCIPvarGetNegationVar(solvar);
         }
      }
   }

   /* perform a binary search for the variable */
   name = SCIPvarGetName(solvar);
   left = 0;
   right = debugsoldata->nsolvals-1;
   while( left <= right )
   {
      middle = (left+right)/2;
      cmp = strcmp(name, debugsoldata->solnames[middle]);
      if( cmp < 0 )
         right = middle-1;
      else if( cmp > 0 )
         left = middle+1;
      else
      {
         *val = scalar * debugsoldata->solvals[middle] + constant;

         if( SCIPsetIsFeasLT(set, *val, SCIPvarGetLbGlobal(var)) || SCIPsetIsFeasGT(set, *val, SCIPvarGetUbGlobal(var)) )
         {
            SCIPmessagePrintWarning(SCIPgetMessagehdlr(set->scip), "invalid solution value %.15g for variable <%s>[%.15g,%.15g]\n",
               *val, SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
         }

         return SCIP_OKAY;
      }
   }
   *val = constant;

   if( SCIPsetIsFeasLT(set, *val, SCIPvarGetLbGlobal(var)) || SCIPsetIsFeasGT(set, *val, SCIPvarGetUbGlobal(var)) )
   {
      SCIPmessagePrintWarning(SCIPgetMessagehdlr(set->scip), "invalid solution value %.15g for variable <%s>[%.15g,%.15g]\n",
         *val, SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   return SCIP_OKAY;
}

/** gets pointer to the debug solution */
SCIP_RETCODE SCIPdebugGetSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol                 /**< buffer to store pointer to the debug solution */
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;

   debugsoldata = SCIPsetGetDebugSolData(scip->set);
   assert(scip != NULL);
   assert(sol != NULL);

   *sol = NULL;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(scip->set) )
      return SCIP_OKAY;

   SCIP_CALL( readSolution(scip->set) );

   if( debugsoldata->debugsol == NULL )
      return SCIP_ERROR;

   *sol = debugsoldata->debugsol;

   return SCIP_OKAY;
}

/** gets value for a variable in the debug solution
 *
 * if no value is stored for the variable, gives 0.0
 */
SCIP_RETCODE SCIPdebugGetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which to get the value */
   SCIP_Real*            val                 /**< buffer to store solution value */
   )
{
   SCIP_CALL( getSolutionValue(scip->set, var, val) );

   return SCIP_OKAY;
}

/** returns whether the debug solution is worse than the best known solution or if the debug solution was found */
static
SCIP_Bool debugSolIsAchieved(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_SOL* bestsol;
   SCIP* scip;
   SCIP_DEBUGSOLDATA* debugsoldata;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   assert(set != NULL);
   debugsoldata = SCIPsetGetDebugSolData(set);

   assert(debugsoldata != NULL);

   if( debugsoldata->solisachieved )
      return TRUE;

   assert(set != NULL);

   scip = set->scip;
   assert(scip != NULL);

   bestsol = SCIPgetBestSol(scip);

   if( bestsol != NULL )
   {
      SCIP_Real solvalue;

      /* don't check solution while in problem creation stage */
      if( SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM )
         return TRUE;

      solvalue = SCIPgetSolOrigObj(scip, bestsol);

      /* make sure a debug solution has been read, so we do not compare against the initial debugsolval == 0 */
      SCIP_CALL( readSolution(set) );

      if( (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE && SCIPsetIsLE(set, solvalue, debugsoldata->debugsolval))
            || (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE && SCIPsetIsGE(set, solvalue, debugsoldata->debugsolval)) )
         debugsoldata->solisachieved = TRUE;
   }

   return debugsoldata->solisachieved;
}

/** returns whether the solution is contained in node's subproblem */
static
SCIP_RETCODE isSolutionInNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< local node where this bound change was applied */
   SCIP_Bool*            solcontained        /**< pointer to store whether the solution is contained in node's subproblem */
   )
{
   SCIP_Bool* boolptr;
   SCIP_DEBUGSOLDATA* debugsoldata;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(solcontained != NULL);

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   debugsoldata = SCIPsetGetDebugSolData(set);
   assert(debugsoldata != NULL);

   if( debugsoldata ->debugsoldisabled )
   {
      *solcontained = FALSE;
      return SCIP_OKAY;
   }

   /* generate the hashmap */
   if( debugsoldata->solinnode == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&debugsoldata->solinnode, blkmem, SCIP_HASHSIZE_DEBUG) );
   }

   /* check, whether we know already whether the solution is contained in the given node */
   boolptr = (SCIP_Bool*)SCIPhashmapGetImage(debugsoldata->solinnode, (void*)node);
   if( boolptr != NULL )
   {
      if( boolptr != &debugsoldata->falseptr && boolptr != &debugsoldata->trueptr )
      {
         SCIPerrorMessage("wrong value in node hashmap\n");
         SCIPABORT();
         return SCIP_ERROR;
      }
      *solcontained = *boolptr;
      return SCIP_OKAY;
   }

   /* if the solution is not contained in the parent of the node, it cannot be contained in the current node */
   *solcontained = TRUE;
   if( node->parent != NULL )
   {
      SCIP_CALL( isSolutionInNode(blkmem, set, node->parent, solcontained) );
   }

   if( *solcontained )
   {
      /* check whether the bound changes at the current node remove the debugging solution from the subproblem */
      if( node->domchg != NULL )
      {
         SCIP_DOMCHGBOUND* domchgbound;
         SCIP_BOUNDCHG* boundchgs;
         int i;

         domchgbound = &node->domchg->domchgbound;
         boundchgs = domchgbound->boundchgs;
         for( i = 0; i < (int)domchgbound->nboundchgs && *solcontained; ++i )
         {
            SCIP_Real varsol;

            /* get solution value of variable */
            SCIP_CALL( getSolutionValue(set, boundchgs[i].var, &varsol) );

            if( varsol != SCIP_UNKNOWN ) /*lint !e777*/
            {
               /* compare the bound change with the solution value */
               if( SCIPboundchgGetBoundtype(&boundchgs[i]) == SCIP_BOUNDTYPE_LOWER )
                  *solcontained = SCIPsetIsFeasGE(set, varsol, boundchgs[i].newbound);
               else
                  *solcontained = SCIPsetIsFeasLE(set, varsol, boundchgs[i].newbound);

               if( !(*solcontained) && SCIPboundchgGetBoundchgtype(&boundchgs[i]) != SCIP_BOUNDCHGTYPE_BRANCHING )
               {
                  SCIPerrorMessage("debugging solution was cut off in local node %p at depth %d by inference <%s>[%.15g] %s %.15g\n",
                     node, SCIPnodeGetDepth(node), SCIPvarGetName(boundchgs[i].var), varsol,
                     SCIPboundchgGetBoundtype(&boundchgs[i]) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", boundchgs[i].newbound);
                  SCIPABORT();
               }
            }
            else if( SCIPboundchgGetBoundchgtype(&boundchgs[i]) == SCIP_BOUNDCHGTYPE_BRANCHING )
            {
               /* we branched on a variable were we don't know the solution: no debugging can be applied in this subtree */
               *solcontained = FALSE;
            }
         }
      }
   }

   /* remember the status of the current node */
   SCIP_CALL( SCIPhashmapSetImage(debugsoldata->solinnode, (void*)node, *solcontained ? (void*)(&debugsoldata->trueptr) : (void*)(&debugsoldata->falseptr)) );

   return SCIP_OKAY;
}

/** frees the debug solution */
SCIP_RETCODE SCIPdebugFreeSol(
   SCIP_SET*             set
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;

   debugsoldata = SCIPsetGetDebugSolData(set);
   assert(debugsoldata != NULL);

   if( debugsoldata->debugsol != NULL && ((SCIPgetStage(set->scip) > SCIP_STAGE_PROBLEM && debugsoldata->debugsolstage > SCIP_STAGE_PROBLEM)
         || (SCIPgetStage(set->scip) <= SCIP_STAGE_PROBLEM && debugsoldata->debugsolstage <= SCIP_STAGE_PROBLEM)) )
   {
      SCIP_CALL( SCIPfreeSol(set->scip, &debugsoldata->debugsol) );
   }

   return SCIP_OKAY;
}

/** resets the data structure after restart */
SCIP_RETCODE SCIPdebugReset(
   SCIP_SET*             set
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;

   assert(set != NULL);

   debugsoldata = SCIPsetGetDebugSolData(set);
   assert(debugsoldata != NULL);

   if( debugsoldata->solinnode != NULL )
   {
      SCIP_CALL( SCIPhashmapRemoveAll(debugsoldata->solinnode) );
   }

   return SCIP_OKAY;
}

/** frees all debugging solution data */
SCIP_RETCODE SCIPdebugFreeDebugData(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int s;

   SCIP_DEBUGSOLDATA* debugsoldata;
   assert(set != NULL);

   debugsoldata = SCIPsetGetDebugSolData(set);
   assert(debugsoldata != NULL);

   for( s = debugsoldata->nsolvals - 1; s >= 0; --s )
      BMSfreeMemoryArrayNull(&(debugsoldata->solnames[s]));

   BMSfreeMemoryArrayNull(&debugsoldata->solnames);
   BMSfreeMemoryArrayNull(&debugsoldata->solvals);

   debugsoldata->nsolvals = 0;
   debugsoldata->debugsolval= 0.0;
   debugsoldata->solisachieved = FALSE;

   if( debugsoldata->solinnode != NULL)
      SCIPhashmapFree(&debugsoldata->solinnode);

   /* free the debug solution */
   SCIP_CALL( SCIPdebugFreeSol(set) );

   BMSfreeMemoryNull(&debugsoldata);

   set->debugsoldata = NULL;

   return SCIP_OKAY;
}

/** checks for validity of the debugging solution in given constraints */
SCIP_RETCODE SCIPdebugCheckConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to check for validity */
   int                   nconss              /**< number of given constraints */
   )
{
   SCIP_RESULT result;
   int c;

   SCIP_DEBUGSOLDATA* debugsoldata;
   assert(scip->set != NULL);

   /* check if we are in the original problem and not in a sub MIP */
   if( !SCIPdebugSolIsEnabled(scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(scip->set) )
      return SCIP_OKAY;

   debugsoldata = SCIPsetGetDebugSolData(scip->set);

   assert(conss != NULL || nconss == 0);
   assert(debugsoldata->debugsol != NULL);

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug
    * solution
    */
   if( debugSolIsAchieved(scip->set) )
      return SCIP_OKAY;

   result = SCIP_FEASIBLE;

   /* checking each given constraint against the debugging solution */
   for( c = nconss - 1; c >= 0; --c )
   {
      assert(conss[c] != NULL);

      if( !SCIPconsIsActive(conss[c]) )
         continue;

      assert(SCIPconsGetActiveDepth(conss[c]) <= SCIPgetDepth(scip));

      /* if the cons is only locally valid, check whether the debugging solution is contained in the local subproblem */
      if( SCIPconsIsLocal(conss[c]) )
      {
         SCIP_Bool solcontained;

         SCIP_CALL( isSolutionInNode(SCIPblkmem(scip), scip->set, SCIPgetCurrentNode(scip), &solcontained) );
         if( !solcontained )
            return SCIP_OKAY;
      }

      SCIP_CALL( SCIPcheckCons(scip, conss[c], debugsoldata->debugsol, TRUE, TRUE, TRUE, &result) );

      SCIPdebugMsg(scip, " -> checking of constraint %s returned result <%d>\n", SCIPconsGetName(conss[c]), result);

      if( result != SCIP_FEASIBLE )
      {
         SCIPerrorMessage("constraint %s violates the debugging solution\n", SCIPconsGetName(conss[c]));
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** checks whether given row is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< row to check for validity */
   )
{
   SCIP_COL** cols;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nnonz;
   int i;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real solval;

   assert(set != NULL);
   assert(row != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(set->scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* if the row is only locally valid, check whether the debugging solution is contained in the local subproblem */
   if( SCIProwIsLocal(row) )
   {
      SCIP_Bool solcontained;

      SCIP_CALL( isSolutionInNode(SCIPblkmem(set->scip), set, SCIPgetCurrentNode(set->scip), &solcontained) );
      if( !solcontained )
         return SCIP_OKAY;
   }

   cols = SCIProwGetCols(row);
   vals = SCIProwGetVals(row);
   nnonz = SCIProwGetNNonz(row);
   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);

   /* calculate row's activity on debugging solution */
   minactivity = SCIProwGetConstant(row);
   maxactivity = minactivity;
   for( i = 0; i < nnonz; ++i )
   {
      SCIP_VAR* var;

      /* get solution value of variable in debugging solution */
      var = SCIPcolGetVar(cols[i]);
      SCIP_CALL( getSolutionValue(set, var, &solval) );

      if( solval != SCIP_UNKNOWN ) /*lint !e777*/
      {
         minactivity += vals[i] * solval;
         maxactivity += vals[i] * solval;
      }
      else if( vals[i] > 0.0 )
      {
         minactivity += vals[i] * SCIPvarGetLbGlobal(var);
         maxactivity += vals[i] * SCIPvarGetUbGlobal(var);
      }
      else if( vals[i] < 0.0 )
      {
         minactivity += vals[i] * SCIPvarGetUbGlobal(var);
         maxactivity += vals[i] * SCIPvarGetLbGlobal(var);
      }
   }
   SCIPsetDebugMsg(set, "debugging solution on row <%s>: %g <= [%g,%g] <= %g\n",
      SCIProwGetName(row), lhs, minactivity, maxactivity, rhs);

   /* check row for violation, using absolute LP feasibility tolerance (as LP solver should do) */
   if( maxactivity + SCIPsetLpfeastol(set) < lhs || minactivity - SCIPsetLpfeastol(set) > rhs )
   {
      printf("***** debug: row <%s> violates debugging solution (lhs=%.15g, rhs=%.15g, activity=[%.15g,%.15g], local=%u, lpfeastol=%g)\n",
         SCIProwGetName(row), lhs, rhs, minactivity, maxactivity, SCIProwIsLocal(row), SCIPsetLpfeastol(set));
      SCIProwPrint(row, SCIPgetMessagehdlr(set->scip), NULL);

      /* output row with solution values */
      printf("\n\n");
      printf("***** debug: violated row <%s>:\n", SCIProwGetName(row));
      printf(" %.15g <= %.15g", lhs, SCIProwGetConstant(row));
      for( i = 0; i < nnonz; ++i )
      {
         /* get solution value of variable in debugging solution */
         SCIP_CALL( getSolutionValue(set, SCIPcolGetVar(cols[i]), &solval) );
         printf(" %+.15g<%s>[%.15g]", vals[i], SCIPvarGetName(SCIPcolGetVar(cols[i])), solval);
      }
      printf(" <= %.15g\n", rhs);

      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given global lower bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lb                  /**< lower bound */
   )
{
   SCIP_Real varsol;

   assert(scip != NULL);
   assert(var != NULL);

   /* check if we are in the original problem and not in a sub MIP */
   if( !SCIPdebugSolIsEnabled(scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(scip->set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(scip->set) )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(scip->set, var, &varsol) );
   SCIPdebugMsg(scip, "debugging solution on lower bound of <%s>[%g] >= %g\n", SCIPvarGetName(var), varsol, lb);

   /* check validity of debugging solution */
   if( varsol != SCIP_UNKNOWN && SCIPisFeasLT(scip, varsol, lb) ) /*lint !e777*/
   {
      SCIPerrorMessage("invalid global lower bound: <%s>[%.15g] >= %.15g\n", SCIPvarGetName(var), varsol, lb);
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given global upper bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   SCIP_Real varsol;

   assert(scip != NULL);
   assert(var != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(scip->set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(scip->set) )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(scip->set, var, &varsol) );
   SCIPdebugMsg(scip, "debugging solution on upper bound of <%s>[%g] <= %g\n", SCIPvarGetName(var), varsol, ub);

   /* check validity of debugging solution */
   if( varsol != SCIP_UNKNOWN && SCIPisFeasGT(scip, varsol, ub) ) /*lint !e777*/
   {
      SCIPerrorMessage("invalid global upper bound: <%s>[%.15g] <= %.15g\n", SCIPvarGetName(var), varsol, ub);
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given local bound implication is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckInference(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< local node where this bound change was applied */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   SCIP_Real varsol;
   SCIP_Bool solcontained;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(var != NULL);

   /* in case we are in probing or diving we have to avoid checking the solution */
   if( SCIPlpDiving(set->scip->lp) || SCIPtreeProbing(set->scip->tree) )
      return SCIP_OKAY;

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(set->scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(set, var, &varsol) );

   /* check validity of debugging solution */
   if( varsol != SCIP_UNKNOWN ) /*lint !e777*/
   {
      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasLT(set, varsol, newbound) )
      {
         SCIPerrorMessage("invalid local lower bound implication: <%s>[%.15g] >= %.15g\n", SCIPvarGetName(var), varsol, newbound);
         SCIPABORT();
      }
      if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasGT(set, varsol, newbound) )
      {
         SCIPerrorMessage("invalid local upper bound implication: <%s>[%.15g] <= %.15g\n", SCIPvarGetName(var), varsol, newbound);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** informs solution debugger, that the given node will be freed */
SCIP_RETCODE SCIPdebugRemoveNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node that will be freed */
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);

   debugsoldata = SCIPsetGetDebugSolData(set);
   assert(debugsoldata != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(set->scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* check if a solution will be cutoff in tree */
   if( SCIPgetStage(set->scip) != SCIP_STAGE_EXITSOLVE && SCIPgetStage(set->scip) != SCIP_STAGE_EXITPRESOLVE
      && SCIPnodeGetType(node) != SCIP_NODETYPE_PROBINGNODE && !SCIPisInRestart(set->scip) )
   {
      SCIP_Bool solisinnode;

      solisinnode = FALSE;

      SCIP_CALL( isSolutionInNode(blkmem, set, node, &solisinnode) );
      /* wrong node will be cutoff */
      if( solisinnode )
      {
         SCIPerrorMessage("debugging solution was cut off in local node #%" SCIP_LONGINT_FORMAT " (%p) at depth %d\n",
            node->number, node, SCIPnodeGetDepth(node));
         SCIPABORT();
      }
   }

   /* remove node from the hash map */
   if( debugsoldata->solinnode != NULL )
   {
      SCIP_CALL( SCIPhashmapRemove(debugsoldata->solinnode, (void*)node) );
   }

   return SCIP_OKAY;
}

/** checks whether given variable bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckVbound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable x in x <= b*z + d  or  x >= b*z + d */
   SCIP_BOUNDTYPE        vbtype,             /**< type of variable bound (LOWER or UPPER) */
   SCIP_VAR*             vbvar,              /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbcoef,             /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbconstant          /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   )
{
   SCIP_Real varsol;
   SCIP_Real vbvarsol;
   SCIP_Real vb;

   assert(set != NULL);
   assert(var != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(set->scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* get solution value of variables */
   SCIP_CALL( getSolutionValue(set, var, &varsol) );
   SCIP_CALL( getSolutionValue(set, vbvar, &vbvarsol) );

   /* check validity of debugging solution */
   if( varsol != SCIP_UNKNOWN && vbvarsol != SCIP_UNKNOWN ) /*lint !e777*/
   {
      vb = vbcoef * vbvarsol + vbconstant;
      if( (vbtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasLT(set, varsol, vb))
         || (vbtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasGT(set, varsol, vb)) )
      {
         SCIPerrorMessage("invalid variable bound: <%s>[%.15g] %s %.15g<%s>[%.15g] %+.15g\n",
            SCIPvarGetName(var), varsol, vbtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", vbcoef,
            SCIPvarGetName(vbvar), vbvarsol, vbconstant);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** checks whether given implication is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckImplic(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound           /**< bound b    in implication y <= b or y >= b */
   )
{
   SCIP_Real solval;

   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(set->scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(set, var, &solval) );
   if( solval == SCIP_UNKNOWN ) /*lint !e777*/
      return SCIP_OKAY;
   assert(SCIPsetIsFeasZero(set, solval) || SCIPsetIsFeasEQ(set, solval, 1.0));

   /* check, whether the implication applies for the debugging solution */
   if( (solval > 0.5) != varfixing )
      return SCIP_OKAY;

   /* get solution value of implied variable */
   SCIP_CALL( getSolutionValue(set, implvar, &solval) );
   if( solval == SCIP_UNKNOWN ) /*lint !e777*/
      return SCIP_OKAY;

   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      if( SCIPsetIsFeasLT(set, solval, implbound) )
      {
         SCIPerrorMessage("invalid implication <%s> == %d -> <%s> >= %.15g (variable has value %.15g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         SCIPABORT();
      }
   }
   else
   {
      if( SCIPsetIsFeasGT(set, solval, implbound) )
      {
         SCIPerrorMessage("invalid implication <%s> == %d -> <%s> <= %.15g (variable has value %.15g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** check whether given clique is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckClique(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< binary variables in the clique: at most one can be set to the given value */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars               /**< number of variables in the clique */
   )
{
   SCIP_Real solval;
   int pos1;
   int pos2;
   int v;

   assert(set != NULL);
   assert(vars != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(set->scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   pos1 = -1;
   pos2 = -1;

   for( v = 0; v < nvars; ++v )
   {
      assert(vars[v] != NULL);
      assert(SCIPvarIsBinary(vars[v]));

      /* get solution value of variable */
      SCIP_CALL( getSolutionValue(set, vars[v], &solval) );

      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         continue;

      assert(SCIPsetIsFeasZero(set, solval) || SCIPsetIsFeasEQ(set, solval, 1.0));

      /* negated solution value if negated variable is in clique */
      if( values != NULL && values[v] == 0 )
         solval = 1.0 - solval;

      if( SCIPsetIsFeasEQ(set, solval, 1.0) )
      {
         if( pos1 == -1 )
            pos1 = v;
         else
         {
            assert(pos2 == -1);
            pos2 = v;
            break;
         }
      }
   }

   /* print debug message if the clique violates the debugging solution */
   if( pos2 != -1 )
   {
      assert(pos1 != -1);
      SCIPerrorMessage("clique violates debugging solution, (at least) variable <%s%s> and variable <%s%s> are both one in the debugging solution\n",
         (values == NULL || values[pos1]) ? "" : "~", SCIPvarGetName(vars[pos1]), (values == NULL || values[pos2]) ? "" : "~", SCIPvarGetName(vars[pos2]));
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** check, whether at least one literals is TRUE in the debugging solution */
static
SCIP_Bool debugCheckBdchginfos(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change informations of the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict, or NULL */
   int                   nbdchginfos         /**< number of bound changes in the conflict set */
   )
{
   SCIP_Real solval;
   int i;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   assert(SCIPdebugSolIsEnabled(set->scip));

   solval = 0.0;
   /* check, whether at least one literals is TRUE in the debugging solution */
   for( i = 0; i < nbdchginfos; ++i )
   {
      SCIP_BDCHGINFO* bdchginfo;
      SCIP_VAR* var;
      SCIP_Real newbound;

      bdchginfo = bdchginfos[i];
      assert(bdchginfo != NULL);

      var = SCIPbdchginfoGetVar(bdchginfo);
      assert(var != NULL);

      if( relaxedbds != NULL )
         newbound = relaxedbds[i];
      else
         newbound = SCIPbdchginfoGetNewbound(bdchginfo);

      SCIP_CALL( getSolutionValue(set, var, &solval) );

      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         return TRUE;

      if( SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER )
      {
         assert(SCIPsetIsLE(set, newbound, SCIPbdchginfoGetNewbound(bdchginfo)));

         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         {
            if( SCIPsetIsLE(set, solval, newbound) )
               return TRUE;
         }
         else
         {
            if( SCIPsetIsLT(set, solval, newbound) )
               return TRUE;
         }
      }
      else
      {
         assert(SCIPsetIsGE(set, newbound, SCIPbdchginfoGetNewbound(bdchginfo)));

         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         {
            if( SCIPsetIsGE(set, solval, newbound) )
               return TRUE;
         }
         else
         {
            if( SCIPsetIsGT(set, solval, newbound) )
               return TRUE;
         }
      }
   }

   return FALSE;
}

/** print bound change information */
static
SCIP_RETCODE printBdchginfo(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO *      bdchginfo,          /**< bound change information */
   SCIP_Real             relaxedbd           /**< array with relaxed bounds which are efficient to create a valid conflict, or NULL */
   )
{
   SCIP_Real solval;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* get solution value within the debug solution */
   SCIP_CALL( getSolutionValue(set, SCIPbdchginfoGetVar(bdchginfo), &solval) );

   printf(" <%s>[%.15g] %s %g(%g)", SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)), solval,
      SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
      SCIPbdchginfoGetNewbound(bdchginfo), relaxedbd);

   return SCIP_OKAY;
}


/** print bound change information */
static
SCIP_RETCODE printBdchginfos(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change information array */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict, or NULL */
   int                   nbdchginfos         /**< number of bound changes in the conflict set */
   )
{
   int i;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   for( i = 0; i < nbdchginfos; ++i )
   {
      SCIP_BDCHGINFO* bdchginfo;

      bdchginfo = bdchginfos[i];
      assert(bdchginfo != NULL);

      printBdchginfo(set, bdchginfo, relaxedbds != NULL ? relaxedbds[i] : SCIPbdchginfoGetNewbound(bdchginfo));
   }

   return SCIP_OKAY;
}

/** checks whether given conflict is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckConflict(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node where the conflict clause is added */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change informations of the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict */
   int                   nbdchginfos         /**< number of bound changes in the conflict set */
   )
{
   SCIP_Bool solcontained;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(nbdchginfos == 0 || bdchginfos != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(set->scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* check, whether at least one literals is TRUE in the debugging solution */
   if( debugCheckBdchginfos(set, bdchginfos, relaxedbds, nbdchginfos) )
      return SCIP_OKAY;

   SCIPerrorMessage("invalid conflict set:");

   /* print bound changes which are already part of the conflict set */
   SCIP_CALL( printBdchginfos(set, bdchginfos, relaxedbds, nbdchginfos) );

   printf("\n");
   SCIPABORT();

   return SCIP_OKAY; /*lint !e527*/
}

/** checks whether given conflict graph frontier is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckConflictFrontier(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node where the conflict clause is added */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change info which got resolved, or NULL */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change informations of the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict */
   int                   nbdchginfos,        /**< number of bound changes in the conflict set */
   SCIP_PQUEUE*          bdchgqueue,         /**< unprocessed conflict bound changes */
   SCIP_PQUEUE*          forcedbdchgqueue    /**< unprocessed conflict bound changes that must be resolved */
   )
{
   SCIP_BDCHGINFO** bdchgqueued;
   SCIP_BDCHGINFO** forcedbdchgqueued;
   SCIP_Bool solcontained;
   int nbdchgqueued;
   int nforcedbdchgqueued;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(nbdchginfos == 0 || bdchginfos != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(set->scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* check, whether one literals is TRUE in the debugging solution */
   if( debugCheckBdchginfos(set, bdchginfos, relaxedbds, nbdchginfos) )
      return SCIP_OKAY;

   /* get the elements of the bound change queue */
   bdchgqueued = (SCIP_BDCHGINFO**)SCIPpqueueElems(bdchgqueue);
   nbdchgqueued = SCIPpqueueNElems(bdchgqueue);

   /* check, whether one literals is TRUE in the debugging solution */
   if( debugCheckBdchginfos(set, bdchgqueued, NULL, nbdchgqueued) )
      return SCIP_OKAY;

   /* get the elements of the bound change queue */
   forcedbdchgqueued = (SCIP_BDCHGINFO**)SCIPpqueueElems(forcedbdchgqueue);
   nforcedbdchgqueued = SCIPpqueueNElems(forcedbdchgqueue);

   /* check, whether one literals is TRUE in the debugging solution */
   if( debugCheckBdchginfos(set, forcedbdchgqueued, NULL, nforcedbdchgqueued) )
      return SCIP_OKAY;

   SCIPerrorMessage("invalid conflict frontier");

   if( bdchginfo != NULL )
   {
      printf(" (after resolving bound change ");
      printBdchginfo(set, bdchginfo, SCIPbdchginfoGetNewbound(bdchginfo));
      printf(")");
   }
   printf(":");

   /* print bound changes which are already part of the conflict set */
   SCIP_CALL( printBdchginfos(set, bdchginfos, relaxedbds, nbdchginfos) );

   /* print bound changes which are queued */
   SCIP_CALL( printBdchginfos(set, bdchgqueued, NULL, nbdchgqueued) );

   /* print bound changes which are queued in the force queue */
   SCIP_CALL( printBdchginfos(set, forcedbdchgqueued, NULL, nforcedbdchgqueued) );

   printf("\n");
   SCIPABORT();

   return SCIP_OKAY; /*lint !e527*/
}

/** check whether the debugging solution is valid in the current node */
SCIP_RETCODE SCIPdebugSolIsValidInSubtree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            isvalidinsubtree    /**< pointer to store whether the solution is valid in the current
                                              *   subtree */
   )
{
   SCIP_Bool solcontained;

   *isvalidinsubtree = FALSE;

   assert(scip->set != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( !SCIPdebugSolIsEnabled(scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(scip->set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(scip->set) )
      return SCIP_OKAY;

   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(SCIPblkmem(scip), scip->set, SCIPgetCurrentNode(scip), &solcontained) );

   if( solcontained )
      *isvalidinsubtree = TRUE;

   return SCIP_OKAY;
}

/** checks whether SCIP data structure is the main SCIP (the one for which debugging is enabled) */
SCIP_Bool SCIPdebugIsMainscip(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return SCIPdebugSolIsEnabled(scip);
}

/** enabling solution debugging mechanism */
void SCIPdebugSolEnable(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;
   assert(scip != NULL);
   assert(scip->set != NULL);

   debugsoldata = SCIPsetGetDebugSolData(scip->set);
   assert(debugsoldata != NULL);

   debugsoldata->debugsoldisabled = FALSE;
}

/** disabling solution debugging mechanism */
void SCIPdebugSolDisable(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;
   assert(scip != NULL);
   assert(scip->set != NULL);

   debugsoldata = SCIPsetGetDebugSolData(scip->set);
   assert(debugsoldata != NULL);

   debugsoldata->debugsoldisabled = TRUE;
}

/** check if solution debugging mechanism is enabled */
SCIP_Bool SCIPdebugSolIsEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;
   assert(scip != NULL);
   assert(scip->set != NULL);

   debugsoldata = SCIPsetGetDebugSolData(scip->set);
   assert(debugsoldata != NULL);

   return (!debugsoldata->debugsoldisabled);
}

/** propagator to force finding the debugging solution */
static
SCIP_DECL_PROPEXEC(propExecDebug)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* check if we are in the original problem and not in a sub MIP */
   if( !SCIPdebugIsMainscip(scip) )
      return SCIP_OKAY;

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(scip->set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(scip->set) )
      return SCIP_OKAY;

#if 1
   /* solve at least one LP */
   if( SCIPgetNLPIterations(scip) == 0 )
      return SCIP_OKAY;
#endif

   vars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      SCIP_CALL( getSolutionValue(scip->set, vars[i], &solval) );
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
      {
         SCIPerrorMessage("original variable without debugging solution value\n");
         SCIPABORT();
      }

      lb = SCIPvarGetLbGlobal(vars[i]);
      ub = SCIPvarGetUbGlobal(vars[i]);
      if( SCIPisLT(scip, solval, lb) || SCIPisGT(scip, solval, ub) )
      {
         SCIPerrorMessage("solution value %.15g of <%s> outside bounds loc=[%.15g,%.15g], glb=[%.15g,%.15g]\n",
            solval, SCIPvarGetName(vars[i]), lb, ub, SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]));
         SCIPABORT();
      }

      SCIP_CALL( SCIPfixVar(scip, vars[i], solval, &infeasible, &fixed) );
      if( infeasible )
         *result = SCIP_CUTOFF;
      else if( fixed )
         *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}

/** creates the debugging propagator and includes it in SCIP */
SCIP_RETCODE SCIPdebugIncludeProp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, "debug", "debugging propagator", 99999999, -1, FALSE,
         SCIP_PROPTIMING_ALWAYS, 99999999, 0, SCIP_PRESOLTIMING_FAST, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, propExecDebug, NULL, NULL) );

   return SCIP_OKAY;
}

/** adds a solution value for a new variable in the transformed problem that has no original counterpart
 * a value can only be set if no value has been set for this variable before
 */
SCIP_RETCODE SCIPdebugAddSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which to add a value */
   SCIP_Real             val                 /**< solution value for variable */
   )
{
   SCIP_DEBUGSOLDATA* debugsoldata;
   SCIP_Real testval;
   const char* varname;
   int i;

   assert(scip != NULL);
   assert(var != NULL);
   assert(scip->set != NULL);

   debugsoldata = SCIPsetGetDebugSolData(scip->set);
   assert(debugsoldata != NULL);

   /* assert that we are in the SCIP instance that we are debugging and not some different (subSCIP,
    * auxiliary CIP, ...)
    */
   if( !SCIPdebugSolIsEnabled(scip) )
      return SCIP_OKAY;

   /* check whether a debug solution is available */
   if( !debugSolutionAvailable(scip->set) )
      return SCIP_OKAY;

   if( debugsoldata->debugsol == NULL )
   {
      /* make sure a debug solution has been read, so we do not compare against the initial debugsolval == 0 */
      SCIP_CALL( readSolution(scip->set) );
   }

   /* allocate memory */
   if( debugsoldata->nsolvals >= debugsoldata->solsize )
   {
      debugsoldata->solsize = MAX(2*debugsoldata->solsize, debugsoldata->nsolvals+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&debugsoldata->solnames, debugsoldata->solsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&debugsoldata->solvals,  debugsoldata->solsize) );
   }
   assert(debugsoldata->nsolvals < debugsoldata->solsize);

   /* store solution value in sorted list */
   varname = SCIPvarGetName(var);
   for( i = debugsoldata->nsolvals; i > 0 && strcmp(varname, debugsoldata->solnames[i-1]) < 0; --i )
   {
      debugsoldata->solnames[i] = debugsoldata->solnames[i-1];
      debugsoldata->solvals[i]  = debugsoldata->solvals[i-1];
   }
   if( i > 0 && strcmp(varname, debugsoldata->solnames[i-1]) == 0 )
   {
      if( REALABS(debugsoldata->solvals[i-1] - val) > 1e-9 )
      {
         SCIPerrorMessage("already have stored different debugging solution value (%g) for variable <%s>, cannot store %g\n", debugsoldata->solvals[i-1], varname, val);
         return SCIP_ERROR;
      }
      else
      {
         SCIPdebugMsg(scip, "already have stored debugging solution value %g for variable <%s>, do not store same value again\n", val, varname);
         for( ; i < debugsoldata->nsolvals; ++i )
         {
            debugsoldata->solnames[i] = debugsoldata->solnames[i+1];
            debugsoldata->solvals[i]  = debugsoldata->solvals[i+1];
         }
         return SCIP_OKAY;
      }
   }

   /* insert new solution value */
   SCIP_ALLOC( BMSduplicateMemoryArray(&(debugsoldata->solnames[i]), varname, strlen(varname)+1) );
   SCIPdebugMsg(scip, "add variable <%s>: value <%g>\n", debugsoldata->solnames[i], val);
   debugsoldata->solvals[i] = val;
   debugsoldata->nsolvals++;

   /* update objective function value of debug solution */
   debugsoldata->debugsolval += debugsoldata->solvals[i] * SCIPvarGetObj(var);
   SCIPdebugMsg(scip, "Debug Solution value is now %g.\n", debugsoldata->debugsolval);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
   {
      /* add values to SCIP debug solution */
      SCIP_CALL( SCIPsetSolVal(scip, debugsoldata->debugsol, var, debugsoldata->solvals[i] ) );
   }

   /* get solution value once to produce warning if solution was cut off */
   SCIPdebugGetSolVal(scip, var, &testval);

   return SCIP_OKAY;
}

#else

/** this is a dummy method to make the SunOS gcc linker happy */
extern void SCIPdummyDebugMethodForSun(void);
void SCIPdummyDebugMethodForSun(void)
{
   return;
}

#endif


/*
 * debug method for LP interface, to check if the LP interface works correct
 */
#ifdef SCIP_DEBUG_LP_INTERFACE

/* check whether coef is the r-th row of the inverse basis matrix B^-1; this is
 * the case if( coef * B ) is the r-th unit vector */
SCIP_RETCODE SCIPdebugCheckBInvRow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< r-th row of the inverse basis matrix */
   )
{
   SCIP_Real vecval;
   SCIP_Real matrixval;
   int* basisind;
   int nrows;
   int idx;
   int i;
   int k;

   assert(scip != NULL);

   nrows = SCIPgetNLPRows(scip);

   /* get basic indices for the basic matrix B */
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );


   /* loop over the columns of B */
   for( k = 0; k < nrows; ++k )
   {
      vecval = 0.0;

      /* indices of basic columns and rows:
       * - index i >= 0 corresponds to column i,
       * - index i < 0 to row -i-1
       */
      idx = basisind[k];

      /* check if we have a slack variable; this is the case if idx < 0 */
      if( idx >= 0 )
      {
         /* loop over the rows to compute the corresponding value in the unit vector */
         for( i = 0; i < nrows; ++i )
         {
            SCIP_CALL( SCIPlpiGetCoef(scip->lp->lpi, i, idx, &matrixval) );
            vecval += coef[i] * matrixval;
         }
      }
      else
      {
         assert( idx < 0 );

         /* retransform idx
          * - index i >= 0 corresponds to column i,
          * - index i < 0 to row -i-1
          */
         idx = -idx - 1;
         assert( idx >= 0 && idx < nrows );

         /* since idx < 0 we are in the case of a slack variable, i.e., the corresponding column
            is the idx-unit vector; note that some LP solver return a -idx-unit vector */
         /*   vecval = REALABS(coef[idx]);*/
         vecval = coef[idx];
      }

      /* check if vecval fits to the r-th unit vector */
      if( k == r && !SCIPisFeasEQ(scip, vecval, 1.0) )
      {
         /* we expected a 1.0 and found something different */
         SCIPmessagePrintWarning(SCIPgetMessagehdlr(scip), "checked SCIPgetLPBInvRow() found value <%g> expected 1.0\n", vecval);
      }
      else if( k != r && !SCIPisFeasZero(scip, vecval) )
      {
         /* we expected a 0.0 and found something different */
         SCIPmessagePrintWarning(SCIPgetMessagehdlr(scip), "checked SCIPgetLPBInvRow() found value <%g> expected 0.0\n", vecval);
      }
   }

   SCIPfreeBufferArray(scip, &basisind);

   return SCIP_OKAY;
}

#endif
