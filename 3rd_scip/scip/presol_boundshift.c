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

/**@file   presol_boundshift.c
 * @brief  presolver that converts variables with domain [a,b] to variables with domain [0,b-a]
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/**@todo test this presolving step to decide whether to turn it in default mode on or off */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_boundshift.h"


#define PRESOL_NAME            "boundshift"
#define PRESOL_DESC            "converts variables with domain [a,b] to variables with domain [0,b-a]"
#define PRESOL_PRIORITY         7900000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              0 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolver (fast, medium, or exhaustive) */

#define MAXABSBOUND             1000.0  /**< maximum absolute variable bounds for aggregation */

/*
 * Default parameter settings
 */

#define DEFAULT_MAXSHIFT      SCIP_LONGINT_MAX  /**< absolute value of maximum shift */
#define DEFAULT_FLIPPING                  TRUE  /**< is flipping allowed? */
#define DEFAULT_INTEGER                   TRUE  /**< are only integer ranges shifted  */

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   SCIP_Longint          maxshift;           /**< absolute value of maximum shift */
   SCIP_Bool             flipping;           /**< is flipping allowed? */
   SCIP_Bool             integer;            /**< shift only integer values? */
};


/*
 * Local methods
 */

/** initializes the presolver data */
static
void initPresoldata(
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(presoldata != NULL);

   presoldata->maxshift = DEFAULT_MAXSHIFT;
   presoldata->flipping = DEFAULT_FLIPPING;
   presoldata->integer = DEFAULT_INTEGER;
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyBoundshift)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolBoundshift(scip) );

   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
/**! [SnippetPresolFreeBoundshift] */
static
SCIP_DECL_PRESOLFREE(presolFreeBoundshift)
{  /*lint --e{715}*/   
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}
/**! [SnippetPresolFreeBoundshift] */


/** presolving execution method */
static
SCIP_DECL_PRESOLEXEC(presolExecBoundshift)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_VAR** scipvars;
   SCIP_VAR** vars;
   int nbinvars;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* get the problem variables */
   scipvars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   nvars = SCIPgetNVars(scip) - nbinvars;

   if( nvars == 0 )
      return SCIP_OKAY;

   if( SCIPdoNotAggr(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* copy the integer variables into an own array, since adding new integer variables affects the left-most slots in
    * the array and thereby interferes with our search loop
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, &scipvars[nbinvars], nvars) );

   /* scan the integer, implicit, and continuous variables for possible conversion */
   for( v = nvars - 1; v >= 0; --v )
   {
      SCIP_VAR* var = vars[v];
      SCIP_Real lb;
      SCIP_Real ub;

      assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);

      /* get current variable's bounds */
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      /* it can happen that the variable bounds of integer variables have not been propagated yet or contain
       * some small noise; this will result in an aggregation that might trigger assertions when updating bounds of
       * aggregated variables (see #1817)
       */
      if( SCIPvarIsIntegral(var) )
      {
         assert(SCIPisIntegral(scip, lb));
         assert(SCIPisIntegral(scip, ub));

         lb = SCIPadjustedVarLb(scip, var, lb);
         ub = SCIPadjustedVarUb(scip, var, ub);
      }

      assert( SCIPisLE(scip, lb, ub) );
      if( SCIPisEQ(scip, lb, ub) )
         continue;
      if( presoldata->integer && !SCIPisIntegral(scip, ub - lb) ) 
         continue;

      /* check if bounds are shiftable */
      if( !SCIPisEQ(scip, lb, 0.0) &&                           /* lower bound != 0.0 */
         SCIPisLT(scip, ub, SCIPinfinity(scip)) &&              /* upper bound != infinity */
         SCIPisGT(scip, lb, -SCIPinfinity(scip)) &&             /* lower bound != -infinity */
#if 0
         SCIPisLT(scip, ub - lb, SCIPinfinity(scip)) &&         /* interval length less than SCIPinfinity(scip) */
#endif
         SCIPisLT(scip, ub - lb, (SCIP_Real) presoldata->maxshift) &&      /* less than max shifting */
         SCIPisLE(scip, REALABS(lb), MAXABSBOUND) &&            /* ensures a small constant in aggregation */
         SCIPisLE(scip, REALABS(ub), MAXABSBOUND) )             /* ensures a small constant in aggregation */
      {
         SCIP_VAR* newvar;
         char newvarname[SCIP_MAXSTRLEN];
         SCIP_Bool infeasible;
         SCIP_Bool redundant;
         SCIP_Bool aggregated;

         SCIPdebugMsg(scip, "convert range <%s>[%g,%g] to [%g,%g]\n", SCIPvarGetName(var), lb, ub, 0.0, (ub - lb) );

         /* create new variable */
         (void) SCIPsnprintf(newvarname, SCIP_MAXSTRLEN, "%s_shift", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateVar(scip, &newvar, newvarname, 0.0, (ub - lb), 0.0, SCIPvarGetType(var),
               SCIPvarIsInitial(var), SCIPvarIsRemovable(var), NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, newvar) );

         /* aggregate old variable with new variable */
         if( presoldata->flipping )
         {
            if( REALABS(ub) < REALABS(lb) )
            {
               SCIP_CALL( SCIPaggregateVars(scip, var, newvar, 1.0, 1.0, ub, &infeasible, &redundant, &aggregated) );
            }
            else
            {
               SCIP_CALL( SCIPaggregateVars(scip, var, newvar, 1.0, -1.0, lb, &infeasible, &redundant, &aggregated) );
            }
         }
         else
         {
            SCIP_CALL( SCIPaggregateVars(scip, var, newvar, 1.0, -1.0, lb, &infeasible, &redundant, &aggregated) );
         }

         assert(!infeasible);
         assert(redundant);
         assert(aggregated);
         SCIPdebugMsg(scip, "var <%s> with bounds [%f,%f] has obj %f\n",
            SCIPvarGetName(newvar),SCIPvarGetLbGlobal(newvar),SCIPvarGetUbGlobal(newvar),SCIPvarGetObj(newvar));

         /* release variable */
         SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

         /* take care of statistic */
         (*naggrvars)++;
         *result = SCIP_SUCCESS;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the boundshift presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolBoundshift(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presolptr;

   /* create boundshift presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );
   initPresoldata(presoldata);

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presolptr, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecBoundshift,
         presoldata) );

   assert(presolptr != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presolptr, presolCopyBoundshift) );
   SCIP_CALL( SCIPsetPresolFree(scip, presolptr, presolFreeBoundshift) );

   /* add probing presolver parameters */
   SCIP_CALL( SCIPaddLongintParam(scip,
         "presolving/boundshift/maxshift", 
         "absolute value of maximum shift",
         &presoldata->maxshift, TRUE, DEFAULT_MAXSHIFT, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/boundshift/flipping", 
         "is flipping allowed (multiplying with -1)?",
         &presoldata->flipping, TRUE, DEFAULT_FLIPPING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/boundshift/integer", 
         "shift only integer ranges?",
         &presoldata->integer, TRUE, DEFAULT_INTEGER, NULL, NULL) );

   return SCIP_OKAY;
}
