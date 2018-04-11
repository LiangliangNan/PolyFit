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

/**@file   presol_convertinttobin.c
 * @ingroup PRESOLVERS
 * @brief  presolver that converts integer variables to binaries
 * @author Michael Winkler
 *
 *  Converts integer variables at the beginning of Presolving into their binary representation. If necessary adds a
 *  bounding knapsack constraint.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_convertinttobin.h"
#include "scip/cons_knapsack.h"
#include "scip/pub_misc.h"


#define PRESOL_NAME            "convertinttobin"
#define PRESOL_DESC            "converts integer variables to binaries"
#define PRESOL_PRIORITY        +6000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              0 /**< maximal number of presolving rounds the presolver participates in (-1: no
					 *   limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_MAXDOMAINSIZE  SCIP_LONGINT_MAX   /**< absolute value of maximum domain size which will be converted */
#define DEFAULT_ONLYPOWERSOFTWO           FALSE   /**< should only integer variables with a domain size of 2^p - 1 be
                                                    *   converted(, there we don't need an knapsack-constraint) */
#define DEFAULT_SAMELOCKSINBOTHDIRECTIONS FALSE   /**< should only integer variables with uplocks equals downlocks be converted */

/** presolver data */
struct SCIP_PresolData
{
   SCIP_Longint          maxdomainsize;      /**< absolute value of maximum domain size */
   SCIP_Bool             onlypoweroftwo;     /**< should only integer variables with a domain size of 2^p - 1 be converted */
   SCIP_Bool             samelocksinbothdirections; /**< should only integer variables with uplocks equals downlocks be converted */
};

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyConvertinttobin)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolConvertinttobin(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeConvertinttobin)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** presolving execution method */
static
SCIP_DECL_PRESOLEXEC(presolExecConvertinttobin)
{  /*lint --e{715}*/
   SCIP_VAR** scipvars;
   SCIP_VAR** vars;
   SCIP_PRESOLDATA* presoldata;
   int nbinvars;
   int nintvars;
   int v;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get the problem variables */
   scipvars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   nintvars = SCIPgetNIntVars(scip);
   if( nintvars == 0 )
      return SCIP_OKAY;

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   *result = SCIP_DIDNOTFIND;

   /* copy the integer variables into an own array, since adding binary variables affects the left-most slots in the
    * array and thereby interferes with our search loop
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, &scipvars[nbinvars], nintvars) );

   /* scan the integer variables for possible conversion into binaries;
    * we have to collect the variables first in an own
    */
   for( v = 0; v < nintvars; ++v )
   {
      SCIP_VAR** newbinvars;
      SCIP_Real* newbinvarcoeffs;
      SCIP_Longint* weights;
      SCIP_CONS* newcons;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Longint domainsize;
      char newbinvarname[SCIP_MAXSTRLEN];
      char newconsname[SCIP_MAXSTRLEN];
      int nnewbinvars;
      int v2;
      SCIP_Longint scalar;
      SCIP_Bool infeasible;
      SCIP_Bool aggregated;
      SCIP_Bool noconsknapsack;

      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER);

      /* skip variables which cannot be multi-aggregated */
      if( SCIPdoNotMultaggrVar(scip, vars[v]) )
	 continue;

      /* check for correct locks */
      if( presoldata->samelocksinbothdirections && SCIPvarGetNLocksUp(vars[v]) != SCIPvarGetNLocksDown(vars[v]) )
         continue;

      /* get variable's bounds */
      lb = SCIPvarGetLbGlobal(vars[v]);
      ub = SCIPvarGetUbGlobal(vars[v]);
      assert( SCIPisIntegral(scip, lb) );
      assert( SCIPisIntegral(scip, ub) );

      if( SCIPisInfinity(scip, ub - lb) )
         domainsize = SCIP_LONGINT_MAX;
      else
         domainsize = (SCIP_Longint) SCIPceil(scip, ub - lb);

      assert(domainsize >= 0);

      /* check for allowed domainsize */
      if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) || domainsize > presoldata->maxdomainsize )
         continue;

      /* check for domainsize is not 2^p - 1 if necessary */
      if( presoldata->onlypoweroftwo )
      {
         /* stop if domainsize is not 2^p - 1*/
         SCIP_Longint tmp;

         assert(domainsize < SCIP_LONGINT_MAX);
         tmp = domainsize + 1;

         while( tmp % 2 == 0 )
            tmp /= 2;
         if( tmp != 1 )
            continue;
      }

      noconsknapsack = FALSE;

      nnewbinvars = (int)SCIPfloor(scip, (log((SCIP_Real) domainsize)/log(2.0))) + 1;

      SCIPdebugMsg(scip, "integer variable <%s> [%g,%g], domainsize %" SCIP_LONGINT_FORMAT "\n, <uplocks = %d, downlocks = %d will be 'binarized' by %d binary variables\n ",
         SCIPvarGetName(vars[v]), lb, ub, domainsize, SCIPvarGetNLocksUp(vars[v]), SCIPvarGetNLocksDown(vars[v]), nnewbinvars);

      assert(nnewbinvars > 0);

      scalar = (SCIP_Longint)pow(2.0, nnewbinvars); /*lint !e747*/
      /* because of rounding errors */
      if( scalar == domainsize )
      {
         scalar *= 2;
         nnewbinvars++;
      }
      else if( scalar == domainsize + 1 )
         noconsknapsack = TRUE;

      assert(scalar > domainsize);

      SCIP_CALL( SCIPallocBufferArray(scip, &newbinvars, nnewbinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newbinvarcoeffs, nnewbinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &weights, nnewbinvars) );

      for( v2 = nnewbinvars - 1; v2 >= 0; --v2 )
      {
         SCIPdebugMsg(scip, "creating for <%s>[%g,%g] %d. binary variable\n", SCIPvarGetName(vars[v]), lb, ub, v2);

         /* create binary variable */
         (void) SCIPsnprintf(newbinvarname, SCIP_MAXSTRLEN, "%s_bin_%d", SCIPvarGetName(vars[v]), v2);
         SCIP_CALL( SCIPcreateVar(scip, &newbinvars[v2], newbinvarname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
               SCIPvarIsInitial(vars[v]), SCIPvarIsRemovable(vars[v]), NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, newbinvars[v2]) );

         scalar /= 2;
         assert(scalar > 0);

         newbinvarcoeffs[v2] = (SCIP_Real)scalar;
         weights[v2] = scalar;
      }

      /* aggregate integer and binary variable */
      SCIP_CALL( SCIPmultiaggregateVar(scip, vars[v], nnewbinvars, newbinvars, (SCIP_Real*)newbinvarcoeffs, lb, &infeasible, &aggregated) );
      assert(!infeasible);
      assert(aggregated);

      (void) SCIPsnprintf(newconsname, SCIP_MAXSTRLEN, "%s_bin_knapsack", SCIPvarGetName(vars[v]));

      if( !noconsknapsack )
      {
         int nodd;
         nodd = 0;
         while( domainsize % 2 == 1 )
         {
            nodd++;
            domainsize = (domainsize - 1) / 2;
         }
         if( nodd > 0 )
         {
            SCIP_Longint divisor;

            divisor = (SCIP_Longint)pow(2.0, nodd); /*lint !e747*/
            assert(divisor >= 2);

            for( v2 = nodd; v2 < nnewbinvars; ++v2 )
            {
               weights[v2] /= divisor;
            }
         }

         SCIP_CALL( SCIPcreateConsKnapsack(scip, &newcons, newconsname, nnewbinvars - nodd, &newbinvars[nodd],
               &weights[nodd], domainsize,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugPrintCons(scip, newcons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      }

      for( v2 = nnewbinvars - 1; v2 >= 0; --v2 )
      {
         /* release binary variable */
         SCIP_CALL( SCIPreleaseVar(scip, &newbinvars[v2]) );
         (*nchgvartypes)++;
      }

      SCIPfreeBufferArray(scip, &newbinvars);
      SCIPfreeBufferArray(scip, &newbinvarcoeffs);
      SCIPfreeBufferArray(scip, &weights);

      if( aggregated ) /*lint !e774*/
         *result = SCIP_SUCCESS;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the convertinttobin presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolConvertinttobin(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presolptr;

   /* create convertinttobin presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   presoldata->maxdomainsize = DEFAULT_MAXDOMAINSIZE;
   presoldata->onlypoweroftwo = DEFAULT_ONLYPOWERSOFTWO;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presolptr, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecConvertinttobin,
         presoldata) );
   assert(presolptr != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presolptr, presolCopyConvertinttobin) );
   SCIP_CALL( SCIPsetPresolFree(scip, presolptr, presolFreeConvertinttobin) );

   /* add convertinttobin presolver parameters */
   SCIP_CALL( SCIPaddLongintParam(scip,
         "presolving/" PRESOL_NAME "/maxdomainsize",
         "absolute value of maximum domain size for converting an integer variable to binaries variables",
         &presoldata->maxdomainsize, TRUE, DEFAULT_MAXDOMAINSIZE, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   /* add convertinttobin presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/onlypoweroftwo",
         "should only integer variables with a domain size of 2^p - 1 be converted(, there we don't need an knapsack-constraint for restricting the sum of the binaries)",
         &presoldata->onlypoweroftwo, TRUE, DEFAULT_ONLYPOWERSOFTWO, NULL, NULL) );

   /* add convertinttobin presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/samelocksinbothdirections",
         "should only integer variables with uplocks equals downlocks be converted",
         &presoldata->samelocksinbothdirections, TRUE, DEFAULT_SAMELOCKSINBOTHDIRECTIONS, NULL, NULL) );

   return SCIP_OKAY;
}
