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

/**@file   sepa_closecuts.c
 * @brief  closecuts meta separator
 * @author Marc Pfetsch
 *
 * This separator generates a convex combination of the current LP solution and either the best
 * primal feasible solution or an interior point of the LP relaxation. If the convex combination is
 * proper, the new point is closer to the convex hull of the feasible points. The separator then
 * calls all other separators to separate this point. The idea is that in this way possibly "deeper"
 * cuts are generated. Note, however, that the new point is not a basic solution, i.e., separators
 * relying basis information, e.g., Gomory cuts, will not work.
 *
 * The other cuts are generated via the sepasol() callbacks in constraints handlers or separators.
 *
 * This separator stops after a certain number (parameter @p maxunsuccessful) of unsuccessful
 * calls. It also inhibits the separation of the ordinary LP solution if it already generated enough
 * (parameter @p sepathreshold) cuts. The convex combination is determined via the parameter @p
 * sepacombvalue.
 *
 * In general, this separator makes sense if it is expected that there will be many separation
 * rounds and many cuts will be again deleted, because they are not active after a certain number of
 * rounds. In particular, branch-and-cut algorithms for combinatorial optimization problems form
 * good candidates.
 *
 * The idea seems to be first proposed in the context of the travelling salesman problem, see@par
 *   The Traveling Salesman Problem: A Computational Study@n
 *   David L. Applegate, Robert E. Bixby, Vasek Chvatal & William J. Cook@n
 *   Princeton University Press  2006@n
 *
 * for more details. See also@par
 *   Acceleration of cutting-plane and column generation algorithms: Applications to network design.@n
 *   Walid Ben-Ameur, Jose Neto@n
 *   Networks 49(1): 3-17 (2007).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_closecuts.h"


#define SEPA_NAME              "closecuts"
#define SEPA_DESC              "closecuts meta separator"
#define SEPA_PRIORITY           1000000
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */


/* default values for parameters */
#define SCIP_DEFAULT_SEPARELINT              TRUE /**< generate close cuts w.r.t. relative interior point (best solution otherwise)? */
#define SCIP_DEFAULT_SEPACOMBVALUE           0.30 /**< convex combination value for close cuts */
#define SCIP_DEFAULT_SEPATHRESHOLD             50 /**< threshold on number of generated cuts below which the ordinary separation is started */
#define SCIP_DEFAULT_INCLOBJCUTOFF          FALSE /**< include the objective cutoff when computing the relative interior? */
#define SCIP_DEFAULT_RECOMPUTERELINT        FALSE /**< recompute relative interior in each separation call? */
#define SCIP_DEFAULT_MAXUNSUCCESSFUL            0 /**< turn off separation in current node after unsuccessful calls (-1 never turn off) */
#define SCIP_DEFAULT_MAXLPITERFACTOR         10.0 /**< factor for maximal LP iterations in relative interior computation compared to node LP iterations */

#define SCIP_MIN_LPITERS                      100 /**< minimum number of allowed LP iterations in relative interior computation */


/** separator data */
struct SCIP_SepaData
{
   SCIP_Bool             separelint;         /**< generate close cuts w.r.t. relative interior point (best solution otherwise)? */
   SCIP_Bool             triedRelint;        /**< tried to compute relative interior */
   SCIP_Real             sepacombvalue;      /**< convex combination value for close cuts */
   int                   sepathreshold;      /**< threshold on number of generated cuts below which the ordinary separation is started */
   SCIP_Bool             inclobjcutoff;      /**< include the objective cutoff when computing the relative interior? */
   SCIP_Bool             recomputerelint;    /**< recompute relative interior in each separation call? */
   int                   maxunsuccessful;    /**< turn off separation in current node after unsuccessful calls (-1 never turn off) */
   SCIP_SOL*             sepasol;            /**< solution that can be used for generating close cuts */
   SCIP_Longint          discardnode;        /**< number of node for which separation is discarded */
   SCIP_Real             maxlpiterfactor;    /**< factor for maximal LP iterations in relative interior computation compared to node LP iterations */
   int                   nunsuccessful;      /**< number of consecutive unsuccessful calls */
};


/** generate point for close cut separation
 *
 *  The constructed point is the convex combination of the point stored in set->closesol and the
 *  current LP solution. The convexity parameter is set->sepa_closecombvalue. If this parameter is
 *  0, the point coincides with the LP solution.
 */
static
SCIP_RETCODE generateCloseCutPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL**            point               /**< point to be generated (or NULL if unsuccessful) */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real alpha;
   SCIP_Real onealpha;
   SCIP_Real lb;
   SCIP_Real ub;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( point != NULL );

   *point = NULL;
   if ( sepadata->sepasol == NULL )
      return SCIP_OKAY;

   alpha = sepadata->sepacombvalue;
   if ( alpha < 0.001 )
      return SCIP_OKAY;
   onealpha = 1.0 - alpha;

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, point, NULL) );

   /* generate convex combination */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for (i = 0; i < nvars; ++i)
   {
      var = vars[i];
      val = alpha * SCIPgetSolVal(scip, sepadata->sepasol, var) + onealpha * SCIPvarGetLPSol(var);

      /* If both the LP relaxation and the base point respect the variable bounds, the computed point will satisfy them
       * as well. However, variables might be fixed (e.g. by branching) since the time of the computation of the base
       * point. Thus, we adapt the value to lie inside the bounds in optimized mode. */
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      val = MAX(val, lb);
      val = MIN(val, ub);

      if ( ! SCIPisZero(scip, val) )
      {
         SCIP_CALL( SCIPsetSolVal(scip, *point, var, val) );
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */


/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyClosecuts)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaClosecuts(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeClosecuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolClosecuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   if ( sepadata->separelint && sepadata->sepasol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &sepadata->sepasol) );
      sepadata->triedRelint = FALSE;
   }

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpClosecuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_Longint currentnodenumber;
   SCIP_SOL* point = NULL;
   SCIP_Bool isroot;

   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* only call separator, if LP has been solved (need LP to compute separation point) */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* exit if we stopped ... */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* get separation data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   /* exit if we already decided to discard the current node */
   currentnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
   if ( sepadata->discardnode == currentnodenumber )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Separation method of closecuts separator.\n");

   /* check whether we have to compute a relative interior point */
   if ( sepadata->separelint )
   {
      if ( sepadata->recomputerelint )
      {
         /* check if previous relative interior point should be forgotten, otherwise it is computed only once and the
          * same point is used for all nodes */
         if ( sepadata->sepasol != NULL )
         {
            SCIP_CALL( SCIPfreeSol(scip, &sepadata->sepasol) );
            sepadata->triedRelint = FALSE;
         }
      }
      else
      {
         /* skip execution, if we unsuccessfully tried to compute a relative interior point */
         if ( sepadata->sepasol == NULL && sepadata->triedRelint )
            return SCIP_OKAY;
      }

      /* if relative interior point is not available ... */
      if ( sepadata->sepasol == NULL )
      {
         SCIP_Longint nlpiters;
         SCIP_Real timelimit;
         int iterlimit;

         /* prepare time limit */
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if ( ! SCIPisInfinity(scip, timelimit) )
            timelimit -= SCIPgetSolvingTime(scip);
         /* exit if no time left */
         if ( timelimit <= 0.0 )
            return SCIP_OKAY;

         /* determine iteration limit */
         if ( sepadata->maxlpiterfactor < 0.0 || SCIPisInfinity(scip, sepadata->maxlpiterfactor) )
            iterlimit = INT_MAX;
         else
         {
            /* determine iteration limit; the number of iterations in the root is only set after its solution, but the
             * total number of LP iterations is always updated. */
            if ( SCIPgetDepth(scip) == 0 )
               nlpiters = SCIPgetNLPIterations(scip);
            else
               nlpiters = SCIPgetNRootLPIterations(scip);
            iterlimit = (int)(sepadata->maxlpiterfactor * nlpiters);
            iterlimit = MAX(iterlimit, SCIP_MIN_LPITERS);
            assert(iterlimit > 0);
         }

         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Computing relative interior point (time limit: %g, iter limit: %d) ...\n", timelimit, iterlimit);
         SCIP_CALL( SCIPcomputeLPRelIntPoint(scip, TRUE, sepadata->inclobjcutoff, timelimit, iterlimit, &sepadata->sepasol) );
         sepadata->triedRelint = TRUE;
      }
   }
   else
   {
      /* get best solution (NULL if not present) */
      sepadata->sepasol = SCIPgetBestSol(scip);
   }

   /* separate close cuts */
   if ( sepadata->sepasol != NULL )
   {
      SCIPdebugMsg(scip, "Generating close cuts ... (combination value: %f)\n", sepadata->sepacombvalue);
      *result = SCIP_DIDNOTFIND;

      /* generate point to be separated */
      SCIP_CALL( generateCloseCutPoint(scip, sepadata, &point) );

      /* apply a separation round to generated point */
      if ( point != NULL )
      {
         int noldcuts;
         SCIP_Bool delayed;
         SCIP_Bool cutoff;

         noldcuts = SCIPgetNCuts(scip);
         isroot = (SCIP_Bool) (SCIPgetDepth(scip) == 0);

         /* separate solution via other separators */
         SCIP_CALL( SCIPseparateSol(scip, point, isroot, TRUE, FALSE, &delayed, &cutoff) );

         SCIP_CALL( SCIPfreeSol(scip, &point) );
         assert( point == NULL );

         /* the cuts might not violated by the current LP if the computed point is strange */
         SCIP_CALL( SCIPremoveInefficaciousCuts(scip) );

         if ( cutoff )
            *result = SCIP_CUTOFF;
         else
         {
            if ( SCIPgetNCuts(scip) - noldcuts > sepadata->sepathreshold )
            {
               sepadata->nunsuccessful = 0;
               *result = SCIP_NEWROUND;
            }
            else
            {
               if ( SCIPgetNCuts(scip) > noldcuts )
               {
                  sepadata->nunsuccessful = 0;
                  *result = SCIP_SEPARATED;
               }
               else
                  ++sepadata->nunsuccessful;
            }
         }

         SCIPdebugMsg(scip, "Separated close cuts: %d (enoughcuts: %d, unsuccessful: %d).\n", SCIPgetNCuts(scip) - noldcuts,
            SCIPgetNCuts(scip) - noldcuts > sepadata->sepathreshold, sepadata->nunsuccessful);

         if ( sepadata->maxunsuccessful >= 0 && sepadata->nunsuccessful > sepadata->maxunsuccessful )
         {
            SCIPdebugMsg(scip, "Turn off close cut separation, because of %d unsuccessful calls.\n", sepadata->nunsuccessful);
            sepadata->discardnode = currentnodenumber;
            sepadata->nunsuccessful = 0;
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the closecuts separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaClosecuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create closecuts separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepadata->sepasol = NULL;
   sepadata->discardnode = -1;
   sepadata->nunsuccessful = 0;
   sepadata->triedRelint = FALSE;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpClosecuts, NULL, sepadata) );

   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyClosecuts) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeClosecuts) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolClosecuts) );

   /* add closecuts separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/closecuts/separelint",
         "generate close cuts w.r.t. relative interior point (best solution otherwise)?",
         &sepadata->separelint, TRUE, SCIP_DEFAULT_SEPARELINT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/closecuts/sepacombvalue",
         "convex combination value for close cuts",
         &sepadata->sepacombvalue, TRUE, SCIP_DEFAULT_SEPACOMBVALUE, 0.0, 1.0,
         NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/closecuts/closethres",
         "threshold on number of generated cuts below which the ordinary separation is started",
         &sepadata->sepathreshold, TRUE, SCIP_DEFAULT_SEPATHRESHOLD, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/closecuts/inclobjcutoff",
         "include an objective cutoff when computing the relative interior?",
         &sepadata->inclobjcutoff, TRUE, SCIP_DEFAULT_INCLOBJCUTOFF, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/closecuts/recomputerelint",
         "recompute relative interior point in each separation call?",
         &sepadata->recomputerelint, TRUE, SCIP_DEFAULT_RECOMPUTERELINT, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/closecuts/maxunsuccessful",
         "turn off separation in current node after unsuccessful calls (-1 never turn off)",
         &sepadata->maxunsuccessful, TRUE, SCIP_DEFAULT_MAXUNSUCCESSFUL, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/closecuts/maxlpiterfactor",
         "factor for maximal LP iterations in relative interior computation compared to node LP iterations (negative for no limit)",
         &sepadata->maxlpiterfactor, TRUE, SCIP_DEFAULT_MAXLPITERFACTOR, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** sets point to be used as base point for computing the point to be separated
 *
 *  The point is only stored if separation of relative interior points is used. The solution is copied.
 */
SCIP_RETCODE SCIPsetBasePointClosecuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< base point solution */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert( scip != NULL );

   /* find separator */
   sepa = SCIPfindSepa(scip, SEPA_NAME);
   if ( sepa == NULL )
   {
      SCIPerrorMessage("Could not find separator <%s>.\n", SEPA_NAME);
      return SCIP_PLUGINNOTFOUND;
   }
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* get sepadata */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   /* store point if we have to separate relative interior points */
   if ( sepadata->separelint )
   {
      /* possibly free solution */
      if ( sepadata->sepasol != NULL )
      {
         SCIP_CALL( SCIPfreeSol(scip, &sepadata->sepasol) );
      }

      /* copy and store solution */
      SCIP_CALL( SCIPcreateSolCopy(scip, &sepadata->sepasol, sol) );
      sepadata->triedRelint = TRUE;
   }

   return SCIP_OKAY;
}
