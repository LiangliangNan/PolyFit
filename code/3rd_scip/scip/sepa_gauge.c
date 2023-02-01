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

/**@file   sepa_gauge.c
 * @brief  gauge separator
 * @author Felipe Serrano
 *
 * @todo should separator only be run when SCIPallColsInLP is true?
 * @todo add SCIPisStopped(scip) to the condition of time consuming loops
 * @todo check if it makes sense to implement the copy callback
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_gauge.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi.h"
#include "nlpi/nlpioracle.h"
#include "nlpi/nlpi_ipopt.h"


#define SEPA_NAME              "gauge"
#define SEPA_DESC              "gauge separator"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define VIOLATIONFAC                100 /* constraints regarded as violated when violation > VIOLATIONFAC*SCIPfeastol */
#define MAX_ITER                     75 /* maximum number of iterations for the line search */

#define DEFAULT_NLPTIMELIMIT        0.0 /**< default time limit of NLP solver; 0.0 for no limit */
#define DEFAULT_NLPITERLIM         1000 /**< default NLP iteration limit */

#define NLPFEASFAC                  1e-1/**< NLP feasibility tolerance = NLPFEASFAC * SCIP's feasibility tolerance */
#define NLPVERBOSITY                  0 /**< NLP solver verbosity */

#define INTERIOROBJVARLB           -100 /**< lower bound of the objective variable when computing interior point */
/*
 * Data structures
 */

/** side that makes a nlrow convex */
enum ConvexSide
{
   LHS = 0,                                  /**< left hand side */
   RHS = 1                                   /**< right hand side */
};
typedef enum ConvexSide CONVEXSIDE;

/** position of a point */
enum Position
{
   INTERIOR = 0,         /**< point is in the interior of the region */
   BOUNDARY = 1,         /**< point is in the boundary of the region */
   EXTERIOR = 2          /**< point is in the exterior of the region */
};
typedef enum Position POSITION;

/** separator data */
struct SCIP_SepaData
{
   SCIP_NLROW**          nlrows;             /**< stores convex nlrows */
   CONVEXSIDE*           convexsides;        /**< which sides make the nlrows convex */
   int*                  nlrowsidx;          /**< indices of nlrows that violate the current lp solution */
   int                   nnlrowsidx;         /**< total number of convex nonlinear nlrows that violate the current lp solution */
   int                   nnlrows;            /**< total number of convex nonlinear nlrows */
   int                   nlrowssize;         /**< memory allocated for nlrows, convexsides and nlrowsidx */

   SCIP_Bool             isintsolavailable;  /**< do we have an interior point available? */
   SCIP_Bool             skipsepa;           /**< whether separator should be skipped */
   SCIP_SOL*             intsol;             /**< stores interior point */

   SCIP_EXPRINT*         exprinterpreter;    /**< expression interpreter to compute gradients */

   int                   ncuts;              /**< number of cuts generated */

   /* parameters */
   SCIP_Real             nlptimelimit;       /**< time limit of NLP solver; 0.0 for no limit */
   int                   nlpiterlimit;       /**< iteration limit of NLP solver; 0 for no limit */
};

/*
 * Local methods
 */

/** stores, from the constraints represented by nlrows, the nonlinear convex ones in sepadata */
static
SCIP_RETCODE storeNonlinearConvexNlrows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW**          nlrows,             /**< nlrows from which to store convex ones */
   int                   nnlrows             /**< number of nlrows */
   )
{
   int i;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(nlrows != NULL);
   assert(nnlrows > 0);

   SCIPdebugMsg(scip, "storing convex nlrows\n");

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->nlrows), nnlrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->convexsides), nnlrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->nlrowsidx), nnlrows) );
   sepadata->nlrowssize = nnlrows;

   sepadata->nnlrows = 0;
   for( i = 0; i < nnlrows; ++i )
   {
      SCIP_NLROW* nlrow;

      nlrow = nlrows[i];
      assert(nlrow != NULL);

      /* linear case */
      if( SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_LINEAR ||
            (SCIPnlrowGetNQuadElems(nlrow) == 0 && SCIPnlrowGetExprtree(nlrow) == NULL) )
         continue;

      /* nonlinear case */
      if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)) && SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONVEX )
      {
         sepadata->convexsides[sepadata->nnlrows] = RHS;
         sepadata->nlrows[sepadata->nnlrows] = nlrow;
         ++(sepadata->nnlrows);
      }
      else if( !SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)) && SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONCAVE )
      {
         sepadata->convexsides[sepadata->nnlrows] = LHS;
         sepadata->nlrows[sepadata->nnlrows] = nlrow;
         ++(sepadata->nnlrows);
      }
   }

   return SCIP_OKAY;
}

/** computes an interior point of a convex NLP relaxation; builds the convex relaxation, modifies it to find an interior
 * point, solves it and frees it; more details in @ref sepa_gauge.h
 *
 * @note the method also counts the number of nonlinear convex constraints and if there are < 2, then the convex
 * relaxation is not interesting and the separator will not run again
 */
static
SCIP_RETCODE computeInteriorPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   SCIP_NLPIORACLE* nlpioracle;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_NLPI* nlpi;
   SCIP_HASHMAP* var2nlpiidx;
   SCIP_Real timelimit;
   SCIP_Real objvarlb;
   SCIP_Real minusone;
   SCIP_Real one;
   int nconvexnlrows;
   int iterlimit;
   int objvaridx;
   int nconss;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(!sepadata->skipsepa);

   SCIPdebugMsg(scip, "Computing interior point\n");

   /* create convex relaxation NLP */
   assert(SCIPgetNNlpis(scip) > 0);

   nlpi = SCIPgetNlpis(scip)[0];
   assert(nlpi != NULL);

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, &nlpiprob, "gauge-interiorpoint-nlp") );
   SCIP_CALL( SCIPhashmapCreate(&var2nlpiidx, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPcreateNlpiProb(scip, nlpi, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip), nlpiprob, var2nlpiidx,
            NULL, SCIPgetCutoffbound(scip), FALSE, TRUE) );

   /* add objective variable; the problem is \min t, s.t. g(x) <= t, l(x) <= 0, where g are nonlinear and l linear */
   objvaridx = nvars;
   objvarlb = INTERIOROBJVARLB;
   one = 1.0;
   SCIP_CALL( SCIPnlpiAddVars(nlpi, nlpiprob, 1, &objvarlb, NULL, NULL) );
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 1, &objvaridx, &one, 0, NULL, NULL, NULL, 0.0) );

   /* add objective variables to constraints; for this we need to get nlpi oracle to have access to number of
    * constraints and which constraints are nonlinear
    */
   /* @todo: this code is only valid when using IPOPT and needs to be changed when new NLP solvers get interfaced */
   assert(strcmp(SCIPnlpiGetName(nlpi), "ipopt") == 0);
   nlpioracle = (SCIP_NLPIORACLE *)SCIPgetNlpiOracleIpopt(nlpiprob);
   assert(nlpioracle != NULL);
   assert(SCIPnlpiOracleGetNVars(nlpioracle) == objvaridx + 1);

   minusone = -1.0;
   nconvexnlrows = 0;
   nconss = SCIPnlpiOracleGetNConstraints(nlpioracle);
   for( i = 0; i < nconss; i++ )
   {
      if( SCIPnlpiOracleGetConstraintDegree(nlpioracle, i) > 1 )
      {
         SCIP_CALL( SCIPnlpiChgLinearCoefs(nlpi, nlpiprob, i, 1, &objvaridx, &minusone) );
         ++nconvexnlrows;
      }
   }
   SCIPdebug( SCIP_CALL( SCIPnlpiOraclePrintProblem(nlpioracle, SCIPgetMessagehdlr(scip), NULL) ) );

   /* check if convex relaxation is interesting */
   if( nconvexnlrows < 2 )
   {
      SCIPdebugMsg(scip, "convex relaxation is not interesting, only %d nonlinear convex rows; abort\n", nconvexnlrows);
      sepadata->skipsepa = TRUE;
      goto CLEANUP;
   }

   /* add linear rows */
   SCIP_CALL( SCIPaddNlpiProbRows(scip, nlpi, nlpiprob, var2nlpiidx, SCIPgetLPRows(scip), SCIPgetNLPRows(scip)) );

   /* set parameters in nlpi; time and iterations limit, tolerance, verbosity; for time limit, get time limit of scip;
    * if scip doesn't have much time left, don't run separator. otherwise, timelimit is the minimum between whats left
    * for scip and the timelimit setting
    */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit <= 1.0 )
      {
         SCIPdebugMsg(scip, "skip NLP solve; no time left\n");
         sepadata->skipsepa = TRUE;
         goto CLEANUP;
      }
   }
   if( sepadata->nlptimelimit > 0.0 )
      timelimit = MIN(sepadata->nlptimelimit, timelimit);
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_TILIM, timelimit) );

   iterlimit = sepadata->nlpiterlimit > 0 ? sepadata->nlpiterlimit : INT_MAX;
   SCIP_CALL( SCIPnlpiSetIntPar(nlpi, nlpiprob, SCIP_NLPPAR_ITLIM, iterlimit) );
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, NLPFEASFAC * SCIPfeastol(scip)) );
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_RELOBJTOL, MAX(SCIPfeastol(scip), SCIPdualfeastol(scip))) ); /*lint !e666*/
   SCIP_CALL( SCIPnlpiSetIntPar(nlpi, nlpiprob, SCIP_NLPPAR_VERBLEVEL, NLPVERBOSITY) );

   /* compute interior point */
   SCIPdebugMsg(scip, "starting interior point computation\n");
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   SCIPdebugMsg(scip, "finish interior point computation\n");

#ifdef SCIP_DEBUG
   {
      SCIP_NLPSTATISTICS* nlpstatistics;

      /* get statistics */
      SCIP_CALL( SCIPnlpStatisticsCreate(SCIPblkmem(scip), &nlpstatistics) );
      SCIP_CALL( SCIPnlpiGetStatistics(nlpi, nlpiprob, nlpstatistics) );

      SCIPdebugMsg(scip, "nlpi took iters %d, time %g searching for an find interior point: solstat %d\n",
            SCIPnlpStatisticsGetNIterations(nlpstatistics), SCIPnlpStatisticsGetTotalTime(nlpstatistics),
            SCIPnlpiGetSolstat(nlpi, nlpiprob));

      SCIPnlpStatisticsFree(SCIPblkmem(scip), &nlpstatistics);
   }
#endif

   if( SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      SCIP_Real* nlpisol;

      SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &nlpisol, NULL, NULL, NULL, NULL) );

      assert(nlpisol != NULL);
      SCIPdebugMsg(scip, "NLP solved: sol found has objvalue = %g\n", nlpisol[objvaridx]);

      /* if we found an interior point store it */
      if( SCIPisFeasNegative(scip, nlpisol[objvaridx]) )
      {
         SCIPdebugMsg(scip, "Interior point found!, storing it\n");
         SCIP_CALL( SCIPcreateSol(scip, &sepadata->intsol, NULL) );
         for( i = 0; i < nvars; i ++ )
         {
            SCIP_VAR* var;

            var = SCIPgetVars(scip)[i];
            assert(SCIPhashmapExists(var2nlpiidx, (void*)var) );

            /* @todo: filter zero? */
            SCIP_CALL( SCIPsetSolVal(scip, sepadata->intsol, var,
                     nlpisol[(int)(size_t)SCIPhashmapGetImage(var2nlpiidx, (void *)var)]) );
         }

         sepadata->isintsolavailable = TRUE;
      }
      else
      {
         SCIPdebugMsg(scip, "We got a feasible point but not interior (objval: %g)\n", nlpisol[objvaridx]);
         sepadata->skipsepa = TRUE;
      }
   }
   else
   {
      SCIPdebugMsg(scip, "We couldn't get an interior point (stat: %d)\n", SCIPnlpiGetSolstat(nlpi, nlpiprob));
      sepadata->skipsepa = TRUE;
   }

CLEANUP:
   /* free memory */
   SCIPhashmapFree(&var2nlpiidx);
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );

   return SCIP_OKAY;
}


/** find whether point is in the interior, at the boundary or in the exterior of the region described by the
 * intersection of nlrows[i] <= rhs if convexsides[i] = RHS or lhs <= nlrows[i] if convexsides[i] = LHS
 * @note: point corresponds to a convex combination between the lp solution and the interior point
 */
static
SCIP_RETCODE findPointPosition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrows,             /**< nlrows defining the region */
   int*                  nlrowsidx,          /**< indices of nlrows defining the region */
   int                   nnlrowsidx,         /**< number of nlrows indices */
   CONVEXSIDE*           convexsides,        /**< sides of the nlrows involved in the region */
   SCIP_SOL*             point,              /**< point for which we want to know its position */
   POSITION*             position            /**< buffer to store position of sol */
   )
{
   int i;

   assert(scip != NULL);
   assert(nlrows != NULL);
   assert(convexsides != NULL);
   assert(nnlrowsidx > 0);
   assert(point != NULL);
   assert(position != NULL);

   *position = INTERIOR;
   for( i = 0; i < nnlrowsidx; i++ )
   {
      SCIP_NLROW* nlrow;
      SCIP_Real activity;
      CONVEXSIDE convexside;

      nlrow = nlrows[nlrowsidx[i]];
      convexside = convexsides[nlrowsidx[i]];

      /* compute activity of nlrow at point */
      SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, point, &activity) );

      if( convexside == RHS )
      {
         assert(!SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)));

         /* if nlrow <= rhs is violated, then we are in the exterior */
         if( SCIPisFeasGT(scip, activity, SCIPnlrowGetRhs(nlrow)) )
         {
            *position = EXTERIOR;
            SCIPdebugMsg(scip, "exterior because cons <%s> has activity %g. rhs: %g\n", SCIPnlrowGetName(nlrow),
                  activity, SCIPnlrowGetRhs(nlrow));
            SCIPdebug( SCIPprintNlRow(scip, nlrow, NULL) );

            return SCIP_OKAY;
         }

         /* if nlrow(point) == rhs, then we are currently at the boundary */
         if( SCIPisFeasEQ(scip, activity, SCIPnlrowGetRhs(nlrow)) )
            *position = BOUNDARY;
      }
      else
      {
         assert(!SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));
         assert(convexside == LHS);

         /* if lhs <= nlrow is violated, then we are in the exterior */
         if( SCIPisFeasLT(scip, activity, SCIPnlrowGetLhs(nlrow)) )
         {
            *position = EXTERIOR;
            return SCIP_OKAY;
         }

         /* if lhs == nlrow(point), then we are currently at the boundary */
         if( SCIPisFeasEQ(scip, activity, SCIPnlrowGetLhs(nlrow)) )
            *position = BOUNDARY;
      }
   }

   return SCIP_OKAY;
}


/** returns, in convexcomb, the convex combination
 * \f$ \lambda \f$ endpoint + (1 - \f$ lambda \f$) startpoint = startpoint + \f$ \lambda \f$ (endpoint - tosepasol) */
static
SCIP_RETCODE buildConvexCombination(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lambda,             /**< convex combination multiplier */
   SCIP_SOL*             startpoint,         /**< point corresponding to \f$ \lambda = 0 \f$ */
   SCIP_SOL*             endpoint,           /**< point corresponding to \f$ \lambda = 1 \f$ */
   SCIP_SOL*             convexcomb          /**< solution to store convex combination of intsol and tosepasol */
   )
{
   SCIP_VAR** vars;
   int        nvars;
   int        i;

   assert(scip != NULL);
   assert(startpoint != NULL);
   assert(endpoint != NULL);
   assert(convexcomb != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real val;
      SCIP_VAR* var;

      var = vars[i];
      val = lambda * SCIPgetSolVal(scip, endpoint, var) + (1.0 - lambda) * SCIPgetSolVal(scip, startpoint, var);

      if( !SCIPisZero(scip, val) )
      {
         SCIP_CALL( SCIPsetSolVal(scip, convexcomb, var, val) );
      }
      else
      {
         SCIP_CALL( SCIPsetSolVal(scip, convexcomb, var, 0.0) );
      }
   }

   return SCIP_OKAY;
}


/** performs binary search to find the point belonging to the segment [intsol, tosepasol] that intersects the boundary
 * of the region described by the intersection of nlrows[i] <= rhs if convexsides[i] = RHS or lhs <= nlrows[i] if not,
 * for i in nlrowsidx
 */
static
SCIP_RETCODE findBoundaryPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrows,             /**< nlrows defining the region */
   int*                  nlrowsidx,          /**< indices of nlrows defining the region */
   int                   nnlrowsidx,         /**< number of nlrows indices */
   CONVEXSIDE*           convexsides,        /**< sides of the nlrows involved in the region */
   SCIP_SOL*             intsol,             /**< point acting as 'interior point' */
   SCIP_SOL*             tosepasol,          /**< solution that should be separated */
   SCIP_SOL*             sol,                /**< convex combination of intsol and lpsol */
   POSITION*             position            /**< buffer to store position of sol */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   int i;

   assert(scip != NULL);
   assert(nlrows != NULL);
   assert(nlrowsidx != NULL);
   assert(convexsides != NULL);
   assert(intsol != NULL);
   assert(tosepasol != NULL);
   assert(sol != NULL);
   assert(position != NULL);

   SCIPdebugMsg(scip, "starting binary search\n");
   lb = 0.0; /* corresponds to intsol */
   ub = 1.0; /* corresponds to tosepasol */
   for( i = 0; i < MAX_ITER; i++ )
   {
      /* sol = (ub+lb)/2 * lpsol + (1 - (ub+lb)/2) * intsol */
      SCIP_CALL( buildConvexCombination(scip, (ub + lb)/2.0, intsol, tosepasol, sol) );

      /* find poisition of point: boundary, interior, exterior */
      SCIP_CALL( findPointPosition(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, sol, position) );
      SCIPdebugMsg(scip, "Position: %d, lambda: %g\n", position, (ub + lb)/2.0);

      switch( *position )
      {
         case BOUNDARY:
            SCIPdebugMsg(scip, "Done\n");
            return SCIP_OKAY;

         case INTERIOR:
            /* want to be closer to tosepasol */
            lb = (ub + lb)/2.0;
            break;

         case EXTERIOR:
            /* want to be closer to intsol */
            ub = (ub + lb)/2.0;
            break;
      }
   }
   SCIPdebugMsg(scip, "Done\n");
   return SCIP_OKAY;
}


/** computes gradient of exprtree at sol */
static
SCIP_RETCODE computeGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expressions interpreter */
   SCIP_SOL*             sol,                /**< point where we compute gradient */
   SCIP_EXPRTREE*        exprtree,           /**< exprtree for which we compute the gradient */
   SCIP_Real*            grad                /**< buffer to store the gradient */
   )
{
   SCIP_Real* x;
   SCIP_Real val;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(exprint != NULL);
   assert(sol != NULL);
   assert(exprtree != NULL);
   assert(grad != NULL);

   nvars = SCIPexprtreeGetNVars(exprtree);
   assert(nvars > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &x, nvars) );

   /* compile expression exprtree, if not done before */
   if( SCIPexprtreeGetInterpreterData(exprtree) == NULL )
   {
      SCIP_CALL( SCIPexprintCompile(exprint, exprtree) );
   }

   for( i = 0; i < nvars; ++i )
   {
      x[i] = SCIPgetSolVal(scip, sol, SCIPexprtreeGetVars(exprtree)[i]);
   }

   SCIP_CALL( SCIPexprintGrad(exprint, exprtree, x, TRUE, &val, grad) );

   /*SCIPdebug( for( i = 0; i < nvars; ++i ) printf("%e [%s]\n", grad[i], SCIPvarGetName(SCIPexprtreeGetVars(exprtree)[i])) );*/

   SCIPfreeBufferArray(scip, &x);

   return SCIP_OKAY;
}


/** computes gradient cut (linearization) of nlrow at sol */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< point used to construct gradient cut (x_0) */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_NLROW*           nlrow,              /**< constraint */
   CONVEXSIDE            convexside,         /**< whether we use rhs or lhs of nlrow */
   SCIP_ROW*             row                 /**< storage for cut */
   )
{
   SCIP_Real activity;
   SCIP_Real gradx0; /* <grad f(x_0), x_0> */
   int i;

   assert(scip != NULL);
   assert(exprint != NULL);
   assert(nlrow != NULL);
   assert(row != NULL);

   gradx0 = 0.0;

   /* an nlrow has a linear part, quadratic part and expression tree; ideally one would just build the gradient but we
    * do not know if the different parts share variables or not, so we can't just build the gradient; for this reason
    * we create the row right away and compute the gradients of each part independently and add them to the row; the
    * row takes care to add coeffs corresponding to the same variable when they appear in different parts of the nlrow
    */

   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   /* linear part */
   for( i = 0; i < SCIPnlrowGetNLinearVars(nlrow); i++ )
   {
      gradx0 += SCIPgetSolVal(scip, sol, SCIPnlrowGetLinearVars(nlrow)[i]) * SCIPnlrowGetLinearCoefs(nlrow)[i];
      SCIP_CALL( SCIPaddVarToRow(scip, row, SCIPnlrowGetLinearVars(nlrow)[i], SCIPnlrowGetLinearCoefs(nlrow)[i]) );
   }

   /* quadratic part */
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); i++ )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real grad1;
      SCIP_Real grad2;
      SCIP_Real solval1;
      SCIP_Real solval2;

      var1    = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx1];
      var2    = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx2];
      solval1 = SCIPgetSolVal(scip, sol, var1);
      solval2 = SCIPgetSolVal(scip, sol, var2);
      grad1   = SCIPnlrowGetQuadElems(nlrow)[i].coef * solval2; /* note that solval2 is correct */
      grad2   = SCIPnlrowGetQuadElems(nlrow)[i].coef * solval1;

      SCIP_CALL( SCIPaddVarToRow(scip, row, var1, grad1) );
      SCIP_CALL( SCIPaddVarToRow(scip, row, var2, grad2) );

      gradx0 += grad1 * solval1 + grad2 * solval2;
   }

   /* expression tree part */
   {
      SCIP_Real* grad;
      SCIP_EXPRTREE* tree;

      tree = SCIPnlrowGetExprtree(nlrow);

      if( tree != NULL && SCIPexprtreeGetNVars(tree) > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &grad, SCIPexprtreeGetNVars(tree)) );

         SCIP_CALL( computeGradient(scip, exprint, sol, tree, grad) );

         for( i = 0; i < SCIPexprtreeGetNVars(tree); i++ )
         {
            gradx0 +=  grad[i] * SCIPgetSolVal(scip, sol, SCIPexprtreeGetVars(tree)[i]);
            SCIP_CALL( SCIPaddVarToRow(scip, row, SCIPexprtreeGetVars(tree)[i], grad[i]) );
         }

         SCIPfreeBufferArray(scip, &grad);
      }
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

#ifdef CUT_DEBUG
   SCIPdebugMsg(scip, "gradient: ");
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
   SCIPdebugMsg(scip, "gradient dot x_0: %g\n", gradx0);
#endif

   /* gradient cut is f(x_0) - <grad f(x_0), x_0> + <grad f(x_0), x> <= rhs or >= lhs */
   SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, sol, &activity) );
   if( convexside == RHS )
   {
      assert(!SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)));
      SCIP_CALL( SCIPchgRowRhs(scip, row, SCIPnlrowGetRhs(nlrow) - activity + gradx0) );
   }
   else
   {
      assert(convexside == LHS);
      assert(!SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));
      SCIP_CALL( SCIPchgRowLhs(scip, row, SCIPnlrowGetLhs(nlrow) - activity + gradx0) );
   }

#ifdef CUT_DEBUG
   SCIPdebugMsg(scip, "gradient cut: ");
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
#endif

   return SCIP_OKAY;
}

/** tries to generate gradient cuts at the point on the segment [intsol, tosepasol] that intersecs the boundary of the
 * convex relaxation
 * -# checks that the relative interior of the segment actually intersects the boundary
 *    (this check is needed since intsol is not necessarily an interior point)
 * -# finds point on the boundary
 * -# generates gradient cut at point on the boundary
 */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SOL*             tosepasol,          /**< solution that should be separated */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_NLROW**   nlrows;
   CONVEXSIDE*    convexsides;
   SCIP_SOL*      sol;
   SCIP_SOL*      intsol;
   POSITION       position;
   int*           nlrowsidx;
   int            nnlrowsidx;
   int            i;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   intsol = sepadata->intsol;
   nlrows = sepadata->nlrows;
   nlrowsidx = sepadata->nlrowsidx;
   nnlrowsidx = sepadata->nnlrowsidx;
   convexsides = sepadata->convexsides;

   assert(intsol != NULL);
   assert(nlrows != NULL);
   assert(nlrowsidx != NULL);
   assert(nnlrowsidx > 0);
   assert(convexsides != NULL);

   /* to evaluate the nlrow one needs a solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* don't separate if, under SCIP tolerances, only a slight perturbation of the interior point in the direction of
    * tosepasol gives a point that is in the exterior */
   SCIP_CALL( buildConvexCombination(scip, VIOLATIONFAC * SCIPfeastol(scip), intsol, tosepasol, sol) );
   SCIP_CALL( findPointPosition(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, sol, &position) );

   if( position == EXTERIOR )
   {
#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "segment joining intsol and tosepasol seems to be contained in the exterior of the region, can't separate\n");
      /* move from intsol in the direction of -tosepasol to check if we are really tangent to the region */
      SCIP_CALL( buildConvexCombination(scip, -1e-3, intsol, tosepasol, sol) );
      SCIP_CALL( findPointPosition(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, sol, &position) );
      if( position == EXTERIOR )
      {
         SCIPdebugMsg(scip, "line through intsol and tosepasol is tangent to region; can't separate\n");
      }
      SCIP_CALL( findPointPosition(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, intsol, &position) );
      printf("Position of intsol is %s\n",
            position == EXTERIOR ? "exterior" : position == INTERIOR ? "interior": "boundary");
      SCIP_CALL( findPointPosition(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, tosepasol, &position) );
      printf("Position of tosepasol is %s\n",
            position == EXTERIOR ? "exterior" : position == INTERIOR ? "interior": "boundary");

      /* slightly move from intsol in the direction of +-tosepasol */
      SCIP_CALL( buildConvexCombination(scip, 1e-5, intsol, tosepasol, sol) );
      SCIP_CALL( findPointPosition(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, sol, &position) );
      printf("Position of intsol + 0.00001(tosepasol - inisol) is %s\n",
            position == EXTERIOR ? "exterior" : position == INTERIOR ? "interior": "boundary");
      SCIPdebug( SCIPprintSol(scip, sol, NULL, FALSE) );

      SCIP_CALL( buildConvexCombination(scip, -1e-5, intsol, tosepasol, sol) );
      SCIP_CALL( findPointPosition(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, sol, &position) );
      printf("Position of intsol - 0.00001(tosepasol - inisol) is %s\n",
            position == EXTERIOR ? "exterior" : position == INTERIOR ? "interior": "boundary");
      SCIPdebug( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
      *result = SCIP_DIDNOTFIND;
      goto CLEANUP;
   }

   /* find point on boundary */
   if( position != BOUNDARY )
   {
      SCIP_CALL( findBoundaryPoint(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, intsol, tosepasol, sol,
               &position) );

      /* if MAX_ITER weren't enough to find a point in the boundary we don't separate */
      if( position != BOUNDARY )
      {
         SCIPdebugMsg(scip, "couldn't find boundary point, don't separate\n");
         goto CLEANUP;
      }
   }

   /* generate cuts at sol */
   for( i = 0; i < nnlrowsidx; i++ )
   {
      SCIP_NLROW* nlrow;
      SCIP_ROW*   row;
      SCIP_Real   activity;
      CONVEXSIDE  convexside;
      char        rowname[SCIP_MAXSTRLEN];

      nlrow = nlrows[nlrowsidx[i]];
      convexside = convexsides[nlrowsidx[i]];

      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%u", SCIPnlrowGetName(nlrow), ++(sepadata->ncuts));

      /* only separate nlrows that are tight at the boundary point */
      SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, sol, &activity) );
      SCIPdebugMsg(scip, "cons <%s> at boundary point has activity: %g\n", SCIPnlrowGetName(nlrow), activity);

      if( (convexside == RHS && !SCIPisFeasEQ(scip, activity, SCIPnlrowGetRhs(nlrow)))
            || (convexside == LHS && !SCIPisFeasEQ(scip, activity, SCIPnlrowGetLhs(nlrow))) )
         continue;

      /* cut is globally valid, since we work on nlrows from the NLP built at the root node, which are globally valid */
      /* @todo: when local nlrows get supported in SCIP, one can think of recomputing the interior point */
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, rowname, -SCIPinfinity(scip), SCIPinfinity(scip),
               FALSE, FALSE , TRUE) );
      SCIP_CALL( generateCut(scip, sol, sepadata->exprinterpreter, nlrow, convexside, row) );

      /* add cut */
      SCIPdebugMsg(scip, "cut <%s> has efficacy %g\n", SCIProwGetName(row), SCIPgetCutEfficacy(scip, NULL, row));
      if( SCIPisCutEfficacious(scip, NULL, row) )
      {
         SCIP_Bool infeasible;

         SCIPdebugMsg(scip, "adding cut\n");
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
            break;
         }
         else
         {
            *result = SCIP_SEPARATED;
         }
      }

      /* release the row */
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

CLEANUP:
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeGauge)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolGauge)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);

   /* free memory and reset data */
   if( sepadata->isintsolavailable )
   {
      SCIPfreeBlockMemoryArray(scip, &sepadata->nlrowsidx, sepadata->nlrowssize);
      SCIPfreeBlockMemoryArray(scip, &sepadata->convexsides, sepadata->nlrowssize);
      SCIPfreeBlockMemoryArray(scip, &sepadata->nlrows, sepadata->nlrowssize);
      SCIP_CALL( SCIPfreeSol(scip, &sepadata->intsol) );
      SCIP_CALL( SCIPexprintFree(&sepadata->exprinterpreter) );

      sepadata->nnlrows = 0;
      sepadata->nnlrowsidx = 0;
      sepadata->nlrowssize = 0;
      sepadata->isintsolavailable = FALSE;
   }
   assert(sepadata->nnlrows == 0);
   assert(sepadata->nnlrowsidx == 0);
   assert(sepadata->nlrowssize == 0);
   assert(sepadata->isintsolavailable == FALSE);

   sepadata->skipsepa = FALSE;

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpGauge)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_SOL* lpsol;
   int i;

   assert(scip != NULL);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* do not run if there is no interesting convex relaxation (with at least two nonlinear convex constraint) */
   if( sepadata->skipsepa )
   {
      SCIPdebugMsg(scip, "not running because convex relaxation is uninteresting\n");
      return SCIP_OKAY;
   }

   /* do not run if SCIP has not constructed an NLP */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "NLP not constructed, skipping gauge separator\n");
      return SCIP_OKAY;
   }

   /* do not run if SCIP has no way of solving nonlinear problems */
   if( SCIPgetNNlpis(scip) == 0 )
   {
      SCIPdebugMsg(scip, "Skip gauge separator: no nlpi and SCIP can't solve nonlinear problems without a nlpi\n");
      return SCIP_OKAY;
   }

   /* if we don't have an interior point compute one; if we fail to compute one, then separator will not be run again;
    * otherwise, we also store the convex nlrows in sepadata
    */
   if( !sepadata->isintsolavailable )
   {
      /* @todo: one could store the convex nonlinear rows inside computeInteriorPoint */
      SCIP_CALL( computeInteriorPoint(scip, sepadata) );
      assert(sepadata->skipsepa || sepadata->isintsolavailable);

      if( sepadata->skipsepa )
         return SCIP_OKAY;

      SCIP_CALL( storeNonlinearConvexNlrows(scip, sepadata, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip)) );
      SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &sepadata->exprinterpreter) );
   }

#ifdef SCIP_DISABLED_CODE
   /* get interior point: try to compute an interior point, otherwise use primal solution, otherwise use NLP solution */
   /* @todo: - decide order:
    *        - we can also use convex combination of solutions; there is a function SCIPvarGetAvgSol!
    *        - can add an event handler to only update when a new solution has been found
    */
   if( !sepadata->isintsolavailable && sepadata->computeintsol )
   {
      SCIP_Bool success;

      success = FALSE;
      SCIP_CALL( computeInteriorPoint(scip, sepa, &success) );

      if( success )
         sepadata->isintsolavailable = TRUE;
   }

   if( !sepadata->isintsolavailable )
   {
      if( SCIPgetNSols(scip) > 0 )
      {
         SCIPdebugMsg(scip, "Using current primal solution as interior point!\n");
         SCIP_CALL( SCIPcreateSolCopy(scip, &sepadata->intsol, SCIPgetBestSol(scip)) );
         sepadata->isintsolavailable = TRUE;
      }
      else if( SCIPhasNLPSolution(scip, NULL) )
      {
         SCIPdebugMsg(scip, "Using NLP solution as interior point!\n");
         SCIP_CALL( SCIPcreateNLPSol(scip, &sepadata->intsol, NULL) );
         sepadata->isintsolavailable = TRUE;
      }
      else
      {
         SCIPdebugMsg(scip, "We couldn't find an interior point, don't have a feasible nor an NLP solution; skip separator\n");
         return SCIP_OKAY;
      }
   }
#endif

   /* store lp sol (or pseudo sol when lp is not solved) to be able to use it to compute nlrows' activities */
   SCIP_CALL( SCIPcreateCurrentSol(scip, &lpsol, NULL) );

   /* store indices of relevant constraints, ie, the ones that violate the lp sol */
   sepadata->nnlrowsidx = 0;
   for( i = 0; i < sepadata->nnlrows; i++ )
   {
      SCIP_NLROW* nlrow;
      SCIP_Real activity;

      nlrow = sepadata->nlrows[i];

      SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, lpsol, &activity) );

      if( sepadata->convexsides[i] == RHS )
      {
         assert(!SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)));

         if( activity - SCIPnlrowGetRhs(nlrow) < VIOLATIONFAC * SCIPfeastol(scip) )
            continue;
      }
      else
      {
         assert(sepadata->convexsides[i] == LHS);
         assert(!SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));

         if( SCIPnlrowGetLhs(nlrow) - activity < VIOLATIONFAC * SCIPfeastol(scip) )
            continue;
      }

      sepadata->nlrowsidx[sepadata->nnlrowsidx] = i;
      ++(sepadata->nnlrowsidx);
   }

   /* separate only if there are violated nlrows */
   SCIPdebugMsg(scip, "there are %d violated nlrows\n", sepadata->nnlrowsidx);
   if( sepadata->nnlrowsidx > 0 )
   {
      SCIP_CALL( separateCuts(scip, sepa, lpsol, result) );
   }

   /* free lpsol */
   SCIP_CALL( SCIPfreeSol(scip, &lpsol) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the gauge separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaGauge(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create gauge separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );

   /* this sets all data in sepadata to 0 */
   BMSclearMemory(sepadata);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpGauge, NULL,
         sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeGauge) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolGauge) );

   /* add gauge separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nlpiterlimit",
         "iteration limit of NLP solver; 0 for no limit",
         &sepadata->nlpiterlimit, TRUE, DEFAULT_NLPITERLIM, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/nlptimelimit",
         "time limit of NLP solver; 0.0 for no limit",
         &sepadata->nlptimelimit, TRUE, DEFAULT_NLPTIMELIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
