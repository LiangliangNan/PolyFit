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

/**@file   sepa_convexproj.c
 * @brief  convexproj separator
 * @author Felipe Serrano
 *
 * @todo should separator only be run when SCIPallColsInLP is true?
 * @todo check if it makes sense to implement the copy callback
 * @todo add SCIPisStopped(scip) to the condition of time consuming loops
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_convexproj.h"
#include "scip/nlp.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi.h"


#define SEPA_NAME              "convexproj"
#define SEPA_DESC              "separate at projection of point onto convex region"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE      /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                 TRUE      /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXDEPTH             -1      /* maximum depth at which the separator is applied; -1 means no limit */
#define DEFAULT_NLPTIMELIMIT        0.0      /**< default time limit of NLP solver; 0.0 for no limit */
#define DEFAULT_NLPITERLIM          250      /**< default NLP iteration limit */

#define VIOLATIONFAC                100      /* points regarded violated if max violation > VIOLATIONFAC*SCIPfeastol */

#define NLPVERBOSITY                  0      /**< NLP solver verbosity */

/*
 * Data structures
 */

/** side that makes an nlrow convex */
enum ConvexSide
{
   LHS = 0,                                  /**< left hand side */
   RHS = 1                                   /**< right hand side */
};
typedef enum ConvexSide CONVEXSIDE;

/** separator data
 * it keeps the nlpi which represents the projection problem (see sepa_convexproj.h); it also keeps the convex nlrows
 * and the side which actually makes them convex; when separating, we use the nlpi to compute the projection and then
 * the convex nlrows to compute the actual gradient cuts */
struct SCIP_SepaData
{
   SCIP_NLPI*            nlpi;               /**< nlpi used to create the nlpi problem */
   SCIP_NLPIPROBLEM*     nlpiprob;           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_VAR**            nlpivars;           /**< array containing all variables of the nlpi */
   SCIP_HASHMAP*         var2nlpiidx;        /**< mapping between variables and nlpi indices */
   int                   nlpinvars;          /**< total number of nlpi variables */

   SCIP_Bool             skipsepa;           /**< should separator be skipped? */

   SCIP_NLROW**          nlrows;             /**< convex nlrows */
   CONVEXSIDE*           convexsides;        /**< which sides make the nlrows convex */
   SCIP_Real*            constraintviolation;/**< array storing the violation of constraint by current solution; 0.0 if it is not violated */
   int                   nnlrows;            /**< total number of nlrows */
   int                   nlrowssize;         /**< memory allocated for nlrows, convexsides and nlrowsidx */

   SCIP_EXPRINT*         exprinterpreter;    /**< expression interpreter to compute gradients */

   /* parameter */
   SCIP_Real             nlptimelimit;       /**< time limit of NLP solver; 0.0 for no limit */
   int                   nlpiterlimit;       /**< iteration limit of NLP solver; 0 for no limit */
   int                   maxdepth;           /**< maximal depth at which the separator is applied */

   int                   ncuts;              /**< number of cuts generated */
};


/*
 * Local methods
 */

/** clears the sepadata data */
static
SCIP_RETCODE sepadataClear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   assert(sepadata != NULL);

   /* nlrowssize gets allocated first and then its decided whether to create the nlpiprob */
   if( sepadata->nlrowssize > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &sepadata->constraintviolation, sepadata->nlrowssize);
      SCIPfreeBlockMemoryArray(scip, &sepadata->convexsides, sepadata->nlrowssize);
      SCIPfreeBlockMemoryArray(scip, &sepadata->nlrows, sepadata->nlrowssize);
      sepadata->nlrowssize = 0;
   }

   if( sepadata->nlpiprob != NULL )
   {
      assert(sepadata->nlpi != NULL);

      SCIPfreeBlockMemoryArray(scip, &sepadata->nlpivars, sepadata->nlpinvars);

      SCIPhashmapFree(&sepadata->var2nlpiidx);
      SCIP_CALL( SCIPnlpiFreeProblem(sepadata->nlpi, &sepadata->nlpiprob) );
      SCIP_CALL( SCIPexprintFree(&sepadata->exprinterpreter) );

      sepadata->nlpinvars = 0;
      sepadata->nnlrows = 0;
   }
   assert(sepadata->nlpinvars == 0);
   assert(sepadata->nnlrows == 0);
   assert(sepadata->nlrowssize == 0);

   sepadata->skipsepa = FALSE;

   return SCIP_OKAY;
}

/** computes gradient of exprtree at projection */
static
SCIP_RETCODE computeGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expressions interpreter */
   SCIP_SOL*             projection,         /**< point where we compute gradient */
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
   assert(projection != NULL);
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
      x[i] = SCIPgetSolVal(scip, projection, SCIPexprtreeGetVars(exprtree)[i]);
   }

   SCIP_CALL( SCIPexprintGrad(exprint, exprtree, x, TRUE, &val, grad) );

   /*SCIPdebug( for( i = 0; i < nvars; ++i ) printf("%e [%s]\n", grad[i], SCIPvarGetName(SCIPexprtreeGetVars(exprtree)[i])) );*/

   SCIPfreeBufferArray(scip, &x);

   return SCIP_OKAY;
}

/** computes gradient cut (linearization) of nlrow at projection */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_SOL*             projection,         /**< point where we compute gradient cut */
   SCIP_NLROW*           nlrow,              /**< constraint for which we generate gradient cut */
   CONVEXSIDE            convexside,         /**< which side makes the nlrow convex */
   SCIP_Real             activity,           /**< activity of constraint at projection */
   SCIP_ROW**            row                 /**< storage for cut */
   )
{
   char rowname[SCIP_MAXSTRLEN];
   SCIP_SEPADATA* sepadata;
   SCIP_Real gradx0; /* <grad f(x_0), x_0> */
   int i;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(exprint != NULL);
   assert(nlrow != NULL);
   assert(row != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);

   gradx0 = 0.0;

   /* an nlrow has a linear part, quadratic part and expression tree; ideally one would just build the gradient but we
    * do not know if the different parts share variables or not, so we can't just build the gradient; for this reason
    * we create the row right away and compute the gradients of each part independently and add them to the row; the
    * row takes care to add coeffs corresponding to the same variable when they appear in different parts of the nlrow
    * NOTE: a gradient cut is globally valid whenever the constraint from which it is deduced is globally valid; since
    *       we build the convex relaxation using only globally valid constraints, the cuts are globally valid
    */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "proj_cut_%s_%u", SCIPnlrowGetName(nlrow), ++(sepadata->ncuts));
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, row, sepa, rowname, -SCIPinfinity(scip), SCIPinfinity(scip), TRUE, FALSE ,
            TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, *row) );

   /* linear part */
   for( i = 0; i < SCIPnlrowGetNLinearVars(nlrow); i++ )
   {
      gradx0 += SCIPgetSolVal(scip, projection, SCIPnlrowGetLinearVars(nlrow)[i]) * SCIPnlrowGetLinearCoefs(nlrow)[i];
      SCIP_CALL( SCIPaddVarToRow(scip, *row, SCIPnlrowGetLinearVars(nlrow)[i], SCIPnlrowGetLinearCoefs(nlrow)[i]) );
   }

   /* quadratic part */
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); i++ )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real grad1;
      SCIP_Real grad2;

      var1  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx1];
      var2  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx2];
      grad1 = SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, projection, var2);
      grad2 = SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, projection, var1);

      SCIP_CALL( SCIPaddVarToRow(scip, *row, var1, grad1) );
      SCIP_CALL( SCIPaddVarToRow(scip, *row, var2, grad2) );

      gradx0 += grad1 * SCIPgetSolVal(scip, projection, var1) + grad2 * SCIPgetSolVal(scip, projection, var2);
   }

   /* expression tree part */
   {
      SCIP_Real* grad;
      SCIP_EXPRTREE* tree;

      tree = SCIPnlrowGetExprtree(nlrow);

      if( tree != NULL && SCIPexprtreeGetNVars(tree) > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &grad, SCIPexprtreeGetNVars(tree)) );

         SCIP_CALL( computeGradient(scip, sepadata->exprinterpreter, projection, tree, grad) );

         for( i = 0; i < SCIPexprtreeGetNVars(tree); i++ )
         {
            gradx0 +=  grad[i] * SCIPgetSolVal(scip, projection, SCIPexprtreeGetVars(tree)[i]);
            SCIP_CALL( SCIPaddVarToRow(scip, *row, SCIPexprtreeGetVars(tree)[i], grad[i]) );
         }

         SCIPfreeBufferArray(scip, &grad);
      }
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, *row) );

   SCIPdebugPrintf("gradient: ");
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *row, NULL) ) );
   SCIPdebugPrintf("gradient dot x_0: %g\n", gradx0);

   /* gradient cut is f(x_0) - <grad f(x_0), x_0> + <grad f(x_0), x> <= rhs or >= lhs */
   if( convexside == RHS )
   {
      assert(!SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)));
      SCIP_CALL( SCIPchgRowRhs(scip, *row, SCIPnlrowGetRhs(nlrow) - activity + gradx0) );
   }
   else
   {
      assert(convexside == LHS);
      assert(!SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));
      SCIP_CALL( SCIPchgRowLhs(scip, *row, SCIPnlrowGetLhs(nlrow) - activity + gradx0) );
   }

   SCIPdebugPrintf("gradient cut: ");
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *row, NULL) ) );

   return SCIP_OKAY;
}

/** set quadratic part of objective function: \f$ \sum_i x_i^2 \f$; the objective function is \f$ ||x - x_0||^2 \f$,
 * where \f$ x_0 \f$ is the point to separate; the only part that changes is the term \f$ -2 \langle x_0, x \rangle \f$
 * which is linear and is set every time we want to separate a point, see separateCuts
 */
static
SCIP_RETCODE setQuadraticObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< the cut separator data */
   )
{
   SCIP_QUADELEM* quadelems;
   int i;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->nlpi != NULL);
   assert(sepadata->nlpiprob != NULL);
   assert(sepadata->var2nlpiidx != NULL);
   assert(sepadata->nlpinvars > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, sepadata->nlpinvars) );
   for( i = 0; i < sepadata->nlpinvars; i++ )
   {
      SCIP_VAR* var;

      var = sepadata->nlpivars[i];
      assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

      quadelems[i].idx1 = (int)(size_t)SCIPhashmapGetImage(sepadata->var2nlpiidx, (void*)var);
      quadelems[i].idx2 = quadelems[i].idx1;
      quadelems[i].coef = 1.0;
   }

   /* set quadratic part of objective function */
   SCIP_CALL( SCIPnlpiSetObjective(sepadata->nlpi, sepadata->nlpiprob,
            0, NULL, NULL, sepadata->nlpinvars, quadelems, NULL, NULL, 0.0) );

   /* free memory */
   SCIPfreeBufferArray(scip, &quadelems);

   return SCIP_OKAY;
}

/** projects sol onto convex relaxation (stored in sepadata) and tries to generate gradient cuts at the projection
 * it generates cuts only for the constraints that were violated by the LP solution and are now active or still
 * violated (in case we don't solve to optimality).
 * @todo: store a feasible solution if one is found to use as warmstart
 */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SOL*             sol,                /**< solution that should be separated */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SOL*      projection;
   SCIP_Real*     linvals;
   SCIP_Real*     nlpisol;
   SCIP_Real      timelimit;
   int            nlpinvars;
   int            i;
   int            iterlimit;
   int*           lininds;
   SCIP_Bool      nlpunstable;

   nlpunstable = FALSE;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(result != NULL);
   assert(sepadata != NULL);
   assert(sepadata->nnlrows > 0);
   assert(sepadata->nlpi != NULL);
   assert(sepadata->nlpinvars > 0);
   assert(sepadata->nlrows != NULL);
   assert(sepadata->nlpiprob != NULL);
   assert(sepadata->var2nlpiidx != NULL);
   assert(sepadata->convexsides != NULL);
   assert(sepadata->constraintviolation != NULL);

   nlpinvars = sepadata->nlpinvars;
   /* set linear part of objective function: \norm(x - x^0)^2 = \norm(x)^2 - \sum 2 * x_i * x^0_i + const
    * we ignore the constant; x0 is `sol`
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nlpinvars) );
   for( i = 0; i < nlpinvars; i++ )
   {
      SCIP_VAR* var;

      var = sepadata->nlpivars[i];
      assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

      lininds[i] = (int)(size_t)SCIPhashmapGetImage(sepadata->var2nlpiidx, (void*)var);
      linvals[i] = - 2.0 * SCIPgetSolVal(scip, sol, var);

      /* if coefficient is too large, don't separate */
      if( SCIPisHugeValue(scip, REALABS(linvals[i])) )
      {
         SCIPdebugMsg(scip, "Don't separate points too close to infinity\n");
         goto CLEANUP;
      }
   }

   /* set linear part of objective function */
   SCIP_CALL( SCIPnlpiChgLinearCoefs(sepadata->nlpi, sepadata->nlpiprob, -1, nlpinvars, lininds, linvals) );

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
         goto CLEANUP;
      }
   }
   if( sepadata->nlptimelimit > 0.0 )
      timelimit = MIN(sepadata->nlptimelimit, timelimit);
   SCIP_CALL( SCIPnlpiSetRealPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_TILIM, timelimit) );

   iterlimit = sepadata->nlpiterlimit > 0 ? sepadata->nlpiterlimit : INT_MAX;
   SCIP_CALL( SCIPnlpiSetIntPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_ITLIM, iterlimit) );
   SCIP_CALL( SCIPnlpiSetRealPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_FEASTOL, SCIPfeastol(scip) / 10.0) ); /* use tighter tolerances for the NLP solver */
   SCIP_CALL( SCIPnlpiSetRealPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_RELOBJTOL, MAX(SCIPfeastol(scip), SCIPdualfeastol(scip))) );  /*lint !e666*/
   SCIP_CALL( SCIPnlpiSetIntPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_VERBLEVEL, NLPVERBOSITY) );

   /* compute the projection onto the convex NLP relaxation */
   SCIP_CALL( SCIPnlpiSolve(sepadata->nlpi, sepadata->nlpiprob) );
   SCIPdebugMsg(scip, "NLP solstat = %d\n", SCIPnlpiGetSolstat(sepadata->nlpi, sepadata->nlpiprob));

   /* if solution is feasible, add cuts */
   switch( SCIPnlpiGetSolstat(sepadata->nlpi, sepadata->nlpiprob) )
   {
      case SCIP_NLPSOLSTAT_GLOBOPT:
      case SCIP_NLPSOLSTAT_LOCOPT:
         /* @todo: if solution is optimal, we might as well add the cut <x - P(x_0), x_0 - P(x_0)> <= 0
          * even though this cut is implied by all the gradient cuts of the rows active at the projection,
          * we do not add them all (only the gradient cuts of constraints that violated the LP solution */
      case SCIP_NLPSOLSTAT_FEASIBLE:
         /* generate cuts for violated constraints (at sol) that are active or still violated at the projection, since
          * a suboptimal solution or numerical issues could give a solution of the projection problem where constraints
          * are not active; if the solution of the projection problem is in the interior of the region, we do nothing
          */

         /* get solution: build SCIP_SOL out of nlpi sol */
         SCIP_CALL( SCIPnlpiGetSolution(sepadata->nlpi, sepadata->nlpiprob, &nlpisol, NULL, NULL, NULL, NULL) );
         assert(nlpisol != NULL);

         SCIP_CALL( SCIPcreateSol(scip, &projection, NULL) );
         for( i = 0; i < nlpinvars; i++ )
         {
            SCIP_VAR* var;

            var = sepadata->nlpivars[i];
            assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

            SCIP_CALL( SCIPsetSolVal(scip, projection, var,
                     nlpisol[(int)(size_t)SCIPhashmapGetImage(sepadata->var2nlpiidx, (void *)var)]) );
         }
         SCIPdebug( SCIPprintSol(scip, projection, NULL, TRUE) );

         /* check for active or violated constraints */
         for( i = 0; i < sepadata->nnlrows; ++i )
         {
            SCIP_NLROW* nlrow;
            CONVEXSIDE convexside;
            SCIP_Real activity;

            /* ignore constraints that are not violated by `sol` */
            if( SCIPisFeasZero(scip, sepadata->constraintviolation[i]) )
               continue;

            convexside = sepadata->convexsides[i];
            nlrow = sepadata->nlrows[i];
            assert(nlrow != NULL);

            /* check for currently active constraints at projected point */
            SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, projection, &activity) );

            SCIPdebugMsg(scip, "NlRow activity at nlpi solution: %g <= %g <= %g\n", SCIPnlrowGetLhs(nlrow), activity,
                  SCIPnlrowGetRhs(nlrow) );

            /* if nlrow is active or violates the projection, build gradient cut at projection */
            if( (convexside == RHS && SCIPisFeasGE(scip, activity, SCIPnlrowGetRhs(nlrow)))
               || (convexside == LHS && SCIPisFeasLE(scip, activity, SCIPnlrowGetLhs(nlrow))) )
            {
               SCIP_ROW* row;

               SCIP_CALL( generateCut(scip, sepa, sepadata->exprinterpreter, projection, nlrow, convexside, activity,
                        &row) );

               SCIPdebugMsg(scip, "active or violated nlrow: (sols vio: %e)\n", sepadata->constraintviolation[i]);
               SCIPdebug( SCIP_CALL( SCIPprintNlRow(scip, nlrow, NULL) ) );
               SCIPdebugMsg(scip, "cut with efficacy %g generated\n", SCIPgetCutEfficacy(scip, sol, row));
               SCIPdebug( SCIPprintRow(scip, row, NULL) );

               /* add cut if it is efficacious for the point we want to separate (sol) */
               if( SCIPisCutEfficacious(scip, sol, row) )
               {
                  SCIP_Bool infeasible;

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
         }

#ifdef SCIP_DEBUG
         {
            SCIP_Real distance;

            /* compute distance between LP sol and its projection (only makes sense when it is optimal) */
            distance = 0.0;
            for( i = 0; i < SCIPgetNNLPVars(scip); ++i )
            {
               SCIP_VAR* var;

               var = SCIPgetNLPVars(scip)[i];
               assert(var != NULL);

               /* assert NLP solution is within the bounds of the variable (only make sense when sol is optimal) */
               if( !SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
                  assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(var), SCIPvarGetNLPSol(var)));
               if( !SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
                  assert(SCIPisFeasLE(scip, SCIPvarGetNLPSol(var), SCIPvarGetUbLocal(var)));

               /*SCIPdebugMsg(scip, "NLP sol (LP sol): %s = %f (%g)\n", SCIPvarGetName(var),
                *     SCIPvarGetNLPSol(var), SCIPgetSolVal(scip, sol, var));
                */

               distance += SQR( SCIPvarGetNLPSol(var) - SCIPgetSolVal(scip, sol, var) );
            }

            SCIPdebugMsg(scip, "NLP objval: %e, distance: %e\n", SCIPgetNLPObjval(scip), distance);
         }
#endif

         /* free solution */
         SCIP_CALL( SCIPfreeSol(scip, &projection) );
         break;

      case SCIP_NLPSOLSTAT_GLOBINFEASIBLE:
      case SCIP_NLPSOLSTAT_LOCINFEASIBLE:
         /* fallthrough;
          * @todo: write what it means to be locinfeasible and why it can't be used to cutoff the node */
      case SCIP_NLPSOLSTAT_UNKNOWN:
         /* unknown... assume numerical issues */
         nlpunstable = TRUE;
         break;

      case SCIP_NLPSOLSTAT_UNBOUNDED:
      default:
         SCIPerrorMessage("Projection NLP is not unbounded by construction, should not get here!\n");
         SCIPABORT();
         nlpunstable = TRUE;
   }


   /* if nlp is detected to be unstable, don't try to separate again */
   if( nlpunstable )
   {
      /* @todo: maybe change objective function to \sum [(x_i - x_i^*)/max(|x_i^*|, 1)]^2
       * or some other scaling when unstable and try again.
       *       maybe free it here */
      sepadata->skipsepa = TRUE;
   }

   /* reset objective */
   BMSclearMemoryArray(linvals, nlpinvars);
   SCIP_CALL( SCIPnlpiChgLinearCoefs(sepadata->nlpi, sepadata->nlpiprob, -1, nlpinvars, lininds, linvals) );

CLEANUP:
   /* free memory */
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);

   return SCIP_OKAY;
}

/** computes the violation and maximum violation of the convex nlrows stored in sepadata wrt sol */
static
SCIP_RETCODE computeMaxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< solution that should be separated */
   SCIP_Real*            maxviolation        /**< buffer to store maximum violation */
   )
{
   SCIP_NLROW*    nlrow;
   int            i;

   assert(sepadata != NULL);
   assert(sepadata->nnlrows > 0);
   assert(sepadata->nlrows != NULL);
   assert(sepadata->convexsides != NULL);
   assert(sepadata->constraintviolation != NULL);

   *maxviolation = 0.0;
   for( i = 0; i < sepadata->nnlrows; i++ )
   {
      SCIP_Real activity;
      SCIP_Real violation;

      nlrow = sepadata->nlrows[i];

      /* get activity of nlrow */
      SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, sol, &activity) );

      /* violation = max{activity - rhs, 0.0} when convex and max{lhs - activity, 0.0} when concave */
      if( sepadata->convexsides[i] == RHS )
      {
         assert(SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONVEX);
         assert(!SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)));

         violation = activity - SCIPnlrowGetRhs(nlrow);
         sepadata->constraintviolation[i] = MAX(violation, 0.0);
      }
      if( sepadata->convexsides[i] == LHS )
      {
         assert(SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONCAVE);
         assert(!SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));

         violation = SCIPnlrowGetLhs(nlrow) - activity;
         sepadata->constraintviolation[i] = MAX(violation, 0.0);
      }

      /* compute maximum */
      if( *maxviolation < sepadata->constraintviolation[i] )
         *maxviolation = sepadata->constraintviolation[i];
   }

   SCIPdebugMsg(scip, "Maximum violation %g\n", *maxviolation);

   return SCIP_OKAY;
}


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

   SCIPdebugMsg(scip, "storing convex nlrows\n");

   sepadata->nlrowssize = nnlrows;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->nlrows), nnlrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->convexsides), nnlrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->constraintviolation), nnlrows) );

   /* count the number of nonlinear convex rows and store them */
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


/*
 * Callback methods of separator
 */


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeConvexproj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( sepadataClear(scip, sepadata) );

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolConvexproj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);

   SCIP_CALL( sepadataClear(scip, sepadata) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpConvexproj)
{  /*lint --e{715}*/

   SCIP_Real maxviolation;
   SCIP_SOL* lpsol;
   SCIP_SEPADATA* sepadata;

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* do not run if there is no interesting convex relaxation (with at least one nonlinear convex constraint),
    * or if we have found it to be numerically unstable
    * @todo: should it be with at least 2 nonlinear convex constraints?
    */
   if( sepadata->skipsepa )
   {
      SCIPdebugMsg(scip, "not running because convex relaxation is uninteresting or numerically unstable\n");
      return SCIP_OKAY;
   }

   /* the separator needs an NLP solver */
   if( SCIPgetNNlpis(scip) == 0 )
      return SCIP_OKAY;

   /* only call separator up to a maximum depth */
   if( sepadata->maxdepth >= 0 && SCIPgetDepth(scip) > sepadata->maxdepth )
      return SCIP_OKAY;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* do not run if SCIP does not have constructed an NLP */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "NLP not constructed, skipping convex projection separator\n");
      return SCIP_OKAY;
   }

   /* recompute convex NLP relaxation if the variable set changed and we are still at the root node
    * @todo: does it make sense to do this??? */
   if( sepadata->nlpiprob != NULL && SCIPgetNVars(scip) != sepadata->nlpinvars  && SCIPgetDepth(scip) == 0 )
   {
      SCIP_CALL( sepadataClear(scip, sepadata) );
      assert(sepadata->nlpiprob == NULL);
   }

   /* create or update convex NLP relaxation */
   if( sepadata->nlpiprob == NULL )
   {
      /* store convex nonlinear constraints */
      SCIP_CALL( storeNonlinearConvexNlrows(scip, sepadata, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip)) );

      /* check that convex NLP relaxation is interesting (more than one nonlinear constraint) */
      if( sepadata->nnlrows < 1 )
      {
         SCIPdebugMsg(scip, "convex relaxation uninteresting, don't run\n");
         sepadata->skipsepa = TRUE;
         return SCIP_OKAY;
      }

      /* create the expression interpreter */
      SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &sepadata->exprinterpreter) );

      sepadata->nlpinvars = SCIPgetNVars(scip);
      sepadata->nlpi = SCIPgetNlpis(scip)[0];
      assert(sepadata->nlpi != NULL);

      SCIP_CALL( SCIPnlpiCreateProblem(sepadata->nlpi, &sepadata->nlpiprob, "convexproj-nlp") );
      SCIP_CALL( SCIPhashmapCreate(&sepadata->var2nlpiidx, SCIPblkmem(scip), sepadata->nlpinvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &sepadata->nlpivars, SCIPgetVars(scip), sepadata->nlpinvars) ); /*lint !e666*/

      SCIP_CALL( SCIPcreateNlpiProb(scip, sepadata->nlpi, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip),
            sepadata->nlpiprob, sepadata->var2nlpiidx, NULL, SCIPgetCutoffbound(scip), FALSE, TRUE) );

      /* add rows of the LP */
      if( SCIPgetDepth(scip) == 0 )
      {
         SCIP_CALL( SCIPaddNlpiProbRows(scip, sepadata->nlpi, sepadata->nlpiprob, sepadata->var2nlpiidx,
                  SCIPgetLPRows(scip), SCIPgetNLPRows(scip)) );
      }

      /* set quadratic part of objective function */
      SCIP_CALL( setQuadraticObj(scip, sepadata) );
   }
   else
   {
      SCIP_CALL( SCIPupdateNlpiProb(scip, sepadata->nlpi, sepadata->nlpiprob, sepadata->var2nlpiidx,
            sepadata->nlpivars, sepadata->nlpinvars, SCIPgetCutoffbound(scip)) );
   }

   /* assert that the lp solution satisfies the cutoff bound; if this fails then we shouldn't have a cutoff bound in the
    * nlpi, since then the projection could be in the interior of the actual convex relaxation */
   assert(SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL ||
         SCIPisLE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)));

   /* get current sol: LP or pseudo solution if LP sol is not available */
   SCIP_CALL( SCIPcreateCurrentSol(scip, &lpsol, NULL) );

   /* do not run if current solution's violation is small */
   SCIP_CALL( computeMaxViolation(scip, sepadata, lpsol, &maxviolation) );
   if( maxviolation < VIOLATIONFAC * SCIPfeastol(scip) )
   {
      SCIPdebugMsg(scip, "solution doesn't violate constraints enough, do not separate\n");
      SCIP_CALL( SCIPfreeSol(scip, &lpsol) );
      return SCIP_OKAY;
   }

   /* run the separator */
   *result = SCIP_DIDNOTFIND;

   /* separateCuts computes the projection and then gradient cuts on each constraint that was originally violated */
   SCIP_CALL( separateCuts(scip, sepa, lpsol, result) );

   /* free memory */
   SCIP_CALL( SCIPfreeSol(scip, &lpsol) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the convexproj separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaConvexproj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create convexproj separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );

   /* this sets all data in sepadata to 0 */
   BMSclearMemory(sepadata);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpConvexproj, NULL,
         sepadata) );
   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeConvexproj) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolConvexproj) );

   /* add convexproj separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxdepth",
         "maximal depth at which the separator is applied (-1: unlimited)",
         &sepadata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nlpiterlimit",
         "iteration limit of NLP solver; 0 for no limit",
         &sepadata->nlpiterlimit, TRUE, DEFAULT_NLPITERLIM, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/nlptimelimit",
         "time limit of NLP solver; 0.0 for no limit",
         &sepadata->nlptimelimit, TRUE, DEFAULT_NLPTIMELIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
