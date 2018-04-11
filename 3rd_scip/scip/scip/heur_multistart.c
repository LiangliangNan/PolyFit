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

/**@file   heur_multistart.c
 * @brief  multistart heuristic for convex and nonconvex MINLPs
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_multistart.h"
#include "scip/heur_subnlp.h"

#include "nlpi/exprinterpret.h"

#define HEUR_NAME             "multistart"
#define HEUR_DESC             "multistart heuristic for convex and nonconvex MINLPs"
#define HEUR_DISPCHAR         'm'
#define HEUR_PRIORITY         -2100000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE           /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_RANDSEED      131            /**< initial random seed */
#define DEFAULT_NRNDPOINTS    100            /**< default number of generated random points per call */
#define DEFAULT_MAXBOUNDSIZE  2e+4           /**< default maximum variable domain size for unbounded variables */
#define DEFAULT_MAXITER       300            /**< default number of iterations to reduce the violation of a point */
#define DEFAULT_MINIMPRFAC    0.05           /**< default minimum required improving factor to proceed in improvement of a point */
#define DEFAULT_MINIMPRITER   10             /**< default number of iteration when checking the minimum improvement */
#define DEFAULT_MAXRELDIST    0.15           /**< default maximum distance between two points in the same cluster */
#define DEFAULT_NLPMINIMPR    0.00           /**< default factor by which heuristic should at least improve the incumbent */
#define DEFAULT_GRADLIMIT     5e+6           /**< default limit for gradient computations for all improvePoint() calls */
#define DEFAULT_MAXNCLUSTER   3              /**< default maximum number of considered clusters per heuristic call */
#define DEFAULT_ONLYNLPS      TRUE           /**< should the heuristic run only on continuous problems? */

#define MINFEAS               -1e+4          /**< minimum feasibility for a point; used for filtering and improving
                                              *   feasibility */
#define MINIMPRFAC            0.95           /**< improvement factor used to discard randomly generated points with a
                                              *   too large objective value */
#define GRADCOSTFAC_LINEAR    1.0            /**< gradient cost factor for the number of linear variables */
#define GRADCOSTFAC_QUAD      2.0            /**< gradient cost factor for the number of quadratic terms */
#define GRADCOSTFAC_NONLINEAR 3.0            /**< gradient cost factor for the number of nodes in nonlinear expression */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_EXPRINT*         exprinterpreter;    /**< expression interpreter to compute gradients */
   int                   nrndpoints;         /**< number of random points generated per execution call */
   SCIP_Real             maxboundsize;       /**< maximum variable domain size for unbounded variables */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_HEUR*            heursubnlp;         /**< sub-NLP heuristic */

   int                   maxiter;            /**< number of iterations to reduce the maximum violation of a point */
   SCIP_Real             minimprfac;         /**< minimum required improving factor to proceed in the improvement of a single point */
   int                   minimpriter;        /**< number of iteration when checking the minimum improvement */

   SCIP_Real             maxreldist;         /**< maximum distance between two points in the same cluster */
   SCIP_Real             nlpminimpr;         /**< factor by which heuristic should at least improve the incumbent */
   SCIP_Real             gradlimit;          /**< limit for gradient computations for all improvePoint() calls (0 for no limit) */
   int                   maxncluster;        /**< maximum number of considered clusters per heuristic call */
   SCIP_Bool             onlynlps;           /**< should the heuristic run only on continuous problems? */
};


/*
 * Local methods
 */


/** returns an unique index of a variable in the range of 0,..,SCIPgetNVars(scip)-1 */
#ifndef NDEBUG
static
int getVarIndex(
   SCIP_HASHMAP*         varindex,           /**< maps variables to indicies between 0,..,SCIPgetNVars(scip)-1 */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert(varindex != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(varindex, (void*)var));

   return (int)(size_t)SCIPhashmapGetImage(varindex, (void*)var);
}
#else
#define getVarIndex(varindex,var) ((int)(size_t)SCIPhashmapGetImage((varindex), (void*)(var)))
#endif

/** samples and stores random points; stores points which have a better objective value than the current incumbent
 *  solution
 */
static
SCIP_RETCODE sampleRandomPoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            rndpoints,          /**< array to store all random points */
   int                   nmaxrndpoints,      /**< maximum number of random points to compute */
   SCIP_Real             maxboundsize,       /**< maximum variable domain size for unbounded variables */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             bestobj,            /**< objective value in the transformed space of the current incumbent */
   int*                  nstored             /**< pointer to store the number of randomly generated points */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* sol;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   int nvars;
   int niter;
   int i;

   assert(scip != NULL);
   assert(rndpoints != NULL);
   assert(nmaxrndpoints > 0);
   assert(maxboundsize > 0.0);
   assert(randnumgen != NULL);
   assert(nstored != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   *nstored = 0;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   for( niter = 0; niter < 3 * nmaxrndpoints && *nstored < nmaxrndpoints; ++niter )
   {
      for( i = 0; i < nvars; ++i )
      {
         lb = MIN(SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])); /*lint !e666*/
         ub = MAX(SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])); /*lint !e666*/

         if( SCIPisFeasEQ(scip, lb, ub) )
            val = (lb + ub) / 2.0;
         /* use a smaller domain for unbounded variables */
         else if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
            val = SCIPrandomGetReal(randnumgen, lb, ub);
         else if( !SCIPisInfinity(scip, -lb) )
            val = lb + SCIPrandomGetReal(randnumgen, 0.0, maxboundsize);
         else if( !SCIPisInfinity(scip, ub) )
            val = ub - SCIPrandomGetReal(randnumgen, 0.0, maxboundsize);
         else
         {
            assert(SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub));
            val = SCIPrandomGetReal(randnumgen, -0.5*maxboundsize, 0.5*maxboundsize);
         }
         assert(SCIPisFeasGE(scip, val, lb) && SCIPisFeasLE(scip, val, ub));

         /* set solution value; round the sampled point for integer variables */
         if( SCIPvarGetType(vars[i]) < SCIP_VARTYPE_CONTINUOUS )
            val = SCIPfeasRound(scip, val);
         SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], val) );
      }

      /* add solution if it is good enough */
      if( SCIPisLE(scip, SCIPgetSolTransObj(scip, sol), bestobj) )
      {
         SCIP_CALL( SCIPcreateSolCopy(scip, &rndpoints[*nstored], sol) );
         ++(*nstored);
      }
   }
   assert(*nstored <= nmaxrndpoints);
   SCIPdebugMsg(scip, "found %d randomly generated points\n", *nstored);

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}

/** computes the minimum feasibility of a given point; a negative value means that there is an infeasibility */
static
SCIP_RETCODE getMinFeas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrows,             /**< array containing all nlrows */
   int                   nnlrows,            /**< total number of nlrows */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Real*            minfeas             /**< buffer to store the minimum feasibility */
   )
{
   SCIP_Real tmp;
   int i;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(minfeas != NULL);
   assert(nlrows != NULL);
   assert(nnlrows > 0);

   *minfeas = SCIPinfinity(scip);

   for( i = 0; i < nnlrows; ++i )
   {
      assert(nlrows[i] != NULL);

      SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrows[i], sol, &tmp) );
      *minfeas = MIN(*minfeas, tmp);
   }

   return SCIP_OKAY;
}

/** computes the gradient for a given point and nonlinear row */
static
SCIP_RETCODE computeGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_EXPRINT*         exprint,            /**< expressions interpreter */
   SCIP_SOL*             sol,                /**< solution to compute the gradient for */
   SCIP_HASHMAP*         varindex,           /**< maps variables to indicies between 0,..,SCIPgetNVars(scip)-1 uniquely */
   SCIP_Real*            grad,               /**< buffer to store the gradient; grad[varindex(i)] corresponds to SCIPgetVars(scip)[i] */
   SCIP_Real*            norm                /**< buffer to store ||grad||^2  */
   )
{
   SCIP_EXPRTREE* tree;
   SCIP_VAR* var;
   int i;

   assert(scip != NULL);
   assert(nlrow != NULL);
   assert(varindex != NULL);
   assert(exprint != NULL);
   assert(sol != NULL);
   assert(norm != NULL);

   BMSclearMemoryArray(grad, SCIPgetNVars(scip));
   *norm = 0.0;

   /* linear part */
   for( i = 0; i < SCIPnlrowGetNLinearVars(nlrow); i++ )
   {
      var = SCIPnlrowGetLinearVars(nlrow)[i];
      assert(var != NULL);
      assert(getVarIndex(varindex, var) >= 0 && getVarIndex(varindex, var) < SCIPgetNVars(scip));

      grad[getVarIndex(varindex, var)] += SCIPnlrowGetLinearCoefs(nlrow)[i];
   }

   /* quadratic part */
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); i++ )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;

      var1  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx1];
      var2  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx2];

      assert(SCIPnlrowGetQuadElems(nlrow)[i].idx1 < SCIPnlrowGetNQuadVars(nlrow));
      assert(SCIPnlrowGetQuadElems(nlrow)[i].idx2 < SCIPnlrowGetNQuadVars(nlrow));
      assert(getVarIndex(varindex, var1) >= 0 && getVarIndex(varindex, var1) < SCIPgetNVars(scip));
      assert(getVarIndex(varindex, var2) >= 0 && getVarIndex(varindex, var2) < SCIPgetNVars(scip));

      grad[getVarIndex(varindex, var1)] += SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, sol, var2);
      grad[getVarIndex(varindex, var2)] += SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, sol, var1);
   }

   /* tree part */
   tree = SCIPnlrowGetExprtree(nlrow);
   if( tree != NULL )
   {
      SCIP_Real* treegrad;
      SCIP_Real* x;
      SCIP_Real val;

      assert(SCIPexprtreeGetNVars(tree) <= SCIPgetNVars(scip));

      SCIP_CALL( SCIPallocBufferArray(scip, &x, SCIPexprtreeGetNVars(tree)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &treegrad, SCIPexprtreeGetNVars(tree)) );

      /* compile expression tree, if not done before */
      if( SCIPexprtreeGetInterpreterData(tree) == NULL )
      {
         SCIP_CALL( SCIPexprintCompile(exprint, tree) );
      }

      /* sets the solution value */
      for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
         x[i] = SCIPgetSolVal(scip, sol, SCIPexprtreeGetVars(tree)[i]);

      SCIP_CALL( SCIPexprintGrad(exprint, tree, x, TRUE, &val, treegrad) );

      /* update corresponding gradient entry */
      for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
      {
         var = SCIPexprtreeGetVars(tree)[i];
         assert(var != NULL);
         assert(getVarIndex(varindex, var) >= 0 && getVarIndex(varindex, var) < SCIPgetNVars(scip));

         grad[getVarIndex(varindex, var)] += treegrad[i];
      }

      SCIPfreeBufferArray(scip, &treegrad);
      SCIPfreeBufferArray(scip, &x);
   }

   /* compute ||grad||^2 */
   for( i = 0; i < SCIPgetNVars(scip); ++i )
      *norm += SQR(grad[i]);

   return SCIP_OKAY;
}

/** use consensus vectors to improve feasibility for a given starting point */
static
SCIP_RETCODE improvePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrows,             /**< array containing all nlrows */
   int                   nnlrows,            /**< total number of nlrows */
   SCIP_HASHMAP*         varindex,           /**< maps variables to indicies between 0,..,SCIPgetNVars(scip)-1 */
   SCIP_EXPRINT*         exprinterpreter,    /**< expression interpreter */
   SCIP_SOL*             point,              /**< random generated point */
   int                   maxiter,            /**< maximum number of iterations */
   SCIP_Real             minimprfac,         /**< minimum required improving factor to proceed in the improvement of a single point */
   int                   minimpriter,        /**< number of iteration when checking the minimum improvement */
   SCIP_Real*            minfeas,            /**< pointer to store the minimum feasibility */
   SCIP_Real*            nlrowgradcosts,     /**< estimated costs for each gradient computation */
   SCIP_Real*            gradcosts           /**< pointer to store the estimated gradient costs */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* grad;
   SCIP_Real* updatevec;
   SCIP_Real lastminfeas;
   int nvars;
   int r;
   int i;

   assert(varindex != NULL);
   assert(exprinterpreter != NULL);
   assert(point != NULL);
   assert(maxiter > 0);
   assert(minfeas != NULL);
   assert(nlrows != NULL);
   assert(nnlrows > 0);
   assert(nlrowgradcosts != NULL);
   assert(gradcosts != NULL);

   *gradcosts = 0.0;

   SCIP_CALL( getMinFeas(scip, nlrows, nnlrows, point, minfeas) );
#ifdef SCIP_DEBUG_IMPROVEPOINT
   printf("start minfeas = %e\n", *minfeas);
#endif

   /* stop since start point is feasible */
   if( !SCIPisFeasLT(scip, *minfeas, 0.0) )
   {
#ifdef SCIP_DEBUG_IMPROVEPOINT
      printf("start point is feasible");
#endif
      return SCIP_OKAY;
   }

   lastminfeas = *minfeas;
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &grad, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &updatevec, nvars) );

   /* main loop */
   for( r = 0; r < maxiter && SCIPisFeasLT(scip, *minfeas, 0.0); ++r )
   {
      SCIP_Real feasibility;
      SCIP_Real activity;
      SCIP_Real nlrownorm;
      SCIP_Real scale;
      int nviolnlrows;

      BMSclearMemoryArray(updatevec, nvars);
      nviolnlrows = 0;

      for( i = 0; i < nnlrows; ++i )
      {
         int j;

         SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrows[i], point, &feasibility) );

         /* do not consider non-violated constraints */
         if( SCIPisFeasGE(scip, feasibility, 0.0) )
            continue;

         /* increase number of violated nlrows */
         ++nviolnlrows;

         SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrows[i], point, &activity) );
         SCIP_CALL( computeGradient(scip, nlrows[i], exprinterpreter, point, varindex, grad, &nlrownorm) );

         /* update estimated costs for computing gradients */
         *gradcosts += nlrowgradcosts[i];

         /* stop if the gradient disappears at the current point */
         if( SCIPisZero(scip, nlrownorm) )
         {
#ifdef SCIP_DEBUG_IMPROVEPOINT
            printf("gradient vanished at current point -> stop\n");
#endif
            goto TERMINATE;
         }

         /* compute -g(x_k) / ||grad(g)(x_k)||^2 for a constraint g(x_k) <= 0 */
         scale = -feasibility / nlrownorm;
         if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrows[i])) && SCIPisGT(scip, activity, SCIPnlrowGetRhs(nlrows[i])) )
            scale *= -1.0;

         /* skip nonliner row if the scaler is too small or too large */
         if( SCIPisEQ(scip, scale, 0.0) || SCIPisHugeValue(scip, REALABS(scale)) )
            continue;

         for( j = 0; j < nvars; ++j )
            updatevec[j] += scale * grad[j];
      }
      assert(nviolnlrows > 0);

      for( i = 0; i < nvars; ++i )
      {
         /* adjust point */
         updatevec[i] = SCIPgetSolVal(scip, point, vars[i]) + updatevec[i] / nviolnlrows;
         updatevec[i] = MIN(updatevec[i], SCIPvarGetUbLocal(vars[i])); /*lint !e666*/
         updatevec[i] = MAX(updatevec[i], SCIPvarGetLbLocal(vars[i])); /*lint !e666*/

         SCIP_CALL( SCIPsetSolVal(scip, point, vars[i], updatevec[i]) );
      }

      /* update feasibility */
      SCIP_CALL( getMinFeas(scip, nlrows, nnlrows, point, minfeas) );

      /* check stopping criterion */
      if( r % minimpriter == 0 && r > 0 )
      {
         if( *minfeas <= MINFEAS
            || (*minfeas-lastminfeas) / MAX(REALABS(*minfeas), REALABS(lastminfeas)) < minimprfac ) /*lint !e666*/
            break;
         lastminfeas = *minfeas;
      }
   }

TERMINATE:
#ifdef SCIP_DEBUG_IMPROVEPOINT
   printf("niter=%d minfeas=%e\n", r, *minfeas);
#endif

   SCIPfreeBufferArray(scip, &grad);
   SCIPfreeBufferArray(scip, &updatevec);

   return SCIP_OKAY;
}

/** sorts points w.r.t their feasibilities; points with a feasibility which is too small (w.r.t. the geometric mean of
 *  all feasibilities) will be filtered out
 */
static
SCIP_RETCODE filterPoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            points,             /**< array containing improved points */
   SCIP_Real*            feasibilities,      /**< array containing feasibility for each point (sorted) */
   int                   npoints,            /**< total number of points */
   int*                  nusefulpoints       /**< pointer to store the total number of useful points */
   )
{
   SCIP_Real minfeas;
   SCIP_Real meanfeas;
   int i;

   assert(points != NULL);
   assert(feasibilities != NULL);
   assert(npoints > 0);
   assert(nusefulpoints != NULL);

   /* sort points w.r.t their feasibilities; non-negative feasibility correspond to feasible points for the NLP */
   SCIPsortDownRealPtr(feasibilities, (void**)points, npoints);
   minfeas = feasibilities[npoints - 1];

   /* check if all points are feasible */
   if( SCIPisFeasGE(scip, minfeas, 0.0) )
   {
      *nusefulpoints = npoints;
      return SCIP_OKAY;
   }

   *nusefulpoints = 0;

   /* compute shifted geometric mean of feasibilities (shift value = 1 - minfeas) */
   meanfeas = 1.0;
   for( i = 0; i < npoints; ++i )
   {
      assert(feasibilities[i] - minfeas + 1.0 > 0.0);
      meanfeas *= pow(feasibilities[i] - minfeas + 1.0, 1.0 / npoints);
   }
   meanfeas += minfeas - 1.0;
   SCIPdebugMsg(scip, "meanfeas = %e\n", meanfeas);

   /* keep all points with which have a feasibility not much below the geometric mean of infeasibilities */
   for( i = 0; i < npoints; ++i )
   {
      if( SCIPisFeasLT(scip, feasibilities[i], 0.0)
         && (feasibilities[i] <= 1.05 * meanfeas || SCIPisLE(scip, feasibilities[i], MINFEAS)) )
         break;

      ++(*nusefulpoints);
   }

   return SCIP_OKAY;
}

/** returns the relative distance between two points; considers a smaller bounded domain for unbounded variables */
static
SCIP_Real getRelDistance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             x,                  /**< first point */
   SCIP_SOL*             y,                  /**< second point */
   SCIP_Real             maxboundsize        /**< maximum variable domain size for unbounded variables */
   )
{
   SCIP_VAR** vars;
   SCIP_Real distance;
   SCIP_Real solx;
   SCIP_Real soly;
   SCIP_Real lb;
   SCIP_Real ub;
   int i;

   assert(x != NULL);
   assert(y != NULL);
   assert(SCIPgetNVars(scip) > 0);

   vars = SCIPgetVars(scip);
   distance = 0.0;

   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      lb = SCIPvarGetLbLocal(vars[i]);
      ub = SCIPvarGetUbLocal(vars[i]);
      solx = SCIPgetSolVal(scip, x, vars[i]);
      soly = SCIPgetSolVal(scip, y, vars[i]);

      /* adjust lower and upper bounds for unbounded variables*/
      if( SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
      {
         lb = -maxboundsize / 2.0;
         ub = +maxboundsize / 2.0;
      }
      else if( SCIPisInfinity(scip, -lb) )
      {
         lb = ub - maxboundsize;
      }
      else if( SCIPisInfinity(scip, ub) )
      {
         ub = lb + maxboundsize;
      }

      /* project solution values to the variable domain */
      solx = MIN(MAX(solx, lb), ub);
      soly = MIN(MAX(soly, lb), ub);

      distance += REALABS(solx - soly) / MAX(1.0, ub - lb);
   }

   return distance / SCIPgetNVars(scip);
}

/** cluster useful points with a greedy algorithm */
static
SCIP_RETCODE clusterPointsGreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            points,             /**< array containing improved points */
   int                   npoints,            /**< total number of points */
   int*                  clusteridx,         /**< array to store for each point the index of the cluster */
   int*                  ncluster,           /**< pointer to store the total number of cluster */
   SCIP_Real             maxboundsize,       /**< maximum variable domain size for unbounded variables */
   SCIP_Real             maxreldist,         /**< maximum relative distance between any two points of the same cluster */
   int                   maxncluster         /**< maximum number of clusters to compute */
   )
{
   int i;

   assert(points != NULL);
   assert(npoints > 0);
   assert(clusteridx != NULL);
   assert(ncluster != NULL);
   assert(maxreldist >= 0.0);
   assert(maxncluster >= 0);

   /* initialize cluster indices */
   for( i = 0; i < npoints; ++i )
      clusteridx[i] = INT_MAX;

   *ncluster = 0;

   for( i = 0; i < npoints && (*ncluster < maxncluster); ++i )
   {
      int j;

      /* point is already assigned to a cluster */
      if( clusteridx[i] != INT_MAX )
         continue;

      /* create a new cluster for i */
      clusteridx[i] = *ncluster;

      for( j = i + 1; j < npoints; ++j )
      {
         if( clusteridx[j] == INT_MAX && getRelDistance(scip, points[i], points[j], maxboundsize) <= maxreldist )
            clusteridx[j] = *ncluster;
      }

      ++(*ncluster);
   }

#ifndef NDEBUG
   for( i = 0; i < npoints; ++i )
   {
      assert(clusteridx[i] >= 0);
      assert(clusteridx[i] < *ncluster || clusteridx[i] == INT_MAX);
   }
#endif

   return SCIP_OKAY;
}

/** calls the sub-NLP heuristic for a given cluster */
static
SCIP_RETCODE solveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< multi-start heuristic */
   SCIP_HEUR*            nlpheur,            /**< pointer to NLP local search heuristics */
   SCIP_SOL**            points,             /**< array containing improved points */
   int                   npoints,            /**< total number of points */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver */
   SCIP_Real             timelimit,          /**< time limit for NLP solver */
   SCIP_Real             minimprove,         /**< desired minimal relative improvement in objective function value */
   SCIP_Bool*            success             /**< pointer to store if we could find a solution */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* refpoint;
   SCIP_RESULT nlpresult;
   SCIP_Real val;
   int nbinvars;
   int nintvars;
   int nvars;
   int i;

   assert(points != NULL);
   assert(npoints > 0);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   *success = FALSE;

   SCIP_CALL( SCIPcreateSol(scip, &refpoint, heur) );

   /* compute reference point */
   for( i = 0; i < nvars; ++i )
   {
      int p;

      val = 0.0;

      for( p = 0; p < npoints; ++p )
      {
         assert(points[p] != NULL);
         val += SCIPgetSolVal(scip, points[p], vars[i]);
      }

      SCIP_CALL( SCIPsetSolVal(scip, refpoint, vars[i], val / npoints) );
   }

   /* round point for sub-NLP heuristic */
   SCIP_CALL( SCIProundSol(scip, refpoint, success) );
   SCIPdebugMsg(scip, "rounding of refpoint successfully? %u\n", *success);

   /* round variables manually if the locks did not allow us to round them */
   if( !(*success) )
   {
      for( i = 0; i < nbinvars + nintvars; ++i )
      {
         val = SCIPgetSolVal(scip, refpoint, vars[i]);

         if( !SCIPisFeasIntegral(scip, val) )
         {
            assert(SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(vars[i])));
            assert(SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(vars[i])));

            /* round and adjust value */
            val = SCIPround(scip, val);
            val = MIN(val, SCIPvarGetUbLocal(vars[i])); /*lint !e666*/
            val = MAX(val, SCIPvarGetLbLocal(vars[i])); /*lint !e666*/
            assert(SCIPisFeasIntegral(scip, val));

            SCIP_CALL( SCIPsetSolVal(scip, refpoint, vars[i], val) );
         }
      }
   }

   /* call sub-NLP heuristic */
   SCIP_CALL( SCIPapplyHeurSubNlp(scip, nlpheur, &nlpresult, refpoint, itercontingent, timelimit, minimprove,
         NULL, NULL) );
   SCIP_CALL( SCIPfreeSol(scip, &refpoint) );

   /* let sub-NLP heuristic decide whether the solution is feasible or not */
   *success = nlpresult == SCIP_FOUNDSOL;

   return SCIP_OKAY;
}

/** recursive helper function to count the number of nodes in a sub-tree */
static
int getExprSize(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int sum;
   int i;

   assert(expr != NULL);

   sum = 0;
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
      sum += getExprSize(child);
   }
   return 1 + sum;
}

/** returns the number of nodes in an expression tree */
static
int getExprtreeSize(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   if( tree == NULL )
      return 0;
   return getExprSize(SCIPexprtreeGetRoot(tree));
}

/** main function of the multi-start heuristic (see @ref heur_multistart.h for more details); it consists of the
 *  following four steps:
 *
 *  1. sampling points in the current domain; for unbounded variables we use a bounded box
 *
 *  2. reduce infeasibility by using a gradient descent method
 *
 *  3. cluster points; filter points with a too large infeasibility
 *
 *  4. compute start point for each cluster and use it in the sub-NLP heuristic (@ref heur_subnlp.h)
 */
static
SCIP_RETCODE applyHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_NLROW** nlrows;
   SCIP_SOL** points;
   SCIP_HASHMAP* varindex;
   SCIP_Real* feasibilities;
   SCIP_Real* nlrowgradcosts;
   int* clusteridx;
   SCIP_Real gradlimit;
   SCIP_Real bestobj;
   int nusefulpoints;
   int nrndpoints;
   int ncluster;
   int nnlrows;
   int npoints;
   int start;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);
   assert(heurdata != NULL);

   SCIPdebugMsg(scip, "call applyHeur()\n");

   nlrows = SCIPgetNLPNlRows(scip);
   nnlrows = SCIPgetNNLPNlRows(scip);
   bestobj = SCIPgetNSols(scip) > 0 ? MINIMPRFAC * SCIPgetSolTransObj(scip, SCIPgetBestSol(scip)) : SCIPinfinity(scip);

   if( heurdata->exprinterpreter == NULL )
   {
      SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &heurdata->exprinterpreter) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &points, heurdata->nrndpoints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlrowgradcosts, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &feasibilities, heurdata->nrndpoints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &clusteridx, heurdata->nrndpoints) );
   SCIP_CALL( SCIPhashmapCreate(&varindex, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   /* create an unique mapping of all variables to 0,..,SCIPgetNVars(scip)-1 */
   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(varindex, (void*)SCIPgetVars(scip)[i], (void*)(size_t)i) );
   }

   /* compute estimated costs of computing a gradient for each nlrow */
   for( i = 0; i < nnlrows; ++i )
   {
      nlrowgradcosts[i] = GRADCOSTFAC_LINEAR * SCIPnlrowGetNLinearVars(nlrows[i])
         + GRADCOSTFAC_QUAD * SCIPnlrowGetNQuadElems(nlrows[i])
         + GRADCOSTFAC_NONLINEAR * getExprtreeSize(SCIPnlrowGetExprtree(nlrows[i]));
   }

   /*
    * 1. sampling points in the current domain; for unbounded variables we use a bounded box
    */
   SCIP_CALL( sampleRandomPoints(scip, points, heurdata->nrndpoints, heurdata->maxboundsize, heurdata->randnumgen,
         bestobj, &nrndpoints) );
   assert(nrndpoints >= 0);

   if( nrndpoints == 0 )
      goto TERMINATE;

   /*
    * 2. improve points via consensus vectors
    */
   gradlimit = heurdata->gradlimit == 0.0 ? SCIPinfinity(scip) : heurdata->gradlimit;
   for( npoints = 0; npoints < nrndpoints && gradlimit >= 0 && !SCIPisStopped(scip); ++npoints )
   {
      SCIP_Real gradcosts;

      SCIP_CALL( improvePoint(scip, nlrows, nnlrows, varindex, heurdata->exprinterpreter, points[npoints],
            heurdata->maxiter, heurdata->minimprfac, heurdata->minimpriter, &feasibilities[npoints], nlrowgradcosts,
            &gradcosts) );

      gradlimit -= gradcosts;
      SCIPdebugMsg(scip, "improve point %d / %d gradlimit = %g\n", npoints, nrndpoints, gradlimit);
   }
   assert(npoints >= 0 && npoints <= nrndpoints);

   if( npoints == 0 )
      goto TERMINATE;

   /*
    * 3. filter and cluster points
    */
   SCIP_CALL( filterPoints(scip, points, feasibilities, npoints, &nusefulpoints) );
   assert(nusefulpoints >= 0);
   SCIPdebugMsg(scip, "nusefulpoints = %d\n", nusefulpoints);

   if( nusefulpoints == 0 )
      goto TERMINATE;

   SCIP_CALL( clusterPointsGreedy(scip, points, nusefulpoints, clusteridx, &ncluster, heurdata->maxboundsize,
         heurdata->maxreldist, heurdata->maxncluster) );
   assert(ncluster >= 0 && ncluster <= heurdata->maxncluster);
   SCIPdebugMsg(scip, "ncluster = %d\n", ncluster);

   SCIPsortIntPtr(clusteridx, (void**)points, nusefulpoints);

   /*
    * 4. compute start point for each cluster and use it in the sub-NLP heuristic (@ref heur_subnlp.h)
    */
   start = 0;
   while( start < nusefulpoints && clusteridx[start] != INT_MAX && !SCIPisStopped(scip) )
   {
      SCIP_Real timelimit;
      SCIP_Bool success;
      int end;

      end = start;
      while( end < nusefulpoints && clusteridx[start] == clusteridx[end] )
         ++end;

      assert(end - start > 0);

      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      if( !SCIPisInfinity(scip, timelimit) )
         timelimit -= SCIPgetSolvingTime(scip);

      /* try to solve sub-NLP if we have enough time left */
      if( timelimit <= 1.0 )
      {
         SCIPdebugMsg(scip, "not enough time left! (%g)\n", timelimit);
         break;
      }

      /* call sub-NLP heuristic */
      SCIP_CALL( solveNLP(scip, heur, heurdata->heursubnlp, &points[start], end - start, -1LL, timelimit,
            heurdata->nlpminimpr, &success) );
      SCIPdebugMsg(scip, "solveNLP result = %d\n", success);

      if( success )
         *result = SCIP_FOUNDSOL;

      /* go to the next cluster */
      start = end;
   }

TERMINATE:
   /* free memory */
   for( i = nrndpoints - 1; i >= 0 ; --i )
   {
      assert(points[i] != NULL);
      SCIP_CALL( SCIPfreeSol(scip, &points[i]) );
   }

   SCIPhashmapFree(&varindex);
   SCIPfreeBufferArray(scip, &clusteridx);
   SCIPfreeBufferArray(scip, &feasibilities);
   SCIPfreeBufferArray(scip, &nlrowgradcosts);
   SCIPfreeBufferArray(scip, &points);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyMultistart)
{  /*lint --e{715}*/
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurMultistart(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMultistart)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);

   if( heurdata->exprinterpreter != NULL )
   {
      SCIP_CALL( SCIPexprintFree(&heurdata->exprinterpreter) );
   }

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitMultistart)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED) );

   /* try to find sub-NLP heuristic */
   heurdata->heursubnlp = SCIPfindHeur(scip, "subnlp");

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitMultistart)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->randnumgen != NULL);

   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMultistart)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check cases for which the heuristic is not applicable */
   if( !SCIPisNLPConstructed(scip) || heurdata->heursubnlp == NULL || SCIPgetNNlpis(scip) <= 0 )
      return SCIP_OKAY;

   /* check whether the heuristic should be applied for a problem containing integer variables */
   if( heurdata->onlynlps && (SCIPgetNBinVars(scip) > 0 || SCIPgetNIntVars(scip) > 0) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( applyHeur(scip, heur, heurdata, result) );

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the multistart primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMultistart(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create multistart primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMultistart, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyMultistart) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMultistart) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitMultistart) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitMultistart) );

   /* add multistart primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nrndpoints",
         "number of random points generated per execution call",
         &heurdata->nrndpoints, FALSE, DEFAULT_NRNDPOINTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxboundsize",
         "maximum variable domain size for unbounded variables",
         &heurdata->maxboundsize, FALSE, DEFAULT_MAXBOUNDSIZE, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxiter",
         "number of iterations to reduce the maximum violation of a point",
         &heurdata->maxiter, FALSE, DEFAULT_MAXITER, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprfac",
         "minimum required improving factor to proceed in improvement of a single point",
         &heurdata->minimprfac, FALSE, DEFAULT_MINIMPRFAC, -SCIPinfinity(scip), SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minimpriter",
         "number of iteration when checking the minimum improvement",
         &heurdata->minimpriter, FALSE, DEFAULT_MINIMPRITER, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxreldist",
         "maximum distance between two points in the same cluster",
         &heurdata->maxreldist, FALSE, DEFAULT_MAXRELDIST, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nlpminimpr",
         "factor by which heuristic should at least improve the incumbent",
         &heurdata->nlpminimpr, FALSE, DEFAULT_NLPMINIMPR, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/gradlimit",
         "limit for gradient computations for all improvePoint() calls (0 for no limit)",
         &heurdata->gradlimit, FALSE, DEFAULT_GRADLIMIT, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxncluster",
         "maximum number of considered clusters per heuristic call",
         &heurdata->maxncluster, FALSE, DEFAULT_MAXNCLUSTER, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/onlynlps",
         "should the heuristic run only on continuous problems?",
         &heurdata->onlynlps, FALSE, DEFAULT_ONLYNLPS, NULL, NULL) );

   return SCIP_OKAY;
}
