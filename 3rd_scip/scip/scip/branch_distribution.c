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

/**@file   branch_distribution.c
 * @ingroup BRANCHINGRULES
 * @brief  probability based branching rule based on an article by J. Pryor and J.W. Chinneck
 * @author Gregor Hendel
 *
 * The distribution branching rule selects a variable based on its impact on row activity distributions. More formally,
 * let \f$ a(x) = a_1 x_1 + \dots + a_n x_n \leq b \f$ be a row of the LP. Let further \f$ l_i, u_i \in R\f$ denote the
 * (finite) lower and upper bound, respectively, of the \f$ i \f$-th variable \f$x_i\f$.
 * Viewing every variable \f$x_i \f$ as (continuously) uniformly distributed within its bounds, we can approximately
 * understand the row activity \f$a(x)\f$ as a Gaussian random variate with mean value \f$ \mu = E[a(x)] = \sum_i a_i\frac{l_i + u_i}{2}\f$
 * and variance \f$ \sigma^2 = \sum_i a_i^2 \sigma_i^2 \f$, with \f$ \sigma_i^2 = \frac{(u_i - l_i)^2}{12}\f$ for
 * continuous and \f$ \sigma_i^2 = \frac{(u_i - l_i + 1)^2 - 1}{12}\f$ for discrete variables.
 * With these two parameters, we can calculate the probability to satisfy the row in terms of the cumulative distribution
 * of the standard normal distribution: \f$ P(a(x) \leq b) = \Phi(\frac{b - \mu}{\sigma})\f$.
 *
 * The impact of a variable on the probability to satisfy a constraint after branching can be estimated by altering
 * the variable contribution to the sums described above. In order to keep the tree size small,
 * variables are preferred which decrease the probability
 * to satisfy a row because it is more likely to create infeasible subproblems.
 *
 * The selection of the branching variable is subject to the parameter @p scoreparam. For both branching directions,
 * an individual score is calculated. Available options for scoring methods are:
 * - @b d: select a branching variable with largest difference in satisfaction probability after branching
 * - @b l: lowest cumulative probability amongst all variables on all rows (after branching).
 * - @b h: highest cumulative probability amongst all variables on all rows (after branching).
 * - @b v: highest number of votes for lowering row probability for all rows a variable appears in.
 * - @b w: highest number of votes for increasing row probability for all rows a variable appears in.
 *
 * If the parameter @p usescipscore is set to @a TRUE, a single branching score is calculated from the respective
 * up and down scores as defined by the SCIP branching score function (see advanced branching parameter @p scorefunc),
 * otherwise, the variable with the single highest score is selected, and the maximizing direction is assigned
 * higher branching priority.
 *
 * The original idea of probability based branching schemes appeared in:
 *
 * J. Pryor and J.W. Chinneck:@n
 * Faster Integer-Feasibility in Mixed-Integer Linear Programs by Branching to Force Change@n
 * Computers and Operations Research, vol. 38, 2011, p. 1143–1152@n
 * (Paper: <a href="http://www.sce.carleton.ca/faculty/chinneck/docs/PryorChinneck.pdf">PDF</a>).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/branch_distribution.h"


#define BRANCHRULE_NAME            "distribution"
#define BRANCHRULE_DESC            "branching rule based on variable influence on cumulative normal distribution of row activities"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define SCOREPARAM_VALUES          "dhlvw"
#define DEFAULT_SCOREPARAM         'v'
#define DEFAULT_PRIORITY           2.0
#define SQRTOFTWO                  1.4142136
#define SQUARED(x) ((x) * (x))
#define DEFAULT_ONLYACTIVEROWS     FALSE     /**< should only rows which are active at the current node be considered? */
#define DEFAULT_USEWEIGHTEDSCORE   FALSE     /**< should the branching score weigh up- and down-scores of a variable */

/* branching rule event handler to catch variable bound changes */
#define EVENTHDLR_NAME             "eventhdlr_distribution"    /**< event handler name */
#define EVENT_DISTRIBUTION         SCIP_EVENTTYPE_BOUNDCHANGED /**< the event tyoe to be handled by this event handler */

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler pointer */
   SCIP_VAR**            updatedvars;        /**< variables to process bound change events for */
   SCIP_Real*            rowmeans;           /**< row activity mean values for all rows */
   SCIP_Real*            rowvariances;       /**< row activity variances for all rows */
   SCIP_Real*            currentubs;         /**< variable upper bounds as currently saved in the row activities */
   SCIP_Real*            currentlbs;         /**< variable lower bounds as currently saved in the row activities */
   int*                  rowinfinitiesdown;  /**< count the number of variables with infinite bounds which allow for always
                                              *   repairing the constraint right hand side */
   int*                  rowinfinitiesup;    /**< count the number of variables with infinite bounds which allow for always
                                              *   repairing the constraint left hand side */
   int*                  varposs;            /**< array of variable positions in the updated variables array */
   int*                  varfilterposs;      /**< array of event filter positions for variable events */
   int                   nupdatedvars;       /**< the current number of variables with pending bound changes */
   int                   memsize;            /**< memory size of current arrays, needed for dynamic reallocation */
   int                   varpossmemsize;     /**< memory size of updated vars and varposs array */
   char                  scoreparam;         /**< parameter how the branch score is calculated */
   SCIP_Bool             onlyactiverows;     /**< should only rows which are active at the current node be considered? */
   SCIP_Bool             usescipscore;       /**< should the branching use SCIP's branching score function */
};

struct SCIP_EventhdlrData
{
   SCIP_BRANCHRULEDATA*  branchruledata;     /**< the branching rule data to access distribution arrays */
};

/*
 * Local methods
 */

/** ensure that maxindex + 1 rows can be represented in data arrays; memory gets reallocated with 10% extra space
 *  to save some time for future allocations */
static
SCIP_RETCODE branchruledataEnsureArraySize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branchruledata */
   int                   maxindex            /**< row index at hand (size must be at least this large) */
   )
{
   int newsize;
   int r;

   /* maxindex fits in current array -> nothing to do */
   if( maxindex < branchruledata->memsize )
      return SCIP_OKAY;

   /* new memory size is the max index + 1 plus 10% additional space */
   newsize = (int)SCIPfeasCeil(scip, (maxindex + 1) * 1.1);
   assert(newsize > branchruledata->memsize);
   assert(branchruledata->memsize >= 0);

   /* alloc memory arrays for row information */
   if( branchruledata->memsize == 0 )
   {
      SCIP_VAR** vars;
      int v;
      int nvars;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->rowinfinitiesdown, newsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->rowinfinitiesup, newsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->rowmeans, newsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->rowvariances, newsize) );

      assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      assert(nvars > 0);

      /* allocate variable update event processing array storage */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->varfilterposs, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->varposs, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->updatedvars, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->currentubs, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->currentlbs, nvars) );

      branchruledata->varpossmemsize = nvars;
      branchruledata->nupdatedvars = 0;

      /* init variable event processing data */
      for( v = 0; v < nvars; ++v )
      {
         assert(SCIPvarIsActive(vars[v]));
         assert(SCIPvarGetProbindex(vars[v]) == v);

         /* set up variable events to catch bound changes */
         SCIP_CALL( SCIPcatchVarEvent(scip, vars[v], EVENT_DISTRIBUTION, branchruledata->eventhdlr, NULL, &(branchruledata->varfilterposs[v])) );
         assert(branchruledata->varfilterposs[v] >= 0);

         branchruledata->varposs[v] = -1;
         branchruledata->updatedvars[v] = NULL;
         branchruledata->currentlbs[v] = SCIP_INVALID;
         branchruledata->currentubs[v] = SCIP_INVALID;
      }

   }
   else
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->rowinfinitiesdown, branchruledata->memsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->rowinfinitiesup, branchruledata->memsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->rowmeans, branchruledata->memsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &branchruledata->rowvariances, branchruledata->memsize, newsize) );
   }

   /* loop over extended arrays and invalidate data to trigger initialization of this row when necessary */
   for( r = branchruledata->memsize; r < newsize; ++r )
   {
      branchruledata->rowmeans[r] = SCIP_INVALID;
      branchruledata->rowvariances[r] = SCIP_INVALID;
      branchruledata->rowinfinitiesdown[r] = 0;
      branchruledata->rowinfinitiesup[r] = 0;
   }

   /* adjust memsize */
   branchruledata->memsize = newsize;

   return SCIP_OKAY;
}

/* update the variables current lower and upper bound */
static
void branchruledataUpdateCurrentBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branchrule data */
   SCIP_VAR*             var                 /**< the variable to update current bounds */
   )
{
   int varindex;
   SCIP_Real lblocal;
   SCIP_Real ublocal;

   assert(var != NULL);

   varindex = SCIPvarGetProbindex(var);
   assert(0 <= varindex && varindex < branchruledata->varpossmemsize);
   lblocal = SCIPvarGetLbLocal(var);
   ublocal = SCIPvarGetUbLocal(var);

   assert(SCIPisFeasLE(scip, lblocal, ublocal));

   branchruledata->currentlbs[varindex] = lblocal;
   branchruledata->currentubs[varindex] = ublocal;
}

/** calculate the variable's distribution parameters (mean and variance) for the bounds specified in the arguments.
 *  special treatment of infinite bounds necessary */
void SCIPvarCalcDistributionParameters(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             varlb,              /**< variable lower bound */
   SCIP_Real             varub,              /**< variable upper bound */
   SCIP_VARTYPE          vartype,            /**< type of the variable */
   SCIP_Real*            mean,               /**< pointer to store mean value */
   SCIP_Real*            variance            /**< pointer to store the variance of the variable uniform distribution */
   )
{
   assert(mean != NULL);
   assert(variance != NULL);

   /* calculate mean and variance of variable uniform distribution before and after branching */
   if( SCIPisInfinity(scip, varub) || SCIPisInfinity(scip, -varlb) )
   {
      /* variables with infinite bounds are not kept in the row activity variance */
      *variance = 0.0;

      if( !SCIPisInfinity(scip, varub) )
         *mean = varub;
      else if( !SCIPisInfinity(scip, -varlb) )
         *mean = varlb;
      else
         *mean = 0.0;
   }
   else
   {
      /* if the variable is continuous, we assume a continuous uniform distribution, otherwise a discrete one */
      if( vartype == SCIP_VARTYPE_CONTINUOUS )
         *variance = SQUARED(varub - varlb);
      else
         *variance = SQUARED(varub - varlb + 1) - 1;
      *variance /= 12.0;
      *mean = (varub + varlb) * .5;
   }

   assert(!SCIPisNegative(scip, *variance));
}

/** calculates the cumulative distribution P(-infinity <= x <= value) that a normally distributed
 *  random variable x takes a value between -infinity and parameter \p value.
 *
 *  The distribution is given by the respective mean and deviation. This implementation
 *  uses the error function SCIPerf().
 */
SCIP_Real SCIPcalcCumulativeDistribution(
   SCIP*                 scip,               /**< current SCIP */
   SCIP_Real             mean,               /**< the mean value of the distribution */
   SCIP_Real             variance,           /**< the square of the deviation of the distribution */
   SCIP_Real             value               /**< the upper limit of the calculated distribution integral */
   )
{
   SCIP_Real normvalue;
   SCIP_Real std;

   /* we need to calculate the standard deviation from the variance */
   assert(!SCIPisNegative(scip, variance));
   if( SCIPisFeasZero(scip, variance) )
      std = 0.0;
   else
      std = sqrt(variance);

   /* special treatment for zero variance */
   if( SCIPisFeasZero(scip, std) )
   {
      if( SCIPisFeasLE(scip, value, mean) )
         return 1.0;
      else
         return 0.0;
   }
   assert( std != 0.0 ); /* for lint */

   /* scale and translate to standard normal distribution. Factor sqrt(2) is needed for SCIPerf() function */
   normvalue = (value - mean)/(std * SQRTOFTWO);

   SCIPdebugMsg(scip, " Normalized value %g = ( %g - %g ) / (%g * 1.4142136)\n", normvalue, value, mean, std);

   /* calculate the cumulative distribution function for normvalue. For negative normvalues, we negate
    * the normvalue and use the oddness of the SCIPerf()-function; special treatment for values close to zero.
    */
   if( SCIPisFeasZero(scip, normvalue) )
      return .5;
   else if( normvalue > 0 )
   {
      SCIP_Real erfresult;

      erfresult = SCIPerf(normvalue);
      return  erfresult / 2.0 + 0.5;
   }
   else
   {
      SCIP_Real erfresult;

      erfresult = SCIPerf(-normvalue);

      return 0.5 - erfresult / 2.0;
   }
}

/** calculates the probability of satisfying an LP-row under the assumption
 *  of uniformly distributed variable values.
 *
 *  For inequalities, we use the cumulative distribution function of the standard normal
 *  distribution PHI(rhs - mu/sqrt(sigma2)) to calculate the probability
 *  for a right hand side row with mean activity mu and variance sigma2 to be satisfied.
 *  Similarly, 1 - PHI(lhs - mu/sqrt(sigma2)) is the probability to satisfy a left hand side row.
 *  For equations (lhs==rhs), we use the centeredness measure p = min(PHI(lhs'), 1-PHI(lhs'))/max(PHI(lhs'), 1 - PHI(lhs')),
 *  where lhs' = lhs - mu / sqrt(sigma2).
 */
SCIP_Real SCIProwCalcProbability(
   SCIP*                 scip,               /**< current scip */
   SCIP_ROW*             row,                /**< the row */
   SCIP_Real             mu,                 /**< the mean value of the row distribution */
   SCIP_Real             sigma2,             /**< the variance of the row distribution */
   int                   rowinfinitiesdown,  /**< the number of variables with infinite bounds to DECREASE activity */
   int                   rowinfinitiesup     /**< the number of variables with infinite bounds to INCREASE activity */
   )
{
   SCIP_Real rowprobability;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real lhsprob;
   SCIP_Real rhsprob;

   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);

   lhsprob = 1.0;
   rhsprob = 1.0;

   /* use the cumulative distribution if the row contains no variable to repair every infeasibility */
   if( !SCIPisInfinity(scip, rhs) && rowinfinitiesdown == 0 )
      rhsprob = SCIPcalcCumulativeDistribution(scip, mu, sigma2, rhs);

   /* use the cumulative distribution if the row contains no variable to repair every infeasibility
    * otherwise the row can always be made feasible by increasing activity far enough
    */
   if( !SCIPisInfinity(scip, -lhs) && rowinfinitiesup == 0 )
      lhsprob = 1.0 - SCIPcalcCumulativeDistribution(scip, mu, sigma2, lhs);

   assert(SCIPisFeasLE(scip, lhsprob, 1.0) && SCIPisFeasGE(scip, lhsprob, 0.0));
   assert(SCIPisFeasLE(scip, rhsprob, 1.0) && SCIPisFeasGE(scip, rhsprob, 0.0));

   /* use centeredness measure for equations; for inequalities, the minimum of the two probabilities is the
    * probability to satisfy the row */
   if( SCIPisFeasEQ(scip, lhs, rhs) )
   {
      SCIP_Real minprobability;
      SCIP_Real maxprobability;

      minprobability = MIN(rhsprob, lhsprob);
      maxprobability = MAX(lhsprob, rhsprob);
      rowprobability = minprobability / maxprobability;
   }
   else
      rowprobability = MIN(rhsprob, lhsprob);

   SCIPdebug( SCIPprintRow(scip, row, NULL) );
   SCIPdebugMsg(scip, " Row %s, mean %g, sigma2 %g, LHS %g, RHS %g has probability %g to be satisfied\n",
      SCIProwGetName(row), mu, sigma2, lhs, rhs, rowprobability);

   assert(SCIPisFeasGE(scip, rowprobability, 0.0) && SCIPisFeasLE(scip, rowprobability, 1.0));

   return rowprobability;
}

/** calculates the initial mean and variance of the row activity normal distribution.
 *
 *  The mean value \f$ \mu \f$ is given by \f$ \mu = \sum_i=1^n c_i * (lb_i +ub_i) / 2 \f$ where
 *  \f$n \f$ is the number of variables, and \f$ c_i, lb_i, ub_i \f$ are the variable coefficient and
 *  bounds, respectively. With the same notation, the variance \f$ \sigma^2 \f$ is given by
 *  \f$ \sigma^2 = \sum_i=1^n c_i^2 * \sigma^2_i \f$, with the variance being
 *  \f$ \sigma^2_i = ((ub_i - lb_i + 1)^2 - 1) / 12 \f$ for integer variables and
 *  \f$ \sigma^2_i = (ub_i - lb_i)^2 / 12 \f$ for continuous variables.
 */
static
void rowCalculateGauss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< the branching rule data */
   SCIP_ROW*             row,                /**< the row for which the gaussian normal distribution has to be calculated */
   SCIP_Real*            mu,                 /**< pointer to store the mean value of the gaussian normal distribution */
   SCIP_Real*            sigma2,             /**< pointer to store the variance value of the gaussian normal distribution */
   int*                  rowinfinitiesdown,  /**< pointer to store the number of variables with infinite bounds to DECREASE activity */
   int*                  rowinfinitiesup     /**< pointer to store the number of variables with infinite bounds to INCREASE activity */
   )
{
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int nrowvals;
   int c;

   assert(scip != NULL);
   assert(row != NULL);
   assert(mu != NULL);
   assert(sigma2 != NULL);
   assert(rowinfinitiesup != NULL);
   assert(rowinfinitiesdown != NULL);

   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   nrowvals = SCIProwGetNNonz(row);

   assert(nrowvals == 0 || rowcols != NULL);
   assert(nrowvals == 0 || rowvals != NULL);

   *mu = SCIProwGetConstant(row);
   *sigma2 = 0.0;
   *rowinfinitiesdown = 0;
   *rowinfinitiesup = 0;

   /* loop over nonzero row coefficients and sum up the variable contributions to mu and sigma2 */
   for( c = 0; c < nrowvals; ++c )
   {
      SCIP_VAR* colvar;
      SCIP_Real colval;
      SCIP_Real colvarlb;
      SCIP_Real colvarub;
      SCIP_Real squarecoeff;
      SCIP_Real varvariance;
      SCIP_Real varmean;
      int varindex;

      assert(rowcols[c] != NULL);
      colvar = SCIPcolGetVar(rowcols[c]);
      assert(colvar != NULL);

      colval = rowvals[c];
      colvarlb = SCIPvarGetLbLocal(colvar);
      colvarub = SCIPvarGetUbLocal(colvar);

      varmean = 0.0;
      varvariance = 0.0;
      varindex = SCIPvarGetProbindex(colvar);
      assert((branchruledata->currentlbs[varindex] == SCIP_INVALID) == (branchruledata->currentubs[varindex] == SCIP_INVALID)); /*lint !e777*/

      /* variable bounds need to be watched from now on */
      if( branchruledata->currentlbs[varindex] == SCIP_INVALID ) /*lint !e777*/
         branchruledataUpdateCurrentBounds(scip, branchruledata, colvar);

      assert(!SCIPisInfinity(scip, colvarlb));
      assert(!SCIPisInfinity(scip, -colvarub));
      assert(SCIPisFeasLE(scip, colvarlb, colvarub));

      /* variables with infinite bounds are skipped for the calculation of the variance; they need to
       * be accounted for by the counters for infinite row activity decrease and increase and they
       * are used to shift the row activity mean in case they have one nonzero, but finite bound */
      if( SCIPisInfinity(scip, -colvarlb) || SCIPisInfinity(scip, colvarub) )
      {
         if( SCIPisInfinity(scip, colvarub) )
         {
         /* an infinite upper bound gives the row an infinite maximum activity or minimum activity, if the coefficient is
          * positive or negative, resp.
          */
            if( colval < 0.0 )
               ++(*rowinfinitiesdown);
            else
               ++(*rowinfinitiesup);
         }

         /* an infinite lower bound gives the row an infinite maximum activity or minimum activity, if the coefficient is
          * negative or positive, resp.
          */
         if( SCIPisInfinity(scip, -colvarlb) )
         {
            if( colval > 0.0 )
               ++(*rowinfinitiesdown);
            else
               ++(*rowinfinitiesup);
         }
      }
      SCIPvarCalcDistributionParameters(scip, colvarlb, colvarub, SCIPvarGetType(colvar), &varmean, &varvariance);

      /* actual values are updated; the contribution of the variable to mu is the arithmetic mean of its bounds */
      *mu += colval * varmean;

      /* the variance contribution of a variable is c^2 * (u - l)^2 / 12.0 for continuous and c^2 * ((u - l + 1)^2 - 1) / 12.0 for integer */
      squarecoeff = SQUARED(colval);
      *sigma2 += squarecoeff * varvariance;

      assert(!SCIPisFeasNegative(scip, *sigma2));
   }

   SCIPdebug( SCIPprintRow(scip, row, NULL) );
   SCIPdebugMsg(scip, "  Row %s has a mean value of %g at a sigma2 of %g \n", SCIProwGetName(row), *mu, *sigma2);
}

/** update the up- and downscore of a single variable after calculating the impact of branching on a
 *  particular row, depending on the chosen score parameter
 */
SCIP_RETCODE SCIPupdateDistributionScore(
   SCIP*                 scip,               /**< current SCIP pointer */
   SCIP_Real             currentprob,        /**< the current probability */
   SCIP_Real             newprobup,          /**< the new probability if branched upwards */
   SCIP_Real             newprobdown,        /**< the new probability if branched downwards */
   SCIP_Real*            upscore,            /**< pointer to store the new score for branching up */
   SCIP_Real*            downscore,          /**< pointer to store the new score for branching down */
   char                  scoreparam          /**< parameter to determine the way the score is calculated */
   )
{
   assert(scip != NULL);
   assert(SCIPisFeasGE(scip, currentprob, 0.0) && SCIPisFeasLE(scip, currentprob, 1.0));
   assert(SCIPisFeasGE(scip, newprobup, 0.0) && SCIPisFeasLE(scip, newprobup, 1.0));
   assert(SCIPisFeasGE(scip, newprobdown, 0.0) && SCIPisFeasLE(scip, newprobdown, 1.0));
   assert(upscore != NULL);
   assert(downscore != NULL);

   /* update up and downscore depending on score parameter */
   switch( scoreparam )
   {
   case 'l' :
      /* 'l'owest cumulative probability */
      if( SCIPisGT(scip, 1.0 - newprobup, *upscore) )
         *upscore = 1.0 - newprobup;
      if( SCIPisGT(scip, 1.0 - newprobdown, *downscore) )
         *downscore = 1.0 - newprobdown;
      break;

   case 'd' :
      /* biggest 'd'ifference currentprob - newprob */
      if( SCIPisGT(scip, currentprob - newprobup, *upscore) )
         *upscore = currentprob - newprobup;
      if( SCIPisGT(scip, currentprob - newprobdown, *downscore) )
         *downscore = currentprob - newprobdown;
      break;

   case 'h' :
      /* 'h'ighest cumulative probability */
      if( SCIPisGT(scip, newprobup, *upscore) )
         *upscore = newprobup;
      if( SCIPisGT(scip, newprobdown, *downscore) )
         *downscore = newprobdown;
      break;

   case 'v' :
      /* 'v'otes lowest cumulative probability */
      if( SCIPisLT(scip, newprobup, newprobdown) )
         *upscore += 1.0;
      else if( SCIPisGT(scip, newprobup, newprobdown) )
         *downscore += 1.0;
      break;

   case 'w' :
      /* votes highest cumulative probability */
      if( SCIPisGT(scip, newprobup, newprobdown) )
         *upscore += 1.0;
      else if( SCIPisLT(scip, newprobup, newprobdown) )
         *downscore += 1.0;
      break;

   default :
      SCIPerrorMessage(" ERROR! No branching scheme selected! Exiting  method.\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

/** calculate the branching score of a variable, depending on the chosen score parameter */
static
SCIP_RETCODE calcBranchScore(
   SCIP*                 scip,               /**< current SCIP */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branch rule data */
   SCIP_VAR*             var,                /**< candidate variable */
   SCIP_Real             lpsolval,           /**< current fractional LP-relaxation solution value  */
   SCIP_Real*            upscore,            /**< pointer to store the variable score when branching on it in upward direction */
   SCIP_Real*            downscore,          /**< pointer to store the variable score when branching on it in downward direction */
   char                  scoreparam          /**< the score parameter of this branching rule */
   )
{
   SCIP_COL* varcol;
   SCIP_ROW** colrows;
   SCIP_Real* rowvals;
   SCIP_Real varlb;
   SCIP_Real varub;
   SCIP_Real squaredbounddiff; /* current squared difference of variable bounds (ub - lb)^2 */
   SCIP_Real newub;            /* new upper bound if branching downwards */
   SCIP_Real newlb;            /* new lower bound if branching upwards */
   SCIP_Real squaredbounddiffup; /* squared difference after branching upwards (ub - lb')^2 */
   SCIP_Real squaredbounddiffdown; /* squared difference after branching downwards (ub' - lb)^2 */
   SCIP_Real currentmean;      /* current mean value of variable uniform distribution */
   SCIP_Real meanup;           /* mean value of variable uniform distribution after branching up */
   SCIP_Real meandown;         /* mean value of variable uniform distribution after branching down*/
   SCIP_VARTYPE vartype;
   int ncolrows;
   int i;

   SCIP_Bool onlyactiverows; /* should only rows which are active at the current node be considered? */

   assert(scip != NULL);
   assert(var != NULL);
   assert(upscore != NULL);
   assert(downscore != NULL);
   assert(!SCIPisIntegral(scip, lpsolval));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   varcol = SCIPvarGetCol(var);
   assert(varcol != NULL);

   colrows = SCIPcolGetRows(varcol);
   rowvals = SCIPcolGetVals(varcol);
   ncolrows = SCIPcolGetNNonz(varcol);
   varlb = SCIPvarGetLbLocal(var);
   varub = SCIPvarGetUbLocal(var);
   assert(SCIPisFeasLT(scip, varlb, varub));
   vartype = SCIPvarGetType(var);

   /* calculate mean and variance of variable uniform distribution before and after branching */
   currentmean = 0.0;
   squaredbounddiff = 0.0;
   SCIPvarCalcDistributionParameters(scip, varlb, varub, vartype, &currentmean, &squaredbounddiff);

   newlb = SCIPfeasCeil(scip, lpsolval);
   newub = SCIPfeasFloor(scip, lpsolval);

   /* calculate the variable's uniform distribution after branching up and down, respectively. */
   squaredbounddiffup = 0.0;
   meanup = 0.0;
   SCIPvarCalcDistributionParameters(scip, newlb, varub, vartype, &meanup, &squaredbounddiffup);

   /* calculate the distribution mean and variance for a variable with finite lower bound */
   squaredbounddiffdown = 0.0;
   meandown = 0.0;
   SCIPvarCalcDistributionParameters(scip, varlb, newub, vartype, &meandown, &squaredbounddiffdown);

   /* initialize the variable's up and down score */
   *upscore = 0.0;
   *downscore = 0.0;

   onlyactiverows = branchruledata->onlyactiverows;

   /* loop over the variable rows and calculate the up and down score */
   for( i = 0; i < ncolrows; ++i )
   {
      SCIP_ROW* row;
      SCIP_Real changedrowmean;
      SCIP_Real rowmean;
      SCIP_Real rowvariance;
      SCIP_Real changedrowvariance;
      SCIP_Real currentrowprob;
      SCIP_Real newrowprobup;
      SCIP_Real newrowprobdown;
      SCIP_Real squaredcoeff;
      SCIP_Real rowval;
      int rowinfinitiesdown;
      int rowinfinitiesup;
      int rowpos;

      row = colrows[i];
      rowval = rowvals[i];
      assert(row != NULL);

      /* we access the rows by their index */
      rowpos = SCIProwGetIndex(row);

      /* skip non-active rows if the user parameter was set this way */
      if( onlyactiverows && SCIPisSumPositive(scip, SCIPgetRowLPFeasibility(scip, row)) )
         continue;

      /* call method to ensure sufficient data capacity */
      SCIP_CALL( branchruledataEnsureArraySize(scip, branchruledata, rowpos) );

      /* calculate row activity distribution if this is the first candidate to appear in this row */
      if( branchruledata->rowmeans[rowpos] == SCIP_INVALID ) /*lint !e777*/
      {
         rowCalculateGauss(scip, branchruledata, row, &branchruledata->rowmeans[rowpos], &branchruledata->rowvariances[rowpos],
               &branchruledata->rowinfinitiesdown[rowpos], &branchruledata->rowinfinitiesup[rowpos]);
      }

      /* retrieve the row distribution parameters from the branch rule data */
      rowmean = branchruledata->rowmeans[rowpos];
      rowvariance = branchruledata->rowvariances[rowpos];
      rowinfinitiesdown = branchruledata->rowinfinitiesdown[rowpos];
      rowinfinitiesup = branchruledata->rowinfinitiesdown[rowpos];
      assert(!SCIPisNegative(scip, rowvariance));

      currentrowprob = SCIProwCalcProbability(scip, row, rowmean, rowvariance,
            rowinfinitiesdown, rowinfinitiesup);

      /* get variable's current expected contribution to row activity */
      squaredcoeff = SQUARED(rowval);

      /* first, get the probability change for the row if the variable is branched on upwards. The probability
       * can only be affected if the variable upper bound is finite
       */
      if( !SCIPisInfinity(scip, varub) )
      {
         int rowinftiesdownafterbranch;
         int rowinftiesupafterbranch;

         /* calculate how branching would affect the row parameters */
         changedrowmean = rowmean + rowval * (meanup - currentmean);
         changedrowvariance = rowvariance + squaredcoeff * (squaredbounddiffup - squaredbounddiff);
         changedrowvariance = MAX(0.0, changedrowvariance);

         rowinftiesdownafterbranch = rowinfinitiesdown;
         rowinftiesupafterbranch = rowinfinitiesup;

         /* account for changes of the row's infinite bound contributions */
         if( SCIPisInfinity(scip, -varlb) && rowval < 0.0 )
            rowinftiesupafterbranch--;
         if( SCIPisInfinity(scip, -varlb) && rowval > 0.0 )
            rowinftiesdownafterbranch--;

         assert(rowinftiesupafterbranch >= 0);
         assert(rowinftiesdownafterbranch >= 0);
         newrowprobup = SCIProwCalcProbability(scip, row, changedrowmean, changedrowvariance, rowinftiesdownafterbranch,
               rowinftiesupafterbranch);
      }
      else
         newrowprobup = currentrowprob;

      /* do the same for the other branching direction */
      if( !SCIPisInfinity(scip, varlb) )
      {
         int rowinftiesdownafterbranch;
         int rowinftiesupafterbranch;

         changedrowmean = rowmean + rowval * (meandown - currentmean);
         changedrowvariance = rowvariance + squaredcoeff * (squaredbounddiffdown - squaredbounddiff);
         changedrowvariance = MAX(0.0, changedrowvariance);

         rowinftiesdownafterbranch = rowinfinitiesdown;
         rowinftiesupafterbranch = rowinfinitiesup;

         /* account for changes of the row's infinite bound contributions */
         if( SCIPisInfinity(scip, varub) && rowval > 0.0 )
            rowinftiesupafterbranch -= 1;
         if( SCIPisInfinity(scip, varub) && rowval < 0.0 )
            rowinftiesdownafterbranch -= 1;

         assert(rowinftiesdownafterbranch >= 0);
         assert(rowinftiesupafterbranch >= 0);
         newrowprobdown = SCIProwCalcProbability(scip, row, changedrowmean, changedrowvariance, rowinftiesdownafterbranch,
               rowinftiesupafterbranch);
      }
      else
         newrowprobdown = currentrowprob;

      /* update the up and down score depending on the chosen scoring parameter */
      SCIP_CALL( SCIPupdateDistributionScore(scip, currentrowprob, newrowprobup, newrowprobdown, upscore, downscore, scoreparam) );

      SCIPdebugMsg(scip, "  Variable %s changes probability of row %s from %g to %g (branch up) or %g;\n",
         SCIPvarGetName(var), SCIProwGetName(row), currentrowprob, newrowprobup, newrowprobdown);
      SCIPdebugMsg(scip, "  -->  new variable score: %g (for branching up), %g (for branching down)\n",
         *upscore, *downscore);
   }

   return SCIP_OKAY;
}

/** free branchrule data */
static
void branchruledataFreeArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   assert(branchruledata->memsize == 0 || branchruledata->rowmeans != NULL);
   assert(branchruledata->memsize >= 0);

   if( branchruledata->memsize > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &branchruledata->rowmeans, branchruledata->memsize);
      SCIPfreeBlockMemoryArray(scip, &branchruledata->rowvariances, branchruledata->memsize);
      SCIPfreeBlockMemoryArray(scip, &branchruledata->rowinfinitiesup, branchruledata->memsize);
      SCIPfreeBlockMemoryArray(scip, &branchruledata->rowinfinitiesdown, branchruledata->memsize);

      branchruledata->memsize = 0;
   }
}

/** add variable to the bound change event queue; skipped if variable is already in there, or if variable has
 *  no row currently watched
 */
static
void branchruledataAddBoundChangeVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branchrule data */
   SCIP_VAR*             var                 /**< the variable whose bound changes need to be processed */
   )
{
   int varindex;
   int varpos;

   assert(var != NULL);

   varindex = SCIPvarGetProbindex(var);
   assert(-1 <= varindex && varindex < branchruledata->varpossmemsize);

   /* if variable is not active, it should not be watched */
   if( varindex == -1 )
      return;
   varpos = branchruledata->varposs[varindex];
   assert(varpos < branchruledata->nupdatedvars);

   /* nothing to do if variable is already in the queue */
   if( varpos >= 0 )
   {
      assert(branchruledata->updatedvars[varpos] == var);

      return;
   }

   /* if none of the variables rows was calculated yet, variable needs not to be watched */
   assert((branchruledata->currentlbs[varindex] == SCIP_INVALID) == (branchruledata->currentubs[varindex] == SCIP_INVALID)); /*lint !e777*/
   if( branchruledata->currentlbs[varindex] == SCIP_INVALID ) /*lint !e777*/
      return;

   /* add the variable to the branch rule data of variables to process updates for */
   assert(branchruledata->varpossmemsize > branchruledata->nupdatedvars);
   varpos = branchruledata->nupdatedvars;
   branchruledata->updatedvars[varpos] = var;
   branchruledata->varposs[varindex] = varpos;
   ++branchruledata->nupdatedvars;
}

/** returns the next unprocessed variable (last in, first out) with pending bound changes, or NULL */
static
SCIP_VAR* branchruledataPopBoundChangeVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branchrule data */
   )
{
   SCIP_VAR* var;
   int varpos;
   int varindex;

   assert(branchruledata->nupdatedvars >= 0);

   /* return if no variable is currently pending */
   if( branchruledata->nupdatedvars == 0 )
      return NULL;

   varpos = branchruledata->nupdatedvars - 1;
   var = branchruledata->updatedvars[varpos];
   assert(var != NULL);
   varindex = SCIPvarGetProbindex(var);
   assert(0 <= varindex && varindex < branchruledata->varpossmemsize);
   assert(varpos == branchruledata->varposs[varindex]);

   branchruledata->varposs[varindex] = -1;
   branchruledata->nupdatedvars--;

   return var;
}

/** process a variable from the queue of changed variables */
static
SCIP_RETCODE varProcessBoundChanges(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branchrule data */
   SCIP_VAR*             var                 /**< the variable whose bound changes need to be processed */
   )
{
   SCIP_ROW** colrows;
   SCIP_COL* varcol;
   SCIP_Real* colvals;
   SCIP_Real oldmean;
   SCIP_Real newmean;
   SCIP_Real oldvariance;
   SCIP_Real newvariance;
   SCIP_Real oldlb;
   SCIP_Real newlb;
   SCIP_Real oldub;
   SCIP_Real newub;
   SCIP_VARTYPE vartype;
   int ncolrows;
   int r;
   int varindex;

   /* skip event execution if SCIP is in Probing mode because these bound changes will be undone anyway before branching
    * rule is called again
    */
   assert(!SCIPinProbing(scip));

   assert(var != NULL);
   varcol = SCIPvarGetCol(var);
   assert(varcol != NULL);
   colrows = SCIPcolGetRows(varcol);
   colvals = SCIPcolGetVals(varcol);
   ncolrows = SCIPcolGetNNonz(varcol);

   varindex = SCIPvarGetProbindex(var);

   oldlb = branchruledata->currentlbs[varindex];
   oldub = branchruledata->currentubs[varindex];

   /* skip update if the variable has never been subject of previously calculated row activities */
   assert((oldlb == SCIP_INVALID) == (oldub == SCIP_INVALID)); /*lint !e777*/
   if( oldlb == SCIP_INVALID ) /*lint !e777*/
      return SCIP_OKAY;

   newlb = SCIPvarGetLbLocal(var);
   newub = SCIPvarGetUbLocal(var);

   /* skip update if the bound change events have cancelled out */
   if( SCIPisFeasEQ(scip, oldlb, newlb) && SCIPisFeasEQ(scip, oldub, newub) )
      return SCIP_OKAY;

   /* calculate old and new variable distribution mean and variance */
   oldvariance = 0.0;
   newvariance = 0.0;
   oldmean = 0.0;
   newmean = 0.0;
   vartype = SCIPvarGetType(var);
   SCIPvarCalcDistributionParameters(scip, oldlb, oldub, vartype, &oldmean, &oldvariance);
   SCIPvarCalcDistributionParameters(scip, newlb, newub, vartype, &newmean, &newvariance);

   /* loop over all rows of this variable and update activity distribution */
   for( r = 0; r < ncolrows; ++r )
   {
      int rowpos;

      assert(colrows[r] != NULL);
      rowpos = SCIProwGetIndex(colrows[r]);
      assert(rowpos >= 0);

      SCIP_CALL( branchruledataEnsureArraySize(scip, branchruledata, rowpos) );

      /* only consider rows for which activity distribution was already calculated */
      if( branchruledata->rowmeans[rowpos] != SCIP_INVALID ) /*lint !e777*/
      {
         SCIP_Real coeff;
         SCIP_Real coeffsquared;
         assert(branchruledata->rowvariances[rowpos] != SCIP_INVALID
               && SCIPisFeasGE(scip, branchruledata->rowvariances[rowpos], 0.0)); /*lint !e777*/

         coeff = colvals[r];
         coeffsquared = SQUARED(coeff);

         /* update variable contribution to row activity distribution */
         branchruledata->rowmeans[rowpos] += coeff * (newmean - oldmean);
         branchruledata->rowvariances[rowpos] += coeffsquared * (newvariance - oldvariance);
         branchruledata->rowvariances[rowpos] = MAX(0.0, branchruledata->rowvariances[rowpos]);

         /* account for changes of the infinite contributions to row activities */
         if( coeff > 0.0 )
         {
            /* if the coefficient is positive, upper bounds affect activity up */
            if( SCIPisInfinity(scip, newub) && !SCIPisInfinity(scip, oldub) )
               ++branchruledata->rowinfinitiesup[rowpos];
            else if( !SCIPisInfinity(scip, newub) && SCIPisInfinity(scip, oldub) )
               --branchruledata->rowinfinitiesup[rowpos];

            if( SCIPisInfinity(scip, newlb) && !SCIPisInfinity(scip, oldlb) )
               ++branchruledata->rowinfinitiesdown[rowpos];
            else if( !SCIPisInfinity(scip, newlb) && SCIPisInfinity(scip, oldlb) )
               --branchruledata->rowinfinitiesdown[rowpos];
         }
         else if( coeff < 0.0 )
         {
            if( SCIPisInfinity(scip, newub) && !SCIPisInfinity(scip, oldub) )
               ++branchruledata->rowinfinitiesdown[rowpos];
            else if( !SCIPisInfinity(scip, newub) && SCIPisInfinity(scip, oldub) )
               --branchruledata->rowinfinitiesdown[rowpos];

            if( SCIPisInfinity(scip, newlb) && !SCIPisInfinity(scip, oldlb) )
               ++branchruledata->rowinfinitiesup[rowpos];
            else if( !SCIPisInfinity(scip, newlb) && SCIPisInfinity(scip, oldlb) )
               --branchruledata->rowinfinitiesup[rowpos];
         }
         assert(branchruledata->rowinfinitiesdown[rowpos] >= 0);
         assert(branchruledata->rowinfinitiesup[rowpos] >= 0);
      }
   }

   /* store the new local bounds in the data */
   branchruledataUpdateCurrentBounds(scip, branchruledata, var);

   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeDistribution)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyDistribution)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeBranchruleDistribution(scip) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolDistribution)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* free row arrays when branch-and-bound data is freed */
   branchruledataFreeArrays(scip, branchruledata);

   /* drop variable events at the end of branch and bound process (cannot be used after restarts, anyway) */
   if( branchruledata->varfilterposs != NULL)
   {
      SCIP_VAR** vars;
      int nvars;
      int v;

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      assert(nvars > 0);
      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, vars[v], EVENT_DISTRIBUTION, branchruledata->eventhdlr, NULL, branchruledata->varfilterposs[v]) );
      }
      SCIPfreeBlockMemoryArray(scip, &(branchruledata->varfilterposs), nvars);
   }
   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeDistribution)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* free internal arrays first */
   branchruledataFreeArrays(scip, branchruledata);
   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpDistribution)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** lpcands;
   SCIP_VAR* bestcand;
   SCIP_NODE* downchild;
   SCIP_NODE* upchild;

   SCIP_Real* lpcandssol;

   SCIP_Real bestscore;
   SCIP_BRANCHDIR bestbranchdir;
   int nlpcands;
   int c;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, NULL, &nlpcands, NULL) );

   if( nlpcands == 0 )
      return SCIP_OKAY;

   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   branchruledata = SCIPbranchruleGetData(branchrule);

   /* if branching rule data arrays were not initialized before (usually the first call of the branching execution),
    * allocate arrays with an initial capacity of the number of LP rows */
   if( branchruledata->memsize == 0 )
   {
      int nlprows;

      /* get LP rows data */
      nlprows = SCIPgetNLPRows(scip);

      /* initialize arrays only if there are actual LP rows to operate on */
      if( nlprows > 0 )
      {
         SCIP_CALL( branchruledataEnsureArraySize(scip, branchruledata, nlprows) );
      }
      else /* if there are no LP rows, branching rule cannot be used */
         return SCIP_OKAY;
   }

   /* process pending bound change events */
   while( branchruledata->nupdatedvars > 0 )
   {
      SCIP_VAR* nextvar;

      /* pop the next variable from the queue and process its bound changes */
      nextvar = branchruledataPopBoundChangeVar(scip, branchruledata);
      assert(nextvar != NULL);
      SCIP_CALL( varProcessBoundChanges(scip, branchruledata, nextvar) );
   }

   bestscore = -1;
   bestbranchdir = SCIP_BRANCHDIR_AUTO;
   bestcand = NULL;

   /* invalidate the current row distribution data which is initialized on the fly when looping over the candidates */

   /* loop over candidate variables and calculate their score in changing the cumulative
    * probability of fulfilling each of their constraints */
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_Real upscore;
      SCIP_Real downscore;
      SCIP_VAR* lpcand;
      int varindex;

      lpcand = lpcands[c];
      assert(lpcand != NULL);

      varindex = SCIPvarGetProbindex(lpcand);

      /* in debug mode, ensure that all bound process events which occurred in the mean time have been captured
       * by the branching rule event system
       */
      assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(lpcand), SCIPvarGetUbLocal(lpcand)));
      assert(0 <= varindex && varindex < branchruledata->varpossmemsize);

      assert((branchruledata->currentlbs[varindex] == SCIP_INVALID) == (branchruledata->currentubs[varindex] == SCIP_INVALID)); /*lint !e777*/
      assert((branchruledata->currentlbs[varindex] == SCIP_INVALID)
            || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(lpcand), branchruledata->currentlbs[varindex])); /*lint !e777*/
      assert((branchruledata->currentubs[varindex] == SCIP_INVALID)
                  || SCIPisFeasEQ(scip, SCIPvarGetUbLocal(lpcand), branchruledata->currentubs[varindex])); /*lint !e777*/

      /* if the branching rule has not captured the variable bounds yet, this can be done now */
      if( branchruledata->currentlbs[varindex] == SCIP_INVALID ) /*lint !e777*/
      {
         branchruledataUpdateCurrentBounds(scip, branchruledata, lpcand);
      }

      upscore = 0.0;
      downscore = 0.0;

      /* loop over candidate rows and determine the candidate up- and down- branching score w.r.t. the score parameter */
      SCIP_CALL( calcBranchScore(scip, branchruledata, lpcand, lpcandssol[c],
            &upscore, &downscore, branchruledata->scoreparam) );

      /* if weighted scoring is enabled, use the branching score method of SCIP to weigh up and down score */
      if( branchruledata->usescipscore )
      {
         SCIP_Real score;

         score = SCIPgetBranchScore(scip, lpcand, downscore, upscore);

         /* select the candidate with the highest branching score */
         if( score > bestscore )
         {
            bestscore = score;
            bestcand = lpcand;
            /* prioritize branching direction with the higher score */
            if( upscore > downscore )
               bestbranchdir = SCIP_BRANCHDIR_UPWARDS;
            else
               bestbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
         }
      }
      else
      {
         /* no weighted score; keep candidate which has the single highest score in one direction */
         if( upscore > bestscore && upscore > downscore )
         {
            bestscore = upscore;
            bestbranchdir = SCIP_BRANCHDIR_UPWARDS;
            bestcand = lpcand;
         }
         else if( downscore > bestscore )
         {
            bestscore = downscore;
            bestbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
            bestcand = lpcand;
         }
      }

      SCIPdebugMsg(scip, "  Candidate %s has score down %g and up %g \n", SCIPvarGetName(lpcand), downscore, upscore);
      SCIPdebugMsg(scip, "  Best candidate: %s, score %g, direction %d\n", SCIPvarGetName(bestcand), bestscore, bestbranchdir);
   }
   assert(!SCIPisFeasIntegral(scip, SCIPvarGetSol(bestcand, TRUE)));
   assert(bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS || bestbranchdir == SCIP_BRANCHDIR_UPWARDS);
   assert(bestcand != NULL);

   SCIPdebugMsg(scip, "  Branching on variable %s with bounds [%g, %g] and solution value <%g>\n", SCIPvarGetName(bestcand),
      SCIPvarGetLbLocal(bestcand), SCIPvarGetUbLocal(bestcand), SCIPvarGetLPSol(bestcand));

   /* branch on the best candidate variable */
   SCIP_CALL( SCIPbranchVar(scip, bestcand, &downchild, NULL, &upchild) );

   assert(downchild != NULL);
   assert(upchild != NULL);

   if( bestbranchdir == SCIP_BRANCHDIR_UPWARDS )
   {
      SCIP_CALL( SCIPchgChildPrio(scip, upchild, DEFAULT_PRIORITY) );
      SCIPdebugMsg(scip, "  Changing node priority of up-child.\n");
   }
   else
   {
      assert(bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS);
      SCIP_CALL( SCIPchgChildPrio(scip, downchild, DEFAULT_PRIORITY) );
      SCIPdebugMsg(scip, "  Changing node priority of down-child.\n");
   }

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** event execution method of distribution branching which handles bound change events of variables */
static
SCIP_DECL_EVENTEXEC(eventExecDistribution)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_VAR* var;

   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   branchruledata = eventhdlrdata->branchruledata;
   var = SCIPeventGetVar(event);

   /* add the variable to the queue of unprocessed variables; method itself ensures that every variable is added at most once */
   branchruledataAddBoundChangeVar(scip, branchruledata, var);

   return SCIP_OKAY;
}


/*
 * branching rule specific interface methods
 */

/** creates the distribution branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleDistribution(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULE* branchrule = NULL;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create distribution branching rule data */
   branchruledata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );

   branchruledata->memsize = 0;
   branchruledata->rowmeans = NULL;
   branchruledata->rowvariances = NULL;
   branchruledata->rowinfinitiesdown = NULL;
   branchruledata->rowinfinitiesup = NULL;
   branchruledata->varfilterposs = NULL;
   branchruledata->currentlbs = NULL;
   branchruledata->currentubs = NULL;

   /* create event handler first to finish branch rule data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhdlrdata) );
   eventhdlrdata->branchruledata = branchruledata;

   branchruledata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &branchruledata->eventhdlr, EVENTHDLR_NAME,
         "event handler for dynamic acitivity distribution updating",
         eventExecDistribution, eventhdlrdata) );
   assert( branchruledata->eventhdlr != NULL);
   SCIP_CALL( SCIPsetEventhdlrFree(scip, branchruledata->eventhdlr, eventFreeDistribution) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
	 BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyDistribution) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeDistribution) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolDistribution) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpDistribution) );

   /* add distribution branching rule parameters */
   SCIP_CALL( SCIPaddCharParam(scip, "branching/" BRANCHRULE_NAME "/scoreparam",
         "the score;largest 'd'ifference, 'l'owest cumulative probability,'h'ighest c.p., 'v'otes lowest c.p., votes highest c.p.('w') ",
         &branchruledata->scoreparam, TRUE, DEFAULT_SCOREPARAM, SCOREPARAM_VALUES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/" BRANCHRULE_NAME "/onlyactiverows",
         "should only rows which are active at the current node be considered?",
         &branchruledata->onlyactiverows, TRUE, DEFAULT_ONLYACTIVEROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/" BRANCHRULE_NAME "/weightedscore",
         "should the branching score weigh up- and down-scores of a variable",
         &branchruledata->usescipscore, TRUE, DEFAULT_USEWEIGHTEDSCORE, NULL, NULL) );

   return SCIP_OKAY;
}
