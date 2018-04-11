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

/**@file   misc.c
 * @brief  miscellaneous methods
 * @author Tobias Achterberg
 * @author Gerald Gamrath
 * @author Stefan Heinz
 * @author Michael Winkler
 * @author Kati Wolter
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "scip/def.h"
#include "scip/pub_message.h"
#include "scip/misc.h"
#include "scip/intervalarith.h"
#include "scip/pub_misc.h"

#ifndef NDEBUG
#include "scip/struct_misc.h"
#endif

/*
 * methods for statistical tests
 */

#define SQRTOFTWO                  1.4142136 /**< the square root of 2 with sufficient precision */

/**< contains all critical values for a one-sided two sample t-test up to 15 degrees of freedom
 *   a critical value represents a threshold for rejecting the null-hypothesis in hypothesis testing at
 *   a certain confidence level;
 *
 *   access through method SCIPstudentTGetCriticalValue()
 *
 *  source: German Wikipedia
 *
 *  for confidence levels
 *  c =
 *  0.75    0.875     0.90      0.95      0.975 (one-sided)
 *  0.50    0.750     0.80      0.90      0.950 (two-sided)
 *
 */
static const SCIP_Real studentt_quartiles[] = {      /* df:*/
   1.000,    2.414,    3.078,    6.314,    12.706,   /*  1 */
   0.816,    1.604,    1.886,    2.920,    4.303,    /*  2 */
   0.765,    1.423,    1.638,    2.353,    3.182,    /*  3 */
   0.741,    1.344,    1.533,    2.132,    2.776,    /*  4 */
   0.727,    1.301,    1.476,    2.015,    2.571,    /*  5 */
   0.718,    1.273,    1.440,    1.943,    2.447,    /*  6 */
   0.711,    1.254,    1.415,    1.895,    2.365,    /*  7 */
   0.706,    1.240,    1.397,    1.860,    2.306,    /*  8 */
   0.703,    1.230,    1.383,    1.833,    2.262,    /*  9 */
   0.700,    1.221,    1.372,    1.812,    2.228,    /* 10 */
   0.697,    1.214,    1.363,    1.796,    2.201,    /* 11 */
   0.695,    1.209,    1.356,    1.782,    2.179,    /* 12 */
   0.694,    1.204,    1.350,    1.771,    2.160,    /* 13 */
   0.692,    1.200,    1.345,    1.761,    2.145,    /* 14 */
   0.691,    1.197,    1.341,    1.753,    2.131     /* 15 */
};

/**< critical values for higher degrees of freedom of Student-T distribution for the same error probabilities; infact,
 *   these are critical values of the standard normal distribution with mean 0 and variance 1
 */
static const SCIP_Real studentt_quartilesabove[] = {
   0.674,    1.150,    1.282,    1.645,    1.960
};

/** the maximum degrees of freedom represented before switching to normal approximation */
static const int studentt_maxdf = sizeof(studentt_quartiles)/(5 * sizeof(SCIP_Real));

/** get critical value of a Student-T distribution for a given number of degrees of freedom at a confidence level */
SCIP_Real SCIPstudentTGetCriticalValue(
   SCIP_CONFIDENCELEVEL  clevel,             /**< (one-sided) confidence level */
   int                   df                  /**< degrees of freedom */
   )
{
   if( df > studentt_maxdf )
      return studentt_quartilesabove[(int)clevel];
   else
      return studentt_quartiles[(int)clevel + 5 * (df - 1)];
}

/** compute a t-value for the hypothesis that x and y are from the same population; Assuming that
 *  x and y represent normally distributed random samples with equal variance, the returned value
 *  comes from a Student-T distribution with countx + county - 2 degrees of freedom; this
 *  value can be compared with a critical value (see also SCIPstudentTGetCriticalValue()) at
 *  a predefined confidence level for checking if x and y significantly differ in location
 */
SCIP_Real SCIPcomputeTwoSampleTTestValue(
   SCIP_Real             meanx,              /**< the mean of the first distribution */
   SCIP_Real             meany,              /**< the mean of the second distribution */
   SCIP_Real             variancex,          /**< the variance of the x-distribution */
   SCIP_Real             variancey,          /**< the variance of the y-distribution */
   SCIP_Real             countx,             /**< number of samples of x */
   SCIP_Real             county              /**< number of samples of y */
   )
{
   SCIP_Real pooledvariance;
   SCIP_Real tresult;

   /* too few samples */
   if( countx < 1.9 || county < 1.9 )
      return SCIP_INVALID;

   /* pooled variance is the weighted average of the two variances */
   pooledvariance = (countx - 1) * variancex + (county - 1) * variancey;
   pooledvariance /= (countx + county - 2);

   /* a variance close to zero means the distributions are basically constant */
   pooledvariance = MAX(pooledvariance, 1e-9);

   /* tresult can be understood as realization of a Student-T distributed variable with
    * countx + county - 2 degrees of freedom
    */
   tresult = (meanx - meany) / pooledvariance;
   tresult *= SQRT(countx * county / (countx + county));

   return tresult;
}

/** returns the value of the Gauss error function evaluated at a given point */
SCIP_Real SCIPerf(
   SCIP_Real             x                   /**< value to evaluate */
   )
{
#if defined(_WIN32) || defined(_WIN64)
   SCIP_Real a1, a2, a3, a4, a5, p, t, y;
   int sign;

   a1 =  0.254829592;
   a2 = -0.284496736;
   a3 =  1.421413741;
   a4 = -1.453152027;
   a5 =  1.061405429;
   p  =  0.3275911;

   sign = (x >= 0) ? 1 : -1;
   x = REALABS(x);

   t = 1.0/(1.0 + p*x);
   y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
   return sign * y;
#else
   return erf(x);
#endif
}

/** get critical value of a standard normal distribution  at a given confidence level */
SCIP_Real SCIPnormalGetCriticalValue(
   SCIP_CONFIDENCELEVEL  clevel              /**< (one-sided) confidence level */
   )
{
   return studentt_quartilesabove[(int)clevel];
}

/** calculates the cumulative distribution P(-infinity <= x <= value) that a normally distributed
 *  random variable x takes a value between -infinity and parameter \p value.
 *
 *  The distribution is given by the respective mean and deviation. This implementation
 *  uses the error function SCIPerf().
 */
SCIP_Real SCIPnormalCDF(
   SCIP_Real             mean,               /**< the mean value of the distribution */
   SCIP_Real             variance,           /**< the square of the deviation of the distribution */
   SCIP_Real             value               /**< the upper limit of the calculated distribution integral */
   )
{
   SCIP_Real normvalue;
   SCIP_Real std;

   /* we need to calculate the standard deviation from the variance */
   assert(variance >= -1e-9);
   if( variance < 1e-9 )
      std = 0.0;
   else
      std = sqrt(variance);

   /* special treatment for zero variance */
   if( std < 1e-9 )
   {
      if( value < mean + 1e-9 )
         return 1.0;
      else
         return 0.0;
   }
   assert( std != 0.0 ); /* for lint */

   /* scale and translate to standard normal distribution. Factor sqrt(2) is needed for SCIPerf() function */
   normvalue = (value - mean)/(std * SQRTOFTWO);

   SCIPdebugMessage(" Normalized value %g = ( %g - %g ) / (%g * 1.4142136)\n", normvalue, value, mean, std);

   /* calculate the cumulative distribution function for normvalue. For negative normvalues, we negate the normvalue and
    * use the oddness of the SCIPerf()-function; special treatment for values close to zero.
    */
   if( normvalue < 1e-9 && normvalue > -1e-9 )
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

/*
 * SCIP regression methods
 */

/** returns the number of observations of this regression */
int SCIPregressionGetNObservations(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   )
{
   assert(regression != NULL);

   return regression->nobservations;
}

/** return the current slope of the regression */
SCIP_Real SCIPregressionGetSlope(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   )
{
   assert(regression != NULL);

   return regression->slope;
}

/** get the current y-intercept of the regression */
SCIP_Real SCIPregressionGetIntercept(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   )
{
   assert(regression != NULL);

   return regression->intercept;
}

/** recomputes regression coefficients from available observation data */
static
void regressionRecompute(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   )
{
   /* regression coefficients require two or more observations and variance in x */
   if( regression->nobservations <= 1 || EPSZ(regression->variancesumx, 1e-9) )
   {
      regression->slope = SCIP_INVALID;
      regression->intercept = SCIP_INVALID;
      regression->corrcoef = SCIP_INVALID;
   }
   else if( EPSZ(regression->variancesumy, 1e-9) )
   {
      /* if there is no variance in the y's (but in the x's), the regression line is horizontal with y-intercept through the mean y */
      regression->slope = 0.0;
      regression->corrcoef = 0.0;
      regression->intercept = regression->meany;
   }
   else
   {
      /* we ruled this case out already, but to please some compilers... */
      assert(regression->variancesumx > 0.0);
      assert(regression->variancesumy > 0.0);

      /* compute slope */
      regression->slope = (regression->sumxy  - regression->nobservations * regression->meanx * regression->meany) / regression->variancesumx;

      /* compute y-intercept */
      regression->intercept = regression->meany - regression->slope * regression->meanx;

      /* compute empirical correlation coefficient */
      regression->corrcoef = (regression->sumxy - regression->nobservations * regression->meanx * regression->meany) /
         sqrt(regression->variancesumx * regression->variancesumy);
   }
}

/* incremental update of statistics describing mean and variance */
static
void incrementalStatsUpdate(
   SCIP_Real             value,              /**< current value to be added to incremental statistics */
   SCIP_Real*            meanptr,            /**< pointer to value of current mean */
   SCIP_Real*            sumvarptr,          /**< pointer to the value of the current variance sum term */
   int                   nobservations,      /**< total number of observations */
   SCIP_Bool             add                 /**< TRUE if the value should be added, FALSE for removing it */
   )
{
   SCIP_Real oldmean;
   SCIP_Real addfactor;
   assert(meanptr != NULL);
   assert(sumvarptr != NULL);
   assert(nobservations > 0 || add);

   addfactor = add ? 1.0 : -1.0;

   oldmean = *meanptr;
   *meanptr = oldmean + addfactor * (value - oldmean)/(SCIP_Real)nobservations;
   *sumvarptr += addfactor * (value - oldmean) * (value - (*meanptr));

   /* it may happen that *sumvarptr is slightly negative, especially after a series of add/removal operations */
   assert(*sumvarptr >= -1e-6);
   *sumvarptr = MAX(0.0, *sumvarptr);
}

/** removes an observation (x,y) from the regression */
void SCIPregressionRemoveObservation(
   SCIP_REGRESSION*      regression,         /**< regression data structure */
   SCIP_Real             x,                  /**< X of observation */
   SCIP_Real             y                   /**< Y of the observation */
   )
{
   assert(regression != NULL);
   assert(regression->nobservations > 0);

   /* simply call the reset function in the case of a single remaining observation to avoid numerical troubles */
   if( regression->nobservations == 1 )
   {
      SCIPregressionReset(regression);
   }
   else
   {
      SCIP_Bool add = FALSE;
      --regression->nobservations;

      /* decrement individual means and variances */
      incrementalStatsUpdate(x, &regression->meanx, &regression->variancesumx, regression->nobservations, add);
      incrementalStatsUpdate(y, &regression->meany, &regression->variancesumy, regression->nobservations, add);

      /* decrement product sum */
      regression->sumxy -= (x * y);
   }

   /* recompute regression parameters */
   regressionRecompute(regression);
}

/** update regression by a new observation (x,y) */
void SCIPregressionAddObservation(
   SCIP_REGRESSION*      regression,         /**< regression data structure */
   SCIP_Real             x,                  /**< X of observation */
   SCIP_Real             y                   /**< Y of the observation */
   )
{
   SCIP_Bool add = TRUE;
   assert(regression != NULL);

   ++(regression->nobservations);
   incrementalStatsUpdate(x, &regression->meanx, &regression->variancesumx, regression->nobservations, add);
   incrementalStatsUpdate(y, &regression->meany, &regression->variancesumy, regression->nobservations, add);

   regression->sumxy += x * y;

   regressionRecompute(regression);
}

/** reset regression data structure */
void SCIPregressionReset(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   )
{
   regression->intercept = SCIP_INVALID;
   regression->slope = SCIP_INVALID;
   regression->corrcoef = SCIP_INVALID;
   regression->meanx = 0;
   regression->variancesumx = 0;
   regression->sumxy = 0;
   regression->meany = 0;
   regression->variancesumy = 0;
   regression->nobservations = 0;
}

/** creates and resets a regression */
SCIP_RETCODE SCIPregressionCreate(
   SCIP_REGRESSION**     regression          /**< regression data structure */
   )
{
   assert(regression != NULL);

   /* allocate necessary memory */
   SCIP_ALLOC (BMSallocMemory(regression) );

   /* reset the regression */
   SCIPregressionReset(*regression);

   return SCIP_OKAY;
}

/** creates and resets a regression */
void SCIPregressionFree(
   SCIP_REGRESSION**     regression          /**< regression data structure */
   )
{
   BMSfreeMemory(regression);
}

/** calculate memory size for dynamically allocated arrays (copied from scip/set.c) */
static
int calcGrowSize(
   int                   initsize,           /**< initial size of array */
   SCIP_Real             growfac,            /**< growing factor of array */
   int                   num                 /**< minimum number of entries to store */
   )
{
   int size;

   assert(initsize >= 0);
   assert(growfac >= 1.0);
   assert(num >= 0);

   if( growfac == 1.0 )
      size = MAX(initsize, num);
   else
   {
      int oldsize;

      /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
      initsize = MAX(initsize, 4);
      size = initsize;
      oldsize = size - 1;

      /* second condition checks against overflow */
      while( size < num && size > oldsize )
      {
         oldsize = size;
         size = (int)(growfac * size + initsize);
      }

      /* if an overflow happened, set the correct value */
      if( size <= oldsize )
         size = num;
   }

   assert(size >= initsize);
   assert(size >= num);

   return size;
}

/*
 * GML graphical printing methods
 * For a detailed format decription see http://docs.yworks.com/yfiles/doc/developers-guide/gml.html
 */

#define GMLNODEWIDTH 120.0
#define GMLNODEHEIGTH 30.0
#define GMLFONTSIZE 13
#define GMLNODETYPE "rectangle"
#define GMLNODEFILLCOLOR "#ff0000"
#define GMLEDGECOLOR "black"
#define GMLNODEBORDERCOLOR "#000000"


/** writes a node section to the given graph file */
void SCIPgmlWriteNode(
   FILE*                 file,               /**< file to write to */
   unsigned int          id,                 /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node, or NULL */
   const char*           fillcolor,          /**< color of the node's interior, or NULL */
   const char*           bordercolor         /**< color of the node's border, or NULL */
   )
{
   assert(file != NULL);
   assert(label != NULL);

   fprintf(file, "  node\n");
   fprintf(file, "  [\n");
   fprintf(file, "    id      %u\n", id);
   fprintf(file, "    label   \"%s\"\n", label);
   fprintf(file, "    graphics\n");
   fprintf(file, "    [\n");
   fprintf(file, "      w       %g\n", GMLNODEWIDTH);
   fprintf(file, "      h       %g\n", GMLNODEHEIGTH);

   if( nodetype != NULL )
      fprintf(file, "      type    \"%s\"\n", nodetype);
   else
      fprintf(file, "      type    \"%s\"\n", GMLNODETYPE);

   if( fillcolor != NULL )
      fprintf(file, "      fill    \"%s\"\n", fillcolor);
   else
      fprintf(file, "      fill    \"%s\"\n", GMLNODEFILLCOLOR);

   if( bordercolor != NULL )
      fprintf(file, "      outline \"%s\"\n", bordercolor);
   else
      fprintf(file, "      outline \"%s\"\n", GMLNODEBORDERCOLOR);

   fprintf(file, "    ]\n");
   fprintf(file, "    LabelGraphics\n");
   fprintf(file, "    [\n");
   fprintf(file, "      text      \"%s\"\n", label);
   fprintf(file, "      fontSize  %d\n", GMLFONTSIZE);
   fprintf(file, "      fontName  \"Dialog\"\n");
   fprintf(file, "      anchor    \"c\"\n");
   fprintf(file, "    ]\n");
   fprintf(file, "  ]\n");
}

/** writes a node section including weight to the given graph file */
void SCIPgmlWriteNodeWeight(
   FILE*                 file,               /**< file to write to */
   unsigned int          id,                 /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node, or NULL */
   const char*           fillcolor,          /**< color of the node's interior, or NULL */
   const char*           bordercolor,        /**< color of the node's border, or NULL */
   SCIP_Real             weight              /**< weight of node */
   )
{
   assert(file != NULL);
   assert(label != NULL);

   fprintf(file, "  node\n");
   fprintf(file, "  [\n");
   fprintf(file, "    id      %u\n", id);
   fprintf(file, "    label   \"%s\"\n", label);
   fprintf(file, "    weight  %g\n", weight);
   fprintf(file, "    graphics\n");
   fprintf(file, "    [\n");
   fprintf(file, "      w       %g\n", GMLNODEWIDTH);
   fprintf(file, "      h       %g\n", GMLNODEHEIGTH);

   if( nodetype != NULL )
      fprintf(file, "      type    \"%s\"\n", nodetype);
   else
      fprintf(file, "      type    \"%s\"\n", GMLNODETYPE);

   if( fillcolor != NULL )
      fprintf(file, "      fill    \"%s\"\n", fillcolor);
   else
      fprintf(file, "      fill    \"%s\"\n", GMLNODEFILLCOLOR);

   if( bordercolor != NULL )
      fprintf(file, "      outline \"%s\"\n", bordercolor);
   else
      fprintf(file, "      outline \"%s\"\n", GMLNODEBORDERCOLOR);

   fprintf(file, "    ]\n");
   fprintf(file, "    LabelGraphics\n");
   fprintf(file, "    [\n");
   fprintf(file, "      text      \"%s\"\n", label);
   fprintf(file, "      fontSize  %d\n", GMLFONTSIZE);
   fprintf(file, "      fontName  \"Dialog\"\n");
   fprintf(file, "      anchor    \"c\"\n");
   fprintf(file, "    ]\n");
   fprintf(file, "  ]\n");
}

/** writes an edge section to the given graph file */
void SCIPgmlWriteEdge(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   )
{
   assert(file != NULL);

   fprintf(file, "  edge\n");
   fprintf(file, "  [\n");
   fprintf(file, "    source  %u\n", source);
   fprintf(file, "    target  %u\n", target);

   if( label != NULL)
      fprintf(file, "    label   \"%s\"\n", label);

   fprintf(file, "    graphics\n");
   fprintf(file, "    [\n");

   if( color != NULL )
      fprintf(file, "      fill    \"%s\"\n", color);
   else
      fprintf(file, "      fill    \"%s\"\n", GMLEDGECOLOR);

   /* fprintf(file, "      arrow     \"both\"\n"); */
   fprintf(file, "    ]\n");

   if( label != NULL)
   {
      fprintf(file, "    LabelGraphics\n");
      fprintf(file, "    [\n");
      fprintf(file, "      text      \"%s\"\n", label);
      fprintf(file, "      fontSize  %d\n", GMLFONTSIZE);
      fprintf(file, "      fontName  \"Dialog\"\n");
      fprintf(file, "      anchor    \"c\"\n");
      fprintf(file, "    ]\n");
   }

   fprintf(file, "  ]\n");
}

/** writes an arc section to the given graph file */
void SCIPgmlWriteArc(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   )
{
   assert(file != NULL);

   fprintf(file, "  edge\n");
   fprintf(file, "  [\n");
   fprintf(file, "    source  %u\n", source);
   fprintf(file, "    target  %u\n", target);

   if( label != NULL)
      fprintf(file, "    label   \"%s\"\n", label);

   fprintf(file, "    graphics\n");
   fprintf(file, "    [\n");

   if( color != NULL )
      fprintf(file, "      fill    \"%s\"\n", color);
   else
      fprintf(file, "      fill    \"%s\"\n", GMLEDGECOLOR);

   fprintf(file, "      targetArrow     \"standard\"\n");
   fprintf(file, "    ]\n");

   if( label != NULL)
   {
      fprintf(file, "    LabelGraphics\n");
      fprintf(file, "    [\n");
      fprintf(file, "      text      \"%s\"\n", label);
      fprintf(file, "      fontSize  %d\n", GMLFONTSIZE);
      fprintf(file, "      fontName  \"Dialog\"\n");
      fprintf(file, "      anchor    \"c\"\n");
      fprintf(file, "    ]\n");
   }

   fprintf(file, "  ]\n");
}

/** writes the starting line to a GML graph file, does not open a file */
void SCIPgmlWriteOpening(
   FILE*                 file,               /**< file to write to */
   SCIP_Bool             directed            /**< is the graph directed */
   )
{
   assert(file != NULL);

   fprintf(file, "graph\n");
   fprintf(file, "[\n");
   fprintf(file, "  hierarchic      1\n");

   if( directed )
      fprintf(file, "  directed        1\n");
}

/** writes the ending lines to a GML graph file, does not close a file */
void SCIPgmlWriteClosing(
   FILE*                 file                /**< file to close */
   )
{
   assert(file != NULL);

   fprintf(file, "]\n");
}


/*
 * Sparse solution
 */

/** creates a sparse solution */
SCIP_RETCODE SCIPsparseSolCreate(
   SCIP_SPARSESOL**      sparsesol,          /**< pointer to store the created sparse solution */
   SCIP_VAR**            vars,               /**< variables in the sparse solution, must not contain continuous
					      *   variables
					      */
   int                   nvars,              /**< number of variables to store, size of the lower and upper bound
					      *   arrays
					      */
   SCIP_Bool             cleared             /**< should the lower and upper bound arrays be cleared (entries set to
					      *	  0)
					      */
   )
{
   assert(sparsesol != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);

   SCIP_ALLOC( BMSallocMemory(sparsesol) );

#ifndef NDEBUG
   {
      int v;

      for( v = nvars - 1; v >= 0; --v )
      {
	 assert(vars[v] != NULL);
	 /* assert(SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS); */
      }
   }
#endif

   /* copy variables */
   SCIP_ALLOC( BMSduplicateMemoryArray(&((*sparsesol)->vars), vars, nvars) );

   /* create bound arrays */
   if( cleared )
   {
      SCIP_ALLOC( BMSallocClearMemoryArray(&((*sparsesol)->lbvalues), nvars) );
      SCIP_ALLOC( BMSallocClearMemoryArray(&((*sparsesol)->ubvalues), nvars) );
   }
   else
   {
      SCIP_ALLOC( BMSallocMemoryArray(&((*sparsesol)->lbvalues), nvars) );
      SCIP_ALLOC( BMSallocMemoryArray(&((*sparsesol)->ubvalues), nvars) );
   }

   (*sparsesol)->nvars = nvars;

   return SCIP_OKAY;
}

/** frees sparse solution */
void SCIPsparseSolFree(
   SCIP_SPARSESOL**      sparsesol           /**< pointer to a sparse solution */
   )
{
   assert(sparsesol != NULL);
   assert(*sparsesol != NULL);

   BMSfreeMemoryArray(&((*sparsesol)->vars));
   BMSfreeMemoryArray(&((*sparsesol)->ubvalues));
   BMSfreeMemoryArray(&((*sparsesol)->lbvalues));
   BMSfreeMemory(sparsesol);
}

/** returns the variables stored in the given sparse solution */
SCIP_VAR** SCIPsparseSolGetVars(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   )
{
   assert(sparsesol != NULL);

   return sparsesol->vars;
}

/** returns the number of variables stored in the given sparse solution */
int SCIPsparseSolGetNVars(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   )
{
   assert(sparsesol != NULL);

   return sparsesol->nvars;
}

/** returns the lower bound array for all variables for a given sparse solution */
SCIP_Longint* SCIPsparseSolGetLbs(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   )
{
   assert(sparsesol != NULL);

   return sparsesol->lbvalues;
}

/** returns the upper bound array for all variables for a given sparse solution */
SCIP_Longint* SCIPsparseSolGetUbs(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   )
{
   assert(sparsesol != NULL);

   return sparsesol->ubvalues;
}

/** constructs the first solution of sparse solution (all variables are set to their lower bound value */
void SCIPsparseSolGetFirstSol(
   SCIP_SPARSESOL*       sparsesol,          /**< sparse solutions */
   SCIP_Longint*         sol,                /**< array to store the first solution */
   int                   nvars               /**< number of variables */
   )
{
   SCIP_Longint* lbvalues;
   int v;

   assert(sparsesol != NULL);
   assert(sol != NULL);
   assert(nvars == SCIPsparseSolGetNVars(sparsesol));

   lbvalues = SCIPsparseSolGetLbs(sparsesol);
   assert(lbvalues != NULL);

   /* copy the lower bounds */
   for( v = 0; v < nvars; ++v )
      sol[v] = lbvalues[v];
}


/** constructs the next solution of the sparse solution and return whether there was one more or not */
SCIP_Bool SCIPsparseSolGetNextSol(
   SCIP_SPARSESOL*       sparsesol,          /**< sparse solutions */
   SCIP_Longint*         sol,                /**< current solution array which get changed to the next solution */
   int                   nvars               /**< number of variables */
   )
{
   SCIP_Longint* lbvalues;
   SCIP_Longint* ubvalues;
   SCIP_Longint lbvalue;
   SCIP_Longint ubvalue;
   SCIP_Bool singular;
   SCIP_Bool carryflag;
   int v;

   assert(sparsesol != NULL);
   assert(sol != NULL);

   if( nvars == 0 )
      return FALSE;

   assert(nvars > 0);
   assert(nvars == SCIPsparseSolGetNVars(sparsesol));

   lbvalues = SCIPsparseSolGetLbs(sparsesol);
   ubvalues = SCIPsparseSolGetUbs(sparsesol);
   assert(lbvalues != NULL);
   assert(ubvalues != NULL);

   singular = TRUE;
   carryflag = FALSE;

   for( v = 0; v < nvars; ++v )
   {
      lbvalue = lbvalues[v];
      ubvalue = ubvalues[v];

      if( lbvalue < ubvalue )
      {
         singular = FALSE;

         if( carryflag == FALSE )
         {
            if( sol[v] < ubvalue )
            {
               sol[v]++;
               break;
            }
            else
            {
               /* in the last solution the variables v was set to its upper bound value */
               assert(sol[v] == ubvalue);
               sol[v] = lbvalue;
               carryflag = TRUE;
            }
         }
         else
         {
            if( sol[v] < ubvalue )
            {
               sol[v]++;
               carryflag = FALSE;
               break;
            }
            else
            {
               assert(sol[v] == ubvalue);
               sol[v] = lbvalue;
            }
         }
      }
   }

   return (!carryflag && !singular);
}


/*
 * Queue
 */

/** resizes element memory to hold at least the given number of elements */
static
SCIP_RETCODE queueResize(
   SCIP_QUEUE*           queue,              /**< pointer to a queue */
   int                   minsize             /**< minimal number of storable elements */
   )
{
   assert(queue != NULL);
   assert(minsize > 0);

   if( minsize <= queue->size )
      return SCIP_OKAY;

   queue->size = MAX(minsize, (int)(queue->size * queue->sizefac));
   SCIP_ALLOC( BMSreallocMemoryArray(&queue->slots, queue->size) );

   return SCIP_OKAY;
}


/** creates a (circular) queue, best used if the size will be fixed or will not be increased that much */
SCIP_RETCODE SCIPqueueCreate(
   SCIP_QUEUE**          queue,              /**< pointer to the new queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac             /**< memory growing factor applied, if more element slots are needed */
   )
{
   assert(queue != NULL);

   initsize = MAX(1, initsize);
   sizefac = MAX(1.0, sizefac);

   SCIP_ALLOC( BMSallocMemory(queue) );
   (*queue)->firstfree = 0;
   (*queue)->firstused = -1;
   (*queue)->size = 0;
   (*queue)->sizefac = sizefac;
   (*queue)->slots = NULL;

   SCIP_CALL( queueResize(*queue, initsize) );

   return SCIP_OKAY;
}

/** frees queue, but not the data elements themselves */
void SCIPqueueFree(
   SCIP_QUEUE**          queue               /**< pointer to a queue */
   )
{
   assert(queue != NULL);

   BMSfreeMemoryArray(&(*queue)->slots);
   BMSfreeMemory(queue);
}

/** clears the queue, but doesn't free the data elements themselves */
void SCIPqueueClear(
   SCIP_QUEUE*           queue               /**< queue */
   )
{
   assert(queue != NULL);

   queue->firstfree = 0;
   queue->firstused = -1;
}

/** inserts element at the end of the queue */
SCIP_RETCODE SCIPqueueInsert(
   SCIP_QUEUE*           queue,              /**< queue */
   void*                 elem                /**< element to be inserted */
   )
{
   assert(queue != NULL);
   assert(queue->slots != NULL);
   assert(queue->firstused >= -1 && queue->firstused < queue->size);
   assert(queue->firstfree >= 0 && queue->firstused < queue->size);
   assert(queue->firstused > -1 || queue->firstfree == 0);
   assert(elem != NULL);

   if( queue->firstfree == queue->firstused )
   {
      int sizediff;
      int oldsize = queue->size;

      SCIP_CALL( queueResize(queue, queue->size+1) );
      assert(oldsize < queue->size);

      sizediff = queue->size - oldsize;

      /* move the used memory at the slots to the end */
      BMSmoveMemoryArray(&(queue->slots[queue->firstused + sizediff]), &(queue->slots[queue->firstused]), oldsize - queue->firstused); /*lint !e866*/
      queue->firstused += sizediff;
   }
   assert(queue->firstfree != queue->firstused);

   /* insert element as leaf in the tree, move it towards the root as long it is better than its parent */
   queue->slots[queue->firstfree] = elem;
   ++(queue->firstfree);

   /* if we saved the value at the last position we need to reset the firstfree position */
   if( queue->firstfree == queue->size )
      queue->firstfree = 0;

   /* if a first element was added, we need to update the firstused counter */
   if( queue->firstused == -1 )
      queue->firstused = 0;

   return SCIP_OKAY;
}

/** removes and returns the first element of the queue */
void* SCIPqueueRemove(
   SCIP_QUEUE*           queue               /**< queue */
   )
{
   int pos;

   assert(queue != NULL);
   assert(queue->firstused >= -1 && queue->firstused < queue->size);
   assert(queue->firstfree >= 0 && queue->firstused < queue->size);
   assert(queue->firstused > -1 || queue->firstfree == 0);

   if( queue->firstused == -1 )
      return NULL;

   assert(queue->slots != NULL);

   pos = queue->firstused;
   ++(queue->firstused);

   /* if we removed the value at the last position we need to reset the firstused position */
   if( queue->firstused == queue->size )
      queue->firstused = 0;

   /* if we reached the first free position we can reset both, firstused and firstused, positions */
   if( queue->firstused == queue->firstfree )
   {
      queue->firstused = -1;
      queue->firstfree = 0; /* this is not necessary but looks better if we have an empty list to reset this value */
   }

   return (queue->slots[pos]);
}

/** returns the first element of the queue without removing it */
void* SCIPqueueFirst(
   SCIP_QUEUE*           queue               /**< queue */
   )
{
   assert(queue != NULL);
   assert(queue->firstused >= -1 && queue->firstused < queue->size);
   assert(queue->firstfree >= 0 && queue->firstused < queue->size);
   assert(queue->firstused > -1 || queue->firstfree == 0);

   if( queue->firstused == -1 )
      return NULL;

   assert(queue->slots != NULL);

   return queue->slots[queue->firstused];
}

/** returns whether the queue is empty */
SCIP_Bool SCIPqueueIsEmpty(
   SCIP_QUEUE*           queue               /**< queue */
   )
{
   assert(queue != NULL);
   assert(queue->firstused >= -1 && queue->firstused < queue->size);
   assert(queue->firstfree >= 0 && queue->firstused < queue->size);
   assert(queue->firstused > -1 || queue->firstfree == 0);

   return (queue->firstused == -1);
}

/** returns the number of elements in the queue */
int SCIPqueueNElems(
   SCIP_QUEUE*           queue               /**< queue */
   )
{
   assert(queue != NULL);
   assert(queue->firstused >= -1 && queue->firstused < queue->size);
   assert(queue->firstfree >= 0 && queue->firstused < queue->size);
   assert(queue->firstused > -1 || queue->firstfree == 0);

   if( queue->firstused == -1 )
      return 0;
   else if( queue->firstused < queue->firstfree )
      return queue->firstfree - queue->firstused;
   else if( queue->firstused == queue->firstfree )
      return queue->size;
   else
      return queue->firstfree + (queue->size - queue->firstused);
}


/*
 * Priority Queue
 */

#define PQ_PARENT(q) (((q)+1)/2-1)
#define PQ_LEFTCHILD(p) (2*(p)+1)
#define PQ_RIGHTCHILD(p) (2*(p)+2)


/** resizes element memory to hold at least the given number of elements */
static
SCIP_RETCODE pqueueResize(
   SCIP_PQUEUE*          pqueue,             /**< pointer to a priority queue */
   int                   minsize             /**< minimal number of storable elements */
   )
{
   assert(pqueue != NULL);

   if( minsize <= pqueue->size )
      return SCIP_OKAY;

   pqueue->size = MAX(minsize, (int)(pqueue->size * pqueue->sizefac));
   SCIP_ALLOC( BMSreallocMemoryArray(&pqueue->slots, pqueue->size) );

   return SCIP_OKAY;
}

/** creates priority queue */
SCIP_RETCODE SCIPpqueueCreate(
   SCIP_PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   assert(pqueue != NULL);
   assert(ptrcomp != NULL);

   initsize = MAX(1, initsize);
   sizefac = MAX(1.0, sizefac);

   SCIP_ALLOC( BMSallocMemory(pqueue) );
   (*pqueue)->len = 0;
   (*pqueue)->size = 0;
   (*pqueue)->sizefac = sizefac;
   (*pqueue)->slots = NULL;
   (*pqueue)->ptrcomp = ptrcomp;
   SCIP_CALL( pqueueResize(*pqueue, initsize) );

   return SCIP_OKAY;
}

/** frees priority queue, but not the data elements themselves */
void SCIPpqueueFree(
   SCIP_PQUEUE**         pqueue              /**< pointer to a priority queue */
   )
{
   assert(pqueue != NULL);

   BMSfreeMemoryArray(&(*pqueue)->slots);
   BMSfreeMemory(pqueue);
}

/** clears the priority queue, but doesn't free the data elements themselves */
void SCIPpqueueClear(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);

   pqueue->len = 0;
}

/** inserts element into priority queue */
SCIP_RETCODE SCIPpqueueInsert(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   )
{
   int pos;

   assert(pqueue != NULL);
   assert(pqueue->len >= 0);
   assert(elem != NULL);

   SCIP_CALL( pqueueResize(pqueue, pqueue->len+1) );

   /* insert element as leaf in the tree, move it towards the root as long it is better than its parent */
   pos = pqueue->len;
   pqueue->len++;
   while( pos > 0 && (*pqueue->ptrcomp)(elem, pqueue->slots[PQ_PARENT(pos)]) < 0 )
   {
      pqueue->slots[pos] = pqueue->slots[PQ_PARENT(pos)];
      pos = PQ_PARENT(pos);
   }
   pqueue->slots[pos] = elem;

   return SCIP_OKAY;
}

/** removes and returns best element from the priority queue */
void* SCIPpqueueRemove(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   void* root;
   void* last;
   int pos;
   int childpos;
   int brotherpos;

   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   if( pqueue->len == 0 )
      return NULL;

   /* remove root element of the tree, move the better child to its parents position until the last element
    * of the queue could be placed in the empty slot
    */
   root = pqueue->slots[0];
   last = pqueue->slots[pqueue->len-1];
   pqueue->len--;
   pos = 0;
   while( pos <= PQ_PARENT(pqueue->len-1) )
   {
      childpos = PQ_LEFTCHILD(pos);
      brotherpos = PQ_RIGHTCHILD(pos);
      if( brotherpos <= pqueue->len && (*pqueue->ptrcomp)(pqueue->slots[brotherpos], pqueue->slots[childpos]) < 0 )
         childpos = brotherpos;
      if( (*pqueue->ptrcomp)(last, pqueue->slots[childpos]) <= 0 )
         break;
      pqueue->slots[pos] = pqueue->slots[childpos];
      pos = childpos;
   }
   assert(pos <= pqueue->len);
   pqueue->slots[pos] = last;

   return root;
}

/** returns the best element of the queue without removing it */
void* SCIPpqueueFirst(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   if( pqueue->len == 0 )
      return NULL;

   return pqueue->slots[0];
}

/** returns the number of elements in the queue */
int SCIPpqueueNElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   return pqueue->len;
}

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
void** SCIPpqueueElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   return pqueue->slots;
}




/*
 * Hash Table
 */

/** table of some prime numbers */
static int primetable[] = {
   2,
   7,
   19,
   31,
   59,
   227,
   617,
   1523,
   3547,
   8011,
   17707,
   38723,
   83833,
   180317,
   385897,
   821411,
   1742369,
   3680893,
   5693959,
   7753849,
   9849703,
   11973277,
   14121853,
   17643961,
   24273817,
   32452843,
   49979687,
   67867967,
   86028121,
   104395301,
   122949823,
   141650939,
   160481183,
   179424673,
   198491317,
   217645177,
   256203161,
   314606869,
   373587883,
   433024223,
   492876847,
   553105243,
   613651349,
   694847533,
   756065159,
   817504243,
   879190747,
   941083981,
   982451653,
   INT_MAX
};
static const int primetablesize = sizeof(primetable)/sizeof(int);

/** simple and fast 2-universal hash function using multiply and shift */
static
uint32_t hashvalue(
   uint64_t              input               /**< key value */
   )
{
   return ( (uint32_t) ((UINT64_C(0x9e3779b97f4a7c15) * input)>>32) ) | 1u;
}

/** returns a reasonable hash table size (a prime number) that is at least as large as the specified value */
int SCIPcalcMultihashSize(
   int                   minsize             /**< minimal size of the hash table */
   )
{
   int pos;

   (void) SCIPsortedvecFindInt(primetable, minsize, primetablesize, &pos);
   assert(pos < primetablesize);

   return primetable[pos];
}

/** appends element to the multihash list */
static
SCIP_RETCODE multihashlistAppend(
   SCIP_MULTIHASHLIST**  multihashlist,      /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 element             /**< element to append to the list */
   )
{
   SCIP_MULTIHASHLIST* newlist;

   assert(multihashlist != NULL);
   assert(blkmem != NULL);
   assert(element != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &newlist) );
   newlist->element = element;
   newlist->next = *multihashlist;
   *multihashlist = newlist;

   return SCIP_OKAY;
}

/** frees a multihash list entry and all its successors */
static
void multihashlistFree(
   SCIP_MULTIHASHLIST**  multihashlist,      /**< pointer to multihash list to free */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_MULTIHASHLIST* list;
   SCIP_MULTIHASHLIST* nextlist;

   assert(multihashlist != NULL);
   assert(blkmem != NULL);

   list = *multihashlist;
   while( list != NULL )
   {
      nextlist = list->next;
      BMSfreeBlockMemory(blkmem, &list);
      list = nextlist;
   }

   *multihashlist = NULL;
}

/** finds multihash list entry pointing to element with given key in the multihash list, returns NULL if not found */
static
SCIP_MULTIHASHLIST* multihashlistFind(
   SCIP_MULTIHASHLIST*   multihashlist,      /**< multihash list */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   uint64_t              keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   uint64_t currentkeyval;
   void* currentkey;

   assert(hashkeyeq != NULL);
   assert(key != NULL);

   while( multihashlist != NULL )
   {
      currentkey = hashgetkey(userptr, multihashlist->element);
      currentkeyval = hashkeyval(userptr, currentkey);
      if( currentkeyval == keyval && hashkeyeq(userptr, currentkey, key) )
         return multihashlist;

      multihashlist = multihashlist->next;
   }

   return NULL;
}

/** retrieves element with given key from the multihash list, or NULL */
static
void* multihashlistRetrieve(
   SCIP_MULTIHASHLIST*   multihashlist,      /**< hash list */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   uint64_t              keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   SCIP_MULTIHASHLIST* h;

   /* find hash list entry */
   h = multihashlistFind(multihashlist, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key);

   /* return element */
   if( h != NULL )
   {
#ifndef NDEBUG
      SCIP_MULTIHASHLIST* h2;

      h2 = multihashlistFind(h->next, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key);

      if( h2 != NULL )
      {
         void* key1;
         void* key2;

         key1 = hashgetkey(userptr, h->element);
         key2 = hashgetkey(userptr, h2->element);
         assert(hashkeyval(userptr, key1) == hashkeyval(userptr, key2));

         if( hashkeyeq(userptr, key1, key2) )
         {
            SCIPerrorMessage("WARNING: hashkey with same value exists multiple times (e.g. duplicate constraint/variable names), so the return value is maybe not correct\n");
         }
      }
#endif

      return h->element;
   }
   else
      return NULL;
}


/** retrieves element with given key from the multihash list, or NULL
 *  returns pointer to multihash table list entry
 */
static
void* multihashlistRetrieveNext(
   SCIP_MULTIHASHLIST**  multihashlist,      /**< on input: hash list to search; on exit: hash list entry corresponding
                                              *   to element after retrieved one, or NULL */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   uint64_t              keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   SCIP_MULTIHASHLIST* h;

   assert(multihashlist != NULL);

   /* find hash list entry */
   h = multihashlistFind(*multihashlist, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key);

   /* return element */
   if( h != NULL )
   {
      *multihashlist = h->next;

      return h->element;
   }

   *multihashlist = NULL;

   return NULL;
}

/** removes element from the multihash list */
static
SCIP_Bool multihashlistRemove(
   SCIP_MULTIHASHLIST**  multihashlist,      /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 element             /**< element to remove from the list */
   )
{
   SCIP_MULTIHASHLIST* nextlist;

   assert(multihashlist != NULL);
   assert(blkmem != NULL);
   assert(element != NULL);

   while( *multihashlist != NULL && (*multihashlist)->element != element )
      multihashlist = &(*multihashlist)->next;

   if( *multihashlist != NULL )
   {
      nextlist = (*multihashlist)->next;
      BMSfreeBlockMemory(blkmem, multihashlist);
      *multihashlist = nextlist;

      return TRUE;
   }

   return FALSE;
}

#define SCIP_MULTIHASH_MAXSIZE 33554431 /* 2^25 - 1*/
#define SCIP_MULTIHASH_RESIZE_PERCENTAGE 65
#define SCIP_MULTIHASH_GROW_FACTOR 1.31

/** resizing(increasing) the given multihash */
static
SCIP_RETCODE multihashResize(
   SCIP_MULTIHASH*       multihash           /**< hash table */
   )
{
   SCIP_MULTIHASHLIST** newlists;
   SCIP_MULTIHASHLIST* multihashlist;
   SCIP_Longint nelements;
   int nnewlists;
   int l;

   assert(multihash != NULL);
   assert(multihash->lists != NULL);
   assert(multihash->nlists > 0);
   assert(multihash->hashgetkey != NULL);
   assert(multihash->hashkeyeq != NULL);
   assert(multihash->hashkeyval != NULL);

   /* get new memeory for hash table lists */
   nnewlists = (int) MIN((unsigned int)(multihash->nlists * SCIP_MULTIHASH_GROW_FACTOR), SCIP_MULTIHASH_MAXSIZE);
   nnewlists = MAX(nnewlists, multihash->nlists);

   SCIPdebugMessage("load = %g, nelements = %" SCIP_LONGINT_FORMAT ", nlists = %d, nnewlist = %d\n", SCIPmultihashGetLoad(multihash), multihash->nelements, multihash->nlists, nnewlists);

   if( nnewlists > multihash->nlists )
   {
      SCIP_Bool onlyone;
      void* key;
      uint64_t keyval;
      unsigned int hashval;

      SCIP_ALLOC( BMSallocClearBlockMemoryArray(multihash->blkmem, &newlists, nnewlists) );

      /* move all lists */
      for( l = multihash->nlists - 1; l >= 0; --l )
      {
         multihashlist = multihash->lists[l];
         onlyone = TRUE;

         /* move all elements frmm the old lists into the new lists */
         while( multihashlist != NULL )
         {
            /* get the hash key and its hash value */
            key = multihash->hashgetkey(multihash->userptr, multihashlist->element);
            keyval = multihash->hashkeyval(multihash->userptr, key);
            hashval = keyval % nnewlists; /*lint !e573*/

            /* if the old hash table list consists of only one entry, we still can use this old memory block instead
             * of creating a new one
             */
            if( multihashlist->next == NULL && onlyone )
            {
               /* the new list is also empty, we can directly copy the entry */
               if( newlists[hashval] == NULL )
                  newlists[hashval] = multihashlist;
               /* the new list is not empty, so we need to find the first empty spot */
               else
               {
                  SCIP_MULTIHASHLIST* lastnext = newlists[hashval];
                  SCIP_MULTIHASHLIST* next = lastnext->next;

                  while( next != NULL )
                  {
                     lastnext = next;
                     next = next->next;
                  }

                  lastnext->next = multihashlist;
               }

               multihash->lists[l] = NULL;
            }
            else
            {
               /* append old element to the list at the hash position */
               SCIP_CALL( multihashlistAppend(&(newlists[hashval]), multihash->blkmem, multihashlist->element) );
            }

            onlyone = FALSE;
            multihashlist = multihashlist->next;
         }
      }

      /* remember number of elements */
      nelements = multihash->nelements;
      /* clear old lists */
      SCIPmultihashRemoveAll(multihash);
      /* free old lists */
      BMSfreeBlockMemoryArray(multihash->blkmem, &(multihash->lists), multihash->nlists);

      /* set new data */
      multihash->lists = newlists;
      multihash->nlists = nnewlists;
      multihash->nelements = nelements;

#ifdef SCIP_MORE_DEBUG
      {
         SCIP_Longint sumslotsize = 0;

         for( l = 0; l < multihash->nlists; ++l )
         {
            multihashlist = multihash->lists[l];
            while( multihashlist != NULL )
            {
               sumslotsize++;
               multihashlist = multihashlist->next;
            }
         }
         assert(sumslotsize == multihash->nelements);
      }
#endif
   }

   return SCIP_OKAY;
}

/** creates a multihash table */
SCIP_RETCODE SCIPmultihashCreate(
   SCIP_MULTIHASH**      multihash,          /**< pointer to store the created multihash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store multihash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr             /**< user pointer */
   )
{
   /* only assert non negative to catch overflow errors
    * but not zeros due to integer divison
    */
   assert(tablesize >= 0);
   assert(multihash != NULL);
   assert(hashgetkey != NULL);
   assert(hashkeyeq != NULL);
   assert(hashkeyval != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, multihash) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &(*multihash)->lists, tablesize) );
   (*multihash)->blkmem = blkmem;
   (*multihash)->nlists = tablesize;
   (*multihash)->hashgetkey = hashgetkey;
   (*multihash)->hashkeyeq = hashkeyeq;
   (*multihash)->hashkeyval = hashkeyval;
   (*multihash)->userptr = userptr;
   (*multihash)->nelements = 0;

   return SCIP_OKAY;
}

/** frees the multihash table */
void SCIPmultihashFree(
   SCIP_MULTIHASH**      multihash           /**< pointer to the multihash table */
   )
{
   int i;
   SCIP_MULTIHASH* table;
   BMS_BLKMEM* blkmem;
   SCIP_MULTIHASHLIST** lists;

   assert(multihash != NULL);
   assert(*multihash != NULL);

   table = (*multihash);
   blkmem = table->blkmem;
   lists = table->lists;

   /* free hash lists */
   for( i = table->nlists - 1; i >= 0; --i )
      multihashlistFree(&lists[i], blkmem);

   /* free main hash table data structure */
   BMSfreeBlockMemoryArray(blkmem, &table->lists, table->nlists);
   BMSfreeBlockMemory(blkmem, multihash);
}


/** inserts element in multihash table (multiple inserts of same element possible)
 *
 *  @note A pointer to a multihashlist returned by SCIPmultihashRetrieveNext() might get invalid when adding an element
 *        to the hash table, due to dynamic resizing.
 */
SCIP_RETCODE SCIPmultihashInsert(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 element             /**< element to insert into the table */
   )
{
   void* key;
   uint64_t keyval;
   unsigned int hashval;

   assert(multihash != NULL);
   assert(multihash->lists != NULL);
   assert(multihash->nlists > 0);
   assert(multihash->hashgetkey != NULL);
   assert(multihash->hashkeyeq != NULL);
   assert(multihash->hashkeyval != NULL);
   assert(element != NULL);

   /* dynamically resizing the hashtables */
   if( SCIPmultihashGetLoad(multihash) > SCIP_MULTIHASH_RESIZE_PERCENTAGE )
   {
      SCIP_CALL( multihashResize(multihash) );
   }

   /* get the hash key and its hash value */
   key = multihash->hashgetkey(multihash->userptr, element);
   keyval = multihash->hashkeyval(multihash->userptr, key);
   hashval = keyval % multihash->nlists; /*lint !e573*/

   /* append element to the list at the hash position */
   SCIP_CALL( multihashlistAppend(&multihash->lists[hashval], multihash->blkmem, element) );

   ++(multihash->nelements);

   return SCIP_OKAY;
}

/** inserts element in multihash table (multiple insertion of same element is checked and results in an error)
 *
 *  @note A pointer to a multihashlist returned by SCIPmultihashRetrieveNext() might get invalid when adding a new
 *        element to the multihash table, due to dynamic resizing.
 */
SCIP_RETCODE SCIPmultihashSafeInsert(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 element             /**< element to insert into the table */
   )
{
   assert(multihash != NULL);
   assert(multihash->hashgetkey != NULL);

   /* check, if key is already existing */
   if( SCIPmultihashRetrieve(multihash, multihash->hashgetkey(multihash->userptr, element)) != NULL )
      return SCIP_KEYALREADYEXISTING;

   /* insert element in hash table */
   SCIP_CALL( SCIPmultihashInsert(multihash, element) );

   return SCIP_OKAY;
}

/** retrieve element with key from multihash table, returns NULL if not existing */
void* SCIPmultihashRetrieve(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 key                 /**< key to retrieve */
   )
{
   uint64_t keyval;
   unsigned int hashval;

   assert(multihash != NULL);
   assert(multihash->lists != NULL);
   assert(multihash->nlists > 0);
   assert(multihash->hashgetkey != NULL);
   assert(multihash->hashkeyeq != NULL);
   assert(multihash->hashkeyval != NULL);
   assert(key != NULL);

   /* get the hash value of the key */
   keyval = multihash->hashkeyval(multihash->userptr, key);
   hashval = keyval % multihash->nlists; /*lint !e573*/

   return multihashlistRetrieve(multihash->lists[hashval], multihash->hashgetkey, multihash->hashkeyeq,
      multihash->hashkeyval, multihash->userptr, keyval, key);
}

/** retrieve element with key from multihash table, returns NULL if not existing
 *  can be used to retrieve all entries with the same key (one-by-one)
 *
 *  @note The returned multimultihashlist pointer might get invalid when adding a new element to the multihash table.
 */
void* SCIPmultihashRetrieveNext(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   SCIP_MULTIHASHLIST**  multihashlist,      /**< input: entry in hash table list from which to start searching, or NULL
                                              *   output: entry in hash table list corresponding to element after
                                              *           retrieved one, or NULL */
   void*                 key                 /**< key to retrieve */
   )
{
   uint64_t keyval;

   assert(multihash != NULL);
   assert(multihash->lists != NULL);
   assert(multihash->nlists > 0);
   assert(multihash->hashgetkey != NULL);
   assert(multihash->hashkeyeq != NULL);
   assert(multihash->hashkeyval != NULL);
   assert(multihashlist != NULL);
   assert(key != NULL);

   keyval = multihash->hashkeyval(multihash->userptr, key);

   if( *multihashlist == NULL )
   {
      unsigned int hashval;

      /* get the hash value of the key */
      hashval = keyval % multihash->nlists; /*lint !e573*/

      *multihashlist = multihash->lists[hashval];
   }

   return multihashlistRetrieveNext(multihashlist, multihash->hashgetkey, multihash->hashkeyeq,
      multihash->hashkeyval, multihash->userptr, keyval, key);
}

/** returns whether the given element exists in the multihash table */
SCIP_Bool SCIPmultihashExists(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 element             /**< element to search in the table */
   )
{
   void* key;
   uint64_t keyval;
   unsigned int hashval;

   assert(multihash != NULL);
   assert(multihash->lists != NULL);
   assert(multihash->nlists > 0);
   assert(multihash->hashgetkey != NULL);
   assert(multihash->hashkeyeq != NULL);
   assert(multihash->hashkeyval != NULL);
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = multihash->hashgetkey(multihash->userptr, element);
   keyval = multihash->hashkeyval(multihash->userptr, key);
   hashval = keyval % multihash->nlists; /*lint !e573*/

   return (multihashlistFind(multihash->lists[hashval], multihash->hashgetkey, multihash->hashkeyeq,
         multihash->hashkeyval, multihash->userptr, keyval, key) != NULL);
}

/** removes element from the multihash table, if it exists */
SCIP_RETCODE SCIPmultihashRemove(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 element             /**< element to remove from the table */
   )
{
   void* key;
   uint64_t keyval;
   unsigned int hashval;

   assert(multihash != NULL);
   assert(multihash->lists != NULL);
   assert(multihash->nlists > 0);
   assert(multihash->hashgetkey != NULL);
   assert(multihash->hashkeyeq != NULL);
   assert(multihash->hashkeyval != NULL);
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = multihash->hashgetkey(multihash->userptr, element);
   keyval = multihash->hashkeyval(multihash->userptr, key);
   hashval = keyval % multihash->nlists; /*lint !e573*/

   /* remove element from the list at the hash position */
   if( multihashlistRemove(&multihash->lists[hashval], multihash->blkmem, element) )
      --(multihash->nelements);

   return SCIP_OKAY;
}

/** removes all elements of the multihash table
 *
 *  @note From a performance point of view you should not fill and clear a hash table too often since the clearing can
 *        be expensive. Clearing is done by looping over all buckets and removing the hash table lists one-by-one.
 */
void SCIPmultihashRemoveAll(
   SCIP_MULTIHASH*       multihash           /**< multihash table */
   )
{
   BMS_BLKMEM* blkmem;
   SCIP_MULTIHASHLIST** lists;
   int i;

   assert(multihash != NULL);

   blkmem = multihash->blkmem;
   lists = multihash->lists;

   /* free hash lists */
   for( i = multihash->nlists - 1; i >= 0; --i )
      multihashlistFree(&lists[i], blkmem);

   multihash->nelements = 0;
}

/** returns number of multihash table elements */
SCIP_Longint SCIPmultihashGetNElements(
   SCIP_MULTIHASH*       multihash           /**< multihash table */
   )
{
   assert(multihash != NULL);

   return multihash->nelements;
}

/** returns the load of the given multihash table in percentage */
SCIP_Real SCIPmultihashGetLoad(
   SCIP_MULTIHASH*       multihash           /**< multihash table */
   )
{
   assert(multihash != NULL);

   return ((SCIP_Real)(multihash->nelements) / (multihash->nlists) * 100.0);
}

/** prints statistics about multihash table usage */
void SCIPmultihashPrintStatistics(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_MULTIHASHLIST* multihashlist;
   int usedslots;
   int maxslotsize;
   int sumslotsize;
   int slotsize;
   int i;

   assert(multihash != NULL);

   usedslots = 0;
   maxslotsize = 0;
   sumslotsize = 0;
   for( i = 0; i < multihash->nlists; ++i )
   {
      multihashlist = multihash->lists[i];
      if( multihashlist != NULL )
      {
         usedslots++;
         slotsize = 0;
         while( multihashlist != NULL )
         {
            slotsize++;
            multihashlist = multihashlist->next;
         }
         maxslotsize = MAX(maxslotsize, slotsize);
         sumslotsize += slotsize;
      }
   }
   assert(sumslotsize == multihash->nelements);

   SCIPmessagePrintInfo(messagehdlr, "%" SCIP_LONGINT_FORMAT " multihash entries, used %d/%d slots (%.1f%%)",
      multihash->nelements, usedslots, multihash->nlists, 100.0*(SCIP_Real)usedslots/(SCIP_Real)(multihash->nlists));
   if( usedslots > 0 )
      SCIPmessagePrintInfo(messagehdlr, ", avg. %.1f entries/used slot, max. %d entries in slot",
         (SCIP_Real)(multihash->nelements)/(SCIP_Real)usedslots, maxslotsize);
   SCIPmessagePrintInfo(messagehdlr, "\n");
}

/** creates a hash table */
SCIP_RETCODE SCIPhashtableCreate(
   SCIP_HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr             /**< user pointer */
   )
{
   unsigned int nslots;

   /* only assert non negative to catch overflow errors
    * but not zeros due to integer divison
    */
   assert(tablesize >= 0);
   assert(hashtable != NULL);
   assert(hashgetkey != NULL);
   assert(hashkeyeq != NULL);
   assert(hashkeyval != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, hashtable) );

   /* dont create too small hashtables, i.e. at least size 32, and increase
    * the given size by divinding it by 0.9, since then no rebuilding will
    * be necessary if the given number of elements are inserted. Finally round
    * to the next power of two.
    */
   (*hashtable)->shift = 32;
   (*hashtable)->shift -= (int)ceil(LOG2(MAX(32.0, tablesize / 0.9)));

   /* compute size from shift */
   nslots = 1u << (32 - (*hashtable)->shift);

   /* compute mask to do a fast modulo by nslots using bitwise and */
   (*hashtable)->mask = nslots - 1;
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*hashtable)->slots, nslots) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &(*hashtable)->hashes, nslots) );
   (*hashtable)->blkmem = blkmem;
   (*hashtable)->hashgetkey = hashgetkey;
   (*hashtable)->hashkeyeq = hashkeyeq;
   (*hashtable)->hashkeyval = hashkeyval;
   (*hashtable)->userptr = userptr;
   (*hashtable)->nelements = 0;

   return SCIP_OKAY;
}

/** frees the hash table */
void SCIPhashtableFree(
   SCIP_HASHTABLE**      hashtable           /**< pointer to the hash table */
   )
{
   uint32_t nslots;
   SCIP_HASHTABLE* table;

   assert(hashtable != NULL);
   assert(*hashtable != NULL);
   table = *hashtable;
   nslots = (*hashtable)->mask + 1;
#ifdef SCIP_DEBUG
   {
      uint32_t maxprobelen = 0;
      uint64_t probelensum = 0;
      uint32_t i;

      assert(table != NULL);

      for( i = 0; i < nslots; ++i )
      {
         if( table->hashes[i] != 0 )
         {
            uint32_t probelen = ((i + table->mask + 1 - (table->hashes[i]>>(table->shift))) & table->mask) + 1;
            probelensum += probelen;
            maxprobelen = MAX(probelen, maxprobelen);
         }
      }

      SCIPdebugMessage("%u hash table entries, used %u/%u slots (%.1f%%)",
                       (unsigned int)table->nelements, (unsigned int)table->nelements, (unsigned int)nslots,
                       100.0*(SCIP_Real)table->nelements/(SCIP_Real)(nslots));
      if( table->nelements > 0 )
         SCIPdebugMessage(", avg. probe length is %.1f, max. probe length is %u",
                              (SCIP_Real)(probelensum)/(SCIP_Real)table->nelements, (unsigned int)maxprobelen);
      SCIPdebugMessage("\n");
   }
#endif

   /* free main hash table data structure */
   BMSfreeBlockMemoryArray((*hashtable)->blkmem, &table->hashes, nslots);
   BMSfreeBlockMemoryArray((*hashtable)->blkmem, &table->slots, nslots);
   BMSfreeBlockMemory((*hashtable)->blkmem, hashtable);
}

/** removes all elements of the hash table
 *
 *  @note From a performance point of view you should not fill and clear a hash table too often since the clearing can
 *        be expensive. Clearing is done by looping over all buckets and removing the hash table lists one-by-one.
 *
 *  @deprecated Please use SCIPhashtableRemoveAll()
 */
void SCIPhashtableClear(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   )
{
   SCIPhashtableRemoveAll(hashtable);
}

/* computes the distance from it's desired position for the element stored at pos */
#define ELEM_DISTANCE(pos) (((pos) + hashtable->mask + 1 - (hashtable->hashes[(pos)]>>(hashtable->shift))) & hashtable->mask)

/** inserts element in hash table (multiple inserts of same element overrides previous one) */
static
SCIP_RETCODE hashtableInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element,            /**< element to insert into the table */
   void*                 key,                /**< key of element */
   uint32_t              hashval,            /**< hash value of element */
   SCIP_Bool             override            /**< should element be overridden or an error be returned if already existing */
   )
{
   uint32_t elemdistance;
   uint32_t pos;
#ifndef NDEBUG
   SCIP_Bool swapped = FALSE;
#endif

   assert(hashtable != NULL);
   assert(hashtable->slots != NULL);
   assert(hashtable->hashes != NULL);
   assert(hashtable->mask > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(element != NULL);

   pos = hashval>>(hashtable->shift);
   elemdistance = 0;
   while( TRUE ) /*lint !e716*/
   {
      uint32_t distance;

      /* if position is empty or key equal insert element */
      if( hashtable->hashes[pos] == 0 )
      {
         hashtable->slots[pos] = element;
         hashtable->hashes[pos] = hashval;
         ++hashtable->nelements;
         return SCIP_OKAY;
      }

      if( hashtable->hashes[pos] == hashval && hashtable->hashkeyeq(hashtable->userptr,
             hashtable->hashgetkey(hashtable->userptr, hashtable->slots[pos]), key) )
      {
         if( override )
         {
#ifndef NDEBUG
            assert(! swapped);
#endif
            hashtable->slots[pos] = element;
            hashtable->hashes[pos] = hashval;
            return SCIP_OKAY;
         }
         else
         {
            return SCIP_KEYALREADYEXISTING;
         }
      }

      /* otherwise check if the current element at this position is closer to its hashvalue */
      distance = ELEM_DISTANCE(pos);
      if( distance < elemdistance )
      {
         uint32_t tmp;

         /* if this is the case we insert the new element here and find a new position for the old one */
         elemdistance = distance;
         SCIPswapPointers(&hashtable->slots[pos], &element);
         tmp = hashval;
         hashval = hashtable->hashes[pos];
         hashtable->hashes[pos] = tmp;
         key = hashtable->hashgetkey(hashtable->userptr, element);

         /* after doing a swap the case that other elements are replaced must not happen anymore */
#ifndef NDEBUG
         swapped = TRUE;
#endif
      }

      /* continue until we have found an empty position */
      pos = (pos + 1) & hashtable->mask;
      ++elemdistance;
   }
}

/** check if the load factor of the hashtable is too high and rebuild if necessary */
static
SCIP_RETCODE hashtableCheckLoad(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   )
{
   assert(hashtable != NULL);
   assert(hashtable->shift < 32);

   /* use integer arithmetic to approximately check if load factor is above 90% */
   if( ((((uint64_t)hashtable->nelements)<<10)>>(32-hashtable->shift) > 921) )
   {
      void** slots;
      uint32_t* hashes;
      uint32_t nslots;
      uint32_t newnslots;
      uint32_t i;

      /* calculate new size (always power of two) */
      nslots = hashtable->mask + 1;
      newnslots = 2*nslots;
      hashtable->mask = newnslots-1;
      --hashtable->shift;

      /* reallocate array */
      SCIP_ALLOC( BMSallocBlockMemoryArray(hashtable->blkmem, &slots, newnslots) );
      SCIP_ALLOC( BMSallocClearBlockMemoryArray(hashtable->blkmem, &hashes, newnslots) );

      SCIPswapPointers((void**) &slots, (void**) &hashtable->slots);
      SCIPswapPointers((void**) &hashes, (void**) &hashtable->hashes);
      hashtable->nelements = 0;

      /* reinsert all elements */
      for( i = 0; i < nslots; ++i )
      {
         /* using SCIP_CALL_ABORT since there are no allocations or duplicates
          * and thus no bad return codes when inserting the elements
          */
         if( hashes[i] != 0 )
         {
            SCIP_CALL_ABORT( hashtableInsert(hashtable, slots[i], hashtable->hashgetkey(hashtable->userptr, slots[i]), hashes[i], FALSE) );
         }
      }
      BMSfreeBlockMemoryArray(hashtable->blkmem, &hashes, nslots);
      BMSfreeBlockMemoryArray(hashtable->blkmem, &slots, nslots);
   }

   return SCIP_OKAY;
}


/** inserts element in hash table
 *
 *  @note multiple inserts of same element overrides previous one
 */
SCIP_RETCODE SCIPhashtableInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   )
{
   void* key;
   uint64_t keyval;
   uint32_t hashval;

   assert(hashtable != NULL);
   assert(hashtable->slots != NULL);
   assert(hashtable->hashes != NULL);
   assert(hashtable->mask > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(element != NULL);

   SCIP_CALL( hashtableCheckLoad(hashtable) );

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = hashvalue(keyval);

   return hashtableInsert(hashtable, element, key, hashval, TRUE);
}

/** inserts element in hash table
 *
 *  @note multiple insertion of same element is checked and results in an error
 */
SCIP_RETCODE SCIPhashtableSafeInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   )
{
   void* key;
   uint64_t keyval;
   uint32_t hashval;

   assert(hashtable != NULL);
   assert(hashtable->slots != NULL);
   assert(hashtable->hashes != NULL);
   assert(hashtable->mask > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(element != NULL);

   SCIP_CALL( hashtableCheckLoad(hashtable) );

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = hashvalue(keyval);

   return hashtableInsert(hashtable, element, key, hashval, FALSE);
}

/** retrieve element with key from hash table, returns NULL if not existing */
void* SCIPhashtableRetrieve(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 key                 /**< key to retrieve */
   )
{
   uint64_t keyval;
   uint32_t hashval;
   uint32_t pos;
   uint32_t elemdistance;

   assert(hashtable != NULL);
   assert(hashtable->slots != NULL);
   assert(hashtable->hashes != NULL);
   assert(hashtable->mask > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(key != NULL);

   /* get the hash value of the key */
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = hashvalue(keyval);

   pos = hashval>>(hashtable->shift);
   elemdistance = 0;

   while( TRUE ) /*lint !e716*/
   {
      uint32_t distance;

      /* slots is empty so element cannot be contained */
      if( hashtable->hashes[pos] == 0 )
         return NULL;

      distance = ELEM_DISTANCE(pos);

      /* element cannot be contained since otherwise we would have swapped it with this one during insert */
      if( elemdistance > distance )
         return NULL;

      /* found element */
      if( hashtable->hashes[pos] == hashval && hashtable->hashkeyeq(hashtable->userptr,
             hashtable->hashgetkey(hashtable->userptr, hashtable->slots[pos]), key) )
         return hashtable->slots[pos];

      pos = (pos + 1) & hashtable->mask;
      ++elemdistance;
   }
}

/** returns whether the given element exists in the table */
SCIP_Bool SCIPhashtableExists(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to search in the table */
   )
{
   assert(hashtable != NULL);
   assert(hashtable->slots != NULL);
   assert(hashtable->hashes != NULL);
   assert(hashtable->mask > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(element != NULL);

   return (SCIPhashtableRetrieve(hashtable, hashtable->hashgetkey(hashtable->userptr, element)) != NULL);
}

/** removes element from the hash table, if it exists */
SCIP_RETCODE SCIPhashtableRemove(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to remove from the table */
   )
{
   void* key;
   uint64_t keyval;
   uint32_t hashval;
   uint32_t elemdistance;
   uint32_t distance;
   uint32_t pos;

   assert(hashtable != NULL);
   assert(hashtable->slots != NULL);
   assert(hashtable->hashes != NULL);
   assert(hashtable->mask > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = hashvalue(keyval);

   elemdistance = 0;
   pos = hashval>>(hashtable->shift);
   while( TRUE ) /*lint !e716*/
   {
      /* slots empty so element not contained */
      if( hashtable->hashes[pos] == 0 )
         return SCIP_OKAY;

      distance = ELEM_DISTANCE(pos);

      /* element can not be contained since otherwise we would have swapped it with this one */
      if( elemdistance > distance )
         return SCIP_OKAY;

      if( hashtable->hashes[pos] == hashval && hashtable->hashkeyeq(hashtable->userptr,
             hashtable->hashgetkey(hashtable->userptr, hashtable->slots[pos]), key) )
      {
         /* element exists at pos so break out of loop */
         break;
      }

      pos = (pos + 1) & hashtable->mask;
      ++elemdistance;
   }

   /* remove element */
   hashtable->hashes[pos] = 0;
   --hashtable->nelements;
   while( TRUE ) /*lint !e716*/
   {
      uint32_t nextpos = (pos + 1) & hashtable->mask;

      /* nothing to do since there is no chain that needs to be moved */
      if( hashtable->hashes[nextpos] == 0 )
         break;

      /* check if the element is the start of a new chain and return if that is the case */
      if( (hashtable->hashes[nextpos]>>(hashtable->shift)) == nextpos )
         break;

      /* element should be moved to the left and next element needs to be checked */
      hashtable->slots[pos] = hashtable->slots[nextpos];
      hashtable->hashes[pos] = hashtable->hashes[nextpos];
      hashtable->hashes[nextpos] = 0;

      pos = nextpos;
   }

   return SCIP_OKAY;
}

/** removes all elements of the hash table */
void SCIPhashtableRemoveAll(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   )
{
   assert(hashtable != NULL);

   BMSclearMemoryArray(hashtable->hashes, hashtable->mask + 1);

   hashtable->nelements = 0;
}

/** returns number of hash table elements */
SCIP_Longint SCIPhashtableGetNElements(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   )
{
   assert(hashtable != NULL);

   return hashtable->nelements;
}

/** gives the number of entries in the internal arrays of a hash table */
int SCIPhashtableGetNEntries(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   )
{
   return (int) hashtable->mask + 1;
}

/** gives the element at the given index or NULL if entry at that index has no element */
void* SCIPhashtableGetEntry(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   int                   entryidx            /**< index of hash table entry */
   )
{
   return hashtable->hashes[entryidx] == 0 ? NULL : hashtable->slots[entryidx];
}

/** returns the load of the given hash table in percentage */
SCIP_Real SCIPhashtableGetLoad(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   )
{
   assert(hashtable != NULL);

   return ((SCIP_Real)(hashtable->nelements) / (hashtable->mask + 1) * 100.0);
}

/** prints statistics about hash table usage */
void SCIPhashtablePrintStatistics(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   uint32_t maxprobelen = 0;
   uint64_t probelensum = 0;
   uint32_t nslots;
   uint32_t i;

   assert(hashtable != NULL);

   nslots = hashtable->mask + 1;

   /* compute the maximum and average probe length */
   for( i = 0; i < nslots; ++i )
   {
      if( hashtable->hashes[i] != 0 )
      {
         uint32_t probelen = ELEM_DISTANCE(i) + 1;
         probelensum += probelen;
         maxprobelen = MAX(probelen, maxprobelen);
      }
   }

   /* print general hash table statistics */
   SCIPmessagePrintInfo(messagehdlr, "%u hash entries, used %u/%u slots (%.1f%%)",
                        (unsigned int)hashtable->nelements, (unsigned int)hashtable->nelements,
                        (unsigned int)nslots, 100.0*(SCIP_Real)hashtable->nelements/(SCIP_Real)(nslots));

   /* if not empty print average and maximum probe length */
   if( hashtable->nelements > 0 )
      SCIPmessagePrintInfo(messagehdlr, ", avg. probe length is %.1f, max. probe length is %u",
         (SCIP_Real)(probelensum)/(SCIP_Real)hashtable->nelements, (unsigned int)maxprobelen);
   SCIPmessagePrintInfo(messagehdlr, "\n");
}

/** returns TRUE iff both keys (i.e. strings) are equal */
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqString)
{  /*lint --e{715}*/
   const char* string1 = (const char*)key1;
   const char* string2 = (const char*)key2;

   return (strcmp(string1, string2) == 0);
}

/** returns the hash value of the key (i.e. string) */
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValString)
{  /*lint --e{715}*/
   const char* str;
   uint64_t hash;

   str = (const char*)key;
   hash = 37;
   while( *str != '\0' )
   {
      hash *= 11;
      hash += (unsigned int)(*str); /*lint !e571*/
      str++;
   }

   return hash;
}


/** gets the element as the key */
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyStandard)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys(pointer) are equal */
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqPtr)
{  /*lint --e{715}*/
   return (key1 == key2);
}

/** returns the hash value of the key */
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValPtr)
{  /*lint --e{715}*/
   /* the key is used as the keyvalue too */
   return (uint64_t) (uintptr_t) key;
}



/*
 * Hash Map
 */

/* redefine ELEM_DISTANCE macro for hashmap */
#undef ELEM_DISTANCE
/* computes the distance from it's desired position for the element stored at pos */
#define ELEM_DISTANCE(pos) (((pos) + hashmap->mask + 1 - (hashmap->hashes[(pos)]>>(hashmap->shift))) & hashmap->mask)


/** inserts element in hash table */
static
SCIP_RETCODE hashmapInsert(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< element to insert into the table */
   SCIP_HASHMAPIMAGE     image,              /**< key of element */
   uint32_t              hashval,            /**< hash value of element */
   SCIP_Bool             override            /**< should element be overridden or error be returned if already existing */
   )
{
   uint32_t elemdistance;
   uint32_t pos;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->hashes != NULL);
   assert(hashmap->mask > 0);
   assert(hashval != 0);

   pos = hashval>>(hashmap->shift);
   elemdistance = 0;
   while( TRUE ) /*lint !e716*/
   {
      uint32_t distance;

      /* if position is empty or key equal insert element */
      if( hashmap->hashes[pos] == 0 )
      {
         hashmap->slots[pos].origin = origin;
         hashmap->slots[pos].image = image;
         hashmap->hashes[pos] = hashval;
         ++hashmap->nelements;
         return SCIP_OKAY;
      }

      if( hashval == hashmap->hashes[pos] && origin == hashmap->slots[pos].origin )
      {
         if( override )
         {
            hashmap->slots[pos].origin = origin;
            hashmap->slots[pos].image = image;
            hashmap->hashes[pos] = hashval;
            return SCIP_OKAY;
         }
         else
         {
            return SCIP_KEYALREADYEXISTING;
         }
      }

      /* otherwise check if the current element at this position is closer to its hashvalue */
      distance = ELEM_DISTANCE(pos);
      if( distance < elemdistance )
      {
         SCIP_HASHMAPIMAGE tmp;
         uint32_t tmphash;

         /* if this is the case we insert the new element here and find a new position for the old one */
         elemdistance = distance;
         tmphash = hashval;
         hashval = hashmap->hashes[pos];
         hashmap->hashes[pos] = tmphash;
         SCIPswapPointers(&hashmap->slots[pos].origin, &origin);
         tmp = image;
         image = hashmap->slots[pos].image;
         hashmap->slots[pos].image = tmp;
      }

      /* continue until we have found an empty position */
      pos = (pos + 1) & hashmap->mask;
      ++elemdistance;
   }
}

/** lookup origin in the hashmap. If element is found returns true and the position of the element,
 *  otherwise returns FALSE.
 */
static
SCIP_Bool hashmapLookup(
   SCIP_HASHMAP*         hashmap,            /**< hash table */
   void*                 origin,             /**< origin to lookup */
   uint32_t*             pos                 /**< pointer to store position of element, if exists */
   )
{
   uint32_t hashval;
   uint32_t elemdistance;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->hashes != NULL);
   assert(hashmap->mask > 0);

   /* get the hash value */
   hashval = hashvalue((size_t)origin);
   assert(hashval != 0);

   *pos = hashval>>(hashmap->shift);
   elemdistance = 0;

   while( TRUE ) /*lint !e716*/
   {
      uint32_t distance;

      /* slots is empty so element cannot be contained */
      if( hashmap->hashes[*pos] == 0 )
         return FALSE;

      distance = ELEM_DISTANCE(*pos);
      /* element can not be contained since otherwise we would have swapped it with this one during insert */
      if( elemdistance > distance )
         return FALSE;

      /* found element */
      if( hashmap->hashes[*pos] == hashval && hashmap->slots[*pos].origin == origin )
         return TRUE;

      *pos = (*pos + 1) & hashmap->mask;
      ++elemdistance;
   }
}

/** check if the load factor of the hashmap is too high and rebuild if necessary */
static
SCIP_RETCODE hashmapCheckLoad(
   SCIP_HASHMAP*         hashmap             /**< hash table */
   )
{
   assert(hashmap != NULL);
   assert(hashmap->shift < 32);

   /* use integer arithmetic to approximately check if load factor is above 90% */
   if( ((((uint64_t)hashmap->nelements)<<10)>>(32-hashmap->shift) > 921) )
   {
      SCIP_HASHMAPENTRY* slots;
      uint32_t* hashes;
      uint32_t nslots;
      uint32_t newnslots;
      uint32_t i;

      /* calculate new size (always power of two) */
      nslots = hashmap->mask + 1;
      --hashmap->shift;
      newnslots = 2*nslots;
      hashmap->mask = newnslots-1;

      /* reallocate array */
      SCIP_ALLOC( BMSallocBlockMemoryArray(hashmap->blkmem, &slots, newnslots) );
      SCIP_ALLOC( BMSallocClearBlockMemoryArray(hashmap->blkmem, &hashes, newnslots) );

      SCIPswapPointers((void**) &slots, (void**) &hashmap->slots);
      SCIPswapPointers((void**) &hashes, (void**) &hashmap->hashes);
      hashmap->nelements = 0;

      /* reinsert all elements */
      for( i = 0; i < nslots; ++i )
      {
         /* using SCIP_CALL_ABORT since there are no allocations or duplicates
          * and thus no bad return codes when inserting the elements
          */
         if( hashes[i] != 0 )
         {
            SCIP_CALL_ABORT( hashmapInsert(hashmap, slots[i].origin, slots[i].image, hashes[i], FALSE) );
         }
      }

      /* free old arrays */
      BMSfreeBlockMemoryArray(hashmap->blkmem, &hashes, nslots);
      BMSfreeBlockMemoryArray(hashmap->blkmem, &slots, nslots);
   }

   return SCIP_OKAY;
}

/** creates a hash map mapping pointers to pointers */
SCIP_RETCODE SCIPhashmapCreate(
   SCIP_HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash map entries */
   int                   mapsize             /**< size of the hash map */
   )
{
   uint32_t nslots;

   assert(hashmap != NULL);
   assert(mapsize >= 0);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, hashmap) );

   /* dont create too small hashtables, i.e. at least size 32, and increase
    * the given size by divinding it by 0.9, since then no rebuilding will
    * be necessary if the given number of elements are inserted. Finally round
    * to the next power of two.
    */
   (*hashmap)->shift = 32;
   (*hashmap)->shift -= (int)ceil(log(MAX(32, mapsize / 0.9)) / log(2.0));
   nslots = 1u << (32 - (*hashmap)->shift);
   (*hashmap)->mask = nslots - 1;
   (*hashmap)->blkmem = blkmem;
   (*hashmap)->nelements = 0;

   SCIP_ALLOC( BMSallocBlockMemoryArray((*hashmap)->blkmem, &(*hashmap)->slots, nslots) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray((*hashmap)->blkmem, &(*hashmap)->hashes, nslots) );

   return SCIP_OKAY;
}

/** frees the hash map */
void SCIPhashmapFree(
   SCIP_HASHMAP**        hashmap             /**< pointer to the hash map */
   )
{
   uint32_t nslots;

   assert(hashmap != NULL);
   assert(*hashmap != NULL);

   nslots = (*hashmap)->mask + 1;
#ifdef SCIP_DEBUG
   {
      uint32_t maxprobelen = 0;
      uint64_t probelensum = 0;
      uint32_t i;

      assert(hashmap != NULL);

      for( i = 0; i < nslots; ++i )
      {
         if( (*hashmap)->hashes[i] != 0 )
         {
            uint32_t probelen = ((i + (*hashmap)->mask + 1 - ((*hashmap)->hashes[i]>>((*hashmap)->shift))) & (*hashmap)->mask) + 1;
            probelensum += probelen;
            maxprobelen = MAX(probelen, maxprobelen);
         }
      }

      SCIPdebugMessage("%u hash map entries, used %u/%u slots (%.1f%%)",
                       (unsigned int)(*hashmap)->nelements, (unsigned int)(*hashmap)->nelements, (unsigned int)nslots,
                       100.0*(SCIP_Real)(*hashmap)->nelements/(SCIP_Real)(nslots));
      if( (*hashmap)->nelements > 0 )
         SCIPdebugMessage(", avg. probe length is %.1f, max. probe length is %u",
                          (SCIP_Real)(probelensum)/(SCIP_Real)(*hashmap)->nelements, (unsigned int)maxprobelen);
      SCIPdebugMessage("\n");
   }
#endif

   /* free main hash map data structure */
   BMSfreeBlockMemoryArray((*hashmap)->blkmem, &(*hashmap)->hashes, nslots);
   BMSfreeBlockMemoryArray((*hashmap)->blkmem, &(*hashmap)->slots, nslots);
   BMSfreeBlockMemory((*hashmap)->blkmem, hashmap);
}

/** inserts new origin->image pair in hash map
 *
 *  @note multiple insertion of same element is checked and results in an error
 */
SCIP_RETCODE SCIPhashmapInsert(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   )
{
   uint32_t hashval;
   SCIP_HASHMAPIMAGE img;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->hashes != NULL);
   assert(hashmap->mask > 0);

   SCIP_CALL( hashmapCheckLoad(hashmap) );

   /* get the hash value */
   hashval = hashvalue((size_t)origin);

   /* append origin->image pair to hash map */
   img.ptr = image;
   SCIP_CALL( hashmapInsert(hashmap, origin, img, hashval, FALSE) );

   return SCIP_OKAY;
}

/** inserts new origin->image pair in hash map
 *
 *  @note multiple insertion of same element is checked and results in an error
 */
SCIP_RETCODE SCIPhashmapInsertReal(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   SCIP_Real             image               /**< new image for origin */
   )
{
   uint32_t hashval;
   SCIP_HASHMAPIMAGE img;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->hashes != NULL);
   assert(hashmap->mask > 0);

   SCIP_CALL( hashmapCheckLoad(hashmap) );

   /* get the hash value */
   hashval = hashvalue((size_t)origin);

   /* append origin->image pair to hash map */
   img.real = image;
   SCIP_CALL( hashmapInsert(hashmap, origin, img, hashval, FALSE) );

   return SCIP_OKAY;
}

/** retrieves image of given origin from the hash map, or NULL if no image exists */
void* SCIPhashmapGetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   )
{
   uint32_t pos;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->hashes != NULL);
   assert(hashmap->mask > 0);

   if( hashmapLookup(hashmap, origin, &pos) )
      return hashmap->slots[pos].image.ptr;

   return NULL;
}

/** retrieves image of given origin from the hash map, or NULL if no image exists */
SCIP_Real SCIPhashmapGetImageReal(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   )
{
   uint32_t pos;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->hashes != NULL);
   assert(hashmap->mask > 0);

   if( hashmapLookup(hashmap, origin, &pos) )
      return hashmap->slots[pos].image.real;

   return SCIP_INVALID;
}

/** sets image for given origin in the hash map, either by modifying existing origin->image pair
 *  or by appending a new origin->image pair
 */
SCIP_RETCODE SCIPhashmapSetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   )
{
   uint32_t hashval;
   SCIP_HASHMAPIMAGE img;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->mask > 0);

   SCIP_CALL( hashmapCheckLoad(hashmap) );

   /* get the hash value */
   hashval = hashvalue((size_t)origin);

   /* append origin->image pair to hash map */
   img.ptr = image;
   SCIP_CALL( hashmapInsert(hashmap, origin, img, hashval, TRUE) );

   return SCIP_OKAY;
}

/** sets image for given origin in the hash map, either by modifying existing origin->image pair
 *  or by appending a new origin->image pair
 */
SCIP_RETCODE SCIPhashmapSetImageReal(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   SCIP_Real             image               /**< new image for origin */
   )
{
   uint32_t hashval;
   SCIP_HASHMAPIMAGE img;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->mask > 0);

   SCIP_CALL( hashmapCheckLoad(hashmap) );

   /* get the hash value */
   hashval = hashvalue((size_t)origin);

   /* append origin->image pair to hash map */
   img.real = image;
   SCIP_CALL( hashmapInsert(hashmap, origin, img, hashval, TRUE) );

   return SCIP_OKAY;
}

/** checks whether an image to the given origin exists in the hash map */
SCIP_Bool SCIPhashmapExists(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to search for */
   )
{
   uint32_t pos;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->hashes != NULL);
   assert(hashmap->mask > 0);

   return hashmapLookup(hashmap, origin, &pos);
}

/** removes origin->image pair from the hash map, if it exists */
SCIP_RETCODE SCIPhashmapRemove(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to remove from the list */
   )
{
   uint32_t pos;

   assert(hashmap != NULL);
   assert(hashmap->slots != NULL);
   assert(hashmap->mask > 0);

   assert(origin != NULL);

   if( hashmapLookup(hashmap, origin, &pos) )
   {
      /* remove element */
      hashmap->hashes[pos] = 0;
      --hashmap->nelements;

      /* move other elements if necessary */
      while( TRUE ) /*lint !e716*/
      {
         uint32_t nextpos = (pos + 1) & hashmap->mask;

         /* nothing to do since there is no chain that needs to be moved */
         if( hashmap->hashes[nextpos] == 0 )
            return SCIP_OKAY;

         /* check if the element is the start of a new chain and return if that is the case */
         if( (hashmap->hashes[nextpos]>>(hashmap->shift)) == nextpos )
            return SCIP_OKAY;

         /* element should be moved to the left and next element needs to be checked */
         hashmap->slots[pos].origin = hashmap->slots[nextpos].origin;
         hashmap->slots[pos].image = hashmap->slots[nextpos].image;
         hashmap->hashes[pos] = hashmap->hashes[nextpos];
         hashmap->hashes[nextpos] = 0;

         pos = nextpos;
      }
   }

   return SCIP_OKAY;
}

/** prints statistics about hash map usage */
void SCIPhashmapPrintStatistics(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   uint32_t maxprobelen = 0;
   uint64_t probelensum = 0;
   uint32_t nslots;
   uint32_t i;

   assert(hashmap != NULL);

   nslots = hashmap->mask + 1;

   /* compute the maximum and average probe length */
   for( i = 0; i < nslots; ++i )
   {
      if( hashmap->hashes[i] != 0 )
      {
         uint32_t probelen = ELEM_DISTANCE(i) + 1;
         probelensum += probelen;
         maxprobelen = MAX(probelen, maxprobelen);
      }
   }

   /* print general hash map statistics */
   SCIPmessagePrintInfo(messagehdlr, "%u hash entries, used %u/%u slots (%.1f%%)",
                        (unsigned int)hashmap->nelements, (unsigned int)hashmap->nelements,
                        (unsigned int)nslots, 100.0*(SCIP_Real)hashmap->nelements/(SCIP_Real)(nslots));

   /* if not empty print average and maximum probe length */
   if( hashmap->nelements > 0 )
      SCIPmessagePrintInfo(messagehdlr, ", avg. probe length is %.1f, max. probe length is %u",
         (SCIP_Real)(probelensum)/(SCIP_Real)hashmap->nelements, (unsigned int)maxprobelen);
   SCIPmessagePrintInfo(messagehdlr, "\n");
}

/** indicates whether a hash map has no entries */
SCIP_Bool SCIPhashmapIsEmpty(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   )
{
   assert(hashmap != NULL);

   return hashmap->nelements == 0;
}

/** gives the number of elements in a hash map */
int SCIPhashmapGetNElements(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   )
{
   return (int) hashmap->nelements;
}

/** gives the number of entries in the internal arrays of a hash map */
int SCIPhashmapGetNEntries(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   )
{
   return (int) hashmap->mask + 1;
}

/** gives the hashmap entry at the given index or NULL if entry is empty */
SCIP_HASHMAPENTRY* SCIPhashmapGetEntry(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   int                   entryidx            /**< index of hash map entry */
   )
{
   assert(hashmap != NULL);

   return hashmap->hashes[entryidx] == 0 ? NULL : &hashmap->slots[entryidx];
}

/** gives the origin of the hashmap entry */
void* SCIPhashmapEntryGetOrigin(
   SCIP_HASHMAPENTRY*    entry               /**< hash map entry */
   )
{
   assert(entry != NULL);

   return entry->origin;
}

/** gives the image of the hashmap entry */
void* SCIPhashmapEntryGetImage(
   SCIP_HASHMAPENTRY*    entry               /**< hash map entry */
   )
{
   assert(entry != NULL);

   return entry->image.ptr;
}

/** gives the image of the hashmap entry */
SCIP_Real SCIPhashmapEntryGetImageReal(
   SCIP_HASHMAPENTRY*    entry               /**< hash map entry */
   )
{
   assert(entry != NULL);

   return entry->image.real;
}

/** sets pointer image of a hashmap entry */
void SCIPhashmapEntrySetImage(
   SCIP_HASHMAPENTRY*    entry,              /**< hash map entry */
   void*                 image               /**< new image */
   )
{
   assert(entry != NULL);

   entry->image.ptr = image;
}

/** sets real image of a hashmap entry */
void SCIPhashmapEntrySetImageReal(
   SCIP_HASHMAPENTRY*    entry,              /**< hash map entry */
   SCIP_Real             image               /**< new image */
   )
{
   assert(entry != NULL);

   entry->image.real = image;
}

/** removes all entries in a hash map. */
SCIP_RETCODE SCIPhashmapRemoveAll(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   )
{
   assert(hashmap != NULL);

   BMSclearMemoryArray(hashmap->hashes, hashmap->mask + 1);

   hashmap->nelements = 0;

   return SCIP_OKAY;
}


/*
 * Hash Set
 */

/* redefine ELEM_DISTANCE macro for hashset */
#undef ELEM_DISTANCE
/* computes the distance from it's desired position for the element stored at pos */
#define ELEM_DISTANCE(pos) (((pos) + nslots - hashSetDesiredPos(hashset, hashset->slots[(pos)])) & mask)

/* calculate desired position of element in hash set */
static
uint32_t hashSetDesiredPos(
   SCIP_HASHSET*         hashset,            /**< the hash set */
   void*                 element             /**< element to calculate position for */
   )
{
   return (uint32_t)((UINT64_C(0x9e3779b97f4a7c15) * (uintptr_t)element)>>(hashset->shift));
}

static
void hashsetInsert(
   SCIP_HASHSET*         hashset,            /**< hash set */
   void*                 element             /**< element to insert */
   )
{
   uint32_t elemdistance;
   uint32_t pos;
   uint32_t nslots;
   uint32_t mask;

   assert(hashset != NULL);
   assert(hashset->slots != NULL);
   assert(element != NULL);

   pos = hashSetDesiredPos(hashset, element);
   nslots = (uint32_t)SCIPhashsetGetNSlots(hashset);
   mask = nslots - 1;

   elemdistance = 0;
   while( TRUE ) /*lint !e716*/
   {
      uint32_t distance;

      /* if position is empty or key equal insert element */
      if( hashset->slots[pos] == NULL )
      {
         hashset->slots[pos] = element;
         ++hashset->nelements;
         return;
      }

      if( hashset->slots[pos] == element )
         return;

      /* otherwise check if the current element at this position is closer to its hashvalue */
      distance = ELEM_DISTANCE(pos);
      if( distance < elemdistance )
      {
         /* if this is the case we insert the new element here and find a new position for the old one */
         elemdistance = distance;
         SCIPswapPointers(&hashset->slots[pos], &element);
      }

      /* continue until we have found an empty position */
      pos = (pos + 1) & mask;
      ++elemdistance;
   }
}

/** check if the load factor of the hash set is too high and rebuild if necessary */
static
SCIP_RETCODE hashsetCheckLoad(
   SCIP_HASHSET*         hashset,            /**< hash set */
   BMS_BLKMEM*           blkmem              /**< block memory used to store hash set entries */
   )
{
   assert(hashset != NULL);
   assert(hashset->shift < 64);

   /* use integer arithmetic to approximately check if load factor is above 90% */
   if( ((((uint64_t)hashset->nelements)<<10)>>(64-hashset->shift) > 921) )
   {
      void** slots;
      uint32_t nslots;
      uint32_t newnslots;
      uint32_t i;

      /* calculate new size (always power of two) */
      nslots = (uint32_t)SCIPhashsetGetNSlots(hashset);
      newnslots = 2*nslots;
      --hashset->shift;

      /* reallocate array */
      SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &slots, newnslots) );

      SCIPswapPointers((void**) &slots, (void**) &hashset->slots);
      hashset->nelements = 0;

      /* reinsert all elements */
      for( i = 0; i < nslots; ++i )
      {
         if( slots[i] != NULL )
            hashsetInsert(hashset, slots[i]);
      }

      BMSfreeBlockMemoryArray(blkmem, &slots, nslots);
   }

   return SCIP_OKAY;
}

/** creates a hash set of pointers */
SCIP_RETCODE SCIPhashsetCreate(
   SCIP_HASHSET**        hashset,            /**< pointer to store the created hash set */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash set entries */
   int                   size                /**< initial size of the hash set; it is guaranteed that the set is not
                                              *   resized if at most that many elements are inserted */
   )
{
   uint32_t nslots;

   assert(hashset != NULL);
   assert(size >= 0);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, hashset) );

   /* dont create too small hashtables, i.e. at least size 32, and increase
    * the given size by divinding it by 0.9, since then no rebuilding will
    * be necessary if the given number of elements are inserted. Finally round
    * to the next power of two.
    */
   (*hashset)->shift = 64;
   (*hashset)->shift -= (int)ceil(log(MAX(8.0, size / 0.9)) / log(2.0));
   nslots = (uint32_t)SCIPhashsetGetNSlots(*hashset);
   (*hashset)->nelements = 0;

   SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &(*hashset)->slots, nslots) );

   return SCIP_OKAY;
}

/** frees the hash set */
void SCIPhashsetFree(
   SCIP_HASHSET**        hashset,            /**< pointer to the hash set */
   BMS_BLKMEM*           blkmem              /**< block memory used to store hash set entries */
   )
{
   BMSfreeBlockMemoryArray(blkmem, &(*hashset)->slots, SCIPhashsetGetNSlots(*hashset));
   BMSfreeBlockMemory(blkmem, hashset);
}

/** inserts new element into the hash set */
SCIP_RETCODE SCIPhashsetInsert(
   SCIP_HASHSET*         hashset,            /**< hash set */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash set entries */
   void*                 element             /**< element to insert */
   )
{
   assert(hashset != NULL);
   assert(hashset->slots != NULL);

   SCIP_CALL( hashsetCheckLoad(hashset, blkmem) );

   hashsetInsert(hashset, element);

   return SCIP_OKAY;
}

/** checks whether an element exists in the hash set */
SCIP_Bool SCIPhashsetExists(
   SCIP_HASHSET*         hashset,            /**< hash set */
   void*                 element             /**< element to search for */
   )
{
   uint32_t pos;
   uint32_t nslots;
   uint32_t mask;
   uint32_t elemdistance;

   assert(hashset != NULL);
   assert(hashset->slots != NULL);

   pos = hashSetDesiredPos(hashset, element);
   nslots = (uint32_t)SCIPhashsetGetNSlots(hashset);
   mask = nslots - 1;
   elemdistance = 0;

   while( TRUE ) /*lint !e716*/
   {
      uint32_t distance;

      /* found element */
      if( hashset->slots[pos] == element )
         return TRUE;

      /* slots is empty so element cannot be contained */
      if( hashset->slots[pos] == NULL )
         return FALSE;

      distance = ELEM_DISTANCE(pos);
      /* element can not be contained since otherwise we would have swapped it with this one during insert */
      if( elemdistance > distance )
         return FALSE;

      pos = (pos + 1) & mask;
      ++elemdistance;
   }
}

/** removes an element from the hash set, if it exists */
SCIP_RETCODE SCIPhashsetRemove(
   SCIP_HASHSET*         hashset,            /**< hash set */
   void*                 element             /**< origin to remove from the list */
   )
{
   uint32_t pos;
   uint32_t nslots;
   uint32_t mask;
   uint32_t elemdistance;

   assert(hashset != NULL);
   assert(hashset->slots != NULL);
   assert(element != NULL);

   pos = hashSetDesiredPos(hashset, element);
   nslots = (uint32_t)SCIPhashsetGetNSlots(hashset);
   mask = nslots - 1;
   elemdistance = 0;

   while( TRUE ) /*lint !e716*/
   {
      uint32_t distance;

      /* found element */
      if( hashset->slots[pos] == element )
         break;

      /* slots is empty so element cannot be contained */
      if( hashset->slots[pos] == NULL )
         return SCIP_OKAY;

      distance = ELEM_DISTANCE(pos);
      /* element can not be contained since otherwise we would have swapped it with this one during insert */
      if( elemdistance > distance )
         return SCIP_OKAY;

      pos = (pos + 1) & mask;
      ++elemdistance;
   }

   assert(hashset->slots[pos] == element);
   assert(SCIPhashsetExists(hashset, element));

   /* remove element */
   --hashset->nelements;

   /* move other elements if necessary */
   while( TRUE ) /*lint !e716*/
   {
      uint32_t nextpos = (pos + 1) & mask;

      /* nothing to do since there is no chain that needs to be moved */
      if( hashset->slots[nextpos] == NULL )
      {
         hashset->slots[pos] = NULL;
         assert(!SCIPhashsetExists(hashset, element));
         return SCIP_OKAY;
      }

      /* check if the element is the start of a new chain and return if that is the case */
      if( hashSetDesiredPos(hashset, hashset->slots[nextpos]) == nextpos )
      {
         hashset->slots[pos] = NULL;
         assert(!SCIPhashsetExists(hashset, element));
         return SCIP_OKAY;
      }

      /* element should be moved to the left and next element needs to be checked */
      hashset->slots[pos] = hashset->slots[nextpos];

      pos = nextpos;
   }
}

/** prints statistics about hash set usage */
void SCIPhashsetPrintStatistics(
   SCIP_HASHSET*         hashset,            /**< hash set */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   uint32_t maxprobelen = 0;
   uint64_t probelensum = 0;
   uint32_t nslots;
   uint32_t mask;
   uint32_t i;

   assert(hashset != NULL);

   nslots = (uint32_t)SCIPhashsetGetNSlots(hashset);
   mask = nslots - 1;

   /* compute the maximum and average probe length */
   for( i = 0; i < nslots; ++i )
   {
      if( hashset->slots[i] != NULL )
      {
         uint32_t probelen = ((hashSetDesiredPos(hashset, hashset->slots[i]) + nslots - i) & mask) + 1;
         probelensum += probelen;
         maxprobelen = MAX(probelen, maxprobelen);
      }
   }

   /* print general hash set statistics */
   SCIPmessagePrintInfo(messagehdlr, "%u hash entries, used %u/%u slots (%.1f%%)",
                        (unsigned int)hashset->nelements, (unsigned int)hashset->nelements,
                        (unsigned int)nslots, 100.0*(SCIP_Real)hashset->nelements/(SCIP_Real)(nslots));

   /* if not empty print average and maximum probe length */
   if( hashset->nelements > 0 )
      SCIPmessagePrintInfo(messagehdlr, ", avg. probe length is %.1f, max. probe length is %u",
         (SCIP_Real)(probelensum)/(SCIP_Real)hashset->nelements, (unsigned int)maxprobelen);
   SCIPmessagePrintInfo(messagehdlr, "\n");
}

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPhashsetIsEmpty
#undef SCIPhashsetGetNElements
#undef SCIPhashsetGetNSlots
#undef SCIPhashsetGetSlots

/** indicates whether a hash set has no entries */
SCIP_Bool SCIPhashsetIsEmpty(
   SCIP_HASHSET*         hashset             /**< hash set */
   )
{
   return hashset->nelements == 0;
}

/** gives the number of elements in a hash set */
int SCIPhashsetGetNElements(
   SCIP_HASHSET*         hashset             /**< hash set */
   )
{
   return (int)hashset->nelements;
}

/** gives the number of slots of a hash set */
int SCIPhashsetGetNSlots(
   SCIP_HASHSET*         hashset             /**< hash set */
   )
{
   return (int) (1u << (64 - hashset->shift));
}

/** gives the array of hash set slots; contains all elements in indetermined order and may contain NULL values */
void** SCIPhashsetGetSlots(
   SCIP_HASHSET*         hashset             /**< hash set */
   )
{
   return hashset->slots;
}

/** removes all entries in a hash set. */
void SCIPhashsetRemoveAll(
   SCIP_HASHSET*         hashset             /**< hash set */
   )
{
   BMSclearMemoryArray(hashset->slots, SCIPhashsetGetNSlots(hashset));

   hashset->nelements = 0;
}

/*
 * Dynamic Arrays
 */

/** creates a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayCreate(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(realarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, realarray) );
   (*realarray)->blkmem = blkmem;
   (*realarray)->vals = NULL;
   (*realarray)->valssize = 0;
   (*realarray)->firstidx = -1;
   (*realarray)->minusedidx = INT_MAX;
   (*realarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayCopy(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REALARRAY*       sourcerealarray     /**< dynamic real array to copy */
   )
{
   assert(realarray != NULL);
   assert(sourcerealarray != NULL);

   SCIP_CALL( SCIPrealarrayCreate(realarray, blkmem) );
   if( sourcerealarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*realarray)->vals, sourcerealarray->vals, \
                     sourcerealarray->valssize) );
   }
   (*realarray)->valssize = sourcerealarray->valssize;
   (*realarray)->firstidx = sourcerealarray->firstidx;
   (*realarray)->minusedidx = sourcerealarray->minusedidx;
   (*realarray)->maxusedidx = sourcerealarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayFree(
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   assert(realarray != NULL);
   assert(*realarray != NULL);

   BMSfreeBlockMemoryArrayNull((*realarray)->blkmem, &(*realarray)->vals, (*realarray)->valssize);
   BMSfreeBlockMemory((*realarray)->blkmem, realarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPrealarrayExtend(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(realarray != NULL);
   assert(realarray->minusedidx == INT_MAX || realarray->firstidx >= 0);
   assert(realarray->maxusedidx == INT_MIN || realarray->firstidx >= 0);
   assert(realarray->minusedidx == INT_MAX || realarray->minusedidx >= realarray->firstidx);
   assert(realarray->maxusedidx == INT_MIN || realarray->maxusedidx < realarray->firstidx + realarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   minidx = MIN(minidx, realarray->minusedidx);
   maxidx = MAX(maxidx, realarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending realarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n",
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > realarray->valssize )
   {
      SCIP_Real* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(realarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( realarray->firstidx != -1 )
      {
         for( i = 0; i < realarray->minusedidx - newfirstidx; ++i )
            newvals[i] = 0.0;

         /* check for possible overflow or negative value */
         assert(realarray->maxusedidx - realarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[realarray->minusedidx - newfirstidx],
            &(realarray->vals[realarray->minusedidx - realarray->firstidx]),
            realarray->maxusedidx - realarray->minusedidx + 1); /*lint !e866 !e776*/
         for( i = realarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = 0.0;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = 0.0;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(realarray->blkmem, &realarray->vals, realarray->valssize);
      realarray->vals = newvals;
      realarray->valssize = newvalssize;
      realarray->firstidx = newfirstidx;
   }
   else if( realarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      realarray->firstidx = minidx - nfree/2;
      assert(realarray->firstidx <= minidx);
      assert(maxidx < realarray->firstidx + realarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < realarray->valssize; ++i )
         assert(realarray->vals[i] == 0.0);
#endif
   }
   else if( minidx < realarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + realarray->valssize);

      if( realarray->minusedidx <= realarray->maxusedidx )
      {
         int shift;

         assert(realarray->firstidx <= realarray->minusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);

         /* shift used part of array to the right */
         shift = realarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = realarray->maxusedidx - realarray->firstidx; i >= realarray->minusedidx - realarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < realarray->valssize);
            realarray->vals[i + shift] = realarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            realarray->vals[realarray->minusedidx - realarray->firstidx + i] = 0.0;
      }
      realarray->firstidx = newfirstidx;
   }
   else if( maxidx >= realarray->firstidx + realarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + realarray->valssize);

      if( realarray->minusedidx <= realarray->maxusedidx )
      {
         int shift;

         assert(realarray->firstidx <= realarray->minusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - realarray->firstidx;
         assert(shift > 0);
         for( i = realarray->minusedidx - realarray->firstidx; i <= realarray->maxusedidx - realarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < realarray->valssize);
            realarray->vals[i - shift] = realarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            realarray->vals[realarray->maxusedidx - realarray->firstidx - i] = 0.0;
      }
      realarray->firstidx = newfirstidx;
   }

   assert(minidx >= realarray->firstidx);
   assert(maxidx < realarray->firstidx + realarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic real array */
SCIP_RETCODE SCIPrealarrayClear(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   SCIPdebugMessage("clearing realarray %p (firstidx=%d, size=%d, range=[%d,%d])\n",
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx);

   if( realarray->minusedidx <= realarray->maxusedidx )
   {
      assert(realarray->firstidx <= realarray->minusedidx);
      assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
      assert(realarray->firstidx != -1);
      assert(realarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&realarray->vals[realarray->minusedidx - realarray->firstidx],
         realarray->maxusedidx - realarray->minusedidx + 1); /*lint !e866*/

      /* mark the array cleared */
      realarray->minusedidx = INT_MAX;
      realarray->maxusedidx = INT_MIN;
   }
   assert(realarray->minusedidx == INT_MAX);
   assert(realarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
SCIP_Real SCIPrealarrayGetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(realarray != NULL);
   assert(idx >= 0);

   if( idx < realarray->minusedidx || idx > realarray->maxusedidx )
      return 0.0;
   else
   {
      assert(realarray->vals != NULL);
      assert(idx - realarray->firstidx >= 0);
      assert(idx - realarray->firstidx < realarray->valssize);

      return realarray->vals[idx - realarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPrealarraySetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   )
{
   assert(realarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting realarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %g\n",
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx, idx, val);

   if( val != 0.0 )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPrealarrayExtend(realarray, arraygrowinit, arraygrowfac, idx, idx) );
      assert(idx >= realarray->firstidx);
      assert(idx < realarray->firstidx + realarray->valssize);

      /* set the array value of the index */
      realarray->vals[idx - realarray->firstidx] = val;

      /* update min/maxusedidx */
      realarray->minusedidx = MIN(realarray->minusedidx, idx);
      realarray->maxusedidx = MAX(realarray->maxusedidx, idx);
   }
   else if( idx >= realarray->firstidx && idx < realarray->firstidx + realarray->valssize )
   {
      /* set the array value of the index to zero */
      realarray->vals[idx - realarray->firstidx] = 0.0;

      /* check, if we can tighten the min/maxusedidx */
      if( idx == realarray->minusedidx )
      {
         assert(realarray->maxusedidx >= 0);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
         do
         {
            realarray->minusedidx++;
         }
         while( realarray->minusedidx <= realarray->maxusedidx
            && realarray->vals[realarray->minusedidx - realarray->firstidx] == 0.0 );

         if( realarray->minusedidx > realarray->maxusedidx )
         {
            realarray->minusedidx = INT_MAX;
            realarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == realarray->maxusedidx )
      {
         assert(realarray->minusedidx >= 0);
         assert(realarray->minusedidx < realarray->maxusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
         do
         {
            realarray->maxusedidx--;
            assert(realarray->minusedidx <= realarray->maxusedidx);
         }
         while( realarray->vals[realarray->maxusedidx - realarray->firstidx] == 0.0 );
      }
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPrealarrayIncVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   )
{
   SCIP_Real oldval;

   oldval = SCIPrealarrayGetVal(realarray, idx);
   if( oldval != SCIP_INVALID ) /*lint !e777*/
      return SCIPrealarraySetVal(realarray, arraygrowinit, arraygrowfac, idx, oldval + incval);
   else
      return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPrealarrayGetMinIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   return realarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPrealarrayGetMaxIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   return realarray->maxusedidx;
}

/** creates a dynamic array of int values */
SCIP_RETCODE SCIPintarrayCreate(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the int array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(intarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, intarray) );
   (*intarray)->blkmem = blkmem;
   (*intarray)->vals = NULL;
   (*intarray)->valssize = 0;
   (*intarray)->firstidx = -1;
   (*intarray)->minusedidx = INT_MAX;
   (*intarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of int values */
SCIP_RETCODE SCIPintarrayCopy(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the copied int array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_INTARRAY*        sourceintarray      /**< dynamic int array to copy */
   )
{
   assert(intarray != NULL);
   assert(sourceintarray != NULL);

   SCIP_CALL( SCIPintarrayCreate(intarray, blkmem) );
   if( sourceintarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*intarray)->vals, sourceintarray->vals, sourceintarray->valssize) );
   }
   (*intarray)->valssize = sourceintarray->valssize;
   (*intarray)->firstidx = sourceintarray->firstidx;
   (*intarray)->minusedidx = sourceintarray->minusedidx;
   (*intarray)->maxusedidx = sourceintarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of int values */
SCIP_RETCODE SCIPintarrayFree(
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   assert(intarray != NULL);
   assert(*intarray != NULL);

   BMSfreeBlockMemoryArrayNull((*intarray)->blkmem, &(*intarray)->vals, (*intarray)->valssize);
   BMSfreeBlockMemory((*intarray)->blkmem, intarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPintarrayExtend(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(intarray != NULL);
   assert(intarray->minusedidx == INT_MAX || intarray->firstidx >= 0);
   assert(intarray->maxusedidx == INT_MIN || intarray->firstidx >= 0);
   assert(intarray->minusedidx == INT_MAX || intarray->minusedidx >= intarray->firstidx);
   assert(intarray->maxusedidx == INT_MIN || intarray->maxusedidx < intarray->firstidx + intarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   minidx = MIN(minidx, intarray->minusedidx);
   maxidx = MAX(maxidx, intarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending intarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n",
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > intarray->valssize )
   {
      int* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(intarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( intarray->firstidx != -1 )
      {
         for( i = 0; i < intarray->minusedidx - newfirstidx; ++i )
            newvals[i] = 0;

         /* check for possible overflow or negative value */
         assert(intarray->maxusedidx - intarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[intarray->minusedidx - newfirstidx],
            &intarray->vals[intarray->minusedidx - intarray->firstidx],
            intarray->maxusedidx - intarray->minusedidx + 1); /*lint !e866 !e776*/
         for( i = intarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = 0;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = 0;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(intarray->blkmem, &intarray->vals, intarray->valssize);
      intarray->vals = newvals;
      intarray->valssize = newvalssize;
      intarray->firstidx = newfirstidx;
   }
   else if( intarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      intarray->firstidx = minidx - nfree/2;
      assert(intarray->firstidx <= minidx);
      assert(maxidx < intarray->firstidx + intarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < intarray->valssize; ++i )
         assert(intarray->vals[i] == 0);
#endif
   }
   else if( minidx < intarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + intarray->valssize);

      if( intarray->minusedidx <= intarray->maxusedidx )
      {
         int shift;

         assert(intarray->firstidx <= intarray->minusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);

         /* shift used part of array to the right */
         shift = intarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = intarray->maxusedidx - intarray->firstidx; i >= intarray->minusedidx - intarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < intarray->valssize);
            intarray->vals[i + shift] = intarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            intarray->vals[intarray->minusedidx - intarray->firstidx + i] = 0;
      }
      intarray->firstidx = newfirstidx;
   }
   else if( maxidx >= intarray->firstidx + intarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + intarray->valssize);

      if( intarray->minusedidx <= intarray->maxusedidx )
      {
         int shift;

         assert(intarray->firstidx <= intarray->minusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - intarray->firstidx;
         assert(shift > 0);
         for( i = intarray->minusedidx - intarray->firstidx; i <= intarray->maxusedidx - intarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < intarray->valssize);
            intarray->vals[i - shift] = intarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            intarray->vals[intarray->maxusedidx - intarray->firstidx - i] = 0;
      }
      intarray->firstidx = newfirstidx;
   }

   assert(minidx >= intarray->firstidx);
   assert(maxidx < intarray->firstidx + intarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic int array */
SCIP_RETCODE SCIPintarrayClear(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   SCIPdebugMessage("clearing intarray %p (firstidx=%d, size=%d, range=[%d,%d])\n",
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx);

   if( intarray->minusedidx <= intarray->maxusedidx )
   {
      assert(intarray->firstidx <= intarray->minusedidx);
      assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
      assert(intarray->firstidx != -1);
      assert(intarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&intarray->vals[intarray->minusedidx - intarray->firstidx],
         intarray->maxusedidx - intarray->minusedidx + 1); /*lint !e866*/

      /* mark the array cleared */
      intarray->minusedidx = INT_MAX;
      intarray->maxusedidx = INT_MIN;
   }
   assert(intarray->minusedidx == INT_MAX);
   assert(intarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
int SCIPintarrayGetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(intarray != NULL);
   assert(idx >= 0);

   if( idx < intarray->minusedidx || idx > intarray->maxusedidx )
      return 0;
   else
   {
      assert(intarray->vals != NULL);
      assert(idx - intarray->firstidx >= 0);
      assert(idx - intarray->firstidx < intarray->valssize);

      return intarray->vals[idx - intarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPintarraySetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   )
{
   assert(intarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting intarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %d\n",
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, idx, val);

   if( val != 0 )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPintarrayExtend(intarray, arraygrowinit, arraygrowfac, idx, idx) );
      assert(idx >= intarray->firstidx);
      assert(idx < intarray->firstidx + intarray->valssize);

      /* set the array value of the index */
      intarray->vals[idx - intarray->firstidx] = val;

      /* update min/maxusedidx */
      intarray->minusedidx = MIN(intarray->minusedidx, idx);
      intarray->maxusedidx = MAX(intarray->maxusedidx, idx);
   }
   else if( idx >= intarray->firstidx && idx < intarray->firstidx + intarray->valssize )
   {
      /* set the array value of the index to zero */
      intarray->vals[idx - intarray->firstidx] = 0;

      /* check, if we can tighten the min/maxusedidx */
      if( idx == intarray->minusedidx )
      {
         assert(intarray->maxusedidx >= 0);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
         do
         {
            intarray->minusedidx++;
         }
         while( intarray->minusedidx <= intarray->maxusedidx
            && intarray->vals[intarray->minusedidx - intarray->firstidx] == 0 );
         if( intarray->minusedidx > intarray->maxusedidx )
         {
            intarray->minusedidx = INT_MAX;
            intarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == intarray->maxusedidx )
      {
         assert(intarray->minusedidx >= 0);
         assert(intarray->minusedidx < intarray->maxusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
         do
         {
            intarray->maxusedidx--;
            assert(intarray->minusedidx <= intarray->maxusedidx);
         }
         while( intarray->vals[intarray->maxusedidx - intarray->firstidx] == 0 );
      }
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPintarrayIncVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   )
{
   return SCIPintarraySetVal(intarray, arraygrowinit, arraygrowfac, idx, SCIPintarrayGetVal(intarray, idx) + incval);
}

/** returns the minimal index of all stored non-zero elements */
int SCIPintarrayGetMinIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   return intarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPintarrayGetMaxIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   return intarray->maxusedidx;
}


/** creates a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayCreate(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the bool array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(boolarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, boolarray) );
   (*boolarray)->blkmem = blkmem;
   (*boolarray)->vals = NULL;
   (*boolarray)->valssize = 0;
   (*boolarray)->firstidx = -1;
   (*boolarray)->minusedidx = INT_MAX;
   (*boolarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayCopy(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the copied bool array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_BOOLARRAY*       sourceboolarray     /**< dynamic bool array to copy */
   )
{
   assert(boolarray != NULL);
   assert(sourceboolarray != NULL);

   SCIP_CALL( SCIPboolarrayCreate(boolarray, blkmem) );
   if( sourceboolarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*boolarray)->vals, sourceboolarray->vals, \
                     sourceboolarray->valssize) );
   }
   (*boolarray)->valssize = sourceboolarray->valssize;
   (*boolarray)->firstidx = sourceboolarray->firstidx;
   (*boolarray)->minusedidx = sourceboolarray->minusedidx;
   (*boolarray)->maxusedidx = sourceboolarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayFree(
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   assert(boolarray != NULL);
   assert(*boolarray != NULL);

   BMSfreeBlockMemoryArrayNull((*boolarray)->blkmem, &(*boolarray)->vals, (*boolarray)->valssize);
   BMSfreeBlockMemory((*boolarray)->blkmem, boolarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPboolarrayExtend(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(boolarray != NULL);
   assert(boolarray->minusedidx == INT_MAX || boolarray->firstidx >= 0);
   assert(boolarray->maxusedidx == INT_MIN || boolarray->firstidx >= 0);
   assert(boolarray->minusedidx == INT_MAX || boolarray->minusedidx >= boolarray->firstidx);
   assert(boolarray->maxusedidx == INT_MIN || boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   minidx = MIN(minidx, boolarray->minusedidx);
   maxidx = MAX(maxidx, boolarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n",
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > boolarray->valssize )
   {
      SCIP_Bool* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(boolarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( boolarray->firstidx != -1 )
      {
         for( i = 0; i < boolarray->minusedidx - newfirstidx; ++i )
            newvals[i] = FALSE;

         /* check for possible overflow or negative value */
         assert(boolarray->maxusedidx - boolarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[boolarray->minusedidx - newfirstidx],
            &boolarray->vals[boolarray->minusedidx - boolarray->firstidx],
            boolarray->maxusedidx - boolarray->minusedidx + 1); /*lint !e866 !e776*/
         for( i = boolarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = FALSE;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = FALSE;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(boolarray->blkmem, &boolarray->vals, boolarray->valssize);
      boolarray->vals = newvals;
      boolarray->valssize = newvalssize;
      boolarray->firstidx = newfirstidx;
   }
   else if( boolarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      boolarray->firstidx = minidx - nfree/2;
      assert(boolarray->firstidx <= minidx);
      assert(maxidx < boolarray->firstidx + boolarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < boolarray->valssize; ++i )
         assert(boolarray->vals[i] == FALSE);
#endif
   }
   else if( minidx < boolarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + boolarray->valssize);

      if( boolarray->minusedidx <= boolarray->maxusedidx )
      {
         int shift;

         assert(boolarray->firstidx <= boolarray->minusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);

         /* shift used part of array to the right */
         shift = boolarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = boolarray->maxusedidx - boolarray->firstidx; i >= boolarray->minusedidx - boolarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < boolarray->valssize);
            boolarray->vals[i + shift] = boolarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            boolarray->vals[boolarray->minusedidx - boolarray->firstidx + i] = FALSE;
      }
      boolarray->firstidx = newfirstidx;
   }
   else if( maxidx >= boolarray->firstidx + boolarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + boolarray->valssize);

      if( boolarray->minusedidx <= boolarray->maxusedidx )
      {
         int shift;

         assert(boolarray->firstidx <= boolarray->minusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - boolarray->firstidx;
         assert(shift > 0);

	 assert(0 <= boolarray->minusedidx - boolarray->firstidx - shift);
         assert(boolarray->maxusedidx - boolarray->firstidx - shift < boolarray->valssize);
	 BMSmoveMemoryArray(&(boolarray->vals[boolarray->minusedidx - boolarray->firstidx - shift]),
            &(boolarray->vals[boolarray->minusedidx - boolarray->firstidx]),
            boolarray->maxusedidx - boolarray->minusedidx + 1); /*lint !e866*/

         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            boolarray->vals[boolarray->maxusedidx - boolarray->firstidx - i] = FALSE;
      }
      boolarray->firstidx = newfirstidx;
   }

   assert(minidx >= boolarray->firstidx);
   assert(maxidx < boolarray->firstidx + boolarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic bool array */
SCIP_RETCODE SCIPboolarrayClear(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   SCIPdebugMessage("clearing boolarray %p (firstidx=%d, size=%d, range=[%d,%d])\n",
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx);

   if( boolarray->minusedidx <= boolarray->maxusedidx )
   {
      assert(boolarray->firstidx <= boolarray->minusedidx);
      assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
      assert(boolarray->firstidx != -1);
      assert(boolarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&boolarray->vals[boolarray->minusedidx - boolarray->firstidx],
         boolarray->maxusedidx - boolarray->minusedidx + 1); /*lint !e866*/

      /* mark the array cleared */
      boolarray->minusedidx = INT_MAX;
      boolarray->maxusedidx = INT_MIN;
   }
   assert(boolarray->minusedidx == INT_MAX);
   assert(boolarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
SCIP_Bool SCIPboolarrayGetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(boolarray != NULL);
   assert(idx >= 0);

   if( idx < boolarray->minusedidx || idx > boolarray->maxusedidx )
      return FALSE;
   else
   {
      assert(boolarray->vals != NULL);
      assert(idx - boolarray->firstidx >= 0);
      assert(idx - boolarray->firstidx < boolarray->valssize);

      return boolarray->vals[idx - boolarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPboolarraySetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   )
{
   assert(boolarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %u\n",
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, idx, val);

   if( val != FALSE )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPboolarrayExtend(boolarray, arraygrowinit, arraygrowfac, idx, idx) );
      assert(idx >= boolarray->firstidx);
      assert(idx < boolarray->firstidx + boolarray->valssize);

      /* set the array value of the index */
      boolarray->vals[idx - boolarray->firstidx] = val;

      /* update min/maxusedidx */
      boolarray->minusedidx = MIN(boolarray->minusedidx, idx);
      boolarray->maxusedidx = MAX(boolarray->maxusedidx, idx);
   }
   else if( idx >= boolarray->firstidx && idx < boolarray->firstidx + boolarray->valssize )
   {
      /* set the array value of the index to zero */
      boolarray->vals[idx - boolarray->firstidx] = FALSE;

      /* check, if we can tighten the min/maxusedidx */
      if( idx == boolarray->minusedidx )
      {
         assert(boolarray->maxusedidx >= 0);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
         do
         {
            boolarray->minusedidx++;
         }
         while( boolarray->minusedidx <= boolarray->maxusedidx
            && boolarray->vals[boolarray->minusedidx - boolarray->firstidx] == FALSE );
         if( boolarray->minusedidx > boolarray->maxusedidx )
         {
            boolarray->minusedidx = INT_MAX;
            boolarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == boolarray->maxusedidx )
      {
         assert(boolarray->minusedidx >= 0);
         assert(boolarray->minusedidx < boolarray->maxusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
         do
         {
            boolarray->maxusedidx--;
            assert(boolarray->minusedidx <= boolarray->maxusedidx);
         }
         while( boolarray->vals[boolarray->maxusedidx - boolarray->firstidx] == FALSE );
      }
   }

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPboolarrayGetMinIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   return boolarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPboolarrayGetMaxIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   return boolarray->maxusedidx;
}


/** creates a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayCreate(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the ptr array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(ptrarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, ptrarray) );
   (*ptrarray)->blkmem = blkmem;
   (*ptrarray)->vals = NULL;
   (*ptrarray)->valssize = 0;
   (*ptrarray)->firstidx = -1;
   (*ptrarray)->minusedidx = INT_MAX;
   (*ptrarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayCopy(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the copied ptr array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PTRARRAY*        sourceptrarray      /**< dynamic ptr array to copy */
   )
{
   assert(ptrarray != NULL);
   assert(sourceptrarray != NULL);

   SCIP_CALL( SCIPptrarrayCreate(ptrarray, blkmem) );
   if( sourceptrarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*ptrarray)->vals, sourceptrarray->vals, sourceptrarray->valssize) );
   }
   (*ptrarray)->valssize = sourceptrarray->valssize;
   (*ptrarray)->firstidx = sourceptrarray->firstidx;
   (*ptrarray)->minusedidx = sourceptrarray->minusedidx;
   (*ptrarray)->maxusedidx = sourceptrarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayFree(
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the ptr array */
   )
{
   assert(ptrarray != NULL);
   assert(*ptrarray != NULL);

   BMSfreeBlockMemoryArrayNull((*ptrarray)->blkmem, &(*ptrarray)->vals, (*ptrarray)->valssize);
   BMSfreeBlockMemory((*ptrarray)->blkmem, ptrarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPptrarrayExtend(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(ptrarray != NULL);
   assert(ptrarray->minusedidx == INT_MAX || ptrarray->firstidx >= 0);
   assert(ptrarray->maxusedidx == INT_MIN || ptrarray->firstidx >= 0);
   assert(ptrarray->minusedidx == INT_MAX || ptrarray->minusedidx >= ptrarray->firstidx);
   assert(ptrarray->maxusedidx == INT_MIN || ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   minidx = MIN(minidx, ptrarray->minusedidx);
   maxidx = MAX(maxidx, ptrarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending ptrarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n",
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > ptrarray->valssize )
   {
      void** newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(ptrarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( ptrarray->firstidx != -1 )
      {
         for( i = 0; i < ptrarray->minusedidx - newfirstidx; ++i )
            newvals[i] = NULL;

         /* check for possible overflow or negative value */
         assert(ptrarray->maxusedidx - ptrarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[ptrarray->minusedidx - newfirstidx],
            &(ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx]),
            ptrarray->maxusedidx - ptrarray->minusedidx + 1); /*lint !e866 !e776*/
         for( i = ptrarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = NULL;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = NULL;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(ptrarray->blkmem, &ptrarray->vals, ptrarray->valssize);
      ptrarray->vals = newvals;
      ptrarray->valssize = newvalssize;
      ptrarray->firstidx = newfirstidx;
   }
   else if( ptrarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      ptrarray->firstidx = minidx - nfree/2;
      assert(ptrarray->firstidx <= minidx);
      assert(maxidx < ptrarray->firstidx + ptrarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < ptrarray->valssize; ++i )
         assert(ptrarray->vals[i] == NULL);
#endif
   }
   else if( minidx < ptrarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + ptrarray->valssize);

      if( ptrarray->minusedidx <= ptrarray->maxusedidx )
      {
         int shift;

         assert(ptrarray->firstidx <= ptrarray->minusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);

         /* shift used part of array to the right */
         shift = ptrarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = ptrarray->maxusedidx - ptrarray->firstidx; i >= ptrarray->minusedidx - ptrarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < ptrarray->valssize);
            ptrarray->vals[i + shift] = ptrarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx + i] = NULL;
      }
      ptrarray->firstidx = newfirstidx;
   }
   else if( maxidx >= ptrarray->firstidx + ptrarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + ptrarray->valssize);

      if( ptrarray->minusedidx <= ptrarray->maxusedidx )
      {
         int shift;

         assert(ptrarray->firstidx <= ptrarray->minusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - ptrarray->firstidx;
         assert(shift > 0);
         for( i = ptrarray->minusedidx - ptrarray->firstidx; i <= ptrarray->maxusedidx - ptrarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < ptrarray->valssize);
            ptrarray->vals[i - shift] = ptrarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            ptrarray->vals[ptrarray->maxusedidx - ptrarray->firstidx - i] = NULL;
      }
      ptrarray->firstidx = newfirstidx;
   }

   assert(minidx >= ptrarray->firstidx);
   assert(maxidx < ptrarray->firstidx + ptrarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic pointer array */
SCIP_RETCODE SCIPptrarrayClear(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   SCIPdebugMessage("clearing ptrarray %p (firstidx=%d, size=%d, range=[%d,%d])\n",
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx);

   if( ptrarray->minusedidx <= ptrarray->maxusedidx )
   {
      assert(ptrarray->firstidx <= ptrarray->minusedidx);
      assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
      assert(ptrarray->firstidx != -1);
      assert(ptrarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx],
         ptrarray->maxusedidx - ptrarray->minusedidx + 1); /*lint !e866*/

      /* mark the array cleared */
      ptrarray->minusedidx = INT_MAX;
      ptrarray->maxusedidx = INT_MIN;
   }
   assert(ptrarray->minusedidx == INT_MAX);
   assert(ptrarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
void* SCIPptrarrayGetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(ptrarray != NULL);
   assert(idx >= 0);

   if( idx < ptrarray->minusedidx || idx > ptrarray->maxusedidx )
      return NULL;
   else
   {
      assert(ptrarray->vals != NULL);
      assert(idx - ptrarray->firstidx >= 0);
      assert(idx - ptrarray->firstidx < ptrarray->valssize);

      return ptrarray->vals[idx - ptrarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPptrarraySetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   )
{
   assert(ptrarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting ptrarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %p\n",
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx, idx, val);

   if( val != NULL )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPptrarrayExtend(ptrarray, arraygrowinit, arraygrowfac, idx, idx) );
      assert(idx >= ptrarray->firstidx);
      assert(idx < ptrarray->firstidx + ptrarray->valssize);

      /* set the array value of the index */
      ptrarray->vals[idx - ptrarray->firstidx] = val;

      /* update min/maxusedidx */
      ptrarray->minusedidx = MIN(ptrarray->minusedidx, idx);
      ptrarray->maxusedidx = MAX(ptrarray->maxusedidx, idx);
   }
   else if( idx >= ptrarray->firstidx && idx < ptrarray->firstidx + ptrarray->valssize )
   {
      /* set the array value of the index to zero */
      ptrarray->vals[idx - ptrarray->firstidx] = NULL;

      /* check, if we can tighten the min/maxusedidx */
      if( idx == ptrarray->minusedidx )
      {
         assert(ptrarray->maxusedidx >= 0);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
         do
         {
            ptrarray->minusedidx++;
         }
         while( ptrarray->minusedidx <= ptrarray->maxusedidx
            && ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx] == NULL );
         if( ptrarray->minusedidx > ptrarray->maxusedidx )
         {
            ptrarray->minusedidx = INT_MAX;
            ptrarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == ptrarray->maxusedidx )
      {
         assert(ptrarray->minusedidx >= 0);
         assert(ptrarray->minusedidx < ptrarray->maxusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
         do
         {
            ptrarray->maxusedidx--;
            assert(ptrarray->minusedidx <= ptrarray->maxusedidx);
         }
         while( ptrarray->vals[ptrarray->maxusedidx - ptrarray->firstidx] == NULL );
      }
   }

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPptrarrayGetMinIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   return ptrarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPptrarrayGetMaxIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   return ptrarray->maxusedidx;
}


/*
 * Sorting algorithms
 */

/** default comparer for integers */
SCIP_DECL_SORTPTRCOMP(SCIPsortCompInt)
{
   int value1;
   int value2;

   value1 = (int)(size_t)elem1;
   value2 = (int)(size_t)elem2;

   if( value1 < value2 )
      return -1;

   if( value2 < value1 )
      return 1;

   return 0;
}

/* first all upwards-sorting methods */

/** sort an indexed element set in non-decreasing order, resulting in a permutation index array */
void SCIPsort(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   )
{
   int pos;

   assert(indcomp != NULL);
   assert(len == 0 || perm != NULL);

   /* create identity permutation */
   for( pos = 0; pos < len; ++pos )
      perm[pos] = pos;

   SCIPsortInd(perm, indcomp, dataptr, len);
}

/* SCIPsortInd(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Ind
#define SORTTPL_KEYTYPE     int
#define SORTTPL_INDCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Ptr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortPtrRealRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrRealRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortPtrRealBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrRealBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrRealIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrRealBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrRealBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrLongInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrLongInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrLongIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrLongIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Real
#define SORTTPL_KEYTYPE     SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealIntLong
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealIntPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealLongRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealLongRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Longint
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealPtrPtrIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealRealBoolBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealBoolBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_FIELD5TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Int
#define SORTTPL_KEYTYPE     int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntLong
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortIntRealLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntRealLong
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtrReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortIntIntIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortIntPtrIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtrIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Long
#define SORTTPL_KEYTYPE     SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtr
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrRealBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrRealBool
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrRealRealBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrRealRealBool
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Real
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrRealRealIntBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrRealRealIntBool
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Real
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_FIELD5TYPE  SCIP_Bool
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrIntInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrBoolInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrBoolInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrIntIntBoolBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrIntIntBoolBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtrIntIntBoolBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtrIntIntBoolBool
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_FIELD5TYPE  SCIP_Bool
#include "scip/sorttpl.c" /*lint !e451*/


/* now all downwards-sorting methods */


/** sort an indexed element set in non-increasing order, resulting in a permutation index array */
void SCIPsortDown(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   )
{
   int pos;

   assert(indcomp != NULL);
   assert(len == 0 || perm != NULL);

   /* create identity permutation */
   for( pos = 0; pos < len; ++pos )
      perm[pos] = pos;

   SCIPsortDownInd(perm, indcomp, dataptr, len);
}


/* SCIPsortDownInd(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownInd
#define SORTTPL_KEYTYPE     int
#define SORTTPL_INDCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownPtrBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrRealBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrRealBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrRealIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrRealBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrRealBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrLongInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrLongInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrLongIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrLongIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownReal
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealIntLong
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealIntPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtrPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownRealRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownRealRealPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealPtrPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealLongRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealLongRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Longint
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownRealPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtrPtrIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealBoolBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealBoolBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_FIELD5TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntLong
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntPtrIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntPtrIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLong
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtr
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrRealBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrRealBool
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrRealRealBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrRealRealBool
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Real
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrRealRealIntBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrRealRealIntBool
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Real
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_FIELD5TYPE  SCIP_Bool
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrIntInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrBoolInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrBoolInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrIntIntBoolBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrIntIntBoolBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntPtrIntIntBoolBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntPtrIntIntBoolBool
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_FIELD5TYPE  SCIP_Bool
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/*
 * Resulting activity
 */

/** create a resource activity */
SCIP_RETCODE SCIPactivityCreate(
   SCIP_RESOURCEACTIVITY** activity,         /**< pointer to store the resource activity */
   SCIP_VAR*             var,                /**< start time variable of the activity */
   int                   duration,           /**< duration of the activity */
   int                   demand              /**< demand of the activity */
   )
{
   assert(activity != NULL);

   SCIP_ALLOC( BMSallocMemory(activity) );

   (*activity)->var = var;
   (*activity)->duration = duration;
   (*activity)->demand = demand;

   return SCIP_OKAY;
}

/** frees a resource activity */
void SCIPactivityFree(
   SCIP_RESOURCEACTIVITY** activity          /**< pointer to the resource activity */
   )
{
   assert(activity != NULL);
   assert(*activity != NULL);

   BMSfreeMemory(activity);
}

/* some simple variable functions implemented as defines */

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPactivityGetVar
#undef SCIPactivityGetDuration
#undef SCIPactivityGetDemand
#undef SCIPactivityGetEnergy

/** returns the start time variable of the resource activity */
SCIP_VAR* SCIPactivityGetVar(
   SCIP_RESOURCEACTIVITY* activity           /**< resource activity */
   )
{
   assert(activity != NULL);

   return activity->var;
}

/** returns the duration of the resource activity */
int SCIPactivityGetDuration(
   SCIP_RESOURCEACTIVITY* activity           /**< resource activity */
   )
{
   assert(activity != NULL);

   return activity->duration;
}

/** returns the demand of the resource activity */
int SCIPactivityGetDemand(
   SCIP_RESOURCEACTIVITY* activity           /**< resource activity */
   )
{
   assert(activity != NULL);

   return activity->demand;
}

/** returns the energy of the resource activity */
int SCIPactivityGetEnergy(
   SCIP_RESOURCEACTIVITY* activity           /**< resource activity */
   )
{
   assert(activity != NULL);

   return activity->duration * activity->demand ;
}

#endif

/*
 * Resource Profile
 */

/** helper method to create a profile */
static
SCIP_RETCODE doProfileCreate(
   SCIP_PROFILE**        profile,            /**< pointer to store the resource profile */
   int                   capacity            /**< resource capacity */
   )
{
   SCIP_ALLOC( BMSallocMemory(profile) );
   BMSclearMemory(*profile);

   (*profile)->arraysize = 10;
   SCIP_ALLOC( BMSallocMemoryArray(&(*profile)->timepoints, (*profile)->arraysize) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*profile)->loads, (*profile)->arraysize) );

   /* setup resource profile for use */
   (*profile)->ntimepoints = 1;
   (*profile)->timepoints[0] = 0;
   (*profile)->loads[0] = 0;
   (*profile)->capacity = capacity;

   return SCIP_OKAY;
}

/** creates resource profile */
SCIP_RETCODE SCIPprofileCreate(
   SCIP_PROFILE**        profile,            /**< pointer to store the resource profile */
   int                   capacity            /**< resource capacity */
   )
{
   assert(profile != NULL);
   assert(capacity > 0);

   SCIP_CALL_FINALLY( doProfileCreate(profile, capacity), SCIPprofileFree(profile) );

   return SCIP_OKAY;
}

/** frees given resource profile */
void SCIPprofileFree(
   SCIP_PROFILE**        profile             /**< pointer to the resource profile */
   )
{
   assert(profile != NULL);

   /* free resource profile */
   if( *profile != NULL )
   {
      BMSfreeMemoryArrayNull(&(*profile)->loads);
      BMSfreeMemoryArrayNull(&(*profile)->timepoints);
      BMSfreeMemory(profile);
   }
}

/** output of the given resource profile */
void SCIPprofilePrint(
   SCIP_PROFILE*         profile,            /**< resource profile to output */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int t;

   SCIPmessageFPrintInfo(messagehdlr, file, "Profile <%p> (capacity %d) --> ", profile, profile->capacity);

   for( t = 0; t < profile->ntimepoints; ++t )
   {
      if( t == 0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "%d:(%d,%d)", t, profile->timepoints[t], profile->loads[t]);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, ", %d:(%d,%d)", t, profile->timepoints[t], profile->loads[t]);
   }

   SCIPmessageFPrintInfo(messagehdlr, file,"\n");
}

/** returns the capacity of the resource profile */
int SCIPprofileGetCapacity(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   )
{
   assert(profile != NULL);

   return profile->capacity;
}

/** returns the number time points of the resource profile */
int SCIPprofileGetNTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   )
{
   assert(profile != NULL);

   return profile->ntimepoints;
}

/** returns the time points of the resource profile */
int* SCIPprofileGetTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   )
{
   assert(profile != NULL);

   return profile->timepoints;
}

/** returns the loads of the resource profile */
int* SCIPprofileGetLoads(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   )
{
   assert(profile != NULL);

   return profile->loads;
}

/** returns the time point for given position of the resource profile */
int SCIPprofileGetTime(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos                 /**< position */
   )
{
   assert(profile != NULL);
   assert(pos >= 0 && pos < profile->ntimepoints);

   return profile->timepoints[pos];
}

/** returns the loads of the resource profile at the given position */
int SCIPprofileGetLoad(
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   pos                 /**< position */
   )
{
   assert(profile != NULL);
   assert(pos >= 0 && pos < profile->ntimepoints);

   return profile->loads[pos];
}

/** returns if the given time point exists in the resource profile and stores the position of the given time point if it
 *  exists; otherwise the position of the next smaller existing time point is stored
 */
SCIP_Bool SCIPprofileFindLeft(
   SCIP_PROFILE*         profile,            /**< resource profile to search */
   int                   timepoint,          /**< time point to search for */
   int*                  pos                 /**< pointer to store the position */
   )
{
   assert(profile != NULL);
   assert(timepoint >= 0);
   assert(profile->ntimepoints > 0);
   assert(profile->timepoints[0] == 0);

   /* find the position of time point in the time points array via binary search */
   if( SCIPsortedvecFindInt(profile->timepoints, timepoint, profile->ntimepoints, pos) )
      return TRUE;

   assert(*pos > 0);
   (*pos)--;

   return FALSE;
}

/* ensures that resource profile arrays is big enough */
static
SCIP_RETCODE ensureProfileSize(
   SCIP_PROFILE*         profile,            /**< resource profile to insert the time point */
   int                   neededsize          /**< needed size */
   )
{
   assert(profile->arraysize > 0);

   /* check whether the arrays are big enough */
   if( neededsize <= profile->arraysize )
      return SCIP_OKAY;

   profile->arraysize *= 2;

   SCIP_ALLOC( BMSreallocMemoryArray(&profile->timepoints, profile->arraysize) );
   SCIP_ALLOC( BMSreallocMemoryArray(&profile->loads, profile->arraysize) );

   return SCIP_OKAY;
}

/** inserts the given time point into the resource profile if it this time point does not exists yet; returns its
 *  position in the time point array
 */
static
SCIP_RETCODE profileInsertTimepoint(
   SCIP_PROFILE*         profile,            /**< resource profile to insert the time point */
   int                   timepoint,          /**< time point to insert */
   int*                  pos                 /**< pointer to store the insert position */
   )
{
   assert(profile != NULL);
   assert(timepoint >= 0);
   assert(profile->arraysize >= profile->ntimepoints);

   /* get the position of the given time point in the resource profile array if it exists; otherwise the position of the
    * next smaller existing time point
    */
   if( !SCIPprofileFindLeft(profile, timepoint, pos) )
   {
      assert(*pos >= 0 && *pos < profile->ntimepoints);
      assert(timepoint >= profile->timepoints[*pos]);

      /* ensure that the arrays are big enough */
      SCIP_CALL( ensureProfileSize(profile, profile->ntimepoints + 1) );
      assert(profile->arraysize > profile->ntimepoints);

      /* insert new time point into the (sorted) resource profile */
      SCIPsortedvecInsertIntInt(profile->timepoints, profile->loads, timepoint, profile->loads[*pos],
         &profile->ntimepoints, pos);
   }

#ifndef NDEBUG
   /* check if the time points are sorted */
   {
      int i;
      for( i = 1; i < profile->ntimepoints; ++i )
         assert(profile->timepoints[i-1] < profile->timepoints[i]);
   }
#endif

   return SCIP_OKAY;
}

/** updates the resource profile due to inserting of a core */
static
SCIP_RETCODE profileUpdate(
   SCIP_PROFILE*         profile,            /**< resource profile to update */
   int                   left,               /**< left side of core interval */
   int                   right,              /**< right side of core interval */
   int                   demand,             /**< demand of the core */
   int*                  pos,                /**< pointer to store the first position were it gets infeasible */
   SCIP_Bool*            infeasible          /**< pointer to store if the update is infeasible */
   )
{
   int startpos;
   int endpos;
   int i;

   assert(profile != NULL);
   assert(profile->arraysize >= profile->ntimepoints);
   assert(left >= 0);
   assert(left < right);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;
   (*pos) = -1;

   /* get position of the starttime in profile */
   SCIP_CALL( profileInsertTimepoint(profile, left, &startpos) );
   assert(profile->timepoints[startpos] == left);

   /* get position of the endtime in profile */
   SCIP_CALL( profileInsertTimepoint(profile, right, &endpos) );
   assert(profile->timepoints[endpos] == right);

   assert(startpos < endpos);
   assert(profile->arraysize >= profile->ntimepoints);

   /* remove/add the given demand from the core */
   for( i = startpos; i < endpos; ++i )
   {
      profile->loads[i] += demand;

      /* check if the core fits */
      if( profile->loads[i] > profile->capacity )
      {
         SCIPdebugMessage("core insertion detected infeasibility (pos %d)\n", i);

         (*infeasible) = TRUE;
         (*pos) = i;

         /* remove the partly inserted core since it does fit completely */
         for( ; i >= startpos; --i ) /*lint !e445*/
            profile->loads[i] -= demand;

         break;
      }
   }

   return SCIP_OKAY;
}

/** insert a core into resource profile; if the core is non-empty the resource profile will be updated otherwise nothing
 *  happens
 */
SCIP_RETCODE SCIPprofileInsertCore(
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   demand,             /**< demand of the core */
   int*                  pos,                /**< pointer to store the first position were it gets infeasible */
   SCIP_Bool*            infeasible          /**< pointer to store if the core does not fit due to capacity */
   )
{
   assert(profile != NULL);
   assert(left < right);
   assert(demand >= 0);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;
   (*pos) = -1;

   /* insert core into the resource profile */
   SCIPdebugMessage("insert core [%d,%d] with demand %d\n", left, right, demand);

   if( demand > 0 )
   {
      /* try to insert core into the resource profile */
      SCIP_CALL( profileUpdate(profile, left, right, demand, pos, infeasible) );
   }

   return SCIP_OKAY;
}

/** subtracts the demand from the resource profile during core time */
SCIP_RETCODE SCIPprofileDeleteCore(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   demand              /**< demand of the core */
   )
{
   SCIP_Bool infeasible;
   int pos;

   assert(left < right);
#ifndef NDEBUG
   {
      /* check if the left and right time points of the core correspond to a time point in the resource profile; this
       * should be the case since we added the core before to the resource profile
       */
      assert(SCIPprofileFindLeft(profile, left, &pos));
      assert(SCIPprofileFindLeft(profile, right, &pos));
   }
#endif

   /* remove the core from the resource profile */
   SCIPdebugMessage("delete core [%d,%d] with demand %d\n", left, right, demand);

   SCIP_CALL( profileUpdate(profile, left, right, -demand, &pos, &infeasible) );
   assert(!infeasible);

   return SCIP_OKAY;
}

/** returns TRUE if the core (given by its demand and during) can be inserted at the given time point; otherwise FALSE */
static
int profileFindFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos,                /**< pointer to store the position in the profile to start the serch */
   int                   lst,                /**< latest start time */
   int                   duration,           /**< duration of the core */
   int                   demand,             /**< demand of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the corer cannot be inserted */
   )
{
   int remainingduration;
   int startpos;

   assert(profile != NULL);
   assert(pos >= 0);
   assert(pos < profile->ntimepoints);
   assert(duration > 0);
   assert(demand > 0);
   assert(profile->loads[profile->ntimepoints-1] == 0);

   remainingduration = duration;
   startpos = pos;
   (*infeasible) = FALSE;

   if( profile->timepoints[startpos] > lst )
   {
      (*infeasible) = TRUE;
      return pos;
   }

   while( pos < profile->ntimepoints - 1 )
   {
      if( profile->loads[pos] + demand > profile->capacity )
      {
         SCIPdebugMessage("profile <%p>: core does not fit at time point %d (pos %d)\n", (void*)profile, profile->timepoints[pos], pos);
         startpos = pos + 1;
         remainingduration = duration;

         if( profile->timepoints[startpos] > lst )
         {
            (*infeasible) = TRUE;
            return pos;
         }
      }
      else
         remainingduration -= profile->timepoints[pos+1] - profile->timepoints[pos];

      if( remainingduration <= 0 )
         break;

      pos++;
   }

   return startpos;
}

/** return the earliest possible starting point within the time interval [lb,ub] for a given core (given by its demand
 *  and duration)
 */
int SCIPprofileGetEarliestFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   est,                /**< earliest starting time of the given core */
   int                   lst,                /**< latest starting time of the given core */
   int                   duration,           /**< duration of the core */
   int                   demand,             /**< demand of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the corer cannot be inserted */
   )
{
   SCIP_Bool found;
   int pos;

   assert(profile != NULL);
   assert(est >= 0);
   assert(est <= lst);
   assert(duration >= 0);
   assert(demand >= 0);
   assert(infeasible != NULL);
   assert(profile->ntimepoints > 0);
   assert(profile->loads[profile->ntimepoints-1] == 0);

   SCIPdebugMessage("profile <%p>: find earliest start time (demad %d, duration %d) [%d,%d]\n", (void*)profile, demand, duration, est, lst);

   if( duration == 0 || demand == 0 )
   {
      *infeasible = FALSE;
      return est;
   }

   found = SCIPprofileFindLeft(profile, est, &pos);
   SCIPdebugMessage("profile <%p>: earliest start time does %s exist as time point (pos %d)\n", (void*)profile, found ? "" : "not", pos);

   /* if the position is the last time point in the profile, the core can be inserted at its earliest start time */
   if( pos == profile->ntimepoints - 1 )
   {
      (*infeasible) = FALSE;
      return est;
   }

   if( found )
   {
      /* if the start time matches a time point in the profile we can just search */
      assert(profile->timepoints[pos] == est);
      pos = profileFindFeasibleStart(profile, pos, lst, duration, demand, infeasible);

      assert(pos < profile->ntimepoints);
      est = profile->timepoints[pos];
   }
   else if( profile->loads[pos] + demand > profile->capacity )
   {
      /* if the the time point left to the start time has not enough free capacity we can just search the profile
       * starting from the next time point
       */
      assert(profile->timepoints[pos] <= est);
      pos = profileFindFeasibleStart(profile, pos+1, lst, duration, demand, infeasible);

      assert(pos < profile->ntimepoints);
      est = profile->timepoints[pos];
   }
   else
   {
      int remainingduration;

      /* check if the core can be placed at its earliest start time */

      assert(pos < profile->ntimepoints - 1);

      remainingduration = duration - (profile->timepoints[pos+1] - est);
      SCIPdebugMessage("remaining duration %d\n", remainingduration);


      if( remainingduration <= 0 )
         (*infeasible) = FALSE;
      else
      {
         pos = profileFindFeasibleStart(profile, pos+1, profile->timepoints[pos+1], remainingduration, demand, infeasible);
         SCIPdebugMessage("remaining duration can%s be processed\n", *infeasible ? "not" : "");

         if( *infeasible )
         {
            pos = profileFindFeasibleStart(profile, pos+1, lst, duration, demand, infeasible);

            assert(pos < profile->ntimepoints);
            est = profile->timepoints[pos];
         }
      }
   }

   return est;
}

/** returns TRUE if the core (given by its demand and during) can be inserted at the given time point; otherwise FALSE */
static
int profileFindDownFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos,                /**< pointer to store the position in the profile to start the search */
   int                   ect,                /**< earliest completion time */
   int                   duration,           /**< duration of the core */
   int                   demand,             /**< demand of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the corer cannot be inserted */
   )
{
   int remainingduration;
   int endpos;

   assert(profile != NULL);
   assert(pos >= 0);
   assert(pos < profile->ntimepoints);
   assert(duration > 0);
   assert(demand > 0);
   assert(profile->ntimepoints > 0);
   assert(profile->loads[profile->ntimepoints-1] == 0);

   remainingduration = duration;
   endpos = pos;
   (*infeasible) = TRUE;

   if( profile->timepoints[endpos] < ect - duration )
      return pos;

   while( pos > 0 )
   {
      if( profile->loads[pos-1] + demand > profile->capacity )
      {
         SCIPdebugMessage("profile <%p>: core does not fit at time point %d (pos %d)\n", (void*)profile, profile->timepoints[pos-1], pos-1);

         endpos = pos - 1;
         remainingduration = duration;

         if( profile->timepoints[endpos] < ect - duration )
            return pos;
      }
      else
         remainingduration -= profile->timepoints[pos] - profile->timepoints[pos-1];

      if( remainingduration <= 0 )
      {
         *infeasible = FALSE;
         break;
      }

      pos--;
   }

   return endpos;
}

/** return the latest possible starting point within the time interval [lb,ub] for a given core (given by its demand and
 *  duration)
 */
int SCIPprofileGetLatestFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   est,                /**< earliest possible start point */
   int                   lst,                /**< latest possible start point */
   int                   duration,           /**< duration of the core */
   int                   demand,             /**< demand of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the core cannot be inserted */
   )
{
   SCIP_Bool found;
   int ect;
   int lct;
   int pos;

   assert(profile != NULL);
   assert(est >= 0);
   assert(est <= lst);
   assert(duration >= 0);
   assert(demand >= 0);
   assert(infeasible != NULL);
   assert(profile->ntimepoints > 0);
   assert(profile->loads[profile->ntimepoints-1] == 0);

   if( duration == 0 || demand == 0 )
   {
      *infeasible = FALSE;
      return lst;
   }

   ect = est + duration;
   lct = lst + duration;

   found = SCIPprofileFindLeft(profile, lct, &pos);
   SCIPdebugMessage("profile <%p>: latest completion time %d does %s exist as time point (pos %d)\n", (void*)profile, lct, found ? "" : "not", pos);

   if( found )
   {
      /* if the start time matches a time point in the profile we can just search */
      assert(profile->timepoints[pos] == lct);
      pos = profileFindDownFeasibleStart(profile, pos, ect, duration, demand, infeasible);

      assert(pos < profile->ntimepoints && pos >= 0);
      lct = profile->timepoints[pos];
   }
   else if( profile->loads[pos] + demand > profile->capacity )
   {
      /* if the time point left to the start time has not enough free capacity we can just search the profile starting
       * from the next time point
       */
      assert(profile->timepoints[pos] < lct);
      pos = profileFindDownFeasibleStart(profile, pos, ect, duration, demand, infeasible);

      assert(pos < profile->ntimepoints && pos >= 0);
      lct = profile->timepoints[pos];
   }
   else
   {
      int remainingduration;

      /* check if the core can be placed at its latest start time */
      assert(profile->timepoints[pos] < lct);

      remainingduration = duration - (lct - profile->timepoints[pos]);

      if( remainingduration <= 0 )
         (*infeasible) = FALSE;
      else
      {
         pos = profileFindDownFeasibleStart(profile, pos, profile->timepoints[pos], remainingduration, demand, infeasible);

         if( *infeasible )
         {
            pos = profileFindDownFeasibleStart(profile, pos, ect, duration, demand, infeasible);

            assert(pos < profile->ntimepoints && pos >= 0);
            lct = profile->timepoints[pos];
         }
      }
   }

   return lct - duration;
}

/*
 * Directed graph
 */

/** creates directed graph structure */
SCIP_RETCODE SCIPdigraphCreate(
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   BMS_BLKMEM*           blkmem,             /**< block memory to store the data */
   int                   nnodes              /**< number of nodes */
   )
{
   assert(digraph != NULL);
   assert(blkmem != NULL);
   assert(nnodes > 0);

   /* allocate memory for the graph and the arrays storing arcs and data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, digraph) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &(*digraph)->successors, nnodes) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &(*digraph)->arcdata, nnodes) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &(*digraph)->successorssize, nnodes) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &(*digraph)->nsuccessors, nnodes) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &(*digraph)->nodedata, nnodes) );

   /* store number of nodes */
   (*digraph)->nnodes = nnodes;

   /* at the beginning, no components are stored */
   (*digraph)->blkmem = blkmem;
   (*digraph)->ncomponents = 0;
   (*digraph)->componentstartsize = 0;
   (*digraph)->components = NULL;
   (*digraph)->componentstarts = NULL;

   return SCIP_OKAY;
}

/** resize directed graph structure */
SCIP_RETCODE SCIPdigraphResize(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   nnodes              /**< new number of nodes */
   )
{
   int n;
   assert(digraph != NULL);
   assert(digraph->blkmem != NULL);

   /* check if the digraph has already a proper size */
   if( nnodes <= digraph->nnodes )
      return SCIP_OKAY;

   /* reallocate memory for increasing the arrays storing arcs and data */
   SCIP_ALLOC( BMSreallocBlockMemoryArray(digraph->blkmem, &digraph->successors, digraph->nnodes, nnodes) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(digraph->blkmem, &digraph->arcdata, digraph->nnodes, nnodes) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(digraph->blkmem, &digraph->successorssize, digraph->nnodes, nnodes) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(digraph->blkmem, &digraph->nsuccessors, digraph->nnodes, nnodes) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(digraph->blkmem, &digraph->nodedata, digraph->nnodes, nnodes) );

   /* initialize the new node data structures */
   for( n = digraph->nnodes; n < nnodes; ++n )
   {
      digraph->nodedata[n] = NULL;
      digraph->arcdata[n] = NULL;
      digraph->successors[n] = NULL;
      digraph->successorssize[n] = 0;
      digraph->nsuccessors[n] = 0;
   }

   /* store the new number of nodes */
   digraph->nnodes = nnodes;

   return SCIP_OKAY;
}

/** copies directed graph structure
 *
 *  @note The data in nodedata is copied verbatim. This possibly has to be adapted by the user.
 */
SCIP_RETCODE SCIPdigraphCopy(
   SCIP_DIGRAPH**        targetdigraph,      /**< pointer to store the copied directed graph */
   SCIP_DIGRAPH*         sourcedigraph,      /**< source directed graph */
   BMS_BLKMEM*           targetblkmem        /**< block memory to store the target block memory, or NULL to use the same
                                              *   the same block memory as used for the \p sourcedigraph */
   )
{
   int ncomponents;
   int nnodes;
   int i;

   assert(sourcedigraph != NULL);
   assert(targetdigraph != NULL);

   /* use the source digraph block memory if not specified otherwise */
   if( targetblkmem == NULL )
      targetblkmem = sourcedigraph->blkmem;

   assert(targetblkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(targetblkmem, targetdigraph) );

   nnodes = sourcedigraph->nnodes;
   ncomponents = sourcedigraph->ncomponents;
   (*targetdigraph)->nnodes = nnodes;
   (*targetdigraph)->ncomponents = ncomponents;
   (*targetdigraph)->blkmem = targetblkmem;

   /* copy arcs and data */
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(targetblkmem, &(*targetdigraph)->successors, nnodes) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(targetblkmem, &(*targetdigraph)->arcdata, nnodes) );
   SCIP_ALLOC( BMSallocClearBlockMemoryArray(targetblkmem, &(*targetdigraph)->nodedata, nnodes) );

   /* copy lists of successors and arc data */
   for( i = 0; i < nnodes; ++i )
   {
      if( sourcedigraph->nsuccessors[i] > 0 )
      {
         assert(sourcedigraph->successors[i] != NULL);
         assert(sourcedigraph->arcdata[i] != NULL);

         SCIP_ALLOC( BMSduplicateBlockMemoryArray(targetblkmem, &((*targetdigraph)->successors[i]),
               sourcedigraph->successors[i], sourcedigraph->nsuccessors[i]) ); /*lint !e866*/
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(targetblkmem, &((*targetdigraph)->arcdata[i]),
               sourcedigraph->arcdata[i], sourcedigraph->nsuccessors[i]) ); /*lint !e866*/
      }
      /* copy node data - careful if these are pointers to some information -> need to be copied by hand */
      (*targetdigraph)->nodedata[i] = sourcedigraph->nodedata[i];
   }

   /* use nsuccessors as size to save memory */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(targetblkmem, &(*targetdigraph)->successorssize, sourcedigraph->nsuccessors, nnodes) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(targetblkmem, &(*targetdigraph)->nsuccessors, sourcedigraph->nsuccessors, nnodes) );

   /* copy component data */
   if( ncomponents > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(targetblkmem, &(*targetdigraph)->components, sourcedigraph->components,
            sourcedigraph->componentstarts[ncomponents]) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(targetblkmem, &(*targetdigraph)->componentstarts,
            sourcedigraph->componentstarts,ncomponents + 1) ); /*lint !e776*/
      (*targetdigraph)->componentstartsize = ncomponents + 1;
   }
   else
   {
      (*targetdigraph)->components = NULL;
      (*targetdigraph)->componentstarts = NULL;
      (*targetdigraph)->componentstartsize = 0;
   }

   return SCIP_OKAY;
}

/** sets the sizes of the successor lists for the nodes in a directed graph and allocates memory for the lists */
SCIP_RETCODE SCIPdigraphSetSizes(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int*                  sizes               /**< sizes of the successor lists */
   )
{
   int i;
   BMS_BLKMEM* blkmem;

   assert(digraph != NULL);
   assert(digraph->nnodes > 0);
   blkmem = digraph->blkmem;

   for( i = 0; i < digraph->nnodes; ++i )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &digraph->successors[i], sizes[i]) ); /*lint !e866*/
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &digraph->arcdata[i], sizes[i]) ); /*lint !e866*/
      digraph->successorssize[i] = sizes[i];
      digraph->nsuccessors[i] = 0;
   }

   return SCIP_OKAY;
}

/** frees given directed graph structure */
void SCIPdigraphFree(
   SCIP_DIGRAPH**        digraph             /**< pointer to the directed graph */
   )
{
   int i;
   BMS_BLKMEM* blkmem;
   SCIP_DIGRAPH* digraphptr;

   assert(digraph != NULL);
   assert(*digraph != NULL);
   assert((*digraph)->blkmem != NULL);

   blkmem = (*digraph)->blkmem;
   digraphptr = *digraph;

   /* free arrays storing the successor nodes and arc data */
   for( i = digraphptr->nnodes - 1; i >= 0; --i )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &digraphptr->successors[i], digraphptr->successorssize[i]);
      BMSfreeBlockMemoryArrayNull(blkmem, &digraphptr->arcdata[i], digraphptr->successorssize[i]);
   }

   /* free components structure */
   SCIPdigraphFreeComponents(digraphptr);
   assert(digraphptr->ncomponents == 0);
   assert(digraphptr->componentstartsize == 0);
   assert(digraphptr->components == NULL);
   assert(digraphptr->componentstarts == NULL);

   /* free directed graph data structure */
   BMSfreeBlockMemoryArray(blkmem, &digraphptr->nodedata, digraphptr->nnodes);
   BMSfreeBlockMemoryArray(blkmem, &digraphptr->successorssize, digraphptr->nnodes);
   BMSfreeBlockMemoryArray(blkmem, &digraphptr->nsuccessors, digraphptr->nnodes);
   BMSfreeBlockMemoryArray(blkmem, &digraphptr->successors, digraphptr->nnodes);
   BMSfreeBlockMemoryArray(blkmem, &digraphptr->arcdata, digraphptr->nnodes);

   BMSfreeBlockMemory(blkmem, digraph);
}

#define STARTSUCCESSORSSIZE 5

/** ensures that successors array of one node in a directed graph is big enough */
static
SCIP_RETCODE ensureSuccessorsSize(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   idx,                /**< index for which the size is ensured */
   int                   newsize             /**< needed size */
   )
{
   BMS_BLKMEM* blkmem;

   assert(digraph != NULL);
   assert(digraph->blkmem != NULL);
   assert(idx >= 0);
   assert(idx < digraph->nnodes);
   assert(newsize > 0);
   assert(digraph->successorssize[idx] == 0 || digraph->successors[idx] != NULL);
   assert(digraph->successorssize[idx] == 0 || digraph->arcdata[idx] != NULL);

   blkmem = digraph->blkmem;

   /* check whether array is big enough, and realloc, if needed */
   if( newsize > digraph->successorssize[idx] )
   {
      if( digraph->successors[idx] == NULL )
      {
         assert(digraph->arcdata[idx] == NULL);
         digraph->successorssize[idx] = STARTSUCCESSORSSIZE;
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &digraph->successors[idx], digraph->successorssize[idx]) ); /*lint !e866*/
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &digraph->arcdata[idx], digraph->successorssize[idx]) ); /*lint !e866*/
      }
      else
      {
         newsize = MAX(newsize, 2 * digraph->successorssize[idx]);
         assert(digraph->arcdata[idx] != NULL);
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &digraph->successors[idx], digraph->successorssize[idx], newsize) ); /*lint !e866*/
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &digraph->arcdata[idx], digraph->successorssize[idx], newsize) ); /*lint !e866*/
         digraph->successorssize[idx] = newsize;
      }
   }

   assert(newsize <= digraph->successorssize[idx]);

   return SCIP_OKAY;
}

/** add (directed) arc and a related data to the directed graph structure
 *
 *  @note if the arc is already contained, it is added a second time
 */
SCIP_RETCODE SCIPdigraphAddArc(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   )
{
   assert(digraph != NULL);
   assert(startnode >= 0);
   assert(endnode >= 0);
   assert(startnode < digraph->nnodes);
   assert(endnode < digraph->nnodes);

   SCIP_CALL( ensureSuccessorsSize(digraph, startnode, digraph->nsuccessors[startnode] + 1) );

   /* add arc */
   digraph->successors[startnode][digraph->nsuccessors[startnode]] = endnode;
   digraph->arcdata[startnode][digraph->nsuccessors[startnode]] = data;
   digraph->nsuccessors[startnode]++;

   return SCIP_OKAY;
}

/** add (directed) arc to the directed graph structure, if it is not contained, yet
 *
 * @note if there already exists an arc from startnode to endnode, the new arc is not added,
 *       even if its data is different
 */
SCIP_RETCODE SCIPdigraphAddArcSafe(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   )
{
   int nsuccessors;
   int i;

   assert(digraph != NULL);
   assert(startnode >= 0);
   assert(endnode >= 0);
   assert(startnode < digraph->nnodes);
   assert(endnode < digraph->nnodes);

   nsuccessors = digraph->nsuccessors[startnode];

   /* search for the arc in existing arcs */
   for( i = 0; i < nsuccessors; ++i )
      if( digraph->successors[startnode][i] == endnode )
         return SCIP_OKAY;

   SCIP_CALL( ensureSuccessorsSize(digraph, startnode, nsuccessors + 1) );

   /* add arc */
   digraph->successors[startnode][nsuccessors] = endnode;
   digraph->arcdata[startnode][nsuccessors] = data;
   ++(digraph->nsuccessors[startnode]);

   return SCIP_OKAY;
}

/** sets the number of successors to a given value */
SCIP_RETCODE SCIPdigraphSetNSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node,               /**< node for which the number of successors has to be changed */
   int                   nsuccessors         /**< new number of successors */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);

   digraph->nsuccessors[node] = nsuccessors;

   return SCIP_OKAY;
}

/** returns the number of nodes of the given digraph */
int SCIPdigraphGetNNodes(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   assert(digraph != NULL);

   return digraph->nnodes;
}

/** returns the node data, or NULL if no data exist */
void* SCIPdigraphGetNodeData(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the node data is returned */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);

   return digraph->nodedata[node];
}

/** sets the node data
 *
 *  @note The old user pointer is not freed. This has to be done by the user
 */
void SCIPdigraphSetNodeData(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   void*                 dataptr,            /**< user node data pointer, or NULL */
   int                   node                /**< node for which the node data is returned */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);

   digraph->nodedata[node] = dataptr;
}

/** returns the total number of arcs in the given digraph */
int SCIPdigraphGetNArcs(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   int i;
   int narcs;

   assert(digraph != NULL);

   /* count number of arcs */
   narcs = 0;
   for( i = 0; i < digraph->nnodes; ++i )
      narcs += digraph->nsuccessors[i];

   return narcs;
}

/** returns the number of successor nodes of the given node */
int SCIPdigraphGetNSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the number of outgoing arcs is returned */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);
   assert(digraph->nsuccessors[node] >= 0);
   assert(digraph->nsuccessors[node] <= digraph->successorssize[node]);

   return digraph->nsuccessors[node];
}

/** returns the array of indices of the successor nodes; this array must not be changed from outside */
int* SCIPdigraphGetSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the array of outgoing arcs is returned */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);
   assert(digraph->nsuccessors[node] >= 0);
   assert(digraph->nsuccessors[node] <= digraph->successorssize[node]);
   assert((digraph->nsuccessors[node] == 0) || (digraph->successors[node] != NULL));

   return digraph->successors[node];
}

/** returns the array of data corresponding to the arcs originating at the given node, or NULL if no data exist; this
 *  array must not be changed from outside
 */
void** SCIPdigraphGetSuccessorsData(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the data corresponding to the outgoing arcs is returned */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);
   assert(digraph->nsuccessors[node] >= 0);
   assert(digraph->nsuccessors[node] <= digraph->successorssize[node]);
   assert(digraph->arcdata != NULL);

   return digraph->arcdata[node];
}

/** performs depth-first-search in the given directed graph from the given start node */
static
void depthFirstSearch(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< node to start the depth-first-search */
   SCIP_Bool*            visited,            /**< array to store for each node, whether it was already visited */
   int*                  dfsstack,           /**< array of size number of nodes to store the stack;
                                              *   only needed for performance reasons */
   int*                  stackadjvisited,    /**< array of size number of nodes to store the number of adjacent nodes already visited
                                              *   for each node on the stack; only needed for performance reasons */
   int*                  dfsnodes,           /**< array of nodes that can be reached starting at startnode, in reverse dfs order */
   int*                  ndfsnodes           /**< pointer to store number of nodes that can be reached starting at startnode */
   )
{
   int stackidx;

   assert(digraph != NULL);
   assert(startnode >= 0);
   assert(startnode < digraph->nnodes);
   assert(visited != NULL);
   assert(visited[startnode] == FALSE);
   assert(dfsstack != NULL);
   assert(dfsnodes != NULL);
   assert(ndfsnodes != NULL);

   /* put start node on the stack */
   dfsstack[0] = startnode;
   stackadjvisited[0] = 0;
   stackidx = 0;

   while( stackidx >= 0 )
   {
      int currnode;
      int sadv;

      /* get next node from stack */
      currnode = dfsstack[stackidx];

      sadv = stackadjvisited[stackidx];
      assert( 0 <= sadv && sadv <= digraph->nsuccessors[currnode] );

      /* mark current node as visited */
      assert( visited[currnode] == (sadv > 0) );
      visited[currnode] = TRUE;

      /* iterate through the successor list until we reach unhandled node */
      while( sadv < digraph->nsuccessors[currnode] && visited[digraph->successors[currnode][sadv]] )
         ++sadv;

      /* the current node was completely handled, remove it from stack */
      if( sadv == digraph->nsuccessors[currnode] )
      {
         --stackidx;

         /* store node in the sorted nodes array */
         dfsnodes[(*ndfsnodes)++] = currnode;
      }
      /* handle next unhandled successor node */
      else
      {
         assert( ! visited[digraph->successors[currnode][sadv]] );

         /* store current stackadjvisted index */
         stackadjvisited[stackidx] = sadv + 1;

         /* put the successor node onto the stack */
         ++stackidx;
         dfsstack[stackidx] = digraph->successors[currnode][sadv];
         stackadjvisited[stackidx] = 0;
         assert( stackidx < digraph->nnodes );
      }
   }
}

/** Compute undirected connected components on the given graph.
 *
 *  @note For each arc, its reverse is added, so the graph does not need to be the directed representation of an
 *        undirected graph.
 */
SCIP_RETCODE SCIPdigraphComputeUndirectedComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   minsize,            /**< all components with less nodes are ignored */
   int*                  components,         /**< array with as many slots as there are nodes in the directed graph
                                              *   to store for each node the component to which it belongs
                                              *   (components are numbered 0 to ncomponents - 1); or NULL, if components
                                              *   are accessed one-by-one using SCIPdigraphGetComponent() */
   int*                  ncomponents         /**< pointer to store the number of components; or NULL, if the
                                              *   number of components is accessed by SCIPdigraphGetNComponents() */
   )
{
   BMS_BLKMEM* blkmem;
   SCIP_Bool* visited;
   int* ndirectedsuccessors;
   int* stackadjvisited;
   int* dfsstack;
   int ndfsnodes;
   int compstart;
   int v;
   int i;
   int j;

   SCIP_RETCODE retcode = SCIP_OKAY;

   assert(digraph != NULL);
   assert(digraph->nnodes > 0);
   assert(digraph->blkmem != NULL);

   blkmem = digraph->blkmem;

   /* first free the old components */
   if( digraph->ncomponents > 0 )
   {
      SCIPdigraphFreeComponents(digraph);
   }

   digraph->ncomponents = 0;
   digraph->componentstartsize = 10;

   /* storage to hold components is stored in block memory */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &digraph->components, digraph->nnodes) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &digraph->componentstarts, digraph->componentstartsize) );

   /* allocate temporary arrays */
   SCIP_ALLOC_TERMINATE( retcode, BMSallocClearMemoryArray(&visited, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&dfsstack, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&stackadjvisited, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&ndirectedsuccessors, digraph->nnodes), TERMINATE );

   digraph->componentstarts[0] = 0;

   /* store the number of directed arcs per node */
   BMScopyMemoryArray(ndirectedsuccessors, digraph->nsuccessors, digraph->nnodes);

   /* add reverse arcs to the graph */
   for( i = digraph->nnodes - 1; i >= 0; --i )
   {
      for( j = 0; j < ndirectedsuccessors[i]; ++j )
      {
         SCIP_CALL_TERMINATE( retcode, SCIPdigraphAddArc(digraph, digraph->successors[i][j], i, NULL), TERMINATE );
      }
   }

   for( v = 0; v < digraph->nnodes; ++v )
   {
      if( visited[v] )
         continue;

      compstart = digraph->componentstarts[digraph->ncomponents];
      ndfsnodes = 0;
      depthFirstSearch(digraph, v, visited, dfsstack, stackadjvisited,
         &digraph->components[compstart], &ndfsnodes);

      /* forget about this component if it is too small */
      if( ndfsnodes >= minsize )
      {
         digraph->ncomponents++;

         /* enlarge componentstartsize array, if needed */
         if( digraph->ncomponents >= digraph->componentstartsize )
         {
            int newsize;

            newsize = 2 * digraph->componentstartsize;
            assert(digraph->ncomponents < newsize);

            SCIP_ALLOC_TERMINATE( retcode, BMSreallocBlockMemoryArray(blkmem, &digraph->componentstarts, digraph->componentstartsize, newsize), TERMINATE );
            digraph->componentstartsize = newsize;
         }
         digraph->componentstarts[digraph->ncomponents] = compstart + ndfsnodes;

         /* store component number for contained nodes if array was given */
         if( components != NULL )
         {
            for( i = digraph->componentstarts[digraph->ncomponents] - 1; i >= compstart; --i )
            {
               components[digraph->components[i]] = digraph->ncomponents - 1;
            }
         }
      }
   }

   /* restore the number of directed arcs per node */
   BMScopyMemoryArray(digraph->nsuccessors, ndirectedsuccessors, digraph->nnodes);
   BMSclearMemoryArray(visited, digraph->nnodes);

   /* return number of components, if the pointer was given */
   if( ncomponents != NULL )
      (*ncomponents) = digraph->ncomponents;

TERMINATE:
   if( retcode != SCIP_OKAY )
   {
      SCIPdigraphFreeComponents(digraph);
   }
   BMSfreeMemoryArrayNull(&ndirectedsuccessors);
   BMSfreeMemoryArrayNull(&stackadjvisited);
   BMSfreeMemoryArrayNull(&dfsstack);
   BMSfreeMemoryArrayNull(&visited);

   return retcode;
}

/** Performs an (almost) topological sort on the undirected components of the given directed graph. The undirected
 *  components should be computed before using SCIPdigraphComputeUndirectedComponents().
 *
 *  @note In general a topological sort is not unique.  Note, that there might be directed cycles, that are randomly
 *        broken, which is the reason for having only almost topologically sorted arrays.
 */
SCIP_RETCODE SCIPdigraphTopoSortComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   SCIP_Bool* visited = NULL;
   int* comps;
   int* compstarts;
   int* stackadjvisited = NULL;
   int* dfsstack = NULL;
   int* dfsnodes = NULL;
   int ndfsnodes;
   int ncomps;
   int i;
   int j;
   int k;
   int endidx;
   SCIP_RETCODE retcode = SCIP_OKAY;

   assert(digraph != NULL);

   ncomps = digraph->ncomponents;
   comps = digraph->components;
   compstarts = digraph->componentstarts;

   SCIP_ALLOC_TERMINATE( retcode, BMSallocClearMemoryArray(&visited, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&dfsnodes, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&dfsstack, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&stackadjvisited, digraph->nnodes), TERMINATE );

   /* sort the components (almost) topologically */
   for( i = 0; i < ncomps; ++i )
   {
      endidx = compstarts[i+1] - 1;
      ndfsnodes = 0;
      for( j = compstarts[i]; j < compstarts[i+1]; ++j )
      {
         if( visited[comps[j]] )
            continue;

         /* perform depth first search, nodes visited in this call are appended to the list dfsnodes in reverse
          * dfs order, after the nodes already contained;
          * so at every point in time, the nodes in dfsnode are in reverse (almost) topological order
          */
         depthFirstSearch(digraph, comps[j], visited, dfsstack, stackadjvisited, dfsnodes, &ndfsnodes);
      }
      assert(endidx - ndfsnodes == compstarts[i] - 1);

      /* copy reverse (almost) topologically sorted array of nodes reached by the dfs searches;
       * reverse their order to get an (almost) topologically sort
       */
      for( k = 0; k < ndfsnodes; ++k )
      {
         digraph->components[endidx - k] = dfsnodes[k];
      }
   }

TERMINATE:
   BMSfreeMemoryArrayNull(&stackadjvisited);
   BMSfreeMemoryArrayNull(&dfsstack);
   BMSfreeMemoryArrayNull(&dfsnodes);
   BMSfreeMemoryArrayNull(&visited);

   return retcode;
}

/** returns the number of previously computed undirected components for the given directed graph */
int SCIPdigraphGetNComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   assert(digraph != NULL);
   assert(digraph->componentstartsize > 0); /* components should have been computed */

   return digraph->ncomponents;
}

/** Returns the previously computed undirected component of the given number for the given directed graph.
 *  If the components were sorted using SCIPdigraphTopoSortComponents(), the component is (almost) topologically sorted.
 */
void SCIPdigraphGetComponent(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   compidx,            /**< number of the component to return */
   int**                 nodes,              /**< pointer to store the nodes in the component; or NULL, if not needed */
   int*                  nnodes              /**< pointer to store the number of nodes in the component;
                                              *   or NULL, if not needed */
   )
{
   assert(digraph != NULL);
   assert(compidx >= 0);
   assert(compidx < digraph->ncomponents);
   assert(nodes != NULL || nnodes != NULL);

   if( nodes != NULL )
      (*nodes) = &(digraph->components[digraph->componentstarts[compidx]]);
   if( nnodes != NULL )
      (*nnodes) = digraph->componentstarts[compidx + 1] - digraph->componentstarts[compidx];
}

/* Performs Tarjan's algorithm for a given directed graph to obtain the strongly connected components
 * which are reachable from a given node.
 */
static
void tarjan(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   v,                  /**< node to start the algorithm */
   int*                  lowlink,            /**< array to store lowlink values */
   int*                  dfsidx,             /**< array to store dfs indices */
   int*                  stack,              /**< array to store a stack */
   int*                  stacksize,          /**< pointer to store the size of the stack */
   SCIP_Bool*            unprocessed,        /**< array to store which node is unprocessed yet */
   SCIP_Bool*            nodeinstack,        /**< array to store which nodes are in the stack */
   int*                  maxdfs,             /**< pointer to store index for DFS */
   int*                  strongcomponents,   /**< array to store for each node the strongly connected
                                              *   component to which it belongs (components are
                                              *   numbered 0 to nstrongcomponents - 1); */
   int*                  nstrongcomponents,  /**< pointer to store the number of computed components so far */
   int*                  strongcompstartidx, /**< array to store the start index of the computed components */
   int*                  nstorednodes        /**< pointer to store the number of already stored nodes */
   )
{
   int i;

   assert(digraph != NULL);
   assert(v >= 0);
   assert(v < digraph->nnodes);
   assert(lowlink != NULL);
   assert(dfsidx != NULL);
   assert(stack != NULL);
   assert(stacksize != NULL);
   assert(*stacksize >= 0);
   assert(*stacksize < digraph->nnodes);
   assert(unprocessed != NULL);
   assert(nodeinstack != NULL);
   assert(maxdfs != NULL);
   assert(strongcomponents != NULL);
   assert(nstrongcomponents != NULL);
   assert(strongcompstartidx != NULL);
   assert(nstorednodes != NULL);
   assert(*nstorednodes >= 0 && *nstorednodes < digraph->nnodes);

   dfsidx[v] = *maxdfs;
   lowlink[v] = *maxdfs;
   *maxdfs += 1;

   /* add v to the stack */
   stack[*stacksize] = v;
   *stacksize += 1;
   nodeinstack[v] = TRUE;

   /* mark v as processed */
   unprocessed[v] = FALSE;

   for( i = 0; i < digraph->nsuccessors[v]; ++i )
   {
      int w;

      /* edge (v,w) */
      w = digraph->successors[v][i];

      if( unprocessed[w] )
      {
         tarjan(digraph, w, lowlink, dfsidx, stack, stacksize, unprocessed, nodeinstack, maxdfs, strongcomponents,
               nstrongcomponents, strongcompstartidx, nstorednodes);

         assert(lowlink[v] >= 0 && lowlink[v] < digraph->nnodes);
         assert(lowlink[w] >= 0 && lowlink[w] < digraph->nnodes);

         /* update lowlink */
         lowlink[v] = MIN(lowlink[v], lowlink[w]);
      }
      else if( nodeinstack[w] )
      {
         assert(lowlink[v] >= 0 && lowlink[v] < digraph->nnodes);
         assert(dfsidx[w] >= 0 && dfsidx[w] < digraph->nnodes);

         /* update lowlink */
         lowlink[v] = MIN(lowlink[v], dfsidx[w]);
      }
   }

   /* found a root of a strong component */
   if( lowlink[v] == dfsidx[v] )
   {
      int w;

      strongcompstartidx[*nstrongcomponents] = *nstorednodes;
      *nstrongcomponents += 1;

      do
      {
         assert(*stacksize > 0);

         /* stack.pop() */
         w = stack[*stacksize - 1];
         *stacksize -= 1;
         nodeinstack[w] = FALSE;

         /* store the node in the corresponding component */
         strongcomponents[*nstorednodes] = w;
         *nstorednodes += 1;
      }
      while( v != w );
   }
}

/** Computes all strongly connected components of an undirected connected component with Tarjan's Algorithm.
 *  The resulting strongly connected components are sorted topologically (starting from the end of the
 *  strongcomponents array).
 *
 *  @note In general a topological sort of the strongly connected components is not unique.
 */
SCIP_RETCODE SCIPdigraphComputeDirectedComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   compidx,            /**< number of the undirected connected component */
   int*                  strongcomponents,   /**< array to store the strongly connected components
                                              *   (length >= size of the component) */
   int*                  strongcompstartidx, /**< array to store the start indices of the strongly connected
                                              *   components (length >= size of the component) */
   int*                  nstrongcomponents   /**< pointer to store the number of strongly connected
                                              *   components */
   )
{
   int* lowlink;
   int* dfsidx;
   int* stack;
   int stacksize;
   SCIP_Bool* unprocessed;
   SCIP_Bool* nodeinstack;
   int maxdfs;
   int nstorednodes;
   int i;
   SCIP_RETCODE retcode;

   assert(digraph != NULL);
   assert(compidx >= 0);
   assert(compidx < digraph->ncomponents);
   assert(strongcomponents != NULL);
   assert(strongcompstartidx != NULL);
   assert(nstrongcomponents != NULL);

   retcode = SCIP_OKAY;

   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&lowlink, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&dfsidx, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&stack, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&unprocessed, digraph->nnodes), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&nodeinstack, digraph->nnodes), TERMINATE );

   for( i = 0; i < digraph->nnodes; ++i )
   {
      lowlink[i] = -1;
      dfsidx[i] = -1;
      stack[i] = -1;
      unprocessed[i] = TRUE;
      nodeinstack[i] = FALSE;
   }

   nstorednodes = 0;
   stacksize = 0;
   maxdfs = 0;
   *nstrongcomponents = 0;

   /* iterate over all nodes in the undirected connected component */
   for( i = digraph->componentstarts[compidx]; i < digraph->componentstarts[compidx + 1]; ++i )
   {
      int v;

      v = digraph->components[i];
      assert(v >= 0 && v < digraph->nnodes);

      /* call Tarjan's algorithm for unprocessed nodes */
      if( unprocessed[v] )
      {
         SCIPdebugMessage("apply Tarjan's algorithm for node %d\n", v);
         tarjan(digraph, v, lowlink, dfsidx, stack, &stacksize, unprocessed, nodeinstack, &maxdfs,
               strongcomponents, nstrongcomponents, strongcompstartidx, &nstorednodes);
      }
   }

   /* we should have stored as many nodes as in the undirected connected component */
   assert(nstorednodes == digraph->componentstarts[compidx + 1] - digraph->componentstarts[compidx]);

   /* to simplify the iteration over all strongly connected components */
   assert(*nstrongcomponents < digraph->nnodes + 1);
   strongcompstartidx[*nstrongcomponents] = nstorednodes;

   assert(retcode == SCIP_OKAY);

TERMINATE:
   BMSfreeMemoryArrayNull(&lowlink);
   BMSfreeMemoryArrayNull(&dfsidx);
   BMSfreeMemoryArrayNull(&stack);
   BMSfreeMemoryArrayNull(&unprocessed);
   BMSfreeMemoryArrayNull(&nodeinstack);

   return retcode;
}

/** frees the component information for the given directed graph */
void SCIPdigraphFreeComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   BMS_BLKMEM* blkmem;

   assert(digraph != NULL);
   assert(digraph->blkmem != NULL);

   blkmem = digraph->blkmem;

   /* free components structure */
   if( digraph->componentstartsize > 0 )
   {
      BMSfreeBlockMemoryArray(blkmem, &digraph->componentstarts, digraph->componentstartsize);
      BMSfreeBlockMemoryArray(blkmem, &digraph->components, digraph->nnodes);
      digraph->components = NULL;
      digraph->componentstarts = NULL;
      digraph->ncomponents = 0;
      digraph->componentstartsize = 0;
   }
#ifndef NDEBUG
   else
   {
      assert(digraph->components == NULL);
      assert(digraph->componentstarts == NULL);
      assert(digraph->ncomponents == 0);
   }
#endif
}

/** output of the given directed graph via the given message handler */
void SCIPdigraphPrint(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int n;

   for( n = 0; n < digraph->nnodes; ++n )
   {
      int* successors;
      int nsuccessors;
      int m;

      nsuccessors = digraph->nsuccessors[n];
      successors = digraph->successors[n];

      SCIPmessageFPrintInfo(messagehdlr, file, "node %d --> ", n);

      for( m = 0; m < nsuccessors ; ++m )
      {
         if( m == 0 )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%d", successors[m]);
         }
         else
         {
            SCIPmessageFPrintInfo(messagehdlr, file, ", %d", successors[m]);
         }
      }
      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }
}

/** prints the given directed graph structure in GML format into the given file */
void SCIPdigraphPrintGml(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   FILE*                 file                /**< file to write to */
   )
{
   int n;

   /* write GML format opening */
   SCIPgmlWriteOpening(file, TRUE);

   /* write all nodes of the graph */
   for( n = 0; n < digraph->nnodes; ++n )
   {
      char label[SCIP_MAXSTRLEN];

      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%d", n);
      SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", NULL, NULL);
   }

   /* write all edges */
   for( n = 0; n < digraph->nnodes; ++n )
   {
      int* successors;
      int nsuccessors;
      int m;

      nsuccessors = digraph->nsuccessors[n];
      successors = digraph->successors[n];

      for( m = 0; m < nsuccessors; ++m )
      {
         SCIPgmlWriteArc(file, (unsigned int)n, (unsigned int)successors[m], NULL, NULL);
      }
   }
   /* write GML format closing */
   SCIPgmlWriteClosing(file);
}

/** output of the given directed graph via the given message handler */
void SCIPdigraphPrintComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int c;
   int i;

   for( c = 0; c < digraph->ncomponents; ++c )
   {
      int start = digraph->componentstarts[c];
      int end =  digraph->componentstarts[c+1];

      SCIPmessageFPrintInfo(messagehdlr, file, "Components %d --> ", c);

      for( i = start; i < end; ++i )
      {
         if( i == start )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%d", digraph->components[i]);
         }
         else
         {
            SCIPmessageFPrintInfo(messagehdlr, file, ", %d", digraph->components[i]);
         }
      }
      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }
}

/*
 * Binary tree
 */

/** creates a node for a binary tree */
static
SCIP_RETCODE btnodeCreateEmpty(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node                /**< pointer to store the created node */
   )
{
   SCIP_ALLOC( BMSallocBlockMemory(tree->blkmem, node) );

   (*node)->parent = NULL;
   (*node)->left = NULL;
   (*node)->right = NULL;
   (*node)->dataptr = NULL;

   return SCIP_OKAY;
}

/** creates a tree node with (optinal) user data */
SCIP_RETCODE SCIPbtnodeCreate(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node,               /**< pointer to store the created node */
   void*                 dataptr             /**< user node data pointer, or NULL */
   )
{
   assert(tree != NULL);
   assert(node != NULL);

   SCIP_CALL( btnodeCreateEmpty(tree, node) );

   assert((*node)->parent == NULL);
   assert((*node)->left == NULL);
   assert((*node)->right == NULL);

   /* initialize user data */
   (*node)->dataptr = dataptr;

   return SCIP_OKAY;
}

/** frees a tree leaf */
static
void btnodeFreeLeaf(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node                /**< pointer to node which has to be freed */
   )
{
   assert(tree != NULL);
   assert(node != NULL);
   assert(*node != NULL);

   assert((*node)->left == NULL);
   assert((*node)->right == NULL);

#if 0
   /* remove reference from parent node */
   if( (*node)->parent != NULL )
   {
      assert(*node != NULL);

      assert((*node)->parent->left == *node || ((*node)->parent->right == *node));

      if( (*node)->parent->left == *node )
      {
         (*node)->parent->left = NULL;
      }
      else
      {
         assert((*node)->parent->right == *node);
         (*node)->parent->right = NULL;
      }
   }
#endif

   assert(*node != NULL);
   BMSfreeBlockMemory(tree->blkmem, node);
   assert(*node == NULL);
}

/** frees the node including the rooted subtree
 *
 *  @note The user pointer (object) is not freed. If needed, it has to be done by the user.
 */
void SCIPbtnodeFree(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node                /**< node to be freed */
   )
{
   assert(tree != NULL);
   assert(node != NULL);
   assert(*node != NULL);

   if( (*node)->left != NULL )
   {
      SCIPbtnodeFree(tree, &(*node)->left);
      assert((*node)->left == NULL);
   }

   if( (*node)->right != NULL )
   {
      SCIPbtnodeFree(tree, &(*node)->right);
      assert((*node)->right == NULL);
   }

   btnodeFreeLeaf(tree, node);
   assert(*node == NULL);
}

/* some simple variable functions implemented as defines */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPbtnodeGetData
#undef SCIPbtnodeGetKey
#undef SCIPbtnodeGetParent
#undef SCIPbtnodeGetLeftchild
#undef SCIPbtnodeGetRightchild
#undef SCIPbtnodeGetSibling
#undef SCIPbtnodeIsRoot
#undef SCIPbtnodeIsLeaf
#undef SCIPbtnodeIsLeftchild
#undef SCIPbtnodeIsRightchild

/** returns the user data pointer stored in that node */
void* SCIPbtnodeGetData(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return node->dataptr;
}

/** returns the parent which can be NULL if the given node is the root */
SCIP_BTNODE* SCIPbtnodeGetParent(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return node->parent;
}

/** returns left child which can be NULL if the given node is a leaf */
SCIP_BTNODE* SCIPbtnodeGetLeftchild(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return node->left;
}

/** returns right child which can be NULL if the given node is a leaf */
SCIP_BTNODE* SCIPbtnodeGetRightchild(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return node->right;
}

/** returns the sibling of the node or NULL if does not exist */
SCIP_BTNODE* SCIPbtnodeGetSibling(
   SCIP_BTNODE*          node                /**< node */
   )
{
   SCIP_BTNODE* parent;

   parent = SCIPbtnodeGetParent(node);

   if( parent == NULL )
      return NULL;

   if( SCIPbtnodeGetLeftchild(parent) == node )
      return SCIPbtnodeGetRightchild(parent);

   assert(SCIPbtnodeGetRightchild(parent) == node);

   return SCIPbtnodeGetLeftchild(parent);
}

/** returns whether the node is a root node */
SCIP_Bool SCIPbtnodeIsRoot(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return (node->parent == NULL);
}

/** returns whether the node is a leaf */
SCIP_Bool SCIPbtnodeIsLeaf(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return (node->left == NULL && node->right == NULL);
}

/** returns TRUE if the given node is left child */
SCIP_Bool SCIPbtnodeIsLeftchild(
   SCIP_BTNODE*          node                /**< node */
   )
{
   SCIP_BTNODE* parent;

   if( SCIPbtnodeIsRoot(node) )
      return FALSE;

   parent = SCIPbtnodeGetParent(node);

   if( SCIPbtnodeGetLeftchild(parent) == node )
      return TRUE;

   return FALSE;
}

/** returns TRUE if the given node is right child */
SCIP_Bool SCIPbtnodeIsRightchild(
   SCIP_BTNODE*          node                /**< node */
   )
{
   SCIP_BTNODE* parent;

   if( SCIPbtnodeIsRoot(node) )
      return FALSE;

   parent = SCIPbtnodeGetParent(node);

   if( SCIPbtnodeGetRightchild(parent) == node )
      return TRUE;

   return FALSE;
}

/** sets the give node data
 *
 *  @note The old user pointer is not freed.
 */
void SCIPbtnodeSetData(
   SCIP_BTNODE*          node,               /**< node */
   void*                 dataptr             /**< node user data pointer */
   )
{
   assert(node != NULL);

   node->dataptr = dataptr;
}

/** sets parent node
 *
 *  @note The old parent including the rooted subtree is not delete.
 */
void SCIPbtnodeSetParent(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          parent              /**< new parent node, or NULL */
   )
{
   assert(node != NULL);

   node->parent = parent;
}

/** sets left child
 *
 *  @note The old left child including the rooted subtree is not delete.
 */
void SCIPbtnodeSetLeftchild(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          left                /**< new left child, or NULL */
   )
{
   assert(node != NULL);

   node->left = left;
}

/** sets right child
 *
 *  @note The old right child including the rooted subtree is not delete.
 */
void SCIPbtnodeSetRightchild(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          right               /**< new right child, or NULL */
   )
{
   assert(node != NULL);

   node->right = right;
}

/** creates an binary tree */
SCIP_RETCODE SCIPbtCreate(
   SCIP_BT**             tree,               /**< pointer to store the created binary tree */
   BMS_BLKMEM*           blkmem              /**< block memory used to createnode */
   )
{
   assert(tree != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, tree) );
   (*tree)->blkmem = blkmem;
   (*tree)->root = NULL;

   return SCIP_OKAY;
}

/** frees binary tree
 *
 *  @note The user pointers (object) of the nodes are not freed. If needed, it has to be done by the user.
 */
void SCIPbtFree(
   SCIP_BT**             tree                /**< pointer to binary tree */
   )
{
   assert(tree != NULL);

   if( (*tree)->root != NULL )
   {
      SCIPbtnodeFree(*tree, &((*tree)->root));
   }

   BMSfreeBlockMemory((*tree)->blkmem, tree);
}

/** prints the rooted subtree of the given binary tree node in GML format into the given file */
static
void btPrintSubtree(
   SCIP_BTNODE*          node,               /**< binary tree node */
   FILE*                 file,               /**< file to write to */
   int*                  nnodes              /**< pointer to count the number of nodes */
   )
{
   SCIP_BTNODE* left;
   SCIP_BTNODE* right;
   char label[SCIP_MAXSTRLEN];

   assert(node != NULL);

   (*nnodes)++;
   (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%d", *nnodes);

   SCIPgmlWriteNode(file, (unsigned int)(size_t)node, label, "circle", NULL, NULL);

   left = SCIPbtnodeGetLeftchild(node);
   right = SCIPbtnodeGetRightchild(node);

   if( left != NULL )
   {
      btPrintSubtree(left, file, nnodes);

      SCIPgmlWriteArc(file, (unsigned int)(size_t)node, (unsigned int)(size_t)left, NULL, NULL);
   }

   if( right != NULL )
   {
      btPrintSubtree(right, file, nnodes);

      SCIPgmlWriteArc(file, (unsigned int)(size_t)node, (unsigned int)(size_t)right, NULL, NULL);
   }
}

/** prints the binary tree in GML format into the given file */
void SCIPbtPrintGml(
   SCIP_BT*              tree,               /**< binary tree */
   FILE*                 file                /**< file to write to */
   )
{
   /* write GML opening */
   SCIPgmlWriteOpening(file, TRUE);

   if( !SCIPbtIsEmpty(tree) )
   {
      SCIP_BTNODE* root;
      int nnodes;

      root = SCIPbtGetRoot(tree);
      assert(root != NULL);

      nnodes = 0;

      btPrintSubtree(root, file, &nnodes);
   }

   /* write GML closing */
   SCIPgmlWriteClosing(file);
}

/* some simple variable functions implemented as defines */
#undef SCIPbtIsEmpty
#undef SCIPbtGetRoot

/** returns whether the binary tree is empty (has no nodes) */
SCIP_Bool SCIPbtIsEmpty(
   SCIP_BT*              tree                /**< binary tree */
   )
{
   assert(tree != NULL);

   return (tree->root == NULL);
}

/** returns the the root node of the binary or NULL if the binary tree is empty */
SCIP_BTNODE* SCIPbtGetRoot(
   SCIP_BT*              tree                /**< tree to be evaluated */
   )
{
   assert(tree != NULL);

   return tree->root;
}

/** sets root node
 *
 *  @note The old root including the rooted subtree is not delete.
 */
void SCIPbtSetRoot(
   SCIP_BT*              tree,               /**< tree to be evaluated */
   SCIP_BTNODE*          root                /**< new root, or NULL */
   )
{
   assert(tree != NULL);

   tree->root = root;
}


/*
 * Numerical methods
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
SCIP_Real SCIPcalcMachineEpsilon(
   void
   )
{
   SCIP_Real eps;
   SCIP_Real lasteps;
   SCIP_Real one;
   SCIP_Real onepluseps;

   one = 1.0;
   eps = 1.0;
   do
   {
      lasteps = eps;
      eps /= 2.0;
      onepluseps = one + eps;
   }
   while( onepluseps > one );

   return lasteps;
}

/** calculates the greatest common divisor of the two given values */
SCIP_Longint SCIPcalcGreComDiv(
   SCIP_Longint          val1,               /**< first value of greatest common devisor calculation */
   SCIP_Longint          val2                /**< second value of greatest common devisor calculation */
   )
{
   int t;

   assert(val1 > 0);
   assert(val2 > 0);

   t = 0;
   /* if val1 is even, divide it by 2 */
   while( !(val1 & 1) )
   {
      val1 >>= 1; /*lint !e704*/

      /* if val2 is even too, divide it by 2 and increase t(=number of e) */
      if( !(val2 & 1) )
      {
         val2 >>= 1; /*lint !e704*/
         ++t;
      }
      /* only val1 can be odd */
      else
      {
         /* while val1 is even, divide it by 2 */
         while( !(val1 & 1) )
            val1 >>= 1; /*lint !e704*/

         break;
      }
   }

   /* while val2 is even, divide it by 2 */
   while( !(val2 & 1) )
      val2 >>= 1; /*lint !e704*/

   /* the following if/else condition is only to make sure that we do not overflow when adding up both values before
    * dividing them by 4 in the following while loop
    */
   if( t == 0 )
   {
      if( val1 > val2 )
      {
         val1 -= val2;

         /* divide val1 by 2 as long as possible  */
         while( !(val1 & 1) )
            val1 >>= 1;   /*lint !e704*/
      }
      else if( val1 < val2 )
      {
         val2 -= val1;

         /* divide val2 by 2 as long as possible  */
         while( !(val2 & 1) )
            val2 >>= 1;   /*lint !e704*/
      }
   }

   /* val1 and val2 are odd */
   while( val1 != val2 )
   {
      if( val1 > val2 )
      {
         /* we can stop if one value reached one */
         if( val2 == 1 )
            return (val2 << t);  /*lint !e647 !e703*/

         /* if ((val1 xor val2) and 2) = 2, then gcd(val1, val2) = gcd((val1 + val2)/4, val2),
          * and otherwise                        gcd(val1, val2) = gcd((val1 − val2)/4, val2)
          */
         if( ((val1 ^ val2) & 2) == 2 )
            val1 += val2;
         else
            val1 -= val2;

         assert((val1 & 3) == 0);
         val1 >>= 2;   /*lint !e704*/

         /* if val1 is still even, divide it by 2  */
         while( !(val1 & 1) )
            val1 >>= 1;   /*lint !e704*/
      }
      else
      {
         /* we can stop if one value reached one */
         if( val1 == 1 )
            return (val1 << t);  /*lint !e647 !e703*/

         /* if ((val2 xor val1) and 2) = 2, then gcd(val2, val1) = gcd((val2 + val1)/4, val1),
          * and otherwise                        gcd(val2, val1) = gcd((val2 − val1)/4, val1)
          */
         if( ((val2 ^ val1) & 2) == 2 )
            val2 += val1;
         else
            val2 -= val1;

         assert((val2 & 3) == 0);
         val2 >>= 2;   /*lint !e704*/

         /* if val2 is still even, divide it by 2  */
         while( !(val2 & 1) )
            val2 >>= 1;   /*lint !e704*/
      }
   }

   return (val1 << t);  /*lint !e703*/
}


/* for the MS compiler, the function nextafter is named _nextafter */
#if defined(_MSC_VER) && !defined(NO_NEXTAFTER)
#define nextafter(x,y) _nextafter(x,y)
#endif

/* on systems where the function nextafter is not defined, we provide an implementation from Sun */
#ifdef NO_NEXTAFTER
/* The following implementation of the routine nextafter() comes with the following license:
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x

static
double nextafter(double x, double y)
{
   int hx;
   int hy;
   int ix;
   int iy;
   unsigned lx;
   unsigned ly;

   /* cppcheck-suppress invalidPointerCast */
   hx = __HI(x);     /* high word of x */
   /* cppcheck-suppress invalidPointerCast */
   lx = __LO(x);     /* low  word of x */
   /* cppcheck-suppress invalidPointerCast */
   hy = __HI(y);     /* high word of y */
   /* cppcheck-suppress invalidPointerCast */
   ly = __LO(y);     /* low  word of y */
   ix = hx&0x7fffffff;     /* |x| */
   iy = hy&0x7fffffff;     /* |y| */

   if( ((ix>=0x7ff00000) && ((ix-0x7ff00000)|lx) != 0 ) ||   /* x is nan */
      ( (iy>=0x7ff00000) && ((iy-0x7ff00000)|ly) != 0 ))     /* y is nan */
      return x + y;

   /* x == y, return x */
   if( x == y )
      return x;

   /* x == 0 */
   if( (ix|lx) == 0 )
   {
      /* return +-minsubnormal */
      /* cppcheck-suppress invalidPointerCast */
      __HI(x) = hy&0x80000000;
      /* cppcheck-suppress invalidPointerCast */
      __LO(x) = 1;
      y = x * x;
      if ( y == x )
         return y;
      else
         return x;  /* raise underflow flag */
   }
   /* x > 0 */
   if( hx >= 0 )
   {
      /* x > y, x -= ulp */
      if( hx > hy || ((hx == hy) && (lx > ly)) )
      {
         if ( lx == 0 )
            hx -= 1;
         lx -= 1;
      }
      else
      {
         /* x < y, x += ulp */
         lx += 1;
         if ( lx == 0 )
            hx += 1;
      }
   }
   else
   {
      /* x < 0 */
      if( hy >= 0 || hx > hy || ((hx == hy) && (lx > ly)) )
      {
         /* x < y, x -= ulp */
         if ( lx == 0 )
            hx -= 1;
         lx -= 1;
      }
      else
      {
         /* x > y, x += ulp */
         lx += 1;
         if( lx == 0 )
            hx += 1;
      }
   }
   hy = hx&0x7ff00000;
   /* overflow  */
   if( hy >= 0x7ff00000 )
      return x + x;
   if( hy < 0x00100000 )
   {
      /* underflow */
      y = x*x;
      if( y != x )
      {
         /* raise underflow flag */
         /* cppcheck-suppress invalidPointerCast */
         __HI(y) = hx;
         /* cppcheck-suppress invalidPointerCast */
         __LO(y) = lx;
         return y;
      }
   }

   /* cppcheck-suppress invalidPointerCast */
   __HI(x) = hx;
   /* cppcheck-suppress invalidPointerCast */
   __LO(x) = lx;
   return x;
}
#endif


/** returns the next representable value of from in the direction of to */
SCIP_Real SCIPnextafter(
   SCIP_Real             from,               /**< value from which the next representable value should be returned */
   SCIP_Real             to                  /**< direction in which the next representable value should be returned */
   )
{
   return nextafter(from, to);
}

/** calculates the smallest common multiple of the two given values */
SCIP_Longint SCIPcalcSmaComMul(
   SCIP_Longint          val1,               /**< first value of smallest common multiple calculation */
   SCIP_Longint          val2                /**< second value of smallest common multiple calculation */
   )
{
   SCIP_Longint gcd;

   assert(val1 > 0);
   assert(val2 > 0);

   gcd = SCIPcalcGreComDiv(val1, val2);

   return val1/gcd * val2;
}

static const SCIP_Real simplednoms[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
                                   17.0, 18.0, 19.0, 25.0, -1.0};

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
SCIP_Bool SCIPrealToRational(
   SCIP_Real             val,                /**< real value r to convert into rational number */
   SCIP_Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   SCIP_Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   )
{
   SCIP_Real a;
   SCIP_Real b;
   SCIP_Real g0;
   SCIP_Real g1;
   SCIP_Real gx;
   SCIP_Real h0;
   SCIP_Real h1;
   SCIP_Real hx;
   SCIP_Real delta0;
   SCIP_Real delta1;
   SCIP_Real epsilon;
   int i;

   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(nominator != NULL);
   assert(denominator != NULL);

   /* try the simple denominators first: each value of the simpledenoms table multiplied by powers of 10
    * is tried as denominator
    */
   for( i = 0; simplednoms[i] > 0.0; ++i )
   {
      SCIP_Real nom;
      SCIP_Real dnom;
      SCIP_Real ratval0;
      SCIP_Real ratval1;

      /* try powers of 10 (including 10^0) */
      dnom = simplednoms[i];
      while( dnom <= maxdnom )
      {
         nom = floor(val * dnom);
         ratval0 = nom/dnom;
         ratval1 = (nom+1.0)/dnom;
         if( mindelta <= val - ratval0 && val - ratval1 <= maxdelta )
         {
            if( val - ratval0 <= maxdelta )
            {
               *nominator = (SCIP_Longint)nom;
               *denominator = (SCIP_Longint)dnom;
               return TRUE;
            }
            if( mindelta <= val - ratval1 )
            {
               *nominator = (SCIP_Longint)(nom+1.0);
               *denominator = (SCIP_Longint)dnom;
               return TRUE;
            }
         }
         dnom *= 10.0;
      }
   }

   /* the simple denominators didn't work: calculate rational representation with arbitrary denominator */
   epsilon = MIN(-mindelta, maxdelta)/2.0;

   b = val;
   a = EPSFLOOR(b, epsilon);
   g0 = a;
   h0 = 1.0;
   g1 = 1.0;
   h1 = 0.0;
   delta0 = val - g0/h0;
   delta1 = (delta0 < 0.0 ? val - (g0-1.0)/h0 : val - (g0+1.0)/h0);

   while( (delta0 < mindelta || delta0 > maxdelta) && (delta1 < mindelta || delta1 > maxdelta) )
   {
      assert(EPSGT(b, a, epsilon));
      assert(h0 >= 0.0);
      assert(h1 >= 0.0);

      b = 1.0 / (b - a);
      a = EPSFLOOR(b, epsilon);

      assert(a >= 0.0);
      gx = g0;
      hx = h0;

      g0 = a * g0 + g1;
      h0 = a * h0 + h1;

      g1 = gx;
      h1 = hx;

      if( h0 > maxdnom )
         return FALSE;

      delta0 = val - g0/h0;
      delta1 = (delta0 < 0.0 ? val - (g0-1.0)/h0 : val - (g0+1.0)/h0);
   }

   if( REALABS(g0) > (SCIP_LONGINT_MAX >> 4) || h0 > (SCIP_LONGINT_MAX >> 4) )
      return FALSE;

   assert(h0 > 0.5);

   if( delta0 < mindelta )
   {
      assert(mindelta <= delta1 && delta1 <= maxdelta);
      *nominator = (SCIP_Longint)(g0 - 1.0);
      *denominator = (SCIP_Longint)h0;
   }
   else if( delta0 > maxdelta )
   {
      assert(mindelta <= delta1 && delta1 <= maxdelta);
      *nominator = (SCIP_Longint)(g0 + 1.0);
      *denominator = (SCIP_Longint)h0;
   }
   else
   {
      *nominator = (SCIP_Longint)g0;
      *denominator = (SCIP_Longint)h0;
   }
   assert(*denominator >= 1);
   assert(val - (SCIP_Real)(*nominator)/(SCIP_Real)(*denominator) >= mindelta);
   assert(val - (SCIP_Real)(*nominator)/(SCIP_Real)(*denominator) <= maxdelta);

   return TRUE;
}

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = floor(sval);
   upval = ceil(sval);

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}

/** additional scalars that are tried in integrality scaling */
static const SCIP_Real scalars[] = {3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};
static const int nscalars = 9;

/** tries to find a value, such that all given values, if scaled with this value become integral in relative allowed
 *  difference in between mindelta and maxdelta
 */
SCIP_RETCODE SCIPcalcIntegralScalar(
   SCIP_Real*            vals,               /**< values to scale */
   int                   nvals,              /**< number of values to scale */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   )
{
   SCIP_Real bestscalar;
   SCIP_Longint gcd;
   SCIP_Longint scm;
   SCIP_Longint nominator;
   SCIP_Longint denominator;
   SCIP_Real val;
   SCIP_Real minval;
   SCIP_Real absval;
   SCIP_Real scaleval;
   SCIP_Bool scalable;
   SCIP_Bool rational;
   int c;
   int s;
   int i;

   assert(vals != NULL);
   assert(nvals >= 0);
   assert(maxdnom >= 1);
   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(success != NULL);

   SCIPdebugMessage("trying to find rational representation for given values\n");

   if( intscalar != NULL )
      *intscalar = SCIP_INVALID;
   *success = FALSE;

   /* get minimal absolute non-zero value */
   minval = SCIP_REAL_MAX;
   for( c = 0; c < nvals; ++c )
   {
      val = vals[c];
      if( val < mindelta || val > maxdelta )
      {
         absval = REALABS(val);
         minval = MIN(minval, absval);
      }
   }

   if( minval == SCIP_REAL_MAX ) /*lint !e777*/
   {
      /* all coefficients are zero (inside tolerances) */
      if( intscalar != NULL )
         *intscalar = 1.0;
      *success = TRUE;
      SCIPdebugMessage(" -> all values are zero (inside tolerances)\n");

      return SCIP_OKAY;
   }
   assert(minval > MIN(-mindelta, maxdelta));

   bestscalar = SCIP_INVALID;

   for( i = 0; i < 2; ++i )
   {
      scalable = TRUE;

      /* try, if values can be made integral multiplying them with the reciprocal of the smallest value and a power of 2 */
      if( i == 0 )
	 scaleval = 1.0/minval;
      /* try, if values can be made integral by multiplying them by a power of 2 */
      else
	 scaleval = 1.0;

      for( c = 0; c < nvals && scalable; ++c )
      {
	 /* check, if the value can be scaled with a simple scalar */
	 val = vals[c];
	 if( val == 0.0 ) /* zeros are allowed in the vals array */
	    continue;

	 absval = REALABS(val);
	 while( scaleval <= maxscale
	    && (absval * scaleval < 0.5 || !isIntegralScalar(val, scaleval, mindelta, maxdelta)) )
	 {
	    for( s = 0; s < nscalars; ++s )
	    {
	       if( isIntegralScalar(val, scaleval * scalars[s], mindelta, maxdelta) )
	       {
		  scaleval *= scalars[s];
		  break;
	       }
	    }
	    if( s >= nscalars )
	       scaleval *= 2.0;
	 }
	 scalable = (scaleval <= maxscale);
	 SCIPdebugMessage(" -> val=%g, scaleval=%g, val*scaleval=%g, scalable=%u\n",
	    val, scaleval, val*scaleval, scalable);
      }
      if( scalable )
      {
	 /* make values integral by dividing them by the smallest value (and multiplying them with a power of 2) */
	 assert(scaleval <= maxscale);

	 /* check if we found a better scaling value */
	 if( scaleval < bestscalar )
	    bestscalar = scaleval;

	 SCIPdebugMessage(" -> integrality could be achieved by scaling with %g\n", scaleval);

	 /* if the scalar is still the reciprocal of the minimal value, all coeffcients are the same and we do not get a better scalar */
	 if( i == 0 && EPSEQ(scaleval, 1.0/minval, SCIP_DEFAULT_EPSILON) )
	 {
	    if( intscalar != NULL )
	       *intscalar = bestscalar;
	    *success = TRUE;

	    return SCIP_OKAY;
	 }
      }
   }

   /* convert each value into a rational number, calculate the greatest common divisor of the nominators
    * and the smallest common multiple of the denominators
    */
   gcd = 1;
   scm = 1;
   rational = TRUE;

   /* first value (to initialize gcd) */
   for( c = 0; c < nvals && rational; ++c )
   {
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;

      rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
      if( rational && nominator != 0 )
      {
         assert(denominator > 0);
         gcd = ABS(nominator);
         scm = denominator;
         rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
         SCIPdebugMessage(" -> c=%d first rational: val: %g == %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ", gcd=%" SCIP_LONGINT_FORMAT ", scm=%" SCIP_LONGINT_FORMAT ", rational=%u\n",
            c, val, nominator, denominator, gcd, scm, rational);
         break;
      }
   }

   /* remaining values */
   for( ++c; c < nvals && rational; ++c )
   {
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;

      rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
      if( rational && nominator != 0 )
      {
         assert(denominator > 0);
         gcd = SCIPcalcGreComDiv(gcd, ABS(nominator));
         scm *= denominator / SCIPcalcGreComDiv(scm, denominator);
         rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
         SCIPdebugMessage(" -> c=%d next rational : val: %g == %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ", gcd=%" SCIP_LONGINT_FORMAT ", scm=%" SCIP_LONGINT_FORMAT ", rational=%u\n",
            c, val, nominator, denominator, gcd, scm, rational);
      }
      else
      {
         SCIPdebugMessage(" -> failed to convert %g into a rational representation\n", val);
      }
   }

   if( rational )
   {
      /* make values integral by multiplying them with the smallest common multiple of the denominators */
      assert((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);

      /* check if we found a better scaling value */
      if( (SCIP_Real)scm/(SCIP_Real)gcd < bestscalar )
	 bestscalar = (SCIP_Real)scm/(SCIP_Real)gcd;

      SCIPdebugMessage(" -> integrality could be achieved by scaling with %g (rational:%" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ")\n",
         (SCIP_Real)scm/(SCIP_Real)gcd, scm, gcd);
   }

   if( bestscalar < SCIP_INVALID )
   {
      if( intscalar != NULL )
         *intscalar = bestscalar;
      *success = TRUE;

      SCIPdebugMessage(" -> smallest value to achieve integrality is %g \n", bestscalar);
   }

   return SCIP_OKAY;
}

/** given a (usually very small) interval, tries to find a rational number with simple denominator (i.e. a small
 *  number, probably multiplied with powers of 10) out of this interval; returns TRUE iff a valid rational
 *  number inside the interval was found
 */
SCIP_Bool SCIPfindSimpleRational(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed for resulting rational number */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   )
{
   SCIP_Real center;
   SCIP_Real delta;

   assert(lb <= ub);

   center = 0.5*(lb+ub);

   /* in order to compute a rational number that is exactly within the bounds (as the user expects),
    * we computed the allowed delta with downward rounding, if available
    */
   if( SCIPintervalHasRoundingControl() )
   {
      SCIP_ROUNDMODE roundmode;

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      delta = 0.5*(ub-lb);

      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      delta = 0.5*(ub-lb);
   }

   return SCIPrealToRational(center, -delta, +delta, maxdnom, nominator, denominator);
}

/** given a (usually very small) interval, selects a value inside this interval; it is tried to select a rational number
 *  with simple denominator (i.e. a small number, probably multiplied with powers of 10);
 *  if no valid rational number inside the interval was found, selects the central value of the interval
 */
SCIP_Real SCIPselectSimpleValue(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom             /**< maximal denominator allowed for resulting rational number */
   )
{
   SCIP_Real val;

   val = 0.5*(lb+ub);
   if( lb < ub )
   {
      SCIP_Longint nominator;
      SCIP_Longint denominator;
      SCIP_Bool success;

      /* try to find a "simple" rational number inside the interval */
      SCIPdebugMessage("simple rational in [%.9f,%.9f]:", lb, ub);
      success = SCIPfindSimpleRational(lb, ub, maxdnom, &nominator, &denominator);
      if( success )
      {
         val = (SCIP_Real)nominator/(SCIP_Real)denominator;
         SCIPdebugPrintf(" %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT " == %.9f\n", nominator, denominator, val);

         if( val - lb < 0.0 || val - ub > 0.0 )
         {
            SCIPdebugPrintf(" value is out of interval bounds by %g -> failed\n", MAX(lb-val, val-ub));
            val = 0.5*(lb+ub);
         }
      }
      else
      {
         SCIPdebugPrintf(" failed\n");
      }
   }

   return val;
}




/*
 * Random Numbers
 */

#if defined(NO_RAND_R) || defined(_WIN32) || defined(_WIN64)

#define SCIP_RAND_MAX 32767
/** returns a random number between 0 and SCIP_RAND_MAX */
static
int getRand(
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   SCIP_Longint nextseed;

   assert(seedp != NULL);

   nextseed = (*seedp) * (SCIP_Longint)1103515245 + 12345;
   *seedp = (unsigned int)nextseed;

   return (int)((unsigned int)(nextseed/(2*(SCIP_RAND_MAX+1))) % (SCIP_RAND_MAX+1));
}

#else

#define SCIP_RAND_MAX RAND_MAX

/** returns a random number between 0 and SCIP_RAND_MAX */
static
int getRand(
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return rand_r(seedp);
}

#endif

/** returns a random integer between minrandval and maxrandval */
static
int getRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   SCIP_Real randnumber;

   randnumber = (SCIP_Real)getRand(seedp)/(SCIP_RAND_MAX+1.0);
   assert(randnumber >= 0.0);
   assert(randnumber < 1.0);

   /* we multiply minrandval and maxrandval separately by randnumber in order to avoid overflow if they are more than INT_MAX
    * apart
    */
   return (int) (minrandval*(1.0 - randnumber) + maxrandval*randnumber + randnumber);
}

/** returns a random real between minrandval and maxrandval */
static
SCIP_Real getRandomReal(
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   SCIP_Real randnumber;

   randnumber = (SCIP_Real)getRand(seedp)/(SCIP_Real)SCIP_RAND_MAX;
   assert(randnumber >= 0.0);
   assert(randnumber <= 1.0);

   /* we multiply minrandval and maxrandval separately by randnumber in order to avoid overflow if they are more than
    * SCIP_REAL_MAX apart
    */
   return minrandval*(1.0 - randnumber) + maxrandval*randnumber;
}

/** returns a random integer between minrandval and maxrandval
 *
 *  @deprecated Please use SCIPrandomGetInt() to request a random integer.
 */
int SCIPgetRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return getRandomInt(minrandval, maxrandval, seedp);
}

/** returns a random real between minrandval and maxrandval
 *
 *  @deprecated Please use SCIPrandomGetReal() to request a random real.
 */
SCIP_Real SCIPgetRandomReal(
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return getRandomReal(minrandval, maxrandval, seedp);
}


/* initial seeds for KISS random number generator */
#define DEFAULT_SEED UINT32_C(123456789)
#define DEFAULT_XOR  UINT32_C(362436000)
#define DEFAULT_MWC  UINT32_C(521288629)
#define DEFAULT_CST  UINT32_C(7654321)


/** initializes a random number generator with a given start seed */
void SCIPrandomSetSeed(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   unsigned int          initseed            /**< initial random seed */
   )
{
   assert(randnumgen != NULL);

   /* use MAX() to avoid zero after over flowing */
   randnumgen->seed = MAX(SCIPhashTwo(DEFAULT_SEED, initseed), 1u);
   randnumgen->xor_seed = MAX(SCIPhashTwo(DEFAULT_XOR, initseed), 1u);
   randnumgen->mwc_seed = MAX(SCIPhashTwo(DEFAULT_MWC, initseed), 1u);
   randnumgen->cst_seed = SCIPhashTwo(DEFAULT_CST, initseed);

   assert(randnumgen->seed > 0);
   assert(randnumgen->xor_seed > 0);
   assert(randnumgen->mwc_seed > 0);
}

/** returns a random number between 0 and UINT32_MAX
 *
 *  implementation of KISS random number generator developed by George Marsaglia.
 *  KISS is combination of three different random number generators:
 *   - Linear congruential generator
 *   - Xorshift
 *   - Lag-1 Multiply-with-carry
 *
 *  KISS has a period of 2^123 and passes all statistical test part of BigCrush-Test of TestU01 [1].
 *
 *  [1] http://dl.acm.org/citation.cfm?doid=1268776.1268777
 */
static
uint32_t randomGetRand(
   SCIP_RANDNUMGEN*      randnumgen          /**< random number generator */
   )
{
   uint64_t t;

   /* linear congruential */
   randnumgen->seed = (uint32_t) (randnumgen->seed * UINT64_C(1103515245) + UINT64_C(12345));

   /* Xorshift */
   randnumgen->xor_seed ^= (randnumgen->xor_seed << 13);
   randnumgen->xor_seed ^= (randnumgen->xor_seed >> 17);
   randnumgen->xor_seed ^= (randnumgen->xor_seed << 5);

   /* Multiply-with-carry */
   t = UINT64_C(698769069) * randnumgen->mwc_seed + randnumgen->cst_seed;
   randnumgen->cst_seed = (uint32_t) (t >> 32);
   randnumgen->mwc_seed = (uint32_t) t;

   return randnumgen->seed + randnumgen->xor_seed + randnumgen->mwc_seed;
}

/** creates and initializes a random number generator */
SCIP_RETCODE SCIPrandomCreate(
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          initialseed         /**< initial random seed */
   )
{
   assert(randnumgen != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, randnumgen) );

   SCIPrandomSetSeed((*randnumgen), initialseed);

   return SCIP_OKAY;
}

/** frees a random number generator */
void SCIPrandomFree(
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(randnumgen != NULL);
   assert((*randnumgen) != NULL);

   BMSfreeBlockMemory(blkmem, randnumgen);

   return;
}



/** returns a random integer between minrandval and maxrandval */
int SCIPrandomGetInt(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval          /**< maximal value to return */
   )
{
   SCIP_Real randnumber;
   SCIP_Longint zeromax;

   randnumber = (SCIP_Real)randomGetRand(randnumgen)/(UINT32_MAX+1.0);
   assert(randnumber >= 0.0);
   assert(randnumber < 1.0);

   /* we need to shift the range to the non-negative integers to handle negative integer values correctly.
    * we use a long integer to avoid overflows.
    */
   zeromax = (SCIP_Longint)maxrandval - (SCIP_Longint)minrandval + 1;

   return (int) ((SCIP_Longint)(zeromax * randnumber) + (SCIP_Longint)minrandval);
}

/** returns a random real between minrandval and maxrandval */
SCIP_Real SCIPrandomGetReal(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval          /**< maximal value to return */
   )
{
   SCIP_Real randnumber;

   randnumber = (SCIP_Real)randomGetRand(randnumgen)/(SCIP_Real)UINT32_MAX;
   assert(randnumber >= 0.0);
   assert(randnumber <= 1.0);

   /* we multiply minrandval and maxrandval separately by randnumber in order to avoid overflow if they are more than
    * SCIP_REAL_MAX apart
    */
   return minrandval*(1.0 - randnumber) + maxrandval*randnumber;
}

/** randomly shuffles parts of an integer array using the Fisher-Yates algorithm */
void SCIPrandomPermuteIntArray(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   int*                  array,              /**< array to be shuffled */
   int                   begin,              /**< first included index that should be subject to shuffling
                                              *   (0 for first array entry)
                                              */
   int                   end                 /**< first excluded index that should not be subject to shuffling
                                              *   (array size for last array entry)
                                              */
   )
{
   int tmp;
   int i;

   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 )
   {
      --end;

      /* get a random position into which the last entry should be shuffled */
      i = SCIPrandomGetInt(randnumgen, begin, end);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}

/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
void SCIPrandomPermuteArray(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first included index that should be subject to shuffling
                                              *   (0 for first array entry)
                                              */
   int                   end                 /**< first excluded index that should not be subject to shuffling
                                              *   (array size for last array entry)
                                              */
   )
{
   void* tmp;
   int i;

   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 )
   {
      end--;

      /* get a random position into which the last entry should be shuffled */
      i = SCIPrandomGetInt(randnumgen, begin, end);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
SCIP_RETCODE SCIPrandomGetSubset(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems           /**< number of elements that should be drawn and stored */
   )
{
   int i;
   int j;

   /* if both sets are of equal size, we just copy the array */
   if( nelems == nsubelems)
   {
      BMScopyMemoryArray(subset,set,nelems);
      return SCIP_OKAY;
   }

   /* abort, if size of subset is too big */
   if( nsubelems > nelems )
   {
      SCIPerrorMessage("Cannot create %d-elementary subset of %d-elementary set.\n", nsubelems, nelems);
      return SCIP_INVALIDDATA;
   }
#ifndef NDEBUG
   for( i = 0; i < nsubelems; i++ )
      for( j = 0; j < i; j++ )
         assert(set[i] != set[j]);
#endif

   /* draw each element individually */
   i = 0;
   while( i < nsubelems )
   {
      int r;

      r = SCIPrandomGetInt(randnumgen, 0, nelems-1);
      subset[i] = set[r];

      /* if we get an element that we already had, we will draw again */
      for( j = 0; j < i; j++ )
      {
         if( subset[i] == subset[j] )
         {
            --i;
            break;
         }
      }
      ++i;
   }
   return SCIP_OKAY;
}

/*
 * Additional math functions
 */

/** calculates a binomial coefficient n over m, choose m elements out of n, maximal value will be 33 over 16 (because
 *  the n=33 is the last line in the Pascal's triangle where each entry fits in a 4 byte value), an error occurs due to
 *  big numbers or an negative value m (and m < n) and -1 will be returned
 */
SCIP_Longint SCIPcalcBinomCoef(
   int                   n,                  /**< number of different elements */
   int                   m                   /**< number to choose out of the above */
   )
{
   if( m == 0 || m >= n )
      return 1;

   if( m < 0 )
      return -1;

   /* symmetry of the binomial coefficient, choose smaller m */
   if( m > n/2 )
      m = n - m;

   /* trivial case m == 1 */
   if( m == 1 )
      return n;

   /* simple case m == 2 */
   if( m == 2 )
   {
      if( ((SCIP_Real)SCIP_LONGINT_MAX) / n >= (n-1) * 2 ) /*lint !e790*/
         return ((SCIP_Longint)n*(n-1)/2); /*lint !e647*/
      else
         return -1;
   }

   /* abort on to big numbers */
   if( m > 16 || n > 33 )
      return -1;

   /* simple case m == 3 */
   if( m == 3 )
      return (n*(n-1)*(n-2)/6); /*lint !e647*/
   else
   {
      /* first half of Pascal's triangle numbers(without the symmetric part) backwards from (33,16) over (32,16),
       * (33,15), (32,15),(31,15, (30,15), (33,14) to (8,4) (rest is calculated directly)
       *
       * due to this order we can extract the right binomial coefficient by (16-m)^2+(16-m)+(33-n)
       */
      static const SCIP_Longint binoms[182] = {
         1166803110, 601080390, 1037158320, 565722720, 300540195, 155117520, 818809200, 471435600, 265182525, 145422675,
         77558760, 40116600, 573166440, 347373600, 206253075, 119759850, 67863915, 37442160, 20058300, 10400600,
         354817320, 225792840, 141120525, 86493225, 51895935, 30421755, 17383860, 9657700, 5200300, 2704156, 193536720,
         129024480, 84672315, 54627300, 34597290, 21474180, 13037895, 7726160, 4457400, 2496144, 1352078, 705432,
         92561040, 64512240, 44352165, 30045015, 20030010, 13123110, 8436285, 5311735, 3268760, 1961256, 1144066,
         646646, 352716, 184756, 38567100, 28048800, 20160075, 14307150, 10015005, 6906900, 4686825, 3124550, 2042975,
         1307504, 817190, 497420, 293930, 167960, 92378, 48620, 13884156, 10518300, 7888725, 5852925, 4292145, 3108105,
         2220075, 1562275, 1081575, 735471, 490314, 319770, 203490, 125970, 75582, 43758, 24310, 12870, 4272048, 3365856,
         2629575, 2035800, 1560780, 1184040, 888030, 657800, 480700, 346104, 245157, 170544, 116280, 77520, 50388, 31824,
         19448, 11440, 6435, 3432, 1107568, 906192, 736281, 593775, 475020, 376740, 296010, 230230, 177100, 134596,
         100947, 74613, 54264, 38760, 27132, 18564, 12376, 8008, 5005, 3003, 1716, 924, 237336, 201376, 169911, 142506,
         118755, 98280, 80730, 65780, 53130, 42504, 33649, 26334, 20349, 15504, 11628, 8568, 6188, 4368, 3003, 2002,
         1287, 792, 462, 252, 40920, 35960, 31465, 27405, 23751, 20475, 17550, 14950, 12650, 10626, 8855, 7315, 5985,
         4845, 3876, 3060, 2380, 1820, 1365, 1001, 715, 495, 330, 210, 126, 70};

      /* m can at most be 16 */
      const int t = 16-m;
      assert(t >= 0);
      assert(n <= 33);

      /* binoms array hast exactly 182 elements */
      assert(t*(t+1)+(33-n) < 182);

      return binoms[t*(t+1)+(33-n)]; /*lint !e662 !e661*/
   }
}

/** negates a number */
SCIP_Real SCIPnegateReal(
   SCIP_Real             x                   /**< value to negate */
   )
{
   return -x;
}

/*
 * Permutations / Shuffling
 */

/** swaps two ints */
void SCIPswapInts(
   int*                  value1,             /**< pointer to first integer */
   int*                  value2              /**< pointer to second integer */
   )
{
   int tmp;

   tmp = *value1;
   *value1 = *value2;
   *value2 = tmp;
}

/** swaps two real values */
void SCIPswapReals(
   SCIP_Real*            value1,             /**< pointer to first real value */
   SCIP_Real*            value2              /**< pointer to second real value */
   )
{
   SCIP_Real tmp;

   tmp = *value1;
   *value1 = *value2;
   *value2 = tmp;
}

/** swaps the addresses of two pointers */
void SCIPswapPointers(
   void**                pointer1,           /**< first pointer */
   void**                pointer2            /**< second pointer */
   )
{
   void* tmp;

   tmp = *pointer1;
   *pointer1 = *pointer2;
   *pointer2 = tmp;
}

/** randomly shuffles parts of an integer array using the Fisher-Yates algorithm
 *
 *  @deprecated Please use SCIPrandomPermuteIntArray()
 */
void SCIPpermuteIntArray(
   int*                  array,              /**< array to be shuffled */
   int                   begin,              /**< first included index that should be subject to shuffling
                                              *   (0 for first array entry)
                                              */
   int                   end,                /**< first excluded index that should not be subject to shuffling
                                              *   (array size for last array entry)
                                              */
   unsigned int*         randseed            /**< seed value for the random generator */
   )
{
   int tmp;
   int i;

   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 )
   {
      --end;

      /* get a random position into which the last entry should be shuffled */
      i = getRandomInt(begin, end, randseed);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}


/** randomly shuffles parts of an array using the Fisher-Yates algorithm
 *
 *  @deprecated Please use SCIPrandomPermuteArray()
 */
void SCIPpermuteArray(
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first included index that should be subject to shuffling
                                              *   (0 for first array entry)
                                              */
   int                   end,                /**< first excluded index that should not be subject to shuffling
                                              *   (array size for last array entry)
                                              */
   unsigned int*         randseed            /**< seed value for the random generator */
   )
{
   void* tmp;
   int i;

   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 )
   {
      end--;

      /* get a random position into which the last entry should be shuffled */
      i = getRandomInt(begin, end, randseed);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 *
 *  @deprecated Please use SCIPrandomGetSubset()
 */
SCIP_RETCODE SCIPgetRandomSubset(
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems,          /**< number of elements that should be drawn and stored */
   unsigned int          randseed            /**< seed value for random generator */
   )
{
   int i;
   int j;

   /* if both sets are of equal size, we just copy the array */
   if( nelems == nsubelems)
   {
      BMScopyMemoryArray(subset,set,nelems);
      return SCIP_OKAY;
   }

   /* abort, if size of subset is too big */
   if( nsubelems > nelems )
   {
      SCIPerrorMessage("Cannot create %d-elementary subset of %d-elementary set.\n", nsubelems, nelems);
      return SCIP_INVALIDDATA;
   }
#ifndef NDEBUG
   for( i = 0; i < nsubelems; i++ )
      for( j = 0; j < i; j++ )
         assert(set[i] != set[j]);
#endif

   /* draw each element individually */
   i = 0;
   while( i < nsubelems )
   {
      int r;

      r = getRandomInt(0, nelems-1, &randseed);
      subset[i] = set[r];

      /* if we get an element that we already had, we will draw again */
      for( j = 0; j < i; j++ )
      {
         if( subset[i] == subset[j] )
         {
            --i;
            break;
         }
      }
      ++i;
   }
   return SCIP_OKAY;
}


/*
 * Arrays
 */

/** computes set intersection (duplicates removed) of two integer arrays that are ordered ascendingly */
SCIP_RETCODE SCIPcomputeArraysIntersection(
   int*                  array1,             /**< first array (in ascending order) */
   int                   narray1,            /**< number of entries of first array */
   int*                  array2,             /**< second array (in ascending order) */
   int                   narray2,            /**< number of entries of second array */
   int*                  intersectarray,     /**< intersection of array1 and array2
                                              *   (note: it is possible to use array1 for this input argument) */
   int*                  nintersectarray     /**< pointer to store number of entries of intersection array
                                              *   (note: it is possible to use narray1 for this input argument) */
   )
{
   int cnt = 0;
   int k = 0;
   int v1;
   int v2;

   assert( array1 != NULL );
   assert( array2 != NULL );
   assert( intersectarray != NULL );
   assert( nintersectarray != NULL );

   /* determine intersection of array1 and array2 */
   for (v1 = 0; v1 < narray1; ++v1)
   {
      assert( v1 == 0 || array1[v1] >= array1[v1-1] );

      /* skip duplicate entries */
      if ( v1+1 < narray1 && array1[v1] == array1[v1+1])
         continue;

      for (v2 = k; v2 < narray2; ++v2)
      {
         assert( v2 == 0 || array2[v2] >= array2[v2-1] );

         if ( array2[v2] > array1[v1] )
         {
            k = v2;
            break;
         }
         else if ( array2[v2] == array1[v1] )
         {
            intersectarray[cnt++] = array2[v2];
            k = v2 + 1;
            break;
         }
      }
   }

   /* store size of intersection array */
   *nintersectarray = cnt;

   return SCIP_OKAY;
}


/** computes set difference (duplicates removed) of two integer arrays that are ordered ascendingly */
SCIP_RETCODE SCIPcomputeArraysSetminus(
   int*                  array1,             /**< first array (in ascending order) */
   int                   narray1,            /**< number of entries of first array */
   int*                  array2,             /**< second array (in ascending order) */
   int                   narray2,            /**< number of entries of second array */
   int*                  setminusarray,      /**< array to store entries of array1 that are not an entry of array2
                                              *   (note: it is possible to use array1 for this input argument) */
   int*                  nsetminusarray      /**< pointer to store number of entries of setminus array
                                              *   (note: it is possible to use narray1 for this input argument) */
   )
{
   int cnt = 0;
   int v1 = 0;
   int v2 = 0;

   assert( array1 != NULL );
   assert( array2 != NULL );
   assert( setminusarray != NULL );
   assert( nsetminusarray != NULL );

   while ( v1 < narray1 )
   {
      int entry1;

      assert( v1 == 0 || array1[v1] >= array1[v1-1] );

      /* skip duplicate entries */
      while ( v1 + 1 < narray1 && array1[v1] == array1[v1 + 1] )
         ++v1;

      entry1 = array1[v1];

      while ( v2 < narray2 && array2[v2] < entry1 )
         ++v2;

      if ( v2 >= narray2 || entry1 < array2[v2] )
         setminusarray[cnt++] = entry1;
      ++v1;
   }

   /* store size of setminus array */
   *nsetminusarray = cnt;

   return SCIP_OKAY;
}


/*
 * Strings
 */


/** copies characters from 'src' to 'dest', copying is stopped when either the 'stop' character is reached or after
 *  'cnt' characters have been copied, whichever comes first.
 *
 *  @note undefined behavior on overlapping arrays
 */
int SCIPmemccpy(
   char*                 dest,               /**< destination pointer to copy to */
   const char*           src,                /**< source pointer to copy from */
   char                  stop,               /**< character when found stop copying */
   unsigned int          cnt                 /**< maximal number of characters to copy */
   )
{
   if( dest == NULL || src == NULL || cnt == 0 )
      return -1;
   else
   {
      char* destination = dest;

      while( cnt-- && (*destination++ = *src++) != stop ); /*lint !e722*/

      return (int)(destination - dest);
   }
}

/** prints an error message containing of the given string followed by a string describing the current system error
 *
 *  Prefers to use the strerror_r method, which is threadsafe. On systems where this method does not exist,
 *  NO_STRERROR_R should be defined (see INSTALL). In this case, strerror is used which is not guaranteed to be
 *  threadsafe (on SUN-systems, it actually is).
 */
void SCIPprintSysError(
   const char*           message             /**< first part of the error message, e.g. the filename */
   )
{
#ifdef NO_STRERROR_R
   SCIPmessagePrintError("%s: %s\n", message, strerror(errno));
#else
   char buf[SCIP_MAXSTRLEN];

#if defined(_WIN32) || defined(_WIN64)
   /* strerror_s returns 0 on success; the string is \0 terminated. */
   if ( strerror_s(buf, SCIP_MAXSTRLEN, errno) != 0 )
      SCIPmessagePrintError("Unknown error number %d or error message too long.\n", errno);
   SCIPmessagePrintError("%s: %s\n", message, buf);
#elif (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && ! defined(_GNU_SOURCE)
   /* We are in the POSIX/XSI case, where strerror_r returns 0 on success; \0 termination is unclear. */
   if ( strerror_r(errno, buf, SCIP_MAXSTRLEN) != 0 )
      SCIPmessagePrintError("Unknown error number %d.\n", errno);
   buf[SCIP_MAXSTRLEN - 1] = '\0';
   SCIPmessagePrintError("%s: %s\n", message, buf);
#else
   /* We are in the GNU case, where strerror_r returns a pointer to the error string. This string is possibly stored
    * in buf and is always \0 terminated.
    * However, if compiling on one system and executing on another system, we might actually call a different
    * variant of the strerror_r function than we had at compile time.
    */
   char* errordescr;
   *buf = '\0';
   errordescr = strerror_r(errno, buf, SCIP_MAXSTRLEN);
   if( *buf != '\0' )
   {
      /* strerror_r wrote into buf */
      SCIPmessagePrintError("%s: %s\n", message, buf);
   }
   else if( errordescr != NULL )
   {
      /* strerror_r returned something non-NULL */
      SCIPmessagePrintError("%s: %s\n", message, errordescr);
   }
   else
   {
      /* strerror_r did return NULL and did not write into buf */
      SCIPmessagePrintError("Could not obtain description for error %d.\n", errno);
   }
#endif
#endif
}

/** extracts tokens from strings - wrapper method for strtok_r() */
char* SCIPstrtok(
   char*                 s,                  /**< string to parse */
   const char*           delim,              /**< delimiters for parsing */
   char**                ptrptr              /**< pointer to working char pointer - must stay the same while parsing */
   )
{
#ifdef NO_STRTOK_R
   return strtok(s, delim);
#else
   return strtok_r(s, delim, ptrptr);
#endif
}

/** translates the given string into a string where symbols ", ', and spaces are escaped with a \ prefix */
void SCIPescapeString(
   char*                 t,                  /**< target buffer to store escaped string */
   int                   bufsize,            /**< size of buffer t */
   const char*           s                   /**< string to transform into escaped string */
   )
{
   int len;
   int i;
   int p;

   assert(t != NULL);
   assert(bufsize > 0);

   len = (int)strlen(s);
   for( p = 0, i = 0; i <= len && p < bufsize; ++i, ++p )
   {
      if( s[i] == ' ' || s[i] == '"' || s[i] == '\'' )
      {
         t[p] = '\\';
         p++;
      }
      if( p < bufsize )
         t[p] = s[i];
   }
   t[bufsize-1] = '\0';
}

/* safe version of snprintf */
int SCIPsnprintf(
   char*                 t,                  /**< target string */
   int                   len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
   )
{
   va_list ap;
   int n;

   assert(t != NULL);
   assert(len > 0);

   va_start(ap, s); /*lint !e826*/

#if defined(_WIN32) || defined(_WIN64)
   n = _vsnprintf(t, (size_t) len, s, ap);
#else
   n = vsnprintf(t, (size_t) len, s, ap); /*lint !e571*/
#endif
   va_end(ap);

   if( n < 0 || n >= len )
   {
#ifndef NDEBUG
      if( n < 0 )
      {
         SCIPerrorMessage("vsnprintf returned %d\n",n);
      }
#endif
      t[len-1] = '\0';
      n = len-1;
   }
   return n;
}

/** extract the next token as a integer value if it is one; in case no value is parsed the endptr is set to @p str
 *
 *  @return Returns TRUE if a value could be extracted, otherwise FALSE
 */
SCIP_Bool SCIPstrToIntValue(
   const char*           str,                /**< string to search */
   int*                  value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   )
{
   assert(str != NULL);
   assert(value != NULL);
   assert(endptr != NULL);

   /* init errno to detect possible errors */
   errno = 0;

   *value = (int) strtol(str, endptr, 10);

   if( *endptr != str && *endptr != NULL )
   {
      SCIPdebugMessage("parsed integer value <%d>\n", *value);
      return TRUE;
   }
   *endptr = (char*)str;

   SCIPdebugMessage("failed parsing integer value <%s>\n", str);

   return FALSE;
}

/** extract the next token as a double value if it is one; in case no value is parsed the endptr is set to @p str
 *
 *  @return Returns TRUE if a value could be extracted, otherwise FALSE
 */
SCIP_Bool SCIPstrToRealValue(
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   )
{
   assert(str != NULL);
   assert(value != NULL);
   assert(endptr != NULL);

   /* init errno to detect possible errors */
   errno = 0;

   *value = strtod(str, endptr);

   if( *endptr != str && *endptr != NULL )
   {
      SCIPdebugMessage("parsed real value <%g>\n", *value);
      return TRUE;
   }
   *endptr = (char*)str;

   SCIPdebugMessage("failed parsing real value <%s>\n", str);

   return FALSE;
}

/** copies the first size characters between a start and end character of str into token, if no error occured endptr
 *  will point to the position after the read part, otherwise it will point to @p str
 */
void SCIPstrCopySection(
   const char*           str,                /**< string to search */
   char                  startchar,          /**< character which defines the beginning */
   char                  endchar,            /**< character which defines the ending */
   char*                 token,              /**< string to store the copy */
   int                   size,               /**< size of the token char array */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   )
{
   const char* copystr;
   int nchars;

   assert(str != NULL);
   assert(token != NULL);
   assert(size > 0);
   assert(endptr != NULL);

   nchars = 0;

   copystr = str;

   /* find starting character */
   while( *str != '\0' && *str != startchar )
      ++str;

   /* did not find start character */
   if( *str == '\0' )
   {
      *endptr = (char*)copystr;
      return;
   }

   /* skip start character */
   ++str;

   /* copy string */
   while( *str != '\0' && *str != endchar && nchars < size-1 )
   {
      assert(nchars < SCIP_MAXSTRLEN);
      token[nchars] = *str;
      nchars++;
      ++str;
   }

   /* add end to token */
   token[nchars] = '\0';

   /* if section was longer than size, we want to reach the end of the parsing section anyway */
   if( nchars == (size-1) )
      while( *str != '\0' && *str != endchar )
         ++str;

   /* did not find end character */
   if( *str == '\0' )
   {
      *endptr = (char*)copystr;
      return;
   }

   /* skip end character */
   ++str;

   SCIPdebugMessage("parsed section <%s>\n", token);

   *endptr = (char*) str;
}

/*
 * File methods
 */

/** returns, whether the given file exists */
SCIP_Bool SCIPfileExists(
   const char*           filename            /**< file name */
   )
{
   FILE* f;

   f = fopen(filename, "r");
   if( f == NULL )
      return FALSE;

   fclose(f);

   return TRUE;
}

/** splits filename into path, name, and extension */
void SCIPsplitFilename(
   char*                 filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**                path,               /**< pointer to store path, or NULL if not needed */
   char**                name,               /**< pointer to store name, or NULL if not needed */
   char**                extension,          /**< pointer to store extension, or NULL if not needed */
   char**                compression         /**< pointer to store compression extension, or NULL if not needed */
   )
{
   char* lastslash;
   char* lastbackslash;
   char* lastdot;

   assert(filename != NULL);

   if( path != NULL )
      *path = NULL;
   if( name != NULL )
      *name = NULL;
   if( extension != NULL )
      *extension = NULL;
   if( compression != NULL )
      *compression = NULL;

   /* treat both slashes '/' and '\' as directory delimiters */
   lastslash = strrchr(filename, '/');
   lastbackslash = strrchr(filename, '\\');
   lastslash = MAX(lastslash, lastbackslash); /*lint !e613*/
   lastdot = strrchr(filename, '.');
   if( lastslash != NULL && lastdot != NULL && lastdot < lastslash ) /* is the last dot belonging to the path? */
      lastdot = NULL;

   /* detect known compression extensions */
#ifdef WITH_ZLIB
   if( lastdot != NULL )
   {
      char* compext;

      compext = lastdot+1;
      if( strcmp(compext, "gz") == 0
        || strcmp(compext, "z") == 0
        || strcmp(compext, "Z") == 0 )
      {
         if( compression != NULL )
            *compression = compext;
         *lastdot = '\0';
      }

      /* find again the last dot in the filename without compression extension */
      lastdot = strrchr(filename, '.');
      if( lastslash != NULL && lastdot != NULL && lastdot < lastslash ) /* is the last dot belonging to the path? */
         lastdot = NULL;
   }
#endif

   if( lastslash == NULL )
   {
      if( name != NULL )
         *name = filename;
   }
   else
   {
      if( path != NULL )
         *path = filename;
      if( name != NULL )
         *name = lastslash+1;
      *lastslash = '\0';
   }

   if( lastdot != NULL )
   {
      if( extension != NULL )
         *extension = lastdot+1;
      *lastdot = '\0';
   }
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPrelDiff

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
SCIP_Real SCIPrelDiff(
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real absval1;
   SCIP_Real absval2;
   SCIP_Real quot;

   absval1 = REALABS(val1);
   absval2 = REALABS(val2);
   quot = MAX3(1.0, absval1, absval2);

   return (val1-val2)/quot;
}


/** computes the gap from the primal and the dual bound */
SCIP_Real SCIPcomputeGap(
   SCIP_Real             eps,                /**< the value treated as zero */
   SCIP_Real             inf,                /**< the value treated as infinity */
   SCIP_Real             primalbound,        /**< the primal bound */
   SCIP_Real             dualbound           /**< the dual bound */
   )
{
   if( EPSEQ(primalbound, dualbound, eps) )
      return 0.0;
   else if( EPSZ(dualbound, eps) ||
            EPSZ(primalbound, eps) ||
            REALABS(primalbound) >= inf ||
            REALABS(dualbound) >= inf ||
            primalbound * dualbound < 0.0 )
      return inf;
   else
   {
      SCIP_Real absdual = REALABS(dualbound);
      SCIP_Real absprimal = REALABS(primalbound);

      return REALABS((primalbound - dualbound)/MIN(absdual, absprimal));
   }
}

/*
 *Union-Find data structure
 */

/** creates a disjoint set (union find) structure \p uf for \p ncomponents many components (of size one) */
SCIP_RETCODE SCIPdisjointsetCreate(
   SCIP_DISJOINTSET**    djset,              /**< disjoint set (union find) data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncomponents         /**< number of components */
   )
{
   assert(djset != NULL);
   assert(blkmem != NULL);

   /* allocate the necessary memory */
   assert(ncomponents > 0);
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, djset) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &((*djset)->parents), ncomponents) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &((*djset)->sizes), ncomponents) );
   (*djset)->size = ncomponents;

   /* clear the data structure */
   SCIPdisjointsetClear(*djset);

   return SCIP_OKAY;
}

/** clears the disjoint set (union find) structure \p uf */
void SCIPdisjointsetClear(
   SCIP_DISJOINTSET*     djset               /**< disjoint set (union find) data structure */
   )
{
   int i;
   djset->componentcount = djset->size;

   /* reset all components to be unconnected */
   for( i = 0; i < djset->componentcount; i++ )
   {
      djset->parents[i] = i;
      djset->sizes[i] = 1;
   }
}


/** finds and returns the component identifier of this \p element */
int SCIPdisjointsetFind(
   SCIP_DISJOINTSET*     djset,              /**< disjoint set (union find) data structure */
   int                   element             /**< element to be found */
   )
{
   int newelement;
   int root = element;
   int* parents = djset->parents;

   /* find root of this element */
   while( root != parents[root] )
   {
      root = parents[root];
   }

   /* compress the path to make future queries faster */
   while( element != root )
   {
      newelement = parents[element];
      parents[element] = root;
      element = newelement;
   }

   return root;
}

/** merges the components containing the elements \p p and \p q */
void SCIPdisjointsetUnion(
   SCIP_DISJOINTSET*     djset,              /**< disjoint set (union find) data structure */
   int                   p,                  /**< first element */
   int                   q,                  /**< second element */
   SCIP_Bool             forcerepofp         /**< force representative of p to be new representative */
   )
{
   int idp;
   int idq;
   int* sizes;
   int* parents;

   assert(djset != NULL);
   assert(0 <= p);
   assert(0 <= q);
   assert(djset->size > p);
   assert(djset->size > q);


   idp = SCIPdisjointsetFind(djset, p);
   idq = SCIPdisjointsetFind(djset, q);

   /* if p and q lie in the same component, there is nothing to be done */
   if( idp == idq )
      return;

   sizes = djset->sizes;
   parents = djset->parents;

   if( forcerepofp )
   {
      parents[idq] = idp;
      sizes[idp] += sizes[idq];
   }
   else
   {
      if( sizes[idp] < sizes[idq] )
      {
         parents[idp] = idq;
         sizes[idq] += sizes[idp];
      }
      else
      {
         parents[idq] = idp;
         sizes[idp] += sizes[idq];
      }
   }
   /* one less component */
   djset->componentcount--;
}

/** frees the disjoint set (union find) data structure */
void SCIPdisjointsetFree(
   SCIP_DISJOINTSET**    djset,              /**< pointer to disjoint set (union find) data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_DISJOINTSET* dsptr;

   assert(djset != NULL);
   assert(*djset != NULL);

   dsptr = *djset;

   BMSfreeBlockMemoryArray(blkmem, &dsptr->sizes, dsptr->size);
   BMSfreeBlockMemoryArray(blkmem, &dsptr->parents, dsptr->size);

   BMSfreeBlockMemory(blkmem, djset);
}

/** returns the number of independent components in this disjoint set (union find) data structure */
int SCIPdisjointsetGetComponentCount(
   SCIP_DISJOINTSET*     djset               /**< disjoint set (union find) data structure */
   )
{
   assert(djset != NULL);

   return djset->componentcount;
}

/** returns the size (number of nodes) of this disjoint set (union find) data structure */
int SCIPdisjointsetGetSize(
   SCIP_DISJOINTSET*     djset               /**< disjoint set (union find) data structure */
   )
{
   assert(djset != NULL);

   return djset->size;
}
