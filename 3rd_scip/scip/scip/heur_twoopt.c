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

/**@file   heur_twoopt.c
 * @brief  primal heuristic to improve incumbent solution by flipping pairs of variables
 * @author Timo Berthold
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/heur_twoopt.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "twoopt"
#define HEUR_DESC             "primal heuristic to improve incumbent solution by flipping pairs of variables"
#define HEUR_DISPCHAR         'B'
#define HEUR_PRIORITY         -20100
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1

#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

/* default parameter values */
#define DEFAULT_INTOPT                FALSE /**< optional integer optimization is applied by default */
#define DEFAULT_WAITINGNODES              0 /**< default number of nodes to wait after current best solution before calling heuristic */
#define DEFAULT_MATCHINGRATE            0.5 /**< default percentage by which two variables have to match in their LP-row set to be
                                             *   associated as pair by heuristic */
#define DEFAULT_MAXNSLAVES              199 /**< default number of slave candidates for a master variable */
#define DEFAULT_ARRAYSIZE                10 /**< the default array size for temporary arrays */
#define DEFAULT_RANDSEED                 37 /**< initial random seed */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;       /**< index of last solution for which heuristic was performed */
   SCIP_Real             matchingrate;       /**< percentage by which two variables have have to match in their LP-row
                                              *   set to be associated as pair by heuristic */
   SCIP_VAR**            binvars;            /**< Array of binary variables which are sorted with respect to their occurrence
                                              *   in the LP-rows */
   int                   nbinvars;           /**< number of binary variables stored in heuristic array */
   int                   waitingnodes;       /**< user parameter to determine number of nodes to wait after last best solution
                                              *   before calling heuristic   */
   SCIP_Bool             presolved;          /**< flag to indicate whether presolving has already been executed */
   int*                  binblockstart;      /**< array to store the start indices of each binary block */
   int*                  binblockend;        /**< array to store the end indices of each binary block */
   int                   nbinblocks;         /**< number of blocks */

   /* integer variable twoopt data */
   SCIP_Bool             intopt;             /**< parameter to determine if integer 2-opt should be applied */
   SCIP_VAR**            intvars;            /**< array to store the integer variables in non-decreasing order
                                              *   with respect to their objective coefficient */
   int                   nintvars;           /**< the number of integer variables stored in array intvars */
   int*                  intblockstart;      /**< array to store the start indices of each binary block */
   int*                  intblockend;        /**< array to store the end indices of each binary block */
   int                   nintblocks;         /**< number of blocks */

   SCIP_Bool             execute;            /**< has presolveTwoOpt detected necessary structure for execution of heuristic? */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   int                   maxnslaves;         /**< delimits the maximum number of slave candidates for a master variable */

#ifdef SCIP_STATISTIC
   /* statistics */
   int                   ntotalbinvars;      /**< total number of binary variables over all runs */
   int                   ntotalintvars;      /**< total number of Integer variables over all runs */
   int                   nruns;              /**< counts the number of runs, i.e. the number of initialized
                                              *   branch and bound processes */
   int                   maxbinblocksize;    /**< maximum size of a binary block */
   int                   maxintblocksize;    /**< maximum size of an integer block */
   int                   binnblockvars;      /**< number of binary variables that appear in blocks  */
   int                   binnblocks;         /**< number of blocks with at least two variables */
   int                   intnblockvars;      /**< number of Integer variables that appear in blocks  */
   int                   intnblocks;         /**< number of blocks with at least two variables */
   int                   binnexchanges;      /**< number of executed changes of binary solution values leading to
                                              *   improvement in objective function */
   int                   intnexchanges;      /**< number of executed changes of Integer solution values leading to improvement in
                                              *   objective function */
#endif
};

/** indicator for optimizing for binaries or integer variables */
enum Opttype
{
   OPTTYPE_BINARY = 1,
   OPTTYPE_INTEGER = 2
};
typedef enum Opttype OPTTYPE;

/** indicator for direction of shifting variables */
enum Direction
{
   DIRECTION_UP = 1,
   DIRECTION_DOWN = -1,
   DIRECTION_NONE = 0
};
typedef enum Direction DIRECTION;

/*
 * Local methods
 */

/** Tries to switch the values of two binary or integer variables and checks feasibility with respect to the LP.
 *
 *  @todo Adapt method not to copy entire activities array, but only the relevant region.
 */
static
SCIP_RETCODE shiftValues(
   SCIP*                 scip,               /**< scip instance */
   SCIP_VAR*             master,             /**< first variable of variable pair */
   SCIP_VAR*             slave,              /**< second variable of pair */
   SCIP_Real             mastersolval,       /**< current value of variable1 in solution */
   DIRECTION             masterdir,          /**< the direction into which the master variable has to be shifted */
   SCIP_Real             slavesolval,        /**< current value of variable2 in solution */
   DIRECTION             slavedir,           /**< the direction into which the slave variable has to be shifted */
   SCIP_Real             shiftval,           /**< the value that variables should be shifted by */
   SCIP_Real*            activities,         /**< the LP-row activities */
   int                   nrows,              /**< size of activities array */
   SCIP_Bool*            feasible            /**< set to true if method has successfully switched the variable values */
   )
{  /*lint --e{715}*/
   SCIP_COL* col;
   SCIP_ROW** masterrows;
   SCIP_ROW** slaverows;
   SCIP_Real* mastercolvals;
   SCIP_Real* slavecolvals;
   int ncolmasterrows;
   int ncolslaverows;
   int i;
   int j;

   assert(scip != NULL);
   assert(master != NULL);
   assert(slave != NULL);
   assert(activities != NULL);
   assert(SCIPisFeasGT(scip, shiftval, 0.0));

   assert(SCIPisFeasGE(scip, mastersolval + (int)masterdir * shiftval, SCIPvarGetLbGlobal(master)));
   assert(SCIPisFeasLE(scip, mastersolval + (int)masterdir * shiftval, SCIPvarGetUbGlobal(master)));

   assert(SCIPisFeasGE(scip, slavesolval + (int)slavedir * shiftval, SCIPvarGetLbGlobal(slave)));
   assert(SCIPisFeasLE(scip, slavesolval + (int)slavedir * shiftval, SCIPvarGetUbGlobal(slave)));

   /* get variable specific rows and coefficients for both master and slave. */
   col = SCIPvarGetCol(master);
   masterrows = SCIPcolGetRows(col);
   mastercolvals = SCIPcolGetVals(col);
   ncolmasterrows = SCIPcolGetNNonz(col);
   assert(ncolmasterrows == 0 || masterrows != NULL);

   col = SCIPvarGetCol(slave);
   slaverows = SCIPcolGetRows(col);
   slavecolvals = SCIPcolGetVals(col);
   ncolslaverows = SCIPcolGetNNonz(col);
   assert(ncolslaverows == 0 || slaverows != NULL);

   /* update the activities of the LP rows of the master variable */
   for( i = 0; i < ncolmasterrows && SCIProwGetLPPos(masterrows[i]) >= 0; ++i )
   {
      int rowpos;

      rowpos = SCIProwGetLPPos(masterrows[i]);
      assert(rowpos < nrows);

      /* skip local rows */
      if( rowpos >= 0 && ! SCIProwIsLocal(masterrows[i]) )
         activities[rowpos] += mastercolvals[i] * (int)masterdir * shiftval;
   }

   /* update the activities of the LP rows of the slave variable */
   for( j = 0; j < ncolslaverows && SCIProwGetLPPos(slaverows[j]) >= 0; ++j )
   {
      int rowpos;

      rowpos = SCIProwGetLPPos(slaverows[j]);
      assert(rowpos < nrows);

      /* skip local rows */
      if( rowpos >= 0 && ! SCIProwIsLocal(slaverows[j]) )
      {
         activities[rowpos] += slavecolvals[j] * (int)slavedir * shiftval;
         assert(SCIPisFeasGE(scip, activities[rowpos], SCIProwGetLhs(slaverows[j])));
         assert(SCIPisFeasLE(scip, activities[rowpos], SCIProwGetRhs(slaverows[j])));
      }
   }

   /* in debug mode, the master rows are checked for feasibility which should be granted by the
    * decision for a shift value */
#ifndef NDEBUG
   for( i = 0; i < ncolmasterrows && SCIProwGetLPPos(masterrows[i]) >= 0; ++i )
   {
      /* local rows can be skipped */
      if( SCIProwIsLocal(masterrows[i]) )
         continue;

      assert(SCIPisFeasGE(scip, activities[SCIProwGetLPPos(masterrows[i])], SCIProwGetLhs(masterrows[i])));
      assert(SCIPisFeasLE(scip, activities[SCIProwGetLPPos(masterrows[i])], SCIProwGetRhs(masterrows[i])));
   }
#endif

   *feasible = TRUE;

   return SCIP_OKAY;
}

/** Compare two variables with respect to their columns.
 *
 *  Columns are treated as {0,1} vector, where every nonzero entry is treated as '1', and compared to each other
 *  lexicographically. I.e. var1 is < var2 if the corresponding column of var2 has the smaller single nonzero index of
 *  the two columns.  This comparison costs O(constraints) in the worst case
 */
static
int varColCompare(
   SCIP_VAR*             var1,               /**< left argument of comparison */
   SCIP_VAR*             var2                /**< right argument of comparison */
   )
{
   SCIP_COL* col1;
   SCIP_COL* col2;
   SCIP_ROW** rows1;
   SCIP_ROW** rows2;
   int nnonzeros1;
   int nnonzeros2;
   int i;

   assert(var1 != NULL);
   assert(var2 != NULL);

   /* get the necessary row and column data */
   col1 = SCIPvarGetCol(var1);
   col2 = SCIPvarGetCol(var2);
   rows1 = SCIPcolGetRows(col1);
   rows2 = SCIPcolGetRows(col2);
   nnonzeros1 = SCIPcolGetNNonz(col1);
   nnonzeros2 = SCIPcolGetNNonz(col2);

   assert(nnonzeros1 == 0 || rows1 != NULL);
   assert(nnonzeros2 == 0 || rows2 != NULL);

   /* loop over the rows, stopped as soon as they differ in one index,
    * or if counter reaches the end of a variables row set */
   for( i = 0; i < nnonzeros1 && i < nnonzeros2; ++i )
   {
      if( SCIProwGetIndex(rows1[i]) != SCIProwGetIndex(rows2[i]) )
         return SCIProwGetIndex(rows1[i]) - SCIProwGetIndex(rows2[i]);
   }

   /* loop is finished, without differing in one of common row indices, due to loop invariant
    * variable i reached either nnonzeros1 or nnonzeros2 or both.
    * one can easily check that the difference of these two numbers always has the desired sign for comparison. */
   return nnonzeros2 - nnonzeros1 ;
}

/** implements a comparator to compare two variables with respect to their column entries */
static
SCIP_DECL_SORTPTRCOMP(SCIPvarcolComp)
{
   return varColCompare((SCIP_VAR*) elem1, (SCIP_VAR*) elem2);
}

/** checks if two given variables are contained in common LP rows,
 *  returns true if variables share the necessary percentage (matchingrate) of rows.
 */
static
SCIP_Bool checkConstraintMatching(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Real             matchingrate        /**< determines the ratio of shared LP rows compared to the total number of
                                              *   LP-rows each variable appears in */
   )
{
   SCIP_COL* col1;
   SCIP_COL* col2;
   SCIP_ROW** rows1;
   SCIP_ROW** rows2;
   int nnonzeros1;
   int nnonzeros2;
   int i;
   int j;
   int nrows1not2;                           /* the number of LP-rows of variable 1 which variable 2 doesn't appear in */
   int nrows2not1;                           /* vice versa */
   int nrowmaximum;
   int nrowabs;

   assert(var1 != NULL);
   assert(var2 != NULL);

   /* get the necessary row and column data */
   col1 = SCIPvarGetCol(var1);
   col2 = SCIPvarGetCol(var2);
   rows1 = SCIPcolGetRows(col1);
   rows2 = SCIPcolGetRows(col2);
   nnonzeros1 = SCIPcolGetNNonz(col1);
   nnonzeros2 = SCIPcolGetNNonz(col2);

   assert(nnonzeros1 == 0 || rows1 != NULL);
   assert(nnonzeros2 == 0 || rows2 != NULL);

   if( nnonzeros1 == 0 && nnonzeros2 == 0 )
      return TRUE;

   /* initialize the counters for the number of rows not shared. */
   nrowmaximum = MAX(nnonzeros1, nnonzeros2);

   nrowabs = ABS(nnonzeros1 - nnonzeros2);
   nrows1not2 = nrowmaximum - nnonzeros2;
   nrows2not1 = nrowmaximum - nnonzeros1;

   /* if the numbers of nonzero rows differs too much, w.r.t.matching ratio, the more expensive check over the rows
    * doesn't have to be applied anymore because the counters for not shared rows can only increase.
    */
   assert(nrowmaximum > 0);

   if( (nrowmaximum - nrowabs) / (SCIP_Real) nrowmaximum < matchingrate )
      return FALSE;

   i = 0;
   j = 0;

   /* loop over all rows and determine number of non-shared rows */
   while( i < nnonzeros1 && j < nnonzeros2 )
   {
      /* variables share a common row */
      if( SCIProwGetIndex(rows1[i]) == SCIProwGetIndex(rows2[j]) )
      {
         ++i;
         ++j;
      }
      /* variable 1 appears in rows1[i], variable 2 doesn't */
      else if( SCIProwGetIndex(rows1[i]) < SCIProwGetIndex(rows2[j]) )
      {
         ++i;
         ++nrows1not2;
      }
      /* variable 2 appears in rows2[j], variable 1 doesn't */
      else
      {
         ++j;
         ++nrows2not1;
      }
   }

   /* now apply the ratio based comparison, that is if the ratio of shared rows is greater or equal the matching rate
    * for each variable */
   return ( SCIPisFeasLE(scip, matchingrate, (nnonzeros1 - nrows1not2) / (SCIP_Real)(nnonzeros1)) ||
      SCIPisFeasLE(scip, matchingrate, (nnonzeros2 - nrows2not1) / (SCIP_Real)(nnonzeros2)) );  /*lint !e795 */
}

/** Determines a bound by which the absolute solution value of two integer variables can be shifted at most.
 *
 *  The criterion is the maintenance of feasibility of any global LP row.
 *  The first implementation only considers shifting proportion 1:1, i.e. if master value is shifted by a certain
 *  integer value k downwards, the value of slave is simultaneously shifted by k upwards.
 */
static
SCIP_Real determineBound(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_SOL*             sol,                /**< current incumbent */
   SCIP_VAR*             master,             /**< current master variable */
   DIRECTION             masterdirection,    /**< the shifting direction of the master variable */
   SCIP_VAR*             slave,              /**< slave variable with same LP_row set as master variable */
   DIRECTION             slavedirection,     /**< the shifting direction of the slave variable */
   SCIP_Real*            activities,         /**< array of LP row activities */
   int                   nrows               /**< the number of rows in LP and the size of the activities array */
   )
{  /*lint --e{715}*/
   SCIP_Real masterbound;
   SCIP_Real slavebound;
   SCIP_Real bound;

   SCIP_COL* col;
   SCIP_ROW** slaverows;
   SCIP_ROW** masterrows;
   SCIP_Real* mastercolvals;
   SCIP_Real* slavecolvals;
   int nslaverows;
   int nmasterrows;
   int i;
   int j;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(master != NULL);
   assert(slave != NULL);
   assert(SCIPvarIsIntegral(master) && SCIPvarIsIntegral(slave));
   assert(masterdirection == DIRECTION_UP || masterdirection == DIRECTION_DOWN);
   assert(slavedirection == DIRECTION_UP || slavedirection == DIRECTION_DOWN);

   /* determine the trivial variable bounds for shift */
   if( masterdirection == DIRECTION_UP )
      masterbound = SCIPvarGetUbGlobal(master) - SCIPgetSolVal(scip, sol, master);
   else
      masterbound = SCIPgetSolVal(scip, sol, master) - SCIPvarGetLbGlobal(master);

   if( slavedirection == DIRECTION_UP )
      slavebound = SCIPvarGetUbGlobal(slave) - SCIPgetSolVal(scip, sol, slave);
   else
      slavebound = SCIPgetSolVal(scip, sol, slave) - SCIPvarGetLbGlobal(slave);

   bound = MIN(slavebound, masterbound);
   assert(!SCIPisInfinity(scip,bound));

   if( bound < 0.5 )
      return 0.0;

   /* get the necessary row and and column data for each variable */
   col = SCIPvarGetCol(slave);
   slaverows = SCIPcolGetRows(col);
   slavecolvals = SCIPcolGetVals(col);
   nslaverows = SCIPcolGetNNonz(col);

   col = SCIPvarGetCol(master);
   mastercolvals = SCIPcolGetVals(col);
   masterrows = SCIPcolGetRows(col);
   nmasterrows = SCIPcolGetNNonz(col);

   assert(nslaverows == 0 || slavecolvals != NULL);
   assert(nmasterrows == 0 || mastercolvals != NULL);

   SCIPdebugMsg(scip, "  Master: %s with direction %d and %d rows, Slave: %s with direction %d and %d rows \n", SCIPvarGetName(master),
      (int)masterdirection, nmasterrows, SCIPvarGetName(slave), (int)slavedirection, nslaverows);

   /* loop over all LP rows and determine the maximum integer bound by which both variables
    * can be shifted without loss of feasibility
    */
   i = 0;
   j = 0;
   while( (i < nslaverows || j < nmasterrows) && SCIPisPositive(scip, bound) )
   {
      SCIP_ROW* row;
      SCIP_Real effect;
      SCIP_Real rhs;
      SCIP_Real lhs;
      SCIP_Real activity;
      int rowpos;
      int masterindex;
      int slaveindex;
      SCIP_Bool slaveincrement;
      SCIP_Bool masterincrement;

      /* check if one pointer already reached the end of the respective array */
      if( i < nslaverows && SCIProwGetLPPos(slaverows[i]) == -1 )
      {
         SCIPdebugMsg(scip, "  Slaverow %s is not in LP (i=%d, j=%d)\n", SCIProwGetName(slaverows[i]), i, j);
         i = nslaverows;
         continue;
      }
      if( j < nmasterrows && SCIProwGetLPPos(masterrows[j]) == -1 )
      {
         SCIPdebugMsg(scip, "  Masterrow %s is not in LP (i=%d, j=%d) \n", SCIProwGetName(masterrows[j]), i, j);
         j = nmasterrows;
         continue;
      }

      slaveincrement = FALSE;
      /* If one counter has already reached its limit, assign a huge number to the corresponding
       * row index to simulate an always greater row position. */
      if( i < nslaverows )
         slaveindex = SCIProwGetIndex(slaverows[i]);
      else
         slaveindex = INT_MAX;

      if( j < nmasterrows )
         masterindex = SCIProwGetIndex(masterrows[j]);
      else
         masterindex = INT_MAX;

      assert(0 <= slaveindex && 0 <= masterindex);

      assert(slaveindex < INT_MAX || masterindex < INT_MAX);

      /* the current row is the one with the smaller index */
      if( slaveindex <= masterindex )
      {
         rowpos = SCIProwGetLPPos(slaverows[i]);
         row = slaverows[i];
         slaveincrement = TRUE;
         masterincrement = (slaveindex == masterindex);
      }
      else
      {
         assert(j < nmasterrows);

         rowpos = SCIProwGetLPPos(masterrows[j]);
         row = masterrows[j];
         masterincrement = TRUE;
      }
      assert(row != NULL);

      /* local rows can be skipped */
      if( !SCIProwIsLocal(row) )
      {
         /* effect is the effect on the row activity by shifting the variables by 1 in the respective directions */
         effect = 0.0;
         if( slaveindex <= masterindex )
            effect += (slavecolvals[i] * (int)slavedirection);
         if( masterindex <= slaveindex )
            effect += (mastercolvals[j] * (int)masterdirection);

         /* get information about the current row */
         if( rowpos >= 0 && !SCIPisFeasZero(scip, effect) )
         {
            /* effect does not equal zero, the bound is determined as minimum positive integer such that
             * feasibility of this constraint is maintained.
             * if constraint is an equality constraint, activity and lhs/rhs should be feasibly equal, which
             * will cause the method to return zero.
             */
            assert(rowpos < nrows);

            activity = activities[rowpos];
            rhs = SCIProwGetRhs(row);
            lhs = SCIProwGetLhs(row);

            /* if the row is an equation, abort because of effect being nonzero */
            if( SCIPisFeasEQ(scip, lhs, rhs) )
               return 0.0;

            assert(SCIPisFeasLE(scip, lhs, activity) && SCIPisFeasLE(scip, activity, rhs));

            SCIPdebugMsg(scip, "   %g <= %g <= %g, bound = %g, effect = %g (%g * %d + %g * %d) (i=%d,j=%d)\n", lhs, activity, rhs, bound, effect,
               slaveindex <= masterindex ? slavecolvals[i] : 0.0, (int)slavedirection, masterindex <= slaveindex ? mastercolvals[j] : 0.0,
               (int)masterdirection, i, j);

            /* if the row has a left hand side, ensure that shifting preserves feasibility of this ">="-constraint */
            if( !SCIPisInfinity(scip, -lhs) && SCIPisFeasLT(scip, activity + (effect * bound), lhs) )
            {
               SCIP_Real newval;

               assert(SCIPisNegative(scip, effect));

               newval = SCIPfeasFloor(scip, (lhs - activity)/effect); /*lint !e414*/
               bound = MIN(bound - 1.0, newval);
            }

            /* if the row has an upper bound, ensure that shifting preserves feasibility of this "<="-constraint */
            if( !SCIPisInfinity(scip, rhs) && SCIPisFeasGT(scip, activity + (effect * bound), rhs) )
            {
               SCIP_Real newval;

               assert(SCIPisPositive(scip, effect));

               newval = SCIPfeasFloor(scip, (rhs - activity)/effect); /*lint !e414*/
               bound = MIN(bound - 1.0, newval);
            }

            assert(SCIPisFeasLE(scip, lhs, activity + effect * bound) && SCIPisFeasGE(scip, rhs, activity + effect * bound));
            assert(SCIPisFeasIntegral(scip, bound));
         }
         else
         {
            SCIPdebugMsg(scip, "  Zero effect on row %s, effect %g, master coeff: %g slave coeff: %g (i=%d, j=%d)\n",
               SCIProwGetName(row), effect, mastercolvals[j], slavecolvals[i], i, j);
         }
      }
      else
      {
         SCIPdebugMsg(scip, "  Row %s is local.\n", SCIProwGetName(row));
      }

      /* increase the counters which belong to the corresponding row. Both counters are increased by
       * 1 iff rowpos1 <= rowpos2 <= rowpos1 */
      if( slaveincrement )
         ++i;
      if( masterincrement )
         ++j;
   }

   /* due to numerical reasons, bound can be negative. A variable shift by a negative bound is not desired by
    * by the heuristic -> Change the return value to zero */
   if( !SCIPisPositive(scip, bound) )
      bound = 0.0;

   return bound;
}

/** Disposes variable with no heuristic relevancy, e.g., due to a fixed solution value, from its neighborhood block.
 *
 *  The affected neighborhood block is reduced by 1.
 */
static
void disposeVariable(
   SCIP_VAR**            vars,               /**< problem variables */
   int*                  blockend,           /**< contains end index of block */
   int                   varindex            /**< variable index */
   )
{
   assert(blockend != NULL);
   assert(varindex <= *blockend);

   vars[varindex] = vars[*blockend];
   --(*blockend);
}

/** realizes the presolve independently from type of variables it's applied to */
static
SCIP_RETCODE innerPresolve(
   SCIP*                 scip,               /**< current scip */
   SCIP_VAR**            vars,               /**< problem vars */
   SCIP_VAR***           varspointer,        /**< pointer to heuristic specific variable memory */
   int                   nvars,              /**< the number of variables */
   int*                  nblocks,            /**< pointer to store the number of detected blocks */
   int*                  maxblocksize,       /**< maximum size of a block */
   int*                  nblockvars,         /**< pointer to store the number of block variables */
   int**                 blockstart,         /**< pointer to store the array of block start indices */
   int**                 blockend,           /**< pointer to store the array of block end indices */
   SCIP_HEUR*            heur,               /**< the heuristic */
   SCIP_HEURDATA*        heurdata            /**< the heuristic data */
   )
{
   int v;
   int startindex;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 2);
   assert(nblocks != NULL);
   assert(nblockvars != NULL);
   assert(blockstart != NULL);
   assert(blockend != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);

   /* allocate the heuristic specific variables */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, varspointer, vars, nvars));

   /* sort the variables with respect to their columns */
   SCIPsortPtr((void**)(*varspointer), SCIPvarcolComp, nvars);

   /* start determining blocks, i.e. a set of at least two variables which share most of their row set.
    * If there is none, heuristic does not need to be executed.
    */
   startindex = 0;
   *nblocks = 0;
   *nblockvars = 0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, blockstart, nvars/2) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, blockend, nvars/2) );

   /* loop over variables and compare neighbors */
   for( v = 1; v < nvars; ++v )
   {
      if( !checkConstraintMatching(scip, (*varspointer)[startindex], (*varspointer)[v], heurdata->matchingrate) )
      {
         /* current block has its last variable at position v-1. If v differs from startindex by at least 2,
          * a block is detected. Update the data correspondingly */
         if( v - startindex >= 2 )
         {
            assert(*nblocks < nvars/2);
            (*nblockvars) += v - startindex;
            (*maxblocksize) = MAX((*maxblocksize), v - startindex);
            (*blockstart)[*nblocks] = startindex;
            (*blockend)[*nblocks] = v - 1;
            (*nblocks)++;
         }
         startindex = v;
      }
      else if( v == nvars - 1 && v - startindex >= 1 )
      {
         assert(*nblocks < nvars/2);
         (*nblockvars) += v - startindex + 1;
         (*maxblocksize) = MAX((*maxblocksize), v - startindex + 1);
         (*blockstart)[*nblocks] = startindex;
         (*blockend)[*nblocks] = v;
         (*nblocks)++;
      }
   }

   /* reallocate memory with respect to the number of found blocks; if there were none, free the memory */
   if( *nblocks > 0 )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, blockstart, nvars/2, *nblocks) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, blockend, nvars/2, *nblocks) );
   }
   else
   {
      SCIPfreeBlockMemoryArray(scip, blockstart, nvars/2);
      SCIPfreeBlockMemoryArray(scip, blockend, nvars/2);

      *blockstart = NULL;
      *blockend = NULL;
   }

   return SCIP_OKAY;
}

/** initializes the required structures for execution of heuristic.
 *
 *  If objective coefficient functions are not all equal, each Binary and Integer variables are sorted
 *  into heuristic-specific arrays with respect to their lexicographical column order,
 *  where every zero in a column is interpreted as zero and every nonzero as '1'.
 *  After the sorting, the variables are compared with respect to user parameter matchingrate and
 *  the heuristic specific blocks are determined.
 */
static
SCIP_RETCODE presolveTwoOpt(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata            /**< the heuristic data */
   )
{
   int nbinvars;
   int nintvars;
   int nvars;
   SCIP_VAR** vars;
   int nbinblockvars = 0;
   int nintblockvars;
   int maxbinblocksize = 0;
   int maxintblocksize;

   assert(scip != NULL);
   assert(heurdata != NULL);

   /* ensure that method is not executed if presolving was already applied once in current branch and bound process */
   if( heurdata->presolved )
      return SCIP_OKAY;

   /* get necessary variable information, i.e. number of binary and integer variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* if number of binary problem variables exceeds 2, they are subject to 2-optimization algorithm, hence heuristic
    * calls innerPresolve method to detect necessary structures. */
   if( nbinvars >= 2 )
   {
      SCIP_CALL( innerPresolve(scip, vars, &(heurdata->binvars), nbinvars, &(heurdata->nbinblocks), &maxbinblocksize,
            &nbinblockvars, &(heurdata->binblockstart), &(heurdata->binblockend), heur, heurdata) );
   }

   heurdata->nbinvars = nbinvars;
   heurdata->execute = nbinvars > 1 && heurdata->nbinblocks > 0;

#ifdef SCIP_STATISTIC
   /* update statistics */
   heurdata->binnblocks += (heurdata->nbinblocks);
   heurdata->binnblockvars += nbinblockvars;
   heurdata->ntotalbinvars += nbinvars;
   heurdata->maxbinblocksize = MAX(maxbinblocksize, heurdata->maxbinblocksize);

   SCIPstatisticMessage("   Twoopt BINARY presolving finished with <%d> blocks, <%d> block variables \n",
      heurdata->nbinblocks, nbinblockvars);
#endif

   if( heurdata->intopt && nintvars > 1 )
   {
      SCIP_CALL( innerPresolve(scip, &(vars[nbinvars]), &(heurdata->intvars), nintvars, &(heurdata->nintblocks), &maxintblocksize,
            &nintblockvars, &(heurdata->intblockstart), &(heurdata->intblockend),
            heur, heurdata) );

      heurdata->execute = heurdata->execute || heurdata->nintblocks > 0;

#ifdef SCIP_STATISTIC
      /* update statistics */
      heurdata->intnblocks += heurdata->nintblocks;
      heurdata->intnblockvars += nintblockvars;
      heurdata->ntotalintvars += nintvars;
      heurdata->maxintblocksize = MAX(maxintblocksize, heurdata->maxintblocksize);
     SCIPstatisticMessage("   Twoopt Integer presolving finished with <%d> blocks, <%d> block variables \n",
         heurdata->nintblocks, nintblockvars);
     SCIPstatisticMessage("   INTEGER coefficients are all equal \n");
#endif
   }
   heurdata->nintvars = nintvars;

   /* presolving is finished, heuristic data is updated*/
   heurdata->presolved = TRUE;
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTwoopt)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTwoopt(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTwoopt)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitTwoopt)
{
   SCIP_HEURDATA* heurdata;
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* heuristic has not run yet, all heuristic specific data is set to initial values */
   heurdata->nbinvars = 0;
   heurdata->nintvars = 0;
   heurdata->lastsolindex = -1;
   heurdata->presolved = FALSE;
   heurdata->nbinblocks = 0;
   heurdata->nintblocks = 0;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED) );

#ifdef SCIP_STATISTIC
   /* initialize statistics */
   heurdata->binnexchanges = 0;
   heurdata->intnexchanges = 0;
   heurdata->binnblockvars = 0;
   heurdata->intnblockvars = 0;
   heurdata->binnblocks = 0;
   heurdata->intnblocks = 0;

   heurdata->maxbinblocksize = 0;
   heurdata->maxintblocksize = 0;

   heurdata->ntotalbinvars = 0;
   heurdata->ntotalintvars = 0;
   heurdata->nruns = 0;
#endif

   /* all pointers are initially set to NULL. Since presolving
    * of the heuristic requires a lot of calculation time on some instances,
    * but might not be needed e.g. if problem is infeasible, presolving is applied
    * when heuristic is executed for the first time
    */
   heurdata->binvars = NULL;
   heurdata->intvars = NULL;
   heurdata->binblockstart = NULL;
   heurdata->binblockend = NULL;
   heurdata->intblockstart = NULL;
   heurdata->intblockend = NULL;

   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}

/* Realizes the 2-optimization algorithm, which tries to improve incumbent solution
 * by shifting pairs of variables which share a common row set.
 */
static
SCIP_RETCODE optimize(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_SOL*             worksol,            /**< working solution */
   SCIP_VAR**            vars,               /**< binary or integer variables */
   int*                  blockstart,         /**< contains start indices of blocks */
   int*                  blockend,           /**< contains end indices of blocks */
   int                   nblocks,            /**< the number of blocks */
   OPTTYPE               opttype,            /**< are binaries or integers optimized */
   SCIP_Real*            activities,         /**< the LP-row activities */
   int                   nrows,              /**< the number of LP rows */
   SCIP_Bool*            improvement,        /**< was there a successful shift? */
   SCIP_Bool*            varboundserr,       /**< has the current incumbent already been cut off */
   SCIP_HEURDATA*        heurdata            /**< the heuristic data */
   )
{  /*lint --e{715}*/
   int b;
   SCIP_Real* objchanges;
   SCIP_VAR** bestmasters;
   SCIP_VAR** bestslaves;
   int* bestdirections;
   int arraysize;
   int npairs = 0;

   assert(scip != NULL);
   assert(nblocks > 0);
   assert(blockstart != NULL && blockend != NULL);
   assert(varboundserr != NULL);
   assert(activities != NULL);
   assert(worksol != NULL);
   assert(improvement != NULL);

   *varboundserr = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &bestmasters, DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestslaves, DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objchanges, DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestdirections, DEFAULT_ARRAYSIZE) );
   arraysize = DEFAULT_ARRAYSIZE;

   /* iterate over blocks */
   for( b = 0; b < nblocks; ++b )
   {
      int m;
      int blocklen;

      blocklen = blockend[b] - blockstart[b] + 1;

      /* iterate over variables in current block */
      for( m = 0; m < blocklen; ++m )
      {
         /* determine the new master variable for heuristic's optimization method */
         SCIP_VAR* master;
         SCIP_Real masterobj;
         SCIP_Real mastersolval;
         SCIP_Real bestimprovement;
         SCIP_Real bestbound;
         int bestslavepos;
         int s;
         int firstslave;
         int nslaves;
         int bestdirection;
         DIRECTION bestmasterdir;
         DIRECTION bestslavedir;

         master = vars[blockstart[b] + m]; /*lint !e679*/
         masterobj = SCIPvarGetObj(master);
         mastersolval = SCIPgetSolVal(scip, worksol, master);

         /* due to cuts or fixings of solution values, worksol might not be feasible w.r.t. its bounds.
          * Exit method in that case. */
         if( SCIPisFeasGT(scip, mastersolval, SCIPvarGetUbGlobal(master)) || SCIPisFeasLT(scip, mastersolval, SCIPvarGetLbGlobal(master)) )
         {
            *varboundserr = TRUE;
            SCIPdebugMsg(scip, "Solution has violated variable bounds for var %s: %g <= %g <= %g \n",
               SCIPvarGetName(master), SCIPvarGetLbGlobal(master), mastersolval, SCIPvarGetUbGlobal(master));
            goto TERMINATE;
         }

         /* if variable has fixed solution value, it is deleted from heuristic array */
         if( SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(master), SCIPvarGetLbGlobal(master)) )
         {
            disposeVariable(vars, &(blockend[b]), blockstart[b] + m);
            --blocklen;
            continue;
         }
         else if( SCIPvarGetStatus(master) != SCIP_VARSTATUS_COLUMN )
            continue;

         assert(SCIPisFeasIntegral(scip, mastersolval));

         assert(opttype == OPTTYPE_INTEGER || (SCIPisFeasLE(scip, mastersolval, 1.0) || SCIPisFeasGE(scip, mastersolval, 0.0)));

         /* initialize the data of the best available shift */
         bestimprovement = 0.0;
         bestslavepos = -1;
         bestbound = 0.0;
         bestmasterdir = DIRECTION_NONE;
         bestslavedir = DIRECTION_NONE;
         bestdirection = -1;

         /* in blocks with more than heurdata->maxnslaves variables, a slave candidate region is chosen */
         if( heurdata->maxnslaves >= 0 && blocklen > heurdata->maxnslaves )
            firstslave = SCIPrandomGetInt(heurdata->randnumgen, blockstart[b] + m, blockend[b]);
         else
            firstslave = blockstart[b] + m + 1;

         nslaves = MIN((heurdata->maxnslaves == -1 ? INT_MAX : heurdata->maxnslaves), blocklen);

         /* Loop over block and determine a slave shift candidate for master variable.
          * If more than one candidate is available, choose the shift which improves objective function
          * the most. */
         for( s = 0; s < nslaves; ++s )
         {
            SCIP_VAR* slave;
            SCIP_Real slaveobj;
            SCIP_Real slavesolval;
            SCIP_Real changedobj;
            SCIP_Real diffdirbound;
            SCIP_Real equaldirbound;
            int directions;
            int slaveindex;

            slaveindex = (firstslave + s - blockstart[b]) % blocklen;
            slaveindex += blockstart[b];

            /* in case of a small block, we do not want to try possible pairings twice */
            if( (blocklen <= heurdata->maxnslaves || heurdata->maxnslaves == -1) && slaveindex < blockstart[b] + m )
               break;
            /* master and slave should not be the same variable */
            if( slaveindex == blockstart[b] + m )
               continue;

            /* get the next slave variable */
            slave = vars[slaveindex];
            slaveobj = SCIPvarGetObj(slave);
            slavesolval = SCIPgetSolVal(scip, worksol, slave);
            changedobj = 0.0;

            assert(SCIPvarGetType(master) == SCIPvarGetType(slave));
            assert(SCIPisFeasIntegral(scip, slavesolval));
            assert(opttype == OPTTYPE_INTEGER || (SCIPisFeasLE(scip, mastersolval, 1.0) || SCIPisFeasGE(scip, mastersolval, 0.0)));

            /* solution is not feasible w.r.t. the variable bounds, stop optimization in this case */
            if( SCIPisFeasGT(scip, slavesolval, SCIPvarGetUbGlobal(slave)) || SCIPisFeasLT(scip, slavesolval, SCIPvarGetLbGlobal(slave)) )
            {
               *varboundserr = TRUE;
               SCIPdebugMsg(scip, "Solution has violated variable bounds for var %s: %g <= %g <= %g \n",
                  SCIPvarGetName(slave), SCIPvarGetLbGlobal(slave), slavesolval, SCIPvarGetUbGlobal(slave));
               goto TERMINATE;
            }

            /* if solution value of the variable is fixed, delete it from the remaining candidates in the block */
            if( SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(slave), SCIPvarGetLbGlobal(slave) ) )
            {
               disposeVariable(vars, &(blockend[b]), slaveindex);
               --blocklen;
               continue;
            }
            else if( SCIPvarGetStatus(master) != SCIP_VARSTATUS_COLUMN )
               continue;

            /* determine the shifting direction to improve the objective function */
            /* assert(SCIPisFeasGT(scip, masterobj, slaveobj)); */

            /* The heuristic chooses the shifting direction and the corresponding maximum nonnegative
             * integer shift value for the two variables which preserves feasibility and improves
             * the objective function value. */
            directions = -1;
            diffdirbound = 0.0;
            equaldirbound = 0.0;

            if( SCIPisFeasLT(scip, masterobj - slaveobj, 0.0) )
            {
               diffdirbound = determineBound(scip, worksol, master, DIRECTION_UP,  slave, DIRECTION_DOWN, activities, nrows);
               directions = 2;
               /* the improvement of objective function is calculated */
               changedobj = (masterobj - slaveobj) * diffdirbound;
            }
            else if( SCIPisFeasGT(scip, masterobj - slaveobj, 0.0) )
            {
               diffdirbound = determineBound(scip, worksol, master, DIRECTION_DOWN,  slave, DIRECTION_UP, activities, nrows);
               directions = 1;
               changedobj = (slaveobj - masterobj) * diffdirbound;
            }

            if( SCIPisFeasLT(scip, masterobj + slaveobj, 0.0) )
            {
               equaldirbound = determineBound(scip, worksol, master, DIRECTION_UP,  slave, DIRECTION_UP, activities, nrows);
               if( SCIPisFeasLT(scip, (slaveobj + masterobj) * equaldirbound, changedobj) )
               {
                  changedobj = (slaveobj + masterobj) * equaldirbound;
                  directions = 3;
               }
            }
            else if( SCIPisFeasGT(scip, masterobj + slaveobj, 0.0) )
            {
               equaldirbound = determineBound(scip, worksol, master, DIRECTION_DOWN,  slave, DIRECTION_DOWN, activities, nrows);
               if( SCIPisFeasLT(scip, -(slaveobj + masterobj) * equaldirbound, changedobj) )
               {
                  changedobj = -(slaveobj + masterobj) * equaldirbound;
                  directions = 0;
               }
            }
            assert(SCIPisFeasIntegral(scip, equaldirbound));
            assert(SCIPisFeasIntegral(scip, diffdirbound));
            assert(SCIPisFeasGE(scip, equaldirbound, 0.0));
            assert(SCIPisFeasGE(scip, diffdirbound, 0.0));

            /* choose the candidate which improves the objective function the most */
            if( (SCIPisFeasGT(scip, equaldirbound, 0.0) || SCIPisFeasGT(scip, diffdirbound, 0.0))
               && SCIPisFeasLT(scip, changedobj, bestimprovement) )
            {
               bestimprovement = changedobj;
               bestslavepos = slaveindex;
               bestdirection = directions;

               /* the most promising shift, i.e., the one which can improve the objective
                * the most, is recorded by the integer 'directions'. It is recovered via the use
                * of a binary representation of the four different combinations for the shifting directions
                * of two variables */
               if( directions / 2 == 1 )
                  bestmasterdir = DIRECTION_UP;
               else
                  bestmasterdir = DIRECTION_DOWN;

               if( directions % 2 == 1 )
                  bestslavedir = DIRECTION_UP;
               else
                  bestslavedir = DIRECTION_DOWN;

               if( bestmasterdir == bestslavedir )
                  bestbound = equaldirbound;
               else
                  bestbound = diffdirbound;
            }
         }

         /* choose the most promising candidate, if one exists */
         if( bestslavepos >= 0 )
         {
            if( npairs == arraysize )
            {
               SCIP_CALL( SCIPreallocBufferArray(scip, &bestmasters, 2 * arraysize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &bestslaves, 2 * arraysize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &objchanges, 2 * arraysize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &bestdirections, 2 * arraysize) );
               arraysize = 2 * arraysize;
            }
            assert( npairs < arraysize );

            bestmasters[npairs] = master;
            bestslaves[npairs] = vars[bestslavepos];
            objchanges[npairs] = ((int)bestslavedir * SCIPvarGetObj(bestslaves[npairs])  + (int)bestmasterdir *  masterobj) * bestbound;
            bestdirections[npairs] = bestdirection;

            assert(objchanges[npairs] < 0);

            SCIPdebugMsg(scip, "  Saved candidate pair {%s=%g, %s=%g} with objectives <%g>, <%g> to be set to {%g, %g} %d\n",
               SCIPvarGetName(master), mastersolval, SCIPvarGetName(bestslaves[npairs]), SCIPgetSolVal(scip, worksol, bestslaves[npairs]) ,
               masterobj, SCIPvarGetObj(bestslaves[npairs]), mastersolval + (int)bestmasterdir * bestbound, SCIPgetSolVal(scip, worksol, bestslaves[npairs])
               + (int)bestslavedir * bestbound, bestdirections[npairs]);

            ++npairs;
         }
      }
   }

   if( npairs == 0 )
      goto TERMINATE;

   SCIPsortRealPtrPtrInt(objchanges, (void**)bestmasters, (void**)bestslaves, bestdirections, npairs);

   for( b = 0; b < npairs; ++b )
   {
      SCIP_VAR* master;
      SCIP_VAR* slave;
      SCIP_Real mastersolval;
      SCIP_Real slavesolval;
      SCIP_Real masterobj;
      SCIP_Real slaveobj;
      SCIP_Real bound;
      DIRECTION masterdir;
      DIRECTION slavedir;

      master = bestmasters[b];
      slave = bestslaves[b];
      mastersolval = SCIPgetSolVal(scip, worksol, master);
      slavesolval = SCIPgetSolVal(scip, worksol, slave);
      masterobj  =SCIPvarGetObj(master);
      slaveobj = SCIPvarGetObj(slave);

      assert(0 <= bestdirections[b] && bestdirections[b] < 4);

      if( bestdirections[b] / 2 == 1 )
         masterdir = DIRECTION_UP;
      else
         masterdir = DIRECTION_DOWN;

      if( bestdirections[b] % 2 == 1 )
         slavedir = DIRECTION_UP;
      else
         slavedir = DIRECTION_DOWN;

      bound = determineBound(scip, worksol, master, masterdir, slave, slavedir, activities, nrows);

      if( !SCIPisZero(scip, bound) )
      {
         SCIP_Bool feasible;
#ifndef NDEBUG
         SCIP_Real changedobj;
#endif

         SCIPdebugMsg(scip, "  Promising candidates {%s=%g, %s=%g} with objectives <%g>, <%g> to be set to {%g, %g}\n",
            SCIPvarGetName(master), mastersolval, SCIPvarGetName(slave), slavesolval,
            masterobj, slaveobj, mastersolval + (int)masterdir * bound, slavesolval + (int)slavedir * bound);

#ifndef NDEBUG
         /* the improvement of objective function is calculated */
         changedobj = ((int)slavedir * slaveobj  + (int)masterdir *  masterobj) * bound;
         assert(SCIPisFeasLT(scip, changedobj, 0.0));
#endif

         assert(SCIPvarGetStatus(master) == SCIP_VARSTATUS_COLUMN && SCIPvarGetStatus(slave) == SCIP_VARSTATUS_COLUMN);
         /* try to change the solution values of the variables */
         feasible = FALSE;
         SCIP_CALL( shiftValues(scip, master, slave, mastersolval, masterdir, slavesolval, slavedir, bound,
               activities, nrows, &feasible) );

         if( feasible )
         {
            /* The variables' solution values were successfully shifted and can hence be updated. */
            assert(SCIPisFeasIntegral(scip, mastersolval + ((int)masterdir * bound)));
            assert(SCIPisFeasIntegral(scip, slavesolval + ((int)slavedir * bound)));

            SCIP_CALL( SCIPsetSolVal(scip, worksol, master, mastersolval + (int)masterdir * bound) );
            SCIP_CALL( SCIPsetSolVal(scip, worksol, slave, slavesolval + (int)slavedir * bound) );
            SCIPdebugMsg(scip, "  Feasible shift: <%s>[%g, %g] (obj: %f)  <%f> --> <%f>\n",
               SCIPvarGetName(master), SCIPvarGetLbGlobal(master), SCIPvarGetUbGlobal(master), masterobj, mastersolval, SCIPgetSolVal(scip, worksol, master));
            SCIPdebugMsg(scip, "                  <%s>[%g, %g] (obj: %f)  <%f> --> <%f>\n",
               SCIPvarGetName(slave), SCIPvarGetLbGlobal(slave), SCIPvarGetUbGlobal(slave), slaveobj, slavesolval, SCIPgetSolVal(scip, worksol, slave));

#ifdef SCIP_STATISTIC
            /* update statistics */
            if( opttype == OPTTYPE_BINARY )
               ++(heurdata->binnexchanges);
            else
               ++(heurdata->intnexchanges);
#endif

            *improvement = TRUE;
         }
      }
   }
 TERMINATE:
   SCIPfreeBufferArray(scip, &bestdirections);
   SCIPfreeBufferArray(scip, &objchanges);
   SCIPfreeBufferArray(scip, &bestslaves);
   SCIPfreeBufferArray(scip, &bestmasters);

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitTwoopt)
{
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /*ensure that initialization was successful */
   assert(heurdata->nbinvars <= 1 || heurdata->binvars != NULL);

#ifdef SCIP_STATISTIC
   /* print relevant statistics to console */
   SCIPstatisticMessage(
      "  Twoopt Binary Statistics  :   "
      "%6.2g   %6.2g   %4.2g   %4.0g %6d (blocks/run, variables/run, varpercentage, avg. block size, max block size) \n",
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->binnblocks/(heurdata->nruns),
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->binnblockvars/(heurdata->nruns),
      heurdata->ntotalbinvars == 0 ? 0.0 : (SCIP_Real)heurdata->binnblockvars/(heurdata->ntotalbinvars) * 100.0,
      heurdata->binnblocks == 0 ? 0.0 : heurdata->binnblockvars/(SCIP_Real)(heurdata->binnblocks),
      heurdata->maxbinblocksize);

   SCIPstatisticMessage(
      "   Twoopt Integer statistics :   "
      "%6.2g   %6.2g   %4.2g   %4.0g %6d (blocks/run, variables/run, varpercentage, avg block size, max block size) \n",
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->intnblocks/(heurdata->nruns),
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->intnblockvars/(heurdata->nruns),
      heurdata->ntotalintvars == 0 ? 0.0 : (SCIP_Real)heurdata->intnblockvars/(heurdata->ntotalintvars) * 100.0,
      heurdata->intnblocks == 0 ? 0.0 : heurdata->intnblockvars/(SCIP_Real)(heurdata->intnblocks),
      heurdata->maxintblocksize);

   SCIPstatisticMessage(
      "  Twoopt results            :   "
      "%6d   %6d   %4d   %4.2g  (runs, binary exchanges, Integer shiftings, matching rate)\n",
      heurdata->nruns,
      heurdata->binnexchanges,
      heurdata->intnexchanges,
      heurdata->matchingrate);

   /* set statistics to initial values*/
   heurdata->binnblockvars = 0;
   heurdata->binnblocks = 0;
   heurdata->intnblocks = 0;
   heurdata->intnblockvars = 0;
   heurdata->binnexchanges = 0;
   heurdata->intnexchanges = 0;
#endif

   /* free the allocated memory for the binary variables */
   if( heurdata->binvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->binvars, heurdata->nbinvars);
   }

   if( heurdata->nbinblocks > 0 )
   {
      assert(heurdata->binblockstart != NULL);
      assert(heurdata->binblockend != NULL);

      SCIPfreeBlockMemoryArray(scip, &heurdata->binblockstart, heurdata->nbinblocks);
      SCIPfreeBlockMemoryArray(scip, &heurdata->binblockend, heurdata->nbinblocks);
   }
   heurdata->nbinvars = 0;
   heurdata->nbinblocks = 0;

   if( heurdata->nintblocks > 0 )
   {
      assert(heurdata->intblockstart != NULL);
      assert(heurdata->intblockend != NULL);

      SCIPfreeBlockMemoryArray(scip, &heurdata->intblockstart, heurdata->nintblocks);
      SCIPfreeBlockMemoryArray(scip, &heurdata->intblockend, heurdata->nintblocks);
   }

   /* free the allocated memory for the integers */
   if( heurdata->intvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->intvars, heurdata->nintvars);
   }

   heurdata->nbinblocks = 0;
   heurdata->nintblocks = 0;
   heurdata->nbinvars = 0;
   heurdata->nintvars = 0;

   assert(heurdata->binvars == NULL);
   assert(heurdata->intvars == NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolTwoopt)
{
   SCIP_HEURDATA* heurdata;
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);
   assert(heurdata->binvars == NULL && heurdata->intvars == NULL);
   assert(heurdata->binblockstart == NULL && heurdata->binblockend == NULL);
   assert(heurdata->intblockstart == NULL && heurdata->intblockend == NULL);

   /* set heuristic data to initial values, but increase the total number of runs */
   heurdata->nbinvars = 0;
   heurdata->nintvars = 0;
   heurdata->lastsolindex = -1;
   heurdata->presolved = FALSE;

#ifdef SCIP_STATISTIC
   ++(heurdata->nruns);
#endif

   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolTwoopt)
{
   SCIP_HEURDATA* heurdata;
   int nbinvars;
   int nintvars;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   nbinvars = heurdata->nbinvars;
   nintvars = heurdata->nintvars;

   /* free the allocated memory for the binary variables */
   if( heurdata->binvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->binvars, nbinvars);
   }
   if( heurdata->binblockstart != NULL )
   {
      assert(heurdata->binblockend != NULL);

      SCIPfreeBlockMemoryArray(scip, &heurdata->binblockstart, heurdata->nbinblocks);
      SCIPfreeBlockMemoryArray(scip, &heurdata->binblockend, heurdata->nbinblocks);
   }
   heurdata->nbinvars = 0;
   heurdata->nbinblocks = 0;

   if( heurdata->intblockstart != NULL )
   {
      assert(heurdata->intblockend != NULL);

      SCIPfreeBlockMemoryArray(scip, &heurdata->intblockstart, heurdata->nintblocks);
      SCIPfreeBlockMemoryArray(scip, &heurdata->intblockend, heurdata->nintblocks);
   }
   heurdata->nintblocks = 0;

   /* free the allocated memory for the integers */
   if( heurdata->intvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->intvars, nintvars);
   }

   heurdata->nintvars = 0;

   assert(heurdata->binvars == NULL && heurdata->intvars == NULL);
   assert(heurdata->binblockstart == NULL && heurdata->binblockend == NULL);
   assert(heurdata->intblockstart == NULL && heurdata->intblockend == NULL);

   /* set heuristic data */
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTwoopt)
{  /*lint --e{715}*/
   SCIP_HEURDATA*  heurdata;
   SCIP_SOL* bestsol;
   SCIP_SOL* worksol;
   SCIP_ROW** lprows;
   SCIP_Real* activities;
   SCIP_COL** cols;
   int ncols;
   int nbinvars;
   int nintvars;
   int ndiscvars;
   int nlprows;
   int i;
   int ncolsforsorting;
   SCIP_Bool improvement;
   SCIP_Bool presolthiscall;
   SCIP_Bool varboundserr;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* we need an LP */
   if( SCIPgetNLPRows(scip) == 0 )
      return SCIP_OKAY;

   bestsol = SCIPgetBestSol(scip);

   /* ensure that heuristic has not already been processed on current incumbent */
   if( bestsol == NULL || heurdata->lastsolindex == SCIPsolGetIndex(bestsol) )
      return SCIP_OKAY;

   heurdata->lastsolindex = SCIPsolGetIndex(bestsol);

   /* we can only work on solutions valid in the transformed space */
   if( SCIPsolIsOriginal(bestsol) )
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintSol(scip, bestsol, NULL, TRUE) );
#endif

   /* ensure that the user defined number of nodes after last best solution has been reached, return otherwise */
   if( (SCIPgetNNodes(scip) - SCIPsolGetNodenum(bestsol)) < heurdata->waitingnodes )
      return SCIP_OKAY;

   presolthiscall = FALSE;
   SCIP_CALL( SCIPgetLPColsData(scip,&cols, &ncols) );
   ndiscvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   ncolsforsorting = MIN(ncols, ndiscvars);

   /* ensure that heuristic specific presolve is applied when heuristic is executed first */
   if( !heurdata->presolved )
   {
      SCIP_CALL( SCIPgetLPColsData(scip,&cols, &ncols) );

      for( i = 0; i < ncolsforsorting; ++i )
         SCIPcolSort(cols[i]);

      SCIP_CALL( presolveTwoOpt(scip, heur, heurdata) );
      presolthiscall = TRUE;
   }

   assert(heurdata->presolved);

   SCIPdebugMsg(scip, "  Twoopt heuristic is %sexecuting.\n", heurdata->execute ? "" : "not ");
   /* ensure that presolve has detected structures in the problem to which the 2-optimization can be applied.
    * That is if variables exist which share a common set of LP-rows. */
   if( !heurdata->execute )
      return SCIP_OKAY;

   nbinvars = heurdata->nbinvars;
   nintvars = heurdata->nintvars;
   ndiscvars = nbinvars + nintvars;

   /* we need to be able to start diving from current node in order to resolve the LP
    * with continuous or implicit integer variables
    */
   if( SCIPgetNVars(scip) > ndiscvars && ( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL ) )
      return SCIP_OKAY;

   /* problem satisfies all necessary conditions for 2-optimization heuristic, execute heuristic! */
   *result = SCIP_DIDNOTFIND;

   /* initialize a working solution as a copy of the current incumbent to be able to store
    * possible improvements obtained by 2-optimization */
   SCIP_CALL( SCIPcreateSolCopy(scip, &worksol, bestsol) );
   SCIPsolSetHeur(worksol, heur);

   /* get the LP row activities from current incumbent bestsol */
   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activities, nlprows) );

   for( i = 0; i < nlprows; i++ )
   {
      SCIP_ROW* row;

      row = lprows[i];
      assert(row != NULL);
      assert(SCIProwGetLPPos(row) == i);
      SCIPdebugMsg(scip, "  Row <%d> is %sin LP: \n", i, SCIProwGetLPPos(row) >= 0 ? "" : "not ");
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
      activities[i] = SCIPgetRowSolActivity(scip, row, bestsol);

      /* Heuristic does not provide infeasibility recovery, thus if any constraint is violated,
       * execution has to be terminated.
       */
      if( !SCIProwIsLocal(row) && (SCIPisFeasLT(scip, activities[i], SCIProwGetLhs(row))
            || SCIPisFeasGT(scip, activities[i], SCIProwGetRhs(row))) )
         goto TERMINATE;
   }

   if( !presolthiscall )
   {
      for( i = 0; i < ncolsforsorting; ++i )
         SCIPcolSort(cols[i]);
   }
   SCIPdebugMsg(scip, "  Twoopt heuristic has initialized activities and sorted rows! \n");

   /* start with binary optimization */
   improvement = FALSE;
   varboundserr = FALSE;

   if( heurdata->nbinblocks > 0 )
   {
      SCIP_CALL( optimize(scip, worksol, heurdata->binvars, heurdata->binblockstart, heurdata->binblockend, heurdata->nbinblocks,
            OPTTYPE_BINARY, activities, nlprows, &improvement, &varboundserr, heurdata) );

      SCIPdebugMsg(scip, "  Binary Optimization finished!\n");
   }

   if( varboundserr )
      goto TERMINATE;

   /* ensure that their are at least two integer variables which do not have the same coefficient
    * in the objective function. In one of these cases, the heuristic will automatically skip the
    * integer variable optimization */
   if( heurdata->nintblocks > 0 )
   {
      assert(heurdata->intopt);
      SCIP_CALL( optimize(scip, worksol, heurdata->intvars, heurdata->intblockstart, heurdata->intblockend, heurdata->nintblocks,
            OPTTYPE_INTEGER, activities, nlprows, &improvement, &varboundserr, heurdata) );

      SCIPdebugMsg(scip, "  Integer Optimization finished!\n");
   }

   if( ! improvement || varboundserr )
      goto TERMINATE;

   if( SCIPgetNVars(scip) == ndiscvars )
   {
      /* the problem is a pure IP, hence, no continuous or implicit variables are left for diving.
       * try if new working solution is feasible in original problem */
      SCIP_Bool success;
#ifndef NDEBUG
      SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
#else
      SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );
#endif

      if( success )
      {
         SCIPdebugMsg(scip, "found feasible shifted solution:\n");
         SCIPdebug( SCIP_CALL( SCIPprintSol(scip, worksol, NULL, FALSE) ) );
         heurdata->lastsolindex = SCIPsolGetIndex(bestsol);
         *result = SCIP_FOUNDSOL;

#ifdef SCIP_STATISTIC
        SCIPstatisticMessage("***Twoopt improved solution found by %10s . \n",
            SCIPsolGetHeur(bestsol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(bestsol)) :"Tree");
#endif
      }
   }
   /* fix the integer variables and start diving to optimize continuous variables with respect to reduced domain */
   else
   {
      SCIP_VAR** allvars;
      SCIP_Bool lperror;
#ifdef NDEBUG
      SCIP_RETCODE retstat;
#endif

      SCIPdebugMsg(scip, "shifted solution should be feasible -> solve LP to fix continuous variables to best values\n");

      allvars = SCIPgetVars(scip);

#ifdef SCIP_DEBUG
      for( i = ndiscvars; i < SCIPgetNVars(scip); ++i )
      {
         SCIPdebugMsg(scip, "  Cont. variable <%s>, status %d with bounds [%g <= %g <= x <= %g <= %g]\n",
            SCIPvarGetName(allvars[i]), SCIPvarGetStatus(allvars[i]), SCIPvarGetLbGlobal(allvars[i]), SCIPvarGetLbLocal(allvars[i]), SCIPvarGetUbLocal(allvars[i]),
            SCIPvarGetUbGlobal(allvars[i]));
      }
#endif

      /* start diving to calculate the LP relaxation */
      SCIP_CALL( SCIPstartDive(scip) );

      /* set the bounds of the variables: fixed for integers, global bounds for continuous */
      for( i = 0; i < SCIPgetNVars(scip); ++i )
      {
         if( SCIPvarGetStatus(allvars[i]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPchgVarLbDive(scip, allvars[i], SCIPvarGetLbGlobal(allvars[i])) );
            SCIP_CALL( SCIPchgVarUbDive(scip, allvars[i], SCIPvarGetUbGlobal(allvars[i])) );
         }
      }

      /* apply this after global bounds to not cause an error with intermediate empty domains */
      for( i = 0; i < ndiscvars; ++i )
      {
         if( SCIPvarGetStatus(allvars[i]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_Real solval;

            solval = SCIPgetSolVal(scip, worksol, allvars[i]);

            assert(SCIPvarGetType(allvars[i]) != SCIP_VARTYPE_CONTINUOUS && SCIPisFeasIntegral(scip, solval));

            SCIP_CALL( SCIPchgVarLbDive(scip, allvars[i], solval) );
            SCIP_CALL( SCIPchgVarUbDive(scip, allvars[i], solval) );
         }
      }
      for( i = 0; i < ndiscvars; ++i )
      {
         assert( SCIPisFeasEQ(scip, SCIPgetVarLbDive(scip, allvars[i]), SCIPgetVarUbDive(scip, allvars[i])) );
      }
      /* solve LP */
      SCIPdebugMsg(scip, " -> old LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));

      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop. */
#ifdef NDEBUG
      retstat = SCIPsolveDiveLP(scip, -1, &lperror, NULL);
      if( retstat != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving LP in Twoopt heuristic; LP solve terminated with code <%d>\n",retstat);
      }
#else
      SCIP_CALL( SCIPsolveDiveLP(scip, -1, &lperror, NULL) );
#endif

      SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
      SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, SCIPgetLPSolstat(scip));

      /* check if this is a feasible solution */
      if( !lperror && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIP_Bool success;

         /* copy the current LP solution to the working solution */
         SCIP_CALL( SCIPlinkLPSol(scip, worksol) );

         /* check solution for feasibility */
#ifndef NDEBUG
         SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
#else
         SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );
#endif

         if( success )
         {
            SCIPdebugMsg(scip, "found feasible shifted solution:\n");
            SCIPdebug( SCIP_CALL( SCIPprintSol(scip, worksol, NULL, FALSE) ) );
            heurdata->lastsolindex = SCIPsolGetIndex(bestsol);
            *result = SCIP_FOUNDSOL;

#ifdef SCIP_STATISTIC
            SCIPstatisticMessage("***   Twoopt improved solution found by %10s . \n",
               SCIPsolGetHeur(bestsol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(bestsol)) :"Tree");
#endif
         }
      }

      /* terminate the diving */
      SCIP_CALL( SCIPendDive(scip) );
   }

 TERMINATE:
   SCIPdebugMsg(scip, "Termination of Twoopt heuristic\n");
   SCIPfreeBufferArray(scip, &activities);
   SCIP_CALL( SCIPfreeSol(scip, &worksol) );

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the twoopt primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTwoopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Twoopt primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTwoopt, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTwoopt) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTwoopt) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitTwoopt) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitTwoopt) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolTwoopt) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolTwoopt) );

   /* include boolean flag intopt */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/twoopt/intopt", " Should Integer-2-Optimization be applied or not?",
         &heurdata->intopt, TRUE, DEFAULT_INTOPT, NULL, NULL) );

   /* include parameter waitingnodes */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/twoopt/waitingnodes", "user parameter to determine number of "
         "nodes to wait after last best solution before calling heuristic",
         &heurdata->waitingnodes, TRUE, DEFAULT_WAITINGNODES, 0, 10000, NULL, NULL));

   /* include parameter maxnslaves */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/twoopt/maxnslaves", "maximum number of slaves for one master variable",
         &heurdata->maxnslaves, TRUE, DEFAULT_MAXNSLAVES, -1, 1000000, NULL, NULL) );

   /* include parameter matchingrate */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/twoopt/matchingrate",
         "parameter to determine the percentage of rows two variables have to share before they are considered equal",
         &heurdata->matchingrate, TRUE, DEFAULT_MATCHINGRATE, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
