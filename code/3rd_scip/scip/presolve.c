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

/**@file   presolve.c
 * @brief  methods for presolving
 * @author Michael Winkler
 */

/* all general presolving methods (not working on any specific kind of data(, e.g. consdata) should be implemented here */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/mem.h"
#include "scip/presolve.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/solve.h"
#include "scip/struct_scip.h"

/*
 * Local methods
 */

/** collect variable bound information for a variable set reduction and global implication; only variable which have the
 *  vartype != SCIP_VARTYPE_BINARY have variable bounds
 */
static
void collectNonBinaryVBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< set variable */
   int                   varidx,             /**< for lower bound set variable index, for upper bound set variable index
                                              *   + number of variables
                                              */
   int                   pos,                /**< variables's position in bdchinfos */
   int                   nredvars,           /**< number of reduced variables so far */
   SCIP_Real*            bounds,             /**< array of bounds where one of them must be fullfilled */
   SCIP_Bool*            boundtypes,         /**< array of bound types */
   SCIP_Real*            newbounds,          /**< array of implied bounds(, size is two times number of variables, first
                                              *   half for implied lower bounds, second for implied upper bounds)
                                              */
   int*                  counts,             /**< array of number of implication on a bound (, size is two times number
                                              *   of variables, first half for implied lower bounds, second for implied
                                              *   upper bounds)
                                              */
   int*                  countnonzeros,      /**< array to store the indices of non-zero entries in the counts array */
   int*                  ncountnonzeros,     /**< pointer to store the number of non-zero entries in the counts array */
   int*                  issetvar,           /**< array containing for set variables the position in the current set, or
                                              *   0 if it is not a set variable or -1, if it is a redundant(i.e. implies
                                              *   another set variable) set variables(, size is two times number of
                                              *   variables, first half for implied lower bounds, second for implied
                                              *   upper bounds) */
   int                   nvars,              /**< number of problem variables */
   int*                  foundbin,           /**< pointer to store the lowest index of a binary implication variable
                                              *   when found
                                              */
   int*                  foundnonbin,        /**< pointer to store the lowest index of a non-binary implication variable
                                              *   when found
                                              */
   int*                  implidx,            /**< array to store the variable indices (for upper bound 'nvars' is added
                                              *   to the index) which are implied
                                              */
   int*                  nimplidx,           /**< pointer to store the number of implied variables */
   SCIP_Real*            lastbounds          /**< array to remember last implied bounds before taken the current
                                              *   variable into account, first nvars for lower bound, second nvars for
                                              *   upper bound
                                              *
                                              *   this array is used when a set variable got redundant, because it
                                              *   implies another set variable, and we need to correct the counts array
                                              */
   )
{
   SCIP_VAR** implvars;
   SCIP_Real* implcoefs;
   SCIP_Real* implconsts;
   int nimpls;
   int idx;
   int w;

   assert(scip != NULL);
   assert(var != NULL);
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
   assert(varidx >= 0);
   assert(pos >= 0);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(newbounds != NULL);
   assert(counts != NULL);
   assert(issetvar != NULL);
   assert(2 * nvars > varidx);
   assert(foundbin != NULL);
   assert(foundnonbin != NULL);
   assert(implidx != NULL);
   assert(nimplidx != NULL);
   assert(lastbounds != NULL);

   /* 1. case: lower bound in set */
   if( !boundtypes[pos] )
   {
      assert(counts[varidx] <= pos - nredvars + 1);

      /* update implication counter of set variable */
      if( counts[varidx] == pos - nredvars )
      {
         ++counts[varidx];

         if( counts[varidx] == 1 )
         {
            assert(*ncountnonzeros < 2*nvars);
            countnonzeros[*ncountnonzeros] = varidx;
            ++(*ncountnonzeros);
            newbounds[varidx] = bounds[pos];
            lastbounds[*nimplidx] = SCIP_INVALID;
         }
         else if( newbounds[varidx] > bounds[pos] )
         {
            lastbounds[*nimplidx] = newbounds[varidx];
            newbounds[varidx] = bounds[pos];
         }
         else
         {
            lastbounds[*nimplidx] = SCIP_INVALID;
         }

         *foundnonbin = MIN(*foundnonbin, varidx);

         implidx[*nimplidx] = varidx;
         ++(*nimplidx);
      }

      nimpls = SCIPvarGetNVubs(var);
      implvars = SCIPvarGetVubVars(var);
      implcoefs = SCIPvarGetVubCoefs(var);
      implconsts = SCIPvarGetVubConstants(var);

      for( w = nimpls - 1; w >= 0; --w )
      {
         assert(!SCIPisZero(scip, implcoefs[w]));
         idx = SCIPvarGetProbindex(implvars[w]);

         /* do not use inactive variables */
         /* @todo if implvars[x] is aggregated, we could transform the variable into the active representation */
         if( idx < 0 || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(implvars[w]), SCIPvarGetUbGlobal(implvars[w])) )
            continue;

         /* the upper bound of implvars[w] is bounding upper bound of var */
         if( implcoefs[w] < 0.0 )
         {
            /* update the counters and implied bounds */
            idx += nvars;

            if( counts[idx] == pos - nredvars )
            {
               if( SCIPvarIsBinary(implvars[w]) )
               {
                  /* do not look at fixed variables */
                  if( SCIPvarGetLbGlobal(implvars[w]) > 0.5 || SCIPvarGetUbGlobal(implvars[w]) < 0.5 )
                     continue;

                  /* (implvars[w] = 1 ===> var <= implcoefs[w] + implconsts[w] and if implcoefs[w] +
                   *  implconsts[w] < bounds[pos]) ===> (because var => bounds[v] ===> implvars[w] = 0)
                   */
                  if( SCIPisFeasLT(scip, implcoefs[w] + implconsts[w], bounds[pos]) )
                  {
                     /* set variable 'var' with bound implies other set variable 'implvars[w]' with corresponding bound
                      * so we can remove the set variable 'var'
                      */
                     if( issetvar[idx] > 0 )
                     {
                        SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g\n",
                           SCIPvarGetName(var), ">=", bounds[pos], SCIPvarGetName(implvars[w]), "<=", 0.0);

                        issetvar[varidx] = -1;
                        break;
                     }

                     ++counts[idx];
                     *foundbin = MIN(*foundbin, idx);

                     if( counts[idx] == 1 )
                     {
                        assert(*ncountnonzeros < 2*nvars);
                        countnonzeros[*ncountnonzeros] = idx;
                        ++(*ncountnonzeros);
                     }

                     implidx[*nimplidx] = idx;
                     ++(*nimplidx);
                  }
               }
               /* if (implvars[w] = ub(implvars[w]) ==> var <= implcoefs[w]*SCIPvarGetUbGlobal(implvars[w]) +
                * implconsts[w]) but (var >= bounds[pos] with bounds[pos] >
                * implcoefs[w]*SCIPvarGetUbGlobal(implvars[w]) + implconsts[w]) it follows (new_ub(var) <
                * ub(var))
                */
               else if( SCIPisFeasLT(scip, implcoefs[w] * SCIPvarGetUbGlobal(implvars[w]) + implconsts[w], bounds[pos]) )
               {
                  SCIP_Real newub;

                  newub = (bounds[pos] - implconsts[w]) / implcoefs[w];

                  /* set variable 'var' with bound implies other set variable 'implvars[w]' with corresponding bound so
                   * we can remove the set variable 'var'
                   */
                  if( issetvar[idx] > 0 && newub <= bounds[issetvar[idx] - 1] )
                  {
                     SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g (%g)\n",
                        SCIPvarGetName(var), ">=", bounds[pos], SCIPvarGetName(implvars[w]), "<=", newub, bounds[issetvar[idx] - 1] );

                     issetvar[varidx] = -1;
                     break;
                  }

                  ++counts[idx];

                  if( counts[idx] == 1 )
                  {
                     assert(*ncountnonzeros < 2*nvars);
                     countnonzeros[*ncountnonzeros] = idx;
                     ++(*ncountnonzeros);
                     newbounds[idx] = newub;
                     lastbounds[*nimplidx] = SCIP_INVALID;
                  }
                  else if( newbounds[idx] < newub )
                  {
                     lastbounds[*nimplidx] = newbounds[idx];
                     newbounds[idx] = newub;
                  }
                  else
                     lastbounds[*nimplidx] = SCIP_INVALID;

                  *foundnonbin = MIN(*foundnonbin, idx);

                  implidx[*nimplidx] = idx;
                  ++(*nimplidx);
               }
            }
         }
         else
         {
            assert(counts[varidx] <= pos - nredvars + 1);

            /* update the counters and implied bounds */
            if( counts[idx] == pos - nredvars )
            {
               if( SCIPvarIsBinary(implvars[w]) )
               {
                  /* do not look at fixed variables */
                  if( SCIPvarGetLbGlobal(implvars[w]) > 0.5 || SCIPvarGetUbGlobal(implvars[w]) < 0.5 )
                     continue;

                  /* (implvars[w] = 0 ===> var <= implconsts[w] and if implconsts[w] < bounds[pos]) ===> (because
                   *  var >= bounds[pos] ===> implvars[w] = 1)
                   */
                  if( SCIPisFeasLT(scip, implconsts[w], bounds[pos]) )
                  {
                     /* set variable 'var' with bound implies other set variable 'implvars[w]' with corresponding bound
                      * so we can remove the set variable 'var'
                      */
                     if( issetvar[idx] > 0 )
                     {
                        SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g\n",
                           SCIPvarGetName(var), ">=", bounds[pos], SCIPvarGetName(implvars[w]), ">=", 1.0);

                        issetvar[varidx] = -1;
                        break;
                     }

                     ++counts[idx];
                     *foundbin = MIN(*foundbin, idx);

                     if( counts[idx] == 1 )
                     {
                        assert(*ncountnonzeros < 2*nvars);
                        countnonzeros[*ncountnonzeros] = idx;
                        ++(*ncountnonzeros);
                     }

                     implidx[*nimplidx] = idx;
                     ++(*nimplidx);
                  }
               }
               /* if (implvars[w] = lb(implvars[w]) => var <= implcoefs[w]*SCIPvarGetLbGlobal(implvars[w]) +
                * implconsts[w]) but (var >= bounds[pos] with bounds[pos] >
                * implcoefs[w]*SCIPvarGetLbGlobal(implvars[w]) + implconsts[w]) it follows (new_lb(var) >
                * lb(var))
                */
               else if( SCIPisFeasLT(scip, implcoefs[w]*SCIPvarGetLbGlobal(implvars[w]) + implconsts[w], bounds[pos]) )
               {
                  SCIP_Real newlb;

                  newlb = (bounds[pos] - implconsts[w]) / implcoefs[w];

                  /* set variable 'var' with bound implies other set variable 'implvars[w]' with corresponding bound so
                   * we can remove the set variable 'var'
                   */
                  if( issetvar[idx] > 0 && newlb >= bounds[issetvar[idx] - 1] )
                  {
                     SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g (%g)\n",
                        SCIPvarGetName(var), ">=", bounds[pos], SCIPvarGetName(implvars[w]), ">=", newlb, bounds[issetvar[idx] - 1] );

                     issetvar[varidx] = -1;
                     break;
                  }

                  ++counts[idx];

                  if( counts[idx] == 1 )
                  {
                     assert(*ncountnonzeros < 2*nvars);
                     countnonzeros[*ncountnonzeros] = idx;
                     ++(*ncountnonzeros);
                     lastbounds[*nimplidx] = SCIP_INVALID;
                     newbounds[idx] = newlb;
                  }
                  else if( newbounds[idx] > newlb )
                  {
                     lastbounds[*nimplidx] = newbounds[idx];
                     newbounds[idx] = newlb;
                  }
                  else
                     lastbounds[*nimplidx] = SCIP_INVALID;


                  *foundnonbin = MIN(*foundnonbin, idx);

                  implidx[*nimplidx] = idx;
                  ++(*nimplidx);
               }
            }
         }
      }
   }
   /* 2.case: upper bound in set */
   else
   {
      assert(boundtypes[pos]);
      assert(counts[varidx] <= pos - nredvars + 1);

      /* update implication counter of set variable */
      if( counts[varidx] == pos - nredvars )
      {
         ++counts[varidx];

         if( counts[varidx] == 1 )
         {
            assert(*ncountnonzeros < 2*nvars);
            countnonzeros[*ncountnonzeros] = varidx;
            ++(*ncountnonzeros);
            newbounds[varidx] = bounds[pos];
            lastbounds[*nimplidx] = SCIP_INVALID;
         }
         else if( newbounds[varidx] < bounds[pos] )
         {
            lastbounds[*nimplidx] = newbounds[varidx];
            newbounds[varidx] = bounds[pos];
         }
         else
         {
            lastbounds[*nimplidx] = SCIP_INVALID;
         }

         *foundnonbin = MIN(*foundnonbin, varidx);

         implidx[*nimplidx] = varidx;
         ++(*nimplidx);
      }

      nimpls = SCIPvarGetNVlbs(var);
      implvars = SCIPvarGetVlbVars(var);
      implcoefs = SCIPvarGetVlbCoefs(var);
      implconsts = SCIPvarGetVlbConstants(var);

      for( w = nimpls - 1; w >= 0; --w )
      {
         assert(!SCIPisZero(scip, implcoefs[w]));
         idx = SCIPvarGetProbindex(implvars[w]);

         /* do not use inactive variables */
         /* @todo if implvars[x] is aggregated, we could transform the variable into the active representation */
         if( idx < 0 || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(implvars[w]), SCIPvarGetUbGlobal(implvars[w])) )
            continue;

         /* the lower bound of implvars[w] is bounding lower bound of var */
         if( implcoefs[w] < 0.0 )
         {
            assert(counts[idx] <= pos - nredvars + 1);

            /* update the counters and implied bounds */
            if( counts[idx] == pos - nredvars )
            {
               if( SCIPvarIsBinary(implvars[w]) )
               {
                  /* do not look at fixed variables */
                  if( SCIPvarGetLbGlobal(implvars[w]) > 0.5 || SCIPvarGetUbGlobal(implvars[w]) < 0.5 )
                     continue;

                  /* (implvars[w] = 0 ===> var >= implconsts[w] and if implconsts[w] > bounds[pos]) ===> (because
                   *  var <= bounds[pos] ===> implvars[w] = 1)
                   */
                  if( SCIPisFeasGT(scip, implconsts[w], bounds[pos]) )
                  {
                     if( issetvar[idx] > 0 )
                     {
                        SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g\n",
                           SCIPvarGetName(var), "<=", bounds[pos], SCIPvarGetName(implvars[w]), ">=", 1.0);

                        issetvar[varidx] = -1;
                        break;
                     }

                     ++counts[idx];
                     *foundbin = MIN(*foundbin, idx);

                     if( counts[idx] == 1 )
                     {
                        assert(*ncountnonzeros < 2*nvars);
                        countnonzeros[*ncountnonzeros] = idx;
                        ++(*ncountnonzeros);
                     }

                     implidx[*nimplidx] = idx;
                     ++(*nimplidx);
                  }
               }
               /* if (implvars[w] = lb(implvars[w]) ==> var <= implcoefs[w]*SCIPvarGetLbGlobal(implvars[w]) +
                * implconsts[w]) but (var <= bounds[pos] with bounds[pos] <
                * implcoefs[w]*SCIPvarGetLbGlobal(implvars[w]) + implconsts[w]) it follows (new_lb(var) >
                * ub(var))
                */
               else if( SCIPisFeasGT(scip, implcoefs[w]*SCIPvarGetLbGlobal(implvars[w]) + implconsts[w], bounds[pos]) )
               {
                  SCIP_Real newlb;

                  newlb = (bounds[pos] - implconsts[w]) / implcoefs[w];

                  if( issetvar[idx] > 0 && newlb >= bounds[issetvar[idx] - 1] )
                  {
                     SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g (%g)\n",
                        SCIPvarGetName(var), "<=", bounds[pos], SCIPvarGetName(implvars[w]), ">=",  newlb, bounds[issetvar[idx] - 1]);

                     issetvar[varidx] = -1;
                     break;
                  }

                  ++counts[idx];

                  if( counts[idx] == 1 )
                  {
                     assert(*ncountnonzeros < 2*nvars);
                     countnonzeros[*ncountnonzeros] = idx;
                     ++(*ncountnonzeros);
                     lastbounds[*nimplidx] = SCIP_INVALID;
                     newbounds[idx] = newlb;
                  }
                  else if( newbounds[idx] > newlb )
                  {
                     lastbounds[*nimplidx] = newbounds[idx];
                     newbounds[idx] = newlb;
                  }
                  else
                     lastbounds[*nimplidx] = SCIP_INVALID;

                  *foundnonbin = MIN(*foundnonbin, idx);

                  implidx[*nimplidx] = idx;
                  ++(*nimplidx);
               }
            }
         }
         else
         {
            /* update the counters and implied bounds */
            idx += nvars;

            assert(counts[idx] <= pos - nredvars + 1);

            if( counts[idx] == pos - nredvars )
            {
               if( SCIPvarIsBinary(implvars[w]) )
               {
                  /* do not look at fixed variables */
                  if( SCIPvarGetLbGlobal(implvars[w]) > 0.5 || SCIPvarGetUbGlobal(implvars[w]) < 0.5 )
                     continue;

                  /* (implvars[w] = 1 ===> var >= implcoefs[w] + implconsts[w] and if implcoefs[w] +
                   *  implconsts[w] > bounds[pos]) ===> (because var <= bounds[pos] ===> implvars[w] = 0)
                   */
                  if( SCIPisFeasGT(scip, implcoefs[w] + implconsts[w], bounds[pos]) )
                  {
                     if( issetvar[idx] > 0 )
                     {
                        SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g\n",
                           SCIPvarGetName(var), "<=", bounds[pos], SCIPvarGetName(implvars[w]), "<=", 0.0);

                        issetvar[varidx] = -1;
                        break;
                     }

                     ++counts[idx];
                     *foundbin = MIN(*foundbin, idx);

                     if( counts[idx] == 1 )
                     {
                        assert(*ncountnonzeros < 2*nvars);
                        countnonzeros[*ncountnonzeros] = idx;
                        ++(*ncountnonzeros);
                     }

                     implidx[*nimplidx] = idx;
                     ++(*nimplidx);
                  }
               }
               /* if (implvars[w] = ub(implvars[w]) => var <= implcoefs[w]*SCIPvarGetUbGlobal(implvars[w]) +
                * implconsts[w]) but (var <= bounds[pos] with bounds[pos] <
                * implcoefs[w]*SCIPvarGetUbGlobal(implvars[w]) + implconsts[w]) it follows (new_ub(var) <
                * ub(var))
                */
               else if( SCIPisFeasGT(scip, implcoefs[w]*SCIPvarGetUbGlobal(implvars[w]) + implconsts[w], bounds[pos]) )
               {
                  SCIP_Real newub;

                  newub = (bounds[pos] - implconsts[w]) / implcoefs[w];

                  if( issetvar[idx] > 0 && newub <= bounds[issetvar[idx] - 1] )
                  {
                     SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g (%g)\n",
                        SCIPvarGetName(var), "<=", bounds[pos], SCIPvarGetName(implvars[w]), "<=",  newub, bounds[issetvar[idx] - 1]);

                     issetvar[varidx] = -1;
                     break;
                  }

                  ++counts[idx];

                  if( counts[idx] == 1 )
                  {
                     assert(*ncountnonzeros < 2*nvars);
                     countnonzeros[*ncountnonzeros] = idx;
                     ++(*ncountnonzeros);
                     lastbounds[*nimplidx] = SCIP_INVALID;
                     newbounds[idx] = newub;
                  }
                  else if( newbounds[idx] < newub )
                  {
                     lastbounds[*nimplidx] = newbounds[idx];
                     newbounds[idx] = newub;
                  }
                  else
                     lastbounds[*nimplidx] = SCIP_INVALID;

                  *foundnonbin = MIN(*foundnonbin, idx);

                  implidx[*nimplidx] = idx;
                  ++(*nimplidx);
               }
            }
         }
      }
   }
}


/** collect non-binary implication data for variable set reduction and global bound implications; only variable which
 *  have the vartype SCIP_VARTYPE_BINARY have implications, otherwise the implications are saved as variable bounds
 */
static
void collectNonBinaryImplicationData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< set variable */
   int                   varidx,             /**< for lower bound set variable index, for upper bound set
                                              *   variable index + number of variables
                                              */
   int                   pos,                /**< variables's position in bdchinfos */
   int                   nredvars,           /**< number of reduced variables so far */
   SCIP_Bool             value,              /**< value used for clique and implication info */
   SCIP_Real*            bounds,             /**< array of bounds where one of them must be fullfilled */
   SCIP_Bool*            boundtypes,         /**< array of bound types */
   SCIP_Real*            newbounds,          /**< array of implied bounds(, size is two times number of variables, first
                                              *   half for implied lower bounds, second for implied upper bounds)
                                              */
   int*                  counts,             /**< array of number of implication on a bound (, size is two times number
                                              *   of variables, first half for implied lower bounds, second for implied
                                              *   upper bounds)
                                              */
   int*                  countnonzeros,      /**< array to store the indices of non-zero entries in the counts array */
   int*                  ncountnonzeros,     /**< pointer to store the number of non-zero entries in the counts array */
   int*                  issetvar,           /**< array containing for set variables the position in the current set, or
                                              *   0 if it is not a set variable or -1, if it is a redundant(i.e. implies
                                              *   another set variable) set variables(, size is two times number of
                                              *   variables, first half for implied lower bounds, second for implied
                                              *   upper bounds) */
   int                   nvars,              /**< number of problem variables */
   int*                  foundbin,           /**< pointer to store the lowest index of a binary implication variable
                                              *   when found
                                              */
   int*                  foundnonbin,        /**< pointer to store the lowest index of a non-binary implication variable
                                              *   when found
                                              */
   int*                  implidx,            /**< array to store the variable indices (for upper bound 'nvars' is added
                                              *   to the index) which are implied
                                              */
   int*                  nimplidx,           /**< pointer to store the number of implied variables */
   SCIP_Real*            lastbounds          /**< array to remember last implied bounds before taken the current
                                              *   variable into account, first nvars for lower bound, second nvars for
                                              *   upper bound
                                              *
                                              *   this array is used when a set variable got redundant, because it
                                              *   implies another set variable, and we need to correct the counts array
                                              */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(varidx >= 0);
   assert(pos >= 0);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(newbounds != NULL);
   assert(counts != NULL);
   assert(issetvar != NULL);
   assert(2 * nvars > varidx);
   assert(foundbin != NULL);
   assert(foundnonbin != NULL);
   assert(implidx != NULL);
   assert(nimplidx != NULL);
   assert(lastbounds != NULL);

   if( issetvar[varidx] > 0 )
   {
      SCIP_VAR** implvars;
      SCIP_Real* implbounds;
      SCIP_BOUNDTYPE* implboundtypes;
      int idx;
      int w;

      /* get implication information */
      implvars = SCIPvarGetImplVars(var, value);
      implboundtypes = SCIPvarGetImplTypes(var, value);
      implbounds = SCIPvarGetImplBounds(var, value);

      for( w = SCIPvarGetNImpls(var, value) - 1; w >= 0; --w )
      {
         assert(implvars != NULL);
         assert(implboundtypes != NULL);

         /* no self implication should exist in the implication data structure */
         assert(implvars[w] != var);

         idx = SCIPvarGetProbindex(implvars[w]);

         /* do not use inactive variables */
         /* @todo if implvars[x] is aggregated, we could transform the variable into the active representation */
         if( idx < 0 || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(implvars[w]), SCIPvarGetUbGlobal(implvars[w])) )
            continue;

         if( implboundtypes[w] == SCIP_BOUNDTYPE_UPPER )
         {
            idx += nvars;

            assert(counts[idx] <= pos - nredvars + 1);

            /* set variable 'var' with bound implies other set variable 'implvars[w]' with a non-worse bound than the
             * bound so we can remove the set variable 'var'
             */
            if( issetvar[idx] > 0 && bounds[issetvar[idx] - 1] >= implbounds[w] )
            {
               SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g (%g)\n",
                  SCIPvarGetName(var), boundtypes[pos] ? "<=" : ">=", bounds[pos], SCIPvarGetName(implvars[w]),
                  "<=", implbounds[w], bounds[issetvar[idx] - 1]);

               issetvar[varidx] = -1;
               break;
            }

            /* update implication counter and implied upper bound */
            if( counts[idx] == pos - nredvars )
            {
               ++counts[idx];

               if( SCIPvarIsBinary(implvars[w]) )
               {
                  /* the implied upper bound on a binary variable should not be trivial, otherwise we might globally fix
                   * this variable to a wrong value
                   *
                   * @note it is possible that the implied bound is lower than zero, when the implied variable has
                   * become binary during the search
                   */
                  assert(SCIPisFeasLE(scip, implbounds[w], 0.0));
                  *foundbin = MIN(*foundbin, idx);

                  if( counts[idx] == 1 )
                  {
                     assert(*ncountnonzeros < 2*nvars);
                     countnonzeros[*ncountnonzeros] = idx;
                     ++(*ncountnonzeros);
                  }
               }
               else
               {
                  *foundnonbin = MIN(*foundnonbin, idx);

                  if( counts[idx] == 1 )
                  {
                     assert(*ncountnonzeros < 2*nvars);
                     countnonzeros[*ncountnonzeros] = idx;
                     ++(*ncountnonzeros);
                     newbounds[idx] = implbounds[w];
                     lastbounds[*nimplidx] = SCIP_INVALID;
                  }
                  else if( newbounds[idx] < implbounds[w] )
                  {
                     lastbounds[*nimplidx] = newbounds[idx];
                     newbounds[idx] = implbounds[w];
                  }
                  else
                     lastbounds[*nimplidx] = SCIP_INVALID;
               }

               implidx[*nimplidx] = idx;
               ++(*nimplidx);
            }
         }
         else
         {
            assert(counts[idx] <= pos - nredvars + 1);

            /* set variable 'var' with bound implies other set variable 'implvars[w]' with a non-worse bound than the
             * bound so we can remove the set variable 'var'
             */
            if( issetvar[idx] > 0 && bounds[issetvar[idx] - 1] <= implbounds[w] )
            {
               SCIPdebugMsg(scip, "set variable <%s> %s %g implies other set variable <%s> %s %g (%g)\n",
                  SCIPvarGetName(var), boundtypes[pos] ? "<=" : ">=", bounds[pos], SCIPvarGetName(implvars[w]),
                  ">=", implbounds[w], bounds[issetvar[idx] - 1]);

               issetvar[varidx] = -1;
               break;
            }

            /* update implication counter */
            if( counts[idx] == pos - nredvars )
            {
               ++counts[idx];

               if( SCIPvarIsBinary(implvars[w]) )
               {
                  /* the implied lower bound on a binary variable should not be trivial, otherwise we might globally fix
                   * this variable to a wrong value
                   *
                   * @note is is possible that the implied bound is greater than one, when the implied variable has
                   * become binary during the search
                   */
                  assert(SCIPisFeasGE(scip, implbounds[w], 1.0));
                  *foundbin = MIN(*foundbin, idx);

                  if( counts[idx] == 1 )
                  {
                     assert(*ncountnonzeros < 2*nvars);
                     countnonzeros[*ncountnonzeros] = idx;
                     ++(*ncountnonzeros);
                  }
               }
               else
               {
                  *foundnonbin = MIN(*foundnonbin, idx);

                  if( counts[idx] == 1 )
                  {
                     assert(*ncountnonzeros < 2*nvars);
                     countnonzeros[*ncountnonzeros] = idx;
                     ++(*ncountnonzeros);
                     newbounds[idx] = implbounds[w];
                     lastbounds[*nimplidx] = SCIP_INVALID;
                  }
                  else if( newbounds[idx] > implbounds[w] )
                  {
                     lastbounds[*nimplidx] = newbounds[idx];
                     newbounds[idx] = implbounds[w];
                  }
                  else
                     lastbounds[*nimplidx] = SCIP_INVALID;
               }

               implidx[*nimplidx] = idx;
               ++(*nimplidx);
            }
         }
      }
   }
}

/** collect clique data on binary variables for variable set reduction and global bound implications */
static
void collectBinaryCliqueData(
   SCIP_VAR*             var,                /**< set variable */
   int                   varidx,             /**< for lower bound set variable index, for upper bound set variable index
                                              *   + number of variables
                                              */
   int                   pos,                /**< variables's position in bdchinfos */
   int                   nredvars,           /**< number of reduced variables so far */
   SCIP_Bool             value,              /**< value used for clique and implication info */
   SCIP_Real*            bounds,             /**< array of bounds where one of them must be fullfilled */
   SCIP_Bool*            boundtypes,         /**< array of bound types */
   SCIP_Real*            newbounds,          /**< array of implied bounds(, size is two times number of variables, first
                                              *   half for implied lower bounds, second for implied upper bounds)
                                              */
   int*                  counts,             /**< array of number of implication on a bound (, size is two times number of
                                              *   variables, first half for implied lower bounds, second for implied upper
                                              *   bounds)
                                              */
   int*                  countnonzeros,      /**< array to store the indices of non-zero entries in the counts array */
   int*                  ncountnonzeros,     /**< pointer to store the number of non-zero entries in the counts array */
   int*                  issetvar,           /**< array containing for set variables the position in the current set, or
                                              *   0 if it is not a set variable, or -1, if it is a redundant (i.e. implies
                                              *   another set variable) set variable
                                              *   (the size of the array is two times the number of variables, first half
                                              *   for implied lower bounds, second for implied upper bounds)
                                              */
   int                   nvars,              /**< number of problem variables */
   int*                  foundbin,           /**< pointer to store the lowest index of a binary implication variable when found */
   int*                  implidx,            /**< array to store the variable indices (for upper bound 'nvars' is added
                                              *   to the index) which are implied
                                              */
   int*                  nimplidx            /**< pointer to store the number of implied variables */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_VAR** clqvars;
   SCIP_Bool* clqvalues;
   int idx;
   int c;
   int w;

   assert(var != NULL);
   assert(SCIPvarIsBinary(var));
   assert(varidx >= 0);
   assert(pos >= 0);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(newbounds != NULL);
   assert(counts != NULL);
   assert(issetvar != NULL);
   assert(2 * nvars > varidx);
   assert(foundbin != NULL);
   assert(implidx != NULL);
   assert(nimplidx != NULL);

   /* implication counter cannot exceed number implication variables */
   assert(counts[varidx] <= pos - nredvars);

   /* if the set variable is not yet redundant we might increase the self implication counter */
   if( issetvar[varidx] > 0 )
   {
      /* update implication counter for set variables */
      if( counts[varidx] == pos - nredvars )
      {
         ++counts[varidx];
         *foundbin = MIN(*foundbin, varidx);

         if( counts[varidx] == 1 )
         {
            assert(*ncountnonzeros < 2*nvars);
            countnonzeros[*ncountnonzeros] = varidx;
            ++(*ncountnonzeros);
         }

         implidx[*nimplidx] = varidx;
         ++(*nimplidx);
      }
   }

   cliques = SCIPvarGetCliques(var, value);

   /* update implication counter on all by cliques implied variables */
   for( c = SCIPvarGetNCliques(var, value) - 1; c >= 0; --c )
   {
      clqvars = SCIPcliqueGetVars(cliques[c]);
      clqvalues = SCIPcliqueGetValues(cliques[c]);

      for( w = SCIPcliqueGetNVars(cliques[c]) - 1; w >= 0; --w )
      {
         /* already handle self-implication and do not look at fixed variables */
         if( clqvars[w] == var || SCIPvarGetLbGlobal(clqvars[w]) > 0.5 || SCIPvarGetUbGlobal(clqvars[w]) < 0.5 )
            continue;

         idx = SCIPvarGetProbindex(clqvars[w]);
         assert(idx >= 0);

         if( clqvalues[w] )
            idx += nvars;

         assert(counts[idx] <= pos - nredvars + 1);

         /* set variable 'var' with bound implies other set variable 'clqvars[w]' with corresponding set bound so we can
          * remove the set variable 'var'
          */
         if( issetvar[idx] > 0 )
         {
            SCIPdebugMessage("set variable <%s> %s %g implies other set variable <%s> %s %g\n",
               SCIPvarGetName(var), boundtypes[pos] ? "<=" : ">=", bounds[pos], SCIPvarGetName(clqvars[w]),
               clqvalues[w] ? "<=" : ">=", clqvalues[w] ? 0.0 : 1.0);

            issetvar[varidx] = -1;
            break;
         }

         /* update implication counter */
         if( counts[idx] == pos - nredvars )
         {
            ++counts[idx];
            *foundbin = MIN(*foundbin, idx);

            if( counts[idx] == 1 )
            {
               assert(*ncountnonzeros < 2*nvars);
               countnonzeros[*ncountnonzeros] = idx;
               ++(*ncountnonzeros);
            }

            implidx[*nimplidx] = idx;
            ++(*nimplidx);
         }
      }
   }
}



/*
 * presolving methods
 */

#define CLEARRATIO 0.8

/** try to reduce the necessary variable in a set of variables with corresponding bounds and boundtypes for which one
 *  must be fulfilled
 *
 *  e.g. a set of logicor or bounddisjunctive constraint variables would be such a set
 *
 *  consider the following set:
 *
 *  x1 >= 1, x2 >= 3, x3 >= 1, x4 <= 0
 *
 *  by (global) implication data (cliques, implications, and variable bounds) we have also the following implications
 *  given:
 *
 *  x1 >= 1 => x3 >= 1
 *  x2 >= 2 => x3 >= 1
 *  x4 <= 0 => x1 >= 1
 *
 *  Because of the last implication x4 is redundant, because x1 >= 1 would also be fulfilled in the variable set, so we
 *  can reduce the set by x4.
 *  Also, the both other implications and x3 >= 1 (in the given variable set) all imply exactly x3 >= 1, so we tighten
 *  the global lower bound of x3 to 1 and the set of variables gets redundant.
 */
SCIP_RETCODE SCIPshrinkDisjunctiveVarSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables array for which at least one must be fulfilled in the
                                              *   following bounds and boundtypes */
   SCIP_Real*            bounds,             /**< bounds array for which at least one must be fulfilled */
   SCIP_Bool*            boundtypes,         /**< boundtypes array (TRUE == SCIP_BOUNDTYPE_UPPER, FALSE == SCIP_BOUNDTYPE_LOWER)
                                              *   for which at least one must be fulfilled */
   SCIP_Bool*            redundants,         /**< array which be filled and then indicate if a variable in the set is redundant */
   int                   nvars,              /**< number of variables */
   int*                  nredvars,           /**< pointer to store how many variables can be removed */
   int*                  nglobalred,         /**< pointer to store number of global reductions on variable bounds found
                                              *   through this set of variables */
   SCIP_Bool*            setredundant,       /**< pointer to store if we found a global reduction on a variable which was part
                                              *   of the given set of variables, this makes this disjunction redundant */
   SCIP_Bool*            glbinfeas,          /**< pointer to store if global infeasibility was detected */
   SCIP_Bool             fullshortening      /**< do we want to try the shortening procedure over the whole set (which might be expensive) */
   )
{
   SCIP_Real* newbounds; /* array saving all overall implied global bounds, first nprobvars for lower bound, second
                          * nprobvars for upper bound
                          */
   SCIP_Real* lastbounds;/* temporary array to remember last implied bounds before taken the current variable into
                          * account, first nprobvars for lower bound, second nprobvars for upper bound
                          *
                          * this array is used when a set variable got redundant, because it implies another set
                          * variable, and we need to correct the counts array
                          */
   int* issetvar;        /* array for mapping from a problem variable to the position in the variable set (,pos + 1 in
                          * the variable set, 0 for no set variable, and -1 if this variable was removed from the set),
                          * first nprobvars for lower bound, second nprobvars for upper bound
                          */
   int* counts;          /* array saving number of implication by set variables, first nprobvars for lower bound, second
                          * nprobvars for upper bound
                          */
   int* implidx;         /* temporary array to remember all indices of implied variables by the current set variable
                          * looked at, first nprobvars for lower bound, second nprobvars for upper bound
                          *
                          * this array is used when a set variable got redundant, because it implies another set
                          * variable, and we need to correct the counts array
                          */
   int* countnonzeros;

   SCIP_VAR* var;
   SCIP_Bool usebin = TRUE;
   SCIP_Bool usenonbin = TRUE;
   SCIP_Bool globalred = TRUE;
   SCIP_Bool reducedset;
   SCIP_Bool value;
   SCIP_Bool implbinvarsexist;
   int start = INT_MAX;
   int nimplidx;
   int foundbin;
   int foundnonbin;
   int varidx;
   int nprobvars;
   int ncountnonzeros;
   int maxcountnonzeros;
   int w;
   int v;

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(redundants != NULL);
   assert(nredvars != NULL);
   assert(nglobalred != NULL);
   assert(setredundant != NULL);
   assert(glbinfeas != NULL);
   assert(scip->transprob != NULL);
   nprobvars = SCIPprobGetNVars(scip->transprob);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &issetvar, 2*nprobvars) ); /*lint !e647*/
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &counts, 2*nprobvars) ); /*lint !e647*/
   SCIP_CALL( SCIPallocBufferArray(scip, &newbounds, 2*nprobvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lastbounds, 2*nprobvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &implidx, 2*nprobvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &countnonzeros, 2*nprobvars) );

   *nredvars = 0;
   *glbinfeas = FALSE;
   ncountnonzeros = 0;

   maxcountnonzeros = (int)(2*nprobvars*CLEARRATIO); /*lint !e790*/

   /* initialize variable indices data */
   for( v = 0; v < nvars; ++v )
   {
      varidx = SCIPvarGetProbindex(vars[v]);
      assert(varidx >= 0);

      if( boundtypes[v] )
         varidx += nprobvars;

      /* initialize issetvar array */
      issetvar[varidx] = v+1;
   }

   /* check if implied binary variables exist, because for these variables the implications can be stored in the
    * variable bounds instead of the 'normal' implications
    */
   implbinvarsexist = (SCIPprobGetNImplBinVars(scip->transprob) > 0);

#if 0
   /* @todo do the cleanup here rather than before calling SCIPshrinkDisjunctiveVarSet()? */
   if( usebin )
   {
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPcleanupCliques(scip, &infeasible) );
   }
#endif

   /* check for same implied binary variables */
   for( v = 0; v < nvars; ++v )
   {
      var =  vars[v];

      foundbin = INT_MAX;
      foundnonbin = INT_MAX;
      reducedset = FALSE;
      nimplidx = 0;

      value = (!boundtypes[v]);

      varidx = SCIPvarGetProbindex(var);
      assert(varidx >= 0);

      if( !value )
         varidx += nprobvars;

      if( usebin )
      {
         /* collect clique data on binary variables */
         if( SCIPvarIsBinary(var) )
         {
            collectBinaryCliqueData(var, varidx, v, *nredvars, value, bounds, boundtypes, newbounds, counts, countnonzeros,
               &ncountnonzeros, issetvar, nprobvars, &foundbin, implidx, &nimplidx);
         }
      }

      /* only variable which have the vartype SCIP_VARTYPE_BINARY have implications, otherwise the implications are
       * saved as variable bounds
       *
       * we only check binary to non-binary implications if we can detect further implications which either lead to
       * global reductions or to redundant set variables
       */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY && ((usebin && implbinvarsexist) || usenonbin) )
      {
         collectNonBinaryImplicationData(scip, var, varidx, v, *nredvars, value, bounds, boundtypes, newbounds, counts,
            countnonzeros, &ncountnonzeros, issetvar, nprobvars, &foundbin, &foundnonbin, implidx, &nimplidx, lastbounds);
      }
      /* only variable which have the vartype != SCIP_VARTYPE_BINARY have variable bounds
       *
       * we only check the variable bounds if we can detect further implications which either lead to global reductions
       * or to redundant set variables
       */
      else if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY && ((usebin && implbinvarsexist) || usenonbin) )
      {
         collectNonBinaryVBoundData(scip, var, varidx, v, *nredvars, bounds, boundtypes, newbounds, counts, countnonzeros,
            &ncountnonzeros, issetvar, nprobvars, &foundbin, &foundnonbin, implidx, &nimplidx, lastbounds);
      }

      /* reduce implication counters on all variables which are implied by a variable now marked as redundant */
      if( issetvar[varidx] < 0 )
      {
         SCIP_VAR** probvars;

         SCIPdebugMsg(scip, "marked variable <%s> as redundant variable in variable set\n", SCIPvarGetName(var));

         probvars = SCIPprobGetVars(scip->transprob);
         assert(probvars != NULL);

         /* correct implication counters and bounds, if the redundant variable implies other variables we need to reduce
          * the counter and get the last bounds before this implication
          */
         for( w = nimplidx - 1; w >= 0; --w )
         {
            assert(implidx[w] < 2 * nprobvars);
            assert(counts[implidx[w]] == v - (*nredvars) + 1);

            --counts[implidx[w]];

            if( implidx[w] == countnonzeros[ncountnonzeros-1] && counts[implidx[w]] == 0 )
               --ncountnonzeros;

            /* only for non-binary variables we need to correct the bounds */
            if( implidx[w] < nprobvars )
            {
               if( !SCIPvarIsBinary(probvars[implidx[w]]) && lastbounds[w] != SCIP_INVALID )/*lint !e777*/
                  newbounds[implidx[w]] = lastbounds[w];
            }
            else
            {
               if( !SCIPvarIsBinary(probvars[implidx[w] - nprobvars]) && lastbounds[w] != SCIP_INVALID )/*lint !e777*/
                  newbounds[implidx[w] - nprobvars] = lastbounds[w];
            }
         }

         reducedset = TRUE;
         ++(*nredvars);
      }

      /* check if we want to shorten the whole set of variables, or terminate early if we did not find any further
       * implication
       */
      if( !fullshortening )
      {
         /* check if it makes sense to look for further binary implications */
         if( foundbin < INT_MAX && !reducedset )
            usebin = FALSE;
         /* check if it makes sense to look for further non-binary implications */
         if( foundnonbin < INT_MAX && !reducedset )
            usenonbin = FALSE;
      }

      /* are global reductions still possible */
      globalred = globalred && (foundbin < INT_MAX || foundnonbin < INT_MAX);

      /* remember the first possible position for a global bound change */
      if( globalred )
      {
         /* get correct variable index(, we added nprobvars for the upper bound implication) */
         if( foundbin < INT_MAX && foundbin >= nprobvars )
            foundbin -= nprobvars;

         /* get correct variable index(, we added nprobvars for the upper bound implication) */
         if( foundnonbin < INT_MAX && foundnonbin >= nprobvars )
            foundnonbin -= nprobvars;

         if( start > foundbin )
            start = foundbin;

         if( start > foundnonbin )
            start = foundnonbin;
      }
      else
         start = INT_MAX;

      /* check if it can find any global implications anymore */
      if( !usebin && !usenonbin )
         break;
   }

   /* remove redundant set variables */
   if( *nredvars > 0 )
   {
#ifndef NDEBUG
      int nreds = 0;
#endif

      for( v = nvars - 1; v >= 0; --v )
      {
         var = vars[v];

         varidx = SCIPvarGetProbindex(var);
         assert(varidx >= 0);

         if( boundtypes[v] )
            varidx += nprobvars;

         /* if set variable was marked to be redundant remove it */
         if( issetvar[varidx] < 0 )
         {
            SCIPdebugMsg(scip, "mark redundant variable <%s> to be removed from variable set\n", SCIPvarGetName(var));

            redundants[v] = TRUE;
#ifndef NDEBUG
            ++nreds;
#endif
         }
      }
      assert((*nredvars) == nreds);
   }

   /* if we found some global boundchanges, we perform then now */
   if( globalred )
   {
      SCIP_VAR** probvars;
      SCIP_VAR* probvar;

      SCIPdebugMsg(scip, "variable set led to global reductions (in %s)\n", SCIPprobGetName(scip->transprob));

      probvars = SCIPprobGetVars(scip->transprob);
      assert(probvars != NULL);

      assert(start < nprobvars);

      /* check for same implied binary variables */
      for( v = start; v < nprobvars; ++v )
      {
         probvar = probvars[v];
         assert(probvar != NULL);

         assert(counts[v] <= nvars);
         assert(counts[nprobvars + v] <= nvars);

         if( counts[v] + (*nredvars) == nvars )
         {
            if( SCIPvarIsBinary(probvar) )
            {
               SCIPdebugMsg(scip, "can fix variable %s [%g, %g] to 1.0\n", SCIPvarGetName(probvar),
                  SCIPvarGetLbGlobal(probvar), SCIPvarGetUbGlobal(probvar));

               if( SCIPvarGetLbGlobal(probvar) < 0.5 )
               {
                  SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
                        scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand,
                        scip->eventqueue, scip->cliquetable, probvar, 1.0, SCIP_BOUNDTYPE_LOWER, FALSE) );

                  assert(SCIPvarGetLbGlobal(probvar) > 0.5 || scip->tree->npendingbdchgs > 0);

                  ++(*nglobalred);

                  if( issetvar[v] > 0 )
                     *setredundant = TRUE;
               }
            }
            else
            {
               SCIPdebugMsg(scip, "can tighten lower bound variable %s [%g, %g] to %g\n", SCIPvarGetName(probvar),
                  SCIPvarGetLbGlobal(probvar), SCIPvarGetUbGlobal(probvar), newbounds[v]);

               /* the new lower bound is greater than the global upper bound => the problem is global infeasible */
               if( SCIPisLT(scip, SCIPvarGetUbGlobal(probvar), newbounds[v]) )
               {
                  SCIPdebugMsg(scip, "-> global infeasibility proven.\n");

                  SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
                  *glbinfeas = TRUE;
                  break;
               }

               if( SCIPisLT(scip, SCIPvarGetLbGlobal(probvar), newbounds[v]) )
               {
                  SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
                        scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand,
                        scip->eventqueue, scip->cliquetable, probvar, newbounds[v], SCIP_BOUNDTYPE_LOWER, FALSE) );

                  ++(*nglobalred);

                  if( issetvar[v] > 0 && newbounds[v] >= bounds[issetvar[v] - 1] )
                     *setredundant = TRUE;
               }
            }
         }
         else if( counts[nprobvars + v] + (*nredvars) == nvars )
         {
            if( SCIPvarIsBinary(probvar) )
            {
               SCIPdebugMsg(scip, "can fix variable %s [%g, %g] to 0.0\n", SCIPvarGetName(probvar),
                  SCIPvarGetLbGlobal(probvar), SCIPvarGetUbGlobal(probvar));

               if( SCIPvarGetUbGlobal(probvar) > 0.5 )
               {
                  SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
                        scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
                        scip->cliquetable, probvar, 0.0, SCIP_BOUNDTYPE_UPPER, FALSE) );

                  assert(SCIPvarGetUbGlobal(probvar) < 0.5 || scip->tree->npendingbdchgs > 0);

                  ++(*nglobalred);

                  if( issetvar[nprobvars + v] > 0 )
                     *setredundant = TRUE;
               }
            }
            else
            {
               int idx = nprobvars + v;

               SCIPdebugMsg(scip, "can tighten upper bound variable %s [%g, %g] to %g\n", SCIPvarGetName(probvar),
                  SCIPvarGetLbGlobal(probvar), SCIPvarGetUbGlobal(probvar), newbounds[idx]);

               /* the new upper bound is small than the global upper bound => the problem is global infeasible */
               if( SCIPisGT(scip, SCIPvarGetLbGlobal(probvar), newbounds[idx]) )
               {
                  SCIPdebugMsg(scip, "-> global infeasibility proven.\n");

                  SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
                  *glbinfeas = TRUE;
                  break;
               }

               if( SCIPisGT(scip, SCIPvarGetUbGlobal(probvar), newbounds[idx]) )
               {
                  SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
                        scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
                        scip->cliquetable, probvar, newbounds[idx], SCIP_BOUNDTYPE_UPPER, FALSE) );

                  ++(*nglobalred);

                  if( issetvar[idx] > 0 && newbounds[idx] <= bounds[issetvar[idx] - 1] )
                     *setredundant = TRUE;
               }
            }
         }
      }
   }

   /* reset issetvar array to 0 */
   for( v = 0; v < nvars; ++v )
   {
      varidx = SCIPvarGetProbindex(vars[v]);
      assert(varidx >= 0);

      if( boundtypes[v] )
         varidx += nprobvars;

      issetvar[varidx] = 0;
   }

   if( ncountnonzeros >= maxcountnonzeros )
   {
      BMSclearMemoryArray(counts, 2*nprobvars);
   }
   else
   {
      while( --ncountnonzeros >= 0 )
         counts[countnonzeros[ncountnonzeros]] = 0;
   }


   SCIPfreeBufferArray(scip, &countnonzeros);
   SCIPfreeBufferArray(scip, &implidx);
   SCIPfreeBufferArray(scip, &lastbounds);
   SCIPfreeBufferArray(scip, &newbounds);
   SCIPfreeCleanBufferArray(scip, &counts);
   SCIPfreeCleanBufferArray(scip, &issetvar);

   return SCIP_OKAY;
}
