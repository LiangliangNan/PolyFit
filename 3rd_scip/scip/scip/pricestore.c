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

/**@file   pricestore.c
 * @brief  methods for storing priced variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/pricestore.h"
#include "scip/pub_message.h"

#include "scip/struct_pricestore.h"



/*
 * dynamic memory arrays
 */

/** resizes vars and score arrays to be able to store at least num entries */
static
SCIP_RETCODE pricestoreEnsureVarsMem(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(pricestore != NULL);
   assert(set != NULL);

   if( num > pricestore->varssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&pricestore->vars, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&pricestore->scores, newsize) );
      pricestore->varssize = newsize;
   }
   assert(num <= pricestore->varssize);

   return SCIP_OKAY;
}

/** resizes bdviolvars arrays to be able to store at least num entries */
static
SCIP_RETCODE pricestoreEnsureBdviolvarsMem(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(pricestore != NULL);
   assert(set != NULL);

   if( num > pricestore->bdviolvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&pricestore->bdviolvars, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&pricestore->bdviolvarslb, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&pricestore->bdviolvarsub, newsize) );
      pricestore->bdviolvarssize = newsize;
   }
   assert(num <= pricestore->bdviolvarssize);

   return SCIP_OKAY;
}


/** creates pricing storage */
SCIP_RETCODE SCIPpricestoreCreate(
   SCIP_PRICESTORE**     pricestore          /**< pointer to store pricing storage */
   )
{
   assert(pricestore != NULL);

   SCIP_ALLOC( BMSallocMemory(pricestore) );

   SCIP_CALL( SCIPclockCreate(&(*pricestore)->probpricingtime, SCIP_CLOCKTYPE_DEFAULT) );
   (*pricestore)->vars = NULL;
   (*pricestore)->scores = NULL;
   (*pricestore)->bdviolvars = NULL;
   (*pricestore)->bdviolvarslb = NULL;
   (*pricestore)->bdviolvarsub = NULL;
   (*pricestore)->varssize = 0;
   (*pricestore)->nvars = 0;
   (*pricestore)->bdviolvarssize = 0;
   (*pricestore)->nbdviolvars = 0;
   (*pricestore)->naddedbdviolvars = 0;
   (*pricestore)->nprobpricings = 0;
   (*pricestore)->nprobvarsfound = 0;
   (*pricestore)->nvarsfound = 0;
   (*pricestore)->nvarsapplied = 0;
   (*pricestore)->initiallp = FALSE;

   return SCIP_OKAY;
}

/** frees pricing storage */
SCIP_RETCODE SCIPpricestoreFree(
   SCIP_PRICESTORE**     pricestore          /**< pointer to store pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(*pricestore != NULL);
   assert((*pricestore)->nvars == 0);
   assert((*pricestore)->nbdviolvars == 0);

   SCIPclockFree(&(*pricestore)->probpricingtime);
   BMSfreeMemoryArrayNull(&(*pricestore)->vars);
   BMSfreeMemoryArrayNull(&(*pricestore)->scores);
   BMSfreeMemoryArrayNull(&(*pricestore)->bdviolvars);
   BMSfreeMemoryArrayNull(&(*pricestore)->bdviolvarslb);
   BMSfreeMemoryArrayNull(&(*pricestore)->bdviolvarsub);
   BMSfreeMemory(pricestore);

   return SCIP_OKAY;
}

/** informs pricing storage, that the setup of the initial LP starts now */
void SCIPpricestoreStartInitialLP(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(!pricestore->initiallp);
   assert(pricestore->nvars == 0);

   pricestore->initiallp = TRUE;
}

/** informs pricing storage, that the setup of the initial LP is now finished */
void SCIPpricestoreEndInitialLP(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->initiallp);
   assert(pricestore->nvars == 0);

   pricestore->initiallp = FALSE;
}

/** adds variable to pricing storage and capture it */
SCIP_RETCODE SCIPpricestoreAddVar(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_VAR*             var,                /**< priced variable */
   SCIP_Real             score,              /**< pricing score of variable (the larger, the better the variable) */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   int maxpricevars;
   int v;

   assert(pricestore != NULL);
   assert(set != NULL);
   assert(var != NULL);

#ifndef NDEBUG
   /* check if we add this variables to the same scip, where we created it */
   if( var->scip != set->scip )
   {
      SCIPerrorMessage("try to add a variable of another scip instance\n");
      return SCIP_INVALIDDATA;
   }
#endif

   SCIPsetDebugMsg(set, "adding variable <%s> (lb=%g, ub=%g) to pricing storage (initiallp=%u)\n",
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), pricestore->initiallp);

   if( pricestore->initiallp )
      maxpricevars = INT_MAX;
   else
   {
      pricestore->nvarsfound++;
      maxpricevars = SCIPsetGetPriceMaxvars(set, root);
   }
   assert(maxpricevars >= 1);
   assert(pricestore->nvars <= maxpricevars);

   /* check, if variable belongs to the best "maxpricevars" pricing variables */
   if( pricestore->nvars < maxpricevars || score > pricestore->scores[maxpricevars-1] )
   {
      /* capture variable */
      SCIPvarCapture(var);

      /* if the array consists of "maxpricevars" variables, release the worst variables */
      if( pricestore->nvars == maxpricevars )
      {
         SCIP_CALL( SCIPvarRelease(&pricestore->vars[pricestore->nvars-1], blkmem, set, eventqueue, lp) );
         pricestore->nvars--;
      }
      assert(pricestore->nvars < maxpricevars);

      /* get enough memory to store additional variable */
      SCIP_CALL( pricestoreEnsureVarsMem(pricestore, set, pricestore->nvars+1) );
      assert(pricestore->nvars <= pricestore->varssize);

      /* insert the variable at the correct position in sorted arrays */
      for( v = pricestore->nvars; v > 0 && score > pricestore->scores[v-1]; --v )
      {
         pricestore->vars[v] = pricestore->vars[v-1];
         pricestore->scores[v] = pricestore->scores[v-1];
      }
      pricestore->vars[v] = var;
      pricestore->scores[v] = score;
      pricestore->nvars++;
   }

   return SCIP_OKAY;
}

/** adds variable where zero violates the bounds to pricing storage, capture it */
SCIP_RETCODE SCIPpricestoreAddBdviolvar(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var                 /**< variable, where zero violates the bounds */
   )
{
   assert(pricestore != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPsetIsPositive(set, SCIPvarGetLbLocal(var)) || SCIPsetIsNegative(set, SCIPvarGetUbLocal(var)));
   assert(pricestore->naddedbdviolvars <= pricestore->nbdviolvars);

   SCIPsetDebugMsg(set, "zero violates bounds of <%s> (lb=%g, ub=%g)\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

   if( !pricestore->initiallp )
      pricestore->nvarsfound++;

   /* get enough memory to store additional variable */
   SCIP_CALL( pricestoreEnsureBdviolvarsMem(pricestore, set, pricestore->nbdviolvars+1) );
   assert(pricestore->nbdviolvars <= pricestore->bdviolvarssize);

   /* capture variable */
   SCIPvarCapture(var);

   /* insert variable in bdviolvars arrays */
   pricestore->bdviolvars[pricestore->nbdviolvars] = var;
   pricestore->bdviolvarslb[pricestore->nbdviolvars] = SCIPvarGetLbLocal(var);
   pricestore->bdviolvarsub[pricestore->nbdviolvars] = SCIPvarGetUbLocal(var);
   pricestore->nbdviolvars++;

   /* Temporarily set bounds, such that zero is feasible, because we don't want to destroy
    * dual feasibility (by adding columns) and primal feasibility (by introducing violated bounds)
    * at the same time.
    * The correct bounds must be reset with a call to SCIPpricestoreResetBounds().
    * The inference information is unimportant for this temporary bound change.
    */
   if( SCIPsetIsPositive(set, SCIPvarGetLbLocal(var)) )
   {
      SCIP_CALL( SCIPvarChgLbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, 0.0) );
   }
   else
   {
      SCIP_CALL( SCIPvarChgUbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, 0.0) );
   }

   return SCIP_OKAY;
}

/** adds given problem variable to pricing storage, if zero is not best bound w.r.t. objective function */
static
SCIP_RETCODE addBoundViolated(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool*            added               /**< pointer to store whether variable was added to pricing storage */
   )
{
   assert(tree != NULL);
   assert(added != NULL);

   *added = FALSE;

   /* add variable, if zero is not feasible within the bounds */
   if( SCIPsetIsPositive(set, SCIPvarGetLbLocal(var)) || SCIPsetIsNegative(set, SCIPvarGetUbLocal(var)) )
   {
      SCIPsetDebugMsg(set, " -> zero violates bounds of <%s> [%g,%g]\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      SCIP_CALL( SCIPpricestoreAddBdviolvar(pricestore, blkmem, set, stat, lp, branchcand, eventqueue, var) );
      *added = TRUE;
   }
   else
   {
      SCIP_Real bestbound;

      /* add variable, if zero is not best bound w.r.t. objective function */
      bestbound = SCIPvarGetBestBoundLocal(var);
      if( !SCIPsetIsZero(set, bestbound) )
      {
         SCIPsetDebugMsg(set, " -> best bound of <%s> [%g,%g] is not zero but %g\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), bestbound);
         SCIP_CALL( SCIPpricestoreAddVar(pricestore, blkmem, set, eventqueue, lp, var, 
               -SCIPvarGetObj(var) * bestbound, (SCIPtreeGetCurrentDepth(tree) == 0)) );
         *added = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** adds problem variables with negative reduced costs to pricing storage */
SCIP_RETCODE SCIPpricestoreAddProbVars(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_VAR* var;
   SCIP_COL* col;
   SCIP_Bool root;
   SCIP_Bool added;
   int v;
   int abortpricevars;
   int maxpricevars;
   int nfoundvars;

   assert(pricestore != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(tree != NULL);
   assert(SCIPtreeHasCurrentNodeLP(tree));
   assert(prob->nvars >= SCIPlpGetNCols(lp));

   /* if all problem variables of status COLUMN are already in the LP, nothing has to be done */
   if( prob->ncolvars == SCIPlpGetNCols(lp) )
      return SCIP_OKAY;

   root = (SCIPtreeGetCurrentDepth(tree) == 0);
   maxpricevars = SCIPsetGetPriceMaxvars(set, root);
   assert(maxpricevars >= 1);
   abortpricevars = (int)(set->price_abortfac * maxpricevars);
   assert(abortpricevars >= maxpricevars);

   /**@todo test pricing: is abortpricevars a good idea? -> like strong branching, lookahead, ... */

   pricestore->nprobpricings++;

   /* start timing */
   SCIPclockStart(pricestore->probpricingtime, set);

   /* price already existing problem variables */
   nfoundvars = 0;
   for( v = 0; v < prob->nvars && nfoundvars < abortpricevars; ++v )
   {
      var = prob->vars[v];
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      {
         col = SCIPvarGetCol(var);
         assert(col != NULL);
         assert(col->var == var);
         assert(col->len >= 0);
         assert(col->lppos >= -1);
         assert(col->lpipos >= -1);
         assert(SCIPcolIsInLP(col) == (col->lpipos >= 0));

         if( !SCIPcolIsInLP(col) )
         {
            SCIPsetDebugMsg(set, "price column variable <%s> in bounds [%g,%g]\n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

            /* add variable to pricing storage, if zero is not best bound w.r.t. objective function */
            SCIP_CALL( addBoundViolated(pricestore, blkmem, set, stat, tree, lp, branchcand, eventqueue, var, &added) );

            if( added )
            {
               pricestore->nprobvarsfound++;
               nfoundvars++;
            }
            else if( SCIPcolGetNNonz(col) > 0 )
            {
               SCIP_Real feasibility;

               /* a column not in LP that doesn't have zero in its bounds was added by bound checking above */
               assert(!SCIPsetIsPositive(set, SCIPvarGetLbLocal(col->var)));
               assert(!SCIPsetIsNegative(set, SCIPvarGetUbLocal(col->var)));

               if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE )
               {
                  /* The LP was proven infeasible, so we have an infeasibility proof by the dual Farkas multipliers y.
                   * The valid inequality  y^T A x >= y^T b  is violated by all x, especially by the (for this
                   * inequality most feasible solution) x' defined by 
                   *    x'_i = ub_i, if y^T A_i > 0
                   *    x'_i = lb_i, if y^T A_i <= 0.
                   * Pricing in this case means to add variables i with positive Farkas value, i.e. y^T A_i x'_i > 0
                   */
                  feasibility = -SCIPcolGetFarkasValue(col, stat, lp);
                  SCIPsetDebugMsg(set, "  <%s> Farkas feasibility: %e\n", SCIPvarGetName(col->var), feasibility);
               }
               else
               {
                  /* The dual LP is feasible, and we have a feasible dual solution. Pricing in this case means to
                   * add variables with negative feasibility, that is
                   *  - positive reduced costs for variables with negative lower bound
                   *  - negative reduced costs for variables with positive upper bound
                   */
                  feasibility = SCIPcolGetFeasibility(col, set, stat, lp);
                  SCIPsetDebugMsg(set, "  <%s> reduced cost feasibility: %e\n", SCIPvarGetName(col->var), feasibility);
               }

               /* the score is -feasibility / (#nonzeros in column + 1) to prefer short columns
                * we must add variables with negative feasibility, but in order to not get a too large lower bound
                * due to missing columns, we better also add variables, that have a very small feasibility
                */
               if( !SCIPsetIsPositive(set, feasibility) )
               {
                  SCIP_CALL( SCIPpricestoreAddVar(pricestore, blkmem, set, eventqueue, lp, var, -feasibility / (col->len+1), root) );
                  pricestore->nprobvarsfound++;
                  nfoundvars++;
               }
            }
         }
      }
   }

   /* stop timing */
   SCIPclockStop(pricestore->probpricingtime, set);

   return SCIP_OKAY;
}

/** adds priced variables to the LP */
SCIP_RETCODE SCIPpricestoreApplyVars(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< LP data */
   )
{
   SCIP_VAR* var;
   SCIP_COL* col;
   int v;

   assert(pricestore != NULL);
   assert(pricestore->naddedbdviolvars <= pricestore->nbdviolvars);
   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(tree != NULL);
   assert(SCIPtreeIsFocusNodeLPConstructed(tree));

   SCIPsetDebugMsg(set, "adding %d variables (%d bound violated and %d priced vars) to %d LP columns\n",
      SCIPpricestoreGetNVars(pricestore), pricestore->nbdviolvars - pricestore->naddedbdviolvars,
      pricestore->nvars, SCIPlpGetNCols(lp));

   /* add the variables with violated bounds to LP */
   for( v = pricestore->naddedbdviolvars; v < pricestore->nbdviolvars; ++v )
   {
      var = pricestore->bdviolvars[v];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);
      assert(var->nuses >= 2); /* at least used in pricing storage and in problem */

      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         /* transform loose variable into column variable */
         SCIP_CALL( SCIPvarColumn(var, blkmem, set, stat, prob, lp) );
      }
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

      col = SCIPvarGetCol(var);
      assert(col != NULL);
      assert(col->lppos == -1);
      SCIPsetDebugMsg(set, "adding bound violated variable <%s> (lb=%g, ub=%g)\n", SCIPvarGetName(var),
         pricestore->bdviolvarslb[v], pricestore->bdviolvarsub[v]);
      SCIP_CALL( SCIPlpAddCol(lp, set, col, SCIPtreeGetCurrentDepth(tree)) );

      if( !pricestore->initiallp )
         pricestore->nvarsapplied++;
   }
   pricestore->naddedbdviolvars = pricestore->nbdviolvars;

   /* add the selected pricing variables to LP */
   for( v = 0; v < pricestore->nvars; ++v )
   {
      var = pricestore->vars[v];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);
      assert(var->nuses >= 2); /* at least used in pricing storage and in problem */

      /* transform variable into column variable, if needed */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         SCIP_CALL( SCIPvarColumn(var, blkmem, set, stat, prob, lp) );
      }
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

      col = SCIPvarGetCol(var);
      assert(col != NULL);
      assert(col->lppos == -1);
      SCIPsetDebugMsg(set, "adding priced variable <%s> (score=%g)\n", SCIPvarGetName(var), pricestore->scores[v]);
      SCIP_CALL( SCIPlpAddCol(lp, set, col, SCIPtreeGetCurrentDepth(tree)) );

      /* release the variable */
      SCIP_CALL( SCIPvarRelease(&pricestore->vars[v], blkmem, set, eventqueue, lp) );

      if( !pricestore->initiallp )
         pricestore->nvarsapplied++;
   }

   /* clear the pricing storage */
   pricestore->nvars = 0;

   return SCIP_OKAY;
}

/** reset variables' bounds violated by zero to its original value */
SCIP_RETCODE SCIPpricestoreResetBounds(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_VAR* var;
   int v;

   assert(pricestore != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(pricestore->nvars == 0);
   assert(pricestore->naddedbdviolvars == pricestore->nbdviolvars);

   /* reset variables' bounds, release them, and clear the boundviolation storage;
    * the inference information is unimportant in these removals of temporary bound changes
    */
   for( v = 0; v < pricestore->nbdviolvars; ++v )
   {
      var = pricestore->bdviolvars[v];
      assert(var != NULL);

      SCIPsetDebugMsg(set, "resetting bounds of <%s> from [%g,%g] to [%g,%g]\n", var->name,
         SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), pricestore->bdviolvarslb[v], pricestore->bdviolvarsub[v]);
      SCIP_CALL( SCIPvarChgLbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, pricestore->bdviolvarslb[v]) );
      SCIP_CALL( SCIPvarChgUbLocal(var, blkmem, set, stat, lp, branchcand, eventqueue, pricestore->bdviolvarsub[v]) );
      SCIP_CALL( SCIPvarRelease(&pricestore->bdviolvars[v], blkmem, set, eventqueue, lp) );
   }
   pricestore->naddedbdviolvars = 0;
   pricestore->nbdviolvars = 0;

   return SCIP_OKAY;
}

/** gets number of variables in pricing storage */
int SCIPpricestoreGetNVars(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->nbdviolvars >= pricestore->naddedbdviolvars);

   return pricestore->nvars + pricestore->nbdviolvars - pricestore->naddedbdviolvars;
}

/** gets number of variables in pricing storage whose bounds must be reset */
int SCIPpricestoreGetNBoundResets(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->nbdviolvars >= pricestore->naddedbdviolvars);

   return pricestore->nbdviolvars - pricestore->naddedbdviolvars;
}

/** gets time needed to price existing problem variables */
SCIP_Real SCIPpricestoreGetProbPricingTime(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return SCIPclockGetTime(pricestore->probpricingtime);
}

/** gets total number of calls to problem variable pricing */
int SCIPpricestoreGetNProbPricings(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->nprobpricings;
}

/** gets total number of times, a problem variable was priced in */
int SCIPpricestoreGetNProbvarsFound(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->nprobvarsfound;
}

/** get total number of variables found so far in pricing */
int SCIPpricestoreGetNVarsFound(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->nvarsfound;
}

/** get total number of variables priced into the LP so far */
int SCIPpricestoreGetNVarsApplied(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->nvarsapplied;
}

