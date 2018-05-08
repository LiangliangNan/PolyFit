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

/**@file   sol.c
 * @brief  methods for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/misc.h"
#include "scip/lp.h"
#include "scip/nlp.h"
#include "scip/relax.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/cons.h"
#include "scip/pub_message.h"

#ifndef NDEBUG
#include "scip/struct_sol.h"
#endif



/** clears solution arrays of primal CIP solution */
static
SCIP_RETCODE solClearArrays(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   SCIP_CALL( SCIPboolarrayClear(sol->valid) );
   sol->hasinfval = FALSE;

   return SCIP_OKAY;
}

/** sets value of variable in the solution's array */
static
SCIP_RETCODE solSetArrayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val                 /**< value to set variable to */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* from now on, variable must not be deleted */
   SCIPvarMarkNotDeletable(var);

   /* mark the variable valid */
   SCIP_CALL( SCIPboolarraySetVal(sol->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

   /* set the value in the solution array */
   SCIP_CALL( SCIPrealarraySetVal(sol->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, val) );

   /* store whether the solution has infinite values assigned to variables */
   if( val != SCIP_UNKNOWN ) /*lint !e777*/
      sol->hasinfval = (sol->hasinfval || SCIPsetIsInfinity(set, val) || SCIPsetIsInfinity(set, -val));

   return SCIP_OKAY;
}

/** increases value of variable in the solution's array */
static
SCIP_RETCODE solIncArrayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             incval              /**< increase of variable's solution value */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* from now on, variable must not be deleted */
   SCIPvarMarkNotDeletable(var);

   /* if the variable was not valid, mark it to be valid and set the value to the incval (it is 0.0 if not valid) */
   if( !SCIPboolarrayGetVal(sol->valid, idx) )
   {
      /* mark the variable valid */
      SCIP_CALL( SCIPboolarraySetVal(sol->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

      /* set the value in the solution array */
      SCIP_CALL( SCIPrealarraySetVal(sol->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, incval) );
   }
   else
   {
      /* increase the value in the solution array */
      SCIP_CALL( SCIPrealarrayIncVal(sol->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, incval) );
   }

   /* store whether the solution has infinite values assigned to variables */
   incval = SCIPrealarrayGetVal(sol->vals, idx);
   if( incval != SCIP_UNKNOWN ) /*lint !e777*/
      sol->hasinfval = (sol->hasinfval || SCIPsetIsInfinity(set, incval) || SCIPsetIsInfinity(set, -incval));

   return SCIP_OKAY;
}

/** returns the value of the variable in the given solution */
static
SCIP_Real solGetArrayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* check, if the variable's value is valid */
   if( SCIPboolarrayGetVal(sol->valid, idx) )
   {
      return SCIPrealarrayGetVal(sol->vals, idx);
   }
   else
   {
      /* return the variable's value corresponding to the origin */
      switch( sol->solorigin )
      {
      case SCIP_SOLORIGIN_ORIGINAL:
      case SCIP_SOLORIGIN_ZERO:
         return 0.0;

      case SCIP_SOLORIGIN_LPSOL:
         return SCIPvarGetLPSol(var);

      case SCIP_SOLORIGIN_NLPSOL:
         return SCIPvarGetNLPSol(var);

      case SCIP_SOLORIGIN_RELAXSOL:
         return SCIPvarGetRelaxSolTransVar(var);

      case SCIP_SOLORIGIN_PSEUDOSOL:
         return SCIPvarGetPseudoSol(var);

      case SCIP_SOLORIGIN_PARTIAL:
      case SCIP_SOLORIGIN_UNKNOWN:
         return SCIP_UNKNOWN;

      default:
         SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
         SCIPABORT();
         return 0.0; /*lint !e527*/
      }
   }
}

/** stores solution value of variable in solution's own array */
static
SCIP_RETCODE solUnlinkVar(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real solval;

   assert(sol != NULL);
   assert(var != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   /* if variable is already valid, nothing has to be done */
   if( SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var)) )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "unlinking solution value of variable <%s>\n", SCIPvarGetName(var));

   /* store the correct solution value into the solution array */
   switch( sol->solorigin )
   {
   case SCIP_SOLORIGIN_ORIGINAL:
      SCIPerrorMessage("cannot unlink solutions of original problem space\n");
      return SCIP_INVALIDDATA;

   case SCIP_SOLORIGIN_ZERO:
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_LPSOL:
      solval = SCIPvarGetLPSol(var);
      SCIP_CALL( solSetArrayVal(sol, set, var, solval) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_NLPSOL:
      solval = SCIPvarGetNLPSol(var);
      SCIP_CALL( solSetArrayVal(sol, set, var, solval) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_RELAXSOL:
      solval = SCIPvarGetRelaxSolTransVar(var);
      SCIP_CALL( solSetArrayVal(sol, set, var, solval) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PSEUDOSOL:
      solval = SCIPvarGetPseudoSol(var);
      SCIP_CALL( solSetArrayVal(sol, set, var, solval) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PARTIAL:
   case SCIP_SOLORIGIN_UNKNOWN:
      SCIP_CALL( solSetArrayVal(sol, set, var, SCIP_UNKNOWN) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
      return SCIP_INVALIDDATA;
   }
}

/** sets the solution time, nodenum, runnum, and depth stamp to the current values */
static
void solStamp(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_Bool             checktime           /**< should the time be updated? */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);

   if( checktime )
   {
      sol->time = SCIPclockGetTime(stat->solvingtime);
#ifndef NDEBUG
      sol->lpcount = stat->lpcount;
#endif
   }
   else
      sol->time = SCIPclockGetLastTime(stat->solvingtime);
   sol->nodenum = stat->nnodes;
   sol->runnum = stat->nruns;
   if( tree == NULL )
      sol->depth = -1;
   else
      sol->depth = SCIPtreeGetCurrentDepth(tree);
}

/** creates primal CIP solution, initialized to zero */
SCIP_RETCODE SCIPsolCreate(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_ZERO;
   (*sol)->obj = 0.0;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   SCIPsolResetViolations(*sol);
   stat->solindex++;
   solStamp(*sol, stat, tree, TRUE);
   SCIPsolResetViolations(*sol);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** creates primal CIP solution in original problem space, initialized to the offset in the original problem */
SCIP_RETCODE SCIPsolCreateOriginal(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_ORIGINAL;
   (*sol)->obj = origprob->objoffset;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   stat->solindex++;
   solStamp(*sol, stat, tree, TRUE);
   SCIPsolResetViolations(*sol);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** creates a copy of a primal CIP solution */
SCIP_RETCODE SCIPsolCopy(
   SCIP_SOL**            sol,                /**< pointer to store the copy of the primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   )
{
   assert(sol != NULL);
   assert(sourcesol != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCopy(&(*sol)->vals, blkmem, sourcesol->vals) );
   SCIP_CALL( SCIPboolarrayCopy(&(*sol)->valid, blkmem, sourcesol->valid) );
   (*sol)->heur = sourcesol->heur;
   (*sol)->obj = sourcesol->obj;
   (*sol)->primalindex = -1;
   (*sol)->time = sourcesol->time;
#ifndef NDEBUG
   (*sol)->lpcount = sourcesol->lpcount;
#endif
   (*sol)->nodenum = sourcesol->nodenum;
   (*sol)->solorigin = sourcesol->solorigin;
   (*sol)->runnum = sourcesol->runnum;
   (*sol)->depth = sourcesol->depth;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = sourcesol->hasinfval;
   stat->solindex++;
   (*sol)->viol.absviolbounds = sourcesol->viol.absviolbounds;
   (*sol)->viol.absviolcons = sourcesol->viol.absviolcons;
   (*sol)->viol.absviolintegrality = sourcesol->viol.absviolintegrality;
   (*sol)->viol.absviollprows = sourcesol->viol.absviollprows;
   (*sol)->viol.relviolbounds = sourcesol->viol.relviolbounds;
   (*sol)->viol.relviolcons = sourcesol->viol.relviolcons;
   (*sol)->viol.relviollprows = sourcesol->viol.relviollprows;

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** transformes given original solution to the transformed space; a corresponding transformed solution has to be given
 *  which is copied into the existing solution and freed afterwards
 */
SCIP_RETCODE SCIPsolTransform(
   SCIP_SOL*             sol,                /**< primal CIP solution to change, living in original space */
   SCIP_SOL**            transsol,           /**< pointer to corresponding transformed primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal              /**< primal data */
   )
{  /*lint --e{715}*/
   SCIP_REALARRAY* tmpvals;
   SCIP_BOOLARRAY* tmpvalid;
   SCIP_SOL* tsol;

   assert(sol != NULL);
   assert(transsol != NULL);
   assert(SCIPsolIsOriginal(sol));
   assert(sol->primalindex > -1);

   tsol = *transsol;
   assert(tsol != NULL);
   assert(!SCIPsolIsOriginal(tsol));

   /* switch vals and valid arrays; the exisiting solution gets the arrays of the transformed solution;
    * the transformed one gets the original arrays, because they have to be freed anyway and freeing the transsol
    * automatically frees its arrays
    */
   tmpvals = sol->vals;
   tmpvalid = sol->valid;
   sol->vals = tsol->vals;
   sol->valid = tsol->valid;
   tsol->vals = tmpvals;
   tsol->valid = tmpvalid;

   /* copy solorigin and objective (should be the same, only to avoid numerical issues);
    * we keep the other statistics of the original solution, since that was the first time that this solution as found
    */
   sol->solorigin = tsol->solorigin;
   sol->obj = tsol->obj;

   SCIP_CALL( SCIPsolFree(transsol, blkmem, primal) );

   return SCIP_OKAY;
}

/** adjusts solution values of implicit integer variables in handed solution. Solution objective value is not
 *  deteriorated by this method.
 */
SCIP_RETCODE SCIPsolAdjustImplicitSolVals(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< either original or transformed problem, depending on sol origin */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             uselprows           /**< should LP row information be considered for none-objective variables */
   )
{
   SCIP_VAR** vars;
   int nimplvars;
   int nbinvars;
   int nintvars;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);

   /* get variable data */
   vars = SCIPprobGetVars(prob);
   nbinvars = SCIPprobGetNBinVars(prob);
   nintvars = SCIPprobGetNIntVars(prob);
   nimplvars = SCIPprobGetNImplVars(prob);

   if( nimplvars == 0 )
      return SCIP_OKAY;

   /* calculate the last array position of implicit integer variables */
   nimplvars = nbinvars + nintvars + nimplvars;

   /* loop over implicit integer variables and round them up or down */
   for( v = nbinvars + nintvars; v < nimplvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real solval;
      SCIP_Real obj;
      SCIP_Real newsolval;
      SCIP_Bool roundup;
      SCIP_Bool rounddown;
      int nuplocks;
      int ndownlocks;

      var = vars[v];

      assert( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT );
      solval = SCIPsolGetVal(sol, set, stat, var);

      /* we do not need to round integral solution values or those of variables which are not column variables */
      if( SCIPsetIsFeasIntegral(set, solval) || SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
         continue;

      nuplocks = SCIPvarGetNLocksUp(var);
      ndownlocks = SCIPvarGetNLocksDown(var);
      obj = SCIPvarGetUnchangedObj(var);

      roundup = FALSE;
      rounddown = FALSE;

      /* in case of a non-zero objective coefficient, there is only one possible rounding direction */
      if( SCIPsetIsFeasNegative(set, obj) )
         roundup = TRUE;
      else if( SCIPsetIsFeasPositive(set, obj) )
         rounddown = TRUE;
      else if( uselprows )
      {
         /* determine rounding direction based on row violations */
         SCIP_COL* col;
         SCIP_ROW** rows;
         SCIP_Real* vals;
         int nrows;
         int r;

         col = SCIPvarGetCol(var);
         vals = SCIPcolGetVals(col);
         rows = SCIPcolGetRows(col);
         nrows = SCIPcolGetNNonz(col);

         /* loop over rows and search for equations whose violation can be decreased by rounding */
         for( r = 0; r < nrows && !(roundup && rounddown); ++r )
         {
            SCIP_ROW* row;
            SCIP_Real activity;
            SCIP_Real rhs;
            SCIP_Real lhs;

            row = rows[r];

            if( SCIProwIsLocal(row) || !SCIProwIsInLP(row) )
               continue;

            rhs = SCIProwGetRhs(row);
            lhs = SCIProwGetLhs(row);

            if( SCIPsetIsInfinity(set, rhs) || SCIPsetIsInfinity(set, -lhs) )
               continue;

            activity = SCIProwGetSolActivity(row, set, stat, sol);
            if( SCIPsetIsFeasLE(set, activity, rhs) && SCIPsetIsFeasLE(set, lhs, activity) )
               continue;

            assert(! SCIPsetIsZero(set, vals[r]));
            if( (SCIPsetIsFeasGT(set, activity, rhs) && SCIPsetIsPositive(set, vals[r]))
                  || (SCIPsetIsFeasLT(set, activity, lhs) && SCIPsetIsNegative(set, vals[r])) )
               rounddown = TRUE;
            else
               roundup = TRUE;
         }
      }

      /* in case of a tie, we select the rounding step based on the number of variable locks */
      if( roundup == rounddown )
      {
         rounddown = ndownlocks <= nuplocks;
         roundup = !rounddown;
      }

      /* round the variable up or down */
      if( roundup )
      {
         newsolval = SCIPsetCeil(set, solval);
         assert(SCIPsetIsFeasLE(set, newsolval, SCIPvarGetUbGlobal(var)));
      }
      else
      {
         assert( rounddown ); /* should be true because of the code above */
         newsolval = SCIPsetFloor(set, solval);
         assert(SCIPsetIsFeasGE(set, newsolval, SCIPvarGetLbGlobal(var)));
      }

      SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, var, newsolval) );
   }

   return SCIP_OKAY;
}
/** creates primal CIP solution, initialized to the current LP solution */
SCIP_RETCODE SCIPsolCreateLPSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(SCIPlpIsSolved(lp));

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkLPSol(*sol, set, stat, prob, tree, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current NLP solution */
SCIP_RETCODE SCIPsolCreateNLPSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(nlp != NULL);

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkNLPSol(*sol, stat, tree, nlp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current relaxation solution */
SCIP_RETCODE SCIPsolCreateRelaxSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(relaxation != NULL);
   assert(SCIPrelaxationIsSolValid(relaxation));

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkRelaxSol(*sol, set, stat, tree, relaxation) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current pseudo solution */
SCIP_RETCODE SCIPsolCreatePseudoSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkPseudoSol(*sol, set, stat, prob, tree, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current solution */
SCIP_RETCODE SCIPsolCreateCurrentSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(tree != NULL);

   if( SCIPtreeHasCurrentNodeLP(tree) )
   {
      SCIP_CALL( SCIPsolCreateLPSol(sol, blkmem, set, stat, prob, primal, tree, lp, heur) );
   }
   else
   {
      SCIP_CALL( SCIPsolCreatePseudoSol(sol, blkmem, set, stat, prob, primal, tree, lp, heur) );
   }

   return SCIP_OKAY;
}

/** creates partial primal CIP solution, initialized to unknown values */
SCIP_RETCODE SCIPsolCreatePartial(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(primal != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_PARTIAL;
   (*sol)->obj = SCIP_UNKNOWN;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   stat->solindex++;
   solStamp(*sol, stat, NULL, TRUE);
   SCIPsolResetViolations(*sol);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to unknown values */
SCIP_RETCODE SCIPsolCreateUnknown(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_UNKNOWN;
   (*sol)->obj = 0.0;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   stat->solindex++;
   solStamp(*sol, stat, tree, TRUE);
   SCIPsolResetViolations(*sol);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** frees primal CIP solution */
SCIP_RETCODE SCIPsolFree(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PRIMAL*          primal              /**< primal data */
   )
{
   assert(sol != NULL);
   assert(*sol != NULL);

   SCIPprimalSolFreed(primal, *sol);

   SCIP_CALL( SCIPrealarrayFree(&(*sol)->vals) );
   SCIP_CALL( SCIPboolarrayFree(&(*sol)->valid) );
   BMSfreeBlockMemory(blkmem, sol);

   return SCIP_OKAY;
}

/** copies current LP solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkLPSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(SCIPlpDiving(lp) || SCIPtreeProbing(tree) || !SCIPlpDivingObjChanged(lp));

   SCIPsetDebugMsg(set, "linking solution to LP\n");

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* link solution to LP solution */
   if( SCIPlpDivingObjChanged(lp) )
   {
      /* the objective value has to be calculated manually, because the LP's value is invalid;
       * use objective values of variables, because columns objective values are changed to dive values
       */
      sol->obj = SCIPlpGetLooseObjval(lp, set, prob);
      if( !SCIPsetIsInfinity(set, -sol->obj) )
      {
         SCIP_VAR* var;
         SCIP_COL** cols;
         int ncols;
         int c;

         cols = SCIPlpGetCols(lp);
         ncols = SCIPlpGetNCols(lp);
         for( c = 0; c < ncols; ++c )
         {
            var = SCIPcolGetVar(cols[c]);
            sol->obj += SCIPvarGetUnchangedObj(var) * cols[c]->primsol;
         }
      }
   }
   else
   {
      /* the objective value in the columns is correct, s.t. the LP's objective value is also correct */
      sol->obj = SCIPlpGetObjval(lp, set, prob);
   }
   sol->solorigin = SCIP_SOLORIGIN_LPSOL;
   solStamp(sol, stat, tree, TRUE);

   SCIPsetDebugMsg(set, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current NLP solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkNLPSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(nlp != NULL);
   assert(SCIPnlpGetSolstat(nlp) <= SCIP_NLPSOLSTAT_FEASIBLE);

   SCIPstatDebugMsg(stat, "linking solution to NLP\n");

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* get objective value of NLP solution */
   if( SCIPnlpIsDivingObjChanged(nlp) )
   {
      /* the objective value has to be calculated manually, because the NLP's value is invalid */

      SCIP_VAR** vars;
      int nvars;
      int v;

      sol->obj = 0.0;

      vars = SCIPnlpGetVars(nlp);
      nvars = SCIPnlpGetNVars(nlp);
      for( v = 0; v < nvars; ++v )
      {
         assert(SCIPvarIsActive(vars[v]));
         sol->obj += SCIPvarGetUnchangedObj(vars[v]) * SCIPvarGetNLPSol(vars[v]);
      }
   }
   else
   {
      sol->obj = SCIPnlpGetObjval(nlp);
   }

   sol->solorigin = SCIP_SOLORIGIN_NLPSOL;
   solStamp(sol, stat, tree, TRUE);

   SCIPstatDebugMsg(stat, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current relaxation solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkRelaxSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   )
{  /*lint --e{715}*/
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(relaxation != NULL);
   assert(SCIPrelaxationIsSolValid(relaxation));

   SCIPsetDebugMsg(set, "linking solution to relaxation\n");

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* the objective value in the columns is correct, s.t. the LP's objective value is also correct */
   sol->obj = SCIPrelaxationGetSolObj(relaxation);
   sol->solorigin = SCIP_SOLORIGIN_RELAXSOL;
   solStamp(sol, stat, tree, TRUE);

   SCIPsetDebugMsg(set, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current pseudo solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkPseudoSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   SCIPsetDebugMsg(set, "linking solution to pseudo solution\n");

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* link solution to pseudo solution */
   sol->obj = SCIPlpGetPseudoObjval(lp, set, prob);
   sol->solorigin = SCIP_SOLORIGIN_PSEUDOSOL;
   solStamp(sol, stat, tree, TRUE);

   SCIPsetDebugMsg(set, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current solution (LP or pseudo solution) into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkCurrentSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);

   SCIPsetDebugMsg(set, "linking solution to current solution\n");

   if( SCIPtreeHasCurrentNodeLP(tree) && SCIPlpIsSolved(lp) )
   {
      SCIP_CALL( SCIPsolLinkLPSol(sol, set, stat, prob, tree, lp) );
   }
   else
   {
      SCIP_CALL( SCIPsolLinkPseudoSol(sol, set, stat, prob, tree, lp) );
   }

   return SCIP_OKAY;
}

/** clears primal CIP solution */
SCIP_RETCODE SCIPsolClear(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(sol != NULL);

   SCIP_CALL( solClearArrays(sol) );
   sol->solorigin = SCIP_SOLORIGIN_ZERO;
   sol->obj = 0.0;
   solStamp(sol, stat, tree, TRUE);

   return SCIP_OKAY;
}

/** declares all entries in the primal CIP solution to be unknown */
SCIP_RETCODE SCIPsolSetUnknown(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(sol != NULL);

   SCIP_CALL( solClearArrays(sol) );
   sol->solorigin = SCIP_SOLORIGIN_UNKNOWN;
   sol->obj = 0.0;
   solStamp(sol, stat, tree, TRUE);

   return SCIP_OKAY;
}

/** stores solution values of variables in solution's own array */
SCIP_RETCODE SCIPsolUnlink(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< transformed problem data */
   )
{
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(prob->nvars == 0 || prob->vars != NULL);

   if( !SCIPsolIsOriginal(sol) && sol->solorigin != SCIP_SOLORIGIN_ZERO
      && sol->solorigin != SCIP_SOLORIGIN_UNKNOWN )
   {
      SCIPsetDebugMsg(set, "completing solution %p\n", (void*)sol);

      for( v = 0; v < prob->nvars; ++v )
      {
         SCIP_CALL( solUnlinkVar(sol, set, prob->vars[v]) );
      }

      sol->solorigin = SCIP_SOLORIGIN_ZERO;
   }

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolSetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Real             val                 /**< solution value of variable */
   )
{
   SCIP_Real oldval;

   assert(sol != NULL);
   assert(stat != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || sol->solorigin == SCIP_SOLORIGIN_PARTIAL
      || sol->solorigin == SCIP_SOLORIGIN_UNKNOWN
      || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(var != NULL);
   assert(SCIPisFinite(val));

   SCIPsetDebugMsg(set, "setting value of <%s> in solution %p to %g\n", SCIPvarGetName(var), (void*)sol, val);

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( SCIPsolIsOriginal(sol) )
      {
         oldval = solGetArrayVal(sol, var);

         if( !SCIPsetIsEQ(set, val, oldval) )
         {
            SCIP_Real obj;
            SCIP_Real objcont;

            SCIP_CALL( solSetArrayVal(sol, set, var, val) );

            /* update the objective value; we do not need to do this for partial solutions */
            if( !SCIPsolIsPartial(sol) )
            {
               /* an unknown solution value does not count towards the objective */
               obj = SCIPvarGetUnchangedObj(var);
               if( oldval != SCIP_UNKNOWN ) /*lint !e777*/
               {
                  objcont = obj * oldval;

                  /* we want to use a clean infinity */
                  if( SCIPsetIsInfinity(set, -objcont) || SCIPsetIsInfinity(set, sol->obj-objcont) )
                     sol->obj = SCIPsetInfinity(set);
                  else if( SCIPsetIsInfinity(set, objcont) || SCIPsetIsInfinity(set, -(sol->obj-objcont)) )
                     sol->obj = -SCIPsetInfinity(set);
                  else
                     sol->obj -= objcont;
               }
               if( val != SCIP_UNKNOWN ) /*lint !e777*/
               {
                  objcont = obj * val;

                  /* we want to use a clean infinity */
                  if( SCIPsetIsInfinity(set, objcont) || SCIPsetIsInfinity(set, sol->obj+objcont) )
                     sol->obj = SCIPsetInfinity(set);
                  else if( SCIPsetIsInfinity(set, objcont) || SCIPsetIsInfinity(set, -(sol->obj+objcont)) )
                     sol->obj = -SCIPsetInfinity(set);
                  else
                     sol->obj += objcont;
               }
            }

            solStamp(sol, stat, tree, FALSE);

         }
         return SCIP_OKAY;
      }
      else
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetTransVar(var), val);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(!SCIPsolIsOriginal(sol));
      assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
         || sol->lpcount == stat->lpcount);
      oldval = solGetArrayVal(sol, var);
      if( !SCIPsetIsEQ(set, val, oldval) )
      {
         SCIP_Real obj;
         SCIP_Real objcont;

         SCIP_CALL( solSetArrayVal(sol, set, var, val) );

         /* update objective: an unknown solution value does not count towards the objective */
         obj = SCIPvarGetUnchangedObj(var);

         if( oldval != SCIP_UNKNOWN ) /*lint !e777*/
         {
            objcont = obj * oldval;

            /* we want to use a clean infinity */
            if( SCIPsetIsInfinity(set, -objcont) || SCIPsetIsInfinity(set, sol->obj-objcont) )
               sol->obj = SCIPsetInfinity(set);
            else if( SCIPsetIsInfinity(set, objcont) || SCIPsetIsInfinity(set, -(sol->obj-objcont)) )
               sol->obj = -SCIPsetInfinity(set);
            else
               sol->obj -= objcont;
         }
         if( val != SCIP_UNKNOWN ) /*lint !e777*/
         {
            objcont = obj * val;

            /* we want to use a clean infinity */
            if( SCIPsetIsInfinity(set, objcont) || SCIPsetIsInfinity(set, sol->obj+objcont) )
               sol->obj = SCIPsetInfinity(set);
            else if( SCIPsetIsInfinity(set, -objcont) || SCIPsetIsInfinity(set, -(sol->obj+objcont)) )
               sol->obj = -SCIPsetInfinity(set);
            else
               sol->obj += objcont;
         }

         solStamp(sol, stat, tree, FALSE);
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(!SCIPsolIsOriginal(sol));
      oldval = SCIPvarGetLbGlobal(var);
      if( !SCIPsetIsEQ(set, val, oldval) )
      {
         SCIPerrorMessage("cannot set solution value for variable <%s> fixed to %.15g to different value %.15g\n",
            SCIPvarGetName(var), oldval, val);
         return SCIP_INVALIDDATA;
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, SCIPvarGetAggrScalar(var)));
      assert(!SCIPsetIsInfinity(set, SCIPvarGetAggrConstant(var)) && !SCIPsetIsInfinity(set, -SCIPvarGetAggrConstant(var)));
      assert(!SCIPsetIsInfinity(set, SCIPvarGetAggrScalar(var)) && !SCIPsetIsInfinity(set, -SCIPvarGetAggrScalar(var)));

      if( val == SCIP_UNKNOWN )/*lint !e777*/
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetAggrVar(var), val);
      if( SCIPsetIsInfinity(set, val) || SCIPsetIsInfinity(set, -val) )
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetAggrVar(var),  SCIPvarGetAggrScalar(var) > 0 ? val : -val);
      else
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetAggrVar(var), (val - SCIPvarGetAggrConstant(var))/SCIPvarGetAggrScalar(var));

   case SCIP_VARSTATUS_MULTAGGR:
      if ( SCIPvarGetMultaggrNVars(var) == 1 )
      {
         SCIP_VAR** multaggrvars;
         SCIP_Real* multaggrscalars;
         SCIP_Real multaggrconstant;

         multaggrvars = SCIPvarGetMultaggrVars(var);
         multaggrscalars = SCIPvarGetMultaggrScalars(var);
         multaggrconstant = SCIPvarGetMultaggrConstant(var);

         if( SCIPsetIsInfinity(set, multaggrconstant) || SCIPsetIsInfinity(set, -multaggrconstant) )
         {
            if( (SCIPsetIsInfinity(set, multaggrconstant) && !SCIPsetIsInfinity(set, val))
               || (SCIPsetIsInfinity(set, -multaggrconstant) && !SCIPsetIsInfinity(set, -val)) )
            {
               SCIPerrorMessage("cannot set solution value for variable <%s> fixed to %.15g to different value %.15g\n",
                  SCIPvarGetName(var), multaggrconstant, val);
               return SCIP_INVALIDDATA;
            }
            return SCIP_OKAY;
         }
         else
         {
            if( SCIPsetIsInfinity(set, val) || SCIPsetIsInfinity(set, -val) )
               return SCIPsolSetVal(sol, set, stat, tree, multaggrvars[0],  multaggrscalars[0] > 0 ? val : -val);
            else
               return SCIPsolSetVal(sol, set, stat, tree, multaggrvars[0], (val - multaggrconstant)/multaggrscalars[0]);
         }
      }
      SCIPerrorMessage("cannot set solution value for multiple aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      assert(!SCIPsetIsInfinity(set, SCIPvarGetNegationConstant(var)) && !SCIPsetIsInfinity(set, -SCIPvarGetNegationConstant(var)));

      if( val == SCIP_UNKNOWN )/*lint !e777*/
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), val);
      else if( SCIPsetIsInfinity(set, val) || SCIPsetIsInfinity(set, -val) )
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), -val);
      else
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), SCIPvarGetNegationConstant(var) - val);

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** increases value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolIncVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Real             incval              /**< increment for solution value of variable */
   )
{
   SCIP_Real oldval;

   assert(sol != NULL);
   assert(stat != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(var != NULL);
   assert(!SCIPsetIsInfinity(set, incval) && !SCIPsetIsInfinity(set, -incval));

   SCIPsetDebugMsg(set, "increasing value of <%s> in solution %p by %g\n", SCIPvarGetName(var), (void*)sol, incval);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
      || sol->lpcount == stat->lpcount);

   oldval = solGetArrayVal(sol, var);
   if( SCIPsetIsInfinity(set, oldval) || SCIPsetIsInfinity(set, -oldval) )
      return SCIP_OKAY;

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   /* @todo: handle strange cases, such as sums that yield infinite values */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( SCIPsolIsOriginal(sol) )
      {
         SCIP_CALL( solIncArrayVal(sol, set, var, incval) );
         sol->obj += SCIPvarGetUnchangedObj(var) * incval;
         solStamp(sol, stat, tree, FALSE);
         return SCIP_OKAY;
      }
      else
         return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetTransVar(var), incval);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(!SCIPsolIsOriginal(sol));
      SCIP_CALL( solIncArrayVal(sol, set, var, incval) );
      sol->obj += SCIPvarGetUnchangedObj(var) * incval;
      solStamp(sol, stat, tree, FALSE);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot increase solution value for fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, SCIPvarGetAggrScalar(var)));
      return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetAggrVar(var), incval/SCIPvarGetAggrScalar(var));

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot increase solution value for multiple aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), -incval);

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** returns value of variable in primal CIP solution */
SCIP_Real SCIPsolGetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* scalars;
   SCIP_Real solval;
   SCIP_Real solvalsum;
   int nvars;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || sol->solorigin == SCIP_SOLORIGIN_PARTIAL
      || sol->solorigin == SCIP_SOLORIGIN_UNKNOWN
      || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(var != NULL);

   /* if the value of a transformed variable in an original solution is requested, we need to project the variable back
    * to the original space, the opposite case is handled below
    */
   if( SCIPsolIsOriginal(sol) && SCIPvarIsTransformed(var) )
   {
      SCIP_RETCODE retcode;
      SCIP_VAR* origvar;
      SCIP_Real scalar;
      SCIP_Real constant;

      /* we cannot get the value of a transformed variable for a solution that lives in the original problem space
       * -> get the corresponding original variable first
       */
      origvar = var;
      scalar = 1.0;
      constant = 0.0;
      retcode = SCIPvarGetOrigvarSum(&origvar, &scalar, &constant);
      if ( retcode != SCIP_OKAY )
         return SCIP_INVALID;
      if( origvar == NULL )
      {
         /* the variable has no original counterpart: in the original solution, it has a value of zero */
         return 0.0;
      }
      assert(!SCIPvarIsTransformed(origvar));
      return scalar * SCIPsolGetVal(sol, set, stat, origvar) + constant;
   }

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( SCIPsolIsOriginal(sol) )
         return solGetArrayVal(sol, var);
      else
         return SCIPsolGetVal(sol, set, stat, SCIPvarGetTransVar(var));

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(!SCIPsolIsOriginal(sol));
      assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
         || sol->lpcount == stat->lpcount);
      return solGetArrayVal(sol, var);

   case SCIP_VARSTATUS_FIXED:
      assert(!SCIPsolIsOriginal(sol));
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetUbGlobal(var)); /*lint !e777*/
      assert(SCIPvarGetLbLocal(var) == SCIPvarGetUbLocal(var)); /*lint !e777*/
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetLbLocal(var)); /*lint !e777*/
      return SCIPvarGetLbGlobal(var);

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      solval = SCIPsolGetVal(sol, set, stat, SCIPvarGetAggrVar(var));
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         return SCIP_UNKNOWN;
      if( SCIPsetIsInfinity(set, solval) || SCIPsetIsInfinity(set, -solval) )
      {
         if( SCIPvarGetAggrScalar(var) * solval > 0.0 )
            return SCIPsetInfinity(set);
         if( SCIPvarGetAggrScalar(var) * solval < 0.0 )
            return -SCIPsetInfinity(set);
      }
      return SCIPvarGetAggrScalar(var) * solval + SCIPvarGetAggrConstant(var);

   case SCIP_VARSTATUS_MULTAGGR:
      nvars = SCIPvarGetMultaggrNVars(var);
      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      solvalsum = SCIPvarGetMultaggrConstant(var);
      for( i = 0; i < nvars; ++i )
      {
         solval = SCIPsolGetVal(sol, set, stat, vars[i]);
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            return SCIP_UNKNOWN;
         if( SCIPsetIsInfinity(set, solval) || SCIPsetIsInfinity(set, -solval) )
         {
            if( scalars[i] * solval > 0.0 )
               return SCIPsetInfinity(set);
            if( scalars[i] * solval < 0.0 )
               return -SCIPsetInfinity(set);
         }
         solvalsum += scalars[i] * solval;
      }
      return solvalsum;

   case SCIP_VARSTATUS_NEGATED:
      solval = SCIPsolGetVal(sol, set, stat, SCIPvarGetNegationVar(var));
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         return SCIP_UNKNOWN;
      if( SCIPsetIsInfinity(set, solval) )
         return -SCIPsetInfinity(set);
      if( SCIPsetIsInfinity(set, -solval) )
         return SCIPsetInfinity(set);
      return SCIPvarGetNegationConstant(var) - solval;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns value of variable in primal ray represented by primal CIP solution */
SCIP_Real SCIPsolGetRayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution, representing a primal ray */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* scalars;
   SCIP_Real solval;
   SCIP_Real solvalsum;
   int nvars;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO);
   assert(var != NULL);

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPsolGetRayVal(sol, set, stat, SCIPvarGetTransVar(var));

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return solGetArrayVal(sol, var);

   case SCIP_VARSTATUS_FIXED:
      assert(!SCIPsolIsOriginal(sol));
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetUbGlobal(var)); /*lint !e777*/
      assert(SCIPvarGetLbLocal(var) == SCIPvarGetUbLocal(var)); /*lint !e777*/
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetLbLocal(var)); /*lint !e777*/
      return 0.0; /* constants are ignored for computing the ray direction */

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      solval = SCIPsolGetRayVal(sol, set, stat, SCIPvarGetAggrVar(var));
      assert(solval != SCIP_UNKNOWN); /*lint !e777*/
      assert(!SCIPsetIsInfinity(set, REALABS(solval)));
      return SCIPvarGetAggrScalar(var) * solval; /* constants are ignored for computing the ray direction */

   case SCIP_VARSTATUS_MULTAGGR:
      nvars = SCIPvarGetMultaggrNVars(var);
      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      solvalsum = 0.0; /* constants are ignored for computing the ray direction */
      for( i = 0; i < nvars; ++i )
      {
         solval = SCIPsolGetRayVal(sol, set, stat, vars[i]);
         assert(solval != SCIP_UNKNOWN ); /*lint !e777*/
         assert(!SCIPsetIsInfinity(set, REALABS(solval)));
         solvalsum += scalars[i] * solval;
      }
      return solvalsum;

   case SCIP_VARSTATUS_NEGATED:
      solval = SCIPsolGetRayVal(sol, set, stat, SCIPvarGetNegationVar(var));
      assert(solval != SCIP_UNKNOWN); /*lint !e777*/
      assert(!SCIPsetIsInfinity(set, REALABS(solval)));
      return -solval; /* constants are ignored for computing the ray direction */

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** gets objective value of primal CIP solution in transformed problem */
SCIP_Real SCIPsolGetObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob            /**< original problem data */
   )
{
   assert(sol != NULL);

   /* for original solutions, sol->obj contains the external objective value */
   if( SCIPsolIsOriginal(sol) )
      return SCIPprobInternObjval(transprob, origprob, set, sol->obj);
   else
      return sol->obj;
}

/** updates primal solutions after a change in a variable's objective value */
void SCIPsolUpdateVarObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldobj,             /**< old objective value */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   SCIP_Real solval;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   solval = solGetArrayVal(sol, var);
   if( solval != SCIP_UNKNOWN ) /*lint !e777*/
      sol->obj += (newobj - oldobj) * solval;
}

/* mark the given solution as partial solution */
SCIP_RETCODE SCIPsolMarkPartial(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars               /**< number of problem variables */
   )
{
   SCIP_Real* vals;
   int v;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL);
   assert(nvars == 0 || vars != NULL);

   if( nvars == 0 )
      return SCIP_OKAY;;

   SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, nvars) );

   /* get values */
   for( v = 0; v < nvars; v++ )
   {
      assert(!SCIPvarIsTransformed(vars[v]));
      vals[v] = SCIPsolGetVal(sol, set, stat, vars[v]);
   }

   /* change origin to partial */
   sol->solorigin = SCIP_SOLORIGIN_PARTIAL;

   /* set values */
   for( v = 0; v < nvars; v++ )
   {
      int idx = SCIPvarGetIndex(vars[v]);

      if( vals[v] != SCIP_UNKNOWN ) /*lint !e777*/
      {
         /* from now on, variable must not be deleted */
         SCIPvarMarkNotDeletable(vars[v]);

         /* mark the variable valid */
         SCIP_CALL( SCIPboolarraySetVal(sol->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

         /* set the value in the solution array */
         SCIP_CALL( SCIPrealarraySetVal(sol->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, vals[v]) );
      }
      else
      {
         /* mark the variable invalid */
         SCIP_CALL( SCIPboolarraySetVal(sol->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, FALSE) );
      }
   }

   /* free buffer */
   SCIPsetFreeBufferArray(set, &vals);

   return SCIP_OKAY;
}

/** checks primal CIP solution for feasibility
 *
 *  @note The difference between SCIPsolCheck() and SCIPcheckSolOrig() is that modifiable constraints are handled
 *        differently. There might be some variables which do not have an original counter part (e.g. in
 *        branch-and-price). Therefore, modifiable constraints can not be double-checked in the original space.
 */
SCIP_RETCODE SCIPsolCheck(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            feasible            /**< stores whether solution is feasible */
   )
{
   SCIP_RESULT result;
   int h;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(set != NULL);
   assert(prob != NULL);
   assert(feasible != NULL);

   SCIPsetDebugMsg(set, "checking solution with objective value %g (nodenum=%" SCIP_LONGINT_FORMAT ", origin=%u)\n",
      sol->obj, sol->nodenum, sol->solorigin);

   *feasible = TRUE;

   SCIPsolResetViolations(sol);

   if( !printreason )
      completely = FALSE;

   /* check whether the solution respects the global bounds of the variables */
   if( checkbounds || sol->hasinfval )
   {
      int v;

      for( v = 0; v < prob->nvars && (*feasible || completely); ++v )
      {
         SCIP_VAR* var;
         SCIP_Real solval;

         var = prob->vars[v];
         solval = SCIPsolGetVal(sol, set, stat, var);

         if( solval != SCIP_UNKNOWN ) /*lint !e777*/
         {
            SCIP_Real lb;
            SCIP_Real ub;

            lb = SCIPvarGetLbGlobal(var);
            ub = SCIPvarGetUbGlobal(var);

            /* if we have to check bound and one of the current bounds is violated */
            if( checkbounds && ((!SCIPsetIsInfinity(set, -lb) && SCIPsetIsFeasLT(set, solval, lb))
                     || (!SCIPsetIsInfinity(set, ub) && SCIPsetIsFeasGT(set, solval, ub))) )
            {
               *feasible = FALSE;

               if( printreason )
               {
                  SCIPmessagePrintInfo(messagehdlr, "solution value %g violates bounds of <%s>[%g,%g] by %g\n", solval, SCIPvarGetName(var),
                        SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), MAX(lb - solval, 0.0) + MAX(solval - ub, 0.0));
               }
#ifdef SCIP_DEBUG
               else
               {
                  SCIPsetDebugMsgPrint(set, "  -> solution value %g violates bounds of <%s>[%g,%g]\n", solval, SCIPvarGetName(var),
                        SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
               }
#endif
            }

            /* check whether there are infinite variable values that lead to an objective value of +infinity */
            if( *feasible && sol->hasinfval )
            {
               *feasible = *feasible && (!SCIPsetIsInfinity(set, solval) || SCIPsetIsLE(set, SCIPvarGetUnchangedObj(var), 0.0) );
               *feasible = *feasible && (!SCIPsetIsInfinity(set, -solval) || SCIPsetIsGE(set, SCIPvarGetUnchangedObj(var), 0.0) );

               if( ((SCIPsetIsInfinity(set, solval) && SCIPsetIsGT(set, SCIPvarGetUnchangedObj(var), 0.0)) || (SCIPsetIsInfinity(set, -solval) && SCIPsetIsLT(set, SCIPvarGetUnchangedObj(var), 0.0))) )
               {
                  if( printreason )
                  {
                     SCIPmessagePrintInfo(messagehdlr, "infinite solution value %g for variable  <%s> with obj %g implies objective value +infinity\n",
                        solval, SCIPvarGetName(var), SCIPvarGetUnchangedObj(var));
                  }
#ifdef SCIP_DEBUG
                  else
                  {
                     SCIPsetDebugMsgPrint(set, "infinite solution value %g for variable  <%s> with obj %g implies objective value +infinity\n",
                        solval, SCIPvarGetName(var), SCIPvarGetUnchangedObj(var));
                  }
#endif
               }
            }
         }
      }
   }

   /* check whether the solution fulfills all constraints */
   for( h = 0; h < set->nconshdlrs && (*feasible || completely); ++h )
   {
      SCIP_CALL( SCIPconshdlrCheck(set->conshdlrs[h], blkmem, set, stat, sol,
            checkintegrality, checklprows, printreason, completely, &result) );
      *feasible = *feasible && (result == SCIP_FEASIBLE);

#ifdef SCIP_DEBUG
      if( !(*feasible) )
      {
         SCIPdebugPrintf("  -> infeasibility detected in constraint handler <%s>\n",
            SCIPconshdlrGetName(set->conshdlrs[h]));
      }
#endif
   }


   return SCIP_OKAY;
}

/** try to round given solution */
SCIP_RETCODE SCIPsolRound(
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool*            success             /**< pointer to store whether rounding was successful */
   )
{
   int nvars;
   int v;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(prob != NULL);
   assert(prob->transformed);
   assert(success != NULL);

   /* round all roundable fractional variables in the corresponding direction as long as no unroundable var was found */
   nvars = prob->nbinvars + prob->nintvars;
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real solval;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;

      var = prob->vars[v];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
         || sol->lpcount == stat->lpcount);
      solval = solGetArrayVal(sol, var);

      /* solutions with unknown entries cannot be rounded */
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         break;

      /* if solution value is already integral with feastol, continue */
      if( SCIPsetIsFeasIntegral(set, solval) )
         continue;

      /* get rounding possibilities */
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);

      /* choose rounding direction */
      if( mayrounddown && mayroundup )
      {
         /* we can round in both directions: round in objective function direction */
         if( SCIPvarGetUnchangedObj(var) >= 0.0 )
            solval = SCIPsetFeasFloor(set, solval);
         else
            solval = SCIPsetFeasCeil(set, solval);
      }
      else if( mayrounddown )
         solval = SCIPsetFeasFloor(set, solval);
      else if( mayroundup )
         solval = SCIPsetFeasCeil(set, solval);
      else
         break;

      /* store new solution value */
      SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, var, solval) );
   }

   /* check, if rounding was successful */
   *success = (v == nvars);

   return SCIP_OKAY;
}

/** updates the solution value sums in variables by adding the value in the given solution */
void SCIPsolUpdateVarsum(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Real             weight              /**< weight of solution in weighted average */
   )
{
   SCIP_Real solval;
   int v;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(0.0 <= weight && weight <= 1.0);

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      solval = SCIPsolGetVal(sol, set, stat, prob->vars[v]);
      if( solval != SCIP_UNKNOWN ) /*lint !e777*/
      {
         prob->vars[v]->primsolavg *= (1.0-weight);
         prob->vars[v]->primsolavg += weight*solval;
      }
   }
}

/** retransforms solution to original problem space */
SCIP_RETCODE SCIPsolRetransform(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_Bool*            hasinfval           /**< pointer to store whether the solution has infinite values */
   )
{
   SCIP_VAR** transvars;
   SCIP_VAR** vars;
   SCIP_VAR** activevars;
   SCIP_Real* solvals;
   SCIP_Real* activevals;
   SCIP_Real* transsolvals;
   SCIP_Real constant;
   int requiredsize;
   int ntransvars;
   int nactivevars;
   int nvars;
   int v;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO);
   assert(origprob != NULL);
   assert(transprob != NULL);
   assert(hasinfval != NULL);
   assert(!origprob->transformed);
   assert(transprob->transformed);

   *hasinfval = FALSE;

   /* This method was a performance bottleneck when retransforming a solution during presolving, before flattening the
    * aggregation graph. In that case, calling SCIPsolGetVal() on the original variable consumed too much
    * time. Therefore, we now first compute the active representation of each original variable using
    * SCIPvarGetActiveRepresentatives(), which is much faster, and sum up the solution values of the active variables by
    * hand for each original variable.
    */
   vars = origprob->vars;
   nvars = origprob->nvars;
   transvars = transprob->vars;
   ntransvars = transprob->nvars;

   /* allocate temporary memory for getting the active representation of the original variables, buffering the solution
    * values of all active variables and storing the original solution values
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &transsolvals, ntransvars + 1) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activevars, ntransvars + 1) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activevals, ntransvars + 1) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &solvals, nvars) );
   assert(transsolvals != NULL); /* for flexelint */
   assert(solvals != NULL); /* for flexelint */

   /* get the solution values of all active variables */
   for( v = 0; v < ntransvars; ++v )
   {
      transsolvals[v] = SCIPsolGetVal(sol, set, stat, transvars[v]);
   }

   /* get the solution in original problem variables */
   for( v = 0; v < nvars; ++v )
   {
      activevars[0] = vars[v];
      activevals[0] = 1.0;
      nactivevars = 1;
      constant = 0.0;

      /* get active representation of the original variable */
      SCIP_CALL( SCIPvarGetActiveRepresentatives(set, activevars, activevals, &nactivevars, ntransvars + 1, &constant,
            &requiredsize, TRUE) );
      assert(requiredsize <= ntransvars);

      /* compute solution value of the original variable */
      solvals[v] = constant;
      for( i = 0; i < nactivevars; ++i )
      {
         assert(0 <= SCIPvarGetProbindex(activevars[i]) && SCIPvarGetProbindex(activevars[i]) < ntransvars);
         assert(!SCIPsetIsInfinity(set, -solvals[v]) || !SCIPsetIsInfinity(set, activevals[i] * transsolvals[SCIPvarGetProbindex(activevars[i])]));
         assert(!SCIPsetIsInfinity(set, solvals[v]) || !SCIPsetIsInfinity(set, -activevals[i] * transsolvals[SCIPvarGetProbindex(activevars[i])]));
         solvals[v] += activevals[i] * transsolvals[SCIPvarGetProbindex(activevars[i])];
      }

      if( SCIPsetIsInfinity(set, solvals[v]) )
      {
         solvals[v] = SCIPsetInfinity(set);
         *hasinfval = TRUE;
      }
      else if( SCIPsetIsInfinity(set, -solvals[v]) )
      {
         solvals[v] = -SCIPsetInfinity(set);
         *hasinfval = TRUE;
      }
   }

   /* clear the solution and convert it into original space */
   SCIP_CALL( solClearArrays(sol) );
   sol->solorigin = SCIP_SOLORIGIN_ORIGINAL;
   sol->obj = origprob->objoffset;

   /* reinsert the values of the original variables */
   for( v = 0; v < nvars; ++v )
   {
      assert(SCIPvarGetUnchangedObj(vars[v]) == SCIPvarGetObj(vars[v])); /*lint !e777*/

      if( !SCIPsetIsZero(set, solvals[v]) )
      {
         SCIP_CALL( solSetArrayVal(sol, set, vars[v], solvals[v]) );
         if( solvals[v] != SCIP_UNKNOWN ) /*lint !e777*/
            sol->obj += SCIPvarGetUnchangedObj(vars[v]) * solvals[v];
      }
   }

   /**@todo remember the variables without original counterpart (priced variables) in the solution */

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &solvals);
   SCIPsetFreeBufferArray(set, &activevals);
   SCIPsetFreeBufferArray(set, &activevars);
   SCIPsetFreeBufferArray(set, &transsolvals);

   return SCIP_OKAY;
}

/** recomputes the objective value of an original solution, e.g., when transferring solutions
 *  from the solution pool (objective coefficients might have changed in the meantime)
 */
void SCIPsolRecomputeObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob            /**< original problem */
   )
{
   SCIP_VAR** vars;
   SCIP_Real solval;
   int nvars;
   int v;

   assert(sol != NULL);
   assert(SCIPsolIsOriginal(sol));
   assert(origprob != NULL);

   vars = origprob->vars;
   nvars = origprob->nvars;

   /* recompute the objective value */
   sol->obj = SCIPprobGetObjoffset(origprob);
   for( v = 0; v < nvars; ++v )
   {
      solval = SCIPsolGetVal(sol, set, stat, vars[v]);
      if( !SCIPsetIsZero(set, solval) && solval != SCIP_UNKNOWN ) /*lint !e777*/
      {
         sol->obj += SCIPvarGetUnchangedObj(vars[v]) * solval;
      }
   }

   if( SCIPsetIsInfinity(set, -sol->obj) )
      sol->obj = -SCIPsetInfinity(set);
}

/** returns whether the given solutions are equal */
SCIP_Bool SCIPsolsAreEqual(
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2,               /**< second primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob           /**< transformed problem after presolve, or NULL if both solution are
                                              *   defined in the original problem space */
   )
{
   SCIP_PROB* prob;
   SCIP_Real obj1;
   SCIP_Real obj2;
   int v;

   assert(sol1 != NULL);
   assert(sol2 != NULL);
   assert((SCIPsolIsOriginal(sol1) && SCIPsolIsOriginal(sol2)) || transprob != NULL);

   /* if both solutions are original or both are transformed, take the objective values stored in the solutions */
   if( SCIPsolIsOriginal(sol1) == SCIPsolIsOriginal(sol2) )
   {
      obj1 = sol1->obj;
      obj2 = sol2->obj;
   }
   /* one solution is original and the other not, so we have to get for both the objective in the transformed problem */
   else
   {
      obj1 = SCIPsolGetObj(sol1, set, transprob, origprob);
      obj2 = SCIPsolGetObj(sol2, set, transprob, origprob);
   }

   /* solutions with different objective values cannot be the same */
   if( !SCIPsetIsEQ(set, obj1, obj2) )
      return FALSE;

   /* if one of the solutions is defined in the original space, the comparison has to be performed in the original
    * space
    */
   prob = transprob;
   if( SCIPsolIsOriginal(sol1) || SCIPsolIsOriginal(sol2) )
      prob = origprob;
   assert(prob != NULL);

   /* compare each variable value */
   for( v = 0; v < prob->nvars; ++v )
   {
      SCIP_Real val1;
      SCIP_Real val2;

      val1 = SCIPsolGetVal(sol1, set, stat, prob->vars[v]);
      val2 = SCIPsolGetVal(sol2, set, stat, prob->vars[v]);
      if( !SCIPsetIsEQ(set, val1, val2) )
         return FALSE;
   }

   return TRUE;
}

/** outputs non-zero elements of solution to file stream */
SCIP_RETCODE SCIPsolPrint(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             mipstart,           /**< should only discrete variables be printed? */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Real solval;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(SCIPsolIsOriginal(sol) || prob->transformed || transprob != NULL);
   assert(!mipstart || !SCIPsolIsPartial(sol));

   /* display variables of problem data */
   for( v = 0; v < prob->nfixedvars; ++v )
   {
      assert(prob->fixedvars[v] != NULL);

      /* skip non-discrete variables in a mip start */
      if( mipstart && !SCIPvarIsIntegral(prob->fixedvars[v]) )
         continue;

      solval = SCIPsolGetVal(sol, set, stat, prob->fixedvars[v]);
      if( printzeros || mipstart
         || (sol->solorigin != SCIP_SOLORIGIN_PARTIAL && !SCIPsetIsZero(set, solval))
         || (sol->solorigin == SCIP_SOLORIGIN_PARTIAL && solval != SCIP_UNKNOWN) ) /*lint !e777*/
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->fixedvars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->fixedvars[v]));
      }
   }

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);

      /* skip non-discrete variables in a mip start */
      if( mipstart && !SCIPvarIsIntegral(prob->vars[v]) )
         continue;

      solval = SCIPsolGetVal(sol, set, stat, prob->vars[v]);
      if( printzeros || mipstart
         || (sol->solorigin != SCIP_SOLORIGIN_PARTIAL && !SCIPsetIsZero(set, solval))
         || (sol->solorigin == SCIP_SOLORIGIN_PARTIAL && solval != SCIP_UNKNOWN) ) /*lint !e777*/
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->vars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->vars[v]));
      }
   }

   /* display additional priced variables (if given problem data is original problem) */
   if( !prob->transformed && !SCIPsolIsOriginal(sol) )
   {
      assert(transprob != NULL);
      for( v = 0; v < transprob->nfixedvars; ++v )
      {
         assert(transprob->fixedvars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->fixedvars[v]) )
            continue;

         /* skip non-discrete variables in a mip start */
         if( mipstart && !SCIPvarIsIntegral(transprob->fixedvars[v]) )
            continue;

         solval = SCIPsolGetVal(sol, set, stat, transprob->fixedvars[v]);
         if( printzeros || mipstart || !SCIPsetIsZero(set, solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->fixedvars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->fixedvars[v]));
         }
      }
      for( v = 0; v < transprob->nvars; ++v )
      {
         assert(transprob->vars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->vars[v]) )
            continue;

         /* skip non-discrete variables in a mip start */
         if( mipstart && !SCIPvarIsIntegral(transprob->vars[v]) )
            continue;

         solval = SCIPsolGetVal(sol, set, stat, transprob->vars[v]);
         if( printzeros || !SCIPsetIsZero(set, solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->vars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->vars[v]));
         }
      }
   }

   return SCIP_OKAY;
}

/** outputs non-zero elements of solution representing a ray to file stream */
SCIP_RETCODE SCIPsolPrintRay(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Real solval;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(SCIPsolIsOriginal(sol) || prob->transformed || transprob != NULL);

   /* display variables of problem data */
   for( v = 0; v < prob->nfixedvars; ++v )
   {
      assert(prob->fixedvars[v] != NULL);
      solval = SCIPsolGetRayVal(sol, set, stat, prob->fixedvars[v]);
      if( printzeros || !SCIPsetIsZero(set, solval) )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->fixedvars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->fixedvars[v]));
      }
   }
   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      solval = SCIPsolGetRayVal(sol, set, stat, prob->vars[v]);
      if( printzeros || !SCIPsetIsZero(set, solval) )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->vars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->vars[v]));
      }
   }

   /* display additional priced variables (if given problem data is original problem) */
   if( !prob->transformed && !SCIPsolIsOriginal(sol) )
   {
      assert(transprob != NULL);
      for( v = 0; v < transprob->nfixedvars; ++v )
      {
         assert(transprob->fixedvars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->fixedvars[v]) )
            continue;

         solval = SCIPsolGetRayVal(sol, set, stat, transprob->fixedvars[v]);
         if( printzeros || !SCIPsetIsZero(set, solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->fixedvars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->fixedvars[v]));
         }
      }
      for( v = 0; v < transprob->nvars; ++v )
      {
         assert(transprob->vars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->vars[v]) )
            continue;

         solval = SCIPsolGetRayVal(sol, set, stat, transprob->vars[v]);
         if( printzeros || !SCIPsetIsZero(set, solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->vars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->vars[v]));
         }
      }
   }

   return SCIP_OKAY;
}

/*
 * methods for accumulated numerical violations of a solution
 */

/** reset violations of a solution */
void SCIPsolResetViolations(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   sol->viol.absviolbounds = 0.0;
   sol->viol.relviolbounds = 0.0;
   sol->viol.absviolintegrality = 0.0;
   sol->viol.absviollprows = 0.0;
   sol->viol.relviollprows = 0.0;
   sol->viol.absviolcons = 0.0;
   sol->viol.relviolcons = 0.0;
}

/** update integrality violation of a solution */
void SCIPsolUpdateIntegralityViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolintegrality  /**< absolute violation of integrality */
   )
{
   assert(sol != NULL);

   sol->viol.absviolintegrality = MAX(sol->viol.absviolintegrality, absviolintegrality);
}

/** update bound violation of a solution */
void SCIPsolUpdateBoundViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolbounds,      /**< absolute violation of bounds */
   SCIP_Real             relviolbounds       /**< relative violation of bounds */
   )
{
   assert(sol != NULL);

   sol->viol.absviolbounds = MAX(sol->viol.absviolbounds, absviolbounds);
   sol->viol.relviolbounds = MAX(sol->viol.relviolbounds, relviolbounds);
}

/** update LP row violation of a solution */
void SCIPsolUpdateLPRowViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviollprows,      /**< absolute violation of LP rows */
   SCIP_Real             relviollprows       /**< relative violation of LP rows */
   )
{
   assert(sol != NULL);

   sol->viol.absviollprows = MAX(sol->viol.absviollprows, absviollprows);
   sol->viol.relviollprows = MAX(sol->viol.relviollprows, relviollprows);
}

/** update constraint violation of a solution */
void SCIPsolUpdateConsViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolcons,        /**< absolute violation of constraint */
   SCIP_Real             relviolcons         /**< relative violation of constraint */
   )
{
   assert(sol != NULL);

   sol->viol.absviolcons = MAX(sol->viol.absviolcons, absviolcons);
   sol->viol.relviolcons = MAX(sol->viol.relviolcons, relviolcons);
}

/** update violation of a constraint that is represented in the LP */
void SCIPsolUpdateLPConsViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation of constraint */
   SCIP_Real             relviol             /**< relative violation of constraint */
   )
{
   assert(sol != NULL);

   SCIPsolUpdateConsViolation(sol, absviol, relviol);
   SCIPsolUpdateLPRowViolation(sol, absviol, relviol);
}

/** get maximum absolute bound violation of solution */
SCIP_Real SCIPsolGetAbsBoundViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.absviolbounds;
}

/** get maximum relative bound violation of solution */
SCIP_Real SCIPsolGetRelBoundViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.relviolbounds;
}

/** get maximum absolute integrality violation of solution */
SCIP_Real SCIPsolGetAbsIntegralityViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.absviolintegrality;
}

/** get maximum absolute LP row violation of solution */
SCIP_Real SCIPsolGetAbsLPRowViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.absviollprows;
}

/** get maximum relative LP row violation of solution */
SCIP_Real SCIPsolGetRelLPRowViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.relviollprows;
}

/** get maximum absolute constraint violation of solution */
SCIP_Real SCIPsolGetAbsConsViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.absviolcons;
}

/** get maximum relative constraint violation of solution */
SCIP_Real SCIPsolGetRelConsViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.relviolcons;
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPsolGetOrigin
#undef SCIPsolIsOriginal
#undef SCIPsolGetOrigObj
#undef SCIPsolGetTime
#undef SCIPsolGetNodenum
#undef SCIPsolGetRunnum
#undef SCIPsolGetDepth
#undef SCIPsolGetHeur
#undef SCIPsolOrigAddObjval
#undef SCIPsolGetPrimalIndex
#undef SCIPsolSetPrimalIndex
#undef SCIPsolGetIndex
#undef SCIPsolSetHeur

/** gets origin of solution */
SCIP_SOLORIGIN SCIPsolGetOrigin(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->solorigin;
}

/** returns whether the given solution is defined on original variables */
SCIP_Bool SCIPsolIsOriginal(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return (sol->solorigin == SCIP_SOLORIGIN_ORIGINAL || sol->solorigin == SCIP_SOLORIGIN_PARTIAL);
}

/** returns whether the given solution is defined on original variables and containes unknown solution values */
SCIP_Bool SCIPsolIsPartial(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return (sol->solorigin == SCIP_SOLORIGIN_PARTIAL);
}

/** gets objective value of primal CIP solution which lives in the original problem space */
SCIP_Real SCIPsolGetOrigObj(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);
   assert(SCIPsolIsOriginal(sol));

   return sol->obj;
}

/** adds value to the objective value of a given original primal CIP solution */
void SCIPsolOrigAddObjval(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             addval              /**< offset value to add */
   )
{
   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL);

   sol->obj += addval;
}

/** gets clock time, when this solution was found */
SCIP_Real SCIPsolGetTime(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->time;
}

/** gets branch and bound run number, where this solution was found */
int SCIPsolGetRunnum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->runnum;
}

/** gets node number, where this solution was found */
SCIP_Longint SCIPsolGetNodenum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->nodenum;
}

/** gets node's depth, where this solution was found */
int SCIPsolGetDepth(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->depth;
}

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
SCIP_HEUR* SCIPsolGetHeur(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->heur;
}

/** gets current position of solution in array of existing solutions of primal data */
int SCIPsolGetPrimalIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->primalindex;
}

/** sets current position of solution in array of existing solutions of primal data */
void SCIPsolSetPrimalIndex(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   primalindex         /**< new primal index of solution */
   )
{
   assert(sol != NULL);

   sol->primalindex = primalindex;
}

/** returns unique index of given solution */
int SCIPsolGetIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->index;
}

/** informs the solution that it now belongs to the given primal heuristic */
void SCIPsolSetHeur(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);

   sol->heur = heur;
}

