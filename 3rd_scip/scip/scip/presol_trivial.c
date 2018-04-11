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

/**@file   presol_trivial.c
 * @brief  trivial presolver: round fractional bounds on integer variables, fix variables with equal bounds
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_trivial.h"


#define PRESOL_NAME            "trivial"
#define PRESOL_DESC            "round fractional bounds on integers, fix variables with equal bounds"
#define PRESOL_PRIORITY        +9000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolver (fast, medium, or exhaustive) */

#ifdef FIXSIMPLEVALUE
#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */
#endif


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyTrivial)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolTrivial(scip) );

   return SCIP_OKAY;
}


/** presolving execution method */
static
SCIP_DECL_PRESOLEXEC(presolExecTrivial)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* get the problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* scan the variables for trivial bound reductions
    * (loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array)
    */
   for( v = nvars-1; v >= 0; --v )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      /* get variable's bounds */
      lb = SCIPvarGetLbGlobal(vars[v]);
      ub = SCIPvarGetUbGlobal(vars[v]);

      /* is variable integral? */
      if( SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real newlb;
         SCIP_Real newub;

         /* round fractional bounds on integer variables */
         newlb = SCIPfeasCeil(scip, lb);
         newub = SCIPfeasFloor(scip, ub);

         /* check bounds on variable for infeasibility */
         if( newlb > newub + 0.5 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
               "problem infeasible: integral variable <%s> has bounds [%.17f,%.17f] rounded to [%.17f,%.17f]\n",
               SCIPvarGetName(vars[v]), lb, ub, newlb, newub);
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* fix variables with equal bounds */
         if( newlb > newub - 0.5 )
         {
            SCIPdebugMsg(scip, "fixing integral variable <%s>: [%.17f,%.17f] -> [%.17f,%.17f]\n", SCIPvarGetName(vars[v]), lb, ub, newlb, newub);
            SCIP_CALL( SCIPfixVar(scip, vars[v], newlb, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, " -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            assert(fixed);
            (*nfixedvars)++;
         }
         else
         {
            /* round fractional bounds */
            if( !SCIPisFeasEQ(scip, lb, newlb) )
            {
               SCIPdebugMsg(scip, "rounding lower bound of integral variable <%s>: [%.17f,%.17f] -> [%.17f,%.17f]\n",
                  SCIPvarGetName(vars[v]), lb, ub, newlb, ub);
               SCIP_CALL( SCIPchgVarLb(scip, vars[v], newlb) );
               (*nchgbds)++;
            }
            if( !SCIPisFeasEQ(scip, ub, newub) )
            {
               SCIPdebugMsg(scip, "rounding upper bound of integral variable <%s>: [%.17f,%.17f] -> [%.17f,%.17f]\n",
                  SCIPvarGetName(vars[v]), newlb, ub, newlb, newub);
               SCIP_CALL( SCIPchgVarUb(scip, vars[v], newub) );
               (*nchgbds)++;
            }
         }
      }
      else
      {
         /* check bounds on continuous variable for infeasibility */
         if( SCIPisFeasGT(scip, lb, ub) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
               "problem infeasible: continuous variable <%s> has bounds [%.17f,%.17f]\n",
               SCIPvarGetName(vars[v]), lb, ub);
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* fix variables with equal bounds */
         if( SCIPisEQ(scip, lb, ub) )
         {
            SCIP_Real fixval;

#ifdef FIXSIMPLEVALUE
            fixval = SCIPselectSimpleValue(lb - 0.9 * SCIPepsilon(scip), ub + 0.9 * SCIPepsilon(scip), MAXDNOM);
#else
            fixval = (lb + ub)/2;
#endif
            SCIPdebugMsg(scip, "fixing continuous variable <%s>[%.17f,%.17f] to %.17f\n", SCIPvarGetName(vars[v]), lb, ub, fixval);
            SCIP_CALL( SCIPfixVar(scip, vars[v], fixval, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, " -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            assert(fixed);
            (*nfixedvars)++;
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the trivial presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolTrivial(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presolptr;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presolptr, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING, presolExecTrivial, NULL) );

   assert(presolptr != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presolptr, presolCopyTrivial) );

   return SCIP_OKAY;
}
