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

/**@file cons_disjunction.c
 * @brief  constraint handler for disjunction constraints
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_disjunction.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "disjunction"
#define CONSHDLR_DESC          "disjunction of constraints (or(cons1, cons2, ..., consn))"
#define CONSHDLR_ENFOPRIORITY   -950000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -900000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in
                                         *   (-1: no limit) */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_FAST
#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP


#define DEFAULT_ALWAYSBRANCH       TRUE /**< alawys perform branching if one of the constraints is violated, otherwise only if all integers are fixed */

/*
 * Data structures
 */

/** constraint data for disjunction constraints */
struct SCIP_ConsData
{
   SCIP_CONS**           conss;              /**< constraints in disjunction */
   SCIP_CONS*            relaxcons;          /**< a conjunction constraint containing the linear relaxation of the
                                              *   disjunction constraint, or NULL
                                              */
   int                   consssize;          /**< size of conss array */
   int                   nconss;             /**< number of constraints in disjunction */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             alwaysbranch;       /**< alawys perform branching if one of the constraints is violated, otherwise only if all integers are fixed */
};

/*
 * Local methods
 */

/** creates disjunction constraint data, captures initial constraints of disjunction */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   SCIP_CONS**           conss,              /**< initial constraint in disjunction */
   int                   nconss,             /**< number of initial constraints in disjunction */
   SCIP_CONS*            relaxcons           /**< a conjuction constraint containing the liner relaxation of the disjunction constraint, or NULL */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   if( nconss > 0 )
   {
      assert(conss != NULL);

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->conss, conss, nconss) );

      (*consdata)->consssize = nconss;
      (*consdata)->nconss = nconss;
      (*consdata)->relaxcons = relaxcons;

      /* we need to capture the constraints to avoid that SCIP deletes them since they are not (yet) added to the
       * problem
       */
      if( SCIPisTransformed(scip) )
      {
         SCIP_CALL( SCIPtransformConss(scip, nconss, (*consdata)->conss, (*consdata)->conss) );

         if( (*consdata)->relaxcons != NULL )
         {
            SCIP_CALL( SCIPtransformCons(scip, (*consdata)->relaxcons, &(*consdata)->relaxcons) );
         }
      }
      else
      {
         int c;

         for( c = 0; c < nconss; ++c )
         {
            assert(conss[c] != NULL);
            SCIP_CALL( SCIPcaptureCons(scip, conss[c]) );
         }

         if( (*consdata)->relaxcons != NULL )
         {
            SCIP_CALL( SCIPcaptureCons(scip, relaxcons) );
         }
      }
   }
   else
   {
      (*consdata)->conss = NULL;
      (*consdata)->consssize = 0;
      (*consdata)->nconss = 0;
      (*consdata)->relaxcons = NULL;
   }

   return SCIP_OKAY;
}

/** frees constraint data and releases all constraints in disjunction */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data */
   )
{
   int c;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release constraints */
   for( c = 0; c < (*consdata)->nconss; ++c )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*consdata)->conss[c]) );
   }

   /* release relaxation constraint */
   if( (*consdata)->relaxcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*consdata)->relaxcons) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->conss, (*consdata)->consssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** adds constraint to disjunction */
static
SCIP_RETCODE consdataAddCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_CONS*            cons                /**< constraint to add to the disjunction */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(cons != NULL);

   /* get memory for additional constraint */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &consdata->conss, &consdata->consssize, consdata->nconss+1) );
   assert(consdata->conss != NULL);
   assert(consdata->nconss < consdata->consssize);

   /* insert constraint in array */
   consdata->conss[consdata->nconss] = cons;
   consdata->nconss++;

   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPtransformCons(scip, consdata->conss[consdata->nconss - 1], &(consdata->conss[consdata->nconss - 1])));
   }
   else
   {
      /* capture constraint */
      SCIP_CALL( SCIPcaptureCons(scip, cons) );
   }

   return SCIP_OKAY;
}

/** branches on disjunctive constraint */
static
SCIP_RETCODE branchCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< active disjunction constraint */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   SCIP_NODE* child;
   SCIP_Real estimate;
   int nconss;
   int i;

   assert(result != NULL);

   /* cannot branch on modifiable constraint */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conss = consdata->conss;
   assert(conss != NULL);

   nconss = consdata->nconss;
   assert(nconss > 0);

   estimate =  SCIPgetLocalTransEstimate(scip);

   /* add all inactive constraints to local subproblem */
   for( i = 0; i < nconss; ++i )
   {
      /* create the branch-and-bound tree child nodes of the current node */
      SCIP_CALL( SCIPcreateChild(scip, &child, 0.0, estimate) );

      /* if disjunctive constraint needs to be checked, the upgraded constraint also needs to be checked */
      if( SCIPconsIsChecked(cons) )
      {
         SCIP_CALL( SCIPsetConsChecked(scip, conss[i], TRUE) );
      }

      /* add constraints to nodes */
      SCIP_CALL( SCIPaddConsNode(scip, child, conss[i], NULL) );

      /* remove disjunction constraint, from child node */
      SCIP_CALL( SCIPdelConsNode(scip, child, cons) );
   }

   SCIPdebugMsg(scip, "disjunction constraint <%s> branched %d childs\n", SCIPconsGetName(cons), nconss);

   /* reset constraint age */
   SCIP_CALL( SCIPresetConsAge(scip, cons) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** checks disjunction constraints if at least one is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< active disjunction constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   int nconss;
   int i;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conss = consdata->conss;
   assert(conss != NULL);

   nconss = consdata->nconss;
   assert(nconss > 0);

   *result = SCIP_INFEASIBLE;

   SCIPdeactivateSolViolationUpdates(scip);

   /* check all constraints */
   for( i = 0; i < nconss && *result != SCIP_FEASIBLE; ++i )
   {
      SCIP_CALL( SCIPcheckCons(scip, conss[i], sol, checkintegrality, checklprows, FALSE, result) );
      assert(*result == SCIP_FEASIBLE || *result == SCIP_INFEASIBLE);
   }

   SCIPactivateSolViolationUpdates(scip);

   if( *result == SCIP_INFEASIBLE )
   {
      if( sol != NULL )
         SCIPupdateSolConsViolation(scip, sol, 1.0, 1.0);

      if( printreason )
      {
         SCIPinfoMessage(scip, NULL, "constraint %s is violated, all sub-constraints in this disjunction are violated by this given solution\n", SCIPconsGetName(cons));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
      }
   }

   return SCIP_OKAY;
}

/** propagation method for disjunction constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< disjunctive constraint */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conss = consdata->conss;
   assert(conss != NULL);

   nconss = consdata->nconss;
   assert(nconss >= 1);

   for( c = 0; c < nconss; ++c )
   {
      /* if a constraint of the disjunction is already active, the disjunction is enforce by this constraint and
       * therefore redundant and can be locally deleted
       */
      if( SCIPconsIsActive(conss[c]) )
      {
         /* if we can globally delete the whole disjunctive constraint, because one constraint is already active, we
          * might need to update the check stage
          */
         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetNNodes(scip) == 0 )
         {
            /* if disjunctive constraint needs to be checked, the upgraded constraint also needs to be checked */
            if( SCIPconsIsChecked(cons) )
            {
             SCIP_CALL( SCIPsetConsChecked(scip, conss[c], TRUE) );
            }
         }

         (*ndelconss)++;
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         break;
      }
      /* if a sub-constraint is globally deleted, it means that this constraint is redundant and always fulfilled and
       * this makes also this disjunction redundant
       */
      else if( SCIPconsIsDeleted(conss[c]) )
      {
         (*ndelconss)++;
         SCIP_CALL( SCIPdelCons(scip, cons) );
         break;
      }
   }

   return SCIP_OKAY;
}

/** helper function to enforce constraints */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool branch;
   int c;

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   branch = SCIPgetNPseudoBranchCands(scip) == 0 || conshdlrdata->alwaysbranch;

   for( c = 0; c < nconss && *result != SCIP_BRANCHED; ++c )
   {
      /* check the disjunction */
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, FALSE, FALSE, result) );

      if( *result == SCIP_INFEASIBLE && branch )
      {
         SCIP_CALL( branchCons(scip, conss[c], result) );
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyDisjunction)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrDisjunction(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeDisjunction)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteDisjunction)
{  /*lint --e{715}*/
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /* get constraint data of source constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->conss, sourcedata->nconss, sourcedata->relaxcons) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   *infeasible = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* if we have a relaxation constraint and it is not active, then we add it locally */
      if( consdata->relaxcons != NULL && !SCIPconsIsActive(consdata->relaxcons) )
      {
         SCIP_CALL( SCIPaddConsLocal(scip, consdata->relaxcons, NULL) );
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpDisjunction)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr,  conss,  nconss,  result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxDisjunction)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr,  conss,  nconss,  result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsDisjunction)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr,  conss,  nconss,  result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckDisjunction)
{  /*lint --e{715}*/
   int c;

   *result = SCIP_FEASIBLE;

   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE || completely); ++c )
   {
      SCIP_RESULT tmpres;

      /* check the disjunction */
      SCIP_CALL( checkCons(scip, conss[c], sol, checkintegrality, checklprows, printreason, &tmpres) );
      assert(tmpres == SCIP_FEASIBLE || tmpres == SCIP_INFEASIBLE);

      if( tmpres == SCIP_INFEASIBLE )
         *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropDisjunction)
{  /*lint --e{715}*/
   int ndelconss;
   int c;

   ndelconss = 0;

   /* in probing mode we do not for deletable constraints */
   if( !SCIPinProbing(scip) )
   {
      for( c = 0; c < nconss; ++c )
      {
         /* propagate constraint */
         SCIP_CALL( propagateCons(scip, conss[c], &ndelconss) );
      }
   }

   /* adjust result code */
   if( ndelconss > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int oldndelconss;
   int c;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;
   oldndelconss = *ndelconss;

   /* all disjunction constraints with one constraint can be replaced with that corresponding constraint */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPconsIsModifiable(conss[c]) && consdata->nconss == 1 )
      {
         /* add constraint to the problem */
         if( !SCIPconsIsActive(consdata->conss[0]) )
         {
            SCIP_CONS* subcons = consdata->conss[0];

            /* if disjunctive constraint needs to be checked, the upgraded constraint also needs to be checked */
            if( SCIPconsIsChecked(conss[c]) )
            {
               SCIP_CALL( SCIPsetConsChecked(scip, subcons, TRUE) );
            }

            SCIP_CALL( SCIPaddCons(scip, subcons) );
         }

         /* remove disjunction constraint */
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );

         *result = SCIP_SUCCESS;

         continue;
      }

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, conss[c], ndelconss) );
   }

   if( *ndelconss > oldndelconss )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock sub constraints */
   for( c = 0; c < consdata->nconss; ++c )
   {
      SCIP_CALL( SCIPaddConsLocks(scip, consdata->conss[c], nlockspos, nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPinfoMessage(scip, file, "disjunction(");

   for( i = 0; i < consdata->nconss; ++i )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIP_CALL( SCIPprintCons(scip, consdata->conss[i], file) );
   }

   /* print relaxation */
   if( consdata->relaxcons != NULL )
   {
      SCIPinfoMessage(scip, file, ",, ");
      SCIP_CALL( SCIPprintCons(scip, consdata->relaxcons, file) );
   }

   SCIPinfoMessage(scip, file, ")");

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseDisjunction)
{  /*lint --e{715}*/
   SCIP_CONS** conss;
   SCIP_Bool relaxed = FALSE;
   int nconss;
   int sconss;
   char* token;
   char* saveptr;
   char* nexttokenstart;
   char* copystr;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);

   SCIPdebugMsg(scip, "parsing disjunction <%s>\n", name);

   *success = TRUE;

   /* allocate memory for constraint in disjunction, initial size is set to 10 */
   nconss = 0;
   sconss = 10;
   SCIP_CALL( SCIPallocBufferArray(scip, &conss, sconss) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &copystr, str, (int)strlen(str)+1) );

   /* find '(' at the beginning, string should start with 'disjunction(' */
   saveptr = strpbrk(copystr, "("); /*lint !e158*/

   if( saveptr == NULL )
   {
      SCIPdebugMsg(scip, "error parsing disjunctive constraint: \"%s\"\n", str);
      *success = FALSE;
      goto TERMINATE;
   }

   /* skip '(' */
   ++saveptr;
   /* remember token start position */
   nexttokenstart = saveptr;

   /* brackets '(' and ')' can exist co we check for them and the constraint delimeter */
   saveptr = strpbrk(saveptr, "(,");

   /* brackets '(' and ')' can exist in the rest of the string so we need to skip them to find the end of the first
    * sub-constraint marked by a ','
    */
   if( saveptr != NULL )
   {
      do
      {
         int bracketcounter = 0;

         if( *saveptr == '(' )
         {
            do
            {
               ++bracketcounter;
               ++saveptr;

               /* find last ending bracket */
               while( bracketcounter > 0 )
               {
                  saveptr = strpbrk(saveptr, "()");

                  if( saveptr != NULL )
                  {
                     if( *saveptr == '(' )
                        ++bracketcounter;
                     else
                        --bracketcounter;

                     ++saveptr;
                  }
                  else
                  {
                     SCIPdebugMsg(scip, "error parsing disjunctive constraint: \"%s\"\n", str);
                     *success = FALSE;
                     goto TERMINATE;
                  }
               }

               saveptr = strpbrk(saveptr, "(,");
            }
            while( saveptr != NULL && *saveptr == '(' );
         }

         /* we found a ',' so the end of the first sub-constraint is determined */
         if( saveptr != NULL )
         {
          assert(*saveptr == ',');

          /* resize constraint array if necessary */
          if( nconss == sconss )
          {
             sconss = SCIPcalcMemGrowSize(scip, nconss+1);
             assert(nconss < sconss);

             SCIP_CALL( SCIPreallocBufferArray(scip, &conss, sconss) );
          }

          assert(saveptr > nexttokenstart);

          /* extract token for parsing */
          SCIP_CALL( SCIPduplicateBufferArray(scip, &token, nexttokenstart, saveptr - nexttokenstart + 1) );
          token[saveptr - nexttokenstart] = '\0';

          SCIPdebugMsg(scip, "disjunctive parsing token(constraint): %s\n", token);

          /* parsing a constraint, part of the disjunction */
          SCIP_CALL( SCIPparseCons(scip, &(conss[nconss]), token, initial, separate, enforce, FALSE, propagate, TRUE, modifiable, dynamic, removable, stickingatnode, success) );

          SCIPfreeBufferArray(scip, &token);

          if( *success )
             ++nconss;
          else
          {
             SCIPdebugMsg(scip, "error parsing disjunctive constraint: \"%s\"\n", str);
             goto TERMINATE;
          }
          /* skip ',' delimeter */
          ++saveptr;
          /* remember token start position */
          nexttokenstart = saveptr;

          /* check if we found the last constraint, which is a conjunctive relaxation of the disjunction, and in the
           * CIP format marked by two consecutive ','
           */
          if( *nexttokenstart == ',' )
          {
             /* remember token start position */
             nexttokenstart = saveptr+1;

             relaxed = TRUE;
             break;
          }

          saveptr = strpbrk(saveptr, "(,");
         }
      }
      while( saveptr != NULL );
   }

   /* find end of disjunction constraint */
   saveptr = strrchr(nexttokenstart, ')');

   if( saveptr == NULL )
   {
      SCIPdebugMsg(scip, "error parsing disjunctive constraint: \"%s\"\n", str);
      *success = FALSE;
      goto TERMINATE;
   }
   /* parse last sub-constraint */
   else
   {
      /* resize constraint array if necessary */
      if( nconss == sconss )
      {
         ++sconss;
         SCIP_CALL( SCIPreallocBufferArray(scip, &conss, sconss) );
      }

      assert(saveptr > nexttokenstart);

      /* extract token for parsing */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &token, nexttokenstart, saveptr - nexttokenstart + 1) );
      token[saveptr - nexttokenstart] = '\0';

      SCIPdebugMsg(scip, "disjunctive parsing token(constraint): %s\n", token);

      /* parsing a constraint, part of the disjunction */
      SCIP_CALL( SCIPparseCons(scip, &(conss[nconss]), token, initial, separate, enforce, FALSE, propagate, TRUE, modifiable, dynamic, removable, stickingatnode, success) );

      if( *success )
         ++nconss;

      SCIPfreeBufferArray(scip, &token);
   }
   assert(nconss > 0 || !(*success));

   /* if parsing sub-constraints was fine, create the disjunctive constraint */
   if( *success )
   {
      /* create disjunctive constraint */
      SCIP_CALL( SCIPcreateConsDisjunction(scip, cons, name, relaxed ? nconss - 1: nconss, conss, relaxed ? conss[nconss - 1] : NULL,
         initial, enforce, check, local, modifiable, dynamic) );
   }

   /* free parsed constraints */
   for( --nconss; nconss >= 0; --nconss )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conss[nconss]) );
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &copystr);
   SCIPfreeBufferArray(scip, &conss);

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONS** sourceconss;
   SCIP_CONS** conss;
   int nconss;
   int c;

   *valid = TRUE;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   nconss = sourcedata->nconss;

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, nconss) );
   sourceconss = sourcedata->conss;

   /* copy each constraint one by one */
   for( c = 0; c < nconss && (*valid); ++c )
   {
      SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourceconss[c], &conss[c], SCIPconsGetHdlr(sourceconss[c]),
            varmap, consmap, SCIPconsGetName(sourceconss[c]),
            SCIPconsIsInitial(sourceconss[c]), SCIPconsIsSeparated(sourceconss[c]), SCIPconsIsEnforced(sourceconss[c]),
            SCIPconsIsChecked(sourceconss[c]), SCIPconsIsPropagated(sourceconss[c]),
            SCIPconsIsLocal(sourceconss[c]), SCIPconsIsModifiable(sourceconss[c]),
            SCIPconsIsDynamic(sourceconss[c]), SCIPconsIsRemovable(sourceconss[c]), SCIPconsIsStickingAtNode(sourceconss[c]),
            global, valid) );
      assert(!(*valid) || conss[c] != NULL);
   }

   if( *valid )
   {
      SCIP_CONS* sourcerelaxcons;
      SCIP_CONS* targetrelaxcons;

      sourcerelaxcons = sourcedata->relaxcons;
      targetrelaxcons = NULL;

      if( sourcerelaxcons != NULL )
      {
         SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourcerelaxcons, &targetrelaxcons, SCIPconsGetHdlr(sourcerelaxcons),
               varmap, consmap, SCIPconsGetName(sourcerelaxcons),
               SCIPconsIsInitial(sourcerelaxcons), SCIPconsIsSeparated(sourcerelaxcons), SCIPconsIsEnforced(sourcerelaxcons),
               SCIPconsIsChecked(sourcerelaxcons), SCIPconsIsPropagated(sourcerelaxcons),
               SCIPconsIsLocal(sourcerelaxcons), SCIPconsIsModifiable(sourcerelaxcons),
               SCIPconsIsDynamic(sourcerelaxcons), SCIPconsIsRemovable(sourcerelaxcons),
               SCIPconsIsStickingAtNode(sourcerelaxcons),
               global, valid) );
      }

      if( *valid )
      {
         if( name == NULL )
         {
            SCIP_CALL( SCIPcreateConsDisjunction(scip, cons, SCIPconsGetName(sourcecons), nconss, conss, targetrelaxcons,
                  initial, enforce, check, local, modifiable, dynamic) );
         }
         else
         {
            SCIP_CALL( SCIPcreateConsDisjunction(scip, cons, name, nconss, conss, targetrelaxcons,
                  initial, enforce, check, local, modifiable, dynamic) );
         }

         if( targetrelaxcons != NULL )
         {
            SCIP_CALL( SCIPreleaseCons(scip, &targetrelaxcons) );
         }
      }
   }

   /* release the copied constraints */
   for( c = (*valid ? c - 1 : c - 2); c >= 0; --c )
   {
      assert(conss[c] != NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &conss[c]) );
   }

   SCIPfreeBufferArray(scip, &conss);

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for disjunction constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrDisjunction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create disjunction constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpDisjunction, consEnfopsDisjunction, consCheckDisjunction, consLockDisjunction,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyDisjunction, consCopyDisjunction) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeDisjunction) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteDisjunction) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpDisjunction) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseDisjunction) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolDisjunction, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintDisjunction) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropDisjunction, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransDisjunction) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxDisjunction) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/alwaysbranch",
         "alawys perform branching if one of the constraints is violated, otherwise only if all integers are fixed",
         &conshdlrdata->alwaysbranch, FALSE, DEFAULT_ALWAYSBRANCH, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a disjunction constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsDisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nconss,             /**< number of initial constraints in disjunction */
   SCIP_CONS**           conss,              /**< initial constraint in disjunction */
   SCIP_CONS*            relaxcons,          /**< a conjunction constraint containing the linear relaxation of the disjunction constraint, or NULL */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic             /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the disjunction constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("disjunction constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conss, nconss, relaxcons) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, FALSE, enforce, check, FALSE,
         local, modifiable, dynamic, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a cumulative constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsDisjunction(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsDisjunction() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicDisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nconss,             /**< number of initial constraints in disjunction */
   SCIP_CONS**           conss,              /**< initial constraint in disjunction */
   SCIP_CONS*            relaxcons           /**< a conjunction constraint containing the linear relaxation of the disjunction constraint, or NULL */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsDisjunction(scip, cons, name, nconss, conss, relaxcons,
         TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


/** adds constraint to the disjunction of constraints */
SCIP_RETCODE SCIPaddConsElemDisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< disjunction constraint */
   SCIP_CONS*            addcons             /**< additional constraint in disjunction */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(addcons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a disjunction constraint\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataAddCons(scip, consdata, addcons) );

   return SCIP_OKAY;

}

