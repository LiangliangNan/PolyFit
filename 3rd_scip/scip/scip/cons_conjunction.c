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

/**@file   cons_conjunction.c
 * @brief  constraint handler for conjunction constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_conjunction.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "conjunction"
#define CONSHDLR_DESC          "conjunction of constraints"
#define CONSHDLR_ENFOPRIORITY   +900000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -900000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_FAST

/*
 * Data structures
 */

/** constraint data for conjunction constraints */
struct SCIP_ConsData
{
   SCIP_CONS**           conss;              /**< constraints in conjunction */
   int                   consssize;          /**< size of conss array */
   int                   nconss;             /**< number of constraints in conjunction */
};


/*
 * Local methods
 */

/** creates conjunction constraint data, captures initial constraints of conjunction */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   SCIP_CONS**           conss,              /**< initial constraint in conjunction */
   int                   nconss              /**< number of initial constraints in conjunction */
   )
{
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   if( nconss > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->conss, conss, nconss) );
      (*consdata)->consssize = nconss;
      (*consdata)->nconss = nconss;

      if( SCIPisTransformed(scip) )
      {
         SCIP_CALL( SCIPtransformConss(scip, nconss, (*consdata)->conss, (*consdata)->conss) );
      }
      else
      {
	 int c;

	 for( c = 0; c < nconss; ++c )
	 {
	    SCIP_CALL( SCIPcaptureCons(scip, conss[c]) );
	 }
      }
   }
   else
   {
      (*consdata)->conss = NULL;
      (*consdata)->consssize = 0;
      (*consdata)->nconss = 0;
   }

   return SCIP_OKAY;
}

/** frees constraint data and releases all constraints in conjunction */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data */
   )
{
   int c;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release constraints */
   for( c = 0; c < (*consdata)->nconss; ++c )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*consdata)->conss[c]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->conss, (*consdata)->consssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** adds constraint to conjunction */
static
SCIP_RETCODE consdataAddCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_CONS*            cons                /**< constraint to add to the conjunction */
   )
{
   assert(consdata != NULL);

   /* get memory for additional constraint */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &consdata->conss, &consdata->consssize, consdata->nconss+1) );
   assert(consdata->conss != NULL);
   assert(consdata->nconss < consdata->consssize);

   /* insert constraint in array */
   consdata->conss[consdata->nconss] = cons;
   consdata->nconss++;

   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPtransformCons(scip, consdata->conss[consdata->nconss - 1], &(consdata->conss[consdata->nconss - 1])) );
   }
   else
   {
      /* capture constraint */
      SCIP_CALL( SCIPcaptureCons(scip, cons) );
   }

   return SCIP_OKAY;
}

/** adds all constraints in conjunction constraints to the problem; disables unmodifiable conjunction constraints */
static
SCIP_RETCODE addAllConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< active conjunction constraints */
   int                   nconss,             /**< number of active conjunction constraints */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   int c;
   int i;

   assert(result != NULL);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* add all inactive constraints to local subproblem */
      for( i = 0; i < consdata->nconss; ++i )
      {
	 /* update check flag for sub constraints when upgrade takes place */
	 if( SCIPconsIsChecked(conss[c]) )
	 {
	    /* make sure, the constraint is checked for feasibility */
	    SCIP_CALL( SCIPsetConsChecked(scip, consdata->conss[i], TRUE) );
	 }

         if( !SCIPconsIsActive(consdata->conss[i]) )
         {
            SCIPdebugMsg(scip, "adding constraint <%s> from add conjunction <%s>\n",
               SCIPconsGetName(consdata->conss[i]), SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPaddConsLocal(scip, consdata->conss[i], NULL) );
            *result = SCIP_CONSADDED;
         }
      }

      /* disable conjunction constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** checks all constraints in conjunction constraints for feasibility */
static
SCIP_RETCODE checkAllConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< active conjunction constraints */
   int                   nconss,             /**< number of active conjunction constraints */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   int c;
   int i;

   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE || completely); ++c )
   {
      SCIP_RESULT subresult = SCIP_FEASIBLE;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* check all constraints */
      for( i = 0; i < consdata->nconss && subresult == SCIP_FEASIBLE; ++i )
      {
         SCIP_CALL( SCIPcheckCons(scip, consdata->conss[i], sol, checkintegrality, checklprows, printreason, &subresult) );
         assert(subresult == SCIP_FEASIBLE || subresult == SCIP_INFEASIBLE);
      }

      if( subresult == SCIP_INFEASIBLE )
      {
         /* mark solution as violated */
         *result = SCIP_INFEASIBLE;
         /* update constraint violation in solution */
         if ( sol != NULL )
            SCIPupdateSolConsViolation(scip, sol, 1.0, 1.0);
         if( printreason )
         {
            assert( 0 < i && i <= consdata->nconss );
            SCIPinfoMessage(scip, NULL, "Conjunction constraint %s is violated, at least the sub-constraint %s is violated by this given solution.\n",
               SCIPconsGetName(conss[c]), SCIPconsGetName(consdata->conss[i-1]));
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */


 /** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyConjunction)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrConjunction(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteConjunction)
{  /*lint --e{715}*/
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransConjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int c;

   /* create constraint data for target constraint */
   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );

   /* get constraint data of source constraint */
   sourcedata = SCIPconsGetData(sourcecons);

   if( sourcedata->nconss > 0 )
   {
      targetdata->consssize = sourcedata->nconss;
      targetdata->nconss = sourcedata->nconss;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &targetdata->conss, targetdata->consssize) );
      for( c = 0; c < sourcedata->nconss; ++c )
      {
         SCIP_CALL( SCIPtransformCons(scip, sourcedata->conss[c], &targetdata->conss[c]) );
      }
   }
   else
   {
      targetdata->conss = NULL;
      targetdata->consssize = 0;
      targetdata->nconss = 0;
   }

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   SCIP_CALL( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   SCIP_CALL( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   SCIP_CALL( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* check all constraints of the conjunction */
   SCIP_CALL( checkAllConss(scip, conss, nconss, sol, checkintegrality, checklprows, printreason, completely, result) );

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolConjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;
   int i;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* all constraints in a conjunction constraint of the global problem can be added directly to the problem and
    * removed from the conjunction constraint;
    * an unmodifiable conjunction constraint can be deleted
    */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* add all inactive constraints to the global problem */
      for( i = 0; i < consdata->nconss; ++i )
      {
	 /* update check flag for sub constraints when upgrade takes place */
	 if( SCIPconsIsChecked(conss[c]) )
	 {
	    /* make sure, the constraint is checked for feasibility */
	    SCIP_CALL( SCIPsetConsChecked(scip, consdata->conss[i], TRUE) );
	 }

         /* add constraint, if it is not active yet */
         if( !SCIPconsIsActive(consdata->conss[i]) )
         {
            SCIPdebugMsg(scip, "adding constraint <%s> from add conjunction <%s>\n",
               SCIPconsGetName(consdata->conss[i]), SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPaddCons(scip, consdata->conss[i]) );
            *result = SCIP_SUCCESS;
         }
         /* release constraint because it will be removed from the conjunction constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &(consdata->conss[i])) );
      }
      /* all constraints where removed, so we need to clear the array */
      consdata->nconss = 0;

      /* delete conjunction constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockConjunction)
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
SCIP_DECL_CONSPRINT(consPrintConjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPinfoMessage(scip, file, "conjunction(");

   for( i = 0; i < consdata->nconss; ++i )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIP_CALL( SCIPprintCons(scip, consdata->conss[i], file) );
   }
   SCIPinfoMessage(scip, file, ")");

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseConjunction)
{  /*lint --e{715}*/
   SCIP_CONS** conss;
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

   SCIPdebugMsg(scip, "parsing conjunction <%s>\n", name);

   *success = TRUE;

   /* allocate memory for constraint in conjunction, initial size is set to 10 */
   nconss = 0;
   sconss = 10;
   SCIP_CALL( SCIPallocBufferArray(scip, &conss, sconss) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &copystr, str, (int)strlen(str)+1) );

   /* find '(' at the beginning, string should start with 'conjunction(' */
   saveptr = strpbrk(copystr, "("); /*lint !e158*/

   if( saveptr == NULL )
   {
      SCIPdebugMsg(scip, "error parsing conjunctive constraint: \"%s\"\n", str);
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
		     SCIPdebugMsg(scip, "error parsing conjunctive constraint: \"%s\"\n", str);
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

	    SCIPdebugMsg(scip, "conjunctive parsing token(constraint): %s\n", token);

	    /* parsing a constraint, part of the conjunction */
	    SCIP_CALL( SCIPparseCons(scip, &(conss[nconss]), token, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, success) );

	    SCIPfreeBufferArray(scip, &token);

	    if( *success )
	       ++nconss;
	    else
	    {
	       SCIPdebugMsg(scip, "error parsing conjunctive constraint: \"%s\"\n", str);
	       goto TERMINATE;
	    }
	    /* skip ',' delimeter */
	    ++saveptr;
	    /* remember token start position */
	    nexttokenstart = saveptr;

	    saveptr = strpbrk(saveptr, "(,");
	 }
      }
      while( saveptr != NULL );
   }

   /* find end of conjunction constraint */
   saveptr = strrchr(nexttokenstart, ')');

   if( saveptr == NULL )
   {
      SCIPdebugMsg(scip, "error parsing conjunctive constraint: \"%s\"\n", str);
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

      SCIPdebugMsg(scip, "conjunctive parsing token(constraint): %s\n", token);

      /* parsing a constraint, part of the conjunction */
      SCIP_CALL( SCIPparseCons(scip, &(conss[nconss]), token, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, success) );

      if( *success )
	 ++nconss;

      SCIPfreeBufferArray(scip, &token);
   }
   assert(nconss > 0 || !(*success));

   /* if parsing sub-constraints was fine, create the conjunctive constraint */
   if( *success )
   {
      /* create conjunctive constraint */
      SCIP_CALL( SCIPcreateConsConjunction(scip, cons, name, nconss, conss,
	    enforce, check, local, modifiable, dynamic) );
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
SCIP_DECL_CONSCOPY(consCopyConjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONS** sourceconss;
   SCIP_CONS** conss;
   int nconss;
   int c;

   *valid = TRUE;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   sourceconss = sourcedata->conss;
   nconss = sourcedata->nconss;

   if( nconss > 0 )
   {
      assert(sourceconss != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &conss, nconss) );

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
         if( name == NULL )
         {
            SCIP_CALL( SCIPcreateConsConjunction(scip, cons, SCIPconsGetName(sourcecons), nconss, conss,
                  enforce, check, local, modifiable, dynamic) );
         }
         else
         {
            SCIP_CALL( SCIPcreateConsConjunction(scip, cons, name, nconss, conss,
                  enforce, check, local, modifiable, dynamic) );
         }
      }

      /* release the copied constraints */
      for( c = (*valid ? c - 1 : c - 2); c >= 0; --c )
      {
         assert(conss[c] != NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &conss[c]) );
      }

      SCIPfreeBufferArray(scip, &conss);
   }

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for conjunction constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrConjunction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpConjunction, consEnfopsConjunction, consCheckConjunction, consLockConjunction,
         NULL) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyConjunction, consCopyConjunction) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteConjunction) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseConjunction) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolConjunction, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintConjunction) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransConjunction) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxConjunction) );

   return SCIP_OKAY;
}

/** creates and captures a conjunction constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsConjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nconss,             /**< number of initial constraints in conjunction */
   SCIP_CONS**           conss,              /**< initial constraint in conjunction */
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

   /* find the conjunction constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("conjunction constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conss, nconss) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, enforce, check, FALSE,
         local, modifiable, dynamic, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures an and constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsConjunction(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsConjunction() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicConjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nconss,             /**< number of initial constraints in conjunction */
   SCIP_CONS**           conss               /**< initial constraint in conjunction */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsConjunction(scip, cons, name, nconss, conss,
         TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** adds constraint to the conjunction of constraints */
SCIP_RETCODE SCIPaddConsElemConjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< conjunction constraint */
   SCIP_CONS*            addcons             /**< additional constraint in conjunction */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(addcons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a conjunction constraint\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataAddCons(scip, consdata, addcons) );

   return SCIP_OKAY;
}
