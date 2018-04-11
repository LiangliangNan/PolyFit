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

/**@file   cons_or.c
 * @brief  Constraint handler for "or" constraints,  \f$r = x_1 \vee x_2 \vee \dots  \vee x_n\f$
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 * This constraint handler deals with "or" constraint. These are constraint of the form:
 *
 * \f[
 *    r = x_1 \vee x_2 \vee \dots  \vee x_n
 * \f]
 *
 * where \f$x_i\f$ is a binary variable for all \f$i\f$. Hence, \f$r\f$ is also of binary type. The variable \f$r\f$ is
 * called resultant and the \f$x\f$'s operators.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_or.h"
#include "scip/cons_and.h"
#include "scip/pub_misc.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "or"
#define CONSHDLR_DESC          "constraint handler for or constraints: r = or(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -850000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -850000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */

#define EVENTHDLR_NAME         "or"
#define EVENTHDLR_DESC         "event handler for or constraints"


/*
 * Data structures
 */

/** constraint data for or constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables in the or operation */
   SCIP_VAR*             resvar;             /**< resultant variable */
   SCIP_ROW**            rows;               /**< rows for linear relaxation of or constraint */
   int                   nvars;              /**< number of variables in or operation */
   int                   varssize;           /**< size of vars array */
   int                   rowssize;           /**< size of rows array */
   int                   watchedvar1;        /**< position of first watched operator variable */
   int                   watchedvar2;        /**< position of second watched operator variable */
   int                   filterpos1;         /**< event filter position of first watched operator variable */
   int                   filterpos2;         /**< event filter position of second watched operator variable */
   unsigned int          propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int          nofixedone:1;       /**< is none of the operator variables fixed to TRUE? */
   unsigned int          impladded:1;        /**< were the implications of the constraint already added? */
   unsigned int          opimpladded:1;      /**< was the implication for 2 operands with fixed resultant added? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for events on watched variables */
};


/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1,                          /**< v_i = TRUE                                   =>  r   = TRUE            */
   PROPRULE_2,                          /**< r   = FALSE                                  =>  v_i = FALSE for all i */
   PROPRULE_3,                          /**< v_i = FALSE for all i                        =>  r   = FALSE           */
   PROPRULE_4,                          /**< r   = TRUE, v_i = FALSE for all i except j   =>  v_j = TRUE            */
   PROPRULE_INVALID                     /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;


/*
 * Local methods
 */

/** installs rounding locks for the given variable in the given or constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< or constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given or constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< or constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** creates constraint handler data */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, conshdlrdata) );

   /* set event handler for catching events on watched variables */
   (*conshdlrdata)->eventhdlr = eventhdlr;

   return SCIP_OKAY;
}

/** frees constraint handler data */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** gets number of LP rows needed for the LP relaxation of the constraint */
static
int consdataGetNRows(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   return consdata->nvars + 1;
}

/** catches events for the watched variable at given position */
static
SCIP_RETCODE consdataCatchWatchedEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< or constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos,                /**< array position of variable to catch bound change events for */
   int*                  filterpos           /**< pointer to store position of event filter entry */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(filterpos != NULL);

   /* catch tightening events for upper bound and relaxed events for lower bounds on watched variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, filterpos) );

   return SCIP_OKAY;
}


/** drops events for the watched variable at given position */
static
SCIP_RETCODE consdataDropWatchedEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< or constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos,                /**< array position of watched variable to drop bound change events for */
   int                   filterpos           /**< position of event filter entry */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(filterpos >= 0);

   /* drop tightening events for upper bound and relaxed events for lower bounds on watched variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, filterpos) );

   return SCIP_OKAY;
}

/** catches needed events on all variables of constraint, except the special ones for watched variables */
static
SCIP_RETCODE consdataCatchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< or constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* catch bound change events for both bounds on resultant variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDCHANGED, 
         eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

   /* catch tightening events for lower bound and relaxed events for upper bounds on operator variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
   }

   return SCIP_OKAY;
}

/** drops events on all variables of constraint, except the special ones for watched variables */
static
SCIP_RETCODE consdataDropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< or constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* drop bound change events for both bounds on resultant variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDCHANGED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* drop tightening events for lower bound and relaxed events for upper bounds on operator variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );
   }

   return SCIP_OKAY;
}

/** stores the given variable numbers as watched variables, and updates the event processing */
static
SCIP_RETCODE consdataSwitchWatchedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< or constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   watchedvar1,        /**< new first watched variable */
   int                   watchedvar2         /**< new second watched variable */
   )
{
   assert(consdata != NULL);
   assert(watchedvar1 == -1 || watchedvar1 != watchedvar2);
   assert(watchedvar1 != -1 || watchedvar2 == -1);
   assert(watchedvar1 == -1 || (0 <= watchedvar1 && watchedvar1 < consdata->nvars));
   assert(watchedvar2 == -1 || (0 <= watchedvar2 && watchedvar2 < consdata->nvars));

   /* if one watched variable is equal to the old other watched variable, just switch positions */
   if( watchedvar1 == consdata->watchedvar2 || watchedvar2 == consdata->watchedvar1 )
   {
      int tmp;

      tmp = consdata->watchedvar1;
      consdata->watchedvar1 = consdata->watchedvar2;
      consdata->watchedvar2 = tmp;
      tmp = consdata->filterpos1;
      consdata->filterpos1 = consdata->filterpos2;
      consdata->filterpos2 = tmp;
   }
   assert(watchedvar1 == -1 || watchedvar1 != consdata->watchedvar2);
   assert(watchedvar2 == -1 || watchedvar2 != consdata->watchedvar1);

   /* drop events on old watched variables */
   if( consdata->watchedvar1 != -1 && consdata->watchedvar1 != watchedvar1 )
   {
      assert(consdata->filterpos1 != -1);
      SCIP_CALL( consdataDropWatchedEvents(scip, consdata, eventhdlr, consdata->watchedvar1, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( consdataDropWatchedEvents(scip, consdata, eventhdlr, consdata->watchedvar2, consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      SCIP_CALL( consdataCatchWatchedEvents(scip, consdata, eventhdlr, watchedvar1, &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      SCIP_CALL( consdataCatchWatchedEvents(scip, consdata, eventhdlr, watchedvar2, &consdata->filterpos2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;

   return SCIP_OKAY;
}

/** ensures, that the vars array can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);

   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}

/** creates constraint data for or constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   nvars,              /**< number of variables in the or operation */
   SCIP_VAR**            vars,               /**< variables in or operation */
   SCIP_VAR*             resvar              /**< resultant variable */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(resvar != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
   (*consdata)->resvar = resvar;
   (*consdata)->rows = NULL;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->rowssize = 0;
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;
   (*consdata)->propagated = FALSE;
   (*consdata)->nofixedone = FALSE;
   (*consdata)->impladded = FALSE;
   (*consdata)->opimpladded = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->resvar, &(*consdata)->resvar) );

      /* catch needed events on variables */
      SCIP_CALL( consdataCatchEvents(scip, *consdata, eventhdlr) );
   }

   return SCIP_OKAY;
}

/** releases LP rows of constraint data and frees rows array */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int r;

   assert(consdata != NULL);

   if( consdata->rows != NULL )
   {
      int nrows;

      nrows = consdataGetNRows(consdata);

      for( r = 0; r < nrows; ++r )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->rows, consdata->rowssize);
   }

   return SCIP_OKAY;
}

/** frees constraint data for or constraint */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to the constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( SCIPisTransformed(scip) )
   {
      /* drop events for watched variables */
      SCIP_CALL( consdataSwitchWatchedvars(scip, *consdata, eventhdlr, -1, -1) );

      /* drop all other events on variables */
      SCIP_CALL( consdataDropEvents(scip, *consdata, eventhdlr) );
   }
   else
   {
      assert((*consdata)->watchedvar1 == -1);
      assert((*consdata)->watchedvar2 == -1);
   }

   /* release and free the rows */
   SCIP_CALL( consdataFreeRows(scip, *consdata) );

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints or constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< or constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(consdata != NULL);

   /* print resultant */
   SCIP_CALL( SCIPwriteVarName(scip, file, consdata->resvar, TRUE) );

   /* start the variable list */
   SCIPinfoMessage(scip, file, " == or(");

   /* print variable list */
   SCIP_CALL( SCIPwriteVarsList(scip, file, consdata->vars, consdata->nvars, TRUE, ',') );

   /* close the variable list */
   SCIPinfoMessage(scip, file, ")");

   return SCIP_OKAY;
}

/** adds coefficient in or constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows == NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1) );
   consdata->vars[consdata->nvars] = var;
   consdata->nvars++;

   /* if we are in transformed problem, catch the variable's events */
   if( transformed )
   {
      /* catch bound change events of variable */
      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
   }

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockRounding(scip, cons, var) );

   /**@todo update LP rows */
   if( consdata->rows != NULL )
   {
      SCIPerrorMessage("cannot add coefficients to or constraint after LP relaxation was created\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from or constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< or constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(consdata->vars[pos]));

   /* remove the rounding locks of the variable */
   SCIP_CALL( unlockRounding(scip, cons, consdata->vars[pos]) );

   if( SCIPconsIsTransformed(cons) )
   {
      /* drop bound change events of variable */
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );
   }

   if( SCIPconsIsTransformed(cons) )
   {
      /* if the position is watched, stop watching the position */
      if( consdata->watchedvar1 == pos )
      {
         SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, consdata->watchedvar2, -1) );
      }
      if( consdata->watchedvar2 == pos )
      {
         SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, consdata->watchedvar1, -1) );
      }
   }
   assert(pos != consdata->watchedvar1);
   assert(pos != consdata->watchedvar2);

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->nvars--;

   /* if the last variable (that moved) was watched, update the watched position */
   if( consdata->watchedvar1 == consdata->nvars )
      consdata->watchedvar1 = pos;
   if( consdata->watchedvar2 == consdata->nvars )
      consdata->watchedvar2 = pos;

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}

/** deletes all zero-fixed variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< or constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else
      {
         SCIP_VAR* repvar;
         SCIP_Bool negated;

         /* get binary representative of variable */
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );

         /* check, if the variable should be replaced with the representative */
         if( repvar != var )
         {
            /* delete old (aggregated) variable */
            SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );

            /* add representative instead */
            SCIP_CALL( addCoef(scip, cons, eventhdlr, repvar) );
         }
         else
            ++v;
      }
   }

   SCIPdebugMsg(scip, "after fixings: ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL)) );
   SCIPdebugMsgPrint(scip, "\n");

   return SCIP_OKAY;
}

/** creates LP rows corresponding to or constraint:
 *   - for each operator variable vi:  resvar - vi            >= 0
 *   - one additional row:             resvar - v1 - ... - vn <= 0
 */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   char rowname[SCIP_MAXSTRLEN];
   int nvars;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows == NULL);

   nvars = consdata->nvars;

   /* get memory for rows */
   consdata->rowssize = consdataGetNRows(consdata);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->rows, consdata->rowssize) );
   assert(consdata->rowssize == nvars+1);

   /* create operator rows */
   for( i = 0; i < nvars; ++i )
   {
      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), i);
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[i], SCIPconsGetHdlr(cons), rowname, 0.0, SCIPinfinity(scip),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[i], consdata->resvar, 1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[i], consdata->vars[i], -1.0) );
   }

   /* create additional row */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_add", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[nvars], SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), 0.0,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[nvars], consdata->resvar, 1.0) );
   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[nvars], nvars, consdata->vars, -1.0) );

   return SCIP_OKAY;
}

/** adds linear relaxation of or constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   SCIP_CONSDATA* consdata;
   int r;
   int nrows;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert( consdata->rows != NULL );

   nrows = consdataGetNRows(consdata);

   for( r = 0; r < nrows && !(*infeasible); ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
      {
         SCIP_CALL( SCIPaddRow(scip, consdata->rows[r], FALSE, infeasible) );
      }
   }

   return SCIP_OKAY;
}

/** checks or constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool mustcheck;
   int r;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   /* check, if we can skip this feasibility check, because all rows are in the LP and doesn't have to be checked */
   mustcheck = checklprows;
   mustcheck = mustcheck || (consdata->rows == NULL);
   if( !mustcheck )
   {
      int nrows;

      assert(consdata->rows != NULL);

      nrows = consdataGetNRows(consdata);

      for( r = 0; r < nrows; ++r )
      {
         mustcheck = !SCIProwIsInLP(consdata->rows[r]);
         if( mustcheck )
            break;
      }         
   }

   /* check feasibility of constraint if necessary */
   if( mustcheck )
   {
      SCIP_Real solval;
      int i;

      /* increase age of constraint; age is reset to zero, if a violation was found only in case we are in
       * enforcement
       */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }

      /* check, if all operator variables are FALSE */
      for( i = 0; i < consdata->nvars; ++i )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
         assert(SCIPisFeasIntegral(scip, solval));
         if( solval > 0.5 )
            break;
      }

      /* if all operator variables are FALSE, the resultant has to be FALSE, otherwise, the resultant has to be TRUE */
      solval = SCIPgetSolVal(scip, sol, consdata->resvar);
      assert(SCIPisFeasIntegral(scip, solval));

      if( (i == consdata->nvars) != (solval < 0.5) )
      {
         *violated = TRUE;

         /* only reset constraint age if we are in enforcement */
         if( sol == NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
         /* update constraint violation in solution */
         else
            SCIPupdateSolConsViolation(scip, sol, 1.0, 1.0);

         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, ";\nviolation:\n");
            if( i == consdata->nvars )
            {
               SCIPinfoMessage(scip, NULL, " all operands are FALSE and resultant <%s> = TRUE\n",
                  SCIPvarGetName(consdata->resvar)); 
            }
            else
            {
               SCIPinfoMessage(scip, NULL, " operand <%s> = TRUE and resultant <%s> = FALSE\n",
                  SCIPvarGetName(consdata->vars[i-1]), SCIPvarGetName(consdata->resvar)); 
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** separates current LP solution */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real feasibility;
   int r;
   int nrows;

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *separated = FALSE;

   /* create all necessary rows for the linear relaxation */
   if( consdata->rows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows != NULL);

   nrows = consdataGetNRows(consdata);

   /* test all rows for feasibility and add infeasible rows */
   for( r = 0; r < nrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
      {
         feasibility = SCIPgetRowSolFeasibility(scip, consdata->rows[r], sol);
         if( SCIPisFeasNegative(scip, feasibility) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPaddRow(scip, consdata->rows[r], FALSE, &infeasible) );
            assert( ! infeasible );
            *separated = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflicting FALSE assignment to resultant of given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflictZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< or constraint that detected the conflict */
   int                   truepos             /**< position of operand that is fixed to TRUE */
   )
{
   SCIP_CONSDATA* consdata;

   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetUbLocal(consdata->resvar) < 0.5);
   assert(0 <= truepos && truepos < consdata->nvars);
   assert(SCIPvarGetLbLocal(consdata->vars[truepos]) > 0.5);

   /* initialize conflict analysis, and add resultant and single operand variable to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[truepos]) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** analyzes conflicting TRUE assignment to resultant of given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflictOne(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< or constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   assert(!SCIPconsIsModifiable(cons));

   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetLbLocal(consdata->resvar) > 0.5);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(SCIPvarGetUbLocal(consdata->vars[v]) < 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagates constraint with the following rules:
 *   (1) v_i = TRUE                                   =>  r   = TRUE
 *   (2) r   = FALSE                                  =>  v_i = FALSE for all i
 *   (3) v_i = FALSE for all i                        =>  r   = FALSE
 *   (4) r   = TRUE, v_i = FALSE for all i except j   =>  v_j = TRUE
 */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< or constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* resvar;
   SCIP_VAR** vars;
   int nvars;
   int watchedvar1;
   int watchedvar2;
   int i;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(cutoff != NULL);
   assert(nfixedvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   resvar = consdata->resvar;
   vars = consdata->vars;
   nvars = consdata->nvars;

   /* don't process the constraint, if none of the operator variables was fixed to TRUE, and if the watched variables
    * and the resultant weren't fixed to any value since last propagation call
    */
   if( consdata->propagated )
   {
      assert(consdata->nofixedone);
      assert(SCIPisFeasEQ(scip, SCIPvarGetUbLocal(resvar), 1.0));
      return SCIP_OKAY;
   }

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   /* if one of the operator variables was fixed to TRUE, the resultant can be fixed to TRUE (rule (1)) */
   if( !consdata->nofixedone )
   {
      for( i = 0; i < nvars && SCIPvarGetLbLocal(vars[i]) < 0.5; ++i ) /* search fixed operator */
      {}
      if( i < nvars )
      {
         SCIPdebugMsg(scip, "constraint <%s>: operator var <%s> fixed to 1.0 -> fix resultant <%s> to 1.0\n",
            SCIPconsGetName(cons), SCIPvarGetName(vars[i]), SCIPvarGetName(resvar));
         SCIP_CALL( SCIPinferBinvarCons(scip, resvar, TRUE, cons, (int)PROPRULE_1, &infeasible, &tightened) );
         if( infeasible )
         {
            /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
            SCIP_CALL( analyzeConflictZero(scip, cons, i) );
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else
         {
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
            if( tightened )
            {
               SCIP_CALL( SCIPresetConsAge(scip, cons) );
               (*nfixedvars)++;
            }
         }

         return SCIP_OKAY;
      }
      else
         consdata->nofixedone = TRUE;
   }
   assert(consdata->nofixedone);

   /* if resultant is fixed to FALSE, all operator variables can be fixed to FALSE (rule (2)) */
   if( SCIPvarGetUbLocal(resvar) < 0.5 )
   {
      for( i = 0; i < nvars && !(*cutoff); ++i )
      {
         SCIPdebugMsg(scip, "constraint <%s>: resultant var <%s> fixed to 0.0 -> fix operator var <%s> to 0.0\n",
            SCIPconsGetName(cons), SCIPvarGetName(resvar), SCIPvarGetName(vars[i]));
         SCIP_CALL( SCIPinferBinvarCons(scip, vars[i], FALSE, cons, (int)PROPRULE_2, &infeasible, &tightened) );
         if( infeasible )
         {
            /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
            SCIP_CALL( analyzeConflictZero(scip, cons, i) );
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( tightened )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            (*nfixedvars)++;
         }
      }

      if( !(*cutoff) )
      {
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }

      return SCIP_OKAY;
   }

   /* rules (3) and (4) can only be applied, if we know all operator variables */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* rules (3) and (4) cannot be applied, if we have at least two unfixed variables left;
    * that means, we only have to watch (i.e. capture events) of two variables, and switch to other variables
    * if these ones get fixed
    */
   watchedvar1 = consdata->watchedvar1;
   watchedvar2 = consdata->watchedvar2;

   /* check, if watched variables are still unfixed */
   if( watchedvar1 != -1 )
   {
      assert(SCIPvarGetLbLocal(vars[watchedvar1]) < 0.5); /* otherwise, rule (1) could be applied */
      if( SCIPvarGetUbLocal(vars[watchedvar1]) < 0.5 )
         watchedvar1 = -1;
   }
   if( watchedvar2 != -1 )
   {
      assert(SCIPvarGetLbLocal(vars[watchedvar2]) < 0.5); /* otherwise, rule (1) could be applied */
      if( SCIPvarGetUbLocal(vars[watchedvar2]) < 0.5 )
         watchedvar2 = -1;
   }

   /* if only one watched variable is still unfixed, make it the first one */
   if( watchedvar1 == -1 )
   {
      watchedvar1 = watchedvar2;
      watchedvar2 = -1;
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if the watched variables are invalid (fixed), find new ones if existing */
   if( watchedvar2 == -1 )
   {
      for( i = 0; i < nvars; ++i )
      {
         assert(SCIPvarGetLbLocal(vars[i]) < 0.5); /* otherwise, rule (1) could be applied */
         if( SCIPvarGetUbLocal(vars[i]) > 0.5 )
         {
            if( watchedvar1 == -1 )
            {
               assert(watchedvar2 == -1);
               watchedvar1 = i;
            }
            else if( watchedvar1 != i )
            {
               watchedvar2 = i;
               break;
            }
         }
      }
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if all variables are fixed to FALSE, the resultant can also be fixed to FALSE (rule (3)) */
   if( watchedvar1 == -1 )
   {
      assert(watchedvar2 == -1);

      SCIPdebugMsg(scip, "constraint <%s>: all operator vars fixed to 0.0 -> fix resultant <%s> to 0.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(resvar));
      SCIP_CALL( SCIPinferBinvarCons(scip, resvar, FALSE, cons, (int)PROPRULE_3, &infeasible, &tightened) );
      if( infeasible )
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflictOne(scip, cons) );
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         if( tightened )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            (*nfixedvars)++;
         }
      }

      return SCIP_OKAY;
   }

   /* if resultant is fixed to TRUE, and only one operator variable is not fixed to FALSE, this operator variable
    * can be fixed to TRUE (rule (4))
    */
   if( SCIPvarGetLbLocal(resvar) > 0.5 && watchedvar2 == -1 )
   {
      assert(watchedvar1 != -1);

      SCIPdebugMsg(scip, "constraint <%s>: resultant <%s> fixed to 1.0, only one unfixed operand -> fix operand <%s> to 1.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(resvar), SCIPvarGetName(vars[watchedvar1]));
      SCIP_CALL( SCIPinferBinvarCons(scip, vars[watchedvar1], TRUE, cons, (int)PROPRULE_4, &infeasible, &tightened) );
      if( infeasible )
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflictOne(scip, cons) );
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         if( tightened )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            (*nfixedvars)++;
         }
      }

      return SCIP_OKAY;
   }

   /* switch to the new watched variables */
   SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, watchedvar1, watchedvar2) );

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** resolves a conflict on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) v_i = TRUE                                   =>  r   = TRUE
 *   (2) r   = FALSE                                  =>  v_i = FALSE for all i
 *   (3) v_i = FALSE for all i                        =>  r   = FALSE
 *   (4) r   = TRUE, v_i = FALSE for all i except j   =>  v_j = TRUE
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   nvars = consdata->nvars;

   switch( proprule )
   {
   case PROPRULE_1:
      /* the resultant was infered to TRUE, because one operand variable was TRUE */
      assert(SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) > 0.5);
      assert(infervar == consdata->resvar);
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, FALSE) > 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
            break;
         }
      }
      assert(i < nvars);
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_2:
      /* the operand variable was infered to FALSE, because the resultant was FALSE */
      assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5);
      assert(SCIPgetVarUbAtIndex(scip, consdata->resvar, bdchgidx, FALSE) < 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_3:
      /* the resultant was infered to FALSE, because all operand variables were FALSE */
      assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5);
      assert(infervar == consdata->resvar);
      for( i = 0; i < nvars; ++i )
      {
         assert(SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, FALSE) < 0.5);
         SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_4:
      /* the operand variable was infered to TRUE, because the resultant was TRUE and all other operands were FALSE */
      assert(SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) > 0.5);
      assert(SCIPgetVarLbAtIndex(scip, consdata->resvar, bdchgidx, FALSE) > 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
      for( i = 0; i < nvars; ++i )
      {
         if( vars[i] != infervar )
         {
            assert(SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, FALSE) < 0.5);
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in or constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** upgrades unmodifiable or constraint into an and constraint on negated variables */
static
SCIP_RETCODE upgradeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   int*                  nupgdconss          /**< pointer to count the number of constraint upgrades */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** negvars;
   SCIP_VAR* negresvar;
   SCIP_CONS* andcons;
   int i;

   assert(nupgdconss != NULL);

   /* we cannot upgrade a modifiable constraint, since we don't know what additional variables to expect */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "upgrading or constraint <%s> into equivalent and constraint on negated variables\n",
      SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get the negated versions of the variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &negvars, consdata->nvars) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPgetNegatedVar(scip, consdata->vars[i], &negvars[i]) );
   }
   SCIP_CALL( SCIPgetNegatedVar(scip, consdata->resvar, &negresvar) );

   /* create and add the and constraint */
   SCIP_CALL( SCIPcreateConsAnd(scip, &andcons, SCIPconsGetName(cons), negresvar, consdata->nvars, negvars,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
         SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, andcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &andcons) );

   /* delete the or constraint */
   SCIP_CALL( SCIPdelCons(scip, cons) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &negvars);

   (*nupgdconss)++;

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyOr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrOr(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeOr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolOr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   /* release and free the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      SCIP_CALL( consdataFreeRows(scip, consdata) );
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransOr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars, sourcedata->vars, sourcedata->resvar) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpOr)
{  /*lint --e{715}*/
   int i;

   *infeasible = FALSE;

   for( i = 0; i < nconss && !(*infeasible); i++ )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i], infeasible) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpOr)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolOr)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpOr)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, FALSE, FALSE, &violated) );
      if( violated )
      {
         SCIP_Bool separated;

         SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
         assert(separated); /* because the solution is integral, the separation always finds a cut */
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxOr)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, FALSE, FALSE, &violated) );
      if( violated )
      {
         SCIP_Bool separated;

         SCIP_CALL( separateCons(scip, conss[i], sol, &separated) );
         assert(separated); /* because the solution is integral, the separation always finds a cut */
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsOr)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, TRUE, FALSE, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckOr)
{  /*lint --e{715}*/
   int i;

   *result = SCIP_FEASIBLE;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss && (*result == SCIP_FEASIBLE || completely); i++ )
   {
      SCIP_Bool violated = FALSE;

      SCIP_CALL( checkCons(scip, conss[i], sol, checklprows, printreason, &violated) );

      if( violated )
         *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int nfixedvars;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   nfixedvars = 0;

   /* propagate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &nfixedvars) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolOr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnupgdconss;
   int c;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnupgdconss = *nupgdconss;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   cutoff = FALSE;
   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->propagated = FALSE;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars) );

      /* remove all variables that are fixed to one */
      SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr) );

      /* transform or constraints into and constraints on the negated variables in order to improve
       * the pairwise constraint presolving possibilities
       */
      SCIP_CALL( upgradeCons(scip, cons, nupgdconss) );

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 1); /* otherwise, propagateCons() has deleted the constraint */

         /* if only one variable is left, the resultant has to be equal to this single variable */
         if( consdata->nvars == 1 )
         {
            SCIPdebugMsg(scip, "or constraint <%s> has only one variable not fixed to 0.0\n", SCIPconsGetName(cons));

            assert(consdata->vars != NULL);
            assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));

            /* aggregate variables: resultant - operand == 0 */
            SCIP_CALL( SCIPaggregateVars(scip, consdata->resvar, consdata->vars[0], 1.0, -1.0, 0.0,
                  &cutoff, &redundant, &aggregated) );
            assert(redundant || SCIPdoNotAggr(scip));

            if( aggregated )
            {
               assert(redundant);
               (*naggrvars)++;
            }

            if( redundant )
            {
               /* delete constraint */
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
         }
         else if( !consdata->impladded )
         {
            int i;

            /* add implications: resultant == 0 -> all operands == 0 */
            for( i = 0; i < consdata->nvars && !cutoff; ++i )
            {
               int nimplbdchgs;

               SCIP_CALL( SCIPaddVarImplication(scip, consdata->resvar, FALSE, consdata->vars[i],
                     SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff, &nimplbdchgs) );
               *nchgbds += nimplbdchgs;
            }
            consdata->impladded = TRUE;
         }

         /* if in r = x or y, the resultant is fixed to one, add implication x = 0 -> y = 1 */
         if( !cutoff && SCIPconsIsActive(cons) && consdata->nvars == 2 && !consdata->opimpladded
            && SCIPvarGetLbGlobal(consdata->resvar) > 0.5 )
         {
            int nimplbdchgs;

            SCIP_CALL( SCIPaddVarImplication(scip, consdata->vars[0], FALSE, consdata->vars[1],
                  SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff, &nimplbdchgs) );
            (*nchgbds) += nimplbdchgs;
            consdata->opimpladded = TRUE;
         }
      }
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nupgdconss > oldnupgdconss )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropOr)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockOr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock resultant variable */
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->resvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );

   /* lock all operand variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintOr)
{  /*lint --e{715}*/

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyOr)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_VAR** vars;
   SCIP_VAR* sourceresvar;
   SCIP_VAR* resvar;
   int nvars;
   int v;

   assert(valid != NULL);
   (*valid) = TRUE;
   resvar = NULL;

   /* get variables that need to be copied */
   sourceresvar = SCIPgetResultantOr(sourcescip, sourcecons); 
   sourcevars = SCIPgetVarsOr(sourcescip, sourcecons);
   nvars = SCIPgetNVarsOr(sourcescip, sourcecons);

   if( nvars == -1 )
      return SCIP_INVALIDCALL;

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   /* map operand variables to active variables of the target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &vars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || vars[v] != NULL);
   }

   /* map resultant to active variable of the target SCIP  */
   if( *valid )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourceresvar, &resvar, varmap, consmap, global, valid) );
      assert(!(*valid) || resvar != NULL);

      if( *valid )
      {
         assert(resvar != NULL);
         SCIP_CALL( SCIPcreateConsOr(scip, cons, SCIPconsGetName(sourcecons), resvar, nvars, vars, 
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      }
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseOr)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR* resvar;
   char* strcopy;
   char* token;
   char* saveptr;
   char* endptr;
   int requiredsize;
   int varssize;
   int nvars;

   SCIPdebugMsg(scip, "parse <%s> as or constraint\n", str);

   /* copy string for truncating it */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &strcopy, str, (int)(strlen(str)+1)));

   /* cutoff "or" form the constraint string */
   token = SCIPstrtok(strcopy, "=", &saveptr ); 

   /* parse variable name */ 
   SCIP_CALL( SCIPparseVarName(scip, token, &resvar, &endptr) );

   if( resvar == NULL )
   {
      SCIPdebugMsg(scip, "resultant variable %s does not exist \n", token);
      *success = FALSE;
   }
   else
   {
      /* cutoff "or" form the constraint string */
      (void) SCIPstrtok(NULL, "(", &saveptr ); 

      /* cutoff ")" form the constraint string */
      token = SCIPstrtok(NULL, ")", &saveptr ); 

      varssize = 100;
      nvars = 0;

      /* allocate buffer array for variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );

      /* parse string */
      SCIP_CALL( SCIPparseVarsList(scip, token, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );

      if( *success )
      {
         /* check if the size of the variable array was great enough */
         if( varssize < requiredsize )
         {
            /* reallocate memory */
            varssize = requiredsize;
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars, varssize) );

            /* parse string again with the correct size of the variable array */
            SCIP_CALL( SCIPparseVarsList(scip, token, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );
         }

         assert(*success);
         assert(varssize >= requiredsize);

         /* create and constraint */
         SCIP_CALL( SCIPcreateConsOr(scip, cons, name, resvar, nvars, vars, 
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      }

      /* free variable buffer */
      SCIPfreeBufferArray(scip, &vars);
   }

   /* free string buffer */
   SCIPfreeBufferArray(scip, &strcopy);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsOr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars + 1 )
      (*success) = FALSE;
   else
   {
      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      vars[consdata->nvars] = consdata->resvar;
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variable (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsOr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars + 1;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecOr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   /* check, if the variable was fixed to one */
   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_LBTIGHTENED )
      consdata->nofixedone = FALSE;

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for or constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;

   /* create event handler for events on variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecOr, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpOr, consEnfopsOr, consCheckOr, consLockOr,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyOr, consCopyOr) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteOr) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolOr) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeOr) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsOr) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsOr) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpOr) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseOr) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolOr, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintOr) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropOr, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropOr) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpOr, consSepasolOr, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY,
         CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransOr) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxOr) );

   return SCIP_OKAY;
}

/** creates and captures an or constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsOr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             resvar,             /**< resultant variable of the operation */
   int                   nvars,              /**< number of operator variables in the constraint */
   SCIP_VAR**            vars,               /**< array with operator variables of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   /* find the or constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("or constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, resvar) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures an or constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicOr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             resvar,             /**< resultant variable of the operation */
   int                   nvars,              /**< number of operator variables in the constraint */
   SCIP_VAR**            vars                /**< array with operator variables of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsOr(scip, cons, name, resvar, nvars, vars, TRUE, TRUE, TRUE, TRUE, TRUE,
         FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** gets number of variables in or constraint */
int SCIPgetNVarsOr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an or constraint\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in or constraint */
SCIP_VAR** SCIPgetVarsOr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an or constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the resultant variable in or constraint */
SCIP_VAR* SCIPgetResultantOr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a or constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->resvar;
}
