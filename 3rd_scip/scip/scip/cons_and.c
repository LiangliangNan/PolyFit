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

/**@file   cons_and.c
 * @brief  Constraint handler for AND-constraints,  \f$r = x_1 \wedge x_2 \wedge \dots  \wedge x_n\f$
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 * This constraint handler deals with AND-constraints. These are constraint of the form:
 *
 * \f[
 *    r = x_1 \wedge x_2 \wedge \dots  \wedge x_n
 * \f]
 *
 * where \f$x_i\f$ is a binary variable for all \f$i\f$. Hence, \f$r\f$ is also of binary type. The variable \f$r\f$ is
 * called resultant and the \f$x\f$'s operators.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_and.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_setppc.h"
#include "scip/cons_pseudoboolean.h"
#include "scip/pub_misc.h"
#include "scip/debug.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "and"
#define CONSHDLR_DESC          "constraint handler for AND-constraints: r = and(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -850100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -850100 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING            (SCIP_PRESOLTIMING_FAST | SCIP_PRESOLTIMING_EXHAUSTIVE)
#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "and"
#define EVENTHDLR_DESC         "bound change event handler for AND-constraints"

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define DEFAULT_LINEARIZE         FALSE /**< should constraint get linearized and removed? */
#define DEFAULT_ENFORCECUTS        TRUE /**< should cuts be separated during LP enforcing? */
#define DEFAULT_AGGRLINEARIZATION FALSE /**< should an aggregated linearization be used? */
#define DEFAULT_UPGRRESULTANT      TRUE /**< should all binary resultant variables be upgraded to implicit binary variables */
#define DEFAULT_DUALPRESOLVING     TRUE /**< should dual presolving be performed? */

#define HASHSIZE_ANDCONS            500 /**< minimal size of hash table in and constraint tables */
#define DEFAULT_PRESOLUSEHASHING   TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define NMINCOMPARISONS          200000 /**< number for minimal pairwise presolving comparisons */
#define MINGAINPERNMINCOMPARISONS 1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */
#define EXPRGRAPHREFORM_PRIORITY 100000 /**< priority of expression graph node reformulation method */

/* @todo maybe use event SCIP_EVENTTYPE_VARUNLOCKED to decide for another dual-presolving run on a constraint */

/*
 * Data structures
 */

/** constraint data for AND-constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables in the AND-constraint */
   SCIP_VAR*             resvar;             /**< resultant variable */
   SCIP_ROW**            rows;               /**< rows for linear relaxation of AND-constraint */
   SCIP_ROW*             aggrrow;            /**< aggregated row for linear relaxation of AND-constraint */
   int                   nvars;              /**< number of variables in AND-constraint */
   int                   varssize;           /**< size of vars array */
   int                   nrows;              /**< number of rows for linear relaxation of AND-constraint */
   int                   watchedvar1;        /**< position of first watched operator variable */
   int                   watchedvar2;        /**< position of second watched operator variable */
   int                   filterpos1;         /**< event filter position of first watched operator variable */
   int                   filterpos2;         /**< event filter position of second watched operator variable */
   unsigned int          propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int          nofixedzero:1;      /**< is none of the operator variables fixed to FALSE? */
   unsigned int          impladded:1;        /**< were the implications of the constraint already added? */
   unsigned int          opimpladded:1;      /**< was the implication for 2 operands with fixed resultant added? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          changed:1;          /**< was constraint changed since last pair preprocessing round? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          checkwhenupgr:1;    /**< if AND-constraint is upgraded to an logicor constraint or the and-
                                              *   constraint is linearized, should the check flag be set to true, even
                                              *   if the AND-constraint has a check flag set to false? */
   unsigned int          notremovablewhenupgr:1;/**< if AND-constraint is upgraded to an logicor constraint or the and-
                                              *   constraint is linearized, should the removable flag be set to false,
                                              *   even if the AND-constraint has a removable flag set to true? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events on watched variables */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
   SCIP_Bool             linearize;          /**< should constraint get linearized and removed? */
   SCIP_Bool             enforcecuts;        /**< should cuts be separated during LP enforcing? */
   SCIP_Bool             aggrlinearization;  /**< should an aggregated linearization be used?  */
   SCIP_Bool             upgrresultant;      /**< upgrade binary resultant variable to an implicit binary variable */
   SCIP_Bool             dualpresolving;     /**< should dual presolving be performed?  */
};


/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_INVALID = 0,                     /**< propagation was applied without a specific propagation rule */
   PROPRULE_1 = 1,                           /**< v_i = FALSE                                  =>  r   = FALSE          */
   PROPRULE_2 = 2,                           /**< r   = TRUE                                   =>  v_i = TRUE for all i */
   PROPRULE_3 = 3,                           /**< v_i = TRUE for all i                         =>  r   = TRUE           */
   PROPRULE_4 = 4                            /**< r   = FALSE, v_i = TRUE for all i except j   =>  v_j = FALSE          */
};
typedef enum Proprule PROPRULE;


/*
 * Local methods
 */

/** installs rounding locks for the given variable in the given AND-constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given AND-constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
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

   /* set event handler for catching bound change events on variables */
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

/** catches events for the watched variable at given position */
static
SCIP_RETCODE consdataCatchWatchedEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< AND-constraint data */
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

   /* catch tightening events for lower bound and relaxed events for upper bounds on watched variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, filterpos) );

   return SCIP_OKAY;
}


/** drops events for the watched variable at given position */
static
SCIP_RETCODE consdataDropWatchedEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< AND-constraint data */
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

   /* drop tightening events for lower bound and relaxed events for upper bounds on watched variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, filterpos) );

   return SCIP_OKAY;
}

/** catches needed events on all variables of constraint, except the special ones for watched variables */
static
SCIP_RETCODE consdataCatchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< AND-constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* catch bound change events for both bounds on resultant variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDCHANGED, 
         eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

   /* catch tightening events for upper bound and relaxed events for lower bounds on operator variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
   }

   return SCIP_OKAY;
}

/** drops events on all variables of constraint, except the special ones for watched variables */
static
SCIP_RETCODE consdataDropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< AND-constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* drop bound change events for both bounds on resultant variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDCHANGED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* drop tightening events for upper bound and relaxed events for lower bounds on operator variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );
   }

   return SCIP_OKAY;
}

/** stores the given variable numbers as watched variables, and updates the event processing */
static
SCIP_RETCODE consdataSwitchWatchedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< AND-constraint data */
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

/** creates constraint data for AND-constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   nvars,              /**< number of variables in the AND-constraint */
   SCIP_VAR**            vars,               /**< variables in AND-constraint */
   SCIP_VAR*             resvar,             /**< resultant variable */
   SCIP_Bool             checkwhenupgr,      /**< should an upgraded constraint be checked despite the fact that this
                                              *   AND-constraint will not be checked
                                              */
   SCIP_Bool             notremovablewhenupgr/**< should an upgraded constraint be  despite the fact that this
                                              *   AND-constraint will not be checked
                                              */
   )
{
   int v;

   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(resvar != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
   (*consdata)->resvar = resvar;
   (*consdata)->rows = NULL;
   (*consdata)->aggrrow = NULL;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->nrows = 0;
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;
   (*consdata)->propagated = FALSE;
   (*consdata)->nofixedzero = FALSE;
   (*consdata)->impladded = FALSE;
   (*consdata)->opimpladded = FALSE;
   (*consdata)->sorted = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->merged = FALSE;
   (*consdata)->checkwhenupgr = checkwhenupgr;
   (*consdata)->notremovablewhenupgr = notremovablewhenupgr;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->resvar, &(*consdata)->resvar) );

      /* catch needed events on variables */
      SCIP_CALL( consdataCatchEvents(scip, *consdata, eventhdlr) );
   }

   assert(SCIPvarIsBinary((*consdata)->resvar));

   /* note: currently, this constraint handler does not handle multiaggregations (e.g. during propagation), hence we forbid
    * multiaggregation from the beginning for the involved variables
    */
   if( SCIPgetStage(scip) <= SCIP_STAGE_EXITPRESOLVE )
   {
      for( v = 0; v < (*consdata)->nvars; ++v )
      {
         assert((*consdata)->vars[v] != NULL);
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars[v]) );
      }
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->resvar) );
   }

   /* capture vars */
   SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->resvar) );
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      assert(SCIPvarIsBinary((*consdata)->vars[v]));
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[v]) );
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
      for( r = 0; r < consdata->nrows; ++r )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->rows, consdata->nrows);

      consdata->nrows = 0;
   }

   if( consdata->aggrrow != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &consdata->aggrrow) );
      consdata->aggrrow = NULL;
   }

   return SCIP_OKAY;
}

/** frees constraint data for AND-constraint */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to the constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int v;

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

   /* release vars */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[v])) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->resvar)) );


   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints AND-constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< AND-constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(consdata != NULL);

   /* print resultant */
   SCIP_CALL( SCIPwriteVarName(scip, file, consdata->resvar, TRUE) );

   /* start the variable list */
   SCIPinfoMessage(scip, file, " == and(");

   /* print variable list */
   SCIP_CALL( SCIPwriteVarsList(scip, file, consdata->vars, consdata->nvars, TRUE, ',') );

   /* close the variable list */
   SCIPinfoMessage(scip, file, ")");

   return SCIP_OKAY;
}

/** adds coefficient to AND-constraint */
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
   consdata->sorted = (consdata->nvars == 1);
   consdata->changed = TRUE;
   consdata->merged = FALSE;

   /* capture variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   /* if we are in transformed problem, catch the variable's events */
   if( transformed )
   {
      /* catch bound change events of variable */
      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
   }

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockRounding(scip, cons, var) );

   /**@todo update LP rows */
   if( consdata->rows != NULL )
   {
      SCIPerrorMessage("cannot add coefficients to AND-constraint after LP relaxation was created\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from AND-constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint */
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
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
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

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &(consdata->vars[pos])) );

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->nvars--;

   /* if the last variable (that moved) was watched, update the watched position */
   if( consdata->watchedvar1 == consdata->nvars )
      consdata->watchedvar1 = pos;
   if( consdata->watchedvar2 == consdata->nvars )
      consdata->watchedvar2 = pos;

   consdata->propagated = FALSE;
   consdata->sorted = FALSE;
   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** sorts AND-constraint's variables by non-decreasing variable index */
static
void consdataSort(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->sorted )
   {
      if( consdata->nvars <= 1 )
	 consdata->sorted = TRUE;
      else
      {
	 SCIP_VAR* var1 = NULL;
	 SCIP_VAR* var2 = NULL;

	 /* remember watch variables */
	 if( consdata->watchedvar1 != -1 )
	 {
	    var1 = consdata->vars[consdata->watchedvar1];
	    assert(var1 != NULL);
	    consdata->watchedvar1 = -1;
	    if( consdata->watchedvar2 != -1 )
	    {
	       var2 = consdata->vars[consdata->watchedvar2];
	       assert(var2 != NULL);
	       consdata->watchedvar2 = -1;
	    }
	 }
	 assert(consdata->watchedvar1 == -1);
	 assert(consdata->watchedvar2 == -1);
	 assert(var1 != NULL || var2 == NULL);

	 /* sort variables after index */
	 SCIPsortPtr((void**)consdata->vars, SCIPvarComp, consdata->nvars);
	 consdata->sorted = TRUE;

	 /* correct watched variables */
	 if( var1 != NULL )
	 {
	    int pos;
#ifndef NDEBUG
	    SCIP_Bool found;

	    found = SCIPsortedvecFindPtr((void**)consdata->vars, SCIPvarComp, (void*)var1, consdata->nvars, &pos);
	    assert(found);
#else
	    SCIPsortedvecFindPtr((void**)consdata->vars, SCIPvarComp, (void*)var1, consdata->nvars, &pos);
#endif
	    assert(pos >= 0 && pos < consdata->nvars);
	    consdata->watchedvar1 = pos;

	    if( var2 != NULL )
	    {
#ifndef NDEBUG
	       found = SCIPsortedvecFindPtr((void**)consdata->vars, SCIPvarComp, (void*)var2, consdata->nvars, &pos);
	       assert(found);
#else
	       SCIPsortedvecFindPtr((void**)consdata->vars, SCIPvarComp, (void*)var2, consdata->nvars, &pos);
#endif
	       assert(pos >= 0 && pos < consdata->nvars);
	       consdata->watchedvar2 = pos;
	    }
	 }
      }
   }

#ifdef SCIP_DEBUG
   /* check sorting */
   {
      int v;

      for( v = 0; v < consdata->nvars; ++v )
      {
         assert(v == consdata->nvars-1 || SCIPvarCompare(consdata->vars[v], consdata->vars[v+1]) <= 0);
      }
   }
#endif
}

/** deletes all one-fixed variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int*                  nchgcoefs           /**< pointer to add up the number of changed coefficients */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);
   assert(nchgcoefs != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs)++;
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

#ifdef SCIP_DISABLED_CODE /* does not work with pseudoboolean constraint handler, need to be fixed */
   /* check, if the resultant should be replaced with the active representative */
   if( !SCIPvarIsActive(consdata->resvar) )
   {
      SCIP_VAR* repvar;
      SCIP_Bool negated;

      /* get binary representative of variable */
      SCIP_CALL( SCIPgetBinvarRepresentative(scip, consdata->resvar, &repvar, &negated) );
      assert(SCIPvarIsBinary(repvar));

      /* check, if the variable should be replaced with the representative */
      if( repvar != consdata->resvar )
      {
	 if( SCIPconsIsTransformed(cons) )
	 {
	    /* drop bound change events of old resultant */
	    SCIP_CALL( SCIPdropVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDCHANGED,
		  eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

	    /* catch bound change events of new resultant */
	    SCIP_CALL( SCIPcatchVarEvent(scip, repvar, SCIP_EVENTTYPE_BOUNDCHANGED,
		  eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
	 }

	 /* release old resultant */
	 SCIP_CALL( SCIPreleaseVar(scip, &(consdata->resvar)) );

	 /* capture new resultant */
	 SCIP_CALL( SCIPcaptureVar(scip, repvar) );

	 consdata->resvar = repvar;
	 consdata->changed = TRUE;
      }
   }
#endif

   SCIPdebugMsg(scip, "after fixings: ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL)) );
   SCIPdebugMsgPrint(scip, "\n");

   return SCIP_OKAY;
}

/** creates a linearization of the AND-constraint */
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
   consdata->nrows = nvars + 1;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->rows, consdata->nrows) );

   /* creates LP rows corresponding to AND-constraint:
    *   - one additional row:             resvar - v1 - ... - vn >= 1-n
    *   - for each operator variable vi:  resvar - vi            <= 0
    */

   /* create additional row */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_add", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[0], SCIPconsGetHdlr(cons), rowname, -consdata->nvars + 1.0, SCIPinfinity(scip),
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[0], consdata->resvar, 1.0) );
   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[0], nvars, consdata->vars, -1.0) );

   /* create operator rows */
   for( i = 0; i < nvars; ++i )
   {
      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), i);
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[i+1], SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), 0.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[i+1], consdata->resvar, 1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[i+1], consdata->vars[i], -1.0) );
   }

   return SCIP_OKAY;
}

/** adds linear relaxation of AND-constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   SCIP_CONSDATA* consdata;

   /* in the root LP we only add the weaker relaxation which consists of two rows:
    *   - one additional row:             resvar - v1 - ... - vn >= 1-n
    *   - aggregated row:               n*resvar - v1 - ... - vn <= 0.0
    *
    * during separation we separate the stronger relaxation which consists of n+1 row:
    *   - one additional row:             resvar - v1 - ... - vn >= 1-n
    *   - for each operator variable vi:  resvar - vi            <= 0.0
    */

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* create the aggregated row */
   if( consdata->aggrrow == NULL )
   {
      char rowname[SCIP_MAXSTRLEN];

      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_operators", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->aggrrow, SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), 0.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->aggrrow, consdata->resvar, (SCIP_Real) consdata->nvars) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->aggrrow, consdata->nvars, consdata->vars, -1.0) );
   }

   /* insert aggregated LP row as cut */
   if( !SCIProwIsInLP(consdata->aggrrow) )
   {
      SCIP_CALL( SCIPaddRow(scip, consdata->aggrrow, FALSE, infeasible) );
   }

   if( !(*infeasible) )
   {
      if( consdata->rows == NULL )
      {
         /* create the n+1 row relaxation */
         SCIP_CALL( createRelaxation(scip, cons) );
      }

      assert(consdata->rows != NULL);

      /* add additional row */
      if( !SCIProwIsInLP(consdata->rows[0]) )
      {
         SCIP_CALL( SCIPaddRow(scip, consdata->rows[0], FALSE, infeasible) );
      }
   }

   return SCIP_OKAY;
}

/** checks AND-constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
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

   /* check whether we can skip this feasibility check, because all rows are in the LP and do not have to be checked */
   mustcheck = checklprows;
   mustcheck = mustcheck || (consdata->rows == NULL);
   if( !mustcheck )
   {
      assert(consdata->rows != NULL);

      for( r = 0; r < consdata->nrows; ++r )
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
      SCIP_Real viol;
      SCIP_Real absviol;
      SCIP_Real relviol;
      int i;

      /* increase age of constraint; age is reset to zero, if a violation was found only in case we are in
       * enforcement
       */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }

      absviol = 0.0;
      relviol = 0.0;

      /* check, if all operator variables are TRUE */
      for( i = 0; i < consdata->nvars; ++i )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[i]);

         viol = REALABS(1 - solval);
         if( absviol < viol )
         {
            absviol = viol;
            relviol = SCIPrelDiff(solval, 1.0);
         }

        /* @todo If "upgraded resultants to varstatus implicit" is fully allowed, than the following assert does not hold
         *       anymore, therefor we need to stop the check and return with the status not violated, because the
         *       integrality condition of this violated operand needs to be enforced by another constraint.
         *
         *       The above should be asserted by marking the constraint handler, for which the result needs to be
         *       SCIP_SEPARATED if the origin was the CONSENFOPS or the CONSENFOLP callback or SCIP_INFEASIBLE if the
         *       origin was CONSCHECK callback.
         *
         */
         assert(SCIPisFeasIntegral(scip, solval));
         if( solval < 0.5 )
            break;
      }

      /* if all operator variables are TRUE, the resultant has to be TRUE, otherwise, the resultant has to be FALSE;
       * in case of an implicit integer resultant variable, we need to ensure the integrality of the solution value
       */
      solval = SCIPgetSolVal(scip, sol, consdata->resvar);
      assert(SCIPvarGetType(consdata->resvar) == SCIP_VARTYPE_IMPLINT || SCIPisFeasIntegral(scip, solval));

      if( !SCIPisFeasIntegral(scip, solval) || (i == consdata->nvars) != (solval > 0.5) )
      {
         *violated = TRUE;
         absviol = 1.0;
         relviol = 1.0;

         /* only reset constraint age if we are in enforcement */
         if( sol == NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }

         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");
            SCIPinfoMessage(scip, NULL, "violation:");
            if( !SCIPisFeasIntegral(scip, solval) )
            {
               SCIPinfoMessage(scip, NULL, " resultant variable <%s> has fractional solution value %" SCIP_REAL_FORMAT "\n",
                     SCIPvarGetName(consdata->resvar), solval);
            }
            else if( i == consdata->nvars )
            {
               SCIPinfoMessage(scip, NULL, " all operands are TRUE and resultant <%s> = FALSE\n",
                  SCIPvarGetName(consdata->resvar));
            }
            else
            {
               SCIPinfoMessage(scip, NULL, " operand <%s> = FALSE and resultant <%s> = TRUE\n",
                  SCIPvarGetName(consdata->vars[i]), SCIPvarGetName(consdata->resvar));
            }
         }
      }
      if( sol != NULL )
         SCIPupdateSolConsViolation(scip, sol, absviol, relviol);
   }

   return SCIP_OKAY;
}

/** separates given primal solution */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            separated,          /**< pointer to store whether a cut was found */
   SCIP_Bool*            cutoff              /**< whether a cutoff has been detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real feasibility;
   int r;

   assert(separated != NULL);
   assert(cutoff != NULL);

   *separated = FALSE;
   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* create all necessary rows for the linear relaxation */
   if( consdata->rows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows != NULL);

   /* test all rows for feasibility and add infeasible rows */
   for( r = 0; r < consdata->nrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
      {
         feasibility = SCIPgetRowSolFeasibility(scip, consdata->rows[r], sol);
         if( SCIPisFeasNegative(scip, feasibility) )
         {
            SCIP_CALL( SCIPaddRow(scip, consdata->rows[r], FALSE, cutoff) );
            if ( *cutoff )
               return SCIP_OKAY;
            *separated = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflicting TRUE assignment to resultant of given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflictOne(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint that detected the conflict */
   int                   falsepos            /**< position of operand that is fixed to FALSE */
   )
{
   SCIP_CONSDATA* consdata;

   /* conflict analysis can only be applied in solving stage and if it turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetLbLocal(consdata->resvar) > 0.5);
   assert(0 <= falsepos && falsepos < consdata->nvars);
   assert(SCIPvarGetUbLocal(consdata->vars[falsepos]) < 0.5);

   /* initialize conflict analysis, and add resultant and single operand variable to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[falsepos]) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** analyzes conflicting FALSE assignment to resultant of given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflictZero(
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
   assert(SCIPvarGetUbLocal(consdata->resvar) < 0.5);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(SCIPvarGetLbLocal(consdata->vars[v]) > 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** tries to fix the given resultant to zero */
static
SCIP_RETCODE consdataFixResultantZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint to be processed */
   SCIP_VAR*             resvar,             /**< resultant variable to fix to zero */
   int                   pos,                /**< position of operand that is fixed to FALSE */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   SCIPdebugMsg(scip, "constraint <%s>: operator %d fixed to 0.0 -> fix resultant <%s> to 0.0\n",
      SCIPconsGetName(cons), pos, SCIPvarGetName(resvar));

   SCIP_CALL( SCIPinferBinvarCons(scip, resvar, FALSE, cons, (int)PROPRULE_1, &infeasible, &tightened) );

   if( infeasible )
   {
      /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
      SCIP_CALL( analyzeConflictOne(scip, cons, pos) );
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      (*cutoff) = TRUE;
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

/** fix all operands to one */
static
SCIP_RETCODE consdataFixOperandsOne(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint to be processed */
   SCIP_VAR**            vars,               /**< array of operands */
   int                   nvars,              /**< number of operands */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   int v;

   for( v = 0; v < nvars && !(*cutoff); ++v )
   {
      SCIPdebugMsg(scip, "constraint <%s>: resultant fixed to 1.0 -> fix operator var <%s> to 1.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(vars[v]));

      SCIP_CALL( SCIPinferBinvarCons(scip, vars[v], TRUE, cons, (int)PROPRULE_2, &infeasible, &tightened) );

      if( infeasible )
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflictOne(scip, cons, v) );
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         (*cutoff) = TRUE;
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

/** linearize AND-constraint due to a globally to zero fixed resultant; that is, creates, adds, and releases a logicor
 *  constraint and remove the AND-constraint globally.
 *
 *  Since the resultant is fixed to zero the AND-constraint collapses to linear constraint of the form:
 *
 *  - \f$\sum_{i=0}^{n-1} v_i \leq n-1\f$
 *
 *  This can be transformed into a logicor constraint of the form
 *
 *  - \f$\sum_{i=0}^{n-1} ~v_i \geq 1\f$
 */
static
SCIP_RETCODE consdataLinearize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint to linearize */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  nupgdconss          /**< pointer to add up the number of upgraded constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_CONS* lincons;
   SCIP_Bool conscreated;
   int nvars;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(!(*cutoff));
   assert(SCIPvarGetUbGlobal(consdata->resvar) < 0.5);

   nvars = consdata->nvars;
   conscreated = FALSE;

   /* allocate memory for variables for updated constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   /* if we only have two variables, we prefer a set packing constraint instead of a logicor constraint */
   if( nvars == 2 && !SCIPconsIsModifiable(cons) )
   {
      SCIP_Bool* negated;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      /* get active representation */
      SCIP_CALL( SCIPallocBufferArray(scip, &negated, nvars) );
      SCIP_CALL( SCIPgetBinvarRepresentatives(scip, nvars, consdata->vars, vars, negated) );
      SCIPfreeBufferArray(scip, &negated);

      /* if one of the two operators is globally fixed to one it follows that the other has to be zero */
      if( SCIPvarGetLbGlobal(vars[0]) > 0.5 )
      {
         SCIP_CALL( SCIPfixVar(scip, vars[1], 0.0, &infeasible, &tightened) );

         if( infeasible )
            *cutoff = TRUE;
         else if( tightened )
            ++(*nfixedvars);
      }
      else if( SCIPvarGetLbGlobal(vars[1]) > 0.5 )
      {
         SCIP_CALL( SCIPfixVar(scip, vars[0], 0.0, &infeasible, &tightened) );

         if( infeasible )
            *cutoff = TRUE;
         else if( tightened )
            ++(*nfixedvars);
      }
      else if( SCIPvarGetUbGlobal(vars[0]) > 0.5 && SCIPvarGetUbGlobal(vars[1]) > 0.5 )
      {
         /* create, add, and release the setppc constraint */
         SCIP_CALL( SCIPcreateConsSetpack(scip, &lincons, SCIPconsGetName(cons), nvars, vars,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               consdata->checkwhenupgr || SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
               !(consdata->notremovablewhenupgr) && SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

         conscreated = TRUE;
      }
   }
   else
   {
      int v;

      /* collect negated variables */
      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, consdata->vars[v], &vars[v]) );
      }

      /* create, add, and release the logicor constraint */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &lincons, SCIPconsGetName(cons), nvars, vars,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            consdata->checkwhenupgr || SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            !(consdata->notremovablewhenupgr) && SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

      conscreated = TRUE;
   }

   if( conscreated )
   {
      /* add and release new constraint */
      SCIPdebugPrintCons(scip, lincons, NULL); /*lint !e644*/
      SCIP_CALL( SCIPaddCons(scip, lincons) ); /*lint !e644*/
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) ); /*lint !e644*/

      ++(*nupgdconss);
   }

   /* remove the AND-constraint globally */
   SCIP_CALL( SCIPdelCons(scip, cons) );

   /* delete temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** the resultant is fixed to zero; in case all except one operator are fixed to TRUE the last operator has to fixed to FALSE */
/** @note consdata->watchedvars might not be the same to the watchedvar parameters, because the update was not yet done */
static
SCIP_RETCODE analyzeZeroResultant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint to be processed */
   int                   watchedvar1,        /**< maybe last unfixed variable position */
   int                   watchedvar2,        /**< second watched position */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetUbLocal(consdata->resvar) < 0.5);

   if( watchedvar2 == -1 )
   {
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      assert(watchedvar1 != -1);

#ifndef NDEBUG
      /* check that all variables regardless of wathcedvar1 are fixed to 1 */
      {
	 int v;

	 for( v = consdata->nvars - 1; v >= 0; --v )
	    if( v != watchedvar1 )
	       assert(SCIPvarGetLbLocal(consdata->vars[v]) > 0.5);
      }
#endif

      SCIPdebugMsg(scip, "constraint <%s>: resultant <%s> fixed to 0.0, only one unfixed operand -> fix operand <%s> to 0.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(consdata->resvar), SCIPvarGetName(consdata->vars[watchedvar1]));

      SCIP_CALL( SCIPinferBinvarCons(scip, consdata->vars[watchedvar1], FALSE, cons, (int)PROPRULE_4, &infeasible, &tightened) );

      if( infeasible )
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflictZero(scip, cons) );
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
   }

   return SCIP_OKAY;
}

/** replaces multiple occurrences of variables */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   unsigned char**       entries,            /**< array to store whether two positions in constraints represent the same variable */
   int*                  nentries,           /**< pointer for array size, if array will be to small it's corrected */
   int*                  nfixedvars,         /**< pointer to store number of fixed variables */
   int*                  nchgcoefs,          /**< pointer to store number of changed coefficients */
   int*                  ndelconss           /**< pointer to store number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_VAR* probvar;
   int probidx;
   int nvars;
   int v;
#ifndef NDEBUG
   int nbinvars;
   int nintvars;
   int nimplvars;
#endif

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);
   assert(*entries != NULL);
   assert(nentries != NULL);
   assert(nfixedvars != NULL);
   assert(nchgcoefs != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged )
      return SCIP_OKAY;

   /* nothing to merge */
   if( consdata->nvars <= 1 )
   {
      consdata->merged = TRUE;
      return SCIP_OKAY;
   }

   vars = consdata->vars;
   nvars = consdata->nvars;

   assert(vars != NULL);
   assert(nvars >= 2);

#ifndef NDEBUG
   nbinvars = SCIPgetNBinVars(scip);
   nintvars = SCIPgetNIntVars(scip);
   nimplvars = SCIPgetNImplVars(scip);
   assert(*nentries >= nbinvars + nintvars + nimplvars);
#endif

   /* initialize entries array */
   for( v = nvars - 1; v >= 0; --v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarIsActive(var) || (SCIPvarIsNegated(var) && SCIPvarIsActive(SCIPvarGetNegatedVar(var))));

      probvar = (SCIPvarIsActive(var) ? var : SCIPvarGetNegatedVar(var));
      assert(probvar != NULL);

      probidx = SCIPvarGetProbindex(probvar);
      assert(0 <= probidx);

      /* check variable type, either pure binary or an integer/implicit integer variable with 0/1 bounds */
      assert((probidx < nbinvars && SCIPvarGetType(probvar) == SCIP_VARTYPE_BINARY)
	 || (SCIPvarIsBinary(probvar) &&
            ((probidx >= nbinvars && probidx < nbinvars + nintvars && SCIPvarGetType(probvar) == SCIP_VARTYPE_INTEGER) ||
               (probidx >= nbinvars + nintvars && probidx < nbinvars + nintvars + nimplvars &&
                  SCIPvarGetType(probvar) == SCIP_VARTYPE_IMPLINT))));

      /* var is not active yet */
      (*entries)[probidx] = 0;
   }

   /* search for multiple variables; scan from back to front because deletion doesn't affect the order of the front
    * variables
    * @note don't reorder variables because we would loose the watched variables and filter position inforamtion
    */
   for( v = nvars - 1; v >= 0; --v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarIsActive(var) || (SCIPvarIsNegated(var) && SCIPvarIsActive(SCIPvarGetNegatedVar(var))));

      probvar = (SCIPvarIsActive(var) ? var : SCIPvarGetNegatedVar(var));
      assert(probvar != NULL);

      probidx = SCIPvarGetProbindex(probvar);
      assert(0 <= probidx && probidx < *nentries);

      /* if var occurs first time in constraint init entries array */
      if( (*entries)[probidx] == 0 )
      {
	 (*entries)[probidx] = (SCIPvarIsActive(var) ? 1 : 2);
      }
      /* if var occurs second time in constraint, first time it was not negated */
      else if( ((*entries)[probidx] == 1 && SCIPvarIsActive(var)) || ((*entries)[probidx] == 2 && !SCIPvarIsActive(var)) )
      {
         /* delete the multiple variable */
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         ++(*nchgcoefs);
      }
      else
      {
	 SCIP_Bool infeasible;
	 SCIP_Bool fixed;

	 assert(((*entries)[probidx] == 1 && !SCIPvarIsActive(var)) || ((*entries)[probidx] == 2 && SCIPvarIsActive(var)));

	 SCIPdebugMsg(scip, "AND-constraint <%s> is redundant: variable <%s> and its negation are present -> fix resultant <%s> = 0\n",
	    SCIPconsGetName(cons), SCIPvarGetName(var), SCIPvarGetName(consdata->resvar));

	 /* negation of the variable is already present in the constraint: fix resultant to zero */
#ifndef NDEBUG
	 {
	    int i;
	    for( i = consdata->nvars - 1; i > v && var != SCIPvarGetNegatedVar(vars[i]); --i )
	    {}
	    assert(i > v);
	 }
#endif

	 SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 0.0, &infeasible, &fixed) );
	 assert(!infeasible);
	 if( fixed )
	    ++(*nfixedvars);

	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 break;
      }
   }

   consdata->merged = TRUE;

   return SCIP_OKAY;
}

/** propagates constraint with the following rules:
 *   (1) v_i = FALSE                                  =>  r   = FALSE
 *   (2) r   = TRUE                                   =>  v_i = TRUE for all i
 *   (3) v_i = TRUE for all i                         =>  r   = TRUE
 *   (4) r   = FALSE, v_i = TRUE for all i except j   =>  v_j = FALSE
 *
 *  additional if the resultant is fixed to zero during presolving or in the root node (globally), then the
 *  AND-constraint is collapsed to a linear (logicor) constraint of the form
 *  -> sum_{i=0}^{n-1} ~v_i >= 1
 */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< AND-constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  nupgdconss          /**< pointer to add up the number of upgraded constraints */
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

   /* don't process the constraint, if none of the operator variables was fixed to FALSE, and if the watched variables
    * and the resultant weren't fixed to any value since last propagation call
    */
   if( consdata->propagated )
   {
      assert(consdata->nofixedzero);
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(resvar), 0.0));
      return SCIP_OKAY;
   }

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   /* if one of the operator variables was fixed to FALSE, the resultant can be fixed to FALSE (rule (1)) */
   if( !consdata->nofixedzero )
   {
      for( i = 0; i < nvars && SCIPvarGetUbLocal(vars[i]) > 0.5; ++i ) /* search for operator fixed to zero */
      {}
      if( i < nvars )
      {
         /* fix resultant to zero */
         SCIP_CALL( consdataFixResultantZero(scip, cons, resvar, i, cutoff, nfixedvars) );
      }
      else
         consdata->nofixedzero = TRUE;
   }

   /* check if resultant variables is globally fixed to zero */
   if( !SCIPinProbing(scip) && SCIPvarGetUbGlobal(resvar) < 0.5 )
   {
      SCIP_CALL( consdataLinearize(scip, cons, cutoff, nfixedvars, nupgdconss) );

      if( *cutoff && SCIPgetDepth(scip) > 0 )
      {
         /* we are done with solving since a global bound change was infeasible */
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
      }

      return SCIP_OKAY;
   }

   /* if the resultant and at least one operand are locally fixed to zero, the constraint is locally redundant */
   if( SCIPvarGetUbLocal(resvar) < 0.5 && !consdata->nofixedzero )
   {
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      return SCIP_OKAY;
   }

   /* if resultant is fixed to TRUE, all operator variables can be fixed to TRUE (rule (2)) */
   if( SCIPvarGetLbLocal(resvar) > 0.5 )
   {
      /* fix operands to one */
      SCIP_CALL( consdataFixOperandsOne(scip, cons, vars, nvars, cutoff, nfixedvars) );

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
      assert(SCIPvarGetUbLocal(vars[watchedvar1]) > 0.5); /* otherwise, rule (1) could be applied */
      if( SCIPvarGetLbLocal(vars[watchedvar1]) > 0.5 )
         watchedvar1 = -1;
   }
   if( watchedvar2 != -1 )
   {
      assert(SCIPvarGetUbLocal(vars[watchedvar2]) > 0.5); /* otherwise, rule (1) could be applied */
      if( SCIPvarGetLbLocal(vars[watchedvar2]) > 0.5 )
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
         assert(SCIPvarGetUbLocal(vars[i]) > 0.5); /* otherwise, rule (1) could be applied */
         if( SCIPvarGetLbLocal(vars[i]) < 0.5 )
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

   /* if all variables are fixed to TRUE, the resultant can also be fixed to TRUE (rule (3)) */
   if( watchedvar1 == -1 )
   {
      assert(watchedvar2 == -1);

      SCIPdebugMsg(scip, "constraint <%s>: all operator vars fixed to 1.0 -> fix resultant <%s> to 1.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(resvar));
      SCIP_CALL( SCIPinferBinvarCons(scip, resvar, TRUE, cons, (int)PROPRULE_3, &infeasible, &tightened) );

      if( infeasible )
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflictZero(scip, cons) );
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

   /* if resultant is fixed to FALSE, and only one operator variable is not fixed to TRUE, this operator variable
    * can be fixed to FALSE (rule (4))
    */
   if( watchedvar2 == -1 && SCIPvarGetUbLocal(resvar) < 0.5 )
   {
      assert(watchedvar1 != -1);

      SCIP_CALL( analyzeZeroResultant(scip, cons, watchedvar1, watchedvar2, cutoff, nfixedvars) );

      return SCIP_OKAY;
   }

   /* switch to the new watched variables */
   SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, watchedvar1, watchedvar2) );

   /* mark the constraint propagated if we have an unfixed resultant or are not in probing, it is necessary that a fixed
    * resulting in probing mode does not lead to a propagated constraint, because the constraint upgrade needs to be performed
    */
   consdata->propagated = (!SCIPinProbing(scip) || (SCIPvarGetLbLocal(consdata->resvar) < 0.5 && SCIPvarGetUbLocal(consdata->resvar) > 0.5));

   return SCIP_OKAY;
}

/** resolves a conflict on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) v_i = FALSE                                  =>  r   = FALSE
 *   (2) r   = TRUE                                   =>  v_i = TRUE for all i
 *   (3) v_i = TRUE for all i                         =>  r   = TRUE
 *   (4) r   = FALSE, v_i = TRUE for all i except j   =>  v_j = FALSE
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
      /* the resultant was infered to FALSE, because one operand variable was FALSE */
      assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5);
      assert(infervar == consdata->resvar);
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, FALSE) < 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
            break;
         }
      }
      assert(i < nvars);
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_2:
      /* the operand variable was infered to TRUE, because the resultant was TRUE */
      assert(SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) > 0.5);
      assert(SCIPgetVarLbAtIndex(scip, consdata->resvar, bdchgidx, FALSE) > 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_3:
      /* the resultant was infered to TRUE, because all operand variables were TRUE */
      assert(SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) > 0.5);
      assert(infervar == consdata->resvar);
      for( i = 0; i < nvars; ++i )
      {
         assert(SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, FALSE) > 0.5);
         SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_4:
      /* the operand variable was infered to FALSE, because the resultant was FALSE and all other operands were TRUE */
      assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5);
      assert(SCIPgetVarUbAtIndex(scip, consdata->resvar, bdchgidx, FALSE) < 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
      for( i = 0; i < nvars; ++i )
      {
         if( vars[i] != infervar )
         {
            assert(SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, FALSE) > 0.5);
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in AND-constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** perform dual presolving on AND-constraints */
static
SCIP_RETCODE dualPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< AND-constraints to perform dual presolving on */
   int                   nconss,             /**< number of AND-constraints */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   unsigned char**       entries,            /**< array to store whether two positions in constraints represent the same variable */
   int*                  nentries,           /**< pointer for array size, if array will be to small it's corrected */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  naggrvars,          /**< pointer to add up the number of aggregated variables */
   int*                  nchgcoefs,          /**< pointer to add up the number of changed coefficients */
   int*                  ndelconss,          /**< pointer to add up the number of deleted constraints */
   int*                  nupgdconss,         /**< pointer to add up the number of upgraded constraints */
   int*                  naddconss           /**< pointer to add up the number of added constraints */
   )
{
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** impoperands;
   SCIP_VAR** vars;
   SCIP_VAR* resvar;
   SCIP_VAR* var;
   int nimpoperands;
   int nvars;
   int size;
   int v;
   int c;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(eventhdlr != NULL);
   assert(*entries != NULL);
   assert(nentries != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgcoefs != NULL);
   assert(ndelconss != NULL);
   assert(nupgdconss != NULL);
   assert(naddconss != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   assert(conss != NULL);

   size = 2 * (SCIPgetNBinVars(scip) + SCIPgetNImplVars(scip));

   SCIP_CALL( SCIPallocBufferArray(scip, &impoperands, size) );

   for( c = nconss - 1; c >= 0 && !(*cutoff); --c )
   {
      cons = conss[c];
      assert(cons != NULL);

      if( !SCIPconsIsActive(cons) || !SCIPconsIsChecked(cons) || SCIPconsIsModifiable(cons) )
	 continue;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, eventhdlr, cutoff, nfixedvars, nupgdconss) );

      if( !SCIPconsIsActive(cons) )
	 continue;

      if( *cutoff )
	 break;

      SCIP_CALL( applyFixings(scip, cons, eventhdlr, nchgcoefs) );

      /* merge multiple occurances of variables or variables with their negated variables */
      SCIP_CALL( mergeMultiples(scip, cons, eventhdlr, entries, nentries, nfixedvars, nchgcoefs, ndelconss) );

      if( !SCIPconsIsActive(cons) )
	 continue;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      vars = consdata->vars;
      nvars = consdata->nvars;
      assert(vars != NULL || nvars == 0);

      if( nvars == 0 )
	 continue;

      assert(vars != NULL);

      resvar = consdata->resvar;
      assert(resvar != NULL);
      /* a fixed resultant needs to be removed, otherwise we might fix operands to a wrong value later on */
      assert(SCIPvarGetLbGlobal(resvar) < 0.5 && SCIPvarGetUbGlobal(resvar) > 0.5);
      assert(SCIPvarGetNLocksUp(resvar) >= 1 && SCIPvarGetNLocksDown(resvar) >= 1);

      if( SCIPvarGetNLocksUp(resvar) == 1 && SCIPvarGetNLocksDown(resvar) == 1 )
      {
	 SCIP_Real resobj;
	 SCIP_Real obj;
	 SCIP_Real posobjsum = 0;
	 SCIP_Real maxobj = -SCIPinfinity(scip);
	 int maxpos = -1;
	 int oldnfixedvars = *nfixedvars;
	 int oldnaggrvars = *naggrvars;

	 nimpoperands = 0;

	 /* collect important operands */
	 for( v = nvars - 1; v >= 0; --v )
	 {
	    var = vars[v];
	    assert(var != NULL);
	    assert(SCIPvarGetNLocksUp(var) >= 1 && SCIPvarGetNLocksDown(var) >= 1);

	    if( SCIPvarGetNLocksUp(var) == 1 && SCIPvarGetNLocksDown(var) == 1 )
	    {
	       impoperands[nimpoperands] = var;
	       ++nimpoperands;

	       /* get aggregated objective value of active variable */
	       SCIP_CALL( SCIPvarGetAggregatedObj(var, &obj) );

	       /* add up all positive objective values of operands which have exactly one lock in both directions */
	       if( obj > 0 )
		  posobjsum += obj;

	       /* memorize maximal objective value of operands and its position */
	       if( obj > maxobj )
	       {
		  maxpos = nimpoperands - 1;
		  maxobj = obj;
	       }
	    }
	 }
	 assert(nimpoperands >= 0 && nimpoperands <= nvars);

	 /* no dual fixable variables found */
	 if( nimpoperands == 0 )
	    continue;

	 /* get aggregated objective value of active variable */
	 SCIP_CALL( SCIPvarGetAggregatedObj(resvar, &resobj) );

	 /* resultant contributes to the objective with a negative value */
	 if( SCIPisLE(scip, resobj, 0.0) )
	 {
	    SCIP_Bool poscontissmall = SCIPisLE(scip, posobjsum, REALABS(resobj));

	    /* if all variables are only locked by this constraint and the resultants contribution more then compensates
	     * the positive contribution, we can fix all variables to 1
	     */
	    if( nimpoperands == nvars && poscontissmall )
	    {
	       SCIPdebugMsg(scip, "dual-fixing all variables in constraint <%s> to 1\n", SCIPconsGetName(cons));

	       SCIP_CALL( SCIPfixVar(scip, resvar, 1.0, &infeasible, &fixed) );

	       *cutoff = *cutoff || infeasible;
	       if( fixed )
		  ++(*nfixedvars);


	       for( v = nvars - 1; v >= 0 && !(*cutoff); --v )
	       {
		  SCIP_CALL( SCIPfixVar(scip, vars[v], 1.0, &infeasible, &fixed) );

		  *cutoff = *cutoff || infeasible;
		  if( fixed )
		     ++(*nfixedvars);
	       }

	       SCIPdebugMsg(scip, "deleting constraint <%s> because all variables are fixed to one\n", SCIPconsGetName(cons));

	       SCIP_CALL( SCIPdelCons(scip, cons) );
	       ++(*ndelconss);
	    }
	    else
	    {
               SCIP_Bool aggregationperformed = FALSE;
               SCIP_Bool zerofix = FALSE;

	       assert(nimpoperands > 0);

	       SCIPdebugMsg(scip, "dual-fixing all variables in constraint <%s> with positive contribution (when together exceeding the negative contribution of the resultant) to 0 and with negative contribution to 1\n", SCIPconsGetName(cons));

	       for( v = nimpoperands - 1; v >= 0 && !(*cutoff); --v )
	       {
		  /* get aggregated objective value of active variable */
		  SCIP_CALL( SCIPvarGetAggregatedObj(impoperands[v], &obj) );

		  if( SCIPisLE(scip, obj, 0.0) )
		  {
		     SCIP_CALL( SCIPfixVar(scip, impoperands[v], 1.0, &infeasible, &fixed) );

		     *cutoff = *cutoff || infeasible;
		     if( fixed )
			++(*nfixedvars);
		  }
		  else if( !poscontissmall )
		  {
		     SCIP_CALL( SCIPfixVar(scip, impoperands[v], 0.0, &infeasible, &fixed) );
		     assert(!infeasible);
		     assert(fixed);

                     ++(*nfixedvars);
                     zerofix = TRUE;
		  }
                  else
                  {
                     SCIP_Bool redundant;
                     SCIP_Bool aggregated;

                     /* aggregate resultant to operand */
                     SCIP_CALL( SCIPaggregateVars(scip, resvar, impoperands[v], 1.0, -1.0, 0.0,
                           &infeasible, &redundant, &aggregated) );
                     assert(!infeasible);

                     if( aggregated )
                     {
                        /* note that we cannot remove the aggregated operand because we do not know the position */
                        ++(*naggrvars);

                        aggregationperformed = TRUE;

                        SCIPdebugMsg(scip, "dual aggregating operand <%s> with 1 up- and downlock to the resultant <%s> in constraint <%s>\n", SCIPvarGetName(impoperands[v]), SCIPvarGetName(resvar), SCIPconsGetName(cons));
                     }
                  }
	       }
	       assert(*nfixedvars - oldnfixedvars + *naggrvars - oldnaggrvars <= nimpoperands);

               /* did we aggregate the resultant, then we can decide the value to fix it on the (aggregated) objective
                * value since it was a independant variable
                */
	       if( aggregationperformed || zerofix )
	       {
                  SCIP_Real fixval;

                  if( zerofix )
                     fixval = 0.0;
                  else
                  {
                     /* get aggregated objective value of active variable, that might be changed */
                     SCIP_CALL( SCIPvarGetAggregatedObj(resvar, &obj) );
                     assert(!SCIPisPositive(scip, obj));

                     fixval = (SCIPisNegative(scip, obj) ? 1.0 : 0.0);
                  }

                  if( fixval < 0.5 || *nfixedvars - oldnfixedvars + *naggrvars - oldnaggrvars == nvars )
                  {
                     SCIPdebugMsg(scip, "constraint <%s> we can fix the resultant <%s> to %g, because the AND-constraint will alwys be fulfilled\n", SCIPconsGetName(cons), SCIPvarGetName(resvar), fixval);

                     SCIP_CALL( SCIPfixVar(scip, resvar, fixval, &infeasible, &fixed) );
                     assert(!infeasible);
                     assert(fixed);

                     ++(*nfixedvars);

                     SCIPdebugMsg(scip, "deleting constraint <%s> because \n", SCIPconsGetName(cons));

                     SCIP_CALL( SCIPdelCons(scip, cons) );
                     ++(*ndelconss);
                  }
	       }
	    }
	 }
	 /* resultant contributes to the objective with a positive value */
	 else
	 {
	    SCIP_Bool zerofix = FALSE;
#ifndef NDEBUG
	    SCIP_Real tmpobj;

	    assert(nimpoperands > 0);
	    assert(maxpos >= 0 && maxpos <= consdata->nvars);
	    assert(!SCIPisInfinity(scip, -maxobj));
	    SCIP_CALL( SCIPvarGetAggregatedObj(impoperands[maxpos], &tmpobj) );
	    assert(SCIPisEQ(scip, tmpobj, maxobj));
#endif

	    /* if the smallest possible contribution is negative, but does not compensate the positive contribution of
	     * the resultant we need to fix this variable to 0
	     */
	    if( nimpoperands == nvars && SCIPisLE(scip, maxobj, 0.0) )
	    {
	       SCIP_Real fixval = (SCIPisLE(scip, REALABS(maxobj), resobj) ? 0.0 : 1.0);

	       SCIPdebugMsg(scip, "dual-fixing variable <%s> in constraint <%s> to %g, because the contribution is%s enough to nullify/exceed the contribution of the resultant \n", SCIPvarGetName(impoperands[maxpos]), SCIPconsGetName(cons), fixval, (fixval < 0.5) ? " not" : "");

	       SCIP_CALL( SCIPfixVar(scip, impoperands[maxpos], fixval, &infeasible, &fixed) );
	       zerofix = (fixval < 0.5);

	       *cutoff = *cutoff || infeasible;
	       if( fixed )
		  ++(*nfixedvars);
	    }

	    SCIPdebugMsg(scip, "dual-fixing all variables, except the variable with the highest contribution to the objective, in constraint <%s> with positive contribution to 0 and with negative contribution to 1\n", SCIPconsGetName(cons));

	    for( v = nimpoperands - 1; v >= 0 && !(*cutoff); --v )
	    {
	       /* get aggregated objective value of active variable */
	       SCIP_CALL( SCIPvarGetAggregatedObj(impoperands[v], &obj) );

	       if( SCIPisLE(scip, obj, 0.0) )
	       {
		  if( v == maxpos )
		     continue;

		  SCIP_CALL( SCIPfixVar(scip, impoperands[v], 1.0, &infeasible, &fixed) );
	       }
	       else
	       {
		  SCIP_CALL( SCIPfixVar(scip, impoperands[v], 0.0, &infeasible, &fixed) );
		  zerofix = TRUE;
	       }

	       *cutoff = *cutoff || infeasible;
	       if( fixed )
		  ++(*nfixedvars);
	    }
	    assert(*nfixedvars - oldnfixedvars <= nimpoperands);
	    /* iff we have fixed all variables, all variables needed to be stored in the impoperands array */
	    assert((*nfixedvars - oldnfixedvars == nvars) == (nimpoperands == nvars));

	    if( *nfixedvars - oldnfixedvars == nvars )
	    {
	       SCIPdebugMsg(scip, "all operands are fixed in constraint <%s> => fix resultant <%s> to %g\n", SCIPconsGetName(cons), SCIPvarGetName(resvar), (zerofix ? 0.0 : 1.0));

	       SCIP_CALL( SCIPfixVar(scip, resvar, zerofix ? 0.0 : 1.0, &infeasible, &fixed) );

	       *cutoff = *cutoff || infeasible;
	       if( fixed )
		  ++(*nfixedvars);

	       SCIPdebugMsg(scip, "deleting constraint <%s> because all variables are fixed\n", SCIPconsGetName(cons));

	       SCIP_CALL( SCIPdelCons(scip, cons) );
	       ++(*ndelconss);
	    }
	 }
      }
      /* resultant is lock by another constraint (handler), check for operands with only one down- and uplock */
      else
      {
	 SCIP_Real maxobj = -SCIPinfinity(scip);
	 SCIP_Real resobj;
	 SCIP_Real obj;
	 SCIP_Bool redundant;
	 SCIP_Bool aggregated;
	 SCIP_Bool resobjispos;
	 SCIP_Bool linearize = FALSE;
	 SCIP_Bool zerofix = FALSE;
#ifndef NDEBUG
	 int oldnchgcoefs = *nchgcoefs;
	 int oldnfixedvars = *nfixedvars;
#endif

	 /* get aggregated objective value of active variable */
	 SCIP_CALL( SCIPvarGetAggregatedObj(resvar, &resobj) );

	 resobjispos = SCIPisGT(scip, resobj, 0.0);

	 /* we can only aggregate when the objective contribution of the resultant is less or equal to 0 */
	 if( !resobjispos )
	 {
            SCIP_Bool goodvarsfound = FALSE;

	    for( v = nvars - 1; v >= 0; --v )
	    {
	       var = vars[v];
	       assert(var != NULL);
	       assert(SCIPvarGetNLocksUp(var) >= 1 && SCIPvarGetNLocksDown(var) >= 1);

	       /* get aggregated objective value of active variable */
	       SCIP_CALL( SCIPvarGetAggregatedObj(var, &obj) );

	       /* all operands which are only locked by this constraint, the objective contribution is greater or equal
		* to 0 can be aggregated to the resultant
		*/
	       if( SCIPvarGetNLocksUp(var) == 1 && SCIPvarGetNLocksDown(var) == 1 )
               {
                  if( !SCIPisNegative(scip, obj) )
                  {
                     /* aggregate resultant to operand */
                     SCIP_CALL( SCIPaggregateVars(scip, resvar, var, 1.0, -1.0, 0.0,
                           &infeasible, &redundant, &aggregated) );

                     if( aggregated )
                     {
                        ++(*naggrvars);

                        linearize = TRUE;

                        /* delete redundant entry from constraint */
                        SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
                        ++(*nchgcoefs);

                        SCIPdebugMsg(scip, "dual aggregating operand <%s> with 1 up- and downlock to the resultant <%s> in constraint <%s>\n", SCIPvarGetName(var), SCIPvarGetName(resvar), SCIPconsGetName(cons));
                     }

                     *cutoff = *cutoff || infeasible;
                  }
                  else
                     goodvarsfound = TRUE;
               }
	    }
	    assert(*nchgcoefs - oldnchgcoefs <= nvars);

            /* if we aggregated an operands with the resultant we can also fix "good" independant operands to 1, since
             * the correctness of "resultant = 0 => at least one operand = 0" in enforced by that aggregation
             * without an aggregation we cannot fix these variables since it might lead to infeasibility, e.g.
             *
             *   obj(x3) = -1
             *   r = x1 * x2 * x3
             *   r = 0
             *   x1 = 1
             *   x2 = 1
             */
            if( !*cutoff && goodvarsfound && linearize )
            {
               /* fix good variables to 1 */
               for( v = consdata->nvars - 1; v >= 0; --v )
               {
                  var = vars[v];
                  assert(var != NULL);

                  if( SCIPvarGetNLocksUp(var) == 1 && SCIPvarGetNLocksDown(var) == 1 )
                  {
#ifndef NDEBUG
                     /* aggregated objective value of active variable need to be negative */
                     SCIP_CALL( SCIPvarGetAggregatedObj(var, &obj) );
                     assert(SCIPisNegative(scip, obj));
#endif
                     SCIPdebugMsg(scip, "dual-fixing variable <%s> in constraint <%s> to 1, because the contribution is negative\n", SCIPvarGetName(var), SCIPconsGetName(cons));

                     SCIP_CALL( SCIPfixVar(scip, var, 1.0, &infeasible, &fixed) );

                     assert(!infeasible);
                     if( fixed )
                        ++(*nfixedvars);
                  }
               }
               assert(*nfixedvars - oldnfixedvars <= consdata->nvars);
            }
	    assert(*nchgcoefs - oldnchgcoefs + *nfixedvars - oldnfixedvars <= nvars);
	 }
	 /* if the downlocks of the resultant are only from this constraint and the objective contribution is positive,
	  * we can try to fix operands
	  */
	 else if( SCIPvarGetNLocksDown(resvar) == 1 )
	 {
	    SCIP_Bool locksareone = TRUE;
	    int maxpos = -1;

	    for( v = nvars - 1; v >= 0; --v )
	    {
	       var = vars[v];
	       assert(var != NULL);
	       assert(SCIPvarGetNLocksUp(var) >= 1 && SCIPvarGetNLocksDown(var) >= 1);

	       /* check if all resultants are only locked by this constraint */
	       locksareone = locksareone && (SCIPvarGetNLocksUp(var) == 1 && SCIPvarGetNLocksDown(var) == 1);

	       /* get aggregated objective value of active variable */
	       SCIP_CALL( SCIPvarGetAggregatedObj(var, &obj) );

	       /* memorize maximal objective value of operands and its position */
	       if( obj > maxobj )
	       {
		  maxpos = v;
		  maxobj = obj;
	       }

	       /* all operands which are only locked by this constraint, the objective contribution is greater or equal
		* to 0, and the absolute value of the contribution of the resultant exceeds can be eliminated and
		* aggregated to the resultant
	        */
	       if( SCIPvarGetNLocksUp(var) == 1 && SCIPvarGetNLocksDown(var) == 1 && SCIPisGE(scip, obj, 0.0) )
	       {
		  SCIPdebugMsg(scip, "dualfix operand <%s> in constraint <%s> to 0\n", SCIPvarGetName(var), SCIPconsGetName(cons));

		  SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );

		  *cutoff = *cutoff || infeasible;
		  if( fixed )
		     ++(*nfixedvars);

		  zerofix = TRUE;
	       }
	    }
	    assert(*nchgcoefs - oldnchgcoefs <= nvars);

	    /* if constraint is still active and all operands are only lock by this constraint, we check if we can fix
	     * the worst (in objective contribution) operand to zero
	     */
	    if( !zerofix && locksareone && SCIPisGE(scip, resobj, REALABS(maxobj)) )
	    {
	       assert(!zerofix);
	       /* objective contribution needs to be negative, otherwise, the variable should already be fixed to 0 */
	       assert(SCIPisLT(scip, maxobj, 0.0));

	       SCIPdebugMsg(scip, "dualfix operand <%s> with worst contribution in constraint <%s> to 0\n", SCIPvarGetName(vars[maxpos]), SCIPconsGetName(cons));

	       SCIP_CALL( SCIPfixVar(scip, vars[maxpos], 0.0, &infeasible, &fixed) );

	       *cutoff = *cutoff || infeasible;
	       if( fixed )
		  ++(*nfixedvars);

	       zerofix = TRUE;
	    }

	    /* fix the resultant if one operand was fixed to zero and delete the constraint */
	    if( zerofix )
	    {
	       SCIPdebugMsg(scip, "fix resultant <%s> in constraint <%s> to 0\n", SCIPvarGetName(resvar), SCIPconsGetName(cons));

	       SCIP_CALL( SCIPfixVar(scip, resvar, 0.0, &infeasible, &fixed) );

	       *cutoff = *cutoff || infeasible;
	       if( fixed )
		  ++(*nfixedvars);

	       SCIPdebugMsg(scip, "deleting constraint <%s> because at least one operand and the resultant is fixed to zero\n", SCIPconsGetName(cons));

	       SCIP_CALL( SCIPdelCons(scip, cons) );
	       ++(*ndelconss);
	    }
	 }

         /* we have to linearize the constraint, otherwise we might get wrong propagations, since due to aggregations a
          * resultant fixed to zero is already fulfilling the constraint, and we must not ensure that some remaining
          * operand needs to be 0
          */
	 if( linearize )
	 {
	    SCIP_CONS* newcons;
	    char consname[SCIP_MAXSTRLEN];
	    SCIP_VAR* consvars[2];
	    SCIP_Real vals[2];

	    assert(SCIPconsIsActive(cons));

	    consvars[0] = consdata->resvar;
            vals[0] = 1.0;
            vals[1] = -1.0;

            /* create operator linear constraints */
            for( v = consdata->nvars - 1; v >= 0; --v )
            {
               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), v);
               consvars[1] = consdata->vars[v];

               SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, 2, consvars, vals, -SCIPinfinity(scip), 0.0,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                     SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                     SCIPconsIsStickingAtNode(cons)) );

               /* add constraint */
               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            }
	    (*naddconss) += consdata->nvars;

	    SCIPdebugMsg(scip, "deleting constraint <%s> because it was linearized\n", SCIPconsGetName(cons));

	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);
	 }
	 /* if only one operand is leftover, aggregate it to the resultant */
	 else if( consdata->nvars == 1 )
	 {
	    SCIPdebugMsg(scip, "aggregating last operand <%s> to the resultant <%s> in constraint <%s>\n", SCIPvarGetName(consdata->vars[0]), SCIPvarGetName(resvar), SCIPconsGetName(cons));

	    /* aggregate resultant to operand */
	    SCIP_CALL( SCIPaggregateVars(scip, resvar, consdata->vars[0], 1.0, -1.0, 0.0,
		  &infeasible, &redundant, &aggregated) );

	    if( aggregated )
	       ++(*naggrvars);

	    *cutoff = *cutoff || infeasible;

	    SCIPdebugMsg(scip, "deleting constraint <%s> because all variables are removed\n", SCIPconsGetName(cons));

	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);
	 }

	 /* if no operand is leftover delete the constraint */
	 if( SCIPconsIsActive(cons) && consdata->nvars == 0 )
	 {
	    SCIPdebugMsg(scip, "deleting constraint <%s> because all variables are removed\n", SCIPconsGetName(cons));

	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);
	 }
      }
   }

   SCIPfreeBufferArray(scip, &impoperands);

   return SCIP_OKAY;
}

/** 1. check if at least two operands or one operand and the resultant are in one clique, if so, we can fix the
 *     resultant to zero and in the former case we can also delete this constraint but we need to extract the clique
 *     information as constraint
 *
 *     x == AND(y, z) and clique(y,z) => x = 0, delete constraint and create y + z <= 1
 *     x == AND(y, z) and clique(x,y) => x = 0
 *
 *     special handled cases are:
 *     - if the resultant is a negation of an operand, in that case we fix the resultant to 0
 *     - if the resultant is equal to an operand, we will linearize this constraint by adding all necessary
 *       set-packing constraints like resultant + ~operand <= 1 and delete the old constraint
 *
 *     x == AND(~x, y) => x = 0
 *     x == AND(x, y)  => add x + ~y <= 1 and delete the constraint
 *
 *  2. check if one operand is in a clique with the negation of all other operands, this means we can aggregate this
 *     operand to the resultant
 *
 *     r == AND(x,y,z) and clique(x,~y) and clique(x,~z) => r == x
 *
 *  3. check if the resultant and the negations of all operands are in a clique
 *
 *     r == AND(x,y) and clique(r, ~x,~y) => upgrade the constraint to a set-partitioning constraint r + ~x + ~y = 1
 *
 *  @note We removed also fixed variables and propagate them, and if only one operand is remaining due to removal, we
 *        will aggregate the resultant with this operand
 */
static
SCIP_RETCODE cliquePresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  naggrvars,          /**< pointer to add up the number of aggregated variables */
   int*                  nchgcoefs,          /**< pointer to add up the number of changed coefficients */
   int*                  ndelconss,          /**< pointer to add up the number of deleted constraints */
   int*                  naddconss           /**< pointer to add up the number of added constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int nvars;
   int vstart;
   int vend;
   int v;
   int v2;
   SCIP_Bool negated;
   SCIP_Bool value1;
   SCIP_Bool value2;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   SCIP_Bool allnegoperandsexist;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgcoefs != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPconsIsActive(cons) || SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   vars = consdata->vars;
   nvars = consdata->nvars;
   assert(vars != NULL || nvars == 0);

   /* remove fixed variables to be able to ask for cliques
    *
    * if an operand is fixed to 0 fix the resultant to 0 and delete the constraint
    * if an operand is fixed to 1 remove it from the constraint
    */
   for( v = nvars - 1; v >= 0; --v )
   {
      assert(vars != NULL);

      if( SCIPvarGetLbGlobal(vars[v]) > 0.5 )
      {
	 SCIPdebugMsg(scip, "In constraint <%s> the operand <%s> is fixed to 1 so remove it from the constraint\n",
	    SCIPconsGetName(cons), SCIPvarGetName(vars[v]));

	 /* because we loop from back to front we can delete the entry in the consdata structure */
	 SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
	 ++(*nchgcoefs);

	 assert(consdata->vars == vars);

	 continue;
      }
      else if( SCIPvarGetUbGlobal(vars[v]) < 0.5 )
      {
	 SCIPdebugMsg(scip, "constraint <%s> redundant: because operand <%s> is fixed to zero so we can fix the resultant <%s> to 0\n",
	    SCIPconsGetName(cons), SCIPvarGetName(vars[v]), SCIPvarGetName(consdata->resvar));

	 SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 0.0, &infeasible, &fixed) );
	 *cutoff = *cutoff || infeasible;
	 if( fixed )
	    ++(*nfixedvars);

	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);

	 return SCIP_OKAY;
      }
   }

   /* if we deleted some operands constraint might be redundant */
   if( consdata->nvars < nvars )
   {
      assert(vars == consdata->vars);

      /* all operands fixed to one were removed, so if no operand is left this means we can fix the resultant to 1
       * too
       */
      if( consdata->nvars == 0 )
      {
	 SCIPdebugMsg(scip, "All operand in constraint <%s> were deleted, so the resultant needs to be fixed to 1\n",
	    SCIPconsGetName(cons));

	 SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 1.0, &infeasible, &fixed) );
	 *cutoff = *cutoff || infeasible;
	 if( fixed )
	    ++(*nfixedvars);

	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);

	 return SCIP_OKAY;
      }
      /* if only one not fixed operand is left, we can aggregate it to the resultant */
      else if( consdata->nvars == 1 )
      {
	 SCIP_Bool redundant;
	 SCIP_Bool aggregated;

	 /* aggregate resultant to last operand */
	 SCIP_CALL( SCIPaggregateVars(scip, consdata->resvar, consdata->vars[0], 1.0, -1.0, 0.0,
	       &infeasible, &redundant, &aggregated) );

	 if( aggregated )
	    ++(*naggrvars);

	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);

	 *cutoff = *cutoff || infeasible;

	 return SCIP_OKAY;
      }

      nvars = consdata->nvars;
   }

   /* @todo when cliques are improved, we only need to collect all clique-ids for all variables and check for doubled
    *       entries
    */
   /* case 1 first part */
   /* check if two operands are in a clique */
   for( v = nvars - 1; v > 0; --v )
   {
      assert(vars != NULL);

      var1 = vars[v];
      assert(var1 != NULL);
      negated = FALSE;

      SCIP_CALL( SCIPvarGetProbvarBinary(&var1, &negated) );
      assert(var1 != NULL);

      if( negated )
	 value1 = FALSE;
      else
	 value1 = TRUE;

      assert(SCIPvarGetStatus(var1) != SCIP_VARSTATUS_FIXED);

      for( v2 = v - 1; v2 >= 0; --v2 )
      {
	 var2 = vars[v2];
	 assert(var2 != NULL);

	 negated = FALSE;
	 SCIP_CALL( SCIPvarGetProbvarBinary(&var2, &negated) );
	 assert(var2 != NULL);

	 if( negated )
	    value2 = FALSE;
	 else
	    value2 = TRUE;

	 assert(SCIPvarGetStatus(var2) != SCIP_VARSTATUS_FIXED);

	 /* if both variables are negated of each other or the same, this will be handled in applyFixings();
	  * @note if both variables are the same, then SCIPvarsHaveCommonClique() will return TRUE, so we better
	  *       continue
	  */
	 if( var1 == var2 )
	    continue;

	 if( SCIPvarsHaveCommonClique(var1, value1, var2, value2, TRUE) )
	 {
	    SCIP_CONS* cliquecons;
	    SCIP_VAR* consvars[2];
	    char name[SCIP_MAXSTRLEN];

	    SCIPdebugMsg(scip, "constraint <%s> redundant: because variable <%s> and variable <%s> are in a clique, the resultant <%s> can be fixed to 0\n",
	       SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2), SCIPvarGetName(consdata->resvar));

	    SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 0.0, &infeasible, &fixed) );
	    *cutoff = *cutoff || infeasible;
	    if( fixed )
	       ++(*nfixedvars);


	    /* create clique constraint which lead to the last fixing */
	    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_clq", SCIPconsGetName(cons), v2);

	    if( value1 )
	       consvars[0] = var1;
	    else
	    {
	       SCIP_CALL( SCIPgetNegatedVar(scip, var1, &(consvars[0])) );
	    }

	    if( value2 )
	       consvars[1] = var2;
	    else
	    {
	       SCIP_CALL( SCIPgetNegatedVar(scip, var2, &(consvars[1])) );
	    }

	    SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, name, 2, consvars,
		  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
		  consdata->checkwhenupgr || SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
		  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
                  !(consdata->notremovablewhenupgr) && SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
	    SCIPdebugMsg(scip, " -> adding clique constraint: ");
	    SCIPdebugPrintCons(scip, cliquecons, NULL);
	    SCIP_CALL( SCIPaddCons(scip, cliquecons) );
	    SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
	    ++(*naddconss);

	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);

	    return SCIP_OKAY;
	 }
      }
   }

   var1 = consdata->resvar;
   assert(var1 != NULL);

   negated = FALSE;
   SCIP_CALL( SCIPvarGetProbvarBinary(&var1, &negated) );
   assert(var1 != NULL);

   /* it may appear that we have a fixed resultant */
   if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_FIXED )
   {
      /* resultant is fixed to 1, so fix all operands to 1 */
      if( SCIPvarGetLbGlobal(consdata->resvar) > 0.5 )
      {
	 SCIPdebugMsg(scip, "In constraint <%s> the resultant <%s> is fixed to 1 so fix all operands to 1\n",
	    SCIPconsGetName(cons), SCIPvarGetName(consdata->resvar));

	 /* fix all operands to 1 */
	 for( v = nvars - 1; v >= 0 && !(*cutoff); --v )
	 {
            assert(vars != NULL);

	    SCIPdebugMsg(scip, "Fixing operand <%s> to 1.\n", SCIPvarGetName(vars[v]));

	    SCIP_CALL( SCIPfixVar(scip, vars[v], 1.0, &infeasible, &fixed) );
	    *cutoff = *cutoff || infeasible;

	    if( fixed )
	       ++(*nfixedvars);
	 }

	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);
      }
      /* the upgrade to a linear constraint because of the to 0 fixed resultant we do in propagateCons() */
      else
	 assert(SCIPvarGetUbGlobal(consdata->resvar) < 0.5);

      return SCIP_OKAY;
   }

   if( negated )
      value1 = FALSE;
   else
      value1 = TRUE;

   /* case 1 second part */
   /* check if one operands is in a clique with the resultant */
   for( v = nvars - 1; v >= 0; --v )
   {
      assert(vars != NULL);

      var2 = vars[v];
      assert(var2 != NULL);

      negated = FALSE;
      SCIP_CALL( SCIPvarGetProbvarBinary(&var2, &negated) );
      assert(var2 != NULL);

      if( negated )
	 value2 = FALSE;
      else
	 value2 = TRUE;

      /* if both variables are negated of each other or the same, this will be handled in applyFixings();
       * @note if both variables are the same, then SCIPvarsHaveCommonClique() will return TRUE, so we better continue
       */
      if( var1 == var2 )
      {
	 /* x1 == AND(~x1, x2 ...) => x1 = 0 */
	 if( value1 != value2 )
	 {
	    SCIPdebugMsg(scip, "In constraint <%s> the resultant <%s> can be fixed to 0 because the negation of it is an operand.\n",
	       SCIPconsGetName(cons), SCIPvarGetName(consdata->resvar));

	    SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 0.0, &infeasible, &fixed) );
	    *cutoff = *cutoff || infeasible;

	    if( fixed )
	       ++(*nfixedvars);

	    return SCIP_OKAY;
	 }
	 /* x1 == AND(x1, x2 ...) => delete constraint and create all set-packing constraints x1 + ~x2 <= 1, x1 + ~... <= 1 */
	 else
	 {
	    SCIP_CONS* cliquecons;
	    SCIP_VAR* consvars[2];
	    char name[SCIP_MAXSTRLEN];

	    assert(value1 == value2);

	    consvars[0] = consdata->resvar;

	    for( v2 = nvars - 1; v2 >= 0; --v2 )
	    {
	       var2 = vars[v2];
	       negated = FALSE;
	       SCIP_CALL( SCIPvarGetProbvarBinary(&var2, &negated) );

	       /* if the active representations of the resultant and an operand are different then we need to extract
		* this as a clique constraint
		*
		* if the active representations of the resultant and an operand are equal then the clique constraint
		* would look like x1 + ~x1 <= 1, which is redundant
		*
		* if the active representations of the resultant and an operand are negated of each other then the
		* clique constraint would look like x1 + x1 <= 1, which will lead to a fixation of the resultant later
		* on
		*/
	       if( var1 == var2 )
	       {
		  if( value1 == negated )
		  {
		     SCIPdebugMsg(scip, "In constraint <%s> the resultant <%s> can be fixed to 0 because the negation of it is an operand.\n",
			SCIPconsGetName(cons), SCIPvarGetName(consdata->resvar));

		     SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 0.0, &infeasible, &fixed) );
		     *cutoff = *cutoff || infeasible;

		     if( fixed )
			++(*nfixedvars);

		     break;
		  }
	       }
	       else
	       {
		  SCIP_CALL( SCIPgetNegatedVar(scip, vars[v2], &consvars[1]) );
		  assert(consvars[1] != NULL);

                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_clq_%d", SCIPconsGetName(cons), v2);

                  SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, name, 2, consvars,
                        SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                        consdata->checkwhenupgr || SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                        SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
                        !(consdata->notremovablewhenupgr) && SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
                  SCIPdebugMsg(scip, " -> adding clique constraint: ");
                  SCIPdebugPrintCons(scip, cliquecons, NULL);
                  SCIP_CALL( SCIPaddCons(scip, cliquecons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
                  ++(*naddconss);
	       }
	    }

	    /* delete old constraint */
	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);

	    return SCIP_OKAY;
	 }
      }

      /* due to SCIPvarsHaveCommonClique() returns on two same variables that they are in a clique, we need to handle
       * it explicitly
       */
      if( var1 == var2 && value1 == value2 )
	 continue;

      /* due to SCIPvarsHaveCommonClique() returns on two negated variables that they are not in a clique, we need to
       * handle it explicitly
       */
      if( (var1 == var2 && value1 != value2) || SCIPvarsHaveCommonClique(var1, value1, var2, value2, TRUE) )
      {
	 SCIPdebugMsg(scip, "in constraint <%s> the resultant <%s> can be fixed to 0 because it is in a clique with operand <%s>\n",
	    SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));

	 SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 0.0, &infeasible, &fixed) );
	 *cutoff = *cutoff || infeasible;
	 if( fixed )
	    ++(*nfixedvars);

	 return SCIP_OKAY;
      }
   }

   if( !SCIPconsIsActive(cons) )
      return SCIP_OKAY;

   v2 = -1;
   /* check which operands have a negated variable */
   for( v = nvars - 1; v >= 0; --v )
   {
      assert(vars != NULL);

      var1 = vars[v];
      assert(var1 != NULL);

      negated = FALSE;
      SCIP_CALL( SCIPvarGetProbvarBinary(&var1, &negated) );
      assert(var1 != NULL);

      if( SCIPvarGetNegatedVar(var1) == NULL )
      {
	 if( v2 >= 0 )
	    break;
	 v2 = v;
      }
   }

   allnegoperandsexist = FALSE;

   /* all operands have a negated variable, so we will check for all possible negated ciques */
   if( v2 == -1 )
   {
      allnegoperandsexist = TRUE;
      vstart = nvars - 1;
      vend = 0;
   }
   /* exactly one operands has no negated variable, so only this variable can be in a clique with all other negations */
   else if( v2 >= 0 && v == -1 )
   {
      vstart = v2;
      vend = v2;
   }
   /* at least two operands have no negated variable, so there is no possible clique with negated variables */
   else
   {
      vstart = -1;
      vend = 0;
   }

   /* case 2 */
   /* check for negated cliques in the operands */
   for( v = vstart; v >= vend; --v )
   {
      assert(vars != NULL);

      var1 = vars[v];
      assert(var1 != NULL);

      negated = FALSE;
      SCIP_CALL( SCIPvarGetProbvarBinary(&var1, &negated) );
      assert(var1 != NULL);

      if( negated )
	 value1 = FALSE;
      else
	 value1 = TRUE;

      for( v2 = nvars - 1; v2 >= 0; --v2 )
      {
	 if( v2 == v )
	    continue;

	 var2 = vars[v2];
	 assert(var2 != NULL);

	 negated = FALSE;
	 SCIP_CALL( SCIPvarGetProbvarBinary(&var2, &negated) );
	 assert(var2 != NULL);

	 if( negated )
	    value2 = FALSE;
	 else
	    value2 = TRUE;

	 assert(SCIPvarGetNegatedVar(var2) != NULL);

	 /* invert flag, because we want to check var 1 against all negations of the other variables */
	 value2 = !value2;

	 /* due to SCIPvarsHaveCommonClique() returns on two same variables that they are in a clique, we need to handle
	  * it explicitly
	  */
	 if( var1 == var2 && value1 == value2 )
	 {
	    SCIPdebugMsg(scip, "in constraint <%s> the resultant <%s> can be fixed to 0 because two operands are negated of each other\n",
	       SCIPconsGetName(cons), SCIPvarGetName(consdata->resvar));

	    SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 0.0, &infeasible, &fixed) );
	    *cutoff = *cutoff || infeasible;
	    if( fixed )
	       ++(*nfixedvars);

	    return SCIP_OKAY;
	 }

	 /* due to SCIPvarsHaveCommonClique() returns on two negated variables that they are not in a clique, we need to
	  * handle it explicitly
	  */
	 if( var1 == var2 && value1 != value2 )
	    continue;

	 if( !SCIPvarsHaveCommonClique(var1, value1, var2, value2, TRUE) )
	    break;
      }

      if( v2 == -1 )
      {
	 SCIP_Bool redundant;
	 SCIP_Bool aggregated;

	 SCIPdebugMsg(scip, "In constraint <%s> the operand <%s> is in a negated clique with all other operands, so we can aggregated this operand to the resultant <%s>.\n",
	    SCIPconsGetName(cons), SCIPvarGetName(vars[v]), SCIPvarGetName(consdata->resvar));

	 SCIP_CALL( SCIPaggregateVars(scip, consdata->resvar, vars[v], 1.0, -1.0, 0.0,
	       &infeasible, &redundant, &aggregated) );
	 *cutoff = *cutoff || infeasible;

	 if( aggregated )
	    ++(*naggrvars);

	 return SCIP_OKAY;
      }
   }

   /* case 3 */
   /* check if the resultant and the negations of the operands are in a clique, then we can upgrade this constraint to a
    * set-partitioning constraint
    */
   if( allnegoperandsexist && SCIPconsIsActive(cons) )
   {
      SCIP_VAR** newvars;
      SCIP_Bool* negations;
      SCIP_Bool upgrade;

      SCIP_CALL( SCIPallocBufferArray(scip, &newvars, nvars + 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &negations, nvars + 1) );
      BMSclearMemoryArray(negations, nvars + 1);

      var1 = consdata->resvar;
      SCIP_CALL( SCIPvarGetProbvarBinary(&var1, &negations[nvars]) );
      assert(var1 != NULL);
      assert(SCIPvarGetStatus(var1) != SCIP_VARSTATUS_FIXED);

      newvars[nvars] = var1;

      /* get active variables */
      for( v = nvars - 1; v >= 0; --v )
      {
	 assert(vars != NULL);

	 var1 = vars[v];
	 SCIP_CALL( SCIPvarGetProbvarBinary(&var1, &negations[v]) );
	 assert(var1 != NULL);
	 assert(SCIPvarGetStatus(var1) != SCIP_VARSTATUS_FIXED);

	 newvars[v] = var1;

	 /* there should be no variable left that is equal or negated to the resultant */
	 assert(newvars[v] != newvars[nvars]);
      }

      upgrade = TRUE;

      /* the resultant is in a clique with the negations of all operands, due to this AND-constraint */
      /* only check if the negations of all operands are in a clique */
      for( v = nvars - 1; v >= 0 && upgrade; --v )
      {
	 for( v2 = v - 1; v2 >= 0; --v2 )
	 {
	    /* the resultant need to be in a clique with the negations of all operands */
	    if( !SCIPvarsHaveCommonClique(newvars[v], negations[v], newvars[v2], negations[v2], TRUE) )
	    {
	       upgrade = FALSE;
	       break;
	    }
	 }
      }

      /* all variables are in a clique, so upgrade thi AND-constraint */
      if( upgrade )
      {
	 SCIP_CONS* cliquecons;
	 char name[SCIP_MAXSTRLEN];

	 /* get new constraint variables */
	 if( negations[nvars] )
	 {
	    /* negation does not need to be existing, so SCIPvarGetNegatedVar() cannot be called
	     * (e.g. resultant = ~x = 1 - x and x = y = newvars[nvars] and negations[nvars] = TRUE,
	     *  then y does not need to have a negated variable, yet)
	     */
	    SCIP_CALL( SCIPgetNegatedVar(scip, newvars[nvars], &(newvars[nvars])) );
	 }
	 assert(newvars[nvars] != NULL);

	 for( v = nvars - 1; v >= 0; --v )
	 {
	    if( !negations[v] )
	    {
	       /* negation does not need to be existing, so SCIPvarGetNegatedVar() cannot be called
		* (e.g. vars[v] = ~x = 1 - x and x = y = newvars[v] and negations[v] = TRUE,
		*  then y does not need to have a negated variable, yet)
		*/
	       SCIP_CALL( SCIPgetNegatedVar(scip, newvars[v], &(newvars[v])) );
	    }
	    assert(newvars[v] != NULL);
	 }

	 (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_clqeq", SCIPconsGetName(cons));

	 SCIP_CALL( SCIPcreateConsSetpart(scip, &cliquecons, name, nvars + 1, newvars,
	       SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
	       consdata->checkwhenupgr || SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
	       SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
               !(consdata->notremovablewhenupgr) && SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
	 SCIPdebugMsg(scip, " -> upgrading AND-constraint <%s> with use of clique information to a set-partitioning constraint: \n", SCIPconsGetName(cons));
	 SCIPdebugPrintCons(scip, cliquecons, NULL);
	 SCIP_CALL( SCIPaddCons(scip, cliquecons) );
	 SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
	 ++(*naddconss);

	 /* delete old constraint */
	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);
      }

      SCIPfreeBufferArray(scip, &negations);
      SCIPfreeBufferArray(scip, &newvars);
   }

   return SCIP_OKAY;
}

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyAndcons)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqAndcons)
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   SCIP_Bool coefsequal;
   int i;
#ifndef NDEBUG
   SCIP* scip;

   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);

   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   /* sorts the constraints */
   consdataSort(consdata1);
   consdataSort(consdata2);
   assert(consdata1->sorted);
   assert(consdata2->sorted);

   coefsequal = TRUE;

   for( i = 0; i < consdata1->nvars ; ++i )
   {
      /* tests if variables are equal */
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 ||
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         coefsequal = FALSE;
         break;
      }
      assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 0);
   }

   return coefsequal;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValAndcons)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int minidx;
   int mididx;
   int maxidx;

   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->sorted);
   assert(consdata->nvars > 0);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   assert(minidx >= 0 && minidx <= maxidx);

   return SCIPhashTwo(SCIPcombineTwoInt(consdata->nvars, minidx),
                      SCIPcombineTwoInt(mididx, maxidx));
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint 
 *  accordingly; in contrast to removeRedundantConstraints(), it uses a hash table 
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  firstchange,        /**< pointer to store first changed constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int c;

   assert(conss != NULL);
   assert(ndelconss != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = nconss;
   hashtablesize = MAX(hashtablesize, HASHSIZE_ANDCONS);
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyAndcons, hashKeyEqAndcons, hashKeyValAndcons, (void*) scip) );

   *cutoff = FALSE;

   /* check all constraints in the given set for redundancy */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata0;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      consdata0 = SCIPconsGetData(cons0);

      /* sort the constraint */
      consdataSort(consdata0);
      assert(consdata0->sorted);

      /* get constraint from current hash table with same variables as cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));

      if( cons1 != NULL )
      {
         SCIP_CONSDATA* consdata1;
         SCIP_Bool redundant;

         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));

         consdata1 = SCIPconsGetData(cons1);

         assert(consdata1 != NULL);
         assert(consdata0->nvars >= 1 && consdata0->nvars == consdata1->nvars);

         assert(consdata0->sorted && consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);

         redundant = FALSE;

         if( consdata0->resvar != consdata1->resvar )
         {
            SCIP_Bool aggregated;

            assert(SCIPvarCompare(consdata0->resvar, consdata1->resvar) != 0);

            /* aggregate resultants */
            SCIP_CALL( SCIPaggregateVars(scip, consdata0->resvar, consdata1->resvar, 1.0, -1.0, 0.0,
                  cutoff, &redundant, &aggregated) );
            assert(redundant || SCIPdoNotAggr(scip));

            if( aggregated )
               ++(*naggrvars);
            if( *cutoff )
               goto TERMINATE;
         }
         else
            redundant = TRUE;

         /* delete consdel */
         if( redundant )
         {
            /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
            SCIP_CALL( SCIPupdateConsFlags(scip, cons1, cons0) );

	    /* also take the check when upgrade flag over if necessary */
	    consdata1->checkwhenupgr = consdata1->checkwhenupgr || consdata0->checkwhenupgr;
	    consdata1->notremovablewhenupgr = consdata1->notremovablewhenupgr || consdata0->notremovablewhenupgr;

            SCIP_CALL( SCIPdelCons(scip, cons0) );
            (*ndelconss)++;
         }

         /* update the first changed constraint to begin the next aggregation round with */
         if( consdata0->changed && SCIPconsGetPos(cons1) < *firstchange )
            *firstchange = SCIPconsGetPos(cons1);

         assert(SCIPconsIsActive(cons1));
      }
      else
      {
         /* no such constraint in current hash table: insert cons0 into hash table */
         SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
      }
   }
 TERMINATE:
   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}

/** helper function to enforce constraints */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool separated;
   SCIP_Bool violated;
   SCIP_Bool cutoff;
   int i;

   separated = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, FALSE, FALSE, &violated) );
      if( violated )
      {
         if( conshdlrdata->enforcecuts )
         {
            SCIP_Bool consseparated;

            SCIP_CALL( separateCons(scip, conss[i], sol, &consseparated, &cutoff) );
            if ( cutoff )
            {
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            separated = separated || consseparated;

            /* following assert is wrong in the case some variables were not in relaxation (dynamic columns),
            *
            * e.g. the resultant, which has a negative objective value, is in the relaxation solution on its upper bound
            * (variables with status loose are in an relaxation solution on it's best bound), but already creating a
            * row, and thereby creating the column, changes the solution value (variable than has status column, and the
            * initialization sets the relaxation solution value) to 0.0, and this already could lead to no violation of
            * the rows, which then are not seperated into the lp
            */
#ifdef SCIP_DISABLED_CODE
            assert(consseparated); /* because the solution is integral, the separation always finds a cut */
#endif
         }
         else
         {
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
      }
   }

   if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** compares constraint with all prior constraints for possible redundancy or aggregation,
 *  and removes or changes constraint accordingly
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  nbdchgs,            /**< pointer to count the number of performed bound changes, or NULL */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   SCIP_Bool cons0changed;
   int c;

   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(cutoff != NULL);
   assert(naggrvars != NULL);
   assert(nbdchgs != NULL);
   assert(ndelconss != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);

   /* sort the constraint */
   consdataSort(consdata0);

   assert(consdata0->nvars >= 1);
   assert(consdata0->sorted);

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;

   if( SCIPconsIsActive(cons0) )
   {
      for( c = (cons0changed ? 0 : firstchange); c < chkind && !(*cutoff); ++c )
      {
         SCIP_CONS* cons1;
         SCIP_CONSDATA* consdata1;
         SCIP_Bool cons0superset;
         SCIP_Bool cons1superset;
         int v0;
         int v1;

         if( c % 1000 == 0 && SCIPisStopped(scip) )
            break;

         cons1 = conss[c];

         /* ignore inactive and modifiable constraints */
         if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
            continue;

         consdata1 = SCIPconsGetData(cons1);
         assert(consdata1 != NULL);

#ifdef SCIP_DISABLED_CODE
         SCIPdebugMsg(scip, "preprocess AND-constraint pair <%s>[chg:%d] and <%s>[chg:%d]\n",
            SCIPconsGetName(cons0), cons0changed, SCIPconsGetName(cons1), consdata1->changed);
#endif

         /* if both constraints were not changed since last round, we can ignore the pair */
         if( !cons0changed && !consdata1->changed )
            continue;

         assert(consdata1->nvars >= 1);

         /* sort the constraint */
         consdataSort(consdata1);
         assert(consdata1->sorted);

         /* check consdata0 against consdata1:
          * - if they consist of the same operands, the resultants can be aggregated
          * - if one operand list is a subset of the other, add implication r0 = 1 -> r1 = 1, or r1 = 1 -> r0 = 1
          */
         v0 = 0;
         v1 = 0;
         cons0superset = TRUE;
         cons1superset = TRUE;
         while( (v0 < consdata0->nvars || v1 < consdata1->nvars) && (cons0superset || cons1superset) )
         {
            int varcmp;

            /* test, if variable appears in only one or in both constraints */
            if( v0 < consdata0->nvars && v1 < consdata1->nvars )
               varcmp = SCIPvarCompare(consdata0->vars[v0], consdata1->vars[v1]);
            else if( v0 < consdata0->nvars )
               varcmp = -1;
            else
               varcmp = +1;

            switch( varcmp )
            {
            case -1:
               /* variable doesn't appear in consdata1 */
               cons1superset = FALSE;
               v0++;
               break;

            case +1:
               /* variable doesn't appear in consdata0 */
               cons0superset = FALSE;
               v1++;
               break;

            case 0:
               /* variable appears in both constraints */
               v0++;
               v1++;
               break;

            default:
               SCIPerrorMessage("invalid comparison result\n");
               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
         }

         /* check for equivalence and domination */
         if( cons0superset && cons1superset )
         {
            SCIP_Bool infeasible;
            SCIP_Bool redundant;
            SCIP_Bool aggregated;

            /* constraints are equivalent */
            SCIPdebugMsg(scip, "equivalent AND-constraints <%s> and <%s>: aggregate resultants <%s> == <%s>\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), SCIPvarGetName(consdata0->resvar),
               SCIPvarGetName(consdata1->resvar));

            /* aggregate resultants */
            SCIP_CALL( SCIPaggregateVars(scip, consdata0->resvar, consdata1->resvar, 1.0, -1.0, 0.0,
                  &infeasible, &redundant, &aggregated) );
            assert(redundant || SCIPdoNotAggr(scip));

            if( aggregated )
            {
               assert(redundant);
               (*naggrvars)++;
            }

            if( redundant )
            {
	       /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
	       SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );

	       /* also take the check when upgrade flag over if necessary */
               consdata0->checkwhenupgr = consdata1->checkwhenupgr || consdata0->checkwhenupgr;
               consdata0->notremovablewhenupgr = consdata1->notremovablewhenupgr || consdata0->notremovablewhenupgr;

               /* delete constraint */
               SCIP_CALL( SCIPdelCons(scip, cons1) );
               (*ndelconss)++;
            }

            *cutoff = *cutoff || infeasible;
         }
         else if( cons0superset )
         {
            SCIP_Bool infeasible;
            int nboundchgs;

            /* the conjunction of cons0 is a superset of the conjunction of cons1 */
            SCIPdebugMsg(scip, "AND-constraint <%s> is superset of <%s>: add implication <%s> = 1 -> <%s> = 1\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), SCIPvarGetName(consdata0->resvar),
               SCIPvarGetName(consdata1->resvar));

            /* add implication */
            SCIP_CALL( SCIPaddVarImplication(scip, consdata0->resvar, TRUE, consdata1->resvar, SCIP_BOUNDTYPE_LOWER, 1.0,
                  &infeasible, &nboundchgs) );
            *cutoff = *cutoff || infeasible;
            (*nbdchgs) += nboundchgs;
         }
         else if( cons1superset )
         {
            SCIP_Bool infeasible;
            int nboundchgs;

            /* the conjunction of cons1 is a superset of the conjunction of cons0 */
            SCIPdebugMsg(scip, "AND-constraint <%s> is superset of <%s>: add implication <%s> = 1 -> <%s> = 1\n",
               SCIPconsGetName(cons1), SCIPconsGetName(cons0), SCIPvarGetName(consdata1->resvar),
               SCIPvarGetName(consdata0->resvar));

            /* add implication */
            SCIP_CALL( SCIPaddVarImplication(scip, consdata1->resvar, TRUE, consdata0->resvar, SCIP_BOUNDTYPE_LOWER, 1.0,
                  &infeasible, &nboundchgs) );
            *cutoff = *cutoff || infeasible;
            (*nbdchgs) += nboundchgs;
         }
      }
   }
   consdata0->changed = FALSE;

   return SCIP_OKAY;
}

/** tries to reformulate an expression graph node that is a product of binary variables via introducing an AND-constraint */
static
SCIP_DECL_EXPRGRAPHNODEREFORM(exprgraphnodeReformAnd)
{
   SCIP_EXPRGRAPHNODE* child;
   char name[SCIP_MAXSTRLEN];
   int nchildren;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int c;

   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(naddcons != NULL);
   assert(reformnode != NULL);

   *reformnode = NULL;

   /* allow only products given as EXPR_PRODUCT or EXPR_POLYNOMIAL with only 1 monomial */
   if( SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_PRODUCT &&
       (SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_POLYNOMIAL || SCIPexprgraphGetNodePolynomialNMonomials(node) > 1)
     )
      return SCIP_OKAY;

   nchildren = SCIPexprgraphGetNodeNChildren(node);

   /* for a polynomial with only one monomial, all children should appear as factors in the monomial
    * since we assume that the factors have been merged, this means that the number of factors in the monomial should equal the number of children of the node
    */
   assert(SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_POLYNOMIAL || SCIPexprGetMonomialNFactors(SCIPexprgraphGetNodePolynomialMonomials(node)[0]) == nchildren);

   /* check only products with at least 3 variables (2 variables are taken of by cons_quadratic) */
   if( nchildren <= 2 )
      return SCIP_OKAY;

   /* check if all factors correspond to binary variables, and if so, setup vars array */
   for( c = 0; c < nchildren; ++c )
   {
      child = SCIPexprgraphGetNodeChildren(node)[c];

      if( SCIPexprgraphGetNodeOperator(child) != SCIP_EXPR_VARIDX )
         return SCIP_OKAY;

      var = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, child);
      if( !SCIPvarIsBinary(var) )
         return SCIP_OKAY;
   }

   /* node corresponds to product of binary variables (maybe with coefficient and constant, if polynomial) */
   SCIPdebugMsg(scip, "reformulate node %p via AND-constraint\n", (void*)node);

   /* collect variables in product */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nchildren) );
   for( c = 0; c < nchildren; ++c )
   {
      child = SCIPexprgraphGetNodeChildren(node)[c];
      vars[c] = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, child);
   }

   /* create variable for resultant
    * cons_and wants to add implications for resultant, which is only possible for binary variables currently
    * so choose binary as vartype, even though implicit integer had been sufficient
    */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%dand", *naddcons);
   SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
      TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, var) );

#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      SCIP_Bool debugval;
      SCIP_Real varval;

      debugval = TRUE;
      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( SCIPdebugGetSolVal(scip, vars[c], &varval) );
         debugval = debugval && (varval > 0.5);
      }
      SCIP_CALL( SCIPdebugAddSolVal(scip, var, debugval ? 1.0 : 0.0) );
   }
#endif

   /* create AND-constraint */
   SCIP_CALL( SCIPcreateConsAnd(scip, &cons, name, var, nchildren, vars,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   ++*naddcons;

   SCIPfreeBufferArray(scip, &vars);

   /* add var to exprgraph */
   SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&var, reformnode) );
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   /* if we have coefficient and constant, then replace reformnode by linear expression in reformnode */
   if( SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_POLYNOMIAL )
   {
      SCIP_Real coef;
      SCIP_Real constant;

      coef = SCIPexprGetMonomialCoef(SCIPexprgraphGetNodePolynomialMonomials(node)[0]);
      constant = SCIPexprgraphGetNodePolynomialConstant(node);

      if( coef != 1.0 || constant != 0.0 )
      {
         SCIP_EXPRGRAPHNODE* linnode;
         SCIP_CALL( SCIPexprgraphCreateNodeLinear(SCIPblkmem(scip), &linnode, 1, &coef, constant) );
         SCIP_CALL( SCIPexprgraphAddNode(exprgraph, linnode, -1, 1, reformnode) );
         *reformnode = linnode;
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyAnd)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrAnd(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( nconss == 0 || conss != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->linearize )
   {
      /* linearize all AND-constraints and remove the AND-constraints */
      SCIP_CONS* newcons;
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      char consname[SCIP_MAXSTRLEN];

      SCIP_VAR** vars;
      SCIP_Real* vals;

      int nvars;
      int c, v;

      /* allocate buffer array */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

      for( c = 0; c < nconss; ++c )
      {
         cons = conss[c];
         assert( cons != NULL );

	 /* only added constraints can be upgraded */
	 if( !SCIPconsIsAdded(cons) )
	    continue;

         consdata = SCIPconsGetData(cons);
         assert( consdata != NULL );
         assert( consdata->resvar != NULL );

         nvars = consdata->nvars;

         if( !conshdlrdata->aggrlinearization )
         {
            vars[0] = consdata->resvar;
            vals[0] = 1.0;
            vals[1] = -1.0;

            /* create operator linear constraints */
            for( v = 0; v < nvars; ++v )
            {
               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), v);
               vars[1] = consdata->vars[v];

               SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, 2, vars, vals, -SCIPinfinity(scip), 0.0,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     consdata->checkwhenupgr || SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                     SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
                     !(consdata->notremovablewhenupgr) && SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );


               /* add constraint */
               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            }
         }

         /* reallocate buffer array */
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, nvars + 1) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &vals, nvars + 1) );

         for( v = 0; v < nvars; ++v )
         {
            vars[v] = consdata->vars[v];
            vals[v] = -1.0;
         }

         vars[nvars] = consdata->resvar;

         if( conshdlrdata->aggrlinearization )
         {
            /* create additional linear constraint */
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_operators", SCIPconsGetName(cons));

            vals[nvars] = (SCIP_Real) nvars;

            SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, nvars + 1, vars, vals, -SCIPinfinity(scip), 0.0,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  consdata->checkwhenupgr || SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
                  !(consdata->notremovablewhenupgr) && SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

            /* add constraint */
            SCIP_CALL( SCIPaddCons(scip, newcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         }

         /* create additional linear constraint */
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_add", SCIPconsGetName(cons));

         vals[nvars] = 1.0;

         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, nvars + 1, vars, vals, -nvars + 1.0, SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               consdata->checkwhenupgr || SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
               !(consdata->notremovablewhenupgr) && SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

         /* add constraint */
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

         /* delete constraint */
         SCIP_CALL( SCIPdelCons(scip, cons) );
      }

      /* free buffer array */
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &vals);
   }

   return SCIP_OKAY;
}


#ifdef GMLGATEPRINTING

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreAnd)
{  /*lint --e{715}*/
   SCIP_HASHMAP* hashmap;
   FILE* gmlfile;
   char fname[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** activeconsvars;
   SCIP_VAR* activevar;
   int* varnodeids;
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int ncontvars;
   int v;
   int c;
   unsigned int resid;
   unsigned int varid;
   unsigned int id = 1;

   /* no AND-constraints available */
   if( nconss == 0 )
      return SCIP_OKAY;

   nvars = SCIPgetNVars(scip);

   /* no variables left anymore */
   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varnodeids, nvars) );
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, &ncontvars) );

   /* open gml file */
   (void) SCIPsnprintf(fname, SCIP_MAXSTRLEN, "and-gates%p.gml", scip);
   gmlfile = fopen(fname, "w");

   if( gmlfile == NULL )
   {
      SCIPerrorMessage("cannot open graph file <%s>\n", fname);
      SCIPABORT();
      return SCIP_WRITEERROR;   /*lint !e527*/
   }

   /* create the variable mapping hash map */
   SCIP_CALL_FINALLY( SCIPhashmapCreate(&hashmap, SCIPblkmem(scip), nvars), fclose(gmlfile) );

   /* write starting of gml file */
   SCIPgmlWriteOpening(gmlfile, TRUE);

   /* walk over all AND-constraints */
   for( c = nconss - 1; c >= 0; --c )
   {
      cons = conss[c];

      /* only handle active constraints */
      if( !SCIPconsIsActive(cons) )
         continue;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* only handle constraints which have operands */
      if( consdata->nvars == 0 )
         continue;

      assert(consdata->vars != NULL);
      assert(consdata->resvar != NULL);

      /* get active variable of resultant */
      activevar = SCIPvarGetProbvar(consdata->resvar);

      /* check if we already found this variables */
      resid = (unsigned int)(size_t) SCIPhashmapGetImage(hashmap, activevar);
      if( resid == 0 )
      {
         resid = id;
         ++id;
         SCIP_CALL( SCIPhashmapInsert(hashmap, (void*)activevar, (void*)(size_t)resid) );

         /* write new gml node for new resultant */
         SCIPgmlWriteNode(gmlfile, resid, SCIPvarGetName(activevar), NULL, NULL, NULL);
      }

      /* copy operands to get problem variables for */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activeconsvars, consdata->vars, consdata->nvars) );

      /* get problem variables of operands */
      SCIPvarsGetProbvar(activeconsvars, consdata->nvars);

      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         /* check if we already found this variables */
         varid = (unsigned int)(size_t) SCIPhashmapGetImage(hashmap, activeconsvars[v]);
         if( varid == 0 )
         {
            varid = id;
            ++id;
            SCIP_CALL( SCIPhashmapInsert(hashmap, (void*)activeconsvars[v], (void*)(size_t)varid) );

            /* write new gml node for new operand */
            SCIPgmlWriteNode(gmlfile, varid, SCIPvarGetName(activeconsvars[v]), NULL, NULL, NULL);
         }
         /* write gml arc between resultant and operand */
         SCIPgmlWriteArc(gmlfile, resid, varid, NULL, NULL);
      }

      /* free temporary memory for active constraint variables */
      SCIPfreeBufferArray(scip, &activeconsvars);
   }

   /* write all remaining variables as nodes */
#ifdef SCIP_DISABLED_CODE
   for( v = nvars - 1; v >= 0; --v )
   {
      activevar = SCIPvarGetProbvar(vars[v]);

      varid = (unsigned int)(size_t) SCIPhashmapGetImage(hashmap, activevar);
      if( varid == 0 )
      {
	 varid = id;
	 ++id;
	 SCIP_CALL( SCIPhashmapInsert(hashmap, (void*)activeconsvars[v], (void*)(size_t)varid) );

	 /* write new gml node for new operand */
	 SCIPgmlWriteNode(gmlfile, varid, SCIPvarGetName(activevar), NULL, NULL, NULL);
      }
   }
#endif

   /* free the variable mapping hash map */
   SCIPhashmapFree(&hashmap);

   SCIPgmlWriteClosing(gmlfile);

   fclose(gmlfile);

   SCIPfreeBufferArray(scip, &varnodeids);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}
#endif

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolAnd)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   /* release and free the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      SCIP_CALL( consdataFreeRows(scip, consdata) );
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransAnd)
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
         sourcedata->nvars, sourcedata->vars, sourcedata->resvar, sourcedata->checkwhenupgr,
         sourcedata->notremovablewhenupgr) );

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
SCIP_DECL_CONSINITLP(consInitlpAnd)
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
SCIP_DECL_CONSSEPALP(consSepalpAnd)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   SCIP_Bool cutoff;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &separated, &cutoff) );
      if ( cutoff )
         *result = SCIP_CUTOFF;
      else if ( separated )
         *result = SCIP_SEPARATED;
   }

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolAnd)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   SCIP_Bool cutoff;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &separated, &cutoff) );
      if ( cutoff )
         *result = SCIP_CUTOFF;
      else if ( separated )
         *result = SCIP_SEPARATED;
   }

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpAnd)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, NULL, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxAnd)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, sol, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsAnd)
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
SCIP_DECL_CONSCHECK(consCheckAnd)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   *result = SCIP_FEASIBLE;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss && (*result == SCIP_FEASIBLE || completely); i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, checklprows, printreason, &violated) );
      if( violated )
         *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int nfixedvars;
   int nupgdconss;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   nfixedvars = 0;
   nupgdconss = 0;

   /* propagate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &nfixedvars, &nupgdconss) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 || nupgdconss > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   unsigned char* entries;
   SCIP_Bool cutoff;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int firstchange;
   int nentries;
   int c;

   assert(result != NULL);

   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;

   nentries = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &entries, nentries) );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   cutoff = FALSE;
   firstchange = INT_MAX;
   for( c = 0; c < nconss && !cutoff && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->propagated = FALSE;

      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars, nupgdconss) );

      /* remove all variables that are fixed to one; merge multiple entries of the same variable;
       * fix resultant to zero if a pair of negated variables is contained in the operand variables
       */
      if( !cutoff && !SCIPconsIsDeleted(cons) )
      {
         SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr, nchgcoefs) );

         /* merge multiple occurances of variables or variables with their negated variables */
         SCIP_CALL( mergeMultiples(scip, cons, conshdlrdata->eventhdlr, &entries, &nentries, nfixedvars, nchgcoefs, ndelconss) );
      }

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 1); /* otherwise, propagateCons() has deleted the constraint */

         /* if only one variable is left, the resultant has to be equal to this single variable */
         if( consdata->nvars == 1 )
         {
            SCIP_Bool redundant;
            SCIP_Bool aggregated;

            SCIPdebugMsg(scip, "AND-constraint <%s> has only one variable not fixed to 1.0\n", SCIPconsGetName(cons));

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

            /* add implications: resultant == 1 -> all operands == 1 */
            for( i = 0; i < consdata->nvars && !cutoff; ++i )
            {
               int nimplbdchgs;

               SCIP_CALL( SCIPaddVarImplication(scip, consdata->resvar, TRUE, consdata->vars[i],
                     SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff, &nimplbdchgs) );
               (*nchgbds) += nimplbdchgs;
            }
            consdata->impladded = TRUE;
         }

         /* if in r = x and y, the resultant is fixed to zero, add implication x = 1 -> y = 0 */
         if( !cutoff && SCIPconsIsActive(cons) && consdata->nvars == 2 && !consdata->opimpladded
            && SCIPvarGetUbGlobal(consdata->resvar) < 0.5 )
         {
            int nimplbdchgs;

            SCIP_CALL( SCIPaddVarImplication(scip, consdata->vars[0], TRUE, consdata->vars[1],
                  SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff, &nimplbdchgs) );
            (*nchgbds) += nimplbdchgs;
            consdata->opimpladded = TRUE;
         }
      }
   }

   /* perform dual presolving on AND-constraints */
   if( conshdlrdata->dualpresolving && !cutoff && !SCIPisStopped(scip) && SCIPallowDualReds(scip))
   {
      SCIP_CALL( dualPresolve(scip, conss, nconss, conshdlrdata->eventhdlr, &entries, &nentries, &cutoff, nfixedvars, naggrvars, nchgcoefs, ndelconss, nupgdconss, naddconss) );
   }

   /* check for cliques inside the AND constraint */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
      {
	 if( SCIPconsIsActive(conss[c]) )
	 {
	    /* check if at least two operands are in one clique */
	    SCIP_CALL( cliquePresolve(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, nfixedvars, naggrvars, nchgcoefs, ndelconss, naddconss) );
	 }
      }
   }

   /* process pairs of constraints: check them for equal operands in order to aggregate resultants;
    * only apply this expensive procedure, if the single constraint preprocessing did not find any reductions
    * (otherwise, we delay the presolving to be called again next time)
    */
   if( !cutoff && conshdlrdata->presolusehashing && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      if( *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars )
      {
         if( firstchange < nconss ) 
         {
            /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
            SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, &cutoff, naggrvars, ndelconss) );
            oldnaggrvars = *naggrvars;
         }
      }
   }

   if( !cutoff && conshdlrdata->presolpairwise && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      if( *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars )
      {
         SCIP_Longint npaircomparisons;
         npaircomparisons = 0;
         oldndelconss = *ndelconss;

         for( c = firstchange; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
         {
            if( SCIPconsIsActive(conss[c]) && !SCIPconsIsModifiable(conss[c]) )
            {
               npaircomparisons += ((SCIPconsGetData(conss[c])->changed) ? (SCIP_Longint) c : ((SCIP_Longint) c - (SCIP_Longint) firstchange));

               SCIP_CALL( preprocessConstraintPairs(scip, conss, firstchange, c, &cutoff, naggrvars, nchgbds,
                     ndelconss) );

               if( npaircomparisons > NMINCOMPARISONS )
               {
                  if( ((*ndelconss - oldndelconss) + (*naggrvars - oldnaggrvars) + (*nchgbds - oldnchgbds)/2.0) / ((SCIP_Real) npaircomparisons) < MINGAINPERNMINCOMPARISONS )
                     break;
                  oldndelconss = *ndelconss;
                  oldnaggrvars = *naggrvars;
                  oldnchgbds = *nchgbds;

                  npaircomparisons = 0;
               }
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &entries);

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds
            || *ndelconss > oldndelconss || *nupgdconss > oldnupgdconss )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropAnd)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockAnd)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* resultant variable */
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->resvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );

   /* operand variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintAnd)
{  /*lint --e{715}*/

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyAnd)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_VAR** vars;
   SCIP_VAR* sourceresvar;
   SCIP_VAR* resvar;
   const char* consname;
   int nvars;
   int v;

   assert(valid != NULL);
   (*valid) = TRUE;

   sourceresvar = SCIPgetResultantAnd(sourcescip, sourcecons);

   /* map resultant to active variable of the target SCIP  */
   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourceresvar, &resvar, varmap, consmap, global, valid) );
   assert(!(*valid) || resvar != NULL);

   /* we do not copy, if a variable is missing */
   if( !(*valid) )
      return SCIP_OKAY;

   /* map operand variables to active variables of the target SCIP  */
   sourcevars = SCIPgetVarsAnd(sourcescip, sourcecons);
   nvars = SCIPgetNVarsAnd(sourcescip, sourcecons);

   if( nvars == -1 )
      return SCIP_INVALIDCALL;

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &vars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || vars[v] != NULL);

      /* we do not copy, if a variable is missing */
      if( !(*valid) )
         goto TERMINATE;
   }

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   /* creates and captures a AND-constraint */
   SCIP_CALL( SCIPcreateConsAnd(scip, cons, consname, resvar, nvars, vars, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

 TERMINATE:   
   /* free buffer array */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseAnd)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR* resvar;
   char* endptr;
   int requiredsize;
   int varssize;
   int nvars;

   SCIPdebugMsg(scip, "parse <%s> as AND-constraint\n", str);

   *success = FALSE;

   /* parse variable name of resultant */
   SCIP_CALL( SCIPparseVarName(scip, str, &resvar, &endptr) );
   str = endptr;

   if( resvar == NULL )
   {
      SCIPdebugMsg(scip, "resultant variable does not exist \n");
   }
   else
   {
      char* strcopy = NULL;
      char* startptr;

      /* cutoff "== and(" form the constraint string */
      startptr = strchr((char*)str, '(');

      if( startptr == NULL )
      {
         SCIPerrorMessage("missing starting character '(' parsing AND-constraint\n");
         return SCIP_OKAY;
      }

      /* skip '(' */
      ++startptr;

      /* find end character ')' */
      endptr = strrchr(startptr, ')');

      if( endptr == NULL )
      {
         SCIPerrorMessage("missing ending character ')' parsing AND-constraint\n");
         return SCIP_OKAY;
      }
      assert(endptr >= startptr);

      if( endptr > startptr )
      {
         /* copy string for parsing */
         SCIP_CALL( SCIPduplicateBufferArray(scip, &strcopy, startptr, (int)(endptr-startptr)) );

         varssize = 100;
         nvars = 0;

         /* allocate buffer array for variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );

         /* parse string */
         SCIP_CALL( SCIPparseVarsList(scip, strcopy, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );

         if( *success )
         {
            /* check if the size of the variable array was great enough */
            if( varssize < requiredsize )
            {
               /* reallocate memory */
               varssize = requiredsize;
               SCIP_CALL( SCIPreallocBufferArray(scip, &vars, varssize) );

               /* parse string again with the correct size of the variable array */
               SCIP_CALL( SCIPparseVarsList(scip, strcopy, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );
            }

            assert(*success);
            assert(varssize >= requiredsize);

            /* create AND-constraint */
            SCIP_CALL( SCIPcreateConsAnd(scip, cons, name, resvar, nvars, vars,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
         }

         /* free variable buffer */
         SCIPfreeBufferArray(scip, &vars);
         SCIPfreeBufferArray(scip, &strcopy);
      }
      else
      {
         if( !modifiable )
         {
            SCIPerrorMessage("cannot create empty AND-constraint\n");
            return SCIP_OKAY;
         }

         /* create empty AND-constraint */
         SCIP_CALL( SCIPcreateConsAnd(scip, cons, name, resvar, 0, NULL,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

         *success = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsAnd)
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
SCIP_DECL_CONSGETNVARS(consGetNVarsAnd)
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
SCIP_DECL_EVENTEXEC(eventExecAnd)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   /* check, if the variable was fixed to zero */
   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_UBTIGHTENED )
      consdata->nofixedzero = FALSE;

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for AND-constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrAnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;

   /* create event handler for events on variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecAnd, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpAnd, consEnfopsAnd, consCheckAnd, consLockAnd,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyAnd, consCopyAnd) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteAnd) );
#ifdef GMLGATEPRINTING
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreAnd) );
#endif
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolAnd) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeAnd) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsAnd) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsAnd) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreAnd) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpAnd) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseAnd) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolAnd, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintAnd) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropAnd, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropAnd) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpAnd, consSepasolAnd, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransAnd) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxAnd) );

   /* add AND-constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/and/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance",
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/linearize",
         "should the AND-constraint get linearized and removed (in presolving)?",
         &conshdlrdata->linearize, TRUE, DEFAULT_LINEARIZE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/enforcecuts",
         "should cuts be separated during LP enforcing?",
         &conshdlrdata->enforcecuts, TRUE, DEFAULT_ENFORCECUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/aggrlinearization",
         "should an aggregated linearization be used?",
         &conshdlrdata->aggrlinearization, TRUE, DEFAULT_AGGRLINEARIZATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/upgraderesultant",
         "should all binary resultant variables be upgraded to implicit binary variables?",
         &conshdlrdata->upgrresultant, TRUE, DEFAULT_UPGRRESULTANT, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/dualpresolving",
         "should dual presolving be performed?",
         &conshdlrdata->dualpresolving, TRUE, DEFAULT_DUALPRESOLVING, NULL, NULL) );

   if( SCIPfindConshdlr(scip, "nonlinear") != NULL )
   {
      /* include the AND-constraint upgrade in the nonlinear constraint handler */
      SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, NULL, exprgraphnodeReformAnd, EXPRGRAPHREFORM_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   return SCIP_OKAY;
}

/** creates and captures a AND-constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsAnd(
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
   SCIP_Bool infeasible;

   /* find the AND-constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("AND-constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* upgrade binary resultant variable to an implicit binary variable */
   /* @todo add implicit upgrade in presolving, improve decision making for upgrade by creating an implication graph */
   if( conshdlrdata->upgrresultant && SCIPvarGetType(resvar) == SCIP_VARTYPE_BINARY
#if 1 /* todo delete following hack,
       *      the following avoids upgrading not artificial variables, for example and-resultants which are generated
       *      from the gate presolver, it seems better to not upgrade these variables
       */
      && strlen(SCIPvarGetName(resvar)) > strlen(ARTIFICIALVARNAMEPREFIX) && strncmp(SCIPvarGetName(resvar), ARTIFICIALVARNAMEPREFIX, strlen(ARTIFICIALVARNAMEPREFIX)) == 0 )
#else
      )
#endif
   {
      SCIP_VAR* activeresvar;
      SCIP_VAR* activevar;
      int v;

      if( SCIPisTransformed(scip) )
         activeresvar = SCIPvarGetProbvar(resvar);
      else
         activeresvar = resvar;

      if( SCIPvarGetType(activeresvar) == SCIP_VARTYPE_BINARY )
      {
         /* check if we can upgrade the variable type of the resultant */
         for( v = nvars - 1; v >= 0; --v )
         {
            if( SCIPisTransformed(scip) )
               activevar = SCIPvarGetProbvar(vars[v]);
            else
               activevar = vars[v];

            if( activevar == activeresvar || SCIPvarGetType(activevar) == SCIP_VARTYPE_IMPLINT )
               break;
         }

         /* upgrade the type of the resultant */
         if( v < 0 )
         {
            SCIP_CALL( SCIPchgVarType(scip, resvar, SCIP_VARTYPE_IMPLINT, &infeasible) );
            assert(!infeasible);
         }
      }
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, resvar, FALSE, FALSE) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures an AND-constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsAnd(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsAnd() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             resvar,             /**< resultant variable of the operation */
   int                   nvars,              /**< number of operator variables in the constraint */
   SCIP_VAR**            vars                /**< array with operator variables of constraint */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsAnd(scip, cons, name, resvar, nvars, vars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


/** gets number of variables in AND-constraint */
int SCIPgetNVarsAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an AND-constraint\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in AND-constraint */
SCIP_VAR** SCIPgetVarsAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an AND-constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}


/** gets the resultant variable in AND-constraint */
SCIP_VAR* SCIPgetResultantAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an AND-constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->resvar;
}

/** return if the variables of the AND-constraint are sorted with respect to their indices */
SCIP_Bool SCIPisAndConsSorted(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an AND-constraint\n");
      SCIPABORT();
      return FALSE;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->sorted;
}

/** sort the variables of the AND-constraint with respect to their indices */
SCIP_RETCODE SCIPsortAndCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an AND-constraint\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdataSort(consdata);
   assert(consdata->sorted);

   return SCIP_OKAY;
}

/** when 'upgrading' the given AND-constraint, should the check flag for the upgraded constraint be set to TRUE, even if
 *  the check flag of this AND-constraint is set to FALSE?
 */
SCIP_RETCODE SCIPchgAndConsCheckFlagWhenUpgr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Bool             flag                /**< should an arising constraint from the given AND-constraint be checked,
                                              *   even if the check flag of the AND-constraint is set to FALSE
                                              */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an AND-constraint\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->checkwhenupgr = flag;

   return SCIP_OKAY;
}

/** when 'upgrading' the given AND-constraint, should the removable flag for the upgraded constraint be set to FALSE,
 *  even if the removable flag of this AND-constraint is set to TRUE?
 */
SCIP_RETCODE SCIPchgAndConsRemovableFlagWhenUpgr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Bool             flag                /**< should an arising constraint from the given AND-constraint be not
                                              *   removable, even if the removable flag of the AND-constraint is set to
                                              *   TRUE
                                              */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an AND-constraint\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->notremovablewhenupgr = flag;

   return SCIP_OKAY;
}
