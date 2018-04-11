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

/**@file   cons_xor.c
 * @brief  Constraint handler for "xor" constraints,  \f$rhs = x_1 \oplus x_2 \oplus \dots  \oplus x_n\f$
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author Michael Winkler
 *
 * This constraint handler deals with "xor" constraint. These are constraint of the form:
 *
 * \f[
 *    rhs = x_1 \oplus x_2 \oplus \dots  \oplus x_n
 * \f]
 *
 * where \f$x_i\f$ is a binary variable for all \f$i\f$ and \f$rhs\f$ is bool. The variables \f$x\f$'s are called
 * operators. This constraint is satisfied if \f$rhs\f$ is TRUE and an odd number of the operators are TRUE or if the
 * \f$rhs\f$ is FALSE and a even number of operators are TRUE. Hence, if the sum of \f$rhs\f$ and operators is even.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/pub_misc.h"
#include "scip/cons_xor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_linear.h"
#include "scip/heur_trysol.h"
#include "scip/debug.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "xor"
#define CONSHDLR_DESC          "constraint handler for xor constraints: r = xor(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850200 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -850200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -850200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_ALWAYS

#define EVENTHDLR_NAME         "xor"
#define EVENTHDLR_DESC         "event handler for xor constraints"

#define LINCONSUPGD_PRIORITY    +600000 /**< priority of the constraint handler for upgrading of linear constraints */

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define DEFAULT_ADDEXTENDEDFORM   FALSE /**< should the extended formulation be added in presolving? */
#define DEFAULT_ADDFLOWEXTENDED   FALSE /**< should the extended flow formulation be added (nonsymmetric formulation otherwise)? */
#define DEFAULT_SEPARATEPARITY    FALSE /**< should parity inequalities be separated? */
#define DEFAULT_GAUSSPROPFREQ         5 /**< frequency for applying the Gauss propagator */
#define HASHSIZE_XORCONS            500 /**< minimal size of hash table in logicor constraint tables */
#define DEFAULT_PRESOLUSEHASHING   TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define NMINCOMPARISONS          200000 /**< number for minimal pairwise presolving comparisons */
#define MINGAINPERNMINCOMPARISONS 1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */
#define MAXXORCONSSSYSTEM          1000 /**< maximal number of active constraints for which checking the system over GF2 is performed */
#define MAXXORVARSSYSTEM           1000 /**< maximal number of variables in xor constraints for which checking the system over GF2 is performed */

#define NROWS 5


/*
 * Data structures
 */

/** type used for matrix entries in function checkGauss() */
typedef unsigned short Type;

/** constraint data for xor constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables in the xor operation */
   SCIP_VAR*             intvar;             /**< internal variable for LP relaxation */
   SCIP_VAR**            extvars;            /**< variables in extended (flow|asymmetric) formulation (order for flow formulation: nn, ns, sn, ss) */
   SCIP_ROW*             rows[NROWS];        /**< rows for linear relaxation of xor constraint */
   int                   nvars;              /**< number of variables in xor operation */
   int                   nextvars;           /**< number of variables in extended flow formulation */
   int                   varssize;           /**< size of vars array */
   int                   extvarssize;        /**< size of extvars array */
   int                   watchedvar1;        /**< position of first watched operator variable */
   int                   watchedvar2;        /**< position of second watched operator variable */
   int                   filterpos1;         /**< event filter position of first watched operator variable */
   int                   filterpos2;         /**< event filter position of second watched operator variable */
   SCIP_Bool             rhs;                /**< right hand side of the constraint */
   unsigned int          deleteintvar:1;     /**< should artificial variable be deleted */
   unsigned int          propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          changed:1;          /**< was constraint changed since last pair preprocessing round? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for events on watched variables */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance? */
   SCIP_Bool             addextendedform;    /**< should the extended formulation be added in presolving? */
   SCIP_Bool             addflowextended;    /**< should the extended flow formulation be added (nonsymmetric formulation otherwise)? */
   SCIP_Bool             separateparity;     /**< should parity inequalities be separated? */
   int                   gausspropfreq;      /**< frequency for applying the Gauss propagator */
};


/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_0,                          /**< all variables are fixed => fix integral variable */
   PROPRULE_1,                          /**< all except one variable fixed  =>  fix remaining variable */
   PROPRULE_INTLB,                      /**< lower bound propagation of integral variable */
   PROPRULE_INTUB,                      /**< upper bound propagation of integral variable */
   PROPRULE_INVALID                     /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;


/*
 * Local methods
 */

/** installs rounding locks for the given variable in the given xor constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given xor constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
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

/** stores the given variable numbers as watched variables, and updates the event processing */
static
SCIP_RETCODE consdataSwitchWatchedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< xor constraint data */
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
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar1], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar2], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[watchedvar1], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[watchedvar2], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, &consdata->filterpos2) );
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

/** creates constraint data for xor constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   SCIP_Bool             rhs,                /**< right hand side of the constraint */
   int                   nvars,              /**< number of variables in the xor operation */
   SCIP_VAR**            vars,               /**< variables in xor operation */
   SCIP_VAR*             intvar              /**< artificial integer variable for linear relaxation */
   )
{
   int r;

   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );

   (*consdata)->rhs = rhs;
   (*consdata)->intvar = intvar;
   for( r = 0; r < NROWS; ++r )
      (*consdata)->rows[r] = NULL;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;
   (*consdata)->deleteintvar = (intvar == NULL);
   (*consdata)->propagated = FALSE;
   (*consdata)->sorted = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->extvars = NULL;
   (*consdata)->nextvars = 0;
   (*consdata)->extvarssize = 0;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

      if( (*consdata)->intvar != NULL )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->intvar, &((*consdata)->intvar)) );
      }

      if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      {
         SCIP_CONSHDLRDATA* conshdlrdata;
         SCIP_CONSHDLR* conshdlr;
         int v;

         conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
         assert(conshdlr != NULL);
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);

         for( v = (*consdata)->nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, (*consdata)->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
                  (SCIP_EVENTDATA*)(*consdata), NULL) );
         }
      }
   }

   if( (*consdata)->intvar != NULL )
   {
      /* capture artificial variable */
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->intvar) );
   }

   return SCIP_OKAY;
}

/** releases LP row of constraint data */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int r;

   assert(consdata != NULL);

   for( r = 0; r < NROWS; ++r )
   {
      if( consdata->rows[r] != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
   }

   return SCIP_OKAY;
}

/** frees constraint data for xor constraint */
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
      int j;

      /* drop events for watched variables */
      SCIP_CALL( consdataSwitchWatchedvars(scip, *consdata, eventhdlr, -1, -1) );

      /* release flow variables */
      if ( (*consdata)->nextvars > 0 )
      {
         assert( (*consdata)->extvars != NULL );
         for (j = 0; j < (*consdata)->extvarssize; ++j)
         {
            if ( (*consdata)->extvars[j] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->extvars[j])) );
            }
         }

         SCIPfreeBlockMemoryArray(scip, &((*consdata)->extvars), (*consdata)->extvarssize);
         (*consdata)->nextvars = 0;
         (*consdata)->extvarssize = 0;
      }
   }
   else
   {
      assert((*consdata)->watchedvar1 == -1);
      assert((*consdata)->watchedvar2 == -1);
   }

   /* release LP row */
   SCIP_CALL( consdataFreeRows(scip, *consdata) );

   /* release internal variable */
   if( (*consdata)->intvar != NULL )
   {
      /* if the constraint is deleted and the integral variable is present, it should be fixed */
      assert( SCIPisEQ(scip, SCIPvarGetLbGlobal((*consdata)->intvar), SCIPvarGetLbGlobal((*consdata)->intvar)) );

      /* We do not delete the integral variable, but leave the handling to SCIP, because it might happen that the
         integral variable is stored in some basis information somewhere. */
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->intvar) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints xor constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< xor constraint data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             endline             /**< should an endline be set? */
   )
{
   assert(consdata != NULL);

   /* start variable list */
   SCIPinfoMessage(scip, file, "xor(");

   /* print variable list */
   SCIP_CALL( SCIPwriteVarsList(scip, file, consdata->vars, consdata->nvars, TRUE, ',') );

   /* close variable list and write right hand side */
   SCIPinfoMessage(scip, file, ") = %d", consdata->rhs);

   /* write integer variable if it exists */
   if( consdata->intvar != NULL )
   {
      SCIPinfoMessage(scip, file, " (intvar = ");
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->intvar, TRUE) );
      SCIPinfoMessage(scip, file, ")");
   }

   if( endline )
      SCIPinfoMessage(scip, file, "\n");

   return SCIP_OKAY;
}

/** sets intvar of an xor constraint */
static
SCIP_RETCODE setIntvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows[0] == NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   /* remove the rounding locks for the old variable and release it */
   if( consdata->intvar != NULL )
   {
      SCIP_CALL( unlockRounding(scip, cons, consdata->intvar) );
      SCIP_CALL( SCIPreleaseVar(scip, &(consdata->intvar)) );
   }

   consdata->intvar = var;
   consdata->changed = TRUE;

   /* install the rounding locks for the new variable and capture it */
   SCIP_CALL( lockRounding(scip, cons, consdata->intvar) );
   SCIP_CALL( SCIPcaptureVar(scip, consdata->intvar) );

   /**@todo update LP rows */
   if( consdata->rows[0] != NULL )
   {
      SCIPerrorMessage("cannot change intvar of xor constraint after LP relaxation was created\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

/** adds coefficient to xor constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows[0] == NULL);

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

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockRounding(scip, cons, var) );

   /* we only catch this event in presolving stages
    * we need to catch this event also during exiting presolving because we call applyFixings to clean up the constraint
    * and this can lead to an insertion of a replacement of variables for which we will try to drop the VARFIXED event.
    */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_INITPRESOLVE
         || SCIPgetStage(scip) == SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_CONSHDLR* conshdlr;

      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
      assert(conshdlr != NULL);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
            (SCIP_EVENTDATA*)consdata, NULL) );
   }

   /**@todo update LP rows */
   if( consdata->rows[0] != NULL )
   {
      SCIPerrorMessage("cannot add coefficients to xor constraint after LP relaxation was created\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from xor constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
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

   /* remove the rounding locks of the deleted variable */
   SCIP_CALL( unlockRounding(scip, cons, consdata->vars[pos]) );

   /* we only catch this event in presolving stage, so we need to only drop it there */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_INITPRESOLVE
         || SCIPgetStage(scip) == SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_VARFIXED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, -1) );
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
   consdata->sorted = FALSE;
   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** sorts and constraint's variables by non-decreasing variable index */
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
         SCIPsortPtr((void**)consdata->vars, SCIPvarCompActiveAndNegated, consdata->nvars);
         consdata->sorted = TRUE;

         /* correct watched variables */
         if( var1 != NULL )
         {
            int v;

            /* since negated variables exist, we need to loop over all variables to find the old variable and cannot use
             * SCIPsortedvecFindPtr()
             */
            for( v = consdata->nvars - 1; v >= 0; --v )
            {
               if( consdata->vars[v] == var1 )
               {
                  consdata->watchedvar1 = v;
                  if( var2 == NULL || consdata->watchedvar2 != -1 )
                     break;
               }
               else if( consdata->vars[v] == var2 )
               {
                  assert(consdata->vars[v] != NULL);
                  consdata->watchedvar2 = v;
                  if( consdata->watchedvar1 != -1 )
                     break;
               }
            }
            assert(consdata->watchedvar1 != -1);
            assert(consdata->watchedvar2 != -1 || var2 == NULL);
            assert(consdata->watchedvar1 < consdata->nvars);
            assert(consdata->watchedvar2 < consdata->nvars);
         }
      }
   }

#ifdef SCIP_DEBUG
   /* check sorting */
   {
      int v;

      for( v = 0; v < consdata->nvars; ++v )
      {
         assert(v == consdata->nvars-1 || SCIPvarCompareActiveAndNegated(consdata->vars[v], consdata->vars[v+1]) <= 0);
      }
   }
#endif
}


/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyXorcons)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqXorcons)
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
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

   for( i = 0; i < consdata1->nvars ; ++i )
   {
      /* tests if variables are equal */
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 ||
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         return FALSE;
      }
      assert(SCIPvarCompareActiveAndNegated(consdata1->vars[i], consdata2->vars[i]) == 0);
   }

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValXorcons)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int minidx;
   int mididx;
   int maxidx;

   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->sorted);
   assert(consdata->nvars > 0);

   /* only active, fixed or negated variables are allowed */
   assert(consdata->vars[0] != NULL);
   assert(consdata->vars[consdata->nvars / 2] != NULL);
   assert(consdata->vars[consdata->nvars - 1] != NULL);
   assert(SCIPvarIsActive(consdata->vars[0]) || SCIPvarGetStatus(consdata->vars[0]) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(consdata->vars[0]) == SCIP_VARSTATUS_FIXED);
   assert(SCIPvarIsActive(consdata->vars[consdata->nvars / 2]) || SCIPvarGetStatus(consdata->vars[consdata->nvars / 2]) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(consdata->vars[consdata->nvars / 2]) == SCIP_VARSTATUS_FIXED);
   assert(SCIPvarIsActive(consdata->vars[consdata->nvars - 1]) || SCIPvarGetStatus(consdata->vars[consdata->nvars - 1]) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(consdata->vars[consdata->nvars - 1]) == SCIP_VARSTATUS_FIXED);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   /* note that for all indices it does not hold that they are sorted, because variables are sorted with
    * SCIPvarCompareActiveAndNegated (see var.c)
    */

   return SCIPhashTwo(SCIPcombineTwoInt(consdata->nvars, minidx),
                      SCIPcombineTwoInt(mididx, maxidx));
}

/** deletes all fixed variables and all pairs of equal variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int*                  nchgcoefs,          /**< pointer to add up the number of changed coefficients */
   int*                  naggrvars,          /**< pointer to add up the number of aggregated variables */
   int*                  naddconss,          /**< pointer to add up the number of added constraints */
   SCIP_Bool*            cutoff              /**< whether a cutoff has been detected */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(nchgcoefs != NULL);

   SCIPdebugMsg(scip, "before fixings: ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   v = 0;
   while( v < consdata->nvars )
   {
      SCIP_VAR* var;

      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs)++;
      }
      else if( SCIPvarGetLbGlobal(var) > 0.5 && consdata->deleteintvar )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         consdata->rhs = !consdata->rhs;
         (*nchgcoefs)++;
      }
      else
      {
         SCIP_VAR* repvar;
         SCIP_Bool negated;

         /* get binary representative of variable */
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );

         /* remove all negations by replacing them with the active variable
          * it holds that xor(x1, ~x2) = 0 <=> xor(x1, x2) = 1
          * @note this can only be done if the integer variable does not exist
          */
         if( negated && consdata->intvar == NULL )
         {
            assert(SCIPvarIsNegated(repvar));

            repvar = SCIPvarGetNegationVar(repvar);
            consdata->rhs = !consdata->rhs;
         }

         /* check, if the variable should be replaced with the representative */
         if( repvar != var )
         {
            /* delete old (aggregated) variable */
            SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );

            /* add representative instead */
            SCIP_CALL( addCoef(scip, cons, repvar) );
         }
         else
            ++v;
      }
   }

   /* sort the variables in the constraint */
   consdataSort(consdata);
   assert(consdata->sorted);

   SCIPdebugMsg(scip, "after sort    : ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   /* delete pairs of equal or negated variables; scan from back to front because deletion doesn't affect the
    * order of the front variables
    */
   v = consdata->nvars-2;
   while ( v >= 0 )
   {
      if( consdata->vars[v] == consdata->vars[v+1] ) /*lint !e679*/
      {
         SCIP_VAR* newvars[3];
         SCIP_Real vals[3];

         newvars[2] = consdata->vars[v];
         vals[2] = 1.0;

         /* delete both variables */
         SCIPdebugMsg(scip, "xor constraint <%s>: deleting pair of equal variables <%s>\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[v]));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v+1) );
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs) += 2;
         v = MIN(v, consdata->nvars-1);

         /* need to update integer variable, consider the following case:
          * xor(x1, x2, x3, x4, x5) = 0  (and x1 == x2) was change above to
          * xor(        x3, x4, x5) = 0
          * assuming we have an integer variable y it needs to be replaced by z with y = x1 + z and z in [lb_y, ub_y]
          */
         if( consdata->intvar != NULL )
         {
            SCIP_CONS* newcons;
            SCIP_Real lb;
            SCIP_Real ub;
            SCIP_VARTYPE vartype;
            char varname[SCIP_MAXSTRLEN];
            char consname[SCIP_MAXSTRLEN];

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "agg_%s", SCIPvarGetName(consdata->intvar));
            lb = SCIPvarGetLbGlobal(consdata->intvar);
            ub = SCIPvarGetUbGlobal(consdata->intvar);
            vartype = SCIPvarGetType(consdata->intvar);

            SCIP_CALL( SCIPcreateVar(scip, &newvars[0], varname, lb, ub, 0.0, vartype,
                  SCIPvarIsInitial(consdata->intvar), SCIPvarIsRemovable(consdata->intvar),
                  NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, newvars[0]) );
            vals[0] = 1.0;

            newvars[1] = consdata->intvar;
            vals[1] = -1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "agg_%s", SCIPconsGetName(cons));

            SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, 3, newvars, vals, 0.0, 0.0,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), TRUE, /*SCIPconsIsEnforced(cons),*/
                  TRUE, TRUE, /*SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),*/
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

            SCIP_CALL( SCIPaddCons(scip, newcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            ++(*naddconss);

            SCIP_CALL( setIntvar(scip, cons, newvars[0]) );
            SCIP_CALL( SCIPreleaseVar(scip, &newvars[0]) );
         }
      }
      else if( consdata->vars[v] == SCIPvarGetNegatedVar(consdata->vars[v+1]) ) /*lint !e679*/
      {
         /* delete both variables and negate the rhs */
         SCIPdebugMsg(scip, "xor constraint <%s>: deleting pair of negated variables <%s> and <%s>\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[v]), SCIPvarGetName(consdata->vars[v+1])); /*lint !e679*/
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v+1) );
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs) += 2;
         consdata->rhs = !consdata->rhs;
         v = MIN(v, consdata->nvars-1);

         /* need to update integer variable, consider the following case:
          * xor(x1, x2, x3, x4, x5) = 0  (and x1 = ~x2) was change above to
          * xor(        x3, x4, x5) = 1
          * assuming we have an integer variable y it needs to be replaced by z with y = 1 + z and z in [lb_y, ub_y - 1]
          */
         if( consdata->rhs && consdata->intvar != NULL )
         {
            SCIP_VAR* newvar;
            SCIP_Real lb;
            SCIP_Real ub;
            SCIP_VARTYPE vartype;
            char varname[SCIP_MAXSTRLEN];
            SCIP_Bool aggregated;
            SCIP_Bool infeasible;
            SCIP_Bool redundant;

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "agg_%s", SCIPvarGetName(consdata->intvar));
            ub = SCIPvarGetUbGlobal(consdata->intvar) - 1;
            lb = MIN(ub, SCIPvarGetLbGlobal(consdata->intvar)); /*lint !e666*/
            vartype = (lb == 0 && ub == 1) ? SCIP_VARTYPE_BINARY : SCIPvarGetType(consdata->intvar);

            SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, lb, ub, 0.0, vartype,
                  SCIPvarIsInitial(consdata->intvar), SCIPvarIsRemovable(consdata->intvar),
                  NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, newvar) );

            SCIP_CALL( SCIPaggregateVars(scip, consdata->intvar, newvar, 1.0, -1.0, 1.0, &infeasible, &redundant, &aggregated) );
            assert(infeasible || redundant || SCIPdoNotAggr(scip));

            if( infeasible )
            {
               *cutoff = TRUE;
               break;
            }

            if( aggregated )
            {
               (*naggrvars)++;

               if( SCIPvarIsActive(newvar) )
               {
                  SCIP_CALL( setIntvar(scip, cons, newvar) );
                  SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
               }
               /* the new variable should only by inactive if it was fixed due to the aggregation, so also the old variable
                * should be fixed now.
                *
                * @todo if newvar is not active we may want to transform the xor into a linear constraint
                */
               else
               {
                  assert(SCIPvarGetStatus(newvar) == SCIP_VARSTATUS_FIXED);
                  assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->intvar), SCIPvarGetUbGlobal(consdata->intvar)));

                  SCIP_CALL( setIntvar(scip, cons, newvar) );
                  SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
               }
            }
            else
            {
               SCIP_CONS* newcons;
               char consname[SCIP_MAXSTRLEN];
               SCIP_VAR* newvars[2];
               SCIP_Real vals[2];

               newvars[0] = consdata->intvar;
               vals[0] = 1.0;
               newvars[1] = newvar;
               vals[1] = -1.0;

               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "agg_%s", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, 2, newvars, vals, 1.0, 1.0,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), TRUE, /*SCIPconsIsEnforced(cons),*/
                     TRUE, TRUE, /*SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),*/
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
               ++(*naddconss);

               SCIP_CALL( setIntvar(scip, cons, newvar) );
               SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
            }
         }
      }
      else
         assert(SCIPvarGetProbvar(consdata->vars[v]) != SCIPvarGetProbvar(consdata->vars[v+1])); /*lint !e679*/
      --v;
   }

   SCIPdebugMsg(scip, "after fixings : ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   return SCIP_OKAY;
}

/** adds extended flow formulation
 *
 *  The extended flow formulation is built as follows: Let \f$x_1, \dots, x_k\f$ be the variables contained in the given
 *  XOR constraint. We construct a two layered flow network. The upper layer is called the north layer and the lower is
 *  called the south layer. For each \f$x_i,\; i = 2, \ldots, k-1\f$, we add arcs that stay in the north and south layer
 *  (denoted by 'nn' and 'ss', respectively), as well as arcs that change the layers (denoted by 'ns' and 'sn'). For
 *  \f$x_1\f$, we only add two arcs from the source to the two layers. The source is located on the north layer. For
 *  \f$x_k\f$, we add two arcs connecting the two layers to the sink. Depending on the rhs of the constraint the sink is
 *  located on the north or south layer. A change in the layers corresponds to a parity change, i.e., the corresponding
 *  variable \f$x_i\f$ is 1 (and 0 otherwise).
 */
static
SCIP_RETCODE addExtendedFlowFormulation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   int*                  naddedconss         /**< number of added constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* consdata;
   SCIP_VAR* varprevnn = NULL;
   SCIP_VAR* varprevns = NULL;
   SCIP_VAR* varprevsn = NULL;
   SCIP_VAR* varprevss = NULL;
   SCIP_VAR* vars[4];
   SCIP_Real vals[4];
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( naddedconss != NULL );
   *naddedconss = 0;

   /* exit if contraints is modifiable */
   if ( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* exit if extended formulation has been added already */
   if ( consdata->extvars != NULL )
      return SCIP_OKAY;

   /* xor constraints with at most 3 variables are handled directly through rows for the convex hull */
   if ( consdata->nvars <= 3 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Add extended formulation for xor constraint <%s> ...\n", SCIPconsGetName(cons));
   assert( consdata->extvars == NULL );
   assert( consdata->nextvars == 0 );
   assert( consdata->extvarssize == 0 );

   /* get storage for auxiliary variables */
   consdata->extvarssize = 4 * (consdata->nvars);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->extvars), consdata->extvarssize) );

   /* pass through components */
   for (i = 0; i < consdata->nvars; ++i)
   {
      /* variables: n - north, s - south */
      SCIP_VAR* varnn = NULL;
      SCIP_VAR* varns = NULL;
      SCIP_VAR* varsn = NULL;
      SCIP_VAR* varss = NULL;
      SCIP_CONS* newcons;
      SCIP_Real rhs = 0.0;
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool redundant = FALSE;
      SCIP_Bool aggregated = FALSE;
      int cnt = 0;

      /* create variables */
      if ( i == 0 )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_nn", SCIPconsGetName(cons), i);
         SCIP_CALL( SCIPcreateVar(scip, &varnn, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, varnn) );

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_ns", SCIPconsGetName(cons), i);
         SCIP_CALL( SCIPcreateVar(scip, &varns, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, varns) );

         /* need to lock variables, because we aggregate them */
         SCIP_CALL( SCIPlockVarCons(scip, varnn, cons, TRUE, TRUE) );
         SCIP_CALL( SCIPlockVarCons(scip, varns, cons, TRUE, TRUE) );

         /* aggregate ns variable with original variable */
         SCIP_CALL( SCIPaggregateVars(scip, varns, consdata->vars[0], 1.0, -1.0, 0.0, &infeasible, &redundant, &aggregated) );
         assert( ! infeasible );
         assert( redundant );
         assert( aggregated );
      }
      else
      {
         if ( i == consdata->nvars-1 )
         {
            if ( consdata->rhs )
            {
               /* if the rhs is 1 (true) the flow goes to the bottom level */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_ns", SCIPconsGetName(cons), i);
               SCIP_CALL( SCIPcreateVar(scip, &varns, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, varns) );

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_ss", SCIPconsGetName(cons), i);
               SCIP_CALL( SCIPcreateVar(scip, &varss, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, varss) );

               /* need to lock variables, because we aggregate them */
               SCIP_CALL( SCIPlockVarCons(scip, varns, cons, TRUE, TRUE) );
               SCIP_CALL( SCIPlockVarCons(scip, varss, cons, TRUE, TRUE) );

               /* aggregate ns variable with original variable */
               SCIP_CALL( SCIPaggregateVars(scip, varns, consdata->vars[i], 1.0, -1.0, 0.0, &infeasible, &redundant, &aggregated) );
               assert( ! infeasible );
               assert( redundant );
               assert( aggregated );
            }
            else
            {
               /* if the rhs is 0 (false) the flow stays on the top level */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_nn", SCIPconsGetName(cons), i);
               SCIP_CALL( SCIPcreateVar(scip, &varnn, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, varnn) );

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_sn", SCIPconsGetName(cons), i);
               SCIP_CALL( SCIPcreateVar(scip, &varsn, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, varsn) );

               /* need to lock variables, because we aggregate them */
               SCIP_CALL( SCIPlockVarCons(scip, varnn, cons, TRUE, TRUE) );
               SCIP_CALL( SCIPlockVarCons(scip, varsn, cons, TRUE, TRUE) );

               /* aggregate sn variable with original variable */
               SCIP_CALL( SCIPaggregateVars(scip, varsn, consdata->vars[i], 1.0, -1.0, 0.0, &infeasible, &redundant, &aggregated) );
               assert( ! infeasible );
               assert( redundant );
               assert( aggregated );
            }
         }
         else
         {
            /* add the four flow variables */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_nn", SCIPconsGetName(cons), i);
            SCIP_CALL( SCIPcreateVar(scip, &varnn, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, varnn) );

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_ns", SCIPconsGetName(cons), i);
            SCIP_CALL( SCIPcreateVar(scip, &varns, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, varns) );

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_sn", SCIPconsGetName(cons), i);
            SCIP_CALL( SCIPcreateVar(scip, &varsn, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, varsn) );

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_ss", SCIPconsGetName(cons), i);
            SCIP_CALL( SCIPcreateVar(scip, &varss, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, varss) );

            SCIP_CALL( SCIPlockVarCons(scip, varnn, cons, TRUE, TRUE) );
            SCIP_CALL( SCIPlockVarCons(scip, varns, cons, TRUE, TRUE) );
            SCIP_CALL( SCIPlockVarCons(scip, varsn, cons, TRUE, TRUE) );
            SCIP_CALL( SCIPlockVarCons(scip, varss, cons, TRUE, TRUE) );

            /* add coupling constraint */
            cnt = 0;
            if ( varns != NULL )
            {
               vars[cnt] = varns;
               vals[cnt++] = 1.0;
            }
            if ( varsn != NULL )
            {
               vars[cnt] = varsn;
               vals[cnt++] = 1.0;
            }
            assert( SCIPvarIsTransformed(consdata->vars[i]) );
            vars[cnt] = consdata->vars[i];
            vals[cnt++] = -1.0;

            assert( cnt >= 2 );
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_couple", SCIPconsGetName(cons));
            /* not initial, separate, do not enforce, do not check, propagate, not local, not modifiable, dynamic, removable, not sticking */
            SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, cnt, vars, vals, 0.0, 0.0,
                  FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, newcons) );
            SCIPdebugPrintCons(scip, newcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            ++(*naddedconss);
         }

         /* add south flow conservation constraint */

         /* incoming variables */
         cnt = 0;
         if ( varprevss != NULL )
         {
            vars[cnt] = varprevss;
            vals[cnt++] = 1.0;
         }
         if ( varprevns != NULL )
         {
            vars[cnt] = varprevns;
            vals[cnt++] = 1.0;
         }

         /* outgoing variables */
         if ( varss != NULL )
         {
            vars[cnt] = varss;
            vals[cnt++] = -1.0;
         }
         if ( varsn != NULL )
         {
            vars[cnt] = varsn;
            vals[cnt++] = -1.0;
         }

         assert( cnt >= 2 );
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_south", SCIPconsGetName(cons));
         /* not initial, separate, do not enforce, do not check, propagate, not local, not modifiable, dynamic, removable, not sticking */
         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, cnt, vars, vals, 0.0, 0.0,
               FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugPrintCons(scip, newcons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         ++(*naddedconss);
      }

      /* add north flow conservation constraint */

      /* incoming variables */
      cnt = 0;
      if ( varprevnn != NULL )
      {
         vars[cnt] = varprevnn;
         vals[cnt++] = 1.0;
      }
      if ( varprevsn != NULL )
      {
         vars[cnt] = varprevsn;
         vals[cnt++] = 1.0;
      }

      /* outgoing variables */
      if ( varnn != NULL )
      {
         vars[cnt] = varnn;
         vals[cnt++] = -1.0;
      }
      if ( varns != NULL )
      {
         vars[cnt] = varns;
         vals[cnt++] = -1.0;
      }

      assert( cnt >= 2 );
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_north", SCIPconsGetName(cons));
      if ( i == 0 )
         rhs = -1.0;
      else
         rhs = 0.0;

      /* not initial, separate, do not enforce, do not check, propagate, not local, not modifiable, dynamic, removable, not sticking */
      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, cnt, vars, vals, rhs, rhs,
            FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIPdebugPrintCons(scip, newcons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      ++(*naddedconss);

      /* store variables */
      consdata->extvars[4*i] = varnn;      /*lint !e679*/
      consdata->extvars[4*i + 1] = varns;  /*lint !e679*/
      consdata->extvars[4*i + 2] = varsn;  /*lint !e679*/
      consdata->extvars[4*i + 3] = varss;  /*lint !e679*/

      if ( varnn != NULL )
         ++(consdata->nextvars);
      if ( varns != NULL )
         ++(consdata->nextvars);
      if ( varsn != NULL )
         ++(consdata->nextvars);
      if ( varss != NULL )
         ++(consdata->nextvars);

      /* store previous variables */
      varprevnn = varnn;
      varprevns = varns;
      varprevsn = varsn;
      varprevss = varss;
   }

   return SCIP_OKAY;
}

/** adds extended asymmetric formulation
 *
 *  The extended asymmetric formulation is constructed as follows: Let \f$x_1, \dots, x_k\f$ be the variables contained
 *  in the given XOR constraint. We introduce variables \f$p_1, \ldots, p_k\f$ with the following constraints: \f$p_1 =
 *  x_1\f$, \f$p_k = 1\f$, and for \f$i = 2, \ldots, k-1\f$:
 *  \f[
 *    \begin{array}{ll}
 *      p_i & \leq p_{i-1} + x_i\\
 *      p_i & \leq 2 - (p_{i-1} + x_i)\\
 *      p_i & \geq p_{i-1} - x_i\\
 *      p_i & \geq x_i - p_{i-1}.
 *    \end{array}
 *  \f]
 *  This formulation is described in
 *
 *  Robert D. Carr and Goran Konjevod@n
 *  Polyhedral combinatorics@n
 *  In Harvey Greenberg, editor, Tutorials on emerging methodologies and applications in Operations Research,@n
 *  Chapter 2, pages (2-1)-(2-48). Springer, 2004.
 */
static
SCIP_RETCODE addExtendedAsymmetricFormulation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   int*                  naddedconss         /**< number of added constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* consdata;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];
   SCIP_VAR* prevvar = NULL;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( naddedconss != NULL );
   *naddedconss = 0;

   /* exit if contraints is modifiable */
   if ( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* exit if extended formulation has been added already */
   if ( consdata->extvars != NULL )
      return SCIP_OKAY;

   /* xor constraints with at most 3 variables are handled directly through rows for the convex hull */
   if ( consdata->nvars <= 3 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Add extended formulation for xor constraint <%s> ...\n", SCIPconsGetName(cons));
   assert( consdata->extvars == NULL );
   assert( consdata->nextvars == 0 );

   /* get storage for auxiliary variables */
   consdata->extvarssize = consdata->nvars;
   consdata->nextvars = consdata->nvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->extvars), consdata->extvarssize ) );

   /* pass through components */
   for (i = 0; i < consdata->nvars; ++i)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool redundant = FALSE;
      SCIP_Bool aggregated = FALSE;
      SCIP_CONS* newcons;
      SCIP_VAR* artvar = NULL;
      SCIP_Real lb = 0.0;
      SCIP_Real ub = 1.0;

      /* determine fixing for last variables */
      if ( i == consdata->nvars-1 )
      {
         if ( consdata->rhs )
         {
            lb = 1.0;
            ub = 1.0;
         }
         else
         {
            lb = 0.0;
            ub = 0.0;
         }
      }

      /* create variable */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "p_%s_%d", SCIPconsGetName(cons), i);
      SCIP_CALL( SCIPcreateVar(scip, &artvar, name, lb, ub, 0.0, SCIP_VARTYPE_IMPLINT, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, artvar) );
      SCIP_CALL( SCIPlockVarCons(scip, artvar, cons, TRUE, TRUE) );

      /* create constraints */
      if ( i == 0 )
      {
         /* aggregate artificial variable with original variable */
         SCIP_CALL( SCIPaggregateVars(scip, artvar, consdata->vars[0], 1.0, -1.0, 0.0, &infeasible, &redundant, &aggregated) );
         assert( ! infeasible );
         assert( redundant );
         assert( aggregated );
      }
      else
      {
         assert( SCIPvarIsTransformed(consdata->vars[i]) );

         /* add first constraint */
         vars[0] = artvar;
         vals[0] = 1.0;
         vars[1] = prevvar;
         vals[1] = -1.0;
         vars[2] = consdata->vars[i];
         vals[2] = -1.0;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_1", SCIPconsGetName(cons), i);
         /* not initial, separate, do not enforce, do not check, propagate, not local, not modifiable, dynamic, removable, not sticking */
         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, 3, vars, vals, -SCIPinfinity(scip), 0.0,
               FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugPrintCons(scip, newcons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         ++(*naddedconss);

         /* add second constraint */
         vars[0] = artvar;
         vals[0] = 1.0;
         vars[1] = prevvar;
         vals[1] = 1.0;
         vars[2] = consdata->vars[i];
         vals[2] = 1.0;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_2", SCIPconsGetName(cons), i);
         /* not initial, separate, do not enforce, do not check, propagate, not local, not modifiable, dynamic, removable, not sticking */
         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, 3, vars, vals, -SCIPinfinity(scip), 2.0,
               FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugPrintCons(scip, newcons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         ++(*naddedconss);

         /* add third constraint */
         vars[0] = artvar;
         vals[0] = -1.0;
         vars[1] = prevvar;
         vals[1] = 1.0;
         vars[2] = consdata->vars[i];
         vals[2] = -1.0;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_3", SCIPconsGetName(cons), i);
         /* not initial, separate, do not enforce, do not check, propagate, not local, not modifiable, dynamic, removable, not sticking */
         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, 3, vars, vals, -SCIPinfinity(scip), 0.0,
               FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugPrintCons(scip, newcons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         ++(*naddedconss);

         /* add fourth constraint */
         vars[0] = artvar;
         vals[0] = -1.0;
         vars[1] = prevvar;
         vals[1] = -1.0;
         vars[2] = consdata->vars[i];
         vals[2] = 1.0;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_4", SCIPconsGetName(cons), i);
         /* not initial, separate, do not enforce, do not check, propagate, not local, not modifiable, dynamic, removable, not sticking */
         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, 3, vars, vals, -SCIPinfinity(scip), 0.0,
               FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugPrintCons(scip, newcons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         ++(*naddedconss);
      }

      /* store variable */
      consdata->extvars[i] = artvar;
      prevvar = artvar;
   }

   return SCIP_OKAY;
}

/** creates LP row corresponding to xor constraint:
 *    x1 + ... + xn - 2q == rhs
 *  with internal integer variable q;
 *  in the special case of 3 variables and c = 0, the following linear system is created:
 *    + x - y - z <= 0
 *    - x + y - z <= 0
 *    - x - y + z <= 0
 *    + x + y + z <= 2
 *  in the special case of 3 variables and c = 1, the following linear system is created:
 *    - x + y + z <= 1
 *    + x - y + z <= 1
 *    + x + y - z <= 1
 *    - x - y - z <= -1
 */
static
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   char varname[SCIP_MAXSTRLEN];

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows[0] == NULL);

   if( SCIPconsIsModifiable(cons) || consdata->nvars != 3 )
   {
      SCIP_Real rhsval;

      /* create internal variable, if not yet existing */
      if( consdata->intvar == NULL )
      {
         int ub;

         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "XOR_artificial_%s_int", SCIPconsGetName(cons));
         ub = consdata->nvars/2;
         SCIP_CALL( SCIPcreateVar(scip, &consdata->intvar, varname, 0.0, (SCIP_Real)ub, 0.0,
               consdata->nvars >= 4 ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_BINARY,
               SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, consdata->intvar) );

#ifdef WITH_DEBUG_SOLUTION
         if( SCIPdebugIsMainscip(scip) )
         {
            SCIP_Real solval;
            int count = 0;
            int v;

            for( v = consdata->nvars - 1; v >= 0; --v )
            {
               SCIP_CALL( SCIPdebugGetSolVal(scip, consdata->vars[v], &solval) );
               count += (solval > 0.5 ? 1 : 0);
            }
            assert((count - consdata->rhs) % 2 == 0);
            solval = (SCIP_Real) ((count - consdata->rhs) / 2);

            /* store debug sol value of artificial integer variable */
            SCIP_CALL( SCIPdebugAddSolVal(scip, consdata->intvar, solval) );
         }
#endif

         /* install the rounding locks for the internal variable */
         SCIP_CALL( lockRounding(scip, cons, consdata->intvar) );
      }

      /* create LP row */
      rhsval = (consdata->rhs ? 1.0 : 0.0);
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[0], SCIPconsGetHdlr(cons), SCIPconsGetName(cons), rhsval, rhsval,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[0], consdata->intvar, -2.0) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[0], consdata->nvars, consdata->vars, 1.0) );
   }
   else if( !consdata->rhs )
   {
      char rowname[SCIP_MAXSTRLEN];
      int r;

      /* create the <= 0 rows with one positive sign */
      for( r = 0; r < 3; ++r )
      {
         int v;

         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), r);
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[r], SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), 0.0,
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         for( v = 0; v < 3; ++v )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[r], consdata->vars[v], v == r ? +1.0 : -1.0) );
         }
      }

      /* create the <= 2 row with all positive signs */
      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_3", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[3], SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), 2.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[3], consdata->nvars, consdata->vars, 1.0) );

      /* create extra LP row if integer variable exists */
      if( consdata->intvar != NULL )
      {
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[4], SCIPconsGetHdlr(cons), SCIPconsGetName(cons), 0.0, 0.0,
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[4], consdata->intvar, -2.0) );
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[4], consdata->nvars, consdata->vars, 1.0) );
      }
   }
   else
   {
      char rowname[SCIP_MAXSTRLEN];
      int r;

      /* create the <= 1 rows with one negative sign */
      for( r = 0; r < 3; ++r )
      {
         int v;

         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), r);
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[r], SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), 1.0,
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         for( v = 0; v < 3; ++v )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[r], consdata->vars[v], v == r ? -1.0 : +1.0) );
         }
      }

      /* create the <= -1 row with all negative signs */
      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_3", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[3], SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), -1.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[3], consdata->nvars, consdata->vars, -1.0) );

      /* create extra LP row if integer variable exists */
      if( consdata->intvar != NULL )
      {
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->rows[4], SCIPconsGetHdlr(cons), SCIPconsGetName(cons), 1.0, 1.0,
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[4], consdata->intvar, -2.0) );
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[4], consdata->nvars, consdata->vars, 1.0) );
      }
   }

   return SCIP_OKAY;
}

/** adds linear relaxation of or constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility was detected */
   )
{
   SCIP_CONSDATA* consdata;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(infeasible != NULL);
   assert(!(*infeasible));

   if( consdata->rows[0] == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows[0] != NULL);
   for( r = 0; r < NROWS && !(*infeasible); ++r )
   {
      if( consdata->rows[r] != NULL && !SCIProwIsInLP(consdata->rows[r]) )
      {
         SCIP_CALL( SCIPaddRow(scip, consdata->rows[r], FALSE, infeasible) );
      }
   }

   return SCIP_OKAY;
}

/** returns whether all rows of the LP relaxation are in the current LP */
static
SCIP_Bool allRowsInLP(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->rows[0] == NULL )      /* LP relaxation does not exist */
      return FALSE;
   else
   {
      int r;
      for( r = 0; r < NROWS; ++r )
      {
         if( consdata->rows[r] != NULL && !SCIProwIsInLP(consdata->rows[r]) )
            return FALSE;
      }
      return TRUE;
   }
}

/** checks xor constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   /* check feasibility of constraint if necessary */
   if( checklprows || !allRowsInLP(consdata) )
   {
      SCIP_Real solval;
      SCIP_Bool odd;
      int ones;
      int i;

      /* increase age of constraint; age is reset to zero, if a violation was found only in case we are in
       * enforcement
       */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }

      /* check, if all variables and the rhs sum up to an even value */
      odd = consdata->rhs;
      ones = 0;
      for( i = 0; i < consdata->nvars; ++i )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
         assert(SCIPisFeasIntegral(scip, solval));
         odd = (odd != (solval > 0.5));

         if( solval > 0.5 )
            ++ones;
      }
      if( odd )
      {
         *violated = TRUE;

         /* only reset constraint age if we are in enforcement */
         if( sol == NULL )
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      else if( consdata->intvar != NULL )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->intvar);
         assert(SCIPisFeasIntegral(scip, solval));

         if( !SCIPisFeasEQ(scip, ones - 2 * solval, (SCIP_Real) consdata->rhs) )
            *violated = TRUE;
      }

      /* only reset constraint age if we are in enforcement */
      if( *violated && sol == NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      /* update constraint violation in solution */
      else if ( *violated && sol != NULL )
         SCIPupdateSolConsViolation(scip, sol, 1.0, 1.0);
   }

   return SCIP_OKAY;
}

/** separates current LP solution
 *
 *  Consider a XOR-constraint
 *  \f[
 *    x_1 \oplus x_2 \oplus \dots \oplus x_n = b
 *  \f]
 *  with \f$b \in \{0,1\}\f$ and a solution \f$x^*\f$ to be cut off. Small XOR constraints are handled by adding the
 *  inequalities of the convex hull.
 *
 *  The separation of larger XOR constraints has been described by @n
 *  Xiaojie Zhang and Paul H. Siegel@n
 *  "Adaptive Cut Generation Algorithm for Improved Linear Programming Decoding of Binary Linear Codes"@n
 *  IEEE Transactions on Information Theory, vol. 58, no. 10, 2012
 *
 *  We separate the inequalities
 *  \f[
 *    \sum_{j \in S} (1 - x_j) + \sum_{j \notin S} x_j \geq 1
 *  \f]
 *  with \f$|S| \equiv (b+1) \mbox{ mod } 2\f$ as follows. That these inequalities are valid can be seen as follows: Let
 *  \f$x\f$ be a feasible solution and suppose that the inequality is violated for some \f$S\f$. Then \f$x_j = 1\f$ for
 *  all \f$j \in S\f$ and \f$x_j = 0\f$ for all \f$j \notin S\f$. Thus we should have
 *  \f[
 *    \oplus_{j \in S} x_j = |S| \mbox{ mod } 2 = b+1 \mbox{ mod } 2,
 *  \f]
 *  which is not equal to \f$b\f$ as required by the XOR-constraint.
 *
 *  Let \f$L= \{j \;:\; x^*_j > \frac{1}{2}\}\f$. Suppose that \f$|L|\f$ has @em not the same parity as \f$b\f$ rhs. Then
 *  \f[
 *    \sum_{j \in L} (1 - x_j) + \sum_{j \notin L} x_j \geq 1
 *  \f]
 *  is the only inequality that can be violated. We rewrite the inequality as
 *  \f[
 *     \sum_{j \in L} x_j - \sum_{j \notin L} x_j \leq |L| - 1.
 *  \f]
 *  These inequalities are added.
 *
 *  Otherwise let \f$k = \mbox{argmin}\{x^*_j \;:\; j \in L\}\f$ and check the inequality for \f$L \setminus \{k\}\f$
 *  and similarly for \f$k = \mbox{argmax}\{x^*_j \;:\; j \in L\}\f$.
 */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool             separateparity,     /**< should parity inequalities be separated? */
   SCIP_Bool*            separated,          /**< pointer to store whether a cut was found */
   SCIP_Bool*            cutoff              /**< whether a cutoff has been detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real feasibility;
   int r;

   assert( separated != NULL );
   assert( cutoff != NULL );
   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *separated = FALSE;

   /* create row for the linear relaxation */
   if( consdata->rows[0] == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows[0] != NULL);

   /* test rows for feasibility and add it, if it is infeasible */
   for( r = 0; r < NROWS; ++r )
   {
      if( consdata->rows[r] != NULL && !SCIProwIsInLP(consdata->rows[r]) )
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

   /* separate parity inequalities if required */
   if ( separateparity && consdata->nvars > 3 )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_Real maxval = -1.0;
      SCIP_Real minval = 2.0;
      SCIP_Real sum = 0.0;
      int maxidx = -1;
      int minidx = -1;
      int ngen = 0;
      int cnt = 0;
      int j;

      SCIPdebugMsg(scip, "separating parity inequalities ...\n");

      /* compute value */
      for (j = 0; j < consdata->nvars; ++j)
      {
         SCIP_Real val;

         val = SCIPgetSolVal(scip, sol, consdata->vars[j]);
         if ( SCIPisFeasGT(scip, val, 0.5) )
         {
            if ( val < minval )
            {
               minval = val;
               minidx = j;
            }
            ++cnt;
            sum += (1.0 - val);
         }
         else
         {
            if ( val > maxval )
            {
               maxval = val;
               maxidx = j;
            }
            sum += val;
         }
      }

      /* if size of set does not have the same parity as rhs (e.g., size is odd if rhs is 0) */
      if ( (cnt - consdata->rhs) % 2 == 1 )
      {
         if ( SCIPisEfficacious(scip, 1.0 - sum) )
         {
            SCIP_ROW* row;

            SCIPdebugMsg(scip, "found violated parity cut (efficiacy: %f)\n", 1.0 - sum);

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "parity#%s", SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), (SCIP_Real) (cnt - 1), FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            /* fill in row */
            for (j = 0; j < consdata->nvars; ++j)
            {
               if ( SCIPisFeasGT(scip, SCIPgetSolVal(scip, sol, consdata->vars[j]), 0.5) )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], 1.0) );
               }
               else
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], -1.0) );
               }
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
            assert( SCIPisGT(scip, SCIPgetRowLPActivity(scip, row), (SCIP_Real) (cnt-1)) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
            ++ngen;
         }
      }
      else
      {
         /* If the parity is equal: check removing the element with smallest value from the set and adding the
          * element with largest value to the set. If we remove the element with smallest value, we have to subtract (1
          * - minval) and add minval to correct the sum. */
         if ( SCIPisEfficacious(scip, 1.0 - (sum - 1.0 + 2.0 * minval)) )
         {
            SCIP_ROW* row;

            SCIPdebugMsg(scip, "found violated parity cut (efficiacy: %f, minval: %f)\n", 1.0 - (sum - 1.0 + 2.0 * minval), minval);

            /* the rhs of the inequality is the corrected set size minus 1 */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "parity#%s", SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), (SCIP_Real) (cnt - 2), FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            /* fill in row */
            for (j = 0; j < consdata->nvars; ++j)
            {
               if ( SCIPisFeasGT(scip, SCIPgetSolVal(scip, sol, consdata->vars[j]), 0.5) )
               {
                  /* if the index corresponds to the smallest element, we reverse the sign */
                  if ( j == minidx )
                     SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], -1.0) );
                  else
                     SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], 1.0) );
               }
               else
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], -1.0) );
               }
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
            assert( SCIPisGT(scip, SCIPgetRowLPActivity(scip, row), (SCIP_Real) (cnt-2)) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
            ++ngen;
         }

         /* If we add the element with largest value, we have to add (1 - maxval) and subtract maxval to get the correct sum. */
         if ( SCIPisEfficacious(scip, 1.0 - (sum + 1.0 - 2.0 * maxval)) )
         {
            SCIP_ROW* row;

            SCIPdebugMsg(scip, "found violated parity cut (efficiacy: %f, maxval: %f)\n", 1.0 - (sum + 1.0 - 2.0 * maxval), maxval);

            /* the rhs of the inequality is the size of the corrected set */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "parity#%s", SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), (SCIP_Real) cnt, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            /* fill in row */
            for (j = 0; j < consdata->nvars; ++j)
            {
               if ( SCIPisFeasGT(scip, SCIPgetSolVal(scip, sol, consdata->vars[j]), 0.5) )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], 1.0) );
               }
               else
               {
                  /* if the index corresponds to the largest element, we reverse the sign */
                  if ( j == maxidx )
                     SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], 1.0) );
                  else
                     SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], -1.0) );
               }
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
            assert( *cutoff || SCIPisGT(scip, SCIPgetRowLPActivity(scip, row), (SCIP_Real)(j-1)) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
            ++ngen;
         }
      }

      SCIPdebugMsg(scip, "separated parity inequalites: %d\n", ngen);
      if ( ngen > 0 )
         *separated = TRUE;
   }

   return SCIP_OKAY;
}

/** Transform linear system \f$A x = b\f$ into row echolon form via the Gauss algorithm with row pivoting over GF2
 *  @returns the rank of @p A
 *
 *  Here, \f$A \in R^{m \times n},\; b \in R^m\f$. On exit, the vector @p p contains a permutation of the row indices
 *  used for pivoting and the function returns the rank @p r of @p A. For each row \f$i = 1, \ldots, r\f$, the entry @p
 *  s[i] contains the column index of the first nonzero in row @p i.
 */
static
int computeRowEcholonGF2(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   m,                  /**< number of rows */
   int                   n,                  /**< number of columns */
   int*                  p,                  /**< row permutation */
   int*                  s,                  /**< steps indicators of the row echolon form */
   Type**                A,                  /**< matrix */
   Type*                 b                   /**< rhs */
   )
{
   int pi;
   int i;
   int j;
   int k;

   assert( A != NULL );
   assert( b != NULL );
   assert( p != NULL );
   assert( s != NULL );

   /* init permutation and step indicators */
   for (i = 0; i < m; ++i)
   {
      p[i] = i;
      s[i] = i;
   }

   /* loop through possible steps in echolon form (stop at min {n, m}) */
   for (i = 0; i < m && i < n; ++i)
   {
      assert( s[i] == i );

      /* init starting column */
      if ( i == 0 )
         j = 0;
      else
         j = s[i-1] + 1;

      /* find pivot row (i.e., first nonzero entry), if all entries in current row are 0 we search the next column */
      do
      {
         /* search in current column j */
         k = i;
         while ( k < m && A[p[k]][j] == 0 )
            ++k;

         /* found pivot */
         if ( k < m )
            break;

         /* otherwise search next column */
         ++j;
      }
      while ( j < n );

      /* if not pivot entry was found (checked all columns), the rank of A is equal to the current index i; in this case
       * all entries in and below row i are 0 */
      if ( j >= n )
         return i;

      /* at this place: we have found a pivot entry (p[k], j) */
      assert( k < m );

      /* store step index */
      s[i] = j;
      assert( A[p[k]][j] != 0 );

      /* swap row indices */
      if ( k != i )
      {
         int h = p[i];
         p[i] = p[k];
         p[k] = h;
      }
      pi = p[i];
      assert( A[pi][s[i]] != 0 );

      /* do elimination */
      for (k = i+1; k < m; ++k)
      {
         int pk = p[k];
         /* if entry in leading column is nonzero (otherwise we already have a 0) */
         if ( A[pk][s[i]] != 0 )
         {
            for (j = s[i]; j < n; ++j)
               A[pk][j] = A[pk][j] ^ A[pi][j];
            b[pk] = b[pk] ^ b[pi];
         }
      }

      /* check stopped (only every 100 rows in order to save time */
      if ( i % 100 == 99 )
      {
         if ( SCIPisStopped(scip) )
            return -1;
      }
   }

   /* at this point we have treated all rows in which a step can occur; the rank is the minimum of the number of rows or
    * columns min {n,m}. */
   if ( n <= m )
      return n;
   return m;
}

/** Construct solution from matrix in row echolon form over GF2
 *
 *  Compute solution of \f$A x = b\f$, which is already in row echolon form (@see computeRowEcholonGF2()) */
static
void solveRowEcholonGF2(
   int                   m,                  /**< number of rows */
   int                   n,                  /**< number of columns */
   int                   r,                  /**< rank of matrix */
   int*                  p,                  /**< row permutation */
   int*                  s,                  /**< steps indicators of the row echolon form */
   Type**                A,                  /**< matrix */
   Type*                 b,                  /**< rhs */
   Type*                 x                   /**< solution vector on exit */
   )
{
   int i;
   int k;

   assert( A != NULL );
   assert( b != NULL );
   assert( s != NULL );
   assert( p != NULL );
   assert( x != NULL );
   assert( r <= m && r <= n );

   /* init solution vector to 0 */
   for (k = 0; k < n; ++k)
      x[k] = 0;

   /* init last entry */
   x[s[r-1]] = b[p[r-1]];

   /* loop backwards through solution vector */
   for (i = r-2; i >= 0; --i)
   {
      Type val;

      assert( i <= s[i] && s[i] <= n );

      /* init val with rhs and then add the contributions of the components of x already computed */
      val = b[p[i]];
      for (k = i+1; k < r; ++k)
      {
         assert( i <= s[k] && s[k] <= n );
         if ( A[p[i]][s[k]] != 0 )
            val = val ^ x[s[k]];
      }

      /* store solution */
      x[s[i]] = val;
   }
}

/** solve equation system over GF 2 by Gauss algorithm and create solution out of it or return cutoff
 *
 *  Collect all information in xor constraints into a linear system over GF2. Then solve the system by computing a row
 *  echolon form. If the system is infeasible, the current node is infeasible. Otherwise, we can compute a solution for
 *  the xor constraints given. We check whether this gives a solution for the whole problem.
 *
 *  We sort the columns with respect to the product of the objective coefficients and 1 minus the current LP solution
 *  value. The idea is that columns that are likely to provide the steps in the row echolon form should appear towards
 *  the front of the matrix. The smaller the product, the more it makes sense to set the variable to 1 (because the
 *  solution value is already close to 1 and the objective function is small).
 *
 *  Note that this function is called from propagation where usually no solution is available. However, the solution is
 *  only used for sorting the columns. Thus, the procedure stays correct even with nonsense solutions.
 */
static
SCIP_RETCODE checkSystemGF2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< xor constraints */
   int                   nconss,             /**< number of xor constraints */
   SCIP_SOL*             currentsol,         /**< current solution (maybe NULL) */
   SCIP_RESULT*          result              /**< result of propagation (possibly cutoff, no change if primal solution has been tried) */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP* varhash;
   SCIP_Bool* xoractive;
   SCIP_Real* xorvals;
   SCIP_VAR** xorvars;
   SCIP_Bool noaggr = TRUE;
   Type** A;
   Type* b;
   int* s;
   int* p;
   int* xoridx;
   int* xorbackidx;
   int nconssactive = 0;
   int nconssmat = 0;
   int nvarsmat = 0;
   int nvars;
   int rank;
   int i;
   int j;

   assert( scip != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   if ( *result == SCIP_CUTOFF )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Checking feasibility via the linear equation system over GF2 using Gauss.\n");

   nvars = SCIPgetNVars(scip);

   /* set up hash map from variable to column index */
   SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xoractive, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xorvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xoridx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xorvals, nvars) );

   /* collect variables */
   for (i = 0; i < nconss; ++i)
   {
      int cnt = 0;

      xoractive[i] = FALSE;

      assert( conss[i] != NULL );
      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );

      /* count nonfixed variables in constraint */
      for (j = 0; j < consdata->nvars; ++j)
      {
         SCIP_VAR* var;

         var = consdata->vars[j];
         assert( var != NULL );
         assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );

         /* replace negated variables */
         if ( SCIPvarIsNegated(var) )
            var = SCIPvarGetNegatedVar(var);
         assert( var != NULL );

         /* consider nonfixed variables */
         if ( SCIPcomputeVarLbLocal(scip, var) < 0.5 && SCIPcomputeVarUbLocal(scip, var) > 0.5 )
         {
            /* consider active variables and collect only new ones */
            if ( SCIPvarIsActive(var) && ! SCIPhashmapExists(varhash, var) )
            {
               /* add variable in map */
               SCIP_CALL( SCIPhashmapInsert(varhash, var, (void*) (size_t) nvarsmat) );
               assert( nvarsmat == (int) (size_t) SCIPhashmapGetImage(varhash, var) );
               xorvals[nvarsmat] = SCIPvarGetObj(var) * (1.0 - SCIPgetSolVal(scip, currentsol, var));
               xorvars[nvarsmat++] = var;
            }
            ++cnt;
         }
      }

      if ( cnt > 0 )
      {
         xoractive[i] = TRUE;
         ++nconssactive;
      }
#ifdef SCIP_DISABLED_CODE
      /* The following can save time, if there are constraints with all variables fixed that are infeasible; this
       * should, however, be detected somewhere else, e.g., in propagateCons(). */
      else
      {
         /* all variables are fixed - check whether constraint is feasible (could be that the constraint is not propagated) */
         assert( cnt == 0 );
         for (j = 0; j < consdata->nvars; ++j)
         {
            /* count variables fixed to 1 */
            if ( SCIPcomputeVarLbLocal(scip, consdata->vars[j]) > 0.5 )
               ++cnt;
            else
               assert( SCIPcomputeVarUbLocal(scip, consdata->vars[j]) < 0.5 );
         }
         if ( ( cnt - consdata->rhs ) % 2 != 0 )
         {
            SCIPdebugMsg(scip, "constraint <%s> with all variables fixed is violated.\n", SCIPconsGetName(conss[i]));
            *result = SCIP_CUTOFF;
            break;
         }
      }
#endif
   }
   assert( nvarsmat <= nvars );
   assert( nconssactive <= nconss );

   if ( nconssactive > MAXXORCONSSSYSTEM || nvarsmat > MAXXORVARSSYSTEM || *result == SCIP_CUTOFF )
   {
      SCIPdebugMsg(scip, "Skip checking the xor system over GF2 (%d conss, %d vars).\n", nconssactive, nvarsmat);
      SCIPfreeBufferArray(scip, &xorvals);
      SCIPfreeBufferArray(scip, &xoridx);
      SCIPfreeBufferArray(scip, &xorvars);
      SCIPfreeBufferArray(scip, &xoractive);
      SCIPhashmapFree(&varhash);
      return SCIP_OKAY;
   }

   /* init index */
   for (j = 0; j < nvarsmat; ++j)
      xoridx[j] = j;

   /* Sort variables non-decreasingly with respect to product of objective and 1 minus the current solution value: the
    * smaller the value the better it would be to set the variable to 1. This is more likely if the variable appears
    * towards the front of the matrix, because only the entries on the steps in the row echolon form will have the
    * chance to be nonzero.
    */
   SCIPsortRealIntPtr(xorvals, xoridx, (void**) xorvars, nvarsmat);
   SCIPfreeBufferArray(scip, &xorvals);

   /* build back index */
   SCIP_CALL( SCIPallocBufferArray(scip, &xorbackidx, nvarsmat) );
   for (j = 0; j < nvarsmat; ++j)
   {
      assert( 0 <= xoridx[j] && xoridx[j] < nvarsmat );
      xorbackidx[xoridx[j]] = j;
   }

   /* init matrix and rhs */
   SCIP_CALL( SCIPallocBufferArray(scip, &b, nconssactive) );
   SCIP_CALL( SCIPallocBufferArray(scip, &A, nconssactive) );
   for (i = 0; i < nconss; ++i)
   {
      if ( ! xoractive[i] )
         continue;

      assert( conss[i] != NULL );
      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );
      assert( consdata->nvars > 0 );

      SCIP_CALL( SCIPallocBufferArray(scip, &(A[nconssmat]), nvarsmat) ); /*lint !e866*/
      BMSclearMemoryArray(A[nconssmat], nvarsmat); /*lint !e866*/

      /* correct rhs w.r.t. to fixed variables and count nonfixed variables in constraint */
      b[nconssmat] = (Type) consdata->rhs;
      for (j = 0; j < consdata->nvars; ++j)
      {
         SCIP_VAR* var;
         int idx;

         var = consdata->vars[j];
         assert( var != NULL );

         /* replace negated variables */
         if ( SCIPvarIsNegated(var) )
         {
            var = SCIPvarGetNegatedVar(var);
            assert( var != NULL );
            b[nconssmat] = ! b[nconssmat];
         }


         /* replace aggregated variables */
         while( SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED )
         {
            SCIP_Real scalar;
            SCIP_Real constant;

            scalar = SCIPvarGetAggrScalar(var);
            constant = SCIPvarGetAggrConstant(var);

            /* the variable resolves to a constant, we just update the rhs */
            if( SCIPisEQ(scip, scalar, 0.0) )
            {
               assert(SCIPisEQ(scip, constant, 0.0) || SCIPisEQ(scip, constant, 1.0));
               if( SCIPisEQ(scip, constant, 1.0) )
                  b[nconssmat] = ! b[nconssmat];
               var = NULL;
               break;
            }
            /* replace aggregated variable by active variable and update rhs, if needed */
            else
            {
               assert(SCIPisEQ(scip, scalar, 1.0) || SCIPisEQ(scip, scalar, -1.0));
               if( SCIPisEQ(scip, constant, 1.0) )
                  b[nconssmat] = ! b[nconssmat];

               var = SCIPvarGetAggrVar(var);
               assert(var != NULL);
            }
         }
         /* variable resolved to a constant */
         if( var == NULL )
            continue;

         /* If the constraint contains multiaggregated variables, the solution might not be valid, since the
          * implications are not represented in the matrix.
          */
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
            noaggr = FALSE;

         if ( SCIPcomputeVarLbLocal(scip, var) > 0.5 )
         {
            /* variable is fixed to 1, invert rhs */
            b[nconssmat] = ! b[nconssmat];
            assert( ! SCIPhashmapExists(varhash, var) );
         }
         else
         {
            assert(SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
               || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
            if ( SCIPvarIsActive(var) && SCIPcomputeVarUbLocal(scip, var) > 0.5 )
            {
               assert( SCIPhashmapExists(varhash, var) );
               idx = (int) (size_t) SCIPhashmapGetImage(varhash, var);
               assert( idx < nvarsmat );
               assert( 0 <= xorbackidx[idx] && xorbackidx[idx] < nvarsmat );
               A[nconssmat][xorbackidx[idx]] = 1;
            }
         }
      }
      ++nconssmat;
   }
   SCIPdebugMsg(scip, "Found %d non-fixed variables in %d nonempty xor constraints.\n", nvarsmat, nconssmat);
   assert( nconssmat == nconssactive );

   /* perform Gauss algorithm */
   SCIP_CALL( SCIPallocBufferArray(scip, &p, nconssmat) );
   SCIP_CALL( SCIPallocBufferArray(scip, &s, nconssmat) );

#ifdef SCIP_OUTPUT
   SCIPinfoMessage(scip, NULL, "Matrix before Gauss (size: %d x %d):\n", nconssmat, nvarsmat);
   for (i = 0; i < nconssmat; ++i)
   {
      for (j = 0; j < nvarsmat; ++j)
         SCIPinfoMessage(scip, NULL, "%d ", A[i][j]);
      SCIPinfoMessage(scip, NULL, " = %d\n", b[i]);
   }
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   rank = -1;
   if ( ! SCIPisStopped(scip) )
   {
      rank = computeRowEcholonGF2(scip, nconssmat, nvarsmat, p, s, A, b);
      assert( rank <= nconssmat && rank <= nvarsmat );
   }

   /* rank is < 0 if the solution process has been stopped */
   if ( rank >= 0 )
   {
#ifdef SCIP_OUTPUT
      SCIPinfoMessage(scip, NULL, "Matrix after Gauss (rank: %d):\n", rank);
      for (i = 0; i < nconssmat; ++i)
      {
         for (j = 0; j < nvarsmat; ++j)
            SCIPinfoMessage(scip, NULL, "%d ", A[p[i]][j]);
         SCIPinfoMessage(scip, NULL, " = %d\n", b[p[i]]);
      }
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* check whether system is feasible */
      for (i = rank; i < nconssmat; ++i)
      {
         if ( b[p[i]] != 0 )
            break;
      }

      /* did not find nonzero entry in b -> equation system is feasible */
      if ( i >= nconssmat )
      {
         SCIPdebugMsg(scip, "System feasible with rank %d (nconss=%d)\n", rank, nconssmat);

         /* matrix has full rank, solution is unique */
         if( rank == nvarsmat && noaggr )
         {
            SCIP_Bool tightened;
            SCIP_Bool infeasible;
            Type* x;

            SCIPdebugMsg(scip, "Found unique solution.\n");

            /* construct solution */
            SCIP_CALL( SCIPallocBufferArray(scip, &x, nvarsmat) );
            solveRowEcholonGF2(nconssmat, nvarsmat, rank, p, s, A, b, x);

#ifdef SCIP_OUTPUT
            SCIPinfoMessage(scip, NULL, "Solution:\n");
            for (j = 0; j < nvarsmat; ++j)
               SCIPinfoMessage(scip, NULL, "%d ", x[j]);
            SCIPinfoMessage(scip, NULL, "\n");
#endif

            /* fix variables according to computed unique solution */
            for( j = 0; j < nvarsmat; ++j )
            {
               assert( (int) (size_t) SCIPhashmapGetImage(varhash, xorvars[j]) < nvars );
               assert( xorbackidx[(int) (size_t) SCIPhashmapGetImage(varhash, xorvars[j])] == j );
               assert( SCIPcomputeVarLbLocal(scip, xorvars[j]) < 0.5 );
               if( x[j] == 0 )
               {
                  SCIP_CALL( SCIPtightenVarUb(scip, xorvars[j], 0.0, FALSE, &infeasible, &tightened) );
                  assert(tightened);
                  assert(!infeasible);
               }
               else
               {
                  assert(x[j] == 1);
                  SCIP_CALL( SCIPtightenVarLb(scip, xorvars[j], 1.0, FALSE, &infeasible, &tightened) );
                  assert(tightened);
                  assert(!infeasible);
               }
            }
            SCIPfreeBufferArray(scip, &x);
         }
         /* matrix does not have full rank, we add the solution, but cannot derive fixings */
         else
         {
            SCIP_HEUR* heurtrysol;

            SCIPdebugMsg(scip, "Found solution.\n");

            /* try solution */
            heurtrysol = SCIPfindHeur(scip, "trysol");

            if ( heurtrysol != NULL )
            {
               SCIP_Bool success;
               SCIP_VAR** vars;
               SCIP_SOL* sol;
               Type* x;

               /* construct solution */
               SCIP_CALL( SCIPallocBufferArray(scip, &x, nvarsmat) );
               solveRowEcholonGF2(nconssmat, nvarsmat, rank, p, s, A, b, x);

#ifdef SCIP_OUTPUT
               SCIPinfoMessage(scip, NULL, "Solution:\n");
               for (j = 0; j < nvarsmat; ++j)
                  SCIPinfoMessage(scip, NULL, "%d ", x[j]);
               SCIPinfoMessage(scip, NULL, "\n");
#endif

               /* create solution */
               SCIP_CALL( SCIPcreateSol(scip, &sol, heurtrysol) );

               /* transfer solution */
               for (j = 0; j < nvarsmat; ++j)
               {
                  if ( x[j] != 0 )
                  {
                     assert( (int) (size_t) SCIPhashmapGetImage(varhash, xorvars[j]) < nvars );
                     assert( xorbackidx[(int) (size_t) SCIPhashmapGetImage(varhash, xorvars[j])] == j );
                     assert( SCIPcomputeVarLbLocal(scip, xorvars[j]) < 0.5 );
                     SCIP_CALL( SCIPsetSolVal(scip, sol, xorvars[j], 1.0) );
                  }
               }
               SCIPfreeBufferArray(scip, &x);

               /* add *all* variables fixed to 1 */
               vars = SCIPgetVars(scip);
               for (j = 0; j < nvars; ++j)
               {
                  if ( SCIPcomputeVarLbLocal(scip, vars[j]) > 0.5 )
                  {
                     SCIP_CALL( SCIPsetSolVal(scip, sol, vars[j], 1.0) );
                     SCIPdebugMsg(scip, "Added fixed variable <%s>.\n", SCIPvarGetName(vars[j]));
                  }
               }

               /* correct integral variables if necessary */
               for (i = 0; i < nconss; ++i)
               {
                  consdata = SCIPconsGetData(conss[i]);
                  assert(consdata != NULL);

                  if ( xoractive[i] && consdata->intvar != NULL )
                  {
                     SCIP_Real val;
                     int nones = 0;

                     for (j = 0; j < consdata->nvars; ++j)
                     {
                        if ( SCIPgetSolVal(scip, sol, consdata->vars[j]) > 0.5 )
                           ++nones;
                     }
                     /* if there are aggregated variables, the solution might not be feasible */
                     assert( ! noaggr || nones % 2 == (int) consdata->rhs );
                     if ( (unsigned int) nones != consdata->rhs )
                     {
                        val = (SCIP_Real) (nones - consdata->rhs)/2;
                        if ( SCIPisGE(scip, val, SCIPvarGetLbGlobal(consdata->intvar)) && SCIPisLE(scip, val, SCIPvarGetUbGlobal(consdata->intvar)) )
                        {
                           SCIP_CALL( SCIPsetSolVal(scip, sol, consdata->intvar, val) );
                        }
                     }
                  }
               }
               SCIPdebug( SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) ) );

               /* check feasibility of new solution and pass it to trysol heuristic */
               SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
               if ( success )
               {
                  SCIP_CALL( SCIPheurPassSolAddSol(scip, heurtrysol, sol) );
                  SCIPdebugMsg(scip, "Creating solution was successful.\n");
               }
#ifdef SCIP_DEBUG
               else
               {
                  /* the solution might not be feasible, because of additional constraints */
                  SCIPdebugMsg(scip, "Creating solution was not successful.\n");
               }
#endif
               SCIP_CALL( SCIPfreeSol(scip, &sol) );
            }
         }
      }
      else
      {
         *result = SCIP_CUTOFF;
         SCIPdebugMsg(scip, "System not feasible.\n");
      }
   }

   /* free storage */
   SCIPfreeBufferArray(scip, &s);
   SCIPfreeBufferArray(scip, &p);
   j = nconssmat - 1;
   for (i = nconss - 1; i >= 0 ; --i)
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if ( consdata->nvars == 0 )
         continue;

      if( !xoractive[i] )
         continue;

      SCIPfreeBufferArray(scip, &(A[j]));
      --j;
   }
   SCIPfreeBufferArray(scip, &A);
   SCIPfreeBufferArray(scip, &b);
   SCIPfreeBufferArray(scip, &xorbackidx);
   SCIPfreeBufferArray(scip, &xoridx);
   SCIPfreeBufferArray(scip, &xorvars);
   SCIPfreeBufferArray(scip, &xoractive);
   SCIPhashmapFree(&varhash);

   return SCIP_OKAY;
}

/** for each variable in the xor constraint, add it to conflict set; for integral variable add corresponding bound */
static
SCIP_RETCODE addConflictBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL (not equal to integral variable) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   PROPRULE              proprule            /**< propagation rule */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   nvars = consdata->nvars;

   switch( proprule )
   {
   case PROPRULE_0:
      assert( infervar == NULL || infervar == consdata->intvar );

      /* the integral variable was fixed, because all variables were fixed */
      for (i = 0; i < nvars; ++i)
      {
         assert( SCIPisEQ(scip, SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, FALSE), SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, FALSE)) );
         SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
      }
      break;

   case PROPRULE_1:
      /* the variable was inferred, because all other variables were fixed */
      for (i = 0; i < nvars; ++i)
      {
         /* add variables that were fixed to 1 before */
         if ( SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, FALSE) > 0.5 )
         {
            assert( SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, TRUE) > 0.5 );
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
         /* add variables that were fixed to 0 */
         else if ( SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, FALSE) < 0.5 )
         {
            assert( SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, TRUE) < 0.5 );
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
         else
         {
            /* check changed variable (changed variable is 0 or 1 afterwards) */
            assert( vars[i] == infervar );
         }
      }
      break;

   case PROPRULE_INTLB:
      assert( consdata->intvar != NULL );

      if( infervar != consdata->intvar )
      {
         /* the variable was fixed, because of the lower bound of the integral variable */
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->intvar, NULL) );
      }
      /* to many and the other fixed variables */
      for (i = 0; i < nvars; ++i)
      {
         /* add variables that were fixed to 0 */
         if ( SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, FALSE) < 0.5 )
         {
            assert( SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, TRUE) < 0.5 );
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
      }
      break;

   case PROPRULE_INTUB:
      assert( consdata->intvar != NULL );

      if( infervar != consdata->intvar )
      {
         /* the variable was fixed, because of upper bound of the integral variable and the other fixed variables */
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->intvar, NULL) );
      }
      for (i = 0; i < nvars; ++i)
      {
         /* add variables that were fixed to 1 */
         if ( SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, FALSE) > 0.5 )
         {
            assert( SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, TRUE) > 0.5 );
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
      }
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in xor constraint <%s>\n", proprule, SCIPconsGetName(cons));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint that detected the conflict */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL (not equal to integral variable) */
   PROPRULE              proprule            /**< propagation rule */
   )
{
   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   /* add bound changes */
   SCIP_CALL( addConflictBounds(scip, cons, infervar, NULL, proprule) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagates constraint with the following rules:
 *   (0) all variables are fixed => can fix integral variable
 *   (1) all except one variable fixed  =>  fix remaining variable and integral variable
 *   (2) depending on the amount of fixed binary variables we can tighten the integral variable
 *   (3) depending on the lower bound of the integral variable one can fix variables to 1
 *   (4) depending on the upper bound of the integral variable one can fix variables to 0
 */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars,         /**< pointer to add up the number of fixed variables */
   int*                  nchgbds             /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Bool odd;
   int nvars;
   int nfixedones;
   int nfixedzeros;
   int watchedvar1;
   int watchedvar2;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(nchgbds != NULL);

   /* propagation can only be applied, if we know all operator variables */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   nvars = consdata->nvars;

   /* don't process the constraint, if the watched variables weren't fixed to any value since last propagation call */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   /* propagation cannot be applied, if we have at least two unfixed variables left;
    * that means, we only have to watch (i.e. capture events) of two variables, and switch to other variables
    * if these ones get fixed
    */
   watchedvar1 = consdata->watchedvar1;
   watchedvar2 = consdata->watchedvar2;

   /* check, if watched variables are still unfixed */
   if( watchedvar1 != -1 )
   {
      if( SCIPvarGetLbLocal(vars[watchedvar1]) > 0.5 || SCIPvarGetUbLocal(vars[watchedvar1]) < 0.5 )
         watchedvar1 = -1;
   }
   if( watchedvar2 != -1 )
   {
      if( SCIPvarGetLbLocal(vars[watchedvar2]) > 0.5 || SCIPvarGetUbLocal(vars[watchedvar2]) < 0.5 )
         watchedvar2 = -1;
   }

   /* if only one watched variable is still unfixed, make it the first one */
   if( watchedvar1 == -1 )
   {
      watchedvar1 = watchedvar2;
      watchedvar2 = -1;
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if the watched variables are invalid (fixed), find new ones if existing; count the parity */
   odd = consdata->rhs;
   nfixedones = 0;
   nfixedzeros = 0;
   if( watchedvar2 == -1 )
   {
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarGetLbLocal(vars[i]) > 0.5 )
         {
            odd = !odd;
            ++nfixedones;
         }
         else if( SCIPvarGetUbLocal(vars[i]) > 0.5 )
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
         else if ( SCIPvarGetUbLocal(vars[i]) < 0.5 )
            ++nfixedzeros;
      }
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if all variables are fixed, we can decide the feasibility of the constraint */
   if( watchedvar1 == -1 )
   {
      assert(watchedvar2 == -1);

      if( odd )
      {
         SCIPdebugMsg(scip, "constraint <%s>: all vars fixed, constraint is infeasible\n", SCIPconsGetName(cons));

         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflict(scip, cons, NULL, PROPRULE_0) );
         SCIP_CALL( SCIPresetConsAge(scip, cons) );

         *cutoff = TRUE;
      }
      else
      {
         /* fix integral variable if present */
         if ( consdata->intvar != NULL && !consdata->deleteintvar )
         {
            int fixval;

            assert( ! *cutoff );
            assert( (nfixedones - consdata->rhs) % 2 == 0 );

            fixval = (nfixedones - consdata->rhs)/2; /*lint !e713*/

            SCIPdebugMsg(scip, "fix integral variable <%s> to %d\n", SCIPvarGetName(consdata->intvar), fixval);

            /* check whether value to fix is outside bounds */
            if ( fixval + 0.5 < SCIPvarGetLbLocal(consdata->intvar) )
            {
               /* cannot fix auxiliary variable (maybe it has been branched on): we are infeasible */
               SCIPdebugMsg(scip, "node infeasible: activity is %d, bounds of integral variable are [%g,%g]\n",
                  fixval, SCIPvarGetLbLocal(consdata->intvar), SCIPvarGetUbLocal(consdata->intvar));

               SCIP_CALL( analyzeConflict(scip, cons, NULL, PROPRULE_INTLB) );
               SCIP_CALL( SCIPresetConsAge(scip, cons) );

               *cutoff = TRUE;
            }
            else if ( fixval - 0.5 > SCIPvarGetUbLocal(consdata->intvar) )
            {
               /* cannot fix auxiliary variable (maybe it has been branched on): we are infeasible */
               SCIPdebugMsg(scip, "node infeasible: activity is %d, bounds of integral variable are [%g,%g]\n",
                  fixval, SCIPvarGetLbLocal(consdata->intvar), SCIPvarGetUbLocal(consdata->intvar));

               SCIP_CALL( analyzeConflict(scip, cons, NULL, PROPRULE_INTUB) );
               SCIP_CALL( SCIPresetConsAge(scip, cons) );

               *cutoff = TRUE;
            }
            else
            {
               if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->intvar), (SCIP_Real) fixval) )
               {
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->intvar, (SCIP_Real) fixval, cons, (int)PROPRULE_0, FALSE, &infeasible, &tightened) );
                  assert( tightened );
                  assert( ! infeasible );
               }

               if ( ! SCIPisEQ(scip, SCIPvarGetUbLocal(consdata->intvar), (SCIP_Real) fixval) )
               {
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->intvar, (SCIP_Real) fixval, cons, (int)PROPRULE_0, FALSE, &infeasible, &tightened) );
                  assert( tightened );
                  assert( ! infeasible );
               }

               ++(*nfixedvars);
            }
         }
         else
         {
            SCIPdebugMsg(scip, "constraint <%s>: all vars fixed, constraint is feasible\n", SCIPconsGetName(cons));
         }
      }
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

      return SCIP_OKAY;
   }

   /* if only one variable is not fixed, this variable can be deduced */
   if( watchedvar2 == -1 )
   {
      assert(watchedvar1 != -1);

      SCIPdebugMsg(scip, "constraint <%s>: only one unfixed variable -> fix <%s> to %u\n",
         SCIPconsGetName(cons), SCIPvarGetName(vars[watchedvar1]), odd);

      SCIP_CALL( SCIPinferBinvarCons(scip, vars[watchedvar1], odd, cons, (int)PROPRULE_1, &infeasible, &tightened) );
      assert(!infeasible);
      assert(tightened);

      (*nfixedvars)++;

      /* fix integral variable if present */
      if ( consdata->intvar != NULL && !consdata->deleteintvar )
      {
         int fixval;

         /* if variable has been fixed to 1, adjust number of fixed variables */
         if ( odd )
            ++nfixedones;

         assert( (nfixedones - consdata->rhs) % 2 == 0 );

         fixval = (nfixedones - consdata->rhs)/2; /*lint !e713*/
         SCIPdebugMsg(scip, "should fix integral variable <%s> to %d\n", SCIPvarGetName(consdata->intvar), fixval);

         /* check whether value to fix is outside bounds */
         if ( fixval + 0.5 < SCIPvarGetLbLocal(consdata->intvar) )
         {
            /* cannot fix auxiliary variable (maybe it has been branched on): we are infeasible */
            SCIPdebugMsg(scip, "node infeasible: activity is %d, bounds of integral variable are [%g,%g]\n",
               fixval, SCIPvarGetLbLocal(consdata->intvar), SCIPvarGetUbLocal(consdata->intvar));

            SCIP_CALL( analyzeConflict(scip, cons, NULL, PROPRULE_INTLB) );
            SCIP_CALL( SCIPresetConsAge(scip, cons) );

            *cutoff = TRUE;
         }
         else if ( fixval - 0.5 > SCIPvarGetUbLocal(consdata->intvar) )
         {
            /* cannot fix auxiliary variable (maybe it has been branched on): we are infeasible */
            SCIPdebugMsg(scip, "node infeasible: activity is %d, bounds of integral variable are [%g,%g]\n",
               fixval, SCIPvarGetLbLocal(consdata->intvar), SCIPvarGetUbLocal(consdata->intvar));

            SCIP_CALL( analyzeConflict(scip, cons, NULL, PROPRULE_INTUB) );
            SCIP_CALL( SCIPresetConsAge(scip, cons) );

            *cutoff = TRUE;
         }
         else
         {
            if( SCIPvarGetLbLocal(consdata->intvar) + 0.5 < (SCIP_Real) fixval )
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, consdata->intvar, (SCIP_Real) fixval, cons, (int)PROPRULE_1, TRUE, &infeasible, &tightened) );
               assert( tightened );
               assert( ! infeasible );
            }

            if( SCIPvarGetUbLocal(consdata->intvar) - 0.5 > (SCIP_Real) fixval )
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, consdata->intvar, (SCIP_Real) fixval, cons, (int)PROPRULE_1, TRUE, &infeasible, &tightened) );
               assert( tightened );
               assert( ! infeasible );
            }
            assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->intvar), SCIPvarGetUbLocal(consdata->intvar)));

            ++(*nfixedvars);
         }
      }

      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

      return SCIP_OKAY;
   }

   /* propagate w.r.t. integral variable */
   if ( consdata->intvar != NULL && !consdata->deleteintvar )
   {
      SCIP_Real newlb;
      SCIP_Real newub;
      int nonesmin;
      int nonesmax;

      assert( nfixedones + nfixedzeros < nvars );

      assert( SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(consdata->intvar)) );
      assert( SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(consdata->intvar)) );

      nonesmin = 2 * (int)(SCIPvarGetLbLocal(consdata->intvar) + 0.5) + consdata->rhs; /*lint !e713*/
      nonesmax = 2 * (int)(SCIPvarGetUbLocal(consdata->intvar) + 0.5) + consdata->rhs; /*lint !e713*/

      /* the number of possible variables that can get value 1 is less than the minimum bound */
      if ( nvars - nfixedzeros < nonesmin )
      {
         SCIPdebugMsg(scip, "constraint <%s>: at most %d variables can take value 1, but there should be at least %d.\n", SCIPconsGetName(cons), nvars - nfixedones, nonesmin);

         SCIP_CALL( analyzeConflict(scip, cons, NULL, PROPRULE_INTLB) );
         SCIP_CALL( SCIPresetConsAge(scip, cons) );

         *cutoff = TRUE;

         return SCIP_OKAY;
      }

      /* the number of variables that are fixed to 1 is larger than the maximum bound */
      if ( nfixedones > nonesmax )
      {
         SCIPdebugMsg(scip, "constraint <%s>: at least %d variables are fixed to 1, but there should be at most %d.\n", SCIPconsGetName(cons), nfixedones, nonesmax);

         SCIP_CALL( analyzeConflict(scip, cons, NULL, PROPRULE_INTUB) );
         SCIP_CALL( SCIPresetConsAge(scip, cons) );

         *cutoff = TRUE;

         return SCIP_OKAY;
      }

      /* compute new bounds on the integral variable */
      newlb = (SCIP_Real)((nfixedones + 1 - consdata->rhs) / 2); /*lint !e653*/
      newub = (SCIP_Real)((nvars - nfixedzeros - consdata->rhs) / 2); /*lint !e653*/

      /* new lower bound is better */
      if( newlb > SCIPvarGetLbLocal(consdata->intvar) + 0.5 )
      {
         SCIPdebugMsg(scip, "constraint <%s>: propagated lower bound of integral variable <%s> to %g\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->intvar), newlb);
         SCIP_CALL( SCIPinferVarLbCons(scip, consdata->intvar, newlb, cons, (int)PROPRULE_INTUB, TRUE, &infeasible, &tightened) );
         assert(tightened);
         assert(!infeasible);

         ++(*nchgbds);

         nonesmin = 2 * (int)(SCIPvarGetLbLocal(consdata->intvar) + 0.5) + consdata->rhs; /*lint !e713*/
      }

      /* new upper bound is better */
      if( newub < SCIPvarGetUbLocal(consdata->intvar) - 0.5 )
      {
         SCIPdebugMsg(scip, "constraint <%s>: propagated upper bound of integral variable <%s> to %g\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->intvar), newub);
         SCIP_CALL( SCIPinferVarUbCons(scip, consdata->intvar, newub, cons, (int)PROPRULE_INTLB, TRUE, &infeasible, &tightened) );
         assert(tightened);
         assert(!infeasible);

         ++(*nchgbds);

         nonesmax = 2 * (int)(SCIPvarGetUbLocal(consdata->intvar) + 0.5) + consdata->rhs; /*lint !e713*/
      }

      assert(nvars - nfixedzeros >= nonesmin);
      assert(nfixedones <= nonesmax);

      /* the number of variables that are free or fixed to 1 is exactly the minimum required -> fix free variables to 1 */
      if ( nvars - nfixedzeros == nonesmin )
      {
         SCIPdebugMsg(scip, "constraint <%s>: fix %d free variables to 1 to reach lower bound of %d\n", SCIPconsGetName(cons), nvars - nfixedzeros - nfixedones, nonesmin);

         for (i = 0; i < nvars; ++i)
         {
            if ( SCIPvarGetLbLocal(vars[i]) < 0.5 && SCIPvarGetUbLocal(vars[i]) > 0.5 )
            {
               SCIP_CALL( SCIPinferBinvarCons(scip, vars[i], TRUE, cons, (int)PROPRULE_INTLB, &infeasible, &tightened) );
               assert( !infeasible );
               assert( tightened );

               ++(*nfixedvars);
            }
         }
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );

         return SCIP_OKAY;
      }

      /* the number of variables that are fixed to 1 is exactly the maximum required -> fix free variables to 0 */
      if ( nfixedones == nonesmax )
      {
         SCIPdebugMsg(scip, "constraint <%s>: fix %d free variables to 0 to guarantee upper bound of %d\n", SCIPconsGetName(cons), nvars - nfixedzeros - nfixedones, nonesmax);

         for (i = 0; i < nvars; ++i)
         {
            if ( SCIPvarGetLbLocal(vars[i]) < 0.5 && SCIPvarGetUbLocal(vars[i]) > 0.5 )
            {
               SCIP_CALL( SCIPinferBinvarCons(scip, vars[i], FALSE, cons, (int)PROPRULE_INTUB, &infeasible, &tightened) );
               assert(!infeasible);
               assert(tightened);
               ++(*nfixedvars);
            }
         }
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );

         return SCIP_OKAY;
      }
   }

   /* switch to the new watched variables */
   SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, watchedvar1, watchedvar2) );

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** resolves a conflict on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rules (see propagateCons())
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
   assert(result != NULL);

   SCIPdebugMsg(scip, "resolving fixations according to rule %d\n", (int) proprule);

   SCIP_CALL( addConflictBounds(scip, cons, infervar, bdchgidx, proprule) );
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** try to use clique information to delete a part of the xor constraint or even fix variables */
static
SCIP_RETCODE cliquePresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  nchgcoefs,          /**< pointer to add up the number of deleted entries */
   int*                  ndelconss,          /**< pointer to add up the number of deleted constraints */
   int*                  naddconss,          /**< pointer to add up the number of added constraints */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if the node can be cut off */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Bool breaked;
   SCIP_Bool restart;
   int posnotinclq1;
   int posnotinclq2;
   int v;
   int v1;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(nchgcoefs != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);
   assert(cutoff != NULL);

   /* propagation can only be applied, if we know all operator variables */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   nvars = consdata->nvars;

   if( nvars < 3 )
      return SCIP_OKAY;

   /* we cannot perform this steps if the integer variables in not artificial */
   if( !consdata->deleteintvar )
      return SCIP_OKAY;

#if 0 /* try to evaluate if clique presolving should only be done multiple times when the constraint changed */
   if( !consdata->changed )
      return SCIP_OKAY;
#endif

   /* @todo: if clique information would have saved the type of the clique, like <= 1, or == 1 we could do more
    *        presolving like:
    *
    *        (xor(x1,x2,x3,x4) = 1 and clique(x1,x2) == 1)  =>  xor(x3,x4) = 0
    *        (xor(x1,x2,x3,x4) = 1 and clique(x1,x2,x3) == 1)  =>  (x4 = 0 and delete xor constraint)
    */

   /* 1. we have only clique information "<=", so we can check if all variables are in the same clique
    *
    * (xor(x1,x2,x3) = 1 and clique(x1,x2,x3) <= 1)  =>  (add set-partioning constraint x1 + x2 + x3 = 1 and delete old
    *                                                     xor-constraint)
    *
    * (xor(x1,x2,x3) = 0 and clique(x1,x2,x3) <= 1)  =>  (fix all variables x1 = x2 = x3 = 0 and delete old xor-
    *                                                     constraint)
    */

   /* 2. we have only clique information "<=", so we can check if all but one variable are in the same clique
    *
    * (xor(x1,x2,x3,x4) = 1 and clique(x1,x2,x3) <= 1)  =>  (add set-partioning constraint x1 + x2 + x3 + x4 = 1 and
    *                                                        delete old xor constraint)
    *
    * (xor(x1,x2,x3,x4) = 0 and clique(x1,x2,x3) <= 1)  =>  (add set-partioning constraint x1 + x2 + x3 + ~x4 = 1 and
    *                                                        delete old xor constraint)
    */

   posnotinclq1 = -1; /* index of variable that is possible not in the clique */
   posnotinclq2 = -1; /* index of variable that is possible not in the clique */
   breaked = FALSE;
   restart = FALSE;

   v = nvars - 2;
   while( v >= 0 )
   {
      SCIP_VAR* var;
      SCIP_VAR* var1;
      SCIP_Bool value;
      SCIP_Bool value1;

      assert(SCIPvarIsActive(vars[v]) || (SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[v]))));

      value = SCIPvarIsActive(vars[v]);

      if( !value )
         var = SCIPvarGetNegationVar(vars[v]);
      else
         var = vars[v];

      if( posnotinclq1 == v )
      {
         --v;
         continue;
      }

      for( v1 = v+1; v1 < nvars; ++v1 )
      {
         if( posnotinclq1 == v1 )
            continue;

         value1 = SCIPvarIsActive(vars[v1]);

         if( !value1 )
            var1 = SCIPvarGetNegationVar(vars[v1]);
         else
            var1 = vars[v1];

         if( !SCIPvarsHaveCommonClique(var, value, var1, value1, TRUE) )
         {
            /* if the position of the variable which is not in the clique with all other variables is not yet
             * initialized, than do now, one of both variables does not fit
             */
            if( posnotinclq1 == -1 )
            {
               posnotinclq1 = v;
               posnotinclq2 = v1;
            }
            else
            {
               /* no clique with exactly nvars-1 variables */
               if( restart || (posnotinclq2 != v && posnotinclq2 != v1) )
               {
                  breaked = TRUE;
                  break;
               }

               /* check the second variables for not fitting into the clique of (nvars - 1) variables */
               posnotinclq1 = posnotinclq2;
               restart = TRUE;
               v = nvars - 1;
            }

            break;
         }
         else
            assert(vars[v] != vars[v1]);
      }

      if( breaked )
         break;

      --v;
   }

   /* at least nvars-1 variables are in one clique */
   if( !breaked )
   {
      /* all variables are in one clique, case 1 */
      if( posnotinclq1 == -1 )
      {
         /* all variables of xor constraints <%s> (with rhs == 1) are in one clique, so create a setpartitioning
          * constraint with all variables and delete this xor-constraint */
         if( consdata->rhs )
         {
            SCIP_CONS* newcons;
            char consname[SCIP_MAXSTRLEN];

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_complete_clq", SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateConsSetpart(scip, &newcons, consname, nvars, vars,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

            SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIPdebugMsg(scip, "added a clique/setppc constraint <%s> \n", SCIPconsGetName(newcons));
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, newcons, NULL) ) );
            ++(*naddconss);

            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         }
         /* all variables of xor constraints <%s> (with rhs == 0) are in one clique, so fixed all variables to 0 */
         else
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            SCIPdebugMsg(scip, "all variables of xor constraints <%s> are in one clique, so fixed all variables to 0\n",
            SCIPconsGetName(cons));
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

            for( v = nvars - 1; v >= 0; --v )
            {
               SCIPdebugMsg(scip, "fixing variable <%s> to 0\n", SCIPvarGetName(vars[v]));
               SCIP_CALL( SCIPfixVar(scip, vars[v], 0.0, &infeasible, &fixed) );

               assert(infeasible || fixed);

               if( infeasible )
               {
                  *cutoff = infeasible;

                  return SCIP_OKAY;
               }
               else
                  ++(*nfixedvars);
            }
         }
      }
      /* all but one variable are in one clique, case 2 */
      else
      {
         SCIP_CONS* newcons;
         char consname[SCIP_MAXSTRLEN];

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_completed_clq", SCIPconsGetName(cons));

         /* complete clique by creating a set partioning constraint over all variables */

         /* if rhs == FALSE we need to exchange the variable not appaering in the clique with the negated variables */
         if( !consdata->rhs )
         {
            SCIP_CALL( SCIPcreateConsSetpart(scip, &newcons, consname, 0, NULL,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

            for( v = 0; v < nvars; ++v )
            {
               if( v == posnotinclq1 )
               {
                  SCIP_VAR* var;

                  SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &var) );
                  assert(var != NULL);

                  SCIP_CALL( SCIPaddCoefSetppc(scip, newcons, var) );
               }
               else
               {
                  SCIP_CALL( SCIPaddCoefSetppc(scip, newcons, vars[v]) );
               }
            }
         }
         /* if rhs == TRUE we can add all variables to the clique constraint directly */
         else
         {
            SCIP_CALL( SCIPcreateConsSetpart(scip, &newcons, consname, nvars, vars,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
         }

         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugMsg(scip, "added a clique/setppc constraint <%s> \n", SCIPconsGetName(newcons));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, newcons, NULL) ) );
         ++(*naddconss);

         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      }

      /* fix integer variable if it exists */
      if( consdata->intvar != NULL )
      {
         SCIP_Bool infeasible;
         SCIP_Bool fixed;

         SCIPdebugMsg(scip, "also fix the integer variable <%s> to 0\n", SCIPvarGetName(consdata->intvar));
         SCIP_CALL( SCIPfixVar(scip, consdata->intvar, 0.0, &infeasible, &fixed) );

         assert(infeasible || fixed || SCIPvarGetStatus(consdata->intvar) == SCIP_VARSTATUS_FIXED);

         if( infeasible )
         {
            *cutoff = infeasible;
            return SCIP_OKAY;
         }
         else if( fixed )
            ++(*nfixedvars);
      }

      /* delete old redundant xor-constraint */
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
   }

   return SCIP_OKAY;
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint
 *  accordingly; in contrast to preprocessConstraintPairs(), it uses a hash table
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  firstchange,        /**< pointer to store first changed constraint */
   int*                  nchgcoefs,          /**< pointer to add up the number of changed coefficients */
   int*                  naggrvars,          /**< pointer to add up the number of aggregated variables */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  naddconss,          /**< pointer to count number of added constraints */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if a cutoff was found */
   )
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int c;

   assert(conss != NULL);
   assert(ndelconss != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = nconss;
   hashtablesize = MAX(hashtablesize, HASHSIZE_XORCONS);

   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyXorcons, hashKeyEqXorcons, hashKeyValXorcons, (void*) scip) );

   /* check all constraints in the given set for redundancy */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata0;
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      /* get constraint handler data */
      conshdlr = SCIPconsGetHdlr(cons0);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* it can happen that during preprocessing some variables got aggregated and a constraint now has not active
       * variables inside so we need to remove them for sorting
       */
      /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
       * merge multiple entries of the same or negated variables
       */
      SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, cutoff) );
      if( *cutoff )
         goto TERMINATE;

      consdata0 = SCIPconsGetData(cons0);

      assert(consdata0 != NULL);

      /* sort the constraint */
      consdataSort(consdata0);
      assert(consdata0->sorted);

      /* get constraint from current hash table with same variables as cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));

      if( cons1 != NULL )
      {
         SCIP_CONSDATA* consdata1;

         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));

         consdata1 = SCIPconsGetData(cons1);

         assert(consdata1 != NULL);
         assert(consdata0->nvars >= 1 && consdata0->nvars == consdata1->nvars);

         assert(consdata0->sorted && consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);

         if( consdata0->rhs != consdata1->rhs )
         {
            *cutoff = TRUE;
            goto TERMINATE;
         }

         /* delete cons0 and update flags of cons1 s.t. nonredundant information doesn't get lost */
         SCIP_CALL( SCIPupdateConsFlags(scip, cons1, cons0) );
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         (*ndelconss)++;

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
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  naddconss,          /**< pointer to count number of added constraints */
   int*                  nchgcoefs           /**< pointer to add up the number of changed coefficients */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   SCIP_Bool cons0changed;
   int c;

   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);

   /* get constraint handler data */
   conshdlr = SCIPconsGetHdlr(cons0);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* it can happen that during preprocessing some variables got aggregated and a constraint now has not active
    * variables inside so we need to remove them for sorting
    */
   /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
    * merge multiple entries of the same or negated variables
    */
   SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   /* sort cons0 */
   consdataSort(consdata0);
   assert(consdata0->sorted);

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;
   consdata0->changed = FALSE;
   for( c = (cons0changed ? 0 : firstchange); c < chkind && !(*cutoff) && SCIPconsIsActive(cons0) && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      SCIP_VAR* singlevar0;
      SCIP_VAR* singlevar1;
      SCIP_Bool parity;
      SCIP_Bool cons0hastwoothervars;
      SCIP_Bool cons1hastwoothervars;
      SCIP_Bool aborted;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      SCIP_Bool redundant;
      SCIP_Bool aggregated;
      int v0;
      int v1;

      cons1 = conss[c];

      /* ignore inactive and modifiable constraints */
      if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
         continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      if( !consdata1->deleteintvar )
         continue;

      /* it can happen that during preprocessing some variables got aggregated and a constraint now has not active
       * variables inside so we need to remove them for sorting
       */
      /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
       * merge multiple entries of the same or negated variables
       */
      SCIP_CALL( applyFixings(scip, cons1, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, cutoff) );
      assert(consdata1 == SCIPconsGetData(cons1));
      if( *cutoff )
         return SCIP_OKAY;

      SCIPdebugMsg(scip, "preprocess xor constraint pair <%s>[chg:%u] and <%s>[chg:%u]\n",
         SCIPconsGetName(cons0), cons0changed, SCIPconsGetName(cons1), consdata1->changed);

      /* if both constraints were not changed since last round, we can ignore the pair */
      if( !cons0changed && !consdata1->changed )
         continue;

      /* applyFixings() led to an empty constraint */
      if( consdata1->nvars == 0 )
      {
         if( consdata1->rhs )
         {
            *cutoff = TRUE;
            break;
         }
         else
         {
            /* delete empty constraint */
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            ++(*ndelconss);

            continue;
         }
      }
      else if( consdata1->nvars == 1 )
      {
         /* fix remaining variable */
         SCIP_CALL( SCIPfixVar(scip, consdata1->vars[0], (SCIP_Real) consdata1->rhs, &infeasible, &fixed) );
         assert(!infeasible);

         if( fixed )
            ++(*nfixedvars);

         SCIP_CALL( SCIPdelCons(scip, cons1) );
         ++(*ndelconss);

         /* check for fixed variable in cons0 and remove it */
         SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, cutoff) );
         assert(!(*cutoff));

         /* sort cons0 */
         consdataSort(consdata0);
         assert(consdata0->sorted);

         continue;
      }
      else if( consdata1->nvars == 2 )
      {
         if( !(consdata1->rhs) )
         {
            /* aggregate var0 == var1 */
            SCIP_CALL( SCIPaggregateVars(scip, consdata1->vars[0], consdata1->vars[1], 1.0, -1.0, 0.0,
                  &infeasible, &redundant, &aggregated) );
         }
         else
         {
            /* aggregate var0 == 1 - var1 */
            SCIP_CALL( SCIPaggregateVars(scip, consdata1->vars[0], consdata1->vars[1], 1.0, 1.0, 1.0,
                  &infeasible, &redundant, &aggregated) );
         }
         assert(!infeasible);
         assert(redundant || SCIPdoNotAggr(scip));

         if( aggregated )
         {
            ++(*naggrvars);

            /* check for aggregated variable in cons0 and remove it */
            SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, cutoff) );
            if( *cutoff )
               return SCIP_OKAY;

            /* sort cons0 */
            consdataSort(consdata0);
            assert(consdata0->sorted);
         }

         if( redundant )
         {
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            ++(*ndelconss);
         }

         continue;
      }
      assert(consdata0->sorted);

      /* sort cons1 */
      consdataSort(consdata1);
      assert(consdata1->sorted);

      /* check whether
       *  (a) one problem variable set is a subset of the other, or
       *  (b) the problem variable sets are almost equal with only one variable in each constraint that is not
       *      member of the other
       */
      aborted = FALSE;
      parity = (consdata0->rhs ^ consdata1->rhs);
      cons0hastwoothervars = FALSE;
      cons1hastwoothervars = FALSE;
      singlevar0 = NULL;
      singlevar1 = NULL;
      v0 = 0;
      v1 = 0;
      while( (v0 < consdata0->nvars || v1 < consdata1->nvars) && !aborted )
      {
         int cmp;

         assert(v0 <= consdata0->nvars);
         assert(v1 <= consdata1->nvars);

         if( v0 == consdata0->nvars )
            cmp = +1;
         else if( v1 == consdata1->nvars )
            cmp = -1;
         else
            cmp = SCIPvarCompareActiveAndNegated(consdata0->vars[v0], consdata1->vars[v1]);

         switch( cmp )
         {
         case -1:
            /* variable doesn't appear in cons1 */
            assert(v0 < consdata0->nvars);
            if( singlevar0 == NULL )
            {
               singlevar0 = consdata0->vars[v0];
               if( cons1hastwoothervars )
                  aborted = TRUE;
            }
            else
            {
               cons0hastwoothervars = TRUE;
               if( singlevar1 != NULL )
                  aborted = TRUE;
            }
            v0++;
            break;

         case +1:
            /* variable doesn't appear in cons0 */
            assert(v1 < consdata1->nvars);
            if( singlevar1 == NULL )
            {
               singlevar1 = consdata1->vars[v1];
               if( cons0hastwoothervars )
                  aborted = TRUE;
            }
            else
            {
               cons1hastwoothervars = TRUE;
               if( singlevar0 != NULL )
                  aborted = TRUE;
            }
            v1++;
            break;

         case 0:
            /* variable appears in both constraints */
            assert(v0 < consdata0->nvars);
            assert(v1 < consdata1->nvars);
            assert(SCIPvarGetProbvar(consdata0->vars[v0]) == SCIPvarGetProbvar(consdata1->vars[v1]));
            if( consdata0->vars[v0] != consdata1->vars[v1] )
            {
               assert(SCIPvarGetNegatedVar(consdata0->vars[v0]) == consdata1->vars[v1]);
               parity = !parity;
            }
            v0++;
            v1++;
            break;

         default:
            SCIPerrorMessage("invalid comparison result\n");
            SCIPABORT();
            return SCIP_INVALIDDATA;  /*lint !e527*/
         }
      }

      /* check if a useful presolving is possible */
      if( (cons0hastwoothervars && singlevar1 != NULL) || (cons1hastwoothervars && singlevar0 != NULL) )
         continue;

      /* check if one problem variable set is a subset of the other */
      if( singlevar0 == NULL && singlevar1 == NULL )
      {
         /* both constraints are equal */
         if( !parity )
         {
            /* even parity: constraints are redundant */
            SCIPdebugMsg(scip, "xor constraints <%s> and <%s> are redundant: delete <%s>\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), SCIPconsGetName(cons1));
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugPrintCons(scip, cons1, NULL);

            /* delete cons1 and update flags of cons0 s.t. nonredundant information doesn't get lost */
            SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            (*ndelconss)++;

            if( consdata1->intvar != NULL )
            {
               /* need to update integer variable, consider the following case:
                * c1: xor(x1, x2, x3) = 1  (intvar1 = y1)
                * c2: xor(x1, x2, x3) = 1  (intvar0 = NULL or intvar0 = y0)
                *
                * if intvar0 = NULL we have to assign intvar0 = y1. otherwise, we have to ensure that y1 = y0 holds.
                * if aggregation is allowed, we can aggregate both variables. otherwise, we have to add a linear
                * constraint y1 - y0 = 0.
                */
               if( consdata0->intvar == NULL )
               {
                  SCIP_CALL( setIntvar(scip, cons0, consdata1->intvar) );
               }
               else
               {
                  /* aggregate integer variables */
                  SCIP_CALL( SCIPaggregateVars(scip, consdata1->intvar, consdata0->intvar, 1.0, -1.0, 0.0,
                        &infeasible, &redundant, &aggregated) );

                  *cutoff = *cutoff || infeasible;

                  if( aggregated )
                  {
                     (*naggrvars)++;
                     assert(SCIPvarIsActive(consdata0->intvar));
                  }
                  else
                  {
                     SCIP_CONS* newcons;
                     char consname[SCIP_MAXSTRLEN];
                     SCIP_VAR* newvars[2];
                     SCIP_Real vals[2];

                     newvars[0] = consdata1->intvar;
                     vals[0] = 1.0;
                     newvars[1] = consdata0->intvar;
                     vals[1] = -1.0;

                     (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "agg_%s", SCIPconsGetName(cons1));

                     SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, 2, newvars, vals, 0.0, 0.0,
                           SCIPconsIsInitial(cons1), SCIPconsIsSeparated(cons1), TRUE, /*SCIPconsIsEnforced(cons),*/
                           TRUE, TRUE, /*SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),*/
                           SCIPconsIsLocal(cons1), SCIPconsIsModifiable(cons1),
                           SCIPconsIsDynamic(cons1), SCIPconsIsRemovable(cons1), SCIPconsIsStickingAtNode(cons1)) );

                     SCIP_CALL( SCIPaddCons(scip, newcons) );
                     SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
                     ++(*naddconss);
                  }
               }
            }
         }
         else
         {
            /* odd parity: constraints are contradicting */
            SCIPdebugMsg(scip, "xor constraints <%s> and <%s> are contradicting\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugPrintCons(scip, cons1, NULL);
            *cutoff = TRUE;
         }
      }
      else if( singlevar1 == NULL )
      {
         /* cons1 is a subset of cons0 */
         if( !cons0hastwoothervars )
         {
            /* only one additional variable in cons0: fix this variable according to the parity */
            SCIPdebugMsg(scip, "xor constraints <%s> and <%s> yield sum %u == <%s>\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity, SCIPvarGetName(singlevar0));
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugPrintCons(scip, cons1, NULL);
            SCIP_CALL( SCIPfixVar(scip, singlevar0, parity ? 1.0 : 0.0, &infeasible, &fixed) );
            *cutoff = *cutoff || infeasible;
            if ( fixed )
               (*nfixedvars)++;

            /* delete cons1 and update flags of cons0 s.t. nonredundant information doesn't get lost */
            SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            (*ndelconss)++;
         }
         else
         {
            int v;

            /* more than one additional variable in cons0: add cons1 to cons0, thus eliminating the equal variables */
            SCIPdebugMsg(scip, "xor constraint <%s> is superset of <%s> with parity %u\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity);
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugPrintCons(scip, cons1, NULL);
            for( v = 0; v < consdata1->nvars; ++v )
            {
               SCIP_CALL( addCoef(scip, cons0, consdata1->vars[v]) );
            }

            SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, cutoff) );
            assert(SCIPconsGetData(cons0) == consdata0);
            assert(consdata0->nvars >= 2); /* at least the two "other" variables should remain in the constraint */
         }

         if( *cutoff )
            return SCIP_OKAY;

         consdataSort(consdata0);
         assert(consdata0->sorted);
      }
      else if( singlevar0 == NULL )
      {
         /* cons0 is a subset of cons1 */
         if( !cons1hastwoothervars )
         {
            /* only one additional variable in cons1: fix this variable according to the parity */
            SCIPdebugMsg(scip, "xor constraints <%s> and <%s> yield sum %u == <%s>\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity, SCIPvarGetName(singlevar1));
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugPrintCons(scip, cons1, NULL);
            SCIP_CALL( SCIPfixVar(scip, singlevar1, parity ? 1.0 : 0.0, &infeasible, &fixed) );
            assert(infeasible || fixed);
            *cutoff = *cutoff || infeasible;
            (*nfixedvars)++;

            /* delete cons1 and update flags of cons0 s.t. nonredundant information doesn't get lost */
            SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            (*ndelconss)++;
         }
         else
         {
            int v;

            /* more than one additional variable in cons1: add cons0 to cons1, thus eliminating the equal variables */
            SCIPdebugMsg(scip, "xor constraint <%s> is subset of <%s> with parity %u\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity);
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugPrintCons(scip, cons1, NULL);
            for( v = 0; v < consdata0->nvars; ++v )
            {
               SCIP_CALL( addCoef(scip, cons1, consdata0->vars[v]) );
            }
            SCIP_CALL( applyFixings(scip, cons1, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, cutoff) );
            assert(SCIPconsGetData(cons1) == consdata1);
            assert(consdata1->nvars >= 2); /* at least the two "other" variables should remain in the constraint */

            if( *cutoff )
               return SCIP_OKAY;

            consdataSort(consdata1);
            assert(consdata1->sorted);
         }
      }
      else
      {
         assert(!cons0hastwoothervars);
         assert(!cons1hastwoothervars);

         /* sum of constraints is parity == singlevar0 xor singlevar1: aggregate variables and delete cons1 */
         SCIPdebugMsg(scip, "xor constraints <%s> and <%s> yield sum %u == xor(<%s>,<%s>)\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity, SCIPvarGetName(singlevar0),
            SCIPvarGetName(singlevar1));
         if( !parity )
         {
            /* aggregate singlevar0 == singlevar1 */
            SCIP_CALL( SCIPaggregateVars(scip, singlevar1, singlevar0, 1.0, -1.0, 0.0,
                  &infeasible, &redundant, &aggregated) );
         }
         else
         {
            /* aggregate singlevar0 == 1-singlevar1 */
            SCIP_CALL( SCIPaggregateVars(scip, singlevar1, singlevar0, 1.0, 1.0, 1.0,
                  &infeasible, &redundant, &aggregated) );
         }
         assert(infeasible || redundant || SCIPdoNotAggr(scip));

         *cutoff = *cutoff || infeasible;
         if( aggregated )
            (*naggrvars)++;

         if( redundant )
         {
            /* delete cons1 and update flags of cons0 s.t. nonredundant information doesn't get lost */
            SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            (*ndelconss)++;

            if( consdata1->intvar != NULL )
            {
               if( consdata0->intvar == NULL )
               {
                  SCIP_CALL( setIntvar(scip, cons0, consdata0->intvar) );
               }
               else
               {
                  /* aggregate integer variables */
                  SCIP_CALL( SCIPaggregateVars(scip, consdata1->intvar, consdata0->intvar, 1.0, -1.0, 0.0,
                        &infeasible, &redundant, &aggregated) );

                  *cutoff = *cutoff || infeasible;
                  if( aggregated )
                     (*naggrvars)++;
               }
            }
         }

         if( !consdata0->sorted )
            consdataSort(consdata0);
         assert(consdata0->sorted);

#if 0
      /* if aggregation in the core of SCIP is not changed we do not need to call applyFixing, this would be the correct
       * way
       */
      /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
       * merge multiple entries of the same or negated variables
       */
      SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;
#endif

      }
   }

   return SCIP_OKAY;
}

/** creates and captures a xor constraint x_0 xor ... xor x_{k-1} = rhs with a given artificial integer variable for the
 *  linear relaxation
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
static
SCIP_RETCODE createConsXorIntvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_Bool             rhs,                /**< right hand side of the constraint */
   int                   nvars,              /**< number of operator variables in the constraint */
   SCIP_VAR**            vars,               /**< array with operator variables of constraint */
   SCIP_VAR*             intvar,             /**< artificial integer variable for linear relaxation */
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
   SCIP_CONSDATA* consdata;

   /* find the xor constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("xor constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, rhs, nvars, vars, intvar) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}



/*
 * Linear constraint upgrading
 */

/** tries to upgrade a linear constraint into an xor constraint
 *
 *  Assuming the only coefficients with an absolute value unequal to one is two in absolute value we can transform:
 *  \f[
 *    \begin{array}{ll}
 *                     & -\sum_{i \in I} x_i + \sum_{j \in J} x_j + a \cdot z = rhs \\
 *     \Leftrightarrow & \sum_{i \in I} ~x_i + \sum_j x_{j \in J} + a \cdot z = rhs + |I| \\
 *     \Leftrightarrow & \sum_{i \in I} ~x_i + \sum_j x_{j \in J} - 2 \cdot y = (rhs + |I|) mod 2 \\
 *                     & z \in [lb_z,ub_z] \Rightarrow y \in (\left\lfloor \frac{rhs + |I|}{2} \right\rfloor + (a = -2\; ?\; lb_z : -ub_z), \left\lfloor \frac{rhs + |I|}{2} \right\rfloor + (a = -2\; ?\; ub_z : -lb_z) ) \\
 *                     & z = (a = -2\; ?\; y - \left\lfloor \frac{rhs + |I|}{2} \right\rfloor\; :\; -y + \left\lfloor \frac{rhs + |I|}{2} \right\rfloor)
 *    \end{array}
 *  \f]
 */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdXor)
{  /*lint --e{715}*/
   assert( upgdcons != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 0 );
   assert( ! SCIPconsIsModifiable(cons) );

   /* check, if linear constraint can be upgraded to xor constraint */
   /* @todo also applicable if the integer variable has a coefficient different from 2, e.g. a coefficient like 0.5 then
    *       we could generate a new integer variable aggregated to the old one, possibly the constraint was then
    *       normalized and all binary variables have coefficients of 2.0, if the coefficient is 4 then we need holes ...
    */
   if( integral && nposcont + nnegcont == 0 && nposbin + nnegbin + nposimplbin + nnegimplbin >= nvars-1 && ncoeffspone + ncoeffsnone == nvars-1 && ncoeffspint + ncoeffsnint == 1 )
   {
      assert( ncoeffspfrac + ncoeffsnfrac == 0 );

      if ( SCIPisEQ(scip, lhs, rhs) && SCIPisIntegral(scip, lhs) )
      {
         SCIP_VAR** xorvars;
         SCIP_VAR* parityvar = NULL;
         SCIP_Bool postwo = FALSE;
         int cnt = 0;
         int j;

         SCIP_CALL( SCIPallocBufferArray(scip, &xorvars, nvars) );

         /* check parity of constraints */
         for( j = nvars - 1; j >= 0; --j )
         {
            if( SCIPisEQ(scip, REALABS(vals[j]), 2.0) )
            {
               parityvar = vars[j];
               postwo = (vals[j] > 0.0);
            }
            else if( !SCIPisEQ(scip, REALABS(vals[j]), 1.0) )
               break;
            else
            {
               /* exit if variable is not binary or implicit binary */
               if ( ! SCIPvarIsBinary(vars[j]) )
               {
                  parityvar = NULL;
                  break;
               }

               /* need negated variables for correct propagation to the integer variable */
               if( vals[j] < 0.0 )
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, vars[j], &(xorvars[cnt])) );
                  assert(xorvars[cnt] != NULL);
               }
               else
                  xorvars[cnt] = vars[j];
               ++cnt;
            }
         }

         if( parityvar != NULL )
         {
            assert(cnt == nvars - 1);

            /* check whether parity variable is present only in this constraint */
            if ( SCIPvarGetNLocksDown(parityvar) <= 1 && SCIPvarGetNLocksUp(parityvar) <= 1 )
            {
               SCIP_VAR* intvar;
               SCIP_Bool rhsparity;
               SCIP_Bool neednew;
               int intrhs;

               SCIPdebugMsg(scip, "upgrading constraint <%s> to an XOR constraint\n", SCIPconsGetName(cons));

               /* adjust the side, since we negated all binary variables with -1.0 as a coefficient */
               rhs += ncoeffsnone;

               intrhs = (int) SCIPfloor(scip, rhs);
               rhsparity = ((SCIP_Bool) (intrhs % 2)); /*lint !e571*/
               neednew = (intrhs != 1 && intrhs != 0);

               /* check if we can use the parity variable as integer variable of the XOR constraint or do we need to
                * create a new variable
                */
               if( neednew )
               {
                  char varname[SCIP_MAXSTRLEN];
                  SCIP_Real lb;
                  SCIP_Real ub;
                  SCIP_Bool isbinary;
                  SCIP_Bool infeasible;
                  SCIP_Bool redundant;
                  SCIP_Bool aggregated;
                  int intrhshalfed;

                  intrhshalfed = intrhs / 2;

                  if( postwo )
                  {
                     lb = intrhshalfed - SCIPvarGetUbGlobal(parityvar);
                     ub = intrhshalfed - SCIPvarGetLbGlobal(parityvar);
                  }
                  else
                  {
                     lb = intrhshalfed + SCIPvarGetLbGlobal(parityvar);
                     ub = intrhshalfed + SCIPvarGetUbGlobal(parityvar);
                  }
                  assert(SCIPisFeasLE(scip, lb, ub));
                  assert(SCIPisFeasIntegral(scip, lb));
                  assert(SCIPisFeasIntegral(scip, ub));

                  /* adjust bounds to be more integral */
                  lb = SCIPfeasFloor(scip, lb);
                  ub = SCIPfeasFloor(scip, ub);

                  isbinary = (SCIPisZero(scip, lb) && SCIPisEQ(scip, ub, 1.0));

                  /* you must not create an artificial integer variable if parity variable is already binary */
                  if( SCIPvarIsBinary(parityvar) && !isbinary )
                  {
                     SCIPfreeBufferArray(scip, &xorvars);
                     return SCIP_OKAY;
                  }

                  (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_xor_upgr", SCIPvarGetName(parityvar));
                  SCIP_CALL( SCIPcreateVar(scip, &intvar, varname, lb, ub, 0.0,
                        isbinary ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER,
                        SCIPvarIsInitial(parityvar), SCIPvarIsRemovable(parityvar), NULL, NULL, NULL, NULL, NULL) );
                  SCIP_CALL( SCIPaddVar(scip, intvar) );

                  SCIPdebugMsg(scip, "created new variable for aggregation:");
                  SCIPdebug( SCIPprintVar(scip, intvar, NULL) );

                  SCIP_CALL( SCIPaggregateVars(scip, parityvar, intvar, 1.0, postwo ? 1.0 : -1.0,
                        (SCIP_Real) (postwo ? intrhshalfed : -intrhshalfed), &infeasible, &redundant, &aggregated) );
                  assert(!infeasible);

                  /* maybe aggregation was forbidden, than we cannot upgrade this constraint */
                  if( !aggregated )
                  {
                     SCIPfreeBufferArray(scip, &xorvars);
                     return SCIP_OKAY;
                  }

                  assert(redundant);
                  assert(SCIPvarIsActive(intvar));

                  SCIPdebugMsg(scip, "aggregated: %s = %g * %s + %g\n", SCIPvarGetName(parityvar),
                     SCIPvarGetAggrScalar(parityvar), SCIPvarGetName(SCIPvarGetAggrVar(parityvar)),
                     SCIPvarGetAggrConstant(parityvar));
               }
               else
                  intvar = parityvar;

               assert(intvar != NULL);

               SCIPdebugMsg(scip, "upgrading constraint <%s> to XOR constraint\n", SCIPconsGetName(cons));

               SCIP_CALL( createConsXorIntvar(scip, upgdcons, SCIPconsGetName(cons), rhsparity, nvars - 1, xorvars, intvar,
                        SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                        SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                        SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                        SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               SCIPdebugPrintCons(scip, *upgdcons, NULL);

               if( neednew )
               {
                  assert(intvar != parityvar);
                  SCIP_CALL( SCIPreleaseVar(scip, &intvar) );
               }
            }
         }

         SCIPfreeBufferArray(scip, &xorvars);
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyXor)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrXor(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeXor)
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
SCIP_DECL_CONSEXITSOL(consExitsolXor)
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
SCIP_DECL_CONSDELETE(consDeleteXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_INITPRESOLVE )
   {
      int v;

      for( v = (*consdata)->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)(*consdata), -1) );
      }
   }

   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->nvars >= 1);
   assert(sourcedata->vars != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->rhs, sourcedata->nvars, sourcedata->vars, sourcedata->intvar) );

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
SCIP_DECL_CONSINITLP(consInitlpXor)
{  /*lint --e{715}*/
   int i;

   assert(infeasible != NULL);

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
SCIP_DECL_CONSSEPALP(consSepalpXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool separated;
   SCIP_Bool cutoff;
   int c;

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, conshdlrdata->separateparity, &separated, &cutoff) );
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
SCIP_DECL_CONSSEPASOL(consSepasolXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool separated;
   SCIP_Bool cutoff;
   int c;

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, conshdlrdata->separateparity, &separated, &cutoff) );
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
SCIP_DECL_CONSENFOLP(consEnfolpXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;
   SCIP_Bool cutoff;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, FALSE, &violated) );
      if( violated )
      {
         SCIP_Bool separated;

         SCIP_CALL( separateCons(scip, conss[i], NULL, conshdlrdata->separateparity, &separated, &cutoff) );
         if ( cutoff )
            *result = SCIP_CUTOFF;
         else
         {
            assert(separated); /* because the solution is integral, the separation always finds a cut */
            *result = SCIP_SEPARATED;
         }
         return SCIP_OKAY;
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;
   SCIP_Bool cutoff;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, FALSE, &violated) );
      if( violated )
      {
         SCIP_Bool separated;

         SCIP_CALL( separateCons(scip, conss[i], sol, conshdlrdata->separateparity, &separated, &cutoff) );
         if ( cutoff )
            *result = SCIP_CUTOFF;
         else
         {
            assert(separated); /* because the solution is integral, the separation always finds a cut */
            *result = SCIP_SEPARATED;
         }
         return SCIP_OKAY;
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsXor)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, TRUE, &violated) );
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
SCIP_DECL_CONSCHECK(consCheckXor)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   *result = SCIP_FEASIBLE;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss && (*result == SCIP_FEASIBLE || completely); i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, checklprows, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            int v;
            int sum = 0;
            SCIP_CONSDATA* consdata;

            consdata = SCIPconsGetData(conss[i]);
            assert( consdata != NULL );

            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );

            for( v = 0; v < consdata->nvars; ++v )
            {
               if( SCIPgetSolVal(scip, sol, consdata->vars[v]) > 0.5 )
                  sum++;
            }

            if( consdata->intvar != NULL )
            {
               SCIPinfoMessage(scip, NULL, ";\nviolation: %d operands are set to TRUE but integer variable has value of %g\n", SCIPgetSolVal(scip, sol, consdata->intvar));
            }
            else
            {
               SCIPinfoMessage(scip, NULL, ";\nviolation: %d operands are set to TRUE\n", sum );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int nfixedvars;
   int nchgbds;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   nfixedvars = 0;
   nchgbds = 0;

   /* propagate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &nfixedvars, &nchgbds) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 || nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
   {
      *result = SCIP_DIDNOTFIND;
      if ( ! SCIPinProbing(scip) )
      {
         int depth;
         int freq;

         depth = SCIPgetDepth(scip);
         freq = conshdlrdata->gausspropfreq;
         if ( (depth == 0 && freq == 0) || (freq > 0 && depth % freq == 0) )
         {
            /* take useful constraints only - might improve success rate to take all */
            SCIP_CALL( checkSystemGF2(scip, conss, nusefulconss, NULL, result) );
         }
      }
   }

   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;
   int v;

   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* catch all variable event for deleted variables, which is only used in presolving */
   for( c = nconss - 1; c >= 0; --c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)consdata, NULL) );
      }
   }

   return SCIP_OKAY;
}

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;
   int v;

   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* drop all variable event for deleted variables, which was only used in presolving */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPconsIsDeleted(conss[c]) )
      {
         for( v = 0; v < consdata->nvars; ++v )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
                  (SCIP_EVENTDATA*)consdata, -1) );
         }
      }
   }

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;
   int oldnfixedvars;
   int oldnchgbds;
   int oldnaggrvars;
   int oldndelconss;
   int oldnchgcoefs;
   int firstchange;
   int c;

   assert(result != NULL);

   oldnfixedvars = *nfixedvars;
   oldnchgbds = *nchgbds;
   oldnaggrvars = *naggrvars;
   oldndelconss = *ndelconss;
   oldnchgcoefs = *nchgcoefs;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   cutoff = FALSE;
   firstchange = INT_MAX;
   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
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

      /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
       * merge multiple entries of the same or negated variables
       */
      SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr, nchgcoefs, naggrvars, naddconss, &cutoff) );

      if( cutoff )
         break;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars, nchgbds) );

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 2); /* otherwise, propagateCons() has deleted the constraint */

         /* if only two variables are left, both have to be equal or opposite, depending on the rhs */
         if( consdata->nvars == 2 )
         {
            SCIPdebugMsg(scip, "xor constraint <%s> has only two unfixed variables, rhs=%u\n",
               SCIPconsGetName(cons), consdata->rhs);

            assert(consdata->vars != NULL);
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[1]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[1]), 1.0));

            if( !consdata->rhs )
            {
               /* aggregate variables: vars[0] - vars[1] == 0 */
               SCIPdebugMsg(scip, " -> aggregate <%s> == <%s>\n", SCIPvarGetName(consdata->vars[0]),
                  SCIPvarGetName(consdata->vars[1]));
               SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], 1.0, -1.0, 0.0,
                     &cutoff, &redundant, &aggregated) );
            }
            else
            {
               /* aggregate variables: vars[0] + vars[1] == 1 */
               SCIPdebugMsg(scip, " -> aggregate <%s> == 1 - <%s>\n", SCIPvarGetName(consdata->vars[0]),
                  SCIPvarGetName(consdata->vars[1]));
               SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], 1.0, 1.0, 1.0,
                     &cutoff, &redundant, &aggregated) );
            }
            assert(redundant || SCIPdoNotAggr(scip));

            if( aggregated )
            {
               assert(redundant);
               (*naggrvars)++;
            }

            /* the constraint can be deleted if the intvar is fixed or NULL */
            if( redundant )
            {
               SCIP_Bool fixedintvar;

               fixedintvar = consdata->intvar == NULL ? TRUE : SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->intvar), SCIPvarGetUbGlobal(consdata->intvar));

               if( fixedintvar )
               {
                  /* if the integer variable is an original variable, i.e.,
                   * consdata->deleteintvar == FALSE then the following
                   * must hold:
                   *
                   *   if consdata->rhs == 1 then the integer variable needs
                   *   to be fixed to zero, otherwise the constraint is
                   *   infeasible,
                   *
                   *   if consdata->rhs == 0 then the integer variable needs
                   *   to be aggregated to one of the binary variables
                   */
                  assert(consdata->deleteintvar || (consdata->rhs && SCIPvarGetLbGlobal(consdata->intvar) < 0.5));

                  /* delete constraint */
                  SCIP_CALL( SCIPdelCons(scip, cons) );
                  (*ndelconss)++;
               }
            }
         }
         else if( (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 )
         {
            /* try to use clique information to upgrade the constraint to a set-partitioning constraint or fix
             * variables
             */
            SCIP_CALL( cliquePresolve(scip, cons, nfixedvars, nchgcoefs, ndelconss, naddconss, &cutoff) );
         }
      }
   }

   /* process pairs of constraints: check them for equal operands;
    * only apply this expensive procedure, if the single constraint preprocessing did not find any reductions
    */
   if( !cutoff && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 && SCIPisPresolveFinished(scip) )
   {
      if( firstchange < nconss && conshdlrdata->presolusehashing )
      {
         /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
         SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, nchgcoefs, naggrvars, ndelconss, naddconss, &cutoff) );
      }
      if( conshdlrdata->presolpairwise )
      {
         SCIP_Longint npaircomparisons;
         int lastndelconss;
         npaircomparisons = 0;
         lastndelconss = *ndelconss;

         for( c = firstchange; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
         {
            if( SCIPconsIsActive(conss[c]) && !SCIPconsIsModifiable(conss[c]) )
            {
               npaircomparisons += (SCIPconsGetData(conss[c])->changed) ? (SCIP_Longint) c : ((SCIP_Longint) c - (SCIP_Longint) firstchange);

               SCIP_CALL( preprocessConstraintPairs(scip, conss, firstchange, c,
                     &cutoff, nfixedvars, naggrvars, ndelconss, naddconss, nchgcoefs) );

               if( npaircomparisons > NMINCOMPARISONS )
               {
                  if( ((SCIP_Real) (*ndelconss - lastndelconss)) / ((SCIP_Real) npaircomparisons) < MINGAINPERNMINCOMPARISONS )
                     break;
                  lastndelconss = *ndelconss;
                  npaircomparisons = 0;
               }
            }
         }
      }
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *nchgbds > oldnchgbds || *naggrvars > oldnaggrvars
      || *ndelconss > oldndelconss || *nchgcoefs > oldnchgcoefs )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   /* add extended formulation at the end of presolving if required */
   if ( conshdlrdata->addextendedform && *result == SCIP_DIDNOTFIND && SCIPisPresolveFinished(scip) )
   {
      for (c = 0; c < nconss && ! SCIPisStopped(scip); ++c)
      {
         int naddedconss = 0;

         cons = conss[c];
         assert(cons != NULL);
         consdata = SCIPconsGetData(cons);
         assert(consdata != NULL);

         if ( consdata->extvars != NULL )
            break;

         if ( conshdlrdata->addflowextended )
         {
            SCIP_CALL( addExtendedFlowFormulation(scip, cons, &naddedconss) );
         }
         else
         {
            SCIP_CALL( addExtendedAsymmetricFormulation(scip, cons, &naddedconss) );
         }
         (*naddconss) += naddedconss;
      }
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropXor)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* external variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   /* internal variable */
   if( consdata->intvar != NULL )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->intvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintXor)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file, FALSE) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** targetvars;
   SCIP_VAR* intvar;
   SCIP_VAR* targetintvar;
   const char* consname;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   (*valid) = TRUE;

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get variables and coefficients of the source constraint */
   sourcevars = sourceconsdata->vars;
   nvars = sourceconsdata->nvars;
   intvar = sourceconsdata->intvar;
   targetintvar = NULL;

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   if( nvars == 0 )
   {
      if( intvar != NULL )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, intvar, &targetintvar, varmap, consmap, global, valid) );
         assert(!(*valid) || targetintvar != NULL);

         SCIPdebugMsg(scip, "Copied integral variable <%s> (bounds: [%g,%g])\n", SCIPvarGetName(targetintvar),
            global ? SCIPvarGetLbGlobal(intvar) : SCIPvarGetLbLocal(intvar),
            global ? SCIPvarGetUbGlobal(intvar) : SCIPvarGetUbLocal(intvar));
      }

      if( *valid )
      {
         SCIP_CALL( createConsXorIntvar(scip, cons, consname, SCIPgetRhsXor(sourcescip, sourcecons), 0, NULL,
               targetintvar, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable,
               stickingatnode) );
      }

      return SCIP_OKAY;
   }

   /* duplicate variable array */
   SCIP_CALL( SCIPallocBufferArray(scip, &targetvars, nvars) );

   /* map variables of the source constraint to variables of the target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &targetvars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || targetvars[v] != NULL);
   }

   /* map artificial relaxation variable of the source constraint to variable of the target SCIP */
   if( *valid && intvar != NULL )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, intvar, &targetintvar, varmap, consmap, global, valid) );
      assert(!(*valid) || targetintvar != NULL);

      SCIPdebugMsg(scip, "Copied integral variable <%s> (bounds: [%g,%g])\n", SCIPvarGetName(targetintvar),
         global ? SCIPvarGetLbGlobal(intvar) : SCIPvarGetLbLocal(intvar),
         global ? SCIPvarGetUbGlobal(intvar) : SCIPvarGetUbLocal(intvar));
   }

   /* only create the target constraints, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( createConsXorIntvar(scip, cons, consname, SCIPgetRhsXor(sourcescip, sourcecons), nvars, targetvars,
            targetintvar, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable,
            stickingatnode) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &targetvars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseXor)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   char* endptr;
   int requiredsize;
   int varssize;
   int nvars;

   SCIPdebugMsg(scip, "parse <%s> as xor constraint\n", str);

   varssize = 100;
   nvars = 0;

   /* allocate buffer array for variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );

   /* parse string */
   SCIP_CALL( SCIPparseVarsList(scip, str, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );

   if( *success )
   {
      SCIP_Real rhs;

      /* check if the size of the variable array was big enough */
      if( varssize < requiredsize )
      {
         /* reallocate memory */
         varssize = requiredsize;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, varssize) );

         /* parse string again with the correct size of the variable array */
         SCIP_CALL( SCIPparseVarsList(scip, str, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );
      }

      assert(*success);
      assert(varssize >= requiredsize);

      SCIPdebugMsg(scip, "successfully parsed %d variables\n", nvars);

      str = endptr;

      /* search for the equal symbol */
      while( *str != '=' && *str != '\0' )
         str++;

      /* if the string end has been reached without finding the '=' */
      if ( *str == '\0' )
      {
         SCIPerrorMessage("Could not find terminating '='.\n");
         *success = FALSE;
      }
      else
      {
         /* skip '=' character */
         ++str;

         if( SCIPstrToRealValue(str, &rhs, &endptr) )
         {
            SCIP_VAR* intvar = NULL;

            assert(SCIPisZero(scip, rhs) || SCIPisEQ(scip, rhs, 1.0));

            str = endptr;

            /* skip white spaces */
            while( *str == ' ' || *str == '\t' )
               str++;

            /* check for integer variable, should look like (intvar = var) */
            if( *str == '(' )
            {
               str++;
               while( *str != '=' && *str != '\0' )
                  str++;

               if( *str != '=' )
               {
                  SCIPerrorMessage("Parsing integer variable of XOR constraint\n");
                  *success = FALSE;
                  goto TERMINATE;
               }

               str++;
               /* skip white spaces */
               while( *str == ' ' || *str == '\t' )
                  str++;

               /* parse variable name */
               SCIP_CALL( SCIPparseVarName(scip, str, &intvar, &endptr) );

               if( intvar == NULL )
               {
                  SCIPdebugMsg(scip, "variable with name <%s> does not exist\n", SCIPvarGetName(intvar));
                  (*success) = FALSE;
                  goto TERMINATE;
               }
               str = endptr;

               /* skip last ')' */
               while( *str != ')' && *str != '\0' )
                  str++;
            }

            if( intvar != NULL )
            {
               /* create or constraint */
               SCIP_CALL( createConsXorIntvar(scip, cons, name, (rhs > 0.5 ? TRUE : FALSE), nvars, vars, intvar,
                     initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
            }
            else
            {
               /* create or constraint */
               SCIP_CALL( SCIPcreateConsXor(scip, cons, name, (rhs > 0.5 ? TRUE : FALSE), nvars, vars,
                     initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
            }

            SCIPdebugPrintCons(scip, *cons, NULL);
         }
         else
            *success = FALSE;
      }
   }

 TERMINATE:
   /* free variable buffer */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int nintvar = 0;
   int cnt;
   int j;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if ( consdata->intvar != NULL )
      nintvar = 1;

   if ( varssize < consdata->nvars + nintvar + consdata->nextvars )
      (*success) = FALSE;
   else
   {
      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);

      if ( consdata->intvar != NULL )
         vars[consdata->nvars] = consdata->intvar;

      if ( consdata->nextvars > 0 )
      {
         assert( consdata->extvars != NULL );
         cnt = consdata->nvars + nintvar;
         for (j = 0; j < consdata->extvarssize; ++j)
         {
            if ( consdata->extvars[j] != NULL )
               vars[cnt++] = consdata->extvars[j];
         }
         assert( cnt == consdata->nvars + nintvar + consdata->nextvars );
      }

      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variable (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->intvar == NULL )
      (*nvars) = consdata->nvars + consdata->nextvars;
   else
      (*nvars) = consdata->nvars + 1 + consdata->nextvars;

   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_VARFIXED )
   {
      /* we only catch this event in presolving stage */
      assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
      assert(SCIPeventGetVar(event) != NULL);

      consdata->sorted = FALSE;
   }

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for xor constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrXor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;

   /* create event handler for events on variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecXor, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpXor, consEnfopsXor, consCheckXor, consLockXor,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyXor, consCopyXor) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteXor) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolXor) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeXor) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsXor) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsXor) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpXor) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseXor) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreXor) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreXor) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolXor, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintXor) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropXor, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropXor) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpXor, consSepasolXor, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransXor) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxXor) );

   if ( SCIPfindConshdlr(scip, "linear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdXor, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }

   /* add xor constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/xor/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/xor/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance?",
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/xor/addextendedform",
         "should the extended formulation be added in presolving?",
         &conshdlrdata->addextendedform, TRUE, DEFAULT_ADDEXTENDEDFORM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/xor/addflowextended",
         "should the extended flow formulation be added (nonsymmetric formulation otherwise)?",
         &conshdlrdata->addflowextended, TRUE, DEFAULT_ADDFLOWEXTENDED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/xor/separateparity",
         "should parity inequalities be separated?",
         &conshdlrdata->separateparity, TRUE, DEFAULT_SEPARATEPARITY, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/xor/gausspropfreq",
         "frequency for applying the Gauss propagator",
         &conshdlrdata->gausspropfreq, TRUE, DEFAULT_GAUSSPROPFREQ, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a xor constraint x_0 xor ... xor x_{k-1} = rhs
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_Bool             rhs,                /**< right hand side of the constraint */
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
   SCIP_CONSDATA* consdata;

   /* find the xor constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("xor constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, rhs, nvars, vars, NULL) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a xor constraint x_0 xor ... xor x_{k-1} = rhs
 *  with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_Bool             rhs,                /**< right hand side of the constraint */
   int                   nvars,              /**< number of operator variables in the constraint */
   SCIP_VAR**            vars                /**< array with operator variables of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsXor(scip,cons, name, rhs, nvars, vars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** gets number of variables in xor constraint */
int SCIPgetNVarsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in xor constraint */
SCIP_VAR** SCIPgetVarsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets integer variable in xor constraint */
SCIP_VAR* SCIPgetIntVarXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->intvar;
}

/** gets the right hand side of the xor constraint */
SCIP_Bool SCIPgetRhsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}
