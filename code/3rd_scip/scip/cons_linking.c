/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_linking.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for linking constraints
 * @author Stefan Heinz
 * @author Jens Schulz
 *
 * The constraints handler stores linking constraints between a linking variable (integer or continuous) and an array of binary variables. Such
 * a linking constraint has the form:
 *
 * linkvar = sum_{i=1}^n {vals[i] * binvars[i]}
 *
 * with the additional side condition that exactly one binary variable has to be one (set partitioning condition).
 *
 * This constraint can be created only with the linking variable if it is an integer variable. In this case the binary variables are only created on
 * demand. That is, whenever someone asks for the binary variables. Therefore, such constraints can be used to get a
 * "binary representation" of the domain of the linking variable which will be dynamically created.
 *
 *
 * @todo add pairwise comparison of constraints in presolving (fast hash table version and complete pairwise comparison)
 * @todo in case the integer variable is set to lower or upper bound it follows that only the corresponding binary
 * variable has a positive value which is one, this can be used to fasten the checking routine
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/cons_linking.h"
#include "scip/cons_setppc.h"
#include "scip/pub_cons.h"
#include "scip/pub_event.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_var.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_cut.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include <ctype.h>
#include <string.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "linking"
#define CONSHDLR_DESC          "linking constraint x = sum_{i=1}^{n} c_i*y_i, y1+...+yn = 1, x real, y's binary"

#define EVENTHDLR_NAME         "linking"
#define EVENTHDLR_DESC         "event handler for linking constraints"

#define CONSHDLR_SEPAPRIORITY    750000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2050000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -750000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation, propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */


#define HASHSIZE_BINVARSCONS        500 /**< minimal size of hash table in linking constraint handler */
#define DEFAULT_LINEARIZE         FALSE /**< should the linking constraint be linearize after the binary variable are created */

/*
 * Data structures
 */

/** constraint data for linking constraints */
struct SCIP_ConsData
{
   SCIP_VAR*             linkvar;            /**< continuous variable which is linked */
   SCIP_VAR**            binvars;            /**< binary variables */
   SCIP_Real*            vals;               /**< coefficients */
   SCIP_ROW*             row1;               /**< LP row for the linking itself */
   SCIP_ROW*             row2;               /**< LP row ensuring the set partitioning condition of the binary variables */
   SCIP_NLROW*           nlrow1;             /**< NLP row for the linking itself */
   SCIP_NLROW*           nlrow2;             /**< NLP row ensuring the set partitioning condition of the binary variables */
   int                   nbinvars;           /**< number of binary variables */
   int                   sizebinvars;        /**< size of the binary variable array */
   int                   nfixedzeros;        /**< current number of variables fixed to zero in the constraint */
   int                   nfixedones;         /**< current number of variables fixed to one in the constraint */
   int                   firstnonfixed;      /**< index of first locally non-fixed binary variable in binvars array */
   int                   lastnonfixed;       /**< index of last locally non-fixed binary variable in binvars array */
   unsigned int          cliqueadded:1;      /**< was the set partitioning condition already added as clique? */
   unsigned int          sorted:1;           /**< are the coefficients of the binary variables are sorted in non-decreasing order */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events on binary variables */
   SCIP_HASHMAP*         varmap;             /**< hash map mapping a linking variable to its linking constraint */
   SCIP_Bool             linearize;          /**< should the linking constraint be linearize after the binary variable are created */
};

/*
 * Local methods
 */

/** returns for a given linking variable the corresponding hash map key */
static
void* getHashmapKey(
   SCIP_VAR*             var                 /**< variable to get the hash map key for */
   )
{
   /* return the unique variable index + 1 */
   return (void*)(size_t)(SCIPvarGetIndex(var) + 1); /*lint !e571 !e776*/
}

/* sort binary variable in non-decreasing order w.r.t. coefficients */
static
void consdataSort(
   SCIP_CONSDATA*        consdata            /**< linking constraint data */
   )
{
   if( consdata->sorted )
      return;

   /* sort binary variable in non-decreasing order w.r.t. coefficients */
   SCIPsortRealPtr(consdata->vals, (void**)consdata->binvars, consdata->nbinvars);

   consdata->sorted = TRUE;
}


/** installs rounding locks for the binary variables in the given linking constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR**            binvars,            /**< binary variables  */
   int                   nbinvars            /**< number of binary variables */
   )
{
   int b;

   for( b = 0; b < nbinvars; ++b )
   {
      SCIP_CALL( SCIPlockVarCons(scip, binvars[b], cons, TRUE, TRUE) );
   }

   return SCIP_OKAY;
}

/** creates constraint handler data for the linking constraint handler */
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

   /* create hash map */
   (*conshdlrdata)->varmap = NULL;

   /* set event handler for bound change events on binary variables */
   (*conshdlrdata)->eventhdlr = eventhdlr;

   return SCIP_OKAY;
}

/** frees constraint handler data for linking constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   /* free hash map */
   if( (*conshdlrdata)->varmap != NULL )
      SCIPhashmapFree(&(*conshdlrdata)->varmap);

   /* free memory of constraint handler data */
   SCIPfreeBlockMemory(scip, conshdlrdata);
}

/** prints linking constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** binvars;
   SCIP_VAR* linkvar;
   int nbinvars;

   assert(scip != NULL);
   assert(consdata != NULL);

   linkvar = consdata->linkvar;
   binvars = consdata->binvars;
   nbinvars = consdata->nbinvars;

   assert(linkvar != NULL);
   assert(binvars != NULL || nbinvars == 0);

   /* print linking variable */
   SCIP_CALL( SCIPwriteVarName(scip, file, linkvar, FALSE) );

   SCIPinfoMessage(scip, file, " = ");

   if( nbinvars == 0 )
   {
      SCIPinfoMessage(scip, file, " no binary variables yet");
   }
   else
   {
      assert(binvars != NULL);

      SCIP_CALL( SCIPwriteVarsLinearsum(scip, file, binvars, consdata->vals, nbinvars, FALSE) );
   }

   return SCIP_OKAY;
}

/** catches events for variable at given position */
static
SCIP_RETCODE catchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_VAR* var;

   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nbinvars);
   assert(consdata->binvars != NULL);

   var = consdata->binvars[pos];
   assert(var != NULL);

   /* catch bound change events on variable */
   /**@todo do we have to add the event SCIP_EVENTTYPE_VARFIXED? */
   SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      consdata->nfixedzeros++;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      consdata->nfixedones++;

   return SCIP_OKAY;
}

/** drops events for variable at given position */
static
SCIP_RETCODE dropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_VAR* var;

   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nbinvars);
   assert(consdata->binvars != NULL);

   var = consdata->binvars[pos];
   assert(var != NULL);

   /* drop events on variable */
   SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      consdata->nfixedzeros--;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      consdata->nfixedones--;

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed linking constraint */
static
SCIP_RETCODE catchAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* author bzfhende
    *
    * TODO should we catch events even in the trivial case of only 1 binary variable
    */

   /* catch event for every single variable */
   for( i = 0; i < consdata->nbinvars; ++i )
   {
      SCIP_CALL( catchEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed linking constraint */
static
SCIP_RETCODE dropAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* author bzfhende
    *
    * TODO drop the events even in the trivial case nbinvars == 1?
    */

   /* drop event of every single variable */
   for( i = 0; i < consdata->nbinvars; ++i )
   {
      SCIP_CALL( dropEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** linearize the given linking constraint into a set partitioning constraint for the binary variables and a linear
 *  constraint for the linking between the linking variable and the binary variables */
static
SCIP_RETCODE consdataLinearize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_CONSDATA*        consdata            /**< linking constraint data */
   )
{
   SCIP_CONS* lincons;
   int b;

   SCIPdebugMsg(scip, "linearized linking constraint <%s>\n", SCIPconsGetName(cons));

   /* create set partitioning constraint for the binary variables */
   SCIP_CALL( SCIPcreateConsSetpart(scip, &lincons, SCIPconsGetName(cons), consdata->nbinvars, consdata->binvars,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

   /* create linear constraint for the linking between the binary variables and the linking variable */
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 0, NULL, NULL, 0.0, 0.0,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

   for( b = 0; b < consdata->nbinvars; ++b )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->binvars[b], consdata->vals[b]) );
   }
   SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->linkvar, -1.0) );

   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

   return SCIP_OKAY;
}

/** creates the binary  variables */
static
SCIP_RETCODE consdataCreateBinvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for bound change events on binary variables */
   SCIP_Bool             linearize           /**< should the linking constraint be linearized */
   )
{
   SCIP_VAR* linkvar;
   SCIP_VAR* binvar;
   int lb;
   int ub;
   char name[SCIP_MAXSTRLEN];
   int nbinvars;
   int b;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nbinvars == 0);
   assert(consdata->binvars == NULL);

   SCIPdebugMsg(scip, "create binary variables for linking variable <%s>\n", SCIPvarGetName(consdata->linkvar));

   /* author bzfhende
    *
    * TODO ensure that this method is only called for integer linking variables, because it does not make sense for continuous linking variables.
    */

   linkvar = consdata->linkvar;
   lb = SCIPconvertRealToInt(scip, SCIPvarGetLbGlobal(linkvar));
   ub = SCIPconvertRealToInt(scip, SCIPvarGetUbGlobal(linkvar));

   nbinvars = ub - lb + 1;
   assert(nbinvars > 0);

   /* allocate block memory for the binary variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->binvars, nbinvars) );
   /* allocate block memory for the binary variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vals, nbinvars) );
   consdata->sizebinvars = nbinvars;

   /* check if the linking variable is fixed */
   if( nbinvars == 1 )
   {
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s[%d]", SCIPvarGetName(linkvar), lb);

      /* creates and captures a fixed binary variables */
      SCIP_CALL( SCIPcreateVar(scip, &binvar, name, 1.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
            FALSE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, binvar) );

      consdata->binvars[0] = binvar;
      consdata->vals[0] = lb;
   }
   else
   {
      for( b = 0; b < nbinvars; ++b)
      {
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s[%d]", SCIPvarGetName(linkvar), lb + b);

         /* creates and captures variables */
         SCIP_CALL( SCIPcreateVar(scip, &binvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
               TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );

         /* add variable to the problem */
         SCIP_CALL( SCIPaddVar(scip, binvar) );
         consdata->binvars[b] = binvar;
         consdata->vals[b] = lb + b;
      }
   }

   consdata->nbinvars = nbinvars;

   assert(consdata->nfixedzeros == 0);
   assert(consdata->nfixedones == 0);

   if( SCIPisTransformed(scip) )
   {
      /* (rounding) lock binary variable */
      SCIP_CALL( lockRounding(scip, cons, consdata->binvars, consdata->nbinvars) );

      /* catch bound change events of variables */
      SCIP_CALL( catchAllEvents(scip, consdata, eventhdlr) );

      if( nbinvars > 1 )
      {
         if( linearize )
         {
            SCIP_CALL( consdataLinearize(scip, cons, consdata) );
         }
         else
         {
            /* enable constraint */
            SCIP_CALL( SCIPenableCons(scip, cons) );
         }
      }
   }

   return SCIP_OKAY;
}

/** creates consdata */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure  */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   SCIP_VAR*             linkvar,            /**< linking variable which is linked */
   SCIP_VAR**            binvars,            /**< binary variables */
   SCIP_Real*            vals,               /**< coefficients of the binary variables */
   int                   nbinvars            /**< number of binary starting variables */
   )
{
   int v;

   assert(scip!= NULL);
   assert(consdata != NULL);
   assert(linkvar != NULL);
   assert(binvars != NULL || nbinvars == 0);
   assert(SCIPvarGetType(linkvar) != SCIP_VARTYPE_CONTINUOUS || nbinvars > 0);

   /* allocate memory for consdata */
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->linkvar = linkvar;
   (*consdata)->nbinvars = nbinvars;
   (*consdata)->sizebinvars = nbinvars;
   (*consdata)->row1 = NULL;
   (*consdata)->row2 = NULL;
   (*consdata)->nlrow1 = NULL;
   (*consdata)->nlrow2 = NULL;
   (*consdata)->cliqueadded = FALSE;

   /* initialize constraint state */
   (*consdata)->sorted = FALSE;
   (*consdata)->firstnonfixed = 0;
   (*consdata)->lastnonfixed = nbinvars - 1;
   (*consdata)->nfixedzeros = 0;
   (*consdata)->nfixedones = 0;

   if( nbinvars == 0 )
   {
      (*consdata)->binvars = NULL;
      (*consdata)->vals = NULL;
   }
   else
   {
      /* copy binary variable array */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->binvars, binvars, nbinvars) );

      /* copy coefficients */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vals, vals, nbinvars) );
   }

   /* get transformed variable, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      if( nbinvars > 0 )
      {
         SCIP_CALL( SCIPgetTransformedVars(scip, nbinvars, (*consdata)->binvars, (*consdata)->binvars) );

         /* catch bound change events of variables */
         SCIP_CALL( catchAllEvents(scip, *consdata, eventhdlr) );
      }

      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->linkvar, &(*consdata)->linkvar) );
   }

   /* author bzfhende
    *
    * TODO do we need to forbid multi-aggregations? This was only needed if we substitute and resubstitute linking
    * variables into linear constraints.
    */

   /* capture variables */
   for( v = 0; v < nbinvars; ++v )
   {
      assert((*consdata)->binvars[v] != NULL);
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->binvars[v]) );
   }
   SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->linkvar) );

   return SCIP_OKAY;
}


/** free consdata */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure  */
   SCIP_CONSDATA**       consdata            /**< pointer to consdata */
   )
{
   int v;

   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->nbinvars == 0 || (*consdata)->binvars != NULL);

   /* release the rows */
   if( (*consdata)->row1 != NULL )
   {
      assert((*consdata)->row2 != NULL);

      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row1) );
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row2) );
   }

   /* release the nlrows */
   if( (*consdata)->nlrow1 != NULL )
   {
      assert((*consdata)->nlrow2 != NULL);

      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow1) );
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow2) );
   }

   /* capture variables */
   for( v = 0; v < (*consdata)->nbinvars; ++v )
   {
      assert((*consdata)->binvars[v] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->binvars[v]) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->linkvar) );

   /* free binary variable array */
   if( (*consdata)->sizebinvars > 0 )
   {
      /* if constraint belongs to transformed problem space, drop bound change events on variables */
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->vals, (*consdata)->sizebinvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->binvars, (*consdata)->sizebinvars);
   }

   /* check if the fixed counters are reset */
   assert((*consdata)->nfixedzeros == 0);
   assert((*consdata)->nfixedones == 0);

   /* free constraint data */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** analyzes conflicting assignment on given constraint where reason comes from the linking variable lower or upper
 *  bound
 */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   SCIP_VAR*             linkvar,            /**< linking variable  */
   SCIP_VAR*             binvar,             /**< binary variable is the reason */
   SCIP_Bool             lblinkvar,          /**< lower bound of linking variable is the reason */
   SCIP_Bool             ublinkvar           /**< upper bound of linking variable is the reason */
   )
{
   assert(scip != NULL);

   /* conflict analysis can only be applied in solving stage and if it is turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   if( lblinkvar )
   {
      assert(linkvar != NULL);
      SCIP_CALL( SCIPaddConflictLb(scip, linkvar, NULL) );
   }

   if( ublinkvar )
   {
      assert(linkvar != NULL);
      SCIP_CALL( SCIPaddConflictUb(scip, linkvar, NULL) );
   }

   if( binvar != NULL )
   {
      SCIP_CALL( SCIPaddConflictBinvar(scip, binvar) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/* author bzfhende
 *
 * TODO check if the method below produces valid results even if the variable is continuous
 */

/** fix linking variable to the value of the binary variable at pos */
static
SCIP_RETCODE consFixLinkvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   int                   pos,                /**< position of binary variable */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if the node can be cut off */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* linkvar;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Real coef;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   linkvar = consdata->linkvar;
   coef = consdata->vals[pos];

   /* change lower bound of the linking variable */
   SCIP_CALL( SCIPinferVarLbCons(scip, linkvar, coef, cons, pos, TRUE, &infeasible, &tightened) );

   if( infeasible )
   {
      assert(coef > SCIPvarGetUbLocal(linkvar));
      assert(coef >= SCIPvarGetLbLocal(linkvar));

      SCIP_CALL( analyzeConflict(scip, cons, linkvar, consdata->binvars[pos], FALSE, TRUE) );

      *cutoff = TRUE;
      return SCIP_OKAY;
   }
   assert(SCIPisFeasLE(scip, coef, SCIPvarGetUbLocal(linkvar)));

   /* change upper bound of the integer variable */
   SCIP_CALL( SCIPinferVarUbCons(scip, linkvar, coef, cons, pos, TRUE, &infeasible, &tightened) );

   if( infeasible )
   {
      assert(coef < SCIPvarGetLbLocal(linkvar));
      assert(coef <= SCIPvarGetUbLocal(linkvar));

      SCIP_CALL( analyzeConflict(scip, cons, linkvar, consdata->binvars[pos], TRUE, FALSE) );

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   assert(SCIPisFeasEQ(scip, SCIPvarGetUbLocal(linkvar), SCIPvarGetLbLocal(linkvar)));

   return SCIP_OKAY;
}

/** checks constraint for violation from the local bound of the linking variable, applies fixings to the binary
 *  variables if possible
 */
static
SCIP_RETCODE processRealBoundChg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nchgbds,            /**< pointer to store the number of changes (foxed) variable bounds */
   SCIP_Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_VAR* linkvar;
   SCIP_Real* vals;
   SCIP_Real lb;
   SCIP_Real ub;
   int nbinvars;
   int b;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* ensure that the binary variables are sorted in non-decreasing order w.r.t. their coefficients */
   consdataSort(consdata);

   nbinvars = consdata->nbinvars;

   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(nbinvars > 1);

   /* if more than one binary variable is fixed to one or at least nbinvars minus one variable are fixed to zero */
   if( consdata->nfixedones > 0  || consdata->nfixedzeros >= nbinvars-1 )
      return  SCIP_OKAY;

   linkvar = consdata->linkvar;
   assert(linkvar != NULL);

   binvars = consdata->binvars;
   vals = consdata->vals;

   lb = SCIPvarGetLbLocal(linkvar);
   ub = SCIPvarGetUbLocal(linkvar);

   assert(lb <= ub);

#ifndef NDEBUG
   /* check that the first variable are locally fixed to zero */
   for( b = 0; b < consdata->firstnonfixed; ++b )
      assert(SCIPvarGetUbLocal(binvars[b]) < 0.5);

   /* check that the last variable are locally fixed to zero */
   for( b = consdata->lastnonfixed + 1; b < nbinvars; ++b )
      assert(SCIPvarGetUbLocal(binvars[b]) < 0.5);
#endif

   for( b = consdata->firstnonfixed; b < nbinvars; ++b )
   {
      if( SCIPisLT(scip, vals[b], lb) )
      {
         SCIP_VAR* var;

         var =  binvars[b];
         assert(var != NULL);

         SCIPdebugMsg(scip, "fix variable <%s> to zero due to the lower bound of the linking variable <%s> [%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetName(linkvar), lb, ub);

         SCIP_CALL( SCIPinferBinvarCons(scip, var, FALSE, cons, -2, &infeasible, &tightened) );

         if( infeasible )
         {
            SCIP_CALL( analyzeConflict(scip, cons, linkvar, var, TRUE, FALSE) );
            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         if( tightened )
            (*nchgbds)++;

         /* adjust constraint state */
         consdata->firstnonfixed++;
      }
      else
         break;
   }

   /* fix binary variables to zero if not yet fixed, from local upper bound + 1*/
   for( b = consdata->lastnonfixed; b >= 0; --b )
   {
      if( SCIPisGT(scip, vals[b], ub) )
      {
         SCIP_VAR* var;

         var = binvars[b];
         assert(var != NULL);

         SCIPdebugMsg(scip, "fix variable <%s> to zero due to the upper bound of the linking variable <%s> [%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetName(linkvar), lb, ub);

         SCIP_CALL( SCIPinferBinvarCons(scip, var, FALSE, cons, -3, &infeasible, &tightened) );

         if( infeasible )
         {
            SCIP_CALL( analyzeConflict(scip, cons, linkvar, var, FALSE, TRUE) );
            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         if( tightened )
            (*nchgbds)++;

         /* adjust constraint state */
         consdata->lastnonfixed--;
      }
      else
         break;
   }

   if( consdata->firstnonfixed > consdata->lastnonfixed )
   {
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   *mustcheck = (*nchgbds) == 0;

   /* if linking variable is fixed, create for the binary variables which have a coefficient equal to the fixed value a
    * set partitioning constraint
    */
   if( SCIPisEQ(scip, lb, ub) )
   {
      if( consdata->firstnonfixed == consdata->lastnonfixed )
      {
         SCIP_VAR* var;

         var = binvars[consdata->firstnonfixed];

         SCIPdebugMsg(scip, "fix variable <%s> to one due to the fixed linking variable <%s> [%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetName(linkvar), lb, ub);

         /* TODO can the forbidden cases be covered more elegantly? */
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
            return SCIP_OKAY;

         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED )
            if( SCIPvarGetStatus(SCIPvarGetAggrVar(var)) == SCIP_VARSTATUS_MULTAGGR ||
                  SCIPvarGetStatus(SCIPvarGetAggrVar(var)) == SCIP_VARSTATUS_AGGREGATED )
               return SCIP_OKAY;

         SCIP_CALL( SCIPinferBinvarCons(scip, var, TRUE, cons, -6, &infeasible, &tightened) );

         if( infeasible )
         {
            SCIP_CALL( analyzeConflict(scip, cons, linkvar, var, TRUE, TRUE) );
            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         if( tightened )
            (*nchgbds)++;

         SCIPdebugMsg(scip, " -> disabling linking constraint <%s>\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );

         *mustcheck = FALSE;
      }
      else if( SCIPgetDepth(scip) <= 0 )
      {
         SCIP_CONS* setppc;
         SCIP_VAR** vars;
         int nvars;

         /* get sub array of variables which have the same coefficient */
         vars = &consdata->binvars[consdata->firstnonfixed];
         nvars = consdata->lastnonfixed - consdata->firstnonfixed + 1;

         SCIP_CALL( SCIPcreateConsSetpart(scip, &setppc, SCIPconsGetName(cons), nvars, vars,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

         SCIP_CALL( SCIPaddCons(scip, setppc) );
         SCIP_CALL( SCIPreleaseCons(scip, &setppc) );

         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from the binary variable array */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_CONS*            cons,               /**< linking constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nbinvars);

   var = consdata->binvars[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* remove the rounding locks for the deleted variable */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );

   /* if we are in transformed problem, delete the event data of the variable */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      SCIP_CALL( dropEvent(scip, consdata, conshdlrdata->eventhdlr, pos) );
   }

   /* move the last variable to the free slot */
   if( pos != consdata->nbinvars - 1 )
   {
      consdata->binvars[pos] = consdata->binvars[consdata->nbinvars-1];
      consdata->vals[pos] = consdata->vals[consdata->nbinvars-1];
      consdata->sorted = FALSE;
   }

   consdata->nbinvars--;

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

/** remove the trailing and leading binary variables that are fixed to zero */
static
SCIP_RETCODE removeFixedBinvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int nbinvars;
   int b;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->sorted);

   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetDepth(scip) <= 0);
   assert(!SCIPinProbing(scip));
   assert(!SCIPinRepropagation(scip));

   nbinvars = consdata->nbinvars;

   for( b = nbinvars - 1; b > consdata->lastnonfixed; --b )
   {
      SCIP_CALL( delCoefPos(scip, eventhdlr, cons, b) );
   }

   for( b = consdata->firstnonfixed - 1; b >= 0; --b )
   {
      SCIP_CALL( delCoefPos(scip, eventhdlr, cons, b) );
   }

   for( b = consdata->nbinvars - 1; b >= 0; --b )
   {
      if( SCIPvarGetUbLocal(consdata->binvars[b]) < 0.5 )
      {
         SCIP_CALL( delCoefPos(scip, eventhdlr, cons, b) );
      }
   }

   /* set the constraint state */
   consdata->firstnonfixed = 0;
   consdata->lastnonfixed = consdata->nbinvars - 1;

   return SCIP_OKAY;
}

/** tightened the linking variable due to binary variables which are fixed to zero */
static
SCIP_RETCODE tightenedLinkvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   SCIP_CONSDATA*        consdata,           /**< linking constraint to be processed */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nchgbds             /**< pointer to store the number of changed variable bounds */
   )
{
   SCIP_VAR** binvars;
   SCIP_VAR* linkvar;
   SCIP_Real* vals;

   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   int nbinvars;
   int b;

   /* if more than one binary variable is fixed to one or at least nbinvars minus one variable are fixed to zero return */
   if( consdata->nfixedones > 1 || consdata->nfixedzeros >= consdata->nbinvars-1 )
      return SCIP_OKAY;

   if( *cutoff )
      return SCIP_OKAY;

   assert(consdata->sorted);

   linkvar = consdata->linkvar;
   binvars = consdata->binvars;
   vals = consdata->vals;
   nbinvars = consdata->nbinvars;

#ifndef NDEBUG
   /* check that the first variable are locally fixed to zero */
   for( b = 0; b < consdata->firstnonfixed; ++b )
      assert(SCIPvarGetUbLocal(binvars[b]) < 0.5);
#endif

   assert(consdata->firstnonfixed < nbinvars);
   assert(consdata->lastnonfixed < nbinvars);

   /* find first non fixed binary variable */
   for( b = consdata->firstnonfixed; b < nbinvars; ++b )
   {
      if( SCIPvarGetUbLocal(binvars[b]) > 0.5 )
         break;

      consdata->firstnonfixed++;
   }

   SCIP_CALL( SCIPinferVarLbCons(scip, linkvar, vals[b], cons, -4, TRUE, &infeasible, &tightened) );

   /* start conflict analysis if infeasible */
   if( infeasible )
   {
      /* analyze the cutoff if if SOLVING stage and conflict analysis is turned on */
      if( (SCIPgetStage(scip) == SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) && SCIPisConflictAnalysisApplicable(scip) )
      {
         SCIPdebugMsg(scip, "conflict at <%s> due to bounds and fixed binvars: [lb,ub] = [%g,%g]; b= %d; coef = %g \n",
            SCIPvarGetName(linkvar), SCIPvarGetLbLocal(linkvar), SCIPvarGetUbLocal(linkvar), b, vals[b]);

         SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

         /* ??????????? use resolve method and only add binvars which are needed to exceed the upper bound */

         /* add conflicting variables */
         SCIP_CALL( SCIPaddConflictUb(scip, linkvar, NULL) );

         for( b = 0;  b < consdata->firstnonfixed; ++b )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, binvars[b]) );
         }

         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   if( tightened )
      (*nchgbds)++;

#ifndef NDEBUG
   /* check that the last variable are locally fixed to zero */
   for( b = consdata->lastnonfixed + 1; b < nbinvars; ++b )
      assert(SCIPvarGetUbLocal(binvars[b]) < 0.5);
#endif

   /* find last non fixed variable */
   for( b = consdata->lastnonfixed; b >= 0; --b )
   {
      if( SCIPvarGetUbLocal(binvars[b]) > 0.5 )
         break;

      consdata->lastnonfixed--;
   }

   if( SCIPvarGetStatus(SCIPvarGetProbvar(linkvar)) != SCIP_VARSTATUS_MULTAGGR )
      SCIP_CALL( SCIPinferVarUbCons(scip, linkvar, (SCIP_Real)vals[b], cons, -5, TRUE, &infeasible, &tightened) );

   if( infeasible )
   {
      /* conflict analysis can only be applied in solving stage and if conflict analysis is turned on */
      if( (SCIPgetStage(scip) == SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) && SCIPisConflictAnalysisApplicable(scip) )
      {
         SCIPdebugMsg(scip, "conflict at <%s> due to bounds and fixed binvars: [lb,ub] = [%g,%g]; b = %d; coef = %g,\n",
            SCIPvarGetName(linkvar), SCIPvarGetLbLocal(linkvar), SCIPvarGetUbLocal(linkvar), b, vals[b]);

         SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

         /* ??????????? use resolve method and only add binvars which are needed to fall below the lower bound */

         /* add conflicting variables */
         SCIP_CALL( SCIPaddConflictLb(scip, linkvar, NULL) );

         for( b = consdata->lastnonfixed + 1; b < nbinvars; ++b )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, binvars[b]) );
         }

         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   if( tightened )
      (*nchgbds)++;

   return SCIP_OKAY;
}

/** checks constraint for violation only looking at the fixed binary variables, applies further fixings if possible */
static
SCIP_RETCODE processBinvarFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nchgbds,            /**< pointer to store the number of changed variable bounds */
   SCIP_Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   SCIP_Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(addcut != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nbinvars == 0 || consdata->binvars != NULL);
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nbinvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nbinvars);

   /* ensure that the binary variables are sorted in non-decreasing order w.r.t. their coefficients */
   consdataSort(consdata);

   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(consdata->nbinvars > 1);

   if( *cutoff )
      return SCIP_OKAY;

   if( consdata->nfixedones == 1 )
   {
      /* exactly one variable is fixed to 1:
       * - all other binary variables in a set partitioning must be zero
       * - linking variable is fixed to that binary variable
       */
      if( consdata->nfixedzeros < consdata->nbinvars - 1 ||
         SCIPisLT(scip, SCIPvarGetLbLocal(consdata->linkvar), SCIPvarGetUbLocal(consdata->linkvar)) )
      {
         SCIP_VAR** vars;
         SCIP_VAR* var;
#ifndef NDEBUG
         SCIP_Bool fixedonefound;
#endif
         int nvars;
         int v;

         SCIPdebugMsg(scip, " -> fixing all other variables to zero due to the set partitioning condition <%s>\n",
            SCIPconsGetName(cons));

         /* unfixed variables exist: fix them to zero;
          * this could result in additional variables fixed to one due to aggregations; in this case, the
          * constraint is infeasible in local bounds
          */
         vars = consdata->binvars;
         nvars = consdata->nbinvars;
#ifndef NDEBUG
         fixedonefound = FALSE;
#endif

         for( v = 0; v < nvars && consdata->nfixedones == 1 && !(*cutoff); ++v )
         {
            var = vars[v];
            assert(SCIPvarIsBinary(var));
            /* TODO can this be handled more elegantly? */
            if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
               continue;

            if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED )
               if( SCIPvarGetStatus(SCIPvarGetAggrVar(var)) == SCIP_VARSTATUS_MULTAGGR ||
                     SCIPvarGetStatus(SCIPvarGetAggrVar(var)) == SCIP_VARSTATUS_AGGREGATED )
                  continue;

            if( SCIPvarGetLbLocal(var) < 0.5 )
            {
               SCIP_CALL( SCIPinferBinvarCons(scip, var, FALSE, cons, -1, &infeasible, &tightened) );
               assert(!infeasible);
               SCIPdebugMsg(scip, "   -> fixed <%s> to zero (tightened=%u)\n", SCIPvarGetName(var), tightened);
            }
            else
            {
#ifndef NDEBUG
               fixedonefound = TRUE;
#endif
               /* fix linking variable */
               /* TODO check if variable status allows fixing (probably in consFixLinkvar) */
               SCIP_CALL( consFixLinkvar(scip, cons, v, cutoff) );
            }
         }
         if( !(*cutoff) )
         {
            /* the fixed to one variable must have been found, and at least one variable must have been fixed */
            assert(consdata->nfixedones >= 1 || fixedonefound);

            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            (*nchgbds)++;
         }
      }

      /* now all other variables are fixed to zero:
       * the constraint is feasible, and if it's not modifiable, it is redundant
       */
      if( !SCIPconsIsModifiable(cons) && consdata->nfixedones == 1 )
      {
         SCIPdebugMsg(scip, " -> disabling set linking constraint <%s>\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }
   }
   else if( consdata->nfixedones >= 2 )
   {
      /* at least two variables are fixed to 1:
       * - the set partitioning condition is violated
       */
      SCIPdebugMsg(scip, " -> conflict on " CONSHDLR_NAME " constraint <%s> due to the set partitioning condition\n", SCIPconsGetName(cons));

      SCIP_CALL( SCIPresetConsAge(scip, cons) );

      /* conflict analysis can only be applied in solving stage and if it is applicable */
      if( (SCIPgetStage(scip) == SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) && SCIPisConflictAnalysisApplicable(scip) )
      {
         SCIP_VAR** vars;
         int nvars;
         int n;
         int v;

         vars = consdata->binvars;
         nvars = consdata->nbinvars;

         /* initialize conflict analysis, and add the two variables assigned to one to conflict candidate queue */
         SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

         n = 0;

         for( v = 0; v < nvars && n < 2; ++v )
         {
            if( SCIPvarGetLbLocal(vars[v]) > 0.5 )
            {
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[v]) );
               n++;
            }
         }
         assert(n == 2);

         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      *cutoff = TRUE;
   }
   else if( consdata->nfixedzeros == consdata->nbinvars )
   {
      /* all variables are fixed to zero:
       * - the set partitioning condition is violated, and if it's unmodifiable, the node
       *   can be cut off -- otherwise, the constraint must be added as a cut and further pricing must
       *   be performed
       */
      assert(consdata->nfixedones == 0);

      SCIPdebugMsg(scip, " -> " CONSHDLR_NAME " constraint <%s> is infeasible due to the set partitioning condition\n",
         SCIPconsGetName(cons));

      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      if( SCIPconsIsModifiable(cons) )
         *addcut = TRUE;
      else
      {
         /* conflict analysis can only be applied in solving stage and if it is applicable */
         if( (SCIPgetStage(scip) == SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) && SCIPisConflictAnalysisApplicable(scip) )
         {
            SCIP_VAR** vars;
            int nvars;
            int v;

            vars = consdata->binvars;
            nvars = consdata->nbinvars;

            /* initialize conflict analysis, add all variables of infeasible constraint to conflict candidate queue */
            SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

            for( v = 0; v < nvars; ++v )
            {
               assert(SCIPvarGetUbLocal(vars[v]) < 0.5);
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[v]) );
            }

            /* analyze the conflict */
            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }
         *cutoff = TRUE;
      }
   }
   else if( consdata->nfixedzeros == consdata->nbinvars - 1 )
   {
      /* all variables except one are fixed to zero:
       * - an unmodifiable set partitioning constraint is feasible and can be disabled after the
       *   remaining variable is fixed to one
       * - a modifiable set partitioning constraint must be checked manually
       */
      assert(consdata->nfixedones == 0);

      if( !SCIPconsIsModifiable(cons) )
      {
         SCIP_VAR** vars;
         SCIP_VAR* var;
         int nvars;
         int v;

         /* search the single variable that can be fixed */
         vars = consdata->binvars;
         nvars = consdata->nbinvars;
         for( v = 0; v < nvars && !(*cutoff); ++v )
         {
            var = vars[v];
            assert(SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)));
            assert(SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) || SCIPisFeasEQ(scip, SCIPvarGetUbLocal(var), 1.0));

            if( SCIPvarGetUbLocal(var) > 0.5 )
            {
               assert(SCIPvarGetLbLocal(var) < 0.5);
               SCIPdebugMsg(scip, " -> fixing remaining binary variable <%s> to one in " CONSHDLR_NAME " constraint <%s>\n",
                  SCIPvarGetName(var), SCIPconsGetName(cons));

               if( SCIPvarGetStatus(SCIPvarGetProbvar(var)) != SCIP_VARSTATUS_MULTAGGR )
               {
                  SCIP_CALL( SCIPinferBinvarCons(scip, var, TRUE, cons, -1, &infeasible, &tightened) );
                  assert(!infeasible);
                  assert(tightened);
               }

               /* fix linking variable */
               /* TODO check if variable status allows fixing (probably in consFixLinkvar)*/
               SCIP_CALL( consFixLinkvar(scip, cons, v, cutoff) );
               break;
            }
         }
         assert(v < nvars);
         assert(consdata->nfixedzeros == consdata->nbinvars - 1);
         assert(consdata->nfixedones == 1);

         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         (*nchgbds)++;
      }
   }
   else
   {
      SCIP_CALL( tightenedLinkvar(scip, cons, consdata, cutoff, nchgbds) );
   }

   *mustcheck = (*nchgbds) == 0;

   assert(consdata->nfixedzeros + consdata->nfixedones <= consdata->nbinvars);

   return SCIP_OKAY;
}

/** returns whether the given solution is feasible for the given linking constraint */
static
SCIP_Bool checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be checked */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo solution */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_Real* vals;
   SCIP_Real solval;
   SCIP_Real linksum;
   SCIP_Real linkvarval;
   SCIP_Real setpartsum;
   SCIP_Real setpartsumbound;
   SCIP_Real absviol;
   SCIP_Real relviol;
   int nbinvars;
   int b;

   assert(scip != NULL);
   assert(cons != NULL);

   SCIPdebugMsg(scip, "checking linking constraint <%s> for feasibility of solution %p\n", SCIPconsGetName(cons), (void*)sol);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->binvars != NULL || consdata->nbinvars == 0);

   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(consdata->nbinvars > 1);

   /* calculate the constraint's activity for the linking part and the set partitioning part */
   binvars = consdata->binvars;
   vals = consdata->vals;
   nbinvars = consdata->nbinvars;

   linksum = 0.0;
   setpartsum = 0.0;
   setpartsumbound = 1.0 + 2*SCIPfeastol(scip);

   for( b = 0; b < nbinvars && setpartsum < setpartsumbound; ++b )  /* if sum >= sumbound, the feasibility is clearly decided */
   {
      assert(SCIPvarIsBinary(binvars[b]));

      solval = SCIPgetSolVal(scip, sol, binvars[b]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));

      linksum += vals[b] * solval;
      setpartsum += solval;
   }

   /* calculate and update absolute and relative violation of the equality constraint */
   linkvarval = SCIPgetSolVal(scip, sol, consdata->linkvar);
   absviol = REALABS(linksum - linkvarval);
   relviol = REALABS(SCIPrelDiff(linksum, linkvarval));
   if( sol != NULL )
      SCIPupdateSolLPConsViolation(scip, sol, absviol, relviol);

   /* calculate and update absolute and relative violation of the set partitioning constraint */
   absviol = REALABS(setpartsum - 1.0);
   relviol = REALABS(SCIPrelDiff(setpartsum, 1.0));
   if( sol != NULL )
      SCIPupdateSolLPConsViolation(scip, sol, absviol, relviol);

   /* check if the fixed binary variable match with the linking variable */
   return SCIPisFeasEQ(scip, linksum, linkvarval) && SCIPisFeasEQ(scip, setpartsum, 1.0);
}

#if 0
/** transfer aggregated integer variables to the corresponding binary variables */
static
SCIP_RETCODE aggregateVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< hash map mapping a integer variables to its linking constraint */
   SCIP_CONS**           conss,              /**< array of linking constraint */
   int                   nconss,             /**< number of linking constraints */
   int*                  naggrvars,          /**< pointer to store the number of aggregate variables */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_CONS* aggrcons;
   SCIP_CONSDATA* aggrconsdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_VAR** aggrbinvars;
   SCIP_VAR* linkvar;
   SCIP_VAR* aggrvar;
   SCIP_Real aggrconst;
   SCIP_Real aggrscalar;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;
   int offset;
   int aggroffset;
   int nbinvars;
   int shift;
   int b;
   int c;

   assert(varmap != NULL);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      linkvar = consdata->linkvar;
      assert(linkvar != NULL);

      if( SCIPvarGetStatus(linkvar) == SCIP_VARSTATUS_AGGREGATED )
      {
         aggrvar =  SCIPvarGetAggrVar(linkvar);
         aggrcons = (SCIP_CONS*) SCIPhashmapGetImage(varmap, getHashmapKey(aggrvar));

         /* check if the aggregate variable belongs to a linking constraint */
         if( aggrcons != NULL )
         {
            aggrconsdata = SCIPconsGetData(aggrcons);
            assert(aggrconsdata != NULL);

            aggrconst = SCIPvarGetAggrConstant(linkvar);
            aggrscalar = SCIPvarGetAggrScalar(linkvar);

            /**@todo extend the aggregation for those cases were the aggrscalar is not equal to 1.0 */
            if( SCIPisEQ(scip, aggrscalar, 1.0 ) )
            {
               /* since both variables are integer variable and the aggrscalar is 1.0 the aggrconst should
                * integral
                */
               assert(SCIPisIntegral(scip, aggrconst));
               shift = SCIPconvertRealToInt(scip, aggrconst);

               offset = consdata->offset;
               binvars = consdata->binvars;
               aggroffset = aggrconsdata->offset;
               aggrbinvars = aggrconsdata->binvars;

               nbinvars = MIN(consdata->nbinvars + offset, aggrconsdata->nbinvars + shift + aggroffset);

               for( b = MAX(offset, aggroffset-shift); b < nbinvars; ++b )
               {
                  assert(b - offset >= 0);
                  assert(b + shift - aggroffset >= 0);
                  assert(b < consdata->nbinvars);
                  assert(b < aggrconsdata->nbinvars - shift);

                  /* add aggregation x - y  = 0.0 */
                  SCIP_CALL( SCIPaggregateVars(scip, binvars[b-offset], aggrbinvars[b+shift-aggroffset], 1.0, -1.0, 0.0,
                        &infeasible, &redundant, &aggregated) );

                  if( infeasible )
                  {
                     (*cutoff) = TRUE;
                     return SCIP_OKAY;
                  }

                  if( aggregated )
                     (*naggrvars)++;
               }
            }
         }
      }
   }
   return SCIP_OKAY;
}
#endif

/** create two rows for the linking constraint
 *
 *  - row1: {sum_{b=1}^n-1 vals[b] * binvars[b]} - linkvar = 0
 *  - row2: {sum_{b=0}^n-1 binvars[b]} = 1.0
 */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;
   char rowname[SCIP_MAXSTRLEN];
   int b;

   assert( cons != NULL);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row1 == NULL);
   assert(consdata->row2 == NULL);
   assert(consdata->nbinvars > 1);

   /* create the LP row which captures the linking between the real and binary variables */
   (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s[link]", SCIPconsGetName(cons));

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->row1, cons, rowname, 0.0, 0.0,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   /* add linking variable to the row */
   assert(consdata->linkvar != NULL);
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->row1, consdata->linkvar, -1.0) );

   /* adding binary variables to the row */
   assert(consdata->binvars != NULL);
   for( b = 0; b < consdata->nbinvars; ++b )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row1, consdata->binvars[b], consdata->vals[b]) );
   }

   /* create the LP row which captures the set partitioning condition of the binary variables */
   (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s[setppc]", SCIPconsGetName(cons));
   assert( consdata->nbinvars > 0 );

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->row2, cons, rowname, 1.0, 1.0,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->row2, consdata->nbinvars, consdata->binvars, 1.0) );

   return SCIP_OKAY;
}


/** adds linking constraint as cut to the LP */
static
SCIP_RETCODE addCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_Bool*            cutoff              /**< whether a cutoff has been detected */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cutoff != NULL );
   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(consdata->nbinvars > 1);

   if( consdata->row1 == NULL )
   {
      assert(consdata->row2 == NULL);

      /* convert linking data into LP rows */
      SCIP_CALL( createRows(scip, cons) );
   }
   assert(consdata->row1 != NULL);
   assert(consdata->row2 != NULL);

   /* insert LP linking row as cut */
   if( !SCIProwIsInLP(consdata->row1) )
   {
      SCIPdebugMsg(scip, "adding linking row of constraint <%s> as cut to the LP\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPaddRow(scip, consdata->row1, TRUE/*FALSE*/, cutoff) );
   }

   /* insert LP set partitioning row as cut */
   if( !SCIProwIsInLP(consdata->row2) )
   {
      SCIPdebugMsg(scip, "adding set partitioning row of constraint <%s> as cut to the LP\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPaddRow(scip, consdata->row2, TRUE/*FALSE*/, cutoff) );
   }

   return SCIP_OKAY;
}

/** adds linking constraint as rows to the NLP, if not added yet */
static
SCIP_RETCODE addNlrow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(SCIPisNLPConstructed(scip));

   /* skip deactivated, redundant, or local constraints (the NLP does not allow for local rows at the moment) */
   if( !SCIPconsIsActive(cons) || !SCIPconsIsChecked(cons) || SCIPconsIsLocal(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow1 == NULL )
   {
      char rowname[SCIP_MAXSTRLEN];
      SCIP_Real* coefs;
      int i;

      assert(consdata->nlrow2 == NULL);

      /* create the NLP row which captures the linking between the real and binary variables */
      (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s[link]", SCIPconsGetName(cons));

      /* create nlrow1 with binary variables */
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow1, rowname,
         0.0, consdata->nbinvars, consdata->binvars, consdata->vals, NULL, 0.0, 0.0, SCIP_EXPRCURV_LINEAR) );
      /* add linking variable to the row */
      SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow1, consdata->linkvar, -1.0) );

      /* create the NLP row which captures the set partitioning condition of the binary variables */
      (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s[setppc]", SCIPconsGetName(cons));

      SCIP_CALL( SCIPallocBufferArray(scip, &coefs, consdata->nbinvars) );
      for( i = 0; i < consdata->nbinvars; ++i )
         coefs[i] = 1.0;

      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow2, rowname,
         0.0, consdata->nbinvars, consdata->binvars, coefs, NULL, 1.0, 1.0, SCIP_EXPRCURV_LINEAR) );

      SCIPfreeBufferArray(scip, &coefs);
   }

   assert(SCIPnlrowIsInNLP(consdata->nlrow1) == SCIPnlrowIsInNLP(consdata->nlrow2));
   if( !SCIPnlrowIsInNLP(consdata->nlrow1) )
   {
      SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow1) );
      SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow2) );
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cuts if possible */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            separated,          /**< pointer to store TRUE, if a cut was found */
   int*                  nchgbds             /**< pointer to store the number of changed variables bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(separated != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(consdata->nbinvars > 1);

   SCIPdebugMsg(scip, "separating constraint <%s>\n", SCIPconsGetName(cons));

   *cutoff = FALSE;
   addcut = FALSE;
   mustcheck = TRUE;

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   if( sol == NULL )
   {
      SCIP_CALL( processRealBoundChg(scip, cons, cutoff, nchgbds, &mustcheck) );
   }

   if( mustcheck && !(*cutoff) )
   {
      /* variable's fixings didn't give us any information -> we have to check the constraint */
      if( sol == NULL && consdata->row1 != NULL )
      {
         SCIP_Real feasibility;
         SCIP_Real tmp;

         assert(consdata->row2 != NULL);

         /* skip constraints already in the LP */
         if( SCIProwIsInLP(consdata->row1) && SCIProwIsInLP(consdata->row2))
            return SCIP_OKAY;

         feasibility = 1.0;

         /* check first row (linking) for feasibility */
         if( !SCIProwIsInLP(consdata->row1) )
         {
            tmp = SCIPgetRowLPFeasibility(scip, consdata->row1);
            feasibility = MIN(feasibility, tmp);
         }

         /* check second row (setppc) for feasibility */
         if( !SCIProwIsInLP(consdata->row2) )
         {
            tmp = SCIPgetRowLPFeasibility(scip, consdata->row2);
            feasibility = MIN(feasibility, tmp);
         }
         addcut = SCIPisFeasNegative(scip, feasibility);
      }
      else
         addcut = !checkCons(scip, cons, sol);

      if( !addcut )
      {
         /* constraint was feasible -> increase age */
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
   }

   if( addcut )
   {
      /* insert LP row as cut */
      assert(!(*cutoff));
      SCIP_CALL( addCuts(scip, cons, cutoff) );
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }

   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
SCIP_RETCODE enforcePseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be separated */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   int*                  nchgbds,            /**< pointer to store the number of changed variable bounds */
   SCIP_Bool*            solvelp             /**< pointer to store TRUE, if the LP has to be solved */
   )
{
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;

   assert(!SCIPhasCurrentNodeLP(scip));
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(nchgbds != NULL);
   assert(solvelp != NULL);

   addcut = FALSE;
   mustcheck = TRUE;

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   SCIP_CALL( processRealBoundChg(scip, cons, cutoff, nchgbds, &mustcheck) );
   SCIP_CALL( processBinvarFixings(scip, cons, cutoff, nchgbds, &addcut, &mustcheck) );

   if( mustcheck )
   {
      assert(!addcut);

      if( checkCons(scip, cons, NULL) )
      {
         /* constraint was feasible -> increase age */
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
      else
      {
         /* constraint was infeasible -> reset age */
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *infeasible = TRUE;
      }
   }

   if( addcut )
   {
      assert(!(*cutoff));
      /* a cut must be added to the LP -> we have to solve the LP immediately */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *solvelp = TRUE;
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
   int                   nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Enforcing %d linking constraints for %s solution\n", nconss, sol == NULL ? "LP" : "relaxation");

   cutoff = FALSE;
   separated = FALSE;
   nchgbds = 0;

   /* check all useful linking constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && nchgbds == 0; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &cutoff, &separated, &nchgbds) );
   }

   /* check all obsolete linking constraints for feasibility */
   for( c = nusefulconss; c < nconss && !cutoff && !separated && nchgbds == 0; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &cutoff, &separated, &nchgbds) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyLinking)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrLinking(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeLinking)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreLinking)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* disable all linking constraints which contain at most one binary variable */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* skip constraints which are not added */
      if( !SCIPconsIsAdded(conss[c]) )
         continue;

      if( consdata->nbinvars <= 1 )
      {
         SCIP_CALL( SCIPdisableCons(scip, conss[c]) );
         assert(consdata->nbinvars == 0 || SCIPvarGetLbGlobal(consdata->binvars[0]) > 0.5);
      }
      else if( conshdlrdata->linearize )
      {
         SCIP_CALL( consdataLinearize(scip, conss[c], consdata) );
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler */
static
SCIP_DECL_CONSINITSOL(consInitsolLinking)
{  /*lint --e{715}*/

   /* add nlrow representations to NLP, if NLP had been constructed */
   if( SCIPisNLPConstructed(scip) )
   {
      int c;
      for( c = 0; c < nconss; ++c )
      {
         SCIP_CALL( addNlrow(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* release the rows and nlrows of all constraints */
      if( consdata->row1 != NULL )
      {
         assert(consdata->row2 != NULL);

         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row1) );
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row2) );
      }

      if( consdata->nlrow1 != NULL )
      {
         assert(consdata->nlrow2 != NULL);

         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow1) );
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow2) );
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteLinking)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* remove linking constraint form variable hash map */
   assert(conshdlrdata->varmap != NULL);
   assert(SCIPhashmapExists(conshdlrdata->varmap, getHashmapKey((*consdata)->linkvar)));
   SCIP_CALL( SCIPhashmapRemove(conshdlrdata->varmap, getHashmapKey((*consdata)->linkvar)) );

   if( (*consdata)->nbinvars > 0 && SCIPisTransformed(scip) )
   {
      SCIP_CALL( dropAllEvents(scip, *consdata, conshdlrdata->eventhdlr) );
   }

   /* free consdata  */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row1 == NULL);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->row2 == NULL);  /* in original problem, there cannot be LP rows */

   SCIPdebugMsg(scip, "transform linking constraint for variable <%s>\n", SCIPvarGetName(sourcedata->linkvar));

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, conshdlrdata->eventhdlr, &targetdata,
         sourcedata->linkvar, sourcedata->binvars, sourcedata->vals, sourcedata->nbinvars) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* insert (transformed) linking constraint into the hash map */
   assert(conshdlrdata->varmap != NULL);
   SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varmap, getHashmapKey(targetdata->linkvar), *targetcons) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   *infeasible = FALSE;

   for( c = 0; c < nconss && !(*infeasible); ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->nbinvars <= 1 )
         continue;

      SCIP_CALL( addCuts(scip, conss[c], infeasible) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpLinking)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "separating %d/%d linking constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   separated = FALSE;
   nchgbds = 0;

   /* check all useful linking constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &separated, &nchgbds) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolLinking)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "separating %d/%d " CONSHDLR_NAME " constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   separated = FALSE;
   nchgbds = 0;

   /* check all useful set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &cutoff, &separated, &nchgbds) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpLinking)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxLinking)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsLinking)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   int nchgbds;
   SCIP_Bool solvelp;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "pseudo enforcing %d " CONSHDLR_NAME " constraints\n", nconss);

   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   cutoff = FALSE;
   infeasible = FALSE;
   nchgbds = 0;
   solvelp = FALSE;

   /* check all linking constraint for domain reductions and feasibility */
   for( c = 0; c < nconss && !cutoff && !solvelp; ++c )
   {
      SCIP_CALL( enforcePseudo(scip, conss[c], &cutoff, &infeasible, &nchgbds, &solvelp) );
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( solvelp )
      *result = SCIP_SOLVELP;
   else if( infeasible )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckLinking)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all linking constraints for feasibility */
   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE || completely); ++c )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( consdata->nbinvars > 1 && (checklprows || consdata->row1 == NULL || !SCIProwIsInLP(consdata->row1)) )
      {
         if( !checkCons(scip, cons, sol) )
         {
            /* constraint is violated */
            *result = SCIP_INFEASIBLE;

            if( printreason )
            {
               int pos;
               int b;

               pos = -1;

#ifndef NDEBUG
               for( b = 0; b < consdata->nbinvars; ++b )
               {
                  assert(consdata->binvars[b] != NULL);
                  assert(SCIPvarIsBinary(consdata->binvars[b]));
               }
#endif

               SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
               SCIPinfoMessage(scip, NULL, ";\n");

               /* check that at most one binary variable is fixed */
               for( b = 0; b < consdata->nbinvars; ++b )
               {
                  assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, consdata->binvars[b])) );

                  /* check if binary variable is fixed */
                  if( SCIPgetSolVal(scip, sol, consdata->binvars[b]) > 0.5 )
                  {
                     if( pos != -1 )
                     {
                        SCIPinfoMessage(scip, NULL, "violation: more than one binary variable is set to one");
                        break;
                     }
                     pos = b ;
                  }
               }

               /* check that at least one binary variable is fixed */
               if( pos == -1 )
               {
                  SCIPinfoMessage(scip, NULL, "violation: none of the binary variables is set to one\n");
               }
               else if( !SCIPisFeasEQ(scip, consdata->vals[pos], SCIPgetSolVal(scip, sol, consdata->linkvar)) )
               {
                  /* check if the fixed binary variable match with the linking variable */
                  SCIPinfoMessage(scip, NULL, "violation: <%s> = <%g> and <%s> is one\n",
                     SCIPvarGetName(consdata->linkvar), SCIPgetSolVal(scip, sol, consdata->linkvar),
                     SCIPvarGetName(consdata->binvars[pos]) );
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropLinking)
{  /*lint --e{715}*/
   SCIP_Bool cutoff = FALSE;
   int nchgbds = 0;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "propagating %d/%d " CONSHDLR_NAME " constraints\n", nusefulconss, nconss);

   /* propagate all useful set partitioning / packing / covering constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_Bool addcut;
      SCIP_Bool mustcheck;

      SCIP_CALL( processRealBoundChg(scip, conss[c], &cutoff, &nchgbds, &mustcheck) );
      SCIP_CALL( processBinvarFixings(scip, conss[c], &cutoff, &nchgbds, &addcut, &mustcheck) );
   } /*lint !e438*/

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolLinking)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int oldnfixedvars;
   int oldnchgbds;
   int oldnaggrvars;
   int oldndelconss;
   int firstchange;
   int firstclique;
   int lastclique;
   int c;
   SCIP_Bool fixed;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_Bool mustcheck;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "presolve %d linking constraints\n", nconss);

   (*result) = SCIP_DIDNOTFIND;

   oldnchgbds = *nchgbds;
   oldnaggrvars = *naggrvars;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;
   cutoff = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   firstchange = INT_MAX;
   firstclique = INT_MAX;
   lastclique = -1;

   /* check for each linking constraint the set partitioning condition */
   for( c = 0; c < nconss && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;

      assert(*result != SCIP_CUTOFF);

      cons = conss[c];
      assert(cons != NULL);
      assert(!SCIPconsIsModifiable(cons));

      SCIPdebugMsg(scip, "presolve linking constraints <%s>\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( !SCIPconsIsEnabled(cons) /* || consdata->nbinvars <= 1 */ )
         continue;

      /* in case there is only at most one binary variables, the constraints should already be disabled */
      assert(consdata->nbinvars > 1);

      /*SCIPdebugMsg(scip, "presolving set partitioning / packing / covering constraint <%s>\n", SCIPconsGetName(cons));*/
      if( consdata->nfixedones >= 2 )
      {
         /* at least two variables are fixed to 1:
          * - a linking constraint is infeasible due to the set partitioning condition
          */
         SCIPdebugMsg(scip, "" CONSHDLR_NAME " constraint <%s> is infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if( consdata->nfixedones == 1 )
      {
         /* exactly one variable is fixed to 1:
          * - all other binary variables must be zero due to the set partitioning condition
          * - linking variable has to be fixed to corresponding binary variable which is fixed to one
          * - if constraint is not modifiable it can be removed
          */
         SCIP_VAR* var;
         int v;

         SCIPdebugMsg(scip, "" CONSHDLR_NAME " constraint <%s> has a binary variable fixed to 1.0\n", SCIPconsGetName(cons));

         for( v = 0; v < consdata->nbinvars; ++v )
         {
            var = consdata->binvars[v];
            assert(var != NULL);

            if( SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5 )
            {
               SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );

               if( infeasible )
               {
                  SCIPdebugMsg(scip, "" CONSHDLR_NAME " constraint <%s>: infeasible fixing <%s> == 0\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var));

                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               assert(fixed);
               (*nfixedvars)++;
            }
            else if( SCIPvarGetLbGlobal(var) > 0.5 )
            {
               /* fix linking variable */
               assert(SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_LOOSE
                   || SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_AGGREGATED
                   || SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_COLUMN
                   || SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_FIXED
                   || SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_NEGATED);
               SCIP_CALL( SCIPfixVar(scip, consdata->linkvar, consdata->vals[v], &infeasible, &fixed) );

               if( infeasible )
               {
                  SCIPdebugMsg(scip, "" CONSHDLR_NAME " constraint <%s>: infeasible fixing <%s> == %g\n",
                     SCIPconsGetName(cons), SCIPvarGetName(consdata->linkvar), consdata->vals[v]);

                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }

               if( fixed )
                  (*nfixedvars)++;
            }
         }

         /* now all other variables are fixed to zero:
          * the constraint is feasible, and if it's not modifiable, it is redundant
          */
         SCIPdebugMsg(scip, "" CONSHDLR_NAME " constraint <%s> is redundant\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
         continue;
      }

      if( consdata->nfixedzeros == consdata->nbinvars )
      {
         /* all variables are fixed to zero:
          * - a linking constraint is infeasible due the set partitioning condition
          */
         assert(consdata->nfixedones == 0);

         SCIPdebugMsg(scip, "linking constraint <%s> is infeasible due to set partitioning condition\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if( consdata->nfixedzeros == consdata->nbinvars - 1 )
      {
         /* all variables except one are fixed to zero:
          * - a linking constraint is feasible due the set partitioning condition
          * - the remaining binary variable can be fixed to one
          * - linking variable has to be fixed to corresponding binary variable which  is fixed  to one
          * - constraint can be deleted since it is not modifiable
          */
         SCIP_VAR* var;
         int v;

         assert(consdata->nfixedones == 0);

         SCIPdebugMsg(scip, "" CONSHDLR_NAME " constraint <%s> has only one binary variable not fixed to zero\n",
            SCIPconsGetName(cons));

         /* search unfixed variable */
         /* intentional empty for loop to increment counter to proper position */
         /* TODO speed up loop by considering only variables between firstnonfixed and lastnonfixed */
         for( v = 0; v < consdata->nbinvars && SCIPvarGetUbGlobal(consdata->binvars[v]) < 0.5; ++v ); /*lint !e722*/
         assert(v < consdata->nbinvars);
         var = consdata->binvars[v];

         /* fix remaining binary variable */
         SCIP_CALL( SCIPfixVar(scip, var, 1.0, &infeasible, &fixed) );
         if( infeasible )
         {
            SCIPdebugMsg(scip, "" CONSHDLR_NAME " constraint <%s>: infeasible fixing <%s> == 1\n",
               SCIPconsGetName(cons), SCIPvarGetName(var));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         assert(fixed);
         (*nfixedvars)++;

         /* fix linking variable */
         assert(SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_LOOSE
             || SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_AGGREGATED
             || SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_COLUMN
             || SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_FIXED
             || SCIPvarGetStatus(consdata->linkvar) == SCIP_VARSTATUS_NEGATED);
         SCIP_CALL( SCIPfixVar(scip, consdata->linkvar, consdata->vals[v], &infeasible, &fixed) );

         if( infeasible )
         {
            SCIPdebugMsg(scip, CONSHDLR_NAME " constraint <%s>: infeasible fixing <%s> == %g\n",
               SCIPconsGetName(cons), SCIPvarGetName(consdata->linkvar), consdata->vals[v]);

            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         assert(!SCIPvarIsActive(consdata->linkvar) || fixed);
         if( fixed )
            (*nfixedvars)++;

         /* delete constraint from  problem */
         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
         continue;
      }

      if( consdata->nfixedzeros == consdata->nbinvars - 2 ) /*lint !e641*/
      {
         SCIP_VAR* var;
         SCIP_VAR* var1;
         SCIP_VAR* var2;
         SCIP_Bool redundant;
         SCIP_Bool aggregated;
         int v;

         /* aggregate variable, if set partitioning condition consists only of two
          * non-fixed variables
          */

         /* search unfixed variable */
         var1 = NULL;
         var2 = NULL;
         for( v = 0; v < consdata->nbinvars && var2 == NULL; ++v )
         {
            var = consdata->binvars[v];
            if( SCIPvarGetUbGlobal(var) > 0.5 )
            {
               if( var1 == NULL )
                  var1 = var;
               else
                  var2 = var;
            }
         }
         assert(var1 != NULL && var2 != NULL);

         /* aggregate binary equality var1 + var2 == 1 */
         SCIPdebugMsg(scip, "" CONSHDLR_NAME " constraint <%s>: aggregate <%s> + <%s> == 1\n",
            SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
         SCIP_CALL( SCIPaggregateVars(scip, var1, var2, 1.0, 1.0, 1.0, &infeasible, &redundant, &aggregated) );

         /* evaluate aggregation result */
         if( infeasible )
         {
            SCIPdebugMsg(scip, "linking constraint <%s>: infeasible aggregation <%s> + <%s> == 1\n",
               SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         if( aggregated )
            (*naggrvars)++;
      }

      /* apply real bound to binary variables */
      SCIP_CALL( processRealBoundChg(scip, cons, &cutoff, nchgbds, &mustcheck) );

      /* tightened linking variable */
      SCIP_CALL( tightenedLinkvar(scip, cons, consdata, &cutoff, nchgbds) );

      /* remove the trailing and leeading binary variable which are fixed to zero */
      SCIP_CALL( removeFixedBinvars(scip, conshdlrdata->eventhdlr, cons) );

      /* fix the linking variable to the only remaining value and the corresponding binary variable to 1.0 */
      if( ! cutoff && consdata->nbinvars == 1 )
      {
         SCIP_VAR* linkvar;
         SCIP_VAR* binvar;
         SCIP_Real val;

         linkvar = consdata->linkvar;
         binvar = consdata->binvars[0];
         val = consdata->vals[0];

         SCIPdebugMsg(scip, "linking constraint <%s>: fix <%s> to %16.9g as only one binary variable remains",
                        SCIPconsGetName(cons), SCIPvarGetName(linkvar), val);

         SCIP_CALL( SCIPfixVar(scip, binvar, 1.0, &infeasible, &fixed) );
         assert(fixed);
         ++(*nfixedvars);

         if( ! infeasible )
         {
            SCIP_CALL( SCIPfixVar(scip, linkvar, val, &infeasible, &fixed) );
            assert(fixed);
            ++(*nfixedvars);
         }
         cutoff = infeasible;

         SCIP_CALL(SCIPdelCons(scip, cons));
         ++(*ndelconss);
      }

      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* remember the first changed constraint to begin the next redundancy round with */
      if( firstchange == INT_MAX )
         firstchange = c;

      /* remember the first and last constraints for which we have to add the clique information */
      if( !consdata->cliqueadded && consdata->nbinvars >= 2 )
      {
         if( firstclique == INT_MAX )
            firstclique = c;
         lastclique = c;
      }
   }

   /* add clique and implication information */
   for( c = firstclique; c < lastclique && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;

      assert(*result != SCIP_CUTOFF);

      cons = conss[c];
      assert(cons != NULL);

      /* ignore deleted constraints */
      if( !SCIPconsIsActive(cons) )
         continue;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( !consdata->cliqueadded && consdata->nbinvars >= 3 )
      {
         /* add set partitioning condition as clique */
         int ncliquebdchgs;

         SCIP_CALL( SCIPaddClique(scip, consdata->binvars, NULL, consdata->nbinvars, TRUE, &infeasible, &ncliquebdchgs) );
         *nchgbds += ncliquebdchgs;

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         consdata->cliqueadded = TRUE;
      }
   }

#if 0
   /* transfer aggregated linking variables to the corresponding binary variables */
   assert(conshdlrdata->varmap != NULL);
   SCIP_CALL( aggregateVariables(scip, conshdlrdata->varmap, conss, nconss, naggrvars, &cutoff) );
#endif

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( oldndelconss < *ndelconss || oldnfixedvars < *nfixedvars || oldnchgbds < *nchgbds || oldnaggrvars < *naggrvars)
      *result = SCIP_SUCCESS;

   return SCIP_OKAY; /*lint !e438*/
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR* linkvar;
   int v;

   SCIPdebugMsg(scip, "conflict resolving method of " CONSHDLR_NAME " constraint handler\n");

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   linkvar = consdata->linkvar;
   assert(linkvar != NULL);

   *result = SCIP_DIDNOTFIND;

   if( inferinfo == -1 )
   {
      /* we have to resolve a fixing of a binary variable which was done due to fixed binary variables */
      assert(SCIPvarIsBinary(infervar));
      assert(SCIPisFeasEQ(scip, SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, FALSE)));
      assert(SCIPisFeasEQ(scip, SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, FALSE)));

      if( boundtype == SCIP_BOUNDTYPE_UPPER )
      {
         /* we fixed the binary variable to zero since one of the other binary variable was fixed to one (set
          * partitioning condition)
          */
         assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5);

         for( v = 0; v < consdata->nbinvars; ++v )
         {
            if( SCIPgetVarLbAtIndex(scip, consdata->binvars[v], bdchgidx, FALSE) > 0.5 )
            {
               SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[v]) );
               break;
            }
         }
         assert(v < consdata->nbinvars);
      }
      else
      {
         /* we fixed the binary variable to one since all other binary variable were fixed to zero */
         assert(boundtype == SCIP_BOUNDTYPE_LOWER);
         assert(SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) > 0.5);

         for( v = 0; v < consdata->nbinvars; ++v )
         {
            if( consdata->binvars[v] != infervar )
            {
               /* the reason variable must be assigned to zero */
               assert(SCIPgetVarUbAtIndex(scip, consdata->binvars[v], bdchgidx, FALSE) < 0.5);
               SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[v]) );
            }
         }
      }
   }
   else if( inferinfo == -2 )
   {
      /* we have to resolve a fixing of a binary variable which was done due to the linking variable lower bound */
      assert(SCIPvarIsBinary(infervar));
      assert(SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5);
      assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5); /*@repair: neu*/
      assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, FALSE) > 0.5); /*@repair: neu*/
      assert( SCIPisFeasEQ(scip, SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, FALSE)) );

      SCIP_CALL( SCIPaddConflictLb(scip, linkvar, bdchgidx) );
   }
   else if( inferinfo == -3 )
   {
      /* we have to resolve a fixing of a binary variable which was done due to the linking variable upper bound */
      assert(SCIPvarIsBinary(infervar));
      assert(SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5);
      assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5);
      assert(SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, FALSE) > 0.5);
      assert( SCIPisFeasEQ(scip, SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, FALSE)) );

      SCIP_CALL( SCIPaddConflictUb(scip, linkvar, bdchgidx) );
   }
   else if( inferinfo == -4 )
   {
      SCIP_VAR** binvars;
      SCIP_Real* vals;
      SCIP_Real lb;
      int nbinvars;
      int b;

      /* we tightened the lower bound of the linking variable due the fixing of the corresponding binary variable to zero */
      assert(infervar == linkvar);
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);

      binvars = consdata->binvars;
      nbinvars = consdata->nbinvars;
      vals = consdata->vals;

      /* get propagated lower bound */
      lb = SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, TRUE);

      for( b = 0;  b < nbinvars; ++b )
      {
         if( vals[b] >= lb )
            break;

         assert(SCIPvarGetUbLocal(binvars[b]) < 0.5);
         SCIP_CALL( SCIPaddConflictBinvar(scip, binvars[b]) );
      }
   }
   else if( inferinfo == -5 )
   {
      SCIP_VAR** binvars;
      SCIP_Real* vals;
      SCIP_Real ub;
      int nbinvars;
      int b;

      /* we tightened the upper bound of the linking variable due the fixing of the corresponding binary variable two zero */

      assert(infervar == linkvar);
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);

      binvars = consdata->binvars;
      nbinvars = consdata->nbinvars;
      vals = consdata->vals;

      /* get old and new upper bound */
      ub = SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, TRUE);

      /* resolve tightening of upper bound of the linking variable by binary variables */
      for( b = nbinvars - 1; b >= 0; --b )
      {
         if( vals[b] <= ub )
            break;

         SCIP_CALL( SCIPaddConflictBinvar(scip, binvars[b]) );
      }
   }
   else if( inferinfo == -6 )
   {
      /* we fixed a binary variable to one since the linking variable was fixed */
      assert(SCIPvarIsBinary(infervar));
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);
      assert( SCIPisFeasEQ(scip, SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPgetVarUbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, TRUE), SCIPgetVarLbAtIndex(scip, linkvar, bdchgidx, FALSE)) );

      assert( !SCIPisFeasEQ(scip, SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE), SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, FALSE))  );

      SCIP_CALL( SCIPaddConflictLb(scip, linkvar, bdchgidx) );
      SCIP_CALL( SCIPaddConflictUb(scip, linkvar, bdchgidx) );
   }
   else
   {
      /* we fixed the linking variable to (vals[inferinfo]) since the corresponding binary variable was fixed to one */
      assert(infervar == linkvar);
      assert(inferinfo >= 0);
      assert(inferinfo < consdata->nbinvars);
      assert(SCIPisEQ(scip, consdata->vals[inferinfo], SCIPgetVarUbAtIndex(scip, consdata->linkvar, bdchgidx, TRUE))
          || SCIPisEQ(scip, consdata->vals[inferinfo], SCIPgetVarLbAtIndex(scip, consdata->linkvar, bdchgidx, TRUE)));

      assert(SCIPgetVarLbAtIndex(scip, consdata->binvars[inferinfo], bdchgidx, FALSE) > 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[inferinfo]) );
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int b;

   assert(locktype == SCIP_LOCKTYPE_MODEL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock linking variable in both directions */
   SCIP_CALL( SCIPaddVarLocksType(scip, consdata->linkvar, locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );

   /* look binary variables in both directions */
   for( b = 0; b < consdata->nbinvars; ++b )
   {
      SCIP_CALL( SCIPaddVarLocksType(scip, consdata->binvars[b], locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}

/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveLinking)
{  /*lint --e{715}*/
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPconsIsTransformed(cons));

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPisNLPConstructed(scip) )
   {
      SCIP_CALL( addNlrow(scip, cons) );
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPconsIsTransformed(cons));

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* remove row from NLP, if still in solving
    * if we are in exitsolve, the whole NLP will be freed anyway
    */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && consdata->nlrow1 != NULL )
   {
      assert(consdata->nlrow2 != NULL);
      SCIP_CALL( SCIPdelNlRow(scip, consdata->nlrow1) );
      SCIP_CALL( SCIPdelNlRow(scip, consdata->nlrow2) );
   }

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableLinking)
{  /*lint --e{715}*/
#if 0
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nbinvars <= 1 )
   {
      SCIP_CALL( SCIPdisableCons(scip, cons) );
      assert(consdata->nbinvars == 0 || SCIPvarGetLbGlobal(consdata->binvars[0]) > 0.5);
   }
   else if( conshdlrdata->linearize )
   {
      SCIP_CALL( consdataLinearize(scip, cons, consdata) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
   }
#endif
   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintLinking)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** binvars;
   SCIP_VAR* linkvar;
   SCIP_Real* vals;
   const char* consname;
   int nbinvars;
   int v;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a linking constraint\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }

   (*valid) = TRUE;

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get number of binary variables, linking variables  */
   nbinvars = sourceconsdata->nbinvars;
   linkvar = sourceconsdata->linkvar;

   /* duplicate variable array */
   if( nbinvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &binvars, sourceconsdata->binvars, nbinvars) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, &vals, sourceconsdata->vals, nbinvars) );
   }
   else
   {
      binvars = NULL;
      vals = NULL;
   }

   /* get copy for the binary variables */
   for( v = 0; v < nbinvars && *valid; ++v )
   {
      assert(binvars != NULL); /* for flexelint */
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, binvars[v], &binvars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || binvars[v] != NULL);
   }

   /* copy the linking variable */
   if( *valid )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, linkvar, &linkvar, varmap, consmap, global, valid) );
      assert(!(*valid) || linkvar != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      if( name != NULL )
         consname = name;
      else
         consname = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsLinking(scip, cons, consname, linkvar, binvars, vals, nbinvars,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free buffer array */
   if( nbinvars > 0 )
   {
      SCIPfreeBufferArrayNull(scip, &vals);
      SCIPfreeBufferArrayNull(scip, &binvars);
   }

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseLinking)
{  /*lint --e{715}*/
   SCIP_VAR** binvars;
   SCIP_VAR* linkvar;
   SCIP_Real* vals;
   char* endptr;
   int varssize;
   int nbinvars;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   *success = TRUE;

   /* parse linking variable */
   SCIP_CALL( SCIPparseVarName(scip, str, &linkvar, &endptr) );

   if( linkvar == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   str = endptr;

   nbinvars = 0;
   varssize = 16;

   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, varssize) );

   while( *str != '=' )
      ++str;

   /* skip '=' */
   ++str;

   /* skip whitespace */
   while( isspace((int)*str) )
      ++str;

   /* check for the string "no binary variables yet" */
   if( strncmp(str, "no binary variables yet", 24) != 0 )
   {
      int requsize;
      int v;

      /* parse linear sum to get variables and coefficients */
      SCIP_CALL( SCIPparseVarsLinearsum(scip, str, binvars, vals, &nbinvars, varssize, &requsize, &endptr, success) );

      if( *success && requsize > varssize )
      {
         /* realloc buffers and try again */
         varssize = requsize;
         SCIP_CALL( SCIPreallocBufferArray(scip, &binvars, varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &vals, varssize) );

         SCIP_CALL( SCIPparseVarsLinearsum(scip, str, binvars, vals, &nbinvars, varssize, &requsize, &endptr, success) );
         assert(!*success || requsize <= varssize); /* if successful, then should have had enough space now */
      }

      /* check coefficients */
      if( *success )
      {
         /* convert SCIP_Real to integer */
         for( v = 0; v < nbinvars;  ++v )
         {
            if( SCIPisIntegral(scip, vals[v]) )
               vals[v] = SCIPconvertRealToInt(scip, vals[v]);
         }
      }
   }

   if( *success )
   {
      SCIP_CALL( SCIPcreateConsLinking(scip, cons, name, linkvar, binvars, vals, nbinvars,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nbinvars + 1)
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->binvars, consdata->nbinvars);
      vars[consdata->nbinvars] = consdata->linkvar;
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nbinvars + 1;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBinvar)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   eventtype = SCIPeventGetType(event);
   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      consdata->nfixedones++;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      consdata->nfixedones--;
      consdata->firstnonfixed = 0;
      consdata->lastnonfixed = consdata->nbinvars - 1;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      consdata->nfixedzeros++;
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      consdata->firstnonfixed = 0;
      consdata->lastnonfixed = consdata->nbinvars - 1;
      consdata->nfixedzeros--;
      break;
   default:
      SCIPerrorMessage("invalid event type\n");
      return SCIP_INVALIDDATA;
   }
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nbinvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nbinvars);

   /*debugMsg(scip, " -> constraint has %d zero-fixed and %d one-fixed of %d variables\n",
     consdata->nfixedzeros, consdata->nfixedones, consdata->nvars);*/

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for linking constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrLinking(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecBinvar, NULL) );

   /* create linking constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpLinking, consEnfopsLinking, consCheckLinking, consLockLinking,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyLinking, consCopyLinking) );
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveLinking) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveLinking) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteLinking) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableLinking) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolLinking) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolLinking) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeLinking) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsLinking) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsLinking) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreLinking) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpLinking) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseLinking) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolLinking, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintLinking) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropLinking, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropLinking) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpLinking, consSepasolLinking, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransLinking) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxLinking) );

   /* include the linear constraint to linking constraint upgrade in the linear constraint handler */
   /* SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdLinking, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) ); */

   /* add linking constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/linearize", "this constraint will not propagate or separate, linear and setppc are used?",
         &conshdlrdata->linearize, FALSE, DEFAULT_LINEARIZE, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a linking constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             linkvar,            /**< linking variable (continuous or integer) which should be linked */
   SCIP_VAR**            binvars,            /**< binary variables */
   SCIP_Real*            vals,               /**< coefficients of the binary variables */
   int                   nbinvars,           /**< number of binary starting variables */
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
   SCIP_CONSHDLRDATA* conshdlrdata;
   int k;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(linkvar)));
   assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(linkvar)));

   /* find the linking constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linking constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMsg(scip, "create linking constraint for variable <%s> with %d binary variables (SCIP stage %d)\n",
      SCIPvarGetName(linkvar), nbinvars, SCIPgetStage(scip));
   for( k = 0; k < nbinvars; k++ )
   {
      SCIPdebugMsg(scip, "Var %d : <%s>\n", k, SCIPvarGetName(binvars[k]));
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->varmap == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->varmap, SCIPblkmem(scip), HASHSIZE_BINVARSCONS) );
   }
   assert(conshdlrdata->varmap != NULL);

   /* check if the linking for the requests linking variable already exists */
   assert(!SCIPhashmapExists(conshdlrdata->varmap, getHashmapKey(linkvar)));

   /* create the constraint specific data */
   SCIP_CALL( consdataCreate(scip, conshdlrdata->eventhdlr, &consdata, linkvar, binvars, vals, nbinvars) );

   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   /* create binary variables for the real domain */
   if( nbinvars == 0 )
   {
      SCIP_CALL( consdataCreateBinvars(scip, *cons, consdata, conshdlrdata->eventhdlr, conshdlrdata->linearize) );
   }

   /* insert linking constraint into the hash map */
   SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varmap, getHashmapKey(linkvar), *cons) );
   assert(SCIPhashmapExists(conshdlrdata->varmap, getHashmapKey(linkvar)));

   return SCIP_OKAY;
}

/** creates and captures a linking constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLinking(); all flags can be set via SCIPsetCons<Flagname>-methods in scip.h
 *
 *  @see SCIPcreateConsLinking() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             linkvar,            /**< linking variable (continuous or integer) which should be linked */
   SCIP_VAR**            binvars,            /**< binary variables, or NULL */
   SCIP_Real*            vals,               /**< coefficients of the binary variables */
   int                   nbinvars            /**< number of binary variables */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsLinking(scip, cons, name, linkvar, binvars, vals, nbinvars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** checks if for the given linking variable (continuous or integer) a linking constraint exists */
SCIP_Bool SCIPexistsConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             linkvar             /**< linking variable (continuous or integer) which should be linked */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return (conshdlrdata->varmap != NULL) && SCIPhashmapExists(conshdlrdata->varmap, getHashmapKey(linkvar));
}

/** returns the linking constraint belonging the given linking variable (continuous or integer) or NULL if it does not exist yet */
SCIP_CONS* SCIPgetConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             linkvar             /**< linking variable (continuous or integer) which should be linked */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->varmap != NULL )
      return (SCIP_CONS*) SCIPhashmapGetImage(conshdlrdata->varmap, getHashmapKey(linkvar));
   else
      return NULL;
}

/** returns the linking variable (continuous or integer) of the linking constraint */
SCIP_VAR* SCIPgetLinkvarLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a " CONSHDLR_NAME " constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->linkvar;
}

/** returns the binary variables of the linking constraint */
SCIP_RETCODE SCIPgetBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR***           binvars,            /**< pointer to store the binary variables array pointer */
   int*                  nbinvars            /**< pointer to store the number of returned binary variables */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a " CONSHDLR_NAME " constraint\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->binvars == NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( consdataCreateBinvars(scip, cons, consdata, conshdlrdata->eventhdlr, conshdlrdata->linearize) );
   }

   assert(consdata->binvars != NULL);

   if( binvars != NULL )
      (*binvars) = consdata->binvars;
   if( nbinvars != NULL )
      (*nbinvars) = consdata->nbinvars;

   return SCIP_OKAY;
}

/** returns the number of binary variables of the linking constraint */
int SCIPgetNBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a " CONSHDLR_NAME " constraint\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nbinvars;
}

/** returns the coefficients of the binary variables */
SCIP_Real* SCIPgetValsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a " CONSHDLR_NAME " constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   consdataSort(consdata);

   return consdata->vals;
}

/** return all binary variable information of the linking constraint */
SCIP_RETCODE SCIPgetBinvarsDataLinking(
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR***           binvars,            /**< pointer to store binary variables, or NULL */
   SCIP_Real**           vals,               /**< pointer to store the binary coefficients, or NULL */
   int*                  nbinvars            /**< pointer to store the number of binary variables, or NULL */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a " CONSHDLR_NAME " constraint\n");
      SCIPABORT();
      return SCIP_ERROR;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdataSort(consdata);

   if( binvars != NULL )
      *binvars = consdata->binvars;
   if( vals != NULL )
      *vals = consdata->vals;
   if( nbinvars != NULL )
      *nbinvars = consdata->nbinvars;

   return SCIP_OKAY;
}
