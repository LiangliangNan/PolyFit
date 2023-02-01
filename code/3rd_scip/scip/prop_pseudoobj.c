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

/**@file   prop_pseudoobj.c
 * @brief  Pseudo objective propagator
 * @author Tobias Achterberg
 * @author Stefan Heinz
 *
 * This propagator propagates the objective function using the cutoff bound and the pseudo objective value. The pseudo
 * objective value can be seen as minimum activity of the linear objective function. Using this, this propagator checks
 * if variables with non-zero objective coefficients can exceed the cutoff bound. If this is the case the corresponding
 * bound can be tightened.
 *
 * @todo use the complete implication to initialize the objective implication data structure
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_pseudoobj.h"


#define PROP_NAME              "pseudoobj"
#define PROP_DESC              "pseudo objective function propagator"
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP | SCIP_PROPTIMING_DURINGLPLOOP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY           3000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_PRESOL_PRIORITY   +6000000 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolving method (fast, medium, or exhaustive) */

#define EVENTHDLR_NAME         "pseudoobj"
#define EVENTHDLR_DESC         "bound change event handler for pseudo objective function propagator"

#define DEFAULT_MINUSELESS          100 /**< minimal number of successive non-binary variable propagator whithout a
                                         *   bound reduction before aborted */
#define DEFAULT_MAXVARSFRAC         0.1 /**< maximal fraction of non-binary variables with non-zero objective
                                         *   without a bound reduction before aborted */
#define DEFAULT_PROPFULLINROOT     TRUE /**< do we want to propagate all non-binary variables if we are propagating the root node? */
#define DEFAULT_PROPCUTOFFBOUND    TRUE /**< propagate new cutoff bound directly globally */
#define DEFAULT_FORCE             FALSE /**< should the propagator be forced even if active pricer are present? Note that
                                         *   can be done if it is known that the pseudo objective activity is given by
                                         *   the zero bound for all variables which are currently not present in the
                                         *   problem */
#define DEFAULT_MAXNEWVARS         1000 /**< number of variable added after the propagator is reinitialized? */
#define DEFAULT_PROPUSEIMPLICS     TRUE /**< use implications to strengthen the propagation of binary variable (increasing the objective change)? */
#define DEFAULT_RESPROPUSEIMPLICS  TRUE /**< use implications to strengthen the resolve propagation of binary variable (increasing the objective change)? */
#define DEFAULT_MAXIMPLVARS       50000 /**< maximum number of binary variables the implications are used if turned on (-1: unlimited)? */


/*
 * Data structures
 */

/** implication data structure for objective contributions of a binary variable */
struct SCIP_ObjImplics
{
   SCIP_VAR**            objvars;            /**< variables y in implications y == 0 or y == 1, first we store the
                                              *   implications by x == 0 and second the implications x == 1 */
   SCIP_Real             maxobjchg;          /**< maximum objective contribution if variables x is fixed to zero or one */
   int                   nlbimpls;           /**< number of all implications result through for x == 0 */
   int                   nubimpls;           /**< number of all implications result through for x == 1 */
   int                   size;               /**< size of the objvars array */
};
typedef struct SCIP_ObjImplics SCIP_OBJIMPLICS; /**< implications in the form x == 0 or x == 1 ==> y == 0 or y == 1 for (x and y binary) */


/** propagator data */
struct SCIP_PropData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for global bound change events */
   SCIP_VAR**            minactvars;         /**< binary variables with non-zero objective contribution w.r.t. minimum activity of the objective function */
   SCIP_OBJIMPLICS**     minactimpls;        /**< implication data structure for the binary variables w.r.t. minimum activity */
   SCIP_VAR**            maxactvars;         /**< binary variables with non-zero objective contribution w.r.t. maximum activity of the objective function */
   SCIP_Real*            maxactchgs;         /**< the maximal potential change of the objective if the binary variable
                                              *   is fixed to its best bound w.r.t. maximum activity of the objective function */

   SCIP_VAR**            objintvars;         /**< non-binary variable with non-zero objective coefficient */
   SCIP_HASHTABLE*       addedvars;          /**< hash table used during resolving of a bound change (conflict analysis) */
   SCIP_Real             lastlowerbound;     /**< last lower bound which was propagated */
   SCIP_Real             cutoffbound;        /**< last cutoff bound used for propagation */
   SCIP_Real             glbpseudoobjval;    /**< last global pseudo objective used in presolving */
   SCIP_Real             maxvarsfrac;        /**< maximal fraction of non-binary variables with non-zero objective
                                              *   without a bound reduction before aborted
                                              */
   SCIP_Real             maxpseudoobjact;    /**< maximal global pseudo objective activity */
   int                   maxpseudoobjactinf; /**< number of coefficients contributing with infinite value to maxpseudoobjact */
   int                   nminactvars;        /**< number of binary variables with non-zero objective contribution w.r.t. minimum activity of the objective function */
   int                   nmaxactvars;        /**< number of binary variables with non-zero objective contribution w.r.t. maximum activity of the objective function */
   int                   nobjintvars;        /**< number of non-binary variables with non-zero objective */
   int                   minuseless;         /**< minimal number of successive non-binary variable propagator whithout
                                              *   a bound reduction before aborted
                                              */
   int                   lastvarnum;         /**< last non-binary variable number that was looked at */
   int                   glbfirstnonfixed;   /**< index of first globally non-fixed binary variable in minactvars array */
   int                   maxactfirstnonfixed;/**< index of first globally non-fixed binary variable in maxctvars array */
   int                   firstnonfixed;      /**< index of first locally non-fixed binary variable in minactvars array */
   int                   nnewvars;           /**< counter for counting number of new variables added */
   int                   maxnewvars;         /**< number of variable added after the propagator is reinitialized? */
   int                   maximplvars;        /**< maximum number of binary variables the implications are used if turned on (-1: unlimited)? */
   int                   minactsize;         /**< size of minactvars and minactimpls array */
   int                   maxactsize;         /**< size of maxactvars and maxactchgs array */
   int                   objintvarssize;     /**< size of objintvars array*/
   SCIP_Bool             glbpropagated;      /**< are global domains propagated */
   SCIP_Bool             propfullinroot;     /**< do we want to propagate all non-binary variables if we are propagating the root node */
   SCIP_Bool             propcutoffbound;    /**< propagate new cutoff bound directly globally */
   SCIP_Bool             force;              /**< should the propagator be forced even if active pricer are present? */
   SCIP_Bool             catchvaradded;      /**< do we catch the variable added event? */
   SCIP_Bool             propuseimplics;     /**< use implications to strengthen the propagation of binary variable (increasing the objective change)? */
   SCIP_Bool             respropuseimplics;  /**< use implications to strengthen the resolve propagation of binary variable (increasing the objective change)? */
   SCIP_Bool             initialized;        /**< is propagator data structure initialized */
};

/*
 * Debug methods
 */

#ifndef NDEBUG
/** check that the implications are applied for a globally fixed variable */
static
void checkImplicsApplied(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to check the implications */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* bounds;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Bool varfixing;
   int nvars;
   int v;

   /* check that the given variable is locally or globally fixed */
   assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);

   /* get fixed value */
   varfixing = SCIPvarGetLbGlobal(var) > 0.5;

   vars = SCIPvarGetImplVars(var, varfixing);
   nvars = SCIPvarGetNImpls(var, varfixing);
   bounds = SCIPvarGetImplBounds(var, varfixing);
   boundtypes = SCIPvarGetImplTypes(var, varfixing);

   /* check that each implication was applied */
   for( v = 0; v < nvars; ++v )
   {
      if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_Real lb;

         lb = SCIPvarGetLbGlobal(vars[v]);
         assert(SCIPisLE(scip, lb, bounds[v]));
      }
      else
      {
         SCIP_Real ub;

         assert(boundtypes[v] == SCIP_BOUNDTYPE_UPPER);

         ub = SCIPvarGetLbGlobal(vars[v]);
         assert(SCIPisGE(scip, ub, bounds[v]));
      }
   }
}

/** check if the global fixed indices are correct */
static
void checkGlbfirstnonfixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR* var;
   int v;

   for( v = 0; v < propdata->glbfirstnonfixed; ++v )
   {
      var = propdata->minactvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5);
   }

   for( v = 0; v < propdata->maxactfirstnonfixed; ++v )
   {
      var = propdata->maxactvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5);
   }
}
#endif /* end of debug methods */

/*
 * Comparer
 */

/** compares objective implications w.r.t. their maximum contribution */
static
SCIP_DECL_SORTPTRCOMP(objimplicsComp)
{
   SCIP_OBJIMPLICS* objimplics1;
   SCIP_OBJIMPLICS* objimplics2;

   objimplics1 = (SCIP_OBJIMPLICS*)elem1;
   objimplics2 = (SCIP_OBJIMPLICS*)elem2;

   if( objimplics1->maxobjchg > objimplics2->maxobjchg )
      return +1;

   if( objimplics1->maxobjchg < objimplics2->maxobjchg )
      return -1;

   return 0;
}

/** compare variables w.r.t.
 *  (i)   the absolute value the objective coefficient;
 *  (ii)  the locks which indicate most effect -- for the variables with a positive (negative) objective coefficient the
 *        down (up) lock is used since this lock indicates that tightened of the upper (lower) bound will triegger
 *        further domain propagations;
 *  (iii) the other locks;
 *  (iv)  variable problem index;
 */
static
SCIP_DECL_SORTPTRCOMP(varCompObj)
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int locks1;
   int locks2;

   var1 = (SCIP_VAR*)elem1;
   var2 = (SCIP_VAR*)elem2;

   assert(SCIPvarGetObj(var1) != 0.0);
   assert(SCIPvarGetObj(var2) != 0.0);

   /* first criteria is the absolute value of objective coefficient */
   if( REALABS(SCIPvarGetObj(var1)) < REALABS(SCIPvarGetObj(var2)) )
      return -1;
   else if( REALABS(SCIPvarGetObj(var1)) > REALABS(SCIPvarGetObj(var2)) )
      return +1;

   /* second criteria the locks which indicate most effect */
   if( SCIPvarGetObj(var1) > 0.0 )
      locks1 = SCIPvarGetNLocksDown(var1);
   else
      locks1 = SCIPvarGetNLocksUp(var1);

   if( SCIPvarGetObj(var2) > 0.0 )
      locks2 = SCIPvarGetNLocksDown(var2);
   else
      locks2 = SCIPvarGetNLocksUp(var2);

   if( locks1 < locks2 )
      return -1;
   if( locks1 > locks2 )
      return 1;

   /* third criteria the other locks */
   if( SCIPvarGetObj(var1) > 0.0 )
      locks1 = SCIPvarGetNLocksUp(var1);
   else
      locks1 = SCIPvarGetNLocksDown(var1);

   if( SCIPvarGetObj(var2) >  0.0 )
      locks2 = SCIPvarGetNLocksUp(var2);
   else
      locks2 = SCIPvarGetNLocksDown(var2);

   if( locks1 < locks2 )
      return -1;
   if( locks1 > locks2 )
      return 1;

   /* forth criteria use the problem index */
   return SCIPvarCompare(var1, var2);
}

/** hash key retrieval function for cliques*/
static
SCIP_DECL_HASHGETKEY(cliqueGetHashkey)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the cliques are equal */
static
SCIP_DECL_HASHKEYEQ(cliqueIsHashkeyEq)
{  /*lint --e{715}*/
   if( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(cliqueGetHashkeyVal)
{  /*lint --e{715}*/
   return SCIPcliqueGetId((SCIP_CLIQUE*) key);
}

/*
 * methods for SCIP_OBJIMPLICS data structure
 */

/** creates an objective implication data structure, fixes (globally) variables which are implied by lower and upper
 *  bound fixing, and clears the collected arrays for lower and upper bound
 */
static
SCIP_RETCODE objimplicsCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJIMPLICS**     objimplics,         /**< pointer to objective implication data structure */
   SCIP_VAR**            objvars,            /**< objective contributor variables, or NULL */
   SCIP_HASHMAP*         binobjvarmap,       /**< hash map mapping binary variables with none-zero objective to position in collected variables arrays, or NULL */
   SCIP_Bool*            collectedlbvars,    /**< temporary buffer to mark collected variables for lower bound fixing, or NULL */
   SCIP_Bool*            collectedubvars,    /**< temporary buffer to mark collected variables for upper bound fixing, or NULL */
   SCIP_Real             maxlbobjchg,        /**< maximum objective contributor if variables id fixed to zero */
   SCIP_Real             maxubobjchg,        /**< maximum objective contributor if variables id fixed to one */
   int                   nlbimpls,           /**< number of variables contributing to to lower bound fix */
   int                   nubimpls            /**< number of variables contributing to to upper bound fix */
   )

{
   assert(scip != NULL);
   assert(objimplics != NULL);
   assert(!SCIPisNegative(scip, maxlbobjchg));
   assert(!SCIPisNegative(scip, maxubobjchg));

   /* allocate block memory for the implication data structure */
   SCIP_CALL( SCIPallocBlockMemory(scip, objimplics) );

   if( nlbimpls + nubimpls == 0 )
   {
      assert(nlbimpls == 0);
      assert(nubimpls == 0);
      (*objimplics)->objvars = NULL;
      (*objimplics)->maxobjchg = 0.0;
      (*objimplics)->nlbimpls = 0;
      (*objimplics)->nubimpls = 0;
      (*objimplics)->size = 0;
   }
   else
   {
      SCIP_VAR* var;
      int nvars;
      int pos;
      int v;

      assert(objvars != NULL);
      assert(binobjvarmap != NULL);
      assert(collectedlbvars != NULL);
      assert(collectedubvars != NULL);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*objimplics)->objvars, nlbimpls + nubimpls) );
      (*objimplics)->size = nlbimpls + nubimpls;

      nvars = 0;

      for( v = 0; v < nlbimpls; ++v )
      {
         var = objvars[v];
         assert(var != NULL);
         assert(!SCIPisZero(scip, SCIPvarGetObj(var)));

         assert(SCIPhashmapExists(binobjvarmap, var));
         pos = (int)(size_t)SCIPhashmapGetImage(binobjvarmap, (void*)var);
         assert(pos > 0);
         assert(collectedlbvars[pos]);

         if( collectedubvars[pos] )
         {
            SCIP_Bool infeasible;
            SCIP_Bool tightened;

            if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
            {
               SCIPdebugMsg(scip, "fix variables <%s> to 1.0 due to implications\n", SCIPvarGetName(var));

               SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, 1.0, FALSE, &infeasible, &tightened) );
               maxlbobjchg -= SCIPvarGetObj(var);
            }
            else
            {
               SCIPdebugMsg(scip, "fix variables <%s> to 0.0 due to implications\n", SCIPvarGetName(var));

               SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, 0.0, FALSE, &infeasible, &tightened) );
               maxlbobjchg += SCIPvarGetObj(var);
            }
            assert(!infeasible);
            assert(tightened);
         }
         else
         {
            (*objimplics)->objvars[nvars] = var;
            nvars++;
         }
         collectedlbvars[pos] = FALSE;
      }
      (*objimplics)->nlbimpls = nvars;

      for( v = 0; v < nubimpls; ++v )
      {
         var = objvars[nlbimpls + v];
         assert(var != NULL);
         assert(!SCIPisZero(scip, SCIPvarGetObj(var)));

         assert(SCIPhashmapExists(binobjvarmap, var));
         pos = (int)(size_t)SCIPhashmapGetImage(binobjvarmap, (void*)var);
         assert(pos > 0);
         assert(collectedubvars[pos]);

         if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
         {
            if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
               maxubobjchg -= SCIPvarGetObj(var);
            else
               maxubobjchg += SCIPvarGetObj(var);
         }
         else
         {
            (*objimplics)->objvars[nvars] = var;
            nvars++;
         }
         collectedubvars[pos] = FALSE;
      }
      (*objimplics)->nubimpls = nvars - (*objimplics)->nlbimpls;

      /* capture the variables */
      for( v = 0; v < nvars; ++v )
      {
         assert(SCIPvarIsBinary((*objimplics)->objvars[v]));
         assert(!SCIPisZero(scip, SCIPvarGetObj((*objimplics)->objvars[v])));
         SCIP_CALL( SCIPcaptureVar(scip, (*objimplics)->objvars[v]) );
      }
   }

   (*objimplics)->maxobjchg = MAX(maxlbobjchg, maxubobjchg);

   return SCIP_OKAY;
}

/** frees an objective implication data structure */
static
SCIP_RETCODE objimplicsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJIMPLICS**     objimplics          /**< pointer to objective implication data structure */
   )
{
   int v;

   assert(scip != NULL);
   assert(objimplics != NULL);
   assert(*objimplics != NULL);

   /* release all variables */
   for( v = 0; v < (*objimplics)->nlbimpls + (*objimplics)->nubimpls; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*objimplics)->objvars[v]) );
   }

   /* free objective variable array */
   SCIPfreeBlockMemoryArrayNull(scip, &(*objimplics)->objvars, (*objimplics)->size);

   /* frees block memory for the implication data structure */
   SCIPfreeBlockMemory(scip, objimplics);

   return SCIP_OKAY;
}

/** remove the given variable at the given pos from the objective implication data structure */
static
SCIP_RETCODE objimplicsDelPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJIMPLICS*      objimplics,         /**< objective implication data structure */
   int                   pos                 /**< position */
   )
{
   assert(0 <= pos);
   assert(pos < objimplics->nlbimpls + objimplics->nubimpls);

   SCIPdebugMsg(scip, "variable <%s> ", SCIPvarGetName(objimplics->objvars[pos]));

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &objimplics->objvars[pos]) );

   /* copy last lower bound variable to that position */
   if( pos < objimplics->nlbimpls )
   {
      objimplics->nlbimpls--;
      assert(objimplics->nlbimpls >= 0);

      /* copy last lower bound variable to that position */
      objimplics->objvars[pos] = objimplics->objvars[objimplics->nlbimpls];

      /* copy last upper bound variable to open slot */
      objimplics->objvars[objimplics->nlbimpls] = objimplics->objvars[objimplics->nlbimpls + objimplics->nubimpls];

      SCIPdebugMsgPrint(scip, "remove lower bound implication\n");
   }
   else
   {
      objimplics->nubimpls--;
      assert(objimplics->nubimpls >= 0);

      /* copy last upper bound variable to that position */
      objimplics->objvars[pos] = objimplics->objvars[objimplics->nlbimpls + objimplics->nubimpls];

      SCIPdebugMsgPrint(scip, "remove upper bound implication\n");
   }

   return SCIP_OKAY;
}

/*
 * Local methods
 */


/** catch bound change events if the variable has a non-zero objective coefficient to check if the maximum activity of
 *  the objective function changed
 */
static
SCIP_RETCODE catchObjEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for global bound change events */
   SCIP_VAR*             var                 /**< variable for which the event should be dropped */
   )
{
   SCIP_Real objval;

   assert(propdata != NULL);
   assert(eventhdlr != NULL);

   objval = SCIPvarGetObj(var);

   if( !SCIPisZero(scip, objval) )
   {
      if( objval > 0.0 )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GUBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
      }
      else
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GLBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
      }
   }

   return SCIP_OKAY;
}

/** drop variable event w.r.t. objective coefficient */
static
SCIP_RETCODE dropObjEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for global bound change events */
   SCIP_VAR*             var                 /**< variable for which the event should be dropped */
   )
{
   SCIP_Real objval;

   assert(propdata != NULL);
   assert(eventhdlr != NULL);

   objval = SCIPvarGetObj(var);

   /* drop bound change event */
   if( !SCIPisZero(scip, objval) )
   {
      if( objval > 0.0 )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GUBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }
      else
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GLBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }
   }
   return SCIP_OKAY;
}

/** drop all variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_VAR* var;
   int k;

   assert(scip != NULL);
   assert(propdata != NULL);

   eventhdlr = propdata->eventhdlr;
   assert(eventhdlr != NULL);

   /* drop all events and release variables */
   for( k = 0; k < propdata->nminactvars; ++k )
   {
      var =  propdata->minactvars[k];
      assert(var != NULL);
      assert(SCIPvarIsBinary(var));

      /* drop bound relax event which is caught for all binary variables which are used for propagation the objective
       * function via the minimum activity of the objective function
       */
      SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDRELAXED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* release variables */
   for( k = 0; k < propdata->nmaxactvars; ++k )
   {
      var = propdata->maxactvars[k];
      assert(var != NULL);
      assert(SCIPvarIsBinary(var));

      /* drop events which are needed for evaluating the maximum activity of the objective function */
      SCIP_CALL( dropObjEvent(scip, propdata, eventhdlr, var) );

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* drop all events and release variables */
   for( k = 0; k < propdata->nobjintvars; ++k )
   {
      var = propdata->objintvars[k];
      assert(var != NULL);

      /* drop events which are needed for evaluating the maximum activity of the objective function */
      SCIP_CALL( dropObjEvent(scip, propdata, eventhdlr, var) );

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   return SCIP_OKAY;
}

/** reset propagatore data structure */
static
void propdataReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   propdata->minactvars = NULL;
   propdata->minactimpls = NULL;
   propdata->maxactvars = NULL;
   propdata->maxactchgs = NULL;
   propdata->objintvars = NULL;
   propdata->nminactvars = 0;
   propdata->nmaxactvars = 0;
   propdata->nobjintvars = 0;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;
   propdata->lastvarnum = -1;
   propdata->glbpropagated = FALSE;
   propdata->cutoffbound = SCIP_INVALID;
   propdata->lastlowerbound = -SCIP_INVALID;
   propdata->glbpseudoobjval = -SCIP_INVALID;
   propdata->glbfirstnonfixed = 0;
   propdata->maxactfirstnonfixed = 0;
   propdata->firstnonfixed = 0;
   propdata->nnewvars = 0;
   propdata->minactsize = 0;
   propdata->maxactsize = 0;
   propdata->objintvarssize = 0;
   propdata->catchvaradded = FALSE;
   propdata->initialized = FALSE;
}

/** free propagator data */
static
SCIP_RETCODE propdataExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int v;

   if( !propdata->initialized )
      return SCIP_OKAY;

   if( propdata->addedvars != NULL )
      SCIPhashtableFree(&propdata->addedvars);

   /* drop events for the variables */
   SCIP_CALL( dropVarEvents(scip, propdata) );

   for( v = 0; v < propdata->nminactvars; ++v )
   {
      SCIP_CALL( objimplicsFree(scip, &(propdata->minactimpls[v])) );
   }

   /* free memory for non-zero objective variables */
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->minactvars, propdata->minactsize);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->minactimpls, propdata->minactsize);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->maxactvars, propdata->maxactsize);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->maxactchgs, propdata->maxactsize);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->objintvars, propdata->objintvarssize);

   /* reset propagator data structure */
   propdataReset(scip, propdata);

   return SCIP_OKAY;
}

/** returns the objective change for the given binary variable */
static
SCIP_Real getVarObjchg(
   SCIP_VAR*             var,                /**< variable to get objective change for */
   SCIP_BOUNDTYPE        boundtype,          /**< bound type to consider */
   SCIP_BOUNDTYPE        bound               /**< fixing bound */
   )
{
   assert(SCIPvarIsBinary(var));
   assert((int)SCIP_BOUNDTYPE_LOWER == 0);
   assert((int)SCIP_BOUNDTYPE_UPPER == 1);

   /* collect contribution of variable itself */
   return (SCIP_Real)((int)bound - (int)(boundtype == SCIP_BOUNDTYPE_UPPER)) * SCIPvarGetObj(var);
}

/** returns the objective change provided by the implication variable by fixing it to the given bound
 *  w.r.t. minimum activity of the objective function; additionally it collects all contributors for that objective
 *  change;
 */
static
SCIP_Real collectMinactImplicVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_HASHMAP*         binobjvarmap,       /**< hash map mapping binary variables with none-zero objective to position in collected variables arrays */
   SCIP_Bool*            collectedvars,      /**< temporary buffer to mark collected variables */
   int                   nbinobjvars,        /**< number of binary variables with non-zero objective coefficient */
   SCIP_VAR**            contributors,       /**< array to store the contributors */
   int*                  ncontributors       /**< pointer to store number of contributor to the objective contribution */
   )
{
   SCIP_Real objval;
   int pos;

   assert(scip != NULL);
   assert(var != NULL);
   assert(binobjvarmap != NULL);
   assert(collectedvars != NULL);
   assert(contributors != NULL);
   assert(ncontributors != NULL);

   /* ignore global fixed variables */
   if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
      return 0.0;

   objval = SCIPvarGetObj(var);

   /* ignore variables with zero objective coefficient */
   if( SCIPisZero(scip, objval) )
      return 0.0;

   assert(SCIPhashmapExists(binobjvarmap, var));
   pos = (int)(size_t)SCIPhashmapGetImage(binobjvarmap, var);
   assert(pos > 0);

   /* check if the variables was already collected through other cliques */
   if( collectedvars[pos] )
      return 0.0;

   /* collect variable */
   assert(*ncontributors < nbinobjvars);
   contributors[*ncontributors] = var;
   (*ncontributors)++;

   /* mark variable to be collected */
   collectedvars[pos] = TRUE;

   /* return the absolute value of the objective coefficient as constriction */
   return REALABS(objval);
}

#define MAX_CLIQUELENGTH 50
/** returns the objective change provided by the implications of the given variable by fixing it to the given bound
 *  w.r.t. minimum activity of the objective function; additionally it collects all contributors for that objective
 *  change;
 *
 *  Let I(0) and I(1) be all implications of the given variable which follow by fixing it to given bound and evaluate to
 *  fixing the implication variable to zero (I(0)) or one (I(1)), respectively. The objective change provided by the
 *  implications are:
 *
 *  \f[
 *  \displaystyle
 *  sum_{x\in I(1)} (1 - \mbox{bestbound}(x)) \cdot \mbox{objval}(x) - sum_{x\in I(1)} \mbox{bestbound}(x) \cdot \mbox{objval}(x)
 *  =
 *  sum_{x\in I(0) \cup I(1)} (\mbox{impliedbound}(x) - \mbox{bestbound}(x)) \cdot \mbox{objval}(x)
 *  \f]
 */
static
SCIP_RETCODE collectMinactImplicVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_HASHMAP*         binobjvarmap,       /**< hash map mapping binary variables with none-zero objective to position in collected variables arrays */
   SCIP_Bool*            collectedvars,      /**< temporary buffer to mark collected variables */
   int                   nbinobjvars,        /**< number of binary variables with non-zero objective coefficient */
   SCIP_VAR**            contributors,       /**< array to store the contributors */
   SCIP_HASHTABLE*       uselesscliques,     /**< hash table to store useless cliques, or NULL */
   int*                  ncontributors,      /**< pointer to store number of contributor to the objective contribution */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_CLIQUE* clique;
   SCIP_VAR** vars;
   SCIP_VAR* implvar;
   SCIP_Bool* values;
   SCIP_Bool varfixing;
   int nbinvars;
   int ncliques;
   int c;
   int v;

   assert(SCIPvarIsBinary(var));
   assert(SCIPvarGetLbGlobal(var) < 0.5);
   assert(SCIPvarGetUbGlobal(var) > 0.5);
   assert(bound == SCIP_BOUNDTYPE_LOWER || bound == SCIP_BOUNDTYPE_UPPER);
   assert(objchg != NULL);
   assert(contributors != NULL);
   assert(ncontributors != NULL);
   assert(*ncontributors == 0);

   assert((SCIP_Bool)SCIP_BOUNDTYPE_LOWER == FALSE);
   assert((SCIP_Bool)SCIP_BOUNDTYPE_UPPER == TRUE);
   varfixing = (SCIP_Bool)bound;

   cliques = SCIPvarGetCliques(var, varfixing);
   ncliques = SCIPvarGetNCliques(var, varfixing);

   if( uselesscliques == NULL )
      return SCIP_INVALIDDATA;

#ifndef NDEBUG
   /* check that the marker array is reset */
   for( c = 0; c < nbinobjvars; ++c )
      assert(collectedvars[c] == FALSE);
#endif

   /* collect all implication which are given via cliques */
   for( c = 0; c < ncliques; ++c )
   {
      SCIP_Bool useless;

      clique = cliques[c];
      assert(clique != NULL);

      /* check if the clique was previously detected to be useless with respect to minimum activity */
      if( SCIPhashtableExists(uselesscliques, (void*)clique) )
         continue;

      nbinvars = SCIPcliqueGetNVars(clique);

      if( nbinvars > MAX_CLIQUELENGTH )
      {
         SCIP_CALL( SCIPhashtableInsert(uselesscliques, (void*)clique) );
         continue;
      }

      vars = SCIPcliqueGetVars(clique);
      values = SCIPcliqueGetValues(clique);
      useless = TRUE;

      for( v = 0; v < nbinvars; ++v )
      {
         implvar = vars[v];
         assert(implvar != NULL);

         if( implvar == var )
         {
            /* check if the clique is useful at all */
            if( useless )
            {
               SCIP_Real objval;

               objval = SCIPvarGetObj(var);

               if( varfixing == (SCIP_Bool)SCIPvarGetBestBoundType(var) && !SCIPisZero(scip, objval) )
                  useless = FALSE;
            }
         }
         else if( values[v] == (SCIP_Bool)SCIPvarGetBestBoundType(implvar) )
         {
            useless = FALSE;
            (*objchg) += collectMinactImplicVar(scip, implvar, binobjvarmap, collectedvars, nbinobjvars, contributors, ncontributors);
         }
      }

      /* if the clique is useless store it in the hash table to skip it later */
      if( useless )
      {
         assert(!SCIPhashtableExists(uselesscliques, (void*)clique));
         SCIP_CALL( SCIPhashtableInsert(uselesscliques, (void*)clique) );
      }
   }

   return SCIP_OKAY;
}

/** returns the objective change provided by the implications of the given variable by fixing it to the given bound
 *  w.r.t. minimum activity of the objective function
 *
 *  Let I(0) and I(1) be all implications of the given variable which follow by fixing it to given bound and evaluate to
 *  fixing the implication variable to zero (I(0)) or one (I(1)), respectively. The objective change provided by the
 *  implications are:
 *
 *  \f[
 *  \displaystyle
 *  sum_{x\in I(1)} (1 - \mbox{bestbound}(x)) \cdot \mbox{objval}(x) - sum_{x\in I(1)} \mbox{bestbound}(x) \cdot \mbox{objval}(x)
 *  =
 *  sum_{x\in I(0) \cup I(1)} (\mbox{impliedbound}(x) - \mbox{bestbound}(x)) \cdot \mbox{objval}(x)
 *  \f]
 *
 *  This can be done w.r.t. global variable bounds (local == FALSE), w.r.t. local variable bounds (local == TRUE &&
 *  bdchgidx == NULL), and w.r.t. given time stamp (local == TRUE && bdchgidx != NULL)
 */
static
SCIP_RETCODE getMinactImplicObjchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_OBJIMPLICS*      objimplics,         /**< objective implication data for the given variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, or NULL */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_Bool             local,              /**< propagate local bounds, otherwise global bounds */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   SCIP_VAR* implvar;
   SCIP_Bool lb;
   SCIP_Bool ub;
   int nbinvars;
   int v;

   assert(SCIPvarIsBinary(var));
   assert(!local || SCIPgetVarLbAtIndex(scip, var, bdchgidx, FALSE) < 0.5);
   assert(!local || SCIPgetVarUbAtIndex(scip, var, bdchgidx, FALSE) > 0.5);
   assert(SCIPvarGetLbGlobal(var) < 0.5);
   assert(SCIPvarGetUbGlobal(var) > 0.5);
   assert(bound == SCIP_BOUNDTYPE_LOWER || bound == SCIP_BOUNDTYPE_UPPER);

   if( bound == SCIP_BOUNDTYPE_LOWER )
   {
      v = 0;
      nbinvars = objimplics->nlbimpls;
   }
   else
   {
      assert(bound == SCIP_BOUNDTYPE_UPPER);
      v = objimplics->nlbimpls;
      nbinvars = objimplics->nlbimpls + objimplics->nubimpls;
   }

   /* loop over all implications */
   while( v < nbinvars )
   {
      implvar = objimplics->objvars[v];
      assert(implvar != NULL);
      assert(!SCIPisZero(scip, SCIPvarGetObj(implvar)));

      if( local )
      {
         lb = SCIPgetVarLbAtIndex(scip, implvar, bdchgidx, FALSE) > 0.5;
         ub = SCIPgetVarUbAtIndex(scip, implvar, bdchgidx, FALSE) > 0.5;

         /* check if variable is fixed */
         if( lb == TRUE || ub == FALSE )
         {
            v++;
            continue;
         }
      }
      else
      {
         lb = SCIPvarGetLbGlobal(implvar) > 0.5;
         ub = SCIPvarGetUbGlobal(implvar) > 0.5;

         /* check if variable is global fixed; if so remove it from the objective implication data structure and
          * continue with the next candidate
          */
         if( lb == TRUE || ub == FALSE )
         {
            SCIP_CALL( objimplicsDelPos(scip, objimplics, v) );
            nbinvars--;
            continue;
         }
      }

      assert(SCIPvarGetObj(implvar) > 0.0 || SCIPvarsHaveCommonClique(var, (SCIP_Bool) bound, implvar, TRUE, TRUE));
      assert(SCIPvarGetObj(implvar) < 0.0 || SCIPvarsHaveCommonClique(var, (SCIP_Bool) bound, implvar, FALSE, TRUE));

      /* add objective change */
      (*objchg) += REALABS(SCIPvarGetObj(implvar));
      v++;
   }

   return SCIP_OKAY;
}

/** computes for the given binary variable the objective contribution by fixing it to given bound w.r.t. minimum
 *  activity of the objective function; additionally it collects all contributors for that objective change;
 */
static
SCIP_RETCODE collectMinactObjchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_HASHMAP*         binobjvarmap,       /**< hash map mapping binary variables with none-zero objective to position in collected variables arrays */
   SCIP_Bool*            collectedvars,      /**< temporary buffer to mark collected variables */
   int                   nbinobjvars,        /**< number of binary variables with non-zero objective coefficient */
   SCIP_VAR**            contributors,       /**< array to store the contriboters */
   SCIP_HASHTABLE*       uselesscliques,     /**< hash table to store useless cliques, or NULL */
   int*                  ncontributors,      /**< pointer to store number of contributor to the objective contribution */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   assert(SCIPvarIsBinary(var));
   assert(contributors != NULL);
   assert(ncontributors != NULL);

   /* collects the contribution of the variable itself w.r.t. the best bound */
   (*objchg) = getVarObjchg(var, SCIPvarGetBestBoundType(var), bound);

   (*ncontributors) = 0;

   /* add the objective contribution from the implication variable */
   SCIP_CALL( collectMinactImplicVars(scip, var, bound, binobjvarmap, collectedvars, nbinobjvars, contributors, uselesscliques, ncontributors, objchg) );

   return SCIP_OKAY;
}

/** computes for the given binary variable the objective contribution by fixing it to given bound w.r.t. minimum
 *  activity of the objective function; this can be done w.r.t. global variable bounds (local == FALSE), w.r.t. local
 *  variable bounds (local == TRUE && bdchgidx == NULL), and w.r.t. given time stamp (local == TRUE && bdchgidx != NULL)
 */
static
SCIP_RETCODE getMinactObjchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_OBJIMPLICS*      objimplics,         /**< objective implication data for the given variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, or NULL */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_Bool             local,              /**< propagate local bounds, otherwise global bounds */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   assert(SCIPvarIsBinary(var));

   /* collects the contribution of the variable itself w.r.t. the best bound */
   (*objchg) = getVarObjchg(var, SCIPvarGetBestBoundType(var), bound);

   /* add the objective contribution from the implication variable */
   SCIP_CALL( getMinactImplicObjchg(scip, var, objimplics, bdchgidx, bound, local, objchg) );

   return SCIP_OKAY;
}

/** returns the global (that means w.r.t. global bounds of the variables) objective change provided by all cliques of
 *  the given variable by fixing it to the given bound w.r.t. maximum activity of the objective function
 *
 *  Let I(0) and I(1) be all implications of the given variable which follow by fixing it to given bound and evaluate to
 *  fixing the implication variable to zero (I(0)) or one (I(1)), respectively. The objective change provided by these
 *  implications are:
 *
 *  \f[
 *  \displaystyle
 *  sum_{x\in I(1)} (1 - \mbox{worstbound}(x)) \cdot \mbox{objval}(x) - sum_{x\in I(1)} \mbox{worst}(x) \cdot \mbox{objval}(x)
 *  =
 *  sum_{x\in I(0) \cup I(1)} (\mbox{impliedbound}(x) - \mbox{worstbound}(x)) \cdot \mbox{objval}(x)
 *  \f]
 */
static
SCIP_RETCODE getMaxactImplicObjchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   SCIP_Bool varfixing;
   int ncliques;
   int nvars;

   assert(scip != NULL);
   assert(SCIPvarIsBinary(var));
   assert(SCIPvarGetLbGlobal(var) < 0.5);
   assert(SCIPvarGetUbGlobal(var) > 0.5);
   assert(bound == SCIP_BOUNDTYPE_LOWER || bound == SCIP_BOUNDTYPE_UPPER);
   assert(objchg != NULL);

   varfixing = (SCIP_Bool)bound;
   assert((SCIP_Bool)SCIP_BOUNDTYPE_LOWER == FALSE);
   assert((SCIP_Bool)SCIP_BOUNDTYPE_UPPER == TRUE);

   *objchg = 0.0;
   ncliques = SCIPvarGetNCliques(var, varfixing);

   if( ncliques > 0 )
   {
      SCIP_CLIQUE** cliques;
      SCIP_CLIQUE* clique;
      SCIP_VAR** clqvars;
      SCIP_VAR** probvars;
      SCIP_VAR* clqvar;
      SCIP_Bool* clqvalues;
      int* entries;
      int* ids;
      SCIP_Real obj;
      int nclqvars;
      int nentries;
      int objmult;
      int nids;
      int id;
      int c;
      int v;

      assert(SCIPisTransformed(scip));

      nentries = SCIPgetNVars(scip) - SCIPgetNContVars(scip) + 1;

      SCIP_CALL( SCIPallocBufferArray(scip, &ids, 2*nentries) );
      nids = 0;
      /* @todo move this memory allocation to SCIP_SET and add a memory list there, to decrease the number of
       *       allocations and clear ups
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &entries, nentries) );
      BMSclearMemoryArray(entries, nentries);

      cliques = SCIPvarGetCliques(var, varfixing);
      assert(cliques != NULL);

      /* iterate over all cliques and determine all importantimplications */
      for( c = ncliques - 1; c >= 0; --c )
      {
         clique = cliques[c];
         clqvars = SCIPcliqueGetVars(clique);
         clqvalues = SCIPcliqueGetValues(clique);
         nclqvars = SCIPcliqueGetNVars(clique);
         assert(nclqvars > 0);
         assert(clqvars != NULL);
         assert(clqvalues != NULL);

         if( nclqvars > MAX_CLIQUELENGTH )
            continue;

         /* iterate over all clique variables */
         for( v = nclqvars - 1; v >= 0; --v )
         {
            clqvar = clqvars[v];
            assert(clqvar != NULL);

            objmult = (int)!clqvalues[v] - (int)SCIPvarGetWorstBoundType(clqvar);
            assert(-1 <= objmult && objmult <= 1);

            /* ignore binary variable which are either fixed and were the objective contribution will not be zero */
            if( clqvar != var && objmult != 0 && SCIPvarIsActive(clqvar) &&
               (SCIPvarGetLbGlobal(clqvar) < 0.5 && SCIPvarGetUbGlobal(clqvar) > 0.5) && !SCIPisZero(scip, SCIPvarGetObj(clqvar)) )
            {
               int probindex = SCIPvarGetProbindex(clqvar) + 1;
               assert(0 < probindex && probindex < nentries);

               /* check that the variable was not yet visited  */
               assert(entries[probindex] == 0 || entries[probindex] == objmult);
               if( entries[probindex] == 0 )
               {
                  /* memorize probindex */
                  ids[nids] = probindex;
                  ++nids;

                  assert(ABS(objmult) == 1);

                  /* mark variable as visited */
                  entries[probindex] = objmult;
               }
            }
         }
      }

      probvars = SCIPgetVars(scip);
      assert(probvars != NULL);

      /* add all implied objective values */
      for( v = nids - 1; v >= 0; --v )
      {
         id = ids[v];
         assert(0 < id && id < nentries);
         assert(entries[id] != 0);

         clqvar = probvars[id - 1];
         assert(clqvar != NULL);
         assert(SCIPvarIsBinary(clqvar));
         assert(SCIPvarIsActive(clqvar));
         assert(SCIPvarGetLbGlobal(clqvar) < 0.5);
         assert(SCIPvarGetUbGlobal(clqvar) > 0.5);

         obj = SCIPvarGetObj(clqvar);
         assert(!SCIPisZero(scip, obj));

         *objchg += entries[id] * obj;
      }

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &entries);
      SCIPfreeBufferArray(scip, &ids);
   }

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(scip, "objective contribution when variable <%s> fixed to %u using cliques is %g\n", SCIPvarGetName(var),
      varfixing, *objchg);
#endif

   /* collect non-binary implication information */
   nvars = SCIPvarGetNImpls(var, varfixing);

   if( nvars > 0 )
   {
      SCIP_VAR** vars;
      SCIP_VAR* implvar;
      SCIP_Real* bounds;
      SCIP_BOUNDTYPE* boundtypes;
      SCIP_Real obj;
      SCIP_Real lb;
      SCIP_Real ub;
      int v;

      vars =  SCIPvarGetImplVars(var, varfixing);
      boundtypes = SCIPvarGetImplTypes(var, varfixing);
      bounds = SCIPvarGetImplBounds(var, varfixing);

      for( v = nvars - 1; v >= 0; --v )
      {
         implvar = vars[v];
         assert(implvar != NULL);

         lb = SCIPvarGetLbLocal(implvar);
         ub = SCIPvarGetUbLocal(implvar);
         obj = SCIPvarGetObj(implvar);

         /* ignore binary variable which are fixed or not of column status */
         if( SCIPisZero(scip, obj) )
            continue;

         /* add up objective change if applicable */
         if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER && SCIPvarGetWorstBoundType(implvar) == SCIP_BOUNDTYPE_LOWER && SCIPisFeasGT(scip, bounds[v], lb) )
            *objchg += (bounds[v] - lb)*obj;
         else if( boundtypes[v] == SCIP_BOUNDTYPE_UPPER && SCIPvarGetWorstBoundType(implvar) == SCIP_BOUNDTYPE_UPPER && SCIPisFeasLT(scip, bounds[v], ub) )
            *objchg += (bounds[v] - ub)*obj;
      }
   }

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(scip, "objective contribution when variable <%s> fixed to %u using cliques and implications is %g\n", SCIPvarGetName(var),
      varfixing, *objchg);
#endif

   return SCIP_OKAY;
}

/** computes for the given binary variable the gloabl (that means w.r.t. global bounds of the variables) objective
 *  contribution by fixing it to given bound w.r.t. maximum activity of the objective function
 */
static
SCIP_RETCODE getMaxactObjchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_Bool             useimplics,         /**< should implications be used */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   assert(scip != NULL);
   assert(SCIPvarIsBinary(var));
   assert(objchg != NULL);

   *objchg = 0;

   /* check if the implications should be used to increase the objective contribution for given variable */
   if( useimplics )
   {
      /* using cliques and @todo other implications */
      SCIP_CALL( getMaxactImplicObjchg(scip, var, bound, objchg) );
   }

   /* collects the contribution of the variable itself w.r.t. the worst bound */
   *objchg += getVarObjchg(var, SCIPvarGetWorstBoundType(var), bound);

   return SCIP_OKAY;
}

/** reset variables array which marks variables which are collected */
static
void resetContributors(
   SCIP_HASHMAP*         binobjvarmap,       /**< hash map mapping binary variables with none-zero objective to position in collected variables arrays */
   SCIP_Bool*            collectedvars,      /**< temporary buffer to mark collected variables which should be reset */
   SCIP_VAR**            contributors,       /**< temporary buffer to use for collecting contributors */
   int                   ncontributors       /**< number of contributors */
   )
{
   SCIP_VAR* var;
   int pos;
   int c;

   for( c = 0; c < ncontributors; ++c )
   {
      var = contributors[c];
      assert(var != NULL);

      assert(SCIPhashmapExists(binobjvarmap, var));
      pos = (int)(size_t)SCIPhashmapGetImage(binobjvarmap, var);
      assert(pos > 0);
      collectedvars[pos] = FALSE;
   }
}

/** check if the given variable should be collected for the minimum activity propagation */
static
SCIP_RETCODE collectMinactVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_OBJIMPLICS**     objimplics,         /**< pointer to store the objective implication data structure w.r.t. minimum activity */
   SCIP_Bool             useimplics,         /**< should implications be used */
   SCIP_HASHMAP*         binobjvarmap,       /**< hash map mapping binary variables with none-zero objective to position in collected variables arrays */
   SCIP_Bool*            collectedlbvars,    /**< temporary buffer to mark collected variables for lower bound fixing */
   SCIP_Bool*            collectedubvars,    /**< temporary buffer to mark collected variables for upper bound fixing */
   int                   nbinobjvars,        /**< number of binary variables with non-zero objective coefficient */
   SCIP_VAR**            contributors,       /**< temporary buffer to use for collecting contributors */
   SCIP_HASHTABLE*       uselesscliques,     /**< hash table to store useless cliques, or NULL */
   SCIP_Bool*            collect             /**< pointer to store if the variable should be stored */
   )
{
   SCIP_Real lbobjchg;
   SCIP_Real ubobjchg;
   SCIP_Real objval;
   int nlbcliques;
   int nubcliques;

   assert(objimplics != NULL);

   objval = SCIPvarGetObj(var);
   (*objimplics) = NULL;

   if( SCIPisZero(scip, objval) )
      (*collect) = FALSE;
   else
      (*collect) = TRUE;

   nlbcliques = SCIPvarGetNCliques(var,  FALSE);
   nubcliques = SCIPvarGetNCliques(var,  TRUE);

   /* check if implications should be used and if implications are existing */
   if( useimplics && nlbcliques + nubcliques > 0 )
   {
      int nlbcontributors;
      int nubcontributors;

      assert((SCIP_Bool)SCIP_BOUNDTYPE_LOWER == FALSE);
      assert((SCIP_Bool)SCIP_BOUNDTYPE_UPPER == TRUE);

      /* get contribution of variable by fixing it to its lower bound w.r.t. minimum activity of the objective function */
      SCIP_CALL( collectMinactObjchg(scip, var, SCIP_BOUNDTYPE_LOWER, binobjvarmap, collectedlbvars, nbinobjvars,
            contributors, uselesscliques, &nlbcontributors, &lbobjchg) );
      assert(!SCIPisNegative(scip, lbobjchg));

      SCIPdebugMsg(scip, "variable <%s> fixed to bound=%d implies %d(%d)\n", SCIPvarGetName(var),
         SCIP_BOUNDTYPE_LOWER, 0, nlbcontributors);

      /* ignore implications if the variable has a zero objective coefficient and implications only one variable, since
       * this is covered by that implied variable
       */
      if( !(*collect) && nlbcontributors == 1 )
      {
         /* reset lower bound contributors */
         resetContributors(binobjvarmap, collectedlbvars, contributors, nlbcontributors);

         assert(SCIPisZero(scip, objval));
         nlbcontributors = 0;
      }

      /* get contribution of variable by fixing it to its upper bound w.r.t. minimum activity of the objective function */
      SCIP_CALL( collectMinactObjchg(scip, var, SCIP_BOUNDTYPE_UPPER, binobjvarmap, collectedubvars, nbinobjvars,
            &contributors[nlbcontributors], uselesscliques, &nubcontributors, &ubobjchg) );
      assert(!SCIPisNegative(scip, ubobjchg));

      SCIPdebugMsg(scip, "variable <%s> fixed to bound=%d implies %d(%d)\n", SCIPvarGetName(var),
         SCIP_BOUNDTYPE_UPPER, 0, nubcontributors);

      /* ignore implications if the variable has a zero objective coefficient and implications only one variable, since
       * this is covered by that implied variable
       */
      if( !(*collect) && nubcontributors == 1 )
      {
         /* reset upper bound contributors */
         resetContributors(binobjvarmap, collectedubvars, &contributors[nlbcontributors], nubcontributors);

         assert(SCIPisZero(scip, objval));
         nubcontributors = 0;
      }

      if( (*collect) || nlbcontributors > 1 || nubcontributors > 1 )
      {
         /* creates an objective implication data structure, fixes (globally) variables which are implied by lower and upper
          * bound fixing, and clears the collected arrays for lower and upper bound
          */
         SCIP_CALL( objimplicsCreate(scip, objimplics, contributors, binobjvarmap, collectedlbvars, collectedubvars, lbobjchg, ubobjchg, nlbcontributors, nubcontributors) );
         (*collect) = TRUE;
      }
      else
      {
         /* reset lower bound contributors */
         resetContributors(binobjvarmap, collectedlbvars, contributors, nlbcontributors);

         /* reset upper bound contributors */
         resetContributors(binobjvarmap, collectedubvars, &contributors[nlbcontributors], nubcontributors);
      }
   }
   else if( (*collect) )
   {
      lbobjchg = getVarObjchg(var, SCIPvarGetBestBoundType(var), SCIP_BOUNDTYPE_LOWER);
      ubobjchg = getVarObjchg(var, SCIPvarGetBestBoundType(var), SCIP_BOUNDTYPE_UPPER);
      assert(!SCIPisZero(scip, lbobjchg) || !SCIPisZero(scip, ubobjchg));
      assert(!SCIPisNegative(scip, lbobjchg));
      assert(!SCIPisNegative(scip, ubobjchg));

      /* creates an "empty" objective implication data structure */
      SCIP_CALL( objimplicsCreate(scip, objimplics, NULL, NULL, NULL, NULL, lbobjchg, ubobjchg, 0, 0) );
   }

   return SCIP_OKAY;
}

/** check if the given variable should be collected for the maximum activity propagation */
static
SCIP_RETCODE collectMaxactVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_Bool             useimplics,         /**< should implications be used */
   SCIP_Real*            objchg,             /**< pointer to store the objective change w.r.t. maximum activity */
   SCIP_Bool*            isnotzero           /**< pointer to store if the objective change is unequal to zero or not */
   )
{
   SCIP_Real lbobjchg;
   SCIP_Real ubobjchg;

   assert(scip != NULL);
   assert(SCIPvarIsBinary(var));
   assert(objchg != NULL);
   assert(isnotzero != NULL);

   /* get contribution of variable by fixing it to its lower bound w.r.t. maximum activity of the objective function */
   SCIP_CALL( getMaxactObjchg(scip, var, SCIP_BOUNDTYPE_LOWER, useimplics, &lbobjchg) );
   assert(!SCIPisPositive(scip, lbobjchg));

   /* get contribution of variable by fixing it to its upper bound w.r.t. maximum activity of the objective function */
   SCIP_CALL( getMaxactObjchg(scip, var, SCIP_BOUNDTYPE_UPPER, useimplics, &ubobjchg) );
   assert(!SCIPisPositive(scip, ubobjchg));

   (*objchg) = MIN(lbobjchg, ubobjchg);

   /* only consider variables with non-zero objective contribution */
   if( SCIPisZero(scip, (*objchg)) )
      *isnotzero = FALSE;
   else
      *isnotzero = TRUE;

   return SCIP_OKAY;
}

/** initializate the propagator */
static
SCIP_RETCODE propdataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_HASHMAP* binobjvarmap;
   int nvars;
   int nbinvars;
   int nintvars;
   int nminactvars;
   int nmaxactvars;
   int nobjintvars;
   int nobjcontvars;
   int nobjvars;
   int nbinobjvars;
   int v;

   assert(scip != NULL);
   assert(propdata != NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nintvars = nvars - SCIPgetNContVars(scip);

   nbinvars = 0;
   nobjvars = 0;
   nbinobjvars = 0;

   SCIP_CALL( SCIPhashmapCreate(&binobjvarmap, SCIPblkmem(scip), SCIPgetNObjVars(scip)) );

   /* count and collect variable problem indices of variables with non-zero objective coefficient */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);

      if( !SCIPisZero(scip, SCIPvarGetObj(var)) )
      {
         nobjvars++;

         if( SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPhashmapInsert(binobjvarmap, (void*)var, (void*)(size_t)(nbinobjvars + 1)) );
            nbinobjvars++;
         }
      }

      if( SCIPvarIsBinary(var) )
         nbinvars++;
   }

   nminactvars = 0;
   nmaxactvars = 0;
   nobjintvars = 0;
   nobjcontvars = 0;

   if( nobjvars > 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;
      SCIP_OBJIMPLICS* objimplics;
      SCIP_HASHTABLE* uselesscliques;
      SCIP_VAR** contributors;
      SCIP_Bool* collectedlbvars;
      SCIP_Bool* collectedubvars;
      SCIP_Bool collect;
      SCIP_Bool useimplics;
      SCIP_Real objval;
      SCIP_Real objchg;

      eventhdlr = propdata->eventhdlr;
      assert(eventhdlr != NULL);

      useimplics = (propdata->propuseimplics && nbinobjvars < propdata->maximplvars);

      /* allocate memory for all arrays */
      propdata->minactsize = nbinvars;
      propdata->maxactsize = nbinvars;
      propdata->objintvarssize = nobjvars - nbinobjvars;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->minactvars, propdata->minactsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->minactimpls, propdata->minactsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->maxactvars, propdata->maxactsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->maxactchgs, propdata->maxactsize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->objintvars, propdata->objintvarssize) );

      if( useimplics )
      {
         int ncliques;

         /* create temporary buffer */
         /* we store both lb and ub contributors in array contributors, and both could be nbinobjvars, we need twice that size */
         SCIP_CALL( SCIPallocBufferArray(scip, &contributors, 2 * nbinobjvars) );
         /* @todo: use SCIPallocCleanBufferArray instead? */
         SCIP_CALL( SCIPallocClearBufferArray(scip, &collectedlbvars, nbinobjvars+1) );
         /* @todo: use SCIPallocCleanBufferArray instead? */
         SCIP_CALL( SCIPallocClearBufferArray(scip, &collectedubvars, nbinobjvars+1) );

         ncliques = SCIPgetNCliques(scip);

         if( ncliques > 0 )
         {
            SCIP_CALL( SCIPhashtableCreate(&uselesscliques, SCIPblkmem(scip), ncliques,
                  cliqueGetHashkey, cliqueIsHashkeyEq, cliqueGetHashkeyVal, NULL) );
         }
         else
            uselesscliques = NULL;
      }
      else
      {
         contributors = NULL;
         collectedlbvars = NULL;
         collectedubvars = NULL;
         uselesscliques = NULL;
      }

      /* collect the variables with non-zero objective contribution and catch global bound tighten events that decrease the
       * maximal pseudo objective activity
       */
      for( v = 0; v < nvars && (nobjintvars == 0 || nobjintvars < propdata->objintvarssize); ++v )
      {
         var = vars[v];
         assert(var != NULL);

         objval = SCIPvarGetObj(var);

         if( SCIPvarIsBinary(var) )
         {
            /* ignore variables which are globally fixed */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
            {
#ifndef NDEBUG
               /* check that the binary implications are applied for binary variables which are globally fixed */
               checkImplicsApplied(scip, var);
#endif
               continue;
            }

            /* check if the variable should be collected for the minimum activity propagation */
            SCIP_CALL( collectMinactVar(scip, var, &objimplics, useimplics, binobjvarmap, collectedlbvars, collectedubvars,
                  nbinobjvars, contributors, uselesscliques, &collect) );

            if( collect )
            {
               assert(nminactvars < nbinvars);
               assert(objimplics != NULL);
               assert(objimplics->nlbimpls + objimplics->nubimpls <= nbinobjvars);

               /* collect the binary variable with non-zero objective contribution */
               propdata->minactvars[nminactvars] = var;
               propdata->minactimpls[nminactvars] = objimplics;
               nminactvars++;

               /* catch bound relax event for the binary variable to handel the firstnonfixed index correctly */
               SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDRELAXED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );

               SCIPdebugMsg(scip, "variable <%s>[obj: <%g>] implicit objective change %g\n",
                  SCIPvarGetName(var), objval, objimplics->maxobjchg);

               /* captures the variable */
               SCIP_CALL( SCIPcaptureVar(scip, var) ) ;
            }
            /* check if the variable should be collected for the maximum activity propagation */
            SCIP_CALL( collectMaxactVar(scip, var, useimplics, &objchg, &collect) );

            if( collect )
            {
               assert(nmaxactvars < nbinvars);

               /* collect the binary variable with non-zero objective contribution */
               propdata->maxactvars[nmaxactvars] = var;
               propdata->maxactchgs[nmaxactvars] = -objchg;
               nmaxactvars++;

               /* catch bound change events if the variable has a non-zero objective coefficient to check if the maximum
                * activity of the objective function changed
                */
               SCIP_CALL( catchObjEvent(scip, propdata, eventhdlr, var) );

               /* captures the variable */
               SCIP_CALL( SCIPcaptureVar(scip, var) ) ;
            }
         }
         else
         {
            /* only consider non-binary variables with a non-zero objective */
            if( SCIPisZero(scip, objval) )
               continue;

            assert(nobjintvars < propdata->objintvarssize);

            propdata->objintvars[nobjintvars] = var;
            nobjintvars++;

            if( v >= nintvars )
               nobjcontvars++;

            /* catch bound change events if the variable has a non-zero objective coefficient to check if the maximum
             * activity of the objective function changed
             */
            SCIP_CALL( catchObjEvent(scip, propdata, eventhdlr, var) );

            /* captures the variable */
            SCIP_CALL( SCIPcaptureVar(scip, var) );
         }
      }

      if( useimplics )
      {
         if( uselesscliques != NULL )
            SCIPhashtableFree(&uselesscliques);

         SCIPfreeBufferArray(scip, &collectedubvars);
         SCIPfreeBufferArray(scip, &collectedlbvars);
         SCIPfreeBufferArray(scip, &contributors);
      }

      if( nminactvars == 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &propdata->minactvars, propdata->minactsize);
         SCIPfreeBlockMemoryArray(scip, &propdata->minactimpls, propdata->minactsize);
         propdata->minactsize = 0;
         propdata->minactvars = NULL;
         propdata->minactimpls = NULL;
      }
      else
      {
         /* sort binary variables with respect to the absolute value of their maximal potential objective contribution for
          * the minimum activity of the objective function
          */
         SCIPsortDownPtrPtr((void**)propdata->minactimpls, (void**)propdata->minactvars, objimplicsComp, nminactvars);

         SCIPdebugMsg(scip, "%d binary variables with non-zero objective contribution w.r.t. the minimum activity of the objective function\n", nminactvars);
      }

      if( nmaxactvars == 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &propdata->maxactvars, propdata->maxactsize);
         SCIPfreeBlockMemoryArray(scip, &propdata->maxactchgs, propdata->maxactsize);
         propdata->maxactsize = 0;
         propdata->maxactvars = NULL;
         propdata->maxactchgs = NULL;
      }
      else
      {
         /* sort binary variables with respect to the absolute value of their maximal potential objective contribution for
          * the maximum activity of the objective function
          */
         SCIPsortDownRealPtr(propdata->maxactchgs, (void**)propdata->maxactvars, nmaxactvars);

         SCIPdebugMsg(scip, "%d binary variables with non-zero objective contribution w.r.t. the maximum activity of the objective function\n", nmaxactvars);
      }

      if( nobjintvars == 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &propdata->objintvars, propdata->objintvarssize);
         propdata->objintvarssize = 0;
         propdata->objintvars = NULL;
      }
      else
      {
         /* sort integer variables with respect to the absolute value of their objective coefficient */
         SCIPsortDownPtr((void**)propdata->objintvars, varCompObj, nobjintvars - nobjcontvars);

         /* sort continuous variables with respect to the absolute value of their objective coefficient */
         SCIPsortDownPtr((void**)(&propdata->objintvars[nobjintvars - nobjcontvars]), varCompObj, nobjcontvars);

         SCIPdebugMsg(scip, "%d integer variables and %d continuous variables with non-zero objective contribution\n",
            nobjintvars - nobjcontvars, nobjcontvars);
      }
   }

   SCIPhashmapFree(&binobjvarmap);

   propdata->nminactvars = nminactvars;
   propdata->nmaxactvars = nmaxactvars;
   propdata->nobjintvars = nobjintvars;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;
   propdata->lastvarnum = -1;
   propdata->glbfirstnonfixed = 0;
   propdata->maxactfirstnonfixed = 0;
   propdata->firstnonfixed = 0;
   propdata->nnewvars = 0;
   propdata->cutoffbound = SCIPinfinity(scip);
   propdata->lastlowerbound = -SCIPinfinity(scip);
   propdata->glbpseudoobjval = -SCIPinfinity(scip);

   propdata->initialized = TRUE;

   /* due to scaling after presolving we need to update the global pseudoactivity and the cutoffbound */
   propdata->glbpropagated = FALSE;
   propdata->glbpseudoobjval = SCIPgetGlobalPseudoObjval(scip);
   propdata->cutoffbound = SCIPgetCutoffbound(scip);
   assert(SCIPgetDepth(scip) > 0 || SCIPisFeasEQ(scip, propdata->glbpseudoobjval, SCIPgetPseudoObjval(scip)));

   /* create hash table which is used for resolving bound changes */
   if( nminactvars > 0 )
   {
      SCIP_CALL( SCIPhashtableCreate(&propdata->addedvars, SCIPblkmem(scip), nvars,
            SCIPvarGetHashkey, SCIPvarIsHashkeyEq, SCIPvarGetHashkeyVal, NULL) );
   }
   else
      propdata->addedvars = NULL;


   return SCIP_OKAY;
}

/** adds for the given non-binary variable a conflict bound depending on its objective contribution */
static
SCIP_RETCODE addConflictBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check for objective contribution */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Real*            reqpseudoobjval     /**< pointer to store the remaining minimum activity which has to be proven */
   )
{
   SCIP_Real objval;

   objval = SCIPvarGetObj(var);
   assert(!SCIPisZero(scip, objval));

   if( objval > 0.0 )
   {
      SCIP_Real loclb;
      SCIP_Real glblb;

      glblb = SCIPvarGetLbGlobal(var);
      loclb = SCIPgetVarLbAtIndex(scip, var, bdchgidx, FALSE);
      assert(SCIPisFeasGE(scip, loclb, glblb));

      /* check if the local lower bound (at time stamp bdchgidx) is larger than the global lower bound */
      if( SCIPisGT(scip, loclb, glblb) )
      {
         SCIPdebugMsg(scip, "  add bound change <%s>[%g] >= <%g>\n", SCIPvarGetName(var), objval, loclb);
         SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );

         /* hard comparison  is enough to make requiredpseudoobjval nonincreasing */
         assert((loclb - glblb) * objval > 0.0);

         (*reqpseudoobjval) -= (loclb - glblb) * objval;
      }
   }
   else
   {
      SCIP_Real locub;
      SCIP_Real glbub;

      glbub = SCIPvarGetUbGlobal(var);
      locub = SCIPgetVarUbAtIndex(scip, var, bdchgidx, FALSE);
      assert(SCIPisFeasLE(scip, locub, glbub));

      /* check if the local upper bound (at time stamp bdchgidx) is smaller than the global upper bound */
      if( SCIPisLT(scip, locub, glbub) )
      {
         SCIPdebugMsg(scip, "  add bound change <%s>[%g] <= <%g>\n", SCIPvarGetName(var), objval, locub);
         SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );

         /* hard comparison  is enough to make requiredpseudoobjval nonincreasing */
         assert((locub - glbub) * objval > 0.0);

         (*reqpseudoobjval) -= (locub - glbub) * objval;
      }
   }

   return SCIP_OKAY;
}

/** check for the given implication variables if they also contribute to the required minimum activity */
static
SCIP_RETCODE getConflictImplics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable to check for objective contribution */
   int                   start,              /**< start index */
   int                   end,                /**< end index */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_HASHTABLE*       addedvars,          /**< hash table containing variables which are already added directly or implicitly due to implications */
   SCIP_Real*            reqpseudoobjval,    /**< pointer to store the remaining minimum activity which has to be proven */
   SCIP_Bool*            foundimplics        /**< pointer to store if an implication is found */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   int v;

   assert(foundimplics != NULL);
   assert(*foundimplics == FALSE);

   for( v = start; v < end; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarIsBinary(var));

      /* we need to take the bounds after the bdchgidx here, since the variable of the bound change may be the implied one;
       * we already counted its contribution before, so we want to see it as fixed here, which it is after the bound change.
       */
      lb = SCIPgetVarLbAtIndex(scip, var, bdchgidx, TRUE);
      ub = SCIPgetVarUbAtIndex(scip, var, bdchgidx, TRUE);

      if( lb < 0.5 && ub > 0.5 && !SCIPhashtableExists(addedvars, (void*)var) )
      {
         (*reqpseudoobjval) -= REALABS(SCIPvarGetObj(var));
         SCIPdebugMsg(scip, "  implicated variables <%s>[%g] bdchgidx [%g,%g] -> remaining <%g>\n", SCIPvarGetName(var), SCIPvarGetObj(var), lb, ub, *reqpseudoobjval);

         SCIP_CALL( SCIPhashtableInsert(addedvars, (void*)var) );
         assert(SCIPhashtableExists(addedvars, (void*)var));
         (*foundimplics) = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** adds for the given binary variable a conflict bound depending on its objective contribution */
static
SCIP_RETCODE addConflictBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check for objective contribution */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_OBJIMPLICS*      objimplics,         /**< objective implication data for the given variable */
   SCIP_HASHTABLE*       addedvars,          /**< hash table containing variables which are already add directly or implicitly due to implications */
   SCIP_Bool             respropuseimplics,  /**< should implications be used */
   SCIP_Real*            reqpseudoobjval     /**< pointer to store the remaining minimum activity which has to be proven */
   )
{
   SCIP_Real objval;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool foundimplics;

   assert(SCIPvarIsBinary(var));

   if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
      return SCIP_OKAY;

   lb = SCIPgetVarLbAtIndex(scip, var, bdchgidx, FALSE);
   ub = SCIPgetVarUbAtIndex(scip, var, bdchgidx, FALSE);

   objval = SCIPvarGetObj(var);
   foundimplics = FALSE;

   /* only consider variables which are fixed */
   if( lb > 0.5 )
   {
      if( respropuseimplics )
      {
         SCIP_CALL( getConflictImplics(scip, objimplics->objvars, objimplics->nlbimpls, objimplics->nlbimpls + objimplics->nubimpls,
               bdchgidx, addedvars, reqpseudoobjval, &foundimplics) );
      }

      /* check if the binary variable has a positive contribution (positive objective coefficient since it is fixed to
       * one) or is needed due a positive contribution of an implied variable
       */
      if( foundimplics || SCIPisPositive(scip, objval) )
      {
         SCIPdebugMsg(scip, "  add bound change <%s>[%g] >= <%g> bdchgidx [%g,%g]\n", SCIPvarGetName(var), objval, lb, lb, ub);
         SCIP_CALL( SCIPaddConflictLb(scip, var, NULL) );

         (*reqpseudoobjval) -= MAX(0.0, objval);

         if( addedvars != NULL )
         {
            assert(!SCIPhashtableExists(addedvars, (void*)var));
            SCIP_CALL( SCIPhashtableInsert(addedvars, (void*)var) );
         }
      }
   }
   else if( ub < 0.5 )
   {
      if( respropuseimplics )
      {
         SCIP_CALL( getConflictImplics(scip, objimplics->objvars, 0, objimplics->nlbimpls,
               bdchgidx, addedvars, reqpseudoobjval, &foundimplics) );
      }

      /* check if the binary variable has a positive contribution (negative objective coefficient since it is fixed to
       * zero) or is needed due a positive contribution of an implied variable
       */
      if( foundimplics || SCIPisNegative(scip, objval) )
      {
         SCIPdebugMsg(scip, "  add bound change <%s>[%g] <= <%g> bdchgidx=[%g,%g]\n", SCIPvarGetName(var), objval, ub, lb, ub);
         SCIP_CALL( SCIPaddConflictUb(scip, var, NULL) );

         (*reqpseudoobjval) +=  MIN(0.0, objval);

         if( addedvars != NULL )
         {
            assert(!SCIPhashtableExists(addedvars, (void*)var));
            SCIP_CALL( SCIPhashtableInsert(addedvars, (void*)var) );
         }
      }
   }

   return SCIP_OKAY;
}


/** resolves a propagation by supplying the variables whose bound changes increased the pseudo objective value above the
 *  cutoff bound
 */
static
SCIP_RETCODE adjustCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable that was deduced */
   int                   inferinfo,          /**< inference information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_HASHTABLE*       addedvars,          /**< hash table which contains variables which are already added or implicitly given as reason for the resolve, or NULL */
   SCIP_Real*            cutoffbound         /**< pointer to store the adjusted cutoff bound */
   )
{
   if( inferinfo != -1 )
   {
      SCIP_OBJIMPLICS* objimplics;
      SCIP_Bool foundimplics;
      int start;
      int end;

      assert(var != NULL);
      assert(SCIPvarIsBinary(var));
      assert(bdchgidx != NULL);
      assert(SCIPisEQ(scip, SCIPgetVarLbAtIndex(scip, var, bdchgidx, TRUE), SCIPgetVarUbAtIndex(scip, var, bdchgidx, TRUE)));
      assert(inferinfo >= 0);
      assert(inferinfo < propdata->nminactvars);
      assert((SCIP_Bool)SCIP_BOUNDTYPE_LOWER == FALSE);
      assert((SCIP_Bool)SCIP_BOUNDTYPE_UPPER == TRUE);

      objimplics = propdata->minactimpls[inferinfo];
      assert(objimplics != NULL);

      /* get the objective contribution if we would fix the binary inference variable to its other bound */
      (*cutoffbound) -= getVarObjchg(var, SCIPvarGetBestBoundType(var), boundtype);
      foundimplics = FALSE;

      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         start = 0;
         end = objimplics->nlbimpls;
      }
      else
      {
         start = objimplics->nlbimpls;
         end = objimplics->nlbimpls + objimplics->nubimpls;
      }

      if( addedvars != NULL )
      {
         SCIP_CALL( getConflictImplics(scip, objimplics->objvars, start, end, bdchgidx, addedvars, cutoffbound, &foundimplics) );
      }
   }
   else
   {
      SCIP_Real glbbound;
      SCIP_Real newbound;
      SCIP_Real objval;

      objval = SCIPvarGetObj(var);

      assert(!SCIPisZero(scip, objval));

      if( objval > 0.0 )
      {
         newbound = SCIPgetVarUbAtIndex(scip, var, bdchgidx, TRUE);
         glbbound = SCIPvarGetLbGlobal(var);
      }
      else
      {
         newbound = SCIPgetVarLbAtIndex(scip, var, bdchgidx, TRUE);
         glbbound = SCIPvarGetUbGlobal(var);
      }

      /* in case the variable is integral we just need to prove the newbound plus/minus (1 - epsilon) since the this bound
       * would be rounded to newbound due to integrability which is global information
       */
      if( SCIPvarIsIntegral(var) )
      {
         if( objval > 0.0 )
            newbound += 1 - 10 * SCIPfeastol(scip);
         else
            newbound -= 1 - 10 * SCIPfeastol(scip);
      }

      /* adjust the cutoff bound by the portion the inference variable contributes to the presudo objective activity
       * (minactivity)
       */
      assert(!SCIPisNegative(scip, objval * (newbound - glbbound)));
      (*cutoffbound) -= objval * (newbound - glbbound);
   }

   return SCIP_OKAY;
}


/** resolves a propagation by supplying the variables whose bound changes increased the pseudo objective value above the
 *  cutoff bound
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Real             cutoffbound,        /**< the global cutoff */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL for conflict analysis initialization */
   int                   inferinfo,          /**< inference information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_HASHTABLE* addedvars;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real glbpseudoobjval;
   SCIP_Real reqpseudoobjval;
   SCIP_Bool infinity;
   int nvars;
   int v;

   infinity = FALSE;
   addedvars = NULL;
   nvars = propdata->nminactvars;

   /* the global pseudo value gives us a global valid minimal activity
    *
    * @note The global pseudo objective activity can be minus infinity. In that case all variable are part of the
    *       reason/explanation
    *
    * @note If the global pseudo objective activity is greater than the required minactivity, the local bound change
    *       which has to explained is actually (now) a global one. That means, the reason/explanation is empty
    */
   glbpseudoobjval = SCIPgetGlobalPseudoObjval(scip);

   if( SCIPisInfinity(scip, -glbpseudoobjval) )
   {
      infinity = TRUE;
      reqpseudoobjval = cutoffbound;
   }
   else
   {
      /* clear hash table for storing variables which are not needed to add the reason due to global implications or
       * already added
       */
      if( nvars > 0 )
      {
         addedvars = propdata->addedvars;
         SCIPhashtableRemoveAll(addedvars);
      }

      if( infervar != NULL )
      {
         SCIP_CALL( adjustCutoffbound(scip, propdata, infervar, inferinfo, boundtype, bdchgidx, addedvars, &cutoffbound) );
      }

      reqpseudoobjval = cutoffbound - glbpseudoobjval;
   }

   SCIPdebugMsg(scip, "resolve propagation global pseudo objective <%g>, cutoff bounda <%g>, required minactivity <%g>\n",
      glbpseudoobjval, cutoffbound, reqpseudoobjval);

   /* the variables responsible for the propagation are the ones with
    *  - obj > 0 and local lb > global lb
    *  - obj < 0 and local ub < global ub
    *
    * collect all variables which contribute positively to the pseudo objective value (minimum activity) until we
    * reached the (adjusted) required minimum activity for the inference bound change
    */

   /* first, consider the binary variables */
   if( nvars > 0 )
   {
      SCIP_OBJIMPLICS** minactimpls;

      vars = propdata->minactvars;
      assert(vars != NULL);

      minactimpls = propdata->minactimpls;
      assert(minactimpls != NULL);

#ifndef NDEBUG
      checkGlbfirstnonfixed(scip, propdata);
#endif

      if( infinity )
      {
         /* if the required minimum activity is minus infinity, we have to add all variables which contribute the local
          * pseudo objective activity
          */

         for( v = propdata->glbfirstnonfixed; v < nvars; ++v )
         {
            var = vars[v];
            assert(var != NULL);

            /* @note binary variables can have a zero objective value */

            if( var == infervar )
               continue;

            SCIP_CALL( addConflictBinvar(scip, var, bdchgidx, NULL, NULL, FALSE, &reqpseudoobjval) );
         }
      }
      else
      {
         assert(addedvars != NULL);

         for( v = propdata->glbfirstnonfixed; v < nvars && SCIPisPositive(scip, reqpseudoobjval); ++v )
         {
            var = vars[v];
            assert(var != NULL);

            /* @note binary variables can have a zero objective value */

            if( var == infervar )
               continue;

            if( SCIPhashtableExists(addedvars, (void*)var) )
               continue;

            SCIP_CALL( addConflictBinvar(scip, var, bdchgidx, minactimpls[v], addedvars, propdata->respropuseimplics, &reqpseudoobjval) );
         }
      }
   }

   vars = propdata->objintvars;
   nvars = propdata->nobjintvars;
   assert(nvars == 0 || vars != NULL);

   /* second, consider the non-binary variables */
   for( v = 0; v < nvars && (infinity || SCIPisPositive(scip, reqpseudoobjval)); ++v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(!SCIPisZero(scip, SCIPvarGetObj(var)));

      if( var == infervar )
         continue;

      SCIP_CALL( addConflictBounds(scip, var, bdchgidx, &reqpseudoobjval) );
   }

   return SCIP_OKAY;
}

/** propagates the given variable against the cutoff bound */
static
SCIP_RETCODE propagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variable to propagate */
   int                   inferinfo,          /**< inference information to store with the bound change */
   SCIP_Real             objchg,             /**< objective change */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   SCIP_Bool             local,              /**< local or global propagation */
   SCIP_Bool*            tightened           /**< pointer to store if the variable domain was tightened */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real newbd;
   SCIP_Bool infeasible;

   assert(!SCIPisZero(scip, objchg));
   assert(!SCIPisInfinity(scip, -pseudoobjval));
   assert(!SCIPisInfinity(scip, cutoffbound));
   assert(SCIPisLT(scip, pseudoobjval, cutoffbound) );
   assert(tightened != NULL);

   *tightened = FALSE;

   /* collect bounds of the variable */
   if( local )
   {
      assert(prop != NULL);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
   }
   else
   {
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
   }

   if( SCIPisFeasEQ(scip, lb, ub) )
      return SCIP_OKAY;

   /* depending on the objective contribution we can try to tighten the lower or upper bound of the variable */
   if( objchg > 0.0 )
   {
      newbd = lb + (cutoffbound - pseudoobjval) / objchg;

      if( local )
      {
         SCIP_CALL( SCIPinferVarUbProp(scip, var, newbd, prop, inferinfo, FALSE, &infeasible, tightened) );
         assert(!infeasible);

         if( *tightened ) /* might not be tightened due to numerical reasons */
         {
            SCIPdebugMsg(scip, " -> new (local) upper bound of variable <%s>[%g,%g]: %g, objective change <%g>, pseudo objective <%g>, cutoff bound <%g>\n",
               SCIPvarGetName(var), lb, ub, newbd, objchg, pseudoobjval, cutoffbound);
         }
      }
      else
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newbd, FALSE, &infeasible, tightened) );
         assert(!infeasible);

         if( *tightened )
         {
            SCIPdebugMsg(scip, " -> new (global) upper bound of variable <%s>[%g,%g]: %g, objective change <%g>, pseudo objective <%g>, cutoff bound <%g>\n",
               SCIPvarGetName(var), lb, ub, newbd, objchg, pseudoobjval, cutoffbound);
         }
      }
   }
   else
   {
      newbd = ub + (cutoffbound - pseudoobjval) / objchg;

      if( local )
      {
         SCIP_CALL( SCIPinferVarLbProp(scip, var, newbd, prop, inferinfo, FALSE, &infeasible, tightened) );
         assert(!infeasible);

         if( *tightened ) /* might not be tightened due to numerical reasons */
         {
            SCIPdebugMsg(scip, " -> new (local) lower bound of variable <%s>[%g,%g]: %g, objective change <%g>, pseudo objective <%g>, cutoff bound <%g>\n",
               SCIPvarGetName(var), lb, ub, newbd, objchg, pseudoobjval, cutoffbound);
         }
      }
      else
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newbd, FALSE, &infeasible, tightened) );
         assert(!infeasible);

         if( *tightened )
         {
            SCIPdebugMsg(scip, " -> new (global) lower bound of variable <%s>[%g,%g]: %g, objective change <%g>, pseudo objective <%g>, cutoff bound <%g>\n",
               SCIPvarGetName(var), lb, ub, newbd, objchg, pseudoobjval, cutoffbound);
         }
      }
   }

   return SCIP_OKAY;
}

/** propagates the given binary variable against the cutoff bound */
static
SCIP_RETCODE propagateCutoffboundBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variable to propagate */
   int                   pos,                /**< position of the variable in the corresponding propdata variable array */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   SCIP_Bool*            tightened,          /**< pointer to store if the variable domain was tightened */
   SCIP_Bool*            cutoff,             /**< pointer to store if a cutoff was detected */
   SCIP_Bool             local               /**< propagate local bounds, otherwise global bounds */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_OBJIMPLICS* objimplics;
   SCIP_Real lbobjchg;
   SCIP_Real ubobjchg;
   SCIP_Real objchg;

   assert(SCIPvarIsBinary(var));

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   objimplics = propdata->minactimpls[pos];
   assert(objimplics != NULL);

   /* get objective change in case of fixing the variable to its lower bound (that is zero) */
   SCIP_CALL( getMinactObjchg(scip, var, objimplics, NULL, SCIP_BOUNDTYPE_LOWER, local, &lbobjchg) );
   assert(!SCIPisNegative(scip, lbobjchg));

   /* get objective change in case of fixing the variable to its upper bound (that is one) */
   SCIP_CALL( getMinactObjchg(scip, var, objimplics, NULL, SCIP_BOUNDTYPE_UPPER, local, &ubobjchg) );
   assert(!SCIPisNegative(scip, ubobjchg));

   (*tightened) = FALSE;

   /* nothing can be done if the objective contribution is zero independently of the bound */
   if( SCIPisZero(scip, lbobjchg) && SCIPisZero(scip, ubobjchg) )
      return SCIP_OKAY;

   /* if the lbobjchg and ubobjchg are both able to fix the variable to its upper (1.0) or lower (0.0) bound,
    * respectively, we detected an cutoff
    *
    * @note There is no need to use SCIPisFeasLT() in case the objective is integral since the cutoff bound in that case
    *       is the upper bound minus 1 plus the SCIPcutoffbounddelta() (which is MIN(100.0 * feastol, 0.0001)). However,
    *       if the objective is not integral we have to check w.r.t. an epsilon to avoid numerical problems.
    */
   if( SCIPisFeasLT(scip, cutoffbound, pseudoobjval + ubobjchg) && SCIPisFeasLT(scip, cutoffbound, pseudoobjval + lbobjchg) )
   {
      /* check if conflict analysis is applicable */
      if( local && SCIPisConflictAnalysisApplicable(scip) )
      {
         assert(SCIPgetDepth(scip) > 0);

         /* initialize conflict analysis */
         SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, TRUE) );

         /* add all variable whose best bound changes increased the pseudo objective value above to cutoff bound */
         SCIP_CALL( resolvePropagation(scip, propdata, pseudoobjval, NULL, -1, SCIP_BOUNDTYPE_UPPER, NULL) );

         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
      }

      (*cutoff) = TRUE;
   }
   else
   {
      /* try to tighten the variable bound use the larger objective contribution */
      if( lbobjchg > ubobjchg )
         objchg = -lbobjchg;
      else
         objchg = ubobjchg;

      SCIP_CALL( propagateCutoffboundVar(scip, prop, var, pos, objchg, cutoffbound, pseudoobjval, local, tightened) );
   }

   return SCIP_OKAY;
}

/** globally propagates if a new cutoff bound or global pseudo objective value (minimum activity of the objective
 *  function) is available
 */
static
SCIP_RETCODE propagateCutoffboundGlobally(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** minactvars;
   SCIP_VAR** objintvars;
   SCIP_VAR* var;
   SCIP_Bool tightened;
   SCIP_Real pseudoobjval;
   SCIP_Real cutoffbound;
   int nminactvars;
   int nobjintvars;
   int v;

   /* this method should not be called in the root node of the search tree since the standard propagation already does
    * the job
    */
   assert(SCIPgetDepth(scip) > 0);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   pseudoobjval = SCIPgetGlobalPseudoObjval(scip);
   cutoffbound = propdata->cutoffbound;

   /* nothing can be done if the global pseudo objective is minus infinity */
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;

   /* check if the global pseudo objective value (minimum activity of the objective function) is greater or equal to
    * the cutoff bound
    */
   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   minactvars = propdata->minactvars;
   objintvars = propdata->objintvars;
   nminactvars = propdata->nminactvars;
   nobjintvars = propdata->nobjintvars;

#ifndef NDEBUG
   checkGlbfirstnonfixed(scip, propdata);
#endif

   *cutoff = FALSE;

   /* always propagate the binary variables completely */
   for( v = propdata->glbfirstnonfixed; v < nminactvars; ++v )
   {
      var = minactvars[v];
      assert(var != NULL);

      /* check if the variables is already globally fixed; if so continue with the potential candidate */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, cutoff, FALSE) );

      /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
       * contribution w.r.t. minimum activity (pseudo objective value) of the objective function; these values are the
       * increase in the pseudo objective activity we would get if we fix the variable to its worse bound; hence, we can
       * stop if for a variable this potential increase is not enough too exceed the cutoff bound;
       */
      if( !tightened )
      {
         SCIPdebugMsg(scip, "interrupt global pseudo objective propagation w.r.t. cutoff bound <%.15g> for binary variables after %d from %d binary variables\n",
            cutoffbound, v, nminactvars);
         break;
      }

      if( *cutoff )
         return SCIP_OKAY;

      /* @note The variable might not be globally fixed right away since this would destroy the local internal
       *       data structure of a search node; the bound change is in that case pending; hence we cannot assert
       *       that the variable is globally fixed
       */
      (*nchgbds)++;
   }
   propdata->glbfirstnonfixed = v;
   propdata->firstnonfixed = MAX(propdata->firstnonfixed, v);

   /* check all binary variables which could potentially be fixed */
   for( ; v < nminactvars && cutoffbound - pseudoobjval <  propdata->minactimpls[v]->maxobjchg; ++v )
   {
      var = minactvars[v];
      assert(var != NULL);

      /* check if the variables is already globally fixed; if so continue with the potential candidate */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, cutoff, FALSE) );

      /* check if the domain of the variable was reduced */
      if( tightened )
         (*nchgbds)++;

      if( *cutoff )
         return SCIP_OKAY;
   }

#if 0 /* might fail, but is not a real error, still need to investigate */
#ifndef NDEBUG
   /* check that the abort criteria for the binary variables works */
   for( ; v < nminactvars; ++v )
   {
      assert(cutoffbound - pseudoobjval >=  propdata->minactimpls[v]->maxobjchg);

      var = minactvars[v];
      assert(var != NULL);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, cutoff, FALSE) );
      assert(!tightened);
      assert(!(*cutoff));
   }
#endif
#endif

   /* propagate the non-binary variables completely */
   for( v = 0; v < nobjintvars; ++v )
   {
      var = objintvars[v];
      assert(var != NULL);
      assert(!SCIPisZero(scip, SCIPvarGetObj(var)));

      /* try to tighten the bound of the variable */
      SCIP_CALL( propagateCutoffboundVar(scip, NULL, var, -1, SCIPvarGetObj(var), cutoffbound, pseudoobjval, FALSE, &tightened) );

      /* check if the domain of the variable was reduced */
      if( tightened )
         (*nchgbds)++;
   }

   propdata->glbpropagated = TRUE;

   return SCIP_OKAY;
}

/** propagates the cutoff bound for binary variables (c*x <= cutoff) */
static
SCIP_RETCODE propagateCutoffboundBinvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** minactvars;
   SCIP_VAR* var;
   SCIP_Bool tightened;
   int nminactvars;
   int v;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   minactvars = propdata->minactvars;
   nminactvars = propdata->nminactvars;
   assert(nminactvars == 0 || minactvars != NULL);

   /* always propagate the binary variables completely; note that the variables before the firstnonfixed indexed are
    * already locally fixed and those before glbfirstnonfixed are already globally fixed
    */

#ifndef NDEBUG
   /* check that the variables before glbfirstnonfixed are globally fixed */
   checkGlbfirstnonfixed(scip, propdata);

   /* check that the variables before firstnonfixed are locally fixed */
   for( v = propdata->glbfirstnonfixed; v < propdata->firstnonfixed; ++v )
   {
      var =  minactvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);
   }
#endif

   (*cutoff) = FALSE;

   for( v = propdata->firstnonfixed; v < nminactvars; ++v )
   {
      var =  minactvars[v];
      assert(var != NULL);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, cutoff, TRUE) );

      /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
       * contribution w.r.t. minimum activity of the objective function; These values are the increase in the pseudo
       * objective activity (minimum activity of the objective function) we would get if we fix the variable to its
       * worse bound; hence, we can stop if for a variable this potential increase is not enough too exceed the cutoff
       * bound;
       */
      if( !tightened )
      {
         SCIPdebugMsg(scip, "interrupt local pseudo objective propagation w.r.t. cutoff bound <%.15g> for binary variables after %d from %d binary variables\n",
            cutoffbound, v, nminactvars);
         break;
      }

      if( *cutoff )
         return SCIP_OKAY;

      /* check that the binary variable is locally fixed */
      assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);
      (*nfixedvars)++;
   }
   propdata->firstnonfixed = v;

   /* check all binary variables which could potentially be fixed */
   for( ; v < nminactvars && cutoffbound - pseudoobjval < propdata->minactimpls[v]->maxobjchg; ++v )
   {
      var =  minactvars[v];
      assert(var != NULL);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, cutoff, TRUE) );

      if( tightened )
      {
         assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);
         (*nfixedvars)++;
      }

      if( *cutoff )
         return SCIP_OKAY;
   }

#if 0 /* might fail, but is not a real error, still need to investigate */
#ifndef NDEBUG
   /* check that the abort criteria for the binary variables works */
   for( ; v < nminactvars; ++v )
   {
      var = minactvars[v];
      assert(var != NULL);

      assert(cutoffbound - pseudoobjval >= propdata->minactimpls[v]->maxobjchg);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, cutoff, TRUE) );
      assert(!tightened);
      assert(!(*cutoff));
   }
#endif
#endif

   return SCIP_OKAY;
}

/** propagates the cutoff bound c*x <= cutoff */
static
SCIP_RETCODE propagateCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_Real pseudoobjval;
   SCIP_Real cutoffbound;
   SCIP_Bool cutoff;
   SCIP_Bool tightened;
   int nchgbds;

   assert(result != NULL);

   (*result) = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* get current pseudo objective value (minimum activity of the objective function) and cutoff bound */
   pseudoobjval = SCIPgetPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;
   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   /* @note A new global pseudo objective value could be used to retrieve global fixings. There is, however, no need to
    *       check if a new global pseudo objective value is available. This is the case since a new (better) global
    *       pseudo activity implies that a global bound change was performed. That causes that the root node of the
    *       search tree gets marked for repropagation. That will result in a propagation call of the pseudo objective
    *       propagator.
    */

   /* check current cutoff bound */
   if( cutoffbound < propdata->cutoffbound )
   {
      propdata->glbpropagated = FALSE;
      propdata->cutoffbound = cutoffbound;
   }

   nchgbds = 0;
   cutoff = FALSE;
   (*result) = SCIP_DIDNOTFIND;

   /* check if we have a new cutoff bound; in that case we globally propagate this new bound
    *
    * @note there is no need to propagate the cutoff bound if we are in the root node since this will be done by the
    *       standard local propagation
    */
   if( propdata->propcutoffbound && !propdata->glbpropagated && SCIPgetDepth(scip) > 0 )
   {
      /* first globally propagate new cutoff bound or pseudo objective activity */
      SCIP_CALL( propagateCutoffboundGlobally(scip, prop, &nchgbds, &cutoff) );

      if( cutoff )
      {
         /* we are done with solving since a global pseudo activity is greater or equal to the cutoff bound */
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );

         (*result) = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* check if the global propagation cut off the active/current node */
      if( SCIPgetCutoffdepth(scip) <= SCIPgetDepth(scip) )
      {
         (*result) = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   /* check if the pseudo objective value (minimum activity of the objective function) is greater or equal to the cutoff
    * bound
    */
   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      SCIPdebugMsg(scip, "pseudo objective value <%g> exceeds cutoff bound <%g>\n", pseudoobjval, cutoffbound);

      /* check if conflict analysis is applicable */
      if( SCIPisConflictAnalysisApplicable(scip) )
      {
         assert(SCIPgetDepth(scip) > 0);

         /* initialize conflict analysis */
         SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, TRUE) );

         /* add all variable whose best bound changes increased the pseudo objective value above the cutoff bound */
         SCIP_CALL( resolvePropagation(scip, propdata, cutoffbound, NULL, -1, SCIP_BOUNDTYPE_UPPER, NULL) );

         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
      }
      (*result) = SCIP_CUTOFF;

      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "propagating pseudo objective function (pseudoobj: %g, cutoffbound: %g)\n", pseudoobjval, cutoffbound);

   /* propagate binary variables */
   SCIP_CALL( propagateCutoffboundBinvars(scip, prop, cutoffbound, pseudoobjval, &nchgbds, &cutoff) );

   if( cutoff )
   {
      (*result) = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* tighten domains of non-binary variables, if they would increase the pseudo objective value above the cutoff
    * bound
    */
   if( propdata->propfullinroot && SCIPgetDepth(scip) == 0 )
   {
      SCIP_VAR** objintvars;
      SCIP_VAR* var;
      SCIP_Real objval;
      int nobjintvars;
      int v;

      objintvars = propdata->objintvars;
      nobjintvars = propdata->nobjintvars;
      assert(nobjintvars == 0 || objintvars != NULL);

      /* propagate all non-binary variables */
      for( v = 0; v < nobjintvars; ++v )
      {
         var = objintvars[v];
         assert(var != NULL);

         objval = SCIPvarGetObj(var);
         assert(!SCIPisZero(scip, objval));

         /* try to tighten the bound of the variable */
         SCIP_CALL( propagateCutoffboundVar(scip, NULL, var, -1, objval, cutoffbound, pseudoobjval, FALSE, &tightened) );

         /* check if the domain of the variable was reduced */
         if( tightened )
            nchgbds++;
      }
   }
   else
   {
      SCIP_VAR** objintvars;
      SCIP_VAR* var;
      SCIP_Real objval;
      int nobjintvars;
      int nmaxuseless;
      int nuseless;
      int c;
      int v;

      objintvars = propdata->objintvars;
      nobjintvars = propdata->nobjintvars;
      assert(nobjintvars == 0 || objintvars != NULL);

      /* compute maximum number of useless propagations before aborting */
      nmaxuseless = MAX(propdata->minuseless, (int)propdata->maxvarsfrac*(nobjintvars));

      nuseless = 0;
      v = propdata->lastvarnum;

      for( c = 0; c < nobjintvars && nuseless < nmaxuseless; ++c )
      {
         v++;
         if( v >= nobjintvars )
            v = 0;

         var = objintvars[v];
         assert(var != NULL);

         objval = SCIPvarGetObj(var);
         assert(!SCIPisZero(scip, objval));

         /* try to tighten the bound of the variable */
         SCIP_CALL( propagateCutoffboundVar(scip, prop, var, -1, objval, cutoffbound, pseudoobjval, TRUE, &tightened) );

         /* check if the domain of the variable was reduced */
         if( tightened )
            nchgbds++;
         else
            nuseless++;
      }
      propdata->lastvarnum = v;
   }

   /* check if we chanced bounds */
   if( nchgbds > 0 )
      (*result) = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}

/** recalculates the maximum objective pseudoactivity */
static
void calcMaxObjPseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR** vars;
   SCIP_Real objval;
   SCIP_Real contrib;
   int nvars;
   int v;

   assert(propdata != NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* calculate current max pseudo activity and largest contribution */
   propdata->maxpseudoobjact = 0.0;
   propdata->maxpseudoobjactinf = 0;

   for( v = 0; v < nvars; ++v )
   {
      objval = SCIPvarGetObj(vars[v]);
      if( SCIPisPositive(scip, objval) )
      {
         contrib = SCIPvarGetUbGlobal(vars[v]);
         if( !SCIPisInfinity(scip, contrib) )
            contrib *= objval;
      }
      else if( SCIPisNegative(scip, objval) )
      {
         contrib = SCIPvarGetLbGlobal(vars[v]);
         if( !SCIPisInfinity(scip, -contrib) )
            contrib *= objval;
         else
            contrib *= -1.0;
      }
      else
         continue;

      if( SCIPisInfinity(scip, contrib) )
         propdata->maxpseudoobjactinf++;
      else
         propdata->maxpseudoobjact += contrib;
   }
}

/** updates the pseudo objective activity if necessary */
static
void updateMaxObjPseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   /* if necessary, calculate the maximum pseudo objective activity */
   if( propdata->maxpseudoobjact == SCIP_INVALID ) /*lint !e777*/
      calcMaxObjPseudoactivity(scip, propdata);
   assert(propdata->maxpseudoobjact != SCIP_INVALID); /*lint !e777*/
}

/** returns the residual pseudo objective activity without the given value */
static
SCIP_Real getMaxObjPseudoactivityResidualValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Real             contrib             /**< value to eliminate from pseudo objective activity */
   )
{
   SCIP_Real residual;

   assert(propdata != NULL);

   /* if necessary, calculate the maximum pseudo objective activity */
   if( propdata->maxpseudoobjact == SCIP_INVALID ) /*lint !e777*/
      calcMaxObjPseudoactivity(scip, propdata);
   assert(propdata->maxpseudoobjact != SCIP_INVALID); /*lint !e777*/

   if( SCIPisInfinity(scip, contrib) )
   {
      assert(propdata->maxpseudoobjactinf >= 1);
      /* check if this variable yields the only infinite contribution */
      if( propdata->maxpseudoobjactinf == 1 )
         residual = propdata->maxpseudoobjact;
      else
         residual = SCIPinfinity(scip);
   }
   else
   {
      /* check if there is an infinite contribution */
      if( propdata->maxpseudoobjactinf >= 1 )
         residual = SCIPinfinity(scip);
      else
         residual = propdata->maxpseudoobjact - contrib;
   }

   return residual;
}

/** returns the residual pseudo objective activity */
static
SCIP_Real getMaxObjPseudoactivityResidual(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var                 /**< variable to get residual activity for */
   )
{
   SCIP_Real objval;
   SCIP_Real contrib;

   assert(propdata != NULL);

   objval = SCIPvarGetObj(var);
   if( SCIPvarGetWorstBoundType(var) == SCIP_BOUNDTYPE_UPPER )
   {
      contrib = SCIPvarGetUbGlobal(var);
      if( !SCIPisInfinity(scip, contrib) )
         contrib *= objval;
   }
   else
   {
      assert(SCIPvarGetWorstBoundType(var) == SCIP_BOUNDTYPE_LOWER);
      contrib = SCIPvarGetLbGlobal(var);
      if( !SCIPisInfinity(scip, -contrib) )
         contrib *= objval;
      else
         contrib *= -1.0;
   }

   return getMaxObjPseudoactivityResidualValue(scip, propdata, contrib);
}

/** returns the maximum pseudo objective activity of the objective function */
static
SCIP_Real getMaxObjPseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   return getMaxObjPseudoactivityResidualValue(scip, propdata, 0.0);
}

/** propagates the global domain of the given binary variable against the lower bound (c*x >= lowerbound) */
static
SCIP_RETCODE propagateLowerboundBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             lowerbound,         /**< lower bound to use */
   SCIP_Real             maxpseudoobjact,    /**< maximum pseudo objective activity */
   SCIP_Bool             useimplics,         /**< should implications be used */
   SCIP_Bool*            infeasible,         /**< pointer to store if the variable domain got empty, infeasible */
   SCIP_Bool*            tightened           /**< pointer to store if the variable domain was tightened */
   )
{
   SCIP_Real lbobjchg;
   SCIP_Real ubobjchg;

   assert(SCIPvarIsBinary(var));
   assert(SCIPisDualfeasLE(scip, lowerbound, maxpseudoobjact));
   assert(!SCIPisInfinity(scip, maxpseudoobjact));

   /*@todo Instead of running always over all implications use SCIP_OBJIMPLICS in the same way as for the propagation of
    *      the minimum activity against the cutoff bound. If so we could use the cliques as well.
    */

   /* get contribution of variable by fixing it to its lower bound w.r.t. maximum activity of the objective function */
   SCIP_CALL( getMaxactObjchg(scip, var, SCIP_BOUNDTYPE_LOWER, useimplics, &lbobjchg) );
   assert(!SCIPisPositive(scip, lbobjchg));

   /* get contribution of variable by fixing it to its upper bound w.r.t. maximum activity of the objective function */
   SCIP_CALL( getMaxactObjchg(scip, var, SCIP_BOUNDTYPE_UPPER, useimplics, &ubobjchg) );
   assert(!SCIPisPositive(scip, ubobjchg));

   (*infeasible) = FALSE;
   (*tightened) = FALSE;

   /* if the maximum activity of the objective function without the contribution of the given variable shrinks below the
    * global lower bound, the contribution of the variable is need; hence, we can fix it to corresponding bound globally
    */
   if( SCIPisFeasLT(scip, maxpseudoobjact + lbobjchg, lowerbound) && SCIPisFeasLT(scip, maxpseudoobjact + ubobjchg, lowerbound) )
   {
      /* fixing the variable to zero or one leads to decreases of the maximum activity below the lower bound, hence, we
       * detected an cutoff
       */
      (*infeasible) = TRUE;
   }
   else if( SCIPisFeasLT(scip, maxpseudoobjact + lbobjchg, lowerbound) )
   {
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, 1.0, FALSE, infeasible, tightened) );
   }
   else if( SCIPisFeasLT(scip, maxpseudoobjact + ubobjchg, lowerbound) )
   {
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, 0.0, FALSE, infeasible, tightened) );
   }

   return SCIP_OKAY;
}

/** propagates the global domains of the given variable with non-zero objective coefficient against the lower bound
 *  (c*x >= lowerbound)
 */
static
SCIP_RETCODE propagateLowerboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             lowerbound,         /**< lower bound to use */
   SCIP_Bool*            infeasible,         /**< pointer to store if the variable domain got empty, infeasible */
   SCIP_Bool*            tightened           /**< pointer to store if the variable domain was tightened */
   )
{
   SCIP_Real residual;
   SCIP_Real newbd;
   SCIP_Real objval;

   objval = SCIPvarGetObj(var);
   assert(!SCIPisZero(scip, objval));

   (*tightened) = FALSE;

   /* get residual pseudo objective activity, that is the pseudo objective activity without the given variable */
   residual = getMaxObjPseudoactivityResidual(scip, propdata, var);

   if( SCIPisInfinity(scip, residual) )
      return SCIP_OKAY;

   /* compute potential mew bound */
   newbd = (lowerbound - residual) / objval;

   /**@note In case the variable is integral we force the update of the new bound */

   if( objval > 0.0 )
   {
      SCIP_Real lb;

      lb = SCIPvarGetLbGlobal(var);

      if( !SCIPvarIsIntegral(var) )
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newbd, FALSE, infeasible, tightened) );
      }
      else if( SCIPisFeasGT(scip, newbd, lb) )
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newbd, TRUE, infeasible, tightened) );
      }
   }
   else
   {
      SCIP_Real ub;

      ub = SCIPvarGetUbGlobal(var);

      if( !SCIPvarIsIntegral(var) )
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newbd, FALSE, infeasible, tightened) );
      }
      else if( SCIPisFeasLT(scip, newbd, ub) )
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newbd, TRUE, infeasible, tightened) );
      }
   }

   return SCIP_OKAY;
}

/** propagates the global lower (dual) bound c*x >= lowerbound */
static
SCIP_RETCODE propagateLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Real lowerbound;
   SCIP_Real maxpseudoobjact;
   SCIP_Bool cutoff;
   int nchgbds;

   assert(result != NULL);

   (*result) = SCIP_DIDNOTRUN;
   cutoff = FALSE;
   nchgbds = 0;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->nminactvars > 0 || propdata->nobjintvars > 0);

   /* check if there is a chance to find a reduction */
   lowerbound = SCIPgetLowerbound(scip);

   if( SCIPisInfinity(scip, -lowerbound) )
      return SCIP_OKAY;

   /* if the lower bound did not change since the last propagation as well as the global bounds of the variables with a
    * non-zero objective coefficient we do nothing since there is no new information available
    */
   if( SCIPisLE(scip, lowerbound, propdata->lastlowerbound) && propdata->maxpseudoobjact < SCIP_INVALID )
      return SCIP_OKAY;

   /* updates the pseudo objective activity if necessary */
   updateMaxObjPseudoactivity(scip, propdata);

   /* if more than one variable contributes an infinity to the maximal pseudo activity we can do nothing */
   assert(propdata->maxpseudoobjact < SCIP_INVALID);
   if( propdata->maxpseudoobjactinf > 1 )
      return SCIP_OKAY;

   maxpseudoobjact = getMaxObjPseudoactivity(scip, propdata);
   assert(!SCIPisInfinity(scip, maxpseudoobjact) || !SCIPisInfinity(scip, lowerbound));

#ifndef NDEBUG
   /* check that the global indices are correct */
   checkGlbfirstnonfixed(scip, propdata);
#endif

   /* if the maximum pseudo objective activity is smaller than the lower bound the problem is infeasible */
   if( SCIPisDualfeasLT(scip, maxpseudoobjact, lowerbound) )
      cutoff = TRUE;
   else
   {
      SCIP_VAR** objintvars;
      SCIP_VAR* var;
      SCIP_Bool tightened;
      int nobjintvars;
      int v;

      if( propdata->maxpseudoobjactinf == 0 && !SCIPisInfinity(scip, maxpseudoobjact) )
      {
         SCIP_VAR** maxactvars;
         int nmaxactvars;

         maxactvars = propdata->maxactvars;
         nmaxactvars = propdata->nmaxactvars;
         assert(nmaxactvars == 0 || maxactvars != NULL);

         for( v = propdata->maxactfirstnonfixed; v < nmaxactvars; ++v )
         {
            var = maxactvars[v];
            assert(var != NULL);
            assert(SCIPvarIsBinary(var));

            /* check if the variables is already globally fixed; if so continue with the next potential candidate */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
               continue;

            /* try to propagate variable domain globally */
            SCIP_CALL( propagateLowerboundBinvar(scip, var, lowerbound, maxpseudoobjact, propdata->propuseimplics, &cutoff, &tightened) );

            /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
             * contribution w.r.t. maximum activity of the objective function; These values are the decrease we would
             * get with the maximum pseudo objective activity if we fix the variable to its best bound; hence, we can
             * stop if for a variable this potential decrease is not enough anymore to fall below the lower bound.
             *
             * @note In case a fixing was performed. The variable might not be globally fixed right away since this would
             *       destroy the local internal data structure of a search node; the bound change is in that case pending;
             *       hence we cannot assert that the variable is globally fixed
             */
            if( !tightened )
            {
               assert(!SCIPisInfinity(scip, propdata->maxpseudoobjact));
               SCIPdebugMsg(scip, "interrupt pseudo objective propagation w.r.t. lower bound <%.15g> for binary variables after %d from %d binary variables\n",
                  lowerbound, v, nmaxactvars);
               break;
            }

            if( cutoff )
               break;

            /* update maximum pseudo activity since the previous global bound change might invalidated the maximum
             * pseudo activity
             */
            maxpseudoobjact = getMaxObjPseudoactivity(scip, propdata);
            nchgbds++;
         }

         /* update globally fixed index if abort criteria was applied */
         propdata->maxactfirstnonfixed = v;

         /* check all binary variables which could potentially be fixed */
         for( ; v < nmaxactvars && maxpseudoobjact - lowerbound < propdata->maxactchgs[v] && !cutoff; ++v )
         {
            var =  maxactvars[v];
            assert(var != NULL);
            assert(SCIPvarIsBinary(var));

            /* check if the variables is already globally fixed; if so continue with the potential candidate */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
               continue;

            /* propagates the cutoff bound for the given binary variable */
            SCIP_CALL( propagateLowerboundBinvar(scip, var, lowerbound, maxpseudoobjact, propdata->propuseimplics, &cutoff, &tightened) );

            if( tightened )
            {
               /* update maximum pseudo activity since the previous global bound change might invalidated the maximum
                * pseudo activity
                */
               maxpseudoobjact = getMaxObjPseudoactivity(scip, propdata);
               nchgbds++;
            }
         }

#if 0 /* might fail, but is not a real error, still need to investigate */
#ifndef NDEBUG
         /* check that the abort criteria for the binary variables works */
         for( ; v < nmaxactvars && !cutoff; ++v )
         {
            var = maxactvars[v];
            assert(var != NULL);
            assert(SCIPvarIsBinary(var));

            /* check if the variables is already globally fixed; if so continue with the next potential candidate */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
               continue;

            /* try to propagate variable domain globally */
            SCIP_CALL( propagateLowerboundBinvar(scip, var, lowerbound, maxpseudoobjact, propdata->propuseimplics, &cutoff, &tightened) );
            assert(!tightened);
            assert(!cutoff);
         }
#endif
#endif
      }

      objintvars = propdata->objintvars;
      nobjintvars = propdata->nobjintvars;
      assert(nobjintvars == 0 || objintvars != NULL);

      /* propagate c*x >= lowerbound for non-binary variables */
      for( v = 0; v < nobjintvars && !cutoff; ++v )
      {
         var = objintvars[v];
         assert(var != NULL);
         assert(!SCIPisZero(scip, SCIPvarGetObj(var)));

         /* try to propagate variable domain globally */
         SCIP_CALL( propagateLowerboundVar(scip, propdata, var, lowerbound, &cutoff, &tightened) );

         if( tightened )
            nchgbds++;
      }
   }

   /* evaluate propagation results */
   if( cutoff )
   {
      /* we are done with solving since a global bound change is infeasible: cutoff root node */
      SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
      (*result) = SCIP_CUTOFF;
   }
   else if( nchgbds > 0 )
      (*result) = SCIP_REDUCEDDOM;

   /* remember the lower bound which we already propagated */
   propdata->lastlowerbound = lowerbound;

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyPseudoobj)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropPseudoobj(scip) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreePseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   SCIPfreeBlockMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolPseudoobj)
{
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* if active pricer are present we want to catch the variable added event */
   if( SCIPgetNActivePricers(scip) > 0 )
   {
      assert(!propdata->catchvaradded);
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, propdata->eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
      propdata->catchvaradded = TRUE;
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   if( propdata->catchvaradded )
   {
      /* drop the variable added event */
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, propdata->eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      propdata->catchvaradded = FALSE;
   }

   /* free propagator data */
   SCIP_CALL( propdataExit(scip, propdata) );

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolPseudoobj)
{  /*lint --e{715}*/

   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_Real cutoffbound;
   SCIP_Real pseudoobjval;
   int oldnchgbds;
   int nvars;
   int v;

   assert(result != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   (*result) = SCIP_DIDNOTRUN;

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* do nothing if objective propagation is not allowed */
   if( !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   pseudoobjval = SCIPgetGlobalPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;

   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      (*result) = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* only propagate if a new cutoff bound or global pseudo objective value is available */
   if( cutoffbound < propdata->cutoffbound || pseudoobjval > propdata->glbpseudoobjval )
   {
      SCIP_Real objval;
      SCIP_Bool tightened;

      (*result) = SCIP_DIDNOTFIND;
      oldnchgbds = *nchgbds;

      /* get the problem variables */
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      /* scan the variables for pseudoobj bound reductions
       * (loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array)
       */
      for( v = nvars - 1; v >= 0; --v )
      {
         objval = SCIPvarGetObj(vars[v]);

         if( SCIPisZero(scip, objval) )
            continue;

         SCIP_CALL( propagateCutoffboundVar(scip, NULL, vars[v], -1, objval, cutoffbound, pseudoobjval, FALSE, &tightened) );

         if( tightened )
            (*nchgbds)++;
      }

      /* evaluate if bound change was detected */
      if( *nchgbds > oldnchgbds )
         (*result) = SCIP_SUCCESS;

      /* store the old values */
      propdata->cutoffbound = cutoffbound;
      propdata->glbpseudoobjval = pseudoobjval;
      propdata->glbpropagated = TRUE;
   }

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   (*result) = SCIP_DIDNOTRUN;

   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   if( proptiming == SCIP_PROPTIMING_DURINGLPLOOP && SCIPgetDepth(scip) != 0 )
      return SCIP_OKAY;

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* do not run if propagation w.r.t. objective is not allowed */
   if( !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   /* check if enough new variable are added (due to column generation to reinitialized the propagator data) */
   if( !propdata->initialized || propdata->nnewvars > propdata->maxnewvars )
   {
      /* free current propdata data */
      SCIP_CALL( propdataExit(scip, propdata) );

      /* initialize propdata data from scratch */
      SCIP_CALL( propdataInit(scip, propdata) );
   }

   /* nothing to do for zero objective */
   if( propdata->nminactvars == 0 && propdata->nmaxactvars == 0 && propdata->nobjintvars == 0 )
      return SCIP_OKAY;

   /* propagate c*x <= cutoff */
   SCIP_CALL( propagateCutoffbound(scip, prop, result) );

   if( (*result) != SCIP_CUTOFF && (propdata->nmaxactvars > 0 || propdata->nobjintvars > 0) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_RESULT dualresult;

      /* propagates the global lower (dual) bound c*x >= lowerbound */
      SCIP_CALL( propagateLowerbound(scip, prop, &dualresult) );

      if( dualresult == SCIP_REDUCEDDOM || dualresult == SCIP_CUTOFF )
         (*result) = dualresult;
   }

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Real cutoffbound;

   assert(!SCIPisEQ(scip, SCIPvarGetLbGlobal(infervar), SCIPvarGetUbGlobal(infervar)));

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   cutoffbound = SCIPgetCutoffbound(scip);
   assert(!SCIPisInfinity(scip, cutoffbound));
   assert(infervar != NULL);

   SCIPdebugMsg(scip, "resolve bound change <%s> %s <%g>(%g), cutoff bound <%g>\n", SCIPvarGetName(infervar),
      boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE),
      SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, FALSE), cutoffbound);

   /* resolve the propagation of the inference variable w.r.t. required minactivity */
   SCIP_CALL( resolvePropagation(scip, propdata, cutoffbound, infervar, inferinfo, boundtype, bdchgidx) );

   (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * Event handler
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_EVENTTYPE eventtype;

   propdata = (SCIP_PROPDATA*)eventdata;
   assert(propdata != NULL);

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   eventtype = SCIPeventGetType(event);

   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if bound got relaxed we need to start up front for trial of bound tightening */
      propdata->firstnonfixed = 0;
      break;
   case SCIP_EVENTTYPE_VARADDED:
      propdata->nnewvars++;
      break;
   default:
      assert(eventtype == SCIP_EVENTTYPE_GLBCHANGED || eventtype == SCIP_EVENTTYPE_GUBCHANGED);

      /* invalidate the maximum pseudo objective activity */
      propdata->maxpseudoobjact = SCIP_INVALID;
      propdata->maxpseudoobjactinf = 0;
   }

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the pseudo objective function propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropPseudoobj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;


   /* create pseudoobj propagator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   /* reset propagator data structure */
   propdataReset(scip, propdata);

   propdata->eventhdlr = NULL;
   /* include event handler for gloabl bound change events and variable added event (in case of pricing) */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &propdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecPseudoobj, NULL) );

   if( propdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for pseudo objective propagator not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecPseudoobj,
         propdata) );
   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyPseudoobj) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreePseudoobj) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolPseudoobj) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolPseudoobj) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolPseudoobj, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropPseudoobj) );

   /* add pseudoobj propagator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/minuseless",
         "minimal number of successive non-binary variable propagator whithout a bound reduction before aborted",
         &propdata->minuseless, TRUE, DEFAULT_MINUSELESS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "propagating/" PROP_NAME "/maxvarsfrac",
         "maximal fraction of non-binary variables with non-zero objective without a bound reduction before aborted",
         &propdata->maxvarsfrac, TRUE, DEFAULT_MAXVARSFRAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/propfullinroot",
         "do we want to propagate all non-binary variables if we are propagating the root node",
         &propdata->propfullinroot, TRUE, DEFAULT_PROPFULLINROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/propcutoffbound",
         "propagate new cutoff bound directly globally",
         &propdata->propcutoffbound, TRUE, DEFAULT_PROPCUTOFFBOUND, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/force",
         "should the propagator be forced even if active pricer are present?",
         &propdata->force, TRUE, DEFAULT_FORCE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxnewvars",
         "number of variable added after the propgatore is reinitialized?",
         &propdata->maxnewvars, TRUE, DEFAULT_MAXNEWVARS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/propuseimplics",
         "use implications to strengthen the propagation of binary variable (increasing the objective change)?",
         &propdata->propuseimplics, TRUE, DEFAULT_PROPUSEIMPLICS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/respropuseimplics",
         "use implications to strengthen the resolve propagation of binary variable (increasing the objective change)?",
         &propdata->respropuseimplics, TRUE, DEFAULT_RESPROPUSEIMPLICS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maximplvars",
         "maximum number of binary variables the implications are used if turned on (-1: unlimited)?",
         &propdata->maximplvars, TRUE, DEFAULT_MAXIMPLVARS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** propagates the cutoff bound for the given variables */
SCIP_RETCODE SCIPpropagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variables to propagate */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   SCIP_Bool*            tightened           /**< pointer to if the domain was tightened */
   )
{
   SCIP_Real objval;

   objval = SCIPvarGetObj(var);

   SCIP_CALL( propagateCutoffboundVar(scip, prop, var, -1, objval, cutoffbound, pseudoobjval, TRUE, tightened) );

   return SCIP_OKAY;
}
