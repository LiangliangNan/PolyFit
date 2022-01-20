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

/**@file    prop_genvbounds.c
 * @ingroup PROPAGATORS
 * @brief   generalized variable bounds propagator
 * @author  Stefan Weltge
 * @author  Ambros Gleixner
 * @author  Benjamin Mueller
 */

/**@todo should we only discard events catched from nodes that are not the current node's ancestors? */
/**@todo improve computation of minactivity */
/**@todo for multaggr vars on left-hand side, create a linear constraint, probably in exitpre */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_genvbounds.h"
#include "scip/debug.h"
#include "scip/cons_linear.h"

#define PROP_NAME                            "genvbounds"
#define PROP_DESC                            "generalized variable bounds propagator"
#define PROP_TIMING   SCIP_PROPTIMING_ALWAYS
#define PROP_PRIORITY                3000000 /**< propagator priority */
#define PROP_FREQ                          1 /**< propagator frequency */
#define PROP_DELAY                     FALSE /**< should propagation method be delayed, if other propagators
                                              *   found reductions? */
#define PROP_PRESOL_PRIORITY        -2000000 /**< priority of the presolving method (>= 0: before, < 0: after
                                              *   constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates
                                              *   in (-1: no limit) */
#define DEFAULT_GLOBAL_PROPAGATION      TRUE /**< apply global propagation? */
#define DEFAULT_PROPAGATE_IN_ROOT_NODE  TRUE /**< apply genvbounds in root node if no new incumbent was found? */
#define DEFAULT_SORT                    TRUE /**< sort genvbounds and wait for bound change events? (otherwise all
                                              *   genvbounds are applied in each node) */
#define DEFAULT_PROPASCONSS            FALSE /**< should genvbounds be transformed to (linear) constraints? */

#define EVENTHDLR_NAME                       "genvbounds"
#define EVENTHDLR_DESC                       "event handler for generalized variable bounds propagator"


/*
 * Data structures
 */

/** GenVBound data */
struct GenVBound
{
   SCIP_VAR**            vars;               /**< pointers to variables x_j occuring in this generalized variable
                                              *   bound */
   SCIP_VAR*             var;                /**< pointer to variable x_i */
   SCIP_Real*            coefs;              /**< coefficients a_j of the variables listed in vars */
   SCIP_Real             constant;           /**< constant term in generalized variable bound */
   SCIP_Real             cutoffcoef;         /**< cutoff bound's coefficient */
   int                   coefssize;          /**< size of coefs array */
   int                   index;              /**< index of this genvbound in genvboundstore array */
   int                   ncoefs;             /**< number of nonzero coefficients a_j */
   SCIP_BOUNDTYPE        boundtype;          /**< type of bound provided by the genvbound, SCIP_BOUNDTYPE_LOWER/UPPER
                                              *   if +/- x_i on left-hand side */
};
typedef struct GenVBound GENVBOUND;

/** starting indices data structure */
struct SCIP_EventData
{
   SCIP_PROP*            prop;               /**< pointer to genvbounds propagator */
   SCIP_VAR*             var;                /**< variable */
   int*                  startindices;       /**< array to store the first indices of genvbounds in components that are
                                              *   impacted by a change of this bound */
   int*                  startcomponents;    /**< array to store the components corresponding to startindices array */
   int                   nstarts;            /**< number of indices stored in startindices array */
   int                   startindicessize;   /**< size of the startindices and startcomponents arrays */
};

/** propagator data */
struct SCIP_PropData
{
   GENVBOUND**           genvboundstore;     /**< array to store genvbounds; fast access is provided by hashmaps
                                              *   lbgenvbounds and ubgenvbounds */
   SCIP_EVENTDATA**      lbevents;           /**< array of lower bound event data */
   SCIP_EVENTDATA**      ubevents;           /**< array of upper bound event data */
   SCIP_EVENTHDLR*       eventhdlr;          /**< genvbounds propagator event handler */
   SCIP_HASHMAP*         lbgenvbounds;       /**< hashmap to provide fast access to lower bound genvbounds in
                                              *   genvboundstore array */
   SCIP_HASHMAP*         ubgenvbounds;       /**< hashmap to provide fast access to upper bound genvbounds in
                                              *   genvboundstore array */
   SCIP_HASHMAP*         lbeventsmap;        /**< hashmap to provide fast access to lbevents array */
   SCIP_HASHMAP*         ubeventsmap;        /**< hashmap to provide fast access to ubevents array */
   SCIP_HASHMAP*         startmap;           /**< hashmap to provide fast access to startindices array */
   SCIP_PROP*            prop;               /**< pointer to genvbounds propagator */
   SCIP_NODE*            lastnodecaught;     /**< last node where events for starting indices were caught */
   SCIP_VAR*             cutoffboundvar;     /**< artificial variable representing primal cutoff bound */
   int*                  componentsstart;    /**< stores the components starting indices in genvboundstore array; the
                                              *   entry componentsstart[ncomponents] is equal to ngenvbounds, which
                                              *   makes it easier to iterate over all components */
   int                   componentsstartsize;/**< size of componentsstart array */
   int*                  startindices;       /**< storing indices of components where local propagation should start */
   int*                  startcomponents;    /**< components corresponding to indices stored in startindices array */
   int                   startindicessize;   /**< size of startindices and startcomponents arrays */
   int*                  gstartindices;      /**< storing indices of components where global propagation, i.e.,
                                              *   propagation of an improved primal bound, should start */
   int*                  gstartcomponents;   /**< components corresponding to indices stored in gstartindices array */
   int                   gstartindicessize;  /**< size of gstartindices and gstartcomponents arrays */
   SCIP_Real             lastcutoff;         /**< cutoff bound's value last time genvbounds propagator was called */
   int                   genvboundstoresize; /**< size of genvboundstore array */
   int                   ngenvbounds;        /**< number of genvbounds stored in genvboundstore array */
   int                   ncomponents;        /**< number of components in genvboundstore array */
   int                   nindices;           /**< number of indices stored in startindices array */
   int                   ngindices;          /**< number of indices stored in gstartindices array */
   int                   nlbevents;          /**< number of data entries in lbevents array */
   int                   nubevents;          /**< number of data entries in ubevents array */
   SCIP_Bool             issorted;           /**< stores wether array genvboundstore is topologically sorted */
   SCIP_Bool             global;             /**< apply global propagation? */
   SCIP_Bool             propinrootnode;     /**< apply genvbounds in root node if no new incumbent was found? */
   SCIP_Bool             sort;               /**< sort genvbounds and wait for bound change events? (otherwise all
                                              *   genvbounds are applied in each node) */
   SCIP_Bool             propasconss;        /**< should genvbounds be transformed to (linear) constraints? */
};


/*
 * Local methods
 */

/** returns correct cutoff bound value */
static
SCIP_Real getCutoffboundGenVBound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIPdebugMsg(scip, "cutoff = %.9g (%.9g + %.9g * %.9g)\n",
      (SCIPgetCutoffbound(scip) + SCIPgetTransObjoffset(scip)) * SCIPgetTransObjscale(scip),
      SCIPgetCutoffbound(scip), SCIPgetTransObjoffset(scip), SCIPgetTransObjscale(scip));

   /* the cutoff bound is valid w.r.t. the current objective function in the transformed problem; during presolving,
    * however, the objective function can change (e.g., when a variable is fixed, its contribution in the objective is
    * subtracted from the cutoff bound and added to the objective offset); we solve this by transforming the
    * contribution of the cutoff bound in the generalized variable bound to the original problem as described in
    * function SCIPgenVBoundAdd()
    */
   return (SCIPgetCutoffbound(scip) + SCIPgetTransObjoffset(scip)) * SCIPgetTransObjscale(scip);
}

/** returns corresponding genvbound in genvboundstore if there is one, NULL otherwise */
static
GENVBOUND* getGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the genvbounds propagator */
   SCIP_VAR*             var,                /**< bounds variable */
   SCIP_BOUNDTYPE        boundtype           /**< bounds type */
   )
{
   SCIP_HASHMAP* hashmap;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(var != NULL);

   hashmap = boundtype == SCIP_BOUNDTYPE_LOWER ? propdata->lbgenvbounds : propdata->ubgenvbounds;

   return (GENVBOUND*) SCIPhashmapGetImage(hashmap, var);
}

#ifdef SCIP_DEBUG
/** prints a genvbound as debug message */
static
void printGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   GENVBOUND*            genvbound           /**< genvbound to be printed */
   )
{
   SCIP_Bool first;
   int i;

   assert(genvbound != NULL);

   if( genvbound->boundtype == SCIP_BOUNDTYPE_UPPER )
   {
      SCIPdebugMsgPrint(scip, "- ");
   }

   SCIPdebugMsgPrint(scip, "<%s> >= ", SCIPvarGetName(genvbound->var));

   first = TRUE;
   for( i = 0; i < genvbound->ncoefs; i++ )
   {
      if( !first )
      {
         SCIPdebugMsgPrint(scip, " + ");
      }

      SCIPdebugMsgPrint(scip, "%g * <%s>", genvbound->coefs[i], SCIPvarGetName(genvbound->vars[i]));

      first = FALSE;
   }

   if( !SCIPisZero(scip, genvbound->cutoffcoef) )
   {
      SCIPdebugMsgPrint(scip, " + %g * cutoff_bound", genvbound->cutoffcoef);
   }

   if( !SCIPisZero(scip, genvbound->constant) )
   {
      SCIPdebugMsgPrint(scip, " + %g", genvbound->constant);
   }

   SCIPdebugMsgPrint(scip, "\n");
}
#endif

/** calculates the minactivity of a linear combination of variables stored in an array */
static
SCIP_Real getGenVBoundsMinActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            coefs,              /**< array of coefficients */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             global              /**< use global variable bounds? */
   )
{
   SCIP_Real minval;
   int i;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(coefs != NULL);
   assert(nvars >= 0);

   minval = 0.0;

   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real bound;

      /* get global or local bound */
      if( global )
         bound = coefs[i] > 0.0 ? SCIPvarGetLbGlobal(vars[i]) : SCIPvarGetUbGlobal(vars[i]);
      else
         bound = coefs[i] > 0.0 ? SCIPvarGetLbLocal(vars[i]) : SCIPvarGetUbLocal(vars[i]);

      /* with infinite bounds we cannot compute a valid minactivity and return minus infinity */
      if( SCIPisInfinity(scip, bound) || SCIPisInfinity(scip, -bound) )
         return -SCIPinfinity(scip);

      /* add contribution to minactivity */
      minval += coefs[i] * bound;
   }

   return minval;
}

/** calculates the minactivity of a linear combination of variables stored in the current conflict set */
static
SCIP_Real getGenVBoundsMinActivityConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            coefs,              /**< array of coefficients */
   int                   nvars,              /**< number of variables */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change at which minactivity should be computed; if NULL use local bounds */
   )
{
   SCIP_Real minval;
   int i;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(coefs != NULL);
   assert(nvars >= 0);

   minval = 0.0;

   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real bound;

      if( coefs[i] > 0.0 )
      {
         /* get bound at current bound change */
         bound = SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, TRUE);

         /* if bdchgidx is NULL, assert that we use local bounds */
         assert(bdchgidx != NULL || SCIPisEQ(scip, bound, SCIPvarGetLbLocal(vars[i])));

         /* if bdchgidx is not NULL, use the possibly tighter bound already in the current conflict set */
         if( bdchgidx != NULL && SCIPgetConflictVarLb(scip, vars[i]) > bound )
            bound = SCIPgetConflictVarLb(scip, vars[i]);
      }
      else
      {
         /* get bound at current bound change */
         bound = SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, TRUE);

         /* if bdchgidx is NULL, assert that we use local bounds */
         assert(bdchgidx != NULL || SCIPisEQ(scip, bound, SCIPvarGetUbLocal(vars[i])));

         /* if bdchgidx is not NULL, use the possibly tighter bound already in the current conflict set */
         if( bdchgidx != NULL && SCIPgetConflictVarUb(scip, vars[i]) < bound )
            bound = SCIPgetConflictVarUb(scip, vars[i]);
      }

      /* with infinite bounds we cannot compute a valid minactivity and return minus infinity */
      if( SCIPisInfinity(scip, bound) || SCIPisInfinity(scip, -bound) )
         return -SCIPinfinity(scip);

      /* add contribution to minactivity */
      minval += coefs[i] * bound;
   }

   return minval;
}

/** returns a valid bound given by a generalized variable bound */
static
SCIP_Real getGenVBoundsBound(
   SCIP*                 scip,               /**< SCIP data structure */
   GENVBOUND*            genvbound,          /**< generalized variable bound */
   SCIP_Bool             global              /**< use global variable bounds? */
   )
{
   SCIP_Real boundval;

   assert(scip != NULL);
   assert(genvbound != NULL);

   boundval = getGenVBoundsMinActivity(scip, genvbound->vars, genvbound->coefs, genvbound->ncoefs, global);

   if( SCIPisInfinity(scip, -boundval) )
      return (genvbound->boundtype == SCIP_BOUNDTYPE_LOWER) ? -SCIPinfinity(scip) : SCIPinfinity(scip);

   if( genvbound->cutoffcoef != 0.0 )
      boundval += genvbound->cutoffcoef * getCutoffboundGenVBound(scip);

   boundval += genvbound->constant;

   if( genvbound->boundtype == SCIP_BOUNDTYPE_UPPER )
      boundval *= -1.0;

   return boundval;
}

#ifdef WITH_DEBUG_SOLUTION
/** checks whether a generalized variable bound violates the debug solution */
static
SCIP_RETCODE checkDebugSolutionGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   GENVBOUND*            genvbound           /**< generalized variable bound */
   )
{
   SCIP_SOL* debugsol;
   SCIP_Real activity;
   SCIP_Real solval;
   int i;

   assert(scip != NULL);
   assert(genvbound != NULL);

   if( !SCIPdebugIsMainscip(scip) )
      return SCIP_OKAY;

   /* the genvbound must be valid for all cutoff bounds greater equal the objective value of the debug solution */
   SCIP_CALL( SCIPdebugGetSol(scip, &debugsol) );

   /* check whether a debug solution is available */
   if( debugsol == NULL )
      return SCIP_OKAY;

   activity = 0.0;
   for( i = 0; i < genvbound->ncoefs; i++ )
   {
      SCIP_CALL( SCIPdebugGetSolVal(scip, genvbound->vars[i], &solval) );
      if( solval != SCIP_UNKNOWN || solval != SCIP_INVALID )
         activity += genvbound->coefs[i] * solval;
      else
         printf("***** debug: ignoring variable with %s value in debug solution\n",
            solval == SCIP_UNKNOWN ? "unknown" : "invalid");
   }

   activity += genvbound->cutoffcoef *
      (SCIPgetSolTransObj(scip, debugsol) + SCIPgetTransObjoffset(scip)) * SCIPgetTransObjscale(scip);
   activity += genvbound->constant;

   SCIP_CALL( SCIPdebugGetSolVal(scip, genvbound->var, &solval) );
   if( solval != SCIP_UNKNOWN || solval != SCIP_INVALID )
   {
      if( genvbound->boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPdebugCheckLbGlobal(scip, genvbound->var, activity) );
      }
      else if( genvbound->boundtype == SCIP_BOUNDTYPE_UPPER )
      {
         SCIP_CALL( SCIPdebugCheckUbGlobal(scip, genvbound->var, -activity) );
      }
   }

   return SCIP_OKAY;
}
#endif

/** allocate local and global startindices, startcomponents and startmap */
static
SCIP_RETCODE createStartingData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   assert(propdata->startcomponents == NULL);
   assert(propdata->startindices == NULL);
   assert(propdata->startmap == NULL);
   assert(propdata->nindices == -1);

   assert(propdata->gstartindices == NULL);
   assert(propdata->gstartcomponents == NULL);
   assert(propdata->ngindices == -1);

   assert(propdata->ngenvbounds >= 1);
   assert(propdata->ncomponents >= 1);

   SCIPdebugMsg(scip, "create starting data\n");

   /* allocate memory for arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->startindices), propdata->ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->startcomponents), propdata->ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->gstartindices), propdata->ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->gstartcomponents), propdata->ncomponents) );
   propdata->startindicessize = propdata->ncomponents;
   propdata->gstartindicessize = propdata->ncomponents;

   /* create hashmap */
   SCIP_CALL( SCIPhashmapCreate(&(propdata->startmap), SCIPblkmem(scip), propdata->ncomponents) );

   propdata->nindices = 0;
   propdata->ngindices = 0;

   return SCIP_OKAY;
}

/** free local and global startindices, startcomponents and startmap */
static
SCIP_RETCODE freeStartingData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   SCIPdebugMsg(scip, "free starting data\n");

   if( propdata->startcomponents != NULL )
   {
      assert(propdata->startindices != NULL);
      assert(propdata->startmap != NULL);
      assert(propdata->nindices >= 0);

      SCIPfreeBlockMemoryArray(scip, &(propdata->startindices), propdata->startindicessize);
      SCIPfreeBlockMemoryArray(scip, &(propdata->startcomponents), propdata->startindicessize);
      propdata->startindicessize = 0;
      SCIPhashmapFree(&(propdata->startmap));
      propdata->nindices = -1;

      assert(propdata->gstartindices != NULL);
      assert(propdata->gstartcomponents != NULL);
      assert(propdata->ngindices >= 0);

      SCIPfreeBlockMemoryArray(scip, &(propdata->gstartindices), propdata->gstartindicessize);
      SCIPfreeBlockMemoryArray(scip, &(propdata->gstartcomponents), propdata->gstartindicessize);
      propdata->gstartindicessize = 0;
      propdata->ngindices = -1;
   }

   assert(propdata->startcomponents == NULL);
   assert(propdata->startindices == NULL);
   assert(propdata->startmap == NULL);
   assert(propdata->nindices == -1);

   assert(propdata->gstartindices == NULL);
   assert(propdata->gstartcomponents == NULL);
   assert(propdata->ngindices == -1);

   return SCIP_OKAY;
}

static
SCIP_RETCODE fillGlobalStartingData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);

   assert(propdata->gstartindices != NULL);
   assert(propdata->gstartcomponents != NULL);
   assert(propdata->ngindices == 0);

   SCIPdebugMsg(scip, "fill global starting data\n");

   for( i = 0; i < propdata->ncomponents; i++ )
   {
      int j;

      for( j = propdata->componentsstart[i]; j < propdata->componentsstart[i+1]; j++ ) /*lint !e679*/
      {
         assert(j < propdata->ngenvbounds);

         if( !SCIPisZero(scip, propdata->genvboundstore[j]->cutoffcoef) )
         {
            assert(SCIPisNegative(scip, propdata->genvboundstore[j]->cutoffcoef));

            propdata->gstartcomponents[propdata->ngindices] = i;
            propdata->gstartindices[propdata->ngindices] = j;

            /* go to next component */
            propdata->ngindices++;
            break;
         }
      }
   }

   /* resize arrays */
   if( propdata->gstartindicessize != propdata->ngindices )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(propdata->gstartindices), propdata->gstartindicessize, \
            propdata->ngindices) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(propdata->gstartcomponents), propdata->gstartindicessize, \
            propdata->ngindices) );
      propdata->gstartindicessize = propdata->ngindices;
   }

   return SCIP_OKAY;
}


/** resets local starting data */
static
SCIP_RETCODE resetLocalStartingData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->startcomponents != NULL);
   assert(propdata->startindices != NULL);
   assert(propdata->startmap != NULL);
   assert(propdata->nindices >= 0);

   SCIP_CALL( SCIPhashmapRemoveAll(propdata->startmap) );
   propdata->nindices = 0;

   return SCIP_OKAY;
}

/** frees sorted components data */
static
SCIP_RETCODE freeComponentsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   SCIPdebugMsg(scip, "free components data\n");

   if( propdata->componentsstart != NULL )
   {
      assert(propdata->ncomponents > 0);

      SCIPfreeBlockMemoryArray(scip, &(propdata->componentsstart), propdata->componentsstartsize);
      propdata->componentsstartsize = 0;
      propdata->ncomponents = -1;
   }

   assert(propdata->componentsstart == NULL);
   assert(propdata->ncomponents == -1);

   return SCIP_OKAY;
}

/** frees memory allocated for a generalized variable bound */
static
SCIP_RETCODE freeGenVBound(
   SCIP*                 scip,
   GENVBOUND*            genvbound
   )
{
   int i;

   assert(scip != NULL);
   assert(genvbound != NULL);
   assert(genvbound->coefs != NULL);
   assert(genvbound->vars != NULL);
   assert(genvbound->var != NULL);

   /* release variables */
   for( i = 0; i < genvbound->ncoefs; ++i )
   {
      assert(genvbound->vars[i] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &(genvbound->vars[i])) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &genvbound->var) );

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &(genvbound->coefs), genvbound->coefssize);
   SCIPfreeBlockMemoryArray(scip, &(genvbound->vars), genvbound->coefssize);
   SCIPfreeBlockMemory(scip, &genvbound);

   return SCIP_OKAY;
}

/** helper function to release all genvbounds */
static
SCIP_RETCODE freeGenVBounds(
   SCIP*                 scip,
   SCIP_PROPDATA*        propdata
   )
{
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);

   if( propdata->genvboundstore != NULL )
   {
      /* free genvbounds */
      for( i = propdata->ngenvbounds - 1; i >= 0; i-- )
      {
         SCIP_CALL( freeGenVBound(scip, propdata->genvboundstore[i]) );
      }

      /* free genvboundstore hashmaps */
      SCIPhashmapFree(&(propdata->lbgenvbounds));
      SCIPhashmapFree(&(propdata->ubgenvbounds));

      /* free genvboundstore array */
      SCIPfreeBlockMemoryArray(scip, &(propdata->genvboundstore), propdata->genvboundstoresize);

      /* set the number of genvbounds to zero */
      propdata->ngenvbounds = 0;

      /* free componentsstart array */
      SCIP_CALL( freeComponentsData(scip, propdata) );

      /* free starting indices data */
      SCIP_CALL( freeStartingData(scip, propdata) );

      /* release the cutoffboundvar and undo the locks */
      if( propdata->cutoffboundvar != NULL )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, propdata->cutoffboundvar, -1, -1) );
         SCIP_CALL( SCIPreleaseVar(scip, &(propdata->cutoffboundvar)) );
         propdata->cutoffboundvar = NULL;
         SCIPdebugMsg(scip, "release cutoffboundvar!\n");
      }
   }

   return SCIP_OKAY;
}

/** resolves propagation of lower bound on +/- left-hand side variable of a generalized variable bound */
static
SCIP_RETCODE resolveGenVBoundPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   GENVBOUND*            genvbound,          /**< genvbound data structure */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Real*            boundval,           /**< pointer to lower bound value on +/- left-hand side variable */
   SCIP_Bool*            success             /**< was the explanation succesful? */
   )
{
   SCIP_VAR* lhsvar;
   SCIP_VAR** vars;
   SCIP_Real minactivity;
   SCIP_Real tmpboundval;
   SCIP_Real slack;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(genvbound != NULL);
   assert(boundval != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* get left-hand side variable */
   lhsvar = genvbound->var;
   assert(lhsvar != NULL);

   /* get right-hand side variables */
   vars =  genvbound->vars;
   nvars = genvbound->ncoefs;
   assert(vars != NULL);

   /* if only the primal bound participates in the propagation, it is globally valid and should not be analyzed */
   assert(nvars > 0);

   /* when resolving a propagation, bdchgidx is not NULL and boundval should be the bound change performed for the
    * left-hand side variable
    */
   assert(bdchgidx == NULL || genvbound->boundtype != SCIP_BOUNDTYPE_LOWER || SCIPisEQ(scip,
         SCIPvarIsIntegral(genvbound->var) ? SCIPfeasCeil(scip, *boundval) : *boundval, SCIPgetVarLbAtIndex(scip, lhsvar, bdchgidx, TRUE)));
   assert(bdchgidx == NULL || genvbound->boundtype != SCIP_BOUNDTYPE_UPPER || SCIPisEQ(scip,
         SCIPvarIsIntegral(genvbound->var) ? SCIPfeasCeil(scip, *boundval) : *boundval, -SCIPgetVarUbAtIndex(scip, lhsvar, bdchgidx, TRUE)));

   /* when creating an initial conflict, bdchgidx is NULL and +/-boundval must exceed the upper/lower bound of the
    * left-hand side variable
    */
   assert(bdchgidx != NULL || genvbound->boundtype != SCIP_BOUNDTYPE_LOWER
      || SCIPisGT(scip, *boundval, SCIPvarGetUbLocal(lhsvar)));
   assert(bdchgidx != NULL || genvbound->boundtype != SCIP_BOUNDTYPE_UPPER
      || SCIPisGT(scip, *boundval, -SCIPvarGetLbLocal(lhsvar)));

   SCIPdebugMsg(scip, "resolving genvbound propagation: lhs=%s<%s> >= boundval=%.15g\n",
      genvbound->boundtype == SCIP_BOUNDTYPE_LOWER ? "+" : "-", SCIPvarGetName(lhsvar), *boundval);

   /* subtract constant terms from bound value */
   tmpboundval = *boundval;
   tmpboundval -= genvbound->cutoffcoef * getCutoffboundGenVBound(scip);
   tmpboundval -= genvbound->constant;

   SCIPdebugMsg(scip, "subtracting constant terms gives boundval=%.15g\n", tmpboundval);

   /* compute minimal activity; if bdchgidx is NULL, we create the initial conflict and use local bounds */
   minactivity = getGenVBoundsMinActivityConflict(scip, genvbound->vars, genvbound->coefs, genvbound->ncoefs, bdchgidx);

   SCIPdebugMsg(scip, "minactivity of right-hand side is minactivity=%.15g\n", minactivity);

   /* a genvbound might have been replaced since the propagation took place, hence we have to check that the current
    * genvbound can explain the propagation at the given bound change index; note that by now, with smaller cutoff
    * bound, we might even perform a stronger propagation
    */
   if( SCIPisLT(scip, minactivity, tmpboundval) )
   {
      SCIPdebugMsg(scip, "minactivity is too small to explain propagation; was genvbound replaced?\n");
      return SCIP_OKAY;
   }

   /* if bdchgidx is NULL, i.e., we create the initial conflict, we should be able to explain the bound change */
   assert(SCIPisGE(scip, minactivity, tmpboundval));

   slack = MAX(minactivity - tmpboundval, 0.0);

   SCIPdebugMsg(scip, "slack=%.15g\n", slack);

   /* add variables on the right-hand side as reasons for propagation */
   for( i = 0; i < nvars; i++ )
   {
      assert(vars[i] != NULL);
      assert(SCIPisEQ(scip, SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, TRUE), SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, FALSE)));
      assert(SCIPisEQ(scip, SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, TRUE), SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, FALSE)));

      /* coefficient is positive */
      if( genvbound->coefs[i] > 0.0 )
      {
         SCIP_Real lbatindex;
         SCIP_Real conflictlb;

         /* get bound at current bound change */
         lbatindex = SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, TRUE);

         /* get bound already enforced by conflict set */
         conflictlb = SCIPgetConflictVarLb(scip, genvbound->vars[i]);
         assert(SCIPisGE(scip, conflictlb, SCIPvarGetLbGlobal(genvbound->vars[i])));

         SCIPdebugMsg(scip, "lower bound of variable <%s> (genvbound->vars[%d]) in conflict set / at index is %.15g / %.15g\n",
            SCIPvarGetName(genvbound->vars[i]), i, conflictlb, lbatindex);

         /* if bound is already enforced by conflict set we do not need to add the bound change; since we used the
          * tighest bound already when computing the initial minactivity, the slack is already correct
          */
         if( SCIPisLE(scip, lbatindex, conflictlb) )
         {
            SCIPdebugMsg(scip, "skipping lower bound of variable <%s> (genvbound->vars[%d]) already enforced in conflict set\n",
               SCIPvarGetName(genvbound->vars[i]), i);
         }
         else
         {
            SCIP_Real relaxedlb;

            /* compute relaxed bound that would suffice to explain the bound change */
            relaxedlb = lbatindex - (slack / genvbound->coefs[i]);
            assert(relaxedlb <= lbatindex);

            /* add variable to conflict set */
            SCIP_CALL( SCIPaddConflictRelaxedLb(scip, genvbound->vars[i], bdchgidx, relaxedlb ) );

            /* get new bound of variable in conflict set; after possible bound widening in SCIPaddConflictLbRelaxed(),
             * it should be between conflictlb and lbatindex
             */
            relaxedlb = SCIPgetConflictVarLb(scip, genvbound->vars[i]);
            assert(SCIPisGE(scip, relaxedlb, conflictlb));
            assert(SCIPisLE(scip, relaxedlb, lbatindex));

            /* update slack and ensure that its nonegative */
            slack -= genvbound->coefs[i] * (lbatindex - relaxedlb);
            slack = MAX(slack, 0.0);

            SCIPdebugMsg(scip, "added lower bound of variable <%s> (genvbound->vars[%d]); new slack=%.15g\n",
               SCIPvarGetName(genvbound->vars[i]), i, slack);
         }
      }
      /* coefficient is negative */
      else
      {
         SCIP_Real ubatindex;
         SCIP_Real conflictub;

         /* get bound at current bound change */
         ubatindex = SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, TRUE);

         /* get bound already enforced by conflict set */
         conflictub = SCIPgetConflictVarUb(scip, genvbound->vars[i]);
         assert(SCIPisLE(scip, conflictub, SCIPvarGetUbGlobal(genvbound->vars[i])));

         SCIPdebugMsg(scip, "upper bound of variable <%s> (genvbound->vars[%d]) in conflict set / at index is %.15g / %.15g\n",
            SCIPvarGetName(genvbound->vars[i]), i, conflictub, ubatindex);

         /* if bound is already enforced by conflict set we do not need to add the bound change; since we used the
          * tighest bound already when computing the initial minactivity, the slack is already correct
          */
         if( SCIPisGE(scip, ubatindex, conflictub) )
         {
            SCIPdebugMsg(scip, "skipping upper bound of variable <%s> (genvbound->vars[%d]) already enforced in conflict set\n",
               SCIPvarGetName(genvbound->vars[i]), i);
         }
         else
         {
            SCIP_Real relaxedub;

            /* compute relaxed bound that would suffice to explain the bound change */
            relaxedub = ubatindex - (slack / genvbound->coefs[i]);
            assert(relaxedub >= ubatindex);

            /* add variable to conflict set */
            SCIP_CALL( SCIPaddConflictRelaxedUb(scip, genvbound->vars[i], bdchgidx, relaxedub ) );

            /* get new bound of variable in conflict set; after possible bound widening in SCIPaddConflictUbRelaxed(),
             * it should be between conflictub and ubatindex
             */
            relaxedub = SCIPgetConflictVarUb(scip, genvbound->vars[i]);
            assert(SCIPisLE(scip, relaxedub, conflictub));
            assert(SCIPisGE(scip, relaxedub, ubatindex));

            /* update slack and ensure that its nonegative */
            slack -= genvbound->coefs[i] * (ubatindex - relaxedub);
            slack = MAX(slack, 0.0);

            SCIPdebugMsg(scip, "added upper bound of variable <%s> (genvbound->vars[%d]); new slack=%.15g\n",
               SCIPvarGetName(genvbound->vars[i]), i, slack);
         }
      }
   }

   /* if slack is positive, return increased boundval */
   if( SCIPisPositive(scip, slack) )
      tmpboundval += slack;

   /* add constant terms again */
   tmpboundval += genvbound->cutoffcoef * getCutoffboundGenVBound(scip);
   tmpboundval += genvbound->constant;

   /* boundval should not have been decreased; if this happened nevertheless, maybe due to numerical errors, we quit
    * without success
    */
   if( SCIPisLT(scip, tmpboundval, *boundval) )
   {
      SCIPdebugMsg(scip, "boundval was reduced from %.15g to %.15g; propagation not resolved\n", *boundval, tmpboundval);
      return SCIP_OKAY;
   }

   /* return widened boundval */
   *boundval = tmpboundval;
   *success = TRUE;

   return SCIP_OKAY;
}

/** create initial conflict */
static
SCIP_RETCODE analyzeGenVBoundConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   GENVBOUND*            genvbound           /**< genvbound data structure */
   )
{
   SCIP_Bool success;

   assert(scip != NULL);
   assert(genvbound != NULL);

   /* check if conflict analysis is applicable */
   if( !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   /* initialize conflict analysis */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, !SCIPisInfinity(scip, REALABS(SCIPgetCutoffbound(scip)))) );

   /* left-hand side variable >= ... */
   if( genvbound->boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_Real infeasthreshold;
      SCIP_Real bound;

      /* get minimal right-hand side bound that leads to infeasibility; first try with a factor of 2 for robustness */
      bound = REALABS(SCIPvarGetUbLocal(genvbound->var));
      infeasthreshold = MAX(bound, 1.0) * 2 * SCIPfeastol(scip);
      bound = SCIPvarGetUbLocal(genvbound->var) + infeasthreshold;

      /* add right-hand side variables that force the lower bound of the left-hand side variable above its upper bound
       * to conflict set
       */
      SCIP_CALL( resolveGenVBoundPropagation(scip, genvbound, NULL, &bound, &success) );
      assert(!success || SCIPisFeasGT(scip, bound, SCIPvarGetUbLocal(genvbound->var)));

      /* if infeasibility cannot be proven with the tighter bound, try with actual bound */
      if( !success )
      {
         bound = REALABS(SCIPvarGetUbLocal(genvbound->var));
         infeasthreshold = MAX(bound, 1.0) * SCIPfeastol(scip);
         bound = SCIPvarGetUbLocal(genvbound->var) + infeasthreshold;

         SCIP_CALL( resolveGenVBoundPropagation(scip, genvbound, NULL, &bound, &success) );
         success = success && SCIPisFeasGT(scip, bound, SCIPvarGetUbLocal(genvbound->var));
      }

      /* compute upper bound on left-hand side variable that leads to infeasibility */
      bound -= infeasthreshold;
      success = success && SCIPisGE(scip, bound, SCIPvarGetUbLocal(genvbound->var));

      /* initial reason could not be constructed, maybe due to numerics; do not apply conflict analysis */
      if( !success )
      {
         SCIPdebugMsg(scip, "strange: could not create initial reason to start conflict analysis\n");
         return SCIP_OKAY;
      }

      /* if bound is already enforced by conflict set we do not have to add it */
      if( SCIPisGE(scip, bound, SCIPgetConflictVarUb(scip, genvbound->var)) )
      {
         SCIPdebugMsg(scip, "skipping upper bound of left-hand side variable <%s> already enforced in conflict set\n",
            SCIPvarGetName(genvbound->var));
      }
      else
      {
         SCIPdebugMsg(scip, "adding upper bound of left-hand side variable <%s>\n", SCIPvarGetName(genvbound->var));

         SCIP_CALL( SCIPaddConflictRelaxedUb(scip, genvbound->var, NULL, bound) );
      }
   }
   /* left-hand side variable <= ..., i.e., - left-hand side variable >= ... */
   else
   {
      SCIP_Real infeasthreshold;
      SCIP_Real bound;

      /* get minimal right-hand side bound that leads to infeasibility; try with a factor of 2 first for robustness */
      bound = REALABS(SCIPvarGetLbLocal(genvbound->var));
      infeasthreshold = MAX(bound, 1.0) * 2 * SCIPfeastol(scip);
      bound = -SCIPvarGetLbLocal(genvbound->var) + infeasthreshold;

      /* add right-hand side variables that force the upper bound of the left-hand side variable below its lower bound
       * to conflict set
       */
      SCIP_CALL( resolveGenVBoundPropagation(scip, genvbound, NULL, &bound, &success) );
      assert(!success || SCIPisFeasLT(scip, -bound, SCIPvarGetLbLocal(genvbound->var)));

      /* if infeasibility cannot be proven with the tighter bound, try with actual bound */
      if( !success )
      {
         bound = REALABS(SCIPvarGetLbLocal(genvbound->var));
         infeasthreshold = MAX(bound, 1.0) * SCIPfeastol(scip);
         bound = -SCIPvarGetLbLocal(genvbound->var) + infeasthreshold;

         SCIP_CALL( resolveGenVBoundPropagation(scip, genvbound, NULL, &bound, &success) );
         success = success && SCIPisFeasLT(scip, -bound, SCIPvarGetLbLocal(genvbound->var));
      }

      /* compute lower bound on left-hand side variable that leads to infeasibility */
      bound = -bound + infeasthreshold;
      success = success && SCIPisLE(scip, bound, SCIPvarGetLbLocal(genvbound->var));

      /* initial reason could not be constructed, maybe due to numerics; do not apply conflict analysis */
      if( !success )
      {
         SCIPdebugMsg(scip, "strange: could not create initial reason to start conflict analysis\n");
         return SCIP_OKAY;
      }

      /* if bound is already enforced by conflict set we do not have to add it */
      if( SCIPisLE(scip, bound, SCIPgetConflictVarLb(scip, genvbound->var)) )
      {
         SCIPdebugMsg(scip, "skipping lower bound of left-hand side variable <%s> already enforced in conflict set\n",
            SCIPvarGetName(genvbound->var));
      }
      else
      {
         SCIPdebugMsg(scip, "adding lower bound of left-hand side variable <%s>\n", SCIPvarGetName(genvbound->var));

         SCIP_CALL( SCIPaddConflictRelaxedLb(scip, genvbound->var, NULL, bound) );
      }
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );

   return SCIP_OKAY;
}

/** apply propagation for one generalized variable bound; also if the left-hand side variable is locally fixed, we
 *  compute the right-hand side minactivity to possibly detect infeasibility
 */
static
SCIP_RETCODE applyGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< genvbounds propagator */
   GENVBOUND*            genvbound,          /**< genvbound data structure */
   SCIP_Bool             global,             /**< apply global bound changes? (global: true, local: false)*/
   SCIP_RESULT*          result,             /**< result pointer */
   int*                  nchgbds             /**< counter to increment if bound was tightened */
   )
{
   SCIP_Real boundval;
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(genvbound != NULL);
   assert(genvbound->var != NULL);
   assert(SCIPvarGetStatus(genvbound->var) != SCIP_VARSTATUS_MULTAGGR);
   assert(result != NULL);
   assert(*result != SCIP_DIDNOTRUN);

   /* get bound value provided by genvbound */
   boundval = getGenVBoundsBound(scip, genvbound, global);

   if( SCIPisInfinity(scip, REALABS(boundval)) )
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   {
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real new_lb;
      SCIP_Real new_ub;

      lb = global ? SCIPvarGetLbGlobal(genvbound->var) : SCIPvarGetLbLocal(genvbound->var);
      ub = global ? SCIPvarGetUbGlobal(genvbound->var) : SCIPvarGetUbLocal(genvbound->var);
      new_lb = genvbound->boundtype == SCIP_BOUNDTYPE_LOWER ? boundval : lb;
      new_ub = genvbound->boundtype == SCIP_BOUNDTYPE_UPPER ? boundval : ub;

      SCIPdebugMsg(scip, "  %s genvbound propagation for <%s>\n", global ?
         "global" : "local", SCIPvarGetName(genvbound->var));
      SCIPdebugMsg(scip, "  genvbound: ");
      printGenVBound(scip, genvbound);
      SCIPdebugMsg(scip, "    [%.15g,%.15g] -> [%.15g,%.15g]\n", lb, ub, new_lb, new_ub);
   }
#endif

   /* tighten bound globally */
   if( global || genvbound->ncoefs <= 0 )
   {
      if( genvbound->boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, genvbound->var, boundval, FALSE, &infeas, &tightened) );
      }
      else
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, genvbound->var, boundval, FALSE, &infeas, &tightened) );
      }
   }
   /* tighten bound locally and start conflict analysis in case of infeasibility; as inferinfo we pass the index of the
    * genvbound that was used for propagation
    */
   else
   {
      if( genvbound->boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPinferVarLbProp(scip, genvbound->var, boundval, prop, genvbound->index, FALSE, &infeas, &tightened) );

         /* initialize conflict analysis if infeasible */
         if( infeas )
         {
            SCIPdebugMsg(scip, " -> lower bound tightening on variable <%s> led to infeasibility\n",
                  SCIPvarGetName(genvbound->var));

            SCIP_CALL( analyzeGenVBoundConflict(scip, genvbound) );
         }
      }
      else
      {
         SCIP_CALL( SCIPinferVarUbProp(scip, genvbound->var, boundval, prop, genvbound->index, FALSE, &infeas, &tightened) );

         /* initialize conflict analysis if infeasible */
         if( infeas )
         {
            SCIPdebugMsg(scip, " -> upper bound tightening on variable <%s> led to infeasibility\n",
                  SCIPvarGetName(genvbound->var));

            SCIP_CALL( analyzeGenVBoundConflict(scip, genvbound) );
         }
      }
   }

   /* handle result */
   if( infeas )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMsg(scip, "    cutoff!\n");
   }
   else if( tightened )
   {
      *result = SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING ? SCIP_SUCCESS : SCIP_REDUCEDDOM;
      if( nchgbds != NULL )
         ++(*nchgbds);
      SCIPdebugMsg(scip, "    tightened!\n");
   }

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** prints event data as debug message */
static
void printEventData(
   SCIP_EVENTDATA*       eventdata,
   SCIP_BOUNDTYPE        boundtype
   )
{
   int i;

   SCIPdebugMessage("event data: %s bound of <%s> tightened ==> start propagating at ",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(eventdata->var));

   /* if there is eventdata it should contain at least one starting index */
   assert(eventdata->nstarts > 0);

   for( i = 0; i < eventdata->nstarts; i++ )
   {
      SCIPdebugPrintf("(component %d, index %d) ", eventdata->startcomponents[i], eventdata->startindices[i]);
   }
   SCIPdebugPrintf("\n");
}
#endif

/** frees event data */
static
SCIP_RETCODE freeEventData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTDATA**      eventdata           /**< event data to be freed */
   )
{
   assert(scip != NULL);
   assert(eventdata != NULL);
   assert(*eventdata != NULL);

   SCIPfreeBlockMemoryArray(scip, &((*eventdata)->startcomponents), (*eventdata)->startindicessize);
   SCIPfreeBlockMemoryArray(scip, &((*eventdata)->startindices), (*eventdata)->startindicessize);

   (*eventdata)->startindicessize = 0;
   (*eventdata)->nstarts = -1;
   (*eventdata)->var = NULL;
   (*eventdata)->prop = NULL;

   SCIPfreeBlockMemory(scip, eventdata);

   return SCIP_OKAY;
}

/** frees all eventdata stored */
static
SCIP_RETCODE freeAllEventData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);

   if( propdata->lbevents != NULL )
   {
      assert(propdata->ubevents != NULL);
      assert(propdata->lbeventsmap != NULL);
      assert(propdata->ubeventsmap != NULL);

      SCIPhashmapFree(&(propdata->lbeventsmap));
      SCIPhashmapFree(&(propdata->ubeventsmap));

      for( i = propdata->nlbevents - 1; i >= 0; i-- )
      {
         SCIP_CALL( freeEventData(scip, &(propdata->lbevents[i])) );
      }

      for( i = propdata->nubevents - 1; i >= 0; i-- )
      {
         SCIP_CALL( freeEventData(scip, &(propdata->ubevents[i])) );
      }

      SCIPfreeBlockMemoryArray(scip, &(propdata->ubevents), propdata->nubevents);
      SCIPfreeBlockMemoryArray(scip, &(propdata->lbevents), propdata->nlbevents);
      propdata->nlbevents = -1;
      propdata->nubevents = -1;
   }

   assert(propdata->lbevents == NULL);
   assert(propdata->ubevents == NULL);
   assert(propdata->lbeventsmap == NULL);
   assert(propdata->ubeventsmap == NULL);
   assert(propdata->nlbevents == -1);
   assert(propdata->nubevents == -1);

   return SCIP_OKAY;
}

/** drops all events caught by genvbounds propagator and frees their data */
static
SCIP_RETCODE dropAndFreeEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   int i;

   SCIPdebugMsg(scip, "drop and free events\n");

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->eventhdlr != NULL);

   if( propdata->lbevents != NULL )
   {
      assert(propdata->ubevents != NULL);
      assert(propdata->nlbevents >= 0);
      assert(propdata->nubevents >= 0);

      for( i = propdata->nlbevents - 1; i >= 0; i-- )
      {
         /* drop event */
         SCIP_CALL( SCIPdropVarEvent(scip, propdata->lbevents[i]->var, SCIP_EVENTTYPE_LBTIGHTENED, propdata->eventhdlr,
               propdata->lbevents[i], -1) );
      }

      for( i = propdata->nubevents - 1; i >= 0; i-- )
      {
         /* drop event */
         SCIP_CALL( SCIPdropVarEvent(scip, propdata->ubevents[i]->var, SCIP_EVENTTYPE_UBTIGHTENED, propdata->eventhdlr,
               propdata->ubevents[i], -1) );
      }

      /* free event data */
      SCIP_CALL( freeAllEventData(scip, propdata) );
   }

   assert(propdata->lbevents == NULL);
   assert(propdata->ubevents == NULL);
   assert(propdata->nlbevents == -1);
   assert(propdata->nubevents == -1);

   return SCIP_OKAY;
}

/** returns the corresponding event data entry in the corresponding array, if there is one; if not: allocates a new
 *  event data entry, stores it in the array and returns its adress
 */
static
SCIP_RETCODE getEventData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the genvbounds propagator */
   SCIP_VAR*             var,                /**< variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound */
   SCIP_EVENTDATA**      eventdata           /**< event data to return */
   )
{
   SCIP_HASHMAP* hashmap;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(var != NULL);

   hashmap = boundtype == SCIP_BOUNDTYPE_LOWER ? propdata->lbeventsmap : propdata->ubeventsmap;

   if( SCIPhashmapExists(hashmap, var) )
      *eventdata = (SCIP_EVENTDATA*) SCIPhashmapGetImage(hashmap, var);
   else
   {
      /* set up new eventdata entry */
      SCIP_CALL( SCIPallocBlockMemory(scip, eventdata) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*eventdata)->startcomponents), propdata->ncomponents) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*eventdata)->startindices), propdata->ncomponents) );
      (*eventdata)->startindicessize = propdata->ncomponents;
      (*eventdata)->nstarts = 0;
      (*eventdata)->var = var;
      (*eventdata)->prop = propdata->prop;

      /* store event data in eventarray */
      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         propdata->lbevents[propdata->nlbevents] = *eventdata;
         propdata->nlbevents++;
      }
      else
      {
         propdata->ubevents[propdata->nubevents] = *eventdata;
         propdata->nubevents++;
      }

      /* store hashmap entry */
      SCIP_CALL( SCIPhashmapInsert(hashmap, var, (*eventdata)) );
   }

   return SCIP_OKAY;
}

/** adds an event to the event array lbevents (if boundtype == SCIP_BOUNDTYPE_LOWER) or ubevents (if boundtype ==
 *  SCIP_BOUNDTYPE_UPPER)
 */
static
SCIP_RETCODE addEventData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the genvbounds propagator */
   SCIP_VAR*             var,                /**< variable thats event to be added */
   int                   startindex,         /**< starting index */
   int                   startcomponent,     /**< starting components index */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound */
   )
{
   SCIP_EVENTDATA* eventdata;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(var != NULL);
   assert(startindex >= 0);
   assert(startcomponent >= 0);

   /* get eventdata entry */
   SCIP_CALL( getEventData(scip, propdata, var, boundtype, &eventdata) );
   assert(eventdata != NULL);

   if( eventdata->nstarts > 0 && eventdata->startcomponents[eventdata->nstarts - 1] == startcomponent )
   {
      /* if there is already a starting index for startcomponent stored at the last entry of eventdata->startindices,
       * it should be smaller; this relies on the implementation of setUpEvents(), calling addEventData() in
       * topological order
       */
      assert(eventdata->startindices[eventdata->nstarts - 1] < startindex);
   }
   else
   {
      /* append starting information */
      eventdata->startcomponents[eventdata->nstarts] = startcomponent;
      eventdata->startindices[eventdata->nstarts] = startindex;

      /* increase counter */
      eventdata->nstarts++;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE setUpEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   int nprobvars;
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->eventhdlr != NULL);
   assert(propdata->lbevents == NULL);
   assert(propdata->ubevents == NULL);
   assert(propdata->issorted);
   assert(propdata->nlbevents == -1);
   assert(propdata->nubevents == -1);

   SCIPdebugMsg(scip, "set up events\n");

   /* allocate lbevents, ubevents, and their hashmaps */
   nprobvars = SCIPgetNVars(scip) + SCIPgetNFixedVars(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->lbevents), nprobvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->ubevents), nprobvars) );
   SCIP_CALL( SCIPhashmapCreate(&(propdata->lbeventsmap), SCIPblkmem(scip), nprobvars) );
   SCIP_CALL( SCIPhashmapCreate(&(propdata->ubeventsmap), SCIPblkmem(scip), nprobvars) );
   propdata->nlbevents = 0;
   propdata->nubevents = 0;

   /* loop over all components of genvboundstore */
   for( i = 0; i < propdata->ncomponents; i++ )
   {
      int j;

      /* loop over all genvbounds in this component */
      for( j = propdata->componentsstart[i]; j < propdata->componentsstart[i+1]; j++ ) /*lint !e679*/
      {
         GENVBOUND* genvbound;
         int k;

         assert(j < propdata->ngenvbounds);

         genvbound = propdata->genvboundstore[j];
         assert(genvbound != NULL);

         /* loop over all coefficients in this genvbound */
         for( k = 0; k < genvbound->ncoefs; k++ )
         {
            if( SCIPisPositive(scip, genvbound->coefs[k]) )
            {
               SCIP_CALL( addEventData(scip, propdata, genvbound->vars[k], j, i, SCIP_BOUNDTYPE_LOWER) );
            }
            else
            {
               SCIP_CALL( addEventData(scip, propdata, genvbound->vars[k], j, i, SCIP_BOUNDTYPE_UPPER) );
            }
         }
      }
   }

   /* resize lbevents and ubevents array */
   assert(propdata->nlbevents <= nprobvars);
   assert(propdata->nubevents <= nprobvars);
   if( propdata->nlbevents < nprobvars )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(propdata->lbevents), nprobvars, propdata->nlbevents) );
   }
   if( propdata->nubevents < nprobvars )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(propdata->ubevents), nprobvars, propdata->nubevents) );
   }

   /* resize and register lower bound events */
   for( i = 0; i < propdata->nlbevents; i++ )
   {
      SCIP_EVENTDATA* eventdata = propdata->lbevents[i];

      assert(eventdata != NULL);
      assert(eventdata->nstarts > 0);
      assert(eventdata->startcomponents != NULL);
      assert(eventdata->startindices != NULL);

      /* resize arrays stored in eventdata */
      if( eventdata->startindicessize != eventdata->nstarts )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(eventdata->startcomponents), eventdata->startindicessize, \
               eventdata->nstarts) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(eventdata->startindices), eventdata->startindicessize, \
               eventdata->nstarts) );
         eventdata->startindicessize = eventdata->nstarts;
      }

      /* register event */
      SCIP_CALL( SCIPcatchVarEvent(scip, eventdata->var, SCIP_EVENTTYPE_LBTIGHTENED, propdata->eventhdlr, eventdata, \
            NULL) );
   }

   /* resize and register upper bound events */
   for( i = 0; i < propdata->nubevents; i++ )
   {
      SCIP_EVENTDATA* eventdata = propdata->ubevents[i];

      assert(eventdata != NULL);
      assert(eventdata->nstarts > 0);
      assert(eventdata->startcomponents != NULL);
      assert(eventdata->startindices != NULL);

      /* resize arrays stored in eventdata */
      if( eventdata->startindicessize != eventdata->nstarts )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(eventdata->startcomponents), eventdata->startindicessize, \
               eventdata->nstarts) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(eventdata->startindices), eventdata->startindicessize, \
               eventdata->nstarts) );
         eventdata->startindicessize = eventdata->nstarts;
      }
      /* register event */
      SCIP_CALL( SCIPcatchVarEvent(scip, eventdata->var, SCIP_EVENTTYPE_UBTIGHTENED, propdata->eventhdlr, eventdata,
            NULL) );
   }

   return SCIP_OKAY;
}

/** performs a topological sort on genvboundstore array
 *
 *  The genvbounds graph is defined as follows: Given two genvbounds
 *
 *    (genvbound1)      c1 * x_i1 >= RHS1
 *
 *  and
 *
 *    (genvbound2)      c2 * x_i2 >= RHS2,
 *
 *  there is an arc from genvbound1 to genvbound2 iff c1 = +1 and x_i1 appears with positive coefficient in RHS2 or
 *  c1 = -1 and x_i1 appears with negative coefficient in RHS2; in this case, a bound change of x_i1 deduced from
 *  genvbound1 improves genvbound2's minactivity in RHS2.
 *
 *  The method computes the strongly connected components and sorts them topologically. The order of the nodes in an
 *  strongly connected component is arbitrary.
 */
static
SCIP_RETCODE sortGenVBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   GENVBOUND** genvboundssorted;            /* array to store the sorted genvbounds */
   SCIP_DIGRAPH* graph;
   int* strongcomponents;
   int* strongcompstartidx;
   int sortedindex;
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->componentsstart == NULL);

   SCIPdebugMsg(scip, "(re-)sort genvbounds topologically\n");

   /* create digraph */
   SCIP_CALL( SCIPcreateDigraph(scip, &graph, propdata->ngenvbounds) );

   /* add outgoing arcs for each genvbound */
   for( i = 0; i < propdata->ngenvbounds; i++ )
   {
      GENVBOUND* genvbound;
      int j;

      assert(i < propdata->ngenvbounds);

      genvbound = propdata->genvboundstore[i];

      for( j = 0; j < genvbound->ncoefs; j++ )
      {
         if( SCIPisPositive(scip, genvbound->coefs[j]) &&
            SCIPhashmapExists(propdata->lbgenvbounds, genvbound->vars[j]) )
         {
            int from = ((GENVBOUND*) SCIPhashmapGetImage(propdata->lbgenvbounds, genvbound->vars[j]))->index;
            SCIP_CALL( SCIPdigraphAddArc(graph, from, i, NULL) );
         }
         else if( SCIPisNegative(scip, genvbound->coefs[j]) &&
            SCIPhashmapExists(propdata->ubgenvbounds, genvbound->vars[j]) )
         {
            int from = ((GENVBOUND*) SCIPhashmapGetImage(propdata->ubgenvbounds, genvbound->vars[j]))->index;
            SCIP_CALL( SCIPdigraphAddArc(graph, from, i, NULL) );
         }
      }
   }

   /* perform the topological sort */
   SCIP_CALL( SCIPdigraphComputeUndirectedComponents(graph, 1, NULL, &(propdata->ncomponents)) );
   SCIP_CALL( SCIPdigraphTopoSortComponents(graph) );
   assert(SCIPdigraphGetNComponents(graph) == propdata->ncomponents);

   /* allocate memory for genvboundssorted and componentsstart array */
   SCIP_CALL( SCIPallocBufferArray(scip, &genvboundssorted, propdata->ngenvbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->componentsstart), propdata->ncomponents + 1) );
   propdata->componentsstartsize = propdata->ncomponents + 1;

   /* allocate memory for strong component arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &strongcomponents, SCIPdigraphGetNNodes(graph)) ); /*lint !e666*/
   SCIP_CALL( SCIPallocBufferArray(scip, &strongcompstartidx, SCIPdigraphGetNNodes(graph) + 1) ); /*lint !e666*/

   /* compute sorted genvbounds array, fill componentsstart array */
   sortedindex = 0;
   propdata->componentsstart[propdata->ncomponents] = propdata->ngenvbounds;
   for( i = 0; i < propdata->ncomponents; i++ )
   {
      int j;
      int *nodes;
      int nnodes;
      int nstrongcomponents;

      SCIPdigraphGetComponent(graph, i, &nodes, &nnodes);
      propdata->componentsstart[i] = sortedindex;

      /* compute the strong components of the i-th undirected component */
      if( nnodes > 2 )
      {
         SCIP_CALL( SCIPdigraphComputeDirectedComponents(graph, i, strongcomponents, strongcompstartidx,
               &nstrongcomponents) );

         for( j = 0; j < nnodes; ++j )
         {
            int node;

            /* take the nodes at the end of the strong components array first to respect the topological
             * order of the different strong components
             */
            node = strongcomponents[nnodes - j - 1];

            assert(node < propdata->ngenvbounds);
            genvboundssorted[sortedindex] = propdata->genvboundstore[node];
            sortedindex++;
         }
      }
      else
      {
         for( j = 0; j < nnodes; j++ )
         {
            assert(nodes[j] < propdata->ngenvbounds);
            genvboundssorted[sortedindex] = propdata->genvboundstore[nodes[j]];
            sortedindex++;
         }
      }
   }
   assert(sortedindex == propdata->ngenvbounds);

   /* copy sorted genvbounds into genvboundstore */
   for( i = 0; i < propdata->ngenvbounds; i++ )
   {
      assert(genvboundssorted[i] != NULL);

      propdata->genvboundstore[i] = genvboundssorted[i];
      propdata->genvboundstore[i]->index = i;
   }
   SCIPfreeBufferArray(scip, &(genvboundssorted));

   /* free strong component arrays */
   SCIPfreeBufferArray(scip, &strongcompstartidx);
   SCIPfreeBufferArray(scip, &strongcomponents);

   /* free digraph */
   SCIPdigraphFree(&graph);

   /* remember genvboundstore as sorted */
   propdata->issorted = TRUE;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "genvbounds got: %d\n", propdata->ngenvbounds);
   for( i = 0; i < propdata->ncomponents; i++ )
   {
      int j;

      SCIPdebugMsg(scip, "{\n");

      for( j = propdata->componentsstart[i]; j < propdata->componentsstart[i+1]; j++ )
      {
         SCIPdebugMsg(scip, "  [%d] ", j);
         printGenVBound(scip, propdata->genvboundstore[j]);
      }

      SCIPdebugMsg(scip, "}\n");
   }
#endif

   return SCIP_OKAY;
}

/** apply propagation of generalized variable bounds */
static
SCIP_RETCODE applyGenVBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< genvbounds propagator */
   SCIP_Bool             global,             /**< use global variable bounds for propagation? */
   SCIP_RESULT*          result,             /**< result pointer */
   int*                  nchgbds             /**< counter to increase by the number of changed bounds */
   )
{
   SCIP_PROPDATA* propdata;
   int* startingcomponents;
   int* startingindices;
   int nindices;
   int i;

   SCIPdebugMsg(scip, "applying %s genvbound propagation in depth %d\n", global ?
      "global" : "local", SCIPgetDepth(scip));

   assert(scip != NULL);
   assert(prop != NULL);
   assert(result != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->genvboundstore != NULL);

   if( *result == SCIP_DIDNOTRUN )
      *result = SCIP_DIDNOTFIND;

   /* if the genvbounds are not sorted, i.e. if root node processing has not been finished, yet, we just propagate in
    * the order in which they have been added to genvboundstore
    */
   if( !propdata->issorted )
   {
      int j;

      assert(!propdata->sort || SCIPinProbing(scip) || SCIPgetDepth(scip) == 0);

      for( j = 0; j < propdata->ngenvbounds && *result != SCIP_CUTOFF; j++ )
      {
         if( ! SCIPvarIsActive(propdata->genvboundstore[j]->var) )
         {
            /**@todo resolve multiaggregation in exitpre */
         }
         else
         {
            SCIPdebugMsg(scip, "applying genvbound with index %d (unsorted mode)\n", j);
            SCIP_CALL( applyGenVBound(scip, prop, propdata->genvboundstore[j], global, result, nchgbds) );
         }
      }

      return SCIP_OKAY;
   }

   /* otherwise, we propagate only components affected by the latest bound changes */
   startingcomponents = global ? propdata->gstartcomponents : propdata->startcomponents;
   startingindices = global ? propdata->gstartindices : propdata->startindices;
   nindices = global ? propdata->ngindices : propdata->nindices;

   for( i = 0; i < nindices && *result != SCIP_CUTOFF; i++ )
   {
      int j;

      SCIPdebugMsg(scip, "starting in component %d at index %d\n", startingcomponents[i], startingindices[i]);
      for( j = startingindices[i]; j < propdata->componentsstart[startingcomponents[i] + 1] &&
         *result != SCIP_CUTOFF; j++ ) /*lint !e679*/
      {
         assert(j < propdata->ngenvbounds);

         if( ! SCIPvarIsActive(propdata->genvboundstore[j]->var) )
         {
            /**@todo resolve multiaggregation in exitpre */
         }
         else
         {
            SCIPdebugMsg(scip, "applying genvbound with index %d, component %d\n", j, startingcomponents[i]);
            SCIP_CALL( applyGenVBound(scip, prop, propdata->genvboundstore[j], global, result, nchgbds) );
         }
      }
   }

   /* we dont want to run again caused by this starting data */
   if( !global )
   {
      SCIP_CALL( resetLocalStartingData(scip, propdata) );
   }

   return SCIP_OKAY;
}

/** initialize propagator data */
static
SCIP_RETCODE initPropdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   int nprobvars;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->eventhdlr != NULL);

   SCIPdebugMsg(scip, "init propdata\n");

   nprobvars = SCIPgetNVars(scip);

   /* init genvboundstore */
   propdata->genvboundstoresize = 2 * nprobvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->genvboundstore), propdata->genvboundstoresize) );
   BMSclearMemoryArray(propdata->genvboundstore, propdata->genvboundstoresize);
   propdata->ngenvbounds = 0;

   /* init genvboundstore hashmaps */
   SCIP_CALL( SCIPhashmapCreate(&(propdata->lbgenvbounds), SCIPblkmem(scip), nprobvars) );
   SCIP_CALL( SCIPhashmapCreate(&(propdata->ubgenvbounds), SCIPblkmem(scip), nprobvars) );

   return SCIP_OKAY;
}

/** adds a new genvbound to genvboundstore array and sets a hashmap entry */
static
SCIP_RETCODE addNewGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the genvbounds propagator */
   GENVBOUND*            genvbound           /**< genvbound to be added */
   )
{
   SCIP_HASHMAP* hashmap;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(genvbound != NULL);
   assert(getGenVBound(scip, propdata, genvbound->var, genvbound->boundtype) == NULL);

   hashmap = genvbound->boundtype == SCIP_BOUNDTYPE_LOWER ? propdata->lbgenvbounds : propdata->ubgenvbounds;

   /* e.g., during presolving after a restart, new variables might have been created; in this case, we need to extend
    * the genvboundstore; the new size may even exceed 2*SCIPgetNVars() if we have genvbounds with nonactive left-hand
    * side variables
    */
   assert(propdata->ngenvbounds <= propdata->genvboundstoresize);
   if( propdata->ngenvbounds == propdata->genvboundstoresize )
   {
      int oldsize = propdata->genvboundstoresize;
      propdata->genvboundstoresize = 2*propdata->genvboundstoresize + 1;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(propdata->genvboundstore), oldsize, propdata->genvboundstoresize) );
   }

   /* new index is propdata->ngenvbounds */
   SCIP_CALL( SCIPhashmapInsert(hashmap, genvbound->var, genvbound) );
   propdata->genvboundstore[propdata->ngenvbounds] = genvbound;
   genvbound->index = propdata->ngenvbounds;
   propdata->ngenvbounds++;

   assert(propdata->ngenvbounds <= propdata->genvboundstoresize);

   return SCIP_OKAY;
}

/** runs propagation routine */
static
SCIP_RETCODE execGenVBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the genvbounds propagator */
   SCIP_RESULT*          result,             /**< result pointer */
   SCIP_Bool             local,              /**< should local propagation be applied? */
   int*                  nchgbds             /**< counter to increase by the number of changed bounds */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->prop != NULL);
   assert(result != NULL);

   /* we only sort after the root node is finished; this avoids having to sort again after adding more genvbounds; if
    * the genvbounds are not sorted, we will simply propagate all of them in the order given
    */
   if( propdata->sort && !SCIPinProbing(scip) && SCIPgetDepth(scip) > 0 )
   {
      if( !propdata->issorted )
      {
         *result = SCIP_DIDNOTFIND;

         SCIPdebugMsg(scip, "genvbounds are not sorted\n");

         /* drop and free old events */
         SCIP_CALL( dropAndFreeEvents(scip, propdata) );

         /* free old starting data */
         SCIP_CALL( freeStartingData(scip, propdata) );

         /* free sorted components data */
         SCIP_CALL( freeComponentsData(scip, propdata) );

         /* sort genvbounds */
         SCIP_CALL( sortGenVBounds(scip, propdata) );

         /* create starting data */
         SCIP_CALL( createStartingData(scip, propdata) );

         /* fill global starting data */
         SCIP_CALL( fillGlobalStartingData(scip, propdata) );
      }

      /* set up new events to catch (if not done so far) */
      if( propdata->lbevents == NULL )
      {
         SCIP_CALL( setUpEvents(scip, propdata) );
      }
   }

   /* apply global propagation if primal bound has improved */
   if( propdata->global && SCIPgetDepth(scip) > 0 && SCIPisFeasLT(scip, SCIPgetCutoffbound(scip), propdata->lastcutoff) )
   {
      if( propdata->ngindices > 0 )
      {
         SCIP_CALL( applyGenVBounds(scip, propdata->prop, TRUE, result, nchgbds) );
         assert(*result != SCIP_DIDNOTRUN);
      }

      /* within the tree, the objective function should not change anymore, hence the cutoff bound should be a stable
       * point of reference
       */
      propdata->lastcutoff = SCIPgetCutoffbound(scip);
   }

   /* apply local propagation if allowed */
   if( local && *result != SCIP_CUTOFF )
   {
      /* check if local propagation in root node is allowed */
      if( SCIPgetDepth(scip) > 0 || propdata->propinrootnode )
      {
         /* if genvbounds are already sorted, check if bound change events were caught; otherwise apply all genvbounds */
         if( !propdata->issorted || ( SCIPgetCurrentNode(scip) == propdata->lastnodecaught && propdata->nindices > 0 ) )
         {
            SCIP_CALL( applyGenVBounds(scip, propdata->prop, FALSE, result, nchgbds) );
            assert(*result != SCIP_DIDNOTRUN);
         }
      }
   }

   return SCIP_OKAY;
}

/* adds all genvbounds in the genvboundstore as constraints to the problem; afterwards clears the genvboundstore */
static
SCIP_RETCODE createConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the genvbounds propagator */
   )
{
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->propasconss);

   /* ensure that the cutoffboundvar is available */
   if( propdata->cutoffboundvar == NULL )
   {
      SCIP_Real ub;
      char name[16];

      /* set the upper bound to the best primal value in the original problem */
      ub = getCutoffboundGenVBound(scip);

      SCIPdebugMsg(scip, "initialize cutoffboundvar with UB = %e\n", ub);

      (void) SCIPsnprintf(name, 16, "cutoffboundvar");
      SCIP_CALL( SCIPcreateVarBasic(scip, &propdata->cutoffboundvar, name, -SCIPinfinity(scip), ub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, propdata->cutoffboundvar) );

      SCIP_CALL( SCIPaddVar(scip, propdata->cutoffboundvar) );

      /* lock the variable because it should not be subject to dual presolving reductions; because we create the
       * linear constraints as non-check constraints, the cutoffboundvar will not be locked by the linear constraint
       * handler
       */
      SCIP_CALL( SCIPaddVarLocks(scip, propdata->cutoffboundvar, 1, 1) );
   }

   assert(propdata->cutoffboundvar != NULL);

   /* now iterate over all genvbounds in the store and construct a linear constraint for each of them */
   for( i = 0; i < propdata->ngenvbounds; ++i )
   {
      GENVBOUND* genvbound;
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      char name[SCIP_MAXSTRLEN];
      int nvars;
      int j;

      genvbound = propdata->genvboundstore[i];
      assert(genvbound != NULL);

      nvars = genvbound->ncoefs + 2;
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );

      SCIPdebugMsgPrint(scip, "add cons: ");

      /* copy the coefs/vars array */
      for( j = 0; j < genvbound->ncoefs; j++ )
      {
         vars[j] = genvbound->vars[j];
         vals[j] = genvbound->coefs[j];
         SCIPdebugMsgPrint(scip, "%e%s + ", vals[j], SCIPvarGetName(vars[j]));
      }

      /* add the variable and the coefficient of the genvbound */
      vars[genvbound->ncoefs] = genvbound->var;
      vals[genvbound->ncoefs] = (genvbound->boundtype == SCIP_BOUNDTYPE_LOWER) ? -1.0 : 1.0;

      SCIPdebugMsgPrint(scip, "%e%s + ", vals[genvbound->ncoefs], SCIPvarGetName(vars[genvbound->ncoefs]));

      /* add cutoffcoef * cutoffboundvar */
      vars[genvbound->ncoefs + 1] = propdata->cutoffboundvar;
      vals[genvbound->ncoefs + 1] = genvbound->cutoffcoef;

      SCIPdebugMsgPrint(scip, "%e%s <= %e\n", vals[genvbound->ncoefs + 1], SCIPvarGetName(vars[genvbound->ncoefs + 1]), -1.0*genvbound->constant);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "genvbound_cons%d", genvbound->index);

      /* create linear constraint with only propagate flag as TRUE */
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nvars, vars, vals, -SCIPinfinity(scip), -genvbound->constant,
            FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* free memory */
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &vals);
   }

   /* now delete all genvbounds in the genvboundstore */
   if( propdata->ngenvbounds > 0 )
   {
      assert(propdata->genvboundstore != NULL);

      for( i = propdata->ngenvbounds - 1; i >= 0; i-- )
      {
         SCIP_CALL( freeGenVBound(scip, propdata->genvboundstore[i]) );
      }

      /* free genvboundstore hashmaps */
      SCIPhashmapFree(&(propdata->lbgenvbounds));
      SCIPhashmapFree(&(propdata->ubgenvbounds));

      /* drop and free all events */
      SCIP_CALL( dropAndFreeEvents(scip, propdata) );

      /* free componentsstart array */
      SCIP_CALL( freeComponentsData(scip, propdata) );

      /* free starting indices data */
      SCIP_CALL( freeStartingData(scip, propdata) );

      SCIPfreeBlockMemoryArray(scip, &(propdata->genvboundstore), propdata->genvboundstoresize);
      propdata->genvboundstore = NULL;
      propdata->genvboundstoresize = 0;
      propdata->ngenvbounds = 0;
   }

   return SCIP_OKAY;
}



/*
 * Public methods
 */

/** adds a generalized variable bound to the genvbounds propagator; if there is already a genvbound for the bound
 *  "boundtype" of variable "var", it will be replaced
 */
SCIP_RETCODE SCIPgenVBoundAdd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            genvboundprop,      /**< genvbound propagator */
   SCIP_VAR**            vars,               /**< array of RHSs variables */
   SCIP_VAR*             var,                /**< LHSs variable */
   SCIP_Real*            coefs,              /**< array of coefficients for the RHSs variables */
   int                   ncoefs,             /**< size of coefs array */
   SCIP_Real             coefcutoffbound,    /**< nonpositive value of the cutoff bounds multiplier */
   SCIP_Real             constant,           /**< constant term */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound provided by the genvbound */
   )
{
   /**@todo in debug mode: check if genvbound is nontrivial */

   SCIP_PROPDATA* propdata;
   GENVBOUND* genvbound;
   SCIP_Bool newgenvbound;
   int i;

   assert(scip != NULL);
   assert(genvboundprop != NULL);
   assert(strcmp(SCIPpropGetName(genvboundprop), PROP_NAME) == 0);
   assert(vars != NULL);
   assert(var != NULL);
   assert(coefs != NULL);
   assert(ncoefs >= 0);
   assert(coefcutoffbound <= 0.0);
   assert(!SCIPisInfinity(scip, -constant));

   if( ncoefs < 0 || coefcutoffbound > 0.0 || SCIPisInfinity(scip, -constant) )
   {
      SCIPerrorMessage("cannot create generalized variable bound from invalid data\n");
      return SCIP_INVALIDDATA;
   }

   propdata = SCIPpropGetData(genvboundprop);
   assert(propdata != NULL);

   /* initialize propdata if not done yet */
   if( propdata->genvboundstore == NULL )
   {
      SCIP_CALL( initPropdata(scip, propdata) );
   }

   genvbound = getGenVBound(scip, propdata, var, boundtype);
   newgenvbound = (genvbound == NULL);

   /* release previous variables */
   if( !newgenvbound )
   {
      for( i = 0; i < genvbound->ncoefs; ++i )
      {
         assert(genvbound->vars[i] != NULL);
         SCIP_CALL( SCIPreleaseVar(scip, &(genvbound->vars[i])) );
      }
   }

   /* check if there already is a genvbound corresponding to this bound, freeing its data and overwriting it */
   if( !newgenvbound && genvbound->ncoefs < ncoefs )
   {
      /* do not realloc since we do not want to keep and possibly copy the old entries */
      SCIPfreeBlockMemoryArray(scip, &(genvbound->coefs), genvbound->coefssize);
      SCIPfreeBlockMemoryArray(scip, &(genvbound->vars), genvbound->coefssize);

      /* allocate and copy arrays in genvbound */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(genvbound->coefs), coefs, ncoefs) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(genvbound->vars), vars, ncoefs) );
      genvbound->coefssize = ncoefs;
   }
   else if( !newgenvbound && genvbound->ncoefs == ncoefs )
   {
      /* just update entries */
      for( i = 0; i < ncoefs; i++ )
      {
         genvbound->coefs[i] = coefs[i];
         genvbound->vars[i] = vars[i];
      }
   }
   else if( !newgenvbound && genvbound->ncoefs > ncoefs )
   {
      /* reallocate memory for arrays in genvbound to free unused memory */
      if( genvbound->coefssize < ncoefs )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(genvbound->coefs), genvbound->coefssize, ncoefs) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(genvbound->vars), genvbound->coefssize, ncoefs) );
         genvbound->coefssize = ncoefs;
      }

      /* update entries */
      for( i = 0; i < ncoefs; i++ )
      {
         genvbound->coefs[i] = coefs[i];
         genvbound->vars[i] = vars[i];
      }
   }
   else if( newgenvbound )
   {
      /* allocate memory for genvbound data */
      SCIP_CALL( SCIPallocBlockMemory(scip, &genvbound) );

      /* allocate and copy arrays in genvbound */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(genvbound->coefs), coefs, ncoefs) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(genvbound->vars), vars, ncoefs) );
      genvbound->coefssize = ncoefs;
   }

   /* set up data for genvbound */
   genvbound->boundtype = boundtype;
   genvbound->var = var;
   genvbound->ncoefs = ncoefs;
   genvbound->constant = constant;

   /* capture variables */
   for( i = 0; i < genvbound->ncoefs; ++i )
   {
      assert(genvbound->vars[i] != NULL);
      SCIP_CALL( SCIPcaptureVar(scip, genvbound->vars[i]) );
   }
   if( newgenvbound )
   {
      assert(genvbound->var != NULL);
      SCIP_CALL( SCIPcaptureVar(scip, genvbound->var) );
   }

   /* the cutoff bound is valid w.r.t. the current objective function in the transformed problem; during presolving,
    * however, the objective function can change (e.g., when a variable is fixed, its contribution in the objective
    * is subtracted from the cutoff bound and added to the objective offset); we solve this by transforming the
    * contribution of the cutoff bound in the generalized variable bound to the original problem as follows:
    *
    *    +/- var >= ... + z * SCIPgetCutoffbound() + constant
    *
    * becomes
    *
    *    +/- var >= ... + (z / SCIPgetTransObjscale()) * origcutoffbound + (constant - z * SCIPgetTransObjoffset())
    *
    * with SCIPgetCutoffbound() = origcutoffbound / SCIPgetTransObjscale() - SCIPgetTransObjoffset(); in the
    * propagation later, we will use (SCIPgetCutoffbound() + SCIPgetTransObjoffset()) * SCIPgetTransObjscale(), see
    * function getCutoffboundGenVBound()
    */
   if( SCIPisNegative(scip, coefcutoffbound) )
   {
      assert(SCIPisPositive(scip, SCIPgetTransObjscale(scip)));
      genvbound->cutoffcoef = coefcutoffbound / SCIPgetTransObjscale(scip);
      genvbound->constant -= (coefcutoffbound * SCIPgetTransObjoffset(scip));
   }
   else
      genvbound->cutoffcoef = 0.0;

   /* if genvbound is not overwritten, create a new entry in genvboundstore */
   if( newgenvbound )
   {
      SCIP_CALL( addNewGenVBound(scip, propdata, genvbound) );
   }

   /* mark genvbounds array to be resorted */
   propdata->issorted = FALSE;

   /* debug message */
   SCIPdebugMsg(scip, "added genvbound ");
   SCIPdebug( printGenVBound(scip, genvbound) );
#ifdef WITH_DEBUG_SOLUTION
   SCIP_CALL( checkDebugSolutionGenVBound(scip, genvbound) );
#endif

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

/** copy method for propagator plugins (called when SCIP copies plugins)
 *
 *  @note The UG framework assumes that all default plug-ins of SCIP implement a copy callback.
 */
static
SCIP_DECL_PROPCOPY(propCopyGenvbounds)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludePropGenvbounds(scip) );

   return SCIP_OKAY;
}

/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitGenvbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->genvboundstore = NULL;
   propdata->genvboundstoresize = 0;
   propdata->lbevents = NULL;
   propdata->ubevents = NULL;
   propdata->lbgenvbounds = NULL;
   propdata->ubgenvbounds = NULL;
   propdata->lbeventsmap = NULL;
   propdata->ubeventsmap = NULL;
   propdata->startmap = NULL;
   propdata->componentsstart = NULL;
   propdata->startindices = NULL;
   propdata->startcomponents = NULL;
   propdata->gstartindices = NULL;
   propdata->gstartcomponents = NULL;
   propdata->lastcutoff = SCIPinfinity(scip);
   propdata->lastnodecaught = NULL;
   propdata->cutoffboundvar = NULL;
   propdata->ngenvbounds = -1;
   propdata->ncomponents = -1;
   propdata->nindices = -1;
   propdata->ngindices = -1;
   propdata->nlbevents = -1;
   propdata->nubevents = -1;
   propdata->issorted = FALSE;

   propdata->prop = prop;

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolGenvbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   *result = SCIP_DIDNOTRUN;

   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPdebugMsg(scip, "proppresol in problem <%s>\n", SCIPgetProbName(scip));

   /* do not run if no genvbounds were added yet */
   if( propdata->ngenvbounds < 1 )
   {
      SCIPdebugMsg(scip, "no bounds were added yet\n");
      return SCIP_OKAY;
   }

   /* propagate */
   SCIP_CALL( execGenVBounds(scip, propdata, result, TRUE, nchgbds) );

   return SCIP_OKAY;
}


/** presolving initialization method of propagator (called when presolving is about to begin) */
static
SCIP_DECL_PROPINITPRE(propInitpreGenvbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* lock the variable because it should not be deleted after a restart */
   if( propdata->cutoffboundvar != NULL )
   {
      SCIPdebugMsg(scip, "propinitpre in problem <%s>: locking cutoffboundvar (current downlocks=%d, uplocks=%d)\n",
         SCIPgetProbName(scip), SCIPvarGetNLocksDown(propdata->cutoffboundvar),
         SCIPvarGetNLocksUp(propdata->cutoffboundvar));

      SCIP_CALL( SCIPaddVarLocks(scip, propdata->cutoffboundvar, 1, 1) );
   }

   return SCIP_OKAY;
}


/** presolving deinitialization method of propagator (called after presolving has been finished) */
static
SCIP_DECL_PROPEXITPRE(propExitpreGenvbounds)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_PROPDATA* propdata;
   int i;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   SCIPdebugMsg(scip, "propexitpre in problem <%s>: removing fixed, aggregated, negated, and multi-aggregated variables from right-hand side\n",
      SCIPgetProbName(scip));

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* there should be no events on the right-hand side variables */
   assert(propdata->lbevents == NULL);
   assert(propdata->ubevents == NULL);

   /* allocate memory to store new variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, SCIPgetNTotalVars(scip)) );

   for( i = 0; i < propdata->ngenvbounds; )
   {
      GENVBOUND* genvbound;
      int requiredsize;
      int nvars;
      int j;

      genvbound = propdata->genvboundstore[i];
      assert(genvbound != NULL);

      /* store variables of the genvbound to release them properly */
      assert(genvbound->ncoefs <= SCIPgetNTotalVars(scip));
      BMScopyMemoryArray(vars, genvbound->vars, genvbound->ncoefs);
      nvars = genvbound->ncoefs;

      /* replace non-active by active variables and update constant; note that this may result in coefficients where
       * SCIPisZero() is true; this should not create any problems
       */
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, genvbound->vars, genvbound->coefs, &genvbound->ncoefs, genvbound->ncoefs, &genvbound->constant, &requiredsize, TRUE) );

      /* if space was not enough we need to resize the buffers */
      if( requiredsize > genvbound->ncoefs )
      {
         /* reallocate memory for arrays in genvbound to free unused memory */
         if( genvbound->coefssize < requiredsize )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(genvbound->coefs), genvbound->coefssize, requiredsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(genvbound->vars), genvbound->coefssize, requiredsize) );
            genvbound->coefssize = requiredsize;
         }

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, genvbound->vars, genvbound->coefs, &genvbound->ncoefs, requiredsize, &genvbound->constant, &requiredsize, TRUE) );
         assert(requiredsize <= genvbound->ncoefs);
      }

      /* capture new and release old variables */
      for( j = 0; j < genvbound->ncoefs; ++j )
      {
         assert(genvbound->vars[j] != NULL);
         SCIP_CALL( SCIPcaptureVar(scip, genvbound->vars[j]) );
      }
      for( j = 0; j < nvars; ++j )
      {
         assert(vars[j] != NULL);
         SCIP_CALL( SCIPreleaseVar(scip, &vars[j]) );
      }

      /* if the resulting genvbound is trivial, remove it */
      /* we remove all genvbounds with an aggregated or multi-aggregated genvbound->var; tightening aggregated variables
       * might lead to some asserts in tree.c if the active variable has been already tightened (see !398);
       *
       * @todo replace aggregated variable by their active part
       */
      if( (genvbound->ncoefs == 0 && SCIPisZero(scip, genvbound->cutoffcoef))
         || SCIPvarGetStatus(genvbound->var) == SCIP_VARSTATUS_MULTAGGR
         || SCIPvarGetStatus(genvbound->var) == SCIP_VARSTATUS_AGGREGATED )
      {
         SCIP_HASHMAP* hashmap;

         hashmap = genvbound->boundtype == SCIP_BOUNDTYPE_LOWER ? propdata->lbgenvbounds : propdata->ubgenvbounds;

         /* remove genvbound from hashmap */
         assert(SCIPhashmapExists(hashmap, genvbound->var));
         SCIP_CALL( SCIPhashmapRemove(hashmap, genvbound->var) );

         /* free genvbound and fill gap */
         SCIP_CALL( freeGenVBound(scip, propdata->genvboundstore[i]) );
         --(propdata->ngenvbounds);

         /* move the last genvbound to the i-th position */
         if( i < propdata->ngenvbounds )
         {
            propdata->genvboundstore[i] = propdata->genvboundstore[propdata->ngenvbounds];
            propdata->genvboundstore[i]->index = i;

            /* mark genvbounds array to be resorted */
            propdata->issorted = FALSE;
         }
      }
      else
         ++i;
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitGenvbounds)
{
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* free remaining genvbounds */
   SCIP_CALL( freeGenVBounds(scip, propdata) );

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecGenvbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   *result = SCIP_DIDNOTRUN;

   /* do not run if propagation w.r.t. current objective is not allowed */
   if( !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* update upper bound of the cutoffboundvar */
   if( propdata->cutoffboundvar != NULL )
   {
      SCIP_Real newub;
      SCIP_Real oldub;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      assert(propdata->propasconss);

      /* compute the primal bound in the original problem */
      newub = getCutoffboundGenVBound(scip);
      oldub = SCIPvarGetUbLocal(propdata->cutoffboundvar);

      if( SCIPisInfinity(scip, newub) == FALSE && SCIPisFeasLT(scip, newub, oldub) )
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, propdata->cutoffboundvar, newub, FALSE, &infeasible, &tightened) );

         if( tightened )
         {
            SCIPdebugMsg(scip, "tightened UB of cutoffboundvar to %e (old: %e, infeas: %u, tightened: %u)\n",
               newub, oldub, infeasible, tightened);
         }

         assert(infeasible == FALSE);
      }
   }

   SCIPdebugMsg(scip, "propexec in problem <%s> at depth %d%s\n", SCIPgetProbName(scip), SCIPgetDepth(scip),
      SCIPinProbing(scip) ? " in probing" : "");

   /* do not run if no genvbounds were added yet */
   if( propdata->ngenvbounds < 1 )
   {
      /**@todo is it really no performance issue to be called each time when there are no genvbounds, e.g., for MIPs? */
      SCIPdebugMsg(scip, "no bounds were added yet\n");
      return SCIP_OKAY;
   }

   /* add the genvbounds in the genvboundstore as constraints to the problem; afterwards clear the genvboundstore */
   if( propdata->propasconss )
   {
      SCIP_CALL( createConstraints(scip, propdata) );
      return SCIP_OKAY;
   }

   /* propagate locally and globally */
   SCIP_CALL( execGenVBounds(scip, propdata, result, !SCIPinProbing(scip), NULL) );

   /* when called in presolving stage the result is set to SCIP_SUCCESS instead of SCIP_REDUCEDDOM, this is corrected
    * here
    */
   if( *result == SCIP_SUCCESS )
      *result = SCIP_REDUCEDDOM;

   SCIPdebugMsg(scip, "end of exec\n");

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropGenvbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   GENVBOUND* genvbound;
   SCIP_Real boundval;
   SCIP_Bool success;

   SCIPdebugMsg(scip, "explain %s bound change of variable <%s>\n",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(infervar));

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->genvboundstore != NULL);

   /* as inferinfo we passed the index of the genvbound that was used for propagation; the genvbound might have been
    * replaced, but also the new genvbound at this position has the same variable on the left-hand side
    */
   assert(inferinfo >= 0);
   assert(inferinfo < propdata->ngenvbounds);

   *result = SCIP_DIDNOTFIND;

   /* check also in optimized mode that inferinfo is correct */
   if( inferinfo >= propdata->ngenvbounds)
   {
      SCIPerrorMessage("generalized variable bounds propagator received inferinfo out of range; propagation not resolved, safe to continue\n");
      return SCIP_OKAY;
   }

   /* get genvbound responsible for the bound change */
   genvbound = propdata->genvboundstore[inferinfo];
   assert(genvbound != NULL);
   assert(genvbound->var == infervar);

   /* check also in optimized mode that inferinfo is correct */
   if( genvbound->var != infervar )
   {
      SCIPerrorMessage("generalized variable bounds propagator received incorrect inferinfo; propagation not resolved, but it's safe to continue\n");
      return SCIP_OKAY;
   }

   /* get value of bound change on left-hand side */
   boundval = genvbound->boundtype == SCIP_BOUNDTYPE_LOWER
      ? SCIPgetVarLbAtIndex(scip, genvbound->var, bdchgidx, TRUE)
      : -SCIPgetVarUbAtIndex(scip, genvbound->var, bdchgidx, TRUE);

   /* if left-hand side variable is integral, it suffices to explain a bound change greater than boundval - 1 */
   if( SCIPvarIsIntegral(genvbound->var) )
   {
      SCIP_Real roundedboundval;

      assert(SCIPisIntegral(scip, boundval));

      roundedboundval = SCIPfeasCeil(scip, boundval - 1.0) + 2 * SCIPfeastol(scip);
      boundval = MIN(boundval, roundedboundval);
   }

   /* resolve propagation */
   SCIP_CALL( resolveGenVBoundPropagation(scip, genvbound, bdchgidx, &boundval, &success) );

   if( success )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolGenvbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   SCIPdebugMsg(scip, "propexitsol in problem <%s>\n", SCIPgetProbName(scip));

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* free all genvbounds if we are not in a restart */
   if( !SCIPisInRestart(scip) )
   {
      SCIP_CALL( freeGenVBounds(scip, propdata) );
   }

   /* drop and free all events */
   SCIP_CALL( dropAndFreeEvents(scip, propdata) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeGenvbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeBlockMemory(scip, &propdata);

   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecGenvbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int i;

   assert(scip != NULL);
   assert(eventdata != NULL);

   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_LBTIGHTENED || SCIPeventGetType(event) ==
      SCIP_EVENTTYPE_UBTIGHTENED);

   assert(eventdata->startcomponents != NULL);
   assert(eventdata->startindices != NULL);
   assert(eventdata->nstarts > 0);
   assert(eventdata->prop != NULL);

   propdata = SCIPpropGetData(eventdata->prop);
   assert(propdata != NULL);

   assert(propdata->startcomponents != NULL);
   assert(propdata->startmap != NULL);
   assert(propdata->startindices != NULL);

   SCIPdebugMsg(scip, "catching eventdata:\n");
   SCIPdebug( printEventData(eventdata, SCIPeventGetType(event) == SCIP_EVENTTYPE_LBTIGHTENED ?
         SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER) );

   /* check if we need to reset old local starting indices data */
   if( SCIPgetCurrentNode(scip) != propdata->lastnodecaught )
   {
      SCIP_CALL( resetLocalStartingData(scip, propdata) );
      propdata->lastnodecaught = SCIPgetCurrentNode(scip);
   }

   for( i = 0; i < eventdata->nstarts; i++ )
   {
      int component;
      int startidx;

      component = eventdata->startcomponents[i];
      assert(component >= 0);
      startidx = eventdata->startindices[i];

      /* there is already an entry for this component */
      if( SCIPhashmapExists(propdata->startmap, (void*)(size_t) (component + 1)) )
      {
         int componentidx;

         /* get its index */
         componentidx = ((int)(size_t) SCIPhashmapGetImage(propdata->startmap, (void*)(size_t) (component + 1))) - 1; /*lint !e776*/
         assert(componentidx >= 0);
         assert(propdata->startcomponents[componentidx] == component);

         if( propdata->startindices[componentidx] > startidx )
            propdata->startindices[componentidx] = startidx;
      }
      else
      {
         /* get a new entry */
         int componentidx;
         componentidx = propdata->nindices;

         /* store index */
         propdata->startcomponents[componentidx] = component;
         propdata->startindices[componentidx] = startidx;

         /* store component in hashmap */
         SCIP_CALL( SCIPhashmapInsert(propdata->startmap, (void*)(size_t) (component + 1),
               (void*)(size_t) (componentidx + 1)) );

         /* increase number of starting indices */
         propdata->nindices++;
      }
   }

   return SCIP_OKAY;
}

/*
 * propagator specific interface methods
 */

/** creates the genvbounds propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropGenvbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create genvbounds propagator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecGenvbounds, propdata) );

   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyGenvbounds) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeGenvbounds) );
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitGenvbounds) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreGenvbounds) );
   SCIP_CALL( SCIPsetPropExitpre(scip, prop, propExitpreGenvbounds) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitGenvbounds) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolGenvbounds) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolGenvbounds, PROP_PRESOL_PRIORITY,
         PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropGenvbounds) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/global",
         "apply global propagation?",
         &propdata->global, TRUE, DEFAULT_GLOBAL_PROPAGATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/propinrootnode",
         "apply genvbounds in root node if no new incumbent was found?",
         &propdata->propinrootnode, TRUE, DEFAULT_PROPAGATE_IN_ROOT_NODE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/sort",
         "sort genvbounds and wait for bound change events?",
         &propdata->sort, TRUE, DEFAULT_SORT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/propasconss",
         "should genvbounds be transformed to (linear) constraints?",
         &propdata->propasconss, TRUE, DEFAULT_PROPASCONSS, NULL, NULL) );

   /* include event handler */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &propdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecGenvbounds, NULL) );

   return SCIP_OKAY;
}
