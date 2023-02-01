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

/**@file   cons_pseudoboolean.c
 * @brief  constraint handler for pseudo Boolean constraints
 * @author Gerald Gamrath
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 *
 * The constraint handler deals with pseudo Boolean constraints. These are constraints of the form
 * \f[
 * \mbox{lhs} \leq \sum_{k=0}^m c_k \cdot x_k  +  \sum_{i=0}^n c_i \cdot \prod_{j \in I_i} x_j \leq \mbox{rhs}
 * \f]
 * where all x are binary and all c are integer
 *
 * @todo Add eventhandling.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_pseudoboolean.h"
#include "scip/cons_and.h"
#include "scip/cons_indicator.h"
#ifdef WITHEQKNAPSACK
#include "scip/cons_eqknapsack.h"
#endif
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_xor.h"
#include "scip/pub_var.h"
#include "scip/debug.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "pseudoboolean"
#define CONSHDLR_DESC          "constraint handler dealing with pseudo Boolean constraints"
#define CONSHDLR_ENFOPRIORITY  -1000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -5000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */

#define DEFAULT_DECOMPOSENORMALPBCONS FALSE /**< decompose all normal pseudo boolean constraint into a "linear" constraint and "and" constraints */
#define DEFAULT_DECOMPOSEINDICATORPBCONS TRUE /**< decompose all indicator pseudo boolean constraint into a "linear" constraint and "and" constraints */

#define DEFAULT_SEPARATENONLINEAR  TRUE /**< if decomposed, should the nonlinear constraints be separated during LP processing */
#define DEFAULT_PROPAGATENONLINEAR TRUE /**< if decomposed, should the nonlinear constraints be propagated during node processing */
#define DEFAULT_REMOVABLENONLINEAR TRUE /**< if decomposed, should the nonlinear constraints be removable */
#define USEINDICATOR               TRUE

/*
 * Data structures
 */
#define HASHSIZE_PSEUDOBOOLEANNONLINEARTERMS 500 /**< minimal size of hash table in and constraint tables */


/* - create special linear(knapsack, setppc, logicor, (eqknapsack)) and and-constraints with check flags FALSE, to
 *   get smaller amount of locks on the term variables, do all presolving ...?! in these constraint handlers
 *
 * - do the checking here, lock and-resultants in both directions and all and-variables according to their
 *   coefficients and sides of the constraint,
 *   @note this only works if the and-resultant has no objective cofficient, otherwise we need to lock variables also in both directions
 *
 * - need to keep and constraint pointer for special propagations like if two ands are due to their variables in
 *   one clique, add this cliques of and-resultants
 *
 * - do special presolving like on instance :
 * check/IP/PseudoBoolean/normalized-PB07/OPT-SMALLINT-NLC/submittedPB07/manquinho/bsg/normalized-bsg_1000_25_1.opb.gz
 *
 *  there exist constraint like:        1 x1 x2 + 1 x1 x3 + 1 x1 x4 + 1 x1 x5 <= 1 ;
 *  which "equals" a linear constraint: 3 x1 + x2 + x3 + x4 + x5 <= 4 ;
 *
 *  in more general terms:                     1 x1 x2 x3 x4 + 1 x1 x2 x5 x6 x7 + 1 x1 x2 x8 x9 <= 1 ;
 *  which "equals" a pseudoboolean constraint: 2 x1 + 2 x2 + 1 x3 x4 + 1 x5 x6 x7 + 1 x8 x9 <= 5 ;
 *
 *  in an even more general terms:             5 x1 x2 x3 x4 + 1 x1 x2 x5 x6 x7 + 1 x1 x2 x8 x9 <= 6 ;
 *            equals(should the knapsack do)   1 x1 x2 x3 x4 + 1 x1 x2 x5 x6 x7 + 1 x1 x2 x8 x9 <= 2 ;
 *  which "equals" a pseudoboolean constraint: 2 x1 + 2 x2 + 1 x3 x4 + 1 x5 x6 x7 + 1 x8 x9 <= 6 ;
 *  (         without knapsack                 7 x1 + 7 x2 + 5 x3 x4 + 1 x5 x6 x7 + 1 x8 x9 <= 20 ; )
 *
 *  another special case :                     1 x1 x2 x3 + 1 x1 x2 x4 + 1 x5 x6 <= 1 ;
 *  which "equals" a pseudoboolean constraint: 2 x1 + 2 x2 + 1 x3 + 1 x4 + 1 x5 x6 <= 5 ;
 *  which "equals" a pseudoboolean constraint: 4 x1 + 4 x2 + 2 x3 + 2 x4 + 1 x5 + 1 x6 <= 10 ;
 *
 *  another special case :                     1 x1 x2 + 1 x1 x3 + 2 x4 x5 <= 3 ;
 *  which "equals" a pseudoboolean constraint: 2 x1 + 1 x2 + 1 x3 + 2 x4 x5 <= 5 ;
 *  which "equals" a pseudoboolean constraint: 2 x1 + 1 x2 + 1 x3 + 1 x4 + 1 x5 <= 5 ;
 */
/* @todo - in and-constraint better count nfixed zeros in both directions and maybe nfixedones for better propagation
 *
 *       - do better conflict analysis by choosing the earliest fixed variable which led to a conflict instead of maybe
 *         best coefficient or create more conflicts by using all to zero fixed variables one by one
 *
 *       - how to make sure that we aggregate in a right way, when aggregating a resultant and a "normal" variable,
 *         maybe add in SCIPaggregateVars a check for original variables, to prefer them if the variable type is the
 *         same; probably it would be better too if we would aggregate two resultants that the one with less variables
 *         inside the and-constraint will stay active
 *
 * @note since product resultants are artificial, we do not care for their solution value, but this can lead to fixation
 *       of the resultant not representing the product, in 'optimization mode' we do not care, but this might make
 *       solution debugging complicated
 */

/** and-constraint data object */
struct ConsAndData
{
   SCIP_CONS*            cons;                /**< pointer to the and-constraint of this 'term' of variables */
   SCIP_CONS*            origcons;            /**< pointer to the original and-constraint of this 'term' of variables
                                               *   only after problem was transformed, NULL otherwise */
   SCIP_VAR**            vars;                /**< all and-constraint variables */
   int                   nvars;               /**< number of all and-constraint variables */
   int                   svars;               /**< size for all and-constraint variables */
   SCIP_VAR**            newvars;             /**< new variables in this presolving round */
   int                   nnewvars;            /**< number of new variables in this presolving round */
   int                   snewvars;            /**< size of new variables in this presolving round */
   int                   noriguses;           /**< how often is this data in used by original constraints */
   int                   nuses;               /**< how often is this data in used by transformed constraints */
   unsigned int          istransformed:1;     /**< is transformed data active */
   unsigned int          isoriginal:1;        /**< is original data active */
};
typedef struct ConsAndData CONSANDDATA;

/** constraint data for pseudoboolean constraints */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of constraint */
   SCIP_Real             rhs;                /**< right hand side of constraint */

   SCIP_CONS*            lincons;            /**< linear constraint which represents this pseudoboolean constraint */
   SCIP_LINEARCONSTYPE   linconstype;        /**< type of linear constraint which represents this pseudoboolean constraint */
   int                   nlinvars;           /**< number of linear variables (without and-resultants) */

   CONSANDDATA**         consanddatas;       /**< array of and-constraints-data-objects sorted after index of
                                              *   and-resultant of corresponding and-constraint */
   SCIP_Real*            andcoefs;           /**< array of coefficients for and-constraints of
                                              *   and-constraints-data-objects
                                              *   (changes during presolving, needs to be updated in every presolving
                                              *   round) */
   SCIP_Bool*            andnegs;            /**< array of negation status for and-constraints of
                                              *   and-constraints-data-objects
                                              *   (changes during presolving, needs to be updated in every presolving
                                              *   round) */
   int                   nconsanddatas;      /**< number of and-constraints-data-objects */
   int                   sconsanddatas;      /**< size of and-constraints-data-objects array */

   SCIP_VAR*             intvar;             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */

   SCIP_VAR*             indvar;             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight;             /**< weight of the soft constraint, if it is one */

   unsigned int          issoftcons:1;       /**< is this a soft constraint */
   unsigned int          changed:1;          /**< was constraint changed? */
   unsigned int          propagated:1;       /**< is constraint already propagated? */
   unsigned int          presolved:1;        /**< is constraint already presolved? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the constraint already extracted? */
   unsigned int          upgradetried:1;     /**< was constraint upgrading already tried */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   CONSANDDATA**         allconsanddatas;    /**< array of all and-constraint data objects inside the whole problem,
                                              *   created via this constraint handler */
   int                   nallconsanddatas;   /**< number of all and-constraint data objects inside the whole problem,
                                              *   created via this constraint handler */
   int                   sallconsanddatas;   /**< size of all and-constraint data objects inside the whole problem,
                                              *   created via this constraint handler */
   SCIP_HASHTABLE*       hashtable;          /**< hash table for all and-constraint data objects */
   int                   hashtablesize;      /**< size for hash table for all and-constraint data objects */

   SCIP_HASHMAP*         hashmap;            /**< hash map for mapping all resultant to and-constraint */
   int                   hashmapsize;        /**< size for hash map for mapping all resultant to and-constraint */

   SCIP_Bool             decomposenormalpbcons;/**< decompose the pseudo boolean constraint into a "linear" constraint and "and" constraints */
   SCIP_Bool             decomposeindicatorpbcons;/**< decompose the indicator pseudo boolean constraint into a "linear" constraint and "and" constraints */
   SCIP_Bool             inithashmapandtable;/**< flag to store if the hashmap and -table is initialized */
   int                   nlinconss;          /**< for counting number of created linear constraints */
   int                   noriguses;          /**< how many consanddata objects are used by original constraints */
};

/*
 * Local methods
 */


/** comparison method for sorting consanddatas according to the index of their corresponding resultant variables, if a
 *  consanddata object is delete it is handled like it has an inactive resultant, so this will be put in front while
 *  sorting
 */
static
SCIP_DECL_SORTPTRCOMP(resvarCompWithInactive)
{
   CONSANDDATA* consanddata1;
   CONSANDDATA* consanddata2;

   consanddata1 = (CONSANDDATA*)elem1;
   consanddata2 = (CONSANDDATA*)elem2;

   /* check if and constraint data object is still valid */
   if( !consanddata1->istransformed )
   {
      if( !consanddata2->istransformed )
      {
         return 0;
      }
      else
         return -1;
   }
   else if( !consanddata2->istransformed )
      return +1;

   assert(consanddata1->cons != NULL);
   assert(consanddata2->cons != NULL);

   /* check if and constraint is still active */
   if( SCIPconsIsDeleted(consanddata1->cons) )
   {
      if( SCIPconsIsDeleted(consanddata2->cons) )
      {
         return 0;
      }
      else
         return -1;
   }
   else if( SCIPconsIsDeleted(consanddata2->cons) )
      return +1;
   else
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;

      /* hack with setting the first pointer to NULL */
      var1 = SCIPgetResultantAnd(NULL, consanddata1->cons);
      var2 = SCIPgetResultantAnd(NULL, consanddata2->cons);

      assert(var1 != NULL);
      assert(var2 != NULL);

      if( SCIPvarGetIndex(var1) < SCIPvarGetIndex(var2) )
         return -1;
      else if( SCIPvarGetIndex(var1) > SCIPvarGetIndex(var2) )
         return +1;
      else
      {
         assert(var1 == var2);
         return 0;
      }
   }
}

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyAndConsDatas)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal; two non-linear terms are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqAndConsDatas)
{
#ifndef NDEBUG
   SCIP* scip;
#endif
   CONSANDDATA* cdata1;
   CONSANDDATA* cdata2;
   int v;

   cdata1 = (CONSANDDATA*)key1;
   cdata2 = (CONSANDDATA*)key2;

#ifndef NDEBUG
   scip = (SCIP*)userptr;
#endif
   assert(scip != NULL);
   assert(cdata1 != NULL);
   assert(cdata2 != NULL);
   assert(cdata1->vars != NULL);
   assert(cdata1->nvars > 1);
   assert(cdata2->vars != NULL);
   assert(cdata2->nvars > 1);

#ifndef NDEBUG
   /* check that cdata1 variables are sorted */
   for( v = cdata1->nvars - 1; v > 0; --v )
      assert(SCIPvarGetIndex(cdata1->vars[v]) >= SCIPvarGetIndex(cdata1->vars[v - 1]));
   /* check that cdata2 variables are sorted */
   for( v = cdata2->nvars - 1; v > 0; --v )
      assert(SCIPvarGetIndex(cdata2->vars[v]) >= SCIPvarGetIndex(cdata2->vars[v - 1]));
#endif

   /* checks trivial case */
   if( cdata1->nvars != cdata2->nvars )
      return FALSE;

   /* checks trivial case */
   if( cdata1->cons != NULL && cdata2->cons != NULL && cdata1->cons != cdata2->cons )
      return FALSE;

   /* check each variable in both cdatas for equality */
   for( v = cdata1->nvars - 1; v >= 0; --v )
   {
      assert(cdata1->vars[v] != NULL);
      assert(cdata2->vars[v] != NULL);

      /* tests if variables are equal */
      if( cdata1->vars[v] != cdata2->vars[v] )
      {
         assert(SCIPvarCompare(cdata1->vars[v], cdata2->vars[v]) == 1 ||
            SCIPvarCompare(cdata1->vars[v], cdata2->vars[v]) == -1);
         return FALSE;
      }
      assert(SCIPvarCompare(cdata1->vars[v], cdata2->vars[v]) == 0);
   }

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValAndConsDatas)
{  /*lint --e{715}*/
   CONSANDDATA* cdata;
   int minidx;
   int mididx;
   int maxidx;

   cdata = (CONSANDDATA*)key;

   assert(cdata != NULL);
   assert(cdata->vars != NULL);
   assert(cdata->nvars > 1);
#ifndef NDEBUG
   {
      /* check that these variables are sorted */
      int v;
      for( v = cdata->nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(cdata->vars[v]) >= SCIPvarGetIndex(cdata->vars[v - 1]));
   }
#endif

   minidx = SCIPvarGetIndex(cdata->vars[0]);
   mididx = SCIPvarGetIndex(cdata->vars[cdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(cdata->vars[cdata->nvars - 1]);
   assert(minidx >= 0 && minidx <= maxidx);

   return SCIPhashTwo(SCIPcombineTwoInt(cdata->nvars, minidx),
                      SCIPcombineTwoInt(mididx, maxidx)); /*lint !e701*/
}

/** initializes the hashmap and -table used in this constraint handler data for artificial variables and specific
 *  and-constraint data objects
 */
static
SCIP_RETCODE inithashmapandtable(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   if( ((*conshdlrdata)->inithashmapandtable) )
   {
      assert((*conshdlrdata)->hashtable != NULL);
      assert((*conshdlrdata)->hashmap != NULL);

      return SCIP_OKAY;
   }

   assert((*conshdlrdata)->hashtable == NULL);
   assert((*conshdlrdata)->hashmap == NULL);

   /* create a hash table for and-constraint data objects */
   (*conshdlrdata)->hashtablesize = HASHSIZE_PSEUDOBOOLEANNONLINEARTERMS;
   SCIP_CALL( SCIPhashtableCreate(&((*conshdlrdata)->hashtable), SCIPblkmem(scip), (*conshdlrdata)->hashtablesize,
         hashGetKeyAndConsDatas, hashKeyEqAndConsDatas, hashKeyValAndConsDatas, (void*) scip) );

   /* create a hash table for and-resultant to and-constraint data objects */
   (*conshdlrdata)->hashmapsize = HASHSIZE_PSEUDOBOOLEANNONLINEARTERMS;
   SCIP_CALL( SCIPhashmapCreate(&((*conshdlrdata)->hashmap), SCIPblkmem(scip), (*conshdlrdata)->hashmapsize) );

   (*conshdlrdata)->inithashmapandtable = TRUE;

   return SCIP_OKAY;
}

/** creates constraint handler data for pseudo boolean constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, conshdlrdata) );

   (*conshdlrdata)->allconsanddatas = NULL;
   (*conshdlrdata)->nallconsanddatas = 0;
   (*conshdlrdata)->sallconsanddatas = 10;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*conshdlrdata)->allconsanddatas), (*conshdlrdata)->sallconsanddatas ) );

   /* set hashmap and -table to NULL, mark them as uninitialized */
   (*conshdlrdata)->inithashmapandtable = FALSE;
   (*conshdlrdata)->hashtable = NULL;
   (*conshdlrdata)->hashtablesize = 0;
   (*conshdlrdata)->hashmap = NULL;
   (*conshdlrdata)->hashmapsize = 0;

   /* for constraint names count number of created constraints */
   (*conshdlrdata)->nlinconss = 0;

   /* initializes how many consanddata objects are used by original constraints */
   (*conshdlrdata)->noriguses = 0;

   return SCIP_OKAY;
}


/** frees constraint handler data for pseudo boolean constraint handler */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);
   assert((*conshdlrdata)->nallconsanddatas == 0);

   /* free hash table if necessary */
   if( (*conshdlrdata)->inithashmapandtable )
   {
      SCIPhashmapFree(&((*conshdlrdata)->hashmap));
      (*conshdlrdata)->hashmapsize = 0;
      SCIPhashtableFree(&((*conshdlrdata)->hashtable));
      (*conshdlrdata)->hashtablesize = 0;
   }
   else
   {
      assert((*conshdlrdata)->hashmap == NULL);
      assert((*conshdlrdata)->hashtable == NULL);
   }
   (*conshdlrdata)->inithashmapandtable = FALSE;

   /* clear array for all consanddata objects */
   SCIPfreeBlockMemoryArray(scip, &((*conshdlrdata)->allconsanddatas), (*conshdlrdata)->sallconsanddatas );

   (*conshdlrdata)->allconsanddatas = NULL;
   (*conshdlrdata)->nallconsanddatas = 0;
   (*conshdlrdata)->sallconsanddatas = 0;

   SCIPfreeBlockMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** gets number of variables in linear constraint */
static
SCIP_RETCODE getLinearConsNVars(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< linear constraint */
   SCIP_LINEARCONSTYPE const constype,       /**< linear constraint type */
   int*const             nvars               /**< pointer to store number variables of linear constraint */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nvars != NULL);

   /* determine for each special linear constranit all variables and coefficients */
   switch( constype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      *nvars = SCIPgetNVarsLinear(scip, cons);
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
      *nvars = SCIPgetNVarsLogicor(scip, cons);
      break;
   case SCIP_LINEARCONSTYPE_KNAPSACK:
      *nvars = SCIPgetNVarsKnapsack(scip, cons);
      break;
   case SCIP_LINEARCONSTYPE_SETPPC:
      *nvars = SCIPgetNVarsSetppc(scip, cons);
      break;
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
      *nvars = SCIPgetNVarsEQKnapsack(scip, cons);
      break;
#endif
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}


/** gets sides of linear constraint */
static
SCIP_RETCODE getLinearConsSides(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< linear constraint */
   SCIP_LINEARCONSTYPE const constype,       /**< linear constraint type */
   SCIP_Real*const       lhs,                /**< pointer to store left hand side of linear constraint */
   SCIP_Real*const       rhs                 /**< pointer to store right hand side of linear constraint */
   )
{
   SCIP_SETPPCTYPE type;

   switch( constype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      *lhs = SCIPgetLhsLinear(scip, cons);
      *rhs = SCIPgetRhsLinear(scip, cons);
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
      *lhs = 1.0;
      *rhs = SCIPinfinity(scip);
      break;
   case SCIP_LINEARCONSTYPE_KNAPSACK:
      *lhs = -SCIPinfinity(scip);
      *rhs = SCIPgetCapacityKnapsack(scip, cons);
      break;
   case SCIP_LINEARCONSTYPE_SETPPC:
      type = SCIPgetTypeSetppc(scip, cons);

      switch( type )
      {
      case SCIP_SETPPCTYPE_PARTITIONING:
         *lhs = 1.0;
         *rhs = 1.0;
         break;
      case SCIP_SETPPCTYPE_PACKING:
         *lhs = -SCIPinfinity(scip);
         *rhs = 1.0;
         break;
      case SCIP_SETPPCTYPE_COVERING:
         *lhs = 1.0;
         *rhs = SCIPinfinity(scip);
         break;
      default:
         SCIPerrorMessage("unknown setppc type\n");
         return SCIP_INVALIDDATA;
      }
      break;
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
      *lhs = SCIPgetCapacityEQKnapsack(scip, cons);
      *rhs = *lhs;
      break;
#endif
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** gets variables and coefficients of linear constraint */
static
SCIP_RETCODE getLinearConsVarsData(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< linear constraint */
   SCIP_LINEARCONSTYPE const constype,       /**< linear constraint type */
   SCIP_VAR**const       vars,               /**< array to store sorted (after indices) variables of linear constraint */
   SCIP_Real*const       coefs,              /**< array to store coefficient of linear constraint, or NULL */
   int*const             nvars               /**< pointer to store number variables of linear constraint */
   )
{
   SCIP_VAR** linvars;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(vars != NULL);
   assert(nvars != NULL);

   /* determine for each special linear constrait all variables and coefficients */
   switch( constype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
   {
      SCIP_Real* lincoefs;

      *nvars = SCIPgetNVarsLinear(scip, cons);
      linvars = SCIPgetVarsLinear(scip, cons);

      if( coefs != NULL )
      {
         lincoefs = SCIPgetValsLinear(scip, cons);

         for( v = 0; v < *nvars; ++v )
         {
            vars[v] = linvars[v];
            coefs[v] = lincoefs[v];
         }
      }
      else
      {
         for( v = 0; v < *nvars; ++v )
            vars[v] = linvars[v];
      }

      break;
   }
   case SCIP_LINEARCONSTYPE_LOGICOR:
      *nvars = SCIPgetNVarsLogicor(scip, cons);
      linvars = SCIPgetVarsLogicor(scip, cons);
      assert( linvars != NULL );

      if( coefs != NULL )
      {
         for( v = 0; v < *nvars; ++v )
         {
            vars[v] = linvars[v];
            coefs[v] = 1.0;
         }
      }
      else
      {
         for( v = 0; v < *nvars; ++v )
            vars[v] = linvars[v];
      }

      break;
   case SCIP_LINEARCONSTYPE_KNAPSACK:
   {
      SCIP_Longint* weights;

      *nvars = SCIPgetNVarsKnapsack(scip, cons);
      linvars = SCIPgetVarsKnapsack(scip, cons);
      assert( linvars != NULL );

      if( coefs != NULL )
      {
         weights = SCIPgetWeightsKnapsack(scip, cons);

         for( v = 0; v < *nvars; ++v )
         {
            vars[v] = linvars[v];
            coefs[v] = (SCIP_Real) weights[v];
         }
      }
      else
      {
         for( v = 0; v < *nvars; ++v )
            vars[v] = linvars[v];
      }

      break;
   }
   case SCIP_LINEARCONSTYPE_SETPPC:
      *nvars = SCIPgetNVarsSetppc(scip, cons);
      linvars = SCIPgetVarsSetppc(scip, cons);
      assert( linvars != NULL );

      if( coefs != NULL )
      {
	 for( v = 0; v < *nvars; ++v )
	 {
	    vars[v] = linvars[v];
	    coefs[v] = 1.0;
	 }
      }
      else
      {
	 for( v = 0; v < *nvars; ++v )
	    vars[v] = linvars[v];
      }

      break;
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
   {
      SCIP_Longint* weights;

      *nvars = SCIPgetNVarsEQKnapsack(scip, cons);
      linvars = SCIPgetVarsEQKnapsack(scip, cons);
      assert( linvars != NULL );

      if( coefs != NULL )
      {
	 weights = SCIPgetWeightsEQKnapsack(scip, cons);

	 for( v = 0; v < *nvars; ++v )
	 {
	    vars[v] = linvars[v];
	    coefs[v] = (SCIP_Real) weights[v];
	 }
      }
      else
      {
	 for( v = 0; v < *nvars; ++v )
	    vars[v] = linvars[v];
      }

      break;
   }
#endif
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   /* sort variables after indices */
   if( coefs != NULL )
   {
      SCIPsortPtrReal((void**)vars, coefs, SCIPvarComp, *nvars);
   }
   else
   {
      SCIPsortPtr((void**)vars, SCIPvarComp, *nvars);
   }

   return SCIP_OKAY;
}

/** calculate all not artificial linear variables and all artificial and-resultants which will be ordered like the
 *  'consanddatas' such that the and-resultant of the and-constraint is the and-resultant in the 'andress' array
 *  afterwards
 */
static
SCIP_RETCODE getLinVarsAndAndRess(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_VAR**const       vars,               /**< all variables of linear constraint */
   SCIP_Real*const       coefs,              /**< all coefficients of linear constraint, or NULL */
   int const             nvars,              /**< number of all variables of linear constraint */
   SCIP_VAR**const       linvars,            /**< array to store not and-resultant variables of linear constraint, or NULL */
   SCIP_Real*const       lincoefs,           /**< array to store coefficients of not and-resultant variables of linear
                                              *   constraint, or NULL */
   int*const             nlinvars,           /**< pointer to store number of not and-resultant variables, or NULL */
   SCIP_VAR**const       andress,            /**< array to store and-resultant variables of linear constraint, or NULL */
   SCIP_Real*const       andcoefs,           /**< array to store coefficients of and-resultant variables of linear
                                              *   constraint, or NULL */
   SCIP_Bool*const       andnegs,            /**< array to store negation status of and-resultant variables of linear
                                              *   constraint, or NULL */
   int*const             nandress            /**< pointer to store number of and-resultant variables, or NULL */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(vars != NULL);
   assert((linvars != NULL) == (nlinvars != NULL));
   assert((andress == NULL) || (nandress != NULL));
   assert((andcoefs != NULL) == (andnegs != NULL));
   assert((coefs != NULL) == ((lincoefs != NULL) || (andcoefs != NULL)));
   assert(linvars != NULL || andress != NULL);

   if( nlinvars != NULL )
      *nlinvars = 0;
   if( nandress != NULL )
      *nandress = 0;

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->hashmap != NULL);

   /* @note it is necessary that the linear constraint is merged (not needed for negated variables) and sorted after
    *       indices
    */

#ifndef NDEBUG
   /* check that old variables are sorted */
   for( v = nvars - 1; v > 0; --v )
      assert(SCIPvarGetIndex(vars[v]) >= SCIPvarGetIndex(vars[v - 1]));
#endif

   /* split variables into original and artificial variables */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Bool hashmapentryexists;
      SCIP_VAR* hashmapvar;

      assert(vars[v] != NULL);

      hashmapentryexists = SCIPhashmapExists(conshdlrdata->hashmap, (void*)(vars[v]));

      if( !hashmapentryexists && SCIPvarIsNegated(vars[v]) )
      {
         hashmapvar = SCIPvarGetNegationVar(vars[v]);
         hashmapentryexists = SCIPhashmapExists(conshdlrdata->hashmap, (void*)(hashmapvar));
      }
      else
         hashmapvar = vars[v];

      /* if and resultant is not a resultant anymore (meaning the corresponding and-constraint was deleted/upgraded),
       * correct the flag and count this variable as normal linear variable
       */
      if( hashmapentryexists )
      {
	 if( !SCIPconsIsOriginal(cons) )
	 {
            CONSANDDATA* consanddata = (CONSANDDATA*) SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)(hashmapvar));
            assert(consanddata != NULL);

	    hashmapentryexists = (consanddata->istransformed);

	    if( hashmapentryexists )
	    {
	       assert(consanddata->cons != NULL);
	       hashmapentryexists = !SCIPconsIsDeleted(consanddata->cons);
	    }
	 }
      }

      if( !hashmapentryexists && linvars != NULL )
      {
         assert(nlinvars != NULL);

         linvars[*nlinvars] = vars[v];
	 if( lincoefs != NULL )
	 {
	    assert(coefs != NULL);
	    lincoefs[*nlinvars] = coefs[v];
	 }
         ++(*nlinvars);
      }
      else if( hashmapentryexists && nandress != NULL )
      {
         if( andress != NULL )
         {
            andress[*nandress] = hashmapvar;

            if( andcoefs != NULL )
            {
               assert(andnegs != NULL);
               assert(coefs != NULL);
               andcoefs[*nandress] = coefs[v];
               andnegs[*nandress] = (vars[v] != hashmapvar);
            }
         }
         ++(*nandress);
      }
   }

   /* @todo try to avoid sorting here */
   if( andress != NULL && nandress != NULL )
   {
      /* sort and resultants by their variable index */
      if( andcoefs != NULL )
      {
         assert(andnegs != NULL);
         SCIPsortPtrRealBool((void**)andress, andcoefs, andnegs, SCIPvarComp, *nandress);
      }
      else
      {
         SCIPsortPtr((void**)andress, SCIPvarComp, *nandress);
      }
   }

   return SCIP_OKAY;
}


#ifdef CHECK_CONSISTENCY
/** check constraint consistency */
static
void checkConsConsistency(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int nvars;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlinvars;
   SCIP_VAR** andress;
   SCIP_Real* andcoefs;
   SCIP_Bool* andnegs;
   int nandress;
   SCIP_Bool* alreadyfound;
   SCIP_VAR* res;
   int c;
   int v;
   SCIP_Real newlhs;
   SCIP_Real newrhs;

   assert(scip != NULL);
   assert(cons != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_FREETRANS )
      return;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check standard pointers and sizes */
   assert(consdata->lincons != NULL);
   assert(!SCIPconsIsDeleted(consdata->lincons));
   assert(consdata->linconstype > SCIP_LINEARCONSTYPE_INVALIDCONS);
   assert(consdata->consanddatas != NULL);
   assert(consdata->nconsanddatas > 0);
   assert(consdata->nconsanddatas <= consdata->sconsanddatas);

   /* get sides of linear constraint */
   SCIP_CALL_ABORT( getLinearConsSides(scip, consdata->lincons, consdata->linconstype, &newlhs, &newrhs) );
   assert(!SCIPisInfinity(scip, newlhs));
   assert(!SCIPisInfinity(scip, -newrhs));
   assert(SCIPisLE(scip, newlhs, newrhs));
   assert(SCIPisEQ(scip, newrhs, consdata->rhs) || SCIPisEQ(scip, newrhs, -consdata->lhs));
   assert(SCIPisEQ(scip, newlhs, consdata->lhs) || SCIPisEQ(scip, newlhs, -consdata->rhs));

   /* check number of linear variables */
   SCIP_CALL_ABORT( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );
   assert(nvars == consdata->nlinvars + consdata->nconsanddatas);

   /* get temporary memory */
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &coefs, nvars) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &linvars, nvars) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &lincoefs, nvars) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &andress, nvars) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &andcoefs, nvars) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &andnegs, nvars) );
   SCIP_CALL_ABORT( SCIPallocClearBufferArray(scip, &alreadyfound, nvars) );

   /* get variables and coefficients */
   SCIP_CALL_ABORT( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );
   assert(nvars == 0 || (coefs != NULL));

   /* calculate all not artificial linear variables and all artificial and-resultants */
   SCIP_CALL_ABORT( getLinVarsAndAndRess(scip, cons, vars, coefs, nvars, linvars, lincoefs, &nlinvars,
         andress, andcoefs, andnegs, &nandress) );
   assert(nlinvars == consdata->nlinvars);
   assert(nandress == consdata->nconsanddatas);

   for( v = nandress - 1; v >= 0; --v )
   {
      SCIP_VAR* andresultant = andress[v];
      int nfound = 0;

      for( c = consdata->nconsanddatas - 1; c >= 0; --c )
      {
         assert(consdata->consanddatas[c] != NULL);
         if( consdata->consanddatas[c]->cons != NULL )
         {
            res = SCIPgetResultantAnd(scip, consdata->consanddatas[c]->cons);
            assert(res != NULL);

            if( res == andresultant && consdata->andnegs[c] == andnegs[v] && consdata->andcoefs[c] == andcoefs[v] )
            {
               /* resultant should be either active or a negated variable of an active one */
               assert(SCIPvarIsActive(res) || (SCIPvarIsNegated(res) && SCIPvarIsActive(SCIPvarGetNegationVar(res))));
               assert(!alreadyfound[c]);

               /* all and-resultants should be merged, so it is only allowed that each variable exists one time */
               alreadyfound[c] = TRUE;
               ++nfound;
               break;
            }
         }
      }
      assert(nfound == 1);
   }

   for( c = consdata->nconsanddatas - 1; c >= 0; --c )
   {
      assert(alreadyfound[c]);
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &alreadyfound);
   SCIPfreeBufferArray(scip, &andnegs);
   SCIPfreeBufferArray(scip, &andcoefs);
   SCIPfreeBufferArray(scip, &andress);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);
}
#else
#define checkConsConsistency(scip, cons) /**/
#endif


/** transforming transformed consanddata object back to original space, if an corresponding original constraint exists,
 *  also clearing all transformed data, i.e. releasing transformed variables
 */
static
SCIP_RETCODE transformToOrig(
   SCIP*const            scip,               /**< SCIP data structure */
   CONSANDDATA*          consanddata,        /**< consanddata object */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   SCIP_VAR** tmpvars;
   SCIP_Bool origdata;
   int ntmpvars;
   int v;

   assert(scip != NULL);
   assert(consanddata != NULL);
   assert(conshdlrdata != NULL);

   origdata = TRUE;

   tmpvars = consanddata->vars;
   ntmpvars = consanddata->nvars;

   /* release all transformed variables */
   for( v = ntmpvars - 1; v >= 0; --v )
   {
      assert(tmpvars[v] != NULL);
      if( SCIPvarIsTransformed(tmpvars[v]) )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &tmpvars[v]) );
         origdata = FALSE;
      }
   }

   tmpvars = consanddata->newvars;
   ntmpvars = consanddata->nnewvars;

   /* release all variables */
   for( v = ntmpvars - 1; v >= 0; --v )
   {
      assert(tmpvars[v] != NULL);
      if( SCIPvarIsTransformed(tmpvars[v]) )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &tmpvars[v]) );
         origdata = FALSE;
      }
   }

   /* reinstall original data */
   if( !origdata || consanddata->nvars == 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(consanddata->vars), consanddata->svars);
      SCIPfreeBlockMemoryArrayNull(scip, &(consanddata->newvars), consanddata->snewvars);

      consanddata->nuses = 0;
      consanddata->nvars = 0;
      consanddata->svars = 0;
      consanddata->nnewvars = 0;
      consanddata->snewvars = 0;
      consanddata->istransformed = FALSE;

      if( consanddata->noriguses > 0 )
      {
         assert(consanddata->origcons != NULL);
         assert(consanddata->isoriginal);

         assert(SCIPgetNVarsAnd(scip, consanddata->origcons) > 0);
         assert(SCIPgetVarsAnd(scip, consanddata->origcons) != NULL);
         consanddata->nvars = SCIPgetNVarsAnd(scip, consanddata->origcons);
         consanddata->svars = consanddata->nvars;

         if( consanddata->nvars > 0 )
         {
            SCIP_VAR** andvars = SCIPgetVarsAnd(scip, consanddata->origcons);

            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consanddata->vars), andvars, consanddata->nvars) );

            /* sort variables */
            SCIPsortPtr((void**)(consanddata->vars), SCIPvarComp, consanddata->nvars);
         }

         /* check that the hash map and tabkle are still having all information */
         if( conshdlrdata->inithashmapandtable )
         {
            assert(conshdlrdata->hashmap != NULL);
            assert(conshdlrdata->hashtable != NULL);
            assert(SCIPgetResultantAnd(scip, consanddata->origcons) != NULL);
            assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddata));
            assert(consanddata == (CONSANDDATA*)(SCIPhashtableRetrieve(conshdlrdata->hashtable, (void*)consanddata)));
            assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddata->origcons)));
            assert(consanddata == (CONSANDDATA*)(SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddata->origcons))));
         }
      }
      else
         assert(consanddata->origcons == NULL);
   }
   else
   {
      assert(consanddata->nuses == 0);
      assert(consanddata->nnewvars == 0);
      assert(consanddata->snewvars == 0);
      assert(consanddata->newvars == NULL);

      consanddata->istransformed = FALSE;

      if( consanddata->noriguses > 0 )
      {
         assert(consanddata->origcons != NULL);
         assert(consanddata->nvars > 0);
         assert(consanddata->svars > 0);
         assert(consanddata->vars != NULL);
         assert(consanddata->isoriginal);

         /* check that the hash map and tabkle are still having all information */
         if( conshdlrdata->inithashmapandtable )
         {
            assert(conshdlrdata->hashmap != NULL);
            assert(conshdlrdata->hashtable != NULL);
            assert(SCIPgetResultantAnd(scip, consanddata->origcons) != NULL);
            assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddata));
            assert(consanddata == (CONSANDDATA*)(SCIPhashtableRetrieve(conshdlrdata->hashtable, (void*)consanddata)));
            assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddata->origcons)));
            assert(consanddata == (CONSANDDATA*)(SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddata->origcons))));
         }
      }
   }

   return SCIP_OKAY;
}



/** creates a pseudo boolean constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*const   conshdlr,           /**< pseudoboolean constraint handler */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   SCIP_CONS*const       lincons,            /**< linear constraint with artificial and-resultants representing this pseudoboolean constraint */
   SCIP_LINEARCONSTYPE const linconstype,    /**< type of linear constraint */
   SCIP_CONS**const      andconss,           /**< array of and-constraints which occur in this pseudoboolean constraint */
   SCIP_Real*const       andcoefs,           /**< coefficients of and-constraints */
   SCIP_Bool*const       andnegs,            /**< negation status of and-constraints (or NULL, if no negated resultants) */
   int const             nandconss,          /**< number of and-constraints */
   SCIP_VAR*const        indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real const       weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool const       issoftcons,         /**< is this a soft constraint */
   SCIP_VAR* const       intvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             check,              /**< is the new constraint a check constraint? */
   SCIP_Bool             transforming        /**< are we called by CONSTRANS */
   )
{
   SCIP_Bool transformed;
   int nvars;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(consdata != NULL);
   assert(lincons != NULL && linconstype > SCIP_LINEARCONSTYPE_INVALIDCONS);
   assert(nandconss == 0 || (andconss != NULL && andcoefs != NULL));
   assert(!issoftcons || (!SCIPisZero(scip, weight) && indvar != NULL));

   /* adjust right hand side */
   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -rhs) )
      rhs = -SCIPinfinity(scip);

   /* adjust left hand side */
   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, lhs) )
      lhs = SCIPinfinity(scip);

   /* check left and right side */
   if( SCIPisGT(scip, lhs, rhs) )
   {
      SCIPerrorMessage("left hand side of pseudo boolean constraint greater than right hand side\n");
      SCIPerrorMessage(" -> lhs=%g, rhs=%g\n", lhs, rhs);
      return SCIP_INVALIDDATA;
   }

   transformed = SCIPisTransformed(scip);

   /* allocate memory for the constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   /* initialize the weights for soft constraints */
   (*consdata)->issoftcons = issoftcons;
   if( issoftcons )
   {
      (*consdata)->weight = weight;
      if( transformed )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, indvar, &((*consdata)->indvar)) );
      }
      else
         (*consdata)->indvar = indvar;
   }
   else
      (*consdata)->indvar = NULL;

   /* copy artificial integer variable if it exist */
   if( intvar != NULL )
   {
      if( transformed )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, intvar, &((*consdata)->intvar)) );
      }
      else
         (*consdata)->intvar = intvar;
   }
   else
      (*consdata)->intvar = NULL;

   /* copy linear constraint */
   (*consdata)->lincons = lincons;
   (*consdata)->linconstype = linconstype;

   /* get transformed linear constraint and capture it if necessary */
   if( transforming )
   {
      /* do not capture the and constraint when scip is in transformed mode; this automatically happens in
       * SCIPtransformCons()
       */
      SCIP_CALL( SCIPtransformCons(scip, (*consdata)->lincons, &((*consdata)->lincons)) );
      assert((*consdata)->lincons != NULL);
   }

   if( transforming || transformed )
   {
      assert(SCIPconsIsTransformed((*consdata)->lincons));

      /* we want to check all necessary transformed linear constraints */
      SCIP_CALL( SCIPsetConsChecked(scip, (*consdata)->lincons, check) );
   }

   /* get number of non-linear terms in pseudoboolean constraint */
   SCIP_CALL( getLinearConsNVars(scip, (*consdata)->lincons, (*consdata)->linconstype, &nvars) );
   (*consdata)->nlinvars = nvars - nandconss;

   /* copy and-constraints */
   if( nandconss > 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_VAR** andress;
      int c;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->consanddatas), nandconss) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &((*consdata)->andcoefs), andcoefs, nandconss) );
      if( andnegs != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &((*consdata)->andnegs), andnegs, nandconss) );
      }
      else
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*consdata)->andnegs), nandconss) );
      }
      (*consdata)->nconsanddatas = nandconss;
      (*consdata)->sconsanddatas = nandconss;

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &andress, nandconss) );

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->hashmap != NULL);

      /* get all and-resultants for sorting */
      for( c = nandconss - 1; c >= 0; --c )
      {
         assert(andconss[c] != NULL);

         andress[c] = SCIPgetResultantAnd(scip, andconss[c]);
         assert(andress[c] != NULL);

         (*consdata)->consanddatas[c] = (CONSANDDATA*) SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)andress[c]);
         assert((*consdata)->consanddatas[c] != NULL);
         assert((*consdata)->consanddatas[c]->origcons == andconss[c] || (*consdata)->consanddatas[c]->cons == andconss[c]);

         if( transforming )
         {
            /* if we perform a new transformation, we need to capture the transformed constraint */
            if( (*consdata)->consanddatas[c]->origcons != NULL && (*consdata)->consanddatas[c]->cons == NULL )
            {
               SCIP_VAR** vars;
               int ncvars;
               int v;

               /* do not capture the and constraint when scip is in transformed mode; this automatically happens in
                * SCIPtransformCons()
                */
               SCIP_CALL( SCIPtransformCons(scip, (*consdata)->consanddatas[c]->origcons, &((*consdata)->consanddatas[c]->cons)) );
               assert((*consdata)->consanddatas[c]->cons != NULL);
               assert((*consdata)->consanddatas[c]->newvars == NULL);
               assert((*consdata)->consanddatas[c]->isoriginal);

               (*consdata)->consanddatas[c]->istransformed = TRUE;

               vars = (*consdata)->consanddatas[c]->vars;
               ncvars = (*consdata)->consanddatas[c]->nvars;
               assert(vars != NULL || ncvars == 0);

               /* get transformed variables */
               SCIP_CALL( SCIPgetTransformedVars(scip, ncvars, vars, vars) );

               /* resort variables in transformed problem, because the order might change while tranforming */
               SCIPsortPtr((void**)vars, SCIPvarComp, ncvars);

               /* capture all transformed variables */
               for( v = ncvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPcaptureVar(scip, vars[v]) ); /*lint !e613*/
               }
            }
            else if( (*consdata)->consanddatas[c]->cons != NULL )
               assert((*consdata)->consanddatas[c]->istransformed);

            ++((*consdata)->consanddatas[c]->nuses);
         }
         else if( transformed )
         {
            assert((*consdata)->consanddatas[c]->cons == andconss[c]);
            assert(SCIPconsIsTransformed(andconss[c]));
            assert((*consdata)->consanddatas[c]->istransformed);
         }
      }

      /* sort and-constraints after indices of corresponding and-resultants */
      SCIPsortPtrPtrRealBool((void**)andress, (void**)((*consdata)->consanddatas), (*consdata)->andcoefs, (*consdata)->andnegs, SCIPvarComp, nandconss);

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &andress);
   }
   else
   {
      (*consdata)->consanddatas = NULL;
      (*consdata)->andcoefs = NULL;
      (*consdata)->andnegs = NULL;
      (*consdata)->nconsanddatas = 0;
      (*consdata)->sconsanddatas = 0;
   }

   /* copy left and right hand side */
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;

   (*consdata)->changed = TRUE;
   (*consdata)->propagated = FALSE;
   (*consdata)->presolved = FALSE;
   (*consdata)->cliquesadded = FALSE;
   (*consdata)->upgradetried = TRUE;

   /* count number of used consanddata objects in original problem */
   if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      conshdlrdata->noriguses += (*consdata)->nconsanddatas;
   }

   return SCIP_OKAY;
}

/** free a pseudo boolean constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   SCIP_Bool             isorig,             /**< are we freeing an original constraint? */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   CONSANDDATA** consanddatas;
   int nconsanddatas;
   int c;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->nconsanddatas == 0 || (*consdata)->consanddatas != NULL);
   assert(conshdlrdata != NULL);

   /* release linear constraint */
   if( (*consdata)->lincons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &((*consdata)->lincons)) );
   }

   nconsanddatas = (*consdata)->nconsanddatas;
   consanddatas = (*consdata)->consanddatas;

   /* count down uses and if necessary release constraints and delete data from hashtable and -map */
   for( c = nconsanddatas - 1; c >= 0; --c )
   {
      assert((consanddatas[c]->origcons == NULL) == (consanddatas[c]->noriguses == 0));
      assert((consanddatas[c]->cons == NULL) == (consanddatas[c]->nuses == 0));
      assert(consanddatas[c]->nuses >= 0);
      assert(consanddatas[c]->noriguses >= 0);
      assert(isorig ? consanddatas[c]->cons == NULL : TRUE);

      /* are we deleteing a transformed constraint */
      if( !isorig && consanddatas[c]->cons != NULL )
      {
         assert(!SCIPconsIsOriginal(consanddatas[c]->cons));

         --(consanddatas[c]->nuses);

         /* if the consanddata is not used anymore, release the constraint and clear the hashmap- and table */
         if( consanddatas[c]->nuses == 0 )
         {
            if( conshdlrdata->inithashmapandtable )
            {
               assert(conshdlrdata->hashmap != NULL);
               assert(conshdlrdata->hashtable != NULL);

               /* remove consanddata from hashtable, if it existed only in transformed space */
               if( consanddatas[c]->origcons == NULL )
               {
                  assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddatas[c]));
                  SCIP_CALL( SCIPhashtableRemove(conshdlrdata->hashtable, (void*)consanddatas[c]) );
               }
               assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddatas[c]->cons)));
               SCIP_CALL( SCIPhashmapRemove(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddatas[c]->cons)) );
            }

            SCIP_CALL( SCIPreleaseCons(scip, &(consanddatas[c]->cons)) );

            /* if the consanddata object was only used in transformed space, delete the memory block */
            if( consanddatas[c]->origcons == NULL )
            {
               int d;

               assert(conshdlrdata->nallconsanddatas > 0);

               for( d = conshdlrdata->nallconsanddatas - 1; d >= 0; --d )
               {
                  if( conshdlrdata->allconsanddatas[d] == consanddatas[c] )
                  {
                     --conshdlrdata->nallconsanddatas;

                     SCIPfreeBlockMemory(scip, &(conshdlrdata->allconsanddatas[d])); /*lint !e866*/

                     conshdlrdata->allconsanddatas[d] = conshdlrdata->allconsanddatas[conshdlrdata->nallconsanddatas];
                     break;
                  }
               }
               assert(d >= 0);
               continue;
            }
         }
      }
      /* are we deleteing an original constraint */
      else if( isorig && consanddatas[c]->origcons != NULL )
      {
         assert(SCIPconsIsOriginal(consanddatas[c]->origcons));
         assert(consanddatas[c]->nuses == 0);
         assert(consanddatas[c]->nnewvars == 0);
         assert(consanddatas[c]->snewvars == 0);
         assert(consanddatas[c]->newvars == NULL);

         --(consanddatas[c]->noriguses);

         /* if the consanddata is not used anymore, release the constraint and clear the hashmap- and table */
         if( consanddatas[c]->noriguses == 0 )
         {
            int d;

            if( conshdlrdata->inithashmapandtable )
            {
               assert(conshdlrdata->hashmap != NULL);
               assert(conshdlrdata->hashtable != NULL);

               assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddatas[c]));
               SCIP_CALL( SCIPhashtableRemove(conshdlrdata->hashtable, (void*)consanddatas[c]) );

               assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddatas[c]->origcons)));
               SCIP_CALL( SCIPhashmapRemove(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddatas[c]->origcons)) );
            }

            if( consanddatas[c]->vars != NULL )
            {
               assert(consanddatas[c]->nvars > 0);
               assert(consanddatas[c]->svars > 0);
               assert(consanddatas[c]->svars >= consanddatas[c]->nvars);

               SCIPfreeBlockMemoryArrayNull(scip, &(consanddatas[c]->vars), consanddatas[c]->svars);
               consanddatas[c]->nvars = 0;
               consanddatas[c]->svars = 0;
            }
            else
            {
               assert(consanddatas[c]->nvars == 0);
               assert(consanddatas[c]->svars == 0);
            }

            SCIP_CALL( SCIPreleaseCons(scip, &(consanddatas[c]->origcons)) );
            assert(consanddatas[c]->origcons == NULL);

            /* delete consanddata object */
            assert(conshdlrdata->nallconsanddatas > 0);
            for( d = conshdlrdata->nallconsanddatas - 1; d >= 0; --d )
            {
               if( conshdlrdata->allconsanddatas[d] == consanddatas[c] )
               {
                  --conshdlrdata->nallconsanddatas;

                  SCIPfreeBlockMemory(scip, &(conshdlrdata->allconsanddatas[d])); /*lint !e866*/

                  conshdlrdata->allconsanddatas[d] = conshdlrdata->allconsanddatas[conshdlrdata->nallconsanddatas];
                  break;
               }
            }
            assert(d >= 0);

            continue;
         }
      }
      else
      {
         assert(!consanddatas[c]->istransformed);
         assert(consanddatas[c]->cons == NULL);
      }

      /* clear and remove capture of transformed consanddata */
      if( consanddatas[c]->nuses == 0 && consanddatas[c]->istransformed )
      {
         SCIP_CALL( transformToOrig(scip, consanddatas[c], conshdlrdata) );
      }
#ifndef NDEBUG
      else if( consanddatas[c]->nuses == 0 )
      {
         SCIP_VAR** tmpvars;
         int ntmpvars;
         int v;

         assert(consanddatas[c]->nnewvars == 0);
         assert(consanddatas[c]->snewvars == 0);
         assert(consanddatas[c]->newvars == NULL);

         tmpvars = consanddatas[c]->vars;
         ntmpvars = consanddatas[c]->nvars;

         /* release all variables */
         for( v = ntmpvars - 1; v >= 0; --v )
         {
            assert(tmpvars[v] != NULL);
            assert(SCIPvarIsOriginal(tmpvars[v]));
         }
      }
#endif

      /* restore original data */
      if( !consanddatas[c]->istransformed && consanddatas[c]->noriguses > 0 )
      {
         assert(consanddatas[c]->origcons != NULL);
         assert(consanddatas[c]->nuses == 0);
         assert(consanddatas[c]->nnewvars == 0);
         assert(consanddatas[c]->snewvars == 0);
         assert(consanddatas[c]->newvars == NULL);
         assert(consanddatas[c]->nvars > 0);
         assert(consanddatas[c]->svars > 0);
         assert(consanddatas[c]->svars >= consanddatas[c]->nvars);
         assert(consanddatas[c]->vars != NULL);
         assert(consanddatas[c]->isoriginal);

         assert(consanddatas[c]->nvars == SCIPgetNVarsAnd(scip, consanddatas[c]->origcons));
         assert(SCIPgetVarsAnd(scip, consanddatas[c]->origcons) != NULL);

         /* check that the hash map and tabkle are still having all information */
         if( conshdlrdata->inithashmapandtable )
         {
            assert(conshdlrdata->hashmap != NULL);
            assert(conshdlrdata->hashtable != NULL);
            assert(SCIPgetResultantAnd(scip, consanddatas[c]->origcons) != NULL);
            assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddatas[c]));
            assert(consanddatas[c] == (CONSANDDATA*)(SCIPhashtableRetrieve(conshdlrdata->hashtable, (void*)consanddatas[c])));
            assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddatas[c]->origcons)));
            assert(consanddatas[c] == (CONSANDDATA*)(SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddatas[c]->origcons))));
         }
      }
   }

   /* free array of and-constraints */
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->andnegs), (*consdata)->sconsanddatas);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->andcoefs), (*consdata)->sconsanddatas);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->consanddatas), (*consdata)->sconsanddatas);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** check the locks of an AND resultant and removes it from all global structures if the resultant is not locked anymore */
static
SCIP_RETCODE checkLocksAndRes(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR*             res                 /**< resultant of AND constraint */
   )
{
   assert(scip != NULL);
   assert(res != NULL);

   /* the resultant has no locks left and might be dual fixed now, we need to delete all its cliques */
   if( SCIPvarIsActive(res) && SCIPvarGetNLocksDown(res) == 0 && SCIPvarGetNLocksUp(res) == 0
      && SCIPgetStage(scip) < SCIP_STAGE_FREETRANS )
   {
      SCIP_CALL( SCIPremoveVarFromGlobalStructures(scip, res) );
   }

   return SCIP_OKAY;
}

/** installs rounding locks for the given and-constraint associated with given coefficient */
static
SCIP_RETCODE lockRoundingAndCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   CONSANDDATA*const     consanddata,        /**< CONSANDDATA object for which we want to add the locks */
   SCIP_Real const       coef,               /**< coefficient which led to old locks */
   SCIP_Real const       lhs,                /**< left hand side */
   SCIP_Real const       rhs                 /**< right hand side */
   )
{
   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR* res;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(consanddata != NULL);
   assert(!SCIPisInfinity(scip, coef) && !SCIPisInfinity(scip, -coef));
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(SCIPisLE(scip, lhs, rhs));

   /* choose correct variable array to add locks for, we only add locks for now valid variables */
   if( consanddata->nnewvars > 0 )
   {
      vars = consanddata->newvars;
      nvars = consanddata->nnewvars;
   }
   else
   {
      vars = consanddata->vars;
      nvars = consanddata->nvars;
   }

   res = SCIPgetResultantAnd(scip, consanddata->cons);
   assert(nvars == 0 || (vars != NULL && res != NULL));

   /* check which sites are infinity */
   haslhs = !SCIPisInfinity(scip, -lhs);
   hasrhs = !SCIPisInfinity(scip, rhs);

   if( SCIPconsIsLocked(cons) )
   {
      /* locking variables */
      if( SCIPisPositive(scip, coef) )
      {
         for( v = nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, haslhs, hasrhs) );
         }
      }
      else
      {
         for( v = nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, hasrhs, haslhs) );
         }
      }
      SCIP_CALL( SCIPlockVarCons(scip, res, cons, TRUE, TRUE) );
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given and-constraint associated with given coefficient */
static
SCIP_RETCODE unlockRoundingAndCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   CONSANDDATA*const     consanddata,        /**< CONSANDDATA object for which we want to delete the locks */
   SCIP_Real const       coef,               /**< coefficient which led to old locks */
   SCIP_Real const       lhs,                /**< left hand side which led to old locks */
   SCIP_Real const       rhs                 /**< right hand side which led to old locks */
   )
{
   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR* res;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(consanddata != NULL);
   assert(!SCIPisInfinity(scip, coef) && !SCIPisInfinity(scip, -coef));
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(SCIPisLE(scip, lhs, rhs));

   vars = consanddata->vars;
   nvars = consanddata->nvars;

   if( consanddata->cons != NULL )
      res = SCIPgetResultantAnd(scip, consanddata->cons);
   else
      res = NULL;
   assert(nvars == 0 || vars != NULL);

   /* check which sites are infinity */
   haslhs = !SCIPisInfinity(scip, -lhs);
   hasrhs = !SCIPisInfinity(scip, rhs);

   if( SCIPconsIsLocked(cons) )
   {
      /* unlock variables */
      if( SCIPisPositive(scip, coef) )
      {
         for( v = nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, haslhs, hasrhs) );
         }
      }
      else
      {
         for( v = nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, hasrhs, haslhs) );
         }
      }

      if( res != NULL )
      {
         SCIP_CALL( SCIPunlockVarCons(scip, res, cons, TRUE, TRUE) );

         SCIP_CALL( checkLocksAndRes(scip, res) );
      }
   }

   return SCIP_OKAY;
}

/** prints pseudoboolean constraint in CIP format to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   FILE*const            file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int nvars;
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlinvars;
   int v;

   SCIP_VAR** andress;
   SCIP_Real* andcoefs;
   SCIP_Bool* andnegs;
   int nandress;

   SCIP_Bool printed;

   assert(scip != NULL);
   assert(cons != NULL);

#ifdef WITHEQKNAPSACK
   if( SCIPconsIsDeleted(cons) )
      return SCIP_OKAY;
#endif

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lincons != NULL);
   /* more than one and-constraint is needed, otherwise this pseudoboolean constraint should be upgraded to a linear constraint */
   assert(consdata->nconsanddatas >= 0);

   /* gets number of variables in linear constraint */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andress, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andnegs, nvars) );

   /* get sides of linear constraint */
   SCIP_CALL( getLinearConsSides(scip, consdata->lincons, consdata->linconstype, &lhs, &rhs) );
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(SCIPisLE(scip, lhs, rhs));

   /* get variables and coefficient of linear constraint */
   SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );
   assert(nvars == 0 || (coefs != NULL));

   /* calculate all not artificial linear variables and all artificial and-resultants which will be ordered like the
    * 'consanddatas' such that the and-resultant of the and-constraint is the and-resultant in the 'andress' array
    * afterwards
    */
   SCIP_CALL( getLinVarsAndAndRess(scip, cons, vars, coefs, nvars, linvars, lincoefs, &nlinvars,
         andress, andcoefs, andnegs, &nandress) );
   assert(consdata->nconsanddatas == nandress);

   /* number of variables should be consistent, number of 'real' linear variables plus number of and-constraints should
    * have to be equal to the number of variables in the linear constraint
    */
   assert(consdata->nlinvars + consdata->nconsanddatas == nvars);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", lhs);

   printed = FALSE;

   /* print coefficients and variables */
   if( nlinvars > 0)
   {
      printed= TRUE;

      /* print linear part of constraint */
      SCIP_CALL( SCIPwriteVarsLinearsum(scip, file, linvars, lincoefs, nlinvars, TRUE) );
   }

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->hashmap != NULL);

   /* print all non-linear terms */
   for( v = nandress - 1; v >= 0; --v )
   {
      CONSANDDATA* consanddata;
      SCIP_CONS* andcons;
      SCIP_VAR** andvars;
      int nandvars;

      if( !SCIPconsIsOriginal(cons) )
      {
         /* if the and resultant was fixed we print a constant */
         if( SCIPvarGetLbLocal(andress[v]) > 0.5 || SCIPvarGetUbLocal(andress[v]) < 0.5 )
         {
            if( SCIPvarGetLbGlobal(andress[v]) > 0.5 )
            {
               printed = TRUE;
               SCIPinfoMessage(scip, file, " %+.15g ", andcoefs[v] * SCIPvarGetLbGlobal(andress[v]));
            }
            continue;
         }
         else if( SCIPvarGetStatus(andress[v]) == SCIP_VARSTATUS_AGGREGATED )
         {
            SCIP_VAR* aggrvar;
            SCIP_Bool negated;

            SCIP_CALL( SCIPgetBinvarRepresentative(scip, andress[v], &aggrvar, &negated) );
            assert(aggrvar != NULL);
            assert(SCIPvarGetType(aggrvar) == SCIP_VARTYPE_BINARY);

            printed = TRUE;
            SCIPinfoMessage(scip, file, " %+.15g %s<%s>[B]", andcoefs[v], negated ? "~" : "", SCIPvarGetName(aggrvar));

            continue;
         }
      }

      consanddata = (CONSANDDATA*) SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)andress[v]);
      assert(consanddata != NULL);

      if( SCIPconsIsOriginal(cons) )
         andcons = consanddata->origcons;
      else
         andcons = consanddata->cons;
      assert(andcons != NULL);

      andvars = SCIPgetVarsAnd(scip, andcons);
      nandvars = SCIPgetNVarsAnd(scip, andcons);
      assert(nandvars == 0 || andvars != NULL);

      if( nandvars > 0 )
      {
         printed = TRUE;
         SCIPinfoMessage(scip, file, " %+.15g %s(", andcoefs[v], andnegs[v] ? "~" : "");

         /* @todo: better write new method SCIPwriteProduct */
         /* print variable list */
         SCIP_CALL( SCIPwriteVarsList(scip, file, andvars, nandvars, TRUE, '*') );

         SCIPinfoMessage(scip, file, ")");
      }
   }

   if( !printed )
   {
      SCIPinfoMessage(scip, file, " 0 ");
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &andnegs);
   SCIPfreeBufferArray(scip, &andcoefs);
   SCIPfreeBufferArray(scip, &andress);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   /* print right hand side */
   if( SCIPisEQ(scip, lhs, rhs) )
      SCIPinfoMessage(scip, file, "== %.15g", rhs);
   else if( !SCIPisInfinity(scip, rhs) )
      SCIPinfoMessage(scip, file, "<= %.15g", rhs);
   else if( !SCIPisInfinity(scip, -lhs) )
      SCIPinfoMessage(scip, file, ">= %.15g", lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}

/** creates and/or adds the resultant for a given term */
static
SCIP_RETCODE createAndAddAndCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*const   conshdlr,           /**< pseudoboolean constraint handler */
   SCIP_VAR**const       vars,               /**< array of variables to get and-constraints for */
   int const             nvars,              /**< number of variables to get and-constraints for */
   SCIP_Bool const       initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool const       enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant
                                              *   constraints. */
   SCIP_Bool const       check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant
                                              *   constraints. */
   SCIP_Bool const       local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching
                                              *   constraints. */
   SCIP_Bool const       modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE
                                              *   if pricing adds coefficients to this constraint. */
   SCIP_Bool const       dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool const       stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent
                                              *   node data. */
   SCIP_CONS**const      andcons             /**< pointer to store and-constraint */
   )
{
   CONSANDDATA* newdata;
   CONSANDDATA* tmpdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char name[SCIP_MAXSTRLEN];
   SCIP_Bool separate;
   SCIP_Bool propagate;
   SCIP_Bool removable;
   SCIP_Bool transformed;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(andcons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->hashtable != NULL);

   transformed = SCIPisTransformed(scip);

   /* allocate memory for a possible new consanddata object */
   SCIP_CALL( SCIPallocBlockMemory(scip, &newdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(newdata->vars), vars, nvars) );
   newdata->nvars = nvars;
   newdata->svars = nvars;
   newdata->newvars = NULL;
   newdata->nnewvars = 0;
   newdata->snewvars = 0;
   newdata->noriguses = 0;
   newdata->nuses = 0;
   newdata->istransformed = transformed;
   newdata->isoriginal = !transformed;
   newdata->cons = NULL;
   newdata->origcons = NULL;

   /* sort variables */
   SCIPsortPtr((void**)(newdata->vars), SCIPvarComp, nvars);

   /* get constraint from current hash table with same variables as cons0 */
   tmpdata = (CONSANDDATA*)(SCIPhashtableRetrieve(conshdlrdata->hashtable, (void*)newdata));

   /* if there is already the same and constraint created use this resultant */
   if( tmpdata != NULL )
   {
#ifndef NDEBUG
      SCIP_VAR* res;
#endif
      if( transformed )
      {
         assert(tmpdata->cons != NULL);
         *andcons = tmpdata->cons;

         assert(tmpdata->nuses > 0);
         /* increase usage of data object */
         ++(tmpdata->nuses);
      }
      else
      {
         assert(tmpdata->origcons != NULL);
         *andcons = tmpdata->origcons;

         assert(tmpdata->noriguses > 0);
         /* increase usage of data object */
         ++(tmpdata->noriguses);
      }
      assert(*andcons != NULL);

#ifndef NDEBUG
      res = SCIPgetResultantAnd(scip, *andcons);
      assert(res != NULL);

      /* check that we already have added this resultant to and-constraint entry */
      assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)res));
#endif
   }
   else
   {
      /* create new and-constraint */
      SCIP_CONS* newcons;
      SCIP_VAR* resultant;

      /* create auxiliary variable */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, ARTIFICIALVARNAMEPREFIX"%d", conshdlrdata->nallconsanddatas);
      SCIP_CALL( SCIPcreateVar(scip, &resultant, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
            TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );

#if 1 /* @todo: check whether we want to branch on artificial variables, the test results show that it is of advantage */
      /* change branching priority of artificial variable to -1 */
      SCIP_CALL( SCIPchgVarBranchPriority(scip, resultant, -1) );
#endif

      /* add auxiliary variable to the problem */
      SCIP_CALL( SCIPaddVar(scip, resultant) );

#if 0 /* does not work for since the value of artificial resultants must not be equal to the value computed by their
       * product, since these variables are irrelevant */
#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real val;
         SCIP_Real debugsolval;
         int v;

         for( v = nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPdebugGetSolVal(scip, vars[v], &val) );
            assert(SCIPisFeasZero(scip, val) || SCIPisFeasEQ(scip, val, 1.0));

            if( val < 0.5 )
               break;
         }
         val = ((val < 0.5) ? 0.0 : 1.0);

         SCIP_CALL( SCIPdebugGetSolVal(scip, resultant, &debugsolval) );
         if( (SCIPvarIsOriginal(resultant) || SCIPvarIsTransformedOrigvar(resultant)) && !SCIPisFeasEQ(scip, debugsolval, val) )
         {
            SCIPerrorMessage("computed solution value %g for resultant <%s> violates debug solution value %g\n", val, SCIPvarGetName(resultant), debugsolval);
            SCIPABORT();
            return SCIP_ERROR; /*lint !e527*/
         }
         else if( !SCIPvarIsOriginal(resultant) && !SCIPvarIsTransformedOrigvar(resultant) )
         {
            SCIP_CALL( SCIPdebugAddSolVal(scip, resultant, val) );
         }
      }
#endif
#endif

      SCIP_CALL( SCIPgetBoolParam(scip, "constraints/" CONSHDLR_NAME "/nlcseparate", &separate) );
      SCIP_CALL( SCIPgetBoolParam(scip, "constraints/" CONSHDLR_NAME "/nlcpropagate", &propagate) );
      SCIP_CALL( SCIPgetBoolParam(scip, "constraints/" CONSHDLR_NAME "/nlcremovable", &removable) );

      /* we do not want to check the and constraints, so the check flag will be FALSE */

      /* create and add "and" constraint for the multiplication of the binary variables */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "andcons_%d", conshdlrdata->nallconsanddatas);
      SCIP_CALL( SCIPcreateConsAnd(scip, &newcons, name, resultant, newdata->nvars, newdata->vars,
            initial, separate, enforce, check && FALSE, propagate,
            local, modifiable, dynamic, removable, stickingatnode) ); /*lint !e506*/
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIPdebugPrintCons(scip, newcons, NULL);

      /* force all deriving constraint from this and constraint to be checked and not removable */
      SCIP_CALL( SCIPchgAndConsCheckFlagWhenUpgr(scip, newcons, TRUE) );
      SCIP_CALL( SCIPchgAndConsRemovableFlagWhenUpgr(scip, newcons, TRUE) );

      *andcons = newcons;
      assert(*andcons != NULL);

      /* resize data for all and-constraints if necessary */
      if( conshdlrdata->nallconsanddatas == conshdlrdata->sallconsanddatas )
      {
         SCIP_CALL( SCIPensureBlockMemoryArray(scip, &(conshdlrdata->allconsanddatas), &(conshdlrdata->sallconsanddatas), SCIPcalcMemGrowSize(scip, conshdlrdata->sallconsanddatas + 1)) );
      }

      /* add new data object to global hash table */
      conshdlrdata->allconsanddatas[conshdlrdata->nallconsanddatas] = newdata;
      ++(conshdlrdata->nallconsanddatas);

      if( transformed )
      {
         int v;

         newdata->cons = newcons;
         SCIP_CALL( SCIPcaptureCons(scip, newdata->cons) );

         /* initialize usage of data object */
         newdata->nuses = 1;

         /* capture all variables */
         for( v = newdata->nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPcaptureVar(scip, newdata->vars[v]) ); /*lint !e613*/
         }
      }
      else
      {
         newdata->origcons = newcons;
         SCIP_CALL( SCIPcaptureCons(scip, newdata->origcons) );

         /* initialize usage of data object */
         newdata->noriguses = 1;
      }

      /* no such and-constraint in current hash table: insert the new object into hash table */
      SCIP_CALL( SCIPhashtableInsert(conshdlrdata->hashtable, (void*)newdata) );

      /* insert new mapping */
      assert(!SCIPhashmapExists(conshdlrdata->hashmap, (void*)resultant));
      SCIP_CALL( SCIPhashmapInsert(conshdlrdata->hashmap, (void*)resultant, (void*)newdata) );

      /* release and-resultant and -constraint */
      SCIP_CALL( SCIPreleaseVar(scip, &resultant) );
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

      return SCIP_OKAY;
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &(newdata->vars), newdata->svars);
   SCIPfreeBlockMemory(scip, &newdata);

   return SCIP_OKAY;
}

/** adds a term to the given pseudoboolean constraint */
static
SCIP_RETCODE addCoefTerm(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_VAR**const       vars,               /**< variables of the nonlinear term */
   int const             nvars,              /**< number of variables of the nonlinear term */
   SCIP_Real const       val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* andcons;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* res;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nvars == 0 || vars != NULL);

   if( nvars == 0 || SCIPisZero(scip, val) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   conshdlrdata =  SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create (and add) and-constraint */
   SCIP_CALL( createAndAddAndCons(scip, conshdlr, vars, nvars,
         SCIPconsIsInitial(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsStickingAtNode(cons),
         &andcons) );
   assert(andcons != NULL);

   /* ensure memory size */
   if( consdata->nconsanddatas == consdata->sconsanddatas )
   {
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, &(consdata->consanddatas), &(consdata->sconsanddatas), consdata->sconsanddatas + 1) );
   }

   res = SCIPgetResultantAnd(scip, andcons);
   assert(res != NULL);
   assert(SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)res) != NULL);

   consdata->consanddatas[consdata->nconsanddatas] = (CONSANDDATA*) SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)res);
   ++(consdata->nconsanddatas);

   /* add auxiliary variables to linear constraint */
   switch( consdata->linconstype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      SCIP_CALL( SCIPaddCoefLinear(scip, consdata->lincons, res, val) );
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
      if( !SCIPisEQ(scip, val, 1.0) )
         return SCIP_INVALIDDATA;

      SCIP_CALL( SCIPaddCoefLogicor(scip, consdata->lincons, res) );
      break;
   case SCIP_LINEARCONSTYPE_KNAPSACK:
      if( !SCIPisIntegral(scip, val) || !SCIPisPositive(scip, val) )
         return SCIP_INVALIDDATA;

      SCIP_CALL( SCIPaddCoefKnapsack(scip, consdata->lincons, res, (SCIP_Longint) val) );
      break;
   case SCIP_LINEARCONSTYPE_SETPPC:
      if( !SCIPisEQ(scip, val, 1.0) )
         return SCIP_INVALIDDATA;

      SCIP_CALL( SCIPaddCoefSetppc(scip, consdata->lincons, res) );
      break;
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
      if( !SCIPisIntegral(scip, val) || !SCIPisPositive(scip, val) )
         return SCIP_INVALIDDATA;

      SCIP_CALL( SCIPaddCoefEQKnapsack(scip, consdata->lincons, res, (SCIP_Longint) val) );
      break;
#endif
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   /* install rounding locks for all new variable */
   SCIP_CALL( lockRoundingAndCons(scip, cons, consdata->consanddatas[consdata->nconsanddatas - 1], val, consdata->lhs, consdata->rhs) );

   /* change flags */
   consdata->changed = TRUE;
   consdata->propagated = FALSE;
   consdata->presolved = FALSE;
   consdata->cliquesadded = FALSE;
   consdata->upgradetried = FALSE;

   return SCIP_OKAY;
}

/** changes left hand side of linear constraint */
static
SCIP_RETCODE chgLhsLinearCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< linear constraint */
   SCIP_LINEARCONSTYPE const constype,       /**< linear constraint type */
   SCIP_Real const       lhs                 /**< new left hand side of linear constraint */
   )
{
   switch( constype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      SCIP_CALL( SCIPchgLhsLinear(scip, cons, lhs) );
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
   case SCIP_LINEARCONSTYPE_KNAPSACK:
   case SCIP_LINEARCONSTYPE_SETPPC:
      SCIPerrorMessage("changing left hand side only allowed on standard lienar constraint \n");
      return SCIP_INVALIDDATA;
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
#endif
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes right hand side of linear constraint */
static
SCIP_RETCODE chgRhsLinearCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< linear constraint */
   SCIP_LINEARCONSTYPE const constype,       /**< linear constraint type */
   SCIP_Real const       rhs                 /**< new right hand side of linear constraint */
   )
{
   switch( constype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      SCIP_CALL( SCIPchgRhsLinear(scip, cons, rhs) );
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
   case SCIP_LINEARCONSTYPE_KNAPSACK:
   case SCIP_LINEARCONSTYPE_SETPPC:
      SCIPerrorMessage("changing left hand side only allowed on standard lienar constraint \n");
      return SCIP_INVALIDDATA;
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
#endif
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** sets left hand side of linear constraint */
static
SCIP_RETCODE chgLhs(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< linear constraint */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int nvars;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlinvars;
   SCIP_VAR** andress;
   SCIP_Real* andcoefs;
   SCIP_Bool* andnegs;
   int nandress;
   SCIP_Real oldlhs;
   SCIP_Real oldrhs;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, lhs));

   /* adjust value to not be smaller than -inf */
   if ( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get sides of linear constraint */
   SCIP_CALL( getLinearConsSides(scip, consdata->lincons, consdata->linconstype, &oldlhs, &oldrhs) );
   assert(!SCIPisInfinity(scip, oldlhs));
   assert(!SCIPisInfinity(scip, -oldrhs));
   assert(SCIPisLE(scip, oldlhs, oldrhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, oldlhs, lhs) )
      return SCIP_OKAY;

   /* gets number of variables in linear constraint */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andress, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andnegs, nvars) );

   /* get variables and coefficient of linear constraint */
   SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );
   assert(nvars == 0 || (coefs != NULL));

   /* calculate all not artificial linear variables and all artificial and-resultants which will be ordered like the
    * 'consanddatas' such that the and-resultant of the and-constraint is the and-resultant in the 'andress' array
    * afterwards
    */
   SCIP_CALL( getLinVarsAndAndRess(scip, cons, vars, coefs, nvars, linvars, lincoefs, &nlinvars, andress, andcoefs, andnegs, &nandress) );
   assert(consdata->nconsanddatas == nandress);

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      SCIP_VAR** andvars;
      int nandvars;
      SCIP_Real val;
      int v;
      int c;

      assert(SCIPconsIsTransformed(cons));

      if( SCIPisInfinity(scip, -oldlhs) && !SCIPisInfinity(scip, -lhs) )
      {
         /* non-linear part */
         for( c = consdata->nconsanddatas - 1; c >= 0; --c )
         {
            CONSANDDATA* consanddata;
            SCIP_CONS* andcons;

            consanddata = consdata->consanddatas[c];
            assert(consanddata != NULL);

            andcons = consanddata->cons;
            assert(andcons != NULL);

            andvars = SCIPgetVarsAnd(scip, andcons);
            nandvars = SCIPgetNVarsAnd(scip, andcons);
            val = andnegs[c] ? -andcoefs[c] : andcoefs[c];

            /* lock variables */
            if( SCIPisPositive(scip, val) )
            {
               for( v = nandvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPlockVarCons(scip, andvars[v], cons, TRUE, FALSE) );
               }
            }
            else
            {
               for( v = nandvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPlockVarCons(scip, andvars[v], cons, FALSE, TRUE) );
               }
            }
         }
      }
      else if( !SCIPisInfinity(scip, -oldlhs) && SCIPisInfinity(scip, -lhs) )
      {
         /* non-linear part */
         for( c = consdata->nconsanddatas - 1; c >= 0; --c )
         {
            CONSANDDATA* consanddata;
            SCIP_CONS* andcons;

            consanddata = consdata->consanddatas[c];
            assert(consanddata != NULL);

            andcons = consanddata->cons;
            assert(andcons != NULL);

            andvars = SCIPgetVarsAnd(scip, andcons);
            nandvars = SCIPgetNVarsAnd(scip, andcons);
            val = andnegs[c] ? -andcoefs[c] : andcoefs[c];

            /* lock variables */
            if( SCIPisPositive(scip, val) )
            {
               for( v = nandvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPunlockVarCons(scip, andvars[v], cons, TRUE, FALSE) );
               }
            }
            else
            {
               for( v = nandvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPunlockVarCons(scip, andvars[v], cons, FALSE, TRUE) );
               }
            }
         }
      }
   }

   /* check whether the left hand side is increased, if and only if that's the case we maybe can propagate, tighten and add more cliques */
   if( SCIPisLT(scip, oldlhs, lhs) )
   {
      consdata->propagated = FALSE;
   }

   /* set new left hand side and update constraint data */
   SCIP_CALL( chgLhsLinearCons(scip, consdata->lincons, consdata->linconstype, lhs) );
   consdata->lhs = lhs;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &andnegs);
   SCIPfreeBufferArray(scip, &andcoefs);
   SCIPfreeBufferArray(scip, &andress);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** sets right hand side of pseudoboolean constraint */
static
SCIP_RETCODE chgRhs(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< linear constraint */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int nvars;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlinvars;
   SCIP_VAR** andress;
   SCIP_Real* andcoefs;
   SCIP_Bool* andnegs;
   int nandress;
   SCIP_Real oldlhs;
   SCIP_Real oldrhs;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, -rhs));

   /* adjust value to not be larger than inf */
   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get sides of linear constraint */
   SCIP_CALL( getLinearConsSides(scip, consdata->lincons, consdata->linconstype, &oldlhs, &oldrhs) );
   assert(!SCIPisInfinity(scip, oldlhs));
   assert(!SCIPisInfinity(scip, -oldrhs));
   assert(SCIPisLE(scip, oldlhs, oldrhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, oldrhs, rhs) )
      return SCIP_OKAY;

   /* gets number of variables in linear constraint */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andress, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andnegs, nvars) );

   /* get variables and coefficient of linear constraint */
   SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );
   assert(nvars == 0 || (coefs != NULL));

   /* calculate all not artificial linear variables and all artificial and-resultants which will be ordered like the
    * 'consanddatas' such that the and-resultant of the and-constraint is the and-resultant in the 'andress' array
    * afterwards
    */
   SCIP_CALL( getLinVarsAndAndRess(scip, cons, vars, coefs, nvars, linvars, lincoefs, &nlinvars, andress, andcoefs, andnegs, &nandress) );
   assert(consdata->nconsanddatas == nandress);

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      SCIP_VAR** andvars;
      int nandvars;
      SCIP_Real val;
      int v;
      int c;

      assert(SCIPconsIsTransformed(cons));

      if( SCIPisInfinity(scip, oldrhs) && !SCIPisInfinity(scip, rhs) )
      {
         /* non-linear part */
         for( c = consdata->nconsanddatas - 1; c >= 0; --c )
         {
            CONSANDDATA* consanddata;
            SCIP_CONS* andcons;

            consanddata = consdata->consanddatas[c];
            assert(consanddata != NULL);

            andcons = consanddata->cons;
            assert(andcons != NULL);

            andvars = SCIPgetVarsAnd(scip, andcons);
            nandvars = SCIPgetNVarsAnd(scip, andcons);
            val = andnegs[c] ? -andcoefs[c] : andcoefs[c];

            /* lock variables */
            if( SCIPisPositive(scip, val) )
            {
               for( v = nandvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPlockVarCons(scip, andvars[v], cons, FALSE, TRUE) );
               }
            }
            else
            {
               for( v = nandvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPlockVarCons(scip, andvars[v], cons, TRUE, FALSE) );
               }
            }
         }
      }
      else if( !SCIPisInfinity(scip, oldrhs) && SCIPisInfinity(scip, rhs) )
      {
         /* non-linear part */
         for( c = consdata->nconsanddatas - 1; c >= 0; --c )
         {
            CONSANDDATA* consanddata;
            SCIP_CONS* andcons;

            consanddata = consdata->consanddatas[c];
            assert(consanddata != NULL);

            andcons = consanddata->cons;
            assert(andcons != NULL);

            andvars = SCIPgetVarsAnd(scip, andcons);
            nandvars = SCIPgetNVarsAnd(scip, andcons);
            val = andnegs[c] ? -andcoefs[c] : andcoefs[c];

            /* lock variables */
            if( SCIPisPositive(scip, val) )
            {
               for( v = nandvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPunlockVarCons(scip, andvars[v], cons, FALSE, TRUE) );
               }
            }
            else
            {
               for( v = nandvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPunlockVarCons(scip, andvars[v], cons, TRUE, FALSE) );
               }
            }
         }
      }
   }

   /* check whether the right hand side is decreased, if and only if that's the case we maybe can propagate, tighten and add more cliques */
   if( SCIPisGT(scip, oldrhs, rhs) )
   {
      consdata->propagated = FALSE;
   }

   /* set new right hand side and update constraint data */
   SCIP_CALL( chgRhsLinearCons(scip, consdata->lincons, consdata->linconstype, rhs) );
   consdata->rhs = rhs;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &andnegs);
   SCIPfreeBufferArray(scip, &andcoefs);
   SCIPfreeBufferArray(scip, &andress);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** create and-constraints and get all and-resultants */
static
SCIP_RETCODE createAndAddAnds(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*const   conshdlr,           /**< pseudoboolean constraint handler */
   SCIP_VAR**const*const terms,              /**< array of term variables to get and-constraints for */
   SCIP_Real*const       termcoefs,          /**< array of coefficients for and-constraints */
   int const             nterms,             /**< number of terms to get and-constraints for */
   int const*const       ntermvars,          /**< array of number of variable in each term */
   SCIP_Bool const       initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool const       enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant
                                              *   constraints. */
   SCIP_Bool const       check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant
                                              *   constraints. */
   SCIP_Bool const       local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching
                                              *   constraints. */
   SCIP_Bool const       modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE
                                              *   if pricing adds coefficients to this constraint. */
   SCIP_Bool const       dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool const       stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent
                                              *   node data. */
   SCIP_CONS**const      andconss,           /**< array to store all created and-constraints for given terms */
   SCIP_Real*const       andvals,            /**< array to store all coefficients of and-constraints */
   SCIP_Bool*const       andnegs,            /**< array to store negation status of and-constraints */
   int*const             nandconss           /**< number of created and constraints */
   )
{
   int t;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nterms == 0 || (terms != NULL && ntermvars != NULL));
   assert(andconss != NULL);
   assert(andvals != NULL);
   assert(nandconss != NULL);

   (*nandconss) = 0;

   if( nterms == 0 )
      return SCIP_OKAY;

   /* loop over all terms and create/get all and constraints */
   for( t = 0; t < nterms; ++t )
   {
      if( !SCIPisZero(scip, termcoefs[t]) && ntermvars[t] > 0 )
      {
         SCIP_CALL( createAndAddAndCons(scip, conshdlr, terms[t], ntermvars[t],
               initial, enforce, check, local, modifiable, dynamic, stickingatnode,
               &(andconss[*nandconss])) );
         assert(andconss[*nandconss] != NULL);
         andvals[*nandconss] = termcoefs[t];
         andnegs[*nandconss] = FALSE;
         ++(*nandconss);
      }
   }

   return SCIP_OKAY;
}

/** created linear constraint of pseudo boolean constraint */
static
SCIP_RETCODE createAndAddLinearCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*const   conshdlr,           /**< pseudoboolean constraint handler */
   SCIP_VAR**const       linvars,            /**< linear variables */
   int const             nlinvars,           /**< number of linear variables */
   SCIP_Real*const       linvals,            /**< linear coefficients */
   SCIP_VAR**const       andress,            /**< and-resultant variables */
   int const             nandress,           /**< number of and-resultant variables */
   SCIP_Real const*const andvals,            /**< and-resultant coefficients */
   SCIP_Bool*const       andnegs,            /**< and-resultant negation status */
   SCIP_Real*const       lhs,                /**< pointer to left hand side of linear constraint */
   SCIP_Real*const       rhs,                /**< pointer to right hand side of linear constraint */
   SCIP_Bool const       initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool const       separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool const       enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant
                                              *   constraints. */
   SCIP_Bool const       check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant
                                              *   constraints. */
   SCIP_Bool const       propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool const       local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching
                                              *   constraints. */
   SCIP_Bool const       modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE
                                              *   if pricing adds coefficients to this constraint. */
   SCIP_Bool const       dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool const       removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user
                                              *   cuts'. */
   SCIP_Bool const       stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent
                                              *   node data. */
   SCIP_CONS**const      lincons,            /**< pointer to store created linear constraint */
   SCIP_LINEARCONSTYPE*const linconstype     /**< pointer to store the type of the linear constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* upgrconshdlr;
   SCIP_CONS* cons;
   char name[SCIP_MAXSTRLEN];
   int v;
   SCIP_Bool created;
   SCIP_Bool integral;
   int nzero;
   int ncoeffspone;
   int ncoeffsnone;
   int ncoeffspint;
   int ncoeffsnint;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nlinvars == 0 || (linvars != NULL && linvals != NULL));
   assert(nandress == 0 || (andress != NULL && andvals != NULL));
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(lincons != NULL);
   assert(linconstype != NULL);
   assert(nlinvars > 0 || nandress > 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   (*linconstype) = SCIP_LINEARCONSTYPE_INVALIDCONS;
   (*lincons) = NULL;
   cons = NULL;

   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "pseudoboolean_linear%d", conshdlrdata->nlinconss);
   ++(conshdlrdata->nlinconss);

   created = FALSE;

   if( !modifiable )
   {
      SCIP_Real val;
      int nvars;

      /* calculate some statistics for upgrading on linear constraint */
      nzero = 0;
      ncoeffspone = 0;
      ncoeffsnone = 0;
      ncoeffspint = 0;
      ncoeffsnint = 0;
      integral = TRUE;
      nvars = nlinvars + nandress;

      /* calculate information over linear part */
      for( v = nlinvars - 1; v >= 0; --v )
      {
         val = linvals[v];

         if( SCIPisZero(scip, val) )
         {
            ++nzero;
            continue;
         }
         if( SCIPisEQ(scip, val, 1.0) )
            ++ncoeffspone;
         else if( SCIPisEQ(scip, val, -1.0) )
            ++ncoeffsnone;
         else if( SCIPisIntegral(scip, val) )
         {
            if( SCIPisPositive(scip, val) )
               ++ncoeffspint;
            else
               ++ncoeffsnint;
         }
         else
         {
            integral = FALSE;
            break;
         }
      }

      if( integral )
      {
         /* calculate information over and-resultants */
         for( v = nandress - 1; v >= 0; --v )
         {
            val = andvals[v];

            if( SCIPisZero(scip, val) )
            {
               ++nzero;
               continue;
            }
            if( SCIPisEQ(scip, val, 1.0) )
               ++ncoeffspone;
            else if( SCIPisEQ(scip, val, -1.0) )
               ++ncoeffsnone;
            else if( SCIPisIntegral(scip, val) )
            {
               if( SCIPisPositive(scip, val) )
                  ++ncoeffspint;
               else
                  ++ncoeffsnint;
            }
            else
            {
               integral = FALSE;
               break;
            }
         }
      }

      SCIPdebugMsg(scip, "While creating the linear constraint of the pseudoboolean constraint we found %d zero coefficients that were removed\n", nzero);

      /* try to upgrade to a special linear constraint */
      if( integral )
      {
         upgrconshdlr = SCIPfindConshdlr(scip, "logicor");

         /* check, if linear constraint can be upgraded to logic or constraint
          * - logic or constraints consist only of binary variables with a
          *   coefficient of +1.0 or -1.0 (variables with -1.0 coefficients can be negated):
          *        lhs     <= x1 + ... + xp - y1 - ... - yn <= rhs
          * - negating all variables y = (1-Y) with negative coefficients gives:
          *        lhs + n <= x1 + ... + xp + Y1 + ... + Yn <= rhs + n
          * - negating all variables x = (1-X) with positive coefficients and multiplying with -1 gives:
          *        p - rhs <= X1 + ... + Xp + y1 + ... + yn <= p - lhs
          * - logic or constraints have left hand side of +1.0, and right hand side of +infinity: x(S) >= 1.0
          *    -> without negations:  (lhs == 1 - n  and  rhs == +inf)  or  (lhs == -inf  and  rhs = p - 1)
          */
         if( upgrconshdlr != NULL && nvars > 2 && ncoeffspone + ncoeffsnone == nvars
            && ((SCIPisEQ(scip, *lhs, 1.0 - ncoeffsnone) && SCIPisInfinity(scip, *rhs))
               || (SCIPisInfinity(scip, -*lhs) && SCIPisEQ(scip, *rhs, ncoeffspone - 1.0))) )
         {
            SCIP_VAR** transvars;
            int mult;

            SCIPdebugMsg(scip, "linear constraint will be logic-or constraint\n");

            /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
            mult = SCIPisInfinity(scip, *rhs) ? +1 : -1;

            /* get temporary memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );

            /* negate positive or negative variables */
            for( v = 0; v < nlinvars; ++v )
            {
               if( mult * linvals[v] > 0.0 )
                  transvars[v] = linvars[v];
               else
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, linvars[v], &transvars[v]) );
               }
               assert(transvars[v] != NULL);
            }

            /* negate positive or negative variables */
            for( v = 0; v < nandress; ++v )
            {
               if( mult * andvals[v] > 0.0 )
                  transvars[nlinvars + v] = andress[v];
               else
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, andress[v], &transvars[nlinvars + v]) );
                  andnegs[v] = TRUE;
               }
               assert(transvars[nlinvars + v] != NULL);
            }

            assert(!modifiable);
            /* create the constraint */
            SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, nvars, transvars,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

            created = TRUE;
            (*linconstype) = SCIP_LINEARCONSTYPE_LOGICOR;

            /* free temporary memory */
            SCIPfreeBufferArray(scip, &transvars);

	    *lhs = 1.0;
	    *rhs = SCIPinfinity(scip);
         }

         upgrconshdlr = SCIPfindConshdlr(scip, "setppc");

         /* check, if linear constraint can be upgraded to set partitioning, packing, or covering constraint
          * - all set partitioning / packing / covering constraints consist only of binary variables with a
          *   coefficient of +1.0 or -1.0 (variables with -1.0 coefficients can be negated):
          *        lhs     <= x1 + ... + xp - y1 - ... - yn <= rhs
          * - negating all variables y = (1-Y) with negative coefficients gives:
          *        lhs + n <= x1 + ... + xp + Y1 + ... + Yn <= rhs + n
          * - negating all variables x = (1-X) with positive coefficients and multiplying with -1 gives:
          *        p - rhs <= X1 + ... + Xp + y1 + ... + yn <= p - lhs
          * - a set partitioning constraint has left hand side of +1.0, and right hand side of +1.0 : x(S) == 1.0
          *    -> without negations:  lhs == rhs == 1 - n  or  lhs == rhs == p - 1
          * - a set packing constraint has left hand side of -infinity, and right hand side of +1.0 : x(S) <= 1.0
          *    -> without negations:  (lhs == -inf  and  rhs == 1 - n)  or  (lhs == p - 1  and  rhs = +inf)
          * - a set covering constraint has left hand side of +1.0, and right hand side of +infinity: x(S) >= 1.0
          *    -> without negations:  (lhs == 1 - n  and  rhs == +inf)  or  (lhs == -inf  and  rhs = p - 1)
          */
         if( upgrconshdlr != NULL && !created && ncoeffspone + ncoeffsnone == nvars )
         {
            SCIP_VAR** transvars;
            int mult;

            if( SCIPisEQ(scip, *lhs, *rhs) && (SCIPisEQ(scip, *lhs, 1.0 - ncoeffsnone) || SCIPisEQ(scip, *lhs, ncoeffspone - 1.0)) )
            {
               SCIPdebugMsg(scip, "linear pseudoboolean constraint will be a set partitioning constraint\n");

               /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
               mult = SCIPisEQ(scip, *lhs, 1.0 - ncoeffsnone) ? +1 : -1;

               /* get temporary memory */
               SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );

               /* negate positive or negative variables for linear variables */
               for( v = 0; v < nlinvars; ++v )
               {
                  if( mult * linvals[v] > 0.0 )
                     transvars[v] = linvars[v];
                  else
                  {
                     SCIP_CALL( SCIPgetNegatedVar(scip, linvars[v], &transvars[v]) );
                  }
                  assert(transvars[v] != NULL);
               }

               /* negate positive or negative variables for and-resultants */
               for( v = 0; v < nandress; ++v )
               {
                  if( mult * andvals[v] > 0.0 )
                     transvars[nlinvars + v] = andress[v];
                  else
                  {
                     SCIP_CALL( SCIPgetNegatedVar(scip, andress[v], &transvars[nlinvars + v]) );
                     andnegs[v] = TRUE;
                  }
                  assert(transvars[nlinvars + v] != NULL);
               }

               /* create the constraint */
               assert(!modifiable);
               SCIP_CALL( SCIPcreateConsSetpart(scip, &cons, name, nvars, transvars,
                     initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

               created = TRUE;
               (*linconstype) = SCIP_LINEARCONSTYPE_SETPPC;

               /* release temporary memory */
               SCIPfreeBufferArray(scip, &transvars);

	       *lhs = 1.0;
	       *rhs = 1.0;
            }
            else if( (SCIPisInfinity(scip, -*lhs) && SCIPisEQ(scip, *rhs, 1.0 - ncoeffsnone))
               || (SCIPisEQ(scip, *lhs, ncoeffspone - 1.0) && SCIPisInfinity(scip, *rhs)) )
            {
               SCIPdebugMsg(scip, "linear pseudoboolean constraint will be a set packing constraint\n");

               /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
               mult = SCIPisInfinity(scip, -*lhs) ? +1 : -1;

               /* get temporary memory */
               SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );

               /* negate positive or negative variables for linear variables */
               for( v = 0; v < nlinvars; ++v )
               {
                  if( mult * linvals[v] > 0.0 )
                     transvars[v] = linvars[v];
                  else
                  {
                     SCIP_CALL( SCIPgetNegatedVar(scip, linvars[v], &transvars[v]) );
                  }
                  assert(transvars[v] != NULL);
               }

               /* negate positive or negative variables for and-resultants*/
               for( v = 0; v < nandress; ++v )
               {
                  if( mult * andvals[v] > 0.0 )
                     transvars[nlinvars + v] = andress[v];
                  else
                  {
                     SCIP_CALL( SCIPgetNegatedVar(scip, andress[v], &transvars[nlinvars + v]) );
                     andnegs[v] = TRUE;
                  }
                  assert(transvars[nlinvars + v] != NULL);
               }

               /* create the constraint */
               assert(!modifiable);
               SCIP_CALL( SCIPcreateConsSetpack(scip, &cons, name, nvars, transvars,
                     initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

               created = TRUE;
               (*linconstype) = SCIP_LINEARCONSTYPE_SETPPC;

               /* release temporary memory */
               SCIPfreeBufferArray(scip, &transvars);

	       *lhs = -SCIPinfinity(scip);
	       *rhs = 1.0;
            }
            else if( (SCIPisEQ(scip, *lhs, 1.0 - ncoeffsnone) && SCIPisInfinity(scip, *rhs))
               || (SCIPisInfinity(scip, -*lhs) && SCIPisEQ(scip, *rhs, ncoeffspone - 1.0)) )
            {
               if( nvars != 1 )
               {
                  if( nvars == 2 )
                  {
                     SCIPwarningMessage(scip, "Does not expect this, because this constraint should be a set packing constraint.\n");
                  }
                  else
                  {
                     SCIPwarningMessage(scip, "Does not expect this, because this constraint should be a logicor constraint.\n");
                  }
               }
               SCIPdebugMsg(scip, "linear pseudoboolean constraint will be a set covering constraint\n");

               /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
               mult = SCIPisInfinity(scip, *rhs) ? +1 : -1;

               /* get temporary memory */
               SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );

               /* negate positive or negative variables for linear variables */
               for( v = 0; v < nlinvars; ++v )
               {
                  if( mult * linvals[v] > 0.0 )
                     transvars[v] = linvars[v];
                  else
                  {
                     SCIP_CALL( SCIPgetNegatedVar(scip, linvars[v], &transvars[v]) );
                  }
                  assert(transvars[v] != NULL);
               }

               /* negate positive or negative variables for and-resultants*/
               for( v = 0; v < nandress; ++v )
               {
                  if( mult * andvals[v] > 0.0 )
                     transvars[nlinvars + v] = andress[v];
                  else
                  {
                     SCIP_CALL( SCIPgetNegatedVar(scip, andress[v], &transvars[nlinvars + v]) );
                     andnegs[v] = TRUE;
                  }
                  assert(transvars[nlinvars + v] != NULL);
               }

               /* create the constraint */
               assert(!modifiable);
               SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, name, nvars, transvars,
                     initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

               created = TRUE;
               (*linconstype) = SCIP_LINEARCONSTYPE_SETPPC;

               /* release temporary memory */
               SCIPfreeBufferArray(scip, &transvars);

	       *lhs = 1.0;
	       *rhs = SCIPinfinity(scip);
            }
         }

         upgrconshdlr = SCIPfindConshdlr(scip, "knapsack");

         /* check, if linear constraint can be upgraded to a knapsack constraint
          * - all variables must be binary
          * - all coefficients must be integral
          * - exactly one of the sides must be infinite
          */
         if( upgrconshdlr != NULL && !created && (ncoeffspone + ncoeffsnone + ncoeffspint + ncoeffsnint == nvars) && (SCIPisInfinity(scip, -*lhs) != SCIPisInfinity(scip, *rhs)) )
         {
            SCIP_VAR** transvars;
            SCIP_Longint* weights;
            SCIP_Longint capacity;
            SCIP_Longint weight;
            int mult;

            SCIPdebugMsg(scip, "linear pseudoboolean constraint will be a knapsack constraint\n");

            /* get temporary memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

            /* if the right hand side is non-infinite, we have to negate all variables with negative coefficient;
             * otherwise, we have to negate all variables with positive coefficient and multiply the row with -1
             */
            if( SCIPisInfinity(scip, *rhs) )
            {
               mult = -1;
               capacity = (SCIP_Longint)SCIPfeasFloor(scip, -*lhs);
            }
            else
            {
               mult = +1;
               capacity = (SCIP_Longint)SCIPfeasFloor(scip, *rhs);
            }

            /* negate positive or negative variables for linear variables */
            for( v = 0; v < nlinvars; ++v )
            {
               assert(SCIPisFeasIntegral(scip, linvals[v]));
               weight = mult * (SCIP_Longint)SCIPfeasFloor(scip, linvals[v]);
               if( weight > 0 )
               {
                  transvars[v] = linvars[v];
                  weights[v] = weight;
               }
               else
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, linvars[v], &transvars[v]) );
                  weights[v] = -weight;
                  capacity -= weight;
               }
               assert(transvars[v] != NULL);
            }
            /* negate positive or negative variables for and-resultants */
            for( v = 0; v < nandress; ++v )
            {
               assert(SCIPisFeasIntegral(scip, andvals[v]));
               weight = mult * (SCIP_Longint)SCIPfeasFloor(scip, andvals[v]);
               if( weight > 0 )
               {
                  transvars[nlinvars + v] = andress[v];
                  weights[nlinvars + v] = weight;
               }
               else
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, andress[v], &transvars[nlinvars + v]) );
                  andnegs[v] = TRUE;
                  weights[nlinvars + v] = -weight;
                  capacity -= weight;
               }
               assert(transvars[nlinvars + v] != NULL);
            }

            /* create the constraint */
            SCIP_CALL( SCIPcreateConsKnapsack(scip, &cons, name, nvars, transvars, weights, capacity,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

            created = TRUE;
            (*linconstype) = SCIP_LINEARCONSTYPE_KNAPSACK;

            /* free temporary memory */
            SCIPfreeBufferArray(scip, &weights);
            SCIPfreeBufferArray(scip, &transvars);

	    *lhs = -SCIPinfinity(scip);
	    *rhs = capacity;
         }
#ifdef WITHEQKNAPSACK

         upgrconshdlr = SCIPfindConshdlr(scip, "eqknapsack");

         /* check, if linear constraint can be upgraded to a knapsack constraint
          * - all variables must be binary
          * - all coefficients must be integral
          * - both sides must be infinite
          */
         if( upgrconshdlr != NULL && !created && (ncoeffspone + ncoeffsnone + ncoeffspint + ncoeffsnint == nvars) && SCIPisEQ(scip, *lhs, *rhs) )
         {
            SCIP_VAR** transvars;
            SCIP_Longint* weights;
            SCIP_Longint capacity;
            SCIP_Longint weight;
            int mult;

            assert(!SCIPisInfinity(scip, *rhs));

            SCIPdebugMsg(scip, "linear pseudoboolean constraint will be a equality-knapsack constraint\n");

            /* get temporary memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

            if( SCIPisPositive(scip, *rhs) )
            {
               mult = +1;
               capacity = (SCIP_Longint)SCIPfeasFloor(scip, *rhs);
            }
            else
            {
               mult = -1;
               capacity = (SCIP_Longint)SCIPfeasFloor(scip, -*rhs);
            }

            /* negate positive or negative variables for linear variables */
            for( v = 0; v < nlinvars; ++v )
            {
               assert(SCIPisFeasIntegral(scip, linvals[v]));
               weight = mult * (SCIP_Longint)SCIPfeasFloor(scip, linvals[v]);
               if( weight > 0 )
               {
                  transvars[v] = linvars[v];
                  weights[v] = weight;
               }
               else
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, linvars[v], &transvars[v]) );
                  weights[v] = -weight;
                  capacity -= weight;
               }
               assert(transvars[v] != NULL);
            }
            /* negate positive or negative variables for and-resultants */
            for( v = 0; v < nandress; ++v )
            {
               assert(SCIPisFeasIntegral(scip, andvals[v]));
               weight = mult * (SCIP_Longint)SCIPfeasFloor(scip, andvals[v]);
               if( weight > 0 )
               {
                  transvars[nlinvars + v] = andress[v];
                  weights[nlinvars + v] = weight;
               }
               else
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, andress[v], &transvars[nlinvars + v]) );
                  andnegs[v] = TRUE;
                  weights[nlinvars + v] = -weight;
                  capacity -= weight;
               }
               assert(transvars[nlinvars + v] != NULL);
            }

            /* create the constraint */
            SCIP_CALL( SCIPcreateConsEqKnapsack(scip, &cons, name, nvars, transvars, weights, capacity,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

            created = TRUE;
            (*linconstype) = SCIP_LINEARCONSTYPE_EQKNAPSACK;

            /* free temporary memory */
            SCIPfreeBufferArray(scip, &weights);
            SCIPfreeBufferArray(scip, &transvars);

	    *lhs = capacity;
	    *rhs = capacity;
         }
#endif
      }
   }

   upgrconshdlr = SCIPfindConshdlr(scip, "linear");
   assert(created || upgrconshdlr != NULL);

   if( !created )
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nlinvars, linvars, linvals, *lhs, *rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

      (*linconstype) = SCIP_LINEARCONSTYPE_LINEAR;

      /* add all and-resultants */
      for( v = 0; v < nandress; ++v )
      {
         assert(andress[v] != NULL);

         /* add auxiliary variables to linear constraint */
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, andress[v], andvals[v]) );
      }
   }

   assert(cons != NULL && *linconstype > SCIP_LINEARCONSTYPE_INVALIDCONS);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugPrintCons(scip, cons, NULL);

   *lincons = cons;
   SCIP_CALL( SCIPcaptureCons(scip, *lincons) );

   /* mark linear constraint not to be upgraded - otherwise we loose control over it */
   SCIPconsAddUpgradeLocks(cons, 1);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** checks one original pseudoboolean constraint for feasibility of given solution */
static
SCIP_RETCODE checkOrigPbCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudo boolean constraint */
   SCIP_SOL*const        sol,                /**< solution to be checked, or NULL for current solution */
   SCIP_Bool*const       violated,           /**< pointer to store whether the constraint is violated */
   SCIP_Bool const       printreason         /**< should violation of constraint be printed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int nvars;
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlinvars;
   int v;

   SCIP_VAR** andress;
   SCIP_Real* andcoefs;
   int nandress;

   SCIP_CONS* andcons;
   SCIP_Real andvalue;
   SCIP_Real activity;
   int c;

   SCIP_Real lhsviol;
   SCIP_Real rhsviol;
   SCIP_Real absviol;
   SCIP_Real relviol;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsOriginal(cons));
   assert(violated != NULL);

   *violated = FALSE;

   SCIPdebugMsg(scip, "checking original pseudo boolean constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lincons != NULL);
   assert(consdata->linconstype > SCIP_LINEARCONSTYPE_INVALIDCONS);
   assert(SCIPconsIsOriginal(consdata->lincons));

   /* gets number of variables in linear constraint */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andress, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andcoefs, nvars) );

   /* get sides of linear constraint */
   SCIP_CALL( getLinearConsSides(scip, consdata->lincons, consdata->linconstype, &lhs, &rhs) );
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(SCIPisLE(scip, lhs, rhs));

   /* get variables and coefficient of linear constraint */
   SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );
   assert(nvars == 0 || (coefs != NULL));

   /* number of variables should be consistent, number of 'real' linear variables plus number of and-constraints should
    * have to be equal to the number of variables in the linear constraint
    */
   assert(consdata->nlinvars + consdata->nconsanddatas == nvars);

   nlinvars = 0;

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->hashmap != NULL);

   nandress = 0;

   activity = 0.0;

   /* split variables into original and artificial variables and compute activity on normal linear variables (without
    * terms)
    */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* hashmapvar;
      SCIP_Bool negated;

      assert(vars[v] != NULL);

      /* negated variables can also exist in the original problem, so we need to check */
      if( !SCIPhashmapExists(conshdlrdata->hashmap, (void*)(vars[v])) && SCIPvarIsNegated(vars[v]) )
      {
         hashmapvar = SCIPvarGetNegationVar(vars[v]);
         negated = TRUE;
      }
      else
      {
         hashmapvar = vars[v];
         negated = FALSE;
      }
      assert(hashmapvar != NULL);

      if( !SCIPhashmapExists(conshdlrdata->hashmap, (void*)(hashmapvar)) )
      {
         assert(!SCIPhashmapExists(conshdlrdata->hashmap, (void*)(vars[v])));

         activity += coefs[v] * SCIPgetSolVal(scip, sol, vars[v]);

         linvars[nlinvars] = vars[v];
         lincoefs[nlinvars] = coefs[v];
         ++nlinvars;
      }
      else
      {
         /* negate coefficient in case of an original negated variable */
         andress[nandress] = hashmapvar;
         if( negated )
         {
            if( !SCIPisInfinity(scip, -lhs) )
               lhs -= coefs[v];
            if( !SCIPisInfinity(scip, rhs) )
               rhs -= coefs[v];
            andcoefs[nandress] = -coefs[v];
         }
         else
            andcoefs[nandress] = coefs[v];
         ++nandress;
      }
   }
   assert(nandress == consdata->nconsanddatas);

   SCIPsortPtrReal((void**)andress, andcoefs, SCIPvarComp, nandress);

   SCIPdebugMsg(scip, "nlinvars = %d, nandress = %d\n", nlinvars, nandress);
   SCIPdebugMsg(scip, "linear activity = %g\n", activity);

   /* compute and add solution values on terms */
   for( c = consdata->nconsanddatas - 1; c >= 0; --c )
   {
      SCIP_VAR** andvars;
      int nandvars;
#ifndef NDEBUG
      SCIP_VAR* res;
#endif
      andcons = consdata->consanddatas[c]->origcons;

      /* if after during or before presolving a solution will be transformed into original space and will be checked
       * there, but origcons was already removed and only the pointer to the transformed and-constraint is existing
       */
      if( andcons == NULL )
      {
         andcons = consdata->consanddatas[c]->cons;
      }
      assert(andcons != NULL);

      andvars = SCIPgetVarsAnd(scip, andcons);
      nandvars = SCIPgetNVarsAnd(scip, andcons);

#ifndef NDEBUG
      res = SCIPgetResultantAnd(scip, andcons);
      assert(nandvars == 0 || (andvars != NULL && res != NULL));
      assert(res == andress[c]);
#endif

      andvalue = 1;
      /* check if the and-constraint is violated */
      for( v = nandvars - 1; v >= 0; --v )
      {
         andvalue *= SCIPgetSolVal(scip, sol, andvars[v]);
         if( SCIPisFeasZero(scip, andvalue) )
            break;
      }
      activity += andvalue * andcoefs[c];
   }
   SCIPdebugMsg(scip, "lhs = %g, overall activity = %g, rhs = %g\n", lhs, activity, rhs);

   /* calculate absolute and relative violation */
   lhsviol = lhs - activity;
   rhsviol = activity - rhs;

   if(lhsviol > rhsviol)
   {
      absviol = lhsviol;
      relviol = SCIPrelDiff(lhs, activity);
   }
   else
   {
      absviol = rhsviol;
      relviol = SCIPrelDiff(activity, rhs);
   }

   /* update absolute and relative violation of the solution */
   if( sol != NULL )
      SCIPupdateSolConsViolation(scip, sol, absviol, relviol);

   /* check left hand side for violation */
   if( SCIPisFeasLT(scip, activity, lhs) )
   {
      if( printreason )
      {
         SCIP_CALL( SCIPprintCons(scip, cons, NULL ) );
	 SCIPinfoMessage(scip, NULL, ";\n");
         SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", lhs - activity);

         /* print linear constraint in SCIP_DEBUG mode too */
         SCIPdebugPrintCons(scip, SCIPconsGetData(cons)->lincons, NULL);
      }

      *violated = TRUE;
   }

   /* check right hand side for violation */
   if( SCIPisFeasGT(scip, activity, rhs) )
   {
      if( printreason )
      {
         SCIP_CALL( SCIPprintCons(scip, cons, NULL ) );
	 SCIPinfoMessage(scip, NULL, ";\n");
         SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", activity - rhs);
      }

      *violated = TRUE;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &andcoefs);
   SCIPfreeBufferArray(scip, &andress);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** checks all and-constraints inside the pseudoboolean constraint handler for feasibility of given solution or current
 *  solution
 */
static
SCIP_RETCODE checkAndConss(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*const   conshdlr,           /**< pseudo boolean constraint handler */
   SCIP_SOL*const        sol,                /**< solution to be checked, or NULL for current solution */
   SCIP_Bool*const       violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* andcons;
   SCIP_VAR** vars;
   SCIP_VAR* res;
   int nvars;
   SCIP_Real andvalue;
   int c;
   int v;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(violated != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *violated = FALSE;

   for( c = conshdlrdata->nallconsanddatas - 1; c >= 0; --c )
   {
      if( !conshdlrdata->allconsanddatas[c]->istransformed )
         continue;

      andcons = conshdlrdata->allconsanddatas[c]->cons;

      /* need to check even locally deleted constraints */
      if( andcons == NULL ) /*|| !SCIPconsIsActive(andcons) )*/
         continue;

      vars = SCIPgetVarsAnd(scip, andcons);
      nvars = SCIPgetNVarsAnd(scip, andcons);
      res = SCIPgetResultantAnd(scip, andcons);
      assert(nvars == 0 || (vars != NULL && res != NULL));

      andvalue = 1;
      /* check if the and-constraint is violated */
      for( v = nvars - 1; v >= 0; --v )
      {
         andvalue *= SCIPgetSolVal(scip, sol, vars[v]);
         if( SCIPisFeasZero(scip, andvalue) )
            break;
      }

      /* check for violation and update aging */
      if( !SCIPisFeasEQ(scip, andvalue, SCIPgetSolVal(scip, sol, res)) )
      {
         /* only reset constraint age if we are in enforcement */
         if( sol == NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, andcons) );
         }

         *violated = TRUE;
         break;
      }
      else if( sol == NULL )
      {
         SCIP_CALL( SCIPincConsAge(scip, andcons) );
      }
   }

   return SCIP_OKAY;
}

/** creates by copying and captures a linear constraint */
static
SCIP_RETCODE copyConsPseudoboolean(
   SCIP*const            targetscip,         /**< target SCIP data structure */
   SCIP_CONS**           targetcons,         /**< pointer to store the created target constraint */
   SCIP*const            sourcescip,         /**< source SCIP data structure */
   SCIP_CONS*const       sourcecons,         /**< source constraint which will be copied */
   const char*           name,               /**< name of constraint */
   SCIP_HASHMAP*const    varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_HASHMAP*const    consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_Bool const       initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool const       separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool const       enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool const       check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool const       propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool const       local,              /**< is constraint only valid locally? */
   SCIP_Bool const       modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool const       dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool const       removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool const       stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool const       global,             /**< create a global or a local copy? */
   SCIP_Bool*const       valid               /**< pointer to store if the copying was valid */
   )
{
   SCIP_CONSDATA* sourceconsdata;
   SCIP_CONS* sourcelincons;

   assert(targetscip != NULL);
   assert(targetcons != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0);
   assert(valid != NULL);

   *valid = TRUE;

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get linear constraint */
   sourcelincons = sourceconsdata->lincons;
   assert(sourcelincons != NULL);

   /* get copied version of linear constraint */
   if( !SCIPconsIsDeleted(sourcelincons) )
   {
      SCIP_CONSHDLR* conshdlrlinear;
      SCIP_CONS* targetlincons;
      SCIP_CONS** targetandconss;
      SCIP_Real* targetandcoefs;
      int ntargetandconss;
      SCIP_LINEARCONSTYPE targetlinconstype;

      targetlinconstype = sourceconsdata->linconstype;

      switch( targetlinconstype )
      {
      case SCIP_LINEARCONSTYPE_LINEAR:
         conshdlrlinear = SCIPfindConshdlr(sourcescip, "linear");
         assert(conshdlrlinear != NULL);
         break;
      case SCIP_LINEARCONSTYPE_LOGICOR:
         conshdlrlinear = SCIPfindConshdlr(sourcescip, "logicor");
         assert(conshdlrlinear != NULL);
         break;
      case SCIP_LINEARCONSTYPE_KNAPSACK:
         conshdlrlinear = SCIPfindConshdlr(sourcescip, "knapsack");
         assert(conshdlrlinear != NULL);
         break;
      case SCIP_LINEARCONSTYPE_SETPPC:
         conshdlrlinear = SCIPfindConshdlr(sourcescip, "setppc");
         assert(conshdlrlinear != NULL);
         break;
#ifdef WITHEQKNAPSACK
      case SCIP_LINEARCONSTYPE_EQKNAPSACK:
         conshdlrlinear = SCIPfindConshdlr(sourcescip, "eqknapsack");
         assert(conshdlrlinear != NULL);
         break;
#endif
      case SCIP_LINEARCONSTYPE_INVALIDCONS:
      default:
         SCIPerrorMessage("unknown linear constraint type\n");
         return SCIP_INVALIDDATA;
      }

      if( conshdlrlinear == NULL ) /*lint !e774*/
      {
         SCIPerrorMessage("linear constraint handler not found\n");
         return SCIP_INVALIDDATA;
      }

      targetlincons = NULL;

      /* copy linear constraint */
      SCIP_CALL( SCIPgetConsCopy(sourcescip, targetscip, sourcelincons, &targetlincons, conshdlrlinear, varmap, consmap, SCIPconsGetName(sourcelincons),
            SCIPconsIsInitial(sourcelincons), SCIPconsIsSeparated(sourcelincons), SCIPconsIsEnforced(sourcelincons), SCIPconsIsChecked(sourcelincons),
            SCIPconsIsPropagated(sourcelincons), SCIPconsIsLocal(sourcelincons), SCIPconsIsModifiable(sourcelincons), SCIPconsIsDynamic(sourcelincons),
            SCIPconsIsRemovable(sourcelincons), SCIPconsIsStickingAtNode(sourcelincons), global, valid) );

      if( *valid )
      {
         assert(targetlincons != NULL);
         assert(SCIPconsGetHdlr(targetlincons) != NULL);
         /* @note  due to copying special linear constraints, now leads only to simple linear constraints, we check that
          *        our target constraint handler is the same as our source constraint handler of the linear constraint,
          *        if not copying was not valid
          */
         if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(targetlincons)), "linear") == 0 )
            targetlinconstype = SCIP_LINEARCONSTYPE_LINEAR;
      }

      targetandconss = NULL;
      targetandcoefs = NULL;
      ntargetandconss = 0;

      if( *valid )
      {
         SCIP_CONSHDLR* conshdlrand;
         int c;
         int nsourceandconss;
         SCIP_HASHTABLE* linconsvarsmap;
         SCIP_VAR** targetlinvars;
         SCIP_Real* targetlincoefs;
         int ntargetlinvars;

         conshdlrand = SCIPfindConshdlr(sourcescip, "and");
         assert(conshdlrand != NULL);

         nsourceandconss = sourceconsdata->nconsanddatas;

         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetandconss, nsourceandconss) );
         SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetandcoefs, nsourceandconss) );

         /* get the number of vars in the copied linear constraint and allocate buffers
          * for the variables and the coefficients
          */
         SCIP_CALL( getLinearConsNVars(targetscip, targetlincons, targetlinconstype, &ntargetlinvars) );
         SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetlinvars, ntargetlinvars) );
         SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetlincoefs, ntargetlinvars) );

         /* retrieve the variables of the copied linear constraint */
         SCIP_CALL( getLinearConsVarsData(targetscip, targetlincons, targetlinconstype,
                                          targetlinvars, targetlincoefs, &ntargetlinvars) );

         /* now create a hashtable and insert the variables into it, so that it
          * can be checked in constant time if a variable was removed due to
          * compressed copying when looping over the and resultants
          */
         SCIP_CALL( SCIPhashtableCreate(&linconsvarsmap, SCIPblkmem(targetscip), ntargetlinvars, SCIPvarGetHashkey,
                                        SCIPvarIsHashkeyEq, SCIPvarGetHashkeyVal, NULL) );

         for( c = 0 ; c < ntargetlinvars; ++c )
         {
            SCIP_CALL( SCIPhashtableInsert(linconsvarsmap, targetlinvars[c]) );
         }

         /* free the buffer arrays that were only required for building the hastable */
         SCIPfreeBufferArray(sourcescip, &targetlincoefs);
         SCIPfreeBufferArray(sourcescip, &targetlinvars);

         for( c = 0 ; c < nsourceandconss; ++c )
         {
            CONSANDDATA* consanddata;
            SCIP_CONS* oldcons;
            SCIP_VAR* targetandresultant;
            SCIP_Bool validand;

            consanddata = sourceconsdata->consanddatas[c];
            assert(consanddata != NULL);

            oldcons = consanddata->cons;
            assert(oldcons != NULL);

            targetandresultant = (SCIP_VAR*) SCIPhashmapGetImage(varmap, SCIPgetResultantAnd(sourcescip, oldcons));
            assert(targetandresultant != NULL);

            /* if compressed copying is active, the resultant might not have been copied by the linear
             * constraint and we don't need to add it to the pseudo boolean constraint in this case
             */
            if( !SCIPhashtableExists(linconsvarsmap, targetandresultant) )
               continue;

            validand = TRUE;

            targetandconss[ntargetandconss] = NULL;

            /* copy and-constraints */
            SCIP_CALL( SCIPgetConsCopy(sourcescip, targetscip, oldcons, &targetandconss[ntargetandconss], conshdlrand, varmap, consmap, SCIPconsGetName(oldcons),
                  SCIPconsIsInitial(oldcons), SCIPconsIsSeparated(oldcons), SCIPconsIsEnforced(oldcons), SCIPconsIsChecked(oldcons),
                  SCIPconsIsPropagated(oldcons), SCIPconsIsLocal(oldcons), SCIPconsIsModifiable(oldcons), SCIPconsIsDynamic(oldcons),
                  SCIPconsIsRemovable(oldcons), SCIPconsIsStickingAtNode(oldcons), global, &validand) );

            *valid &= validand;

            if( validand )
            {
               targetandcoefs[ntargetandconss] = sourceconsdata->andcoefs[c];
               ++ntargetandconss;
            }
         }

         SCIPhashtableFree(&linconsvarsmap);
         assert(ntargetandconss <= ntargetlinvars);
      }

      /* no correct pseudoboolean constraint */
      if( ntargetandconss == 0 )
      {
         SCIPdebugMsg(sourcescip, "no and-constraints copied for pseudoboolean constraint <%s>\n", SCIPconsGetName(sourcecons));
         *valid = FALSE;
      }

      if( *valid )
      {
         SCIP_Real targetrhs;
         SCIP_Real targetlhs;

         SCIP_VAR* intvar;
         SCIP_VAR* indvar;
         const char* consname;

         /* third the indicator and artificial integer variable part */
         assert(sourceconsdata->issoftcons == (sourceconsdata->indvar != NULL));
         indvar = sourceconsdata->indvar;
         intvar = sourceconsdata->intvar;

         /* copy indicator variable */
         if( indvar != NULL )
         {
            assert(*valid);
            SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, indvar, &indvar, varmap, consmap, global, valid) );
            assert(!(*valid) || indvar != NULL);
         }
         /* copy artificial integer variable */
         if( intvar != NULL && *valid )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, intvar, &intvar, varmap, consmap, global, valid) );
            assert(!(*valid) || intvar != NULL);
         }

         if( name != NULL )
            consname = name;
         else
            consname = SCIPconsGetName(sourcecons);

         /* get new left and right hand sides of copied linear constraint since
          * they might have changed if compressed copying is used
          */
         SCIP_CALL( getLinearConsSides(targetscip, targetlincons, targetlinconstype, &targetlhs, &targetrhs) );

         /* create new pseudoboolean constraint */
         SCIP_CALL( SCIPcreateConsPseudobooleanWithConss(targetscip, targetcons, consname,
               targetlincons, targetlinconstype, targetandconss, targetandcoefs, ntargetandconss,
               indvar, sourceconsdata->weight, sourceconsdata->issoftcons, intvar, targetlhs, targetrhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      }
      else if( !SCIPisConsCompressionEnabled(sourcescip) )
      {
         SCIPverbMessage(sourcescip, SCIP_VERBLEVEL_MINIMAL, NULL, "could not copy constraint <%s>\n", SCIPconsGetName(sourcecons));
      }

      /* release copied linear constraint */
      if( targetlincons != NULL )
      {
         SCIP_CALL( SCIPreleaseCons(targetscip, &targetlincons) );
      }

      /* release copied and constraint */
      if( targetandconss != NULL )
      {
         int c;

         assert(ntargetandconss <= sourceconsdata->nconsanddatas);

         for( c = 0 ; c < ntargetandconss; ++c )
         {
            if( targetandconss[c] != NULL )
            {
               SCIP_CALL( SCIPreleaseCons(targetscip, &targetandconss[c]) );
            }
         }
      }

      /* free temporary memory */
      SCIPfreeBufferArrayNull(sourcescip, &targetandcoefs);
      SCIPfreeBufferArrayNull(sourcescip, &targetandconss);
   }
   else
      *valid = FALSE;

   return SCIP_OKAY;
}

/** compute all changes in consanddatas array */
static
SCIP_RETCODE computeConsAndDataChanges(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*const conshdlrdata      /**< pseudoboolean constraint handler data */
   )
{
   CONSANDDATA** allconsanddatas;
   CONSANDDATA* consanddata;
   int c;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   allconsanddatas = conshdlrdata->allconsanddatas;
   assert(allconsanddatas != NULL);
   assert(conshdlrdata->nallconsanddatas > 0);
   assert(conshdlrdata->nallconsanddatas <= conshdlrdata->sallconsanddatas);

   for( c = conshdlrdata->nallconsanddatas - 1; c >= 0; --c )
   {
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      int nvars;
      SCIP_VAR** newvars;
      int nnewvars;
      int v;

      consanddata = allconsanddatas[c];

      if( !consanddata->istransformed )
         continue;

      if( consanddata->nuses == 0 )
         continue;

      vars = consanddata->vars;
      nvars = consanddata->nvars;
      assert(nvars == 0 || vars != NULL);
      assert(consanddata->nnewvars == 0 && ((consanddata->snewvars > 0) == (consanddata->newvars != NULL)));

      if( nvars == 0 )
      {
#ifndef NDEBUG
         /* if an old consanddata-object has no variables left there should be no new variables */
         if( consanddata->cons != NULL )
            assert(SCIPgetNVarsAnd(scip, consanddata->cons) == 0);
#endif
         continue;
      }

      cons = consanddata->cons;
      assert(cons != NULL);

      if( SCIPconsIsDeleted(cons) )
         continue;

      /* sort and-variables */
      if( !SCIPisAndConsSorted(scip, consanddata->cons) )
      {
         SCIP_CALL( SCIPsortAndCons(scip, consanddata->cons) );
         assert(SCIPisAndConsSorted(scip, consanddata->cons));
      }

      /* get new and-variables */
      nnewvars = SCIPgetNVarsAnd(scip, consanddata->cons);
      newvars = SCIPgetVarsAnd(scip, consanddata->cons);

      /* stop if the constraint has no variables or there was an error (coverity issue) */
      if( nnewvars <= 0 )
         continue;

#ifndef NDEBUG
      /* check that old variables are sorted */
      for( v = nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(vars[v]) >= SCIPvarGetIndex(vars[v - 1]));
      /* check that new variables are sorted */
      for( v = nnewvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(newvars[v]) >= SCIPvarGetIndex(newvars[v - 1]));
#endif

      /* check for changings, if and-constraint did not change we do not need to copy all variables */
      if( nvars == nnewvars )
      {
         SCIP_Bool changed;

         changed = FALSE;

         /* check each variable */
         for( v = nvars - 1; v >= 0; --v )
         {
            if( vars[v] != newvars[v] )
            {
               changed = TRUE;
               break;
            }
         }

         if( !changed )
            continue;
      }

      /* resize newvars array if necessary */
      if( nnewvars > consanddata->snewvars )
      {
         SCIP_CALL( SCIPensureBlockMemoryArray(scip, &(consanddata->newvars), &(consanddata->snewvars), nnewvars) );
      }

      /* copy all variables */
      BMScopyMemoryArray(consanddata->newvars, newvars, nnewvars);
      consanddata->nnewvars = nnewvars;

      /* capture all variables */
      for( v = consanddata->nnewvars - 1; v >= 0; --v )
      {
         /* in original problem the variables was already deleted */
         assert(consanddata->newvars[v] != NULL);
         SCIP_CALL( SCIPcaptureVar(scip, consanddata->newvars[v]) );
      }
   }

   return SCIP_OKAY;
}

/** remove old locks */
static
SCIP_RETCODE removeOldLocks(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   CONSANDDATA*const     consanddata,        /**< CONSANDDATA object for which we want to delete the locks and the
                                              *   capture of the corresponding and-constraint */
   SCIP_Real const       coef,               /**< coefficient which led to old locks */
   SCIP_Real const       lhs,                /**< left hand side which led to old locks */
   SCIP_Real const       rhs                 /**< right hand side which led to old locks */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(consanddata != NULL);
   assert(!SCIPisInfinity(scip, coef) && !SCIPisInfinity(scip, -coef));
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(SCIPisLE(scip, lhs, rhs));

   /* remove rounding locks */
   SCIP_CALL( unlockRoundingAndCons(scip, cons, consanddata, coef, lhs, rhs) );

   assert(consanddata->cons != NULL);

   return SCIP_OKAY;
}

/** add new locks */
static
SCIP_RETCODE addNewLocks(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   CONSANDDATA*const     consanddata,        /**< CONSANDDATA object for which we want to delete the locks and the
                                              *   capture of the corresponding and-constraint */
   SCIP_Real const       coef,               /**< coefficient which lead to new locks */
   SCIP_Real const       lhs,                /**< left hand side which lead to new locks */
   SCIP_Real const       rhs                 /**< right hand side which lead to new locks */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(consanddata != NULL);
   assert(!SCIPisInfinity(scip, coef) && !SCIPisInfinity(scip, -coef));
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(SCIPisLE(scip, lhs, rhs));

   /* add rounding locks due to old variables in consanddata object */
   SCIP_CALL( lockRoundingAndCons(scip, cons, consanddata, coef, lhs, rhs) );

   assert(consanddata->cons != NULL);

   return SCIP_OKAY;
}

/** update all locks inside this constraint and all captures on all and-constraints */
static
SCIP_RETCODE correctLocksAndCaptures(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONSHDLRDATA*const conshdlrdata,     /**< pseudoboolean constraint handler data */
   SCIP_Real const       newlhs,             /**< new left hand side of pseudoboolean constraint */
   SCIP_Real const       newrhs,             /**< new right hand side of pseudoboolean constraint */
   SCIP_VAR**const       andress,            /**< current and-resultants in pseudoboolean constraint */
   SCIP_Real*const       andcoefs,           /**< current and-resultants-coeffcients in pseudoboolean constraint */
   SCIP_Bool*const       andnegs,            /**< current negation status of and-resultants in pseudoboolean constraint */
   int const             nandress            /**< number of current and-resultants in pseudoboolean constraint */
   )
{
   CONSANDDATA** newconsanddatas;
   int nnewconsanddatas;
   int snewconsanddatas;
   SCIP_Real* newandcoefs;
   SCIP_Real* oldandcoefs;
   SCIP_Bool* newandnegs;
   SCIP_Bool* oldandnegs;
   CONSANDDATA** consanddatas;
   int nconsanddatas;
   SCIP_CONSDATA* consdata;
   int oldnvars;
   int c;
   int c1;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->hashmap != NULL);
   assert(nandress == 0 || (andress != NULL && andcoefs != NULL));
   assert(!SCIPisInfinity(scip, newlhs));
   assert(!SCIPisInfinity(scip, -newrhs));
   assert(SCIPisLE(scip, newlhs, newrhs));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* sort and-constraints after indices of corresponding and-resultants */
   SCIPsortPtrRealBool((void**)(consdata->consanddatas), consdata->andcoefs, consdata->andnegs, resvarCompWithInactive, consdata->nconsanddatas);

   consanddatas = consdata->consanddatas;
   oldandcoefs = consdata->andcoefs;
   oldandnegs = consdata->andnegs;
   nconsanddatas = consdata->nconsanddatas;
   assert(nconsanddatas == 0 || (consanddatas != NULL && oldandcoefs != NULL));

#ifndef NDEBUG
   /* check that and-resultants are sorted, and coefficents are not zero */
   for( c = nandress - 1; c > 0; --c )
   {
      assert(!SCIPisZero(scip, andcoefs[c]));
      assert(SCIPvarGetIndex(andress[c]) > SCIPvarGetIndex(andress[c - 1]));
   }
   /* check that consanddata objects are sorted due to the index of the corresponding resultants, and coefficents are
    * not zero
    */
   for( c = nconsanddatas - 1; c > 0; --c )
   {
      SCIP_VAR* res1;
      SCIP_VAR* res2;

      assert(consanddatas[c] != NULL);

      if( !consanddatas[c]->istransformed )
         continue;

      assert(!SCIPisZero(scip, oldandcoefs[c]));
      assert(consanddatas[c - 1] != NULL);

      if( !consanddatas[c - 1]->istransformed )
         continue;

      assert(!SCIPisZero(scip, oldandcoefs[c - 1]));

      if( SCIPconsIsDeleted(consanddatas[c]->cons) || SCIPconsIsDeleted(consanddatas[c - 1]->cons) )
         continue;

      assert(consanddatas[c]->cons != NULL);
      res1 = SCIPgetResultantAnd(scip, consanddatas[c]->cons);
      assert(res1 != NULL);
      assert(consanddatas[c - 1]->cons != NULL);
      res2 = SCIPgetResultantAnd(scip, consanddatas[c - 1]->cons);
      assert(res2 != NULL);

      assert(SCIPvarGetIndex(res1) >= SCIPvarGetIndex(res2));
   }
#endif

   snewconsanddatas = nconsanddatas + nandress;

   /* allocate new block memory arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newconsanddatas, snewconsanddatas) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newandcoefs, snewconsanddatas) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newandnegs, snewconsanddatas) );

   nnewconsanddatas = 0;

   /* collect new consanddata objects and update locks and captures */
   for( c = 0, c1 = 0; c < nconsanddatas && c1 < nandress; )
   {
      SCIP_CONS* andcons;
      SCIP_VAR* res1;
      SCIP_VAR* res2;

      assert(consanddatas[c] != NULL);

      /* consanddata object could have been deleted in the last presolving round */
      if( !consanddatas[c]->istransformed )
      {
         ++c;
         consdata->changed = TRUE;
         consdata->upgradetried = FALSE;
         continue;
      }

      andcons = consanddatas[c]->cons;
      assert(andcons != NULL);

      if( andcons == NULL ) /*lint !e774*/
      {
         ++c;
         consdata->changed = TRUE;
         consdata->upgradetried = FALSE;
         continue;
      }
      else if( SCIPconsIsDeleted(andcons) )
      {
         /* remove rounding locks, because the and constraint was deleted  */
         SCIP_CALL( unlockRoundingAndCons(scip, cons, consanddatas[c],
               oldandnegs[c] ? -oldandcoefs[c] : oldandcoefs[c], consdata->lhs, consdata->rhs) );
         ++c;
         consdata->changed = TRUE;
         consdata->upgradetried = FALSE;
         continue;
      }
      assert(andcons != NULL);

      /* get and-resultants of consanddata object in constraint data */
      res1 = SCIPgetResultantAnd(scip, andcons);
      assert(res1 != NULL);
      assert(SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)res1) == consanddatas[c]);

      /* get and-resultants in new corresponding linear constraint */
      res2 = andress[c1];
      assert(res2 != NULL);
      assert(SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)res2) != NULL);

      /* collect new consanddata objects in sorted order due to the variable index of corresponding and-resultants */
      if( SCIPvarGetIndex(res1) < SCIPvarGetIndex(res2) )
      {
	 assert(consanddatas[c]->nuses > 0);
	 --(consanddatas[c]->nuses);

         /* remove old locks */
         SCIP_CALL( removeOldLocks(scip, cons, consanddatas[c], oldandnegs[c] ? -oldandcoefs[c] : oldandcoefs[c],
               consdata->lhs, consdata->rhs) );
         ++c;
         consdata->changed = TRUE;
         consdata->upgradetried = FALSE;
	 consdata->propagated = FALSE;
	 consdata->presolved = FALSE;
      }
      else if( SCIPvarGetIndex(res1) > SCIPvarGetIndex(res2) )
      {
         assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)res2));
         newconsanddatas[nnewconsanddatas] = (CONSANDDATA*) SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)res2);
         newandcoefs[nnewconsanddatas] = andcoefs[c1];
         newandnegs[nnewconsanddatas] = andnegs[c1];
	 ++(newconsanddatas[nnewconsanddatas]->nuses);

         /* add new locks */
         SCIP_CALL( addNewLocks(scip, cons, newconsanddatas[nnewconsanddatas], newandnegs[nnewconsanddatas] ?
               -newandcoefs[nnewconsanddatas] : newandcoefs[nnewconsanddatas], newlhs, newrhs) );
         ++c1;
         consdata->changed = TRUE;
         consdata->upgradetried = FALSE;
         consdata->cliquesadded = FALSE;
	 consdata->propagated = FALSE;
	 consdata->presolved = FALSE;

         ++nnewconsanddatas;
      }
      else
      {
         SCIP_Bool coefsignchanged;
         SCIP_Bool lhschanged;
         SCIP_Bool rhschanged;

         assert(SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)res2) == consanddatas[c]);

         /* copy old consanddata object and new coefficent */
         newconsanddatas[nnewconsanddatas] = consanddatas[c];

         newandcoefs[nnewconsanddatas] = andcoefs[c1];
         newandnegs[nnewconsanddatas] = andnegs[c1];

         if( ((oldandnegs[c] == andnegs[c1]) && !SCIPisEQ(scip, oldandcoefs[c], newandcoefs[c1]))
            || ((oldandnegs[c] != newandnegs[c1]) && !SCIPisEQ(scip, oldandcoefs[c], -newandcoefs[c1])) )
            consdata->upgradetried = FALSE;

         coefsignchanged = (oldandnegs[c] == andnegs[c1]) &&
            ((oldandcoefs[c] < 0 && andcoefs[c1] > 0) || (oldandcoefs[c] > 0 && andcoefs[c1] < 0));
         coefsignchanged = coefsignchanged || ((oldandnegs[c] != andnegs[c1]) &&
            ((oldandcoefs[c] < 0 && andcoefs[c1] < 0) || (oldandcoefs[c] > 0 && andcoefs[c1] > 0)));
         lhschanged = (SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, -newlhs)) || (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, -newlhs))
            || (consdata->lhs < 0 && newlhs > 0) || (consdata->lhs > 0 && newlhs < 0);
         rhschanged = (SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, newrhs)) || (!SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, newrhs))
            || (consdata->rhs < 0 && newrhs > 0) || (consdata->rhs > 0 && newrhs < 0);

         /* update or renew locks */
         if( coefsignchanged || lhschanged || rhschanged || newconsanddatas[nnewconsanddatas]->nnewvars > 0)
         {
            /* renew locks */
            SCIP_CALL( removeOldLocks(scip, cons, newconsanddatas[nnewconsanddatas], oldandnegs[c] ?
                  -oldandcoefs[c] : oldandcoefs[c], consdata->lhs, consdata->rhs) );
            SCIP_CALL( addNewLocks(scip, cons, newconsanddatas[nnewconsanddatas], newandnegs[nnewconsanddatas] ?
                  -newandcoefs[nnewconsanddatas] : newandcoefs[nnewconsanddatas], newlhs, newrhs) );

            consdata->changed = TRUE;
            consdata->upgradetried = FALSE;
            consdata->cliquesadded = FALSE;
            consdata->propagated = FALSE;
            consdata->presolved = FALSE;
         }

         ++c;
         ++c1;
         ++nnewconsanddatas;
      }
   }

   /* add all remaining consanddatas and update locks and captures */
   if( c < nconsanddatas )
   {
      assert(c1 == nandress);

      for( ; c < nconsanddatas; ++c )
      {
         SCIP_CONS* andcons;
#ifndef NDEBUG
         SCIP_VAR* res1;

         assert(consanddatas[c] != NULL);
#endif
         andcons = consanddatas[c]->cons;
#ifndef NDEBUG
         if( andcons != NULL )
         {
            res1 = SCIPgetResultantAnd(scip, andcons);
            assert(res1 != NULL);
            assert(SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)res1) == consanddatas[c]);
         }
#endif
         if( andcons == NULL )
         {
            consdata->changed = TRUE;
            consdata->upgradetried = FALSE;
            continue;
         }

	 assert(consanddatas[c]->nuses > 0);
	 --(consanddatas[c]->nuses);

         /* remove old locks */
         SCIP_CALL( removeOldLocks(scip, cons, consanddatas[c], oldandnegs[c] ? -oldandcoefs[c] : oldandcoefs[c],
               consdata->lhs, consdata->rhs) );
         consdata->changed = TRUE;
         consdata->upgradetried = FALSE;
	 consdata->propagated = FALSE;
	 consdata->presolved = FALSE;
      }
   }
   else if( c1 < nandress )
   {
      for( ; c1 < nandress; ++c1 )
      {
         SCIP_VAR* res2;

         res2 = andress[c1];
         assert(res2 != NULL);
         assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)res2));
         newconsanddatas[nnewconsanddatas] = (CONSANDDATA*) SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)res2);
         newandcoefs[nnewconsanddatas] = andcoefs[c1];
         newandnegs[nnewconsanddatas] = andnegs[c1];
	 ++(newconsanddatas[nnewconsanddatas]->nuses);

         /* add new locks */
         SCIP_CALL( addNewLocks(scip, cons, newconsanddatas[nnewconsanddatas], newandnegs[nnewconsanddatas] ?
               -newandcoefs[nnewconsanddatas] : newandcoefs[nnewconsanddatas], newlhs, newrhs) );

         ++nnewconsanddatas;
         consdata->changed = TRUE;
         consdata->upgradetried = FALSE;
	 consdata->cliquesadded = FALSE;
	 consdata->propagated = FALSE;
	 consdata->presolved = FALSE;
      }
   }
   assert(c == nconsanddatas && c1 == nandress);

   /* delete old and-coefficients and consanddata objects */
   SCIPfreeBlockMemoryArray(scip, &(consdata->andcoefs), consdata->sconsanddatas);
   SCIPfreeBlockMemoryArray(scip, &(consdata->andnegs), consdata->sconsanddatas);
   SCIPfreeBlockMemoryArray(scip, &(consdata->consanddatas), consdata->sconsanddatas);

   if( !SCIPisEQ(scip, consdata->lhs, newlhs) || !SCIPisEQ(scip, consdata->rhs, newrhs) )
   {
      consdata->upgradetried = FALSE;
      consdata->lhs = newlhs;
      consdata->rhs = newrhs;
   }

   consdata->consanddatas = newconsanddatas;
   consdata->andcoefs = newandcoefs;
   consdata->andnegs = newandnegs;
   consdata->nconsanddatas = nnewconsanddatas;
   consdata->sconsanddatas = snewconsanddatas;

   oldnvars = consdata->nlinvars;
   /* update number of linear variables without and-resultants */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &(consdata->nlinvars)) );
   consdata->nlinvars -= nnewconsanddatas;

   if( oldnvars != consdata->nlinvars )
   {
      consdata->changed = TRUE;
      consdata->upgradetried = FALSE;
      consdata->cliquesadded = FALSE;
      consdata->propagated = FALSE;
      consdata->presolved = FALSE;
   }

   /* we need to re-sort and-constraints after indices of corresponding and-resultants, since we might have replaced
    * negated variables
    */
   SCIPsortPtrRealBool((void**)(consdata->consanddatas), consdata->andcoefs, consdata->andnegs, resvarCompWithInactive, consdata->nconsanddatas);

#ifndef NDEBUG
   consanddatas = consdata->consanddatas;
   nconsanddatas = consdata->nconsanddatas;
   assert(nconsanddatas == 0 || consanddatas != NULL);

   /* check that consanddata objects are sorted with respect to the index of the corresponding resultants */
   for( c = nconsanddatas - 1; c > 0; --c )
   {
      SCIP_VAR* res1;
      SCIP_VAR* res2;

      assert(consanddatas[c] != NULL);
      assert(consanddatas[c]->cons != NULL);
      res1 = SCIPgetResultantAnd(scip, consanddatas[c]->cons);
      assert(res1 != NULL);
      assert(consanddatas[c - 1] != NULL);
      assert(consanddatas[c - 1]->cons != NULL);
      res2 = SCIPgetResultantAnd(scip, consanddatas[c - 1]->cons);
      assert(res2 != NULL);

      assert(SCIPvarGetIndex(res1) > SCIPvarGetIndex(res2));
   }
#endif

   return SCIP_OKAY;
}

/** adds cliques of the pseudoboolean constraint to the global clique table */
static
SCIP_RETCODE addCliques(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_Bool*const       cutoff,             /**< pointer to store whether the node can be cut off */
   int*const             naggrvars,          /**< pointer to count the number of aggregated variables */
   int*const             nchgbds             /**< pointer to count the number of performed bound changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR** linvars;
   SCIP_VAR* andres;
   SCIP_VAR* andres2;
   int nlinvars;
   int nandress;
   int c;
   int v2;
   int v1;
   int nchgbdslocal;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(naggrvars != NULL);
   assert(nchgbds != NULL);
   assert(SCIPconsIsActive(cons));

   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   /* if we have no and-constraints left, we should not be here and this constraint should be deleted (only the linaer should survive) */
   assert(consdata->nconsanddatas > 0);

   /* check whether the cliques have already been added */
   if( consdata->cliquesadded )
      return SCIP_OKAY;

   consdata->cliquesadded = TRUE;

   checkConsConsistency(scip, cons);

   /* check standard pointers and sizes */
   assert(consdata->lincons != NULL);
   assert(SCIPconsIsActive(consdata->lincons));
   assert(consdata->linconstype > SCIP_LINEARCONSTYPE_INVALIDCONS);
   assert(consdata->consanddatas != NULL);
   assert(consdata->nconsanddatas > 0);
   assert(consdata->nconsanddatas <= consdata->sconsanddatas);

   /* check number of linear variables */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );
   assert(nvars == consdata->nlinvars + consdata->nconsanddatas);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nvars) );

   /* get variables and coefficients */
   SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, NULL, &nvars) );

   /* calculate all not artificial linear variables and all artificial and-resultants which will be ordered like the
    * 'consanddatas' such that the and-resultant of the and-constraint is the and-resultant in the 'andress' array
    * afterwards
    * @todo should we take into accout the negation status of the cliques?
    */
   SCIP_CALL( getLinVarsAndAndRess(scip, cons, vars, NULL, nvars, linvars, NULL, &nlinvars,
         NULL, NULL, NULL, &nandress) );

   assert(nandress == consdata->nconsanddatas);
   assert(consdata->consanddatas != NULL);

   /* find cliques from linear variable to and-resultant */
   for( c = nandress - 1; c >= 0; --c )
   {
      CONSANDDATA* consanddata;
      SCIP_VAR** andvars;
      int nandvars;

      consanddata = consdata->consanddatas[c];
      assert(consanddata != NULL);

      andres = SCIPgetResultantAnd(scip, consanddata->cons);

      /* choose correct variable array */
      if( consanddata->nnewvars > 0 )
      {
         andvars = consanddata->newvars;
         nandvars = consanddata->nnewvars;
      }
      else
      {
         andvars = consanddata->vars;
         nandvars = consanddata->nvars;
      }

      for( v1 = nandvars - 1; v1 >= 0; --v1 )
      {
         SCIP_VAR* var1;
         SCIP_Bool values[2];

         var1 = andvars[v1];
         if( !SCIPvarIsActive(var1) && (!SCIPvarIsNegated(var1) || !SCIPvarIsActive(SCIPvarGetNegationVar(var1))) )
            continue;

         /* get active counterpart to check for common cliques */
         if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED )
         {
            var1 = SCIPvarGetNegationVar(var1);
            values[0] = FALSE;
         }
         else
            values[0] = TRUE;

         for( v2 = nlinvars - 1; v2 >= 0; --v2 )
         {
            SCIP_VAR* var2;

            var2 = linvars[v2];
            if( !SCIPvarIsActive(var2) && (!SCIPvarIsNegated(var2) || !SCIPvarIsActive(SCIPvarGetNegationVar(var2))) )
               continue;

            /* get active counterpart to check for common cliques */
            if( SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED )
            {
               var2 = SCIPvarGetNegationVar(var2);
               values[1] = FALSE;
            }
            else
               values[1] = TRUE;

            /* if variable in and-constraint1 is the negated variable of a normal linear variable, than we can add a
             * clique between the and-resultant and the normal linear variable, negated variables are not save in
             * cliquetables
             *
             * set r_1 = var1 * z; (z is some product)
             * var1 == ~var2
             *
             * if:
             * var1 + ~var1 <= 1;          r_1
             *    0 +     1 <= 1             0   \
             *    1 +     0 <= 1   ==>  1 or 0    >   ==>    r_1 + var2 <= 1
             *    0 +     0 <= 1             0   /
             */
            if( values[0] != values[1] && var1 == var2 )
            {
               SCIP_CONS* newcons;
               SCIP_VAR* clqvars[2];
               char consname[SCIP_MAXSTRLEN];

               clqvars[0] = andres;
               clqvars[1] = values[1] ? var2 : SCIPvarGetNegatedVar(var2);
               assert(clqvars[1] != NULL);

               /* @todo: check whether it is better to only add the clique or to add the setppc constraint or do both */

               /* add clique */
               SCIP_CALL( SCIPaddClique(scip, clqvars, NULL, 2, FALSE, cutoff, &nchgbdslocal) );
               if( *cutoff )
                  goto TERMINATE;

	       *nchgbds += nchgbdslocal;

               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_clq_%s_%s", SCIPconsGetName(cons), SCIPvarGetName(clqvars[0]), SCIPvarGetName(clqvars[1]) );
               SCIP_CALL( SCIPcreateConsSetpack(scip, &newcons, consname, 2, clqvars,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     FALSE, SCIPconsIsPropagated(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIPdebugMsg(scip, "added a clique/setppc constraint <%s> \n", SCIPconsGetName(newcons));
               SCIPdebugPrintCons(scip, newcons, NULL);

               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            }
            /* if a variable in an and-constraint is in a clique with another normal linear variable, we can add the
             * clique between the linear variable and the and-resultant
             *
             * set r_1 = var1 * z; (z is some product)
             *
             * if:
             * var1 + var2 <= 1;          r_1
             *    0 +    1 <= 1             0   \
             *    1 +    0 <= 1   ==>  1 or 0    >   ==>    r_1 + var2 <= 1
             *    0 +    0 <= 1             0   /
             */
            if( (var1 != var2) && SCIPvarsHaveCommonClique(var1, values[0], var2, values[1], TRUE) )
            {
               SCIP_CONS* newcons;
               SCIP_VAR* clqvars[2];
               char consname[SCIP_MAXSTRLEN];

               clqvars[0] = andres;
               clqvars[1] = values[1] ? var2 : SCIPvarGetNegatedVar(var2);
               assert(clqvars[1] != NULL);

               /* @todo: check whether it is better to only add the clique or to add the setppc constraint or do both */

               /* add clique */
               SCIP_CALL( SCIPaddClique(scip, clqvars, NULL, 2, FALSE, cutoff, &nchgbdslocal) );
               if( *cutoff )
                  goto TERMINATE;

	       *nchgbds += nchgbdslocal;

               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_clq_%s_%s", SCIPconsGetName(cons), SCIPvarGetName(clqvars[0]), SCIPvarGetName(clqvars[1]) );
               SCIP_CALL( SCIPcreateConsSetpack(scip, &newcons, consname, 2, clqvars,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     FALSE, SCIPconsIsPropagated(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIPdebugMsg(scip, "added a clique/setppc constraint <%s> \n", SCIPconsGetName(newcons));
               SCIPdebugPrintCons(scip, newcons, NULL);

               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            }
         }
      }
   }

   /* find cliques over variables which are in different and-constraints */
   for( c = nandress - 1; c > 0; --c )
   {
      CONSANDDATA* consanddata1;
      CONSANDDATA* consanddata2;
      SCIP_VAR** andvars1;
      int nandvars1;
      SCIP_VAR** andvars2;
      int nandvars2;

      consanddata1 = consdata->consanddatas[c];
      assert(consanddata1 != NULL);
      consanddata2 = consdata->consanddatas[c - 1];
      assert(consanddata2 != NULL);

      andres = SCIPgetResultantAnd(scip, consanddata1->cons);
      andres2 = SCIPgetResultantAnd(scip, consanddata2->cons);

      /* choose correct variable array of consanddata object 1 */
      if( consanddata1->nnewvars > 0 )
      {
         andvars1 = consanddata1->newvars;
         nandvars1 = consanddata1->nnewvars;
      }
      else
      {
         andvars1 = consanddata1->vars;
         nandvars1 = consanddata1->nvars;
      }

      /* choose correct variable array of consanddata object 2 */
      if( consanddata2->nnewvars > 0 )
      {
         andvars2 = consanddata2->newvars;
         nandvars2 = consanddata2->nnewvars;
      }
      else
      {
         andvars2 = consanddata2->vars;
         nandvars2 = consanddata2->nvars;
      }

      /* compare both terms for finding new aggregated variables and new cliques */
      for( v1 = nandvars1 - 1; v1 >= 0; --v1 )
      {
         SCIP_VAR* var1;
         SCIP_Bool values[2];

         var1 = andvars1[v1];
         if( !SCIPvarIsActive(var1) && (!SCIPvarIsNegated(var1) || !SCIPvarIsActive(SCIPvarGetNegationVar(var1))) )
            continue;

         /* get active counterpart to check for common cliques */
         if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED )
         {
            var1 = SCIPvarGetNegationVar(var1);
            values[0] = FALSE;
         }
         else
            values[0] = TRUE;

         for( v2 = nandvars2 - 1; v2 >= 0; --v2 )
         {
            SCIP_VAR* var2;

            var2 = andvars2[v2];
            if( !SCIPvarIsActive(var2) && (!SCIPvarIsNegated(var2) || !SCIPvarIsActive(SCIPvarGetNegationVar(var2))) )
               continue;

            /* get active counterpart to check for common cliques */
            if( SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED )
            {
               var2 = SCIPvarGetNegationVar(var2);
               values[1] = FALSE;
            }
            else
               values[1] = TRUE;

            /* if a variable in and-constraint1 is the negated variable of a variable in and-constraint2, than we can
             * add a clique between both and-resultants, negated variables are not save in cliquetables
             *
             * set r_1 = var1 * z_1; (z_1 is some product)
             * set r_2 = var2 * z_2; (z_2 is some product)
             * var1 == ~var2
             *
             * if:
             * var1 + ~var1 <= 1;          r_1     r_2
             *    0 +     1 <= 1             0  1 or 0   \
             *    1 +     0 <= 1   ==>  1 or 0       0    >   ==>    r_1 + r_2 <= 1
             *    0 +     0 <= 1             0       0   /
             */
            if( values[0] != values[1] && var1 == var2 )
            {
               SCIP_CONS* newcons;
               SCIP_VAR* clqvars[2];
               char consname[SCIP_MAXSTRLEN];

               clqvars[0] = andres;
               clqvars[1] = andres2;

               /* @todo: check whether it is better to only add the clique or to add the setppc constraint or do both */

               /* add clique */
               SCIP_CALL( SCIPaddClique(scip, clqvars, NULL, 2, FALSE, cutoff, &nchgbdslocal) );
               if( *cutoff )
                  goto TERMINATE;

	       *nchgbds += nchgbdslocal;

               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_clq_%s_%s", SCIPconsGetName(cons), SCIPvarGetName(clqvars[0]), SCIPvarGetName(clqvars[1]) );
               SCIP_CALL( SCIPcreateConsSetpack(scip, &newcons, consname, 2, clqvars,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     FALSE, SCIPconsIsPropagated(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIPdebugMsg(scip, "added a clique/setppc constraint <%s> \n", SCIPconsGetName(newcons));
               SCIPdebugPrintCons(scip, newcons, NULL);

               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            }
            /* if a variable in an and-constraint is in a clique with a variable in another and-constraint, we can add
             * the clique between both and-resultant
             *
             * let r_1 = var1 * z_1; (z_1 is some product)
             * let r_2 = var2 * z_2; (z_2 is some product)
             *
             * if:
             * var1 + var2 <= 1;          r_1     r_2
             *    0 +    1 <= 1             0  1 or 0   \
             *    1 +    0 <= 1   ==>  1 or 0       0    >   ==>    r_1 + r_2 <= 1
             *    0 +    0 <= 1             0       0   /
             */
            else if( SCIPvarsHaveCommonClique(var1, values[0], var2, values[1], TRUE) && (var1 != var2) )
            {
               SCIP_CONS* newcons;
               SCIP_VAR* clqvars[2];
               char consname[SCIP_MAXSTRLEN];

               clqvars[0] = andres;
               clqvars[1] = values[1] ? var2 : SCIPvarGetNegatedVar(var2);
               assert(clqvars[1] != NULL);

               /* @todo: check whether it is better to only add the clique or to add the setppc constraint or do both */

               /* add clique */
               SCIP_CALL( SCIPaddClique(scip, clqvars, NULL, 2, FALSE, cutoff, &nchgbdslocal) );
               if( *cutoff )
                  goto TERMINATE;

	       *nchgbds += nchgbdslocal;

               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_clq_%s_%s", SCIPconsGetName(cons), SCIPvarGetName(clqvars[0]), SCIPvarGetName(clqvars[1]) );
               SCIP_CALL( SCIPcreateConsSetpack(scip, &newcons, consname, 2, clqvars,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     FALSE, SCIPconsIsPropagated(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIPdebugMsg(scip, "added a clique/setppc constraint <%s> \n", SCIPconsGetName(newcons));
               SCIPdebugPrintCons(scip, newcons, NULL);

               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            }
         }
      }
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** propagation method for pseudoboolean constraints */
static
SCIP_RETCODE propagateCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< knapsack constraint */
   SCIP_Bool*const       cutoff,             /**< pointer to store whether the node can be cut off */
   int*const             ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);

   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lincons != NULL);

   /* if linear constraint is redundant, than pseudoboolean constraint is redundant too */
   if( SCIPconsIsDeleted(consdata->lincons) )
   {
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      ++(*ndelconss);
   }

   /* check if the constraint was already propagated */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** update and-constraint flags due to pseudoboolean constraint flags */
static
SCIP_RETCODE updateAndConss(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   )
{
   CONSANDDATA** consanddatas;
   int nconsanddatas;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consanddatas = consdata->consanddatas;
   nconsanddatas = consdata->nconsanddatas;
   assert(nconsanddatas == 0 || consanddatas != NULL);

   if( !SCIPconsIsActive(cons) )
      return SCIP_OKAY;

   /* release and-constraints and change check flag of and-constraint */
   for( c = nconsanddatas - 1; c >= 0; --c )
   {
      SCIP_CONS* andcons;

      assert(consanddatas[c] != NULL);

      if( !consanddatas[c]->istransformed )
         continue;

      andcons = consanddatas[c]->cons;
      assert(andcons != NULL);

      SCIP_CALL( SCIPsetConsChecked(scip, andcons, SCIPconsIsChecked(cons)) );
   }

   return SCIP_OKAY;
}

/** delete unused information in constraint handler data */
static
SCIP_RETCODE correctConshdlrdata(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*const conshdlrdata,     /**< pseudoboolean constraint handler data */
   int*const             ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   CONSANDDATA** allconsanddatas;
   CONSANDDATA* consanddata;
   int c;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(ndelconss != NULL);

   allconsanddatas = conshdlrdata->allconsanddatas;
   assert(allconsanddatas != NULL);
   assert(conshdlrdata->nallconsanddatas > 0);
   assert(conshdlrdata->nallconsanddatas <= conshdlrdata->sallconsanddatas);

   for( c = conshdlrdata->nallconsanddatas - 1; c >= 0; --c )
   {
      SCIP_VAR** tmpvars;
      int stmpvars;
      SCIP_CONS* cons;
      int v;

      consanddata = allconsanddatas[c];

      assert(consanddata->nvars == 0 || (consanddata->vars != NULL && consanddata->svars > 0));
      assert(consanddata->nnewvars == 0 || (consanddata->newvars != NULL && consanddata->snewvars > 0));

      if( !consanddata->istransformed )
      {
         assert(consanddata->vars == NULL || consanddata->origcons != NULL);
         assert(consanddata->nvars == 0 || consanddata->origcons != NULL);
         assert(consanddata->svars == 0 || consanddata->origcons != NULL);
         assert(consanddata->newvars == NULL);
         assert(consanddata->nnewvars == 0);
         assert(consanddata->snewvars == 0);

         continue;
      }

      /* if no variables are left, delete variables arrays */
      if( consanddata->nvars == 0 )
      {
         SCIP_VAR* resvar = SCIPgetResultantAnd(scip, consanddata->cons);

         /* if we have no old variables, than also no new variables */
         assert(consanddata->nnewvars == 0);
         assert(consanddata->nuses > 0);
         assert(resvar != NULL);

         /* delete and-constraint */
         SCIP_CALL( SCIPdelCons(scip, consanddata->cons) );
         ++(*ndelconss);

         SCIP_CALL( transformToOrig(scip, consanddata, conshdlrdata) );

         /* release and-constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &consanddata->cons) );
         consanddata->nuses = 0;

         /* remove consanddata from hashtable, if it existed only in transformed space */
         if( consanddata->origcons == NULL )
         {
            assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddata));
            SCIP_CALL( SCIPhashtableRemove(conshdlrdata->hashtable, (void*)consanddata) );
         }
         assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)resvar));
         SCIP_CALL( SCIPhashmapRemove(conshdlrdata->hashmap, (void*)resvar) );

         continue;
      }

      /* the consanddata object is not used anymore, so extract the and constraint and delete other data */
      if( consanddata->nuses == 0 )
      {
	 SCIP_Bool looseorcolumn;
	 SCIP_VARSTATUS varstatus;

         if( consanddata->cons == NULL )
         {
            assert(!consanddata->istransformed || consanddata->noriguses > 0);
            assert((consanddata->noriguses > 0) == (consanddata->origcons != NULL));
            assert(consanddata->vars == NULL || consanddata->origcons != NULL);
            assert(consanddata->nvars == 0 || consanddata->origcons != NULL);
            assert(consanddata->svars == 0 || consanddata->origcons != NULL);
            assert(consanddata->newvars == NULL);
            assert(consanddata->nnewvars == 0);
            assert(consanddata->snewvars == 0);

            continue;
         }

         SCIP_CALL( transformToOrig(scip, consanddata, conshdlrdata) );

	 varstatus = SCIPvarGetStatus(SCIPgetResultantAnd(scip, consanddata->cons));
	 looseorcolumn = (varstatus == SCIP_VARSTATUS_LOOSE || varstatus == SCIP_VARSTATUS_COLUMN);

#if 1
	 /* @note  due to aggregations or fixings the resultant may need to be propagated later on, so we can only
	  *        delete the and-constraint if the resultant is of column or loose status
	  *        and is not an active variable of another (multi-)aggregated/negated variable
	  */
	 if( looseorcolumn )
	 {
	    SCIP_Bool del = TRUE;
            int nfixedvars = SCIPgetNFixedVars(scip);

	    if( nfixedvars > 0 )
	    {
	       SCIP_VAR** fixedvars;
	       SCIP_VAR** scipfixedvars;
	       SCIP_VAR** activevars = NULL;
	       SCIP_Real* activescalars = NULL;
	       SCIP_Real activeconstant;
	       int nactivevars;
	       int requiredsize;
	       int pos;
	       int w;

	       scipfixedvars = SCIPgetFixedVars(scip);
	       SCIP_CALL( SCIPduplicateBufferArray(scip, &fixedvars, scipfixedvars, nfixedvars) );

	       SCIPvarsGetProbvar(fixedvars, nfixedvars);

               /* all inactive variables have a loose, column, fixed or multi-aggregated variable as counterpart,
                * for multi-aggregated variables, we need to check all active representatives
                * @todo move this outside of the consanddata loop
	        */
	       for( w = nfixedvars - 1; w >= 0; --w )
	       {
                  if( SCIPvarGetStatus(fixedvars[w]) == SCIP_VARSTATUS_MULTAGGR )
                  {
                     if( activevars == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &activevars, SCIPgetNVars(scip)) );
                        SCIP_CALL( SCIPallocBufferArray(scip, &activescalars, SCIPgetNVars(scip)) );
                     }
                     assert(activevars != NULL);
                     assert(activescalars != NULL);

                     activevars[0] = fixedvars[w];
                     activescalars[0] = 1.0;
                     activeconstant = 0.0;
                     nactivevars = 1;

                     SCIP_CALL( SCIPgetProbvarLinearSum(scip, activevars, activescalars, &nactivevars, SCIPgetNVars(scip),
                           &activeconstant, &requiredsize, TRUE) );
                     assert(requiredsize <= SCIPgetNVars(scip));

                     if( nactivevars == 0 )
                     {
                        --nfixedvars;
                        fixedvars[w] = fixedvars[nfixedvars];
                     }
                     else
                     {
                        fixedvars[w] = activevars[0];

                        if( nactivevars > 1 )
                        {
                           int i;

                           SCIP_CALL( SCIPreallocBufferArray(scip, &fixedvars, nfixedvars + nactivevars - 1) );
                           for( i = 1; i < nactivevars; ++i )
                           {
                              assert(SCIPvarGetStatus(activevars[i]) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(activevars[i]) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(activevars[i]) == SCIP_VARSTATUS_FIXED);
                              fixedvars[nfixedvars] = activevars[i];
                              ++nfixedvars;
                           }
                        }
                     }
                  }

		  assert(SCIPvarGetStatus(fixedvars[w]) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(fixedvars[w]) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(fixedvars[w]) == SCIP_VARSTATUS_FIXED);
	       }

               if( activevars != NULL )
               {
                  SCIPfreeBufferArray(scip, &activevars);
                  SCIPfreeBufferArray(scip, &activescalars);
               }

	       SCIPsortPtr((void**)fixedvars, SCIPvarComp, nfixedvars);

	       if( SCIPsortedvecFindPtr((void**)fixedvars, SCIPvarComp, SCIPgetResultantAnd(scip, consanddata->cons), nfixedvars, &pos) )
		  del = FALSE;

	       SCIPfreeBufferArray(scip, &fixedvars);
	    }

	    if( del )
	    {
	       SCIP_CALL( SCIPdelCons(scip, consanddata->cons) );
	    }
	 }
#endif

	 if( !SCIPconsIsDeleted(consanddata->cons) )
	 {
	    /* change flags */
	    if( !looseorcolumn )
	    {
	       SCIP_CALL( SCIPsetConsInitial(scip, consanddata->cons, FALSE) );
#if 0
	       SCIP_CALL( SCIPsetConsSeparated(scip, consanddata->cons, FALSE) );
#endif
	    }
	    SCIP_CALL( SCIPsetConsChecked(scip, consanddata->cons, TRUE) );
	 }

         /* remove consanddata from hashtable, if it existed only in transformed space */
         if( consanddata->origcons == NULL )
         {
            assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddata));
            SCIP_CALL( SCIPhashtableRemove(conshdlrdata->hashtable, (void*)consanddata) );
         }
         assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddata->cons)));
         SCIP_CALL( SCIPhashmapRemove(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddata->cons)) );

         SCIP_CALL( SCIPreleaseCons(scip, &(consanddata->cons)) );
         ++(*ndelconss);

         continue;
      }

      cons = consanddata->cons;
      assert(cons != NULL);

      /* if and-constraint is deleted, delete variables arrays */
      if( SCIPconsIsDeleted(cons) )
      {
         SCIP_VAR* resvar = SCIPgetResultantAnd(scip, consanddata->cons);

         assert(consanddata->nuses > 0);
         assert(resvar != NULL);

         SCIP_CALL( transformToOrig(scip, consanddata, conshdlrdata) );

         /* release and-constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &consanddata->cons) );
         consanddata->nuses = 0;

         /* remove consanddata from hashtable, if it existed only in transformed space */
         if( consanddata->origcons == NULL )
         {
            assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddata));
            SCIP_CALL( SCIPhashtableRemove(conshdlrdata->hashtable, (void*)consanddata) );
         }
         assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)resvar));
         SCIP_CALL( SCIPhashmapRemove(conshdlrdata->hashmap, (void*)resvar) );

         continue;
      }

      /* if no new variables exist, we do not need to do anything here */
      if( consanddata->nnewvars == 0 )
         continue;

      tmpvars = consanddata->vars;
      /* release all variables */
      for( v = consanddata->nvars - 1; v >= 0; --v )
      {
         /* in original problem the variables was already deleted */
         assert(tmpvars[v] != NULL);
         SCIP_CALL( SCIPreleaseVar(scip, &tmpvars[v]) );
      }

      /* exchange newvars with old vars array */
      tmpvars = consanddata->vars;
      stmpvars = consanddata->svars;
      consanddata->vars = consanddata->newvars;
      consanddata->svars = consanddata->snewvars;
      consanddata->nvars = consanddata->nnewvars;
      consanddata->newvars = tmpvars;
      consanddata->snewvars = stmpvars;
      /* reset number of variables in newvars array */
      consanddata->nnewvars = 0;
   }

   return SCIP_OKAY;
}

/** update the uses counter of consandata objects which are used in pseudoboolean constraint, that were deleted and
 *  probably delete and-constraints
 */
static
SCIP_RETCODE updateConsanddataUses(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONSHDLRDATA*const conshdlrdata,     /**< pseudoboolean constraint handler data */
   int*const             ndelconss           /**< pointer to store number of deleted constraints */
   )
{
   CONSANDDATA** consanddatas;
   int nconsanddatas;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);
   assert(ndelconss != NULL);

   /* can only be called when constraint was deleted */
   assert(SCIPconsIsDeleted(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consanddatas = consdata->consanddatas;
   nconsanddatas = consdata->nconsanddatas;
   assert(nconsanddatas > 0 && consanddatas != NULL);

   /* remove old locks */
   if( nconsanddatas > 0 )
   {
      assert(consdata->andcoefs != NULL);

      for( c = nconsanddatas - 1; c >= 0; --c )
      {
         CONSANDDATA* consanddata;

         consanddata = consanddatas[c];
         assert(consanddata != NULL);

         if( !consanddata->istransformed )
            continue;

         SCIP_CALL( removeOldLocks(scip, cons, consanddata, consdata->andcoefs[c], consdata->lhs, consdata->rhs) );
      }
   }

   /* correct consandata usage counters and data */
   for( c = nconsanddatas - 1; c >= 0; --c )
   {
      CONSANDDATA* consanddata;

      consanddata = consanddatas[c];
      assert(consanddata != NULL);
      assert(consanddatas[c]->istransformed);

      assert(consanddata->nuses > 0);

      if( consanddata->nuses > 0 )
         --(consanddata->nuses);

      /* if data object is not used anymore, delete it */
      if( consanddata->nuses == 0 )
      {
         SCIP_VAR* resvar;
	 SCIP_VARSTATUS varstatus;
	 SCIP_Bool looseorcolumn;

         SCIP_CALL( transformToOrig(scip, consanddata, conshdlrdata) );

         resvar = SCIPgetResultantAnd(scip, consanddata->cons);
         assert(resvar != NULL);

	 varstatus = SCIPvarGetStatus(resvar);
	 looseorcolumn = (varstatus == SCIP_VARSTATUS_LOOSE || varstatus == SCIP_VARSTATUS_COLUMN);

#if 1
	 /* @note  due to aggregations or fixings the resultant may need to be propagated later on, so we can only
	  *        delete the and-constraint if the resultant is of column or loose status
          *        and is not an active variable of another (multi-)aggregated/negated variable
	  */
	 if( looseorcolumn )
	 {
	    SCIP_Bool delcons = TRUE;
#if 0
	    const int nfixedvars = SCIPgetNFixedVars(scip);

	    if( nfixedvars > 0 )
	    {
	       SCIP_VAR** fixedvars;
               SCIP_Bool foundmultiaggrvar = FALSE; /* workaround for multi-aggregated variables */
	       int pos;
	       int w;

	       SCIP_CALL( SCIPduplicateBufferArray(scip, &fixedvars, SCIPgetFixedVars(scip), nfixedvars) );

	       SCIPvarsGetProbvar(fixedvars, nfixedvars);

	       /* all inactive variables have a loose, column, fixed or multi-aggregated variable as counterpart, but
		* because we have only binary variables (in pseudobbolean contest) there should also be no
		* multi-aggregated variable
		*
		* @todo for multi-aggregated variables check also all active representatives for this resultant
	        */
	       for( w = nfixedvars - 1; w >= 0; --w )
	       {
                  if( SCIPvarGetStatus(fixedvars[w]) == SCIP_VARSTATUS_MULTAGGR )
                     foundmultiaggrvar = TRUE;
                  else
                     assert(SCIPvarGetStatus(fixedvars[w]) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(fixedvars[w]) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(fixedvars[w]) == SCIP_VARSTATUS_FIXED);
	       }

	       SCIPsortPtr((void**)fixedvars, SCIPvarComp, nfixedvars);

               if( foundmultiaggrvar )
		  delcons = FALSE;
	       else if( SCIPsortedvecFindPtr((void**)fixedvars, SCIPvarComp, resvar, nfixedvars, &pos) )
		  delcons = FALSE;

	       SCIPfreeBufferArray(scip, &fixedvars);
	    }
#endif
            /* we can only delete and constraints if the resultant is an artificial variable and also active, because
             * then the assigned value is not of interest and the artificial and constraint does not need to be
             * fulfilled
             *
             * if this variable is not such an artificial variable we need the IRRELEVANT vartype which should be the
             * correct way to fix this
             */
	    if( delcons
#if 0
               && strlen(SCIPvarGetName(resvar)) > strlen(ARTIFICIALVARNAMEPREFIX) &&
               strncmp(SCIPvarGetName(resvar)+2, ARTIFICIALVARNAMEPREFIX, strlen(ARTIFICIALVARNAMEPREFIX)) == 0
#endif
               ) /*lint !e774*/
	    {
               assert(!SCIPconsIsChecked(consanddata->cons));
	       SCIP_CALL( SCIPdelCons(scip, consanddata->cons) );
	    }
	 }
#endif

#if 0
	 /* @note  due to aggregations or fixings the resultant may need to be propagated later on, so we can only
	  *        delete the and-constraint if the resultant is of column or loose status
	  *        and is not an active variable of another (multi-)aggregated/negated variable
	  */
	 if( looseorcolumn )
	 {
	    SCIP_CALL( SCIPdelCons(scip, consanddata->cons) );
	 }
#endif

	 if( !SCIPconsIsDeleted(consanddata->cons) )
	 {
	    /* change flags */
	    if( !looseorcolumn )
	    {
	       SCIP_CALL( SCIPsetConsInitial(scip, consanddata->cons, FALSE) );
#if 0
	       SCIP_CALL( SCIPsetConsSeparated(scip, consanddata->cons, FALSE) );
#endif
	    }
	    SCIP_CALL( SCIPsetConsChecked(scip, consanddata->cons, TRUE) );
	 }

         /* remove consanddata from hashtable, if it existed only in transformed space */
         if( consanddata->origcons == NULL )
         {
            assert(SCIPhashtableExists(conshdlrdata->hashtable, (void*)consanddata));
            SCIP_CALL( SCIPhashtableRemove(conshdlrdata->hashtable, (void*)consanddata) );
         }
         assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddata->cons)));
         SCIP_CALL( SCIPhashmapRemove(conshdlrdata->hashmap, (void*)SCIPgetResultantAnd(scip, consanddata->cons)) );

         SCIP_CALL( SCIPreleaseCons(scip, &(consanddata->cons)) );
         ++(*ndelconss);
      }
   }

   consdata->nconsanddatas = 0;

   return SCIP_OKAY;
}


/* maximal number to enumerate solutions for one pseudoboolean constraint to check for an upgrade to an XOR constraint */
#define MAXNVARS 10 /* note that this cannot be bigger than 31 */

/** calculate result for a given pseudoboolean constraint with given values, this is used to decide whether a
 *  pseudoboolean constraint can be upgrade to an XOR constraint
 */
static
SCIP_RETCODE checkSolution(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< all variables which occur */
   int const             nvars,              /**< number of all variables which appear in the pseudoboolean
					      *   constraint
					      */
   SCIP_Bool*const       values,             /**< values of all variables which appear in the pseudoboolean
					      *   constraint
					      */
   SCIP_VAR**const       linvars,            /**< linear variables */
   SCIP_Real*const       lincoefs,           /**< linear coefficients */
   int const             nlinvars,           /**< number of linear variables */
   SCIP_Real const       constant,           /**< offset to the linear part */
   SCIP_Real const       side,               /**< side of pseudoboolean constraint */
   CONSANDDATA**const    consanddatas,       /**< all consanddata objects in a constraint */
   SCIP_Real*const       consanddatacoefs,   /**< nonlinear coefficients */
   SCIP_Bool*const       consanddatanegs,    /**< negation status of and resultants in pseudo-boolean constraint */
   int const             nconsanddatas,      /**< number of all consanddata objects */
   int const             cnt,                /**< number of variables set to 1 */
   int*const             xortype             /**< pointer to save the possible xor type if a solution was valid and does
					      *   not violate the old xortype
					      */
   )
{
   CONSANDDATA* consanddata;
   SCIP_VAR** termvars;
   SCIP_VAR** repvars;
   int ntermvars;
   SCIP_Bool* negated;
   SCIP_Real value;
   int pos;
   int v;
   int c;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(values != NULL);
   assert(linvars != NULL || nlinvars == 0);
   assert(lincoefs != NULL || nlinvars == 0);
   assert(nvars >= nlinvars);
   assert(SCIPisEQ(scip, side, 1.0) || SCIPisZero(scip, side));
   assert(consanddatas != NULL);
   assert(consanddatacoefs != NULL);
   assert(nconsanddatas > 0);
   assert(*xortype >= -1 && *xortype <= 1);

   /* order the variables after index, to compare them easier */
   SCIPsortPtr((void**)linvars, SCIPvarCompActiveAndNegated, nlinvars);
   SCIPsortPtr((void**)vars, SCIPvarCompActiveAndNegated, nvars);

   value = constant;
   for( v = nlinvars - 1; v >= 0; --v )
   {
      if( SCIPsortedvecFindPtr((void**)vars, SCIPvarCompActiveAndNegated, linvars[v], nvars, &pos) ) /*lint !e613*/
      {
	 if( values[pos] )
	    value += lincoefs[v]; /*lint !e613*/
      }
      else
      {
	 /* this cannot happen, all linear variables should be a part of 'vars' */
         SCIPABORT();

	 *xortype = -1;  /*lint !e527*/
	 return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &repvars, MAXNVARS) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negated, MAXNVARS) );

   for( c = nconsanddatas - 1; c >= 0; --c )
   {
      SCIP_Bool val = TRUE;

      consanddata = consanddatas[c];
      assert(consanddata != NULL);
      assert(consanddata->istransformed);

      /* choose correct variable array to add locks for, we only add locks for now valid variables */
      if( consanddata->nnewvars > 0 )
      {
	 termvars = consanddata->newvars;
	 ntermvars = consanddata->nnewvars;
      }
      else
      {
	 termvars = consanddata->vars;
	 ntermvars = consanddata->nvars;
      }
      assert(ntermvars > 0 && termvars != NULL);

      BMSclearMemoryArray(negated, MAXNVARS);

      /* get linear active representation */
      SCIP_CALL( SCIPgetBinvarRepresentatives(scip, ntermvars, termvars, repvars, negated) );
      SCIPsortPtrBool((void**)repvars, negated, SCIPvarCompActiveAndNegated, ntermvars);

      for( v = ntermvars - 1; v >= 0; --v )
      {
	 SCIP_VAR* var;

	 assert(!negated[v] || (SCIPvarIsNegated(repvars[v]) && SCIPvarGetNegatedVar(repvars[v]) != NULL));

	 var = ( negated[v] ? SCIPvarGetNegationVar(repvars[v]) : repvars[v]);
	 if( SCIPsortedvecFindPtr((void**)vars, SCIPvarCompActiveAndNegated, var, nvars, &pos) )
	 {
	    if( (negated[v] && values[pos]) || (!negated[v] && !values[pos]) )
	    {
	       val = FALSE;
	       break;
	    }
	 }
	 else
	 {
	    /* this cannot happen, all non-linear variables should be a part of 'vars' */
	    SCIPABORT();

	    *xortype = -1; /*lint !e527*/
	    goto TERMINATE;
	 }
      }

      if( val != consanddatanegs[c] )
	 value += consanddatacoefs[c];
   }

   if( SCIPisEQ(scip, value, side) )
   {
      /* first solution is checked, so determine the possible xor upgrade */
      if( *xortype == -1 )
      {
	 if( cnt % 2 == 0 )
	    *xortype = 0;
	 else
	    *xortype = 1;
      }
      /* check if this solution does not fit in all possible xor solutions */
      else if( *xortype == 1 && cnt % 2 == 0 )
	 *xortype = -1;
      else if( *xortype == 0 && cnt % 2 == 1 )
	 *xortype = -1;
   }
   else
   {
      /* first not-solution is checked, so determine the possible xor upgrade */
      if( *xortype == -1 )
      {
	 if( cnt % 2 == 0 )
	    *xortype = 1;
	 else
	    *xortype = 0;
      }
      /* check if this had to be a solution for an upgrade to an xor */
      else if( *xortype == 1 && cnt % 2 == 1 )
	 *xortype = -1;
      else if( *xortype == 0 && cnt % 2 == 0 )
	 *xortype = -1;
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &negated);
   SCIPfreeBufferArray(scip, &repvars);

   return SCIP_OKAY;
}

/** try upgrading pseudoboolean linear constraint to an XOR constraint and/or remove possible and-constraints
 *
 *  @note An XOR(x_1,..,x_n) = 1 <=> XOR(x1,..,~x_j,..,x_n) = 0, for j in {1,..,n}, which is not yet checked while
 *  trying to upgrade
 */
static
SCIP_RETCODE tryUpgradingXor(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONSHDLRDATA*const conshdlrdata,     /**< pseudoboolean constraint handler data */
   int*const             ndelconss,          /**< pointer to store number of deleted constraints */
   int*const             naddconss,          /**< pointer to count number of added constraints */
   int*const             nfixedvars,         /**< pointer to store number of fixed variables */
   int*const             nchgcoefs,          /**< pointer to store number of changed coefficients constraints */
   int*const             nchgsides,          /**< pointer to store number of changed sides constraints */
   SCIP_Bool*const       cutoff              /**< pointer to store if a cutoff happened */
   )
{
   SCIP_CONSDATA* consdata;
   CONSANDDATA** consanddatas;
   int nconsanddatas;
   CONSANDDATA* consanddata;
   SCIP_VAR** allvars;
   SCIP_Real* allcoefs;
   int nallvars;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlinvars;
   SCIP_Real* andcoefs;
   SCIP_Bool* andnegs;
   int nandress;
   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR** repvars;
   SCIP_Bool* negated;
   SCIP_VAR** activelinvars;
   SCIP_Bool* values;
   SCIP_CONS* lincons;
   SCIP_CONS* newcons;
   char newname[SCIP_MAXSTRLEN];
   SCIP_Real constant;
   int requiredsize;
   int firstnlinvars;
   int oldnlinvars;
   int xortype;
   int v;
   int v1;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);
   assert(ndelconss != NULL);
   assert(nfixedvars != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(cutoff != NULL);
   assert(SCIPconsIsActive(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consanddatas = consdata->consanddatas;
   andcoefs = consdata->andcoefs;
   andnegs = consdata->andnegs;
   nconsanddatas = consdata->nconsanddatas;
   assert(nconsanddatas > 0 && consanddatas != NULL);

   assert(consdata->lincons != NULL);
   assert(consdata->linconstype == SCIP_LINEARCONSTYPE_LINEAR || consdata->linconstype == SCIP_LINEARCONSTYPE_SETPPC);

   /* only equations can be updated */
   if( !SCIPisEQ(scip, consdata->lhs, consdata->rhs) || (!SCIPisEQ(scip, consdata->lhs, 1.0) && !SCIPisZero(scip, consdata->lhs)) )
      return SCIP_OKAY;

   assert(consanddatas[0] != NULL);
   assert(consanddatas[0]->cons != NULL);

   lincons = consdata->lincons;

   /* check number of linear variables */
   SCIP_CALL( getLinearConsNVars(scip, lincons, consdata->linconstype, &nallvars) );
   assert(nallvars - nconsanddatas == consdata->nlinvars);
   nlinvars = consdata->nlinvars;

   if( nlinvars > MAXNVARS )
      return SCIP_OKAY;

   checkConsConsistency(scip, cons);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &allvars, nallvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &allcoefs, nallvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, MAXNVARS) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, MAXNVARS) );
   SCIP_CALL( SCIPallocBufferArray(scip, &repvars, MAXNVARS) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negated, MAXNVARS) );

   /* get variables and coefficients */
   SCIP_CALL( getLinearConsVarsData(scip, lincons, consdata->linconstype, allvars, allcoefs, &nallvars) );
   assert(nallvars > 0);

   /* calculate all not artificial linear variables */
   SCIP_CALL( getLinVarsAndAndRess(scip, cons, allvars, allcoefs, nallvars, linvars, lincoefs, &nlinvars,
         NULL, NULL, NULL, &nandress) );
   assert(nlinvars == consdata->nlinvars);
   assert(nandress == nallvars-nlinvars);

   constant = 0;

   /* get linear active representation */
   SCIP_CALL( SCIPgetProbvarLinearSum(scip, linvars, lincoefs, &nlinvars, MAXNVARS, &constant, &requiredsize, TRUE) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activelinvars, linvars, nlinvars) );

   if( requiredsize > MAXNVARS )
      goto TERMINATE;

   firstnlinvars = nlinvars;

   /* order the variables after index, to compare them easier */
   SCIPsortPtr((void**)linvars, SCIPvarCompActiveAndNegated, nlinvars);

   for( c = nconsanddatas - 1; c >= 0; --c )
   {
      consanddata = consanddatas[c];
      assert(consanddata != NULL);
      assert(consanddata->istransformed);

      /* choose correct variable array */
      if( consanddata->nnewvars > 0 )
      {
         vars = consanddata->newvars;
         nvars = consanddata->nnewvars;
      }
      else
      {
         vars = consanddata->vars;
         nvars = consanddata->nvars;
      }
      assert(nvars > 0 && vars != NULL);

      if( nvars > MAXNVARS )
         goto TERMINATE;

      BMSclearMemoryArray(negated, MAXNVARS);

      /* get linear active representation */
      SCIP_CALL( SCIPgetBinvarRepresentatives(scip, nvars, vars, repvars, negated) );
      SCIPsortPtr((void**)repvars, SCIPvarCompActiveAndNegated, nvars);

      oldnlinvars = nlinvars;

      /* determine all different variables over the linear variables and all variables in all and constraints */
      for( v = nvars - 1, v1 = nlinvars - 1; v >= 0 && v1 >= 0; )
      {
         SCIP_VAR* var;

         /* it appears that some fixed variables were not yet deleted */
         if( SCIPvarGetLbGlobal(repvars[v]) > 0.5 || SCIPvarGetUbGlobal(repvars[v]) < 0.5 )
            goto TERMINATE;

         assert(SCIPvarIsActive(linvars[v1]));
         assert(SCIPvarIsActive(repvars[v]) || (SCIPvarIsNegated(repvars[v]) && SCIPvarIsActive(SCIPvarGetNegationVar(repvars[v]))));

         if( SCIPvarIsActive(repvars[v]) )
            var = repvars[v];
         else
            var = SCIPvarGetNegationVar(repvars[v]);

         if( SCIPvarGetIndex(var) > SCIPvarGetIndex(linvars[v1]) )
         {
            if( nlinvars + 1 < MAXNVARS )
            {
               linvars[nlinvars] = var;
               ++nlinvars;
            }
            else
               goto TERMINATE;

            --v;
         }
         else if( SCIPvarGetIndex(var) < SCIPvarGetIndex(linvars[v1]) )
            --v1;
         else
         {
            --v;
            --v1;
         }
      }

      /* add the rest of variables */
      if( v >= 0 )
      {
         SCIP_VAR* var;

         for( ; v >= 0; --v )
         {
            /* it appears that some fixed variables were not yet deleted */
            if( SCIPvarGetLbGlobal(repvars[v]) > 0.5 || SCIPvarGetUbGlobal(repvars[v]) < 0.5 )
               goto TERMINATE;

            assert(SCIPvarIsActive(repvars[v]) || (SCIPvarIsNegated(repvars[v]) && SCIPvarIsActive(SCIPvarGetNegationVar(repvars[v]))));

            if( SCIPvarIsActive(repvars[v]) )
               var = repvars[v];
            else
               var = SCIPvarGetNegationVar(repvars[v]);

            if( nlinvars + 1 < MAXNVARS )
            {
               linvars[nlinvars] = var;
               ++nlinvars;
            }
            else
               goto TERMINATE;
         }
      }

      /* if some new variables were inserted we need to reorder the array */
      if( nlinvars > oldnlinvars )
      {
         /* order the variables after index, to compare them easier */
         SCIPsortPtr((void**)linvars, SCIPvarCompActiveAndNegated, nlinvars);
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &values, nlinvars) );
   xortype = -1;

   /* check values for variables which result in solutions which in the end lead to an XOR upgrade */
   for( v = (1 << nlinvars) - 1; v >= 0; --v ) /*lint !e701*/
   {
      int cnt = 0;
      for( v1 = nlinvars - 1; v1 >= 0; --v1 )
         if( v & (1 << v1) ) /*lint !e701*/
         {
            values[v1] = TRUE;
            ++cnt;
         }
         else
            values[v1] = FALSE;

      /* at maximum nlinvars values could be set to TRUE */
      assert(cnt <= nlinvars);

      SCIP_CALL( checkSolution(scip, linvars, nlinvars, values, activelinvars, lincoefs, firstnlinvars, constant,
            consdata->lhs, consanddatas, andcoefs, andnegs, nconsanddatas, cnt, &xortype) );
      if( xortype == -1 )
         break;
   }

   SCIPfreeBufferArray(scip, &values);

   assert(xortype >= -1 && xortype <= 1);

   if( xortype >= 0 )
   {
      (void) SCIPsnprintf(newname, SCIP_MAXSTRLEN, "%s_upgraded", SCIPconsGetName(lincons));

      SCIP_CALL( SCIPcreateConsXor(scip, &newcons, newname, (unsigned int) xortype, nlinvars, linvars,
            SCIPconsIsInitial(lincons), SCIPconsIsSeparated(lincons), SCIPconsIsEnforced(lincons), SCIPconsIsChecked(lincons),
            SCIPconsIsPropagated(lincons), SCIPconsIsLocal(lincons), SCIPconsIsModifiable(lincons),
            SCIPconsIsDynamic(lincons), SCIPconsIsRemovable(lincons), SCIPconsIsStickingAtNode(lincons)) );

      /* add and release new constraint */
      SCIP_CALL( SCIPaddCons(scip, newcons) );

      SCIPdebugMsg(scip, "created upgraded XOR constraint:\n");
      SCIPdebugMsg(scip, "old -> ");
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIPdebugMsg(scip, "new -> ");
      SCIPdebugPrintCons(scip, newcons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      ++(*naddconss);

      /* delete old constraints */
      SCIP_CALL( SCIPdelCons(scip, lincons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss) += 2;
   }

 TERMINATE:
   /* delete temporary memory */
   SCIPfreeBufferArray(scip, &activelinvars);
   SCIPfreeBufferArray(scip, &negated);
   SCIPfreeBufferArray(scip, &repvars);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &allcoefs);
   SCIPfreeBufferArray(scip, &allvars);

   return SCIP_OKAY;
}

/** try upgrading pseudoboolean logicor constraint to a linear constraint and/or remove possible and-constraints */
static
SCIP_RETCODE tryUpgradingLogicor(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONSHDLRDATA*const conshdlrdata,     /**< pseudoboolean constraint handler data */
   int*const             ndelconss,          /**< pointer to store number of deleted constraints */
   int*const             naddconss,          /**< pointer to count number of added constraints */
   int*const             nfixedvars,         /**< pointer to store number of fixed variables */
   int*const             nchgcoefs,          /**< pointer to store number of changed coefficients constraints */
   int*const             nchgsides,          /**< pointer to store number of changed sides constraints */
   SCIP_Bool*const       cutoff              /**< pointer to store if a cutoff happened */
   )
{
   CONSANDDATA** consanddatas;
   int nconsanddatas;
   SCIP_CONSDATA* consdata;
   int c;
   int v;
   int v2;
   SCIP_VAR** eqvars;
   int neqvars;
   int nminvars;
   int nmaxvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);
   assert(ndelconss != NULL);
   assert(nfixedvars != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(cutoff != NULL);
   assert(SCIPconsIsActive(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consanddatas = consdata->consanddatas;
   nconsanddatas = consdata->nconsanddatas;
   assert(nconsanddatas > 0 && consanddatas != NULL);

   assert(consdata->lincons != NULL);
   assert(consdata->linconstype == SCIP_LINEARCONSTYPE_LOGICOR);

   assert(consanddatas[0] != NULL);
   assert(consanddatas[0]->cons != NULL);

   if( nconsanddatas == 1 )
   {
      CONSANDDATA* consanddata;
      SCIP_VAR** allvars;
      SCIP_Real* allcoefs;
      int nallvars;
      SCIP_VAR** linvars;
      SCIP_Real* lincoefs;
      int nlinvars;
      SCIP_VAR** vars;
      int nvars;
      SCIP_CONS* lincons;
      SCIP_CONS* newcons;
      char newname[SCIP_MAXSTRLEN];
      SCIP_Real lhs;
      SCIP_Real rhs;

      /* if we have only one term left in the logicor constraint, the presolving should be done by the logicor
       * constraint handler
       */
      if( consdata->nlinvars == 0 )
      {
         return SCIP_OKAY;
      }

      /* for every old logicor constraint: sum_i (x_i) + res >= 1 , with an and-constraint of res as the resultant,
       *     which looks like 'res = y_1 * ... * y_n' => sum_i (n * x_i) + sum_j=1^n y_j >= n
       *
       * i.e. x_1 +  x_2 +  x_3 + x_4 * x_5 * x_6 >= 1
       *  => 3x_1 + 3x_2 + 3x_3 + x_4 + x_5 + x_6 >= 3
       */

      lincons = consdata->lincons;

      consanddata = consanddatas[0];
      assert(consanddata != NULL);
      assert(consanddata->istransformed);

      /* choose correct variable array to add locks for, we only add locks for now valid variables */
      if( consanddata->nnewvars > 0 )
      {
         vars = consanddata->newvars;
         nvars = consanddata->nnewvars;
      }
      else
      {
         vars = consanddata->vars;
         nvars = consanddata->nvars;
      }
      assert(nvars > 0 && vars != NULL);

      lhs = nvars;
      rhs = SCIPinfinity(scip);

      (void) SCIPsnprintf(newname, SCIP_MAXSTRLEN, "%s_upgraded", SCIPconsGetName(lincons));

      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, newname, 0, NULL, NULL, lhs, rhs,
            SCIPconsIsInitial(lincons), SCIPconsIsSeparated(lincons), SCIPconsIsEnforced(lincons), SCIPconsIsChecked(lincons),
            SCIPconsIsPropagated(lincons), SCIPconsIsLocal(lincons), SCIPconsIsModifiable(lincons),
            SCIPconsIsDynamic(lincons), SCIPconsIsRemovable(lincons), SCIPconsIsStickingAtNode(lincons)) );

      /* check number of linear variables */
      SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nallvars) );
      assert(nallvars == consdata->nlinvars + 1);

      nlinvars = consdata->nlinvars;

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &allvars, nallvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &allcoefs, nallvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nlinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nlinvars) );

      /* get variables and coefficients */
      SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, allvars, allcoefs, &nallvars) );
      assert(allcoefs != NULL);

      /* calculate all not artificial linear variables */
      SCIP_CALL( getLinVarsAndAndRess(scip, cons, allvars, allcoefs, nallvars, linvars, lincoefs, &nlinvars,
            NULL, NULL, NULL, NULL) );
      assert(nlinvars == consdata->nlinvars);

      /* add linear part to new constraint */
      for( v = 0; v < nlinvars; ++v )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, newcons, linvars[v], (SCIP_Real) nvars) );
      }

      /* add non-linear part to new constraint */
      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, newcons, vars[v], 1.0) );
      }

      /* add and release new constraint */
      SCIP_CALL( SCIPaddCons(scip, newcons) );

      SCIPdebugMsg(scip, "created upgraded linear constraint:\n");
      SCIPdebugMsg(scip, "old -> ");
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIPdebugMsg(scip, "new -> ");
      SCIPdebugPrintCons(scip, newcons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      ++(*naddconss);

      /* delete old constraints */
      SCIP_CALL( SCIPdelCons(scip, lincons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss) += 2;

      /* delete temporary memory */
      SCIPfreeBufferArray(scip, &lincoefs);
      SCIPfreeBufferArray(scip, &linvars);
      SCIPfreeBufferArray(scip, &allcoefs);
      SCIPfreeBufferArray(scip, &allvars);

      return SCIP_OKAY;
   }

   /* initializing array for variables which can appear in all consanddata objects */
   c = nconsanddatas - 1;
   assert(consanddatas[c]->istransformed);

   /* choose correct variable array */
   if( consanddatas[c]->nnewvars > 0 )
   {
      neqvars = consanddatas[c]->nnewvars;
      /* allocate temporary memory */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &eqvars, consanddatas[c]->newvars, neqvars) );
   }
   else
   {
      neqvars = consanddatas[c]->nvars;
      /* allocate temporary memory */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &eqvars, consanddatas[c]->vars, neqvars) );
   }
   nminvars = neqvars;
   nmaxvars = neqvars;
   assert(neqvars > 0 && eqvars != NULL);

#ifndef NDEBUG
   /* check that variables are sorted */
   for( v = neqvars - 1; v > 0; --v )
      assert(SCIPvarGetIndex(eqvars[v]) > SCIPvarGetIndex(eqvars[v - 1]));
#endif
   /* computing all variables which appear in all consanddata objects */
   for( --c ; c >= 0; --c )
   {
      CONSANDDATA* consanddata;
      SCIP_VAR** vars;
      int nvars;
      int nneweqvars;

      consanddata = consanddatas[c];
      assert(consanddata != NULL);
      assert(consanddatas[c]->istransformed);

      /* choose correct variable array to add locks for, we only add locks for now valid variables */
      if( consanddata->nnewvars > 0 )
      {
         vars = consanddata->newvars;
         nvars = consanddata->nnewvars;
      }
      else
      {
         vars = consanddata->vars;
         nvars = consanddata->nvars;
      }
      assert(nvars > 0 && vars != NULL);

#ifndef NDEBUG
      /* check that variables are sorted */
      for( v = nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(vars[v]) > SCIPvarGetIndex(vars[v - 1]));
#endif

      /* update minimal number of variables in and-constraint */
      if( nvars < nminvars )
         nminvars = nvars;
      /* update maximal number of variables in and-constraint */
      else if( nvars > nmaxvars )
         nmaxvars = nvars;
      assert(nminvars > 0);
      assert(nminvars <= nmaxvars);

      nneweqvars = 0;
      for( v = 0, v2 = 0; v < neqvars && v2 < nvars; )
      {
         int index1;
         int index2;

         assert(eqvars[v] != NULL);
         assert(vars[v2] != NULL);
         index1 = SCIPvarGetIndex(eqvars[v]);
         index2 = SCIPvarGetIndex(vars[v2]);

         /* check which variables are still in all and-constraints */
         if( index1 < index2 )
            ++v;
         else if( index1 > index2 )
            ++v2;
         else
         {
            assert(index1 == index2);
            assert(nneweqvars <= v);

            if( nneweqvars < v )
               eqvars[nneweqvars] = eqvars[v];
            ++nneweqvars;
            ++v;
            ++v2;
         }
      }
      neqvars = nneweqvars;

      /* now we only want to handle the easy case where nminvars == neqvars + 1
       * @todo: implement for the othercase too
       */
      if( nminvars > neqvars + 1 )
         break;

      /* if no variables overlap we have to stop */
      if( neqvars == 0 )
         break;
   }

   /* if all and-constraints in pseudoboolean constraint have some equal variables we can extract them and create a new
    * linear constraint; iff the number of equal variables is equal to the number of variables - 1 in all consanddata
    * objects then the new constraint will not contain any products; if no normal linear variables exist we can fix all
    * equal variables to 1
    *
    * e.g. x1 * x2 + x1 * x3 + x1 * x4 >= 1
    * =>   x1 = 1 /\ x2 + x3 + x4 >= 1
    *
    * e.g. x1 * x2 * x3 + x1 * x2 * x4 + x5 >= 1
    * =>  2x1 + 2x2 + x3 + x4 + 5x5 >= 5
    *
    * e.g. x1 * x2 * x3 + x1 * x4 >= 1
    * =>   x1 = 1 /\ x2 * x3 + x4 >= 1 (constraint is created indirectly, caused by the fixing of x1)
    *
    * @todo: implement the next cases
    *
    * e.g. x1 * x2 * x3 + x1 * x4 + x5 >= 1
    * =>  2x1 + x2 * x3 + x4 + 3x5 >= 3 (x2 * x3 will be a new and-constraint)
    *
    * e.g. x1 * x2 + x1 * x2 * x3 + x4 >= 1
    * =>  x1 + x2 + 2x4 >= 2
    *
    * e.g. x1 * x2 + x1 * x3 + x2 * x3 + sum_i x_i >= 1
    * =>  x1 + x2 + x3 + 2 * sum_i x_i >= 2
    *
    */

   /* Extract additional information ???
    *
    * e.g. x1 * x2 * x4 + x1 * x3 * x5 + x2 * x3 * x6 >= 1
    * =>  extract x1 + x2 + x3 >= 2
    */

   /* if we have no normal linear variable in the logicor constraint, we can fix all equal variables */
   if( neqvars > 0 && consdata->nlinvars == 0 )
   {
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      /* fix all equal variable in logicor constraints which have to be one to fulfill the constraint */
      for( v = 0; v < neqvars; ++v )
      {
         /* fix the variable which cannot be one */
         SCIP_CALL( SCIPfixVar(scip, eqvars[v], 1.0, &infeasible, &fixed) );
         if( infeasible )
         {
            SCIPdebugMsg(scip, " -> infeasible fixing\n");
            *cutoff = TRUE;
            goto TERMINATE;
         }
         if( fixed )
            ++(*nfixedvars);
      }

      /* if a complete consanddata object have all variables in common with all other consanddata objects, than we can
       * delete this constraint after fixing all equal variables
       */
      if( nminvars == neqvars )
      {
         /* delete old constraints */
         SCIP_CALL( SCIPdelCons(scip, consdata->lincons) );
         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss) += 2;

         goto TERMINATE;
      }
   }

   /* now the following condition grant us that we can linearize the whole constraint */
   if( neqvars > 0 && nminvars == nmaxvars && nminvars == neqvars + 1 )
   {
      SCIP_CONS* lincons;
      SCIP_CONS* newcons;
      char newname[SCIP_MAXSTRLEN];
      SCIP_Real lhs;
      SCIP_Real rhs;

      lhs = 1.0;
      rhs = SCIPinfinity(scip);

      lincons = consdata->lincons;

      (void) SCIPsnprintf(newname, SCIP_MAXSTRLEN, "%s_upgraded", SCIPconsGetName(lincons));

      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, newname, 0, NULL, NULL, lhs, rhs,
            SCIPconsIsInitial(lincons), SCIPconsIsSeparated(lincons), SCIPconsIsEnforced(lincons), SCIPconsIsChecked(lincons),
            SCIPconsIsPropagated(lincons), SCIPconsIsLocal(lincons), SCIPconsIsModifiable(lincons),
            SCIPconsIsDynamic(lincons), SCIPconsIsRemovable(lincons), SCIPconsIsStickingAtNode(lincons)) );

      /* if createcons == TRUE add all variables which are not in the eqvars array to the new constraint with
       * coefficient 1.0
       */
      for( c = nconsanddatas - 1; c >= 0; --c )
      {
         CONSANDDATA* consanddata;
         SCIP_VAR** vars;
         int nvars;

         consanddata = consanddatas[c];
         assert(consanddata != NULL);
         assert(consanddatas[c]->istransformed);

         /* choose correct variable array to add locks for, we only add locks for now valid variables */
         if( consanddata->nnewvars > 0 )
         {
            vars = consanddata->newvars;
            nvars = consanddata->nnewvars;
         }
         else
         {
            vars = consanddata->vars;
            nvars = consanddata->nvars;
         }
         assert(nvars > 0 && vars != NULL);

         for( v = 0, v2 = 0; v < neqvars && v2 < nvars; )
         {
            int index1;
            int index2;

            assert(eqvars[v] != NULL);
            assert(vars[v2] != NULL);
            index1 = SCIPvarGetIndex(eqvars[v]);
            index2 = SCIPvarGetIndex(vars[v2]);

            /* all variables in eqvars array must exist in all and-constraints */
            assert(index1 >= index2);

            if( index1 > index2 )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, newcons, vars[v2], 1.0) );
               ++v2;
            }
            else
            {
               assert(index1 == index2);
               ++v;
               ++v2;
            }
         }

         /* if we did not loop over all variables in the and-constraint, go on and fix variables */
         if( v2 < nvars )
         {
            assert(v == neqvars);
            for( ; v2 < nvars; ++v2)
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, newcons, vars[v2], 1.0) );
            }
         }
         assert(v == neqvars && v2 == nvars);
      }

      /* if we have normal linear variable in the logicor constraint, we did not fix all equal variables and we have to
       * add them with a coefficient of 'nconsanddatas'
       * we have to add also all normal linear variables with a coefficient of 'nconsanddatas * neqvars + 1'
       */
      if( consdata->nlinvars > 0 )
      {
         SCIP_VAR** vars;
         SCIP_Real* coefs;
         int nvars;
         SCIP_VAR** linvars;
         SCIP_Real* lincoefs;
         int nlinvars;

         /* add all equal variables */
         for( v = 0; v < neqvars; ++v )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, newcons, eqvars[v], (SCIP_Real)nconsanddatas) );
         }

         /* check number of linear variables */
         SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );
         assert(nvars == consdata->nlinvars + consdata->nconsanddatas);

         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nvars) );

         /* get variables and coefficients */
         SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );
         assert(nvars == 0 || (coefs != NULL));

#ifndef NDEBUG
         /* all coefficients have to be 1 */
         for( v = 0; v < nvars; ++v )
            assert(SCIPisEQ(scip, coefs[v], 1.0));
#endif
         /* calculate all not artificial linear variables */
         SCIP_CALL( getLinVarsAndAndRess(scip, cons, vars, coefs, nvars, linvars, lincoefs, &nlinvars,
               NULL, NULL, NULL, NULL) );
         assert(nlinvars == consdata->nlinvars);

         /* add all old normal linear variables */
         for( v = 0; v < nlinvars; ++v )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, newcons, linvars[v], (SCIP_Real)(nconsanddatas * neqvars + 1)) ); /*lint !e732 !e790*/
         }

         /* reset left hand side to correct value */
         SCIP_CALL( SCIPchgLhsLinear(scip, newcons, (SCIP_Real)(nconsanddatas * neqvars + 1)) ); /*lint !e732 !e790*/

         /* free temporary memory */
         SCIPfreeBufferArray(scip, &lincoefs);
         SCIPfreeBufferArray(scip, &linvars);
         SCIPfreeBufferArray(scip, &coefs);
         SCIPfreeBufferArray(scip, &vars);
      }

      /* add and release new constraint */
      SCIP_CALL( SCIPaddCons(scip, newcons) );

      SCIPdebugMsg(scip, "created upgraded linear constraint:\n");
      SCIPdebugMsg(scip, "old -> ");
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIPdebugMsg(scip, "new -> ");
      SCIPdebugPrintCons(scip, newcons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      ++(*naddconss);

      /* delete old constraints */
      SCIP_CALL( SCIPdelCons(scip, lincons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss) += 2;
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &eqvars);

   return SCIP_OKAY;
}

/** try upgrading pseudoboolean setppc constraint to a linear constraint and/or remove possible and-constraints */
static
SCIP_RETCODE tryUpgradingSetppc(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONSHDLRDATA*const conshdlrdata,     /**< pseudoboolean constraint handler data */
   int*const             ndelconss,          /**< pointer to store number of deleted constraints */
   int*const             naddconss,          /**< pointer to count number of added constraints */
   int*const             nfixedvars,         /**< pointer to store number of fixed variables */
   int*const             nchgcoefs,          /**< pointer to store number of changed coefficients constraints */
   int*const             nchgsides,          /**< pointer to store number of changed sides constraints */
   SCIP_Bool*const       cutoff              /**< pointer to store if a cutoff happened */
   )
{
   CONSANDDATA** consanddatas;
   int nconsanddatas;
   SCIP_CONSDATA* consdata;
   SCIP_SETPPCTYPE type;
   int c;
   int v;
   int v2;
   SCIP_VAR** eqvars;
   int neqvars;
   int nminvars;
   int nmaxvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);
   assert(ndelconss != NULL);
   assert(nfixedvars != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(cutoff != NULL);
   assert(SCIPconsIsActive(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consanddatas = consdata->consanddatas;
   nconsanddatas = consdata->nconsanddatas;
   assert(nconsanddatas > 0 && consanddatas != NULL);

   assert(consdata->lincons != NULL);
   assert(consdata->linconstype == SCIP_LINEARCONSTYPE_SETPPC);

   type = SCIPgetTypeSetppc(scip, consdata->lincons);

   switch( type )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
   case SCIP_SETPPCTYPE_PACKING:
      break;
   case SCIP_SETPPCTYPE_COVERING:
      return SCIP_OKAY;
   default:
      SCIPerrorMessage("unknown setppc type\n");
      return SCIP_INVALIDDATA;
   }

   assert(consanddatas[0] != NULL);
   assert(consanddatas[0]->cons != NULL);

   if( nconsanddatas == 1 )
   {
      /* if we have only one term left in the setppc constraint, the presolving should be done by the setppc constraint handler */
      if( consdata->nlinvars == 0 )
      {
         return SCIP_OKAY;
      }

      /* @todo: implement the following */

      /* for each set packing constraint:
       *     sum_i (x_i) + res <= 1 , with and-constraint of res as the resultant like res = y_1 * ... * y_n
       *  => sum_i (n * x_i) + sum_j=1^n y_j <= n + n-1
       *
       * i.e. x_1 + x_2 + x_3 + x_4*x_5*x_6 <= 1
       *  => 3x_1 + 3x_2 + 3x_3 + x_4 + x_5 + x_6 <= 5
       */

      /* for each set partitioning constraint:
       *     sum_i (x_i) + res = 1 , with the corresponding and-constraint of res like
       *                             res = y_1 * ... * y_n
       *
       *  => n <= sum_i (n * x_i) + sum_j=1^n y_j <= 2 * n - 1
       *
       * i.e. x_1 + x_2 + x_3 + x_4*x_5*x_6 = 1
       *  => 3 <= 3x_1 + 3x_2 + 3x_3 + x_4 + x_5 + x_6 <= 5
       *
       */

      return SCIP_OKAY;
   }

   if( consdata->nlinvars > 0 )
   {
      /* @todo: */
      return SCIP_OKAY;
   }
   assert(consdata->nlinvars == 0 && nconsanddatas > 1);

   c = nconsanddatas - 1;
   assert(consanddatas[c]->istransformed);

   /* initializing array for variables which can appear in all consanddata objects */
   if( consanddatas[c]->nnewvars > 0 )
   {
      neqvars = consanddatas[c]->nnewvars;
      /* allocate temporary memory */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &eqvars, consanddatas[c]->newvars, neqvars) );
   }
   else
   {
      neqvars = consanddatas[c]->nvars;
      /* allocate temporary memory */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &eqvars, consanddatas[c]->vars, neqvars) );
   }
   nminvars = neqvars;
   nmaxvars = neqvars;
   assert(neqvars > 0 && eqvars != NULL);

#ifndef NDEBUG
   /* check that variables are sorted */
   for( v = neqvars - 1; v > 0; --v )
      assert(SCIPvarGetIndex(eqvars[v]) > SCIPvarGetIndex(eqvars[v - 1]));
#endif

   for( --c ; c >= 0; --c )
   {
      CONSANDDATA* consanddata;
      SCIP_VAR** vars;
      int nvars;
      int nneweqvars;

      consanddata = consanddatas[c];
      assert(consanddata != NULL);
      assert(consanddatas[c]->istransformed);

      /* choose correct variable array to add locks for, we only add locks for now valid variables */
      if( consanddata->nnewvars > 0 )
      {
         vars = consanddata->newvars;
         nvars = consanddata->nnewvars;
      }
      else
      {
         vars = consanddata->vars;
         nvars = consanddata->nvars;
      }
      assert(nvars > 0 && vars != NULL);

#ifndef NDEBUG
      /* check that variables are sorted */
      for( v = nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(vars[v]) > SCIPvarGetIndex(vars[v - 1]));
#endif

      /* update minimal number of variables in and-constraint */
      if( nvars < nminvars )
         nminvars = nvars;
      /* update maximal number of variables in and-constraint */
      else if( nvars > nmaxvars )
         nmaxvars = nvars;
      assert(nminvars > 0);
      assert(nminvars <= nmaxvars);

      nneweqvars = 0;
      for( v = 0, v2 = 0; v < neqvars && v2 < nvars; )
      {
         int index1;
         int index2;

         assert(eqvars[v] != NULL);
         assert(vars[v2] != NULL);
         index1 = SCIPvarGetIndex(eqvars[v]);
         index2 = SCIPvarGetIndex(vars[v2]);

         /* check which variables are still in all and-constraints */
         if( index1 < index2 )
            ++v;
         else if( index1 > index2 )
            ++v2;
         else
         {
            assert(index1 == index2);
            assert(nneweqvars <= v);

            if( nneweqvars < v )
               eqvars[nneweqvars] = eqvars[v];
            ++nneweqvars;
            ++v;
            ++v2;
         }
      }
      neqvars = nneweqvars;

      /* now we only want to handle the easy case where nminvars == neqvars + 1
       * @todo: implement for the othercase too
       */
      if( nminvars > neqvars + 1 && type != SCIP_SETPPCTYPE_PARTITIONING)
         break;

      if( neqvars == 0 )
         break;
   }

   /* if all and-constraints in pseudoboolean constraint have the same length and some equal variables we can upgrade
    * the linear constraint and fix some variables in setpartitioning case
    *
    * e.g. x1 * x2 + x1 * x3 + x1 * x4 <= 1
    * =>  3x1 + x2 + x3 + x4 <= 4
    *
    * e.g. x1 * x2 * x3 + x1 * x2 * x4 <= 1
    * =>  2x1 + 2x2 + x3 + x4 <= 5
    *
    * e.g. x1 * x2 + x1 * x2 * x3 + x1 * x2 * x4 <= 1
    * =>  3x1 + 3x2 + x3 + x4 <= 6
    *
    * e.g. x1 * x2 + x1 * x3 == 1
    * =>   x1 = 1 /\ x2 + x3 == 1
    *
    * e.g. x1 * x2 * x3 + x1 * x4 == 1
    * =>   x1 = 1 /\ x2 * x3 + x4 == 1 (constraint is created indirectly, caused by the fixing of x1)
    *
    * e.g. x1 * x2 + x1 * x2 * x3 + x1 * x2 * x4 == 1
    * =>   x1 = 1, x2 = 1, x3 = 0, x4 = 0
    *
    * e.g. x1 * x2 + x1 * x2 * x3 + x1 * x2 * x4 * x5 == 1
    * =>   x1 = 1, x2 = 1, x3 = 0 /\ x4 * x5 == 0
    *
    * @todo: implement the next cases
    *
    * e.g. x1 * x2 * x3 + x1 * x2 * x4 + x5 <= 1
    * =>  2x1 + 2x2 + x3 + x4 + x5 <= 5
    *
    */
   if( neqvars > 0 && ((nminvars == nmaxvars && nminvars == neqvars + 1) || (nminvars == neqvars) || (type == SCIP_SETPPCTYPE_PARTITIONING)) )
   {
      SCIP_CONS* lincons;
      SCIP_CONS* newcons;
      char newname[SCIP_MAXSTRLEN];
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      SCIP_Bool createcons;
      SCIP_Bool deletecons;

      newcons = NULL;

      /* determine new sides of linear constraint */
      if( type == SCIP_SETPPCTYPE_PARTITIONING )
      {
         lhs = 1.0;
         rhs = 1.0;
      }
      else
      {
         assert(type == SCIP_SETPPCTYPE_PACKING);
         lhs = -SCIPinfinity(scip);
         rhs = 1.0;
      }

      /* if one and-constraint was completely contained in all other and-constraints, we have to reduced the right hand
       * side by 1
       */
      if( neqvars == nminvars )
         rhs -= 1.0;

      createcons = (SCIPisLE(scip, lhs, rhs) && ((nminvars == nmaxvars && nminvars == neqvars + 1) || (nminvars == neqvars)));
      assert(createcons || type == SCIP_SETPPCTYPE_PARTITIONING);

      deletecons = (type == SCIP_SETPPCTYPE_PARTITIONING && nminvars == neqvars);

      lincons = consdata->lincons;

      if( createcons )
      {
         (void) SCIPsnprintf(newname, SCIP_MAXSTRLEN, "%s_upgraded", SCIPconsGetName(lincons));

         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, newname, 0, NULL, NULL, lhs, rhs,
               SCIPconsIsInitial(lincons), SCIPconsIsSeparated(lincons), SCIPconsIsEnforced(lincons), SCIPconsIsChecked(lincons),
               SCIPconsIsPropagated(lincons), SCIPconsIsLocal(lincons), SCIPconsIsModifiable(lincons),
               SCIPconsIsDynamic(lincons), SCIPconsIsRemovable(lincons), SCIPconsIsStickingAtNode(lincons)) );
      }

      /* if createcons == TRUE add all variables which are not in the eqvars array to the new constraint with
       * coefficient 1.0
       *
       * otherwise (if createcons == FALSE) fix all variables to zero which are not in the eqvars array and if we have a
       * set partitioning constraint
       */
      for( c = nconsanddatas - 1; c >= 0; --c )
      {
         CONSANDDATA* consanddata;
         SCIP_VAR** vars;
         int nvars;

         consanddata = consanddatas[c];
         assert(consanddata != NULL);
         assert(consanddatas[c]->istransformed);

         /* choose correct variable array to add locks for, we only add locks for now valid variables */
         if( consanddata->nnewvars > 0 )
         {
            vars = consanddata->newvars;
            nvars = consanddata->nnewvars;
         }
         else
         {
            vars = consanddata->vars;
            nvars = consanddata->nvars;
         }
         assert(nvars > 0 && vars != NULL);

         /* if the consanddata object has at least two more different variables then the equal variables we have to fix the resultant to zero */
         if( deletecons && neqvars + 1 < nvars )
         {
            assert(SCIPgetResultantAnd(scip, consanddata->cons) != NULL);

            /* fix the resultant variable which have to be zero */
            SCIP_CALL( SCIPfixVar(scip, SCIPgetResultantAnd(scip, consanddata->cons), 0.0, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, " -> infeasible fixing\n");
               *cutoff = TRUE;
               goto TERMINATE;
            }
            if( fixed )
               ++(*nfixedvars);

            continue;
         }

         /* if the consanddata object has at exactly one more different variable then the equal variables we have to fix it to zero */
         for( v = 0, v2 = 0; v < neqvars && v2 < nvars; )
         {
            int index1;
            int index2;

            assert(eqvars[v] != NULL);
            assert(vars[v2] != NULL);
            index1 = SCIPvarGetIndex(eqvars[v]);
            index2 = SCIPvarGetIndex(vars[v2]);

            /* all variables in eqvars array must exist in all and-constraints */
            assert(index1 >= index2);

            if( index1 > index2 )
            {
               if( createcons )
               {
                  assert(newcons != NULL);
                  SCIP_CALL( SCIPaddCoefLinear(scip, newcons, vars[v2], 1.0) );
               }
               else if( deletecons )
               {
                  /* fix the variable which cannot be one */
                  SCIP_CALL( SCIPfixVar(scip, vars[v2], 0.0, &infeasible, &fixed) );
                  if( infeasible )
                  {
                     SCIPdebugMsg(scip, " -> infeasible fixing\n");
                     *cutoff = TRUE;
                     goto TERMINATE;
                  }
                  if( fixed )
                     ++(*nfixedvars);
               }
               ++v2;
            }
            else
            {
               assert(index1 == index2);

               ++v;
               ++v2;
            }
         }

         /* if we did not loop over all variables in the and-constraint, go on and fix variables */
         if( v2 < nvars )
         {
            assert(v == neqvars);
            for( ; v2 < nvars; ++v2)
            {
               if( createcons )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, newcons, vars[v2], 1.0) );
               }
               else if( deletecons )
               {
                  /* fix the variable which cannot be one */
                  SCIP_CALL( SCIPfixVar(scip, vars[v2], 0.0, &infeasible, &fixed) );
                  if( infeasible )
                  {
                     SCIPdebugMsg(scip, " -> infeasible fixing\n");
                     *cutoff = TRUE;
                     goto TERMINATE;
                  }
                  if( fixed )
                     ++(*nfixedvars);
               }
            }
         }
         assert(v == neqvars && v2 == nvars);
      }

      /* fix all equal variable in set-partitioning constraints which have to be one, in set-packing constraint we have
       * to add these variable with a coeffcient as big as (nconsanddatas - 1)
       */
      for( v = 0; v < neqvars; ++v )
      {
         if( type == SCIP_SETPPCTYPE_PARTITIONING )
         {
            /* fix the variable which have to be one */
            SCIP_CALL( SCIPfixVar(scip, eqvars[v], 1.0, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, " -> infeasible fixing\n");
               *cutoff = TRUE;
               goto TERMINATE;
            }
            if( fixed )
               ++(*nfixedvars);
         }
         else
         {
            assert(type == SCIP_SETPPCTYPE_PACKING);
            SCIP_CALL( SCIPaddCoefLinear(scip, newcons, eqvars[v], (SCIP_Real)(nconsanddatas - 1)) );
         }
      }

      /* correct right hand side for set packing constraint */
      if( type == SCIP_SETPPCTYPE_PACKING )
      {
         assert(createcons);
         assert(newcons != NULL);

         SCIP_CALL( SCIPchgRhsLinear(scip, newcons, rhs + (SCIP_Real)((nconsanddatas - 1) * neqvars)) ); /*lint !e790*/
      }

      /* add and release new constraint */
      if( createcons )
      {
         SCIP_CALL( SCIPaddCons(scip, newcons) );

         SCIPdebugMsg(scip, "created upgraded linear constraint:\n");
         SCIPdebugMsg(scip, "old -> ");
         SCIPdebugPrintCons(scip, lincons, NULL);
         SCIPdebugMsg(scip, "new -> ");
         SCIPdebugPrintCons(scip, newcons, NULL);

         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         ++(*naddconss);

         assert(!deletecons);
         deletecons = TRUE;
      }

      if( deletecons )
      {
         /* delete old constraints */
         SCIP_CALL( SCIPdelCons(scip, lincons) );
         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss) += 2;
      }
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &eqvars);

   return SCIP_OKAY;
}

/** try upgrading pseudoboolean constraint to a linear constraint and/or remove possible and-constraints */
static
SCIP_RETCODE tryUpgrading(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONSHDLRDATA*const conshdlrdata,     /**< pseudoboolean constraint handler data */
   int*const             ndelconss,          /**< pointer to store number of upgraded constraints */
   int*const             naddconss,          /**< pointer to count number of added constraints */
   int*const             nfixedvars,         /**< pointer to store number of fixed variables */
   int*const             nchgcoefs,          /**< pointer to store number of changed coefficients constraints */
   int*const             nchgsides,          /**< pointer to store number of changed sides constraints */
   SCIP_Bool*const       cutoff              /**< pointer to store if a cutoff happened */
   )
{
#ifndef NDEBUG
   CONSANDDATA** consanddatas;
#endif
   SCIP_CONSDATA* consdata;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);
   assert(ndelconss != NULL);
   assert(nfixedvars != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(cutoff != NULL);
   assert(SCIPconsIsActive(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lincons != NULL);

#ifndef NDEBUG
   consanddatas = consdata->consanddatas;
   assert(consdata->nconsanddatas == 0 || consanddatas != NULL);
#endif

   /* if no consanddata-objects in pseudoboolean constraint are left, create the corresponding linear constraint */
   if( consdata->nconsanddatas == 0 )
   {
      SCIPconsAddUpgradeLocks(consdata->lincons, -1);
      assert(SCIPconsGetNUpgradeLocks(consdata->lincons) == 0);

      /* @TODO: maybe it is better to create everytime a standard linear constraint instead of letting the special
       *        linear constraint stay
       */
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);

      return SCIP_OKAY;
   }

   /* check number of linear variables */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );
   assert(consdata->nlinvars + consdata->nconsanddatas == nvars);

   switch( consdata->linconstype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      SCIP_CALL( tryUpgradingXor(scip, cons, conshdlrdata, ndelconss, naddconss, nfixedvars, nchgcoefs, nchgsides, cutoff) );
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
      SCIP_CALL( tryUpgradingLogicor(scip, cons, conshdlrdata, ndelconss, naddconss, nfixedvars, nchgcoefs, nchgsides, cutoff) );
      break;
   case SCIP_LINEARCONSTYPE_KNAPSACK:
      break;
   case SCIP_LINEARCONSTYPE_SETPPC:
      SCIP_CALL( tryUpgradingSetppc(scip, cons, conshdlrdata, ndelconss, naddconss, nfixedvars, nchgcoefs, nchgsides, cutoff) );
      if( !SCIPconsIsDeleted(cons) )
      {
	 SCIP_CALL( tryUpgradingXor(scip, cons, conshdlrdata, ndelconss, naddconss, nfixedvars, nchgcoefs, nchgsides, cutoff) );
      }
      break;
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
      SCIP_CALL( tryUpgradingXor(scip, cons, conshdlrdata, ndelconss, naddconss, nfixedvars, nchgcoefs, nchgsides, cutoff) );
#endif
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPconsIsDeleted(cons) )
   {
      /* update the uses counter of consandata objects which are used in pseudoboolean constraint, which was deleted and
       * probably delete and-constraints
       */
      SCIP_CALL( updateConsanddataUses(scip, cons, conshdlrdata, ndelconss) );
   }

   consdata->upgradetried = TRUE;

   return SCIP_OKAY;
}

/** check if we can aggregated some variables */
static
SCIP_RETCODE findAggregation(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONSHDLRDATA*const conshdlrdata,     /**< pseudoboolean constraint handler data */
   int*const             ndelconss,          /**< pointer to store number of upgraded constraints */
   int*const             naggrvars,          /**< pointer to store number of aggregated variables */
   SCIP_Bool*const       cutoff              /**< pointer to store if a cutoff happened */
   )
{
   CONSANDDATA** consanddatas;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** allvars;
   int* varcount[2];
   SCIP_VAR** repvars;
   SCIP_Bool* negated;
   SCIP_VAR** vars;
   int nconsanddatas;
   int nvars;
   int zerocount;
   int onecount;
   int twocount;
   int othercount;
   int c;
   int v;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);
   assert(ndelconss != NULL);
   assert(naggrvars != NULL);
   assert(cutoff != NULL);
   assert(SCIPconsIsActive(cons));

   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lincons != NULL);

   consanddatas = consdata->consanddatas;
   nconsanddatas = consdata->nconsanddatas;
   assert(nconsanddatas == 0 || consanddatas != NULL);

   /* we have only one special case for aggregations, a set-partinioning constraint */
   if( consdata->linconstype != SCIP_LINEARCONSTYPE_SETPPC || SCIPgetTypeSetppc(scip, consdata->lincons) != SCIP_SETPPCTYPE_PARTITIONING )
      return SCIP_OKAY;

   assert(SCIPisEQ(scip, consdata->rhs, consdata->lhs));
   assert(SCIPisEQ(scip, consdata->rhs, 1.0));

   if( nconsanddatas < 2 || nconsanddatas > 3 )
      return SCIP_OKAY;

#ifndef NDEBUG
   /* check number of linear variables */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );
   assert(consdata->nlinvars + nconsanddatas == nvars);
#endif

   if( consdata->nlinvars != 1 )
      return SCIP_OKAY;

   /* check valid number of variables */
   if( consanddatas[0]->nnewvars > 0 )
      nvars = consanddatas[0]->nnewvars;
   else
      nvars = consanddatas[0]->nvars;

   if( consanddatas[1]->nnewvars > 0 )
   {
      if( nvars != consanddatas[1]->nnewvars )
	 return SCIP_OKAY;
   }
   else if( nvars != consanddatas[1]->nvars )
      return SCIP_OKAY;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &allvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(varcount[0]), nvars) );
   BMSclearMemoryArray(varcount[0], nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &(varcount[1]), nvars) );
   BMSclearMemoryArray(varcount[1], nvars);

   SCIP_CALL( SCIPallocBufferArray(scip, &repvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negated, nvars) );
   BMSclearMemoryArray(negated, nvars);

   /* get valid variables */
   if( consanddatas[nconsanddatas - 1]->nnewvars > 0 )
      vars = consanddatas[nconsanddatas - 1]->newvars;
   else
      vars = consanddatas[nconsanddatas - 1]->vars;

   /* get linear active representation */
   SCIP_CALL( SCIPgetBinvarRepresentatives(scip, nvars, vars, repvars, negated) );
   SCIPsortPtrBool((void**)repvars, negated, SCIPvarCompActiveAndNegated, nvars);

#ifndef NDEBUG
   /* and-constraints have to be merged in order to check for aggregation */
   for( v = 1; v < nvars; ++v )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;

      /* it appears that some fixed variables were not yet deleted */
      if( SCIPvarGetLbGlobal(repvars[v-1]) > 0.5 || SCIPvarGetUbGlobal(repvars[v-1]) < 0.5 )
	 goto TERMINATE;

      /* it appears that some fixed variables were not yet deleted */
      if( SCIPvarGetLbGlobal(repvars[v]) > 0.5 || SCIPvarGetUbGlobal(repvars[v]) < 0.5 )
	 goto TERMINATE;

      assert(SCIPvarIsActive(repvars[v]) || (SCIPvarIsNegated(repvars[v]) && SCIPvarIsActive(SCIPvarGetNegationVar(repvars[v]))));
      assert(SCIPvarIsActive(repvars[v-1]) || (SCIPvarIsNegated(repvars[v-1]) && SCIPvarIsActive(SCIPvarGetNegationVar(repvars[v-1]))));
      assert(SCIPvarIsActive(repvars[v]) != negated[v]);
      assert(SCIPvarIsActive(repvars[v-1]) != negated[v-1]);

      var1 = (negated[v-1] ? SCIPvarGetNegationVar(repvars[v-1]) : repvars[v-1]);
      var2 = (negated[v] ? SCIPvarGetNegationVar(repvars[v]) : repvars[v]);
      assert(var1 != var2);
   }
#endif

   /* initializing the statuses of all appearing variables */
   for( v = nvars - 1; v >= 0; --v )
   {
      /* it appears that some fixed variables were not yet deleted */
      if( SCIPvarGetLbGlobal(repvars[v]) > 0.5 || SCIPvarGetUbGlobal(repvars[v]) < 0.5 )
	 goto TERMINATE;

      assert(SCIPvarIsActive(repvars[v]) || (SCIPvarIsNegated(repvars[v]) && SCIPvarIsActive(SCIPvarGetNegationVar(repvars[v]))));
      assert(SCIPvarIsActive(repvars[v]) != negated[v]);

      allvars[v] = negated[v] ? SCIPvarGetNegationVar(repvars[v]) : repvars[v];

      ++(varcount[negated[v]][v]);
   }

   for( c = nconsanddatas - 2; c >= 0; --c )
   {
      int pos = -1;

      /* get valid variables */
      if( consanddatas[nconsanddatas - 1]->nnewvars > 0 )
	 vars = consanddatas[c]->newvars;
      else
	 vars = consanddatas[c]->vars;

      /* need to reset the negated flags */
      BMSclearMemoryArray(negated, nvars);

      /* get linear active representation */
      SCIP_CALL( SCIPgetBinvarRepresentatives(scip, nvars, vars, repvars, negated) );
      SCIPsortPtrBool((void**)repvars, negated, SCIPvarCompActiveAndNegated, nvars);

#ifndef NDEBUG
      /* and-constraints have to be merged in order to check for aggregation */
      for( v = 1; v < nvars; ++v )
      {
	 SCIP_VAR* var1;
	 SCIP_VAR* var2;

	 /* it appears that some fixed variables were not yet deleted */
	 if( SCIPvarGetLbGlobal(repvars[v-1]) > 0.5 || SCIPvarGetUbGlobal(repvars[v-1]) < 0.5 )
	    goto TERMINATE;

	 /* it appears that some fixed variables were not yet deleted */
	 if( SCIPvarGetLbGlobal(repvars[v]) > 0.5 || SCIPvarGetUbGlobal(repvars[v]) < 0.5 )
	    goto TERMINATE;

	 assert(SCIPvarIsActive(repvars[v]) || (SCIPvarIsNegated(repvars[v]) && SCIPvarIsActive(SCIPvarGetNegationVar(repvars[v]))));
	 assert(SCIPvarIsActive(repvars[v-1]) || (SCIPvarIsNegated(repvars[v-1]) && SCIPvarIsActive(SCIPvarGetNegationVar(repvars[v-1]))));
	 assert(SCIPvarIsActive(repvars[v]) != negated[v]);
	 assert(SCIPvarIsActive(repvars[v-1]) != negated[v-1]);

	 var1 = (negated[v-1] ? SCIPvarGetNegationVar(repvars[v-1]) : repvars[v-1]);
	 var2 = (negated[v] ? SCIPvarGetNegationVar(repvars[v]) : repvars[v]);
	 assert(var1 != var2);
      }
#endif

      /* update the statuses of all appearing variables */
      for( v = nvars - 1; v >= 0; --v )
      {
	 /* it appears that some fixed variables were not yet deleted */
	 if( SCIPvarGetLbGlobal(repvars[v]) > 0.5 || SCIPvarGetUbGlobal(repvars[v]) < 0.5 )
	    goto TERMINATE;

	 assert(SCIPvarIsActive(repvars[v]) || (SCIPvarIsNegated(repvars[v]) && SCIPvarIsActive(SCIPvarGetNegationVar(repvars[v]))));
	 assert(SCIPvarIsActive(repvars[v]) != negated[v]);

	 /* we can only find an aggregation if all and constraints have the same variables */
	 if( SCIPsortedvecFindPtr((void**)allvars, SCIPvarCompActiveAndNegated, repvars[v], nvars, &pos) )
	 {
	    assert(pos >= 0 && pos < nvars);

	    ++(varcount[negated[v]][pos]);
	 }
	 else
	    goto TERMINATE;
      }
   }

   zerocount = 0;
   onecount = 0;
   twocount = 0;
   othercount = 0;

   /* count number of multiple appearances of a variable */
   for( i = 1; i >= 0; --i )
   {
      for( v = nvars - 1; v >= 0; --v )
      {
	 assert(SCIPvarIsActive(allvars[v]));

	 if( varcount[i][v] == 0 )
	    ++zerocount;
	 else if( varcount[i][v] == 1 )
	    ++onecount;
	 else if( varcount[i][v] == 2 )
	    ++twocount;
	 else
	    ++othercount;
      }
   }

   /* exactly one variable in all and-constraints appears as active and as negated variable */
   if( othercount == 0 )
   {
      /* we have a constraint in the form of: x1 + x2 * x3 * ... * x_n + ~x2 * x3 * ... * x_n == 1
       * this leads to the aggregation x1 = 1 - x3 * ... * x_n
       */
      if( nconsanddatas == 2 && twocount == nvars - 1 && onecount == 2 && zerocount == 1 )
      {
	 SCIP_VAR** consvars;
	 SCIP_Real* conscoefs;
	 int nconsvars;
	 SCIP_VAR* linvar;
	 SCIP_Real lincoef;
	 int nlinvars;

	 /* allocate temporary memory */
	 SCIP_CALL( SCIPallocBufferArray(scip, &consvars, consdata->nlinvars + nconsanddatas) );
	 SCIP_CALL( SCIPallocBufferArray(scip, &conscoefs, consdata->nlinvars + nconsanddatas) );

	 /* get variables and coefficients */
	 SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, consvars, conscoefs, &nconsvars) );
	 assert(nconsvars == consdata->nlinvars + nconsanddatas);
	 assert(conscoefs != NULL);

#ifndef NDEBUG
	 /* all coefficients have to be 1 */
	 for( v = 0; v < nconsvars; ++v )
	    assert(SCIPisEQ(scip, conscoefs[v], 1.0));
#endif
	 linvar = NULL;

	 /* calculate all not artificial linear variables */
	 SCIP_CALL( getLinVarsAndAndRess(scip, cons, consvars, conscoefs, nconsvars, &linvar, &lincoef, &nlinvars,
               NULL, NULL, NULL, NULL) );
	 assert(nlinvars == 1);
	 assert(linvar != NULL);

	 SCIPfreeBufferArray(scip, &conscoefs);
	 SCIPfreeBufferArray(scip, &consvars);

	 /* if all and-constraints have exactly two variables */
	 if( nvars == 2 )
	 {
	    SCIP_VAR* var;
	    SCIP_Bool breaked;
	    SCIP_Bool redundant;
	    SCIP_Bool infeasible;
	    SCIP_Bool aggregated;

	    var = NULL;
	    breaked = FALSE;

	    /* find necessary variables, which only occur once */
	    for( i = 1; i >= 0; --i )
	    {
	       for( v = nvars - 1; v >= 0; --v )
	       {
		  assert(i == 1 || SCIPvarGetNegatedVar(allvars[v]) != NULL);
		  if( varcount[i][v] == 2 )
		  {
		     var = i ? SCIPvarGetNegatedVar(allvars[v]) : allvars[v];

		     breaked = TRUE;
		     break;
		  }
	       }

	       if( breaked )
		  break;
	    }
	    assert(var != NULL);

	    SCIPdebugMsg(scip, "aggregating variables <%s> == 1 - <%s> in pseudoboolean <%s>\n", SCIPvarGetName(linvar), SCIPvarGetName(var), SCIPconsGetName(cons));

	    SCIP_CALL( SCIPaggregateVars(scip, linvar, var, 1.0, 1.0, 1.0, &infeasible, &redundant, &aggregated) );

	    SCIPdebugPrintCons(scip, cons, NULL);
	    SCIPdebugMsg(scip, "aggregation of variables: <%s> == 1 - <%s>, infeasible = %u, aggregated = %u\n", SCIPvarGetName(linvar), SCIPvarGetName(var), infeasible, aggregated);

	    if( infeasible )
	       *cutoff = TRUE;
	    else
	    {
	       if( aggregated )
		  ++(*naggrvars);

	       /* delete old constraints */
	       SCIP_CALL( SCIPdelCons(scip, consdata->lincons) );
	       SCIP_CALL( SCIPdelCons(scip, cons) );
	       (*ndelconss) += 2;
	    }
	 }
#if 0
	 else
	 {
	    /* @todo */
	    /* delete allvars[samepos] from all and-constraints which appear in this pseudoboolean constraint, and delete
	     * all but one of the remaining and-constraint
	     *
	     * it is the same like aggregating linvar with the resultant of the product, which is the same in all and-
	     * constraints without allvars[samepos]
	     *
	     * e.g. x1 + x2*x_3*...x_n + ~x2*x_3*...x_n = 1 => x1 = 1 - x_3*...x_n
	     */
	 }
#endif
      }
      /* we have a constraint in the form of: x1 + x2 * x3 + ~x2 * x3 + ~x2 * ~x3 == 1
       * this leads to the aggregation x1 = x2 * ~x3
       *
       * @todo: implement more general step, that one combination of the variables in the and constraints is missing in
       *        the pseudoboolean constraint, which leads to the same result, that the only linear variable is the
       *        resultant of the missing and-constraint
       */
      else if( nvars == 2 && nconsanddatas == 3 && twocount == 2 && onecount == 2 && zerocount == 0)
      {
	 SCIP_VAR** consvars;
	 SCIP_Real* conscoefs;
	 int nconsvars;
	 SCIP_VAR* linvar;
	 SCIP_Real lincoef;
	 int nlinvars;
	 SCIP_VAR* newandvars[2];
	 SCIP_Bool breaked;
	 SCIP_CONS* newcons;
	 char name[SCIP_MAXSTRLEN];

	 /* allocate temporary memory */
	 SCIP_CALL( SCIPallocBufferArray(scip, &consvars, consdata->nlinvars + nconsanddatas) );
	 SCIP_CALL( SCIPallocBufferArray(scip, &conscoefs, consdata->nlinvars + nconsanddatas) );

	 /* get variables and coefficients */
	 SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, consvars, conscoefs, &nconsvars) );
	 assert(nconsvars == consdata->nlinvars + nconsanddatas);
	 assert(conscoefs != NULL);

#ifndef NDEBUG
	 /* all coefficients have to be 1 */
	 for( v = 0; v < nconsvars; ++v )
	    assert(SCIPisEQ(scip, conscoefs[v], 1.0));
#endif
	 linvar = NULL;

	 /* calculate all not artificial linear variables */
	 SCIP_CALL( getLinVarsAndAndRess(scip, cons, consvars, conscoefs, nconsvars, &linvar, &lincoef, &nlinvars,
               NULL, NULL, NULL, NULL) );
	 assert(nlinvars == 1);
	 assert(linvar != NULL);

	 SCIPfreeBufferArray(scip, &conscoefs);
	 SCIPfreeBufferArray(scip, &consvars);

	 newandvars[0] = NULL;
	 newandvars[1] = NULL;
	 breaked = FALSE;

	 /* find necessary variables, which only occur once */
	 for( i = 1; i >= 0; --i )
	 {
	    for( v = nvars - 1; v >= 0; --v )
	    {
	       assert(i == 1 || SCIPvarGetNegatedVar(allvars[v]) != NULL);
	       if( varcount[i][v] == 1 )
	       {
		  if( newandvars[0] == NULL )
		     newandvars[0] = i ? SCIPvarGetNegatedVar(allvars[v]) : allvars[v];
		  else
		  {
		     assert(newandvars[1] == NULL);
		     newandvars[1] = i ? SCIPvarGetNegatedVar(allvars[v]) : allvars[v];

		     breaked = TRUE;
		     break;
		  }
	       }
	    }

	    if( breaked )
	       break;
	 }
	 assert(newandvars[0] != NULL && newandvars[1] != NULL);

	 (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "andcons_%s_%s", SCIPconsGetName(cons), SCIPvarGetName(linvar));
	 SCIP_CALL( SCIPcreateConsAnd(scip, &newcons, name, linvar, nvars, newandvars,
	       TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	 SCIP_CALL( SCIPaddCons(scip, newcons) );
	 SCIPdebugPrintCons(scip, newcons, NULL);
	 SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

	 /* delete old constraints */
	 SCIP_CALL( SCIPdelCons(scip, consdata->lincons) );
	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 (*ndelconss) += 2;
      }
   }

   if( SCIPconsIsDeleted(cons) )
   {
      /* update the uses counter of consandata objects which are used in pseudoboolean constraint, which was deleted and
       * probably delete and-constraints
       */
      SCIP_CALL( updateConsanddataUses(scip, cons, conshdlrdata, ndelconss) );
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &negated);
   SCIPfreeBufferArray(scip, &repvars);
   SCIPfreeBufferArray(scip, &(varcount[1]));
   SCIPfreeBufferArray(scip, &(varcount[0]));
   SCIPfreeBufferArray(scip, &allvars);


   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyPseudoboolean)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrPseudoboolean(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreePseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check each constraint and get transformed constraints */
   for( c = conshdlrdata->nallconsanddatas - 1; c >= 0; --c )
   {
      SCIP_CONS* andcons;
      SCIP_VAR* resultant;
#ifndef NDEBUG
      SCIP_VAR** vars;
      int nvars;
      int v;

      assert(conshdlrdata->allconsanddatas[c] != NULL);
      assert(conshdlrdata->allconsanddatas[c]->newvars == NULL);

      vars = conshdlrdata->allconsanddatas[c]->vars;
      nvars = conshdlrdata->allconsanddatas[c]->nvars;
      assert(vars != NULL || nvars == 0);

      /* check for correct variables data */
      for( v = nvars - 1; v > 0; --v )
      {
         assert(SCIPvarIsTransformed(vars[v])); /*lint !e613*/
         assert(SCIPvarGetIndex(vars[v]) >= SCIPvarGetIndex(vars[v-1])); /*lint !e613*/
      }
      assert(nvars == 0 || SCIPvarIsTransformed(vars[0])); /*lint !e613*/
#endif

      andcons = conshdlrdata->allconsanddatas[c]->cons;
      assert(andcons != NULL);

      assert(SCIPconsIsTransformed(andcons));

      resultant = SCIPgetResultantAnd(scip, andcons);
      /* insert new mapping */
      assert(!SCIPhashmapExists(conshdlrdata->hashmap, (void*)resultant));
      SCIP_CALL( SCIPhashmapInsert(conshdlrdata->hashmap, (void*)resultant, (void*)(conshdlrdata->allconsanddatas[c])) );

      SCIPdebugMsg(scip, "insert into hashmap <%s> (%p) -> <%s> (%p/%p)\n", SCIPvarGetName(resultant), (void*)resultant,
         SCIPconsGetName(conshdlrdata->allconsanddatas[c]->cons), (void*)(conshdlrdata->allconsanddatas[c]),
         (void*)(conshdlrdata->allconsanddatas[c]->cons));
   }

   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitprePseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* decompose all pseudo boolean constraints into a "linear" constraint and "and" constraints */
   if( conshdlrdata->decomposeindicatorpbcons || conshdlrdata->decomposenormalpbcons )
   {
      for( c = 0; c < nconss; ++c )
      {
         SCIP_CONS* cons;
         SCIP_CONSDATA* consdata;
         SCIP_VAR** vars;
         SCIP_Real* coefs;
         int nvars;

         cons = conss[c];
         assert(cons != NULL);

	 /* only added constraints can be upgraded */
	 if( !SCIPconsIsAdded(cons) )
	    continue;

         consdata = SCIPconsGetData(cons);
         assert(consdata != NULL);

         /* gets number of variables in linear constraint */
         SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );

         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );

         /* get variables and coefficient of linear constraint */
         SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );
         assert(nvars == 0 || (coefs != NULL));

         if( consdata->issoftcons && conshdlrdata->decomposeindicatorpbcons )
         {
            SCIP_VAR* negindvar;
            char name[SCIP_MAXSTRLEN];
            SCIP_Real lhs;
            SCIP_Real rhs;
            SCIP_Bool initial;
            SCIP_Bool updateandconss;
            int v;
#if USEINDICATOR == FALSE
            SCIP_CONS* lincons;
            SCIP_Real maxact;
            SCIP_Real minact;
            SCIP_Real lb;
            SCIP_Real ub;
#else
            SCIP_CONS* indcons;
#endif

            assert(consdata->weight != 0);
            assert(consdata->indvar != NULL);

            /* if it is a soft constraint, there should be no integer variable */
            assert(consdata->intvar == NULL);

            /* get negation of indicator variable */
            SCIP_CALL( SCIPgetNegatedVar(scip, consdata->indvar, &negindvar) );
            assert(negindvar != NULL);

            /* get sides of linear constraint */
            SCIP_CALL( getLinearConsSides(scip, consdata->lincons, consdata->linconstype, &lhs, &rhs) );
            assert(!SCIPisInfinity(scip, lhs));
            assert(!SCIPisInfinity(scip, -rhs));
            assert(SCIPisLE(scip, lhs, rhs));

            updateandconss = FALSE;

#if USEINDICATOR == FALSE
            maxact = 0.0;
            minact = 0.0;

            /* adding all linear coefficients up */
            for( v = nvars - 1; v >= 0; --v )
               if( coefs[v] > 0 )
                  maxact += coefs[v];
               else
                  minact += coefs[v];

            if( SCIPisInfinity(scip, maxact) )
            {
               SCIPwarningMessage(scip, "maxactivity = %g exceed infinity value.\n", maxact);
            }
            if( SCIPisInfinity(scip, -minact) )
            {
               SCIPwarningMessage(scip, "minactivity = %g exceed -infinity value.\n", minact);
            }

            /* @todo check whether it's better to set the initial flag to false */
            initial = SCIPconsIsInitial(cons); /* FALSE; */

            /* first soft constraints for lhs */
            if( !SCIPisInfinity(scip, -lhs) )
            {
               /* first we are modelling the feasibility of the soft constraint by adding a slack variable */
               /* we ensure that if indvar == 1 => (a^T*x + ub*indvar >= lhs) */
               ub = lhs - minact;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_lhs_part1", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, nvars, vars, coefs, lhs, SCIPinfinity(scip),
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               /* update and constraint flags */
               SCIP_CALL( updateAndConss(scip, cons) );
               updateandconss = TRUE;

               /* add artificial indicator variable */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->indvar, ub) );

               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIPdebugPrintCons(scip, lincons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

               /* second we are modelling the implication that if the slack variable is on( negation is off), the constraint
                * is disabled, so only the cost arise if the slack variable is necessary */
               /* indvar == 1 => (a^T*x (+ ub * negindvar) <= lhs - 1) */
               ub = lhs - maxact - 1;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_lhs_part2", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, nvars, vars, coefs, -SCIPinfinity(scip), lhs - 1,
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               /* add artificial indicator variable */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, negindvar, ub) );

               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIPdebugPrintCons(scip, lincons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
            }

            /* second soft constraints for rhs */
            if( !SCIPisInfinity(scip, rhs) )
            {
               /* first we are modelling the feasibility of the soft-constraint by adding a slack variable */
               /* indvar == 1 => (a^T*x + lb * indvar <= rhs) */
               lb = rhs - maxact;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_rhs_part1", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, nvars, vars, coefs, -SCIPinfinity(scip), rhs,
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               if( !updateandconss )
               {
                  /* update and constraint flags */
                  SCIP_CALL( updateAndConss(scip, cons) );
               }

               /* add artificial indicator variable */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->indvar, lb) );

               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIPdebugPrintCons(scip, lincons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

               /* second we are modelling the implication that if the slack variable is on( negation is off), the constraint
                * is disabled, so only the cost arise if the slack variable is necessary */
               /* indvar == 1 => (a^T*x (+ lb * negindvar) >= rhs + 1) */
               lb = rhs - minact + 1;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_rhs_part2", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, nvars, vars, coefs, rhs + 1, SCIPinfinity(scip),
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               /* add artificial indicator variable */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, negindvar, lb) );

               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIPdebugPrintCons(scip, lincons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
            }
#else /* with indicator */
            /* @todo check whether it's better to set the initial flag to false */
            initial = SCIPconsIsInitial(cons); /* FALSE; */

            if( !SCIPisInfinity(scip, rhs) )
            {
               /* first we are modelling the implication that if the negation of the indicator variable is on, the constraint
                * is enabled */
               /* indvar == 0 => a^T*x <= rhs */

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_rhs_ind", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsIndicator(scip, &indcons, name, negindvar, nvars, vars, coefs, rhs,
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               /* update and constraint flags */
               SCIP_CALL( updateAndConss(scip, cons) );
               updateandconss = TRUE;

               SCIP_CALL( SCIPaddCons(scip, indcons) );
               SCIPdebugPrintCons(scip, indcons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &indcons) );
            }

            if( !SCIPisInfinity(scip, -lhs) )
            {
               /* second we are modelling the implication that if the negation of the indicator variable is on, the constraint
                * is enabled */
               /* change the a^T*x >= lhs to -a^Tx<= -lhs, for indicator constraint */

               for( v = nvars - 1; v >= 0; --v )
                  coefs[v] *= -1;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_lhs_ind", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsIndicator(scip, &indcons, name, negindvar, nvars, vars, coefs, -lhs,
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               if( !updateandconss )
               {
                  /* update and constraint flags */
                  SCIP_CALL( updateAndConss(scip, cons) );
               }

               SCIP_CALL( SCIPaddCons(scip, indcons) );
               SCIPdebugPrintCons(scip, indcons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &indcons) );
            }
#endif
            /* remove pseudo boolean and corresponding linear constraint, new linear constraints were created,
             * and-constraints still active
             */
            SCIP_CALL( SCIPdelCons(scip, consdata->lincons) );
            SCIP_CALL( SCIPdelCons(scip, cons) );
         }
         /* no soft constraint */
         else if( !consdata->issoftcons && conshdlrdata->decomposenormalpbcons )
         {
            /* todo: maybe better create a new linear constraint and let scip do the upgrade */

            /* mark linear constraint not to be upgraded - otherwise we loose control over it */
            SCIPconsAddUpgradeLocks(consdata->lincons, 1);

            /* update and constraint flags */
            SCIP_CALL( updateAndConss(scip, cons) );

#if 0 /* not implemented correctly */
            if( consdata->intvar != NULL )
            {
               /* add auxiliary integer variables to linear constraint */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->intvar, -1.0) );
            }
#endif
            /* remove pseudo boolean constraint, old linear constraint is still active, and-constraints too */
            SCIP_CALL( SCIPdelCons(scip, cons) );
         }

         /* free temporary memory */
         SCIPfreeBufferArray(scip, &coefs);
         SCIPfreeBufferArray(scip, &vars);
      }
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeletePseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool isorig;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   isorig = SCIPconsIsOriginal(cons);

   /* count number of used consanddata objects in original problem */
   if( isorig )
   {
#ifndef NDEBUG
      int c;
      assert((*consdata)->lincons == NULL || SCIPconsIsOriginal((*consdata)->lincons));

      for( c = (*consdata)->nconsanddatas - 1; c >= 0; --c )
      {
         assert((*consdata)->consanddatas[c]->nuses == 0);
         assert((*consdata)->consanddatas[c]->cons == NULL);
         assert((*consdata)->consanddatas[c]->noriguses == 0 || ((*consdata)->consanddatas[c]->origcons != NULL && SCIPconsIsOriginal((*consdata)->consanddatas[c]->origcons)));
      }
#endif
      conshdlrdata->noriguses -= (*consdata)->nconsanddatas;
   }
   assert(conshdlrdata->noriguses >= 0);

   /* free pseudo boolean constraint */
   SCIP_CALL( consdataFree(scip, consdata, isorig, conshdlrdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   SCIP_CONS** andconss;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   assert(sourcedata->nconsanddatas == 0 || sourcedata->consanddatas != NULL);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &andconss, sourcedata->nconsanddatas) );

   /* copy and-constraints */
   for( c = sourcedata->nconsanddatas - 1; c >= 0; --c )
   {
      assert(sourcedata->consanddatas[c] != NULL);
      andconss[c] = sourcedata->consanddatas[c]->origcons;
      assert(andconss[c] != NULL);
      assert(SCIPconsIsOriginal(andconss[c]));
   }

   /* create pseudoboolean constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, conshdlr, &targetdata, sourcedata->lincons, sourcedata->linconstype,
         andconss, sourcedata->andcoefs, sourcedata->andnegs, sourcedata->nconsanddatas, sourcedata->indvar, sourcedata->weight,
         sourcedata->issoftcons, sourcedata->intvar, sourcedata->lhs, sourcedata->rhs, SCIPconsIsChecked(sourcecons),
         TRUE) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &andconss);

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
SCIP_DECL_CONSENFOLP(consEnfolpPseudoboolean)
{  /*lint --e{715}*/
   SCIP_Bool violated;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   violated = FALSE;

   /* check all and-constraints */
   SCIP_CALL( checkAndConss(scip, conshdlr, NULL, &violated) );

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxPseudoboolean)
{  /*lint --e{715}*/
   SCIP_Bool violated;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   violated = FALSE;

   /* check all and-constraints */
   SCIP_CALL( checkAndConss(scip, conshdlr, sol, &violated) );

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsPseudoboolean)
{  /*lint --e{715}*/
   SCIP_Bool violated;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   violated = FALSE;

   /* check all and-constraints */
   SCIP_CALL( checkAndConss(scip, conshdlr, NULL, &violated) );

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckPseudoboolean)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(sol != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( nconss > 0 )
   {
      if( SCIPconsIsOriginal(conss[0]) )
      {
         SCIP_CONSDATA* consdata;

         for( c = nconss - 1; c >= 0 && (*result == SCIP_FEASIBLE || completely); --c )
         {
            consdata = SCIPconsGetData(conss[c]);
            assert(consdata != NULL);

            if( consdata->issoftcons )
            {
               assert(consdata->indvar != NULL);

               if( SCIPisEQ(scip, SCIPgetSolVal(scip, sol, consdata->indvar), 1.0) )
                  continue;
            }

            SCIP_CALL( checkOrigPbCons(scip, conss[c], sol, &violated, printreason) );
            if( violated )
               *result = SCIP_INFEASIBLE;
         }
      }
      else
      {
         /* check all and-constraints */
         SCIP_CALL( checkAndConss(scip, conshdlr, sol, &violated) );
         if( violated )
            *result = SCIP_INFEASIBLE;
      }
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int firstchange;
   int firstupgradetry;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* remember old preprocessing counters */
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* compute all changes in consanddata objects */
   SCIP_CALL( computeConsAndDataChanges(scip, conshdlrdata) );

   firstchange = INT_MAX;
   firstupgradetry = INT_MAX;
   cutoff = FALSE;

   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      SCIP_VAR** vars;
      SCIP_Real* coefs;
      int nvars;
      SCIP_VAR** linvars;
      SCIP_Real* lincoefs;
      int nlinvars;
      SCIP_VAR** andress;
      SCIP_Real* andcoefs;
      SCIP_Bool* andnegs;
      int nandress;
      SCIP_Real newlhs;
      SCIP_Real newrhs;

      cons = conss[c];
      assert(cons != NULL);
      assert(SCIPconsIsActive(cons));

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      assert(consdata->lincons != NULL);

      /* if linear constraint is redundant, than pseudoboolean constraint is redundant too */
      if( SCIPconsIsDeleted(consdata->lincons) )
      {
         /* update and constraint flags */
         SCIP_CALL( updateAndConss(scip, cons) );

         SCIP_CALL( SCIPdelCons(scip, cons) );
         ++(*ndelconss);
         continue;
      }

      /* get sides of linear constraint */
      SCIP_CALL( getLinearConsSides(scip, consdata->lincons, consdata->linconstype, &newlhs, &newrhs) );
      assert(!SCIPisInfinity(scip, newlhs));
      assert(!SCIPisInfinity(scip, -newrhs));
      assert(SCIPisLE(scip, newlhs, newrhs));

      /* gets number of variables in linear constraint */
      SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &andress, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &andcoefs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &andnegs, nvars) );

      /* get variables and coefficient of linear constraint */
      SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );
      assert(nvars == 0 || (coefs != NULL));

      /* calculate all not artificial linear variables and all artificial and-resultants which will be ordered like the
       * 'consanddatas' such that the and-resultant of the and-constraint is the and-resultant in the 'andress' array
       * afterwards
       */
      SCIP_CALL( getLinVarsAndAndRess(scip, cons, vars, coefs, nvars, linvars, lincoefs, &nlinvars,
            andress, andcoefs, andnegs, &nandress) );

      /* update all locks inside this constraint and all captures on all and-constraints */
      SCIP_CALL( correctLocksAndCaptures(scip, cons, conshdlrdata, newlhs, newrhs, andress, andcoefs, andnegs, nandress) );

      /* we can only presolve pseudoboolean constraints, that are not modifiable */
      if( SCIPconsIsModifiable(cons) )
         goto CONTTERMINATE;

      SCIPdebugMsg(scip, "presolving pseudoboolean constraint <%s>\n", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);

      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      if( consdata->changed && !SCIPisStopped(scip) )
      {
	 /* check if we can aggregated some variables */
	 SCIP_CALL( findAggregation(scip, cons, conshdlrdata, ndelconss, naggrvars, &cutoff) );
      }

      /* if aggregation also deleted the constraint we can go to the next */
      if( !SCIPconsIsActive(cons) )
         goto CONTTERMINATE;

      if( consdata->changed )
      {
         /* try upgrading pseudoboolean constraint to a linear constraint and/or remove possible and-constraints */
         SCIP_CALL( tryUpgrading(scip, cons, conshdlrdata, ndelconss, naddconss, nfixedvars, nchgcoefs, nchgsides, &cutoff) );
         if( cutoff )
            goto CONTTERMINATE;
      }

      /* if upgrading deleted the pseudoboolean constraint we go on */
      if( !SCIPconsIsActive(cons) )
         goto CONTTERMINATE;

      /* remember the first constraint that was not yet tried to be upgraded, to begin the next upgrading round with */
      if( firstupgradetry == INT_MAX && !consdata->upgradetried )
         firstupgradetry = c;

      while( !consdata->presolved && !SCIPisStopped(scip) )
      {
         /* mark constraint being presolved and propagated */
         consdata->presolved = TRUE;

         /* add cliques to the clique table */
         SCIP_CALL( addCliques(scip, cons, &cutoff, naggrvars, nchgbds) );
         if( cutoff )
            break;

         /* propagate constraint */
         SCIP_CALL( propagateCons(scip, cons, &cutoff, ndelconss) );
         if( cutoff )
            break;
      }

   CONTTERMINATE:

      /* reset changed flag */
      if( SCIPconsIsActive(cons) )
      {
	 consdata->changed = FALSE;
      }

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &andnegs);
      SCIPfreeBufferArray(scip, &andcoefs);
      SCIPfreeBufferArray(scip, &andress);
      SCIPfreeBufferArray(scip, &lincoefs);
      SCIPfreeBufferArray(scip, &linvars);
      SCIPfreeBufferArray(scip, &coefs);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* delete unused information in constraint handler data */
   SCIP_CALL( correctConshdlrdata(scip, conshdlrdata, ndelconss) );

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds || *ndelconss > oldndelconss
      || *nupgdconss > oldnupgdconss || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int v;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhs = consdata->lhs;
   rhs = consdata->rhs;
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(SCIPisLE(scip, lhs, rhs));

   haslhs = !SCIPisInfinity(scip, -lhs);
   hasrhs = !SCIPisInfinity(scip, rhs);

   SCIPdebugMsg(scip, "%socking constraint <%s> by [%d;%d].\n", (nlocksneg < 0) || (nlockspos < 0) ? "Unl" : "L", SCIPconsGetName(cons), nlocksneg, nlockspos);

   /* update rounding locks of every single variable corresponding to the and-constraints */
   for( c = consdata->nconsanddatas - 1; c >= 0; --c )
   {
      SCIP_VAR* andres;
      SCIP_VAR** andvars;
      SCIP_Real val;
      int nandvars;
      SCIP_CONS* andcons;
      CONSANDDATA* consanddata;

      consanddata = consdata->consanddatas[c];
      assert( consanddata != NULL );

      if( !consanddata->istransformed )
         continue;

      andcons = consanddata->cons;

      if( andcons == NULL )
      {
         /* we should have no new variables */
         assert(consanddata->nnewvars == 0);
         assert(consanddata->nvars == 0);

         SCIPfreeBlockMemoryArrayNull(scip, &(consanddata->vars), consanddata->svars);
         SCIPfreeBlockMemoryArrayNull(scip, &(consanddata->newvars), consanddata->snewvars);

         consanddata->nvars = 0;
         consanddata->svars = 0;
         consanddata->nnewvars = 0;
         consanddata->snewvars = 0;
         consanddata->istransformed = FALSE;

         continue;
      }
      assert(andcons != NULL);
      if( consanddata->nnewvars > 0 )
      {
         andvars = consanddata->newvars;
         nandvars = consanddata->nnewvars;
      }
      else
      {
         andvars = consanddata->vars;
         nandvars = consanddata->nvars;
      }

      /* probably we need to store the resultant too, now it's not possible to remove the resultant from the and-constraint */
      andres = SCIPgetResultantAnd(scip, andcons);
      assert(nandvars == 0 || andvars != NULL);
      assert(andres != NULL);
      val = consdata->andnegs[c] ? -consdata->andcoefs[c] : consdata->andcoefs[c];

      /* lock variables */
      if( SCIPisPositive(scip, val) )
      {
         if( haslhs )
         {
            for( v = nandvars - 1; v >= 0; --v )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, andvars[v], nlockspos, nlocksneg) );
            }
            SCIP_CALL( SCIPaddVarLocks(scip, andres, nlocksneg + nlockspos, nlocksneg + nlockspos) );

            SCIP_CALL( checkLocksAndRes(scip, andres) );
         }
         if( hasrhs )
         {
            for( v = nandvars - 1; v >= 0; --v )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, andvars[v], nlocksneg, nlockspos) );
            }
            /* don't double the locks on the and-resultant */
            if( !haslhs )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, andres, nlocksneg + nlockspos, nlocksneg + nlockspos) );

               SCIP_CALL( checkLocksAndRes(scip, andres) );
            }
         }
      }
      else
      {
         if( haslhs )
         {
            for( v = nandvars - 1; v >= 0; --v )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, andvars[v], nlocksneg, nlockspos) );
            }
            SCIP_CALL( SCIPaddVarLocks(scip, andres, nlocksneg + nlockspos, nlocksneg + nlockspos) );

            SCIP_CALL( checkLocksAndRes(scip, andres) );
         }
         if( hasrhs )
         {
            for( v = nandvars - 1; v >= 0; --v )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, andvars[v], nlockspos, nlocksneg) );
            }
            /* don't double the locks on the and-resultant */
            if( !haslhs )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, andres, nlocksneg + nlockspos, nlocksneg + nlockspos) );

               SCIP_CALL( checkLocksAndRes(scip, andres) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintPseudoboolean)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   SCIP_CALL( consdataPrint(scip, cons, file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyPseudoboolean)
{  /*lint --e{715}*/
   const char* consname;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIP_CALL( copyConsPseudoboolean(scip, cons, sourcescip, sourcecons, consname, varmap, consmap,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global,
         valid) );
   assert(cons != NULL || *valid == FALSE);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   CONSANDDATA* consanddata;
   SCIP_VAR** linconsvars;
   SCIP_VAR** linvars;
   SCIP_VAR** andress;
   int nlinconsvars;
   int nlinvars;
   int nandress;
   SCIP_Bool transformed;
   int nvars;
   int r;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(vars != NULL);
   assert(varssize >= 0);
   assert(success != NULL);

   if( varssize < 0 )
      return SCIP_INVALIDDATA;

   (*success) = TRUE;

   /* pseudoboolean constraint is already deleted */
   if( SCIPconsIsDeleted(cons) )
   {
      vars = NULL;

      return SCIP_OKAY; /*lint !e438*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lincons != NULL);

   /* linear constraint of pseudoboolean is already deleted */
   if( SCIPconsIsDeleted(consdata->lincons) )
   {
      vars = NULL;

      return SCIP_OKAY; /*lint !e438*/
   }

   /* gets number of variables in linear constraint */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nlinconsvars) );
   assert(nlinconsvars >= 0);

   /* no variables exist */
   if( nlinconsvars == 0 )
   {
      vars = NULL;

      return SCIP_OKAY; /*lint !e438*/
   }
   /* not enough space in the variables array */
   else if( varssize < nlinconsvars )
   {
      (*success) = FALSE;

      return SCIP_OKAY;
   }

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &linconsvars, nlinconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nlinconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andress, nlinconsvars) );

   /* get variables and coefficient of linear constraint */
   SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, linconsvars, NULL, &nlinconsvars) );

   /* calculate all non-artificial linear variables and all artificial and-resultants which will be ordered like the
    * 'consanddatas' such that the and-resultant of the and-constraint is the and-resultant in the 'andress' array
    * afterwards
    */
   SCIP_CALL( getLinVarsAndAndRess(scip, cons, linconsvars, NULL, nlinconsvars, linvars, NULL, &nlinvars,
         andress, NULL, NULL, &nandress) );
   assert(nlinconsvars == nlinvars + nandress);

   nvars = nlinvars;

   if( nlinvars > 0 )
   {
      assert(linvars != NULL);
      BMScopyMemoryArray(vars, linvars, nvars);
   }

   if( nandress == 0 )
      goto TERMINATE;

   assert(andress != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->hashmap != NULL);

   transformed = SCIPconsIsTransformed(cons);

   for( r = nandress - 1; r >= 0; --r )
   {
      SCIP_CONS* andcons;

      assert(andress[r] != NULL);

      consanddata = (CONSANDDATA*) SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)andress[r]);

      assert(consanddata != NULL);
      assert(consanddata->istransformed);

      if( transformed )
         andcons = consanddata->cons;
      else
         andcons = consanddata->origcons;

      assert(andcons != NULL);

      /* not enough space for all variables */
      if( varssize <= nvars )
      {
	 (*success) = FALSE;

	 goto TERMINATE;
      }

      /* add the resultant */
      vars[nvars] = andress[r];
      ++nvars;

      /* add all and-operands and the resultant */
      if( !SCIPconsIsDeleted(andcons) )
      {
	 int noperands = SCIPgetNVarsAnd(scip, andcons);

	 assert(noperands >= 0);

	 /* not enough space for all variables */
	 if( varssize < nvars + noperands )
	 {
	    (*success) = FALSE;

	    goto TERMINATE;
	 }

	 /* copy operands */
	 if( noperands > 0 )
	 {
	    assert(SCIPgetVarsAnd(scip, andcons) != NULL);
	    BMScopyMemoryArray(&(vars[nvars]), SCIPgetVarsAnd(scip, andcons), noperands); /*lint !e866*/
	    nvars += noperands;
	 }
      }
   }

 TERMINATE:

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &andress);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &linconsvars);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   CONSANDDATA* consanddata;
   SCIP_VAR** linconsvars;
   SCIP_VAR** linvars;
   SCIP_VAR** andress;
   int nlinconsvars;
   int nlinvars;
   int nandress;
   SCIP_Bool transformed;
   int r;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(nvars != NULL);
   assert(success != NULL);

   (*success) = TRUE;

   /* pseudoboolean constraint is already deleted */
   if( SCIPconsIsDeleted(cons) )
   {
      *nvars = 0;

      return SCIP_OKAY;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lincons != NULL);

   /* linear constraint of pseudoboolean is already deleted */
   if( SCIPconsIsDeleted(consdata->lincons) )
   {
      *nvars = 0;

      return SCIP_OKAY;
   }

   /* gets number of variables in linear constraint */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nlinconsvars) );
   assert(nlinconsvars >= 0);

   /* no variables exist */
   if( nlinconsvars == 0 )
   {
      *nvars = 0;

      return SCIP_OKAY;
   }

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &linconsvars, nlinconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nlinconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andress, nlinconsvars) );

   /* get variables and coefficient of linear constraint */
   SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, linconsvars, NULL, &nlinconsvars) );

   /* calculate all non-artificial linear variables and all artificial and-resultants which will be ordered like the
    * 'consanddatas' such that the and-resultant of the and-constraint is the and-resultant in the 'andress' array
    * afterwards
    */
   SCIP_CALL( getLinVarsAndAndRess(scip, cons, linconsvars, NULL, nlinconsvars, linvars, NULL, &nlinvars,
         andress, NULL, NULL, &nandress) );
   assert(nlinconsvars == nlinvars + nandress);

   *nvars = nlinvars;

   if( nandress == 0 )
      goto TERMINATE;

   assert(andress != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->hashmap != NULL);

   transformed = SCIPconsIsTransformed(cons);

   for( r = nandress - 1; r >= 0; --r )
   {
      SCIP_CONS* andcons;

      assert(andress[r] != NULL);

      consanddata = (CONSANDDATA*) SCIPhashmapGetImage(conshdlrdata->hashmap, (void*)andress[r]);

      assert(consanddata != NULL);
      assert(consanddata->istransformed);

      if( transformed )
         andcons = consanddata->cons;
      else
         andcons = consanddata->origcons;

      assert(andcons != NULL);

      if( SCIPconsIsDeleted(andcons) )
      {
	 /* only add one for the resultant */
	 ++(*nvars);
      }
      else
      {
	 /* add all and-operands and one for the resultant */
	 *nvars += SCIPgetNVarsAnd(scip, andcons) + 1;
      }
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &andress);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &linconsvars);

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for pseudoboolean constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrPseudoboolean(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create pseudoboolean constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpPseudoboolean, consEnfopsPseudoboolean, consCheckPseudoboolean, consLockPseudoboolean,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyPseudoboolean, consCopyPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeletePseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreePseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitprePseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolPseudoboolean, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxPseudoboolean) );

   /* add pseudoboolean constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/decomposenormal",
         "decompose all normal pseudo boolean constraint into a \"linear\" constraint and \"and\" constraints",
         &conshdlrdata->decomposenormalpbcons, TRUE, DEFAULT_DECOMPOSENORMALPBCONS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/decomposeindicator",
         "decompose all indicator pseudo boolean constraint into a \"linear\" constraint and \"and\" constraints",
         &conshdlrdata->decomposeindicatorpbcons, TRUE, DEFAULT_DECOMPOSEINDICATORPBCONS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/nlcseparate", "should the nonlinear constraints be separated during LP processing?",
         NULL, TRUE, DEFAULT_SEPARATENONLINEAR, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/nlcpropagate", "should the nonlinear constraints be propagated during node processing?",
         NULL, TRUE, DEFAULT_PROPAGATENONLINEAR, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/nlcremovable", "should the nonlinear constraints be removable?",
         NULL, TRUE, DEFAULT_REMOVABLENONLINEAR, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a pseudoboolean constraint, with given linear and and-constraints */
SCIP_RETCODE SCIPcreateConsPseudobooleanWithConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONS*            lincons,            /**< associated linear constraint */
   SCIP_LINEARCONSTYPE   linconstype,        /**< linear constraint type of associated linear constraint */
   SCIP_CONS**           andconss,           /**< associated and-constraints */
   SCIP_Real*            andcoefs,           /**< associated coefficients of and-constraints */
   int                   nandconss,          /**< number of associated and-constraints */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< an artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
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
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   CONSANDDATA* newdata;
   CONSANDDATA* tmpdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR* res;
   SCIP_Bool memisinvalid;
   SCIP_Bool transformed;
   int nvars;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(lincons != NULL);
   assert(linconstype > SCIP_LINEARCONSTYPE_INVALIDCONS);
   assert(andconss != NULL);
   assert(andcoefs != NULL);
   assert(nandconss >= 1);
   assert(issoftcons == (indvar != NULL));

   /* find the pseudoboolean constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("pseudo boolean constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* initial hashmap and -table */
   SCIP_CALL( inithashmapandtable(scip, &conshdlrdata) );

   assert(conshdlrdata->hashmap != NULL);
   assert(conshdlrdata->hashtable != NULL);
   assert(conshdlrdata->allconsanddatas != NULL);
   assert(conshdlrdata->nallconsanddatas <= conshdlrdata->sallconsanddatas);

   memisinvalid = TRUE;
   newdata = NULL;

   transformed = SCIPconsIsTransformed(lincons);

   /* create hash map and hash table entries */
   for( c = nandconss - 1; c >= 0; --c )
   {
      assert(andconss[c] != NULL);
      res = SCIPgetResultantAnd(scip, andconss[c]);
      vars = SCIPgetVarsAnd(scip, andconss[c]);
      nvars = SCIPgetNVarsAnd(scip, andconss[c]);
      assert(vars != NULL && nvars > 0);
      assert(res != NULL);

      /* stop if the constraint has 0 variables or an error occurred (coverity issue) */
      if( nvars <= 0 )
         continue;

      /* if allocated memory in this for loop was already used, allocate a new block, otherwise we only need to copy the variables */
      if( memisinvalid )
      {
         /* allocate memory for a possible new consanddata object */
         SCIP_CALL( SCIPallocBlockMemory(scip, &newdata) );
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(newdata->vars), vars, nvars) );
         newdata->svars = nvars;
         newdata->newvars = NULL;
         newdata->nnewvars = 0;
         newdata->snewvars = 0;
         newdata->istransformed = transformed;
         newdata->isoriginal = !transformed;
         newdata->noriguses = 0;
         newdata->nuses = 0;
         newdata->cons = NULL;
         newdata->origcons = NULL;
      }
      else
      {
         assert(newdata != NULL);
         /* resize variable array if necessary */
         if( newdata->svars < nvars )
         {
            SCIP_CALL( SCIPensureBlockMemoryArray(scip, &(newdata->vars), &(newdata->svars), nvars) );
         }

         /* copy variables in already allocated array */
         BMScopyMemoryArray(newdata->vars, vars, nvars);
      }

      /* sort variables */
      SCIPsortPtr((void**)(newdata->vars), SCIPvarComp, nvars);

      newdata->nvars = nvars;
      assert(newdata->vars != NULL && newdata->nvars > 0);

      if( SCIPconsIsTransformed(andconss[c]) )
      {
         int v;

         /* either all constraints are transformed or none */
         assert(transformed);
         newdata->cons = andconss[c];

         /* capture all variables */
         for( v = newdata->nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( SCIPcaptureVar(scip, newdata->vars[v]) ); /*lint !e613*/
         }
      }
      else
      {
         /* either all constraints are transformed or none */
         assert(!transformed);
         newdata->origcons = andconss[c];
      }

      /* get constraint from current hash table with same variables as andconss[c] */
      tmpdata = (CONSANDDATA*)(SCIPhashtableRetrieve(conshdlrdata->hashtable, (void*)newdata));
      assert(tmpdata == NULL || tmpdata->cons != NULL || tmpdata->origcons != NULL);

      if( tmpdata == NULL || (tmpdata->cons != andconss[c] && tmpdata->origcons != andconss[c]))
      {
         if( tmpdata != NULL && (tmpdata->cons != NULL || tmpdata->origcons != NULL) )
         {
            SCIPwarningMessage(scip, "Another and-constraint with the same variables but different and-resultant is added to the global and-constraint hashtable of pseudoboolean constraint handler.\n");
         }

         /* resize data for all and-constraints if necessary */
         if( conshdlrdata->nallconsanddatas == conshdlrdata->sallconsanddatas )
         {
            SCIP_CALL( SCIPensureBlockMemoryArray(scip, &(conshdlrdata->allconsanddatas), &(conshdlrdata->sallconsanddatas), SCIPcalcMemGrowSize(scip, conshdlrdata->sallconsanddatas + 1)) );
         }

         conshdlrdata->allconsanddatas[conshdlrdata->nallconsanddatas] = newdata;
         ++(conshdlrdata->nallconsanddatas);

         /* no such and-constraint in current hash table: insert the new object into hash table */
         SCIP_CALL( SCIPhashtableInsert(conshdlrdata->hashtable, (void*)newdata) );

         /* if newdata object was new we want to allocate new memory in next loop iteration */
         memisinvalid = TRUE;
         assert(!SCIPhashmapExists(conshdlrdata->hashmap, (void*)res));

         /* capture and-constraint */
         if( transformed )
         {
            SCIP_CALL( SCIPcaptureCons(scip, newdata->cons) );

            /* initialize usage of data object */
            newdata->nuses = 1;
         }
         else
         {
            SCIP_CALL( SCIPcaptureCons(scip, newdata->origcons) );

            /* initialize usage of data object */
            newdata->noriguses = 1;
         }

         /* insert new mapping */
         assert(!SCIPhashmapExists(conshdlrdata->hashmap, (void*)res));
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->hashmap, (void*)res, (void*)newdata) );
      }
      else
      {
         assert(SCIPhashmapExists(conshdlrdata->hashmap, (void*)res));
         memisinvalid = FALSE;

         if( transformed )
         {
            assert(tmpdata->nuses > 0);

            /* increase usage of data object */
            ++(tmpdata->nuses);
         }
         else
         {
            assert(tmpdata->noriguses > 0);

            /* increase usage of data object */
            ++(tmpdata->noriguses);
         }
      }
   }

   if( !memisinvalid )
   {
      assert(newdata != NULL);

      /* free temporary memory */
      SCIPfreeBlockMemoryArray(scip, &(newdata->vars), newdata->svars);
      SCIPfreeBlockMemory(scip, &newdata);
   }

   /* adjust right hand side */
   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -rhs) )
      rhs = -SCIPinfinity(scip);

   /* capture linear constraint */
   SCIP_CALL( SCIPcaptureCons(scip, lincons) );

   /* todo: make the constraint upgrade flag global, now it works only for the common linear constraint */
   /* mark linear constraint not to be upgraded - otherwise we loose control over it */
   SCIPconsAddUpgradeLocks(lincons, 1);

   /* create constraint data */
   /* checking for and-constraints will be FALSE, we check all information in this constraint handler */
   SCIP_CALL( consdataCreate(scip, conshdlr, &consdata, lincons, linconstype, andconss, andcoefs, NULL, nandconss,
         indvar, weight, issoftcons, intvar, lhs, rhs, check, FALSE) );
   assert(consdata != NULL);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a pseudoboolean constraint
 *
 *  @note linear and nonlinear terms can be added using SCIPaddCoefPseudoboolean() and SCIPaddTermPseudoboolean(),
 *        respectively
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            linvars,            /**< variables of the linear part, or NULL */
   int                   nlinvars,           /**< number of variables of the linear part */
   SCIP_Real*            linvals,            /**< coefficients of linear part, or NULL */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   SCIP_Real*            termvals,           /**< coefficients of nonlinear parts, or NULL */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< an artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
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
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** andress;
   SCIP_CONS** andconss;
   SCIP_Real* andcoefs;
   SCIP_Bool* andnegs;
   int nandconss;
   SCIP_CONS* lincons;
   SCIP_LINEARCONSTYPE linconstype;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nlinvars == 0 || (linvars != NULL && linvals != NULL));
   assert(nterms == 0 || (terms != NULL && termvals != NULL && ntermvars != NULL));
   assert(issoftcons == (indvar != NULL));

   /* find the pseudoboolean constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("pseudo boolean constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

#if USEINDICATOR == TRUE
   if( issoftcons && modifiable )
   {
      SCIPerrorMessage("Indicator constraint handler can't work with modifiable constraints\n");
      return SCIP_INVALIDDATA;
   }
#endif

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* initial hashmap and -table */
   SCIP_CALL( inithashmapandtable(scip, &conshdlrdata) );

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &andconss, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andress, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andcoefs, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &andnegs, nterms) );

   nandconss = 0;
   /* create and-constraints */
   SCIP_CALL( createAndAddAnds(scip, conshdlr, terms, termvals, nterms, ntermvars,
         initial, enforce, check, local, modifiable, dynamic, stickingatnode,
         andconss, andcoefs, andnegs, &nandconss) );
   assert(nterms >= nandconss);

   /* get all and-resultants for linear constraint */
   for( c = nandconss - 1; c >= 0; --c )
   {
      assert(andconss[c] != NULL);
      andress[c] = SCIPgetResultantAnd(scip, andconss[c]);
   }

   linconstype = SCIP_LINEARCONSTYPE_INVALIDCONS;

   /* adjust right hand side */
   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -rhs) )
      rhs = -SCIPinfinity(scip);

   /* create and add linear constraint */
   /* checking for original linear constraint will be FALSE, transformed linear constraints get the check flag like this
    * pseudoboolean constraint, in this constraint handler we only will check all and-constraints
    */
   SCIP_CALL( createAndAddLinearCons(scip, conshdlr, linvars, nlinvars, linvals, andress, nandconss, andcoefs, andnegs,
         &lhs, &rhs, initial, separate, enforce, FALSE/*check*/, propagate, local, modifiable, dynamic, removable,
         stickingatnode, &lincons, &linconstype) );
   assert(lincons != NULL);
   assert(linconstype > SCIP_LINEARCONSTYPE_INVALIDCONS);

   /* create constraint data */
   /* checking for and-constraints will be FALSE, we check all information in this constraint handler */
   SCIP_CALL( consdataCreate(scip, conshdlr, &consdata, lincons, linconstype, andconss, andcoefs, andnegs, nandconss,
         indvar, weight, issoftcons, intvar, lhs, rhs, check, FALSE) );
   assert(consdata != NULL);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &andnegs);
   SCIPfreeBufferArray(scip, &andcoefs);
   SCIPfreeBufferArray(scip, &andress);
   SCIPfreeBufferArray(scip, &andconss);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a pseudoboolean constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            linvars,            /**< variables of the linear part, or NULL */
   int                   nlinvars,           /**< number of variables of the linear part */
   SCIP_Real*            linvals,            /**< coefficients of linear part, or NULL */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   SCIP_Real*            termvals,           /**< coefficients of nonlinear parts, or NULL */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsPseudoboolean(scip, cons, name, linvars, nlinvars, linvals,
         terms, nterms, ntermvars, termvals, indvar, weight, issoftcons, intvar, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** adds a variable to the pseudo boolean constraint (if it is not zero)
 *
 * @note  you can only add a coefficient if the special type of linear constraint won't changed
 *
 * @todo  if adding a coefficient would change the type of the special linear constraint, we need to erase it and
 *         create a new linear constraint
 */
SCIP_RETCODE SCIPaddCoefPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< constraint data */
   SCIP_VAR*const        var,                /**< variable of constraint entry */
   SCIP_Real const       val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   if( SCIPisZero(scip, val) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->linconstype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      SCIP_CALL( SCIPaddCoefLinear(scip, consdata->lincons, var, val) );
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
      if( !SCIPisEQ(scip, val, 1.0) )
         return SCIP_INVALIDDATA;

      SCIP_CALL( SCIPaddCoefLogicor(scip, consdata->lincons, var) );
      break;
   case SCIP_LINEARCONSTYPE_KNAPSACK:
      if( !SCIPisIntegral(scip, val) || !SCIPisPositive(scip, val) )
         return SCIP_INVALIDDATA;

      SCIP_CALL( SCIPaddCoefKnapsack(scip, consdata->lincons, var, (SCIP_Longint) val) );
      break;
   case SCIP_LINEARCONSTYPE_SETPPC:
      if( !SCIPisEQ(scip, val, 1.0) )
         return SCIP_INVALIDDATA;

      SCIP_CALL( SCIPaddCoefSetppc(scip, consdata->lincons, var) );
      break;
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
      if( !SCIPisIntegral(scip, val) || !SCIPisPositive(scip, val) )
         return SCIP_INVALIDDATA;

      SCIP_CALL( SCIPaddCoefEQKnapsack(scip, consdata->lincons, var, (SCIP_Longint) val) );
      break;
#endif
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   consdata->propagated = FALSE;
   consdata->presolved = FALSE;
   consdata->cliquesadded = FALSE;

   return SCIP_OKAY;
}

/** adds nonlinear term to pseudo boolean constraint (if it is not zero)
 *
 * @note  you can only add a coefficient if the special type of linear constraint won't changed
 *
 * @todo if adding a coefficient would change the type of the special linear constraint, we need to erase it and
 *         create a new linear constraint
 */
SCIP_RETCODE SCIPaddTermPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_VAR**const       vars,               /**< variables of the nonlinear term */
   int const             nvars,              /**< number of variables of the nonlinear term */
   SCIP_Real const       val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nvars == 0 || vars != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   SCIP_CALL( addCoefTerm(scip, cons, vars, nvars, val) );

   return SCIP_OKAY;
}

/** gets indicator variable of pseudoboolean constraint, or NULL if there is no */
SCIP_VAR* SCIPgetIndVarPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->indvar;
}

/** gets linear constraint of pseudoboolean constraint */
SCIP_CONS* SCIPgetLinearConsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lincons;
}

/** gets type of linear constraint of pseudoboolean constraint */
SCIP_LINEARCONSTYPE SCIPgetLinearConsTypePseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return SCIP_LINEARCONSTYPE_INVALIDCONS; /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->linconstype;
}

/** gets number of linear variables without artificial terms variables of pseudoboolean constraint */
int SCIPgetNLinVarsWithoutAndPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   checkConsConsistency(scip, cons);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nlinvars;
}

/** gets linear constraint of pseudoboolean constraint */
SCIP_RETCODE SCIPgetLinDatasWithoutAndPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_VAR**const       linvars,            /**< array to store and-constraints */
   SCIP_Real*const       lincoefs,           /**< array to store and-coefficients */
   int*const             nlinvars            /**< pointer to store the required array size for and-constraints, have to
                                              *   be initialized with size of given array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nlinvars != NULL);
   assert(*nlinvars == 0 || linvars != NULL);
   assert(*nlinvars == 0 || lincoefs != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   checkConsConsistency(scip, cons);

   if( *nlinvars < consdata->nlinvars )
   {
      *nlinvars = consdata->nlinvars;
      return SCIP_OKAY;
   }

   /* gets number of variables in linear constraint */
   SCIP_CALL( getLinearConsNVars(scip, consdata->lincons, consdata->linconstype, &nvars) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );

   /* get variables and coefficient of linear constraint */
   SCIP_CALL( getLinearConsVarsData(scip, consdata->lincons, consdata->linconstype, vars, coefs, &nvars) );

   /* calculate all not artificial linear variables */
   SCIP_CALL( getLinVarsAndAndRess(scip, cons, vars, coefs, nvars, linvars, lincoefs, nlinvars, NULL, NULL, NULL, NULL) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** gets and-constraints of pseudoboolean constraint */
SCIP_RETCODE SCIPgetAndDatasPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONS**const      andconss,           /**< array to store and-constraints */
   SCIP_Real*const       andcoefs,           /**< array to store and-coefficients */
   int*const             nandconss           /**< pointer to store the required array size for and-constraints, have to
                                              *   be initialized with size of given array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool isorig;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nandconss != NULL);
   assert(*nandconss == 0 || andconss != NULL);
   assert(*nandconss == 0 || andcoefs != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   checkConsConsistency(scip, cons);

   if( *nandconss < consdata->nconsanddatas )
   {
      *nandconss = consdata->nconsanddatas;
      return SCIP_OKAY;
   }

   *nandconss = consdata->nconsanddatas;
   assert(*nandconss == 0 || consdata->consanddatas != NULL);

   isorig = SCIPconsIsOriginal(cons);

   for( c = *nandconss - 1; c >= 0; --c )
   {
      assert(consdata->consanddatas[c] != NULL);
      assert(consdata->consanddatas[c]->istransformed ? (consdata->consanddatas[c]->cons != NULL) : TRUE);
      assert(consdata->consanddatas[c]->isoriginal ? (consdata->consanddatas[c]->origcons != NULL) : TRUE);
      assert(consdata->consanddatas[c]->cons != NULL || consdata->consanddatas[c]->origcons != NULL);
      assert(isorig ? consdata->consanddatas[c]->origcons != NULL : consdata->consanddatas[c]->cons != NULL);

      andconss[c] = (isorig ? consdata->consanddatas[c]->origcons : consdata->consanddatas[c]->cons);
      assert(andconss[c] != NULL);

      andcoefs[c] = consdata->andcoefs[c];
   }

   return SCIP_OKAY;
}

/** gets number of and constraints of pseudoboolean constraint */
int SCIPgetNAndsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   checkConsConsistency(scip, cons);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nconsanddatas;
}

/** changes left hand side of pseudoboolean constraint
 *
 * @note you can only change the left hand side if the special type of linear constraint won't changed
 *
 * @todo if changing the left hand side would change the type of the special linear constraint, we need to erase it
 *       and create a new linear constraint
 */
SCIP_RETCODE SCIPchgLhsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< constraint data */
   SCIP_Real const       lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      return SCIP_INVALIDDATA;
   }

   checkConsConsistency(scip, cons);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->linconstype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      SCIP_CALL( chgLhs(scip, cons, lhs) );
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
   case SCIP_LINEARCONSTYPE_KNAPSACK:
   case SCIP_LINEARCONSTYPE_SETPPC:
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
#endif
      SCIPerrorMessage("changing left hand side only allowed on standard linear constraint \n");
      return SCIP_INVALIDDATA;
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes right hand side of pseudoboolean constraint
 *
 * @note you can only change the right hand side if the special type of linear constraint won't changed
 *
 * @todo if changing the right hand side would change the type of the special linear constraint, we need to erase it
 *       and create a new linear constraint
 */
SCIP_RETCODE SCIPchgRhsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< constraint data */
   SCIP_Real const       rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      return SCIP_INVALIDDATA;
   }

   checkConsConsistency(scip, cons);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->linconstype )
   {
   case SCIP_LINEARCONSTYPE_LINEAR:
      SCIP_CALL( chgRhs(scip, cons, rhs) );
      break;
   case SCIP_LINEARCONSTYPE_LOGICOR:
   case SCIP_LINEARCONSTYPE_KNAPSACK:
   case SCIP_LINEARCONSTYPE_SETPPC:
#ifdef WITHEQKNAPSACK
   case SCIP_LINEARCONSTYPE_EQKNAPSACK:
#endif
      SCIPerrorMessage("changing right hand side only allowed on standard linear constraint \n");
      return SCIP_INVALIDDATA;
   case SCIP_LINEARCONSTYPE_INVALIDCONS:
   default:
      SCIPerrorMessage("unknown linear constraint type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** get left hand side of pseudoboolean constraint */
SCIP_Real SCIPgetLhsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }

   checkConsConsistency(scip, cons);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** get right hand side of pseudoboolean constraint */
SCIP_Real SCIPgetRhsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }

   checkConsConsistency(scip, cons);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}
