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

/**@file   heur_crossover.c
 * @brief  crossover primal heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/heur_crossover.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "crossover"
#define HEUR_DESC             "LNS heuristic that fixes all variables that are identic in a couple of solutions"
#define HEUR_DISPCHAR         'C'
#define HEUR_PRIORITY         -1104000
#define HEUR_FREQ             30
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL         /* maximum number of nodes to regard in the subproblem                   */
#define DEFAULT_MINIMPROVE    0.01           /* factor by which Crossover should at least improve the incumbent       */
#define DEFAULT_MINNODES      50LL           /* minimum number of nodes to regard in the subproblem                   */
#define DEFAULT_MINFIXINGRATE 0.666          /* minimum percentage of integer variables that have to be fixed         */
#define DEFAULT_NODESOFS      500LL          /* number of nodes added to the contingent of the total nodes            */
#define DEFAULT_NODESQUOT     0.1            /* subproblem nodes in relation to nodes of the original problem         */
#define DEFAULT_LPLIMFAC      2.0            /* factor by which the limit on the number of LP depends on the node limit */
#define DEFAULT_NUSEDSOLS     3              /* number of solutions that will be taken into account                   */
#define DEFAULT_NWAITINGNODES 200LL          /* number of nodes without incumbent change heuristic should wait        */
#define DEFAULT_RANDOMIZATION TRUE           /* should the choice which sols to take be randomized?                   */
#define DEFAULT_DONTWAITATROOT FALSE         /* should the nwaitingnodes parameter be ignored at the root node?       */
#define DEFAULT_USELPROWS     FALSE          /* should subproblem be created out of the rows in the LP rows,
                                              * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE           /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the
                                              * cutpool of the original scip be copied to constraints of the subscip
                                              */
#define DEFAULT_PERMUTE       FALSE          /* should the subproblem be permuted to increase diversification?        */
#define HASHSIZE_SOLS         500            /* size of hash table for solution tuples in crossover heuristic         */
#define DEFAULT_BESTSOLLIMIT   -1            /* limit on number of improving incumbent solutions in sub-CIP           */
#define DEFAULT_USEUCT        FALSE          /* should uct node selection be used at the beginning of the search?     */
#define DEFAULT_RANDSEED         7           /* initial random seed                                                   */

/* event handler properties */
#define EVENTHDLR_NAME         "Crossover"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"

/*
 * Data structures
 */

typedef struct SolTuple SOLTUPLE;

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             prevlastsol;        /**< worst solution taken into account during the previous run         */
   SCIP_SOL*             prevbestsol;        /**< best solution during the previous run                             */
   int                   prevnsols;          /**< number of all solutions during the previous run                   */

   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem               */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem               */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes        */
   SCIP_Longint          usednodes;          /**< nodes already used by crossover in earlier calls                  */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem     */

   int                   nusedsols;          /**< number of solutions that will be taken into account               */
   SCIP_Longint          nwaitingnodes;      /**< number of nodes without incumbent change heuristic should wait    */
   unsigned int          nfailures;          /**< number of failures since last successful call                     */
   SCIP_Longint          nextnodenumber;     /**< number of nodes at which crossover should be called the next time */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed     */
   SCIP_Real             minimprove;         /**< factor by which Crossover should at least improve the incumbent   */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   SCIP_Bool             randomization;      /**< should the choice which sols to take be randomized?               */
   SCIP_Bool             dontwaitatroot;     /**< should the nwaitingnodes parameter be ignored at the root node?   */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator                                           */
   SCIP_HASHTABLE*       hashtable;          /**< hashtable used to store the solution tuples already used          */
   SOLTUPLE*             lasttuple;          /**< last tuple of solutions created by crossover                      */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?      */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?                                     */
   SCIP_Bool             permute;            /**< should the subproblem be permuted to increase diversification?    */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP       */
   SCIP_Bool             useuct;             /**< should uct node selection be used at the beginning of the search?  */
};

/** n-tuple of solutions and their hashkey */
struct SolTuple
{
   int*                  indices;            /**< sorted array of solution indices                                 */
   int                   size;               /**< size of the array (should be heurdata->nusedsols)                */
   unsigned int          key;                /**< hashkey of the tuple                                             */
   SOLTUPLE*             prev;               /**< previous solution tuple created                                  */
};


/*
 * Local methods
 */

/** gets the hash key of a solution tuple */
static
SCIP_DECL_HASHGETKEY(hashGetKeySols)
{  /*lint --e{715}*/
   return elem;
}


/** returns TRUE iff both solution tuples are identical */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqSols)
{  /*lint --e{715}*/
   int i;
   int size;

   int* indices1;
   int* indices2;

   indices1 = ((SOLTUPLE*)key1)->indices;
   indices2 = ((SOLTUPLE*)key2)->indices;

   /* there should be two nonempty arrays of the same size */
   assert(indices1 != NULL);
   assert(indices2 != NULL);
   assert(((SOLTUPLE*)key1)->size == ((SOLTUPLE*)key2)->size);

   size = ((SOLTUPLE*)key1)->size;

   /* compare arrays by components, return TRUE, iff equal */
   for( i = 0; i < size; i++ )
   {
      if( indices1[i] != indices2[i] )
         return FALSE;
   }

   return TRUE;
}


/** returns hashkey of a solution tuple */
static
SCIP_DECL_HASHKEYVAL(hashKeyValSols)
{  /*lint --e{715}*/
   return ((SOLTUPLE*)key)->key;
}


/** calculates a hash key for a given tuple of solution indices */
static
unsigned int calculateHashKey(
   int*                  indices,            /**< indices of solutions */
   int                   size                /**< number of solutions */
   )
{
   int i;
   unsigned int hashkey;

   /* hashkey should be (x1+1) * (x2+1) * ... * (xn+1) + x1 + x2 + ... + xn */
   hashkey = 1;
   for( i = 0; i < size; i++ )
      hashkey *= indices[i] + 1;
   for( i = 0; i < size; i++ )
      hashkey += indices[i];

   return hashkey;
}


/** insertion sort for a small int array */
static void sortArray(
   int*                  a,                  /**< array to be sorted */
   int                   size                /**< size of array */
   )
{
   int i;
   int j;
   int tmp;

   /* simple insertion sort algorithm */
   for( i = 1; i < size; i++ )
   {
      tmp = a[i];
      j = i-1;
      while( j >= 0 && a[j] > tmp )
      {
         a[j+1] = a[j]; /*lint !e679*/
         j = j-1;
      }
      a[j+1] = tmp; /*lint !e679*/
   }
}


/** creates a new tuple of solutions */
static
SCIP_RETCODE createSolTuple(
   SCIP*                 scip,               /**< original SCIP data structure */
   SOLTUPLE**            elem,               /**< tuple of solutions which should be created */
   int*                  indices,            /**< indices of solutions */
   int                   size,               /**< number of solutions */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   /* memory allocation */
   SCIP_CALL( SCIPallocBlockMemory(scip, elem) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*elem)->indices, size) );
   BMScopyMemoryArray((*elem)->indices, indices, size);

   /* data input */
   sortArray(indices,size);
   (*elem)->size = size;
   (*elem)->key = calculateHashKey((*elem)->indices, (*elem)->size);
   (*elem)->prev = heurdata->lasttuple;

   /* update heurdata */
   heurdata->lasttuple = *elem;
   return SCIP_OKAY;
}


/** checks whether the new solution was found at the same node by the same heuristic as an already selected one */
static
SCIP_Bool solHasNewSource(
   SCIP_SOL**            sols,               /**< feasible SCIP solutions */
   int*                  selection,          /**< pool of solutions crossover uses */
   int                   selectionsize,      /**< size of solution pool */
   int                   newsol              /**< candidate solution */
   )
{
   int i;

   for( i = 0; i < selectionsize; i++ )
   {
      if( SCIPsolGetHeur(sols[selection[i]]) == SCIPsolGetHeur(sols[newsol])
         && SCIPsolGetNodenum(sols[selection[i]]) == SCIPsolGetNodenum(sols[newsol]) )
         return FALSE;
   }

   return TRUE;
}

/** randomly selects the solutions crossover will use from the pool of all solutions found so far */
static
SCIP_RETCODE selectSolsRandomized(
   SCIP*                 scip,               /**< original SCIP data structure */
   int*                  selection,          /**< pool of solutions crossover uses  */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            success             /**< pointer to store whether the process was successful */
   )
{
   int i;
   int j;
   int lastsol;          /* the worst solution possible to choose */
   int nusedsols;        /* number of solutions which will be chosen */

   SOLTUPLE* elem;
   SCIP_SOL** sols;

   assert( success != NULL );

   /* initialization */
   nusedsols = heurdata->nusedsols;
   lastsol = SCIPgetNSols(scip);
   sols = SCIPgetSols(scip);
   assert(nusedsols < lastsol);

   i = 0;
   *success = FALSE;

   /* perform at maximum 10 restarts and stop as soon as a new set of solutions is found */
   while( !*success && i < 10 )
   {
      SCIP_Bool validtuple;

      validtuple = TRUE;
      for( j = 0; j < nusedsols && validtuple; j++ )
      {
         int k;
         k = SCIPrandomGetInt(heurdata->randnumgen, nusedsols-j-1, lastsol-1);

         /* ensure that the solution does not have a similar source as the others */
         while( k >= nusedsols-j-1 && !solHasNewSource(sols, selection, j, k) )
            k--;

         validtuple = (k >= nusedsols-j-1);
         selection[j] = k;
         lastsol = k;
      }

      if( validtuple )
      {
         /* creates an object ready to be inserted into the hashtable */
         SCIP_CALL( createSolTuple(scip, &elem, selection, nusedsols, heurdata) );

         /* check whether the randomized set is already in the hashtable, if not, insert it */
         if( !SCIPhashtableExists(heurdata->hashtable, elem) )
         {
            SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
            *success = TRUE;
         }
      }
      i++;
   }

   return SCIP_OKAY;
}


/** determines the fixings for the CROSSOVER subproblem and checks whether enough fixings were found */
static
SCIP_RETCODE fixVariables(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR**            fixedvars,          /**< array to store source SCIP variables whose copies should be fixed in the sub-SCIP */
   SCIP_Real*            fixedvals,          /**< array to store solution values for variable fixing */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int                   fixedvarssize,      /**< size of the arrays to store fixing variables */
   int*                  selection,          /**< pool of solutions crossover will use */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                */
   SCIP_SOL** sols;                          /* pool of solutions                      */
   SCIP_Real fixingrate;                     /* percentage of variables that are fixed */

   int nvars;
   int nbinvars;
   int nintvars;

   int i;
   int j;

   sols = SCIPgetSols(scip);
   assert(sols != NULL);
   assert(fixedvars != NULL);
   assert(fixedvals != NULL);

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   assert(fixedvarssize >= nbinvars + nintvars);

   *nfixedvars = 0;

   /* go through discrete variables and collect potential fixings */
   for( i = 0; i < nbinvars + nintvars; i++ )
   {
      SCIP_Real solval;
      SCIP_Bool fixable;

      fixable = TRUE;
      solval = SCIPgetSolVal(scip, sols[selection[0]], vars[i]);

      /* check, whether variable's value is identical for each selected solution */
      for( j = 1; j < heurdata->nusedsols; j++ )
      {
         SCIP_Real varsolval;
         varsolval = SCIPgetSolVal(scip, sols[selection[j]], vars[i]);
         if( REALABS(solval - varsolval) > 0.5 )
         {
            fixable = FALSE;
            break;
         }
      }

      /* original solval can be outside transformed global bounds */
      fixable = fixable && SCIPvarGetLbGlobal(vars[i]) <= solval && solval <= SCIPvarGetUbGlobal(vars[i]);

      /* if solutions' values are equal, variable should be fixed in the subproblem */
      if( fixable )
      {
         fixedvars[(*nfixedvars)] = vars[i];
         fixedvals[(*nfixedvars)] = solval;
         (*nfixedvars)++;
      }
   }

   fixingrate = (SCIP_Real)(*nfixedvars) / (SCIP_Real)(MAX(nbinvars + nintvars, 1));

   /* if all variables would be fixed or amount of fixed variables is insufficient,
    * skip subproblem creation and abort immediately
    */
   *success = (*nfixedvars) < nbinvars + nintvars && fixingrate >= heurdata->minfixingrate;

   return SCIP_OKAY;
}

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE determineVariableFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR**            fixedvars,          /**< array to store source SCIP variables whose copies should be fixed in the sub-SCIP */
   SCIP_Real*            fixedvals,          /**< array to store solution values for variable fixing */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int                   fixedvarssize,      /**< size of the arrays to store fixing variables */
   int*                  selection,          /**< pool of solutions crossover will use */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_SOL** sols;                         /* array of all solutions found so far         */
   int nsols;                               /* number of all solutions found so far        */
   int nusedsols;                           /* number of solutions to use in crossover     */
   int i;

   /* get solutions' data */
   nsols = SCIPgetNSols(scip);
   sols = SCIPgetSols(scip);
   nusedsols = heurdata->nusedsols;

   assert(nusedsols > 1);
   assert(nsols >= nusedsols);

   /* use nusedsols best solutions if randomization is deactivated or there are only nusedsols solutions at hand
    * or a good new solution was found since last call */
   if( !heurdata->randomization || nsols == nusedsols || heurdata->prevlastsol != sols[nusedsols-1] )
   {
      SOLTUPLE* elem;
      SCIP_HEUR* solheur;
      SCIP_Longint solnodenum;
      SCIP_Bool allsame;

      for( i = 0; i < nusedsols; i++ )
         selection[i] = i;
      SCIP_CALL( createSolTuple(scip, &elem, selection, nusedsols, heurdata) );

      solheur = SCIPsolGetHeur(sols[0]);
      solnodenum = SCIPsolGetNodenum(sols[0]);
      allsame = TRUE;

      /* check, whether all solutions have been found by the same heuristic at the same node; in this case we do not run
       * crossover, since it would probably just optimize over the same space as the other heuristic
       */
      for( i = 1; i < nusedsols; i++ )
      {
         if( SCIPsolGetHeur(sols[i]) != solheur || SCIPsolGetNodenum(sols[i]) != solnodenum )
            allsame = FALSE;
      }
      *success = !allsame && !SCIPhashtableExists(heurdata->hashtable, elem);

      /* check, whether solution tuple has already been tried */
      if( !SCIPhashtableExists(heurdata->hashtable, elem) )
      {
         SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
      }

      /* if solution tuple has already been tried, randomization is allowed and enough solutions are at hand, try
       * to randomize another tuple. E.g., this can happen if the last crossover solution was among the best ones */
      if( !(*success) && heurdata->randomization && nsols > nusedsols )
      {
         SCIP_CALL( selectSolsRandomized(scip, selection, heurdata, success) );
      }

   }
   /* otherwise randomize the set of solutions */
   else
   {
      SCIP_CALL( selectSolsRandomized(scip, selection, heurdata, success) );
   }

   /* no acceptable solution tuple could be created */
   if( !(*success) )
      return SCIP_OKAY;

   /* set up the variables of the subproblem */
   SCIP_CALL( fixVariables(scip, fixedvars, fixedvals, nfixedvars, fixedvarssize, selection, heurdata, success) );

   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEUR*            heur,               /**< crossover heuristic structure */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   int*                  solindex,           /**< index of the solution */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );
   *solindex = SCIPsolGetIndex(newsol);

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** updates heurdata after a run of crossover */
static
void updateFailureStatistic(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   /* increase number of failures, calculate next node at which crossover should be called and update actual solutions */
   heurdata->nfailures++;
   heurdata->nextnodenumber = (heurdata->nfailures <= 25
      ? SCIPgetNNodes(scip) + 100*(2LL << heurdata->nfailures) /*lint !e703*/
      : SCIP_LONGINT_MAX);
}

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecCrossover)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetNLPs(scip) > heurdata->lplimfac * heurdata->nodelimit )
   {
      SCIPdebugMsg(scip, "interrupt after  %" SCIP_LONGINT_FORMAT " LPs\n", SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyCrossover)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurCrossover(scip) );

   return SCIP_OKAY;
}

/** setup and solve the subproblem and catch the return code */
static
SCIP_RETCODE setupAndSolveSubscipCrossover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_HEUR*            heur,               /**< mutation heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristics data */
   SCIP_VAR**            vars,               /**< SCIP variables */
   SCIP_VAR**            fixedvars,          /**< array to store the variables that should be fixed in the subproblem */
   SCIP_Real*            fixedvals,          /**< array to store the fixing values to fix variables in the subproblem */
   SCIP_Longint          nstallnodes,        /**< node limit for the subproblem */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  selection,          /**< pool of solutions crossover uses */
   int                   nvars,              /**< number of original problem's variables */
   int                   nfixedvars,         /**< the number of variables that should be fixed */
   int                   nusedsols           /**< number of solutions which will be chosen */
   )
{
   SCIP_EVENTHDLR* eventhdlr;                /* event handler for LP events                     */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR** subvars;                       /* subproblem's variables                          */
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem             */
   SCIP_Real upperbound;
   SCIP_Bool success;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
   success = FALSE;

   /* create a copy of the transformed problem to be used by the heuristic */
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "crossover", fixedvars, fixedvals, nfixedvars,
         heurdata->uselprows, heurdata->copycuts, &success, NULL) );

   eventhdlr = NULL;
   /* create event handler for LP events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecCrossover, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* store copied variables in the order in which they appear in the main SCIP */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   for( i = 0; i < nvars; i++ )
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsollimit) );

   /* set limits for the subproblem */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   heurdata->nodelimit = nstallnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nstallnodes) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* activate uct node selection at the top of the tree */
   if( heurdata->useuct && SCIPfindNodesel(subscip, "uct") != NULL && !SCIPisParamFixed(subscip, "nodeselection/uct/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/uct/stdpriority", INT_MAX/2) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* enable conflict analysis, disable analysis of boundexceeding LPs, and restrict conflict pool */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetCharParam(subscip, "conflict/useboundlp", 'o') );
   }
   if( !SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 500) );
   }

   /* add an objective cutoff */
   assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));

   upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
   if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
   {
      cutoff = (1-heurdata->minimprove)*SCIPgetUpperbound(scip) + heurdata->minimprove*SCIPgetLowerbound(scip);
   }
   else
   {
      if( SCIPgetUpperbound ( scip ) >= 0 )
         cutoff = ( 1 - heurdata->minimprove ) * SCIPgetUpperbound ( scip );
      else
         cutoff = ( 1 + heurdata->minimprove ) * SCIPgetUpperbound ( scip );
   }
   cutoff = MIN(upperbound, cutoff );
   SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );

   /* permute the subproblem to increase diversification */
   if( heurdata->permute )
   {
      SCIP_CALL( SCIPpermuteProb(subscip, SCIPinitializeRandomSeed(scip, (int) SCIPheurGetNCalls(heur)), TRUE, TRUE, TRUE, TRUE, TRUE) );
   }

   /* catch LP events of sub-SCIP */
   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   /* this code can be enabled whenever the subproblem should be written out */
#ifdef SCIP_DISABLED_CODE
   SCIPdebug( SCIP_CALL( SCIPwriteOrigProblem(subscip, "crossoverprob.cip", "cip", FALSE) ) );
#endif
   /* solve the subproblem */
   SCIPdebugMsg(scip, "Solve Crossover subMIP\n");

   /* Errors in solving the subproblem should not kill the overall solving process.
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop. */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* drop LP events of sub-SCIP */
   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   heurdata->usednodes += SCIPgetNNodes(subscip);

   /* merge histories of the subscip-variables to the SCIP variables. */
   SCIP_CALL( SCIPmergeVariableStatistics(subscip, scip, subvars, vars, nvars) );
   SCIPdebugMsg(scip, "Transferring variable histories complete\n");

   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;
      int solindex;                             /* index of the solution created by crossover          */

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      solindex = -1;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &solindex, &success) );
      }

      if( success )
      {
         int tmp;

         assert(solindex != -1);

         *result = SCIP_FOUNDSOL;

         /* insert all crossings of the new solution and (nusedsols-1) of its parents into the hashtable
          * in order to avoid incest ;)
          */
         for( i = 0; i < nusedsols; i++ )
         {
            SOLTUPLE* elem;
            tmp = selection[i];
            selection[i] = solindex;

            SCIP_CALL( createSolTuple(scip, &elem, selection, nusedsols, heurdata) );
            SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
            selection[i] = tmp;
         }

         /* if solution was among the best ones, crossover should not be called until another good solution was found */
         if( !heurdata->randomization )
         {
            heurdata->prevbestsol = SCIPgetBestSol(scip);
            heurdata->prevlastsol = SCIPgetSols(scip)[heurdata->nusedsols-1];
         }
      }

      /* if solution is not better than incumbent or could not be added to problem => run is counted as a failure */
      if( !success || solindex != SCIPsolGetIndex(SCIPgetBestSol(scip)) )
         updateFailureStatistic(scip, heurdata);
   }
   else
   {
      /* if no new solution was found, run was a failure */
      updateFailureStatistic(scip, heurdata);
   }

   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeCrossover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitCrossover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->prevlastsol = NULL;
   heurdata->prevbestsol = NULL;
   heurdata->lasttuple = NULL;
   heurdata->nfailures = 0;
   heurdata->prevnsols = 0;
   heurdata->nextnodenumber = 0;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, DEFAULT_RANDSEED) );

   /* initialize hash table */
   SCIP_CALL( SCIPhashtableCreate(&heurdata->hashtable, SCIPblkmem(scip), HASHSIZE_SOLS,
         hashGetKeySols, hashKeyEqSols, hashKeyValSols, NULL) );
   assert(heurdata->hashtable != NULL);

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitCrossover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SOLTUPLE* soltuple;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   soltuple = heurdata->lasttuple;

   /* free all soltuples iteratively */
   while( soltuple != NULL )
   {
      SOLTUPLE* tmp;
      tmp = soltuple->prev;
      SCIPfreeBlockMemoryArray(scip, &soltuple->indices, soltuple->size);
      SCIPfreeBlockMemory(scip, &soltuple);
      soltuple = tmp;
   }

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   /* free hash table */
   assert(heurdata->hashtable != NULL);
   SCIPhashtableFree(&heurdata->hashtable);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecCrossover)
{  /*lint --e{715}*/
   SCIP* subscip;                            /* the subproblem created by crossover                 */
   SCIP_HEURDATA* heurdata;                  /* primal heuristic data                               */
   SCIP_VAR** vars;                          /* original problem's variables                        */
   SCIP_VAR** fixedvars;
   SCIP_SOL** sols;
   SCIP_RETCODE retcode;
   SCIP_Longint nstallnodes;                 /* node limit for the subproblem                       */
   SCIP_Bool success;
   SCIP_Real* fixedvals;
   int* selection;                           /* pool of solutions crossover uses                    */
   int nvars;                                /* number of original problem's variables              */
   int nbinvars;
   int nintvars;
   int nusedsols;
   int nfixedvars;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   nusedsols = heurdata->nusedsols;

   *result = SCIP_DELAYED;

   /* only call heuristic, if enough solutions are at hand */
   if( SCIPgetNSols(scip) < nusedsols  )
      return SCIP_OKAY;

   sols = SCIPgetSols(scip);
   assert(sols != NULL);

   /* if one good solution was found, heuristic should not be delayed any longer */
   if( sols[nusedsols-1] != heurdata->prevlastsol )
   {
      heurdata->nextnodenumber = SCIPgetNNodes(scip);
      if( sols[0] != heurdata->prevbestsol )
         heurdata->nfailures = 0;
   }
   /* in nonrandomized mode: only recall heuristic, if at least one new good solution was found in the meantime */
   else if( !heurdata->randomization )
      return SCIP_OKAY;

   /* if heuristic should be delayed, wait until certain number of nodes is reached */
   if( SCIPgetNNodes(scip) < heurdata->nextnodenumber )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip, SCIPgetBestSol(scip)) < heurdata->nwaitingnodes
      && (SCIPgetDepth(scip) > 0 || !heurdata->dontwaitatroot) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward Crossover if it succeeded often */
   nstallnodes = (SCIP_Longint)
      (nstallnodes * (1.0 + 2.0*(SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur)+1.0)));

   /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
      return SCIP_OKAY;

   /* consider time and memory limits of the main SCIP */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if( !success )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
     return SCIP_OKAY;

   /* get variable information */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   assert(nvars > 0);

   /* check whether discrete variables are available */
   if( nbinvars == 0 && nintvars == 0 )
      return SCIP_OKAY;

   /* allocate necessary buffer storage for selection of variable fixings */
   SCIP_CALL( SCIPallocBufferArray(scip, &selection, nusedsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, nbinvars + nintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, nbinvars + nintvars) );

   success = FALSE;
   nfixedvars = 0;
   /* determine fixings of variables with same value in a certain set of solutions */
   SCIP_CALL( determineVariableFixings(scip, fixedvars, fixedvals, &nfixedvars, nbinvars + nintvars, selection, heurdata, &success) );

   heurdata->prevbestsol = SCIPgetBestSol(scip);
   heurdata->prevlastsol = sols[heurdata->nusedsols-1];

   /* if creation of sub-SCIP was aborted (e.g. due to number of fixings), free sub-SCIP and abort */
   if( !success )
   {
      /* this run will be counted as a failure since no new solution tuple could be generated or the neighborhood of the
       * solution was not fruitful in the sense that it was too big
       */
      updateFailureStatistic(scip, heurdata);

      goto TERMINATE;
   }

   *result = SCIP_DIDNOTFIND;
   /* initializing the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* setup and solve the subproblem and catch the return code */
   retcode = setupAndSolveSubscipCrossover(scip, subscip, heur, heurdata, vars,
               fixedvars, fixedvals, nstallnodes, result, selection, nvars, nfixedvars, nusedsols);

   /* free the subscip in any case */
   SCIP_CALL( SCIPfree(&subscip) );
   SCIP_CALL( retcode );

TERMINATE:
   /* free buffer storage for variable fixings */
   SCIPfreeBufferArray(scip, &fixedvals);
   SCIPfreeBufferArray(scip, &fixedvars);
   SCIPfreeBufferArray(scip, &selection);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the crossover primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCrossover(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Crossover primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecCrossover, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyCrossover) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeCrossover) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitCrossover) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitCrossover) );

   /* add crossover primal heuristic parameters */

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nusedsols",
         "number of solutions to be taken into account",
         &heurdata->nusedsols, FALSE, DEFAULT_NUSEDSOLS, 2, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "minimum percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which Crossover should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/randomization",
         "should the choice which sols to take be randomized?",
         &heurdata->randomization, TRUE, DEFAULT_RANDOMIZATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/dontwaitatroot",
         "should the nwaitingnodes parameter be ignored at the root node?",
         &heurdata->dontwaitatroot, TRUE, DEFAULT_DONTWAITATROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/permute",
         "should the subproblem be permuted to increase diversification?",
         &heurdata->permute, TRUE, DEFAULT_PERMUTE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/bestsollimit",
         "limit on number of improving incumbent solutions in sub-CIP",
         &heurdata->bestsollimit, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useuct",
         "should uct node selection be used at the beginning of the search?",
         &heurdata->useuct, TRUE, DEFAULT_USEUCT, NULL, NULL) );
   return SCIP_OKAY;
}
