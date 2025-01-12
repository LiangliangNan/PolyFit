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

/**@file   prop_symmetry.c
 * @ingroup DEFPLUGINS_PROP
 * @brief  propagator for handling symmetries
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 * @author Fabian Wegscheider
 *
 * This propagator combines the following symmetry handling functionalities:
 * - It allows to compute symmetries of the problem and to store this information in adequate form. The symmetry
 *   information can be accessed through external functions.
 * - It allows to add the following symmetry breaking constraints:
 *    - symresack constraints, which separate minimal cover inequalities
 *    - orbitope constraints, if special symmetry group structures are detected
 * - It allows to apply orbital fixing.
 *
 *
 * @section SYMCOMP Symmetry Computation
 *
 * The following comments apply to symmetry computation.
 *
 * - The generic functionality of the compute_symmetry.h interface is used.
 * - We treat implicit integer variables as if they were continuous/real variables. The reason is that there is currently
 *   no distinction between implicit integer and implicit binary. Moreover, currently implicit integer variables hurt
 *   our code more than continuous/real variables (we basically do not handle integral variables at all).
 * - We do not copy symmetry information, since it is not clear how this information transfers. Moreover, copying
 *   symmetry might inhibit heuristics. But note that solving a sub-SCIP might then happen without symmetry information!
 *
 *
 * @section SYMBREAK Symmetry Handling Constraints
 *
 * The following comments apply to adding symmetry handling constraints.
 *
 * - The code automatically detects whether symmetry substructures like symresacks or orbitopes are present and possibly
 *   adds the corresponding constraints.
 * - If orbital fixing is active, only orbitopes are added (if present) and no symresacks.
 * - We try to compute symmetry as late as possible and then add constraints based on this information.
 * - Currently, we only allocate memory for pointers to symresack constraints for group generators. If further
 *   constraints are considered, we have to reallocate memory.
 *
 *
 * @section OF Orbital Fixing
 *
 * Orbital fixing is implemented as introduced by@n
 * F. Margot: Exploiting orbits in symmetric ILP. Math. Program., 98(1-3):3–21, 2003.
 *
 * The method computes orbits of variables with respect to the subgroup of the symmetry group that stabilizes the
 * variables globally fixed or branched to 1. Then one can fix all variables in an orbit to 0 or 1 if one of the other
 * variables in the orbit is fixed to 0 or 1, respectively. Different from Margot, the subgroup is obtained by filtering
 * out generators that do not individually stabilize the variables branched to 1.
 *
 * @pre All variable fixings applied by other components are required to be strict, i.e., if one variable is fixed to
 *      a certain value v, all other variables in the same variable orbit can be fixed to v as well, c.f.@n
 *      F. Margot: Symmetry in integer linear programming. 50 Years of Integer Programming, 647-686, Springer 2010.
 *
 * To illustrate this, consider the example \f$\max\{x_1 + x_2 : x_1 + x_2 \leq 1, Ay \leq b,
 * (x,y) \in \{0,1\}^{2 + n}\} \f$. Since \f$x_1\f$ and \f$x_2\f$ are independent from the remaining problem, the
 * setppc constraint handler may fix \f$(x_1,x_2) = (1,0)\f$. However, since both variables are symmetric, this setting
 * is not strict (if it was strict, both variables would have been set to the same value) and orbital fixing would
 * declare this subsolution as infeasible (there exists an orbit of non-branching variables that are fixed to different
 * values). To avoid this situation, we have to assume that all non-strict settings fix variables globally, i.e., we
 * can take care of it by taking variables into account that have been globally fixed to 1. In fact, it suffices to
 * consider one kind of global fixings since stabilizing one kind prevents an orbit to contain variables that have
 * been fixed globally to different values.
 *
 * @pre All non-strict settings are global settings, since otherwise, we cannot (efficiently) take care of them.
 *
 * @pre No non-strict setting algorithm is interrupted early (e.g., by a time or iteration limit), since this may lead to
 * wrong decisions by orbital fixing as well. For example, if cons_setppc in the above toy example starts by fixing
 * \f$x_2 = 0\f$ and is interrupted afterwards, orbital fixing detects that the orbit \f$\{x_1, x_2\}\f$ contains
 * one variable that is fixed to 0, and thus, it fixes \f$x_1\f$ to 0 as well. Thus, after these reductions, every
 * feasible solution has objective 0 which is not optimal. This situation would not occur if the non-strict setting is
 * complete, because then \f$x_1\f$ is globally fixed to 1, and thus, is stabilized in orbital fixing.
 *
 * Note that orbital fixing might lead to wrong results if it is called in repropagation of a node, because the path
 * from the node to the root might have been changed. Thus, the stabilizers of global 1-fixing and 1-branchings of the
 * initial propagation and repropagation might differ, which may cause conflicts. For this reason, orbital fixing cannot
 * be called in repropagation.
 *
 * @note If, besides orbital fixing, also symmetry handling constraints shall be added, orbital fixing is only applied
 *       to symmetry components that are not handled by orbitope constraints.
 *
 *
 * @section SST Cuts derived from the Schreier Sims table
 *
 * SST cuts have been introduced by@n
 * D. Salvagnin: Symmetry Breaking Inequalities from the Schreier-Sims table. CPAIOR 2018 Proceedings, 521-529, 2018.
 *
 * These inequalities are computed as follows. Throughout these procedure a set of so-called leaders is maintained.
 * Initially the set of leaders is empty. In a first step, select a variable \f$x_i\f$ and compute its orbit w.r.t.
 * the symmetry group of the mixed-integer program. For each variable \f$x_j\f$ in the orbit of \f$x_i\f$, the
 * inequality \f$x_i \geq x_j\f$ is a valid symmetry handling inequality, which can be added to the mixed-integer
 * program. We call \f$x_i\f$ the leader of this inequality. Add the leader \f$x_i\f$ to the set of leaders and
 * compute the pointwise stabilizer of the leader set. In the next step, select a new variable, compute its orbit
 * w.r.t. the stabilizer group of the leaders, add the inequalities based on this orbit, and add the new leader
 * to the set of leaders. This procedure is iterated until the pointwise stabilizer group of the leaders has become
 * trivial.
 *
 * @todo Possibly turn off propagator in subtrees.
 * @todo Check application of conflict resolution.
 * @todo Check whether one should switch the role of 0 and 1
 * @todo Implement stablizer computation?
 * @todo Implement isomorphism pruning?
 * @todo Implement particular preprocessing rules
 * @todo Separate permuted cuts (first experiments not successful)
 * @todo Allow the computation of local symmetries
 * @todo Order rows of orbitopes (in particular packing/partitioning) w.r.t. cliques in conflict graph.
 */
/* #define SCIP_OUTPUT */
/* #define SCIP_OUTPUT_COMPONENT */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/cons_linear.h>
#include <scip/cons_knapsack.h>
#include <scip/cons_varbound.h>
#include <scip/cons_setppc.h>
#include <scip/cons_and.h>
#include <scip/cons_logicor.h>
#include <scip/cons_or.h>
#include <scip/cons_orbitope.h>
#include <scip/cons_symresack.h>
#include <scip/cons_xor.h>
#include <scip/cons_linking.h>
#include <scip/cons_bounddisjunction.h>
#include <scip/cons_nonlinear.h>
#include <scip/pub_expr.h>
#include <scip/misc.h>
#include <scip/scip_datastructures.h>

#include <scip/prop_symmetry.h>
#include <symmetry/compute_symmetry.h>
#include <scip/symmetry.h>

#include <string.h>

/* propagator properties */
#define PROP_NAME            "symmetry"
#define PROP_DESC            "propagator for handling symmetry"
#define PROP_TIMING    SCIP_PROPTIMING_BEFORELP   /**< propagation timing mask */
#define PROP_PRIORITY          -1000000           /**< propagator priority */
#define PROP_FREQ                     1           /**< propagator frequency */
#define PROP_DELAY                FALSE           /**< should propagation method be delayed, if other propagators found reductions? */

#define PROP_PRESOL_PRIORITY  -10000000           /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers) */
#define PROP_PRESOLTIMING   SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1           /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */


/* default parameter values for symmetry computation */
#define DEFAULT_MAXGENERATORS        1500    /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
#define DEFAULT_CHECKSYMMETRIES     FALSE    /**< Should all symmetries be checked after computation? */
#define DEFAULT_DISPLAYNORBITVARS   FALSE    /**< Should the number of variables affected by some symmetry be displayed? */
#define DEFAULT_USECOLUMNSPARSITY   FALSE    /**< Should the number of conss a variable is contained in be exploited in symmetry detection? */
#define DEFAULT_DOUBLEEQUATIONS     FALSE    /**< Double equations to positive/negative version? */
#define DEFAULT_COMPRESSSYMMETRIES   TRUE    /**< Should non-affected variables be removed from permutation to save memory? */
#define DEFAULT_COMPRESSTHRESHOLD     0.5    /**< Compression is used if percentage of moved vars is at most the threshold. */
#define DEFAULT_SYMFIXNONBINARYVARS FALSE    /**< Whether all non-binary variables shall be not affected by symmetries if OF is active? */
#define DEFAULT_ONLYBINARYSYMMETRY   TRUE    /**< Is only symmetry on binary variables used? */

/* default parameters for linear symmetry constraints */
#define DEFAULT_CONSSADDLP           TRUE    /**< Should the symmetry breaking constraints be added to the LP? */
#define DEFAULT_ADDSYMRESACKS       FALSE    /**< Add inequalities for symresacks for each generator? */
#define DEFAULT_DETECTORBITOPES      TRUE    /**< Should we check whether the components of the symmetry group can be handled by orbitopes? */
#define DEFAULT_DETECTSUBGROUPS      TRUE    /**< Should we try to detect orbitopes in subgroups of the symmetry group? */
#define DEFAULT_ADDWEAKSBCS          TRUE    /**< Should we add weak SBCs for enclosing orbit of symmetric subgroups? */
#define DEFAULT_ADDSTRONGSBCS       FALSE    /**< Should we add strong SBCs for enclosing orbit of symmetric subgroups if orbitopes are not used? */
#define DEFAULT_ADDCONSSTIMING          2    /**< timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving) */
#define DEFAULT_MAXNCONSSSUBGROUP  500000    /**< Maximum number of constraints up to which subgroup structures are detected */
#define DEFAULT_USEDYNAMICPROP       TRUE    /**< whether dynamic propagation should be used for full orbitopes */
#define DEFAULT_PREFERLESSROWS       TRUE    /**< Shall orbitopes with less rows be preferred in detection? */

/* default parameters for orbital fixing */
#define DEFAULT_OFSYMCOMPTIMING         2    /**< timing of symmetry computation for orbital fixing (0 = before presolving, 1 = during presolving, 2 = at first call) */
#define DEFAULT_PERFORMPRESOLVING   FALSE    /**< Run orbital fixing during presolving? */
#define DEFAULT_RECOMPUTERESTART        0    /**< Recompute symmetries after a restart has occurred? (0 = never, 1 = always, 2 = if OF found reduction) */

/* default parameters for Schreier Sims constraints */
#define DEFAULT_SSTTIEBREAKRULE   1          /**< index of tie break rule for selecting orbit for Schreier Sims constraints? */
#define DEFAULT_SSTLEADERRULE     0          /**< index of rule for selecting leader variables for Schreier Sims constraints? */
#define DEFAULT_SSTLEADERVARTYPE 14          /**< bitset encoding which variable types can be leaders (1: binary; 2: integer; 4: impl. int; 8: continuous);
                                              *   if multiple types are allowed, take the one with most affected vars */
#define DEFAULT_ADDCONFLICTCUTS       TRUE   /**< Should Schreier Sims constraints be added if we use a conflict based rule? */
#define DEFAULT_SSTADDCUTS            TRUE   /**< Should Schreier Sims constraints be added? */
#define DEFAULT_SSTMIXEDCOMPONENTS    TRUE   /**< Should Schreier Sims constraints be added if a symmetry component contains variables of different types? */

/* event handler properties */
#define EVENTHDLR_SYMMETRY_NAME    "symmetry"
#define EVENTHDLR_SYMMETRY_DESC    "filter global variable fixing event handler for orbital fixing"

/* output table properties */
#define TABLE_NAME_ORBITALFIXING        "orbitalfixing"
#define TABLE_DESC_ORBITALFIXING        "orbital fixing statistics"
#define TABLE_POSITION_ORBITALFIXING    7001                    /**< the position of the statistics table */
#define TABLE_EARLIEST_ORBITALFIXING    SCIP_STAGE_SOLVING      /**< output of the statistics table is only printed from this stage onwards */


/* other defines */
#define MAXGENNUMERATOR          64000000    /**< determine maximal number of generators by dividing this number by the number of variables */
#define SCIP_SPECIALVAL 1.12345678912345e+19 /**< special floating point value for handling zeros in bound disjunctions */
#define COMPRESSNVARSLB             25000    /**< lower bound on the number of variables above which compression could be performed */

/* macros for getting activeness of symmetry handling methods */
#define ISSYMRETOPESACTIVE(x)      (((unsigned) x & SYM_HANDLETYPE_SYMBREAK) != 0)
#define ISORBITALFIXINGACTIVE(x)   (((unsigned) x & SYM_HANDLETYPE_ORBITALFIXING) != 0)
#define ISSSTACTIVE(x)             (((unsigned) x & SYM_HANDLETYPE_SST) != 0)

#define ISSSTBINACTIVE(x)          (((unsigned) x & SCIP_SSTTYPE_BINARY) != 0)
#define ISSSTINTACTIVE(x)          (((unsigned) x & SCIP_SSTTYPE_INTEGER) != 0)
#define ISSSTIMPLINTACTIVE(x)      (((unsigned) x & SCIP_SSTTYPE_IMPLINT) != 0)
#define ISSSTCONTACTIVE(x)         (((unsigned) x & SCIP_SSTTYPE_CONTINUOUS) != 0)


/** propagator data */
struct SCIP_PropData
{
   /* symmetry group information */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nbinpermvars;       /**< number of binary variables for permuations */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
#ifndef NDEBUG
   SCIP_Real*            permvarsobj;        /**< objective values of permuted variables (for debugging) */
#endif
   int                   nperms;             /**< number of permutations */
   int                   nmaxperms;          /**< maximal number of permutations (needed for freeing storage) */
   int**                 perms;              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   int**                 permstrans;         /**< pointer to store transposed permutation generators as (npermvars x nperms) matrix */
   SCIP_HASHMAP*         permvarmap;         /**< map of variables to indices in permvars array */
   int                   nmovedpermvars;     /**< number of variables moved by any permutation */
   int                   nmovedbinpermvars;  /**< number of binary variables moved by any permutation */
   int                   nmovedintpermvars;  /**< number of integer variables moved by any permutation */
   int                   nmovedimplintpermvars; /**< number of implicitly integer variables moved by any permutation */
   int                   nmovedcontpermvars; /**< number of continuous variables moved by any permutation */
   SCIP_Shortbool*       nonbinpermvarcaptured; /**< array to store which non-binary variables have been captured
                                                 *   (only necessary for SST cuts) */

   /* components of symmetry group */
   int                   ncomponents;        /**< number of components of symmetry group */
   int                   ncompblocked;       /**< number of components that have been blocked */
   int*                  components;         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins;    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int*                  vartocomponent;     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   unsigned*             componentblocked;   /**< array to store which symmetry methods have been applied to a component using
                                              *   the same bitset as for misc/usesymmetry */

   /* further symmetry information */
   int                   nmovedvars;         /**< number of variables moved by some permutation */
   SCIP_Real             log10groupsize;     /**< log10 of size of symmetry group */
   SCIP_Bool             binvaraffected;     /**< whether binary variables are affected by some symmetry */

   /* for symmetry computation */
   int                   maxgenerators;      /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
   SCIP_Bool             checksymmetries;    /**< Should all symmetries be checked after computation? */
   SCIP_Bool             displaynorbitvars;  /**< Whether the number of variables in non-trivial orbits shall be computed */
   SCIP_Bool             compresssymmetries; /**< Should non-affected variables be removed from permutation to save memory? */
   SCIP_Real             compressthreshold;  /**< Compression is used if percentage of moved vars is at most the threshold. */
   SCIP_Bool             compressed;         /**< Whether symmetry data has been compressed */
   SCIP_Bool             computedsymmetry;   /**< Have we already tried to compute symmetries? */
   int                   usesymmetry;        /**< encoding of active symmetry handling methods (for debugging) */
   SCIP_Bool             usecolumnsparsity;  /**< Should the number of conss a variable is contained in be exploited in symmetry detection? */
   SCIP_Bool             doubleequations;    /**< Double equations to positive/negative version? */
   SCIP_Bool             symfixnonbinaryvars; /**< Whether all non-binary variables shall be not affected by symmetries if OF is active? */
   SCIP_Bool             onlybinarysymmetry; /**< Whether only symmetry on binary variables is used */

   /* for symmetry constraints */
   SCIP_Bool             symconsenabled;     /**< Should symmetry constraints be added? */
   SCIP_Bool             triedaddconss;      /**< whether we already tried to add symmetry breaking constraints */
   SCIP_Bool             conssaddlp;         /**< Should the symmetry breaking constraints be added to the LP? */
   SCIP_Bool             addsymresacks;      /**< Add symresack constraints for each generator? */
   int                   addconsstiming;     /**< timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving) */
   SCIP_CONS**           genorbconss;        /**< list of generated orbitope/orbisack/symresack constraints */
   SCIP_CONS**           genlinconss;        /**< list of generated linear constraints */
   int                   ngenorbconss;       /**< number of generated orbitope/orbisack/symresack constraints */
   int                   ngenlinconss;       /**< number of generated linear constraints */
   int                   genlinconsssize;    /**< size of linear constraints array */
   int                   nsymresacks;        /**< number of symresack constraints */
   SCIP_Bool             detectorbitopes;    /**< Should we check whether the components of the symmetry group can be handled by orbitopes? */
   SCIP_Bool             detectsubgroups;    /**< Should we try to detect orbitopes in subgroups of the symmetry group? */
   SCIP_Bool             addweaksbcs;        /**< Should we add weak SBCs for enclosing orbit of symmetric subgroups? */
   SCIP_Bool             addstrongsbcs;      /**< Should we add strong SBCs for enclosing orbit of symmetric subgroups if orbitopes are not used? */
   int                   norbitopes;         /**< number of orbitope constraints */
   SCIP_Bool*            isnonlinvar;        /**< array indicating whether variables apper non-linearly */
   SCIP_CONSHDLR*        conshdlr_nonlinear; /**< nonlinear constraint handler */
   int                   maxnconsssubgroup;  /**< maximum number of constraints up to which subgroup structures are detected */
   SCIP_Bool             usedynamicprop;     /**< whether dynamic propagation should be used for full orbitopes */
   SCIP_Bool             preferlessrows;     /**< Shall orbitopes with less rows be preferred in detection? */

   /* data necessary for orbital fixing */
   SCIP_Bool             ofenabled;          /**< Run orbital fixing? */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for handling global variable fixings */
   SCIP_Shortbool*       bg0;                /**< bitset to store variables globally fixed to 0 */
   int*                  bg0list;            /**< list of variables globally fixed to 0 */
   int                   nbg0;               /**< number of variables in bg0 and bg0list */
   SCIP_Shortbool*       bg1;                /**< bitset to store variables globally fixed or branched to 1 */
   int*                  bg1list;            /**< list of variables globally fixed or branched to 1 */
   int                   nbg1;               /**< number of variables in bg1 and bg1list */
   int*                  permvarsevents;     /**< stores events caught for permvars */
   SCIP_Shortbool*       inactiveperms;      /**< array to store whether permutations are inactive */
   SCIP_Bool             performpresolving;  /**< Run orbital fixing during presolving? */
   int                   recomputerestart;   /**< Recompute symmetries after a restart has occured? (0 = never, 1 = always, 2 = if OF found reduction) */
   int                   ofsymcomptiming;    /**< timing of orbital fixing (0 = before presolving, 1 = during presolving, 2 = at first call) */
   int                   lastrestart;        /**< last restart for which symmetries have been computed */
   int                   nfixedzero;         /**< number of variables fixed to 0 */
   int                   nfixedone;          /**< number of variables fixed to 1 */
   SCIP_Longint          nodenumber;         /**< number of node where propagation has been last applied */
   SCIP_Bool             offoundreduction;   /**< whether orbital fixing has found a reduction since the last time computing symmetries */

   /* data necessary for Schreier Sims constraints */
   SCIP_Bool             sstenabled;         /**< Use Schreier Sims constraints? */
   SCIP_CONS**           sstconss;           /**< list of generated schreier sims conss */
   int                   nsstconss;          /**< number of generated schreier sims conss */
   int                   maxnsstconss;       /**< maximum number of conss in sstconss */
   int                   sstleaderrule;      /**< rule to select leader  */
   int                   ssttiebreakrule;    /**< tie break rule for leader selection */
   int                   sstleadervartype;   /**< bitset encoding which variable types can be leaders;
                                              *   if multiple types are allowed, take the one with most affected vars */
   int*                  leaders;            /**< index of orbit leaders in permvars */
   int                   nleaders;           /**< number of orbit leaders in leaders array */
   int                   maxnleaders;        /**< maximum number of leaders in leaders array */
   SCIP_Bool             addconflictcuts;    /**< Should Schreier Sims constraints be added if we use a conflict based rule? */
   SCIP_Bool             sstaddcuts;         /**< Should Schreier Sims constraints be added? */
   SCIP_Bool             sstmixedcomponents; /**< Should Schreier Sims constraints be added if a symmetry component contains variables of different types? */
};

/** node data of a given node in the conflict graph */
struct SCIP_NodeData
{
   SCIP_VAR*             var;                /**< variable belonging to node */
   int                   orbitidx;           /**< orbit of variable w.r.t. current stabilizer subgroup
                                              *   or -1 if not affected by symmetry */
   int                   nconflictinorbit;   /**< number of variables the node's var is in conflict with */
   int                   orbitsize;          /**< size of the variable's orbit */
   int                   posinorbit;         /**< position of variable in its orbit */
   SCIP_Bool             active;             /**< whether variable has not been fixed by Schreier Sims code */
};
typedef struct SCIP_NodeData SCIP_NODEDATA;


/*
 * Event handler callback methods
 */

/** exec the event handler for handling global variable bound changes (necessary for orbital fixing)
 *
 *  Global variable fixings during the solving process might arise because parts of the tree are pruned or if certain
 *  preprocessing steps are performed that do not correspond to strict setting algorithms. Since these fixings might be
 *  caused by or be in conflict with orbital fixing, they can be in conflict with the symmetry handling decisions of
 *  orbital fixing in the part of the tree that is not pruned. Thus, we have to take global fixings into account when
 *  filtering out symmetries.
 */
static
SCIP_DECL_EVENTEXEC(eventExecSymmetry)
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR* var;
   int varidx;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_SYMMETRY_NAME) == 0 );
   assert( event != NULL );

   propdata = (SCIP_PROPDATA*) eventdata;
   assert( propdata != NULL );
   assert( propdata->permvarmap != NULL );
   assert( propdata->permstrans != NULL );
   assert( propdata->nperms > 0 );
   assert( propdata->permvars != NULL );
   assert( propdata->npermvars > 0 );

   /* get fixed variable */
   var = SCIPeventGetVar(event);
   assert( var != NULL );
   assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );

   if ( ! SCIPhashmapExists(propdata->permvarmap, (void*) var) )
   {
      SCIPerrorMessage("Invalid variable.\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }
   varidx = SCIPhashmapGetImageInt(propdata->permvarmap, (void*) var);
   assert( 0 <= varidx && varidx < propdata->npermvars );

   if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_GUBCHANGED )
   {
      assert( SCIPisEQ(scip, SCIPeventGetNewbound(event), 0.0) );
      assert( SCIPisEQ(scip, SCIPeventGetOldbound(event), 1.0) );

      SCIPdebugMsg(scip, "Mark variable <%s> as globally fixed to 0.\n", SCIPvarGetName(var));
      assert( ! propdata->bg0[varidx] );
      propdata->bg0[varidx] = TRUE;
      propdata->bg0list[propdata->nbg0++] = varidx;
      assert( propdata->nbg0 <= propdata->npermvars );
   }

   if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_GLBCHANGED )
   {
      assert( SCIPisEQ(scip, SCIPeventGetNewbound(event), 1.0) );
      assert( SCIPisEQ(scip, SCIPeventGetOldbound(event), 0.0) );

      SCIPdebugMsg(scip, "Mark variable <%s> as globally fixed to 1.\n", SCIPvarGetName(var));
      assert( ! propdata->bg1[varidx] );
      propdata->bg1[varidx] = TRUE;
      propdata->bg1list[propdata->nbg1++] = varidx;
      assert( propdata->nbg1 <= propdata->npermvars );
   }

   return SCIP_OKAY;
}




/*
 * Table callback methods
 */

/** table data */
struct SCIP_TableData
{
   SCIP_PROPDATA*        propdata;           /** pass data of propagator for table output function */
};


/** output method of orbital fixing propagator statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputOrbitalfixing)
{
   SCIP_TABLEDATA* tabledata;

   assert( scip != NULL );
   assert( table != NULL );

   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );
   assert( tabledata->propdata != NULL );

   if ( tabledata->propdata->nperms > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "Orbital fixing     :\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  vars fixed to 0  :%11d\n", tabledata->propdata->nfixedzero);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  vars fixed to 1  :%11d\n", tabledata->propdata->nfixedone);
   }

   return SCIP_OKAY;
}


/** destructor of statistics table to free user data (called when SCIP is exiting) */
static
SCIP_DECL_TABLEFREE(tableFreeOrbitalfixing)
{
   SCIP_TABLEDATA* tabledata;
   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );

   SCIPfreeBlockMemory(scip, &tabledata);

   return SCIP_OKAY;
}



/*
 * local data structures
 */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(SYMhashGetKeyVartype)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal
 *
 *  Compare the types of two variables according to objective, lower and upper bound, variable type, and column sparsity.
 */
static
SCIP_DECL_HASHKEYEQ(SYMhashKeyEQVartype)
{
   SCIP* scip;
   SYM_VARTYPE* k1;
   SYM_VARTYPE* k2;

   scip = (SCIP*) userptr;
   k1 = (SYM_VARTYPE*) key1;
   k2 = (SYM_VARTYPE*) key2;

   /* first check objective coefficients */
   if ( ! SCIPisEQ(scip, k1->obj, k2->obj) )
      return FALSE;

   /* if still undecided, take lower bound */
   if ( ! SCIPisEQ(scip, k1->lb, k2->lb) )
      return FALSE;

   /* if still undecided, take upper bound */
   if ( ! SCIPisEQ(scip, k1->ub, k2->ub) )
      return FALSE;

   /* if still undecided, take variable type */
   if ( k1->type != k2->type )
      return FALSE;

   /* if still undecided, take number of conss var is contained in */
   if ( k1->nconss != k2->nconss )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(SYMhashKeyValVartype)
{  /*lint --e{715}*/
   SYM_VARTYPE* k;

   k = (SYM_VARTYPE*) key;

   return SCIPhashFour(SCIPrealHashCode(k->obj), SCIPrealHashCode(k->lb), SCIPrealHashCode((double) k->nconss), SCIPrealHashCode(k->ub));
}

/** data structure to store arrays used for sorting rhs types */
struct SYM_Sortrhstype
{
   SCIP_Real*            vals;               /**< array of values */
   SYM_RHSSENSE*         senses;             /**< array of senses of rhs */
   int                   nrhscoef;           /**< size of arrays (for debugging) */
};
typedef struct SYM_Sortrhstype SYM_SORTRHSTYPE;

/** data structure to store arrays used for sorting colored component types */
struct SYM_Sortgraphcompvars
{
   int*                  components;         /**< array of components */
   int*                  colors;             /**< array of colors */
};
typedef struct SYM_Sortgraphcompvars SYM_SORTGRAPHCOMPVARS;

/** sorts rhs types - first by sense, then by value
 *
 *  Due to numerical issues, we first sort by sense, then by value.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortRhsTypes)
{
   SYM_SORTRHSTYPE* data;
   SCIP_Real diffvals;

   data = (SYM_SORTRHSTYPE*) dataptr;
   assert( 0 <= ind1 && ind1 < data->nrhscoef );
   assert( 0 <= ind2 && ind2 < data->nrhscoef );

   /* first sort by senses */
   if ( data->senses[ind1] < data->senses[ind2] )
      return -1;
   else if ( data->senses[ind1] > data->senses[ind2] )
      return 1;

   /* senses are equal, use values */
   diffvals = data->vals[ind1] - data->vals[ind2];

   if ( diffvals < 0.0 )
      return -1;
   else if ( diffvals > 0.0 )
      return 1;

   return 0;
}

/** sorts matrix coefficients
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortMatCoef)
{
   SCIP_Real diffvals;
   SCIP_Real* vals;

   vals = (SCIP_Real*) dataptr;
   diffvals = vals[ind1] - vals[ind2];

   if ( diffvals < 0.0 )
      return -1;
   else if ( diffvals > 0.0 )
      return 1;

   return 0;
}


/** sorts variable indices according to their corresponding component in the graph
 *
 *  Variables are sorted first by the color of their component and then by the component index.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortGraphCompVars)
{
   SYM_SORTGRAPHCOMPVARS* data;

   data = (SYM_SORTGRAPHCOMPVARS*) dataptr;

   if ( data->colors[ind1] < data->colors[ind2] )
      return -1;
   else if ( data->colors[ind1] > data->colors[ind2] )
      return 1;

   if ( data->components[ind1] < data->components[ind2] )
      return -1;
   if ( data->components[ind1] > data->components[ind2] )
      return 1;

   return 0;
}



/*
 * Local methods
 */

#ifndef NDEBUG
/** checks that symmetry data is all freed */
static
SCIP_Bool checkSymmetryDataFree(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert( propdata->permvarmap == NULL );
   assert( propdata->permvarsevents == NULL );
   assert( propdata->bg0list == NULL );
   assert( propdata->bg0 == NULL );
   assert( propdata->bg1list == NULL );
   assert( propdata->bg1 == NULL );
   assert( propdata->nbg0 == 0 );
   assert( propdata->nbg1 == 0 );
   assert( propdata->genorbconss == NULL );
   assert( propdata->genlinconss == NULL );
   assert( propdata->sstconss == NULL );
   assert( propdata->leaders == NULL );

   assert( propdata->permvars == NULL );
   assert( propdata->permvarsobj == NULL );
   assert( propdata->inactiveperms == NULL );
   assert( propdata->perms == NULL );
   assert( propdata->permstrans == NULL );
   assert( propdata->nonbinpermvarcaptured == NULL );
   assert( propdata->npermvars == 0 );
   assert( propdata->nbinpermvars == 0 );
   assert( propdata->nperms == -1 || propdata->nperms == 0 );
   assert( propdata->nmaxperms == 0 );
   assert( propdata->nmovedpermvars == -1 );
   assert( propdata->nmovedbinpermvars == 0 );
   assert( propdata->nmovedintpermvars == 0 );
   assert( propdata->nmovedimplintpermvars == 0 );
   assert( propdata->nmovedcontpermvars == 0 );
   assert( propdata->nmovedvars == -1 );
   assert( propdata->binvaraffected == FALSE );
   assert( propdata->isnonlinvar == NULL );

   assert( propdata->componentblocked == NULL );
   assert( propdata->componentbegins == NULL );
   assert( propdata->components == NULL );
   assert( propdata->ncomponents == -1 );
   assert( propdata->ncompblocked == 0 );

   return TRUE;
}
#endif


/** checks whether a variable has a type compatible with the leader vartype */
static
SCIP_Bool isLeadervartypeCompatible(
   SCIP_VAR*             var,                /**< variable to check */
   int                   leadervartype       /**< bit set encoding possible leader variable types */
   )
{
   SCIP_VARTYPE vartype;
   unsigned int vartypeencoding;

   assert( var != NULL );
   assert( leadervartype >= 0 );
   assert( leadervartype <= 15 );

   vartype = SCIPvarGetType(var);

   if ( vartype == SCIP_VARTYPE_BINARY )
      vartypeencoding = 1;
   else if ( vartype == SCIP_VARTYPE_INTEGER )
      vartypeencoding = 2;
   else if ( vartype == SCIP_VARTYPE_IMPLINT )
      vartypeencoding = 4;
   else
      vartypeencoding = 8;

   return (SCIP_Bool) (vartypeencoding & (unsigned) leadervartype);
}


/** returns whether a recomputation of symmetries is required */
static
SCIP_Bool isSymmetryRecomputationRequired(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{ /*lint --e{641}*/
   assert( scip != NULL );
   assert( propdata != NULL );

   if ( propdata->recomputerestart == SCIP_RECOMPUTESYM_NEVER )
      return FALSE;

   /* we do not need to recompute symmetries if no restart has occured */
   if ( SCIPgetNRuns(scip) == propdata->lastrestart || propdata->lastrestart == 0 || SCIPgetNRuns(scip) == 1 )
      return FALSE;

   if ( propdata->recomputerestart == SCIP_RECOMPUTESYM_ALWAYS )
      return TRUE;

   /* recompute symmetries if OF found a reduction */
   assert( propdata->recomputerestart == SCIP_RECOMPUTESYM_OFFOUNDRED );
   if ( propdata->offoundreduction )
      return TRUE;

   return FALSE;
}


/** sets in propdata which symmetry handling methods are active */
static
SCIP_RETCODE setSymmetryMethodEnabledValues(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert( propdata != NULL );

   if ( ISSYMRETOPESACTIVE(propdata->usesymmetry) )
      propdata->symconsenabled = TRUE;
   else
      propdata->symconsenabled = FALSE;

   if ( ISORBITALFIXINGACTIVE(propdata->usesymmetry) )
      propdata->ofenabled = TRUE;
   else
      propdata->ofenabled = FALSE;

   if ( ISSSTACTIVE(propdata->usesymmetry) )
      propdata->sstenabled = TRUE;
   else
      propdata->sstenabled = FALSE;

   return SCIP_OKAY;
}


/** frees symmetry data */
static
SCIP_RETCODE freeSymmetryData(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );

   if ( propdata->permvarmap != NULL )
   {
      SCIPhashmapFree(&propdata->permvarmap);
   }

   /* drop events */
   if ( propdata->permvarsevents != NULL )
   {
      assert( propdata->permvars != NULL );
      assert( propdata->npermvars > 0 );

      for (i = 0; i < propdata->npermvars; ++i)
      {
         if ( SCIPvarGetType(propdata->permvars[i]) == SCIP_VARTYPE_BINARY )
         {
            /* If symmetry is computed before presolving, it might happen that some variables are turned into binary
             * variables, for which no event has been catched. Since there currently is no way of checking whether a var
             * event has been caught for a particular variable, we use the stored eventfilter positions. */
            if ( propdata->permvarsevents[i] >= 0 )
            {
               SCIP_CALL( SCIPdropVarEvent(scip, propdata->permvars[i], SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
                     propdata->eventhdlr, (SCIP_EVENTDATA*) propdata, propdata->permvarsevents[i]) );
            }
         }
      }
      SCIPfreeBlockMemoryArray(scip, &propdata->permvarsevents, propdata->npermvars);
   }

   /* release variables */
   if ( propdata->nonbinpermvarcaptured != NULL )
   {
      int cnt;

      /* memory should have been allocated only if the leader type is not binary */
      assert( propdata->sstenabled && propdata->sstleadervartype != (int) SCIP_SSTTYPE_BINARY );

      for (i = propdata->nbinpermvars, cnt = 0; i < propdata->npermvars; ++i, ++cnt)
      {
         /* release captured non-binary variables
          * (cannot use isLeadervartypeCompatible(), because vartype may have changed in between)
          */
         if ( propdata->nonbinpermvarcaptured[cnt] )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &propdata->permvars[i]) );
         }
      }
      SCIPfreeBlockMemoryArray(scip, &propdata->nonbinpermvarcaptured, propdata->npermvars - propdata->nbinpermvars);
      propdata->nonbinpermvarcaptured = NULL;
   }

   if ( propdata->binvaraffected )
   {
      for (i = 0; i < propdata->nbinpermvars; ++i)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &propdata->permvars[i]) );
      }
   }

   /* free lists for orbitopal fixing */
   if ( propdata->bg0list != NULL )
   {
      assert( propdata->bg0 != NULL );
      assert( propdata->bg1list != NULL );
      assert( propdata->bg1 != NULL );

      SCIPfreeBlockMemoryArray(scip, &propdata->bg0list, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->bg0, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->bg1list, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->bg1, propdata->npermvars);

      propdata->nbg0 = 0;
      propdata->nbg1 = 0;
   }

   /* other data */
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->inactiveperms, propdata->nperms);

   /* free permstrans matrix*/
   if ( propdata->permstrans != NULL )
   {
      assert( propdata->nperms > 0 );
      assert( propdata->permvars != NULL );
      assert( propdata->npermvars > 0 );
      assert( propdata->nmaxperms > 0 );

      for (i = 0; i < propdata->npermvars; ++i)
      {
         SCIPfreeBlockMemoryArray(scip, &propdata->permstrans[i], propdata->nmaxperms);
      }
      SCIPfreeBlockMemoryArray(scip, &propdata->permstrans, propdata->npermvars);
   }

   /* free data of added orbitope/orbisack/symresack constraints */
   if ( propdata->genorbconss != NULL )
   {
      assert( propdata->ngenorbconss + propdata->ngenlinconss > 0
         || (ISORBITALFIXINGACTIVE(propdata->usesymmetry) && propdata->norbitopes == 0) );

      /* release constraints */
      for (i = 0; i < propdata->ngenorbconss; ++i)
      {
         assert( propdata->genorbconss[i] != NULL );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->genorbconss[i]) );
      }

      /* free pointers to symmetry group and binary variables */
      SCIPfreeBlockMemoryArray(scip, &propdata->genorbconss, propdata->nperms);
      propdata->ngenorbconss = 0;
   }

   /* free data of added constraints */
   if ( propdata->genlinconss != NULL )
   {
      /* release constraints */
      for (i = 0; i < propdata->ngenlinconss; ++i)
      {
         assert( propdata->genlinconss[i] != NULL );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->genlinconss[i]) );
      }

      /* free pointers to symmetry group and binary variables */
      SCIPfreeBlockMemoryArray(scip, &propdata->genlinconss, propdata->genlinconsssize);
      propdata->ngenlinconss = 0;
      propdata->genlinconsssize = 0;
   }

   if ( propdata->sstconss != NULL )
   {
      assert( propdata->nsstconss > 0 );

      /* release constraints */
      for (i = 0; i < propdata->nsstconss; ++i)
      {
         assert( propdata->sstconss[i] != NULL );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->sstconss[i]) );
      }

      /* free pointers to symmetry group and binary variables */
      SCIPfreeBlockMemoryArray(scip, &propdata->sstconss, propdata->maxnsstconss);
      propdata->sstconss = NULL;
      propdata->nsstconss = 0;
      propdata->maxnsstconss = 0;
   }

   if ( propdata->leaders != NULL )
   {
      assert( propdata->maxnleaders > 0 );

      SCIPfreeBlockMemoryArray(scip, &propdata->leaders, propdata->maxnleaders);
      propdata->maxnleaders = 0;
      propdata->leaders = NULL;
      propdata->nleaders = 0;
   }

   /* free components */
   if ( propdata->ncomponents > 0 )
   {
      assert( propdata->componentblocked != NULL );
      assert( propdata->vartocomponent != NULL );
      assert( propdata->componentbegins != NULL );
      assert( propdata->components != NULL );

      SCIPfreeBlockMemoryArray(scip, &propdata->componentblocked, propdata->ncomponents);
      SCIPfreeBlockMemoryArray(scip, &propdata->vartocomponent, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->componentbegins, propdata->ncomponents + 1);
      SCIPfreeBlockMemoryArray(scip, &propdata->components, propdata->nperms);

      propdata->ncomponents = -1;
      propdata->ncompblocked = 0;
   }

   /* free main symmetry data */
   if ( propdata->nperms > 0 )
   {
      assert( propdata->permvars != NULL );

      SCIPfreeBlockMemoryArray(scip, &propdata->permvars, propdata->npermvars);

      /* if orbital fixing runs exclusively, propdata->perms was already freed in determineSymmetry() */
      if ( propdata->perms != NULL )
      {
         for (i = 0; i < propdata->nperms; ++i)
         {
            SCIPfreeBlockMemoryArray(scip, &propdata->perms[i], propdata->npermvars);
         }
         SCIPfreeBlockMemoryArray(scip, &propdata->perms, propdata->nmaxperms);
      }

#ifndef NDEBUG
      SCIPfreeBlockMemoryArrayNull(scip, &propdata->permvarsobj, propdata->npermvars);
#endif

      SCIPfreeBlockMemoryArrayNull(scip, &propdata->isnonlinvar, propdata->npermvars);

      propdata->npermvars = 0;
      propdata->nbinpermvars = 0;
      propdata->nperms = -1;
      propdata->nmaxperms = 0;
      propdata->nmovedpermvars = -1;
      propdata->nmovedbinpermvars = 0;
      propdata->nmovedintpermvars = 0;
      propdata->nmovedimplintpermvars = 0;
      propdata->nmovedcontpermvars = 0;
      propdata->nmovedvars = -1;
      propdata->log10groupsize = -1.0;
      propdata->binvaraffected = FALSE;
      propdata->isnonlinvar = NULL;
   }
   propdata->nperms = -1;

   assert( checkSymmetryDataFree(propdata) );

   propdata->computedsymmetry = FALSE;
   propdata->compressed = FALSE;

   return SCIP_OKAY;
}


/** deletes symmetry handling constraints */
static
SCIP_RETCODE delSymConss(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );

   if ( propdata->ngenorbconss == 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &propdata->genorbconss, propdata->nperms);
   }
   else
   {
      assert( propdata->genorbconss != NULL );
      assert( propdata->nperms > 0 );
      assert( propdata->nperms >= propdata->ngenorbconss );

      for (i = 0; i < propdata->ngenorbconss; ++i)
      {
         assert( propdata->genorbconss[i] != NULL );

         SCIP_CALL( SCIPdelCons(scip, propdata->genorbconss[i]) );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->genorbconss[i]) );
      }

      SCIPfreeBlockMemoryArray(scip, &propdata->genorbconss, propdata->nperms);
      propdata->ngenorbconss = 0;
   }

   /* free Schreier Sims data */
   if ( propdata->nsstconss > 0 )
   {
      for (i = 0; i < propdata->nsstconss; ++i)
      {
         assert( propdata->sstconss[i] != NULL );

         SCIP_CALL( SCIPdelCons(scip, propdata->sstconss[i]) );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->sstconss[i]) );
      }

      SCIPfreeBlockMemoryArray(scip, &propdata->sstconss, propdata->maxnsstconss);
      propdata->nsstconss = 0;
      propdata->maxnsstconss = 0;
   }

   if ( propdata->ngenlinconss == 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &propdata->genlinconss, propdata->genlinconsssize);
   }
   else
   {
      assert( propdata->genlinconss != NULL );
      assert( propdata->nperms > 0 );

      for (i = 0; i < propdata->ngenlinconss; ++i)
      {
         assert( propdata->genlinconss[i] != NULL );

         SCIP_CALL( SCIPdelCons(scip, propdata->genlinconss[i]) );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->genlinconss[i]) );
      }

      SCIPfreeBlockMemoryArray(scip, &propdata->genlinconss, propdata->genlinconsssize);
      propdata->ngenlinconss = 0;
   }

   /* free pointers to symmetry group and binary variables */
   assert( propdata->nperms > 0 || propdata->genorbconss == NULL );
   assert( propdata->nperms >= propdata->ngenorbconss || propdata->genorbconss == NULL );
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->genorbconss, propdata->nperms);
   propdata->ngenorbconss = 0;
   propdata->triedaddconss = FALSE;

   return SCIP_OKAY;
}


/** determines whether variable should be fixed by permutations */
static
SCIP_Bool SymmetryFixVar(
   SYM_SPEC              fixedtype,          /**< bitset of variable types that should be fixed */
   SCIP_VAR*             var                 /**< variable to be considered */
   )
{
   if ( (fixedtype & SYM_SPEC_INTEGER) && SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
      return TRUE;
   if ( (fixedtype & SYM_SPEC_BINARY) && SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      return TRUE;
   if ( (fixedtype & SYM_SPEC_REAL) &&
      (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT) )
      return TRUE;
   return FALSE;
}


/** Transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant.
 *
 *  @note @p constant needs to be initialized!
 */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( scalars != NULL );
   assert( *vars != NULL );
   assert( *scalars != NULL );
   assert( nvars != NULL );
   assert( constant != NULL );

   if ( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if ( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, scalars, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for (v = 0; v < *nvars; ++v)
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&(*vars)[v], &(*scalars)[v], constant) );
      }
   }
   return SCIP_OKAY;
}


/** fills in matrix elements into coefficient arrays */
static
SCIP_RETCODE collectCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             doubleequations,    /**< Double equations to positive/negative version? */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficients values (or NULL if all linear coefficient values are 1) */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             istransformed,      /**< whether the constraint is transformed */
   SYM_RHSSENSE          rhssense,           /**< identifier of constraint type */
   SYM_MATRIXDATA*       matrixdata,         /**< matrix data to be filled in */
   int*                  nconssforvar        /**< pointer to array to store for each var the number of conss */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real constant = 0.0;
   int nrhscoef;
   int nmatcoef;
   int nvars;
   int j;

   assert( scip != NULL );
   assert( nlinvars == 0 || linvars != NULL );
   assert( lhs <= rhs );

   /* do nothing if constraint is empty */
   if ( nlinvars == 0 )
      return SCIP_OKAY;

   /* ignore redundant constraints */
   if ( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   /* duplicate variable and value array */
   nvars = nlinvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, linvars, nvars) );
   if ( linvals != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &vals, linvals, nvars) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
      for (j = 0; j < nvars; ++j)
         vals[j] = 1.0;
   }
   assert( vars != NULL );
   assert( vals != NULL );

   /* get active variables */
   SCIP_CALL( getActiveVariables(scip, &vars, &vals, &nvars, &constant, istransformed) );

   /* check whether constraint is empty after transformation to active variables */
   if ( nvars <= 0 )
   {
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
      return SCIP_OKAY;
   }

   /* handle constant */
   if ( ! SCIPisInfinity(scip, -lhs) )
      lhs -= constant;
   if ( ! SCIPisInfinity(scip, rhs) )
      rhs -= constant;

   /* check whether we have to resize; note that we have to add 2 * nvars since two inequalities may be added */
   if ( matrixdata->nmatcoef + 2 * nvars > matrixdata->nmaxmatcoef )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, matrixdata->nmatcoef + 2 * nvars);
      assert( newsize >= 0 );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(matrixdata->matidx), matrixdata->nmaxmatcoef, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(matrixdata->matrhsidx), matrixdata->nmaxmatcoef, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(matrixdata->matvaridx), matrixdata->nmaxmatcoef, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(matrixdata->matcoef), matrixdata->nmaxmatcoef, newsize) );
      SCIPdebugMsg(scip, "Resized matrix coefficients from %d to %d.\n", matrixdata->nmaxmatcoef, newsize);
      matrixdata->nmaxmatcoef = newsize;
   }

   nrhscoef = matrixdata->nrhscoef;
   nmatcoef = matrixdata->nmatcoef;

   /* check lhs/rhs */
   if ( SCIPisEQ(scip, lhs, rhs) )
   {
      SCIP_Bool poscoef = FALSE;
      SCIP_Bool negcoef = FALSE;

      assert( ! SCIPisInfinity(scip, rhs) );

      /* equality constraint */
      matrixdata->rhscoef[nrhscoef] = rhs;

      /* if we deal with special constraints */
      if ( rhssense >= SYM_SENSE_XOR )
         matrixdata->rhssense[nrhscoef] = rhssense;
      else
         matrixdata->rhssense[nrhscoef] = SYM_SENSE_EQUATION;
      matrixdata->rhsidx[nrhscoef] = nrhscoef;

      for (j = 0; j < nvars; ++j)
      {
         assert( nmatcoef < matrixdata->nmaxmatcoef );

         matrixdata->matidx[nmatcoef] = nmatcoef;
         matrixdata->matrhsidx[nmatcoef] = nrhscoef;

         assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

         if ( nconssforvar != NULL )
            nconssforvar[SCIPvarGetProbindex(vars[j])] += 1;
         matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);
         matrixdata->matcoef[nmatcoef++] = vals[j];
         if ( SCIPisPositive(scip, vals[j]) )
            poscoef = TRUE;
         else
            negcoef = TRUE;
      }
      nrhscoef++;

      /* add negative of equation; increases chance to detect symmetry, but might increase time to compute symmetry. */
      if ( doubleequations && poscoef && negcoef )
      {
         for (j = 0; j < nvars; ++j)
         {
            assert( nmatcoef < matrixdata->nmaxmatcoef );
            assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

            matrixdata->matidx[nmatcoef] = nmatcoef;
            matrixdata->matrhsidx[nmatcoef] = nrhscoef;
            matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);
            matrixdata->matcoef[nmatcoef++] = -vals[j];
         }
         matrixdata->rhssense[nrhscoef] = SYM_SENSE_EQUATION;
         matrixdata->rhsidx[nrhscoef] = nrhscoef;
         matrixdata->rhscoef[nrhscoef++] = -rhs;
      }
   }
   else
   {
#ifndef NDEBUG
      if ( rhssense == SYM_SENSE_BOUNDIS_TYPE_2 )
      {
         assert( ! SCIPisInfinity(scip, -lhs) );
         assert( ! SCIPisInfinity(scip, rhs) );
      }
#endif

      if ( ! SCIPisInfinity(scip, -lhs) )
      {
         matrixdata->rhscoef[nrhscoef] = -lhs;
         if ( rhssense >= SYM_SENSE_XOR )
         {
            assert( rhssense == SYM_SENSE_BOUNDIS_TYPE_2 );
            matrixdata->rhssense[nrhscoef] = rhssense;
         }
         else
            matrixdata->rhssense[nrhscoef] = SYM_SENSE_INEQUALITY;

         matrixdata->rhsidx[nrhscoef] = nrhscoef;

         for (j = 0; j < nvars; ++j)
         {
            assert( nmatcoef < matrixdata->nmaxmatcoef );
            matrixdata->matidx[nmatcoef] = nmatcoef;
            matrixdata->matrhsidx[nmatcoef] = nrhscoef;
            matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);

            assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

            if ( nconssforvar != NULL )
               nconssforvar[SCIPvarGetProbindex(vars[j])] += 1;

            matrixdata->matcoef[nmatcoef++] = -vals[j];
         }
         nrhscoef++;
      }

      if ( ! SCIPisInfinity(scip, rhs) )
      {
         matrixdata->rhscoef[nrhscoef] = rhs;
         if ( rhssense >= SYM_SENSE_XOR )
         {
            assert( rhssense == SYM_SENSE_BOUNDIS_TYPE_2 );
            matrixdata->rhssense[nrhscoef] = rhssense;
         }
         else
            matrixdata->rhssense[nrhscoef] = SYM_SENSE_INEQUALITY;

         matrixdata->rhsidx[nrhscoef] = nrhscoef;

         for (j = 0; j < nvars; ++j)
         {
            assert( nmatcoef < matrixdata->nmaxmatcoef );
            matrixdata->matidx[nmatcoef] = nmatcoef;
            matrixdata->matrhsidx[nmatcoef] = nrhscoef;

            assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

            if ( nconssforvar != NULL )
               nconssforvar[SCIPvarGetProbindex(vars[j])] += 1;

            matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);
            matrixdata->matcoef[nmatcoef++] = vals[j];
         }
         nrhscoef++;
      }
   }
   matrixdata->nrhscoef = nrhscoef;
   matrixdata->nmatcoef = nmatcoef;

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** checks whether given permutations form a symmetry of a MIP
 *
 *  We need the matrix and rhs in the original order in order to speed up the comparison process. The matrix is needed
 *  in the right order to easily check rows. The rhs is used because of cache effects.
 */
static
SCIP_RETCODE checkSymmetriesAreSymmetries(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SPEC              fixedtype,          /**< variable types that must be fixed by symmetries */
   SYM_MATRIXDATA*       matrixdata,         /**< matrix data */
   int                   nperms,             /**< number of permutations */
   int**                 perms               /**< permutations */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_HASHMAP* varmap;
   SCIP_VAR** occuringvars;
   SCIP_Real* permrow = 0;
   SCIP_Bool success;
   int* rhsmatbeg = 0;
   int nconss;
   int noccuringvars;
   int oldrhs;
   int i;
   int j;
   int p;

   SCIPdebugMsg(scip, "Checking whether symmetries are symmetries (generators: %d).\n", nperms);

   /* set up dense row for permuted row */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &permrow, matrixdata->npermvars) );

   /* set up map between rows and first entry in matcoef array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rhsmatbeg, matrixdata->nrhscoef) );
   for (j = 0; j < matrixdata->nrhscoef; ++j)
      rhsmatbeg[j] = -1;

   /* get info for non-linear part */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   nconss = conshdlr != NULL ? SCIPconshdlrGetNConss(conshdlr) : 0;

   /* create hashmaps for variable permutation and constraints in non-linear part array for occuring variables */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), matrixdata->npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &occuringvars, matrixdata->npermvars) );

   /* build map from rhs into matrix */
   oldrhs = -1;
   for (j = 0; j < matrixdata->nmatcoef; ++j)
   {
      int rhs;

      rhs = matrixdata->matrhsidx[j];
      if ( rhs != oldrhs )
      {
         assert( 0 <= rhs && rhs < matrixdata->nrhscoef );
         rhsmatbeg[rhs] = j;
         oldrhs = rhs;
      }
   }

   /* create row */
   for (j = 0; j < matrixdata->npermvars; ++j)
      permrow[j] = 0.0;

   /* check all generators */
   for (p = 0; p < nperms; ++p)
   {
      int* P;
      int r1;
      int r2;

      SCIPdebugMsg(scip, "Verifying automorphism group generator #%d for linear part ...\n", p);
      P = perms[p];
      assert( P != NULL );

      for (j = 0; j < matrixdata->npermvars; ++j)
      {
         if ( SymmetryFixVar(fixedtype, matrixdata->permvars[j]) && P[j] != j )
         {
            SCIPdebugMsg(scip, "Permutation does not fix types %u, moving variable %d.\n", fixedtype, j);
            return SCIP_ERROR;
         }
      }

      /*
       *  linear part
       */

      /* check all linear constraints == rhs */
      for (r1 = 0; r1 < matrixdata->nrhscoef; ++r1)
      {
         int npermuted = 0;

         /* fill row into permrow (dense) */
         j = rhsmatbeg[r1];
         assert( 0 <= j && j < matrixdata->nmatcoef );
         assert( matrixdata->matrhsidx[j] == r1 ); /* note: row cannot be empty by construction */

         /* loop through row */
         while ( j < matrixdata->nmatcoef && matrixdata->matrhsidx[j] == r1 )
         {
            int varidx;

            assert( matrixdata->matvaridx[j] < matrixdata->npermvars );
            varidx = P[matrixdata->matvaridx[j]];
            assert( 0 <= varidx && varidx < matrixdata->npermvars );
            if ( varidx != matrixdata->matvaridx[j] )
               ++npermuted;
            assert( SCIPisZero(scip, permrow[varidx]) );
            permrow[varidx] = matrixdata->matcoef[j];
            ++j;
         }

         /* if row is not affected by permutation, we do not have to check it */
         if ( npermuted > 0 )
         {
            /* check other rows (sparse) */
            SCIP_Bool found = FALSE;
            for (r2 = 0; r2 < matrixdata->nrhscoef; ++r2)
            {
               /* a permutation must map constraints of the same type and respect rhs coefficients */
               if ( matrixdata->rhssense[r1] == matrixdata->rhssense[r2] && SCIPisEQ(scip, matrixdata->rhscoef[r1], matrixdata->rhscoef[r2]) )
               {
                  j = rhsmatbeg[r2];
                  assert( 0 <= j && j < matrixdata->nmatcoef );
                  assert( matrixdata->matrhsidx[j] == r2 );
                  assert( matrixdata->matvaridx[j] < matrixdata->npermvars );

                  /* loop through row r2 and check whether it is equal to permuted row r */
                  while ( j < matrixdata->nmatcoef && matrixdata->matrhsidx[j] == r2 && SCIPisEQ(scip, permrow[matrixdata->matvaridx[j]], matrixdata->matcoef[j] ) )
                     ++j;

                  /* check whether rows are completely equal */
                  if ( j >= matrixdata->nmatcoef || matrixdata->matrhsidx[j] != r2 )
                  {
                     /* perm[p] is indeed a symmetry */
                     found = TRUE;
                     break;
                  }
               }
            }

            assert( found );
            if ( ! found ) /*lint !e774*/
            {
               SCIPerrorMessage("Found permutation that is not a symmetry.\n");
               return SCIP_ERROR;
            }
         }

         /* reset permrow */
         j = rhsmatbeg[r1];
         while ( j < matrixdata->nmatcoef && matrixdata->matrhsidx[j] == r1 )
         {
            int varidx;
            varidx = P[matrixdata->matvaridx[j]];
            permrow[varidx] = 0.0;
            ++j;
         }
      }

      /*
       *  non-linear part
       */

      SCIPdebugMsg(scip, "Verifying automorphism group generator #%d for non-linear part ...\n", p);

      /* fill hashmap according to permutation */
      for (j = 0; j < matrixdata->npermvars; ++j)
      {
         SCIP_CALL( SCIPhashmapInsert(varmap, matrixdata->permvars[j], matrixdata->permvars[P[j]]) );
      }

      /* check all non-linear constraints */
      for (i = 0; i < nconss; ++i)
      {
         SCIP_CONS* cons1;
         SCIP_Bool permuted = FALSE;

         cons1 = SCIPconshdlrGetConss(conshdlr)[i];

         SCIP_CALL( SCIPgetConsVars(scip, cons1, occuringvars, matrixdata->npermvars, &success) );
         assert(success);
         SCIP_CALL( SCIPgetConsNVars(scip, cons1, &noccuringvars, &success) );
         assert(success);

         /* count number of affected variables in this constraint */
         for (j = 0; j < noccuringvars && ! permuted; ++j)
         {
            int varidx;

            varidx = SCIPvarGetProbindex(occuringvars[j]);
            assert( varidx >= 0 && varidx < matrixdata->npermvars );

            if ( P[varidx] != varidx )
               permuted = TRUE;
         }

         /* if constraint is not affected by permutation, we do not have to check it */
         if ( permuted )
         {
            SCIP_CONS* permutedcons = NULL;
            SCIP_EXPR* permutedexpr;
            SCIP_Bool found = FALSE;
            SCIP_Bool infeasible;

            /* copy contraints but exchange variables according to hashmap */
            SCIP_CALL( SCIPgetConsCopy(scip, scip, cons1, &permutedcons, conshdlr, varmap, NULL, NULL,
                  SCIPconsIsInitial(cons1), SCIPconsIsSeparated(cons1), SCIPconsIsEnforced(cons1),
                  SCIPconsIsChecked(cons1), SCIPconsIsPropagated(cons1), SCIPconsIsLocal(cons1),
                  SCIPconsIsModifiable(cons1), SCIPconsIsDynamic(cons1), SCIPconsIsRemovable(cons1),
                  SCIPconsIsStickingAtNode(cons1), FALSE, &success) );
            assert(success);
            assert(permutedcons != NULL);

            /* simplify permuted expr in order to guarantee sorted variables */
            permutedexpr = SCIPgetExprNonlinear(permutedcons);
            SCIP_CALL( SCIPsimplifyExpr(scip, permutedexpr, &permutedexpr, &success, &infeasible, NULL, NULL) );
            assert( !infeasible );

            /* look for a constraint with same lhs, rhs and expression */
            for (j = 0; j < nconss; ++j)
            {
               SCIP_CONS* cons2;

               cons2 = SCIPconshdlrGetConss(conshdlr)[j];

               if ( SCIPisEQ(scip, SCIPgetRhsNonlinear(cons2), SCIPgetRhsNonlinear(permutedcons))
                  && SCIPisEQ(scip, SCIPgetLhsNonlinear(cons2), SCIPgetLhsNonlinear(permutedcons))
                  && (SCIPcompareExpr(scip, SCIPgetExprNonlinear(cons2), permutedexpr) == 0) )
               {
                  found = TRUE;
                  break;
               }
            }

            /* release copied constraint and expression because simplify captures it */
            SCIP_CALL( SCIPreleaseExpr(scip, &permutedexpr) );
            SCIP_CALL( SCIPreleaseCons(scip, &permutedcons) );

            assert( found );
            if( !found ) /*lint !e774*/
            {
               SCIPerrorMessage("Found permutation that is not a symmetry.\n");
               return SCIP_ERROR;
            }
         }
      }

      /* reset varmap */
      SCIP_CALL( SCIPhashmapRemoveAll(varmap) );
   }

   SCIPhashmapFree(&varmap);
   SCIPfreeBufferArray(scip, &occuringvars);
   SCIPfreeBlockMemoryArray(scip, &rhsmatbeg, matrixdata->nrhscoef);
   SCIPfreeBlockMemoryArray(scip, &permrow, matrixdata->npermvars);

   return SCIP_OKAY;
}


/** returns the number of active constraints that can be handled by symmetry */
static
int getNSymhandableConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONSHDLR*        conshdlr_nonlinear  /**< nonlinear constraint handler, if included */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int nhandleconss = 0;

   assert( scip != NULL );

   conshdlr = SCIPfindConshdlr(scip, "linear");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "linking");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "setppc");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "xor");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "and");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "or");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "logicor");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "knapsack");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "varbound");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "bounddisjunction");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   if( conshdlr_nonlinear != NULL )
      nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr_nonlinear);

   return nhandleconss;
}

/** returns whether there are any active nonlinear constraints */
static
SCIP_Bool hasNonlinearConstraints(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   return propdata->conshdlr_nonlinear != NULL && SCIPconshdlrGetNActiveConss(propdata->conshdlr_nonlinear) > 0;
}

/** set symmetry data */
static
SCIP_RETCODE setSymmetryData(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR**            vars,               /**< vars present at time of symmetry computation */
   int                   nvars,              /**< number of vars present at time of symmetry computation */
   int                   nbinvars,           /**< number of binary vars present at time of symmetry computation */
   SCIP_VAR***           permvars,           /**< pointer to permvars array */
   int*                  npermvars,          /**< pointer to store number of permvars */
   int*                  nbinpermvars,       /**< pointer to store number of binary permvars */
   int**                 perms,              /**< permutations matrix (nperms x nvars) */
   int                   nperms,             /**< number of permutations */
   int*                  nmovedvars,         /**< pointer to store number of vars affected by symmetry (if usecompression) or NULL */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether a binary variable is affected by symmetry */
   SCIP_Bool             usecompression,     /**< whether symmetry data shall be compressed */
   SCIP_Real             compressthreshold,  /**< if percentage of moved vars is at most threshold, compression is done */
   SCIP_Bool*            compressed          /**< pointer to store whether compression has been performed */
   )
{
   int i;
   int p;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( permvars != NULL );
   assert( npermvars != NULL );
   assert( nbinpermvars != NULL );
   assert( perms != NULL );
   assert( nperms > 0 );
   assert( binvaraffected != NULL );
   assert( SCIPisGE(scip, compressthreshold, 0.0) );
   assert( SCIPisLE(scip, compressthreshold, 1.0) );
   assert( compressed != NULL );

   /* set default return values */
   *permvars = vars;
   *npermvars = nvars;
   *nbinpermvars = nbinvars;
   *binvaraffected = FALSE;
   *compressed = FALSE;

   /* if we possibly perform compression */
   if ( usecompression && SCIPgetNVars(scip) >= COMPRESSNVARSLB )
   {
      SCIP_Real percentagemovedvars;
      int* labelmovedvars;
      int* labeltopermvaridx;
      int nbinvarsaffected = 0;

      assert( nmovedvars != NULL );

      *nmovedvars = 0;

      /* detect number of moved vars and label moved vars */
      SCIP_CALL( SCIPallocBufferArray(scip, &labelmovedvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &labeltopermvaridx, nvars) );
      for (i = 0; i < nvars; ++i)
      {
         labelmovedvars[i] = -1;

         for (p = 0; p < nperms; ++p)
         {
            if ( perms[p][i] != i )
            {
               labeltopermvaridx[*nmovedvars] = i;
               labelmovedvars[i] = (*nmovedvars)++;

               if ( SCIPvarIsBinary(vars[i]) )
                  ++nbinvarsaffected;
               break;
            }
         }
      }

      if ( nbinvarsaffected > 0 )
         *binvaraffected = TRUE;

      /* check whether compression should be performed */
      percentagemovedvars = (SCIP_Real) *nmovedvars / (SCIP_Real) nvars;
      if ( *nmovedvars > 0 && SCIPisLE(scip, percentagemovedvars, compressthreshold) )
      {
         /* remove variables from permutations that are not affected by any permutation */
         for (p = 0; p < nperms; ++p)
         {
            /* iterate over labels and adapt permutation */
            for (i = 0; i < *nmovedvars; ++i)
            {
               assert( i <= labeltopermvaridx[i] );
               perms[p][i] = labelmovedvars[perms[p][labeltopermvaridx[i]]];
            }

            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &perms[p], nvars, *nmovedvars) );
         }

         /* remove variables from permvars array that are not affected by any symmetry */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, permvars, *nmovedvars) );
         for (i = 0; i < *nmovedvars; ++i)
         {
            (*permvars)[i] = vars[labeltopermvaridx[i]];
         }
         *npermvars = *nmovedvars;
         *nbinpermvars = nbinvarsaffected;
         *compressed = TRUE;

         SCIPfreeBlockMemoryArray(scip, &vars, nvars);
      }
      SCIPfreeBufferArray(scip, &labeltopermvaridx);
      SCIPfreeBufferArray(scip, &labelmovedvars);
   }
   else
   {
      /* detect whether binary variable is affected by symmetry and count number of binary permvars */
      for (i = 0; i < nbinvars; ++i)
      {
         for (p = 0; p < nperms && ! *binvaraffected; ++p)
         {
            if ( perms[p][i] != i )
            {
               if ( SCIPvarIsBinary(vars[i]) )
                  *binvaraffected = TRUE;
               break;
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** computes symmetry group of a MIP */
static
SCIP_RETCODE computeSymmetryGroup(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_Bool             doubleequations,    /**< Double equations to positive/negative version? */
   SCIP_Bool             compresssymmetries, /**< Should non-affected variables be removed from permutation to save memory? */
   SCIP_Real             compressthreshold,  /**< Compression is used if percentage of moved vars is at most the threshold. */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_SPEC              fixedtype,          /**< variable types that must be fixed by symmetries */
   SCIP_Bool             local,              /**< Use local variable bounds? */
   SCIP_Bool             checksymmetries,    /**< Should all symmetries be checked after computation? */
   SCIP_Bool             usecolumnsparsity,  /**< Should the number of conss a variable is contained in be exploited in symmetry detection? */
   SCIP_CONSHDLR*        conshdlr_nonlinear, /**< Nonlinear constraint handler, if included */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   int*                  nbinpermvars,       /**< pointer to store number of binary variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   int*                  nmovedvars,         /**< pointer to store number of moved vars */
   SCIP_Bool**           isnonlinvar,        /**< pointer to store which variables appear nonlinearly */
   SCIP_Bool*            binvaraffected,     /**< pointer to store wether a binary variable is affected by symmetry */
   SCIP_Bool*            compressed,         /**< pointer to store whether compression has been performed */
   SCIP_Bool*            success             /**< pointer to store whether symmetry computation was successful */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SYM_MATRIXDATA matrixdata;
   SYM_EXPRDATA exprdata;
   SCIP_HASHTABLE* vartypemap;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_EXPRITER* it = NULL;
   SCIP_HASHSET* auxvars = NULL;
   SYM_VARTYPE* uniquevararray;
   SYM_RHSSENSE oldsense = SYM_SENSE_UNKOWN;
   SYM_SORTRHSTYPE sortrhstype;
   SCIP_Real oldcoef = SCIP_INVALID;
   SCIP_Real val;
   int* nconssforvar = NULL;
   int nuniquevararray = 0;
   int nhandleconss;
   int nactiveconss;
   int nnlconss;
   int nconss;
   int nvars;
   int nbinvars;
   int nvarsorig;
   int nallvars;
   int c;
   int j;

   assert( scip != NULL );
   assert( npermvars != NULL );
   assert( nbinpermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );
   assert( binvaraffected != NULL );
   assert( compressed != NULL );
   assert( success != NULL );
   assert( isnonlinvar != NULL );
   assert( SYMcanComputeSymmetry() );

   /* init */
   *npermvars = 0;
   *nbinpermvars = 0;
   *permvars = NULL;
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;
   *nmovedvars = -1;
   *binvaraffected = FALSE;
   *compressed = FALSE;
   *success = FALSE;

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   nvarsorig = nvars;

   /* exit if no constraints or no variables are available */
   if ( nconss == 0 || nvars == 0 )
   {
      *success = TRUE;
      return SCIP_OKAY;
   }

   conss = SCIPgetConss(scip);
   assert( conss != NULL );

   /* compute the number of active constraints */
   nactiveconss = SCIPgetNActiveConss(scip);
   nnlconss = conshdlr_nonlinear != NULL ? SCIPconshdlrGetNActiveConss(conshdlr_nonlinear) : 0;

   /* exit if no active constraints are available */
   if ( nactiveconss == 0 )
   {
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* before we set up the matrix, check whether we can handle all constraints */
   nhandleconss = getNSymhandableConss(scip, conshdlr_nonlinear);
   assert( nhandleconss <= nactiveconss );
   if ( nhandleconss < nactiveconss )
   {
      /* In this case we found unkown constraints and we exit, since we cannot handle them. */
      *success = FALSE;
      *nperms = -1;
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "Detecting %ssymmetry on %d variables and %d constraints.\n", local ? "local " : "", nvars, nactiveconss);

   /* copy variables */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vars, SCIPgetVars(scip), nvars) ); /*lint !e666*/
   assert( vars != NULL );

   /* fill matrixdata */

   /* use a staggered scheme for allocating space for non-zeros of constraint matrix since it can become large */
   if ( nvars <= 100000 )
      matrixdata.nmaxmatcoef = 100 * nvars;
   else if ( nvars <= 1000000 )
      matrixdata.nmaxmatcoef = 32 * nvars;
   else if ( nvars <= 16700000 )
      matrixdata.nmaxmatcoef = 16 * nvars;
   else
      matrixdata.nmaxmatcoef = INT_MAX / 10;

   matrixdata.nmatcoef = 0;
   matrixdata.nrhscoef = 0;
   matrixdata.nuniquemat = 0;
   matrixdata.nuniquevars = 0;
   matrixdata.nuniquerhs = 0;
   matrixdata.npermvars = nvars;
   matrixdata.permvars = vars;
   matrixdata.permvarcolors = NULL;
   matrixdata.matcoefcolors = NULL;
   matrixdata.rhscoefcolors = NULL;

   /* fill exprdata */
   exprdata.nuniqueoperators = 0;
   exprdata.nuniquecoefs = 0;
   exprdata.nuniqueconstants = 0;

   /* prepare matrix data (use block memory, since this can become large) */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhscoef, 2 * nactiveconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhssense, 2 * nactiveconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhsidx, 2 * nactiveconss) );

   /* prepare temporary constraint data (use block memory, since this can become large);
    * also allocate memory for fixed vars since some vars might have been deactivated meanwhile */
   nallvars = nvars + SCIPgetNFixedVars(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consvars, nallvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consvals, nallvars) );

   /* create hashset for auxvars and iterator for nonlinear constraints */
   if( nnlconss > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, isnonlinvar, nvars) );
      SCIP_CALL( SCIPhashsetCreate(&auxvars, SCIPblkmem(scip), nnlconss) );
      SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   }
   else
      *isnonlinvar = NULL;

   /* allocate memory for getting the number of constraints that contain a variable */
   if ( usecolumnsparsity )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &nconssforvar, nvars) );
   }

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      const char* conshdlrname;
      SCIP_CONS* cons;
      SCIP_VAR** linvars;
      int nconsvars;

      /* get constraint */
      cons = conss[c];
      assert( cons != NULL );

      /* skip non-active constraints */
      if ( ! SCIPconsIsActive(cons) )
         continue;

      /* Skip conflict constraints if we are late in the solving process */
      if ( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPconsIsConflict(cons) )
         continue;

      /* get constraint handler */
      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( conshdlrname != NULL );

      /* check type of constraint */
      if ( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( collectCoefficients(scip, doubleequations, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
               SCIPgetNVarsLinear(scip, cons), SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons),
               SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "linking") == 0 )
      {
         SCIP_VAR** curconsvars;
         SCIP_Real* curconsvals;
         int i;

         /* get constraint variables and their coefficients */
         curconsvals = SCIPgetValsLinking(scip, cons);
         SCIP_CALL( SCIPgetBinvarsLinking(scip, cons, &curconsvars, &nconsvars) );
         /* SCIPgetBinVarsLinking returns the number of binary variables, but we also need the integer variable */
         nconsvars++;

         /* copy vars and vals for binary variables */
         for (i = 0; i < nconsvars - 1; i++)
         {
            consvars[i] = curconsvars[i];
            consvals[i] = (SCIP_Real) curconsvals[i];
         }

         /* set final entry of vars and vals to the linking variable and its coefficient, respectively */
         consvars[nconsvars - 1] = SCIPgetLinkvarLinking(scip, cons);
         consvals[nconsvars - 1] = -1.0;

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nconsvars, 0.0, 0.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata, nconssforvar) );
         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, NULL, nconsvars - 1, 1.0, 1.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "setppc") == 0 )
      {
         linvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch ( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( collectCoefficients(scip, doubleequations, linvars, 0, nconsvars, 1.0, 1.0, SCIPconsIsTransformed(cons), SYM_SENSE_EQUATION, &matrixdata, nconssforvar) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( collectCoefficients(scip, doubleequations, linvars, 0, nconsvars, -SCIPinfinity(scip), 1.0, SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( collectCoefficients(scip, doubleequations, linvars, 0, nconsvars, 1.0, SCIPinfinity(scip), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
            break;
         default:
            SCIPerrorMessage("Unknown setppc type %d.\n", SCIPgetTypeSetppc(scip, cons));
            return SCIP_ERROR;
         }
      }
      else if ( strcmp(conshdlrname, "xor") == 0 )
      {
         SCIP_VAR** curconsvars;
         SCIP_VAR* var;

         /* get number of variables of XOR constraint (without integer variable) */
         nconsvars = SCIPgetNVarsXor(scip, cons);

         /* get variables of XOR constraint */
         curconsvars = SCIPgetVarsXor(scip, cons);
         for (j = 0; j < nconsvars; ++j)
         {
            assert( curconsvars[j] != NULL );
            consvars[j] = curconsvars[j];
            consvals[j] = 1.0;
         }

         /* intVar of xor constraint might have been removed */
         var = SCIPgetIntVarXor(scip, cons);
         if ( var != NULL )
         {
            consvars[nconsvars] = var;
            consvals[nconsvars++] = 2.0;
         }
         assert( nconsvars <= nallvars );

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nconsvars, (SCIP_Real) SCIPgetRhsXor(scip, cons),
               (SCIP_Real) SCIPgetRhsXor(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_XOR, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "and") == 0 )
      {
         SCIP_VAR** curconsvars;

         /* get number of variables of AND constraint (without resultant) */
         nconsvars = SCIPgetNVarsAnd(scip, cons);

         /* get variables of AND constraint */
         curconsvars = SCIPgetVarsAnd(scip, cons);

         for (j = 0; j < nconsvars; ++j)
         {
            assert( curconsvars[j] != NULL );
            consvars[j] = curconsvars[j];
            consvals[j] = 1.0;
         }

         assert( SCIPgetResultantAnd(scip, cons) != NULL );
         consvars[nconsvars] = SCIPgetResultantAnd(scip, cons);
         consvals[nconsvars++] = 2.0;
         assert( nconsvars <= nallvars );

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nconsvars, 0.0, 0.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_AND, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "or") == 0 )
      {
         SCIP_VAR** curconsvars;

         /* get number of variables of OR constraint (without resultant) */
         nconsvars = SCIPgetNVarsOr(scip, cons);

         /* get variables of OR constraint */
         curconsvars = SCIPgetVarsOr(scip, cons);

         for (j = 0; j < nconsvars; ++j)
         {
            assert( curconsvars[j] != NULL );
            consvars[j] = curconsvars[j];
            consvals[j] = 1.0;
         }

         assert( SCIPgetResultantOr(scip, cons) != NULL );
         consvars[nconsvars] = SCIPgetResultantOr(scip, cons);
         consvals[nconsvars++] = 2.0;
         assert( nconsvars <= nallvars );

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nconsvars, 0.0, 0.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_OR, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( collectCoefficients(scip, doubleequations, SCIPgetVarsLogicor(scip, cons), 0, SCIPgetNVarsLogicor(scip, cons),
               1.0, SCIPinfinity(scip), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;

         nconsvars = SCIPgetNVarsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array and get active variables of constraint */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         for (j = 0; j < nconsvars; ++j)
            consvals[j] = (SCIP_Real) weights[j];
         assert( nconsvars <= nallvars );

         SCIP_CALL( collectCoefficients(scip, doubleequations, SCIPgetVarsKnapsack(scip, cons), consvals, nconsvars, -SCIPinfinity(scip),
               (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "varbound") == 0 )
      {
         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvals[0] = 1.0;

         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
               SCIPgetRhsVarbound(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "bounddisjunction") == 0 )
      {
         /* To model bound disjunctions, we normalize each constraint
          * \f[
          *   (x_1 \{\leq,\geq\} b_1) \vee \ldots \vee (x_n \{\leq,\geq\} b_n)
          * \f]
          * to a constraint of type
          * \f[
          *   (x_1 \leq b'_1 \vee \ldots \vee (x_n \leq b'_n).
          * \f]
          *
          * If no variable appears twice in such a normalized constraint, we say this bound disjunction
          * is of type 1. If the bound disjunction has length two and both disjunctions contain the same variable,
          * we say the bound disjunction is of type 2. Further bound disjunctions are possible, but can currently
          * not be handled.
          *
          * Bound disjunctions of type 1 are modeled as the linear constraint
          * \f[
          *    b'_1 \cdot x_1 + \ldots +  b'_n \cdot x_n = 0
          * \f]
          * and bound disjunctions of type 2 are modeled as the linear constraint
          * \f[
          *    \min\{b'_1, b'_2\} \leq x_1 \leq \max\{b'_1, b'_2\}.
          * \f]
          * Note that problems arise if \fb'_i = 0\f for some variable \fx_i\f, because its coefficient in the
          * linear constraint is 0. To avoid this, we replace 0 by a special number.
          */
         SCIP_VAR** bounddisjvars;
         SCIP_BOUNDTYPE* boundtypes;
         SCIP_Real* bounds;
         SCIP_Bool repetition = FALSE;
         int nbounddisjvars;
         int k;

         /* collect coefficients for normalized constraint */
         nbounddisjvars = SCIPgetNVarsBounddisjunction(scip, cons);
         bounddisjvars = SCIPgetVarsBounddisjunction(scip, cons);
         boundtypes = SCIPgetBoundtypesBounddisjunction(scip, cons);
         bounds = SCIPgetBoundsBounddisjunction(scip, cons);

         /* copy data */
         for (j = 0; j < nbounddisjvars; ++j)
         {
            consvars[j] = bounddisjvars[j];

            /* normalize bounddisjunctions to SCIP_BOUNDTYPE_LOWER */
            if ( boundtypes[j] == SCIP_BOUNDTYPE_LOWER )
               consvals[j] = - bounds[j];
            else
               consvals[j] = bounds[j];

            /* special treatment of 0 values */
            if ( SCIPisZero(scip, consvals[j]) )
               consvals[j] = SCIP_SPECIALVAL;

            /* detect whether a variable appears in two literals */
            for (k = 0; k < j && ! repetition; ++k)
            {
               if ( consvars[j] == consvars[k] )
                  repetition = TRUE;
            }

            /* stop, we cannot handle bounddisjunctions with more than two variables that contain a variable twice */
            if ( repetition && nbounddisjvars > 2 )
            {
               *success = FALSE;

               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                  "   Deactivated symmetry handling methods, there exist constraints that cannot be handled by symmetry methods.\n");

               if ( usecolumnsparsity )
                  SCIPfreeBlockMemoryArrayNull(scip, &nconssforvar, nvars);

               SCIPfreeBlockMemoryArrayNull(scip, &consvals, nallvars);
               SCIPfreeBlockMemoryArrayNull(scip, &consvars, nallvars);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);
               SCIPfreeBlockMemoryArrayNull(scip, &vars, nvars);

               return SCIP_OKAY;
            }
         }
         assert( ! repetition || nbounddisjvars == 2 );

         /* if no variable appears twice */
         if ( ! repetition )
         {
            /* add information for bounddisjunction of type 1 */
            SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nbounddisjvars, 0.0, 0.0,
                  SCIPconsIsTransformed(cons), SYM_SENSE_BOUNDIS_TYPE_1, &matrixdata, nconssforvar) );
         }
         else
         {
            /* add information for bounddisjunction of type 2 */
            SCIP_Real lhs;
            SCIP_Real rhs;

            lhs = MIN(consvals[0], consvals[1]);
            rhs = MAX(consvals[0], consvals[1]);

            consvals[0] = 1.0;

            SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, 1, lhs, rhs,
                  SCIPconsIsTransformed(cons), SYM_SENSE_BOUNDIS_TYPE_2, &matrixdata, nconssforvar) );
         }
      }
      else if ( strcmp(conshdlrname, "nonlinear") == 0 )
      {
         SCIP_EXPR* expr;
         SCIP_EXPR* rootexpr;

         rootexpr = SCIPgetExprNonlinear(cons);
         assert(rootexpr != NULL);

         /* for nonlinear constraints, only collect auxiliary variables for now */
         SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
         SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR);

         for (expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it)) /*lint !e441*/ /*lint !e440*/
         {
            assert( SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_ENTEREXPR );

            /* for variables, we check whether they appear nonlinearly and store the result in the resp. array */
            if ( SCIPisExprVar(scip, expr) )
            {
               assert(*isnonlinvar != NULL);
               (*isnonlinvar)[SCIPvarGetProbindex(SCIPgetVarExprVar(expr))] = (SCIPexpriterGetParentDFS(it) != rootexpr || !SCIPisExprSum(scip, rootexpr));
            }
            else
            {
               SCIP_VAR* auxvar = SCIPgetExprAuxVarNonlinear(expr);

               if ( auxvar != NULL && !SCIPhashsetExists(auxvars, (void*) auxvar) )
               {
                  SCIP_CALL( SCIPhashsetInsert(auxvars, SCIPblkmem(scip), (void*) auxvar) );
               }

               if ( SCIPisExprValue(scip, expr) )
                  ++exprdata.nuniqueconstants;
               else if ( SCIPisExprSum(scip, expr) )
               {
                  ++exprdata.nuniqueoperators;
                  ++exprdata.nuniqueconstants;
                  exprdata.nuniquecoefs += SCIPexprGetNChildren(expr);
               }
               else
                  ++exprdata.nuniqueoperators;
            }
         }
      }
      else
      {
         /* if constraint is not one of the previous types, it cannot be handled */
         SCIPerrorMessage("Cannot determine symmetries for constraint <%s> of constraint handler <%s>.\n",
            SCIPconsGetName(cons), SCIPconshdlrGetName(conshdlr) );
         return SCIP_ERROR;
      }
   }
   assert( matrixdata.nrhscoef <= 2 * (nactiveconss - nnlconss) );
   assert( matrixdata.nrhscoef >= 0 );

   SCIPfreeBlockMemoryArray(scip, &consvals, nallvars);
   SCIPfreeBlockMemoryArray(scip, &consvars, nallvars);

   /* if no active constraint contains active variables */
   if ( nnlconss == 0 && matrixdata.nrhscoef == 0 )
   {
      *success = TRUE;

      if ( usecolumnsparsity )
         SCIPfreeBlockMemoryArrayNull(scip, &nconssforvar, nvars);

      /* free matrix data */
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);

      SCIPfreeBlockMemoryArray(scip, &vars, nvars);

      return SCIP_OKAY;
   }

   /* sort matrix coefficients (leave matrix array intact) */
   SCIPsort(matrixdata.matidx, SYMsortMatCoef, (void*) matrixdata.matcoef, matrixdata.nmatcoef);

   /* sort rhs types (first by sense, then by value, leave rhscoef intact) */
   sortrhstype.vals = matrixdata.rhscoef;
   sortrhstype.senses = matrixdata.rhssense;
   sortrhstype.nrhscoef = matrixdata.nrhscoef;
   SCIPsort(matrixdata.rhsidx, SYMsortRhsTypes, (void*) &sortrhstype, matrixdata.nrhscoef);

   /* create map for variables to indices */
   SCIP_CALL( SCIPhashtableCreate(&vartypemap, SCIPblkmem(scip), 5 * nvars, SYMhashGetKeyVartype, SYMhashKeyEQVartype, SYMhashKeyValVartype, (void*) scip) );
   assert( vartypemap != NULL );

   /* allocate space for mappings to colors */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.permvarcolors, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matcoefcolors, matrixdata.nmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhscoefcolors, matrixdata.nrhscoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &uniquevararray, nvars) );

   /* determine number of different coefficients */

   /* find non-equivalent variables: same objective, lower and upper bounds, and variable type */
   for (j = 0; j < nvars; ++j)
   {
      SCIP_VAR* var;

      var = vars[j];
      assert( var != NULL );

      /* if the variable type should be fixed, just increase the color */
      if ( SymmetryFixVar(fixedtype, var) || (nnlconss > 0 && SCIPhashsetExists(auxvars, (void*) var)) )
      {
         matrixdata.permvarcolors[j] = matrixdata.nuniquevars++;
#ifdef SCIP_OUTPUT
         SCIPdebugMsg(scip, "Detected variable <%s> of fixed type %d - color %d.\n", SCIPvarGetName(var), SCIPvarGetType(var), matrixdata.nuniquevars - 1);
#endif
      }
      else
      {
         SYM_VARTYPE* vt;

         vt = &uniquevararray[nuniquevararray];
         assert( nuniquevararray <= matrixdata.nuniquevars );

         vt->obj = SCIPvarGetObj(var);
         if ( local )
         {
            vt->lb = SCIPvarGetLbLocal(var);
            vt->ub = SCIPvarGetUbLocal(var);
         }
         else
         {
            vt->lb = SCIPvarGetLbGlobal(var);
            vt->ub = SCIPvarGetUbGlobal(var);
         }
         vt->type = SCIPvarGetType(var);
         vt->nconss = usecolumnsparsity ? nconssforvar[j] : 0; /*lint !e613*/

         if ( ! SCIPhashtableExists(vartypemap, (void*) vt) )
         {
            SCIP_CALL( SCIPhashtableInsert(vartypemap, (void*) vt) );
            vt->color = matrixdata.nuniquevars;
            matrixdata.permvarcolors[j] = matrixdata.nuniquevars++;
            ++nuniquevararray;
#ifdef SCIP_OUTPUT
            SCIPdebugMsg(scip, "Detected variable <%s> of new type (probindex: %d, obj: %g, lb: %g, ub: %g, type: %d) - color %d.\n",
               SCIPvarGetName(var), SCIPvarGetProbindex(var), vt->obj, vt->lb, vt->ub, vt->type, matrixdata.nuniquevars - 1);
#endif
         }
         else
         {
            SYM_VARTYPE* vtr;

            vtr = (SYM_VARTYPE*) SCIPhashtableRetrieve(vartypemap, (void*) vt);
            matrixdata.permvarcolors[j] = vtr->color;
         }
      }
   }

   /* If every variable is unique, terminate. -> no symmetries can be present */
   if ( matrixdata.nuniquevars == nvars )
   {
      *success = TRUE;

      /* free matrix data */
      SCIPfreeBlockMemoryArray(scip, &uniquevararray, nvars);

      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoefcolors, matrixdata.nrhscoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoefcolors, matrixdata.nmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.permvarcolors, nvars);
      SCIPhashtableFree(&vartypemap);

      if ( usecolumnsparsity )
         SCIPfreeBlockMemoryArrayNull(scip, &nconssforvar, nvars);

      if ( nnlconss > 0 )
      {
         SCIPfreeExpriter(&it);
         SCIPhashsetFree(&auxvars, SCIPblkmem(scip));
         SCIPfreeBlockMemoryArrayNull(scip, isnonlinvar, nvars);
      }

      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);

      SCIPfreeBlockMemoryArray(scip, &vars, nvars);

      return SCIP_OKAY;
   }

   /* find non-equivalent matrix entries (use sorting to avoid too many map calls) */
   for (j = 0; j < matrixdata.nmatcoef; ++j)
   {
      int idx;

      idx = matrixdata.matidx[j];
      assert( 0 <= idx && idx < matrixdata.nmatcoef );

      val = matrixdata.matcoef[idx];
      assert( oldcoef == SCIP_INVALID || oldcoef <= val ); /*lint !e777*/

      if ( ! SCIPisEQ(scip, val, oldcoef) )
      {
#ifdef SCIP_OUTPUT
         SCIPdebugMsg(scip, "Detected new matrix entry type %f - color: %d\n.", val, matrixdata.nuniquemat);
#endif
         matrixdata.matcoefcolors[idx] = matrixdata.nuniquemat++;
         oldcoef = val;
      }
      else
      {
         assert( matrixdata.nuniquemat > 0 );
         matrixdata.matcoefcolors[idx] = matrixdata.nuniquemat - 1;
      }
   }

   /* find non-equivalent rhs */
   oldcoef = SCIP_INVALID;
   for (j = 0; j < matrixdata.nrhscoef; ++j)
   {
      SYM_RHSSENSE sense;
      int idx;

      idx = matrixdata.rhsidx[j];
      assert( 0 <= idx && idx < matrixdata.nrhscoef );
      sense = matrixdata.rhssense[idx];
      val = matrixdata.rhscoef[idx];

      /* make sure that new senses are treated with new color */
      if ( sense != oldsense )
         oldcoef = SCIP_INVALID;
      oldsense = sense;
      assert( oldcoef == SCIP_INVALID || oldcoef <= val ); /*lint !e777*/

      /* assign new color to new type */
      if ( ! SCIPisEQ(scip, val, oldcoef) )
      {
#ifdef SCIP_OUTPUT
         SCIPdebugMsg(scip, "Detected new rhs type %f, type: %u - color: %d\n", val, sense, matrixdata.nuniquerhs);
#endif
         matrixdata.rhscoefcolors[idx] = matrixdata.nuniquerhs++;
         oldcoef = val;
      }
      else
      {
         assert( matrixdata.nuniquerhs > 0 );
         matrixdata.rhscoefcolors[idx] = matrixdata.nuniquerhs - 1;
      }
   }
   assert( 0 < matrixdata.nuniquevars && matrixdata.nuniquevars <= nvars );
   assert( 0 <= matrixdata.nuniquerhs && matrixdata.nuniquerhs <= matrixdata.nrhscoef );
   assert( 0 <= matrixdata.nuniquemat && matrixdata.nuniquemat <= matrixdata.nmatcoef );

   SCIPdebugMsg(scip, "Number of detected different variables: %d (total: %d).\n", matrixdata.nuniquevars, nvars);
   SCIPdebugMsg(scip, "Number of detected different rhs types: %d (total: %d).\n", matrixdata.nuniquerhs, matrixdata.nrhscoef);
   SCIPdebugMsg(scip, "Number of detected different matrix coefficients: %d (total: %d).\n", matrixdata.nuniquemat, matrixdata.nmatcoef);

   /* do not compute symmetry if all variables are non-equivalent (unique) or if all matrix coefficients are different */
   if ( matrixdata.nuniquevars < nvars && (matrixdata.nuniquemat == 0 || matrixdata.nuniquemat < matrixdata.nmatcoef) )
   {
      /* determine generators */
      SCIP_CALL( SYMcomputeSymmetryGenerators(scip, maxgenerators, &matrixdata, &exprdata, nperms, nmaxperms,
            perms, log10groupsize) );
      assert( *nperms <= *nmaxperms );

      /* SCIPisStopped() might call SCIPgetGap() which is only available after initpresolve */
      if ( checksymmetries && SCIPgetStage(scip) > SCIP_STAGE_INITPRESOLVE && ! SCIPisStopped(scip) )
      {
         SCIP_CALL( checkSymmetriesAreSymmetries(scip, fixedtype, &matrixdata, *nperms, *perms) );
      }

      if ( *nperms > 0 )
      {
         SCIP_CALL( setSymmetryData(scip, vars, nvars, nbinvars, permvars, npermvars, nbinpermvars, *perms, *nperms,
               nmovedvars, binvaraffected, compresssymmetries, compressthreshold, compressed) );
      }
      else
      {
         SCIPfreeBlockMemoryArrayNull(scip, isnonlinvar, nvars);
         SCIPfreeBlockMemoryArray(scip, &vars, nvars);
      }
   }
   else
   {
      SCIPfreeBlockMemoryArrayNull(scip, isnonlinvar, nvars);
      SCIPfreeBlockMemoryArray(scip, &vars, nvars);
   }
   *success = TRUE;

   /* free matrix data */
   SCIPfreeBlockMemoryArray(scip, &uniquevararray, nvarsorig);

   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoefcolors, matrixdata.nrhscoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoefcolors, matrixdata.nmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.permvarcolors, nvarsorig);
   SCIPhashtableFree(&vartypemap);

   if ( usecolumnsparsity )
      SCIPfreeBlockMemoryArrayNull(scip, &nconssforvar, nvarsorig);

   /* free cons expr specific data */
   if ( nnlconss > 0 )
   {
      assert( it != NULL );
      assert( auxvars != NULL );

      SCIPfreeExpriter(&it);
      SCIPhashsetFree(&auxvars, SCIPblkmem(scip));
   }

   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);

   return SCIP_OKAY;
}


/** determines symmetry */
static
SCIP_RETCODE determineSymmetry(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SYM_SPEC              symspecrequire,     /**< symmetry specification for which we need to compute symmetries */
   SYM_SPEC              symspecrequirefixed /**< symmetry specification of variables which must be fixed by symmetries */
   )
{ /*lint --e{641}*/
   SCIP_Bool successful;
   int maxgenerators;
   int nhandleconss;
   int nconss;
   unsigned int type = 0;
   int nvars;
   int j;
   int p;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );
   assert( propdata->ofenabled || propdata->symconsenabled || propdata->sstenabled );

   /* do not compute symmetry if reoptimization is enabled */
   if ( SCIPisReoptEnabled(scip) )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;
      return SCIP_OKAY;
   }

   /* do not compute symmetry if Benders decomposition enabled */
   if ( SCIPgetNActiveBenders(scip) > 0 )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;
      return SCIP_OKAY;
   }

   /* skip symmetry computation if no graph automorphism code was linked */
   if ( ! SYMcanComputeSymmetry() )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;

      nconss = SCIPgetNActiveConss(scip);
      nhandleconss = getNSymhandableConss(scip, propdata->conshdlr_nonlinear);

      /* print verbMessage only if problem consists of symmetry handable constraints */
      assert( nhandleconss <=  nconss );
      if ( nhandleconss < nconss )
         return SCIP_OKAY;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   Deactivated symmetry handling methods, since SCIP was built without symmetry detector (SYM=none).\n");

      return SCIP_OKAY;
   }

   /* do not compute symmetry if there are active pricers */
   if ( SCIPgetNActivePricers(scip) > 0 )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;

      return SCIP_OKAY;
   }

   /* avoid trivial cases */
   nvars = SCIPgetNVars(scip);
   if ( nvars <= 0 )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;

      return SCIP_OKAY;
   }

   /* do not compute symmetry if there are no binary variables and non-binary variables cannot be handled, but only binary variables would be used */
   if ( propdata->onlybinarysymmetry && SCIPgetNBinVars(scip) == 0 )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      /* terminate if only Schreier Sims for binary variables is selected */
      if ( propdata->sstenabled )
      {
         if ( ! ((ISSSTINTACTIVE(propdata->sstleadervartype) && SCIPgetNIntVars(scip) > 0)
               || (ISSSTIMPLINTACTIVE(propdata->sstleadervartype) && SCIPgetNImplVars(scip) > 0)
               || (ISSSTCONTACTIVE(propdata->sstleadervartype) && SCIPgetNContVars(scip) > 0)) )
            return SCIP_OKAY;
      }
      else
         return SCIP_OKAY;
   }

   /* determine symmetry specification */
   if ( SCIPgetNBinVars(scip) > 0 )
      type |= (int) SYM_SPEC_BINARY;
   if ( SCIPgetNIntVars(scip) > 0 )
      type |= (int) SYM_SPEC_INTEGER;
   /* count implicit integer variables as real variables, since we cannot currently handle integral variables well */
   if ( SCIPgetNContVars(scip) > 0 || SCIPgetNImplVars(scip) > 0 )
      type |= (int) SYM_SPEC_REAL;

   /* skip symmetry computation if required variables are not present */
   if ( ! (type & symspecrequire) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   (%.1fs) symmetry computation skipped: type (bin %c, int %c, cont %c) does not match requirements (bin %c, int %c, cont %c).\n",
         SCIPgetSolvingTime(scip),
         SCIPgetNBinVars(scip) > 0 ? '+' : '-',
         SCIPgetNIntVars(scip) > 0  ? '+' : '-',
         SCIPgetNContVars(scip) + SCIPgetNImplVars(scip) > 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_REAL) != 0 ? '+' : '-');

      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;

      return SCIP_OKAY;
   }

   /* if a restart occured, possibly prepare symmetry data to be recomputed */
   if ( isSymmetryRecomputationRequired(scip, propdata) )
   {
      /* reset symmetry information */
      SCIP_CALL( delSymConss(scip, propdata) );
      SCIP_CALL( freeSymmetryData(scip, propdata) );

      propdata->lastrestart = SCIPgetNRuns(scip);
      propdata->offoundreduction = FALSE;
   }

   /* skip computation if symmetry has already been computed */
   if ( propdata->computedsymmetry )
      return SCIP_OKAY;

   /* skip symmetry computation if there are constraints that cannot be handled by symmetry */
   nconss = SCIPgetNActiveConss(scip);
   nhandleconss = getNSymhandableConss(scip, propdata->conshdlr_nonlinear);
   if ( nhandleconss < nconss )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   (%.1fs) symmetry computation skipped: there exist constraints that cannot be handled by symmetry methods.\n",
         SCIPgetSolvingTime(scip));

      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;

      return SCIP_OKAY;
   }

   assert( propdata->npermvars == 0 );
   assert( propdata->permvars == NULL );
#ifndef NDEBUG
   assert( propdata->permvarsobj == NULL );
#endif
   assert( propdata->nperms < 0 );
   assert( propdata->nmaxperms == 0 );
   assert( propdata->perms == NULL );

   /* output message */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
      "   (%.1fs) symmetry computation started: requiring (bin %c, int %c, cont %c), (fixed: bin %c, int %c, cont %c)\n",
      SCIPgetSolvingTime(scip),
      (symspecrequire & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
      (symspecrequire & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
      (symspecrequire & (int) SYM_SPEC_REAL) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_REAL) != 0 ? '+' : '-');

   /* output warning if we want to fix certain symmetry parts that we also want to compute */
   if ( symspecrequire & symspecrequirefixed )
      SCIPwarningMessage(scip, "Warning: some required symmetries must be fixed.\n");

   /* determine maximal number of generators depending on the number of variables */
   maxgenerators = propdata->maxgenerators;
   maxgenerators = MIN(maxgenerators, MAXGENNUMERATOR / nvars);

   /* actually compute (global) symmetry */
   SCIP_CALL( computeSymmetryGroup(scip, propdata->doubleequations, propdata->compresssymmetries, propdata->compressthreshold,
	 maxgenerators, symspecrequirefixed, FALSE, propdata->checksymmetries, propdata->usecolumnsparsity, propdata->conshdlr_nonlinear,
         &propdata->npermvars, &propdata->nbinpermvars, &propdata->permvars, &propdata->nperms, &propdata->nmaxperms,
         &propdata->perms, &propdata->log10groupsize, &propdata->nmovedvars, &propdata->isnonlinvar,
         &propdata->binvaraffected, &propdata->compressed, &successful) );

   /* mark that we have computed the symmetry group */
   propdata->computedsymmetry = TRUE;

   /* store restart level */
   propdata->lastrestart = SCIPgetNRuns(scip);

   /* return if not successful */
   if ( ! successful )
   {
      assert( checkSymmetryDataFree(propdata) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) could not compute symmetry\n", SCIPgetSolvingTime(scip));

      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;

      return SCIP_OKAY;
   }

   /* return if no symmetries found */
   if ( propdata->nperms == 0 )
   {
      assert( checkSymmetryDataFree(propdata) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) no symmetry present\n", SCIPgetSolvingTime(scip));

      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      propdata->sstenabled = FALSE;

      return SCIP_OKAY;
   }

   /* display statistics */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) symmetry computation finished: %d generators found (max: ",
      SCIPgetSolvingTime(scip), propdata->nperms);

   /* display statistics: maximum number of generators */
   if ( maxgenerators == 0 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "-");
   else
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%d", maxgenerators);

   /* display statistics: log10 group size, number of affected vars*/
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", log10 of symmetry group size: %.1f", propdata->log10groupsize);

   if ( propdata->displaynorbitvars )
   {
      if ( propdata->nmovedvars == -1 )
      {
         SCIP_CALL( SCIPdetermineNVarsAffectedSym(scip, propdata->perms, propdata->nperms, propdata->permvars,
               propdata->npermvars, &(propdata->nmovedvars)) );
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", number of affected variables: %d)\n", propdata->nmovedvars);
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ")\n");

   /* exit if no binary variables are affected by symmetry and we cannot handle non-binary symmetries */
   if ( ! propdata->binvaraffected )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) no symmetry on binary variables present.\n", SCIPgetSolvingTime(scip));

      /* disable OF and symmetry handling constraints based on symretopes */
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      if ( ! propdata->sstenabled ||
         ! ( ISSSTINTACTIVE(propdata->sstleadervartype) || ISSSTIMPLINTACTIVE(propdata->sstleadervartype)
            || ISSSTCONTACTIVE(propdata->sstleadervartype) ) )
      {
         propdata->sstenabled = FALSE;

         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) -> no handable symmetry found, free symmetry data.\n",
            SCIPgetSolvingTime(scip));

         /* free data and exit */
         SCIP_CALL( freeSymmetryData(scip, propdata) );

         return SCIP_OKAY;
      }
   }

   assert( propdata->nperms > 0 );
   assert( propdata->npermvars > 0 );

   /* compute components */
   assert( propdata->components == NULL );
   assert( propdata->componentbegins == NULL );
   assert( propdata->vartocomponent == NULL );

#ifdef SCIP_OUTPUT_COMPONENT
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) component computation started\n", SCIPgetSolvingTime(scip));
#endif

   /* we only need the components for orbital fixing, orbitope and subgroup detection, and Schreier Sims constraints */
   if ( propdata->ofenabled || ( propdata->symconsenabled && propdata->detectorbitopes )
      || propdata->detectsubgroups || propdata->sstenabled )
   {
      SCIP_CALL( SCIPcomputeComponentsSym(scip, propdata->perms, propdata->nperms, propdata->permvars,
            propdata->npermvars, FALSE, &propdata->components, &propdata->componentbegins,
            &propdata->vartocomponent, &propdata->componentblocked, &propdata->ncomponents) );
   }

#ifdef SCIP_OUTPUT_COMPONENT
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) component computation finished\n", SCIPgetSolvingTime(scip));
#endif

   /* set up data for OF */
   if ( propdata->ofenabled )
   {
      int componentidx;
      int v;

      /* transpose symmetries matrix here if necessary */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->permstrans, propdata->npermvars) );
      for (v = 0; v < propdata->npermvars; ++v)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->permstrans[v]), propdata->nmaxperms) );
         for (p = 0; p < propdata->nperms; ++p)
         {
            if ( SCIPvarIsBinary(propdata->permvars[v]) || propdata->sstenabled )
               propdata->permstrans[v][p] = propdata->perms[p][v];
            else
               propdata->permstrans[v][p] = v; /* ignore symmetry information on non-binary variables */
         }
      }

      /* prepare array for active permutations */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->inactiveperms, propdata->nperms) );
      for (v = 0; v < propdata->nperms; ++v)
         propdata->inactiveperms[v] = FALSE;

      /* create hashmap for storing the indices of variables */
      assert( propdata->permvarmap == NULL );
      SCIP_CALL( SCIPhashmapCreate(&propdata->permvarmap, SCIPblkmem(scip), propdata->npermvars) );

      /* prepare data structures */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->permvarsevents, propdata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->bg0, propdata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->bg0list, propdata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->bg1, propdata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->bg1list, propdata->npermvars) );

      /* insert variables into hashmap  */
      assert( propdata->nmovedpermvars == -1 );
      propdata->nmovedpermvars = 0;
      for (v = 0; v < propdata->npermvars; ++v)
      {
         SCIP_CALL( SCIPhashmapInsertInt(propdata->permvarmap, propdata->permvars[v], v) );

         propdata->bg0[v] = FALSE;
         propdata->bg1[v] = FALSE;
         propdata->permvarsevents[v] = -1;

         /* collect number of moved permvars */
         componentidx = propdata->vartocomponent[v];
         if ( componentidx > -1 && ! propdata->componentblocked[componentidx] )
         {
            propdata->nmovedpermvars += 1;

            if ( SCIPvarGetType(propdata->permvars[v]) == SCIP_VARTYPE_BINARY )
               ++propdata->nmovedbinpermvars;
            else if ( SCIPvarGetType(propdata->permvars[v]) == SCIP_VARTYPE_INTEGER )
               ++propdata->nmovedintpermvars;
            else if ( SCIPvarGetType(propdata->permvars[v]) == SCIP_VARTYPE_IMPLINT )
               ++propdata->nmovedimplintpermvars;
            else
               ++propdata->nmovedcontpermvars;
         }

         /* Only catch binary variables, since integer variables should be fixed pointwise; implicit integer variables
          * are not branched on. */
         if ( SCIPvarGetType(propdata->permvars[v]) == SCIP_VARTYPE_BINARY )
         {
            /* catch whether binary variables are globally fixed; also store filter position */
            SCIP_CALL( SCIPcatchVarEvent(scip, propdata->permvars[v], SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
                  propdata->eventhdlr, (SCIP_EVENTDATA*) propdata, &propdata->permvarsevents[v]) );
         }
      }
      assert( propdata->nbg1 == 0 );
   }

   /* set up data for Schreier Sims constraints or subgroup detection */
   if ( (propdata->sstenabled || propdata->detectsubgroups) && ! propdata->ofenabled )
   {
      int v;

      /* transpose symmetries matrix here if necessary */
      assert( propdata->permstrans == NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->permstrans, propdata->npermvars) );
      for (v = 0; v < propdata->npermvars; ++v)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->permstrans[v]), propdata->nmaxperms) );
         for (p = 0; p < propdata->nperms; ++p)
            propdata->permstrans[v][p] = propdata->perms[p][v];
      }

      /* create hashmap for storing the indices of variables */
      assert( propdata->permvarmap == NULL );
      SCIP_CALL( SCIPhashmapCreate(&propdata->permvarmap, SCIPblkmem(scip), propdata->npermvars) );

      /* insert variables into hashmap  */
      for (v = 0; v < propdata->npermvars; ++v)
      {
         SCIP_CALL( SCIPhashmapInsertInt(propdata->permvarmap, propdata->permvars[v], v) );
      }
   }

   /* handle several general aspects */
#ifndef NDEBUG
   /* store objective coefficients for debug purposes */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->permvarsobj, propdata->npermvars) );
   for (j = 0; j < propdata->npermvars; ++j)
      propdata->permvarsobj[j] = SCIPvarGetObj(propdata->permvars[j]);
#endif

   /* capture symmetric variables and forbid multi aggregation */

   /* binary symmetries are always handled
    *
    * note: binary variables are in the beginning of permvars
    */
   if ( propdata->binvaraffected )
   {
      for (j = 0; j < propdata->nbinpermvars; ++j)
      {
         SCIP_CALL( SCIPcaptureVar(scip, propdata->permvars[j]) );

         if ( propdata->compressed )
         {
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, propdata->permvars[j]) );
         }
         else
         {
            for (p = 0; p < propdata->nperms; ++p)
            {
               if ( propdata->perms[p][j] != j )
               {
                  SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, propdata->permvars[j]) );
                  break;
               }
            }
         }
      }
   }

   /* if Schreier-Sims constraints are enabled, also capture symmetric variables and forbid multi aggregation of handable vars */
   if ( propdata->sstenabled && propdata->sstleadervartype != (int) SCIP_SSTTYPE_BINARY )
   {
      int cnt;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->nonbinpermvarcaptured,
            propdata->npermvars - propdata->nbinpermvars) );
      for (j = propdata->nbinpermvars, cnt = 0; j < propdata->npermvars; ++j, ++cnt)
      {
         if ( ! isLeadervartypeCompatible(propdata->permvars[j], propdata->sstleadervartype) )
         {
            propdata->nonbinpermvarcaptured[cnt] = FALSE;
            continue;
         }

         SCIP_CALL( SCIPcaptureVar(scip, propdata->permvars[j]) );
         propdata->nonbinpermvarcaptured[cnt] = TRUE;

         if ( propdata->compressed )
         {
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, propdata->permvars[j]) );
         }
         else
         {
            for (p = 0; p < propdata->nperms; ++p)
            {
               if ( propdata->perms[p][j] != j )
               {
                  SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, propdata->permvars[j]) );
                  break;
               }
            }
         }
      }
   }

   /* free original perms matrix if no symmetry constraints are added */
   if ( ! propdata->symconsenabled && ! propdata->sstenabled )
   {
      for (p = 0; p < propdata->nperms; ++p)
      {
         SCIPfreeBlockMemoryArray(scip, &(propdata->perms)[p], propdata->npermvars);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &propdata->perms, propdata->nmaxperms);
   }

   return SCIP_OKAY;
}


/*
 * Functions for symmetry constraints
 */


/** Checks whether given set of 2-cycle permutations forms an orbitope and if so, builds the variable index matrix.
 *
 *  If @p activevars == NULL, then the function assumes all permutations of the component are active and therefore all
 *  moved vars are considered.
 *
 *  We need to keep track of the number of generatored columns, because we might not be able to detect all orbitopes.
 *  For example (1,2), (2,3), (3,4), (3,5) defines the symmetric group on {1,2,3,4,5}, but the generators we expect
 *  in our construction need shape (1,2), (2,3), (3,4), (4,5).
 *
 *  @pre @p orbitopevaridx has to be an initialized 2D array of size @p ntwocycles x @p nperms
 *  @pre @p columnorder has to be an initialized array of size nperms
 *  @pre @p nusedelems has to be an initialized array of size npermvars
 */
static
SCIP_RETCODE checkTwoCyclePermsAreOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< array of all permutation variables */
   int                   npermvars,          /**< number of permutation variables */
   int**                 perms,              /**< array of all permutations of the symmety group */
   int*                  activeperms,        /**< indices of the relevant permutations in perms */
   int                   ntwocycles,         /**< number of 2-cycles in the permutations */
   int                   nactiveperms,       /**< number of active permutations */
   int**                 orbitopevaridx,     /**< pointer to store variable index matrix */
   int*                  columnorder,        /**< pointer to store column order */
   int*                  nusedelems,         /**< pointer to store how often each element was used */
   int*                  nusedcols,          /**< pointer to store numer of columns used in orbitope (or NULL) */
   SCIP_Shortbool*       rowisbinary,        /**< pointer to store which rows are binary (or NULL) */
   SCIP_Bool*            isorbitope,         /**< buffer to store result */
   SCIP_Shortbool*       activevars          /**< bitset to store whether a variable is active (or NULL) */
   )
{  /*lint --e{571}*/
   SCIP_Bool* usedperm;
   SCIP_Bool foundperm = FALSE;
   int nusedperms = 0;
   int nfilledcols;
   int coltoextend;
   int ntestedperms = 0;
   int row = 0;
   int j;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( perms != NULL );
   assert( activeperms != NULL );
   assert( orbitopevaridx != NULL );
   assert( columnorder != NULL );
   assert( nusedelems != NULL );
   assert( isorbitope != NULL );
   assert( nactiveperms > 0 );
   assert( ntwocycles > 0 );
   assert( npermvars > 0 );
   assert( activevars == NULL || (0 <= nactiveperms && nactiveperms < npermvars) );

   *isorbitope = TRUE;
   if ( nusedcols != NULL )
      *nusedcols = 0;

   /* whether a permutation was considered to contribute to orbitope */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &usedperm, nactiveperms) );

   /* fill first two columns of orbitopevaridx matrix */

   /* look for the first active permutation which moves an active variable */
   while ( ! foundperm )
   {
      int permidx;

      assert( ntestedperms < nactiveperms );

      permidx = activeperms[ntestedperms];

      for (j = 0; j < npermvars; ++j)
      {
         if ( activevars != NULL && ! activevars[j] )
            continue;

         assert( activevars == NULL || activevars[perms[permidx][j]] );

         /* avoid adding the same 2-cycle twice */
         if ( perms[permidx][j] > j )
         {
            assert( SCIPvarIsBinary(permvars[j]) == SCIPvarIsBinary(permvars[perms[permidx][j]]) );

            if ( rowisbinary != NULL && SCIPvarIsBinary(permvars[j]) )
               rowisbinary[row] = TRUE;

            orbitopevaridx[row][0] = j;
            orbitopevaridx[row++][1] = perms[permidx][j];
            ++(nusedelems[j]);
            ++(nusedelems[perms[permidx][j]]);

            foundperm = TRUE;
         }

         if ( row == ntwocycles )
            break;
      }

      ++ntestedperms;
   }

   /* in the subgroup case it might happen that a generator has less than ntwocycles many 2-cyles */
   if ( row < ntwocycles )
   {
      *isorbitope = FALSE;
      SCIPfreeBufferArray(scip, &usedperm);
      return SCIP_OKAY;
   }

   usedperm[ntestedperms - 1] = TRUE;
   ++nusedperms;
   columnorder[0] = 0;
   columnorder[1] = 1;
   nfilledcols = 2;

   /* extend orbitopevaridx matrix to the left, i.e., iteratively find new permutations that
    * intersect the last added left column in each row in exactly one entry, starting with
    * column 0 */
   coltoextend = 0;
   for (j = ntestedperms; j < nactiveperms; ++j)
   {  /* lint --e{850} */
      SCIP_Bool success = FALSE;
      SCIP_Bool infeasible = FALSE;

      if ( nusedperms == nactiveperms )
         break;

      if ( usedperm[j] )
         continue;

      SCIP_CALL( SCIPextendSubOrbitope(orbitopevaridx, ntwocycles, nfilledcols, coltoextend,
            perms[activeperms[j]], TRUE, &nusedelems, permvars, NULL, &success, &infeasible) );

      if ( infeasible )
      {
         *isorbitope = FALSE;
         break;
      }
      else if ( success )
      {
         usedperm[j] = TRUE;
         ++nusedperms;
         coltoextend = nfilledcols;
         columnorder[nfilledcols++] = -1; /* mark column to be filled from the left */
         j = 0; /*lint !e850*/ /* reset j since previous permutations can now intersect with the latest added column */
      }
   }

   if ( ! *isorbitope ) /*lint !e850*/
   {
      SCIPfreeBufferArray(scip, &usedperm);
      return SCIP_OKAY;
   }

   coltoextend = 1;
   for (j = ntestedperms; j < nactiveperms; ++j)
   {  /*lint --e(850)*/
      SCIP_Bool success = FALSE;
      SCIP_Bool infeasible = FALSE;

      if ( nusedperms == nactiveperms )
         break;

      if ( usedperm[j] )
         continue;

      SCIP_CALL( SCIPextendSubOrbitope(orbitopevaridx, ntwocycles, nfilledcols, coltoextend,
            perms[activeperms[j]], FALSE, &nusedelems, permvars, NULL, &success, &infeasible) );

      if ( infeasible )
      {
         *isorbitope = FALSE;
         break;
      }
      else if ( success )
      {
         usedperm[j] = TRUE;
         ++nusedperms;
         coltoextend = nfilledcols;
         columnorder[nfilledcols] = 1; /* mark column to be filled from the right */
         ++nfilledcols;
         j = 0; /*lint !e850*/ /* reset j since previous permutations can now intersect with the latest added column */
      }
   }

   if ( activevars == NULL && nusedperms < nactiveperms ) /*lint !e850*/
      *isorbitope = FALSE;

   if ( nusedcols != NULL )
      *nusedcols = nfilledcols;
   assert( ! *isorbitope || activevars == NULL || nusedperms < nfilledcols );

   SCIPfreeBufferArray(scip, &usedperm);

   return SCIP_OKAY;
}

/** choose an order in which the generators should be added for subgroup detection */
static
SCIP_RETCODE chooseOrderOfGenerators(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int                   compidx,            /**< index of component */
   int**                 genorder,           /**< (initialized) buffer to store the resulting order of generator */
   int*                  ntwocycleperms      /**< pointer to store the number of 2-cycle permutations in component compidx */
   )
{
   int** perms;
   int* components;
   int* componentbegins;
   int* ntwocycles;
   int npermvars;
   int npermsincomp;
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( compidx >= 0 );
   assert( compidx < propdata->ncomponents );
   assert( genorder != NULL );
   assert( *genorder != NULL );
   assert( ntwocycleperms != NULL );
   assert( propdata->computedsymmetry );
   assert( propdata->nperms > 0 );
   assert( propdata->perms != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->ncomponents > 0 );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );

   perms = propdata->perms;
   npermvars = propdata->npermvars;
   components = propdata->components;
   componentbegins = propdata->componentbegins;
   npermsincomp = componentbegins[compidx + 1] - componentbegins[compidx];
   *ntwocycleperms = npermsincomp;

   SCIP_CALL( SCIPallocBufferArray(scip, &ntwocycles, npermsincomp) );

   for (i = 0; i < npermsincomp; ++i)
   {
      int* perm;
      int nbincycles;

      perm = perms[components[componentbegins[compidx] + i]];

      SCIP_CALL( SCIPisInvolutionPerm(perm, propdata->permvars, npermvars, &(ntwocycles[i]), &nbincycles, FALSE) );

      /* we skip permutations which do not purely consist of 2-cycles */
      if ( ntwocycles[i] == 0 )
      {
         /* we change the number of two cycles for this perm so that it will be sorted to the end */
         if ( propdata->preferlessrows )
            ntwocycles[i] = npermvars;
         else
            ntwocycles[i] = 0;
         --(*ntwocycleperms);
      }
      else if ( ! propdata->preferlessrows )
         ntwocycles[i] = - ntwocycles[i];
   }

   SCIPsortIntInt(ntwocycles, *genorder, npermsincomp);

   SCIPfreeBufferArray(scip, &ntwocycles);

   return SCIP_OKAY;
}


/** builds the graph for symmetric subgroup detection from the given permutation of generators
 *
 *  After execution, @p graphcomponents contains all permvars sorted by their color and component,
 *  @p graphcompbegins points to the indices where new components in @p graphcomponents start and
 *  @p compcolorbegins points to the indices where new colors in @p graphcompbegins start.
*/
static
SCIP_RETCODE buildSubgroupGraph(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  genorder,           /**< order in which the generators should be considered */
   int                   ntwocycleperms,     /**< number of 2-cycle permutations in this component */
   int                   compidx,            /**< index of the component */
   int**                 graphcomponents,    /**< buffer to store the components of the graph (ordered var indices) */
   int**                 graphcompbegins,    /**< buffer to store the indices of each new graph component */
   int**                 compcolorbegins,    /**< buffer to store at which indices a new color begins */
   int*                  ngraphcomponents,   /**< pointer to store the number of graph components */
   int*                  ncompcolors,        /**< pointer to store the number of different colors */
   int**                 usedperms,          /**< buffer to store the indices of permutations that were used */
   int*                  nusedperms,         /**< pointer to store the number of used permutations in the graph */
   int                   usedpermssize,      /**< initial size of usedperms */
   SCIP_Shortbool*       permused            /**< initialized buffer to store which permutations have been used
                                              *   (identified by index in component) */
   )
{
   SCIP_DISJOINTSET* vartocomponent;
   SCIP_DISJOINTSET* comptocolor;
   int** perms;
   int* components;
   int* componentbegins;
   int* componentslastperm;
   SYM_SORTGRAPHCOMPVARS graphcompvartype;
   int npermvars;
   int nextcolor;
   int nextcomp;
   int j;
   int k;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( graphcomponents != NULL );
   assert( graphcompbegins != NULL );
   assert( compcolorbegins != NULL );
   assert( ngraphcomponents != NULL );
   assert( ncompcolors != NULL );
   assert( genorder != NULL );
   assert( usedperms != NULL );
   assert( nusedperms != NULL );
   assert( usedpermssize > 0 );
   assert( permused != NULL );
   assert( ntwocycleperms >= 0 );
   assert( compidx >= 0 );
   assert( compidx < propdata->ncomponents );
   assert( propdata->computedsymmetry );
   assert( propdata->nperms > 0 );
   assert( propdata->perms != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->ncomponents > 0 );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );
   assert( !propdata->componentblocked[compidx] );

   perms = propdata->perms;
   npermvars = propdata->npermvars;
   components = propdata->components;
   componentbegins = propdata->componentbegins;
   *nusedperms = 0;

   assert( ntwocycleperms <= componentbegins[compidx + 1] - componentbegins[compidx] );

   SCIP_CALL( SCIPcreateDisjointset(scip, &vartocomponent, npermvars) );
   SCIP_CALL( SCIPcreateDisjointset(scip, &comptocolor, npermvars) );
   SCIP_CALL( SCIPallocBufferArray( scip, &componentslastperm, npermvars) );

   for (k = 0; k < npermvars; ++k)
      componentslastperm[k] = -1;

   for (j = 0; j < ntwocycleperms; ++j)
   {
      int* perm;
      int firstcolor = -1;

      /* use given order of generators */
      perm = perms[components[componentbegins[compidx] + genorder[j]]];
      assert( perm != NULL );

      /* iteratively handle each swap of perm until an invalid one is found or all edges have been added */
      for (k = 0; k < npermvars; ++k)
      {
         int comp1;
         int comp2;
         int color1;
         int color2;
         int img;

         img = perm[k];
         assert( perm[img] == k );

         if ( img <= k )
            continue;

         comp1 = SCIPdisjointsetFind(vartocomponent, k);
         comp2 = SCIPdisjointsetFind(vartocomponent, img);

         if ( comp1 == comp2 )
         {
            /* another permutation has already merged these variables into one component; store its color */
            if ( firstcolor < 0 )
            {
               assert( SCIPdisjointsetFind(comptocolor, comp1) == SCIPdisjointsetFind(comptocolor, comp2) );
               firstcolor = SCIPdisjointsetFind(comptocolor, comp1);
            }
            componentslastperm[comp1] = j;
            continue;
         }

         /* if it is the second time that the component is used for this generator,
          * it is not guaranteed that the group acts like the symmetric group, so skip it
          */
         if ( componentslastperm[comp1] == j || componentslastperm[comp2] == j )
            break;

         color1 = SCIPdisjointsetFind(comptocolor, comp1);
         color2 = SCIPdisjointsetFind(comptocolor, comp2);

         /* a generator is not allowed to connect two components of the same color, since they depend on each other */
         if ( color1 == color2 )
            break;

         componentslastperm[comp1] = j;
         componentslastperm[comp2] = j;

         if ( firstcolor < 0 )
            firstcolor = color1;
      }

      /* if the generator is invalid, delete the newly added edges, go to next generator */
      if ( k < npermvars )
         continue;

      /* if the generator only acts on already existing components, we don't have to store it */
      if ( firstcolor == -1 )
         continue;

      /* check whether we need to resize */
      if ( *nusedperms >= usedpermssize )
      {
         int newsize = SCIPcalcMemGrowSize(scip, (*nusedperms) + 1);
         assert( newsize > usedpermssize );

         SCIP_CALL( SCIPreallocBufferArray(scip, usedperms, newsize) );

         usedpermssize = newsize;
      }

      (*usedperms)[*nusedperms] = components[componentbegins[compidx] + genorder[j]];
      ++(*nusedperms);
      permused[genorder[j]] = TRUE;

      /* if the generator can be added, update the datastructures for graph components and colors */
      for (k = 0; k < npermvars; ++k)
      {
         int comp1;
         int comp2;
         int color1;
         int color2;
         int img;

         img = perm[k];
         assert( perm[img] == k );

         if ( img <= k )
            continue;

         comp1 = SCIPdisjointsetFind(vartocomponent, k);
         comp2 = SCIPdisjointsetFind(vartocomponent, img);

         /* components and colors don't have to be updated if the components are the same */
         if ( comp1 == comp2 )
            continue;

         color1 = SCIPdisjointsetFind(comptocolor, comp1);
         color2 = SCIPdisjointsetFind(comptocolor, comp2);

         if ( color1 != color2 )
         {
            SCIPdisjointsetUnion(comptocolor, firstcolor, color1, TRUE);
            SCIPdisjointsetUnion(comptocolor, firstcolor, color2, TRUE);
         }

         SCIPdisjointsetUnion(vartocomponent, comp1, comp2, FALSE);

         assert( SCIPdisjointsetFind(vartocomponent, k) == SCIPdisjointsetFind(vartocomponent, img) );
         assert( SCIPdisjointsetFind(comptocolor, SCIPdisjointsetFind(vartocomponent, k)) == firstcolor );
         assert( SCIPdisjointsetFind(comptocolor, SCIPdisjointsetFind(vartocomponent, img)) == firstcolor );
      }
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, graphcomponents, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(graphcompvartype.components), npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(graphcompvartype.colors), npermvars) );

   /*
    * At this point, we have built the colored graph. Now we transform the information in the
    * disjoint sets to the arrays graphcomponents, graphcompbegins, and compcolorbegins (see above).
    */

   /* build the struct graphcompvartype which is used to sort the graphcomponents array */
   for (j = 0; j < npermvars; ++j)
   {
      int comp;

      comp = SCIPdisjointsetFind(vartocomponent, j);

      graphcompvartype.components[j] = comp;
      graphcompvartype.colors[j] = SCIPdisjointsetFind(comptocolor, comp);

      (*graphcomponents)[j] = j;
   }

   /* sort graphcomponents first by color, then by component */
   SCIPsort(*graphcomponents, SYMsortGraphCompVars, (void*) &graphcompvartype, npermvars);

   *ngraphcomponents = SCIPdisjointsetGetComponentCount(vartocomponent);
   *ncompcolors = SCIPdisjointsetGetComponentCount(comptocolor);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, graphcompbegins, (*ngraphcomponents) + 1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, compcolorbegins, (*ncompcolors) + 1) );

   nextcolor = 1;
   nextcomp = 1;
   (*graphcompbegins)[0] = 0;
   (*compcolorbegins)[0] = 0;

   /* find the starting indices of new components and new colors */
   for (j = 1; j < npermvars; ++j)
   {
      int idx1;
      int idx2;

      idx1 = (*graphcomponents)[j];
      idx2 = (*graphcomponents)[j-1];

      assert( graphcompvartype.colors[idx1] >= graphcompvartype.colors[idx2] );

      if ( graphcompvartype.components[idx1] != graphcompvartype.components[idx2] )
      {
         (*graphcompbegins)[nextcomp] = j;

         if ( graphcompvartype.colors[idx1] > graphcompvartype.colors[idx2] )
         {
            (*compcolorbegins)[nextcolor] = nextcomp;
            ++nextcolor;
         }

         ++nextcomp;
      }
   }
   assert( nextcomp == *ngraphcomponents );
   assert( nextcolor == *ncompcolors );

   (*compcolorbegins)[nextcolor] = *ngraphcomponents;
   (*graphcompbegins)[nextcomp] = npermvars;

   SCIPfreeBufferArray(scip, &(graphcompvartype.colors));
   SCIPfreeBufferArray(scip, &(graphcompvartype.components));
   SCIPfreeBufferArray(scip, &componentslastperm);
   SCIPfreeDisjointset(scip, &comptocolor);
   SCIPfreeDisjointset(scip, &vartocomponent);

   return SCIP_OKAY;
}

/** adds an orbitope constraint for a suitable color of the subgroup graph */
static
SCIP_RETCODE addOrbitopeSubgroup(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  usedperms,          /**< array of the permutations that build the orbitope */
   int                   nusedperms,         /**< number of permutations in usedperms */
   int*                  compcolorbegins,    /**< array indicating where a new graphcolor begins */
   int*                  graphcompbegins,    /**< array indicating where a new graphcomponent begins */
   int*                  graphcomponents,    /**< array of all variable indices sorted by color and comp */
   int                   graphcoloridx,      /**< index of the graph color */
   int                   nrows,              /**< number of rows in the orbitope  */
   int                   ncols,              /**< number of columns in the orbitope  */
   int*                  firstvaridx,        /**< buffer to store the index of the largest variable (or NULL) */
   int*                  compidxfirstrow,    /**< buffer to store the comp index for the first row (or NULL) */
   int**                 lexorder,           /**< pointer to array storing lexicographic order defined by sub orbitopes */
   int*                  nvarslexorder,      /**< number of variables in lexicographic order */
   int*                  maxnvarslexorder,   /**< maximum number of variables in lexicographic order */
   SCIP_Bool             mayinteract,        /**< whether orbitope's symmetries might interact with other symmetries */
   SCIP_Bool*            success             /**< whether the orbitpe could be added */
   )
{  /*lint --e{571}*/
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR*** orbitopevarmatrix;
   SCIP_Shortbool* activevars;
   int** orbitopevaridx;
   int* columnorder;
   int* nusedelems;
   SCIP_CONS* cons;
   SCIP_Bool isorbitope;
   SCIP_Bool infeasible = FALSE;
   int nactivevars = 0;
   int ngencols = 0;
   int k;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( usedperms != NULL );
   assert( compcolorbegins != NULL );
   assert( graphcompbegins != NULL );
   assert( graphcomponents != NULL );
   assert( nusedperms > 0 );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( lexorder != NULL );
   assert( nvarslexorder != NULL );
   assert( maxnvarslexorder != NULL );

   *success = FALSE;

   /* create hashset to mark variables */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &activevars, propdata->npermvars) );

   /* orbitope matrix for indices of variables in permvars array */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx, nrows) );
   for (k = 0; k < nrows; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx[k], ncols) ); /*lint !e866*/
   }

   /* order of columns of orbitopevaridx */
   SCIP_CALL( SCIPallocBufferArray(scip, &columnorder, ncols) );
   for (k = 0; k < ncols; ++k)
      columnorder[k] = ncols + 1;

   /* count how often an element was used in the potential orbitope */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nusedelems, propdata->npermvars) );

   /* mark variables in this subgroup orbitope */
   for (k = compcolorbegins[graphcoloridx]; k < compcolorbegins[graphcoloridx+1]; ++k)
   {
      SCIP_VAR* firstvar;
      int compstart;
      int l;

      compstart = graphcompbegins[k];
      firstvar = propdata->permvars[graphcomponents[compstart]];

      if ( ! SCIPvarIsBinary(firstvar) )
         continue;

      for (l = 0; l < ncols; ++l)
      {
         int varidx;

         varidx = graphcomponents[compstart + l];
         assert( ! activevars[varidx] );

         activevars[varidx] = TRUE;
         ++nactivevars;
      }
   }
   assert( nactivevars == nrows * ncols );

   /* build the variable index matrix for the orbitope
    *
    * It is possible that we find an orbitope, but not using all possible columns. For example
    * (1,2), (2,3), (3,4), (3,5) defines the symmetric group on {1,2,3,4,5}, but the generators
    * we expect in our construction need shape (1,2), (2,3), (3,4), (4,5). For this reason,
    * we need to store how many columns have been generated.
    *
    * @todo ensure compatibility with more general generators
    */
   SCIP_CALL( checkTwoCyclePermsAreOrbitope(scip, propdata->permvars, propdata->npermvars,
         propdata->perms, usedperms, nrows, nusedperms, orbitopevaridx, columnorder,
         nusedelems, &ngencols, NULL, &isorbitope, activevars) );

   /* it might happen that we cannot detect the orbitope if it is generated by permutations with different
    *  number of 2-cycles.
    */
   if ( ! isorbitope )
   {
      SCIPfreeBufferArray(scip, &nusedelems);
      SCIPfreeBufferArray(scip, &columnorder);
      SCIPfreeBufferArray(scip, &orbitopevaridx);
      SCIPfreeBufferArray(scip, &activevars);

      return SCIP_OKAY;
   }

   /* There are three possibilities for the structure of columnorder:
    * 1)  [0, 1, -1, -1, ..., -1]
    * 2)  [0, 1, 1, 1, ..., 1]
    * 3)  [0, 1, -1, -1, ...., -1, 1, 1, ..., 1]
    *
    * The '1'-columns will be added to the matrix first and in the last 2
    * cases the method starts from the right. So to store the variable index
    * that will be in the upper-left corner, we need either the entryin the
    * second column (case 1) or the entry in the last column (cases 2 and 3).
    */
   if ( firstvaridx != NULL )
   {
      if ( columnorder[ngencols-1] > -1 )
         *firstvaridx = orbitopevaridx[0][ngencols-1];
      else
         *firstvaridx = orbitopevaridx[0][1];
   }

   /* find corresponding graphcomponent of first variable (needed for weak sbcs) */
   if ( compidxfirstrow != NULL && firstvaridx != NULL )
   {
      *compidxfirstrow = -1;

      for (k = compcolorbegins[graphcoloridx]; k < compcolorbegins[graphcoloridx+1] && (*compidxfirstrow) < 0; ++k)
      {
         SCIP_VAR* firstvar;
         int compstart;
         int l;

         compstart = graphcompbegins[k];
         firstvar = propdata->permvars[graphcomponents[compstart]];

         if ( ! SCIPvarIsBinary(firstvar) )
            continue;

         /* iterate over all columns (elements in orbit), because we cannot see from ngencols which columns
          * have been left out
          */
         for (l = 0; l < ncols; ++l)
         {
            if ( graphcomponents[compstart + l] == *firstvaridx )
            {
               *compidxfirstrow = k;
               break;
            }
         }
      }
      assert( *compidxfirstrow > -1 );
   }

   /* prepare orbitope variable matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix, nrows) );
   for (k = 0; k < nrows; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix[k], ngencols) );
   }

   /* build the matrix containing the actual variables of the orbitope */
   SCIP_CALL( SCIPgenerateOrbitopeVarsMatrix(scip, &orbitopevarmatrix, nrows, ngencols,
         propdata->permvars, propdata->npermvars, orbitopevaridx, columnorder,
         nusedelems, NULL, &infeasible, TRUE, lexorder, nvarslexorder, maxnvarslexorder) );

   assert( ! infeasible );
   assert( firstvaridx == NULL || propdata->permvars[*firstvaridx] == orbitopevarmatrix[0][0] );

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "suborbitope_%d_%d", graphcoloridx, propdata->norbitopes);

   SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, orbitopevarmatrix,
         SCIP_ORBITOPETYPE_FULL, nrows, ngencols, FALSE, mayinteract, FALSE, FALSE, propdata->conssaddlp,
         TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   *success = TRUE;

   /* do not release constraint here - will be done later */
   propdata->genorbconss[propdata->ngenorbconss++] = cons;
   ++propdata->norbitopes;

   for (k = nrows - 1; k >= 0; --k)
      SCIPfreeBufferArray(scip, &orbitopevarmatrix[k]);
   SCIPfreeBufferArray(scip, &orbitopevarmatrix);
   SCIPfreeBufferArray(scip, &nusedelems);
   SCIPfreeBufferArray(scip, &columnorder);
   for (k = nrows - 1; k >= 0; --k)
      SCIPfreeBufferArray(scip, &orbitopevaridx[k]);
   SCIPfreeBufferArray(scip, &orbitopevaridx);
   SCIPfreeBufferArray(scip, &activevars);

   return SCIP_OKAY;
}

/** adds strong SBCs for a suitable color of the subgroup graph */
static
SCIP_RETCODE addStrongSBCsSubgroup(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  graphcompbegins,    /**< array indicating where a new graphcomponent begins */
   int*                  graphcomponents,    /**< array of all variable indices sorted by color and comp */
   int                   graphcompidx,       /**< index of the graph component */
   SCIP_Bool             storelexorder,      /**< whether the lexicographic order induced by the orbitope shall be stored */
   int**                 lexorder,           /**< pointer to array storing lexicographic order defined by sub orbitopes */
   int*                  nvarsorder,         /**< number of variables in lexicographic order */
   int*                  maxnvarsorder       /**< maximum number of variables in lexicographic order */
   )
{
   int k;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( graphcompbegins != NULL );
   assert( graphcomponents != NULL );
   assert( graphcompidx >= 0 );
   assert( ! storelexorder || lexorder != NULL );
   assert( ! storelexorder || nvarsorder != NULL );
   assert( ! storelexorder || maxnvarsorder != NULL );

   /* possibly store lexicographic order defined by strong SBCs */
   if ( storelexorder )
   {
      if ( *maxnvarsorder == 0 )
      {
         *maxnvarsorder = graphcompbegins[graphcompidx + 1] - graphcompbegins[graphcompidx + 1];
         *nvarsorder = 0;

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, lexorder, *maxnvarsorder) );
      }
      else
      {
         assert( *nvarsorder == *maxnvarsorder );

         *maxnvarsorder += graphcompbegins[graphcompidx + 1] - graphcompbegins[graphcompidx + 1];

         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, lexorder, *nvarsorder, *maxnvarsorder) );
      }

      (*lexorder)[*nvarsorder++] = graphcomponents[graphcompbegins[graphcompidx]];
   }

   /* add strong SBCs (lex-max order) for chosen graph component */
   for (k = graphcompbegins[graphcompidx]+1; k < graphcompbegins[graphcompidx+1]; ++k)
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_CONS* cons;
      SCIP_VAR* vars[2];
      SCIP_Real vals[2] = {1, -1};

      vars[0] = propdata->permvars[graphcomponents[k-1]];
      vars[1] = propdata->permvars[graphcomponents[k]];

      if ( storelexorder )
         (*lexorder)[*nvarsorder++] = graphcomponents[k];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "strong_sbcs_%s_%s", SCIPvarGetName(vars[0]), SCIPvarGetName(vars[1]));

      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, 0.0,
            SCIPinfinity(scip), propdata->conssaddlp, propdata->conssaddlp, TRUE,
            FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* check whether we need to resize */
      if ( propdata->ngenlinconss >= propdata->genlinconsssize )
      {
         int newsize;

         newsize = SCIPcalcMemGrowSize(scip, propdata->ngenlinconss + 1);
         assert( newsize > propdata->ngenlinconss );

         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &propdata->genlinconss, propdata->genlinconsssize, newsize) );

         propdata->genlinconsssize = newsize;
      }

      propdata->genlinconss[propdata->ngenlinconss] = cons;
      ++propdata->ngenlinconss;
   }

   return SCIP_OKAY;
}

/** adds weak SBCs for a suitable color of the subgroup graph */
static
SCIP_RETCODE addWeakSBCsSubgroup(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  compcolorbegins,    /**< array indicating where a new graphcolor begins */
   int*                  graphcompbegins,    /**< array indicating where a new graphcomponent begins */
   int*                  graphcomponents,    /**< array of all variable indices sorted by color and comp */
   int                   ncompcolors,        /**< number of colors in the graph  */
   int*                  chosencomppercolor, /**< array indicating which comp was handled per color */
   int*                  firstvaridxpercolor,/**< array indicating the largest variable per color */
   int                   symgrpcompidx,      /**< index of the component of the symmetry group */
   int*                  naddedconss,        /**< buffer to store the number of added constraints */
   SCIP_Bool             storelexorder,      /**< whether the lexicographic order induced by the orbitope shall be stored */
   int**                 lexorder,           /**< pointer to array storing lexicographic order defined by sub orbitopes */
   int*                  nvarsorder,         /**< number of variables in lexicographic order */
   int*                  maxnvarsorder       /**< maximum number of variables in lexicographic order */
   )
{  /*lint --e{571}*/
   SCIP_HASHMAP* varsinlexorder;
   SCIP_Shortbool* usedvars;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2] = {1, -1};
   SCIP_Shortbool* varfound;
   int* orbit[2];
   int orbitsize[2] = {1, 1};
   int activeorb = 0;
   int chosencolor = -1;
   int j;
   int k;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( compcolorbegins != NULL );
   assert( graphcompbegins != NULL );
   assert( graphcomponents != NULL );
   assert( firstvaridxpercolor != NULL );
   assert( chosencomppercolor != NULL );
   assert( naddedconss != NULL );
   assert( symgrpcompidx >= 0 );
   assert( symgrpcompidx < propdata->ncomponents );
   assert( ! storelexorder || lexorder != NULL );
   assert( ! storelexorder || nvarsorder != NULL );
   assert( ! storelexorder || maxnvarsorder != NULL );

   *naddedconss = 0;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &usedvars, propdata->npermvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &varfound, propdata->npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbit[0], propdata->npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbit[1], propdata->npermvars) );

   /* Store the entries in lexorder in a hashmap, for fast lookups. */
   if ( lexorder == NULL || *lexorder == NULL )
   {
      /* Lexorder does not exist, so do not create hashmap. */
      varsinlexorder = NULL;
   }
   else
   {
      assert( *maxnvarsorder >= 0 );
      assert( *nvarsorder >= 0 );

      SCIP_CALL( SCIPhashmapCreate(&varsinlexorder, SCIPblkmem(scip), *maxnvarsorder) );

      for (k = 0; k < *nvarsorder; ++k)
      {
         /* add element from lexorder to hashmap.
          * Use insert, as duplicate entries in lexorder is not permitted. */
         assert( ! SCIPhashmapExists(varsinlexorder, (void*) (long) (*lexorder)[k]) ); /* Use int as pointer */
         SCIP_CALL( SCIPhashmapInsertInt(varsinlexorder, (void*) (long) (*lexorder)[k], k) );
      }
   }

   /* We will store the newest and the largest orbit and activeorb will be used to mark at which entry of the array
    * orbit the newly computed one will be stored. */
   for (j = 0; j < ncompcolors; ++j)
   {
      int graphcomp;
      int graphcompsize;
      int varidx;

      /* skip color for which we did not add anything */
      if ( chosencomppercolor[j] < 0 )
         continue;

      assert( firstvaridxpercolor[j] >= 0 );

      graphcomp = chosencomppercolor[j];
      graphcompsize = graphcompbegins[graphcomp+1] - graphcompbegins[graphcomp];
      varidx = firstvaridxpercolor[j];

      /* if the first variable was already contained in another orbit or if there are no variables left anyway, skip the
       * component */
      if ( varfound[varidx] || graphcompsize == propdata->npermvars )
         continue;

      /* If varidx is in lexorder, then it must be the first entry of lexorder. */
      if ( varsinlexorder != NULL
         && SCIPhashmapExists(varsinlexorder, (void*) (long) varidx)
         && lexorder != NULL && *lexorder != NULL && *maxnvarsorder > 0 && *nvarsorder > 0
         && (*lexorder)[0] != varidx )
         continue;

      /* mark all variables that have been used in strong SBCs */
      for (k = graphcompbegins[graphcomp]; k < graphcompbegins[graphcomp+1]; ++k)
      {
         assert( 0 <= graphcomponents[k] && graphcomponents[k] < propdata->npermvars );

         usedvars[graphcomponents[k]] = TRUE;
      }

      SCIP_CALL( SCIPcomputeOrbitVar(scip, propdata->npermvars, propdata->perms,
            propdata->permstrans, propdata->components, propdata->componentbegins,
            usedvars, varfound, varidx, symgrpcompidx,
            orbit[activeorb], &orbitsize[activeorb]) );

      assert( orbit[activeorb][0] ==  varidx );

      if ( orbitsize[activeorb] > orbitsize[1 - activeorb] ) /*lint !e514*/
      {
         /* if the new orbit is larger then the old largest one, flip activeorb */
         activeorb = 1 - activeorb;
         chosencolor = j;
      }

      /* reset array */
      for (k = graphcompbegins[graphcomp]; k < graphcompbegins[graphcomp+1]; ++k)
         usedvars[graphcomponents[k]] = FALSE;
   }

   /* check if we have found at least one non-empty orbit */
   if ( chosencolor > -1 )
   {
      /* flip activeorb again to avoid confusion, it is then at the largest orbit */
      activeorb = 1 - activeorb;

      assert( orbit[activeorb][0] == firstvaridxpercolor[chosencolor] );
      vars[0] = propdata->permvars[orbit[activeorb][0]];

      assert( chosencolor > -1 );
      SCIPdebugMsg(scip, "    adding %d weak sbcs for enclosing orbit of color %d.\n", orbitsize[activeorb]-1, chosencolor);

      *naddedconss = orbitsize[activeorb] - 1;

      /* add weak SBCs for rest of enclosing orbit */
      for (j = 1; j < orbitsize[activeorb]; ++j)
      {
         SCIP_CONS* cons;
         char name[SCIP_MAXSTRLEN];

         vars[1] = propdata->permvars[orbit[activeorb][j]];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "weak_sbcs_%d_%s_%s", symgrpcompidx, SCIPvarGetName(vars[0]), SCIPvarGetName(vars[1]));

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, 0.0,
               SCIPinfinity(scip), propdata->conssaddlp, propdata->conssaddlp, TRUE,
               FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         /* check whether we need to resize */
         if ( propdata->ngenlinconss >= propdata->genlinconsssize )
         {
            int newsize;

            newsize = SCIPcalcMemGrowSize(scip, propdata->ngenlinconss + 1);
            assert( newsize > propdata->ngenlinconss );

            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &propdata->genlinconss, propdata->genlinconsssize, newsize) );

            propdata->genlinconsssize = newsize;
         }

         propdata->genlinconss[propdata->ngenlinconss] = cons;
         ++propdata->ngenlinconss;
      }

      /* possibly store lexicographic order defined by weak SBCs */
      if ( storelexorder )
      {
         int varidx;

         varidx = orbit[activeorb][0];

         if ( *maxnvarsorder == 0 )
         {
            *maxnvarsorder = 1;
            *nvarsorder = 0;

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, lexorder, *maxnvarsorder) );
            (*lexorder)[(*nvarsorder)++] = varidx;
         }
         else
         {
            assert( *nvarsorder == *maxnvarsorder );
            assert( varsinlexorder != NULL );
            assert( lexorder != NULL );
            assert( *lexorder != NULL );

            /* the leader of the weak inequalities has to be the first element in the lexicographic order */
            if ( varidx == (*lexorder)[0] )
            {
               /* lexorder is already ok!! */
               assert( SCIPhashmapExists(varsinlexorder, (void*) (long) varidx) );
            }
            else
            {
               /* Then varidx must not be in the lexorder,
                * We must add it at the front of the array, and maintain the current order. */
               assert( ! SCIPhashmapExists(varsinlexorder, (void*) (long) varidx) );

               ++(*maxnvarsorder);
               ++(*nvarsorder);

               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, lexorder, *nvarsorder, *maxnvarsorder) );

               /* Shift array by one position to the right */
               for (k = *maxnvarsorder - 1; k >= 1; --k)
                  (*lexorder)[k] = (*lexorder)[k - 1];

               (*lexorder)[0] = varidx;
            }
         }
      }
   }
   else
      SCIPdebugMsg(scip, "  no further weak sbcs are valid\n");

   SCIPfreeBufferArray(scip, &orbit[1]);
   SCIPfreeBufferArray(scip, &orbit[0]);
   if ( varsinlexorder != NULL )
      SCIPhashmapFree(&varsinlexorder);
   SCIPfreeBufferArray(scip, &varfound);
   SCIPfreeCleanBufferArray(scip, &usedvars);

   return SCIP_OKAY;
}


/** temporarily adapt symmetry data to new variable order given by Schreier Sims */
static
SCIP_RETCODE adaptSymmetryDataSST(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 origperms,          /**< permutation matrix w.r.t. original variable ordering */
   int**                 modifiedperms,      /**< memory for permutation matrix w.r.t. new variable ordering */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            origpermvars,       /**< array of permutation vars w.r.t. original variable ordering */
   SCIP_VAR**            modifiedpermvars,   /**< memory for array of permutation vars w.r.t. new variable ordering */
   int                   npermvars,          /**< length or modifiedpermvars array */
   int*                  leaders,            /**< leaders of Schreier Sims constraints */
   int                   nleaders            /**< number of leaders */
   )
{
   int* permvaridx;
   int* posinpermvar;
   int leader;
   int curposleader;
   int varidx;
   int lidx;
   int i;
   int l;
   int p;

   assert( scip != NULL );
   assert( origperms != NULL );
   assert( modifiedperms != NULL );
   assert( nperms > 0 );
   assert( origpermvars != NULL );
   assert( modifiedpermvars != NULL );
   assert( npermvars > 0 );
   assert( leaders != NULL );
   assert( nleaders > 0 );

   /* initialize map from position in lexicographic order to index of original permvar */
   SCIP_CALL( SCIPallocBufferArray(scip, &permvaridx, npermvars) );
   for (i = 0; i < npermvars; ++i)
      permvaridx[i] = i;

   /* initialize map from permvaridx to its current position in the reordered permvars array */
   SCIP_CALL( SCIPallocBufferArray(scip, &posinpermvar, npermvars) );
   for (i = 0; i < npermvars; ++i)
      posinpermvar[i] = i;

   /* Iterate over leaders and put the l-th leader to the l-th position of the lexicographic order.
    * We do this by swapping the l-th leader with the element at position l of the current permvars array. */
   for (l = 0; l < nleaders; ++l)
   {
      leader = leaders[l];
      curposleader = posinpermvar[leader];
      varidx = permvaridx[curposleader];
      lidx = permvaridx[l];

      /* swap the permvar at position l with the l-th leader */
      permvaridx[curposleader] = lidx;
      permvaridx[l] = varidx;

      /* update the position map */
      posinpermvar[lidx] = curposleader;
      posinpermvar[leader] = l;
   }

   /* update the permvars array to new variable order */
   for (i = 0; i < npermvars; ++i)
      modifiedpermvars[i] = origpermvars[permvaridx[i]];

   /* update the permutation to the new variable order */
   for (p = 0; p < nperms; ++p)
   {
      for (i = 0; i < npermvars; ++i)
         modifiedperms[p][i] = posinpermvar[origperms[p][permvaridx[i]]];
   }

   SCIPfreeBufferArray(scip, &permvaridx);
   SCIPfreeBufferArray(scip, &posinpermvar);

   return SCIP_OKAY;
}


/* returns the number of found orbitopes with at least three columns per graph component or 0
 * if the found orbitopes do not satisfy certain criteria for being used
 */
static
int getNOrbitopesInComp(
   SCIP_VAR**            permvars,           /**< array of variables affected by symmetry */
   int*                  graphcomponents,    /**< array of graph components */
   int*                  graphcompbegins,    /**< array indicating starting position of graph components */
   int*                  compcolorbegins,    /**< array indicating starting positions of potential orbitopes */
   int                   ncompcolors,        /**< number of components encoded in compcolorbegins */
   int                   symcompsize         /**< size of symmetry component for that we detect suborbitopes */
   )
{
   SCIP_Bool oneorbitopecriterion = FALSE;
   SCIP_Bool multorbitopecriterion = FALSE;
   int norbitopes = 0;
   int j;

   assert( graphcompbegins != NULL );
   assert( compcolorbegins != NULL );
   assert( ncompcolors >= 0 );
   assert( symcompsize > 0 );

   for (j = 0; j < ncompcolors; ++j)
   {
      SCIP_VAR* firstvar;
      int largestcompsize = 0;
      int nbinrows= 0;
      int k;

      /* skip trivial components */
      if ( graphcompbegins[compcolorbegins[j+1]] - graphcompbegins[compcolorbegins[j]] < 2 )
         continue;

      /* check whether components of this color build an orbitope (with > 2 columns) */
      for (k = compcolorbegins[j]; k < compcolorbegins[j+1]; ++k)
      {
         int compsize;

         compsize = graphcompbegins[k+1] - graphcompbegins[k];

         /* the first component that we are looking at for this color */
         if ( largestcompsize < 1 )
         {
            if ( compsize < 3 )
               break;

            largestcompsize = compsize;
         }
         else if ( compsize != largestcompsize )
            break;

         firstvar = permvars[graphcomponents[graphcompbegins[k]]];

         /* count number of binary orbits (comps) */
         if ( SCIPvarIsBinary(firstvar) )
            ++nbinrows;
      }

      /* we have found an orbitope */
      if ( k == compcolorbegins[j+1] )
      {
         SCIP_Real threshold;
         int ncols;

         ++norbitopes;
         ncols = graphcompbegins[compcolorbegins[j] + 1] - graphcompbegins[compcolorbegins[j]];

         threshold = 0.7 * (SCIP_Real) symcompsize;

         /* check whether criteria for adding orbitopes are satisfied */
         if ( nbinrows <= 2 * ncols || (nbinrows <= 8 * ncols && nbinrows < 100) )
            multorbitopecriterion = TRUE;
         else if ( nbinrows <= 3 * ncols || (SCIP_Real) nbinrows * ncols >= threshold )
            oneorbitopecriterion = TRUE;
      }
   }

   if ( (norbitopes == 1 && oneorbitopecriterion) || (norbitopes >= 2 && multorbitopecriterion) )
      return norbitopes;

   return 0;
}


/** checks whether subgroups of the components are symmetric groups and adds SBCs for them */
static
SCIP_RETCODE detectAndHandleSubgroups(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata            /**< pointer to data of symmetry propagator */
   )
{
   int* genorder;
   int i;
#ifdef SCIP_DEBUG
   int norbitopes = 0;
   int nstrongsbcs = 0;
   int nweaksbcs = 0;
#endif
   int** modifiedperms;
   SCIP_VAR** modifiedpermvars;
   int* nvarsincomponent;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->computedsymmetry );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );
   assert( propdata->ncomponents > 0 );
   assert( propdata->nperms >= 0 );

   /* exit if no symmetry is present */
   if ( propdata->nperms == 0 )
      return SCIP_OKAY;

   /* exit if instance is too large */
   if ( SCIPgetNConss(scip) > propdata->maxnconsssubgroup )
      return SCIP_OKAY;

   assert( propdata->nperms > 0 );
   assert( propdata->perms != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->permvars != NULL );

   /* create array for permutation order */
   SCIP_CALL( SCIPallocBufferArray(scip, &genorder, propdata->nperms) );

   /* create arrays for modified permutations in case we adapt the lexicographic order because of suborbitopes */
   SCIP_CALL( SCIPallocBufferArray(scip, &modifiedperms, propdata->nperms) );
   for (i = 0; i < propdata->nperms; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &modifiedperms[i], propdata->npermvars) );
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &modifiedpermvars, propdata->npermvars) );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarsincomponent, propdata->npermvars) );
   for (i = 0; i < propdata->npermvars; ++i)
   {
      if ( propdata->vartocomponent[i] >= 0 )
         ++nvarsincomponent[propdata->vartocomponent[i]];
   }

   SCIPdebugMsg(scip, "starting subgroup detection routine for %d components\n", propdata->ncomponents);

   /* iterate over components */
   for (i = 0; i < propdata->ncomponents; ++i)
   {
      int* graphcomponents;
      int* graphcompbegins;
      int* compcolorbegins;
      int* chosencomppercolor = NULL;
      int* firstvaridxpercolor = NULL;
      int* usedperms;
      int usedpermssize;
      int ngraphcomponents;
      int ncompcolors;
      int ntwocycleperms;
      int npermsincomp;
      int nusedperms;
      int ntrivialcolors = 0;
      int j;
      int* lexorder = NULL;
      int nvarslexorder = 0;
      int maxnvarslexorder = 0;
      SCIP_Shortbool* permused;
      SCIP_Bool allpermsused = FALSE;
      SCIP_Bool handlednonbinarysymmetry = FALSE;
      int norbitopesincomp;

      /* if component is blocked, skip it */
      if ( propdata->componentblocked[i] )
      {
         SCIPdebugMsg(scip, "component %d has already been handled and will be skipped\n", i);
         continue;
      }

      npermsincomp = propdata->componentbegins[i + 1] - propdata->componentbegins[i];

      /* set the first npermsincomp entries of genorder; the others are not used for this component */
      for (j = 0; j < npermsincomp; ++j)
         genorder[j] = j;

      SCIP_CALL( chooseOrderOfGenerators(scip, propdata, i, &genorder, &ntwocycleperms) );

      assert( ntwocycleperms >= 0 );
      assert( ntwocycleperms <= npermsincomp );

      SCIPdebugMsg(scip, "component %d has %d permutations consisting of 2-cycles\n", i, ntwocycleperms);

#ifdef SCIP_MORE_DEBUG
      SCIP_Bool* used;
      int perm;
      int p;
      int k;

      SCIP_CALL( SCIPallocBufferArray(scip, &used, propdata->npermvars) );
      for (p = propdata->componentbegins[i]; p < propdata->componentbegins[i+1]; ++p)
      {
         perm = propdata->components[p];

         for (k = 0; k < propdata->npermvars; ++k)
            used[k] = FALSE;

         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "permutation %d\n", perm);

         for (k = 0; k < propdata->npermvars; ++k)
         {
            if ( used[k] )
               continue;

            j = propdata->perms[perm][k];

            if ( k == j )
               continue;

            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "(%s,", SCIPvarGetName(propdata->permvars[k]));
            used[k] = TRUE;
            while (j != k)
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s,", SCIPvarGetName(propdata->permvars[j]));
               used[j] = TRUE;

               j = propdata->perms[perm][j];
            }
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ")");
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");
      }

      SCIPfreeBufferArray(scip, &used);
#endif

      if ( ntwocycleperms < 2 )
      {
         SCIPdebugMsg(scip, "  --> skip\n");
         continue;
      }

      usedpermssize = ntwocycleperms / 2;
      SCIP_CALL( SCIPallocBufferArray(scip, &usedperms, usedpermssize) );
      SCIP_CALL( SCIPallocClearBufferArray(scip, &permused, npermsincomp) );

      SCIP_CALL( buildSubgroupGraph(scip, propdata, genorder, ntwocycleperms, i,
            &graphcomponents, &graphcompbegins, &compcolorbegins, &ngraphcomponents,
            &ncompcolors, &usedperms, &nusedperms, usedpermssize, permused) );

      SCIPdebugMsg(scip, "  created subgroup detection graph using %d of the permutations\n", nusedperms);

      if ( nusedperms == npermsincomp )
         allpermsused = TRUE;

      assert( graphcomponents != NULL );
      assert( graphcompbegins != NULL );
      assert( compcolorbegins != NULL );
      assert( ngraphcomponents > 0 );
      assert( ncompcolors > 0 );
      assert( nusedperms <= ntwocycleperms );
      assert( ncompcolors < propdata->npermvars );

      if ( nusedperms == 0 )
      {
         SCIPdebugMsg(scip, "  -> skipping component, since less no permutation was used\n");

         SCIPfreeBufferArray(scip, &permused);
         SCIPfreeBufferArray(scip, &usedperms);

         continue;
      }

      SCIPdebugMsg(scip, "  number of different colors in the graph: %d\n", ncompcolors);

      if ( propdata->addstrongsbcs || propdata->addweaksbcs )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &chosencomppercolor, ncompcolors) );
         SCIP_CALL( SCIPallocBufferArray(scip, &firstvaridxpercolor, ncompcolors) );

         /* Initialize the arrays with -1 to encode that we have not added orbitopes/strong SBCs
          * yet. In case we do not modify this entry, no weak inequalities are added based on
          * this component.
          */
         for (j = 0; j < ncompcolors; ++j)
         {
            chosencomppercolor[j] = -1;
            firstvaridxpercolor[j] = -1;
         }
      }

      norbitopesincomp = getNOrbitopesInComp(propdata->permvars, graphcomponents, graphcompbegins, compcolorbegins,
         ncompcolors, nvarsincomponent[i]);

      /* if there is just one orbitope satisfying the requirements, handle the full component by symresacks */
      if ( norbitopesincomp == 1 )
      {
         int k;

         for (k = 0; k < npermsincomp; ++k)
         {
            SCIP_CONS* cons;
            char name[SCIP_MAXSTRLEN];

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_comp%d_perm%d", i, k);

            SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, propdata->perms[propdata->components[propdata->componentbegins[i] + k]],
                  propdata->permvars, propdata->npermvars, FALSE,
                  propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons));

            /* do not release constraint here - will be done later */
            propdata->genorbconss[propdata->ngenorbconss++] = cons;
            ++propdata->nsymresacks;

            if ( ! propdata->componentblocked[i] )
            {
               propdata->componentblocked[i] |= SYM_HANDLETYPE_SYMBREAK;
               ++propdata->ncompblocked;
            }

            SCIPdebugMsg(scip, "  add symresack for permutation %d of component %d\n", k, i);
         }

         goto FREELOOPMEMORY;
      }

      for (j = 0; j < ncompcolors; ++j)
      {
         int nbinarycomps = 0;
         int largestcolorcomp = -1;
         int largestcompsize = 0;
         int k;
         SCIP_Bool isorbitope = TRUE;
         SCIP_Bool orbitopeadded = FALSE;
         SCIP_Bool useorbitope;
#ifdef SCIP_DEBUG
         SCIP_Bool binaffected = FALSE;
         SCIP_Bool intaffected = FALSE;
         SCIP_Bool contaffected = FALSE;
#endif

         /* skip trivial components */
         if ( graphcompbegins[compcolorbegins[j+1]] - graphcompbegins[compcolorbegins[j]] < 2 )
         {
            if( chosencomppercolor != NULL )
               chosencomppercolor[j] = -1;

            ++ntrivialcolors;
            continue;
         }

         SCIPdebugMsg(scip, "    color %d has %d components with overall %d variables\n", j, compcolorbegins[j+1] - compcolorbegins[j],
            graphcompbegins[compcolorbegins[j+1]] - graphcompbegins[compcolorbegins[j]]);

         /* check whether components of this color might build an orbitope (with > 2 columns) */
         for (k = compcolorbegins[j]; k < compcolorbegins[j+1]; ++k)
         {
            SCIP_VAR* firstvar;
            int compsize;

            compsize = graphcompbegins[k+1] - graphcompbegins[k];

            /* the first component that we are looking at for this color */
            if ( largestcompsize < 1 )
            {
               if ( compsize < 3 )
               {
                  isorbitope = FALSE;
                  break;
               }

               largestcompsize = compsize;
               largestcolorcomp = k;
            }
            else if ( compsize != largestcompsize )
            {
               /* variable orbits (compsize) have not the same size, cannot define orbitope */
               isorbitope = FALSE;
               break;
            }

            firstvar = propdata->permvars[graphcomponents[graphcompbegins[k]]];

            /* count number of binary orbits (comps) */
            if ( SCIPvarIsBinary(firstvar) )
               ++nbinarycomps;
         }

#ifdef SCIP_DEBUG
         for (k = compcolorbegins[j]; k < compcolorbegins[j+1]; ++k)
         {
            SCIP_VAR* firstvar;

            firstvar = propdata->permvars[graphcomponents[graphcompbegins[k]]];

            if ( SCIPvarIsBinary(firstvar) )
               binaffected = TRUE;
            else if (SCIPvarIsIntegral(firstvar) )
               intaffected = TRUE;
            else
               contaffected = TRUE;
         }

         SCIPdebugMsg(scip, "      affected types (bin,int,cont): (%d,%d,%d)\n", binaffected, intaffected, contaffected);
#endif

         /* only use the orbitope if there are binary rows */
         useorbitope = FALSE;
         if ( norbitopesincomp > 0 && nbinarycomps > 0 )
            useorbitope = TRUE;

         if ( isorbitope && useorbitope )
         {
            int firstvaridx;
            int chosencomp;

            SCIPdebugMsg(scip, "      detected an orbitope with %d rows and %d columns\n", nbinarycomps, largestcompsize);

            assert( nbinarycomps > 0 );
            assert( largestcompsize > 2 );

            /* add the orbitope constraint for this color
             *
             * It might happen that we cannot generate the orbitope matrix if the orbitope is not generated by permutations
             * all having the same number of 2-cycles, e.g., the orbitope generated by (1,2)(4,5), (2,3), (5,6).
             */
            SCIP_CALL( addOrbitopeSubgroup(scip, propdata, usedperms, nusedperms, compcolorbegins,
                  graphcompbegins, graphcomponents, j, nbinarycomps, largestcompsize, &firstvaridx, &chosencomp,
                  &lexorder, &nvarslexorder, &maxnvarslexorder, allpermsused, &orbitopeadded) );

            if ( orbitopeadded )
            {
               if ( propdata->addstrongsbcs || propdata->addweaksbcs )
               {
                  assert( chosencomppercolor != NULL );
                  assert( firstvaridxpercolor != NULL );

                  /* adapt the first variable per color to be compatible with the created orbiope (upper left variable) */
                  assert( compcolorbegins[j] <= chosencomp && chosencomp < compcolorbegins[j+1] );
                  assert( 0 <= firstvaridx && firstvaridx < propdata->npermvars );

                  chosencomppercolor[j] = chosencomp;
                  firstvaridxpercolor[j] = firstvaridx;
               }

               if ( ! propdata->componentblocked[i] )
               {
                  propdata->componentblocked[i] |= SYM_HANDLETYPE_SYMBREAK;
                  ++propdata->ncompblocked;
               }

#ifdef SCIP_DEBUG
               ++norbitopes;
#endif
            }
         }

         /* if no (useable) orbitope was found, possibly add strong SBCs */
         if ( propdata->addstrongsbcs && ! orbitopeadded )
         {
            assert( largestcolorcomp >= 0 );
            assert( largestcolorcomp < ngraphcomponents );
            assert( largestcompsize > 0 );

            if( propdata->addweaksbcs )
            {
               assert( chosencomppercolor != NULL );
               assert( firstvaridxpercolor != NULL );

               chosencomppercolor[j] = largestcolorcomp;
               firstvaridxpercolor[j] = graphcomponents[graphcompbegins[largestcolorcomp]];
            }

            SCIPdebugMsg(scip, "      choosing component %d with %d variables and adding strong SBCs\n",
               largestcolorcomp, graphcompbegins[largestcolorcomp+1] - graphcompbegins[largestcolorcomp]);

            /* add the strong SBCs for the corresponding component */
            SCIP_CALL( addStrongSBCsSubgroup(scip, propdata, graphcompbegins, graphcomponents, largestcolorcomp,
                  propdata->addsymresacks, &lexorder, &nvarslexorder, &maxnvarslexorder) );

            /* store whether symmetries on non-binary symmetries have been handled */
            if ( ! SCIPvarIsBinary(propdata->permvars[graphcomponents[graphcompbegins[largestcolorcomp]]]) )
               handlednonbinarysymmetry = TRUE;

            if ( ! propdata->componentblocked[i] )
            {
               propdata->componentblocked[i] |= SYM_HANDLETYPE_SYMBREAK;
               ++propdata->ncompblocked;
            }

#ifdef SCIP_DEBUG
            nstrongsbcs += graphcompbegins[largestcolorcomp+1] - graphcompbegins[largestcolorcomp] - 1;
#endif
         }
         else if ( ! orbitopeadded )
         {
            SCIPdebugMsg(scip, "      no useable orbitope found and no SBCs added\n");

            /* mark the color as not handled */
            if ( propdata->addweaksbcs )
            {
               assert( chosencomppercolor != NULL );
               chosencomppercolor[j] = -1; /*lint !e613*/
            }
         }
      }

      SCIPdebugMsg(scip, "    skipped %d trivial colors\n", ntrivialcolors);

      /* possibly add weak SBCs for enclosing orbit of first component */
      if ( propdata->addweaksbcs && propdata->componentblocked[i] && nusedperms < npermsincomp )
      {
         int naddedconss;

         assert( firstvaridxpercolor != NULL );
         assert( chosencomppercolor != NULL );

         SCIP_CALL( addWeakSBCsSubgroup(scip, propdata, compcolorbegins, graphcompbegins,
               graphcomponents, ncompcolors, chosencomppercolor, firstvaridxpercolor,
               i, &naddedconss, propdata->addsymresacks, &lexorder, &nvarslexorder, &maxnvarslexorder) );

         assert( naddedconss < propdata->npermvars );

#ifdef SCIP_DEBUG
         nweaksbcs += naddedconss;
#endif
      }
      else
         SCIPdebugMsg(scip, "  don't add weak sbcs because all generators were used or the settings forbid it\n");

      /* if suborbitopes or strong group actions have been found, potentially add symresacks adapted to
       * variable order given by lexorder if no symmetries on non-binary variables have been handled
       */
      if ( nvarslexorder > 0 && propdata->addsymresacks && ! handlednonbinarysymmetry )
      {
         int k;

         SCIP_CALL( adaptSymmetryDataSST(scip, propdata->perms, modifiedperms, propdata->nperms,
               propdata->permvars, modifiedpermvars, propdata->npermvars, lexorder, nvarslexorder) );

         for (k = 0; k < npermsincomp; ++k)
         {
            SCIP_CONS* cons;
            char name[SCIP_MAXSTRLEN];

            /* skip permutations that have been used to build an orbitope */
            if ( permused[k] )
               continue;

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_comp%d_perm%d", i, k);

            SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, modifiedperms[propdata->components[propdata->componentbegins[i] + k]],
                  modifiedpermvars, propdata->npermvars, FALSE,
                  propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons));

            /* do not release constraint here - will be done later */
            propdata->genorbconss[propdata->ngenorbconss++] = cons;
            ++propdata->nsymresacks;

            if ( ! propdata->componentblocked[i] )
            {
               propdata->componentblocked[i] |= SYM_HANDLETYPE_SYMBREAK;
               ++propdata->ncompblocked;
            }

            SCIPdebugMsg(scip, "  add symresack for permutation %d of component %d adapted to suborbitope lexorder\n", k, i);
         }
      }

   FREELOOPMEMORY:
      SCIPfreeBlockMemoryArrayNull(scip, &lexorder, maxnvarslexorder);

      SCIPfreeBufferArrayNull(scip, &firstvaridxpercolor);
      SCIPfreeBufferArrayNull(scip, &chosencomppercolor);
      SCIPfreeBlockMemoryArrayNull(scip, &compcolorbegins, ncompcolors + 1);
      SCIPfreeBlockMemoryArrayNull(scip, &graphcompbegins, ngraphcomponents + 1);
      SCIPfreeBlockMemoryArrayNull(scip, &graphcomponents, propdata->npermvars);
      SCIPfreeBufferArrayNull(scip, &permused);
      SCIPfreeBufferArrayNull(scip, &usedperms);
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "total number of added (sub-)orbitopes: %d\n", norbitopes);
   SCIPdebugMsg(scip, "total number of added strong sbcs: %d\n", nstrongsbcs);
   SCIPdebugMsg(scip, "total number of added weak sbcs: %d\n", nweaksbcs);
#endif

   SCIPfreeBufferArray(scip, &nvarsincomponent);

   SCIPfreeBufferArray(scip, &modifiedpermvars);
   for (i = propdata->nperms - 1; i >= 0; --i)
   {
      SCIPfreeBufferArray(scip, &modifiedperms[i]);
   }
   SCIPfreeBufferArray(scip, &modifiedperms);
   SCIPfreeBufferArray(scip, &genorder);

   return SCIP_OKAY;
}


/*
 * Functions for symmetry constraints
 */


/** sorts orbitope vars matrix such that rows are sorted increasingly w.r.t. minimum variable index in row;
 *  columns are sorted such that first row is sorted increasingly w.r.t. variable indices
 */
static
SCIP_RETCODE SCIPsortOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 orbitopevaridx,     /**< variable index matrix of orbitope */
   SCIP_VAR***           vars,               /**< variable matrix of orbitope */
   int                   nrows,              /**< number of binary rows of orbitope */
   int                   ncols               /**< number of columns of orbitope */
   )
{
   SCIP_VAR** sortedrow;
   int* colorder;
   int* idcs;
   int arrlen;
   int minrowidx = INT_MAX;
   int minrow = INT_MAX;
   int i;
   int j;

   assert( scip != NULL );
   assert( orbitopevaridx != NULL );
   assert( vars != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );

   arrlen = MAX(nrows, ncols);
   SCIP_CALL( SCIPallocBufferArray(scip, &idcs, arrlen) );

   /* detect minimum index per row */
   for (i = 0; i < nrows; ++i)
   {
      int idx;

      idcs[i] = INT_MAX;

      for (j = 0; j < ncols; ++j)
      {
         idx = orbitopevaridx[i][j];

         if ( idx < idcs[i] )
            idcs[i] = idx;

         if ( idx < minrowidx )
         {
            minrowidx = idx;
            minrow = i;
         }
      }
   }

   /* sort rows increasingly w.r.t. minimum variable indices */
   SCIPsortIntPtr(idcs, (void**) vars, nrows);

   /* sort columns increasingly w.r.t. variable indices of first row */
   SCIP_CALL( SCIPallocBufferArray(scip, &colorder, ncols) );
   for (j = 0; j < ncols; ++j)
   {
      idcs[j] = orbitopevaridx[minrow][j];
      colorder[j] = j;
   }

   /* sort columns of first row and store new column order */
   SCIPsortIntIntPtr(idcs, colorder, (void**) vars[0], ncols);

   /* adapt rows 1, ..., nrows - 1 to new column order*/
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedrow, ncols) );
   for (i = 1; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
         sortedrow[j] = vars[i][colorder[j]];
      for (j = 0; j < ncols; ++j)
         vars[i][j] = sortedrow[j];
   }

   SCIPfreeBufferArray(scip, &sortedrow);
   SCIPfreeBufferArray(scip, &colorder);
   SCIPfreeBufferArray(scip, &idcs);

   return SCIP_OKAY;
}


/** checks whether components of the symmetry group can be completely handled by orbitopes */
static
SCIP_RETCODE detectOrbitopes(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  components,         /**< array containing components of symmetry group */
   int*                  componentbegins,    /**< array containing begin positions of components in components array */
   int                   ncomponents         /**< number of components */
   )
{
   SCIP_VAR** permvars;
   int** perms;
   int npermvars;
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( ncomponents > 0 );
   assert( propdata->nperms >= 0 );

   /* exit if no symmetry is present */
   if ( propdata->nperms == 0 )
      return SCIP_OKAY;

   assert( propdata->nperms > 0 );
   assert( propdata->perms != NULL );
   assert( propdata->nbinpermvars >= 0 );
   assert( propdata->npermvars >= 0 );
   assert( propdata->permvars != NULL );

   /* exit if no symmetry on binary variables is present */
   if ( propdata->nbinpermvars == 0 )
   {
      assert( ! propdata->binvaraffected );
      return SCIP_OKAY;
   }

   perms = propdata->perms;
   npermvars = propdata->npermvars;
   permvars = propdata->permvars;

   /* iterate over components */
   for (i = 0; i < ncomponents; ++i)
   {
      SCIP_VAR*** vars;
      SCIP_VAR*** varsallocorder;
      SCIP_CONS* cons;
      SCIP_Shortbool* rowisbinary;
      SCIP_Bool isorbitope = TRUE;
      SCIP_Bool infeasibleorbitope;
      int** orbitopevaridx;
      int* columnorder;
      int npermsincomponent;
      int ntwocyclescomp = INT_MAX;
      int nbincyclescomp = INT_MAX;
      int* nusedelems;
      int j;
      int cnt;

      /* orbitopes are detected first, so no component should be blocked */
      assert( ! propdata->componentblocked[i] );

      /* get properties of permutations */
      npermsincomponent = componentbegins[i + 1] - componentbegins[i];
      assert( npermsincomponent > 0 );
      for (j = componentbegins[i]; j < componentbegins[i + 1]; ++j)
      {
         int ntwocyclesperm = 0;
         int nbincyclesperm = 0;

         SCIP_CALL( SCIPisInvolutionPerm(perms[components[j]], permvars, npermvars,
               &ntwocyclesperm, &nbincyclesperm, FALSE) );

         if ( ntwocyclescomp == 0 )
         {
            isorbitope = FALSE;
            break;
         }

         /* if we are checking the first permutation */
         if ( ntwocyclescomp == INT_MAX )
         {
            ntwocyclescomp = ntwocyclesperm;
            nbincyclescomp = nbincyclesperm;

            /* if there are no binary rows */
            if ( nbincyclescomp == 0 )
            {
               isorbitope = FALSE;
               break;
            }
         }

         /* no or different number of 2-cycles or not all vars binary: permutations cannot generate orbitope */
         if ( ntwocyclescomp != ntwocyclesperm || nbincyclesperm != nbincyclescomp )
         {
            isorbitope = FALSE;
            break;
         }
      }

      /* if no orbitope was detected */
      if ( ! isorbitope )
         continue;
      assert( ntwocyclescomp > 0 );
      assert( ntwocyclescomp < INT_MAX );

      /* iterate over permutations and check whether for each permutation there exists
       * another permutation whose 2-cycles intersect pairwise in exactly one element */

      /* orbitope matrix for indices of variables in permvars array */
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx, ntwocyclescomp) );
      for (j = 0; j < ntwocyclescomp; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx[j], npermsincomponent + 1) ); /*lint !e866*/
      }

      /* order of columns of orbitopevaridx */
      SCIP_CALL( SCIPallocBufferArray(scip, &columnorder, npermsincomponent + 1) );
      for (j = 0; j < npermsincomponent + 1; ++j)
         columnorder[j] = npermsincomponent + 2;

      /* count how often an element was used in the potential orbitope */
      SCIP_CALL( SCIPallocClearBufferArray(scip, &nusedelems, npermvars) );

      /* store whether a row of the potential orbitope contains only binary variables */
      SCIP_CALL( SCIPallocClearBufferArray(scip, &rowisbinary, ntwocyclescomp) );

      /* check if the permutations fulfill properties of an orbitope */
      SCIP_CALL( checkTwoCyclePermsAreOrbitope(scip, permvars, npermvars, perms,
            &(components[componentbegins[i]]), ntwocyclescomp, npermsincomponent,
            orbitopevaridx, columnorder, nusedelems, NULL, rowisbinary, &isorbitope, NULL) );

      if ( ! isorbitope )
         goto FREEDATASTRUCTURES;

      /* we have found a potential orbitope, prepare data for orbitope conshdlr */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbincyclescomp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varsallocorder, nbincyclescomp) );
      cnt = 0;
      for (j = 0; j < ntwocyclescomp; ++j)
      {
         if ( ! rowisbinary[j] )
            continue;

         SCIP_CALL( SCIPallocBufferArray(scip, &vars[cnt], npermsincomponent + 1) ); /*lint !e866*/
         varsallocorder[cnt] = vars[cnt]; /* to ensure that we can free the buffer in reverse order */
         ++cnt;
      }
      assert( cnt == nbincyclescomp );

      /* prepare variable matrix (reorder columns of orbitopevaridx) */
      infeasibleorbitope = FALSE;
      SCIP_CALL( SCIPgenerateOrbitopeVarsMatrix(scip, &vars, ntwocyclescomp, npermsincomponent + 1, permvars,
            npermvars, orbitopevaridx, columnorder, nusedelems, rowisbinary, &infeasibleorbitope, FALSE, NULL, NULL, NULL) );

      if ( ! infeasibleorbitope )
      {
         char name[SCIP_MAXSTRLEN];

         SCIPdebugMsg(scip, "found an orbitope of size %d x %d in component %d\n", ntwocyclescomp, npermsincomponent + 1, i);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "orbitope_component%d", i);

         /* to ensure same orbitope is added if different sets of generators are found */
         SCIP_CALL( SCIPsortOrbitope(scip, orbitopevaridx, vars, nbincyclescomp, npermsincomponent + 1) );

         SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, vars, SCIP_ORBITOPETYPE_FULL,
               nbincyclescomp, npermsincomponent + 1, propdata->usedynamicprop, FALSE, FALSE, FALSE,
               propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         propdata->genorbconss[propdata->ngenorbconss++] = cons;
         ++propdata->norbitopes;

         propdata->componentblocked[i] |= SYM_HANDLETYPE_SYMBREAK;
      }

      /* free data structures */
      for (j = nbincyclescomp - 1; j >= 0; --j)
      {
         SCIPfreeBufferArray(scip, &varsallocorder[j]);
      }
      SCIPfreeBufferArray(scip, &varsallocorder);
      SCIPfreeBufferArray(scip, &vars);

   FREEDATASTRUCTURES:
      SCIPfreeBufferArray(scip, &rowisbinary);
      SCIPfreeBufferArray(scip, &nusedelems);
      SCIPfreeBufferArray(scip, &columnorder);
      for (j = ntwocyclescomp - 1; j >= 0; --j)
      {
         SCIPfreeBufferArray(scip, &orbitopevaridx[j]);
      }
      SCIPfreeBufferArray(scip, &orbitopevaridx);
   }

   return SCIP_OKAY;
}


/** update symmetry information of conflict graph */
static
SCIP_RETCODE updateSymInfoConflictGraphSST(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_VAR**            graphvars,          /**< variables encoded in conflict graph (either all vars or permvars) */
   int                   ngraphvars,         /**< number of nodes/vars in conflict graph */
   SCIP_VAR**            permvars,           /**< variables considered in permutations */
   int                   npermvars,          /**< number of permvars */
   SCIP_Bool             onlypermvars,       /**< whether conflict graph contains only permvars */
   SCIP_HASHMAP*         varmap,             /**< map from graphvar to node label in conflict graph
                                              *   (or NULL if onlypermvars == TRUE) */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int                   norbits             /**< number of non-trivial orbits */
   )
{
   int i;
   int j;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( graphvars != NULL );
   assert( ngraphvars > 0 );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( onlypermvars || varmap != NULL );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( norbits >= 0 );

   /* initialize/reset variable information of nodes in conflict graph */
   for (i = 0; i < ngraphvars; ++i)
   {
      SCIP_NODEDATA* nodedata;

      nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, i);

      /* possibly create node data */
      if ( nodedata == NULL )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &nodedata) );
         nodedata->var = graphvars[i];
         nodedata->active = TRUE;
      }

      /* (re-)set node data */
      nodedata->orbitidx = -1;
      nodedata->nconflictinorbit = 0;
      nodedata->orbitsize = -1;
      nodedata->posinorbit = -1;

      /* set node data */
      SCIPdigraphSetNodeData(conflictgraph, (void*) nodedata, i);
   }

   /* add orbit information to nodes of conflict graph */
   for (j = 0; j < norbits; ++j)
   {
      int posinorbit = 0;
      int orbitsize;

      orbitsize = orbitbegins[j + 1] - orbitbegins[j];
      assert( orbitsize >= 0 );

      for (i = orbitbegins[j]; i < orbitbegins[j + 1]; ++i)
      {
         SCIP_NODEDATA* nodedata;
         SCIP_VAR* var;
         int pos;

         /* get variable and position in conflict graph */
         if ( onlypermvars )
         {
            pos = orbits[i];
            var = permvars[pos];
         }
         else
         {
            var = permvars[orbits[i]];
            assert( var != NULL );

            assert( SCIPhashmapExists(varmap, var) );
            pos = SCIPhashmapGetImageInt(varmap, var);
            assert( pos != INT_MAX );
         }

         /* get node data */
         nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, pos);
         assert( nodedata != NULL );
         assert( nodedata->var == var );

         nodedata->orbitidx = j;
         nodedata->orbitsize = orbitsize;
         nodedata->posinorbit = posinorbit++;
      }
   }

   /* add information on number of conflicts within orbit to conflict graph */
   for (i = 0; i < ngraphvars; ++i)
   {
      SCIP_NODEDATA* nodedata;
      SCIP_NODEDATA* nodedataconflict;
      int* conflictvaridx;
      int nconflictinorbit = 0;
      int curorbit;

      nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, i);
      conflictvaridx = SCIPdigraphGetSuccessors(conflictgraph, i);

      assert( nodedata != NULL );
      assert( nodedata->nconflictinorbit == 0 );
      assert( conflictvaridx != NULL || SCIPdigraphGetNSuccessors(conflictgraph, i) == 0 );

      curorbit = nodedata->orbitidx;

      /* i-th variable is fixed by all permutations */
      if ( curorbit == -1 )
      {
         nodedata->nconflictinorbit = 0;
         continue;
      }

      /* get conflicts in orbit by couting the active neighbors of i in the same orbit */
      for (j = 0; j < SCIPdigraphGetNSuccessors(conflictgraph, i); ++j)
      {
         assert( conflictvaridx != NULL );
         nodedataconflict = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, conflictvaridx[j]);
         assert( nodedataconflict != NULL );

         if ( nodedataconflict->active && nodedataconflict->orbitidx == curorbit )
            ++nconflictinorbit;
      }

      nodedata->nconflictinorbit = nconflictinorbit;
   }

   return SCIP_OKAY;
}


/** create conflict graph either for symmetric or for all variables
 *
 *  This routine just creates the graph, but does not add (symmetry) information to its nodes.
 *  This has to be done separately by the routine updateSymInfoConflictGraphSST().
 */
static
SCIP_RETCODE createConflictGraphSST(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_DIGRAPH**        conflictgraph,      /**< pointer to store conflict graph */
   SCIP_VAR**            graphvars,          /**< variables encoded in conflict graph */
   int                   ngraphvars,         /**< number of vars encoded in conflict graph */
   SCIP_Bool             onlypermvars,       /**< whether conflict graph contains only permvars */
   SCIP_HASHMAP*         permvarmap,         /**< map of variables to indices in permvars array (or NULL) */
   SCIP_Bool*            success             /**< pointer to store whether conflict graph could be created successfully */
   )
{
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   SCIP_VAR** setppcvars;
   SCIP_CONS* cons;
   int nsetppcconss;
   int nsetppcvars;
   int nodei;
   int nodej;
   int c;
   int i;
   int j;
   int nedges = 0;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( graphvars != NULL );
   assert( ngraphvars > 0 );
   assert( success != NULL );

   *success = FALSE;

   /* get setppcconss for creating conflict graph */
   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   if ( setppcconshdlr == NULL )
   {
      SCIPdebugMsg(scip, "Could not find setppc conshdlr --> construction of conflict graph aborted.\n");
      return SCIP_OKAY;
   }
   assert( setppcconshdlr != NULL );

   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);
   if ( nsetppcconss == 0 )
   {
      SCIPdebugMsg(scip, "No setppc constraints present --> construction of conflict graph aborted.\n");
      return SCIP_OKAY;
   }

   /* construct conflict graph */
   SCIP_CALL( SCIPcreateDigraph(scip, conflictgraph, ngraphvars) );
   *success = TRUE;

   SCIPdebugMsg(scip, "Construction of conflict graph:\n");

   for (c = 0; c < nsetppcconss; ++c)
   {
      cons = setppcconss[c];
      assert( cons != NULL );

      /* skip covering constraints */
      if ( SCIPgetTypeSetppc(scip, cons) == SCIP_SETPPCTYPE_COVERING )
         continue;

      setppcvars = SCIPgetVarsSetppc(scip, cons);
      nsetppcvars = SCIPgetNVarsSetppc(scip, cons);
      assert( setppcvars != NULL );
      assert( nsetppcvars > 0 );

      SCIPdebugMsg(scip, "\tAdd edges for constraint %s.\n", SCIPconsGetName(cons));

      /* iterate over pairs of variables in constraint and add bidirected arc
       * if both are affected by a symmetry or active */
      for (i = 0; i < nsetppcvars; ++i)
      {
         if ( onlypermvars )
         {
            nodei = SCIPhashmapGetImageInt(permvarmap, setppcvars[i]);

            /* skip variables that are not affected by symmetry */
            if ( nodei == INT_MAX )
               continue;
         }
         else
         {
            nodei = SCIPvarGetProbindex(setppcvars[i]);

            /* skip inactive variables */
            if ( nodei < 0 )
               continue;
         }

         for (j = i + 1; j < nsetppcvars; ++j)
         {
            if ( onlypermvars )
            {
               nodej = SCIPhashmapGetImageInt(permvarmap, setppcvars[j]);

               /* skip variables that are not affected by symmetyr */
               if ( nodej == INT_MAX )
                  continue;
            }
            else
            {
               nodej = SCIPvarGetProbindex(setppcvars[j]);

               /* skip inactive variables */
               if ( nodej < 0 )
                  continue;
            }

            SCIP_CALL( SCIPdigraphAddArcSafe(*conflictgraph, nodei, nodej, NULL) );
            SCIP_CALL( SCIPdigraphAddArcSafe(*conflictgraph, nodej, nodei, NULL) );
            ++nedges;
         }
      }
   }
   SCIPdebugMsg(scip, "Construction of conflict graph terminated; %d conflicts detected.\n", nedges);

   return SCIP_OKAY;
}


/** frees conflict graph */
static
SCIP_RETCODE freeConflictGraphSST(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_DIGRAPH**        conflictgraph,      /**< conflict graph */
   int                   nnodes              /**< number of nodes in conflict graph */
   )
{
   int i;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( *conflictgraph != NULL );
   assert( nnodes > 0 );

   /* free node data */
   for (i = 0; i < nnodes; ++i)
   {
      SCIP_NODEDATA* nodedata;

      /* get node data */
      nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(*conflictgraph, i);

      /* free node data (might not have been allocated if all components are already blocked) */
      if ( nodedata != NULL )
      {
         SCIPfreeBlockMemory(scip, &nodedata);
      }
   }

   /* free conflict graph */
   SCIPdigraphFree(conflictgraph);

   return SCIP_OKAY;
}


/** adds symresack constraints */
static
SCIP_RETCODE addSymresackConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< symmetry breaking propagator */
   int*                  components,         /**< array containing components of symmetry group */
   int*                  componentbegins,    /**< array containing begin positions of components in components array */
   int                   ncomponents         /**< number of components */
   )
{ /*lint --e{641}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** permvars;
   SCIP_Bool conssaddlp;
   int** modifiedperms = NULL;
   SCIP_VAR** modifiedpermvars = NULL;
   int** perms;
   int nsymresackcons = 0;
   int npermvars;
   int nperms;
   int i;
   int p;

   assert( scip != NULL );
   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->npermvars >= 0 );
   assert( propdata->nbinpermvars >= 0 );

   /* if no symmetries on binary variables are present */
   if ( propdata->nbinpermvars == 0 )
   {
      assert( propdata->binvaraffected == 0 );
      return SCIP_OKAY;
   }

   perms = propdata->perms;
   nperms = propdata->nperms;
   permvars = propdata->permvars;
   npermvars = propdata->npermvars;
   conssaddlp = propdata->conssaddlp;

   assert( nperms <= 0 || perms != NULL );
   assert( permvars != NULL );
   assert( npermvars > 0 );

   /* adapt natural variable order to a variable order that is compatible with Schreier Sims constraints */
   if ( propdata->nleaders > 0 && ISSSTBINACTIVE(propdata->sstleadervartype) )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &modifiedperms, nperms) );
      for (p = 0; p < nperms; ++p)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &modifiedperms[p], npermvars) );
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &modifiedpermvars, npermvars) );

      for (i = 0; i < npermvars; ++i)
         modifiedpermvars[i] = permvars[i];

      SCIP_CALL( adaptSymmetryDataSST(scip, perms, modifiedperms, nperms, permvars, modifiedpermvars, npermvars,
            propdata->leaders, propdata->nleaders) );
   }

   /* if components have not been computed */
   if ( ncomponents == -1 )
   {
      assert( ! propdata->ofenabled );
      assert( ! propdata->detectorbitopes );
      assert( ! propdata->sstenabled );

      /* loop through perms and add symresack constraints */
      for (p = 0; p < propdata->nperms; ++p)
      {
         SCIP_CONS* cons;
         char name[SCIP_MAXSTRLEN];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symbreakcons_perm%d", p);

         SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, perms[p], permvars, npermvars, FALSE,
                  conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         propdata->genorbconss[propdata->ngenorbconss++] = cons;
         ++propdata->nsymresacks;
         ++nsymresackcons;
      }
   }
   else
   {
      /* loop through components */
      for (i = 0; i < ncomponents; ++i)
      {
         SCIP_Bool sstcompatible = TRUE;

         if ( ISSSTINTACTIVE(propdata->sstleadervartype)
            || ISSSTIMPLINTACTIVE(propdata->sstleadervartype)
            || ISSSTCONTACTIVE(propdata->sstleadervartype) )
            sstcompatible = FALSE;

         /* skip components that were treated by incompatible symmetry handling techniques */
         if ( (propdata->componentblocked[i] & SYM_HANDLETYPE_SYMBREAK) != 0
            || (propdata->componentblocked[i] & SYM_HANDLETYPE_ORBITALFIXING) != 0
            || ((propdata->componentblocked[i] & SYM_HANDLETYPE_SST) != 0 && ! sstcompatible) )
            continue;

         /* loop through perms in component i and add symresack constraints */
         for (p = componentbegins[i]; p < componentbegins[i + 1]; ++p)
         {
            SCIP_CONS* cons;
            int permidx;
            char name[SCIP_MAXSTRLEN];

            permidx = components[p];

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symbreakcons_component%d_perm%d", i, permidx);

            /* adapt permutation to leader */
            if ( propdata->nleaders > 0 && ISSSTBINACTIVE(propdata->sstleadervartype) )
            {
               assert( (propdata->componentblocked[i] & SYM_HANDLETYPE_SST) != 0 );
               assert( modifiedperms != NULL );
               assert( modifiedpermvars != NULL );

               SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, modifiedperms[permidx], modifiedpermvars, npermvars, FALSE,
                     conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, perms[permidx], permvars, npermvars, FALSE,
                     conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            }
            propdata->componentblocked[i] |= SYM_HANDLETYPE_SYMBREAK;
            SCIP_CALL( SCIPaddCons(scip, cons) );

            /* do not release constraint here - will be done later */
            propdata->genorbconss[propdata->ngenorbconss++] = cons;
            ++propdata->nsymresacks;
            ++nsymresackcons;
         }
      }
   }

   if ( propdata->nleaders > 0 && ISSSTBINACTIVE(propdata->sstleadervartype) )
   {
      assert( modifiedperms != NULL );
      assert( modifiedpermvars != NULL );

      SCIPfreeBufferArray(scip, &modifiedpermvars);
      for (p = nperms - 1; p >= 0; --p)
      {
         SCIPfreeBufferArray(scip, &modifiedperms[p]);
      }
      SCIPfreeBufferArray(scip, &modifiedperms);
   }

   SCIPdebugMsg(scip, "Added %d symresack constraints.\n", nsymresackcons);

   return SCIP_OKAY;
}


/** add Schreier Sims constraints for a specific orbit and update Schreier Sims table */
static
SCIP_RETCODE addSSTConssOrbitAndUpdateSST(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph or NULL if useconflictgraph == FALSE */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   SCIP_VAR**            permvars,           /**< permvars array */
   int*                  orbits,             /**< symmetry orbits */
   int*                  orbitbegins,        /**< array storing begin position for each orbit */
   int                   orbitidx,           /**< index of orbit for Schreier Sims constraints */
   int                   orbitleaderidx,     /**< index of leader variable for Schreier Sims constraints */
   SCIP_Shortbool*       orbitvarinconflict, /**< indicator whether orbitvar is in conflict with orbit leader */
   int                   norbitvarinconflict, /**< number of variables in conflict with orbit leader */
   int*                  nchgbds,            /**< pointer to store number of bound changes (or NULL) */
   SCIP_Bool             useconflictgraph    /**< whether conflict graph shall be used */
   )
{ /*lint --e{613,641}*/
   SCIP_CONS* cons;
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   int orbitsize;
   int posleader;
   int poscur;
   int ncuts = 0;
   SCIP_Bool addcuts = FALSE;
   int i;
#ifndef NDEBUG
   int j;
#endif

   assert( scip != NULL );
   assert( conflictgraph != NULL || ! useconflictgraph );
   assert( propdata != NULL );
   assert( permvars != NULL );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( orbitidx >= 0 );
   assert( orbitleaderidx >= 0 );
   assert( orbitvarinconflict != NULL || ! useconflictgraph );
   assert( norbitvarinconflict >= 0 );
   assert( nchgbds != NULL );

   orbitsize = orbitbegins[orbitidx + 1] - orbitbegins[orbitidx];

   /* variables in conflict with leader are fixed and not treated by a cut; trailing -1 to not count the leader */
   if ( propdata->sstaddcuts )
      addcuts = TRUE;
   else if ( propdata->sstleaderrule == SCIP_LEADERRULE_MAXCONFLICTSINORBIT
      || propdata->sstleaderrule == SCIP_LEADERRULE_MAXCONFLICTS
      || propdata->ssttiebreakrule == SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT )
      addcuts = propdata->addconflictcuts;

   if ( addcuts )
      ncuts = orbitsize - norbitvarinconflict - 1;

   /* (re-)allocate memory for Schreier Sims constraints and leaders */
   if ( ncuts > 0 )
   {
      if ( propdata->nsstconss == 0 )
      {
         assert( propdata->sstconss == NULL );
         assert( propdata->maxnsstconss == 0 );
         propdata->maxnsstconss = 2 * ncuts;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->sstconss), propdata->maxnsstconss) );
      }
      else if ( propdata->nsstconss + ncuts > propdata->maxnsstconss )
      {
         int newsize;

         newsize = SCIPcalcMemGrowSize(scip, propdata->maxnsstconss + 2 * ncuts);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(propdata->sstconss),
               propdata->maxnsstconss, newsize) );
         propdata->maxnsstconss = newsize;
      }
   }

   if ( propdata->nleaders == 0 )
   {
      propdata->maxnleaders = MIN(propdata->nperms, propdata->npermvars);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->leaders), propdata->maxnleaders) );
   }
   assert( propdata->nleaders < propdata->maxnleaders );

   /* add Schreier Sims constraints vars[0] >= vars[1], where vars[0] is always the leader */
   posleader = orbitbegins[orbitidx] + orbitleaderidx;
   vars[0] = permvars[orbits[posleader]];
   vals[0] = -1.0;
   vals[1] = 1.0;
   propdata->leaders[propdata->nleaders++] = orbits[posleader];
   *nchgbds = 0;
   for (i = 0, poscur = orbitbegins[orbitidx]; i < orbitsize; ++i, ++poscur)
   {
      if ( i == orbitleaderidx )
      {
         assert( orbitvarinconflict == NULL || ! orbitvarinconflict[i] );
         continue;
      }

      vars[1] = permvars[orbits[poscur]];
#ifndef NDEBUG
      for (j = 0; j < propdata->nleaders - 1; ++j)
      {
         assert( propdata->leaders[j] != orbits[poscur] );
      }
#endif

      /* if the i-th variable in the orbit is in a conflict with the leader, fix it to 0 */
      if ( useconflictgraph )
      {
         if ( orbitvarinconflict[i] )
         {
            assert( SCIPvarIsBinary(vars[1]) );
            assert( SCIPvarGetLbLocal(vars[1]) < 0.5 );
            assert( useconflictgraph );

            /* if variable is fixed */
            if ( SCIPvarGetUbLocal(vars[1]) > 0.5 )
            {
               SCIP_NODEDATA* nodedata;

               SCIP_CALL( SCIPchgVarUb(scip, vars[1], 0.0) );
               ++(*nchgbds);

               /* deactivate the fixed variable (cannot contribute to a conflict anymore) */
               nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, orbits[poscur]);
               assert( nodedata != NULL );
               assert( nodedata->active );

               nodedata->active = FALSE;
            }

            /* reset value */
            orbitvarinconflict[i] = FALSE;
         }
         else if ( addcuts )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "SSTcut_%d_%d", orbits[posleader], orbits[poscur]);
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, - SCIPinfinity(scip), 0.0,
                  FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
            propdata->sstconss[propdata->nsstconss++] = cons;
         }
      }
      else if ( addcuts )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "SSTcut_%d_%d", orbits[posleader], orbits[poscur]);
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, - SCIPinfinity(scip), 0.0,
               FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         propdata->sstconss[propdata->nsstconss++] = cons;
      }
   }

   return SCIP_OKAY;
}


/** selection rule of next orbit/leader in orbit for Schreier Sims constraints */
static
SCIP_RETCODE selectOrbitLeaderSSTConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph or NULL if useconflictgraph == FALSE */
   SCIP_VAR**            graphvars,          /**< variables encoded in conflict graph */
   int                   ngraphvars,         /**< number of variables encoded in conflict graph */
   SCIP_HASHMAP*         varmap,             /**< map from variable to node label in conflict graph or NULL if useconflictgraph == FALSE  */
   SCIP_VAR**            permvars,           /**< vars encoded in a permutation */
   int                   npermvars,          /**< number of vars in a permutation */
   int*                  orbits,             /**< orbits of stabilizer subgroup */
   int*                  orbitbegins,        /**< array storing the begin position of each orbit in orbits */
   int                   norbits,            /**< number of orbits */
   int                   leaderrule,         /**< rule to select leader */
   int                   tiebreakrule,       /**< tie break rule to select leader */
   SCIP_VARTYPE          leadervartype,      /**< variable type of leader */
   int*                  orbitidx,           /**< pointer to index of selected orbit */
   int*                  leaderidx,          /**< pointer to leader in orbit */
   SCIP_Shortbool*       orbitvarinconflict, /**< array to store whether a var in the orbit is conflicting with leader */
   int*                  norbitvarinconflict, /**< pointer to store number of vars in the orbit in conflict with leader */
   SCIP_Bool             useconflictgraph,   /**< whether conflict graph shall be used */
   SCIP_Bool*            success             /**< pointer to store whether orbit cut be selected successfully */
   )
{
   SCIP_NODEDATA* nodedata;
   int* conflictvars;
   int nconflictvars = 0;
   int varidx;
   int orbitcriterion;
   int curcriterion = INT_MIN;
   int orbitsize;
   int i;
   SCIP_NODEDATA* neighbordata;
   int leader = -1;
   int j;

   assert( scip != NULL );
   assert( conflictgraph != NULL || ! useconflictgraph );
   assert( graphvars != NULL );
   assert( ngraphvars > 0 );
   assert( varmap != NULL || ! useconflictgraph );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( norbits > 0 );
   assert( orbitidx != NULL );
   assert( leaderidx != NULL );
   assert( orbitvarinconflict != NULL || ! useconflictgraph );
   assert( norbitvarinconflict != NULL );
   assert( success != NULL );

   *orbitidx = 0;
   *leaderidx = 0;
   *norbitvarinconflict = 0;
   *success = FALSE;

   /* terminate if leader or tiebreak rule cannot be checked */
   if ( ! useconflictgraph && (leaderrule == (int) SCIP_LEADERRULE_MAXCONFLICTS
         || leaderrule == (int) SCIP_LEADERRULE_MAXCONFLICTSINORBIT
         || tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT) )
      return SCIP_OKAY;

   /* select the leader and its orbit */
   if ( leaderrule == (int) SCIP_LEADERRULE_FIRSTINORBIT || leaderrule == (int) SCIP_LEADERRULE_LASTINORBIT )
   {
      orbitcriterion = INT_MIN;

      /* iterate over orbits and select the first one that meets the tiebreak rule */
      for (i = 0; i < norbits; ++i)
      {
         /* skip orbits containing vars different to the leader's vartype */
         if ( SCIPvarGetType(permvars[orbits[orbitbegins[i]]]) != leadervartype )
            continue;

         if ( tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MINORBIT )
            curcriterion = orbitbegins[i] - orbitbegins[i + 1];
         else if ( tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MAXORBIT )
            curcriterion = orbitbegins[i + 1] - orbitbegins[i];
         else
         {
            /* get first or last active variable in orbit */
            if ( leaderrule == (int) SCIP_LEADERRULE_FIRSTINORBIT )
            {
               int cnt = orbitbegins[i];

               do
               {
                  varidx = SCIPvarGetProbindex(permvars[orbits[cnt++]]);
               }
               while ( varidx == -1 && cnt < orbitbegins[i + 1]);
            }
            else
            {
               int cnt = orbitbegins[i + 1] - 1;

               do
               {
                  varidx = SCIPvarGetProbindex(permvars[orbits[cnt--]]);
               }
               while ( varidx == -1 && cnt >= orbitbegins[i]);
            }

            /* skip inactive variables */
            if ( varidx == -1 )
               continue;

            nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, varidx);
            assert( nodedata != NULL );
            assert( nodedata->orbitidx == i );

            if ( nodedata->nconflictinorbit > 0 )
               curcriterion = nodedata->nconflictinorbit;
         }

         /* update selected orbit */
         if ( curcriterion > orbitcriterion )
         {
            orbitcriterion = curcriterion;
            *orbitidx = i;
            *success = TRUE;

            if ( leaderrule == (int) SCIP_LEADERRULE_FIRSTINORBIT )
               *leaderidx = 0;
            else
               *leaderidx = orbitbegins[i + 1] - orbitbegins[i] - 1;
         }
      }

      /* store variables in conflict with leader */
      if ( useconflictgraph )
      {
         leader = SCIPhashmapGetImageInt(varmap, permvars[orbits[orbitbegins[*orbitidx] + *leaderidx]]);
         assert( leader < SCIPdigraphGetNNodes(conflictgraph) );

         nconflictvars = SCIPdigraphGetNSuccessors(conflictgraph, leader);
      }

      if ( *success && tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT && nconflictvars > 0 )
      {
         SCIP_VAR* var;
         orbitsize = orbitbegins[*orbitidx + 1] - orbitbegins[*orbitidx];
         assert( useconflictgraph );
         assert( leader >= 0 && leader < npermvars );

         conflictvars = SCIPdigraphGetSuccessors(conflictgraph, leader);
         assert( conflictvars != NULL );
         assert( orbitvarinconflict != NULL );

         for (i = 0; i < orbitsize; ++i)
         {
            /* skip the leader */
            if ( i == *leaderidx )
               continue;

            var = permvars[orbits[orbitbegins[*orbitidx] + i]];

            for (j = 0; j < nconflictvars; ++j)
            {
               neighbordata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, conflictvars[j]);

               assert( neighbordata != NULL );

               if ( neighbordata->var == var && neighbordata->active )
               {
                  orbitvarinconflict[i] = TRUE;
                  *norbitvarinconflict += 1;
                  break;
               }
            }
         }
      }
   }
   else if ( useconflictgraph )
   {
      orbitcriterion = 0;

      /* iterate over variables and select the first one that meets the tiebreak rule */
      for (i = 0; i < ngraphvars; ++i)
      {
         /* skip vars different to the leader's vartype */
         if ( SCIPvarGetType(graphvars[i]) != leadervartype )
            continue;

         nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, i);
         assert( nodedata != NULL );

         /* skip variables not affected by symmetry */
         if ( nodedata->orbitidx == -1 )
            continue;

         if ( leaderrule == (int) SCIP_LEADERRULE_MAXCONFLICTSINORBIT )
            curcriterion = nodedata->nconflictinorbit;
         else
            curcriterion = SCIPdigraphGetNSuccessors(conflictgraph, i);

         if ( curcriterion > orbitcriterion )
         {
            orbitcriterion = curcriterion;
            *orbitidx = nodedata->orbitidx;
            *leaderidx = nodedata->posinorbit;
            *success = TRUE;
         }
      }

      /* store variables in conflict with leader */
      leader = SCIPhashmapGetImageInt(varmap, permvars[orbits[orbitbegins[*orbitidx] + *leaderidx]]);
      assert( leader < SCIPdigraphGetNNodes(conflictgraph) );
      assert( norbitvarinconflict != NULL );

      nconflictvars = SCIPdigraphGetNSuccessors(conflictgraph, leader);
      if ( *success && nconflictvars > 0 )
      {
         SCIP_VAR* var;
         assert( orbitvarinconflict != NULL );

         orbitsize = orbitbegins[*orbitidx + 1] - orbitbegins[*orbitidx];

         conflictvars = SCIPdigraphGetSuccessors(conflictgraph, leader);
         assert( conflictvars != NULL );

         for (i = 0; i < orbitsize; ++i)
         {
            /* skip the leader */
            if ( i == *leaderidx )
               continue;

            var = permvars[orbits[orbitbegins[*orbitidx] + i]];

            for (j = 0; j < nconflictvars; ++j)
            {
               neighbordata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(conflictgraph, conflictvars[j]);
               assert( neighbordata != NULL );

               if ( neighbordata->var == var && neighbordata->active )
               {
                  orbitvarinconflict[i] = TRUE;
                  *norbitvarinconflict += 1;
                  break;
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** add Schreier Sims constraints to the problem */
static
SCIP_RETCODE addSSTConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< datas of symmetry propagator */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{ /*lint --e{641}*/
   SCIP_DIGRAPH* conflictgraph = NULL;
   SCIP_HASHMAP* varmap = NULL;
   SCIP_VAR** vars;
   int nvars;

   SCIP_HASHMAP* permvarmap;
   SCIP_VAR** permvars;
   int** permstrans;
   int npermvars;
   int nmovedpermvars;
   int nmovedbinpermvars;
   int nmovedintpermvars;
   int nmovedimplintpermvars;
   int nmovedcontpermvars;
   int nperms;

   int* orbits;
   int* orbitbegins;
   int norbits;
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int ncomponents;
   unsigned* componentblocked;

   int orbitidx;
   int orbitleaderidx;
   SCIP_Shortbool* orbitvarinconflict = NULL;
   int norbitvarinconflict;
   SCIP_Shortbool* inactiveperms;
   int ninactiveperms;
   int posleader;
   int leaderrule;
   int tiebreakrule;
   int leadervartype;
   SCIP_VARTYPE selectedtype = SCIP_VARTYPE_CONTINUOUS;
   int nvarsselectedtype;
   SCIP_Bool conflictgraphcreated = FALSE;
   SCIP_Bool mixedcomponents;
   int* norbitleadercomponent;

   int c;
   int v;
   int p;

   assert( scip != NULL );
   assert( propdata != NULL );

   permvars = propdata->permvars;
   npermvars = propdata->npermvars;
   permvarmap = propdata->permvarmap;
   permstrans = propdata->permstrans;
   nperms = propdata->nperms;
   components = propdata->components;
   componentbegins = propdata->componentbegins;
   componentblocked = propdata->componentblocked;
   vartocomponent = propdata->vartocomponent;
   ncomponents = propdata->ncomponents;
   nmovedpermvars = propdata->nmovedpermvars;
   nmovedbinpermvars = propdata->nmovedbinpermvars;
   nmovedintpermvars = propdata->nmovedintpermvars;
   nmovedimplintpermvars = propdata->nmovedimplintpermvars;
   nmovedcontpermvars = propdata->nmovedcontpermvars;

   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( permvarmap != NULL );
   assert( permstrans != NULL );
   assert( nperms > 0 );
   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( vartocomponent != NULL );
   assert( ncomponents > 0 );
   assert( nmovedpermvars > 0 || ! propdata->ofenabled );
   assert( nmovedbinpermvars > 0 || ! propdata->ofenabled );

   leaderrule = propdata->sstleaderrule;
   tiebreakrule = propdata->ssttiebreakrule;
   leadervartype = propdata->sstleadervartype;
   mixedcomponents = propdata->sstmixedcomponents;

   /* if not already computed, get number of affected vars */
   if ( nmovedpermvars == -1 )
   {
      nmovedpermvars = 0;

      for (v = 0; v < npermvars; ++v)
      {
         for (p = 0; p < nperms; ++p)
         {
            if ( permstrans[v][p] != v )
            {
               ++nmovedpermvars;

               switch ( SCIPvarGetType(permvars[v]) )
               {
               case SCIP_VARTYPE_BINARY:
                  ++nmovedbinpermvars;
                  break;
               case SCIP_VARTYPE_INTEGER:
                  ++nmovedintpermvars;
                  break;
               case SCIP_VARTYPE_IMPLINT:
                  ++nmovedimplintpermvars;
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
               default:
                  ++nmovedcontpermvars;
               }
            }
         }
      }
   }
   propdata->nmovedbinpermvars = nmovedbinpermvars;
   propdata->nmovedintpermvars = nmovedintpermvars;
   propdata->nmovedimplintpermvars = nmovedimplintpermvars;
   propdata->nmovedcontpermvars = nmovedcontpermvars;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* determine the leader's vartype */
   nvarsselectedtype = 0;
   if ( ISSSTBINACTIVE(leadervartype) && nmovedbinpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_BINARY;
      nvarsselectedtype = nmovedbinpermvars;
   }

   if ( ISSSTINTACTIVE(leadervartype) && nmovedintpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_INTEGER;
      nvarsselectedtype = nmovedintpermvars;
   }

   if ( ISSSTIMPLINTACTIVE(leadervartype) && nmovedimplintpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_IMPLINT;
      nvarsselectedtype = nmovedimplintpermvars;
   }

   if ( ISSSTCONTACTIVE(leadervartype) && nmovedcontpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_CONTINUOUS;
      nvarsselectedtype = nmovedcontpermvars;
   }

   /* terminate if no variables of a possible leader type is affected */
   if ( nvarsselectedtype == 0 )
      return SCIP_OKAY;

   /* possibly create conflict graph; graph is not created if no setppc conss are present */
   if ( selectedtype == SCIP_VARTYPE_BINARY && (leaderrule == SCIP_LEADERRULE_MAXCONFLICTSINORBIT
         || leaderrule == SCIP_LEADERRULE_MAXCONFLICTS
         || tiebreakrule == SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT) )
   {
      SCIP_CALL( createConflictGraphSST(scip, &conflictgraph, vars, nvars, FALSE,
            permvarmap, &conflictgraphcreated) );
   }

   /* allocate data structures necessary for orbit computations and conflict graph */
   SCIP_CALL( SCIPallocBufferArray(scip, &inactiveperms, nperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );

   if ( conflictgraphcreated )
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &orbitvarinconflict, npermvars) );
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), nvars) );
      for (v = 0; v < nvars; ++v)
      {
         assert( ! SCIPhashmapExists(varmap, vars[v]) );
         SCIP_CALL( SCIPhashmapInsertInt(varmap, vars[v], v) );
      }
   }

   SCIPdebugMsg(scip, "Start selection of orbits and leaders for Schreier Sims constraints.\n");
   SCIPdebugMsg(scip, "orbitidx\tleaderidx\torbitsize\n");

   if ( nchgbds != NULL )
      *nchgbds = 0;

   /* initialize array indicating whether permutations shall not be considered for orbit permutations */
   for (p = 0; p < nperms; ++p)
      inactiveperms[p] = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &norbitleadercomponent, ncomponents) );
   for (c = 0; c < ncomponents; ++c)
      norbitleadercomponent[c] = 0;

   /* iterate over components and compute orbits */
   for (c = 0; c < ncomponents; ++c)
   {
      SCIP_Bool success = TRUE;

      if ( componentblocked[c] )
         continue;

      for (p = componentbegins[c]; p < componentbegins[c + 1]; ++p)
         inactiveperms[components[p]] = FALSE;
      ninactiveperms = nperms - componentbegins[c + 1] + componentbegins[c];

      /* as long as the stabilizer is non-trivial, add Schreier Sims constraints */
      while ( ninactiveperms < nperms )
      {
         int nchanges = 0;

         /* compute orbits w.r.t. active perms */
         SCIP_CALL( SCIPcomputeOrbitsFilterSym(scip, npermvars, permstrans, nperms, inactiveperms,
               orbits, orbitbegins, &norbits, components, componentbegins, vartocomponent,
               componentblocked, ncomponents, nmovedpermvars) );

         /* stop if we require pure components and a component contains variables of different types */
         if ( ! mixedcomponents )
         {
            for (p = 0; p < norbits; ++p)
            {
               /* stop if the first element of an orbits has the wrong vartype */
               if ( SCIPvarGetType(permvars[orbits[orbitbegins[p]]]) != selectedtype )
               {
                  success = FALSE;
                  break;
               }
            }
         }

         if ( ! success )
            break;

         /* update symmetry information of conflict graph */
         if ( conflictgraphcreated )
         {
            assert( conflictgraph != NULL );
            SCIP_CALL( updateSymInfoConflictGraphSST(scip, conflictgraph, vars, nvars, permvars, npermvars, FALSE,
                  varmap, orbits, orbitbegins, norbits) );
         }

         /* possibly adapt the leader and tie-break rule */
         if ( (leaderrule == SCIP_LEADERRULE_MAXCONFLICTSINORBIT || leaderrule == SCIP_LEADERRULE_MAXCONFLICTS)
            && ! conflictgraphcreated )
            leaderrule = SCIP_LEADERRULE_FIRSTINORBIT;
         if ( (leaderrule == SCIP_LEADERRULE_MAXCONFLICTSINORBIT || leaderrule == SCIP_LEADERRULE_MAXCONFLICTS)
            && selectedtype != SCIP_VARTYPE_BINARY )
            leaderrule = SCIP_LEADERRULE_FIRSTINORBIT;
         if ( tiebreakrule == SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT && ! conflictgraphcreated )
            tiebreakrule = SCIP_LEADERTIEBREAKRULE_MAXORBIT;
         if ( tiebreakrule == SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT && selectedtype != SCIP_VARTYPE_BINARY )
            tiebreakrule = SCIP_LEADERTIEBREAKRULE_MAXORBIT;

         /* select orbit and leader */
         SCIP_CALL( selectOrbitLeaderSSTConss(scip, conflictgraph, vars, nvars, varmap,
               permvars, npermvars, orbits, orbitbegins, norbits, propdata->sstleaderrule, propdata->ssttiebreakrule, selectedtype,
               &orbitidx, &orbitleaderidx, orbitvarinconflict, &norbitvarinconflict, conflictgraphcreated, &success) );

         if ( ! success )
            break;

         assert( 0 <= orbitidx && orbitidx < norbits );
         assert( 0 <= orbitleaderidx && orbitleaderidx < orbitbegins[orbitidx + 1] - orbitbegins[orbitidx] );
         SCIPdebugMsg(scip, "%d\t\t%d\t\t%d\n", orbitidx, orbitleaderidx, orbitbegins[orbitidx + 1] - orbitbegins[orbitidx]);

         /* add Schreier Sims constraints for the selected orbit and update Schreier Sims table */
         SCIP_CALL( addSSTConssOrbitAndUpdateSST(scip, conflictgraph, propdata, permvars,
               orbits, orbitbegins, orbitidx, orbitleaderidx, orbitvarinconflict, norbitvarinconflict, &nchanges, conflictgraphcreated) );

         ++norbitleadercomponent[propdata->vartocomponent[orbits[orbitbegins[orbitidx] + orbitleaderidx]]];

         if ( nchgbds != NULL )
            *nchgbds += nchanges;

         /* deactivate permutations that move the orbit leader */
         posleader = orbits[orbitbegins[orbitidx] + orbitleaderidx];
         for (p = 0; p < nperms; ++p)
         {
            if ( inactiveperms[p] )
               continue;

            if ( permstrans[posleader][p] != posleader )
            {
               inactiveperms[p] = TRUE;
               ++ninactiveperms;
            }
         }
      }

      for (p = componentbegins[c]; p < componentbegins[c + 1]; ++p)
         inactiveperms[components[p]] = TRUE;
   }

   /* if Schreier Sims constraints have been added, store that Schreier Sims has been used for this component */
   for (c = 0; c < ncomponents; ++c)
   {
      if ( norbitleadercomponent[c] > 0 )
         componentblocked[c] |= SYM_HANDLETYPE_SST;
   }
   SCIPfreeBufferArray(scip, &norbitleadercomponent);

   if ( conflictgraphcreated )
   {
      SCIPhashmapFree(&varmap);
      SCIPfreeBufferArray(scip, &orbitvarinconflict);
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);
   if ( conflictgraphcreated )
   {
      assert( conflictgraph != NULL );
      SCIP_CALL( freeConflictGraphSST(scip, &conflictgraph, nvars) );
   }
   SCIPfreeBufferArray(scip, &inactiveperms);

   return SCIP_OKAY;
}


/** finds problem symmetries */
static
SCIP_RETCODE tryAddSymmetryHandlingConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< symmetry breaking propagator */
   int*                  nchgbds,            /**< pointer to store number of bound changes (or NULL)*/
   SCIP_Bool*            earlyterm           /**< pointer to store whether we terminated early (or NULL) */
   )
{
   SCIP_PROPDATA* propdata;

   assert( prop != NULL );
   assert( scip != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->symconsenabled || propdata->sstenabled );

   /* if constraints have already been added */
   if ( propdata->triedaddconss && isSymmetryRecomputationRequired(scip, propdata) )
   {
      /* remove symmetry handling constraints to be prepared for a recomputation */
      SCIP_CALL( delSymConss(scip, propdata) );
      SCIP_CALL( freeSymmetryData(scip, propdata) );

      propdata->lastrestart = SCIPgetNRuns(scip);
      propdata->offoundreduction = FALSE;
   }
   else if ( propdata->triedaddconss )
   {
      assert( propdata->nperms > 0 );

      if ( earlyterm != NULL )
         *earlyterm = TRUE;

      return SCIP_OKAY;
   }

   /* possibly compute symmetry */
   if ( propdata->ofenabled && SCIPgetNBinVars(scip) > 1 )
   {
      SCIP_Bool oldsymconsenabled;

      oldsymconsenabled = propdata->symconsenabled;

      /* in the nonlinear case, all non-binary variables have to be fixed
         (fix non-binary potential branching variables)
      */
      if ( hasNonlinearConstraints(propdata) || propdata->symfixnonbinaryvars )
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY, SYM_SPEC_INTEGER | SYM_SPEC_REAL) );
      }
      else
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_REAL, SYM_SPEC_INTEGER) );
      }

      /* if there is no symmetry compatible with OF, check whether there are more general symmetries */
      if ( propdata->nperms == 0 && SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip) > 1 )
      {
         SCIP_CALL( freeSymmetryData(scip, propdata) );
         propdata->symconsenabled = oldsymconsenabled;
         propdata->ofenabled = FALSE;
         propdata->sstenabled = FALSE;

         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_INTEGER | SYM_SPEC_REAL, 0) );
      }
   }
   else
   {
      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_INTEGER | SYM_SPEC_REAL, 0) );
   }
   assert( propdata->binvaraffected || ! propdata->ofenabled || ! propdata->symconsenabled );

   if ( propdata->nperms <= 0 || (! propdata->symconsenabled && ! propdata->sstenabled) )
      return SCIP_OKAY;

   if ( ! propdata->binvaraffected )
   {
      SCIPdebugMsg(scip, "Symmetry propagator: problem is linear and no symmetry on binary variables has been found, turning symretope constraints off.\n");
      propdata->symconsenabled = FALSE;
   }
   assert( propdata->nperms > 0 );
   assert( hasNonlinearConstraints(propdata) || propdata->binvaraffected || propdata->sstenabled );

   propdata->triedaddconss = TRUE;

   if ( propdata->symconsenabled )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->genorbconss, propdata->nperms) );

      if ( propdata->detectorbitopes )
      {
         SCIP_CALL( detectOrbitopes(scip, propdata, propdata->components, propdata->componentbegins, propdata->ncomponents) );
      }
   }

   /* disable orbital fixing if all components are handled by orbitopes */
   if ( propdata->ncomponents == propdata->norbitopes )
      propdata->ofenabled = FALSE;

   /* possibly stop */
   if ( SCIPisStopped(scip) )
   {
      if ( propdata->ngenorbconss == 0 )
      {
         SCIPfreeBlockMemoryArrayNull(scip, &propdata->genorbconss, propdata->nperms);
      }
      return SCIP_OKAY;
   }

   if ( propdata->ncompblocked < propdata->ncomponents && propdata->detectsubgroups && propdata->symconsenabled )
   {
      /* @TODO: create array only when needed */
      propdata->genlinconsssize = propdata->nperms;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->genlinconss, propdata->genlinconsssize) );

      SCIP_CALL( detectAndHandleSubgroups(scip, propdata) );
   }

   if ( propdata->sstenabled )
   {
      SCIP_CALL( addSSTConss(scip, propdata, nchgbds) );
   }

   /* possibly stop */
   if ( SCIPisStopped(scip) || ! propdata->symconsenabled )
      return SCIP_OKAY;

   /* add symmetry breaking constraints if orbital fixing is not used outside orbitopes */
   if ( ! propdata->ofenabled )
   {
      /* exit if no or only trivial symmetry group is available */
      if ( propdata->nperms < 1 || ! propdata->binvaraffected )
         return SCIP_OKAY;

      if ( propdata->addsymresacks )
      {
         SCIP_CALL( addSymresackConss(scip, prop, propdata->components, propdata->componentbegins, propdata->ncomponents) );
      }

      /* free symmetry conss if no orbitope/symresack constraints have been found (may happen if Schreier-Sims constraints are active) */
      if ( propdata->ngenorbconss == 0 )
         SCIPfreeBlockMemoryArrayNull(scip, &propdata->genorbconss, propdata->nperms);
   }

   return SCIP_OKAY;
}



/*
 * Local methods for orbital fixing
 */


/** performs orbital fixing
 *
 *  Note that we do not have to distinguish between variables that have been fixed or branched to 1, since the
 *  stabilizer is with respect to the variables that have been branched to 1. Thus, if an orbit contains a variable that
 *  has been branched to 1, the whole orbit only contains variables that have been branched to 1 - and nothing can be
 *  fixed.
 */
static
SCIP_RETCODE performOrbitalFixing(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR**            permvars,           /**< variables */
   int                   npermvars,          /**< number of variables */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int                   norbits,            /**< number of orbits */
   SCIP_Bool*            infeasible,         /**< pointer to store whether problem is infeasible */
   int*                  nfixedzero,         /**< pointer to store number of variables fixed to 0 */
   int*                  nfixedone           /**< pointer to store number of variables fixed to 1 */
   )
{
   SCIP_Bool tightened;
   int i;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( infeasible != NULL );
   assert( nfixedzero != NULL );
   assert( nfixedone != NULL );
   assert( norbits > 0 );
   assert( orbitbegins[0] == 0 );

   *infeasible = FALSE;
   *nfixedzero = 0;
   *nfixedone = 0;

   /* check all orbits */
   for (i = 0; i < norbits; ++i)
   {
      SCIP_Bool havefixedone = FALSE;
      SCIP_Bool havefixedzero = FALSE;
      SCIP_VAR* var;
      int j;

      /* we only have nontrivial orbits */
      assert( orbitbegins[i+1] - orbitbegins[i] >= 2 );

      /* check all variables in the orbit */
      for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
      {
         assert( 0 <= orbits[j] && orbits[j] < npermvars );
         var = permvars[orbits[j]];
         assert( var != NULL );

         /* check whether variable is not binary (and not implicit integer!) */
         if ( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
         {
            /* skip orbit if there are non-binary variables */
            havefixedone = FALSE;
            havefixedzero = FALSE;
            break;
         }

         /* if variable is fixed to 1 -> can fix all variables in orbit to 1 */
         if ( SCIPvarGetLbLocal(var) > 0.5 )
            havefixedone = TRUE;

         /* check for zero-fixed variables */
         if ( SCIPvarGetUbLocal(var) < 0.5 )
            havefixedzero = TRUE;
      }

      /* check consistency */
      if ( havefixedone && havefixedzero )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }

      /* fix all variables to 0 if there is one variable fixed to 0 */
      if ( havefixedzero )
      {
         assert( ! havefixedone );

         for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         {
            assert( 0 <= orbits[j] && orbits[j] < npermvars );
            var = permvars[orbits[j]];
            assert( var != NULL );

            /* only variables that are not yet fixed to 0 */
            if ( SCIPvarGetUbLocal(var) > 0.5 )
            {
               SCIPdebugMsg(scip, "can fix <%s> (index %d) to 0.\n", SCIPvarGetName(var), orbits[j]);
               assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
               /* due to aggregation, var might already be fixed to 1, so do not put assert here */

               /* do not use SCIPinferBinvarProp(), since conflict analysis is not valid */
               SCIP_CALL( SCIPtightenVarUb(scip, var, 0.0, FALSE, infeasible, &tightened) );
               if ( *infeasible )
                  return SCIP_OKAY;
               if ( tightened )
                  ++(*nfixedzero);
            }
         }
      }

      /* fix all variables to 1 if there is one variable fixed to 1 */
      if ( havefixedone )
      {
         assert( ! havefixedzero );

         for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         {
            assert( 0 <= orbits[j] && orbits[j] < npermvars );
            var = permvars[orbits[j]];
            assert( var != NULL );

            /* only variables that are not yet fixed to 1 */
            if ( SCIPvarGetLbLocal(var) < 0.5)
            {
               SCIPdebugMsg(scip, "can fix <%s> (index %d) to 1.\n", SCIPvarGetName(var), orbits[j]);
               assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
               /* due to aggregation, var might already be fixed to 0, so do not put assert here */

               /* do not use SCIPinferBinvarProp(), since conflict analysis is not valid */
               SCIP_CALL( SCIPtightenVarLb(scip, var, 1.0, FALSE, infeasible, &tightened) );
               if ( *infeasible )
                  return SCIP_OKAY;
               if ( tightened )
                  ++(*nfixedone);
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** Gets branching variables on the path to root
 *
 *  The variables are added to bg1 and bg1list, which are prefilled with the variables globally fixed to 1.
 */
static
SCIP_RETCODE computeBranchingVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   nvars,              /**< number of variables */
   SCIP_HASHMAP*         varmap,             /**< map of variables to indices in vars array */
   SCIP_Shortbool*       bg1,                /**< bitset marking the variables globally fixed or branched to 1 */
   int*                  bg1list,            /**< array to store the variable indices globally fixed or branched to 1 */
   int*                  nbg1                /**< pointer to store the number of variables in bg1 and bg1list */
   )
{
   SCIP_NODE* node;

   assert( scip != NULL );
   assert( varmap != NULL );
   assert( bg1 != NULL );
   assert( bg1list != NULL );
   assert( nbg1 != NULL );
   assert( *nbg1 >= 0 );

   /* get current node */
   node = SCIPgetCurrentNode(scip);

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPprintNodeRootPath(scip, node, NULL) );
#endif

   /* follow path to the root (in the root no domains were changed due to branching) */
   while ( SCIPnodeGetDepth(node) != 0 )
   {
      SCIP_BOUNDCHG* boundchg;
      SCIP_DOMCHG* domchg;
      SCIP_VAR* branchvar;
      int nboundchgs;
      int i;

      /* get domain changes of current node */
      domchg = SCIPnodeGetDomchg(node);

      /* If we stopped due to a solving limit, it might happen that a non-root node has no domain changes, in all other
       * cases domchg should not be NULL. */
      if ( domchg != NULL )
      {
         /* loop through all bound changes */
         nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
         for (i = 0; i < nboundchgs; ++i)
         {
            /* get bound change info */
            boundchg = SCIPdomchgGetBoundchg(domchg, i);
            assert( boundchg != NULL );

            /* branching decisions have to be in the beginning of the bound change array */
            if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
               break;

            /* get corresponding branching variable */
            branchvar = SCIPboundchgGetVar(boundchg);

            /* we only consider binary variables */
            if ( SCIPvarGetType(branchvar) == SCIP_VARTYPE_BINARY )
            {
               /* if branching variable is not known (may have been created meanwhile,
                * e.g., by prop_inttobinary; may have been removed from symmetry data
                * due to compression), continue with parent node */
               if ( ! SCIPhashmapExists(varmap, (void*) branchvar) )
                  break;

               if ( SCIPvarGetLbLocal(branchvar) > 0.5 )
               {
                  int branchvaridx;

                  branchvaridx = SCIPhashmapGetImageInt(varmap, (void*) branchvar);
                  assert( branchvaridx < nvars );

                  /* the variable might already be fixed to 1 */
                  if ( ! bg1[branchvaridx] )
                  {
                     bg1[branchvaridx] = TRUE;
                     bg1list[(*nbg1)++] = branchvaridx;
                  }
               }
            }
         }
      }

      node = SCIPnodeGetParent(node);
   }

   return SCIP_OKAY;
}


/** propagates orbital fixing */
static
SCIP_RETCODE propagateOrbitalFixing(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry breaking propagator */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the node is detected to be infeasible */
   int*                  nprop               /**< pointer to store the number of propagations */
   )
{
   SCIP_Shortbool* inactiveperms;
   SCIP_Shortbool* bg0;
   SCIP_Shortbool* bg1;
   SCIP_VAR** permvars;
   int* orbitbegins;
   int* orbits;
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int ncomponents;
   int* bg0list;
   int nbg0;
   int* bg1list;
   int nbg1;
   int nactiveperms;
   int norbits;
   int npermvars;
   int nbinpermvars;
   int** permstrans;
   int nperms;
   int p;
   int v;
   int j;
   int componentidx;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->ofenabled );
   assert( infeasible != NULL );
   assert( nprop != NULL );

   *infeasible = FALSE;
   *nprop = 0;

   /* possibly compute symmetry; fix non-binary potential branching variables */
   if ( hasNonlinearConstraints(propdata) || propdata->symfixnonbinaryvars )
   {
      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY, SYM_SPEC_INTEGER | SYM_SPEC_REAL) );
   }
   else
   {
      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_REAL, SYM_SPEC_INTEGER) );
   }
   assert( hasNonlinearConstraints(propdata) || propdata->binvaraffected || ! propdata->ofenabled );

   /* return if there is no symmetry available */
   nperms = propdata->nperms;
   if ( nperms <= 0 || ! propdata->ofenabled )
      return SCIP_OKAY;

   assert( propdata->permvars != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->permvarmap != NULL );
   assert( propdata->permstrans != NULL );
   assert( propdata->inactiveperms != NULL );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );
   assert( propdata->vartocomponent != NULL );
   assert( propdata->ncomponents > 0 );

   permvars = propdata->permvars;
   npermvars = propdata->npermvars;
   nbinpermvars = propdata->nbinpermvars;
   permstrans = propdata->permstrans;
   inactiveperms = propdata->inactiveperms;
   components = propdata->components;
   componentbegins = propdata->componentbegins;
   vartocomponent = propdata->vartocomponent;
   ncomponents = propdata->ncomponents;

   /* init bitset for marking variables (globally fixed or) branched to 1 */
   assert( propdata->bg1 != NULL );
   assert( propdata->bg1list != NULL );
   assert( propdata->nbg1 >= 0 );
   assert( propdata->nbg1 <= npermvars );

   bg1 = propdata->bg1;
   bg1list = propdata->bg1list;
   nbg1 = propdata->nbg1;

   /* get branching variables */
   SCIP_CALL( computeBranchingVariables(scip, npermvars, propdata->permvarmap, bg1, bg1list, &nbg1) );
   assert( nbg1 >= propdata->nbg1 );

   /* reset inactive permutations */
   nactiveperms = nperms;
   for (p = 0; p < nperms; ++p)
      propdata->inactiveperms[p] = FALSE;

   /* get pointers for bg0 */
   assert( propdata->bg0 != NULL );
   assert( propdata->bg0list != NULL );
   assert( propdata->nbg0 >= 0 );
   assert( propdata->nbg0 <= npermvars );

   bg0 = propdata->bg0;
   bg0list = propdata->bg0list;
   nbg0 = propdata->nbg0;

   /* filter out permutations that move variables that are fixed to 0 */
   for (j = 0; j < nbg0 && nactiveperms > 0; ++j)
   {
      int* pt;

      v = bg0list[j];
      assert( 0 <= v && v < npermvars );
      assert( bg0[v] );

      componentidx = vartocomponent[v];

      /* skip unaffected variables and blocked components */
      if ( componentidx < 0 || propdata->componentblocked[componentidx] )
         continue;

      pt = permstrans[v];
      assert( pt != NULL );

      for (p = componentbegins[componentidx]; p < componentbegins[componentidx + 1]; ++p)
      {
         int img;
         int perm;

         perm = components[p];

         /* skip inactive permutations */
         if ( inactiveperms[perm] )
            continue;

         img = pt[perm];

         if ( img != v )
         {
#ifndef NDEBUG
            SCIP_VAR* varv = permvars[v];
            SCIP_VAR* varimg = permvars[img];

            /* check whether moved variables have the same type (might have been aggregated in the meanwhile) */
            assert( SCIPvarGetType(varv) == SCIPvarGetType(varimg) ||
               (SCIPvarIsBinary(varv) && SCIPvarIsBinary(varimg)) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_IMPLINT && SCIPvarGetType(varimg) == SCIP_VARTYPE_CONTINUOUS &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(varimg) == SCIP_VARTYPE_IMPLINT &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) );
            assert( SCIPisEQ(scip, propdata->permvarsobj[v], propdata->permvarsobj[img]) );
#endif

            /* we are moving a variable globally fixed to 0 to a variable not of this type */
            if ( ! bg0[img] )
            {
               inactiveperms[perm] = TRUE; /* mark as inactive */
               --nactiveperms;
            }
         }
      }
   }

   /* filter out permutations that move variables that are fixed to different values */
   for (j = 0; j < nbg1 && nactiveperms > 0; ++j)
   {
      int* pt;

      v = bg1list[j];
      assert( 0 <= v && v < npermvars );
      assert( bg1[v] );

      componentidx = vartocomponent[v];

      /* skip unaffected variables and blocked components */
      if ( componentidx < 0 || propdata->componentblocked[componentidx] )
         continue;

      pt = permstrans[v];
      assert( pt != NULL );

      for (p = componentbegins[componentidx]; p < componentbegins[componentidx + 1]; ++p)
      {
         int img;
         int perm;

         perm = components[p];

         /* skip inactive permutations */
         if ( inactiveperms[perm] )
            continue;

         img = pt[perm];

         if ( img != v )
         {
#ifndef NDEBUG
            SCIP_VAR* varv = permvars[v];
            SCIP_VAR* varimg = permvars[img];

            /* check whether moved variables have the same type (might have been aggregated in the meanwhile) */
            assert( SCIPvarGetType(varv) == SCIPvarGetType(varimg) ||
               (SCIPvarIsBinary(varv) && SCIPvarIsBinary(varimg)) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_IMPLINT && SCIPvarGetType(varimg) == SCIP_VARTYPE_CONTINUOUS &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(varimg) == SCIP_VARTYPE_IMPLINT &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) );
            assert( SCIPisEQ(scip, propdata->permvarsobj[v], propdata->permvarsobj[img]) );
#endif

            /* we are moving a variable globally fixed or branched to 1 to a variable not of this type */
            if ( ! bg1[img] )
            {
               inactiveperms[perm] = TRUE; /* mark as inactive */
               --nactiveperms;
            }
         }
      }
   }

   /* Clean bg1 list - need to do this after the main loop! (Not needed for bg0.)
    * Note that variables globally fixed to 1 are not resetted, since the loop starts at propdata->nbg1. */
   for (j = propdata->nbg1; j < nbg1; ++j)
      bg1[bg1list[j]] = FALSE;

   /* exit if no active permuations left */
   if ( nactiveperms == 0 )
      return SCIP_OKAY;

   /* compute orbits of binary variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, nbinpermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, nbinpermvars) );
   SCIP_CALL( SCIPcomputeOrbitsFilterSym(scip, nbinpermvars, permstrans, nperms, inactiveperms,
         orbits, orbitbegins, &norbits, components, componentbegins, vartocomponent, propdata->componentblocked, ncomponents, propdata->nmovedpermvars) );

   if ( norbits > 0 )
   {
      int nfixedzero = 0;
      int nfixedone = 0;

      SCIPdebugMsg(scip, "Perform orbital fixing on %d orbits (%d active perms).\n", norbits, nactiveperms);
      SCIP_CALL( performOrbitalFixing(scip, permvars, nbinpermvars, orbits, orbitbegins, norbits, infeasible, &nfixedzero, &nfixedone) );

      propdata->nfixedzero += nfixedzero;
      propdata->nfixedone += nfixedone;
      *nprop = nfixedzero + nfixedone;

      SCIPdebugMsg(scip, "Orbital fixings: %d 0s, %d 1s.\n", nfixedzero, nfixedone);
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   return SCIP_OKAY;
}



/*
 * Callback methods of propagator
 */

/** presolving initialization method of propagator (called when presolving is about to begin) */
static
SCIP_DECL_PROPINITPRE(propInitpreSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* get nonlinear conshdlr for future checks on whether there are nonlinear constraints */
   propdata->conshdlr_nonlinear = SCIPfindConshdlr(scip, "nonlinear");

   /* check whether we should run */
   if ( propdata->usesymmetry < 0 )
   {
      SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &propdata->usesymmetry) );
      SCIP_CALL( setSymmetryMethodEnabledValues(propdata) );
   }
   else if ( SCIPgetNRuns(scip) > propdata->lastrestart && isSymmetryRecomputationRequired(scip, propdata) )
   {
      assert( SCIPgetNRuns(scip) > 1 );

      SCIP_CALL( setSymmetryMethodEnabledValues(propdata) );
   }

   /* add symmetry handling constraints if required  */
   if ( (propdata->symconsenabled || propdata->sstenabled) && propdata->addconsstiming == 0 )
   {
      SCIPdebugMsg(scip, "Try to add symmetry handling constraints before presolving.");

      SCIP_CALL( tryAddSymmetryHandlingConss(scip, prop, NULL, NULL) );
   }
   else if ( propdata->ofenabled && propdata->ofsymcomptiming == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Symmetry computation before presolving:\n");

      /* otherwise compute symmetry if timing requests it; fix non-binary potential branching variables */
      if ( hasNonlinearConstraints(propdata) || propdata->symfixnonbinaryvars )
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY, SYM_SPEC_INTEGER | SYM_SPEC_REAL) );
      }
      else
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_REAL, SYM_SPEC_INTEGER) );
      }
      assert( propdata->binvaraffected || ! propdata->ofenabled );
   }

   return SCIP_OKAY;
}


/** presolving deinitialization method of propagator (called after presolving has been finished) */
static
SCIP_DECL_PROPEXITPRE(propExitpreSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   SCIPdebugMsg(scip, "Exitpre method of propagator <%s> ...\n", PROP_NAME);

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );

   /* guarantee that symmetries are computed (and handled) if the solving process has not been interrupted
    * and even if presolving has been disabled */
   if ( (propdata->symconsenabled || propdata->sstenabled) && SCIPgetStatus(scip) == SCIP_STATUS_UNKNOWN )
   {
      SCIP_CALL( tryAddSymmetryHandlingConss(scip, prop, NULL, NULL) );
   }

   /* if timing requests it, guarantee that symmetries are computed even if presolving is disabled */
   if ( propdata->ofenabled && propdata->ofsymcomptiming <= 1 && SCIPgetStatus(scip) == SCIP_STATUS_UNKNOWN )
   {
      /* fix non-binary potential branching variables */
      if ( hasNonlinearConstraints(propdata) || propdata->symfixnonbinaryvars )
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY, SYM_SPEC_INTEGER | SYM_SPEC_REAL) );
      }
      else
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_REAL, SYM_SPEC_INTEGER) );
      }
      assert( propdata->binvaraffected || ! propdata->ofenabled );
   }

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int i;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );

   /* possibly create symmetry handling constraints */
   if ( propdata->symconsenabled || propdata->sstenabled )
   {
      int noldngenconns;
      int nchanges = 0;
      SCIP_Bool earlyterm = FALSE;

      /* skip presolving if we are not at the end if addconsstiming == 2 */
      assert( 0 <= propdata->addconsstiming && propdata->addconsstiming <= SYM_COMPUTETIMING_AFTERPRESOL );
      if ( propdata->addconsstiming > SYM_COMPUTETIMING_DURINGPRESOL && ! SCIPisPresolveFinished(scip) )
         return SCIP_OKAY;

      /* possibly stop */
      if ( SCIPisStopped(scip) )
         return SCIP_OKAY;

      noldngenconns = propdata->ngenorbconss + propdata->nsstconss + propdata->ngenlinconss;

      SCIP_CALL( tryAddSymmetryHandlingConss(scip, prop, &nchanges, &earlyterm) );

      /* if we actually tried to add symmetry handling constraints */
      if ( ! earlyterm ) /*lint !e774*/
      {
         *result = SCIP_DIDNOTFIND;

         if ( nchanges > 0 )
         {
            *result = SCIP_SUCCESS;
            *nchgbds += nchanges;
         }

         /* if symmetry handling constraints have been added, presolve each */
         if ( propdata->ngenorbconss > 0 || propdata->ngenlinconss > 0 || propdata->nsstconss > 0 )
         {
            /* at this point, the symmetry group should be computed and nontrivial */
            assert( propdata->nperms > 0 );
            assert( propdata->triedaddconss );

            /* we have added at least one symmetry handling constraints, i.e., we were successful */
            *result = SCIP_SUCCESS;

            *naddconss += propdata->ngenorbconss + propdata->ngenlinconss + propdata->nsstconss - noldngenconns;
            SCIPdebugMsg(scip, "Added symmetry breaking constraints: %d.\n", *naddconss);

            /* if constraints have been added, loop through generated constraints and presolve each */
            for (i = 0; i < propdata->ngenorbconss; ++i)
            {
               SCIP_CALL( SCIPpresolCons(scip, propdata->genorbconss[i], nrounds, SCIP_PROPTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
                     nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
                     nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

               /* exit if cutoff or unboundedness has been detected */
               if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
               {
                  SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(propdata->genorbconss[i]));
                  return SCIP_OKAY;
               }
            }

            for (i = 0; i < propdata->ngenlinconss; ++i)
            {
               SCIP_CALL( SCIPpresolCons(scip, propdata->genlinconss[i], nrounds, SCIP_PROPTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
                     nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
                     nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

               /* exit if cutoff or unboundedness has been detected */
               if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
               {
                  SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(propdata->genlinconss[i]));
                  return SCIP_OKAY;
               }
            }
            SCIPdebugMsg(scip, "Presolved %d generated constraints.\n",
               propdata->ngenorbconss + propdata->ngenlinconss);

            for (i = 0; i < propdata->nsstconss; ++i)
            {
               SCIP_CALL( SCIPpresolCons(scip, propdata->sstconss[i], nrounds, SCIP_PROPTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
                     nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
                     nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

               /* exit if cutoff or unboundedness has been detected */
               if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
               {
                  SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(propdata->sstconss[i]));
                  return SCIP_OKAY;
               }
            }
            SCIPdebugMsg(scip, "Presolved %d generated Schreier Sims constraints.\n", propdata->nsstconss);
         }
      }
   }

   /* run OF presolving */
   assert( 0 <= propdata->ofsymcomptiming && propdata->ofsymcomptiming <= SYM_COMPUTETIMING_AFTERPRESOL );
   if ( propdata->ofenabled && propdata->performpresolving && propdata->ofsymcomptiming <= SYM_COMPUTETIMING_DURINGPRESOL )
   {
      SCIP_Bool infeasible;
      int nprop;

      /* if we have not tried to add symmetry handling constraints */
      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;

      SCIPdebugMsg(scip, "Presolving <%s>.\n", PROP_NAME);

      SCIP_CALL( propagateOrbitalFixing(scip, propdata, &infeasible, &nprop) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         propdata->offoundreduction = TRUE;
      }
      else if ( nprop > 0 )
      {
         *result = SCIP_SUCCESS;
         *nfixedvars += nprop;
         propdata->offoundreduction = TRUE;
      }
   }
   else if ( propdata->ofenabled && propdata->ofsymcomptiming == SYM_COMPUTETIMING_DURINGPRESOL )
   {
      /* otherwise compute symmetry early if timing requests it; fix non-binary potential branching variables */
      if ( hasNonlinearConstraints(propdata) || propdata->symfixnonbinaryvars )
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY, SYM_SPEC_INTEGER | SYM_SPEC_REAL) );
      }
      else
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_REAL, SYM_SPEC_INTEGER) );
      }
      assert( propdata->binvaraffected || ! propdata->ofenabled );
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Longint nodenumber;
   int nprop = 0;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if ( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* do nothing if we are in a probing node */
   if ( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* do not run again in repropagation, since the path to the root might have changed */
   if ( SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* if usesymmetry has not been read so far */
   if ( propdata->usesymmetry < 0 )
   {
      SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &propdata->usesymmetry) );
      SCIP_CALL( setSymmetryMethodEnabledValues(propdata) );
   }
   else if ( SCIPgetNRuns(scip) > propdata->lastrestart && isSymmetryRecomputationRequired(scip, propdata) )
   {
      assert( SCIPgetNRuns(scip) > 1 );

      SCIP_CALL( setSymmetryMethodEnabledValues(propdata) );
   }

   /* do not propagate if orbital fixing is not enabled */
   if ( ! propdata->ofenabled )
      return SCIP_OKAY;

   /* return if there is no symmetry available */
   if ( propdata->nperms == 0 )
      return SCIP_OKAY;

   /* return if we already ran in this node */
   nodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
   if ( nodenumber == propdata->nodenumber )
      return SCIP_OKAY;
   propdata->nodenumber = nodenumber;

   /* propagate */
   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "Propagating <%s>.\n", SCIPpropGetName(prop));

   SCIP_CALL( propagateOrbitalFixing(scip, propdata, &infeasible, &nprop) );

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      propdata->offoundreduction = TRUE;
   }
   else if ( nprop > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      propdata->offoundreduction = TRUE;
   }

   return SCIP_OKAY;
}


/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitSymmetry)
{
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   SCIPdebugMsg(scip, "Exiting propagator <%s>.\n", PROP_NAME);

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( freeSymmetryData(scip, propdata) );

   /* reset basic data */
   propdata->usesymmetry = -1;
   propdata->symconsenabled = FALSE;
   propdata->triedaddconss = FALSE;
   propdata->nsymresacks = 0;
   propdata->norbitopes = 0;
   propdata->ofenabled = FALSE;
   propdata->sstenabled = FALSE;
   propdata->lastrestart = 0;
   propdata->nfixedzero = 0;
   propdata->nfixedone = 0;
   propdata->nodenumber = -1;
   propdata->offoundreduction = FALSE;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator
 *
 *  @todo Implement reverse propagation.
 *
 *  Note that this is relatively difficult to obtain: One needs to include all bounds of variables that are responsible
 *  for creating the orbit in which the variables that was propagated lies. This includes all variables that are moved
 *  by the permutations which are involved in creating the orbit.
 */
static
SCIP_DECL_PROPRESPROP(propRespropSymmetry)
{  /*lint --e{715,818}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   SCIPdebugMsg(scip, "Freeing symmetry propagator.\n");

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}


/*
 * External methods
 */

/** include symmetry propagator */
SCIP_RETCODE SCIPincludePropSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLEDATA* tabledata;
   SCIP_PROPDATA* propdata = NULL;
   SCIP_PROP* prop = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   assert( propdata != NULL );

   propdata->npermvars = 0;
   propdata->nbinpermvars = 0;
   propdata->permvars = NULL;
#ifndef NDEBUG
   propdata->permvarsobj = NULL;
#endif
   propdata->nperms = -1;
   propdata->nmaxperms = 0;
   propdata->perms = NULL;
   propdata->permstrans = NULL;
   propdata->permvarmap = NULL;
   propdata->nonbinpermvarcaptured = NULL;

   propdata->ncomponents = -1;
   propdata->ncompblocked = 0;
   propdata->components = NULL;
   propdata->componentbegins = NULL;
   propdata->vartocomponent = NULL;
   propdata->componentblocked = NULL;

   propdata->log10groupsize = -1.0;
   propdata->nmovedvars = -1;
   propdata->binvaraffected = FALSE;
   propdata->computedsymmetry = FALSE;
   propdata->conshdlr_nonlinear = NULL;

   propdata->usesymmetry = -1;
   propdata->symconsenabled = FALSE;
   propdata->triedaddconss = FALSE;
   propdata->genorbconss = NULL;
   propdata->genlinconss = NULL;
   propdata->ngenorbconss = 0;
   propdata->ngenlinconss = 0;
   propdata->genlinconsssize = 0;
   propdata->nsymresacks = 0;
   propdata->norbitopes = 0;
   propdata->isnonlinvar = NULL;

   propdata->ofenabled = FALSE;
   propdata->bg0 = NULL;
   propdata->bg0list = NULL;
   propdata->nbg0 = 0;
   propdata->bg1 = NULL;
   propdata->bg1list = NULL;
   propdata->nbg1 = 0;
   propdata->permvarsevents = NULL;
   propdata->inactiveperms = NULL;
   propdata->nmovedpermvars = -1;
   propdata->nmovedbinpermvars = 0;
   propdata->nmovedintpermvars = 0;
   propdata->nmovedimplintpermvars = 0;
   propdata->nmovedcontpermvars = 0;
   propdata->lastrestart = 0;
   propdata->nfixedzero = 0;
   propdata->nfixedone = 0;
   propdata->nodenumber = -1;
   propdata->offoundreduction = FALSE;

   propdata->sstenabled = FALSE;
   propdata->sstconss = NULL;
   propdata->nsstconss = 0;
   propdata->maxnsstconss = 0;
   propdata->leaders = NULL;
   propdata->nleaders = 0;
   propdata->maxnleaders = 0;

   /* create event handler */
   propdata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(propdata->eventhdlr), EVENTHDLR_SYMMETRY_NAME, EVENTHDLR_SYMMETRY_DESC,
         eventExecSymmetry, NULL) );
   assert( propdata->eventhdlr != NULL );

   /* include constraint handler */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC,
         PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, propExecSymmetry, propdata) );
   assert( prop != NULL );

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSymmetry) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitSymmetry) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreSymmetry) );
   SCIP_CALL( SCIPsetPropExitpre(scip, prop, propExitpreSymmetry) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropSymmetry) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolSymmetry, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );

   /* include table */
   SCIP_CALL( SCIPallocBlockMemory(scip, &tabledata) );
   tabledata->propdata = propdata;
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_ORBITALFIXING, TABLE_DESC_ORBITALFIXING, TRUE,
         NULL, tableFreeOrbitalfixing, NULL, NULL, NULL, NULL, tableOutputOrbitalfixing,
         tabledata, TABLE_POSITION_ORBITALFIXING, TABLE_EARLIEST_ORBITALFIXING) );

   /* add parameters for computing symmetry */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxgenerators",
         "limit on the number of generators that should be produced within symmetry detection (0 = no limit)",
         &propdata->maxgenerators, TRUE, DEFAULT_MAXGENERATORS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/checksymmetries",
         "Should all symmetries be checked after computation?",
         &propdata->checksymmetries, TRUE, DEFAULT_CHECKSYMMETRIES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/displaynorbitvars",
         "Should the number of variables affected by some symmetry be displayed?",
         &propdata->displaynorbitvars, TRUE, DEFAULT_DISPLAYNORBITVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/doubleequations",
         "Double equations to positive/negative version?",
         &propdata->doubleequations, TRUE, DEFAULT_DOUBLEEQUATIONS, NULL, NULL) );

   /* add parameters for adding symmetry handling constraints */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/conssaddlp",
         "Should the symmetry breaking constraints be added to the LP?",
         &propdata->conssaddlp, TRUE, DEFAULT_CONSSADDLP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addsymresacks",
         "Add inequalities for symresacks for each generator?",
         &propdata->addsymresacks, TRUE, DEFAULT_ADDSYMRESACKS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/detectorbitopes",
         "Should we check whether the components of the symmetry group can be handled by orbitopes?",
         &propdata->detectorbitopes, TRUE, DEFAULT_DETECTORBITOPES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/detectsubgroups",
         "Should we try to detect symmetric subgroups of the symmetry group on binary variables?",
         &propdata->detectsubgroups, TRUE, DEFAULT_DETECTSUBGROUPS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addweaksbcs",
         "Should we add weak SBCs for enclosing orbit of symmetric subgroups?",
         &propdata->addweaksbcs, TRUE, DEFAULT_ADDWEAKSBCS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/addconsstiming",
         "timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving)",
         &propdata->addconsstiming, TRUE, DEFAULT_ADDCONSSTIMING, 0, 2, NULL, NULL) );

   /* add parameters for orbital fixing */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/ofsymcomptiming",
         "timing of symmetry computation for orbital fixing (0 = before presolving, 1 = during presolving, 2 = at first call)",
         &propdata->ofsymcomptiming, TRUE, DEFAULT_OFSYMCOMPTIMING, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/performpresolving",
         "run orbital fixing during presolving?",
         &propdata->performpresolving, TRUE, DEFAULT_PERFORMPRESOLVING, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/recomputerestart",
         "recompute symmetries after a restart has occured? (0 = never, 1 = always, 2 = if OF found reduction)",
         &propdata->recomputerestart, TRUE, DEFAULT_RECOMPUTERESTART, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/compresssymmetries",
         "Should non-affected variables be removed from permutation to save memory?",
         &propdata->compresssymmetries, TRUE, DEFAULT_COMPRESSSYMMETRIES, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "propagating/" PROP_NAME "/compressthreshold",
         "Compression is used if percentage of moved vars is at most the threshold.",
         &propdata->compressthreshold, TRUE, DEFAULT_COMPRESSTHRESHOLD, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/usecolumnsparsity",
         "Should the number of conss a variable is contained in be exploited in symmetry detection?",
         &propdata->usecolumnsparsity, TRUE, DEFAULT_USECOLUMNSPARSITY, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxnconsssubgroup",
         "maximum number of constraints up to which subgroup structures are detected",
         &propdata->maxnconsssubgroup, TRUE, DEFAULT_MAXNCONSSSUBGROUP, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/usedynamicprop",
         "whether dynamic propagation should be used for full orbitopes",
         &propdata->usedynamicprop, TRUE, DEFAULT_USEDYNAMICPROP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addstrongsbcs",
         "Should strong SBCs for enclosing orbit of symmetric subgroups be added if orbitopes are not used?",
         &propdata->addstrongsbcs, TRUE, DEFAULT_ADDSTRONGSBCS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/ssttiebreakrule",
         "rule to select the orbit in Schreier Sims inequalities (variable in 0: minimum size orbit; 1: maximum size orbit; 2: orbit with most variables in conflict with leader)",
         &propdata->ssttiebreakrule, TRUE, DEFAULT_SSTTIEBREAKRULE, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/sstleaderrule",
         "rule to select the leader in an orbit (0: first var; 1: last var; 2: var having most conflicting vars in orbit; 3: var having most conflicting vars in problem)",
         &propdata->sstleaderrule, TRUE, DEFAULT_SSTLEADERRULE, 0, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/sstleadervartype",
         "bitset encoding which variable types can be leaders (1: binary; 2: integer; 4: impl. int; 8: continuous);" \
         "if multiple types are allowed, take the one with most affected vars",
         &propdata->sstleadervartype, TRUE, DEFAULT_SSTLEADERVARTYPE, 1, 15, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addconflictcuts",
         "Should Schreier Sims constraints be added if we use a conflict based rule?",
         &propdata->addconflictcuts, TRUE, DEFAULT_ADDCONFLICTCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/sstaddcuts",
         "Should Schreier Sims constraints be added?",
         &propdata->sstaddcuts, TRUE, DEFAULT_SSTADDCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/sstmixedcomponents",
         "Should Schreier Sims constraints be added if a symmetry component contains variables of different types?",
         &propdata->sstmixedcomponents, TRUE, DEFAULT_SSTMIXEDCOMPONENTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/symfixnonbinaryvars",
         "Whether all non-binary variables shall be not affected by symmetries if OF is active?",
         &propdata->symfixnonbinaryvars, TRUE, DEFAULT_SYMFIXNONBINARYVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/onlybinarysymmetry",
         "Is only symmetry on binary variables used?",
         &propdata->onlybinarysymmetry, TRUE, DEFAULT_ONLYBINARYSYMMETRY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/preferlessrows",
         "Shall orbitopes with less rows be preferred in detection?",
         &propdata->preferlessrows, TRUE, DEFAULT_PREFERLESSROWS, NULL, NULL) );

   /* possibly add description */
   if ( SYMcanComputeSymmetry() )
   {
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SYMsymmetryGetName(), SYMsymmetryGetDesc()) );
   }

   return SCIP_OKAY;
}


/** return currently available symmetry group information */
SCIP_RETCODE SCIPgetSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   SCIP_HASHMAP**        permvarmap,         /**< pointer to store hash map of permvars (or NULL) */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix (or NULL)*/
   int***                permstrans,         /**< pointer to store permutation generators as (npermvars x nperms) matrix (or NULL)*/
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of group size (or NULL) */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether binary variables are affected (or NULL) */
   int**                 components,         /**< pointer to store components of symmetry group (or NULL) */
   int**                 componentbegins,    /**< pointer to store begin positions of components in components array (or NULL) */
   int**                 vartocomponent,     /**< pointer to store assignment from variable to its component (or NULL) */
   int*                  ncomponents         /**< pointer to store number of components (or NULL) */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   assert( scip != NULL );
   assert( npermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( perms != NULL || permstrans != NULL );
   assert( ncomponents != NULL || (components == NULL && componentbegins == NULL && vartocomponent == NULL) );

   /* find symmetry propagator */
   prop = SCIPfindProp(scip, "symmetry");
   if ( prop == NULL )
   {
      SCIPerrorMessage("Could not find symmetry propagator.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   *npermvars = propdata->npermvars;
   *permvars = propdata->permvars;

   if ( permvarmap != NULL )
      *permvarmap = propdata->permvarmap;

   *nperms = propdata->nperms;
   if ( perms != NULL )
   {
      *perms = propdata->perms;
      assert( *perms != NULL || *nperms <= 0 );
   }

   if ( permstrans != NULL )
   {
      *permstrans = propdata->permstrans;
      assert( *permstrans != NULL || *nperms <= 0 );
   }

   if ( log10groupsize != NULL )
      *log10groupsize = propdata->log10groupsize;

   if ( binvaraffected != NULL )
      *binvaraffected = propdata->binvaraffected;

   if ( components != NULL )
      *components = propdata->components;

   if ( componentbegins != NULL )
      *componentbegins = propdata->componentbegins;

   if ( vartocomponent )
      *vartocomponent = propdata->vartocomponent;

   if ( ncomponents )
      *ncomponents = propdata->ncomponents;

   return SCIP_OKAY;
}

/** return whether orbital fixing is enabled */
SCIP_Bool SCIPisOrbitalfixingEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );

   prop = SCIPfindProp(scip, PROP_NAME);
   if ( prop == NULL )
      return FALSE;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   return propdata->ofenabled;
}

/** return number of the symmetry group's generators */
int SCIPgetSymmetryNGenerators(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );

   prop = SCIPfindProp(scip, PROP_NAME);
   if ( prop == NULL )
      return 0;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   if ( propdata->nperms < 0 )
      return 0;
   else
      return propdata->nperms;
}
