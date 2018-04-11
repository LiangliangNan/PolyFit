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

/**@file   cons_sos1.c
 * @brief  constraint handler for SOS type 1 constraints
 * @author Tobias Fischer
 * @author Marc Pfetsch
 *
 * A specially ordered set of type 1 (SOS1) is a sequence of variables such that at most one
 * variable is nonzero. The special case of two variables arises, for instance, from equilibrium or
 * complementary conditions like \f$x \cdot y = 0\f$. Note that it is in principle allowed that a
 * variables appears twice, but it then can be fixed to 0.
 *
 * This implementation of this constraint handler is based on classical ideas, see e.g.@n
 *  "Special Facilities in General Mathematical Programming System for
 *  Non-Convex Problems Using Ordered Sets of Variables"@n
 *  E. Beale and J. Tomlin, Proc. 5th IFORS Conference, 447-454 (1970)
 *
 *
 * The order of the variables is determined as follows:
 *
 * - If the constraint is created with SCIPcreateConsSOS1() and weights are given, the weights
 *   determine the order (decreasing weights). Additional variables can be added with
 *   SCIPaddVarSOS1(), which adds a variable with given weight.
 *
 * - If an empty constraint is created and then variables are added with SCIPaddVarSOS1(), weights
 *   are needed and stored.
 *
 * - All other calls ignore the weights, i.e., if a nonempty constraint is created or variables are
 *   added with SCIPappendVarSOS1().
 *
 * The validity of the SOS1 constraints can be enforced by different branching rules:
 *
 * - If classical SOS branching is used, branching is performed on only one SOS1 constraint.
 *   Depending on the parameters, there are two ways to choose this branching constraint. Either
 *   the constraint with the most number of nonzeros or the one with the largest nonzero-variable
 *   weight. The later version allows the user to specify an order for the branching importance of
 *   the constraints. Constraint branching can also be turned off.
 *
 * - Another way is to branch on the neighborhood of a single variable @p i, i.e., in one branch
 *   \f$x_i\f$ is fixed to zero and in the other its neighbors from the conflict graph.
 *
 * - If bipartite branching is used, then we branch using complete bipartite subgraphs of the
 *   conflict graph, i.e., in one branch fix the variables from the first bipartite partition and
 *   the variables from the second bipartite partition in the other.
 *
 * - In addition to variable domain fixings, it is sometimes also possible to add new SOS1
 *   constraints to the branching nodes. This results in a nonstatic conflict graph, which may
 *   change dynamically with every branching node.
 *
 *
 * @todo Possibly allow to generate local cuts via strengthened local cuts (would need to modified coefficients of rows).
 *
 * @todo Check whether we can avoid turning off multi-aggregation (it is sometimes possible to fix a multi-aggregated
 * variable to 0 by fixing the aggregating variables to 0).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_sos1.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/pub_misc.h"
#include "scip/misc.h"
#include "scip/struct_misc.h"
#include "tclique/tclique.h"
#include <string.h>
#include <ctype.h>


/* constraint handler properties */
#define CONSHDLR_NAME          "SOS1"
#define CONSHDLR_DESC          "SOS1 constraint handler"
#define CONSHDLR_SEPAPRIORITY      1000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            10 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_MEDIUM

/* adjacency matrix */
#define DEFAULT_MAXSOSADJACENCY   10000 /**< do not create an adjacency matrix if number of SOS1 variables is larger than predefined value
                                         *   (-1: no limit) */

/* presolving */
#define DEFAULT_MAXEXTENSIONS         1 /**< maximal number of extensions that will be computed for each SOS1 constraint */
#define DEFAULT_MAXTIGHTENBDS         5 /**< maximal number of bound tightening rounds per presolving round (-1: no limit) */
#define DEFAULT_PERFIMPLANALYSIS  FALSE /**< if TRUE then perform implication graph analysis (might add additional SOS1 constraints) */
#define DEFAULT_DEPTHIMPLANALYSIS    -1 /**< number of recursive calls of implication graph analysis (-1: no limit) */

/* propagation */
#define DEFAULT_CONFLICTPROP      TRUE /**< whether to use conflict graph propagation */
#define DEFAULT_IMPLPROP          TRUE /**< whether to use implication graph propagation */
#define DEFAULT_SOSCONSPROP      FALSE /**< whether to use SOS1 constraint propagation */

/* branching rules */
#define DEFAULT_BRANCHSTRATEGIES  "nbs" /**< possible branching strategies (see parameter DEFAULT_BRANCHINGRULE) */
#define DEFAULT_BRANCHINGRULE       'n' /**< which branching rule should be applied ? ('n': neighborhood, 'b': bipartite, 's': SOS1/clique)
                                         *   (note: in some cases an automatic switching to SOS1 branching is possible) */
#define DEFAULT_AUTOSOS1BRANCH     TRUE /**< if TRUE then automatically switch to SOS1 branching if the SOS1 constraints do not overlap */
#define DEFAULT_FIXNONZERO        FALSE /**< if neighborhood branching is used, then fix the branching variable (if positive in sign) to the value of the
                                         *   feasibility tolerance */
#define DEFAULT_ADDCOMPS          FALSE /**< if TRUE then add complementarity constraints to the branching nodes (can be used in combination with
                                         *   neighborhood or bipartite branching) */
#define DEFAULT_MAXADDCOMPS          -1 /**< maximal number of complementarity constraints added per branching node (-1: no limit) */
#define DEFAULT_ADDCOMPSDEPTH        30 /**< only add complementarity constraints to branching nodes for predefined depth (-1: no limit) */
#define DEFAULT_ADDCOMPSFEAS       -0.6 /**< minimal feasibility value for complementarity constraints in order to be added to the branching node */
#define DEFAULT_ADDBDSFEAS          1.0 /**< minimal feasibility value for bound inequalities in order to be added to the branching node */
#define DEFAULT_ADDEXTENDEDBDS     TRUE /**< should added complementarity constraints be extended to SOS1 constraints to get tighter bound inequalities */

/* selection rules */
#define DEFAULT_NSTRONGROUNDS         0 /**< maximal number of strong branching rounds to perform for each node (-1: auto)
                                         *   (only available for neighborhood and bipartite branching) */
#define DEFAULT_NSTRONGITER       10000 /**< maximal number LP iterations to perform for each strong branching round (-2: auto, -1: no limit) */

/* separation */
#define DEFAULT_BOUNDCUTSFROMSOS1 FALSE /**< if TRUE separate bound inequalities from SOS1 constraints */
#define DEFAULT_BOUNDCUTSFROMGRAPH TRUE /**< if TRUE separate bound inequalities from the conflict graph */
#define DEFAULT_AUTOCUTSFROMSOS1   TRUE /**< if TRUE then automatically switch to separating from SOS1 constraints if the SOS1 constraints do not overlap */
#define DEFAULT_BOUNDCUTSFREQ        10 /**< frequency for separating bound cuts; zero means to separate only in the root node */
#define DEFAULT_BOUNDCUTSDEPTH       40 /**< node depth of separating bound cuts (-1: no limit) */
#define DEFAULT_MAXBOUNDCUTS         50 /**< maximal number of bound cuts separated per branching node */
#define DEFAULT_MAXBOUNDCUTSROOT    150 /**< maximal number of bound cuts separated per iteration in the root node */
#define DEFAULT_STRTHENBOUNDCUTS   TRUE /**< if TRUE then bound cuts are strengthened in case bound variables are available */
#define DEFAULT_IMPLCUTSFREQ          0 /**< frequency for separating implied bound cuts; zero means to separate only in the root node */
#define DEFAULT_IMPLCUTSDEPTH        40 /**< node depth of separating implied bound cuts (-1: no limit) */
#define DEFAULT_MAXIMPLCUTS          50 /**< maximal number of implied bound cuts separated per branching node */
#define DEFAULT_MAXIMPLCUTSROOT     150 /**< maximal number of implied bound cuts separated per iteration in the root node */

/* event handler properties */
#define EVENTHDLR_NAME         "SOS1"
#define EVENTHDLR_DESC         "bound change event handler for SOS1 constraints"

/* defines */
#define DIVINGCUTOFFVALUE     1e6


/** constraint data for SOS1 constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables in the constraint */
   int                   maxvars;            /**< maximal number of variables (= size of storage) */
   int                   nfixednonzeros;     /**< number of variables fixed to be nonzero */
   SCIP_Bool             local;              /**< TRUE if constraint is only valid locally */
   SCIP_VAR**            vars;               /**< variables in constraint */
   SCIP_ROW*             rowlb;              /**< row corresponding to lower bounds, or NULL if not yet created */
   SCIP_ROW*             rowub;              /**< row corresponding to upper bounds, or NULL if not yet created */
   SCIP_Real*            weights;            /**< weights determining the order (ascending), or NULL if not used */
};


/** node data of a given node in the conflict graph */
struct SCIP_NodeData
{
   SCIP_VAR*             var;                /**< variable belonging to node */
   SCIP_VAR*             lbboundvar;         /**< bound variable @p z from constraint \f$x \geq \mu \cdot z\f$ (or NULL if not existent) */
   SCIP_VAR*             ubboundvar;         /**< bound variable @p z from constraint \f$x \leq \mu \cdot z\f$ (or NULL if not existent) */
   SCIP_Real             lbboundcoef;        /**< value \f$\mu\f$ from constraint \f$x \geq \mu z \f$ (0.0 if not existent) */
   SCIP_Real             ubboundcoef;        /**< value \f$\mu\f$ from constraint \f$x \leq \mu z \f$ (0.0 if not existent) */
   SCIP_Bool             lbboundcomp;        /**< TRUE if the nodes from the connected component of the conflict graph the given node belongs to
                                              *   all have the same lower bound variable */
   SCIP_Bool             ubboundcomp;        /**< TRUE if the nodes from the connected component of the conflict graph the given node belongs to
                                              *   all have the same lower bound variable */
};
typedef struct SCIP_NodeData SCIP_NODEDATA;


/** successor data of a given nodes successor in the implication graph */
struct SCIP_SuccData
{
   SCIP_Real             lbimpl;             /**< lower bound implication */
   SCIP_Real             ubimpl;             /**< upper bound implication */
};
typedef struct SCIP_SuccData SCIP_SUCCDATA;


/** tclique data for bound cut generation */
struct TCLIQUE_Data
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr;           /**< SOS1 constraint handler */
   SCIP_DIGRAPH*         conflictgraph;      /**< conflict graph */
   SCIP_SOL*             sol;                /**< LP solution to be separated (or NULL) */
   SCIP_Real             scaleval;           /**< factor for scaling weights */
   SCIP_Bool             cutoff;             /**< whether a cutoff occurred */
   int                   ncuts;              /**< number of bound cuts found in this iteration */
   int                   nboundcuts;         /**< number of bound cuts found so far */
   int                   maxboundcuts;       /**< maximal number of clique cuts separated per separation round (-1: no limit) */
   SCIP_Bool             strthenboundcuts;   /**< if TRUE then bound cuts are strengthened in case bound variables are available */
};


/** SOS1 constraint handler data */
struct SCIP_ConshdlrData
{
   /* conflict graph */
   SCIP_DIGRAPH*         conflictgraph;      /**< conflict graph */
   SCIP_DIGRAPH*         localconflicts;     /**< local conflicts */
   SCIP_Bool             isconflocal;        /**< if TRUE then local conflicts are present and conflict graph has to be updated for each node */
   SCIP_HASHMAP*         varhash;            /**< hash map from variable to node in the conflict graph */
   int                   nsos1vars;          /**< number of problem variables that are part of the SOS1 conflict graph */
   /* adjacency matrix */
   int                   maxsosadjacency;    /**< do not create an adjacency matrix if number of SOS1 variables is larger than predefined
                                              *   value (-1: no limit) */
   /* implication graph */
   SCIP_DIGRAPH*         implgraph;          /**< implication graph (@p j is successor of @p i if and only if \f$ x_i\not = 0 \Rightarrow x_j\not = 0\f$) */
   int                   nimplnodes;         /**< number of nodes in the implication graph */
   /* tclique graph */
   TCLIQUE_GRAPH*        tcliquegraph;       /**< tclique graph data structure */
   TCLIQUE_DATA*         tcliquedata;        /**< tclique data */
   /* event handler */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_VAR**            fixnonzerovars;     /**< stack of variables fixed to nonzero marked by event handler */
   int                   maxnfixnonzerovars; /**< size of stack fixnonzerovars */
   int                   nfixnonzerovars;    /**< number of variables fixed to nonzero marked by event handler */
   /* presolving */
   int                   cntextsos1;         /**< counts number of extended SOS1 constraints */
   int                   maxextensions;      /**< maximal number of extensions that will be computed for each SOS1 constraint */
   int                   maxtightenbds;      /**< maximal number of bound tightening rounds per presolving round (-1: no limit) */
   SCIP_Bool             perfimplanalysis;   /**< if TRUE then perform implication graph analysis (might add additional SOS1 constraints) */
   int                   depthimplanalysis;  /**< number of recursive calls of implication graph analysis (-1: no limit) */
   /* propagation */
   SCIP_Bool             conflictprop;       /**< whether to use conflict graph propagation */
   SCIP_Bool             implprop;           /**< whether to use implication graph propagation */
   SCIP_Bool             sosconsprop;        /**< whether to use SOS1 constraint propagation */
   /* branching */
   char                  branchingrule;      /**< which branching rule should be applied ? ('n': neighborhood, 'b': bipartite, 's': SOS1/clique)
                                              *   (note: in some cases an automatic switching to SOS1 branching is possible) */
   SCIP_Bool             autosos1branch;     /**< if TRUE then automatically switch to SOS1 branching if the SOS1 constraints do not overlap */
   SCIP_Bool             fixnonzero;         /**< if neighborhood branching is used, then fix the branching variable (if positive in sign) to the value of the
                                              *   feasibility tolerance */
   SCIP_Bool             addcomps;           /**< if TRUE then add complementarity constraints to the branching nodes additionally to domain fixings
                                              *   (can be used in combination with neighborhood or bipartite branching) */
   int                   maxaddcomps;        /**< maximal number of complementarity cons. and cor. bound ineq. added per branching node (-1: no limit) */
   int                   addcompsdepth;      /**< only add complementarity constraints to branching nodes for predefined depth (-1: no limit) */
   SCIP_Real             addcompsfeas;       /**< minimal feasibility value for complementarity constraints in order to be added to the branching node */
   SCIP_Real             addbdsfeas;         /**< minimal feasibility value for bound inequalities in order to be added to the branching node */
   SCIP_Bool             addextendedbds;     /**< should added complementarity constraints be extended to SOS1 constraints to get tighter bound inequalities */
   SCIP_Bool             branchsos;          /**< Branch on SOS condition in enforcing? This value can only be set to false if all SOS1 variables are binary */
   SCIP_Bool             branchnonzeros;     /**< Branch on SOS cons. with most number of nonzeros? */
   SCIP_Bool             branchweight;       /**< Branch on SOS cons. with highest nonzero-variable weight for branching - needs branchnonzeros to be false */
   SCIP_Bool             switchsos1branch;   /**< whether to switch to SOS1 branching */
   /* selection rules */
   int                   nstrongrounds;      /**< maximal number of strong branching rounds to perform for each node (-1: auto)
                                              *   (only available for neighborhood and bipartite branching) */
   int                   nstrongiter;        /**< maximal number LP iterations to perform for each strong branching round (-2: auto, -1: no limit) */
   /* separation */
   SCIP_Bool             boundcutsfromsos1;  /**< if TRUE separate bound inequalities from SOS1 constraints */
   SCIP_Bool             boundcutsfromgraph; /**< if TRUE separate bound inequalities from the conflict graph */
   SCIP_Bool             autocutsfromsos1;   /**< if TRUE then automatically switch to separating SOS1 constraints if the SOS1 constraints do not overlap */
   SCIP_Bool             switchcutsfromsos1; /**< whether to switch to separate bound inequalities from SOS1 constraints */
   int                   boundcutsfreq;      /**< frequency for separating bound cuts; zero means to separate only in the root node */
   int                   boundcutsdepth;     /**< node depth of separating bound cuts (-1: no limit) */
   int                   maxboundcuts;       /**< maximal number of bound cuts separated per branching node */
   int                   maxboundcutsroot;   /**< maximal number of bound cuts separated per iteration in the root node */
   int                   nboundcuts;         /**< number of bound cuts found so far */
   SCIP_Bool             strthenboundcuts;   /**< if TRUE then bound cuts are strengthened in case bound variables are available */
   int                   implcutsfreq;       /**< frequency for separating implied bound cuts; zero means to separate only in the root node */
   int                   implcutsdepth;      /**< node depth of separating implied bound cuts (-1: no limit) */
   int                   maximplcuts;        /**< maximal number of implied bound cuts separated per branching node */
   int                   maximplcutsroot;    /**< maximal number of implied bound cuts separated per iteration in the root node */
};



/*
 * local methods
 */

/** returns whether two vertices are adjacent in the conflict graph */
static
SCIP_Bool isConnectedSOS1(
   SCIP_Bool**           adjacencymatrix,    /**< adjacency matrix of conflict graph (lower half) (or NULL if an adjacencymatrix is not at hand) */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph (or NULL if an adjacencymatrix is at hand) */
   int                   vertex1,            /**< first vertex */
   int                   vertex2             /**< second vertex */
   )
{
   assert( adjacencymatrix != NULL || conflictgraph != NULL );

   /* we do not allow self-loops */
   if ( vertex1 == vertex2 )
      return FALSE;

   /* for debugging */
   if ( adjacencymatrix == NULL )
   {
      int succvertex;
      int* succ;
      int nsucc1;
      int nsucc2;
      int j;

      nsucc1 = SCIPdigraphGetNSuccessors(conflictgraph, vertex1);
      nsucc2 = SCIPdigraphGetNSuccessors(conflictgraph, vertex2);

      if ( nsucc1 < 1 || nsucc2 < 1 )
         return FALSE;

      if ( nsucc1 > nsucc2 )
      {
         SCIPswapInts(&vertex1, &vertex2);
         SCIPswapInts(&nsucc1, &nsucc2);
      }

      succ = SCIPdigraphGetSuccessors(conflictgraph, vertex1);
      SCIPsortInt(succ, nsucc1);

      for (j = 0; j < nsucc1; ++j)
      {
         succvertex = succ[j];
         if ( succvertex == vertex2 )
            return TRUE;
         else if ( succvertex > vertex2 )
            return FALSE;
      }
   }
   else
   {
      if ( vertex1 < vertex2 )
         return adjacencymatrix[vertex2][vertex1];
      else
         return adjacencymatrix[vertex1][vertex2];
   }

   return FALSE;
}


/** checks whether a variable violates an SOS1 constraint w.r.t. sol together with at least one other variable */
static
SCIP_Bool isViolatedSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph (or NULL if an adjacencymatrix is at hand) */
   int                   node,               /**< node of variable in the conflict graph */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   )
{
   SCIP_Real solval;
   SCIP_VAR* var;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( node >= 0 );

   var = SCIPnodeGetVarSOS1(conflictgraph, node);
   assert( var != NULL );
   solval = SCIPgetSolVal(scip, sol, var);

   /* check whether variable is nonzero w.r.t. sol and the bounds have not been fixed to zero by propagation */
   if ( ! SCIPisFeasZero(scip, solval) && ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) ) )
   {
      int* succ;
      int nsucc;
      int s;

      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, node);
      succ = SCIPdigraphGetSuccessors(conflictgraph, node);

      /* check whether a neighbor variable is nonzero w.r.t. sol */
      for (s = 0; s < nsucc; ++s)
      {
         var = SCIPnodeGetVarSOS1(conflictgraph, succ[s]);
         assert( var != NULL );
         solval = SCIPgetSolVal(scip, sol, var);
         if ( ! SCIPisFeasZero(scip, solval) && ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) ) )
            return TRUE;
      }
   }

   return FALSE;
}


/** returns solution value of imaginary binary big-M variable of a given node from the conflict graph */
static
SCIP_Real nodeGetSolvalBinaryBigMSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   node                /**< node of the conflict graph */
   )
{
   SCIP_Real bound;
   SCIP_VAR* var;
   SCIP_Real val;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( node >= 0 && node < SCIPdigraphGetNNodes(conflictgraph) );

   var = SCIPnodeGetVarSOS1(conflictgraph, node);
   val = SCIPgetSolVal(scip, sol, var);

   if ( SCIPisFeasNegative(scip, val) )
   {
      bound = SCIPvarGetLbLocal(var);
      assert( SCIPisFeasNegative(scip, bound) );

      if ( SCIPisInfinity(scip, -val) )
         return 1.0;
      else if ( SCIPisInfinity(scip, -bound) )
         return 0.0;
      else
         return (val/bound);
   }
   else if ( SCIPisFeasPositive(scip, val) )
   {
      bound = SCIPvarGetUbLocal(var);
      assert( SCIPisFeasPositive(scip, bound) );
      assert( SCIPisFeasPositive(scip, val) );

      if ( SCIPisInfinity(scip, val) )
         return 1.0;
      else if ( SCIPisInfinity(scip, bound) )
         return 0.0;
      else
         return (val/bound);
   }
   else
      return 0.0;
}


/** gets (variable) lower bound value of current LP relaxation solution for a given node from the conflict graph */
static
SCIP_Real nodeGetSolvalVarboundLbSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   node                /**< node of the conflict graph */
   )
{
   SCIP_NODEDATA* nodedata;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( node >= 0 && node < SCIPdigraphGetNNodes(conflictgraph) );

   /* get node data */
   nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, node);
   assert( nodedata != NULL );

   /* if variable is not involved in a variable upper bound constraint */
   if ( nodedata->lbboundvar == NULL || ! nodedata->lbboundcomp )
      return SCIPvarGetLbLocal(nodedata->var);

   return nodedata->lbboundcoef * SCIPgetSolVal(scip, sol, nodedata->lbboundvar);
}


/** gets (variable) upper bound value of current LP relaxation solution for a given node from the conflict graph */
static
SCIP_Real nodeGetSolvalVarboundUbSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   node                /**< node of the conflict graph */
   )
{
   SCIP_NODEDATA* nodedata;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( node >= 0 && node < SCIPdigraphGetNNodes(conflictgraph) );

   /* get node data */
   nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, node);
   assert( nodedata != NULL );

   /* if variable is not involved in a variable upper bound constraint */
   if ( nodedata->ubboundvar == NULL || ! nodedata->ubboundcomp )
      return SCIPvarGetUbLocal(nodedata->var);

   return nodedata->ubboundcoef * SCIPgetSolVal(scip, sol, nodedata->ubboundvar);
}


/** returns whether variable is part of the SOS1 conflict graph */
static
SCIP_Bool varIsSOS1(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert( conshdlrdata != NULL );
   assert( var != NULL );

   if ( conshdlrdata->varhash == NULL || ! SCIPhashmapExists(conshdlrdata->varhash, var) )
      return FALSE;

   return TRUE;
}


/** returns SOS1 index of variable or -1 if variable is not part of the SOS1 conflict graph */
static
int varGetNodeSOS1(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert( conshdlrdata != NULL );
   assert( var != NULL );
   assert( conshdlrdata->varhash != NULL );

   if ( ! SCIPhashmapExists(conshdlrdata->varhash, var) )
      return -1;

   return (int) (size_t) SCIPhashmapGetImage(conshdlrdata->varhash, var);
}


/** fix variable in given node to 0 or add constraint if variable is multi-aggregated
 *
 *  @todo Try to handle multi-aggregated variables as in fixVariableZero() below.
 */
static
SCIP_RETCODE fixVariableZeroNode(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0*/
   SCIP_NODE*            node,               /**< node */
   SCIP_Bool*            infeasible          /**< if fixing is infeasible */
   )
{
   /* if variable cannot be nonzero */
   *infeasible = FALSE;
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* if variable is multi-aggregated */
   if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_CONS* cons;
      SCIP_Real val;

      val = 1.0;

      if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMsg(scip, "creating constraint to force multi-aggregated variable <%s> to 0.\n", SCIPvarGetName(var));
         /* we have to insert a local constraint var = 0 */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &var, &val, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   else
   {
      if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) )
         SCIP_CALL( SCIPchgVarLbNode(scip, node, var, 0.0) );
      if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
         SCIP_CALL( SCIPchgVarUbNode(scip, node, var, 0.0) );
   }

   return SCIP_OKAY;
}


/** try to fix variable to 0
 *
 *  Try to treat fixing by special consideration of multiaggregated variables. For a multi-aggregation
 *  \f[
 *  x = \sum_{i=1}^n \alpha_i x_i + c,
 *  \f]
 *  we can express the fixing \f$x = 0\f$ by fixing all \f$x_i\f$ to 0 if \f$c = 0\f$ and the lower bounds of \f$x_i\f$
 *  are nonnegative if \f$\alpha_i > 0\f$ or the upper bounds are nonpositive if \f$\alpha_i < 0\f$.
 */
static
SCIP_RETCODE fixVariableZero(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0*/
   SCIP_Bool*            infeasible,         /**< if fixing is infeasible */
   SCIP_Bool*            tightened           /**< if fixing was performed */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   assert( infeasible != NULL );
   assert( tightened != NULL );

   *infeasible = FALSE;
   *tightened = FALSE;

   if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Real aggrconst;

      /* if constant is 0 */
      aggrconst = SCIPvarGetMultaggrConstant(var);
      if ( SCIPisZero(scip, aggrconst) )
      {
         SCIP_VAR** aggrvars;
         SCIP_Real* aggrvals;
         SCIP_Bool allnonnegative = TRUE;
         int naggrvars;
         int i;

         SCIP_CALL( SCIPflattenVarAggregationGraph(scip, var) );

         /* check whether all variables are "nonnegative" */
         naggrvars = SCIPvarGetMultaggrNVars(var);
         aggrvars = SCIPvarGetMultaggrVars(var);
         aggrvals = SCIPvarGetMultaggrScalars(var);
         for (i = 0; i < naggrvars; ++i)
         {
            if ( (SCIPisPositive(scip, aggrvals[i]) && SCIPisNegative(scip, SCIPvarGetLbLocal(aggrvars[i]))) ||
                 (SCIPisNegative(scip, aggrvals[i]) && SCIPisPositive(scip, SCIPvarGetUbLocal(aggrvars[i]))) )
            {
               allnonnegative = FALSE;
               break;
            }
         }

         if ( allnonnegative )
         {
            /* all variables are nonnegative -> fix variables */
            for (i = 0; i < naggrvars; ++i)
            {
               SCIP_Bool fixed;
               SCIP_CALL( SCIPfixVar(scip, aggrvars[i], 0.0, infeasible, &fixed) );
               if ( *infeasible )
                  return SCIP_OKAY;
               *tightened = *tightened || fixed;
            }
         }
      }
   }
   else
   {
      SCIP_CALL( SCIPfixVar(scip, var, 0.0, infeasible, tightened) );
   }

   return SCIP_OKAY;
}


/** fix variable in local node to 0, and return whether the operation was feasible
 *
 *  @note We do not add a linear constraint if the variable is multi-aggregated as in
 *  fixVariableZeroNode(), since this would be too time consuming.
 */
static
SCIP_RETCODE inferVariableZero(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0*/
   SCIP_CONS*            cons,               /**< constraint */
   int                   inferinfo,          /**< info for reverse prop. */
   SCIP_Bool*            infeasible,         /**< if fixing is infeasible */
   SCIP_Bool*            tightened,          /**< if fixing was performed */
   SCIP_Bool*            success             /**< whether fixing was successful, i.e., variable is not multi-aggregated */
   )
{
   *infeasible = FALSE;
   *tightened = FALSE;
   *success = FALSE;

   /* if variable cannot be nonzero */
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* directly fix variable if it is not multi-aggregated */
   if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Bool tighten;

      /* fix lower bound */
      SCIP_CALL( SCIPinferVarLbCons(scip, var, 0.0, cons, inferinfo, FALSE, infeasible, &tighten) );
      *tightened = *tightened || tighten;

      /* fix upper bound */
      SCIP_CALL( SCIPinferVarUbCons(scip, var, 0.0, cons, inferinfo, FALSE, infeasible, &tighten) );
      *tightened = *tightened || tighten;

      *success = TRUE;
   }

   return SCIP_OKAY;
}


/** add lock on variable */
static
SCIP_RETCODE lockVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( var != NULL );

   /* rounding down == bad if lb < 0, rounding up == bad if ub > 0 */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)), SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var))) );

   return SCIP_OKAY;
}


/** remove lock on variable */
static
SCIP_RETCODE unlockVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( var != NULL );

   /* rounding down == bad if lb < 0, rounding up == bad if ub > 0 */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)), SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var))) );

   return SCIP_OKAY;
}


/** ensures that the vars and weights array can store at least num entries */
static
SCIP_RETCODE consdataEnsurevarsSizeSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   num,                /**< minimum number of entries to store */
   SCIP_Bool             reserveWeights      /**< whether the weights array is handled */
   )
{
   assert( consdata != NULL );
   assert( consdata->nvars <= consdata->maxvars );

   if ( num > consdata->maxvars )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->maxvars, newsize) );
      if ( reserveWeights )
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->weights, consdata->maxvars, newsize) );
      consdata->maxvars = newsize;
   }
   assert( num <= consdata->maxvars );

   return SCIP_OKAY;
}


/** handle new variable */
static
SCIP_RETCODE handleNewVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Bool             transformed         /**< whether original variable was transformed */
   )
{
   SCIP_DIGRAPH* conflictgraph;
   int node;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( conshdlrdata != NULL );
   assert( var != NULL );

   /* if we are in transformed problem, catch the variable's events */
   if ( transformed )
   {
      assert( conshdlrdata->eventhdlr != NULL );

      /* catch bound change events of variable */
      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
            (SCIP_EVENTDATA*)cons, NULL) ); /*lint !e740*/

      /* if the variable if fixed to nonzero */
      assert( consdata->nfixednonzeros >= 0 );
      if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
         ++consdata->nfixednonzeros;
   }

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockVariableSOS1(scip, cons, var) );

   /* branching on multiaggregated variables does not seem to work well, so avoid it */
   SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) );

   /* add the new coefficient to the upper bound LP row, if necessary */
   if ( consdata->rowub != NULL && ! SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) && ! SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rowub, var, 1.0/SCIPvarGetUbGlobal(var)) );
   }

   /* add the new coefficient to the lower bound LP row, if necessary */
   if ( consdata->rowlb != NULL && ! SCIPisInfinity(scip, SCIPvarGetLbGlobal(var)) && ! SCIPisZero(scip, SCIPvarGetLbGlobal(var)) )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rowlb, var, 1.0/SCIPvarGetLbGlobal(var)) );
   }

   /* return if the conflict graph has not been created yet */
   conflictgraph = conshdlrdata->conflictgraph;
   if ( conflictgraph == NULL )
      return SCIP_OKAY;

   /* get node of variable in the conflict graph (or -1) */
   node = varGetNodeSOS1(conshdlrdata, var);
   assert( node < conshdlrdata->nsos1vars );

   /* if the variable is not already a node of the conflict graph */
   if ( node < 0 )
   {
      /* variable does not appear in the conflict graph: switch to SOS1 branching rule, which does not make use of a conflict graph
       * @todo: maybe recompute the conflict graph, implication graph and varhash instead */
      SCIPdebugMsg(scip, "Switched to SOS1 branching rule, since conflict graph could be infeasible.\n");
      conshdlrdata->switchsos1branch = TRUE;
      return SCIP_OKAY;
   }

   /* if the constraint is local, then there is no need to act, since local constraints are handled by the local conflict graph in the
    * function enforceConflictgraph() */
   if ( ! consdata->local )
   {
      SCIP_VAR** vars;
      int nvars;
      int v;

      vars = consdata->vars;
      nvars = consdata->nvars;

      for (v = 0; v < nvars; ++v)
      {
         int nodev;

         if ( var == vars[v] )
            continue;

         /* get node of variable in the conflict graph (or -1) */
         nodev = varGetNodeSOS1(conshdlrdata, vars[v]);
         assert( nodev < conshdlrdata->nsos1vars );

         /* if the variable is already a node of the conflict graph */
         if ( nodev >= 0 )
         {
            int nsucc;
            int nsuccv;

            nsucc = SCIPdigraphGetNSuccessors(conflictgraph, node);
            nsuccv = SCIPdigraphGetNSuccessors(conflictgraph, nodev);

            /* add arcs if not existent */
            SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraph, nodev, node, NULL) );
            SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraph, node, nodev, NULL) );

            /* in case of new arcs: sort successors in ascending order */
            if ( nsucc < SCIPdigraphGetNSuccessors(conflictgraph, node) )
            {
               SCIPdebugMsg(scip, "Added new conflict graph arc from variable %s to variable %s.\n", SCIPvarGetName(var), SCIPvarGetName(vars[v]));
               SCIPsortInt(SCIPdigraphGetSuccessors(conflictgraph, node), SCIPdigraphGetNSuccessors(conflictgraph, node));
            }

            if ( nsuccv < SCIPdigraphGetNSuccessors(conflictgraph, nodev) )
            {
               SCIPdebugMsg(scip, "Added new conflict graph arc from variable %s to variable %s.\n", SCIPvarGetName(vars[v]), SCIPvarGetName(var));
               SCIPsortInt(SCIPdigraphGetSuccessors(conflictgraph, nodev), SCIPdigraphGetNSuccessors(conflictgraph, nodev));
            }
         }
         else
         {
            /* variable does not appear in the conflict graph: switch to SOS1 branching rule, which does not make use of a conflict graph
             * @todo: maybe recompute the conflict graph, implication graph and varhash instead */
            SCIPdebugMsg(scip, "Switched to SOS1 branching rule, since conflict graph could be infeasible.\n");
            conshdlrdata->switchsos1branch = TRUE;
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** adds a variable to an SOS1 constraint, at position given by weight - ascending order */
static
SCIP_RETCODE addVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_Real             weight              /**< weight to determine position */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;
   int pos;
   int j;

   assert( var != NULL );
   assert( cons != NULL );
   assert( conshdlrdata != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( consdata->weights == NULL && consdata->maxvars > 0 )
   {
      SCIPerrorMessage("cannot add variable to SOS1 constraint <%s> that does not contain weights.\n", SCIPconsGetName(cons));
      return SCIP_INVALIDCALL;
   }

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if ( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert( var != NULL );
   assert( transformed == SCIPvarIsTransformed(var) );

   SCIP_CALL( consdataEnsurevarsSizeSOS1(scip, consdata, consdata->nvars + 1, TRUE) );
   assert( consdata->weights != NULL );
   assert( consdata->maxvars >= consdata->nvars+1 );

   /* find variable position */
   for (pos = 0; pos < consdata->nvars; ++pos)
   {
      if ( consdata->weights[pos] > weight )
         break;
   }
   assert( 0 <= pos && pos <= consdata->nvars );

   /* move other variables, if necessary */
   for (j = consdata->nvars; j > pos; --j)
   {
      consdata->vars[j] = consdata->vars[j-1];
      consdata->weights[j] = consdata->weights[j-1];
   }

   /* insert variable */
   consdata->vars[pos] = var;
   consdata->weights[pos] = weight;
   ++consdata->nvars;

   /* handle the new variable */
   SCIP_CALL( handleNewVariableSOS1(scip, cons, consdata, conshdlrdata, var, transformed) );

   return SCIP_OKAY;
}


/** appends a variable to an SOS1 constraint */
static
SCIP_RETCODE appendVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert( var != NULL );
   assert( cons != NULL );
   assert( conshdlrdata != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if ( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert( var != NULL );
   assert( transformed == SCIPvarIsTransformed(var) );

   SCIP_CALL( consdataEnsurevarsSizeSOS1(scip, consdata, consdata->nvars + 1, FALSE) );

   /* insert variable */
   consdata->vars[consdata->nvars] = var;
   assert( consdata->weights != NULL || consdata->nvars > 0 );
   if ( consdata->weights != NULL && consdata->nvars > 0 )
      consdata->weights[consdata->nvars] = consdata->weights[consdata->nvars-1] + 1.0;
   ++consdata->nvars;

   /* handle the new variable */
   SCIP_CALL( handleNewVariableSOS1(scip, cons, consdata, conshdlrdata, var, transformed) );

   return SCIP_OKAY;
}


/** deletes a variable of an SOS1 constraint */
static
SCIP_RETCODE deleteVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< corresponding event handler */
   int                   pos                 /**< position of variable in array */
   )
{
   int j;

   assert( 0 <= pos && pos < consdata->nvars );

   /* remove lock of variable */
   SCIP_CALL( unlockVariableSOS1(scip, cons, consdata->vars[pos]) );

   /* drop events on variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)cons, -1) ); /*lint !e740*/

   /* delete variable - need to copy since order is important */
   for (j = pos; j < consdata->nvars-1; ++j)
   {
      consdata->vars[j] = consdata->vars[j+1]; /*lint !e679*/
      if ( consdata->weights != NULL )
         consdata->weights[j] = consdata->weights[j+1]; /*lint !e679*/
   }
   --consdata->nvars;

   return SCIP_OKAY;
}


/* ----------------------------- presolving --------------------------------------*/

/** extends a given clique of the conflict graph
 *
 *  Implementation of the Bron-Kerbosch Algorithm from the paper:
 *  Algorithm 457: Finding all Cliques of an Undirected Graph, Bron & Kerbosch, Commun. ACM, 1973
 */
static
SCIP_RETCODE extensionOperatorSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_Bool**           adjacencymatrix,    /**< adjacencymatrix of the conflict graph (only lower half filled) */
   SCIP_DIGRAPH*         vertexcliquegraph,  /**< graph that contains the information which cliques contain a given vertex
                                              *   vertices of variables = 0, ..., nsos1vars-1; vertices of cliques = nsos1vars, ..., nsos1vars+ncliques-1*/
   int                   nsos1vars,          /**< number of SOS1 variables */
   int                   nconss,             /**< number of SOS1 constraints */
   SCIP_CONS*            cons,               /**< constraint to be extended */
   SCIP_VAR**            vars,               /**< variables of extended clique */
   SCIP_Real*            weights,            /**< weights of extended clique */
   SCIP_Bool             firstcall,          /**< whether this is the first call of extension operator */
   SCIP_Bool             usebacktrack,       /**< whether backtracking is needed for the computation */
   int**                 cliques,            /**< all cliques found so far */
   int*                  ncliques,           /**< number of clique found so far */
   int*                  cliquesizes,        /**< number of variables of current clique */
   int*                  newclique,          /**< clique we want to extended*/
   int*                  workingset,         /**< set of vertices that already served as extension and set of candidates that probably will lead to an extension */
   int                   nworkingset,        /**< length of array workingset */
   int                   nexts,              /**< number of vertices that already served as extension */
   int                   pos,                /**< position of potential candidate */
   int*                  maxextensions,      /**< maximal number of extensions */
   int*                  naddconss,          /**< number of added constraints */
   SCIP_Bool*            success             /**< pointer to store if at least one new clique was found */
   )
{
   int* workingsetnew = NULL;
   int nextsnew;
   int nworkingsetnew;
   int mincands;
   int btriter = 0; /* backtrack iterator */
   int selvertex;
   int selpos = -1;
   int fixvertex = -1;
   int i;
   int j;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( adjacencymatrix != NULL );
   assert( vertexcliquegraph != NULL );
   assert( cons != NULL );
   assert( cliques != NULL );
   assert( cliquesizes != NULL );
   assert( newclique != NULL );
   assert( workingset != NULL );
   assert( maxextensions != NULL );
   assert( naddconss != NULL );
   assert( success != NULL );

   if ( firstcall )
      *success = FALSE;

   mincands = nworkingset;
   if ( mincands < 1 )
      return SCIP_OKAY;

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &workingsetnew, nworkingset) );

#ifdef SCIP_DEBUG
   for (i = 0; i < nexts; ++i)
   {
      int vertex = workingset[i];
      for (j = nexts; j < nworkingset; ++j)
      {
         assert( isConnectedSOS1(adjacencymatrix, NULL, vertex, workingset[j]) );
      }
   }
#endif

   /* determine candidate with minimum number of disconnections */
   for (i = 0; i < nworkingset; ++i)
   {
      int vertex;
      int cnt = 0;

      vertex = workingset[i];

      /* count disconnections */
      for (j = nexts; j < nworkingset && cnt < mincands; ++j)
      {
         if ( vertex != workingset[j] && ! isConnectedSOS1(adjacencymatrix, NULL, vertex, workingset[j]) )
         {
            cnt++;

            /* save position of potential candidate */
            pos = j;
         }
      }

      /* check whether a new minimum was found */
      if ( cnt < mincands )
      {
         fixvertex = vertex;
         mincands = cnt;
         if ( i < nexts )
         {
            assert( pos >= 0 );
            selpos = pos;
         }
         else
         {
            selpos = i;

            /* preincrement */
            btriter = 1;
         }
      }
   }

   /* If fixed point is initially chosen from candidates then number of disconnections will be preincreased by one. */

   /* backtrackcycle */
   for (btriter = mincands + btriter; btriter >= 1; --btriter)
   {
      assert( selpos >= 0);
      assert( fixvertex >= 0);

      /* interchange */
      selvertex = workingset[selpos];
      workingset[selpos] = workingset[nexts];
      workingset[nexts] = selvertex;

      /* create new workingset */
      nextsnew = 0;
      for (j = 0 ; j < nexts; ++j)
      {
         if ( isConnectedSOS1(adjacencymatrix, NULL, selvertex, workingset[j]) )
            workingsetnew[nextsnew++] = workingset[j];
      }
      nworkingsetnew = nextsnew;
      for (j = nexts + 1; j < nworkingset; ++j)
      {
         if ( isConnectedSOS1(adjacencymatrix, NULL, selvertex, workingset[j]) )
            workingsetnew[nworkingsetnew++] = workingset[j];
      }

      newclique[cliquesizes[*ncliques]++] = selvertex;

      /* if we found a new clique */
      if ( nworkingsetnew == 0 )
      {
         char consname[SCIP_MAXSTRLEN];
         SCIP_CONSDATA* consdata;
         SCIP_CONS* newcons;
         int cliqueind;

         cliqueind = nsos1vars + *ncliques; /* index of clique in the vertex-clique graph */

         /* save new clique */
         assert( cliquesizes[*ncliques] >= 0 && cliquesizes[*ncliques] <= nsos1vars );
         assert( *ncliques < MAX(1, conshdlrdata->maxextensions) * nconss );
         SCIP_CALL( SCIPallocBufferArray(scip, &(cliques[*ncliques]), cliquesizes[*ncliques]) );/*lint !e866*/
         for (j = 0 ; j < cliquesizes[*ncliques]; ++j)
         {
            vars[j] = SCIPnodeGetVarSOS1(conshdlrdata->conflictgraph, newclique[j]);
            weights[j] = j+1;
            cliques[*ncliques][j] = newclique[j];
         }

         SCIPsortInt(cliques[*ncliques], cliquesizes[*ncliques]);

         /* create new constraint */
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "extsos1_%" SCIP_LONGINT_FORMAT, conshdlrdata->cntextsos1, conshdlrdata->cntextsos1);

         SCIP_CALL( SCIPcreateConsSOS1(scip, &newcons, consname, cliquesizes[*ncliques], vars, weights,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
               SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
               SCIPconsIsDynamic(cons),
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

         consdata = SCIPconsGetData(newcons);

         /* add directed edges to the vertex-clique graph */
         for (j = 0; j < consdata->nvars; ++j)
         {
            /* add arc from clique vertex to clique (needed in presolRoundConssSOS1() to delete redundand cliques) */
            SCIP_CALL( SCIPdigraphAddArcSafe(vertexcliquegraph, cliques[*ncliques][j], cliqueind, NULL) );
         }

         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

         ++(*naddconss);
         ++(conshdlrdata->cntextsos1);
         ++(*ncliques);
         cliquesizes[*ncliques] = cliquesizes[*ncliques-1]; /* cliquesizes[*ncliques] = size of newclique */

         *success = TRUE;

         --(*maxextensions);

         if ( *maxextensions <= 0 )
         {
            SCIPfreeBufferArray(scip, &workingsetnew);
            return SCIP_OKAY;
         }
      }
      else if ( nextsnew < nworkingsetnew ) /* else if the number of of candidates equals zero */
      {
         /* if backtracking is used, it is necessary to keep the memory for 'workingsetnew' */
         if ( usebacktrack )
         {
            SCIP_CALL( extensionOperatorSOS1(scip, conshdlrdata, adjacencymatrix, vertexcliquegraph, nsos1vars, nconss, cons, vars, weights, FALSE, usebacktrack,
                  cliques, ncliques, cliquesizes, newclique, workingsetnew, nworkingsetnew, nextsnew, pos, maxextensions, naddconss, success) );
            if ( *maxextensions <= 0 )
            {
               SCIPfreeBufferArrayNull(scip, &workingsetnew);
               return SCIP_OKAY;
            }
         }
         else
         {
            int w;

            assert( nworkingset >= nworkingsetnew );
            for (w = 0; w < nworkingsetnew; ++w)
               workingset[w] = workingsetnew[w];
            nworkingset = nworkingsetnew;

            SCIPfreeBufferArrayNull(scip, &workingsetnew);

            SCIP_CALL( extensionOperatorSOS1(scip, conshdlrdata, adjacencymatrix, vertexcliquegraph, nsos1vars, nconss, cons, vars, weights, FALSE, usebacktrack,
                  cliques, ncliques, cliquesizes, newclique, workingset, nworkingset, nextsnew, pos, maxextensions, naddconss, success) );
            assert( *maxextensions <= 0 );
            return SCIP_OKAY;
         }
      }
      assert( workingsetnew != NULL );
      assert( workingset != NULL );

      /* remove selvertex from clique */
      --cliquesizes[*ncliques];

      /* add selvertex to the set of vertices that already served as extension */
      ++nexts;

      if ( btriter > 1 )
      {
         /* select a candidate that is not connected to the fixed vertex */
         for (j = nexts; j < nworkingset; ++j)
         {
            assert( fixvertex != workingset[j] );
            if ( ! isConnectedSOS1(adjacencymatrix, NULL, fixvertex, workingset[j]) )
            {
               selpos = j;
               break;
            }
         }
      }
   }

   SCIPfreeBufferArrayNull(scip, &workingsetnew);

   return SCIP_OKAY;
}


/** generates conflict graph that is induced by the variables of a linear constraint */
static
SCIP_RETCODE genConflictgraphLinearCons(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraphlin,   /**< conflict graph of linear constraint (nodes: 1, ..., nlinvars) */
   SCIP_DIGRAPH*         conflictgraphorig,  /**< original conflict graph (nodes: 1, ..., nsos1vars) */
   SCIP_VAR**            linvars,            /**< linear variables in linear constraint */
   int                   nlinvars,           /**< number of linear variables in linear constraint */
   int*                  posinlinvars        /**< posinlinvars[i] = position (index) of SOS1 variable i in linear constraint,
                                              *   posinlinvars[i]= -1 if @p i is not a SOS1 variable or not a variable of the linear constraint */
   )
{
   int indexinsosvars;
   int indexinlinvars;
   int* succ;
   int nsucc;
   int v;
   int s;

   assert( conflictgraphlin != NULL );
   assert( conflictgraphorig != NULL );
   assert( linvars != NULL );
   assert( posinlinvars != NULL );

   for (v = 1; v < nlinvars; ++v) /* we start with v = 1, since "indexinlinvars < v" (see below) is never fulfilled for v = 0 */
   {
      indexinsosvars = varGetNodeSOS1(conshdlrdata, linvars[v]);

      /* if linvars[v] is contained in at least one SOS1 constraint */
      if ( indexinsosvars >= 0 )
      {
         succ = SCIPdigraphGetSuccessors(conflictgraphorig, indexinsosvars);
         nsucc = SCIPdigraphGetNSuccessors(conflictgraphorig, indexinsosvars);

         for (s = 0; s < nsucc; ++s)
         {
            assert( succ[s] >= 0 );
            indexinlinvars = posinlinvars[succ[s]];
            assert( indexinlinvars < nlinvars );

            if ( indexinlinvars >= 0 && indexinlinvars < v )
            {
               SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraphlin, v, indexinlinvars, NULL) );
               SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraphlin, indexinlinvars, v, NULL) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** determine the common successors of the vertices from the considered clique */
static
SCIP_RETCODE cliqueGetCommonSuccessorsSOS1(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int*                  clique,             /**< current clique */
   SCIP_VAR**            vars,               /**< clique variables */
   int                   nvars,              /**< number of clique variables */
   int*                  comsucc,            /**< pointer to store common successors of clique vertices (size = nvars) */
   int*                  ncomsucc            /**< pointer to store number common successors of clique vertices */
   )
{
   int nsucc;
   int* succ;
   int ind;
   int k = 0;
   int v;
   int i;
   int j;

   assert( conflictgraph != NULL );
   assert( clique != NULL );
   assert( vars != NULL );
   assert( comsucc != NULL );
   assert( ncomsucc != NULL );

   *ncomsucc = 0;

   /* determine the common successors of the vertices from the considered clique */

   /* determine successors of variable var[0] that are not in the clique */
   assert(vars[0] != NULL );
   ind =  varGetNodeSOS1(conshdlrdata, vars[0]);

   if( ind == -1 )
      return SCIP_INVALIDDATA;

   assert( ind < SCIPdigraphGetNNodes(conflictgraph) );
   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, ind);
   succ = SCIPdigraphGetSuccessors(conflictgraph, ind);

   for (j = 0; j < nvars; ++j)
   {
      for (i = k; i < nsucc; ++i)
      {
         if ( succ[i] > clique[j] )
         {
            k = i;
            break;
         }
         else if ( succ[i] == clique[j] )
         {
            k = i + 1;
            break;
         }
         else
            comsucc[(*ncomsucc)++] = succ[i];
      }
   }

   /* for all variables except the first one */
   for (v = 1; v < nvars; ++v)
   {
      int ncomsuccsave = 0;
      k = 0;

      assert(vars[v] != NULL );
      ind =  varGetNodeSOS1(conshdlrdata, vars[v]);
      assert( ind >= 0 && ind < SCIPdigraphGetNNodes(conflictgraph) );

      if ( ind >= 0 )
      {
         nsucc = SCIPdigraphGetNSuccessors(conflictgraph, ind);
         succ = SCIPdigraphGetSuccessors(conflictgraph, ind);

         /* determine successors that are in comsucc */
         for (j = 0; j < *ncomsucc; ++j)
         {
            for (i = k; i < nsucc; ++i)
            {
               if ( succ[i] > comsucc[j] )
               {
                  k = i;
                  break;
               }
               else if ( succ[i] == comsucc[j] )
               {
                  comsucc[ncomsuccsave++] = succ[i];
                  k = i + 1;
                  break;
               }
            }
         }
         *ncomsucc = ncomsuccsave;
      }
   }

   return SCIP_OKAY;
}


/** get nodes whose corresponding SOS1 variables are nonzero if an SOS1 variable of a given node is nonzero */
static
SCIP_RETCODE getSOS1Implications(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR**            vars,               /**< problem and SOS1 variables */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph (@p j is successor of @p i if and only if \f$ x_i\not = 0 \Rightarrow x_j\not = 0\f$) */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_Bool*            implnodes,          /**< implnodes[i] = TRUE if the SOS1 variable corresponding to node i in the implication graph is implied to be nonzero */
   int                   node                /**< node of the implication graph */
   )
{
   SCIP_SUCCDATA** succdatas;
   int sos1node;
   int* succ;
   int nsucc;
   int s;

   assert( scip != NULL );
   assert( implgraph != NULL );
   assert( implnodes != NULL );
   assert( node >= 0 );
   assert( vars[node] != NULL );
   assert( (int) (size_t) SCIPhashmapGetImage(implhash, vars[node]) == node );

   /* get node of variable in the conflict graph (-1 if variable is no SOS1 variable) */
   sos1node = varGetNodeSOS1(conshdlrdata, vars[node]);
   if ( sos1node < 0 )
      return SCIP_OKAY;

   succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, node);
   nsucc = SCIPdigraphGetNSuccessors(implgraph, node);
   succ = SCIPdigraphGetSuccessors(implgraph, node);

   for (s = 0; s < nsucc; ++s)
   {
      SCIP_SUCCDATA* data;
      int succnode;
      succnode = succ[s];
      data = succdatas[s];
      sos1node = varGetNodeSOS1(conshdlrdata, vars[succnode]);

      /* if node is SOS1 and the corresponding variable is implied to be nonzero */
      assert( succdatas[s] != NULL );
      if ( sos1node >= 0 && ! implnodes[sos1node] && ( SCIPisFeasPositive(scip, data->lbimpl) || SCIPisFeasNegative(scip, data->ubimpl) ) )
      {
         assert( sos1node == succnode );
         implnodes[sos1node] = TRUE;
         SCIP_CALL( getSOS1Implications(scip, conshdlrdata, vars, implgraph, implhash, implnodes, succnode) );
      }
   }

   return SCIP_OKAY;
}


/** perform one presolving round for a single SOS1 constraint
 *
 *  We perform the following presolving steps.
 *
 *  - If the bounds of some variable force it to be nonzero, we can
 *    fix all other variables to zero and remove the SOS1 constraints
 *    that contain it.
 *  - If a variable is fixed to zero, we can remove the variable.
 *  - If a variable appears twice, it can be fixed to 0.
 *  - We substitute appregated variables.
 */
static
SCIP_RETCODE presolRoundConsSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_Bool*            substituted,        /**< whether a variable was substituted */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   SCIP_Bool*            success,            /**< whether we performed a successful reduction */
   int*                  ndelconss,          /**< number of deleted constraints */
   int*                  nupgdconss,         /**< number of upgraded constraints */
   int*                  nfixedvars,         /**< number of fixed variables */
   int*                  nremovedvars        /**< number of variables removed */
   )
{
   SCIP_VAR** vars;
   SCIP_Bool allvarsbinary;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   int nfixednonzeros;
   int lastFixedNonzero;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( eventhdlr != NULL );
   assert( cutoff != NULL );
   assert( success != NULL );
   assert( ndelconss != NULL );
   assert( nfixedvars != NULL );
   assert( nremovedvars != NULL );

   *substituted = FALSE;
   *cutoff = FALSE;
   *success = FALSE;

   SCIPdebugMsg(scip, "Presolving SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

   j = 0;
   nfixednonzeros = 0;
   lastFixedNonzero = -1;
   allvarsbinary = TRUE;
   vars = consdata->vars;

   /* check for variables fixed to 0 and bounds that fix a variable to be nonzero */
   while ( j < consdata->nvars )
   {
      int l;
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real scalar;
      SCIP_Real constant;

      scalar = 1.0;
      constant = 0.0;

      /* check for aggregation: if the constant is zero the variable is zero iff the aggregated
       * variable is 0 */
      var = vars[j];
      SCIP_CALL( SCIPgetProbvarSum(scip, &var, &scalar, &constant) );

      /* if constant is zero and we get a different variable, substitute variable */
      if ( SCIPisZero(scip, constant) && ! SCIPisZero(scip, scalar) && var != vars[j] )
      {
         SCIPdebugMsg(scip, "substituted variable <%s> by <%s>.\n", SCIPvarGetName(vars[j]), SCIPvarGetName(var));
         SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)cons, -1) ); /*lint !e740*/
         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)cons, NULL) ); /*lint !e740*/

         /* change the rounding locks */
         SCIP_CALL( unlockVariableSOS1(scip, cons, consdata->vars[j]) );
         SCIP_CALL( lockVariableSOS1(scip, cons, var) );

         vars[j] = var;
         *substituted = TRUE;
      }

      /* check whether the variable appears again later */
      for (l = j+1; l < consdata->nvars; ++l)
      {
         /* if variable appeared before, we can fix it to 0 and remove it */
         if ( vars[j] == vars[l] )
         {
            SCIPdebugMsg(scip, "variable <%s> appears twice in constraint, fixing it to 0.\n", SCIPvarGetName(vars[j]));
            SCIP_CALL( SCIPfixVar(scip, vars[j], 0.0, &infeasible, &fixed) );

            if ( infeasible )
            {
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if ( fixed )
               ++(*nfixedvars);
         }
      }

      /* get bounds */
      lb = SCIPvarGetLbLocal(vars[j]);
      ub = SCIPvarGetUbLocal(vars[j]);

      /* if the variable if fixed to nonzero */
      if ( SCIPisFeasPositive(scip, lb) || SCIPisFeasNegative(scip, ub) )
      {
         ++nfixednonzeros;
         lastFixedNonzero = j;
      }

      /* if the variable is fixed to 0 */
      if ( SCIPisFeasZero(scip, lb) && SCIPisFeasZero(scip, ub) )
      {
         SCIPdebugMsg(scip, "deleting variable <%s> fixed to 0.\n", SCIPvarGetName(vars[j]));
         SCIP_CALL( deleteVarSOS1(scip, cons, consdata, eventhdlr, j) );
         ++(*nremovedvars);
      }
      else
      {
         /* check whether all variables are binary */
         if ( ! SCIPvarIsBinary(vars[j]) )
            allvarsbinary = FALSE;

         ++j;
      }
   }

   /* if the number of variables is less than 2 */
   if ( consdata->nvars < 2 )
   {
      SCIPdebugMsg(scip, "Deleting SOS1 constraint <%s> with < 2 variables.\n", SCIPconsGetName(cons));

      /* delete constraint */
      assert( ! SCIPconsIsModifiable(cons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* if more than one variable are fixed to be nonzero, we are infeasible */
   if ( nfixednonzeros > 1 )
   {
      SCIPdebugMsg(scip, "The problem is infeasible: more than one variable has bounds that keep it from being 0.\n");
      assert( lastFixedNonzero >= 0 );
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if there is exactly one fixed nonzero variable */
   if ( nfixednonzeros == 1 )
   {
      assert( lastFixedNonzero >= 0 );

      /* fix all other variables to zero */
      for (j = 0; j < consdata->nvars; ++j)
      {
         if ( j != lastFixedNonzero )
         {
            SCIP_CALL( fixVariableZero(scip, vars[j], &infeasible, &fixed) );
            if ( infeasible )
            {
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if ( fixed )
               ++(*nfixedvars);
         }
      }

      SCIPdebugMsg(scip, "Deleting redundant SOS1 constraint <%s> with one variable.\n", SCIPconsGetName(cons));

      /* delete original constraint */
      assert( ! SCIPconsIsModifiable(cons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
      *success = TRUE;
   }
   /* note: there is no need to update consdata->nfixednonzeros, since the constraint is deleted as soon nfixednonzeros > 0. */
   else
   {
      /* if all variables are binary create a set packing constraint */
      if ( allvarsbinary && SCIPfindConshdlr(scip, "setppc") != NULL )
      {
         SCIP_CONS* setpackcons;

         /* create, add, and release the logicor constraint */
         SCIP_CALL( SCIPcreateConsSetpack(scip, &setpackcons, SCIPconsGetName(cons), consdata->nvars, consdata->vars,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), 
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, setpackcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &setpackcons) );

         SCIPdebugMsg(scip, "Upgrading SOS1 constraint <%s> to set packing constraint.\n", SCIPconsGetName(cons));

         /* remove the SOS1 constraint globally */
         assert( ! SCIPconsIsModifiable(cons) );
         SCIP_CALL( SCIPdelCons(scip, cons) );
         ++(*nupgdconss);
         *success = TRUE;
      }
   }

   return SCIP_OKAY;
}



/** perform one presolving round for all SOS1 constraints
 *
 *  We perform the following presolving steps.
 *
 *  - If the bounds of some variable force it to be nonzero, we can
 *    fix all other variables to zero and remove the SOS1 constraints
 *    that contain it.
 *  - If a variable is fixed to zero, we can remove the variable.
 *  - If a variable appears twice, it can be fixed to 0.
 *  - We substitute appregated variables.
 *  - Remove redundant SOS1 constraints
 *
 *  If the adjacency matrix of the conflict graph is present, then
 *  we perform the following additional presolving steps
 *
 *  - Search for larger SOS1 constraints in the conflict graph
 */
static
SCIP_RETCODE presolRoundConssSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_Bool**           adjacencymatrix,    /**< adjacency matrix of conflict graph (or NULL) */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   int                   nconss,             /**< number of SOS1 constraints */
   int                   nsos1vars,          /**< number of SOS1 variables */
   int*                  naddconss,          /**< number of added constraints */
   int*                  ndelconss,          /**< number of deleted constraints */
   int*                  nupgdconss,         /**< number of upgraded constraints */
   int*                  nfixedvars,         /**< number of fixed variables */
   int*                  nremovedvars,       /**< number of variables removed */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_DIGRAPH* vertexcliquegraph;
   SCIP_VAR** consvars;
   SCIP_Real* consweights;
   int** cliques = NULL;
   int ncliques = 0;
   int* cliquesizes = NULL;
   int* newclique = NULL;
   int* indconss = NULL;
   int* lengthconss = NULL;
   int* comsucc = NULL;
   int csize;
   int iter;
   int c;

   assert( scip != NULL );
   assert( eventhdlr != NULL );
   assert( conshdlrdata != NULL );
   assert( conflictgraph != NULL );
   assert( conss != NULL );
   assert( naddconss != NULL );
   assert( ndelconss != NULL );
   assert( nupgdconss != NULL );
   assert( nfixedvars != NULL );
   assert( nremovedvars != NULL );
   assert( result != NULL );

   /* create digraph whose nodes represent variables and cliques in the conflict graph */
   csize = MAX(1, conshdlrdata->maxextensions) * nconss;
   SCIP_CALL( SCIPcreateDigraph(scip, &vertexcliquegraph, nsos1vars + csize) );

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consweights, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquesizes, csize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newclique, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indconss, csize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lengthconss, csize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comsucc, MAX(nsos1vars, csize)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliques, csize) );

   /* get constraint indices and sort them in descending order of their lengths */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      indconss[c] = c;
      lengthconss[c] = consdata->nvars;
   }
   SCIPsortDownIntInt(lengthconss, indconss, nconss);

   /* check each constraint */
   for (iter = 0; iter < nconss; ++iter)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_Bool substituted;
      SCIP_Bool success;
      SCIP_Bool cutoff;
      int savennupgdconss;
      int savendelconss;

      SCIP_VAR** vars;
      int nvars;

      c = indconss[iter];

      assert( conss != NULL );
      assert( conss[c] != NULL );
      cons = conss[c];
      consdata = SCIPconsGetData(cons);

      assert( consdata != NULL );
      assert( consdata->nvars >= 0 );
      assert( consdata->nvars <= consdata->maxvars );
      assert( ! SCIPconsIsModifiable(cons) );
      assert( ncliques < csize );

      savendelconss = *ndelconss;
      savennupgdconss = *nupgdconss;

      /* perform one presolving round for SOS1 constraint */
      SCIP_CALL( presolRoundConsSOS1(scip, cons, consdata, eventhdlr, &substituted, &cutoff, &success, ndelconss, nupgdconss, nfixedvars, nremovedvars) );

      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( *ndelconss > savendelconss || *nupgdconss > savennupgdconss || substituted )
      {
         *result = SCIP_SUCCESS;
         continue;
      }

      if ( success )
         *result = SCIP_SUCCESS;

      /* get number of variables of constraint */
      nvars = consdata->nvars;

      /* get variables of constraint */
      vars = consdata->vars;

      if ( nvars > 1 && conshdlrdata->maxextensions != 0 )
      {
         SCIP_Bool extended = FALSE;
         int cliquesize = 0;
         int ncomsucc = 0;
         int varprobind;
         int j;

         /* get clique and size of clique */
         for (j = 0; j < nvars; ++j)
         {
            varprobind = varGetNodeSOS1(conshdlrdata, vars[j]);

            if ( varprobind >= 0 )
               newclique[cliquesize++] = varprobind;
         }

         if ( cliquesize > 1 )
         {
            cliquesizes[ncliques] = cliquesize;

            /* sort clique vertices */
            SCIPsortInt(newclique, cliquesizes[ncliques]);

            /* check if clique is contained in an already known clique */
            if ( ncliques > 0 )
            {
               int* succ;
               int nsucc;
               int v;

               varprobind = newclique[0];
               ncomsucc = SCIPdigraphGetNSuccessors(vertexcliquegraph, varprobind);
               succ = SCIPdigraphGetSuccessors(vertexcliquegraph, varprobind);

               /* get all (already processed) cliques that contain 'varpropind' */
               for (j = 0; j < ncomsucc; ++j)
               {
                  /* successors should have been sorted in a former step of the algorithm */
                  assert( j == 0 || succ[j] > succ[j-1] );
                  comsucc[j] = succ[j];
               }

               /* loop through remaining nodes of clique (case v = 0 already processed) */
               for (v = 1; v < cliquesize && ncomsucc > 0; ++v)
               {
                  varprobind = newclique[v];

                  /* get all (already processed) cliques that contain 'varpropind' */
                  nsucc = SCIPdigraphGetNSuccessors(vertexcliquegraph, varprobind);
                  succ = SCIPdigraphGetSuccessors(vertexcliquegraph, varprobind);
                  assert( succ != NULL || nsucc == 0 );

                  if ( nsucc < 1 )
                  {
                     ncomsucc = 0;
                     break;
                  }

                  /* get intersection with comsucc */
                  SCIP_CALL( SCIPcomputeArraysIntersection(comsucc, ncomsucc, succ, nsucc, comsucc, &ncomsucc) );
               }
            }

            /* if constraint is redundand then delete it */
            if ( ncomsucc > 0 )
            {
               assert( ! SCIPconsIsModifiable(cons) );
               SCIP_CALL( SCIPdelCons(scip, cons) );
               ++(*ndelconss);
               *result = SCIP_SUCCESS;
               continue;
            }

            if ( conshdlrdata->maxextensions != 0 && adjacencymatrix != NULL )
            {
               int maxextensions;
               ncomsucc = 0;

               /* determine the common successors of the vertices from the considered clique */
               SCIP_CALL( cliqueGetCommonSuccessorsSOS1(conshdlrdata, conflictgraph, newclique, vars, nvars, comsucc, &ncomsucc) );

               /* find extensions for the clique */
               maxextensions = conshdlrdata->maxextensions;
               extended = FALSE;
               SCIP_CALL( extensionOperatorSOS1(scip, conshdlrdata, adjacencymatrix, vertexcliquegraph, nsos1vars, nconss, cons, consvars, consweights,
                     TRUE, (maxextensions <= 1) ? FALSE : TRUE, cliques, &ncliques, cliquesizes, newclique, comsucc, ncomsucc, 0, -1, &maxextensions,
                     naddconss, &extended) );
            }

            /* if an extension was found for the current clique then free the old SOS1 constraint */
            if ( extended )
            {
               assert( ! SCIPconsIsModifiable(cons) );
               SCIP_CALL( SCIPdelCons(scip, cons) );
               ++(*ndelconss);
               *result = SCIP_SUCCESS;
            }
            else /* if we keep the constraint */
            {
               int cliqueind;

               cliqueind = nsos1vars + ncliques; /* index of clique in vertex-clique graph */

               /* add directed edges to the vertex-clique graph */
               assert( cliquesize >= 0 && cliquesize <= nsos1vars );
               assert( ncliques < csize );
               SCIP_CALL( SCIPallocBufferArray(scip, &cliques[ncliques], cliquesize) );/*lint !e866*/
               for (j = 0; j < cliquesize; ++j)
               {
                  cliques[ncliques][j] = newclique[j];
                  SCIP_CALL( SCIPdigraphAddArcSafe(vertexcliquegraph, cliques[ncliques][j], cliqueind, NULL) );
               }

               /* update number of maximal cliques */
               ++ncliques;
            }
         }
      }
   }

   /* free buffer arrays */
   for (c = ncliques-1; c >= 0; --c)
      SCIPfreeBufferArrayNull(scip, &cliques[c]);
   SCIPfreeBufferArrayNull(scip, &cliques);
   SCIPfreeBufferArrayNull(scip, &comsucc);
   SCIPfreeBufferArrayNull(scip, &lengthconss);
   SCIPfreeBufferArrayNull(scip, &indconss);
   SCIPfreeBufferArrayNull(scip, &newclique);
   SCIPfreeBufferArrayNull(scip, &cliquesizes);
   SCIPfreeBufferArrayNull(scip, &consweights);
   SCIPfreeBufferArrayNull(scip, &consvars);
   SCIPdigraphFree(&vertexcliquegraph);

   return SCIP_OKAY;
}


/** performs implication graph analysis
 *
 *  Tentatively fixes a variable to nonzeero and extracts consequences from it:
 *  - adds (possibly new) complementarity constraints to the problem if variables are implied to be zero
 *  - returns that the subproblem is infeasible if the domain of a variable turns out to be empty
 */
static
SCIP_RETCODE performImplicationGraphAnalysis(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_VAR**            totalvars,          /**< problem and SOS1 variables */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph (@p j is successor of @p i if and only if \f$ x_i\not = 0 \Rightarrow x_j\not = 0\f$) */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_Bool**           adjacencymatrix,    /**< adjacencymatrix of the conflict graph (only lower half filled) */
   int                   givennode,          /**< node of the conflict graph */
   int                   nonznode,           /**< node of the conflict graph that is implied to be nonzero if given node is nonzero */
   SCIP_Real*            impllbs,            /**< current lower variable bounds if given node is nonzero (update possible) */
   SCIP_Real*            implubs,            /**< current upper variable bounds if given node is nonzero (update possible) */
   SCIP_Bool*            implnodes,          /**< indicates which variables are currently implied to be nonzero if given node is nonzero (update possible) */
   int*                  naddconss,          /**< pointer to store number of added SOS1 constraints */
   int*                  probingdepth,       /**< pointer to store current probing depth */
   SCIP_Bool*            infeasible          /**< pointer to store whether the subproblem gets infeasible if variable to 'nonznode' is nonzero */
   )
{
   SCIP_SUCCDATA** succdatas;
   int succnode;
   int* succ;
   int nsucc;
   int s;

   assert( nonznode >= 0 && nonznode < SCIPdigraphGetNNodes(conflictgraph) );

   /* check probing depth */
   if ( conshdlrdata->depthimplanalysis >= 0 && *probingdepth >= conshdlrdata->depthimplanalysis )
      return SCIP_OKAY;
   ++(*probingdepth);

   /* get successors of 'nonznode' in the conflict graph */
   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, nonznode);
   succ = SCIPdigraphGetSuccessors(conflictgraph, nonznode);

   /* loop through neighbors of 'nonznode' in the conflict graph; these variables are implied to be zero */
   for (s = 0; s < nsucc; ++s)
   {
      succnode = succ[s];

      /* if the current variable domain of the successor node does not contain the value zero then return that the problem is infeasible
       * else if 'succnode' is not already complementary to 'givennode' then add a new complementarity constraint */
      if ( givennode == succnode || SCIPisFeasPositive(scip, impllbs[succnode]) || SCIPisFeasNegative(scip, implubs[succnode]) )
      {
	*infeasible = TRUE;
	return SCIP_OKAY;
      }
      else if ( ! isConnectedSOS1(adjacencymatrix, NULL, givennode, succnode) )
      {
	char namesos[SCIP_MAXSTRLEN];
	SCIP_CONS* soscons = NULL;
	SCIP_VAR* var1;
	SCIP_VAR* var2;

	/* update implied bounds of succnode */
	impllbs[succnode] = 0;
	implubs[succnode] = 0;

	/* add arcs to the conflict graph */
	SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraph, givennode, succnode, NULL) );
	SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraph, succnode, givennode, NULL) );

	/* resort successors */
	SCIPsortInt(SCIPdigraphGetSuccessors(conflictgraph, givennode), SCIPdigraphGetNSuccessors(conflictgraph, givennode));
	SCIPsortInt(SCIPdigraphGetSuccessors(conflictgraph, succnode), SCIPdigraphGetNSuccessors(conflictgraph, succnode));

	/* update adjacencymatrix */
	if ( givennode > succnode )
	  adjacencymatrix[givennode][succnode] = 1;
	else
	  adjacencymatrix[succnode][givennode] = 1;

	var1 = SCIPnodeGetVarSOS1(conflictgraph, givennode);
	var2 = SCIPnodeGetVarSOS1(conflictgraph, succnode);

	/* create SOS1 constraint */
	assert( SCIPgetDepth(scip) == 0 );
	(void) SCIPsnprintf(namesos, SCIP_MAXSTRLEN, "presolved_sos1_%s_%s", SCIPvarGetName(var1), SCIPvarGetName(var2) );
	SCIP_CALL( SCIPcreateConsSOS1(scip, &soscons, namesos, 0, NULL, NULL, TRUE, TRUE, TRUE, FALSE, TRUE,
				      FALSE, FALSE, FALSE, FALSE) );

	/* add variables to SOS1 constraint */
	SCIP_CALL( addVarSOS1(scip, soscons, conshdlrdata, var1, 1.0) );
	SCIP_CALL( addVarSOS1(scip, soscons, conshdlrdata, var2, 2.0) );

	/* add constraint */
	SCIP_CALL( SCIPaddCons(scip, soscons) );

	/* release constraint */
	SCIP_CALL( SCIPreleaseCons(scip, &soscons) );

	++(*naddconss);
      }
   }

   /* by construction: nodes of SOS1 variables are equal for conflict graph and implication graph */
   assert( nonznode == (int) (size_t) SCIPhashmapGetImage(implhash, SCIPnodeGetVarSOS1(conflictgraph, nonznode)) );
   succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, nonznode);
   nsucc = SCIPdigraphGetNSuccessors(implgraph, nonznode);
   succ = SCIPdigraphGetSuccessors(implgraph, nonznode);

   /* go further in implication graph */
   for (s = 0; s < nsucc; ++s)
   {
     SCIP_SUCCDATA* data;
     int oldprobingdepth;

     succnode = succ[s];
     data = succdatas[s];
     oldprobingdepth = *probingdepth;

     /* if current lower bound is smaller than implied lower bound */
     if ( SCIPisFeasLT(scip, impllbs[succnode], data->lbimpl) )
     {
	impllbs[succnode] = data->lbimpl;

	/* if node is SOS1 and implied to be nonzero for the first time, then this recursively may imply further bound changes */
	if ( varGetNodeSOS1(conshdlrdata, totalvars[succnode]) >= 0 && ! implnodes[succnode] && SCIPisFeasPositive(scip, data->lbimpl) )
	{
	   /* by construction: nodes of SOS1 variables are equal for conflict graph and implication graph */
	   assert( succnode == (int) (size_t) SCIPhashmapGetImage(implhash, SCIPnodeGetVarSOS1(conflictgraph, succnode)) );
	   implnodes[succnode] = TRUE; /* in order to avoid cycling */
	   SCIP_CALL( performImplicationGraphAnalysis(scip, conshdlrdata, conflictgraph, totalvars, implgraph, implhash, adjacencymatrix, givennode, succnode, impllbs, implubs, implnodes, naddconss, probingdepth, infeasible) );
	   *probingdepth = oldprobingdepth;

	   /* return if the subproblem is known to be infeasible */
	   if ( *infeasible )
	      return SCIP_OKAY;
	}
     }

     /* if current upper bound is larger than implied upper bound */
     if ( SCIPisFeasGT(scip, implubs[succnode], data->ubimpl) )
     {
	implubs[succnode] = data->ubimpl;

	/* if node is SOS1 and implied to be nonzero for the first time, then this recursively may imply further bound changes */
	if ( varGetNodeSOS1(conshdlrdata, totalvars[succnode]) >= 0 && ! implnodes[succnode] && SCIPisFeasNegative(scip, data->ubimpl) )
	{
	   /* by construction: nodes of SOS1 variables are equal for conflict graph and implication graph */
	   assert( succnode == (int) (size_t) SCIPhashmapGetImage(implhash, SCIPnodeGetVarSOS1(conflictgraph, succnode)) );
	   implnodes[succnode] = TRUE; /* in order to avoid cycling */
	   SCIP_CALL( performImplicationGraphAnalysis(scip, conshdlrdata, conflictgraph, totalvars, implgraph, implhash, adjacencymatrix, givennode, succnode, impllbs, implubs, implnodes, naddconss, probingdepth, infeasible) );
	   *probingdepth = oldprobingdepth;

	   /* return if the subproblem is known to be infeasible */
	   if ( *infeasible )
	      return SCIP_OKAY;
	}
     }
   }

   return SCIP_OKAY;
}


/** returns whether node is implied to be zero; this information is taken from the input array 'implnodes' */
static
SCIP_Bool isImpliedZero(
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_Bool*            implnodes,          /**< implnodes[i] = TRUE if the SOS1 variable corresponding to node i in the implication graph is implied to be nonzero */
   int                   node                /**< node of the conflict graph (or -1) */
   )
{
   int* succ;
   int nsucc;
   int s;

   if ( node < 0 )
      return FALSE;

   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, node);
   succ = SCIPdigraphGetSuccessors(conflictgraph, node);

   /* check whether any successor is implied to be nonzero */
   for (s = 0; s < nsucc; ++s)
   {
      if ( implnodes[succ[s]] )
         return TRUE;
   }

   return FALSE;
}


/** updates arc data of implication graph */
static
SCIP_RETCODE updateArcData(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_VAR**            totalvars,          /**< problem and SOS1 variables */
   SCIP_VAR*             varv,               /**< variable that is assumed to be nonzero */
   SCIP_VAR*             varw,               /**< implication variable */
   SCIP_Real             lb,                 /**< old lower bound of \f$x_w\f$ */
   SCIP_Real             ub,                 /**< old upper bound of \f$x_w\f$ */
   SCIP_Real             newbound,           /**< new bound of \f$x_w\f$ */
   SCIP_Bool             lower,              /**< whether to consider lower bound implication (otherwise upper bound) */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   SCIP_Bool*            update,             /**< pointer to store whether implication graph has been updated */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility has been detected */
   )
{
   SCIP_SUCCDATA** succdatas;
   SCIP_SUCCDATA* data = NULL;
   int nsucc;
   int* succ;
   int indv;
   int indw;
   int s;

   assert( scip != NULL );
   assert( implgraph != NULL );
   assert( implhash != NULL );
   assert( totalvars != NULL );
   assert( varv != NULL );
   assert( varw != NULL );

   /* if x_v != 0 turns out to be infeasible then fix x_v = 0 */
   if ( ( lower && SCIPisFeasLT(scip, ub, newbound) ) || ( ! lower && SCIPisFeasGT(scip, lb, newbound) ) )
   {
      SCIP_Bool infeasible1;
      SCIP_Bool infeasible2;
      SCIP_Bool tightened1;
      SCIP_Bool tightened2;

      SCIP_CALL( SCIPtightenVarLb(scip, varv, 0.0, FALSE, &infeasible1, &tightened1) );
      SCIP_CALL( SCIPtightenVarUb(scip, varv, 0.0, FALSE, &infeasible2, &tightened2) );

      if ( infeasible1 || infeasible2 )
      {
         SCIPdebugMsg(scip, "detected infeasibility while trying to fix variable <%s> to zero\n", SCIPvarGetName(varv));
         *infeasible = TRUE;
      }

      if ( tightened1 || tightened2 )
      {
         SCIPdebugMsg(scip, "fixed variable %s from lb = %f and ub = %f to 0.0 \n", SCIPvarGetName(varv), lb, ub);
         ++(*nchgbds);
      }
   }

   /* get successor information */
   indv = (int) (size_t) SCIPhashmapGetImage(implhash, varv); /* get index of x_v in implication graph */
   assert( (int) (size_t) SCIPhashmapGetImage(implhash, totalvars[indv]) == indv );
   succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, indv);
   nsucc = SCIPdigraphGetNSuccessors(implgraph, indv);
   succ = SCIPdigraphGetSuccessors(implgraph, indv);

   /* search for nodew in existing successors. If this is the case then check whether the lower implication bound may be updated ... */
   indw = (int) (size_t) SCIPhashmapGetImage(implhash, varw);
   assert( (int) (size_t) SCIPhashmapGetImage(implhash, totalvars[indw]) == indw );
   for (s = 0; s < nsucc; ++s)
   {
      if ( succ[s] == indw )
      {
         data = succdatas[s];
         assert( data != NULL );
         if ( lower && SCIPisFeasLT(scip, data->lbimpl, newbound) )
         {
            if ( SCIPvarIsIntegral(varw) )
               data->lbimpl = SCIPceil(scip, newbound);
            else
               data->lbimpl = newbound;

            *update = TRUE;
            SCIPdebugMsg(scip, "updated to implication %s != 0 -> %s >= %f\n", SCIPvarGetName(varv), SCIPvarGetName(varw), newbound);
         }
         else if ( ! lower && SCIPisFeasGT(scip, data->ubimpl, newbound) )
         {
            if ( SCIPvarIsIntegral(varw) )
               data->ubimpl = SCIPfloor(scip, newbound);
            else
               data->ubimpl = newbound;

            *update = TRUE;
            SCIPdebugMsg(scip, "updated to implication %s != 0 -> %s >= %f\n", SCIPvarGetName(varv), SCIPvarGetName(varw), newbound);
         }
         break;
      }
   }

   /* ..., otherwise if there does not exist an arc between indv and indw already, then create one and add implication */
   if ( s == nsucc )
   {
      assert( data == NULL );
      SCIP_CALL( SCIPallocBlockMemory(scip, &data) );
      if ( lower )
      {
         data->lbimpl = newbound;
         data->ubimpl = ub;
         SCIPdebugMsg(scip, "add implication %s != 0 -> %s >= %f\n", SCIPvarGetName(varv), SCIPvarGetName(varw), newbound);
      }
      else
      {
         data->lbimpl = lb;
         data->ubimpl = newbound;
         SCIPdebugMsg(scip, "add implication %s != 0 -> %s <= %f\n", SCIPvarGetName(varv), SCIPvarGetName(varw), newbound);
      }
      SCIP_CALL( SCIPdigraphAddArc(implgraph, indv, indw, (void*)data) );
      *update = TRUE;
   }

   return SCIP_OKAY;
}


/** updates implication graph
 *
 *  Assume the variable from the input is nonzero. If this implies that some other variable is also nonzero, then
 *  store this information in an implication graph
 */
static
SCIP_RETCODE updateImplicationGraphSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_Bool**           adjacencymatrix,    /**< adjacency matrix of conflict graph (lower half) */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph (@p j is successor of @p i if and only if \f$ x_i\not = 0 \Rightarrow x_j\not = 0\f$) */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_Bool*            implnodes,          /**< implnodes[i] = TRUE if the SOS1 variable corresponding to node i in the implication graph is implied to be nonzero */
   SCIP_VAR**            totalvars,          /**< problem and SOS1 variables */
   int**                 cliquecovers,       /**< clique covers of linear constraint */
   int*                  cliquecoversizes,   /**< size of clique covers */
   int*                  varincover,         /**< array with varincover[i] = cover of SOS1 index @p i */
   SCIP_VAR**            vars,               /**< variables to be checked */
   SCIP_Real*            coefs,              /**< coefficients of variables in linear constraint */
   int                   nvars,              /**< number of variables to be checked */
   SCIP_Real*            bounds,             /**< bounds of variables */
   SCIP_VAR*             var,                /**< variable that is assumed to be nonzero */
   SCIP_Real             bound,              /**< bound of variable */
   SCIP_Real             boundnonzero,       /**< bound of variable if it is known to be nonzero if infinity values are not summarized */
   int                   ninftynonzero,      /**< number of times infinity/-infinity has to be summarized to boundnonzero */
   SCIP_Bool             lower,              /**< TRUE if lower bounds are consideres; FALSE for upper bounds */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   SCIP_Bool*            update,             /**< pointer to store whether implication graph has been updated */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility has been detected */
   )
{
   int nodev;
   int w;

   assert( update != NULL );

   /* update implication graph if possible */
   *update = FALSE;
   *infeasible = FALSE;
   nodev = varGetNodeSOS1(conshdlrdata, var); /* possibly -1 if var is not involved in an SOS1 constraint */

   /* if nodev is an index of an SOS1 variable and at least one lower bound of a variable that is not x_v is infinity */
   if ( nodev < 0 || SCIPisInfinity(scip, REALABS(bound)) || ninftynonzero > 1 )
      return SCIP_OKAY;

   /* for every variable x_w: compute upper bound of a_w * x_w if x_v is known to be nonzero */
   for (w = 0; w < nvars; ++w)
   {
      int newninftynonzero;
      SCIP_Bool implinfty = FALSE;
      int nodew;

      /* get node of x_w in conflict graph: nodew = -1 if it is no SOS1 variable */
      nodew = varGetNodeSOS1(conshdlrdata, vars[w]);

      newninftynonzero = ninftynonzero;

      /* variable should not be fixed to be already zero (note x_v is fixed to be nonzero by assumption) */
      if ( nodew < 0 || ( nodev != nodew && ! isConnectedSOS1(adjacencymatrix, NULL, nodev, nodew) && ! isImpliedZero(conflictgraph, implnodes, nodew) ) )
      {
         SCIP_Real implbound;
         SCIP_Bool implcoverw;
         int nodecliq;
         int indcliq;
         int ind;
         int j;

         /* boundnonzero is the bound of x_v if x_v is nonzero we use this information to get a bound of x_w if x_v is
          * nonzero; therefore, we have to perform some recomputations */
         implbound = boundnonzero - bound;
         ind = varincover[w];
         assert( cliquecoversizes[ind] > 0 );

         implcoverw = FALSE;
         for (j = 0; j < cliquecoversizes[ind]; ++j)
         {
            indcliq = cliquecovers[ind][j];
            assert( 0 <= indcliq && indcliq < nvars );

            nodecliq = varGetNodeSOS1(conshdlrdata, vars[indcliq]); /* possibly -1 if variable is not involved in an SOS1 constraint */

            /* if nodecliq is not a member of an SOS1 constraint or the variable corresponding to nodecliq is not implied to be zero if x_v != 0  */
            if ( nodecliq < 0 || (! isConnectedSOS1(adjacencymatrix, NULL, nodev, nodecliq) && ! isImpliedZero(conflictgraph, implnodes, nodecliq) ) )
            {
               if ( indcliq == w )
               {
                  if ( !SCIPisInfinity(scip, REALABS(bounds[w])) && !SCIPisInfinity(scip, REALABS(implbound + bounds[w])) )
                     implbound += bounds[w];
                  else
                     --newninftynonzero;
                  implcoverw = TRUE;
               }
               else if ( implcoverw )
               {
                  if ( SCIPisInfinity(scip, REALABS(bounds[indcliq])) || SCIPisInfinity(scip, REALABS(implbound - bounds[indcliq])) )
                     implinfty = TRUE;
                  else
                     implbound -= bounds[indcliq];
                  break;
               }
               else
               {
                  if ( SCIPisInfinity(scip, REALABS(bounds[indcliq])) )
                     implinfty = TRUE;
                  break;
               }
            }
         }

         /* check whether x_v != 0 implies a bound change of x_w */
         if ( ! implinfty && newninftynonzero == 0 )
         {
            SCIP_Real newbound;
            SCIP_Real coef;
            SCIP_Real lb;
            SCIP_Real ub;

            lb = SCIPvarGetLbLocal(vars[w]);
            ub = SCIPvarGetUbLocal(vars[w]);
            coef = coefs[w];

            if ( SCIPisFeasZero(scip, coef) )
               continue;

            newbound = implbound / coef;

            /* check if an implication can be added/updated or assumption x_v != 0 is infeasible */
            if ( lower )
            {
               if ( SCIPisFeasPositive(scip, coef) && SCIPisFeasLT(scip, lb, newbound) )
               {
                  SCIP_CALL( updateArcData(scip, implgraph, implhash, totalvars, var, vars[w], lb, ub, newbound, TRUE, nchgbds, update, infeasible) );
               }
               else if ( SCIPisFeasNegative(scip, coef) && SCIPisFeasGT(scip, ub, newbound) )
               {
                  SCIP_CALL( updateArcData(scip, implgraph, implhash, totalvars, var, vars[w], lb, ub, newbound, FALSE, nchgbds, update, infeasible) );
               }
            }
            else
            {
               if ( SCIPisFeasPositive(scip, coef) && SCIPisFeasGT(scip, ub, newbound) )
               {
                  SCIP_CALL( updateArcData(scip, implgraph, implhash, totalvars, var, vars[w], lb, ub, newbound, FALSE, nchgbds, update, infeasible) );
               }
               else if ( SCIPisFeasNegative(scip, coef) && SCIPisFeasLT(scip, lb, newbound) )
               {
                  SCIP_CALL( updateArcData(scip, implgraph, implhash, totalvars, var, vars[w], lb, ub, newbound, TRUE, nchgbds, update, infeasible) );
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** search new disjoint clique that covers given node
 *
 *  For a given vertex @p v search for a clique of the conflict graph induced by the variables of a linear constraint that
 *  - covers @p v and
 *  - has an an empty intersection with already computed clique cover.
 */
static
SCIP_RETCODE computeVarsCoverSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraphroot,  /**< conflict graph of the root node (nodes: 1, ..., @p nsos1vars) */
   SCIP_DIGRAPH*         conflictgraphlin,   /**< conflict graph of linear constraint (nodes: 1, ..., @p nlinvars) */
   SCIP_VAR**            linvars,            /**< variables in linear constraint */
   SCIP_Bool*            coveredvars,        /**< states which variables of the linear constraint are currently covered by a clique */
   int*                  clique,             /**< array to store new clique in cover */
   int*                  cliquesize,         /**< pointer to store the size of @p clique */
   int                   v,                  /**< position of variable in linear constraint that should be covered */
   SCIP_Bool             considersolvals     /**< TRUE if largest auxiliary bigM values of variables should be prefered */
   )
{
   int nsucc;
   int s;

   assert( conflictgraphlin != NULL );
   assert( linvars != NULL );
   assert( coveredvars != NULL );
   assert( clique != NULL );
   assert( cliquesize != NULL );

   assert( ! coveredvars[v] );  /* we should produce a new clique */

   /* add index 'v' to the clique cover */
   clique[0] = v;
   *cliquesize = 1;

   nsucc = SCIPdigraphGetNSuccessors(conflictgraphlin, v);
   if ( nsucc > 0 )
   {
      int* extensions;
      int nextensions = 0;
      int nextensionsnew;
      int succnode;
      int* succ;

      /* allocate buffer array */
      SCIP_CALL( SCIPallocBufferArray(scip, &extensions, nsucc) );

      succ = SCIPdigraphGetSuccessors(conflictgraphlin, v);

      /* compute possible extensions for the clique cover */
      for (s = 0; s < nsucc; ++s)
      {
         succnode = succ[s];
         if ( ! coveredvars[succnode] )
            extensions[nextensions++] = succ[s];
      }

      /* while there exist possible extensions for the clique cover */
      while ( nextensions > 0 )
      {
         int bestindex = -1;

         if ( considersolvals )
         {
            SCIP_Real bestbigMval;
            SCIP_Real bigMval;

            bestbigMval = -SCIPinfinity(scip);

            /* search for the extension with the largest absolute value of its LP relaxation solution value */
            for (s = 0; s < nextensions; ++s)
            {
               bigMval = nodeGetSolvalBinaryBigMSOS1(scip, conflictgraphroot, NULL, extensions[s]);
               if ( SCIPisFeasLT(scip, bestbigMval, bigMval) )
               {
                  bestbigMval = bigMval;
                  bestindex = extensions[s];
               }
            }
         }
         else
            bestindex = extensions[0];

         assert( bestindex != -1 );

         /* add bestindex to the clique cover */
         clique[(*cliquesize)++] = bestindex;

         /* compute new 'extensions' array */
         nextensionsnew = 0;
         for (s = 0; s < nextensions; ++s)
         {
            if ( s != bestindex && isConnectedSOS1(NULL, conflictgraphlin, bestindex, extensions[s]) )
               extensions[nextensionsnew++] = extensions[s];
         }
         nextensions = nextensionsnew;
      }

      /* free buffer array */
      SCIPfreeBufferArray(scip, &extensions);
   }

   /* mark covered indices */
   for (s = 0; s < *cliquesize; ++s)
   {
      int ind;

      ind = clique[s];
      assert( 0 <= ind );
      assert( ! coveredvars[ind] );
      coveredvars[ind] = TRUE;
   }

   return SCIP_OKAY;
}


/** try to tighten upper and lower bounds for variables */
static
SCIP_RETCODE tightenVarsBoundsSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph (@p j is successor of @p i if and only if \f$ x_i\not = 0 \f$ implies a new lower/upper bound for \f$ x_j\f$) */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_Bool**           adjacencymatrix,    /**< adjacencymatrix of conflict graph */
   SCIP_VAR**            totalvars,          /**< problem and SOS1 vars */
   int                   ntotalvars,         /**< number of problem and SOS1 variables*/
   int                   nsos1vars,          /**< number of SOS1 variables */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   SCIP_Bool*            implupdate,         /**< pointer to store whether the implication graph has been updated in this function call */
   SCIP_Bool*            cutoff              /**< pointer to store if current nodes LP is infeasible  */
   )
{
   SCIP_CONSHDLR* conshdlrlinear;
   SCIP_CONS** linearconss;
   int nlinearconss;

   SCIP_Bool* implnodes = NULL;     /* implnodes[i] = TRUE if the SOS1 variable corresponding to node i in the implication graph is implied to be nonzero */
   SCIP_Bool* coveredvars = NULL;   /* coveredvars[i] = TRUE if variable with index i is covered by the clique cover */
   int* varindincons = NULL;        /* varindincons[i] = position of SOS1 index i in linear constraint (-1 if x_i is not involved in linear constraint) */

   SCIP_VAR** trafolinvars = NULL;  /* variables of transformed linear constraints without (multi)aggregated variables */
   int ntrafolinvars = 0;
   SCIP_Real* trafolinvals = NULL;
   SCIP_Real* trafoubs = NULL;
   SCIP_Real* trafolbs = NULL;
   SCIP_Real traforhs;
   SCIP_Real trafolhs;

   SCIP_VAR** sos1linvars = NULL;  /* variables that are not contained in linear constraint, but are in conflict with a variable from the linear constraint */
   int nsos1linvars;
   int c;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( adjacencymatrix != NULL );
   assert( nchgbds != NULL );
   assert( cutoff != NULL );

   *cutoff = FALSE;
   *implupdate = FALSE;

   /* get constraint handler data of linear constraints */
   conshdlrlinear = SCIPfindConshdlr(scip, "linear");
   if ( conshdlrlinear == NULL )
      return SCIP_OKAY;

   /* get linear constraints and number of linear constraints */
   nlinearconss = SCIPconshdlrGetNConss(conshdlrlinear);
   linearconss = SCIPconshdlrGetConss(conshdlrlinear);

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &sos1linvars, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &implnodes, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varindincons, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coveredvars, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &trafoubs, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &trafolbs, ntotalvars) );

   /* for every linear constraint and every SOS1 variable */
   for (c = 0; c < nlinearconss + nsos1vars && ! (*cutoff); ++c)
   {
      SCIP_DIGRAPH* conflictgraphlin;
      int** cliquecovers = NULL;           /* clique covers of indices of variables in linear constraint */
      int* cliquecoversizes = NULL;        /* size of each cover */
      int ncliquecovers;
      SCIP_Real* cliquecovervals = NULL;
      int* varincover = NULL;              /* varincover[i] = cover of SOS1 index i */

      int v;
      int i;
      int j;

      /* get transformed linear constraints (without aggregated variables) */
      if ( c < nlinearconss )
      {
         SCIP_VAR** origlinvars;
         int noriglinvars;
         SCIP_Real* origlinvals;
         SCIP_Real origrhs;
         SCIP_Real origlhs;
         SCIP_Real constant;
         int requiredsize;

         /* get data of linear constraint */
         noriglinvars = SCIPgetNVarsLinear(scip, linearconss[c]);
         origlinvars = SCIPgetVarsLinear(scip, linearconss[c]);
         origlinvals = SCIPgetValsLinear(scip, linearconss[c]);
         origrhs = SCIPgetRhsLinear(scip, linearconss[c]);
         origlhs = SCIPgetLhsLinear(scip, linearconss[c]);

         if ( noriglinvars < 1 )
            continue;
         assert( origlinvars != NULL );
         assert( origlinvals != NULL );

         /* copy variables and coefficients of linear constraint */
         SCIP_CALL( SCIPduplicateBufferArray(scip, &trafolinvars, origlinvars, noriglinvars) );
         SCIP_CALL( SCIPduplicateBufferArray(scip, &trafolinvals, origlinvals, noriglinvars) );
         ntrafolinvars = noriglinvars;

         /* transform linear constraint */
         constant = 0.0;
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, trafolinvars, trafolinvals, &ntrafolinvars, noriglinvars, &constant, &requiredsize, TRUE) );
         if( requiredsize > ntrafolinvars )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &trafolinvars, requiredsize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &trafolinvals, requiredsize) );

            SCIP_CALL( SCIPgetProbvarLinearSum(scip, trafolinvars, trafolinvals, &ntrafolinvars, requiredsize, &constant, &requiredsize, TRUE) );
            assert( requiredsize <= ntrafolinvars );
         }
         trafolhs = origlhs - constant;
         traforhs = origrhs - constant;
      }
      else
      {
         SCIP_VAR* var;

         var = SCIPnodeGetVarSOS1(conflictgraph, c-nlinearconss);

         if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED )
         {
            SCIP_Real constant;

            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvals, 2) );

            constant = SCIPvarGetAggrConstant(var);
            trafolinvars[0] = SCIPvarGetAggrVar(var);
            trafolinvals[0] = SCIPvarGetAggrScalar(var);
            trafolinvars[1] = var;
            trafolinvals[1] = -1.0;
            trafolhs = -constant;
            traforhs = -constant;
            ntrafolinvars = 2;
         }
         else if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
         {
            SCIP_Real* scalars;
            SCIP_VAR** agrvars;
            SCIP_Real constant;
            int nagrvars;

            nagrvars = SCIPvarGetMultaggrNVars(var);

            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvars, nagrvars+1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvals, nagrvars+1) );

            agrvars = SCIPvarGetMultaggrVars(var);
            scalars = SCIPvarGetMultaggrScalars(var);
            constant = SCIPvarGetMultaggrConstant(var);

            for (v = 0; v < nagrvars; ++v)
            {
               trafolinvars[v] = agrvars[v];
               trafolinvals[v] = scalars[v];
            }
            trafolinvars[nagrvars] = var;
            trafolinvals[nagrvars] = -1.0;
            trafolhs = -constant;
            traforhs = -constant;
            ntrafolinvars = nagrvars + 1;
         }
         else if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
         {
            SCIP_VAR* negvar;
            SCIP_Real negcons;

            /* get negation variable and negation offset */
            negvar = SCIPvarGetNegationVar(var);
            negcons = SCIPvarGetNegationConstant(var);

            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvals, 2) );

            trafolinvars[0] = negvar;
            trafolinvars[1] = var;
            trafolinvals[0] = 1.0;
            trafolinvals[1] = 1.0;
            trafolhs = negcons;
            traforhs = negcons;
            ntrafolinvars = 2;
         }
         else
            continue;
      }

      if ( ntrafolinvars == 0 )
      {
         SCIPfreeBufferArray(scip, &trafolinvars);
         SCIPfreeBufferArray(scip, &trafolinvals);
         continue;
      }

      /* compute lower and upper bounds of each term a_i * x_i of transformed constraint */
      for (v = 0; v < ntrafolinvars; ++v)
      {
         SCIP_Real lb = SCIPvarGetLbLocal(trafolinvars[v]);
         SCIP_Real ub = SCIPvarGetUbLocal(trafolinvars[v]);

         if ( trafolinvals[v] < 0.0 )
         {
            SCIP_Real temp;

            temp = lb;
            lb = ub;
            ub = temp;
         }

         assert(!SCIPisInfinity(scip, REALABS(trafolinvals[v])));

         if ( SCIPisInfinity(scip, REALABS(lb)) || SCIPisInfinity(scip, REALABS(lb * trafolinvals[v])) )
            trafolbs[v] = -SCIPinfinity(scip);
         else
            trafolbs[v] = lb * trafolinvals[v];

         if ( SCIPisInfinity(scip, REALABS(ub)) || SCIPisInfinity(scip, REALABS(ub * trafolinvals[v])) )
            trafoubs[v] = SCIPinfinity(scip);
         else
            trafoubs[v] = ub * trafolinvals[v];
      }

      /* initialization: mark all the SOS1 variables as 'not a member of the linear constraint' */
      for (v = 0; v < nsos1vars; ++v)
         varindincons[v] = -1;

      /* save position of SOS1 variables in linear constraint */
      for (v = 0; v < ntrafolinvars; ++v)
      {
         int node;

         node = varGetNodeSOS1(conshdlrdata, trafolinvars[v]);

         if ( node >= 0 )
            varindincons[node] = v;
      }

      /* create conflict graph of linear constraint */
      SCIP_CALL( SCIPcreateDigraph(scip, &conflictgraphlin, ntrafolinvars) );
      SCIP_CALL( genConflictgraphLinearCons(conshdlrdata, conflictgraphlin, conflictgraph, trafolinvars, ntrafolinvars, varindincons) );

      /* mark all the variables as 'not covered by some clique cover' */
      for (i = 0; i < ntrafolinvars; ++i)
         coveredvars[i] = FALSE;

      /* allocate buffer array */
      SCIP_CALL( SCIPallocBufferArray(scip, &cliquecovervals, ntrafolinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cliquecoversizes, ntrafolinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cliquecovers, ntrafolinvars) );

      /* compute distinct cliques that cover all the variables of the linear constraint */
      ncliquecovers = 0;
      for (v = 0; v < ntrafolinvars; ++v)
      {
         /* if variable is not already covered by an already known clique cover */
         if ( ! coveredvars[v] )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(cliquecovers[ncliquecovers]), ntrafolinvars) ); /*lint !e866*/
            SCIP_CALL( computeVarsCoverSOS1(scip, conflictgraph, conflictgraphlin, trafolinvars, coveredvars, cliquecovers[ncliquecovers], &(cliquecoversizes[ncliquecovers]), v, FALSE) );
            ++ncliquecovers;
         }
      }

      /* free conflictgraph */
      SCIPdigraphFree(&conflictgraphlin);

      /* compute variables that are not contained in transformed linear constraint, but are in conflict with a variable from the transformed linear constraint */
      nsos1linvars = 0;
      for (v = 0; v < ntrafolinvars; ++v)
      {
         int nodev;

         nodev = varGetNodeSOS1(conshdlrdata, trafolinvars[v]);

         /* if variable is an SOS1 variable */
         if ( nodev >= 0 )
         {
            int succnode;
            int nsucc;
            int* succ;
            int s;

            succ = SCIPdigraphGetSuccessors(conflictgraph, nodev);
            nsucc = SCIPdigraphGetNSuccessors(conflictgraph, nodev);

            for (s = 0; s < nsucc; ++s)
            {
               succnode = succ[s];

               /* if variable is not a member of linear constraint and not already listed in the array sos1linvars */
               if ( varindincons[succnode] == -1 )
               {
                  sos1linvars[nsos1linvars] = SCIPnodeGetVarSOS1(conflictgraph, succnode);
                  varindincons[succnode] = -2; /* mark variable as listed in array sos1linvars */
                  ++nsos1linvars;
               }
            }
         }
      }


      /* try to tighten lower bounds */

      /* sort each cliquecover array in ascending order of the lower bounds of a_i * x_i; fill vector varincover */
      SCIP_CALL( SCIPallocBufferArray(scip, &varincover, ntrafolinvars) );
      for (i = 0; i < ncliquecovers; ++i)
      {
         for (j = 0; j < cliquecoversizes[i]; ++j)
         {
            int ind = cliquecovers[i][j];

            varincover[ind] = i;
            cliquecovervals[j] = trafoubs[ind];
         }
         SCIPsortDownRealInt(cliquecovervals, cliquecovers[i], cliquecoversizes[i]);
      }

      /* for every variable in transformed constraint: try lower bound tightening */
      for (v = 0; v < ntrafolinvars + nsos1linvars; ++v)
      {
         SCIP_Real newboundnonzero; /* new bound of a_v * x_v if we assume that x_v != 0 */
         SCIP_Real newboundnores;   /* new bound of a_v * x_v if we assume that x_v = 0 is possible */
         SCIP_Real newbound;        /* resulting new bound of x_v */
         SCIP_VAR* var;
         SCIP_Real trafoubv;
         SCIP_Real linval;
         SCIP_Real ub;
         SCIP_Real lb;
         SCIP_Bool tightened;
         SCIP_Bool infeasible;
         SCIP_Bool inftynores = FALSE;
         SCIP_Bool update;
         int ninftynonzero = 0;
         int nodev;
         int w;

         if ( v < ntrafolinvars )
         {
            var = trafolinvars[v];
            trafoubv = trafoubs[v];
         }
         else
         {
            assert( v >= ntrafolinvars );
            var = sos1linvars[v-ntrafolinvars];/*lint !e679*/
            trafoubv = 0.0;
         }

         ub = SCIPvarGetUbLocal(var);
         lb = SCIPvarGetLbLocal(var);

         if ( SCIPisInfinity(scip, -trafolhs) || SCIPisZero(scip, ub - lb) )
            continue;

         newboundnonzero = trafolhs;
         newboundnores = trafolhs;
         nodev = varGetNodeSOS1(conshdlrdata, var); /* possibly -1 if var is not involved in an SOS1 constraint */
         assert( nodev < nsos1vars );

         /* determine incidence vector of implication variables */
         for (w = 0; w < nsos1vars; ++w)
            implnodes[w] = FALSE;
         SCIP_CALL( getSOS1Implications(scip, conshdlrdata, totalvars, implgraph, implhash, implnodes, (int) (size_t) SCIPhashmapGetImage(implhash, var)) );

         /* compute new bound */
         for (i = 0; i < ncliquecovers; ++i)
         {
            int indcliq;
            int nodecliq;

            assert( cliquecoversizes[i] > 0 );

            indcliq = cliquecovers[i][0];
            assert( 0 <= indcliq && indcliq < ntrafolinvars );

            /* determine maximum without index v (note that the array 'cliquecovers' is sorted by the values of trafoub in non-increasing order) */
            if ( v != indcliq )
            {
               if ( SCIPisInfinity(scip, trafoubs[indcliq]) || SCIPisInfinity(scip, REALABS(newboundnores - trafoubs[indcliq])) )
                  inftynores = TRUE;
               else
                  newboundnores -= trafoubs[indcliq];
            }
            else if ( cliquecoversizes[i] > 1 )
            {
               assert( 0 <= cliquecovers[i][1] && cliquecovers[i][1] < ntrafolinvars );
               if ( SCIPisInfinity(scip, trafoubs[cliquecovers[i][1]]) || SCIPisInfinity(scip, REALABS(newboundnores - trafoubs[cliquecovers[i][1]])) )
                  inftynores = TRUE;
               else
                  newboundnores -= trafoubs[cliquecovers[i][1]];/*lint --e{679}*/
            }

            /* determine maximum without index v and if x_v is nonzero (note that the array 'cliquecovers' is sorted by the values of trafoub in non-increasing order) */
            for (j = 0; j < cliquecoversizes[i]; ++j)
            {
               indcliq = cliquecovers[i][j];
               assert( 0 <= indcliq && indcliq < ntrafolinvars );

               nodecliq = varGetNodeSOS1(conshdlrdata, trafolinvars[indcliq]); /* possibly -1 if variable is not involved in an SOS1 constraint */
               assert( nodecliq < nsos1vars );

               if ( v != indcliq )
               {
                  /* if nodev or nodecliq are not a member of an SOS1 constraint or the variable corresponding to nodecliq is not implied to be zero if x_v != 0  */
                  if ( nodev < 0 || nodecliq < 0 || (! isConnectedSOS1(adjacencymatrix, NULL, nodev, nodecliq) && ! isImpliedZero(conflictgraph, implnodes, nodecliq) ) )
                  {
                     if ( SCIPisInfinity(scip, trafoubs[indcliq]) || SCIPisInfinity(scip, REALABS(newboundnonzero - trafoubs[indcliq])) )
                        ++ninftynonzero;
                     else
                        newboundnonzero -= trafoubs[indcliq];
                     break; /* break since we are only interested in the maximum upper bound among the variables in the clique cover;
                             * the variables in the clique cover form an SOS1 constraint, thus only one of them can be nonzero */
                  }
               }
            }
         }
         assert( ninftynonzero == 0 || inftynores );

         /* if computed upper bound is not infinity and variable is contained in linear constraint */
         if ( ninftynonzero == 0 && v < ntrafolinvars )
         {
            linval = trafolinvals[v];

            if ( SCIPisFeasZero(scip, linval) )
               continue;

            /* compute new bound */
            if ( SCIPisFeasPositive(scip, newboundnores) && ! inftynores )
               newbound = newboundnonzero;
            else
               newbound = MIN(0, newboundnonzero);
            newbound /= linval;

            /* check if new bound is tighter than the old one or problem is infeasible */
            if ( SCIPisFeasPositive(scip, linval) && SCIPisFeasLT(scip, lb, newbound) )
            {
               if ( SCIPisFeasLT(scip, ub, newbound) )
               {
                  *cutoff = TRUE;
                  break;
               }

               if ( SCIPvarIsIntegral(var) )
                  newbound = SCIPceil(scip, newbound);

               SCIP_CALL( SCIPtightenVarLb(scip, var, newbound, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );

               if ( tightened )
               {
                  SCIPdebugMsg(scip, "changed lower bound of variable %s from %f to %f \n", SCIPvarGetName(var), lb, newbound);
                  ++(*nchgbds);
               }
            }
            else if ( SCIPisFeasNegative(scip, linval) && SCIPisFeasGT(scip, ub, newbound) )
            {
               /* if assumption a_i * x_i != 0 was not correct */
               if ( SCIPisFeasGT(scip, SCIPvarGetLbLocal(var), newbound) )
               {
                  *cutoff = TRUE;
                  break;
               }

               if ( SCIPvarIsIntegral(var) )
                  newbound = SCIPfloor(scip, newbound);

               SCIP_CALL( SCIPtightenVarUb(scip, var, newbound, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );

               if ( tightened )
               {
                  SCIPdebugMsg(scip, "changed upper bound of variable %s from %f to %f \n", SCIPvarGetName(var), ub, newbound);
                  ++(*nchgbds);
               }
            }
         }

         /* update implication graph if possible */
         SCIP_CALL( updateImplicationGraphSOS1(scip, conshdlrdata, conflictgraph, adjacencymatrix, implgraph, implhash, implnodes, totalvars, cliquecovers, cliquecoversizes, varincover,
               trafolinvars, trafolinvals, ntrafolinvars, trafoubs, var, trafoubv, newboundnonzero, ninftynonzero, TRUE, nchgbds, &update, &infeasible) );
         if ( infeasible )
            *cutoff = TRUE;
         else if ( update )
            *implupdate = TRUE;
      }

      if ( *cutoff == TRUE )
      {
         /* free memory */
         SCIPfreeBufferArrayNull(scip, &varincover);
         for (j = ncliquecovers-1; j >= 0; --j)
            SCIPfreeBufferArrayNull(scip, &cliquecovers[j]);
         SCIPfreeBufferArrayNull(scip, &cliquecovers);
         SCIPfreeBufferArrayNull(scip, &cliquecoversizes);
         SCIPfreeBufferArrayNull(scip, &cliquecovervals);
         SCIPfreeBufferArrayNull(scip, &trafolinvals);
         SCIPfreeBufferArrayNull(scip, &trafolinvars);
         break;
      }


      /* try to tighten upper bounds */

      /* sort each cliquecover array in ascending order of the lower bounds of a_i * x_i; fill vector varincover */
      for (i = 0; i < ncliquecovers; ++i)
      {
         for (j = 0; j < cliquecoversizes[i]; ++j)
         {
            int ind = cliquecovers[i][j];

            varincover[ind] = i;
            cliquecovervals[j] = trafolbs[ind];
         }
         SCIPsortRealInt(cliquecovervals, cliquecovers[i], cliquecoversizes[i]);
      }

      /* for every variable that is in transformed constraint or every variable that is in conflict with some variable from trans. cons.:
         try upper bound tightening */
      for (v = 0; v < ntrafolinvars + nsos1linvars; ++v)
      {
         SCIP_Real newboundnonzero; /* new bound of a_v*x_v if we assume that x_v != 0 */
         SCIP_Real newboundnores;   /* new bound of a_v*x_v if there are no restrictions */
         SCIP_Real newbound;        /* resulting new bound of x_v */
         SCIP_VAR* var;
         SCIP_Real linval;
         SCIP_Real trafolbv;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Bool tightened;
         SCIP_Bool infeasible;
         SCIP_Bool inftynores = FALSE;
         SCIP_Bool update;
         int ninftynonzero = 0;
         int nodev;
         int w;

         if ( v < ntrafolinvars )
         {
            var = trafolinvars[v];
            trafolbv = trafolbs[v];
         }
         else
         {
            assert( v-ntrafolinvars >= 0 );
            var = sos1linvars[v-ntrafolinvars];/*lint !e679*/
            trafolbv = 0.0; /* since variable is not a member of linear constraint */
         }
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         if ( SCIPisInfinity(scip, traforhs) || SCIPisEQ(scip, lb, ub) )
            continue;

         newboundnonzero = traforhs;
         newboundnores = traforhs;
         nodev = varGetNodeSOS1(conshdlrdata, var); /* possibly -1 if var is not involved in an SOS1 constraint */
         assert( nodev < nsos1vars );

         /* determine incidence vector of implication variables (i.e., which SOS1 variables are nonzero if x_v is nonzero) */
         for (w = 0; w < nsos1vars; ++w)
            implnodes[w] = FALSE;
         SCIP_CALL( getSOS1Implications(scip, conshdlrdata, totalvars, implgraph, implhash, implnodes, (int) (size_t) SCIPhashmapGetImage(implhash, var)) );

         /* compute new bound */
         for (i = 0; i < ncliquecovers; ++i)
         {
            int indcliq;
            int nodecliq;

            assert( cliquecoversizes[i] > 0 );

            indcliq = cliquecovers[i][0];
            assert( 0 <= indcliq && indcliq < ntrafolinvars );

            /* determine minimum without index v (note that the array 'cliquecovers' is sorted by the values of trafolb in increasing order) */
            if ( v != indcliq )
            {
               /* if bound would be infinity */
               if ( SCIPisInfinity(scip, -trafolbs[indcliq]) || SCIPisInfinity(scip, REALABS(newboundnores - trafolbs[indcliq])) )
                  inftynores = TRUE;
               else
                  newboundnores -= trafolbs[indcliq];
            }
            else if ( cliquecoversizes[i] > 1 )
            {
               assert( 0 <= cliquecovers[i][1] && cliquecovers[i][1] < ntrafolinvars );
               if ( SCIPisInfinity(scip, -trafolbs[cliquecovers[i][1]]) || SCIPisInfinity(scip, REALABS(newboundnores - trafolbs[cliquecovers[i][1]])) )
                  inftynores = TRUE;
               else
                  newboundnores -= trafolbs[cliquecovers[i][1]]; /*lint --e{679}*/
            }

            /* determine minimum without index v and if x_v is nonzero (note that the array 'cliquecovers' is sorted by the values of trafolb in increasing order) */
            for (j = 0; j < cliquecoversizes[i]; ++j)
            {
               indcliq = cliquecovers[i][j];
               assert( 0 <= indcliq && indcliq < ntrafolinvars );

               nodecliq = varGetNodeSOS1(conshdlrdata, trafolinvars[indcliq]); /* possibly -1 if variable is not involved in an SOS1 constraint */
               assert( nodecliq < nsos1vars );

               if ( v != indcliq )
               {
                  /* if nodev or nodecliq are not a member of an SOS1 constraint or the variable corresponding to nodecliq is not implied to be zero if x_v != 0  */
                  if ( nodev < 0 || nodecliq < 0 || (! isConnectedSOS1(adjacencymatrix, NULL, nodev, nodecliq) && ! isImpliedZero(conflictgraph, implnodes, nodecliq) ) )
                  {
                     /* if bound would be infinity */
                     if ( SCIPisInfinity(scip, -trafolbs[indcliq]) || SCIPisInfinity(scip, REALABS(newboundnonzero - trafolbs[indcliq])) )
                        ++ninftynonzero;
                     else
                        newboundnonzero -= trafolbs[indcliq];
                     break; /* break since we are only interested in the minimum lower bound among the variables in the clique cover;
                             * the variables in the clique cover form an SOS1 constraint, thus only one of them can be nonzero */
                  }
               }
            }
         }
         assert( ninftynonzero == 0 || inftynores );


         /* if computed bound is not infinity and variable is contained in linear constraint */
         if ( ninftynonzero == 0 && v < ntrafolinvars )
         {
            linval = trafolinvals[v];

            if ( SCIPisFeasZero(scip, linval) )
               continue;

            /* compute new bound */
            if ( SCIPisFeasNegative(scip, newboundnores) && ! inftynores )
               newbound = newboundnonzero;
            else
               newbound = MAX(0, newboundnonzero);
            newbound /= linval;

            /* check if new bound is tighter than the old one or problem is infeasible */
            if ( SCIPisFeasPositive(scip, linval) && SCIPisFeasGT(scip, ub, newbound) )
            {
               /* if new upper bound is smaller than the lower bound, we are infeasible */
               if ( SCIPisFeasGT(scip, lb, newbound) )
               {
                  *cutoff = TRUE;
                  break;
               }

               if ( SCIPvarIsIntegral(var) )
                  newbound = SCIPfloor(scip, newbound);

               SCIP_CALL( SCIPtightenVarUb(scip, var, newbound, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );

               if ( tightened )
               {
                  SCIPdebugMsg(scip, "changed upper bound of variable %s from %f to %f \n", SCIPvarGetName(var), ub, newbound);
                  ++(*nchgbds);
               }
            }
            else if ( SCIPisFeasNegative(scip, linval) && SCIPisFeasLT(scip, lb, newbound) )
            {
               /* if assumption a_i * x_i != 0 was not correct */
               if ( SCIPisFeasLT(scip, ub, newbound) )
               {
                  *cutoff = TRUE;
                  break;
               }

               if ( SCIPvarIsIntegral(var) )
                  newbound = SCIPceil(scip, newbound);

               SCIP_CALL( SCIPtightenVarLb(scip, var, newbound, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );

               if ( tightened )
               {
                  SCIPdebugMsg(scip, "changed lower bound of variable %s from %f to %f \n", SCIPvarGetName(var), lb, newbound);
                  ++(*nchgbds);
               }
            }
         }

         /* update implication graph if possible */
         SCIP_CALL( updateImplicationGraphSOS1(scip, conshdlrdata, conflictgraph, adjacencymatrix, implgraph, implhash, implnodes, totalvars, cliquecovers, cliquecoversizes, varincover,
               trafolinvars, trafolinvals, ntrafolinvars, trafolbs, var, trafolbv, newboundnonzero, ninftynonzero, FALSE, nchgbds, &update, &infeasible) );
         if ( infeasible )
            *cutoff = TRUE;
         else if ( update )
            *implupdate = TRUE;
      }

      /* free memory */
      SCIPfreeBufferArrayNull(scip, &varincover);
      for (j = ncliquecovers-1; j >= 0; --j)
         SCIPfreeBufferArrayNull(scip, &cliquecovers[j]);
      SCIPfreeBufferArrayNull(scip, &cliquecovers);
      SCIPfreeBufferArrayNull(scip, &cliquecoversizes);
      SCIPfreeBufferArrayNull(scip, &cliquecovervals);
      SCIPfreeBufferArrayNull(scip, &trafolinvals);
      SCIPfreeBufferArrayNull(scip, &trafolinvars);

      if ( *cutoff == TRUE )
         break;
   } /* end for every linear constraint */

   /* free buffer arrays */
   SCIPfreeBufferArrayNull(scip, &sos1linvars);
   SCIPfreeBufferArrayNull(scip, &trafolbs);
   SCIPfreeBufferArrayNull(scip, &trafoubs);
   SCIPfreeBufferArrayNull(scip, &coveredvars);
   SCIPfreeBufferArrayNull(scip, &varindincons);
   SCIPfreeBufferArrayNull(scip, &implnodes);

   return SCIP_OKAY;
}


/** perform one presolving round for variables
 *
 *  We perform the following presolving steps:
 *  - Tighten the bounds of the variables
 *  - Update conflict graph based on bound implications of the variables
 */
static
SCIP_RETCODE presolRoundVarsSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_Bool**           adjacencymatrix,    /**< adjacencymatrix of conflict graph */
   int                   nsos1vars,          /**< number of SOS1 variables */
   int*                  nfixedvars,         /**< pointer to store number of fixed variables */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   int*                  naddconss,          /**< pointer to store number of addded constraints */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_DIGRAPH* implgraph;
   SCIP_HASHMAP* implhash;

   SCIP_Bool cutoff = FALSE;
   SCIP_Bool updateconfl;

   SCIP_VAR** totalvars;
   SCIP_VAR** probvars;
   int ntotalvars = 0;
   int nprobvars;
   int i;
   int j;

   /* determine totalvars (union of SOS1 and problem variables) */
   probvars = SCIPgetVars(scip);
   nprobvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPhashmapCreate(&implhash, SCIPblkmem(scip), nsos1vars + nprobvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &totalvars, nsos1vars + nprobvars) );

   for (i = 0; i < nsos1vars; ++i)
   {
      SCIP_VAR* var;
      var = SCIPnodeGetVarSOS1(conflictgraph, i);

      /* insert node number to hash map */
      assert( ! SCIPhashmapExists(implhash, var) );
      SCIP_CALL( SCIPhashmapInsert(implhash, var, (void*) (size_t) ntotalvars) );/*lint !e571*/
      assert( ntotalvars == (int) (size_t) SCIPhashmapGetImage(implhash, var) );
      totalvars[ntotalvars++] = var;
   }

   for (i = 0; i < nprobvars; ++i)
   {
      SCIP_VAR* var;
      var = probvars[i];

      /* insert node number to hash map if not existent */
      if ( ! SCIPhashmapExists(implhash, var) )
      {
         SCIP_CALL( SCIPhashmapInsert(implhash, var, (void*) (size_t) ntotalvars) );/*lint !e571*/
         assert( ntotalvars == (int) (size_t) SCIPhashmapGetImage(implhash, var) );
         totalvars[ntotalvars++] = var;
      }
   }

   /* create implication graph */
   SCIP_CALL( SCIPcreateDigraph(scip, &implgraph, ntotalvars) );

   /* try to tighten the lower and upper bounds of the variables */
   updateconfl = FALSE;
   for (j = 0; (j < conshdlrdata->maxtightenbds || conshdlrdata->maxtightenbds == -1 ) && ! cutoff; ++j)
   {
      SCIP_Bool implupdate;
      int nchgbdssave;

      nchgbdssave = *nchgbds;

      assert( ntotalvars > 0 );
      SCIP_CALL( tightenVarsBoundsSOS1(scip, conshdlrdata, conflictgraph, implgraph, implhash, adjacencymatrix, totalvars, ntotalvars, nsos1vars, nchgbds, &implupdate, &cutoff) );
      if ( *nchgbds > nchgbdssave )
      {
         *result = SCIP_SUCCESS;
         if ( implupdate )
            updateconfl = TRUE;
      }
      else if ( implupdate )
         updateconfl = TRUE;
      else
         break;
   }

   /* perform implication graph analysis */
   if ( updateconfl && conshdlrdata->perfimplanalysis && ! cutoff )
   {
      SCIP_Real* implubs;
      SCIP_Real* impllbs;
      SCIP_Bool* implnodes;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      int naddconsssave;
      int probingdepth;

      /* allocate buffer arrays */
      SCIP_CALL( SCIPallocBufferArray(scip, &implnodes, nsos1vars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &impllbs, ntotalvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &implubs, ntotalvars) );

      naddconsssave = *naddconss;
      for (i = 0; i < nsos1vars; ++i)
      {
	 /* initialize data for implication graph analysis */
	 infeasible = FALSE;
	 probingdepth = 0;
         for (j = 0; j < nsos1vars; ++j)
            implnodes[j] = FALSE;
	 for (j = 0; j < ntotalvars; ++j)
	 {
	    impllbs[j] = SCIPvarGetLbLocal(totalvars[j]);
	    implubs[j] = SCIPvarGetUbLocal(totalvars[j]);
	 }

	 /* try to update the conflict graph based on the information of the implication graph */
	 SCIP_CALL( performImplicationGraphAnalysis(scip, conshdlrdata, conflictgraph, totalvars, implgraph, implhash, adjacencymatrix, i, i, impllbs, implubs, implnodes, naddconss, &probingdepth, &infeasible) );

	 /* if the subproblem turned out to be infeasible then fix variable to zero */
	 if ( infeasible )
	 {
	    SCIP_CALL( SCIPfixVar(scip, totalvars[i], 0.0, &infeasible, &fixed) );

	    if ( fixed )
	    {
	       SCIPdebugMsg(scip, "fixed variable %s with lower bound %f and upper bound %f to zero\n",
				SCIPvarGetName(totalvars[i]), SCIPvarGetLbLocal(totalvars[i]), SCIPvarGetUbLocal(totalvars[i]));
	       ++(*nfixedvars);
	    }

	    if ( infeasible )
	       cutoff = TRUE;
	 }
      }

      if ( *naddconss > naddconsssave )
         *result = SCIP_SUCCESS;

      /* free buffer arrays */
      SCIPfreeBufferArrayNull(scip, &implubs);
      SCIPfreeBufferArrayNull(scip, &impllbs);
      SCIPfreeBufferArrayNull(scip, &implnodes);
   }

   /* if an infeasibility has been detected */
   if ( cutoff )
   {
      SCIPdebugMsg(scip, "cutoff \n");
      *result = SCIP_CUTOFF;
   }

   /* free memory */;
   for (j = ntotalvars-1; j >= 0; --j)
   {
      SCIP_SUCCDATA** succdatas;
      int nsucc;
      int s;

      succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, j);
      nsucc = SCIPdigraphGetNSuccessors(implgraph, j);

      for (s = nsucc-1; s >= 0; --s)
         SCIPfreeBlockMemory(scip, &succdatas[s]);/*lint !e866*/
   }
   SCIPdigraphFree(&implgraph);
   SCIPfreeBufferArrayNull(scip, &totalvars);
   SCIPhashmapFree(&implhash);

   return SCIP_OKAY;
}


/* ----------------------------- propagation -------------------------------------*/

/** propagate variables of SOS1 constraint */
static
SCIP_RETCODE propConsSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   int*                  ngen                /**< number of domain changes */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( cutoff != NULL );
   assert( ngen != NULL );

   *cutoff = FALSE;

   /* if more than one variable is fixed to be nonzero */
   if ( consdata->nfixednonzeros > 1 )
   {
      SCIPdebugMsg(scip, "the node is infeasible, more than 1 variable is fixed to be nonzero.\n");
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if exactly one variable is fixed to be nonzero */
   if ( consdata->nfixednonzeros == 1 )
   {
      SCIP_VAR** vars;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;
      SCIP_Bool success;
      SCIP_Bool allVarFixed;
      int firstFixedNonzero;
      int nvars;
      int j;

      firstFixedNonzero = -1;
      nvars = consdata->nvars;
      vars = consdata->vars;
      assert( vars != NULL );

      /* search nonzero variable - is needed for propinfo */
      for (j = 0; j < nvars; ++j)
      {
         if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(vars[j])) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(vars[j])) )
         {
            firstFixedNonzero = j;
            break;
         }
      }
      assert( firstFixedNonzero >= 0 );

      SCIPdebugMsg(scip, "variable <%s> is fixed nonzero, fixing other variables to 0.\n", SCIPvarGetName(vars[firstFixedNonzero]));

      /* fix variables before firstFixedNonzero to 0 */
      allVarFixed = TRUE;
      for (j = 0; j < firstFixedNonzero; ++j)
      {
         /* fix variable */
         SCIP_CALL( inferVariableZero(scip, vars[j], cons, firstFixedNonzero, &infeasible, &tightened, &success) );
         assert( ! infeasible );
         allVarFixed = allVarFixed && success;
         if ( tightened )
            ++(*ngen);
      }

      /* fix variables after firstFixedNonzero to 0 */
      for (j = firstFixedNonzero+1; j < nvars; ++j)
      {
         /* fix variable */
         SCIP_CALL( inferVariableZero(scip, vars[j], cons, firstFixedNonzero, &infeasible, &tightened, &success) );
         assert( ! infeasible ); /* there should be no variables after firstFixedNonzero that are fixed to be nonzero */
         allVarFixed = allVarFixed && success;
         if ( tightened )
            ++(*ngen);
      }

      /* reset constraint age counter */
      if ( *ngen > 0 )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }

      /* delete constraint locally */
      if ( allVarFixed )
      {
         assert( !SCIPconsIsModifiable(cons) );
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }
   }

   return SCIP_OKAY;
}


/** propagate a variable that is known to be nonzero */
static
SCIP_RETCODE propVariableNonzero(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph */
   SCIP_CONS*            cons,               /**< some arbitrary SOS1 constraint */
   int                   node,               /**< conflict graph node of variable that is known to be nonzero */
   SCIP_Bool             implprop,           /**< whether implication graph propagation shall be applied */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   int*                  ngen                /**< number of domain changes */
   )
{
   int inferinfo;
   int* succ;
   int nsucc;
   int s;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( cutoff != NULL );
   assert( ngen != NULL );
   assert( node >= 0 );

   *cutoff = FALSE;
   inferinfo = -node - 1;

   /* by assumption zero is outside the domain of variable */
   assert( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(SCIPnodeGetVarSOS1(conflictgraph, node))) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(SCIPnodeGetVarSOS1(conflictgraph, node))) );

   /* apply conflict graph propagation (fix all neighbors in the conflict graph to zero) */
   succ = SCIPdigraphGetSuccessors(conflictgraph, node);
   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, node);
   for (s = 0; s < nsucc; ++s)
   {
      SCIP_VAR* succvar;
      SCIP_Real lb;
      SCIP_Real ub;

      succvar = SCIPnodeGetVarSOS1(conflictgraph, succ[s]);
      lb = SCIPvarGetLbLocal(succvar);
      ub = SCIPvarGetUbLocal(succvar);

      if ( ! SCIPisFeasZero(scip, lb) || ! SCIPisFeasZero(scip, ub) )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Bool success;

         /* fix variable if it is not multi-aggregated */
         SCIP_CALL( inferVariableZero(scip, succvar, cons, inferinfo, &infeasible, &tightened, &success) );

         if ( infeasible )
         {
            /* variable cannot be nonzero */
            *cutoff = TRUE;
            return SCIP_OKAY;
         }
         if ( tightened )
            ++(*ngen);
         assert( success || SCIPvarGetStatus(succvar) == SCIP_VARSTATUS_MULTAGGR );
      }
   }


   /* apply implication graph propagation */
   if ( implprop && implgraph != NULL )
   {
      SCIP_SUCCDATA** succdatas;

#ifndef NDEBUG
      SCIP_NODEDATA* nodedbgdata;
      nodedbgdata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(implgraph, node);
      assert( SCIPvarCompare(nodedbgdata->var, SCIPnodeGetVarSOS1(conflictgraph, node)) == 0 );
#endif

      /* get successor datas */
      succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, node);

      if ( succdatas != NULL )
      {
         succ = SCIPdigraphGetSuccessors(implgraph, node);
         nsucc = SCIPdigraphGetNSuccessors(implgraph, node);
         for (s = 0; s < nsucc; ++s)
         {
            SCIP_SUCCDATA* succdata;
            SCIP_NODEDATA* nodedata;
            SCIP_VAR* var;

            nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(implgraph, succ[s]);
            assert( nodedata != NULL );
            succdata = succdatas[s];
            assert( succdata != NULL );
            var = nodedata->var;
            assert( var != NULL );

            /* tighten variable if it is not multi-aggregated */
            if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
            {
               /* check for lower bound implication */
               if ( SCIPisFeasLT(scip, SCIPvarGetLbLocal(var), succdata->lbimpl) )
               {
                  SCIP_Bool infeasible;
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPinferVarLbCons(scip, var, succdata->lbimpl, cons, inferinfo, FALSE, &infeasible, &tightened) );
                  if ( infeasible )
                  {
                     *cutoff = TRUE;
                     return SCIP_OKAY;
                  }
                  if ( tightened )
                     ++(*ngen);
               }

               /* check for upper bound implication */
               if ( SCIPisFeasGT(scip, SCIPvarGetUbLocal(var), succdata->ubimpl) )
               {
                  SCIP_Bool infeasible;
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPinferVarUbCons(scip, var, succdata->ubimpl, cons, inferinfo, FALSE, &infeasible, &tightened) );
                  if ( infeasible )
                  {
                     *cutoff = TRUE;
                     return SCIP_OKAY;
                  }
                  if ( tightened )
                     ++(*ngen);
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** initialize implication graph
 *
 *  @p j is successor of @p i if and only if \f$ x_i\not = 0 \Rightarrow x_j\not = 0\f$
 *
 *  @note By construction the implication graph is globally valid.
 */
static
SCIP_RETCODE initImplGraphSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   nsos1vars,          /**< number of SOS1 variables */
   int                   maxrounds,          /**< maximal number of propagation rounds for generating implications */
   int*                  nchgbds,            /**< pointer to store number of bound changes */
   SCIP_Bool*            cutoff,             /**< pointer to store  whether a cutoff occurred */
   SCIP_Bool*            success             /**< whether initialization was successful */
   )
{
   SCIP_HASHMAP* implhash = NULL;
   SCIP_Bool** adjacencymatrix = NULL;
   SCIP_Bool* implnodes = NULL;
   SCIP_VAR** implvars = NULL;
   SCIP_VAR** probvars;
   int nimplnodes;
   int nprobvars;
   int i;
   int j;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( conflictgraph != NULL );
   assert( conshdlrdata->implgraph == NULL );
   assert( conshdlrdata->nimplnodes == 0 );
   assert( cutoff != NULL );
   assert( nchgbds != NULL );

   *nchgbds = 0;
   *cutoff = FALSE;

   /* we do not create the adjacency matrix of the conflict graph if the number of SOS1 variables is larger than a predefined value */
   if ( conshdlrdata->maxsosadjacency != -1 && nsos1vars > conshdlrdata->maxsosadjacency )
   {
      *success = FALSE;
      SCIPdebugMsg(scip, "Implication graph was not created since number of SOS1 variables (%d) is larger than %d.\n", nsos1vars, conshdlrdata->maxsosadjacency);

      return SCIP_OKAY;
   }
   *success = TRUE;

   /* only add globally valid implications to implication graph */
   assert ( SCIPgetDepth(scip) == 0 );

   probvars = SCIPgetVars(scip);
   nprobvars = SCIPgetNVars(scip);
   nimplnodes = 0;

   /* create implication graph */
   SCIP_CALL( SCIPdigraphCreate(&conshdlrdata->implgraph, SCIPblkmem(scip), nsos1vars + nprobvars) );

   /* create hashmap */
   SCIP_CALL( SCIPhashmapCreate(&implhash, SCIPblkmem(scip), nsos1vars + nprobvars) );

   /* determine implvars (union of SOS1 and problem variables)
    * Note: For separation of implied bound cuts it is important that SOS1 variables are enumerated first
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &implvars, nsos1vars + nprobvars) );
   for (i = 0; i < nsos1vars; ++i)
   {
      SCIP_VAR* var;
      var = SCIPnodeGetVarSOS1(conflictgraph, i);

      /* insert node number to hash map */
      assert( ! SCIPhashmapExists(implhash, var) );
      SCIP_CALL( SCIPhashmapInsert(implhash, var, (void*) (size_t) nimplnodes) );/*lint !e571*/
      assert( nimplnodes == (int) (size_t) SCIPhashmapGetImage(implhash, var) );
      implvars[nimplnodes++] = var;
   }

   for (i = 0; i < nprobvars; ++i)
   {
      SCIP_VAR* var;
      var = probvars[i];

      /* insert node number to hash map if not existent */
      if ( ! SCIPhashmapExists(implhash, var) )
      {
         SCIP_CALL( SCIPhashmapInsert(implhash, var, (void*) (size_t) nimplnodes) );/*lint !e571*/
         assert( nimplnodes == (int) (size_t) SCIPhashmapGetImage(implhash, var) );
         implvars[nimplnodes++] = var;
      }
   }
   conshdlrdata->nimplnodes = nimplnodes;

   /* add variables to nodes of implication graph */
   for (i = 0; i < nimplnodes; ++i)
   {
      SCIP_NODEDATA* nodedata = NULL;

      /* create node data */
      SCIP_CALL( SCIPallocBlockMemory(scip, &nodedata) );
      nodedata->var = implvars[i];

      /* set node data */
      SCIPdigraphSetNodeData(conshdlrdata->implgraph, (void*) nodedata, i);
   }

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &implnodes, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &adjacencymatrix, nsos1vars) );

   for (i = 0; i < nsos1vars; ++i)
      SCIP_CALL( SCIPallocBufferArray(scip, &adjacencymatrix[i], i+1) ); /*lint !e866*/

   /* create adjacency matrix */
   for (i = 0; i < nsos1vars; ++i)
   {
      for (j = 0; j < i+1; ++j)
         adjacencymatrix[i][j] = 0;
   }

   for (i = 0; i < nsos1vars; ++i)
   {
      int* succ;
      int nsucc;
      succ = SCIPdigraphGetSuccessors(conflictgraph, i);
      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, i);

      for (j = 0; j < nsucc; ++j)
      {
         if ( i > succ[j] )
            adjacencymatrix[i][succ[j]] = 1;
      }
   }

   assert( SCIPgetDepth(scip) == 0 );

   /* compute SOS1 implications from linear constraints and tighten bounds of variables */
   for (j = 0; (j < maxrounds || maxrounds == -1 ); ++j)
   {
      SCIP_Bool implupdate;
      int nchgbdssave;

      nchgbdssave = *nchgbds;

      assert( nimplnodes > 0 );
      SCIP_CALL( tightenVarsBoundsSOS1(scip, conshdlrdata, conflictgraph, conshdlrdata->implgraph, implhash, adjacencymatrix, implvars, nimplnodes, nsos1vars, nchgbds, &implupdate, cutoff) );
      if ( *cutoff || ( ! implupdate && ! ( *nchgbds > nchgbdssave ) ) )
         break;
   }

   /* free memory */
   for (i = nsos1vars-1; i >= 0; --i)
      SCIPfreeBufferArrayNull(scip, &adjacencymatrix[i]);
   SCIPfreeBufferArrayNull(scip, &adjacencymatrix);
   SCIPfreeBufferArrayNull(scip, &implnodes);
   SCIPfreeBufferArrayNull(scip, &implvars);
   SCIPhashmapFree(&implhash);

#ifdef SCIP_DEBUG
   /* evaluate results */
   if ( cutoff )
   {
      SCIPdebugMsg(scip, "cutoff \n");
   }
   else if ( *nchgbds > 0 )
   {
      SCIPdebugMsg(scip, "found %d bound changes\n", *nchgbds);
   }
#endif

   assert( conshdlrdata->implgraph != NULL );

   return SCIP_OKAY;
}


/** deinitialize implication graph */
static
SCIP_RETCODE freeImplGraphSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   int j;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );

   /* free whole memory of implication graph */
   if ( conshdlrdata->implgraph == NULL )
   {
      assert( conshdlrdata->nimplnodes == 0 );
      return SCIP_OKAY;
   }

   /* free arc data */
   for (j = conshdlrdata->nimplnodes-1; j >= 0; --j)
   {
      SCIP_SUCCDATA** succdatas;
      int nsucc;
      int s;

      succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(conshdlrdata->implgraph, j);
      nsucc = SCIPdigraphGetNSuccessors(conshdlrdata->implgraph, j);

      for (s = nsucc-1; s >= 0; --s)
      {
         assert( succdatas[s] != NULL );
         SCIPfreeBlockMemory(scip, &succdatas[s]);/*lint !e866*/
      }
   }

   /* free node data */
   for (j = conshdlrdata->nimplnodes-1; j >= 0; --j)
   {
      SCIP_NODEDATA* nodedata;
      nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conshdlrdata->implgraph, j);
      assert( nodedata != NULL );
      SCIPfreeBlockMemory(scip, &nodedata);
      SCIPdigraphSetNodeData(conshdlrdata->implgraph, NULL, j);
   }

   /* free implication graph */
   SCIPdigraphFree(&conshdlrdata->implgraph);
   conshdlrdata->nimplnodes = 0;

   return SCIP_OKAY;
}


/* ----------------------------- branching -------------------------------------*/

/** get the vertices whose neighbor set covers a subset of the neighbor set of a given other vertex.
 *
 *  This function can be used to compute sets of variables to branch on.
 */
static
SCIP_RETCODE getCoverVertices(
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_Bool*            verticesarefixed,   /**< array that indicates which variables are currently fixed to zero */
   int                   vertex,             /**< vertex (-1 if not needed) */
   int*                  neightocover,       /**< neighbors of given vertex to be covered (or NULL if all neighbors shall be covered) */
   int                   nneightocover,      /**< number of entries of neightocover (or 0 if all neighbors shall be covered )*/
   int*                  coververtices,      /**< array to store the vertices whose neighbor set covers the neighbor set of the given vertex */
   int*                  ncoververtices      /**< pointer to store size of coververtices */
   )
{
   int* succ1;
   int nsucc1;
   int s;

   assert( conflictgraph != NULL );
   assert( verticesarefixed != NULL );
   assert( coververtices != NULL );
   assert( ncoververtices != NULL );

   *ncoververtices = 0;

   /* if all the neighbors shall be covered */
   if ( neightocover == NULL )
   {
      assert( nneightocover == 0 );
      nsucc1 = SCIPdigraphGetNSuccessors(conflictgraph, vertex);
      succ1 = SCIPdigraphGetSuccessors(conflictgraph, vertex);
   }
   else
   {
      nsucc1 = nneightocover;
      succ1 = neightocover;
   }

   /* determine all the successors of the first unfixed successor */
   for (s = 0; s < nsucc1; ++s)
   {
      int succvertex1 = succ1[s];

      if ( ! verticesarefixed[succvertex1] )
      {
         int succvertex2;
         int* succ2;
         int nsucc2;
         int j;

         nsucc2 = SCIPdigraphGetNSuccessors(conflictgraph, succvertex1);
         succ2 = SCIPdigraphGetSuccessors(conflictgraph, succvertex1);

         /* for the first unfixed vertex */
         if ( *ncoververtices == 0 )
         {
            for (j = 0; j < nsucc2; ++j)
            {
               succvertex2 = succ2[j];
               if ( ! verticesarefixed[succvertex2] )
                  coververtices[(*ncoververtices)++] = succvertex2;
            }
         }
         else
         {
            int vv = 0;
            int k = 0;
            int v;

            /* determine all the successors that are in the set "coververtices" */
            for (v = 0; v < *ncoververtices; ++v)
            {
               assert( vv <= v );
               for (j = k; j < nsucc2; ++j)
               {
                  succvertex2 = succ2[j];
                  if ( succvertex2 > coververtices[v] )
                  {
                     /* coververtices[v] does not appear in succ2 list, go to next vertex in coververtices */
                     k = j;
                     break;
                  }
                  else if ( succvertex2 == coververtices[v] )
                  {
                     /* vertices are equal, copy to free position vv */
                     coververtices[vv++] = succvertex2;
                     k = j + 1;
                     break;
                  }
               }
            }
            /* store new size of coververtices */
            *ncoververtices = vv;
         }
      }
   }

#ifdef SCIP_DEBUG
   /* check sorting */
   for (s = 0; s < *ncoververtices; ++s)
   {
      assert( *ncoververtices <= 1 || coververtices[*ncoververtices - 1] > coververtices[*ncoververtices - 2] );
   }
#endif

   return SCIP_OKAY;
}


/** get vertices of variables that will be fixed to zero for each node */
static
SCIP_RETCODE getBranchingVerticesSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< solution to be enforced (NULL for LP solution) */
   SCIP_Bool*            verticesarefixed,   /**< vector that indicates which variables are currently fixed to zero */
   SCIP_Bool             bipbranch,          /**< TRUE if bipartite branching method should be used */
   int                   branchvertex,       /**< branching vertex */
   int*                  fixingsnode1,       /**< vertices of variables that will be fixed to zero for the first node */
   int*                  nfixingsnode1,      /**< pointer to store number of fixed variables for the first node */
   int*                  fixingsnode2,       /**< vertices of variables that will be fixed to zero for the second node */
   int*                  nfixingsnode2       /**< pointer to store number of fixed variables for the second node */
   )
{
   SCIP_Bool takeallsucc; /* whether to set fixingsnode1 = neighbors of 'branchvertex' in the conflict graph */
   int* succ;
   int nsucc;
   int j;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( verticesarefixed != NULL );
   assert( ! verticesarefixed[branchvertex] );
   assert( fixingsnode1 != NULL );
   assert( fixingsnode2 != NULL );
   assert( nfixingsnode1 != NULL );
   assert( nfixingsnode2 != NULL );

   *nfixingsnode1 = 0;
   *nfixingsnode2 = 0;
   takeallsucc = TRUE;

   /* get successors and number of successors of branching vertex */
   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, branchvertex);
   succ = SCIPdigraphGetSuccessors(conflictgraph, branchvertex);

   /* if bipartite branching method is turned on */
   if ( bipbranch )
   {
      SCIP_Real solval;
      int cnt = 0;

      /* get all the neighbors of the variable with index 'branchvertex' whose solution value is nonzero */
      for (j = 0; j < nsucc; ++j)
      {
         if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, SCIPnodeGetVarSOS1(conflictgraph, succ[j]))) )
         {
            assert( ! verticesarefixed[succ[j]] );
            fixingsnode1[(*nfixingsnode1)++] = succ[j];
         }
      }

      /* if one of the sets fixingsnode1 or fixingsnode2 contains only one variable with a nonzero LP value we perform standard neighborhood branching */
      if ( *nfixingsnode1 > 0 )
      {
         /* get the vertices whose neighbor set cover the selected subset of the neighbors of the given branching vertex */
         SCIP_CALL( getCoverVertices(conflictgraph, verticesarefixed, branchvertex, fixingsnode1, *nfixingsnode1, fixingsnode2, nfixingsnode2) );

         /* determine the intersection of the neighbors of branchvertex with the intersection of all the neighbors of fixingsnode2 */
         SCIP_CALL( getCoverVertices(conflictgraph, verticesarefixed, branchvertex, fixingsnode2, *nfixingsnode2, fixingsnode1, nfixingsnode1) );

         for (j = 0; j < *nfixingsnode2; ++j)
         {
            solval = SCIPgetSolVal(scip, sol, SCIPnodeGetVarSOS1(conflictgraph, fixingsnode2[j]));
            if( ! SCIPisFeasZero(scip, solval) )
               ++cnt;
         }

         /* we decide whether to use all successors if one partition of complete bipartite subgraph has only one node */
         if ( cnt >= 2 )
         {
            cnt = 0;
            for (j = 0; j < *nfixingsnode1; ++j)
            {
               solval = SCIPgetSolVal(scip, sol, SCIPnodeGetVarSOS1(conflictgraph, fixingsnode1[j]));
               if( ! SCIPisFeasZero(scip, solval) )
                  ++cnt;
            }

            if ( cnt >= 2 )
               takeallsucc = FALSE;
         }
      }
   }

   if ( takeallsucc )
   {
      /* get all the unfixed neighbors of the branching vertex */
      *nfixingsnode1 = 0;
      for (j = 0; j < nsucc; ++j)
      {
         if ( ! verticesarefixed[succ[j]] )
            fixingsnode1[(*nfixingsnode1)++] = succ[j];
      }

      if ( bipbranch )
      {
         /* get the vertices whose neighbor set covers the neighbor set of a given branching vertex */
         SCIP_CALL( getCoverVertices(conflictgraph, verticesarefixed, branchvertex, fixingsnode1, *nfixingsnode1, fixingsnode2, nfixingsnode2) );
      }
      else
      {
         /* use neighborhood branching, i.e, for the second node only the branching vertex can be fixed */
         fixingsnode2[0] = branchvertex;
         *nfixingsnode2 = 1;
      }
   }

   return SCIP_OKAY;
}


/** gets branching priorities for SOS1 variables and applies 'most infeasible selection' rule to determine a vertex for the next branching decision */
static
SCIP_RETCODE getBranchingPrioritiesSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< solution to be enforced (NULL for LP solution) */
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_Bool*            verticesarefixed,   /**< vector that indicates which variables are currently fixed to zero */
   SCIP_Bool             bipbranch,          /**< TRUE if bipartite branching method should be used */
   int*                  fixingsnode1,       /**< vertices of variables that will be fixed to zero for the first node (size = nsos1vars) */
   int*                  fixingsnode2,       /**< vertices of variables that will be fixed to zero for the second node (size = nsos1vars) */
   SCIP_Real*            branchpriors,       /**< pointer to store branching priorities (size = nsos1vars) or NULL if not needed */
   int*                  vertexbestprior,    /**< pointer to store vertex with the best branching priority or NULL if not needed */
   SCIP_Bool*            relsolfeas          /**< pointer to store if LP relaxation solution is feasible */
   )
{
   SCIP_Real bestprior;
   int i;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( conflictgraph != NULL );
   assert( verticesarefixed != NULL );
   assert( fixingsnode1 != NULL );
   assert( fixingsnode2 != NULL );
   assert( relsolfeas != NULL );

   bestprior = -SCIPinfinity(scip);

   for (i = 0; i < nsos1vars; ++i)
   {
      SCIP_Real prior;
      SCIP_Real solval;
      SCIP_Real sum1;
      SCIP_Real sum2;
      int nfixingsnode1;
      int nfixingsnode2;
      int nsucc;
      int j;

      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, i);

      if ( nsucc == 0 || SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, SCIPnodeGetVarSOS1(conflictgraph, i))) || verticesarefixed[i] )
         prior = -SCIPinfinity(scip);
      else
      {
         SCIP_Bool iszero1 = TRUE;
         SCIP_Bool iszero2 = TRUE;

         /* get vertices of variables that will be fixed to zero for each strong branching execution */
         assert( ! verticesarefixed[i] );
         SCIP_CALL( getBranchingVerticesSOS1(scip, conflictgraph, sol, verticesarefixed, bipbranch, i, fixingsnode1, &nfixingsnode1, fixingsnode2, &nfixingsnode2) );

         sum1 = 0.0;
         for (j = 0; j < nfixingsnode1; ++j)
         {
            solval = SCIPgetSolVal(scip, sol, SCIPnodeGetVarSOS1(conflictgraph, fixingsnode1[j]));
            if ( ! SCIPisFeasZero(scip, solval) )
            {
               sum1 += REALABS( solval );
               iszero1 = FALSE;
            }
         }

         sum2 = 0.0;
         for (j = 0; j < nfixingsnode2; ++j)
         {
            solval = SCIPgetSolVal(scip, sol, SCIPnodeGetVarSOS1(conflictgraph, fixingsnode2[j]));
            if ( ! SCIPisFeasZero(scip, solval) )
            {
               sum2 += REALABS( solval );
               iszero2 = FALSE;
            }
         }

         if ( iszero1 || iszero2 )
            prior = -SCIPinfinity(scip);
         else
            prior = sum1 * sum2;
      }

      if ( branchpriors != NULL )
         branchpriors[i] = prior;
      if ( bestprior < prior )
      {
         bestprior = prior;

         if ( vertexbestprior != NULL )
            *vertexbestprior = i;
      }
   }

   if ( SCIPisInfinity(scip, -bestprior) )
      *relsolfeas = TRUE;
   else
      *relsolfeas = FALSE;

   return SCIP_OKAY;
}


/** performs strong branching with given domain fixings */
static
SCIP_RETCODE performStrongbranchSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int*                  fixingsexec,        /**< vertices of variables to be fixed to zero for this strong branching execution */
   int                   nfixingsexec,       /**< number of vertices of variables to be fixed to zero for this strong branching execution */
   int*                  fixingsop,          /**< vertices of variables to be fixed to zero for the opposite strong branching execution */
   int                   nfixingsop,         /**< number of vertices of variables to be fixed to zero for the opposite strong branching execution */
   int                   inititer,           /**< maximal number of LP iterations to perform */
   SCIP_Bool             fixnonzero,         /**< shall opposite variable (if positive in sign) fixed to the feasibility tolerance
                                              *   (only possible if nfixingsop = 1) */
   int*                  domainfixings,      /**< vertices that can be used to reduce the domain (should have size equal to number of variables) */
   int*                  ndomainfixings,     /**< pointer to store number of vertices that can be used to reduce the domain, could be filled by earlier calls */
   SCIP_Bool*            infeasible,         /**< pointer to store whether branch is infeasible */
   SCIP_Real*            objval,             /**< pointer to store objective value of LP with fixed variables (SCIP_INVALID if reddomain = TRUE or lperror = TRUE) */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error or a strange solution status occurred */
   )
{
   SCIP_LPSOLSTAT solstat;
   int i;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( fixingsexec != NULL );
   assert( nfixingsop > 0 );
   assert( fixingsop != NULL );
   assert( nfixingsop > 0 );
   assert( inititer >= -1 );
   assert( domainfixings != NULL );
   assert( ndomainfixings != NULL );
   assert( *ndomainfixings >= 0 );
   assert( infeasible != NULL );
   assert( objval != NULL );
   assert( lperror != NULL );

   *objval = SCIP_INVALID; /* for debugging */
   *lperror = FALSE;
   *infeasible = FALSE;

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* perform domain fixings */
   if ( fixnonzero && nfixingsop == 1 )
   {
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;

      var = SCIPnodeGetVarSOS1(conflictgraph, fixingsop[0]);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
      {
         if ( SCIPisZero(scip, lb) )
         {
            /* fix variable to some very small, but positive number or to 1.0 if variable is integral */
            if (SCIPvarIsIntegral(var) )
            {
               SCIP_CALL( SCIPchgVarLbProbing(scip, var, 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarLbProbing(scip, var, 1.5 * SCIPfeastol(scip)) );
            }
         }
         else if ( SCIPisZero(scip, ub) )
         {
            /* fix variable to some negative number with small absolute value or to -1.0 if variable is integral */
            if (SCIPvarIsIntegral(var) )
            {
               SCIP_CALL( SCIPchgVarUbProbing(scip, var, -1.0) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarUbProbing(scip, var, -1.5 * SCIPfeastol(scip)) );
            }
         }
      }
   }

   /* injects variable fixings into current probing node */
   for (i = 0; i < nfixingsexec && ! *infeasible; ++i)
   {
      SCIP_VAR* var;

      var = SCIPnodeGetVarSOS1(conflictgraph, fixingsexec[i]);
      if ( SCIPisFeasGT(scip, SCIPvarGetLbLocal(var), 0.0) || SCIPisFeasLT(scip, SCIPvarGetUbLocal(var), 0.0) )
         *infeasible = TRUE;
      else
      {
         SCIP_CALL( SCIPfixVarProbing(scip, var, 0.0) );
      }
   }

   /* apply domain propagation */
   if ( ! *infeasible )
   {
      SCIP_CALL( SCIPpropagateProbing(scip, 0, infeasible, NULL) );
   }

   if ( *infeasible )
      solstat = SCIP_LPSOLSTAT_INFEASIBLE;
   else
   {
      /* solve the probing LP */
      SCIP_CALL( SCIPsolveProbingLP(scip, inititer, lperror, NULL) );
      if ( *lperror )
      {
         SCIP_CALL( SCIPendProbing(scip) );
         return SCIP_OKAY;
      }

      /* get solution status */
      solstat = SCIPgetLPSolstat(scip);
   }

   /* if objective limit was reached, then the domain can be reduced */
   if ( solstat == SCIP_LPSOLSTAT_OBJLIMIT || solstat == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      *infeasible = TRUE;

      for (i = 0; i < nfixingsop; ++i)
         domainfixings[(*ndomainfixings)++] = fixingsop[i];
   }
   else if ( solstat == SCIP_LPSOLSTAT_OPTIMAL || solstat == SCIP_LPSOLSTAT_TIMELIMIT || solstat == SCIP_LPSOLSTAT_ITERLIMIT )
   {
      /* get objective value of probing LP */
      *objval = SCIPgetLPObjval(scip);
   }
   else
      *lperror = TRUE;

   /* end probing */
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}


/** apply strong branching to determine the vertex for the next branching decision */
static
SCIP_RETCODE getBranchingDecisionStrongbranchSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< solution to be enforced (NULL for LP solution) */
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_Real             lpobjval,           /**< current LP relaxation solution */
   SCIP_Bool             bipbranch,          /**< TRUE if bipartite branching method should be used */
   int                   nstrongrounds,      /**< number of strong branching rounds */
   SCIP_Bool*            verticesarefixed,   /**< vector that indicates which variables are currently fixed to zero */
   int*                  fixingsnode1,       /**< pointer to store vertices of variables that will be fixed to zero for the first node (size = nsos1vars) */
   int*                  fixingsnode2,       /**< pointer to store vertices of variables that will be fixed to zero for the second node (size = nsos1vars) */
   int*                  vertexbestprior,    /**< pointer to store vertex with the best strong branching priority */
   SCIP_Real*            bestobjval1,        /**< pointer to store LP objective for left child node of branching decision with best priority */
   SCIP_Real*            bestobjval2,        /**< pointer to store LP objective for right child node of branching decision with best priority */
   SCIP_RESULT*          result              /**< pointer to store result of strong branching */
   )
{
   SCIP_Real* branchpriors = NULL;
   int* indsos1vars = NULL;
   int* domainfixings = NULL;
   int ndomainfixings;
   int nfixingsnode1;
   int nfixingsnode2;

   SCIP_Bool relsolfeas;
   SCIP_Real bestscore;
   int lastscorechange;
   int maxfailures;

   SCIP_Longint nlpiterations;
   SCIP_Longint nlps;
   int inititer;
   int j;
   int i;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( conflictgraph != NULL );
   assert( verticesarefixed != NULL );
   assert( fixingsnode1 != NULL );
   assert( fixingsnode2 != NULL );
   assert( vertexbestprior != NULL );
   assert( result != NULL );

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchpriors, nsos1vars) );

   /* get branching priorities */
   SCIP_CALL( getBranchingPrioritiesSOS1(scip, conshdlrdata, conflictgraph, sol, nsos1vars, verticesarefixed,
         bipbranch, fixingsnode1, fixingsnode2, branchpriors, NULL, &relsolfeas) );

   /* if LP relaxation solution is feasible */
   if ( relsolfeas )
   {
      SCIPdebugMsg(scip, "all the SOS1 constraints are feasible.\n");
      *result = SCIP_FEASIBLE;

      /* free memory */
      SCIPfreeBufferArrayNull(scip, &branchpriors);

      return SCIP_OKAY;
   }

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &indsos1vars, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &domainfixings, nsos1vars) );

   /* sort branching priorities (descending order) */
   for (j = 0; j < nsos1vars; ++j)
      indsos1vars[j] = j;
   SCIPsortDownRealInt(branchpriors, indsos1vars, nsos1vars);

   /* determine the number of LP iterations to perform in each strong branch */
   nlpiterations =  SCIPgetNDualResolveLPIterations(scip);
   nlps = SCIPgetNDualResolveLPs(scip);
   if ( nlps == 0 )
   {
      nlpiterations = SCIPgetNNodeInitLPIterations(scip);
      nlps = SCIPgetNNodeInitLPs(scip);
      if ( nlps == 0 )
      {
         nlpiterations = 1000;
         nlps = 1;
      }
   }
   assert(nlps >= 1);

   /* compute number of LP iterations performed per strong branching iteration */
   if ( conshdlrdata->nstrongiter == -2 )
   {
      inititer = (int)(2*nlpiterations / nlps);
      inititer = (int)((SCIP_Real)inititer * (1.0 + 20.0/SCIPgetNNodes(scip)));
      inititer = MAX(inititer, 10);
      inititer = MIN(inititer, 500);
   }
   else
      inititer = conshdlrdata->nstrongiter;

   /* get current LP relaxation solution */
   lpobjval = SCIPgetLPObjval(scip);

   /* determine branching variable by strong branching or reduce domain */
   ndomainfixings = 0;
   lastscorechange = -1;
   *vertexbestprior = indsos1vars[0]; /* for the case that nstrongrounds = 0 */
   bestscore = -SCIPinfinity(scip);
   *bestobjval1 = -SCIPinfinity(scip);
   *bestobjval2 = -SCIPinfinity(scip);
   maxfailures = nstrongrounds;

   /* for each strong branching round */
   for (j = 0; j < nstrongrounds; ++j)
   {
      int testvertex;

      /* get branching vertex for the current strong branching iteration */
      testvertex = indsos1vars[j];

      /* if variable with index 'vertex' does not violate any complementarity in its neighborhood for the current LP relaxation solution */
      if ( SCIPisPositive(scip, branchpriors[j]) )
      {
         SCIP_Bool infeasible1;
         SCIP_Bool infeasible2;
         SCIP_Bool lperror;
         SCIP_Real objval1;
         SCIP_Real objval2;
         SCIP_Real score;

         /* get vertices of variables that will be fixed to zero for each strong branching execution */
         assert( ! verticesarefixed[testvertex] );
         SCIP_CALL( getBranchingVerticesSOS1(scip, conflictgraph, sol, verticesarefixed, bipbranch, testvertex, fixingsnode1, &nfixingsnode1, fixingsnode2, &nfixingsnode2) );

         /* get information for first strong branching execution */
         SCIP_CALL( performStrongbranchSOS1(scip, conflictgraph, fixingsnode1, nfixingsnode1, fixingsnode2, nfixingsnode2,
               inititer, conshdlrdata->fixnonzero, domainfixings, &ndomainfixings, &infeasible1, &objval1, &lperror) );
         if ( lperror )
            continue;

         /* get information for second strong branching execution */
         SCIP_CALL( performStrongbranchSOS1(scip, conflictgraph, fixingsnode2, nfixingsnode2, fixingsnode1, nfixingsnode1,
               inititer, FALSE, domainfixings, &ndomainfixings, &infeasible2, &objval2, &lperror) );
         if ( lperror )
            continue;

         /* if both subproblems are infeasible */
         if ( infeasible1 && infeasible2 )
         {
            SCIPdebugMsg(scip, "detected cutoff.\n");

            /* update result */
            *result = SCIP_CUTOFF;

            /* free memory */
            SCIPfreeBufferArrayNull(scip, &domainfixings);
            SCIPfreeBufferArrayNull(scip, &indsos1vars);
            SCIPfreeBufferArrayNull(scip, &branchpriors);

            return SCIP_OKAY;
         }
         else if ( ! infeasible1 && ! infeasible2 ) /* both subproblems are feasible */
         {
            /* if domain has not been reduced in this for-loop */
            if ( ndomainfixings == 0 )
            {
               score = MAX( REALABS( objval1 - lpobjval ), SCIPfeastol(scip) ) * MAX( REALABS( objval2 - lpobjval ), SCIPfeastol(scip) );/*lint !e666*/

               if ( SCIPisPositive(scip, score - bestscore) )
               {
                  bestscore = score;
                  *vertexbestprior = testvertex;
                  *bestobjval1 = objval1;
                  *bestobjval2 = objval2;

                  lastscorechange = j;
               }
               else if ( j - lastscorechange > maxfailures )
                  break;
            }
         }
      }
   }

   /* if variable fixings have been detected by probing, then reduce domain */
   if ( ndomainfixings > 0 )
   {
      SCIP_NODE* node = SCIPgetCurrentNode(scip);
      SCIP_Bool infeasible;

      for (i = 0; i < ndomainfixings; ++i)
      {
         SCIP_CALL( fixVariableZeroNode(scip, SCIPnodeGetVarSOS1(conflictgraph, domainfixings[i]), node, &infeasible) );
         assert( ! infeasible );
      }

      SCIPdebugMsg(scip, "found %d domain fixings.\n", ndomainfixings);

      /* update result */
      *result = SCIP_REDUCEDDOM;
   }

   /* free buffer arrays */
   SCIPfreeBufferArrayNull(scip, &domainfixings);
   SCIPfreeBufferArrayNull(scip, &indsos1vars);
   SCIPfreeBufferArrayNull(scip, &branchpriors);

   return SCIP_OKAY;
}


/** for two given vertices @p v1 and @p v2 search for a clique in the conflict graph that contains these vertices. From
 *  this clique, we create a bound constraint.
 */
static
SCIP_RETCODE getBoundConsFromVertices(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< solution to be enforced (NULL for LP solution) */
   int                   v1,                 /**< first vertex that shall be contained in bound constraint */
   int                   v2,                 /**< second vertex that shall be contained in bound constraint */
   SCIP_VAR*             boundvar,           /**< bound variable of @p v1 and @p v2 (or NULL if not existent) */
   SCIP_Bool             extend,             /**< should @p v1 and @p v2 be greedily extended to a clique of larger size */
   SCIP_CONS*            cons,               /**< bound constraint */
   SCIP_Real*            feas                /**< feasibility value of bound constraint */
   )
{
   SCIP_NODEDATA* nodedata;
   SCIP_Bool addv2 = TRUE;
   SCIP_Real solval;
   SCIP_VAR* var;
   SCIP_Real coef = 0.0;
   int nsucc;
   int s;

   int* extensions = NULL;
   int nextensions = 0;
   int nextensionsnew;
   int* succ;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( cons != NULL );
   assert( feas != NULL );

   *feas = 0.0;

   /* add index 'v1' to the clique */
   nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, v1);
   var = nodedata->var;
   assert( boundvar == NULL || SCIPvarCompare(boundvar, nodedata->ubboundvar) == 0 );
   solval = SCIPgetSolVal(scip, sol, var);

   /* if 'v1' and 'v2' have the same bound variable then the bound cut can be strengthened */
   if ( boundvar == NULL )
   {
      if ( SCIPisFeasPositive(scip, solval) )
      {
         SCIP_Real ub;
         ub = SCIPvarGetUbLocal(var);
         assert( SCIPisFeasPositive(scip, ub));

         if ( ! SCIPisInfinity(scip, ub) )
            coef = 1.0/ub;
      }
      else if ( SCIPisFeasNegative(scip, solval) )
      {
         SCIP_Real lb;
         lb = SCIPvarGetLbLocal(var);
         assert( SCIPisFeasNegative(scip, lb) );
         if ( ! SCIPisInfinity(scip, -lb) )
            coef = 1.0/lb;
      }
   }
   else if ( boundvar == nodedata->ubboundvar )
   {
      if ( SCIPisFeasPositive(scip, solval) )
      {
         SCIP_Real ub;

         ub = nodedata->ubboundcoef;
         assert( SCIPisFeasPositive(scip, ub) );
         if ( ! SCIPisInfinity(scip, ub) )
            coef = 1.0/ub;
      }
      else if ( SCIPisFeasNegative(scip, solval) )
      {
         SCIP_Real lb;

         lb = nodedata->lbboundcoef;
         assert( SCIPisFeasPositive(scip, lb) );
         if ( ! SCIPisInfinity(scip, lb) )
            coef = 1.0/lb;
      }
   }

   if ( ! SCIPisZero(scip, coef) )
   {
      *feas += coef * solval;
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, coef) );
   }

   /* if clique shall be greedily extended to a clique of larger size */
   if ( extend )
   {
      /* get successors */
      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, v1);
      succ = SCIPdigraphGetSuccessors(conflictgraph, v1);
      assert( nsucc > 0 );

      /* allocate buffer array */
      SCIP_CALL( SCIPallocBufferArray(scip, &extensions, nsucc) );

      /* get possible extensions for the clique cover */
      for (s = 0; s < nsucc; ++s)
         extensions[s] = succ[s];
      nextensions = nsucc;
   }
   else
      nextensions = 1;

   /* while there exist possible extensions for the clique cover */
   while ( nextensions > 0 )
   {
      SCIP_Real bestbigMval;
      SCIP_Real bigMval;
      int bestindex = -1;
      int ext;

      bestbigMval = -SCIPinfinity(scip);

      /* if v2 has not been added to clique already */
      if ( addv2 )
      {
         bestindex = v2;
         addv2 = FALSE;
      }
      else /* search for the extension with the largest absolute value of its LP relaxation solution value */
      {
         assert( extensions != NULL );
         for (s = 0; s < nextensions; ++s)
         {
            ext = extensions[s];
            bigMval = nodeGetSolvalBinaryBigMSOS1(scip, conflictgraph, sol, ext);
            if ( SCIPisFeasLT(scip, bestbigMval, bigMval) )
            {
               bestbigMval = bigMval;
               bestindex = ext;
            }
         }
      }
      assert( bestindex != -1 );

      /* add bestindex variable to the constraint */
      nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, bestindex);
      var = nodedata->var;
      solval = SCIPgetSolVal(scip, sol, var);
      coef = 0.0;
      if ( boundvar == NULL )
      {
         if ( SCIPisFeasPositive(scip, solval) )
         {
            SCIP_Real ub;
            ub = SCIPvarGetUbLocal(var);
            assert( SCIPisFeasPositive(scip, ub));

            if ( ! SCIPisInfinity(scip, ub) )
               coef = 1.0/ub;
         }
         else if ( SCIPisFeasNegative(scip, solval) )
         {
            SCIP_Real lb;
            lb = SCIPvarGetLbLocal(var);
            assert( SCIPisFeasNegative(scip, lb) );
            if ( ! SCIPisInfinity(scip, -lb) )
               coef = 1.0/lb;
         }
      }
      else if ( boundvar == nodedata->ubboundvar )
      {
         if ( SCIPisFeasPositive(scip, solval) )
         {
            SCIP_Real ub;

            ub = nodedata->ubboundcoef;
            assert( SCIPisFeasPositive(scip, ub) );
            if ( ! SCIPisInfinity(scip, ub) )
               coef = 1.0/ub;
         }
         else if ( SCIPisFeasNegative(scip, solval) )
         {
            SCIP_Real lb;

            lb = nodedata->lbboundcoef;
            assert( SCIPisFeasPositive(scip, lb) );
            if ( ! SCIPisInfinity(scip, -lb) )
               coef = 1.0/lb;
         }
      }
      if ( ! SCIPisZero(scip, coef) )
      {
         *feas += coef * solval;
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, coef) );
      }

      if ( extend )
      {
         assert( extensions != NULL );
         /* compute new 'extensions' array */
         nextensionsnew = 0;
         for (s = 0; s < nextensions; ++s)
         {
            if ( s != bestindex && isConnectedSOS1(NULL, conflictgraph, bestindex, extensions[s]) )
               extensions[nextensionsnew++] = extensions[s];
         }
         nextensions = nextensionsnew;
      }
      else
         nextensions = 0;
   }

   /* free buffer array */
   if ( extend )
      SCIPfreeBufferArray(scip, &extensions);

   /* subtract rhs of constraint from feasibility value or add bound variable if existent */
   if ( boundvar == NULL )
      *feas -= 1.0;
   else
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, boundvar, -1.0) );
      *feas -= SCIPgetSolVal(scip, sol, boundvar);
   }

   return SCIP_OKAY;
}


/** tries to add feasible complementarity constraints to a given child branching node.
 *
 *  @note In this function the conflict graph is updated to the conflict graph of the considered child branching node.
 */
static
SCIP_RETCODE addBranchingComplementaritiesSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_NODE*            node,               /**< branching node */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph of the current node */
   SCIP_DIGRAPH*         localconflicts,     /**< local conflicts (updates to local conflicts of child node) */
   SCIP_SOL*             sol,                /**< solution to be enforced (NULL for LP solution) */
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_Bool*            verticesarefixed,   /**< vector that indicates which variables are currently fixed to zerox */
   int*                  fixingsnode1,       /**< vertices of variables that will be fixed to zero for the branching node in the input of this function */
   int                   nfixingsnode1,      /**< number of entries of array nfixingsnode1 */
   int*                  fixingsnode2,       /**< vertices of variables that will be fixed to zero for the other branching node */
   int                   nfixingsnode2,      /**< number of entries of array nfixingsnode2 */
   int*                  naddedconss,        /**< pointer to store the number of added SOS1 constraints */
   SCIP_Bool             onlyviolsos1        /**< should only SOS1 constraints be added that are violated by the LP solution */
   )
{
   assert( scip != NULL );
   assert( node != NULL );
   assert( conshdlrdata != NULL );
   assert( conflictgraph != NULL );
   assert( verticesarefixed != NULL );
   assert( fixingsnode1 != NULL );
   assert( fixingsnode2 != NULL );
   assert( naddedconss != NULL );

   *naddedconss = 0;

   if ( nfixingsnode2 > 1 )
   {
      int* fixingsnode21; /* first partition of fixingsnode2 */
      int* fixingsnode22; /* second partition of fixingsnode2 */
      int nfixingsnode21;
      int nfixingsnode22;

      int* coverarray; /* vertices, not in fixingsnode1 that cover all the vertices in array fixingsnode22 */
      int ncoverarray;

      SCIP_Bool* mark;
      int* succarray;
      int nsuccarray;
      int* succ;
      int nsucc;

      int i;
      int s;

      /* allocate buffer arrays */
      SCIP_CALL( SCIPallocBufferArray(scip, &succarray, nsos1vars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &mark, nsos1vars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &fixingsnode21, nfixingsnode2) );
      SCIP_CALL( SCIPallocBufferArray(scip, &fixingsnode22, nfixingsnode2) );

      /* mark all the unfixed vertices with FALSE */
      for (i = 0; i < nsos1vars; ++i)
         mark[i] = (verticesarefixed[i]);

      /* mark all the vertices that are in the set fixingsnode1 */
      for (i = 0; i < nfixingsnode1; ++i)
      {
         assert( nfixingsnode1 <= 1 || (fixingsnode1[nfixingsnode1 - 1] > fixingsnode1[nfixingsnode1 - 2]) ); /* test: vertices are sorted */
         mark[fixingsnode1[i]] = TRUE;
      }

      /* mark all the vertices that are in the set fixingsnode2 */
      for (i = 0; i < nfixingsnode2; ++i)
      {
         assert( nfixingsnode2 <= 1 || (fixingsnode2[nfixingsnode2 - 1] > fixingsnode2[nfixingsnode2 - 2]) ); /* test: vertices are sorted */
         mark[fixingsnode2[i]] = TRUE;
      }

      /* compute the set of vertices that have a neighbor in the set fixingsnode2, but are not in the set fixingsnode1 or fixingsnode2 and are not already fixed */
      nsuccarray = 0;
      for (i = 0; i < nfixingsnode2; ++i)
      {
         nsucc = SCIPdigraphGetNSuccessors(conflictgraph, fixingsnode2[i]);
         succ = SCIPdigraphGetSuccessors(conflictgraph, fixingsnode2[i]);

         for (s = 0; s < nsucc; ++s)
         {
            int succnode = succ[s];

            if ( ! mark[succnode] )
            {
               mark[succnode] = TRUE;
               succarray[nsuccarray++] = succnode;
            }
         }
      }

      /* allocate buffer array */
      SCIP_CALL( SCIPallocBufferArray(scip, &coverarray, nsos1vars) );

      /* mark all the vertices with FALSE */
      for (i = 0; i < nsos1vars; ++i)
         mark[i] = FALSE;

      /* mark all the vertices that are in the set fixingsnode2 */
      for (i = 0; i < nfixingsnode2; ++i)
         mark[fixingsnode2[i]] = TRUE;

      /* for every node in succarray */
      for (i = 0; i < nsuccarray; ++i)
      {
         SCIP_Real solval1;
         SCIP_VAR* var1;
         int vertex1;
         int j;

         vertex1 = succarray[i];
         var1 = SCIPnodeGetVarSOS1(conflictgraph, vertex1);
         solval1 = SCIPgetSolVal(scip, sol, var1);

         /* we only add complementarity constraints if they are violated by the current LP solution */
         if ( ! onlyviolsos1 || ! SCIPisFeasZero(scip, solval1) )
         {
            /* compute first partition of fixingsnode2 that is the intersection of the neighbors of 'vertex1' with the set fixingsnode2 */
            nsucc = SCIPdigraphGetNSuccessors(conflictgraph, vertex1);
            succ = SCIPdigraphGetSuccessors(conflictgraph, vertex1);
            nfixingsnode21 = 0;

            for (s = 0; s < nsucc; ++s)
            {
               if ( mark[succ[s]] )
               {
                  fixingsnode21[nfixingsnode21++] = succ[s];
                  assert( nfixingsnode21 == 1 || (fixingsnode21[nfixingsnode21 - 1] > fixingsnode21[nfixingsnode21 - 2]) ); /* test: successor vertices are sorted */
               }
            }

            /* if variable can be fixed to zero */
            if ( nfixingsnode21 == nfixingsnode2 )
            {
               SCIP_Bool infeasible;

               SCIP_CALL( fixVariableZeroNode(scip, var1, node, &infeasible) );
               assert( ! infeasible );
               continue;
            }

            /* compute second partition of fixingsnode2 (that is fixingsnode2 \setminus fixingsnode21 ) */
            SCIP_CALL( SCIPcomputeArraysSetminus(fixingsnode2, nfixingsnode2, fixingsnode21, nfixingsnode21, fixingsnode22, &nfixingsnode22) );
            assert ( nfixingsnode22 + nfixingsnode21 == nfixingsnode2 );

            /* compute cover set (that are all the vertices not in fixingsnode1 and fixingsnode21, whose neighborhood covers all the vertices of fixingsnode22) */
            SCIP_CALL( getCoverVertices(conflictgraph, verticesarefixed, -1, fixingsnode22, nfixingsnode22, coverarray, &ncoverarray) );
            SCIP_CALL( SCIPcomputeArraysSetminus(coverarray, ncoverarray, fixingsnode1, nfixingsnode1, coverarray, &ncoverarray) );
            SCIP_CALL( SCIPcomputeArraysSetminus(coverarray, ncoverarray, fixingsnode21, nfixingsnode21, coverarray, &ncoverarray) );

            for (j = 0; j < ncoverarray; ++j)
            {
               int vertex2;

               vertex2 = coverarray[j];
               assert( vertex2 != vertex1 );

               /* prevent double enumeration */
               if ( vertex2 < vertex1 )
               {
                  SCIP_VAR* var2;
                  SCIP_Real solval2;

                  var2 = SCIPnodeGetVarSOS1(conflictgraph, vertex2);
                  solval2 = SCIPgetSolVal(scip, sol, var2);

                  if ( onlyviolsos1 && ( SCIPisFeasZero(scip, solval1) || SCIPisFeasZero(scip, solval2) ) )
                     continue;

                  if ( ! isConnectedSOS1(NULL, conflictgraph, vertex1, vertex2) )
                  {
                     char name[SCIP_MAXSTRLEN];
                     SCIP_CONS* conssos1 = NULL;
                     SCIP_Bool takebound = FALSE;
                     SCIP_Real feas;

                     SCIP_NODEDATA* nodedata;
                     SCIP_Real lbboundcoef1;
                     SCIP_Real lbboundcoef2;
                     SCIP_Real ubboundcoef1;
                     SCIP_Real ubboundcoef2;
                     SCIP_VAR* boundvar1;
                     SCIP_VAR* boundvar2;

                     /* get bound variables if available */
                     nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, vertex1);
                     assert( nodedata != NULL );
                     boundvar1 = nodedata->ubboundvar;
                     lbboundcoef1 = nodedata->lbboundcoef;
                     ubboundcoef1 = nodedata->ubboundcoef;
                     nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, vertex2);
                     assert( nodedata != NULL );
                     boundvar2 = nodedata->ubboundvar;
                     lbboundcoef2 = nodedata->lbboundcoef;
                     ubboundcoef2 = nodedata->ubboundcoef;

                     if ( boundvar1 != NULL && boundvar2 != NULL && SCIPvarCompare(boundvar1, boundvar2) == 0 )
                        takebound = TRUE;

                     /* add new arc to local conflicts in order to generate tighter bound inequalities */
                     if ( conshdlrdata->addextendedbds )
                     {
                        if ( localconflicts == NULL )
                        {
                           SCIP_CALL( SCIPdigraphCreate(&conshdlrdata->localconflicts, SCIPblkmem(scip), nsos1vars) );
                           localconflicts = conshdlrdata->localconflicts;
                        }
                        SCIP_CALL( SCIPdigraphAddArc(localconflicts, vertex1, vertex2, NULL) );
                        SCIP_CALL( SCIPdigraphAddArc(localconflicts, vertex2, vertex1, NULL) );
                        SCIP_CALL( SCIPdigraphAddArc(conflictgraph, vertex1, vertex2, NULL) );
                        SCIP_CALL( SCIPdigraphAddArc(conflictgraph, vertex2, vertex1, NULL) );

                        /* can sort successors in place - do not use arcdata */
                        SCIPsortInt(SCIPdigraphGetSuccessors(localconflicts, vertex1), SCIPdigraphGetNSuccessors(localconflicts, vertex1));
                        SCIPsortInt(SCIPdigraphGetSuccessors(localconflicts, vertex2), SCIPdigraphGetNSuccessors(localconflicts, vertex2));
                        SCIPsortInt(SCIPdigraphGetSuccessors(conflictgraph, vertex1), SCIPdigraphGetNSuccessors(conflictgraph, vertex1));
                        SCIPsortInt(SCIPdigraphGetSuccessors(conflictgraph, vertex2), SCIPdigraphGetNSuccessors(conflictgraph, vertex2));

                        /* mark conflictgraph as not local such that the new arcs are deleted after currents node processing */
                        conshdlrdata->isconflocal = TRUE;
                     }

                     /* measure feasibility of complementarity between var1 and var2 */
                     if ( ! takebound )
                     {
                        feas = -1.0;
                        if ( SCIPisFeasPositive(scip, solval1) )
                        {
                           assert( SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var1)));
                           if ( ! SCIPisInfinity(scip, SCIPvarGetUbLocal(var1)) )
                              feas += solval1/SCIPvarGetUbLocal(var1);
                        }
                        else if ( SCIPisFeasNegative(scip, solval1) )
                        {
                           assert( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var1)));
                           if ( ! SCIPisInfinity(scip, -SCIPvarGetLbLocal(var1)) )
                              feas += solval1/SCIPvarGetLbLocal(var1);
                        }

                        if ( SCIPisFeasPositive(scip, solval2) )
                        {
                           assert( SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var2)));
                           if ( ! SCIPisInfinity(scip, SCIPvarGetUbLocal(var2)) )
                              feas += solval2/SCIPvarGetUbLocal(var2);
                        }
                        else if ( SCIPisFeasNegative(scip, solval2) )
                        {
                           assert( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var2)));
                           if ( ! SCIPisInfinity(scip, -SCIPvarGetLbLocal(var2)) )
                              feas += solval2/SCIPvarGetLbLocal(var2);
                        }
                     }
                     else
                     {
                        feas = -SCIPgetSolVal(scip, sol, boundvar1);
                        if ( SCIPisFeasPositive(scip, solval1) )
                        {
                           assert( SCIPisFeasPositive(scip, ubboundcoef1));
                           if ( ! SCIPisInfinity(scip, ubboundcoef1) )
                              feas += solval1/ubboundcoef1;
                        }
                        else if ( SCIPisFeasNegative(scip, solval1) )
                        {
                           assert( SCIPisFeasPositive(scip, lbboundcoef1));
                           if ( ! SCIPisInfinity(scip, -lbboundcoef1) )
                              feas += solval1/lbboundcoef1;
                        }

                        if ( SCIPisFeasPositive(scip, solval2) )
                        {
                           assert( SCIPisFeasPositive(scip, ubboundcoef2));
                           if ( ! SCIPisInfinity(scip, ubboundcoef2) )
                              feas += solval2/ubboundcoef2;
                        }
                        else if ( SCIPisFeasNegative(scip, solval2) )
                        {
                           assert( SCIPisFeasPositive(scip, lbboundcoef2));
                           if ( ! SCIPisInfinity(scip, -lbboundcoef2) )
                              feas += solval2/lbboundcoef2;
                        }
                        assert( ! SCIPisFeasNegative(scip, solval2) );
                     }

                     if ( SCIPisGT(scip, feas, conshdlrdata->addcompsfeas) )
                     {
                        /* create SOS1 constraint */
                        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sos1_branchnode_%i_no_%i", SCIPnodeGetNumber(node), *naddedconss);
                        SCIP_CALL( SCIPcreateConsSOS1(scip, &conssos1, name, 0, NULL, NULL, TRUE, TRUE, TRUE, FALSE, TRUE,
                              TRUE, FALSE, FALSE, FALSE) );

                        /* add variables to SOS1 constraint */
                        SCIP_CALL( addVarSOS1(scip, conssos1, conshdlrdata, var1, 1.0) );
                        SCIP_CALL( addVarSOS1(scip, conssos1, conshdlrdata, var2, 2.0) );

                        /* add SOS1 constraint to the branching node */
                        SCIP_CALL( SCIPaddConsNode(scip, node, conssos1, NULL) );
                        ++(*naddedconss);

                        /* release constraint */
                        SCIP_CALL( SCIPreleaseCons(scip, &conssos1) );
                     }


                     /* add bound inequality*/
                     if ( ! SCIPisFeasZero(scip, solval1) && ! SCIPisFeasZero(scip, solval2) )
                     {
                        /* possibly create linear constraint of the form x_i/u_i + x_j/u_j <= t if a bound variable t with x_i <= u_i * t and x_j <= u_j * t exists.
                         * Otherwise try to create a constraint of the form x_i/u_i + x_j/u_j <= 1. Try the same for the lower bounds. */
                        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "boundcons_branchnode_%i_no_%i", SCIPnodeGetNumber(node), *naddedconss);
                        if ( takebound )
                        {
                           /* create constraint with right hand side = 0.0 */
                           SCIP_CALL( SCIPcreateConsLinear(scip, &conssos1, name, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, TRUE, FALSE, TRUE, FALSE, FALSE,
                                 TRUE, FALSE, FALSE, FALSE, FALSE) );

                           /* add variables */
                           SCIP_CALL( getBoundConsFromVertices(scip, conflictgraph, sol, vertex1, vertex2, boundvar1, conshdlrdata->addextendedbds, conssos1, &feas) );
                        }
                        else
                        {
                           /* create constraint with right hand side = 1.0 */
                           SCIP_CALL( SCIPcreateConsLinear(scip, &conssos1, name, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, FALSE, TRUE, FALSE, FALSE,
                                 TRUE, FALSE, FALSE, FALSE, FALSE) );

                           /* add variables */
                           SCIP_CALL( getBoundConsFromVertices(scip, conflictgraph, sol, vertex1, vertex2, NULL, conshdlrdata->addextendedbds, conssos1, &feas) );
                        }

                        /* add linear constraint to the branching node if usefull */
                        if ( SCIPisGT(scip, feas, conshdlrdata->addbdsfeas ) )
                        {
                           SCIP_CALL( SCIPaddConsNode(scip, node, conssos1, NULL) );
                           ++(*naddedconss);
                        }

                        /* release constraint */
                        SCIP_CALL( SCIPreleaseCons(scip, &conssos1) );
                     }

                     /* break if number of added constraints exceeds a predefined value */
                     if ( conshdlrdata->maxaddcomps >= 0 && *naddedconss > conshdlrdata->maxaddcomps )
                        break;
                  }
               }
            }
         }

         /* break if number of added constraints exceeds a predefined value */
         if ( conshdlrdata->maxaddcomps >= 0 && *naddedconss > conshdlrdata->maxaddcomps )
            break;
      }

      /* free buffer array */
      SCIPfreeBufferArray(scip, &coverarray);
      SCIPfreeBufferArray(scip, &fixingsnode22);
      SCIPfreeBufferArray(scip, &fixingsnode21);
      SCIPfreeBufferArray(scip, &mark);
      SCIPfreeBufferArray(scip, &succarray);
   }

   return SCIP_OKAY;
}


/** resets local conflict graph to the conflict graph of the root node */
static
SCIP_RETCODE resetConflictgraphSOS1(
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph of root node */
   SCIP_DIGRAPH*         localconflicts,     /**< local conflicts that should be removed from conflict graph */
   int                   nsos1vars           /**< number of SOS1 variables */
   )
{
   int j;

   for (j = 0; j < nsos1vars; ++j)
   {
      int nsuccloc;

      nsuccloc = SCIPdigraphGetNSuccessors(localconflicts, j);
      if ( nsuccloc > 0 )
      {
         int* succloc;
         int* succ;
         int nsucc;
         int k = 0;

         succloc = SCIPdigraphGetSuccessors(localconflicts, j);
         succ = SCIPdigraphGetSuccessors(conflictgraph, j);
         nsucc = SCIPdigraphGetNSuccessors(conflictgraph, j);

         /* reset number of successors */
         SCIP_CALL( SCIPcomputeArraysSetminus(succ, nsucc, succloc, nsuccloc, succ, &k) );
         SCIP_CALL( SCIPdigraphSetNSuccessors(conflictgraph, j, k) );
         SCIP_CALL( SCIPdigraphSetNSuccessors(localconflicts, j, 0) );
      }
   }

   return SCIP_OKAY;
}


/** Conflict graph enforcement method
 *
 *  The conflict graph can be enforced by different branching rules:
 *
 *  - Branch on the neighborhood of a single variable @p i, i.e., in one branch \f$x_i\f$ is fixed to zero and in the
 *    other its neighbors from the conflict graph.
 *
 *  - Branch on complete bipartite subgraphs of the conflict graph, i.e., in one branch fix the variables from the first
 *    bipartite partition and the variables from the second bipartite partition in the other.
 *
 *  - In addition to variable domain fixings, it is sometimes also possible to add new SOS1 constraints to the branching
 *    nodes. This results in a nonstatic conflict graph, which may change dynamically with every branching node.
 *
 *  We make use of different selection rules that define on which system of SOS1 variables to branch next:
 *
 *  - Most infeasible branching: Branch on the system of SOS1 variables with largest violation.
 *
 *  - Strong branching: Here, the LP-relaxation is partially solved for each branching decision among a candidate list.
 *    Then the decision with best progress is chosen.
 */
static
SCIP_RETCODE enforceConflictgraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   SCIP_SOL*             sol,                /**< solution to be enforced (NULL for LP solution) */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_DIGRAPH* conflictgraph;
   int nsos1vars;

   SCIP_Bool* verticesarefixed = NULL;
   int* fixingsnode1 = NULL;
   int* fixingsnode2 = NULL;
   int nfixingsnode1;
   int nfixingsnode2;

   SCIP_Real bestobjval1 = -SCIPinfinity(scip);
   SCIP_Real bestobjval2 = -SCIPinfinity(scip);
   SCIP_Real lpobjval = -SCIPinfinity(scip);

   SCIP_Bool infeasible;
   SCIP_Bool bipbranch = FALSE;
   int nstrongrounds;

   int branchvertex;
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   SCIP_Real nodeselest;
   SCIP_Real objest;

   int i;
   int j;
   int c;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing SOS1 conflict graph <%s>.\n", SCIPconshdlrGetName(conshdlr) );
   *result = SCIP_DIDNOTRUN;

   /* get number of SOS1 variables */
   nsos1vars = conshdlrdata->nsos1vars;

   /* get conflict graph */
   conflictgraph = conshdlrdata->conflictgraph;
   assert( ! conshdlrdata->isconflocal ); /* conflictgraph should be the one of the root node */

   /* check each constraint and update conflict graph if necessary */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_Bool cutoff;int ngen;

      cons = conss[c];
      assert( cons != NULL );
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      /* do nothing if there are not enough variables - this is usually eliminated by preprocessing */
      if ( consdata->nvars < 2 )
         continue;

      /* first perform propagation (it might happen that standard propagation is turned off) */
      ngen = 0;
      SCIP_CALL( propConsSOS1(scip, cons, consdata, &cutoff, &ngen) );
      SCIPdebugMsg(scip, "propagating <%s> in enforcing (cutoff: %u, domain reductions: %d).\n", SCIPconsGetName(cons), cutoff, ngen);
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
	 break;
      }
      if ( ngen > 0 )
      {
         *result = SCIP_REDUCEDDOM;
	 break;
      }
      assert( ngen == 0 );

      /* add local conflicts to conflict graph and save them in 'localconflicts' */
      if ( consdata->local )
      {
         SCIP_VAR** vars;
         int nvars;
         int indi;
         int indj;

         if ( conshdlrdata->localconflicts == NULL )
         {
            SCIP_CALL( SCIPdigraphCreate(&conshdlrdata->localconflicts, SCIPblkmem(scip), nsos1vars ) );
         }

         vars = consdata->vars;
         nvars = consdata->nvars;
         for (i = 0; i < nvars-1; ++i)
         {
            SCIP_VAR* var;

            var = vars[i];
            indi = varGetNodeSOS1(conshdlrdata, var);

            if( indi == -1 )
               return SCIP_INVALIDDATA;

            if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) )
            {
               for (j = i+1; j < nvars; ++j)
               {
                  var = vars[j];
                  indj = varGetNodeSOS1(conshdlrdata, var);

                  if( indj == -1 )
                     return SCIP_INVALIDDATA;

                  if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) )
                  {
                     if ( ! isConnectedSOS1(NULL, conflictgraph, indi, indj) )
                     {
                        SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraph, indi, indj, NULL) );
                        SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraph, indj, indi, NULL) );

                        SCIP_CALL( SCIPdigraphAddArcSafe(conshdlrdata->localconflicts, indi, indj, NULL) );
                        SCIP_CALL( SCIPdigraphAddArcSafe(conshdlrdata->localconflicts, indj, indi, NULL) );

                        conshdlrdata->isconflocal = TRUE;
                     }
                  }
               }
            }
         }
      }
   }

   /* sort successor list of conflict graph if necessary */
   if ( conshdlrdata->isconflocal )
   {
      for (j = 0; j < nsos1vars; ++j)
      {
         int nsuccloc;

         nsuccloc = SCIPdigraphGetNSuccessors(conshdlrdata->localconflicts, j);
         if ( nsuccloc > 0 )
         {
            SCIPsortInt(SCIPdigraphGetSuccessors(conflictgraph, j), SCIPdigraphGetNSuccessors(conflictgraph, j));
            SCIPsortInt(SCIPdigraphGetSuccessors(conshdlrdata->localconflicts, j), nsuccloc);
         }
      }
   }

   if ( *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM )
   {
      /* remove local conflicts from conflict graph */
      if ( conshdlrdata->isconflocal )
      {
	 SCIP_CALL( resetConflictgraphSOS1(conflictgraph, conshdlrdata->localconflicts, nsos1vars) );
	 conshdlrdata->isconflocal = FALSE;
      }
      return SCIP_OKAY;
   }


   /* detect fixed variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &verticesarefixed, nsos1vars) );
   for (j = 0; j < nsos1vars; ++j)
   {
      SCIP_VAR* var;
      SCIP_Real ub;
      SCIP_Real lb;

      var = SCIPnodeGetVarSOS1(conflictgraph, j);
      ub = SCIPvarGetUbLocal(var);
      lb = SCIPvarGetLbLocal(var);
      if ( SCIPisFeasZero(scip, ub) && SCIPisFeasZero(scip, lb) )
         verticesarefixed[j] = TRUE;
      else
         verticesarefixed[j] = FALSE;
   }

   /* should bipartite branching be used? */
   if ( conshdlrdata->branchingrule == 'b' )
      bipbranch = TRUE;

   /* determine number of strong branching iterations */
   if ( conshdlrdata->nstrongrounds >= 0 )
      nstrongrounds = MIN(conshdlrdata->nstrongrounds, nsos1vars);
   else
   {
      /* determine number depending on depth, based on heuristical considerations */
      if ( SCIPgetDepth(scip) <= 10 )
         nstrongrounds = MAX(10, (int)SCIPfloor(scip, pow(log((SCIP_Real)nsos1vars), 1.0)));/*lint !e666*/
      else if ( SCIPgetDepth(scip) <= 20 )
         nstrongrounds = MAX(5, (int)SCIPfloor(scip, pow(log((SCIP_Real)nsos1vars), 0.7)));/*lint !e666*/
      else
         nstrongrounds = 0;
      nstrongrounds = MIN(nsos1vars, nstrongrounds);
   }


   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &fixingsnode1, nsos1vars) );
   if ( bipbranch )
      SCIP_CALL( SCIPallocBufferArray(scip, &fixingsnode2, nsos1vars) );
   else
      SCIP_CALL( SCIPallocBufferArray(scip, &fixingsnode2, 1) );


   /* if strongbranching is turned off: use most infeasible branching */
   if ( nstrongrounds == 0 )
   {
      SCIP_Bool relsolfeas;

      /* get branching vertex using most infeasible branching */
      SCIP_CALL( getBranchingPrioritiesSOS1(scip, conshdlrdata, conflictgraph, sol, nsos1vars, verticesarefixed, bipbranch, fixingsnode1, fixingsnode2, NULL, &branchvertex, &relsolfeas) );

      /* if LP relaxation solution is feasible */
      if ( relsolfeas )
      {
         SCIPdebugMsg(scip, "all the SOS1 constraints are feasible.\n");

         /* update result */
         *result = SCIP_FEASIBLE;

         /* remove local conflicts from conflict graph */
         if ( conshdlrdata->isconflocal )
         {
            SCIP_CALL( resetConflictgraphSOS1(conflictgraph, conshdlrdata->localconflicts, nsos1vars) );
            conshdlrdata->isconflocal = FALSE;
         }

         /* free memory */
         SCIPfreeBufferArrayNull(scip, &fixingsnode2);
         SCIPfreeBufferArrayNull(scip, &fixingsnode1);
         SCIPfreeBufferArrayNull(scip, &verticesarefixed);

         return SCIP_OKAY;
      }
   }
   else
   {
      /* get branching vertex using strong branching */
      SCIP_CALL( getBranchingDecisionStrongbranchSOS1(scip, conshdlrdata, conflictgraph, sol, nsos1vars, lpobjval, bipbranch, nstrongrounds, verticesarefixed,
            fixingsnode1, fixingsnode2, &branchvertex, &bestobjval1, &bestobjval2, result) );

      if ( *result == SCIP_CUTOFF || *result == SCIP_FEASIBLE || *result == SCIP_REDUCEDDOM )
      {
         /* remove local conflicts from conflict graph */
         if ( conshdlrdata->isconflocal )
         {
            SCIP_CALL( resetConflictgraphSOS1(conflictgraph, conshdlrdata->localconflicts, nsos1vars) );
            conshdlrdata->isconflocal = FALSE;
         }

         /* free memory */
         SCIPfreeBufferArrayNull(scip, &fixingsnode2);
         SCIPfreeBufferArrayNull(scip, &fixingsnode1);
         SCIPfreeBufferArrayNull(scip, &verticesarefixed);

         return SCIP_OKAY;
      }
   }

   /* if we shouldleave branching decision to branching rules */
   if ( ! conshdlrdata->branchsos )
   {
      /* remove local conflicts from conflict graph */
      if ( conshdlrdata->isconflocal )
      {
	 SCIP_CALL( resetConflictgraphSOS1(conflictgraph, conshdlrdata->localconflicts, nsos1vars) );
	 conshdlrdata->isconflocal = FALSE;
      }

      if ( SCIPvarIsBinary(SCIPnodeGetVarSOS1(conflictgraph, branchvertex)) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
      else
      {
         SCIPerrorMessage("Incompatible parameter setting: branchsos can only be set to false if all SOS1 variables are binary.\n");
         return SCIP_PARAMETERWRONGVAL;
      }
   }

   /* create branching nodes */

   /* get vertices of variables that will be fixed to zero for each node */
   assert( branchvertex >= 0 && branchvertex < nsos1vars );
   assert( ! verticesarefixed[branchvertex] );
   SCIP_CALL( getBranchingVerticesSOS1(scip, conflictgraph, sol, verticesarefixed, bipbranch, branchvertex, fixingsnode1, &nfixingsnode1, fixingsnode2, &nfixingsnode2) );

   /* calculate node selection and objective estimate for node 1 */
   nodeselest = 0.0;
   objest = 0.0;
   for (j = 0; j < nfixingsnode1; ++j)
   {
      SCIP_VAR* var;

      var = SCIPnodeGetVarSOS1(conflictgraph, fixingsnode1[j]);
      nodeselest += SCIPcalcNodeselPriority(scip, var, SCIP_BRANCHDIR_DOWNWARDS, 0.0);
      objest += SCIPcalcChildEstimate(scip, var, 0.0);
   }
   /* take the average of the individual estimates */
   objest = objest/((SCIP_Real) nfixingsnode1);

   /* create node 1 */
   SCIP_CALL( SCIPcreateChild(scip, &node1, nodeselest, objest) );

   /* fix variables for the first node */
   if ( conshdlrdata->fixnonzero && nfixingsnode2 == 1 )
   {
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;

      var = SCIPnodeGetVarSOS1(conflictgraph, fixingsnode2[0]);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
      {
         if ( SCIPisZero(scip, lb) )
         {
            /* fix variable to some very small, but positive number or to 1.0 if variable is integral */
            if (SCIPvarIsIntegral(var) )
            {
               SCIP_CALL( SCIPchgVarLbNode(scip, node1, var, 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarLbNode(scip, node1, var, 1.5 * SCIPfeastol(scip)) );
            }
         }
         else if ( SCIPisZero(scip, ub) )
         {
            if (SCIPvarIsIntegral(var) )
            {
               /* fix variable to some negative number with small absolute value to -1.0 if variable is integral */
               SCIP_CALL( SCIPchgVarUbNode(scip, node1, var, -1.0) );
            }
            else
            {
               /* fix variable to some negative number with small absolute value to -1.0 if variable is integral */
               SCIP_CALL( SCIPchgVarUbNode(scip, node1, var, -1.5 * SCIPfeastol(scip)) );
            }
         }
      }
   }
   for (j = 0; j < nfixingsnode1; ++j)
   {
      /* fix variable to zero */
      SCIP_CALL( fixVariableZeroNode(scip, SCIPnodeGetVarSOS1(conflictgraph, fixingsnode1[j]), node1, &infeasible) );
      assert( ! infeasible );
   }

   /* calculate node selection and objective estimate for node 2 */
   nodeselest = 0.0;
   objest = 0.0;
   for (j = 0; j < nfixingsnode2; ++j)
   {
      SCIP_VAR* var;

      var = SCIPnodeGetVarSOS1(conflictgraph, fixingsnode2[j]);
      nodeselest += SCIPcalcNodeselPriority(scip, var, SCIP_BRANCHDIR_DOWNWARDS, 0.0);
      objest += SCIPcalcChildEstimate(scip, var, 0.0);
   }

   /* take the average of the individual estimates */
   objest = objest/((SCIP_Real) nfixingsnode2);

   /* create node 2 */
   SCIP_CALL( SCIPcreateChild(scip, &node2, nodeselest, objest) );

   /* fix variables to zero */
   for (j = 0; j < nfixingsnode2; ++j)
   {
      SCIP_CALL( fixVariableZeroNode(scip, SCIPnodeGetVarSOS1(conflictgraph, fixingsnode2[j]), node2, &infeasible) );
      assert( ! infeasible );
   }


   /* add complementarity constraints to the branching nodes */
   if ( conshdlrdata->addcomps && ( conshdlrdata->addcompsdepth == -1 || conshdlrdata->addcompsdepth >= SCIPgetDepth(scip) ) )
   {
      int naddedconss;

      assert( ! conshdlrdata->fixnonzero );

      /* add complementarity constraints to the left branching node */
      SCIP_CALL( addBranchingComplementaritiesSOS1(scip, node1, conshdlrdata, conflictgraph, conshdlrdata->localconflicts, sol,
               nsos1vars, verticesarefixed, fixingsnode1, nfixingsnode1, fixingsnode2, nfixingsnode2, &naddedconss, TRUE) );

      if ( naddedconss == 0 )
      {
         /* add complementarity constraints to the right branching node */
         SCIP_CALL( addBranchingComplementaritiesSOS1(scip, node2, conshdlrdata, conflictgraph, conshdlrdata->localconflicts, sol,
               nsos1vars, verticesarefixed, fixingsnode2, nfixingsnode2, fixingsnode1, nfixingsnode1, &naddedconss, TRUE) );
      }
   }

   /* sets node's lower bound to the best known value */
   if ( nstrongrounds > 0 )
   {
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, node1, MAX(lpobjval, bestobjval1) ) );
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, node2, MAX(lpobjval, bestobjval2) ) );
   }

   /* remove local conflicts from conflict graph */
   if ( conshdlrdata->isconflocal )
   {
      SCIP_CALL( resetConflictgraphSOS1(conflictgraph, conshdlrdata->localconflicts, nsos1vars) );
      conshdlrdata->isconflocal = FALSE;
   }

   /* free buffer arrays */
   SCIPfreeBufferArrayNull(scip, &fixingsnode2);
   SCIPfreeBufferArrayNull(scip, &fixingsnode1);
   SCIPfreeBufferArrayNull(scip, &verticesarefixed );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** SOS1 branching enforcement method
 *
 *  We check whether the current solution is feasible, i.e., contains at most one nonzero
 *  variable. If not, we branch along the lines indicated by Beale and Tomlin:
 *
 *  We first compute \f$W = \sum_{j=1}^n |x_i|\f$ and \f$w = \sum_{j=1}^n j\, |x_i|\f$. Then we
 *  search for the index \f$k\f$ that satisfies
 *  \f[
 *        k \leq \frac{w}{W} < k+1.
 *  \f]
 *  The branches are then
 *  \f[
 *        x_1 = 0, \ldots, x_k = 0 \qquad \mbox{and}\qquad x_{k+1} = 0, \ldots, x_n = 0.
 *  \f]
 *
 *  If the constraint contains two variables, the branching of course simplifies.
 *
 *  Depending on the parameters (@c branchnonzeros, @c branchweight) there are three ways to choose
 *  the branching constraint.
 *
 *  <TABLE>
 *  <TR><TD>@c branchnonzeros</TD><TD>@c branchweight</TD><TD>constraint chosen</TD></TR>
 *  <TR><TD>@c true          </TD><TD> ?             </TD><TD>most number of nonzeros</TD></TR>
 *  <TR><TD>@c false         </TD><TD> @c true       </TD><TD>maximal weight corresponding to nonzero variable</TD></TR>
 *  <TR><TD>@c false         </TD><TD> @c true       </TD><TD>largest sum of variable values</TD></TR>
 *  </TABLE>
 *
 *  @c branchnonzeros = @c false, @c branchweight = @c true allows the user to specify an order for
 *  the branching importance of the constraints (setting the weights accordingly).
 *
 *  Constraint branching can also be turned off using parameter @c branchsos.
 */
static
SCIP_RETCODE enforceConssSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< indicator constraints */
   SCIP_SOL*             sol,                /**< solution to be enforced (NULL for LP solution) */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   SCIP_CONS* branchCons;
   SCIP_Real maxWeight;
   SCIP_VAR** vars;
   int nvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   maxWeight = -SCIP_REAL_MAX;
   branchCons = NULL;

   SCIPdebugMsg(scip, "Enforcing SOS1 constraints <%s>.\n", SCIPconshdlrGetName(conshdlr) );
   *result = SCIP_FEASIBLE;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons;
      SCIP_Bool cutoff;
      SCIP_Real weight;
      int ngen;
      int cnt;
      int j;

      cons = conss[c];
      assert( cons != NULL );
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      ngen = 0;
      cnt = 0;
      nvars = consdata->nvars;
      vars = consdata->vars;

      /* do nothing if there are not enough variables - this is usually eliminated by preprocessing */
      if ( nvars < 2 )
         continue;

      /* first perform propagation (it might happen that standard propagation is turned off) */
      SCIP_CALL( propConsSOS1(scip, cons, consdata, &cutoff, &ngen) );
      SCIPdebugMsg(scip, "propagating <%s> in enforcing (cutoff: %u, domain reductions: %d).\n", SCIPconsGetName(cons), cutoff, ngen);
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if ( ngen > 0 )
      {
         *result = SCIP_REDUCEDDOM;
         return SCIP_OKAY;
      }
      assert( ngen == 0 );

      /* check constraint */
      weight = 0.0;
      for (j = 0; j < nvars; ++j)
      {
         SCIP_Real val = REALABS(SCIPgetSolVal(scip, sol, vars[j]));

         if ( ! SCIPisFeasZero(scip, val) )
         {
            if ( conshdlrdata->branchnonzeros )
               weight += 1.0;
            else
            {
               if ( conshdlrdata->branchweight )
               {
                  /* choose maximum nonzero-variable weight */
                  if ( consdata->weights[j] > weight )
                     weight = consdata->weights[j];
               }
               else
                  weight += val;
            }
            ++cnt;
         }
      }
      /* if constraint is violated */
      if ( cnt > 1 && weight > maxWeight )
      {
         maxWeight = weight;
         branchCons = cons;
      }
   }

   /* if all constraints are feasible */
   if ( branchCons == NULL )
   {
      SCIPdebugMsg(scip, "All SOS1 constraints are feasible.\n");
      return SCIP_OKAY;
   }

   /* if we should leave branching decision to branching rules */
   if ( ! conshdlrdata->branchsos )
   {
      int j;

      consdata = SCIPconsGetData(branchCons);
      for (j = 0; j < consdata->nvars; ++j)
      {
         if ( ! SCIPvarIsBinary(consdata->vars[j]) )
            break;
      }

      if ( j == consdata->nvars )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
      else
      {
         SCIPerrorMessage("Incompatible parameter setting: branchsos can only be set to false if all SOS1 variables are binary.\n");
         return SCIP_PARAMETERWRONGVAL;
      }
   }

   /* otherwise create branches */
   SCIPdebugMsg(scip, "Branching on constraint <%s> (weight: %f).\n", SCIPconsGetName(branchCons), maxWeight);
   consdata = SCIPconsGetData(branchCons);
   assert( consdata != NULL );
   nvars = consdata->nvars;
   vars = consdata->vars;

   if ( nvars == 2 )
   {
      SCIP_Bool infeasible;

      /* constraint is infeasible: */
      assert( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, vars[0])) && ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, vars[1])) );

      /* create branches */
      SCIPdebugMsg(scip, "Creating two branches.\n");

      SCIP_CALL( SCIPcreateChild(scip, &node1, SCIPcalcNodeselPriority(scip, vars[0], SCIP_BRANCHDIR_DOWNWARDS, 0.0), SCIPcalcChildEstimate(scip, vars[0], 0.0) ) );
      SCIP_CALL( fixVariableZeroNode(scip, vars[0], node1, &infeasible) );
      assert( ! infeasible );

      SCIP_CALL( SCIPcreateChild(scip, &node2, SCIPcalcNodeselPriority(scip, vars[1], SCIP_BRANCHDIR_DOWNWARDS, 0.0), SCIPcalcChildEstimate(scip, vars[1], 0.0) ) );
      SCIP_CALL( fixVariableZeroNode(scip, vars[1], node2, &infeasible) );
      assert( ! infeasible );
   }
   else
   {
      SCIP_Bool infeasible;
      SCIP_Real weight1;
      SCIP_Real weight2;
      SCIP_Real nodeselest;
      SCIP_Real objest;
      SCIP_Real w;
      int j;
      int ind;
      int cnt;

      cnt = 0;

      weight1 = 0.0;
      weight2 = 0.0;

      /* compute weight */
      for (j = 0; j < nvars; ++j)
      {
         SCIP_Real val = REALABS(SCIPgetSolVal(scip, sol, vars[j]));
         weight1 += val * (SCIP_Real) j;
         weight2 += val;

         if ( ! SCIPisFeasZero(scip, val) )
            ++cnt;
      }

      assert( cnt >= 2 );
      assert( !SCIPisFeasZero(scip, weight2) );
      w = weight1/weight2;  /*lint !e795*/

      ind = (int) SCIPfloor(scip, w);
      assert( 0 <= ind && ind < nvars-1 );

      /* branch on variable ind: either all variables up to ind or all variables after ind are zero */
      SCIPdebugMsg(scip, "Branching on variable <%s>.\n", SCIPvarGetName(vars[ind]));

      /* calculate node selection and objective estimate for node 1 */
      nodeselest = 0.0;
      objest = 0.0;
      for (j = 0; j <= ind; ++j)
      {
         nodeselest += SCIPcalcNodeselPriority(scip, vars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
         objest += SCIPcalcChildEstimate(scip, vars[j], 0.0);
      }
      /* take the average of the individual estimates */
      objest = objest/(SCIP_Real)(ind + 1.0);

      /* create node 1 */
      SCIP_CALL( SCIPcreateChild(scip, &node1, nodeselest, objest) );
      for (j = 0; j <= ind; ++j)
      {
         SCIP_CALL( fixVariableZeroNode(scip, vars[j], node1, &infeasible) );
         assert( ! infeasible );
      }

      /* calculate node selection and objective estimate for node 1 */
      nodeselest = 0.0;
      objest = 0.0;
      for (j = ind+1; j < nvars; ++j)
      {
         nodeselest += SCIPcalcNodeselPriority(scip, vars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
         objest += SCIPcalcChildEstimate(scip, vars[j], 0.0);
      }
      /* take the average of the individual estimates */
      objest = objest/((SCIP_Real) (nvars - ind - 1));

      /* create node 2 */
      SCIP_CALL( SCIPcreateChild(scip, &node2, nodeselest, objest) );
      for (j = ind+1; j < nvars; ++j)
      {
         SCIP_CALL( fixVariableZeroNode(scip, vars[j], node2, &infeasible) );
         assert( ! infeasible );
      }
   }
   SCIP_CALL( SCIPresetConsAge(scip, branchCons) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler */
static
SCIP_RETCODE enforceSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< indicator constraints */
   SCIP_SOL*             sol,                /**< solution to be enforced (NULL for LP solution) */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );
   
   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->addcomps && conshdlrdata->fixnonzero )
   {
      SCIPerrorMessage("Incompatible parameter setting: addcomps = TRUE and fixnonzero = TRUE.\n");
      return SCIP_PARAMETERWRONGVAL;
   }

   if ( conshdlrdata->fixnonzero && ( conshdlrdata->branchingrule == 'b' || conshdlrdata->branchingrule == 's' ) )
   {
      SCIPerrorMessage("Incompatible parameter setting: nonzero fixing is not compatible with bipartite or sos1 branching.\n");
      return SCIP_PARAMETERWRONGVAL;
   }

   if ( conshdlrdata->branchingrule == 's' && conshdlrdata->nstrongrounds != 0 )
   {
      SCIPerrorMessage("Strong branching is not available for SOS1 branching.\n");
      return SCIP_PARAMETERWRONGVAL;
   }

   if ( conshdlrdata->branchingrule == 's' || conshdlrdata->switchsos1branch )
   {
      /* enforce SOS1 constraints */
      SCIP_CALL( enforceConssSOS1(scip, conshdlr, nconss, conss, sol, result) );
   }
   else
   {
      if ( conshdlrdata->branchingrule != 'n' && conshdlrdata->branchingrule != 'b' )
      {
         SCIPerrorMessage("branching rule %c unknown\n", conshdlrdata->branchingrule);
         return SCIP_PARAMETERWRONGVAL;
      }

      /* enforce conflict graph */
      SCIP_CALL( enforceConflictgraph(scip, conshdlrdata, conshdlr, nconss, conss, sol, result) );
   }

   return SCIP_OKAY;
}


/* ----------------------------- separation ------------------------------------*/

/** initialitze tclique graph and create clique data */
static
SCIP_RETCODE initTCliquegraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   nsos1vars           /**< number of SOS1 variables */
   )
{
   TCLIQUE_DATA* tcliquedata;
   int j;

   /* try to generate bound cuts */
   if ( ! tcliqueCreate(&conshdlrdata->tcliquegraph) )
      return SCIP_NOMEMORY;

   /* add nodes */
   for (j = 0; j < nsos1vars; ++j)
   {
      if ( ! tcliqueAddNode(conshdlrdata->tcliquegraph, j, 0 ) )
         return SCIP_NOMEMORY;
   }

   /* add edges */
   for (j = 0; j < nsos1vars; ++j)
   {
      int* succ;
      int nsucc;
      int succnode;
      int i;

      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, j);
      succ = SCIPdigraphGetSuccessors(conflictgraph, j);

      for (i = 0; i < nsucc; ++i)
      {
         succnode = succ[i];

         if ( succnode > j && SCIPvarIsActive(SCIPnodeGetVarSOS1(conflictgraph, succnode)) )
         {
            if ( ! tcliqueAddEdge(conshdlrdata->tcliquegraph, j, succnode) )
               return SCIP_NOMEMORY;
         }
      }
   }
   if ( ! tcliqueFlush(conshdlrdata->tcliquegraph) )
      return SCIP_NOMEMORY;


   /* allocate clique data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata->tcliquedata) );
   tcliquedata = conshdlrdata->tcliquedata;

   /* initialize clique data */
   tcliquedata->scip = scip;
   tcliquedata->sol = NULL;
   tcliquedata->conshdlr = conshdlr;
   tcliquedata->conflictgraph = conflictgraph;
   tcliquedata->scaleval = 1000.0;
   tcliquedata->ncuts = 0;
   tcliquedata->nboundcuts = conshdlrdata->nboundcuts;
   tcliquedata->strthenboundcuts = conshdlrdata->strthenboundcuts;
   tcliquedata->maxboundcuts = conshdlrdata->maxboundcutsroot;

   return SCIP_OKAY;
}


/** update weights of tclique graph */
static
SCIP_RETCODE updateWeightsTCliquegraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   TCLIQUE_DATA*         tcliquedata,        /**< tclique data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   int                   nsos1vars           /**< number of SOS1 variables */
   )
{
   SCIP_Real scaleval;
   int j;

   scaleval = tcliquedata->scaleval;

   for (j = 0; j < nsos1vars; ++j)
   {
      SCIP_Real solval;
      SCIP_Real bound;
      SCIP_VAR* var;

      var = SCIPnodeGetVarSOS1(conflictgraph, j);
      solval = SCIPgetSolVal(scip, sol, var);

      if ( SCIPisFeasPositive(scip, solval) )
      {
         if ( conshdlrdata->strthenboundcuts )
            bound = REALABS( nodeGetSolvalVarboundUbSOS1(scip, conflictgraph, sol, j) );
         else
            bound = REALABS( SCIPvarGetUbLocal(var) );
      }
      else if ( SCIPisFeasNegative(scip, solval) )
      {
         if ( conshdlrdata->strthenboundcuts )
            bound = REALABS( nodeGetSolvalVarboundLbSOS1(scip, conflictgraph, sol, j) );
         else
            bound = REALABS( SCIPvarGetLbLocal(var) );
      }
      else
         bound = 0.0;

      solval = REALABS( solval );

      if ( ! SCIPisFeasZero(scip, bound) && ! SCIPisInfinity(scip, bound) )
      {
         SCIP_Real nodeweight;
         nodeweight = REALABS( solval/bound ) * scaleval;/*lint !e414*/
         tcliqueChangeWeight(conshdlrdata->tcliquegraph, j, (int)nodeweight);
      }
      else
      {
         tcliqueChangeWeight(conshdlrdata->tcliquegraph, j, 0);
      }
   }

   return SCIP_OKAY;
}


/** adds bound cut(s) to separation storage */
static
SCIP_RETCODE addBoundCutSepa(
   SCIP*                 scip,               /**< SCIP pointer */
   TCLIQUE_DATA*         tcliquedata,        /**< clique data */
   SCIP_ROW*             rowlb,              /**< row for lower bounds (or NULL) */
   SCIP_ROW*             rowub,              /**< row for upper bounds (or NULL) */
   SCIP_Bool*            success,            /**< pointer to store if bound cut was added */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff occurred */
   )
{
   assert( scip != NULL );
   assert( tcliquedata != NULL );
   assert( success != NULL);
   assert( cutoff != NULL );

   *success = FALSE;
   *cutoff = FALSE;

   /* add cut for lower bounds */
   if ( rowlb != NULL )
   {
      if ( ! SCIProwIsInLP(rowlb) && SCIPisCutEfficacious(scip, NULL, rowlb) )
      {
         SCIP_Bool infeasible;

         SCIP_CALL( SCIPaddRow(scip, rowlb, FALSE, &infeasible) );
         if ( infeasible )
            *cutoff = TRUE;
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rowlb, NULL) ) );
         ++tcliquedata->nboundcuts;
         ++tcliquedata->ncuts;
         *success = TRUE;
      }
   }

   /* add cut for upper bounds */
   if ( rowub != NULL )
   {
      if ( ! SCIProwIsInLP(rowub) && SCIPisCutEfficacious(scip, NULL, rowub) )
      {
         SCIP_Bool infeasible;

         SCIP_CALL( SCIPaddRow(scip, rowub, FALSE, &infeasible) );
         if ( infeasible )
            *cutoff = TRUE;
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rowub, NULL) ) );
         ++tcliquedata->nboundcuts;
         ++tcliquedata->ncuts;
         *success = TRUE;
      }
   }

   return SCIP_OKAY;
}


/** Generate bound constraint
 *
 *  We generate the row corresponding to the following simple valid inequalities:
 *  \f[
 *         \frac{x_1}{u_1} + \ldots + \frac{x_n}{u_n} \leq 1\qquad\mbox{and}\qquad
 *         \frac{x_1}{\ell_1} + \ldots + \frac{x_n}{\ell_1} \leq 1,
 *  \f]
 *  where \f$\ell_1, \ldots, \ell_n\f$ and \f$u_1, \ldots, u_n\f$ are the nonzero and finite lower and upper bounds of
 *  the variables \f$x_1, \ldots, x_n\f$. If an upper bound < 0 or a lower bound > 0, the constraint itself is
 *  redundant, so the cut is not applied (lower bounds > 0 and upper bounds < 0 are usually detected in presolving or
 *  propagation). Infinite bounds and zero are skipped. Thus \f$\ell_1, \ldots, \ell_n\f$ are all negative, which
 *  results in the \f$\leq\f$ inequality. In case of the presence of variable upper bounds, the bound inequality can
 *  be further strengthened.
 *
 *  Note that in fact, any mixture of nonzero finite lower and upper bounds would lead to a valid inequality as
 *  above. However, usually either the lower or upper bound is nonzero. Thus, the above inequalities are the most
 *  interesting.
 */
static
SCIP_RETCODE generateBoundInequalityFromSOS1Nodes(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int*                  nodes,              /**< conflict graph nodes for bound constraint */
   int                   nnodes,             /**< number of conflict graph nodes for bound constraint */
   SCIP_Real             rhs,                /**< right hand side of bound constraint */
   SCIP_Bool             local,              /**< in any case produce a local cut (even if local bounds of variables are valid globally) */
   SCIP_Bool             global,             /**< in any case produce a global cut */
   SCIP_Bool             strengthen,         /**< whether trying to strengthen bound constraint */
   SCIP_Bool             removable,          /**< should the inequality be removed from the LP due to aging or cleanup? */
   const char*           nameext,            /**< part of name of bound constraints */
   SCIP_ROW**            rowlb,              /**< output: row for lower bounds (or NULL if not needed) */
   SCIP_ROW**            rowub               /**< output: row for upper bounds (or NULL if not needed) */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* lbboundvar = NULL;
   SCIP_VAR* ubboundvar = NULL;
   SCIP_Bool locallbs;
   SCIP_Bool localubs;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conflictgraph != NULL );
   assert( ! local || ! global );
   assert( nodes != NULL );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nnodes+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nnodes+1) );

   /* take care of upper bounds */
   if ( rowub != NULL )
   {
      SCIP_Bool useboundvar;
      int cnt = 0;
      int j;

      /* Loop through all variables. We check whether all bound variables (if existent) are equal; if this is the
       * case then the bound constraint can be strengthened */
      localubs = local;
      useboundvar = strengthen;
      for (j = 0; j < nnodes; ++j)
      {
         SCIP_NODEDATA* nodedata;
         SCIP_VAR* var;
         SCIP_Real val;

         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, nodes[j]);
         assert( nodedata != NULL );
         var = nodedata->var;
         assert( var != NULL );

         /* if variable is not involved in a variable bound constraint */
         if ( ! useboundvar || nodedata->ubboundvar == NULL )
         {
            useboundvar = FALSE;
            if ( localubs )
            {
               assert( ! global );
               val = SCIPvarGetUbLocal(var);
            }
            else
            {
               val = SCIPvarGetUbGlobal(var);
               if ( ! global && ! SCIPisFeasEQ(scip, val, SCIPvarGetUbLocal(var)) )
               {
                  localubs = TRUE;
                  val = SCIPvarGetUbLocal(var);
               }
            }
         }
         else
         {
            /* in this case the cut is always valid globally */

            /* if we have a bound variable for the first time */
            if ( ubboundvar == NULL )
            {
               ubboundvar = nodedata->ubboundvar;
               val = nodedata->ubboundcoef;
            }
            /* else if the bound variable equals the stored bound variable */
            else if ( ubboundvar == nodedata->ubboundvar )
               val = nodedata->ubboundcoef;
            else /* else use bounds on the variables */
            {
               useboundvar = FALSE;

               /* restart 'for'-loop */
               j = -1;
               cnt = 0;
               continue;
            }
         }

         /* should not apply the cut if a variable is fixed to be negative -> constraint is redundant */
         if ( SCIPisNegative(scip, val) )
            break;

         /* store variable if relevant for bound inequality */
         if ( ! SCIPisInfinity(scip, val) && ! SCIPisZero(scip, val) )
         {
            vars[cnt] = var;

            /* if only two nodes then we scale the cut differently */
            if ( nnodes == 2 )
               vals[cnt++] = val;
            else
               vals[cnt++] = 1.0/val;
         }
      }

      /* if cut is meaningful */
      if ( j == nnodes && cnt >= 2 )/*lint !e850*/
      {
         /* if only two nodes then we scale the cut differently */
         if ( nnodes == 2 )
         {
            SCIP_Real save;

            save = vals[0];
            vals[0] = vals[1];
            vals[1] = save;
            rhs = rhs * vals[0] * vals[1];
            assert( (! useboundvar && cnt == 2 ) || (useboundvar && cnt == 3 ) );
         }

         if ( useboundvar )
         {
            /* add bound variable to array */
            vars[cnt] = ubboundvar;
            vals[cnt++] = -rhs;
            assert(ubboundvar != NULL );

            /* create upper bound inequality if at least two of the bounds are finite and nonzero */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sosub#%s", nameext);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowub, conshdlr, name, -SCIPinfinity(scip), 0.0, localubs, FALSE, removable) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *rowub, cnt, vars, vals) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowub, NULL) ) );
         }
         else
         {
            /* create upper bound inequality if at least two of the bounds are finite and nonzero */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sosub#%s", nameext);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowub, conshdlr, name, -SCIPinfinity(scip), rhs, localubs, FALSE, removable) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *rowub, cnt, vars, vals) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowub, NULL) ) );
         }
      }
   }


   /* take care of lower bounds */
   if ( rowlb != NULL )
   {
      SCIP_Bool useboundvar;
      int cnt;
      int j;

      /* loop through all variables. We check whether all bound variables (if existent) are equal; if this is the
       * case then the bound constraint can be strengthened */
      cnt = 0;
      locallbs = local;
      useboundvar = strengthen;
      for (j = 0; j < nnodes; ++j)
      {
         SCIP_NODEDATA* nodedata;
         SCIP_VAR* var;
         SCIP_Real val;

         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, nodes[j]);
         assert( nodedata != NULL );
         var = nodedata->var;
         assert( var != NULL );

         /* if variable is not involved in a variable bound constraint */
         if ( ! useboundvar || nodedata->lbboundvar == NULL )
         {
            useboundvar = FALSE;
            if ( locallbs )
            {
               assert( ! global );
               val = SCIPvarGetLbLocal(var);
            }
            else
            {
               val = SCIPvarGetLbGlobal(var);
               if ( ! global && ! SCIPisFeasEQ(scip, val, SCIPvarGetLbLocal(var)) )
               {
                  locallbs = TRUE;
                  val = SCIPvarGetUbLocal(var);
               }
            }
         }
         else
         {
            /* in this case the cut is always valid globally */

            /* if we have a bound variable for the first time */
            if ( lbboundvar == NULL )
            {
               lbboundvar = nodedata->lbboundvar;
               val = nodedata->lbboundcoef;
            }
            /* else if the bound variable equals the stored bound variable */
            else if ( SCIPvarCompare(lbboundvar, nodedata->lbboundvar) == 0 )
            {
               val = nodedata->lbboundcoef;
            }
            else /* else use bounds on the variables */
            {
               useboundvar = FALSE;

               /* restart 'for'-loop */
               j = -1;
               cnt = 0;
               continue;
            }
         }

         /* should not apply the cut if a variable is fixed to be positive -> constraint is redundant */
         if ( SCIPisPositive(scip, val) )
            break;

         /* store variable if relevant for bound inequality */
         if ( ! SCIPisInfinity(scip, -val) && ! SCIPisZero(scip, val) )
         {
            vars[cnt] = var;

            /* if only two nodes then we scale the cut differently */
            if ( nnodes == 2 )
               vals[cnt++] = val;
            else
               vals[cnt++] = 1.0/val;
         }
      }

      /* if cut is meaningful */
      if ( j == nnodes && cnt >= 2 )/*lint !e850*/
      {
         /* if only two nodes then we scale the cut differently */
         if ( nnodes == 2 )
         {
            SCIP_Real save;

            save = vals[0];
            vals[0] = vals[1];
            vals[1] = save;
            rhs = rhs * vals[0] * vals[1];
            assert( (! useboundvar && cnt == 2 ) || (useboundvar && cnt == 3 ) );
         }

         if ( useboundvar )
         {
            /* add bound variable to array */
            vars[cnt] = lbboundvar;
            vals[cnt++] = -rhs;
            assert(lbboundvar != NULL );

            /* create upper bound inequality if at least two of the bounds are finite and nonzero */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "soslb#%s", nameext);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowlb, conshdlr, name, -SCIPinfinity(scip), 0.0, locallbs, FALSE, TRUE) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *rowlb, cnt, vars, vals) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowlb, NULL) ) );
         }
         else
         {
            /* create upper bound inequality if at least two of the bounds are finite and nonzero */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "soslb#%s", nameext);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowlb, conshdlr, name, -SCIPinfinity(scip), rhs, locallbs, FALSE, TRUE) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *rowlb, cnt, vars, vals) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowlb, NULL) ) );
         }
      }
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** generates bound cuts using a clique found by algorithm for maximum weight clique
 *  and decides whether to stop generating cliques with the algorithm for maximum weight clique
 */
static
TCLIQUE_NEWSOL(tcliqueNewsolClique)
{
   TCLIQUE_WEIGHT minweightinc;

   assert( acceptsol != NULL );
   assert( stopsolving != NULL );
   assert( tcliquedata != NULL );

   /* we don't accept the solution as new incumbent, because we want to find many violated clique inequalities */
   *acceptsol = FALSE;
   *stopsolving = FALSE;

   /* slightly increase the minimal weight for additional cliques */
   minweightinc = (cliqueweight - *minweight)/10;
   minweightinc = MAX(minweightinc, 1);
   *minweight += minweightinc;

   /* adds cut if weight of the clique is greater than 1 */
   if( cliqueweight > tcliquedata->scaleval )
   {
      SCIP* scip;
      SCIP_SOL* sol;
      SCIP_Real unscaledweight;
      SCIP_Real solval;
      SCIP_Real bound;
      SCIP_VAR* var;
      int node;
      int i;

      scip = tcliquedata->scip;
      sol = tcliquedata->sol;
      assert( scip != NULL );

      /* calculate the weight of the clique in unscaled fractional variable space */
      unscaledweight = 0.0;
      for( i = 0; i < ncliquenodes; i++ )
      {
         node = cliquenodes[i];
         var = SCIPnodeGetVarSOS1(tcliquedata->conflictgraph, node);
         solval = SCIPgetSolVal(scip, sol, var);

         if ( SCIPisFeasPositive(scip, solval) )
         {
            if ( tcliquedata->strthenboundcuts )
               bound = REALABS( nodeGetSolvalVarboundUbSOS1(scip, tcliquedata->conflictgraph, sol, node) );
            else
               bound = REALABS( SCIPvarGetUbLocal(var) );
         }
         else if ( SCIPisFeasNegative(scip, solval) )
         {
            if ( tcliquedata->strthenboundcuts )
               bound = REALABS( nodeGetSolvalVarboundLbSOS1(scip, tcliquedata->conflictgraph, sol, node) );
            else
               bound = REALABS( SCIPvarGetLbLocal(var) );
         }
         else
            bound = 0.0;

         solval = REALABS( solval );

         if ( ! SCIPisFeasZero(scip, bound) && ! SCIPisInfinity(scip, bound) )
            unscaledweight += REALABS( solval/bound );/*lint !e414*/
      }

      if ( SCIPisEfficacious(scip, unscaledweight - 1.0) )
      {
         char nameext[SCIP_MAXSTRLEN];
         SCIP_ROW* rowlb = NULL;
         SCIP_ROW* rowub = NULL;
         SCIP_Bool success;
         SCIP_Bool cutoff;

         /* generate bound inequalities for lower and upper bound case
          * NOTE: tests have shown that non-removable rows give the best results */
         (void) SCIPsnprintf(nameext, SCIP_MAXSTRLEN, "%d", tcliquedata->nboundcuts);
         if ( generateBoundInequalityFromSOS1Nodes(scip, tcliquedata->conshdlr, tcliquedata->conflictgraph,
               cliquenodes, ncliquenodes, 1.0, FALSE, FALSE, tcliquedata->strthenboundcuts, FALSE, nameext, &rowlb, &rowub) != SCIP_OKAY )
         {
            SCIPerrorMessage("Unexpected error in bound cut creation.\n");
            SCIPABORT();
            return;   /*lint !e527*/
         }

         /* add bound cut(s) to separation storage if existent */
         if ( addBoundCutSepa(scip, tcliquedata, rowlb, rowub, &success, &cutoff) != SCIP_OKAY )
         {
            SCIPerrorMessage("Unexpected error in bound cut creation.\n");
            SCIPABORT();
            return;   /*lint !e527*/
         }

         if ( rowlb != NULL )
         {
            if ( SCIPreleaseRow(scip, &rowlb) != SCIP_OKAY )
            {
               SCIPerrorMessage("Cannot release row,\n");
               SCIPABORT();
               return;   /*lint !e527*/
            }
         }
         if ( rowub != NULL )
         {
            if ( SCIPreleaseRow(scip, &rowub) != SCIP_OKAY )
            {
               SCIPerrorMessage("Cannot release row,\n");
               SCIPABORT();
               return;   /*lint !e527*/
            }
         }

         /* if at least one cut has been added */
         if ( success )
         {
            SCIPdebugMsg(scip, " -> found bound cut corresponding to clique (act=%g)\n", unscaledweight);

            /* if we found more than half the cuts we are allowed to generate, we accept the clique as new incumbent,
             * such that only more violated cuts are generated afterwards
             */
            if( tcliquedata->maxboundcuts >= 0 )
            {
               if ( tcliquedata->ncuts > tcliquedata->maxboundcuts/2 )
                  *acceptsol = TRUE;
               if ( tcliquedata->ncuts >= tcliquedata->maxboundcuts )
                  *stopsolving = TRUE;
            }
         }
         else
            *stopsolving = TRUE;
      }
   }
}


/** separate bound inequalities from conflict graph */
static
SCIP_RETCODE sepaBoundInequalitiesFromGraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   int                   maxboundcuts,       /**< maximal number of bound cuts separated per separation round (-1: no limit) */
   int*                  ngen,               /**< pointer to store number of cuts generated */
   SCIP_Bool*            cutoff              /**< pointer whether a cutoff occurred */
   )
{
   SCIP_DIGRAPH* conflictgraph;
   TCLIQUE_DATA* tcliquedata;
   TCLIQUE_WEIGHT cliqueweight;
   TCLIQUE_STATUS tcliquestatus;
   int nsos1vars;

   SCIP_Real scaleval = 1000.0;                  /* factor for scaling weights */
   int maxtreenodes = 10000;                     /* maximal number of nodes of b&b tree */
   int maxzeroextensions = 1000;                 /* maximal number of zero-valued variables extending the clique (-1: no limit) */
   int backtrackfreq = 1000;                     /* frequency for premature backtracking up to tree level 1 (0: no backtracking) */
   int ntreenodes;
   int* cliquenodes;
   int ncliquenodes;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conshdlrdata != NULL );
   assert( ngen != NULL );

   /* get conflict graph */
   conflictgraph = SCIPgetConflictgraphSOS1(conshdlr);
   assert( conflictgraph != NULL );

   /* get number of SOS1 variables */
   nsos1vars = SCIPgetNSOS1Vars(conshdlr);

   /* initialize data of tclique graph*/
   tcliquedata = conshdlrdata->tcliquedata;
   tcliquedata->scaleval = scaleval;
   tcliquedata->maxboundcuts = maxboundcuts;
   tcliquedata->sol = sol;
   tcliquedata->ncuts = 0;
   tcliquedata->cutoff = FALSE;

   /* update the weights of the tclique graph */
   SCIP_CALL( updateWeightsTCliquegraph(scip, conshdlrdata, tcliquedata, conflictgraph, sol, nsos1vars) );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquenodes, nsos1vars) );

   /* start algorithm to find maximum weight cliques and use them to generate bound cuts */
   tcliqueMaxClique(tcliqueGetNNodes, tcliqueGetWeights, tcliqueIsEdge, tcliqueSelectAdjnodes,
      conshdlrdata->tcliquegraph, tcliqueNewsolClique, tcliquedata,
      cliquenodes, &ncliquenodes, &cliqueweight, (int)scaleval-1, (int)scaleval+1,
      maxtreenodes, backtrackfreq, maxzeroextensions, -1, &ntreenodes, &tcliquestatus);

   /* free buffer array */
   SCIPfreeBufferArray(scip, &cliquenodes);

   /* get number of cuts of current separation round */
   *ngen = tcliquedata->ncuts;

   /* store whether a cutoff occurred */
   *cutoff = tcliquedata->cutoff;

   /* update number of bound cuts in separator data */
   conshdlrdata->nboundcuts = tcliquedata->nboundcuts;

   return SCIP_OKAY;
}


/** Generate a bound constraint from the variables of an SOS1 constraint (see generateBoundInequalityFromSOS1Nodes() for more information) */
static
SCIP_RETCODE generateBoundInequalityFromSOS1Cons(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< SOS1 constraint */
   SCIP_Bool             local,              /**< in any case produce a local cut (even if local bounds of variables are valid globally) */
   SCIP_Bool             global,             /**< in any case produce a global cut */
   SCIP_Bool             strengthen,         /**< whether trying to strengthen bound constraint */
   SCIP_Bool             removable,          /**< should the inequality be removed from the LP due to aging or cleanup? */
   SCIP_ROW**            rowlb,              /**< output: row for lower bounds (or NULL if not needed) */
   SCIP_ROW**            rowub               /**< output: row for upper bounds (or NULL if not needed) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int* nodes;
   int nvars;
   int cnt = 0;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   nvars = consdata->nvars;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nvars) );

   /* get nodes in the conflict graph */
   for (j = 0; j < nvars; ++j)
   {
      if ( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(consdata->vars[j])) || SCIPisFeasPositive(scip, SCIPvarGetUbLocal(consdata->vars[j])) )
      {
         assert( varGetNodeSOS1(conshdlrdata, consdata->vars[j]) >= 0 );
         nodes[cnt++] = varGetNodeSOS1(conshdlrdata, consdata->vars[j]);
      }
   }

   /* generate bound constraint from conflict graph nodes */
   if ( cnt > 0 )
   {
      SCIP_CALL( generateBoundInequalityFromSOS1Nodes(scip, conshdlr, conshdlrdata->conflictgraph, nodes, cnt, 1.0, local, global,
            strengthen, removable, SCIPconsGetName(cons), rowlb, rowub) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &nodes);

   return SCIP_OKAY;
}


/** initialize or separate bound inequalities from SOS1 constraints */
static
SCIP_RETCODE initsepaBoundInequalityFromSOS1Cons(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   int                   nconss,             /**< number of SOS1 constraints */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   SCIP_Bool             solvedinitlp,       /**< TRUE if initial LP relaxation at a node is solved */
   int                   maxboundcuts,       /**< maximal number of bound cuts separated per separation round (-1: no limit) */
   int*                  ngen,               /**< pointer to store number of cuts generated (or NULL) */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff occurred */
   )
{
   int cnt = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( conss != NULL );

   *cutoff = FALSE;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_ROW* rowub = NULL;
      SCIP_ROW* rowlb = NULL;
      SCIP_Bool release = FALSE;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      if ( solvedinitlp )
      {
         SCIPdebugMsg(scip, "Separating inequalities for SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );
      }
      else
      {
         SCIPdebugMsg(scip, "Checking for initial rows for SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );
      }

      /* in case that the SOS1 constraint is local, we always generate new rows - the former rows might be invalid;
       * otherwise if the SOS1 constraint is global, we only generate rows if not yet done */
      if ( consdata->local )
      {
         SCIP_CALL( generateBoundInequalityFromSOS1Cons(scip, conshdlr, conss[c], TRUE, FALSE, TRUE, FALSE, &rowlb, &rowub) );
         release = TRUE;
      }
      else
      {
         if ( consdata->rowub == NULL || consdata->rowlb == NULL )
         {
            SCIP_CALL( generateBoundInequalityFromSOS1Cons(scip, conshdlr, conss[c], FALSE, TRUE, TRUE, FALSE,
                  (consdata->rowlb == NULL) ? &consdata->rowlb : NULL,
                  (consdata->rowub == NULL) ? &consdata->rowub : NULL) ); /*lint !e826*/
         }
         rowub = consdata->rowub;
         rowlb = consdata->rowlb;
      }

      /* put corresponding rows into LP */
      if ( rowub != NULL && ! SCIProwIsInLP(rowub) && ( solvedinitlp || SCIPisCutEfficacious(scip, sol, rowub) ) )
      {
         SCIP_CALL( SCIPaddRow(scip, rowub, FALSE, cutoff) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rowub, NULL) ) );

         if ( solvedinitlp )
         {
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
         }
         ++cnt;
      }

      if ( ! (*cutoff) && rowlb != NULL && ! SCIProwIsInLP(rowlb) && ( solvedinitlp || SCIPisCutEfficacious(scip, sol, rowlb) ) )
      {
         SCIP_CALL( SCIPaddRow(scip, rowlb, FALSE, cutoff) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rowlb, NULL) ) );

         if ( solvedinitlp )
         {
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
         }
         ++cnt;
      }

      /* release rows if they are local */
      if ( release )
      {
         if ( rowlb != NULL )
         {
            SCIP_CALL( SCIPreleaseRow(scip, &rowlb) );
         }
         if ( rowub != NULL )
         {
            SCIP_CALL( SCIPreleaseRow(scip, &rowub) );
         }
      }

      if ( *cutoff || ( maxboundcuts >= 0 && cnt >= maxboundcuts ) )
         break;
   }

   /* store number of generated cuts */
   if ( ngen != NULL )
      *ngen = cnt;

   return SCIP_OKAY;
}


/** separate implied bound cuts */
static
SCIP_RETCODE sepaImplBoundCutsSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   int                   maxcuts,            /**< maximal number of implied bound cuts separated per separation round (-1: no limit) */
   int*                  ngen,               /**< pointer to store number of cuts generated */
   SCIP_Bool*            cutoff              /**< pointer whether a cutoff occurred */
   )
{
   SCIP_DIGRAPH* implgraph;
   SCIP_Bool genbreak;
   int nimplnodes;
   int i;

   assert( scip != NULL);
   assert( conshdlrdata != NULL);
   assert( conshdlr != NULL);
   assert( ngen != NULL);
   assert( cutoff != NULL);

   *cutoff = FALSE;
   *ngen = 0;

   /* return if conflict graph is not available */
   if ( conshdlrdata->conflictgraph == NULL )
      return SCIP_OKAY;

   /* get implication graph  */
   implgraph = conshdlrdata->implgraph;

   /* create implication graph if not done already */
   if ( implgraph == NULL )
   {
      int nchbds;

      if ( SCIPgetDepth(scip) == 0 )
      {
         SCIP_Bool success;
         SCIP_CALL( initImplGraphSOS1(scip, conshdlrdata, conshdlrdata->conflictgraph, conshdlrdata->nsos1vars, conshdlrdata->maxtightenbds, &nchbds, cutoff, &success) );
         if ( *cutoff || ! success )
            return SCIP_OKAY;
         implgraph = conshdlrdata->implgraph;
      }
      else
      {
         return SCIP_OKAY;
      }
   }
   nimplnodes = conshdlrdata->nimplnodes;
   assert( implgraph != NULL );
   assert( nimplnodes > 0);

   /* exit if implication graph has no arcs between its nodes */
   if ( SCIPdigraphGetNArcs(implgraph) < 1 )
      return SCIP_OKAY;

   /* loop through all nodes of the implication graph */
   genbreak = FALSE;
   for (i = 0; i < nimplnodes && ! genbreak; ++i)
   {
      SCIP_SUCCDATA** succdatas;
      SCIP_NODEDATA* nodedata;
      SCIP_Real solval;
      SCIP_VAR* var;
      int* succ;
      int nsucc;
      int s;

      succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, i);
      nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(implgraph, i);
      assert( nodedata != NULL );
      var = nodedata->var;
      assert( var != NULL );
      solval = SCIPgetSolVal(scip, sol, var);

      if ( succdatas != NULL && ! SCIPisFeasZero(scip, solval) )
      {
         succ = SCIPdigraphGetSuccessors(implgraph, i);
         nsucc = SCIPdigraphGetNSuccessors(implgraph, i);

         for (s = 0; s < nsucc && ! genbreak; ++s)
         {
            SCIP_SUCCDATA* succdata;
            SCIP_VAR* succvar;
            SCIP_ROW* cut = NULL;
            SCIP_Bool bound1lower;
            SCIP_Bool bound2lower;
            SCIP_Real solvalsucc;
            SCIP_Real bound1;
            SCIP_Real bound2;
            SCIP_Real lhsrhs;
            SCIP_Real impl;
            int k;

            nodedata = (SCIP_NODEDATA*) SCIPdigraphGetNodeData(implgraph, succ[s]);
            succdata = succdatas[s];
            assert( nodedata != NULL && succdata != NULL && nodedata->var != NULL );
            succvar = nodedata->var;
            solvalsucc = SCIPgetSolVal(scip, sol, succvar);

            /* determine coefficients for bound inequality */
            assert( ! SCIPisFeasZero(scip, solval) );
            if ( SCIPisFeasNegative(scip, solval) )
            {
               bound1lower = TRUE;
               bound1 = SCIPvarGetLbGlobal(var);
            }
            else
            {
               bound1lower = FALSE;
               bound1 = SCIPvarGetUbGlobal(var);
            }

            /* handle lower bound upper bound implications */
            for (k = 0; k < 2; ++k)
            {
               if ( k == 0 )
               {
                  SCIP_Real lbsucc;
                  lbsucc = SCIPvarGetLbGlobal(succvar);
                  if ( SCIPisFeasLT(scip, lbsucc, succdata->lbimpl) )
                  {
                     impl = succdata->lbimpl;
                     bound2 = lbsucc;
                  }
                  else
                     continue;
               }
               else
               {
                  SCIP_Real ubsucc;
                  ubsucc = SCIPvarGetUbGlobal(succvar);
                  if ( SCIPisFeasGT(scip, ubsucc, succdata->ubimpl) )
                  {
                     impl = succdata->ubimpl;
                     bound2 = ubsucc;
                  }
                  else
                     continue;
               }

               if ( SCIPisInfinity(scip, REALABS(bound1)) || SCIPisInfinity(scip, REALABS(bound2)) )
                  continue;
               assert( ! SCIPisInfinity(scip, REALABS(impl)) );

               if ( SCIPisFeasNegative(scip, bound2-impl) )
                  bound2lower = TRUE;
               else
                  bound2lower = FALSE;

               /* determine left/right hand side of bound inequality */
               lhsrhs = bound1 * bound2;

               /* create cut */
               if ( bound1lower == bound2lower )
               {
                  if ( SCIPisFeasGT(scip, solval * (bound2-impl) + solvalsucc * bound1, lhsrhs) )
                  {
                     SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "", -SCIPinfinity(scip), lhsrhs, FALSE, FALSE, TRUE) );
                  }
                  else
                     continue;
               }
               else
               {
                  if ( SCIPisFeasLT(scip, solval * (bound2-impl) + solvalsucc * bound1, lhsrhs) )
                  {
                     SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "", lhsrhs, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
                  }
                  else
                     continue;
               }

               /* add coefficients of variables */
               SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
               SCIP_CALL( SCIPaddVarToRow(scip, cut, var, bound2-impl) );
               SCIP_CALL( SCIPaddVarToRow(scip, cut, succvar, bound1) );
               SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

               /* add cut if useful */
               if ( ! SCIProwIsInLP(cut) && SCIPisCutEfficacious(scip, NULL, cut) )
               {
                  SCIP_Bool infeasible;
                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );
                  if ( infeasible )
                  {
                     genbreak = TRUE;
                     *cutoff = TRUE;
                     break;
                  }
                  SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );
#ifdef SCIP_DEBUG
                  if ( k == 0 )
                  {
                     SCIPdebugMsg(scip, "added cut for implication %s != 0 -> %s >= %f \n", SCIPvarGetName(var), SCIPvarGetName(succvar), succdata->lbimpl);
                  }
                  else
                  {
                     SCIPdebugMsg(scip, "added cut for implication %s != 0 -> %s <= %f \n", SCIPvarGetName(var), SCIPvarGetName(succvar), succdata->ubimpl);
                  }
#endif

                  ++(*ngen);
               }

               if ( maxcuts >= 0 && *ngen > maxcuts )
               {
                  genbreak = TRUE;
                  break;
               }
            }

            if ( cut != NULL )
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
      }
   }

   return SCIP_OKAY;
}


/** separates SOS1 constraints for arbitrary solutions */
static
SCIP_RETCODE separateSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL) */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int depth;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( nconss == 0 )
      return SCIP_OKAY;

   /* only separate cuts if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* get node depth */
   depth = SCIPgetDepth(scip);


   /* separate bound (clique) inequalities */
   if ( conshdlrdata->boundcutsfreq >= 0 && ( (conshdlrdata->boundcutsfreq == 0 && depth == 0) || (conshdlrdata->boundcutsfreq > 0 && depth % conshdlrdata->boundcutsfreq == 0)) )
   {
      int maxboundcuts;
      int ngen = 0;

      /* determine maximal number of cuts*/
      if ( depth == 0 )
         maxboundcuts = conshdlrdata->maxboundcutsroot;
      else
         maxboundcuts = conshdlrdata->maxboundcuts;

      if ( maxboundcuts >= 1 )
      {
         /* separate bound inequalities from SOS1 constraints */
         if( conshdlrdata->boundcutsfromsos1 || conshdlrdata->switchcutsfromsos1 )
         {
            SCIP_Bool cutoff;

            SCIP_CALL( initsepaBoundInequalityFromSOS1Cons(scip, conshdlr, conshdlrdata, conss, nconss, sol, TRUE, maxboundcuts, &ngen, &cutoff) );
            if ( cutoff )
            {
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
         }

         /* separate bound inequalities from the conflict graph */
         if( conshdlrdata->boundcutsfromgraph && ! conshdlrdata->switchcutsfromsos1 )
         {
            SCIP_Bool cutoff;
            SCIP_CALL( sepaBoundInequalitiesFromGraph(scip, conshdlr, conshdlrdata, sol, maxboundcuts, &ngen, &cutoff) );
            if ( cutoff )
            {
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
         }
      }

      /* evaluate results */
      if ( ngen > 0 )
         *result = SCIP_SEPARATED;
      SCIPdebugMsg(scip, "Separated %d bound (clique) inequalities.\n", ngen);
   }


   /* separate implied bound inequalities */
   if ( conshdlrdata->implcutsfreq >= 0 && ( (conshdlrdata->implcutsfreq == 0 && depth == 0) || (conshdlrdata->implcutsfreq > 0 && depth % conshdlrdata->implcutsfreq == 0)) )
   {
      int maximplcuts;
      int ngen = 0;

      /* determine maximal number of cuts*/
      if ( depth == 0 )
         maximplcuts = conshdlrdata->maximplcutsroot;
      else
         maximplcuts = conshdlrdata->maximplcuts;

      /* call separator for implied bound cuts */
      if ( maximplcuts >= 1 )
      {
         SCIP_Bool cutoff;
         SCIP_CALL( sepaImplBoundCutsSOS1(scip, conshdlr, conshdlrdata, sol, maximplcuts, &ngen, &cutoff) );
         if ( cutoff )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      }

      /* evaluate results */
      if ( ngen > 0 )
         *result = SCIP_SEPARATED;
      SCIPdebugMsg(scip, "Separated %d implied bound inequalities.\n", ngen);
   }

   return SCIP_OKAY;
}


/* -------------------------- heuristic methods --------------------------------*/

/** gets weights determining an order of the variables in a heuristic for the maximum weighted independent set problem */
static
SCIP_RETCODE getVectorOfWeights(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SOL*             sol,                /**< primal solution or NULL for current LP solution */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_Bool*            indicatorzero,      /**< vector that indicates which variables are currently fixed to zero */
   SCIP_Real*            weights             /**< pointer to store weights determining the order of the variables (length = nsos1vars) */
   )
{
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real sum;
   int nviols;
   int* succ;
   int nsucc;
   int i;
   int j;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( indicatorzero != NULL );
   assert( weights != NULL );

   for (i = 0; i < nsos1vars; ++i)
   {
      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, i);

      if( nsucc == 0 || indicatorzero[i] )
         weights[i] = 0.0;
      else
      {
         var = SCIPnodeGetVarSOS1(conflictgraph, i);
         val = REALABS( SCIPgetSolVal(scip, sol, var) );
         if ( SCIPisFeasZero(scip, val) )
            weights[i] = 0.0;
         else
         {
            succ = SCIPdigraphGetSuccessors(conflictgraph, i);

            nviols = 0;
            sum = 0.0;
            for (j = 0; j < nsucc; ++j)
            {
               SCIP_Real valsucc;

               valsucc = REALABS( SCIPgetSolVal(scip, sol, SCIPnodeGetVarSOS1(conflictgraph, succ[j])) );
               if( ! SCIPisFeasZero(scip, valsucc) )
               {
                  sum += MIN(10E05, valsucc);
                  ++nviols;
               }
            }

            if ( nviols == 0 )
               weights[i] = 0.0;
            else
            {
               assert( SCIPisFeasPositive(scip, sum * (SCIP_Real)nviols));
               val = MIN(1e6, val);
               weights[i] = ( val + SCIPsumepsilon(scip) ) / ( sum * (SCIP_Real)nviols + SCIPsumepsilon(scip) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/* marks neighbors of a given node as not a member of the maximal independent set */
static
SCIP_RETCODE markNeighborsMWISHeuristic(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   node,               /**< node of the conflict graph */
   SCIP_Bool*            mark,               /**< indicator vector of processed nodes */
   SCIP_Bool*            indset,             /**< indicator vector of current independent */
   int*                  cnt,                /**< pointer to store number of marked nodes */
   SCIP_Bool*            cutoff              /**< pointer to store whether operation is infeasible */
   )
{
   int nsucc;
   int* succ;
   int j;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( mark != NULL );
   assert( indset != NULL );
   assert( cutoff != NULL );
   assert( cnt != NULL );

   *cutoff = FALSE;

   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, node);
   succ = SCIPdigraphGetSuccessors(conflictgraph, node);

   /* for all successors */
   for (j = 0; j < nsucc && !(*cutoff); ++j)
   {
      int succj;

      succj = succ[j];
      assert( indset[succj] == 0 );
      if( ! mark[succj] )
      {
         SCIP_VARSTATUS varstatus;
         SCIP_VAR* var;

         /* mark node as processed */
         mark[succj] = TRUE;
         ++(*cnt);

         /* get variable and variable status corresponding to successor node */
         var = SCIPnodeGetVarSOS1(conflictgraph, succj);
         varstatus = SCIPvarGetStatus(var);

         /* if variable is aggregated */
         if ( varstatus == SCIP_VARSTATUS_AGGREGATED )
         {
            int aggrnode;

            aggrnode = SCIPvarGetNodeSOS1(conshdlr, SCIPvarGetAggrVar(var));

            /* if aggregated variable is an SOS1 variable */
            if ( aggrnode >= 0 )
            {
               /* if aggregated variable is implied to be zero */
               if ( SCIPisFeasZero(scip, SCIPvarGetAggrConstant(var)) )
               {
                  if ( ! mark[aggrnode] )
                  {
                     mark[aggrnode] = TRUE;
                     ++(*cnt);
                  }
                  else if ( indset[aggrnode] == 1 )
                  {
                     *cutoff = TRUE;
                     return SCIP_OKAY;
                  }
               }
               else
               {
                  /* if aggregated variable is not already a member of the maximal independent set */
                  if ( indset[aggrnode] == 0 )
                  {
                     /* if variable is already marked */
                     if ( mark[aggrnode] )
                     {
                        *cutoff = TRUE;
                        return SCIP_OKAY;
                     }
                     else
                     {
                        indset[aggrnode] = 1;
                        mark[aggrnode] = TRUE;
                        ++(*cnt);
                     }

                     /* mark neighbors of aggregated variable */
                     SCIP_CALL( markNeighborsMWISHeuristic(scip, conshdlr, conflictgraph, aggrnode, mark, indset, cnt, cutoff) );
                  }
               }
            }
         }
         else if ( varstatus == SCIP_VARSTATUS_NEGATED )
         {
            int negnode;

            negnode = SCIPvarGetNodeSOS1(conshdlr, SCIPvarGetNegationVar(var));

            /* if negated variable is an SOS1 variable */
            if ( negnode >= 0 )
            {
               if ( SCIPisFeasZero(scip, SCIPvarGetNegationConstant(var) ) )
               {
                  if ( indset[negnode] == 1 )
                  {
                     *cutoff = TRUE;
                     return SCIP_OKAY;
                  }
                  else if ( ! mark[negnode] )
                  {
                     mark[negnode] = TRUE;
                     ++(*cnt);
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** calls greedy algorithm for the maximum weighted independent set problem (MWIS)
 *
 * We compute a feasible solution to
 * \f[
 *  \begin{array}{ll}
 *  \min\limits_{z} & {x^*}^T z \\
 *                  & z_i + z_j \leq 1, \qquad (i,j)\in E \\
 *                  &       z_i \in  \{0,1\}, \qquad\quad  i\in V
 * \end{array}
 * \f]
 * by the algorithm GGWMIN of Shuichi Sakai, Mitsunori Togasaki and Koichi Yamazaki in "A note on greedy algorithms for the
 * maximum weighted independent set problem", Discrete Applied Mathematics. Here \f$x^*\f$ denotes the current LP
 * relaxation solution. Note that the solution of the MWIS is the indicator vector of an independent set.
 */
static
SCIP_RETCODE maxWeightIndSetHeuristic(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SOL*             sol,                /**< primal solution or NULL for current LP solution */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_Bool*            indicatorzero,      /**< vector that indicates which variables are currently fixed to zero */
   SCIP_Bool*            indset              /**< pointer to store indicator vector of an independent set */
   )
{
   SCIP_Bool* mark = NULL;
   SCIP_Real* weights = NULL;
   int* indscipvars = NULL;
   int ind;
   int nsucc;
   int i;
   int k;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( indicatorzero != NULL );
   assert( indset != NULL );

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &mark, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indscipvars, nsos1vars) );

   /* sort SOS1 variables in nonincreasing order of weights */
   for (i = 0; i < nsos1vars; ++i)
      indscipvars[i] = i;

   SCIP_CALL( getVectorOfWeights(scip, sol, conflictgraph, nsos1vars, indicatorzero, weights) );
   SCIPsortDownRealInt(weights, indscipvars, nsos1vars);

   /* mark fixed variables and variables without any neighbors in the conflict graph */
   k = 0;
   for (i = 0; i < nsos1vars; ++i)
   {
      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, i);

      if ( indset[i] == 0 )
      {
         if( indicatorzero[i] )
         {
            mark[i] = TRUE;
            ++k;
         }
         else if ( nsucc == 0 )
         {
            indset[i] = 1;
            mark[i] = TRUE;
            ++k;
         }
         else
            mark[i] = FALSE;
      }
      else
      {
         SCIP_Bool cutoff;

         ++k;
         mark[i] = TRUE;

         SCIP_CALL( markNeighborsMWISHeuristic(scip, conshdlr, conflictgraph, i, mark, indset, &k, &cutoff) );
         assert( ! cutoff );
      }
   }

   /* mark vertices in the order of their largest weight */
   for (i = 0; k < nsos1vars; ++i) /*lint !e440*/
   {
      assert( i < nsos1vars );

      ind = indscipvars[i];

      if ( ! mark[ind] )
      {
         SCIP_Bool cutoff;

         /* mark ind */
         indset[ind] = 1;
         mark[ind] = TRUE;
         ++k;

         SCIP_CALL( markNeighborsMWISHeuristic(scip, conshdlr, conflictgraph, ind, mark, indset, &k, &cutoff) );
         if ( cutoff )
            indset[ind] = 0;
      }
   }
   assert( k == nsos1vars );

   /* free buffer arrays */
   SCIPfreeBufferArrayNull(scip, &indscipvars);
   SCIPfreeBufferArrayNull(scip, &weights);
   SCIPfreeBufferArrayNull(scip, &mark);

   return SCIP_OKAY;
}


/** based on solution values of the variables, fixes variables of the conflict graph to zero to turn all SOS1 constraints feasible
 *
 *  if the SOS1 constraints do not overlap, the method makeSOS1constraintsFeasible() may be faster
 */
static
SCIP_RETCODE makeSOS1conflictgraphFeasible(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            changed,            /**< pointer to store whether the solution has been changed */
   SCIP_Bool*            allroundable        /**< pointer to store whether all variables are roundable */
   )
{
   SCIP_DIGRAPH* conflictgraph;  /* conflict graph for SOS1 constraints */
   SCIP_Bool* indicatorzero;     /* indicates which solution values are zero */
   SCIP_Bool* indset;            /* indicator vector of feasible solution; i.e., an independent set */
   int nsos1vars;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( sol != NULL );
   assert( changed != NULL );

   *allroundable = TRUE;
   *changed = FALSE;

   /* get number of SOS1 variables */
   nsos1vars = SCIPgetNSOS1Vars(conshdlr);
   assert( nsos1vars >= 0 );

   /* get conflict graph */
   conflictgraph = SCIPgetConflictgraphSOS1(conshdlr);
   assert( conflictgraph != NULL );

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &indset, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indicatorzero, nsos1vars) );

   /* determine if variables with nonzero solution value are roundable */
   for (j = 0; j < nsos1vars; ++j)
   {
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;

      var = SCIPnodeGetVarSOS1(conflictgraph, j);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      indset[j] = 0;

      /* if solution value of variable is zero */
      if ( SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, var)) )
         indicatorzero[j] = TRUE;
      else
      {
         indicatorzero[j] = FALSE;

         /* if variable is not roundable */
         if ( ! SCIPvarMayRoundDown(var) && ! SCIPvarMayRoundUp(var) )
         {
            *allroundable = FALSE;
            break;
         }

         /* if bounds of variable are fixed to zero */
         if ( SCIPisFeasZero(scip, ub) && SCIPisFeasZero(scip, lb) )
            indicatorzero[j] = TRUE;
         else if ( SCIPisFeasPositive(scip, lb) || SCIPisFeasNegative(scip, ub) ) /* if variable is fixed to be nonzero */
            indset[j] = 1;
      }
   }

   /* return if at least one SOS1 variable is not roundable */
   if ( ! (*allroundable) )
   {
      SCIPfreeBufferArray(scip, &indicatorzero);
      SCIPfreeBufferArray(scip, &indset);
      return SCIP_OKAY;
   }

   /* call greedy algorithm for the maximum weighted independent set problem */
   SCIP_CALL( maxWeightIndSetHeuristic(scip, sol, conshdlr, conflictgraph, nsos1vars, indicatorzero, indset) );

   /* make solution feasible */
   for (j = 0; j < nsos1vars; ++j)
   {
      if ( indset[j] == 0 )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPnodeGetVarSOS1(conflictgraph, j), 0.0) );
         *changed = TRUE;
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &indicatorzero);
   SCIPfreeBufferArray(scip, &indset);

#ifdef SCIP_NDEBUG
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS** conss;
      int nconss;
      int c;

      conss = SCIPconshdlrGetConss(conshdlr);
      nconss = SCIPconshdlrGetNConss(conshdlr);
      for (c = 0; c < nconss; ++c)
      {
         int cnt = 0;
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );

         for (j = 0; j < consdata->nvars; ++j)
         {
            if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[j])) )
            {
               ++cnt;
            }
         }
         assert( cnt < 2 );
      }
   }
#endif

   return SCIP_OKAY;
}


/** based on solution values of the variables, fixes variables of the SOS1 constraints to zero to turn these constraints feasible
 *
 *  if the SOS1 constraints overlap, the method makeSOS1constraintsFeasible() may result in better primal solutions
 */
static
SCIP_RETCODE makeSOS1constraintsFeasible(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            changed,            /**< pointer to store whether the solution has been changed */
   SCIP_Bool*            allroundable        /**< pointer to store whether all variables are roundable */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( sol != NULL );
   assert( changed != NULL );

   *allroundable = TRUE;
   *changed = FALSE;

   /* get SOS1 constraints and number of SOS1 constraints */
   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);
   assert( nconss > 0 );

   /* loop through all SOS1 constraints */
   for (c = 0; c < nconss && *allroundable; ++c)
   {
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_Bool varisfixed = FALSE;
      SCIP_Real maxval = 0.0;
      int pos = -1;
      int nvars;
      int j;

      cons = conss[c];
      assert( cons != NULL );
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      nvars = consdata->nvars;
      vars = consdata->vars;

      /* search for maximum solution value */
      for (j = 0; j < nvars; ++j)
      {
         SCIP_VAR* var;

         var = vars[j];

         if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, var)) )
         {
            SCIP_Real lb;
            SCIP_Real ub;

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);

            /* if variable is not roundable */
            if ( ! SCIPvarMayRoundDown(var) && ! SCIPvarMayRoundUp(var) )
            {
               *allroundable = FALSE;
               break;
            }

            /* it is possible that the bounds were proagated to zero although the current solution value is nonzero
             * in this case fix the solution value to zero */
            if ( SCIPisFeasZero(scip, ub) && SCIPisFeasZero(scip, lb) )
            {
               SCIP_CALL( SCIPsetSolVal(scip, sol, var, 0.0) );
               *changed = TRUE;
            }
            else if ( SCIPisFeasPositive(scip, lb) || SCIPisFeasNegative(scip, ub) ) /* if variable is fixed to be nonzero */
            {
               assert( ! varisfixed );
               varisfixed = TRUE;
               maxval = SCIPgetSolVal(scip, sol, var);
               pos = j;
            }
            else if ( ! varisfixed && SCIPisFeasGT(scip, REALABS(SCIPgetSolVal(scip, sol, var)), REALABS(maxval)) ) /* search for variable with maximum solution value */
            {
               maxval = SCIPgetSolVal(scip, sol, var);
               pos = j;
            }

            /* fix variable to zero; the solution value of the variable with maximum solution value
             * will be restored in a later step */
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, 0.0) );
            *changed = TRUE;
         }
      }

      if ( ! (*allroundable) )
         break;
      else if ( pos >= 0 ) /* restore solution of variable with maximum solution value */
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, vars[pos], maxval) );
      }
   }

#ifdef SCIP_NDEBUG
   if ( *allroundable )
   {
      for (c = 0; c < nconss; ++c)
      {
         int cnt = 0;
         int j;

         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );

         for (j = 0; j < consdata->nvars; ++j)
         {
            if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[j])) )
            {
               ++cnt;
            }
         }
         assert( cnt < 2 );
      }
   }
#endif

   return SCIP_OKAY;
}


/** determine a diving variables and boundchanges of diving variables by analyzing the conflict graph
 *
 *  if the SOS1 constraints do not overlap, the method getDiveBdChgsSOS1constraints() may be faster
 */
static
SCIP_RETCODE getDiveBdChgsSOS1conflictgraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            success             /**< pointer to store */
   )
{
   SCIP_DIGRAPH* conflictgraph;
   SCIP_VAR* bestvar = NULL;
   SCIP_Bool bestvarfixneigh = FALSE;
   SCIP_Real bestscore = SCIP_REAL_MIN;
   int bestnode = -1;
   int nsos1vars;
   int v;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( diveset != NULL );
   assert( success != NULL );

   *success = FALSE;

   /* get number of SOS1 variables */
   nsos1vars = SCIPgetNSOS1Vars(conshdlr);

   /* get conflict graph of SOS1 constraints */
   conflictgraph = SCIPgetConflictgraphSOS1(conshdlr);

   /* loop over SOS1 variables  */
   for (v = 0; v < nsos1vars; ++v)
   {
      /* check whether the variable violates an SOS1 constraint together with at least one other variable */
      if ( isViolatedSOS1(scip, conflictgraph, v, sol) )
      {
         SCIP_VAR* var;
         SCIP_Real solval;
         SCIP_Real score;
         SCIP_Real bound;
         SCIP_Real fracval;
         SCIP_Bool fixneigh;

         var = SCIPnodeGetVarSOS1(conflictgraph, v);
         solval = SCIPgetSolVal(scip, sol, var);

         /* compute (variable) bound of candidate */
         if ( SCIPisFeasNegative(scip, solval) )
            bound = nodeGetSolvalVarboundLbSOS1(scip, conflictgraph, sol, v);
         else
            bound = nodeGetSolvalVarboundUbSOS1(scip, conflictgraph, sol, v);

         /* bound may have changed in propagation; ensure that fracval <= 1 */
         if ( SCIPisFeasLT(scip, REALABS(bound), REALABS(solval)) )
            bound = solval;

         /* ensure finiteness */
         bound = MIN(DIVINGCUTOFFVALUE, REALABS(bound)); /*lint !e666*/
         fracval = MIN(DIVINGCUTOFFVALUE, REALABS(solval)); /*lint !e666*/
         assert( ! SCIPisInfinity(scip, bound) );
         assert( ! SCIPisInfinity(scip, fracval) );
         assert( SCIPisPositive(scip, bound) );

         /* get fractionality of candidate */
         fracval /= (bound + SCIPsumepsilon(scip));

         /* should SOS1 variables be scored by the diving heuristics specific score function;
          *  otherwise use the score function of the SOS1 constraint handler */
         if ( SCIPdivesetSupportsType(diveset, SCIP_DIVETYPE_SOS1VARIABLE) )
         {
            SCIP_Bool roundup;

            SCIP_CALL( SCIPgetDivesetScore(scip, diveset, SCIP_DIVETYPE_SOS1VARIABLE, var, solval, fracval,
                  &score, &roundup) );

            fixneigh = roundup;
            if ( SCIPisFeasNegative(scip, solval) )
               fixneigh = !fixneigh;
         }
         else
         {
            /* we always fix the candidates neighbors in the conflict graph to zero */
            fixneigh = TRUE;

            /* score fractionality of candidate */
            score = fracval;
         }

         /* best candidate maximizes the score */
         if ( score > bestscore )
         {
            bestscore = score;

            *success = TRUE;
            bestvar = var;
            bestnode = v;
            bestvarfixneigh = fixneigh;
         }
      }
   }
   assert( !(*success) || bestvar != NULL );

   if ( *success )
   {
      int* succ;
      int nsucc;
      int s;

      assert( bestnode >= 0 && bestnode < nsos1vars );

      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, bestnode);
      succ = SCIPdigraphGetSuccessors(conflictgraph, bestnode);

      /* if the diving score voted for fixing the best variable to 0.0, we add this as the preferred bound change;
       * otherwise, fixing the neighbors in the conflict graph to 0.0 is the preferred bound change.
       */
      assert( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(bestvar)) || SCIPisFeasPositive(scip, SCIPvarGetUbLocal(bestvar)) );
      SCIP_CALL( SCIPaddDiveBoundChange(scip, bestvar, SCIP_BRANCHDIR_FIXED, 0.0, !bestvarfixneigh) );
      for (s = 0; s < nsucc; ++s)
      {
         SCIP_VAR* var;

         var = SCIPnodeGetVarSOS1(conflictgraph, succ[s]);

         /* if variable is not already fixed */
         if ( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddDiveBoundChange(scip, var, SCIP_BRANCHDIR_FIXED, 0.0, bestvarfixneigh) );
         }
      }
   }

   return SCIP_OKAY;
}


/** determine a diving variables and boundchanges of diving variables by analyzing the SOS1 constraints
 *
 *  if the SOS1 constraints overlap, the method getDiveBdChgsSOS1conflictgraph() may produce better results (e.g., due to more
 *  diving candidates)
 */
static
SCIP_RETCODE getDiveBdChgsSOS1constraints(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            success             /**< pointer to store */
   )
{
   SCIP_VAR* bestvar = NULL;
   SCIP_Bool bestvarfixcomp = FALSE;
   SCIP_Real bestscore = SCIP_REAL_MIN;
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   int nconss;
   int bestcons = -1;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( diveset != NULL );
   assert( success != NULL );

   *success = FALSE;

   /* get SOS1 constraints and number of SOS1 constraints */
   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   /* loop through all SOS1 constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_VAR** vars;
      int nvars;
      int cnt = 0;
      int j;

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      nvars = consdata->nvars;
      vars = consdata->vars;

      /* check whether SOS1 constraint is violated */
      for (j = 0; j < nvars && cnt < 2; ++j)
      {
         SCIP_VAR* var;

         var = vars[j];

         /* check whether variable is nonzero w.r.t. sol and the bounds have not been fixed to zero by propagation */
         if ( !SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, var))
               && (!SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || !SCIPisFeasZero(scip, SCIPvarGetUbLocal(var))) )
            ++cnt;
      }

      /* if SOS1 constraint is not violated then continue with the next SOS1 constraint */
      if ( cnt < 2 )
         continue;

      /* get diving score of every variable in constraint */
      for (j = 0; j < nvars; ++j)
      {
         SCIP_VAR* var;
         SCIP_Real solval;
         SCIP_Real score;
         SCIP_Real bound;
         SCIP_Real fracval;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Bool fixcomp;  /* whether to fix the complementary variables of the candidate in the SOS1 constraint to zero */

         var = vars[j];
         solval = SCIPgetSolVal(scip, sol, var);
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);

         /* check whether variable is nonzero w.r.t. sol and the bounds have not been fixed to zero by propagation */
         if ( ! SCIPisFeasZero(scip, solval) && ( ! SCIPisFeasZero(scip, lb) || ! SCIPisFeasZero(scip, ub) ) )
         {
            /* compute (variable) bound of candidate */
            if ( SCIPisFeasNegative(scip, solval) )
               bound = lb;
            else
               bound = ub;

            /* bound may have changed in propagation; ensure that fracval <= 1 */
            if ( SCIPisFeasLT(scip, REALABS(bound), REALABS(solval)) )
               bound = solval;

            /* ensure finiteness */
            bound = MIN(DIVINGCUTOFFVALUE, REALABS(bound)); /*lint !e666*/
            fracval = MIN(DIVINGCUTOFFVALUE, REALABS(solval)); /*lint !e666*/
            assert( ! SCIPisInfinity(scip, bound) );
            assert( ! SCIPisInfinity(scip, fracval) );
            assert( SCIPisPositive(scip, bound) );

            /* get fractionality of candidate */
            fracval /= (bound + SCIPsumepsilon(scip));

            /* should SOS1 variables be scored by the diving heuristics specific score function;
             *  otherwise use the score function of the SOS1 constraint handler
             */
            if ( SCIPdivesetSupportsType(diveset, SCIP_DIVETYPE_SOS1VARIABLE) )
            {
               SCIP_Bool roundup;

               SCIP_CALL( SCIPgetDivesetScore(scip, diveset, SCIP_DIVETYPE_SOS1VARIABLE, var, solval, fracval,
                &score, &roundup) );

               fixcomp = roundup;
               if ( SCIPisFeasNegative(scip, solval) )
                  fixcomp = !fixcomp;
            }
            else
            {
               /* we always fix the complementary variables of the candidate in the SOS1 constraint to zero */
               fixcomp = TRUE;

               /* score fractionality of candidate */
               score = fracval;
            }

            /* best candidate maximizes the score */
            if ( score > bestscore )
            {
               bestscore = score;

               *success = TRUE;
               bestvar = var;
               bestcons = c;
               bestvarfixcomp = fixcomp;
            }
         }
      }
   }
   assert( !(*success) || bestvar != NULL );

   if ( *success )
   {
      SCIP_VAR** vars;
      int nvars;
      int j;

      consdata = SCIPconsGetData(conss[bestcons]);
      assert( consdata != NULL );

      nvars = consdata->nvars;
      vars = consdata->vars;

      assert( bestcons >= 0 && bestcons < nconss );

      /* if the diving score voted for fixing the best variable to 0.0, we add this as the preferred bound change;
       * otherwise, fixing the complementary variables of the candidate in the SOS1 constraint to 0.0 is the preferred bound change.
       */
      assert( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(bestvar)) || SCIPisFeasPositive(scip, SCIPvarGetUbLocal(bestvar)) );

      SCIP_CALL( SCIPaddDiveBoundChange(scip, bestvar, SCIP_BRANCHDIR_FIXED, 0.0, !bestvarfixcomp) );
      for (j = 0; j < nvars; ++j)
      {
         SCIP_VAR* var;

         var = vars[j];

         /* if variable is not already fixed and is not the candidate variable */
         if ( var != bestvar && ( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var)) ) )
         {
            SCIP_CALL( SCIPaddDiveBoundChange(scip, var, SCIP_BRANCHDIR_FIXED, 0.0, bestvarfixcomp) );
         }
      }
   }

   return SCIP_OKAY;
}


/* --------------------initialization/deinitialization ------------------------*/

/** check whether \f$x_1\f$ is a bound variable of \f$x_0\f$; i.e., \f$x_0 \leq c\cdot x_1\f$ or \f$x_0 \geq d\cdot x_1\f$
 *  for positive values \f$c, d\f$. If true, then add this information to the node data of the conflict graph.
 */
static
SCIP_RETCODE detectVarboundSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler data */
   SCIP_VAR*             var0,               /**< first variable */
   SCIP_VAR*             var1,               /**< second variable */
   SCIP_Real             val0,               /**< first coefficient */
   SCIP_Real             val1                /**< second coefficient */
   )
{
   int node0;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( var0 != NULL && var1 != NULL );

   /* get nodes of variable in the conflict graph (node = -1 if no SOS1 variable) */
   node0 = varGetNodeSOS1(conshdlrdata, var0);

   /* if var0 is an SOS1 variable */
   if ( node0 >= 0 )
   {
      SCIP_Real val;

      assert( ! SCIPisFeasZero(scip, val0) );
      val = -val1/val0;

      /* check variable bound relation of variables */

      /* handle lower bound case */
      if ( SCIPisFeasNegative(scip, val0) && SCIPisFeasNegative(scip, val) )
      {
         SCIP_NODEDATA* nodedata;

         /* get node data of the conflict graph */
         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conshdlrdata->conflictgraph, node0);

         /* @todo: maybe save multiple variable bounds for each SOS1 variable */
         if ( nodedata->lbboundvar == NULL )
         {
            /* add variable bound information to node data */
            nodedata->lbboundvar = var1;
            nodedata->lbboundcoef = val;

            SCIPdebugMsg(scip, "detected variable bound constraint %s >= %f %s.\n", SCIPvarGetName(var0), val, SCIPvarGetName(var1));
         }
      }
      /* handle upper bound case */
      else if ( SCIPisFeasPositive(scip, val0) && SCIPisFeasPositive(scip, val) )
      {
         SCIP_NODEDATA* nodedata;
         assert( SCIPisFeasPositive(scip, val0) );

         /* get node data of the conflict graph */
         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conshdlrdata->conflictgraph, node0);

         if ( nodedata->ubboundvar == NULL )
         {
            /* add variable bound information to node data */
            nodedata->ubboundvar = var1;
            nodedata->ubboundcoef = val;

            SCIPdebugMsg(scip, "detected variable bound constraint %s <= %f %s.\n", SCIPvarGetName(var0), val, SCIPvarGetName(var1));
         }
      }
   }

   return SCIP_OKAY;
}


/** pass connected component \f$C\f$ of the conflict graph and check whether all the variables correspond to a unique variable upper bound variable \f$z\f$,
 *  i.e., \f$x_i \leq u_i z\f$ for every \f$i\in C\f$.
 *
 *  @note if the bound variable is unique, then bound inequalities can be strengthened.
 */
static
SCIP_RETCODE passConComponentVarbound(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   node,               /**< current node of connected component */
   SCIP_VAR*             boundvar,           /**< bound variable of connected component */
   SCIP_Bool             checklb,            /**< whether to check lower bound variable (else upper bound variable) */
   SCIP_Bool*            processed,          /**< states for each variable whether it has been processed */
   int*                  concomp,            /**< current connected component */
   int*                  nconcomp,           /**< pointer to store number of elements of connected component */
   SCIP_Bool*            unique              /**< pointer to store whether bound variable is unique */
   )
{
   int* succ;
   int nsucc;
   int s;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( processed != NULL );
   assert( concomp != NULL );
   assert( nconcomp != NULL );
   assert( unique != NULL );

   processed[node] = TRUE;/*lint !e737*/
   concomp[(*nconcomp)++] = node;

   /* if bound variable of connected component without new node is unique */
   if ( *unique )
   {
      SCIP_NODEDATA* nodedata;
      SCIP_VAR* comparevar;
      nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, node);
      assert( nodedata != NULL );

      if ( checklb )
         comparevar = nodedata->lbboundvar;
      else
         comparevar = nodedata->ubboundvar;

      /* check whether bound variable is unique for connected component without new node */
      if ( boundvar == NULL )
      {
         if ( comparevar != NULL )
            *unique = FALSE;
      }
      else
      {
         if ( comparevar == NULL )
            *unique = FALSE;
         else if ( SCIPvarCompare(boundvar, comparevar) != 0 )
            *unique = FALSE;
      }
   }

   /* pass through successor variables */
   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, node);
   succ = SCIPdigraphGetSuccessors(conflictgraph, node);
   for (s = 0; s < nsucc; ++s)
   {
      if ( ! processed[succ[s]] )
         SCIP_CALL( passConComponentVarbound(scip, conflictgraph, succ[s], boundvar, checklb, processed, concomp, nconcomp, unique) );
   }

   return SCIP_OKAY;
}


/** for each connected component \f$C\f$ of the conflict graph check whether all the variables correspond to a unique variable upper bound variable \f$z\f$
 *  (e.g., for the upper bound case this means that \f$x_i \leq u_i z\f$ for every \f$i\in C\f$).
 *
 *  @note if the bound variable is unique, then bound inequalities can be strengthened.
 */
static
SCIP_RETCODE checkConComponentsVarbound(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_Bool             checklb             /**< whether to check lower bound variable (else check upper bound variable) */
   )
{
   SCIP_Bool* processed;  /* states for each variable whether it has been processed */
   int* concomp;          /* current connected component */
   int nconcomp;
   int j;

   assert( scip != NULL );
   assert( conflictgraph != NULL );

   /* allocate buffer arrays and initialize 'processed' array */
   SCIP_CALL( SCIPallocBufferArray(scip, &processed, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &concomp, nsos1vars) );
   for (j = 0; j < nsos1vars; ++j)
      processed[j] = FALSE;

   /* run through all SOS1 variables */
   for (j = 0; j < nsos1vars; ++j)
   {
      /* if variable belongs to a connected component that has not been processed so far */
      if ( ! processed[j] )
      {
         SCIP_NODEDATA* nodedata;
         SCIP_VAR* boundvar;
         SCIP_Bool unique;
         int* succ;
         int nsucc;
         int s;

         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, j);
         assert( nodedata != NULL );

         if ( checklb )
            boundvar = nodedata->lbboundvar;
         else
            boundvar = nodedata->ubboundvar;
         unique = TRUE;

         processed[j] = TRUE;
         concomp[0] = j;
         nconcomp = 1;

         /* pass through successor variables */
         nsucc = SCIPdigraphGetNSuccessors(conflictgraph, j);
         succ = SCIPdigraphGetSuccessors(conflictgraph, j);
         for (s = 0; s < nsucc; ++s)
         {
            if ( ! processed[succ[s]] )
            {
               SCIP_CALL( passConComponentVarbound(scip, conflictgraph, succ[s], boundvar, checklb, processed, concomp, &nconcomp, &unique) );
            }
         }

         /* if the connected component has a unique bound variable */
         if ( unique && boundvar != NULL )
         {
            for (s = 0; s < nconcomp; ++s)
            {
               nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, concomp[s]);
               assert( processed[concomp[s]] == TRUE );
               assert( nodedata != NULL );

               if ( checklb )
                  nodedata->lbboundcomp = TRUE;
               else
                  nodedata->ubboundcomp = TRUE;
            }
            SCIPdebugMsg(scip, "Found a connected component of size <%i> with unique bound variable.\n", nconcomp);
         }
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &concomp);
   SCIPfreeBufferArray(scip, &processed);

   return SCIP_OKAY;
}


/** check all linear constraints for variable bound constraints of the form \f$c\cdot z \leq x \leq d\cdot z\f$, where @p x is some SOS1
 *  variable and @p z is some arbitrary variable (not necessarily binary)
 */
static
SCIP_RETCODE checkLinearConssVarboundSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler data */
   SCIP_CONS**           linconss,           /**< linear constraints */
   int                   nlinconss           /**< number of linear constraints */
   )
{
   int c;

   /* loop through linear constraints */
   for (c = 0; c < nlinconss; ++c)
   {
      SCIP_CONS* lincons;
      int nvars;

      lincons = linconss[c];

      /* variable bound constraints only contain two variables */
      nvars = SCIPgetNVarsLinear(scip, lincons);
      if ( nvars == 2 )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         SCIP_VAR* var0;
         SCIP_VAR* var1;
         SCIP_Real lhs;
         SCIP_Real rhs;

         /* get constraint data */
         vars = SCIPgetVarsLinear(scip, lincons);
         vals = SCIPgetValsLinear(scip, lincons);
         lhs = SCIPgetLhsLinear(scip, lincons);
         rhs = SCIPgetRhsLinear(scip, lincons);

         var0 = vars[0];
         var1 = vars[1];
         assert( var0 != NULL && var1 != NULL );

         /* at least one variable should be an SOS1 variable */
         if ( varIsSOS1(conshdlrdata, var0) || varIsSOS1(conshdlrdata, var1) )
         {
            SCIP_Real val0;
            SCIP_Real val1;

            /* check whether right hand side or left hand side of constraint is zero */
            if ( SCIPisFeasZero(scip, lhs) )
            {
               val0 = -vals[0];
               val1 = -vals[1];

               /* check whether the two variables are in a variable bound relation */
               SCIP_CALL( detectVarboundSOS1(scip, conshdlrdata, var0, var1, val0, val1) );
               SCIP_CALL( detectVarboundSOS1(scip, conshdlrdata, var1, var0, val1, val0) );
            }
            else if( SCIPisFeasZero(scip, rhs) )
            {
               val0 = vals[0];
               val1 = vals[1];

               /* check whether the two variables are in a variable bound relation */
               SCIP_CALL( detectVarboundSOS1(scip, conshdlrdata, var0, var1, val0, val1) );
               SCIP_CALL( detectVarboundSOS1(scip, conshdlrdata, var1, var0, val1, val0) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** switch to SOS1 branching and separating bound iniqualities from SOS1 constraints if the SOS1 constraints do not overlap */
static
SCIP_RETCODE checkSwitchNonoverlappingSOS1Methods(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   int                   nconss              /**< number of SOS1 constraints */
   )
{
   SCIP_Bool nonoverlap = TRUE;
   int c;

   /* loop through all SOS1 constraints */
   if ( conshdlrdata->nsos1vars > 0 )
   {
      for (c = 0; c < nconss && nonoverlap; ++c)
      {
         SCIP_CONSDATA* consdata;
         SCIP_VAR** vars;
         int notfixed = 0;
         int nvars;
         int i;

         assert( conss[c] != NULL );

         /* get constraint data field of the constraint */
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );

         /* get variables and number of variables of constraint */
         nvars = consdata->nvars;
         vars = consdata->vars;

         /* get number of variables of SOS1 constraint that are not fixed to zero */
         for (i = 0; i < nvars; ++i)
         {
            if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(vars[i])) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i])) )
               ++notfixed;
         }

         /* check variables of SOS1 constraint */
         for (i = 0; i < nvars; ++i)
         {
            int node;

            assert( vars[i] != NULL );

            node = varGetNodeSOS1(conshdlrdata, vars[i]);
            assert( node >= 0 || ( SCIPisFeasZero(scip, SCIPvarGetLbLocal(vars[i])) && SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i]))) );
            assert( node < conshdlrdata->nsos1vars );
            assert( node < 0 || SCIPdigraphGetNSuccessors(conflictgraph, node) >= notfixed-1 );
            if ( node >= 0 && SCIPdigraphGetNSuccessors(conflictgraph, node) > notfixed-1 )
            {
               nonoverlap = FALSE;
               break;
            }
         }
      }
   }

   /* if the SOS1 constraints do not overlap */
   if ( nonoverlap )
   {
      if ( conshdlrdata->autosos1branch )
      {
         conshdlrdata->switchsos1branch = TRUE;
         SCIPdebugMsg(scip, "Switched to SOS1 branching, since the SOS1 constraints do not overlap\n");
      }

      if ( conshdlrdata->autocutsfromsos1 )
      {
         conshdlrdata->switchcutsfromsos1 = TRUE;
         SCIPdebugMsg(scip, "Switched to separating bound cuts from SOS1 constraints (and not from the conflict graph), since the SOS1 constraints do not overlap\n");
      }
   }

   return SCIP_OKAY;
}


/** sets node data of conflict graph nodes */
static
SCIP_RETCODE computeNodeDataSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler data */
   int                   nsos1vars           /**< number of SOS1 variables */
   )
{
   SCIP_CONSHDLR* linconshdlr;
   SCIP_CONS** linconss;
   int nlinconss;

   /* if no SOS1 variables exist -> exit */
   if ( nsos1vars == 0 )
      return SCIP_OKAY;

   /* get constraint handler data of linear constraints */
   linconshdlr = SCIPfindConshdlr(scip, "linear");
   if ( linconshdlr == NULL )
      return SCIP_OKAY;

   /* get linear constraints and number of linear constraints */
   nlinconss = SCIPconshdlrGetNConss(linconshdlr);
   linconss = SCIPconshdlrGetConss(linconshdlr);

   /* check linear constraints for variable bound constraints */
   SCIP_CALL( checkLinearConssVarboundSOS1(scip, conshdlrdata, linconss, nlinconss) );

   /* for each connected component of the conflict graph check whether all the variables correspond to a unique variable
    * upper bound variable */
   SCIP_CALL( checkConComponentsVarbound(scip, conshdlrdata->conflictgraph, conshdlrdata->nsos1vars, TRUE) );
   SCIP_CALL( checkConComponentsVarbound(scip, conshdlrdata->conflictgraph, conshdlrdata->nsos1vars, FALSE) );

   return SCIP_OKAY;
}


/** initialize conflictgraph and create hashmap for SOS1 variables */
static
SCIP_RETCODE initConflictgraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   int                   nconss              /**< number of SOS1 constraints */
   )
{
   SCIP_Bool* nodecreated; /* nodecreated[i] = TRUE if a node in the conflict graph is already created for index i
                            * (with i index of the original variables) */
   int* nodeorig;          /* nodeorig[i] = node of original variable x_i in the conflict graph */
   int ntotalvars;
   int cntsos;
   int i;
   int j;
   int c;

   assert( conshdlrdata != NULL );
   assert( nconss == 0 || conss != NULL );

   /* get the number of original problem variables */
   ntotalvars = SCIPgetNTotalVars(scip);

   /* initialize vector 'nodecreated' */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeorig, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodecreated, ntotalvars) );
   for (i = 0; i < ntotalvars; ++i)
      nodecreated[i] = FALSE;

   /* compute number of SOS1 variables */
   cntsos = 0;
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR** vars;
      int nvars;

      assert( conss[c] != NULL );

      /* get constraint data field of the constraint */
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* get variables and number of variables of constraint */
      nvars = consdata->nvars;
      vars = consdata->vars;

      /* update number of SOS1 variables */
      for (i = 0; i < nvars; ++i)
      {
         SCIP_VAR* var;

         var = vars[i];

         /* if the variable is not fixed to zero */
         if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
         {
            int ind;

            ind = SCIPvarGetIndex(var);
            assert( ind >= 0 && ind < ntotalvars );
            if ( ! nodecreated[ind] )
            {
               nodecreated[ind] = TRUE; /* mark node as counted */
               nodeorig[ind] = cntsos;
               ++cntsos;
            }
         }
      }
   }
   if ( cntsos <= 0 )
   {
      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &nodecreated);
      SCIPfreeBufferArray(scip, &nodeorig);
      conshdlrdata->nsos1vars = 0;
      return SCIP_OKAY;
   }

   /* reinitialize vector 'nodecreated' */
   for (i = 0; i < ntotalvars; ++i)
      nodecreated[i] = FALSE;

   /* create conflict graph */
   SCIP_CALL( SCIPdigraphCreate(&conshdlrdata->conflictgraph, SCIPblkmem(scip), cntsos) );

   /* set up hash map */
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->varhash, SCIPblkmem(scip), cntsos) );

   /* for every SOS1 constraint */
   cntsos = 0;
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR** vars;
      int nvars;

      assert( conss[c] != NULL );

      /* get constraint data field of the constraint */
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* get variables and number of variables of constraint */
      nvars = consdata->nvars;
      vars = consdata->vars;

      /* add edges to the conflict graph and create node data for each of its nodes */
      for (i = 0; i < nvars; ++i)
      {
         SCIP_VAR* var;

         var = vars[i];

         /* if the variable is not fixed to zero */
         if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
         {
            int indi;

            indi = SCIPvarGetIndex(var);

            if ( ! nodecreated[indi] )
            {
               SCIP_NODEDATA* nodedata = NULL;

               /* insert node number to hash map */
               assert( ! SCIPhashmapExists(conshdlrdata->varhash, var) );
               SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varhash, var, (void*) (size_t) cntsos) );/*lint !e571*/
               assert( cntsos == (int) (size_t) SCIPhashmapGetImage(conshdlrdata->varhash, var) );
               assert( SCIPhashmapExists(conshdlrdata->varhash, var) );

               /* create node data */
               SCIP_CALL( SCIPallocBlockMemory(scip, &nodedata) );
               nodedata->var = var;
               nodedata->lbboundvar = NULL;
               nodedata->ubboundvar = NULL;
               nodedata->lbboundcoef = 0.0;
               nodedata->ubboundcoef = 0.0;
               nodedata->lbboundcomp = FALSE;
               nodedata->ubboundcomp = FALSE;

               /* set node data */
               SCIPdigraphSetNodeData(conshdlrdata->conflictgraph, (void*)nodedata, cntsos);

               /* mark node and var data of node as created and update SOS1 counter */
               nodecreated[indi] = TRUE;
               ++cntsos;
            }

            /* add edges to the conflict graph */
            for (j = i+1; j < nvars; ++j)
            {
               var = vars[j];

               /* if the variable is not fixed to zero */
               if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
               {
                  int indj;

                  indj = SCIPvarGetIndex(var);

                  /* in case indi = indj the variable will be deleted in the presolving step */
                  if ( indi != indj )
                  {
                     /* arcs have to be added 'safe' */
                     SCIP_CALL( SCIPdigraphAddArcSafe(conshdlrdata->conflictgraph, nodeorig[indi], nodeorig[indj], NULL) );
                     SCIP_CALL( SCIPdigraphAddArcSafe(conshdlrdata->conflictgraph, nodeorig[indj], nodeorig[indi], NULL) );
                  }
               }
            }
         }
      }
   }

   /* set number of problem variables that are contained in at least one SOS1 constraint */
   conshdlrdata->nsos1vars = cntsos;

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &nodecreated);
   SCIPfreeBufferArray(scip, &nodeorig);

   /* sort successors in ascending order */
   for (j = 0; j < conshdlrdata->nsos1vars; ++j)
   {
      int nsucc;

      nsucc = SCIPdigraphGetNSuccessors(conshdlrdata->conflictgraph, j);
      SCIPsortInt(SCIPdigraphGetSuccessors(conshdlrdata->conflictgraph, j), nsucc);
   }

   return SCIP_OKAY;
}


/** free conflict graph, nodedata and hashmap */
static
SCIP_RETCODE freeConflictgraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   int j;

   if ( conshdlrdata->conflictgraph == NULL )
   {
      assert( conshdlrdata->nsos1vars == 0 );
      return SCIP_OKAY;
   }

   /* for every SOS1 variable */
   assert( conshdlrdata->nsos1vars > 0 );
   for (j = 0; j < conshdlrdata->nsos1vars; ++j)
   {
      SCIP_NODEDATA* nodedata;

      /* get node data */
      assert( conshdlrdata->conflictgraph != NULL );
      nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conshdlrdata->conflictgraph, j);
      assert( nodedata != NULL );

      /* free node data */
      SCIPfreeBlockMemory(scip, &nodedata);
      SCIPdigraphSetNodeData(conshdlrdata->conflictgraph, NULL, j);
   }

   /* free conflict graph and hash map */
   assert( conshdlrdata->varhash != NULL );
   SCIPhashmapFree(&conshdlrdata->varhash);
   SCIPdigraphFree(&conshdlrdata->conflictgraph);
   conshdlrdata->nsos1vars = 0;

   assert( conshdlrdata->varhash == NULL );
   assert( conshdlrdata->conflictgraph == NULL );

   return SCIP_OKAY;
}


/* ---------------------------- constraint handler callback methods ----------------------*/

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSOS1(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSOS1)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free stack of variables fixed to nonzero (usually already freed in consExitsolSOS1 unless instance was solved during presolving) */
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->fixnonzerovars, conshdlrdata->maxnfixnonzerovars); /*lint !e737*/

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolSOS1)
{  /*lint --e{715}*/
    SCIP_CONSHDLRDATA* conshdlrdata;

    assert( scip != NULL );
    assert( conshdlr != NULL );
    assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert( conshdlrdata != NULL );

    conshdlrdata->nsos1vars = 0;
    conshdlrdata->varhash = NULL;

    if ( nconss > 0 )
    {
       /* initialize conflict graph and hashmap for SOS1 variables */
       SCIP_CALL( initConflictgraph(scip, conshdlrdata, conss, nconss) );

       /* add data to conflict graph nodes */
       SCIP_CALL( computeNodeDataSOS1(scip, conshdlrdata, conshdlrdata->nsos1vars) );

       if ( ( conshdlrdata->autosos1branch || conshdlrdata->autocutsfromsos1 )
          && ( ! conshdlrdata->switchsos1branch || ! conshdlrdata->switchcutsfromsos1 )
          )
       {
          /* switch to nonoverlapping methods if the SOS1 constraints do not overlap */
          SCIP_CALL( checkSwitchNonoverlappingSOS1Methods(scip, conshdlrdata, conshdlrdata->conflictgraph, conss, nconss) );
       }

       /* initialize tclique graph */
       SCIP_CALL( initTCliquegraph(scip, conshdlr, conshdlrdata, conshdlrdata->conflictgraph, conshdlrdata->nsos1vars) );

       /* create local conflict graph if needed */
       if ( conshdlrdata->addcomps )
       {
          SCIP_CALL( SCIPdigraphCreate(&conshdlrdata->localconflicts, SCIPblkmem(scip), conshdlrdata->nsos1vars) );
       }

       /* initialize stack of variables fixed to nonzero (memory may be already allocated in consTransSOS1()) */
       if ( conshdlrdata->fixnonzerovars == NULL )
       {
          conshdlrdata->maxnfixnonzerovars = conshdlrdata->nsos1vars;
          SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->fixnonzerovars, conshdlrdata->maxnfixnonzerovars) );
       }
    }

    return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolSOS1)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      SCIPdebugMsg(scip, "Exiting SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* free rows */
      if ( consdata->rowub != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rowub) );
      }

      if ( consdata->rowlb != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rowlb) );
      }
   }

   /* free implication graph */
   if ( conshdlrdata->implgraph != NULL )
   {
      SCIP_CALL( freeImplGraphSOS1(scip, conshdlrdata) );
   }
   assert( conshdlrdata->implgraph == NULL );

   /* free tclique graph and tclique data */
   if ( conshdlrdata->tcliquegraph != NULL )
   {
      assert( conshdlrdata->tcliquedata != NULL );
      SCIPfreeBlockMemory(scip, &conshdlrdata->tcliquedata);
      tcliqueFree(&conshdlrdata->tcliquegraph);
   }
   assert(conshdlrdata->tcliquegraph == NULL);
   assert(conshdlrdata->tcliquedata == NULL);

   /* free stack of variables fixed to nonzero */
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->fixnonzerovars, conshdlrdata->maxnfixnonzerovars); /*lint !e737*/
   conshdlrdata->nfixnonzerovars = 0;
   conshdlrdata->maxnfixnonzerovars = 0;

   /* free graph for storing local conflicts */
   if ( conshdlrdata->localconflicts != NULL )
      SCIPdigraphFree(&conshdlrdata->localconflicts);
   assert( conshdlrdata->localconflicts == NULL );

   /* free conflict graph  */
   SCIP_CALL( freeConflictgraph(scip, conshdlrdata) );
   assert( conshdlrdata->conflictgraph == NULL );

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSOS1)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMsg(scip, "Deleting SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

   /* drop events on transformed variables */
   if ( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int j;

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      assert( conshdlrdata->eventhdlr != NULL );

      for (j = 0; j < (*consdata)->nvars; ++j)
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)cons, -1) ); /*lint !e737 !e740*/
      }
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->maxvars);
   if ( (*consdata)->weights != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, (*consdata)->maxvars);
   }

   /* free rows */
   if ( (*consdata)->rowub != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->rowub) );
   }
   if ( (*consdata)->rowlb != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->rowlb) );
   }
   assert( (*consdata)->rowub == NULL );
   assert( (*consdata)->rowlb == NULL );

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSOS1)
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   char s[SCIP_MAXSTRLEN];
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( conshdlrdata->eventhdlr != NULL );

   SCIPdebugMsg(scip, "Transforming SOS1 constraint: <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->nvars > 0 );
   assert( sourcedata->nvars <= sourcedata->maxvars );

   /* initialize stack of variables fixed to nonzero */
   if ( conshdlrdata->fixnonzerovars == NULL )
   {
      conshdlrdata->maxnfixnonzerovars = SCIPgetNTotalVars(scip);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->fixnonzerovars, conshdlrdata->maxnfixnonzerovars) );
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nvars = sourcedata->nvars;
   consdata->maxvars = sourcedata->nvars;
   consdata->rowub = NULL;
   consdata->rowlb = NULL;
   consdata->nfixednonzeros = 0;
   consdata->local = sourcedata->local;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, consdata->nvars) );
   /* if weights were used */
   if ( sourcedata->weights != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->weights, sourcedata->weights, consdata->nvars) );
   }
   else
      consdata->weights = NULL;

   for (j = 0; j < sourcedata->nvars; ++j)
   {
      assert( sourcedata->vars[j] != 0 );
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->vars[j], &(consdata->vars[j])) );

      /* if variable is fixed to be nonzero */
      if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(consdata->vars[j])) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(consdata->vars[j])) )
         ++(consdata->nfixednonzeros);
   }

   /* create transformed constraint with the same flags */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch bound change events on variable */
   for (j = 0; j < consdata->nvars; ++j)
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
            (SCIP_EVENTDATA*)*targetcons, NULL) ); /*lint !e740*/
   }

#ifdef SCIP_DEBUG
   if ( consdata->nfixednonzeros > 0 )
   {
      SCIPdebugMsg(scip, "constraint <%s> has %d variables fixed to be nonzero.\n", SCIPconsGetName(*targetcons),
         consdata->nfixednonzeros );
   }
#endif

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSOS1)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   /* cppcheck-suppress unassignedVariable */
   int oldnfixedvars;
   /* cppcheck-suppress unassignedVariable */
   int oldnchgbds;
   /* cppcheck-suppress unassignedVariable */
   int oldndelconss;
   /* cppcheck-suppress unassignedVariable */
   int oldnupgdconss;
   int nremovedvars;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPdebugMsg(scip, "Presolving SOS1 constraints.\n");

   *result = SCIP_DIDNOTRUN;

   SCIPdebug( oldnfixedvars = *nfixedvars; )
   SCIPdebug( oldnchgbds = *nchgbds; )
   SCIPdebug( oldndelconss = *ndelconss; )
   SCIPdebug( oldnupgdconss = *nupgdconss; )
   nremovedvars = 0;

   /* only run if success if possible */
   if( nconss > 0 && ( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 || nnewchgbds > 0 ) )
   {
      SCIP_Bool** adjacencymatrix = NULL;
      SCIP_DIGRAPH* conflictgraph;
      SCIP_EVENTHDLR* eventhdlr;
      int nsos1vars;
      int i;
      int j;

      *result = SCIP_DIDNOTFIND;

      /* get constraint handler data */
      assert( SCIPconshdlrGetData(conshdlr) != NULL );
      eventhdlr = SCIPconshdlrGetData(conshdlr)->eventhdlr;
      assert( eventhdlr != NULL );

      /* initialize conflict graph */
      SCIP_CALL( initConflictgraph(scip, conshdlrdata, conss, nconss));

      /* get conflict graph and number of SOS1 variables */
      conflictgraph = conshdlrdata->conflictgraph;
      nsos1vars = conshdlrdata->nsos1vars;
      if ( nsos1vars < 2 )
      {
         SCIP_CALL( freeConflictgraph(scip, conshdlrdata));
         return SCIP_OKAY;
      }

      /* we do not create the adjacency matrix of the conflict graph if the number of SOS1 variables is larger than a predefined value */
      if ( conshdlrdata->maxsosadjacency == -1 || nsos1vars <= conshdlrdata->maxsosadjacency )
      {
         /* allocate buffer arrays for adjacency matrix  */
         SCIP_CALL( SCIPallocBufferArray(scip, &adjacencymatrix, nsos1vars) );
         for (i = 0; i < nsos1vars; ++i)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &adjacencymatrix[i], i+1) );/*lint !e866*/
         }

         /* create adjacency matrix */
         for (i = 0; i < nsos1vars; ++i)
         {
            for (j = 0; j < i+1; ++j)
               adjacencymatrix[i][j] = 0;
         }
         for (i = 0; i < nsos1vars; ++i)
         {
            int* succ;
            int nsucc;

            succ = SCIPdigraphGetSuccessors(conflictgraph, i);
            nsucc = SCIPdigraphGetNSuccessors(conflictgraph, i);

            for (j = 0; j < nsucc; ++j)
            {
               if ( i > succ[j] )
                  adjacencymatrix[i][succ[j]] = 1;
            }
         }
      }
      else
      {
         SCIPdebugMsg(scip, "Adjacency matrix was not created since number of SOS1 variables (%d) is larger than %d.\n", nsos1vars, conshdlrdata->maxsosadjacency);
      }

      /* perform one presolving round for SOS1 constraints */
      SCIP_CALL( presolRoundConssSOS1(scip, eventhdlr, conshdlrdata, conflictgraph, adjacencymatrix, conss, nconss, nsos1vars, naddconss, ndelconss, nupgdconss, nfixedvars, &nremovedvars, result) );

      if ( adjacencymatrix != NULL )
      {
         /* perform one presolving round for SOS1 variables */
         if ( conshdlrdata->maxtightenbds != 0 && *result != SCIP_CUTOFF )
         {
            SCIP_CALL( presolRoundVarsSOS1(scip, conshdlrdata, conflictgraph, adjacencymatrix, nsos1vars, nfixedvars, nchgbds, naddconss, result) );
         }

         /* free adjacency matrix */
         for (j = nsos1vars-1; j >= 0; --j)
            SCIPfreeBufferArrayNull(scip, &adjacencymatrix[j]);
         SCIPfreeBufferArrayNull(scip, &adjacencymatrix);
      }

      /* free memory allocated in function initConflictgraph() */
      SCIP_CALL( freeConflictgraph(scip, conshdlrdata));
   }
   (*nchgcoefs) += nremovedvars;

   SCIPdebugMsg(scip, "presolving fixed %d variables, changed %d bounds, removed %d variables, deleted %d constraints, and upgraded %d constraints.\n",
                *nfixedvars - oldnfixedvars, *nchgbds - oldnchgbds, nremovedvars, *ndelconss - oldndelconss, *nupgdconss - oldnupgdconss);

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpSOS1)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   *infeasible = FALSE;

   /* checking for initial rows for SOS1 constraints */
   if( conshdlrdata->boundcutsfromsos1 || conshdlrdata->switchcutsfromsos1 )
   {
      SCIP_CALL( initsepaBoundInequalityFromSOS1Cons(scip, conshdlr, conshdlrdata, conss, nconss, NULL, FALSE, -1, NULL, infeasible) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpSOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( separateSOS1(scip, conshdlr, NULL, nconss, conss, result) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolSOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( separateSOS1(scip, conshdlr, sol, nconss, conss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceSOS1(scip, conshdlr, nconss, conss, NULL, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceSOS1(scip, conshdlr, nconss, conss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceSOS1(scip, conshdlr, nconss, conss, NULL, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions
 *
 *  We simply check whether at most one variable is nonzero in the given solution.
 */
static
SCIP_DECL_CONSCHECK(consCheckSOS1)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* check each constraint */
   for (c = 0; c < nconss && (*result == SCIP_FEASIBLE || completely); ++c)
   {
      SCIP_CONSDATA* consdata;
      int j;
      int cnt;

      cnt = 0;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMsg(scip, "Checking SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]));

      /* check all variables */
      for (j = 0; j < consdata->nvars; ++j)
      {
         /* if variable is nonzero */
         if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[j])) )
         {
            ++cnt;

            /* if more than one variable is nonzero */
            if ( cnt > 1 )
            {
               SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
               *result = SCIP_INFEASIBLE;

               /* update constraint violation in solution */
               if ( sol != NULL )
                  SCIPupdateSolConsViolation(scip, sol, 1.0, 1.0);

               if ( printreason )
               {
                  int l;

                  SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
                  SCIPinfoMessage(scip, NULL, ";\nviolation: ");

                  for (l = 0; l < consdata->nvars; ++l)
                  {
                     /* if variable is nonzero */
                     if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[l])) )
                     {
                        SCIPinfoMessage(scip, NULL, "<%s> = %.15g ",
                           SCIPvarGetName(consdata->vars[l]), SCIPgetSolVal(scip, sol, consdata->vars[l]));
                     }
                  }
                  SCIPinfoMessage(scip, NULL, "\n");
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSOS1)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_DIGRAPH* conflictgraph;
   SCIP_DIGRAPH* implgraph;
   int ngen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );
   assert( SCIPisTransformed(scip) );

   /* return if number of SOS1 constraints is zero */
   if ( nconss < 1 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   *result = SCIP_DIDNOTFIND;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* get conflict graph */
   conflictgraph = conshdlrdata->conflictgraph;

   /* get/initialize implication graph */
   implgraph = conshdlrdata->implgraph;
   if ( implgraph == NULL && conshdlrdata->implprop && conflictgraph != NULL )
   {
      if ( SCIPgetDepth(scip) == 0 )
      {
         SCIP_Bool success;
         SCIP_Bool cutoff;
         int nchbds;

         SCIP_CALL( initImplGraphSOS1(scip, conshdlrdata, conflictgraph, conshdlrdata->nsos1vars, conshdlrdata->maxtightenbds, &nchbds, &cutoff, &success) );
         if ( ! success )
            conshdlrdata->implprop = FALSE;

         if ( cutoff )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( nchbds > 0 )
            *result = SCIP_REDUCEDDOM;
         implgraph = conshdlrdata->implgraph;
      }
      else
         conshdlrdata->implprop = FALSE;
   }

   /* if conflict graph propagation shall be used */
   if ( conshdlrdata->conflictprop && conflictgraph != NULL )
   {
      SCIP_VAR** fixnonzerovars;
      int nfixnonzerovars;
      int j;

      assert( nconss > 0 );

      /* stack of variables fixed to nonzero */
      nfixnonzerovars = conshdlrdata->nfixnonzerovars;
      fixnonzerovars = conshdlrdata->fixnonzerovars;
      assert( fixnonzerovars != NULL );

      /* check each variable from stack */
      for (j = 0; j < nfixnonzerovars; ++j)
      {
         SCIP_VAR* var;

         var = fixnonzerovars[j];
         if ( var != NULL )
         {
            int node;
            node = varGetNodeSOS1(conshdlrdata, var);

            /* if variable is involved in an SOS1 constraint */
            if ( node >= 0 )
            {
               assert( varGetNodeSOS1(conshdlrdata, var) < conshdlrdata->nsos1vars );
               SCIPdebugMsg(scip, "Propagating SOS1 variable <%s>.\n", SCIPvarGetName(var) );

               /* if zero is outside the domain of variable */
               if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
               {
                  SCIP_Bool cutoff;

                  SCIP_CALL( propVariableNonzero(scip, conflictgraph, implgraph, conss[0], node, conshdlrdata->implprop, &cutoff, &ngen) );
                  if ( cutoff )
                  {
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
               }
            }
         }
      }
   }
   conshdlrdata->nfixnonzerovars = 0;

   /* if SOS1 constraint propagation shall be used */
   if ( conshdlrdata->sosconsprop || conflictgraph == NULL )
   {
      int c;

      /* check each constraint */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_CONS* cons;
         SCIP_CONSDATA* consdata;
         SCIP_Bool cutoff;

         assert( conss[c] != NULL );
         cons = conss[c];
         consdata = SCIPconsGetData(cons);
         assert( consdata != NULL );
         SCIPdebugMsg(scip, "Propagating SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

         SCIP_CALL( propConsSOS1(scip, cons, consdata, &cutoff, &ngen) );
         if ( cutoff )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      }
   }

   SCIPdebugMsg(scip, "Propagated %d domains.\n", ngen);
   if ( ngen > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler
 *
 *  We check which bound changes were the reason for infeasibility. We
 *  use that @a inferinfo stores the index of the variable that has
 *  bounds that fix it to be nonzero (these bounds are the reason). */
static
SCIP_DECL_CONSRESPROP(consRespropSOS1)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;
   SCIPdebugMsg(scip, "Propagation resolution method of SOS1 constraint <%s>.\n", SCIPconsGetName(cons));

   /* check whether conflict was detected in variable propagation or constraint propagation */
   if ( inferinfo < 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      assert( conshdlr != NULL );

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      assert( conshdlrdata->conflictgraph != NULL );
      assert( inferinfo >= -conshdlrdata->maxnfixnonzerovars );
      assert( inferinfo >= -conshdlrdata->nsos1vars );
      assert( inferinfo <= -1 );

      var = SCIPnodeGetVarSOS1(conshdlrdata->conflictgraph, -inferinfo - 1);
   }
   else
   {
      SCIP_CONSDATA* consdata;

      /* get constraint data */
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( inferinfo < consdata->nvars );

      var = consdata->vars[inferinfo];
   }
   assert( var != NULL );
   assert( var != infervar );

   /* check if lower bound of var was the reason */
   if ( SCIPisFeasPositive(scip, SCIPgetVarLbAtIndex(scip, var, bdchgidx, FALSE)) )
   {
      SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   /* check if upper bound of var was the reason */
   if ( SCIPisFeasNegative(scip, SCIPgetVarUbAtIndex(scip, var, bdchgidx, FALSE)) )
   {
      SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler
 *
 *  Let lb and ub be the lower and upper bounds of a
 *  variable. Preprocessing usually makes sure that lb <= 0 <= ub.
 *
 *  - If lb < 0 then rounding down may violate the constraint.
 *  - If ub > 0 then rounding up may violated the constraint.
 *  - If lb > 0 or ub < 0 then the constraint is infeasible and we do
 *    not have to deal with it here.
 *  - If lb == 0 then rounding down does not violate the constraint.
 *  - If ub == 0 then rounding up does not violate the constraint.
 */
static
SCIP_DECL_CONSLOCK(consLockSOS1)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPdebugMsg(scip, "Locking constraint <%s>.\n", SCIPconsGetName(cons));

   vars = consdata->vars;
   nvars = consdata->nvars;
   assert( vars != NULL );

   for (j = 0; j < nvars; ++j)
   {
      SCIP_VAR* var;
      var = vars[j];

      /* if lower bound is negative, rounding down may violate constraint */
      if ( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)) )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, var, nlockspos, nlocksneg) );
      }

      /* additionally: if upper bound is positive, rounding up may violate constraint */
      if ( SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, var, nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   for (j = 0; j < consdata->nvars; ++j)
   {
      if ( j > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[j], FALSE) );
      if ( consdata->weights == NULL )
         SCIPinfoMessage(scip, file, " (%d)", j+1);
      else
         SCIPinfoMessage(scip, file, " (%3.2f)", consdata->weights[j]);
   }

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** targetvars;
   SCIP_Real* sourceweights;
   SCIP_Real* targetweights;
   const char* consname;
   int nvars;
   int v;

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );

   *valid = TRUE;

   if ( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIPdebugMsg(scip, "Copying SOS1 constraint <%s> ...\n", consname);

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert( sourceconsdata != NULL );

   /* get variables and weights of the source constraint */
   nvars = sourceconsdata->nvars;

   if ( nvars == 0 )
      return SCIP_OKAY;

   sourcevars = sourceconsdata->vars;
   assert( sourcevars != NULL );
   sourceweights = sourceconsdata->weights;
   assert( sourceweights != NULL );

   /* duplicate variable array */
   SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetvars, nvars) );
   SCIP_CALL( SCIPduplicateBufferArray(sourcescip, &targetweights, sourceweights, nvars) );

   /* get copied variables in target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &(targetvars[v]), varmap, consmap, global, valid) );
   }

    /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsSOS1(scip, cons, consname, nvars, targetvars, targetweights,
            initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(sourcescip, &targetweights);
   SCIPfreeBufferArray(sourcescip, &targetvars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSOS1)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_Real weight;
   const char* s;
   char* t;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert(cons != NULL);
   assert(success != NULL);

   *success = TRUE;
   s = str;

   /* create empty SOS1 constraint */
   SCIP_CALL( SCIPcreateConsSOS1(scip, cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );

   /* loop through string */
   do
   {
      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, s, &var, &t) );
      s = t;

      /* skip until beginning of weight */
      while ( *s != '\0' && *s != '(' )
         ++s;

      if ( *s == '\0' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected weight at input: %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      /* skip '(' */
      ++s;

      /* find weight */
      weight = strtod(s, &t);
      if ( t == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error during parsing of the weight: %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      s = t;

      /* skip white space, ',', and ')' */
      while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' || *s == ')' ) )
         ++s;

      /* add variable */
      SCIP_CALL( SCIPaddVarSOS1(scip, *cons, var, weight) );
   }
   while ( *s != '\0' );

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/** exec the event handler
 *
 *  We update the number of variables fixed to be nonzero
 */
static
SCIP_DECL_EVENTEXEC(eventExecSOS1)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   SCIP_Real oldbound;
   SCIP_Real newbound;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );
   assert( event != NULL );

   cons = (SCIP_CONS*)eventdata;
   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( 0 <= consdata->nfixednonzeros && consdata->nfixednonzeros <= consdata->nvars );

   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   eventtype = SCIPeventGetType(event);
   switch ( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasPositive(scip, oldbound) && SCIPisFeasPositive(scip, newbound) )
      {
         conshdlr = SCIPconsGetHdlr(cons);
         assert( conshdlr != NULL );
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert( conshdlrdata != NULL );

         /* store variable fixed to be nonzero on stack */
         assert( 0 <= conshdlrdata->nfixnonzerovars && conshdlrdata->nfixnonzerovars <= SCIPgetNTotalVars(scip) );
         if ( conshdlrdata->nfixnonzerovars < conshdlrdata->maxnfixnonzerovars )
         {
            assert( conshdlrdata->fixnonzerovars != NULL );
            assert( SCIPeventGetVar(event) != NULL );
            conshdlrdata->fixnonzerovars[conshdlrdata->nfixnonzerovars++] = SCIPeventGetVar(event);
         }

         ++(consdata->nfixednonzeros);
      }
      break;

   case SCIP_EVENTTYPE_UBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasNegative(scip, oldbound) && SCIPisFeasNegative(scip, newbound) )
      {
         conshdlr = SCIPconsGetHdlr(cons);
         assert( conshdlr != NULL );
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert( conshdlrdata != NULL );

         /* store variable fixed to be nonzero on stack */
         assert( 0 <= conshdlrdata->nfixnonzerovars && conshdlrdata->nfixnonzerovars <= SCIPgetNTotalVars(scip) );
         if ( conshdlrdata->nfixnonzerovars < conshdlrdata->maxnfixnonzerovars )
         {
            assert( conshdlrdata->fixnonzerovars != NULL );
            assert( SCIPeventGetVar(event) != NULL );
            conshdlrdata->fixnonzerovars[conshdlrdata->nfixnonzerovars++] = SCIPeventGetVar(event);
         }

         ++(consdata->nfixednonzeros);
      }
      break;

   case SCIP_EVENTTYPE_LBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasPositive(scip, oldbound) && ! SCIPisFeasPositive(scip, newbound) )
         --(consdata->nfixednonzeros);
      break;

   case SCIP_EVENTTYPE_UBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasNegative(scip, oldbound) && ! SCIPisFeasNegative(scip, newbound) )
         --(consdata->nfixednonzeros);
      break;

   default:
      SCIPerrorMessage("invalid event type.\n");
      return SCIP_INVALIDDATA;
   }
   assert( 0 <= consdata->nfixednonzeros && consdata->nfixednonzeros <= consdata->nvars );

   SCIPdebugMsg(scip, "changed bound of variable <%s> from %f to %f (nfixednonzeros: %d).\n", SCIPvarGetName(SCIPeventGetVar(event)),
                    oldbound, newbound, consdata->nfixednonzeros);

   return SCIP_OKAY;
}


/** constraint handler method to determine a diving variable by assigning a variable and two values for diving */
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsSOS1)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( diveset != NULL );
   assert( success != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;
   *success = FALSE;

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("not an SOS1 constraint handler.\n");
      return SCIP_INVALIDDATA;
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* if the SOS1 constraints do not overlap, we apply a faster method getDiveBdChgsSOS1constraints() that does not make use of the conflict graph;
    * for overlapping SOS1 constraints we apply the method getDiveBdChgsSOS1conflictgraph(), which then may produce better results (e.g. due to more
    * diving candidates) */
   if ( conshdlrdata->switchsos1branch )
   {
      SCIP_CALL( getDiveBdChgsSOS1constraints(scip, conshdlr, diveset, sol, success) );
   }
   else
   {
      SCIP_CALL( getDiveBdChgsSOS1conflictgraph(scip, conshdlr, diveset, sol, success) );
   }

   return SCIP_OKAY;
}


/* ---------------- Constraint specific interface methods ---------------- */

/** creates the handler for SOS1 constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSOS1(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->branchsos = TRUE;
   conshdlrdata->switchsos1branch = FALSE;
   conshdlrdata->switchcutsfromsos1 = FALSE;
   conshdlrdata->eventhdlr = NULL;
   conshdlrdata->fixnonzerovars = NULL;
   conshdlrdata->maxnfixnonzerovars = 0;
   conshdlrdata->nfixnonzerovars = 0;
   conshdlrdata->conflictgraph = NULL;
   conshdlrdata->localconflicts = NULL;
   conshdlrdata->isconflocal = FALSE;
   conshdlrdata->implgraph = NULL;
   conshdlrdata->nimplnodes = 0;
   conshdlrdata->nboundcuts = 0;
   conshdlrdata->tcliquegraph = NULL;
   conshdlrdata->tcliquedata = NULL;
   conshdlrdata->cntextsos1 = -1;
   conshdlrdata->varhash = NULL;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSOS1, NULL) );
   if ( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for SOS1 constraints not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSOS1, consEnfopsSOS1, consCheckSOS1, consLockSOS1, conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySOS1, consCopySOS1) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSOS1) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsSOS1) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolSOS1) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolSOS1) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSOS1) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSOS1) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSOS1) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSOS1) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSOS1) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSOS1, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSOS1) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSOS1, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSOS1) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSOS1, consSepasolSOS1, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSOS1) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSOS1) );

   /* add SOS1 constraint handler parameters */

   /* adjacency matrix parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxsosadjacency",
         "do not create an adjacency matrix if number of SOS1 variables is larger than predefined value (-1: no limit)",
         &conshdlrdata->maxsosadjacency, TRUE, DEFAULT_MAXSOSADJACENCY, -1, INT_MAX, NULL, NULL) );

   /* presolving parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxextensions",
         "maximal number of extensions that will be computed for each SOS1 constraint  (-1: no limit)",
         &conshdlrdata->maxextensions, TRUE, DEFAULT_MAXEXTENSIONS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxtightenbds",
         "maximal number of bound tightening rounds per presolving round (-1: no limit)",
         &conshdlrdata->maxtightenbds, TRUE, DEFAULT_MAXTIGHTENBDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/perfimplanalysis",
         "if TRUE then perform implication graph analysis (might add additional SOS1 constraints)",
         &conshdlrdata->perfimplanalysis, TRUE, DEFAULT_PERFIMPLANALYSIS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/depthimplanalysis",
         "number of recursive calls of implication graph analysis (-1: no limit)",
         &conshdlrdata->depthimplanalysis, TRUE, DEFAULT_DEPTHIMPLANALYSIS, -1, INT_MAX, NULL, NULL) );

   /* propagation parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/conflictprop",
         "whether to use conflict graph propagation",
         &conshdlrdata->conflictprop, TRUE, DEFAULT_CONFLICTPROP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/implprop",
         "whether to use implication graph propagation",
         &conshdlrdata->implprop, TRUE, DEFAULT_IMPLPROP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/sosconsprop",
         "whether to use SOS1 constraint propagation",
         &conshdlrdata->sosconsprop, TRUE, DEFAULT_SOSCONSPROP, NULL, NULL) );

   /* branching parameters */
   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/branchingrule",
         "which branching rule should be applied ? ('n': neighborhood, 'b': bipartite, 's': SOS1/clique) (note: in some cases an automatic switching to SOS1 branching is possible)",
         &conshdlrdata->branchingrule, TRUE, DEFAULT_BRANCHINGRULE, DEFAULT_BRANCHSTRATEGIES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/autosos1branch",
         "if TRUE then automatically switch to SOS1 branching if the SOS1 constraints do not overlap",
         &conshdlrdata->autosos1branch, TRUE, DEFAULT_AUTOSOS1BRANCH, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/fixnonzero",
         "if neighborhood branching is used, then fix the branching variable (if positive in sign) to the value of the feasibility tolerance",
         &conshdlrdata->fixnonzero, TRUE, DEFAULT_FIXNONZERO, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/addcomps",
         "if TRUE then add complementarity constraints to the branching nodes (can be used in combination with neighborhood or bipartite branching)",
         &conshdlrdata->addcomps, TRUE, DEFAULT_ADDCOMPS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxaddcomps",
         "maximal number of complementarity constraints added per branching node (-1: no limit)",
         &conshdlrdata->maxaddcomps, TRUE, DEFAULT_MAXADDCOMPS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/addcompsfeas",
         "minimal feasibility value for complementarity constraints in order to be added to the branching node",
         &conshdlrdata->addcompsfeas, TRUE, DEFAULT_ADDCOMPSFEAS, -SCIP_REAL_MAX, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/addbdsfeas",
         "minimal feasibility value for bound inequalities in order to be added to the branching node",
         &conshdlrdata->addbdsfeas, TRUE, DEFAULT_ADDBDSFEAS, -SCIP_REAL_MAX, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/addextendedbds",
         "should added complementarity constraints be extended to SOS1 constraints to get tighter bound inequalities",
         &conshdlrdata->addextendedbds, TRUE, DEFAULT_ADDEXTENDEDBDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/branchsos",
         "Use SOS1 branching in enforcing (otherwise leave decision to branching rules)? This value can only be set to false if all SOS1 variables are binary",
         &conshdlrdata->branchsos, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/branchnonzeros",
         "Branch on SOS constraint with most number of nonzeros?",
         &conshdlrdata->branchnonzeros, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/branchweight",
         "Branch on SOS cons. with highest nonzero-variable weight for branching (needs branchnonzeros = false)?",
         &conshdlrdata->branchweight, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/addcompsdepth",
         "only add complementarity constraints to branching nodes for predefined depth (-1: no limit)",
         &conshdlrdata->addcompsdepth, TRUE, DEFAULT_ADDCOMPSDEPTH, -1, INT_MAX, NULL, NULL) );

   /* selection rule parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/nstrongrounds",
         "maximal number of strong branching rounds to perform for each node (-1: auto); only available for neighborhood and bipartite branching",
         &conshdlrdata->nstrongrounds, TRUE, DEFAULT_NSTRONGROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/nstrongiter",
         "maximal number LP iterations to perform for each strong branching round (-2: auto, -1: no limit)",
         &conshdlrdata->nstrongiter, TRUE, DEFAULT_NSTRONGITER, -2, INT_MAX, NULL, NULL) );

   /* separation parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/boundcutsfromsos1",
         "if TRUE separate bound inequalities from initial SOS1 constraints",
         &conshdlrdata->boundcutsfromsos1, TRUE, DEFAULT_BOUNDCUTSFROMSOS1, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/boundcutsfromgraph",
         "if TRUE separate bound inequalities from the conflict graph",
         &conshdlrdata->boundcutsfromgraph, TRUE, DEFAULT_BOUNDCUTSFROMGRAPH, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/autocutsfromsos1",
         "if TRUE then automatically switch to separating initial SOS1 constraints if the SOS1 constraints do not overlap",
         &conshdlrdata->autocutsfromsos1, TRUE, DEFAULT_AUTOCUTSFROMSOS1, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/boundcutsfreq",
         "frequency for separating bound cuts; zero means to separate only in the root node",
         &conshdlrdata->boundcutsfreq, TRUE, DEFAULT_BOUNDCUTSFREQ, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/boundcutsdepth",
         "node depth of separating bound cuts (-1: no limit)",
         &conshdlrdata->boundcutsdepth, TRUE, DEFAULT_BOUNDCUTSDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxboundcuts",
         "maximal number of bound cuts separated per branching node",
         &conshdlrdata->maxboundcuts, TRUE, DEFAULT_MAXBOUNDCUTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxboundcutsroot",
         "maximal number of bound cuts separated per iteration in the root node",
         &conshdlrdata->maxboundcutsroot, TRUE, DEFAULT_MAXBOUNDCUTSROOT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/strthenboundcuts",
         "if TRUE then bound cuts are strengthened in case bound variables are available",
         &conshdlrdata->strthenboundcuts, TRUE, DEFAULT_STRTHENBOUNDCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/implcutsfreq",
         "frequency for separating implied bound cuts; zero means to separate only in the root node",
         &conshdlrdata->implcutsfreq, TRUE, DEFAULT_IMPLCUTSFREQ, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/implcutsdepth",
         "node depth of separating implied bound cuts (-1: no limit)",
         &conshdlrdata->implcutsdepth, TRUE, DEFAULT_IMPLCUTSDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maximplcuts",
         "maximal number of implied bound cuts separated per branching node",
         &conshdlrdata->maximplcuts, TRUE, DEFAULT_MAXIMPLCUTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maximplcutsroot",
         "maximal number of implied bound cuts separated per iteration in the root node",
         &conshdlrdata->maximplcutsroot, TRUE, DEFAULT_MAXIMPLCUTSROOT, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}


/** creates and captures a SOS1 constraint
 *
 *  We set the constraint to not be modifable. If the weights are non NULL, the variables are ordered according to these
 *  weights (in ascending order).
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights,            /**< weights determining the variable order, or NULL if natural order should be used */
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
   SCIP_Bool modifiable;
   SCIP_Bool transformed;
   int v;

   modifiable = FALSE;

   /* find the SOS1 constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("<%s> constraint handler not found\n", CONSHDLR_NAME);
      return SCIP_PLUGINNOTFOUND;
   }

   /* are we in the transformed problem? */
   transformed = SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED;

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->vars = NULL;
   consdata->nvars = nvars;
   consdata->maxvars = nvars;
   consdata->rowub = NULL;
   consdata->rowlb = NULL;
   consdata->nfixednonzeros = transformed ? 0 : -1;
   consdata->weights = NULL;
   consdata->local = local;

   if ( nvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->vars, vars, nvars) );

      /* check weights */
      if ( weights != NULL )
      {
         /* store weights */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->weights, weights, nvars) );

         /* sort variables - ascending order */
         SCIPsortRealPtr(consdata->weights, (void**)consdata->vars, nvars);
      }
   }
   else
   {
      assert( weights == NULL );
   }

   /* branching on multiaggregated variables does not seem to work well, so avoid it */
   for (v = 0; v < nvars; ++v)
   {
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, consdata->vars[v]) );
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );
   assert( transformed == SCIPconsIsTransformed(*cons) );

   /* replace original variables by transformed variables in transformed constraint, add locks, and catch events */
   for (v = nvars - 1; v >= 0; --v)
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* always use transformed variables in transformed constraints */
      if ( transformed )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, consdata->vars[v], &(consdata->vars[v])) );
      }
      assert( consdata->vars[v] != NULL );
      assert( transformed == SCIPvarIsTransformed(consdata->vars[v]) );

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      /* handle the new variable */
      SCIP_CALL( handleNewVariableSOS1(scip, *cons, consdata, conshdlrdata, consdata->vars[v], transformed) );
   }

   return SCIP_OKAY;
}


/** creates and captures a SOS1 constraint with all constraint flags set to their default values.
 *
 *  @warning Do NOT set the constraint to be modifiable manually, because this might lead
 *  to wrong results as the variable array will not be resorted
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights             /**< weights determining the variable order, or NULL if natural order should be used */
   )
{
   SCIP_CALL( SCIPcreateConsSOS1( scip, cons, name, nvars, vars, weights, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


/** adds variable to SOS1 constraint, the position is determined by the given weight */
SCIP_RETCODE SCIPaddVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_Real             weight              /**< weight determining position of variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   assert( scip != NULL );
   assert( var != NULL );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "adding variable <%s> to constraint <%s> with weight %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), weight);

   conshdlr = SCIPconsGetHdlr(cons);
   assert( conshdlr != NULL );
   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      return SCIP_INVALIDDATA;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIP_CALL( addVarSOS1(scip, cons, conshdlrdata, var, weight) );

   return SCIP_OKAY;
}


/** appends variable to SOS1 constraint */
SCIP_RETCODE SCIPappendVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   assert( scip != NULL );
   assert( var != NULL );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "appending variable <%s> to constraint <%s>\n", SCIPvarGetName(var), SCIPconsGetName(cons));

   conshdlr = SCIPconsGetHdlr(cons);
   assert( conshdlr != NULL );
   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      return SCIP_INVALIDDATA;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIP_CALL( appendVarSOS1(scip, cons, conshdlrdata, var) );

   return SCIP_OKAY;
}


/** gets number of variables in SOS1 constraint */
int SCIPgetNVarsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->nvars;
}


/** gets array of variables in SOS1 constraint */
SCIP_VAR** SCIPgetVarsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->vars;
}


/** gets array of weights in SOS1 constraint (or NULL if not existent) */
SCIP_Real* SCIPgetWeightsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->weights;
}


/** gets conflict graph of SOS1 constraints (or NULL if not existent)
 *
 *  @note The conflict graph is globally valid; local changes are not taken into account.
 */
SCIP_DIGRAPH* SCIPgetConflictgraphSOS1(
   SCIP_CONSHDLR*        conshdlr            /**< SOS1 constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("not an SOS1 constraint handler.\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   return conshdlrdata->conflictgraph;
}


/** gets number of problem variables that are part of the SOS1 conflict graph */
int SCIPgetNSOS1Vars(
   SCIP_CONSHDLR*        conshdlr            /**< SOS1 constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("not an SOS1 constraint handler.\n");
      SCIPABORT();
      return -1; /*lint !e527*/
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   return conshdlrdata->nsos1vars;
}


/** returns whether variable is part of the SOS1 conflict graph */
SCIP_Bool SCIPvarIsSOS1(
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( var != NULL );
   assert( conshdlr != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("not an SOS1 constraint handler.\n");
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   return varIsSOS1(conshdlrdata, var);
}


/** returns SOS1 index of variable or -1 if variable is not part of the SOS1 conflict graph */
int SCIPvarGetNodeSOS1(
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );
   assert( var != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("Not an SOS1 constraint handler.\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->varhash == NULL )
   {
      SCIPerrorMessage("Hashmap not yet initialized.\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   return varGetNodeSOS1(conshdlrdata, var);
}


/** returns variable that belongs to a given node from the conflict graph */
SCIP_VAR* SCIPnodeGetVarSOS1(
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   node                /**< node from the conflict graph */
   )
{
   SCIP_NODEDATA* nodedata;

   assert( conflictgraph != NULL );
   assert( node >= 0 && node < SCIPdigraphGetNNodes(conflictgraph) );

   /* get node data */
   nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, node);

   if ( nodedata == NULL )
   {
      SCIPerrorMessage("variable is not assigned to an index.\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   return nodedata->var;
}


/** based on solution values of the variables, fixes variables to zero to turn all SOS1 constraints feasible */
SCIP_RETCODE SCIPmakeSOS1sFeasible(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            changed,            /**< pointer to store whether the solution has been changed */
   SCIP_Bool*            success             /**< pointer to store whether SOS1 constraints have been turned feasible and
                                              *   solution was good enough */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real roundobjval;
   SCIP_Bool allroundable;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( sol != NULL );
   assert( changed != NULL );
   assert( success != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("Not an SOS1 constraint handler.\n");
      return SCIP_PARAMETERWRONGVAL;
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   *changed = FALSE;
   *success = FALSE;
   allroundable = FALSE;

   /* check number of SOS1 constraints */
   if ( SCIPconshdlrGetNConss(conshdlr) < 1 )
   {
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* if the SOS1 constraints do not overlap, we apply a faster method makeSOS1constraintsFeasible() that does not make use of the conflict graph;
    * for overlapping SOS1 constraints we apply the method makeSOS1conflictgraphFeasible(), which then may produce better feasible solutions */
   if ( conshdlrdata->switchsos1branch )
   {
      SCIP_CALL( makeSOS1constraintsFeasible(scip, conshdlr, sol, changed, &allroundable) );
   }
   else
   {
      SCIP_CALL( makeSOS1conflictgraphFeasible(scip, conshdlr, sol, changed, &allroundable) );
   }

   if ( ! allroundable )
      return SCIP_OKAY;

   /* check whether objective value of rounded solution is good enough */
   roundobjval = SCIPgetSolOrigObj(scip, sol);
   if ( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE )
      roundobjval *= -1;

   if ( SCIPisLT(scip, roundobjval, SCIPgetUpperbound(scip) ) )
      *success = TRUE;

   return SCIP_OKAY;
}
