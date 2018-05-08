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

/**@file   struct_var.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_VAR_H__
#define __SCIP_STRUCT_VAR_H__


#include "scip/def.h"
#include "scip/type_history.h"
#include "scip/type_event.h"
#include "scip/type_var.h"
#include "scip/type_implics.h"
#include "scip/type_cons.h"
#include "scip/type_prop.h"
#include "scip/type_lp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** hole in a domain */
struct SCIP_Hole
{
   SCIP_Real             left;               /**< left bound of open interval defining the hole (left,right) */
   SCIP_Real             right;              /**< right bound of open interval defining the hole (left,right) */
};

/** list of domain holes */
struct SCIP_Holelist
{
   SCIP_HOLE             hole;               /**< this hole */
   SCIP_HOLELIST*        next;               /**< next hole in list */
};

/** change in a hole list */
struct SCIP_HoleChg
{
   SCIP_HOLELIST**       ptr;                /**< changed list pointer */
   SCIP_HOLELIST*        newlist;            /**< new value of list pointer */
   SCIP_HOLELIST*        oldlist;            /**< old value of list pointer */
};

/** data for branching decision bound changes */
struct SCIP_BranchingData
{
   SCIP_Real             lpsolval;           /**< sol val of var in last LP prior to bound change, or SCIP_INVALID if unknown */
};

/** data for infered bound changes */
struct SCIP_InferenceData
{
   SCIP_VAR*             var;                /**< variable that was changed (parent of var, or var itself) */
   union
   {
      SCIP_CONS*         cons;               /**< constraint that infered this bound change, or NULL */
      SCIP_PROP*         prop;               /**< propagator that infered this bound change, or NULL */
   } reason;
   int                   info;               /**< user information for inference to help resolving the conflict */
};

/** change in one bound of a variable */
struct SCIP_BoundChg
{
   SCIP_Real             newbound;           /**< new value for bound */
   union
   {
      SCIP_BRANCHINGDATA branchingdata;      /**< data for branching decisions */
      SCIP_INFERENCEDATA inferencedata;      /**< data for infered bound changes */
   } data;
   SCIP_VAR*             var;                /**< active variable to change the bounds for */
   unsigned int          boundchgtype:2;     /**< bound change type: branching decision or infered bound change */
   unsigned int          boundtype:1;        /**< type of bound for var: lower or upper bound */
   unsigned int          inferboundtype:1;   /**< type of bound for inference var (see inference data): lower or upper bound */
   unsigned int          applied:1;          /**< was this bound change applied at least once? */
   unsigned int          redundant:1;        /**< is this bound change redundant? */
};

/** bound change index representing the time of the bound change in path from root to current node */
struct SCIP_BdChgIdx
{
   int                   depth;              /**< depth of node where the bound change was created */
   int                   pos;                /**< position of bound change in node's domchg array */
};

/** bound change information to track bound changes from root node to current node */
struct SCIP_BdChgInfo
{
   SCIP_Real             oldbound;           /**< old value for bound */
   SCIP_Real             newbound;           /**< new value for bound */
   SCIP_VAR*             var;                /**< active variable that changed the bounds */
   SCIP_INFERENCEDATA    inferencedata;      /**< data for infered bound changes */
   SCIP_BDCHGIDX         bdchgidx;           /**< bound change index in path from root to current node */
   unsigned int          pos:27;             /**< position in the variable domain change array */
   unsigned int          boundchgtype:2;     /**< bound change type: branching decision or infered bound change */
   unsigned int          boundtype:1;        /**< type of bound for var: lower or upper bound */
   unsigned int          inferboundtype:1;   /**< type of bound for inference var (see inference data): lower or upper bound */
   unsigned int          redundant:1;        /**< does the bound change info belong to a redundant bound change? */
};

/** tracks changes of the variables' domains (static arrays, bound changes only) */
struct SCIP_DomChgBound
{
   unsigned int          nboundchgs:30;      /**< number of bound changes (must be first structure entry!) */
   unsigned int          domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   SCIP_BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
};

/** tracks changes of the variables' domains (static arrays, bound and hole changes) */
struct SCIP_DomChgBoth
{
   unsigned int          nboundchgs:30;      /**< number of bound changes (must be first structure entry!) */
   unsigned int          domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   SCIP_BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
   SCIP_HOLECHG*         holechgs;           /**< array with changes in hole lists */
   int                   nholechgs;          /**< number of hole list changes */
};

/** tracks changes of the variables' domains (dynamic arrays) */
struct SCIP_DomChgDyn
{
   unsigned int          nboundchgs:30;      /**< number of bound changes (must be first structure entry!) */
   unsigned int          domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   SCIP_BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
   SCIP_HOLECHG*         holechgs;           /**< array with changes in hole lists */
   int                   nholechgs;          /**< number of hole list changes */
   int                   boundchgssize;      /**< size of bound changes array */
   int                   holechgssize;       /**< size of hole changes array */
};

/** tracks changes of the variables' domains */
union SCIP_DomChg
{
   SCIP_DOMCHGBOUND      domchgbound;        /**< bound changes */
   SCIP_DOMCHGBOTH       domchgboth;         /**< bound and hole changes */
   SCIP_DOMCHGDYN        domchgdyn;          /**< bound and hole changes with dynamic arrays */
};

/** domain of a variable */
struct SCIP_Dom
{
   SCIP_Real             lb;                 /**< lower bounds of variables */
   SCIP_Real             ub;                 /**< upper bounds of variables */
   SCIP_HOLELIST*        holelist;           /**< list of holes */
};

/** original variable information */
struct SCIP_Original
{
   SCIP_DOM              origdom;            /**< domain of variable in original problem */
   SCIP_VAR*             transvar;           /**< pointer to representing transformed variable */
};

/** aggregation information: x = a*y + c */
struct SCIP_Aggregate
{
   SCIP_Real             scalar;             /**< multiplier a in aggregation */
   SCIP_Real             constant;           /**< constant shift c in aggregation */
   SCIP_VAR*             var;                /**< variable y in aggregation */
};

/** multiple aggregation information: x = a_1*y_1 + ... + a_k*y_k + c */
struct SCIP_Multaggr
{
   SCIP_Real             constant;           /**< constant shift c in multiple aggregation */
   SCIP_Real*            scalars;            /**< multipliers a in multiple aggregation */
   SCIP_VAR**            vars;               /**< variables y in multiple aggregation */
   int                   nvars;              /**< number of variables in aggregation */
   int                   varssize;           /**< size of vars and scalars arrays */
};

/** negation information: x' = c - x */
struct SCIP_Negate
{
   SCIP_Real             constant;           /**< constant shift c in negation */
};

/** variable of the problem */
struct SCIP_Var
{
#ifndef NDEBUG
   SCIP*                 scip;               /**< SCIP data structure */
#endif
   SCIP_Real             obj;                /**< objective function value of variable (might be changed temporarily in probing mode)*/
   SCIP_Real             unchangedobj;       /**< unchanged objective function value of variable (ignoring temporary changes in probing mode) */
   SCIP_Real             branchfactor;       /**< factor to weigh variable's branching score with */
   SCIP_Real             rootsol;            /**< last primal solution of variable in root node, or zero */
   SCIP_Real             bestrootsol;        /**< best primal solution of variable in root node, or zero, w.r.t. root LP value and root reduced cost */
   SCIP_Real             bestrootredcost;    /**< best reduced costs of variable in root node, or zero, w.r.t. root LP value and root solution value */
   SCIP_Real             bestrootlpobjval;   /**< best root LP objective value, or SCIP_INVALID, w.r.t. root solution value and root reduced cost */
   SCIP_Real             relaxsol;           /**< primal solution of variable in current relaxation solution, or SCIP_INVALID */
   SCIP_Real             nlpsol;             /**< primal solution of variable in current NLP solution, or SCIP_INVALID */
   SCIP_Real             primsolavg;         /**< weighted average of all values of variable in primal feasible solutions */
   SCIP_Real             conflictlb;         /**< maximal lower bound of variable in the current conflict */
   SCIP_Real             conflictub;         /**< minimal upper bound of variable in the current conflict */
   SCIP_Real             conflictrelaxedlb;  /**< minimal relaxed lower bound of variable in the current conflict (conflictrelqxlb <= conflictlb) */
   SCIP_Real             conflictrelaxedub;  /**< minimal release upper bound of variable in the current conflict (conflictrelqxlb <= conflictlb) */
   SCIP_Real             lazylb;             /**< global lower bound that is ensured by constraints and has not to be added to the LP */
   SCIP_Real             lazyub;             /**< global upper bound that is ensured by constraints and has not to be added to the LP */
   SCIP_DOM              glbdom;             /**< domain of variable in global problem */
   SCIP_DOM              locdom;             /**< domain of variable in current subproblem */
   union
   {
      SCIP_ORIGINAL      original;           /**< original variable information */
      SCIP_COL*          col;                /**< LP column (for column variables) */
      SCIP_AGGREGATE     aggregate;          /**< aggregation information (for aggregated variables) */
      SCIP_MULTAGGR      multaggr;           /**< multiple aggregation information (for multiple aggregated variables) */
      SCIP_NEGATE        negate;             /**< negation information (for negated variables) */
   } data;
   char*                 name;               /**< name of the variable */
   SCIP_DECL_VARCOPY     ((*varcopy));       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig));    /**< frees user data of original variable */
   SCIP_DECL_VARTRANS    ((*vartrans));      /**< creates transformed user data by transforming original user data */
   SCIP_DECL_VARDELTRANS ((*vardeltrans));   /**< frees user data of transformed variable */
   SCIP_VARDATA*         vardata;            /**< user data for this specific variable */
   SCIP_VAR**            parentvars;         /**< parent variables in the aggregation tree */
   SCIP_VAR*             negatedvar;         /**< pointer to the variables negation: x' = lb + ub - x, or NULL if not created */
   SCIP_VBOUNDS*         vlbs;               /**< variable lower bounds x >= b*y + d */
   SCIP_VBOUNDS*         vubs;               /**< variable upper bounds x <= b*y + d */
   SCIP_IMPLICS*         implics;            /**< implications y >=/<= b following from x <= 0 and x >= 1 (x binary), or NULL if x is not binary */
   SCIP_CLIQUELIST*      cliquelist;         /**< list of cliques the variable and its negation is member of */
   SCIP_EVENTFILTER*     eventfilter;        /**< event filter for events concerning this variable; not for ORIGINAL vars */
   SCIP_BDCHGINFO*       lbchginfos;         /**< bound change informations for lower bound changes from root to current node */
   SCIP_BDCHGINFO*       ubchginfos;         /**< bound change informations for upper bound changes from root to current node */
   SCIP_HISTORY*         history;            /**< branching and inference history information */
   SCIP_HISTORY*         historycrun;        /**< branching and inference history information for current run */
   SCIP_VALUEHISTORY*    valuehistory;       /**< branching and inference history information which are value based, or NULL if not used */
   SCIP_Longint          closestvblpcount;   /**< LP count for which the closestvlbidx/closestvubidx entries are valid */
   int                   index;              /**< consecutively numbered variable identifier */
   int                   probindex;          /**< array position in problems vars array, or -1 if not assigned to a problem */
   int                   pseudocandindex;    /**< array position in pseudo branching candidates array, or -1 */
   int                   eventqueueindexobj; /**< array position in event queue of objective change event, or -1 */
   int                   eventqueueindexlb;  /**< array position in event queue of lower bound change event, or -1 */
   int                   eventqueueindexub;  /**< array position in event queue of upper bound change event, or -1 */
   int                   parentvarssize;     /**< available slots in parentvars array */
   int                   nparentvars;        /**< number of parent variables in aggregation tree (used slots of parentvars) */
   int                   nuses;              /**< number of times, this variable is referenced */
   int                   nlocksdown;         /**< number of locks for rounding down; if zero, rounding down is always feasible */
   int                   nlocksup;           /**< number of locks for rounding up; if zero, rounding up is always feasible */
   int                   branchpriority;     /**< priority of the variable for branching */
   int                   lbchginfossize;     /**< available slots in lbchginfos array */
   int                   nlbchginfos;        /**< number of lower bound changes from root node to current node */
   int                   ubchginfossize;     /**< available slots in ubchginfos array */
   int                   nubchginfos;        /**< number of upper bound changes from root node to current node */
   int                   conflictlbcount;    /**< number of last conflict, the lower bound was member of */
   int                   conflictubcount;    /**< number of last conflict, the upper bound was member of */
   int                   closestvlbidx;      /**< index of closest VLB variable in current LP solution, or -1 */
   int                   closestvubidx;      /**< index of closest VUB variable in current LP solution, or -1 */
   unsigned int          initial:1;          /**< TRUE iff var's column should be present in the initial root LP */
   unsigned int          removable:1;        /**< TRUE iff var's column is removable from the LP (due to aging or cleanup) */
   unsigned int          deletable:1;        /**< TRUE iff the variable is removable from the problem */
   unsigned int          deleted:1;          /**< TRUE iff variable was marked for deletion from the problem */
   unsigned int          donotmultaggr:1;    /**< TRUE iff variable is not allowed to be multi-aggregated */
   unsigned int          vartype:2;          /**< type of variable: binary, integer, implicit integer, continuous */
   unsigned int          varstatus:3;        /**< status of variable: original, loose, column, fixed, aggregated, multiaggregated, negated */
   unsigned int          pseudocostflag:2;   /**< temporary flag used in pseudo cost update */
   unsigned int          branchdirection:2;  /**< preferred branching direction of the variable (downwards, upwards, auto) */
   unsigned int          eventqueueimpl:1;   /**< is an IMPLADDED event on this variable currently in the event queue? */
   unsigned int          delglobalstructs:1; /**< is variable marked to be removed from global structures (cliques etc.)? */
};

#ifdef __cplusplus
}
#endif

#endif
