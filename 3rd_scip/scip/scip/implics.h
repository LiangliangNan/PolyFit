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

/**@file   implics.h
 * @ingroup INTERNALAPI
 * @brief  methods for implications, variable bounds, and cliques
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_IMPLICS_H__
#define __SCIP_IMPLICS_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_implics.h"
#include "scip/type_branch.h"
#include "scip/pub_implics.h"

#ifdef NDEBUG
#include "scip/struct_implics.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Methods for Variable Bounds
 */

/** frees a variable bounds data structure */
extern
void SCIPvboundsFree(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** adds a variable bound to the variable bounds data structure */
extern
SCIP_RETCODE SCIPvboundsAdd(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        vboundtype,         /**< type of variable bound (LOWER or UPPER) */
   SCIP_VAR*             var,                /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             coef,               /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             constant,           /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Bool*            added               /**< pointer to store whether the variable bound was added */
   );

/** removes from variable x a variable bound x >=/<= b*z + d with binary or integer z */
extern
SCIP_RETCODE SCIPvboundsDel(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             vbdvar,             /**< variable z    in x >=/<= b*z + d */
   SCIP_Bool             negativecoef        /**< is coefficient b negative? */
   );

/** reduces the number of variable bounds stored in the given variable bounds data structure */
void SCIPvboundsShrink(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   newnvbds            /**< new number of variable bounds */
   );


/** gets number of variable bounds contained in given variable bounds data structure */
extern
int SCIPvboundsGetNVbds(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   );

/** gets array of variables contained in given variable bounds data structure */
extern
SCIP_VAR** SCIPvboundsGetVars(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   );

/** gets array of coefficients contained in given variable bounds data structure */
extern
SCIP_Real* SCIPvboundsGetCoefs(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   );

/** gets array of constants contained in given variable bounds data structure */
extern
SCIP_Real* SCIPvboundsGetConstants(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPvboundsGetNVbds(vbounds)     ((vbounds) != NULL ? (vbounds)->len : 0)
#define SCIPvboundsGetVars(vbounds)      ((vbounds) != NULL ? (vbounds)->vars : NULL)
#define SCIPvboundsGetCoefs(vbounds)     ((vbounds) != NULL ? (vbounds)->coefs : NULL)
#define SCIPvboundsGetConstants(vbounds) ((vbounds) != NULL ? (vbounds)->constants : NULL)

#endif




/*
 * Methods for Implications
 */

/** frees an implications data structure */
extern
void SCIPimplicsFree(
   SCIP_IMPLICS**        implics,            /**< pointer of implications data structure to free */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** adds an implication x == 0/1 -> y <= b or y >= b to the implications data structure;
 *  the implication must be non-redundant
 */
extern
SCIP_RETCODE SCIPimplicsAdd(
   SCIP_IMPLICS**        implics,            /**< pointer to implications data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             varfixing,          /**< FALSE if implication for x == 0 has to be added, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool             isshortcut,         /**< is the implication a shortcut, i.e., added as part of the transitive closure of another implication? */
   SCIP_Bool*            conflict,           /**< pointer to store whether implication causes a conflict for variable x */
   SCIP_Bool*            added               /**< pointer to store whether the implication was added */
   );

/** removes the implication  x <= 0 or x >= 1  ==>  y <= b  or  y >= b  from the implications data structure */
extern
SCIP_RETCODE SCIPimplicsDel(
   SCIP_IMPLICS**        implics,            /**< pointer to implications data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             varfixing,          /**< FALSE if y should be removed from implications for x <= 0, TRUE for x >= 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype            /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   );

/** returns which implications on given variable y are contained in implications for x == 0 or x == 1 */
extern
void SCIPimplicsGetVarImplics(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_Bool*            haslowerimplic,     /**< pointer to store whether there exists an implication y >= l */
   SCIP_Bool*            hasupperimplic      /**< pointer to store whether there exists an implication y <= u */
   );

/** returns whether an implication y <= b or y >= b is contained in implications for x == 0 or x == 1 */
extern
SCIP_Bool SCIPimplicsContainsImpl(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_BOUNDTYPE        impltype            /**< type of implication y <=/>= b to search for */
   );


/** gets number of implications for a given binary variable fixing */
extern
int SCIPimplicsGetNImpls(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** gets array with implied variables for a given binary variable fixing */
extern
SCIP_VAR** SCIPimplicsGetVars(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** gets array with implication types for a given binary variable fixing */
extern
SCIP_BOUNDTYPE* SCIPimplicsGetTypes(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** gets array with implication bounds for a given binary variable fixing */
extern
SCIP_Real* SCIPimplicsGetBounds(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** Gets array with unique implication identifiers for a given binary variable fixing.
 *  If an implication is a shortcut, i.e., it was added as part of the transitive closure of another implication,
 *  its id is negative, otherwise it is nonnegative.
 */
extern
int* SCIPimplicsGetIds(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPimplicsGetNImpls(implics, varfixing)       ((implics) != NULL ? (implics)->nimpls[varfixing] : 0)
#define SCIPimplicsGetNBinImpls(implics, varfixing)    ((implics) != NULL ? (implics)->nbinimpls[varfixing] : 0)
#define SCIPimplicsGetVars(implics, varfixing)         ((implics) != NULL ? (implics)->vars[varfixing] : NULL)
#define SCIPimplicsGetTypes(implics, varfixing)        ((implics) != NULL ? (implics)->types[varfixing] : NULL)
#define SCIPimplicsGetBounds(implics, varfixing)       ((implics) != NULL ? (implics)->bounds[varfixing] : NULL)
#define SCIPimplicsGetIds(implics, varfixing)          ((implics) != NULL ? (implics)->ids[varfixing] : NULL)

#endif




/*
 * methods for cliques
 */

/** adds a single variable to the given clique */
extern
SCIP_RETCODE SCIPcliqueAddVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to add to the clique */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_Bool*            doubleentry,        /**< pointer to store whether the variable and value occurs twice in the clique */
   SCIP_Bool*            oppositeentry       /**< pointer to store whether the variable with opposite value is in the clique */
   );

/** removes a single variable from the given clique */
extern
void SCIPcliqueDelVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< variable to remove from the clique */
   SCIP_Bool             value               /**< value of the variable in the clique */
   );

/** frees a clique list data structure */
extern
void SCIPcliquelistFree(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** adds a clique to the clique list */
extern
SCIP_RETCODE SCIPcliquelistAdd(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             value,              /**< value of the variable for which the clique list should be extended */
   SCIP_CLIQUE*          clique              /**< clique that should be added to the clique list */
   );

/** removes a clique from the clique list */
extern
SCIP_RETCODE SCIPcliquelistDel(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             value,              /**< value of the variable for which the clique list should be reduced */
   SCIP_CLIQUE*          clique              /**< clique that should be deleted from the clique list */
   );

/** returns whether the given clique lists have a non-empty intersection, i.e. whether there is a clique that appears
 *  in both lists
 */
extern
SCIP_Bool SCIPcliquelistsHaveCommonClique(
   SCIP_CLIQUELIST*      cliquelist1,        /**< first clique list data structure */
   SCIP_Bool             value1,             /**< value of first variable */
   SCIP_CLIQUELIST*      cliquelist2,        /**< second clique list data structure */
   SCIP_Bool             value2              /**< value of second variable */
   );

/** removes all listed entries from the cliques */
extern
void SCIPcliquelistRemoveFromCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< active problem variable the clique list belongs to */
   SCIP_Bool             irrelevantvar       /**< has the variable become irrelevant, meaning that equality
                                              *   cliques need to be relaxed? */
   );

/** creates a clique table data structure */
extern
SCIP_RETCODE SCIPcliquetableCreate(
   SCIP_CLIQUETABLE**    cliquetable,        /**< pointer to store clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** frees a clique table data structure */
extern
SCIP_RETCODE SCIPcliquetableFree(
   SCIP_CLIQUETABLE**    cliquetable,        /**< pointer to store clique table data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** adds a clique to the clique table, using the given values for the given variables;
 *  performs implications if the clique contains the same variable twice
 */
extern
SCIP_RETCODE SCIPcliquetableAdd(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR**            vars,               /**< binary variables in the clique: at most one can be set to the given value */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars,              /**< number of variables in the clique */
   SCIP_Bool             isequation,         /**< is the clique an equation or an inequality? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   );

/** removes all empty and single variable cliques from the clique table; removes double entries from the clique table
 *
 * @note cliques can be processed several times by this method
 */
extern
SCIP_RETCODE SCIPcliquetableCleanup(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int*                  nchgbds,            /**< pointer to store number of fixed variables */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   );

/** computes connected components of the clique graph
 *
 *  use depth-first search similarly to the components presolver/constraint handler, representing a clique as a
 *  path to reduce memory usage, but leaving the connected components the same
 *
 *  an update becomes necessary if a clique gets added with variables from different components
 */
extern
SCIP_RETCODE SCIPcliquetableComputeCliqueComponents(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< array of problem variables, sorted by variable type */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars,           /**< number of integer variables */
   int                   nimplvars           /**< number of implicit integer variables */
   );

/** returns the index of the connected component of the clique graph that the variable belongs to, or -1  */
extern
int SCIPcliquetableGetVarComponentIdx(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the number of cliques stored in the clique list */
extern
int SCIPcliquelistGetNCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_Bool             value               /**< value of the variable for which the cliques should be returned */
   );

/** returns the cliques stored in the clique list, or NULL if the clique list is empty */
extern
SCIP_CLIQUE** SCIPcliquelistGetCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_Bool             value               /**< value of the variable for which the cliques should be returned */
   );

/** checks whether variable is contained in all cliques of the cliquelist */
extern
void SCIPcliquelistCheck(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_VAR*             var                 /**< variable, the clique list belongs to */
   );

/** gets the number of cliques stored in the clique table */
extern
int SCIPcliquetableGetNCliques(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** gets the number of cliques created so far by the clique table */
extern
int SCIPcliquetableGetNCliquesCreated(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** gets the array of cliques stored in the clique table */
extern
SCIP_CLIQUE** SCIPcliquetableGetCliques(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** gets the number of entries in the whole clique table */
extern
SCIP_Longint SCIPcliquetableGetNEntries(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** returns the number of clique components, or -1 if update is necessary first */
extern
int SCIPcliquetableGetNCliqueComponents(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** returns TRUE iff the connected clique components need an update (because new cliques were added) */
extern
SCIP_Bool SCIPcliquetableNeedsComponentUpdate(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcliquelistGetNCliques(cliquelist, value) ((cliquelist) != NULL ? (cliquelist)->ncliques[value] : 0)
#define SCIPcliquelistGetCliques(cliquelist, value)  ((cliquelist) != NULL ? (cliquelist)->cliques[value] : NULL)
#define SCIPcliquelistCheck(cliquelist, var)         /**/
#define SCIPcliquetableGetNCliques(cliquetable)      ((cliquetable)->ncliques)
#define SCIPcliquetableGetCliques(cliquetable)       ((cliquetable)->cliques)
#define SCIPcliquetableGetNEntries(cliquetable)      ((cliquetable)->nentries)
#define SCIPcliquetableGetNCliqueComponents(cliquetable) (cliquetable->compsfromscratch ? -1 : cliquetable->ncliquecomponents)
#define SCIPcliquetableNeedsComponentUpdate(cliquetable) (cliquetable->compsfromscratch || cliquetable->djset == NULL)
#endif

#ifdef __cplusplus
}
#endif

#endif
