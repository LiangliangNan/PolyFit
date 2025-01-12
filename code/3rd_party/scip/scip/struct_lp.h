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

/**@file   struct_lp.h
 * @ingroup INTERNALAPI
 * @brief  data structures for LP management
 * @author Tobias Achterberg
 *
 *  In SCIP, the LP is defined as follows:
 *
 *     min       obj * x
 *        lhs <=   A * x + const <= rhs
 *        lb  <=       x         <= ub
 *
 *  The row activities are defined as activity = A * x + const and must
 *  therefore be in the range of [lhs,rhs].
 *
 *  Mathematically, each range constraint would account for two dual
 *  variables, one for each inequality. Since in an optimal solution (at
 *  least) one of them may be chosen to be zero, we may define one dual
 *  multiplier for each row as the difference of those two.
 *
 *  Let y be the vector of dual multipliers for the rows, then the reduced
 *  costs are defined as
 *
 *     redcost = obj - A^T * y.
 *
 *  In an optimal solution, y must be
 *
 *     - nonnegative, if the corresponding row activity is not tight at its rhs
 *     - nonpositive, if the corresponding row activity is not tight at its lhs
 *     - zero, if the corresponding row activity is not at any of its sides
 *
 *  and the reduced costs must be
 *
 *     - nonnegative, if the corresponding variable is not tight at its ub
 *     - nonpositive, if the corresponding variable is not tight at its lb
 *     - zero, if the corresponding variable is not at any of its bounds.
 *
 *  The main datastructures for storing an LP are the rows and the columns.
 *  A row can live on its own (if it was created by a separator), or as SCIP_LP
 *  relaxation of a constraint. Thus, it has a nuses-counter, and is
 *  deleted, if not needed any more.
 *  A column cannot live on its own. It is always connected to a problem
 *  variable. Because pricing is always problem specific, it cannot create
 *  LP columns without introducing new variables. Thus, each column is
 *  connected to exactly one variable, and is deleted, if the variable
 *  is deleted.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_LP_H__
#define __SCIP_STRUCT_LP_H__


#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_event.h"
#include "lpi/type_lpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** collected values of a column which depend on the LP solution
 *  We store these values in each column to recover the LP solution at start of diving or probing mode, say, without
 *  having to resolve the LP.  Note that we do not store the farkascoef value since we do expect a node with infeasible
 *  LP to be pruned anyway.
 */
struct SCIP_ColSolVals
{
   SCIP_Real             primsol;            /**< primal solution value in LP, is 0 if col is not in LP */
   SCIP_Real             redcost;            /**< reduced cost value in LP, or SCIP_INVALID if not yet calculated */
   unsigned int          basisstatus:2;      /**< basis status of column in last LP solution, invalid for non-LP columns */
};

/** collected values of a row which depend on the LP solution
 *  We store these values in each row to recover the LP solution at start of diving or probing mode, say, without having
 *  to resolve the LP.  We do not store the dualfarkas value since we expect a node with infeasible LP to be pruned
 *  anyway. In this unlikely case, we have to resolve the LP.
 */
struct SCIP_RowSolVals
{
   SCIP_Real             dualsol;            /**< dual solution value in LP, is 0 if row is not in LP */
   SCIP_Real             activity;           /**< row activity value in LP, or SCIP_INVALID if not yet calculated */
   unsigned int          basisstatus:2;      /**< basis status of row in last LP solution, invalid for non-LP rows */
};

/** collected values of the LP data which depend on the LP solution
 *  We store these values to recover the LP solution at start of diving or probing mode, say, without having to resolve
 *  the LP.
 */
struct SCIP_LpSolVals
{
   SCIP_LPSOLSTAT        lpsolstat;          /**< solution status of last LP solution */
   SCIP_Real             lpobjval;           /**< objective value of LP without loose variables, or SCIP_INVALID */
   SCIP_Bool             primalfeasible;     /**< is current LP solution primal feasible? */
   SCIP_Bool             primalchecked;      /**< was current LP solution checked for primal feasibility? */
   SCIP_Bool             dualfeasible;       /**< is current LP solution dual feasible? */
   SCIP_Bool             dualchecked;        /**< was current LP solution checked for primal feasibility? */
   SCIP_Bool             solisbasic;         /**< is current LP solution a basic solution? */
   SCIP_Bool             lpissolved;         /**< is current LP solved? */
};

/** LP column;
 *  The row vector of the LP column is partitioned into two parts: The first col->nlprows rows in the rows array
 *  are the ones that belong to the current LP (col->rows[j]->lppos >= 0) and that are linked to the column
 *  (col->linkpos[j] >= 0). The remaining col->len - col->nlprows rows in the rows array are the ones that
 *  don't belong to the current LP (col->rows[j]->lppos == -1) or that are not linked to the column
 *  (col->linkpos[j] == -1).
 */
struct SCIP_Col
{
   SCIP_Real             obj;                /**< current objective value of column in LP (might be changed in diving or probing) */
   SCIP_Real             lb;                 /**< current lower bound of column in LP */
   SCIP_Real             ub;                 /**< current upper bound of column in LP */
   SCIP_Real             unchangedobj;       /**< unchanged objective value of column (ignoring diving or probing changes) */
   SCIP_Real             lazylb;             /**< lazy lower bound of the column; if the current lower bound is not greater than 
                                              *   the lazy lower bound, then the lower bound has not to be added to the LP */
   SCIP_Real             lazyub;             /**< lazy upper bound of the column; if the current upper bound is not smaller than 
                                              *   the lazy upper bound, then the upper bound has not to be added to the LP */
   SCIP_Real             flushedobj;         /**< objective value of column already flushed to the LP solver */
   SCIP_Real             flushedlb;          /**< lower bound of column already flushed to the LP solver */
   SCIP_Real             flushedub;          /**< upper bound of column already flushed to the LP solver */
   SCIP_Real             primsol;            /**< primal solution value in LP, is 0 if col is not in LP */
   SCIP_Real             redcost;            /**< reduced cost value in LP, or SCIP_INVALID if not yet calculated */
   SCIP_Real             farkascoef;         /**< coefficient in dual Farkas infeasibility proof (== dualfarkas^T A_c) */
   SCIP_Real             minprimsol;         /**< minimal LP solution value, this column ever assumed */
   SCIP_Real             maxprimsol;         /**< maximal LP solution value, this column ever assumed */
   SCIP_Real             sbdown;             /**< strong branching information for downwards branching */
   SCIP_Real             sbup;               /**< strong branching information for upwards branching */
   SCIP_Real             sbsolval;           /**< LP solution value of column at last strong branching call */
   SCIP_Real             sblpobjval;         /**< LP objective value at last strong branching call on the column */
   SCIP_Longint          sbnode;             /**< node number of the last strong branching call on this column */
   SCIP_Longint          obsoletenode;       /**< last node where this column was removed due to aging */
   SCIP_COLSOLVALS*      storedsolvals;      /**< values stored before entering diving or probing mode */
   SCIP_VAR*             var;                /**< variable, this column represents; there cannot be a column without variable */
   SCIP_ROW**            rows;               /**< rows of column entries, that may have a nonzero dual solution value */
   SCIP_Real*            vals;               /**< coefficients of column entries */
   SCIP_Longint          validredcostlp;     /**< LP number for which reduced cost value is valid */
   SCIP_Longint          validfarkaslp;      /**< LP number for which Farkas coefficient is valid */
   SCIP_Longint          validsblp;          /**< LP number for which strong branching values are valid */
   int*                  linkpos;            /**< position of col in col vector of the row, or -1 if not yet linked */
   int                   index;              /**< consecutively numbered column identifier */
   int                   size;               /**< size of the row- and val-arrays */
   int                   len;                /**< number of nonzeros in column */
   int                   nlprows;            /**< number of linked rows in column, that belong to the current LP */
   int                   nunlinked;          /**< number of column entries, where the rows don't know about the column */
   int                   lppos;              /**< column position number in current LP, or -1 if not in current LP */
   int                   lpipos;             /**< column position number in LP solver, or -1 if not in LP solver */
   int                   lpdepth;            /**< depth level at which column entered the LP, or -1 if not in current LP */
   int                   sbitlim;            /**< strong branching iteration limit used to get strong branching values, or -1 */
   int                   nsbcalls;           /**< number of times, strong branching was applied on the column */
   int                   age;                /**< number of successive times this variable was in LP and was 0.0 in solution */
   int                   var_probindex;      /**< copy of var->probindex for avoiding expensive dereferencing */
   unsigned int          basisstatus:2;      /**< basis status of column in last LP solution, invalid for non-LP columns */
   unsigned int          lprowssorted:1;     /**< are the linked LP rows in the rows array sorted by non-decreasing index? */
   unsigned int          nonlprowssorted:1;  /**< are the non-LP/not linked rows sorted by non-decreasing index? */
   unsigned int          objchanged:1;       /**< has objective value changed, and has data of LP solver to be updated? */
   unsigned int          lbchanged:1;        /**< has lower bound changed, and has data of LP solver to be updated? */
   unsigned int          ubchanged:1;        /**< has upper bound changed, and has data of LP solver to be updated? */
   unsigned int          coefchanged:1;      /**< has the coefficient vector changed, and has LP solver to be updated? */
   unsigned int          integral:1;         /**< is associated variable of integral type? */
   unsigned int          removable:1;        /**< is column removable from the LP (due to aging or cleanup)? */
   unsigned int          sbdownvalid:1;      /**< stores whether the stored strong branching down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   unsigned int          sbupvalid:1;        /**< stores whether the stored strong branching up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
};

/** LP row
 *  The column vector of the LP row is partitioned into two parts: The first row->nlpcols columns in the cols array
 *  are the ones that belong to the current LP (row->cols[j]->lppos >= 0) and that are linked to the row
 *  (row->linkpos[j] >= 0). The remaining row->len - row->nlpcols columns in the cols array are the ones that
 *  don't belong to the current LP (row->cols[j]->lppos == -1) or that are not linked to the row
 *  (row->linkpos[j] == -1).
 */
struct SCIP_Row
{
   SCIP_Real             constant;           /**< constant shift c in row lhs <= ax + c <= rhs */
   SCIP_Real             lhs;                /**< left hand side of row */
   SCIP_Real             rhs;                /**< right hand side of row */
   SCIP_Real             flushedlhs;         /**< left hand side minus constant of row already flushed to the LP solver */
   SCIP_Real             flushedrhs;         /**< right hand side minus constant of row already flushed to the LP solver */
   SCIP_Real             sqrnorm;            /**< squared Euclidean norm of row vector */
   SCIP_Real             sumnorm;            /**< sum norm of row vector (sum of absolute values of coefficients) */
   SCIP_Real             objprod;            /**< scalar product of row vector with objective function */
   SCIP_Real             maxval;             /**< maximal absolute value of row vector, only valid if nummaxval > 0 */
   SCIP_Real             minval;             /**< minimal absolute non-zero value of row vector, only valid if numminval > 0 */
   SCIP_Real             dualsol;            /**< dual solution value in LP, is 0 if row is not in LP */
   SCIP_Real             activity;           /**< row activity value in LP, or SCIP_INVALID if not yet calculated */
   SCIP_Real             dualfarkas;         /**< multiplier value in dual Farkas infeasibility proof */
   SCIP_Real             pseudoactivity;     /**< row activity value in pseudo solution, or SCIP_INVALID if not yet calculated */
   SCIP_Real             minactivity;        /**< minimal activity value w.r.t. the column's bounds, or SCIP_INVALID */
   SCIP_Real             maxactivity;        /**< maximal activity value w.r.t. the column's bounds, or SCIP_INVALID */
   SCIP_Longint          validpsactivitydomchg; /**< domain change number for which pseudo activity value is valid */
   SCIP_Longint          validactivitybdsdomchg;/**< domain change number for which activity bound values are valid */
   SCIP_Longint          obsoletenode;       /**< last node where this row was removed due to aging */
   SCIP_Longint          activeinlpcounter;  /**< counter for the number of times this row was active in an optimal LP solution */
   SCIP_Longint          nlpsaftercreation;  /**< counter for the number of LPs after the row has been created */
   SCIP_ROWSOLVALS*      storedsolvals;      /**< values stored before entering diving or probing mode */
   void*                 origin;             /**< pointer to constraint handler or separator who created the row (NULL if unknown) */
   char*                 name;               /**< name of the row */
   SCIP_COL**            cols;               /**< columns of row entries, that may have a nonzero primal solution value */
   int*                  cols_index;         /**< copy of cols[i]->index for avoiding expensive dereferencing */
   SCIP_Real*            vals;               /**< coefficients of row entries */
   int*                  linkpos;            /**< position of row in row vector of the column, or -1 if not yet linked */
   SCIP_EVENTFILTER*     eventfilter;        /**< event filter for events concerning this row */
   SCIP_Longint          validactivitylp;    /**< LP number for which activity value is valid */
   int                   index;              /**< consecutively numbered row identifier */
   int                   size;               /**< size of the col- and val-arrays */
   int                   len;                /**< number of nonzeros in row */
   int                   nlpcols;            /**< number of linked columns in row, that belong to the current LP */
   int                   nunlinked;          /**< number of row entries, where the columns don't know about the row */
   int                   nuses;              /**< number of times, this row is referenced */
   int                   lppos;              /**< row position number in current LP, or -1 if not in current LP */
   int                   lpipos;             /**< row position number in LP solver, or -1 if not in LP solver */
   int                   lpdepth;            /**< depth level at which row entered the LP, or -1 if not in current LP */
   int                   minidx;             /**< minimal column index of row entries */
   int                   maxidx;             /**< maximal column index of row entries */
   int                   numintcols;         /**< number of integral columns */
   int                   nummaxval;          /**< number of coefs with absolute value equal to maxval, zero if maxval invalid */
   int                   numminval;          /**< number of coefs with absolute value equal to minval, zero if minval invalid */
   int                   age;                /**< number of successive times this row was in LP and was not sharp in solution */
   int                   rank;               /**< rank of the row (upper bound, to be precise) */
   unsigned int          fromcutpool:1;      /**< added from cutpool to sepastore */
   unsigned int          basisstatus:2;      /**< basis status of row in last LP solution, invalid for non-LP rows */
   unsigned int          lpcolssorted:1;     /**< are the linked LP columns in the cols array sorted by non-decreasing index? */
   unsigned int          nonlpcolssorted:1;  /**< are the non-LP/not linked columns sorted by non-decreasing index? */
   unsigned int          delaysort:1;        /**< should the row sorting be delayed and done in a lazy fashion? */
   unsigned int          validminmaxidx:1;   /**< are minimal and maximal column index valid? */
   unsigned int          lhschanged:1;       /**< was left hand side or constant changed, and has LP solver to be updated? */
   unsigned int          rhschanged:1;       /**< was right hand side or constant changed, and has LP solver to be updated? */
   unsigned int          coefchanged:1;      /**< was the coefficient vector changed, and has LP solver to be updated? */
   unsigned int          integral:1;         /**< is activity (without constant) of row always integral in feasible solution? */
   unsigned int          local:1;            /**< is row only valid locally? */
   unsigned int          modifiable:1;       /**< is row modifiable during node processing (subject to column generation)? */
   unsigned int          removable:1;        /**< is row removable from the LP (due to aging or cleanup)? */
   unsigned int          inglobalcutpool:1;  /**< is row contained in the global cut pool? */
   unsigned int          normunreliable:1;   /**< is the objective product of the row unreliable? */
   unsigned int          nlocks:13;          /**< number of sealed locks of an unmodifiable row */
   unsigned int          origintype:3;       /**< origin of row (0: unknown, 1: constraint handler, 2: constraint, 3: separator, 4: reoptimization) */
};

/** current LP data */
struct SCIP_Lp
{
   SCIP_Real             lpobjval;           /**< objective value of LP without loose variables, or SCIP_INVALID */
   SCIP_Real             looseobjval;        /**< current solution value of all loose variables set to their best bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Real             rellooseobjval;     /**< last reliable solution value of all loose variables set to their best bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Real             glbpseudoobjval;    /**< global pseudo solution value with all variables set to their best global bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Real             relglbpseudoobjval; /**< last reliable global pseudo solution value */
   SCIP_Real             pseudoobjval;       /**< current pseudo solution value with all variables set to their best bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Real             relpseudoobjval;    /**< last reliable pseudo solution value */
   SCIP_Real             rootlpobjval;       /**< objective value of root LP without loose variables, or SCIP_INVALID */
   SCIP_Real             rootlooseobjval;    /**< objective value of loose variables in root node, or SCIP_INVALID */
   SCIP_Real             cutoffbound;        /**< upper objective limit of LP (copy of primal->cutoffbound) */
   SCIP_Real             feastol;            /**< current feasibility tolerance */
   SCIP_Real             lpiobjlim;          /**< current objective limit in LPI */
   SCIP_Real             lpifeastol;         /**< current feasibility tolerance in LPI */
   SCIP_Real             lpidualfeastol;     /**< current reduced costs feasibility tolerance in LPI */
   SCIP_Real             lpibarrierconvtol;  /**< current convergence tolerance used in barrier algorithm in LPI */
   SCIP_Real             lpiconditionlimit;  /**< current condition number limit in LPI */
   SCIP_Real             lpimarkowitz;       /**< current markowitz threshold */
   SCIP_Real             objsqrnorm;         /**< squared Euclidean norm of objective function vector of problem variables */
   SCIP_Real             objsumnorm;         /**< sum norm of objective function vector of problem variables */
   SCIP_Real             degeneracy;         /**< share of degenerate non-basic variables in the current LP */
   SCIP_Real             varconsratio;       /**< variable-constraint ratio of the optimal face */
   SCIP_LPI*             lpi;                /**< LP solver interface */
   SCIP_COL**            lpicols;            /**< array with columns currently stored in the LP solver */
   SCIP_ROW**            lpirows;            /**< array with rows currently stored in the LP solver */
   SCIP_COL**            chgcols;            /**< array of changed columns not yet applied to the LP solver */
   SCIP_ROW**            chgrows;            /**< array of changed rows not yet applied to the LP solver */
   SCIP_COL**            cols;               /**< array with current LP columns in correct order */
   SCIP_COL**            lazycols;           /**< array with current LP lazy columns */
   SCIP_ROW**            rows;               /**< array with current LP rows in correct order */
   SCIP_Real*            soldirection;       /**< normalized vector in direction of primal solution from current LP solution */
   SCIP_LPISTATE*        divelpistate;       /**< stores LPI state (basis information) before diving starts */
   SCIP_Real*            divechgsides;       /**< stores the lhs/rhs changed in the current diving */
   SCIP_SIDETYPE*        divechgsidetypes;   /**< stores the side type of the changes done in the current diving */
   SCIP_ROW**            divechgrows;        /**< stores the rows changed in the current diving */
   SCIP_LPSOLVALS*       storedsolvals;      /**< collected values of the LP data which depend on the LP solution */
   SCIP_SOL*             validsoldirsol;     /**< primal solution for which the currently stored solution direction vector is valid */
   SCIP_Longint          validsollp;         /**< LP number for which the currently stored solution values are valid */
   SCIP_Longint          validfarkaslp;      /**< LP number for which the currently stored Farkas row multipliers are valid */
   SCIP_Longint          validsoldirlp;      /**< LP number for which the currently stored solution direction vector is valid */
   SCIP_Longint          validdegeneracylp;  /**< LP number for which the currently stored degeneracy information is valid */
   SCIP_Longint          divenolddomchgs;    /**< number of domain changes before diving has started */
   int                   lpicolssize;        /**< available slots in lpicols vector */
   int                   nlpicols;           /**< number of columns in the LP solver */
   int                   lpifirstchgcol;     /**< first column of the LP which differs from the column in the LP solver */
   int                   lpirowssize;        /**< available slots in lpirows vector */
   int                   nlpirows;           /**< number of rows in the LP solver */
   int                   lpifirstchgrow;     /**< first row of the LP which differs from the row in the LP solver */
   int                   chgcolssize;        /**< available slots in chgcols vector */
   int                   nchgcols;           /**< current number of chgcols (number of used slots in chgcols vector) */
   int                   chgrowssize;        /**< available slots in chgrows vector */
   int                   nchgrows;           /**< current number of chgrows (number of used slots in chgrows vector) */
   int                   colssize;           /**< available slots in cols vector */
   int                   soldirectionsize;   /**< available slots in soldirection vector */
   int                   ncols;              /**< current number of LP columns (number of used slots in cols vector) */
   int                   lazycolssize;       /**< available slots in lazycols vector */
   int                   nlazycols;          /**< current number of LP lazy columns (number of used slots in lazycols vector) */
   int                   nremovablecols;     /**< number of removable columns in the LP */
   int                   firstnewcol;        /**< first column added at the current node */
   int                   rowssize;           /**< available slots in rows vector */
   int                   nrows;              /**< current number of LP rows (number of used slots in rows vector) */
   int                   nremovablerows;     /**< number of removable rows in the LP */
   int                   firstnewrow;        /**< first row added at the current node */
   int                   looseobjvalinf;     /**< number of loose variables with infinite best bound in current solution */
   int                   nloosevars;         /**< number of loose variables in LP */
   int                   glbpseudoobjvalinf; /**< number of variables with infinite best bound in global pseudo solution */
   int                   pseudoobjvalinf;    /**< number of variables with infinite best bound in current pseudo solution */
   int                   ndivingrows;        /**< number of rows when entering diving mode */
   int                   ndivechgsides;      /**< number of side changes in current diving */
   int                   divechgsidessize;   /**< size of the arrays */
   int                   divinglpiitlim;     /**< LPI iteration limit when entering diving mode */
   int                   lpiitlim;           /**< current iteration limit setting in LPI */
   int                   lpifastmip;         /**< current FASTMIP setting in LPI */
   int                   lpithreads;         /**< current THREADS setting in LPI */
   int                   lpitiming;          /**< current timing type in LPI */
   int                   lpirandomseed;      /**< current initial random seed in LPI */
   int                   lpiscaling;         /**< current SCALING setting in LPI */
   int                   lpirefactorinterval;/**< current refactorization interval */
   SCIP_PRICING          lpipricing;         /**< current pricing setting in LPI */
   SCIP_LPSOLSTAT        lpsolstat;          /**< solution status of last LP solution */
   SCIP_LPALGO           lastlpalgo;         /**< algorithm used for last LP solve */
   SCIP_Bool             objsqrnormunreliable;/**< is squared Euclidean norm of objective function vector of problem
                                               *   variables unreliable and need recalculation? */
   SCIP_Bool             lpisolutionpolishing;/**< LP solution polishing method (0: disabled, 1: enabled) */
   SCIP_Bool             looseobjvalid;      /**< is the loose objective value valid or should it be recomputed from scratch? */
   SCIP_Bool             glbpseudoobjvalid;  /**< is the global pseudo solution value valid or should it be recomputed from scratch? */
   SCIP_Bool             pseudoobjvalid;     /**< is the pseudo solution value valid or should it be recomputed from scratch? */
   SCIP_Bool             flushdeletedcols;   /**< have LPI-columns been deleted in the last lpFlush() call? */
   SCIP_Bool             flushaddedcols;     /**< have LPI-columns been added in the last lpFlush() call? */
   SCIP_Bool             flushdeletedrows;   /**< have LPI-rows been deleted in the last lpFlush() call? */
   SCIP_Bool             flushaddedrows;     /**< have LPI-rows been added in the last lpFlush() call? */
   SCIP_Bool             updateintegrality;  /**< does integrality information need to be updated? */
   SCIP_Bool             flushed;            /**< are all cached changes applied to the LP solver? */
   SCIP_Bool             solved;             /**< is current LP solved? */
   SCIP_Bool             primalfeasible;     /**< is current LP solution (rather LPI state) primal feasible? */
   SCIP_Bool             primalchecked;      /**< was current LP solution checked for primal feasibility?? */
   SCIP_Bool             dualfeasible;       /**< is current LP solution (rather LPI state) dual feasible? */
   SCIP_Bool             dualchecked;        /**< was current LP solution checked for primal feasibility?? */
   SCIP_Bool             solisbasic;         /**< is current LP solution a basic solution? */
   SCIP_Bool             rootlpisrelax;      /**< is root LP a relaxation of the problem and its solution value a valid global lower bound? */
   SCIP_Bool             isrelax;            /**< is the current LP a relaxation of the problem for which it has been solved and its 
                                              *   solution value a valid local lower bound? */
   SCIP_Bool             installing;         /**< whether the solution process is in stalling */
   SCIP_Bool             strongbranching;    /**< whether the lp is used for strong branching */
   SCIP_Bool             probing;            /**< are we currently in probing mode? */
   SCIP_Bool             strongbranchprobing;/**< are we currently in probing mode for strong branching? */
   SCIP_Bool             diving;             /**< LP is used for diving: col bounds and obj don't correspond to variables */
   SCIP_Bool             divingobjchg;       /**< objective values were changed in diving or probing: LP objective is invalid */
   SCIP_Bool             divinglazyapplied;  /**< lazy bounds were applied to the LP during diving */
   SCIP_Bool             resolvelperror;     /**< an error occurred during resolving the LP after diving or probing */
   SCIP_Bool             adjustlpval;        /**< does an infinite LP objective value has been adjusted so far? */
   SCIP_Bool             lpifromscratch;     /**< current FROMSCRATCH setting in LPI */
   SCIP_Bool             lpipresolving;      /**< current PRESOLVING setting in LPI */
   SCIP_Bool             lpilpinfo;          /**< current LPINFO setting in LPI */
   SCIP_Bool             lpihasfeastol;      /**< does the LPI support the FEASTOL parameter? */
   SCIP_Bool             lpihasdualfeastol;  /**< does the LPI support the DUALFEASTOL parameter? */
   SCIP_Bool             lpihasbarrierconvtol;/**< does the LPI support the BARRIERCONVTOL parameter? */
   SCIP_Bool             lpihasfastmip;      /**< does the LPI support the FASTMIP parameter? */
   SCIP_Bool             lpihasscaling;      /**< does the LPI support the SCALING parameter? */
   SCIP_Bool             lpihaspresolving;   /**< does the LPI support the PRESOLVING parameter? */
   SCIP_Bool             lpihasrowrep;       /**< does the LPI support row representation of a simplex basis? */
   SCIP_Bool             lpihaspolishing;    /**< does the LPI support solution polishing? */
   SCIP_Bool             lpihasrefactor;     /**< does the LPI support changing the refactorization interval? */
   SCIP_Real             lpirowrepswitch;    /**< simplex algorithm shall use row representation of the basis
                                              *   if number of rows divided by number of columns exceeds this value */
   SCIP_Bool             divelpwasprimfeas;  /**< primal feasibility when diving started */
   SCIP_Bool             divelpwasprimchecked;/**< primal feasibility was checked when diving started */
   SCIP_Bool             divelpwasdualfeas;  /**< dual feasibility when diving started */
   SCIP_Bool             divelpwasdualchecked;/**< dual feasibility was checked when diving started */
};

#ifdef __cplusplus
}
#endif

#endif
