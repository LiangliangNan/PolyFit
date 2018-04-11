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

/**@file   heur_shiftandpropagate.c
 * @brief  shiftandpropagate primal heuristic
 * @author Timo Berthold
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/pub_misc.h"
#include "scip/heur_shiftandpropagate.h"

#define HEUR_NAME             "shiftandpropagate"
#define HEUR_DESC             "Pre-root heuristic to expand an auxiliary branch-and-bound tree and apply propagation techniques"
#define HEUR_DISPCHAR         'T'
#define HEUR_PRIORITY         1000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE     /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_WEIGHT_INEQUALITY   1   /**< the heuristic row weight for inequalities */
#define DEFAULT_WEIGHT_EQUALITY     3   /**< the heuristic row weight for equations */
#define DEFAULT_RELAX            TRUE   /**< Should continuous variables be relaxed from the problem? */
#define DEFAULT_PROBING          TRUE   /**< Is propagation of solution values enabled? */
#define DEFAULT_ONLYWITHOUTSOL   TRUE   /**< Should heuristic only be executed if no primal solution was found, yet? */
#define DEFAULT_NPROPROUNDS        10   /**< The default number of propagation rounds for each propagation used */
#define DEFAULT_PROPBREAKER     65000   /**< fixed maximum number of propagations */
#define DEFAULT_CUTOFFBREAKER      15   /**< fixed maximum number of allowed cutoffs before the heuristic stops */
#define DEFAULT_RANDSEED           29   /**< the default random seed for random number generation */
#define DEFAULT_SORTKEY            'v'  /**< the default key for variable sorting */
#define DEFAULT_SORTVARS         TRUE   /**< should variables be processed in sorted order? */
#define DEFAULT_COLLECTSTATS     TRUE   /**< should variable statistics be collected during probing? */
#define DEFAULT_STOPAFTERFEASIBLE TRUE  /**< Should the heuristic stop calculating optimal shift values when no more rows are violated? */
#define DEFAULT_PREFERBINARIES   TRUE   /**< Should binary variables be shifted first? */
#define DEFAULT_SELECTBEST       FALSE  /**< should the heuristic choose the best candidate in every round? (set to FALSE for static order)? */
#define DEFAULT_MAXCUTOFFQUOT      0.0  /**< maximum percentage of allowed cutoffs before stopping the heuristic */
#define SORTKEYS                 "nrtuv"/**< options sorting key: (n)orms down, norms (u)p, (v)iolated rows decreasing,
                                         *   viola(t)ed rows increasing, or (r)andom */
#define DEFAULT_NOZEROFIXING      FALSE /**< should variables with a zero shifting value be delayed instead of being fixed? */
#define DEFAULT_FIXBINLOCKS       TRUE  /**< should binary variables with no locks in one direction be fixed to that direction? */
#define DEFAULT_BINLOCKSFIRST     FALSE  /**< should binary variables with no locks be preferred in the ordering? */
#define DEFAULT_NORMALIZE         TRUE  /**< should coefficients and left/right hand sides be normalized by max row coeff? */
#define DEFAULT_UPDATEWEIGHTS     FALSE /**< should row weight be increased every time the row is violated? */
#define DEFAULT_IMPLISCONTINUOUS   TRUE /**< should implicit integer variables be treated as continuous variables? */

#define EVENTHDLR_NAME         "eventhdlrshiftandpropagate"
#define EVENTHDLR_DESC         "event handler to catch bound changes"
#define EVENTTYPE_SHIFTANDPROPAGATE (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_GBDCHANGED)


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_COL**            lpcols;             /**< stores lp columns with discrete variables before cont. variables */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generation */
   int*                  rowweights;         /**< row weight storage */
   SCIP_Bool             relax;              /**< should continuous variables be relaxed from the problem */
   SCIP_Bool             probing;            /**< should probing be executed? */
   SCIP_Bool             onlywithoutsol;     /**< Should heuristic only be executed if no primal solution was found, yet? */
   int                   nlpcols;            /**< the number of lp columns */
   int                   nproprounds;        /**< The default number of propagation rounds for each propagation used */
   int                   cutoffbreaker;      /**< the number of cutoffs before heuristic execution is stopped, or -1 for no
                                               * limit */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler to register and process variable bound changes */

   SCIP_Real             maxcutoffquot;      /**< maximum percentage of allowed cutoffs before stopping the heuristic */
   char                  sortkey;            /**< the key by which variables are sorted */
   SCIP_Bool             sortvars;           /**< should variables be processed in sorted order? */
   SCIP_Bool             collectstats;       /**< should variable statistics be collected during probing? */
   SCIP_Bool             stopafterfeasible;  /**< Should the heuristic stop calculating optimal shift values when no
                                              *   more rows are violated? */
   SCIP_Bool             preferbinaries;     /**< Should binary variables be shifted first? */
   SCIP_Bool             nozerofixing;       /**< should variables with a zero shifting value be delayed instead of being fixed? */
   SCIP_Bool             fixbinlocks;        /**< should binary variables with no locks in one direction be fixed to that direction? */
   SCIP_Bool             binlocksfirst;      /**< should binary variables with no locks be preferred in the ordering? */
   SCIP_Bool             normalize;          /**< should coefficients and left/right hand sides be normalized by max row coeff? */
   SCIP_Bool             updateweights;      /**< should row weight be increased every time the row is violated? */
   SCIP_Bool             impliscontinuous;   /**< should implicit integer variables be treated as continuous variables? */
   SCIP_Bool             selectbest;         /**< should the heuristic choose the best candidate in every round? (set to FALSE for static order)? */
   SCIPstatistic(
      SCIP_LPSOLSTAT     lpsolstat;          /**< the probing status after probing */
      SCIP_Longint       ntotaldomredsfound; /**< the total number of domain reductions during heuristic */
      SCIP_Longint       nlpiters;           /**< number of LP iterations which the heuristic needed */
      int                nremainingviols;    /**< the number of remaining violations */
      int                nprobings;          /**< how many probings has the heuristic executed? */
      int                ncutoffs;           /**< has the probing node been cutoff? */
      )
};

/** status of a variable in heuristic transformation */
enum TransformStatus
{
   TRANSFORMSTATUS_NONE = 0,            /**< variable has not been transformed yet */
   TRANSFORMSTATUS_LB   = 1,            /**< variable has been shifted by using lower bound (x-lb) */
   TRANSFORMSTATUS_NEG  = 2,            /**< variable has been negated by using upper bound (ub-x) */
   TRANSFORMSTATUS_FREE = 3             /**< variable does not have to be shifted */
};
typedef enum TransformStatus TRANSFORMSTATUS;

/** information about the matrix after its heuristic transformation */
struct ConstraintMatrix
{
   SCIP_Real*            rowmatvals;         /**< matrix coefficients row by row */
   int*                  rowmatind;          /**< the indices of the corresponding variables */
   int*                  rowmatbegin;        /**< the starting indices of each row */
   SCIP_Real*            colmatvals;         /**< matrix coefficients column by column */
   int*                  colmatind;          /**< the indices of the corresponding rows for each coefficient */
   int*                  colmatbegin;        /**< the starting indices of each column */
   int*                  violrows;           /**< the number of violated rows for every variable */
   TRANSFORMSTATUS*      transformstatus;    /**< information about transform status of every discrete variable */
   SCIP_Real*            lhs;                /**< left hand side vector after normalization */
   SCIP_Real*            rhs;                /**< right hand side vector after normalization */
   SCIP_Real*            colnorms;           /**< vector norms of all discrete problem variables after normalization */
   SCIP_Real*            upperbounds;        /**< the upper bounds of every non-continuous variable after transformation*/
   SCIP_Real*            transformshiftvals; /**< values by which original discrete variable bounds were shifted */
   int                   nnonzs;             /**< number of nonzero column entries */
   int                   nrows;              /**< number of rows of matrix */
   int                   ncols;              /**< the number of columns in matrix (including continuous vars) */
   int                   ndiscvars;          /**< number of discrete problem variables */
   SCIP_Bool             normalized;         /**< indicates if the matrix data has already been normalized */
};
typedef struct ConstraintMatrix CONSTRAINTMATRIX;

struct SCIP_EventhdlrData
{
   CONSTRAINTMATRIX*    matrix;              /**< the constraint matrix of the heuristic */
   SCIP_HEURDATA*       heurdata;            /**< heuristic data */
   int*                 violatedrows;        /**< all currently violated LP rows */
   int*                 violatedrowpos;      /**< position in violatedrows array for every row */
   int*                 nviolatedrows;       /**< pointer to the total number of currently violated rows */
};

struct SCIP_EventData
{
   int                  colpos;              /**< column position of the event-related variable */
};
/*
 * Local methods
 */

/** returns whether a given variable is counted as discrete, depending on the parameter impliscontinuous */
static
SCIP_Bool varIsDiscrete(
   SCIP_VAR*             var,                /**< variable to check for discreteness */
   SCIP_Bool             impliscontinuous    /**< should implicit integer variables be counted as continuous? */
   )
{
   return SCIPvarIsIntegral(var) && (SCIPvarGetType(var) != SCIP_VARTYPE_IMPLINT || !impliscontinuous);
}

/** returns whether a given column is counted as discrete, depending on the parameter impliscontinuous */
static
SCIP_Bool colIsDiscrete(
   SCIP_COL*             col,                /**< column to check for discreteness */
   SCIP_Bool             impliscontinuous    /**< should implicit integer variables be counted as continuous? */
   )
{
   return SCIPcolIsIntegral(col) && (!impliscontinuous || SCIPvarGetType(SCIPcolGetVar(col)) != SCIP_VARTYPE_IMPLINT);
}

/** returns nonzero values and corresponding columns of given row */
static
void getRowData(
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   rowindex,           /**< index of the desired row */
   SCIP_Real**           valpointer,         /**< pointer to store the nonzero coefficients of the row */
   SCIP_Real*            lhs,                /**< lhs of the row */
   SCIP_Real*            rhs,                /**< rhs of the row */
   int**                 indexpointer,       /**< pointer to store column indices which belong to the nonzeros */
   int*                  nrowvals            /**< pointer to store number of nonzeros in the desired row (or NULL) */
   )
{
   int arrayposition;

   assert(matrix != NULL);
   assert(0 <= rowindex && rowindex < matrix->nrows);

   arrayposition = matrix->rowmatbegin[rowindex];

   if ( nrowvals != NULL )
   {
      if( rowindex == matrix->nrows - 1 )
         *nrowvals = matrix->nnonzs - arrayposition;
      else
         *nrowvals = matrix->rowmatbegin[rowindex + 1] - arrayposition; /*lint !e679*/
   }

   if( valpointer != NULL )
      *valpointer = &(matrix->rowmatvals[arrayposition]);
   if( indexpointer != NULL )
      *indexpointer = &(matrix->rowmatind[arrayposition]);

   if( lhs != NULL )
      *lhs = matrix->lhs[rowindex];

   if( rhs != NULL )
      *rhs = matrix->rhs[rowindex];
}

/** returns nonzero values and corresponding rows of given column */
static
void getColumnData(
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   colindex,           /**< the index of the desired column */
   SCIP_Real**           valpointer,         /**< pointer to store the nonzero coefficients of the column */
   int**                 indexpointer,       /**< pointer to store row indices which belong to the nonzeros */
   int*                  ncolvals            /**< pointer to store number of nonzeros in the desired column */
   )
{
   int arrayposition;

   assert(matrix != NULL);
   assert(0 <= colindex && colindex < matrix->ncols);

   arrayposition = matrix->colmatbegin[colindex];

   if( ncolvals != NULL )
   {
      if( colindex == matrix->ncols - 1 )
         *ncolvals = matrix->nnonzs - arrayposition;
      else
         *ncolvals = matrix->colmatbegin[colindex + 1] - arrayposition; /*lint !e679*/
   }
   if( valpointer != NULL )
      *valpointer = &(matrix->colmatvals[arrayposition]);

   if( indexpointer != NULL )
      *indexpointer = &(matrix->colmatind[arrayposition]);
}

/** relaxes a continuous variable from all its rows, which has influence
 *  on both the left and right hand side of the constraint.
 */
static
void relaxVar(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_VAR*             var,                /**< variable which is relaxed from the problem */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   SCIP_Bool             normalize           /**< should coefficients and be normalized by rows maximum norms? */
   )
{
   SCIP_ROW** colrows;
   SCIP_COL* varcol;
   SCIP_Real* colvals;
   SCIP_Real ub;
   SCIP_Real lb;
   int ncolvals;
   int r;

   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   varcol = SCIPvarGetCol(var);
   assert(varcol != NULL);

   /* get nonzero values and corresponding rows of variable */
   colvals = SCIPcolGetVals(varcol);
   ncolvals = SCIPcolGetNLPNonz(varcol);
   colrows = SCIPcolGetRows(varcol);

   ub = SCIPvarGetUbGlobal(var);
   lb = SCIPvarGetLbGlobal(var);

   assert(colvals != NULL || ncolvals == 0);

   SCIPdebugMsg(scip, "Relaxing variable <%s> with lb <%g> and ub <%g>\n",
      SCIPvarGetName(var), lb, ub);

   assert(matrix->normalized);
   /* relax variable from all its constraints */
   for( r = 0; r < ncolvals; ++r )
   {
      SCIP_ROW* colrow;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real lhsvarbound;
      SCIP_Real rhsvarbound;
      SCIP_Real rowabs;
      SCIP_Real colval;
      int rowindex;

      colrow = colrows[r];
      rowindex = SCIProwGetLPPos(colrow);

      if( rowindex == -1 )
         break;

      rowabs = SCIPgetRowMaxCoef(scip, colrow);
      assert(colvals != NULL); /* to please flexelint */
      colval = colvals[r];
      if( normalize && SCIPisFeasGT(scip, rowabs, 0.0) )
         colval /= rowabs;

      assert(0 <= rowindex && rowindex < matrix->nrows);
      getRowData(matrix, rowindex, NULL, &lhs, &rhs, NULL, NULL);
      /* variables bound influence the lhs and rhs of current row depending on the sign
       * of the variables coefficient.
       */
      if( SCIPisFeasPositive(scip, colval) )
      {
         lhsvarbound = ub;
         rhsvarbound = lb;
      }
      else if( SCIPisFeasNegative(scip, colval) )
      {
         lhsvarbound = lb;
         rhsvarbound = ub;
      }
      else
         continue;

      /* relax variable from the current row */
      if( !SCIPisInfinity(scip, -matrix->lhs[rowindex]) && !SCIPisInfinity(scip, ABS(lhsvarbound)) )
         matrix->lhs[rowindex] -= colval * lhsvarbound;
      else
         matrix->lhs[rowindex] = -SCIPinfinity(scip);

      if( !SCIPisInfinity(scip, matrix->rhs[rowindex]) && !SCIPisInfinity(scip, ABS(rhsvarbound)) )
         matrix->rhs[rowindex] -= colval * rhsvarbound;
      else
         matrix->rhs[rowindex] = SCIPinfinity(scip);

      SCIPdebugMsg(scip, "Row <%s> changed:Coefficient <%g>, LHS <%g> --> <%g>, RHS <%g> --> <%g>\n",
         SCIProwGetName(colrow), colval, lhs, matrix->lhs[rowindex], rhs, matrix->rhs[rowindex]);
   }
}

/** transforms bounds of a given variable s.t. its lower bound equals zero afterwards.
 *  If the variable already has lower bound zero, the variable is not transformed,
 *  if not, the variable's bounds are changed w.r.t. the smaller absolute value of its
 *  bounds in order to avoid numerical inaccuracies. If both lower and upper bound
 *  of the variable differ from infinity, there are two cases. If |lb| <= |ub|,
 *  the bounds are shifted by -lb, else a new variable ub - x replaces x.
 *  The transformation is memorized by the transform status of the variable s.t.
 *  retransformation is possible.
 */
static
void transformVariable(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   colpos              /**< position of variable column in matrix */
   )
{
   SCIP_COL* col;
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;

   SCIP_Bool negatecoeffs; /* do the row coefficients need to be negated? */
   SCIP_Real deltashift;   /* difference from previous transformation */

   assert(matrix != NULL);
   assert(0 <= colpos && colpos < heurdata->nlpcols);
   col = heurdata->lpcols[colpos];
   assert(col != NULL);
   assert(SCIPcolIsInLP(col));

   var = SCIPcolGetVar(col);
   assert(var != NULL);
   assert(SCIPvarIsIntegral(var));
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   negatecoeffs = FALSE;
   /* if both lower and upper bound are -infinity and infinity, resp., this is reflected by a free transform status.
    * If the lower bound is already zero, this is reflected by identity transform status. In both cases, none of the
    * corresponding rows needs to be modified.
    */
   if( SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
   {
      if( matrix->transformstatus[colpos] == TRANSFORMSTATUS_NEG )
         negatecoeffs = TRUE;

      deltashift = matrix->transformshiftvals[colpos];
      matrix->transformshiftvals[colpos] = 0.0;
      matrix->transformstatus[colpos] = TRANSFORMSTATUS_FREE;
   }
   else if( SCIPisFeasLE(scip, ABS(lb), ABS(ub)) )
   {
      assert(!SCIPisInfinity(scip, lb));
      matrix->transformstatus[colpos] = TRANSFORMSTATUS_LB;
      deltashift = lb;
      matrix->transformshiftvals[colpos] = lb;
   }
   else
   {
      assert(!SCIPisInfinity(scip, ub));
      if( matrix->transformstatus[colpos] != TRANSFORMSTATUS_NEG )
         negatecoeffs = TRUE;
      matrix->transformstatus[colpos] = TRANSFORMSTATUS_NEG;
      deltashift = ub;
      matrix->transformshiftvals[colpos] = ub;
   }

   /* determine the upper bound for this variable in heuristic transformation (lower bound is implicit; always 0) */
   if( !SCIPisInfinity(scip, ub) && !SCIPisInfinity(scip, lb) )
      matrix->upperbounds[colpos] = ub - lb;
   else
      matrix->upperbounds[colpos] = SCIPinfinity(scip);

   /* a real transformation is necessary. The variable x is either shifted by -lb or
    * replaced by ub - x, depending on the smaller absolute of lb and ub.
    */
   if( !SCIPisFeasZero(scip, deltashift) || negatecoeffs )
   {
      SCIP_Real* vals;
      int* rows;
      int nrows;
      int i;

      assert(!SCIPisInfinity(scip, deltashift));

      /* get nonzero values and corresponding rows of column */
      getColumnData(matrix, colpos, &vals, &rows, &nrows);
      assert(nrows == 0 ||(vals != NULL && rows != NULL));

      /* go through rows and modify its lhs, rhs and the variable coefficient, if necessary */
      for( i = 0; i < nrows; ++i )
      {
         int rowpos = rows[i];
         assert(rowpos >= 0);
         assert(rowpos < matrix->nrows);

         if( !SCIPisInfinity(scip, -(matrix->lhs[rowpos])) )
            matrix->lhs[rowpos] -= (vals[i]) * deltashift;

         if( !SCIPisInfinity(scip, matrix->rhs[rowpos]) )
            matrix->rhs[rowpos] -= (vals[i]) * deltashift;

         if( negatecoeffs )
           (vals[i]) = -(vals[i]);

         assert(SCIPisFeasLE(scip, matrix->lhs[rowpos], matrix->rhs[rowpos]));
      }
   }
   SCIPdebugMsg(scip, "Variable <%s> at colpos %d transformed. LB <%g> --> <%g>, UB <%g> --> <%g>\n",
      SCIPvarGetName(var), colpos, lb, 0.0, ub, matrix->upperbounds[colpos]);
}

/** initializes copy of the original coefficient matrix and applies heuristic specific adjustments: normalizing row
 *  vectors, transforming variable domains such that lower bound is zero, and relaxing continuous variables.
 */
static
SCIP_RETCODE initMatrix(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object to be initialized */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  colposs,            /**< position of columns according to variable type sorting */
   SCIP_Bool             normalize,          /**< should coefficients and be normalized by rows maximum norms? */
   int*                  nmaxrows,           /**< maximum number of rows a variable appears in */
   SCIP_Bool             relax,              /**< should continuous variables be relaxed from the problem? */
   SCIP_Bool*            initialized,        /**< was the initialization successful? */
   SCIP_Bool*            infeasible          /**< is the problem infeasible? */
   )
{
   SCIP_ROW** lprows;
   SCIP_COL** lpcols;
   SCIP_Bool impliscontinuous;
   int i;
   int j;
   int currentpointer;

   int nrows;
   int ncols;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(initialized!= NULL);
   assert(infeasible != NULL);
   assert(nmaxrows != NULL);

   SCIPdebugMsg(scip, "entering Matrix Initialization method of SHIFTANDPROPAGATE heuristic!\n");

   /* get LP row data; column data is already initialized in heurdata */
   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nrows) );
   lpcols = heurdata->lpcols;
   ncols = heurdata->nlpcols;

   matrix->nrows = nrows;
   matrix->nnonzs = 0;
   matrix->normalized = FALSE;
   matrix->ndiscvars = 0;
   *nmaxrows = 0;
   impliscontinuous = heurdata->impliscontinuous;

   /* count the number of nonzeros of the LP constraint matrix */
   for( j = 0; j < ncols; ++j )
   {
      assert(lpcols[j] != NULL);
      assert(SCIPcolGetLPPos(lpcols[j]) >= 0);

      if( colIsDiscrete(lpcols[j], impliscontinuous) )
      {
         matrix->nnonzs += SCIPcolGetNLPNonz(lpcols[j]);
         ++matrix->ndiscvars;
      }
   }

   matrix->ncols = matrix->ndiscvars;

   if( matrix->nnonzs == 0 )
   {
      SCIPdebugMsg(scip, "No matrix entries - Terminating initialization of matrix.\n");

      *initialized = FALSE;

      return SCIP_OKAY;
   }

   /* allocate memory for the members of heuristic matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatvals, matrix->nnonzs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatind, matrix->nnonzs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatvals, matrix->nnonzs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatind, matrix->nnonzs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatbegin, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatbegin, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lhs, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rhs, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colnorms, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->violrows, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->transformstatus, matrix->ndiscvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->upperbounds, matrix->ndiscvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->transformshiftvals, matrix->ndiscvars) );

   /* set transform status of variables */
   for( j = 0; j < matrix->ndiscvars; ++j )
      matrix->transformstatus[j] = TRANSFORMSTATUS_NONE;

   currentpointer = 0;
   *infeasible = FALSE;

   /* initialize the rows vector of the heuristic matrix together with its corresponding
    * lhs, rhs.
    */
   for( i = 0; i < nrows; ++i )
   {
      SCIP_COL** cols;
      SCIP_ROW* row;
      SCIP_Real* rowvals;
      SCIP_Real constant;
      SCIP_Real maxval;
      int nrowlpnonz;

      /* get LP row information */
      row = lprows[i];
      rowvals = SCIProwGetVals(row);
      nrowlpnonz = SCIProwGetNLPNonz(row);
      maxval = SCIPgetRowMaxCoef(scip, row);
      cols = SCIProwGetCols(row);
      constant = SCIProwGetConstant(row);

      SCIPdebugMsg(scip, " %s : lhs=%g, rhs=%g, maxval=%g \n", SCIProwGetName(row), matrix->lhs[i], matrix->rhs[i], maxval);
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
      assert(!SCIPisInfinity(scip, constant));

      matrix->rowmatbegin[i] = currentpointer;

      /* modify the lhs and rhs w.r.t to the rows constant and normalize by 1-norm, i.e divide the lhs and rhs by the
       * maximum absolute value of the row
       */
      if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
         matrix->lhs[i] = SCIProwGetLhs(row) - constant;
      else
         matrix->lhs[i] = -SCIPinfinity(scip);

      if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) )
         matrix->rhs[i] = SCIProwGetRhs(row) - constant;
      else
         matrix->rhs[i] = SCIPinfinity(scip);

      /* make sure that maxval is larger than zero before normalization.
       * Maxval may be zero if the constraint contains no variables but is modifiable, hence not redundant
       */
      if( normalize && !SCIPisFeasZero(scip, maxval) )
      {
         if( !SCIPisInfinity(scip, -matrix->lhs[i]) )
            matrix->lhs[i] /= maxval;
         if( !SCIPisInfinity(scip, matrix->rhs[i]) )
            matrix->rhs[i] /= maxval;
      }


      /* in case of empty rows with a 0 < lhs <= 0.0 or 0.0 <= rhs < 0 we deduce the infeasibility of the problem */
      if( nrowlpnonz == 0 && (SCIPisFeasPositive(scip, matrix->lhs[i]) || SCIPisFeasNegative(scip, matrix->rhs[i])) )
      {
         *infeasible = TRUE;
         SCIPdebugMsg(scip, "  Matrix initialization stopped because of row infeasibility! \n");
         break;
      }

      /* row coefficients are normalized and copied to heuristic matrix */
      for( j = 0; j < nrowlpnonz; ++j )
      {
         if( !colIsDiscrete(cols[j], impliscontinuous) )
            continue;
         assert(SCIPcolGetLPPos(cols[j]) >= 0);
         assert(currentpointer < matrix->nnonzs);

         matrix->rowmatvals[currentpointer] = rowvals[j];
         if( normalize && SCIPisFeasGT(scip, maxval, 0.0) )
            matrix->rowmatvals[currentpointer] /= maxval;

         matrix->rowmatind[currentpointer] = colposs[SCIPcolGetLPPos(cols[j])];

         ++currentpointer;
      }
   }

   matrix->normalized = TRUE;

   if( *infeasible )
      return SCIP_OKAY;

   assert(currentpointer == matrix->nnonzs);

   currentpointer = 0;

   /* copy the nonzero coefficient data column by column to heuristic matrix */
   for( j = 0; j < matrix->ncols; ++j )
   {
      SCIP_COL* currentcol;
      SCIP_ROW** rows;
      SCIP_Real* colvals;
      int ncolnonz;


      assert(SCIPcolGetLPPos(lpcols[j]) >= 0);

      currentcol = lpcols[j];
      assert(colIsDiscrete(currentcol, impliscontinuous));

      colvals = SCIPcolGetVals(currentcol);
      rows = SCIPcolGetRows(currentcol);
      ncolnonz = SCIPcolGetNLPNonz(currentcol);
      matrix->colnorms[j] = ncolnonz;

      *nmaxrows = MAX(*nmaxrows, ncolnonz);

      /* loop over all rows with nonzero coefficients in the column, transform them and add them to the heuristic matrix */
      matrix->colmatbegin[j] = currentpointer;

      for( i = 0; i < ncolnonz; ++i )
      {
         SCIP_Real maxval;

         assert(rows[i] != NULL);
         assert(0 <= SCIProwGetLPPos(rows[i]));
         assert(SCIProwGetLPPos(rows[i]) < nrows);
         assert(currentpointer < matrix->nnonzs);

         /* rows are normalized by maximum norm */
         maxval = SCIPgetRowMaxCoef(scip, rows[i]);

         assert(maxval > 0);

         matrix->colmatvals[currentpointer] = colvals[i];
         if( normalize && SCIPisFeasGT(scip, maxval, 0.0) )
            matrix->colmatvals[currentpointer] /= maxval;

         matrix->colmatind[currentpointer] = SCIProwGetLPPos(rows[i]);

         /* update the column norm */
         matrix->colnorms[j] += ABS(matrix->colmatvals[currentpointer]);
         ++currentpointer;
      }
   }
   assert(currentpointer == matrix->nnonzs);

   /* each variable is either transformed, if it supposed to be integral, or relaxed */
   for( j = 0; j < (relax ? ncols : matrix->ndiscvars); ++j )
   {
      SCIP_COL* col;

      col = lpcols[j];
      if( colIsDiscrete(col, impliscontinuous) )
      {
         matrix->transformshiftvals[j] = 0.0;
         transformVariable(scip, matrix, heurdata, j);
      }
      else
      {
         SCIP_VAR* var;
         var = SCIPcolGetVar(col);
         assert(!varIsDiscrete(var, impliscontinuous));
         relaxVar(scip, var, matrix, normalize);
      }
   }
   *initialized = TRUE;

   SCIPdebugMsg(scip, "Matrix initialized for %d discrete variables with %d cols, %d rows and %d nonzero entries\n",
      matrix->ndiscvars, matrix->ncols, matrix->nrows, matrix->nnonzs);
   return SCIP_OKAY;
}

/** frees all members of the heuristic matrix */
static
void freeMatrix(
   SCIP*                 scip,               /**< current SCIP instance */
   CONSTRAINTMATRIX**    matrix              /**< constraint matrix object */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);

   /* all fields are only allocated, if problem is not empty  */
   if( (*matrix)->nnonzs > 0 )
   {
      assert((*matrix) != NULL);
      assert((*matrix)->rowmatbegin != NULL);
      assert((*matrix)->rowmatvals != NULL);
      assert((*matrix)->rowmatind != NULL);
      assert((*matrix)->colmatbegin != NULL);
      assert((*matrix)->colmatvals!= NULL);
      assert((*matrix)->colmatind != NULL);
      assert((*matrix)->lhs != NULL);
      assert((*matrix)->rhs != NULL);
      assert((*matrix)->transformstatus != NULL);
      assert((*matrix)->transformshiftvals != NULL);

      /* free all fields */
      SCIPfreeBufferArray(scip, &((*matrix)->transformshiftvals));
      SCIPfreeBufferArray(scip, &((*matrix)->upperbounds));
      SCIPfreeBufferArray(scip, &((*matrix)->transformstatus));
      SCIPfreeBufferArray(scip, &((*matrix)->violrows));
      SCIPfreeBufferArray(scip, &((*matrix)->colnorms));
      SCIPfreeBufferArray(scip, &((*matrix)->rhs));
      SCIPfreeBufferArray(scip, &((*matrix)->lhs));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatbegin));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatvals));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatvals));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatbegin));

     (*matrix)->nrows = 0;
     (*matrix)->ncols = 0;
   }

   /* free matrix */
   SCIPfreeBuffer(scip, matrix);
}

/** updates the information about a row whenever violation status changes */
static
void checkRowViolation(
   SCIP*                 scip,               /**< current SCIP instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   rowindex,           /**< index of the row */
   int*                  violatedrows,       /**< contains all violated rows */
   int*                  violatedrowpos,     /**< positions of rows in the violatedrows array */
   int*                  nviolatedrows,      /**< pointer to update total number of violated rows */
   int*                  rowweights,         /**< row weight storage */
   SCIP_Bool             updateweights       /**< should row weight be increased every time the row is violated? */
   )
{
   int* cols;
   int ncols;
   int c;
   int violadd;
   assert(matrix != NULL);
   assert(violatedrows != NULL);
   assert(violatedrowpos != NULL);
   assert(nviolatedrows != NULL);

   getRowData(matrix, rowindex, NULL, NULL, NULL, &cols, &ncols);
   violadd = 0;

   /* row is now violated. Enqueue it in the set of violated rows. */
   if( violatedrowpos[rowindex] == -1 && (SCIPisFeasGT(scip, matrix->lhs[rowindex], 0.0) || SCIPisFeasLT(scip, matrix->rhs[rowindex], 0.0)) )
   {
      assert(*nviolatedrows < matrix->nrows);

      violatedrows[*nviolatedrows] = rowindex;
      violatedrowpos[rowindex] = *nviolatedrows;
      ++(*nviolatedrows);
      if( updateweights )
         ++rowweights[rowindex];

      violadd = 1;
   }
   /* row is now feasible. Remove it from the set of violated rows. */
   else if( violatedrowpos[rowindex] >= 0 && SCIPisFeasLE(scip, matrix->lhs[rowindex], 0.0) && SCIPisFeasGE(scip, matrix->rhs[rowindex], 0.0) )
   {
      /* swap the row with last violated row */
      if( violatedrowpos[rowindex] != *nviolatedrows - 1  )
      {
         assert(*nviolatedrows - 1 >= 0);
         violatedrows[violatedrowpos[rowindex]] = violatedrows[*nviolatedrows - 1];
         violatedrowpos[violatedrows[*nviolatedrows - 1]] = violatedrowpos[rowindex];
      }

      /* unlink the row from its position in the array and decrease number of violated rows */
      violatedrowpos[rowindex] = -1;
      --(*nviolatedrows);
      violadd = -1;
   }

   /* increase or decrease the column violation counter */
   for( c = 0; c < ncols; ++c )
   {
      matrix->violrows[cols[c]] += violadd;
      assert(matrix->violrows[cols[c]] >= 0);
   }
}

/** collects the necessary information about row violations for the zero-solution. That is,
 *  all solution values in heuristic transformation are zero.
 */
static
void checkViolations(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   colidx,             /**< column index for specific column, or -1 for all rows */
   int*                  violatedrows,       /**< violated rows */
   int*                  violatedrowpos,     /**< row positions of violated rows */
   int*                  nviolatedrows,      /**< pointer to store the number of violated rows */
   int*                  rowweights,         /**< weight array for every row */
   SCIP_Bool             updateweights       /**< should row weight be increased every time the row is violated? */
   )
{
   int nrows;
   int* rowindices;
   int i;

   assert(matrix != NULL);
   assert(violatedrows != NULL);
   assert(violatedrowpos != NULL);
   assert(nviolatedrows != NULL);
   assert(-1 <= colidx && colidx < matrix->ncols);

   /* check if we requested an update for a single variable, or if we want to (re)-initialize the whole violation info */
   if( colidx >= 0 )
      getColumnData(matrix, colidx, NULL, &rowindices, &nrows);
   else
   {
      nrows = matrix->nrows;
      rowindices = NULL;
      *nviolatedrows = 0;

      /* reinitialize the violated rows */
      for( i = 0; i < nrows; ++i )
         violatedrowpos[i] = -1;

      /* clear the violated row counters for all variables */
      BMSclearMemoryArray(matrix->violrows, matrix->ndiscvars);
   }

   assert(colidx < 0 || *nviolatedrows >= 0);
   SCIPdebugMsg(scip, "Entering violation check for %d rows! \n", nrows);
   /* loop over rows and check if it is violated */
   for( i = 0; i < nrows; ++i )
   {
      int rowpos;
      if( colidx >= 0 )
      {
         assert(rowindices != NULL);
         rowpos = rowindices[i];
      }
      else
         rowpos = i;
      /* check, if zero solution violates this row */
      checkRowViolation(scip, matrix, rowpos, violatedrows, violatedrowpos, nviolatedrows, rowweights, updateweights);

      assert((violatedrowpos[rowpos] == -1 && SCIPisFeasGE(scip, matrix->rhs[rowpos], 0.0) && SCIPisFeasLE(scip, matrix->lhs[rowpos], 0.0))
         || (violatedrowpos[rowpos] >= 0 &&(SCIPisFeasLT(scip, matrix->rhs[rowpos], 0.0) || SCIPisFeasGT(scip, matrix->lhs[rowpos], 0.0))));
   }
}

/** retransforms solution values of variables according to their transformation status */
static
SCIP_Real retransformVariable(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   SCIP_VAR*             var,                /**< variable whose solution value has to be retransformed */
   int                   varindex,           /**< permutation of variable indices according to sorting */
   SCIP_Real             solvalue            /**< solution value of the variable */
   )
{
   TRANSFORMSTATUS status;

   assert(matrix != NULL);
   assert(var != NULL);

   status = matrix->transformstatus[varindex];
   assert(status != TRANSFORMSTATUS_NONE);

   /* check if original variable has different bounds and transform solution value correspondingly */
   if( status == TRANSFORMSTATUS_LB )
   {
      assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));

      return solvalue + matrix->transformshiftvals[varindex];
   }
   else if( status == TRANSFORMSTATUS_NEG )
   {
      assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
      return matrix->transformshiftvals[varindex] - solvalue;
   }
   return solvalue;
}

/** determines the best shifting value of a variable
 *  @todo if there is already an incumbent solution, try considering the objective cutoff as additional constraint */
static
SCIP_RETCODE getOptimalShiftingValue(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   varindex,           /**< index of variable which should be shifted */
   int                   direction,          /**< the direction for this variable */
   int*                  rowweights,         /**< weighting of rows for best shift calculation */
   SCIP_Real*            steps,              /**< buffer array to store the individual steps for individual rows */
   int*                  violationchange,    /**< buffer array to store the individual change of feasibility of row */
   SCIP_Real*            beststep,           /**< pointer to store optimal shifting step */
   int*                  rowviolations       /**< pointer to store new weighted sum of row violations, i.e, v - f */
   )
{
   SCIP_Real* vals;
   int* rows;

   SCIP_Real slacksurplus;
   SCIP_Real upperbound;

   int nrows;
   int sum;
   int i;

   SCIP_Bool allzero;

   assert(beststep != NULL);
   assert(rowviolations != NULL);
   assert(rowweights != NULL);
   assert(steps != NULL);
   assert(violationchange != NULL);
   assert(direction == 1 || direction == -1);

   upperbound = matrix->upperbounds[varindex];

   /* get nonzero values and corresponding rows of variable */
   getColumnData(matrix, varindex, &vals, &rows, &nrows);

   /* loop over rows and calculate, which is the minimum shift to make this row feasible
    * or the minimum shift to violate this row
    */
   allzero = TRUE;
   slacksurplus = 0.0;
   for( i = 0; i < nrows; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real val;
      int rowpos;
      SCIP_Bool rowisviolated;
      int rowweight;

      /* get the row data */
      rowpos = rows[i];
      assert(rowpos >= 0);
      lhs = matrix->lhs[rowpos];
      rhs = matrix->rhs[rowpos];
      rowweight = rowweights[rowpos];
      val = direction * vals[i];

      /* determine if current row is violated or not */
      rowisviolated =(SCIPisFeasLT(scip, rhs, 0.0) || SCIPisFeasLT(scip, -lhs, 0.0));

      /* for a feasible row, determine the minimum integer value within the bounds of the variable by which it has to be
       * shifted to make row infeasible.
       */
      if( !rowisviolated )
      {
         SCIP_Real maxfeasshift;

         maxfeasshift = SCIPinfinity(scip);

         /* feasibility can only be violated if the variable has a lock in the corresponding direction,
          * i.e. a positive coefficient for a "<="-constraint, a negative coefficient for a ">="-constraint.
          */
         if( SCIPisFeasGT(scip, val, 0.0) && !SCIPisInfinity(scip, rhs) )
            maxfeasshift = SCIPfeasFloor(scip, rhs/val);
         else if( SCIPisFeasLT(scip, val, 0.0) && !SCIPisInfinity(scip, -lhs) )
            maxfeasshift = SCIPfeasFloor(scip, lhs/val);

         /* if the variable has no lock in the current row, it can still help to increase the slack of this row;
          * we measure slack increase for shifting by one
          */
         if( SCIPisFeasGT(scip, val, 0.0) && SCIPisInfinity(scip, rhs) )
            slacksurplus += val;
         if( SCIPisFeasLT(scip, val, 0.0) && SCIPisInfinity(scip, -lhs) )
            slacksurplus -= val;

         /* check if the least violating shift lies within variable bounds and set corresponding array values */
         if( !SCIPisInfinity(scip, maxfeasshift) && SCIPisFeasLE(scip, maxfeasshift + 1.0, upperbound) )
         {
            steps[i] = maxfeasshift + 1.0;
            violationchange[i] = rowweight;
            allzero = FALSE;
         }
         else
         {
            steps[i] = upperbound;
            violationchange[i] = 0;
         }
      }
      /* for a violated row, determine the minimum integral value within the bounds of the variable by which it has to be
       * shifted to make row feasible.
       */
      else
      {
         SCIP_Real minfeasshift;

         minfeasshift = SCIPinfinity(scip);

         /* if coefficient has the right sign to make row feasible, determine the minimum integer to shift variable
          * to obtain feasibility
          */
         if( SCIPisFeasLT(scip, -lhs, 0.0) && SCIPisFeasGT(scip, val, 0.0) )
            minfeasshift = SCIPfeasCeil(scip, lhs/val);
         else if( SCIPisFeasLT(scip, rhs,0.0) && SCIPisFeasLT(scip, val, 0.0) )
            minfeasshift = SCIPfeasCeil(scip, rhs/val);

         /* check if the minimum feasibility recovery shift lies within variable bounds and set corresponding array
          * values
          */
         if( !SCIPisInfinity(scip, minfeasshift) && SCIPisFeasLE(scip, minfeasshift, upperbound) )
         {
            steps[i] = minfeasshift;
            violationchange[i] = -rowweight;
            allzero = FALSE;
         }
         else
         {
            steps[i] = upperbound;
            violationchange[i] = 0;
         }
      }
   }

   /* in case that the variable cannot affect the feasibility of any row, in particular it cannot violate
    * a single row, but we can add slack to already feasible rows, we will do this
    */
   if( allzero )
   {
      if( ! SCIPisInfinity(scip, upperbound) && SCIPisGT(scip, slacksurplus, 0.0) )
         *beststep = direction * upperbound;
      else
         *beststep = 0.0;

      return SCIP_OKAY;
   }

   /* sorts rows by increasing value of steps */
   SCIPsortRealInt(steps, violationchange, nrows);

   *beststep = 0.0;
   *rowviolations = 0;
   sum = 0;

   /* best shifting step is calculated by summing up the violation changes for each relevant step and
    * taking the one which leads to the minimum sum. This sum measures the balance of feasibility recovering and
    * violating changes which will be obtained by shifting the variable by this step
    * note, the sums for smaller steps have to be taken into account for all bigger steps, i.e., the sums can be
    * computed iteratively
    */
   for( i = 0; i < nrows && !SCIPisInfinity(scip, steps[i]); ++i )
   {
      sum += violationchange[i];

      /* if we reached the last entry for the current step value, we have finished computing its sum and
       * update the step defining the minimum sum
       */
      if( (i == nrows-1 || steps[i+1] > steps[i]) && sum < *rowviolations ) /*lint !e679*/
      {
         *rowviolations = sum;
         *beststep = direction * steps[i];
      }
   }
   assert(*rowviolations <= 0);
   assert(!SCIPisInfinity(scip, *beststep));

   return SCIP_OKAY;
}

/** updates transformation of a given variable by taking into account current local bounds. if the bounds have changed
 *  since last update, updating the heuristic specific upper bound of the variable, its current transformed solution value
 *  and all affected rows is necessary.
 */
static
SCIP_RETCODE updateTransformation(
   SCIP*                 scip,               /**< current scip */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   varindex,           /**< index of variable in matrix */
   SCIP_Real             lb,                 /**< local lower bound of the variable */
   SCIP_Real             ub,                 /**< local upper bound of the variable */
   int*                  violatedrows,       /**< violated rows */
   int*                  violatedrowpos,     /**< violated row positions */
   int*                  nviolatedrows       /**< pointer to store number of violated rows */
   )
{
   TRANSFORMSTATUS status;
   SCIP_Real deltashift;
   SCIP_Bool checkviolations;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= varindex && varindex < matrix->ndiscvars);

   /* deltashift is the difference between the old and new transformation value. */
   deltashift = 0.0;
   status = matrix->transformstatus[varindex];

   SCIPdebugMsg(scip, "  Variable <%d> [%g,%g], status %d(%g), ub %g \n", varindex, lb, ub, status,
      matrix->transformshiftvals[varindex], matrix->upperbounds[varindex]);

   checkviolations = FALSE;
   /* depending on the variable status, deltashift is calculated differently. */
   switch( status )
   {
   case TRANSFORMSTATUS_LB:
      if( SCIPisInfinity(scip, -lb) )
      {
         transformVariable(scip, matrix, heurdata, varindex);
         checkviolations = TRUE;
      }
      else
      {
         deltashift = lb - (matrix->transformshiftvals[varindex]);
         matrix->transformshiftvals[varindex] = lb;
         if( !SCIPisInfinity(scip, ub) )
            matrix->upperbounds[varindex] = ub - lb;
         else
            matrix->upperbounds[varindex] = SCIPinfinity(scip);
      }
      break;
   case TRANSFORMSTATUS_NEG:
      if( SCIPisInfinity(scip, ub) )
      {
         transformVariable(scip, matrix, heurdata, varindex);
         checkviolations = TRUE;
      }
      else
      {
         deltashift = (matrix->transformshiftvals[varindex]) - ub;
         matrix->transformshiftvals[varindex] = ub;

         if( !SCIPisInfinity(scip, -lb) )
            matrix->upperbounds[varindex] = ub - lb;
         else
            matrix->upperbounds[varindex] = SCIPinfinity(scip);
      }
      break;
   case TRANSFORMSTATUS_FREE:
      /* in case of a free transform status, if one of the bounds has become finite, we want
       * to transform this variable to a variable with a lowerbound or a negated transform status */
      if( !SCIPisInfinity(scip, -lb) || !SCIPisInfinity(scip, ub) )
      {
         transformVariable(scip, matrix, heurdata, varindex);

         /* violations have to be rechecked for rows in which variable appears */
         checkviolations = TRUE;

         assert(matrix->transformstatus[varindex] == TRANSFORMSTATUS_LB || TRANSFORMSTATUS_NEG);
         assert(SCIPisFeasLE(scip, ABS(lb), ABS(ub)) || matrix->transformstatus[varindex] == TRANSFORMSTATUS_NEG);
      }
      break;

   case TRANSFORMSTATUS_NONE:
   default:
      SCIPerrorMessage("Error: Invalid variable status <%d> in shift and propagagate heuristic, aborting!\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;   /*lint !e527*/
   }
   /* if the bound, by which the variable was shifted, has changed, deltashift is different from zero, which requires
    * an update of all affected rows
    */
   if( !SCIPisFeasZero(scip, deltashift) )
   {
      int i;
      int* rows;
      SCIP_Real* vals;
      int nrows;

      /* get nonzero values and corresponding rows of variable */
      getColumnData(matrix, varindex, &vals, &rows, &nrows);

      /* go through rows, update the rows w.r.t. the influence of the changed transformation of the variable */
      for( i = 0; i < nrows; ++i )
      {
         SCIPdebugMsg(scip, "  update slacks of row<%d>:  coefficient <%g>, %g <= 0 <= %g \n",
            rows[i], vals[i], matrix->lhs[rows[i]], matrix->rhs[rows[i]]);

         if( !SCIPisInfinity(scip, -(matrix->lhs[rows[i]])) )
            matrix->lhs[rows[i]] -= (vals[i]) * deltashift;

         if( !SCIPisInfinity(scip, matrix->rhs[rows[i]]) )
            matrix->rhs[rows[i]] -= (vals[i]) * deltashift;
      }
      checkviolations = TRUE;
   }

   /* check and update information about violated rows, if necessary */
   if( checkviolations )
      checkViolations(scip, matrix, varindex, violatedrows, violatedrowpos, nviolatedrows, heurdata->rowweights, heurdata->updateweights);

   SCIPdebugMsg(scip, "  Variable <%d> [%g,%g], status %d(%g), ub %g \n", varindex, lb, ub, status,
      matrix->transformshiftvals[varindex], matrix->upperbounds[varindex]);

   return SCIP_OKAY;
}

/** comparison method for columns; binary < integer < implicit < continuous variables */
static
SCIP_DECL_SORTPTRCOMP(heurSortColsShiftandpropagate)
{
   SCIP_COL* col1;
   SCIP_COL* col2;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VARTYPE vartype1;
   SCIP_VARTYPE vartype2;
   int key1;
   int key2;

   col1 = (SCIP_COL*)elem1;
   col2 = (SCIP_COL*)elem2;
   var1 = SCIPcolGetVar(col1);
   var2 = SCIPcolGetVar(col2);
   assert(var1 != NULL && var2 != NULL);

   vartype1 = SCIPvarGetType(var1);
   vartype2 = SCIPvarGetType(var2);

   switch (vartype1)
   {
      case SCIP_VARTYPE_BINARY:
         key1 = 1;
         break;
      case SCIP_VARTYPE_INTEGER:
         key1 = 2;
         break;
      case SCIP_VARTYPE_IMPLINT:
         key1 = 3;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         key1 = 4;
         break;
      default:
         key1 = -1;
         SCIPerrorMessage("unknown variable type\n");
         SCIPABORT();
         break;
   }
   switch (vartype2)
   {
      case SCIP_VARTYPE_BINARY:
         key2 = 1;
         break;
      case SCIP_VARTYPE_INTEGER:
         key2 = 2;
         break;
      case SCIP_VARTYPE_IMPLINT:
         key2 = 3;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         key2 = 4;
         break;
      default:
         key2 = -1;
         SCIPerrorMessage("unknown variable type\n");
         SCIPABORT();
         break;
   }
   return key1 - key2;
}

/*
 * Callback methods of primal heuristic
 */

/** deinitialization method of primal heuristic(called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitShiftandpropagate)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   /* if statistic mode is enabled, statistics are printed to console */
   SCIPstatistic(
      SCIPstatisticMessage(
         "  DETAILS                    :  %d violations left, %d probing status\n",
         heurdata->nremainingviols,
         heurdata->lpsolstat
         );
      SCIPstatisticMessage(
         "  SHIFTANDPROPAGATE PROBING  :  %d probings, %" SCIP_LONGINT_FORMAT " domain reductions,  ncutoffs: %d ,  LP iterations: %" SCIP_LONGINT_FORMAT " \n ",
         heurdata->nprobings,
         heurdata->ntotaldomredsfound,
         heurdata->ncutoffs,
         heurdata->nlpiters);
      );

   return SCIP_OKAY;
}

/** initialization method of primal heuristic(called after problem was transformed). We only need this method for
 *  statistic mode of heuristic.
 */
static
SCIP_DECL_HEURINIT(heurInitShiftandpropagate)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED) );

   SCIPstatistic(
      heurdata->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      heurdata->nremainingviols = 0;
      heurdata->nprobings = 0;
      heurdata->ntotaldomredsfound = 0;
      heurdata->ncutoffs = 0;
      heurdata->nlpiters = 0;
   )
   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data(called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeShiftandpropagate)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   eventhdlr = heurdata->eventhdlr;
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   SCIPfreeBlockMemoryNull(scip, &eventhdlrdata);

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);

   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** copy method for primal heuristic plugins(called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyShiftandpropagate)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurShiftandpropagate(scip) );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecShiftandpropagate)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;       /* heuristic data */
   SCIP_EVENTHDLR* eventhdlr;     /* shiftandpropagate event handler */
   SCIP_EVENTHDLRDATA* eventhdlrdata; /* event handler data */
   SCIP_EVENTDATA**  eventdatas;  /* event data for every variable */

   CONSTRAINTMATRIX* matrix;      /* constraint matrix object */
   SCIP_COL** lpcols;             /* lp columns */
   SCIP_SOL* sol;                 /* solution pointer */
   SCIP_Real* colnorms;           /* contains Euclidean norms of column vectors */

   SCIP_Real* steps;              /* buffer arrays for best shift selection in main loop */
   int* violationchange;

   int* violatedrows;             /* the violated rows */
   int* violatedrowpos;           /* the array position of a violated row, or -1 */
   int* permutation;              /* reflects the position of the variables after sorting */
   int* violatedvarrows;          /* number of violated rows for each variable */
   int* colposs;                  /* position of columns according to variable type sorting */
   int nlpcols;                   /* number of lp columns */
   int nviolatedrows;             /* number of violated rows */
   int ndiscvars;                 /* number of non-continuous variables of the problem */
   int lastindexofsusp;           /* last variable which has been swapped due to a cutoff */
   int nbinvars;                  /* number of binary variables */
   int nintvars;                  /* number of integer variables */
   int i;
   int r;
   int v;
   int c;
   int ncutoffs;                  /* counts the number of cutoffs for this execution */
   int nprobings;                 /* counts the number of probings */
   int nlprows;                   /* the number LP rows */
   int nmaxrows;                  /* maximum number of LP rows of a variable */

   SCIP_Bool initialized;         /* has the matrix been initialized? */
   SCIP_Bool cutoff;              /* has current probing node been cutoff? */
   SCIP_Bool probing;             /* should probing be applied or not? */
   SCIP_Bool infeasible;          /* FALSE as long as currently infeasible rows have variables left */
   SCIP_Bool impliscontinuous;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   eventhdlr = heurdata->eventhdlr;
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   *result = SCIP_DIDNOTRUN;
   SCIPdebugMsg(scip, "entering execution method of shift and propagate heuristic\n");

   /* heuristic is obsolete if there are only continuous variables */
   if( SCIPgetNVars(scip) - SCIPgetNContVars(scip) == 0 )
      return SCIP_OKAY;

   /* stop execution method if there is already a primarily feasible solution at hand */
   if( SCIPgetBestSol(scip) != NULL && heurdata->onlywithoutsol )
      return SCIP_OKAY;

   /* stop if there is no LP available */
   if ( ! SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   if( !SCIPisLPConstructed(scip) )
   {
      /* @note this call can have the side effect that variables are created */
      SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

      /* manually cut off the node if the LP construction detected infeasibility (heuristics cannot return such a result) */
      if( cutoff )
      {
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetCurrentNode(scip)) );
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPflushLP(scip) );
   }

   SCIPstatistic( heurdata->nlpiters = SCIPgetNLPIterations(scip) );

   nlprows = SCIPgetNLPRows(scip);

   SCIP_CALL( SCIPgetLPColsData(scip, &lpcols, &nlpcols) );
   assert(nlpcols == 0 || lpcols != NULL);

   /* we need an LP */
   if( nlprows == 0 || nlpcols == 0 )
      return SCIP_OKAY;


   *result = SCIP_DIDNOTFIND;
   initialized = FALSE;

   /* allocate lp column array */
   SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->lpcols, nlpcols) );
   heurdata->nlpcols = nlpcols;

   impliscontinuous = heurdata->impliscontinuous;

#ifndef NDEBUG
   BMSclearMemoryArray(heurdata->lpcols, nlpcols);
#endif

   /* copy and sort the columns by their variable types (binary before integer before implicit integer before continuous) */
   BMScopyMemoryArray(heurdata->lpcols, lpcols, nlpcols);

   SCIPsortPtr((void**)heurdata->lpcols, heurSortColsShiftandpropagate, nlpcols);

   SCIP_CALL( SCIPallocBufferArray(scip, &colposs, nlpcols) );

   /* we have to collect the number of different variable types before we start probing since during probing variable
    * can be created (e.g., cons_xor.c)
    */
   ndiscvars = 0;
   nbinvars = 0;
   nintvars = 0;
   for( c = 0; c < nlpcols; ++c )
   {
      SCIP_COL* col;
      SCIP_VAR* colvar;

      col = heurdata->lpcols[c];
      assert(col != NULL);
      colvar = SCIPcolGetVar(col);
      assert(colvar != NULL);

      if( varIsDiscrete(colvar, impliscontinuous) )
         ++ndiscvars;
      if( SCIPvarGetType(colvar) == SCIP_VARTYPE_BINARY )
         ++nbinvars;
      else if( SCIPvarGetType(colvar) == SCIP_VARTYPE_INTEGER )
         ++nintvars;

      /* save the position of this column in the array such that it can be accessed as the "true" column position */
      assert(SCIPcolGetLPPos(col) >= 0);
      colposs[SCIPcolGetLPPos(col)] = c;
   }
   assert(nbinvars + nintvars <= ndiscvars);

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   if( heurdata->collectstats )
      SCIPenableVarHistory(scip);
   else
      SCIPdisableVarHistory(scip);

   /* this should always be fulfilled becase we perform shift and propagate only at the root node */
   assert(SCIP_MAXTREEDEPTH > SCIPgetDepth(scip));

   /* @todo check if this node is necessary (I don't think so) */
   SCIP_CALL( SCIPnewProbingNode(scip) );
   ncutoffs = 0;
   nprobings = 0;
   nmaxrows = 0;
   infeasible = FALSE;

   /* initialize heuristic matrix and working solution */
   SCIP_CALL( SCIPallocBuffer(scip, &matrix) );
   SCIP_CALL( initMatrix(scip, matrix, heurdata, colposs, heurdata->normalize, &nmaxrows, heurdata->relax, &initialized, &infeasible) );

   /* the column positions are not needed anymore */
   SCIPfreeBufferArray(scip, &colposs);

   /* could not initialize matrix */
   if( !initialized || infeasible )
   {
      SCIPdebugMsg(scip, " MATRIX not initialized -> Execution of heuristic stopped! \n");
      goto TERMINATE;
   }

   /* the number of discrete LP column variables can be less than the actual number of variables, if, e.g., there
    * are nonlinearities in the problem. The heuristic execution can be terminated in that case.
    */
   if( matrix->ndiscvars < ndiscvars )
   {
      SCIPdebugMsg(scip, "Not all discrete variables are in the current LP. Shiftandpropagate execution terminated.\n");
      goto TERMINATE;
   }

   assert(nmaxrows > 0);

   eventhdlrdata->matrix = matrix;
   eventhdlrdata->heurdata = heurdata;

   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   SCIPsolSetHeur(sol, heur);

   /* allocate arrays for execution method */
   SCIP_CALL( SCIPallocBufferArray(scip, &permutation, ndiscvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->rowweights, matrix->nrows) );

   /* allocate necessary memory for best shift search */
   SCIP_CALL( SCIPallocBufferArray(scip, &steps, nmaxrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violationchange, nmaxrows) );

   /* allocate arrays to store information about infeasible rows */
   SCIP_CALL( SCIPallocBufferArray(scip, &violatedrows, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violatedrowpos, matrix->nrows) );

   eventhdlrdata->violatedrows = violatedrows;
   eventhdlrdata->violatedrowpos = violatedrowpos;
   eventhdlrdata->nviolatedrows = &nviolatedrows;



   /* initialize arrays. Before sorting, permutation is the identity permutation */
   for( i = 0; i < ndiscvars; ++i )
      permutation[i] = i;

   /* initialize row weights */
   for( r = 0; r < matrix->nrows; ++r )
   {
      if( !SCIPisInfinity(scip, -(matrix->lhs[r])) && !SCIPisInfinity(scip, matrix->rhs[r]) )
         heurdata->rowweights[r] = DEFAULT_WEIGHT_EQUALITY;
      else
         heurdata->rowweights[r] = DEFAULT_WEIGHT_INEQUALITY;

   }
   colnorms = matrix->colnorms;

   assert(nbinvars >= 0);
   assert(nintvars >= 0);

   /* check rows for infeasibility */
   checkViolations(scip, matrix, -1, violatedrows, violatedrowpos, &nviolatedrows, heurdata->rowweights, heurdata->updateweights);

   /* allocate memory for violatedvarrows array only if variable ordering relies on it */
   if( heurdata->sortvars && (heurdata->sortkey == 't' || heurdata->sortkey == 'v') )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &violatedvarrows, ndiscvars) );
      BMScopyMemoryArray(violatedvarrows, matrix->violrows, ndiscvars);
   }
   else
      violatedvarrows = NULL;

   /* sort variables w.r.t. the sorting key parameter. Sorting is indirect, all matrix column data
    * stays in place, but permutation array gives access to the sorted order of variables
    */
   if( heurdata->sortvars )
   {
      switch (heurdata->sortkey)
      {
         case 'n':
            /* variable ordering w.r.t. column norms nonincreasing */
            if( heurdata->preferbinaries )
            {
               if( nbinvars > 0 )
                  SCIPsortDownRealInt(colnorms, permutation, nbinvars);
               if( nbinvars < ndiscvars )
                  SCIPsortDownRealInt(&colnorms[nbinvars], &permutation[nbinvars], ndiscvars - nbinvars);
            }
            else
            {
               SCIPsortDownRealInt(colnorms, permutation, ndiscvars);
            }
            SCIPdebugMsg(scip, "Variables sorted down w.r.t their normalized columns!\n");
            break;
         case 'u':
            /* variable ordering w.r.t. column norms nondecreasing */
            if( heurdata->preferbinaries )
            {
               if( nbinvars > 0 )
                  SCIPsortRealInt(colnorms, permutation, nbinvars);
               if( nbinvars < ndiscvars )
                  SCIPsortRealInt(&colnorms[nbinvars], &permutation[nbinvars], ndiscvars - nbinvars);
            }
            else
            {
               SCIPsortRealInt(colnorms, permutation, ndiscvars);
            }
            SCIPdebugMsg(scip, "Variables sorted w.r.t their normalized columns!\n");
            break;
         case 'v':
            /* variable ordering w.r.t. nonincreasing number of violated rows */
            assert(violatedvarrows != NULL);
            if( heurdata->preferbinaries )
            {
               if( nbinvars > 0 )
                  SCIPsortDownIntInt(violatedvarrows, permutation, nbinvars);
               if( nbinvars < ndiscvars )
                  SCIPsortDownIntInt(&violatedvarrows[nbinvars], &permutation[nbinvars], ndiscvars - nbinvars);
            }
            else
            {
               SCIPsortDownIntInt(violatedvarrows, permutation, ndiscvars);
            }

            SCIPdebugMsg(scip, "Variables sorted down w.r.t their number of currently infeasible rows!\n");
            break;
         case 't':
            /* variable ordering w.r.t. nondecreasing number of violated rows */
            assert(violatedvarrows != NULL);
            if( heurdata->preferbinaries )
            {
               if( nbinvars > 0 )
                  SCIPsortIntInt(violatedvarrows, permutation, nbinvars);
               if( nbinvars < ndiscvars )
                  SCIPsortIntInt(&violatedvarrows[nbinvars], &permutation[nbinvars], ndiscvars - nbinvars);
            }
            else
            {
               SCIPsortIntInt(violatedvarrows, permutation, ndiscvars);
            }

            SCIPdebugMsg(scip, "Variables sorted (upwards) w.r.t their number of currently infeasible rows!\n");
            break;
         case 'r':
            /* random sorting */
            if( heurdata->preferbinaries )
            {
               if( nbinvars > 0 )
                  SCIPrandomPermuteIntArray(heurdata->randnumgen, permutation, 0, nbinvars - 1);
               if( nbinvars < ndiscvars )
                  SCIPrandomPermuteIntArray(heurdata->randnumgen, &permutation[nbinvars], nbinvars - 1,
                        ndiscvars - nbinvars - 1);
            }
            else
            {
               SCIPrandomPermuteIntArray(heurdata->randnumgen, permutation, 0, ndiscvars - 1);
            }
            SCIPdebugMsg(scip, "Variables permuted randomly!\n");
            break;
         default:
            SCIPdebugMsg(scip, "No variable permutation applied\n");
            break;
      }
   }

   /* should binary variables without locks be treated first? */
   if( heurdata->binlocksfirst )
   {
      SCIP_VAR* var;
      int nbinwithoutlocks = 0;

      /* count number of binaries without locks */
      if( heurdata->preferbinaries )
      {
         for( c = 0; c < nbinvars; ++c )
         {
            var = SCIPcolGetVar(heurdata->lpcols[permutation[c]]);
            if( SCIPvarGetNLocksUp(var) == 0 || SCIPvarGetNLocksDown(var) == 0 )
               ++nbinwithoutlocks;
         }
      }
      else
      {
         for( c = 0; c < ndiscvars; ++c )
         {
            var = SCIPcolGetVar(heurdata->lpcols[permutation[c]]);
            if( SCIPvarIsBinary(var) )
            {
               if( SCIPvarGetNLocksUp(var) == 0 || SCIPvarGetNLocksDown(var) == 0 )
                  ++nbinwithoutlocks;
            }
         }
      }

      if( nbinwithoutlocks > 0 )
      {
         SCIP_VAR* binvar;
         int b = 1;
         int tmp;
         c = 0;

         /* if c reaches nbinwithoutlocks, then all binary variables without locks were sorted to the beginning of the array */
         while( c < nbinwithoutlocks && b < ndiscvars )
         {
            assert(c < b);
            assert(c < ndiscvars);
            assert(b < ndiscvars);
            var = SCIPcolGetVar(heurdata->lpcols[permutation[c]]);
            binvar = SCIPcolGetVar(heurdata->lpcols[permutation[b]]);

            /* search for next variable which is not a binary variable without locks */
            while( SCIPvarIsBinary(var) && (SCIPvarGetNLocksUp(var) == 0 || SCIPvarGetNLocksDown(var) == 0) )
            {
               ++c;
               if( c >= nbinwithoutlocks )
                  break;
               var = SCIPcolGetVar(heurdata->lpcols[permutation[c]]);
            }
            if( c >= nbinwithoutlocks )
               break;

            /* search for next binary variable without locks (with position > c) */
            if( b <= c )
            {
               b = c + 1;
               binvar = SCIPcolGetVar(heurdata->lpcols[permutation[b]]);
            }
            while( !SCIPvarIsBinary(binvar) || (SCIPvarGetNLocksUp(binvar) > 0 && SCIPvarGetNLocksDown(binvar) > 0) )
            {
               ++b;
               assert(b < ndiscvars);
               binvar = SCIPcolGetVar(heurdata->lpcols[permutation[b]]);
            }

            /* swap the two variables */
            tmp = permutation[b];
            permutation[b] = permutation[c];
            permutation[c] = tmp;

            /* increase counters */
            ++c;
            ++b;
         }
      }

#ifndef NDEBUG
      for( c = 0; c < ndiscvars; ++c )
      {
         assert((c < nbinwithoutlocks) == (SCIPvarIsBinary(SCIPcolGetVar(heurdata->lpcols[permutation[c]]))
               && (SCIPvarGetNLocksUp(SCIPcolGetVar(heurdata->lpcols[permutation[c]])) == 0
                  || SCIPvarGetNLocksDown(SCIPcolGetVar(heurdata->lpcols[permutation[c]])) == 0)));
      }
#endif
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &eventdatas, matrix->ndiscvars) );
   BMSclearMemoryArray(eventdatas, matrix->ndiscvars);

   /* initialize variable events to catch bound changes during propagation */
   for( c = 0; c < matrix->ndiscvars; ++c )
   {
      SCIP_VAR* var;

      var = SCIPcolGetVar(heurdata->lpcols[c]);
      assert(var != NULL);
      assert(SCIPvarIsIntegral(var));
      assert(eventdatas[c] == NULL);

      SCIP_CALL( SCIPallocBuffer(scip, &(eventdatas[c])) ); /*lint !e866*/

      eventdatas[c]->colpos = c;

      SCIP_CALL( SCIPcatchVarEvent(scip, var, EVENTTYPE_SHIFTANDPROPAGATE, eventhdlr, eventdatas[c], NULL) );
   }

   cutoff = FALSE;

   lastindexofsusp = -1;
   probing = heurdata->probing;
   infeasible = FALSE;

   SCIPdebugMsg(scip, "SHIFT_AND_PROPAGATE heuristic starts main loop with %d violations and %d remaining variables!\n",
      nviolatedrows, ndiscvars);

   assert(matrix->ndiscvars == ndiscvars);

   /* loop over variables, shift them according to shifting criteria and try to reduce the global infeasibility */
   for( c = 0; c < ndiscvars; ++c )
   {
      SCIP_VAR* var;
      SCIP_Longint ndomredsfound;
      SCIP_Real optimalshiftvalue;
      SCIP_Real origsolval;
      SCIP_Real lb;
      SCIP_Real ub;
      int nviolations;
      int permutedvarindex;
      int j;
      SCIP_Bool marksuspicious;

      if( heurdata->selectbest )
      {  /* search for best candidate */
         j = c + 1;
         while( j < ndiscvars )
         {
            /* run through remaining variables and search for best candidate */
            if( matrix->violrows[permutation[c]] < matrix->violrows[permutation[j]] )
            {
               int tmp;
               tmp = permutation[c];
               permutation[c] = permutation[j];
               permutation[j] = tmp;
            }
            ++j;
         }
      }
      permutedvarindex = permutation[c];
      optimalshiftvalue = 0.0;
      nviolations = 0;
      var = SCIPcolGetVar(heurdata->lpcols[permutedvarindex]);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      assert(SCIPcolGetLPPos(SCIPvarGetCol(var)) >= 0);
      assert(SCIPvarIsIntegral(var));

      /* check whether we hit some limit, e.g. the time limit, in between
       * since the check itself consumes some time, we only do it every tenth iteration
       */
      if( c % 10 == 0 && SCIPisStopped(scip) )
         goto TERMINATE2;

      /* if propagation is enabled, check if propagation has changed the variables bounds
       * and update the transformed upper bound correspondingly
       * @todo this should not be necessary
       */
      if( heurdata->probing )
         SCIP_CALL( updateTransformation(scip, matrix, heurdata, permutedvarindex,lb, ub, violatedrows, violatedrowpos,
               &nviolatedrows) );

      SCIPdebugMsg(scip, "Variable %s with local bounds [%g,%g], status <%d>, matrix bound <%g>\n",
         SCIPvarGetName(var), lb, ub, matrix->transformstatus[permutedvarindex], matrix->upperbounds[permutedvarindex]);

      /* ignore variable if propagation fixed it (lb and ub will be zero) */
      if( SCIPisFeasZero(scip, matrix->upperbounds[permutedvarindex]) )
      {
         assert(!SCIPisInfinity(scip, ub));
         assert(SCIPisFeasEQ(scip, lb, ub));

         SCIP_CALL( SCIPsetSolVal(scip, sol, var, ub) );

         continue;
      }

      marksuspicious = FALSE;

      /* check whether the variable is binary and has no locks in one direction, so that we want to fix it to the
       * respective bound (only enabled by parameter)
       */
      if( heurdata->fixbinlocks && SCIPvarIsBinary(var) && (SCIPvarGetNLocksUp(var) == 0 || SCIPvarGetNLocksDown(var) == 0) )
      {
         if( SCIPvarGetNLocksUp(var) == 0 )
            origsolval = SCIPvarGetUbLocal(var);
         else
         {
            assert(SCIPvarGetNLocksDown(var) == 0);
            origsolval = SCIPvarGetLbLocal(var);
         }
      }
      else
      {
         /* only apply the computationally expensive best shift selection, if there is a violated row left */
         if( !heurdata->stopafterfeasible || nviolatedrows > 0 )
         {
            /* compute optimal shift value for variable */
            SCIP_CALL( getOptimalShiftingValue(scip, matrix, permutedvarindex, 1, heurdata->rowweights, steps, violationchange,
                  &optimalshiftvalue, &nviolations) );
            assert(SCIPisFeasGE(scip, optimalshiftvalue, 0.0));

            /* Variables with FREE transform have to be dealt with twice */
            if( matrix->transformstatus[permutedvarindex] == TRANSFORMSTATUS_FREE )
            {
               SCIP_Real downshiftvalue;
               int ndownviolations;

               downshiftvalue = 0.0;
               ndownviolations = 0;
               SCIP_CALL( getOptimalShiftingValue(scip, matrix, permutedvarindex, -1, heurdata->rowweights, steps, violationchange,
                     &downshiftvalue, &ndownviolations) );

               assert(SCIPisLE(scip, downshiftvalue, 0.0));

               /* compare to positive direction and select the direction which makes more rows feasible */
               if( ndownviolations < nviolations )
               {
                  optimalshiftvalue = downshiftvalue;
               }
            }
         }
         else
            optimalshiftvalue = 0.0;

         /* if zero optimal shift values are forbidden by the user parameter, delay the variable by marking it suspicious */
         if( heurdata->nozerofixing && nviolations > 0 && SCIPisFeasZero(scip, optimalshiftvalue) )
            marksuspicious = TRUE;

         /* retransform the solution value from the heuristic transformation space */
         assert(varIsDiscrete(var, impliscontinuous));
         origsolval = retransformVariable(scip, matrix, var, permutedvarindex, optimalshiftvalue);
      }
      assert(SCIPisFeasGE(scip, origsolval, lb) && SCIPisFeasLE(scip, origsolval, ub));

      /* check if propagation should still be performed
       * @todo do we need the hard coded value? we could use SCIP_MAXTREEDEPTH
       */
      if( nprobings > DEFAULT_PROPBREAKER )
         probing = FALSE;

      /* if propagation is enabled, fix the variable to the new solution value and propagate the fixation
       * (to fix other variables and to find out early whether solution is already infeasible)
       */
      if( !marksuspicious && probing )
      {
         /* this assert should be always fulfilled because we run this heuristic at the root node only and do not
          * perform probing if nprobings is less than DEFAULT_PROPBREAKER (currently: 65000)
          */
         assert(SCIP_MAXTREEDEPTH > SCIPgetDepth(scip));

         SCIP_CALL( SCIPnewProbingNode(scip) );
         SCIP_CALL( SCIPfixVarProbing(scip, var, origsolval) );
         ndomredsfound = 0;

         SCIPdebugMsg(scip, "  Shift %g(%g originally) is optimal, propagate solution\n", optimalshiftvalue, origsolval);
         SCIP_CALL( SCIPpropagateProbing(scip, heurdata->nproprounds, &cutoff, &ndomredsfound) );

         ++nprobings;
         SCIPstatistic( heurdata->ntotaldomredsfound += ndomredsfound );
         SCIPdebugMsg(scip, "Propagation finished! <%" SCIP_LONGINT_FORMAT "> domain reductions %s, <%d> probing depth\n", ndomredsfound, cutoff ? "CUTOFF" : "",
            SCIPgetProbingDepth(scip));
      }
      assert(!cutoff || probing);

      /* propagation led to an empty domain, hence we backtrack and postpone the variable */
      if( cutoff )
      {
         assert(probing);

         ++ncutoffs;

         /* only continue heuristic if number of cutoffs occured so far is reasonably small */
         if( heurdata->cutoffbreaker >= 0 && ncutoffs >= ((heurdata->maxcutoffquot * SCIPgetProbingDepth(scip)) + heurdata->cutoffbreaker) )
            break;

         cutoff = FALSE;

         /* backtrack to the parent of the current node */
         assert(SCIPgetProbingDepth(scip) >= 1);
         SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );




         /* this assert should be always fulfilled because we run this heuristic at the root node only and do not
          * perform probing if nprobings is less than DEFAULT_PROPBREAKER (currently: 65000)
          */
         assert(SCIP_MAXTREEDEPTH > SCIPgetDepth(scip));

         /* if the variable upper and lower bound are equal to the solution value to which we tried to fix the variable,
          * we are trapped at an infeasible node and break; this can only happen due to an intermediate global bound change of the variable,
          * I guess
          */
         if( SCIPisFeasEQ(scip, SCIPvarGetUbLocal(var), origsolval) && SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), origsolval) )
         {
            cutoff = TRUE;
            break;
         }
         else if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), origsolval) )
         {
            /* if the variable were to be set to one of its bounds, repropagate by tightening this bound by 1.0
             * into the direction of the other bound, if possible */
            assert(SCIPisFeasGE(scip, SCIPvarGetUbLocal(var), origsolval + 1.0));

            ndomredsfound = 0;
            SCIP_CALL( SCIPnewProbingNode(scip) );
            SCIP_CALL( SCIPchgVarLbProbing(scip, var, origsolval + 1.0) );
            SCIP_CALL( SCIPpropagateProbing(scip, heurdata->nproprounds, &cutoff, &ndomredsfound) );

            SCIPstatistic( heurdata->ntotaldomredsfound += ndomredsfound );
         }
         else if( SCIPisFeasEQ(scip, SCIPvarGetUbLocal(var), origsolval) )
         {
            /* if the variable were to be set to one of its bounds, repropagate by tightening this bound by 1.0
             * into the direction of the other bound, if possible */
            assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(var), origsolval - 1.0));

            ndomredsfound = 0;

            SCIP_CALL( SCIPnewProbingNode(scip) );
            SCIP_CALL( SCIPchgVarUbProbing(scip, var, origsolval - 1.0) );
            SCIP_CALL( SCIPpropagateProbing(scip, heurdata->nproprounds, &cutoff, &ndomredsfound) );

            SCIPstatistic( heurdata->ntotaldomredsfound += ndomredsfound );

         }
         /* if the tightened bound again leads to a cutoff, both subproblems are proven infeasible and the heuristic
          * can be stopped */
         if( cutoff )
         {
            break;
         }
         else
         {
            /* since repropagation was successful, we indicate that this variable led to a cutoff in one direction */
            marksuspicious = TRUE;
         }
      }

      if( marksuspicious )
      {
         /* mark the variable as suspicious */
         assert(permutedvarindex == permutation[c]);

         ++lastindexofsusp;
         assert(lastindexofsusp >= 0 && lastindexofsusp <= c);

         permutation[c] = permutation[lastindexofsusp];
         permutation[lastindexofsusp] = permutedvarindex;

         SCIPdebugMsg(scip, "  Suspicious variable! Postponed from pos <%d> to position <%d>\n", c, lastindexofsusp);
      }
      else
      {
         SCIPdebugMsg(scip, "Variable <%d><%s> successfully shifted by value <%g>!\n", permutedvarindex,
            SCIPvarGetName(var), optimalshiftvalue);

         /* update solution */
         SCIP_CALL( SCIPsetSolVal(scip, sol, var, origsolval) );

         /* only to ensure that some assertions can be made later on */
         if( !probing )
         {
            SCIP_CALL( SCIPfixVarProbing(scip, var, origsolval) );
         }
      }
   }
   SCIPdebugMsg(scip, "Heuristic finished with %d remaining violations and %d remaining variables!\n",
      nviolatedrows, lastindexofsusp + 1);

   /* if constructed solution might be feasible, go through the queue of suspicious variables and set the solution
    * values
    */
   if( nviolatedrows == 0 && !cutoff )
   {
      SCIP_Bool stored;
      SCIP_Bool trysol;

      for( v = 0; v <= lastindexofsusp; ++v )
      {
         SCIP_VAR* var;
         SCIP_Real origsolval;
         int permutedvarindex;

         /* get the column position of the variable */
         permutedvarindex = permutation[v];
         var = SCIPcolGetVar(heurdata->lpcols[permutedvarindex]);
         assert(varIsDiscrete(var, impliscontinuous));

         /* update the transformation of the variable, since the bound might have changed after the last update. */
         if( heurdata->probing )
            SCIP_CALL( updateTransformation(scip, matrix, heurdata, permutedvarindex, SCIPvarGetLbLocal(var),
                  SCIPvarGetUbLocal(var), violatedrows, violatedrowpos, &nviolatedrows) );

         /* retransform the solution value from the heuristic transformed space, set the solution value accordingly */
         assert(varIsDiscrete(var, impliscontinuous));
         origsolval = retransformVariable(scip, matrix, var, permutedvarindex, 0.0);
         assert(SCIPisFeasGE(scip, origsolval, SCIPvarGetLbLocal(var))
            && SCIPisFeasLE(scip, origsolval, SCIPvarGetUbLocal(var)));
         SCIP_CALL( SCIPsetSolVal(scip, sol, var, origsolval) );
         SCIP_CALL( SCIPfixVarProbing(scip, var, origsolval) ); /* only to ensure that some assertions can be made later */

         SCIPdebugMsg(scip, "  Remaining variable <%s> set to <%g>; %d Violations\n", SCIPvarGetName(var), origsolval,
            nviolatedrows);
      }

      /* Fixing of remaining variables led to infeasibility */
      if( nviolatedrows > 0 )
         goto TERMINATE2;

      trysol = TRUE;

      /* if the constructed solution might still be extendable to a feasible solution, try this by
       * solving the remaining LP
       */
      if( nlpcols != matrix->ndiscvars )
      {
         /* case that remaining LP has to be solved */
         SCIP_Bool lperror;

#ifndef NDEBUG
         {
            SCIP_VAR** vars;

            vars = SCIPgetVars(scip);
            assert(vars != NULL);
            /* ensure that all discrete variables in the remaining LP are fixed */
            for( v = 0; v < ndiscvars; ++v )
            {
               if( SCIPvarIsInLP(vars[v]) )
                  assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v])));

            }
         }
#endif

         SCIPdebugMsg(scip, " -> old LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));

#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPwriteLP(scip, "shiftandpropagatelp.mps") );
#endif
         /* solve LP;
          * errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
          * hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
          */
#ifdef NDEBUG
         {
            SCIP_RETCODE retstat;
            retstat = SCIPsolveProbingLP(scip, -1, &lperror, NULL);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving LP in SHIFTANDPROPAGATE heuristic; LP solve terminated with code <%d>\n",
                     retstat);
            }
         }
#else
         SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
#endif

         SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
         SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, SCIPgetLPSolstat(scip));

         /* check if this is a feasible solution */
         if( !lperror && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            /* copy the current LP solution to the working solution */
            SCIP_CALL( SCIPlinkLPSol(scip, sol) );
         }
         else
            trysol = FALSE;

         SCIPstatistic( heurdata->lpsolstat = SCIPgetLPSolstat(scip) );
      }

      /* check solution for feasibility, and add it to solution store if possible.
       * None of integrality, feasibility of LP rows, variable bounds have to be checked, because they
       * are guaranteed by the heuristic at this stage.
       */
      if( trysol )
      {
         SCIP_Bool printreason;
         SCIP_Bool completely;
#ifdef SCIP_DEBUG
         printreason = TRUE;
#else
         printreason = FALSE;
#endif
#ifndef NDEBUG
         completely = TRUE; /*lint !e838*/
#else
         completely = FALSE;
#endif

         /* we once also checked the variable bounds which should not be necessary */
         SCIP_CALL( SCIPtrySol(scip, sol, printreason, completely, FALSE, FALSE, FALSE, &stored) );

         if( stored )
         {
            SCIPdebugMsg(scip, "found feasible shifted solution:\n");
            SCIPdebug( SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) ) );
            *result = SCIP_FOUNDSOL;

            SCIPstatisticMessage("  Shiftandpropagate solution value: %16.9g \n", SCIPgetSolOrigObj(scip, sol));
         }
      }
   }
   else
   {
      SCIPdebugMsg(scip, "Solution constructed by heuristic is already known to be infeasible\n");
   }

   SCIPstatistic( heurdata->nremainingviols = nviolatedrows; );

 TERMINATE2:
   /* free allocated memory in reverse order of allocation */
   for( c = matrix->ndiscvars - 1; c >= 0; --c )
   {
      SCIP_VAR* var;

      var = SCIPcolGetVar(heurdata->lpcols[c]);
      assert(var != NULL);
      assert(eventdatas[c] != NULL);

      SCIP_CALL( SCIPdropVarEvent(scip, var, EVENTTYPE_SHIFTANDPROPAGATE, eventhdlr, eventdatas[c], -1) );
      SCIPfreeBuffer(scip, &(eventdatas[c]));
   }
   SCIPfreeBufferArray(scip, &eventdatas);

   if( violatedvarrows != NULL )
   {
      assert(heurdata->sortkey == 'v' || heurdata->sortkey == 't');
      SCIPfreeBufferArray(scip, &violatedvarrows);
   }
   /* free all allocated memory */
   SCIPfreeBufferArray(scip, &violatedrowpos);
   SCIPfreeBufferArray(scip, &violatedrows);
   SCIPfreeBufferArray(scip, &violationchange);
   SCIPfreeBufferArray(scip, &steps);
   SCIPfreeBufferArray(scip, &heurdata->rowweights);
   SCIPfreeBufferArray(scip, &permutation);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   eventhdlrdata->nviolatedrows = NULL;
   eventhdlrdata->violatedrowpos = NULL;
   eventhdlrdata->violatedrows = NULL;

 TERMINATE:
   /* terminate probing mode and free the remaining memory */
   SCIPstatistic(
      heurdata->ncutoffs += ncutoffs;
      heurdata->nprobings += nprobings;
      heurdata->nlpiters = SCIPgetNLPIterations(scip) - heurdata->nlpiters;
      );

   SCIP_CALL( SCIPendProbing(scip) );
   SCIPfreeBufferArray(scip, &heurdata->lpcols);
   freeMatrix(scip, &matrix);
   eventhdlrdata->matrix = NULL;

   return SCIP_OKAY;
}

/** event handler execution method for the heuristic which catches all
 *  events in which a lower or upper bound were tightened */
static
SCIP_DECL_EVENTEXEC(eventExecShiftandpropagate)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_VAR* var;
   SCIP_COL* col;
   SCIP_Real lb;
   SCIP_Real ub;
   int colpos;
   CONSTRAINTMATRIX* matrix;
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(EVENTHDLR_NAME, SCIPeventhdlrGetName(eventhdlr)) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   matrix = eventhdlrdata->matrix;

   heurdata = eventhdlrdata->heurdata;
   assert(heurdata != NULL && heurdata->lpcols != NULL);

   colpos = eventdata->colpos;

   assert(0 <= colpos && colpos < matrix->ndiscvars);

   col = heurdata->lpcols[colpos];
   var = SCIPcolGetVar(col);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   SCIP_CALL( updateTransformation(scip, matrix, eventhdlrdata->heurdata, colpos, lb, ub, eventhdlrdata->violatedrows,
         eventhdlrdata->violatedrowpos, eventhdlrdata->nviolatedrows) );

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the shiftandpropagate primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurShiftandpropagate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;
   SCIP_EVENTHDLRDATA* eventhandlerdata;
   SCIP_EVENTHDLR* eventhdlr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhandlerdata) );
   eventhandlerdata->matrix = NULL;

   eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecShiftandpropagate, eventhandlerdata) );
   assert(eventhdlr != NULL);

   /* create Shiftandpropagate primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   heurdata->rowweights = NULL;
   heurdata->nlpcols = 0;
   heurdata->eventhdlr = eventhdlr;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecShiftandpropagate, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyShiftandpropagate) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeShiftandpropagate) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitShiftandpropagate) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitShiftandpropagate) );


   /* add shiftandpropagate primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nproprounds",
         "The number of propagation rounds used for each propagation",
         &heurdata->nproprounds, TRUE, DEFAULT_NPROPROUNDS, -1, 1000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/relax", "Should continuous variables be relaxed?",
         &heurdata->relax, TRUE, DEFAULT_RELAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/probing", "Should domains be reduced by probing?",
         &heurdata->probing, TRUE, DEFAULT_PROBING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/onlywithoutsol",
         "Should heuristic only be executed if no primal solution was found, yet?",
         &heurdata->onlywithoutsol, TRUE, DEFAULT_ONLYWITHOUTSOL, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/cutoffbreaker", "The number of cutoffs before heuristic stops",
         &heurdata->cutoffbreaker, TRUE, DEFAULT_CUTOFFBREAKER, -1, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/sortkey",
         "the key for variable sorting: (n)orms down, norms (u)p, (v)iolations down, viola(t)ions up, or (r)andom",
         &heurdata->sortkey, TRUE, DEFAULT_SORTKEY, SORTKEYS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/sortvars", "Should variables be sorted for the heuristic?",
         &heurdata->sortvars, TRUE, DEFAULT_SORTVARS, NULL, NULL));
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/collectstats", "should variable statistics be collected during probing?",
         &heurdata->collectstats, TRUE, DEFAULT_COLLECTSTATS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/stopafterfeasible",
         "Should the heuristic stop calculating optimal shift values when no more rows are violated?",
         &heurdata->stopafterfeasible, TRUE, DEFAULT_STOPAFTERFEASIBLE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/preferbinaries",
         "Should binary variables be shifted first?",
         &heurdata->preferbinaries, TRUE, DEFAULT_PREFERBINARIES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/nozerofixing",
         "should variables with a zero shifting value be delayed instead of being fixed?",
         &heurdata->nozerofixing, TRUE, DEFAULT_NOZEROFIXING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/fixbinlocks",
         "should binary variables with no locks in one direction be fixed to that direction?",
         &heurdata->fixbinlocks, TRUE, DEFAULT_FIXBINLOCKS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/binlocksfirst",
         "should binary variables with no locks be preferred in the ordering?",
         &heurdata->binlocksfirst, TRUE, DEFAULT_BINLOCKSFIRST, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/normalize",
         "should coefficients and left/right hand sides be normalized by max row coeff?",
         &heurdata->normalize, TRUE, DEFAULT_NORMALIZE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/updateweights",
         "should row weight be increased every time the row is violated?",
         &heurdata->updateweights, TRUE, DEFAULT_UPDATEWEIGHTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/impliscontinuous",
         "should implicit integer variables be treated as continuous variables?",
         &heurdata->impliscontinuous, TRUE, DEFAULT_IMPLISCONTINUOUS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/selectbest",
         "should the heuristic choose the best candidate in every round? (set to FALSE for static order)?",
         &heurdata->selectbest, TRUE, DEFAULT_SELECTBEST, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxcutoffquot",
         "maximum percentage of allowed cutoffs before stopping the heuristic",
         &heurdata->maxcutoffquot, TRUE, DEFAULT_MAXCUTOFFQUOT, 0.0, 2.0, NULL, NULL) );

   return SCIP_OKAY;
}
