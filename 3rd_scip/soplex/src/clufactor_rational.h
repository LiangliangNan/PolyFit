/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  clufactor_rational.h
 * @brief Implementation of sparse LU factorization with Rational precision.
 */
#ifndef _CLUFACTOR_RATIONAL_H_
#define _CLUFACTOR_RATIONAL_H_

#include "spxdefines.h"
#include "slinsolver_rational.h"
#include "timer.h"
#include "svector.h"
#include "rational.h"
#include "basevectors.h"

#define WITH_L_ROWS 1

namespace soplex
{
/**@brief   Implementation of sparse LU factorization with Rational precision.
 * @ingroup Algo
 *
 * This class implements a sparse LU factorization with either
 * FOREST-TOMLIN or ETA updates, using dynamic Markowitz pivoting.
 */
class CLUFactorRational
{
public:

   //----------------------------------------
   /**@name Public types */
   //@{
   /** Doubly linked ring structure for garbage collection of column or
    *  row file in working matrix.
    */
   struct Dring
   {
      Dring* next;
      Dring* prev;
      int    idx;
   };

   /// Pivot Ring
   class Pring
   {
   public:
      Pring* next;                ///<
      Pring* prev;                ///<
      int    idx;                 ///< index of pivot row
      int    pos;                 ///< position of pivot column in row
      int    mkwtz;               ///< markowitz number of pivot

      Pring()                     ///< constructor
         : next(0)
         , prev(0)
         , idx(0)
         , pos(0)
         , mkwtz(0)
      {}

   private:
      Pring(const Pring&);             ///< blocked copy constructor
      Pring& operator= (const Pring&); ///< blocked assignment operator
   };
   //@}

protected:

   //----------------------------------------
   /**@name Protected types */
   //@{
   /// Temporary data structures.
   class Temp
   {
   public:
      int*    s_mark;       ///< marker
      DVectorRational s_max;     ///< maximum absolute value per row (or -1)
      int*    s_cact;       ///< lengths of columns of active submatrix
      int     stage;        ///< stage of the structure
      Pring   pivots;       ///< ring of selected pivot rows
      Pring*  pivot_col;    ///< column index handlers for Real linked list
      Pring*  pivot_colNZ;  ///< lists for columns to number of nonzeros
      Pring*  pivot_row;    ///< row index handlers for Real linked list
      Pring*  pivot_rowNZ;  ///< lists for rows to number of nonzeros

      Temp();               ///< constructor
      ~Temp();              ///< destructor
      void init(int p_dim); ///< initialization
      void clear();         ///< clears the structure

   private:
      Temp( const Temp& );             ///< blocked copy constructor
      Temp& operator= ( const Temp& ); ///< blocked assignment operator
   };

   /// Data structures for saving the row and column permutations.
   struct Perm
   {
      int* orig;          ///< orig[p] original index from p
      int* perm;          ///< perm[i] permuted index from i
   };

   /// Data structures for saving the working matrix and U factor.
   struct U
   {
      ///
      struct Row
      {
         Dring list;         /*!< \brief Double linked ringlist of vector
                               indices in the order they appear
                               in the row file                      */
         Dring* elem;        ///< %Array of ring elements.
         int    used;        ///< used entries of arrays idx and val
         DVectorRational val;      ///< hold nonzero values
         int*   idx;         ///< array of length val.dim() to hold column indices of nonzeros in val
         int*   start;       ///< starting positions in val and idx
         int*   len;         ///< used nonzeros per row vectors
         int*   max;         /*!< \brief maximum available nonzeros per row:
                               start[i] + max[i] == start[elem[i].next->idx]
                               len[i] <= max[i].                    */
      } row;

      ///
      struct Col
      {
         Dring  list;        /*!< \brief Double linked ringlist of vector
                                indices in the order they appear
                                in the column file                  */
         Dring* elem;        ///< %Array of ring elements.
         int    size;        ///< size of array idx
         int    used;        ///< used entries of array idx
         int*   idx;         ///< hold row indices of nonzeros
         DVectorRational val;      /*!< \brief hold nonzero values: this is only initialized
                                in the end of the factorization with DEFAULT
                                updates.                            */
         int*   start;       ///< starting positions in val and idx
         int*   len;         ///< used nonzeros per column vector
         int*   max;         /*!< \brief maximum available nonzeros per colunn:
                               start[i] + max[i] == start[elem[i].next->idx]
                               len[i] <= max[i].                    */
      } col;
   };


   /// Data structures for saving the working matrix and L factor.
   struct L
   {
      DVectorRational val;  ///< values of L vectors
      int*   idx;           ///< array of size val.dim() storing indices of L vectors
      int    startSize;     ///< size of array start
      int    firstUpdate;   ///< number of first update L vector
      int    firstUnused;   ///< number of first unused L vector
      int*   start;         ///< starting positions in val and idx
      int*   row;           ///< column indices of L vectors
      int    updateType;    ///< type of updates to be used.

      /* The following arrays have length |firstUpdate|, since they keep
       * rows of the L-vectors occuring during the factorization (without
       * updates), only:
       */
      DVectorRational rval; ///< values of rows of L
      int*   ridx;          ///< indices of rows of L
      int*   rbeg;          ///< start of rows in rval and ridx
      int*   rorig;         ///< original row permutation
      int*   rperm;         ///< original row permutation
   };
   //@}

   //----------------------------------------
   /**@name Protected data */
   //@{
   SLinSolverRational::Status stat;   ///< Status indicator.

   int     thedim;            ///< dimension of factorized matrix
   int     nzCnt;             ///< number of nonzeros in U
   Rational   initMaxabs;     ///< maximum abs number in initail Matrix
   Rational   maxabs;         ///< maximum abs number in L and U

   Real    rowMemMult;        ///< factor of minimum Memory * number of nonzeros
   Real    colMemMult;        ///< factor of minimum Memory * number of nonzeros
   Real    lMemMult;          ///< factor of minimum Memory * number of nonzeros

   Perm    row;               ///< row permutation matrices
   Perm    col;               ///< column permutation matrices

   L       l;                 ///< L matrix
   DVectorRational diag;      ///< Array of pivot elements
   U       u;                 ///< U matrix

   Rational*   work;          ///< Working array: must always be left as 0!

   Timer*  factorTime;        ///< Time spent in factorizations
   int     factorCount;       ///< Number of factorizations
   Real    timeLimit;         ///< Time limit on factorization or solves
   //@}

private:

   //----------------------------------------
   /**@name Private data */
   //@{
   Temp    temp;              ///< Temporary storage
   //@}

   //----------------------------------------
   /**@name Solving
      These helper methods are used during the factorization process.
      The solve*-methods solve lower and upper triangular systems from
      the left or from the right, respectively  The methods with '2' in
      the end solve two systems at the same time.  The methods with
      "Eps" in the end consider elements smaller then the passed epsilon
      as zero.
   */
   //@{
   ///
   void solveUright(Rational* wrk, Rational* vec);
   ///
   int  solveUrightEps(Rational* vec, int* nonz, Rational* rhs);
   ///
   void solveUright2(Rational* work1, Rational* vec1, Rational* work2, Rational* vec2);
   ///
   int  solveUright2eps(Rational* work1, Rational* vec1, Rational* work2, Rational* vec2, int* nonz);
   ///
   void solveLright2(Rational* vec1, Rational* vec2);
   ///
   void solveUpdateRight(Rational* vec);
   ///
   void solveUpdateRight2(Rational* vec1, Rational* vec2);
   ///
   void solveUleft(Rational* work, Rational* vec);
   ///
   void solveUleft2(Rational* work1, Rational* vec1, Rational* work2, Rational* vec2);
   ///
   int solveLleft2forest(Rational* vec1, int* /* nonz */, Rational* vec2);
   ///
   void solveLleft2(Rational* vec1, int* /* nonz */, Rational* vec2);
   ///
   int solveLleftForest(Rational* vec, int* /* nonz */);
   ///
   void solveLleft(Rational* vec);
   ///
   int solveLleftEps(Rational* vec, int* nonz);
   ///
   void solveUpdateLeft(Rational* vec);
   ///
   void solveUpdateLeft2(Rational* vec1, Rational* vec2);

   ///
   int vSolveLright(Rational* vec, int* ridx, int rn);
   ///
   void vSolveLright2(Rational* vec, int* ridx, int* rnptr,
      Rational* vec2, int* ridx2, int* rn2ptr);
   ///
   void vSolveLright3(Rational* vec, int* ridx, int* rnptr,
      Rational* vec2, int* ridx2, int* rn2ptr,
      Rational* vec3, int* ridx3, int* rn3ptr);
   ///
   int vSolveUright(Rational* vec, int* vidx, Rational* rhs, int* ridx, int rn);
   ///
   void vSolveUrightNoNZ(Rational* vec, Rational* rhs, int* ridx, int rn);
   ///
   int vSolveUright2(Rational* vec, int* vidx, Rational* rhs, int* ridx, int rn,
      Rational* vec2, Rational* rhs2, int* ridx2, int rn2);
   ///
   int vSolveUpdateRight(Rational* vec, int* ridx, int n);
   ///
   void vSolveUpdateRightNoNZ(Rational* vec);
   ///
   int solveUleft(Rational* vec, int* vecidx, Rational* rhs, int* rhsidx, int rhsn);
   ///
   void solveUleftNoNZ(Rational* vec, Rational* rhs, int* rhsidx, int rhsn);
   ///
   int solveLleftForest(Rational* vec, int* nonz, int n);
   ///
   void solveLleftForestNoNZ(Rational* vec);
   ///
   int solveLleft(Rational* vec, int* nonz, int rn);
   ///
   void solveLleftNoNZ(Rational* vec);
   ///
   int solveUpdateLeft(Rational* vec, int* nonz, int n);

   ///
   void forestPackColumns();
   ///
   void forestMinColMem(int size);
   ///
   void forestReMaxCol(int col, int len);

   ///
   void initPerm();
   ///
   void initFactorMatrix(const SVectorRational** vec );
   ///
   void minLMem(int size);
   ///
   void setPivot(const int p_stage, const int p_col, const int p_row, const Rational& val);
   ///
   void colSingletons();
   ///
   void rowSingletons();

   ///
   void initFactorRings();
   ///
   void freeFactorRings();

   ///
   int setupColVals();
   ///
   void setupRowVals();

   ///
   void eliminateRowSingletons();
   ///
   void eliminateColSingletons();
   ///
   void selectPivots(const Rational& threshold);
   ///
   int updateRow(int r, int lv, int prow, int pcol, const Rational& pval);

   ///
   void eliminatePivot(int prow, int pos);
   ///
   void eliminateNucleus(const Rational& threshold);
   ///
   void minRowMem(int size);
   ///
   void minColMem(int size);
   ///
   void remaxCol(int p_col, int len);
   ///
   void packRows();
   ///
   void packColumns();
   ///
   void remaxRow(int p_row, int len);
   ///
   int makeLvec(int p_len, int p_row);
   ///
   bool timeLimitReached()
   {
      if( timeLimit >= 0.0 && factorTime->time() >= timeLimit )
      {
         stat = SLinSolverRational::TIME;
         return true;
      }
      else
         return false;
   }
   //@}

protected:

   //----------------------------------------
   /**@name Solver methods */
   //@{
   ///
   void solveLright(Rational* vec);
   ///
   int  solveRight4update(Rational* vec, int* nonz, Rational* rhs,
      Rational* forest, int* forestNum, int* forestIdx);
   ///
   void solveRight(Rational* vec, Rational* rhs);
   ///
   int  solveRight2update(Rational* vec1, Rational* vec2, Rational* rhs1,
      Rational* rhs2, int* nonz, Rational* forest, int* forestNum, int* forestIdx);
   ///
   void solveRight2(Rational* vec1, Rational* vec2, Rational* rhs1, Rational* rhs2);
   ///
   void solveLeft(Rational* vec, Rational* rhs);
   ///
   int solveLeftEps(Rational* vec, Rational* rhs, int* nonz);
   ///
   int solveLeft2(Rational* vec1, int* nonz, Rational* vec2, Rational* rhs1, Rational* rhs2);

   ///
   int vSolveRight4update(
      Rational* vec, int* idx,               /* result       */
      Rational* rhs, int* ridx, int rn,      /* rhs & Forest */
      Rational* forest, int* forestNum, int* forestIdx);
   ///
   int vSolveRight4update2(
      Rational* vec, int* idx,              /* result1 */
      Rational* rhs, int* ridx, int rn,     /* rhs1    */
      Rational* vec2,                       /* result2 */
      Rational* rhs2, int* ridx2, int rn2,  /* rhs2    */
      Rational* forest, int* forestNum, int* forestIdx);
   ///
   int vSolveRight4update3(
      Rational* vec, int* idx,              /* result1 */
      Rational* rhs, int* ridx, int rn,     /* rhs1    */
      Rational* vec2,                       /* result2 */
      Rational* rhs2, int* ridx2, int rn2,  /* rhs2    */
      Rational* vec3,                       /* result3 */
      Rational* rhs3, int* ridx3, int rn3,  /* rhs3    */
      Rational* forest, int* forestNum, int* forestIdx);
   ///
   void vSolveRightNoNZ(Rational* vec2,              /* result2 */
      Rational* rhs2, int* ridx2, int rn2);          /* rhs2    */
   ///
   int vSolveLeft(
      Rational* vec, int* idx,                      /* result */
      Rational* rhs, int* ridx, int rn);            /* rhs    */
   ///
   void vSolveLeftNoNZ(
      Rational* vec,                           /* result */
      Rational* rhs, int* ridx, int rn);       /* rhs    */
   ///
   int vSolveLeft2(
      Rational* vec, int* idx,                     /* result */
      Rational* rhs, int* ridx, int rn,            /* rhs    */
      Rational* vec2,                              /* result2 */
      Rational* rhs2, int* ridx2, int rn2);        /* rhs2    */
   ///
   int vSolveLeft3(Rational* vec, int* idx,                     /* result */
                   Rational* rhs, int* ridx, int rn,            /* rhs    */
                   Rational* vec2,                              /* result2 */
                   Rational* rhs2, int* ridx2, int rn2,         /* rhs2    */
                   Rational* vec3,                              /* result3 */
                   Rational* rhs3, int* ridx3, int rn3);        /* rhs3    */

   void forestUpdate(int col, Rational* work, int num, int *nonz);

   void update(int p_col, Rational* p_work, const int* p_idx, int num);
   void updateNoClear(int p_col, const Rational* p_work, const int* p_idx, int num);

   ///
   void factor(const SVectorRational** vec,   ///< Array of column vector pointers
               const Rational& threshold);           ///< pivoting threshold
   //@}

   //----------------------------------------
   /**@name Debugging */
   //@{
   ///
   void dump() const;

   ///
   bool isConsistent() const;
   //@}
};

} // namespace soplex
#endif // _CLUFACTOR_RATIONAL_H_
