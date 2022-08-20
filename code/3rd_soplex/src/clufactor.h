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

/**@file  clufactor.h
 * @brief Implementation of sparse LU factorization.
 */
#ifndef _CLUFACTOR_H_
#define _CLUFACTOR_H_

#include "spxdefines.h"
#include "slinsolver.h"
#include "timer.h"
#include "svector.h"

#define WITH_L_ROWS 1

namespace soplex
{
/**@brief   Implementation of sparse LU factorization.
 * @ingroup Algo
 * 
 * This class implements a sparse LU factorization with either
 * FOREST-TOMLIN or ETA updates, using dynamic Markowitz pivoting.
 */
class CLUFactor
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

      Pring() : next(0), prev(0)  ///< constructor
      {
         mkwtz = -1;
         idx = -1;
         pos = -1;
      }

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
      Real*   s_max;        ///< maximum absolute value per row (or -1) 
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
         int    size;        ///< size of arrays val and idx
         int    used;        ///< used entries of arrays idx and val
         Real*  val;         ///< hold nonzero values
         int*   idx;         ///< hold column indices of nonzeros 
         int*   start;       ///< starting positions in val and idx
         int*   len;         ///< used nonzeros per row vectors
         int*   max;         /*!< \brief maximum available nonzeros per row:
                               start[i] + max[i] == start[elem[i].next->idx] 
                               len[i] <= max[i].                    */
      } row;

      ///
      struct Col
      {
         Dring list;         /*!< \brief Double linked ringlist of vector
                                indices in the order they appear
                                in the column file                  */
         Dring *elem;        ///< %Array of ring elements.
         int size;           ///< size of array idx
         int used;           ///< used entries of array idx
         int *idx;           ///< hold row indices of nonzeros
         Real *val;          /*!< \brief hold nonzero values: this is only initialized
                                in the end of the factorization with DEFAULT
                                updates.                            */
         int *start;         ///< starting positions in val and idx
         int *len;           ///< used nonzeros per column vector
         int *max;           /*!< \brief maximum available nonzeros per colunn:
                               start[i] + max[i] == start[elem[i].next->idx] 
                               len[i] <= max[i].                    */
      } col;
   };


   /// Data structures for saving the working matrix and L factor.
   struct L
   {
      int  size;           ///< size of arrays val and idx
      Real *val;           ///< values of L vectors
      int  *idx;           ///< indices of L vectors
      int  startSize;      ///< size of array start
      int  firstUpdate;    ///< number of first update L vector
      int  firstUnused;    ///< number of first unused L vector
      int  *start;         ///< starting positions in val and idx
      int  *row;           ///< column indices of L vectors
      int  updateType;     ///< type of updates to be used.

      /* The following arrays have length |firstUpdate|, since they keep
       * rows of the L-vectors occuring during the factorization (without
       * updates), only:
       */
      Real *rval;          ///< values of rows of L
      int  *ridx;          ///< indices of rows of L
      int  *rbeg;          ///< start of rows in rval and ridx
      int  *rorig;         ///< original row permutation
      int  *rperm;         ///< original row permutation
   };
   //@}

   //----------------------------------------
   /**@name Protected data */
   //@{
   SLinSolver::Status stat;   ///< Status indicator.

   int     thedim;            ///< dimension of factorized matrix   
   int     nzCnt;             ///< number of nonzeros in U      
   Real    initMaxabs;        ///< maximum abs number in initail Matrix 
   Real    maxabs;            ///< maximum abs number in L and U        

   Real    rowMemMult;        ///< factor of minimum Memory * number of nonzeros 
   Real    colMemMult;        ///< factor of minimum Memory * number of nonzeros 
   Real    lMemMult;          ///< factor of minimum Memory * number of nonzeros 

   Perm    row;               ///< row permutation matrices 
   Perm    col;               ///< column permutation matrices 

   L       l;                 ///< L matrix 
   Real*   diag;              ///< Array of pivot elements          
   U       u;                 ///< U matrix 

   Real*   work;              ///< Working array: must always be left as 0! 

   Timer*  factorTime;        ///< Time spent in factorizations
   int     factorCount;       ///< Number of factorizations
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
   void solveUright(Real* wrk, Real* vec) const;
   ///
   int  solveUrightEps(Real* vec, int* nonz, Real eps, Real* rhs);
   ///
   void solveUright2(Real* work1, Real* vec1, Real* work2, Real* vec2);
   ///
   int  solveUright2eps(Real* work1, Real* vec1, Real* work2, Real* vec2, int* nonz, Real eps);
   ///
   void solveLright2(Real* vec1, Real* vec2);
   ///
   void solveUpdateRight(Real* vec);
   ///
   void solveUpdateRight2(Real* vec1, Real* vec2);
   ///
   void solveUleft(Real* work, Real* vec);
   ///
   void solveUleft2(Real* work1, Real* vec1, Real* work2, Real* vec2);
   ///
   int solveLleft2forest(Real* vec1, int* /* nonz */, Real* vec2, Real /* eps */);
   ///
   void solveLleft2(Real* vec1, int* /* nonz */, Real* vec2, Real /* eps */);
   ///
   int solveLleftForest(Real* vec, int* /* nonz */, Real /* eps */);
   ///
   void solveLleft(Real* vec) const;
   ///
   int solveLleftEps(Real* vec, int* nonz, Real eps);
   ///
   void solveUpdateLeft(Real* vec);
   ///
   void solveUpdateLeft2(Real* vec1, Real* vec2);

   void inline updateSolutionVectorLright(Real change, int j, Real& vec, int* idx, int& nnz);
   ///
   void vSolveLright(Real* vec, int* ridx, int& rn, Real eps);
   ///
   void vSolveLright2(Real* vec, int* ridx, int& rn, Real eps,
      Real* vec2, int* ridx2, int& rn2, Real eps2);
   ///
   void vSolveLright3(Real* vec, int* ridx, int& rn, Real eps,
      Real* vec2, int* ridx2, int& rn2, Real eps2,
      Real* vec3, int* ridx3, int& rn3, Real eps3);
   ///
   int vSolveUright(Real* vec, int* vidx, Real* rhs, int* ridx, int rn, Real eps);
   ///
   void vSolveUrightNoNZ(Real* vec, Real* rhs, int* ridx, int rn, Real eps);
   ///
   int vSolveUright2(Real* vec, int* vidx, Real* rhs, int* ridx, int rn, Real eps,
      Real* vec2, Real* rhs2, int* ridx2, int rn2, Real eps2);
   ///
   int vSolveUpdateRight(Real* vec, int* ridx, int n, Real eps);
   ///
   void vSolveUpdateRightNoNZ(Real* vec, Real /*eps*/);
   ///
   int solveUleft(Real eps, Real* vec, int* vecidx, Real* rhs, int* rhsidx, int rhsn);
   ///
   void solveUleftNoNZ(Real eps, Real* vec, Real* rhs, int* rhsidx, int rhsn);
   ///
   int solveLleftForest(Real eps, Real* vec, int* nonz, int n);
   ///
   void solveLleftForestNoNZ(Real* vec);
   ///
   int solveLleft(Real eps, Real* vec, int* nonz, int rn);
   ///
   void solveLleftNoNZ(Real* vec);
   ///
   int solveUpdateLeft(Real eps, Real* vec, int* nonz, int n);

   ///
   void forestPackColumns();
   ///
   void forestMinColMem(int size);
   ///
   void forestReMaxCol(int col, int len);

   ///
   void initPerm();
   ///
   void initFactorMatrix(const SVector** vec, const Real eps );
   ///
   void minLMem(int size);
   ///
   void setPivot(const int p_stage, const int p_col, const int p_row, const Real val);
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
   void selectPivots(Real threshold);
   ///
   int updateRow(int r, int lv, int prow, int pcol, Real pval, Real eps);

   ///
   void eliminatePivot(int prow, int pos, Real eps);
   ///
   void eliminateNucleus(const Real eps, const Real threshold);
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
   //@}

protected:

   //----------------------------------------
   /**@name Solver methods */
   //@{
   ///
   void solveLright(Real* vec);
   ///
   int  solveRight4update(Real* vec, int* nonz, Real eps, Real* rhs,
      Real* forest, int* forestNum, int* forestIdx);
   ///
   void solveRight(Real* vec, Real* rhs);
   ///
   int  solveRight2update(Real* vec1, Real* vec2, Real* rhs1,
      Real* rhs2, int* nonz, Real eps, Real* forest, int* forestNum, int* forestIdx);
   ///
   void solveRight2(Real* vec1, Real* vec2, Real* rhs1, Real* rhs2);
   ///
   void solveLeft(Real* vec, Real* rhs);
   ///
   int solveLeftEps(Real* vec, Real* rhs, int* nonz, Real eps);
   ///
   int solveLeft2(Real* vec1, int* nonz, Real* vec2, Real eps, Real* rhs1, Real* rhs2);

   ///
   int vSolveRight4update(Real eps, 
      Real* vec, int* idx,               /* result       */
      Real* rhs, int* ridx, int rn,      /* rhs & Forest */
      Real* forest, int* forestNum, int* forestIdx);
   ///
   int vSolveRight4update2(Real eps,
      Real* vec, int* idx,              /* result1 */
      Real* rhs, int* ridx, int rn,     /* rhs1    */
      Real* vec2, Real eps2,            /* result2 */
      Real* rhs2, int* ridx2, int rn2,  /* rhs2    */
      Real* forest, int* forestNum, int* forestIdx);
   /// sparse version of above method
   void vSolveRight4update2sparse(
      Real eps, Real* vec, int* idx,    /* result1 */
      Real* rhs, int* ridx, int& rn,    /* rhs1    */
      Real eps2, Real* vec2, int* idx2, /* result2 */
      Real* rhs2, int* ridx2, int& rn2, /* rhs2    */
      Real* forest, int* forestNum, int* forestIdx);
   ///
   int vSolveRight4update3(Real eps,
      Real* vec, int* idx,              /* result1 */
      Real* rhs, int* ridx, int rn,     /* rhs1    */
      Real* vec2, Real eps2,            /* result2 */
      Real* rhs2, int* ridx2, int rn2,  /* rhs2    */
      Real* vec3, Real eps3,            /* result3 */
      Real* rhs3, int* ridx3, int rn3,  /* rhs3    */
      Real* forest, int* forestNum, int* forestIdx);
   /// sparse version of above method
   void vSolveRight4update3sparse(
      Real eps, Real* vec, int* idx,    /* result1 */
      Real* rhs, int* ridx, int& rn,    /* rhs1    */
      Real eps2, Real* vec2, int* idx2, /* result2 */
      Real* rhs2, int* ridx2, int& rn2, /* rhs2    */
      Real eps3, Real* vec3, int* idx3, /* result3 */
      Real* rhs3, int* ridx3, int& rn3, /* rhs3    */
      Real* forest, int* forestNum, int* forestIdx);
   ///
   void vSolveRightNoNZ(Real* vec, Real eps,    /* result */
      Real* rhs, int* ridx, int rn);            /* rhs    */
   ///
   int vSolveLeft(Real eps,
      Real* vec, int* idx,                      /* result */
      Real* rhs, int* ridx, int rn);            /* rhs    */
   ///
   void vSolveLeftNoNZ(Real eps,
      Real* vec,                           /* result */
      Real* rhs, int* ridx, int rn);       /* rhs    */
   ///
   int vSolveLeft2(Real eps,
      Real* vec, int* idx,                     /* result */
      Real* rhs, int* ridx, int rn,            /* rhs    */
      Real* vec2,                              /* result2 */
      Real* rhs2, int* ridx2, int rn2);        /* rhs2    */
   /// sparse version of solving 2 systems of equations
   void vSolveLeft2sparse(Real eps,
                          Real* vec, int* idx,                     /* result */
                          Real* rhs, int* ridx, int& rn,           /* rhs    */
                          Real* vec2, int* idx2,                   /* result2 */
                          Real* rhs2, int* ridx2, int& rn2);       /* rhs2    */
   ///
   int vSolveLeft3(Real eps,
                   Real* vec, int* idx,                     /* result */
                   Real* rhs, int* ridx, int rn,            /* rhs    */
                   Real* vec2,                              /* result2 */
                   Real* rhs2, int* ridx2, int rn2,         /* rhs2    */
                   Real* vec3,                              /* result3 */
                   Real* rhs3, int* ridx3, int rn3);        /* rhs3    */
   /// sparse version of solving 3 systems of equations
   void vSolveLeft3sparse(Real eps,
                          Real* vec, int* idx,                     /* result */
                          Real* rhs, int* ridx, int& rn,           /* rhs    */
                          Real* vec2, int* idx2,                   /* result2 */
                          Real* rhs2, int* ridx2, int& rn2,        /* rhs2    */
                          Real* vec3, int* idx3,                   /* result2 */
                          Real* rhs3, int* ridx3, int& rn3);       /* rhs2    */

   void forestUpdate(int col, Real* work, int num, int *nonz);

   void update(int p_col, Real* p_work, const int* p_idx, int num);
   void updateNoClear(int p_col, const Real* p_work, const int* p_idx, int num);

   ///
   void factor(const SVector** vec,   ///< Array of column vector pointers  
               Real threshold,    ///< pivoting threshold                
               Real eps);         ///< epsilon for zero detection        
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
#endif // _CLUFACTOR_H_
