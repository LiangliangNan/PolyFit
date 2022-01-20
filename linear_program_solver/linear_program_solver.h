/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/ - liangliang.nan@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef SOLVER_SUITE_LINEAR_PROGRAM_SOLVER_H
#define SOLVER_SUITE_LINEAR_PROGRAM_SOLVER_H

#include "linear_program_solver_common.h"

#include <vector>

/***********************************************************************
* I have collected seven different LP/MIP solvers: 
* GUROBI, CBC, SCIP, GLPK, CLP, SOPLEX, LPSOLVE, among which CLP and SOPLEX
* are for continuous problem only (others can solve both LP and MIP problems).
* 
* I have tried my best to tune their options/parameters to achieve their best
* performance. Below is a short report on two problems: 
*    -) the times are in milliseconds and include model construction.
*    -) '-' means too much wait and program was terminated.
* 
**************************************************************
\       Program name: Fig4g (of my PolyFit paper in ICCV 2017)
\       Program type: Binary
\       #Variables: 21556
\               -Binary: 21556
\               -Integer: 0
\               -Continuous: 0
\       #Constraints: 44191
\       #Objective sense: Minimize
\       #Objective offset: 0
**************************************************************
Solved use the MIP solvers (ms):
----------------------------------------------------------------------------
  GUROBI   |   CBC  | SCIP (with CLP) | SCIP (with SoPlex) |  GLPK | LPSOLVE |
-----------|--------|----------------------------------------------|----------
   1520    |   4612 |       5754      |        8357        |   -   |     -   |
----------------------------------------------------------------------------
Solved as LP problem using the LP solvers (ms):
-----------------
   CLP |  SOPLEX |
-----------------|
   156 |   272   |
-----------------

**************************************************************
\       Program name: building
\       Program type: Binary
\       #Variables: 1284
\               -Binary: 1284
\               -Integer: 0
\               -Continuous: 0
\       #Constraints: 2739
\       #Objective sense: Minimize
\       #Objective offset: 0
**************************************************************
Solved use the MIP solvers (ms):
----------------------------------------------------------------------------
  GUROBI   |   CBC  | SCIP (with CLP) | SCIP (with SoPlex) |  GLPK | LPSOLVE |
-----------|--------|----------------------------------------------|----------
    33     |   140  |        87       |         89         |  302  |     -   |
----------------------------------------------------------------------------
Solved as LP problem using the LP solvers (ms):
-----------------
   CLP |  SOPLEX |
-----------------|
    7  |   13    |
-----------------

************************************************************************/
// ToDo:
//   - enable switching between different lp solvers (e.g., clp, soplex) used by SCIP.
//   - check if blas, cholmod, etc. are available for Clp and SCIP
//   - add supports for SOC, SOS1, SOS2, quadratic, NLP constraints and quadratic objectives. 
//	   Check here: https://github.com/SCIP-Interfaces/CSIP
//   - support lazy constraints by implementing a single callback function.
//     


class LinearProgram;

class SOLVER_API LinearProgramSolver
{
public:
    LinearProgramSolver();
    virtual ~LinearProgramSolver();

    /* Solves the problem and returns false if failed.
     *
     * Continuous linear programming solvers.
     *       LP_CLP		// Recommended default value.
     *       LP_SOPLEX
     * (Mixed) Integer programming solvers.
     *       CBC:       // Recommended default value.
     *       SCIP,
     *       GLPK,
     *       LPSOLVE,
     *       GUROBI: best performance, but it requires a license :-(
     *
     * NOTE: Both CBC and SCIP are (approximately 5-10 times) slower than Gurobi.
     *       GLPK and LPSOLVE are too slow and may fail for problems with more
     *       than three hundred variables. If you have a really LARGE problem,
     *       please consider using Gurobi.
     */
    virtual bool solve(const LinearProgram* program) = 0;

    /* Returns the result.
     * The result can also be retrieved using Variable::solution_value().
     * NOTE: (1) result is valid only if the solver succeeded.
     *       (2) each entry in the result corresponds to the variable with the
     * 			 same index in the linear program.
     */
	const std::vector<double>& solution() const { return result_; }

    /*  Returns the objective value.
     *  NOTE: (1) result is valid only if the solver succeeded.
     * 		  (2) if a variable is integer, then rounded values are used.
     *        (3) the constant term is also included.
     */
	double objective_value() const { return objective_value_; }

protected:
	// write the solution values to the program
	void upload_solution(const LinearProgram* program);

protected:
	std::vector<double> result_;
	double				objective_value_;

private:
    //copying disabled
    LinearProgramSolver(const LinearProgramSolver&);
    LinearProgramSolver& operator=(const LinearProgramSolver&);
};

#endif  // SOLVER_SUITE_LINEAR_PROGRAM_SOLVER_H
