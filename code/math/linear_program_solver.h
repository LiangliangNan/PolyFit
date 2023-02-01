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

#ifndef _MATH_LINEAR_PROGRAM_SOLVER_H_
#define _MATH_LINEAR_PROGRAM_SOLVER_H_

#include "math_common.h"
#include "linear_program.h"

#include <vector>


class MATH_API LinearProgramSolver
{
public:
	enum SolverName { 
#ifdef HAS_GUROBI	// Gurobi is commercial and requires license :-(
		GUROBI,	
#endif
		SCIP,		// Recommended default value.
		GLPK,
		LPSOLVE,
	};

public:
	LinearProgramSolver() {}
	~LinearProgramSolver() {}

	// Solves the problem and returns false if fails.
	// NOTE: The SCIP and CBC solvers are slower than Gurobi but acceptable. 
	//       GLPK and LPSOLVE may be too slow or even fail for large problems.
	//       If you have a really LARGE problem, you may consider using Gurobi.
    bool solve(const LinearProgram* program, SolverName solver);

	// Returns the result. 
	// The result can also be retrieved using Variable::solution_value().
	// NOTE: (1) result is valid only if the solver succeeded.
	//       (2) each entry in the result corresponds to the variable with the
	//			 same index in the linear program.
	const std::vector<double>& solution() const { return result_; }

	// Returns the objective value.
	// NOTE: (1) result is valid only if the solver succeeded.
	//       (2) the constant term is not included.
	double objective_value() const { return objective_value_; }

private:
	bool check_program(const LinearProgram* program) const;
	void upload_solution(const LinearProgram* program);

private:
#ifdef HAS_GUROBI
	bool _solve_GUROBI(const LinearProgram* program);
#endif
	bool _solve_SCIP(const LinearProgram* program);
	bool _solve_GLPK(const LinearProgram* program);
	bool _solve_LPSOLVE(const LinearProgram* program);

private:
	std::vector<double> result_;
	double				objective_value_;
};

#endif
