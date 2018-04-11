#pragma once

#include "math_common.h"
#include "linear_program.h"

#include <vector>

class MATH_API LinearProgramSolver
{
public:
	enum SolverName { GUROBI, SCIP, GLPK, LPSOLVE };

	typedef LinearProgram<double>	LinearProgram;

public:
	LinearProgramSolver() {}
	~LinearProgramSolver() {}

	// Solve the problem; returns false if fails
	// NOTE: Gurobi solver recommended;
	//		 The SCIP solver is about x10 slower than Gurobi;
	//       LPSOLVE and GLPK may be too slow or even fail.
	bool solve(const LinearProgram* program, SolverName solver = GUROBI);

	// returns the objective value
	// NOTE: result is valid only if the solver succeeded
	double get_objective_value() const { return objective_value_; }

	// returns the result
	// NOTE: (1) result is valid only if the solver succeeded
	//       (2) the result includes all auxiliary variables
	const std::vector<double>& get_result() const { 
		return result_; 
	}

private:
	bool _solve_GUROBI(const LinearProgram* program);
	bool _solve_SCIP(const LinearProgram* program);
	bool _solve_LPSOLVE(const LinearProgram* program);
	bool _solve_GLPK(const LinearProgram* program);

private:
	std::vector<double> result_;
	double				objective_value_;
};

