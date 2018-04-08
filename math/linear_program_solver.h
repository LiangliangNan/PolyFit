#pragma once

#include "math_common.h"
#include "linear_program.h"

#include <vector>

class MATH_API LinearProgramSolver
{
public:
	enum LP_Solver { GUROBI, SCIP, LPSOLVE, GLPK };

	typedef LinearProgram<double>	LinearProgram;

public:
	LinearProgramSolver() {}
	~LinearProgramSolver() {}

	// Solve the problem; returns false if fails
	// NOTE: Gurobi solver recommended;
	//		 The SCIP solver is about x10 slower than Gurobi;
	//       LPSOLVE and GLPK may be too slow or even fail.
	bool solve(const LinearProgram* program, LP_Solver solver = GUROBI);

	// returns the result
	// NOTE: (1) result is valid only if the solver succeeded
	//       (2) the result includes all auxiliary variables
	const std::vector<double>& get_result() const { 
		return result_; 
	}

private:
	std::vector<double> result_;
};

