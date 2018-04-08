#include "linear_program_solver.h"
#include "../basic/logger.h"


bool LinearProgramSolver::solve(const LinearProgram* program, LP_Solver solver /* = GUROBI */) {
	switch (solver) {
	case GUROBI:
		return _solve_GUROBI(program);
	case SCIP:
		return _solve_SCIP(program);
	case LPSOLVE:
		return _solve_LPSOLVE(program);
	case GLPK:
		return _solve_GLPK(program);
	default:
		Logger::warn("-") << "no such solver" << std::endl;
		return false;
	}
}

