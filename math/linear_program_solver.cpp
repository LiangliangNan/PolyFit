#include "linear_program_solver.h"


bool LinearProgramSolver::solve(const LinearProgram* program, SolverName solver /* = GUROBI */) {
	result_.clear();

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
		std::cerr << "no such solver" << std::endl;
		return false;
	}
}

