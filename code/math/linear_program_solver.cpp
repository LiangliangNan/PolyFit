/* ---------------------------------------------------------------------------
 * Copyright (C) 2017 Liangliang Nan <liangliang.nan@gmail.com>
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of PolyFit. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 *
 *     Liangliang Nan and Peter Wonka.
 *     PolyFit: Polygonal Surface Reconstruction from Point Clouds.
 *     ICCV 2017.
 *
 *  For more information:
 *  https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html
 * ---------------------------------------------------------------------------
 */


#include <math/linear_program_solver.h>

#include <iostream>
#include <cmath>


bool LinearProgramSolver::check_program(const LinearProgram* program) const {
	if (program->objective()->sense() == LinearObjective::UNDEFINED) {
		std::cerr << "incomplete objective: undefined objective sense." << std::endl;
		return false;
	}

	const std::vector<Variable*>& variables = program->variables();
	if (variables.empty()) {
		std::cerr << "variable set is empty" << std::endl;
		return false;
	}

	// TODO: check if multiple variables have the same name or index

	// TODO: check if multiple constraints have the same name or index

	return true;
}


void LinearProgramSolver::upload_solution(const LinearProgram* program) {
	std::vector<Variable*>& variables = const_cast<LinearProgram*>(program)->variables();
	for (std::size_t i = 0; i < variables.size(); ++i) {
		Variable* v = variables[i];
		v->set_solution_value(result_[i]);
		if (v->variable_type() != Variable::CONTINUOUS)
			result_[i] = static_cast<int>(std::round(result_[i]));
	}
}


bool LinearProgramSolver::solve(const LinearProgram* program, SolverName solver) {
	switch (solver) {
#ifdef HAS_GUROBI
	case GUROBI:
        return _solve_GUROBI(program);
#endif
//	case GLPK:
//        return _solve_GLPK(program);
//	case LPSOLVE:
//        return _solve_LPSOLVE(program);
	case SCIP:
        return _solve_SCIP(program);
	}
    return false;
}

