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


#include "linear_program_solver.h"

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
	case GLPK:
        return _solve_GLPK(program);
	case LPSOLVE:
        return _solve_LPSOLVE(program);
	case SCIP:
        return _solve_SCIP(program);
	}
    return false;
}

