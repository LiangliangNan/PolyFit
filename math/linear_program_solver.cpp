/*
Copyright (C) 2017  Liangliang Nan
http://web.siat.ac.cn/~liangliang/ - liangliang.nan@gmail.com

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


bool LinearProgramSolver::solve(const LinearProgram* program, SolverName solver) {
	if (program->objective_sense() == LinearProgram::UNDEFINED) {
		std::cerr << "incomplete objective: undefined objective sense." << std::endl;
		return false;
	}

	result_.clear();

	switch (solver) {
#ifdef HAS_GUROBI_SOLVER
	case GUROBI:
		return _solve_GUROBI(program);
#endif
#ifdef HAS_CBC_SOLVER
	case CBC:
		return _solve_CBC(program);
#endif
	case GLPK:
		return _solve_GLPK(program);
	case LPSOLVE:
		return _solve_LPSOLVE(program);
#ifdef HAS_SCIP_SOLVER
	case SCIP:
	default: // use SCIP
		return _solve_SCIP(program);
#endif
	}

    std::cerr << "no such solver doesn't exist" << std::endl;
    return false;
}

