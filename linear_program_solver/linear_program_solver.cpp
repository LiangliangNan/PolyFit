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


#include <linear_program_solver/linear_program_solver.h>

#include <cmath>    // for std::round()

#include <linear_program_solver/linear_program.h>


LinearProgramSolver::LinearProgramSolver()
{

}


LinearProgramSolver::~LinearProgramSolver()
{

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
