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

#define WIN32_LEAN_AND_MEAN 1

#include "linear_program_solver.h"
#include "../basic/basic_types.h"


// Liangliang: CBC has many dependences. I am lazy, so I am using the 
//             google ortools' CBC wrapper class. 

#ifdef HAS_CBC_SOLVER

#include "F:/3_code/SolverSuite/ortools/linear_solver/linear_solver.h"

#ifdef _DEBUG
#pragma comment(lib, "F:/3_code/SolverSuite/coin/build/x64/Debug/linear_solver.lib")
#else
#pragma comment(lib, "F:/3_code/SolverSuite/coin/build/x64/Release/linear_solver.lib")
#endif

using namespace operations_research;

bool LinearProgramSolver::_solve_CBC(const LinearProgram* program) {
	try {
		typedef Variable<double>			Variable;
		typedef LinearExpression<double>	Objective;
		typedef LinearConstraint<double>	Constraint;

		const std::vector<Variable>& variables = program->variables();
		if (variables.empty()) {
			std::cerr << "variable set is empty" << std::endl;
			return false;
		}

		MPSolver solver("CBC_Solver", MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);

		// create variables
		std::vector<MPVariable*> X(variables.size());
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable& var = variables[i];
			const std::string& name = "x" + std::to_string(i);

			double lb, ub;
			var.get_double_bounds(lb, ub);

			switch (var.variable_type())
			{
			case Variable::CONTINUOUS:
				X[i] = solver.MakeNumVar(lb, ub, name.data());
				break;
			case Variable::INTEGER:
				X[i] = solver.MakeIntVar(lb, ub, name.data());
				break;
			case Variable::BINARY:
				X[i] = solver.MakeBoolVar(name.data());
				break;
			}
		}

		MPObjective* const cbc_objective = solver.MutableObjective();
		cbc_objective->SetOptimizationDirection(program->objective_sense() == LinearProgram::MAXIMIZE);

		const Objective& objective = program->objective();
		const std::unordered_map<std::size_t, double>& obj_coeffs = objective.coefficients();
		std::unordered_map<std::size_t, double>::const_iterator it = obj_coeffs.begin();
		for (; it != obj_coeffs.end(); ++it) {
			std::size_t var_idx = it->first;
			double coeff = it->second;
			cbc_objective->SetCoefficient(X[var_idx], coeff);
		}

		// Add constraints
		const std::vector<Constraint>& constraints = program->constraints();
		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const Constraint& cstr = constraints[i];
			const std::string& name = "cstr" + std::to_string(i);
			MPConstraint* const c = solver.MakeRowConstraint();

			const std::unordered_map<std::size_t, double>& cstr_coeffs = cstr.coefficients();
			std::unordered_map<std::size_t, double>::const_iterator cur = cstr_coeffs.begin();
			for (; cur != cstr_coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first;
				double coeff = cur->second;
				c->SetCoefficient(X[var_idx], coeff);
			}

			switch (cstr.bound_type())
			{
			case Constraint::FIXED:
				c->SetBounds(cstr.get_single_bound(), cstr.get_single_bound());
				break;
			case Constraint::LOWER:
				c->SetLB(cstr.get_single_bound());
				break;
			case Constraint::UPPER:
				c->SetUB(cstr.get_single_bound());
				break;
			case Constraint::DOUBLE: {
				double lb, ub;
				cstr.get_double_bounds(lb, ub);
				c->SetBounds(lb, ub);
				break;
				}
			default:
				break;
			}
		}

		// Optimize model
		const MPSolver::ResultStatus result_status = solver.Solve();

		// Check that the problem has an optimal solution.
		if (result_status == MPSolver::OPTIMAL) {
			objective_value_ = cbc_objective->Value();
			result_.resize(variables.size());
			for (std::size_t i = 0; i < variables.size(); ++i) {
				result_[i] = X[i]->solution_value();
			}
		}
		else {
			std::cerr << "The problem does not have an optimal solution!" << std::endl;
		}

		return (result_status == MPSolver::OPTIMAL);
	}
	catch (std::exception e) {
		std::cerr << "Error code = " << e.what() << std::endl;
	}
	catch (...) {
		std::cerr << "Exception during optimization" << std::endl;
	}

	return false;
}

#endif