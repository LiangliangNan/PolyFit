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
#include "../3rd_lpsolve/lp_lib.h"
#include "../basic/logger.h"

#include <iostream>


bool LinearProgramSolver::_solve_LPSOLVE(const LinearProgram* program) {
	try {
		if (!check_program(program))
			return false;

		// Note that there are several restrictions with its row entry mode:
		// Only use set_add_rowmode() after a make_lp call. Also, if this function is used, first add the objective function 
		// and after that add the constraints. Don't call other API functions while in row entry mode. No other data matrix 
		// access is allowed while in row entry mode. After adding the constraints, turn row entry mode back off. Once turned 
		// off, you cannot switch back to row entry mode. So in short:
		//	- turn row entry mode on
		//	- set the objective function
		//	- create the constraints
		//	- turn row entry mode off

		// Create a new LP model
		const std::vector<Variable*>& variables = program->variables();
		const std::vector<LinearConstraint*>& constraints = program->constraints();
		lprec* lp = make_lp(constraints.size(), variables.size());
		if (!lp) {
			std::cerr << "error in creating a LP model" << std::endl;
			return false;
		}

		set_lp_name(lp, const_cast<char*>(program->name().c_str()));
		set_verbose(lp, SEVERE);
		set_presolve(lp, PRESOLVE_ROWS | PRESOLVE_COLS | PRESOLVE_LINDEP, get_presolveloops(lp));

		// create variables
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable* var = variables[i];

			std::size_t var_idx = i + 1;	// The LP_SOLVE manual says the first element is ignored
			if (var->variable_type() == Variable::INTEGER)
				set_int(lp, var_idx, TRUE);
			else if (var->variable_type() == Variable::BINARY)
				set_binary(lp, var_idx, TRUE);
			else {
				// continuous variable
			}

			switch (var->bound_type())
			{
			case Variable::FIXED: { // value known, actually not a variable 
				double bd = var->get_bound();
				set_bounds(lp, var_idx, bd, bd);
				break;
			}
			case Variable::LOWER:
				set_lowbo(lp, var_idx, var->get_bound());
				break;
			case Variable::UPPER:
				set_upbo(lp, var_idx, var->get_bound());
				break;
			case Variable::DOUBLE: {
				double lb, ub;
				var->get_bounds(lb, ub);
				set_bounds(lp, var_idx, lb, ub);
				break;
			}
			case Variable::FREE:
			default:
				set_unbounded(lp, var_idx);
				break;
			}
		}

		// turn row entry mode on
		set_add_rowmode(lp, TRUE);

		// set objective 

		// The LP_SOLVE manual says the first element is ignored, so any value is OK.
		std::vector<double> row(variables.size() + 1, 0);

		// determine the coefficient of each variable in the objective function
		const LinearObjective* objective = program->objective();
		const std::unordered_map<int, double>& obj_coeffs = objective->coefficients();
		std::unordered_map<int, double>::const_iterator cur = obj_coeffs.begin();
		for (; cur != obj_coeffs.end(); ++cur) {
			std::size_t var_idx = cur->first;
			double coeff = cur->second;
			row[var_idx + 1] = coeff;	// The LP_SOLVE manual says the first element is ignored
		}
		set_obj_fn(lp, row.data());

		// Set objective function sense
		set_sense(lp, objective->sense() == LinearObjective::MAXIMIZE); // true for maximize

		// Add constraints

		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const LinearConstraint* c = constraints[i];
			const std::unordered_map<int, double>& coeffs = c->coefficients();
			std::unordered_map<int, double>::const_iterator cur = coeffs.begin();

			// Liangliang: Annoying LPSOLVE: some functions read an array from 0 but some from 1!!!
			// set_rowex() is one of the functions that read arrays forom 0.
			// The LP_SOLVE manual: In contrary to set_row(), set_rowex() reads the arrays starting from element 0.
			std::vector<int>	indices(coeffs.size(), 0);			
			std::vector<double> coefficients(coeffs.size(), 0.0);
			std::size_t idx = 0; // lp_solve uses 1-based arrays
			for (; cur != coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first;
				double coeff = cur->second;

				indices[idx] = var_idx + 1;	 // The LP_SOLVE manual says the first element is ignored
				coefficients[idx] = coeff;
				++idx;
			}

			// set the coefficients
			std::size_t row_idx = i + 1;	// The LP_SOLVE manual says the first element is ignored
			set_rowex(lp, row_idx, static_cast<int>(indices.size()), coefficients.data(), indices.data());
			switch (c->bound_type())
			{
			case LinearConstraint::FIXED:
				set_constr_type(lp, row_idx, EQ);
				set_rh(lp, row_idx, c->get_bound());
				break;
			case LinearConstraint::LOWER:
				set_constr_type(lp, row_idx, GE);
				set_rh(lp, row_idx, c->get_bound());
				break;
			case LinearConstraint::UPPER:
				set_constr_type(lp, row_idx, LE);
				set_rh(lp, row_idx, c->get_bound());
				break;
			case LinearConstraint::DOUBLE: {
				double lb, ub;
				c->get_bounds(lb, ub);
				set_constr_type(lp, row_idx, GE); // I choose GE and I will set the range using set_rh_range()
				set_rh(lp, row_idx, lb);
				set_rh_range(lp, row_idx, ub - lb);
				break;
				}
			default:
				break;
			}
		}

		// turn row entry mode off
		set_add_rowmode(lp, FALSE);

		Logger::out("-") << "using the LPSOLVE solver" << std::endl;
		int status = ::solve(lp);
		switch (status) {
		case 0: {
			objective_value_ = get_objective(lp);
			result_.resize(variables.size());
			get_variables(lp, result_.data());
			upload_solution(program);
			break;
		}
		case -2:
			std::cerr << "Out of memory" << std::endl;
			break;
		case 1:
			std::cerr << "The model is sub-optimal. Only happens if there are integer variables and there is already an integer solution found. The solution is not guaranteed the most optimal one." << std::endl;
			break;
		case 2:
			std::cerr << "The model is infeasible" << std::endl;
			break;
		case 3:
			std::cerr << "The model is unbounded" << std::endl;
			break;
		case 4:
			std::cerr << "The model is degenerative" << std::endl;
			break;
		case 5:
			std::cerr << "Numerical failure encountered" << std::endl;
			break;
		case 6:
			std::cerr << "The abort() routine was called" << std::endl;
			break;
		case 7:
			std::cerr << "A timeout occurred" << std::endl;
			break;
		case 9:
			std::cerr << "The model could be solved by presolve. This can only happen if presolve is active via set_presolve()" << std::endl;
			break;
		case 25:
			std::cerr << "Accuracy error encountered" << std::endl;
			break;
		default:
			std::cerr << "optimization was stopped with status = " << status << std::endl;
			break;
		}

		delete_lp(lp);

		return (status == 0);
	}
	catch (std::exception e) {
		std::cerr << "Error code = " << e.what() << std::endl;
	}
	catch (...) {
		std::cerr << "Exception during optimization" << std::endl;
	}

	return false;
}

