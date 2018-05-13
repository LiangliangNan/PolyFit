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
#include "../basic/basic_types.h"
#include "../3rd_lpsolve/lp_lib.h"

#include <iostream>


bool LinearProgramSolver::_solve_LPSOLVE(const LinearProgram* program) {
	try {
		typedef Variable<double>			Variable;
		typedef LinearExpression<double>	Objective;
		typedef LinearConstraint<double>	Constraint;

		const std::vector<Variable>& variables = program->variables();
		if (variables.empty()) {
			std::cerr << "variable set is empty" << std::endl;
			return false;
		}

		// Create a new LP model
		lprec* lp = make_lp(0, variables.size());
		if (!lp) {
			std::cerr << "error in creating a LP model" << std::endl;
			return false;
		}
		set_verbose(lp, SEVERE);

		// set object sense
		set_sense(lp, program->objective_sense() == LinearProgram::MAXIMIZE); // true for maximize

		// create variables
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable& var = variables[i];

			if (var.variable_type() == Variable::INTEGER)
				set_int(lp, i+1, TRUE);		// lp_solve uses 1-based arrays
			else if (var.variable_type() == Variable::BINARY)
				set_binary(lp, i+1, TRUE);	// lp_solve uses 1-based arrays
			else {
				// continuous variable
			}

			switch (var.bound_type())
			{
			case Variable::FIXED: { // value known, actually not a variable 
				double bd = var.get_single_bound();
				set_bounds(lp, i + 1, bd, bd);
				break;
			}
			case Variable::LOWER:
				set_lowbo(lp, i + 1, var.get_single_bound());
				break;
			case Variable::UPPER:
				set_upbo(lp, i + 1, var.get_single_bound());
				break;
			case Variable::DOUBLE: {
				double lb = -DBL_MAX;
				double ub = DBL_MAX;
				var.get_double_bounds(lb, ub);
				set_bounds(lp, i + 1, lb, ub);
				break;
			}
			case Variable::FREE:
			default:
				set_unbounded(lp, i + 1);
				break;
			}
		}

		// set objective 

		// lp_solve uses 1-based arrays
		// The LP_SOLVE manual says the first element is ignored, so any value is OK.
		std::vector<double> row(variables.size() + 1, 0);
									
		const Objective& objective = program->objective();
		const std::map<std::size_t, double>& obj_coeffs = objective.coefficients();
		std::map<std::size_t, double>::const_iterator it = obj_coeffs.begin();
		for (; it != obj_coeffs.end(); ++it) {
			std::size_t var_idx = it->first;
			double coeff = it->second;	
			row[var_idx + 1] = coeff;	// lp_solve uses 1-based arrays
		}
		set_obj_fn(lp, row.data());

		// Add constraints

		const std::vector<Constraint>& constraints = program->constraints();

		/* num constraints will be added, so allocate memory for it in advance to make things faster */
		resize_lp(lp, constraints.size(), get_Ncolumns(lp));

		set_add_rowmode(lp, TRUE);

		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const Constraint& cstr = constraints[i];
			const std::map<std::size_t, double>& cstr_coeffs = cstr.coefficients();
			std::map<std::size_t, double>::const_iterator cur = cstr_coeffs.begin();

			std::vector<int>	colno(cstr_coeffs.size() + 1, 0);		// lp_solve uses 1-based arrays
			std::vector<double> sparserow(cstr_coeffs.size() + 1, 0.0); // lp_solve uses 1-based arrays
			std::size_t idx = 1; // lp_solve uses 1-based arrays
			for (; cur != cstr_coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first;
				double coeff = cur->second;

				colno[idx] = var_idx + 1;	 // lp_solve uses 1-based arrays
				sparserow[idx] = coeff;
				++idx;
			}

			switch (cstr.bound_type())
			{
			case Constraint::FIXED:
				add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), EQ, cstr.get_single_bound());
				break;
			case Constraint::LOWER:
				add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), GE, cstr.get_single_bound());
				break;
			case Constraint::UPPER:
				add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), LE, cstr.get_single_bound());
				break;
			case Constraint::DOUBLE: {
				double lb = -DBL_MAX;
				double ub = DBL_MAX;
				cstr.get_double_bounds(lb, ub);
				add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), GE, lb);
				add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), LE, ub);
				break;
				}
			default:
				break;
			}
		}

		set_add_rowmode(lp, FALSE);
		int status = ::solve(lp);
		switch (status) {
		case 0: {
			objective_value_ = get_objective(lp);
			result_.resize(variables.size());
			get_variables(lp, result_.data());
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

