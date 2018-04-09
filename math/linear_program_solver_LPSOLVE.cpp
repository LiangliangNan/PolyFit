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

		/* Create a new LP model */
		lprec* lp = make_lp(0, variables.size());
		if (!lp) {
			std::cerr << "error in creating a LP model" << std::endl;
			return false;
		}
		set_verbose(lp, SEVERE);

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
			case Variable::FIXED: // value known, actually not a variable 
				set_bounds(lp, i + 1, var.get_bound(), var.get_bound());
				break;

			case Variable::LOWER:
				set_lowbo(lp, i + 1, var.get_bound());
				break;

			case Variable::UPPER:
				set_upbo(lp, i + 1, var.get_bound());
				break;

			case Variable::DOUBLE: {
				double lb = -DBL_MAX;
				double ub = DBL_MAX;
				var.get_bound(lb, ub);
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
		const std::unordered_map<std::size_t, double>& obj_coeffs = objective.coefficients();
		std::unordered_map<std::size_t, double>::const_iterator it = obj_coeffs.begin();
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
			std::vector<int> colno;
			std::vector<double> sparserow;

			const Constraint& cstr = constraints[i];
			const std::unordered_map<std::size_t, double>& cstr_coeffs = cstr.coefficients();
			std::unordered_map<std::size_t, double>::const_iterator cur = cstr_coeffs.begin();
			for (; cur != cstr_coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first;
				double coeff = cur->second;
				colno.push_back(var_idx + 1);	 // lp_solve uses 1-based arrays
				sparserow.push_back(coeff);
			}

			switch (cstr.bound_type())
			{
			case Constraint::FIXED:
				add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), EQ, cstr.get_bound());
				break;
			case Constraint::LOWER:
				add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), GE, cstr.get_bound());
				break;
			case Constraint::UPPER:
				add_constraintex(lp, static_cast<int>(colno.size()), sparserow.data(), colno.data(), LE, cstr.get_bound());
				break;
			case Constraint::DOUBLE: {
				double lb = -DBL_MAX;
				double ub = DBL_MAX;
				cstr.get_bound(lb, ub);
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
			//double objval = get_objective(lp);
			//std::cout << "objective: " << objval << std::endl;
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

