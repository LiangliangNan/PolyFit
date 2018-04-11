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
#include "../3rd_glpk/glpk.h"

#include <iostream>


bool LinearProgramSolver::_solve_GLPK(const LinearProgram* program) {
	try {
		typedef Variable<double>			Variable;
		typedef LinearExpression<double>	Objective;
		typedef LinearConstraint<double>	Constraint;

		const std::vector<Variable>& variables = program->variables();
		if (variables.empty()) {
			std::cerr << "variable set is empty" << std::endl;
			return false;
		}

		glp_prob* lp = glp_create_prob();
		if (!lp) {
			std::cerr << "error in creating a LP model" << std::endl;
			return false;
		}
		glp_set_prob_name(lp, "unknown");
		glp_set_obj_dir(lp, GLP_MIN);

		std::size_t num_integer_variables = 0;

		// create variables
		glp_add_cols(lp, variables.size());
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable& var = variables[i];
			const std::string& name = "x" + std::to_string(i + 1);
			glp_set_col_name(lp, i + 1, name.data());

			if (var.variable_type() == Variable::INTEGER) {
				glp_set_col_kind(lp, i + 1, GLP_IV);	// glpk uses 1-based arrays
				++num_integer_variables;
			}
			else if (var.variable_type() == Variable::BINARY) {
				glp_set_col_kind(lp, i + 1, GLP_BV);	// glpk uses 1-based arrays
				++num_integer_variables;
			}
			else {
				glp_set_col_kind(lp, i + 1, GLP_CV);	// continuous variable
			}

			int bound_type = GLP_FR;
			switch (var.bound_type())
			{
			case Variable::FIXED:  bound_type = GLP_FX; break;
			case Variable::LOWER:  bound_type = GLP_LO; break;
			case Variable::UPPER:  bound_type = GLP_UP; break;
			case Variable::DOUBLE: bound_type = GLP_DB; break;
			case Variable::FREE:
			default:
				break;
			}

			double lb, ub;
			var.get_double_bounds(lb, ub);
			glp_set_col_bnds(lp, i + 1, bound_type, lb, ub);
		}

		// set objective 

		const Objective& objective = program->objective();
		const std::unordered_map<std::size_t, double>& obj_coeffs = objective.coefficients();
		std::unordered_map<std::size_t, double>::const_iterator it = obj_coeffs.begin();
		for (; it != obj_coeffs.end(); ++it) {
			std::size_t var_idx = it->first;
			double coeff = it->second;
			glp_set_obj_coef(lp, var_idx + 1, coeff); // glpk uses 1-based arrays
		}

		// Add constraints

		const std::vector<Constraint>& constraints = program->constraints();
		glp_add_rows(lp, constraints.size());

		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const Constraint& cstr = constraints[i];
			const std::unordered_map<std::size_t, double>& cstr_coeffs = cstr.coefficients();
			std::unordered_map<std::size_t, double>::const_iterator cur = cstr_coeffs.begin();

			std::vector<int>	colno(cstr_coeffs.size() + 1, 0);		// glpk uses 1-based arrays
			std::vector<double> sparserow(cstr_coeffs.size() + 1, 0.0); // glpk uses 1-based arrays
			std::size_t idx = 1; // glpk uses 1-based arrays
			for (; cur != cstr_coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first;
				double coeff = cur->second;

				colno[idx] = var_idx + 1;	 // glpk uses 1-based arrays
				sparserow[idx] = coeff;
				++idx;
			}

			glp_set_mat_row(lp, i + 1, static_cast<int>(cstr_coeffs.size()), colno.data(), sparserow.data());

			int bound_type = GLP_FR;
			switch (cstr.bound_type())
			{
			case Constraint::FIXED:  bound_type = GLP_FX; break;
			case Constraint::LOWER:  bound_type = GLP_LO; break;
			case Constraint::UPPER:  bound_type = GLP_UP; break;
			case Constraint::DOUBLE: bound_type = GLP_DB; break;
			case Constraint::FREE:
			default:
				break;
			}

			double lb, ub;
			cstr.get_double_bounds(lb, ub);
			glp_set_row_bnds(lp, i + 1, bound_type, lb, ub);
		}

		int msg_level = GLP_MSG_ERR;
		int status = -1;
		if (num_integer_variables == 0) { // continuous problem
			glp_smcp parm;
			glp_init_smcp(&parm);
			parm.msg_lev = msg_level;
			status = glp_simplex(lp, &parm);
		}
		else { // solve as MIP problem
			glp_iocp parm;
			glp_init_iocp(&parm);
			parm.msg_lev = msg_level;
			parm.presolve = GLP_ON;
			// The routine glp_intopt is a driver to the MIP solver based on the branch-and-cut method,
			// which is a hybrid of branch-and-bound and cutting plane methods.
			status = glp_intopt(lp, &parm);	
		}

		switch (status) {
		case 0: {
			if (num_integer_variables == 0) { // continuous problem
				objective_value_ = glp_get_obj_val(lp);
				assert(variables.size() == glp_get_num_cols(lp));
				result_.resize(variables.size());
				for (std::size_t i = 0; i < variables.size(); ++i) {
					result_[i] = glp_get_col_prim(lp, i + 1);
				}
			}
			else { // MIP problem
				objective_value_ = glp_mip_obj_val(lp);
				assert(variables.size() == glp_get_num_cols(lp));
				result_.resize(variables.size());
				for (std::size_t i = 0; i < variables.size(); ++i) {
					result_[i] = glp_mip_col_val(lp, i + 1);
				}
			}
			break;
		}

		case GLP_EBOUND:
			std::cerr << 
				"Unable to start the search, because some double-bounded variables have incorrect"
				"bounds or some integer variables have non - integer(fractional) bounds." << std::endl;
			break;

		case GLP_EROOT:
			std::cerr << 
				"Unable to start the search, because optimal basis for initial LP relaxation is not"
				"provided. (This code may appear only if the presolver is disabled.)" << std::endl;
			break;

		case GLP_ENOPFS:
			std::cerr << 
				"Unable to start the search, because LP relaxation of the MIP problem instance has"
				"no primal feasible solution. (This code may appear only if the presolver is enabled.)" << std::endl;
			break;

		case GLP_ENODFS:
			std::cerr << 
				"Unable to start the search, because LP relaxation of the MIP problem instance has"
				"no dual feasible solution.In other word, this code means that if the LP relaxation"
				"has at least one primal feasible solution, its optimal solution is unbounded, so if the"
				"MIP problem has at least one integer feasible solution, its(integer) optimal solution"
				"is also unbounded. (This code may appear only if the presolver is enabled.)" << std::endl;
			break;

		case GLP_EFAIL:
			std::cerr << "The search was prematurely terminated due to the solver failure." << std::endl;
			break;

		case GLP_EMIPGAP:
			std::cerr << 
				"The search was prematurely terminated, because the relative mip gap tolerance has been reached." << std::endl;
			break;

		case GLP_ETMLIM:
			std::cerr << "The search was prematurely terminated, because the time limit has been exceeded." << std::endl;
			break;

		case GLP_ESTOP:
			std::cerr << 
				"The search was prematurely terminated by application. (This code may appear only"
				"if the advanced solver interface is used.)" << std::endl;
			break;

		default:
			std::cerr << "optimization was stopped with status = " << status << std::endl;
			break;
		}

		glp_delete_prob(lp);

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