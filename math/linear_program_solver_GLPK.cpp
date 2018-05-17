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
#include "../basic/logger.h"
#include "../3rd_glpk/glpk.h"

#include <iostream>


bool LinearProgramSolver::_solve_GLPK(const LinearProgram* program) {
	try {
		if (!check_program(program))
			return false;

		glp_prob* lp = glp_create_prob();
		if (!lp) {
			std::cerr << "error in creating a LP model" << std::endl;
			return false;
		}
		glp_set_prob_name(lp, program->name().c_str());

		std::size_t num_integer_variables = 0;

		// create variables
		const std::vector<Variable*>& variables = program->variables();
		glp_add_cols(lp, variables.size());
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable* var = variables[i];
			glp_set_col_name(lp, i + 1, var->name().c_str());

			if (var->variable_type() == Variable::INTEGER) {
				glp_set_col_kind(lp, i + 1, GLP_IV);	// glpk uses 1-based arrays
				++num_integer_variables;
			}
			else if (var->variable_type() == Variable::BINARY) {
				glp_set_col_kind(lp, i + 1, GLP_BV);	// glpk uses 1-based arrays
				++num_integer_variables;
			}
			else {
				glp_set_col_kind(lp, i + 1, GLP_CV);	// continuous variable
			}

			int bound_type = GLP_FR;
			switch (var->bound_type())
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
			var->get_bounds(lb, ub);
			glp_set_col_bnds(lp, i + 1, bound_type, lb, ub);
		}

		// Add constraints

		const std::vector<LinearConstraint*>& constraints = program->constraints();
		glp_add_rows(lp, constraints.size());

		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const LinearConstraint* c = constraints[i];
			const std::unordered_map<int, double>& coeffs = c->coefficients();
			std::unordered_map<int, double>::const_iterator cur = coeffs.begin();

			std::vector<int>	indices(coeffs.size() + 1, 0);		// glpk uses 1-based arrays
			std::vector<double> coefficients(coeffs.size() + 1, 0.0);  // glpk uses 1-based arrays
			std::size_t idx = 1; // glpk uses 1-based arrays
			for (; cur != coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first;
				double coeff = cur->second;

				indices[idx] = var_idx + 1;	 // glpk uses 1-based arrays
				coefficients[idx] = coeff;
				++idx;
			}

			glp_set_mat_row(lp, i + 1, static_cast<int>(coeffs.size()), indices.data(), coefficients.data());

			int bound_type = GLP_FR;
			switch (c->bound_type())
			{
			case LinearConstraint::FIXED:  bound_type = GLP_FX; break;
			case LinearConstraint::LOWER:  bound_type = GLP_LO; break;
			case LinearConstraint::UPPER:  bound_type = GLP_UP; break;
			case LinearConstraint::DOUBLE: bound_type = GLP_DB; break;
			case LinearConstraint::FREE:
			default:
				break;
			}

			double lb, ub;
			c->get_bounds(lb, ub);
			glp_set_row_bnds(lp, i + 1, bound_type, lb, ub);

			glp_set_row_name(lp, i + 1, c->name().c_str());
		}

		// set objective 

		// determine the coefficient of each variable in the objective function
		const LinearObjective* objective = program->objective();
		const std::unordered_map<int, double>& obj_coeffs = objective->coefficients();
		std::unordered_map<int, double>::const_iterator cur = obj_coeffs.begin();
		for (; cur != obj_coeffs.end(); ++cur) {
			std::size_t var_idx = cur->first;
			double coeff = cur->second;
			glp_set_obj_coef(lp, var_idx + 1, coeff); // glpk uses 1-based arrays
		}

		Logger::out("-") << "using the GLPK solver" << std::endl;

		// Set objective function sense
		bool minimize = (objective->sense() == LinearObjective::MINIMIZE);
		glp_set_obj_dir(lp, minimize ? GLP_MIN : GLP_MAX);
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
				result_.resize(variables.size());
				for (std::size_t i = 0; i < variables.size(); ++i) {
					result_[i] = glp_get_col_prim(lp, i + 1);	 // glpk uses 1-based arrays
				}
			}
			else { // MIP problem
				objective_value_ = glp_mip_obj_val(lp);
				result_.resize(variables.size());
				for (std::size_t i = 0; i < variables.size(); ++i) {
					result_[i] = glp_mip_col_val(lp, i + 1);	 // glpk uses 1-based arrays
				}
			}
			upload_solution(program);
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