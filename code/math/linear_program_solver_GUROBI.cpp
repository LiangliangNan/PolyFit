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


#ifdef HAS_GUROBI

#include <gurobi_c++.h>

#ifdef WIN32
#if (_MSC_VER == 1900) // vs 2015
#pragma comment(lib, "gurobi80.lib")
#ifdef _DEBUG
#pragma comment(lib, "gurobi_c++mdd2015.lib")
#else
#pragma comment(lib, "gurobi_c++md2015.lib")
#endif
#elif (_MSC_VER >= 1910 && _MSC_VER <= 1915) // vs 2017
#pragma comment(lib, "gurobi81.lib")
#ifdef _DEBUG
#pragma comment(lib, "gurobi_c++mdd2017.lib")
#else
#pragma comment(lib, "gurobi_c++md2017.lib")
#endif
#endif
#endif


bool LinearProgramSolver::_solve_GUROBI(const LinearProgram* program) {
	try {
		if (!check_program(program))
			return false;

		// I am using an academic license of Gurobi. Each time when a Gurobi environment is created, it pops up a
		// notice "Academic license - for non-commercial use only".
		// It is not possible to suppress this notice completely, but we can get the Gurobi environment only once and
		// reuse it later on.
		static GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);

		GRBModel model = GRBModel(env);

		// create variables
		const std::vector<Variable*>& variables = program->variables();
		std::vector<GRBVar> X(variables.size());
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable* var = variables[i];

			double lb, ub;
			var->get_bounds(lb, ub);

			char vtype = GRB_CONTINUOUS;
			if (var->variable_type() == Variable::INTEGER)
				vtype = GRB_INTEGER;
			else if (var->variable_type() == Variable::BINARY)
				vtype = GRB_BINARY;

			X[i] = model.addVar(lb, ub, 0.0, vtype);
		}

		// Integrate new variables
		model.update();

		// Add constraints
		const std::vector<LinearConstraint*>& constraints = program->constraints();
		for (std::size_t i = 0; i < constraints.size(); ++i) {
			GRBLinExpr expr;
			const LinearConstraint* c = constraints[i];
			const std::unordered_map<int, double>& coeffs = c->coefficients();
			std::unordered_map<int, double>::const_iterator cur = coeffs.begin();
			for (; cur != coeffs.end(); ++cur) {
                std::size_t var_idx = static_cast<std::size_t>(cur->first);
				double coeff = cur->second;
				expr += coeff * X[var_idx];
			}

			switch (c->bound_type())
			{
			case LinearConstraint::FIXED:
				model.addConstr(expr == c->get_bound());
				break;
			case LinearConstraint::LOWER:
				model.addConstr(expr >= c->get_bound());
				break;
			case LinearConstraint::UPPER:
				model.addConstr(expr <= c->get_bound());
				break;
			case LinearConstraint::DOUBLE: {
				double lb, ub;
				c->get_bounds(lb, ub);
				model.addConstr(expr >= lb);
				model.addConstr(expr <= ub);
				break;
				}
			default:
				break;
			}
		}

		// Set objective
		GRBLinExpr obj;

		const LinearObjective* objective = program->objective();
		const std::unordered_map<int, double>& obj_coeffs = objective->coefficients();
		std::unordered_map<int, double>::const_iterator it = obj_coeffs.begin();
		for (; it != obj_coeffs.end(); ++it) {
            std::size_t var_idx = static_cast<std::size_t>(it->first);
			double coeff = it->second;
			obj += coeff * X[var_idx];
		}
		// Set objective function sense
		bool minimize = (objective->sense() == LinearObjective::MINIMIZE);
		model.setObjective(obj, minimize ? GRB_MINIMIZE : GRB_MAXIMIZE);

		// Optimize model
        Logger::out("-") << "using the GUROBI solver (version " << GRB_VERSION_MAJOR << "." << GRB_VERSION_MINOR << ")." << std::endl;
		model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        switch (status) {
		case GRB_OPTIMAL: {
			objective_value_ = model.get(GRB_DoubleAttr_ObjVal);
			result_.resize(variables.size());
			for (std::size_t i = 0; i < variables.size(); ++i) {
				result_[i] = X[i].get(GRB_DoubleAttr_X);
			}
			upload_solution(program);
			break;
		}
		
		case GRB_INF_OR_UNBD:
			std::cerr << "model is infeasible or unbounded" << std::endl;
			break;

		case GRB_INFEASIBLE:
			std::cerr << "model is infeasible" << std::endl;
			break;

		case GRB_UNBOUNDED:
			std::cerr << "model is unbounded" << std::endl;
			break;

		default:
			std::cerr << "optimization was stopped with status = " << status << std::endl;
			break;
		}

		return (status == GRB_OPTIMAL);
	}
	catch (GRBException e) {
        Logger::err("-") << e.getMessage() << " (error code: " << e.getErrorCode() << ")." << std::endl;
        if (e.getErrorCode() == GRB_ERROR_NO_LICENSE) {
            Logger::warn("-") << "Gurobi installed but license is missing or expired. Please choose another solver, e.g., SCIP." << std::endl;
        }
	}
	catch (...) {
		std::cerr << "Exception during optimization" << std::endl;
	}

	return false;
}

#endif
