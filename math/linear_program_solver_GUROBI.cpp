#include "linear_program_solver.h"
#include "../basic/basic_types.h"

#include <gurobi_c++.h>

#ifdef WIN32
#if (_MSC_VER == 1800)  // vs2013
#pragma comment(lib, "gurobi70.lib")
#ifdef _DEBUG
#pragma comment(lib, "gurobi_c++mdd2013.lib")
#else
#pragma comment(lib, "gurobi_c++md2013.lib")
#endif
#elif (_MSC_VER == 1900) // vs 2015
#pragma comment(lib, "gurobi75.lib")
#ifdef _DEBUG
#pragma comment(lib, "gurobi_c++mdd2015.lib")
#else
#pragma comment(lib, "gurobi_c++md2015.lib")
#endif
#elif (_MSC_VER >= 1910 && _MSC_VER <= 1913) // vs 2017
#pragma comment(lib, "gurobi75.lib")
#ifdef _DEBUG
#pragma comment(lib, "gurobi_c++mdd2017.lib")
#else
#pragma comment(lib, "gurobi_c++md2017.lib")
#endif
#endif
#endif


bool LinearProgramSolver::_solve_GUROBI(const LinearProgram* program) {
	try {
		typedef Variable<double>			Variable;
		typedef LinearExpression<double>	Objective;
		typedef LinearConstraint<double>	Constraint;

		const std::vector<Variable>& variables = program->variables();
		if (variables.empty()) {
			std::cerr << "variable set is empty" << std::endl;
			return false;
		}

		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);

		GRBModel model = GRBModel(env);

		// create variables
		std::vector<GRBVar> X(variables.size());
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable& var = variables[i];

			char vtype = GRB_CONTINUOUS;
			if (var.variable_type() == Variable::INTEGER)
				vtype = GRB_INTEGER;
			else if (var.variable_type() == Variable::BINARY)
				vtype = GRB_BINARY;

			double lb = -DBL_MAX;
			double ub = DBL_MAX;
			if (var.bound_type() == Variable::DOUBLE)
				var.get_bound(lb, ub);
			X[i] = model.addVar(lb, ub, 0.0, vtype);
		}

		// Integrate new variables
		model.update();

		// Set objective
		GRBLinExpr obj;

		const Objective& objective = program->objective();
		const std::unordered_map<std::size_t, double>& obj_coeffs = objective.coefficients();
		std::unordered_map<std::size_t, double>::const_iterator it = obj_coeffs.begin();
		for (; it != obj_coeffs.end(); ++it) {
			std::size_t var_idx = it->first;
			double coeff = it->second;
			obj += coeff * X[var_idx];
		}
		model.setObjective(obj, GRB_MINIMIZE);

		// Add constraints
		const std::vector<Constraint>& constraints = program->constraints();
		for (std::size_t i = 0; i < constraints.size(); ++i) {
			GRBLinExpr expr;
			const Constraint& cstr = constraints[i];
			const std::unordered_map<std::size_t, double>& cstr_coeffs = cstr.coefficients();
			std::unordered_map<std::size_t, double>::const_iterator cur = cstr_coeffs.begin();
			for (; cur != cstr_coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first;
				double coeff = cur->second;
				expr += coeff * X[var_idx];
			}

			switch (cstr.bound_type())
			{
			case Constraint::FIXED:
				model.addConstr(expr == cstr.get_bound());
				break;
			case Constraint::LOWER:
				model.addConstr(expr >= cstr.get_bound());
				break;
			case Constraint::UPPER:
				model.addConstr(expr <= cstr.get_bound());
				break;
			case Constraint::DOUBLE: {
				double lb, ub;
				cstr.get_bound(lb, ub);
				model.addConstr(expr >= lb);
				model.addConstr(expr <= ub);
				break;
				}
			default:
				break;
			}
		}

		// Optimize model
		model.optimize();

		int optimstatus = model.get(GRB_IntAttr_Status);
		if (optimstatus == GRB_OPTIMAL) {
			double objval = model.get(GRB_DoubleAttr_ObjVal);
			//Logger::out("-") << "done. objective: " << truncate_digits(objval, 3) << std::endl;

			result_.resize(variables.size());
			for (std::size_t i=0; i<variables.size(); ++i) {
				result_[i] = X[i].get(GRB_DoubleAttr_X);
			}
		}
		else if (optimstatus == GRB_INF_OR_UNBD) 
			std::cerr << "model is infeasible or unbounded" << std::endl;
		else if (optimstatus == GRB_INFEASIBLE) 
			std::cerr << "model is infeasible" << std::endl;
		else if (optimstatus == GRB_UNBOUNDED) 
			std::cerr << "model is unbounded" << std::endl;
		else 
			std::cerr << "optimization was stopped with status = " << optimstatus << std::endl;

		return (optimstatus == GRB_OPTIMAL);
	}
	catch (GRBException e) {
		std::cerr << "Error code = " << e.getErrorCode() << std::endl;
		std::cerr << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cerr << "Exception during optimization" << std::endl;
	}

	return false;
}

