#include "binary_program.h"
#include "../basic/logger.h"
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


bool BinaryProgram::solve(Solver solver /* = GUROBI */) const {
	try {

		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);

		GRBModel model = GRBModel(env);

		// create variables
		std::vector<GRBVar> X(num_variables_);
		for (std::size_t i = 0; i < num_variables_; ++i) {
			X[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		// Integrate new variables
		model.update();

		// Set objective
		GRBLinExpr obj;

		const std::unordered_map<std::size_t, double>& obj_coeffs = objective_.coefficients();
		std::unordered_map<std::size_t, double>::const_iterator it = obj_coeffs.begin();
		for (; it != obj_coeffs.end(); ++it) {
			std::size_t var_idx = it->first;
			double coeff = it->second;
			obj += coeff * X[var_idx];
		}
		model.setObjective(obj, GRB_MINIMIZE);

		// Add constraints
		for (std::size_t i = 0; i < constraints_.size(); ++i) {
			GRBLinExpr expr;
			const Constraint& cstr = constraints_[i];
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
				model.addConstr(expr == cstr.get_fixed_bound());
				break;
			case Constraint::LOWER:
				model.addConstr(expr >= cstr.get_lower_bound());
				break;
			case Constraint::UPPER:
				model.addConstr(expr <= cstr.get_upper_bound());
				break;
			case Constraint::DOUBLE: {
				double lb, ub;
				cstr.get_double_bound(lb, ub);
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

			std::vector<int>& output = const_cast<BinaryProgram*>(this)->result_;
			output.resize(num_variables_);
			for (std::size_t i=0; i<num_variables_; ++i) {
				output[i] = static_cast<int>(X[i].get(GRB_DoubleAttr_X));
			}
			return true;
		}
		else if (optimstatus == GRB_INF_OR_UNBD) {
			Logger::err("-") << "model is infeasible or unbounded" << std::endl;
			return false;
		}
		else if (optimstatus == GRB_INFEASIBLE) {
			Logger::err("-") << "model is infeasible" << std::endl;
			return false;
		}
		else if (optimstatus == GRB_UNBOUNDED) {
			Logger::err("-") << "model is unbounded" << std::endl;
			return false;
		}
		else {
			Logger::err("-") << "optimization was stopped with status = " << optimstatus << std::endl;
			return false;
		}
	}
	catch (GRBException e) {
		Logger::err("-") << "Error code = " << e.getErrorCode() << std::endl;
		Logger::err("-") << e.getMessage() << std::endl;
	}
	catch (...) {
		Logger::err("-") << "Exception during optimization" << std::endl;
	}

	return false;
}

