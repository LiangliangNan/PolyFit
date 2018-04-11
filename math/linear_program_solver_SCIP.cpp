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
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <iostream>


bool LinearProgramSolver::_solve_SCIP(const LinearProgram* program) {
	try {
		typedef Variable<double>			Variable;
		typedef LinearExpression<double>	Objective;
		typedef LinearConstraint<double>	Constraint;

		const std::vector<Variable>& variables = program->variables();
		if (variables.empty()) {
			std::cerr << "variable set is empty" << std::endl;
			return false;
		}

		Scip* scip = 0;
		SCIP_CALL(SCIPcreate(&scip));
		SCIP_CALL(SCIPincludeDefaultPlugins(scip));

		// disable scip output to stdout
		SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scip), TRUE);

		// use wall clock time because getting CPU user seconds
		// involves calling times() which is very expensive
		SCIP_CALL(SCIPsetIntParam(scip, "timing/clocktype", SCIP_CLOCKTYPE_WALL));

		// create empty problem 
		SCIP_CALL(SCIPcreateProbBasic(scip, "PolyFit"));

		// set the objective sense to maximize, default is minimize
		SCIP_CALL(SCIPfreeTransform(scip));
		SCIP_CALL(SCIPsetObjsense(scip, program->objective_sense() == LinearProgram::MINIMIZE ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE));

		// create variables
		std::vector<SCIP_VAR*> scip_variables;
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable& var = variables[i];
			const std::string& name = "x" + std::to_string(i + 1);
			SCIP_VAR* v = 0;

			double lb, ub;
			var.get_double_bounds(lb, ub);

			SCIP_CALL(SCIPfreeTransform(scip));
			// The true objective coefficient will be set later in ExtractObjective.
			double tmp_obj_coef = 0.0;
			switch (var.variable_type())
			{
			case Variable::CONTINUOUS:
				SCIP_CALL(SCIPcreateVar(scip, &v, name.data(), lb, ub, tmp_obj_coef, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0));
				break;
			case Variable::INTEGER:
				SCIP_CALL(SCIPcreateVar(scip, &v, name.data(), lb, ub, tmp_obj_coef, SCIP_VARTYPE_INTEGER, TRUE, FALSE, 0, 0, 0, 0, 0));
				break;
			case Variable::BINARY:
				SCIP_CALL(SCIPcreateVar(scip, &v, name.data(), 0, 1, tmp_obj_coef, SCIP_VARTYPE_BINARY, TRUE, FALSE, 0, 0, 0, 0, 0));
				break;
			}
			// add the SCIP_VAR object to the scip problem
			SCIP_CALL(SCIPaddVar(scip, v));

			// storing the SCIP_VAR pointer for later access
			scip_variables.push_back(v);
		}

		// Add constraints

		std::vector<SCIP_CONS*> scip_constraints;
		const std::vector<Constraint>& constraints = program->constraints();
		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const Constraint& cstr = constraints[i];
			const std::unordered_map<std::size_t, double>& cstr_coeffs = cstr.coefficients();
			std::unordered_map<std::size_t, double>::const_iterator cur = cstr_coeffs.begin();

			std::vector<SCIP_VAR*>	cstr_variables(cstr_coeffs.size());
			std::vector<double>		cstr_values(cstr_coeffs.size());
			std::size_t idx = 0;
			for (; cur != cstr_coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first;
				double coeff = cur->second;
				cstr_variables[idx] = scip_variables[var_idx];
				cstr_values[idx] = coeff;
				++idx;
			}

			// create SCIP_CONS object
			SCIP_CONS* cons = 0;
			const std::string& name = "cstr" + std::to_string(i + 1);

			double lb, ub;
			cstr.get_double_bounds(lb, ub);

			SCIP_CALL(SCIPfreeTransform(scip));
			SCIP_CALL(SCIPcreateConsLinear(scip, &cons, name.data(), cstr_coeffs.size(), cstr_variables.data(), cstr_values.data(), lb, ub, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
			SCIP_CALL(SCIPaddCons(scip, cons));			// add the constraint to scip

			// store the constraint for later on
			scip_constraints.push_back(cons);
		}

		// set objective

		// determine the coefficient of each variable in the objective function
		const Objective& objective = program->objective();
		const std::unordered_map<std::size_t, double>& obj_coeffs = objective.coefficients();
		std::unordered_map<std::size_t, double>::const_iterator it = obj_coeffs.begin();
		for (; it != obj_coeffs.end(); ++it) {
			std::size_t var_idx = it->first;
			double coeff = it->second;
			SCIP_CALL(SCIPfreeTransform(scip));
			SCIP_CALL(SCIPchgVarObj(scip, scip_variables[var_idx], coeff));
		}

		// set SCIP parameters
		double tolerance = 1e-7;
		SCIP_CALL(SCIPsetRealParam(scip, "numerics/feastol", tolerance));
		SCIP_CALL(SCIPsetRealParam(scip, "numerics/dualfeastol", tolerance));
		SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", -1));  // enable presolve
		double MIP_gap = 1e-4;
		SCIP_CALL(SCIPsetRealParam(scip, "limits/gap", MIP_gap));

		bool status = false;
		// this tells scip to start the solution process
		if (SCIPsolve(scip) == SCIP_OKAY) {
			// get the best found solution from scip
			SCIP_SOL* sol = SCIPgetBestSol(scip);
			if (sol) {
				// If optimal or feasible solution is found.
				objective_value_ = SCIPgetSolOrigObj(scip, sol);
				result_.resize(variables.size());
				for (std::size_t i = 0; i < variables.size(); ++i) {
					result_[i] = SCIPgetSolVal(scip, sol, scip_variables[i]);
				}
				status = true;
			}
		}

		// report the status: optimal, infeasible, etc.
		SCIP_STATUS scip_status = SCIPgetStatus(scip);
		switch (scip_status) {
		case SCIP_STATUS_OPTIMAL:
			// provides info only if fails.
			break;
		case SCIP_STATUS_GAPLIMIT:
			// To be consistent with the other solvers.
			// provides info only if fails.
			break;
		case SCIP_STATUS_INFEASIBLE:
			std::cerr << "model was infeasible" << std::endl;
			break;
		case SCIP_STATUS_UNBOUNDED:
			std::cerr << "model was unbounded" << std::endl;
			break;
		case SCIP_STATUS_INFORUNBD:
			std::cerr << "model was either infeasible or unbounded" << std::endl;
			break;
		default:
			if (scip_status == SCIP_STATUS_TIMELIMIT) 
				std::cerr << "time limit reached" << std::endl;
			break;
		}

		SCIP_CALL(SCIPresetParams(scip));

		// since the SCIPcreateVar captures all variables, we have to release them now
		for (std::size_t i = 0; i < scip_variables.size(); ++i)
			SCIP_CALL(SCIPreleaseVar(scip, &scip_variables[i]));
		scip_variables.clear();

		// the same for the constraints
		for (std::size_t i = 0; i < scip_constraints.size(); ++i) 
			SCIP_CALL(SCIPreleaseCons(scip, &scip_constraints[i]));
		scip_constraints.clear();

		// after releasing all vars and cons we can free the scip problem
	   // remember this has always to be the last call to scip
		SCIP_CALL(SCIPfree(&scip));

		return status;
	}
	catch (std::exception e) {
		std::cerr << "Error code = " << e.what() << std::endl;
	}
	catch (...) {
		std::cerr << "Exception during optimization" << std::endl;
	}
	return false;
}