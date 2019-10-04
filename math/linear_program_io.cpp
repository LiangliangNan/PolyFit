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


#include "linear_program.h"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "scip.h"
#include "scip/scipdefplugins.h"


// sort variables in an increasing order according to their indices
#define SORT_VARIABLES_IN_INCREASING_ORDER

namespace details {

	static const char * const PATH_SEPARATORS = "/\\";

    bool is_file(const std::string& file_name) {
        std::ifstream input(file_name.c_str());
        if (input.fail())
            return false;
		return true;
	}

	std::string extension(const std::string& file_name) {
		std::string::size_type dot = file_name.find_last_of('.');
		std::string::size_type slash = file_name.find_last_of(PATH_SEPARATORS);
		if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
			return std::string("");

		return std::string(file_name.begin() + dot + 1, file_name.end());
	}

	std::string name_without_path(const std::string& file_name) {
		std::string::size_type slash = file_name.find_last_of(PATH_SEPARATORS);
		if (slash == std::string::npos)
			return file_name;
		else
			return std::string(file_name.begin() + slash + 1, file_name.end());
	}

	// strip one level of extension from the filename.
	std::string name_less_extension(const std::string& file_name)
	{
		std::string::size_type dot = file_name.find_last_of('.');
		std::string::size_type slash = file_name.find_last_of(PATH_SEPARATORS);        // Finds forward slash *or* back slash
		if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
			return file_name;

		return std::string(file_name.begin(), file_name.begin() + dot);
	}

	std::string base_name(const std::string& file_path) {
		std::string name = name_without_path(file_path);
		return name_less_extension(name);
	}

#ifdef SORT_VARIABLES_IN_INCREASING_ORDER

	std::vector<std::pair<int, double>> ordered_coefficients(const LinearExpression* expr, bool increasing_order /* = true*/) {
		std::vector<std::pair<int, double>> coeffs(expr->coefficients().begin(), expr->coefficients().end());

		// it is possible to make the order an parameter when constructing SortObj, but this is a bit more efficient 
		if (increasing_order) {
			struct SortIncreasing {
				bool operator()(const std::pair<int, double>& a, const std::pair<int, double>& b) {
					return a.first < b.first;
				}
			};
			std::sort(coeffs.begin(), coeffs.end(), SortIncreasing());
		}
		else {
			struct SortDecreasing {
				bool operator()(const std::pair<int, double>& a, const std::pair<int, double>& b) {
					return a.first > b.first;
				}
			};
			std::sort(coeffs.begin(), coeffs.end(), SortDecreasing());
		}

		return coeffs;
	}
#endif

}


bool LinearProgram::save(const std::string& file_name, bool simple_name /* = false*/) const {
	std::ofstream output(file_name.c_str());
	if (output.fail()) {
		std::cerr << "could not create/open file to save:\'" << file_name << "\'" << std::endl;
		output.close();
		return false;
	}

	const std::string& ext = details::extension(file_name);
	if (ext != "lp" && ext != "mps" && ext != "cip") {
		std::cerr << "unsupported format: \'" << ext << "\'" << std::endl;
		return false;
	}

	Scip* scip = 0;
	SCIP_CALL(SCIPcreate(&scip));
	SCIP_CALL(SCIPincludeDefaultPlugins(scip));

	// disable scip output to stdout
	SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scip), TRUE);

	// create empty problem 
	SCIP_CALL(SCIPcreateProbBasic(scip, name_.c_str()));

	// set the objective sense to maximize, default is minimize
//		SCIP_CALL(SCIPfreeTransform(scip));
	SCIP_CALL(SCIPsetObjsense(scip, objective_->sense() == LinearObjective::MINIMIZE ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE));

	// create variables
	std::vector<SCIP_VAR*> scip_variables;
	for (std::size_t i = 0; i < variables_.size(); ++i) {
		const Variable* var = variables_[i];

		std::string name = var->name();
		if (simple_name)
			name = "x" + std::to_string(i);

		double lb, ub;
		var->get_bounds(lb, ub);

		SCIP_VAR* v = 0;
		//			SCIP_CALL(SCIPfreeTransform(scip));
		switch (var->variable_type())
		{
		case Variable::CONTINUOUS:
			SCIP_CALL(SCIPcreateVar(scip, &v, name.c_str(), lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0));
			break;
		case Variable::INTEGER:
			SCIP_CALL(SCIPcreateVar(scip, &v, name.c_str(), lb, ub, 0.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE, 0, 0, 0, 0, 0));
			break;
		case Variable::BINARY:
			SCIP_CALL(SCIPcreateVar(scip, &v, name.c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, 0, 0, 0, 0, 0));
			break;
		}
		// add the SCIP_VAR object to the scip problem
		SCIP_CALL(SCIPaddVar(scip, v));

		// storing the SCIP_VAR pointer for later access
		scip_variables.push_back(v);
	}

	// Add constraints

	std::vector<SCIP_CONS*> scip_constraints;
	for (std::size_t i = 0; i < constraints_.size(); ++i) {
		const LinearConstraint* c = constraints_[i];
#ifdef SORT_VARIABLES_IN_INCREASING_ORDER
		const std::vector< std::pair<int, double> >& coeffs = details::ordered_coefficients(c, true);
		std::vector< std::pair<int, double> >::const_iterator cur = coeffs.begin();
#else
		const std::unordered_map<const Variable*, double>& coeffs = c->coefficients();
		std::unordered_map<const Variable*, double>::const_iterator cur = coeffs.begin();
#endif
		std::vector<SCIP_VAR*>	cstr_variables(coeffs.size());
		std::vector<double>		cstr_values(coeffs.size());
		std::size_t idx = 0;
		for (; cur != coeffs.end(); ++cur) {
			std::size_t var_idx = cur->first;
			double coeff = cur->second;
			cstr_variables[idx] = scip_variables[var_idx];
			cstr_values[idx] = coeff;
			++idx;
		}

		// create SCIP_CONS object
		SCIP_CONS* cons = 0;
		std::string name = c->name();
		if (simple_name)
			name = "c" + std::to_string(i);

		double lb, ub;
		c->get_bounds(lb, ub);

		//			SCIP_CALL(SCIPfreeTransform(scip));
		SCIP_CALL(SCIPcreateConsLinear(scip, &cons, name.c_str(), coeffs.size(), cstr_variables.data(), cstr_values.data(), lb, ub, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
		SCIP_CALL(SCIPaddCons(scip, cons));			// add the constraint to scip

		// store the constraint for later on
		scip_constraints.push_back(cons);
	}

	// set objective

	// determine the coefficient of each variable in the objective function
#ifdef SORT_VARIABLES_IN_INCREASING_ORDER
	const std::vector< std::pair<int, double> >& obj_coeffs = details::ordered_coefficients(objective_, true);
	std::vector< std::pair<int, double> >::const_iterator cur = obj_coeffs.begin();
#else
	const std::unordered_map<const Variable*, double>& obj_coeffs = objective()->coefficients();
	std::unordered_map<const Variable*, double>::const_iterator cur = obj_coeffs.begin();
#endif
	for (; cur != obj_coeffs.end(); ++cur) {
		std::size_t var_idx = cur->first;
		double coeff = cur->second;
		//			SCIP_CALL(SCIPfreeTransform(scip));
		SCIP_CALL(SCIPchgVarObj(scip, scip_variables[var_idx], coeff));
	}

	//		SCIP_CALL(SCIPfreeTransform(scip));
	SCIP_RETCODE status = SCIPwriteOrigProblem(scip, file_name.c_str(), 0, 0);

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

	return  status == SCIP_OKAY;
}


bool LinearProgram::load(const std::string& file_name) {
	if (!details::is_file(file_name)) {
		std::cerr << "file does not exist: \'" << file_name << "\'" << std::endl;
		return false;
	}

	clear();

	const std::string& ext = details::extension(file_name);
	if (ext != "lp" && ext != "mps" && ext != "cip") {
		std::cerr << "unsupported format: \'" << ext << "\'" << std::endl;
		return false;
	}

	Scip* scip = 0;
	SCIP_CALL(SCIPcreate(&scip));
	SCIP_CALL(SCIPincludeDefaultPlugins(scip));

	// disable scip output to stdout
	SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scip), TRUE);

	SCIP_CALL(SCIPreadProb(scip, file_name.c_str(), 0));

	// SCIPgetProbName() returns the original file name
	name_ = details::base_name(SCIPgetProbName(scip));
	LinearObjective::Sense s = SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? LinearObjective::MINIMIZE : LinearObjective::MAXIMIZE;
	objective_->set_sense(s);

	const double infinity = Bound::infinity();

	// create variables
	int num_var = SCIPgetNVars(scip);
	SCIP_VAR** scip_variables = SCIPgetVars(scip);
	const std::vector<Variable*>& variables = create_n_variables(num_var);
	for (std::size_t i = 0; i < num_var; ++i) {
		SCIP_VAR* v = scip_variables[i];
		const std::size_t idx = SCIPvarGetIndex(v);
		Variable* var = variables[idx];

		const char* name = SCIPvarGetName(v);
		var->set_name(name);

		double lb = SCIPvarGetLbGlobal(v);
		double ub = SCIPvarGetUbGlobal(v);
		var->set_bounds(lb, ub);

		switch (SCIPvarGetType(v))
		{
		case SCIP_VARTYPE_BINARY:
			var->set_variable_type(Variable::BINARY);
			break;
		case SCIP_VARTYPE_INTEGER:
			var->set_variable_type(Variable::INTEGER);
			break;
		case SCIP_VARTYPE_CONTINUOUS:
		default:
			var->set_variable_type(Variable::CONTINUOUS);
			break;
		}
	}

	// Add constraints

	int num_cons = SCIPgetNConss(scip);
	SCIP_CONS** scip_constraints = SCIPgetConss(scip);
	const std::vector<LinearConstraint*>& constraints = create_n_constraints(num_cons);
	for (std::size_t i = 0; i < num_cons; ++i) {
		SCIP_CONS* cons = scip_constraints[i];
		const int idx = i;// SCIPconsGetPos(cons); // Liangliang: why the idx is -1? 
		LinearConstraint* c = constraints[idx];

		const char* name = SCIPconsGetName(cons);
		c->set_name(name);

		SCIP_VAR** vars = SCIPgetVarsLinear(scip, cons);
		const double* coeffs = SCIPgetValsLinear(scip, cons);

		const int nv = SCIPgetNVarsLinear(scip, cons);
		for (std::size_t j = 0; j < nv; ++j) {
			SCIP_VAR* v = vars[j];
			std::size_t idx = SCIPvarGetIndex(v);
			double coeff = coeffs[j];
			c->add_coefficient(idx, coeff);
		}

		double lb = SCIPgetLhsLinear(scip, cons);
		double ub = SCIPgetRhsLinear(scip, cons);
		c->set_bounds(lb, ub);
	}

	// set objective

	for (std::size_t i = 0; i < num_var; ++i) {
		SCIP_VAR* v = scip_variables[i];
		std::size_t idx = SCIPvarGetIndex(v);
		double coeff = SCIPvarGetObj(v);
		objective_->add_coefficient(idx, coeff);
	}

	// after releasing all vars and cons we can free the scip problem
	// remember this has always to be the last call to scip
	SCIP_CALL(SCIPfree(&scip));

	return  true;
}
