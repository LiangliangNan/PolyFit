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

#include "linear_program.h"
#include "../basic/file_utils.h"

#include <fstream>

bool save(const LinearProgram<double>& program, const std::string& file) {
	std::ofstream output(file.c_str());
	if (output.fail()) {
		std::cerr << "failed to create file \'" << file.c_str() << "\'" << std::endl;
		return false;
	}

	typedef Variable<double>			Variable;
	typedef LinearExpression<double>	Objective;
	typedef LinearConstraint<double>	Constraint;

	const std::string& name = program.name();
	output << "#problem: " << (name.empty() ? FileUtils::base_name(file) : name) << std::endl;
	output << "num_variables: " << program.num_variables() << std::endl;
	output << "num_constraints: " << program.num_constraints() << std::endl;
	output << "optimization_direction: " << program.objective_sense() << std::endl;

	// variables
	output << "#variables: variable_index, variable_type, bound_type, lower_bound, upper_bound" << std::endl;
	const std::vector<Variable>& variables = program.variables();

	for (std::size_t i = 0; i < program.num_variables(); ++i) {
		const Variable& var = variables[i];
		double lb, ub;
		var.get_double_bounds(lb, ub);
		output << "variable: " << i << " " << var.variable_type() << " " << var.bound_type() << " " << lb << " " << ub << std::endl;
	}

	// constraints
	output << "#constraints: constraint_index, bound_type, lower_bound, upper_bound, num_coefficients, variable_index, coefficient, ..." << std::endl;
	const std::vector<Constraint>& constraints = program.constraints();
	for (std::size_t i = 0; i < constraints.size(); ++i) {
		const Constraint& cstr = constraints[i];
		const std::unordered_map<std::size_t, double>& coeffs = cstr.coefficients();
		std::unordered_map<std::size_t, double>::const_iterator cur = coeffs.begin();

		double lb, ub;
		cstr.get_double_bounds(lb, ub);
		output << "constraint: " << i << " " << cstr.bound_type() << " " << lb << " " << ub << " " << coeffs.size() << " ";
		for (; cur != coeffs.end(); ++cur)
			output << cur->first << " " << cur->second << " ";
		output << std::endl;
	}

	// objective
	output << "#objective: num_coefficients, variable_index, coefficient, ..." << std::endl;
	const Objective& obj = program.objective();
	const std::unordered_map<std::size_t, double>& obj_coeffs = obj.coefficients();
	std::unordered_map<std::size_t, double>::const_iterator cur = obj_coeffs.begin();
	output << obj_coeffs.size() << std::endl;
	for (; cur != obj_coeffs.end(); ++cur)
		output << cur->first << " " << cur->second << std::endl;

	return true;
}

