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

#include <linear_program_solver/linear_program.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>    // for std::round()


namespace details {
	/**
	* Converts an integer v to a string of specified 'width' by
	* filling with character 'fill'
	*/
	template <class Int>
	inline std::string from_integer(Int v, int width, char fill) {
		std::ostringstream string_stream;
		string_stream << std::setfill(fill) << std::setw(width) << v;
		return string_stream.str();
	}
}

//double Bounded::infinity_ = std::numeric_limits<double>::max();
double Bound::infinity_ = 1e20;		// in SCIP, values larger than 1e20 are considered infinity


Bound::Bound(double lb /* = -infinity() */, double ub /* = +infinity() */)
	: lower_bound_(lb)
	, upper_bound_(ub)
{
}

double Bound::infinity() {
	return infinity_;
}

Bound::BoundType Bound::bound_type() const {
	if (lower_bound_ <= -infinity_ && upper_bound_ >= infinity_) // free variable
		return FREE;

	else if (lower_bound_ > -infinity_ && upper_bound_ >= infinity_)
		return LOWER;                                                                                                                                                

	else if (lower_bound_ <= -infinity_ && upper_bound_ < infinity_) 
		return UPPER;

	else { // lower_bound_ > -infinity_ && upper_bound_ < infinity_
		if (lower_bound_ == upper_bound_)
			return FIXED;
		else
			return DOUBLE;
	}
}


void Bound::set_bounds(double lb, double ub) {
	lower_bound_ = lb;
	upper_bound_ = ub;
}

void Bound::get_bounds(double& lb, double& ub) const {
	lb = lower_bound_;
	ub = upper_bound_;
}

//////////////////////////////////////////////////////////////////////////


Variable::Variable(
	LinearProgram* program, 
	VariableType type /* = CONTINUOUS */, 
	double lb /* = -infinity() */, 
	double ub /* = +infinity() */, 
	const std::string& name /* = "" */,
	int idx /* = 0*/
)
	: ProgramEntry(program, name, idx)
	, Bound(lb, ub)
	, variable_type_(type)
	, solution_value_(0.0)
{
	if (type == BINARY)
		Bound::set_bounds(0.0, 1.0);
}


void Variable::set_variable_type(VariableType type) {
	variable_type_ = type;
	if (type == BINARY)
		Bound::set_bounds(0.0, 1.0);
}


double Variable::solution_value(bool rounded /* = false*/) const { 
	if (rounded && variable_type_ != CONTINUOUS)
		return std::round(solution_value_);
	else
		return solution_value_;
}


//////////////////////////////////////////////////////////////////////////

LinearExpression::LinearExpression(LinearProgram* program, const std::string& name, int idx)
	: ProgramEntry(program, name, idx)
	, offset_(0.0)
{
}


void LinearExpression::add_coefficient(const Variable* var, double coeff) {
	if (!program()->has_variable(var)) {
		std::cerr << "program does not own variable " << var->name() << " (index " << var->index() << ")" << std::endl;
		return;
	}

	if (coefficients_.find(var) == coefficients_.end())
		coefficients_[var] = coeff;
	else
		coefficients_[var] += coeff;
}


double LinearExpression::get_coefficient(const Variable* var) const {
	if (!program()->has_variable(var)) {
		std::cerr << "program does not own variable " << var->name() << " (index " << var->index() << ")" << std::endl;
		return 0.0;
	}

	std::unordered_map<const Variable*, double>::const_iterator pos = coefficients_.find(var);
	if (pos != coefficients_.end())
		return pos->second;
	else {
		std::cerr << "linear expression does not own variable " << var->name() << " (index " << var->index() << ")" << std::endl;
		return 0.0;
	}
}


double LinearExpression::solution_value(bool rounded /* = false*/) const {
	double solution = offset_;

	std::unordered_map<const Variable*, double>::const_iterator it = coefficients_.begin();
	for (; it != coefficients_.end(); ++it) {
		const Variable* var = it->first;
		double coeff = it->second;
		solution += var->solution_value(rounded) * coeff;
	}
	return solution;
}


//////////////////////////////////////////////////////////////////////////


void LinearObjective::clear() {
	LinearExpression::clear();
	set_name("");
	set_index(0);
	set_offset(0.0);
}


//////////////////////////////////////////////////////////////////////////


LinearConstraint::LinearConstraint(
	LinearProgram* program, 
	double lb /* = -infinity() */, 
	double ub /* = +infinity() */, 
	const std::string& name/* = "" */,
	int idx /* = 0*/
)
	: LinearExpression(program, name, idx)
	, Bound(lb, ub) 
{
}


//////////////////////////////////////////////////////////////////////////


LinearObjective::LinearObjective(LinearProgram* program, Sense sense) 
	: LinearExpression(program)
	, sense_(sense) 
{
}

//////////////////////////////////////////////////////////////////////////


LinearProgram::LinearProgram(const std::string& name /* = "no_name"*/)
	: name_(name)
{
	// intentionally set the objective to UNDEFINED, so it will allow me to warn
	// the user if he/she forgot to set the objective sense.
	objective_ = new LinearObjective(this, LinearObjective::UNDEFINED);
}


LinearProgram::~LinearProgram() 
{
	clear();

	delete objective_;
}


void LinearProgram::clear() {
	for (std::size_t i = 0; i < variables_.size(); ++i)
		delete variables_[i];
	variables_.clear();

	for (std::size_t i = 0; i < constraints_.size(); ++i)
		delete constraints_[i];
	constraints_.clear();

	objective_->clear();
}


Variable* LinearProgram::create_variable(
	Variable::VariableType type /* = Variable::CONTINUOUS */, 
	double lb /* = -Variable::infinity() */, 
	double ub /* = +Variable::infinity() */, 
	const std::string& name /* = "" */)
{
	Variable* v = new Variable(this, type, lb, ub);

	std::size_t idx = variables_.size();
	v->set_index(idx);

	const std::string& fixed_name = name.empty() ? "x" + details::from_integer(idx, 9, '0') : name;
	v->set_name(fixed_name);

	variables_.push_back(v);
	return v;
}


std::vector<Variable*> LinearProgram::create_n_variables(std::size_t n) {
	std::vector<Variable*> variables;
	for (std::size_t i = 0; i < n; ++i) {
		Variable* v = create_variable();
		variables.push_back(v);
	}
	return variables;
}


LinearConstraint* LinearProgram::create_constraint( 
	double lb /* = -Variable::infinity() */, 
	double ub /* = +Variable::infinity() */, 
	const std::string& name /* = "" */ )
{
	LinearConstraint* c = new LinearConstraint(this, lb, ub);

	std::size_t idx = constraints_.size();
	c->set_index(idx);

	const std::string& fixed_name = name.empty() ? "c" + details::from_integer(idx, 9, '0') : name;
	c->set_name(fixed_name);

	constraints_.push_back(c);
	return c;
}


std::vector<LinearConstraint*> LinearProgram::create_n_constraints(std::size_t n) {
	std::vector<LinearConstraint*> constraints;
	for (std::size_t i = 0; i < n; ++i) {
		LinearConstraint* v = create_constraint();
		constraints.push_back(v);
	}
	return constraints;
}


LinearObjective* LinearProgram::create_objective(LinearObjective::Sense sense /* = LinearObjective::MINIMIZE*/) {
	if (objective_)
		delete objective_;

	objective_ = new LinearObjective(this, sense);
	return objective_;
}


bool LinearProgram::has_variable(const Variable* var) const {
	if (var == nullptr) 
		return false;

	if (var->index() >= 0 && var->index() < variables_.size()) {
		// Then, verify that the variable with this index has the same address.
		return variables_[var->index()] == var;
	}
	return false;
}


bool LinearProgram::has_constraint(const LinearConstraint* cons) const {
	if (cons == nullptr)
		return false;

	if (cons->index() >= 0 && cons->index() < constraints_.size()) {
		// Then, verify that the constraint with this index has the same address.
		return constraints_[cons->index()] == cons;
	}
	return false;
}


const LinearObjective* LinearProgram::objective() const { 
	if (!objective_)
		std::cerr << "please call \'create_objective()\' to create an objective first" << std::endl;

	return objective_; 
}


LinearObjective* LinearProgram::objective() { 
	if (!objective_)
		std::cerr << "please call \'create_objective()\' to create an objective first" << std::endl;

	return objective_; 
}


std::size_t LinearProgram::num_continuous_variables() const {
	std::size_t num_continuous_var = 0;
	for (std::size_t i = 0; i < variables_.size(); ++i) {
		const Variable* v = variables_[i];
		if (v->variable_type() == Variable::CONTINUOUS)
			++num_continuous_var;
	}
	return num_continuous_var;
}


std::size_t LinearProgram::num_integer_variables() const {
	std::size_t num_iteger_var = 0;
	for (std::size_t i = 0; i < variables_.size(); ++i) {
		const Variable* v = variables_[i];
		if (v->variable_type() == Variable::INTEGER/* || v->variable_type() == Variable::BINARY*/)
			++num_iteger_var;
	}
	return num_iteger_var;
}


std::size_t LinearProgram::num_binary_variables() const {
	std::size_t num_binary_var = 0;
	for (std::size_t i = 0; i < variables_.size(); ++i) {
		const Variable* v = variables_[i];
		if (v->variable_type() == Variable::BINARY)
			++num_binary_var;
	}
	return num_binary_var;
}


// returns true if all variables are continuous
bool LinearProgram::is_continuous() const {
	std::size_t num = num_continuous_variables();
	return (num > 0) && (num == variables_.size());
}


// returns true if mixed inter program
bool LinearProgram::is_mix_integer_program() const {
	std::size_t num = num_continuous_variables();
	return (num > 0) && (num < variables_.size());
}


// returns true if inter program
bool LinearProgram::is_integer_program() const {
	std::size_t num = num_integer_variables();
	return (num > 0) && (num == variables_.size()); // or (num_continuous_variables() == 0)
}


// returns true if binary program
bool LinearProgram::is_binary_proram() const {
	std::size_t num = num_binary_variables();
	return (num > 0) && (num == variables_.size());
}


// print statistics of the program
void LinearProgram::print_statistics(std::ostream& output) const {
	output << "**************************************************************" << std::endl;
	output << "\\\tProgram name: " << name() << std::endl;
	output << "\\\tProgram type: ";
	if (is_binary_proram())
		output << "Binary" << std::endl;
	else if (is_integer_program())
		output << "Integer" << std::endl;
	else if (is_mix_integer_program())
		output << "Mixed Integer" << std::endl;
	else
		output << "Continuous" << std::endl;
	output << "\\\t#Variables: " << num_variables() << std::endl;
	output << "\\\t\t-Binary: " << num_binary_variables() << std::endl;
	output << "\\\t\t-Integer: " << num_integer_variables() << std::endl;
	output << "\\\t\t-Continuous: " << num_continuous_variables() << std::endl;
	output << "\\\t#Constraints: " << num_constraints() << std::endl;
	output << "\\\t#Objective sense: " << ((objective()->sense() == LinearObjective::MINIMIZE) ? "Minimize" : "Maximize") << std::endl;
	output << "\\\t#Objective offset: " << objective()->offset() << std::endl;
	output << "**************************************************************" << std::endl;
}


bool LinearProgram::is_valid(std::ostream* output) const {
	bool valid = true;

	if (objective()->sense() == LinearObjective::UNDEFINED) {
		valid = false;
		if (output)
			*output << "incomplete objective: undefined objective sense." << std::endl;
	}

	if (variables_.empty()) {
		valid = false;
		if (output)
			*output << "variable set is empty" << std::endl;
	}

	// check if multiple variables have the same name or index
	std::unordered_map<int, const Variable*> vindices;
	std::unordered_map<std::string, const Variable*> vnames;
	for (std::size_t i = 0; i < variables_.size(); ++i) {
		const Variable* v = variables_[i];

		const std::string& name = v->name();
		int idx = v->index();

		if (v->program() != this) {
			valid = false;
			if (output)
				*output << "variable " << v->name() << " (index " << idx << ") is not owned by this program" << std::endl;
		}

		if (idx != i) {
			valid = false;
			if (output)
				*output << "variable " << v->name() << " (index " << idx << ") has an inconsistent index" << std::endl;
		}

		if (v->lower_bound() > v->upper_bound()) {
			valid = false;
			if (output)
				*output << "variable " << v->name() << " (index " << idx << ") has contradictory bounds (i.e. lb > ub): "
				<< " lb = " << v->lower_bound() << ", ub = " << v->upper_bound() << std::endl;
		}

		std::unordered_map<int, const Variable*>::const_iterator vpos_id = vindices.find(idx);
		if (vpos_id == vindices.end())
			vindices[idx] = v;
		else {
			valid = false;
			if (output) {
				if (vpos_id->second == v)
					*output << "duplicated variable " << name << " (index " << idx << ")" << std::endl;
				else
					*output << "variables " << vpos_id->second->name() << " and " << v->name() << " have the same index " << idx << std::endl;
			}
		}

		std::unordered_map<std::string, const Variable*>::const_iterator vpos_name = vnames.find(name);
		if (vpos_name == vnames.end())
			vnames[name] = v;
		else {
			valid = false;
			if (output) {
				if (vpos_name->second == v)
					*output << "duplicated variable " << name << " (index " << idx << ")" << std::endl;
				else
					*output << "variables " << vpos_name->second->index() << " and " << v->index() << " have the same name: " << vpos_name->first << std::endl;
			}
		}
	}

	// check if multiple constraints have the same name or index
	std::unordered_map<int, const LinearConstraint*> cindices;
	std::unordered_map<std::string, const LinearConstraint*> cnames;
	for (std::size_t i = 0; i < constraints_.size(); ++i) {
		const LinearConstraint* c = constraints_[i];
		const std::string& name = c->name();
		int idx = c->index();

		if (c->program() != this) {
			valid = false;
			if (output)
				*output << "constraint " << c->name() << " (index " << idx << ") is not owned by this program" << std::endl;
		}

		if (idx != i) {
			valid = false;
			if (output)
				*output << "constraint indices are not consistent" << std::endl;
		}

		if (c->lower_bound() > c->upper_bound()) {
			valid = false;
			if (output)
				*output << "constraint " << c->name() << " (index " << idx << ") has contradictory bounds (i.e. lb > ub): "
				<< " lb = " << c->lower_bound() << ", ub = " << c->upper_bound() << std::endl;
		}

		std::unordered_map<int, const LinearConstraint*>::const_iterator cpos_id = cindices.find(idx);
		if (cpos_id == cindices.end())
			cindices[idx] = c;
		else {
			valid = false;
			if (output) {
				if (cpos_id->second == c)
					*output << "duplicated constraint " << name << " (index " << idx << ")" << std::endl;
				else
					*output << "constraints " << cpos_id->second->name() << " and " << name << " have the same index " << idx << std::endl;
			}
		}

		std::unordered_map<std::string, const LinearConstraint*>::const_iterator cpos_name = cnames.find(name);
		if (cpos_name == cnames.end())
			cnames[name] = c;
		else {
			valid = false;
			if (output) {
				if (cpos_name->second == c)
					*output << "duplicated constraint " << name << " (index " << idx << ")" << std::endl;
				else
					*output << "constraints " << cpos_name->second->index() << " and " << c->index() << " have the same name: " << cpos_name->first << std::endl;
			}
		}
	}

	return valid;
}
