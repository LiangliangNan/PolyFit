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
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <cmath>


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


//double Bound::infinity_ = std::numeric_limits<double>::max();
double Bound::infinity_ = 1e20;		// in SCIP, values larger than 1e20 are considered infinity


Bound::Bound(
	BoundType type /* = FREE */,
	double lb /* = -infinity() */, 
	double ub /* = +infinity() */
)
{
	set_bounds(type, lb, ub);
}

double Bound::infinity() {
	return infinity_;
}

void Bound::set_bounds(BoundType type, double lb, double ub) {
	switch (type)
	{
	case FIXED:
		if (std::abs(lb - ub) > 1e-10)
			std::cerr << "lower/upper bounds must be equal for FIXED bound" << std::endl;
		lower_bound_ = upper_bound_ = lb;
		break;
	case LOWER:		
		lower_bound_ = lb;	
		break;
	case UPPER:	
		upper_bound_ = ub;	
		break;
	case DOUBLE:
		lower_bound_ = lb;
		upper_bound_ = ub;
		break;
	case FREE:
		lower_bound_ = -infinity_;	// adjust the bound value
		upper_bound_ = +infinity_;	// adjust the bound value
		break;
	default:
		std::cerr << "fatal error: unrecognized bound type" << std::endl;
		break;
	}
	bound_type_ = type;
}

void Bound::set_bound(BoundType type, double value) 
{
	switch (type)
	{
	case FIXED:
		bound_type_ = FIXED;
		lower_bound_ = upper_bound_ = value;
		break;

	case LOWER:	
		bound_type_ = LOWER;
		lower_bound_ = value;
		break;
	case UPPER:	
		bound_type_ = UPPER;
		upper_bound_ = value;	
		break;

	case DOUBLE:
		std::cerr << "Warning: please use \'set_bounds()\' to set DOUBLE bound" << std::endl;
		break;

	case FREE:
		lower_bound_ = -infinity_;	// adjust the bound value
		upper_bound_ = +infinity_;	// adjust the bound value
		bound_type_ = FREE;
		break;
	default:
		std::cerr << "fatal error: unrecognized bound type" << std::endl;
		break;
	}
}

// query the single bound according to its bound type
double Bound::get_bound() const {
	switch (bound_type_)
	{
	case FIXED:
	case LOWER:
		return lower_bound_;
	case UPPER:
		return upper_bound_;
	case DOUBLE:
		std::cerr << "warning: please use \'get_bounds()\' to query DOUBLE bound" << std::endl;
		return lower_bound_;
	case FREE:
	default:
		std::cerr << "warning: no bound exists (this is a FREE variable)" << std::endl;
		return lower_bound_;
	}
}

void Bound::set_bounds(double lb, double ub) {
	if (lb <= -infinity_ && ub >= infinity_) { // free variable
		std::cerr << "variable with bounds (" << lb << ", " << ub << ") should be a FREE variable" << std::endl;
		bound_type_ = FREE;
	}
	else if (lb > -infinity_ && ub < infinity_) {
		if (lb == ub) 
			bound_type_ = FIXED;
		else 
			bound_type_ = DOUBLE;
	}
	else if (lb > -infinity_ && ub >= infinity_) 
		bound_type_ = LOWER;
	else if (lb <= -infinity_ && ub < infinity_) 
		bound_type_ = UPPER;
	else
		std::cerr << "fatal error: should not reach here" << std::endl;

	lower_bound_ = lb;
	upper_bound_ = ub;
}

void Bound::get_bounds(double& lb, double& ub) const {
	lb = lower_bound_;
	ub = upper_bound_;
}

//////////////////////////////////////////////////////////////////////////


Variable::Variable(LinearProgram* program, VariableType t /* = CONTINUOUS*/)
	: ProgramElement(program)
	, variable_type_(t)
	, solution_value_(0.0)
{
	if (t == BINARY)
		Bound::set_bounds(0.0, 1.0);
}


void Variable::set_variable_type(VariableType t) {
	variable_type_ = t;
	if (t == BINARY)
		Bound::set_bounds(0.0, 1.0);
}


double Variable::solution_value(bool rounded /* = false*/) const { 
	if (rounded && variable_type_ != CONTINUOUS)
		return std::round(solution_value_);
	else
		return solution_value_;
}


//////////////////////////////////////////////////////////////////////////

LinearExpression::LinearExpression(LinearProgram* program)
	: ProgramElement(program)
{
}


void LinearExpression::add_coefficient(int var_index, double coeff) {
	if (coefficients_.find(var_index) == coefficients_.end())
		coefficients_[var_index] = coeff;
	else
		coefficients_[var_index] += coeff;
}

double LinearExpression::solution_value(bool rounded /* = false*/) const {
	double solution = 0.0;

	const std::vector<Variable*>& variables = program()->variables();
	std::unordered_map<int, double>::const_iterator it = coefficients_.begin();
	for (; it != coefficients_.end(); ++it) {
		int var_index = it->first;
		solution += variables[var_index]->solution_value(rounded) * it->second;
	}
	return solution;
}


//////////////////////////////////////////////////////////////////////////


LinearConstraint::LinearConstraint(LinearProgram* program, LinearConstraint::BoundType bt, double lb, double ub) 
	: LinearExpression(program)
	, Bound(bt, lb, ub) 
{
}


//////////////////////////////////////////////////////////////////////////


LinearObjective::LinearObjective(LinearProgram* program, Sense s) 
	: LinearExpression(program)
	, sense_(s) 
{
}

//////////////////////////////////////////////////////////////////////////


LinearProgram::LinearProgram()
	: name_("unknown")
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
	Variable::VariableType vt /* = Variable::CONTINUOUS */, 
	Variable::BoundType bt /* = Variable::FREE */, 
	double lb /* = -Variable::infinity() */, 
	double ub /* = +Variable::infinity() */, 
	const std::string& name /* = "" */)
{
	Variable* v = new Variable(this, vt);
	v->set_bounds(bt, lb, ub);

	std::size_t idx = variables_.size();
	v->set_index(idx);

	const std::string& fixed_name = name.empty() ? "x" + from_integer(idx, 9, '0') : name;
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
	LinearConstraint::BoundType bt /* = Variable::FREE */, 
	double lb /* = -Variable::infinity() */, 
	double ub /* = +Variable::infinity() */, 
	const std::string& name /* = "" */ )
{
	LinearConstraint* c = new LinearConstraint(this, bt, lb, ub);

	std::size_t idx = constraints_.size();
	c->set_index(idx);

	const std::string& fixed_name = name.empty() ? "c" + from_integer(idx, 9, '0') : name;
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
		if (v->variable_type() == Variable::INTEGER || v->variable_type() == Variable::BINARY)
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
