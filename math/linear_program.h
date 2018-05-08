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

#ifndef _MATH_LINEAR_PROGRAM_H_
#define _MATH_LINEAR_PROGRAM_H_

#include "math_common.h"

#include <vector>
#include <cassert>
#include <unordered_map>
#include <iostream>
#include <cfloat>



template <class FT>
class Bound
{
public:
	// bound type for variables and expressions
enum BoundType { FIXED, LOWER, UPPER, DOUBLE, FREE };

public:
	Bound()
		: bound_type_(FREE)
		, lower_bound_(-DBL_MAX)
		, upper_bound_(DBL_MAX)
	{
	}

	BoundType bound_type() const { return bound_type_; }

	void set_bounds(BoundType type, FT lb, FT ub) {
		if (bound_type_ != FREE) {
			// In general, bound(s) once set should not be changed.
			// Easier life: print a message if you want to change bound(s)
			std::cerr << "Warning: are you sure you want to changed the bound(s)" << std::endl;
		}

		switch (type)
		{
		case FIXED:
			if (std::abs(lb - ub) > 1e-10)
				std::cerr << "lower/upper bounds must be equal for FIXED bounds" << std::endl;
			lower_bound_ = upper_bound_ = lb;	
			break;
		case LOWER:		lower_bound_ = lb;	break;
		case UPPER:		upper_bound_ = ub;	break;
		case DOUBLE:	
			lower_bound_ = lb; 
			upper_bound_ = ub; 
			break;
		case FREE:
		default:
			std::cerr << "no FREE bound(s)" << std::endl;
			break;
		}
		bound_type_ = type;
	}

	// query the single bound according to its bound type
	FT   get_single_bound() const { 
		switch (bound_type_)
		{
		case FIXED:
		case LOWER:
			return lower_bound_;
		case UPPER:
			return upper_bound_;
		case DOUBLE:
			std::cerr << "please use get_double_bounds() to get double bound(s)" << std::endl;
			return lower_bound_;
		case FREE:
		default:
			std::cerr << "no bound specified" << std::endl;
			return lower_bound_;
		}
	}

	void get_double_bounds(FT& lb, FT& ub) const { lb = lower_bound_; ub = upper_bound_; }

private:
	BoundType	bound_type_;
	FT			lower_bound_;
	FT			upper_bound_;	
};


template <class FT>
class Variable : public Bound<FT> {
public:
enum VariableType { CONTINUOUS, INTEGER, BINARY };

public:
	Variable(VariableType type) : variable_type_(type) {
		if (type == BINARY)
			Bound<FT>::set_bounds(Bound<FT>::DOUBLE, 0.0, 1.0);
	}

	VariableType variable_type() const { return variable_type_; }

private:
	VariableType variable_type_;
};


template <class FT>
class LinearExpression {
public:
	LinearExpression() {}

	void add_coefficient(std::size_t var_index, FT coeff) {
		if (coefficients_.find(var_index) == coefficients_.end())
			coefficients_[var_index] = coeff;
		else 
			coefficients_[var_index] += coeff;
	}

	const std::unordered_map<std::size_t, FT>& coefficients() const {
		return coefficients_;
	}

private:
	std::unordered_map<std::size_t, FT>	coefficients_;
};


template <class FT>
class LinearConstraint : public LinearExpression<FT>, public Bound<FT> 
{
public:
	LinearConstraint() {}
};


template <class FT>
class LinearProgram
{
public:
	enum Solver { GUROBI, SCIP, LPSOLVE, GLPK };

	enum Sense  { MINIMIZE, MAXIMIZE, UNDEFINED };

	typedef Variable<FT>			Variable;
	typedef LinearExpression<FT>	Objective;
	typedef LinearConstraint<FT>	Constraint;

public:
	LinearProgram() : objective_sense_(UNDEFINED) {}
	~LinearProgram() {}

	void add_variable(const Variable& var) { variables_.push_back(var);	}
	void add_variables(const std::vector<Variable>& vars) {	variables_.insert(variables_.end(), vars.begin(), vars.end()); }
	const std::vector<Variable>& variables() const { return variables_;	}

	void set_objective(const Objective& obj, Sense sense) {
		objective_ = obj; 
		objective_sense_ = sense;
	}

	const Objective& objective() const { return objective_; }
	Sense objective_sense() const { return objective_sense_; }

	// add an existing constraint to the current problem
	void add_constraint(const Constraint& cstr) { 
		if (cstr.bound_type() == Constraint::FREE) {
			std::cerr << "incomplete constraint: no bound(s) specified. Constraint ignored." << std::endl;
			return;
		}
		constraints_.push_back(cstr); 
	}

	// add a set of existing constraints to the current problem
	void add_constraints(const std::vector<Constraint>& cstrs) {
		for (std::size_t i = 0; i < cstrs.size(); ++i)
			add_constraint(cstrs[i]);
	}
	const std::vector<Constraint>& constraints() const {
		return constraints_; 
	}

	std::size_t num_constraints() const { return constraints_.size(); }

	std::size_t num_continuous_variables() const {
		std::size_t num_continuous_var = 0;
		for (std::size_t i = 0; i < variables_.size(); ++i) {
			if (variables_[i].variable_type() == Variable::CONTINUOUS)
				++num_continuous_var;
		}
		return num_continuous_var;
	}	
	
	std::size_t num_integer_variables() const {
		std::size_t num_iteger_var = 0;
		for (std::size_t i = 0; i < variables_.size(); ++i) {
			if (variables_[i].variable_type() == Variable::INTEGER)
				++num_iteger_var;
		}
		return num_iteger_var;
	}

	std::size_t num_binary_variables() const {
		std::size_t num_binary_var = 0;
		for (std::size_t i = 0; i < variables_.size(); ++i) {
			if (variables_[i].variable_type() == Variable::BINARY)
				++num_binary_var;
		}
		return num_binary_var;
	}

	// returns true if mixed inter program
	bool is_mix_integer_program() const {
		std::size_t num = num_continuous_variables();
		return (num > 0) && (num < variables_.size());
	}

	// returns true if inter program
	bool is_integer_program() const {
		std::size_t num = num_integer_variables();
		return (num > 0) && (num == variables_.size());
	}

	// returns true if binary program
	bool is_binary_proram() const {
		std::size_t num = num_binary_variables();
		return (num > 0) && (num == variables_.size());
	}

private:
	std::vector<Variable>	variables_;
	std::vector<Constraint>	constraints_;

	Objective	objective_;
	Sense		objective_sense_;
};


#endif
