#pragma once

#include "math_common.h"

#include <vector>
#include <cassert>
#include <unordered_map>
#include <iostream>



template <class FT>
class Bound
{
public:
	// bound type for variables and expression
enum BoundType { FIXED, LOWER, UPPER, DOUBLE, FREE };

public:
	Bound()
		: bound_type_(FREE)
		, bound_(0.0)
		, bound2_(1.0)
	{
	}

	void set_bound(BoundType type, FT bd) {
		bound_type_ = type;
		bound_ = bd;
	}
	void set_bound(FT lb, FT ub) {
		bound_type_ = DOUBLE;
		bound_ = lb;
		bound2_ = ub;
	}

	BoundType bound_type() const { return bound_type_; }

	// you can only call the one corresponding to the correct bound type
	FT get_bound() const {
		if (bound_type_ == FREE)
			std::cerr << "invalid bound(s): not specified" << std::endl;
		else if (bound_type_ == DOUBLE)
			std::cerr << "wrong function call for double bounded expression" << std::endl;
		return bound_;
	}
	void get_bound(FT& lb, FT& ub) const {
		if (bound_type_ != DOUBLE)
			std::cerr << "wrong function call for single bounded expression" << std::endl;
		lb = bound_;
		ub = bound2_;
	}

private:
	BoundType	bound_type_;
	FT			bound_;
	FT			bound2_;	// used only for double bounds
};


template <class FT>
class Variable : public Bound<FT> {
public:
enum VariableType { CONTINUOUS, INTEGER, BINARY };

public:
	Variable(VariableType type) : variable_type_(type) {
		if (type == BINARY)
			set_bound(0.0, 1.0);
	}

	void set_variable_type(VariableType type) { variable_type_ = type; }
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
	typedef Variable<FT>			Variable;
	typedef LinearExpression<FT>	Objective;
	typedef LinearConstraint<FT>	Constraint;

	enum Solver { GUROBI, SCIP, LPSOLVE, GLPK };

public:
	LinearProgram() {}
	~LinearProgram() {}

	// num of binary variables
	void add_variable(const Variable& var) { 
		variables_.push_back(var);
	}
	void add_variables(const std::vector<Variable>& vars) {
		variables_.insert(variables_.end(), vars.begin(), vars.end());
	}
	const std::vector<Variable>& variables() const {
		return variables_;
	}

	void set_objective(const Objective& obj) { 
		objective_ = obj; 
	}
	const Objective& objective() const {
		return objective_;
	}

	void add_constraint(const Constraint& cstr) { 
		if (cstr.bound_type() == Constraint::FREE) {
			std::cerr << "incomplete constraint: no bound(s) specified" << std::endl;
			return;
		}
		constraints_.push_back(cstr); 
	}
	void add_constraints(const std::vector<Constraint>& cstrs) {
		constraints_.insert(constraints_.end(), cstrs.begin(), cstrs.end());
	}
	const std::vector<Constraint>& constraints() const {
		return constraints_; 
	}

private:
	std::vector<Variable>	variables_;
	std::vector<Constraint>	constraints_;

	Objective	objective_;
};

