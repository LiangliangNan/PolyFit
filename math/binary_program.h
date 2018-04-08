#pragma once

#include "math_common.h"
#include <vector>
#include <unordered_map>
#include <cassert>


template < class FT>
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
class LinearConstraint : public LinearExpression<FT> 
{
public:
	enum BoundType {
		FIXED,  // constant value
		LOWER,
		UPPER,
		DOUBLE,	// both lower and upper bounds
		UNKNOWN
	};

public:
	LinearConstraint() 
		: bound_type_(UNKNOWN)
	{
	}

	// once set, bound type cannot be changed
	void set_fixed_bound(FT bd) { assert(bound_type_ == UNKNOWN);	bound_type_ = FIXED;	bound_ = bd; }
	void set_lower_bound(FT bd) { assert(bound_type_ == UNKNOWN);	bound_type_ = LOWER;	bound_ = bd; }
	void set_upper_bound(FT bd) { assert(bound_type_ == UNKNOWN);	bound_type_ = UPPER;	bound2_ = bd; }
	void set_double_bound(FT lb, FT ub) { assert(bound_type_ == UNKNOWN);	bound_type_ = DOUBLE;	bound_ = lb; bound2_ = ub; }

	BoundType bound_type() const { return bound_type_; }

	// you can only call the one corresponding to the correct bound type
	FT get_fixed_bound() const { assert(bound_type_ == FIXED);	return bound_; }
	FT get_lower_bound() const { assert(bound_type_ == LOWER);	return bound_; }
	FT get_upper_bound() const { assert(bound_type_ == UPPER);	return bound_; }
	void get_double_bound(FT& lb, FT& ub) const { assert(bound_type_ == DOUBLE);	lb = bound_; ub = bound2_; }

private:
	BoundType	bound_type_;
	FT			bound_;
	FT			bound2_;
};


class MATH_API BinaryProgram
{
public:
	typedef LinearExpression<double>	Objective;
	typedef LinearConstraint<double>	Constraint;
	typedef std::vector<Constraint>		Constraints;

	enum Solver { GUROBI, SCIP, LPSOLVE, GLPK };

public:
	BinaryProgram() {}
	~BinaryProgram() {}

	// num of binary variables
	std::size_t num_variables() const { 
		return num_variables_; 
	}

	void set_objective(const Objective& obj, std::size_t num_var) { 
		objective_ = obj; 
		num_variables_ = num_var;
	}

	const Objective& objective() const {
		return objective_;
	}

	void add_constraint(const Constraint& cstr) { 
		assert(cstr.bound_type() != Constraint::UNKNOWN); 
		constraints_.push_back(cstr); 
	}

	const Constraints& constraints() const { 
		return constraints_; 
	}

	// Solve the problem; returns false if fails
	// NOTE: Gurobi solver recommended;
	//		 The SCIP solver is about x10 slower than Gurobi;
	//       LPSOLVE and GLPK may be too slow or even fail.
	bool solve(Solver solver = GUROBI) const;

	// NOTE:result is valid only if the solver succeeded
	const std::vector<int>& get_result() const { 
		return result_; 
	}

private:
	std::size_t	num_variables_;
	Objective	objective_;
	Constraints	constraints_;

	std::vector<int> result_;
};

