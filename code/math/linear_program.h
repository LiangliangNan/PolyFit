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

#ifndef _MATH_LINEAR_PROGRAM_H_
#define _MATH_LINEAR_PROGRAM_H_

#include "math_common.h"

#include <string>
#include <vector>
#include <unordered_map>


class LinearProgram;


class MATH_API ProgramElement
{
public:
	// A program element cannot belong to multiple models.
	// "program" is the program the owns this element.
	ProgramElement(LinearProgram* program, const std::string& name = "", int idx = 0) : program_(program), name_(name), index_(idx) {}

	const std::string& name() const { return name_; }
	void set_name(const std::string& n) { name_ = n; }

	int index() const { return index_; }
	void set_index(int idx) { index_ = idx; }

	// the program the owns this element
	const LinearProgram* program() const { return program_; }
	LinearProgram* program() { return program_; }

private:
	LinearProgram * program_; // the program the owns this element
	std::string		name_;
	int				index_;
};


class MATH_API Bound
{
public:
	// bound type for variables and expressions
enum BoundType { FIXED, LOWER, UPPER, DOUBLE, FREE };

public:
	Bound(BoundType bt = FREE, double lb = -infinity(), double ub = +infinity());

	BoundType bound_type() const { return bound_type_; }

	// if you are sure about the bound type
	void set_bounds(BoundType type, double lb, double ub);
	void set_bound(BoundType type, double value); // for FIXED, LOWER, UPPER

	// bound type will be automatically determined
	void set_bounds(double lb, double ub);

	double get_bound() const;	// query the single bound according to its bound type
	void   get_bounds(double& lb, double& ub) const;

	static double infinity();

private:
	BoundType	bound_type_;
	double		lower_bound_;
	double		upper_bound_;

	static double infinity_;
};


class MATH_API Variable : public Bound, public ProgramElement
{
public:
enum VariableType { CONTINUOUS, INTEGER, BINARY };

public:
	// A variable cannot belong to several models.
	// "program" is the program the owns this variable.
	Variable(LinearProgram* program, VariableType t = CONTINUOUS);

	VariableType variable_type() const { return variable_type_; }
	void set_variable_type(VariableType t);

	// Returns the value of the variable in the current solution.
	// Note: (1) valid only if the problem was successfully solved.
	//       (2) if the variable is integer and rounded == true, then the 
	//           value will be rounded to the nearest integer.
	double solution_value(bool rounded = false) const;

protected:
	friend class LinearProgramSolver;
	void set_solution_value(double value) { solution_value_ = value; }

private:
	VariableType variable_type_;
	double		 solution_value_;
};


class MATH_API LinearExpression : public ProgramElement
{
public:
	// An expression cannot belong to several models.
	// "program" is the program the owns this expression.
	LinearExpression(LinearProgram* program);

	void add_coefficient(int var_index, double coeff);	// coefficients can accumulate

	const std::unordered_map<int, double>& coefficients() const { return coefficients_; }

	// Evaluates the value of this expression at the solution found.
	// Note: (1) valid only if the problem was successfully solved.
	//       (2) if a variable is integer and rounded == true, then the 
	//           variable value will be rounded to the nearest integer.
	double solution_value(bool rounded = false) const;
	
	void clear() { coefficients_.clear(); }

private:
	std::unordered_map<int, double>	coefficients_;
};


class MATH_API LinearConstraint : public LinearExpression, public Bound
{
public:
	// A constraint cannot belong to several models.
	// "program" is the program the owns this constraint.
	LinearConstraint(LinearProgram* program, LinearConstraint::BoundType bt, double lb, double ub);
};


class MATH_API LinearObjective : public LinearExpression
{
public:
	enum Sense { MINIMIZE, MAXIMIZE, UNDEFINED };

public:
	// An objective cannot belong to several models.
	// "program" is the program the owns this objective.
	LinearObjective(LinearProgram* program, Sense s);

	void set_sense(Sense s) { sense_ = s; }
	Sense sense() const { return sense_; }

private:
	Sense sense_;
};


class MATH_API LinearProgram
{
public:
	LinearProgram();
	~LinearProgram();

	const std::string& name() const { return name_; }
	void set_name(const std::string& name) { name_ = name; }

	//////////////////////////////////////////////////////////////////////////

	// create a single variable, it to the program, and returns the pointer.
	// Note: if name is empty or not provided, a default name (e.g., x0, x1...) will be given.
	Variable* create_variable(
		Variable::VariableType vt = Variable::CONTINUOUS,
		Variable::BoundType bt = Variable::FREE,
		double lb = -Variable::infinity(),
		double ub = +Variable::infinity(),
		const std::string& name = ""
		);

	// create a set of variables and add them to the program.
	// Note: variables with be given default names, e.g., x0, x1...
	std::vector<Variable*> create_n_variables(std::size_t n);

	// create a single linear constraint, add it to the program, and returns the pointer.
	// Note: if name is empty or not provided, a default name (e.g., c0, c1...) will be given.
	LinearConstraint* create_constraint(
		LinearConstraint::BoundType bt = Variable::FREE,
		double lb = -Variable::infinity(),
		double ub = +Variable::infinity(),
		const std::string& name = ""
	);

	// create a set of linear constraints and add them to the program.	
	// Note: constraints with be given default names, e.g., c0, c1...
	std::vector<LinearConstraint*> create_n_constraints(std::size_t n);

	// create the objective function and returns the pointer.
	LinearObjective* create_objective(LinearObjective::Sense sense = LinearObjective::MINIMIZE);

	//////////////////////////////////////////////////////////////////////////

	std::size_t num_variables() const { return variables_.size(); }
	const std::vector<Variable*>& variables() const { return variables_; }
	std::vector<Variable*>& variables() { return variables_; }

	std::size_t num_constraints() const { return constraints_.size(); }
	const std::vector<LinearConstraint*>& constraints() const { return constraints_; }
	std::vector<LinearConstraint*>& constraints() { return constraints_; }

	const LinearObjective* objective() const;
	LinearObjective* objective();

	//////////////////////////////////////////////////////////////////////////

	std::size_t num_continuous_variables() const;
	std::size_t num_integer_variables() const;
	std::size_t num_binary_variables() const;

	bool is_continuous() const;				// returns true if all variables are continuous
	bool is_mix_integer_program() const;	// returns true if mixed inter program
	bool is_integer_program() const;		// returns true if inter program
	bool is_binary_proram() const;			// returns true if binary program
	
	//////////////////////////////////////////////////////////////////////////

	// clear all variables, constraints, and the objective.
	void clear();

	//////////////////////////////////////////////////////////////////////////

		// read/write linear program from/to a file. Format determined by file extension:
	//  - "lp":   CPLEX LP format. CPLEX .lp files with linear and quadratic constraints 
	//			  and objective, special ordered sets of type 1 and 2, indicators on linear 
	//			  constraints, and semi-continuous variables. For writing, linear (general 
	//			  and specialized), indicator, quadratic, second order cone, and special 
	//			  ordered set constraints are supported.
	//  - "mps":  MPS format. Allows to parse and write MPS files with linear and quadratic 
	//			  constraints and objective, special ordered sets of type 1 and 2, indicators 
	//			  on linear constraints, and semi-continuous variables. For writing, linear 
	//			  (general and specialized), indicator, quadratic, second order cone, and 
	//			  special ordered set constraints are supported.
	//			  See http://en.wikipedia.org/wiki/MPS_%28format%29 for a description.
	//  - "cip":  SCIP cip format. The CIP format consists of information written by the 
	//            individual constraints. Thus, the format is defined within the constraint 
	//            handlers. The CIP format is the only format within SCIP that allows to write 
	//            and read all constraints; all other file formats are restricted to some 
	//            particular sub-class of constraint integer programs.
	bool load(const std::string& file_name);

	// the parameter "use_simple_name" provides an option to save the variables/constrains' 
	// original names or simple names like x0, x1... and c0, c1...
	bool save(const std::string& file_name, bool use_simple_name = false) const;
	
private:
	std::string			name_;
	LinearObjective*	objective_;

	std::vector<Variable*>			variables_;
	std::vector<LinearConstraint*>	constraints_;
};


#endif
