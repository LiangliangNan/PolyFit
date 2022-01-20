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

#ifndef SOLVER_SUITE_LINEAR_PROGRAM_H
#define SOLVER_SUITE_LINEAR_PROGRAM_H

#include "linear_program_solver_common.h"

#include <string>
#include <vector>
#include <unordered_map>


class LinearProgram;


// The base class of linear program entries (e.g., variable, constraint, objective)

class SOLVER_API ProgramEntry
{
private:
	// A program entry cannot belong to multiple programs.
	// "program" owns this entry.
    ProgramEntry(LinearProgram* program, const std::string& name = "", int idx = 0)
        : program_(program)
        , name_(name)
        , index_(idx)
    {}

public:
	const std::string& name() const { return name_; }
	void set_name(const std::string& n) { name_ = n; }

	int  index() const { return index_; }
	void set_index(int idx) { index_ = idx; }

	// the program that owns this entry
	const LinearProgram* program() const { return program_; }
	LinearProgram* program() { return program_; }

private:
	LinearProgram * program_; // the program that owns this entry
	std::string		name_;
	int				index_;

private:
    friend class Variable;
	friend class LinearExpression;
	friend class LinearProgramSolver;
	friend class LinearProgram;
};


class SOLVER_API Bound
{
public:
	// bound type for variables and expressions
	enum BoundType { FIXED, LOWER, UPPER, DOUBLE, FREE };

private:
	Bound(double lb = -infinity(), double ub = +infinity());

public:
	// bound type is determined by the bound values
	BoundType bound_type() const;

	void   set_bounds(double lb, double ub);
	void   set_lower_bound(double lb) { lower_bound_ = lb; }
	void   set_upper_bound(double ub) { upper_bound_ = ub; }

	void   get_bounds(double& lb, double& ub) const;
	double lower_bound() const { return lower_bound_; }
	double upper_bound() const { return upper_bound_; }

	static double infinity();

private:
	double		lower_bound_;
	double		upper_bound_;

	static double infinity_;

private:
	friend class Variable;
	friend class LinearConstraint;
	friend class LinearProgramSolver;
	friend class LinearProgram;
};


class SOLVER_API Variable : public ProgramEntry, public Bound
{
public:
	enum VariableType { CONTINUOUS, INTEGER, BINARY };

private:
	// A variable cannot belong to several programs.
	// "program" owns this variable.
	Variable(LinearProgram* program, VariableType type = CONTINUOUS, double lb = -infinity(), double ub = +infinity(), const std::string& name = "", int idx = 0);

public:
	VariableType variable_type() const { return variable_type_; }
	void set_variable_type(VariableType t);

	// Returns the value of the variable in the current solution.
	// Note: (1) valid only if the problem was successfully solved.
	//       (2) if the variable is integer and rounded == true, then the 
	//           value will be rounded to the nearest integer.
	double solution_value(bool rounded = false) const;

private:
	void set_solution_value(double value) { solution_value_ = value; }

private:
	VariableType variable_type_;
	double		 solution_value_;

private:
    //copying disabled
    Variable(const Variable&);
    Variable& operator=(const Variable&);

    friend class LinearProgramSolver;
    friend class LinearProgram;
};


class SOLVER_API LinearExpression : public ProgramEntry
{
private:
	// An expression cannot belong to several programs.
	// "program" owns this expression.
	LinearExpression(LinearProgram* program, const std::string& name = "", int idx = 0);
    virtual ~LinearExpression() {}

public:
	// we use "add": coefficients can accumulate
	void add_coefficient(const Variable* var, double coeff);	

	double get_coefficient(const Variable* var) const;

	const std::unordered_map<const Variable*, double>& coefficients() const { return coefficients_; }

	// the constant term
	void set_offset(double value) { offset_ = value; }
	double offset() const { return offset_; }

	// evaluates the value of this expression at the solution found.
	// Note: (1) valid only if the problem was successfully solved.
	//       (2) if a variable is integer and rounded == true, then the 
	//           variable value will be rounded to the nearest integer.
	double solution_value(bool rounded = false) const;

	virtual void clear() { coefficients_.clear(); offset_ = 0.0; }

private:
	std::unordered_map<const Variable*, double>	coefficients_;
    double	offset_;

private:
    //copying disabled
    LinearExpression(const LinearExpression&);
    LinearExpression& operator=(const LinearExpression&);

    friend class LinearConstraint;
    friend class LinearObjective;
    friend class LinearProgramSolver;
    friend class LinearProgram;
};


class SOLVER_API LinearConstraint : public LinearExpression, public Bound
{
private:
	// A constraint cannot belong to several programs.
	// "program" owns this constraint.
	LinearConstraint(LinearProgram* program, double lb = -infinity(), double ub = +infinity(), const std::string& name = "", int idx = 0);
    virtual ~LinearConstraint() {}

private:
    //copying disabled
    LinearConstraint(const LinearConstraint&);
    LinearConstraint& operator=(const LinearConstraint&);

    friend class LinearProgramSolver;
    friend class LinearProgram;
};


class SOLVER_API LinearObjective : public LinearExpression
{
public:
	enum Sense { MINIMIZE, MAXIMIZE, UNDEFINED };

private:
	// An objective cannot belong to several programs.
	// "program" owns this objective.
	LinearObjective(LinearProgram* program, Sense sense);

public:
	void  set_sense(Sense sense) { sense_ = sense; }
	Sense sense() const { return sense_; }

	void clear();

private:
	Sense sense_;

private:
    //copying disabled
    LinearObjective(const LinearObjective&);
    LinearObjective& operator=(const LinearObjective&);

    friend class LinearProgramSolver;
    friend class LinearProgram;
};


class SOLVER_API LinearProgram
{
public:
    LinearProgram(const std::string& name = "no_name");
    ~LinearProgram();

	const std::string& name() const { return name_; }
	void set_name(const std::string& name) { name_ = name; }

	//////////////////////////////////////////////////////////////////////////

	// create a single variable, add it to the program, and return the pointer.
	// Note: if name is empty or not provided, a default name (e.g., x0, x1...) will be given.
	Variable* create_variable(
		Variable::VariableType type = Variable::CONTINUOUS,
		double lb = -Variable::infinity(),
		double ub = +Variable::infinity(),
		const std::string& name = ""
	);

	// create a set of variables and add them to the program.
	// Note: variables will be given default names, e.g., x0, x1...
	std::vector<Variable*> create_n_variables(std::size_t n);

	// create a single linear constraint, add it to the program, and returns the pointer.
	// Note: if name is empty or not provided, a default name (e.g., c0, c1...) will be given.
	LinearConstraint* create_constraint(
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

	bool has_variable(const Variable* var) const;
	bool has_constraint(const LinearConstraint* cons) const;

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

    // print statistics of the program to the stream.
    // the stream can be "std::cout"
	void print_statistics(std::ostream& output) const;

	// checks the program if there are issues and reports them to the output stream.
	// The issues that will be checked are:
	//  -) variables, constraints, and/or objective not owned by the program
	//  -) duplicated variables or constraints
	//  -) variables have the same name or index
	//  -) constraints have the same name or index
	//  -) variables with infeasible bounds (e.g., lb > ub)
	//  -) constraints with infeasible bounds (e.g., lb > ub)
    bool is_valid(std::ostream* output = nullptr) const;

	//////////////////////////////////////////////////////////////////////////

	// clear all variables, constraints, and the objective.
	void clear();

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

private:
    //copying disabled
    LinearProgram(const LinearProgram&);
    LinearProgram& operator=(const LinearProgram&);
};


#endif  // SOLVER_SUITE_LINEAR_PROGRAM_H
