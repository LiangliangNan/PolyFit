#include "linear_program_solver.h"
#include "../basic/logger.h"
#include "../basic/basic_types.h"

bool LinearProgramSolver::_solve_GLPK(const LinearProgram* program) {
	Logger::warn("-") << "GLPK solver is in preparation and will be available soon..." << std::endl;
	return false;
}