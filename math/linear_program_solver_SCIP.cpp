#include "linear_program_solver.h"
#include "../basic/logger.h"
#include "../basic/basic_types.h"

bool LinearProgramSolver::_solve_SCIP(const LinearProgram* program) {
	Logger::warn("-") << "SCIP solver is in preparation and will be available soon..." << std::endl;
	return false;
}