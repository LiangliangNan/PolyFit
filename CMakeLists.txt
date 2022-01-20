set(POLYFIT_ROOT           ${CMAKE_CURRENT_LIST_DIR})
set(POLYFIT_3RD_PARTY_DIR  ${POLYFIT_ROOT}/3rd_party)
set(POLYFIT_scip_DIR       ${POLYFIT_3RD_PARTY_DIR}/3rd_scip)
set(POLYFIT_soplex_DIR     ${POLYFIT_3RD_PARTY_DIR}/3rd_soplex)

add_subdirectory(3rd_party)
add_subdirectory(linear_program_solver)
add_subdirectory(polyfit)
add_subdirectory(examples)