cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(lpsolve_HEADERS
    fortify.h
    ini.h
    lp_bit.h
    lp_crash.h
    lp_explicit.h
    lp_Hash.h
    lp_lib.h
    lp_matrix.h
    lp_MDO.h
    lp_mipbb.h
    lp_MPS.h
    lp_presolve.h
    lp_price.h
    lp_pricePSE.h
    lp_report.h
    lp_rlp.h
    lp_scale.h
    lp_simplex.h
    lp_SOS.h
    lp_types.h
    lp_utils.h
    lp_wlp.h
    lpkit.h
    yacc_read.h
    shared/commonlib.h
    shared/mmio.h
    shared/myblas.h
    bfp/lp_BFP.h
    bfp/bfp_LUSOL/LUSOL/lusol.h
    bfp/bfp_LUSOL/bfp_LUSOL.h
    bfp/bfp_LUSOL/lp_LUSOL.h
    colamd/colamd.h
    )

set(lpsolve_SOURCES
    fortify.c
    ini.c
    lp_crash.c
    lp_Hash.c
    lp_lib.c
    lp_matrix.c
    lp_MDO.c
    lp_mipbb.c
    lp_MPS.c
    lp_params.c
    lp_presolve.c
    lp_price.c
    lp_pricePSE.c
    lp_report.c
    lp_rlp.c
    lp_scale.c
    lp_simplex.c
    lp_SOS.c
    lp_utils.c
    lp_wlp.c
    yacc_read.c
    shared/commonlib.c
    shared/mmio.c
    shared/myblas.c
    bfp/bfp_LUSOL/LUSOL/lusol.c
#    bfp/bfp_LUSOL/LUSOL/lusol1.c
#    bfp/bfp_LUSOL/LUSOL/lusol2.c
#    bfp/bfp_LUSOL/LUSOL/lusol6a.c
#    bfp/bfp_LUSOL/LUSOL/lusol6l0.c
#    bfp/bfp_LUSOL/LUSOL/lusol6u.c
#    bfp/bfp_LUSOL/LUSOL/lusol7a.c
#    bfp/bfp_LUSOL/LUSOL/lusol8a.c
    bfp/bfp_LUSOL/lp_LUSOL.c
    colamd/colamd.c
    )


add_library(${PROJECT_NAME} STATIC ${lpsolve_SOURCES} ${lpsolve_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
                      FOLDER "3rd_party")

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${POLYFIT_lpsolve_DIR}
    ${POLYFIT_lpsolve_DIR}/shared
    ${POLYFIT_lpsolve_DIR}/bfp
    ${POLYFIT_lpsolve_DIR}/bfp/bfp_LUSOL
    ${POLYFIT_lpsolve_DIR}/bfp/bfp_LUSOL/LUSOL
    ${POLYFIT_lpsolve_DIR}/colamd
    )

target_compile_definitions(${PROJECT_NAME} PRIVATE CHECK_SOLUTION YY_NEVER_INTERACTIVE PARSER_LP INVERSE_ACTIVE=INVERSE_LUSOL RoleIsExternalInvEngine)

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()
