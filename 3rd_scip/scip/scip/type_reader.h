/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_reader.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_READER_H__
#define __SCIP_TYPE_READER_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Reader SCIP_READER;               /**< reader data structure */
typedef struct SCIP_ReaderData SCIP_READERDATA;       /**< reader specific data */


/** copy method for reader plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - reader          : the reader itself
 */
#define SCIP_DECL_READERCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_READER* reader)


/** destructor of reader to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - reader          : the reader itself
 */
#define SCIP_DECL_READERFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_READER* reader)

/** problem reading method of reader
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - reader          : the reader itself
 *  - filename        : full path and name of file to read, or NULL if stdin should be used
 *  - result          : pointer to store the result of the file reading call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropriate problem
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERROR or SCIP_NOFILE.
 */
#define SCIP_DECL_READERREAD(x) SCIP_RETCODE x (SCIP* scip, SCIP_READER* reader, const char* filename, SCIP_RESULT* result)

/** problem writing method of reader; NOTE: if the parameter "genericnames" is TRUE, then
 *  SCIP already set all variable and constraint names to generic names; therefore, this
 *  method should always use SCIPvarGetName() and SCIPconsGetName(); 
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - reader          : the reader itself
 *  - file            : output file, or NULL if standard output should be used
 *  - name            : problem name
 *  - probdata        : user problem data set by the reader
 *  - transformed     : TRUE iff problem is the transformed problem
 *  - objsense        : objective sense
 *  - objscale        : scalar applied to objective function; external objective value is
                        extobj = objsense * objscale * (intobj + objoffset)
 *  - objoffset       : objective offset from bound shifting and fixing 
 *  - vars            : array with active variables ordered binary, integer, implicit, continuous 
 *  - nvars           : number of active variables in the problem
 *  - nbinvars        : number of binary variables
 *  - nintvars        : number of general integer variables
 *  - nimplvars       : number of implicit integer variables 
 *  - ncontvars;      : number of continuous variables
 *  - fixedvars       : array with fixed and aggregated variables
 *  - nfixedvars      : number of fixed and aggregated variables in the problem
 *  - startnvars      : number of variables existing when problem solving started
 *  - conss           : array with constraints of the problem
 *  - nconss          : number of constraints in the problem
 *  - maxnconss       : maximum number of constraints existing at the same time 
 *  - startnconss     : number of constraints existing when problem solving started
 *  - genericnames    : using generic variable and constraint names?
 *  - result          : pointer to store the result of the file reading call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader wrote the file correctly
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error while writing the output file, it should return with RETCODE SCIP_WRITEERROR 
 */
#define SCIP_DECL_READERWRITE(x) SCIP_RETCODE x (SCIP* scip, SCIP_READER* reader, FILE* file, \
      const char* name, SCIP_PROBDATA* probdata, SCIP_Bool transformed, \
      SCIP_OBJSENSE objsense, SCIP_Real objscale, SCIP_Real objoffset,  \
      SCIP_VAR** vars, int nvars, int nbinvars, int nintvars, int nimplvars, int ncontvars, \
      SCIP_VAR** fixedvars, int nfixedvars, int startnvars, \
      SCIP_CONS** conss, int nconss, int maxnconss, int startnconss, \
      SCIP_Bool genericnames, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
