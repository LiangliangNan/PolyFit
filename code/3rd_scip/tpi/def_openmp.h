/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   def_openmp.h
 * @ingroup TASKINTERFACE
 * @brief  wrappers for OpenMP defines
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DEF_OPENMP_H__
#define __DEF_OPENMP_H__

#define STR(x)                #x
#define STRINGIFY(x)          STR(x)
#define CONCATENATE(x, y)     x y
#define CONCATPARENTH(x, y)   x ( y )

#define TPI_NULL  NULL

#ifdef TPI_OMP

#define TPI_PRAGMA_CLAUSE(directive, clause)    _Pragma( STRINGIFY( CONCATENATE( directive, clause ) ) )
#define TPI_PRAGMA(directive)               _Pragma( STRINGIFY( directive ) )
#define TPI_PRAGMA_PARENTH(directive, var)  _Pragma( STRINGIFY( CONCATPARENTH( directive, var ) ) )

#else

#define TPI_PRAGMA_CLAUSE(directive, clause)
#define TPI_PRAGMA(directive)
#define TPI_PRAGMA_PARENTH(directive, var)

#endif


#define TPI_FOR_CLAUSE(clause)      TPI_PRAGMA_CLAUSE( TPI_DIR_FOR, clause )
#define TPI_FOR                     TPI_PRAGMA( TPI_DIR_FOR )
#define TPI_PARA_CLAUSE(clause)     TPI_PRAGMA_CLAUSE( TPI_DIR_PARA, clause )

#define TPI_PARA_CLAUSE_SHARED(priv, clause)       TPI_PRAGMA_CLAUSE( TPI_DIR_PARA,                           \
                                                      TPI_CLAUSE_DEFAULT( shared )                            \
                                                      TPI_CLAUSE_PRIVATE( (priv) ) clause )

#define TPI_PARA_SHARED                            TPI_PRAGMA_CLAUSE( TPI_DIR_PARA,                           \
                                                      TPI_CLAUSE_DEFAULT( shared ) )

#define TPI_PARA_SHARED_PRIVATE(priv)              TPI_PRAGMA_CLAUSE( TPI_DIR_PARA,                           \
                                                      TPI_CLAUSE_DEFAULT( shared )                            \
                                                      TPI_CLAUSE_PRIVATE( ( priv ) ) )

#define TPI_PARA_CLAUSE_NONE(share, priv, clause)  TPI_PRAGMA_CLAUSE( TPI_DIR_PARA,                           \
                                                      TPI_CLAUSE_DEFAULT( (none) )                            \
                                                      TPI_CLAUSE_SHARED( (share) )                            \
                                                      TPI_CLAUSE_PRIVATE( (priv) ) clause )

#define TPI_PARA                    TPI_PRAGMA( TPI_DIR_PARA )
#define TPI_CRITICAL(var)           TPI_PRAGMA_PARENTH( TPI_DIR_CRITICAL, var)
#define TPI_MASTER                  TPI_PRAGMA( TPI_DIR_MASTER )
#define TPI_WAIT                    TPI_PRAGMA( TPI_DIR_WAIT )
#define TPI_ORDERED                 TPI_PRAGMA( TPI_DIR_ORDERED )
#define TPI_SINGLE                  TPI_PRAGMA( TPI_DIR_SINGLE )
#define TPI_CLAUSE_SINGLE(clause)   TPI_PRAGMA_CLAUSE( TPI_DIR_SINGLE, clause )
#define TPI_TASK                    TPI_PRAGMA( TPI_DIR_TASK )
#define TPI_TASK_SHARED             TPI_PRAGMA_CLAUSE( TPI_DIR_TASK,                                           \
                                       TPI_CLAUSE_DEFAULT(shared) )
#define TPI_CLAUSE_TASK(clause)     TPI_PRAGMA_CLAUSE( TPI_DIR_TASK, clause )
#define TPI_TASKWAIT                TPI_PRAGMA( TPI_DIR_TASKWAIT )


/* OpenMP pragma directives */
#define TPI_DIR_PARA             omp parallel
#define TPI_DIR_FOR              omp for
#define TPI_DIR_CRITICAL         omp critical
#define TPI_DIR_MASTER           omp master
#define TPI_DIR_WAIT             omp barrier
#define TPI_DIR_ORDERED          omp ordered
#define TPI_DIR_TASK             omp task
#define TPI_DIR_SINGLE           omp single
#define TPI_DIR_TASKWAIT         omp taskwait


/* OpenMP clauses */
#define TPI_CLAUSE_PRIVATE(var)                 CONCATENATE( private, var )
#define TPI_CLAUSE_FSTPRIVATE(var)              CONCATENATE( firstprivate, var )
#define TPI_CLAUSE_LSTPRIVATE(var)              CONCATENATE( lastprivate, var )
#define TPI_CLAUSE_CPYPRIVATE(var)              CONCATENATE( copyprivate, var )
#define TPI_CLAUSE_NOWAIT                       nowait
#define TPI_CLAUSE_SHARED(var)                  CONCATENATE( shared, var )
#define TPI_CLAUSE_DEFAULT(var)                 CONCATPARENTH( default, var )
/* The reduce clause requires op as either an operator or intrinsic procedure.
 * Operators: +, *, .and., .or., .eqv., .neqv.
 * intrinsic procedures: max, min, iand, ior, or ieor*/
#define TPI_CLAUSE_REDUCE(op, var)              CONCATENATE( reduction, CONCATENATE( CONCATENATE( op, : ), var ) )
#define TPI_CLAUSE_ORDERED                      ordered
#define TPI_CLAUSE_IF(var)                      CONCATENATE( if, var )
#define TPI_CLAUSE_NUMTHREADS(var)              CONCATENATE( num_threads, var )
#define TPI_CLAUSE_SCHEDULE(type)               CONCATENATE( schedule, type )
#define TPI_CLAUSE_SCHEDULE_CHUNK(type, chunk)  CONCATENATE( schedule, CONCATPARENTH( type, chunk ) )
#define TPI_CLAUSE_COPYIN(var)                  CONCATENATE( copyin, var )
#define TPI_CLAUSE_FINAL(var)                   CONCATENATE( final, var )
#define TPI_CLAUSE_UNTIED                       untied
#define TPI_CLAUSE_MERGEABLE                    mergeable
#define TPI_CLAUSE_DEPEND(type, var)            CONCATENATE( depend, CONCATENATE( CONCATENATE( type, : ), var ) )
#define TPI_CLAUSE_PRIORITY(var)                CONCATENATE( priority, var )



#define TPI_SHARED_DATA(name, members)       struct TPI_Shared_Data {                                         \
                                                members                                                       \
                                             } name;

#define TPI_PRIVATE_DATA(name, members)      struct TPI_Private_Data {                                        \
                                                members                                                       \
                                             } name;

#define TPI_FSTPRIVATE_DATA(name, members)   struct TPI_FirstPrivate_Data {                                   \
                                                members                                                       \
                                             } name;

#define TPI_LSTPRIVATE_DATA(name, members)   struct TPI_LastPrivate_Data {                                    \
                                                members                                                       \
                                             } name;

#define TPI_COPYIN_DATA(name, members)       struct TPI_CopyIn_Data {                                         \
                                                members                                                       \
                                             } name;


#endif
