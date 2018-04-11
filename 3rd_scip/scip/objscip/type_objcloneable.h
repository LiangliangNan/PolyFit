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

/**@file   type_objcloneable.h
 * @ingroup TYPEDEFINITIONS
 * @brief  function type definitions for clonable classes
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_OBJPCLONEABLE_H__
#define __SCIP_TYPE_OBJPCLONEABLE_H__

/** clone method, used to copy plugins which are not constraint handlers or variable pricer plugins
 *
 *  input:
 *  - scip            : SCIP main data structure
 */
#define SCIP_DECL_OBJCLONEABLECLONE(x) x (SCIP* scip) const
#define SCIP_DECL_BRANCHCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_DIALOGCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_DISPCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_EVENTCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_HEURCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_NODESELCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_PRESOLCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_PROPCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_READERCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_RELAXCLONE(x) x (SCIP* scip) const
#define SCIP_DECL_SEPACLONE(x) x (SCIP* scip) const


/** returns whether the plugin object is copyable
 *
 *  return value      : whether object is copyable
 */
#define SCIP_DECL_OBJCLONEABLEISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_BRANCHISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_DIALOGISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_DISPISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_EVENTISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_HEURISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_NODESELISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_PRESOLISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_PROPISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_READERISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_RELAXISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_SEPAISCLONEABLE(x) SCIP_Bool x (void) const

#endif
