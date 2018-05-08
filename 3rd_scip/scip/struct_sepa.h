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

/**@file   struct_sepa.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SEPA_H__
#define __SCIP_STRUCT_SEPA_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_sepa.h"

#ifdef __cplusplus
extern "C" {
#endif

/** separators data */
struct SCIP_Sepa
{
   SCIP_Longint          lastsepanode;       /**< last (total) node where this separator was called */
   SCIP_Longint          ncalls;             /**< number of times, this separator was called */
   SCIP_Longint          ncutoffs;           /**< number of cutoffs found so far by this separator */
   SCIP_Longint          ncutsfound;         /**< number of cutting planes found so far by this separator */
   SCIP_Longint          ncutsapplied;       /**< number of cutting planes applied to LP */
   SCIP_Longint          nconssfound;        /**< number of additional constraints added by this separator */
   SCIP_Longint          ndomredsfound;      /**< number of domain reductions found so far by this separator */
   SCIP_Real             maxbounddist;       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   char*                 name;               /**< name of separator */
   char*                 desc;               /**< description of separator */
   SCIP_DECL_SEPACOPY    ((*sepacopy));      /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_SEPAFREE    ((*sepafree));      /**< destructor of separator */
   SCIP_DECL_SEPAINIT    ((*sepainit));      /**< initialize separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit));      /**< deinitialize separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol));   /**< solving process initialization method of separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol));   /**< solving process deinitialization method of separator */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp));    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol));   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata;           /**< separators local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this separator for the next stages */
   SCIP_CLOCK*           sepaclock;          /**< separation time */
   int                   priority;           /**< priority of the separator */
   int                   freq;               /**< frequency for calling separator */
   int                   ncallsatnode;       /**< number of times, this separator was called at the current node */
   int                   ncutsfoundatnode;   /**< number of cutting planes found at the current node */
   int                   expbackoff;         /**< base for exponential increase of frequency at which the separator is called */
   SCIP_Bool             usessubscip;        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay;              /**< should separator be delayed, if other separators found cuts? */
   SCIP_Bool             lpwasdelayed;       /**< was the LP separation delayed at the last call? */
   SCIP_Bool             solwasdelayed;      /**< was the solution separation delayed at the last call? */
   SCIP_Bool             initialized;        /**< is separator initialized? */
};

#ifdef __cplusplus
}
#endif

#endif
