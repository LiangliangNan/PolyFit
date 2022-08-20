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

/**@file   struct_pricestore.h
 * @ingroup INTERNALAPI
 * @brief  data structures for storing priced variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PRICESTORE_H__
#define __SCIP_STRUCT_PRICESTORE_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_var.h"
#include "scip/type_pricestore.h"

#ifdef __cplusplus
extern "C" {
#endif

/** storage for priced variables */
struct SCIP_Pricestore
{
   SCIP_CLOCK*           probpricingtime;    /**< time needed to price existing problem variables */
   SCIP_VAR**            vars;               /**< array with priced variables with violated reduced costs sorted by score */
   SCIP_Real*            scores;             /**< score for each priced variable (e.g. |redcost|/no. of nonzeros) */
   SCIP_VAR**            bdviolvars;         /**< variables where zero violates the bounds */
   SCIP_Real*            bdviolvarslb;       /**< lower bounds of bdviolvars */
   SCIP_Real*            bdviolvarsub;       /**< upper bounds of bdbiolvars */
   int                   varssize;           /**< size of vars and score arrays */
   int                   nvars;              /**< number of priced variables (max. is set->price_maxvars) */
   int                   bdviolvarssize;     /**< size of bdviolvars, bdviolvarslb, and bdviolvarsub arrays */
   int                   nbdviolvars;        /**< number of variables, where zero violates the bounds */
   int                   naddedbdviolvars;   /**< number of bound violated variables already added to the LP */
   int                   nprobpricings;      /**< total number of calls to problem variable pricing */
   int                   nprobvarsfound;     /**< total number of problem variables, that were added (and possibly thrown away) */
   int                   nvarsfound;         /**< total number of variables, that were added (and possibly thrown away) */
   int                   nvarsapplied;       /**< total number of variables, that were added to the LP */
   SCIP_Bool             initiallp;          /**< is the pricing storage currently being filled with the initial LP columns? */
};

#ifdef __cplusplus
}
#endif

#endif
