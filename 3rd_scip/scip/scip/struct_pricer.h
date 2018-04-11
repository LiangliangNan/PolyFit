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

/**@file   struct_pricer.h
 * @ingroup INTERNALAPI
 * @brief  data structures for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PRICER_H__
#define __SCIP_STRUCT_PRICER_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_pricer.h"

#ifdef __cplusplus
extern "C" {
#endif

/** variable pricers data */
struct SCIP_Pricer
{
   char*                 name;               /**< name of variable pricer */
   char*                 desc;               /**< description of variable pricer */
   SCIP_DECL_PRICERCOPY  ((*pricercopy));    /**< copy method of pricer or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRICERFREE  ((*pricerfree));    /**< destructor of variable pricer */
   SCIP_DECL_PRICERINIT  ((*pricerinit));    /**< initialize variable pricer */
   SCIP_DECL_PRICEREXIT  ((*pricerexit));    /**< deinitialize variable pricer */
   SCIP_DECL_PRICERINITSOL((*pricerinitsol));/**< solving process initialization method of variable pricer */
   SCIP_DECL_PRICEREXITSOL((*pricerexitsol));/**< solving process deinitialization method of variable pricer */
   SCIP_DECL_PRICERREDCOST((*pricerredcost));/**< reduced cost pricing method of variable pricer for feasible LPs */
   SCIP_DECL_PRICERFARKAS((*pricerfarkas));  /**< Farkas pricing method of variable pricer for infeasible LPs */
   SCIP_PRICERDATA*      pricerdata;         /**< variable pricers local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this pricer for the next stages */
   SCIP_CLOCK*           pricerclock;        /**< pricer execution time */
   int                   priority;           /**< priority of the variable pricer */
   int                   ncalls;             /**< number of times, this pricer was called */
   int                   nvarsfound;         /**< number of variables priced in found so far by this pricer */
   SCIP_Bool             delay;              /**< should the pricer be delayed until no other pricers or already existing
                                              *   problem variables with negative reduced costs are found */
   SCIP_Bool             active;             /**< is variable pricer in use for the current problem? */
   SCIP_Bool             initialized;        /**< is variable pricer initialized? */
};

#ifdef __cplusplus
}
#endif

#endif
