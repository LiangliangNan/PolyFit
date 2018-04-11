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

/**@file   message.h
 * @ingroup INTERNALAPI
 * @brief  message output methods
 * @author Tobias Achterberg
 *
 * Because the message functions are implemented as defines with more than one
 * function call, they shouldn't be used as a single statement like in:
 *   if( error )
 *      SCIPerrorMessage("an error occured");
 * because this would produce the following macro extension:
 *   if( error )
 *      printf(("[%s:%d] ERROR: ", __FILE__, __LINE__);
 *   printf(("an error occured");
 * Instead, they should be protected with brackets:
 *   if( error )
 *   {
 *      SCIPerrorMessage("an error occured");
 *   }
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MESSAGE_H__
#define __SCIP_MESSAGE_H__


#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
}
#endif

#endif

