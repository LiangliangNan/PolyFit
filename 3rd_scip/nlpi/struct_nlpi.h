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

/**@file   struct_nlpi.h
 * @brief  data definitions for an NLP solver interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_NLPI_H__
#define __SCIP_STRUCT_NLPI_H__

#include "scip/def.h"
#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** NLP interface data */
struct SCIP_Nlpi
{
   char*                           name;                        /**< name of NLP solver */
   char*                           description;                 /**< description of NLP solver */
   int                             priority;                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy));               /**< copy an NLPI */
   SCIP_DECL_NLPIFREE              ((*nlpifree));               /**< free NLPI user data */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer));   /**< get solver pointer */
   SCIP_DECL_NLPICREATEPROBLEM     ((*nlpicreateproblem));      /**< create a new problem instance */
   SCIP_DECL_NLPIFREEPROBLEM       ((*nlpifreeproblem));        /**< free a problem instance */
   SCIP_DECL_NLPIGETPROBLEMPOINTER ((*nlpigetproblempointer));  /**< get problem pointer */
   SCIP_DECL_NLPIADDVARS           ((*nlpiaddvars));            /**< add variables to a problem */
   SCIP_DECL_NLPIADDCONSTRAINTS    ((*nlpiaddconstraints));     /**< add constraints to a problem  */
   SCIP_DECL_NLPISETOBJECTIVE      ((*nlpisetobjective));       /**< set objective of a problem  */
   SCIP_DECL_NLPICHGVARBOUNDS      ((*nlpichgvarbounds));       /**< change variable bounds in a problem  */
   SCIP_DECL_NLPICHGCONSSIDES      ((*nlpichgconssides));       /**< change constraint sides in a problem  */
   SCIP_DECL_NLPIDELVARSET         ((*nlpidelvarset));          /**< delete a set of variables from a problem  */
   SCIP_DECL_NLPIDELCONSSET        ((*nlpidelconsset));         /**< delete a set of constraints from a problem */
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoefs));     /**< change coefficients in linear part of a constraint or objective */
   SCIP_DECL_NLPICHGQUADCOEFS      ((*nlpichgquadcoefs));       /**< change coefficients in quadratic part of a constraint or objective */
   SCIP_DECL_NLPICHGEXPRTREE       ((*nlpichgexprtree));        /**< change nonlinear expression a constraint or objective */
   SCIP_DECL_NLPICHGNONLINCOEF     ((*nlpichgnonlincoef));      /**< change one parameter in nonlinear expressions of a constraint or objective */
   SCIP_DECL_NLPICHGOBJCONSTANT    ((*nlpichgobjconstant));     /**< change the constant offset in the objective */
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess));    /**< set initial guess for primal variables in a problem  */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve));              /**< solve a problem */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat));         /**< get solution status for a problem  */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat));        /**< get termination status for a problem  */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution));        /**< get solution of a problem  */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics));      /**< get solve statistics for a problem  */
   SCIP_DECL_NLPIGETWARMSTARTSIZE  ((*nlpigetwarmstartsize));   /**< get size for warmstart object buffer for a problem  */
   SCIP_DECL_NLPIGETWARMSTARTMEMO  ((*nlpigetwarmstartmemo));   /**< get warmstart object for a problem  */
   SCIP_DECL_NLPISETWARMSTARTMEMO  ((*nlpisetwarmstartmemo));   /**< set warmstart object for a problem  */
   SCIP_DECL_NLPIGETINTPAR         ((*nlpigetintpar));          /**< get value of integer parameter in a problem  */
   SCIP_DECL_NLPISETINTPAR         ((*nlpisetintpar));          /**< set value of integer parameter in a problem  */
   SCIP_DECL_NLPIGETREALPAR        ((*nlpigetrealpar));         /**< get value of floating point parameter in a problem  */
   SCIP_DECL_NLPISETREALPAR        ((*nlpisetrealpar));         /**< set value of floating point parameter in a problem  */
   SCIP_DECL_NLPIGETSTRINGPAR      ((*nlpigetstringpar));       /**< get value of string parameter in a problem  */
   SCIP_DECL_NLPISETSTRINGPAR      ((*nlpisetstringpar));       /**< set value of string parameter in a problem  */
   SCIP_DECL_NLPISETMESSAGEHDLR    ((*nlpisetmessagehdlr));     /**< set message handler  */
   SCIP_NLPIDATA*                  nlpidata;                    /**< NLP interface local data */
};

/** Statistics from an NLP solve */
struct SCIP_NlpStatistics
{
   int       niterations;   /**< number of iterations the NLP solver spend in the last solve command */
   SCIP_Real totaltime;     /**< total time in CPU sections the NLP solver spend in the last solve command */ 
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_NLPI_H__ */
