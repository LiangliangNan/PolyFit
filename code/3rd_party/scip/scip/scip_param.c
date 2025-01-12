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

/**@file   scip_param.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for SCIP parameter handling
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/paramset.h"
#include "scip/pub_message.h"
#include "scip/scip_param.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"

/** creates a SCIP_Bool parameter, sets it to its default value, and adds it to the parameter set
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPaddBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Bool             defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPsetAddBoolParam(scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, valueptr, isadvanced,
         defaultvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a int parameter, sets it to its default value, and adds it to the parameter set
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPaddIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   int*                  valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   int                   defaultvalue,       /**< default value of the parameter */
   int                   minvalue,           /**< minimum value for parameter */
   int                   maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPsetAddIntParam(scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, valueptr, isadvanced,
         defaultvalue, minvalue, maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a SCIP_Longint parameter, sets it to its default value, and adds it to the parameter set
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPaddLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Longint          defaultvalue,       /**< default value of the parameter */
   SCIP_Longint          minvalue,           /**< minimum value for parameter */
   SCIP_Longint          maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPsetAddLongintParam(scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, valueptr, isadvanced,
         defaultvalue, minvalue, maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a SCIP_Real parameter, sets it to its default value, and adds it to the parameter set
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPaddRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Real             defaultvalue,       /**< default value of the parameter */
   SCIP_Real             minvalue,           /**< minimum value for parameter */
   SCIP_Real             maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPsetAddRealParam(scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, valueptr, isadvanced,
         defaultvalue, minvalue, maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a char parameter, sets it to its default value, and adds it to the parameter set
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPaddCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char*                 valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   char                  defaultvalue,       /**< default value of the parameter */
   const char*           allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPsetAddCharParam(scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, valueptr, isadvanced,
         defaultvalue, allowedvalues, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a string(char*) parameter, sets it to its default value, and adds it to the parameter set
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPaddStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char**                valueptr,           /**< pointer to store the current parameter value, or NULL; if not NULL then *valueptr should be NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   const char*           defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPsetAddStringParam(scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, valueptr, isadvanced,
         defaultvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** gets the fixing status of an existing parameter
 *
 *  @return TRUE if the parameter is fixed to a value, otherwise FALSE.
 */
SCIP_Bool SCIPisParamFixed(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsParamFixed(scip->set, name);
}

/** returns the pointer to the SCIP parameter with the given name
 *
 *  @return pointer to the parameter with the given name
 */
SCIP_PARAM* SCIPgetParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetGetParam(scip->set, name);
}

/** gets the value of an existing SCIP_Bool parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPgetBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool*            value               /**< pointer to store the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetGetBoolParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing int parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPgetIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   int*                  value               /**< pointer to store the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetGetIntParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Longint parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPgetLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint*         value               /**< pointer to store the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetGetLongintParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Real parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPgetRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Real*            value               /**< pointer to store the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetGetRealParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing char parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPgetCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char*                 value               /**< pointer to store the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetGetCharParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing string(char*) parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPgetStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char**                value               /**< pointer to store the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetGetStringParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** fixes the value of an existing parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note: Be careful with this method! Some general settings, e.g., the time or node limit, should not be fixed because
 *         they have to be changed for sub-SCIPs.
 */
SCIP_RETCODE SCIPfixParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetChgParamFixed(scip->set, name, TRUE) );

   return SCIP_OKAY;
}

/** unfixes the value of an existing parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPunfixParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetChgParamFixed(scip->set, name, FALSE) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Bool parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(scip->set != NULL);

   retcode = SCIPsetChgBoolParam(scip->set, scip->messagehdlr, param, value);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing SCIP_Bool parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetSetBoolParam(scip->set, scip->messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** checks whether the value of an existing SCIP_Bool parameter is valid */
SCIP_Bool SCIPisBoolParamValid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Bool             value               /**< value to check */
   )
{
   assert(scip != NULL);
   assert(param != NULL);

   return SCIPparamIsValidBool(param, value);
}

/** changes the value of an existing int parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   int                   value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(scip->set != NULL);

   retcode = SCIPsetChgIntParam(scip->set, scip->messagehdlr, param, value);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing int parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   int                   value               /**< new value of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetSetIntParam(scip->set, scip->messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** checks whether parameter value of an existing int paramter is valid */
SCIP_Bool SCIPisIntParamValid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   int                   value               /**< value to check */
   )
{
   assert(scip != NULL);
   assert(param != NULL);

   return SCIPparamIsValidInt(param, value);
}

/** changes the value of an existing SCIP_Longint parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(scip->set != NULL);

   retcode = SCIPsetChgLongintParam(scip->set, scip->messagehdlr, param, value);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing SCIP_Longint parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetSetLongintParam(scip->set, scip->messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** checks whether parameter value of an existing SCIP_Longint paramter is valid */
SCIP_Bool SCIPisLongintParamValid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Longint          value               /**< value to check */
   )
{
   assert(scip != NULL);
   assert(param != NULL);

   return SCIPparamIsValidLongint(param, value);
}

/** changes the value of an existing SCIP_Real parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Real             value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(scip->set != NULL);

   retcode = SCIPsetChgRealParam(scip->set, scip->messagehdlr, param, value);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing SCIP_Real parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Real             value               /**< new value of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetSetRealParam(scip->set, scip->messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** checks whether parameter value of an existing SCIP_Real paramter is valid */
SCIP_Bool SCIPisRealParamValid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Real             value               /**< value to check */
   )
{
   assert(scip != NULL);
   assert(param != NULL);

   return SCIPparamIsValidReal(param, value);
}

/** changes the value of an existing char parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   char                  value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(scip->set != NULL);

   retcode = SCIPsetChgCharParam(scip->set, scip->messagehdlr, param, value);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing char parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char                  value               /**< new value of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetSetCharParam(scip->set, scip->messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** checks whether parameter value for a given SCIP_Real parameter is valid */
SCIP_Bool SCIPisCharParamValid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   const char            value               /**< value to check */
   )
{
   assert(scip != NULL);
   assert(param != NULL);

   return SCIPparamIsValidChar(param, value);
}

/** changes the value of an existing string(char*) parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   const char*           value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(scip->set != NULL);

   retcode = SCIPsetChgStringParam(scip->set, scip->messagehdlr, param, value);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing string(char*) parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           value               /**< new value of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetSetStringParam(scip->set, scip->messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** checks whether parameter value for a given string parameter is valid */
SCIP_Bool SCIPisStringParamValid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   const char*           value               /**< value to check */
   )
{
   assert(scip != NULL);
   assert(param != NULL);

   return SCIPparamIsValidString(param, value);
}

/** reads parameters from a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPreadParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetReadParams(scip->set, scip->messagehdlr, filename) );

   return SCIP_OKAY;
}

/** writes a single parameter to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPwriteParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   const char*           filename,           /**< file name, or NULL for stdout */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only those parameters be written that are changed from their
                                              *   default value?
                                              */
   )
{
   assert(scip != NULL);
   assert(param != NULL);

   SCIP_CALL( SCIPparamWrite(param, scip->messagehdlr, filename, comments, onlychanged) );

   return SCIP_OKAY;
}

/** writes all parameters in the parameter set to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPwriteParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< file name, or NULL for stdout */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only those parameters be written that are changed from their
                                              *   default value?
                                              */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetWriteParams(scip->set, scip->messagehdlr, filename, comments, onlychanged) );

   return SCIP_OKAY;
}

/** resets a single parameter to its default value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPresetParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the parameter */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetResetParam(scip->set, scip->messagehdlr, name) );

   return SCIP_OKAY;
}

/** resets all parameters to their default values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPresetParams(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetResetParams(scip->set, scip->messagehdlr) );

   return SCIP_OKAY;
}

/** sets parameters to
 *
 *  - \ref SCIP_PARAMEMPHASIS_DEFAULT to use default values (see also SCIPresetParams())
 *  - \ref SCIP_PARAMEMPHASIS_COUNTER to get feasible and "fast" counting process
 *  - \ref SCIP_PARAMEMPHASIS_CPSOLVER to get CP like search (e.g. no LP relaxation)
 *  - \ref SCIP_PARAMEMPHASIS_EASYCIP to solve easy problems fast
 *  - \ref SCIP_PARAMEMPHASIS_FEASIBILITY to detect feasibility fast
 *  - \ref SCIP_PARAMEMPHASIS_HARDLP to be capable to handle hard LPs
 *  - \ref SCIP_PARAMEMPHASIS_OPTIMALITY to prove optimality fast
 *  - \ref SCIP_PARAMEMPHASIS_PHASEFEAS to find feasible solutions during a 3 phase solution process
 *  - \ref SCIP_PARAMEMPHASIS_PHASEIMPROVE to find improved solutions during a 3 phase solution process
 *  - \ref SCIP_PARAMEMPHASIS_PHASEPROOF to proof optimality during a 3 phase solution process
 *  - \ref SCIP_PARAMEMPHASIS_NUMERICS to solve problems which cause numerical issues
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetEmphasis(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMEMPHASIS    paramemphasis,      /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetSetEmphasis(scip->set, scip->messagehdlr, paramemphasis, quiet) );

   return SCIP_OKAY;
}

/** sets parameters to deactivate separators and heuristics that use auxiliary SCIP instances; should be called for
 *  auxiliary SCIP instances to avoid recursion
 *
 *  @note only deactivates plugins which could cause recursion, some plugins which use sub-SCIPs stay activated
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetSubscipsOff(
   SCIP*                 scip,               /**< (auxiliary) SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPsetSetSubscipsOff(scip->set, scip->messagehdlr, quiet) );

   return SCIP_OKAY;
}

/** sets heuristic parameters values to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spent on heuristics is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristics are called more aggressively
 *  - SCIP_PARAMSETTING_OFF which turn off all heuristics
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetHeuristics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(paramsetting == SCIP_PARAMSETTING_DEFAULT || paramsetting == SCIP_PARAMSETTING_FAST
      || paramsetting == SCIP_PARAMSETTING_AGGRESSIVE || paramsetting == SCIP_PARAMSETTING_OFF);

   SCIP_CALL( SCIPsetSetHeuristics(scip->set, scip->messagehdlr, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** sets presolving parameters to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all presolving parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spent on presolving is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the presolving is more aggressive
 *  - SCIP_PARAMSETTING_OFF which turn off all presolving
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(paramsetting == SCIP_PARAMSETTING_DEFAULT || paramsetting == SCIP_PARAMSETTING_FAST
      || paramsetting == SCIP_PARAMSETTING_AGGRESSIVE || paramsetting == SCIP_PARAMSETTING_OFF);

   SCIP_CALL( SCIPsetSetPresolving(scip->set, scip->messagehdlr, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** sets separating parameters to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separating parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spent on separating is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that separating is more aggressive
 *  - SCIP_PARAMSETTING_OFF which turn off all separating
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetSeparating(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(paramsetting == SCIP_PARAMSETTING_DEFAULT || paramsetting == SCIP_PARAMSETTING_FAST
      || paramsetting == SCIP_PARAMSETTING_AGGRESSIVE || paramsetting == SCIP_PARAMSETTING_OFF);

   SCIP_CALL( SCIPsetSetSeparating(scip->set, scip->messagehdlr, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** returns the array of all available SCIP parameters
 *
 *  @return SCIP_PARAM* array, containing all SCIP parameters.
 */
SCIP_PARAM** SCIPgetParams(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetGetParams(scip->set);
}

/** returns the total number of all available SCIP parameters
 *
 *  @return number of all SCIP parameters.
 */
int SCIPgetNParams(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetGetNParams(scip->set);
}

/** returns whether plugins with sub-SCIPs that could cause recursion have been disabled
 *
 *  @return the value of the variable set->subscipsoff
 */
SCIP_Bool SCIPgetSubscipsOff(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetGetSubscipsOff(scip->set);
}
