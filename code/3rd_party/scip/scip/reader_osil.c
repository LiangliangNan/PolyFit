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

/**@file   reader_osil.c
 * @ingroup DEFPLUGINS_READER
 * @brief  OS instance language (OSiL) format file reader
 * @author Stefan Vigerske
 * @author Ingmar Vierhaus
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI and M_E on Windows */  /*lint !750 */
#include "blockmemshell/memory.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_linear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/expr_abs.h"
#include "scip/expr_erf.h"
#include "scip/expr_exp.h"
#include "scip/expr_log.h"
#include "scip/expr_pow.h"
#include "scip/expr_product.h"
#include "scip/expr_sum.h"
#include "scip/expr_trig.h"
#include "scip/expr_value.h"
#include "scip/expr_var.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_nlp.h"
#include "scip/pub_var.h"
#include "scip/reader_osil.h"
#include "scip/scip_cons.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"
#include "scip/scip_var.h"
#include <stdlib.h>
#include <string.h>
#include "xml/xml.h"

#define READER_NAME             "osilreader"
#define READER_DESC             "file reader for OS instance language (OSiL) format"
#define READER_EXTENSION        "osil"

/*
 * Local methods
 */

/** create variables with bounds and type according to xml data */
static
SCIP_RETCODE readVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR***           vars,               /**< buffer to store pointer to variable array */
   int*                  nvars,              /**< buffer to store number of variables */
   SCIP_Bool             initialconss,       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss,       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamiccols,        /**< should columns be added and removed dynamically to the LP? */
   SCIP_Bool             dynamicrows,        /**< should rows be added and removed dynamically to the LP? */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* variables;
   const XML_NODE* varnode;
   const char* attrval;
   int varssize;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL);
   assert(nvars != NULL);
   assert(doingfine != NULL);

   *vars = NULL;
   *nvars = 0;

   variables = xmlFindNodeMaxdepth(datanode, "variables", 0, 1);

   if( variables == NULL )
   {
      /* no variables: strange but ok so far */
      return SCIP_OKAY;
   }

   /* get number of variables */
   attrval = xmlGetAttrval(variables, "numberOfVariables");
   if( attrval == NULL )
   {
      SCIPerrorMessage("Attribute \"numberOfVariables\" not found in <variables> node.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   varssize = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || varssize < 0 )
   {
      SCIPerrorMessage("Invalid value '%s' for \"numberOfVariables\" attribute.\n", xmlGetAttrval(variables, "numberOfVariables"));
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(varssize >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, vars, varssize) );

   /* parse variable nodes, create SCIP vars and add to problem
    * create bounddisjunction constraints for semicontinuous/semiinteger variables
    */
   for( varnode = xmlFirstChild(variables); varnode != NULL; varnode = xmlNextSibl(varnode) )
   {
      const char* varname;
      SCIP_VARTYPE vartype;
      SCIP_Real varlb;
      SCIP_Real varub;
      SCIP_Real semibound;

      if( varssize == *nvars )
      {
         SCIPerrorMessage("Expected %d variables, got at least %d many.\n", varssize, *nvars+1);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* find variable name */
      varname = xmlGetAttrval(varnode, "name");

      /* check for mult attribute */
      attrval = xmlGetAttrval(varnode, "mult");
      if( attrval != NULL && strcmp(attrval, "1") != 0 )
      {
         SCIPerrorMessage("Variable attribute 'mult' not supported (while parsing variable <%s>)\n", varname);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* find variable lower bound (default is 0.0 !) */
      attrval = xmlGetAttrval(varnode, "lb");
      if( attrval == NULL )
         varlb = 0.0;
      else if( strcmp(attrval, "-INF") == 0 )
         varlb = -SCIPinfinity(scip);
      else if( strcmp(attrval, "INF") == 0 )
         varlb = SCIPinfinity(scip);
      else
      {
         varlb = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("Error parsing variable lower bound '%s' for variable <%s>\n", attrval, varname);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }

      /* find variable upper bound (default is infinity) */
      attrval = xmlGetAttrval(varnode, "ub");
      if( attrval == NULL )
         varub = SCIPinfinity(scip);
      else if( strcmp(attrval, "-INF") == 0 )
         varub = -SCIPinfinity(scip);
      else if( strcmp(attrval, "INF") == 0 )
         varub = SCIPinfinity(scip);
      else
      {
         varub = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("Error parsing variable upper bound '%s' for variable <%s>\n", attrval, varname);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }

      semibound = SCIP_INVALID;

      /* find variable type (default is continuous)
       * adjust variable lower bound for semicontinuous variables
       */
      attrval = xmlGetAttrval(varnode, "type");
      if( attrval == NULL )
         vartype = SCIP_VARTYPE_CONTINUOUS;
      else switch( *attrval )
      {
      case 'C':
         vartype = SCIP_VARTYPE_CONTINUOUS;
         break;
      case 'B':
         vartype = SCIP_VARTYPE_BINARY;
         if( varub > 1.0 )
            varub = 1.0;
         break;
      case 'I':
         vartype = SCIP_VARTYPE_INTEGER;
         break;
      case 'D':
         vartype = SCIP_VARTYPE_CONTINUOUS;
         if( varlb > 0.0 )
            semibound = varlb;
         varlb = 0.0;
         break;
      case 'J':
         vartype = SCIP_VARTYPE_INTEGER;
         if( varlb > 0.0 )
            semibound = varlb;
         varlb = 0.0;
         break;
      default:
         SCIPerrorMessage("Unsupported variable type '%s' for variable <%s>\n", attrval, varname);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      if( vartype != SCIP_VARTYPE_CONTINUOUS )
      {
         varlb = SCIPceil(scip, varlb);
         varub = SCIPfloor(scip, varub);
      }

      /* create SCIP variable */
      SCIP_CALL( SCIPcreateVar(scip, &(*vars)[*nvars], varname, varlb, varub, 0.0, vartype, !dynamiccols, dynamiccols, NULL, NULL, NULL, NULL, NULL) );
      assert((*vars)[*nvars] != NULL);

      /* add variable to problem */
      SCIP_CALL( SCIPaddVar(scip, (*vars)[*nvars]) );

      /* if variable is actually semicontinuous or semiintegral, create bounddisjunction constraint (var <= 0.0 || var >= semibound) */
      if( semibound != SCIP_INVALID )  /*lint !e777*/
      {
         SCIP_CONS* cons;
         SCIP_VAR* consvars[2];
         SCIP_BOUNDTYPE boundtypes[2];
         SCIP_Real bounds[2];
         char name[SCIP_MAXSTRLEN];

         consvars[0] = (*vars)[*nvars];
         consvars[1] = (*vars)[*nvars];

         boundtypes[0] = SCIP_BOUNDTYPE_UPPER;
         boundtypes[1] = SCIP_BOUNDTYPE_LOWER;

         bounds[0] = 0.0;
         bounds[1] = semibound;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_semibound", SCIPvarGetName((*vars)[*nvars]));

         SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, 2, consvars, boundtypes, bounds,
               initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }

      ++*nvars;
   }
   if( *nvars < varssize )
   {
      SCIPerrorMessage("Expected %d variables, but got only %d many.\n", varssize, *nvars);
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** setup linear coefficients and constant of objective and objective sense */
static
SCIP_RETCODE readObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             dynamiccols,        /**< should columns be added and removed dynamically to the LP? */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* objective;
   const XML_NODE* coefnode;
   const char* attrval;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(doingfine != NULL);

   /* check for first objective */
   objective = xmlFindNodeMaxdepth(datanode, "obj", 0, 2);

   /* if no objective, then nothing to do here */
   if( objective == NULL )
      return SCIP_OKAY;

   /* check for mult attribute */
   attrval = xmlGetAttrval(objective, "mult");
   if( attrval != NULL && strcmp(attrval, "1") != 0 )
   {
      SCIPerrorMessage("Objective attribute 'mult' not supported.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   /* objective sense */
   attrval = xmlGetAttrval(objective, "maxOrMin");
   if( attrval == NULL )
   {
      SCIPerrorMessage("Objective sense missing.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   else if( strcmp(attrval, "min") == 0 )
   {
      SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   }
   else if( strcmp(attrval, "max") == 0 )
   {
      SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   }
   else
   {
      SCIPerrorMessage("Cannot parse objective sense '%s'.\n", attrval);
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   /* objective coefficients */
   for( coefnode = xmlFirstChild(objective); coefnode != NULL; coefnode = xmlNextSibl(coefnode) )
   {
      SCIP_Real val;
      int idx;

      /* get variable index */
      attrval = xmlGetAttrval(coefnode, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing \"idx\" attribute in objective coefficient.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      idx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' )
      {
         SCIPerrorMessage("Error parsing variable index '%s' of objective coefficient.\n", xmlGetAttrval(coefnode, "idx"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      if( idx < 0 || idx >= nvars )
      {
         SCIPerrorMessage("Invalid variable index '%d' of objective coefficient.\n", idx);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* get coefficient value */
      if( xmlFirstChild(coefnode) == NULL || xmlGetData(xmlFirstChild(coefnode)) == NULL )
      {
         SCIPerrorMessage("No objective coefficient stored for %d'th variable (<%s>).\n", idx, SCIPvarGetName(vars[idx]));  /*lint !e613*/
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      attrval = xmlGetData(xmlFirstChild(coefnode));
      val = strtod(attrval, (char**)&attrval);
      if( *attrval != '\0' )
      {
         SCIPerrorMessage("Error parsing objective coefficient value '%s' for %d'th variable (<%s>).\n", xmlGetData(xmlFirstChild(coefnode)), idx, SCIPvarGetName(vars[idx]));  /*lint !e613*/
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* change objective coefficient of SCIP variable */
      SCIP_CALL( SCIPchgVarObj(scip, vars[idx], val) );  /*lint !e613*/
   }

   /* objective constant: model as fixed variable, if nonzero */
   attrval = xmlGetAttrval(objective, "constant");
   if( attrval != NULL )
   {
      SCIP_Real objconst;

      objconst = strtod(attrval, (char**)&attrval);
      if( *attrval != '\0' )
      {
         SCIPerrorMessage("Error parsing objective constant '%s'\n", xmlGetAttrval(objective, "constant"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      if( objconst != 0.0 )
      {
         SCIP_VAR* objconstvar;

         SCIP_CALL( SCIPcreateVar(scip, &objconstvar, "objconstvar", objconst, objconst, 1.0, SCIP_VARTYPE_CONTINUOUS, !dynamiccols, dynamiccols, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, objconstvar) );
         SCIP_CALL( SCIPreleaseVar(scip, &objconstvar) );
      }
   }

   if( xmlNextSibl(objective) != NULL )
   {
      SCIPerrorMessage("Multiple objectives not supported by SCIP.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** helper method to get the total number of constraints */
static
SCIP_RETCODE readNConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   int*                  nconss,             /**< pointer to store the total number of constraints */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* constraints;
   const char* attrval;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(nconss != NULL);
   assert(doingfine != NULL);

   *nconss = 0;

   constraints = xmlFindNodeMaxdepth(datanode, "constraints", 0, 1);

   /* if no constraints, then nothing to do here */
   if( constraints == NULL )
      return SCIP_OKAY;

   /* read number of constraints */
   attrval = xmlGetAttrval(constraints, "numberOfConstraints");
   if( attrval == NULL )
   {
      SCIPerrorMessage("Attribute \"numberOfConstraints\" not found in <constraints> node.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   *nconss = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || *nconss < 0 )
   {
      SCIPerrorMessage("Invalid value '%s' for \"numberOfConstraints\" attribute.\n", xmlGetAttrval(constraints, "numberOfConstraints"));
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(*nconss >= 0);

   return SCIP_OKAY;
}

/** helper method to create and add a constraint (or a nonlinear objective constraint) */
static
SCIP_RETCODE createConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            linvars,            /**< array containing the linear variables (might be NULL) */
   SCIP_Real*            lincoefs,           /**< array containing the coefficients of the linear variables (might be NULL) */
   int                   nlinvars,           /**< the total number of linear variables */
   SCIP_VAR**            quadvars1,          /**< array containing the first variables of the quadratic terms (might be NULL) */
   SCIP_VAR**            quadvars2,          /**< array containing the second variables of the quadratic terms (might be NULL) */
   SCIP_Real*            quadcoefs,          /**< array containing the coefficients of the quadratic terms (might be NULL) */
   int                   nquadterms,         /**< the total number of quadratic terms */
   SCIP_EXPR*            nlexpr,             /**< the nonlinear part (might be NULL) */
   SCIP_Real             lhs,                /**< left-hand side */
   SCIP_Real             rhs,                /**< right-hand side */
   const char*           name,               /**< name of the constraint */
   SCIP_Bool             objcons,            /**< whether to add an objective constraints */
   SCIP_Bool             initialconss,       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss,       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamicrows         /**< should rows be added and removed dynamically to the LP? */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* objvar = NULL;

   assert(nlinvars >= 0);
   assert(nquadterms >= 0);

   /* create objective variable, if requested */
   if( objcons )
   {
      SCIP_CALL( SCIPcreateVar(scip, &objvar, "nlobjvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, objvar) );
   }

   /* linear constraint (can be empty) */
   if( nquadterms == 0 && nlexpr == NULL )
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name,
         nlinvars, linvars, lincoefs, lhs, rhs, initialconss,
         TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows, FALSE) );

      /* add objective variable, if requested */
      if( objcons )
      {
         assert(objvar != NULL);
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, objvar, -1.0) );
      }
   }
   /* nonlinear constraint */
   else
   {
      SCIP_EXPR* expr = NULL;
      SCIP_EXPR* varexpr = NULL;

      /* create variable expression for objvar */
      if( objcons )
      {
         SCIP_CALL( SCIPcreateExprVar(scip, &varexpr, objvar, NULL, NULL) );
      }

      /* check whether there is a quadratic part */
      if( nlinvars > 0 || nquadterms > 0 )
      {
         /* create quadratic expression; note that this is always a sum */
         SCIP_CALL( SCIPcreateExprQuadratic(scip, &expr, nlinvars, linvars, lincoefs,
            nquadterms, quadvars1, quadvars2, quadcoefs, NULL, NULL) );
         assert(SCIPisExprSum(scip, expr));

         /* add nonlinear expression as a child to expr */
         if( nlexpr != NULL )
         {
            SCIP_CALL( SCIPappendExprSumExpr(scip, expr, nlexpr, 1.0) );
         }

         /* add expression that represents the objective variable as a child to expr */
         if( varexpr != NULL )
         {
            SCIP_CALL( SCIPappendExprSumExpr(scip, expr, varexpr, -1.0) );
         }

         /* create nonlinear constraint */
         SCIP_CALL( SCIPcreateConsNonlinear(scip, &cons, name, expr, lhs, rhs,
            initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows) );

         /* release created expression */
         SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
      }

      /* there is no quadratic part but we might need to take care of the objective variable */
      else
      {
         assert(nlexpr != NULL);

         if( objcons )
         {
            SCIP_EXPR* sumexpr;
            SCIP_EXPR* children[2] = {nlexpr, varexpr};
            SCIP_Real coefs[2] = {1.0, -1.0};

            assert(varexpr != NULL);

            /* create sum expression */
            SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 2, children, coefs, 0.0, NULL, NULL) );

            /* create nonlinear constraint */
            SCIP_CALL( SCIPcreateConsNonlinear(scip, &cons, name, sumexpr, lhs, rhs,
               initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows) );

            /* release sum expression */
            SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
         }
         else
         {
            /* create nonlinear constraint */
            SCIP_CALL( SCIPcreateConsNonlinear(scip, &cons, name, nlexpr, lhs, rhs,
               initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows) );
         }
      }

      /* release variable expression */
      if( objcons )
      {
         assert(varexpr != NULL);
         SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );
      }
   }

   /* add and release constraint */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* release objective variable */
   if( objcons )
   {
      assert(objvar != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &objvar) );
   }

   return SCIP_OKAY;
}


/** reads constraint-specific information; creates and adds linear and nonlinear constraints based on the
 * information that have been collected by @ref readLinearCoefs, @ref readQuadraticCoefs, and @ref readNonlinearExprs
 */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   int                   nconss,             /**< total number of constraints */
   SCIP_VAR***           linvars,            /**< array containing for each constraint the linear variables */
   SCIP_Real**           lincoefs,           /**< array containing for each constraint the coefficients of the linear variables */
   int*                  nlinvars,           /**< array containing for each constraint the total number of linear variables */
   SCIP_VAR***           quadvars1,          /**< array containing for each constraint the first variables of the quadratic terms */
   SCIP_VAR***           quadvars2,          /**< array containing for each constraint the second variables of the quadratic terms */
   SCIP_Real**           quadcoefs,          /**< array containing for each constraint the coefficients of the quadratic terms */
   int*                  nquadterms,         /**< array containing for each constraint the total number of quadratic terms */
   SCIP_EXPR**           nlexprs,            /**< array containing for each constraint the nonlinear part */
   SCIP_Bool             initialconss,       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss,       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamicrows,        /**< should rows be added and removed dynamically to the LP? */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* constraints;
   const XML_NODE* consnode;
   const char* attrval;
   char name[SCIP_MAXSTRLEN];
   int c = 0;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(doingfine != NULL);
   assert(linvars != NULL);
   assert(lincoefs != NULL);
   assert(nlinvars != NULL);
   assert(quadvars1 != NULL);
   assert(quadvars2 != NULL);
   assert(quadcoefs != NULL);
   assert(nquadterms != NULL);
   assert(nlexprs != NULL);

   constraints = xmlFindNodeMaxdepth(datanode, "constraints", 0, 1);

   /* if no constraints, then nothing to do here */
   if( constraints == NULL )
      return SCIP_OKAY;

   /* read constraint names, lhs, rhs, constant */
   for( consnode = xmlFirstChild(constraints); consnode != NULL; consnode = xmlNextSibl(consnode) )
   {
      const char* consname;
      SCIP_Real conslhs;
      SCIP_Real consrhs;

      if( c == nconss )
      {
         SCIPerrorMessage("Expected %d constraints, but got at least %d many.\n", nconss, c+1);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* find constraint name */
      consname = xmlGetAttrval(consnode, "name");
      if( consname == NULL )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cons%d", c);
         consname = name;
      }

      /* check for mult attribute */
      attrval = xmlGetAttrval(consnode, "mult");
      if( attrval != NULL && strcmp(attrval, "1") != 0 )
      {
         SCIPerrorMessage("Constraint attribute 'mult' not supported (while parsing constraint <%s>).\n", consname);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* find constraint lower bound (=lhs) (default is -infinity) */
      attrval = xmlGetAttrval(consnode, "lb");
      if( attrval == NULL )
         conslhs = -SCIPinfinity(scip);
      else if( strcmp(attrval, "-INF") == 0 )
         conslhs = -SCIPinfinity(scip);
      else if( strcmp(attrval, "INF") == 0 )
         conslhs = SCIPinfinity(scip);
      else
      {
         conslhs = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("Error parsing constraint lower bound '%s' for constraint <%s>.\n", attrval, consname);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }

      /* find constraint upper bound (=rhs) (default is +infinity) */
      attrval = xmlGetAttrval(consnode, "ub");
      if( attrval == NULL )
         consrhs = SCIPinfinity(scip);
      else if( strcmp(attrval, "-INF") == 0 )
         consrhs = -SCIPinfinity(scip);
      else if( strcmp(attrval, "INF") == 0 )
         consrhs = SCIPinfinity(scip);
      else
      {
         consrhs = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("Error parsing constraint upper bound '%s' for constraint <%s>.\n", attrval, consname);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }

      /* find constraint constant (default is 0.0) and substract from lhs/rhs */
      attrval = xmlGetAttrval(consnode, "constant");
      if( attrval != NULL )
      {
         SCIP_Real consconstant;

         consconstant = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("Error parsing constraint constant '%s' for constraint <%s>.\n", attrval, consname);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
         if( conslhs > -SCIPinfinity(scip) )
            conslhs -= consconstant;
         if( consrhs <  SCIPinfinity(scip) )
            consrhs -= consconstant;
      }

      /* create, add, and release constraint */
      SCIP_CALL( createConstraint(scip, linvars[c], lincoefs[c], nlinvars[c],
         quadvars1[c], quadvars2[c], quadcoefs[c], nquadterms[c], nlexprs[c],
         conslhs, consrhs, consname, FALSE, initialconss, dynamicconss, dynamicrows) );

      ++c;
   }

   if( c != nconss )
   {
      SCIPerrorMessage("Got %d constraints, but expected %d many.\n", c, nconss);
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** reads mult and incr attributes of an OSiL node
 *
 *  if mult attribute is not present, then returns mult=1
 *  if incr attribute is not present, then returns incrint=0 and incrreal=0
 */
static
void readMultIncr(
   const XML_NODE*       node,               /**< XML node to read attributes from */
   int*                  mult,               /**< buffer to store mult */
   int*                  incrint,            /**< buffer to store incr as int, or NULL if no int expected */
   SCIP_Real*            incrreal,           /**< buffer to store incr as real, or NULL if no real expected */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const char* attrval;

   assert(node != NULL);
   assert(mult != NULL);
   assert(doingfine != NULL);

   *mult = 1;
   if( incrint != NULL )
      *incrint = 0;
   if( incrreal != NULL )
      *incrreal = 0.0;

   attrval = xmlGetAttrval(node, "mult");
   if( attrval == NULL )
      return;

   /* read "mult" attribute */
   *mult = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || *mult < 1 )
   {
      SCIPerrorMessage("Invalid value '%s' in \"mult\" attribute of node.\n", xmlGetAttrval(node, "mult"));
      *doingfine = FALSE;
      return;
   }

   if( *mult == 1 )
      return;

   /* read "incr" attribute */
   attrval = xmlGetAttrval(node, "incr");
   if( attrval == NULL )
      return;

   if( incrint != NULL )
   {
      *incrint = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' )
      {
         SCIPerrorMessage("Invalid value '%s' in \"incr\" attribute of node.\n", xmlGetAttrval(node, "incr"));
         *doingfine = FALSE;
         return;
      }
   }

   if( incrreal != NULL )
   {
      *incrreal = strtod(attrval, (char**)&attrval);
      if( *attrval != '\0' || !SCIPisFinite(*incrreal) )
      {
         SCIPerrorMessage("Invalid value '%s' in \"incr\" attribute of node.\n", xmlGetAttrval(node, "incr"));
         *doingfine = FALSE;
         return;
      }
   }
}

/** parse linear coefficients of constraints */
static
SCIP_RETCODE readLinearCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   int                   nconss,             /**< number of constraints */
   SCIP_VAR***           linvars,            /**< array to store for each constraint the linear variables */
   SCIP_Real**           lincoefs,           /**< array to store for each constraint the coefficients of the linear variables */
   int*                  nlinvars,           /**< array to store for each constraint the total number of linear variables */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* lincoef;
   const XML_NODE* startnode;
   const XML_NODE* idxnode;
   const XML_NODE* valnode;
   const XML_NODE* elnode;
   const char* attrval;
   SCIP_Bool rowmajor;
   int* start;
   int* idx;
   SCIP_Real* val;
   int nnz;
   int count;
   int mult;
   int incrint;
   SCIP_Real incrreal;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(doingfine != NULL);

   lincoef = xmlFindNodeMaxdepth(datanode, "linearConstraintCoefficients", 0, 1);

   if( lincoef == NULL )
      return SCIP_OKAY;

   /* get number of linear constraint coefficients */
   attrval = xmlGetAttrval(lincoef, "numberOfValues");
   if( attrval == NULL )
   {
      SCIPerrorMessage("Attribute \"numberOfValues\" not found for <linearConstraintCoefficients> node.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   nnz = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || nnz < 0 )
   {
      SCIPerrorMessage("Invalid value '%s' for \"numberOfValues\" attribute in <linearConstraintCoefficients> node.\n", xmlGetAttrval(lincoef, "numberOfValues"));
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(nnz >= 0);

   /* check for start, rowIdx, colIdx, and value nodes */
   startnode = xmlFindNodeMaxdepth(lincoef, "start", 0, 1);
   if( startnode == NULL )
   {
      SCIPerrorMessage("Node <start> not found inside <linearConstraintCoefficients> node.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   idxnode = xmlFindNodeMaxdepth(lincoef, "rowIdx", 0, 1);
   if( idxnode != NULL )
   {
      if( xmlFindNodeMaxdepth(lincoef, "colIdx", 0, 1) != NULL )
      {
         SCIPerrorMessage("Both <rowIdx> and <colIdx> found under <linearConstraintCoefficients> node.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      rowmajor = FALSE;
   }
   else
   {
      idxnode = xmlFindNodeMaxdepth(lincoef, "colIdx", 0, 1);
      if( idxnode == NULL )
      {
         SCIPerrorMessage("Both <rowIdx> and <colIdx> not found under <linearConstraintCoefficients> node.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      rowmajor = TRUE;
   }

   valnode = xmlFindNodeMaxdepth(lincoef, "value", 0, 1);
   if( valnode == NULL )
   {
      SCIPerrorMessage("<value> node not found under <linearConstraintCoefficients> node.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   start = NULL;
   idx = NULL;
   val = NULL;

   /* read row or column start indices */
   SCIP_CALL( SCIPallocBufferArray(scip, &start, (rowmajor ? nconss : nvars) + 1) );

   count = 0;
   for( elnode = xmlFirstChild(startnode); elnode != NULL; elnode = xmlNextSibl(elnode), ++count )
   {
      /* check for <el> node and read it's data */
      if( strcmp(xmlGetName(elnode), "el") != 0 )
      {
         SCIPerrorMessage("Expected <el> node under <start> node in <linearConstraintCoefficients>, but got '%s'.\n", xmlGetName(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( count >= (rowmajor ? nconss : nvars) + 1 )
      {
         SCIPerrorMessage("Too many elements under <start> node in <linearConstraintCoefficients>, expected %d many, got at least %d.\n", (rowmajor ? nconss : nvars) + 1, count + 1);
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( xmlFirstChild(elnode) == NULL || xmlGetData(xmlFirstChild(elnode)) == NULL )
      {
         SCIPerrorMessage("No data in <el> node in <linearConstraintCoefficients>.\n");
         *doingfine = FALSE;
         goto CLEANUP;
      }

      start[count] = (int)strtol(xmlGetData(xmlFirstChild(elnode)), (char**)&attrval, 10);

      if( *attrval != '\0' || start[count] < 0 || (start[count] > nnz) )
      {
         SCIPerrorMessage("Invalid value '%s' in <el> node under <start> node in <linearConstraintCoefficients>.\n", xmlGetData(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }

      /* add additional start-indices according to mult and incr attributes */
      readMultIncr(elnode, &mult, &incrint, NULL, doingfine);
      if( !*doingfine )
         goto CLEANUP;

      for( --mult; mult > 0; --mult )
      {
         ++count;
         if( count >= (rowmajor ? nconss : nvars) + 1 )
         {
            SCIPerrorMessage("Too many elements under <start> node in <linearConstraintCoefficients>, expected %d many, got at least %d.\n", (rowmajor ? nconss : nvars) + 1, count + 1);
            *doingfine = FALSE;
            goto CLEANUP;
         }
         start[count] = start[count-1] + incrint;
      }
   }
   if( count != (rowmajor ? nconss : nvars) + 1 )
   {
      SCIPerrorMessage("Got only %d <start> entries in <linearConstraintCoefficients>, but expected %d many.\n", count, (rowmajor ? nconss : nvars) + 1);
      *doingfine = FALSE;
      goto CLEANUP;
   }

   /* read row or column indices */
   SCIP_CALL( SCIPallocBufferArray(scip, &idx, nnz) );

   count = 0;
   for( elnode = xmlFirstChild(idxnode); elnode != NULL; elnode = xmlNextSibl(elnode), ++count )
   {
      /* check for <el> node and read it's data */
      if( strcmp(xmlGetName(elnode), "el") != 0 )
      {
         SCIPerrorMessage("Expected <el> node under <%s> node in <linearConstraintCoefficients>, but got '%s'.\n", rowmajor ? "colIdx" : "rowIdx", xmlGetName(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( count >= nnz )
      {
         SCIPerrorMessage("Too many elements under <%s> node in <linearConstraintCoefficients>, expected %d many, but got at least %d.\n", rowmajor ? "colIdx" : "rowIdx", nnz, count + 1);
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( xmlFirstChild(elnode) == NULL || xmlGetData(xmlFirstChild(elnode)) == NULL )
      {
         SCIPerrorMessage("No data in <el> node under <%s> node in <linearConstraintCoefficients>.\n", rowmajor ? "colIdx" : "rowIdx");
         *doingfine = FALSE;
         goto CLEANUP;
      }

      idx[count] = (int)strtol(xmlGetData(xmlFirstChild(elnode)), (char**)&attrval, 10);

      if( *attrval != '\0' || idx[count] < 0 || (idx[count] >= (rowmajor ? nvars : nconss)) )
      {
         SCIPerrorMessage("Invalid value '%s' in <el> node under <%s> node in <linearConstraintCoefficients>.\n", xmlGetData(elnode), rowmajor ? "colIdx" : "rowIdx");
         *doingfine = FALSE;
         goto CLEANUP;
      }

      /* add additional indices according to mult and incr attributes */
      readMultIncr(elnode, &mult, &incrint, NULL, doingfine);
      if( !*doingfine )
         goto CLEANUP;

      for( --mult; mult > 0; --mult )
      {
         ++count;
         if( count >= nnz )
         {
            SCIPerrorMessage("Too many elements under <%s> node in <linearConstraintCoefficients>, expected %d many, got at least %d.\n", rowmajor ? "colIdx" : "rowIdx", nnz, count + 1);
            *doingfine = FALSE;
            goto CLEANUP;
         }
         idx[count] = idx[count-1] + incrint;
      }
   }
   if( count != nnz )
   {
      SCIPerrorMessage("Got only %d entries in <%s> node in <linearConstraintCoefficients>, expected %d many.\n", count, rowmajor ? "colIdx" : "rowIdx", nnz);
      *doingfine = FALSE;
      goto CLEANUP;
   }

   /* read coefficient values */
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnz) );

   count = 0;
   for( elnode = xmlFirstChild(valnode); elnode != NULL; elnode = xmlNextSibl(elnode), ++count )
   {
      /* check for <el> node and read it's data */
      if( strcmp(xmlGetName(elnode), "el") != 0 )
      {
         SCIPerrorMessage("Expected <el> node under <value> node in <linearConstraintCoefficients>, but got '%s'.\n", xmlGetName(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( count >= nnz )
      {
         SCIPerrorMessage("Too many elements under <value> node in <linearConstraintCoefficients>, expected %d many, got at least %d.\n", nnz, count + 1);
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( xmlFirstChild(elnode) == NULL || xmlGetData(xmlFirstChild(elnode)) == NULL )
      {
         SCIPerrorMessage("No data in <el> node under <value> node in <linearConstraintCoefficients>.\n");
         *doingfine = FALSE;
         goto CLEANUP;
      }

      val[count] = strtod(xmlGetData(xmlFirstChild(elnode)), (char**)&attrval);

      if( *attrval != '\0' || !SCIPisFinite(val[count]) )
      {
         SCIPerrorMessage("Invalid value '%s' in <el> node under <value> node in <linearConstraintCoefficients>.\n", xmlGetData(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }

      /* add additional values according to mult and incr attributes */
      readMultIncr(elnode, &mult, NULL, &incrreal, doingfine);
      if( !*doingfine )
         goto CLEANUP;

      for( --mult; mult > 0; --mult )
      {
         ++count;
         if( count >= nnz )
         {
            SCIPerrorMessage("Too many elements under <value> node in <linearConstraintCoefficients>, expected %d many, got at least %d.\n", nnz, count + 1);
            *doingfine = FALSE;
            goto CLEANUP;
         }
         val[count] = val[count-1] + incrreal;
      }
   }
   if( count != nnz )
   {
      SCIPerrorMessage("Got only %d entries under <value> node in <linearConstraintCoefficients>, expected %d many.\n", count, nnz);
      *doingfine = FALSE;
      goto CLEANUP;
   }

   /* add coefficients to linear constraints */
   if( rowmajor )
   {
      int row;
      int pos;
      for( row = 0; row < nconss; ++row )
      {
         int nterms;

         /* these asserts were checked above */
         assert(start[row] >= 0);
         assert(start[row+1] >= 0);
         assert(start[row] <= nnz);
         assert(start[row+1] <= nnz);

         assert(linvars[row] == NULL);
         assert(lincoefs[row] == NULL);
         assert(nlinvars[row] == 0);

         nterms = start[row+1] - start[row];
         SCIP_CALL( SCIPallocBufferArray(scip, &linvars[row], nterms) );
         SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs[row], nterms) );

         for( pos = start[row]; pos < start[row+1]; ++pos )
         {
            /* these asserts were checked above */
            assert(pos >= 0);
            assert(pos < nnz);
            assert(idx[pos] >= 0);
            assert(idx[pos] < nvars);

            linvars[row][nlinvars[row]] = vars[idx[pos]];
            lincoefs[row][nlinvars[row]] = val[pos];
            ++(nlinvars[row]);
         }
         assert(nlinvars[row] == nterms);
      }
   }
   else
   {
      int col;
      int pos;
      int k;

      /* allocate memory for the coefficients in iteration k=0; in k=1 fill in the data */
      for( k = 0; k < 2; ++k )
      {
         for( col = 0; col < nvars; ++col )
         {
            /* these asserts were checked above */
            assert(start[col] >= 0);
            assert(start[col+1] >= 0);
            assert(start[col] <= nnz);
            assert(start[col+1] <= nnz);
            for( pos = start[col]; pos < start[col+1]; ++pos )
            {
               int considx = idx[pos];

               /* these asserts were checked above */
               assert(pos >= 0);
               assert(pos < nnz);
               assert(considx >= 0);
               assert(considx < nconss);

               if( k == 0 )
               {
                  ++(nlinvars[considx]);
               }
               else
               {
                  linvars[considx][nlinvars[considx]] = vars[col];
                  lincoefs[considx][nlinvars[considx]] = val[pos];
                  ++(nlinvars[considx]);
               }
            }
         }

         /* allocate memory to store the linear coefficients for each constraint after the first iteration */
         if( k == 0 )
         {
            int c;

            for( c = 0; c < nconss; ++c )
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &linvars[c], nlinvars[c]) );
               SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs[c], nlinvars[c]) );

               /* reset nlinvars[c] so it can be used for iteration k=1 */
               nlinvars[c] = 0;
            }
         }
      }
   }

 CLEANUP:
   SCIPfreeBufferArrayNull(scip, &val);
   SCIPfreeBufferArrayNull(scip, &idx);
   SCIPfreeBufferArrayNull(scip, &start);

   return SCIP_OKAY;
}

/** read quadratic coefficients of constraints and objective */
static
SCIP_RETCODE readQuadraticCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   int                   nconss,             /**< number of constraints */
   SCIP_VAR***           quadvars1,          /**< array to store for each constraint the first variables of the quadratic terms */
   SCIP_VAR***           quadvars2,          /**< array to store for each constraint the second variables of the quadratic terms */
   SCIP_Real**           quadcoefs,          /**< array to store for each constraint the coefficients of the quadratic terms */
   int*                  nquadterms,         /**< array to store for each constraint the total number of quadratic terms */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* quadcoef;
   const XML_NODE* qterm;
   const char* attrval;
   int* termssize;
   SCIP_Real coef;
   int nqterms;
   int count;
   int considx;
   int varidx1;
   int varidx2;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(quadvars1 != NULL);
   assert(quadvars2 != NULL);
   assert(quadcoefs != NULL);
   assert(nquadterms != NULL);
   assert(doingfine != NULL);

   quadcoef = xmlFindNodeMaxdepth(datanode, "quadraticCoefficients", 0, 1);

   if( quadcoef == NULL )
      return SCIP_OKAY;

   /* read number of quadratic terms */
   attrval = xmlGetAttrval(quadcoef, "numberOfQuadraticTerms");
   if( attrval == NULL )
   {
      SCIPerrorMessage("Attribute \"numberOfQuadraticTerms\" not found for <quadraticCoefficients> node.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   nqterms = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || nqterms < 0 )
   {
      SCIPerrorMessage("Invalid value '%s' for \"numberOfQuadraticTerms\" attribute of <quadraticCoefficients> node.\n", xmlGetAttrval(quadcoef, "numberOfQuadraticTerms"));
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(nqterms >= 0);

   if( nqterms == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   SCIP_CALL( SCIPallocClearBufferArray(scip, &termssize, nconss + 1) );

   count = 0;
   for( qterm = xmlFirstChild(quadcoef); qterm != NULL; qterm = xmlNextSibl(qterm), ++count )
   {
      /* check for qterm node */
      if( strcmp(xmlGetName(qterm), "qTerm") != 0 )
      {
         SCIPerrorMessage("Expected <qTerm> node under <quadraticCoefficients> node, but got <%s>\n", xmlGetName(qterm));
         *doingfine = FALSE;
         goto TERMINATE;
      }
      if( count >= nqterms )
      {
         SCIPerrorMessage("Too many quadratic terms under <quadraticCoefficients> node, expected %d many, but got at least %d.\n", nqterms, count + 1);
         *doingfine = FALSE;
         goto TERMINATE;
      }

      /* get constraint index, or -1 for objective */
      attrval = xmlGetAttrval(qterm, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing \"idx\" attribute in %d'th <qTerm> node under <quadraticCoefficients> node.\n", count);
         *doingfine = FALSE;
         goto TERMINATE;
      }

      considx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || considx < -1 || considx >= nconss )
      {
         SCIPerrorMessage("Invalid value '%s' in \"idx\" attribute of %d'th <qTerm> node under <quadraticCoefficients> node.\n", xmlGetAttrval(qterm, "idx"), count);
         *doingfine = FALSE;
         goto TERMINATE;
      }

      /* get index of first variable */
      attrval = xmlGetAttrval(qterm, "idxOne");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing \"idxOne\" attribute in %d'th <qTerm> node under <quadraticCoefficients> node.\n", count);
         *doingfine = FALSE;
         goto TERMINATE;
      }

      varidx1 = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || varidx1 < 0 || varidx1 >= nvars )
      {
         SCIPerrorMessage("Invalid value '%s' in \"idxOne\" attribute of %d'th <qTerm> node under <quadraticCoefficients> node.\n", xmlGetAttrval(qterm, "idxOne"), count);
         *doingfine = FALSE;
         goto TERMINATE;
      }

      /* get index of second variable */
      attrval = xmlGetAttrval(qterm, "idxTwo");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing \"idxTwo\" attribute in %d'th <qTerm> node under <quadraticCoefficients> node.\n", count);
         *doingfine = FALSE;
         goto TERMINATE;
      }

      varidx2 = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || varidx2 < 0 || varidx2 >= nvars )
      {
         SCIPerrorMessage("Invalid value '%s' in \"idxTwo\" attribute of %d'th <qTerm> node under <quadraticCoefficients> node.\n", xmlGetAttrval(qterm, "idxTwo"), count);
         *doingfine = FALSE;
         goto TERMINATE;
      }

      /* get (optional) coefficient of quadratic term */
      attrval = xmlGetAttrval(qterm, "coef");
      if( attrval != NULL )
      {
         coef = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' || (coef != coef) )  /*lint !e777*/
         {
            SCIPerrorMessage("Invalid value '%s' in \"coef\" attribute of %d'th <qTerm> node under <quadraticCoefficients> node.\n", xmlGetAttrval(qterm, "coef"), count);
            *doingfine = FALSE;
            goto TERMINATE;
         }
      }
      else
      {
         /* default is 1.0 according to specification */
         coef = 1.0;
      }

      /* skip zero coefficients */
      if( coef == 0.0 )
         continue;

      /* put objective at end of array */
      if( considx == -1 )
         considx = nconss;

      if( nquadterms[considx] + 1 > termssize[considx] )
      {
         termssize[considx] = SCIPcalcMemGrowSize(scip, nquadterms[considx] + 1);
         SCIP_CALL( SCIPreallocBufferArray(scip, &quadvars1[considx], termssize[considx]) );  /*lint !e866*/
         SCIP_CALL( SCIPreallocBufferArray(scip, &quadvars2[considx], termssize[considx]) );  /*lint !e866*/
         SCIP_CALL( SCIPreallocBufferArray(scip, &quadcoefs[considx], termssize[considx]) );  /*lint !e866*/
      }

      quadvars1[considx][nquadterms[considx]] = vars[varidx1];
      quadvars2[considx][nquadterms[considx]] = vars[varidx2];
      quadcoefs[considx][nquadterms[considx]] = coef;
      ++nquadterms[considx];
   }

   if( count != nqterms )
   {
      SCIPerrorMessage("Got only %d quadratic terms under <quadraticCoefficients> node, but expected %d many.\n", count, nqterms);
      *doingfine = FALSE;
      goto TERMINATE;
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &termssize);

   return SCIP_OKAY;
}

/** transforms OSnL expression tree into SCIP expression */
static
SCIP_RETCODE readExpression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to created expression */
   const XML_NODE*       node,               /**< root node of expression to be read */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< total number of variables in problem */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const char* exprname;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(node != NULL);
   assert(vars != NULL);
   assert(doingfine != NULL);

   exprname = xmlGetName(node);
   assert(exprname != NULL);

   *expr = NULL;

   /* zero argument operands */
   if( strcmp(exprname, "variable") == 0 )
   {
      const char* attrval;
      SCIP_Real coef;
      int idx;

      /* read variable index */
      attrval = xmlGetAttrval(node, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Attribute \"idx\" required for <variable> node in nonlinear expression\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      idx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || idx < 0 || idx >= nvars )
      {
         SCIPerrorMessage("Invalid value '%s' in \"idx\" attribute of <variable> node in nonlinear expression.\n", xmlGetAttrval(node, "idx"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* read variable coefficient */
      attrval = xmlGetAttrval(node, "coef");
      if( attrval != NULL )
      {
         coef = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' || !SCIPisFinite(coef) )
         {
            SCIPerrorMessage("Invalid value '%s' in \"coef\" attribute of <variable> node in nonlinear expression.\n", xmlGetAttrval(node, "coef"));
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         coef = 1.0;
      }

      /* create variable expression */
      SCIP_CALL( SCIPcreateExprVar(scip, expr, vars[idx], NULL, NULL) );

      /* create a sum if the coefficient != 1 */
      if( coef != 1.0 )
      {
         SCIP_EXPR* sumexpr;

         SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 1, expr, &coef, 0.0, NULL, NULL) );

         /* release the variable expression and store the sum */
         SCIP_CALL( SCIPreleaseExpr(scip, expr) );
         *expr = sumexpr;
      }

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "number") == 0 )
   {
      const char* attrval;
      SCIP_Real val;

      attrval = xmlGetAttrval(node, "type");
      if( attrval != NULL && (strcmp(attrval, "real") != 0) )
      {
         SCIPerrorMessage("Type '%s' for <number> node in nonlinear expression not supported.\n", attrval);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      attrval = xmlGetAttrval(node, "value");
      if( attrval != NULL )
      {
         val = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' || !SCIPisFinite(val) )
         {
            SCIPerrorMessage("Invalid value '%s' in \"value\" attribute of <number> node in nonlinear expression.\n", xmlGetAttrval(node, "value"));
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         /* according to OSnL.xsd, the value attribute is optional
          * I guess the default is the empty string, which should correspond to 0.0
          */
         val = 0.0;
      }

      /* create constant expression */
      SCIP_CALL( SCIPcreateExprValue(scip, expr, val, NULL, NULL) );

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "PI") == 0 )
   {
      /* create constant expression with PI value */
      SCIP_CALL( SCIPcreateExprValue(scip, expr, M_PI, NULL, NULL) );

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "E") == 0 )
   {
      /* create constant expression with PI value */
      SCIP_CALL( SCIPcreateExprValue(scip, expr, M_E, NULL, NULL) );

      return SCIP_OKAY;
   }

   /* single argument operands */
   if( strcmp(exprname, "negate") == 0 ||
      strcmp(exprname, "abs") == 0 ||
      strcmp(exprname, "squareRoot") == 0 ||
      strcmp(exprname, "sqrt") == 0 ||
      strcmp(exprname, "square") == 0 ||
      strcmp(exprname, "exp") == 0 ||
      strcmp(exprname, "ln") == 0 ||
      strcmp(exprname, "log10") == 0 ||
      strcmp(exprname, "sin") == 0 ||
      strcmp(exprname, "cos") == 0 ||
      strcmp(exprname, "erf") == 0
      )
   {
      SCIP_EXPR* arg;

      /* check number of children */
      if( xmlFirstChild(node) == NULL || xmlNextSibl(xmlFirstChild(node)) != NULL )
      {
         SCIPerrorMessage("Expected exactly one child in <%s> node in nonlinear expression\n", exprname);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* read child expression */
      SCIP_CALL( readExpression(scip, &arg, xmlFirstChild(node), vars, nvars, doingfine) );

      if( !*doingfine )
         return SCIP_OKAY;
      assert(arg != NULL);

      /* create SCIP expression according to expression name */
      if( strcmp(exprname, "negate") == 0 )
      {
         SCIP_Real minusone;

         minusone = -1.0;

         SCIP_CALL( SCIPcreateExprSum(scip, expr, 1, &arg, &minusone, 0.0, NULL, NULL) );
      }
      else if( strcmp(exprname, "abs") == 0 )
      {
         SCIP_CALL( SCIPcreateExprAbs(scip, expr, arg, NULL, NULL) );
      }
      else if( strcmp(exprname, "squareRoot") == 0 || strcmp(exprname, "sqrt") == 0 )
      {
         SCIP_CALL( SCIPcreateExprPow(scip, expr, arg, 0.5, NULL, NULL) );
      }
      else if( strcmp(exprname, "square") == 0 )
      {
         SCIP_CALL( SCIPcreateExprPow(scip, expr, arg, 2.0, NULL, NULL) );
      }
      else if( strcmp(exprname, "exp") == 0 )
      {
         SCIP_CALL( SCIPcreateExprExp(scip, expr, arg, NULL, NULL) );
      }
      else if( strcmp(exprname, "ln") == 0 )
      {
         SCIP_CALL( SCIPcreateExprLog(scip, expr, arg, NULL, NULL) );
      }
      else if( strcmp(exprname, "log10") == 0 )
      {
         SCIP_EXPR* logexpr;
         SCIP_Real coef = 1.0/log(10.0);

         SCIP_CALL( SCIPcreateExprLog(scip, &logexpr, arg, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, expr, 1, &logexpr, &coef, 0.0, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &logexpr) );
      }
      else if( strcmp(exprname, "sin") == 0 )
      {
         SCIP_CALL( SCIPcreateExprSin(scip, expr, arg, NULL, NULL) );
      }
      else if( strcmp(exprname, "cos") == 0 )
      {
         SCIP_CALL( SCIPcreateExprCos(scip, expr, arg, NULL, NULL) );
      }
      else if( strcmp(exprname, "erf") == 0 )
      {
         SCIPwarningMessage(scip, "Danger! You're entering a construction area. Implementation of support for 'erf' is incomplete.\n");
         SCIP_CALL( SCIPcreateExprErf(scip, expr, arg, NULL, NULL) );
      }

      /* release argument expression */
      SCIP_CALL( SCIPreleaseExpr(scip, &arg) );

      return SCIP_OKAY;
   }

   /* two argument operands */
   if( strcmp(exprname, "plus") == 0 ||
      strcmp(exprname, "minus") == 0 ||
      strcmp(exprname, "times") == 0 ||
      strcmp(exprname, "divide") == 0 ||
      strcmp(exprname, "power") == 0 ||
      strcmp(exprname, "signpower") == 0 ||
      strcmp(exprname, "log") == 0
     )
   {
      SCIP_EXPR* args[2] = {NULL, NULL};

      /* check number of children */
      if( xmlFirstChild(node) == NULL ||
         xmlNextSibl(xmlFirstChild(node)) == NULL ||
         xmlNextSibl(xmlNextSibl(xmlFirstChild(node))) != NULL )
      {
         SCIPerrorMessage("Expected exactly two children in <%s> node in nonlinear expression.\n", exprname);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* read first child expression */
      SCIP_CALL( readExpression(scip, &args[0], xmlFirstChild(node), vars, nvars, doingfine) );
      if( !*doingfine )
         goto TERMINATE_TWO_ARGS;
      assert(args[0] != NULL);

      /* read second child expression */
      SCIP_CALL( readExpression(scip, &args[1], xmlNextSibl(xmlFirstChild(node)), vars, nvars, doingfine) );
      if( !*doingfine )
         goto TERMINATE_TWO_ARGS;
      assert(args[1] != NULL);

      if( strcmp(exprname, "plus") == 0 )
      {
         SCIP_CALL( SCIPcreateExprSum(scip, expr, 2, args, NULL, 0.0, NULL, NULL) );
      }
      else if( strcmp(exprname, "minus") == 0 )
      {
         SCIP_Real coefs[2] = {1.0, -1.0};
         SCIP_CALL( SCIPcreateExprSum(scip, expr, 2, args, coefs, 0.0, NULL, NULL) );
      }
      else if( strcmp(exprname, "times") == 0 )
      {
         SCIP_CALL( SCIPcreateExprProduct(scip, expr, 2, args, 1.0, NULL, NULL) );
      }
      else if( strcmp(exprname, "divide") == 0 )
      {
         SCIP_EXPR* tmp[2];
         SCIP_EXPR* powexpr;

         SCIP_CALL( SCIPcreateExprPow(scip, &powexpr, args[1], -1.0, NULL, NULL) );
         tmp[0] = args[0];
         tmp[1] = powexpr;
         SCIP_CALL( SCIPcreateExprProduct(scip, expr, 2, tmp, 1.0, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &powexpr) );
      }
      else if( strcmp(exprname, "power") == 0 )
      {
         /* case 1: expr^number */
         if( SCIPisExprValue(scip, args[1]) )
         {
            SCIP_CALL( SCIPcreateExprPow(scip, expr, args[0], SCIPgetValueExprValue(args[1]), NULL, NULL) );
         }
         /* case 2: number^expr = exp(arg2 * ln(number)) */
         else if( SCIPisExprValue(scip, args[0]) )
         {
            SCIP_Real value = SCIPgetValueExprValue(args[0]);

            if( value <= 0.0 )
            {
               SCIPerrorMessage("Negative base in <power> node with nonconstant exponent not allowed in nonlinear expression.\n");
               *doingfine = FALSE;
               goto TERMINATE_TWO_ARGS;
            }
            else
            {
               SCIP_EXPR* sumexpr;
               SCIP_Real coef = log(value);

               SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 1, &args[1], &coef, 0.0, NULL, NULL) );
               SCIP_CALL( SCIPcreateExprExp(scip, expr, sumexpr, NULL, NULL) );
               SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
            }
         }
         /* case 3: arg1^arg2 is exp(arg2 * ln(arg1)) */
         else
         {
            SCIP_EXPR* logexpr;
            SCIP_EXPR* prodexpr;
            SCIP_EXPR* tmp[2];

            SCIP_CALL( SCIPcreateExprLog(scip, &logexpr, args[0], NULL, NULL) );
            tmp[0] = args[1];
            tmp[1] = logexpr;
            SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr, 2, tmp, 1.0, NULL, NULL) );
            SCIP_CALL( SCIPcreateExprExp(scip, expr, prodexpr, NULL, NULL) );
            SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
            SCIP_CALL( SCIPreleaseExpr(scip, &logexpr) );
         }
      }
      else if( strcmp(exprname, "signpower") == 0 )
      {
         /* signpower(expr,number) with number > 1 is the only one we can handle */
         if( !SCIPisExprValue(scip, args[1]) )
         {
            SCIPerrorMessage("Signpower only supported for constant exponents.\n");
            *doingfine = FALSE;
            goto TERMINATE_TWO_ARGS;
         }
         if( SCIPgetValueExprValue(args[1]) <= 1.0 )
         {
            SCIPerrorMessage("Signpower only supported for exponents > 1, but got %g.\n",
               SCIPgetValueExprValue(args[1]));
            *doingfine = FALSE;
            goto TERMINATE_TWO_ARGS;
         }

         SCIP_CALL( SCIPcreateExprSignpower(scip, expr, args[0], SCIPgetValueExprValue(args[1]), NULL, NULL) );
      }
      /* logarithm of arg2 w.r.t. base arg1 = ln(arg2) / ln(arg1) */
      else if( strcmp(exprname, "log") == 0 )
      {
         SCIP_EXPR* logexpr0;
         SCIP_EXPR* logexpr1;
         SCIP_EXPR* powexpr;
         SCIP_EXPR* tmp[2];

         /* logarithm of arg2 w.r.t. base arg1 = ln(arg2) / ln(arg1) = ln(arg2) * pow(ln(arg1),-1) */
         SCIP_CALL( SCIPcreateExprLog(scip, &logexpr0, args[0], NULL, NULL) );
         SCIP_CALL( SCIPcreateExprLog(scip, &logexpr1, args[1], NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &powexpr, logexpr0, -1.0, NULL, NULL) );
         tmp[0] = logexpr1;
         tmp[1] = powexpr;
         SCIP_CALL( SCIPcreateExprProduct(scip, expr, 2, tmp, 1.0, NULL, NULL) );

         SCIP_CALL( SCIPreleaseExpr(scip, &powexpr) );
         SCIP_CALL( SCIPreleaseExpr(scip, &logexpr1) );
         SCIP_CALL( SCIPreleaseExpr(scip, &logexpr0) );
      }
      else if( strcmp(exprname, "min") == 0 )
      {
         SCIPerrorMessage("min expressions are not supported\n");
         *doingfine = FALSE;
         goto TERMINATE_TWO_ARGS;
      }
      else /* if( strcmp(exprname, "max") == 0 ) */
      {
         assert(strcmp(exprname, "max") == 0);

         SCIPerrorMessage("max expressions are not supported\n");
         *doingfine = FALSE;
         goto TERMINATE_TWO_ARGS;
      }

TERMINATE_TWO_ARGS:

      /* release first and second argument expression */
      if( args[0] != NULL )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &args[0]) );
      }
      if( args[1] != NULL )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &args[1]) );
      }

      return SCIP_OKAY;
   }

   /* arbitrary argument operands */
   if( strcmp(exprname, "sum") == 0 || strcmp(exprname, "product") == 0 )
   {
      const XML_NODE* argnode;
      SCIP_EXPR** args;
      int nargs;
      int argssize;
      int i;

      /* a sum or product w.r.t. 0 arguments is constant */
      if( xmlFirstChild(node) == NULL )
      {
         SCIP_CALL( SCIPcreateExprValue(scip, expr, (strcmp(exprname, "sum") == 0) ? 0.0 : 1.0, NULL, NULL) );

         return SCIP_OKAY;
      }

      /* read all child expressions */
      argssize = 5;
      SCIP_CALL( SCIPallocBufferArray(scip, &args, argssize) );

      nargs = 0;
      for( argnode = xmlFirstChild(node); argnode != NULL; argnode = xmlNextSibl(argnode), ++nargs )
      {
         if( nargs >= argssize )
         {
            argssize = SCIPcalcMemGrowSize(scip, nargs + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &args, argssize) );
         }
         assert(nargs < argssize);

         SCIP_CALL( readExpression(scip, &args[nargs], argnode, vars, nvars, doingfine) );
         if( !*doingfine )
         {
            assert(args[nargs] == NULL);
            break;
         }
      }

      if( *doingfine )
      {
         switch( nargs )
         {
            case 0:
            {
               SCIP_CALL( SCIPcreateExprValue(scip, expr, (strcmp(exprname, "sum") == 0) ? 0.0 : 1.0, NULL, NULL) );
               break;
            }
            case 1:
            {
               *expr = args[0];
               /* capture expression here because args[0] will be released at the end */
               SCIPcaptureExpr(*expr);
               break;
            }

            default:
            {
               /* create sum or product expression */
               if( strcmp(exprname, "sum") == 0 )
               {
                  SCIP_CALL( SCIPcreateExprSum(scip, expr, nargs, args, NULL, 0.0, NULL, NULL) );
               }
               else
               {
                  SCIP_CALL( SCIPcreateExprProduct(scip, expr, nargs, args, 1.0, NULL, NULL) );
               }

               break;
            }
         }
      }

      /* release argument expressions */
      for( i = 0; i < nargs; ++i )
      {
         assert(args[i] != NULL);
         SCIP_CALL( SCIPreleaseExpr(scip, &args[i]) );
      }

      SCIPfreeBufferArray(scip, &args);

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "min") == 0 || strcmp(exprname, "max") == 0 )
   {
      SCIPerrorMessage("min or max expressions are not supported\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   if( strcmp(exprname, "quadratic") == 0 )
   {
      const char* attrval;
      const XML_NODE* qterm;
      SCIP_VAR** quadvars1;
      SCIP_VAR** quadvars2;
      SCIP_Real* quadcoefs;
      int nquadelems;
      int quadelemssize;
      int idx;

      quadelemssize = 5;
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, quadelemssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, quadelemssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, quadelemssize) );
      nquadelems = 0;

      /* read quadratic terms */
      for( qterm = xmlFirstChild(node); qterm != NULL; qterm = xmlNextSibl(qterm), ++nquadelems )
      {
         /* check for qpTerm node */
         if( strcmp(xmlGetName(qterm), "qpTerm") != 0 )
         {
            SCIPerrorMessage("Unexpected <%s> node under <quadratic> node in nonlinear expression, expected <qpTerm>.\n", xmlGetName(qterm));
            *doingfine = FALSE;
            return SCIP_OKAY;
         }

         if( nquadelems >= quadelemssize )
         {
            quadelemssize = SCIPcalcMemGrowSize(scip, nquadelems + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &quadvars1, quadelemssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &quadvars2, quadelemssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &quadcoefs, quadelemssize) );
         }
         assert(quadelemssize > nquadelems);

         /* get index of first variable */
         attrval = xmlGetAttrval(qterm, "idxOne");
         if( attrval == NULL )
         {
            SCIPerrorMessage("Missing \"idxOne\" attribute in %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", nquadelems);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }

         idx = (int)strtol(attrval, (char**)&attrval, 10);
         if( *attrval != '\0' || idx < 0 || idx >= nvars )
         {
            SCIPerrorMessage("Invalid value '%s' for \"idxOne\" attribute of %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", xmlGetAttrval(qterm, "idxOne"), nquadelems);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
         quadvars1[nquadelems] = vars[idx];

         /* get index of second variable */
         attrval = xmlGetAttrval(qterm, "idxTwo");
         if( attrval == NULL )
         {
            SCIPerrorMessage("Missing \"idxTwo\" attribute in %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", nquadelems);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }

         idx = (int)strtol(attrval, (char**)&attrval, 10);
         if( *attrval != '\0' || idx < 0 || idx >= nvars )
         {
            SCIPerrorMessage("Invalid value '%s' for \"idxTwo\" attribute of %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", xmlGetAttrval(qterm, "idxTwo"), nquadelems);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
         quadvars2[nquadelems] = vars[idx];

         /* get coefficient */
         attrval = xmlGetAttrval(qterm, "coef");
         if( attrval != NULL )
         {
            quadcoefs[nquadelems] = strtod(attrval, (char**)&attrval);
            if( *attrval != '\0' || !SCIPisFinite(quadcoefs[nquadelems]) )  /*lint !e777*/
            {
               SCIPerrorMessage("Invalid value '%s' for \"coef\" attribute of %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", xmlGetAttrval(qterm, "coef"), nquadelems);
               *doingfine = FALSE;
               return SCIP_OKAY;
            }
         }
         else
         {
            quadcoefs[nquadelems] = 1.0;
         }
      }

      /* create quadratic expression */
      SCIP_CALL( SCIPcreateExprQuadratic(scip, expr, 0, NULL, NULL, nquadelems, quadvars1, quadvars2, quadcoefs, NULL, NULL) );
   }

   SCIPerrorMessage("Expression operand <%s> in nonlinear expression not supported by SCIP so far.\n", exprname);
   *doingfine = FALSE;

   return SCIP_OKAY;
}


/** read nonlinear expressions of constraints and objective */
static
SCIP_RETCODE readNonlinearExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   int                   nconss,             /**< number of constraints */
   SCIP_EXPR**           exprs,              /**< array to store for each constraint a nonlinear expression */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* nlexprs;
   const XML_NODE* nlexpr;
   const char* attrval;
   int nnlexprs;
   int count;
   int considx;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(exprs != NULL);
   assert(doingfine != NULL);

   nlexprs = xmlFindNodeMaxdepth(datanode, "nonlinearExpressions", 0, 1);

   if( nlexprs == NULL )
      return SCIP_OKAY;

   /* get number of nonlinear expressions */
   attrval = xmlGetAttrval(nlexprs, "numberOfNonlinearExpressions");
   if( attrval == NULL )
   {
      SCIPerrorMessage("Attribute \"numberOfNonlinearExpressions\" in <nonlinearExpressions> node not found.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   nnlexprs = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || nnlexprs < 0 )
   {
      SCIPerrorMessage("Invalid value '%s' for \"numberOfNonlinearExpressions\" attribute in <nonlinearExpressions>.\n", xmlGetAttrval(nlexprs, "numberOfNonlinearExpressions"));
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(nnlexprs >= 0);

   /* read nonlinear expressions and store in constraints */
   count = 0;
   for( nlexpr = xmlFirstChild(nlexprs); nlexpr != NULL; nlexpr = xmlNextSibl(nlexpr), ++count )
   {
      if( strcmp(xmlGetName(nlexpr), "nl") != 0 )
      {
         SCIPerrorMessage("Expected <nl> node under <nonlinearExpressions> node, but got '%s'.\n", xmlGetName(nlexpr));
         *doingfine = FALSE;
         break;
      }
      if( count >= nnlexprs )
      {
         SCIPerrorMessage("Too many nonlinear expressions under <nonlinearExpressions> node, expected %d many, but got at least %d.\n", nnlexprs, count + 1);
         *doingfine = FALSE;
         break;
      }

      /* treat empty expression as 0.0 and continue */
      if( xmlFirstChild(nlexprs) == NULL )
         continue;

      /* get constraint index, or -1 for objective */
      attrval = xmlGetAttrval(nlexpr, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing \"idx\" attribute in %d'th <nl> node under <nonlinearExpressions> node.\n", count);
         *doingfine = FALSE;
         break;
      }

      considx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || considx < -1 || considx >= nconss )
      {
         SCIPerrorMessage("Invalid value '%s' in \"idx\" attribute of %d'th <nl> node under <nonlinearExpressions> node.\n", xmlGetAttrval(nlexpr, "idx"), count);
         *doingfine = FALSE;
         break;
      }

      /* turn OSiL expression into SCIP expression and assign indices to variables; store a nonlinear objective at position nconss */
      SCIP_CALL( readExpression(scip, considx == -1 ? &exprs[nconss] : &exprs[considx],
         xmlFirstChild(nlexpr), vars, nvars, doingfine) );
      if( !*doingfine )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** read sos1 and sos2 constraints
 *
 *  sos constraints are expected to be given as a node of \<instanceData\> in the following way:
 *    @code
 *    <specialOrderedSets numberOfSpecialOrderedSets="1">
 *       <sos numberOfVar="2" order="2">
 *           <var idx="1"></var>
 *           <var idx="2"></var>
 *       </sos>
 *    </specialOrderedSets>
 *    @endcode
 * Weights are determined by the order in which the variables are given
 *
 */
static
SCIP_RETCODE readSOScons(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             initialconss,       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss,       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamicrows,        /**< should rows be added and removed dynamically to the LP? */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* soscons;
   const XML_NODE* sosvar;
   const char* attrval;
   int nsoscons;
   int nsosvars;
   int sosorder;
   int type;
   int count;
   int varcount;
   int idx;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   char name[SCIP_MAXSTRLEN];

   /* standard settings for SOS constraints: */
   initial = initialconss;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   dynamic = dynamicconss;
   removable = dynamicrows;

   soscons= xmlFindNodeMaxdepth(datanode, "specialOrderedSets", 0, 1);

   if( soscons== NULL )
      return SCIP_OKAY;

   /* get number of sos constraints */
   attrval = xmlGetAttrval(soscons, "numberOfSOS");
   if( attrval == NULL )
   {
      SCIPerrorMessage("Attribute \"numberOfSOS in <specialOrderedSets> node not found.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   nsoscons = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || nsoscons < 0 )
   {
      SCIPerrorMessage("Invalid value '%s' for \"numberOfSOS\" attribute in <specialOrderedSets>.\n", xmlGetAttrval(soscons, "numberOfSOS"));
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(nsoscons >= 0);

   /* read sos constraints and create corresponding constraint */
   count = 0;
   for( soscons = xmlFirstChild(soscons); soscons != NULL; soscons = xmlNextSibl(soscons), ++count )
   {
      SCIP_CONS* cons;

      /* Make sure we get a sos node and not more then announced*/
      if( strcmp(xmlGetName(soscons), "sos") != 0 )
      {
         SCIPerrorMessage("Expected <sos> node under <specialOrderedSet> node, but got '%s'.\n", xmlGetName(soscons));
         *doingfine = FALSE;
         break;
      }

      if( count >= nsoscons)
      {
         SCIPerrorMessage("Too many sos under <specialOrderedSets> node, expected %d many, but got at least %d.\n", nsoscons, count + 1);
         *doingfine = FALSE;
         break;
      }

      /* get number of variables in this sos constraint */
      attrval = xmlGetAttrval(soscons, "numberOfVar");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Attribute \"numberOfVar in <sos> node not found.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      nsosvars = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || nsosvars < 0 )
      {
         SCIPerrorMessage("Invalid value '%s' for \"numberOfVar\" attribute in <sos>.\n", xmlGetAttrval(soscons, "numberOfVar"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      assert(nsosvars >= 0);

      /* get order of this sos constraint */
      attrval = xmlGetAttrval(soscons, "type");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Attribute \"order\" in <sos> node not found.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      sosorder = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || sosorder < 0 || sosorder > 2 )
      {
         SCIPerrorMessage("Invalid/unsupported value '%s' for \"order\" attribute in <sos>.\n", xmlGetAttrval(soscons, "order"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      assert(sosorder == 1 || sosorder == 2);
      type = sosorder;

      /* set artificial name for sos constraint*/
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "SOS%d_%d", type, count);

      /* Create sos constraint */
      switch( type )
      {
      case 1:
         SCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
            local, dynamic, removable, FALSE) );
         break;
      case 2:
         SCIP_CALL( SCIPcreateConsSOS2(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
            local, dynamic, removable, FALSE) );
         break;
      default:
         SCIPerrorMessage("unknown SOS type: <%d>\n", type); /* should not happen */
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      }

      varcount = 0;
      for( sosvar = xmlFirstChild(soscons); sosvar!= NULL; sosvar = xmlNextSibl(sosvar), ++varcount )
      {
         /* get variable id*/
          attrval = xmlGetAttrval(sosvar, "idx");
          if( attrval == NULL )
          {
             SCIPerrorMessage("Attribute \"idx\" in <var> node below <specialOrderedSets> node not found.\n");
             *doingfine = FALSE;
             return SCIP_OKAY;
          }

          idx = (int)strtol(attrval, (char**)&attrval, 10);
          if( *attrval != '\0' || idx < 0  || idx > nvars - 1 )
          {
             SCIPerrorMessage("Invalid value '%s' for \"idx\" attribute in <var>.\n", xmlGetAttrval(sosvar, "idx"));
             *doingfine = FALSE;
             return SCIP_OKAY;
          }
          assert(idx >= 0);

          /* we now know that we have a variable/weight pair -> add variable*/
          switch( type )
          {
          case 1:
             SCIP_CALL( SCIPaddVarSOS1(scip, cons, vars[idx], (SCIP_Real) (nsosvars - varcount)) );
             break;
          case 2:
             SCIP_CALL( SCIPaddVarSOS2(scip, cons, vars[idx], (SCIP_Real) (nsosvars - varcount)) );
             break;
          /* coverity[dead_error_begin] */
          default:
             SCIPerrorMessage("unknown SOS type: <%d>\n", type); /* should not happen */
             SCIPABORT();
             return SCIP_INVALIDDATA;  /*lint !e527*/
          }
      } /* Close loop over variables in sos constraint */

      /* add the SOS constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

 /*
 * Callback methods of reader
 */


/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyOsil)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeReaderOsil(scip) );

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadOsil)
{  /*lint --e{715}*/
   const XML_NODE* header;
   const XML_NODE* data;
   XML_NODE* start;
   SCIP_VAR** vars;
   const char* name;
   SCIP_RETCODE retcode;
   SCIP_Bool doingfine;
   SCIP_Bool initialconss;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamiccols;
   SCIP_Bool dynamicrows;
   int nconss;
   int nvars;
   int c;
   int i;

   /* linear parts */
   SCIP_VAR*** linvars = NULL;
   SCIP_Real** lincoefs = NULL;
   int* nlinvars = NULL;

   /* quadratic parts */
   SCIP_VAR*** quadvars1 = NULL;
   SCIP_VAR*** quadvars2 = NULL;
   SCIP_Real** quadcoefs = NULL;
   int* nquadterms = NULL;

   /* nonlinear parts */
   SCIP_EXPR** nlexprs = NULL;

   assert(scip != NULL);
   assert(reader != NULL);
   assert(result != NULL);
   assert(filename != NULL);

   *result = SCIP_DIDNOTRUN;
   retcode = SCIP_READERROR;
   doingfine = TRUE;
   vars = NULL;
   nvars = 0;
   nconss = -1;

   /* read OSiL xml file */
   start = xmlProcess(filename);

   if( start == NULL )
   {
      SCIPerrorMessage("Some error occurred when parsing the OSiL XML file '%s'.\n", filename);
      goto CLEANUP;
   }

   SCIPdebug( xmlShowNode(start) );

   /* parse header to get problem name */
   name = filename;
   header = xmlFindNodeMaxdepth(start, "instanceHeader", 0, 2);
   if( header != NULL )
   {
      const XML_NODE* namenode;

      namenode = xmlFindNodeMaxdepth(header, "name", 0, 2);

      if( namenode != NULL && xmlFirstChild(namenode) != NULL )
         name = xmlGetData(xmlFirstChild(namenode));
      else
      {
         namenode = xmlFindNodeMaxdepth(header, "description", 0, 2);

         if( namenode != NULL && xmlFirstChild(namenode) != NULL )
            name = xmlGetData(xmlFirstChild(namenode));
      }
   }

   /* create SCIP problem */
   SCIP_CALL( SCIPcreateProb(scip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* process instance data */
   data = xmlFindNodeMaxdepth(start, "instanceData", 0, 2);
   if( data == NULL )
   {
      SCIPerrorMessage("Node <instanceData> not found.\n");
      goto CLEANUP;
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/initialconss", &initialconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamiccols", &dynamiccols) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicrows", &dynamicrows) );

   /* read variables */
   SCIP_CALL_TERMINATE( retcode, readVariables(scip, data, &vars, &nvars, initialconss, dynamicconss, dynamiccols, dynamicrows, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;
   assert(vars != NULL || nvars == 0);

   /* read objective sense, coefficient, and constant */
   SCIP_CALL_TERMINATE( retcode, readObjective(scip, data, vars, nvars, dynamiccols, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* read total number of constraints */
   SCIP_CALL_TERMINATE( retcode, readNConstraints(scip, data, &nconss, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* allocate memory to store constraint information */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &linvars, nconss + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &lincoefs, nconss + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nlinvars, nconss + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &quadvars1, nconss + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &quadvars2, nconss + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &quadcoefs, nconss + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nquadterms, nconss + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nlexprs, nconss + 1) );

   /* read linear coefficients matrix */
   SCIP_CALL_TERMINATE( retcode, readLinearCoefs(scip, data, vars, nvars, nconss, linvars, lincoefs, nlinvars, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* read quadratic coefficients */
   SCIP_CALL_TERMINATE( retcode, readQuadraticCoefs(scip, data, vars, nvars, nconss, quadvars1, quadvars2, quadcoefs,
      nquadterms, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* read nonlinear expressions */
   SCIP_CALL_TERMINATE( retcode, readNonlinearExprs(scip, data, vars, nvars, nconss, nlexprs, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* read constraint data; generate constraints */
   SCIP_CALL_TERMINATE( retcode, readConstraints(scip, data, nconss, linvars, lincoefs, nlinvars,
      quadvars1, quadvars2, quadcoefs, nquadterms, nlexprs, initialconss, dynamicconss, dynamicrows, &doingfine),
      CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* add nonlinear objective constraint */
   if( nlinvars[nconss] > 0 || nquadterms[nconss] > 0 || nlexprs[nconss] != NULL )
   {
      SCIP_CALL( createConstraint(scip, linvars[nconss], lincoefs[nconss], nlinvars[nconss],
         quadvars1[nconss], quadvars2[nconss], quadcoefs[nconss], nquadterms[nconss], nlexprs[nconss],
         SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
         SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0,
         "objcons", TRUE, TRUE, FALSE, FALSE) );
   }

   /* read sos2 constraints and add to problem */
   SCIP_CALL_TERMINATE( retcode, readSOScons(scip, data, vars, nvars, initialconss, dynamicconss, dynamicrows, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   *result = SCIP_SUCCESS;
   retcode = SCIP_OKAY;

 CLEANUP:
   /* free xml data */
   if( start != NULL )
      xmlFreeNode(start);

   /* free memory for constraint information (position nconss belongs to the nonlinear objective function) */
   for( c = nconss; c >= 0; --c )
   {
      /* free nonlinear parts */
      if( nlexprs != NULL && nlexprs[c] != NULL )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &nlexprs[c]) );
      }

      /* free quadratic parts */
      SCIPfreeBufferArrayNull(scip, &quadcoefs[c]);
      SCIPfreeBufferArrayNull(scip, &quadvars1[c]);
      SCIPfreeBufferArrayNull(scip, &quadvars2[c]);

      /* free linear parts */
      SCIPfreeBufferArrayNull(scip, &lincoefs[c]);
      SCIPfreeBufferArrayNull(scip, &linvars[c]);
   }
   SCIPfreeBufferArrayNull(scip, &nlexprs);
   SCIPfreeBufferArrayNull(scip, &nquadterms);
   SCIPfreeBufferArrayNull(scip, &quadcoefs);
   SCIPfreeBufferArrayNull(scip, &quadvars2);
   SCIPfreeBufferArrayNull(scip, &quadvars1);
   SCIPfreeBufferArrayNull(scip, &nlinvars);
   SCIPfreeBufferArrayNull(scip, &lincoefs);
   SCIPfreeBufferArrayNull(scip, &linvars);

   /* free variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );  /*lint !e613*/
   }
   SCIPfreeBufferArrayNull(scip, &vars);

   /* return read error retcode if something went wrong */
   if( !doingfine )
      return SCIP_READERROR;

   if( retcode == SCIP_PLUGINNOTFOUND )
      retcode = SCIP_READERROR;

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the osil file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderOsil(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include osil reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyOsil) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadOsil) );

   return SCIP_OKAY;
}
