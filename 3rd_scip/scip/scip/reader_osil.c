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

/**@file   reader_osil.c
 * @brief  OS instance language (OSiL) format file reader
 * @author Stefan Vigerske
 * @author Ingmar Vierhaus
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI and M_E on Windows */  /*lint !750 */

#include <assert.h>
#include <string.h>
#include <math.h>

#include "scip/reader_osil.h"
#include "scip/scip.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_linear.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "xml/xml.h"


#define READER_NAME             "osilreader"
#define READER_DESC             "file reader for OS instance language (OSiL) format"
#define READER_EXTENSION        "osil"

/*
 * Data structures
 */

/** type of constraint */
typedef enum
{
   LINEAR,               /**< linear constraint */
   QUADRATIC,            /**< quadratic constraint */
   NONLINEAR             /**< general nonlinear constraint */
} CONSTYPE;


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
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
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
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
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

/** setup constraint sides as linear constraints
 *
 * constraints are not added to the problem yet
 */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_CONS***          conss,              /**< buffer to store array of (linear) constraints */
   CONSTYPE**            constypes,          /**< buffer to store type of constraints (will be all LINEAR) */
   int*                  nconss,             /**< buffer to store number of constraints */
   SCIP_Bool             initialconss,       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss,       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamicrows,        /**< should rows be added and removed dynamically to the LP? */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
   )
{
   const XML_NODE* constraints;
   const XML_NODE* consnode;
   const char* attrval;
   int consssize;
   char name[20];

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(conss != NULL);
   assert(constypes != NULL);
   assert(nconss != NULL);
   assert(doingfine != NULL);

   *conss = NULL;
   *constypes = NULL;
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

   consssize = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || consssize < 0 )
   {
      SCIPerrorMessage("Invalid value '%s' for \"numberOfConstraints\" attribute.\n", xmlGetAttrval(constraints, "numberOfConstraints"));
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(consssize >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, conss, consssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, constypes, consssize) );

   /* read constraint names, lhs, rhs, constant */
   for( consnode = xmlFirstChild(constraints); consnode != NULL; consnode = xmlNextSibl(consnode) )
   {
      const char* consname;
      SCIP_Real conslhs;
      SCIP_Real consrhs;

      if( consssize == *nconss )
      {
         SCIPerrorMessage("Expected %d constraints, but got at least %d many.\n", consssize, *nconss+1);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* find constraint name */
      consname = xmlGetAttrval(consnode, "name");
      if( consname == NULL )
      {
         (void) SCIPsnprintf(name, 20, "cons%d", *nconss);
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

      /* create SCIP linear constraint */
      SCIP_CALL( SCIPcreateConsLinear(scip, &(*conss)[*nconss], consname, 0, NULL, NULL, conslhs, consrhs,
            initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows, FALSE) );
      assert((*conss)[*nconss] != NULL);

      (*constypes)[*nconss] = LINEAR;

      ++*nconss;
   }

   if( *nconss < consssize )
   {
      SCIPerrorMessage("Got %d constraints, but expected %d many.\n", *nconss, consssize);
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
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
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
   SCIP_CONS**           conss,              /**< constraints in order of OSiL indices */
   CONSTYPE*             constypes,          /**< type of constraints (assumed to be LINEAR) */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
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
   assert(conss != NULL || nconss == 0);
   assert(constypes != NULL || nconss == 0);
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
         /* these asserts were checked above */
         assert(start[row] >= 0);
         assert(start[row+1] >= 0);
         assert(start[row] <= nnz);
         assert(start[row+1] <= nnz);
         for( pos = start[row]; pos < start[row+1]; ++pos )
         {
            /* these asserts were checked above */
            assert(pos >= 0);
            assert(pos < nnz);
            assert(idx[pos] >= 0);
            assert(idx[pos] < nvars);

            assert(constypes[row] == LINEAR);  /*lint !e613*/

            SCIP_CALL( SCIPaddCoefLinear(scip, conss[row], vars[idx[pos]], val[pos]) );  /*lint !e613*/
         }
      }
   }
   else
   {
      int col;
      int pos;
      for( col = 0; col < nvars; ++col )
      {
         /* these asserts were checked above */
         assert(start[col] >= 0);
         assert(start[col+1] >= 0);
         assert(start[col] <= nnz);
         assert(start[col+1] <= nnz);
         for( pos = start[col]; pos < start[col+1]; ++pos )
         {
            /* these asserts were checked above */
            assert(pos >= 0);
            assert(pos < nnz);
            assert(idx[pos] >= 0);
            assert(idx[pos] < nconss);

            assert(constypes[idx[pos]] == LINEAR);  /*lint !e613*/

            SCIP_CALL( SCIPaddCoefLinear(scip, conss[idx[pos]], vars[col], val[pos]) );  /*lint !e613*/
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
   SCIP_CONS**           conss,              /**< constraints in order of OSiL indices */
   CONSTYPE*             constypes,          /**< type of constraints (assumed to be LINEAR) */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           objcons,            /**< buffer to store constraint for nonlinear part of objective function, or to add to if already existing */
   CONSTYPE*             objconstype,        /**< buffer to store type of objective constraint, if created (should be QUADRATIC) */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
   )
{
   const XML_NODE* quadcoef;
   const XML_NODE* qterm;
   const char* attrval;
   SCIP_CONS* cons;
   int nqterms;
   int count;
   int considx;
   int varidx1;
   int varidx2;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(conss != NULL || nconss == 0);
   assert(constypes != NULL || nconss == 0);
   assert(objcons != NULL);
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

   count = 0;
   for( qterm = xmlFirstChild(quadcoef); qterm != NULL; qterm = xmlNextSibl(qterm), ++count )
   {
      /* check for qterm node */
      if( strcmp(xmlGetName(qterm), "qTerm") != 0 )
      {
         SCIPerrorMessage("Expected <qTerm> node under <quadraticCoefficients> node, but got <%s>\n", xmlGetName(qterm));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      if( count >= nqterms )
      {
         SCIPerrorMessage("Too many quadratic terms under <quadraticCoefficients> node, expected %d many, but got at least %d.\n", nqterms, count + 1);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* get constraint index, or -1 for objective */
      attrval = xmlGetAttrval(qterm, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing \"idx\" attribute in %d'th <qTerm> node under <quadraticCoefficients> node.\n", count);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      considx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || considx < -1 || considx >= nconss )
      {
         SCIPerrorMessage("Invalid value '%s' in \"idx\" attribute of %d'th <qTerm> node under <quadraticCoefficients> node.\n", xmlGetAttrval(qterm, "idx"), count);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* get index of first variable */
      attrval = xmlGetAttrval(qterm, "idxOne");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing \"idxOne\" attribute in %d'th <qTerm> node under <quadraticCoefficients> node.\n", count);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      varidx1 = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || varidx1 < 0 || varidx1 >= nvars )
      {
         SCIPerrorMessage("Invalid value '%s' in \"idxOne\" attribute of %d'th <qTerm> node under <quadraticCoefficients> node.\n", xmlGetAttrval(qterm, "idxOne"), count);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* get index of second variable */
      attrval = xmlGetAttrval(qterm, "idxTwo");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing \"idxTwo\" attribute in %d'th <qTerm> node under <quadraticCoefficients> node.\n", count);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      varidx2 = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || varidx2 < 0 || varidx2 >= nvars )
      {
         SCIPerrorMessage("Invalid value '%s' in \"idxTwo\" attribute of %d'th <qTerm> node under <quadraticCoefficients> node.\n", xmlGetAttrval(qterm, "idxTwo"), count);
         *doingfine = FALSE;
         return SCIP_OKAY;
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
            return SCIP_OKAY;
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

      if( considx == -1 )
      {
         if( *objcons == NULL )
         {
            /* create constraint to hold quadratic part of objective; note that
             * reading/{initialconss,dynamicconss,dynamicrows,dynamiccols} apply only to model constraints and
             * variables, not to an auxiliary objective constraint (otherwise it can happen that an auxiliary objective
             * variable is loose with infinite best bound, triggering the problem that an LP that is unbounded because
             * of loose variables with infinite best bound cannot be solved)
             */

            SCIP_VAR* objvar;
            SCIP_Real minusone;

            SCIP_CALL( SCIPcreateVar(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, objvar) );

            minusone = -1.0;
            SCIP_CALL( SCIPcreateConsQuadratic(scip, objcons, "objcons", 1, &objvar, &minusone, 0, NULL, NULL, NULL,
                  SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
                  SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
            *objconstype = QUADRATIC;

            SCIP_CALL( SCIPreleaseVar(scip, &objvar) );
         }
         cons = *objcons;
         assert(*objconstype == QUADRATIC);
      }
      else if( constypes[considx] == LINEAR )  /*lint !e613*/
      {
         /* replace linear constraint by quadratic constraint */
         cons = conss[considx];  /*lint !e613*/

         SCIP_CALL( SCIPcreateConsQuadratic(scip, &cons, SCIPconsGetName(cons),
            SCIPgetNVarsLinear(scip, cons), SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
            0, NULL, NULL, NULL,
            SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );

         SCIP_CALL( SCIPreleaseCons(scip, &conss[considx]) );  /*lint !e613*/

         conss[considx] = cons;  /*lint !e613*/
         constypes[considx] = QUADRATIC;  /*lint !e613*/
      }
      else
      {
         cons = conss[considx];  /*lint !e613*/
         assert(constypes[considx] == QUADRATIC);  /*lint !e613*/
      }

      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, vars[varidx1], vars[varidx2], coef) );  /*lint !e613*/
   }

   if( count != nqterms )
   {
      SCIPerrorMessage("Got only %d quadratic terms under <quadraticCoefficients> node, but expected %d many.\n", count, nqterms);
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** transforms OSnL expression tree into SCIP expression */
static
SCIP_RETCODE readExpression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to created expression */
   const XML_NODE*       node,               /**< root node of expression to be read */
   int*                  exprvaridx,         /**< array with index of problem variables in expression graph */
   int*                  nexprvars,          /**< number of variables in currently processed expression so far */
   int                   nvars,              /**< total number of variables in problem (and length of exprvaridx array) */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const char* exprname;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(node != NULL);
   assert(exprvaridx != NULL || nvars == 0);
   assert(nexprvars != NULL);
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

      /* assign index to variable, if we see it the first time */
      if( exprvaridx[idx] == -1 )  /*lint !e613*/
      {
         exprvaridx[idx] = *nexprvars;  /*lint !e613*/
         ++*nexprvars;
      }

      /* create VARIDX expression, put into LINEAR expression if we have coefficient != 1 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_VARIDX, exprvaridx[idx]) );  /*lint !e613*/
      if( coef != 1.0 )
      {
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), expr, 1, expr, &coef, 0.0) );
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

      /* create CONST expression */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, val) );

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "PI") == 0 )
   {
      /* create CONST expression with PI value*/
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, M_PI) );

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "E") == 0 )
   {
      /* create CONST expression with E value*/
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, M_E) );

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
      strcmp(exprname, "log10") == 0
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
      SCIP_CALL( readExpression(scip, &arg, xmlFirstChild(node), exprvaridx, nexprvars, nvars, doingfine) );
      if( !*doingfine )
         return SCIP_OKAY;

      /* create SCIP expression according to expression name */
      if( strcmp(exprname, "negate") == 0 )
      {
         SCIP_Real minusone;

         minusone = -1.0;
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), expr, 1, &arg, &minusone, 0.0) );
      }
      else if( strcmp(exprname, "abs") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_ABS, arg) );
      }
      else if( strcmp(exprname, "squareRoot") == 0 || strcmp(exprname, "sqrt") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_SQRT, arg) );
      }
      else if( strcmp(exprname, "square") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_SQUARE, arg) );
      }
      else if( strcmp(exprname, "exp") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_EXP, arg) );
      }
      else if( strcmp(exprname, "ln") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_LOG, arg) );
      }
      else /* if( strcmp(exprname, "log10") == 0 ) */
      {
         /* log10(expr) = ln(expr)*1/ln(10) */
         SCIP_EXPR* tmp;

         assert(strcmp(exprname, "log10") == 0);

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, 1.0/log(10.0)) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg, SCIP_EXPR_LOG, arg) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MUL, arg, tmp) );
      }

      return SCIP_OKAY;
   }

   /* two argument operands */
   if( strcmp(exprname, "plus") == 0 ||
      strcmp(exprname, "minus") == 0 ||
      strcmp(exprname, "times") == 0 ||
      strcmp(exprname, "divide") == 0 ||
      strcmp(exprname, "power") == 0 ||
      strcmp(exprname, "log") == 0
     )
   {
      SCIP_EXPR* arg1;
      SCIP_EXPR* arg2;

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
      SCIP_CALL( readExpression(scip, &arg1, xmlFirstChild(node), exprvaridx, nexprvars, nvars, doingfine) );
      if( !*doingfine )
         return SCIP_OKAY;

      /* read second child expression */
      SCIP_CALL( readExpression(scip, &arg2, xmlNextSibl(xmlFirstChild(node)), exprvaridx, nexprvars, nvars, doingfine) );
      if( !*doingfine )
      {
         SCIPexprFreeDeep(SCIPblkmem(scip), &arg1);
         return SCIP_OKAY;
      }

      if( strcmp(exprname, "plus") == 0 )
      {
         SCIP_CALL( SCIPexprAdd(SCIPblkmem(scip), expr, 1.0, arg1, 1.0, arg2, 0.0) );
         /* SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_PLUS, arg1, arg2) ); */
      }
      else if( strcmp(exprname, "minus") == 0 )
      {
         SCIP_CALL( SCIPexprAdd(SCIPblkmem(scip), expr, 1.0, arg1, -1.0, arg2, 0.0) );
         /* SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MINUS, arg1, arg2) ); */
      }
      else if( strcmp(exprname, "times") == 0 )
      {
         if( SCIPexprGetOperator(arg1) == SCIP_EXPR_CONST )
         {
            SCIP_CALL( SCIPexprMulConstant(SCIPblkmem(scip), expr, arg2, SCIPexprGetOpReal(arg1)) );
            SCIPexprFreeDeep(SCIPblkmem(scip), &arg1);
         }
         else if( SCIPexprGetOperator(arg2) == SCIP_EXPR_CONST )
         {
            SCIP_CALL( SCIPexprMulConstant(SCIPblkmem(scip), expr, arg1, SCIPexprGetOpReal(arg2)) );
            SCIPexprFreeDeep(SCIPblkmem(scip), &arg2);
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MUL, arg1, arg2) );
         }
      }
      else if( strcmp(exprname, "divide") == 0 )
      {
         if( SCIPexprGetOperator(arg2) == SCIP_EXPR_CONST )
         {
            assert(SCIPexprGetOpReal(arg2) != 0.0);
            SCIP_CALL( SCIPexprMulConstant(SCIPblkmem(scip), expr, arg1, 1.0/SCIPexprGetOpReal(arg2)) );
            SCIPexprFreeDeep(SCIPblkmem(scip), &arg2);
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_DIV, arg1, arg2) );
         }
      }
      else if( strcmp(exprname, "power") == 0 )
      {
         if( SCIPexprGetOperator(arg2) == SCIP_EXPR_CONST )
         {
            /* expr^number is intpower or realpower */
            if( SCIPisIntegral(scip, SCIPexprGetOpReal(arg2)) )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_INTPOWER, arg1, (int)SCIPround(scip, SCIPexprGetOpReal(arg2))) );
            }
            else
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_REALPOWER, arg1, SCIPexprGetOpReal(arg2)) );
            }
            SCIPexprFreeDeep(SCIPblkmem(scip), &arg2);
         }
         else if( SCIPexprGetOperator(arg1) == SCIP_EXPR_CONST )
         {
            /* number^arg2 is exp(arg2 * ln(number)) */
            if( SCIPexprGetOpReal(arg1) < 0.0 )
            {
               SCIPerrorMessage("Negative base in <power> node with nonconstant exponent not allowed in nonlinear expression.\n");
               SCIPexprFreeDeep(SCIPblkmem(scip), &arg1);
               SCIPexprFreeDeep(SCIPblkmem(scip), &arg2);
               *doingfine = FALSE;
               return SCIP_OKAY;
            }
            else
            {
               SCIP_EXPR* tmp;

               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, log(SCIPexprGetOpReal(arg1))) );
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_MUL, tmp, arg2) );
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_EXP, tmp) );
               SCIPexprFreeDeep(SCIPblkmem(scip), &arg1);
            }
         }
         else
         {
            /* arg1^arg2 is exp(arg2 * ln(arg1)) */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg1, SCIP_EXPR_LOG, arg1) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg2, SCIP_EXPR_MUL, arg1, arg2) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_EXP, arg2) );
         }
      }
      else if( strcmp(exprname, "log") == 0 )
      {
         /* logarithm of arg2 w.r.t. base arg1 = ln(arg2) / ln(arg1) */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg1, SCIP_EXPR_LOG, arg1) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg2, SCIP_EXPR_LOG, arg2) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_DIV, arg2, arg1) );
      }
      else if( strcmp(exprname, "min") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MIN, arg1, arg2) );
      }
      else /* if( strcmp(exprname, "max") == 0 ) */
      {
         assert(strcmp(exprname, "max") == 0);

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MAX, arg1, arg2) );
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

      /* a sum or product w.r.t. 0 arguments is constant */
      if( xmlFirstChild(node) == NULL )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, (strcmp(exprname, "sum") == 0) ? 0.0 : 1.0) );

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

         SCIP_CALL( readExpression(scip, &args[nargs], argnode, exprvaridx, nexprvars, nvars, doingfine) );
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
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, (strcmp(exprname, "sum") == 0) ? 0.0 : 1.0) );
               break;
            }
            case 1:
            {
               *expr = args[0];
               break;
            }
            case 2:
            {
               if( strcmp(exprname, "sum") == 0 )
               {
                  SCIP_CALL( SCIPexprAdd(SCIPblkmem(scip), expr, 1.0, args[0], 1.0, args[1], 0.0) );
               }
               else
               {
                  SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MUL, args[0], args[1]) );
               }
               break;
            }
            default:
            {
               /* create sum or product expression */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, (strcmp(exprname, "sum") == 0) ? SCIP_EXPR_SUM : SCIP_EXPR_PRODUCT, nargs, args) );
               break;
            }
         }
      }
      else
      {
         /* cleanup if parsing error */
         for( ; nargs > 0; --nargs )
            SCIPexprFreeDeep(SCIPblkmem(scip), &args[nargs-1]);
      }

      SCIPfreeBufferArray(scip, &args);

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "min") == 0 || strcmp(exprname, "max") == 0 )
   {
      const XML_NODE* argnode;
      SCIP_EXPROP exprop;
      SCIP_EXPR* arg2;

      /* check that we have children */
      if( xmlFirstChild(node) == NULL )
      {
         SCIPerrorMessage("Expected at least one child in <%s> node of nonlinear expression.\n", exprname);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* read expression corresponding to first child and store in expr */
      argnode = xmlFirstChild(node);
      SCIP_CALL( readExpression(scip, expr, argnode, exprvaridx, nexprvars, nvars, doingfine) );
      if( !*doingfine )
      {
         assert(*expr == NULL);
         return SCIP_OKAY;
      }
      arg2 = NULL;

      exprop = (strcmp(exprname, "min") == 0) ? SCIP_EXPR_MIN : SCIP_EXPR_MAX;

      /* read expressions corresponding to other children in arg and store exprop(expr, arg) in expr */
      for( argnode = xmlNextSibl(argnode); argnode != NULL; argnode = xmlNextSibl(argnode) )
      {
         assert(arg2 == NULL);
         SCIP_CALL( readExpression(scip, &arg2, argnode, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
         {
            assert(arg2 == NULL);
            break;
         }

         assert(*expr != NULL);
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, exprop, *expr, arg2) );
         arg2 = NULL;
      }

      if( !*doingfine )
      {
         /* cleanup if failure */
         SCIPexprFreeDeep(SCIPblkmem(scip), expr);
      }
      assert(arg2 == NULL);

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "quadratic") == 0 )
   {
      const char* attrval;
      const XML_NODE* qterm;
      SCIP_QUADELEM* quadelems;
      int nquadelems;
      int quadelemssize;
      int* quadvarsidxs;
      int nquadvars;
      int i;

      quadelemssize = 5;
      SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, quadelemssize) );
      nquadelems = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &quadvarsidxs, nvars) );
      for( i = 0; i < nvars; ++i )
         quadvarsidxs[i] = -1;
      nquadvars = 0;

      /* read quadratic terms */
      for( qterm = xmlFirstChild(node); qterm != NULL; qterm = xmlNextSibl(qterm), ++nquadelems )
      {
         /* check for qpTerm node */
         if( strcmp(xmlGetName(qterm), "qpTerm") != 0 )
         {
            SCIPerrorMessage("Unexpected <%s> node under <quadratic> node in nonlinear expression, expected <qpTerm>.\n", xmlGetName(qterm));
            *doingfine = FALSE;
            break;
         }

         if( nquadelems >= quadelemssize )
         {
            quadelemssize = SCIPcalcMemGrowSize(scip, nquadelems + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &quadelems, quadelemssize) );
         }
         assert(quadelemssize > nquadelems);

         /* get index of first variable */
         attrval = xmlGetAttrval(qterm, "idxOne");
         if( attrval == NULL )
         {
            SCIPerrorMessage("Missing \"idxOne\" attribute in %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", nquadelems);
            *doingfine = FALSE;
            break;
         }

         quadelems[nquadelems].idx1 = (int)strtol(attrval, (char**)&attrval, 10);
         if( *attrval != '\0' || quadelems[nquadelems].idx1 < 0 || quadelems[nquadelems].idx1 >= nvars )
         {
            SCIPerrorMessage("Invalid value '%s' for \"idxOne\" attribute of %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", xmlGetAttrval(qterm, "idxOne"), nquadelems);
            *doingfine = FALSE;
            break;
         }

         /* get index of second variable */
         attrval = xmlGetAttrval(qterm, "idxTwo");
         if( attrval == NULL )
         {
            SCIPerrorMessage("Missing \"idxTwo\" attribute in %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", nquadelems);
            *doingfine = FALSE;
            break;
         }

         quadelems[nquadelems].idx2 = (int)strtol(attrval, (char**)&attrval, 10);
         if( *attrval != '\0' || quadelems[nquadelems].idx2 < 0 || quadelems[nquadelems].idx2 >= nvars )
         {
            SCIPerrorMessage("Invalid value '%s' for \"idxTwo\" attribute of %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", xmlGetAttrval(qterm, "idxTwo"), nquadelems);
            *doingfine = FALSE;
            break;
         }

         /* get coefficient */
         attrval = xmlGetAttrval(qterm, "coef");
         if( attrval != NULL )
         {
            quadelems[nquadelems].coef = strtod(attrval, (char**)&attrval);
            if( *attrval != '\0' || (quadelems[nquadelems].coef != quadelems[nquadelems].coef) )  /*lint !e777*/
            {
               SCIPerrorMessage("Invalid value '%s' for \"coef\" attribute of %d'th <qpTerm> node under <quadratic> node in nonlinear expression.\n", xmlGetAttrval(qterm, "coef"), nquadelems);
               *doingfine = FALSE;
               break;
            }
         }
         else
         {
            quadelems[nquadelems].coef = 1.0;
         }

         /* get index for first variable in quadratic element */
         if( quadvarsidxs[quadelems[nquadelems].idx1] < 0 )
         {
            quadvarsidxs[quadelems[nquadelems].idx1] = nquadvars;
            quadelems[nquadelems].idx1 = nquadvars;

            ++nquadvars;
         }
         else
         {
            quadelems[nquadelems].idx1 = quadvarsidxs[quadelems[nquadelems].idx1];
         }

         /* get index for second variable in quadratic element */
         if( quadvarsidxs[quadelems[nquadelems].idx2] < 0 )
         {
            quadvarsidxs[quadelems[nquadelems].idx2] = nquadvars;
            quadelems[nquadelems].idx2 = nquadvars;

            ++nquadvars;
         }
         else
         {
            quadelems[nquadelems].idx2 = quadvarsidxs[quadelems[nquadelems].idx2];
         }

         /* swap indices if in wrong order */
         if( quadelems[nquadelems].idx1 > quadelems[nquadelems].idx2 )
         {
            int tmp;

            tmp = quadelems[nquadelems].idx1;
            quadelems[nquadelems].idx1 = quadelems[nquadelems].idx2;
            quadelems[nquadelems].idx2 = tmp;
         }
      }

      if( *doingfine )
      {
         SCIP_EXPR** children;

         /* setup array with children expressions corresponding to variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &children, nquadvars) );
         for( i = 0; i < nvars; ++i )
         {
            if( quadvarsidxs[i] == -1 )
               continue;

            /* assign new index to variable, if we see it the first time in this exprtree */
            if( exprvaridx[i] == -1 )  /*lint !e613*/
            {
               exprvaridx[i] = *nexprvars;  /*lint !e613*/
               ++*nexprvars;
            }

            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[quadvarsidxs[i]], SCIP_EXPR_VARIDX, exprvaridx[i]) );  /*lint !e613*/
         }

         /* create quadratic expression */
         SCIP_CALL( SCIPexprCreateQuadratic(SCIPblkmem(scip), expr, nquadvars, children, 0.0, NULL, nquadelems, quadelems) );

         SCIPfreeBufferArray(scip, &children);
      }

      SCIPfreeBufferArray(scip, &quadelems);
      SCIPfreeBufferArray(scip, &quadvarsidxs);
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
   SCIP_CONS**           conss,              /**< constraints in order of OSiL indices */
   CONSTYPE*             constypes,          /**< type of constraints (assumed to be LINEAR) */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           objcons,            /**< buffer to store constraint for nonlinear part of objective function, or to add to if already existing */
   CONSTYPE*             objconstype,        /**< buffer to store type of objective constraint, if created (should be QUADRATIC) */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const XML_NODE* nlexprs;
   const XML_NODE* nlexpr;
   const char* attrval;
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* expr;
   SCIP_VAR** exprvars;
   int* exprvaridx;
   SCIP_RETCODE retcode;
   int nexprvars;
   int nnlexprs;
   int count;
   int considx;
   int i;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(conss != NULL || nconss == 0);
   assert(constypes != NULL || nconss == 0);
   assert(objcons != NULL);
   assert(doingfine != NULL);

   retcode = SCIP_OKAY;

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

   /* buffer array to store index of variable in expression graph, or -1 if not present */
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvars, nvars) );

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

      expr = NULL;
      nexprvars = 0;
      for( i = 0; i < nvars; ++i )
         exprvaridx[i] = -1;

      /* turn OSiL expression into SCIP expression and assign indices to variables */
      SCIP_CALL( readExpression(scip, &expr, xmlFirstChild(nlexpr), exprvaridx, &nexprvars, nvars, doingfine) );
      if( !*doingfine )
      {
         assert(expr == NULL);
         break;
      }

      /* assemble array with SCIP_VAR*'s */
      for( i = 0; i < nvars; ++i )
      {
         assert(exprvaridx[i] < nexprvars );

         if( exprvaridx[i] >= 0 )
            exprvars[exprvaridx[i]] = vars[i];  /*lint !e613*/
      }

      /* create expression tree */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, nexprvars, 0, NULL) );
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, nexprvars, exprvars) );

      /* add expression tree to objective or constraint */
      if( considx == -1 && *objcons == NULL )
      {
         /* create constraint to hold nonlinear part of objective; note that
          * reading/{initialconss,dynamicconss,dynamicrows,dynamiccols} apply only to model constraints and variables,
          * not to an auxiliary objective constraint (otherwise it can happen that an auxiliary objective variable is
          * loose with infinite best bound, triggering the problem that an LP that is unbounded because of loose
          * variables with infinite best bound cannot be solved)
          */

         SCIP_VAR* objvar;
         SCIP_Real minusone;
         SCIP_Real one;

         SCIP_CALL( SCIPcreateVar(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, objvar) );

         minusone = -1.0;
         one = 1.0;
         SCIP_CALL_TERMINATE( retcode, SCIPcreateConsNonlinear(scip, objcons, "objcons", 1, &objvar, &minusone, 1, &exprtree, &one,
               SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
               SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), TERMINATE );
         *objconstype = NONLINEAR;

         SCIP_CALL( SCIPreleaseVar(scip, &objvar) );
      }
      else
      {
         SCIP_CONS** cons;
         SCIP_CONS* oldcons;
         CONSTYPE* constype;

         if( considx == -1 )
         {
            cons = objcons;
            constype = objconstype;
         }
         else
         {
            cons = &conss[considx];  /*lint !e613*/
            constype = &constypes[considx];  /*lint !e613*/
         }
         oldcons = *cons;

         /* replace cons by nonlinear constraint or add to already existing nonlinear constraint */
         switch( *constype )
         {
         case LINEAR:
         {
            SCIP_Real one;

            one = 1.0;
            SCIP_CALL_TERMINATE( retcode, SCIPcreateConsNonlinear(scip, cons, SCIPconsGetName(*cons),
                  SCIPgetNVarsLinear(scip, *cons), SCIPgetVarsLinear(scip, *cons), SCIPgetValsLinear(scip, *cons),
                  1, &exprtree, &one,
                  SCIPgetLhsLinear(scip, *cons), SCIPgetRhsLinear(scip, *cons),
                  SCIPconsIsInitial(*cons), SCIPconsIsSeparated(*cons), SCIPconsIsEnforced(*cons),
                  SCIPconsIsChecked(*cons), SCIPconsIsPropagated(*cons), SCIPconsIsLocal(*cons),
                  SCIPconsIsModifiable(*cons), SCIPconsIsDynamic(*cons), SCIPconsIsRemovable(*cons),
                  SCIPconsIsStickingAtNode(*cons)), TERMINATE );

            SCIP_CALL( SCIPreleaseCons(scip, &oldcons) );

            break;
         }

         case QUADRATIC:
         {
            SCIP_EXPRTREE* exprtrees[2];
            SCIP_Real exprcoefs[2];

            SCIP_EXPR* quadexpr;
            SCIP_QUADELEM* quadelems;
            SCIP_Real* lincoefs;
            SCIP_EXPR** children;
            SCIP_QUADVARTERM* quadvarterms;
            SCIP_BILINTERM* bilinterms;
            int nquadelems;
            int nquadvars;
            int nbilin;
            int j;

            exprtrees[0] = exprtree;
            exprcoefs[0] = 1.0;

            /* turn quadratic part into expression tree */
            SCIP_CALL( SCIPsortQuadVarTermsQuadratic(scip, *cons) );

            quadvarterms = SCIPgetQuadVarTermsQuadratic(scip, *cons);
            nquadvars = SCIPgetNQuadVarTermsQuadratic(scip, *cons);
            bilinterms = SCIPgetBilinTermsQuadratic(scip, *cons);
            nbilin = SCIPgetNBilinTermsQuadratic(scip, *cons);

            SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nquadvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &children, nquadvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nbilin + nquadvars) );
            nquadelems = 0;

            for( i = 0; i < nquadvars; ++i )
            {
               lincoefs[i] = quadvarterms[i].lincoef;
               exprvars[i] = quadvarterms[i].var;

               if( quadvarterms[i].sqrcoef != 0.0 )
               {
                  quadelems[nquadelems].idx1 = i;
                  quadelems[nquadelems].idx2 = i;
                  quadelems[nquadelems].coef = quadvarterms[i].sqrcoef;
                  ++nquadelems;
               }

               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[i], SCIP_EXPR_VARIDX, i) );

               for( j = 0; j < quadvarterms[i].nadjbilin; ++j )
               {
                  if( bilinterms[quadvarterms[i].adjbilin[j]].var1 == quadvarterms[i].var )
                  {
                     int otheridx;

                     assert(bilinterms[quadvarterms[i].adjbilin[j]].var2 != quadvarterms[i].var);

                     SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, *cons, bilinterms[quadvarterms[i].adjbilin[j]].var2, &otheridx) );
                     assert(otheridx >= 0);
                     assert(otheridx < nquadvars);

                     quadelems[nquadelems].idx1 = MIN(i, otheridx);
                     quadelems[nquadelems].idx2 = MAX(i, otheridx);
                     quadelems[nquadelems].coef = bilinterms[quadvarterms[i].adjbilin[j]].coef;
                     ++nquadelems;
                  }
               }
            }

            SCIP_CALL( SCIPexprCreateQuadratic(SCIPblkmem(scip), &quadexpr, nquadvars, children, 0.0, lincoefs, nquadelems, quadelems) );
            SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtrees[1], quadexpr, nquadvars, 0, NULL) );
            SCIP_CALL( SCIPexprtreeSetVars(exprtrees[1], nquadvars, exprvars) );
            exprcoefs[1] = 1.0;

            SCIPfreeBufferArray(scip, &lincoefs);
            SCIPfreeBufferArray(scip, &children);
            SCIPfreeBufferArray(scip, &quadelems);

            SCIP_CALL_TERMINATE( retcode, SCIPcreateConsNonlinear(scip, cons, SCIPconsGetName(*cons),
                  SCIPgetNLinearVarsNonlinear(scip, *cons), SCIPgetLinearVarsNonlinear(scip, *cons),
                  SCIPgetLinearCoefsNonlinear(scip, *cons), 2, exprtrees, exprcoefs,
                  SCIPgetLhsNonlinear(scip, *cons), SCIPgetRhsNonlinear(scip, *cons),
                  SCIPconsIsInitial(*cons), SCIPconsIsSeparated(*cons), SCIPconsIsEnforced(*cons),
                  SCIPconsIsChecked(*cons), SCIPconsIsPropagated(*cons), SCIPconsIsLocal(*cons),
                  SCIPconsIsModifiable(*cons), SCIPconsIsDynamic(*cons), SCIPconsIsRemovable(*cons),
                  SCIPconsIsStickingAtNode(*cons)), TERMINATE );

            SCIP_CALL( SCIPreleaseCons(scip, &oldcons) );

            break;
         }

         case NONLINEAR:
         {
            SCIP_Real one;

            one = 1.0;
            SCIP_CALL( SCIPaddExprtreesNonlinear(scip, *cons, 1, &exprtree, &one) );
            break;
         }
         }

         *constype = NONLINEAR;
      }

   TERMINATE:
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );

      if( retcode != SCIP_OKAY )
         break;
   }

   SCIPfreeBufferArray(scip, &exprvars);
   SCIPfreeBufferArray(scip, &exprvaridx);

   SCIP_CALL( retcode );

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
   const char* name;
   XML_NODE* start;
   const XML_NODE* header;
   const XML_NODE* data;
   SCIP_RETCODE retcode;
   SCIP_Bool doingfine;
   SCIP_Bool initialconss;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamiccols;
   SCIP_Bool dynamicrows;
   SCIP_VAR** vars;
   int nvars;
   SCIP_CONS** conss;
   CONSTYPE* constypes;
   int nconss;
   SCIP_CONS* objcons;
   CONSTYPE objconstype;
   int i;

   assert(scip != NULL);
   assert(reader != NULL);
   assert(result != NULL);
   assert(filename != NULL);

   *result = SCIP_DIDNOTRUN;
   retcode = SCIP_READERROR;
   doingfine = TRUE;
   vars = NULL;
   nvars = 0;
   conss = NULL;
   constypes = NULL;
   nconss = 0;
   objcons = NULL;

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

   /* read constraint data (names, constants, lhs/rhs) */
   SCIP_CALL_TERMINATE( retcode, readConstraints(scip, data, &conss, &constypes, &nconss, initialconss, dynamicconss, dynamicrows, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;
   assert(conss != NULL || nconss == 0);

   /* read linear coefficients matrix */
   SCIP_CALL_TERMINATE( retcode, readLinearCoefs(scip, data, vars, nvars, conss, constypes, nconss, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* read quadratic coefficients (turns linear constraints into quadratic ones, may create objcons) */
   SCIP_CALL_TERMINATE( retcode, readQuadraticCoefs(scip, data, vars, nvars, conss, constypes, nconss, &objcons, &objconstype, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* read nonlinear expressions (turns constraints into nonlinear ones, may create objcons) */
   SCIP_CALL_TERMINATE( retcode, readNonlinearExprs(scip, data, vars, nvars, conss, constypes, nconss, &objcons, &objconstype, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;

   /* add constraints to problem */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);  /*lint !e613*/
      SCIP_CALL( SCIPaddCons(scip, conss[i]) );  /*lint !e613*/
   }
   if( objcons != NULL )
   {
      SCIP_CALL( SCIPaddCons(scip, objcons) );
   }

   /* read sos2 constraints  and add to problem*/
   SCIP_CALL_TERMINATE( retcode, readSOScons(scip, data, vars, nvars, initialconss, dynamicconss, dynamicrows, &doingfine), CLEANUP );
   if( !doingfine )
      goto CLEANUP;


   *result = SCIP_SUCCESS;
   retcode = SCIP_OKAY;

 CLEANUP:
   /* free xml data */
   if( start != NULL )
      xmlFreeNode(start);

   /* free variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );  /*lint !e613*/
   }
   SCIPfreeBufferArrayNull(scip, &vars);

   /* free constraints */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conss[i]) );  /*lint !e613*/
   }
   SCIPfreeBufferArrayNull(scip, &conss);
   SCIPfreeBufferArrayNull(scip, &constypes);

   if( objcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &objcons) );
   }

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
