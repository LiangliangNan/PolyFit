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

/**@file   reader_cnf.c
 * @brief  CNF file reader
 * @author Thorsten Koch
 * @author Tobias Achterberg
 *
 * The DIMACS CNF (conjunctive normal form) is a file format used for example for SAT problems. For a detailed description of
 * this format see http://people.sc.fsu.edu/~jburkardt/data/cnf/cnf.html .
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/reader_cnf.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/pub_misc.h"


#define READER_NAME             "cnfreader"
#define READER_DESC             "file reader for SAT problems in conjunctive normal form"
#define READER_EXTENSION        "cnf"

#define MAXLINELEN       65536


/*
 * cnf reader internal methods
 */

static
void readError(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   linecount,          /**< line number of error */
   const char*           errormsg            /**< error message */
   )
{
   SCIPerrorMessage("read error in line <%d>: %s\n", linecount, errormsg);
}

static
void readWarning(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   linecount,          /**< line number of error */
   const char*           warningmsg          /**< warning message */
   )
{
   SCIPwarningMessage(scip, "Line <%d>: %s\n", linecount, warningmsg);
}

/** reads the next non-empty non-comment line of a cnf file */
static
SCIP_RETCODE readCnfLine(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            file,               /**< input file */
   char*                 buffer,             /**< buffer for storing the input line */
   int                   size,               /**< size of the buffer */
   int*                  linecount           /**< pointer to the line number counter */
   )
{
   char* line;
   int linelen;

   assert(file != NULL);
   assert(buffer != NULL);
   assert(size >= 2);
   assert(linecount != NULL);

   do
   {
      (*linecount)++;
      line = SCIPfgets(buffer, size, file);
      if( line != NULL )
      {
         linelen = (int)strlen(line);
         if( linelen == size-1 )
         {
            char s[SCIP_MAXSTRLEN];
            (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "line too long (exceeds %d characters)", size-2);
            readError(scip, *linecount, s);
            return SCIP_READERROR;
         }
      }
      else
         linelen = 0;
   }
   while( line != NULL && (*line == 'c' || *line == '\n') );

   if( line != NULL && linelen >= 2 && line[linelen-2] == '\n' )
      line[linelen-2] = '\0';
   else if( linelen == 0 )
      *buffer = '\0';

   assert((line == NULL) == (*buffer == '\0'));

   return SCIP_OKAY;
}

/* Read SAT formula in "CNF File Format".
 * 
 *  The specification is taken from the
 *
 *  Satisfiability Suggested Format
 *
 *  Online available at http://www.intellektik.informatik.tu-darmstadt.de/SATLIB/Benchmarks/SAT/satformat.ps
 *
 *  The method reads all files of CNF format. Other formats (SAT, SATX, SATE) are not supported.
 */  
static
SCIP_RETCODE readCnf(
   SCIP*                 scip,               /**< SCIP data structure */   
   SCIP_FILE*            file                /**< input file */
   )
{
   SCIP_RETCODE retcode;
   SCIP_VAR** vars;
   SCIP_VAR** clausevars;
   SCIP_CONS* cons;
   int* varsign;
   char* tok;
   char* nexttok;
   char line[MAXLINELEN];
   char format[SCIP_MAXSTRLEN];
   char varname[SCIP_MAXSTRLEN];
   char s[SCIP_MAXSTRLEN];
   SCIP_Bool initialconss;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamiccols;
   SCIP_Bool dynamicrows;
   SCIP_Bool useobj;
   int linecount;
   int clauselen;
   int clausenum;
   int nvars;
   int nclauses;
   int varnum;
   int v;

   assert(scip != NULL);
   assert(file != NULL);

   linecount = 0;

   /* read header */
   SCIP_CALL( readCnfLine(scip, file, line, (int) sizeof(line), &linecount) );
   if( *line != 'p' )
   {
      readError(scip, linecount, "problem declaration line expected");
      return SCIP_READERROR;
   }
   /* cppcheck-suppress invalidScanfFormatWidth_smaller */
   if( sscanf(line, "p %8s %d %d", format, &nvars, &nclauses) != 3 )
   {
      readError(scip, linecount, "invalid problem declaration (must be 'p cnf <nvars> <nclauses>')");
      return SCIP_READERROR;
   }
   if( strcmp(format, "cnf") != 0 )
   {
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "invalid format tag <%s> (must be 'cnf')", format);
      readError(scip, linecount, s);
      return SCIP_READERROR;
   }
   if( nvars <= 0 )
   {
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "invalid number of variables <%d> (must be positive)", nvars);
      readError(scip, linecount, s);
      return SCIP_READERROR;
   }
   if( nclauses <= 0 )
   {
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "invalid number of clauses <%d> (must be positive)", nclauses);
      readError(scip, linecount, s);
      return SCIP_READERROR;
   }

   /* get parameter values */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/initialconss", &initialconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamiccols", &dynamiccols) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicrows", &dynamicrows) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/cnfreader/useobj", &useobj) );

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &clausevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, nvars) );

   /* create the variables */
   for( v = 0; v < nvars; ++v )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x%d", v+1);
      SCIP_CALL( SCIPcreateVar(scip, &vars[v], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, !dynamiccols, dynamiccols,
            NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, vars[v]) );
      varsign[v] = 0;
   }

   /* read clauses */
   clausenum = 0;
   clauselen = 0;
   do
   {
      retcode = readCnfLine(scip, file, line, (int) sizeof(line), &linecount);
      if( retcode != SCIP_OKAY )
         goto TERMINATE;

      if( *line != '\0' && *line != '%' )
      {
         tok = SCIPstrtok(line, " \f\n\r\t", &nexttok);
         while( tok != NULL )
         {
            /* parse literal and check for errors */
            if( sscanf(tok, "%d", &v) != 1 )
            {
               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "invalid literal <%s>", tok);
               readError(scip, linecount, s);
               retcode = SCIP_READERROR;
               goto TERMINATE;
            }

            /* interpret literal number: v == 0: end of clause, v < 0: negated literal, v > 0: positive literal */
            if( v == 0 )
            {
               /* end of clause: construct clause and add it to SCIP */
               if( clauselen == 0 )
                  readWarning(scip, linecount, "empty clause detected in line -- problem infeasible");

               clausenum++;
               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "c%d", clausenum);

               if( SCIPfindConshdlr(scip, "logicor") != NULL )
               {
                  /* if the constraint handler logicor exit create a logicor constraint */
                  SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, s, clauselen, clausevars,
                        initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows, FALSE) );
               }
               else if( SCIPfindConshdlr(scip, "setppc") != NULL )
               {
                  /* if the constraint handler logicor does not exit but constraint
                   *  handler setppc create a setppc constraint */
                  SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, s, clauselen, clausevars,
                        initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows, FALSE) );
               }
               else
               {
                  /* if none of the previous constraint handler exits create a linear
                   * constraint */
                  SCIP_Real* vals;
                  int i;

                  SCIP_CALL( SCIPallocBufferArray(scip, &vals, clauselen) );

                  for( i = 0; i < clauselen; ++i )
                     vals[i] = 1.0;

                  SCIP_CALL( SCIPcreateConsLinear(scip, &cons, s, clauselen, clausevars, vals, 1.0, SCIPinfinity(scip),
                        initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows, FALSE) );

                  SCIPfreeBufferArray(scip, &vals);
               }

               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               clauselen = 0;
            }
            else if( v >= -nvars && v <= nvars )
            {
               if( clauselen >= nvars )
               {
                  readError(scip, linecount, "too many literals in clause");
                  retcode = SCIP_READERROR;
                  goto TERMINATE;
               }

               /* add literal to clause */
               varnum = ABS(v)-1;
               if( v < 0 )
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, vars[varnum], &clausevars[clauselen]) );
                  varsign[varnum]--;
               }
               else
               {
                  clausevars[clauselen] = vars[varnum];
                  varsign[varnum]++;
               }
               clauselen++;
            }
            else
            {
               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "invalid variable number <%d>", ABS(v));
               readError(scip, linecount, s);
               retcode = SCIP_READERROR;
               goto TERMINATE;
            }

            /* get next token */
            tok = SCIPstrtok(NULL, " \f\n\r\t", &nexttok);
         }
      }
   }
   while( *line != '\0' && *line != '%' );

   /* check for additional literals */
   if( clauselen > 0 )
   {
      SCIPwarningMessage(scip, "found %d additional literals after last clause\n", clauselen);
   }

   /* check number of clauses */
   if( clausenum != nclauses )
   {
      SCIPwarningMessage(scip, "expected %d clauses, but found %d\n", nclauses, clausenum);
   }

 TERMINATE:
   /* change objective values and release variables */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   for( v = 0; v < nvars; ++v )
   {
      if( useobj )
      {
         SCIP_CALL( SCIPchgVarObj(scip, vars[v], (SCIP_Real)varsign[v]) );
      }
      SCIP_CALL( SCIPreleaseVar(scip, &vars[v]) );
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &varsign);
   SCIPfreeBufferArray(scip, &clausevars);
   SCIPfreeBufferArray(scip, &vars);

   return retcode;
}


/*
 * Callback methods
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCnf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderCnf(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCnf)
{  /*lint --e{715}*/
   SCIP_FILE* f;
   SCIP_RETCODE retcode;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(filename != NULL);
   assert(result != NULL);

   /* open file */
   f = SCIPfopen(filename, "r");
   if( f == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* read cnf file */
   retcode = readCnf(scip, f);

   /* close file */
   SCIPfclose(f);

   *result = SCIP_SUCCESS;

   return retcode;
}


/*
 * cnf file reader specific interface methods
 */

/** includes the cnf file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCnf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCnf) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCnf) );

   /* add cnf reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/cnfreader/useobj", "should an artificial objective, depending on the number of clauses a variable appears in, be used?",
         NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}

