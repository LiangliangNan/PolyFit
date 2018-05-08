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

/**@file   scipshell.c
 * @brief  SCIP command line interface
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"
#include "scip/message_default.h"

/*
 * Message Handler
 */

static
SCIP_RETCODE readParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< parameter file name */
   )
{
   if( SCIPfileExists(filename) )
   {
      SCIPinfoMessage(scip, NULL, "reading user parameter file <%s>\n", filename);
      SCIP_CALL( SCIPreadParams(scip, filename) );
   }
   else
      SCIPinfoMessage(scip, NULL, "user parameter file <%s> not found - using default parameters\n", filename);

   return SCIP_OKAY;
}

static
SCIP_RETCODE fromCommandLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< input file name */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Bool outputorigsol = FALSE;

   /********************
    * Problem Creation *
    ********************/

   /** @note The message handler should be only fed line by line such the message has the chance to add string in front
    *        of each message
    */
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "read problem <%s>\n", filename);
   SCIPinfoMessage(scip, NULL, "============\n");
   SCIPinfoMessage(scip, NULL, "\n");


   retcode = SCIPreadProb(scip, filename, NULL);

   switch( retcode )
   {
   case SCIP_NOFILE:
      SCIPinfoMessage(scip, NULL, "file <%s> not found\n", filename);
      return SCIP_OKAY;
   case SCIP_PLUGINNOTFOUND:
      SCIPinfoMessage(scip, NULL, "no reader for input file <%s> available\n", filename);
      return SCIP_OKAY;
   case SCIP_READERROR:
      SCIPinfoMessage(scip, NULL, "error reading file <%s>\n", filename);
      return SCIP_OKAY;
   default:
      SCIP_CALL( retcode );
   } /*lint !e788*/

   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   SCIPinfoMessage(scip, NULL, "\nsolve problem\n");
   SCIPinfoMessage(scip, NULL, "=============\n\n");

   SCIP_CALL( SCIPsolve(scip) );

   /*******************
    * Solution Output *
    *******************/

   SCIP_CALL( SCIPgetBoolParam(scip, "misc/outputorigsol", &outputorigsol) );
   if ( outputorigsol )
   {
      SCIP_SOL* bestsol;

      SCIPinfoMessage(scip, NULL, "\nprimal solution (original space):\n");
      SCIPinfoMessage(scip, NULL, "=================================\n\n");

      bestsol = SCIPgetBestSol(scip);
      if ( bestsol == NULL )
         SCIPinfoMessage(scip, NULL, "no solution available\n");
      else
      {
         SCIP_SOL* origsol;

         SCIP_CALL( SCIPcreateSolCopy(scip, &origsol, bestsol) );
         SCIP_CALL( SCIPretransformSol(scip, origsol) );
         SCIP_CALL( SCIPprintSol(scip, origsol, NULL, FALSE) );
         SCIP_CALL( SCIPfreeSol(scip, &origsol) );
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "\nprimal solution (transformed space):\n");
      SCIPinfoMessage(scip, NULL, "====================================\n\n");

      SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
   }


   /**************
    * Statistics *
    **************/

   SCIPinfoMessage(scip, NULL, "\nStatistics\n");
   SCIPinfoMessage(scip, NULL, "==========\n\n");

   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

/** evaluates command line parameters and runs SCIP appropriately in the given SCIP instance */
SCIP_RETCODE SCIPprocessShellArguments(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{  /*lint --e{850}*/
   char* probname = NULL;
   char* settingsname = NULL;
   char* logname = NULL;
   int randomseed;
   SCIP_Bool randomseedread;
   SCIP_Bool quiet;
   SCIP_Bool paramerror;
   SCIP_Bool interactive;
   SCIP_Bool onlyversion;
   SCIP_Real primalreference = SCIP_UNKNOWN;
   SCIP_Real dualreference = SCIP_UNKNOWN;
   const char* dualrefstring;
   const char* primalrefstring;
   int i;

   /********************
    * Parse parameters *
    ********************/

   quiet = FALSE;
   paramerror = FALSE;
   interactive = FALSE;
   onlyversion = FALSE;
   randomseedread = FALSE;
   randomseed = 0;
   primalrefstring = NULL;
   dualrefstring = NULL;

   for( i = 1; i < argc; ++i )
   {
      if( strcmp(argv[i], "-l") == 0 )
      {
         i++;
         if( i < argc )
            logname = argv[i];
         else
         {
            printf("missing log filename after parameter '-l'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-q") == 0 )
         quiet = TRUE;
      else if( strcmp(argv[i], "-v") == 0 )
         onlyversion = TRUE;
      else if( strcmp(argv[i], "--version") == 0 )
         onlyversion = TRUE;
      else if( strcmp(argv[i], "-s") == 0 )
      {
         i++;
         if( i < argc )
            settingsname = argv[i];
         else
         {
            printf("missing settings filename after parameter '-s'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-f") == 0 )
      {
         i++;
         if( i < argc )
            probname = argv[i];
         else
         {
            printf("missing problem filename after parameter '-f'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-c") == 0 )
      {
         i++;
         if( i < argc )
         {
            SCIP_CALL( SCIPaddDialogInputLine(scip, argv[i]) );
            interactive = TRUE;
         }
         else
         {
            printf("missing command line after parameter '-c'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-b") == 0 )
      {
         i++;
         if( i < argc )
         {
            SCIP_FILE* file;

            file = SCIPfopen(argv[i], "r");
            if( file == NULL )
            {
               printf("cannot read command batch file <%s>\n", argv[i]);
               SCIPprintSysError(argv[i]);
               paramerror = TRUE;
            }
            else
            {
               while( !SCIPfeof(file) )
               {
                  char buffer[SCIP_MAXSTRLEN];

                  (void)SCIPfgets(buffer, (int) sizeof(buffer), file);
                  if( buffer[0] != '\0' )
                  {
                     SCIP_CALL_FINALLY( SCIPaddDialogInputLine(scip, buffer), SCIPfclose(file) );
                  }
               }
               SCIPfclose(file);
               interactive = TRUE;
            }
         }
         else
         {
            printf("missing command batch filename after parameter '-b'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-r") == 0 )
      {
         /*read a random seed from the command line */
         i++;
         if( i < argc && isdigit(argv[i][0]) )
         {
            randomseed = atoi(argv[i]);
            randomseedread = TRUE;
         }
         else
         {
            printf("Random seed parameter '-r' followed by something that is not an integer\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-o") == 0 )
      {
         if( i >= argc - 2 )
         {
            printf("wrong usage of reference objective parameter '-o': -o <primref> <dualref>\n");
            paramerror = TRUE;
         }
         else
         {
            /* do not parse the strings directly, the settings could still influence the value of +-infinity */
            primalrefstring = argv[i + 1];
            dualrefstring = argv[i+2];
         }
         i += 2;
      }
      else
      {
         printf("invalid parameter <%s>\n", argv[i]);
         paramerror = TRUE;
      }
   }

   if( interactive && probname != NULL )
   {
      printf("cannot mix batch mode '-c' and '-b' with file mode '-f'\n");
      paramerror = TRUE;
   }

   if( !paramerror )
   {
      /***********************************
       * create log file message handler *
       ***********************************/

      if( quiet )
      {
         SCIPsetMessagehdlrQuiet(scip, quiet);
      }

      if( logname != NULL )
      {
         SCIPsetMessagehdlrLogfile(scip, logname);
      }

      /***********************************
       * Version and library information *
       ***********************************/

      SCIPprintVersion(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      SCIPprintExternalCodes(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      if( onlyversion )
      {
         SCIPprintBuildOptions(scip, NULL);
         SCIPinfoMessage(scip, NULL, "\n");
         return SCIP_OKAY;
      }

      /*****************
       * Load settings *
       *****************/

      if( settingsname != NULL )
      {
         SCIP_CALL( readParams(scip, settingsname) );
      }
      else if( defaultsetname != NULL )
      {
         SCIP_CALL( readParams(scip, defaultsetname) );
      }

      /************************************
       * Change random seed, if specified *
       ***********************************/
      if( randomseedread )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "randomization/randomseedshift", randomseed) );
      }

      /**************
       * Start SCIP *
       **************/

      if( probname != NULL )
      {
         SCIP_Bool validatesolve = FALSE;

         if( primalrefstring != NULL && dualrefstring != NULL )
         {
            char *endptr;
            if( ! SCIPparseReal(scip, primalrefstring, &primalreference, &endptr) ||
                     ! SCIPparseReal(scip, dualrefstring, &dualreference, &endptr) )
            {
               printf("error parsing primal and dual reference values for validation: %s %s\n", primalrefstring, dualrefstring);
               return SCIP_ERROR;
            }
            else
               validatesolve = TRUE;
         }
         SCIP_CALL( fromCommandLine(scip, probname) );

         /* validate the solve */
         if( validatesolve )
         {
            SCIP_CALL( SCIPvalidateSolve(scip, primalreference, dualreference, SCIPfeastol(scip), FALSE, NULL, NULL, NULL) );
         }
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "\n");
         SCIP_CALL( SCIPstartInteraction(scip) );
      }
   }
   else
   {
      printf("\nsyntax: %s [-l <logfile>] [-q] [-s <settings>] [-r <randseed>] [-f <problem>] [-b <batchfile>] [-c \"command\"]\n"
         "  -v, --version : print version and build options\n"
         "  -l <logfile>  : copy output into log file\n"
         "  -q            : suppress screen messages\n"
         "  -s <settings> : load parameter settings (.set) file\n"
         "  -f <problem>  : load and solve problem file\n"
         "  -o <primref> <dualref> : pass primal and dual objective reference values for validation at the end of the solve\n"
         "  -b <batchfile>: load and execute dialog command batch file (can be used multiple times)\n"
         "  -r <randseed> : nonnegative integer to be used as random seed. "
         "Has priority over random seed specified through parameter settings (.set) file\n"
         "  -c \"command\"  : execute single line of dialog commands (can be used multiple times)\n\n",
         argv[0]);
   }

   return SCIP_OKAY;
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
SCIP_RETCODE SCIPrunShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;
   SCIP_RETCODE retcode = SCIP_OKAY;
   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   /* include default SCIP plugins */
   SCIP_CALL_TERMINATE( retcode, SCIPincludeDefaultPlugins(scip), TERMINATE );

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL_TERMINATE( retcode, SCIPprocessShellArguments(scip, argc, argv, defaultsetname), TERMINATE );


   /********************
    * Deinitialization *
    ********************/
TERMINATE:
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return retcode;
}
