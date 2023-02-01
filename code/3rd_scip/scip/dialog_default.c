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

/**@file   dialog_default.c
 * @brief  default user interface dialog
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/dialog_default.h"
#include "nlpi/nlpi.h"
#include "scip/pub_cons.h"
#include "scip/type_cons.h"
#include "scip/cons_linear.h"



/** executes a menu dialog */
static
SCIP_RETCODE dialogExecMenu(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog,             /**< dialog menu */
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG**         nextdialog          /**< pointer to store next dialog to execute */
   )
{
   char* command;
   SCIP_Bool again;
   SCIP_Bool endoffile;
   int nfound;

   do
   {
      again = FALSE;

      /* get the next word of the command string */
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, NULL, &command, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }

      /* exit to the root dialog, if command is empty */
      if( command[0] == '\0' )
      {
         *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
         return SCIP_OKAY;
      }
      else if( strcmp(command, "..") == 0 )
      {
         *nextdialog = SCIPdialogGetParent(dialog);
         if( *nextdialog == NULL )
            *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
         return SCIP_OKAY;
      }

      /* find command in dialog */
      nfound = SCIPdialogFindEntry(dialog, command, nextdialog);

      /* check result */
      if( nfound == 0 )
      {
         SCIPdialogMessage(scip, NULL, "command <%s> not available\n", command);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
         *nextdialog = dialog;
      }
      else if( nfound >= 2 )
      {
         SCIPdialogMessage(scip, NULL, "\npossible completions:\n");
         SCIP_CALL( SCIPdialogDisplayCompletions(dialog, scip, command) );
         SCIPdialogMessage(scip, NULL, "\n");
         SCIPdialoghdlrClearBuffer(dialoghdlr);
         again = TRUE;
      }
   }
   while( again );

   return SCIP_OKAY;
}


/* parse the given string to detect a Boolean value and returns it */
static
SCIP_Bool parseBoolValue(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           valuestr,           /**< string to parse */
   SCIP_Bool*            error               /**< pointer to store the error result */
   )
{
   assert( scip  != NULL );
   assert( valuestr != NULL );
   assert( error != NULL );

   *error = FALSE;

   switch( valuestr[0] )
   {
   case 'f':
   case 'F':
   case '0':
   case 'n':
   case 'N':
      return FALSE;
   case 't':
   case 'T':
   case '1':
   case 'y':
   case 'Y':
      return TRUE;
   default:
      *error = TRUE;
      break;
   }

   return FALSE;
}


/* display the reader information */
static
void displayReaders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             reader,             /**< display reader which can read */
   SCIP_Bool             writer              /**< display reader which can write */
   )
{
   SCIP_READER** readers;
   int nreaders;
   int r;

   assert( scip != NULL );

   readers = SCIPgetReaders(scip);
   nreaders = SCIPgetNReaders(scip);

   /* display list of readers */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " file reader          extension  description\n");
   SCIPdialogMessage(scip, NULL, " -----------          ---------  -----------\n");
   for( r = 0; r < nreaders; ++r )
   {
      if( (reader && SCIPreaderCanRead(readers[r])) || (writer && SCIPreaderCanWrite(readers[r])) )
      {
         SCIPdialogMessage(scip, NULL, " %-20s ", SCIPreaderGetName(readers[r]));
         if( strlen(SCIPreaderGetName(readers[r])) > 20 )
            SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
         SCIPdialogMessage(scip, NULL, "%9s  ", SCIPreaderGetExtension(readers[r]));
         SCIPdialogMessage(scip, NULL, "%s", SCIPreaderGetDesc(readers[r]));
         SCIPdialogMessage(scip, NULL, "\n");
      }
   }
   SCIPdialogMessage(scip, NULL, "\n");
}


/* writes problem to file */
static
SCIP_RETCODE writeProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog,             /**< dialog menu */
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG**         nextdialog,         /**< pointer to store next dialog to execute */
   SCIP_Bool             transformed,        /**< output the transformed problem? */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   char* filename;
   SCIP_Bool endoffile;
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   if( filename[0] != '\0' )
   {
      char* tmpfilename;
      char* extension;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      /* copy filename */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpfilename, filename, (int)strlen(filename)+1) );
      extension = NULL;

      do
      {
         if( transformed )
            retcode = SCIPwriteTransProblem(scip, tmpfilename, extension, genericnames);
         else
            retcode = SCIPwriteOrigProblem(scip, tmpfilename, extension, genericnames);

         if( retcode == SCIP_FILECREATEERROR )
         {
            SCIPdialogMessage(scip, NULL, "error creating the file <%s>\n", filename);
            SCIPdialoghdlrClearBuffer(dialoghdlr);
            break;
         }         
         else if(retcode == SCIP_WRITEERROR )
         {
            SCIPdialogMessage(scip, NULL, "error writing file <%s>\n", filename);
            SCIPdialoghdlrClearBuffer(dialoghdlr);
            break;
         }
         else if( retcode == SCIP_PLUGINNOTFOUND )
         {
            /* ask user once for a suitable reader */
            if( extension == NULL )
            {
               SCIPdialogMessage(scip, NULL, "no reader for requested output format\n");

               SCIPdialogMessage(scip, NULL, "The following readers are available for writing:\n");
               displayReaders(scip, FALSE, TRUE);

               SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, 
                     "select a suitable reader by extension (or return): ", &extension, &endoffile) );

               if( extension[0] == '\0' )
                  break;
            }
            else
            {
               SCIPdialogMessage(scip, NULL, "no reader for output in <%s> format\n", extension);
               extension = NULL;
            }
         }
         else
         {
            /* check for unexpected errors */
            SCIP_CALL( retcode );

            /* print result message if writing was successful */
            if( transformed )
               SCIPdialogMessage(scip, NULL, "written transformed problem to file <%s>\n", tmpfilename);
            else
               SCIPdialogMessage(scip, NULL, "written original problem to file <%s>\n", tmpfilename);
            break;
         }
      }
      while( extension != NULL );

      SCIPfreeBufferArray(scip, &tmpfilename);
   }

   return SCIP_OKAY;
}

/** copy method for dialog plugins (called when SCIP copies plugins) */
static
SCIP_DECL_DIALOGCOPY(dialogCopyDefault)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(dialog != NULL);

   /* call inclusion method of dialog */
   SCIP_CALL( SCIPincludeDialogDefault(scip) );

   return SCIP_OKAY;
}

/** standard menu dialog execution method, that displays it's help screen if the remaining command line is empty */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecMenu)
{  /*lint --e{715}*/
   /* if remaining command string is empty, display menu of available options */
   if( SCIPdialoghdlrIsBufferEmpty(dialoghdlr) )
   {
      SCIPdialogMessage(scip, NULL, "\n");
      SCIP_CALL( SCIPdialogDisplayMenu(dialog, scip) );
      SCIPdialogMessage(scip, NULL, "\n");
   }

   SCIP_CALL( dialogExecMenu(scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}

/** standard menu dialog execution method, that doesn't display it's help screen */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecMenuLazy)
{  /*lint --e{715}*/
   SCIP_CALL( dialogExecMenu(scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}

/** dialog execution method for the change add constraint */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeAddCons)
{  /*lint --e{715}*/

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM )
      SCIPdialogMessage(scip, NULL, "cannot call method after problem was transformed\n");
   else if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
      SCIPdialogMessage(scip, NULL, "cannot call method before problem was created\n");
   else
   {
      SCIP_CONS* cons;
      SCIP_Bool endoffile;
      char* str;

      cons = NULL;

      SCIP_CALL( SCIPdialoghdlrGetLine(dialoghdlr, dialog, "write constraint in <cip> format\n", &str, &endoffile) );

      if( str[0] != '\0' )
      {
         SCIP_Bool success;

         printf("<%s>\n", str);

         SCIP_CALL( SCIPparseCons(scip, &cons, str, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

         if( success )
         {
            char consstr[SCIP_MAXSTRLEN];

            /* add and release constraint */
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            SCIPdialogMessage(scip, NULL, "successfully added constraint\n"); 
            SCIPescapeString(consstr, SCIP_MAXSTRLEN, str);

            SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, consstr, FALSE) );
         }
         else
         {
            SCIPdialogMessage(scip, NULL, "constraint was not recognizable\n");
         }
      }
   }

   /* set root dialog as next dialog */
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the change bounds command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeBounds)
{  /*lint --e{715}*/

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM )
      SCIPdialogMessage(scip, NULL, "cannot call method after problem was transformed\n");
   else if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
      SCIPdialogMessage(scip, NULL, "cannot call method before problem was created\n");
   else
   {
      SCIP_VAR* var;
      SCIP_Bool endoffile;
      char* varname;

      var = NULL;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

      do
      {
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter variable name: ", &varname, &endoffile) );

         /* if we get a return or we reached the end of the file, then we stop */
         if( varname[0] == '\0' || endoffile )
            break;

         var = SCIPfindVar(scip, varname);

         if( var == NULL )
            SCIPdialogMessage(scip, NULL, "variable <%s> does not exist\n", varname);
      }
      while( var == NULL );

      if( var != NULL )
      {
         do
         {
            char* boundstr;
            char message[SCIP_MAXSTRLEN];
            SCIP_Real bound;

            SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, varname, FALSE) );

            (void)SCIPsnprintf(message, SCIP_MAXSTRLEN, "current lower bound <%.15g> (Return to skip): ", SCIPvarGetLbGlobal(var));
            SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, message, &boundstr, &endoffile) );

            /* if we reached the end of the file, then we stop */
            if( endoffile )
               break;

            if( boundstr[0] != '\0' )
            {
               char* endptr;

               bound = strtod(boundstr, &endptr);
               if( endptr == boundstr || *endptr != '\0' )
               {
                  printf("<%s> <%s>\n", endptr, boundstr);
                  SCIPdialogMessage(scip, NULL, "ignore none value string\n");
               }
               else if( SCIPisGT(scip, bound, SCIPvarGetUbGlobal(var)) )
               {
                  SCIPdialogMessage(scip, NULL, "ignore lower bound <%.15g> since it is larger than the current upper bound <%.15g>\n",
                     bound, SCIPvarGetUbGlobal(var));
               }
               else
               {
                  SCIP_CALL( SCIPchgVarLbGlobal(scip, var, bound) );
               }
            }

            (void)SCIPsnprintf(message, SCIP_MAXSTRLEN, "current upper bound <%.15g> (Return to skip): ", SCIPvarGetUbGlobal(var));
            SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, message, &boundstr, &endoffile) );

            /* if we reached the end of the file, then we stop */
            if( endoffile )
               break;

            if( boundstr[0] != '\0' )
            {
               char* endptr;

               bound = strtod(boundstr, &endptr);
               if( endptr == boundstr || *endptr != '\0' )
               {
                  SCIPdialogMessage(scip, NULL, "ignore none value string\n");
               }
               else if( SCIPisLT(scip, bound, SCIPvarGetLbGlobal(var)) )
               {
                  SCIPdialogMessage(scip, NULL, "ignore new upper bound <%.15g> since it is smaller than the current lower bound <%.15g>\n",
                     bound, SCIPvarGetLbGlobal(var));
               }
               else
               {
                  SCIP_CALL( SCIPchgVarUbGlobal(scip, var, bound) );
               }
            }
         }
         while( FALSE);

         SCIPdialogMessage(scip, NULL, "variable <%s> global bounds [%.15g,%.15g]\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      }
   }

   /* set root dialog as next dialog */
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the freetransproblem command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeFreetransproblem)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   /* free transformed problem */
   SCIP_CALL( SCIPfreeTransform(scip) );

   /* set root dialog as next dialog */
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the changing the objective sense */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeObjSense)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM )
      SCIPdialogMessage(scip, NULL, "cannot call method after problem was transformed\n");
   else if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
      SCIPdialogMessage(scip, NULL, "cannot call method before problem was created\n");
   else
   {
      SCIP_Bool endoffile;
      char* objsense;

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "new objective sense {min,max}: ", &objsense, &endoffile) );

      /* if we get a return or we reached the end of the file, then we stop */
      if( objsense[0] != '\0' && !endoffile )
      {
         if( strncmp(objsense, "max", 3) == 0 )
         {
            SCIP_CALL( SCIPsetObjsense(scip,  SCIP_OBJSENSE_MAXIMIZE) );
         }
         else if( strncmp(objsense , "min", 3) == 0 )
         {
            SCIP_CALL( SCIPsetObjsense(scip,  SCIP_OBJSENSE_MINIMIZE) );
         }
         else
         {
            SCIPdialogMessage(scip, NULL, "invalid argument <%s>\n", objsense);
         }
      }
   }

   /* set root dialog as next dialog */
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the checksol command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChecksol)
{  /*lint --e{715}*/
   SCIP_SOL* sol;
   SCIP_Bool feasible;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
      sol = SCIPgetBestSol(scip);
   else
      sol = NULL;

   if( sol == NULL )
      SCIPdialogMessage(scip, NULL, "no feasible solution available\n");
   else
   {
      SCIP_Real oldfeastol;
      SCIP_Real checkfeastolfac;
      SCIP_Bool dispallviols;

      oldfeastol = SCIPfeastol(scip);
      SCIP_CALL( SCIPgetRealParam(scip, "numerics/checkfeastolfac", &checkfeastolfac) );
      SCIP_CALL( SCIPgetBoolParam(scip, "display/allviols", &dispallviols) );

      /* scale feasibility tolerance by set->num_checkfeastolfac */
      if( !SCIPisEQ(scip, checkfeastolfac, 1.0) )
      {
         SCIP_CALL( SCIPchgFeastol(scip, oldfeastol * checkfeastolfac) );
      }

      SCIPinfoMessage(scip, NULL, "check best solution\n");
      SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, dispallviols) );

      /* restore old feasibilty tolerance */
      if( !SCIPisEQ(scip, checkfeastolfac, 1.0) )
      {
         SCIP_CALL( SCIPchgFeastol(scip, oldfeastol) );
      }

      if( feasible )
         SCIPdialogMessage(scip, NULL, "solution is feasible in original problem\n");

      SCIPdialogMessage(scip, NULL, "%-19s: %11s %11s\n", "Violation", "absolute", "relative");
      SCIPdialogMessage(scip, NULL, "%-19s: %11.5e %11.5e\n", "  bounds", SCIPsolGetAbsBoundViolation(sol), SCIPsolGetRelBoundViolation(sol));
      SCIPdialogMessage(scip, NULL, "%-19s: %11.5e %11s\n", "  integrality", SCIPsolGetAbsIntegralityViolation(sol), "-");
      SCIPdialogMessage(scip, NULL, "%-19s: %11.5e %11.5e\n", "  LP rows", SCIPsolGetAbsLPRowViolation(sol), SCIPsolGetRelLPRowViolation(sol));
      SCIPdialogMessage(scip, NULL, "%-19s: %11.5e %11.5e\n", "  constraints", SCIPsolGetAbsConsViolation(sol), SCIPsolGetRelConsViolation(sol));
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the cliquegraph command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCliquegraph)
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;
   SCIP_Bool endoffile;
   char* filename;

   assert(nextdialog != NULL);

   *nextdialog = NULL;

   if( !SCIPisTransformed(scip) )
   {
      SCIPdialogMessage(scip, NULL, "cannot call method before problem was transformed\n");
      SCIPdialoghdlrClearBuffer(dialoghdlr);
   }
   else
   {
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }

      if( filename[0] != '\0' )
      {
         SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

         retcode = SCIPwriteCliqueGraph(scip, filename, FALSE);
         if( retcode == SCIP_FILECREATEERROR )
            SCIPdialogMessage(scip, NULL, "error creating file <%s>\n", filename);
         else
         {
            SCIP_CALL( retcode );
         }
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display branching command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayBranching)
{  /*lint --e{715}*/
   SCIP_BRANCHRULE** branchrules;
   SCIP_BRANCHRULE** sorted;
   int nbranchrules;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   branchrules = SCIPgetBranchrules(scip);
   nbranchrules = SCIPgetNBranchrules(scip);

   /* copy branchrules array into temporary memory for sorting */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sorted, branchrules, nbranchrules) );

   /* sort the branching rules */
   SCIPsortPtr((void**)sorted, SCIPbranchruleComp, nbranchrules);

   /* display sorted list of branching rules */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " branching rule       priority maxdepth maxbddist  description\n");
   SCIPdialogMessage(scip, NULL, " --------------       -------- -------- ---------  -----------\n");
   for( i = 0; i < nbranchrules; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPbranchruleGetName(sorted[i]));
      if( strlen(SCIPbranchruleGetName(sorted[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d %8d %8.1f%%  ", SCIPbranchruleGetPriority(sorted[i]),
         SCIPbranchruleGetMaxdepth(sorted[i]), 100.0 * SCIPbranchruleGetMaxbounddist(sorted[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPbranchruleGetDesc(sorted[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &sorted);

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display relaxators command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayRelaxators)
{  /*lint --e{715}*/
   SCIP_RELAX** relaxs;
   SCIP_RELAX** sorted;
   int nrelaxs;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   relaxs = SCIPgetRelaxs(scip);
   nrelaxs = SCIPgetNRelaxs(scip);

   /* copy relaxs array into temporary memory for sorting */
   if( nrelaxs != 0 )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &sorted, relaxs, nrelaxs) );
   }
   else
      sorted = NULL;

   /* sort the relaxators */
   SCIPsortPtr((void**)sorted, SCIPrelaxComp, nrelaxs);

   /* display sorted list of relaxators */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " relaxator            priority freq  description\n");
   SCIPdialogMessage(scip, NULL, " --------------       -------- ----  -----------\n");
   for( i = 0; i < nrelaxs; ++i )
   {
      assert(sorted != NULL); /* for flexelint */
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPrelaxGetName(sorted[i]));
      if( strlen(SCIPrelaxGetName(sorted[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d %4d  ", SCIPrelaxGetPriority(sorted[i]),
         SCIPrelaxGetFreq(sorted[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPrelaxGetDesc(sorted[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   /* free temporary memory */
   SCIPfreeBufferArrayNull(scip, &sorted);

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display conflict command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayConflict)
{  /*lint --e{715}*/
   SCIP_CONFLICTHDLR** conflicthdlrs;
   SCIP_CONFLICTHDLR** sorted;
   int nconflicthdlrs;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   conflicthdlrs = SCIPgetConflicthdlrs(scip);
   nconflicthdlrs = SCIPgetNConflicthdlrs(scip);

   /* copy conflicthdlrs array into temporary memory for sorting */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sorted, conflicthdlrs, nconflicthdlrs) );

   /* sort the conflict handlers */
   SCIPsortPtr((void**)sorted, SCIPconflicthdlrComp, nconflicthdlrs);

   /* display sorted list of conflict handlers */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " conflict handler     priority  description\n");
   SCIPdialogMessage(scip, NULL, " ----------------     --------  -----------\n");
   for( i = 0; i < nconflicthdlrs; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPconflicthdlrGetName(sorted[i]));
      if( strlen(SCIPconflicthdlrGetName(sorted[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d  ", SCIPconflicthdlrGetPriority(sorted[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPconflicthdlrGetDesc(sorted[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &sorted);

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display conshdlrs command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayConshdlrs)
{  /*lint --e{715}*/
   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* display list of constraint handlers */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " Legend:\n");
   SCIPdialogMessage(scip, NULL, "  prestim (presolve timing): 'f'ast, 'm'edium, 'e'xhaustive\n\n");
   SCIPdialogMessage(scip, NULL, " constraint handler   chckprio enfoprio sepaprio sepaf propf eager prestim description\n");
   SCIPdialogMessage(scip, NULL, " ------------------   -------- -------- -------- ----- ----- ----- ------- -----------\n");
   for( i = 0; i < nconshdlrs; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPconshdlrGetName(conshdlrs[i]));
      if( strlen(SCIPconshdlrGetName(conshdlrs[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d %8d %8d %5d %5d %5d  ",
         SCIPconshdlrGetCheckPriority(conshdlrs[i]),
         SCIPconshdlrGetEnfoPriority(conshdlrs[i]),
         SCIPconshdlrGetSepaPriority(conshdlrs[i]),
         SCIPconshdlrGetSepaFreq(conshdlrs[i]),
         SCIPconshdlrGetPropFreq(conshdlrs[i]),
         SCIPconshdlrGetEagerFreq(conshdlrs[i]));
      SCIPdialogMessage(scip, NULL, "   %c", (SCIPconshdlrGetPresolTiming(conshdlrs[i]) & SCIP_PRESOLTIMING_FAST) ? 'f' : ' ');
      SCIPdialogMessage(scip, NULL, "%c", (SCIPconshdlrGetPresolTiming(conshdlrs[i]) & SCIP_PRESOLTIMING_MEDIUM) ? 'm' : ' ');
      SCIPdialogMessage(scip, NULL, "%c  ", (SCIPconshdlrGetPresolTiming(conshdlrs[i]) & SCIP_PRESOLTIMING_EXHAUSTIVE) ? 'e' : ' ');
      SCIPdialogMessage(scip, NULL, "%s", SCIPconshdlrGetDesc(conshdlrs[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display displaycols command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayDisplaycols)
{  /*lint --e{715}*/
   SCIP_DISP** disps;
   int ndisps;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   disps = SCIPgetDisps(scip);
   ndisps = SCIPgetNDisps(scip);

   /* display list of display columns */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " display column       header           position width priority status  description\n");
   SCIPdialogMessage(scip, NULL, " --------------       ------           -------- ----- -------- ------  -----------\n");
   for( i = 0; i < ndisps; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPdispGetName(disps[i]));
      if( strlen(SCIPdispGetName(disps[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%-16s ", SCIPdispGetHeader(disps[i]));
      if( strlen(SCIPdispGetHeader(disps[i])) > 16 )
         SCIPdialogMessage(scip, NULL, "\n %20s %16s ", "", "-->");
      SCIPdialogMessage(scip, NULL, "%8d ", SCIPdispGetPosition(disps[i]));
      SCIPdialogMessage(scip, NULL, "%5d ", SCIPdispGetWidth(disps[i]));
      SCIPdialogMessage(scip, NULL, "%8d ", SCIPdispGetPriority(disps[i]));
      switch( SCIPdispGetStatus(disps[i]) )
      {
      case SCIP_DISPSTATUS_OFF:
         SCIPdialogMessage(scip, NULL, "%6s  ", "off");
         break;
      case SCIP_DISPSTATUS_AUTO:
         SCIPdialogMessage(scip, NULL, "%6s  ", "auto");
         break;
      case SCIP_DISPSTATUS_ON:
         SCIPdialogMessage(scip, NULL, "%6s  ", "on");
         break;
      default:
         SCIPdialogMessage(scip, NULL, "%6s  ", "?");
         break;
      }
      SCIPdialogMessage(scip, NULL, "%s", SCIPdispGetDesc(disps[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display heuristics command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayHeuristics)
{  /*lint --e{715}*/
   SCIP_HEUR** heurs;
   int nheurs;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   heurs = SCIPgetHeurs(scip);
   nheurs = SCIPgetNHeurs(scip);

   /* display list of primal heuristics */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " primal heuristic     c priority freq ofs  description\n");
   SCIPdialogMessage(scip, NULL, " ----------------     - -------- ---- ---  -----------\n");
   for( i = 0; i < nheurs; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPheurGetName(heurs[i]));
      if( strlen(SCIPheurGetName(heurs[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%c ", SCIPheurGetDispchar(heurs[i]));
      SCIPdialogMessage(scip, NULL, "%8d ", SCIPheurGetPriority(heurs[i]));
      SCIPdialogMessage(scip, NULL, "%4d ", SCIPheurGetFreq(heurs[i]));
      SCIPdialogMessage(scip, NULL, "%3d  ", SCIPheurGetFreqofs(heurs[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPheurGetDesc(heurs[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display memory command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayMemory)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   SCIPprintMemoryDiagnostic(scip);
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display nlpi command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayNlpi)
{  /*lint --e{715}*/
   SCIP_NLPI** nlpis;
   SCIP_NLPI** sorted;
   int nnlpis;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   nlpis  = SCIPgetNlpis(scip);
   nnlpis = SCIPgetNNlpis(scip);

   /* copy nlpis array into temporary memory for sorting */
   if( nnlpis != 0 )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &sorted, nlpis, nnlpis) );
   }
   else
      sorted = NULL;

   /* sort the branching rules */
   SCIPsortPtr((void**)sorted, SCIPnlpiComp, nnlpis);

   /* display sorted list of branching rules */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " NLP interface        priority description\n");
   SCIPdialogMessage(scip, NULL, " -------------        -------- -----------\n");
   for( i = 0; i < nnlpis; ++i )
   {
      assert(sorted != NULL);
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPnlpiGetName(sorted[i]));
      if( strlen(SCIPnlpiGetName(sorted[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d ", SCIPnlpiGetPriority(sorted[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPnlpiGetDesc(sorted[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   /* free temporary memory */
   if( nnlpis != 0 )
   {
      SCIPfreeBufferArray(scip, &sorted);
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display nodeselectors command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayNodeselectors)
{  /*lint --e{715}*/
   SCIP_NODESEL** nodesels;
   int nnodesels;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   nodesels = SCIPgetNodesels(scip);
   nnodesels = SCIPgetNNodesels(scip);

   /* display list of node selectors */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " node selector        std priority memsave prio  description\n");
   SCIPdialogMessage(scip, NULL, " -------------        ------------ ------------  -----------\n");
   for( i = 0; i < nnodesels; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPnodeselGetName(nodesels[i]));
      if( strlen(SCIPnodeselGetName(nodesels[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%12d ", SCIPnodeselGetStdPriority(nodesels[i]));
      SCIPdialogMessage(scip, NULL, "%12d  ", SCIPnodeselGetMemsavePriority(nodesels[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPnodeselGetDesc(nodesels[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display parameters command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayParameters)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, "number of parameters = %d\n", SCIPgetNParams(scip));
   SCIPdialogMessage(scip, NULL, "non-default parameter settings:\n");
   SCIP_CALL( SCIPwriteParams(scip, NULL, FALSE, TRUE) );
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display presolvers command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPresolvers)
{  /*lint --e{715}*/
   SCIP_PRESOL** presols;
   int npresols;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   presols = SCIPgetPresols(scip);
   npresols = SCIPgetNPresols(scip);

   /* display list of presolvers */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " Legend:\n");
   SCIPdialogMessage(scip, NULL, "  priority:  presolver called before constraint handlers iff priority > 0\n");
   SCIPdialogMessage(scip, NULL, "  timing:    'f'ast, 'm'edium, 'e'xhaustive\n\n");
   SCIPdialogMessage(scip, NULL, "  maxrounds: -1: no limit, 0: off, >0: limited number of rounds\n\n");
   SCIPdialogMessage(scip, NULL, " presolver            priority  timing  maxrounds  description\n");
   SCIPdialogMessage(scip, NULL, " ---------            --------  ------  ---------  -----------\n");
   for( i = 0; i < npresols; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPpresolGetName(presols[i]));
      if( strlen(SCIPpresolGetName(presols[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d  ", SCIPpresolGetPriority(presols[i]));
      SCIPdialogMessage(scip, NULL, "   %c", (SCIPpresolGetTiming(presols[i]) & SCIP_PRESOLTIMING_FAST) ? 'f' : ' ');
      SCIPdialogMessage(scip, NULL, "%c", (SCIPpresolGetTiming(presols[i]) & SCIP_PRESOLTIMING_MEDIUM) ? 'm' : ' ');
      SCIPdialogMessage(scip, NULL, "%c  ", (SCIPpresolGetTiming(presols[i]) & SCIP_PRESOLTIMING_EXHAUSTIVE) ? 'e' : ' ');
      SCIPdialogMessage(scip, NULL, "%9d  ", SCIPpresolGetMaxrounds(presols[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPpresolGetDesc(presols[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display pricer command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPricers)
{  /*lint --e{715}*/
   SCIP_PRICER** pricers;
   int npricers;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   pricers = SCIPgetPricers(scip);
   npricers = SCIPgetNPricers(scip);

   /* display list of pricers */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " pricer               priority  description\n");
   SCIPdialogMessage(scip, NULL, " ----------           --------  -----------\n");
   for( i = 0; i < npricers; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPpricerGetName(pricers[i]));
      if( strlen(SCIPpricerGetName(pricers[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d%c ", SCIPpricerGetPriority(pricers[i]), SCIPpricerIsDelayed(pricers[i]) ? 'd' : ' ');
      SCIPdialogMessage(scip, NULL, "%s", SCIPpricerGetDesc(pricers[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display problem command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayProblem)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM )
   {
      SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );
   }
   else
      SCIPdialogMessage(scip, NULL, "no problem available\n");

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display propagators command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPropagators)
{  /*lint --e{715}*/
   SCIP_PROP** props;
   int nprops;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   props = SCIPgetProps(scip);
   nprops = SCIPgetNProps(scip);

   /* display list of propagators */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " Legend:\n");
   SCIPdialogMessage(scip, NULL, "  presprio: propagator presolving called before constraint handlers iff presprio > 0\n");
   SCIPdialogMessage(scip, NULL, "  prestim (presolve timing): 'f'ast, 'm'edium, 'e'xhaustive\n\n");

   SCIPdialogMessage(scip, NULL, " propagator           propprio  freq  presprio  prestim   description\n");
   SCIPdialogMessage(scip, NULL, " ----------           --------  ----  --------  -------  -----------\n");
   for( i = 0; i < nprops; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPpropGetName(props[i]));
      if( strlen(SCIPpropGetName(props[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d%c ", SCIPpropGetPriority(props[i]), SCIPpropIsDelayed(props[i]) ? 'd' : ' ');
      SCIPdialogMessage(scip, NULL, "%4d  ", SCIPpropGetFreq(props[i]));
      SCIPdialogMessage(scip, NULL, "%8d  ", SCIPpropGetPresolPriority(props[i]));
      SCIPdialogMessage(scip, NULL, "    %c", (SCIPpropGetPresolTiming(props[i]) & SCIP_PRESOLTIMING_FAST) ? 'f' : ' ');
      SCIPdialogMessage(scip, NULL, "%c", (SCIPpropGetPresolTiming(props[i]) & SCIP_PRESOLTIMING_MEDIUM) ? 'm' : ' ');
      SCIPdialogMessage(scip, NULL, "%c  ", (SCIPpropGetPresolTiming(props[i]) & SCIP_PRESOLTIMING_EXHAUSTIVE) ? 'e' : ' ');
      SCIPdialogMessage(scip, NULL, "%s", SCIPpropGetDesc(props[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display readers command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayReaders)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   /* print reader information */
   displayReaders(scip, TRUE, TRUE);

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display separators command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySeparators)
{  /*lint --e{715}*/
   SCIP_SEPA** sepas;
   int nsepas;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   sepas = SCIPgetSepas(scip);
   nsepas = SCIPgetNSepas(scip);

   /* display list of separators */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " separator            priority  freq bddist  description\n");
   SCIPdialogMessage(scip, NULL, " ---------            --------  ---- ------  -----------\n");
   for( i = 0; i < nsepas; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-20s ", SCIPsepaGetName(sepas[i]));
      if( strlen(SCIPsepaGetName(sepas[i])) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d%c ", SCIPsepaGetPriority(sepas[i]), SCIPsepaIsDelayed(sepas[i]) ? 'd' : ' ');
      SCIPdialogMessage(scip, NULL, "%4d ", SCIPsepaGetFreq(sepas[i]));
      SCIPdialogMessage(scip, NULL, "%6.2f  ", SCIPsepaGetMaxbounddist(sepas[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPsepaGetDesc(sepas[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display solution command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySolution)
{  /*lint --e{715}*/
   SCIP_VAR** fixedvars;
   SCIP_VAR* var;
   SCIP_Bool printzeros;
   int nfixedvars;
   int v;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
      SCIPdialogMessage(scip, NULL, "No problem exists. Read (and solve) problem first.\n");
   else
   {
      SCIP_CALL( SCIPgetBoolParam(scip, "write/printzeros", &printzeros) );

      SCIPdialogMessage(scip, NULL, "\n");
      SCIP_CALL( SCIPprintBestSol(scip, NULL, printzeros) );
      SCIPdialogMessage(scip, NULL, "\n");

      /* check if there are infinite fixings and print a reference to 'display finitesolution', if needed */
      fixedvars = SCIPgetFixedVars(scip);
      nfixedvars = SCIPgetNFixedVars(scip);
      assert(fixedvars != NULL || nfixedvars == 0);

      /* check whether there are variables fixed to an infinite value */
      for( v = 0; v < nfixedvars; ++v )
      {
         var = fixedvars[v]; /*lint !e613*/

         /* skip (multi-)aggregated variables */
         if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED )
            continue;

         if( (SCIPisInfinity(scip, SCIPvarGetLbGlobal(var)) || SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var))) )
         {
            SCIPdialogMessage(scip, NULL, "The primal solution contains variables fixed to infinite values.\n\
If you want SCIP to display an optimal solution without infinite values, use 'display finitesolution'.\n");
            SCIPdialogMessage(scip, NULL, "\n");
            break;
         }
      }
   }
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display finitesolution command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayFiniteSolution)
{  /*lint --e{715}*/
   SCIP_SOL* bestsol = SCIPgetBestSol(scip);

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   if( bestsol != NULL )
   {
      SCIP_SOL* sol;
      SCIP_Bool success;
      SCIP_RETCODE retcode;

      /* create copy of solution with finite values */
      retcode = SCIPcreateFiniteSolCopy(scip, &sol, bestsol, &success);

      if( retcode == SCIP_OKAY && success )
      {
         SCIP_Bool printzeros;

         SCIP_CALL( SCIPgetBoolParam(scip, "write/printzeros", &printzeros) );
         retcode = SCIPprintSol(scip, sol, NULL, printzeros);
         SCIPdialogMessage(scip, NULL, "\n");
      }
      else
      {
         SCIPdialogMessage(scip, NULL, "error while creating finite solution\n");
      }

      /* free solution copy */
      if( retcode == SCIP_OKAY && sol != NULL )
      {
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      }
   }
   else
   {
      SCIP_Bool printzeros;

      SCIP_CALL( SCIPgetBoolParam(scip, "write/printzeros", &printzeros) );
      SCIP_CALL( SCIPprintBestSol(scip, NULL, printzeros) );
      SCIPdialogMessage(scip, NULL, "\n");
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display dual solution command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayDualSolution)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPprintDualSol(scip, NULL, FALSE) );
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}


/** dialog execution method for the display of solutions in the pool command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySolutionPool)
{  /*lint --e{715}*/
   char prompt[SCIP_MAXSTRLEN];
   SCIP_Bool endoffile;
   SCIP_SOL** sols;
   char* idxstr;
   char* endstr;
   int nsols;
   int idx;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
   SCIPdialogMessage(scip, NULL, "\n");

   if ( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPdialogMessage(scip, NULL, "No solution available.\n\n");
      return SCIP_OKAY;
   }

   nsols = SCIPgetNSols(scip);
   if ( nsols == 0 )
   {
      SCIPdialogMessage(scip, NULL, "No solution available.\n\n");
      return SCIP_OKAY;
   }

   /* parse solution number */
   (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "index of solution [0-%d]: ", nsols-1);

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &idxstr, &endoffile) );

   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   if ( SCIPstrToIntValue(idxstr, &idx, &endstr) )
   {
      SCIP_Bool printzeros;

      if ( idx < 0 || idx >= nsols )
      {
         SCIPdialogMessage(scip, NULL, "Solution index out of bounds [0-%d].\n", nsols-1);
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPgetBoolParam(scip, "write/printzeros", &printzeros) );

      sols = SCIPgetSols(scip);
      assert( sols[idx] != NULL );
      SCIP_CALL( SCIPprintSol(scip, sols[idx], NULL, FALSE) );
   }
   SCIPdialogMessage(scip, NULL, "\n");

   return SCIP_OKAY;
}

/** dialog execution method for the display statistics command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayStatistics)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display reoptstatistics command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayReoptStatistics)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPprintReoptStatistics(scip, NULL) );
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display compression command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayCompression)
{  /*lint --e{715}*/
   SCIP_COMPR** comprs;
   SCIP_COMPR** sorted;
   int ncomprs;
   int i;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   comprs = SCIPgetComprs(scip);
   ncomprs = SCIPgetNCompr(scip);

   /* copy compression array into temporary memory for sorting */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sorted, comprs, ncomprs) );

   /* sort the compression t */
   SCIPsortPtr((void**)sorted, SCIPcomprComp, ncomprs);

   /* display sorted list of branching rules */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " compression method       priority minnodes  description\n");
   SCIPdialogMessage(scip, NULL, " ------------------       -------- --------  -----------\n");
   for( i = 0; i < ncomprs; ++i )
   {
      SCIPdialogMessage(scip, NULL, " %-24s ", SCIPcomprGetName(sorted[i]));
      if( strlen(SCIPcomprGetName(sorted[i])) > 24 )
         SCIPdialogMessage(scip, NULL, "\n %24s ", "-->");
      SCIPdialogMessage(scip, NULL, "%8d %8d  ", SCIPcomprGetPriority(sorted[i]), SCIPcomprGetMinNodes(sorted[i]));
      SCIPdialogMessage(scip, NULL, "%s", SCIPcomprGetDesc(sorted[i]));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &sorted);

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display transproblem command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayTransproblem)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   if(SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED)
   {
      SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
   }
   else
      SCIPdialogMessage(scip, NULL, "no transformed problem available\n");

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display value command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayValue)
{  /*lint --e{715}*/
   SCIP_SOL* sol;
   SCIP_VAR* var;
   char* varname;
   SCIP_Real solval;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
      sol = SCIPgetBestSol(scip);
   else
      sol = NULL;

   if( sol == NULL )
   {
      SCIPdialogMessage(scip, NULL, "no feasible solution available\n");
      SCIPdialoghdlrClearBuffer(dialoghdlr);
   }
   else
   {
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter variable name: ", &varname, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }

      if( varname[0] != '\0' )
      {
         SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, varname, TRUE) );

         var = SCIPfindVar(scip, varname);
         if( var == NULL )
            SCIPdialogMessage(scip, NULL, "variable <%s> not found\n", varname);
         else
         {
            solval = SCIPgetSolVal(scip, sol, var);
            SCIPdialogMessage(scip, NULL, "%-32s", SCIPvarGetName(var));
            if( SCIPisInfinity(scip, solval) )
               SCIPdialogMessage(scip, NULL, " +infinity");
            else if( SCIPisInfinity(scip, -solval) )
               SCIPdialogMessage(scip, NULL, " -infinity");
            else
               SCIPdialogMessage(scip, NULL, " %20.15g", solval);
            SCIPdialogMessage(scip, NULL, " \t(obj:%.15g)\n", SCIPvarGetObj(var));
         }
      }
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display varbranchstatistics command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayVarbranchstatistics)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPprintBranchingStatistics(scip, NULL) );
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display LP solution quality command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayLPSolutionQuality)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPprintLPSolutionQuality(scip, NULL) );
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the help command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecHelp)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPdialogDisplayMenu(SCIPdialogGetParent(dialog), scip) );
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the display transsolution command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayTranssolution)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      if( SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
      {
         SCIPdialogMessage(scip, NULL, "best solution exists only in original problem space\n");
      }
      else
      {
         SCIP_CALL( SCIPprintBestTransSol(scip, NULL, FALSE) );
      }
   }
   else
      SCIPdialogMessage(scip, NULL, "no solution available\n");
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the free command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecFree)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIP_CALL( SCIPfreeProb(scip) );

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the newstart command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecNewstart)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIP_CALL( SCIPfreeSolve(scip, TRUE) );

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the transform command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecTransform)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPtransformProb(scip) );
      break;

   case SCIP_STAGE_TRANSFORMED:
      SCIPdialogMessage(scip, NULL, "problem is already transformed\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the concurrentopt command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecConcurrentOpt)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPsolveParallel(scip) );
      break;

   case SCIP_STAGE_SOLVED:
      SCIPdialogMessage(scip, NULL, "problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the optimize command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecOptimize)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPsolve(scip) );
      break;

   case SCIP_STAGE_SOLVED:
      SCIPdialogMessage(scip, NULL, "problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the presolve command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecPresolve)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      SCIP_CALL( SCIPpresolve(scip) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIPdialogMessage(scip, NULL, "problem is already presolved\n");
      break;

   case SCIP_STAGE_SOLVED:
      SCIPdialogMessage(scip, NULL, "problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the quit command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecQuit)
{  /*lint --e{715}*/
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = NULL;

   return SCIP_OKAY;
}

/** dialog execution method for the read command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecRead)
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;
   char* filename;
   SCIP_Bool endoffile;

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   if( filename[0] != '\0' )
   {
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      if( SCIPfileExists(filename) )
      {
         char* tmpfilename;
         char* extension;

         /* copy filename */
         SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpfilename, filename, (int)strlen(filename)+1) );
         extension = NULL;

         SCIPinfoMessage(scip, NULL, "\n");
         SCIPinfoMessage(scip, NULL, "read problem <%s>\n", filename);
         SCIPinfoMessage(scip, NULL, "============\n");
         SCIPinfoMessage(scip, NULL, "\n");

         do
         {
            retcode = SCIPreadProb(scip, tmpfilename, extension);
            if( retcode == SCIP_READERROR || retcode == SCIP_NOFILE )
            {
               if( extension == NULL )
                  SCIPdialogMessage(scip, NULL, "error reading file <%s>\n", tmpfilename);
               else
                  SCIPdialogMessage(scip, NULL, "error reading file <%s> using <%s> file format\n",
                     tmpfilename, extension);

               SCIP_CALL( SCIPfreeProb(scip) );
               break;
            }
            else if( retcode == SCIP_PLUGINNOTFOUND )
            {
               /* ask user once for a suitable reader */
               if( extension == NULL )
               {
                  SCIPdialogMessage(scip, NULL, "no reader for input file <%s> available\n", tmpfilename);

                  SCIPdialogMessage(scip, NULL, "The following readers are available for reading:\n");
                  displayReaders(scip, TRUE, FALSE);

                  SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
                        "select a suitable reader by extension (or return): ", &extension, &endoffile) );

                  if( extension[0] == '\0' )
                     break;
               }
               else
               {
                  SCIPdialogMessage(scip, NULL, "no reader for file extension <%s> available\n", extension);
                  extension = NULL;
               }
            }
            else
            {
               /* check if an unexpected error occurred during the reading process */
               SCIP_CALL( retcode );
               break;
            }
         }
         while( extension != NULL );

         /* free buffer array */
         SCIPfreeBufferArray(scip, &tmpfilename);
      }
      else
      {
         SCIPdialogMessage(scip, NULL, "file <%s> not found\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set default command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDefault)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIP_CALL( SCIPresetParams(scip) );
   SCIPdialogMessage(scip, NULL, "reset parameters to their default values\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set load command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetLoad)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   if( filename[0] != '\0' )
   {
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      if( SCIPfileExists(filename) )
      {
         SCIP_CALL( SCIPreadParams(scip, filename) );
         SCIPdialogMessage(scip, NULL, "loaded parameter file <%s>\n", filename);
      }
      else
      {
         SCIPdialogMessage(scip, NULL, "file <%s> not found\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set save command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSave)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   if( filename[0] != '\0' )
   {
      SCIP_RETCODE retcode;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      retcode =  SCIPwriteParams(scip, filename, TRUE, FALSE);

      if( retcode == SCIP_FILECREATEERROR )
      {
         SCIPdialogMessage(scip, NULL, "error creating file  <%s>\n", filename);
      }
      else
      {
         SCIP_CALL( retcode );
         SCIPdialogMessage(scip, NULL, "saved parameter file <%s>\n", filename);
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set diffsave command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDiffsave)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   if( filename[0] != '\0' )
   {
      SCIP_RETCODE retcode;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      retcode = SCIPwriteParams(scip, filename, TRUE, TRUE);

      if( retcode == SCIP_FILECREATEERROR )
      {
         SCIPdialogMessage(scip, NULL, "error creating file  <%s>\n", filename);
      }
      else
      {
         SCIP_CALL( retcode );
         SCIPdialogMessage(scip, NULL, "saved non-default parameter settings to file <%s>\n", filename);
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set parameter command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetParam)
{  /*lint --e{715}*/
   SCIP_PARAM* param;
   char prompt[SCIP_MAXSTRLEN];
   char* valuestr;
   SCIP_Bool boolval;
   int intval;
   SCIP_Longint longintval;
   SCIP_Real realval;
   char charval;
   SCIP_Bool endoffile;
   SCIP_Bool error;

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* get the parameter to set */
   param = (SCIP_PARAM*)SCIPdialogGetData(dialog);

   /* depending on the parameter type, request a user input */
   switch( SCIPparamGetType(param) )
   {
   case SCIP_PARAMTYPE_BOOL:
      (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: %s, new value (TRUE/FALSE): ",
         SCIPparamGetBool(param) ? "TRUE" : "FALSE");
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      boolval = parseBoolValue(scip, valuestr, &error);

      if( error )
      {
         SCIPdialogMessage(scip, NULL, "\nInvalid value <%s> for bool parameter <%s>. Must be <0>, <1>, <FALSE>, or <TRUE>.\n\n",
            valuestr, SCIPparamGetName(param));
         SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr, TRUE) );
      }
      else
      {
         assert(SCIPisBoolParamValid(scip, param, boolval));

         SCIP_CALL( SCIPchgBoolParam(scip, param, boolval) );
         SCIPdialogMessage(scip, NULL, "%s = %s\n", SCIPparamGetName(param), boolval ? "TRUE" : "FALSE");
         SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, boolval ? "TRUE" : "FALSE", TRUE) );

      }

      break;

   case SCIP_PARAMTYPE_INT:
      (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: %d, new value [%d,%d]: ",
         SCIPparamGetInt(param), SCIPparamGetIntMin(param), SCIPparamGetIntMax(param));
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr, TRUE) );

      if( sscanf(valuestr, "%d", &intval) != 1 || !SCIPisIntParamValid(scip, param, intval) )
      {
         SCIPdialogMessage(scip, NULL, "\nInvalid value <%s> for int parameter <%s>. Must be integral in range [%d,%d].\n\n",
            valuestr, SCIPparamGetName(param), SCIPparamGetIntMin(param), SCIPparamGetIntMax(param));
      }
      else
      {
         SCIP_CALL( SCIPchgIntParam(scip, param, intval) );
         SCIPdialogMessage(scip, NULL, "%s = %d\n", SCIPparamGetName(param), intval);
      }

      break;

   case SCIP_PARAMTYPE_LONGINT:
      (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: %" SCIP_LONGINT_FORMAT ", new value [%" SCIP_LONGINT_FORMAT ",%" SCIP_LONGINT_FORMAT "]: ",
         SCIPparamGetLongint(param), SCIPparamGetLongintMin(param), SCIPparamGetLongintMax(param));
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr, TRUE) );

      if( sscanf(valuestr, "%" SCIP_LONGINT_FORMAT, &longintval) != 1 || !SCIPisLongintParamValid(scip, param, longintval) )
      {
         SCIPdialogMessage(scip, NULL, "\nInvalid value <%s> for longint parameter <%s>. Must be integral in range [%" SCIP_LONGINT_FORMAT ",%" SCIP_LONGINT_FORMAT "].\n\n",
            valuestr, SCIPparamGetName(param), SCIPparamGetLongintMin(param), SCIPparamGetLongintMax(param));
      }
      else
      {
         SCIP_CALL( SCIPchgLongintParam(scip, param, longintval) );
         SCIPdialogMessage(scip, NULL, "%s = %" SCIP_LONGINT_FORMAT "\n", SCIPparamGetName(param), longintval);
      }
      break;

   case SCIP_PARAMTYPE_REAL:
      (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: %.15g, new value [%.15g,%.15g]: ",
         SCIPparamGetReal(param), SCIPparamGetRealMin(param), SCIPparamGetRealMax(param));
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr, TRUE) );

      if( sscanf(valuestr, "%" SCIP_REAL_FORMAT, &realval) != 1 || !SCIPisRealParamValid(scip, param, realval) )
      {
         SCIPdialogMessage(scip, NULL, "\nInvalid real parameter value <%s> for parameter <%s>. Must be in range [%.15g,%.15g].\n\n",
            valuestr, SCIPparamGetName(param), SCIPparamGetRealMin(param), SCIPparamGetRealMax(param));
      }
      else
      {
         SCIP_CALL( SCIPchgRealParam(scip, param, realval) );
         SCIPdialogMessage(scip, NULL, "%s = %.15g\n", SCIPparamGetName(param), realval);
      }
      break;

   case SCIP_PARAMTYPE_CHAR:
      (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: <%c>, new value: ", SCIPparamGetChar(param));
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr, TRUE) );

      if( sscanf(valuestr, "%c", &charval) != 1 || !SCIPisCharParamValid(scip, param, charval) )
      {
         SCIPdialogMessage(scip, NULL, "\nInvalid char parameter value <%s>. Must be in set {%s}.\n\n",
            valuestr, SCIPparamGetCharAllowedValues(param));
      }
      else
      {
         SCIP_CALL( SCIPchgCharParam(scip, param, charval) );
         SCIPdialogMessage(scip, NULL, "%s = %c\n", SCIPparamGetName(param), charval);
      }
      break;

   case SCIP_PARAMTYPE_STRING:
      (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: <%s>, new value: ", SCIPparamGetString(param));
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr, TRUE) );

      if ( !SCIPisStringParamValid(scip, param, valuestr) )
      {
         SCIPdialogMessage(scip, NULL, "\nInvalid character in string parameter.\n\n");
      }
      else
      {
         SCIP_CALL( SCIPchgStringParam(scip, param, valuestr) );
         SCIPdialogMessage(scip, NULL, "%s = %s\n", SCIPparamGetName(param), valuestr);
      }
      break;

   default:
      SCIPerrorMessage("invalid parameter type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** dialog description method for the set parameter command */
SCIP_DECL_DIALOGDESC(SCIPdialogDescSetParam)
{  /*lint --e{715}*/
   SCIP_PARAM* param;
   char valuestr[SCIP_MAXSTRLEN];

   /* get the parameter to set */
   param = (SCIP_PARAM*)SCIPdialogGetData(dialog);

   /* retrieve parameter's current value */
   switch( SCIPparamGetType(param) )
   {
   case SCIP_PARAMTYPE_BOOL:
      if( SCIPparamGetBool(param) )
         (void) SCIPsnprintf(valuestr, SCIP_MAXSTRLEN, "TRUE");
      else
         (void) SCIPsnprintf(valuestr, SCIP_MAXSTRLEN, "FALSE");
      break;

   case SCIP_PARAMTYPE_INT:
      (void) SCIPsnprintf(valuestr, SCIP_MAXSTRLEN, "%d", SCIPparamGetInt(param));
      break;

   case SCIP_PARAMTYPE_LONGINT:
      (void) SCIPsnprintf(valuestr, SCIP_MAXSTRLEN, "%" SCIP_LONGINT_FORMAT, SCIPparamGetLongint(param));
      break;

   case SCIP_PARAMTYPE_REAL:
      (void) SCIPsnprintf(valuestr, SCIP_MAXSTRLEN, "%.15g", SCIPparamGetReal(param));
      if( strchr(valuestr, '.') == NULL && strchr(valuestr, 'e') == NULL )
         (void) SCIPsnprintf(valuestr, SCIP_MAXSTRLEN, "%.1f", SCIPparamGetReal(param));
      break;

   case SCIP_PARAMTYPE_CHAR:
      (void) SCIPsnprintf(valuestr, SCIP_MAXSTRLEN, "%c", SCIPparamGetChar(param));
      break;

   case SCIP_PARAMTYPE_STRING:
      (void) SCIPsnprintf(valuestr, SCIP_MAXSTRLEN, "%s", SCIPparamGetString(param));
      break;

   default:
      SCIPerrorMessage("invalid parameter type\n");
      return SCIP_INVALIDDATA;
   }
   valuestr[SCIP_MAXSTRLEN-1] = '\0';

   /* display parameter's description */
   SCIPdialogMessage(scip, NULL, "%s", SCIPparamGetDesc(param));

   /* display parameter's current value */
   SCIPdialogMessage(scip, NULL, " [%s]", valuestr);

   return SCIP_OKAY;
}

/** dialog execution method for the fix parameter command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecFixParam)
{  /*lint --e{715}*/
   SCIP_PARAM* param;
   char prompt[SCIP_MAXSTRLEN];
   char* valuestr;
   SCIP_Bool fix;
   SCIP_Bool endoffile;
   SCIP_Bool error;

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* get the parameter to fix */
   param = (SCIP_PARAM*)SCIPdialogGetData(dialog);

   (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current fixing status: %s, new value (TRUE/FALSE): ",
         SCIPparamIsFixed(param) ? "TRUE" : "FALSE");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   fix = parseBoolValue(scip, valuestr, &error);

   if( !error )
   {
      SCIPparamSetFixed(param, fix);
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, (fix ? "TRUE" : "FALSE"), TRUE) );
      SCIPdialogMessage(scip, NULL, "<%s> %s\n", SCIPparamGetName(param), (fix ? "fixed" : "unfixed"));
   }
   else
   {
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr, TRUE) );
      SCIPdialogMessage(scip, NULL, "\nInvalid value <%s> for fixing status. Must be <0>, <1>, <FALSE>, or <TRUE>.\n\n",
         valuestr);
   }

   return SCIP_OKAY;
}

/** dialog description method for the fix parameter command */
SCIP_DECL_DIALOGDESC(SCIPdialogDescFixParam)
{  /*lint --e{715}*/
   SCIP_PARAM* param;

   /* get the parameter to set */
   param = (SCIP_PARAM*)SCIPdialogGetData(dialog);

   /* display parameter's description */
   SCIPdialogMessage(scip, NULL, "%s", SCIPparamGetDesc(param));

   /* display parameter's current fixing status */
   if( SCIPparamIsFixed(param) )
      SCIPdialogMessage(scip, NULL, " [fixed]");
   else
      SCIPdialogMessage(scip, NULL, " [not fixed]");

   return SCIP_OKAY;
}

/** dialog execution method for the set branching direction command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetBranchingDirection)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   char prompt[SCIP_MAXSTRLEN];
   char* valuestr;
   int direction;
   SCIP_Bool endoffile;

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* branching priorities cannot be set, if no problem was created */
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPdialogMessage(scip, NULL, "cannot set branching directions before problem was created\n");
      return SCIP_OKAY;
   }

   /* get variable name from user */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "variable name: ", &valuestr, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   /* find variable */
   var = SCIPfindVar(scip, valuestr);
   if( var == NULL )
   {
      SCIPdialogMessage(scip, NULL, "variable <%s> does not exist in problem\n", valuestr);
      return SCIP_OKAY;
   }

   /* get new branching direction from user */
   switch( SCIPvarGetBranchDirection(var) )
   {
   case SCIP_BRANCHDIR_DOWNWARDS:
      direction = -1;
      break;
   case SCIP_BRANCHDIR_AUTO:
      direction = 0;
      break;
   case SCIP_BRANCHDIR_UPWARDS:
      direction = +1;
      break;
   case SCIP_BRANCHDIR_FIXED:
   default:
      SCIPerrorMessage("invalid preferred branching direction <%d> of variable <%s>\n",
         SCIPvarGetBranchDirection(var), SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }
   (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: %d, new value: ", direction);
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   SCIPescapeString(prompt, SCIP_MAXSTRLEN, SCIPvarGetName(var));
   (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "%s %s", prompt, valuestr);
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, prompt, FALSE) );

   if( sscanf(valuestr, "%d", &direction) != 1 )
   {
      SCIPdialogMessage(scip, NULL, "\ninvalid input <%s>\n\n", valuestr);
      return SCIP_OKAY;
   }
   if( direction < -1 || direction > +1 )
   {
      SCIPdialogMessage(scip, NULL, "\ninvalid input <%d>: direction must be -1, 0, or +1\n\n", direction);
      return SCIP_OKAY;
   }

   /* set new branching direction */
   if( direction == -1 )
      SCIP_CALL( SCIPchgVarBranchDirection(scip, var, SCIP_BRANCHDIR_DOWNWARDS) );
   else if( direction == 0 )
      SCIP_CALL( SCIPchgVarBranchDirection(scip, var, SCIP_BRANCHDIR_AUTO) );
   else
      SCIP_CALL( SCIPchgVarBranchDirection(scip, var, SCIP_BRANCHDIR_UPWARDS) );

   SCIPdialogMessage(scip, NULL, "branching direction of variable <%s> set to %d\n", SCIPvarGetName(var), direction);

   return SCIP_OKAY;
}

/** dialog execution method for the set branching priority command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetBranchingPriority)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   char prompt[SCIP_MAXSTRLEN];
   char* valuestr;
   int priority;
   SCIP_Bool endoffile;

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* branching priorities cannot be set, if no problem was created */
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPdialogMessage(scip, NULL, "cannot set branching priorities before problem was created\n");
      return SCIP_OKAY;
   }

   /* get variable name from user */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "variable name: ", &valuestr, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   /* find variable */
   var = SCIPfindVar(scip, valuestr);
   if( var == NULL )
   {
      SCIPdialogMessage(scip, NULL, "variable <%s> does not exist in problem\n", valuestr);
      return SCIP_OKAY;
   }

   /* get new branching priority from user */
   (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: %d, new value: ", SCIPvarGetBranchPriority(var));
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   SCIPescapeString(prompt, SCIP_MAXSTRLEN, SCIPvarGetName(var));
   (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "%s %s", prompt, valuestr);
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, prompt, FALSE) );

   if( sscanf(valuestr, "%d", &priority) != 1 )
   {
      SCIPdialogMessage(scip, NULL, "\ninvalid input <%s>\n\n", valuestr);
      return SCIP_OKAY;
   }

   /* set new branching priority */
   SCIP_CALL( SCIPchgVarBranchPriority(scip, var, priority) );
   SCIPdialogMessage(scip, NULL, "branching priority of variable <%s> set to %d\n", SCIPvarGetName(var), SCIPvarGetBranchPriority(var));

   return SCIP_OKAY;
}

/** dialog execution method for the set heuristics aggressive command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsAggressive)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set heuristics default command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsDefault)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_DEFAULT, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set heuristics fast command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsFast)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set heuristics off command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsOff)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set presolving aggressive command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingAggressive)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_AGGRESSIVE, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set presolving default command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingDefault)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_DEFAULT, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set presolving fast command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingFast)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set presolving off command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingOff)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set separating aggressive command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingAggressive)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_AGGRESSIVE, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set separating default command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingDefault)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_DEFAULT, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set separating fast command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingFast)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_FAST, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set separating off command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingOff)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set emphasis counter command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisCounter)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* reset SCIP parameters */
   SCIP_CALL( SCIPresetParams(scip) );

   /* set parameters for counting problems */
   SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_COUNTER, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set emphasis cpsolver command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisCpsolver)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* reset SCIP parameters */
   SCIP_CALL( SCIPresetParams(scip) );

   /* set parameters for CP like search problems */
   SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_CPSOLVER, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set emphasis easy CIP command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisEasycip)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* reset SCIP parameters */
   SCIP_CALL( SCIPresetParams(scip) );

   /* set parameters for easy CIP problems */
   SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_EASYCIP, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set emphasis feasibility command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisFeasibility)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* reset SCIP parameters */
   SCIP_CALL( SCIPresetParams(scip) );

   /* set parameters for feasibility problems */
   SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_FEASIBILITY, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set emphasis hard LP command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisHardlp)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* reset SCIP parameters */
   SCIP_CALL( SCIPresetParams(scip) );

   /* set parameters for problems with hard LP */
   SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_HARDLP, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set emphasis optimality command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisOptimality)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* reset SCIP parameters */
   SCIP_CALL( SCIPresetParams(scip) );

   /* set parameters for problems to prove optimality fast */
   SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_OPTIMALITY, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set limits objective command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetLimitsObjective)
{  /*lint --e{715}*/
   char prompt[SCIP_MAXSTRLEN];
   char* valuestr;
   SCIP_Real objlim;
   SCIP_Bool endoffile;

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* objective limit cannot be set, if no problem was created */
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPdialogMessage(scip, NULL, "cannot set objective limit before problem was created\n");
      return SCIP_OKAY;
   }

   /* get new objective limit from user */
   (void) SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "current value: %.15g, new value: ", SCIPgetObjlimit(scip));
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &valuestr, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr, TRUE) );

   if( sscanf(valuestr, "%" SCIP_REAL_FORMAT, &objlim) != 1 )
   {
      SCIPdialogMessage(scip, NULL, "\ninvalid input <%s>\n\n", valuestr);
      return SCIP_OKAY;
   }

   /* check, if new objective limit is valid */
   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM
      && SCIPtransformObj(scip, objlim) > SCIPtransformObj(scip, SCIPgetObjlimit(scip)) )
   {
      SCIPdialogMessage(scip, NULL, "\ncannot relax objective limit from %.15g to %.15g after problem was transformed\n\n",
         SCIPgetObjlimit(scip), objlim);
      return SCIP_OKAY;
   }

   /* set new objective limit */
   SCIP_CALL( SCIPsetObjlimit(scip, objlim) );
   SCIPdialogMessage(scip, NULL, "objective value limit set to %.15g\n", SCIPgetObjlimit(scip));

   return SCIP_OKAY;
}

/** dialog execution method for the write LP command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteLp)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   /* node relaxations only exist in solving & solved stage */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
   {
      SCIPdialogMessage(scip, NULL, "There is no node LP relaxation before solving starts\n");
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
   {
      SCIPdialogMessage(scip, NULL, "There is no node LP relaxation after problem was solved\n");
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      SCIP_RETCODE retcode;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );
      retcode =  SCIPwriteLP(scip, filename);

      if( retcode == SCIP_FILECREATEERROR )
      {
         SCIPdialogMessage(scip, NULL, "error not creating file  <%s>\n", filename);
      }
      else
      {
         SCIP_CALL( retcode );

         SCIPdialogMessage(scip, NULL, "written node LP relaxation to file <%s>\n", filename);
      }
   }

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write MIP command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteMip)
{  /*lint --e{715}*/
   char command[SCIP_MAXSTRLEN];
   char filename[SCIP_MAXSTRLEN];
   SCIP_Bool endoffile;
   char* valuestr;
   SCIP_Bool offset;
   SCIP_Bool generic;
   SCIP_Bool lazyconss;
   SCIP_Bool error;

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* node relaxations only exist in solving & solved stage */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
   {
      SCIPdialogMessage(scip, NULL, "There is no node MIP relaxation before solving starts\n");
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
   {
      SCIPdialogMessage(scip, NULL, "There is no node MIP relaxation after problem was solved\n");
      return SCIP_OKAY;
   }

   /* first get file name */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &valuestr, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   (void)strncpy(filename, valuestr, SCIP_MAXSTRLEN-1);

   /* second ask for generic variable and row names */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "using generic variable and row names (TRUE/FALSE): ",
         &valuestr, &endoffile) );

   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   generic = parseBoolValue(scip, valuestr, &error);

   if( error )
   {
      SCIPdialogMessage(scip, NULL, "\nInvalid value <%s>. Must be <0>, <1>, <FALSE>, or <TRUE>.\n\n",
         valuestr);

      return SCIP_OKAY;
   }

   /* adjust command and add to the history */
   SCIPescapeString(command, SCIP_MAXSTRLEN, filename);
   (void) SCIPsnprintf(command, SCIP_MAXSTRLEN, "%s %s", command, generic ? "TRUE" : "FALSE");

   /* third ask if for adjusting the objective offset */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "using original objective function (TRUE/FALSE): ",
         &valuestr, &endoffile) );

   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   offset = parseBoolValue(scip, valuestr, &error);

   if( error )
   {
      SCIPdialogMessage(scip, NULL, "\nInvalid value <%s>. Must be <0>, <1>, <FALSE>, or <TRUE>.\n\n",
         valuestr);

      return SCIP_OKAY;
   }

   (void) SCIPsnprintf(command, SCIP_MAXSTRLEN, "%s %s", command, offset ? "TRUE" : "FALSE");
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, command, FALSE) );

   /* fourth ask for lazy constraints */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "output removable rows as lazy constraints (TRUE/FALSE): ",
         &valuestr, &endoffile) );

   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   lazyconss = parseBoolValue(scip, valuestr, &error);

   if( error )
   {
      SCIPdialogMessage(scip, NULL, "\nInvalid value <%s>. Must be <0>, <1>, <FALSE>, or <TRUE>.\n\n",
         valuestr);

      return SCIP_OKAY;
   }

   /* adjust command and add to the history */
   SCIPescapeString(command, SCIP_MAXSTRLEN, filename);
   (void) SCIPsnprintf(command, SCIP_MAXSTRLEN, "%s %s", command, lazyconss ? "TRUE" : "FALSE");

   /* execute command */
   SCIP_CALL( SCIPwriteMIP(scip, filename, generic, offset, lazyconss) );
   SCIPdialogMessage(scip, NULL, "written node MIP relaxation to file <%s>\n", filename);

   SCIPdialogMessage(scip, NULL, "\n");

   return SCIP_OKAY;
}


/** dialog execution method for the write NLP command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteNlp)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   /* node relaxations only exist in solving & solved stage */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
   {
      SCIPdialogMessage(scip, NULL, "There is no node NLP relaxation before solving starts\n");
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
   {
      SCIPdialogMessage(scip, NULL, "There is no node NLP relaxation after problem was solved\n");
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      return SCIP_OKAY;
   }
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdialogMessage(scip, NULL, "There has been no node NLP relaxation constructed\n");
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      SCIP_RETCODE retcode;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );
      retcode =  SCIPwriteNLP(scip, filename);

      if( retcode == SCIP_FILECREATEERROR )
      {
         SCIPdialogMessage(scip, NULL, "error not creating file  <%s>\n", filename);
      }
      else
      {
         SCIP_CALL( retcode );

         SCIPdialogMessage(scip, NULL, "written node NLP relaxation to file <%s>\n", filename);
      }
   }

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write problem command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteProblem)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM )
   {
      SCIP_CALL( writeProblem(scip, dialog, dialoghdlr, nextdialog, FALSE, FALSE) );
   }
   else
      SCIPdialogMessage(scip, NULL, "no problem available\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write generic problem command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteGenProblem)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM )
   {
      SCIP_CALL( writeProblem(scip, dialog, dialoghdlr, nextdialog, FALSE, TRUE) );
   }
   else
      SCIPdialogMessage(scip, NULL, "no problem available\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write solution command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteSolution)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      FILE* file;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPdialogMessage(scip, NULL, "error creating file <%s>\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
      else
      {
         SCIP_Bool printzeros;

         SCIPinfoMessage(scip, file, "solution status: ");
         SCIP_CALL_FINALLY( SCIPprintStatus(scip, file), fclose(file) );

         SCIP_CALL_FINALLY( SCIPgetBoolParam(scip, "write/printzeros", &printzeros), fclose(file) );

         SCIPinfoMessage(scip, file, "\n");
         SCIP_CALL_FINALLY( SCIPprintBestSol(scip, file, printzeros), fclose(file) );

         SCIPdialogMessage(scip, NULL, "written solution information to file <%s>\n", filename);
         fclose(file);
      }
   }

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write mipstart command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteMIPStart)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      FILE* file;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPdialogMessage(scip, NULL, "error creating file <%s>\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
      else
      {
         SCIP_SOL* sol;

         SCIPinfoMessage(scip, file, "\n");

         sol = SCIPgetBestSol(scip);

         if( sol == NULL )
         {
            SCIPdialogMessage(scip, NULL, "no mip start available\n");
            fclose(file);
         }
         else
         {
            SCIP_CALL_FINALLY( SCIPprintMIPStart(scip, sol, file), fclose(file) );

            SCIPdialogMessage(scip, NULL, "written mip start information to file <%s>\n", filename);
            fclose(file);
         }
      }
   }

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for writing command line history */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteCommandHistory)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      SCIP_RETCODE retcode;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      retcode = SCIPdialogWriteHistory(filename);

      if( retcode != SCIP_OKAY )
      {
         SCIPdialogMessage(scip, NULL, "error writing to file <%s>\n"
               "check that the directory exists and that you have correct permissions\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
      else
      {
         SCIPdialogMessage(scip, NULL, "wrote available command line history to <%s>\n", filename);
      }
   }

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write finitesolution command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteFiniteSolution)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      FILE* file;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPdialogMessage(scip, NULL, "error creating file <%s>\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
      else
      {
         SCIP_SOL* bestsol = SCIPgetBestSol(scip);
         SCIP_Bool printzeros;

         SCIPinfoMessage(scip, file, "solution status: ");

         SCIP_CALL_FINALLY( SCIPprintStatus(scip, file), fclose(file) );

         SCIPinfoMessage(scip, file, "\n");

         if( bestsol != NULL )
         {
            SCIP_SOL* sol;
            SCIP_Bool success;

            SCIP_CALL_FINALLY( SCIPcreateFiniteSolCopy(scip, &sol, bestsol, &success), fclose(file) );

            SCIP_CALL_FINALLY( SCIPgetBoolParam(scip, "write/printzeros", &printzeros), fclose(file) );

            if( sol != NULL )
            {
               SCIP_CALL_FINALLY( SCIPprintSol(scip, sol, file, printzeros), fclose(file) );

               SCIPdialogMessage(scip, NULL, "written solution information to file <%s>\n", filename);

               SCIP_CALL_FINALLY( SCIPfreeSol(scip, &sol), fclose(file) );
            }
            else
            {
               SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "finite solution could not be created\n");
               SCIPdialogMessage(scip, NULL, "finite solution could not be created\n", filename);
            }
         }
         else
         {
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "no solution available\n");
            SCIPdialogMessage(scip, NULL, "no solution available\n", filename);
         }

         fclose(file);
      }
   }

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write statistics command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteStatistics)
{  /*lint --e{715}*/
   char* filename;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      FILE* file;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPdialogMessage(scip, NULL, "error creating file <%s>\n", filename);
         SCIPprintSysError(filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
      else
      {
         SCIP_CALL_FINALLY( SCIPprintStatistics(scip, file), fclose(file) );

         SCIPdialogMessage(scip, NULL, "written statistics to file <%s>\n", filename);
         fclose(file);
      }
   }

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write transproblem command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteTransproblem)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( writeProblem(scip, dialog, dialoghdlr, nextdialog, TRUE, FALSE) );
   }
   else
      SCIPdialogMessage(scip, NULL, "no transformed problem available\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the write generic transproblem command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteGenTransproblem)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( writeProblem(scip, dialog, dialoghdlr, nextdialog, TRUE, TRUE) );
   }
   else
      SCIPdialogMessage(scip, NULL, "no transformed problem available\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for solution validation */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecValidateSolve)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPdialogMessage(scip, NULL, "\nNo problem available for validation\n");
   }
   else
   {
      char *refstrs[2];
      SCIP_Real refvals[2] = {SCIP_INVALID, SCIP_INVALID};
      const char* primaldual[] = {"primal", "dual"};
      char prompt[SCIP_MAXSTRLEN];
      int i;

      /* read in primal and dual reference values */
      for( i = 0; i < 2; ++i )
      {
         char * endptr;
         SCIP_Bool endoffile;

         (void)SCIPsnprintf(prompt, SCIP_MAXSTRLEN, "Please enter %s validation reference bound (or use +/-infinity) :", primaldual[i]);
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt, &(refstrs[i]), &endoffile) );

         /* treat no input as SCIP_UNKNOWN */
         if( endoffile || strncmp(refstrs[i], "\0", 1) == 0 ) /*lint !e840*/
         {
            refvals[i] = SCIP_UNKNOWN;
         }
         else if( strncmp(refstrs[i], "q", 1) == 0 )
            break;
         else if( ! SCIPparseReal(scip, refstrs[i], &refvals[i], &endptr) )
         {
            SCIPdialogMessage(scip, NULL, "Could not parse value '%s', please try again or type 'q' to quit\n", refstrs[i]);
            --i;
         }
      }

      /* check if the loop finished by checking the value of 'i'. Do not validate if user input is missing */
      if( i == 2 ) /*lint !e850*/
      {
         assert(refvals[0] != SCIP_INVALID); /*lint !e777*/
         assert(refvals[1] != SCIP_INVALID); /*lint !e777*/
         SCIP_CALL( SCIPvalidateSolve(scip, refvals[0], refvals[1], SCIPfeastol(scip), FALSE, NULL, NULL, NULL) );
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for linear constraint type classification */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayLinearConsClassification)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
      SCIPdialogMessage(scip, NULL, "\nNo problem available for classification\n");
   else
   {
      SCIP_LINCONSSTATS* linconsstats;

      SCIP_CALL( SCIPlinConsStatsCreate(scip, &linconsstats) );

      /* call linear constraint classification and print the statistics to standard out */
      SCIP_CALL( SCIPclassifyConstraintTypesLinear(scip, linconsstats) );

      SCIPprintLinConsStats(scip, NULL, linconsstats);

      SCIPlinConsStatsFree(scip, &linconsstats);
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** creates a root dialog */
SCIP_RETCODE SCIPcreateRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         root                /**< pointer to store the root dialog */
   )
{
   SCIP_CALL( SCIPincludeDialog(scip, root, 
         dialogCopyDefault,
         SCIPdialogExecMenuLazy, NULL, NULL,
         "SCIP", "SCIP's main menu", TRUE, NULL) );

   SCIP_CALL( SCIPsetRootDialog(scip, *root) );
   SCIP_CALL( SCIPreleaseDialog(scip, root) );
   *root = SCIPgetRootDialog(scip);

   return SCIP_OKAY;
}


/** includes or updates the default dialog menus in SCIP */
SCIP_RETCODE SCIPincludeDialogDefault(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DIALOG* root;
   SCIP_DIALOG* submenu;
   SCIP_DIALOG* dialog;

   /* root menu */
   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      SCIP_CALL( SCIPcreateRootDialog(scip, &root) );
   }

   /* change */
   if( !SCIPdialogHasEntry(root, "change") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu, 
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "change", "change the problem", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(root, "change", &submenu) != 1 )
   {
      SCIPerrorMessage("change sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* change add */
   if( !SCIPdialogHasEntry(submenu, "add") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecChangeAddCons, NULL, NULL,
            "add", "add constraint", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* change bounds */
   if( !SCIPdialogHasEntry(submenu, "bounds") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecChangeBounds, NULL, NULL,
            "bounds", "change bounds of a variable", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* free transformed problem */
   if( !SCIPdialogHasEntry(submenu, "freetransproblem") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecChangeFreetransproblem, NULL, NULL,
            "freetransproblem", "free transformed problem", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* change objective sense */
   if( !SCIPdialogHasEntry(submenu, "objsense") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecChangeObjSense, NULL, NULL,
            "objsense", "change objective sense", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* checksol */
   if( !SCIPdialogHasEntry(root, "checksol") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, 
            NULL,
            SCIPdialogExecChecksol, NULL, NULL,
            "checksol", "double checks best solution w.r.t. original problem", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display */
   if( !SCIPdialogHasEntry(root, "display") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "display", "display information", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(root, "display", &submenu) != 1 )
   {
      SCIPerrorMessage("display sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* display branching */
   if( !SCIPdialogHasEntry(submenu, "branching") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayBranching, NULL, NULL,
            "branching", "display branching rules", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display compressions */
   if( !SCIPdialogHasEntry(submenu, "compression") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayCompression, NULL, NULL,
            "compression", "display compression techniques", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display conflict */
   if( !SCIPdialogHasEntry(submenu, "conflict") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayConflict, NULL, NULL,
            "conflict", "display conflict handlers", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display conshdlrs */
   if( !SCIPdialogHasEntry(submenu, "conshdlrs") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayConshdlrs, NULL, NULL,
            "conshdlrs", "display constraint handlers", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display displaycols */
   if( !SCIPdialogHasEntry(submenu, "displaycols") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayDisplaycols, NULL, NULL,
            "displaycols", "display display columns", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display heuristics */
   if( !SCIPdialogHasEntry(submenu, "heuristics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayHeuristics, NULL, NULL,
            "heuristics", "display primal heuristics", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display memory */
   if( !SCIPdialogHasEntry(submenu, "memory") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayMemory, NULL, NULL,
            "memory", "display memory diagnostics", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display nlpi */
   if( !SCIPdialogHasEntry(submenu, "nlpis") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayNlpi, NULL, NULL,
            "nlpis", "display NLP solver interfaces", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display nodeselectors */
   if( !SCIPdialogHasEntry(submenu, "nodeselectors") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayNodeselectors, NULL, NULL,
            "nodeselectors", "display node selectors", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display parameters */
   if( !SCIPdialogHasEntry(submenu, "parameters") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayParameters, NULL, NULL,
            "parameters", "display non-default parameter settings", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display presolvers */
   if( !SCIPdialogHasEntry(submenu, "presolvers") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayPresolvers, NULL, NULL,
            "presolvers", "display presolvers", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display pricers */
   if( !SCIPdialogHasEntry(submenu, "pricers") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayPricers, NULL, NULL,
            "pricers", "display pricers", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display problem */
   if( !SCIPdialogHasEntry(submenu, "problem") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayProblem, NULL, NULL,
            "problem", "display original problem", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display propagators */
   if( !SCIPdialogHasEntry(submenu, "propagators") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayPropagators, NULL, NULL,
            "propagators", "display propagators", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display readers */
   if( !SCIPdialogHasEntry(submenu, "readers") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayReaders, NULL, NULL,
            "readers", "display file readers", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display relaxing */
   if( !SCIPdialogHasEntry(submenu, "relaxators") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayRelaxators, NULL, NULL,
            "relaxators", "display relaxators", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display separators */
   if( !SCIPdialogHasEntry(submenu, "separators") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplaySeparators, NULL, NULL,
            "separators", "display cut separators", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display solution */
   if( !SCIPdialogHasEntry(submenu, "solution") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplaySolution, NULL, NULL,
            "solution", "display best primal solution", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display finite solution */
   if( !SCIPdialogHasEntry(submenu, "finitesolution") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayFiniteSolution, NULL, NULL,
            "finitesolution", "display best primal solution (try to make solution values finite, first)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display solution */
   if( !SCIPdialogHasEntry(submenu, "dualsolution") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
                                   NULL,
                                   SCIPdialogExecDisplayDualSolution, NULL, NULL,
                                   "dualsolution", "display dual solution vector (LP only, without presolving)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display solution */
   if( !SCIPdialogHasEntry(submenu, "sols") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplaySolutionPool, NULL, NULL,
            "sols", "display solutions from pool", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display statistics */
   if( !SCIPdialogHasEntry(submenu, "statistics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayStatistics, NULL, NULL,
            "statistics", "display problem and optimization statistics", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display reoptimization statistics */
   if( !SCIPdialogHasEntry(submenu, "reoptstatistics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayReoptStatistics, NULL, NULL,
            "reoptstatistics", "display reoptimitazion statistics", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display transproblem */
   if( !SCIPdialogHasEntry(submenu, "transproblem") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayTransproblem, NULL, NULL,
            "transproblem", "display current node transformed problem", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display value */
   if( !SCIPdialogHasEntry(submenu, "value") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayValue, NULL, NULL,
            "value", "display value of single variable in best primal solution", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display varbranchstatistics */
   if( !SCIPdialogHasEntry(submenu, "varbranchstatistics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayVarbranchstatistics, NULL, NULL,
            "varbranchstatistics", "display statistics for branching on variables", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display varbranchstatistics */
   if( !SCIPdialogHasEntry(submenu, "lpsolquality") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayLPSolutionQuality, NULL, NULL,
            "lpsolquality", "display quality of the current LP solution, if available", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display transsolution */
   if( !SCIPdialogHasEntry(submenu, "transsolution") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayTranssolution, NULL, NULL,
            "transsolution", "display best primal solution in transformed variables", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display linear constraint type classification */
   if( !SCIPdialogHasEntry(submenu, "linclass") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecDisplayLinearConsClassification, NULL, NULL,
            "linclass", "linear constraint classification as used for MIPLIB", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* free */
   if( !SCIPdialogHasEntry(root, "free") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecFree, NULL, NULL,
            "free", "free current problem from memory", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* help */
   if( !SCIPdialogHasEntry(root, "help") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecHelp, NULL, NULL,
            "help", "display this help", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* newstart */
   if( !SCIPdialogHasEntry(root, "newstart") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecNewstart, NULL, NULL,
            "newstart", "reset branch and bound tree to start again from root", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

#ifndef NDEBUG
   /* transform problem (for debugging) */
   if( !SCIPdialogHasEntry(root, "transform") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecTransform, NULL, NULL,
            "transform", "transforms problem from original state", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }
#endif

   /* optimize */
   if( !SCIPdialogHasEntry(root, "optimize") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecOptimize, NULL, NULL,
            "optimize", "solve the problem", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* optimize */
   if( !SCIPdialogHasEntry(root, "concurrentopt") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
                                   NULL,
                                   SCIPdialogExecConcurrentOpt, NULL, NULL,
                                   "concurrentopt", "solve the problem using concurrent solvers", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* presolve */
   if( !SCIPdialogHasEntry(root, "presolve") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecPresolve, NULL, NULL,
            "presolve", "solve the problem, but stop after presolving stage", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* quit */
   if( !SCIPdialogHasEntry(root, "quit") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecQuit, NULL, NULL,
            "quit", "leave SCIP", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* read */
   if( !SCIPdialogHasEntry(root, "read") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecRead, NULL, NULL,
            "read", "read a problem", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set */
   SCIP_CALL( SCIPincludeDialogDefaultSet(scip) );

   /* fix */
   SCIP_CALL( SCIPincludeDialogDefaultFix(scip) );

   /* write */
   if( !SCIPdialogHasEntry(root, "write") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "write", "write information to file", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(root, "write", &submenu) != 1 )
   {
      SCIPerrorMessage("write sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* write LP */
   if( !SCIPdialogHasEntry(submenu, "lp") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteLp, NULL, NULL,
            "lp", "write current node LP relaxation in LP format to file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write MIP */
   if( !SCIPdialogHasEntry(submenu, "mip") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteMip, NULL, NULL,
            "mip", "write current node MIP relaxation in LP format to file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write NLP */
   if( !SCIPdialogHasEntry(submenu, "nlp") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteNlp, NULL, NULL,
            "nlp", "write current node NLP relaxation to file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write problem */
   if( !SCIPdialogHasEntry(submenu, "problem") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteProblem, NULL, NULL,
            "problem",
            "write original problem to file (format is given by file extension, e.g., orig.{lp,rlp,cip,mps})",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write generic problem */
   if( !SCIPdialogHasEntry(submenu, "genproblem") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteGenProblem, NULL, NULL,
            "genproblem",
            "write original problem with generic names to file (format is given by file extension, e.g., orig.{lp,rlp,cip,mps})",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write solution */
   if( !SCIPdialogHasEntry(submenu, "solution") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteSolution, NULL, NULL,
            "solution", "write best primal solution to file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write finite solution */
   if( !SCIPdialogHasEntry(submenu, "finitesolution") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteFiniteSolution, NULL, NULL,
            "finitesolution", "write best primal solution to file (try to make solution values finite, first)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write mip start */
   if( !SCIPdialogHasEntry(submenu, "mipstart") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteMIPStart, NULL, NULL,
            "mipstart", "write mip start to file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write statistics */
   if( !SCIPdialogHasEntry(submenu, "statistics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteStatistics, NULL, NULL,
            "statistics", "write statistics to file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write transproblem */
   if( !SCIPdialogHasEntry(submenu, "transproblem") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteTransproblem, NULL, NULL,
            "transproblem",
            "write current node transformed problem to file (format is given by file extension, e.g., trans.{lp,rlp,cip,mps})",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write transproblem with generic names */
   if( !SCIPdialogHasEntry(submenu, "gentransproblem") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteGenTransproblem, NULL, NULL,
            "gentransproblem",
            "write current node transformed problem with generic names to file (format is given by file extension, e.g., trans.{lp,rlp,cip,mps})",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write cliquegraph */
   if( !SCIPdialogHasEntry(submenu, "cliquegraph") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecCliquegraph, NULL, NULL,
            "cliquegraph",
            "write graph of cliques and implications of binary variables to GML file (better call after presolving)",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write command line history */
   if( !SCIPdialogHasEntry(submenu, "history") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteCommandHistory, NULL, NULL,
            "history",
            "write command line history to a file (only works if SCIP was compiled with 'readline')",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* validate solve */
   if( !SCIPdialogHasEntry(root, "validatesolve") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecValidateSolve, NULL, NULL,
               "validatesolve",
               "validate the solution against external objective reference interval",
               FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   return SCIP_OKAY;
}

/** if a '/' occurs in the parameter's name, adds a sub menu dialog to the given menu and inserts the parameter dialog
 *  recursively in the sub menu; if no '/' occurs in the name, adds a parameter change dialog into the given dialog menu
 */
static
SCIP_RETCODE addSetParamDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          menu,               /**< dialog menu to insert the parameter into */
   SCIP_PARAM*           param,              /**< parameter to add a dialog for */
   char*                 paramname           /**< parameter name to parse */
   )
{
   char* slash;
   char* dirname;

   assert(paramname != NULL);

   /* check for a '/' */
   slash = strchr(paramname, '/');

   if( slash == NULL )
   {
      /* check, if the corresponding dialog already exists */
      if( !SCIPdialogHasEntry(menu, paramname) )
      {
         SCIP_DIALOG* paramdialog;

         if( SCIPparamIsAdvanced(param) )
         {
            SCIP_DIALOG* advmenu;

            if( !SCIPdialogHasEntry(menu, "advanced") )
            {
               /* if not yet existing, create an advanced sub menu */
               char desc[SCIP_MAXSTRLEN];

               (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "advanced parameters");
               SCIP_CALL( SCIPincludeDialog(scip, &advmenu,
                     NULL,
                     SCIPdialogExecMenu, NULL, NULL, "advanced", desc, TRUE, NULL) );
               SCIP_CALL( SCIPaddDialogEntry(scip, menu, advmenu) );
               SCIP_CALL( SCIPreleaseDialog(scip, &advmenu) );
            }

            /* find the corresponding sub menu */
            (void)SCIPdialogFindEntry(menu, "advanced", &advmenu);
            if( advmenu == NULL )
            {
               SCIPerrorMessage("dialog sub menu not found\n");
               return SCIP_PLUGINNOTFOUND;
            }

            if( !SCIPdialogHasEntry(advmenu, paramname) )
            {
               /* create a parameter change dialog */
               SCIP_CALL( SCIPincludeDialog(scip, &paramdialog,
                     NULL,
                     SCIPdialogExecSetParam, SCIPdialogDescSetParam, NULL,
                     paramname, SCIPparamGetDesc(param), FALSE, (SCIP_DIALOGDATA*)param) );
               SCIP_CALL( SCIPaddDialogEntry(scip, advmenu, paramdialog) );
               SCIP_CALL( SCIPreleaseDialog(scip, &paramdialog) );
            }
         }
         else
         {
            /* create a parameter change dialog */
            SCIP_CALL( SCIPincludeDialog(scip, &paramdialog,
                  NULL,
                  SCIPdialogExecSetParam, SCIPdialogDescSetParam, NULL,
                  paramname, SCIPparamGetDesc(param), FALSE, (SCIP_DIALOGDATA*)param) );
            SCIP_CALL( SCIPaddDialogEntry(scip, menu, paramdialog) );
            SCIP_CALL( SCIPreleaseDialog(scip, &paramdialog) );
         }
      }
   }
   else
   {
      SCIP_DIALOG* submenu;

      /* split the parameter name into dirname and parameter name */
      dirname = paramname;
      paramname = slash+1;
      *slash = '\0';

      /* if not yet existing, create a corresponding sub menu */
      if( !SCIPdialogHasEntry(menu, dirname) )
      {
         char desc[SCIP_MAXSTRLEN];

         (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "parameters for <%s>", dirname);
         SCIP_CALL( SCIPincludeDialog(scip, &submenu,
               NULL,
               SCIPdialogExecMenu, NULL, NULL, dirname, desc, TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, menu, submenu) );
         SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
      }

      /* find the corresponding sub menu */
      (void)SCIPdialogFindEntry(menu, dirname, &submenu);
      if( submenu == NULL )
      {
         SCIPerrorMessage("dialog sub menu not found\n");
         return SCIP_PLUGINNOTFOUND;
      }

      /* recursively call add parameter method */
      SCIP_CALL( addSetParamDialog(scip, submenu, param, paramname) );
   }

   return SCIP_OKAY;
}

/** if a '/' occurs in the parameter's name, adds a sub menu dialog to the given menu and inserts the parameter dialog
 *  recursively in the sub menu; if no '/' occurs in the name, adds a fix parameter dialog into the given dialog menu
 */
static
SCIP_RETCODE addFixParamDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          menu,               /**< dialog menu to insert the parameter into */
   SCIP_PARAM*           param,              /**< parameter to add a dialog for */
   char*                 paramname           /**< parameter name to parse */
   )
{
   char* slash;
   char* dirname;

   assert(paramname != NULL);

   /* check for a '/' */
   slash = strchr(paramname, '/');

   if( slash == NULL )
   {
      /* check, if the corresponding dialog already exists */
      if( !SCIPdialogHasEntry(menu, paramname) )
      {
         SCIP_DIALOG* paramdialog;

         if( SCIPparamIsAdvanced(param) )
         {
            SCIP_DIALOG* advmenu;

            if( !SCIPdialogHasEntry(menu, "advanced") )
            {
               /* if not yet existing, create an advanced sub menu */
               char desc[SCIP_MAXSTRLEN];

               (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "advanced parameters");
               SCIP_CALL( SCIPincludeDialog(scip, &advmenu,
                     NULL,
                     SCIPdialogExecMenu, NULL, NULL, "advanced", desc, TRUE, NULL) );
               SCIP_CALL( SCIPaddDialogEntry(scip, menu, advmenu) );
               SCIP_CALL( SCIPreleaseDialog(scip, &advmenu) );
            }

            /* find the corresponding sub menu */
            (void)SCIPdialogFindEntry(menu, "advanced", &advmenu);
            if( advmenu == NULL )
            {
               SCIPerrorMessage("dialog sub menu not found\n");
               return SCIP_PLUGINNOTFOUND;
            }

            if( !SCIPdialogHasEntry(advmenu, paramname) )
            {
               /* create a fix parameter dialog */
               SCIP_CALL( SCIPincludeDialog(scip, &paramdialog,
                     NULL,
                     SCIPdialogExecFixParam, SCIPdialogDescFixParam, NULL,
                     paramname, SCIPparamGetDesc(param), FALSE, (SCIP_DIALOGDATA*)param) );
               SCIP_CALL( SCIPaddDialogEntry(scip, advmenu, paramdialog) );
               SCIP_CALL( SCIPreleaseDialog(scip, &paramdialog) );
            }
         }
         else
         {
            /* create a fix parameter dialog */
            SCIP_CALL( SCIPincludeDialog(scip, &paramdialog,
                  NULL,
                  SCIPdialogExecFixParam, SCIPdialogDescFixParam, NULL,
                  paramname, SCIPparamGetDesc(param), FALSE, (SCIP_DIALOGDATA*)param) );
            SCIP_CALL( SCIPaddDialogEntry(scip, menu, paramdialog) );
            SCIP_CALL( SCIPreleaseDialog(scip, &paramdialog) );
         }
      }
   }
   else
   {
      SCIP_DIALOG* submenu;

      /* split the parameter name into dirname and parameter name */
      dirname = paramname;
      paramname = slash+1;
      *slash = '\0';

      /* if not yet existing, create a corresponding sub menu */
      if( !SCIPdialogHasEntry(menu, dirname) )
      {
         char desc[SCIP_MAXSTRLEN];

         (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "parameters for <%s>", dirname);
         SCIP_CALL( SCIPincludeDialog(scip, &submenu,
               NULL,
               SCIPdialogExecMenu, NULL, NULL, dirname, desc, TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, menu, submenu) );
         SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
      }

      /* find the corresponding sub menu */
      (void)SCIPdialogFindEntry(menu, dirname, &submenu);
      if( submenu == NULL )
      {
         SCIPerrorMessage("dialog sub menu not found\n");
         return SCIP_PLUGINNOTFOUND;
      }

      /* recursively call add parameter method */
      SCIP_CALL( addFixParamDialog(scip, submenu, param, paramname) );
   }

   return SCIP_OKAY;
}

/** create a "emphasis" sub menu */
static
SCIP_RETCODE createEmphasisSubmenu(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          root,               /**< the menu to add the empty sub menu */
   SCIP_DIALOG**         submenu             /**< pointer to store the created emphasis sub menu */
   )
{
   if( !SCIPdialogHasEntry(root, "emphasis") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, submenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "emphasis", "predefined parameter settings", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, *submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, submenu) );
   }
   else if( SCIPdialogFindEntry(root, "emphasis", submenu) != 1 )
   {
      SCIPerrorMessage("emphasis sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert(*submenu != NULL);

   return SCIP_OKAY;
}


/** includes or updates the "set" menu for each available parameter setting */
SCIP_RETCODE SCIPincludeDialogDefaultSet(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DIALOG* root;
   SCIP_DIALOG* setmenu;
   SCIP_DIALOG* emphasismenu;
   SCIP_DIALOG* submenu;
   SCIP_DIALOG* dialog;
   SCIP_PARAM** params;
   char* paramname;
   int nparams;
   int i;

   SCIP_BRANCHRULE** branchrules; 
   SCIP_CONFLICTHDLR** conflicthdlrs;
   SCIP_CONSHDLR** conshdlrs;
   SCIP_DISP** disps;
   SCIP_HEUR** heurs;
   SCIP_NLPI** nlpis;
   SCIP_NODESEL** nodesels;
   SCIP_PRESOL** presols;
   SCIP_PRICER** pricers;
   SCIP_READER** readers;
   SCIP_SEPA** sepas;
   int nbranchrules;
   int nconflicthdlrs;
   int nconshdlrs;
   int ndisps;
   int nheurs;
   int nnlpis;
   int nnodesels;
   int npresols;
   int npricers;
   int nreaders;
   int nsepas;

   /* get root dialog */
   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      SCIPerrorMessage("root dialog not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* find (or create) the "set" menu of the root dialog */
   if( !SCIPdialogHasEntry(root, "set") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &setmenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "set", "load/save/change parameters", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, setmenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &setmenu) );
   }
   if( SCIPdialogFindEntry(root, "set", &setmenu) != 1 )
   {
      SCIPerrorMessage("set sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* set default */
   if( !SCIPdialogHasEntry(setmenu, "default") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecSetDefault, NULL, NULL,
            "default", "reset parameter settings to their default values", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set load */
   if( !SCIPdialogHasEntry(setmenu, "load") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecSetLoad, NULL, NULL,
            "load", "load parameter settings from a file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set save */
   if( !SCIPdialogHasEntry(setmenu, "save") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecSetSave, NULL, NULL,
            "save", "save parameter settings to a file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set diffsave */
   if( !SCIPdialogHasEntry(setmenu, "diffsave") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecSetDiffsave, NULL, NULL,
            "diffsave", "save non-default parameter settings to a file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set branching */
   if( !SCIPdialogHasEntry(setmenu, "branching") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "branching", "change parameters for branching rules", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "branching", &submenu) != 1 )
   {
      SCIPerrorMessage("branching sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nbranchrules = SCIPgetNBranchrules(scip);
   branchrules = SCIPgetBranchrules(scip);

   for( i = 0; i < nbranchrules; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPbranchruleGetName(branchrules[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPbranchruleGetName(branchrules[i]), SCIPbranchruleGetDesc(branchrules[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set branching priority */
   if( !SCIPdialogHasEntry(submenu, "priority") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecSetBranchingPriority, NULL, NULL,
            "priority", "change branching priority of a single variable", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set branching direction */
   if( !SCIPdialogHasEntry(submenu, "direction") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecSetBranchingDirection, NULL, NULL,
            "direction", "change preferred branching direction of a single variable (-1:down, 0:auto, +1:up)",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set conflict */
   if( !SCIPdialogHasEntry(setmenu, "conflict") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "conflict", "change parameters for conflict handlers", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "conflict", &submenu) != 1 )
   {
      SCIPerrorMessage("conflict sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nconflicthdlrs = SCIPgetNConflicthdlrs(scip);
   conflicthdlrs = SCIPgetConflicthdlrs(scip);

   for( i = 0; i < nconflicthdlrs; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPconflicthdlrGetName(conflicthdlrs[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPconflicthdlrGetName(conflicthdlrs[i]), SCIPconflicthdlrGetDesc(conflicthdlrs[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set constraints */
   if( !SCIPdialogHasEntry(setmenu, "constraints") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "constraints", "change parameters for constraint handlers", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "constraints", &submenu) != 1 )
   {
      SCIPerrorMessage("constraints sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);

   for( i = 0; i < nconshdlrs; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPconshdlrGetName(conshdlrs[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPconshdlrGetName(conshdlrs[i]), SCIPconshdlrGetDesc(conshdlrs[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set display */
   if( !SCIPdialogHasEntry(setmenu, "display") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "display", "change parameters for display columns", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "display", &submenu) != 1 )
   {
      SCIPerrorMessage("display sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   ndisps = SCIPgetNDisps(scip);
   disps = SCIPgetDisps(scip);

   for( i = 0; i < ndisps; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPdispGetName(disps[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPdispGetName(disps[i]), SCIPdispGetDesc(disps[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set heuristics */
   if( !SCIPdialogHasEntry(setmenu, "heuristics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "heuristics", "change parameters for primal heuristics", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "heuristics", &submenu) != 1 )
   {
      SCIPerrorMessage("heuristics sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nheurs = SCIPgetNHeurs(scip);
   heurs = SCIPgetHeurs(scip);

   for( i = 0; i < nheurs; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPheurGetName(heurs[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPheurGetName(heurs[i]), SCIPheurGetDesc(heurs[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* create set heuristics emphasis */
   SCIP_CALL( createEmphasisSubmenu(scip, submenu, &emphasismenu) );
   assert(emphasismenu != NULL);

   /* set heuristics emphasis aggressive */
   if( !SCIPdialogHasEntry(emphasismenu, "aggressive") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetHeuristicsAggressive, NULL, NULL,
            "aggressive", "sets heuristics <aggressive>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set heuristics emphasis default */
   if( !SCIPdialogHasEntry(emphasismenu, "default") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetHeuristicsDefault, NULL, NULL,
            "default", "sets heuristics settings to <default> ", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set heuristics emphasis fast */
   if( !SCIPdialogHasEntry(emphasismenu, "fast") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetHeuristicsFast, NULL, NULL,
            "fast", "sets heuristics <fast>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set heuristics emphasis off */
   if( !SCIPdialogHasEntry(emphasismenu, "off") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetHeuristicsOff, NULL, NULL,
            "off", "turns <off> all heuristics", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set limits */
   if( !SCIPdialogHasEntry(setmenu, "limits") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "limits", "change parameters for time, memory, objective value, and other limits", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );

      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecSetLimitsObjective, NULL, NULL,
            "objective", "set limit on objective function, such that only solutions better than this limit are accepted", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );

      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set LP */
   if( !SCIPdialogHasEntry(setmenu, "lp") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "lp", "change parameters for linear programming relaxations", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set NLP */
   if( !SCIPdialogHasEntry(setmenu, "nlp") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "nlp", "change parameters for nonlinear programming relaxations", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set memory */
   if( !SCIPdialogHasEntry(setmenu, "memory") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "memory", "change parameters for memory management", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set misc */
   if( !SCIPdialogHasEntry(setmenu, "misc") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "misc", "change parameters for miscellaneous stuff", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set nlpi */
   if( !SCIPdialogHasEntry(setmenu, "nlpi") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "nlpi", "change parameters for NLP solver interfaces", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "nlpi", &submenu) != 1 )
   {
      SCIPerrorMessage("nlpi sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nnlpis = SCIPgetNNlpis(scip);
   nlpis = SCIPgetNlpis(scip);

   for( i = 0; i < nnlpis; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPnlpiGetName(nlpis[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPnlpiGetName(nlpis[i]), SCIPnlpiGetDesc(nlpis[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set nodeselection */
   if( !SCIPdialogHasEntry(setmenu, "nodeselection") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "nodeselection", "change parameters for node selectors", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "nodeselection", &submenu) != 1 )
   {
      SCIPerrorMessage("nodeselection sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nnodesels = SCIPgetNNodesels(scip);
   nodesels = SCIPgetNodesels(scip);

   for( i = 0; i < nnodesels; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPnodeselGetName(nodesels[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPnodeselGetName(nodesels[i]), SCIPnodeselGetDesc(nodesels[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set numerics */
   if( !SCIPdialogHasEntry(setmenu, "numerics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "numerics", "change parameters for numerical values", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set parallel */
   if( !SCIPdialogHasEntry(setmenu, "parallel") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "parallel", "change parameters for parallel implementation", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set presolving */
   if( !SCIPdialogHasEntry(setmenu, "presolving") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "presolving", "change parameters for presolving", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "presolving", &submenu) != 1 )
   {
      SCIPerrorMessage("presolving sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   npresols = SCIPgetNPresols(scip);
   presols = SCIPgetPresols(scip);

   for( i = 0; i < npresols; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPpresolGetName(presols[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL, SCIPdialogExecMenu, NULL, NULL,
               SCIPpresolGetName(presols[i]), SCIPpresolGetDesc(presols[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* create set presolving emphasis */
   SCIP_CALL( createEmphasisSubmenu(scip, submenu, &emphasismenu) );
   assert(emphasismenu != NULL);

   /* set presolving emphasis aggressive */
   if( !SCIPdialogHasEntry(emphasismenu, "aggressive") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetPresolvingAggressive, NULL, NULL,
            "aggressive", "sets presolving <aggressive>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set presolving emphasis default */
   if( !SCIPdialogHasEntry(emphasismenu, "default") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetPresolvingDefault, NULL, NULL,
            "default", "sets presolving settings to <default>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set presolving emphasis fast */
   if( !SCIPdialogHasEntry(emphasismenu, "fast") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetPresolvingFast, NULL, NULL,
            "fast", "sets presolving <fast>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set presolving emphasis off */
   if( !SCIPdialogHasEntry(emphasismenu, "off") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetPresolvingOff, NULL, NULL,
            "off", "turns <off> all presolving", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set pricing */
   if( !SCIPdialogHasEntry(setmenu, "pricing") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "pricing", "change parameters for pricing variables", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "pricing", &submenu) != 1 )
   {
      SCIPerrorMessage("pricing sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   npricers = SCIPgetNPricers(scip);
   pricers = SCIPgetPricers(scip);

   for( i = 0; i < npricers; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPpricerGetName(pricers[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPpricerGetName(pricers[i]), SCIPpricerGetDesc(pricers[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set propagation */
   if( !SCIPdialogHasEntry(setmenu, "propagating") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "propagating", "change parameters for constraint propagation", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set reading */
   if( !SCIPdialogHasEntry(setmenu, "reading") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "reading", "change parameters for problem file readers", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "reading", &submenu) != 1 )
   {
      SCIPerrorMessage("reading sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nreaders = SCIPgetNReaders(scip);
   readers = SCIPgetReaders(scip);

   for( i = 0; i < nreaders; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPreaderGetName(readers[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPreaderGetName(readers[i]), SCIPreaderGetDesc(readers[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set separating */
   if( !SCIPdialogHasEntry(setmenu, "separating") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "separating", "change parameters for cut separators", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "separating", &submenu) != 1 )
   {
      SCIPerrorMessage("separating sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nsepas = SCIPgetNSepas(scip);
   sepas = SCIPgetSepas(scip);

   for( i = 0; i < nsepas; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPsepaGetName(sepas[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL, SCIPdialogExecMenu, NULL, NULL,
               SCIPsepaGetName(sepas[i]), SCIPsepaGetDesc(sepas[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* create set separating emphasis */
   SCIP_CALL( createEmphasisSubmenu(scip, submenu, &emphasismenu) );
   assert(emphasismenu != NULL);

   /* set separating emphasis aggressive */
   if( !SCIPdialogHasEntry(emphasismenu, "aggressive") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetSeparatingAggressive, NULL, NULL,
            "aggressive", "sets separating <aggressive>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set separating emphasis default */
   if( !SCIPdialogHasEntry(emphasismenu, "default") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetSeparatingDefault, NULL, NULL,
            "default", "sets separating settings to <default>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set separating emphasis fast */
   if( !SCIPdialogHasEntry(emphasismenu, "fast") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetSeparatingFast, NULL, NULL,
            "fast", "sets separating <fast>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set separating emphasis off */
   if( !SCIPdialogHasEntry(emphasismenu, "off") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetSeparatingOff, NULL, NULL,
            "off", "turns <off> all separation", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set timing */
   if( !SCIPdialogHasEntry(setmenu, "timing") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "timing", "change parameters for timing issues", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set visualization */
   if( !SCIPdialogHasEntry(setmenu, "visual") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "visual", "change parameters for visualization output", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set emphasis */
   SCIP_CALL( createEmphasisSubmenu(scip, setmenu, &submenu) );

   /* get SCIP's parameters */
   params = SCIPgetParams(scip);
   nparams = SCIPgetNParams(scip);

   /* insert each parameter into the set menu */
   for( i = 0; i < nparams; ++i )
   {
      const char* pname;

      pname = SCIPparamGetName(params[i]);
      SCIP_ALLOC( BMSduplicateMemoryArray(&paramname, pname, strlen(pname)+1) );
      SCIP_CALL( addSetParamDialog(scip, setmenu, params[i], paramname) );
      BMSfreeMemoryArray(&paramname);
   }

   /* set emphasis feasibility */
   /* add "counter" dialog to "set/emphasis" sub menu */
   if( !SCIPdialogHasEntry(submenu, "counter") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecSetEmphasisCounter, NULL, NULL,
            "counter", "predefined parameter settings for a \"feasible\" and \"fast\" counting process", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* add "cpsolver" dialog to "set/emphasis" sub menu */
   if( !SCIPdialogHasEntry(submenu, "cpsolver") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecSetEmphasisCpsolver, NULL, NULL,
            "cpsolver", "predefined parameter settings for CP like search", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* add "easycip" dialog to "set/emphasis" sub menu */
   if( !SCIPdialogHasEntry(submenu, "easycip") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecSetEmphasisEasycip, NULL, NULL,
            "easycip", "predefined parameter settings for easy problems", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* add "feasibility" dialog to "set/emphasis" sub menu */
   if( !SCIPdialogHasEntry(submenu, "feasibility") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecSetEmphasisFeasibility, NULL, NULL,
            "feasibility", "predefined parameter settings for feasibility problems", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* add "hardlp" dialog to "set/emphasis" sub menu */
   if( !SCIPdialogHasEntry(submenu, "hardlp") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecSetEmphasisHardlp, NULL, NULL,
            "hardlp", "predefined parameter settings for problems with a hard LP", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* add "optimality" dialog to "set/emphasis" sub menu */
   if( !SCIPdialogHasEntry(submenu, "optimality") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecSetEmphasisOptimality, NULL, NULL,
            "optimality", "predefined parameter settings for proving optimality fast", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   return SCIP_OKAY;
}

/** includes or updates the "fix" menu for each available parameter setting */
SCIP_RETCODE SCIPincludeDialogDefaultFix(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DIALOG* root;
   SCIP_DIALOG* fixmenu;
   SCIP_DIALOG* submenu;
   SCIP_DIALOG* dialog;
   SCIP_PARAM** params;
   char* paramname;
   int nparams;
   int i;

   SCIP_BRANCHRULE** branchrules;
   SCIP_CONFLICTHDLR** conflicthdlrs;
   SCIP_CONSHDLR** conshdlrs;
   SCIP_DISP** disps;
   SCIP_HEUR** heurs;
   SCIP_NLPI** nlpis;
   SCIP_NODESEL** nodesels;
   SCIP_PRESOL** presols;
   SCIP_PRICER** pricers;
   SCIP_READER** readers;
   SCIP_SEPA** sepas;
   int nbranchrules;
   int nconflicthdlrs;
   int nconshdlrs;
   int ndisps;
   int nheurs;
   int nnlpis;
   int nnodesels;
   int npresols;
   int npricers;
   int nreaders;
   int nsepas;

   /* get root dialog */
   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      SCIPerrorMessage("root dialog not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* find (or create) the "fix" menu of the root dialog */
   if( !SCIPdialogHasEntry(root, "fix") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &fixmenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "fix", "fix/unfix parameters", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, fixmenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &fixmenu) );
   }
   if( SCIPdialogFindEntry(root, "fix", &fixmenu) != 1 )
   {
      SCIPerrorMessage("fix sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* fix branching */
   if( !SCIPdialogHasEntry(fixmenu, "branching") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "branching", "fix parameters for branching rules", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "branching", &submenu) != 1 )
   {
      SCIPerrorMessage("branching sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nbranchrules = SCIPgetNBranchrules(scip);
   branchrules = SCIPgetBranchrules(scip);

   for( i = 0; i < nbranchrules; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPbranchruleGetName(branchrules[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPbranchruleGetName(branchrules[i]), SCIPbranchruleGetDesc(branchrules[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix conflict */
   if( !SCIPdialogHasEntry(fixmenu, "conflict") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "conflict", "fix parameters for conflict handlers", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "conflict", &submenu) != 1 )
   {
      SCIPerrorMessage("conflict sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nconflicthdlrs = SCIPgetNConflicthdlrs(scip);
   conflicthdlrs = SCIPgetConflicthdlrs(scip);

   for( i = 0; i < nconflicthdlrs; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPconflicthdlrGetName(conflicthdlrs[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPconflicthdlrGetName(conflicthdlrs[i]), SCIPconflicthdlrGetDesc(conflicthdlrs[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix constraints */
   if( !SCIPdialogHasEntry(fixmenu, "constraints") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "constraints", "fix parameters for constraint handlers", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "constraints", &submenu) != 1 )
   {
      SCIPerrorMessage("constraints sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);

   for( i = 0; i < nconshdlrs; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPconshdlrGetName(conshdlrs[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPconshdlrGetName(conshdlrs[i]), SCIPconshdlrGetDesc(conshdlrs[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix display */
   if( !SCIPdialogHasEntry(fixmenu, "display") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "display", "fix parameters for display columns", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "display", &submenu) != 1 )
   {
      SCIPerrorMessage("display sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   ndisps = SCIPgetNDisps(scip);
   disps = SCIPgetDisps(scip);

   for( i = 0; i < ndisps; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPdispGetName(disps[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPdispGetName(disps[i]), SCIPdispGetDesc(disps[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix heuristics */
   if( !SCIPdialogHasEntry(fixmenu, "heuristics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "heuristics", "fix parameters for primal heuristics", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "heuristics", &submenu) != 1 )
   {
      SCIPerrorMessage("heuristics sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nheurs = SCIPgetNHeurs(scip);
   heurs = SCIPgetHeurs(scip);

   for( i = 0; i < nheurs; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPheurGetName(heurs[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPheurGetName(heurs[i]), SCIPheurGetDesc(heurs[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix limits */
   if( !SCIPdialogHasEntry(fixmenu, "limits") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "limits", "fix parameters for time, memory, objective value, and other limits", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );

      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* fix LP */
   if( !SCIPdialogHasEntry(fixmenu, "lp") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "lp", "fix parameters for linear programming relaxations", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* fix NLP */
   if( !SCIPdialogHasEntry(fixmenu, "nlp") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "nlp", "fix parameters for nonlinear programming relaxations", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* fix memory */
   if( !SCIPdialogHasEntry(fixmenu, "memory") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "memory", "fix parameters for memory management", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* fix misc */
   if( !SCIPdialogHasEntry(fixmenu, "misc") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "misc", "fix parameters for miscellaneous stuff", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* fix nlpi */
   if( !SCIPdialogHasEntry(fixmenu, "nlpi") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "nlpi", "fix parameters for NLP solver interfaces", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "nlpi", &submenu) != 1 )
   {
      SCIPerrorMessage("nlpi sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nnlpis = SCIPgetNNlpis(scip);
   nlpis = SCIPgetNlpis(scip);

   for( i = 0; i < nnlpis; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPnlpiGetName(nlpis[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPnlpiGetName(nlpis[i]), SCIPnlpiGetDesc(nlpis[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix nodeselection */
   if( !SCIPdialogHasEntry(fixmenu, "nodeselection") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "nodeselection", "fix parameters for node selectors", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "nodeselection", &submenu) != 1 )
   {
      SCIPerrorMessage("nodeselection sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nnodesels = SCIPgetNNodesels(scip);
   nodesels = SCIPgetNodesels(scip);

   for( i = 0; i < nnodesels; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPnodeselGetName(nodesels[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPnodeselGetName(nodesels[i]), SCIPnodeselGetDesc(nodesels[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix numerics */
   if( !SCIPdialogHasEntry(fixmenu, "numerics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "numerics", "fix parameters for numerical values", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* fix presolving */
   if( !SCIPdialogHasEntry(fixmenu, "presolving") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "presolving", "fix parameters for presolving", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "presolving", &submenu) != 1 )
   {
      SCIPerrorMessage("presolving sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   npresols = SCIPgetNPresols(scip);
   presols = SCIPgetPresols(scip);

   for( i = 0; i < npresols; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPpresolGetName(presols[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL, SCIPdialogExecMenu, NULL, NULL,
               SCIPpresolGetName(presols[i]), SCIPpresolGetDesc(presols[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix pricing */
   if( !SCIPdialogHasEntry(fixmenu, "pricing") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "pricing", "fix parameters for pricing variables", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "pricing", &submenu) != 1 )
   {
      SCIPerrorMessage("pricing sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   npricers = SCIPgetNPricers(scip);
   pricers = SCIPgetPricers(scip);

   for( i = 0; i < npricers; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPpricerGetName(pricers[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPpricerGetName(pricers[i]), SCIPpricerGetDesc(pricers[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix propagation */
   if( !SCIPdialogHasEntry(fixmenu, "propagating") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "propagating", "fix parameters for constraint propagation", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* fix reading */
   if( !SCIPdialogHasEntry(fixmenu, "reading") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "reading", "fix parameters for problem file readers", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "reading", &submenu) != 1 )
   {
      SCIPerrorMessage("reading sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nreaders = SCIPgetNReaders(scip);
   readers = SCIPgetReaders(scip);

   for( i = 0; i < nreaders; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPreaderGetName(readers[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL,
               SCIPdialogExecMenu, NULL, NULL,
               SCIPreaderGetName(readers[i]), SCIPreaderGetDesc(readers[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix separating */
   if( !SCIPdialogHasEntry(fixmenu, "separating") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "separating", "fix parameters for cut separators", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(fixmenu, "separating", &submenu) != 1 )
   {
      SCIPerrorMessage("separating sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   nsepas = SCIPgetNSepas(scip);
   sepas = SCIPgetSepas(scip);

   for( i = 0; i < nsepas; ++i )
   {
      if( !SCIPdialogHasEntry(submenu, SCIPsepaGetName(sepas[i])) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
               NULL, SCIPdialogExecMenu, NULL, NULL,
               SCIPsepaGetName(sepas[i]), SCIPsepaGetDesc(sepas[i]), TRUE, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* fix timing */
   if( !SCIPdialogHasEntry(fixmenu, "timing") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "timing", "fix parameters for timing issues", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, fixmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }

   /* get SCIP's parameters */
   params = SCIPgetParams(scip);
   nparams = SCIPgetNParams(scip);

   /* insert each parameter into the fix menu */
   for( i = 0; i < nparams; ++i )
   {
      const char* pname;

      pname = SCIPparamGetName(params[i]);
      SCIP_ALLOC( BMSduplicateMemoryArray(&paramname, pname, strlen(pname)+1) );
      SCIP_CALL( addFixParamDialog(scip, fixmenu, params[i], paramname) );
      BMSfreeMemoryArray(&paramname);
   }

   return SCIP_OKAY;
}
