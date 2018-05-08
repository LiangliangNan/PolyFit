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

/**@file   visual.c
 * @brief  methods for creating output for visualization tools (VBC, BAK)
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 *
 * Output can be generated for the following visualization tools:
 *
 * - VBCTOOL - a graphical interface for Visualization of Branch Cut algorithms @n
 *   See <a href="http://www.informatik.uni-koeln.de/ls_juenger/research/vbctool">VBCTOOL</a>.
 * - BAK: Branch-and-bound Analysis Kit @n
 *   BAK is available through COIN-OR, see <a href="https://projects.coin-or.org/CoinBazaar/wiki/Projects/BAK">BAK</a>.
 *   A description is <a href="http://www.optimization-online.org/DB_HTML/2007/09/1785.html">available</a> as well.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>

#include "blockmemshell/memory.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/var.h"
#include "scip/tree.h"
#include "scip/visual.h"
#include "scip/struct_visual.h"


/** returns the branching variable of the node, or NULL */
static
void getBranchInfo(
   SCIP_NODE*            node,               /**< node */
   SCIP_VAR**            var,                /**< pointer to store the branching variable */
   SCIP_BOUNDTYPE*       boundtype,          /**< pointer to store the branching type: lower or upper bound */
   SCIP_Real*            bound               /**< pointer to store the new bound of the branching variable */
   )
{
   SCIP_DOMCHGBOUND* domchgbound;

   (*var) = NULL;
   (*bound) = 0.0;
   (*boundtype) = SCIP_BOUNDTYPE_LOWER;

   assert(node != NULL);
   if( node->domchg == NULL )
      return;

   domchgbound = &node->domchg->domchgbound;
   if( domchgbound->nboundchgs == 0 )
      return;

   (*var) = domchgbound->boundchgs[0].var;
   (*bound) = domchgbound->boundchgs[0].newbound;
   (*boundtype) = (SCIP_BOUNDTYPE) domchgbound->boundchgs[0].boundtype;
}

/** creates visualization data structure */
SCIP_RETCODE SCIPvisualCreate(
   SCIP_VISUAL**         visual,             /**< pointer to store visualization information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_ALLOC( BMSallocMemory(visual) );

   (*visual)->vbcfile = NULL;
   (*visual)->bakfile = NULL;
   (*visual)->messagehdlr = messagehdlr;
   (*visual)->nodenum = NULL;
   (*visual)->timestep = 0;
   (*visual)->lastnode = NULL;
   (*visual)->lastcolor = SCIP_VBCCOLOR_NONE;
   (*visual)->userealtime = FALSE;

   return SCIP_OKAY;
}

/** frees visualization data structure */
void SCIPvisualFree(
   SCIP_VISUAL**         visual              /**< pointer to store visualization information */
   )
{
   assert( visual != NULL );
   assert( *visual != NULL );
   assert( (*visual)->vbcfile == NULL );
   assert( (*visual)->bakfile == NULL );
   assert( (*visual)->nodenum == NULL );

   BMSfreeMemory(visual);
}

/** initializes visualization information and creates a file for visualization output */
SCIP_RETCODE SCIPvisualInit(
   SCIP_VISUAL*          visual,             /**< visualization information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert( visual != NULL );
   assert( set != NULL );
   assert( set->visual_vbcfilename != NULL );
   assert( set->visual_bakfilename != NULL );
   assert( visual->nodenum == NULL );

   /* check whether we should initialize VBC output */
   if ( set->visual_vbcfilename[0] != '-' || set->visual_vbcfilename[1] != '\0' )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
         "storing VBC information in file <%s>\n", set->visual_vbcfilename);
      visual->vbcfile = fopen(set->visual_vbcfilename, "w");
      visual->timestep = 0;
      visual->lastnode = NULL;
      visual->lastcolor = SCIP_VBCCOLOR_NONE;
      visual->userealtime = set->visual_realtime;

      if( visual->vbcfile == NULL )
      {
         SCIPerrorMessage("error creating file <%s>\n", set->visual_vbcfilename);
         SCIPprintSysError(set->visual_vbcfilename);
         return SCIP_FILECREATEERROR;
      }

      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#TYPE: COMPLETE TREE\n");
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#TIME: SET\n");
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#BOUNDS: SET\n");
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#INFORMATION: STANDARD\n");
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#NODE_NUMBER: NONE\n");
   }

   /* check whether we should initialize BAK output */
   if ( set->visual_bakfilename[0] != '-' || set->visual_bakfilename[1] != '\0' )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
         "storing BAK information in file <%s>\n", set->visual_bakfilename);
      visual->bakfile = fopen(set->visual_bakfilename, "w");
      visual->timestep = 0;
      visual->lastnode = NULL;
      visual->lastcolor = SCIP_VBCCOLOR_NONE;
      visual->userealtime = set->visual_realtime;

      if ( visual->bakfile == NULL )
      {
         SCIPerrorMessage("error creating file <%s>\n", set->visual_bakfilename);
         SCIPprintSysError(set->visual_bakfilename);
         return SCIP_FILECREATEERROR;
      }
   }

   /* possibly init hashmap for nodes */
   if ( visual->vbcfile != NULL || visual->bakfile != NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&visual->nodenum, blkmem, SCIP_HASHSIZE_VBC) );
   }

   return SCIP_OKAY;
}

/** closes the visualization output file */
void SCIPvisualExit(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert( visual != NULL );
   assert( set != NULL );

   if ( visual->vbcfile != NULL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL, "closing VBC information file\n");

      fclose(visual->vbcfile);
      visual->vbcfile = NULL;
   }

   if ( visual->bakfile != NULL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL, "closing BAK information file\n");

      fclose(visual->bakfile);
      visual->bakfile = NULL;
   }

   if ( visual->nodenum )
      SCIPhashmapFree(&visual->nodenum);
}

/** prints current solution time to visualization output file */
static
void printTime(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             vbc                 /**< whether we use vbc output (bak otherwise) */
   )
{
   SCIP_Longint step;
   int hours;
   int mins;
   int secs;
   int hunds;

   assert( visual != NULL );
   assert( stat != NULL );

   if( visual->userealtime )
   {
      double time;
      time = SCIPclockGetTime(stat->solvingtime);
      step = (SCIP_Longint)(time * 100.0);
   }
   else
   {
      step = visual->timestep;
      visual->timestep++;
   }

   if ( vbc )
   {
      hours = (int)(step / (60*60*100));
      step %= 60*60*100;
      mins = (int)(step / (60*100));
      step %= 60*100;
      secs = (int)(step / 100);
      step %= 100;
      hunds = (int)step;

      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "%02d:%02d:%02d.%02d ", hours, mins, secs, hunds);
   }
   else
   {
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "%f ", (SCIP_Real) step/100.0);
   }
}

/** creates a new node entry in the visualization output file */
SCIP_RETCODE SCIPvisualNewChild(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   SCIP_Real lowerbound;
   size_t parentnodenum;
   size_t nodenum;

   assert( visual != NULL );
   assert( stat != NULL );
   assert( node != NULL );

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* check whether output should be created */
   if ( visual->vbcfile == NULL && visual->bakfile == NULL )
      return SCIP_OKAY;

   /* insert mapping node -> nodenum into hash map */
   if( stat->ncreatednodesrun >= (SCIP_Longint)INT_MAX )
   {
      SCIPerrorMessage("too many nodes to store in the visualization file\n");
      return SCIP_INVALIDDATA;
   }

   nodenum = (size_t)stat->ncreatednodesrun;
   assert(nodenum > 0);
   SCIP_CALL( SCIPhashmapSetImage(visual->nodenum, node, (void*)nodenum) );

   /* get nodenum of parent node from hash map */
   parentnodenum = (node->parent != NULL ? (size_t)SCIPhashmapGetImage(visual->nodenum, node->parent) : 0);
   assert(node->parent == NULL || parentnodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   /* determine lower bound */
   if ( set->visual_objextern )
      lowerbound = SCIPretransformObj(set->scip, SCIPnodeGetLowerbound(node));
   else
      lowerbound = SCIPnodeGetLowerbound(node);

   if ( visual->vbcfile != NULL )
   {
      printTime(visual, stat, TRUE);
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "N %d %d %d\n", (int)parentnodenum, (int)nodenum, SCIP_VBCCOLOR_UNSOLVED);
      printTime(visual, stat, TRUE);
      if( branchvar != NULL )
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
            SCIPvarGetName(branchvar), SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
            branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, lowerbound);
      }
      else
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), lowerbound);
      }
   }

   /* For BAK, not all available information is available here. Use SCIPvisualUpdateChild() instead */

   return SCIP_OKAY;
}

/** updates a node entry in the visualization output file */
SCIP_RETCODE SCIPvisualUpdateChild(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   SCIP_Real lowerbound;
   size_t nodenum;

   assert( visual != NULL );
   assert( stat != NULL );
   assert( node != NULL );

   /* check whether output should be created */
   if ( visual->vbcfile == NULL && visual->bakfile == NULL )
      return SCIP_OKAY;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   /* determine lower bound */
   if ( set->visual_objextern )
      lowerbound = SCIPretransformObj(set->scip, SCIPnodeGetLowerbound(node));
   else
      lowerbound = SCIPnodeGetLowerbound(node);

   if ( visual->vbcfile != NULL )
   {
      printTime(visual, stat, TRUE);
      if( branchvar != NULL )
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
            SCIPvarGetName(branchvar), SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
            branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, lowerbound);
      }
      else
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), lowerbound);
      }
   }

   if ( visual->bakfile != NULL )
   {
      size_t parentnodenum;
      SCIP_Real* lpcandsfrac;
      SCIP_Real sum = 0.0;
      int nlpcands = 0;
      char t = 'M';
      const char* nodeinfo;
      int j;

      /* determine branching type */
      if ( branchvar != NULL )
         t = (branchtype == SCIP_BOUNDTYPE_LOWER ? 'R' : 'L');

      /* get nodenum of parent node from hash map */
      parentnodenum = (node->parent != NULL ? (size_t)SCIPhashmapGetImage(visual->nodenum, node->parent) : 0);
      assert(node->parent == NULL || parentnodenum > 0);

      /* update info depending on the node type */
      switch( SCIPnodeGetType(node) )
      {
      case SCIP_NODETYPE_CHILD:
         /* the child is a new candidate */
         nodeinfo = "candidate";
         break;
      case SCIP_NODETYPE_FOCUSNODE:
         /* the focus node is updated to a branch node */
         nodeinfo = "branched";

         /* calculate infeasibility information only if the LP was solved to optimality */
         if( SCIPgetLPSolstat(set->scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIP_CALL( SCIPgetLPBranchCands(set->scip, NULL, NULL, &lpcandsfrac, &nlpcands, NULL, NULL) );
            for( j = 0; j < nlpcands; ++j )
               sum += lpcandsfrac[j];
         }

         break;
      default:
         SCIPerrorMessage("Error: Unexpected node type <%d> in Update Child Method", SCIPnodeGetType(node));
         return SCIP_INVALIDDATA;
      } /*lint !e788*/
      /* append new status line with updated node information to the bakfile */
      printTime(visual, stat, FALSE);
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "%s %d %d %c %f %f %d\n", nodeinfo, (int)nodenum, (int)parentnodenum, t,
            lowerbound, sum, nlpcands);
   }

   return SCIP_OKAY;
}

/** changes the color of the node to the given color */
static
void vbcSetColor(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node to change color for */
   SCIP_VBCCOLOR         color               /**< new color of node, or SCIP_VBCCOLOR_NONE */
   )
{
   assert( visual != NULL );
   assert( node != NULL );

   if( visual->vbcfile != NULL && color != SCIP_VBCCOLOR_NONE && (node != visual->lastnode || color != visual->lastcolor) )
   {
      size_t nodenum;

      nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
      assert(nodenum > 0);
      printTime(visual, stat, TRUE);
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "P %d %d\n", (int)nodenum, color);
      visual->lastnode = node;
      visual->lastcolor = color;
   }
}

/** marks node as solved in visualization output file */
void SCIPvisualSolvedNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was solved */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   SCIP_Real lowerbound;
   size_t nodenum;

   assert( visual != NULL );
   assert( stat != NULL );
   assert( node != NULL );

   /* check whether output should be created */
   if ( visual->vbcfile == NULL && visual->bakfile == NULL )
      return;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   /* determine lower bound */
   if ( set->visual_objextern )
      lowerbound = SCIPretransformObj(set->scip, SCIPnodeGetLowerbound(node));
   else
      lowerbound = SCIPnodeGetLowerbound(node);

   if ( visual->vbcfile != NULL )
   {
      printTime(visual, stat, TRUE);
      if( branchvar != NULL )
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\\nnr:\\t%" SCIP_LONGINT_FORMAT "\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
            SCIPvarGetName(branchvar),  SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
            branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, lowerbound, stat->nnodes);
      }
      else
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\\nnr:\\t%" SCIP_LONGINT_FORMAT "\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), lowerbound, stat->nnodes);
      }
      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_SOLVED);
   }

   /* do nothing for BAK */
}

/** changes the color of the node to the color of cutoff nodes */
void SCIPvisualCutoffNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node, that was cut off */
   SCIP_Bool             infeasible          /**< whether the node is infeasible (otherwise exceeded the cutoff bound) */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   SCIP_Real lowerbound;
   size_t nodenum;

   assert( visual != NULL );
   assert( stat != NULL );
   assert( node != NULL );

   /* check whether output should be created */
   if ( visual->vbcfile == NULL && visual->bakfile == NULL )
      return;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   /* determine lower bound */
   if ( set->visual_objextern )
      lowerbound = SCIPretransformObj(set->scip, SCIPnodeGetLowerbound(node));
   else
      lowerbound = SCIPnodeGetLowerbound(node);

   if ( visual->vbcfile != NULL )
   {
      printTime(visual, stat, TRUE);
      if( branchvar != NULL )
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\\nnr:\\t%" SCIP_LONGINT_FORMAT "\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
            SCIPvarGetName(branchvar),  SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
            branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, lowerbound, stat->nnodes);
      }
      else
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\\nnr:\\t%" SCIP_LONGINT_FORMAT "\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), lowerbound, stat->nnodes);
      }
      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_CUTOFF);
   }

   if ( visual->bakfile != NULL )
   {
      size_t parentnodenum;
      char t = 'M';

      /* determine branching type */
      if ( branchvar != NULL )
         t = (branchtype == SCIP_BOUNDTYPE_LOWER ? 'R' : 'L');

      /* get nodenum of parent node from hash map */
      parentnodenum = (node->parent != NULL ? (size_t)SCIPhashmapGetImage(visual->nodenum, node->parent) : 0);
      assert(node->parent == NULL || parentnodenum > 0);

      printTime(visual, stat, FALSE);
      if ( infeasible )
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "infeasible %d %d %c\n", (int)nodenum, (int)parentnodenum, t);
      else
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "fathomed %d %d %c\n", (int)nodenum, (int)parentnodenum, t);
   }
}

/** changes the color of the node to the color of nodes where a conflict constraint was found */
void SCIPvisualFoundConflict(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, where the conflict was found */
   )
{
   assert(node != NULL);

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_CONFLICT);

   /* do nothing for BAK */
}

/** changes the color of the node to the color of nodes that were marked to be repropagated */
void SCIPvisualMarkedRepropagateNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was marked to be repropagated */
   )
{
   assert(node != NULL);

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* if the node number is zero, then SCIP is currently in probing and wants to mark a probing node; however this node
    * is not part of the search tree */
   if( SCIPnodeGetNumber(node) > 0 )
      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_MARKREPROP);

   /* do nothing for BAK */
}

/** changes the color of the node to the color of repropagated nodes */
void SCIPvisualRepropagatedNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was repropagated */
   )
{
   assert(node != NULL);

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_REPROP);

   /* do nothing for BAK */
}

/** changes the color of the node to the color of nodes with a primal solution */
void SCIPvisualFoundSolution(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node where the solution was found, or NULL */
   SCIP_Bool             bettersol,          /**< the solution was better than the previous ones */
   SCIP_SOL*             sol                 /**< solution that has been found */
   )
{
   if( node == NULL || ! set->visual_dispsols )
      return;

   if( visual->vbcfile != NULL )
   {
      SCIP_Real obj;
      size_t nodenum;

      /* if we are in probing, determine original parent node */
      while ( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
         node = SCIPnodeGetParent(node);

      /* get node num from hash map */
      assert(node != NULL);
      nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
      assert(nodenum > 0);

      /* get objective of solution */
      if( set->visual_objextern )
         obj = SCIPgetSolOrigObj(set->scip, sol);
      else
         obj = SCIPgetSolTransObj(set->scip, sol);

      printTime(visual, stat, TRUE);
      if( bettersol )
      {
         /* note that this output is in addition to the one by SCIPvisualUpperbound() */
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "A %d \\nfound better solution: %f\n", (int)nodenum, obj);
      }
      else
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "A %d \\nfound solution: %f\n", (int)nodenum, obj);

      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_SOLUTION);
   }

   if( visual->bakfile != NULL && bettersol )
   {
      SCIP_Real obj;

      if( set->visual_objextern )
         obj = SCIPgetSolOrigObj(set->scip, sol);
      else
         obj = SCIPgetSolTransObj(set->scip, sol);

      if( SCIPsolGetHeur(sol) == NULL )
      {
         /* if LP solution was feasible ... */
         SCIP_VAR* branchvar;
         SCIP_BOUNDTYPE branchtype;
         SCIP_Real branchbound;
         SCIP_NODE *pnode;
         size_t parentnodenum;
         size_t nodenum;
         char t = 'M';

         /* find first parent that is not a probing node */
         assert(node != NULL);
         pnode = node;
         while( pnode != NULL && SCIPnodeGetType(pnode) == SCIP_NODETYPE_PROBINGNODE )
            pnode = pnode->parent;

         if( pnode != NULL )
         {
            /* get node num from hash map */
            nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, pnode);

            /* get nodenum of parent node from hash map */
            parentnodenum = (pnode->parent != NULL ? (size_t)SCIPhashmapGetImage(visual->nodenum, pnode->parent) : 0);
            assert(pnode->parent == NULL || parentnodenum > 0);

            /* get branching information */
            getBranchInfo(pnode, &branchvar, &branchtype, &branchbound);

            /* determine branching type */
            if( branchvar != NULL )
               t = (branchtype == SCIP_BOUNDTYPE_LOWER ? 'R' : 'L');

            printTime(visual, stat, FALSE);
            SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "integer %d %d %c %f\n", (int)nodenum, (int)parentnodenum, t, obj);
         }
      }
      else
      {
         printTime(visual, stat, FALSE);
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "heuristic %f\n", obj);
      }
   }
}

/** outputs a new global lower bound to the visualization output file */
void SCIPvisualLowerbound(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             lowerbound          /**< new lower bound */
   )
{
   assert(visual != NULL);

   /* check, if VBC output should be created */
   if( visual->vbcfile == NULL )
      return;

   /* determine external lower bound */
   if( set->visual_objextern )
      lowerbound = SCIPretransformObj(set->scip, lowerbound);

   printTime(visual, stat, TRUE);
   if( SCIPgetObjsense(set->scip) == SCIP_OBJSENSE_MINIMIZE )
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "L %f\n", lowerbound);
   else
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "U %f\n", lowerbound);

   /* do nothing for BAK */
}

/** outputs a new global upper bound to the visualization output file */
void SCIPvisualUpperbound(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             upperbound          /**< new upper bound */
   )
{
   assert(visual != NULL);

   /* check, if VBC output should be created */
   if( visual->vbcfile == NULL )
      return;

   /* determine external upper bound */
   if( set->visual_objextern )
      upperbound = SCIPretransformObj(set->scip, upperbound);

   printTime(visual, stat, TRUE);
   if( SCIPgetObjsense(set->scip) == SCIP_OBJSENSE_MINIMIZE )
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "U %f\n", upperbound);
   else
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "L %f\n", upperbound);

   /* do nothing for BAK */
}
