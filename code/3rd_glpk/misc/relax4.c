/* relax4.c (relaxation method of Bertsekas and Tseng) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  THIS CODE IS THE RESULT OF TRANSLATION OF THE FORTRAN CODE RELAX4.
*
*  THE TRANSLATION HAS BEEN DONE WITH THE PERMISSION OF THE AUTHOR OF
*  THE ORIGINAL FORTRAN CODE PROF. DIMITRI P. BERTSEKAS, MASSACHUSETTS
*  INSTITUTE OF TECHNOLOGY, CAMBRIDGE, MASSACHUSETTS, USA.
*
*  The translation was made by Andrew Makhorin <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include "relax4.h"

/***********************************************************************
*  WARNING
*
*  A serious bug was *tentatively* fixed in this code (see #if/#endif
*  marked by 'mao').
*
*  This bug is inherited from the original Fortran version of the
*  RELAX-IV code. Unfortunately, the code is very intricate, so this
*  bug is still under investigation. Thanks to Sylvain Fournier for bug
*  report.
*
*  RELAX-IV bug details
*  --------------------
*  In the original RELAX-IV code there are four similar fragments in
*  subroutines ascnt1 and ascnt2 like this:
*
*  C
*  C     DECREASE THE PRICES OF THE SCANNED NODES BY DELPRC.
*  C     ADJUST FLOW TO MAINTAIN COMPLEMENTARY SLACKNESS WITH
*  C     THE PRICES.
*  C
*        NB = 0
*        DO 6 I=1,NSAVE
*        . . .
*           IF (RC(ARC).EQ.0) THEN
*              DELX=DELX+U(ARC)
*              NB = NB + 1
*              PRDCSR(NB) = ARC
*           END IF
*        . . .
*
*  On some instances the variable NB becomes greater than N (the number
*  of nodes) that leads to indexing error, because the array PRDCSR is
*  declared as array of N elements (more precisely, as array of MAXNN
*  elements, however, NB becomes even much greater than MAXNN).
***********************************************************************/

#define false 0
#define true  1

/***********************************************************************
*  NAME
*
*  RELAX-IV (version of October 1994)
*
*  PURPOSE
*
*  This routine implements the relaxation method of Bertsekas and Tseng
*  (see [1], [2]) for linear cost ordinary network flow problems.
*
*  [1] Bertsekas, D. P., "A Unified Framework for Primal-Dual Methods"
*      Mathematical Programming, Vol. 32, 1985, pp. 125-145.
*  [2] Bertsekas, D. P., and Tseng, P., "Relaxation Methods for
*      Minimum Cost" Operations Research, Vol. 26, 1988, pp. 93-114.
*
*  The relaxation method is also described in the books:
*
*  [3] Bertsekas, D. P., "Linear Network Optimization: Algorithms and
*      Codes" MIT Press, 1991.
*  [4] Bertsekas, D. P. and Tsitsiklis, J. N., "Parallel and Distributed
*      Computation: Numerical Methods", Prentice-Hall, 1989.
*  [5] Bertsekas, D. P., "Network Optimization: Continuous and Discrete
*      Models", Athena Scientific, 1998.
*
*  RELEASE NOTE
*
*  This version of relaxation code has option for a special crash
*  procedure for the initial price-flow pair. This is recommended for
*  difficult problems where the default initialization results in long
*  running times. crash = 1 corresponds to an auction/shortest path
*  method
*
*  These initializations are recommended in the absence of any prior
*  information on a favorable initial flow-price vector pair that
*  satisfies complementary slackness.
*
*  The relaxation portion of the code differs from the code RELAXT-III
*  and other earlier relaxation codes in that it maintains the set of
*  nodes with nonzero deficit in a fifo queue. Like its predecessor
*  RELAXT-III, this code maintains a linked list of balanced (i.e., of
*  zero reduced cost) arcs so to reduce the work in labeling and
*  scanning. Unlike RELAXT-III, it does not use selectively shortest
*  path iterations for initialization.
*
*  SOURCE
*
*  The original Fortran code was written by Dimitri P. Bertsekas and
*  Paul Tseng, with a contribution by Jonathan Eckstein in the phase II
*  initialization. The original Fortran routine AUCTION was written by
*  Dimitri P. Bertsekas and is based on the method described in the
*  paper:
*
*  [6] Bertsekas, D. P., "An Auction/Sequential Shortest Path Algorithm
*      for the Minimum Cost Flow Problem", LIDS Report P-2146, MIT,
*      Nov. 1992.
*
*  For inquiries about the original Fortran code, please contact:
*
*  Dimitri P. Bertsekas
*  Laboratory for information and decision systems
*  Massachusetts Institute of Technology
*  Cambridge, MA 02139
*  (617) 253-7267, dimitrib@mit.edu
*
*  This code is the result of translation of the original Fortran code.
*  The translation was made by Andrew Makhorin <mao@gnu.org>.
*
*  USER GUIDELINES
*
*  This routine is in the public domain to be used only for research
*  purposes. It cannot be used as part of a commercial product, or to
*  satisfy in any part commercial delivery requirements to government
*  or industry, without prior agreement with the authors. Users are
*  requested to acknowledge the authorship of the code, and the
*  relaxation method.
*
*  No modification should be made to this code other than the minimal
*  necessary to make it compatible with specific platforms.
*
*  INPUT PARAMETERS (see notes 1, 2, 4)
*
*  n         = number of nodes
*  na        = number of arcs
*  large     = a very large integer to represent infinity
*              (see note 3)
*  repeat    = true if initialization is to be skipped
*              (false otherwise)
*  crash     = 0 if default initialization is used
*              1 if auction initialization is used
*  startn[j] = starting node for arc j, j = 1,...,na
*  endn[j]   = ending node for arc j, j = 1,...,na
*  fou[i]    = first arc out of node i, i = 1,...,n
*  nxtou[j]  = next arc out of the starting node of arc j, j = 1,...,na
*  fin[i]    = first arc into node i, i = 1,...,n
*  nxtin[j]  = next arc into the ending node of arc j, j = 1,...,na
*
*  UPDATED PARAMETERS (see notes 1, 3, 4)
*
*  rc[j]     = reduced cost of arc j, j = 1,...,na
*  u[j]      = capacity of arc j on input
*              and (capacity of arc j) - x[j] on output, j = 1,...,na
*  dfct[i]   = demand at node i on input
*              and zero on output, i = 1,...,n
*
*  OUTPUT PARAMETERS (see notes 1, 3, 4)
*
*  x[j]      = flow on arc j, j = 1,...,na
*  nmultinode = number of multinode relaxation iterations in RELAX4
*  iter      = number of relaxation iterations in RELAX4
*  num_augm  = number of flow augmentation steps in RELAX4
*  num_ascnt = number of multinode ascent steps in RELAX4
*  nsp       = number of auction/shortest path iterations
*
*  WORKING PARAMETERS (see notes 1, 4, 5)
*
*  label[1+n], prdcsr[1+n], save[1+na], tfstou[1+n], tnxtou[1+na],
*  tfstin[1+n], tnxtin[1+na], nxtqueue[1+n], scan[1+n], mark[1+n],
*  extend_arc[1+n], sb_level[1+n], sb_arc[1+n]
*
*  RETURNS
*
*  0         = normal return
*  1,...,8   = problem is found to be infeasible
*
*  NOTE 1
*
*  To run in limited memory systems, declare the arrays startn, endn,
*  nxtin, nxtou, fin, fou, label, prdcsr, save, tfstou, tnxtou, tfstin,
*  tnxtin, ddpos, ddneg, nxtqueue as short instead.
*
*  NOTE 2
*
*  This routine makes no effort to initialize with a favorable x from
*  amongst those flow vectors that satisfy complementary slackness with
*  the initial reduced cost vector rc. If a favorable x is known, then
*  it can be passed, together with the corresponding arrays u and dfct,
*  to this routine directly. This, however, requires that the capacity
*  tightening portion and the flow initialization portion of this
*  routine (up to line labeled 90) be skipped.
*
*  NOTE 3
*
*  All problem data should be less than large in magnitude, and large
*  should be less than, say, 1/4 the largest int of the machine used.
*  This will guard primarily against overflow in uncapacitated problems
*  where the arc capacities are taken finite but very large. Note,
*  however, that as in all codes operating with integers, overflow may
*  occur if some of the problem data takes very large values.
*
*  NOTE 4
*
*  [This note being specific to Fortran was removed.-A.M.]
*
*  NOTE 5
*
*  ddpos and ddneg are arrays that give the directional derivatives for
*  all positive and negative single-node price changes. These are used
*  only in phase II of the initialization procedure, before the linked
*  list of balanced arcs comes to play. Therefore, to reduce storage,
*  they are equivalence to tfstou and tfstin, which are of the same size
*  (number of nodes) and are used only after the tree comes into use. */

static void ascnt1(struct relax4_csa *csa, int dm, int *delx,
      int *nlabel, int *feasbl, int *svitch, int nscan, int curnode,
      int *prevnode);

static void ascnt2(struct relax4_csa *csa, int dm, int *delx,
      int *nlabel, int *feasbl, int *svitch, int nscan, int curnode,
      int *prevnode);

static int auction(struct relax4_csa *csa);

int relax4(struct relax4_csa *csa)
{     /* input parameters */
      int n = csa->n;
      int na = csa->na;
      int large = csa->large;
      int repeat = csa->repeat;
      int crash = csa->crash;
      int *startn = csa->startn;
      int *endn = csa->endn;
      int *fou = csa->fou;
      int *nxtou = csa->nxtou;
      int *fin = csa->fin;
      int *nxtin = csa->nxtin;
      /* updated parameters */
      int *rc = csa->rc;
      int *u = csa->u;
      int *dfct = csa->dfct;
      /* output parameters */
      int *x = csa->x;
#     define nmultinode (csa->nmultinode)
#     define iter (csa->iter)
#     define num_augm (csa->num_augm)
#     define num_ascnt (csa->num_ascnt)
#     define nsp (csa->nsp)
      /* working parameters */
      int *label = csa->label;
      int *prdcsr = csa->prdcsr;
      int *save = csa->save;
      int *tfstou = csa->tfstou;
      int *tnxtou = csa->tnxtou;
      int *tfstin = csa->tfstin;
      int *tnxtin = csa->tnxtin;
      int *nxtqueue = csa->nxtqueue;
      char *scan = csa->scan;
      char *mark = csa->mark;
      int *ddpos = tfstou;
      int *ddneg = tfstin;
      /* local variables */
      int arc, augnod, capin, capout, defcit, delprc, delx, dm, dp,
         dx, feasbl, i, ib, indef, j, lastqueue, maxcap, narc, nb,
         nlabel, node, node2, node_def, naugnod, nscan, num_passes,
         numnz, numnz_new, numpasses, nxtarc, nxtbrk, nxtnode, passes,
         pchange, posit, prevnode, prvarc, quit, rdcost, scapin,
         scapou, svitch, t, t1, t2, tmparc, tp, trc, ts;
      /*--------------------------------------------------------------*/
      /* Initialization phase I */
      /* In this phase, we reduce the arc capacities by as much as
       * possible without changing the problem; then we set the initial
       * flow array x, together with the corresponding arrays u and
       * dfct. */
      /* This phase and phase II (from here up to line labeled 90) can
       * be skipped (by setting repeat to true) if the calling program
       * places in common user-chosen values for the arc flows, the
       * residual arc capacities, and the nodal deficits. When this is
       * done, it is critical that the flow and the reduced cost for
       * each arc satisfy complementary slackness and the dfct array
       * properly correspond to the initial arc/flows. */
      if (repeat)
         goto L90;
      for (node = 1; node <= n; node++)
      {  node_def = dfct[node];
         ddpos[node] = node_def;
         ddneg[node] = -node_def;
         maxcap = 0;
         scapou = 0;
         for (arc = fou[node]; arc > 0; arc = nxtou[arc])
         {  if (scapou <= large - u[arc])
               scapou += u[arc];
            else
               goto L10;
         }
         if (scapou <= large - node_def)
            capout = scapou + node_def;
         else
            goto L10;
         if (capout < 0)
         {  /* problem is infeasible */
            /* exogenous flow into node exceeds out capacity */
            return 1;
         }
         scapin = 0;
         for (arc = fin[node]; arc > 0; arc = nxtin[arc])
         {  if (u[arc] > capout)
               u[arc] = capout;
            if (maxcap < u[arc])
               maxcap = u[arc];
            if (scapin <= large - u[arc])
               scapin += u[arc];
            else
               goto L10;
         }
         if (scapin <= large + node_def)
            capin = scapin - node_def;
         else
            goto L10;
         if (capin < 0)
         {  /* problem is infeasible */
            /* exogenous flow out of node exceeds in capacity */
            return 2;
         }
         for (arc = fou[node]; arc > 0; arc = nxtou[arc])
         {  if (u[arc] > capin)
               u[arc] = capin;
         }
L10:     ;
      }
      /*--------------------------------------------------------------*/
      /* Initialization phase II */
      /* In this phase, we initialize the prices and flows by either
       * calling the routine auction or by performing only single node
       * (coordinate) relaxation iterations. */
      if (crash == 1)
      {  nsp = 0;
         if (auction(csa) != 0)
         {  /* problem is found to be infeasible */
            return 3;
         }
         goto L70;
      }
      /* Initialize the arc flows to satisfy complementary slackness
       * with the prices. u[arc] is the residual capacity of arc, and
       * x[arc] is the flow. These two always add up to the total
       * capacity for arc. Also compute the directional derivatives for
       * each coordinate and compute the actual deficits. */
      for (arc = 1; arc <= na; arc++)
      {  x[arc] = 0;
         if (rc[arc] <= 0)
         {  t = u[arc];
            t1 = startn[arc];
            t2 = endn[arc];
            ddpos[t1] += t;
            ddneg[t2] += t;
            if (rc[arc] < 0)
            {  x[arc] = t;
               u[arc] = 0;
               dfct[t1] += t;
               dfct[t2] -= t;
               ddneg[t1] -= t;
               ddpos[t2] -= t;
            }
         }
      }
      /* Make 2 or 3 passes through all nodes, performing only single
       * node relaxation iterations. The number of passes depends on the
       * density of the network. */
      if (na > n * 10)
         numpasses = 2;
      else
         numpasses = 3;
      for (passes = 1; passes <= numpasses; passes++)
      for (node = 1; node <= n; node++)
      {  if (dfct[node] == 0)
            continue;
         if (ddpos[node] <= 0)
         {  /* Compute delprc, the stepsize to the next breakpoint in
             * the dual cost as the price of node is increased.
             * [Since the reduced cost of all outgoing (resp., incoming)
             * arcs will decrease (resp., increase) as the price of node
             * is increased, the next breakpoint is the minimum of the
             * positive reduced cost on outgoing arcs and of the
             * negative reduced cost on incoming arcs.] */
            delprc = large;
            for (arc = fou[node]; arc > 0; arc = nxtou[arc])
            {  trc = rc[arc];
               if ((trc > 0) && (trc < delprc))
                  delprc = trc;
            }
            for (arc = fin[node]; arc > 0; arc = nxtin[arc])
            {  trc = rc[arc];
               if ((trc < 0) && (trc > -delprc))
                  delprc = -trc;
            }
            /* If no breakpoint is left and dual ascent is still
             * possible, the problem is infeasible. */
            if (delprc >= large)
            {  if (ddpos[node] == 0)
                  continue;
               return 4;
            }
            /* delprc is the stepsize to next breakpoint. Increase
             * price of node by delprc and compute the stepsize to the
             * next breakpoint in the dual cost. */
L53:        nxtbrk = large;
            /* Look at all arcs out of node. */
            for (arc = fou[node]; arc > 0; arc = nxtou[arc])
            {  trc = rc[arc];
               if (trc == 0)
               {  t1 = endn[arc];
                  t = u[arc];
                  if (t > 0)
                  {  dfct[node] += t;
                     dfct[t1] -= t;
                     x[arc] = t;
                     u[arc] = 0;
                  }
                  else
                     t = x[arc];
                  ddneg[node] -= t;
                  ddpos[t1] -= t;
               }
               /* Decrease the reduced costs on all outgoing arcs. */
               trc -= delprc;
               if ((trc > 0) && (trc < nxtbrk))
                  nxtbrk = trc;
               else if (trc == 0)
               {  /* Arc goes from inactive to balanced. Update the rate
                   * of dual ascent at node and at its neighbor. */
                  ddpos[node] += u[arc];
                  ddneg[endn[arc]] += u[arc];
               }
               rc[arc] = trc;
            }
            /* Look at all arcs into node. */
            for (arc = fin[node]; arc > 0; arc = nxtin[arc])
            {  trc = rc[arc];
               if (trc == 0)
               {  t1 = startn[arc];
                  t = x[arc];
                  if (t > 0)
                  {  dfct[node] += t;
                     dfct[t1] -= t;
                     u[arc] = t;
                     x[arc] = 0;
                  }
                  else
                     t = u[arc];
                  ddpos[t1] -= t;
                  ddneg[node] -= t;
               }
               /* Increase the reduced cost on all incoming arcs. */
               trc += delprc;
               if ((trc < 0) && (trc > -nxtbrk))
                  nxtbrk = -trc;
               else if (trc == 0)
               {  /* Arc goes from active to balanced. Update the rate
                   * of dual ascent at node and at its neighbor. */
                  ddneg[startn[arc]] += x[arc];
                  ddpos[node] += x[arc];
               }
               rc[arc] = trc;
            }
            /* If price of node can be increased further without
             * decreasing the dual cost (even the dual cost doesn't
             * increase), return to increase the price further. */
            if ((ddpos[node] <= 0) && (nxtbrk < large))
            {  delprc = nxtbrk;
               goto L53;
            }
         }
         else if (ddneg[node] <= 0)
         {  /* Compute delprc, the stepsize to the next breakpoint in
             * the dual cost as the price of node is decreased.
             * [Since the reduced cost of all outgoing (resp., incoming)
             * arcs will increase (resp., decrease) as the price of node
             * is decreased, the next breakpoint is the minimum of the
             * negative reduced cost on outgoing arcs and of the
             * positive reduced cost on incoming arcs.] */
            delprc = large;
            for (arc = fou[node]; arc > 0; arc = nxtou[arc])
            {  trc = rc[arc];
               if ((trc < 0) && (trc > -delprc))
                  delprc = -trc;
            }
            for (arc = fin[node]; arc > 0; arc = nxtin[arc])
            {  trc = rc[arc];
               if ((trc > 0) && (trc < delprc))
                  delprc = trc;
            }
            /* If no breakpoint is left and dual ascent is still
             * possible, the problem is infeasible. */
            if (delprc == large)
            {  if (ddneg[node] == 0)
                  continue;
               return 5;
            }
            /* delprc is the stepsize to next breakpoint. Decrease
             * price of node by delprc and compute the stepsize to the
             * next breakpoint in the dual cost. */
L63:        nxtbrk = large;
            /* Look at all arcs out of node. */
            for (arc = fou[node]; arc > 0; arc = nxtou[arc])
            {  trc = rc[arc];
               if (trc == 0)
               {  t1 = endn[arc];
                  t = x[arc];
                  if (t > 0)
                  {  dfct[node] -= t;
                     dfct[t1] += t;
                     u[arc] = t;
                     x[arc] = 0;
                  }
                  else
                     t = u[arc];
                  ddpos[node] -= t;
                  ddneg[t1] -= t;
               }
               /* Increase the reduced cost on all outgoing arcs. */
               trc += delprc;
               if ((trc < 0) && (trc > -nxtbrk))
                  nxtbrk = -trc;
               else if (trc == 0)
               {  /* Arc goes from active to balanced. Update the rate
                   * of dual ascent at node and at its neighbor. */
                  ddneg[node] += x[arc];
                  ddpos[endn[arc]] += x[arc];
               }
               rc[arc] = trc;
            }
            /* Look at all arcs into node. */
            for (arc = fin[node]; arc > 0; arc = nxtin[arc])
            {  trc = rc[arc];
               if (trc == 0)
               {  t1 = startn[arc];
                  t = u[arc];
                  if (t > 0)
                  {  dfct[node] -= t;
                     dfct[t1] += t;
                     x[arc] = t;
                     u[arc] = 0;
                  }
                  else
                     t = x[arc];
                  ddneg[t1] -= t;
                  ddpos[node] -= t;
               }
               /* Decrease the reduced cost on all incoming arcs. */
               trc -= delprc;
               if ((trc > 0) && (trc < nxtbrk))
                  nxtbrk = trc;
               else if (trc == 0)
               {  /* Arc goes from inactive to balanced. Update the rate
                   * of dual ascent at node and at its neighbor. */
                  ddpos[startn[arc]] += u[arc];
                  ddneg[node] += u[arc];
               }
               rc[arc] = trc;
            }
            /* If price of node can be decreased further without
             * decreasing the dual cost (even the dual cost doesn't
             * increase), return to decrease the price further. */
            if ((ddneg[node] <= 0) && (nxtbrk < large))
            {  delprc = nxtbrk;
               goto L63;
            }
         }
      }
      /*--------------------------------------------------------------*/
L70:  /* Initialize tree data structure. */
      for (i = 1; i <= n; i++)
         tfstou[i] = tfstin[i] = 0;
      for (i = 1; i <= na; i++)
      {  tnxtin[i] = tnxtou[i] = -1;
         if (rc[i] == 0)
         {  tnxtou[i] = tfstou[startn[i]];
            tfstou[startn[i]] = i;
            tnxtin[i] = tfstin[endn[i]];
            tfstin[endn[i]] = i;
         }
      }
L90:  /* Initialize other variables. */
      feasbl = true;
      iter = 0;
      nmultinode = 0;
      num_augm = 0;
      num_ascnt = 0;
      num_passes = 0;
      numnz = n;
      numnz_new = 0;
      svitch = false;
      for (i = 1; i <= n; i++)
         mark[i] = scan[i] = false;
      nlabel = 0;
      /* RELAX4 uses an adaptive strategy to decide whether to continue
       * the scanning process after a multinode price change.
       * The threshold parameter tp and ts that control this strategy
       * are set in the next two lines. */
      tp = 10;
      ts = n / 15;
      /* Initialize the queue of nodes with nonzero deficit. */
      for (node = 1; node <= n - 1; node++)
         nxtqueue[node] = node + 1;
      nxtqueue[n] = 1;
      node = lastqueue = n;
      /*--------------------------------------------------------------*/
      /* Start the relaxation algorithm. */
L100: /* Code for advancing the queue of nonzero deficit nodes. */
      prevnode = node;
      node = nxtqueue[node];
      defcit = dfct[node];
      if (node == lastqueue)
      {  numnz = numnz_new;
         numnz_new = 0;
         lastqueue = prevnode;
         num_passes++;
      }
      /* Code for deleting a node from the queue. */
      if (defcit == 0)
      {  nxtnode = nxtqueue[node];
         if (node == nxtnode)
            return 0;
         else
         {  nxtqueue[prevnode] = nxtnode;
            nxtqueue[node] = 0;
            node = nxtnode;
            goto L100;
         }
      }
      else
         posit = (defcit > 0);
      iter++;
      numnz_new++;
      if (posit)
      {  /* Attempt a single node iteration from node with positive
          * deficit. */
         pchange = false;
         indef = defcit;
         delx = 0;
         nb = 0;
         /* Check outgoing (probably) balanced arcs from node. */
         for (arc = tfstou[node]; arc > 0; arc = tnxtou[arc])
         {  if ((rc[arc] == 0) && (x[arc] > 0))
            {  delx += x[arc];
               nb++;
               save[nb] = arc;
            }
         }
         /* Check incoming arcs. */
         for (arc = tfstin[node]; arc > 0; arc = tnxtin[arc])
         {  if ((rc[arc] == 0) && (u[arc] > 0))
            {  delx += u[arc];
               nb++;
               save[nb] = -arc;
            }
         }
         /* End of initial node scan. */
L4018:   /* If no price change is possible, exit. */
         if (delx > defcit)
         {  quit = (defcit < indef);
            goto L4016;
         }
         /* RELAX4 searches along the ascent direction for the best
          * price by checking the slope of the dual cost at successive
          * break points. First, we compute the distance to the next
          * break point. */
         delprc = large;
         for (arc = fou[node]; arc > 0; arc = nxtou[arc])
         {  rdcost = rc[arc];
            if ((rdcost < 0) && (rdcost > -delprc))
               delprc = -rdcost;
         }
         for (arc = fin[node]; arc > 0; arc = nxtin[arc])
         {  rdcost = rc[arc];
            if ((rdcost > 0) && (rdcost < delprc))
               delprc = rdcost;
         }
         /* Check if problem is infeasible. */
         if ((delx < defcit) && (delprc == large))
         {  /* The dual cost can be decreased without bound. */
            return 6;
         }
         /* Skip flow adjustment if there is no flow to modify. */
         if (delx == 0)
            goto L4014;
         /* Adjust the flow on the balanced arcs incident to node to
          * maintain complementary slackness after the price change. */
         for (j = 1; j <= nb; j++)
         {  arc = save[j];
            if (arc > 0)
            {  node2 = endn[arc];
               t1 = x[arc];
               dfct[node2] += t1;
               if (nxtqueue[node2] == 0)
               {  nxtqueue[prevnode] = node2;
                  nxtqueue[node2] = node;
                  prevnode = node2;
               }
               u[arc] += t1;
               x[arc] = 0;
            }
            else
            {  narc = -arc;
               node2 = startn[narc];
               t1 = u[narc];
               dfct[node2] += t1;
               if (nxtqueue[node2] == 0)
               {  nxtqueue[prevnode] = node2;
                  nxtqueue[node2] = node;
                  prevnode = node2;
               }
               x[narc] += t1;
               u[narc] = 0;
            }
         }
         defcit -= delx;
L4014:   if (delprc == large)
         {  quit = true;
            goto L4019;
         }
         /* Node corresponds to a dual ascent direction. Decrease the
          * price of node by delprc and compute the stepsize to the next
          * breakpoint in the dual cost. */
         nb = 0;
         pchange = true;
         dp = delprc;
         delprc = large;
         delx = 0;
         for (arc = fou[node]; arc > 0; arc = nxtou[arc])
         {  rdcost = rc[arc] + dp;
            rc[arc] = rdcost;
            if (rdcost == 0)
            {  nb++;
               save[nb] = arc;
               delx += x[arc];
            }
            if ((rdcost < 0) && (rdcost > -delprc))
               delprc = -rdcost;
         }
         for (arc = fin[node]; arc > 0; arc = nxtin[arc])
         {  rdcost = rc[arc] - dp;
            rc[arc] = rdcost;
            if (rdcost == 0)
            {  nb++;
               save[nb] = -arc;
               delx += u[arc];
            }
            if ((rdcost > 0) && (rdcost < delprc))
               delprc = rdcost;
         }
         /* Return to check if another price change is possible. */
         goto L4018;
L4016:   /* Perform flow augmentation at node. */
         for (j = 1; j <= nb; j++)
         {  arc = save[j];
            if (arc > 0)
            {  /* arc is an outgoing arc from node. */
               node2 = endn[arc];
               t1 = dfct[node2];
               if (t1 < 0)
               {  /* Decrease the total deficit by decreasing flow of
                   * arc. */
                  quit = true;
                  t2 = x[arc];
                  dx = defcit;
                  if (dx > -t1) dx = -t1;
                  if (dx > t2) dx = t2;
                  defcit -= dx;
                  dfct[node2] = t1 + dx;
                  if (nxtqueue[node2] == 0)
                  {  nxtqueue[prevnode] = node2;
                     nxtqueue[node2] = node;
                     prevnode = node2;
                  }
                  x[arc] = t2 - dx;
                  u[arc] += dx;
                  if (defcit == 0)
                     break;
               }
            }
            else
            {  /* -arc is an incoming arc to node. */
               narc = -arc;
               node2 = startn[narc];
               t1 = dfct[node2];
               if (t1 < 0)
               {  /* Decrease the total deficit by increasing flow of
                   * -arc. */
                  quit = true;
                  t2 = u[narc];
                  dx = defcit;
                  if (dx > -t1) dx = -t1;
                  if (dx > t2) dx = t2;
                  defcit -= dx;
                  dfct[node2] = t1 + dx;
                  if (nxtqueue[node2] == 0)
                  {  nxtqueue[prevnode] = node2;
                     nxtqueue[node2] = node;
                     prevnode = node2;
                  }
                  x[narc] += dx;
                  u[narc] = t2 - dx;
                  if (defcit == 0)
                     break;
               }
            }
         }
L4019:   dfct[node] = defcit;
         /* Reconstruct the linked list of balance arcs incident to this
          * node. For each adjacent node, we add any newly balanced arcs
          * to the list, but do not bother removing formerly balanced
          * ones (they will be removed the next time each adjacent node
          * is scanned). */
         if (pchange)
         {  arc = tfstou[node];
            tfstou[node] = 0;
            while (arc > 0)
            {  nxtarc = tnxtou[arc];
               tnxtou[arc] = -1;
               arc = nxtarc;
            }
            arc = tfstin[node];
            tfstin[node] = 0;
            while (arc > 0)
            {  nxtarc = tnxtin[arc];
               tnxtin[arc] = -1;
               arc = nxtarc;
            }
            /* Now add the currently balanced arcs to the list for this
             * node (which is now empty), and the appropriate adjacent
             * ones. */
            for (j = 1; j <= nb; j++)
            {  arc = save[j];
               if (arc < 0)
                  arc = -arc;
               if (tnxtou[arc] < 0)
               {  tnxtou[arc] = tfstou[startn[arc]];
                  tfstou[startn[arc]] = arc;
               }
               if (tnxtin[arc] < 0)
               {  tnxtin[arc] = tfstin[endn[arc]];
                  tfstin[endn[arc]] = arc;
               }
            }
         }
         /* End of single node iteration for positive deficit node. */
      }
      else
      {  /* Attempt a single node iteration from node with negative
          * deficit. */
         pchange = false;
         defcit = -defcit;
         indef = defcit;
         delx = 0;
         nb = 0;
         for (arc = tfstin[node]; arc > 0; arc = tnxtin[arc])
         {  if ((rc[arc] == 0) && (x[arc] > 0))
            {  delx += x[arc];
               nb++;
               save[nb] = arc;
            }
         }
         for (arc = tfstou[node]; arc > 0; arc = tnxtou[arc])
         {  if ((rc[arc] == 0) && (u[arc] > 0))
            {  delx += u[arc];
               nb++;
               save[nb] = -arc;
            }
         }
L4028:   if (delx >= defcit)
         {  quit = (defcit < indef);
            goto L4026;
         }
         /* Compute distance to next breakpoint. */
         delprc = large;
         for (arc = fin[node]; arc > 0; arc = nxtin[arc])
         {  rdcost = rc[arc];
            if ((rdcost < 0) && (rdcost > -delprc))
               delprc = -rdcost;
         }
         for (arc = fou[node]; arc > 0; arc = nxtou[arc])
         {  rdcost = rc[arc];
            if ((rdcost > 0) && (rdcost < delprc))
               delprc = rdcost;
         }
         /* Check if problem is infeasible. */
         if ((delx < defcit) && (delprc == large))
            return 7;
         if (delx == 0)
            goto L4024;
         /* Flow augmentation is possible. */
         for (j = 1; j <= nb; j++)
         {  arc = save[j];
            if (arc > 0)
            {  node2 = startn[arc];
               t1 = x[arc];
               dfct[node2] -= t1;
               if (nxtqueue[node2] == 0)
               {  nxtqueue[prevnode] = node2;
                  nxtqueue[node2] = node;
                  prevnode = node2;
               }
               u[arc] += t1;
               x[arc] = 0;
            }
            else
            {  narc = -arc;
               node2 = endn[narc];
               t1 = u[narc];
               dfct[node2] -= t1;
               if (nxtqueue[node2] == 0)
               {  nxtqueue[prevnode] = node2;
                  nxtqueue[node2] = node;
                  prevnode = node2;
               }
               x[narc] += t1;
               u[narc] = 0;
            }
         }
         defcit -= delx;
L4024:   if (delprc == large)
         {  quit = true;
            goto L4029;
         }
         /* Price increase at node is possible. */
         nb = 0;
         pchange = true;
         dp = delprc;
         delprc = large;
         delx = 0;
         for (arc = fin[node]; arc > 0; arc = nxtin[arc])
         {  rdcost = rc[arc] + dp;
            rc[arc] = rdcost;
            if (rdcost == 0)
            {  nb++;
               save[nb] = arc;
               delx += x[arc];
            }
            if ((rdcost < 0) && (rdcost > -delprc))
               delprc = -rdcost;
         }
         for (arc = fou[node]; arc > 0; arc = nxtou[arc])
         {  rdcost = rc[arc] - dp;
            rc[arc] = rdcost;
            if (rdcost == 0)
            {  nb++;
               save[nb] = -arc;
               delx += u[arc];
            }
            if ((rdcost > 0) && (rdcost < delprc))
               delprc = rdcost;
         }
         goto L4028;
L4026:   /* Perform flow augmentation at node. */
         for (j = 1; j <= nb; j++)
         {  arc = save[j];
            if (arc > 0)
            {  /* arc is an incoming arc to node. */
               node2 = startn[arc];
               t1 = dfct[node2];
               if (t1 > 0)
               {  quit = true;
                  t2 = x[arc];
                  dx = defcit;
                  if (dx > t1) dx = t1;
                  if (dx > t2) dx = t2;
                  defcit -= dx;
                  dfct[node2] = t1 - dx;
                  if (nxtqueue[node2] == 0)
                  {  nxtqueue[prevnode] = node2;
                     nxtqueue[node2] = node;
                     prevnode = node2;
                  }
                  x[arc] = t2 - dx;
                  u[arc] += dx;
                  if (defcit == 0)
                     break;
               }
            }
            else
            {  /* -arc is an outgoing arc from node. */
               narc = -arc;
               node2 = endn[narc];
               t1 = dfct[node2];
               if (t1 > 0)
               {  quit = true;
                  t2 = u[narc];
                  dx = defcit;
                  if (dx > t1) dx = t1;
                  if (dx > t2) dx = t2;
                  defcit -= dx;
                  dfct[node2] = t1 - dx;
                  if (nxtqueue[node2] == 0)
                  {  nxtqueue[prevnode] = node2;
                     nxtqueue[node2] = node;
                     prevnode = node2;
                  }
                  x[narc] += dx;
                  u[narc] = t2 - dx;
                  if (defcit == 0)
                     break;
               }
            }
         }
L4029:   dfct[node] = -defcit;
         /* Reconstruct the list of balanced arcs incident to node. */
         if (pchange)
         {  arc = tfstou[node];
            tfstou[node] = 0;
            while (arc > 0)
            {  nxtarc = tnxtou[arc];
               tnxtou[arc] = -1;
               arc = nxtarc;
            }
            arc = tfstin[node];
            tfstin[node] = 0;
            while (arc > 0)
            {  nxtarc = tnxtin[arc];
               tnxtin[arc] = -1;
               arc = nxtarc;
            }
            /* Now add the currently balanced arcs to the list for this
             * node (which is now empty), and the appropriate adjacent
             * ones. */
            for (j = 1; j <= nb; j++)
            {  arc = save[j];
               if (arc <= 0)
                  arc = -arc;
               if (tnxtou[arc] < 0)
               {  tnxtou[arc] = tfstou[startn[arc]];
                  tfstou[startn[arc]] = arc;
               }
               if (tnxtin[arc] < 0)
               {  tnxtin[arc] = tfstin[endn[arc]];
                  tfstin[endn[arc]] = arc;
               }
            }
         }
         /* End of single node iteration for a negative deficit node. */
      }
      if (quit || (num_passes <= 3))
         goto L100;
      /* Do a multinode iteration from node. */
      nmultinode++;
      /* If number of nonzero deficit nodes is small, continue labeling
       * until a flow augmentation is done. */
      svitch = (numnz < tp);
      /* Unmark nodes labeled earlier. */
      for (j = 1; j <= nlabel; j++)
      {  node2 = label[j];
         mark[node2] = scan[node2] = false;
      }
      /* Initialize labeling. */
      nlabel = 1;
      label[1] = node;
      mark[node] = true;
      prdcsr[node] = 0;
      /* Scan starting node. */
      scan[node] = true;
      nscan = 1;
      dm = dfct[node];
      delx = 0;
      for (j = 1; j <= nb; j++)
      {  arc = save[j];
         if (arc > 0)
         {  if (posit)
               node2 = endn[arc];
            else
               node2 = startn[arc];
            if (!mark[node2])
            {  nlabel++;
               label[nlabel] = node2;
               prdcsr[node2] = arc;
               mark[node2] = true;
               delx += x[arc];
            }
         }
         else
         {  narc = -arc;
            if (posit)
               node2 = startn[narc];
            else
               node2 = endn[narc];
            if (!mark[node2])
            {  nlabel++;
               label[nlabel] = node2;
               prdcsr[node2] = arc;
               mark[node2] = true;
               delx += u[narc];
            }
         }
      }
L4120:/* Start scanning a labeled but unscanned node. */
      nscan++;
      /* Check to see if switch needs to be set to true so to continue
       * scanning even after a price change. */
      svitch = svitch || ((nscan > ts) && (numnz < ts));
      /* Scanning will continue until either an overestimate of the
       * residual capacity across the cut corresponding to the scanned
       * set of nodes (called delx) exceeds the absolute value of the
       * total deficit of the scanned nodes (called dm), or else an
       * augmenting path is found. Arcs that are in the tree but are not
       * balanced are removed as part of the scanning process. */
      i = label[nscan];
      scan[i] = true;
      naugnod = 0;
      if (posit)
      {  /* Scanning node i in case of positive deficit. */
         prvarc = 0;
         arc = tfstou[i];
         while (arc > 0)
         {  /* arc is an outgoing arc from node. */
            if (rc[arc] == 0)
            {  if (x[arc] > 0)
               {  node2 = endn[arc];
                  if (!mark[node2])
                  {  /* node2 is not labeled, so add node2 to the
                        labeled set. */
                     prdcsr[node2] = arc;
                     if (dfct[node2] < 0)
                     {  naugnod++;
                        save[naugnod] = node2;
                     }
                     nlabel++;
                     label[nlabel] = node2;
                     mark[node2] = true;
                     delx += x[arc];
                  }
               }
               prvarc = arc;
               arc = tnxtou[arc];
            }
            else
            {  tmparc = arc;
               arc = tnxtou[arc];
               tnxtou[tmparc] = -1;
               if (prvarc == 0)
                  tfstou[i] = arc;
               else
                  tnxtou[prvarc] = arc;
            }
         }
         prvarc = 0;
         arc = tfstin[i];
         while (arc > 0)
         {  /* arc is an incoming arc into node. */
            if (rc[arc] == 0)
            {  if (u[arc] > 0)
               {  node2 = startn[arc];
                  if (!mark[node2])
                  {  /* node2 is not labeled, so add node2 to the
                      * labeled set. */
                     prdcsr[node2] = -arc;
                     if (dfct[node2] < 0)
                     {  naugnod++;
                        save[naugnod] = node2;
                     }
                     nlabel++;
                     label[nlabel] = node2;
                     mark[node2] = true;
                     delx += u[arc];
                  }
               }
               prvarc = arc;
               arc = tnxtin[arc];
            }
            else
            {  tmparc = arc;
               arc = tnxtin[arc];
               tnxtin[tmparc] = -1;
               if (prvarc == 0)
                  tfstin[i] = arc;
               else
                  tnxtin[prvarc] = arc;
            }
         }
         /* Correct the residual capacity of the scanned node cut. */
         arc = prdcsr[i];
         if (arc > 0)
            delx -= x[arc];
         else
            delx -= u[-arc];
         /* End of scanning of node i for positive deficit case. */
      }
      else
      {  /* Scanning node i for negative deficit case. */
         prvarc = 0;
         arc = tfstin[i];
         while (arc > 0)
         {  if (rc[arc] == 0)
            {  if (x[arc] > 0)
               {  node2 = startn[arc];
                  if (!mark[node2])
                  {  prdcsr[node2] = arc;
                     if (dfct[node2] > 0)
                     {  naugnod++;
                        save[naugnod] = node2;
                     }
                     nlabel++;
                     label[nlabel] = node2;
                     mark[node2] = true;
                     delx += x[arc];
                  }
               }
               prvarc = arc;
               arc = tnxtin[arc];
            }
            else
            {  tmparc = arc;
               arc = tnxtin[arc];
               tnxtin[tmparc] = -1;
               if (prvarc == 0)
                  tfstin[i] = arc;
               else
                  tnxtin[prvarc] = arc;
            }
         }
         prvarc = 0;
         arc = tfstou[i];
         while (arc > 0)
         {  if (rc[arc] == 0)
            {  if (u[arc] > 0)
               {  node2 = endn[arc];
                  if (!mark[node2])
                  {  prdcsr[node2] = -arc;
                     if (dfct[node2] > 0)
                     {  naugnod++;
                        save[naugnod] = node2;
                     }
                     nlabel++;
                     label[nlabel] = node2;
                     mark[node2] = true;
                     delx += u[arc];
                  }
               }
               prvarc = arc;
               arc = tnxtou[arc];
            }
            else
            {  tmparc = arc;
               arc = tnxtou[arc];
               tnxtou[tmparc] = -1;
               if (prvarc == 0)
                  tfstou[i] = arc;
               else
                  tnxtou[prvarc] = arc;
            }
         }
         arc = prdcsr[i];
         if (arc > 0)
            delx -= x[arc];
         else
            delx -= u[-arc];
      }
      /* Add deficit of node scanned to dm. */
      dm += dfct[i];
      /* Check if the set of scanned nodes correspond to a dual ascent
       * direction; if yes, perform a price adjustment step, otherwise
       * continue labeling. */
      if (nscan < nlabel)
      {  if (svitch)
            goto L4210;
         if ((delx >= dm) && (delx >= -dm))
            goto L4210;
      }
      /* Try a price change.
       * [Note that since delx - abs(dm) is an overestimate of ascent
       * slope, we may occasionally try a direction that is not an
       * ascent direction. In this case the ascnt routines return with
       * quit = false, so we continue labeling nodes.] */
      if (posit)
      {  ascnt1(csa, dm, &delx, &nlabel, &feasbl, &svitch, nscan, node,
            &prevnode);
         num_ascnt++;
      }
      else
      {  ascnt2(csa, dm, &delx, &nlabel, &feasbl, &svitch, nscan, node,
            &prevnode);
         num_ascnt++;
      }
      if (!feasbl)
         return 8;
      if (!svitch)
         goto L100;
      /* Store those newly labeled nodes to which flow augmentation is
       * possible. */
      naugnod = 0;
      for (j = nscan + 1; j <= nlabel; j++)
      {  node2 = label[j];
         if (posit && (dfct[node2] < 0))
         {  naugnod++;
            save[naugnod] = node2;
         }
         else if ((!posit) && (dfct[node2] > 0))
         {  naugnod++;
            save[naugnod] = node2;
         }
      }
L4210:/* Check if flow augmentation is possible. If not, return to scan
       * another node. */
      if (naugnod == 0)
         goto L4120;
      for (j = 1; j <= naugnod; j++)
      {  num_augm++;
         augnod = save[j];
         if (posit)
         {  /* Do the augmentation from node with positive deficit. */
            dx = -dfct[augnod];
            ib = augnod;
            while (ib != node)
            {  arc = prdcsr[ib];
               if (arc > 0)
               {  if (dx > x[arc]) dx = x[arc];
                  ib = startn[arc];
               }
               else
               {  if (dx > u[-arc]) dx = u[-arc];
                  ib = endn[-arc];
               }
            }
            if (dx > dfct[node]) dx = dfct[node];
            if (dx > 0)
            {  /* Increase (decrease) the flow of all forward (backward)
                * arcs in the flow augmenting path. Adjust node deficit
                * accordingly. */
               if (nxtqueue[augnod] == 0)
               {  nxtqueue[prevnode] = augnod;
                  nxtqueue[augnod] = node;
                  prevnode = augnod;
               }
               dfct[augnod] += dx;
               dfct[node] -= dx;
               ib = augnod;
               while (ib != node)
               {  arc = prdcsr[ib];
                  if (arc > 0)
                  {  x[arc] -= dx;
                     u[arc] += dx;
                     ib = startn[arc];
                  }
                  else
                  {  narc = -arc;
                     x[narc] += dx;
                     u[narc] -= dx;
                     ib = endn[narc];
                  }
               }
            }
         }
         else
         {  /* Do the augmentation from node with negative deficit. */
            dx = dfct[augnod];
            ib = augnod;
            while (ib != node)
            {  arc = prdcsr[ib];
               if (arc > 0)
               {  if (dx > x[arc]) dx = x[arc];
                  ib = endn[arc];
               }
               else
               {  if (dx > u[-arc]) dx = u[-arc];
                  ib = startn[-arc];
               }
            }
            if (dx > -dfct[node]) dx = -dfct[node];
            if (dx > 0)
            {  /* Update the flow and deficits. */
               if (nxtqueue[augnod] == 0)
               {  nxtqueue[prevnode] = augnod;
                  nxtqueue[augnod] = node;
                  prevnode = augnod;
               }
               dfct[augnod] -= dx;
               dfct[node] += dx;
               ib = augnod;
               while (ib != node)
               {  arc = prdcsr[ib];
                  if (arc > 0)
                  {  x[arc] -= dx;
                     u[arc] += dx;
                     ib = endn[arc];
                  }
                  else
                  {  narc = -arc;
                     x[narc] += dx;
                     u[narc] -= dx;
                     ib = startn[narc];
                  }
               }
            }
         }
         if (dfct[node] == 0)
            goto L100;
         if (dfct[augnod] != 0)
            svitch = false;
      }
      /* If node still has nonzero deficit and all newly labeled nodes
       * have same sign for their deficit as node, we can continue
       * labeling. In this case, continue labeling only when flow
       * augmentation is done relatively infrequently. */
      if (svitch && (iter > 8 * num_augm))
         goto L4120;
      /* Return to do another relaxation iteration. */
      goto L100;
#     undef nmultinode
#     undef iter
#     undef num_augm
#     undef num_ascnt
#     undef nsp
}

/***********************************************************************
*  NAME
*
*  relax4_inidat - construct linked lists for network topology
*
*  PURPOSE
*
*  This routine constructs two linked lists for the network topology:
*  one list (given by fou, nxtou) for the outgoing arcs of nodes and
*  one list (given by fin, nxtin) for the incoming arcs of nodes. These
*  two lists are required by RELAX4.
*
*  INPUT PARAMETERS
*
*  n         = number of nodes
*  na        = number of arcs
*  startn[j] = starting node for arc j, j = 1,...,na
*  endn[j]   = ending node for arc j, j = 1,...,na
*
*  OUTPUT PARAMETERS
*
*  fou[i]    = first arc out of node i, i = 1,...,n
*  nxtou[j]  = next arc out of the starting node of arc j, j = 1,...,na
*  fin[i]    = first arc into node i, i = 1,...,n
*  nxtin[j]  = next arc into the ending node of arc j, j = 1,...,na
*
*  WORKING PARAMETERS
*
*  tempin[1+n], tempou[1+n] */

void relax4_inidat(struct relax4_csa *csa)
{     /* input parameters */
      int n = csa->n;
      int na = csa->na;
      int *startn = csa->startn;
      int *endn = csa->endn;
      /* output parameters */
      int *fou = csa->fou;
      int *nxtou = csa->nxtou;
      int *fin = csa->fin;
      int *nxtin = csa->nxtin;
      /* working parameters */
      int *tempin = csa->label;
      int *tempou = csa->prdcsr;
      /* local variables */
      int i, i1, i2;
      for (i = 1; i <= n; i++)
      {  fin[i] = fou[i] = 0;
         tempin[i] = tempou[i] = 0;
      }
      for (i = 1; i <= na; i++)
      {  nxtin[i] = nxtou[i] = 0;
         i1 = startn[i];
         i2 = endn[i];
         if (fou[i1] != 0)
            nxtou[tempou[i1]] = i;
         else
            fou[i1] = i;
         tempou[i1] = i;
         if (fin[i2] != 0)
            nxtin[tempin[i2]] = i;
         else
            fin[i2] = i;
         tempin[i2] = i;
      }
      return;
}

/***********************************************************************
*  NAME
*
*  ascnt1 - multi-node price adjustment for positive deficit case
*
*  PURPOSE
*
*  This subroutine performs the multi-node price adjustment step for
*  the case where the scanned nodes have positive deficit. It first
*  checks if decreasing the price of the scanned nodes increases the
*  dual cost. If yes, then it decreases the price of all scanned nodes.
*  There are two possibilities for price decrease: if switch = true,
*  then the set of scanned nodes corresponds to an elementary direction
*  of maximal rate of ascent, in which case the price of all scanned
*  nodes are decreased until the next breakpoint in the dual cost is
*  encountered. At this point, some arc becomes balanced and more
*  node(s) are added to the labeled set and the subroutine is exited.
*  If switch = false, then the price of all scanned nodes are decreased
*  until the rate of ascent becomes negative (this corresponds to the
*  price adjustment step in which both the line search and the
*  degenerate ascent iteration are implemented).
*
*  INPUT PARAMETERS
*
*  dm        = total deficit of scanned nodes
*  switch    = true if labeling is to continue after price change
*  nscan     = number of scanned nodes
*  curnode   = most recently scanned node
*  n         = number of nodes
*  na        = number of arcs
*  large     = a very large integer to represent infinity (see note 3)
*  startn[i] = starting node for the i-th arc, i = 1,...,na
*  endn[i]   = ending node for the i-th arc, i = 1,...,na
*  fou[i]    = first arc leaving i-th node, i = 1,...,n
*  nxtou[i]  = next arc leaving the starting node of j-th arc,
*              i = 1,...,na
*  fin[i]    = first arc entering i-th node, i = 1,...,n
*  nxtin[i]  = next arc entering the ending node of j-th arc,
*              i = 1,...,na
*
*  UPDATED PARAMETERS
*
*  delx      = a lower estimate of the total flow on balanced arcs in
*              the scanned-nodes cut
*  nlabel    = number of labeled nodes
*  feasbl    = false if problem is found to be infeasible
*  prevnode  = the node before curnode in queue
*  rc[j]     = reduced cost of arc j, j = 1,...,na
*  u[j]      = residual capacity of arc j, j = 1,...,na
*  x[j]      = flow on arc j, j = 1,...,na
*  dfct[i]   = deficit at node i, i = 1,...,n
*  label[k]  = k-th node labeled, k = 1,...,nlabel
*  prdcsr[i] = predecessor of node i in tree of labeled nodes (0 if i
*              is unlabeled), i = 1,...,n
*  tfstou[i] = first balanced arc out of node i, i = 1,...,n
*  tnxtou[j] = next balanced arc out of the starting node of arc j,
*              j = 1,...,na
*  tfstin[i] = first balanced arc into node i, i = 1,...,n
*  tnxtin[j] = next balanced arc into the ending node of arc j,
*              j = 1,...,na
*  nxtqueue[i] = node following node i in the fifo queue (0 if node is
*              not in the queue), i = 1,...,n
*  scan[i]   = true if node i is scanned, i = 1,...,n
*  mark[i]   = true if node i is labeled, i = 1,...,n
*
*  WORKING PARAMETERS
*
*  save[1+na] */

static void ascnt1(struct relax4_csa *csa, int dm, int *delx,
      int *nlabel, int *feasbl, int *svitch, int nscan, int curnode,
      int *prevnode)
{     /* input parameters */
      int n = csa->n;
      /* int na = csa->na; */
      int large = csa->large;
      int *startn = csa->startn;
      int *endn = csa->endn;
      int *fou = csa->fou;
      int *nxtou = csa->nxtou;
      int *fin = csa->fin;
      int *nxtin = csa->nxtin;
      /* updated parameters */
#     define delx (*delx)
#     define nlabel (*nlabel)
#     define feasbl (*feasbl)
#     define svitch (*svitch)
#     define prevnode (*prevnode)
      int *rc = csa->rc;
      int *u = csa->u;
      int *x = csa->x;
      int *dfct = csa->dfct;
      int *label = csa->label;
      int *prdcsr = csa->prdcsr;
      int *tfstou = csa->tfstou;
      int *tnxtou = csa->tnxtou;
      int *tfstin = csa->tfstin;
      int *tnxtin = csa->tnxtin;
      int *nxtqueue = csa->nxtqueue;
      char *scan = csa->scan;
      char *mark = csa->mark;
      int *save = csa->save;
      /* local variables */
      int arc, delprc, dlx, i, j, nb, node, node2, nsave, rdcost, t1,
         t2, t3;
      /* Store the arcs between the set of scanned nodes and its
       * complement in save and compute delprc, the stepsize to the next
       * breakpoint in the dual cost in the direction of decreasing
       * prices of the scanned nodes.
       * [The arcs are stored into save by looking at the arcs incident
       * to either the set of scanned nodes or its complement, depending
       * on whether nscan > n/2 or not. This improves the efficiency of
       * storing.] */
      delprc = large;
      dlx = 0;
      nsave = 0;
      if (nscan <= n / 2)
      {  for (i = 1; i <= nscan; i++)
         {  node = label[i];
            for (arc = fou[node]; arc > 0; arc = nxtou[arc])
            {  /* arc points from scanned node to an unscanned node. */
               node2 = endn[arc];
               if (!scan[node2])
               {  nsave++;
                  save[nsave] = arc;
                  rdcost = rc[arc];
                  if ((rdcost == 0) && (prdcsr[node2] != arc))
                     dlx += x[arc];
                  if ((rdcost < 0) && (rdcost > -delprc))
                     delprc = -rdcost;
               }
            }
            for (arc = fin[node]; arc > 0; arc = nxtin[arc])
            {  /* arc points from unscanned node to scanned node. */
               node2 = startn[arc];
               if (!scan[node2])
               {  nsave++;
                  save[nsave] = -arc;
                  rdcost = rc[arc];
                  if ((rdcost == 0) && (prdcsr[node2] != -arc))
                     dlx += u[arc];
                  if ((rdcost > 0) && (rdcost < delprc))
                     delprc = rdcost;
               }
            }
         }
      }
      else
      {  for (node = 1; node <= n; node++)
         {  if (scan[node])
               continue;
            for (arc = fin[node]; arc > 0; arc = nxtin[arc])
            {  node2 = startn[arc];
               if (scan[node2])
               {  nsave++;
                  save[nsave] = arc;
                  rdcost = rc[arc];
                  if ((rdcost == 0) && (prdcsr[node] != arc))
                     dlx += x[arc];
                  if ((rdcost < 0) && (rdcost > -delprc))
                     delprc = -rdcost;
               }
            }
            for (arc = fou[node]; arc > 0; arc = nxtou[arc])
            {  node2 = endn[arc];
               if (scan[node2])
               {  nsave++;
                  save[nsave] = -arc;
                  rdcost = rc[arc];
                  if ((rdcost == 0) && (prdcsr[node] != -arc))
                     dlx += u[arc];
                  if ((rdcost > 0) && (rdcost < delprc))
                     delprc = rdcost;
               }
            }
         }
      }
      /* Check if the set of scanned nodes truly corresponds to a dual
       * ascent direction. [Here delx + dlx is the exact sum of the flow
       * on arcs from the scanned set to the unscanned set plus the
       * (capacity - flow) on arcs from the unscanned set to the scanned
       * set.] If this were not the case, set switch to true and exit
       * subroutine. */
      if ((!svitch) && (delx + dlx >= dm))
      {  svitch = true;
         return;
      }
      delx += dlx;
L4:   /* Check that the problem is feasible. */
      if (delprc == large)
      {  /* We can increase the dual cost without bound, so the primal
          * problem is infeasible. */
         feasbl = false;
         return;
      }
      /* Decrease the prices of the scanned nodes, add more nodes to
       * the labeled set and check if a newly labeled node has negative
       * deficit. */
      if (svitch)
      {  for (i = 1; i <= nsave; i++)
         {  arc = save[i];
            if (arc > 0)
            {  rc[arc] += delprc;
               if (rc[arc] == 0)
               {  node2 = endn[arc];
                  if (tnxtou[arc] < 0)
                  {  tnxtou[arc] = tfstou[startn[arc]];
                     tfstou[startn[arc]] = arc;
                  }
                  if (tnxtin[arc] < 0)
                  {  tnxtin[arc] = tfstin[node2];
                     tfstin[node2] = arc;
                  }
                  if (!mark[node2])
                  {  prdcsr[node2] = arc;
                     nlabel++;
                     label[nlabel] = node2;
                     mark[node2] = true;
                  }
               }
            }
            else
            {  arc = -arc;
               rc[arc] -= delprc;
               if (rc[arc] == 0)
               {  node2 = startn[arc];
                  if (tnxtou[arc] < 0)
                  {  tnxtou[arc] = tfstou[node2];
                     tfstou[node2] = arc;
                  }
                  if (tnxtin[arc] < 0)
                  {  tnxtin[arc] = tfstin[endn[arc]];
                     tfstin[endn[arc]] = arc;
                  }
                  if (!mark[node2])
                  {  prdcsr[node2] = -arc;
                     nlabel++;
                     label[nlabel] = node2;
                     mark[node2] = true;
                  }
               }
            }
         }
         return;
      }
      else
      {  /* Decrease the prices of the scanned nodes by delprc. Adjust
          * flow to maintain complementary slackness with the prices. */
         nb = 0;
         for (i = 1; i <= nsave; i++)
         {  arc = save[i];
            if (arc > 0)
            {  t1 = rc[arc];
               if (t1 == 0)
               {  t2 = x[arc];
                  t3 = startn[arc];
                  dfct[t3] -= t2;
                  if (nxtqueue[t3] == 0)
                  {  nxtqueue[prevnode] = t3;
                     nxtqueue[t3] = curnode;
                     prevnode = t3;
                  }
                  t3 = endn[arc];
                  dfct[t3] += t2;
                  if (nxtqueue[t3] == 0)
                  {  nxtqueue[prevnode] = t3;
                     nxtqueue[t3] = curnode;
                     prevnode = t3;
                  }
                  u[arc] += t2;
                  x[arc] = 0;
               }
               rc[arc] = t1 + delprc;
#if 0 /* by mao; 26/IV-2013 */
               if (rc[arc] == 0)
#else
               if (rc[arc] == 0 && nb < n)
#endif
               {  delx += x[arc];
                  nb++;
                  prdcsr[nb] = arc;
               }
            }
            else
            {  arc = -arc;
               t1 = rc[arc];
               if (t1 == 0)
               {  t2 = u[arc];
                  t3 = startn[arc];
                  dfct[t3] += t2;
                  if (nxtqueue[t3] == 0)
                  {  nxtqueue[prevnode] = t3;
                     nxtqueue[t3] = curnode;
                     prevnode = t3;
                  }
                  t3 = endn[arc];
                  dfct[t3] -= t2;
                  if (nxtqueue[t3] == 0)
                  {  nxtqueue[prevnode] = t3;
                     nxtqueue[t3] = curnode;
                     prevnode = t3;
                  }
                  x[arc] += t2;
                  u[arc] = 0;
               }
               rc[arc] = t1 - delprc;
#if 0 /* by mao; 26/IV-2013 */
               if (rc[arc] == 0)
#else
               if (rc[arc] == 0 && nb < n)
#endif
               {  delx += u[arc];
                  nb++;
                  prdcsr[nb] = arc;
               }
            }
         }
      }
      if (delx <= dm)
      {  /* The set of scanned nodes still corresponds to a dual
          * (possibly degenerate) ascent direction. Compute the stepsize
          * delprc to the next breakpoint in the dual cost. */
         delprc = large;
         for (i = 1; i <= nsave; i++)
         {  arc = save[i];
            if (arc > 0)
            {  rdcost = rc[arc];
               if ((rdcost < 0) && (rdcost > -delprc))
                  delprc = -rdcost;
            }
            else
            {  arc = -arc;
               rdcost = rc[arc];
               if ((rdcost > 0) && (rdcost < delprc))
                  delprc = rdcost;
            }
         }
         if ((delprc != large) || (delx < dm))
            goto L4;
      }
      /* Add new balanced arcs to the superset of balanced arcs. */
      for (i = 1; i <= nb; i++)
      {  arc = prdcsr[i];
         if (tnxtin[arc] == -1)
         {  j = endn[arc];
            tnxtin[arc] = tfstin[j];
            tfstin[j] = arc;
         }
         if (tnxtou[arc] == -1)
         {  j = startn[arc];
            tnxtou[arc] = tfstou[j];
            tfstou[j] = arc;
         }
      }
      return;
#     undef delx
#     undef nlabel
#     undef feasbl
#     undef svitch
#     undef prevnode
}

/***********************************************************************
*  NAME
*
*  ascnt2 - multi-node price adjustment for negative deficit case
*
*  PURPOSE
*
*  This routine is analogous to ascnt1 but for the case where the
*  scanned nodes have negative deficit. */

static void ascnt2(struct relax4_csa *csa, int dm, int *delx,
      int *nlabel, int *feasbl, int *svitch, int nscan, int curnode,
      int *prevnode)
{     /* input parameters */
      int n = csa->n;
      /* int na = csa->na; */
      int large = csa->large;
      int *startn = csa->startn;
      int *endn = csa->endn;
      int *fou = csa->fou;
      int *nxtou = csa->nxtou;
      int *fin = csa->fin;
      int *nxtin = csa->nxtin;
      /* updated parameters */
#     define delx (*delx)
#     define nlabel (*nlabel)
#     define feasbl (*feasbl)
#     define svitch (*svitch)
#     define prevnode (*prevnode)
      int *rc = csa->rc;
      int *u = csa->u;
      int *x = csa->x;
      int *dfct = csa->dfct;
      int *label = csa->label;
      int *prdcsr = csa->prdcsr;
      int *tfstou = csa->tfstou;
      int *tnxtou = csa->tnxtou;
      int *tfstin = csa->tfstin;
      int *tnxtin = csa->tnxtin;
      int *nxtqueue = csa->nxtqueue;
      char *scan = csa->scan;
      char *mark = csa->mark;
      int *save = csa->save;
      /* local variables */
      int arc, delprc, dlx, i, j, nb, node, node2, nsave, rdcost, t1,
         t2, t3;
      /* Store the arcs between the set of scanned nodes and its
       * complement in save and compute delprc, the stepsize to the next
       * breakpoint in the dual cost in the direction of increasing
       * prices of the scanned nodes. */
      delprc = large;
      dlx = 0;
      nsave = 0;
      if (nscan <= n / 2)
      {  for (i = 1; i <= nscan; i++)
         {  node = label[i];
            for (arc = fin[node]; arc > 0; arc = nxtin[arc])
            {  node2 = startn[arc];
               if (!scan[node2])
               {  nsave++;
                  save[nsave] = arc;
                  rdcost = rc[arc];
                  if ((rdcost == 0) && (prdcsr[node2] != arc))
                     dlx += x[arc];
                  if ((rdcost < 0) && (rdcost > -delprc))
                     delprc = -rdcost;
               }
            }
            for (arc = fou[node]; arc > 0; arc = nxtou[arc])
            {  node2 = endn[arc];
               if (!scan[node2])
               {  nsave++;
                  save[nsave] = -arc;
                  rdcost = rc[arc];
                  if ((rdcost == 0) && (prdcsr[node2] != -arc))
                     dlx += u[arc];
                  if ((rdcost > 0) && (rdcost < delprc))
                     delprc = rdcost;
               }
            }
         }
      }
      else
      {  for (node = 1; node <= n; node++)
         {  if (scan[node])
               continue;
            for (arc = fou[node]; arc > 0; arc = nxtou[arc])
            {  node2 = endn[arc];
               if (scan[node2])
               {  nsave++;
                  save[nsave] = arc;
                  rdcost = rc[arc];
                  if ((rdcost == 0) && (prdcsr[node] != arc))
                     dlx += x[arc];
                  if ((rdcost < 0) && (rdcost > -delprc))
                     delprc = -rdcost;
               }
            }
            for (arc = fin[node]; arc > 0; arc = nxtin[arc])
            {  node2 = startn[arc];
               if (scan[node2])
               {  nsave++;
                  save[nsave] = -arc;
                  rdcost = rc[arc];
                  if ((rdcost == 0) && (prdcsr[node] != -arc))
                     dlx += u[arc];
                  if ((rdcost > 0) && (rdcost < delprc))
                     delprc = rdcost;
               }
            }
         }
      }
      if ((!svitch) && (delx + dlx >= -dm))
      {  svitch = true;
         return;
      }
      delx += dlx;
      /* Check that the problem is feasible. */
L4:   if (delprc == large)
      {  feasbl = false;
         return;
      }
      /* Increase the prices of the scanned nodes, add more nodes to
       * the labeled set and check if a newly labeled node has positive
       * deficit. */
      if (svitch)
      {  for (i = 1; i <= nsave; i++)
         {  arc = save[i];
            if (arc > 0)
            {  rc[arc] += delprc;
               if (rc[arc] == 0)
               {  node2 = startn[arc];
                  if (tnxtou[arc] < 0)
                  {  tnxtou[arc] = tfstou[node2];
                     tfstou[node2] = arc;
                  }
                  if (tnxtin[arc] < 0)
                  {  tnxtin[arc] = tfstin[endn[arc]];
                     tfstin[endn[arc]] = arc;
                  }
                  if (!mark[node2])
                  {  prdcsr[node2] = arc;
                     nlabel++;
                     label[nlabel] = node2;
                     mark[node2] = true;
                  }
               }
            }
            else
            {  arc = -arc;
               rc[arc] -= delprc;
               if (rc[arc] == 0)
               {  node2 = endn[arc];
                  if (tnxtou[arc] < 0)
                  {  tnxtou[arc] = tfstou[startn[arc]];
                     tfstou[startn[arc]] = arc;
                  }
                  if (tnxtin[arc] < 0)
                  {  tnxtin[arc] = tfstin[node2];
                     tfstin[node2] = arc;
                  }
                  if (!mark[node2])
                  {  prdcsr[node2] = -arc;
                     nlabel++;
                     label[nlabel] = node2;
                     mark[node2] = true;
                  }
               }
            }
         }
         return;
      }
      else
      {  nb = 0;
         for (i = 1; i <= nsave; i++)
         {  arc = save[i];
            if (arc > 0)
            {  t1 = rc[arc];
               if (t1 == 0)
               {  t2 = x[arc];
                  t3 = startn[arc];
                  dfct[t3] -= t2;
                  if (nxtqueue[t3] == 0)
                  {  nxtqueue[prevnode] = t3;
                     nxtqueue[t3] = curnode;
                     prevnode = t3;
                  }
                  t3 = endn[arc];
                  dfct[t3] += t2;
                  if (nxtqueue[t3] == 0)
                  {  nxtqueue[prevnode] = t3;
                     nxtqueue[t3] = curnode;
                     prevnode = t3;
                  }
                  u[arc] += t2;
                  x[arc] = 0;
               }
               rc[arc] = t1 + delprc;
#if 0 /* by mao; 26/IV-2013 */
               if (rc[arc] == 0)
#else
               if (rc[arc] == 0 && nb < n)
#endif
               {  delx += x[arc];
                  nb++;
                  prdcsr[nb] = arc;
               }
            }
            else
            {  arc = -arc;
               t1 = rc[arc];
               if (t1 == 0)
               {  t2 = u[arc];
                  t3 = startn[arc];
                  dfct[t3] += t2;
                  if (nxtqueue[t3] == 0)
                  {  nxtqueue[prevnode] = t3;
                     nxtqueue[t3] = curnode;
                     prevnode = t3;
                  }
                  t3 = endn[arc];
                  dfct[t3] -= t2;
                  if (nxtqueue[t3] == 0)
                  {  nxtqueue[prevnode] = t3;
                     nxtqueue[t3] = curnode;
                     prevnode = t3;
                  }
                  x[arc] += t2;
                  u[arc] = 0;
               }
               rc[arc] = t1 - delprc;
#if 0 /* by mao; 26/IV-2013 */
               if (rc[arc] == 0)
#else
               if (rc[arc] == 0 && nb < n)
#endif
               {  delx += u[arc];
                  nb++;
                  prdcsr[nb] = arc;
               }
            }
         }
      }
      if (delx <= -dm)
      {  delprc = large;
         for (i = 1; i <= nsave; i++)
         {  arc = save[i];
            if (arc > 0)
            {  rdcost = rc[arc];
               if ((rdcost < 0) && (rdcost > -delprc))
                  delprc = -rdcost;
            }
            else
            {  arc = -arc;
               rdcost = rc[arc];
               if ((rdcost > 0) && (rdcost < delprc))
                  delprc = rdcost;
            }
         }
         if ((delprc != large) || (delx < -dm))
            goto L4;
      }
      /* Add new balanced arcs to the superset of balanced arcs. */
      for (i = 1; i <= nb; i++)
      {  arc = prdcsr[i];
         if (tnxtin[arc] == -1)
         {  j = endn[arc];
            tnxtin[arc] = tfstin[j];
            tfstin[j] = arc;
         }
         if (tnxtou[arc] == -1)
         {  j = startn[arc];
            tnxtou[arc] = tfstou[j];
            tfstou[j] = arc;
         }
      }
      return;
#     undef delx
#     undef nlabel
#     undef feasbl
#     undef svitch
#     undef prevnode
}

/***********************************************************************
*  NAME
*
*  auction - compute good initial flow and prices
*
*  PURPOSE
*
*  This subroutine uses a version of the auction algorithm for min
*  cost network flow to compute a good initial flow and prices for the
*  problem.
*
*  INPUT PARAMETERS
*
*  n         = number of nodes
*  na        = number of arcs
*  large     = a very large integer to represent infinity (see note 3)
*  startn[i] = starting node for the i-th arc, i = 1,...,na
*  endn[i]   = ending node for the i-th arc, i = 1,...,na
*  fou[i]    = first arc leaving i-th node, i = 1,...,n
*  nxtou[i]  = next arc leaving the starting node of j-th arc,
*              i = 1,...,na
*  fin[i]    = first arc entering i-th node, i = 1,...,n
*  nxtin[i]  = next arc entering the ending node of j-th arc,
*              i = 1,...,na
*
*  UPDATED PARAMETERS
*
*  rc[j]     = reduced cost of arc j, j = 1,...,na
*  u[j]      = residual capacity of arc j, j = 1,...,na
*  x[j]      = flow on arc j, j = 1,...,na
*  dfct[i]   = deficit at node i, i = 1,...,n
*
*  OUTPUT PARAMETERS
*
*  nsp       = number of auction/shortest path iterations
*
*  WORKING PARAMETERS
*
*  p[1+n], prdcsr[1+n], save[1+na], fpushf[1+n], nxtpushf[1+na],
*  fpushb[1+n], nxtpushb[1+na], nxtqueue[1+n], extend_arc[1+n],
*  sb_level[1+n], sb_arc[1+n], path_id[1+n]
*
*  RETURNS
*
*  0         = normal return
*  1         = problem is found to be infeasible */

static int auction(struct relax4_csa *csa)
{     /* input parameters */
      int n = csa->n;
      int na = csa->na;
      int large = csa->large;
      int *startn = csa->startn;
      int *endn = csa->endn;
      int *fou = csa->fou;
      int *nxtou = csa->nxtou;
      int *fin = csa->fin;
      int *nxtin = csa->nxtin;
      /* updated parameters */
#     define crash (csa->crash)
      int *rc = csa->rc;
      int *u = csa->u;
      int *x = csa->x;
      int *dfct = csa->dfct;
      /* output parameters */
#     define nsp (csa->nsp)
      /* working parameters */
      int *p = csa->label;
      int *prdcsr = csa->prdcsr;
      int *save = csa->save;
      int *fpushf = csa->tfstou;
      int *nxtpushf = csa->tnxtou;
      int *fpushb = csa->tfstin;
      int *nxtpushb = csa->tnxtin;
      int *nxtqueue = csa->nxtqueue;
      int *extend_arc = csa->extend_arc;
      int *sb_level = csa->sb_level;
      int *sb_arc = csa->sb_arc;
      char *path_id = csa->mark;
      /* local variables */
      int arc, bstlevel, end, eps, extarc, factor, flow, i, incr,
         last, lastqueue, maxcost, mincost, nas, naug, new_level, node,
         nolist, num_passes, nxtnode, pass, pend, pr_term, prd,
         prevarc, prevlevel, prevnode, pstart, pterm, rdcost, red_cost,
         resid, root, secarc, seclevel, start, term, thresh_dfct;
      /* start initialization using auction */
      naug = 0;
      pass = 0;
      thresh_dfct = 0;
      /* factor determines by how much epsilon is reduced at each
       * minimization */
      factor = 3;
      /* num_passes determines how many auction scaling phases are
       * performed */
      num_passes = 1;
      /* set arc flows to satisfy cs and calculate maxcost and
       * mincost */
      maxcost = -large;
      mincost = large;
      for (arc = 1; arc <= na; arc++)
      {  start = startn[arc];
         end = endn[arc];
         rdcost = rc[arc];
         if (maxcost < rdcost)
            maxcost = rdcost;
         if (mincost > rdcost)
            mincost = rdcost;
         if (rdcost < 0)
         {  dfct[start] += u[arc];
            dfct[end] -= u[arc];
            x[arc] = u[arc];
            u[arc] = 0;
         }
         else
            x[arc] = 0;
      }
      /* set initial epsilon */
      if ((maxcost - mincost) >= 8)
         eps = (maxcost - mincost) / 8;
      else
         eps = 1;
      /* set initial prices to zero */
      for (node = 1; node <= n; node++)
         p[node] = 0;
      /* Initialization using auction/shortest paths. */
L100: /* Start of the first scaling phase. */
      pass++;
      if ((pass == num_passes) || (eps == 1))
         crash = 0;
      nolist = 0;
      /* construct list of positive surplus nodes and queue of negative
       * surplus nodes */
      for (node = 1; node <= n; node++)
      {  prdcsr[node] = 0;
         path_id[node] = false;
         extend_arc[node] = 0;
         sb_level[node] = -large;
         nxtqueue[node] = node + 1;
         if (dfct[node] > 0)
         {  nolist++;
            save[nolist] = node;
         }
      }
      nxtqueue[n] = 1;
      root = 1;
      prevnode = lastqueue = n;
      /* initialization with down iterations for negative surplus
       * nodes */
      for (i = 1; i <= nolist; i++)
      {  node = save[i];
         nsp++;
         /* build the list of arcs w/ room for pushing flow and find
          * proper price for down iteration */
         bstlevel = -large;
         fpushf[node] = 0;
         for (arc = fou[node]; arc > 0; arc = nxtou[arc])
         {  if (u[arc] > 0)
            {  if (fpushf[node] == 0)
               {  fpushf[node] = arc;
                  nxtpushf[arc] = 0;
                  last = arc;
               }
               else
               {  nxtpushf[last] = arc;
                  nxtpushf[arc] = 0;
                  last = arc;
               }
            }
            if (x[arc] > 0)
            {  new_level = p[endn[arc]] + rc[arc];
               if (new_level > bstlevel)
               {  bstlevel = new_level;
                  extarc = arc;
               }
            }
         }
         fpushb[node] = 0;
         for (arc = fin[node]; arc > 0; arc = nxtin[arc])
         {  if (x[arc] > 0)
            {  if (fpushb[node] == 0)
               {  fpushb[node] = arc;
                  nxtpushb[arc] = 0;
                  last = arc;
               }
               else
               {  nxtpushb[last] = arc;
                  nxtpushb[arc] = 0;
                  last = arc;
               }
            }
            if (u[arc] > 0)
            {  new_level = p[startn[arc]] - rc[arc];
               if (new_level > bstlevel)
               {  bstlevel = new_level;
                  extarc = -arc;
               }
            }
         }
         extend_arc[node] = extarc;
         p[node] = bstlevel - eps;
      }
L200: /* Start the augmentation cycles of the new scaling phase. */
      if (dfct[root] >= thresh_dfct)
         goto L3000;
      term = root;
      path_id[root] = true;
L500: /* Main forward algorithm with root as origin. */
      /* start of a new forward iteration */
      pterm = p[term];
      extarc = extend_arc[term];
      if (extarc == 0)
      {  /* build the list of arcs w/ room for pushing flow */
         fpushf[term] = 0;
         for (arc = fou[term]; arc > 0; arc = nxtou[arc])
         {  if (u[arc] > 0)
            {  if (fpushf[term] == 0)
               {  fpushf[term] = arc;
                  nxtpushf[arc] = 0;
                  last = arc;
               }
               else
               {  nxtpushf[last] = arc;
                  nxtpushf[arc] = 0;
                  last = arc;
               }
            }
         }
         fpushb[term] = 0;
         for (arc = fin[term]; arc > 0; arc = nxtin[arc])
         {  if (x[arc] > 0)
            {  if (fpushb[term] == 0)
               {  fpushb[term] = arc;
                  nxtpushb[arc] = 0;
                  last = arc;
               }
               else
               {  nxtpushb[last] = arc;
                  nxtpushb[arc] = 0;
                  last = arc;
               }
            }
         }
         goto L600;
      }
      /* speculative path extension attempt */
      /* note: arc > 0 means that arc is oriented from the root to the
       * destinations
       * arc < 0 means that arc is oriented from the destinations to the
       * root
       * extarc = 0 or prdarc = 0, means the extension arc or the
       * predecessor arc, respectively, has not been established */
      if (extarc > 0)
      {  if (u[extarc] == 0)
         {  seclevel = sb_level[term];
            goto L580;
         }
         end = endn[extarc];
         bstlevel = p[end] + rc[extarc];
         if (pterm >= bstlevel)
         {  if (path_id[end])
               goto L1200;
            term = end;
            prdcsr[term] = extarc;
            path_id[term] = true;
            /* if negative surplus node is found, do an augmentation */
            if (dfct[term] > 0)
               goto L2000;
            /* return for another iteration */
            goto L500;
         }
      }
      else
      {  extarc = -extarc;
         if (x[extarc] == 0)
         {  seclevel = sb_level[term];
            goto L580;
         }
         start = startn[extarc];
         bstlevel = p[start] - rc[extarc];
         if (pterm >= bstlevel)
         {  if (path_id[start])
               goto L1200;
            term = start;
            prdcsr[term] = -extarc;
            path_id[term] = true;
            /* if negative surplus node is found, do an augmentation */
            if (dfct[term] > 0)
               goto L2000;
            /* return for another iteration */
            goto L500;
         }
      }
L550: /* second best logic test applied to save a full node scan
       * if old best level continues to be best go for another
       * contraction */
      seclevel = sb_level[term];
      if (bstlevel <= seclevel)
         goto L800;
L580: /* if second best can be used do either a contraction or start
       * over with a speculative extension */
      if (seclevel > -large)
      {  extarc = sb_arc[term];
         if (extarc > 0)
         {  if (u[extarc] == 0)
               goto L600;
            bstlevel = p[endn[extarc]] + rc[extarc];
         }
         else
         {  if (x[-extarc] == 0)
               goto L600;
            bstlevel = p[startn[-extarc]] - rc[-extarc];
         }
         if (bstlevel == seclevel)
         {  sb_level[term] = -large;
            extend_arc[term] = extarc;
            goto L800;
         }
      }
L600: /* extension/contraction attempt was unsuccessful, so scan
       * terminal node */
      nsp++;
      bstlevel = seclevel = large;
      for (arc = fpushf[term]; arc > 0; arc = nxtpushf[arc])
      {  new_level = p[endn[arc]] + rc[arc];
         if (new_level < seclevel)
         {  if (new_level < bstlevel)
            {  seclevel = bstlevel;
               bstlevel = new_level;
               secarc = extarc;
               extarc = arc;
            }
            else
            {  seclevel = new_level;
               secarc = arc;
            }
         }
      }
      for (arc = fpushb[term]; arc > 0; arc = nxtpushb[arc])
      {  new_level = p[startn[arc]] - rc[arc];
         if (new_level < seclevel)
         {  if (new_level < bstlevel)
            {  seclevel = bstlevel;
               bstlevel = new_level;
               secarc = extarc;
               extarc = -arc;
            }
            else
            {  seclevel = new_level;
               secarc = -arc;
            }
         }
      }
      sb_level[term] = seclevel;
      sb_arc[term] = secarc;
      extend_arc[term] = extarc;
L800: /* End of node scan. */
      /* if the terminal node is the root, adjust its price and change
       * root */
      if (term == root)
      {  p[term] = bstlevel + eps;
         if (p[term] >= large)
         {  /* no path to the destination */
            /* problem is found to be infeasible */
            return 1;
         }
         path_id[root] = false;
         prevnode = root;
         root = nxtqueue[root];
         goto L200;
      }
      /* check whether extension or contraction */
      prd = prdcsr[term];
      if (prd > 0)
      {  pr_term = startn[prd];
         prevlevel = p[pr_term] - rc[prd];
      }
      else
      {  pr_term = endn[-prd];
         prevlevel = p[pr_term] + rc[-prd];
      }
      if (prevlevel > bstlevel)
      {  /* path extension */
         if (prevlevel >= bstlevel + eps)
            p[term] = bstlevel + eps;
         else
            p[term] = prevlevel;
         if (extarc > 0)
         {  end = endn[extarc];
            if (path_id[end])
               goto L1200;
            term = end;
         }
         else
         {  start = startn[-extarc];
            if (path_id[start])
               goto L1200;
            term = start;
         }
         prdcsr[term] = extarc;
         path_id[term] = true;
         /* if negative surplus node is found, do an augmentation */
         if (dfct[term] > 0)
            goto L2000;
         /* return for another iteration */
         goto L500;
      }
      else
      {  /* path contraction */
         p[term] = bstlevel + eps;
         path_id[term] = false;
         term = pr_term;
         if (pr_term != root)
         {  if (bstlevel <= pterm + eps)
               goto L2000;
         }
         pterm = p[term];
         extarc = prd;
         if (prd > 0)
            bstlevel += eps + rc[prd];
         else
            bstlevel += eps - rc[-prd];
         /* do a second best test and if that fails, do a full node
          * scan */
         goto L550;
      }
L1200:/* A cycle is about to form; do a retreat sequence. */
      node = term;
L1600:if (node != root)
      {  path_id[node] = false;
         prd = prdcsr[node];
         if (prd > 0)
         {  pr_term = startn[prd];
            if (p[pr_term] == p[node] + rc[prd] + eps)
            {  node = pr_term;
               goto L1600;
            }
         }
         else
         {  pr_term = endn[-prd];
            if (p[pr_term] == p[node] - rc[-prd] + eps)
            {  node = pr_term;
               goto L1600;
            }
         }
         /* do a full scan and price rise at pr_term */
         nsp++;
         bstlevel = seclevel = large;
         for (arc = fpushf[pr_term]; arc > 0; arc = nxtpushf[arc])
         {  new_level = p[endn[arc]] + rc[arc];
            if (new_level < seclevel)
            {  if (new_level < bstlevel)
               {  seclevel = bstlevel;
                  bstlevel = new_level;
                  secarc = extarc;
                  extarc = arc;
               }
               else
               {  seclevel = new_level;
                  secarc = arc;
               }
            }
         }
         for (arc = fpushb[pr_term]; arc > 0; arc = nxtpushb[arc])
         {  new_level = p[startn[arc]] - rc[arc];
            if (new_level < seclevel)
            {  if (new_level < bstlevel)
               {  seclevel = bstlevel;
                  bstlevel = new_level;
                  secarc = extarc;
                  extarc = -arc;
               }
               else
               {  seclevel = new_level;
                  secarc = -arc;
               }
            }
         }
         sb_level[pr_term] = seclevel;
         sb_arc[pr_term] = secarc;
         extend_arc[pr_term] = extarc;
         p[pr_term] = bstlevel + eps;
         if (pr_term == root)
         {  prevnode = root;
            path_id[root] = false;
            root = nxtqueue[root];
            goto L200;
         }
         path_id[pr_term] = false;
         prd = prdcsr[pr_term];
         if (prd > 0)
            term = startn[prd];
         else
            term = endn[-prd];
         if (term == root)
         {  prevnode = root;
            path_id[root] = false;
            root = nxtqueue[root];
            goto L200;
         }
         else
            goto L2000;
      }
L2000:/* End of auction/shortest path routine. */
      /* do augmentation from root and correct the push lists */
      incr = -dfct[root];
      for (node = root;;)
      {  extarc = extend_arc[node];
         path_id[node] = false;
         if (extarc > 0)
         {  node = endn[extarc];
            if (incr > u[extarc])
               incr = u[extarc];
         }
         else
         {  node = startn[-extarc];
            if (incr > x[-extarc])
               incr = x[-extarc];
         }
         if (node == term)
            break;
      }
      path_id[term] = false;
      if (dfct[term] > 0)
      {  if (incr > dfct[term])
            incr = dfct[term];
      }
      for (node = root;;)
      {  extarc = extend_arc[node];
         if (extarc > 0)
         {  end = endn[extarc];
            /* add arc to the reduced graph */
            if (x[extarc] == 0)
            {  nxtpushb[extarc] = fpushb[end];
               fpushb[end] = extarc;
               new_level = p[node] - rc[extarc];
               if (sb_level[end] > new_level)
               {  sb_level[end] = new_level;
                  sb_arc[end] = -extarc;
               }
            }
            x[extarc] += incr;
            u[extarc] -= incr;
            /* remove arc from the reduced graph */
            if (u[extarc] == 0)
            {  nas++;
               arc = fpushf[node];
               if (arc == extarc)
                  fpushf[node] = nxtpushf[arc];
               else
               {  prevarc = arc;
                  arc = nxtpushf[arc];
                  while (arc > 0)
                  {  if (arc == extarc)
                     {  nxtpushf[prevarc] = nxtpushf[arc];
                        break;
                     }
                     prevarc = arc;
                     arc = nxtpushf[arc];
                  }
               }
            }
            node = end;
         }
         else
         {  extarc = -extarc;
            start = startn[extarc];
            /* add arc to the reduced graph */
            if (u[extarc] == 0)
            {  nxtpushf[extarc] = fpushf[start];
               fpushf[start] = extarc;
               new_level = p[node] + rc[extarc];
               if (sb_level[start] > new_level)
               {  sb_level[start] = new_level;
                  sb_arc[start] = extarc;
               }
            }
            u[extarc] += incr;
            x[extarc] -= incr;
            /* remove arc from the reduced graph */
            if (x[extarc] == 0)
            {  nas++;
               arc = fpushb[node];
               if (arc == extarc)
                  fpushb[node] = nxtpushb[arc];
               else
               {  prevarc = arc;
                  arc = nxtpushb[arc];
                  while (arc > 0)
                  {  if (arc == extarc)
                     {  nxtpushb[prevarc] = nxtpushb[arc];
                        break;
                     }
                     prevarc = arc;
                     arc = nxtpushb[arc];
                  }
               }
            }
            node = start;
         }
         if (node == term)
            break;
      }
      dfct[term] -= incr;
      dfct[root] += incr;
      /* insert term in the queue if it has a large enough surplus */
      if (dfct[term] < thresh_dfct)
      {  if (nxtqueue[term] == 0)
         {  nxtnode = nxtqueue[root];
            if ((p[term] >= p[nxtnode]) && (root != nxtnode))
            {  nxtqueue[root] = term;
               nxtqueue[term] = nxtnode;
            }
            else
            {  nxtqueue[prevnode] = term;
               nxtqueue[term] = root;
               prevnode = term;
            }
         }
      }
      /* if root has a large enough surplus, keep it in the queue and
       * return for another iteration */
      if (dfct[root] < thresh_dfct)
      {  prevnode = root;
         root = nxtqueue[root];
         goto L200;
      }
L3000:/* end of augmentation cycle */
      /* Check for termination of scaling phase. If scaling phase is not
       * finished, advance the queue and return to take another node. */
      nxtnode = nxtqueue[root];
      if (root != nxtnode)
      {  nxtqueue[root] = 0;
         nxtqueue[prevnode] = nxtnode;
         root = nxtnode;
         goto L200;
      }
      /* End of subproblem (scaling phase). */
      /* Reduce epsilon. */
      eps /= factor;
      if (eps < 1) eps = 1;
      thresh_dfct /= factor;
      if (eps == 1) thresh_dfct = 0;
      /* if another auction scaling phase remains, reset the flows &
       * the push lists; else reset arc flows to satisfy cs and compute
       * reduced costs */
      if (crash == 1)
      {  for (arc = 1; arc <= na; arc++)
         {  start = startn[arc];
            end = endn[arc];
            pstart = p[start];
            pend = p[end];
            if (pstart > pend + eps + rc[arc])
            {  resid = u[arc];
               if (resid > 0)
               {  dfct[start] += resid;
                  dfct[end] -= resid;
                  x[arc] += resid;
                  u[arc] = 0;
               }
            }
            else if (pstart < pend - eps + rc[arc])
            {  flow = x[arc];
               if (flow > 0)
               {  dfct[start] -= flow;
                  dfct[end] += flow;
                  x[arc] = 0;
                  u[arc] += flow;
               }
            }
         }
         /* return for another phase */
         goto L100;
      }
      else
      {  crash = 1;
         for (arc = 1; arc <= na; arc++)
         {  start = startn[arc];
            end = endn[arc];
            red_cost = rc[arc] + p[end] - p[start];
            if (red_cost < 0)
            {  resid = u[arc];
               if (resid > 0)
               {  dfct[start] += resid;
                  dfct[end] -= resid;
                  x[arc] += resid;
                  u[arc] = 0;
               }
            }
            else if (red_cost > 0)
            {  flow = x[arc];
               if (flow > 0)
               {  dfct[start] -= flow;
                  dfct[end] += flow;
                  x[arc] = 0;
                  u[arc] += flow;
               }
            }
            rc[arc] = red_cost;
         }
      }
      return 0;
#     undef crash
#     undef nsp
}

/* eof */
