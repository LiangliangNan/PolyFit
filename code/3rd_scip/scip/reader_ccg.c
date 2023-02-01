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

/**@file   reader_ccg.c
 * @brief  Graph file reader (actually, only a writer)
 * @author Marc Pfetsch
 *
 * Write a weighted column/variable graph, i.e., the nodes correspond to the columns (variables) of
 * the constraint matrix. Two nodes are adjacent if the corresponding columns/variables appear
 * in a common row/constraint (with nonzero coefficient).  The weight is obtained by summing for
 * each row that produces an edge the absolute values of coefficients in the row; hence, we avoid
 * parallel edges.
 *
 * This graph gives an indication of the connectivity structure of the constraint matrix.
 *
 * The graph is output in DIMACS graph format.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "scip/reader_ccg.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"

#define READER_NAME             "ccgreader"
#define READER_DESC             "file writer for column connectivity graph file format"
#define READER_EXTENSION        "ccg"

/*
 * Data structures
 */


/* graph data structure */
struct sparseGraph
{
   unsigned int          n;                  /**< number of nodes */
   unsigned int          m;                  /**< number of edges */
   int**                 A;                  /**< adjacency list (= adjacent nodes) for each node (-1 for end of list) */
   SCIP_Real**           W;                  /**< weights for each edge */
   unsigned int*         deg;                /**< degree each node */
   unsigned int*         size;               /**< size of A/w for each node */
};

typedef struct sparseGraph SparseGraph;


/*
 * Local methods (for writing)
 */

/** initialize graph */
static
SCIP_RETCODE initGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SparseGraph*          G,                  /**< graph to free */
   unsigned int          nNodes,             /**< number of nodes */
   unsigned int          initSize            /**< initial size of lists */
   )
{
   unsigned int i;

   G->n = nNodes;
   G->m = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &G->deg, (int) nNodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &G->size, (int) nNodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &G->A, (int) nNodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &G->W, (int) nNodes) );

   for( i = 0; i < nNodes; ++i )
   {
      G->deg[i] = 0;
      G->size[i] = initSize;

      SCIP_CALL( SCIPallocBufferArray(scip, &(G->A[i]), (int) initSize) );   /*lint !e866 */
      SCIP_CALL( SCIPallocBufferArray(scip, &(G->W[i]), (int) initSize) );   /*lint !e866 */

      G->A[i][0] = -1;
   }

   return SCIP_OKAY;
}


/** frees graph */
static
void freeGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SparseGraph*          G                   /**< graph to free */
   )
{
   unsigned int i;

   for( i = 0; i < G->n; ++i )
   {
      SCIPfreeBufferArray(scip, &G->A[i]);
      SCIPfreeBufferArray(scip, &G->W[i]);
   }

   SCIPfreeBufferArray(scip, &G->W);
   SCIPfreeBufferArray(scip, &G->A);
   SCIPfreeBufferArray(scip, &G->size);
   SCIPfreeBufferArray(scip, &G->deg);
}


/** check whether there is enough capacity for one additional edge in the given adjacency list */
static
SCIP_RETCODE ensureEdgeCapacity(
   SCIP*                 scip,               /**< SCIP data structure */
   SparseGraph*          G,                  /**< graph */
   unsigned int          node                /**< list for node */
   )
{
   if( G->deg[node] + 2 > G->size[node] )
   {
      unsigned int newSize;
      newSize = G->size[node] * 2;
      SCIP_CALL( SCIPreallocBufferArray(scip, &(G->A[node]), (int) newSize) );  /*lint !e866 */
      SCIP_CALL( SCIPreallocBufferArray(scip, &(G->W[node]), (int) newSize) );  /*lint !e866 */
      G->size[node] = newSize;
   }

   return SCIP_OKAY;
}


/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n inrc/scip/reader_graph.c linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( scalars != NULL );
   assert( nvars != NULL );
   assert( constant != NULL );

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
   }
   return SCIP_OKAY;
}


/* Generate edges from given row
 *
 * We avoid parallel edges. Each row generates a clique in the graph.
 */
static
SCIP_RETCODE createEdgesFromRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of constraint variables */
   SCIP_Real*            vals,               /**< array of constraint values */
   int                   nvars,              /**< number of constraint variables */
   SparseGraph*          G                   /**< graph */
   )
{
   int i, j;
   SCIP_Real w;

   assert( scip != NULL );
   assert( nvars > 0 );

   /* compute weight */
   w = 0;
   for( i = 0; i < nvars; ++i )
      w += ABS(vals[i]);

   /* generate edges */
   for( i = 0; i < nvars; ++i )
   {
      int s;
      s = SCIPvarGetProbindex(vars[i]);
      assert( s >= 0 );

      for( j = i+1; j < nvars; ++j )
      {
         unsigned int k;
         int t;
         int a;
         t = SCIPvarGetProbindex(vars[j]);
         assert( t >= 0 );

         /* search whether edge is already present */
         k = 0;
         a = G->A[s][k];
         while( a >= 0 )
         {
            /* if we found edge, add weight */
            if( a == t )
            {
               G->W[s][k] += w;
               break;
            }
            a = G->A[s][++k];
            assert( k <= G->size[s] );
         }

         /* add new edge */
         if( a < 0 )
         {
            /* forward edge */
            SCIP_CALL( ensureEdgeCapacity(scip, G, (unsigned int) s) );
            k = G->deg[s];
            assert( G->A[s][k] == -1 );

            G->A[s][k] = t;
            G->W[s][k] = w;

            G->A[s][k+1] = -1; /*lint !e679*/
            ++G->deg[s];

            /* backward edge */
            SCIP_CALL( ensureEdgeCapacity(scip, G, (unsigned int) t) );
            k = G->deg[t];
            assert( G->A[t][k] == -1 );

            G->A[t][k] = s;
            G->W[t][k] = w;

            G->A[t][k+1] = -1; /*lint !e679*/
            ++G->deg[t];

            /* increase number of edges */
            ++G->m;
         }
      }
   }

   return SCIP_OKAY;
}


/** handle given linear constraint information */
static
SCIP_RETCODE handleLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SparseGraph*          G                   /**< graph */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( nvars > 0 );

   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );
   if( vals != NULL )
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

   /* print constraint */
   SCIP_CALL( createEdgesFromRow(scip, activevars, activevals, nactivevars, G) );

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCcg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderCcg(scip) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCcg)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPwriteCcg(scip, file, name, transformed, vars, nvars, conss, nconss, result) );

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the ccg file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   assert(reader != NULL);

   /* set non-fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCcg) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteCcg) );

   return SCIP_OKAY;
}


/** writes problem to file */
SCIP_RETCODE SCIPwriteCcg(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_VAR**            vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{  /*lint --e{715}*/
   int c;
   int v;
   int i;

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SparseGraph G;

   assert( scip != NULL );
   assert( nvars >= 0 );

   /* initialize graph */
   SCIP_CALL( initGraph(scip, &G, (unsigned int) nvars, 10) );

   /* check all constraints */
   for( c = 0; c < nconss; ++c)
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed is written only constraint are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {  
         consvars = SCIPgetVarsLinear(scip, cons);
         nconsvars = SCIPgetNVarsLinear(scip, cons);
         assert( consvars != NULL || nconsvars == 0 );

         if( nconsvars > 0 ) 
         { 
            SCIP_CALL( handleLinearCons(scip, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
                  SCIPgetNVarsLinear(scip, cons), transformed, &G) );
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);
         assert( consvars != NULL || nconsvars == 0 );

         if( nconsvars > 0 ) 
         {
            SCIP_CALL( handleLinearCons(scip, consvars, NULL, nconsvars, transformed, &G) );
         }
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {  
         consvars = SCIPgetVarsLogicor(scip, cons);
         nconsvars = SCIPgetNVarsLogicor(scip, cons);
         assert( consvars != NULL || nconsvars == 0 );

         if( nconsvars > 0 ) 
         { 
            SCIP_CALL( handleLinearCons(scip, SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons), transformed, &G) );
         }
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* w;

         consvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);
         assert( consvars != NULL || nconsvars == 0 );

         /* copy Longint array to SCIP_Real array */
         w = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
         for( v = 0; v < nconsvars; ++v )
            consvals[v] = (SCIP_Real)w[v];

         if( nconsvars > 0 ) 
         { 
            SCIP_CALL( handleLinearCons(scip, consvars, consvals, nconsvars, transformed, &G) );
         }
         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( handleLinearCons(scip, consvars, consvals, 2, transformed, &G) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else
      {
         SCIPwarningMessage(scip, "constraint handler <%s> cannot print requested format\n", conshdlrname );
         SCIPinfoMessage(scip, file, "\\ ");
         SCIP_CALL( SCIPprintCons(scip, cons, file) );
         SCIPinfoMessage(scip, file, ";\n");
      }
   }

   /* output graph */
   SCIPinfoMessage(scip, file, "c graph generated from %s\n", name);
   SCIPinfoMessage(scip, file, "p edge %d %d\n", nvars, G.m);

   for( i = 0; i < nvars; ++i )
   {
      unsigned int k;
      int a;

      k = 0;
      a = G.A[i][k];
      while( a >= 0 )
      {
         /* only output edges from lower to higher number */
         if( i < a )
         {
            /* note: node numbers start with 1 in the DIMACS format */
            SCIPinfoMessage(scip, file, "e %d %d %f\n", i+1, a+1, G.W[i][k]);
         }

         a = G.A[i][++k];
         assert( k <= G.size[i] );
      }
      assert( k == G.deg[i] );
   }

   freeGraph(scip, &G);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
