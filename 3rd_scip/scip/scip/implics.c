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

/**@file   implics.c
 * @brief  methods for implications, variable bounds, and clique tables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/event.h"
#include "scip/var.h"
#include "scip/implics.h"
#include "scip/misc.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/debug.h"
#include "scip/misc.h"

#ifndef NDEBUG
#include "scip/struct_implics.h"
#endif



/*
 * methods for variable bounds
 */

/** creates a variable bounds data structure */
static
SCIP_RETCODE vboundsCreate(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(vbounds != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, vbounds) );
   (*vbounds)->vars = NULL;
   (*vbounds)->coefs = NULL;
   (*vbounds)->constants = NULL;
   (*vbounds)->len = 0;
   (*vbounds)->size = 0;

   return SCIP_OKAY;
}

/** frees a variable bounds data structure */
void SCIPvboundsFree(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(vbounds != NULL);

   if( *vbounds != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*vbounds)->vars, (*vbounds)->size);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*vbounds)->coefs, (*vbounds)->size);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*vbounds)->constants, (*vbounds)->size);
      BMSfreeBlockMemory(blkmem, vbounds);
   }
}

/** ensures, that variable bounds arrays can store at least num entries */
static
SCIP_RETCODE vboundsEnsureSize(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(vbounds != NULL);

   /* create variable bounds data structure, if not yet existing */
   if( *vbounds == NULL )
   {
      SCIP_CALL( vboundsCreate(vbounds, blkmem) );
   }
   assert(*vbounds != NULL);
   assert((*vbounds)->len <= (*vbounds)->size);

   if( num > (*vbounds)->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*vbounds)->vars, (*vbounds)->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*vbounds)->coefs, (*vbounds)->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*vbounds)->constants, (*vbounds)->size, newsize) );
      (*vbounds)->size = newsize;
   }
   assert(num <= (*vbounds)->size);

   return SCIP_OKAY;
}

/** binary searches the insertion position of the given variable in the vbounds data structure */
static
SCIP_RETCODE vboundsSearchPos(
   SCIP_VBOUNDS*         vbounds,            /**< variable bounds data structure, or NULL */
   SCIP_VAR*             var,                /**< variable to search in vbounds data structure */
   SCIP_Bool             negativecoef,       /**< is coefficient b negative? */
   int*                  insertpos,          /**< pointer to store position where to insert new entry */
   SCIP_Bool*            found               /**< pointer to store whether the same variable was found at the returned pos */
   )
{
   SCIP_Bool exists;
   int pos;

   assert(insertpos != NULL);
   assert(found != NULL);

   /* check for empty vbounds data */
   if( vbounds == NULL )
   {
      *insertpos = 0;
      *found = FALSE;
      return SCIP_OKAY;
   }
   assert(vbounds->len >= 0);

   /* binary search for the given variable */
   exists = SCIPsortedvecFindPtr((void**)vbounds->vars, SCIPvarComp, var, vbounds->len, &pos);

   if( exists )
   {
      /* we found the variable: check if the sign of the coefficient matches */
      assert(var == vbounds->vars[pos]);
      if( (vbounds->coefs[pos] < 0.0) == negativecoef )
      {
         /* the variable exists with the same sign at the current position */
         *insertpos = pos;
         *found = TRUE;
      }
      else if( negativecoef )
      {
         assert(vbounds->coefs[pos] > 0.0);
         if( pos+1 < vbounds->len && vbounds->vars[pos+1] == var )
         {
            /* the variable exists with the desired sign at the next position */
            assert(vbounds->coefs[pos+1] < 0.0);
            *insertpos = pos+1;
            *found = TRUE;
         }
         else
         {
            /* the negative coefficient should be inserted to the right of the positive one */
            *insertpos = pos+1;
            *found = FALSE;
         }
      }
      else
      {
         assert(vbounds->coefs[pos] < 0.0);
         if( pos-1 >= 0 && vbounds->vars[pos-1] == var )
         {
            /* the variable exists with the desired sign at the previous position */
            assert(vbounds->coefs[pos-1] > 0.0);
            *insertpos = pos-1;
            *found = TRUE;
         }
         else
         {
            /* the positive coefficient should be inserted to the left of the negative one */
            *insertpos = pos;
            *found = FALSE;
         }
      }
   }
   else
   {
      *insertpos = pos;
      *found = FALSE;
   }

   return SCIP_OKAY;
}

/** adds a variable bound to the variable bounds data structure */
SCIP_RETCODE SCIPvboundsAdd(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        vboundtype,         /**< type of variable bound (LOWER or UPPER) */
   SCIP_VAR*             var,                /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             coef,               /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             constant,           /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Bool*            added               /**< pointer to store whether the variable bound was added */
   )
{
   int insertpos;
   SCIP_Bool found;

   assert(vbounds != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
   assert(added != NULL);
   assert(!SCIPsetIsZero(set, coef));

   *added = FALSE;

   /* identify insertion position of variable */
   SCIP_CALL( vboundsSearchPos(*vbounds, var, (coef < 0.0), &insertpos, &found) );
   if( found )
   {
      /* the same variable already exists in the vbounds data structure: use the better vbound */
      assert(*vbounds != NULL);
      assert(0 <= insertpos && insertpos < (*vbounds)->len);
      assert((*vbounds)->vars[insertpos] == var);
      assert(((*vbounds)->coefs[insertpos] < 0.0) == (coef < 0.0));

      if( vboundtype == SCIP_BOUNDTYPE_UPPER )
      {
         if( constant + MIN(coef, 0.0) < (*vbounds)->constants[insertpos] + MIN((*vbounds)->coefs[insertpos], 0.0) )
         {
            (*vbounds)->coefs[insertpos] = coef;
            (*vbounds)->constants[insertpos] = constant;
            *added = TRUE;
         }
      }
      else
      {
         if( constant + MAX(coef, 0.0) > (*vbounds)->constants[insertpos] + MAX((*vbounds)->coefs[insertpos], 0.0) )
         {
            (*vbounds)->coefs[insertpos] = coef;
            (*vbounds)->constants[insertpos] = constant;
            *added = TRUE;
         }
      }
   }
   else
   {
      int i;

      /* the given variable does not yet exist in the vbounds */
      SCIP_CALL( vboundsEnsureSize(vbounds, blkmem, set, *vbounds != NULL ? (*vbounds)->len+1 : 1) );
      assert(*vbounds != NULL);
      assert(0 <= insertpos && insertpos <= (*vbounds)->len);
      assert(0 <= insertpos && insertpos < (*vbounds)->size);

      /* insert variable at the correct position */
      for( i = (*vbounds)->len; i > insertpos; --i )
      {
         assert(!SCIPsetIsZero(set, (*vbounds)->coefs[i-1]));
         (*vbounds)->vars[i] = (*vbounds)->vars[i-1];
         (*vbounds)->coefs[i] = (*vbounds)->coefs[i-1];
         (*vbounds)->constants[i] = (*vbounds)->constants[i-1];
      }
      assert(!SCIPsetIsZero(set, coef));
      (*vbounds)->vars[insertpos] = var;
      (*vbounds)->coefs[insertpos] = coef;
      (*vbounds)->constants[insertpos] = constant;
      (*vbounds)->len++;
      *added = TRUE;
   }

   return SCIP_OKAY;
}

/** removes from variable x a variable bound x >=/<= b*z + d with binary or integer z */
SCIP_RETCODE SCIPvboundsDel(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             vbdvar,             /**< variable z    in x >=/<= b*z + d */
   SCIP_Bool             negativecoef        /**< is coefficient b negative? */
   )
{
   SCIP_Bool found;
   int pos;
   int i;

   assert(vbounds != NULL);
   assert(*vbounds != NULL);

   /* searches for variable z in variable bounds of x */
   SCIP_CALL( vboundsSearchPos(*vbounds, vbdvar, negativecoef, &pos, &found) );
   if( !found )
      return SCIP_OKAY;

   assert(0 <= pos && pos < (*vbounds)->len);
   assert((*vbounds)->vars[pos] == vbdvar);
   assert(((*vbounds)->coefs[pos] < 0.0) == negativecoef);

   /* removes z from variable bounds of x */
   for( i = pos; i < (*vbounds)->len - 1; i++ )
   {
      (*vbounds)->vars[i] = (*vbounds)->vars[i+1];
      (*vbounds)->coefs[i] = (*vbounds)->coefs[i+1];
      (*vbounds)->constants[i] = (*vbounds)->constants[i+1];
   }
   (*vbounds)->len--;

#ifndef NDEBUG
   SCIP_CALL( vboundsSearchPos(*vbounds, vbdvar, negativecoef, &pos, &found) );
   assert(!found);
#endif

   /* free vbounds data structure, if it is empty */
   if( (*vbounds)->len == 0 )
      SCIPvboundsFree(vbounds, blkmem);

   return SCIP_OKAY;
}

/** reduces the number of variable bounds stored in the given variable bounds data structure */
void SCIPvboundsShrink(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   newnvbds            /**< new number of variable bounds */
   )
{
   assert(vbounds != NULL);
   assert(*vbounds != NULL);
   assert(newnvbds <= (*vbounds)->len);

   if( newnvbds == 0 )
      SCIPvboundsFree(vbounds, blkmem);
   else
      (*vbounds)->len = newnvbds;
}




/*
 * methods for implications
 */

#ifndef NDEBUG
/** comparator function for implication variables in the implication data structure */
static
SCIP_DECL_SORTPTRCOMP(compVars)
{  /*lint --e{715}*/
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int var1idx;
   int var2idx;

   var1 = (SCIP_VAR*)elem1;
   var2 = (SCIP_VAR*)elem2;
   assert(var1 != NULL);
   assert(var2 != NULL);
   var1idx = SCIPvarGetIndex(var1);
   var2idx = SCIPvarGetIndex(var2);

   if( var1idx < var2idx )
      return -1;
   else if( var1idx > var2idx )
      return +1;
   else
      return 0;
}

/** performs integrity check on implications data structure */
static
void checkImplics(
   SCIP_IMPLICS*         implics             /**< implications data structure */
   )
{
   SCIP_Bool varfixing;

   if( implics == NULL )
      return;

   varfixing = FALSE;
   do
   {
      SCIP_VAR** vars;
      SCIP_BOUNDTYPE* types;
      int nimpls;
      int i;

      vars = implics->vars[varfixing];
      types = implics->types[varfixing];
      nimpls = implics->nimpls[varfixing];

      assert(0 <= nimpls && nimpls <= implics->size[varfixing]);
      for( i = 1; i < nimpls; ++i )
      {
         int cmp;

         cmp = compVars(vars[i-1], vars[i]);
         assert(cmp <= 0);
         assert((cmp == 0) == (vars[i-1] == vars[i]));
         assert(cmp < 0 || (types[i-1] == SCIP_BOUNDTYPE_LOWER && types[i] == SCIP_BOUNDTYPE_UPPER));
      }

      varfixing = !varfixing;
   }
   while( varfixing == TRUE );
}
#else
#define checkImplics(implics) /**/
#endif

/** creates an implications data structure */
static
SCIP_RETCODE implicsCreate(
   SCIP_IMPLICS**        implics,            /**< pointer to store implications data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(implics != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, implics) );

   (*implics)->vars[0] = NULL;
   (*implics)->types[0] = NULL;
   (*implics)->bounds[0] = NULL;
   (*implics)->ids[0] = NULL;
   (*implics)->size[0] = 0;
   (*implics)->nimpls[0] = 0;
   (*implics)->vars[1] = NULL;
   (*implics)->types[1] = NULL;
   (*implics)->bounds[1] = NULL;
   (*implics)->ids[1] = NULL;
   (*implics)->size[1] = 0;
   (*implics)->nimpls[1] = 0;

   return SCIP_OKAY;
}

/** frees an implications data structure */
void SCIPimplicsFree(
   SCIP_IMPLICS**        implics,            /**< pointer of implications data structure to free */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(implics != NULL);

   if( *implics != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->vars[0], (*implics)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->types[0], (*implics)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->bounds[0], (*implics)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->ids[0], (*implics)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->vars[1], (*implics)->size[1]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->types[1], (*implics)->size[1]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->bounds[1], (*implics)->size[1]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->ids[1], (*implics)->size[1]);
      BMSfreeBlockMemory(blkmem, implics);
   }
}

/** ensures, that arrays for x == 0 or x == 1 in implications data structure can store at least num entries */
static
SCIP_RETCODE implicsEnsureSize(
   SCIP_IMPLICS**        implics,            /**< pointer to implications data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             varfixing,          /**< FALSE if size of arrays for x == 0 has to be ensured, TRUE for x == 1 */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(implics != NULL);

   /* create implications data structure, if not yet existing */
   if( *implics == NULL )
   {
      SCIP_CALL( implicsCreate(implics, blkmem) );
   }
   assert(*implics != NULL);
   assert((*implics)->nimpls[varfixing] <= (*implics)->size[varfixing]);

   if( num > (*implics)->size[varfixing] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*implics)->vars[varfixing], (*implics)->size[varfixing], newsize) ); /*lint !e866*/
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*implics)->types[varfixing], (*implics)->size[varfixing], newsize) ); /*lint !e866*/
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*implics)->bounds[varfixing], (*implics)->size[varfixing], newsize) ); /*lint !e866*/
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*implics)->ids[varfixing], (*implics)->size[varfixing], newsize) ); /*lint !e866*/
      (*implics)->size[varfixing] = newsize;
   }
   assert(num <= (*implics)->size[varfixing]);

   return SCIP_OKAY;
}

/** gets the positions of the implications y >= l and y <= u in the implications data structure;
 *  if no lower or upper bound implication for y was found, -1 is returned
 */
static
void implicsSearchVar(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y is searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   int*                  poslower,           /**< pointer to store position of y_lower (-1 if not found) */
   int*                  posupper,           /**< pointer to store position of y_upper (-1 if not found) */
   int*                  posadd              /**< pointer to store position of first y entry, or where a new y entry
                                              *   should be placed */
   )
{
   SCIP_Bool found;
   int right;
   int pos;

   assert(implics != NULL);
   assert(poslower != NULL);
   assert(posupper != NULL);
   assert(posadd != NULL);

   if( implics->nimpls[varfixing] == 0 )
   {
      /* there are no implications with non-binary variable y */
      *posadd = 0;
      *poslower = -1;
      *posupper = -1;
      return;
   }
   right = implics->nimpls[varfixing] - 1;
   assert(0 <= right);

   /* search for the position in the sorted array (via binary search) */
   found = SCIPsortedvecFindPtr((void**)(&(implics->vars[varfixing][0])), SCIPvarComp, (void*)implvar, right + 1, &pos);

   if( !found )
   {
      /* y was not found */
      assert(pos >= right || compVars((void*)implics->vars[varfixing][pos], (void*)implvar) > 0);
      assert(pos == 0 || compVars((void*)implics->vars[varfixing][pos-1], (void*)implvar) < 0);
      *poslower = -1;
      *posupper = -1;
      *posadd = pos;
   }
   else
   {
      /* y was found */
      assert(implvar == implics->vars[varfixing][pos]);

      /* set poslower and posupper */
      if( implics->types[varfixing][pos] == SCIP_BOUNDTYPE_LOWER )
      {
         /* y was found as y_lower (on position middle) */
         *poslower = pos;
         if( pos + 1 < implics->nimpls[varfixing] && implics->vars[varfixing][pos+1] == implvar )
         {  
            assert(implics->types[varfixing][pos+1] == SCIP_BOUNDTYPE_UPPER);
            *posupper = pos + 1;
         }
         else
            *posupper = -1;
         *posadd = pos;
      }
      else
      {
         /* y was found as y_upper (on position pos) */
         *posupper = pos;
         if( pos - 1 >= 0 && implics->vars[varfixing][pos-1] == implvar )
         {  
            assert(implics->types[varfixing][pos-1] == SCIP_BOUNDTYPE_LOWER);
            *poslower = pos - 1;
            *posadd = pos - 1;
         }
         else
         {
            *poslower = -1;
            *posadd = pos;
         }
      }
   }
}

/** returns whether variable y is already contained in implications for x == 0 or x == 1 with the given impltype
 *  y can be contained in structure with y >= b (y_lower) and y <= b (y_upper) 
 */
static
SCIP_Bool implicsSearchImplic(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y is searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_BOUNDTYPE        impltype,           /**< type of implication y <=/>= b to search for */
   int*                  poslower,           /**< pointer to store position of y_lower (inf if not found) */
   int*                  posupper,           /**< pointer to store position of y_upper (inf if not found) */
   int*                  posadd              /**< pointer to store correct position (with respect to impltype) to add y */
   )
{
   assert(implics != NULL);
   assert(poslower != NULL);
   assert(posupper != NULL);
   assert(posadd != NULL);

   implicsSearchVar(implics, varfixing, implvar, poslower, posupper, posadd);
   assert(*poslower == -1 || *posupper == -1 || *posupper == (*poslower)+1);
   assert(*poslower == -1 || *posadd == *poslower);
   assert(*poslower >= 0 || *posupper == -1 || *posadd == *posupper);
   assert(0 <= *posadd && *posadd <= implics->nimpls[varfixing]);

   if( impltype == SCIP_BOUNDTYPE_LOWER )
      return (*poslower >= 0);
   else
   {
      if( *poslower >= 0 )
      {
         assert(*posadd == *poslower);
         (*posadd)++;
      }
      return (*posupper >= 0);
   }
}

/** adds an implication x == 0/1 -> y <= b or y >= b to the implications data structure;
 *  the implication must be non-redundant
 */
SCIP_RETCODE SCIPimplicsAdd(
   SCIP_IMPLICS**        implics,            /**< pointer to implications data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             varfixing,          /**< FALSE if implication for x == 0 has to be added, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool             isshortcut,         /**< is the implication a shortcut, i.e., added as part of the transitive closure of another implication? */
   SCIP_Bool*            conflict,           /**< pointer to store whether implication causes a conflict for variable x */
   SCIP_Bool*            added               /**< pointer to store whether the implication was added */
   )
{
   int poslower;
   int posupper;
   int posadd;
   SCIP_Bool found;
#ifndef NDEBUG
   int k;
#endif

   assert(implics != NULL);
   assert(*implics == NULL || 0 <= (*implics)->nimpls[varfixing]);
   assert(stat != NULL);
   assert(SCIPvarIsActive(implvar));
   assert(SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_LOOSE); 
   assert((impltype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasGT(set, implbound, SCIPvarGetLbGlobal(implvar)))
      || (impltype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasLT(set, implbound, SCIPvarGetUbGlobal(implvar))));
   assert(conflict != NULL);
   assert(added != NULL);

   SCIPsetDebugMsg(set, "adding implication to implics %p [%u]: <%s> %s %g\n",
      (void*)*implics, varfixing, SCIPvarGetName(implvar), impltype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", implbound);

   checkImplics(*implics);

   *conflict = FALSE;
   *added = FALSE;

   /* check if variable is already contained in implications data structure */
   if( *implics != NULL )
   {
      found = implicsSearchImplic(*implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
      assert(-1 <= poslower && poslower < (*implics)->nimpls[varfixing]);
      assert(-1 <= posupper && posupper < (*implics)->nimpls[varfixing]);
      assert(0 <= posadd && posadd <= (*implics)->nimpls[varfixing]);
      assert(poslower == -1 || (*implics)->types[varfixing][poslower] == SCIP_BOUNDTYPE_LOWER);
      assert(posupper == -1 || (*implics)->types[varfixing][posupper] == SCIP_BOUNDTYPE_UPPER);
   }
   else
   {
      found = FALSE;
      poslower = -1;
      posupper = -1;
      posadd = 0;
   }

   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      assert(found == (poslower >= 0));

      /* check if y >= b is redundant */
      if( poslower >= 0 && SCIPsetIsFeasLE(set, implbound, (*implics)->bounds[varfixing][poslower]) )
      {
         SCIPsetDebugMsg(set, " -> implication is redundant to <%s> >= %g\n", SCIPvarGetName(implvar), (*implics)->bounds[varfixing][poslower]);
         return SCIP_OKAY;
      }

      /* check if y >= b causes conflict for x (i.e. y <= a (with a < b) is also valid) */
      if( posupper >= 0 && SCIPsetIsFeasGT(set, implbound, (*implics)->bounds[varfixing][posupper]) )
      {
         SCIPsetDebugMsg(set, " -> implication is conflicting to <%s> <= %g\n", SCIPvarGetName(implvar), (*implics)->bounds[varfixing][posupper]);
         *conflict = TRUE;
         return SCIP_OKAY;
      }

      *added = TRUE;

      /* check if entry of the same type already exists */
      if( found )
      {
         assert(poslower >= 0);
         assert(posadd == poslower);

         /* add y >= b by changing old entry on poslower */
         assert((*implics)->vars[varfixing][poslower] == implvar);
         assert(SCIPsetIsFeasGT(set, implbound, (*implics)->bounds[varfixing][poslower]));
         (*implics)->bounds[varfixing][poslower] = implbound;
      }
      else
      {
         assert(poslower == -1);

         /* add y >= b by creating a new entry on posadd */
         SCIP_CALL( implicsEnsureSize(implics, blkmem, set, varfixing,
               *implics != NULL ? (*implics)->nimpls[varfixing]+1 : 1) );
         assert(*implics != NULL);

	 if( (*implics)->nimpls[varfixing] - posadd > 0 )
	 {
	    int amount = ((*implics)->nimpls[varfixing] - posadd);

#ifndef NDEBUG
	    for( k = (*implics)->nimpls[varfixing]; k > posadd; k-- )
	    {
	       assert(compVars((void*)(*implics)->vars[varfixing][k-1], (void*)implvar) >= 0);
	    }
#endif
	    BMSmoveMemoryArray(&((*implics)->types[varfixing][posadd+1]), &((*implics)->types[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->ids[varfixing][posadd+1]), &((*implics)->ids[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->vars[varfixing][posadd+1]), &((*implics)->vars[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->bounds[varfixing][posadd+1]), &((*implics)->bounds[varfixing][posadd]), amount); /*lint !e866*/
	 }

         (*implics)->vars[varfixing][posadd] = implvar;
         (*implics)->types[varfixing][posadd] = impltype;
         (*implics)->bounds[varfixing][posadd] = implbound;
         (*implics)->ids[varfixing][posadd] = (isshortcut ? -stat->nimplications : stat->nimplications);
         (*implics)->nimpls[varfixing]++;
#ifndef NDEBUG
         for( k = posadd-1; k >= 0; k-- )
            assert(compVars((void*)(*implics)->vars[varfixing][k], (void*)implvar) <= 0);
#endif
         stat->nimplications++;
      }
   }
   else
   {
      assert(found == (posupper >= 0));

      /* check if y <= b is redundant */
      if( posupper >= 0 && SCIPsetIsFeasGE(set, implbound, (*implics)->bounds[varfixing][posupper]) )
      {
         SCIPsetDebugMsg(set, " -> implication is redundant to <%s> <= %g\n", SCIPvarGetName(implvar), (*implics)->bounds[varfixing][posupper]);
         return SCIP_OKAY;
      }

      /* check if y <= b causes conflict for x (i.e. y >= a (with a > b) is also valid) */
      if( poslower >= 0 && SCIPsetIsFeasLT(set, implbound, (*implics)->bounds[varfixing][poslower]) )
      {
         SCIPsetDebugMsg(set, " -> implication is conflicting to <%s> >= %g\n", SCIPvarGetName(implvar), (*implics)->bounds[varfixing][poslower]);
         *conflict = TRUE;
         return SCIP_OKAY;
      }

      *added = TRUE;

      /* check if entry of the same type already exists */
      if( found )
      {
         assert(posupper >= 0);
         assert(posadd == posupper);

         /* add y <= b by changing old entry on posupper */
         assert((*implics)->vars[varfixing][posupper] == implvar);
         assert(SCIPsetIsFeasLT(set, implbound,(*implics)->bounds[varfixing][posupper]));
         (*implics)->bounds[varfixing][posupper] = implbound;
      }
      else
      {
         /* add y <= b by creating a new entry on posadd */
         assert(posupper == -1);

         SCIP_CALL( implicsEnsureSize(implics, blkmem, set, varfixing,
               *implics != NULL ? (*implics)->nimpls[varfixing]+1 : 1) );
         assert(*implics != NULL);

	 if( (*implics)->nimpls[varfixing] - posadd > 0 )
	 {
	    int amount = ((*implics)->nimpls[varfixing] - posadd);

#ifndef NDEBUG
	    for( k = (*implics)->nimpls[varfixing]; k > posadd; k-- )
	    {
	       assert(compVars((void*)(*implics)->vars[varfixing][k-1], (void*)implvar) >= 0);
	    }
#endif
	    BMSmoveMemoryArray(&((*implics)->types[varfixing][posadd+1]), &((*implics)->types[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->ids[varfixing][posadd+1]), &((*implics)->ids[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->vars[varfixing][posadd+1]), &((*implics)->vars[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->bounds[varfixing][posadd+1]), &((*implics)->bounds[varfixing][posadd]), amount); /*lint !e866*/
	 }

         (*implics)->vars[varfixing][posadd] = implvar;
         (*implics)->types[varfixing][posadd] = impltype;
         (*implics)->bounds[varfixing][posadd] = implbound;
         (*implics)->ids[varfixing][posadd] = (isshortcut ? -stat->nimplications : stat->nimplications);
         (*implics)->nimpls[varfixing]++;
#ifndef NDEBUG
         for( k = posadd-1; k >= 0; k-- )
            assert(compVars((void*)(*implics)->vars[varfixing][k], (void*)implvar) <= 0);
#endif
         stat->nimplications++;
      }
   }

   checkImplics(*implics);

   return SCIP_OKAY;
}

/** removes the implication  x <= 0 or x >= 1  ==>  y <= b  or  y >= b  from the implications data structure */
SCIP_RETCODE SCIPimplicsDel(
   SCIP_IMPLICS**        implics,            /**< pointer to implications data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             varfixing,          /**< FALSE if y should be removed from implications for x <= 0, TRUE for x >= 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype            /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   )
{  /*lint --e{715}*/
   int poslower;
   int posupper; 
   int posadd;
   SCIP_Bool found;

   assert(implics != NULL);
   assert(*implics != NULL);
   assert(implvar != NULL);

   SCIPsetDebugMsg(set, "deleting implication from implics %p [%u]: <%s> %s x\n",
      (void*)*implics, varfixing, SCIPvarGetName(implvar), impltype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=");

   checkImplics(*implics);

   /* searches for y in implications of x */
   found = implicsSearchImplic(*implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
   if( !found )
   {
      SCIPsetDebugMsg(set, " -> implication was not found\n");
      return SCIP_OKAY;
   }

   assert((impltype == SCIP_BOUNDTYPE_LOWER && poslower >= 0 && posadd == poslower) 
      || (impltype == SCIP_BOUNDTYPE_UPPER && posupper >= 0 && posadd == posupper));
   assert(0 <= posadd && posadd < (*implics)->nimpls[varfixing]);
   assert((*implics)->vars[varfixing][posadd] == implvar);
   assert((*implics)->types[varfixing][posadd] == impltype);

   /* removes y from implications of x */
   if( (*implics)->nimpls[varfixing] - posadd > 1 )
   {
      int amount = ((*implics)->nimpls[varfixing] - posadd - 1);

      BMSmoveMemoryArray(&((*implics)->types[varfixing][posadd]), &((*implics)->types[varfixing][posadd+1]), amount); /*lint !e866*/
      BMSmoveMemoryArray(&((*implics)->vars[varfixing][posadd]), &((*implics)->vars[varfixing][posadd+1]), amount); /*lint !e866*/
      BMSmoveMemoryArray(&((*implics)->bounds[varfixing][posadd]), &((*implics)->bounds[varfixing][posadd+1]), amount); /*lint !e866*/
   }
   (*implics)->nimpls[varfixing]--;

   /* free implics data structure, if it is empty */
   if( (*implics)->nimpls[0] == 0 && (*implics)->nimpls[1] == 0 )
      SCIPimplicsFree(implics, blkmem);

   checkImplics(*implics);

   return SCIP_OKAY;
}

/** returns which implications on given variable y are contained in implications for x == 0 or x == 1 */
void SCIPimplicsGetVarImplics(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_Bool*            haslowerimplic,     /**< pointer to store whether there exists an implication y >= l */
   SCIP_Bool*            hasupperimplic      /**< pointer to store whether there exists an implication y <= u */
   )
{
   int poslower;
   int posupper;
   int posadd;

   assert(haslowerimplic != NULL);
   assert(hasupperimplic != NULL);

   implicsSearchVar(implics, varfixing, implvar, &poslower, &posupper, &posadd);

   *haslowerimplic = (poslower >= 0);
   *hasupperimplic = (posupper >= 0);
}

/** returns whether an implication y <= b or y >= b is contained in implications for x == 0 or x == 1 */
SCIP_Bool SCIPimplicsContainsImpl(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_BOUNDTYPE        impltype            /**< type of implication y <=/>= b to search for */
   )
{
   int poslower;
   int posupper;
   int posadd;

   return implicsSearchImplic(implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
}




/*
 * methods for cliques
 */

/* swaps cliques at positions first and second in cliques array of clique table */
static
void cliquetableSwapCliques(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   int                   first,              /**< first index */
   int                   second              /**< second index */
   )
{
   /* nothing to do if clique happens to be at the pole position of clean cliques */
   if( first != second )
   {
      SCIP_CLIQUE* tmp;

      tmp = cliquetable->cliques[first];
      assert(tmp->index == first);
      assert(cliquetable->cliques[second]->index == second);

      cliquetable->cliques[first] = cliquetable->cliques[second];
      cliquetable->cliques[second] = tmp;

      /* change the indices accordingly */
      tmp->index = second;
      cliquetable->cliques[first]->index = first;
   }
}

/* moves clique to the last position of currently dirty cliques */
static
void cliquetableMarkCliqueForCleanup(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(cliquetable->ndirtycliques <= clique->index);
   assert(cliquetable->cliques[clique->index] == clique);
   assert(cliquetable->ndirtycliques < cliquetable->ncliques);
   assert(SCIPcliqueIsCleanedUp(cliquetable->cliques[cliquetable->ndirtycliques]));

   /* nothing to do if clique happens to be at the pole position of clean cliques */
   if( clique->index > cliquetable->ndirtycliques )
   {
      cliquetableSwapCliques(cliquetable, clique->index, cliquetable->ndirtycliques);
   }

   ++cliquetable->ndirtycliques;
}


/** creates a clique data structure with already created variables and values arrays in the size of 'size' */
static
SCIP_RETCODE cliqueCreateWithData(
   SCIP_CLIQUE**         clique,             /**< pointer to store clique data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   size,               /**< initial size of clique */
   SCIP_VAR**            vars,               /**< binary variables in the clique: at most one can be set to the given
                                              *   value */
   SCIP_Bool*            values,             /**< values of the variables in the clique */
   int                   nvars,              /**< number of variables in the clique */
   int                   id,                 /**< unique identifier of the clique */
   SCIP_Bool             isequation          /**< is the clique an equation or an inequality? */
   )
{
   assert(clique != NULL);
   assert(blkmem != NULL);
   assert(size >= nvars && nvars > 0);
   assert(vars != NULL);
   assert(values != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, clique) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*clique)->vars, vars, size) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*clique)->values, values, size) );
   (*clique)->nvars = nvars;
   (*clique)->size = size;
   (*clique)->startcleanup = -1;
   (*clique)->id = (unsigned int)id;
   (*clique)->eventsissued = FALSE;
   (*clique)->equation = isequation;
   (*clique)->index = -1;

   return SCIP_OKAY;
}

/** frees a clique data structure */
static
void cliqueFree(
   SCIP_CLIQUE**         clique,             /**< pointer to store clique data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(clique != NULL);

   if( *clique != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*clique)->vars, (*clique)->size);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*clique)->values, (*clique)->size);
      BMSfreeBlockMemory(blkmem, clique);
   }
}

/** ensures, that clique arrays can store at least num entries */
static
SCIP_RETCODE cliqueEnsureSize(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(clique != NULL);

   if( num > clique->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &clique->vars, clique->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &clique->values, clique->size, newsize) );
      clique->size = newsize;
   }
   assert(num <= clique->size);

   return SCIP_OKAY;
}

/** returns the position of the given variable/value pair in the clique; returns -1 if variable/value pair is not member
 *  of the clique
 */
int SCIPcliqueSearchVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_VAR*             var,                /**< variable to search for */
   SCIP_Bool             value               /**< value of the variable in the clique */
   )
{
   int varidx;
   int left;
   int right;

   assert(clique != NULL);

   varidx = SCIPvarGetIndex(var);
   left = -1;
   right = clique->nvars;
   while( left < right-1 )
   {
      int middle;
      int idx;

      middle = (left+right)/2;
      idx = SCIPvarGetIndex(clique->vars[middle]);
      assert(idx >= 0);
      if( varidx < idx )
         right = middle;
      else if( varidx > idx )
         left = middle;
      else
      {
         assert(var == clique->vars[middle]);

         /* now watch out for the correct value */
         if( clique->values[middle] < value )
         {
            int i;
            for( i = middle+1; i < clique->nvars && clique->vars[i] == var; ++i )
            {
               if( clique->values[i] == value )
                  return i;
            }
            return -1;
         }
         if( clique->values[middle] > value )
         {
            int i;
            for( i = middle-1; i >= 0 && clique->vars[i] == var; --i )
            {
               if( clique->values[i] == value )
                  return i;
            }
            return -1;
         }
         return middle;
      }
   }

   return -1;
}

/** returns whether the given variable/value pair is member of the given clique */
SCIP_Bool SCIPcliqueHasVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_VAR*             var,                /**< variable to remove from the clique */
   SCIP_Bool             value               /**< value of the variable in the clique */
   )
{
   return (SCIPcliqueSearchVar(clique, var, value) >= 0);
}

/** adds a single variable to the given clique */
SCIP_RETCODE SCIPcliqueAddVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to add to the clique */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_Bool*            doubleentry,        /**< pointer to store whether the variable and value occurs twice in the clique */
   SCIP_Bool*            oppositeentry       /**< pointer to store whether the variable with opposite value is in the clique */
   )
{
   int pos;
   int i;

   assert(clique != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarIsBinary(var));
   assert(doubleentry != NULL);
   assert(oppositeentry != NULL);

   SCIPsetDebugMsg(set, "adding variable <%s> == %u to clique %u\n", SCIPvarGetName(var), value, clique->id);

   *doubleentry = FALSE;
   *oppositeentry = FALSE;

   /* allocate memory */
   SCIP_CALL( cliqueEnsureSize(clique, blkmem, set, clique->nvars+1) );

   /* search for insertion position by binary variable, note that first the entries are order after variable index and
    * second after the bool value of the corresponding variable
    */
   (void) SCIPsortedvecFindPtr((void**) clique->vars, SCIPvarComp, var, clique->nvars, &pos);

   assert(pos >= 0 && pos <= clique->nvars);
   /* remember insertion position for later, pos might change */
   i = pos;

   if( pos < clique->nvars )
   {
      const int amount = clique->nvars - pos;

      /* moving elements to correct position */
      BMSmoveMemoryArray(&(clique->vars[pos+1]), &(clique->vars[pos]), amount); /*lint !e866*/
      BMSmoveMemoryArray(&(clique->values[pos+1]), &(clique->values[pos]), amount); /*lint !e866*/
      clique->nvars++;

      /* insertion for a variable with cliquevalue FALSE */
      if( !value )
      {
	 /* find last entry with the same variable and value behind the insertion position */
	 for( ; pos < clique->nvars - 1 && clique->vars[pos + 1] == var && clique->values[pos + 1] == value; ++pos ); /*lint !e722*/

	 /* check if the same variable with other value also exists */
	 if( pos < clique->nvars - 1 && clique->vars[pos + 1] == var )
	 {
	    assert(clique->values[pos + 1] != value);
	    *oppositeentry = TRUE;
	 }

	 /* check if we found the same variable with the same value more than once */
	 if( pos != i )
	    *doubleentry = TRUE;
	 else
	 {
	    /* find last entry with the same variable and different value before the insertion position */
	    for( ; pos > 0 && clique->vars[pos - 1] == var && clique->values[pos - 1] != value; --pos ); /*lint !e722*/

	    /* check if we found the same variable with the same value more than once */
	    if( pos > 0 && clique->vars[pos - 1] == var )
	    {
	       assert(clique->values[pos - 1] == value);

	       *doubleentry = TRUE;
	    }
	    /* if we found the same variable with different value, we need to order them correctly */
	    if( pos != i )
	    {
	       assert(clique->vars[pos] == var);
	       assert(clique->values[pos] != value);

	       clique->values[pos] = value;
	       value = !value;
	    }
	 }
      }
      /* insertion for a variable with cliquevalue TRUE */
      else
      {
	 /* find last entry with the same variable and different value behind the insertion position */
	 for( ; pos < clique->nvars - 1 && clique->vars[pos + 1] == var && clique->values[pos + 1] != value; ++pos ); /*lint !e722*/

	 /* check if the same variable with other value also exists */
	 if( pos < clique->nvars - 1 && clique->vars[pos + 1] == var )
	 {
	    assert(clique->values[pos + 1] == value);
	    *doubleentry = TRUE;
	 }

	 /* check if we found the same variable with different value */
	 if( pos != i )
	 {
	    *oppositeentry = TRUE;

	    /* if we found the same variable with different value, we need to order them correctly */
	    assert(clique->vars[pos] == var);
	    assert(clique->values[pos] != value);

	    clique->values[pos] = value;
	    value = !value;
	 }
	 else
	 {
	    /* find last entry with the same variable and value before the insertion position */
	    for( ; pos > 0 && clique->vars[pos - 1] == var && clique->values[pos - 1] == value; --pos ); /*lint !e722*/

	    if( pos != i )
	       *doubleentry = TRUE;

	    /* check if we found the same variable with different value up front */
	    if( pos > 0 && clique->vars[pos - 1] == var && clique->values[pos - 1] != value )
	       *oppositeentry = TRUE;
	 }
      }
   }
   else
      clique->nvars++;

   clique->vars[i] = var;
   clique->values[i] = value;
   clique->eventsissued = FALSE;

   return SCIP_OKAY;
}

/** removes a single variable from the given clique */
void SCIPcliqueDelVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< variable to remove from the clique */
   SCIP_Bool             value               /**< value of the variable in the clique */
   )
{
   int pos;

   assert(clique != NULL);
   assert(SCIPvarIsBinary(var));
   assert(cliquetable != NULL);

   /* if the clique is the leading clique during the cleanup step, we do not need to insert it again */
   if( cliquetable->incleanup && clique->index == 0 )
      return;

   SCIPdebugMessage("marking variable <%s> == %u from clique %u for deletion\n", SCIPvarGetName(var), value, clique->id);

   /* find variable in clique */
   pos = SCIPcliqueSearchVar(clique, var, value);

   assert(0 <= pos && pos < clique->nvars);
   assert(clique->vars[pos] == var);
   assert(clique->values[pos] == value);

   /* inform the clique table that this clique should be cleaned up */
   if( clique->startcleanup == -1 )
      cliquetableMarkCliqueForCleanup(cliquetable, clique);

   if( clique->startcleanup == -1 || pos < clique->startcleanup )
      clique->startcleanup = pos;

#ifdef SCIP_MORE_DEBUG
   {
      int v;
      /* all variables prior to the one marked for startcleanup should be unfixed */
      for( v = clique->startcleanup - 1; v >= 0; --v )
      {
         assert(SCIPvarGetUbGlobal(clique->vars[v]) > 0.5);
         assert(SCIPvarGetLbGlobal(clique->vars[v]) < 0.5);
      }
   }
#endif
}

/** gets the position of the given clique in the cliques array; returns -1 if clique is not member of cliques array */
static
int cliquesSearchClique(
   SCIP_CLIQUE**         cliques,            /**< array of cliques */
   int                   ncliques,           /**< number of cliques in the cliques array */
   SCIP_CLIQUE*          clique              /**< clique to search for */
   )
{
   unsigned int cliqueid;
   int left;
   int right;

   assert(cliques != NULL || ncliques == 0);
   assert(clique != NULL);

   cliqueid = clique->id; /*lint !e732*/
   left = -1;
   right = ncliques;
   while( left < right-1 )
   {
      unsigned int id;
      int middle;

      assert(cliques != NULL);
      middle = (left+right)/2;
      id = cliques[middle]->id; /*lint !e732*/
      if( cliqueid < id )
         right = middle;
      else if( cliqueid > id )
         left = middle;
      else
      {
         assert(clique == cliques[middle]);
         return middle;
      }
   }

   return -1;
}

#ifdef SCIP_MORE_DEBUG
/** checks whether clique appears in all clique lists of the involved variables */
static
void cliqueCheck(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   int i;

   assert(clique != NULL);

   for( i = 0; i < clique->nvars; ++i )
   {
      SCIP_CLIQUE** cliques;
      int ncliques;
      int pos;
      SCIP_VAR* var;

      var = clique->vars[i];
      assert(i == 0 || SCIPvarGetIndex(clique->vars[i-1]) <= SCIPvarGetIndex(var));
      assert(i == 0 || clique->vars[i-1] != var || clique->values[i-1] <= clique->values[i]);
      ncliques = SCIPvarGetNCliques(var, clique->values[i]);

      assert(SCIPvarIsActive(var) || ncliques == 0);

      /* cliquelist of inactive variables are already destroyed */
      if( ncliques == 0 )
         continue;

      cliques = SCIPvarGetCliques(var, clique->values[i]);
      pos = cliquesSearchClique(cliques, ncliques, clique);

      /* assert that the clique is correctly listed in all clique lists of unfixed variables. For fixed variables,
       * we require that a clean up has been correctly scheduled, but not yet been processed
       */
      if( SCIPvarGetUbGlobal(var) - SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(0 <= pos && pos < ncliques);
         assert(cliques[pos] == clique);
      }
      else
         assert(0 <= clique->startcleanup && clique->startcleanup <= i);
   }
   assert(clique->index >= 0);
}
#else
#define cliqueCheck(clique) /**/
#endif

/** creates a clique list data structure */
static
SCIP_RETCODE cliquelistCreate(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to store clique list data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(cliquelist != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, cliquelist) );
   (*cliquelist)->cliques[0] = NULL;
   (*cliquelist)->cliques[1] = NULL;
   (*cliquelist)->ncliques[0] = 0;
   (*cliquelist)->ncliques[1] = 0;
   (*cliquelist)->size[0] = 0;
   (*cliquelist)->size[1] = 0;

   return SCIP_OKAY;
}

/** frees a clique list data structure */
void SCIPcliquelistFree(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(cliquelist != NULL);

   if( *cliquelist != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*cliquelist)->cliques[0], (*cliquelist)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*cliquelist)->cliques[1], (*cliquelist)->size[1]);
      BMSfreeBlockMemory(blkmem, cliquelist);
   }
}

/** ensures, that clique list arrays can store at least num entries */
static
SCIP_RETCODE cliquelistEnsureSize(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             value,              /**< value of the variable for which the clique list should be extended */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(cliquelist != NULL);

   if( num > cliquelist->size[value] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &cliquelist->cliques[value], cliquelist->size[value], newsize) ); /*lint !e866*/
      cliquelist->size[value] = newsize;
   }
   assert(num <= cliquelist->size[value]);

   return SCIP_OKAY;
}

/** adds a clique to the clique list */
SCIP_RETCODE SCIPcliquelistAdd(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             value,              /**< value of the variable for which the clique list should be extended */
   SCIP_CLIQUE*          clique              /**< clique that should be added to the clique list */
   )
{
   unsigned int id;
   int i = 0;

   assert(cliquelist != NULL);

   /* insert clique into list, sorted by clique id */
   id = clique->id; /*lint !e732*/

   /* allocate memory */
   if( *cliquelist == NULL )
   {
      SCIP_CALL( cliquelistCreate(cliquelist, blkmem) );
   }
   else
   {
      if( (*cliquelist)->cliques[value] != NULL )
      {
         for( i = (*cliquelist)->ncliques[value]; i > 0 && (*cliquelist)->cliques[value][i - 1]->id > id; --i ); /*lint !e574 !e722*/
         /* do not put the same clique twice in the cliquelist */
         if( i > 0 && (*cliquelist)->cliques[value][i - 1]->id == id )
            return SCIP_OKAY;
      }
   }
   SCIP_CALL( cliquelistEnsureSize(*cliquelist, blkmem, set, value, (*cliquelist)->ncliques[value] + 1) );

   SCIPsetDebugMsg(set, "adding clique %u to cliquelist %p value %u (length: %d)\n",
      clique->id, (void*)*cliquelist, value, (*cliquelist)->ncliques[value]);

   BMSmoveMemoryArray(&((*cliquelist)->cliques[value][i+1]), &((*cliquelist)->cliques[value][i]), (*cliquelist)->ncliques[value] - i); /*lint !e866*/

   (*cliquelist)->cliques[value][i] = clique;
   (*cliquelist)->ncliques[value]++;

   return SCIP_OKAY;
}

/** removes a clique from the clique list */
SCIP_RETCODE SCIPcliquelistDel(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             value,              /**< value of the variable for which the clique list should be reduced */
   SCIP_CLIQUE*          clique              /**< clique that should be deleted from the clique list */
   )
{
   int pos;

   assert(cliquelist != NULL);

   /* if a variable appeared twice in its last clique and is now removed both times, the cliquelist is already cleaned
      up after the first removal call of this "doubleton" */
   if( *cliquelist == NULL )
      return SCIP_OKAY;

   SCIPdebugMessage("deleting clique %u from cliquelist %p value %u (length: %d)\n", 
      clique->id, (void*)*cliquelist, value, (*cliquelist)->ncliques[value]);

   pos = cliquesSearchClique((*cliquelist)->cliques[value], (*cliquelist)->ncliques[value], clique);

   /* clique does not exist in cliquelist, the clique should contain multiple entries of the same variable */
   if( pos < 0 )
   {
#ifdef SCIP_MORE_DEBUG
      SCIP_VAR** clqvars;
      SCIP_Bool* clqvalues;
      int nclqvars = SCIPcliqueGetNVars(clique);
      int v;

      assert(nclqvars >= 2);
      assert(SCIPcliqueGetVars(clique) != NULL);
      assert(SCIPcliqueGetValues(clique) != NULL);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &clqvars, SCIPcliqueGetVars(clique), nclqvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &clqvalues, SCIPcliqueGetValues(clique), nclqvars) );

      /* sort variables and corresponding clique values after variables indices */
      SCIPsortPtrBool((void**) clqvars, clqvalues, SCIPvarComp, nclqvars);

      for( v = nclqvars - 1; v > 0; --v )
      {
         if( clqvars[v] == clqvars[v - 1] )
         {
            if( clqvalues[v] == clqvalues[v - 1] || (v > 1 && clqvars[v] == clqvars[v - 2]) )
               break;
         }
      }
      assert(v > 0);

      BMSfreeBlockMemoryArray(blkmem, &clqvalues, nclqvars);
      BMSfreeBlockMemoryArray(blkmem, &clqvars, nclqvars);
#endif
      return SCIP_OKAY;
   }

   assert(0 <= pos && pos < (*cliquelist)->ncliques[value]);
   assert((*cliquelist)->cliques[value][pos] == clique);

   /* remove clique from list */
   /* @todo maybe buffered deletion */
   (*cliquelist)->ncliques[value]--;
   if( pos < (*cliquelist)->ncliques[value] )
   {
      BMSmoveMemoryArray(&((*cliquelist)->cliques[value][pos]), &((*cliquelist)->cliques[value][pos+1]),
         (*cliquelist)->ncliques[value] - pos); /*lint !e866*/
   }

   /* free cliquelist if it is empty */
   if( (*cliquelist)->ncliques[0] == 0 && (*cliquelist)->ncliques[1] == 0 )
      SCIPcliquelistFree(cliquelist, blkmem);

   return SCIP_OKAY;
}

/** returns whether the given clique lists have a non-empty intersection, i.e. whether there is a clique that appears
 *  in both lists
 */
SCIP_Bool SCIPcliquelistsHaveCommonClique(
   SCIP_CLIQUELIST*      cliquelist1,        /**< first clique list data structure */
   SCIP_Bool             value1,             /**< value of first variable */
   SCIP_CLIQUELIST*      cliquelist2,        /**< second clique list data structure */
   SCIP_Bool             value2              /**< value of second variable */
   )
{
   SCIP_CLIQUE** cliques1;
   SCIP_CLIQUE** cliques2;
   int ncliques1;
   int ncliques2;
   int i1;
   int i2;

   if( cliquelist1 == NULL || cliquelist2 == NULL )
      return FALSE;

   ncliques1 = cliquelist1->ncliques[value1];
   cliques1 = cliquelist1->cliques[value1];
   ncliques2 = cliquelist2->ncliques[value2];
   cliques2 = cliquelist2->cliques[value2];

   i1 = 0;
   i2 = 0;

   if( i1 < ncliques1 && i2 < ncliques2 )
   {
      unsigned int cliqueid;

      /* make the bigger clique the first one */
      if( ncliques2 > ncliques1 )
      {
         SCIP_CLIQUE** tmpc;
         int tmpi;

         tmpc = cliques1;
         tmpi = ncliques1;
         cliques1 = cliques2;
         ncliques1 = ncliques2;
         cliques2 = tmpc;
         ncliques2 = tmpi;
      }

      /* check whether both clique lists have a same clique */
      while( TRUE )  /*lint !e716*/
      {
         cliqueid = SCIPcliqueGetId(cliques2[i2]);

         /* if last item in clique1 has a smaller index than the actual clique in clique2, than cause of increasing order
          * there will be no same item and we can stop */
         if( SCIPcliqueGetId(cliques1[ncliques1 - 1]) < cliqueid )
            break;

         while( SCIPcliqueGetId(cliques1[i1]) < cliqueid )
         {
            ++i1;
            assert(i1 < ncliques1);
         }
         cliqueid = SCIPcliqueGetId(cliques1[i1]);

         /* if last item in clique2 has a smaller index than the actual clique in clique1, than cause of increasing order
          * there will be no same item and we can stop */
         if( SCIPcliqueGetId(cliques2[ncliques2 - 1]) < cliqueid )
            break;

         while( SCIPcliqueGetId(cliques2[i2]) < cliqueid )
         {
            ++i2;
            assert(i2 < ncliques2);
         }
         if( SCIPcliqueGetId(cliques2[i2]) == cliqueid )
            return TRUE;
      }
   }
   return FALSE;
}

/** removes all listed entries from the cliques */
void SCIPcliquelistRemoveFromCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< active problem variable the clique list belongs to */
   SCIP_Bool             irrelevantvar       /**< has the variable become irrelevant, meaning that equality
                                              *   cliques need to be relaxed? */
   )
{
   assert(SCIPvarIsBinary(var));

   if( cliquelist != NULL )
   {
      int value;

      SCIPdebugMessage("removing variable <%s> from cliques (%d with value 0, %d with value 1)\n",
         SCIPvarGetName(var), cliquelist->ncliques[0], cliquelist->ncliques[1]);

      for( value = 0; value < 2; ++value )
      {
         int i;

         assert(SCIPvarGetCliques(var, (SCIP_Bool)value) == cliquelist->cliques[value]);
         assert(SCIPvarGetNCliques(var, (SCIP_Bool)value) == cliquelist->ncliques[value]);

         /* it is important to iterate from the end of the array because otherwise, each removal causes
          * a memory move of the entire array
          */
         for( i = cliquelist->ncliques[value] - 1; i >= 0; --i )
         {
            SCIP_CLIQUE* clique;

            clique = cliquelist->cliques[value][i];
            assert(clique != NULL);

            SCIPdebugMessage(" -> removing variable <%s> == %d from clique %u (size %d)\n",
               SCIPvarGetName(var), value, clique->id, clique->nvars);

            SCIPcliqueDelVar(clique, cliquetable, var, (SCIP_Bool)value);

            if( irrelevantvar )
               clique->equation = FALSE;

#ifdef SCIP_MORE_DEBUG
            /* during the cleanup step, we skip the consistency check because clique may be temporarily inconsistent */
            if( ! cliquetable->incleanup || clique->index > 0 )
            {
               cliqueCheck(clique);
            }
#endif
         }
      }
   }
}

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashgetkeyClique)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(hashkeyeqClique)
{  /*lint --e{715}*/
   SCIP_CLIQUE* clique1;
   SCIP_CLIQUE* clique2;
   int i;

   clique1 = (SCIP_CLIQUE*)key1;
   clique2 = (SCIP_CLIQUE*)key2;
   assert(clique1 != NULL);
   assert(clique2 != NULL);

   if( clique1->nvars != clique2->nvars )
      return FALSE;

   /* the variables are sorted: we can simply check the equality of each pair of variable/values */
   for( i = 0; i < clique1->nvars; ++i )
   {
      if( clique1->vars[i] != clique2->vars[i] || clique1->values[i] != clique2->values[i] )
         return FALSE;
   }

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashkeyvalClique)
{  /*lint --e{715}*/
   SCIP_CLIQUE* clique;

   clique = (SCIP_CLIQUE*)key;

   return clique->nvars == 0 ? 0 : SCIPhashTwo(SCIPcombineTwoInt(SCIPvarGetIndex(clique->vars[0]), \
                                                                 SCIPvarGetIndex(clique->vars[clique->nvars-1])), \
                                               SCIPcombineThreeInt(clique->nvars, \
                                                                   clique->values[0], \
                                                                   clique->values[clique->nvars-1]));
}

#define HASHTABLE_CLIQUETABLE_SIZE 100

/** creates a clique table data structure */
SCIP_RETCODE SCIPcliquetableCreate(
   SCIP_CLIQUETABLE**    cliquetable,        /**< pointer to store clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int hashtablesize;

   assert(cliquetable != NULL);

   SCIP_ALLOC( BMSallocMemory(cliquetable) );

   /* create hash table to test for multiple cliques */
   hashtablesize = HASHTABLE_CLIQUETABLE_SIZE;
   hashtablesize = MAX(hashtablesize, (set->misc_usesmalltables ? SCIP_HASHSIZE_CLIQUES_SMALL : SCIP_HASHSIZE_CLIQUES));
   SCIP_CALL( SCIPhashtableCreate(&((*cliquetable)->hashtable), blkmem, hashtablesize,
         hashgetkeyClique, hashkeyeqClique, hashkeyvalClique, NULL) );

   (*cliquetable)->varidxtable = NULL;
   (*cliquetable)->djset = NULL;
   (*cliquetable)->cliques = NULL;
   (*cliquetable)->ncliques = 0;
   (*cliquetable)->size = 0;
   (*cliquetable)->ncreatedcliques = 0;
   (*cliquetable)->ncleanupfixedvars = 0;
   (*cliquetable)->ncleanupaggrvars = 0;
   (*cliquetable)->ndirtycliques = 0;
   (*cliquetable)->nentries = 0;
   (*cliquetable)->incleanup = FALSE;
   (*cliquetable)->compsfromscratch = FALSE;
   (*cliquetable)->ncliquecomponents = -1;

   return SCIP_OKAY;
}

/** frees a clique table data structure */
SCIP_RETCODE SCIPcliquetableFree(
   SCIP_CLIQUETABLE**    cliquetable,        /**< pointer to store clique table data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;

   assert(cliquetable != NULL);
   assert(*cliquetable != NULL);

   /* free all cliques */
   for( i = (*cliquetable)->ncliques - 1; i >= 0; --i )
   {
      cliqueFree(&(*cliquetable)->cliques[i], blkmem);
   }

   /* free disjoint set (union find) data structure */
   if( (*cliquetable)->djset != NULL )
      SCIPdisjointsetFree(&(*cliquetable)->djset, blkmem);

   /* free hash table for variable indices */
   if( (*cliquetable)->varidxtable != NULL )
      SCIPhashmapFree(&(*cliquetable)->varidxtable);

   /* free clique table data */
   BMSfreeMemoryArrayNull(&(*cliquetable)->cliques);

   /* free hash table */
   SCIPhashtableFree(&((*cliquetable)->hashtable));

   BMSfreeMemory(cliquetable);

   return SCIP_OKAY;
}

/** ensures, that clique table arrays can store at least num entries */
static
SCIP_RETCODE cliquetableEnsureSize(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(cliquetable != NULL);

   if( num > cliquetable->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&cliquetable->cliques, newsize) );
      cliquetable->size = newsize;
   }
   assert(num <= cliquetable->size);

   return SCIP_OKAY;
}

/** sort variables regarding their index and remove multiple entries of the same variable */
static
SCIP_RETCODE sortAndMergeClique(
   SCIP_VAR**            clqvars,            /**< variables of a clique */
   SCIP_Bool*            clqvalues,          /**< clique values, active or negated, for the variables in a clique */
   int*                  nclqvars,           /**< number of clique variables */
   SCIP_Bool*            isequation,         /**< do we have an equation clique at hand? */
   SCIP_CLIQUE*          clique,             /**< clique data structure or NULL */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   int*                  nbdchgs,            /**< pointer to store number of fixed variables */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   SCIP_VAR* var;
   int noldbdchgs;
   int startidx;

   assert(nclqvars != NULL);

   SCIPsetDebugMsg(set, "starting merging %d variables in clique %d\n", *nclqvars, (clique == NULL) ? -1 : (int) clique->id);

   if( *nclqvars <= 1 )
      return SCIP_OKAY;

   assert(clqvars != NULL);
   assert(clqvalues != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(eventqueue != NULL);
   assert(nbdchgs != NULL);
   assert(infeasible != NULL);

   /* sort variables and corresponding clique values regarding variables indices before merging multiples */
   SCIPsortPtrBool((void**) clqvars, clqvalues, SCIPvarComp, *nclqvars);

   var = NULL;
   noldbdchgs = *nbdchgs;
   /* check for multiple occurences or pairs of negations in the variable array, this should be very rare when creating a
    * new clique */
   startidx = *nclqvars - 1;
   while( startidx >= 0 )
   {
      int nones;
      int nzeros;
      int curr;

      var = clqvars[startidx];
       /* only column and loose variables can exist now */
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
      assert(SCIPvarIsBinary(var));

      /* the counters denote the occurrences of the variable var with value TRUE and FALSE as clqvalues, respectively */
      if( clqvalues[startidx] )
      {
         nones = 1;
         nzeros = 0;
      }
      else
      {
         nzeros = 1;
         nones = 0;
      }
      curr = startidx - 1;

      /* loop and increase the counter based on the clique value */
      while( curr >= 0 )
      {
         if( clqvars[curr] == var )
         {
            if( clqvalues[curr] )
               nones++;
            else
               nzeros++;
         }
         else
            /* var is different from variable at index curr */
            break;

         --curr;
      }
      assert(nzeros + nones <= *nclqvars);

      /* single occurrence of the variable */
      if( nones + nzeros == 1 )
      {
         startidx = curr;
         continue;
      }
      /* at least two occurrences of both the variable and its negation; infeasible */
      if( nones >= 2 && nzeros >= 2 )
      {
         *infeasible = TRUE;
         break;
      }

      /* do we have the same variable with the same clique value twice? */
      if( nones >= 2 || nzeros >= 2 )
      {
         SCIP_Bool fixvalue;

         /* other cases should be handled as infeasible */
         assert(nones <= 1 || nzeros <= 1);

         /* ensure that the variable is removed from both clique lists before fixing it */
         if( clique != NULL )
         {
            if( nones >= 1 )
            {
               SCIP_CALL( SCIPvarDelCliqueFromList(var, blkmem, TRUE, clique) );
            }
            if( nzeros >= 1 )
            {
               SCIP_CALL( SCIPvarDelCliqueFromList(var, blkmem, FALSE, clique) );
            }
         }

         /* we fix the variable into the direction of fewer occurrences */
         fixvalue = nones >= 2 ? FALSE : TRUE;

         /* a variable multiple times in one clique forces this variable to be zero */
         SCIP_CALL( SCIPvarFixBinary(var, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               cliquetable, fixvalue, infeasible, nbdchgs) );

         SCIPsetDebugMsg(set, "same var %s twice in a clique with value %d fixed to %d (was %s)\n", SCIPvarGetName(var), fixvalue ? 0 : 1,
               fixvalue ? 1 : 0, *infeasible ? "infeasible" : "feasible");

         if( *infeasible )
            break;
      }

      /* do we have a pair of negated variables? */
      if( nones >= 1 && nzeros >= 1 )
      {
         SCIPsetDebugMsg(set, "var %s is paired with its negation in one clique -> fix all other variables\n", SCIPvarGetName(var));

         /* a pair of negated variables in one clique forces all other variables in the clique to be zero */
         if( nzeros + nones < *nclqvars )
         {
            int w;

            w = *nclqvars - 1;
            while( w >= 0 )
            {
               /* skip all occurrences of variable var itself, these are between curr and startidx */
               if( w == startidx )
               {
                  if( curr == -1 )
                     break;

                  w = curr;
               }
               assert(clqvars[w] != var);

               /* if one of the other variables occurs multiple times, we can
                * 1) deduce infeasibility if it occurs with different values
                * 2) wait for the last occurence to fix it
                */
               while( w > 0 && clqvars[w-1] == clqvars[w] )
               {
                  if( clqvalues[w-1] != clqvalues[w] )
                  {
                     *infeasible = TRUE;
                     break;
                  }
                  --w;
               }

               if( *infeasible )
                  break;

               if( clique != NULL )
               {
                  SCIP_CALL( SCIPvarDelCliqueFromList(clqvars[w], blkmem, clqvalues[w], clique) );
               }
               SCIP_CALL( SCIPvarFixBinary(clqvars[w], blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     branchcand, eventqueue, cliquetable, !clqvalues[w], infeasible, nbdchgs) );

               SCIPsetDebugMsg(set, "fixed var %s with value %u to %d (was %s)\n", SCIPvarGetName(clqvars[w]), clqvalues[w],
                  clqvalues[w] ? 0 : 1, *infeasible ? "infeasible" : "feasible");

               if( *infeasible )
                  break;

               --w;
            }
         }
         /* all variables except for var were fixed */
         if( clique != NULL )
         {
            SCIP_CALL( SCIPvarDelCliqueFromList(var, blkmem, TRUE, clique) );
            SCIP_CALL( SCIPvarDelCliqueFromList(var, blkmem, FALSE, clique) );
         }
         *nclqvars = 0;
         *isequation = FALSE;

         /* break main loop */
         break;
      }

      /* otherwise, we would have an endless loop */
      assert(curr < startidx);
      startidx = curr;
   }


   /* we might have found an infeasibility or reduced the clique to size 0 */
   if( *infeasible || *nclqvars == 0 )
      return SCIP_OKAY;

   /* we fixed a variable because it appears at least twice, now we need to remove the fixings
    *
    * @note that if we are in probing or solving stage, the fixation on the variable might not be carried out yet,
    *       because it was contradicting a local bound
    */
   if( *nbdchgs > noldbdchgs )
   {
      SCIP_VAR* onefixedvar;
      SCIP_Bool onefixedvalue;
      int w;

      /* we initialize the former value only to please gcc */
      onefixedvalue = TRUE;
      onefixedvar = NULL;

      /* remove all inactive variables, thereby checking for a variable that is fixed to one */
      startidx = 0;
      w = 0;
      while( startidx < *nclqvars )
      {
         var = clqvars[startidx];

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
            || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

         /* check if we have a variable fixed to one for this clique */
         if( (clqvalues[startidx] && SCIPvarGetLbGlobal(var) > 0.5) || (!clqvalues[startidx] && SCIPvarGetUbGlobal(var) < 0.5) )
         {
            if( onefixedvar != NULL )
            {
               *infeasible = TRUE;

               SCIPsetDebugMsg(set, "two variables in clique %d fixed to one %s%s and %s%s\n", (clique != NULL) ? (int) clique->id : -1,
                     onefixedvalue ? "" : "~", SCIPvarGetName(onefixedvar), clqvalues[startidx] ? "" : "~", SCIPvarGetName(clqvars[startidx]));
               return SCIP_OKAY;
            }
            onefixedvar = var;
            onefixedvalue = clqvalues[startidx];
         }
         /* only keep active variables */
         else if( SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5 )
         {
            /* we remove all fixed variables */
            if( startidx > w  )
            {
               clqvars[w] = clqvars[startidx];
               clqvalues[w] = clqvalues[startidx];
            }
            ++w;
         }
         else
         {
            /* can we have some variable fixed to one? */
            assert((SCIPvarGetUbGlobal(var) < 0.5 && clqvalues[startidx]) || (SCIPvarGetLbGlobal(var) > 0.5 && !clqvalues[startidx]));

            if( clique != NULL )
            {
               SCIP_CALL( SCIPvarDelCliqueFromList(var, blkmem, clqvalues[startidx], clique) );
            }
         }
         ++startidx;
      }

      /* the clique has been reduced to w active (unfixed) variables */
      *nclqvars = w;

      /* in rare cases, a variable fixed to one is only detected in the merging step */
      if( onefixedvar != NULL )
      {
         SCIPsetDebugMsg(set, "variable %s%s in clique %d fixed to one, fixing all other variables to zero\n",
            onefixedvalue ? "" : "~", SCIPvarGetName(onefixedvar), (clique != NULL) ? (int) clique->id : -1);

         /* handle all active variables by fixing them */
         for( startidx = *nclqvars; startidx >= 0; --startidx )
         {
            assert(SCIPvarGetStatus(clqvars[startidx]) == SCIP_VARSTATUS_COLUMN
                  || SCIPvarGetStatus(clqvars[startidx]) == SCIP_VARSTATUS_LOOSE);

            /* delete the variable from its clique lists and fix it */
            if( onefixedvalue != clqvalues[startidx] || clqvars[startidx] != onefixedvar )
            {
               if( clique != NULL )
               {
                  SCIP_CALL( SCIPvarDelCliqueFromList(clqvars[startidx], blkmem, clqvalues[startidx], clique) );
               }

               /* if one of the other variables occurs multiple times, we can
                * 1) deduce infeasibility if it occurs with different values
                * 2) wait for the last occurence to fix it
                */
               while( startidx > 0 && clqvars[startidx - 1] == clqvars[startidx] )
               {
                  if( clqvalues[startidx - 1] != clqvalues[startidx] )
                  {
                     *infeasible = TRUE;
                     return SCIP_OKAY;
                  }
                  --startidx;
               }

               SCIPsetDebugMsg(set, "fixing variable %s in clique %d to %d\n", SCIPvarGetName(clqvars[startidx]), (clique != NULL) ? (int) clique->id : -1,
                     clqvalues[startidx] ? 0 : 1);

               /* note that the variable status will remain unchanged */
               SCIP_CALL( SCIPvarFixBinary(clqvars[startidx], blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     branchcand, eventqueue, cliquetable, !clqvalues[startidx], infeasible, nbdchgs) );

               if( *infeasible )
                  return SCIP_OKAY;
            }
         }

         /* delete clique from onefixedvars clique list */
         if( clique != NULL ) /*lint !e850*/
         {
            SCIP_CALL( SCIPvarDelCliqueFromList(onefixedvar, blkmem, onefixedvalue, clique) );
         }

         *nclqvars = 0;
         *isequation = FALSE;
      }
   }

   assert(!(*infeasible));

   if( *isequation )
   {
      if( *nclqvars == 0 )
      {
         SCIPsetDebugMsg(set, "empty equation clique left over -> infeasible\n");

         *infeasible = TRUE;
      }
      else if( *nclqvars == 1 )
      {
         assert(SCIPvarGetStatus(clqvars[0]) == SCIP_VARSTATUS_COLUMN
            || SCIPvarGetStatus(clqvars[0]) == SCIP_VARSTATUS_LOOSE);
         /* clearing data and removing variable from its clique list */
         if( clique != NULL )
         {
            SCIP_CALL( SCIPvarDelCliqueFromList(clqvars[0], blkmem, clqvalues[0], clique) );
         }
         SCIP_CALL( SCIPvarFixBinary(clqvars[0], blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
               eventqueue, cliquetable, clqvalues[0], infeasible, nbdchgs) );

         SCIPsetDebugMsg(set, "fixed last clique var %s with value %u to %d (was %s)\n", SCIPvarGetName(clqvars[0]), clqvalues[0],
            clqvalues[0] ? 1 : 0, *infeasible ? "infeasible" : "feasible");

         *nclqvars = 0;
         *isequation = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** helper function that returns the graph node index for a variable during connected component detection */
static
int cliquetableGetNodeIndexBinvar(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             binvar              /**< binary (or binary integer or implicit binary) variable */
   )
{
   SCIP_VAR* activevar;
   int nodeindex;

   assert(binvar != NULL);
   assert(SCIPvarIsBinary(binvar));

   if( cliquetable->varidxtable == NULL )
      return -1;

   /* get active representative of variable */
   if( !SCIPvarIsActive(binvar) )
   {
      activevar = SCIPvarGetProbvar(binvar);
      if( !SCIPvarIsActive(activevar) )
         return -1;
   }
   else
      activevar = binvar;


   assert(SCIPvarIsBinary(activevar));

   /* if the map does not contain an index for this variable, return -1 and mark that the components must be
    * recomputed from scratch
    */
   if( SCIPhashmapExists(cliquetable->varidxtable, (void*)activevar) )
      nodeindex = (int)(size_t)SCIPhashmapGetImage(cliquetable->varidxtable, (void *)activevar);
   else
   {
      nodeindex = -1;
      cliquetable->compsfromscratch = TRUE;
   }

   return nodeindex;
}


/** updates connectedness information for the \p clique */
static
void cliquetableUpdateConnectednessClique(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_CLIQUE*          clique              /**< clique that should be added */
   )
{
   int i;
   int lastnode;
   SCIP_VAR** clqvars;
   int nclqvars;

   assert(cliquetable != NULL);
   assert(clique != NULL);

   /* check if disjoint set (union find) data structure has not been allocated yet */
   if( cliquetable->djset == NULL )
      cliquetable->compsfromscratch = TRUE;

   /* return immediately if a component update is already pending */
   if( cliquetable->compsfromscratch )
      return;

   lastnode = -1;
   clqvars = clique->vars;
   nclqvars = clique->nvars;

   /* loop over variables in the clique and connect the corresponding components */
   for( i = 0; i < nclqvars && !cliquetable->compsfromscratch; ++i )
   {

      /* this method may also detect that the clique table must entirely recompute connectedness */
      int currnode = cliquetableGetNodeIndexBinvar(cliquetable, clqvars[i]);

      /* skip inactive variables */
      if( currnode == - 1 )
         continue;

      /* connect this and the last active variable */
      if( lastnode != -1 )
         SCIPdisjointsetUnion(cliquetable->djset, lastnode, currnode, FALSE);

      lastnode = currnode;
   }
}

/** returns the index of the connected component of the clique graph that the variable belongs to, or -1  */
int SCIPcliquetableGetVarComponentIdx(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   int cmpidx = -1;
   assert(var != NULL);
   assert(cliquetable->djset != NULL);

   /* only binary variables can have a component index
    *(both integer and implicit integer variables can be binary)
    */
   if( SCIPvarIsBinary(var) )
   {
      int nodeindex = cliquetableGetNodeIndexBinvar(cliquetable, var);
      assert(nodeindex < SCIPdisjointsetGetSize(cliquetable->djset));

      if( nodeindex >= 0 )
         cmpidx = SCIPdisjointsetFind(cliquetable->djset, nodeindex);
   }

   return cmpidx;
}


/** adds a clique to the clique table, using the given values for the given variables;
 *  performs implications if the clique contains the same variable twice
 */
SCIP_RETCODE SCIPcliquetableAdd(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR**            vars,               /**< binary variables in the clique: at most one can be set to the given value */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars,              /**< number of variables in the clique */
   SCIP_Bool             isequation,         /**< is the clique an equation or an inequality? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   )
{
   SCIP_VAR** clqvars;
   SCIP_Bool* clqvalues;
   SCIP_CLIQUE* clique;
   SCIP_VAR* var;
   int size;
   int nlocalbdchgs = 0;
   int v;
   int w;

   assert(cliquetable != NULL);
   assert(vars != NULL);

   SCIPsetDebugMsg(set, "trying to add clique %d with %d vars to clique table\n", cliquetable->ncliques, nvars);

   /* check clique on debugging solution */
   SCIP_CALL( SCIPdebugCheckClique(set, vars, values, nvars) ); /*lint !e506 !e774*/

   *infeasible = FALSE;
   size = nvars;

   /* copy clique data */
   if( values == NULL )
   {
      SCIP_CALL( SCIPsetAllocBufferArray(set, &clqvalues, size) );

      /* initialize clique values data */
      for( v = nvars - 1; v >= 0; --v )
         clqvalues[v] = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPsetDuplicateBufferArray(set, &clqvalues, values, size) );
   }
   SCIP_CALL( SCIPsetDuplicateBufferArray(set, &clqvars, vars, size) );

   /* get active variables */
   SCIP_CALL( SCIPvarsGetProbvarBinary(&clqvars, &clqvalues, nvars) );

   /* remove all inactive vars */
   for( v = nvars - 1; v >= 0; --v )
   {
      var = clqvars[v];

      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
      assert(SCIPvarIsBinary(var));

      /* if we have a variables already fixed to one in the clique, fix all other to zero */
      if( ! SCIPvarIsMarkedDeleteGlobalStructures(var) &&
            ((clqvalues[v] && SCIPvarGetLbGlobal(var) > 0.5) || (!clqvalues[v] && SCIPvarGetUbGlobal(var) < 0.5)) )
      {
         SCIPsetDebugMsg(set, "in a clique var %s with value %u is fixed to %d -> fix the rest\n", SCIPvarGetName(var), clqvalues[v], clqvalues[v] ? 1 : 0);

         for( w = nvars - 1; w >= 0; --w )
         {
            SCIP_VAR* clqvar = clqvars[w];

            /* the rare case where the same variable appears several times is handled later during sort and merge */
            if( clqvar != var )
            {
               SCIP_Bool clqval = clqvalues[w];

               /* check if variable is fixed already and terminate with infeasible if this fixing contradicts the clique info */
               if( SCIPvarGetLbGlobal(clqvar) > SCIPvarGetUbGlobal(clqvar) - 0.5 )
               {

                  /* check if fixing contradicts clique constraint */
                  if( (clqval && SCIPvarGetLbGlobal(clqvar) > 0.5)
                     || (! clqval && SCIPvarGetUbGlobal(clqvar) < 0.5) )
                  {
                     *infeasible = TRUE;

                     goto FREEMEM;
                  }

                  continue;
               }

/* the following piece of code is disabled because it assumes the clqvars to be sorted which they aren't until the method
 * sortAndMergeClique() is called
 */
#ifdef SCIP_DISABLED_CODE
               /* if one of the other variables occurs multiple times, we can
                * 1) deduce infeasibility if it occurs with different values
                * 2) wait for the last occurence to fix it
                */
               while( w > 0 && clqvars[w - 1] == clqvars[w] )
               {
                  if( clqvalues[w - 1] != clqvalues[w] )
                  {
                     *infeasible = TRUE;
                     return SCIP_OKAY;
                  }
                  --w;
               }
#endif

               /* fix the variable */
               SCIP_CALL( SCIPvarFixBinary(clqvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                     branchcand, eventqueue, cliquetable, !clqval, infeasible, &nlocalbdchgs) );

               SCIPsetDebugMsg(set, "fixed var %s with value %u to %d (was %s)\n", SCIPvarGetName(clqvar), clqval,
                  clqval ? 0 : 1, *infeasible ? "infeasible" : "feasible");

               if( *infeasible )
                  break;
            }
         }

         /* all variables are fixed so we can stop */
         break; /*lint !e850*/
      }

      /* only unfixed column and loose variables may be member of a clique */
      if( SCIPvarGetUbGlobal(var) - SCIPvarGetLbGlobal(var) < 0.5 ||
            (SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN && SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE) ||
            SCIPvarIsMarkedDeleteGlobalStructures(var) )
      {
         --nvars;
         clqvars[v] = clqvars[nvars];
         clqvalues[v] = clqvalues[nvars];
         isequation = isequation && !(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR) && ! SCIPvarIsMarkedDeleteGlobalStructures(var);
      }
   }

   if( nbdchgs != NULL )
      *nbdchgs += nlocalbdchgs;

   /* did we fix all variables or are we infeasible? */
   if( v >= 0 )
      goto FREEMEM;

   assert(!*infeasible);

   /* check if only one or less variables are left */
   if( v < 0 && nvars <= 1)
   {
      if( isequation )
      {
         if( nvars == 1 )
         {
            nlocalbdchgs = 0;

            SCIP_CALL( SCIPvarFixBinary(clqvars[0], blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                  branchcand, eventqueue, cliquetable, clqvalues[0], infeasible, &nlocalbdchgs) );

            SCIPsetDebugMsg(set, "fixed last clique var %s with value %u to %d (was %s)\n", SCIPvarGetName(clqvars[0]), clqvalues[0],
               clqvalues[0] ? 1 : 0, *infeasible ? "infeasible" : "feasible");

            if( nbdchgs != NULL )
               *nbdchgs += nlocalbdchgs;
         }
         else if( nvars == 0 )
         {
            SCIPsetDebugMsg(set, "empty equation clique left over -> infeasible\n");

            *infeasible = TRUE;
         }
      }

      goto FREEMEM;
   }

   nlocalbdchgs = 0;

   /* sort the variables and values and remove multiple entries of the same variable */
   SCIP_CALL( sortAndMergeClique(clqvars, clqvalues, &nvars, &isequation, NULL, blkmem, set, stat, transprob, origprob, tree,
         reopt, lp, branchcand, eventqueue, cliquetable, &nlocalbdchgs, infeasible) );

   if( nbdchgs != NULL )
      *nbdchgs += nlocalbdchgs;

   /* did we stop early due to a pair of negated variables? */
   if( nvars == 0 || *infeasible )
      goto FREEMEM;

   /* if less than two variables are left over, the clique is redundant */
   if( nvars > 1 )
   {
      SCIP_CLIQUE* sameclique;

      /* @todo check if we can aggregate variables if( clique->equation && clique->nvars == 2 ) */

      /* create the clique data structure */
      SCIP_CALL( cliqueCreateWithData(&clique, blkmem, nvars, clqvars, clqvalues, nvars, cliquetable->ncreatedcliques, isequation) );

      sameclique = (SCIP_CLIQUE*)SCIPhashtableRetrieve(cliquetable->hashtable, (void*)clique);

      if( sameclique == NULL )
      {
         SCIPsetDebugMsg(set, "adding clique %u with %d vars to clique table\n", clique->id, nvars);

         cliquetable->ncreatedcliques++;

         /* add clique to clique table */
         SCIP_CALL( cliquetableEnsureSize(cliquetable, set, cliquetable->ncliques+1) );
         cliquetable->cliques[cliquetable->ncliques] = clique;
         clique->index = cliquetable->ncliques;
         cliquetable->ncliques++;
         cliquetable->nentries += nvars;
         cliquetableUpdateConnectednessClique(cliquetable, clique);

         SCIP_CALL( SCIPhashtableInsert(cliquetable->hashtable, (void*)clique) );

         /* add filled clique to the cliquelists of all corresponding variables */
         SCIP_CALL( SCIPvarsAddClique(clqvars, clqvalues, nvars, blkmem, set, clique) );
      }
      else
      {
         SCIPsetDebugMsg(set, "same clique %p (id: %d) already found in cliquetable -> discard new one\n", (void*) sameclique, sameclique->id);

         /* update equation status of clique */
         /* @note if we might change the order of the clique, e.g. put the equations up front, we should rather remove
          *       the sameclique from the hashmap while adding the new clique
          */
         if( !sameclique->equation && clique->equation )
            sameclique->equation = TRUE;

         cliqueFree(&clique, blkmem);

         goto FREEMEM;
      }
   }
   else
   {
      assert(!isequation);
      assert(nvars == 1);

      goto FREEMEM;
   }
   cliqueCheck(clique);

  FREEMEM:
   SCIPsetFreeBufferArray(set, &clqvalues);
   SCIPsetFreeBufferArray(set, &clqvars);

   return SCIP_OKAY;
}

/** clean up given clique by removing fixed variables */
static
SCIP_RETCODE cliqueCleanup(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   int*                  nchgbds,            /**< pointer to store number of fixed variables */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   assert(clique != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(nchgbds != NULL);
   assert(infeasible != NULL);

   /* do we need to clean up fixed variables? */
   if( !SCIPcliqueIsCleanedUp(clique) )
   {
      SCIP_VAR* onefixedvar = NULL;
      SCIP_Bool onefixedvalue;
      SCIP_Bool needsorting = FALSE;
      int v;
      int w;

      w = clique->startcleanup;

      SCIPsetDebugMsg(set, "Starting clean up of clique %u (size %d) from position %d\n", clique->id, clique->nvars, w);

      /* exchange inactive by active variables */
      for( v = w; v < clique->nvars; ++v )
      {
         SCIP_Bool addvartoclique; /* has the variable status changed such that it needs to be replaced by an active representative? */

         addvartoclique = FALSE;
         if( SCIPvarGetStatus(clique->vars[v]) == SCIP_VARSTATUS_AGGREGATED
            || SCIPvarGetStatus(clique->vars[v]) == SCIP_VARSTATUS_NEGATED
            || SCIPvarGetStatus(clique->vars[v]) == SCIP_VARSTATUS_MULTAGGR )
         {
            needsorting = TRUE;

            SCIP_CALL( SCIPvarGetProbvarBinary(&(clique->vars[v]), &(clique->values[v])) );
            if( SCIPvarGetStatus(clique->vars[v]) == SCIP_VARSTATUS_NEGATED )
            {
               clique->vars[v] = SCIPvarGetNegationVar(clique->vars[v]);
               clique->values[v] = !clique->values[v];
            }
            else if( SCIPvarGetStatus(clique->vars[v]) == SCIP_VARSTATUS_MULTAGGR )
            {
               clique->equation = FALSE;
               continue;
            }

            addvartoclique = TRUE;
         }

         assert(SCIPvarGetStatus(clique->vars[v]) == SCIP_VARSTATUS_COLUMN
               || SCIPvarGetStatus(clique->vars[v]) == SCIP_VARSTATUS_LOOSE
               || SCIPvarGetStatus(clique->vars[v]) == SCIP_VARSTATUS_FIXED);

         /* check for variables that are either fixed to zero or marked for deletion from global structures */
         if( (clique->values[v] && SCIPvarGetUbGlobal(clique->vars[v]) < 0.5) ||
               (!clique->values[v] && SCIPvarGetLbGlobal(clique->vars[v]) > 0.5) ||
               SCIPvarIsMarkedDeleteGlobalStructures(clique->vars[v]) )
         {
            if( clique->equation && SCIPvarIsMarkedDeleteGlobalStructures(clique->vars[v]) )
               clique->equation = FALSE;

            /* the variable will be overwritten by subsequent active variables */
            continue;
         }

         /* check for a variable fixed to one in the clique */
         else if( (clique->values[v] && SCIPvarGetLbGlobal(clique->vars[v]) > 0.5)
               || (!clique->values[v] && SCIPvarGetUbGlobal(clique->vars[v]) < 0.5) )
         {
            if( onefixedvar != NULL )
            {
               *infeasible = TRUE;

               SCIPsetDebugMsg(set, "two variables in clique %u fixed to one %s%s and %s%s\n", clique->id,
                     onefixedvalue ? "" : "~", SCIPvarGetName(onefixedvar), clique->values[v] ? "" : "~",
                           SCIPvarGetName(clique->vars[v])); /*lint !e530*/
               return SCIP_OKAY;
            }
            onefixedvar = clique->vars[v];
            onefixedvalue = clique->values[v];
         }
         else
         {
            assert(SCIPvarGetStatus(clique->vars[v]) != SCIP_VARSTATUS_FIXED);
            assert(w <= v);

            if( w < v )
            {
               clique->vars[w] = clique->vars[v];
               clique->values[w] = clique->values[v];
            }

            /* add clique to active variable if it replaced a aggregated or negated var */
            if( addvartoclique )
            {
               SCIP_CALL( SCIPvarAddCliqueToList(clique->vars[w], blkmem, set, clique->values[w], clique) );
            }
            /* increase indexer of last active, i.e. unfixed, variable in clique */
            ++w;
         }

      }
      clique->nvars = w;

      if( onefixedvar != NULL )
      {
         SCIPsetDebugMsg(set, "variable %s%s in clique %u fixed to one, fixing all other variables to zero\n",
            onefixedvalue ? "" : "~", SCIPvarGetName(onefixedvar), clique->id);

         for( v = 0; v < clique->nvars ; ++v )
         {
            SCIP_VAR* clqvar = clique->vars[v];
            SCIP_Bool clqval = clique->values[v];

            assert(SCIPvarGetStatus(clqvar) == SCIP_VARSTATUS_COLUMN
               || SCIPvarGetStatus(clqvar) == SCIP_VARSTATUS_LOOSE);

            if( onefixedvalue != clqval || clqvar != onefixedvar )
            {
               /* the variable could have been fixed already because it appears more than once in the clique */
               if( SCIPvarGetLbGlobal(clqvar) > SCIPvarGetUbGlobal(clqvar) - 0.5 )
               {
                  /* the fixing must have occurred previously inside this loop. It may happen that
                   * the variable occurs together with its negation in that clique, which is usually
                   * handled during the merge step, but leads to a direct infeasibility because it
                   * contradicts any other variable fixed to one in this clique
                   */
                  if( (clqval && SCIPvarGetLbGlobal(clqvar) > 0.5)
                     || (! clqval && SCIPvarGetUbGlobal(clqvar) < 0.5) )
                  {
                     /* impossible because we would have detected this in the previous cleanup loop */
                     assert(clqvar != onefixedvar);
                     *infeasible = TRUE;

                     return SCIP_OKAY;
                  }
                  continue;
               }

               SCIP_CALL( SCIPvarDelCliqueFromList(clqvar, blkmem, clqval, clique) );

/* the following piece of code is wrong at this point because it assumes sorted variables. can be enabled if sorting
 * of variables is performed earlier than it is now
 */
#ifdef SCIP_DISABLED_CODE
               /* if one of the other variables occurs multiple times, we can
                * 1) deduce infeasibility if it occurs with different values
                * 2) wait for the last occurence to fix it
                */
               while( v < clique->nvars - 1 && clique->vars[v + 1] == clqvar )
               {
                  if( clique->values[v + 1] != clique->values[v] )
                  {
                     *infeasible = TRUE;
                     return SCIP_OKAY;
                  }
                  ++v;
               }
#endif
               SCIPsetDebugMsg(set, "fixing variable %s in clique %u to %d\n", SCIPvarGetName(clqvar), clique->id, clqval ? 0 : 1);

               SCIP_CALL( SCIPvarFixBinary(clqvar, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
                     eventqueue, cliquetable, !clqval, infeasible, nchgbds) );

               if( *infeasible )
                  return SCIP_OKAY;
            }
         }

         if( SCIPvarGetStatus(onefixedvar) == SCIP_VARSTATUS_COLUMN /*lint !e850*/
            || SCIPvarGetStatus(onefixedvar) == SCIP_VARSTATUS_LOOSE )
         {
            SCIP_CALL( SCIPvarDelCliqueFromList(onefixedvar, blkmem, onefixedvalue, clique) );
         }

         clique->nvars = 0;
         clique->equation = FALSE;
         clique->startcleanup = -1;

         return SCIP_OKAY;
      }

      if( clique->equation )
      {
         if( clique->nvars == 0 )
         {
            *infeasible = TRUE;
            return SCIP_OKAY;
         }
         else if( clique->nvars == 1 )
         {
            assert(SCIPvarGetStatus(clique->vars[0]) == SCIP_VARSTATUS_COLUMN
               || SCIPvarGetStatus(clique->vars[0]) == SCIP_VARSTATUS_LOOSE);

            /* clearing data and removing variable from its clique list */
            SCIP_CALL( SCIPvarDelCliqueFromList(clique->vars[0], blkmem, clique->values[0], clique) );

            SCIPsetDebugMsg(set, "fixing last variable %s in clique %u to %d\n", SCIPvarGetName(clique->vars[0]), clique->id,
               clique->values[0] ? 1 : 0);

            SCIP_CALL( SCIPvarFixBinary(clique->vars[0], blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                  branchcand, eventqueue, cliquetable, clique->values[0], infeasible, nchgbds) );

            clique->nvars = 0;
            clique->equation = FALSE;
            clique->startcleanup = -1;

            return SCIP_OKAY;
         }
      }

      if( needsorting )
      {
         SCIP_Bool isequation = clique->equation;

         /* remove multiple entries of the same variable */
         SCIP_CALL( sortAndMergeClique(clique->vars, clique->values, &(clique->nvars), &isequation, clique, blkmem, set, stat,
               transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable, nchgbds, infeasible) );

         clique->equation = isequation;
      }

      /* @todo check if we can aggregate variables if( clique->equation && clique->nvars == 2 ) */

      clique->startcleanup = -1;
   }
   assert(SCIPcliqueIsCleanedUp(clique));

   return SCIP_OKAY;
}

#ifdef SCIP_MORE_DEBUG
/** checks whether clique appears in all clique lists of the involved variables */
static
SCIP_Bool checkNEntries(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   SCIP_Longint nentries = 0;
   int c;

   assert(cliquetable != NULL);

   for( c = cliquetable->ncliques - 1; c >= 0; --c )
      nentries += cliquetable->cliques[c]->nvars;

   return (nentries == cliquetable->nentries);
}
#else
#define checkNEntries(cliquetable) TRUE
#endif

/** removes all empty and single variable cliques from the clique table; removes double entries from the clique table
 *
 * @note cliques can be processed several times by this method
 *
 * @todo try to detect infeasible implications, e.g., x = 1 => (y = 0 && y = 1)
 */
SCIP_RETCODE SCIPcliquetableCleanup(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int*                  nchgbds,            /**< pointer to store number of fixed variables */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   assert(cliquetable != NULL);
   assert(stat != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* check if we have anything to do */
   if( stat->npresolfixedvars == cliquetable->ncleanupfixedvars
      && stat->npresolaggrvars == cliquetable->ncleanupaggrvars
      && cliquetable->ndirtycliques == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "cleaning up clique table with %d cliques (with %" SCIP_LONGINT_FORMAT " entries)\n", cliquetable->ncliques, cliquetable->nentries);

   /* delay events */
   SCIP_CALL( SCIPeventqueueDelay(eventqueue) );

   cliquetable->incleanup = TRUE;
   while( cliquetable->ndirtycliques > 0 && !(*infeasible) )
   {
      SCIP_CLIQUE* clique;
      SCIP_CLIQUE* sameclique;

      clique = cliquetable->cliques[0];

      assert(!SCIPcliqueIsCleanedUp(clique));

      /* remove not clean up clique from hastable */
      SCIP_CALL( SCIPhashtableRemove(cliquetable->hashtable, (void*)clique) );
      cliquetable->nentries -= clique->nvars;
      assert(cliquetable->nentries >= 0);

      SCIP_CALL( cliqueCleanup(clique, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
            cliquetable, nchgbds, infeasible) );

      if( *infeasible )
         break;

      assert(SCIPcliqueIsCleanedUp(clique));

      /* swap freshly cleaned clique with last dirty clique */
      cliquetable->ndirtycliques--;
      cliquetableSwapCliques(cliquetable, 0, cliquetable->ndirtycliques);
      cliqueCheck(clique);

      /* @todo check if we can aggregate variables if( clique->equation && clique->nvars == 2 && SCIPsetGetStage(set) == SCIP_STAGE_PRESOLVING */
#if 0
      if( clique->nvars == 2 && clique->equation && SCIPsetGetStage(set) == SCIP_STAGE_PRESOLVING )
      {
         SCIP_VAR* var0;
         SCIP_VAR* var1;
         SCIP_Real scalarx;
         SCIP_Real scalary;
         SCIP_Real rhs = 1.0;
         SCIP_Bool aggregated;

         printf("aggr vars, clique %d\n", clique->id);

         if( SCIPvarGetType(clique->vars[0]) >= SCIPvarGetType(clique->vars[1]) )
         {
            var0 = clique->vars[0];
            var1 = clique->vars[1];

            if( !clique->values[0] )
            {
               scalarx = -1.0;
               rhs -= 1.0;
            }
            else
               scalarx = 1.0;

            if( !clique->values[1] )
            {
               scalary = -1.0;
               rhs -= 1.0;
            }
            else
               scalary = 1.0;
         }
         else
         {
            var0 = clique->vars[1];
            var1 = clique->vars[0];

            if( !clique->values[0] )
            {
               scalary = -1.0;
               rhs -= 1.0;
            }
            else
               scalary = 1.0;

            if( !clique->values[1] )
            {
               scalarx = -1.0;
               rhs -= 1.0;
            }
            else
               scalarx = 1.0;
         }

         assert((SCIPvarGetStatus(var0) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var0) == SCIP_VARSTATUS_COLUMN)
            && (SCIPvarGetStatus(var1) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_COLUMN));

         /* aggregate the variable */
         SCIP_CALL( SCIPvarTryAggregateVars(set, blkmem, stat, transprob, origprob, primal,
	    tree, lp, cliquetable, branchcand, eventfilter, eventqueue,
	    var0, var1, scalarx, scalary, rhs, infeasible, &aggregated) );

         assert(aggregated || *infeasible);
      }
#endif

      sameclique = (SCIP_CLIQUE*)SCIPhashtableRetrieve(cliquetable->hashtable, (void*)clique);

      /* check if the clique is already contained in the clique table, or if it is redundant (too small) */
      if( clique->nvars <= 1 || sameclique != NULL )
      {
         int j;

         /* infeasible or fixing should be performed already on trivial clique */
         assert(clique->nvars > 1 || !clique->equation);

         /* if the clique which is already in the hashtable is an inequality and the current clique is an equation, we
          * update the equation status of the old one
          */
         if( clique->nvars > 1 && clique->equation && !sameclique->equation )
         {
            assert(sameclique->nvars >= 2);

            /* @note if we might change the order of the clique, e.g. put the equations up front, we should rather remove
             *       the sameclique from the hashmap while adding the new clique
             */
            sameclique->equation = TRUE;
         }

         /* delete the clique from the variables' clique lists */
         for( j = 0; j < clique->nvars; ++j )
         {
            SCIP_CALL( SCIPvarDelCliqueFromList(clique->vars[j], blkmem, clique->values[j], clique) );
         }

         /* free clique and remove it from clique table */
         cliqueFree(&clique, blkmem);
         cliquetable->ncliques--;
         /* insert a clean clique from the end of the array */
         if( cliquetable->ndirtycliques < cliquetable->ncliques )
         {
            assert(SCIPcliqueIsCleanedUp(cliquetable->cliques[cliquetable->ncliques]));
            cliquetable->cliques[cliquetable->ndirtycliques] = cliquetable->cliques[cliquetable->ncliques];
            cliquetable->cliques[cliquetable->ndirtycliques]->index = cliquetable->ndirtycliques;
         }
      }
      else
      {
         cliquetable->nentries += clique->nvars;

         SCIP_CALL( SCIPhashtableInsert(cliquetable->hashtable, (void*)clique) );
         if( !clique->eventsissued )
         {
            int j;

            /* issue IMPLADDED event on each variable in the clique */
            for( j = 0; j < clique->nvars; ++j )
            {
               SCIP_EVENT* event;

               SCIP_CALL( SCIPeventCreateImplAdded(&event, blkmem, clique->vars[j]) );
               SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, NULL, &event) );
            }
            clique->eventsissued = TRUE;
         }
      }
   }
   cliquetable->incleanup = FALSE;

   /* remember the number of fixed variables and cliques in order to avoid unnecessary cleanups */
   cliquetable->ncleanupfixedvars = stat->npresolfixedvars;
   cliquetable->ncleanupaggrvars = stat->npresolaggrvars;
   assert(*infeasible || cliquetable->ndirtycliques == 0);

   assert(*infeasible || checkNEntries(cliquetable));

   SCIPsetDebugMsg(set, "cleaned up clique table has %d cliques left (with %" SCIP_LONGINT_FORMAT " entries)\n", cliquetable->ncliques, cliquetable->nentries);

   /* process events */
   SCIP_CALL( SCIPeventqueueProcess(eventqueue, blkmem, set, NULL, lp, branchcand, NULL) );

   return SCIP_OKAY;
}

/** computes connected components of the clique table
 *
 *  an update becomes necessary if a clique gets added with variables from different components
 */
SCIP_RETCODE SCIPcliquetableComputeCliqueComponents(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< array of problem variables, sorted by variable type */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars,           /**< number of integer variables */
   int                   nimplvars           /**< number of implicit integer variables */
   )
{
   SCIP_DISJOINTSET* djset;
   int nimplbinvars;
   int v;
   int c;
   int nbinvarstotal;
   int ndiscvars;
   int nnonbinvars;

   SCIP_CLIQUE** cliques;

   assert(cliquetable != NULL);
   assert(vars != NULL);

   nimplbinvars = 0;
   cliquetable->compsfromscratch = FALSE;
   ndiscvars = nbinvars + nintvars + nimplvars;

   /* detect integer and implicit integer variables with bounds {0,1} because they might appear in cliques, as well */
   for( v = nbinvars; v < ndiscvars; ++v )
   {
      if( SCIPvarIsBinary(vars[v]) )
         ++nimplbinvars;
   }

   /* count binary and implicit binary variables */
   nbinvarstotal = nbinvars + nimplbinvars;

   /* if there are no binary variables, return */
   if( nbinvarstotal == 0 )
   {
      SCIPsetDebugMsg(set, "0 binary variables in total --> 0 connected components in the clique table");
      cliquetable->ncliquecomponents = 0;
      return SCIP_OKAY;
   }

   /* no cliques are present */
   if( cliquetable->ncliques == 0 )
   {
      SCIPsetDebugMsg(set, "0 cliques --> Clique table has %d isolated nodes", nbinvarstotal);
      cliquetable->ncliquecomponents = nbinvarstotal;
      return SCIP_OKAY;
   }

   /* create or clear the variable index mapping */
   if( cliquetable->varidxtable == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&(cliquetable)->varidxtable, blkmem, ndiscvars) );
   }
   else
   {
      SCIP_CALL( SCIPhashmapRemoveAll(cliquetable->varidxtable) );
   }

   /* loop through variables and store their respective positions in the hash map if they are binary */
   for( v = 0; v < ndiscvars; ++v )
   {
      SCIP_VAR* var = vars[v];
      if( SCIPvarIsBinary(var) )
      {
         /* consider only active representatives */
         if( SCIPvarIsActive(var) )
         {
            SCIP_CALL( SCIPhashmapInsert(cliquetable->varidxtable, (void*)var, (void*)(size_t)v) );
         }
         else
         {
            var = SCIPvarGetProbvar(var);
            if( SCIPvarIsActive(var) )
            {
               SCIP_CALL( SCIPhashmapInsert(cliquetable->varidxtable, (void*)var, (void*)(size_t)v) );
            }
         }
      }
   }

   /* free previous disjoint set (union find) data structure which has become outdated if new variables are present */
   if( cliquetable->djset != NULL )
      SCIPdisjointsetFree(&cliquetable->djset, blkmem);

   /* we need to consider integer and implicit integer variables for which SCIPvarIsBinary() returns TRUE.
    * These may be scattered across the ninteger + nimplvars implicit integer variables.
    * For simplicity, we add all integer and implicit integer variables as nodes to the graph, and subtract
    * the amount of nonbinary integer and implicit integer variables afterwards.
    */
   SCIP_CALL( SCIPdisjointsetCreate(&cliquetable->djset, blkmem, ndiscvars) );
   djset = cliquetable->djset;

   /* subtract all (implicit) integer for which SCIPvarIsBinary() returns FALSE */
   nnonbinvars = (nintvars + nimplvars) - nimplbinvars;

   cliques = cliquetable->cliques;

   /* for every clique, we connect clique variable components, treating a clique as a path
    *
    * if the graph turns out to be completely connected (except for the nonbinary variables), we terminate */
   for( c = 0; c < cliquetable->ncliques && SCIPdisjointsetGetComponentCount(djset) > 1 + nnonbinvars; ++c )
   {
      SCIP_CLIQUE* clique;

      clique = cliques[c];
      cliquetableUpdateConnectednessClique(cliquetable, clique);
   }

   /* subtract superfluous integer and implicit integer variables added to the auxiliary graph */
   cliquetable->ncliquecomponents = SCIPdisjointsetGetComponentCount(djset) - nnonbinvars;
   assert(cliquetable->ncliquecomponents >= 0);
   assert(cliquetable->ncliquecomponents <= nbinvarstotal);

   SCIPsetDebugMsg(set, "connected components detection: %d comps (%d clqs, %d vars)", cliquetable->ncliquecomponents, cliquetable->ncliques, nbinvarstotal);

   return SCIP_OKAY;
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPvboundsGetNVbds
#undef SCIPvboundsGetVars
#undef SCIPvboundsGetCoefs
#undef SCIPvboundsGetConstants
#undef SCIPimplicsGetNImpls
#undef SCIPimplicsGetVars
#undef SCIPimplicsGetTypes
#undef SCIPimplicsGetBounds
#undef SCIPimplicsGetIds
#undef SCIPcliqueGetNVars
#undef SCIPcliqueGetVars
#undef SCIPcliqueGetValues
#undef SCIPcliqueGetId
#undef SCIPcliqueGetIndex
#undef SCIPcliqueIsCleanedUp
#undef SCIPcliqueIsEquation
#undef SCIPcliquelistGetNCliques
#undef SCIPcliquelistGetCliques
#undef SCIPcliquelistCheck
#undef SCIPcliquetableGetNCliques
#undef SCIPcliquetableGetCliques
#undef SCIPcliquetableGetNEntries
#undef SCIPcliquetableGetNCliqueComponents
#undef SCIPcliquetableNeedsComponentUpdate

/** gets number of variable bounds contained in given variable bounds data structure */
int SCIPvboundsGetNVbds(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   return vbounds != NULL ? vbounds->len : 0;
}

/** gets array of variables contained in given variable bounds data structure */
SCIP_VAR** SCIPvboundsGetVars(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   return vbounds != NULL ? vbounds->vars : NULL;
}

/** gets array of coefficients contained in given variable bounds data structure */
SCIP_Real* SCIPvboundsGetCoefs(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   return vbounds != NULL ? vbounds->coefs : NULL;
}

/** gets array of constants contained in given variable bounds data structure */
SCIP_Real* SCIPvboundsGetConstants(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   return vbounds != NULL ? vbounds->constants : NULL;
}

/** gets number of implications for a given binary variable fixing */
int SCIPimplicsGetNImpls(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->nimpls[varfixing] : 0;
}

/** gets array with implied variables for a given binary variable fixing */
SCIP_VAR** SCIPimplicsGetVars(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->vars[varfixing] : NULL;
}

/** gets array with implication types for a given binary variable fixing */
SCIP_BOUNDTYPE* SCIPimplicsGetTypes(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->types[varfixing] : NULL;
}

/** gets array with implication bounds for a given binary variable fixing */
SCIP_Real* SCIPimplicsGetBounds(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->bounds[varfixing] : NULL;
}

/** Gets array with unique implication identifiers for a given binary variable fixing.
 *  If an implication is a shortcut, i.e., it was added as part of the transitive closure of another implication,
 *  its id is negative, otherwise it is nonnegative.
 */
int* SCIPimplicsGetIds(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->ids[varfixing] : NULL;
}

/** gets number of variables in the cliques */
int SCIPcliqueGetNVars(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->nvars;
}

/** gets array of active problem variables in the cliques */
SCIP_VAR** SCIPcliqueGetVars(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->vars;
}

/** gets array of values of active problem variables in the cliques, i.e. whether the variable is fixed to FALSE or
 *  to TRUE in the clique
 */
SCIP_Bool* SCIPcliqueGetValues(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->values;
}

/** gets unique identifier of the clique */
unsigned int SCIPcliqueGetId(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   int id;

   assert(clique != NULL);

   id = clique->id;

   return (unsigned int)id;
}

/** gets index of the clique in the clique table */
int SCIPcliqueGetIndex(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->index;
}

/** gets unique identifier of the clique */
SCIP_Bool SCIPcliqueIsCleanedUp(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return (clique->startcleanup == -1);
}

/** return whether the given clique is an equation */
SCIP_Bool SCIPcliqueIsEquation(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return (SCIP_Bool)(clique->equation); /*lint !e571*/
}

/** returns the number of cliques stored in the clique list */
int SCIPcliquelistGetNCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_Bool             value               /**< value of the variable for which the cliques should be returned */
   )
{
   return cliquelist != NULL ? cliquelist->ncliques[value] : 0;
}

/** returns the cliques stored in the clique list, or NULL if the clique list is empty */
SCIP_CLIQUE** SCIPcliquelistGetCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_Bool             value               /**< value of the variable for which the cliques should be returned */
   )
{
   return cliquelist != NULL ? cliquelist->cliques[value] : NULL;
}

/** checks whether variable is contained in all cliques of the cliquelist */
void SCIPcliquelistCheck(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_VAR*             var                 /**< variable, the clique list belongs to */
   )
{
   /* @todo might need to change ifndef NDEBUG to ifdef SCIP_MORE_DEBUG because it can take at lot of time to check for
    *       correctness
    */
#ifndef NDEBUG
   int value;

   assert(SCIPvarGetNCliques(var, FALSE) == SCIPcliquelistGetNCliques(cliquelist, FALSE));
   assert(SCIPvarGetCliques(var, FALSE) == SCIPcliquelistGetCliques(cliquelist, FALSE));
   assert(SCIPvarGetNCliques(var, TRUE) == SCIPcliquelistGetNCliques(cliquelist, TRUE));
   assert(SCIPvarGetCliques(var, TRUE) == SCIPcliquelistGetCliques(cliquelist, TRUE));

   for( value = 0; value < 2; ++value )
   {
      SCIP_CLIQUE** cliques;
      int ncliques;
      int i;

      ncliques = SCIPcliquelistGetNCliques(cliquelist, (SCIP_Bool)value);
      cliques = SCIPcliquelistGetCliques(cliquelist, (SCIP_Bool)value);
      for( i = 0; i < ncliques; ++i )
      {
         SCIP_CLIQUE* clique;
         int pos;

         clique = cliques[i];
         assert(clique != NULL);

         pos = SCIPcliqueSearchVar(clique, var, (SCIP_Bool)value);
         assert(0 <= pos && pos < clique->nvars);
         assert(clique->vars[pos] == var);
         assert(clique->values[pos] == (SCIP_Bool)value);
      }
   }
#endif
}

/** gets the number of cliques stored in the clique table */
int SCIPcliquetableGetNCliques(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(cliquetable != NULL);

   return cliquetable->ncliques;
}

/** gets the number of cliques created so far by the clique table */
int SCIPcliquetableGetNCliquesCreated(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(cliquetable != NULL);

   return cliquetable->ncreatedcliques;
}

/** gets the array of cliques stored in the clique table */
SCIP_CLIQUE** SCIPcliquetableGetCliques(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(cliquetable != NULL);

   return cliquetable->cliques;
}

/** gets the number of entries in the whole clique table */
SCIP_Longint SCIPcliquetableGetNEntries(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(cliquetable != NULL);

   return cliquetable->nentries;
}

/** returns the number of clique components, or -1 if update is necessary first */
int SCIPcliquetableGetNCliqueComponents(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   return cliquetable->compsfromscratch ? -1 : cliquetable->ncliquecomponents;
}

/** returns TRUE iff the connected clique components need an update (because new cliques were added) */
SCIP_Bool SCIPcliquetableNeedsComponentUpdate(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   return cliquetable->compsfromscratch || cliquetable->djset == NULL;
}
