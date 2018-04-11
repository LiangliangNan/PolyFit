/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SOPLEX_LEGACY

#include <iostream>
#include <assert.h>

#include "spxdefines.h"
#include "ratrecon.h"

namespace soplex
{
#ifdef SOPLEX_WITH_GMP
   static void GMPv_init(mpz_t* vect, int dim)
   {
      for( int i = 0; i < dim; i++ )
         mpz_init(vect[i]);
   }

   static void GMPv_clear(mpz_t* vect, int dim)
   {
      for( int i = 0; i < dim; i++ )
         mpz_clear(vect[i]);
   }

#ifdef SOPLEX_DEBUG
   /** print integer to stream */
   static std::ostream& operator<<(std::ostream& os, const mpz_t* number)
   {
      os << mpz_get_str(0, 10, *number);
      return os;
   }
#endif

   /** this reconstruction routine will set x equal to the mpq vector where each component is the best rational
    *  approximation of xnum / denom with where the GCD of denominators of x is at most Dbound; it will return true on
    *  success and false if more accuracy is required: specifically if componentwise rational reconstruction does not
    *  produce such a vector
    */
   static int Reconstruct(VectorRational& resvec, mpz_t* xnum, mpz_t denom, int dim, const Rational& denomBoundSquared, const DIdxSet* indexSet = 0)
   {
      bool rval = true;
      int done = 0;

      /* denominator must be positive */
      assert(mpz_sgn(denom) > 0);
      assert(denomBoundSquared > 0);

      mpz_t temp;
      mpz_t td;
      mpz_t tn;
      mpz_t Dbound;
      mpz_t gcd;

      mpz_init(gcd); /* stores the gcd of denominators; abort if too big */
      mpz_set_ui(gcd, 1);
      mpz_init(temp);
      mpz_init(td);
      mpz_init(tn);
      mpz_init(Dbound);

#if 1
      mpz_set_q(Dbound, denomBoundSquared.getMpqRef()); /* this is the working bound on the denominator size */
#else
      mpz_set(Dbound, denom); /* this is the working bound on the denominator size */
#endif

      mpz_sqrt(Dbound, Dbound);

      MSG_DEBUG( std::cout << "reconstructing " << dim << " dimensional vector with denominator bound " << mpz_get_str(0, 10, Dbound) << "\n" );

      /* if Dbound is below 2^24 increase it to this value, this avoids changing input vectors that have low denominator
       * because they are floating point representable
       */
      if( mpz_cmp_ui(Dbound,16777216) < 0 )
         mpz_set_ui(Dbound,16777216);

      /* The following represent a_i, the cont frac representation and p_i/q_i, the convergents */
      mpz_t a0;
      mpz_t ai;
      mpz_init(a0);
      mpz_init(ai);

      /* here we use p[2]=pk, p[1]=pk-1,p[0]=pk-2 and same for q */
      mpz_t p[3];
      GMPv_init(p, 3);
      mpz_t q[3];
      GMPv_init(q, 3);

      for( int c = 0; (indexSet == 0 && c < dim) || (indexSet != 0 && c < indexSet->size()); c++ )
      {
         int j = (indexSet == 0 ? c : indexSet->index(c));

         assert(j >= 0);
         assert(j < dim);

         MSG_DEBUG( std::cout << "  --> component " << j << " = " << &xnum[j] << " / denom\n" );

         /* if xnum =0 , then just leave x[j] as zero */
         if( mpz_sgn(xnum[j]) != 0 )
         {
            /* setup n and d for computing a_i the cont. frac. rep */
            mpz_set(tn,xnum[j]);
            mpz_set(td,denom);

            /* divide tn and td by gcd */
            mpz_gcd(temp,tn,td);
            mpz_divexact(tn,tn,temp);
            mpz_divexact(td,td,temp);

            if(mpz_cmp(td,Dbound)<=0)
            {
               MSG_DEBUG( std::cout << "marker 1\n" );

               mpq_set_num(resvec[j].getMpqRef_w(),tn);
               mpq_set_den(resvec[j].getMpqRef_w(),td);
            }
            else
            {
               MSG_DEBUG( std::cout << "marker 2\n" );

               mpz_set_ui(temp,1);

               mpz_fdiv_q(a0,tn,td);
               mpz_fdiv_r(temp,tn,td);
               mpz_set(tn,td);
               mpz_set(td,temp);
               mpz_fdiv_q(ai,tn,td);
               mpz_fdiv_r(temp,tn,td);
               mpz_set(tn,td);
               mpz_set(td,temp);

               mpz_set(p[1],a0);
               mpz_set_ui(p[2],1);
               mpz_addmul(p[2],a0,ai);

               mpz_set_ui(q[1],1);
               mpz_set(q[2],ai);

               done = 0;

               /* if q is already big, skip loop */
               if( mpz_cmp(q[2],Dbound) > 0 )
               {
                  MSG_DEBUG( std::cout << "marker 3\n" );
                  done = 1;
               }

               int cfcnt = 2;
               while( !done && mpz_cmp_ui(td,0) )
               {
                  /* update everything: compute next ai, then update convergents */

                  /* update ai */
                  mpz_fdiv_q(ai, tn, td);
                  mpz_fdiv_r(temp, tn, td);
                  mpz_set(tn, td);
                  mpz_set(td, temp);

                  /* shift p,q */
                  mpz_set(q[0], q[1]);
                  mpz_set(q[1], q[2]);
                  mpz_set(p[0], p[1]);
                  mpz_set(p[1], p[2]);

                  /* compute next p,q */
                  mpz_set(p[2], p[0]);
                  mpz_addmul(p[2], p[1], ai);
                  mpz_set(q[2], q[0]);
                  mpz_addmul(q[2], q[1], ai);

                  if( mpz_cmp(q[2], Dbound) > 0 )
                     done = 1;
                  cfcnt++;

                  MSG_DEBUG( std::cout << "  --> convergent denominator = " << &q[2] << "\n" );
               }

               /* Assign the values */
               mpq_set_num(resvec[j].getMpqRef_w(), p[1]);
               mpq_set_den(resvec[j].getMpqRef_w(), q[1]);
               mpq_canonicalize(resvec[j].getMpqRef_w());
               mpz_gcd(temp, gcd, mpq_denref(resvec[j].getMpqRef()));
               mpz_mul(gcd, gcd, temp);

               if( mpz_cmp(gcd, Dbound) > 0 )
               {
                  MSG_DEBUG( std::cout << "terminating with gcd " << &gcd << " exceeding Dbound " << &Dbound << "\n" );
                  rval = false;
                  goto CLEANUP;
               }
            }
         }
      }

   CLEANUP:
      GMPv_clear(q, 3);
      GMPv_clear(p, 3);
      mpz_clear(td);
      mpz_clear(tn);
      mpz_clear(a0);
      mpz_clear(ai);
      mpz_clear(temp);
      mpz_clear(Dbound);
      mpz_clear(gcd);

      return rval;
   }
#endif



   /** reconstruct a rational vector */
   bool reconstructVector(VectorRational& input, const Rational& denomBoundSquared, const DIdxSet* indexSet)
   {
#ifdef SOPLEX_WITH_GMP
      mpz_t* xnum = 0; /* numerator of input vector */
      mpz_t denom; /* common denominator of input vector */
      int rval = true;
      int dim;

      dim = input.dim();

      /* convert vector to mpz format */
      spx_alloc(xnum, dim);
      GMPv_init(xnum, dim);
      mpz_init_set_ui(denom, 1);

      /* find common denominator */
      if( indexSet == 0 )
      {
         for( int i = 0; i < dim; i++ )
            mpz_lcm(denom, denom, mpq_denref(input[i].getMpqRef()));

         for( int i = 0; i < dim; i++ )
         {
            mpz_mul(xnum[i], denom, mpq_numref(input[i].getMpqRef()));
            mpz_divexact(xnum[i], xnum[i], mpq_denref(input[i].getMpqRef()));
         }
      }
      else
      {
         for( int i = 0; i < indexSet->size(); i++ )
         {
            assert(indexSet->index(i) >= 0);
            assert(indexSet->index(i) < input.dim());
            mpz_lcm(denom, denom, mpq_denref(input[indexSet->index(i)].getMpqRef()));
         }

         for( int i = 0; i < indexSet->size(); i++ )
         {
            int k = indexSet->index(i);
            assert(k >= 0);
            assert(k < input.dim());
            mpz_mul(xnum[k], denom, mpq_numref(input[k].getMpqRef()));
            mpz_divexact(xnum[k], xnum[k], mpq_denref(input[k].getMpqRef()));
         }
      }

      MSG_DEBUG( std::cout << "LCM = " << mpz_get_str(0, 10, denom) << "\n" );

      /* reconstruct */
      rval = Reconstruct(input, xnum, denom, dim, denomBoundSquared, indexSet);

      mpz_clear(denom);
      GMPv_clear(xnum, dim);
      spx_free(xnum);

      return rval;
#else
      return false;
#endif
   }



   /** reconstruct a rational solution */
   /**@todo make this a method of class SoPlex */
   bool reconstructSol(SolRational& solution)
   {
#if 0
      DVectorRational buffer;

      if( solution.hasPrimal() )
      {
         buffer.reDim((solution._primal).dim());
         solution.getPrimal(buffer);
         reconstructVector(buffer);
         solution._primal = buffer;

         buffer.reDim((solution._slacks).dim());
         solution.getSlacks(buffer);
         reconstructVector(buffer);
         solution._slacks = buffer;
      }

      if( solution.hasPrimalray() )
      {
         buffer.reDim((solution._primalray).dim());
         solution.getPrimalray(buffer);
         reconstructVector(buffer);
         solution._primalray = buffer;
      }

      if( solution.hasDual() )
      {
         buffer.reDim((solution._dual).dim());
         solution.getDual(buffer);
         reconstructVector(buffer);
         solution._dual = buffer;

         buffer.reDim((solution._redcost).dim());
         solution.getRedcost(buffer);
         reconstructVector(buffer);
         solution._redcost = buffer;
      }

      if( solution.hasDualfarkas() )
      {
         buffer.reDim((solution._dualfarkas).dim());
         solution.getDualfarkas(buffer);
         reconstructVector(buffer);
         solution._dualfarkas = buffer;
      }
#endif
      return true;
   }
} // namespace soplex

#endif
