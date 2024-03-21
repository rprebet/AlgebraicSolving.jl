/**
Author: Mohab Safey El Din



 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <getopt.h>
#include <omp.h>

#include "flint/flint.h"
#include "flint/ulong_extras.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#include "flint/nmod_poly.h"

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

#define MPN_NORM(a, an)                         \
  do {                                          \
    while ((an) != 0 && (a)[(an) - 1] == 0)     \
      (an)--;                                   \
  } while (0)


void fmpz_poly_random_dense(fmpz_poly_t pol, slong deg,   mp_bitcnt_t bits, flint_rand_t state){

  fmpz_t coeff;
  fmpz_init(coeff);

  //  flint_randseed(state, (ulong) time(NULL), (ulong) clock());

  long i;
  flint_randseed(state, (ulong) clock(), (ulong)(10000*clock()));

  for(i=0;i<=deg;i++){
    fmpz_randbits(coeff, state, bits);
    fmpz_poly_set_coeff_fmpz(pol, i, coeff);
    flint_randseed(state, (ulong) (i), (ulong)(10000*clock()*(i+1)));
  }
  fmpz_clear(coeff);
}

void fmpz_poly_random_dense_old(fmpz_poly_t pol, slong deg,   mp_bitcnt_t bits){

  fmpz_t coeff;
  fmpz_init(coeff);

  flint_rand_t state;
  flint_randinit(state);
  //  flint_randseed(state, (ulong) time(NULL), (ulong) clock());

  long i;

  for(i=0;i<=deg;i++){
    fmpz_randbits(coeff, state, bits);
    fmpz_poly_set_coeff_fmpz(pol, i, coeff);
    //    flint_randseed(state, (ulong) (i), (ulong)(10000*time(NULL)*(i+1)));
  }
  fmpz_clear(coeff);
}

void make_mod_poly(nmod_poly_t a_mod, fmpz_poly_t A, mp_limb_t prime){
  for(unsigned long int i=0; i<A->length; i++){
    nmod_poly_set_coeff_ui(a_mod, i, fmpz_fdiv_ui(A->coeffs+i, prime));
  }
}


void run_res_flint(unsigned long int deg, flint_rand_t state){

  unsigned long int newdeg = deg;

  puts("FLINT resultant functions (nbits = deg / 2)");

  while(newdeg<=8090){
    unsigned long int deg1 = newdeg, deg2 = deg1;
    unsigned int nbits = deg1 / 2;

    fmpz_poly_t A, B;
    fmpz_poly_init2(A, deg1+1);
    fmpz_poly_init2(B, deg2+1);
    fmpz_poly_random_dense(A, deg1, nbits, state);
    fmpz_poly_random_dense(B, deg2, nbits, state);

    fmpz_t res_normal;
    fmpz_init(res_normal);

    double e_normal = 0, t;
    t = omp_get_wtime();
    fmpz_poly_resultant(res_normal, A, B);
    e_normal = omp_get_wtime() - t;

    fprintf(stdout, "Degree: %lu\n", newdeg);
    fprintf(stdout, "Elapsed time (normal): %f\n\n", e_normal);

    if(fmpz_cmp_ui(res_normal, 0)==0){
      puts("Resultant was 0\n");
    }
    newdeg = 2*newdeg;

    fmpz_poly_clear(A);
    fmpz_poly_clear(B);

  }
}

void myresultantold(fmpz_t res, fmpz_poly_t A, fmpz_poly_t B, int nthreads){

  mp_limb_t myprime;

  unsigned long bound ; //=   (deg1+deg2) * nbits + 2 * LOG2(  (deg1 + deg2) ) ;
  bound = (A->length + B->length - 1)*FLINT_BIT_COUNT((10*(A->length + B->length - 1) + 26)/27) + 3;

  mp_bitcnt_t bits1 = FLINT_ABS(_fmpz_vec_max_bits(A->coeffs, A->length)); 
  mp_bitcnt_t bits2 = FLINT_ABS(_fmpz_vec_max_bits(B->coeffs, B->length));

  /* Upper bound Hadamard bound */
  bound += (A->length - 1)*bits1 + (B->length - 1)*bits2;

  unsigned int PRIMESIZE = 32; 

  myprime = n_nextprime(1152921504606847009, 0);
  //  myprime = n_nextprime(4611686018427388039, 0);
  fprintf(stdout,"LOG2(prime) = %d\n", LOG2(myprime));

  unsigned int nbits_prime = 0; //LOG2(myprime);
  long count = 1;
  long count_primes = 0;

  double e_crt = 0, e_modular  = 0, t = omp_get_wtime();
  long size_table = (bound / PRIMESIZE) + 2;
  mp_limb_t * primes_table = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * size_table)) ;
  mp_limb_t * res_mod_table = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * size_table)) ;

  while(nbits_prime < bound){
    myprime = n_nextprime(myprime, 0);

    if(fmpz_divisible_si(fmpz_poly_lead(A), myprime)==0 &&
       fmpz_divisible_si(fmpz_poly_lead(B), myprime)==0){
      nbits_prime += LOG2(myprime);
      primes_table[count_primes] = myprime;
      count_primes++;
    }
    else{
      fprintf(stdout, "! %ld\n", count);
    }
    count++;
  }
  fprintf(stdout, "count = %ld\n", count);
  fflush(stdout);

  nmod_poly_t a, b;
  nmod_poly_init2(a, myprime, A->length);
  nmod_poly_init2(b, myprime, B->length);

  double t_in = 0, e_in = 0;

  omp_set_num_threads(nthreads);
#pragma omp parallel for num_threads(nthreads)
  for(long i = 0; i < count_primes; i++){
    mp_limb_t res_mod;
    mp_limb_t prime = primes_table[i];
    a->mod.n = prime;
    a->mod.ninv = n_preinvert_limb(prime);
    count_leading_zeros(a->mod.norm, prime);

    b->mod.n = myprime;
    b->mod.ninv = n_preinvert_limb(prime);
    count_leading_zeros(b->mod.norm, prime);

    make_mod_poly(a, A, prime);
    make_mod_poly(b, B, prime);
    t_in = omp_get_wtime();
    res_mod = _nmod_poly_resultant(a->coeffs, A->length, b->coeffs, B->length, a->mod);
    e_in += omp_get_wtime() - t_in;
    res_mod_table[i] = res_mod;
  }
  e_modular = omp_get_wtime() - t;

  //utilisation de fmpz_CRT pour recuperer le resultant
  fmpz_t prod;
  fmpz_init(prod);
  fmpz_zero(res);
  fmpz_one(prod);

  puts("CRT starts");
  t = omp_get_wtime() ;
  for(long i = 0; i < count_primes; i++){
    fmpz_CRT_ui(res, res, prod, res_mod_table[i], primes_table[i], 1);
    fmpz_mul_ui(prod, prod, primes_table[i]);
  }
  e_crt = omp_get_wtime() - t;


  fprintf(stdout, "Elapsed time (multi-mod resultant): %f (%f + %f)\n", e_modular + e_crt, e_modular, e_crt);
  fprintf(stdout, "Elapsed time MM in %f\n", e_in);
  fprintf(stdout, "Elapsed time per mod comp. %f\n", e_in / count);

  if (fmpz_cmp_ui(res, 0)==0){
    puts("\nResultant was 0\n");
  }
  fmpz_clear(prod);

  free(primes_table);
  nmod_poly_clear(a);
  nmod_poly_clear(b);

}



//ce resultant est parallelisable mais consomme trop de memoire (et ca nuit a ses performances)
void myresultant_mem_not_eff(fmpz_t res, fmpz_poly_t A, fmpz_poly_t B, int nthreads){

  mp_limb_t myprime;

  unsigned long bound ; //=   (deg1+deg2) * nbits + 2 * LOG2(  (deg1 + deg2) ) ;
  bound = (A->length + B->length - 1)*FLINT_BIT_COUNT((10*(A->length + B->length - 1) + 26)/27) + 3;

  mp_bitcnt_t bits1 = FLINT_ABS(_fmpz_vec_max_bits(A->coeffs, A->length)); 
  mp_bitcnt_t bits2 = FLINT_ABS(_fmpz_vec_max_bits(B->coeffs, B->length));

  /* Upper bound Hadamard bound */
  bound += (A->length - 1)*bits1 + (B->length - 1)*bits2;

  unsigned int PRIMESIZE = 32; 

  myprime = n_nextprime(1152921504606847009, 0);
  //  myprime = n_nextprime(4611686018427388039, 0);
  fprintf(stdout,"LOG2(prime) = %d\n", LOG2(myprime));

  unsigned int nbits_prime = 0; //LOG2(myprime);
  long count = 1;
  long count_primes = 0;

  double e_crt = 0, e_modular  = 0, t = omp_get_wtime();
  long size_table = (bound / PRIMESIZE) + 2;
  mp_limb_t * primes_table = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * size_table)) ;
  mp_limb_t * res_mod_table = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * size_table)) ;

  fmpz * leadA = fmpz_poly_lead(A);
  fmpz * leadB = fmpz_poly_lead(B);
  while(nbits_prime < bound){
    myprime = n_nextprime(myprime, 0);

    if(fmpz_divisible_si(leadA, myprime)==0 &&
       fmpz_divisible_si(leadB, myprime)==0){
      nbits_prime += LOG2(myprime);
      primes_table[count_primes] = myprime;
      count_primes++;
    }
    else{
      fprintf(stdout, "! %ld\n", count);
    }
    count++;
  }
  fprintf(stdout, "count = %ld\n", count);
  fflush(stdout);

  nmod_poly_t *amod_table = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * count_primes));
  nmod_poly_t *bmod_table = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * count_primes));

  for(long i = 0 ; i < count_primes; i++){
    mp_limb_t prime = primes_table[i];
    nmod_poly_init2(amod_table[i], prime, A->length);
    nmod_poly_init2(bmod_table[i], prime, A->length);

    make_mod_poly(amod_table[i], A, prime);
    make_mod_poly(bmod_table[i], B, prime);

  }
  double t_in = 0, e_in = 0;

  omp_set_num_threads(nthreads);
#pragma omp parallel for num_threads(nthreads)
  for(long i =0 ; i < count_primes ; i++){

    t_in = omp_get_wtime();
    mp_limb_t res_mod = _nmod_poly_resultant(amod_table[i]->coeffs, A->length, bmod_table[i]->coeffs, B->length, amod_table[i]->mod);
    e_in += omp_get_wtime() - t_in;
    res_mod_table[i] = res_mod;

  }

  e_modular = omp_get_wtime() - t;

  //utilisation de fmpz_CRT pour recuperer le resultant
  fmpz_t prod;
  fmpz_init(prod);
  fmpz_zero(res);
  fmpz_one(prod);

  t = omp_get_wtime() ;
  for(long i = 0; i < count_primes; i++){
    fmpz_CRT_ui(res, res, prod, res_mod_table[i], primes_table[i], 1);
    fmpz_mul_ui(prod, prod, primes_table[i]);
  }
  e_crt = omp_get_wtime() - t;


  fprintf(stdout, "Elapsed time (multi-mod resultant): %f (%f + %f)\n", e_modular + e_crt, e_modular, e_crt);
  fprintf(stdout, "Elapsed time MM in %f\n", e_in);
  fprintf(stdout, "Elapsed time per mod comp. %f\n", e_in / count);

  if (fmpz_cmp_ui(res, 0)==0){
    puts("\nResultant was 0\n");
  }
  fmpz_clear(prod);

  free(primes_table);
  //nettoyage memoire a terminer

}

unsigned long int compute_bsize_norm(fmpz *l, unsigned long int s){
  fmpz_t sum, sq;
  fmpz_init_set_ui(sum, 0);
  fmpz_init(sq);
  for(unsigned long int i = 0; i < s; i++){
    fmpz_pow_ui(sq, l+i, 2);
    fmpz_add(sum, sum, sq);
  }
  return fmpz_bits(sum) / 2;
}

//l1 et l2 sont 2 tableaux de coefficients de taille respectives l1 et l2
unsigned long int compute_bsize_bound_resultant(fmpz *l1, fmpz *l2, unsigned long int s1, unsigned long int s2){
  unsigned long int b1 = compute_bsize_norm(l1, s1);
  unsigned long int b2 = compute_bsize_norm(l2, s2);
  return (unsigned long int) ( (s1-1) * b2 + (s2-1) * b1 );
}



static mp_limb_t
_nmod_poly_dth_resultant(mp_srcptr poly1, slong len1,
                         mp_srcptr poly2, slong len2,
                         nmod_t mod,
                         mp_ptr head,
                         slong sdeg
                         )
{
  if (poly1 == poly2)
    {
        return 0;
    }
    else if (len2 == 1)
    {
        if (len1 == 1)
        {
            return 1;
        }
        else if (len1 == 2)
        {
            return poly2[0];
        }
        else
        {
            return n_powmod2_ui_preinv(poly2[0], len1 - 1, mod.n, mod.ninv);
        }
    }
    else  /* len1 >= len2 >= 2 */
    {
        mp_limb_t res = 1;

        mp_ptr u, v, r, t, w;
        slong l0, l1, l2;
        mp_limb_t lc;
        w = _nmod_vec_init(3 * len1);

        u = w;
        v = w + len1;
        r = v + len1;
        //u == poly1
        _nmod_vec_set(u, poly1, len1);
        //v == poly2
        _nmod_vec_set(v, poly2, len2);
        l1 = len1;
        l2 = len2;

        do
        {
            l0 = l1;
            l1 = l2;
            //coefficient dominant de v
            lc = v[l1 - 1];

            /* if(DEBUG>3){ */
            /*   fprintf(stderr, "Affichage de v\n"); */
            /*   nmod_vec_print(stderr, v, l1); */
            /* } */

            //u est de longueur l0, v est de longueur l1
            //r = remainder(u, v)
            _nmod_poly_rem(r, u, l0, v, l1, mod);
            l2 = l1 - 1;
            //determine la longueur de r (decremente l2 tant que le top coeff == 0)
            //r est de degre l2 - 1
            MPN_NORM(r, l2);
            //mise a jour
            {
                t = u;
                u = v;
                v = r;
                r = t;
            }

            if (l2 >= 1) 
            {
                lc  = n_powmod2_preinv(lc, l0 - l2, mod.n, mod.ninv);
                res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);

                if (((l0 | l1) & 1) == 0)
                {
                    res = nmod_neg(res, mod);
                }
                if(l2==sdeg+1){
                  mp_limb_t hi, lo;
                  //                  umul_ppmm(hi, lo, v[0], res);
                  //                  NMOD_RED2(tail[0], hi, lo, mod); /* hi already reduced mod n */
                  umul_ppmm(hi, lo, v[sdeg], res);
                  NMOD_RED2(head[0], hi, lo, mod); /* hi already reduced mod n */
                  _nmod_vec_clear(w);
                  return res;
                }
            }
            else 
            {
                if (l1 == 1)
                {
                    lc  = n_powmod2_preinv(lc, l0 - 1, mod.n, mod.ninv);
                    res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);
                }
                else
                {
                    res = 0;
                }
            }
        }
        while (l2 > 0);
        _nmod_vec_clear(w);
        return res;
    }
}



//Calcul du resultant de A et B
//nthreads est le nombre de threads utilises
void my_dth_subresultant(fmpz_t res, fmpz_poly_t A, fmpz_poly_t B,
                         slong sdeg, int nthreads){
  if(A->length<B->length){
    fprintf(stderr, "Warning: polynomials have been given with increasing degree ; possible sign variation\n");
    my_dth_subresultant(res, B, A, sdeg, nthreads);
    return;
  }

  mp_limb_t myprime;

  unsigned long bound ; //=   (deg1+deg2) * nbits + 2 * LOG2(  (deg1 + deg2) ) ;
  bound = (A->length + B->length - 1)*FLINT_BIT_COUNT((10*(A->length + B->length - 1) + 26)/27) + 3;

  mp_bitcnt_t bits1 = FLINT_ABS(_fmpz_vec_max_bits(A->coeffs, A->length)); 
  mp_bitcnt_t bits2 = FLINT_ABS(_fmpz_vec_max_bits(B->coeffs, B->length));

  /* Upper bound Hadamard bound */
  bound += (A->length - 1)*bits1 + (B->length - 1)*bits2;

  //  fprintf(stderr, "Old bound %lu\n", bound);

  bound = compute_bsize_bound_resultant(A->coeffs, B->coeffs, A->length, B->length);

  fprintf(stderr, "Output bit size bound %lu\n", bound);

  //4611686018427388039
  //1152921504606847009

  myprime = n_nextprime(4611686018427388039, 0);
  unsigned int PRIMESIZE = LOG2(myprime); 
  //  myprime = n_nextprime(4611686018427388039, 0);
  fprintf(stderr,"LOG2(prime) = %d\n", LOG2(myprime));

  unsigned int nbits_prime = 0; //LOG2(myprime);
  long count = 1;
  long count_primes = 0;

  double e_end, e_start; 
  long size_table = (bound / PRIMESIZE) + 2;
  mp_limb_t * primes_table_init = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * size_table)) ;
  mp_limb_t * primes_table = primes_table_init;
  mp_limb_t * res_mod_table = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * nthreads)) ;

  fmpz * leadA = fmpz_poly_lead(A);
  fmpz * leadB = fmpz_poly_lead(B);
  e_start = omp_get_wtime();
  while(nbits_prime < bound){
    myprime = n_nextprime(myprime, 0);
    if(fmpz_divisible_si(leadA, myprime)==0 &&
       fmpz_divisible_si(leadB, myprime)==0){
      nbits_prime += LOG2(myprime);
      primes_table_init[count_primes] = myprime;
      count_primes++;
    }
    else{
      fprintf(stderr, "! %ld\n", count);
    }
    count++;
  }
  fprintf(stderr, "count = %ld / count_primes = %ld\n", count, count_primes);
  //  fflush(stdout);

  nmod_poly_t *amod_table = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * nthreads));
  nmod_poly_t *bmod_table = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * nthreads));

  //Il faut les initialiser
  for(int i =0; i < nthreads; i++){
    mp_limb_t prime = primes_table[i];
    nmod_poly_init2(amod_table[i], prime, A->length);
    nmod_poly_init2(bmod_table[i], prime, A->length);
  }
  mp_ptr head = malloc(sizeof(mp_limb_t)*size_table);
  fprintf(stderr, "First initialization done\n\n");
  unsigned long int nblocks = count_primes / nthreads;

  //  double t_in = 0, e_in = 0;
  fmpz_t prod;
  fmpz_init(prod);
  fmpz_zero(res);
  fmpz_one(prod);
  int index = 0;
  omp_set_num_threads(nthreads);
  double e_in = 0;

  for(unsigned long int block = 0; block < nblocks; block++){
    //    unsigned long int index = block*nthreads;
#pragma omp parallel for num_threads(nthreads)
    for(int i =0; i < nthreads; i++){

      mp_limb_t prime = primes_table[i]; 
      //      fprintf(stderr, "prime = %lu\n", prime);
      amod_table[i]->mod.n = prime;
      amod_table[i]->mod.ninv = n_preinvert_limb(prime);
      count_leading_zeros(amod_table[i]->mod.norm, prime);

      bmod_table[i]->mod.n = prime;
      bmod_table[i]->mod.ninv = n_preinvert_limb(prime);
      count_leading_zeros(bmod_table[i]->mod.norm, prime);

      make_mod_poly(amod_table[i], A, prime);
      make_mod_poly(bmod_table[i], B, prime);

      double t_in = omp_get_wtime();
      _nmod_poly_dth_resultant(amod_table[i]->coeffs,
                               A->length,
                               bmod_table[i]->coeffs,
                               B->length,
                               amod_table[i]->mod,
                               head+i+index,
                               sdeg);
      if(e_in==0){
        e_in = omp_get_wtime() - t_in;
        fprintf(stderr, "\nRun time modp comp. = %f\n",e_in);
      }
      res_mod_table[i+index] = head[i+index];

    }

    if(fmpz_fdiv_ui(res, primes_table[0])==res_mod_table[index]){
      //      fprintf(stderr, "Probabilistic stop\n");
      break;
    }
    else{
      for(int i =0; i < nthreads; i++){
        //        unsigned long int index = block*nthreads;
        //        fprintf(stderr, "Ici: %lu (%d)\n", primes_table[i+index], i);
        fmpz_CRT_ui(res, res, prod, res_mod_table[i+index], primes_table[i+index], 1);
        //        fprintf(stderr, "La?\n");
        fmpz_mul_ui(prod, prod, primes_table[i+index]);
        //        fprintf(stderr, "On incremente\n");
        //        fprintf(stderr, "%lu, %lu, %lu\n", primes_table_init[0], primes_table_init[1], primes_table_init[2]);
      }
      primes_table+=nthreads;
    }
  }
  /* //  fprintf(stderr, "after loops\n"); */
  /* //  unsigned long int index = nblocks*nthreads; */
  /* for(unsigned long int i = 0; i < count_primes - index; i++){ */
  /*   mp_limb_t prime = primes_table[0];  */

  /*   amod_table[i]->mod.n = prime; */
  /*   amod_table[i]->mod.ninv = n_preinvert_limb(prime); */
  /*   count_leading_zeros(amod_table[i]->mod.norm, prime); */

  /*   bmod_table[i]->mod.n = prime; */
  /*   bmod_table[i]->mod.ninv = n_preinvert_limb(prime); */
  /*   count_leading_zeros(bmod_table[i]->mod.norm, prime); */

  /*   make_mod_poly(amod_table[i], A, prime); */
  /*   make_mod_poly(bmod_table[i], B, prime); */


  /*   //    t_in = omp_get_wtime(); */
  /*   mp_limb_t res_mod = _nmod_poly_dth_resultant( */
  /*                                            amod_table[i]->coeffs, */
  /*                                            A->length, */
  /*                                            bmod_table[i]->coeffs, */
  /*                                            B->length, */
  /*                                            amod_table[i]->mod, */
  /*                                            head+i+index, */
  /*                                            sdeg); */
  /*   //    e_in += omp_get_wtime() - t_in; */

  /*   res_mod_table[i+index] = head[i+index]; */

  /*   if(fmpz_fdiv_ui(res, primes_table[0])==res_mod_table[i+index]){ */
  /*     break; */
  /*   } */
  /*   else{ */
  /*     fmpz_CRT_ui(res, res, prod, res_mod_table[i+index], primes_table[i+index], 1); */
  /*     fmpz_mul_ui(prod, prod, primes_table[i+index]); */
  /*     primes_table++; */
  /*   } */
  /* } */
  e_end = omp_get_wtime() - e_start;


  //utilisation de fmpz_CRT pour recuperer le resultant
  //Obsolete

  /* t = omp_get_wtime() ; */
  /* for(long i = 0; i < count_primes; i++){ */
  /*   fmpz_CRT_ui(res, res, prod, res_mod_table[i], primes_table[i], 1); */
  /*   fmpz_mul_ui(prod, prod, primes_table[i]); */
  /* } */
  /* e_crt = omp_get_wtime() - t; */
  fprintf(stderr, "Bit size of the result = %ld\n", fmpz_bits(res));
  fprintf(stderr, "Elapsed time (multi-mod resultant): %f\n", e_end);
  //  fprintf(stderr, "Elapsed time MM in %f\n", e_in);
  //  fprintf(stderr, "Elapsed time per mod comp. %f\n", e_in / count);

  if (fmpz_cmp_ui(res, 0)==0){
    fprintf(stderr,"\nResultant was 0\n");
  }
  fmpz_clear(prod);

  free(primes_table_init);

  for(int i = 0; i < nthreads; i++){
    nmod_poly_clear(amod_table[i]);
    nmod_poly_clear(bmod_table[i]);
  }

  free(res_mod_table);
  free(amod_table);
  free(bmod_table);

  //nettoyage memoire a terminer

}


//Calcul du resultant de A et B
//nthreads est le nombre de threads utilises
void myresultant(fmpz_t res, fmpz_poly_t A, fmpz_poly_t B, int nthreads){
  if(A->length<B->length){
    fprintf(stderr, "Warning: polynomials have been given with increasing degree ; possible sign variation\n");
    myresultant(res, B, A, nthreads);
    return;
  }

  mp_limb_t myprime;

  unsigned long bound ; //=   (deg1+deg2) * nbits + 2 * LOG2(  (deg1 + deg2) ) ;
  bound = (A->length + B->length - 1)*FLINT_BIT_COUNT((10*(A->length + B->length - 1) + 26)/27) + 3;

  mp_bitcnt_t bits1 = FLINT_ABS(_fmpz_vec_max_bits(A->coeffs, A->length)); 
  mp_bitcnt_t bits2 = FLINT_ABS(_fmpz_vec_max_bits(B->coeffs, B->length));

  /* Upper bound Hadamard bound */
  bound += (A->length - 1)*bits1 + (B->length - 1)*bits2;

  bound = compute_bsize_bound_resultant(A->coeffs, B->coeffs, A->length, B->length);

  fprintf(stderr, "Output bit size bound %lu\n", bound);

  //4611686018427388039
  myprime = n_nextprime(1152921504606847009, 0);
  unsigned int PRIMESIZE = LOG2(myprime); 
  //  myprime = n_nextprime(4611686018427388039, 0);
  fprintf(stderr,"LOG2(prime) = %d\n", LOG2(myprime));

  unsigned int nbits_prime = 0; //LOG2(myprime);
  long count = 1;
  long count_primes = 0;

  double e_end, e_start; //= omp_get_wtime();
  long size_table = (bound / PRIMESIZE) + 2;
  mp_limb_t * primes_table = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * size_table)) ;
  mp_limb_t * res_mod_table = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * size_table)) ;

  fmpz * leadA = fmpz_poly_lead(A);
  fmpz * leadB = fmpz_poly_lead(B);
  e_start = omp_get_wtime();
  while(nbits_prime < bound){
    myprime = n_nextprime(myprime, 0);
    if(fmpz_divisible_si(leadA, myprime)==0 &&
       fmpz_divisible_si(leadB, myprime)==0){
      nbits_prime += LOG2(myprime);
      primes_table[count_primes] = myprime;
      count_primes++;
    }
    else{
      fprintf(stderr, "! %ld\n", count);
    }
    count++;
  }
  fprintf(stderr, "count = %ld / count_primes = %ld\n", count, count_primes);
  //  fflush(stdout);

  nmod_poly_t *amod_table = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * nthreads));
  nmod_poly_t *bmod_table = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * nthreads));

  //Il faut les initialiser
  for(int i =0; i < nthreads; i++){
    mp_limb_t prime = primes_table[i];
    nmod_poly_init2(amod_table[i], prime, A->length);
    nmod_poly_init2(bmod_table[i], prime, A->length);
  }
  fprintf(stderr, "First initialization done\n\n");
  unsigned long int nblocks = count_primes / nthreads;

  //  double t_in = 0, e_in = 0;
  fmpz_t prod;
  fmpz_init(prod);
  fmpz_zero(res);
  fmpz_one(prod);

  omp_set_num_threads(nthreads);
  for(unsigned long int block = 0; block < nblocks; block++){
    unsigned long int index = block*nthreads;
#pragma omp parallel for num_threads(nthreads)
    for(int i =0; i < nthreads; i++){

      mp_limb_t prime = primes_table[i+index];

      amod_table[i]->mod.n = prime;
      amod_table[i]->mod.ninv = n_preinvert_limb(prime);
      count_leading_zeros(amod_table[i]->mod.norm, prime);

      bmod_table[i]->mod.n = prime;
      bmod_table[i]->mod.ninv = n_preinvert_limb(prime);
      count_leading_zeros(bmod_table[i]->mod.norm, prime);

      make_mod_poly(amod_table[i], A, prime);
      make_mod_poly(bmod_table[i], B, prime);


      //      t_in = omp_get_wtime();
      mp_limb_t res_mod = _nmod_poly_resultant(
                                               amod_table[i]->coeffs,
                                               A->length,
                                               bmod_table[i]->coeffs,
                                               B->length,
                                               amod_table[i]->mod);
      //      e_in += omp_get_wtime() - t_in;

      res_mod_table[i+index] = res_mod;

    }
    if(fmpz_fdiv_ui(res, primes_table[index])==res_mod_table[index]){
      //      fprintf(stderr, "Probabilistic stop\n");
      break;
    }
    else{
      for(int i =0; i < nthreads; i++){
        unsigned long int index = block*nthreads;

        //        mp_limb_t prime = primes_table[i+index];
        fmpz_CRT_ui(res, res, prod, res_mod_table[i+index], primes_table[i+index], 1);
        fmpz_mul_ui(prod, prod, primes_table[i+index]);
      }
    }
  }
  //  fprintf(stderr, "after loops\n");
  unsigned long int index = nblocks*nthreads;
  for(unsigned long int i = 0; i < count_primes - index; i++){
    mp_limb_t prime = primes_table[i+index];

    amod_table[i]->mod.n = prime;
    amod_table[i]->mod.ninv = n_preinvert_limb(prime);
    count_leading_zeros(amod_table[i]->mod.norm, prime);

    bmod_table[i]->mod.n = prime;
    bmod_table[i]->mod.ninv = n_preinvert_limb(prime);
    count_leading_zeros(bmod_table[i]->mod.norm, prime);

    make_mod_poly(amod_table[i], A, prime);
    make_mod_poly(bmod_table[i], B, prime);


    //    t_in = omp_get_wtime();
    mp_limb_t res_mod = _nmod_poly_resultant(
                                             amod_table[i]->coeffs,
                                             A->length,
                                             bmod_table[i]->coeffs,
                                             B->length,
                                             amod_table[i]->mod);
    //    e_in += omp_get_wtime() - t_in;

    res_mod_table[i+index] = res_mod;

    if(fmpz_fdiv_ui(res, primes_table[i+index])==res_mod_table[i+index]){
      break;
    }
    else{
      fmpz_CRT_ui(res, res, prod, res_mod_table[i+index], primes_table[i+index], 1);
      fmpz_mul_ui(prod, prod, primes_table[i+index]);
    }
  }
  e_end = omp_get_wtime() - e_start;


  //utilisation de fmpz_CRT pour recuperer le resultant
  //Obsolete

  /* t = omp_get_wtime() ; */
  /* for(long i = 0; i < count_primes; i++){ */
  /*   fmpz_CRT_ui(res, res, prod, res_mod_table[i], primes_table[i], 1); */
  /*   fmpz_mul_ui(prod, prod, primes_table[i]); */
  /* } */
  /* e_crt = omp_get_wtime() - t; */
  fprintf(stderr, "Bit size of the result = %ld\n", fmpz_bits(res));
  fprintf(stderr, "Elapsed time (multi-mod resultant): %f \n", e_end);
  //  fprintf(stderr, "Elapsed time MM in %f\n", e_in);
  //  fprintf(stderr, "Elapsed time per mod comp. %f\n", e_in / count);

  if (fmpz_cmp_ui(res, 0)==0){
    fprintf(stderr,"\nResultant was 0\n");
  }
  fmpz_clear(prod);

  free(primes_table);
  for(int i = 0; i < nthreads; i++){
    nmod_poly_clear(amod_table[i]);
    nmod_poly_clear(bmod_table[i]);
  }

  free(res_mod_table);
  free(amod_table);
  free(bmod_table);

  //nettoyage memoire a terminer

}

//Compares FLINT resultant with a hand-made multi-mod resultant
//(still using FLINT basics)
void run_compare(unsigned long int deg, flint_rand_t state){

  unsigned long int deg1 = deg, deg2 = deg1;
  unsigned int nbits = deg1 ;

  fmpz_poly_t A, B;
  fmpz_poly_init2(A, deg1+1);
  fmpz_poly_init2(B, deg2+1);
  fmpz_poly_random_dense(A, deg1, nbits, state);
  fmpz_poly_random_dense(B, deg2, nbits, state);

  /* fmpz_poly_random_dense_old(A, deg1, nbits); */
  /* fmpz_poly_random_dense_old(B, deg2, nbits); */

  fmpz_t res_normal, res;
  fmpz_init(res_normal);
  fmpz_init(res);

  double e_normal = 0, t;
  t = omp_get_wtime();
  fmpz_poly_resultant(res_normal, A, B);
  e_normal = omp_get_wtime() - t;
  fprintf(stderr, "Elapsed time (FLINT resultant): %f\n\n", e_normal);

  t = omp_get_wtime();
  myresultant(res, A, B, 1);
  e_normal = omp_get_wtime() - t;
  fprintf(stderr, "Elapsed time (my resultant): %f\n", e_normal);

  if(fmpz_cmp(res, res_normal)!=0){
    puts("\n\t BUG BUG BUG\n");
  }

  puts("");

  /* if(fmpz_cmp_ui(res, 0)==1){ */
  /*   puts("Positive"); */
  /* } */
  /* else{ */
  /*   puts("Negative"); */
  /* } */
  fmpz_clear(res);
  fmpz_clear(res_normal);

}

//runs modular FLINT resultant modulo 32 bits prime
void run_modular16(flint_rand_t state, int ntimes){

  unsigned long int deg = 64;
  mp_limb_t myprime  = 65537;

  puts("\nFLINT resultant modulo 16 bits prime\n");
  while(deg<=2*8192){

    fmpz_poly_t A, B;
    fmpz_poly_init2(A, deg+1);
    fmpz_poly_init2(B, deg+1);
    fmpz_poly_random_dense(A, deg, deg / 2, state);
    fmpz_poly_random_dense(B, deg, deg / 2, state);

    nmod_poly_t a_mod,b_mod;
    nmod_poly_init2(a_mod,myprime, deg+1);
    nmod_poly_init2(b_mod,myprime, deg+1);
    //attention a la gestion des coefficients negatifs
    //  nmod_poly_set_coeff_ui(x,3,9002);
    //  nmod_poly_set_coeff_ui(x,0,6);

    make_mod_poly(a_mod, A, myprime);
    make_mod_poly(b_mod, B, myprime);

    fprintf(stdout,"Degree = %lu\n", deg);

    //  nmod_poly_print(x);flint_printf("\n");
    //  nmod_poly_print(y);flint_printf("\n");
    double e = 0, t = omp_get_wtime();
    int j;
    for(int i = 1; i<=ntimes; i++){
      //      res = nmod_poly_resultant_hgcd(a_mod, b_mod);
      nmod_poly_resultant(a_mod, b_mod); //combine gcd et euclide en fonction du degre
      j = j+i;
    }
    e = omp_get_wtime() - t;
    //    fprintf(stdout,"%lu \n", res);
    fprintf(stdout, "Elapsed time Half gcd resultant = %f\n", e);

    double e2 = 0, t2 = omp_get_wtime();
    for(int i =1 ; i <=ntimes; i++){
      nmod_poly_resultant_euclidean(a_mod, b_mod);
      j = j+i;
    }
    e2 = omp_get_wtime() - t2;
    fprintf(stdout, "Elapsed time Euclidean resultant = %f\n", e2);
    puts("");

    nmod_poly_clear(a_mod);
    nmod_poly_clear(b_mod);

    fmpz_poly_clear(A);
    fmpz_poly_clear(B);

    deg = 2*deg;
  }
}

//runs modular FLINT resultant modulo 32 bits prime
void run_modular32(flint_rand_t state, int ntimes){

  unsigned long int deg = 30;
  mp_limb_t myprime  = 4294967311;

  puts("\nFLINT resultant modulo 32 bits prime\n");
  while(deg<=2*8192){

    fmpz_poly_t A, B;
    fmpz_poly_init2(A, deg+1);
    fmpz_poly_init2(B, deg+1);
    fmpz_poly_random_dense(A, deg, deg / 2, state);
    fmpz_poly_random_dense(B, deg, deg / 2, state);

    nmod_poly_t a_mod,b_mod;
    nmod_poly_init2(a_mod,myprime, deg+1);
    nmod_poly_init2(b_mod,myprime, deg+1);
    //attention a la gestion des coefficients negatifs
    //  nmod_poly_set_coeff_ui(x,3,9002);
    //  nmod_poly_set_coeff_ui(x,0,6);

    make_mod_poly(a_mod, A, myprime);
    make_mod_poly(b_mod, B, myprime);

    fprintf(stdout,"Degree = %lu\n", deg);

    double e = 0, t = omp_get_wtime();
    int j;
    for(int i = 1; i<=ntimes; i++){
      //      res = nmod_poly_resultant_hgcd(a_mod, b_mod);
      nmod_poly_resultant_hgcd(a_mod, b_mod); //combine gcd et euclide en fonction du degre
      j = j+i;
    }
    e = omp_get_wtime() - t;
    fprintf(stdout, "Elapsed time Half gcd resultant = %f\n", e);

    double e2 = 0, t2 = omp_get_wtime();
    for(int i =1 ; i <=ntimes; i++){
      nmod_poly_resultant_euclidean(a_mod, b_mod);
      j = j+i;
    }
    e2 = omp_get_wtime() - t2;
    fprintf(stdout, "Elapsed time Euclidean resultant = %f\n", e2);
    puts("");

    nmod_poly_clear(a_mod);
    nmod_poly_clear(b_mod);

    fmpz_poly_clear(A);
    fmpz_poly_clear(B);

    deg = 2*deg;
  }
}

//runs modular FLINT resultant modulo 60 bits fftprime
void run_modular60fftprime(flint_rand_t state, int ntimes){

  unsigned long int deg = 64;
  mp_limb_t myprime  = 1945555039024054273;

  puts("\nFLINT resultant modulo 60 bits FFT prime\n");
  while(deg<=8192){

    fmpz_poly_t A, B;
    fmpz_poly_init2(A, deg+1);
    fmpz_poly_init2(B, deg+1);
    fmpz_poly_random_dense(A, deg, deg / 2, state);
    fmpz_poly_random_dense(B, deg, deg / 2, state);

    nmod_poly_t a_mod,b_mod;
    nmod_poly_init2(a_mod,myprime, deg+1);
    nmod_poly_init2(b_mod,myprime, deg+1);
    //attention a la gestion des coefficients negatifs
    //  nmod_poly_set_coeff_ui(x,3,9002);
    //  nmod_poly_set_coeff_ui(x,0,6);

    make_mod_poly(a_mod, A, myprime);
    make_mod_poly(b_mod, B, myprime);

    fprintf(stdout,"Degree = %lu\n", deg);

    double e = 0, t = omp_get_wtime();
    int j;
    for(int i = 1; i<=ntimes; i++){
      //      res = nmod_poly_resultant_hgcd(a_mod, b_mod);
      nmod_poly_resultant(a_mod, b_mod); //combine gcd et euclide en fonction du degre
      j = j+i;
    }
    e = omp_get_wtime() - t;

    fprintf(stdout, "Elapsed time Half gcd resultant = %f\n", e);

    double e2 = 0, t2 = omp_get_wtime();
    for(int i =1 ; i <=ntimes; i++){
      nmod_poly_resultant_euclidean(a_mod, b_mod);
      j = j+i;
    }
    e2 = omp_get_wtime() - t2;
    fprintf(stdout, "Elapsed time Euclidean resultant = %f\n", e2);
    puts("");

    nmod_poly_clear(a_mod);
    nmod_poly_clear(b_mod);

    fmpz_poly_clear(A);
    fmpz_poly_clear(B);

    deg = 2*deg;
  }
}

//runs modular FLINT resultant modulo 60 bits prime
void run_modular60prime(flint_rand_t state, int ntimes){

  unsigned long int deg = 256;
  mp_limb_t myprime  = 1152921504606847009;

  puts("\nFLINT resultant modulo 60 bits prime\n");
  while(deg<=8192){

    fmpz_poly_t A, B;
    fmpz_poly_init2(A, deg+1);
    fmpz_poly_init2(B, deg+1);
    fmpz_poly_random_dense(A, deg, deg / 2, state);
    fmpz_poly_random_dense(B, deg, deg / 2, state);

    nmod_poly_t a_mod,b_mod;
    nmod_poly_init2(a_mod,myprime, deg+1);
    nmod_poly_init2(b_mod,myprime, deg+1);
    //attention a la gestion des coefficients negatifs
    //  nmod_poly_set_coeff_ui(x,3,9002);
    //  nmod_poly_set_coeff_ui(x,0,6);

    make_mod_poly(a_mod, A, myprime);
    make_mod_poly(b_mod, B, myprime);

    fprintf(stdout,"Degree = %lu\n", deg);

    double e = 0, t = omp_get_wtime();
    int j;
    for(int i = 1; i<=ntimes; i++){
      //      res = nmod_poly_resultant_hgcd(a_mod, b_mod);
      nmod_poly_resultant_hgcd(a_mod, b_mod); //combine gcd et euclide en fonction du degre
      j = j+i;
    }
    e = omp_get_wtime() - t;

    fprintf(stdout, "Elapsed time Half gcd resultant = %f\n", e/ntimes);

    double e2 = 0, t2 = omp_get_wtime();
    for(int i =1 ; i <=ntimes; i++){
      nmod_poly_resultant_euclidean(a_mod, b_mod);
      j = j+i;
    }
    e2 = omp_get_wtime() - t2;
    fprintf(stdout, "Elapsed time Euclidean resultant = %f\n", e2/ntimes);
    puts("");

    nmod_poly_clear(a_mod);
    nmod_poly_clear(b_mod);

    fmpz_poly_clear(A);
    fmpz_poly_clear(B);

    deg = deg * 2;
  }
}

void display_help(){
  fprintf(stdout, "Basic Flint based C program for computing resultant and subresultant in Z[x]\n");
  fprintf(stdout, "Implemented by M. Safey El Din (Sorbonne Univ./CNRS/Inria), mohab.safey@lip6.fr\n");
  fprintf(stdout, "Provided without any guarantee for ** academic use only ** by students, teachers and researchers hosted in high-schools, universities or academic (non-industtrial) research institutes\n");
  fprintf(stdout, "Any kind of other use requires formal agreement from the author\n\n");

  fprintf(stdout, "-h \t diplays this help\n");
  fprintf(stdout, "-r <d>\t generates random dense polynomials of degree d and runs the resultant computation\n");
  fprintf(stdout, "-f <filename>\t reads polynomials in filename and computes the resultant of those polynomials\n");
  fprintf(stdout, "\n Format file: \nFiles provide 2 polynomials encoded as follows \ndegree\ncoef_0\n...\ncoef_degree\n");

  fprintf(stdout, "\n\nTODO:\n");
  fprintf(stdout, "Main difference with FLINT functions are in memory management+multi-threading\n");
  fprintf(stdout, "Lifting is not that efficient - can be improved\n");
  fprintf(stdout, "Use of fast reductions still needs to be implemented\n");
}



static char *getoptions(int argc, char **argv, char ** input_file,
                        ulong *deg, slong *sdeg, int *nthreads){
  int optch, errflag = 0, fflag = 1;
	extern int opterr;

	char options[] = "hf:r:d:t:";
  char *filename = NULL;
	opterr = 1;

	while ((optch = getopt(argc, argv, options)) != -1)
    switch (optch) {
		case 'f':
      fflag = 0;
      *input_file = optarg;
      optind--;
			break;
    case 'h':
      display_help();
      exit(1);
    case 'r':
      *deg = strtol(optarg, NULL, 10);
      optind--;
      break;
    case 'd':
      *sdeg = strtol(optarg, NULL, 10);
      optind--;
      break;
    case 't':
      *nthreads = strtol(optarg, NULL, 10);
      optind--;
      break;
    default:
      errflag++;
      break;
    }
  if(fflag==1){
    fprintf(stderr, "No given file\n");
    *input_file = NULL;
  }
  if(errflag){
    fprintf(stderr,"Invalid usage\n");
    fprintf(stderr,"%s -f file_name\n", argv[0]);
    exit(1);
  }
  return filename;
}

void get_input_polynomials(FILE *file, ulong *degA, ulong *degB,
                           fmpz_poly_t A, fmpz_poly_t B){
  fmpz_t c;
  fmpz_init(c);
  if(fscanf(file, "%lu", degA)){
    fprintf(stderr, "Degree %lu\n", *degA);
    for(ulong i=0; i <= *degA; i++){
      fmpz_fread(file, c);
      fmpz_poly_set_coeff_fmpz(A, i, c);
    }
  }
  if(fscanf(file, "%lu", degB)){
    fprintf(stderr, "Degree %lu\n", *degB);
    for(ulong i=0; i <= *degB; i++){
      fmpz_fread(file, c);
      fmpz_poly_set_coeff_fmpz(B, i, c);
    }
  }
  fmpz_clear(c);
}

int main(int argc, char *argv[]){
  char **input_file = malloc(sizeof(char**));
  ulong deg = 0;
  slong sdeg = -1;
  int nthreads = 1;
  getoptions(argc, argv, input_file, &deg, &sdeg, &nthreads);

  if(*input_file!=NULL){
    fprintf(stderr, "Reading file...\n");
    FILE *file = fopen(*input_file,"r");
    ulong degA = 0, degB = 0;
    fmpz_poly_t A;
    fmpz_poly_t B;
    fmpz_poly_init(A);
    fmpz_poly_init(B);
    get_input_polynomials(file, &degA, &degB, A, B);
    fprintf(stderr, "Done.\n\n");
    fmpz_t res;
    fmpz_init(res);
    if(sdeg>0){
      my_dth_subresultant(res, A, B, sdeg, nthreads);
    }
    else{
      myresultant(res, A, B, nthreads);
    }
    //    fmpz_zero(res);
    fmpz_fprint(stdout, res);
    fmpz_poly_clear(A);
    fmpz_poly_clear(B);
    fmpz_clear(res);
    free(input_file);
    return 0;
  }

  fprintf(stdout, "\nRunning resultant computations for random dense polynomials in Z[x]\n\n");
  //  run_modular();
  flint_rand_t state;
  flint_randinit(state);
  run_compare(deg, state);

  //FLINT resultant sur les entiers
  //    run_res_flint(64, state);

  //FLINT resultant modulo premier 16 bits
  //    run_modular16(state, 100);

  //FLINT resultant modulo premier 32 bits
  //    run_modular32(state, 10000);

  //FLINT resultant modulo premier 60 bits
  run_modular60prime(state, 400);

  //FLINT resultant modulo premier fft 60 bits
  //    run_modular60fftprime(state, 100);

  flint_randclear(state);
  return 0;
}
