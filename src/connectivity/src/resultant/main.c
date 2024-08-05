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
#include "flint/mpn_extras.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#include "flint/nmod_vec.h"
#include "flint/nmod_poly.h"
#include "flint/fmpz_mpoly.h"

#define nthreads 1

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define DEBUG 2

#include "data.c"
#include "ratreconstruct.c"
#include "io.c"
#include "flint_euclide.c"

//taken from FLINT
typedef unsigned int UDItype __attribute__ ((mode(DI)));
#define umul_ppmm2(w1, w0, u, v)                \
  __asm__ ("mul{q} %3"                          \
           : "=a" ((UDItype) (w0)),             \
             "=d" ((UDItype) (w1))              \
           : "%0" ((UDItype) (u)),              \
             "rm" ((UDItype) (v)))

#define NMOD_DIVREM_DIVCONQUER_CUTOFF  300


/*
Input: Les coeffs des polynomes sont stockes du plus petit au plus grand degre

pol devient un polynome de degre deg, coeffs de taille bits aleatoires

*/
static void fmpz_poly_random_dense(fmpz_poly_t pol, slong deg, mp_bitcnt_t bits, flint_rand_t state){

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

/*
bivarpol est un tableau de taille degmain + 1 et dont les elements sont
des pointeurs sur tableau de taille degvar + 1 

Les coefficients sont donnes du plus petit au plus grand degre
*/
static void fmpz_bpoly_random_dense(fmpz_poly_t *bivarpol,
                                  slong degmain,
                                  mp_bitcnt_t bits, flint_rand_t state){
  for(slong i = 0; i <= degmain; i++){
    fmpz_poly_random_dense(bivarpol[i], degmain - i, bits, state); 
  }
}

/*

A est un polynome a coefficients dans Z et prime est un nbre premier (qui tient sur un mot-machine)

a_mod contiendra A mod prime

*/
static void make_nmod_poly(nmod_poly_t a_mod, fmpz_poly_t A, mp_limb_t prime){
  a_mod->mod.n = prime;
  a_mod->mod.ninv = n_preinvert_limb(prime);
  flint_clz(a_mod->mod.norm);
  a_mod->length = 0;
  for(unsigned long int i=0; i<A->length; i++){
    nmod_poly_set_coeff_ui(a_mod, i, fmpz_fdiv_ui(A->coeffs+i, prime));
  }
}

//On calcule A mod prime qu'on met dans a (a est suppose deja alloue)
static void make_nmod_bpoly(nmod_poly_t *a, fmpz_poly_t *A, unsigned long int deg, mp_limb_t prime){
  for(unsigned long int i = 0; i <= deg; i++){
    make_nmod_poly(a[i], A[i], prime);
  }
}


static void nmod_bpoly_fast_eval_nmod_2_pols(mp_ptr *a_s_ptr, mp_ptr *b_s_ptr,
                                             nmod_poly_t *a, const slong degA, 
                                             nmod_poly_t *b, const slong degB, 
                                             const mp_ptr eval_points, const ulong n,
                                             mp_ptr * tree){
  nmod_t mod = a[0]->mod;
  _nmod_poly_tree_build(tree, eval_points, n, mod);

  for(ulong i = 0; i <= degA; i++){
    _nmod_poly_evaluate_nmod_vec_fast_precomp(a_s_ptr[i], a[i]->coeffs, a[i]->length, tree, n, mod);
  }
  for(ulong i = 0; i <= degB; i++){
    _nmod_poly_evaluate_nmod_vec_fast_precomp(b_s_ptr[i], b[i]->coeffs, b[i]->length, tree, n, mod);
  }

}


static void nmod_bpoly_eval_nmod(nmod_poly_t a_s,
                                 nmod_poly_t *a, const slong degmain, 
                                 const mp_limb_t eval_point){
  for(ulong i = 0; i <=degmain; i++){
    a_s->coeffs[i] = nmod_poly_evaluate_nmod(a[i], eval_point);
  }
  a_s->length = degmain + 1;
}


/*
res contiendra le resultat
a et b sont des pointeurs sur tableaux de nmod_poly_t de longueur degmain+1
(degmain = degre de la variable principale)

Chaque nmod_poly_t est de degre degvar

numpoints est le nombre de points d'interpolation dont on a besoin ; il doit etre precalcule
 */
static void nmod_bpoly_resultant(nmod_poly_t res, nmod_poly_t *a, const slong degA, nmod_poly_t *b,
                                 const slong degB, const ulong numpoints,
                                 nmod_poly_t *a_s_ptr, nmod_poly_t *b_s_ptr,
                                 mp_ptr restable, mp_ptr eval_points){

 ulong count_points = 0, i = 0;
  while(count_points < numpoints){
    if(nmod_poly_evaluate_nmod(a[degA], i) != 0){
      eval_points[count_points] = i;
      count_points++;
      i++;
    }
    else{
      i++;
    }
  }
  mp_limb_t prime = a[0]->mod.n;

  mp_limb_t ninv = n_preinvert_limb(prime);

  a_s_ptr[0]->mod.n = prime;
  a_s_ptr[0]->mod.ninv = ninv;
  flint_clz(a_s_ptr[0]->mod.norm);
  a_s_ptr[0]->length = 0;

  b_s_ptr[0]->mod.n = prime;
  b_s_ptr[0]->mod.ninv = ninv;
  flint_clz(b_s_ptr[0]->mod.norm);
  b_s_ptr[0]->length = 0;

  //  nmod_poly_init2(a_s, prime, degA+1);
  //  nmod_poly_init2(b_s, prime, degB+1);
  //  mp_limb_t * restable = (mp_limb_t *)(malloc(sizeof(mp_limb_t)*numpoints));
  //  nmod_poly_evaluate_nmod_vec();

  for(i = 0; i < numpoints; i++){
    nmod_bpoly_eval_nmod(a_s_ptr[0], a, degA, eval_points[i]);
    nmod_bpoly_eval_nmod(b_s_ptr[0], b, degB, eval_points[i]);
    restable[i] = nmod_poly_resultant(a_s_ptr[0], b_s_ptr[0]);
  }

  nmod_poly_interpolate_nmod_vec(res, eval_points, restable, numpoints);
#ifdef DEBUG
  char str[1]= {"x"};
  nmod_poly_fprint_pretty(stderr, res, str);fprintf(stderr, "\n");
#endif
  return;
}


static inline void get_coeffs_from_eval(nmod_poly_t a_s, mp_ptr *a_s_ptr, ulong deg, slong n){
  for(ulong i=0; i <= deg; i++){
    a_s->coeffs[i] = a_s_ptr[i][n];
  }
}



//calcul du resultant par evaluation interpolation.
//on suppose que le tableau des points d'evaluation ne contient pas que des 0 et est deja alloue
//les donnees necessaires au calcul sont deja allouees egalement.
static void nmod_bpoly_resultant_allocated_feval(data_heap_t *heap,
                                                 slong mdegA, slong mdegB,
                                                 int pos){
  //a ameliorer
  ulong count_points = 0, i = 0;
  while(count_points < heap->npts){
    if(nmod_poly_evaluate_nmod(heap->a_mod[mdegA], i) != 0 && nmod_poly_evaluate_nmod(heap->b_mod[mdegB], i) != 0){
      heap->ev_pts[count_points] = i;
      count_points++;
      i++;
    }
    else{
      i++;
    }
  }

  nmod_bpoly_fast_eval_nmod_2_pols(heap->a_mod_eval_all, heap->b_mod_eval_all,
                                   heap->a_mod, mdegA, heap->b_mod, mdegB,
                                   heap->ev_pts, heap->npts, heap->tree);

  for(ulong i = 0; i < heap->npts; i++){
    //On ne peut pas paralleliser a cause de w qui est partage. Il faudra modifier ca.

    get_coeffs_from_eval(heap->a_mod_eval, heap->a_mod_eval_all, mdegA, i);
    get_coeffs_from_eval(heap->b_mod_eval, heap->b_mod_eval_all, mdegB, i);
    heap->res_mod_eval[i] = _nmod_poly_resultant_euclidean_allocated(heap->a_mod_eval->coeffs,
                                                           heap->a_mod_eval->length,
                                                           heap->b_mod_eval->coeffs,
                                                           heap->b_mod_eval->length,
                                                           heap->a_mod_eval->mod,
                                                           heap->nmod_coeffs_array);
  }

  nmod_poly_interpolate_nmod_vec(heap->res_mod_array[pos], heap->ev_pts, heap->res_mod_eval, heap->npts);

  return;
}

//Fonction qui renvoie la partie square-free d'un polynome pol

static void square_free_part_nmod_poly(nmod_poly_t pol, nmod_poly_factor_t sqf){

  sqf->num = 0;
  mp_limb_t lc = pol->coeffs[pol->length-1];

  nmod_poly_factor_squarefree(sqf, pol);
  nmod_poly_set(pol, sqf->p);

  for(ulong count = 1; count < sqf->num; count++){
    nmod_poly_mul(pol, pol, sqf->p+count);
  }

  nmod_poly_scalar_mul_nmod(pol, pol, lc);
  sqf->num = 0;
}

static inline void nmod_poly_init_copy(nmod_poly_t *poly1, nmod_poly_t *poly2){
  (*poly1)->length = (*poly2)->length;
  (*poly1)->alloc = (*poly2)->alloc;
  (*poly1)->coeffs = (*poly2)->coeffs;
}

static inline void nmod_poly_copy(nmod_poly_t poly1, nmod_poly_t poly2){
  (poly1)->length = (poly2)->length;
  for(ulong i = 0; i < (poly1)->length; i++){
    (poly1)->coeffs[i] = (poly2)->coeffs[i];
  }
}

//calcul d'une param rationnelle par evaluation interpolation.
//on suppose que le tableau des points d'evaluation ne contient pas que des 0 et est deja alloue
//les donnees necessaires au calcul sont deja allouees egalement.
//On renvoie 0 si il y a un probleme
//Notamment, on suppose que le resultant est square-free.
//Du coup, une seule parametrisation est construite dans nmod_all_params_table+pos
static int nmod_bpoly_param_allocated_feval_sqfree(data_heap_t *heap,
                                                   slong mdegA, slong mdegB,
                                                   int pos){

  ulong count_points = 0, i = 0;
  while(count_points < heap->npts){
    if(nmod_poly_evaluate_nmod(heap->a_mod[mdegA], i) != 0 &&
       nmod_poly_evaluate_nmod(heap->b_mod[mdegB], i) != 0){
      heap->ev_pts[count_points] = i;
      count_points++;
      i++;
    }
    else{
      i++;
    }
  }

  nmod_bpoly_fast_eval_nmod_2_pols(heap->a_mod_eval_all,
                                   heap->b_mod_eval_all,
                                   heap->a_mod, mdegA,
                                   heap->b_mod, mdegB,
                                   heap->ev_pts,
                                   heap->npts, heap->tree);
  for(ulong i = 0; i < heap->npts; i++){
    //On ne peut pas paralleliser a cause de w qui est partage. Il faudra modifier ca.

    get_coeffs_from_eval(heap->a_mod_eval, heap->a_mod_eval_all, mdegA, i);
    get_coeffs_from_eval(heap->b_mod_eval, heap->b_mod_eval_all, mdegB, i);

    if(DEBUG > 3){
      fprintf(stderr, "ev_pt = %ld\n", heap->ev_pts[i]);
      char str[1] = {"x"};
      nmod_poly_fprint_pretty(stderr, heap->a_mod_eval, str);
      fprintf(stderr, "\n");
      nmod_poly_fprint_pretty(stderr, heap->b_mod_eval, str);
      fprintf(stderr, "\n");
    }
    heap->res_mod_eval[i] =
      _nmod_poly_param_euclidean_allocated(heap->a_mod_eval->coeffs,
                                           heap->a_mod_eval->length,
                                           heap->b_mod_eval->coeffs,
                                           heap->b_mod_eval->length,
                                           heap->a_mod_eval->mod,
                                           heap->nmod_coeffs_array,
                                           heap->par_h_mod_eval + i,
                                           heap->par_t_mod_eval + i);
    if(DEBUG > 3){
      fprintf(stderr, "Computed resultant = %ld\n",
              heap->res_mod_eval[i]);
      fprintf(stderr, "head -> %ld\n",
              heap->par_h_mod_eval[i]);
      fprintf(stderr, "tail -> %ld\n",
              heap->par_t_mod_eval[i]);
    }
  }

  nmod_poly_interpolate_nmod_vec(heap->nmod_lparams_array[pos]->params[0]->elim,
                                 heap->ev_pts, heap->res_mod_eval, heap->npts);

  nmod_poly_interpolate_nmod_vec(heap->nmod_lparams_array[pos]->params[0]->denom,
                                 heap->ev_pts, heap->par_h_mod_eval, heap->npts);

  nmod_poly_interpolate_nmod_vec(heap->nmod_lparams_array[pos]->params[0]->coords[0],
                                 heap->ev_pts, heap->par_t_mod_eval, heap->npts);
  //On suppose que le resultant est square-free.

  if(DEBUG > 3){
    for(ulong k = 0; k < heap->npts; k++){
      fprintf(stderr, "%lu ", heap->par_t_mod_eval[k]);
    }
    fprintf(stderr, "\n");
    for(ulong k = 0; k < heap->npts; k++){
      fprintf(stderr, "%lu ", heap->ev_pts[k]);
    }
    fprintf(stderr, "\n");
    char str[1] = {"x"};
    fprintf(stderr, "degree = %ld\n", heap->par_t_mod_array[pos]->length - 1);
    nmod_poly_fprint_pretty(stderr, heap->par_t_mod_array[pos], str);
    fprintf(stderr,"\n");
    fprintf(stderr, "degree = %ld\n", heap->par_h_mod_array[pos]->length - 1);
    nmod_poly_fprint_pretty(stderr, heap->par_h_mod_array[pos], str);
    fprintf(stderr,"\n\n");
    nmod_poly_fprint_pretty(stderr, heap->nmod_lparams_array[pos]->params[0]->coords[0], str);
    fprintf(stderr,"\n");
  }

  if(DEBUG > 4){
    fprintf(stderr, "\nParametrization mod %ld\n", heap->primes_array[pos]);
    print_nmod_param(stderr, heap->nmod_lparams_array[pos]->params[0]);
    fprintf(stderr, "prime in param struct = \n%ld, \n%ld\n",
            heap->nmod_lparams_array[pos]->params[0]->denom->mod.n,
            heap->nmod_lparams_array[pos]->params[0]->coords[0]->mod.n);
  }

  return 1;
}


//calcul d'une param rationnelle par evaluation interpolation.
//on suppose que le tableau des points d'evaluation ne contient pas que des 0 et est deja alloue
//les donnees necessaires au calcul sont deja allouees egalement.
//On renvoie 0 si il y a un probleme
//Ici, on suppose que le resultant n'est pas square-free.

static int nmod_bpoly_param_allocated_feval_nonsqfree(data_heap_t *heap,
                                                      slong mdegA, slong mdegB,
                                                      int pos){

  ulong count_points = 0, i = 0;
  while(count_points < heap->npts){
    if(nmod_poly_evaluate_nmod(heap->a_mod[mdegA], i) != 0 &&
       nmod_poly_evaluate_nmod(heap->b_mod[mdegB], i) != 0){
      heap->ev_pts[count_points] = i;
      count_points++;
      i++;
    }
    else{
      i++;
    }
  }

  nmod_bpoly_fast_eval_nmod_2_pols(heap->a_mod_eval_all, heap->b_mod_eval_all,
                                   heap->a_mod, mdegA, heap->b_mod, mdegB,
                                   heap->ev_pts, heap->npts, heap->tree);
  for(ulong i = 0; i < heap->npts; i++){
    //On ne peut pas paralleliser a cause de w qui est partage. Il faudra modifier ca.

    get_coeffs_from_eval(heap->a_mod_eval, heap->a_mod_eval_all, mdegA, i);
    get_coeffs_from_eval(heap->b_mod_eval, heap->b_mod_eval_all, mdegB, i);

    heap->res_mod_eval[i] =
      _nmod_poly_param_euclidean_allocated(heap->a_mod_eval->coeffs,
                                           heap->a_mod_eval->length,
                                           heap->b_mod_eval->coeffs,
                                           heap->b_mod_eval->length,
                                           heap->a_mod_eval->mod,
                                           heap->nmod_coeffs_array,
                                           heap->par_h_mod_eval + i,
                                           heap->par_t_mod_eval + i);

  }

  nmod_poly_interpolate_nmod_vec(heap->res_mod_array[pos], heap->ev_pts, heap->res_mod_eval, heap->npts);
  nmod_poly_interpolate_nmod_vec(heap->par_h_mod_array[pos],
                                 heap->ev_pts, heap->par_h_mod_eval, heap->npts);
  nmod_poly_interpolate_nmod_vec(heap->par_t_mod_array[pos],
                                 heap->ev_pts, heap->par_t_mod_eval, heap->npts);

  heap->sqf_part->num = 0;
  nmod_poly_factor_squarefree(heap->sqf_part, heap->res_mod_array[pos]);

  for(ulong i = 0; i < heap->sqf_part->num; i++){
    nmod_poly_copy(heap->nmod_lparams_array[pos]->params[i]->elim,
                     heap->sqf_part->p+i);

    if(nmod_poly_invmod(heap->GCDINV[pos],
                        heap->par_h_mod_array[pos],
                        heap->nmod_lparams_array[pos]->params[i]->elim)){

      nmod_poly_derivative(heap->nmod_lparams_array[pos]->params[i]->denom,
                     heap->nmod_lparams_array[pos]->params[i]->elim);

      nmod_poly_mul(heap->nmod_lparams_array[pos]->params[i]->coords[0],
                    heap->par_t_mod_array[pos],
                    heap->GCDINV[pos]);

      nmod_poly_rem(heap->nmod_lparams_array[pos]->params[i]->coords[0],
                    heap->nmod_lparams_array[pos]->params[i]->coords[0],
                    heap->nmod_lparams_array[pos]->params[i]->elim);
      nmod_poly_mul(heap->nmod_lparams_array[pos]->params[i]->coords[0],
                    heap->nmod_lparams_array[pos]->params[i]->coords[0],
                    heap->nmod_lparams_array[pos]->params[i]->denom
                    );

      nmod_poly_rem(heap->nmod_lparams_array[pos]->params[i]->coords[0],
                    heap->nmod_lparams_array[pos]->params[i]->coords[0],
                    heap->nmod_lparams_array[pos]->params[i]->elim);
    }
    else{
      fprintf(stderr, "Pas d'inversion possible\n");
      return 0;
    }
  }

  if(DEBUG>2){
    fprintf(stderr, "Prime = %ld\n", heap->primes_array[pos]);
    print_nmod_all_params_pretty(stderr, heap->nmod_lparams_array[pos]);
  }

  return 1;
}

static void fmpz_bpoly_resultant_deterministic(fmpz_poly_t res,
                                        fmpz_poly_t *A, slong degA,
                                        fmpz_poly_t *B, slong degB){
  slong totaldegA = get_total_degree(A, degA);
  slong totaldegB = get_total_degree(B, degB);

  fprintf(stderr, "Total degrees : %ld %ld\n", totaldegA, totaldegB);
  //maximum degree in the variable that is not eliminated
  slong degAx = get_partial_degree(A, degA);
  slong degBx = get_partial_degree(B, degB);

  ulong numpoints = MIN(totaldegA * totaldegB + 1, degB * degAx + degA * degBx + 1);
  ulong nbitsA = nbits_bound_fmpz_bvar_poly(A, degA);
  ulong nbitsB = nbits_bound_fmpz_bvar_poly(B, degB);
  fprintf(stderr, "Max nbits : %lu %lu\n", nbitsA, nbitsB);
  ulong i;

  //Hadamard-type bound (Goldstein-Graham)
  ulong bound_nbits = nbitsA*(degB) + nbitsB*(degA);
  fprintf(stderr, "Bound (nbits) = %lu\n", bound_nbits);

  mp_limb_t prime = 1073741827;
  prime = n_nextprime(1152921504606847009, 0);
  //  prime = n_nextprime(prime, 0);
  ulong nbits_prime = LOG2(prime);

  nmod_poly_t *a = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * (degA + 1)));
  nmod_poly_t *b = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * (degB + 1)));

  for(i = 0; i <= degA; i++){
    nmod_poly_init(a[i], prime);
    make_nmod_poly(a[i], A[i], prime);
  }
  for(i = 0; i <= degB; i++){
    nmod_poly_init(b[i], prime);
    make_nmod_poly(b[i], B[i], prime);
  }

  //data for calling nmod_bpoly_resultant
  //a_s_ptr contiendra la specialisation des coefficients de a en un point. 
  nmod_poly_t *a_s_ptr = malloc(sizeof(nmod_poly_t));
  nmod_poly_init2(a_s_ptr[0], prime, degA+1);

  nmod_poly_t *b_s_ptr = malloc(sizeof(nmod_poly_t));
  nmod_poly_init2(b_s_ptr[0], prime, degB+1);

  //contiendra les valeurs des resultants modulo des premiers
  mp_limb_t * restable = (mp_limb_t *)(malloc(sizeof(mp_limb_t)*numpoints));
  //contiendra les points d'evaluation
  mp_limb_t * eval_points = (mp_limb_t *)(malloc(sizeof(mp_limb_t)*numpoints));
  //////////////////////////////////

  ulong count_nprimes = 1;

  nmod_poly_t *res_mod_table = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * (bound_nbits / nbits_prime + 1) ));
  mp_limb_t *primes_table = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * (bound_nbits / nbits_prime + 1)));
  double e = omp_get_wtime();
  nbits_prime = 0;

  fmpz_t *prod_mod = malloc(sizeof(fmpz_t));
  fmpz_init(prod_mod[0]);
  while(nbits_prime < bound_nbits){
    for(i = 0; i <= degA; i++){
      make_nmod_poly(a[i], A[i], prime);
    }
    for(i = 0; i <= degB; i++){
      make_nmod_poly(b[i], B[i], prime);
    }
    if(a[degA]->length>0 && b[degB]->length>0){
      nmod_poly_init(res_mod_table[count_nprimes-1], prime);
      primes_table[count_nprimes-1] = prime;
      nmod_bpoly_resultant(res_mod_table[count_nprimes-1], a, degA, b, degB, numpoints,
                           a_s_ptr, b_s_ptr, restable, eval_points);
      prime = n_nextprime(prime, 0);
      nbits_prime += LOG2(prime);
      count_nprimes++;
    }
    else{
      prime = n_nextprime(prime, 0);
    }
  }
  count_nprimes--;
  fprintf(stderr, "Number of used primes = %lu\n", count_nprimes);
  fprintf(stderr, "Elapsed time (mod comp.): %f\n", omp_get_wtime() - e);

  e=omp_get_wtime();
  crt_lift_resultant(res, res_mod_table, primes_table, count_nprimes, prod_mod);
  fprintf(stderr, "Elapsed time (CRT comp.): %f\n", omp_get_wtime() - e);
  fprintf(stderr, "Final size coeffs: %ld\n", FLINT_ABS(_fmpz_vec_max_bits(res->coeffs, res->length)));
  fprintf(stderr, "Degree of the output: %ld\n", res->length-1);

  fmpz_clear(prod_mod[0]);
  free(prod_mod);
  nmod_poly_clear(a_s_ptr[0]);
  nmod_poly_clear(b_s_ptr[0]);
  free(a_s_ptr);
  free(b_s_ptr);
  free(eval_points);

  for(i = 0; i< count_nprimes; i++){
    nmod_poly_clear(res_mod_table[i]);
  }
  for(i = 0; i <= degA; i++){
    nmod_poly_clear(a[i]);
  }
  for(i = 0; i <= degB; i++){
    nmod_poly_clear(b[i]);
  }
  free(res_mod_table);
  free(restable);
  free(primes_table);
  free(a);
  free(b);
  return;
}

/*
  Input: dans res_mod_table on a N = size_table images modulaires de resultant
  (modulo les premiers de primes_table).

  On vérifie que res modulo ces premiers correspond aux elements de res_mod_table
*/
/* static int check_modular_images_old(fmpz_poly_t res, nmod_poly_t res_mod, nmod_poly_t *res_mod_table, */
/*                          mp_limb_t *primes_table, int size_table, slong *max_degree){ */
/*   for(int i = 0; i < size_table; i++){ */
/*     if(*max_degree < res_mod_table[i]->length - 1){ */
/*       *max_degree = res_mod_table[i]->length - 1; */
/*     } */
/*   } */
/*   for(int i = 0; i < size_table; i++){ */
/*     make_nmod_poly(res_mod, res, primes_table[i]); */
/*     if(!nmod_poly_equal(res_mod, res_mod_table[i])){ */
/*       return 0; */
/*     } */
/*   } */
/*   return 1; */
/* } */

/*
  Input: dans res_mod_table on a N = size_table images modulaires de resultant
  (modulo les premiers de primes_table).

  On vérifie que res modulo ces premiers correspond aux elements de res_mod_table
*/
static int check_modular_images(fmpz_poly_t res,
                                data_heap_t *heap,
                                int size_table, slong *max_degree){
  for(int i = 0; i < size_table; i++){
    if(*max_degree < heap->res_mod_array[i]->length - 1){
      *max_degree = heap->res_mod_array[i]->length - 1;
    }
  }
  for(int i = 0; i < size_table; i++){
    make_nmod_poly(heap->res_mod[0], res, heap->primes_array[i]);
    if(!nmod_poly_equal(heap->res_mod[0], heap->res_mod_array[i])){
      return 0;
    }
  }
  return 1;
}


static void update_poly_tables_with_prime(data_heap_t *heap,
                                          int pos, mp_limb_t prime){
  mp_limb_t p = n_preinvert_limb(prime);
  heap->primes_array[pos] = prime;

  heap->GCDINV[pos]->mod.n = prime;
  heap->GCDINV[pos]->mod.ninv = p;
  flint_clz(heap->GCDINV[pos]->mod.norm);

  heap->res_mod_array[pos]->mod.n = prime;
  heap->res_mod_array[pos]->mod.ninv = p;
  flint_clz(heap->res_mod_array[pos]->mod.norm);

  heap->par_t_mod_array[pos]->mod.n = prime;
  heap->par_t_mod_array[pos]->mod.ninv = p;
  flint_clz(heap->par_t_mod_array[pos]->mod.norm);

  heap->par_h_mod_array[pos]->mod.n = prime;
  heap->par_h_mod_array[pos]->mod.ninv = p;
  flint_clz(heap->par_h_mod_array[pos]->mod.norm);

  heap->a_mod_eval->mod.n = prime;
  heap->a_mod_eval->mod.ninv = p;
  flint_clz(heap->a_mod_eval->mod.norm);

  heap->b_mod_eval->mod.n = prime;
  heap->b_mod_eval->mod.ninv = p;
  flint_clz(heap->b_mod_eval->mod.norm);

}

static void update_nmod_param_with_prime(nmod_all_params_t *lparams,
                                         mp_limb_t prime){
  mp_limb_t p = n_preinvert_limb(prime);

  for(int i = 0; i <lparams->nparams; i++){
    nmod_param_t *param = lparams->params[i];
    param->denom->mod.n = prime;
    param->denom->mod.ninv = p;
    flint_clz(param->denom->mod.norm);

    param->elim->mod.n = prime;
    param->elim->mod.ninv = p;
    flint_clz(param->elim->mod.norm);

    for(int j = 0; j < param->nvars; j++){
      param->coords[j]->mod.n = prime;
      param->coords[j]->mod.ninv = p;
      flint_clz(param->coords[j]->mod.norm);
    }
  }
}


//Fonction qui teste si pol est square-free
//renvoie 0 si pol est square-free, un entier >= 1 sinon
static slong is_non_square_free(nmod_poly_t pol, nmod_poly_factor_t sqf){
  nmod_poly_factor_squarefree(sqf, pol);
  int nf =   sqf->num ;
  sqf->num = 0;
  return nf - 1;
}


static void get_degrees_and_nbits(fmpz_poly_t *A, slong mdegA,
                                  fmpz_poly_t *B, slong mdegB,
                                  slong *tdegA, slong *tdegB,
                                  slong *sdegA, slong *sdegB,
                                  ulong *nbitsA, ulong *nbitsB){
  //total degrees
  *tdegA = get_total_degree(A, mdegA);
  *tdegB = get_total_degree(B, mdegB);

  //maximum degree in the variable that is not eliminated
  *sdegA = get_partial_degree(A, mdegA);
  *sdegB = get_partial_degree(B, mdegB);

  *nbitsA = nbits_bound_fmpz_bvar_poly(A, mdegA);
  *nbitsB = nbits_bound_fmpz_bvar_poly(B, mdegB);
}

//renvoie 0 si on a encore besoin de lifter les coefficients
static inline int check_param_modular_images(fmpz_param_t *param,
                                             data_heap_t *heap, int size){
  for(ulong i = 0; i < size; i++){
    //on verifie le polynome eliminant
    make_nmod_poly(heap->res_mod[0], param->elim, heap->primes_array[i]);
    if(!nmod_poly_equal(heap->res_mod[0],
                        heap->nmod_lparams_array[i]->params[0]->elim)){
      return 0;
    }
    make_nmod_poly(heap->res_mod[0], param->denom, heap->primes_array[i]);
    //on verifie le polynome denominateur
    if(!nmod_poly_equal(heap->res_mod[0],
                        heap->nmod_lparams_array[i]->params[0]->denom)){
      return 0;
    }
    //on verifie les coordonnees
    make_nmod_poly(heap->res_mod[0], param->coords[0], heap->primes_array[i]);
    if(!nmod_poly_equal(heap->res_mod[0],
                        heap->nmod_lparams_array[i]->params[0]->coords[0])){
      return 0;
    }
  }
  return 1;
}


//renvoie 0 si on a encore besoin de lifter les coefficients
static inline int check_all_params_modular_image(fmpz_all_params_t *all_params,
                                             data_heap_t *heap, int size){
  int nparams = all_params->nparams;
  for(int i = 0; i < nparams; i++){
    if(!check_param_modular_images(all_params->params[i], heap, size)){
      return 0;
    }
  }
  return 1;
}

static inline void integer_param_crt_lifting_all_params(fmpz_all_params_t *all_params,
                                         data_heap_t *heap, int size){
  slong maxdegree = -1;
  fmpz_param_t *param = all_params->params[0];
  for(ulong i = 0; i < size; i++){
    if(heap->nmod_lparams_array[0]->params[0]->elim->length - 1 > maxdegree){
      maxdegree = heap->nmod_lparams_array[0]->params[0]->elim->length - 1;
    }
  }
  if((param)->elim->length-1!=maxdegree){
    fprintf(stderr, "Error with degree estimates in integer_param_crt_lifting\n");
    fprintf(stderr, "maxdegree = %ld\n", maxdegree);
    fprintf(stderr, "First param of degree = %ld\n", (param)->elim->length-1);
    exit(1);
  }

  if(DEBUG>=5){
    for(int i = 0; i < size; i++){
      fprintf(stderr, "prime = %ld\n", heap->primes_array[i]);
      print_nmod_param(stderr, (heap->nmod_lparams_array[i]->params[0]));
      fprintf(stderr, "\n");
    }
  }

  for(ulong j = 0; j <= maxdegree; j++){
    for(ulong i = 0; i < size; i++){

      fmpz_CRT_ui(((param)->elim->coeffs)+j, ((param)->elim->coeffs)+j,
                  heap->prod_primes[i],
                  (heap->nmod_lparams_array[i]->params[0]->elim->coeffs)[j],
                  heap->primes_array[i], 1);

      fmpz_CRT_ui(((param)->denom->coeffs)+j, ((param)->denom->coeffs)+j,
                  heap->prod_primes[i],
                  (heap->nmod_lparams_array[i]->params[0]->denom->coeffs)[j],
                  heap->primes_array[i], 1);

      fmpz_CRT_ui(((param)->coords[0]->coeffs)+j, ((param)->coords[0]->coeffs)+j,
                  heap->prod_primes[i],
                  (heap->nmod_lparams_array[i]->params[0]->coords[0]->coeffs)[j],
                  heap->primes_array[i], 1);

    }
  }

}

static inline int rational_param_lifting(fmpz_param_t *param,
                                          data_heap_t *heap,
                                          int size, int mul){
  int b = 1;

  ulong deg = heap->nmod_lparams_array[0]->params[mul]->elim->length - 1;
  if(param->rec_elim < deg){
    for(ulong j = 0; j <= deg; j++){
      for(ulong i = 0; i < size; i++){
        fmpz_CRT_ui(((param)->elim->coeffs)+j, ((param)->elim->coeffs)+j,
                    heap->prod_primes[i],
                    (heap->nmod_lparams_array[i]->params[mul]->elim->coeffs)[j],
                    heap->primes_array[i], 1);
      }
    }
    fprintf(stderr, "Lifting the elimination polynomial\n");
    if(rat_recon(param->elim, param->elim_num, param->elim_den, param->elim_lcm,
                 heap, deg, &(param->rec_elim)) == 0){
      b = 0;
    }
  }

  deg = (heap->nmod_lparams_array[0]->params[mul]->coords[0]->length) - 1;
  if(param->rec_coords[0] < deg){
    for(ulong j = 0; j <= deg; j++){
      for(ulong i = 0; i < size; i++){
        fmpz_CRT_ui(((param)->coords[0]->coeffs)+j, ((param)->coords[0]->coeffs)+j,
                    heap->prod_primes[i],
                    (heap->nmod_lparams_array[i]->params[mul]->coords[0]->coeffs)[j],
                    heap->primes_array[i], 1);
      }
    }
    if(b==0){
      return b;
    }
    fprintf(stderr, "Lifting of the parametrization\n");
    if(rat_recon(param->coords[0], param->coords_num[0], param->coords_den[0],
                 param->coords_lcm[0],
                 heap, deg, param->rec_coords) == 0){
      b = 0;
    }
  }
  fprintf(stderr, "LCM = ");
  fmpz_fprint(stderr, param->coords_lcm[0]);
  fprintf(stderr, "\n");

  return b;
}

//attention : ici le denominateur sera de tout de maniere la derivee du polynome
//eliminant
static inline int rational_param_lifting_all_params(fmpz_all_params_t *all_params,
                                         data_heap_t *heap, int size){
  int b = 1;
  for(ulong i = 0; i < all_params->nparams; i++){
    fmpz_param_t *param = all_params->params[i];
    if(param->recons == 0){
      if(rational_param_lifting(param, heap, size, i)==0){
        b = 0;
      }
      else{
        param->recons = 1;
      }
    }
  }
  return b;
}


static inline void copy_coeffs_poly(fmpz_poly_t poly1, fmpz_poly_t poly2){

  for(ulong i = 0; i < poly2->length; i++){
    fmpz_poly_set_coeff_fmpz(poly1, i, poly2->coeffs + i);
  }
}

static inline bsolve_fmpz_param_t *copy_to_final_parametrization(fmpz_param_t *param){
  bsolve_fmpz_param_t *bs_param = malloc(sizeof(bsolve_fmpz_param_t));
  bs_param->nvars = param->nvars;
  bs_param->sepelem = malloc(sizeof(mp_limb_t) * 2);
  bs_param->cfs = malloc(sizeof(fmpz_t) * param->nvars);

  for(int i = 0; i < 2; i++){
    bs_param->sepelem[i] = param->sepelem[i];
  }
  fmpz_poly_init2(bs_param->elim, param->elim_num->length);

  copy_coeffs_poly(bs_param->elim, param->elim_num);

  bs_param->coords = malloc(sizeof(fmpz_poly_t) * param->nvars);

  for(int i = 0; i < param->nvars; i++){
    fmpz_poly_init2(bs_param->coords[i], param->coords_num[i]->length);
    copy_coeffs_poly(bs_param->coords[i], param->coords_num[i]);
    fmpz_init(bs_param->cfs[i]);
    fmpz_set(bs_param->cfs[i], param->coords_lcm[i]);
    fprintf(stderr, "Warning: changing cfs\n");
    fmpz_set_ui(bs_param->cfs[i], 1);
  }
  fmpz_poly_init2(bs_param->denom, param->elim_num->length - 1);

  for(ulong i = 0; i < param->elim_num->length - 1; i++){
    fmpz_poly_set_coeff_ui(bs_param->denom, i, 1);
    fmpz_mul_ui(bs_param->denom->coeffs + i, bs_param->elim->coeffs + (i+1), i+1);
  }
  return bs_param;
}

static inline void copy_to_final_all_params(bsolve_fmpz_all_params_t *bs_params,
                                            fmpz_all_params_t *params){
  bs_params->nparams = params->nparams;
  fprintf(stderr, "Number of parametrizations = %d\n", params->nparams);
  bs_params->params = malloc(sizeof(bsolve_fmpz_param_t *) * params->nparams);
  for(int i = 0; i < params->nparams; i++){
    bs_params->params[i] = copy_to_final_parametrization(params->params[i]);
  }
}

//retourne 1 si il y a une param. rat. des solutions
//sinon retourne 0

/* output parametrization will be stored in bs_param */
static int _fmpz_bpoly_param_probabilistic(bsolve_fmpz_all_params_t *bs_param,
                                           fmpz_poly_t *A, slong mdegA,
                                           fmpz_poly_t *B, slong mdegB){
  slong tdegA, tdegB, sdegA, sdegB;
  ulong nbitsA, nbitsB;
  int size_table = 8;
  /*
   tdegA and degB will be the total degrees of A and B

   sdegA and sdegB will be the partial degrees of A and B
   in the variable that is kept.

   nbitsA and nbitsB are the max bit size coeffs of A and B
   */
  get_degrees_and_nbits(A, mdegA,
                        B, mdegB,
                        &tdegA, &tdegB,
                        &sdegA, &sdegB,
                        &nbitsA, &nbitsB);

  data_heap_t *heap = allocate_data_heap(tdegA, tdegB,
                                         mdegA, mdegB,
                                         sdegA, sdegB,
                                         size_table);

  fprintf(stderr, "Total degrees : %ld %ld\n", tdegA, tdegB);
  fprintf(stderr, "Main degrees : %ld %ld\n", mdegA, mdegB);
  fprintf(stderr, "Bit sizes : %lu %lu\n", nbitsA, nbitsB);
  fprintf(stderr, "numpoints = %lu\n", heap->npts);

  /*Hadamard-type (Goldstein-Graham) to obtain a bound on the
  bit size of the resultant
  */
  ulong bound_nbits = nbitsA*(mdegB) + nbitsB*(mdegA);

  fprintf(stderr, "Bound (nbits) = %lu\n", bound_nbits);

  mp_limb_t prime = 1073741827;
  prime = n_nextprime(1152921504606847009, 0);
  ulong nbits_prime = LOG2(prime);

  ulong count_nprimes = 1;

  double e = omp_get_wtime();
  nbits_prime = 0;


  double e_mod = 0;

  //On teste si le resultant est square-free
  make_nmod_bpoly(heap->a_mod, A, mdegA, prime);
  make_nmod_bpoly(heap->b_mod, B, mdegB, prime);

  update_poly_tables_with_prime(heap, 0, prime);
  nmod_bpoly_param_allocated_feval_sqfree(heap,
                                          mdegA, mdegB, 0);

  fprintf(stderr, "Degree of the resultant is %ld\n", heap->nmod_lparams_array[0]->params[0]->elim->length - 1);

  /* resultant is a non-zero constant -> no solution */
  if(heap->nmod_lparams_array[0]->params[0]->elim->length == 1){
#if(DEBUG > 3)
    char str[1] = {"x"};
    fprintf(stderr, "elim = ");
    nmod_poly_fprint_pretty(stderr,
                            heap->nmod_lparams_array[0]->params[0]->elim,
                            str);
    fprintf(stderr, "\n");
    fprintf(stderr, "denom = ");
    nmod_poly_fprint_pretty(stderr,
                            heap->nmod_lparams_array[0]->params[0]->denom,
                            str);
    fprintf(stderr, "\n");
#endif
    bs_param->nparams = 1;
    bs_param->params = malloc(sizeof(bsolve_fmpz_param_t *));
    bs_param->params[0] = malloc(sizeof(bsolve_fmpz_param_t));
    bs_param->params[0]->nvars = 1;

    bs_param->params[0]->sepelem = malloc(sizeof(mp_limb_t) * 2);
    bs_param->params[0]->sepelem[0] = 0;
    bs_param->params[0]->sepelem[1] = 1;

    bs_param->params[0]->coords = malloc(sizeof(fmpz_poly_t) );

    fmpz_poly_init2(bs_param->params[0]->elim, 1);
    fmpz_poly_set_ui(bs_param->params[0]->elim, 1);
    fmpz_poly_init2(bs_param->params[0]->coords[0], 1);
    fmpz_poly_set_ui(bs_param->params[0]->coords[0], 0);
    return 1;
  }
  if(DEBUG>3){
    char str[1] = {"x"};
    nmod_poly_fprint_pretty(stderr, *heap->res_mod_array, str);
    fprintf(stderr, "\n");
    fprintf(stderr, "RESULTANT = ");
    nmod_poly_fprint_pretty(stderr, heap->nmod_lparams_array[0]->params[0]->elim, str);
    fprintf(stderr, "\n");
  }

  heap->sqf_part->num = 0;
  nmod_poly_factor_squarefree(heap->sqf_part, heap->nmod_lparams_array[0]->params[0]->elim);

  int nb_sqfree = heap->sqf_part->num;
  int is_sqf = 0;
  if(nb_sqfree==1 && (heap->sqf_part->p)->length == heap->nmod_lparams_array[0]->params[0]->elim->length)
    {
    fprintf(stderr, "Resultant is square-free\n");
    is_sqf = 1;
    }
  else{
    fprintf(stderr, "Resultant is not square-free\n");
    fprintf(stderr, "Array of degree/multiplicity is [");
    for(int i = 0; i < nb_sqfree-1; i++){
      fprintf(stderr, "[%lu, %lu], ", (heap->sqf_part->p+i)->length - 1, (heap->sqf_part->exp[i]));
    }
    fprintf(stderr, "[%lu, %lu]]\n", (heap->sqf_part->p+(nb_sqfree-1))->length - 1, (heap->sqf_part->exp[nb_sqfree-1]));
  }

  if(nmod_poly_invmod(heap->GCDINV[0],
                      heap->nmod_lparams_array[0]->params[0]->denom,
                      heap->nmod_lparams_array[0]->params[0]->elim)
     == 0){

    fprintf(stderr, "The input system is not in generic position\n");

    return 0;
  }
  mp_ptr sepelem = malloc(sizeof(mp_limb_t) * 2);
  sepelem[0]=0;
  sepelem[1]=1;

  for(int i = 0; i < size_table; i++){
    heap->nmod_lparams_array[i] = allocate_degrees_sepelem_nmod_all_params(heap, sepelem,
                                                                           nb_sqfree, 1,
                                                                           prime);
  }

  fmpz_all_params_t *param = malloc(sizeof(fmpz_all_params_t));
  allocate_degrees_sepelem_fmpz_all_params(param, heap, sepelem,
                                           nb_sqfree, 1);
  fprintf(stderr, "Allocation done\n");
  fprintf(stderr, "First param of degree = %ld\n", param->params[0]->elim->length-1);

  prime = n_nextprime(prime, 0);
  int boo = 1;
  while(boo){
    int pos = 0;

    while(pos < size_table){

      make_nmod_bpoly(heap->a_mod, A, mdegA, prime);
      make_nmod_bpoly(heap->b_mod, B, mdegB, prime);

      if(heap->a_mod[mdegA]->length>0 && heap->b_mod[mdegB]->length>0){

        double start_mod = omp_get_wtime();
        update_poly_tables_with_prime(heap, pos, prime);
        //Ici on met a jour nmod_all_primes_table avec prime
        update_nmod_param_with_prime(heap->nmod_lparams_array[pos], prime);

        if(is_sqf==1){
          nmod_bpoly_param_allocated_feval_sqfree(heap,
                                                  mdegA, mdegB, pos);
        }
        else{
          if(nmod_bpoly_param_allocated_feval_nonsqfree(heap,
                                                        mdegA, mdegB, pos)==0){
            return 0;
          }
        }

#if(DEBUG>2)
        char str[1] = {"x"};
        fprintf(stderr, "prime = %ld\n", prime);
        fprintf(stderr, "RESULTANT = ");
        nmod_poly_fprint_pretty(stderr, heap->nmod_lparams_array[pos]->params[0]->elim, str);
        fprintf(stderr, "\n");
#endif
        e_mod += (omp_get_wtime() - start_mod);
        pos++;
        fmpz_mul_ui(heap->prod_primes[pos], heap->prod_primes[pos - 1], prime);

        nbits_prime += LOG2(prime);
        prime = n_nextprime(prime, 0);
        count_nprimes++;
      }
      else{
        prime = n_nextprime(prime, 0);
      }
    }

    if(nb_sqfree==1){
      if(check_all_params_modular_image(param, heap, size_table) == 0){
        //Il faut continuer la construction
        integer_param_crt_lifting_all_params(param, heap,
                                             size_table);
        fmpz_set(heap->prod_primes[0], heap->prod_primes[size_table]);
        if(DEBUG>1){
          if(check_all_params_modular_image(param, heap, size_table)==0){
            fprintf(stderr, "Error in lifting\n");
            exit(1);
          }
        }
      }
      else{
        boo = 0;
      }
    }
    else{
      if(rational_param_lifting_all_params(param, heap, size_table) == 0){
        fmpz_set(heap->prod_primes[0], heap->prod_primes[size_table]);
      }
      else{
        boo = 0;
      }
    }
  }
  copy_to_final_all_params(bs_param, param);
  count_nprimes--;
  fprintf(stderr, "Number of used primes = %lu\n", count_nprimes);
  fprintf(stderr, "Elapsed time (modp): %f\n", e_mod);
  fprintf(stderr, "Elapsed time (total): %f\n", omp_get_wtime() - e);

  free(sepelem);
//  free_fmpz_all_params(all_params);
  free_data_heap(mdegA, mdegB, size_table, heap);

  return 1;
}

static int fmpz_bpoly_param(bsolve_fmpz_all_params_t *bs_param, 
                            fmpz_poly_t *A, slong mdegA,
                            fmpz_poly_t *B, slong mdegB){
  if(mdegA < mdegB){
    return _fmpz_bpoly_param_probabilistic(bs_param,
                                           B, mdegB,
                                           A, mdegA);
  }
  else{
    return _fmpz_bpoly_param_probabilistic(bs_param,
                                           A, mdegA,
                                           B, mdegB);
  }
}



static void _fmpz_bpoly_resultant_probabilistic(fmpz_poly_t res,
                                                fmpz_poly_t *A, slong mdegA,
                                                fmpz_poly_t *B, slong mdegB){
  slong tdegA, tdegB, sdegA, sdegB;
  ulong nbitsA, nbitsB;
  int size_table = 4;
  fprintf(stderr, "ici on est la\n");
  get_degrees_and_nbits(A, mdegA,
                        B, mdegB,
                        &tdegA, &tdegB,
                        &sdegA, &sdegB,
                        &nbitsA, &nbitsB);
  fprintf(stderr, "allocation starts\n");
  data_heap_t *heap = allocate_data_heap(tdegA, tdegB,
                                       mdegA, mdegB,
                                       sdegA, sdegB,
                                       size_table);
  fprintf(stderr, "allocation done\n");

  fprintf(stderr, "Total degrees : %ld %ld\n", tdegA, tdegB);
  fprintf(stderr, "numpoints = %lu\n", heap->npts);

  //Hadamard-type (Goldstein-Graham)
  ulong bound_nbits = nbitsA*(mdegB) + nbitsB*(mdegA);

  fprintf(stderr, "Bound (nbits) = %lu\n", bound_nbits);

  mp_limb_t prime = 1073741827;
  prime = n_nextprime(1152921504606847009, 0);
  //  prime = 4611686018427388039;  //  prime = n_nextprime(prime, 0);
  ulong nbits_prime = LOG2(prime);

  ulong count_nprimes = 1;

  //max degree reached during multi-mod computation
  slong maxdegree = -1;


  double e = omp_get_wtime();
  nbits_prime = 0;

  fmpz_t coeff;
  fmpz_init(coeff);

  double e_mod = 0;

  //On teste si le resultant est square-free
  make_nmod_bpoly(heap->a_mod, A, mdegA, prime);
  make_nmod_bpoly(heap->b_mod, B, mdegB, prime);

  update_poly_tables_with_prime(heap, 0, prime);
  nmod_bpoly_resultant_allocated_feval(heap, mdegA, mdegB, 0);

  fprintf(stderr, "prime = %lu\n", prime);
  fprintf(stderr, "Degree of the resultant is %ld\n", heap->res_mod_array[0]->length - 1);
  prime = n_nextprime(prime, 0);
  double e_res = 0;
  int boo = 1;
  while(nbits_prime < bound_nbits && boo){
    int pos = 0;

    while(pos < size_table){
      make_nmod_bpoly(heap->a_mod, A, mdegA, prime);
      make_nmod_bpoly(heap->b_mod, B, mdegB, prime);

      if(heap->a_mod[mdegA]->length>0 && heap->b_mod[mdegB]->length>0){

        double start_mod = omp_get_wtime();
        update_poly_tables_with_prime(heap, pos, prime);
        double t = omp_get_wtime();
        nmod_bpoly_resultant_allocated_feval(heap, mdegA, mdegB, pos);
        e_res += (omp_get_wtime() - t);
        e_mod += (omp_get_wtime() - start_mod);
        pos++;
        fmpz_mul_ui(heap->prod_primes[pos], heap->prod_primes[pos - 1], prime);

        nbits_prime += LOG2(prime);
        prime = n_nextprime(prime, 0);
        count_nprimes++;
      }
      else{
        prime = n_nextprime(prime, 0);
      }
    }

      //on verifie si les images modulaires collent au resultant qu'on a calcule
    if(check_modular_images(res, heap, size_table, &maxdegree)){
      fprintf(stderr, "\n\nfinal prime = %lu\n\n", prime);
      boo = 0;
    }
    else{//ca colle pas, on continue la reconstruction. 
      fmpz_poly_fit_length(res, maxdegree+1);
      fmpz_one(coeff);
      for(ulong j = 0; j <= maxdegree; j++){
        fmpz_set(coeff, res->coeffs + j);
        for(int i = 0; i < size_table; i++){
          myfmpz_CRT_ui(coeff, coeff,
                      heap->prod_primes[i], (heap->res_mod_array[i]->coeffs)[j], heap->primes_array[i], 1, heap->prod_mod);
        }
        fmpz_poly_set_coeff_fmpz(res, j, coeff);
      }
      fmpz_set(heap->prod_primes[0], heap->prod_primes[size_table]);
    }
  }

  fprintf(stderr, "Wall time spent in modular resultant = %f\n", e_res);
  count_nprimes--;
  fprintf(stderr, "Number of used primes = %lu\n", count_nprimes);
  fprintf(stderr, "Elapsed time (modp): %f\n", e_mod);
  fprintf(stderr, "Elapsed time (total): %f\n", omp_get_wtime() - e);

  fprintf(stderr, "Final size coeffs: %ld\n", FLINT_ABS(_fmpz_vec_max_bits(res->coeffs, res->length)));
  fprintf(stderr, "Degree of the output: %ld\n", res->length-1);

  fmpz_clear(coeff);
  free_data_heap(mdegA, mdegB, size_table, heap);
  return;
}


static void _fmpz_bpoly_sqfree_resultant_probabilistic(fmpz_poly_t res,
                                                fmpz_poly_t *A, slong mdegA,
                                                fmpz_poly_t *B, slong mdegB){

  slong tdegA, tdegB, sdegA, sdegB;
  ulong nbitsA, nbitsB;
  get_degrees_and_nbits(A, mdegA,
                        B, mdegB,
                        &tdegA, &tdegB,
                        &sdegA, &sdegB,
                        &nbitsA, &nbitsB);
  int size_table = 4;
  data_heap_t *heap = allocate_data_heap(tdegA, tdegB,
                                       mdegA, mdegB,
                                       sdegA, sdegB,
                                       size_table);

  fprintf(stderr, "Total degrees : %ld %ld\n", tdegA, tdegB);
  fprintf(stderr, "numpoints = %lu\n", heap->npts);

  //Hadamard-type (Goldstein-Graham)
  ulong bound_nbits = nbitsA*(mdegB) + nbitsB*(mdegA);

  fprintf(stderr, "Bound (nbits) = %lu\n", bound_nbits);

  mp_limb_t prime = 1073741827;
  //On devrait mettre un peu d'alea
  prime = n_nextprime(1152921504606847009, 0);
  ulong nbits_prime = LOG2(prime);

  ulong count_nprimes = 1;

  //max degree reached during multi-mod computation
  slong maxdegree = -1;

  double e = omp_get_wtime();
  nbits_prime = 0;

  fmpz_t coeff;
  fmpz_init(coeff);

  double e_mod = 0;

  //On teste si le resultant est square-free
  make_nmod_bpoly(heap->a_mod, A, mdegA, prime);
  make_nmod_bpoly(heap->b_mod, B, mdegB, prime);

  update_poly_tables_with_prime(heap, 0, prime);
  nmod_bpoly_resultant_allocated_feval(heap, mdegA, mdegB, 0);

  fprintf(stderr, "prime = %lu\n", prime);
  fprintf(stderr, "Degree of the resultant is %ld\n", heap->res_mod_array[0]->length - 1);
  if(!is_non_square_free(heap->res_mod_array[0], heap->sqf_part)){
    fprintf(stderr, "The resultant is square-free\n");

    free_data_heap(mdegA, mdegB, size_table, heap);
    fprintf(stderr, "free heap done\n");
    _fmpz_bpoly_resultant_probabilistic(res,
                                        A, mdegA,
                                        B, mdegB);
    return;
  }
  fprintf(stderr, "The resultant is not square-free\n");
  prime = n_nextprime(prime, 0);

  while(nbits_prime < bound_nbits || 0 == 0){
    int pos = 0;

    while(pos < size_table){
      make_nmod_bpoly(heap->a_mod, A, mdegA, prime);
      make_nmod_bpoly(heap->b_mod, B, mdegB, prime);

      if(heap->a_mod[mdegA]->length>0 && heap->b_mod[mdegB]->length>0){

        double start_mod = omp_get_wtime();
        update_poly_tables_with_prime(heap, pos, prime);
        nmod_bpoly_resultant_allocated_feval(heap, mdegA, mdegB, pos);

        square_free_part_nmod_poly(heap->res_mod_array[pos], heap->sqf_part);

        e_mod += (omp_get_wtime() - start_mod);
        pos++;
        fmpz_mul_ui(heap->prod_primes[pos], heap->prod_primes[pos - 1], prime); 
        nbits_prime += LOG2(prime);
        prime = n_nextprime(prime, 0);
        count_nprimes++;
      }
      else{
        prime = n_nextprime(prime, 0);
      }
    }

    //on verifie si les images modulaires collent au resultant qu'on a calcule
    if(check_modular_images(res, heap, size_table, &maxdegree)){
      fprintf(stderr, "\n\nfinal prime = %lu\n\n", prime);
      break;
    }
    else{//ca colle pas, on continue la reconstruction. 
      fmpz_poly_fit_length(res, maxdegree+1);
      fmpz_one(coeff);
      for(ulong j = 0; j <= maxdegree; j++){
        fmpz_set(coeff, res->coeffs + j);
        for(int i = 0; i < size_table; i++){
          myfmpz_CRT_ui(coeff, coeff,
                      heap->prod_primes[i], (heap->res_mod_array[i]->coeffs)[j], heap->primes_array[i], 1, heap->prod_mod);
        }
        fmpz_poly_set_coeff_fmpz(res, j, coeff);
        //        fmpz_print(coeff);fprintf(stdout, "\n");
      }
      fmpz_set(heap->prod_primes[0], heap->prod_primes[size_table]);
    }
  }

  count_nprimes--;
  fprintf(stderr, "Number of used primes = %lu\n", count_nprimes);
  fprintf(stderr, "Elapsed time (modp): %f\n", e_mod);
  fprintf(stderr, "Elapsed time (total): %f\n", omp_get_wtime() - e);

  fprintf(stderr, "Final size coeffs: %ld\n", FLINT_ABS(_fmpz_vec_max_bits(res->coeffs, res->length)));
  fprintf(stderr, "Degree of the output: %ld\n", res->length-1);

  fmpz_clear(coeff);
  free_data_heap(mdegA, mdegB, size_table, heap);
  return;
}

/* /\* */
/* Input: A et B sont des pointeurs sur tableaux de  fmpz_poly_t de */
/* tailles respectives degA + 1 et degB + 1 */

/* *\/ */


static void fmpz_bpoly_sqfree_resultant_probabilistic(fmpz_poly_t res,
                                                      fmpz_poly_t *A, slong degA,
                                                      fmpz_poly_t *B, slong degB){
  if(degA < degB){
    _fmpz_bpoly_sqfree_resultant_probabilistic(res, B, degB, A, degA);
  }
  else{
    _fmpz_bpoly_sqfree_resultant_probabilistic(res,
                                        A, degA,
                                        B, degB);
  }
}


static void fmpz_bpoly_resultant_probabilistic(fmpz_poly_t res,
                                               fmpz_poly_t *A, slong degA,
                                               fmpz_poly_t *B, slong degB){
  if(degA < degB){
    _fmpz_bpoly_resultant_probabilistic(res, B, degB, A, degA);
  }
  else{
    _fmpz_bpoly_resultant_probabilistic(res,
                                        A, degA,
                                        B, degB);
  }
}


int main(int argc, char **argv){

  bsolve_flags * flags = allocate_initialize_flags(); 

  int binary_format = 0;
  getoptions(argc, argv, flags);

  FILE *file = fopen(flags->input_file,"r");

  if(file==NULL){
    fprintf(stderr, "No such file\n");
    return 1;
  }
  slong degA = 0, degB = 0;
  //Array est un tableau de fmpz_poly_t de lgr degA + degB + 1
  fmpz_poly_t *Array = get_input_polynomials(file, &degA, &degB);
  fmpz_poly_t *bA = Array;
  fmpz_poly_t *bB = Array  + (degA + 1);
  fclose(file);

  if(DEBUG >= 2){

    fprintf(stderr, "\nA := ");
    fmpz_bpoly_fprint_pretty(stderr, bA, degA);
    fprintf(stderr, "\n");
    fflush(stderr);
    fprintf(stderr, "B := ");
    fmpz_bpoly_fprint_pretty(stderr, bB, degB);
    fprintf(stderr, "\n");
    fflush(stderr);
  }



  if(0==1){
    flint_rand_t state;
    flint_randinit(state);
    fmpz_bpoly_random_dense(bA, degA, degA, state);
    fmpz_bpoly_random_dense(bB, degB, degB, state);
    flint_randclear(state);
  }

  fprintf(stderr, "flags->task = %d\n", flags->task);
  if(flags->task == 1){
    fmpz_poly_t res;
    fmpz_poly_init(res);
    fmpz_bpoly_resultant_probabilistic(res, bA, degA, bB, degB);
    display_resutant(res, flags, binary_format);
    fmpz_poly_clear(res);
  }
  if(flags->task == 2){
    fmpz_poly_t res;
    fmpz_poly_init(res);
    fmpz_bpoly_sqfree_resultant_probabilistic(res, bA, degA, bB, degB);
    display_resutant(res, flags, binary_format);
    fmpz_poly_clear(res);
  }
  if(DEBUG>10){
    fmpz_poly_t res;
    fmpz_poly_init(res);
    fmpz_bpoly_resultant_deterministic(res, bA, degA, bB, degB);
    fmpz_poly_clear(res);
  }
  if(flags->task==0){
    bsolve_fmpz_all_params_t *bs_param = malloc(sizeof(bsolve_fmpz_all_params_t));
    if(fmpz_bpoly_param(bs_param, bA, degA, bB, degB)){
      if(flags->output_file != NULL){
        FILE *ofile = fopen(flags->output_file, "w");
        if(flags->binary){
          print_bs_fmpz_all_params_bin(ofile, bs_param);
        }
        else{
          print_bs_fmpz_all_params_pretty(ofile, bs_param);
        }
        fclose(ofile);
      }
      else{
        if(flags->binary){
          print_bs_fmpz_all_params_bin(stdout, bs_param);
        }
        else{
          print_bs_fmpz_all_params_pretty(stdout, bs_param);
        }
      }

      free_bsolve_fmpz_all_params(bs_param);
      return 0;
    }
    free(bs_param);
    return 0;
  }

  for(ulong i = 0; i <= degA + degB + 1; i++){
    fmpz_poly_clear(Array[i]);
  }
  free(Array);
  free_flags(flags);
  return 0;
}
