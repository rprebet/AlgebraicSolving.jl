#include "data.h"

//alloue ce qu'il faut, connaissant le degree deg de la param.
static fmpz_param_t *allocate_degree_sepelem_fmpz_param(ulong deg,
                                                        mp_ptr sepelem,
                                                        int nvars){
  fmpz_param_t *param = malloc(sizeof(fmpz_param_t));

  param->recons = 0;
  param->nvars = nvars;
  param->degree = deg;
  param->sepelem = sepelem;
  param->coords = malloc(sizeof(fmpz_poly_t) * param->nvars);
  param->coords_num = malloc(sizeof(fmpz_poly_t) * param->nvars);
  param->coords_den = malloc(sizeof(fmpz_poly_t) * param->nvars);
  param->coords_lcm = malloc(sizeof(fmpz_t) * param->nvars);
  for(ulong i = 0; i < param->nvars; i++){
    fmpz_poly_init2(param->coords[i], deg);
    fmpz_poly_init2(param->coords_num[i], deg);
    fmpz_poly_init2(param->coords_den[i], deg);
    fmpz_init(param->coords_lcm[i]);
    fmpz_one(param->coords_lcm[i]);
    for(ulong j = 0; j < deg; j++){
      fmpz_poly_set_coeff_ui(param->coords[i], j, 1);
      fmpz_poly_set_coeff_ui(param->coords_num[i], j, 1);
      fmpz_poly_set_coeff_ui(param->coords_den[i], j, 1);
    }
  }

  fmpz_poly_init2(param->denom, deg + 1);
  fmpz_poly_init2(param->elim, deg + 1);

  fmpz_poly_init2(param->elim_num, deg + 1);
  fmpz_poly_init2(param->denom_num, deg + 1);

  fmpz_poly_init2(param->elim_den, deg + 1);
  fmpz_poly_init2(param->denom_den, deg + 1);

  fmpz_init(param->elim_lcm);
  fmpz_init(param->denom_lcm);

  fmpz_one(param->elim_lcm);
  fmpz_one(param->denom_lcm);

  for(ulong j = 0; j <= deg; j++){
    fmpz_poly_set_coeff_ui(param->denom, j, 1);
    fmpz_poly_set_coeff_ui(param->elim, j, 1);

    fmpz_poly_set_coeff_ui(param->denom_num, j, 1);
    fmpz_poly_set_coeff_ui(param->elim_num, j, 1);

    fmpz_poly_set_coeff_ui(param->denom_den, j, 1);
    fmpz_poly_set_coeff_ui(param->elim_den, j, 1);
  }
  param->rec_elim = 0;
  param->rec_denom = 0;
  param->rec_coords = malloc(sizeof(ulong)*nvars);
  for(int i = 0; i < nvars; i++){
    param->rec_coords[i] = 0;
  }
  return param;
}


static void allocate_degrees_sepelem_fmpz_all_params(fmpz_all_params_t * fmpz_all_params,
                                                     data_heap_t *heap,
                                                     mp_ptr sepelem,
                                                     int nb_params, int nvars){

  fmpz_all_params->nparams = nb_params;
  fmpz_all_params->params = malloc(sizeof(fmpz_param_t *) * fmpz_all_params->nparams);
  for(int i = 0; i < fmpz_all_params->nparams; i++){
    fmpz_all_params->params[i] =
      allocate_degree_sepelem_fmpz_param((heap->sqf_part->p + i)->length - 1, sepelem, nvars);
  }

}


static nmod_param_t *allocate_nmod_param(int nvars, mp_limb_t prime){
  nmod_param_t *param = malloc(sizeof(nmod_param_t));
  param->nvars = nvars;
  param->degree = 0;
  param->sepelem = NULL;
  param->coords = malloc(sizeof(nmod_poly_t) * param->nvars);
  for(int i = 0; i < param->nvars; i++){
    nmod_poly_init(param->coords[i], prime);
  }
  nmod_poly_init(param->denom, prime);
  nmod_poly_init(param->elim, prime);

  return param;
}

static nmod_param_t *allocate_degree_sepelem_nmod_param(ulong deg, mp_ptr sepelem,
                                                        int nvars, mp_limb_t prime){
  nmod_param_t *param = malloc(sizeof(nmod_param_t));
  param->nvars = nvars;
  param->degree = deg;
  param->sepelem = sepelem;
  param->coords = malloc(sizeof(nmod_poly_t) * param->nvars);
  for(int i = 0; i < param->nvars; i++){
    nmod_poly_init2(param->coords[i], prime, deg);
  }
  nmod_poly_init2(param->denom, prime, deg);
  nmod_poly_init2(param->elim, prime, deg + 1);

  return param;
}


static nmod_all_params_t *allocate_nmod_all_params(int nb_params, int nvars,
                                                   mp_limb_t prime){
  nmod_all_params_t *nmod_all_params = malloc(sizeof(nmod_all_params_t));
  nmod_all_params->nparams = nb_params;
  nmod_all_params->params = malloc(sizeof(nmod_param_t *) * nmod_all_params->nparams);
  for(int i = 0; i < nmod_all_params->nparams; i++){
    nmod_all_params->params[i] = allocate_nmod_param(nvars, prime);
  }
  return nmod_all_params;
}

static nmod_all_params_t *allocate_degrees_sepelem_nmod_all_params(data_heap_t *heap,
                                                           mp_ptr sepelem,
                                                           int nb_params, int nvars,
                                                           mp_limb_t prime){
  nmod_all_params_t *nmod_all_params = malloc(sizeof(nmod_all_params_t));
  nmod_all_params->nparams = nb_params;
  nmod_all_params->params = malloc(sizeof(nmod_param_t *) * nmod_all_params->nparams);
  for(int i = 0; i < nmod_all_params->nparams; i++){
    nmod_all_params->params[i] =
      allocate_degree_sepelem_nmod_param((heap->sqf_part->p+i)->length - 1,
                                         sepelem,
                                         nvars, prime);
  }
  return nmod_all_params;
}

static void free_bsolve_fmpz_param(bsolve_fmpz_param_t *param){
  fmpz_poly_clear(param->elim);
  fmpz_poly_clear(param->denom);
  for(int i = 0; i < param->nvars; i++){
    fmpz_poly_clear(param->coords[i]);
    fmpz_clear(param->cfs[i]);
  }
  free(param->sepelem);
  free(param->coords);
  free(param->cfs);
  free(param);
}

static void free_bsolve_fmpz_all_params(bsolve_fmpz_all_params_t *allpars){
  for(int i = 0; i < allpars->nparams; i++ ){
    free_bsolve_fmpz_param(allpars->params[i]);
  }
  free(allpars);
}

static void free_nmod_param(nmod_param_t *param){
  for(int i = 0; i < param->nvars; i++){
    nmod_poly_clear(param->coords[i]);
  }
  nmod_poly_clear(param->denom);
  nmod_poly_clear(param->elim);

  free(param);
}

static void free_nmod_all_params(nmod_all_params_t *nmod_all_params){
  for(int i = 0; i < nmod_all_params->nparams; i++){
    free_nmod_param(nmod_all_params->params[i]);
  }
  free(nmod_all_params);
}

static void free_fmpz_param(fmpz_param_t *param){
  for(int i = 0; i < param->nvars; i++){
    fmpz_poly_clear(param->coords[i]);
  }
  fmpz_poly_clear(param->denom);
  fmpz_poly_clear(param->elim);
  free(param->rec_coords);
  free(param);
}

static void free_fmpz_all_params(fmpz_all_params_t *fmpz_all_params){
  for(int i = 0; i < fmpz_all_params->nparams; i++){
    free_fmpz_param(fmpz_all_params->params[i]);
  }
  free(fmpz_all_params);
}

/*

  tdeg = total degree
  mdeg = main degree
  sdeg = degree in the other variable

*/
static data_heap_t *allocate_data_heap(slong tdegA, slong tdegB,
                               slong mdegA, slong mdegB,
                               slong sdegA, slong sdegB,
                               int size_mod_array){
  if(mdegA < mdegB){
    fprintf(stderr, "Bad call to allocator (error in degrees)\n");
    exit(1);
  }

  data_heap_t * heap_bsolve = (data_heap_t *)(malloc(sizeof(data_heap_t)));
  heap_bsolve->size_mod_array = size_mod_array;
  heap_bsolve->nmod_coeffs_array = _nmod_vec_init(3 * (mdegA + 1));

  heap_bsolve->npts = MIN(tdegA * tdegB + 1, mdegB * sdegA + mdegA * sdegB + 1);

  heap_bsolve->primes_array = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * size_mod_array));

  mp_limb_t   prime = n_nextprime(1152921504606847009, 0);
  heap_bsolve->a_mod_eval_all = malloc(sizeof(mp_ptr) * (mdegA + 1));
  for(ulong i = 0; i<=mdegA; i++){
    heap_bsolve->a_mod_eval_all[i] = malloc(sizeof(mp_limb_t) * heap_bsolve->npts);
  }

  heap_bsolve->b_mod_eval_all = malloc(sizeof(mp_ptr) * (mdegB + 1));
  for(ulong i = 0; i<=mdegB; i++){
    heap_bsolve->b_mod_eval_all[i] = malloc(sizeof(mp_limb_t) * heap_bsolve->npts);
  }


  nmod_poly_init2(heap_bsolve->a_mod_eval, prime, mdegA + 1);
  heap_bsolve->a_mod_eval->length = mdegA + 1;
  nmod_poly_init2(heap_bsolve->b_mod_eval, prime, mdegB + 1);
  heap_bsolve->b_mod_eval->length = mdegB + 1;

  heap_bsolve->res_mod_eval = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * 2 * heap_bsolve->npts));
  heap_bsolve->ev_pts = heap_bsolve->res_mod_eval + heap_bsolve->npts;

  heap_bsolve->par_h_mod_eval = (mp_limb_t *)(malloc(sizeof(mp_limb_t) * 2 * heap_bsolve->npts));
  heap_bsolve->par_t_mod_eval = heap_bsolve->par_h_mod_eval + heap_bsolve->npts;

  heap_bsolve->tree = _nmod_poly_tree_alloc(heap_bsolve->npts);

  heap_bsolve->a_mod = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * (mdegA + 1 + mdegB + 1)));
  heap_bsolve->b_mod = heap_bsolve->a_mod + mdegA + 1;

  for(ulong i = 0; i <= mdegA; i++){
    nmod_poly_init(heap_bsolve->a_mod[i], prime);
  }
  for(ulong i = 0; i <= mdegB; i++){
    nmod_poly_init(heap_bsolve->b_mod[i], prime);
  }
  heap_bsolve->res_mod_array = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * ( 1 + size_mod_array )));
  for(ulong i = 0; i <= size_mod_array; i++){
    nmod_poly_init(heap_bsolve->res_mod_array[i], prime);
  }
  heap_bsolve->res_mod = heap_bsolve->res_mod_array + size_mod_array;

  heap_bsolve->res_mod_prime_array = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * ( size_mod_array )));
  for(ulong i = 0; i < size_mod_array; i++){
    nmod_poly_init(heap_bsolve->res_mod_prime_array[i], prime);
  }

  heap_bsolve->par_h_mod_array = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * (  size_mod_array )));
  for(ulong i = 0; i < size_mod_array; i++){
    nmod_poly_init(heap_bsolve->par_h_mod_array[i], prime);
  }

  heap_bsolve->par_t_mod_array = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * ( size_mod_array )));
  for(ulong i = 0; i < size_mod_array; i++){
    nmod_poly_init(heap_bsolve->par_t_mod_array[i], prime);
  }

  heap_bsolve->GCDINV = (nmod_poly_t *)(malloc(sizeof(nmod_poly_t) * ( size_mod_array )));
  for(ulong i = 0; i < size_mod_array; i++){
    nmod_poly_init(heap_bsolve->GCDINV[i], prime);
  }

  heap_bsolve->prod_primes = (fmpz_t *)(malloc(sizeof(fmpz_t)*(size_mod_array + 1)));

  for(int i  = 0; i <= size_mod_array; i++){
    fmpz_init(heap_bsolve->prod_primes[i]);
    fmpz_one(heap_bsolve->prod_primes[i]);
  }

  nmod_poly_factor_init(heap_bsolve->sqf_part);

  heap_bsolve->nmod_lparams_array = malloc(sizeof(nmod_all_params_t *) * size_mod_array);
  for(int i = 0; i < size_mod_array; i++){
    heap_bsolve->nmod_lparams_array[i] = allocate_nmod_all_params(1, 2, prime);
  }

  heap_bsolve->coef = malloc(sizeof(mpz_rat_t));
  fmpz_init(heap_bsolve->coef->numer);
  fmpz_init(heap_bsolve->coef->denom);

  fmpz_init(heap_bsolve->h_prod);

  fmpz_init(heap_bsolve->lcm);
  fmpz_one(heap_bsolve->lcm);

  heap_bsolve->prod_mod = malloc(sizeof(fmpz_t));
  fmpz_init(heap_bsolve->prod_mod[0]);
  //allocation si besoin dans la fonction de calcul de parametrisation
  //square-free (quand le resultant ne l'est pas)
  return heap_bsolve;
}


static void free_data_heap(slong mdegA, slong mdegB,
                           int size_mod_array, data_heap_t *heap_bsolve){

  _nmod_vec_clear(heap_bsolve->nmod_coeffs_array);
  nmod_poly_factor_clear(heap_bsolve->sqf_part);
  _nmod_poly_tree_free(heap_bsolve->tree, heap_bsolve->npts);

  free(heap_bsolve->primes_array);

  for(ulong i = 0; i<=mdegA; i++){
    free(heap_bsolve->a_mod_eval_all[i]);
  }
  free(heap_bsolve->a_mod_eval_all);

  for(ulong i = 0; i<=mdegB; i++){
    free(heap_bsolve->b_mod_eval_all[i]);
  }
  free(heap_bsolve->b_mod_eval_all);

  nmod_poly_clear(heap_bsolve->a_mod_eval);
  nmod_poly_clear(heap_bsolve->b_mod_eval);

  free(heap_bsolve->res_mod_eval);
  free(heap_bsolve->par_h_mod_eval);

  for(ulong i = 0; i <= mdegA; i++){
    nmod_poly_clear(heap_bsolve->a_mod[i]);
  }
  for(ulong i = 0; i <= mdegB; i++){
    nmod_poly_clear(heap_bsolve->b_mod[i]);
  }
  free(heap_bsolve->a_mod);

  for(ulong i = 0; i <= size_mod_array; i++){
    nmod_poly_clear(heap_bsolve->res_mod_array[i]);
  }
  free(heap_bsolve->res_mod_array);

  for(ulong i = 0; i < size_mod_array; i++){
    nmod_poly_clear(heap_bsolve->res_mod_prime_array[i]);
  }
  free(heap_bsolve->res_mod_prime_array);

  for(ulong i = 0; i < size_mod_array; i++){
    nmod_poly_clear(heap_bsolve->par_h_mod_array[i]);
  }
  free(heap_bsolve->par_h_mod_array);

  for(ulong i = 0; i < size_mod_array; i++){
    nmod_poly_clear(heap_bsolve->par_t_mod_array[i]);
  }
  free(heap_bsolve->par_t_mod_array);

  for(ulong i = 0; i < size_mod_array; i++){
    nmod_poly_clear(heap_bsolve->GCDINV[i]);
  }
  free(heap_bsolve->GCDINV);

  for(int i  = 0; i <= size_mod_array; i++){
    fmpz_clear(heap_bsolve->prod_primes[i]);
  }

  fmpz_clear(heap_bsolve->coef->numer);
  fmpz_clear(heap_bsolve->coef->denom);
  free(heap_bsolve->coef);

  fmpz_clear(heap_bsolve->h_prod);
  fmpz_clear(heap_bsolve->lcm);
  fmpz_clear(heap_bsolve->prod_mod[0]);
  free(heap_bsolve->prod_mod);

  //attention : ecrire une fonction speciale pour desallouer ces zones car elles
  //ne sont allouees que quand c'est necessaire (resultant non square-free et
  //calcul de parametrisation)

  /* for(int i = 0; i < size_mod_array; i++){ */
  /*   free_nmod_all_params(heap_bsolve->nmod_lparams_array[i]); */
  /* } */
  free(heap_bsolve->nmod_lparams_array);
}
