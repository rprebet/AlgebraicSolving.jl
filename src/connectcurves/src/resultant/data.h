typedef struct{
  fmpz_t numer;
  fmpz_t denom;
} mpz_rat_t;


typedef struct{
  int nvars;
  mp_ptr sepelem;
  fmpz_poly_t * coords;
  fmpz_poly_t denom;
  fmpz_poly_t elim;
  fmpz_t *cfs;
} bsolve_fmpz_param_t;


typedef struct{
  int nparams;
  bsolve_fmpz_param_t **params; 
} bsolve_fmpz_all_params_t; 

typedef struct{
  int recons;
  int nvars;
  int degree;
  mp_ptr sepelem;
  fmpz_poly_t * coords; // numerateurs des parametrisations
  fmpz_poly_t denom; //denominateur
  fmpz_poly_t elim; //polynome eliminant
  fmpz_poly_t * coords_num; // numerateurs des parametrisations
  fmpz_poly_t denom_num; //denominateur
  fmpz_poly_t elim_num; //polynome eliminant
  fmpz_poly_t * coords_den; // numerateurs des parametrisations
  fmpz_poly_t denom_den; //denominateur
  fmpz_poly_t elim_den; //polynome eliminant
  fmpz_t * coords_lcm; // numerateurs des parametrisations
  fmpz_t denom_lcm; //denominateur
  fmpz_t elim_lcm; //polynome eliminant
  ulong *rec_coords;
  ulong rec_denom;
  ulong rec_elim;
} fmpz_param_t;


typedef struct{
  int nvars; //nvars est le nbre de variables initiales. 
  ulong degree; // degre de la parametrisation
  mp_ptr sepelem;
  nmod_poly_t * coords; // numerateurs des parametrisations
  //on en a autant que de coordonnees initiales
  nmod_poly_t denom; // polynome denomninateur
  nmod_poly_t elim; // polynome eliminant
} nmod_param_t;
// donc au total on a exactement nvars + 1 polynomes


typedef struct{
  int nparams; // nombre de parametrisations
  fmpz_param_t **params; // tableau de parametrisations
} fmpz_all_params_t;


typedef struct{
  int nparams;
  nmod_param_t **params;
} nmod_all_params_t;


typedef struct{
  int verb ;
  int task ; // 0 -> parametrisation
    // 1 -> square-free resultant
    // 2 -> resultant
  int binary ; // 1 -> output file in binary format
  char *input_file ;
  char *output_file ;
} bsolve_flags;

/*

  Zone memoire pour donnees intermediaires.

  Le but est de calculer le resultant de A et B qui sont des polynomes bivaries.
  (ou bien de calculer une parametrisation rationnelle des solutions)

*/
typedef struct
{
  ulong size_mod_array;
  mp_ptr nmod_coeffs_array; // tableau pour suite des restes euclidiens
  ulong npts; // nbre de points pour l'evaluation interpolation

  mp_ptr primes_array; // tableau de premiers utilises pour les calculs
                       // multi-modulaires
  mp_ptr ev_pts; // tableau contenant les points d'evaluation utilises pour
                 // l'interpolation
  mp_ptr *tree;  // arbre pour l'evaluation rapide multi-point
  mp_ptr *a_mod_eval_all; // tableau qui recevra les coeffs de A mod prime
                               // dont on a specialise une variable aux entrees
                               // de ev_pts
  mp_ptr *b_mod_eval_all; // tableau qui recevra les coeffs de B mod prime
                               // dont on a specialise une variable aux entrees
                               // de ev_pts
  nmod_poly_t a_mod_eval; // evaluation de a_mod
  nmod_poly_t b_mod_eval; // evaluation de b_mod

  nmod_poly_t * a_mod; // image modulare de A
  nmod_poly_t * b_mod; // image modulaire de B


  mp_ptr res_mod_eval; // tableau recevant les specialisations du resultant de 
                       // A mod prime et B mod prime

  nmod_poly_t * res_mod_array; // tableau contenant les images modulaires du resultant
  nmod_poly_t * res_mod_prime_array; // tableau contenant les images modulaires
                                     // de la derivee du resultant
  nmod_poly_t * res_mod; // tableau contenant une image modulaire du resultant ;
                         // sera utilise pour verifier qu'on a fini le calcul

  // la longueur de ces tableaux est la meme que celle de res_mod_eval
  mp_ptr par_h_mod_eval; // tableau recevant les coeffs dominants de
                         // sous-resultants (univaries) de degre 1
                         // ici on a evalue une variable
  mp_ptr par_t_mod_eval; // tableau recevant les coeffs constants de
                         // sous-resultants (univaries) de degre 1
                         // ici on a evalue une variable

  nmod_poly_t * par_h_mod_array; // tableau contenant les images modulaires du
                                 // coeff dominant du sous-resultant de degre 1
  nmod_poly_t * par_t_mod_array; // tableau contenant les images modulaires du
                                 // coeff constant du sous-resultant de degre 1

  nmod_poly_t * GCDINV; // quand le GCD du coeff de dom. du
  // sous-resultant de degre 1 avec le resultant est = 1 on inverse le coeff
  // dominant modulo le resultant

  nmod_all_params_t **nmod_lparams_array; //tableau contenant des images
                                          //modulaires des parametrisations

  fmpz_t *prod_primes; // contiendra des produits de nbre premiers
  nmod_poly_factor_t sqf_part; // contiendra la partie square-free

  //data pour les reconstructions.
  fmpz_t *prod_mod;
  mpz_rat_t *coef;
  fmpz_t h_prod;
  fmpz_t lcm;
  fmpz_t *coef_num_array;
  fmpz_t *coef_den_array;
} data_heap_t;
