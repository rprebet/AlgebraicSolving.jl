static void print_nmod_param(FILE * file, nmod_param_t *param){
  char str[2]={"x"};
  fprintf(file, "%d\n", param->nvars);
  fprintf(file, "[");

  for(int i = 0; i < param->nvars; i++){
    nmod_poly_fprint_pretty(file, param->coords[i], str);
  }
  fprintf(file, "]\n");
  nmod_poly_fprint_pretty(file, param->denom, str);
  fprintf(file, "\n");
  nmod_poly_fprint_pretty(file, param->elim, str);
  fprintf(stderr, "\n");
}


static inline void print_nmod_all_params_pretty(FILE *file, nmod_all_params_t *all_params){
  fprintf(file, "\n");
  for(int i = 0; i < all_params->nparams; i++){
    print_nmod_param(file, (all_params->params)[i]);
  }
}

static inline void fmpz_poly_fprint_bin(FILE *file, const fmpz_poly_t poly){
  slong len = poly->length;
  fprintf(file, "%ld\n", len);
  fprintf(stderr, "LENGTH = %ld\n", len);
  for(slong i = 0; i < len; i++){
    fmpz_out_raw(file, poly->coeffs + i);
  }
}

static inline void print_fmpz_param_bin(FILE *file, fmpz_param_t *param){
  fprintf(file, "%d\n", param->nvars);
  for(int i = 0; i < param->nvars-1; i++){
    fmpz_poly_fprint_bin(file, param->coords[i]);
  }
  fmpz_poly_fprint_bin(file, param->denom);
  fmpz_poly_fprint_bin(file, param->elim);
}

static inline void print_fmpz_param_pretty(FILE *file, fmpz_param_t *param){
  char stra[2] = {"x"};
  fprintf(file, "%d\n", param->nvars);
  fprintf(file, "[");

  for(int i = 0; i < param->nvars; i++){
    fmpz_poly_fprint_pretty(file, param->coords[i], stra);
    fprintf(file, ",\n");
  }
  fmpz_poly_fprint_pretty(file, param->coords[param->nvars-1], stra);
  fprintf(file, "]\n");
  fmpz_poly_fprint_pretty(file, param->denom, stra);
  fprintf(file, "\n");
  fmpz_poly_fprint_pretty(file, param->elim, stra);
  fprintf(stderr, "\n");
}

static inline void print_bs_fmpz_param_bin(FILE *file, bsolve_fmpz_param_t *param){
  fprintf(stderr, "ELIM -> ");
  fmpz_poly_fprint_bin(file, param->elim);
  fprintf(stderr, "DENOM -> ");
  fmpz_poly_fprint_bin(file, param->denom);
  int nv = param->nvars + 1;
  fprintf(stderr, "NVARS = %d\n", nv);
  fprintf(file, "%d\n", nv);
  for(int i = 0; i < param->nvars; i++){
    fprintf(stderr, "COORD -> ");
    fmpz_poly_fprint_bin(file, param->coords[i]);
    fmpz_out_raw(file, param->cfs[i]);
  }
}

static inline void print_bs_fmpz_param_pretty(FILE *file, bsolve_fmpz_param_t *param){
  char stra[2] = {"x"};
  fprintf(file, "[");
  fmpz_poly_fprint_pretty(file, param->elim, stra);
  fprintf(file, ",\n");
  fmpz_poly_fprint_pretty(file, param->denom, stra);
  fprintf(file, ",\n");
  for(int i = 0; i < param->nvars-1; i++){
    fmpz_poly_fprint_pretty(file, param->coords[i], stra);
    fprintf(file, ",\n");
  }
  fmpz_poly_fprint_pretty(file, param->coords[param->nvars-1], stra);
  fprintf(file, "]");
}

static inline void print_fmpz_all_params_bin(FILE *file, fmpz_all_params_t *all_params){
  fprintf(file, "%d\n", all_params->nparams);

  for(int i = 0; i < all_params->nparams; i++){
    print_fmpz_param_bin(file, all_params->params[i]);
  }
}


static inline void print_fmpz_all_params_pretty(FILE *file, fmpz_all_params_t *all_params){
  fprintf(file, "\n");
  for(int i = 0; i < all_params->nparams; i++){
    print_fmpz_param_pretty(file, all_params->params[i]);
  }
}

static inline void print_bs_fmpz_all_params_bin(FILE *file,
                                                bsolve_fmpz_all_params_t *all_params){
  fprintf(file, "%d\n", all_params->nparams);

  for(int i = 0; i < all_params->nparams-1; i++){
    fprintf(stderr, "Display PARAM %d\n", i);
    print_bs_fmpz_param_bin(file, all_params->params[i]);
    fprintf(file, "\n");
  }
  fprintf(stderr, "Display PARAM %d\n", all_params->nparams-1);
  print_bs_fmpz_param_bin(file, all_params->params[all_params->nparams-1]);
  fprintf(file, "\n");
}

static inline void print_bs_fmpz_all_params_pretty(FILE *file,
                                                   bsolve_fmpz_all_params_t *all_params){
  fprintf(file, "[\n");
  for(int i = 0; i < all_params->nparams-1; i++){
    print_bs_fmpz_param_pretty(file, all_params->params[i]);
    fprintf(file, ",\n");
  }
  print_bs_fmpz_param_pretty(file, all_params->params[all_params->nparams-1]);
  fprintf(file, "\n];\n");
}

static void nmod_vec_print(FILE *file, mp_srcptr vec, slong len){
  if(len==1) {
    fprintf(file, "%lu ; \n", vec[0]);
    return;
  }
  fprintf(file, "%lu + ", vec[0]);
  for(slong i = 1; i < len - 1; i++){
    fprintf(file, "%lu * x^%ld + ", vec[i], i);
  }
  fprintf(file, "%lu * x^%ld; \n", vec[len-1], len-1);
}

/*
  Input:
*/

static void fmpz_bpoly_fprint_pretty(FILE *file, fmpz_poly_t *A, ulong deg){
  char strx[1]={"x"};
  for(slong i=deg; i >=1 ; i--){
    fprintf(file, "(");
    fmpz_poly_fprint_pretty(file, A[i], strx);
    fprintf(file, ")*y^%lu+", i);
  }
  fprintf(file, "(");
  fmpz_poly_fprint_pretty(file, A[0], strx);
  fprintf(file, ");");
  fprintf(file, "\n");
}

static ulong nbits_bound_fmpz_bvar_poly(fmpz_poly_t *A, ulong deg){
  ulong nbits = 0;
  fmpz_t sum, sumtmp, tmp;
  fmpz_init(sum);
  fmpz_init(sumtmp);
  fmpz_init(tmp);
  fmpz_zero(sum);
  for(ulong i = 0; i <= deg; i++){
    fmpz_zero(sumtmp);
    for(ulong j = 0; j < A[i]->length; j++){
      fmpz_abs(tmp, A[i]->coeffs+j); 
      fmpz_add(sumtmp, sumtmp, tmp);
    }
    fmpz_mul(sumtmp, sumtmp, sumtmp);
    fmpz_add(sum, sum, sumtmp);
  }
  nbits = fmpz_bits(sum) + 1;
  fmpz_clear(sumtmp);
  fmpz_clear(sum);
  fmpz_clear(tmp);
  return nbits / 2;
}

static void crt_lift_resultant(fmpz_poly_t res, nmod_poly_t *res_mod_table, mp_limb_t *primes_table, ulong size, fmpz_t *prod_mod){
  slong maxdegree = -1;
  ulong i;
  for(i = 0; i < size; i++){
    if(res_mod_table[i]->length - 1 > maxdegree){
      maxdegree = res_mod_table[i]->length - 1;
    }
  }
  if(res->length-1!=maxdegree){
    res->coeffs = malloc(sizeof(fmpz_t) * (maxdegree + 1));
    for(i = 0; i <= maxdegree; i++){
      fmpz_poly_set_coeff_ui(res, i, 1);
    }
  }

  fmpz_t prod_primes;
  fmpz_init(prod_primes);
  fmpz_one(prod_primes);

  for(i = 0; i < size; i++){
    for(ulong j = 0; j <= maxdegree; j++){
      myfmpz_CRT_ui((res->coeffs)+j, (res->coeffs)+j,
                      prod_primes, (res_mod_table[i]->coeffs)[j],
                    primes_table[i], 1,
                    prod_mod);
    }
    fmpz_mul_ui(prod_primes, prod_primes, primes_table[i]);
  }
  fmpz_clear(prod_primes);
}



static slong get_total_degree(fmpz_poly_t *A, slong degA){
  slong total_degree = -1;

  for(slong i = 0; i <= degA; i++){
    if(A[i]->length - 1 + i > total_degree){
      total_degree = A[i]->length - 1 + i; 
    }
  }
  return total_degree;
}

static slong get_partial_degree(fmpz_poly_t *A, slong degA){
  slong degAx = -1;
  for(slong i = 0; i <= degA ; i++){
    if(A[i]->length - 1 > degAx){
      degAx = A[i]->length -1;
    }
  }
  return degAx;
}


static void get_exponent(FILE *file, ulong *exp, ulong nvars){
  for(ulong i = 0; i < nvars-1; i++){
    if(fscanf(file, "%lu ", exp+i)){
    }
  }
  if(fscanf(file, "%lu\n", exp+nvars-1)){
  }
}

static void get_mpoly_from_file(FILE *file, fmpz_mpoly_t A, fmpz_mpoly_ctx_t ctx){
  ulong nterms, nvars;
  if(fscanf(file, "%lu", &nterms)){
  }
  else{
    fprintf(stderr, "Error when reading file\n");
  }
  if(fscanf(file, "%lu", &nvars)){
  }
  else{
    fprintf(stderr, "Error when reading file\n");
  }
  if(nvars!=2){
    fprintf(stderr, "Bad number of variables\n");
    exit(1);
  }
  fmpz_mpoly_fit_length(A, nterms, ctx);
  fmpz_t coeff;
  fmpz_init(coeff);
  ulong * exp = malloc(sizeof(ulong)*nvars);
  for(ulong i=0; i < nterms; i++){
    fmpz_fread(file, coeff);
    get_exponent(file, exp, nvars);
    //    _fmpz_mpoly_emplacebackterm_fmpz_ui(A, coeff, exp, ctx);
    fmpz_mpoly_push_term_fmpz_ui(A, coeff, exp, ctx);
  }
  fmpz_mpoly_sort_terms(A, ctx);
  fmpz_mpoly_combine_like_terms(A, ctx);
  free(exp);
  fmpz_clear(coeff);
}


/*

Input: A et B dependent de deux variables (x1,x2)
L'ordre utilise est l'ordre lex x1 > x2

On renvoie un pointeur sur les coefficients de A suivis de
ceux de B vus comme un polynome de Z[x2][x1]

main_degree_A et main_degree_B sont les degres en x1 de A et B

*/
static fmpz_poly_t *make_bpoly_from_mpoly2(fmpz_mpoly_t A, fmpz_mpoly_t B, fmpz_mpoly_ctx_t ctx,
                                    slong *main_degree_A, slong *main_degree_B){

  *main_degree_A = fmpz_mpoly_degree_si(A, 0, ctx);
  *main_degree_B = fmpz_mpoly_degree_si(B, 0, ctx);

  fmpz_poly_t *A_univ = (fmpz_poly_t *)(malloc(sizeof(fmpz_poly_t) *
                                               (*main_degree_A + 1 + *main_degree_B + 1)));
  for(ulong i = 0; i <= (*main_degree_A) + (*main_degree_B) + 1; i++){
    fmpz_poly_init(A_univ[i]);
    fmpz_poly_zero(A_univ[i]);
  }
  fmpz_mpoly_t Term;
  fmpz_mpoly_init(Term, ctx);
  slong deg_mainvar_term = 0, old_deg_mainvar_term = -1;
  slong deg_var2;
  fmpz_t c;
  fmpz_init(c);
  for(ulong i = 0; i < A->length; i++){
    fmpz_mpoly_get_term(Term, A, i, ctx);
    fmpz_mpoly_get_term_coeff_fmpz(c, Term, 0, ctx);

    deg_mainvar_term = fmpz_mpoly_degree_si(Term, 0, ctx);
    deg_var2 = fmpz_mpoly_degree_si(Term, 1, ctx);
    if(old_deg_mainvar_term != deg_mainvar_term){
      fmpz_poly_fit_length(A_univ[deg_mainvar_term], deg_var2 + 1);
    }
    fmpz_poly_set_coeff_fmpz(A_univ[deg_mainvar_term], deg_var2, c);
    old_deg_mainvar_term = deg_mainvar_term;
  }
  for(ulong i = 0; i < B->length; i++){
    fmpz_mpoly_get_term(Term, B, i, ctx);
    fmpz_mpoly_get_term_coeff_fmpz(c, Term, 0, ctx);

    deg_mainvar_term = fmpz_mpoly_degree_si(Term, 0, ctx);
    deg_var2 = fmpz_mpoly_degree_si(Term, 1, ctx);
    if(old_deg_mainvar_term!=deg_mainvar_term){
      fmpz_poly_fit_length(A_univ[*main_degree_A + 1 + deg_mainvar_term], deg_var2 + 1);
    }
    fmpz_poly_set_coeff_fmpz(A_univ[*main_degree_A + 1 + deg_mainvar_term], deg_var2, c);
    old_deg_mainvar_term = deg_mainvar_term;
  }

  fmpz_mpoly_clear(Term, ctx);
  fmpz_clear(c);
  return A_univ;
}

#ifdef TESTDENSE
static void test_dense(slong deg, ulong nbits, flint_rand_t state){
  fmpz_poly_t *A = (fmpz_poly_t *)(malloc(sizeof(fmpz_poly_t) * (deg + 1)));
  fmpz_poly_t *B = (fmpz_poly_t *)(malloc(sizeof(fmpz_poly_t) * (deg + 1)));
  for(ulong i = 0; i <= deg; i++){
    fmpz_poly_init(A[i]);
    fmpz_poly_init(B[i]);
  }

  fmpz_bpoly_random_dense(A, deg, nbits, state);
  fmpz_bpoly_random_dense(B, deg, nbits, state);

  fmpz_poly_t res;
  fmpz_poly_init(res);
  fmpz_poly_t res2;
  fmpz_poly_init(res2);
  fmpz_poly_t res3;
  fmpz_poly_init(res3);

  fmpz_bpoly_resultant_deterministic(res, A, deg, B, deg);
  fprintf(stderr, "\n\nStarts probabilistic implem.\n\n");
  fmpz_bpoly_resultant_probabilistic(res2, A, deg, B, deg, 0);
  //  fmpz_bpoly_sqfree_resultant_probabilistic(res2, A, deg, B, deg, 0);
  fprintf(stderr, "\n\nStarts SQRFREE probabilistic implem.\n\n");
  fmpz_bpoly_sqfree_resultant_probabilistic(res2, A, deg, B, deg);


  /* fprintf(stdout,"R<x>:=PolynomialRing(Rationals());\n"); */
  /* fprintf(stdout,"P<y>:=PolynomialRing(R);\n"); */
  /* fprintf(stdout, "A := "); */
  /* fmpz_bpoly_fprint_pretty(A, deg); */
  /* fprintf(stdout, "\nB := "); */
  /* fmpz_bpoly_fprint_pretty(B, deg); */
  /* char str[1]={"x"}; */
  /* fprintf(stdout, "\nmyres := "); */
  /* fmpz_poly_fprint_pretty(stdout, res2, str); */
  /* fprintf(stdout, ";\n"); */
  /* fprintf(stdout, "\ntime res := Resultant(A, B);\n"); */
  /* fprintf(stdout, "res-myres;"); */
  /* fprintf(stdout, "\n"); */

  for(ulong i = 0; i <= deg; i++){
    fmpz_poly_clear(A[i]);
    fmpz_poly_clear(B[i]);
  }
  fmpz_poly_clear(res);
  fmpz_poly_clear(res2);
  fmpz_poly_clear(res3);
  free(A);
  free(B);
}
#endif

static void display_help(char **argv){
  fprintf(stdout, "Basic usage\n\n%s -f file_in\n\n", argv[0]);
  fprintf(stdout, "Whenever it is possible, by default, it returns\n");
  fprintf(stdout, "a finite sequence of zero-dimensional parametrizations\n");
  fprintf(stdout, "of the set of common complex solutions to the input polynomials\n\n");

  fprintf(stdout, "Each parametrization corresponds to distinct multiplicities \n");
  fprintf(stdout, "of the resultant of the input polynomials w.r.t. the first variable\n\n");

  fprintf(stdout, "Additional options: \n");
  fprintf(stdout, "-r \t Computes the resultant\n");
  fprintf(stdout, "-q \t Computes the square-free part of the resultant\n");
  fprintf(stdout, "-b \t Output file will be in binary format\n");
  fprintf(stdout, "-o \t output file \n");
  fprintf(stdout, "-h \t displays this help\n\n");

  fprintf(stdout, "\n\n file_in must have the following format: \n");
  fprintf(stdout, "First line looks like\nnt n\nwhere nt is the number of terms and n the number of variables\n");
  fprintf(stdout, "Next lines look like \ncoeff e1 e2\n where coeff is an integer coefficient to the monomial x1^e1...xn^en\n");
  fprintf(stdout, "And a next polynomial should be given with the same format\n");
  fprintf(stdout, "\nThe variable which is eliminated is the first one (the one whose exponent comes at first)\n");
}


static fmpz_poly_t *get_input_polynomials(FILE *file, slong *degA, slong *degB){
  fmpz_mpoly_ctx_t ctx;
  ordering_t ord = ORD_LEX;
  fmpz_mpoly_t A, B;
  fmpz_mpoly_ctx_init(ctx, 2, ord);
  fmpz_mpoly_init(A, ctx);
  fmpz_mpoly_init(B, ctx);

  get_mpoly_from_file(file, A, ctx);
  get_mpoly_from_file(file, B, ctx);

  fmpz_poly_t *bA = make_bpoly_from_mpoly2(A, B, ctx, degA, degB);

  fmpz_mpoly_clear(A, ctx);
  fmpz_mpoly_clear(B, ctx);
  fmpz_mpoly_ctx_clear(ctx);
  return bA;
}

static void display_fmpz_poly_to_usolve(FILE *file, fmpz_poly_t res){
  fprintf(file, "%ld\n", res->length-1);
  fmpz_t c;
  fmpz_init(c);
  if(res->length==0){
    fprintf(file,"%d\n", 0);
    return;
  }
  for(ulong i=0; i < res->length; i++){
    fmpz_poly_get_coeff_fmpz(c, res, i);
    fmpz_fprint(file, c);
    fprintf(file, "\n");
  }
  fmpz_clear(c);
}

static void display_fmpz_poly_binary_to_usolve(FILE *file, fmpz_poly_t res){
  return;
  /* const long length = res->length; */
  /* if(res->length==0){ */
  /*   mpz_t * tab = (mpz_t *)(malloc(sizeof(mpz_t) )); */
  /*   mpz_init(tab[0]); */
  /*   mpz_set_ui(tab[0], 0); */
  /*   mpz_t p; */
  /*   binary_write_INT(file, p);//cette fonction fait un mpz_clear(p) */
  /*   binary_write(file, tab, length); */
  /*   mpz_clear(tab[0]); */
  /*   free(tab); */
  /*   return; */
  /* } */
  /* mpz_t * tab = (mpz_t *)(malloc(sizeof(mpz_t) * length)); */
  /* for(long i=0; i < length; i++){ */
  /*   mpz_init(tab[i]); */
  /*   fmpz_get_mpz(tab[i], res->coeffs+i); */
  /* } */
  /* mpz_t p; */
  /* binary_write_INT(file, p);//cette fonction fait un mpz_clear(p) */
  /* binary_write(file, tab, length); */
  /* for(long i=0; i < length; i++){ */
  /*   mpz_clear(tab[i]); */
  /* } */
  /* free(tab); */
}


static char *getoptions(int argc, char **argv, bsolve_flags * flags){
  int optch, errflag = 0, fflag = 1; 
	extern int opterr;

	char options[] = "rqhbf:o:v:";
  char *filename = NULL;
	opterr = 1;

	while ((optch = getopt(argc, argv, options)) != -1)
    switch (optch) {
    case 'b':
      flags->binary = 1;
      break;
		case 'f':
      fflag = 0;
      flags->input_file = optarg;
      optind--;
			break;
    case 'h':
      display_help(argv);
      exit(1);
    case 'o':
      flags->output_file = optarg;
      optind--;
      break;
    case 'q':
      flags->task = 2;
      break;
    case 'r':
      flags->task = 1;
      break;
    case 'v':
      flags->verb = strtol(optarg, NULL, 10);
      optind--;
      break;
    default:
      errflag++;
      break;
    }
  if(fflag){
    fprintf(stderr,"No given file\n");
    fprintf(stderr,"%s -f file_name\n", argv[0]);
    exit(1);
  }
  if(errflag){
    fprintf(stderr,"Invalid usage\n");
    fprintf(stderr,"%s -f file_name\n", argv[0]);
    exit(1);
  }
  return filename;
}

static bsolve_flags *allocate_initialize_flags(){
  bsolve_flags *flags = malloc(sizeof(bsolve_flags));
  flags->verb=0;
  flags->task = 0;
  flags->binary = 0;
  flags->input_file = NULL;
  flags->output_file = NULL;

  return flags;
}

static void free_flags(bsolve_flags * flags){
  free(flags);
}

static inline void display_resutant(fmpz_poly_t res, bsolve_flags *flags,
                                    int binary_format){
  if(flags->output_file!=NULL){
    FILE *outputfile = fopen(flags->output_file,"w");
    if(binary_format==0){
      display_fmpz_poly_to_usolve(outputfile, res);
    }
    else{
      display_fmpz_poly_binary_to_usolve(outputfile, res);
    }
    fclose(outputfile);
  }
  else{
    display_fmpz_poly_to_usolve(stdout, res);
  }
}

