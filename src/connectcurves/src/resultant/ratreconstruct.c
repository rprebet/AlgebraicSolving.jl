#define ROT(u,v,t)   \
    do { fmpz _t = *u; *u = *v; *v = *t; *t = _t; } while (0);

void myfmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
                   ulong r2, ulong m2, int sign, fmpz_t *m1m2)
{
  mp_limb_t c;
  //  fmpz_t m1m2;

  c = fmpz_fdiv_ui(m1, m2);
  c = n_invmod(c, m2);

  if (c == 0)
    {
      flint_printf("Exception (fmpz_CRT_ui). m1 not invertible modulo m2.\n");
      flint_abort();
    }

  //  fmpz_init(m1m2);
  fmpz_mul_ui(m1m2[0], m1, m2);

  _fmpz_CRT_ui_precomp(out, r1, m1, r2, m2, n_preinvert_limb(m2),
                       m1m2[0], c, sign);

  //  fmpz_clear(m1m2);
}

int
_myfmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d,
    const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
{
    fmpz_t q, r, s, t;
    int success = 0;

    /* Quickly identify small integers */
    if (fmpz_cmp(a, N) <= 0)
    {
        fmpz_set(n, a);
        fmpz_one(d);
        return 1;
    }
    fmpz_sub(n, a, m);
    if(fmpz_cmp_ui(n, 0)>=0){
      if(fmpz_cmp(n, N) <= 0){
        fmpz_one(d);
        return 1;
      }
    }
    else{
      fmpz_neg(n, n);
      if(fmpz_cmp(n, N) <= 0){
        fmpz_one(d);
        fmpz_neg(n, n);
        return 1;
      }
    }
    /* if (fmpz_cmpabs(n, N) <= 0) */
    /* { */
    /*     fmpz_one(d); */
    /*     fprintf(stderr, "ici?\n"); */
    /*     return 1; */
    /* } */

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(s);
    fmpz_init(t);

    fmpz_set(r, m); fmpz_zero(s);
    fmpz_set(n, a); fmpz_one(d);

    while (fmpz_cmpabs(n, N) > 0)
    {
        fmpz_fdiv_q(q, r, n);
        fmpz_mul(t, q, n); fmpz_sub(t, r, t); ROT(r, n, t);
        fmpz_mul(t, q, d); fmpz_sub(t, s, t); ROT(s, d, t);
    }

    if (fmpz_sgn(d) < 0)
    {
        fmpz_neg(n, n);
        fmpz_neg(d, d);
    }

    if (fmpz_cmp(d, D) <= 0)
    {
        fmpz_gcd(t, n, d);
        success = fmpz_is_one(t);
    }

    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(s);
    fmpz_clear(t);
    return success;
}

int
myfmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m,
                        const fmpz_t N, const fmpz_t D)
{
  return _myfmpq_reconstruct_fmpz_2(fmpq_numref(res),
                                  fmpq_denref(res), a, m, N, D);
}

int
_myfmpq_reconstruct_fmpz(fmpz_t n, fmpz_t d,
                       const fmpz_t a, const fmpz_t m)
{
  fmpz_t N;
  int result;

  fmpz_init(N);
  fmpz_fdiv_q_2exp(N, m, 1);
  fmpz_sqrt(N, N);
  result = _myfmpq_reconstruct_fmpz_2(n, d, a, m, N, N);
  fmpz_clear(N);

  return result;
}

int
myfmpq_reconstruct_fmpz(mpz_rat_t *res, fmpz_t a, const fmpz_t m)
{
  if(fmpz_cmp_ui(a, 0)>=0){
    return _myfmpq_reconstruct_fmpz(res->numer,
                                    res->denom, a, m);
  }
  else{
    while(fmpz_cmp_ui(a, 0) < 0){
      //      fmpz_fprint(stderr, a); fprintf(stderr, "\n");
      fmpz_add(a, a, m);
    }
    int b = _myfmpq_reconstruct_fmpz(res->numer,
                                     res->denom, a, m);
    return b;
  }
}

/* returns 0 if the rational reconstruction failed */
int rat_recon(fmpz_poly_t poly, fmpz_poly_t poly_num, fmpz_poly_t poly_den, fmpz_t lcm,
              data_heap_t *heap, ulong deg, ulong *rec_poly){
  fprintf(stderr, "[%lu, %lu] ", (*rec_poly), deg);
  ulong size = heap->size_mod_array;
  fmpz_mul_ui(heap->h_prod, heap->prod_primes[size-1], heap->primes_array[size-1]);
  int b = myfmpq_reconstruct_fmpz(heap->coef, poly->coeffs+(*rec_poly), heap->h_prod);
  fmpz_lcm(lcm, lcm, heap->coef->denom);
  fmpz_set(poly_num->coeffs + (*rec_poly), heap->coef->numer);
  fmpz_set(poly_den->coeffs + (*rec_poly), heap->coef->denom);
  if(b==0){
    (*rec_poly) = 0;
    fmpz_one(lcm);
    fprintf(stderr, "--> [%lu, %lu]\n", (*rec_poly), deg);
    return b;
  }
  for(ulong i = (*rec_poly) + 1; i <= deg; i++){
    b = myfmpq_reconstruct_fmpz(heap->coef, poly->coeffs+i, heap->h_prod);
    if(b == 0){
      (*rec_poly) = i-1;
      fprintf(stderr, "--> [%lu, %lu]\n", (*rec_poly), deg);
      return b;
    }
    fmpz_lcm(lcm, lcm, heap->coef->denom);
    fmpz_set(poly_num->coeffs + i, heap->coef->numer);
    fmpz_set(poly_den->coeffs + i, heap->coef->denom);
  }
  (*rec_poly) = deg;
  for(ulong i = 0; i <= deg; i++){
    fmpz_divexact(heap->coef->numer, lcm, poly_den->coeffs + i);
    fmpz_mul(poly_num->coeffs + i, poly_num->coeffs + i, heap->coef->numer);
  }
  fprintf(stderr, "--> [%lu, %lu]\n", (*rec_poly), deg);
  return b;
}


