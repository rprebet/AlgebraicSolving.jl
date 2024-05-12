/*
  Functions below are very slight modifications of FLINT file
  resultant_euclidean.c which is part of FLINT.

  Copyright (C) 2007, 2008 William Hart
  Copyright (C) 2011 Sebastian Pancratz
  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

//taken from FLINT
#define MYMPN_NORM(a, an)                       \
  do {                                          \
    while ((an) != 0 && (a)[(an) - 1] == 0)     \
      (an)--;                                   \
  } while (0)

static void _my_nmod_vec_scalar_mul_nmod_large(mp_ptr res, mp_srcptr vec, 
                                  slong len, mp_limb_t c, nmod_t mod){
  for (slong i = 0; i < len; i++)
    {
      mp_limb_t hi, lo;
      umul_ppmm(hi, lo, vec[i], c);
      NMOD_RED2(res[i], hi, lo, mod); /* hi already reduced mod n */
    }
}


static void _my_nmod_poly_rem_q1(mp_ptr R,
                       mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                       nmod_t mod)
{
  const mp_limb_t invL = (B[lenB-1] == 1) ? 1 : n_invmod(B[lenB-1], mod.n);

  if (lenB > 1)
    {
      mp_limb_t t, q0, q1;

      q1 = n_mulmod2_preinv(A[lenA-1], invL, mod.n, mod.ninv);
      t  = n_mulmod2_preinv(q1, B[lenB-2], mod.n, mod.ninv);
      t  = n_submod(A[lenA-2], t, mod.n);
      q0 = n_mulmod2_preinv(t, invL, mod.n, mod.ninv);

      if (FLINT_BITS + 2 <= 2 * mod.norm)
        {
          mpn_mul_1(R, B, lenB - 1, q0);
          if (lenB > 2)
            mpn_addmul_1(R + 1, B, lenB - 2, q1);
          _nmod_vec_reduce(R, R, lenB - 1, mod);
        }
      else
        {
          _my_nmod_vec_scalar_mul_nmod_large(R, B, lenB - 1, q0, mod);
          if (lenB > 2)
            _nmod_vec_scalar_addmul_nmod(R + 1, B, lenB - 2, q1, mod);
        }

      _nmod_vec_sub(R, A, R, lenB - 1, mod);
    }
}

static void _my_nmod_poly_rem(mp_ptr R, mp_srcptr A, slong lenA, 
                    mp_srcptr B, slong lenB, nmod_t mod)
{
  TMP_INIT;

  if (lenA - lenB == 1)
    {
      _my_nmod_poly_rem_q1(R, A, lenA, B, lenB, mod);
    }
  else if (lenA < FMPZ_POLY_SQRTREM_DIVCONQUER_CUTOFF)
    {
      mp_ptr W;

      TMP_START;
      W = TMP_ALLOC(NMOD_DIVREM_BC_ITCH(lenA, lenB, mod)*sizeof(mp_limb_t));

      _nmod_poly_rem_basecase(R, W, A, lenA, B, lenB, mod);
      TMP_END;
    }
  else
    {
      mp_ptr Q = _nmod_vec_init(lenA - lenB + 1);

      _nmod_poly_divrem(Q, R, A, lenA, B, lenB, mod);
      _nmod_vec_clear(Q);
    }
}

//Les coefficients sont ranges dans le tableau par degres croissants.
//dans param on se retrouvera avec les coefficients de la parametrisation 
//attention si il n'y a pas de polynome de degre 1 dans la suite des sous-resultants
//le contenu de param reste inchange (c'est pour ca qu'il faut l'initialiser a 0)
static mp_limb_t
_nmod_poly_param_euclidean_allocated(mp_srcptr poly1, slong len1,
                                     mp_srcptr poly2, slong len2,
                                     nmod_t mod, mp_ptr w,
                                     mp_ptr head,
                                     mp_ptr tail)
{
  if(DEBUG>3){
    fprintf(stderr, "lengths = [%ld, %ld], prime = %ld\n", len1, len2, mod.n);
  }
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
          head[0] = poly1[1];
          tail[0] = poly1[0];
          return poly2[0];
        }
        else
        {
            return n_powmod2_ui_preinv(poly2[0], len1 - 1, mod.n, mod.ninv);
        }
    }
    else{  /* len1 >= len2 >= 2 */
      if(len2 == 2){
        head[0] = poly2[1];
        tail[0] = poly2[0];
      }
        mp_limb_t res = 1;

        mp_ptr u, v, r, t;
        slong l0, l1, l2;
        mp_limb_t lc;

        u = w;
        v = w + len1;
        r = v + len1;
        //u == poly1
        _nmod_vec_set(u, poly1, len1);
        //v == poly2
        _nmod_vec_set(v, poly2, len2);
        l1 = len1;
        l2 = len2;

        if(DEBUG>5){
          fprintf(stderr, "Affichage de u\n");
          nmod_vec_print(stderr, u, len1);
          fprintf(stderr, "Affichage de v\n");
          nmod_vec_print(stderr, v, len2);
        }

        do
        {
            l0 = l1;
            l1 = l2;
            //coefficient dominant de v
            lc = v[l1 - 1];

            if(DEBUG>5){
              fprintf(stderr, "u = ");
              nmod_vec_print(stderr, u, l0);
              fprintf(stderr, "v = ");
              nmod_vec_print(stderr, v, l1);
            }
            //u est de longueur l0, v est de longueur l1
            //r = remainder(u, v)
            _my_nmod_poly_rem(r, u, l0, v, l1, mod);
            //            _nmod_poly_rem(r, u, l0, v, l1, mod);
            l2 = l1 - 1;
            //determine la longueur de r (decremente l2 tant que le top coeff == 0)
            //r est de degre l2 - 1
            MYMPN_NORM(r, l2);
            if(DEBUG>5){
              fprintf(stderr, "r = ");
              nmod_vec_print(stderr, r, l2);
              fprintf(stderr, "l0 = %ld\n", l0);
              fprintf(stderr, "l1 = %ld\n", l1);
              fprintf(stderr, "l2 = %ld\n", l2);
              fprintf(stderr, "lc = %ld\n", lc);
              fprintf(stderr, "res = %ld\n", res);
            }
            //mise a jour
            {
                t = u;
                u = v;
                v = r;
                r = t;
            }
            if (l2 >= 1){
                lc  = n_powmod2_preinv(lc, l0 - l2, mod.n, mod.ninv);
                res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);

                if (((l0 | l1) & 1) == 0)
                {
                  if(DEBUG>5){
                    fprintf(stderr, "[%ld, %ld] ", l0, l1);
                  }
                    res = nmod_neg(res, mod);
                }
                if(l2==2){ /* degree of the remainder -- which is now stored in v -- is 1 */
                  mp_limb_t hi, lo; //, c;
                  umul_ppmm(hi, lo, v[0], res);
                  NMOD_RED2(tail[0], hi, lo, mod); /* hi already reduced mod n */
                  hi = 0; lo = 0;
                  umul_ppmm(hi, lo, v[1], res);
                  NMOD_RED2(head[0], hi, lo, mod); /* hi already reduced mod n */
                }
                if(l2==1){/* degree of the remainder -- which is now stored in v -- is 0 */
                  /* head[0] = u[1]; */
                  /* tail[0] = u[0]; */
                  /* fprintf(stderr, "--> head = %ld\n", head[0]); */
                  /* fprintf(stderr, "--> tail = %ld\n", tail[0]); */
                }
            }
            else /* l2 == 0 */
            {
                if (l1 == 1)
                {
                    lc  = n_powmod2_preinv(lc, l0 - 1, mod.n, mod.ninv);
                    res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);
                    /* head[0] = u[1]; */
                    /* tail[0] = u[0]; */
                }
                else
                {
                  if(l1==2){
                    fprintf(stderr, "l1 == 2 et lc = %ld\n", lc);
                    lc  = n_powmod2_preinv(lc, l0 - 1, mod.n, mod.ninv);
                    fprintf(stderr, "et maintenant lc = %ld\n", lc);
                    res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);
                    fprintf(stderr, "et res = %ld\n", res);
                    head[0] = u[1];
                    tail[0] = u[0];
                    /* mp_limb_t hi, lo; //, c; */
                    /* umul_ppmm(hi, lo, u[0], res); */
                    /* NMOD_RED2(tail[0], hi, lo, mod); /\* hi already reduced mod n *\/ */
                    /* hi = 0; lo = 0; */
                    /* umul_ppmm(hi, lo, u[1], res); */
                    /* NMOD_RED2(head[0], hi, lo, mod); /\* hi already reduced mod n *\/ */
                  }
                  else{
                    head[0] = 0;
                    tail[0] = 0;
                  }
                  res = 0;
                }
            }
        }
        while (l2 > 0);

        return res;
    }
}


//Les coefficients sont ranges dans le tableau par degre croissant.
static mp_limb_t
_nmod_poly_resultant_euclidean_allocated(mp_srcptr poly1, slong len1, 
                                         mp_srcptr poly2, slong len2, nmod_t mod, mp_ptr w)
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

        mp_ptr u, v, r, t;
        slong l0, l1, l2;
        mp_limb_t lc;

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

            //u est de longueur l0, v est de longueur l1
            //r = remainder(u, v)
            _my_nmod_poly_rem(r, u, l0, v, l1, mod);
            l2 = l1 - 1;
            //determine la longueur de r (decremente l2 tant que le top coeff == 0)
            //r est de degre l2 - 1
            MYMPN_NORM(r, l2);
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

        return res;
    }
}


/*
  END OF FLINT modified resultant_euclidean.c
 */
