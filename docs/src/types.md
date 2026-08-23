```@meta
CurrentModule = AlgebraicSolving
DocTestSetup = quote
  using AlgebraicSolving
end
```

```@setup algebraicsolving
using AlgebraicSolving
```

```@contents
Pages = ["types.md"]
```

# Data Types

## Introduction

AlgebraicSolving implements functionality for ideals in multivariate polynomial rings with coefficients in 

- a finite field of cardinality $p$ where $p$ is a prime number $<2^{31}$,
- the field of rational numbers,
- the field of rational fractions of a polynomial ring with base field the rational numbers.

## Polynomial Rings

We use [Nemo](https://www.nemocas.org/)'s multivariate polynomial 
ring structures:

```@repl
using AlgebraicSolving
R, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"], internal_ordering=:degrevlex)
```
The above example defines a multivariate polynomial ring in three variables `x`, 
`y`, and `z` over the rationals using the dgree reverse lexicographical ordering 
for printing polynomials in the following. One can also define polynomial rings 
over finite fields:

```@repl
using AlgebraicSolving
R, (x,y,z) = polynomial_ring(GF(101), ["x", "y", "z"], internal_ordering=:degrevlex)
```
In order to define polynomial rings over rational fraction fields one can proceed as follows:
```@repl
using AlgebraicSolving
R, (t, ) = polynomial_ring(QQ, ["t"])
K = fraction_field(R)
R, (x,y,z) = polynomial_ring(K, ["x", "y", "z"])
```

## Ideals

Ideals can be constructed by giving an array of generators. Ideals cache varies 
data structures connected to ideals in order to make computational algebra more 
effective:

```@repl
using AlgebraicSolving
R, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"], internal_ordering=:degrevlex)
I = Ideal([x+y+1, y*z^2-13*y^2])
```

Ideals in polynomial rings over rational fraction fields are implemented using 
the special type `ParametricIdeal`:
```@repl
using AlgebraicSolving
R, (t, ) = polynomial_ring(QQ, ["t"])
K = fraction_field(R)
R, (x,y,z) = polynomial_ring(K, ["x", "y", "z"])
I = ParametricIdeal([t*x^2 - 1, t*x*y - z + 4])
```
