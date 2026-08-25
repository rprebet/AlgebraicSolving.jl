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
Pages = ["groebner-bases.md"]
```

# Gröbner bases

## Introduction

AlgebraicSolving allows to compute Gröbner bases for input generators over prime
fields of characteristic smaller $2^{31}$, over the rationals w.r.t. the degree
reverse lexicographical monomial order and over fields of rational fractions of
polynomial rings w.r.t. the degree reverse lexicographical monomial order.

At the moment different variants of Faugère's F4 Algorithm are implemented as
well as a signature based algorithm to compute Gröbner bases.

## Functionality

```@docs
    groebner_basis(
        I::Ideal{T} where T <: MPolyRingElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        eliminate::Int=0,
        intersect::Bool=true,
        complete_reduction::Bool=true,
        info_level::Int=0
        )
```

```@docs
	groebner_basis(
		I::ParametricIdeal{K};
		retry::Int=10,
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		) where {K<:FracFieldElem}
```

The engine supports the elimination of one block of variables considering the
product monomial ordering of two blocks, both ordered w.r.t. the degree
reverse lexicographical order. One can either directly add the number of
variables of the first block via the `eliminate` parameter in the
`groebner_basis` call. By using `intersect=false` it is possible to only 
use block ordering without intersecting. We have also implemented an alias 
for this call:

```@docs
    eliminate(
        I::Ideal{T} where T <: MPolyRingElem,
        eliminate::Int;
        intersect::Bool=true,
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        complete_reduction::Bool=true,
        info_level::Int=0
        )
```

To compute signature Gröbner bases use

```@docs
    sig_groebner_basis(sys::Vector{T}; info_level::Int = 0, degbound::Int = 0, mod_ord::Symbol=:POT) where {T <: MPolyRingElem}
```

AlgebraicSolving also allows the computation of multiplication
matrices of polynomials modulo polynomial ideals. For a polynomial
ideal $I$ in a polynomial ring $R$ and a polynomial $g \in R$, this is
a matrix representing the linear map sending an element $f\in R / I$
to $gf \in R/I$.

```@docs
	multiplication_matrix(gb::Vector{T}, b::Vector{T}, g::T) where {T<:MPolyRingElem}
```

```@docs
	multiplication_matrix(
		I::ParametricIdeal{K},
		g::MPoly{K};
		retry::Int=10,
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		) where {K<:FracFieldElem}
```
