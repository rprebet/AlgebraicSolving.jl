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
Pages = ["real-algebraic-geometry.md"]
```

# Real Algebraic Geometry

## Introduction

AlgebraicSolving allows to work with real semi-algebraic sets, i.e. subsets of real affine space
defined by polynomial equations and inequalities.

## Functionality

AlgebraicSolving implements finding sample points per connected
components of semi-algebraic sets:

```@docs
    points_per_components(
		eqs::Vector{QQMPolyRingElem},
		pos::Vector{QQMPolyRingElem},
		ineqs::Vector{QQMPolyRingElem};
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		info_level::Int=0
		)
```

AlgebraicSolving allows to compute a roadmap for the real trace of the zero-set
of the ideal spanned by given input generators over the rationals.

It assumes that the underlying algebraic set is **smooth**, and its real trace is **bounded**.

The underlying engine is provided by msolve.

```@docs
    roadmap(
      I::Ideal{P} where P <: QQMPolyRingElem;
      C::Vector{Vector{Vector{QQFieldElem}}}=Vector{Vector{QQFieldElem}}[],
      info_level::Int=0
    )
```

In addition, AlgebraicSolving can compute equations definition critical loci of polynomial maps over the given algebraic set.

```@docs
   computepolar(
        J::Union{Vector{Int},UnitRange{Int}},
        V::Ideal{P};
        phi::Vector{P} = P[],
        dimproj = length(J)-1,
        only_mins = false
    ) where (P <: MPolyRingElem)
```


AlgebraicSolving also implements functionality for real root classification:

```@docs
    real_root_classification(
		I::ParametricIdeal{K} where {K<:FracFieldElem},
		g::Vector{MPoly{K}} where {K<:FracFieldElem};
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		info_level::Int=0,
		show_progress::Bool=info_level >= 1,
		output_form::Symbol=:sign,
		ignore_no_real_roots::Bool=false,
		)
```

This is based on the computation of Hermite matrices:

```@docs
    hermite_matrix(gb::Vector{T}, b::Vector{T}, g::T) where {T<:MPolyRingElem}
```

```@docs
    hermite_matrix(
		I::ParametricIdeal{K},
		g::MPoly{K};
		retry::Int=10,
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		) where {K<:FracFieldElem}
```

```@docs
	hermite_matrix(
		g::Vector{MPoly{K}},
		α::Vector{Int},
		I::ParametricIdeal{K};
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		) where {K<:FracFieldElem}
```

```@docs
	hermite_matrix(
		g::Vector{MPoly{K}},
		α::Vector{Int},
		I::ParametricIdeal{K},
		vals::Vector{QQFieldElem};
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		) where {K<:FracFieldElem}
```
