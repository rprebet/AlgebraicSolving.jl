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
Pages = ["solvers.md"]
```

# Algebraic Systems Solving

## Introduction

AlgebraicSolving allows to compute solutions of polynomial systems
with input generators over the rationals.  In addition, the
computation of rational parametrizations of solution sets is
implemented for polynomial systems with input generators over finite
fields of order $p$ where $p$ is a prime number $<2^{31}$, over the
rationals and over fraction fields of polynomial rings.

The underlying engine is provided by msolve.

## Functionality

```@docs
    rational_parametrization(
        I::Ideal{T} where T <: MPolyRingElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        info_level::Int=0,
        precision::Int=32

        )
		
	rational_parametrization(
		I::ParametricIdeal{K} where K<:FracFieldElem;
		retry::Int=10,
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		)

    real_solutions(
        I::Ideal{T} where T <: MPolyRingElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        info_level::Int=0,
        precision::Int=32
        )
    rational_solutions(
        I::Ideal{T} where T <: MPolyRingElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        info_level::Int=0,
        precision::Int=32
        )

    curve_rational_parametrization(
        I::Ideal{P} where P<:QQMPolyRingElem;
        info_level::Int=0,
        use_lfs::Bool = false,
        cfs_lfs::Vector{Vector{ZZRingElem}} = Vector{ZZRingElem}[],
        nr_thrds::Int=1,
        check_gen::Bool = true
    )
```

