# DECUHR.jl

Pure-Julia port of the **DECUHR** algorithm (Espelid & Genz, 1994) for
automatic adaptive integration of functions with **vertex singularities**
over hyper-rectangular regions.

Exposed as a pluggable algorithm for the
[Integrals.jl](https://docs.sciml.ai/Integrals/stable/) solver stack via
the `SciMLBase.AbstractIntegralAlgorithm` interface.

## Features

- 2-D and 3-D integration on hyper-rectangles.
- Vertex singularity handling with user-supplied or auto-estimated
  exponent `α` (the singular strength).
- Logarithmic singularities (`logf = k` for `(log)^k` weights).
- Vector-valued integrands (any `NUMFUN`).
- Richardson extrapolation on sub-region averages.
- Reports a `retcode` compatible with Integrals.jl (`Success`,
  `MaxIters`, …).

## Installation

While the package is private in the MicMacTools organisation:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/MicMacTools/DECUHR.jl.git")
```

Once made public and registered, `Pkg.add("DECUHR")` will suffice.

## Quick start

```julia
using Integrals, DECUHR

# ∫₀¹∫₀¹ (x·y)^(-0.5) dx dy = 4
f = (u, _) -> (u[1] * u[2])^(-0.5)
prob = IntegralProblem(f, zeros(2), ones(2))
sol  = solve(prob, DecuhrAlgorithm(singul = 2, alpha = -0.5); abstol = 1e-8)

@show sol.u              # ≈ 4.0
@show sol.retcode        # Success
```

A broader set of examples is in
[`examples/basic_usage.jl`](./examples/basic_usage.jl):

```shell
julia --project=. examples/basic_usage.jl
```

## Tests

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Documentation

```shell
julia --project=docs -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

HTML output is placed in `docs/build/`.

## References

- T.O. Espelid and A. Genz, *DECUHR: An Algorithm for Automatic
  Integration of Singular Functions over a Hyperrectangular Region*,
  Numerical Algorithms **8** (1994), 201–220.
- T.O. Espelid, *On integrating Vertex Singularities using
  Extrapolation*, BIT **34** (1994), 62–79.

## Licence

MIT — see [LICENSE](./LICENSE).

## Author

Jean-François Barthélémy — Cerema / UMR MCD.
